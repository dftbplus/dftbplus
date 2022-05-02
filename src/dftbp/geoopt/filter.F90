!--------------------------------------------------------------------------------------------------!
!  DFTB+: general package for performing fast atomistic simulations                                !
!  Copyright (C) 2006 - 2022  DFTB+ developers group                                               !
!                                                                                                  !
!  See the LICENSE file for terms of usage and distribution.                                       !
!--------------------------------------------------------------------------------------------------!

#:include 'error.fypp'

!> Geometry transformation filter to transform derivatives and stress tensors to
!> a common gradient vector and the apply the resulting displacement to a geometry.
module dftbp_geoopt_filter
  use dftbp_common_accuracy, only : dp
  use dftbp_math_simplealgebra, only : invert33, determinant33
  implicit none
  private
  public :: TFilter, TFilterInput, TFilter_init


  !> Input for the cartesian constraints
  type :: TFilterInput

    !> Transform lattice
    logical :: lattice = .false.

    !> Allow shearing of lattice
    logical :: fixAngles = .false.

    !> Allow compression / stretching of lattice
    logical :: fixLength(3) = .false.

    !> Allow only isotropic deformations
    logical :: isotropic = .false.

    !> Selected atomic indices
    integer, allocatable :: indMovedAtom(:)

  end type TFilterInput


  !> Cartesian transformation constraints
  type :: TFilter

    !> Number of resulting variables
    integer :: nvar

    !> Transform lattice
    logical :: lattice

    !> Allow shearing of lattice
    logical :: fixAngles

    !> Allow compression / stretching of lattice
    logical :: fixLength(3)

    !> Allow only isotropic deformations
    logical :: isotropic

    !> Selected atomic indices
    logical, allocatable :: mask(:)

  contains

    !> Get number of variables defined by this filter
    procedure :: getDimension

    !> Apply displacement to actual geometry
    procedure :: transformStructure

    !> Transform derivatives and stress tensor to common gradient vector
    procedure :: transformDerivative

  end type TFilter

contains


  !> Create new transformation filter
  subroutine TFilter_init(this, input, coord0, latVec, stat)

    !> Instance of the transformation filter
    type(TFilter), intent(out) :: this

    !> Input to create transformation filter
    type(TFilterInput), intent(in) :: input

    !> Initial coordinates in the central cell
    real(dp), intent(in) :: coord0(:, :)

    !> Initial lattice vectors
    real(dp), intent(in), optional :: latVec(:, :)

    !> Status of operation
    integer, intent(out), optional :: stat

    integer :: ii, iat

    if (allocated(input%indMovedAtom)) then
      allocate(this%mask(size(coord0)))
      do ii = 1, size(coord0)
        iat = (ii - 1)/3 + 1
        this%mask(ii) = any(input%indMovedAtom == iat)
      end do
      this%nvar = count(this%mask)
      if (count(this%mask) == size(coord0)) then
        deallocate(this%mask)
      end if
    else
      this%nvar = size(coord0)
    end if
    this%lattice = input%lattice .and. present(latVec)
    this%fixAngles = input%fixAngles
    this%fixLength = input%fixLength
    this%isotropic = input%isotropic
    if (this%lattice .and. present(latVec)) then
      if (size(coord0) > this%nvar) then
        @:ERROR_HANDLING(stat, 1, &
            & "Subset of optimising atoms not currently possible with lattice optimisation.")
      end if
      this%nvar = this%nvar + size(latVec)
    end if
  end subroutine TFilter_init


  !> Get number of variables defined by this filter
  pure function getDimension(this) result(n)

    !> Instance of the transformation filter
    class(TFilter), intent(in) :: this

    !> Number of variables
    integer :: n

    n = this%nvar
  end function getDimension


  !> Apply displacement to actual geometry
  subroutine transformStructure(this, coord0, latVec, displacement)

    !> Instance of the transformation filter
    class(TFilter), intent(in) :: this

    !> Cartesian coordinates in central cell
    real(dp), intent(inout) :: coord0(:, :)

    !> Lattice vectors
    real(dp), intent(inout), optional :: latVec(:, :)

    !> Displacement vector
    real(dp), intent(in) :: displacement(:)

    real(dp), parameter :: zvec(3) = 0.0_dp

    if (allocated(this%mask)) then
      coord0(:, :) = coord0 + unpack(displacement(:count(this%mask)), &
          & reshape(this%mask, shape(coord0)), spread(zvec, 2, size(coord0, 2)))
    else
      coord0(:, :) = coord0 + reshape(displacement(:size(coord0)), shape(coord0))
      if (this%lattice) then
        latVec(:, :) = latVec + reshape(displacement(size(coord0)+1:), shape(latVec))
      end if
    end if
  end subroutine transformStructure


  !> Transform derivatives and stress tensor to common gradient vector
  subroutine transformDerivative(this, coord0, latVec, gradient, stress, deriv)

    !> Instance of the transformation filter
    class(TFilter), intent(in) :: this

    !> Cartesian coordinates in central cell
    real(dp), intent(in) :: coord0(:, :)

    !> Lattice vectors
    real(dp), intent(in), optional :: latVec(:, :)

    !> Current gradient vector
    real(dp), intent(in) :: gradient(:, :)

    !> Current stress tensor
    real(dp), intent(in) :: stress(:, :)

    !> Common gradient vector containing both cartesian and lattice derivatives
    real(dp), intent(out) :: deriv(:)

    real(dp) :: inv_lat(3, 3), lat_grad(3, 3), vol, sigma(3, 3), pressure
    logical :: mask(3, 3)
    logical, parameter :: diagonal(3, 3) = reshape(&
       & [.true., .false., .false., .false., .true., .false., .false., .false., .true.], &
       & shape(diagonal))

    if (allocated(this%mask)) then
      deriv(:) = pack(gradient, reshape(this%mask, shape(gradient)))
    else
      deriv(:size(gradient)) = reshape(gradient, [size(gradient)])
      if (this%lattice) then
        vol = abs(determinant33(latVec))
        sigma(:, :) = -vol*stress - matmul(gradient, transpose(coord0))
        mask(:, :) = .not.(spread(this%fixLength, 1, 3) .and. spread(this%fixLength, 2, 3))
        if (this%fixAngles) mask = mask .and. diagonal
        where(.not.mask)
          sigma = 0.0_dp
        end where
        if (this%isotropic) then
          pressure = sum(pack(sigma, diagonal)) / 3.0_dp
          where(diagonal)
            sigma = pressure
          elsewhere
            sigma = 0.0_dp
          end where
        endif
        call invert33(inv_lat, latVec)
        inv_lat(:, :) = transpose(inv_lat)
        lat_grad(:, :) = matmul(sigma, inv_lat)
        deriv(size(gradient)+1:) = reshape(lat_grad, [size(lat_grad)])
      end if
    end if
  end subroutine transformDerivative


end module dftbp_geoopt_filter
