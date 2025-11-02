
!  DFTB+: general package for performing fast atomistic simulations                                !
!  Copyright (C) 2006 - 2025  DFTB+ developers group                                               !
!                                                                                                  !
!  See the LICENSE file for terms of usage and distribution.                                       !
!--------------------------------------------------------------------------------------------------!

#:include 'common.fypp'
#:set pure = "" if defined('WITH_ASSERT') else "pure"
!> Holds TRadialTableOrbital, a concrete implementation of TOrbital utilising a 1d interpolated lookup table.
!> This allows for arbitrary radial functions.
module dftbp_wavegrid_basis_lut
  use dftbp_common_accuracy, only : dp
  use dftbp_io_message, only : error
  use dftbp_wavegrid_basis_orbital, only : TOrbital
  implicit none

  private

  public :: TRadialTableOrbital, TRadialTableOrbital_initFromArray, TRadialTableOrbital_initFromOrbital
  
  !> Concrete class for an orbital represented by a linearly
  !! interpolated lookup table.
  type, extends(TOrbital) :: TRadialTableOrbital
    !> Grid spacing (resolution)
    real(dp) :: gridDist
    !> Inverse of the grid spacing
    real(dp) :: invLutStep
    !> Orbital values on the grid
    real(dp), allocatable :: gridValue(:)
  contains
    procedure :: getRadial => TRadialTableOrbital_getRadial
    procedure, pass(lhs) :: assign => TRadialTableOrbital_assign
  end type TRadialTableOrbital

contains

  !> Initialises using a verbatim array of values.
  subroutine TRadialTableOrbital_initFromArray(this, gridValue, gridDist, angMom)
    type(TRadialTableOrbital), intent(out) :: this
    real(dp), intent(in) :: gridValue(:)
    real(dp), intent(in) :: gridDist
    integer, intent(in) :: angMom
    real(dp) :: cutoff

    this%angMom = angMom
    this%gridDist = gridDist
    this%invLutStep = 1.0_dp / gridDist

    cutoff = this%gridDist * real(size(gridValue) - 1, dp)
    this%cutoffSq = cutoff**2

    allocate(this%gridValue, source=gridValue)

  end subroutine TRadialTableOrbital_initFromArray

  !> Resamples another orbital onto a LUT with given resolution.
  subroutine TRadialTableOrbital_initFromOrbital(this, other, resolution, newCutoff)
    type(TRadialTableOrbital), intent(out) :: this
    class(TOrbital), intent(in) :: other
    real(dp), intent(in) :: resolution
    !> New cutoff, required to be larger than original cutoff.
    real(dp), intent(in), optional :: newCutoff

    integer :: iGrid
    real(dp) :: r, cutoff, norm

    @:ASSERT(resolution > 0.0_dp)

    ! Optionally enlarge cutoff
    cutoff = sqrt(other%cutoffSq)
    if (present(newCutoff)) then
      @:ASSERT(newCutoff >= cutoff)
      cutoff = newCutoff
    end if

    ! Set parameters
    this%angMom = other%angMom
    this%cutoffSq = other%cutoffSq
    this%gridDist = resolution
    this%invLutStep = 1.0_dp / resolution
    
    ! Allocate LUT grid
    allocate(this%gridValue(floor(cutoff / resolution) + 2))

    ! Populate LUT by sampling the other orbital
    do iGrid = 1, size(this%gridValue)
      r = real(iGrid - 1, dp) * resolution
      this%gridValue(iGrid) = other%getRadial(r)
    end do

  end subroutine TRadialTableOrbital_initFromOrbital


  !> Returns the value of the RadialFunction at a given point.
  !! Builds a 1d cache grid across which the result is interpolated
  !! in order to speed up evaluation for subsequent calls.
  ${pure}$ function TRadialTableOrbital_getRadial(this, r) result(sto)

    !> RadialTable instance
    class(TRadialTableOrbital), intent(in) :: this

    !> Distance, where value should be calculated
    real(dp), intent(in) :: r

    !> Contains the value of the function on return
    real(dp) :: sto

    integer :: ind
    real(dp) :: frac, posOnGrid

    @:ASSERT(r >= 0.0_dp)

    ! ind = 1 means zero distance as r = (ind - 1) * gridDist
    posOnGrid = r * this%invLutStep
    ind = floor(posOnGrid) + 1
    if (ind < size(this%gridValue)) then
      frac = posOnGrid - real(ind - 1, dp)
      sto = (1.0_dp - frac) * this%gridValue(ind) + frac * this%gridValue(ind+1)
    else
      sto = 0.0_dp
    end if

  end function TRadialTableOrbital_getRadial

  !> Assignment operator for TRadialTableOrbital
  subroutine TRadialTableOrbital_assign(lhs, rhs)
    class(TRadialTableOrbital), intent(out) :: lhs
    class(TOrbital), intent(in) :: rhs

    select type (rhs)
      type is (TRadialTableOrbital)
        lhs%angMom = rhs%angMom
        lhs%cutoffSq = rhs%cutoffSq
        lhs%gridDist = rhs%gridDist
        lhs%invLutStep = rhs%invLutStep

        if (allocated(lhs%gridValue)) deallocate(lhs%gridValue)
        allocate(lhs%gridValue, source=rhs%gridValue)
      class default
        call error("Implicit conversion to TRadialTableOrbital not allowed, use initFromOrbital instead.")
    end select

  end subroutine TRadialTableOrbital_assign


end module dftbp_wavegrid_basis_lut
