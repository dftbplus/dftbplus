!--------------------------------------------------------------------------------------------------!
!  DFTB+: general package for performing fast atomistic simulations                                !
!  Copyright (C) 2006 - 2025  DFTB+ developers group                                               !
!                                                                                                  !
!  See the LICENSE file for terms of usage and distribution.                                       !
!--------------------------------------------------------------------------------------------------!

#:include 'common.fypp'
#:set pure = "" if defined('WITH_ASSERT') else "pure"
!> Routines to calculate the radial part of a Slater type orbital (STO).
module dftbp_wavegrid_slater
  use dftbp_common_accuracy, only : dp
  implicit none

  private

  public :: TSlaterOrbital
  

  !> Data type for STOs.
  type TSlaterOrbital
    !> Angular momentum of the orbital
    integer :: angMom

    !> Square of the Cutoff, after which the orbital is assumed to be zero
    real(dp) :: cutoffSq

    !> Whether to use cached values instead of direct calculation.
    logical :: useRadialLut

    ! ##  LUT parameters ##
    !> Grid distance (resolution)
    real(dp) :: gridDist

    !> Inverse of the grid distance
    real(dp) :: invLutStep

    !> Number of grid points
    integer :: nGrid

    !> STO values on the distance grid
    real(dp), allocatable :: gridValue(:)

    ! ## Direct calculation STO parameters ##
    !> Maximum power of the radial distance
    integer :: nPow

    !> Number of exponential coefficients
    integer :: nAlpha

    !> Summation coefficients. Shape: [nPow, nAlpha]
    real(dp), allocatable :: aa(:,:)

    !> Exponential coefficients
    real(dp), allocatable :: alpha(:)

  contains

    !> Initialises a SlaterOrbital.
    procedure :: init => TSlaterOrbital_init

    !> Initialises using a LUT for radial values.
    procedure :: initFromLut => TSlaterOrbital_initFromLut

    !> Returns the value of the Lut in a given point.
    procedure :: getRadialCached => TSlaterOrbital_getRadialValueCached
  
    !> Non-interpolated version using direct STO calculation.
    procedure :: getRadialDirect => TSlaterOrbital_getRadialValueDirect

    !> Dispatch to cached or direct version based on this%useRadialLut
    procedure :: getRadial
  
  end type TSlaterOrbital


contains

  !> Initialises using a LUT for radial values.
  subroutine TSlaterOrbital_initFromLut(this, gridValue, gridDist, angMom)
    class(TSlaterOrbital), intent(inout) :: this
    real(dp), intent(in) :: gridValue(:)
    real(dp), intent(in) :: gridDist
    integer, intent(in) :: angMom
    real(dp) :: cutoff

    this%useRadialLut = .true.
    this%angMom = angMom
    this%gridDist = gridDist
    this%invLutStep = 1.0_dp / gridDist
    this%nGrid = size(gridValue)

    cutoff = this%gridDist * real(this%nGrid - 1, dp)
    this%cutoffSq = cutoff**2

    allocate(this%gridValue(this%nGrid))
    this%gridValue(:) = gridValue(:)

  end subroutine TSlaterOrbital_initFromLut

  !> Initialises a SlaterOrbital.
  subroutine TSlaterOrbital_init(this, aa, alpha, ll, resolution, cutoff, useRadialLut)

    !> SlaterOrbital instance to initialise
    class(TSlaterOrbital), intent(inout) :: this

    !> Summation coefficients (nCoeffPerAlpha, nAlpha)
    real(dp), intent(in) :: aa(:,:)

    !> Exponential coefficients
    real(dp), intent(in) :: alpha(:)

    !> Angular momentum of the orbital
    integer, intent(in) :: ll

    !> Grid distance for the orbital
    real(dp), intent(in) :: resolution

    !> Cutoff, after which orbital is assumed to be zero
    real(dp), intent(in) :: cutoff

    !> Whether to use the cached grid for evaluation
    logical, intent(in), optional :: useRadialLut

    integer :: iGrid
    real(dp) :: r

    @:ASSERT(cutoff > 0.0_dp)

    this%angMom = ll
    this%cutoffSq = cutoff ** 2

    ! Store parameters for direct calculation
    this%nAlpha = size(alpha)
    @:ASSERT(size(aa, dim=2) == this%nAlpha)
    this%nPow = size(aa, dim=1)
    allocate(this%aa(this%nPow, this%nAlpha))
    allocate(this%alpha(this%nAlpha))
    this%aa(:,:) = aa

    ! The STO formula uses exp(-alpha * r).
    ! Directly store -alpha to avoid repated negation when calculating.
    this%alpha(:) = -1.0_dp * alpha

    if (present(useRadialLut)) then
      this%useRadialLut = useRadialLut
    else
      this%useRadialLut = .false.
    end if

    if (this%useRadialLut) then
      ! Obtain STO on a grid
      @:ASSERT(resolution > 0.0_dp)
      this%nGrid = floor(cutoff / resolution) + 2
      this%gridDist = resolution
      this%invLutStep = 1.0_dp / resolution

      allocate(this%gridValue(this%nGrid))
      do iGrid = 1, this%nGrid
        r = real(iGrid - 1, dp) * resolution
        this%gridValue(iGrid) = this%getRadialDirect(r)
      end do
    end if

  end subroutine TSlaterOrbital_init


  ${pure}$ function getRadial(this, r) result(sto)
    class(TSlaterOrbital), intent(in) :: this
    real(dp), intent(in) :: r
    real(dp) :: sto

    if (this%useRadialLut) then
      sto = this%getRadialCached(r)
    else
      sto = this%getRadialDirect(r)
    end if

  end function getRadial

  !> Returns the value of the SlaterOrbital at a given point.
  !! Builds a 1d cache grid across which the result is interpolated
  !! in order to speed up evaluation for subsequent calls.
  ${pure}$ function TSlaterOrbital_getRadialValueCached(this, r) result(sto)

    !> SlaterOrbital instance
    class(TSlaterOrbital), intent(in) :: this

    !> Distance, where STO should be calculated
    real(dp), intent(in) :: r

    !> Contains the value of the function on return
    real(dp) :: sto

    integer :: ind
    real(dp) :: frac, posOnGrid

    @:ASSERT(r >= 0.0_dp)

    ! ind = 1 means zero distance as r = (ind - 1) * gridDist
    posOnGrid = r * this%invLutStep
    ind = floor(posOnGrid) + 1
    if (ind < this%nGrid) then
      frac = posOnGrid - real(ind - 1, dp)
      sto = (1.0_dp - frac) * this%gridValue(ind) + frac * this%gridValue(ind+1)
    else
      sto = 0.0_dp
    end if

  end function TSlaterOrbital_getRadialValueCached


  !> Calculates the value of an STO analytically.
  ${pure}$ function TSlaterOrbital_getRadialValueDirect(this, r) result(sto)

    !> SlaterOrbital instance
    class(TSlaterOrbital), intent(in) :: this

    !> Distance, where the STO should be calculated
    real(dp), intent(in) :: r

    !> Value of the STO on return
    real(dp) :: sto

    real(dp) :: pows(this%nPow)
    real(dp) :: rTmp
    integer :: ii, jj

    ! Avoid 0.0**0 as it may lead to arithmetic exception on pre-2008 compilers
    if (this%angMom == 0 .and. r < epsilon(1.0_dp)) then
      rTmp = 1.0_dp
    else
      rTmp = r**this%angMom
    end if

    do ii = 1, this%nPow
      pows(ii) = rTmp
      rTmp = rTmp * r
    end do
    sto = 0.0_dp
    do ii = 1, this%nAlpha
      rTmp = 0.0_dp
      do jj = 1, this%nPow
        rTmp = rTmp + this%aa(jj, ii) * pows(jj)
      end do
      sto = sto + rTmp * exp(this%alpha(ii) * r)
    end do

  end function TSlaterOrbital_getRadialValueDirect


end module dftbp_wavegrid_slater
