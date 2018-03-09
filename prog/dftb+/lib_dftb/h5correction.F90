!--------------------------------------------------------------------------------------------------!
!  DFTB+: general package for performing fast atomistic simulations                                !
!  Copyright (C) 2017  DFTB+ developers group                                                      !
!                                                                                                  !
!  See the LICENSE file for terms of usage and distribution.                                       !
!--------------------------------------------------------------------------------------------------!

#:include 'common.fypp'

!> H5 H-bond correction.
module h5correction
  use accuracy
  implicit none
  private

  public :: H5Corr, H5Corr_init

  !> Internal data of the H5 correction
  type :: H5Corr
    !> Global parameters
    real(dp) :: rScale, wScale
    !> Nr. of rpecies
    integer :: nSpecies
    !> Species names for H5 correction
    character(mc), allocatable :: speciesName(:)
    !> Elementwise parameters for species
    real(dp), allocatable :: elementPara(:)

  contains

    procedure :: printH5Setup
    procedure :: getParams
    procedure :: scaleShortGamma
    procedure :: scaleShortGammaDeriv

  end type H5Corr

contains

  !> Initialization of a H5Corr instance.
  subroutine H5Corr_init(this)
    !> Initialised instance at return.
    type(H5Corr), intent(out) :: this
  end subroutine H5Corr_init

  !> Print all the parameters for debugging
  subroutine printH5Setup(this)
    ! Arguments
    class(H5Corr), intent(in) :: this
    ! Local variables
    integer :: iSp1

     write(37,*) "H5 setup:"
     write(37,*) "   rScale = ", this%rScale
     write(37,*) "   wScale = ", this%wScale
     write(37,*) "H5 species parameters:"
     do iSp1 = 1, this%nSpecies
       write(37,*) "   ", this%speciesName(iSp1), " = ", this%elementPara(iSp1)
     end do
  end subroutine printH5Setup

  !> Get H5 parameters for a pair of species
  subroutine getParams(this, iSp1, iSp2, applyCorrection, h5Scaling, sumVDW)
    !> Returns a flag whether the correction is applied to the pair,
    !> and the species-specific parameters

    ! Arguments
    class(H5Corr), intent(in) :: this
    !> Input: specie names
    integer, intent(in) :: iSp1, iSp2
    !> Output: flag indicating if the correction is applied to that pair
    logical, intent(out) :: applyCorrection
    !> Output: pair-specific scaling factor
    real(dp), intent(out) :: h5Scaling
    !> Output: Sum of vdW radii of the pair
    real(dp), intent(out) :: sumVDW

    ! Local variables
    character(mc) :: spName1, spName2
    integer :: iSpHeavy
    character(mc) :: spNameHeavy

    spName1 = this%speciesName(iSp1)
    spName2 = this%speciesName(iSp2)
    applyCorrection = .false.

    ! If there is one hydrogen and one other atom, save the heavy atom,
    ! otherwise return with applyCorrection = false
    if (spName1 .ne. "H" .and. spName2 == "H") then
      iSpHeavy = iSp1
      spNameHeavy = spName1
    else if (spName2 .ne. "H" .and. spName1 == "H") then
      iSpHeavy = iSp2
      spNameHeavy = spName2
    else
      return
    end if
   
    ! For each species the correction is applied to,
    ! the correction is enabled and corresponding
    ! parameters are returned
    if (spNameHeavy == "O" ) then
      ! Correction for OH
      applyCorrection = .true.
      ! Parameters
      h5Scaling = this%elementPara(iSpHeavy)
      sumVDW = 2.72_dp
      return
    end if

    if (spNameHeavy == "N" ) then
      ! Correction for NH
      applyCorrection = .true.
      ! Parameters
      h5Scaling = this%elementPara(iSpHeavy)
      sumVDW = 2.75_dp
      return
    end if

    if (spNameHeavy == "S" ) then
      ! Correction for SH
      applyCorrection = .true.
      ! Parameters
      h5Scaling = this%elementPara(iSpHeavy)
      sumVDW = 3.00_dp
      return
    end if
  end subroutine getParams

  !> Apply the correction to the short-range part of gamma function
  subroutine scaleShortGamma(this, shortGamma, iSp1, iSp2, rab)
    !> Returns modified shortGamma

    ! Arguments
    class(H5Corr), intent(in) :: this
    real(dp), intent(inout) :: shortGamma
    integer, intent(in) :: iSp1, iSp2
    real(dp), intent(in) :: rab
    
    ! Local variables
    real(dp) :: h5Scaling, gauss, sumVDW, fwhm, r0, c
    logical :: applyCorrection

    ! Get parameters for current pair of species
    call this%getParams(iSp1, iSp2, applyCorrection, h5Scaling, sumVDW)

    ! If applicable to the current pair, modify the gamma
    if (applyCorrection) then
      ! Gaussian calculation
      fwhm = this%wScale * sumVDW
      r0 = this%rScale * sumVDW
      c = fwhm / 2.35482_dp
      gauss = exp(-1.0_dp * ((rab*0.5291772083_dp)-r0)**2 / 2.0_dp / c**2) * h5Scaling
      ! Apply the correction to original gamma
      shortGamma = shortGamma * (1.0_dp + gauss) - gauss / rab
    end if
  end subroutine scaleShortGamma

  !> Apply the correction to the derivative of the short-range part of gamma function
  subroutine scaleShortGammaDeriv(this, shortGamma, shortGammaDeriv, iSp1, iSp2, rab)
    !> Returns modified shortGamma derivative

    ! Arguments
    class(H5Corr), intent(in) :: this
    real(dp), intent(in) :: shortGamma
    real(dp), intent(inout) :: shortGammaDeriv
    integer, intent(in) :: iSp1, iSp2
    real(dp), intent(in) :: rab
    
    ! Local variables
    real(dp) :: h5Scaling, gauss, sumVDW, fwhm, r0, c, dgauss, deriv1, deriv2
    logical :: applyCorrection

    ! Get parameters for current pair of species
    call this%getParams(iSp1, iSp2, applyCorrection, h5Scaling, sumVDW)

    ! If applicable to the current pair, modify the gamma
    if (applyCorrection) then
      ! Gaussian calculation
      fwhm = this%wScale * sumVDW
      r0 = this%rScale * sumVDW
      c = fwhm / 2.35482_dp
      gauss = exp(-1.0_dp * ((rab*0.5291772083_dp)-r0)**2 / 2.0_dp / c**2) * h5Scaling
      ! Derivative calculation
      dgauss = -1.0 * (0.5291772083 * ( 0.5291772083 * rab - r0 ) ) / c**2 * gauss
      deriv1 = shortGamma * dgauss + shortGammaDeriv * (1.0+gauss)
      deriv2 = dgauss/rab - (h5Scaling*exp(-1.0*(0.5*(0.5291772083*rab - r0)**2)/c**2))/rab**2
      shortGammaDeriv = deriv1 - deriv2
    end if
  end subroutine scaleShortGammaDeriv

end module h5correction
