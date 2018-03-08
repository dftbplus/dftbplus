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
    real(dp) :: rscale, wscale
    !> Nr. of rpecies
    integer :: nSpecies
    !> Species names for H5 correction
    character(mc), allocatable :: species_name(:)
    !> Elementwise parameters for species
    real(dp), allocatable :: elementPara(:)

  contains

    procedure :: printH5Setup
    procedure :: getParams
    procedure :: scaleShortGamma
    procedure :: scaleShortGammaDeriv

  end type H5Corr

contains

  !> Inits a H5Corr instance.
  subroutine H5Corr_init(this)
    !> Initialised instance at return.
    type(H5Corr), intent(out) :: this

  end subroutine H5Corr_init

  subroutine printH5Setup(this)
     ! Arguments
     class(H5Corr), intent(in) :: this
     ! Local variables
     integer :: iSp1

      write(37,*) "H5 setup:"
      write(37,*) "   rscale = ", this%rscale
      write(37,*) "   wscale = ", this%wscale
      write(37,*) "H5 species parameters:"
      do iSp1 = 1, this%nSpecies
        write(37,*) "   ", this%species_name(iSp1), " = ", this%elementPara(iSp1)
      end do

  end subroutine printH5Setup

  ! This method gets H5 parameters for a pair of species.
  ! It also returns a logical value do_corr if the correction should be applied
  ! to this pair.
  subroutine getParams(this, iSp1, iSp2, do_corr, h5scaling, sumvdw)
     ! Arguments
     class(H5Corr), intent(in) :: this
     integer, intent(in) :: iSp1, iSp2
     logical, intent(out) :: do_corr
     real(dp), intent(out) :: h5scaling, sumvdw

     ! Local variables
     character(mc) :: spname1, spname2
     integer :: iSpHeavy
     character(mc) :: spnameHeavy

     spname1 = this%species_name(iSp1)
     spname2 = this%species_name(iSp2)
     do_corr = .false.

     ! If there is one hydrogen and one other atom, save the heavy atom,
     ! otherwise return with do_corr = false
     if (spname1 .ne. "H" .and. spname2 == "H") then
             iSpHeavy = iSp1
             spnameHeavy = spname1
     else if (spname2 .ne. "H" .and. spname1 == "H") then
             iSpHeavy = iSp2
             spnameHeavy = spname2
     else
             return
     end if
   
     ! For each species the correction is applied to,
     ! the correction is enabled and corresponding
     ! parameters are returned
     if (spnameHeavy == "O" ) then
             ! Correction for OH
             do_corr = .true.
             ! Parameters
             h5scaling = this%elementPara(iSpHeavy)
             sumvdw = 2.72_dp
             return
     end if

     if (spnameHeavy == "N" ) then
             ! Correction for NH
             do_corr = .true.
             ! Parameters
             h5scaling = this%elementPara(iSpHeavy)
             sumvdw = 2.75_dp
             return
     end if

     if (spnameHeavy == "S" ) then
             ! Correction for SH
             do_corr = .true.
             ! Parameters
             h5scaling = this%elementPara(iSpHeavy)
             sumvdw = 3.00_dp
             return
     end if
  end subroutine getParams

  subroutine scaleShortGamma(this, shortGamma, iSp1, iSp2, rab)
     ! Arguments
     class(H5Corr), intent(in) :: this
     real(dp), intent(inout) :: shortGamma
     integer, intent(in) :: iSp1, iSp2
     real(dp), intent(in) :: rab
     
     ! Local variables
     real(dp) :: h5scaling, gauss, sumvdw, fwhm, r0, c
     logical :: do_corr

     ! Get parameters for current pair of species
     call this%getParams(iSp1, iSp2, do_corr, h5scaling, sumvdw)

     ! If applicable to the current pair, modify the gamma
     if (do_corr) then
             ! Gaussian calculation
             fwhm = this%wscale * sumvdw
             r0 = this%rscale * sumvdw
             c = fwhm / 2.35482_dp
             gauss = exp(-1.0_dp * ((rab*0.5291772083_dp)-r0)**2 / 2.0_dp / c**2) * h5scaling
             ! Apply the correction to original gamma
             ShortGamma = ShortGamma * (1.0_dp + gauss) - gauss / rab
     end if
  end subroutine scaleShortGamma

  subroutine scaleShortGammaDeriv(this, shortGamma, shortGammaDeriv, iSp1, iSp2, rab)
     ! Arguments
     class(H5Corr), intent(in) :: this
     real(dp), intent(in) :: shortGamma
     real(dp), intent(inout) :: shortGammaDeriv
     integer, intent(in) :: iSp1, iSp2
     real(dp), intent(in) :: rab
     
     ! Local variables
     real(dp) :: h5scaling, gauss, sumvdw, fwhm, r0, c, dgauss, deriv1, deriv2
     logical :: do_corr

     ! Get parameters for current pair of species
     call this%getParams(iSp1, iSp2, do_corr, h5scaling, sumvdw)

     ! If applicable to the current pair, modify the gamma
     if (do_corr) then
             ! Gaussian calculation
             fwhm = this%wscale * sumvdw
             r0 = this%rscale * sumvdw
             c = fwhm / 2.35482_dp
             gauss = exp(-1.0_dp * ((rab*0.5291772083_dp)-r0)**2 / 2.0_dp / c**2) * h5scaling
             ! Derivative calculation
             dgauss = -1.0 * (0.5291772083 * ( 0.5291772083 * rab - r0 ) ) / c**2 * gauss

             deriv1 = shortGamma * dgauss + shortGammaDeriv * (1.0+gauss)
             deriv2 = dgauss/rab - (h5scaling*exp(-1.0*(0.5*(0.5291772083*rab - r0)**2)/c**2))/rab**2
             shortGammaDeriv = deriv1 - deriv2
     end if
   end subroutine scaleShortGammaDeriv


end module h5correction
