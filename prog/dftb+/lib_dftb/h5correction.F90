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
  contains

    procedure :: printH5Setup
    procedure :: getParams
    procedure :: scaleShortGamma

  end type H5Corr

contains

  !> Inits a H5Corr instance.
  subroutine H5Corr_init(this)
    !> Initialised instance at return.
    type(H5Corr), intent(out) :: this

  end subroutine H5Corr_init

  subroutine printH5Setup(this)
     class(H5Corr), intent(inout) :: this
      write(37,*) "H5 setup:"
      write(37,*) "   rscale = ", this%rscale
      write(37,*) "   wscale = ", this%wscale
  end subroutine printH5Setup

  ! This method gets H5 parameters for a pair of species.
  ! It also returns a logical value do_corr if the correction should be applied
  ! to this pair.
  subroutine getParams(this, iSp1, iSp2, do_corr, h5scaling, sumvdw)
     ! Arguments
     class(H5Corr), intent(inout) :: this
     integer, intent(in) :: iSp1, iSp2
     logical, intent(out) :: do_corr
     real(dp), intent(out) :: h5scaling, sumvdw

     ! Local variables
     character(mc) :: spname1, spname2

     spname1 = this%species_name(iSp1)
     spname2 = this%species_name(iSp2)
     do_corr = .false.
     
     if ((spname1 == "O" .and. spname2 == "H") .or. (spname1 == "H" .and. spname2 == "O")) then
             ! Correction for OH
             do_corr = .true.
             ! Parameters
             h5scaling = 0.06_dp
             sumvdw = 2.72_dp
     end if

     if ((spname1 == "N" .and. spname2 == "H") .or. (spname1 == "H" .and. spname2 == "N")) then
             ! Correction for NH
             do_corr = .true.
             ! Parameters
             h5scaling = 0.18_dp
             sumvdw = 2.75_dp
     end if

     if ((spname1 == "S" .and. spname2 == "H") .or. (spname1 == "H" .and. spname2 == "S")) then
             ! Correction for SH
             do_corr = .true.
             ! Parameters
             h5scaling = 0.21_dp
             sumvdw = 3.00_dp
     end if
  end subroutine getParams

  subroutine scaleShortGamma(this, shortGamma, iSp1, iSp2, rab)
     ! Arguments
     class(H5Corr), intent(inout) :: this
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
             gauss = exp(-1.0_dp * ((rab*0.5291772083_dp)-r0)**2 / 2.0_dp / c**2)
             gauss = gauss * h5scaling
             ! Apply the correction to original gamma
             ShortGamma = ShortGamma * (1.0_dp + gauss)
             ShortGamma = ShortGamma - gauss / rab
     end if
  end subroutine scaleShortGamma


end module h5correction
