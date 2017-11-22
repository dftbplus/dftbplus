!--------------------------------------------------------------------------------------------------!
!  DFTB+: general package for performing fast atomistic simulations                                !
!  Copyright (C) 2017  DFTB+ developers group                                                      !
!                                                                                                  !
!  See the LICENSE file for terms of usage and distribution.                                       !
!--------------------------------------------------------------------------------------------------!

#:include 'common.fypp'

!> H5 H-bond correction dispersion model.
module h5correction
  use accuracy
  implicit none
  private

  public :: H5Inp, H5Corr, H5Corr_init

  !> Input structure for the H5 correction
  type :: H5Inp

    !> Global parameters
    real(dp) :: rscale, wscale

    !> Species names for H5 correction
    character(mc), allocatable :: species_name(:)

  end type H5Inp


  !> Internal data of the H5 correction
  type :: H5Corr

    !> Global parameters
    real(dp) :: rscale, wscale

    !> Species names for H5 correction
    character(mc), allocatable :: species_name(:)

  contains

    procedure :: printH5Setup

  end type H5Corr

contains

  !> Inits a H5Corr instance.
  subroutine H5Corr_init(this, inp)

    !> Initialised instance at return.
    type(H5Corr), intent(out) :: this

    !> Specific input parameters for H5 correction
    type(H5Inp), intent(in) :: inp

    this%rscale = inp%rscale
    this%wscale = inp%wscale

  end subroutine H5Corr_init

  subroutine printH5Setup(this)
     class(H5Corr), intent(inout) :: this
      write(37,*) "H5 setup:"
      write(37,*) "   rscale = ", this%rscale
      write(37,*) "   wscale = ", this%wscale
  end subroutine printH5Setup

end module h5correction
