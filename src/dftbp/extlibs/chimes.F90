!--------------------------------------------------------------------------------------------------!
!  DFTB+: general package for performing fast atomistic simulations                                !
!  Copyright (C) 2006 - 2022  DFTB+ developers group                                               !
!                                                                                                  !
!  See the LICENSE file for terms of usage and distribution.                                       !
!--------------------------------------------------------------------------------------------------!

#:include "common.fypp"

!> Exports the entities used from the ChIMES calculator library
module dftbp_extlibs_chimes
#:if WITH_CHIMES
  use chimescalc_serial08, only : TChimesCalc => ChimesCalc, TChimesCalc_init => ChimesCalc_init
#:endif
  implicit none

  private
  public :: withChimes
  public :: TChimesCalc, TChimesCalc_init


  !> Whether the code was built with ChIMES support
  logical, parameter :: withChimes = ${FORTRAN_LOGICAL(WITH_CHIMES)}$


#:if not WITH_CHIMES

  !> Dummy placeholder type
  type :: TChimesCalc
  end type TChimesCalc

#:endif

contains

#:if not WITH_CHIMES

  subroutine TChimesCalc_init()
  end subroutine TChimesCalc_init

#:endif

end module dftbp_extlibs_chimes
