!--------------------------------------------------------------------------------------------------!
!  DFTB+: general package for performing fast atomistic simulations                                !
!  Copyright (C) 2006 - 2020  DFTB+ developers group                                               !
!                                                                                                  !
!  See the LICENSE file for terms of usage and distribution.                                       !
!--------------------------------------------------------------------------------------------------!

#:include "common.fypp"

!> Exports the entities used from the ChIMES calculator library
module dftbp_chimes
#:if WITH_CHIMES
  use chimes_serial08, only : TChimesCalc => ChimesCalc, TChimesCalc_init => ChimesCalc_init
#:endif
  implicit none

  private
  public :: withChimes
  public :: TChimesCalc
#:if WITH_CHIMES
  public :: TChimesCalc_init
#:endif


  !> Whether the code was built with ChIMES support
  logical, parameter :: withChimes = ${FORTRAN_LOGICAL(WITH_CHIMES)}$


#:if not WITH_CHIMES

  !> Dummy placeholder type
  type :: TChimesCalc
  end type TChimesCalc

#:endif

end module dftbp_chimes
