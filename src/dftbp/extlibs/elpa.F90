!--------------------------------------------------------------------------------------------------!
!  DFTB+: general package for performing fast atomistic simulations                                !
!  Copyright (C) 2006 - 2025  DFTB+ developers group                                               !
!                                                                                                  !
!  See the LICENSE file for terms of usage and distribution.                                       !
!--------------------------------------------------------------------------------------------------!

#:include 'common.fypp'


!> Interface wrapper for the ELPA library
!>
module dftbp_extlibs_elpa
#:if WITH_ELPA
  use elpa, only : elpa_allocate, elpa_autotune_deallocate, elpa_autotune_t, elpa_deallocate,&
      & elpa_init, elpa_t, elpa_uninit, ELPA_AUTOTUNE_DOMAIN_COMPLEX, ELPA_AUTOTUNE_DOMAIN_REAL,&
      & ELPA_AUTOTUNE_FAST, ELPA_AUTOTUNE_MEDIUM, ELPA_OK, ELPA_SOLVER_1STAGE, ELPA_SOLVER_2STAGE
#:endif
  implicit none
  private

  public :: elpa_t, withElpa
  public :: elpa_autotune_t
#:if WITH_ELPA
  public :: elpa_allocate, elpa_deallocate
  public :: elpa_init, elpa_uninit
  public :: elpa_autotune_deallocate
  public :: ELPA_OK
  public :: ELPA_SOLVER_1STAGE, ELPA_SOLVER_2STAGE
  public :: ELPA_AUTOTUNE_FAST, ELPA_AUTOTUNE_MEDIUM
  public :: ELPA_AUTOTUNE_DOMAIN_REAL, ELPA_AUTOTUNE_DOMAIN_COMPLEX
#:else
  type elpa_t
  end type
  type elpa_autotune_t
  end type
#:endif

  !> Whether code was built with ELPA support
  logical, parameter :: withElpa = #{if WITH_ELPA}# .true. #{else}# .false. #{endif}#

end module dftbp_extlibs_elpa
