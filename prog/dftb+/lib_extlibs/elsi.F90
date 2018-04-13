!--------------------------------------------------------------------------------------------------!
!  DFTB+: general package for performing fast atomistic simulations                                !
!  Copyright (C) 2017  DFTB+ developers group                                                      !
!                                                                                                  !
!  See the LICENSE file for terms of usage and distribution.                                       !
!--------------------------------------------------------------------------------------------------!

#:include 'common.fypp'

!> Interface wrapper for the ELSI library
!>
module elsiInterface
  implicit none

  public

#:if WITH_ELSI

  !> Whether code was built with ELSI support
  logical, parameter :: withELSI = .true.

#:if WITH_PEXSI

  !> Whether code was built with PEXSI support
  logical, parameter :: withPEXSI = .true.

#:else

  !> Whether code was built with PEXSI support
  logical, parameter :: withPEXSI = .false.

#:endif

#:else

  !> Whether code was built with ELSI support
  logical, parameter :: withELSI = .false.

  !> Whether code was built with PEXSI support
  logical, parameter :: withPEXSI = .false.

#:endif

end module elsiInterface
