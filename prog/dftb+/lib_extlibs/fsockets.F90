!--------------------------------------------------------------------------------------------------!
!  DFTB+: general package for performing fast atomistic simulations                                !
!  Copyright (C) 2006 - 2020  DFTB+ developers group                                               !
!                                                                                                  !
!  See the LICENSE file for terms of usage and distribution.                                       !
!--------------------------------------------------------------------------------------------------!

!> Socket interface library
module dftbp_fsockets
  use fsockets
  implicit none

#:if WITH_SOCKETS

  !> Whether code was built with socket support
  logical, parameter :: withSockets = .true.

#:else

  !> Whether code was built with socket support
  logical, parameter :: withSockets = .false.

#:endif

end module dftbp_fsockets
