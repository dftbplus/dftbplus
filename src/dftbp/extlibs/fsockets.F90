!--------------------------------------------------------------------------------------------------!
!  DFTB+: general package for performing fast atomistic simulations                                !
!  Copyright (C) 2006 - 2025  DFTB+ developers group                                               !
!                                                                                                  !
!  See the LICENSE file for terms of usage and distribution.                                       !
!--------------------------------------------------------------------------------------------------!

!> Socket interface library
module dftbp_extlibs_fsockets
  use fsockets, only : close_socket, connect_inet_socket, connect_unix_socket, readbuffer,&
      & writebuffer
  implicit none

#:if WITH_SOCKETS

  !> Whether code was built with socket support
  logical, parameter :: withSockets = .true.

#:else

  !> Whether code was built with socket support
  logical, parameter :: withSockets = .false.

#:endif

end module dftbp_extlibs_fsockets
