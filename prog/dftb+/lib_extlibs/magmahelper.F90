!--------------------------------------------------------------------------------------------------!
!  DFTB+: general package for performing fast atomistic simulations                                !
!  Copyright (C) 2018  DFTB+ developers group                                                      !
!                                                                                                  !
!  See the LICENSE file for terms of usage and distribution.                                       !
!--------------------------------------------------------------------------------------------------!

!> MAGMA GPU interface library
module magmahelper
  use device_info
  implicit none

#:if WITH_GPU

  !> Whether code was built with GPU support
  logical, parameter :: withGPU = .true.

#:else

  !> Whether code was built with GPU support
  logical, parameter :: withGPU = .false.

#:endif

end module magmahelper
