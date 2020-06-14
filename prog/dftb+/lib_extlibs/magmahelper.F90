!--------------------------------------------------------------------------------------------------!
!  DFTB+: general package for performing fast atomistic simulations                                !
!  Copyright (C) 2006 - 2020  DFTB+ developers group                                               !
!                                                                                                  !
!  See the LICENSE file for terms of usage and distribution.                                       !
!--------------------------------------------------------------------------------------------------!

#:include 'common.fypp'

!> MAGMA GPU interface library
module dftbp_magmahelper
#:if WITH_GPU
  use device_info
#:endif
  implicit none

#:if WITH_GPU

  !> Whether code was built with GPU support
  logical, parameter :: withGPU = .true.

#:else

  !> Whether code was built with GPU support
  logical, parameter :: withGPU = .false.

#:endif

end module dftbp_magmahelper
