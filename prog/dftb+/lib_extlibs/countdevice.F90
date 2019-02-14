
!--------------------------------------------------------------------------------------------------!
!  DFTB+: general package for performing fast atomistic simulations                                !
!  Copyright (C) 2018  DFTB+ developers group                                                      !
!                                                                                                  !
!  See the LICENSE file for terms of usage and distribution.                                       !
!--------------------------------------------------------------------------------------------------!

#:include 'common.fypp'

!> Exporting device functionality if compiled with gpu support, otherwise empty.
module countdevice
#:if WITH_GPU
  use count_device_module
#:endif
  implicit none
  public

end module countdevice
