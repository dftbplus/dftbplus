!--------------------------------------------------------------------------------------------------!
!  DFTB+: general package for performing fast atomistic simulations                                !
!  Copyright (C) 2006 - 2020  DFTB+ developers group                                               !
!                                                                                                  !
!  See the LICENSE file for terms of usage and distribution.                                       !
!--------------------------------------------------------------------------------------------------!

#:include 'common.fypp'

!> Exports scalapackfx functionality if compiled with scalapack support, otherwise empty.
module dftbp_scalapackfx
#:if WITH_SCALAPACK
  use libscalapackfx_module
#:endif
  implicit none
  public

end module dftbp_scalapackfx
