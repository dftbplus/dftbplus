!--------------------------------------------------------------------------------------------------!
!  DFTB+: general package for performing fast atomistic simulations                                !
!  Copyright (C) 2006 - 2025  DFTB+ developers group                                               !
!                                                                                                  !
!  See the LICENSE file for terms of usage and distribution.                                       !
!--------------------------------------------------------------------------------------------------!

!> Exporting the functionality we use from the library mbd.
module dftbp_extlibs_mbd
  use mbd, only : mbd_calc_t, TDispMbdInp => mbd_input_t
  implicit none

  private
  public :: TDispMbdInp, mbd_calc_t

end module dftbp_extlibs_mbd
