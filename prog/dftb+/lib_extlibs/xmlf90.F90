!--------------------------------------------------------------------------------------------------!
!  DFTB+: general package for performing fast atomistic simulations                                !
!  Copyright (C) 2018  DFTB+ developers group                                                      !
!                                                                                                  !
!  See the LICENSE file for terms of usage and distribution.                                       !
!--------------------------------------------------------------------------------------------------!

!> Exporting the functionality we use from the library XMLF90.
module xmlf90
  use m_strings
  use flib_wxml
  use flib_dom
  implicit none
  public

  public :: xmlf_t

end module xmlf90
