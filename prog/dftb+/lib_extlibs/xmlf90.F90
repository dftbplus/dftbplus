!--------------------------------------------------------------------------------------------------!
!  DFTB+: general package for performing fast atomistic simulations                                !
!  Copyright (C) 2006 - 2020  DFTB+ developers group                                               !
!                                                                                                  !
!  See the LICENSE file for terms of usage and distribution.                                       !
!--------------------------------------------------------------------------------------------------!

!> Exporting the functionality we use from the library XMLF90.
module dftbp_xmlf90
  use xmlf90_strings
  use xmlf90_flib_wxml
  use xmlf90_flib_dom
  implicit none
  public

  public :: xmlf_t

end module dftbp_xmlf90
