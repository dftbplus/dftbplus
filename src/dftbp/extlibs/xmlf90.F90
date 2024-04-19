!--------------------------------------------------------------------------------------------------!
!  DFTB+: general package for performing fast atomistic simulations                                !
!  Copyright (C) 2006 - 2023  DFTB+ developers group                                               !
!                                                                                                  !
!  See the LICENSE file for terms of usage and distribution.                                       !
!--------------------------------------------------------------------------------------------------!

!> Exporting the functionality we use from the library XMLF90.
module dftbp_extlibs_xmlf90
  use xmlf90_flib_dom
  use xmlf90_flib_wxml
  use xmlf90_strings
  implicit none
  public

  public :: xmlf_t

end module dftbp_extlibs_xmlf90
