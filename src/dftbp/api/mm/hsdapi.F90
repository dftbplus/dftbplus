!--------------------------------------------------------------------------------------------------!
!  DFTB+: general package for performing fast atomistic simulations                                !
!  Copyright (C) 2006 - 2022  DFTB+ developers group                                               !
!                                                                                                  !
!  See the LICENSE file for terms of usage and distribution.                                       !
!--------------------------------------------------------------------------------------------------!

#:include 'common.fypp'

!> Proviedes access to HSD manipulation functions
module dftbp_hsdapi
  use dftbp_extlibs_xmlf90, only : fnode, fNodeList
  use dftbp_io_hsdparser, only : dumpHsd
  use dftbp_io_hsdutils, only : getChild, getChildValue, getChildren, setChild, setChildValue
  implicit none
  private

  public :: fnode, fnodeList
  public :: getChild, getChildren, setChild, getChildValue, setChildValue
  public :: dumpHsd

end module dftbp_hsdapi
