!--------------------------------------------------------------------------------------------------!
!  DFTB+: general package for performing fast atomistic simulations                                !
!  Copyright (C) 2006 - 2020  DFTB+ developers group                                               !
!                                                                                                  !
!  See the LICENSE file for terms of usage and distribution.                                       !
!--------------------------------------------------------------------------------------------------!

#:include 'common.fypp'

!> Proviedes access to HSD manipulation functions
module dftbp_hsdapi
  use dftbp_hsdparser, only : dumpHsd
  use dftbp_hsdutils
  use dftbp_xmlf90
  implicit none
  private

  public :: fnode, fnodeList
  public :: getChild, getChildren, setChild, getChildValue, setChildValue
  public :: dumpHsd

end module dftbp_hsdapi
