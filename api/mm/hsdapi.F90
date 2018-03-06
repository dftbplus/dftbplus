!--------------------------------------------------------------------------------------------------!
!  DFTB+: general package for performing fast atomistic simulations                                !
!  Copyright (C) 2018  DFTB+ developers group                                                      !
!                                                                                                  !
!  See the LICENSE file for terms of usage and distribution.                                       !
!--------------------------------------------------------------------------------------------------!

#:include 'common.fypp'

!> Proviedes access to HSD manipulation functions
module hsdapi
  use hsdparser, only : dumpHsd
  use hsdutils
  use xmlf90
  implicit none
  private

  public :: fnode, fnodeList
  public :: getChild, getChildren, setChild, getChildValue, setChildValue
  public :: dumpHsd

end module hsdapi
