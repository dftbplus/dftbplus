!--------------------------------------------------------------------------------------------------!
!  DFTB+: general package for performing fast atomistic simulations                                !
!  Copyright (C) 2006 - 2022  DFTB+ developers group                                               !
!                                                                                                  !
!  See the LICENSE file for terms of usage and distribution.                                       !
!--------------------------------------------------------------------------------------------------!

#:include "allocatablelist.fypp"

module alloclisthelper
  implicit none

  @:declare_allocatable_list(TIntList, integer)

contains

  @:implement_allocatable_list(TIntList, integer)

end module alloclisthelper
