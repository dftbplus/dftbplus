!--------------------------------------------------------------------------------------------------!
!  DFTB+: general package for performing fast atomistic simulations                                !
!  Copyright (C) 2006 - 2022  DFTB+ developers group                                               !
!                                                                                                  !
!  See the LICENSE file for terms of usage and distribution.                                       !
!--------------------------------------------------------------------------------------------------!

#:include "pointerlist.fypp"

module ptrlisthelper
  implicit none

  type :: TInt
    integer :: value = 0
  end type TInt

  @:declare_pointer_list(TIntList, type(TInt))

contains

  @:implement_pointer_list(TIntList, type(TInt))

end module ptrlisthelper
