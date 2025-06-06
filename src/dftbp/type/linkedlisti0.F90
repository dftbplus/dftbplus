!--------------------------------------------------------------------------------------------------!
!  DFTB+: general package for performing fast atomistic simulations                                !
!  Copyright (C) 2006 - 2025  DFTB+ developers group                                               !
!                                                                                                  !
!  See the LICENSE file for terms of usage and distribution.                                       !
!--------------------------------------------------------------------------------------------------!

#:include 'linkedlist.fypp'

!> Linked list for single integers
module dftbp_type_linkedlisti0
  implicit none

  private

  $:define_list(&
      & TYPE_NAME='TListInt',&
      & ITEM_TYPE='integer',&
      & PADDING='0')

end module dftbp_type_linkedlisti0
