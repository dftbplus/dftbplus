!--------------------------------------------------------------------------------------------------!
!  DFTB+: general package for performing fast atomistic simulations                                !
!  Copyright (C) 2018  DFTB+ developers group                                                      !
!                                                                                                  !
!  See the LICENSE file for terms of usage and distribution.                                       !
!--------------------------------------------------------------------------------------------------!

#:include 'linkedlist.fypp'

!> Linked list for single integers
module dftbp_linkedlisti0
  use dftbp_assert
  implicit none
  private

  $:define_list(&
      & TYPE_NAME='listInt',&
      & ITEM_TYPE='integer',&
      & PADDING='0')

end module dftbp_linkedlisti0
