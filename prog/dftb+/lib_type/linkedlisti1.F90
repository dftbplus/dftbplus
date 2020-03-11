!--------------------------------------------------------------------------------------------------!
!  DFTB+: general package for performing fast atomistic simulations                                !
!  Copyright (C) 2006 - 2020  DFTB+ developers group                                               !
!                                                                                                  !
!  See the LICENSE file for terms of usage and distribution.                                       !
!--------------------------------------------------------------------------------------------------!

#:include 'linkedlist.fypp'

!> Linked list for single real vectors
module dftbp_linkedlisti1
  use dftbp_assert
  implicit none
  private

  $:define_list(&
      & TYPE_NAME='TListIntR1',&
      & ITEM_TYPE='integer',&
      & ITEM_RANK=1,&
      & PADDING='0')

end module dftbp_linkedlisti1
