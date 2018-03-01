!--------------------------------------------------------------------------------------------------!
!  DFTB+: general package for performing fast atomistic simulations                                !
!  Copyright (C) 2018  DFTB+ developers group                                                      !
!                                                                                                  !
!  See the LICENSE file for terms of usage and distribution.                                       !
!--------------------------------------------------------------------------------------------------!

#:include 'linkedlist.fypp'

!> Linked list for single real vectors
module linkedlisti1
  use assert
  implicit none
  private

  $:define_list(&
      & TYPE_NAME='listIntR1',&
      & ITEM_TYPE='integer',&
      & ITEM_RANK=1,&
      & PADDING='0')

end module linkedlisti1
