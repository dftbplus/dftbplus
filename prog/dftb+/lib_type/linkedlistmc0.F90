!--------------------------------------------------------------------------------------------------!
!  DFTB+: general package for performing fast atomistic simulations                                !
!  Copyright (C) 2018  DFTB+ developers group                                                      !
!                                                                                                  !
!  See the LICENSE file for terms of usage and distribution.                                       !
!--------------------------------------------------------------------------------------------------!

#:include 'linkedlist.fypp'

!> Linked list for single charcter strings
module linkedlistmc0
  use accuracy, only : mc
  use assert
  implicit none
  private

  $:define_list(&
      & TYPE_NAME='listCharMc',&
      & ITEM_TYPE='character(mc)',&
      & PADDING='""')

end module linkedlistmc0
