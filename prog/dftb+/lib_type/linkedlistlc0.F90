!--------------------------------------------------------------------------------------------------!
!  DFTB+: general package for performing fast atomistic simulations                                !
!  Copyright (C) 2018  DFTB+ developers group                                                      !
!                                                                                                  !
!  See the LICENSE file for terms of usage and distribution.                                       !
!--------------------------------------------------------------------------------------------------!

#:include 'linkedlist.fypp'

!> Linked list for single strings
module linkedlistlc0
  use accuracy, only : lc
  use assert
  implicit none
  private

  $:define_list(&
      & TYPE_NAME='listCharLc',&
      & ITEM_TYPE='character(lc)',&
      & PADDING='""')

end module linkedlistlc0
