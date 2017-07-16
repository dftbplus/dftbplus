!--------------------------------------------------------------------------------------------------!
!  DFTB+: general package for performing fast atomistic simulations                                !
!  Copyright (C) 2017  DFTB+ developers group                                                      !
!                                                                                                  !
!  See the LICENSE file for terms of usage and distribution.                                       !
!--------------------------------------------------------------------------------------------------!

#:include 'linkedlist.fypp'

module linkedlistlc0
  use accuracy, only : lc
  implicit none
  private

  $:define_list(&
      & TYPE_NAME='listCharLc',&
      & ITEM_TYPE='character(lc)',&
      & PADDING='""')

end module linkedlistlc0
