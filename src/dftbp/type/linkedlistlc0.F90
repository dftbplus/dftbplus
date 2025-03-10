!--------------------------------------------------------------------------------------------------!
!  DFTB+: general package for performing fast atomistic simulations                                !
!  Copyright (C) 2006 - 2025  DFTB+ developers group                                               !
!                                                                                                  !
!  See the LICENSE file for terms of usage and distribution.                                       !
!--------------------------------------------------------------------------------------------------!

#:include 'linkedlist.fypp'

!> Linked list for single strings
module dftbp_type_linkedlistlc0
  use dftbp_common_accuracy, only : lc
  implicit none

  private

  $:define_list(&
      & TYPE_NAME='TListCharLc',&
      & ITEM_TYPE='character(lc)',&
      & PADDING='""')

end module dftbp_type_linkedlistlc0
