!--------------------------------------------------------------------------------------------------!
!  DFTB+: general package for performing fast atomistic simulations                                !
!  Copyright (C) 2006 - 2020  DFTB+ developers group                                               !
!                                                                                                  !
!  See the LICENSE file for terms of usage and distribution.                                       !
!--------------------------------------------------------------------------------------------------!

#:include 'linkedlist.fypp'

!> Linked list for single charcter strings
module dftbp_linkedlistmc0
  use dftbp_accuracy, only : mc
  use dftbp_assert
  implicit none
  private

  $:define_list(&
      & TYPE_NAME='TListCharMc',&
      & ITEM_TYPE='character(mc)',&
      & PADDING='""')

end module dftbp_linkedlistmc0
