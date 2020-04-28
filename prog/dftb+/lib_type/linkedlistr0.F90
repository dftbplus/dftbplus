!--------------------------------------------------------------------------------------------------!
!  DFTB+: general package for performing fast atomistic simulations                                !
!  Copyright (C) 2006 - 2020  DFTB+ developers group                                               !
!                                                                                                  !
!  See the LICENSE file for terms of usage and distribution.                                       !
!--------------------------------------------------------------------------------------------------!

#:include 'linkedlist.fypp'

!> Linked list for single real values
module dftbp_linkedlistr0
  use dftbp_accuracy, only : dp
  use dftbp_assert
  implicit none
  private

  $:define_list(&
      & TYPE_NAME='TListReal',&
      & ITEM_TYPE='real(dp)',&
      & PADDING='0.0_dp')

end module dftbp_linkedlistr0
