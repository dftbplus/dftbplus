!--------------------------------------------------------------------------------------------------!
!  DFTB+: general package for performing fast atomistic simulations                                !
!  Copyright (C) 2018  DFTB+ developers group                                                      !
!                                                                                                  !
!  See the LICENSE file for terms of usage and distribution.                                       !
!--------------------------------------------------------------------------------------------------!

#:include 'linkedlist.fypp'

!> Linked list for single real values
module linkedlistr0
  use accuracy, only : dp
  use assert
  implicit none
  private

  $:define_list(&
      & TYPE_NAME='listReal',&
      & ITEM_TYPE='real(dp)',&
      & PADDING='0.0_dp')

end module linkedlistr0
