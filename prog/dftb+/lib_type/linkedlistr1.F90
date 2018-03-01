!--------------------------------------------------------------------------------------------------!
!  DFTB+: general package for performing fast atomistic simulations                                !
!  Copyright (C) 2018  DFTB+ developers group                                                      !
!                                                                                                  !
!  See the LICENSE file for terms of usage and distribution.                                       !
!--------------------------------------------------------------------------------------------------!

#:include 'linkedlist.fypp'

!> Linked list for real vectors
module linkedlistr1
  use accuracy, only : dp
  use assert
  implicit none
  private

  $:define_list(&
      & TYPE_NAME='listRealR1',&
      & ITEM_TYPE='real(dp)',&
      & ITEM_RANK=1,&
      & PADDING='0.0_dp')

end module linkedlistr1
