!--------------------------------------------------------------------------------------------------!
!  DFTB+: general package for performing fast atomistic simulations                                !
!  Copyright (C) 2006 - 2024  DFTB+ developers group                                               !
!                                                                                                  !
!  See the LICENSE file for terms of usage and distribution.                                       !
!--------------------------------------------------------------------------------------------------!

#:include 'linkedlist.fypp'

!> Linked list for real vectors
module dftbp_type_linkedlistr1
  use dftbp_common_accuracy, only : dp
  implicit none

  private

  $:define_list(&
      & TYPE_NAME='TListRealR1',&
      & ITEM_TYPE='real(dp)',&
      & ITEM_RANK=1,&
      & PADDING='0.0_dp')

end module dftbp_type_linkedlistr1
