!--------------------------------------------------------------------------------------------------!
!  DFTB+: general package for performing fast atomistic simulations                                !
!  Copyright (C) 2006 - 2025  DFTB+ developers group                                               !
!                                                                                                  !
!  See the LICENSE file for terms of usage and distribution.                                       !
!--------------------------------------------------------------------------------------------------!

#:include 'linkedlist.fypp'

!> Linked list for real arrays
module dftbp_type_linkedlistr2
  use dftbp_common_accuracy, only : dp
  implicit none

  private

  $:define_list(&
      & TYPE_NAME='TListRealR2',&
      & ITEM_TYPE='real(dp)',&
      & ITEM_RANK=2,&
      & PADDING='0.0_dp')

end module dftbp_type_linkedlistr2
