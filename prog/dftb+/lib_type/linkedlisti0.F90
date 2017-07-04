!--------------------------------------------------------------------------------------------------!
!  DFTB+: general package for performing fast atomistic simulations                                !
!  Copyright (C) 2017  DFTB+ developers group                                                      !
!                                                                                                  !
!  See the LICENSE file for terms of usage and distribution.                                       !
!--------------------------------------------------------------------------------------------------!

#:include 'linkedlist.fypp'

module linkedlisti0
  implicit none
  private

  $:define_list(&
      & TYPE_NAME='listInt',&
      & ITEM_TYPE='integer',&
      & PADDING='0')

end module linkedlisti0
