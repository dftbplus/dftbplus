!--------------------------------------------------------------------------------------------------!
!  DFTB+: general package for performing fast atomistic simulations                                !
!  Copyright (C) 2018  DFTB+ developers group                                                      !
!                                                                                                  !
!  See the LICENSE file for terms of usage and distribution.                                       !
!--------------------------------------------------------------------------------------------------!

#:include 'linkedlist.fypp'

!> Linked list of single strings
module linkedlists0
  use assert
  use xmlf90
  implicit none
  private

  $:define_list(&
      & TYPE_NAME='listString',&
      & ITEM_TYPE='character(len=*)',&
      & NODE_TYPE='type(string)',&
      & PADDING="''", )

end module linkedlists0
