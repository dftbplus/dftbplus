!--------------------------------------------------------------------------------------------------!
!  DFTB+: general package for performing fast atomistic simulations                                !
!  Copyright (C) 2006 - 2021  DFTB+ developers group                                               !
!                                                                                                  !
!  See the LICENSE file for terms of usage and distribution.                                       !
!--------------------------------------------------------------------------------------------------!

#:include 'linkedlist.fypp'

!> Linked list of single strings
module dftbp_type_linkedlists0
  use dftbp_extlibs_xmlf90, only : string, assignment(=), operator(==), len
  implicit none

  private

  $:define_list(&
      & TYPE_NAME='TListString',&
      & ITEM_TYPE='character(len=*)',&
      & NODE_TYPE='type(string)',&
      & PADDING="''", )

end module dftbp_type_linkedlists0
