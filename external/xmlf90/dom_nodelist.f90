module xmlf90_dom_nodelist

use xmlf90_dom_types

implicit none

private

public :: item, getItem1
public :: getLength
public :: append

interface append
   module procedure append_nl
end interface

interface item
   module procedure item_nl
end interface

interface getItem1
  module procedure getitem1_nl
end interface

interface getLength
   module procedure getLength_nl
end interface

CONTAINS

  !-----------------------------------------------------------
  ! METHODS FOR NODELISTS
  !-----------------------------------------------------------
  function item_nl(nodeList, i)
    
    integer, intent(in)             :: i
    type(fnodeList), pointer        :: nodeList
    type(fnode), pointer            :: item_nl
    
    type(flistNode), pointer :: lp
    integer :: n

    item_nl => null()            ! In case there is no such item
    if (.not. associated(nodeList) .or. i >= nodeList%length) then
      return
    end if

    lp => nodeList%head
    n = -1
    do 
       if (.not. associated(lp))  exit
       n = n + 1
       if (n == i) then
          item_nl => lp%node
          exit
       endif
       lp => lp%next
    enddo

  end function item_nl


  !!* Subroutine version for item_nl with shifted indices and cache (B.A.)
  subroutine getitem1_nl(nodeList, i, item)
    type(fnodeList), pointer        :: nodeList
    integer, intent(in)             :: i
    type(fnode), pointer            :: item
    
    type(flistNode), pointer :: lp
    integer :: n
    
    item => null()            ! In case there is no such item
    if (.not. associated(nodeList) .or. i > nodeList%length) then
      return
    end if

    if (associated(nodeList%pCache) .and. nodeList%iCache <= i) then
      lp => nodeList%pCache
      n = nodeList%iCache
    else
      lp => nodeList%head
      n = 1
    end if
    do 
      if (.not. associated(lp)) then
        exit
      end if
      if (n == i) then
        item => lp%node
        nodeList%pCache => lp
        nodeList%iCache = n
        exit
      end if
      n = n + 1
      lp => lp%next
    end do

  end subroutine getitem1_nl


  !----------------------------------------------------------- 

  function getLength_nl(nodeList)
    
    type(fnodeList), pointer :: nodeList
    integer                  :: getLength_nl

    if (.not. associated(nodeList)) then
       getLength_nl = 0
    else
       getLength_nl = nodeList % length
    endif

  end function getLength_nl

  subroutine append_nl(nodeList,node)
    type(fnodeList), pointer :: nodeList
    type(fnode), pointer :: node

    if (.not. associated(nodeList)) then
       allocate(nodeList)
       nodelist%length = 1
       allocate(nodelist%head)
       nodelist%head%node => node
       nodelist%head%next => null()
       nodelist%tail => nodelist%head
    else
       allocate(nodelist%tail%next)
       nodelist%tail%next%node => node
       nodelist%tail%next%next => null()
       nodelist%tail => nodelist%tail%next
       nodelist%length = nodelist%length + 1
    endif

  end subroutine append_nl

end module xmlf90_dom_nodelist

