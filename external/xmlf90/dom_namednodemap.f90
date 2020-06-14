module xmlf90_dom_namednodemap
!
! This is basically a dictionary module, but written with the
! DOM node structure in mind.
!
use xmlf90_dom_types
use xmlf90_strings

implicit none

private
  !-------------------------------------------------------  
  ! METHODS FOR NAMEDNODEMAPS
  !-------------------------------------------------------   
  public :: getNamedItem
  public :: setNamedItem
  public :: removeNamedItem

  public :: item
  public :: getLength
  public :: append

  interface append
     module procedure append_nnm
  end interface

  interface item
     module procedure item_nnm
  end interface 

  interface getLength
     module procedure getLength_nnm
  end interface

CONTAINS

  function item_nnm(namedNodeMap, i)
    
    integer, intent(in)             :: i
    type(fnamedNodeMap), pointer    :: namedNodeMap
    type(fnode), pointer            :: item_nnm
    
    type(fnamedNode), pointer :: nnp

    integer :: n

    item_nnm => null()            ! In case there is no such item
    if (.not. associated(namedNodeMap)) RETURN

    nnp => namedNodeMap%head
    n = -1
    do 
       if (.not. associated(nnp))  exit
       n = n + 1
       if (n == i) then
          item_nnm => nnp%node
          exit
       endif
       nnp => nnp%next
    enddo

  end function item_nnm

  !----------------------------------------------------------- 
  
  function getLength_nnm(namedNodeMap)
  
    type(fnamedNodeMap), pointer :: namedNodeMap
    integer :: getLength_nnm

    getLength_nnm = 0
    if (.not. associated(namedNodeMap)) return

    getLength_nnm = namedNodeMap % length    
    
  end function getLength_nnm

  !----------------------------------------------------------- 


  subroutine append_nnm(nodeMap,node)
    type(fnamednodeMap), pointer :: nodeMap
    type(fnode), pointer :: node

    if (.not. associated(nodeMap)) then
       allocate(nodeMap)
       nodeMap%length = 1
       allocate(nodeMap%head)
       nodeMap%head%name = node%nodeName
       nodeMap%head%node => node
       nodeMap%tail => nodeMap%head
    else
      allocate(nodeMap%tail%next)
       nodeMap%tail%next%node => node
       nodeMap%tail%next%name =  node%nodeName
       nodeMap%tail => nodeMap%tail%next
       nodeMap%length = nodeMap%length + 1
    endif

  end subroutine append_nnm

  !----------------------------------------------------------- 

  function getNamedItem(namedNodeMap, name)
    
    type(fnamedNodeMap), pointer    :: namedNodeMap
    character(len=*), intent(in)    :: name
    type(fnode), pointer            :: getNamedItem

    type(fnamedNode), pointer :: nnp

    getNamedItem => null()
    if (.not. associated(namedNodeMap)) return 

    nnp => namedNodeMap%head
    do while (associated(nnp)) 
       if (nnp%name == name) then
          getNamedItem => nnp%node
          exit                 ! one or zero nodes with a given name
       endif
       nnp => nnp%next
    enddo

  end function getNamedItem

  
  function setNamedItem(namedNodeMap, node)

!!AG: Do we need to clone the node ?
    
    type(fnamedNodeMap), pointer    :: namedNodeMap
    type(fnode), pointer            :: node
    type(fnode), pointer            :: setNamedItem

    type(fnamedNode), pointer :: nnp

    if (.not. associated(namedNodeMap)) then

       call append(namedNodeMap,node)
       setNamedItem => node
      
    else

       nnp => namedNodeMap%head
       do while (associated(nnp)) 
          if (nnp%name == node%nodeName) then
            !setNamedItem => nnp%node
            call destroyNode(nnp%node)
             nnp%node => node
             setNamedItem => node
             return
          endif
          nnp => nnp%next
       enddo

       !   If not found, insert it at the end of the linked list

       call append(namedNodeMap,node)
       setNamedItem => node
    endif

  end function setNamedItem

!------------------------------------------------------------
   function removeNamedItem(namedNodeMap, name)
    
    type(fnamedNodeMap), pointer   :: namedNodeMap
    character(len=*), intent(in)   :: name
    type(fnode), pointer           :: removeNamedItem

    type(fnamedNode), pointer :: nnp, previous

    removeNamedItem => null()
    if (.not. associated(namedNodeMap)) return  

    previous => null()
    nnp => namedNodeMap%head
    do while (associated(nnp)) 
       if (nnp%name == name) then
          removeNamedItem => nnp%node
          if (associated(nnp,namedNodeMap%head)) then
             ! we remove the first fnamedNode in the chain...
             namedNodeMap%head => nnp%next
          else if (.not. associated(nnp%next)) then
             ! we remove the last fnamedNode in the chain
             previous%next => null()
             namedNodeMap%tail => previous
          else
             ! we remove a link in the middle of the chain
             previous%next => nnp%next
          endif
          namedNodeMap%length =  namedNodeMap%length - 1
          deallocate(nnp)
          EXIT                 ! one or zero nodes with a given name
       endif
       previous => nnp
       nnp => nnp%next
    enddo
    !! Deallocate the entire dictionary, if it is empty
    if (namedNodeMap%length == 0) then
      deallocate(namedNodeMap)
    end if

  end function removeNamedItem


end module xmlf90_dom_namednodemap
