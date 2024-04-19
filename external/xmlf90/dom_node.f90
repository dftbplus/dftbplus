module xmlf90_dom_node

use xmlf90_dom_types
use xmlf90_dom_nodelist
use xmlf90_dom_namednodemap
use xmlf90_dom_debug
use xmlf90_dom_error

use xmlf90_strings

implicit none

private

  !-------------------------------------------------------   
  ! METHODS FOR NODES
  !-------------------------------------------------------   

  public :: getNodeName
  public :: getNodevalue
  public :: getNodeType
  public :: hasChildNodes
  public :: hasAttributes
  public :: getParentNode
  public :: getFirstChild
  public :: getLastChild
  public :: getNextSibling
  public :: getPreviousSibling
  public :: getOwnerDocument
  public :: getAttributes
  public :: getChildNodes
  public :: setNodeValue
  public :: appendChild
  public :: removeChild
  public :: replaceChild
  public :: cloneNode  
  public :: isSameNode
  public :: insertBefore

  private :: name_len, value_len

CONTAINS

  pure function name_len(node)
    type(fnode), pointer :: node
    integer :: name_len

    name_len = len_trim(node % nodeName)

  end function name_len

  pure function value_len(node)
    type(fnode), pointer :: node
    integer :: value_len

    value_len = len_trim(node % nodeValue)

  end function value_len

  !-----------------------------------------------------------
  !  METHODS FOR NODES
  !-----------------------------------------------------------
  subroutine getNodeName(node, nodeName)
    type(fnode), pointer :: node
    type(string), intent(inout)  :: nodeName

    if (.not. associated(node))  &
        call dom_error("getNodeName",0,"Node not allocated")
    nodeName = node%nodeName

  end subroutine getNodeName

  !-----------------------------------------------------------

  subroutine getNodeValue(node, nodeValue)
    type(fnode), pointer :: node
    type(string), intent(inout)  :: nodeValue

    if (.not. associated(node))  &
        call dom_error("getNodeValue",0,"Node not allocated")
    nodeValue = node%nodeValue
    
  end subroutine getNodeValue

  !-----------------------------------------------------------

  function getNodeType(node)

    type(fnode), pointer :: node
    integer :: getNodeType

    if (.not. associated(node)) call dom_error("getNodeType",0,"Node not allocated")
    getNodeType = node % nodeType

  end function getNodeTYpe

  !-----------------------------------------------------------

  function hasChildNodes(node)

    type(fnode), pointer :: node
    logical :: hasChildNodes

    if (.not. associated(node)) call dom_error("hasChildNodes",0,"Node not allocated")
    hasChildNodes = associated(node % firstChild)

  end function hasChildNodes

  !-----------------------------------------------------------

  function hasAttributes(node)

    type(fnode), pointer    :: node
    logical                 :: hasAttributes

    hasAttributes = .false.
    if (.not. associated(node)) call dom_error("hasAttributes",0,"Node not allocated")
    if (node % nodeType /= ELEMENT_NODE) RETURN
    if ( getLength(node%attributes) > 0) hasAttributes = .true.

  end function hasAttributes

  !-----------------------------------------------------------

  function getParentNode(node)

    type(fnode), pointer    :: node
    type(fnode), pointer    :: getParentNode

    if (.not. associated(node)) call dom_error("getParentNode",0,"Node not allocated")
    getParentNode => node % parentNode
    
  end function getParentNode
  
  !-----------------------------------------------------------

  function getFirstChild(node)

    type(fnode), pointer    :: node
    type(fnode), pointer    :: getFirstChild

    if (.not. associated(node)) call dom_error("getFirstChild",0,"Node not allocated")
    getFirstChild => node % firstChild

  end function getFirstChild

  !-----------------------------------------------------------

  function getLastChild(node)

    type(fnode), pointer :: node
    type(fnode), pointer    :: getLastChild

    if (.not. associated(node)) call dom_error("getLastChild",0,"Node not allocated")
    getLastChild => node % lastChild

  end function getLastChild

  !-----------------------------------------------------------

  function getNextSibling(node)

    type(fnode), pointer :: node
    type(fnode), pointer    :: getNextSibling

    if (.not. associated(node)) call dom_error("getNextSibling",0,"Node not allocated")
    getNextSibling => node % nextSibling

  end function getNextSibling

  !-----------------------------------------------------------

  function getPreviousSibling(node)

    type(fnode), pointer     :: node
    type(fnode), pointer    :: getPreviousSibling

    if (.not. associated(node)) call dom_error("getPreviousSibling",0,"Node not allocated")
    getPreviousSibling => node % previousSibling

  end function getPreviousSibling

  !-----------------------------------------------------------

  function getOwnerDocument(node)

    type(fnode), pointer    :: node
    type(fnode), pointer    :: getOwnerDocument

    if (.not. associated(node)) call dom_error("getOwnerDocument",0,"Node not allocated")
    getOwnerDocument => node % ownerDocument

  end function getOwnerDocument

  !----------------------------------------------------------- 

  function getChildNodes(node) result(nodelist)
    
    type(fnode), pointer        :: node
    type(fnodeList), pointer    :: nodelist      !!! NB nodeList

    type(fnode), pointer        :: np

    if (.not. associated(node)) call dom_error("getChildNodes",0,"Node not allocated")
    nodelist => null()
    np => node%firstChild
    do 
       if (.not. associated(np)) exit
       call append(nodelist,np)
       np => np%nextSibling
    enddo

  end function getChildNodes

  !----------------------------------------------------------- 

  function getAttributes(node)

    type(fnode), pointer         :: node
    type(fnamedNodeMap), pointer :: getAttributes       !!! NB namedNodeMap
    
    if (.not. associated(node))  &
        call dom_error("getAttributes",0,"Node not allocated")
    getAttributes => node % attributes

  end function getAttributes

  !----------------------------------------------------------- 

  subroutine setNodeValue(node, value)

    type(fnode), pointer :: node
    character(len=*), intent(in) :: value
    
    if (.not. associated(node))  &
               call dom_error("setNodeValue",0,"Node not allocated")

    select case(node % nodeType)

    case(ATTRIBUTE_NODE)
       node % nodeValue = trim(value)    !!AG: use just value ??

    case(COMMENT_NODE)
       node % nodeValue = value

    case(TEXT_NODE)
       node % nodeValue = value

    case(PROCESSING_INSTRUCTION_NODE)
       node % nodeValue = value

    case(CDATA_SECTION_NODE)
       node % nodeValue = value

    end select

  end subroutine setNodeValue

  !-----------------------------------------------------------
  
  function appendChild(node, newChild)
    type(fnode), pointer :: node
    type(fnode), pointer :: newChild
    type(fnode), pointer :: appendChild
    
    if (.not. associated(node))  & 
               call dom_error("appendChild",0,"Node not allocated")

    if ((node%nodeType /= ELEMENT_NODE) .and. &
        (node%nodeType /= DOCUMENT_NODE)) &
    call dom_error("appendChild",HIERARCHY_REQUEST_ERR, &
           "this node cannot have children")

    if (.not.(associated(node % firstChild))) then
       node % firstChild => newChild
    else 
       newChild % previousSibling   => node % lastChild
       node % lastChild % nextSibling => newChild 
    endif

    node % lastChild               => newChild
    newChild % parentNode          => node
    newChild % ownerDocument       => node % ownerDocument
    node%nc  = node%nc + 1

    appendChild => newChild

  end function appendChild

  !-----------------------------------------------------------
  
  function removeChild(node, oldChild)

    type(fnode), pointer :: removeChild
    type(fnode), pointer :: node
    type(fnode), pointer :: oldChild
    type(fnode), pointer :: np
    
    if (.not. associated(node)) call dom_error("removeChild",0,"Node not allocated")
    np => node % firstChild

    do while (associated(np))
       if (associated(np, oldChild)) then   ! Two argument form 
                                              !  of associated()
          if (associated(np,node%firstChild)) then
             node%firstChild => np%nextSibling
             if (associated(np % nextSibling)) then
                np%nextSibling % previousSibling => null()
             else
                node%lastChild => null()    ! there was just 1 node
             endif
          else if (associated(np,node%lastChild)) then
             ! one-node-only case covered above
             node%lastChild => np%previousSibling
             np%previousSibling%nextSibling => null()
          else
             np % previousSibling % nextSibling => np % nextSibling
             np % nextSibling % previousSibling => np % previousSibling
          endif
          node%nc = node%nc -1
          np % previousSibling => null()    ! Are these necessary?
          np % nextSibling => null()
          np % parentNode => null()
          removeChild => oldChild
          RETURN
       endif
       np => np % nextSibling
    enddo

    call dom_error("removeChild",NOT_FOUND_ERR,"oldChild not found")

  end function removeChild

 !-----------------------------------------------------------
  
  function replaceChild(node, newChild, oldChild)

    type(fnode), pointer :: replaceChild
    type(fnode), pointer :: node
    type(fnode), pointer :: newChild
    type(fnode), pointer :: oldChild

    type(fnode), pointer :: np
    
    if (.not. associated(node)) call dom_error("replaceChild",0,"Node not allocated")
    if ((node%nodeType /= ELEMENT_NODE) .and. &
        (node%nodeType /= DOCUMENT_NODE)) &
    call dom_error("replaceChild",HIERARCHY_REQUEST_ERR, &
           "this node cannot have children")

    np => node % firstChild

    do while (associated(np))    
       if (associated(np, oldChild)) then
          if (associated(np,node%firstChild)) then
             node%firstChild => newChild
             if (associated(np % nextSibling)) then
                oldChild%nextSibling % previousSibling => newChild
             else
                node%lastChild => newChild    ! there was just 1 node
             endif
          else if (associated(np,node%lastChild)) then
             ! one-node-only case covered above
             node%lastChild => newChild
             oldChild%previousSibling%nextSibling => newChild
          else
             oldChild % previousSibling % nextSibling => newChild
             oldChild % nextSibling % previousSibling => newChild
          endif

          newChild % parentNode      => oldChild % parentNode
          newChild % nextSibling     => oldChild % nextSibling
          newChild % previousSibling => oldChild % previousSibling
          replaceChild => oldChild
          np % previousSibling => null()    ! Are these necessary?
          np % nextSibling => null()
          np % parentNode => null()
          RETURN
       endif
       np => np % nextSibling
    enddo

    call dom_error("replaceChild",NOT_FOUND_ERR,"oldChild not found")

  end function replaceChild

  !-----------------------------------------------------------

  function cloneNode(node, deep)             
    type(fnode), pointer :: cloneNode
    type(fnode), pointer :: node

    logical, intent(in), optional :: deep
    logical           :: do_children

    type(fnode), pointer :: original
    type(fnode), pointer :: parent_clone
    
    if (.not. associated(node)) call dom_error("cloneNode",0,"Node not allocated")

    do_children = .false.
    if (present(deep)) then
       do_children = deep
    endif
    
    original => node             ! Keep node
    cloneNode => null()
    parent_clone => null()
    call recursive_clone(original, cloneNode)
    cloneNode%parentNode => null()     ! as per specs   , superfÃluous
 
  Contains

    recursive subroutine recursive_clone(original, cloneNode)
      type(fnode), pointer :: original        ! node to clone
      type(fnode), pointer :: cloneNode       ! new node

      type(fnode), pointer :: np, clone
      type(fnode), pointer :: previous_clone, attr, newattr
      type(string)         :: name
      logical :: first_sibling
      integer :: i

      np => original
      previous_clone => null()
      first_sibling = .true.
      do 

         ! Keep going across siblings
         ! (2nd and lower levels only)

         if (.not.(associated(np))) EXIT


         !----------------------------------------------------!
         clone => createNode()
         if (first_sibling) then
            cloneNode => clone       ! Rest of siblings are chained
                                     ! automatically, but must not
                                     ! be aliases of cloneNode !!
            first_sibling = .false.
         endif
         clone % nodeName    = np % nodeName
         name = np%nodeName
         if (dom_debug) print *, "Cloning ", char(name)
         clone % nodeValue   = np % nodeValue
         clone % nodeType    = np % nodeType
         clone % ownerDocument => np % ownerDocument
         clone % parentNode  => parent_clone
         !
         ! always deep copy attributes, as per specs
         ! Note that this will not work for "deep" attributes, with
         ! hanging entity nodes, etc
         if (associated(np % attributes)) then
            do i = 0, getLength(np%attributes) - 1
               attr => item(np%attributes,i)
               newattr => createNode()
               call getNodeName(attr, newattr%nodeName)
               call getNodeValue(attr, newattr%nodeValue)
               newattr%nodeType = ATTRIBUTE_NODE
               call append(clone%attributes, newattr)
            enddo
         endif

         ! Deal with first sibling
         if (associated(previous_clone)) then
            if (dom_debug) print *, "linking to previous sibling"
            previous_clone%nextSibling => clone
            clone%previousSibling => previous_clone
         else
            if (dom_debug) print *, "marking as first child of parent"
            if (associated(parent_clone))  &
                           parent_clone%firstChild => clone
         endif

         ! Deal with last sibling
         if (.not. associated(np%nextSibling)) then
            if (dom_debug) print *, "this is the last sibling"
            if (associated(parent_clone)) then
               if (dom_debug) print *, "marking as last child of parent"
               parent_clone%lastChild => clone
            endif
         endif
            
         if (do_children .and. associated(np%firstChild)) then
            parent_clone => clone
            if (dom_debug) print *, ".... going for its children"
            call recursive_clone(np%firstChild,clone%firstChild)
            parent_clone => clone%parentNode
         endif

         if (associated(np,node)) then
            if (dom_debug) print *, "No more siblings of ", char(name)
            EXIT  ! no siblings of main node
         endif
         np => np % nextSibling
         previous_clone => clone

      enddo

    end subroutine recursive_clone

  end function cloneNode

  !-----------------------------------------------------------

  function isSameNode(node1, node2)    ! DOM 3.0
    type(fnode), pointer :: node1
    type(fnode), pointer :: node2
    logical :: isSameNode

    isSameNode = associated(node1, node2)

  end function isSameNode

  !-----------------------------------------------------------

  function insertBefore(node, newChild, refChild)
    type(fnode), pointer :: insertBefore
    type(fnode), pointer :: node
    type(fnode), pointer :: newChild
    type(fnode), pointer :: refChild
    type(fnode), pointer :: np

    if (.not. associated(node)) call dom_error("insertBefore",0,"Node not allocated")
    if ((node%nodeType /= ELEMENT_NODE) .and. &
        (node%nodeType /= DOCUMENT_NODE)) &
    call dom_error("insertBefore",HIERARCHY_REQUEST_ERR, &
           "cannot insert node here")

    if (.not.associated(refChild)) then
       insertBefore => appendChild(node, newChild)
       RETURN
    endif

    np => node % firstChild
    do while (associated(np))
       if (associated(np, refChild)) then
          if (associated(np,node%firstChild)) then
             node%firstChild => newChild
          else
             refChild%previousSibling%nextSibling => newChild
          endif

          refChild % previousSibling => newChild
          newChild % nextSibling => refChild
          newChild % parentNode => node
          newChild % ownerDocument => refChild % ownerDocument
          insertBefore => newChild
          RETURN
       endif
       np => np % nextSibling
    enddo

    call dom_error("insertBefore",NOT_FOUND_ERR,"refChild not found")

  end function insertBefore

!----------------------------------------------------------------------

end module xmlf90_dom_node

