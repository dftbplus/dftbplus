module xmlf90_dom_element

use xmlf90_dom_types
use xmlf90_dom_namednodemap
use xmlf90_dom_nodelist
use xmlf90_dom_attribute
use xmlf90_dom_document
use xmlf90_dom_debug
use xmlf90_dom_node
use xmlf90_strings
implicit none

private

  !-------------------------------------------------------   
  ! METHODS FOR ELEMENT NODES
  !-------------------------------------------------------   
  public :: getTagName
  public :: getElementsByTagName
  public :: getAttribute
  public :: getAttributeNode
  public :: setAttribute
  public :: setAttributeNode
  public :: removeAttribute
  public :: normalize       !--- combines adjacent text nodes ---!
  public :: setTagName

CONTAINS

  !-----------------------------------------------------------
  !  METHODS FOR ELEMENT NODES
  !-----------------------------------------------------------
  subroutine getTagName(element, tagName)
    type(fnode), intent(in) :: element   
    type(string), intent(inout) :: tagName

    if (element % nodeType == ELEMENT_NODE) then
      tagName = element % nodeName 
    else
      tagName = ''
    endif

  end subroutine getTagName

  !-----------------------------------------------------------
  function getElementsByTagName(element, tag) result(nodelist)
    type(fnode), pointer         :: element
    character(len=*), intent(in) :: tag
    type(fnodeList), pointer     :: nodelist 

    type(fnode), pointer        :: np

    nodelist => null()

    np => element
    if (dom_debug) print *, "Going into search for tag: ", trim(tag)
    call search(np)

    CONTAINS

    recursive subroutine search(np)
    type(fnode), pointer        :: np

    type(string)                :: name

    !
    ! Could replace the calls to helper methods by direct lookups of node 
    ! components to make it faster.
    ! 
    do
       if (.not. associated(np)) exit
       select case(np%nodeType)

          case(DOCUMENT_NODE) 
             ! special case ... search its children 
             if (hasChildNodes(np)) call search(getFirstChild(np))
             ! will exit for lack of siblings
          case(ELEMENT_NODE)

            call getNodeName(np, name)
             if (dom_debug) print *, "exploring node: ", char(name)
             if ((tag == "*") .or. (tag == name)) then
                call append(nodelist,np)
                if (dom_debug) print *, "found match ", nodelist%length
             endif
             if (hasChildNodes(np)) call search(getFirstChild(np))

          case default
             
             ! do nothing

        end select

        if (associated(np,element)) exit  ! no siblings of element...
        np => getNextSibling(np)

     enddo

    end subroutine search

  end function getElementsByTagName


  !-----------------------------------------------------------

  subroutine getAttribute(element, name, attribute)
    type(fnode), intent(in) :: element
    character(len=*), intent(in) :: name
    type(string), intent(inout) :: attribute

    type(fnode), pointer :: nn

    attribute = ""  ! as per specs, if not found
    if (element % nodeType /= ELEMENT_NODE) RETURN
    nn => getNamedItem(element%attributes,name)
    if (.not. associated(nn)) RETURN
    attribute = nn%nodeValue
        
  end subroutine getAttribute
  

  !-----------------------------------------------------------

  function getAttributeNode(element, name)
    
    type(fnode), intent(in) :: element
    type(fnode), pointer    :: getAttributeNode
    character(len=*), intent(in) :: name

    getAttributeNode => null()     ! as per specs, if not found
    if (element % nodeType /= ELEMENT_NODE) RETURN
    getAttributeNode => getNamedItem(element%attributes,name)

  end function getAttributeNode
  
  !-----------------------------------------------------------

  subroutine setAttributeNode(element, newattr)
    type(fnode), pointer :: element
    type(fnode), pointer :: newattr

    type(fnode), pointer :: dummy

    if (element % nodeType /= ELEMENT_NODE) then
       if (dom_debug) print *, "not an element node in setAttributeNode..."
       RETURN
    endif

    dummy => setNamedItem(element%attributes,newattr)
     
  end subroutine setAttributeNode

!-------------------------------------------------------------------
  subroutine setAttribute(element, name, value)
    type(fnode), pointer :: element
    character(len=*), intent(in) :: name
    character(len=*), intent(in) :: value

    type(fnode), pointer      :: newattr

    newattr => createAttribute(name)
    call setValue(newattr,value)
    call setAttributeNode(element,newattr)

  end subroutine setAttribute

  !-----------------------------------------------------------

  subroutine removeAttribute(element, name)
    type(fnode), pointer :: element
    character(len=*), intent(in) :: name

    type(fnode), pointer :: dummy

    if (element % nodeType /= ELEMENT_NODE) RETURN
    if (.not. associated(element%attributes)) RETURN

    dummy => removeNamedItem(element%attributes,name)
     
  end subroutine removeAttribute


  !-----------------------------------------------------------
  recursive subroutine normalize(element)
    type(fnode), pointer         :: element

    type(fnode), pointer        :: np
    type(fnode), pointer        :: head
    integer :: length, nTextNode


    if (dom_debug) print *, "Normalizing: ", trim(element%nodeName)
    np => element%firstChild
    nTextNode = 0
    do while (associated(np))

      select case(np%nodeType)

      case(TEXT_NODE)

        if (nTextNode == 0) then
          if (dom_debug) print *, "normalize: found first in chain"
          head => np
          length = len(np%nodeValue)
          nTextNode = 1
        else
          if (dom_debug) print *, "normalize: found second in chain"
          nTextNode = nTextNode + 1
          length = length + len(np%nodeValue)
        endif

      case(ELEMENT_NODE)
        if (nTextNode > 1) then
          call concatenateText(head, nTextNode, length)
          nTextNode = 0
        end if
        if (dom_debug) print *, "element sibling: ", trim(np%nodeName)
        if (hasChildNodes(np)) call normalize(np)

      case default
        ! do nothing, just mark that we break the chain of text nodes
        if (dom_debug) print *, "other sibling: ", trim(np%nodeName)
        if (nTextNode > 1) then
          call concatenateText(head, nTextNode, length)
          nTextNode = 0
        end if

      end select

      np => getNextSibling(np)

    end do
    if (nTextNode > 1) then
      call concatenateText(head, nTextNode, length)
    end if

  end subroutine normalize



  !!* Concatenates text nodes.
  !!* @param head      The first the node
  !!* @param nTextNode Number of all text nodes (must be greater as 1)
  !!* @param length    Total length of the text contained in the text nodes
  subroutine concatenateText(head, nTextNode, length)
    type(fnode), pointer :: head
    integer, intent(in) :: nTextNode
    integer, intent(in) :: length

    type(fnode), pointer :: np, ghost
    type(string) :: buffer
    integer :: ii

    call resize_string(buffer, length)
    call append_to_string(buffer, head%nodeValue)
    np => getNextSibling(head)
    do ii = 2, nTextNode - 1
      call append_to_string(buffer, np%nodeValue)
      ghost => np
      np => getNextSibling(np)
      call destroyNode(ghost)
    end do
    call append_to_string(buffer, np%nodeValue)
    head%nextSibling => np%nextSibling
    if (associated(np,np%parentNode%lastChild)) then
      np%parentNode%lastChild => head
      head%nextSibling => null()
    else
      np%nextSibling%previousSibling => head
    endif
    call destroyNode(np)
    head%nodeValue = buffer
    
  end subroutine concatenateText

  

  !!* Sets a new name for a given tag.
  !!* @param element Tag to rename.
  !!* @param tagName New name
  subroutine setTagName(element, tagName)
    type(fnode), intent(inout) :: element   
    character(len=*), intent(in) :: tagName
    
    if (element%nodeType == ELEMENT_NODE) then
      element%nodeName = tagName
    end if
    
  end subroutine setTagName




end module xmlf90_dom_element
