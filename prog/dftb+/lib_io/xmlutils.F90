!--------------------------------------------------------------------------------------------------!
!  DFTB+: general package for performing fast atomistic simulations                                !
!  Copyright (C) 2017  DFTB+ developers group                                                      !
!                                                                                                  !
!  See the LICENSE file for terms of usage and distribution.                                       !
!--------------------------------------------------------------------------------------------------!

#:include 'common.fypp'

module xmlutils
  use assert
  use charmanip
  use xmlf90
  implicit none

  private

  public :: getFirstChildByName, getLastChildByName, removeSpace
  public :: getTagsWithoutAttribute, removeChildNodes, removeNodes
  public :: getChildrenByName


contains

  !!* Returns first child with a certain name.
  !!* @param node Parent node containing children
  !!* @param name Child name to look for (empty string returns the node itself)
  !!* @return Pointer to child with the specified name or null pointer if not
  !!*   found. If the name parameter was empty, a pointer to the node itself
  !!*   will be returned.
  function getFirstChildByName(node, name) result(child)
    type(fnode), pointer :: node
    character(len=*), intent(in) :: name
    type(fnode), pointer :: child

    type(string) :: buffer

    @:ASSERT(associated(node))

    if (len(name) == 0) then
      child => node
    else
      child => getFirstChild(node)
      do while (associated(child))
        call getNodeName(child, buffer)
        if (buffer == name) then
          exit
        end if
        child => getNextSibling(child)
      end do
    end if

  end function getFirstChildByName



  !!* Returns last child with a certain name.
  !!* @param node Parent node containing children
  !!* @param name Child name to look for (empty string returns the node itself)
  !!* @return Pointer to child with the specified name or null pointer if not
  !!*   found. If the name parameter was empty, a pointer to the node itself
  !!*   will be returned.
  function getLastChildByName(node, name) result(child)
    type(fnode), pointer :: node
    character(len=*), intent(in) :: name
    type(fnode), pointer :: child

    type(string) :: buffer

    @:ASSERT(associated(node))

    if (len(name) == 0) then
      child => node
    else
      child => getLastChild(node)
      do while (associated(child))
        call getNodeName(child, buffer)
        if (buffer == name) then
          exit
        end if
        child => getPreviousSibling(child)
      end do
    end if

  end function getLastChildByName



  !!* Returns a list of children with a specific node name.
  !!* @param node      Parent node to investigate.
  !!* @param name      Name of the children to look for.
  !!* @param childList List of the children on return.
  function getChildrenByName(node, name) result(childList)
    type(fnode), pointer :: node
    character(len=*), intent(in) :: name
    type(fnodeList), pointer :: childList

    type(fnode), pointer :: child
    type(string) :: buffer

    @:ASSERT(associated(node))

    nullify(childList)
    if (len(name) == 0) then
      return
    end if

    child => getFirstChild(node)
    do while (associated(child))
      call getNodeName(child, buffer)
      if (buffer == name) then
        call append(childList, child)
      end if
      child => getNextSibling(child)
    end do

  end function getChildrenByName




  !!* Remove text nodes with only whitespace characters from node and children.
  !!* @param node Node to investigate
  recursive subroutine removeSpace(node)
    type(fnode), pointer :: node

    type(fnode), pointer :: child, child2, dummy
    type(string) :: buffer

    @:ASSERT(associated(node))

    child => getFirstChild(node)
    do while (associated(child))
      child2 => getNextSibling(child)
      if (getNodeType(child) == TEXT_NODE) then
        call getNodeValue(child, buffer)
        if (len_trim2(char(buffer)) == 0) then
          dummy => removeChild(node, child)
          call destroyNode(child)
        end if
      else
        call removeSpace(child)
      end if
      child => child2
    end do

  end subroutine removeSpace



  !!* Collects nodes in a tree without a specific attribute.
  !!* @param node     Tree to investigate
  !!* @param name     Name of the attribute to look for
  !!* @param rootOnly Should children of a found attributeless node be ignored?
  !!* @return List of the nodes without the specified attribute
  function getTagsWithoutAttribute(node, name, rootOnly) result(nodeList)
    type(fnode), pointer :: node
    character(len=*), intent(in) :: name
    logical, intent(in), optional :: rootOnly
    type(fnodeList), pointer :: nodeList

    logical :: tRootOnly

    @:ASSERT(associated(node))
    @:ASSERT(len(name) > 0)

    if (present(rootOnly)) then
      tRootOnly = rootOnly
    else
      tRootOnly = .true.
    end if
    nodeList => null()
    call getTagsWithoutAttr_recursive(node, name, tRootOnly, nodeList)

  end function getTagsWithoutAttribute



  !!* Recursive working subroutine for the getTagsWithoutAttribute routine
  !!* @param node     Tree to investigate
  !!* @param name     Name of the attribute to look for
  !!* @param rootOnly Should children of a found attributeless node be ignored?
  !!* @param nodeList List of the nodes without the specified attribute
  recursive subroutine getTagsWithoutAttr_recursive(node, name, rootOnly, &
      &nodeList)
    type(fnode), pointer :: node
    character(len=*), intent(in) :: name
    logical, intent(in) :: rootOnly
    type(fnodeList), pointer :: nodeList

    type(fnode), pointer :: attr, child

    @:ASSERT(associated(node))
    @:ASSERT(len(name) > 0)

    attr => getAttributeNode(node, name)
    if (.not. associated(attr)) then
      call append(nodeList, node)
      if (rootOnly) then
        return
      end if
    end if
    child => getFirstChild(node)
    do while (associated(child))
      if (getNodeType(child) /= TEXT_NODE) then
        call getTagsWithoutAttr_recursive(child, name, rootOnly, nodeList)
      end if
      child => getNextSibling(child)
    end do

  end subroutine getTagsWithoutAttr_recursive



  !!* Remove and destroy all children of a node.
  !!* @param node Node to process
  subroutine removeChildNodes(node)
    type(fnode), pointer :: node

    type(fnode), pointer :: child, child2

    @:ASSERT(associated(node))

    child => getFirstChild(node)
    do while (associated(child))
      child2 => getNextSibling(child)
      child => removeChild(node, child)
      call destroyNode(child)
      child => child2
    end do

  end subroutine removeChildNodes



  !!* Removes nodes from a tree and destroys them.
  !!* @param nodeList Contains the nodes to remove
  !!* @caveat The nodes must not be each others children.
  subroutine removeNodes(nodeList)
    type(fnodeList), pointer :: nodeList

    type(fnode), pointer :: child
    integer :: ii

    do ii = 0, getLength(nodeList) - 1
      child => item(nodeList, ii)
      child => removeChild(getParentNode(child), child)
      call destroyNode(child)
    end do

  end subroutine removeNodes



end module xmlutils
