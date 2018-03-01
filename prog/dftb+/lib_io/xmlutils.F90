!--------------------------------------------------------------------------------------------------!
!  DFTB+: general package for performing fast atomistic simulations                                !
!  Copyright (C) 2018  DFTB+ developers group                                                      !
!                                                                                                  !
!  See the LICENSE file for terms of usage and distribution.                                       !
!--------------------------------------------------------------------------------------------------!

#:include 'common.fypp'

!> Utilities for processing an XML tree
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


  !> Returns first child with the specified name.
  function getFirstChildByName(node, name) result(child)

    !> Parent node containing children
    type(fnode), pointer :: node

    !> Child name to look for (empty string returns the node itself)
    character(len=*), intent(in) :: name

    !> Pointer to child with the specified name or null pointer if not found. If the name parameter
    !> was empty, a pointer to the node itself will be returned.
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


  !> Returns last child with the specified name.
  function getLastChildByName(node, name) result(child)

    !> Parent node containing children
    type(fnode), pointer :: node

    !> Child name to look for (empty string returns the node itself)
    character(len=*), intent(in) :: name

    !> Pointer to child with the specified name or null pointer if not found. If the name parameter
    !> was empty, a pointer to the node itself will be returned.
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


  !> Returns a list of children with the specified node name.
  function getChildrenByName(node, name) result(childList)

    !> Parent node to investigate.
    type(fnode), pointer :: node

    !> Name of the children to look for.
    character(len=*), intent(in) :: name

    !> List of the children on return.
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


  !> Remove text nodes with only whitespace characters from node and children.
  recursive subroutine removeSpace(node)

    !> Node to investigate
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


  !> Collects nodes in a tree that are without a specific attribute.
  function getTagsWithoutAttribute(node, name, rootOnly) result(nodeList)

    !> Tree to investigate
    type(fnode), pointer :: node

    !> Name of the attribute to look for
    character(len=*), intent(in) :: name

    !> Should children of a found attribute-less node be ignored?
    logical, intent(in), optional :: rootOnly

    !> List of the nodes without the specified attribute
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


  !> Recursive working subroutine for the getTagsWithoutAttribute routine
  recursive subroutine getTagsWithoutAttr_recursive(node, name, rootOnly, &
      &nodeList)

    !> Tree to investigate
    type(fnode), pointer :: node

    !> Name of the attribute to look for
    character(len=*), intent(in) :: name

    !> Should children of a found attribute-less node be ignored?
    logical, intent(in) :: rootOnly

    !> List of the nodes without the specified attribute
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


  !> Remove and destroy all children of a node.
  subroutine removeChildNodes(node)

    !> Node to process
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


  !> Removes nodes from a tree and destroys them.
  !> Caveat: The nodes must not be children of each other.
  subroutine removeNodes(nodeList)

    !> Contains the nodes to remove
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
