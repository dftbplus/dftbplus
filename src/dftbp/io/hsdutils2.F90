!--------------------------------------------------------------------------------------------------!
!  DFTB+: general package for performing fast atomistic simulations                                !
!  Copyright (C) 2006 - 2022  DFTB+ developers group                                               !
!                                                                                                  !
!  See the LICENSE file for terms of usage and distribution.                                       !
!--------------------------------------------------------------------------------------------------!

#:include 'common.fypp'

#:set UNIT_CONVERSION_RANKS = [0, 1, 2]

!> Contains more high level functions for converting the values in a XML/HSD DOM-tree to Fortran
!> intrinsic types.
module dftbp_io_hsdutils2
  use dftbp_common_accuracy, only : dp
  use dftbp_common_unitconversion, only : TUnit, unitConvStat => statusCodes, convertUnit
  use dftbp_extlibs_xmlf90, only : fnode, fnodeList, string, trim, len, assignment(=), parsefile,&
      & getLength, item, char, removeAttribute, getAttribute, setAttribute, setTagName,&
      & normalize, append_to_string, destroyNodeList, removeAttribute
  use dftbp_io_charmanip, only : newline, tolower, i2c
  use dftbp_io_hsdparser, only : attrName, attrModifier
  use dftbp_io_hsdutils, only : attrProcessed, getChild, setChildValue, detailedError,&
      & appendPathAndLine
  use dftbp_io_message, only : error, warning
  use dftbp_io_xmlutils, only : getTagsWithoutAttribute, removeNodes, removeSpace
  implicit none

  private
  public :: getUnprocessedNodes, warnUnprocessedNodes
  public :: readHSDAsXML
  public :: getNodeName2, setNodeName, removeModifier, splitModifier
  public :: setUnprocessed, getDescendant
  public :: convertUnitHsd


  !> Converts according to passed modifier and array of possible units by multplicating the provided
  !> value with the appropriate conversion factor.
  interface convertUnitHsd
    #:for RANK in UNIT_CONVERSION_RANKS
      module procedure convertUnitHsdR${RANK}$
    #:endfor
  end interface convertUnitHsd


  !> Separator for modifiers
  character, parameter :: sepModifier = ","


  !> common error message within this module
  character(len=*), parameter :: MSG_INVALID_MODIFIER = "Invalid modifier: "

contains


  !> Removes the processed flag from node
  subroutine setUnprocessed(node)

    !> The node to process
    type(fnode), pointer :: node

    if (associated(node)) then
      call removeAttribute(node, attrProcessed)
    end if

  end subroutine setUnprocessed


  !> Prints a warning message about unprocessed nodes

  subroutine getUnprocessedNodes(node, nodeList)

    !> Root element of the tree to investigate
    type(fnode), pointer :: node

    !> Containst the list of unprocessed nodes.
    type(fnodeList), pointer :: nodeList

    nodeList => getTagsWithoutAttribute(node, attrProcessed)

  end subroutine getUnprocessedNodes


  !> Prints a warning message about unprocessed nodes
  subroutine warnUnprocessedNodes(node, tIgnoreUnprocessed, nodeList)

    !> Root element of the tree to investigate
    type(fnode), pointer :: node

    !> if anything left after processing should be flagged
    logical, intent(in), optional :: tIgnoreUnprocessed

    !> list of left over nodes (if present)
    type(fnodeList), pointer, optional :: nodeList

    type(fnodeList), pointer :: list
    type(fnode), pointer :: child
    type(string) :: msg
    integer :: ii, ll
    logical :: tIgnoreUnprocessed0

    call getUnprocessedNodes(node, list)
    ll = getLength(list)
    if (ll > 0) then
      msg = "The following " // i2c(ll) &
          &// " node(s) have been ignored by the parser:" // newline
      do ii = 0, ll-1
        child => item(list, ii)
        call append_to_string(msg, "(" // i2c(ii+1) // ")")
        call appendPathAndLine(child, msg)
        call append_to_string(msg, newline)
      end do
      call warning(char(msg))
    end if
    if (present(nodeList)) then
      nodeList => list
    else
      call removeNodes(list)
      call destroyNodeList(list)
    end if

    if (present(tIgnoreUnprocessed)) then
      tIgnoreUnprocessed0 = tIgnoreUnprocessed
    else
      tIgnoreUnprocessed0 = .false.
    end if
    if (.not. tIgnoreUnprocessed0 .and. (ll > 0)) then
      call error("Code halting due to the presence of errors in dftb_in file.")
    end if

  end subroutine warnUnprocessedNodes


  !> Reads a HSD tree from an xml-file
  subroutine readHSDAsXML(fileName, fp)

    !> file to read
    character(len=*), intent(in) :: fileName

    !> to the tree
    type(fnode), pointer :: fp

    fp => parsefile(fileName)
    call removeSpace(fp)
    call normalize(fp)

  end subroutine readHSDAsXML


  !> Returns the name of a node, with empty string for unassociated nodes.
  subroutine getNodeName2(node, nodeName)

    !> Node to get the name from
    type(fnode), pointer :: node

    !> Contains the node name for an associated node or empty string for an unassociated one.
    type(string), intent(inout) :: nodeName

    if (.not. associated(node)) then
      nodeName = ""
    else
      nodeName = node%nodeName
    end if

  end subroutine getNodeName2


  !> Changes the name of a given node.
  !>
  !> Returns if node is not associated.
  subroutine setNodeName(node, name)

    !> Node to change.
    type(fnode), pointer :: node

    !> New name of the node.
    character(len=*), intent(in) :: name

    type(string) :: buffer

    @:ASSERT(associated(node))

    call getAttribute(node, attrName, buffer)
    if (len(buffer) > 0) then
      call removeAttribute(node, attrName)
      call setAttribute(node, attrName, name)
    end if
    call setTagName(node, tolower(name))

  end subroutine setNodeName


  !> Removes the modifier attribute from a given node.
  subroutine removeModifier(node)

    !> The node to process.
    type(fnode), pointer :: node

    if (associated(node)) then
      call removeAttribute(node, attrModifier)
    end if

  end subroutine removeModifier


  !> Splits a modifier containing coma separated list of modifiers into components.
  !>
  !> Note: if the number of the modifiers found differs from the size of the modifiers array, the
  !> program stops with error.
  subroutine splitModifier(modifier, child, modifiers)

    !> The list of modifers as a string.
    character(len=*), intent(in) :: modifier

    !>  The child which carries this modifier (for error messages)
    type(fnode), pointer :: child

    !>  Array of the modifiers, occurring in modifer.
    type(string), intent(inout) :: modifiers(:)

    integer :: nModif
    integer :: ii, iStart, iEnd

    nModif = size(modifiers)
    iStart = 1
    do ii = 1, nModif - 1
      iEnd = index(modifier(iStart:), sepModifier)
      if (iEnd == 0) then
        call detailedError(child, "Invalid number of specified modifiers (" &
            &// i2c(ii) // " instead of " // i2c(nModif) // ").")
      end if
      iEnd = iStart + iEnd - 1
      modifiers(ii) = trim(adjustl(modifier(iStart:iEnd-1)))
      iStart = iEnd + 1
    end do
    if (index(modifier(iStart:), sepModifier) /= 0) then
      call detailedError(child, "Invalid number of specified modifiers (&
          &more than " // i2c(nModif) // ").")
    end if
    modifiers(nModif) = trim(adjustl(modifier(iStart:)))

  end subroutine splitModifier


#:for RANK in UNIT_CONVERSION_RANKS

  !> Implementation of convertUnitHsd for given rank
  !>
  !> Stops by calling detailedError(), if modifier is not found among the passed units.
  !>
  subroutine convertUnitHsdR${RANK}$(modifier, units, child, convertValue, replace, changed)

    !> Modifier (name of the unit to use)
    character(len=*), intent(in) :: modifier

    !> Array of the possible units
    type(TUnit), intent(in) :: units(:)

    !> The child, which carries the modifier.
    type(fnode), pointer :: child

    !> Value to convert, converted value on return.
    real(dp), intent(inout) :: convertValue${FORTRAN_ARG_DIM_SUFFIX(RANK)}$

    !> If childs value should replaced by the new value (default: .false.)
    logical, intent(in), optional :: replace

    !> Contains flag on return, if childs value was changed.
    logical, intent(out), optional :: changed


    logical :: replace_, changed_
    integer :: status

    if (present(replace)) then
      replace_ = replace
    else
      replace_ = .false.
    end if

    changed_ = len(modifier) > 0
    if (changed_) then
      call convertUnit(units, modifier, convertValue, status)
      if (status /= unitConvStat%ok) then
        call detailedError(child, MSG_INVALID_MODIFIER // modifier)
      end if
      if (replace_) then
        call setChildValue(child, "", convertValue, .true.)
      end if
    end if

    if (present(changed)) changed = changed_

  end subroutine convertUnitHsdR${RANK}$

#:endfor



  !> Returns a descendant of a given node.
  subroutine getDescendant(root, path, child, requested, processed, parent)

    !> Node to seek the descendants of
    type(fnode), pointer :: root

    !> Path to the descendant. Parents are separated by "/" from their children
    !> (e.g. node1/node2/node3)
    character(len=*), intent(in) :: path

    !> Pointer to the child on return or null pointer if not found
    type(fnode), pointer :: child

    !> Should the program stop, if specified descendant is not present (default: .false.)
    logical, intent(in), optional :: requested

    !> Should elements along the path marked as processed? (default: .false.)
    logical, intent(in), optional :: processed

    !>If provided, contains parent node of the child, or the last associated node, if the child was
    !>not found.
    type(fnode), pointer, optional :: parent

    character(len=*), parameter :: pathSep = "/"

    logical :: tRequested, tUnprocessed
    type(fnode), pointer :: par
    integer :: iStart, iPos

    if (present(requested)) then
      tRequested = requested
    else
      tRequested = .false.
    end if
    if (present(processed)) then
      tUnprocessed = .not. processed
    else
      tUnprocessed = .true.
    end if

    iStart = 1
    par => null()
    child => root
    iPos = index(path, pathSep)
    do while (iPos /= 0 .and. associated(child))
      par => child
      call getChild(par, path(iStart:iStart+iPos-2), child, &
          &requested=tRequested)
      if (.not. associated(child)) then
        exit
      end if
      if (tUnprocessed) then
        call setUnprocessed(child)
      end if
      iStart = iStart + iPos
      iPos = index(path(iStart:), pathSep)
    end do
    if (associated(child)) then
      par => child
      call getChild(par, path(iStart:), child, requested=tRequested)
      if (associated(child) .and. tUnprocessed) then
        call setUnprocessed(child)
      end if
    end if

    if (present(parent)) then
      parent => par
    end if

  end subroutine getDescendant

end module dftbp_io_hsdutils2
