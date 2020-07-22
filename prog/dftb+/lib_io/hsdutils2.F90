!--------------------------------------------------------------------------------------------------!
!  DFTB+: general package for performing fast atomistic simulations                                !
!  Copyright (C) 2006 - 2020  DFTB+ developers group                                               !
!                                                                                                  !
!  See the LICENSE file for terms of usage and distribution.                                       !
!--------------------------------------------------------------------------------------------------!

#:include 'common.fypp'

!> Contains more high level functions for converting the values in a XML/HSD DOM-tree to Fortran
!> intrinsic types.
module dftbp_hsdutils2
  use dftbp_assert
  use dftbp_accuracy
  use dftbp_hsdutils
  use dftbp_hsdparser
  use dftbp_xmlutils
  use dftbp_unitconversion, only : unit
  use dftbp_message
  use dftbp_charmanip
  use dftbp_xmlf90
  implicit none
  private

  public :: getUnprocessedNodes, warnUnprocessedNodes,  getModifierIndex
  public :: readHSDAsXML
  public :: getNodeName2, setNodeName, removeModifier, splitModifier
  public :: setUnprocessed, getDescendant
  public :: convertByMul


  !> Converts according to passed modifier and array of possible units by multplicating the provided
  !> value with the appropriate conversion factor.
  interface convertByMul
    module procedure convertByMul_real
    module procedure convertByMul_realR1
    module procedure convertByMul_realR2
  end interface convertByMul


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


  !> Gets the index of a modifier from an array of possible modifier names.
  function getModifierIndex(modifier, modifiers, node, requested) result(ind)

    !> String containing the parsed modifier
    character(len=*), intent(in) :: modifier

    !> Array containing the names of the possible modifiers
    type(unit), intent(in) :: modifiers(:)

    !> Node for which the modifier was obtained (for errors)
    type(fnode), pointer :: node

    !> Should an error be raised, if the modifier is not found?
    logical, intent(in), optional :: requested

    !> Index of the modifer (zero if not found)
    integer :: ind

    character(len=len(modifier)) :: modifierLo
    logical :: mandatory
    integer :: ii

    if (present(requested)) then
      mandatory = requested
    else
      mandatory = .true.
    end if
    modifierLo = tolower(modifier)

    ind = 0
    do ii = 1, size(modifiers)
      if (trim(modifiers(ii)%name) == modifierLo) then
        ind = ii
        exit
      end if
    end do

    if (ind == 0 .and. mandatory) then
      call detailedError(node, MSG_INVALID_MODIFIER // modifier)
    end if

  end function getModifierIndex


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

    !>  Array of the modifiers, occuring in modifer.
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


  !> Implementation of convertByMul for real scalar.
  subroutine convertByMul_real(modifier, units, child, convertValue, replace, changed)

    !> Modifier (name of the unit to use)
    character(len=*), intent(in) :: modifier

    !> Array of the possible units
    type(unit), intent(in) :: units(:)

    !> The child, which carries the modifier.
    type(fnode), pointer :: child

    !> Value to convert, converted value on return.
    real(dp), intent(inout) :: convertValue

    !> If childs value should replaced by the new value (default: .false.)
    logical, intent(in), optional :: replace

    !> Contains flag on return, if childs value was changed.
    logical, intent(out), optional :: changed

    logical :: tReplace, tChanged
    integer :: ind

    if (present(replace)) then
      tReplace = replace
    else
      tReplace = .false.
    end if

    if (len(modifier) > 0) then
      tChanged = .true.
      ind = getModifierIndex(modifier, units, child)
      convertValue = convertValue * units(ind)%convertValue
      if (tReplace) then
        call setChildValue(child, "", convertValue, .true.)
      end if
    else
      tChanged = .false.
    end if

    if (present(changed)) then
      changed = tChanged
    end if

  end subroutine convertByMul_real


  !> Implementation of convertByMul for real rank one array.
  subroutine convertByMul_realR1(modifier, units, child, convertValue, replace, changed)

    !> Modifier (name of the unit to use)
    character(len=*), intent(in) :: modifier

    !> Array of the possible units
    type(unit), intent(in) :: units(:)

    !> The child, which carries the modifier.
    type(fnode), pointer :: child

    !> Value to convert, converted value on return.
    real(dp), intent(inout) :: convertValue(:)

    !> If childs value should replaced by the new value (default: .false.)
    logical, intent(in), optional :: replace

    !> Contains flag on return, if childs value was changed.
    logical, intent(out), optional :: changed

    logical :: tReplace, tChanged
    integer :: ind

    if (present(replace)) then
      tReplace = replace
    else
      tReplace = .false.
    end if

    if (len(modifier) > 0) then
      tChanged = .true.
      ind = getModifierIndex(modifier, units, child)
      convertValue = convertValue * units(ind)%convertValue
      if (tReplace) then
        call setChildValue(child, "", convertValue, .true.)
      end if
    else
      tChanged = .false.
    end if

    if (present(changed)) then
      changed = tChanged
    end if

  end subroutine convertByMul_realR1


  !> Implementation of convertByMul for real rank two array.
  subroutine convertByMul_realR2(modifier, units, child, convertValue, replace, changed)

    !> Modifier (name of the unit to use)
    character(len=*), intent(in) :: modifier

    !> Array of the possible units
    type(unit), intent(in) :: units(:)

    !> The child, which carries the modifier.
    type(fnode), pointer :: child

    !> Value to convert, converted value on return.
    real(dp), intent(inout) :: convertValue(:,:)

    !> If childs value should replaced by the new value (default: .false.)
    logical, intent(in), optional :: replace

    !> Contains flag on return, if childs value was changed.
    logical, intent(out), optional :: changed

    logical :: tReplace, tChanged
    integer :: ind

    if (present(replace)) then
      tReplace = replace
    else
      tReplace = .false.
    end if

    if (len(modifier) > 0) then
      tChanged = .true.
      ind = getModifierIndex(modifier, units, child)
      convertValue = convertValue * units(ind)%convertValue
      if (tReplace) then
        call setChildValue(child, "", convertValue, .true.)
      end if
    else
      tChanged = .false.
    end if

    if (present(changed)) then
      changed = tChanged
    end if

  end subroutine convertByMul_realR2


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

end module dftbp_hsdutils2
