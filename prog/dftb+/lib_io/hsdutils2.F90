!--------------------------------------------------------------------------------------------------!
!  DFTB+: general package for performing fast atomistic simulations                                !
!  Copyright (C) 2017  DFTB+ developers group                                                      !
!                                                                                                  !
!  See the LICENSE file for terms of usage and distribution.                                       !
!--------------------------------------------------------------------------------------------------!

#:include 'common.fypp'

!!* Contains more high level functions for converting the values in a XML/HSD
!!* DOM-tree to Fortran intrinsic types.
module hsdutils2
  use assert
  use accuracy
  use hsdutils
  use hsdparser
  use xmlutils
  use unitconversion, only : unit
  use message
  use charmanip
  use xmlf90
  implicit none
  private

  public :: getUnprocessedNodes, warnUnprocessedNodes,  getModifierIndex
  public :: readHSDAsXML, readHSDOrXML
  public :: getNodeName2, setNodeName, removeModifier, splitModifier
  public :: setUnprocessed, getDescendant
  public :: convertByMul

  !!* Converts according to passed modifier and array of possible units by
  !!* multplicating the provided value with the appropriate conversion factor.
  !!* @param modifier Modifier (name of the unit to use)
  !!* @param units Array of the possible units
  !!* @param child The child, which carries the modifier.
  !!* @param value Value to convert, converted value on return.
  !!* @param replace If childs value should replaced by the new value
  !!*   (default: .false.)
  !!* @param changed Contains flag on return, if childs value was changed.
  interface convertByMul
    module procedure convertByMul_real
    module procedure convertByMul_realR1
    module procedure convertByMul_realR2
  end interface


  !!* Separator for modifiers
  character, parameter :: sepModifier = ","


  character(len=*), parameter :: MSG_INVALID_MODIFIER = "Invalid modifier: "



contains


  !!* Removes the processed flag from node
  !!* @param node The node to process
  subroutine setUnprocessed(node)
    type(fnode), pointer :: node

    if (associated(node)) then
      call removeAttribute(node, attrProcessed)
    end if

  end subroutine setUnprocessed



  !!* Gets the index of a modifier from an array of possible modifier names.
  !!* @param modifier  String containing the parsed modifier
  !!* @param modifiers Array containing the names of the possible modifiers
  !!* @param node      Node for which the modifier was obtained (for errors)
  !!* @param requested Should an error be raised, if the modifier is not found?
  !!* @return Index of the modifer (zero if not found)
  function getModifierIndex(modifier, modifiers, node, requested) result(ind)
    character(len=*), intent(in) :: modifier
    type(unit), intent(in) :: modifiers(:)
    type(fnode), pointer :: node
    logical, intent(in), optional :: requested
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



  !!* Prints a warning message about unprocessed nodes
  !!* @param node      Root element of the tree to investigate
  !!* @param nodesList Containst the list of unprocessed nodes.
  subroutine getUnprocessedNodes(node, nodeList)
    type(fnode), pointer :: node
    type(fnodeList), pointer :: nodeList

    nodeList => getTagsWithoutAttribute(node, attrProcessed)

  end subroutine getUnprocessedNodes



  !!* Prints a warning message about unprocessed nodes
  !!* @param node Root element of the tree to investigate
  subroutine warnUnprocessedNodes(node, tIgnoreUnprocessed, nodeList)
    type(fnode), pointer               :: node
    logical, intent(in), optional      :: tIgnoreUnprocessed
    type(fnodeList), pointer, optional :: nodeList


    type(fnodeList), pointer :: list
    type(fnode), pointer     :: child
    type(string)             :: msg
    integer                  :: ii, ll
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
    if (.not. tIgnoreUnprocessed .and. (ll > 0)) then
      call error("Code halting due to the presence of errors in dftb_in file.")
    end if

  end subroutine warnUnprocessedNodes



  !!* Reads a HSD tree from an xml-file
  !!* @param fileName file to read
  !!* @return Pointer to the tree
  subroutine readHSDAsXML(fileName, fp)
    character(len=*), intent(in) :: fileName
    type(fnode), pointer :: fp

    fp => parsefile(fileName)
    call removeSpace(fp)
    call normalize(fp)

  end subroutine readHSDAsXML



  !!* Reads HSD from an HSD-file or from an xml-file, but stop if input is
  !!* ambiguous or missing.
  !!* @param hsdFile  File containing the HSD input
  !!* @param xmlFile  File containing the XML input
  !!* @param rootTag  Name of the root tag
  !!* @param fp       Parsed tree on retunr
  !!* @param hsdInput True, if tree was read from an HSD file
  !!* @param missing  True, if neither xml nor hsd input was found
  !!* @note If optional parameter missing is not passed, the subroutine
  !!*   stops if input is not found. If ambigous input is found, the subroutine
  !!*   stops anyway.
  subroutine readHSDOrXML(hsdFile, xmlFile, rootTag, fp, hsdInput, missing)
    character(len=*), intent(in) :: hsdFile
    character(len=*), intent(in) :: xmlFile
    character(len=*), intent(in) :: rootTag
    type(fnode), pointer :: fp
    logical, intent(out), optional :: hsdInput
    logical, intent(out), optional :: missing

    logical :: tExist
    logical :: tHSD, tXML

    if (present(missing)) then
      missing = .false.
    end if
    inquire(file=hsdFile, exist=tExist)
    if (tExist) then
      tHSD = .true.
    else
      tHSD = .false.
    end if
    inquire(file=xmlFile, exist=tExist)
    if (tExist) then
      tXML = .true.
    else
      tXML = .false.
    end if
    if (tXML .and. tHSD) then
      call error("Ambiguous input: '" // hsdFile // "' and '" // xmlFile &
          &// "' both present")
    elseif (tXML) then
      call readHSDAsXML(xmlFile, fp)
    elseif (tHSD) then
      call parseHSD(rootTag, hsdFile, fp)
    elseif (present(missing)) then
      missing = .true.
    else
      call error("Missing input: nor '" // hsdFile // "' nor '" // xmlFile &
          &// "' present")
    end if
    if (present(hsdInput)) then
      hsdInput = tHSD
    end if

  end subroutine readHSDOrXML



  !!* Returns the name of a node, with empty string for unassociated nodes.
  !!* @param node      Node to get the name from
  !!* @param nodeName  Contains the node name for an associated node or empty
  !!*                  string for an unassociated one.
  subroutine getNodeName2(node, nodeName)
    type(fnode), pointer :: node
    type(string), intent(inout)  :: nodeName

    if (.not. associated(node)) then
      nodeName = ""
    else
      nodeName = node%nodeName
    end if

  end subroutine getNodeName2



  !!* Changes the name of a given node.
  !!* @param node Node to change.
  !!* @param name New name of the node.
  !!* @note Returns if node is not associated.
  subroutine setNodeName(node, name)
    type(fnode), pointer :: node
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



  !!* Removes the modifier attribute from a given node.
  !!* @param node The node to process.
  subroutine removeModifier(node)
    type(fnode), pointer :: node

    if (associated(node)) then
      call removeAttribute(node, attrModifier)
    end if

  end subroutine removeModifier



  !!* Splits a modifier containing coma separated list of modifiers into
  !!* components.
  !!* @param modifier The list of modifers as string.
  !!* @param child The child which carries this modifier (for error messages)
  !!* @param modifiers Array of the modifiers, occuring in modifer.
  !!* @note If the number of the modifiers found differs from the size of the
  !!*   modifiers array, the program stops with error.
  subroutine splitModifier(modifier, child, modifiers)
    character(len=*), intent(in) :: modifier
    type(fnode), pointer :: child
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



!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!! convertByMul
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


  !!* Implementation of convertByMul for real scalar.
  subroutine convertByMul_real(modifier, units, child, value, replace, changed)
    character(len=*), intent(in) :: modifier
    type(unit), intent(in) :: units(:)
    type(fnode), pointer :: child
    real(dp), intent(inout) :: value
    logical, intent(in), optional :: replace
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
      value = value * units(ind)%value
      if (tReplace) then
        call setChildValue(child, "", value, .true.)
      end if
    else
      tChanged = .false.
    end if

    if (present(changed)) then
      changed = tChanged
    end if

  end subroutine convertByMul_real



  !!* Implementation of convertByMul for real rank one array.
  subroutine convertByMul_realR1(modifier, units, child, value, &
      &replace, changed)
    character(len=*), intent(in) :: modifier
    type(unit), intent(in) :: units(:)
    type(fnode), pointer :: child
    real(dp), intent(inout) :: value(:)
    logical, intent(in), optional :: replace
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
      value = value * units(ind)%value
      if (tReplace) then
        call setChildValue(child, "", value, .true.)
      end if
    else
      tChanged = .false.
    end if

    if (present(changed)) then
      changed = tChanged
    end if

  end subroutine convertByMul_realR1

  !!* Implementation of convertByMul for real rank two array.
  subroutine convertByMul_realR2(modifier, units, child, value, &
      &replace, changed)
    character(len=*), intent(in) :: modifier
    type(unit), intent(in) :: units(:)
    type(fnode), pointer :: child
    real(dp), intent(inout) :: value(:,:)
    logical, intent(in), optional :: replace
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
      value = value * units(ind)%value
      if (tReplace) then
        call setChildValue(child, "", value, .true.)
      end if
    else
      tChanged = .false.
    end if

    if (present(changed)) then
      changed = tChanged
    end if

  end subroutine convertByMul_realR2



  !!* Returns a descendant of a given node.
  !!* @param root Node to seek the descendants of
  !!* @param path Path to the descendant. Parents are separated by "/" from
  !!*   their children (e.g. node1/node2/node3)
  !!* @param child Pointer to the child on return or null pointer if not found
  !!* @param requested Should the program stop, if specified descendant is
  !!*   not present (default: .false.)
  !!* @param processed Should elements along the path marked as processed?
  !!*  (default: .false.)
  !!* @param parent If provided, contains parent node of the child, or the last
  !!*  associated node, if the child was not found.
  subroutine getDescendant(root, path, child, requested, processed, parent)
    type(fnode), pointer :: root
    character(len=*), intent(in) :: path
    type(fnode), pointer :: child
    logical, intent(in), optional :: requested
    logical, intent(in), optional :: processed
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



end module hsdutils2

