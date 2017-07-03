!--------------------------------------------------------------------------------------------------!
!  DFTB+: general package for performing fast atomistic simulations                                !
!  Copyright (C) 2017  DFTB+ developers group                                                      !
!                                                                                                  !
!  See the LICENSE file for terms of usage and distribution.                                       !
!--------------------------------------------------------------------------------------------------!

#:include 'common.fypp'

!!* Contains high level functions for converting the values in a XML/HSD
!!* DOM-tree to Fortran intrinsic types.
!!* @todo Some more routines for complex numbers?
module hsdutils
  use assert
  use xmlf90
  use tokenreader
  use hsdparser
  use xmlutils
  use charmanip
  use message
  use linkedlist
  use accuracy
  implicit none
  private

  public :: checkError, detailedError, detailedWarning
  public :: getFirstTextChild, getChildValue, setChildValue
  public :: writeChildValue, getAsString
  public :: convAtomRangeToInt, convRangeToInt, appendPathAndLine
  public :: getChild, getChildren, setChild
  public :: attrProcessed

  !!* Returns the value (the child) of a child node identified by its name.
  !!* @desc
  !!*   These routines investigate the provided node and look for a child with
  !!*   the supplied name. If this child found, its child (which should be a
  !!*   single text node or a usual node if the value argument is of type node)
  !!*   is returned as value converted to the appropriate type. If the child is
  !!*   not found, an error is raised, unless a default value was specified.In
  !!*   that case, a child is created with the provided name and is appended to
  !!*   the node. Furthermore a text node containing the string converted
  !!*   default value is appended to the child node. If default value is
  !!*   provided, it must be also indicated, if the created child is only
  !!*   allowed to have one further child or not. (This corresponds to an
  !!*   assignment with '=' in the HSD input.) If the child (identified by the
  !!*   provided name) is allowed to have a modifier, an argument for the
  !!*   modifier must be provided to contain the the parsed value on return. If
  !!*   the argument for the modifier is missing, but a modifier is found, the
  !!*   program raises an error. The pointer to the found (or created) child can
  !!*   be queried through an appropriate argument. If the name of the child to
  !!*   look for is an empty string, the passed node itself is treated as if it
  !!*   would be the child, which had been found.
  interface getChildValue
    module procedure getChVal_logical
    module procedure getChVal_logicalR1
    module procedure getChVal_node
    module procedure getChVal_string
    module procedure getChVal_lString
    module procedure getChVal_lReal
    module procedure getChVal_lRealR1
    module procedure getChVal_lInt
    module procedure getChVal_lIntR1
    module procedure getChVal_real
    module procedure getChVal_realR1
    module procedure getChVal_realR2
    module procedure getChVal_int
    module procedure getChVal_intR1
    module procedure getChVal_intR2
    module procedure getChVal_lIntR1RealR1
    module procedure getChVal_lStringIntR1RealR1
  end interface


  !!* Sets the value (the child) of a child node identified by its name
  !!* @desc
  !!*   Those functions are the inverse of the getChildValue functions. They
  !!*   create a child with the provided name and append to that child a
  !!*   text node (or a normal node, if the provided value is of type node)
  !!*   containing the provided value. It must be indicated, if the
  !!*   created child is allowed to have only one single further child. If a
  !!*   child with the specified name already exists, the program raises an
  !!*   error, unless replacement flag is set on .true.. In that case, the
  !!*   the existing child is replaced. If the name of the child is the empty
  !!*   string, the current node is treated as if it would be the child, which
  !!*   had been found.
  interface setChildValue
    module procedure setChVal_logical
    module procedure setChVal_logicalR1
    module procedure setChVal_node
    module procedure setChVal_char
    module procedure setChVal_charR1
    module procedure setChVal_real
    module procedure setChVal_realR1
    module procedure setChVal_realR2
    module procedure setChVal_int
    module procedure setChVal_intR1
    module procedure setChVal_intR2
    module procedure setChVal_intR2RealR2
    module procedure setChVal_charR1intR2RealR2
  end interface


  !!* Writes a child and its value to an xml-write stream
  interface writeChildValue
    module procedure writeChVal_logical
    module procedure writeChVal_logicalR1
    module procedure writeChVal_real
    module procedure writeChVal_realR1
    module procedure writeChVal_realR2
    module procedure writeChVal_int
    module procedure writeChVal_intR1
    module procedure writeChVal_intR2
    module procedure writeChVal_intR2RealR2
    module procedure writeChVal_charR1
    module procedure writeChVal_charR1IntR2RealR2
  end interface


  !!* Returns a string representation of an object
  interface getAsString
    module procedure getAsString_logical
    module procedure getAsString_logicalR1
    module procedure getAsString_real
    module procedure getAsString_realR1
    module procedure getAsString_realR2
    module procedure getAsString_int
    module procedure getAsString_intR1
    module procedure getAsString_intR2
    module procedure getAsString_intR2RealR2
    module procedure getAsString_charR1
    module procedure getAsString_charR1IntR2RealR2
  end interface


  !! Error messages
  character(len=*), parameter :: MSG_MISSING_FIELD = "Missing child: "
  character(len=*), parameter :: MSG_EXISTING_CHILD = "Already existing child: "
  character(len=*), parameter :: MSG_NOMODIFIER = "Entity is not allowed to &
      &have a modifier"
  character(len=*), parameter :: MSG_MISSING_VALUES = "Not enough values &
      &provided."

  !!* Length of a line (for wrapping long lines when writing values)
  integer, parameter :: lineLength = 80

  !!* Maximal number of characters needed to represent an integer
  integer, parameter :: nCharInt = 50

  !!* Maximal number of characters needed to represent a real number
  integer, parameter :: nCharReal = 50

  !!* Maximal number of characters needed to represent a logical value
  integer, parameter :: nCharLogical = 4

  !!* Attribute signalising that a tag was processed
  character(len=*), parameter :: attrProcessed = "proc"

  !!* Preallocateated size for temporary buffer strings
  integer, parameter :: preAllocSize = 1024


contains

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!! getChildValue and getNode routines
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !!* Returns the value (the child) of a child node as logical.
  !!* @param node     The node to investigate.
  !!* @param name     Name of the child to look for
  !!* @param value    Value on return
  !!* @param default  Default value for the child, if child is not found
  !!* @param modifier Modifier of the child on return
  !!* @param child    Pointer to the child node (with the spec. name) on return
  subroutine getChVal_logical(node, name, value, default, modifier, child)
    type(fnode), pointer :: node
    character(len=*), intent(in) :: name
    logical, intent(out) :: value
    logical, intent(in), optional :: default
    type(string), intent(inout), optional :: modifier
    type(fnode), pointer, optional :: child

    type(string) :: text, modif
    integer :: iStart, iErr
    type(fnode), pointer :: child2

    @:ASSERT(associated(node))

    child2 => getFirstChildByName(node, tolower(name))
    if (associated(child2)) then
      call getAttribute(child2, attrModifier, modif)
      if (present(modifier)) then
        modifier = modif
      elseif (len(modif) > 0) then
        call detailedError(child2, MSG_NOMODIFIER)
      end if
      iStart = 1
      call getFirstTextChild(child2, text)
      call getNextToken(char(text), value, iStart, iErr)
      call checkError(child2, iErr, "Invalid logical value")
      call checkNoData(child2, char(text), iStart)
      call setAttribute(child2, attrProcessed, "")
    elseif (present(default)) then
      value = default
      if (present(modifier)) then
        modifier = ""
      end if
      call setChildValue(node, name, value, .false., child=child2)
    else
      call detailedError(node, MSG_MISSING_FIELD // name)
    end if
    if (present(child)) then
      child => child2
    end if

  end subroutine getChVal_logical

  !!* Returns the value (the child) of a child node as logical.
  !!* @param node     The node to investigate.
  !!* @param name     Name of the child to look for
  !!* @param value    Value on return
  !!* @param default  Default value for the child, if child is not found
  !!* @param nItem    Nr. of read items. If this argument is not passed, and
  !!*   the nr. of read items is less than the size of the array, the
  !!*   subroutine raises an error.
  !!* @param modifier Modifier of the child on return
  !!* @param child    Pointer to the child node (with the spec. name) on return
  subroutine getChVal_logicalR1(node, name, value, default, nItem, modifier, &
      & child)
    type(fnode), pointer :: node
    character(len=*), intent(in) :: name
    logical, intent(out) :: value(:)
    logical, intent(in), optional :: default(:)
    integer, intent(out), optional :: nItem
    type(string), intent(inout), optional :: modifier
    type(fnode), pointer, optional :: child

    type(string) :: text, modif
    integer :: iStart, iErr, nReadItem
    type(fnode), pointer :: child2

    @:ASSERT(associated(node))
  #:call ASSERT_CODE
    if (present(default)) then
      @:ASSERT(all(shape(default) == shape(value)))
    end if
  #:endcall ASSERT_CODE

    if (present(nItem)) then
      nItem = 0
    end if
    child2 => getFirstChildByName(node, tolower(name))
    if (associated(child2)) then
      call getAttribute(child2, attrModifier, modif)
      if (present(modifier)) then
        modifier = modif
      elseif (len(modif) > 0) then
        call detailedError(child2, MSG_NOMODIFIER)
      end if
      iStart = 1
      call getFirstTextChild(child2, text)
      call getNextToken(char(text), value, iStart, iErr, nReadItem)
      call checkError(child2, iErr, "Invalid logical value")
      call checkNoData(child2, char(text), iStart)
      if (present(nItem)) then
        nItem = nReadItem
      elseif (nReadItem /= size(value)) then
        call detailedError(node, MSG_MISSING_VALUES)
      end if
      call setAttribute(child2, attrProcessed, "")
    elseif (present(default)) then
      value = default
      if (present(modifier)) then
        modifier = ""
      end if
      call setChildValue(node, name, value, .false., child=child2)
    else
      call detailedError(node, MSG_MISSING_FIELD // name)
    end if
    call setAttribute(child2, attrProcessed, "")
    if (present(child)) then
      child => child2
    end if

  end subroutine getChVal_logicalR1

  !!* Returns the value (the child) of a child node as string.
  !!* @param node     The node to investigate.
  !!* @param name     Name of the child to look for
  !!* @param value    Value on return
  !!* @param default  Default value for the child, if child is not found
  !!* @param modifier Modifier of the child on return
  !!* @param child    Pointer to the child node (with the spec. name) on return
  !!* @param child    If true, string contains as many tokens as possible, not
  !!*                 just one (with spaces between the tokens).
  subroutine getChVal_string(node, name, value, default, modifier, child, &
      &multiple)
    type(fnode), pointer :: node
    character(len=*), intent(in) :: name
    type(string), intent(inout) :: value
    character(len=*), intent(in), optional :: default
    type(string), intent(inout), optional :: modifier
    type(fnode), pointer, optional :: child
    logical, intent(in), optional :: multiple

    type(string) :: text, modif
    integer :: iStart, iErr
    type(fnode), pointer :: child2
    logical :: tMultiple

    @:ASSERT(associated(node))

    if (present(multiple)) then
      tMultiple = multiple
    else
      tMultiple = .false.
    end if

    child2 => getFirstChildByName(node, tolower(name))
    if (associated(child2)) then
      call getAttribute(child2, attrModifier, modif)
      if (present(modifier)) then
        modifier = modif
      elseif (len(modif) > 0) then
        call detailedError(child2, MSG_NOMODIFIER)
      end if
      call getFirstTextChild(child2, text)
      if (tMultiple) then
        value = unquote(trim(adjustl(char(text))))
      else
        iStart = 1
        call getNextToken(char(text), value, iStart, iErr)
        call checkError(child2, iErr, "Invalid string value")
        call checkNoData(child2, char(text), iStart)
      end if
      call setAttribute(child2, attrProcessed, "")
    elseif (present(default)) then
      value = default
      if (present(modifier)) then
        modifier = ""
      end if
      call setChildValue(node, name, default, .false., child=child2)
    else
      call detailedError(node, MSG_MISSING_FIELD // name)
    end if
    if (present(child)) then
      child => child2
    end if

  end subroutine getChVal_string



  !!* Returns the value (the child) of a child node as real.
  !!* @param node     The node to investigate.
  !!* @param name     Name of the child to look for
  !!* @param value    Value on return
  !!* @param default  Default value for the child, if child is not found
  !!* @param modifier Modifier of the child on return
  !!* @param child    Pointer to the child node (with the spec. name) on return
  subroutine getChVal_real(node, name, value, default, modifier, child)
    type(fnode), pointer :: node
    character(len=*), intent(in) :: name
    real(dp), intent(out) :: value
    real(dp), intent(in), optional :: default
    type(string), intent(inout), optional :: modifier
    type(fnode), pointer, optional :: child

    type(string) :: text, modif
    integer :: iStart, iErr
    type(fnode), pointer :: child2

    @:ASSERT(associated(node))

    child2 => getFirstChildByName(node, tolower(name))
    if (associated(child2)) then
      call getAttribute(child2, attrModifier, modif)
      if (present(modifier)) then
        modifier = modif
      elseif (len(modif) > 0) then
        call detailedError(child2, MSG_NOMODIFIER)
      end if
      iStart = 1
      call getFirstTextChild(child2, text)
      call getNextToken(char(text), value, iStart, iErr)
      call checkError(child2, iErr, "Invalid real value")
      call checkNoData(child2, char(text), iStart)
      call setAttribute(child2, attrProcessed, "")
    elseif (present(default)) then
      value = default
      if (present(modifier)) then
        modifier = ""
      end if
      call setChildValue(node, name, value, .false., child=child2)
    else
      call detailedError(node, MSG_MISSING_FIELD // name)
    end if
    call setAttribute(child2, attrProcessed, "")
    if (present(child)) then
      child => child2
    end if

  end subroutine getChVal_real



  !!* Returns the value (the child) of a child node as a rank one real array.
  !!* @param node     The node to investigate.
  !!* @param name     Name of the child to look for
  !!* @param value    Value on return
  !!* @param default  Default value for the child, if child is not found
  !!* @param nItem    Nr. of read items. If this argument is not passed, and
  !!*   the nr. of read items is less than the size of the array, the
  !!*   subroutine raises an error.
  !!* @param modifier Modifier of the child on return
  !!* @param child    Pointer to the child node (with the spec. name) on return
  subroutine getChVal_realR1(node, name, value, default, nItem, modifier, &
      &child)
    type(fnode), pointer :: node
    character(len=*), intent(in) :: name
    real(dp), intent(out) :: value(:)
    real(dp), intent(in), optional :: default(:)
    integer, intent(out), optional :: nItem
    type(string), intent(inout), optional :: modifier
    type(fnode), pointer, optional :: child

    type(string) :: text, modif
    integer :: iStart, iErr, nReadItem
    type(fnode), pointer :: child2

    @:ASSERT(associated(node))
  #:call ASSERT_CODE
    if (present(default)) then
      @:ASSERT(all(shape(default) == shape(value)))
    end if
  #:endcall ASSERT_CODE

    if (present(nItem)) then
      nItem = 0
    end if
    child2 => getFirstChildByName(node, tolower(name))
    if (associated(child2)) then
      call getAttribute(child2, attrModifier, modif)
      if (present(modifier)) then
        modifier = modif
      elseif (len(modif) > 0) then
        call detailedError(child2, MSG_NOMODIFIER)
      end if
      iStart = 1
      call getFirstTextChild(child2, text)
      call getNextToken(char(text), value, iStart, iErr, nReadItem)
      call checkError(child2, iErr, "Invalid real value")
      call checkNoData(child2, char(text), iStart)
      if (present(nItem)) then
        nItem = nReadItem
      elseif (nReadItem /= size(value)) then
        call detailedError(node, MSG_MISSING_VALUES)
      end if
      call setAttribute(child2, attrProcessed, "")
    elseif (present(default)) then
      value = default
      if (present(modifier)) then
        modifier = ""
      end if
      call setChildValue(node, name, value, .false., child=child2)
    else
      call detailedError(node, MSG_MISSING_FIELD // name)
    end if
    call setAttribute(child2, attrProcessed, "")
    if (present(child)) then
      child => child2
    end if

  end subroutine getChVal_realR1



  !!* Returns the value (the child) of a child node as a rank two real array.
  !!* @param node     The node to investigate.
  !!* @param name     Name of the child to look for
  !!* @param value    Value on return
  !!* @param default  Default value for the child, if child is not found
  !!* @param nItem    Nr. of read items. If this argument is not passed, and
  !!*   the nr. of read items is less than the size of the array, the
  !!*   subroutine raises an error.
  !!* @param modifier Modifier of the child on return
  !!* @param child    Pointer to the child node (with the spec. name) on return
  !!* @note This is just a wrapper around the rank one version, to make sure
  !!*   that two dimensional arrays are pretty printed. For higher ranked
  !!*   arrays the rank one version should be used with some reshaping after.
  subroutine getChVal_realR2(node, name, value, default, nItem, modifier, &
      &child)
    type(fnode), pointer :: node
    character(len=*), intent(in) :: name
    real(dp), intent(out) :: value(:,:)
    real(dp), intent(in), optional :: default(:,:)
    integer, intent(out), optional :: nItem
    type(string), intent(inout), optional :: modifier
    type(fnode), pointer, optional :: child

    real(dp) :: buffer(size(value))
    integer :: nReadItem
    type(string) :: modif
    type(fnode), pointer :: child2

    @:ASSERT(associated(node))
  #:call ASSERT_CODE
    if (present(default)) then
      @:ASSERT(all(shape(default) == shape(value)))
    end if
  #:endcall ASSERT_CODE

    nReadItem = 0
    value = 0.0_dp
    if (present(default)) then
      call getChildValue(node, name, buffer, reshape(default, shape(buffer)), &
          &nReadItem, modifier=modif, child=child2)
    else
      call getChildValue(node, name, buffer, nItem=nReadItem, modifier=modif, &
          &child=child2)
    end if
    if (present(nItem)) then
      nItem = nReadItem
    elseif (nReadItem /= size(value)) then
      call detailedError(node, MSG_MISSING_VALUES)
    end if
    if (present(modifier)) then
      modifier = modif
    elseif (len(modif) > 0) then
      call detailedError(child2, MSG_NOMODIFIER)
    end if
    value(:,:) = reshape(buffer, shape(value))
    if (present(child)) then
      child => child2
    end if

  end subroutine getChVal_realR2



  !!* Returns the value (the child) of a child node as integer.
  !!* @param node     The node to investigate.
  !!* @param name     Name of the child to look for
  !!* @param value    Value on return
  !!* @param default  Default value for the child, if child is not found
  !!* @param modifier Modifier of the child on return
  !!* @param child    Pointer to the child node (with the spec. name) on return
  subroutine getChVal_int(node, name, value, default, modifier, child)
    type(fnode), pointer :: node
    character(len=*), intent(in) :: name
    integer, intent(out) :: value
    integer, intent(in), optional :: default
    type(string), intent(inout), optional :: modifier
    type(fnode), pointer, optional :: child

    type(string) :: text, modif
    integer :: iStart, iErr
    type(fnode), pointer :: child2

    @:ASSERT(associated(node))

    child2 => getFirstChildByName(node, tolower(name))
    if (associated(child2)) then
      call getAttribute(child2, attrModifier, modif)
      if (present(modifier)) then
        modifier = modif
      elseif (len(modif) > 0) then
        call detailedError(child2, MSG_NOMODIFIER)
      end if
      iStart = 1
      call getFirstTextChild(child2, text)
      call getNextToken(char(text), value, iStart, iErr)
      call checkError(child2, iErr, "Invalid integer value")
      call checkNoData(child2, char(text), iStart)
      call setAttribute(child2, attrProcessed, "")
    elseif (present(default)) then
      value = default
      if (present(modifier)) then
        modifier = ""
      end if
      call setChildValue(node, name, value, .false., child=child2)
    else
      call detailedError(node, MSG_MISSING_FIELD // name)
    end if
    if (present(child)) then
      child => child2
    end if

  end subroutine getChVal_int



  !!* Returns the value (the child) of a child node as a rank one integer array.
  !!* @param node     The node to investigate.
  !!* @param name     Name of the child to look for
  !!* @param value    Value on return
  !!* @param default  Default value for the child, if child is not found
  !!* @param nItem    Nr. of read items. If this argument is not passed, and
  !!*   the nr. of read items is less than the size of the array, the
  !!*   subroutine raises an error.
  !!* @param modifier Modifier of the child on return
  !!* @param child    Pointer to the child node (with the spec. name) on return
  subroutine getChVal_intR1(node, name, value, default, nItem, modifier, &
      &child)
    type(fnode), pointer :: node
    character(len=*), intent(in) :: name
    integer, intent(out) :: value(:)
    integer, intent(in), optional :: default(:)
    integer, intent(out), optional :: nItem
    type(string), intent(inout), optional :: modifier
    type(fnode), pointer, optional :: child

    type(string) :: text, modif
    integer :: iStart, iErr, nReadItem
    type(fnode), pointer :: child2

    @:ASSERT(associated(node))
  #:call ASSERT_CODE
    if (present(default)) then
      @:ASSERT(all(shape(default) == shape(value)))
    end if
  #:endcall ASSERT_CODE

    if (present(nItem)) then
      nItem = 0
    end if
    child2 => getFirstChildByName(node, tolower(name))
    if (associated(child2)) then
      call getAttribute(child2, attrModifier, modif)
      if (present(modifier)) then
        modifier = modif
      elseif (len(modif) > 0) then
        call detailedError(child2, MSG_NOMODIFIER)
      end if
      iStart = 1
      call getFirstTextChild(child2, text)
      call getNextToken(char(text), value, iStart, iErr, nReadItem)
      call checkError(child2, iErr, "Invalid integer value")
      call checkNoData(child2, char(text), iStart)
      if (present(nItem)) then
        nItem = nReadItem
      elseif (nReadItem /= size(value)) then
        call detailedError(node, MSG_MISSING_VALUES)
      end if
      call setAttribute(child2, attrProcessed, "")
    elseif (present(default)) then
      value = default
      if (present(modifier)) then
        modifier = ""
      end if
      call setChildValue(node, name, value, .false., child=child2)
    else
      call detailedError(node, MSG_MISSING_FIELD // name)
    end if
    if (present(child)) then
      child => child2
    end if

  end subroutine getChVal_intR1



  !!* Returns the value (the child) of a child node as a rank two integer array.
  !!* @param node     The node to investigate.
  !!* @param name     Name of the child to look for
  !!* @param value    Value on return
  !!* @param default  Default value for the child, if child is not found
  !!* @param nItem    Nr. of read items. If this argument is not passed, and
  !!*   the nr. of read items is less than the size of the array, the
  !!*   subroutine raises an error.
  !!* @param modifier Modifier of the child on return
  !!* @param child    Pointer to the child node (with the spec. name) on return
  !!* @note This is just a wrapper around the rank one version, to make sure
  !!*   that two dimensional arrays are pretty printed. For higher ranked
  !!*   arrays the rank one version should be used with some reshaping after.
  subroutine getChVal_intR2(node, name, value, default, nItem, modifier, &
      &child)
    type(fnode), pointer :: node
    character(len=*), intent(in) :: name
    integer, intent(out) :: value(:,:)
    integer, intent(in), optional :: default(:,:)
    integer, intent(out), optional :: nItem
    type(string), intent(inout), optional :: modifier
    type(fnode), pointer, optional :: child

    integer :: buffer(size(value))
    integer :: nReadItem
    type(string) :: modif
    type(fnode), pointer :: child2

    @:ASSERT(associated(node))
  #:call ASSERT_CODE
    if (present(default)) then
      @:ASSERT(all(shape(default) == shape(value)))
    end if
  #:endcall ASSERT_CODE

    nReadItem = 0
    if (present(default)) then
      call getChildValue(node, name, buffer, reshape(default, shape(buffer)), &
          &nReadItem, modif, child=child2)
    else
      call getChildValue(node, name, buffer, nItem=nReadItem, modifier=modif, &
          &child=child2)
    end if
    if (present(nItem)) then
      nItem = nReadItem
    elseif (nReadItem /= size(value)) then
      call detailedError(node, MSG_MISSING_VALUES)
    end if
    if (present(modifier)) then
      modifier = modif
    elseif (len(modif) > 0) then
      call detailedError(child2, MSG_NOMODIFIER)
    end if
    value(:,:) = reshape(buffer, shape(value))
    if (present(child)) then
      child => child2
    end if

  end subroutine getChVal_intR2



  !!* Returns the value (the child) of a child node as a linked list of strings.
  !!* @param node     The node to investigate.
  !!* @param name     Name of the child to look for
  !!* @param value    Value on return
  !!* @param modifier Modifier of the child on return
  !!* @param child    Pointer to the child node (with the spec. name) on return
  !!* @note In order to prevent a double packaging (from array to linked list
  !!*   and then from linked list to array), the setting of defaults for list
  !!*   types is not allowed. The presence of the child must be explicitely
  !!*   queried in the caller routine and an eventual default setting must be
  !!*   set with an explicit setChildValue call.
  subroutine getChVal_lString(node, name, value, modifier, child)
    type(fnode), pointer :: node
    character(len=*), intent(in) :: name
    type(listString), intent(inout) :: value
    type(string), intent(inout), optional :: modifier
    type(fnode), pointer, optional :: child

    type(string) :: text, modif
    type(fnode), pointer :: child2

    @:ASSERT(associated(node))
    child2 => getFirstChildByName(node, tolower(name))
    if (associated(child2)) then
      call getAttribute(child2, attrModifier, modif)
      if (present(modifier)) then
        modifier = modif
      elseif (len(modif) > 0) then
        call detailedError(child2, MSG_NOMODIFIER)
      end if
      call getFirstTextChild(child2, text)
      call getChVal_lString_h(char(text), value, child2)
      call setAttribute(child2, attrProcessed, "")
    else
      call detailedError(node, MSG_MISSING_FIELD // name)
    end if
    if (present(child)) then
      child => child2
    end if

  end subroutine getChVal_lString



  !!* Helper function for getChVal_lString to avoid string to character
  !!* conversion in the do-loop.
  !!* @param text  Text to parse
  !!* @param value Contains the value of the parsed text
  subroutine getChVal_lString_h(text, value, node)
    character(len=*), intent(in) :: text
    type(listString), intent(inout) :: value
    type(fnode), pointer :: node

    integer :: iStart, iErr
    type(string) :: token

    iStart = 1
    call getNextToken(text, token, iStart, iErr)
    do while (iErr == TOKEN_OK)
      call append(value, trim(unquote(char(token))))
      call getNextToken(text, token, iStart, iErr)
    end do
    if (iErr == TOKEN_ERROR) then
      call detailedError(node, "Invalid string")
    end if

  end subroutine getChVal_lString_h



  !!* Returns the value (the child) of a child node as a linked list of reals.
  !!* @param node     The node to investigate.
  !!* @param name     Name of the child to look for
  !!* @param value    Value on return
  !!* @param modifier Modifier of the child on return
  !!* @param child    Pointer to the child node (with the spec. name) on return
  !!* @note In order to prevent a double packaging (from array to linked list
  !!*   and then from linked list to array), the setting of defaults for list
  !!*   types is not allowed. The presence of the child must be explicitely
  !!*   queried in the caller routine and an eventual default setting must be
  !!*   set with an explicit setChildValue call.
  subroutine getChVal_lReal(node, name, value, modifier, child)
    type(fnode), pointer :: node
    character(len=*), intent(in) :: name
    type(listReal), intent(inout) :: value
    type(string), intent(inout), optional :: modifier
    type(fnode), pointer, optional :: child

    type(string) :: text, modif
    type(fnode), pointer :: child2

    @:ASSERT(associated(node))

    child2 => getFirstChildByName(node, tolower(name))
    if (associated(child2)) then
      call getAttribute(child2, attrModifier, modif)
      if (present(modifier)) then
        modifier = modif
      elseif (len(modif) > 0) then
        call detailedError(child2, MSG_NOMODIFIER)
      end if
      call getFirstTextChild(child2, text)
      call getChVal_lReal_h(char(text), value, child2)
      call setAttribute(child2, attrProcessed, "")
    else
      call detailedError(node, MSG_MISSING_FIELD // name)
    end if
    if (present(child)) then
      child => child2
    end if

  end subroutine getChVal_lReal



  !!* Helper function for getChVal_lReal to avoid string to character
  !!* conversion in the do-loop.
  !!* @param text  Text to parse
  !!* @param value Contains the value of the parsed text
  subroutine getChVal_lReal_h(text, value, node)
    character(len=*), intent(in) :: text
    type(listReal), intent(inout) :: value
    type(fnode), pointer :: node

    integer :: iStart, iErr
    real(dp) :: buffer

    iStart = 1
    call getNextToken(text, buffer, iStart, iErr)
    do while (iErr == TOKEN_OK)
      call append(value, buffer)
      call getNextToken(text, buffer, iStart, iErr)
    end do
    if (iErr == TOKEN_ERROR) then
      call detailedError(node, "Invalid real value")
    end if

  end subroutine getChVal_lReal_h



  !!* Returns the value (the child) of a child node as a linked list of rank one
  !!* real arrays.
  !!* @param node     The node to investigate.
  !!* @param name     Name of the child to look for
  !!* @param dim      Dimension of the arrays
  !!* @param value    Value on return
  !!* @param modifier Modifier of the child on return
  !!* @param child    Pointer to the child node (with the spec. name) on return
  !!* @note In order to prevent a double packaging (from array to linked list
  !!*   and then from linked list to array), the setting of defaults for list
  !!*   types is not allowed. The presence of the child must be explicitely
  !!*   queried in the caller routine and an eventual default setting must be
  !!*   set with an explicit setChildValue call.
  subroutine getChVal_lRealR1(node, name, dim, value, modifier, child)
    type(fnode), pointer :: node
    character(len=*), intent(in) :: name
    integer, intent(in) :: dim
    type(listRealR1), intent(inout) :: value
    type(string), intent(inout), optional :: modifier
    type(fnode), pointer, optional :: child

    type(string) :: text, modif
    type(fnode), pointer :: child2

    @:ASSERT(associated(node))

    child2 => getFirstChildByName(node, tolower(name))
    if (associated(child2)) then
      call getAttribute(child2, attrModifier, modif)
      if (present(modifier)) then
        modifier = modif
      elseif (len(modif) > 0) then
        call detailedError(child2, MSG_NOMODIFIER)
      end if
      call getFirstTextChild(child2, text)
      call getChVal_lRealR1_h(char(text), dim, value, child2)
      call setAttribute(child2, attrProcessed, "")
    else
      call detailedError(node, MSG_MISSING_FIELD // name)
    end if
    if (present(child)) then
      child => child2
    end if

  end subroutine getChVal_lRealR1



  !!* Helper function for getChVal_lReal to avoid string to character
  !!* conversion in the do-loop.
  !!* @param text  Text to parse
  !!* @param value Contains the value of the parsed text
  subroutine getChVal_lRealR1_h(text, dim, value, node)
    character(len=*), intent(in) :: text
    integer, intent(in) :: dim
    type(listRealR1), intent(inout) :: value
    type(fnode), pointer :: node

    integer :: iStart, iErr
    real(dp) :: buffer(dim)
    integer :: nItem

    iStart = 1
    call getNextToken(text, buffer, iStart, iErr, nItem)
    do while (iErr == TOKEN_OK)
      call append(value, buffer)
      call getNextToken(text, buffer, iStart, iErr, nItem)
    end do
    if (iErr == TOKEN_ERROR) then
      call detailedError(node, "Invalid real value")
    elseif (iErr == TOKEN_EOS .and. nItem /= 0) then
      call detailedError(node, "Unexpected end of data")
    end if

  end subroutine getChVal_lRealR1_h



  !!* Returns the value (the child) of a child node as linked list of integers.
  !!* @param node     The node to investigate.
  !!* @param name     Name of the child to look for
  !!* @param value    Value on return
  !!* @param modifier Modifier of the child on return
  !!* @param child    Pointer to the child node (with the spec. name) on return
  !!* @note In order to prevent a double packaging (from array to linked list
  !!*   and then from linked list to array), the setting of defaults for list
  !!*   types is not allowed. The presence of the child must be explicitely
  !!*   queried in the caller routine and an eventual default setting must be
  !!*   set with an explicit setChildValue call.
  subroutine getChVal_lInt(node, name, value, modifier, child)
    type(fnode), pointer :: node
    character(len=*), intent(in) :: name
    type(listInt), intent(inout) :: value
    type(string), intent(inout), optional :: modifier
    type(fnode), pointer, optional :: child

    type(string) :: text, modif
    type(fnode), pointer :: child2

    @:ASSERT(associated(node))

    child2 => getFirstChildByName(node, tolower(name))
    if (associated(child2)) then
      call getAttribute(child2, attrModifier, modif)
      if (present(modifier)) then
        modifier = modif
      elseif (len(modif) > 0) then
        call detailedError(child2, MSG_NOMODIFIER)
      end if
      call getFirstTextChild(child2, text)
      call getChVal_lInt_h(char(text), value, child2)
      call setAttribute(child2, attrProcessed, "")
    else
      call detailedError(node, MSG_MISSING_FIELD // name)
    end if
    if (present(child)) then
      child => child2
    end if

  end subroutine getChVal_lInt



  !!* Helper function for getChVal_lReal to avoid string to character
  !!* conversion in the do-loop.
  !!* @param text  Text to parse
  !!* @param value Contains the value of the parsed text
  subroutine getChVal_lInt_h(text, value, node)
    character(len=*), intent(in) :: text
    type(listInt), intent(inout) :: value
    type(fnode), pointer :: node

    integer :: iStart, iErr
    integer :: buffer

    iStart = 1
    call getNextToken(text, buffer, iStart, iErr)
    do while (iErr == TOKEN_OK)
      call append(value, buffer)
      call getNextToken(text, buffer, iStart, iErr)
    end do
    if (iErr == TOKEN_ERROR) then
      call detailedError(node, "Invalid real value")
    end if

  end subroutine getChVal_lInt_h



  !!* Returns the value (the child) of a child node as linked list of rank one
  !!* integer arrays.
  !!* @param node     The node to investigate.
  !!* @param name     Name of the child to look for
  !!* @param value    Value on return
  !!* @param modifier Modifier of the child on return
  !!* @param child    Pointer to the child node (with the spec. name) on return
  !!* @note In order to prevent a double packaging (from array to linked list
  !!*   and then from linked list to array), the setting of defaults for list
  !!*   types is not allowed. The presence of the child must be explicitely
  !!*   queried in the caller routine and an eventual default setting must be
  !!*   set with an explicit setChildValue call.
  subroutine getChVal_lIntR1(node, name, dim, value, modifier, child)
    type(fnode), pointer :: node
    character(len=*), intent(in) :: name
    integer, intent(in) :: dim
    type(listIntR1), intent(inout) :: value
    type(string), intent(inout), optional :: modifier
    type(fnode), pointer, optional :: child

    type(string) :: text, modif
    type(fnode), pointer :: child2

    @:ASSERT(associated(node))

    child2 => getFirstChildByName(node, tolower(name))
    if (associated(child2)) then
      call getAttribute(child2, attrModifier, modif)
      if (present(modifier)) then
        modifier = modif
      elseif (len(modif) > 0) then
        call detailedError(child2, MSG_NOMODIFIER)
      end if
      call getFirstTextChild(child2, text)
      call getChVal_lIntR1_h(char(text), dim, value, child2)
      call setAttribute(child2, attrProcessed, "")
    else
      call detailedError(node, MSG_MISSING_FIELD // name)
    end if
    if (present(child)) then
      child => child2
    end if

  end subroutine getChVal_lIntR1



  !!* Helper function for getChVal_lReal to avoid string to character
  !!* conversion in the do-loop.
  !!* @param text  Text to parse
  !!* @param value Contains the value of the parsed text
  subroutine getChVal_lIntR1_h(text, dim, value, node)
    character(len=*), intent(in) :: text
    integer, intent(in) :: dim
    type(listIntR1), intent(inout) :: value
    type(fnode), pointer :: node

    integer :: iStart, iErr
    integer :: buffer(dim)
    integer :: nItem

    iStart = 1
    call getNextToken(text, buffer, iStart, iErr, nItem)
    do while (iErr == TOKEN_OK)
      call append(value, buffer)
      call getNextToken(text, buffer, iStart, iErr, nItem)
    end do
    if (iErr == TOKEN_ERROR) then
      call detailedError(node, "Invalid real value")
    elseif (iErr == TOKEN_EOS .and. nItem /= 0) then
      call detailedError(node, "Unexpected end of data")
    end if

  end subroutine getChVal_lIntR1_h



  !!* Returns the value (the child) of a child node as a linked list rank one
  !!* integer and rank one real arrays.
  !!* @param node      The node to investigate.
  !!* @param name      Name of the child to look for
  !!* @param dimInt    Dimension of the integer arrays in the list
  !!* @param valueInt  List of integer arrays on return
  !!* @param dimReal   Dimensio of the real arrays in the list
  !!* @param valueReal List of real array on return
  !!* @param modifier  Modifier of the child on return
  !!* @param child     Pointer to the child on return
  !!* @note In order to prevent a double packaging (from array to linked list
  !!*   and then from linked list to array), the setting of defaults for list
  !!*   types is not allowed. The presence of the child must be explicitely
  !!*   queried in the caller routine and an eventual default setting must be
  !!*   set with an explicit setChildValue call.
  subroutine getChVal_lIntR1RealR1(node, name, dimInt, valueInt, &
      &dimReal, valueReal, modifier, child)
    type(fnode), pointer :: node
    character(len=*), intent(in) :: name
    integer, intent(in) :: dimInt
    integer, intent(in) :: dimReal
    type(listIntR1), intent(inout) :: valueInt
    type(listRealR1), intent(inout) :: valueReal
    type(string), intent(inout), optional :: modifier
    type(fnode), pointer, optional :: child

    type(string) :: text, modif
    type(fnode), pointer :: child2

    @:ASSERT(associated(node))
    @:ASSERT(dimInt > 0)
    @:ASSERT(dimReal > 0)

    child2 => getFirstChildByName(node, tolower(name))
    if (associated(child2)) then
      call getAttribute(child2, attrModifier, modif)
      if (present(modifier)) then
        modifier = modif
      elseif (len(modif) > 0) then
        call detailedError(child2, MSG_NOMODIFIER)
      end if
      call getFirstTextChild(child2, text)
      call getChVal_lIntR1RealR1_h(char(text), dimInt, valueInt, &
          &dimReal, valueReal, child2)
      if (len(valueInt) /= len(valueReal)) then
        call detailedError(node, "Unexpected end of data")
      end if
      call setAttribute(child2, attrProcessed, "")
    else
      call detailedError(node, MSG_MISSING_FIELD // name)
    end if
    if (present(child)) then
      child => child2
    end if

  end subroutine getChVal_lIntR1RealR1



  !!* Helper function for getChVal_lIntR1RealR1 to avoid string to char
  !!* conversion in the do-loop.
  !!* @param text  Text to parse
  !!* @param value Contains the value of the parsed text
  subroutine getChVal_lIntR1RealR1_h(text, dimInt, valueInt, &
      &dimReal, valueReal, node)
    character(len=*), intent(in) :: text
    integer, intent(in) :: dimInt
    type(listIntR1), intent(inout) :: valueInt
    integer, intent(in) :: dimReal
    type(listRealR1), intent(inout) :: valueReal
    type(fnode), pointer :: node

    integer :: iStart, iErr
    real(dp) :: bufferReal(dimReal)
    integer :: bufferInt(dimInt)
    integer :: nItem

    iErr = TOKEN_OK
    iStart = 1
    do while (iErr == TOKEN_OK)
      call getNextToken(text, bufferInt, iStart, iErr, nItem)
      if (iErr == TOKEN_ERROR) then
        call detailedError(node, "Invalid integer")
      elseif (iErr == TOKEN_EOS .and. nItem /= 0) then
        call detailedError(node, "Unexpected end of data")
      end if
      if (iErr == TOKEN_OK) then
        call append(valueInt, bufferInt)
        call getNextToken(text, bufferReal, iStart, iErr, nItem)
        call checkError(node, iErr, "Invalid real")
        if (iErr == TOKEN_OK) then
          call append(valueReal, bufferReal)
        end if
      end if
    end do

  end subroutine getChVal_lIntR1RealR1_h



  !!* Returns the value (the child) of a child node as a linked list of string,
  !!* rank one integer and rank one real arrays.
  !!* @param node      The node to investigate.
  !!* @param name      Name of the child to look for
  !!* @param valueStr  List of strings on return.
  !!* @param dimInt    Dimension of the integer arrays in the list
  !!* @param valueInt  List of integer arrays on return
  !!* @param dimReal   Dimension of the real arrays in the list
  !!* @param valueReal List of real array on return
  !!* @param modifier  Modifier of the child on return
  !!* @param child     Pointer to the child on return
  !!* @note In order to prevent a double packaging (from array to linked list
  !!*   and then from linked list to array), the setting of defaults for list
  !!*   types is not allowed. The presence of the child must be explicitely
  !!*   queried in the caller routine and an eventual default setting must be
  !!*   set with an explicit setChildValue call.
  subroutine getChVal_lStringIntR1RealR1(node, name, valueStr, dimInt, &
      &valueInt, dimReal, valueReal, modifier, child)
    type(fnode), pointer :: node
    character(len=*), intent(in) :: name
    type(listString), intent(inout) :: valueStr
    integer, intent(in) :: dimInt
    type(listIntR1), intent(inout) :: valueInt
    integer, intent(in) :: dimReal
    type(listRealR1), intent(inout) :: valueReal
    type(string), intent(inout), optional :: modifier
    type(fnode), pointer, optional :: child

    type(string) :: text, modif
    type(fnode), pointer :: child2

    @:ASSERT(associated(node))
    @:ASSERT(dimInt > 0)
    @:ASSERT(dimReal > 0)

    child2 => getFirstChildByName(node, tolower(name))
    if (associated(child2)) then
      call getAttribute(child2, attrModifier, modif)
      if (present(modifier)) then
        modifier = modif
      elseif (len(modif) > 0) then
        call detailedError(child2, MSG_NOMODIFIER)
      end if
      call getFirstTextChild(child2, text)
      call getChVal_lStringIntR1RealR1_h(char(text), valueStr, &
          &dimInt, valueInt, dimReal, valueReal, child2)
      if (len(valueStr) /= len(valueInt) &
          &.or. len(valueInt) /= len(valueReal)) then
        call detailedError(node, "Unexpected end of data")
      end if
      call setAttribute(child2, attrProcessed, "")
    else
      call detailedError(node, MSG_MISSING_FIELD // name)
    end if
    if (present(child)) then
      child => child2
    end if

  end subroutine getChVal_lStringIntR1RealR1



  !!* Helper function for getChVal_lIntR1RealR1 to avoid string to char
  !!* conversion in the do-loop.
  !!* @param text  Text to parse
  !!* @param value Contains the value of the parsed text
  subroutine getChVal_lStringIntR1RealR1_h(text, valueStr, dimInt, &
      &valueInt, dimReal, valueReal, node)
    character(len=*), intent(in) :: text
    type(listString), intent(inout) :: valueStr
    integer, intent(in) :: dimInt
    type(listIntR1), intent(inout) :: valueInt
    integer, intent(in) :: dimReal
    type(listRealR1), intent(inout) :: valueReal
    type(fnode), pointer :: node

    integer :: iStart, iErr
    real(dp) :: bufferReal(dimReal)
    integer :: bufferInt(dimInt)
    integer :: nItem
    type(string) :: bufferStr

    iErr = TOKEN_OK
    iStart = 1
    do while (iErr == TOKEN_OK)
      call getNextToken(text, bufferStr, iStart, iErr)
      if (iErr == TOKEN_ERROR) then
        call detailedError(node, "Invalid string")
      elseif (iErr == TOKEN_EOS) then
        exit
      end if
      call append(valueStr, char(bufferStr))

      call getNextToken(text, bufferInt, iStart, iErr, nItem)
      call checkError(node, iErr, "Invalid integer")
      call append(valueInt, bufferInt)

      call getNextToken(text, bufferReal, iStart, iErr, nItem)
      call checkError(node, iErr, "Invalid real")
      call append(valueReal, bufferReal)
    end do

  end subroutine getChVal_lStringIntR1RealR1_h



  !!* Returns the value (the child) of a child node as a node.
  !!* @param node     The node to investigate.
  !!* @param name     Name of the child to look for
  !!* @param value    Value on return
  !!* @param default  Default value for the child, if child is not found. If
  !!*   the empty string is passed as default value, the child is created
  !!*   but no value is added to it. The returned value pointer will be
  !!*   unassociated. (allowEmptyValue must be explicitely set to .true.)
  !!* @param modifier Modifier of the child on return
  !!* @param child    Pointer to the child node (with the spec. name) on return
  !!* @param list     If the node created as default should be tagged as list.
  !!* @param allowEmptyValue If the child is allowed to have an empty value.
  !!* @param dummyValue If true, the value is not marked as processed.
  !!* @caveat If allowEmptyValue is set to .true. and the child has no subnodes
  !!*   (empty value) the returned value is an unassociated pointer!
  subroutine getChVal_node(node, name, value, default, modifier, child, &
      &list, allowEmptyValue, dummyValue)
    type(fnode), pointer :: node
    character(len=*), intent(in) :: name
    type(fnode), pointer :: value
    character(len=*), intent(in), optional :: default
    type(string), intent(inout), optional :: modifier
    type(fnode), pointer, optional :: child
    logical, intent(in), optional :: list
    logical, intent(in), optional :: allowEmptyValue
    logical, intent(in), optional :: dummyValue

    type(string) :: modif
    type(fnode), pointer :: child2
    logical :: tList, tAllowEmptyVal, tDummyValue

    @:ASSERT(associated(node))
  #:call ASSERT_CODE
    if (present(default)) then
      if (len(default) == 0) then
        @:ASSERT(present(allowEmptyValue))
        @:ASSERT(allowEmptyValue)
      end if
    end if
  #:endcall ASSERT_CODE

    if (present(list)) then
      tList = list
    else
      tList = .false.
    end if
    if (present(allowEmptyValue)) then
      tAllowEmptyVal = allowEmptyValue
    else
      tAllowEmptyVal = .false.
    end if
    if (present(dummyValue)) then
      tDummyValue = dummyValue
    else
      tDummyValue = .false.
    end if

    child2 => getFirstChildByName(node, tolower(name))
    if (associated(child2)) then
      call getAttribute(child2, attrModifier, modif)
      if (present(modifier)) then
        modifier = modif
      elseif (len(modif) > 0) then
        call detailedError(child2, MSG_NOMODIFIER)
      end if
      value => getFirstChild(child2)
      if ((.not. associated(value)) .and. (.not. tAllowEmptyVal)) then
        call detailedError(child2, "Missing value")
      end if
      call setAttribute(child2, attrProcessed, "")
    elseif (present(default)) then
      if (present(modifier)) then
        modifier = ""
      end if
      if (len(default) > 0) then
        value => createElement(tolower(default))
        call setChildValue(node, name, value, .false., child=child2, list=tList)
        call setAttribute(value, attrName, default)
      else
        nullify(value)
        call setChild(node, name, child2, .false., list=tList)
      end if
    else
      call detailedError(node, MSG_MISSING_FIELD // name)
    end if
    if (associated(value) .and. .not. tDummyValue) then
      if (getNodeType(value) == ELEMENT_NODE) then
        call setAttribute(value, attrProcessed, "")
      end if
    end if
    if (present(child)) then
      child => child2
    end if

  end subroutine getChVal_node


  !!* Converts a string containing atom indices, ranges and species names to
  !!* a list of atom indeces.
  !!* @param str  String to convert
  !!* @param speciesNames  Contains the valid species names.
  !!* @param species  Contains for every atom its species index
  !!* @param node  Master node for detailed errors.
  !!* @param val  Integer list of atom indices on return.
  subroutine convAtomRangeToInt(str, speciesNames, species, node, val)
    character(len=*), intent(in) :: str
    character(len=*), intent(in) :: speciesNames(:)
    integer, intent(in) :: species(:)
    type(fnode), pointer :: node
    integer, allocatable, intent(out) :: val(:)

    type(string) :: buffer
    type(ListInt) :: li
    integer :: nAtom, iStart, iostat

    nAtom = size(species)
    call init(li)
    iStart = 1
    call getNextToken(str, buffer, iStart, iostat)
    do while (iostat == TOKEN_OK)
      call process(char(buffer), speciesNames, species, nAtom, node, li)
      call getNextToken(str, buffer, iStart, iostat)
    end do
    allocate(val(len(li)))
    if (len(li) > 0) then
      call asArray(li, val)
    end if
    call destruct(li)

  contains

    ! Helper routine.
    subroutine process(cbuffer, speciesNames, species, nAtom, node, li)
      character(len=*), intent(in) :: cbuffer
      character(len=*), intent(in) :: speciesNames(:)
      integer, intent(in) :: species(:)
      integer, intent(in) :: nAtom
      type(fnode), pointer :: node
      type(ListInt), intent(inout) :: li

      integer :: iPos, bounds(2), iSp, ii
      integer :: iStart1, iStart2, iost(2)

      if ((cbuffer(1:1) >= "0" .and. cbuffer(1:1) <= "9") &
          &.or. cbuffer(1:1) == "-") then
        iPos = scan(cbuffer, ":")
        if (iPos /= 0) then
          iStart1 = 1
          iStart2 = iPos + 1
          call getNextToken(cbuffer(1:iPos-1), bounds(1), iStart1, iost(1))
          call getNextToken(cbuffer, bounds(2), iStart2, iost(2))
          if (any(iost /= TOKEN_OK)) then
            call detailedError(node, "Invalid range specification '" &
                &// trim(cbuffer) // "'")
          end if
          if (any(bounds > nAtom) .or. any(bounds < -nAtom) &
              &.or. any(bounds == 0)) then
            call detailedError(node, "Specified number out of range in '" &
                &// trim(cbuffer) // "'")
          end if
          bounds = modulo(bounds, nAtom + 1)
          if (bounds(1) > bounds(2)) then
            call detailedError(node, "Negative range '" // trim(cbuffer) &
                &// "'")
          end if
          do ii = bounds(1), bounds(2)
            call append(li, ii)
          end do
        else
          iStart1 = 1
          call getNextToken(cbuffer, ii, iStart1, iost(1))
          if (iost(1) /= TOKEN_OK) then
            call detailedError(node, "Invalid integer '" // trim(cbuffer) &
                &// "'")
          end if
          if (ii > nAtom .or. ii < -nAtom .or. ii == 0) then
            call detailedError(node, "Specified number (" // trim(cbuffer) // &
                &") out of range.")
          end if
          ii = modulo(ii, nAtom + 1)
          call append(li, ii)
        end if
      else
        ! Try to interprete it as a species name
        iPos = 0
        do iSp = 1, size(speciesNames)
          if (speciesNames(iSp) == cbuffer) then
            iPos = iSp
            exit
          end if
        end do
        if (iPos == 0) then
          call detailedError(node, "Invalid species name '" // trim(cbuffer) &
              &// "'")
        end if
        do ii = 1, size(species)
          if (species(ii) == iPos) then
            call append(li, ii)
          end if
        end do
      end if

    end subroutine process

  end subroutine convAtomRangeToInt


  !!* Converts a string containing indices and ranges to a list of indices.
  !!* @param str  String to convert
  !!* @param node  Master node for detailed errors.
  !!* @param val  Integer list of atom indices on return.
  !!* @param nMax Maximum number for an index
  subroutine convRangeToInt(str, node, val, nMax)
    character(len=*), intent(in) :: str
    type(fnode), pointer :: node
    integer, allocatable, intent(out) :: val(:)
    integer, intent(in) :: nMax

    type(string) :: buffer
    type(ListInt) :: li
    integer :: iStart, iostat

    call init(li)
    iStart = 1
    call getNextToken(str, buffer, iStart, iostat)
    do while (iostat == TOKEN_OK)
      call process(char(buffer), nMax, node, li)
      call getNextToken(str, buffer, iStart, iostat)
    end do
    allocate(val(len(li)))
    if (len(li) > 0) then
      call asArray(li, val)
    end if
    call destruct(li)

  contains

    ! Helper routine.
    subroutine process(cbuffer, nMax, node, li)
      character(len=*), intent(in) :: cbuffer
      integer, intent(in) :: nMax
      type(fnode), pointer :: node
      type(ListInt), intent(inout) :: li

      integer :: iPos, bounds(2), ii
      integer :: iStart1, iStart2, iost(2)

      if ((cbuffer(1:1) >= "0" .and. cbuffer(1:1) <= "9") &
          &.or. cbuffer(1:1) == "-") then
        iPos = scan(cbuffer, ":")
        if (iPos /= 0) then
          iStart1 = 1
          iStart2 = iPos + 1
          call getNextToken(cbuffer(1:iPos-1), bounds(1), iStart1, iost(1))
          call getNextToken(cbuffer, bounds(2), iStart2, iost(2))
          if (any(iost /= TOKEN_OK)) then
            call detailedError(node, "Invalid range specification '" &
                &// trim(cbuffer) // "'")
          end if
          if (any(bounds > nMax) .or. any(bounds < -nMax) &
              &.or. any(bounds == 0)) then
            call detailedError(node, "Specified number out of range in '" &
                &// trim(cbuffer) // "'")
          end if
          bounds = modulo(bounds, nMax + 1)
          if (bounds(1) > bounds(2)) then
            call detailedError(node, "Negative range '" // trim(cbuffer) &
                &// "'")
          end if
          do ii = bounds(1), bounds(2)
            call append(li, ii)
          end do
        else
          iStart1 = 1
          call getNextToken(cbuffer, ii, iStart1, iost(1))
          if (iost(1) /= TOKEN_OK) then
            call detailedError(node, "Invalid integer '" // trim(cbuffer) &
                &// "'")
          end if
          if (ii > nMax .or. ii < -nMax .or. ii == 0) then
            call detailedError(node, "Specified number (" // trim(cbuffer) // &
                &") out of range.")
          end if
          call append(li, ii)
        end if
      else
        call detailedError(node, "Invalid range '" // trim(cbuffer) // "'")
      end if

    end subroutine process

  end subroutine convRangeToInt




  !!* Returns a child node with a specified name
  !!* @param node      Node to investigate
  !!* @param name      Name of the child node to look for
  !!* @param child     Contains a pointer to the child on return
  !!* @param requested If true and child not found, error is issued
  !!* @param modifier  Contains modifier on exit.
  subroutine getChild(node, name, child, requested, modifier)
    type(fnode), pointer :: node
    character(len=*), intent(in) :: name
    type(fnode), pointer :: child
    logical, intent(in), optional :: requested
    type(string), intent(inout), optional :: modifier

    logical :: tRequested
    type(string) :: modif

    @:ASSERT(associated(node))

    if (present(requested)) then
      tRequested = requested
    else
      tRequested = .true.
    end if

    child => getFirstChildByName(node, tolower(name))
    if (associated(child)) then
      call getAttribute(child, attrModifier, modif)
      if (present(modifier)) then
        modifier = modif
      elseif (len(modif) > 0) then
        call detailedError(child, MSG_NOMODIFIER)
      end if
      call setAttribute(child, attrProcessed, "")
    elseif (tRequested) then
      call detailedError(node, MSG_MISSING_FIELD // name)
    end if

  end subroutine getChild



  !!* Returns a list of children with the specified name.
  !!* @param node     Parent node to investigate
  !!* @param name     Name of the children to look for
  !!* @param children List of the children.
  subroutine getChildren(node, name, children)
    type(fnode), pointer :: node
    character(len=*), intent(in) :: name
    type(fnodeList), pointer :: children

    type(fnode), pointer :: child
    integer :: ii

    children => getChildrenByName(node, tolower(name))
    do ii = 1, getLength(children)
      call getItem1(children, ii, child)
      call setAttribute(child, attrProcessed, "")
    end do

  end subroutine getChildren



!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!! setChildValue and setNode routines
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !!* Sets the value (child) of a child with given name.
  !!* @param node    The node to investigate
  !!* @param name    Name of the child to look for
  !!* @param value   Value to set
  !!* @param replace Replace if child with same name already exists
  !!* @param child   Pointer to the child node (with the provided name)
  !!* @param modifier Optional modifier for the child
  subroutine setChVal_logical(node, name, value, replace, child, modifier)
    type(fnode), pointer :: node
    character(len=*), intent(in) :: name
    logical, intent(in) :: value
    logical, intent(in), optional :: replace
    type(fnode), pointer, optional :: child
    character(len=*), optional, intent(in) :: modifier

    type(string) :: strBuffer
    type(fnode), pointer :: child2
    logical :: tReplace

    if (present(replace)) then
      tReplace = replace
    else
      tReplace = .false.
    end if

    call getAsString(value, strBuffer)
    call createChild_local(node, name, .false., tReplace, child2, &
        &value=char(strBuffer))
    if (present(child)) then
      child => child2
    end if
    if (present(modifier)) then
      call setAttribute(child2, attrModifier, modifier)
    end if

  end subroutine setChVal_logical


  !!* Sets the value (child) of a child with given name.
  !!* @param node    The node to investigate
  !!* @param name    Name of the child to look for
  !!* @param value   Value to set
  !!* @param replace Replace if child with same name already exists
  !!* @param child   Pointer to the child node (with the provided name)
  !!* @param modifier Optional modifier for the child
  subroutine setChVal_logicalR1(node, name, value, replace, child, modifier)
    type(fnode), pointer :: node
    character(len=*), intent(in) :: name
    logical, intent(in) :: value(:)
    logical, intent(in), optional :: replace
    type(fnode), pointer, optional :: child
    character(len=*), optional, intent(in) :: modifier

    type(string) :: strBuffer
    type(fnode), pointer :: child2
    logical :: tReplace

    if (present(replace)) then
      tReplace = replace
    else
      tReplace = .false.
    end if

    call getAsString(value, strBuffer)
    call createChild_local(node, name, .false., tReplace, child2, &
        &value=char(strBuffer))
    if (present(child)) then
      child => child2
    end if
    if (present(modifier)) then
      call setAttribute(child2, attrModifier, modifier)
    end if

  end subroutine setChVal_logicalR1



  !!* Writes the text representation of a node and its value to an xmlwriter.
  !!* @param xf    Xmlwriter stream
  !!* @param name  Name of the node
  !!* @param value Value of the node
  subroutine writeChVal_logical(xf, name, value)
    type(xmlf_t), intent(inout) :: xf
    character(len=*), intent(in) :: name
    logical, intent(in) :: value

    type(string) :: strBuffer

    call getAsString(value, strBuffer)
    call writeChild_local(xf, name, char(strBuffer))

  end subroutine writeChVal_logical

  !!* Writes the text representation of a node and its value to an xmlwriter.
  !!* @param xf    Xmlwriter stream
  !!* @param name  Name of the node
  !!* @param value Value of the node
  subroutine writeChVal_logicalR1(xf, name, value)
    type(xmlf_t), intent(inout) :: xf
    character(len=*), intent(in) :: name
    logical, intent(in) :: value(:)

    type(string) :: strBuffer

    call getAsString(value, strBuffer)
    call writeChild_local(xf, name, char(strBuffer))

  end subroutine writeChVal_logicalR1

  !!* Returns the text representation of the passed object
  !!* @param value     Value to represent
  !!* @param strBuffer Text representation on exit
  subroutine getAsString_logical(value, strBuffer)
    logical, intent(in) :: value
    type(string), intent(inout) :: strBuffer

    if (value) then
      strBuffer = LOGICAL_TRUE
    else
      strBuffer = LOGICAL_FALSE
    end if

  end subroutine getAsString_logical


  !!* Returns the text representation of the passed object
  !!* @param value     Value to represent
  !!* @param strBuffer Text representation on exit
  subroutine getAsString_logicalR1(value, strBuffer)
    logical, intent(in) :: value(:)
    type(string), intent(inout) :: strBuffer

    character(len=nCharLogical) :: buffer
    integer :: buffLen, len
    integer :: ii

    call resize_string(strBuffer, preAllocSize)
    len = 0
    do ii = 1, size(value)
      if (value(ii)) then
        write (buffer, *)LOGICAL_TRUE
      else
        write (buffer, *)LOGICAL_FALSE
      end if
      buffer = adjustl(buffer)
      buffLen = len_trim(buffer)
      len = len + buffLen
      if (len > lineLength) then
        call append_to_string(strBuffer, newline // trim(buffer))
        len = buffLen
      else
        call append_to_string(strBuffer, space // trim(buffer))
      end if
    end do

  end subroutine getAsString_logicalR1



  !!* Sets the value (child) of a child with given name.
  !!* @param node    The node to investigate
  !!* @param name    Name of the child to look for
  !!* @param value   Value to set
  !!* @param replace Replace if child with same name already exists
  !!* @param child   Pointer to the child node (with the provided name)
  !!* @param modifier Optional modifier for the child
  !!* @caveat This subroutines assumes, that a real can be represented as text
  !!*   with less than nCharReal characters.
  subroutine setChVal_real(node, name, value, replace, child, modifier)
    type(fnode), pointer :: node
    character(len=*), intent(in) :: name
    real(dp), intent(in) :: value
    logical, intent(in), optional :: replace
    type(fnode), pointer, optional :: child
    character(len=*), optional, intent(in) :: modifier

    type(string) :: strBuffer
    type(fnode), pointer :: child2
    logical :: tReplace

    if (present(replace)) then
      tReplace = replace
    else
      tReplace = .false.
    end if

    call getAsString(value, strBuffer)
    call createChild_local(node, name, .false., tReplace, child2, &
        &value=char(strBuffer))
    if (present(child)) then
      child => child2
    end if
    if (present(modifier)) then
      call setAttribute(child2, attrModifier, modifier)
    end if

  end subroutine setChVal_real



  !!* Writes the text representation of a node and its value to an xmlwriter.
  !!* @param xf    Xmlwriter stream
  !!* @param name  Name of the node
  !!* @param value Value of the node
  subroutine writeChVal_real(xf, name, value)
    type(xmlf_t), intent(inout) :: xf
    character(len=*), intent(in) :: name
    real(dp), intent(in) :: value

    type(string) :: strBuffer

    call getAsString(value, strBuffer)
    call writeChild_local(xf, name, char(strBuffer))

  end subroutine writeChVal_real



  !!* Returns the text representation of the passed object
  !!* @param value     Value to represent
  !!* @param strBuffer Text representation on exit
  subroutine getAsString_real(value, strBuffer)
    real(dp), intent(in) :: value
    type(string), intent(inout) :: strBuffer

    character(len=nCharReal) :: buffer

    write (buffer, *) value
    strBuffer = trim(adjustl(buffer))

  end subroutine getAsString_real



  !!* Sets the value (child) of a child with given name.
  !!* @param node    The node to investigate
  !!* @param name    Name of the child to look for
  !!* @param value   Value to set
  !!* @param replace Replace if child with same name already exists
  !!* @param child   Pointer to the child node (with the provided name)
  !!* @param modifier Optional modifier for the child
  !!* @caveat This subroutines assumes, that a real can be represented as text
  !!*   with less than nCharReal characters.
  subroutine setChVal_realR1(node, name, value, replace, child, modifier)
    type(fnode), pointer :: node
    character(len=*), intent(in) :: name
    real(dp), intent(in) :: value(:)
    logical, intent(in), optional :: replace
    type(fnode), pointer, optional :: child
    character(len=*), optional, intent(in) :: modifier

    type(string) :: strBuffer
    type(fnode), pointer :: child2
    logical :: tReplace

    if (present(replace)) then
      tReplace = replace
    else
      tReplace = .false.
    end if
    call getAsString(value, strBuffer)
    call createChild_local(node, name, .true., tReplace, child2, &
        &value=char(strBuffer))
    if (present(child)) then
      child => child2
    end if
    if (present(modifier)) then
      call setAttribute(child2, attrModifier, modifier)
    end if

  end subroutine setChVal_realR1



  !!* Writes the text representation of a node and its value to an xmlwriter.
  !!* @param xf    Xmlwriter stream
  !!* @param name  Name of the node
  !!* @param value Value of the node
  subroutine writeChVal_realR1(xf, name, value)
    type(xmlf_t), intent(inout) :: xf
    character(len=*), intent(in) :: name
    real(dp), intent(in) :: value(:)

    type(string) :: strBuffer

    call getAsString(value, strBuffer)
    call writeChild_local(xf, name, char(strBuffer))

  end subroutine writeChVal_realR1



  !!* Returns the text representation of the passed object
  !!* @param value     Value to represent
  !!* @param strBuffer Text representation on exit
  subroutine getAsString_realR1(value, strBuffer)
    real(dp), intent(in) :: value(:)
    type(string), intent(inout) :: strBuffer

    character(len=nCharReal) :: buffer
    integer :: buffLen, len
    integer :: ii

    call resize_string(strBuffer, preAllocSize)
    len = 0
    do ii = 1, size(value)
      write (buffer, *) value(ii)
      buffer = adjustl(buffer)
      buffLen = len_trim(buffer)
      len = len + buffLen
      if (len > lineLength) then
        call append_to_string(strBuffer, newline // trim(buffer))
        len = buffLen
      else
        call append_to_string(strBuffer, space // trim(buffer))
      end if
    end do

  end subroutine getAsString_realR1



  !!* Sets the value (child) of a child with given name.
  !!* @param node    The node to investigate
  !!* @param name    Name of the child to look for
  !!* @param value   Value to set
  !!* @param replace Replace if child with same name already exists
  !!* @caveat This subroutines assumes, that a real can be represented as text
  !!*   with less than nCharReal characters.
  !!* @param child   Pointer to the child node (with the provided name)
  !!* @param modifier Optional modifier for the child
  !!* @note This is just a wrapper around the rank one version, to make sure
  !!*   that two dimensional arrays are pretty printed. For higher ranked
  !!*   arrays the rank one version should be used with some reshaping before.
  subroutine setChVal_realR2(node, name, value, replace, child, modifier)
    type(fnode), pointer :: node
    character(len=*), intent(in) :: name
    real(dp), intent(in) :: value(:,:)
    logical, intent(in), optional :: replace
    type(fnode), pointer, optional :: child
    character(len=*), intent(in), optional :: modifier

    type(fnode), pointer :: child2
    type(string) :: strBuffer
    logical :: tReplace

    if (present(replace)) then
      tReplace = replace
    else
      tReplace = .false.
    end if

    call getAsString(value, strBuffer)
    call createChild_local(node, name, .true., tReplace, child2, &
        &value=char(strBuffer))
    if (present(child)) then
      child => child2
    end if
    if (present(modifier)) then
      call setAttribute(child2, attrModifier, modifier)
    end if

  end subroutine setChVal_realR2



  !!* Writes the text representation of a node and its value to an xmlwriter.
  !!* @param xf    Xmlwriter stream
  !!* @param name  Name of the node
  !!* @param value Value of the node
  subroutine writeChVal_realR2(xf, name, value)
    type(xmlf_t), intent(inout) :: xf
    character(len=*), intent(in) :: name
    real(dp), intent(in) :: value(:,:)

    type(string) :: strBuffer

    call getAsString(value, strBuffer)
    call writeChild_local(xf, name, char(strBuffer))

  end subroutine writeChVal_realR2



  !!* Returns the text representation of the passed object
  !!* @param value     Value to represent
  !!* @param strBuffer Text representation on exit
  subroutine getAsString_realR2(value, strBuffer)
    real(dp), intent(in) :: value(:,:)
    type(string), intent(inout) :: strBuffer

    character(len=nCharReal) :: buffer
    integer :: ii, jj

    call resize_string(strBuffer, preAllocSize)
    do ii = 1, size(value, dim=2)
      do jj = 1, size(value, dim=1)
        write (buffer, *) value(jj, ii)
        buffer = adjustl(buffer)
        call append_to_string(strBuffer, space // trim(buffer))
      end do
      call append_to_string(strBuffer, newline)
    end do

  end subroutine getAsString_realR2



  !!* Sets the value (child) of a child with given name.
  !!* @param node    The node to investigate
  !!* @param name    Name of the child to look for
  !!* @param value   Value to set
  !!* @param replace Replace if child with same name already exists
  !!* @param child   Pointer to the child node (with the provided name)
  !!* @param modifier Optional modifier for the child
  !!* @caveat This subroutines assumes, that an integer can be represented
  !!*   as text with less than nCharInt characters.
  subroutine setChVal_int(node, name, value, replace, child, modifier)
    type(fnode), pointer :: node
    character(len=*), intent(in) :: name
    integer, intent(in) :: value
    logical, intent(in), optional :: replace
    type(fnode), pointer, optional :: child
    character(len=*), optional, intent(in) :: modifier

    type(fnode), pointer :: child2
    type(string) :: strBuffer
    logical :: tReplace

    if (present(replace)) then
      tReplace = replace
    else
      tReplace = .false.
    end if
    call getAsString(value, strBuffer)
    call createChild_local(node, name, .false., tReplace, child2, &
        &value=char(strBuffer))
    if (present(child)) then
      child => child2
    end if
    if (present(modifier)) then
      call setAttribute(child2, attrModifier, modifier)
    end if

  end subroutine setChVal_int



  !!* Writes the text representation of a node and its value to an xmlwriter.
  !!* @param xf    Xmlwriter stream
  !!* @param name  Name of the node
  !!* @param value Value of the node
  subroutine writeChVal_int(xf, name, value)
    type(xmlf_t), intent(inout) :: xf
    character(len=*), intent(in) :: name
    integer, intent(in) :: value

    type(string) :: strBuffer

    call getAsString(value, strBuffer)
    call writeChild_local(xf, name, char(strBuffer))

  end subroutine writeChVal_int



  !!* Returns the text representation of the passed object
  !!* @param value     Value to represent
  !!* @param strBuffer Text representation on exit
  subroutine getAsString_int(value, strBuffer)
    integer, intent(in) :: value
    type(string), intent(inout) :: strBuffer

    character(len=nCharInt) :: buffer

    write (buffer, *) value
    strBuffer = trim(adjustl(buffer))

  end subroutine getAsString_int



  !!* Sets the value (child) of a child with given name.
  !!* @param node    The node to investigate
  !!* @param name    Name of the child to look for
  !!* @param value   Value to set
  !!* @param replace Replace if child with same name already exists
  !!* @param modifier Optional modifier for the child
  !!* @caveat This subroutines assumes, that an integer can be represented
  !!*   as text with less than nCharInt characters.
  subroutine setChVal_intR1(node, name, value, replace, child, modifier)
    type(fnode), pointer :: node
    character(len=*), intent(in) :: name
    integer, intent(in) :: value(:)
    logical, intent(in), optional :: replace
    type(fnode), pointer, optional :: child
    character(len=*), optional, intent(in) :: modifier

    type(fnode), pointer :: child2
    type(string) :: strBuffer
    logical :: tReplace

    if (present(replace)) then
      tReplace = replace
    else
      tReplace = .false.
    end if
    call getAsString(value, strBuffer)
    call createChild_local(node, name, .true., tReplace, child2, &
        &value=char(strBuffer))
    if (present(child)) then
      child => child2
    end if
    if (present(modifier)) then
      call setAttribute(child2, attrModifier, modifier)
    end if

  end subroutine setChVal_intR1



  !!* Writes the text representation of a node and its value to an xmlwriter.
  !!* @param xf    Xmlwriter stream
  !!* @param name  Name of the node
  !!* @param value Value of the node
  subroutine writeChVal_intR1(xf, name, value)
    type(xmlf_t), intent(inout) :: xf
    character(len=*), intent(in) :: name
    integer, intent(in) :: value(:)

    type(string) :: strBuffer

    call getAsString(value, strBuffer)
    call writeChild_local(xf, name, char(strBuffer))

  end subroutine writeChVal_intR1



  !!* Returns the text representation of the passed object
  !!* @param value     Value to represent
  !!* @param strBuffer Text representation on exit
  subroutine getAsString_intR1(value, strBuffer)
    integer, intent(in) :: value(:)
    type(string), intent(inout) :: strBuffer

    character(len=nCharInt) :: buffer
    integer :: buffLen, len
    integer :: ii

    call resize_string(strBuffer, preAllocSize)
    len = 0
    do ii = 1, size(value)
      write (buffer, *) value(ii)
      buffer = adjustl(buffer)
      buffLen = len_trim(buffer)
      len = len + buffLen
      if (len > lineLength) then
        call append_to_string(strBuffer, newline // trim(buffer))
        len = buffLen
      else
        call append_to_string(strBuffer, space // trim(buffer))
      end if
    end do

  end subroutine getAsString_intR1



  !!* Sets the value (child) of a child with given name.
  !!* @param node    The node to investigate
  !!* @param name    Name of the child to look for
  !!* @param value   Value to set
  !!* @param replace Replace if child with same name already exists
  !!* @param child   Pointer to the child node (with the provided name)
  !!* @caveat This subroutines assumes, that an integer can be represented as
  !!*   text with less than nCharInt characters.
  !!* @param modifier Optional modifier for the child
  !!* @note This is just a wrapper around the rank one version, to make sure
  !!*   that two dimensional arrays are pretty printed. For higher ranked
  !!*   arrays the rank one version should be used with some reshaping before.
  subroutine setChVal_intR2(node, name, value, replace, child, modifier)
    type(fnode), pointer :: node
    character(len=*), intent(in) :: name
    integer, intent(in) :: value(:,:)
    logical, intent(in), optional :: replace
    type(fnode), pointer, optional :: child
    character(len=*), optional, intent(in) :: modifier

    type(fnode), pointer :: child2
    type(string) :: strBuffer
    logical :: tReplace

    if (present(replace)) then
      tReplace = replace
    else
      tReplace = .false.
    end if
    call getAsString(value, strBuffer)
    call createChild_local(node, name, .true., tReplace, child2, &
        &value=char(strBuffer))
    if (present(child)) then
      child => child2
    end if
    if (present(modifier)) then
      call setAttribute(child2, attrModifier, modifier)
    end if

  end subroutine setChVal_intR2



  !!* Writes the text representation of a node and its value to an xmlwriter.
  !!* @param xf    Xmlwriter stream
  !!* @param name  Name of the node
  !!* @param value Value of the node
  subroutine writeChVal_intR2(xf, name, value)
    type(xmlf_t), intent(inout) :: xf
    character(len=*), intent(in) :: name
    integer, intent(in) :: value(:,:)

    type(string) :: strBuffer

    call getAsString(value, strBuffer)
    call writeChild_local(xf, name, char(strBuffer))

  end subroutine writeChVal_intR2



  !!* Returns the text representation of the passed object
  !!* @param value     Value to represent
  !!* @param strBuffer Text representation on exit
  subroutine getAsString_intR2(value, strBuffer)
    integer, intent(in) :: value(:,:)
    type(string), intent(inout) :: strBuffer

    character(len=nCharInt) :: buffer
    integer :: ii, jj

    call resize_string(strBuffer, preAllocSize)
    do ii = 1, size(value, dim=2)
      do jj = 1, size(value, dim=1)
        write (buffer, *) value(jj, ii)
        buffer = adjustl(buffer)
        call append_to_string(strBuffer, space // trim(buffer))
      end do
      call append_to_string(strBuffer, newline)
    end do

  end subroutine getAsString_intR2



  !!* Sets the value (child) of a child with given name.
  !!* @param node    The node to investigate
  !!* @param name    Name of the child to look for
  !!* @param value   Value to set
  !!* @param replace Replace if child with same name already exists
  !!* @param child   Pointer to the child node (with the provided name)
  !!* @param omitQuotes If qoutes around the string should be omitted
  !!* @param modifier Optional modifier for the child
  subroutine setChVal_char(node, name, value, replace, child, omitQuotes, &
      &modifier)
    type(fnode), pointer :: node
    character(len=*), intent(in) :: name
    character(len=*), intent(in) :: value
    logical, intent(in), optional :: replace
    type(fnode), pointer, optional :: child
    logical, intent(in), optional :: omitQuotes
    character(len=*), optional, intent(in) :: modifier

    type(fnode), pointer :: child2
    logical :: tReplace, tQuotes

    if (present(replace)) then
      tReplace = replace
    else
      tReplace = .false.
    end if
    if (present(omitQuotes)) then
      tQuotes = .not. omitQuotes
    else
      tQuotes = .true.
    end if
    if (tQuotes) then
      call createChild_local(node, name, .false., tReplace, child2, &
          &value='"'//value//'"')
    else
      call createChild_local(node, name, .false., tReplace, child2, value=value)
    end if

    if (present(child)) then
      child => child2
    end if
    if (present(modifier)) then
      call setAttribute(child2, attrModifier, modifier)
    end if

  end subroutine setChVal_char



  !!* Sets the value (child) of a child with given name.
  !!* @param node    The node to investigate
  !!* @param name    Name of the child to look for
  !!* @param value   Value to set
  !!* @param replace Replace if child with same name already exists
  !!* @param child   Pointer to the child node (with the provided name)
  !!* @param modifier Optional modifier for the child
  subroutine setChVal_charR1(node, name, value, replace, child, modifier)
    type(fnode), pointer :: node
    character(len=*), intent(in) :: name
    character(len=*), intent(in) :: value(:)
    logical, intent(in), optional :: replace
    type(fnode), pointer, optional :: child
    character(len=*), optional, intent(in) :: modifier

    type(string) :: strBuffer
    type(fnode), pointer :: child2
    logical :: tReplace

    if (present(replace)) then
      tReplace = replace
    else
      tReplace = .false.
    end if
    call getAsString(value, strBuffer)
    call createChild_local(node, name, .true., tReplace, child2, &
        &value=char(strBuffer))
    if (present(child)) then
      child => child2
    end if
    if (present(modifier)) then
      call setAttribute(child2, attrModifier, modifier)
    end if

  end subroutine setChVal_charR1



  !!* Writes the text representation of a node and its value to an xmlwriter.
  !!* @param xf    Xmlwriter stream
  !!* @param name  Name of the node
  !!* @param value Value of the node
  subroutine writeChVal_charR1(xf, name, value)
    type(xmlf_t), intent(inout) :: xf
    character(len=*), intent(in) :: name
    character(len=*), intent(in) :: value(:)

    type(string) :: strBuffer

    call getAsString(value, strBuffer)
    call writeChild_local(xf, name, char(strBuffer))

  end subroutine writeChVal_charR1



  !!* Returns the text representation of the passed object
  !!* @param value     Value to represent
  !!* @param strBuffer Text representation on exit
  subroutine getAsString_charR1(value, strBuffer)
    character(len=*), intent(in) :: value(:)
    type(string), intent(inout) :: strBuffer

    integer :: buffLen, len
    integer :: ii

    call resize_string(strBuffer, preAllocSize)
    len = 0
    do ii = 1, size(value)
      buffLen = len_trim(value(ii))
      len = len + buffLen
      if (len > lineLength) then
        call append_to_string(strBuffer, newline // '"'//trim(value(ii))//'"')
        len = buffLen
      else
        call append_to_string(strBuffer, space // '"'//trim(value(ii))//'"')
      end if
    end do

  end subroutine getAsString_charR1



  !!* Sets the value (child) of a child with given name.
  !!* @param node      The node to investigate
  !!* @param name      Name of the child to look for
  !!* @param intValue  Value for the integers
  !!* @param realValue Value for the reals
  !!* @param replace Replace if child with same name already exists
  !!* @param child     Pointer to the child node (with the provided name)
  !!* @param modifier Optional modifier for the child
  subroutine setChVal_intR2RealR2(node, name, intValue, realValue, &
      &replace, child, modifier)
    type(fnode), pointer :: node
    character(len=*), intent(in) :: name
    integer, intent(in) :: intValue(:,:)
    real(dp), intent(in) :: realValue(:,:)
    logical, intent(in), optional :: replace
    type(fnode), pointer, optional :: child
    character(len=*), optional, intent(in) :: modifier

    type(fnode), pointer :: child2
    type(string) :: strBuffer
    logical :: tReplace

    if (present(replace)) then
      tReplace = replace
    else
      tReplace = .false.
    end if
    call getAsString(intValue, realValue, strBuffer)
    call createChild_local(node, name, .true., tReplace, child2, &
        &value=char(strBuffer))
    if (present(child)) then
      child => child2
    end if
    if (present(modifier)) then
      call setAttribute(child2, attrModifier, modifier)
    end if

  end subroutine setChVal_intR2RealR2



  !!* Writes the text representation of a node and its value to an xmlwriter.
  !!* @param xf    Xmlwriter stream
  !!* @param name  Name of the node
  !!* @param value Value of the node
  subroutine writeChVal_intR2RealR2(xf, name, intValue, realValue)
    type(xmlf_t), intent(inout) :: xf
    character(len=*), intent(in) :: name
    integer, intent(in) :: intValue(:,:)
    real(dp), intent(in) :: realValue(:,:)

    type(string) :: strBuffer

    call getAsString(intValue, realValue, strBuffer)
    call writeChild_local(xf, name, char(strBuffer))

  end subroutine writeChVal_intR2RealR2



  !!* Returns the text representation of the passed object
  !!* @param intValue
  !!* @param realValue
  !!* @param strBuffer Text representation on exit
  subroutine getAsString_intR2RealR2(intValue, realValue, strBuffer)
    integer, intent(in) :: intValue(:,:)
    real(dp), intent(in) :: realValue(:,:)
    type(string), intent(inout) :: strBuffer

    character(len=100) :: buffer
    integer :: nRow, nCol1, nCol2
    integer :: ii, jj

    nRow = size(intValue, dim=2)
    @:ASSERT(size(realValue, dim=2) == nRow)

    nCol1 = size(intValue, dim=1)
    nCol2 = size(realValue, dim=1)
    call resize_string(strBuffer, preAllocSize)
    do ii = 1, nRow
      do jj = 1, nCol1
        write (buffer, *) intValue(jj, ii)
        buffer = adjustl(buffer)
        call append_to_string(strBuffer, space // trim(buffer))
      end do
      do jj = 1, nCol2
        write (buffer, *) realValue(jj, ii)
        buffer = adjustl(buffer)
        call append_to_string(strBuffer, space // trim(buffer))
      end do
      call append_to_string(strBuffer, newline)
    end do

  end subroutine getAsString_intR2RealR2



  !!* Sets the value (child) of a child with given name.
  !!* @param node      The node to investigate
  !!* @param name      Name of the child to look for
  !!* @param charValue Value for the characters
  !!* @param intValue  Value for the integers
  !!* @param realValue Value for the reals
  !!* @param replace Replace if child with same name already exists
  !!* @param child     Pointer to the child node (with the provided name)
  !!* @param modifier Optional modifier for the child
  subroutine setChVal_charR1IntR2RealR2(node, name, charValue, intValue, &
      &realValue, replace, child, modifier)
    type(fnode), pointer :: node
    character(len=*), intent(in) :: name
    character(len=*), intent(in) :: charValue(:)
    integer, intent(in) :: intValue(:,:)
    real(dp), intent(in) :: realValue(:,:)
    logical, intent(in), optional :: replace
    type(fnode), pointer, optional :: child
    character(len=*), optional, intent(in) :: modifier

    type(fnode), pointer :: child2
    type(string) :: strBuffer
    logical :: tReplace

    if (present(replace)) then
      tReplace = replace
    else
      tReplace = .false.
    end if
    call getAsString(charValue, intValue, realValue, strBuffer)
    call createChild_local(node, name, .true., tReplace, child2, &
        &value=char(strBuffer))
    if (present(child)) then
      child => child2
    end if
    if (present(modifier)) then
      call setAttribute(child2, attrModifier, modifier)
    end if

  end subroutine setChVal_charR1IntR2RealR2



  !!* Writes the text representation of a node and its value to an xmlwriter.
  !!* @param xf    Xmlwriter stream
  !!* @param name  Name of the node
  !!* @param value Value of the node
  subroutine writeChVal_charR1IntR2RealR2(xf, name, charValue, intValue, &
      &realValue)
    type(xmlf_t), intent(inout) :: xf
    character(len=*), intent(in) :: name
    character(len=*), intent(in) :: charValue(:)
    integer, intent(in) :: intValue(:,:)
    real(dp), intent(in) :: realValue(:,:)

    type(string) :: strBuffer

    call getAsString(charValue, intValue, realValue, strBuffer)
    call writeChild_local(xf, name, char(strBuffer))

  end subroutine writeChVal_charR1IntR2RealR2



  !!* Returns the text representation of the passed object
  !!* @param charValue
  !!* @param intValue
  !!* @param realValue
  !!* @param strBuffer Text representation on exit
  subroutine getAsString_charR1IntR2RealR2(charValue, intValue, realValue, &
      &strBuffer)
    character(len=*), intent(in) :: charValue(:)
    integer, intent(in) :: intValue(:,:)
    real(dp), intent(in) :: realValue(:,:)
    type(string), intent(inout) :: strBuffer

    character(len=100) :: buffer
    integer :: nRow, nCol1, nCol2
    integer :: ii, jj

    nRow = size(charValue)
    @:ASSERT(size(intValue, dim=2) == nRow)
    @:ASSERT(size(realValue, dim=2) == nRow)

    nCol1 = size(intValue, dim=1)
    nCol2 = size(realValue, dim=1)
    call resize_string(strBuffer, preAllocSize)
    do ii = 1, nRow
      call append_to_string(strBuffer, charValue(ii))
      do jj = 1, nCol1
        write (buffer, *) intValue(jj, ii)
        buffer = adjustl(buffer)
        call append_to_string(strBuffer, space // trim(buffer))
      end do
      do jj = 1, nCol2
        write (buffer, *) realValue(jj, ii)
        buffer = adjustl(buffer)
        call append_to_string(strBuffer, space // trim(buffer))
      end do
      call append_to_string(strBuffer, newline)
    end do

  end subroutine getAsString_charR1IntR2RealR2



  !!* Sets the value (child) of a child with given name.
  !!* @param node    The node to investigate
  !!* @param name    Name of the child to look for
  !!* @param value   Value to set
  !!* @param replace Replace if child with same name already exists
  !!* @param child   Pointer to the child node (with the provided name)
  !!* @param modifier Optional modifier for the child
  !!* @param list    If created child should be marked as a list.
  subroutine setChVal_node(node, name, value, replace, child, modifier, &
      &list)
    type(fnode), pointer :: node
    character(len=*), intent(in) :: name
    type(fnode), pointer :: value
    logical, intent(in), optional :: replace
    type(fnode), pointer, optional :: child
    character(len=*), optional, intent(in) :: modifier
    logical, optional, intent(in) :: list

    type(fnode), pointer :: child2, dummy
    logical :: tReplace, tList

    if (present(replace)) then
      tReplace = replace
    else
      tReplace = .false.
    end if
    if (present(list)) then
      tList = list
    else
      tList = .false.
    end if
    call createChild_local(node, name, tList, tReplace, child2)
    if (associated(value)) then
      dummy => appendChild(child2, value)
    end if
    if (present(child)) then
      child => child2
    end if
    if (present(modifier)) then
      call setAttribute(child2, attrModifier, modifier)
    end if

  end subroutine setChVal_node



  !!* Working horse for the setChildValue routines
  !!* @param node    The node to investigate
  !!* @param name    Name of the child to create
  !!* @param list    True, if child should be signed as a list
  !!* @param replace Replace if child with same name already exists
  !!* @param child   Pointer to the created child on return
  !!* @param value   Value to set (if empty, no child is appended to the created
  !!*   child)
  !!* @note If the empty string is provided as child name, no child is created,
  !!*   and the current node is replace instead. The pointer "node" becames
  !!*   associated with the new node, since the old instance will be destroyed.
  subroutine createChild_local(node, name, list, replace, child, value)
    type(fnode), pointer :: node
    character(len=*), intent(in) :: name
    logical, intent(in) :: list
    logical, intent(in) :: replace
    type(fnode), pointer :: child
    character(len=*), intent(in), optional :: value

    type(fnode), pointer :: parent, oldChild, child2, text, dummy
    character(len=len(name)) :: loName
    type(string) :: newName, parentname

    if (replace) then
      if (len(name) == 0) then
        call getNodeHSDName(node, newName)
        parent => getParentNode(node)
        oldChild => node
        child2 => createElement(tolower(char(newName)))
        node => child2
      else
        newName = name
        parent => node
        loName = tolower(name)
        oldChild => getFirstChildByName(node, loName)
        child2 => createElement(loName)
      end if
    else
      newName = name
      parent => node
      oldChild => null()
      child2 => createElement(tolower(name))
    end if

    ! If parent is a text mode, no subnodes should be allowed.
    dummy => getFirstChild(parent)
    if (associated(dummy)) then
      call getNodeName(dummy, parentname)
      if (char(parentname) == textNodeName) then
        call detailedError(node, "Node contains superfluous free text: '"&
            & // trim(dummy%nodeValue) // "'")
      end if
    end if

    if (associated(oldChild)) then
      dummy => replaceChild(parent, child2, oldChild)
      call destroyNode(oldChild)
    else
      dummy => appendChild(parent, child2)
    end if

    if (len(newName) > 0) then
      call setAttribute(child2, attrName, char(newName))
    end if
    if (list) then
      call setAttribute(child2, attrList, "")
    end if

    child => child2
    call setAttribute(child, attrProcessed, "")
    if (present(value)) then
      text => createTextNode(value)
      dummy => appendChild(child, text)
    end if

  end subroutine createChild_local



  subroutine writeChild_local(xf, name, value)
    type(xmlf_t), intent(inout) :: xf
    character(len=*), intent(in) :: name
    character(len=*), intent(in) :: value

    call xml_NewElement(xf, name)
    call xml_AddPCData(xf, value)
    call xml_EndElement(xf, name)

  end subroutine writeChild_local



  !!* Creates a child with the given name
  !!* @param node    Node to append the child to
  !!* @param name    Name of the child node to append
  !!* @param child   Contains the pointer to the added child node on return
  !!* @param replace If an already existing child with the same name should be
  !!*   replaced
  !!* @param list    If child should be signed as a list tag
  !!* @param modifier Optional modifier for the child
  subroutine setChild(node, name, child, replace, list, modifier)
    type(fnode), pointer :: node
    character(len=*), intent(in) :: name
    type(fnode), pointer :: child
    logical, intent(in), optional :: replace
    logical, intent(in), optional :: list
    character(len=*), optional, intent(in) :: modifier

    logical :: tReplace, tList
    type(fnode), pointer :: dummy

    if (present(replace)) then
      tReplace = replace
    else
      tReplace = .false.
    end if
    if (present(list)) then
      tList = list
    else
      tList = .false.
    end if

    child => getFirstChildByName(node, tolower(name))
    if (associated(child)) then
      if (tReplace) then
        dummy => removeChild(node, child)
        call destroyNode(child)
      else
        call detailedError(node, MSG_EXISTING_CHILD // name)
      end if
    end if
    child => createElement(tolower(name))
    dummy => appendChild(node, child)
    call setAttribute(child, attrName, name)
    call setAttribute(child, attrProcessed, "")
    if (tList) then
      call setAttribute(child, attrList, "")
    end if
    if (present(modifier)) then
      call setAttribute(child, attrModifier, modifier)
    end if

  end subroutine setChild


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!! Miscancellous routines
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !!* Returns the content of the first TEXT_NODE child of a given node or
  !!* empty string, if such node does not exist.
  !!* @param node The node to investigate.
  !!* @return String representation of the TEXT_NODE.
  !!* @note If the document tree is normalized, every node has only one
  !!*   TEXT_NODE child.
  subroutine getFirstTextChild(node, str)
    type(fnode), pointer :: node
    type(string), intent(out) :: str

    type(fnode), pointer :: child

    child => getFirstChild(node)
    if (.not. associated(child)) then
      str = ""
    elseif (getNodeType(child) /= TEXT_NODE) then
      call detailedError(child, "Invalid node type.")
    else
      call getNodeValue(child, str)
    end if

  end subroutine getFirstTextChild



  !!* Checks if error flag signalizes an error. If yes, raises error.
  !!* @param node Node which the error flag was set for
  !!* @param iErr Content of the error flag.
  !!* @param msg  Message to print, if error occured
  subroutine checkError(node, iErr, msg)
    type(fnode), pointer :: node
    integer, intent(in) :: iErr
    character(len=*), intent(in) :: msg

    if (iErr == TOKEN_ERROR) then
      call detailedError(node, msg)
    elseif (iErr == TOKEN_EOS) then
      call detailedError(node, "Unexpected end of data")
    end if

  end subroutine checkError


  !!* Issues an error, if the string from a given position contains
  !!* non-whitespace characters.
  !!* @param node  Node which is being processed (for error message)
  !!* @param str String content of the child.
  !!* @param start Starting position, after which the string should not contain
  !!*   any whitespace characters.
  subroutine checkNoData(node, str, start)
    type(fnode), pointer :: node
    character(len=*), intent(in) :: str
    integer, intent(in) :: start

    if (complementaryScan(str(start:), whiteSpaces) > 0) then
      call detailedError(node, "Superfluous data found.")
    end if

  end subroutine checkNoData



  !!* Prints detailed error, including line number and path
  !!* @param node Node where the error occured.
  !!* @param msg  Message to print
  subroutine detailedError(node, msg)
    type(fnode), pointer :: node
    character(len=*), intent(in) :: msg

    type(string) :: str

    str = msg
    call appendPathAndLine(node, str)
    call error(char(str) // newline)

  end subroutine detailedError



  !!* Prints detailed warning, including line number and path
  !!* @param node Node where the error occured.
  !!* @param msg  Message to print
  subroutine detailedWarning(node, msg)
    type(fnode), pointer :: node
    character(len=*), intent(in) :: msg

    type(string) :: str

    str = msg
    call appendPathAndLine(node, str)
    call warning(char(str) // newline)

  end subroutine detailedWarning



  !!* Appends path and line information to a string.
  !!* @param node Node, for which path and line should be added
  !!* @param str  String prepending the path and line information
  subroutine appendPathAndLine(node, str)
    type(fnode), pointer :: node
    type(string), intent(inout) :: str

    type(string) :: str2, str3

    call append_to_string(str, newline // "Path: ")
    call getHSDPath(node, str2, excludeRoot=.true.)
    call append_to_string(str, str2)
    call getAttribute(node, attrStart, str2)
    call getAttribute(node, attrEnd, str3)
    if (len(str2) /= 0) then
      call append_to_string(str, newline // "Line: ")
      call append_to_string(str, str2)
      if (len(str3) /= 0) then
        call append_to_string(str, "-")
        call append_to_string(str, str3)
      end if
    end if
    call getAttribute(node, attrFile, str2)
    if (len(str2) /= 0) then
      call append_to_string(str, " (File: ")
      call append_to_string(str, str2)
      call append_to_string(str, ")")
    end if

  end subroutine appendPathAndLine


end module hsdutils
