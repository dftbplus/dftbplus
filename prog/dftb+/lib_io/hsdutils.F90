!--------------------------------------------------------------------------------------------------!
!  DFTB+: general package for performing fast atomistic simulations                                !
!  Copyright (C) 2006 - 2020  DFTB+ developers group                                               !
!                                                                                                  !
!  See the LICENSE file for terms of usage and distribution.                                       !
!--------------------------------------------------------------------------------------------------!

#:include 'common.fypp'

!> Contains high level functions for converting the values in a XML/HSD DOM-tree to Fortran
!> intrinsic types.
!> Todo: Some more routines for complex numbers?
module dftbp_hsdutils
  use dftbp_assert
  use dftbp_xmlf90
  use dftbp_tokenreader
  use dftbp_hsdparser
  use dftbp_xmlutils
  use dftbp_charmanip
  use dftbp_message
  use dftbp_linkedlist
  use dftbp_accuracy
  implicit none
  private

  public :: checkError, detailedError, detailedWarning
  public :: getFirstTextChild, getChildValue, setChildValue
  public :: writeChildValue, getAsString
  public :: convAtomRangeToInt, convRangeToInt, appendPathAndLine
  public :: getChild, getChildren, setChild
  public :: attrProcessed


  !> Returns the value (the child) of a child node identified by its name.
  !>
  !> These routines investigate the provided node and look for a child with the supplied name. If
  !> this child found, its child (which should be a single text node or a usual node if the value
  !> argument is of type node) is returned as value converted to the appropriate type. If the child
  !> is not found, an error is raised, unless a default value was specified.In that case, a child is
  !> created with the provided name and is appended to the node. Furthermore a text node containing
  !> the string converted default value is appended to the child node. If default value is provided,
  !> it must be also indicated, if the created child is only allowed to have one further child or
  !> not. (This corresponds to an assignment with '=' in the HSD input.) If the child (identified by
  !> the provided name) is allowed to have a modifier, an argument for the modifier must be provided
  !> to contain the the parsed value on return. If the argument for the modifier is missing, but a
  !> modifier is found, the program raises an error. The pointer to the found (or created) child can
  !> be queried through an appropriate argument. If the name of the child to look for is an empty
  !> string, the passed node itself is treated as if it would be the child, which had been found.
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
  end interface getChildValue


  !> Sets the value (the child) of a child node identified by its name
  !>
  !> Those functions are the inverse of the getChildValue functions. They create a child with the
  !> provided name and append to that child a text node (or a normal node, if the provided value is
  !> of type node) containing the provided value. It must be indicated, if the created child is
  !> allowed to have only one single further child. If a child with the specified name already
  !> exists, the program raises an error, unless replacement flag is set on .true.. In that case,
  !> the the existing child is replaced. If the name of the child is the empty string, the current
  !> node is treated as if it would be the child, which had been found.
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
  end interface setChildValue


  !> Writes a child and its value to an xml-write stream
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
  end interface writeChildValue


  !> Returns a string representation of an object
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
  end interface getAsString


  !> Error messages
  character(len=*), parameter :: MSG_MISSING_FIELD = "Missing child: "
  character(len=*), parameter :: MSG_EXISTING_CHILD = "Already existing child: "
  character(len=*), parameter :: MSG_NOMODIFIER = "Entity is not allowed to &
      &have a modifier"
  character(len=*), parameter :: MSG_MISSING_VALUES = "Not enough values &
      &provided."


  !> Length of a line (for wrapping long lines when writing values)
  integer, parameter :: lineLength = 80


  !> Maximal number of characters needed to represent an integer
  integer, parameter :: nCharInt = 50


  !> Maximal number of characters needed to represent a real number
  integer, parameter :: nCharReal = 50


  !> Maximal number of characters needed to represent a logical value
  integer, parameter :: nCharLogical = 4


  !> Attribute signals that a tag was processed
  character(len=*), parameter :: attrProcessed = "proc"


  !> Preallocateated size for temporary buffer strings
  integer, parameter :: preAllocSize = 1024

contains


  !> Returns the value (the child) of a child node as logical.
  subroutine getChVal_logical(node, name, variableValue, default, modifier, child)

    !> The node to investigate.
    type(fnode), pointer :: node

    !> Name of the child to look for
    character(len=*), intent(in) :: name

    !> Value on return
    logical, intent(out) :: variableValue

    !> Default value for the child, if child is not found
    logical, intent(in), optional :: default

    !> Modifier of the child on return
    type(string), intent(inout), optional :: modifier

    !> Pointer to the child node (with the spec. name) on return
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
      call getNextToken(char(text), variableValue, iStart, iErr)
      call checkError(child2, iErr, "Invalid logical value")
      call checkNoData(child2, char(text), iStart)
      call setAttribute(child2, attrProcessed, "")
    elseif (present(default)) then
      variableValue = default
      if (present(modifier)) then
        modifier = ""
      end if
      call setChildValue(node, name, variableValue, .false., child=child2)
    else
      call detailedError(node, MSG_MISSING_FIELD // name)
    end if
    if (present(child)) then
      child => child2
    end if

  end subroutine getChVal_logical


  !> Returns the value (the child) of a child node as logical.
  subroutine getChVal_logicalR1(node, name, variableValue, default, nItem, modifier, child)

    !> The node to investigate.
    type(fnode), pointer :: node

    !> Name of the child to look for
    character(len=*), intent(in) :: name

    !> Value on return
    logical, intent(out) :: variableValue(:)

    !> Default value for the child, if child is not found
    logical, intent(in), optional :: default(:)

    !> Nr. of read items. If this argument is not passed, and the nr. of read items is less than the
    !> size of the array, the subroutine raises an error.
    integer, intent(out), optional :: nItem

    !> Modifier of the child on return
    type(string), intent(inout), optional :: modifier

    !> Pointer to the child node (with the spec. name) on return
    type(fnode), pointer, optional :: child

    type(string) :: text, modif
    integer :: iStart, iErr, nReadItem
    type(fnode), pointer :: child2

    @:ASSERT(associated(node))
  #:block DEBUG_CODE
    if (present(default)) then
      @:ASSERT(all(shape(default) == shape(variableValue)))
    end if
  #:endblock DEBUG_CODE

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
      call getNextToken(char(text), variableValue, iStart, iErr, nReadItem)
      call checkError(child2, iErr, "Invalid logical value")
      call checkNoData(child2, char(text), iStart)
      if (present(nItem)) then
        nItem = nReadItem
      elseif (nReadItem /= size(variableValue)) then
        call detailedError(node, MSG_MISSING_VALUES)
      end if
      call setAttribute(child2, attrProcessed, "")
    elseif (present(default)) then
      variableValue = default
      if (present(nItem)) then
        nItem = size(default)
      end if
      if (present(modifier)) then
        modifier = ""
      end if
      call setChildValue(node, name, variableValue, .false., child=child2)
    else
      call detailedError(node, MSG_MISSING_FIELD // name)
    end if
    call setAttribute(child2, attrProcessed, "")
    if (present(child)) then
      child => child2
    end if

  end subroutine getChVal_logicalR1


  !> Returns the value (the child) of a child node as string.
  subroutine getChVal_string(node, name, variableValue, default, modifier, child, multiple)

    !> The node to investigate.
    type(fnode), pointer :: node

    !> Name of the child to look for
    character(len=*), intent(in) :: name

    !> Value on return
    type(string), intent(inout) :: variableValue

    !> Default value for the child, if child is not found
    character(len=*), intent(in), optional :: default

    !> Modifier of the child on return
    type(string), intent(inout), optional :: modifier

    !> Pointer to the child node (with the spec. name) on return
    type(fnode), pointer, optional :: child

    !> If true, string contains as many tokens as possible, not just one (with spaces between the
    !> tokens).
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
        variableValue = unquote(trim(adjustl(char(text))))
      else
        iStart = 1
        call getNextToken(char(text), variableValue, iStart, iErr)
        call checkError(child2, iErr, "Invalid string value")
        call checkNoData(child2, char(text), iStart)
      end if
      call setAttribute(child2, attrProcessed, "")
    elseif (present(default)) then
      variableValue = default
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


  !> Returns the value (the child) of a child node as real.
  subroutine getChVal_real(node, name, variableValue, default, modifier, child)

    !> The node to investigate.
    type(fnode), pointer :: node

    !> Name of the child to look for
    character(len=*), intent(in) :: name

    !> Value on return
    real(dp), intent(out) :: variableValue

    !> Default value for the child, if child is not found
    real(dp), intent(in), optional :: default

    !> Modifier of the child on return
    type(string), intent(inout), optional :: modifier

    !> Pointer to the child node (with the spec. name) on return
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
      call getNextToken(char(text), variableValue, iStart, iErr)
      call checkError(child2, iErr, "Invalid real value")
      call checkNoData(child2, char(text), iStart)
      call setAttribute(child2, attrProcessed, "")
    elseif (present(default)) then
      variableValue = default
      if (present(modifier)) then
        modifier = ""
      end if
      call setChildValue(node, name, variableValue, .false., child=child2)
    else
      call detailedError(node, MSG_MISSING_FIELD // name)
    end if
    call setAttribute(child2, attrProcessed, "")
    if (present(child)) then
      child => child2
    end if

  end subroutine getChVal_real


  !> Returns the value (the child) of a child node as a rank one real array.
  subroutine getChVal_realR1(node, name, variableValue, default, nItem, modifier, child)

    !> The node to investigate.
    type(fnode), pointer :: node

    !> Name of the child to look for
    character(len=*), intent(in) :: name

    !> Value on return
    real(dp), intent(out) :: variableValue(:)

    !> Default value for the child, if child is not found
    real(dp), intent(in), optional :: default(:)

    !> Nr. of read items. If this argument is not passed, and the nr. of read items is less than the
    !> size of the array, the subroutine raises an error.
    integer, intent(out), optional :: nItem

    !> Modifier of the child on return
    type(string), intent(inout), optional :: modifier

    !> Pointer to the child node (with the spec. name) on return
    type(fnode), pointer, optional :: child

    type(string) :: text, modif
    integer :: iStart, iErr, nReadItem
    type(fnode), pointer :: child2

    @:ASSERT(associated(node))
  #:block DEBUG_CODE
    if (present(default)) then
      @:ASSERT(all(shape(default) == shape(variableValue)))
    end if
  #:endblock DEBUG_CODE

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
      call getNextToken(char(text), variableValue, iStart, iErr, nReadItem)
      call checkError(child2, iErr, "Invalid real value")
      call checkNoData(child2, char(text), iStart)
      if (present(nItem)) then
        nItem = nReadItem
      elseif (nReadItem /= size(variableValue)) then
        call detailedError(node, MSG_MISSING_VALUES)
      end if
      call setAttribute(child2, attrProcessed, "")
    elseif (present(default)) then
      variableValue = default
      if (present(nItem)) then
        nItem = size(default)
      end if
      if (present(modifier)) then
        modifier = ""
      end if
      call setChildValue(node, name, variableValue, .false., child=child2)
    else
      call detailedError(node, MSG_MISSING_FIELD // name)
    end if
    call setAttribute(child2, attrProcessed, "")
    if (present(child)) then
      child => child2
    end if

  end subroutine getChVal_realR1


  !> Returns the value (the child) of a child node as a rank two real array.
  !>
  !> This is just a wrapper around the rank one version, to make sure that two dimensional arrays
  !> are pretty printed. For higher ranked arrays the rank one version should be used with some
  !> reshaping after.
  subroutine getChVal_realR2(node, name, variableValue, default, nItem, modifier, child)

    !> The node to investigate.
    type(fnode), pointer :: node

    !> Name of the child to look for
    character(len=*), intent(in) :: name

    !> Value on return
    real(dp), intent(out) :: variableValue(:,:)

    !> Default value for the child, if child is not found
    real(dp), intent(in), optional :: default(:,:)

    !> Nr. of read items. If this argument is not passed, and the nr. of read items is less than the
    !> size of the array, the subroutine raises an error.
    integer, intent(out), optional :: nItem

    !> Modifier of the child on return
    type(string), intent(inout), optional :: modifier

    !> Pointer to the child node (with the spec. name) on return
    type(fnode), pointer, optional :: child

    real(dp) :: buffer(size(variableValue))
    integer :: nReadItem
    type(string) :: modif
    type(fnode), pointer :: child2

    @:ASSERT(associated(node))
  #:block DEBUG_CODE
    if (present(default)) then
      @:ASSERT(all(shape(default) == shape(variableValue)))
    end if
  #:endblock DEBUG_CODE

    nReadItem = 0
    variableValue = 0.0_dp
    if (present(default)) then
      call getChildValue(node, name, buffer, reshape(default, shape(buffer)), &
          &nReadItem, modifier=modif, child=child2)
    else
      call getChildValue(node, name, buffer, nItem=nReadItem, modifier=modif, &
          &child=child2)
    end if
    if (present(nItem)) then
      nItem = nReadItem
    elseif (nReadItem /= size(variableValue)) then
      call detailedError(node, MSG_MISSING_VALUES)
    end if
    if (present(modifier)) then
      modifier = modif
    elseif (len(modif) > 0) then
      call detailedError(child2, MSG_NOMODIFIER)
    end if
    variableValue(:,:) = reshape(buffer, shape(variableValue))
    if (present(child)) then
      child => child2
    end if

  end subroutine getChVal_realR2


  !> Returns the value (the child) of a child node as integer.
  subroutine getChVal_int(node, name, variableValue, default, modifier, child)

    !> The node to investigate.
    type(fnode), pointer :: node

    !> Name of the child to look for
    character(len=*), intent(in) :: name

    !> Value on return
    integer, intent(out) :: variableValue

    !> Default value for the child, if child is not found
    integer, intent(in), optional :: default

    !> Modifier of the child on return
    type(string), intent(inout), optional :: modifier

    !> Pointer to the child node (with the spec. name) on return
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
      call getNextToken(char(text), variableValue, iStart, iErr)
      call checkError(child2, iErr, "Invalid integer variableValue")
      call checkNoData(child2, char(text), iStart)
      call setAttribute(child2, attrProcessed, "")
    elseif (present(default)) then
      variableValue = default
      if (present(modifier)) then
        modifier = ""
      end if
      call setChildValue(node, name, variableValue, .false., child=child2)
    else
      call detailedError(node, MSG_MISSING_FIELD // name)
    end if
    if (present(child)) then
      child => child2
    end if

  end subroutine getChVal_int


  !> Returns the value (the child) of a child node as a rank one integer array.
  subroutine getChVal_intR1(node, name, variableValue, default, nItem, modifier, child)

    !> The node to investigate.
    type(fnode), pointer :: node

    !> Name of the child to look for
    character(len=*), intent(in) :: name

    !> Value on return
    integer, intent(out) :: variableValue(:)

    !> Default value for the child, if child is not found
    integer, intent(in), optional :: default(:)

    !> Nr. of read items. If this argument is not passed, and the nr. of read items is less than the
    !> size of the array, the subroutine raises an error.
    integer, intent(out), optional :: nItem

    !> Modifier of the child on return
    type(string), intent(inout), optional :: modifier

    !> Pointer to the child node (with the spec. name) on return
    type(fnode), pointer, optional :: child

    type(string) :: text, modif
    integer :: iStart, iErr, nReadItem
    type(fnode), pointer :: child2

    @:ASSERT(associated(node))
  #:block DEBUG_CODE
    if (present(default)) then
      @:ASSERT(all(shape(default) == shape(variableValue)))
    end if
  #:endblock DEBUG_CODE

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
      call getNextToken(char(text), variableValue, iStart, iErr, nReadItem)
      call checkError(child2, iErr, "Invalid integer value")
      call checkNoData(child2, char(text), iStart)
      if (present(nItem)) then
        nItem = nReadItem
      elseif (nReadItem /= size(variableValue)) then
        call detailedError(node, MSG_MISSING_VALUES)
      end if
      call setAttribute(child2, attrProcessed, "")
    elseif (present(default)) then
      variableValue = default
      if (present(nItem)) then
        nItem = size(default)
      end if
      if (present(modifier)) then
        modifier = ""
      end if
      call setChildValue(node, name, variableValue, .false., child=child2)
    else
      call detailedError(node, MSG_MISSING_FIELD // name)
    end if
    if (present(child)) then
      child => child2
    end if

  end subroutine getChVal_intR1


  !> Returns the value (the child) of a child node as a rank two integer array.
  !>
  !> This is just a wrapper around the rank one version, to make sure that two dimensional arrays
  !> are pretty printed. For higher ranked arrays the rank one version should be used with some
  !> reshaping after.
  subroutine getChVal_intR2(node, name, variableValue, default, nItem, modifier, child)

    !> The node to investigate.
    type(fnode), pointer :: node

    !> Name of the child to look for
    character(len=*), intent(in) :: name

    !> Value on return
    integer, intent(out) :: variableValue(:,:)

    !> Default value for the child, if child is not found
    integer, intent(in), optional :: default(:,:)

    !> Nr. of read items. If this argument is not passed, and the nr. of read items is less than the
    !> size of the array, the subroutine raises an error.
    integer, intent(out), optional :: nItem

    !> Modifier of the child on return
    type(string), intent(inout), optional :: modifier

    !> Pointer to the child node (with the spec. name) on return
    type(fnode), pointer, optional :: child

    integer :: buffer(size(variableValue))
    integer :: nReadItem
    type(string) :: modif
    type(fnode), pointer :: child2

    @:ASSERT(associated(node))
  #:block DEBUG_CODE
    if (present(default)) then
      @:ASSERT(all(shape(default) == shape(variableValue)))
    end if
  #:endblock DEBUG_CODE

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
    elseif (nReadItem /= size(variableValue)) then
      call detailedError(node, MSG_MISSING_VALUES)
    end if
    if (present(modifier)) then
      modifier = modif
    elseif (len(modif) > 0) then
      call detailedError(child2, MSG_NOMODIFIER)
    end if
    variableValue(:,:) = reshape(buffer, shape(variableValue))
    if (present(child)) then
      child => child2
    end if

  end subroutine getChVal_intR2


  !> Returns the value (the child) of a child node as a linked list of strings.
  !>
  !> In order to prevent a double packaging (from array to linked list and then from linked list to
  !> array), the setting of defaults for list types is not allowed. The presence of the child must
  !> be explicitely queried in the caller routine and an eventual default setting must be set with
  !> an explicit setChildValue call.
  subroutine getChVal_lString(node, name, variableValue, modifier, child)

    !> The node to investigate.
    type(fnode), pointer :: node

    !> Name of the child to look for
    character(len=*), intent(in) :: name

    !> Value on return
    type(TListString), intent(inout) :: variableValue

    !> Modifier of the child on return
    type(string), intent(inout), optional :: modifier

    !> Pointer to the child node (with the spec. name) on return
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
      call getChVal_lString_h(char(text), variableValue, child2)
      call setAttribute(child2, attrProcessed, "")
    else
      call detailedError(node, MSG_MISSING_FIELD // name)
    end if
    if (present(child)) then
      child => child2
    end if

  end subroutine getChVal_lString


  !> Helper function for getChVal_lString to avoid string to character conversion in the do-loop.
  subroutine getChVal_lString_h(text, variableValue, node)

    !> Text to parse
    character(len=*), intent(in) :: text

    !> Contains the value of the parsed text
    type(TListString), intent(inout) :: variableValue

    !> node for error handling
    type(fnode), pointer :: node

    integer :: iStart, iErr
    type(string) :: token

    iStart = 1
    call getNextToken(text, token, iStart, iErr)
    do while (iErr == TOKEN_OK)
      call append(variableValue, trim(unquote(char(token))))
      call getNextToken(text, token, iStart, iErr)
    end do
    if (iErr == TOKEN_ERROR) then
      call detailedError(node, "Invalid string")
    end if

  end subroutine getChVal_lString_h


  !> Returns the value (the child) of a child node as a linked list of reals.
  !>
  !> In order to prevent a double packaging (from array to linked list and then from linked list to
  !> array), the setting of defaults for list types is not allowed. The presence of the child must
  !> be explicitely queried in the caller routine and an eventual default setting must be set with
  !> an explicit setChildValue call.
  subroutine getChVal_lReal(node, name, variableValue, modifier, child)

    !> The node to investigate.
    type(fnode), pointer :: node

    !> Name of the child to look for
    character(len=*), intent(in) :: name

    !> Value on return
    type(TListReal), intent(inout) :: variableValue

    !> Modifier of the child on return
    type(string), intent(inout), optional :: modifier

    !> Pointer to the child node (with the spec. name) on return
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
      call getChVal_lReal_h(char(text), variableValue, child2)
      call setAttribute(child2, attrProcessed, "")
    else
      call detailedError(node, MSG_MISSING_FIELD // name)
    end if
    if (present(child)) then
      child => child2
    end if

  end subroutine getChVal_lReal


  !> Helper function for getChVal_lReal to avoid string to character conversion in the do-loop.
  subroutine getChVal_lReal_h(text, variableValue, node)

    !> text  Text to parse
    character(len=*), intent(in) :: text

    !> value Contains the value of the parsed text
    type(TListReal), intent(inout) :: variableValue
    type(fnode), pointer :: node

    integer :: iStart, iErr
    real(dp) :: buffer

    iStart = 1
    call getNextToken(text, buffer, iStart, iErr)
    do while (iErr == TOKEN_OK)
      call append(variableValue, buffer)
      call getNextToken(text, buffer, iStart, iErr)
    end do
    if (iErr == TOKEN_ERROR) then
      call detailedError(node, "Invalid real value")
    end if

  end subroutine getChVal_lReal_h


  !> Returns the value (the child) of a child node as a linked list of rank one real arrays.
  !>
  !> In order to prevent a double packaging (from array to linked list and then from linked list to
  !> array), the setting of defaults for list types is not allowed. The presence of the child must
  !> be explicitely queried in the caller routine and an eventual default setting must be set with
  !> an explicit setChildValue call.
  subroutine getChVal_lRealR1(node, name, dim, variableValue, modifier, child)

    !> The node to investigate.
    type(fnode), pointer :: node

    !> Name of the child to look for
    character(len=*), intent(in) :: name

    !> Dimension of the arrays
    integer, intent(in) :: dim

    !> Value on return
    type(TListRealR1), intent(inout) :: variableValue

    !> Modifier of the child on return
    type(string), intent(inout), optional :: modifier

    !> Pointer to the child node (with the spec. name) on return
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
      call getChVal_lRealR1_h(char(text), dim, variableValue, child2)
      call setAttribute(child2, attrProcessed, "")
    else
      call detailedError(node, MSG_MISSING_FIELD // name)
    end if
    if (present(child)) then
      child => child2
    end if

  end subroutine getChVal_lRealR1


  !> Helper function for getChVal_lReal to avoid string to character conversion in the do-loop.
  subroutine getChVal_lRealR1_h(text, dim, variableValue, node)

    !> Text to parse
    character(len=*), intent(in) :: text

    !> buffer sizing
    integer, intent(in) :: dim

    !> Contains the value of the parsed text
    type(TListRealR1), intent(inout) :: variableValue

    !> nodes for error handling
    type(fnode), pointer :: node

    integer :: iStart, iErr
    real(dp) :: buffer(dim)
    integer :: nItem

    iStart = 1
    call getNextToken(text, buffer, iStart, iErr, nItem)
    do while (iErr == TOKEN_OK)
      call append(variableValue, buffer)
      call getNextToken(text, buffer, iStart, iErr, nItem)
    end do
    if (iErr == TOKEN_ERROR) then
      call detailedError(node, "Invalid real value")
    elseif (iErr == TOKEN_EOS .and. nItem /= 0) then
      call detailedError(node, "Unexpected end of data")
    end if

  end subroutine getChVal_lRealR1_h


  !> Returns the value (the child) of a child node as linked list of integers.
  !>
  !> In order to prevent a double packaging (from array to linked list and then from linked list to
  !> array), the setting of defaults for list types is not allowed. The presence of the child must
  !> be explicitely queried in the caller routine and an eventual default setting must be set with
  !> an explicit setChildValue call.
  subroutine getChVal_lInt(node, name, variableValue, modifier, child)

    !> The node to investigate.
    type(fnode), pointer :: node

    !> Name of the child to look for
    character(len=*), intent(in) :: name

    !> Value on return
    type(TListInt), intent(inout) :: variableValue

    !> Modifier of the child on return
    type(string), intent(inout), optional :: modifier

    !> Pointer to the child node (with the spec. name) on return
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
      call getChVal_lInt_h(char(text), variableValue, child2)
      call setAttribute(child2, attrProcessed, "")
    else
      call detailedError(node, MSG_MISSING_FIELD // name)
    end if
    if (present(child)) then
      child => child2
    end if

  end subroutine getChVal_lInt


  !> Helper function for getChVal_lReal to avoid string to character conversion in the do-loop.
  subroutine getChVal_lInt_h(text, variableValue, node)

    !> Text to parse
    character(len=*), intent(in) :: text

    !> Contains the value of the parsed text
    type(TListInt), intent(inout) :: variableValue

    !> node for error handling
    type(fnode), pointer :: node

    integer :: iStart, iErr
    integer :: buffer

    iStart = 1
    call getNextToken(text, buffer, iStart, iErr)
    do while (iErr == TOKEN_OK)
      call append(variableValue, buffer)
      call getNextToken(text, buffer, iStart, iErr)
    end do
    if (iErr == TOKEN_ERROR) then
      call detailedError(node, "Invalid real value")
    end if

  end subroutine getChVal_lInt_h


  !> Returns the value (the child) of a child node as linked list of rank one integer arrays.
  !>
  !> In order to prevent a double packaging (from array to linked list and then from linked list to
  !> array), the setting of defaults for list types is not allowed. The presence of the child must
  !> be explicitely queried in the caller routine and an eventual default setting must be set with
  !> an explicit setChildValue call.
  subroutine getChVal_lIntR1(node, name, dim, variableValue, modifier, child)

    !> The node to investigate.
    type(fnode), pointer :: node

    !> Name of the child to look for
    character(len=*), intent(in) :: name

    !> Value on return
    integer, intent(in) :: dim

    !> Modifier of the child on return
    type(TListIntR1), intent(inout) :: variableValue

    !> Pointer to the child node (with the spec. name) on return
    type(string), intent(inout), optional :: modifier

    !> the child itself
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
      call getChVal_lIntR1_h(char(text), dim, variableValue, child2)
      call setAttribute(child2, attrProcessed, "")
    else
      call detailedError(node, MSG_MISSING_FIELD // name)
    end if
    if (present(child)) then
      child => child2
    end if

  end subroutine getChVal_lIntR1


  !> Helper function for getChVal_lReal to avoid string to character conversion in the do-loop.
  subroutine getChVal_lIntR1_h(text, dim, variableValue, node)

    !> Text to parse
    character(len=*), intent(in) :: text

    !> buffer sizing
    integer, intent(in) :: dim

    !> Contains the value of the parsed text
    type(TListIntR1), intent(inout) :: variableValue

    !> node for error handling
    type(fnode), pointer :: node

    integer :: iStart, iErr
    integer :: buffer(dim)
    integer :: nItem

    iStart = 1
    call getNextToken(text, buffer, iStart, iErr, nItem)
    do while (iErr == TOKEN_OK)
      call append(variableValue, buffer)
      call getNextToken(text, buffer, iStart, iErr, nItem)
    end do
    if (iErr == TOKEN_ERROR) then
      call detailedError(node, "Invalid real value")
    elseif (iErr == TOKEN_EOS .and. nItem /= 0) then
      call detailedError(node, "Unexpected end of data")
    end if

  end subroutine getChVal_lIntR1_h


  !> Returns the value (the child) of a child node as a linked list rank one integer and rank one
  !> real arrays.
  !>
  !> In order to prevent a double packaging (from array to linked list and then from linked list to
  !> array), the setting of defaults for list types is not allowed. The presence of the child must
  !> be explicitely queried in the caller routine and an eventual default setting must be set with
  !> an explicit setChildValue call.
  subroutine getChVal_lIntR1RealR1(node, name, dimInt, valueInt, dimReal, valueReal, modifier, &
      & child)

    !> The node to investigate.
    type(fnode), pointer :: node

    !> Name of the child to look for
    character(len=*), intent(in) :: name

    !> Dimension of the integer arrays in the list
    integer, intent(in) :: dimInt

    !> List of integer arrays on return
    integer, intent(in) :: dimReal

    !> Dimensio of the real arrays in the list
    type(TListIntR1), intent(inout) :: valueInt

    !> List of real array on return
    type(TListRealR1), intent(inout) :: valueReal

    !> Modifier of the child on return
    type(string), intent(inout), optional :: modifier

    !> Pointer to the child on return
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


  !> Helper function for getChVal_lIntR1RealR1 to avoid string to char conversion in the do-loop.
  subroutine getChVal_lIntR1RealR1_h(text, dimInt, valueInt, dimReal, valueReal, node)

    !> Text to parse
    character(len=*), intent(in) :: text

    !> integer buffer dimensioning
    integer, intent(in) :: dimInt

    !> Contains the value of the integer in the parsed text
    type(TListIntR1), intent(inout) :: valueInt

    !> real buffer dimensioning
    integer, intent(in) :: dimReal

    !> Contains the value of the real in the parsed text
    type(TListRealR1), intent(inout) :: valueReal

    !> for error handling
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


  !> Returns the value (the child) of a child node as a linked list of string, rank one integer and
  !> rank one real arrays.
  !>
  !> In order to prevent a double packaging (from array to linked list and then from linked list to
  !> array), the setting of defaults for list types is not allowed. The presence of the child must
  !> be explicitely queried in the caller routine and an eventual default setting must be set with
  !> an explicit setChildValue call.
  subroutine getChVal_lStringIntR1RealR1(node, name, valueStr, dimInt, valueInt, dimReal, &
      & valueReal, modifier, child)

    !> The node to investigate.
    type(fnode), pointer :: node

    !> Name of the child to look for
    character(len=*), intent(in) :: name

    !> List of strings on return.
    type(TListString), intent(inout) :: valueStr

    !> Dimension of the integer arrays in the list
    integer, intent(in) :: dimInt

    !> List of integer arrays on return
    type(TListIntR1), intent(inout) :: valueInt

    !> Dimension of the real arrays in the list
    integer, intent(in) :: dimReal

    !> List of real array on return
    type(TListRealR1), intent(inout) :: valueReal

    !> Modifier of the child on return
    type(string), intent(inout), optional :: modifier

    !> Pointer to the child on return
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


  !> Helper function for getChVal_lIntR1RealR1 to avoid string to char conversion in the do-loop.
  subroutine getChVal_lStringIntR1RealR1_h(text, valueStr, dimInt, valueInt, dimReal, valueReal, &
      & node)

    !> Text to parse
    character(len=*), intent(in) :: text

    !> Contains the string part of the parsed text
    type(TListString), intent(inout) :: valueStr

    !> integer buffer dimensioning
    integer, intent(in) :: dimInt

    !> Contains the integer part of the parsed text
    type(TListIntR1), intent(inout) :: valueInt

    !> integer buffer dimensioning
    integer, intent(in) :: dimReal

    !> Contains the real value part of the parsed text
    type(TListRealR1), intent(inout) :: valueReal

    !> for error handling
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


  !> Returns the value (the child) of a child node as a node.
  !>
  !> Caveat: If allowEmptyValue is set to .true. and the child has no subnodes (empty value) then
  !> the returned value is an unassociated pointer
  subroutine getChVal_node(node, name, variableValue, default, modifier, child, list, &
      & allowEmptyValue, dummyValue)

    !> The node to investigate.
    type(fnode), pointer :: node

    !> Name of the child to look for
    character(len=*), intent(in) :: name

    !> Value on return
    type(fnode), pointer :: variableValue

    !> Default value for the child, if child is not found. If the empty string is passed as default
    !> value, the child is created but no value is added to it. The returned value pointer will be

    !> unassociated. (allowEmptyValue must be explicitely set to .true.)
    character(len=*), intent(in), optional :: default

    !> Modifier of the child on return
    type(string), intent(inout), optional :: modifier

    !> Pointer to the child node (with the spec. name) on return
    type(fnode), pointer, optional :: child

    !> If the node created as default should be tagged as list.
    logical, intent(in), optional :: list

    !> If the child is allowed to have an empty value.
    logical, intent(in), optional :: allowEmptyValue

    !> If true, the value is not marked as processed.
    logical, intent(in), optional :: dummyValue

    type(string) :: modif
    type(fnode), pointer :: child2
    logical :: tList, tAllowEmptyVal, tDummyValue

    @:ASSERT(associated(node))
  #:block DEBUG_CODE
    if (present(default)) then
      if (len(default) == 0) then
        @:ASSERT(present(allowEmptyValue))
        @:ASSERT(allowEmptyValue)
      end if
    end if
  #:endblock DEBUG_CODE

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
      variableValue => getFirstChild(child2)
      if ((.not. associated(variableValue)) .and. (.not. tAllowEmptyVal)) then
        call detailedError(child2, "Missing value")
      end if
      call setAttribute(child2, attrProcessed, "")
    elseif (present(default)) then
      if (present(modifier)) then
        modifier = ""
      end if
      if (len(default) > 0) then
        variableValue => createElement(tolower(default))
        call setChildValue(node, name, variableValue, .false., child=child2, list=tList)
        call setAttribute(variableValue, attrName, default)
      else
        nullify(variableValue)
        call setChild(node, name, child2, .false., list=tList)
      end if
    else
      call detailedError(node, MSG_MISSING_FIELD // name)
    end if
    if (associated(variableValue) .and. .not. tDummyValue) then
      if (getNodeType(variableValue) == ELEMENT_NODE) then
        call setAttribute(variableValue, attrProcessed, "")
      end if
    end if
    if (present(child)) then
      child => child2
    end if

  end subroutine getChVal_node


  !> Converts a string containing atom indices, ranges and species names to a list of atom indices.
  subroutine convAtomRangeToInt(str, speciesNames, species, node, val, ishift, maxRange)

    !> String to convert
    character(len=*), intent(in) :: str

    !> Contains the valid species names.
    character(len=*), intent(in) :: speciesNames(:)

    !> Contains for every atom its species index
    integer, intent(in) :: species(:)

    !> Top node for detailed errors.
    type(fnode), pointer :: node

    !> Integer list of atom indices on return.
    integer, allocatable, intent(out) :: val(:)

    !> Shift to be applied to provided atomic indices
    integer, intent(in), optional :: ishift

    !> Upper range of atoms
    integer, intent(in), optional :: maxRange

    type(string) :: buffer
    type(TListInt) :: li
    integer :: nAtom, iStart, iostat, shift

    shift = 0
    if (present(ishift)) then
      shift = ishift
    end if
    if (present(maxRange)) then
      nAtom = maxRange
    else
      nAtom = size(species)
    end if
    call init(li)
    iStart = 1
    call getNextToken(str, buffer, iStart, iostat)
    do while (iostat == TOKEN_OK)
      call convAtomRangeToIntProcess(char(buffer), speciesNames, species, nAtom, node, li, shift)
      call getNextToken(str, buffer, iStart, iostat)
    end do
    allocate(val(len(li)))
    if (len(li) > 0) then
      call asArray(li, val)
    end if
    call destruct(li)

  end subroutine convAtomRangeToInt


  !> Helper routine.
  subroutine convAtomRangeToIntProcess(cbuffer, speciesNames, species, nAtom, node, li, shift)

    !> Chunk of the specified atoms
    character(len=*), intent(in) :: cbuffer

    !> Name of chemical species
    character(len=*), intent(in) :: speciesNames(:)

    !> Chemical species of atoms
    integer, intent(in) :: species(:)

    !> Upper limit on range of atoms
    integer, intent(in) :: nAtom

    !> Top node for detailed errors.
    type(fnode), pointer :: node

    !> List of the converted atom numbers
    type(TListInt), intent(inout) :: li

    !> Shift in lower range of index
    integer, intent(in) :: shift 

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
        bounds = bounds + shift
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
        ii = ii + shift
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
      do ii = 1, nAtom
        if (species(ii) == iPos) then
          call append(li, ii)
        end if
      end do
    end if

  end subroutine convAtomRangeToIntProcess


  !> Converts a string containing indices and ranges to a list of indices.
  subroutine convRangeToInt(str, node, val, nMax)

    !> String to convert
    character(len=*), intent(in) :: str

    !> Top node for detailed errors.
    type(fnode), pointer :: node

    !> Integer list of atom indices on return.
    integer, allocatable, intent(out) :: val(:)

    !> Maximum number for an index
    integer, intent(in) :: nMax

    type(string) :: buffer
    type(TListInt) :: li
    integer :: iStart, iostat

    call init(li)
    iStart = 1
    call getNextToken(str, buffer, iStart, iostat)
    do while (iostat == TOKEN_OK)
      call convRangeToIntProcess(char(buffer), nMax, node, li)
      call getNextToken(str, buffer, iStart, iostat)
    end do
    allocate(val(len(li)))
    if (len(li) > 0) then
      call asArray(li, val)
    end if
    call destruct(li)

  end subroutine convRangeToInt

  !> Helper routine.
  subroutine convRangeToIntProcess(cbuffer, nMax, node, li)
    character(len=*), intent(in) :: cbuffer
    integer, intent(in) :: nMax
    type(fnode), pointer :: node
    type(TListInt), intent(inout) :: li

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

  end subroutine convRangeToIntProcess

  !> Returns a child node with a specified name
  subroutine getChild(node, name, child, requested, modifier)

    !> Node to investigate
    type(fnode), pointer :: node

    !> Name of the child node to look for
    character(len=*), intent(in) :: name

    !> Contains a pointer to the child on return
    type(fnode), pointer :: child

    !> If true and child not found, error is issued
    logical, intent(in), optional :: requested

    !> Contains modifier on exit.
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


  !> Returns a list of children with the specified name.
  subroutine getChildren(node, name, children)

    !> Parent node to investigate
    type(fnode), pointer :: node

    !> Name of the children to look for
    character(len=*), intent(in) :: name

    !> List of the children.
    type(fnodeList), pointer :: children

    type(fnode), pointer :: child
    integer :: ii

    children => getChildrenByName(node, tolower(name))
    do ii = 1, getLength(children)
      call getItem1(children, ii, child)
      call setAttribute(child, attrProcessed, "")
    end do

  end subroutine getChildren


  !> Sets the value (child) of a child with given name.
  subroutine setChVal_logical(node, name, variableValue, replace, child, modifier)

    !> The node to investigate
    type(fnode), pointer :: node

    !> Name of the child to look for
    character(len=*), intent(in) :: name

    !> Value to set
    logical, intent(in) :: variableValue

    !> Replace if child with same name already exists
    logical, intent(in), optional :: replace

    !> Pointer to the child node (with the provided name)
    type(fnode), pointer, optional :: child

    !> Optional modifier for the child
    character(len=*), optional, intent(in) :: modifier

    type(string) :: strBuffer
    type(fnode), pointer :: child2
    logical :: tReplace

    if (present(replace)) then
      tReplace = replace
    else
      tReplace = .false.
    end if

    call getAsString(variableValue, strBuffer)
    call createChild_local(node, name, .false., tReplace, child2, &
        &variableValue=char(strBuffer))
    if (present(child)) then
      child => child2
    end if
    if (present(modifier)) then
      call setAttribute(child2, attrModifier, modifier)
    end if

  end subroutine setChVal_logical


  !> Sets the value (child) of a child with given name.
  subroutine setChVal_logicalR1(node, name, variableValue, replace, child, modifier)

    !> The node to investigate
    type(fnode), pointer :: node

    !> Name of the child to look for
    character(len=*), intent(in) :: name

    !> Value to set
    logical, intent(in) :: variableValue(:)

    !> Replace if child with same name already exists
    logical, intent(in), optional :: replace

    !> Pointer to the child node (with the provided name)
    type(fnode), pointer, optional :: child

    !> Optional modifier for the child
    character(len=*), optional, intent(in) :: modifier

    type(string) :: strBuffer
    type(fnode), pointer :: child2
    logical :: tReplace

    if (present(replace)) then
      tReplace = replace
    else
      tReplace = .false.
    end if

    call getAsString(variableValue, strBuffer)
    call createChild_local(node, name, .false., tReplace, child2, &
        &variableValue=char(strBuffer))
    if (present(child)) then
      child => child2
    end if
    if (present(modifier)) then
      call setAttribute(child2, attrModifier, modifier)
    end if

  end subroutine setChVal_logicalR1


  !> Writes the text representation of a node and its value to an xmlwriter.
  subroutine writeChVal_logical(xf, name, variableValue)

    !> Xmlwriter stream
    type(xmlf_t), intent(inout) :: xf

    !> Name of the node
    character(len=*), intent(in) :: name

    !> Value of the node
    logical, intent(in) :: variableValue

    type(string) :: strBuffer

    call getAsString(variableValue, strBuffer)
    call writeChild_local(xf, name, char(strBuffer))

  end subroutine writeChVal_logical


  !> Writes the text representation of a node and its value to an xmlwriter.
  subroutine writeChVal_logicalR1(xf, name, variableValue)

    !> Xmlwriter stream
    type(xmlf_t), intent(inout) :: xf

    !> Name of the node
    character(len=*), intent(in) :: name

    !> Value of the node
    logical, intent(in) :: variableValue(:)

    type(string) :: strBuffer

    call getAsString(variableValue, strBuffer)
    call writeChild_local(xf, name, char(strBuffer))

  end subroutine writeChVal_logicalR1


  !> Returns the text representation of the passed object
  subroutine getAsString_logical(variableValue, strBuffer)

    !> Value to represent
    logical, intent(in) :: variableValue

    !> Text representation on exit
    type(string), intent(inout) :: strBuffer

    if (variableValue) then
      strBuffer = LOGICAL_TRUE
    else
      strBuffer = LOGICAL_FALSE
    end if

  end subroutine getAsString_logical


  !> Returns the text representation of the passed object
  subroutine getAsString_logicalR1(variableValue, strBuffer)

    !> Value to represent
    logical, intent(in) :: variableValue(:)

    !> Text representation on exit
    type(string), intent(inout) :: strBuffer

    character(len=nCharLogical) :: buffer
    integer :: buffLen, len
    integer :: ii

    call resize_string(strBuffer, preAllocSize)
    len = 0
    do ii = 1, size(variableValue)
      if (variableValue(ii)) then
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


  !> Sets the value (child) of a child with given name.
  !>
  !> Caveat: This subroutines assumes, that a real can be represented as text with less than
  !> nCharReal characters.
  subroutine setChVal_real(node, name, variableValue, replace, child, modifier)

    !> The node to investigate
    type(fnode), pointer :: node

    !> Name of the child to look for
    character(len=*), intent(in) :: name

    !> Value to set
    real(dp), intent(in) :: variableValue

    !> Replace if child with same name already exists
    logical, intent(in), optional :: replace

    !> Pointer to the child node (with the provided name)
    type(fnode), pointer, optional :: child

    !> Optional modifier for the child
    character(len=*), optional, intent(in) :: modifier

    type(string) :: strBuffer
    type(fnode), pointer :: child2
    logical :: tReplace

    if (present(replace)) then
      tReplace = replace
    else
      tReplace = .false.
    end if

    call getAsString(variableValue, strBuffer)
    call createChild_local(node, name, .false., tReplace, child2, &
        &variableValue=char(strBuffer))
    if (present(child)) then
      child => child2
    end if
    if (present(modifier)) then
      call setAttribute(child2, attrModifier, modifier)
    end if

  end subroutine setChVal_real


  !> Writes the text representation of a node and its value to an xmlwriter.
  subroutine writeChVal_real(xf, name, variableValue)

    !> Xmlwriter stream
    type(xmlf_t), intent(inout) :: xf

    !> Name of the node
    character(len=*), intent(in) :: name

    !> Value of the node
    real(dp), intent(in) :: variableValue

    type(string) :: strBuffer

    call getAsString(variableValue, strBuffer)
    call writeChild_local(xf, name, char(strBuffer))

  end subroutine writeChVal_real


  !> Returns the text representation of the passed object
  subroutine getAsString_real(variableValue, strBuffer)

    !> Value to represent
    real(dp), intent(in) :: variableValue

    !> Text representation on exit
    type(string), intent(inout) :: strBuffer

    character(len=nCharReal) :: buffer

    write (buffer, *) variableValue
    strBuffer = trim(adjustl(buffer))

  end subroutine getAsString_real


  !> Sets the value (child) of a child with given name.
  !>
  !> Caveat: This subroutines assumes, that a real can be represented as text with less than
  !> nCharReal characters.
  subroutine setChVal_realR1(node, name, variableValue, replace, child, modifier)

    !> The node to investigate
    type(fnode), pointer :: node

    !> Name of the child to look for
    character(len=*), intent(in) :: name

    !> Value to set
    real(dp), intent(in) :: variableValue(:)

    !> Replace if child with same name already exists
    logical, intent(in), optional :: replace

    !> Pointer to the child node (with the provided name)
    type(fnode), pointer, optional :: child

    !> Optional modifier for the child
    character(len=*), optional, intent(in) :: modifier

    type(string) :: strBuffer
    type(fnode), pointer :: child2
    logical :: tReplace

    if (present(replace)) then
      tReplace = replace
    else
      tReplace = .false.
    end if
    call getAsString(variableValue, strBuffer)
    call createChild_local(node, name, .true., tReplace, child2, &
        &variableValue=char(strBuffer))
    if (present(child)) then
      child => child2
    end if
    if (present(modifier)) then
      call setAttribute(child2, attrModifier, modifier)
    end if

  end subroutine setChVal_realR1


  !> Writes the text representation of a node and its value to an xmlwriter.
  subroutine writeChVal_realR1(xf, name, variableValue)

    !> Xmlwriter stream
    type(xmlf_t), intent(inout) :: xf

    !> Name of the node
    character(len=*), intent(in) :: name

    !> Value of the node
    real(dp), intent(in) :: variableValue(:)

    type(string) :: strBuffer

    call getAsString(variableValue, strBuffer)
    call writeChild_local(xf, name, char(strBuffer))

  end subroutine writeChVal_realR1


  !> Returns the text representation of the passed object
  subroutine getAsString_realR1(variableValue, strBuffer)

    !> Value to represent
    real(dp), intent(in) :: variableValue(:)

    !> Text representation on exit
    type(string), intent(inout) :: strBuffer

    character(len=nCharReal) :: buffer
    integer :: buffLen, len
    integer :: ii

    call resize_string(strBuffer, preAllocSize)
    len = 0
    do ii = 1, size(variableValue)
      write (buffer, *) variableValue(ii)
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


  !> Sets the value (child) of a child with given name.  The node to investigate
  !>
  !> This is just a wrapper around the rank one version, to make sure that two dimensional arrays
  !> are pretty printed. For higher ranked arrays the rank one version should be used with some
  !> reshaping before.
  !>
  !> This subroutines assumes, that a real can be represented as text with less than nCharReal
  !> characters.
  subroutine setChVal_realR2(node, name, variableValue, replace, child, modifier)

    !> node to process from
    type(fnode), pointer :: node

    !> Name of the child to look for
    character(len=*), intent(in) :: name

    !> Value to set
    real(dp), intent(in) :: variableValue(:,:)

    !> Replace if child with same name already exists
    logical, intent(in), optional :: replace

    !> Pointer to the child node (with the provided name)
    type(fnode), pointer, optional :: child

    !> Optional modifier for the child
    character(len=*), intent(in), optional :: modifier

    type(fnode), pointer :: child2
    type(string) :: strBuffer
    logical :: tReplace

    if (present(replace)) then
      tReplace = replace
    else
      tReplace = .false.
    end if

    call getAsString(variableValue, strBuffer)
    call createChild_local(node, name, .true., tReplace, child2, &
        &variableValue=char(strBuffer))
    if (present(child)) then
      child => child2
    end if
    if (present(modifier)) then
      call setAttribute(child2, attrModifier, modifier)
    end if

  end subroutine setChVal_realR2


  !> Writes the text representation of a node and its value to an xmlwriter.
  subroutine writeChVal_realR2(xf, name, variableValue)

    !> Xmlwriter stream
    type(xmlf_t), intent(inout) :: xf

    !> Name of the node
    character(len=*), intent(in) :: name

    !> Value of the node
    real(dp), intent(in) :: variableValue(:,:)

    type(string) :: strBuffer

    call getAsString(variableValue, strBuffer)
    call writeChild_local(xf, name, char(strBuffer))

  end subroutine writeChVal_realR2


  !> Returns the text representation of the passed object
  subroutine getAsString_realR2(variableValue, strBuffer)

    !> Value to represent
    real(dp), intent(in) :: variableValue(:,:)

    !> Text representation on exit
    type(string), intent(inout) :: strBuffer

    character(len=nCharReal) :: buffer
    integer :: ii, jj

    call resize_string(strBuffer, preAllocSize)
    do ii = 1, size(variableValue, dim=2)
      do jj = 1, size(variableValue, dim=1)
        write (buffer, *) variableValue(jj, ii)
        buffer = adjustl(buffer)
        call append_to_string(strBuffer, space // trim(buffer))
      end do
      call append_to_string(strBuffer, newline)
    end do

  end subroutine getAsString_realR2


  !> Sets the value (child) of a child with given name.
  !>
  !> Caveat: This subroutines assumes, that an integer can be represented as text with less than
  !> nCharInt characters.
  subroutine setChVal_int(node, name, variableValue, replace, child, modifier)

    !> The node to investigate
    type(fnode), pointer :: node

    !> Name of the child to look for
    character(len=*), intent(in) :: name

    !> Value to set
    integer, intent(in) :: variableValue

    !> Replace if child with same name already exists
    logical, intent(in), optional :: replace

    !> Pointer to the child node (with the provided name)
    type(fnode), pointer, optional :: child

    !> Optional modifier for the child
    character(len=*), optional, intent(in) :: modifier

    type(fnode), pointer :: child2
    type(string) :: strBuffer
    logical :: tReplace

    if (present(replace)) then
      tReplace = replace
    else
      tReplace = .false.
    end if
    call getAsString(variableValue, strBuffer)
    call createChild_local(node, name, .false., tReplace, child2, &
        &variableValue=char(strBuffer))
    if (present(child)) then
      child => child2
    end if
    if (present(modifier)) then
      call setAttribute(child2, attrModifier, modifier)
    end if

  end subroutine setChVal_int


  !> Writes the text representation of a node and its value to an xmlwriter.
  subroutine writeChVal_int(xf, name, variableValue)

    !> Xmlwriter stream
    type(xmlf_t), intent(inout) :: xf

    !> Name of the node
    character(len=*), intent(in) :: name

    !> Value of the node
    integer, intent(in) :: variableValue

    type(string) :: strBuffer

    call getAsString(variableValue, strBuffer)
    call writeChild_local(xf, name, char(strBuffer))

  end subroutine writeChVal_int


  !> Returns the text representation of the passed object
  subroutine getAsString_int(variableValue, strBuffer)

    !> Value to represent
    integer, intent(in) :: variableValue

    !> Text representation on exit
    type(string), intent(inout) :: strBuffer

    character(len=nCharInt) :: buffer

    write (buffer, *) variableValue
    strBuffer = trim(adjustl(buffer))

  end subroutine getAsString_int


  !> Sets the value (child) of a child with given name.
  !>
  !> Caveat: This subroutines assumes, that an integer can be represented as text with less than
  !> nCharInt characters.
  subroutine setChVal_intR1(node, name, variableValue, replace, child, modifier)

    !> The node to investigate
    type(fnode), pointer :: node

    !> Name of the child to look for
    character(len=*), intent(in) :: name

    !> Value to set
    integer, intent(in) :: variableValue(:)

    !> Replace if child with same name already exists
    logical, intent(in), optional :: replace

    !> Optional modifier for the child
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
    call getAsString(variableValue, strBuffer)
    call createChild_local(node, name, .true., tReplace, child2, &
        &variableValue=char(strBuffer))
    if (present(child)) then
      child => child2
    end if
    if (present(modifier)) then
      call setAttribute(child2, attrModifier, modifier)
    end if

  end subroutine setChVal_intR1


  !> Writes the text representation of a node and its value to an xmlwriter.
  subroutine writeChVal_intR1(xf, name, variableValue)

    !> Xmlwriter stream
    type(xmlf_t), intent(inout) :: xf

    !> Name of the node
    character(len=*), intent(in) :: name

    !> Value of the node
    integer, intent(in) :: variableValue(:)

    type(string) :: strBuffer

    call getAsString(variableValue, strBuffer)
    call writeChild_local(xf, name, char(strBuffer))

  end subroutine writeChVal_intR1


  !> Returns the text representation of the passed object
  subroutine getAsString_intR1(variableValue, strBuffer)

    !> Value to represent
    integer, intent(in) :: variableValue(:)

    !> Text representation on exit
    type(string), intent(inout) :: strBuffer

    character(len=nCharInt) :: buffer
    integer :: buffLen, len
    integer :: ii

    call resize_string(strBuffer, preAllocSize)
    len = 0
    do ii = 1, size(variableValue)
      write (buffer, *) variableValue(ii)
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


  !> Sets the value (child) of a child with given name.
  !>
  !> This is just a wrapper around the rank one version, to make sure that two dimensional arrays
  !> are pretty printed. For higher ranked arrays the rank one version should be used with some
  !> reshaping beforehand.
  !>
  !> Caveat: This subroutines assumes, that an integer can be represented as text with less than
  !> nCharInt characters.
  subroutine setChVal_intR2(node, name, variableValue, replace, child, modifier)

    !> The node to investigate
    type(fnode), pointer :: node

    !> Name of the child to look for
    character(len=*), intent(in) :: name

    !> Value to set
    integer, intent(in) :: variableValue(:,:)

    !> Replace if child with same name already exists
    logical, intent(in), optional :: replace

    !> Pointer to the child node (with the provided name)
    type(fnode), pointer, optional :: child

    !> Optional modifier for the child
    character(len=*), optional, intent(in) :: modifier

    type(fnode), pointer :: child2
    type(string) :: strBuffer
    logical :: tReplace

    if (present(replace)) then
      tReplace = replace
    else
      tReplace = .false.
    end if
    call getAsString(variableValue, strBuffer)
    call createChild_local(node, name, .true., tReplace, child2, &
        &variableValue=char(strBuffer))
    if (present(child)) then
      child => child2
    end if
    if (present(modifier)) then
      call setAttribute(child2, attrModifier, modifier)
    end if

  end subroutine setChVal_intR2


  !> Writes the text representation of a node and its value to an xmlwriter.
  subroutine writeChVal_intR2(xf, name, variableValue)

    !> Xmlwriter stream
    type(xmlf_t), intent(inout) :: xf

    !> Name of the node
    character(len=*), intent(in) :: name

    !> Value of the node
    integer, intent(in) :: variableValue(:,:)

    type(string) :: strBuffer

    call getAsString(variableValue, strBuffer)
    call writeChild_local(xf, name, char(strBuffer))

  end subroutine writeChVal_intR2


  !> Returns the text representation of the passed object
  subroutine getAsString_intR2(variableValue, strBuffer)

    !> Value to represent
    integer, intent(in) :: variableValue(:,:)

    !> Text representation on exit
    type(string), intent(inout) :: strBuffer

    character(len=nCharInt) :: buffer
    integer :: ii, jj

    call resize_string(strBuffer, preAllocSize)
    do ii = 1, size(variableValue, dim=2)
      do jj = 1, size(variableValue, dim=1)
        write (buffer, *) variableValue(jj, ii)
        buffer = adjustl(buffer)
        call append_to_string(strBuffer, space // trim(buffer))
      end do
      call append_to_string(strBuffer, newline)
    end do

  end subroutine getAsString_intR2


  !> Sets the value (child) of a child with given name.
  subroutine setChVal_char(node, name, variableValue, replace, child, omitQuotes, modifier)

    !> The node to investigate
    type(fnode), pointer :: node

    !> Name of the child to look for
    character(len=*), intent(in) :: name

    !> Value to set
    character(len=*), intent(in) :: variableValue

    !> Replace if child with same name already exists
    logical, intent(in), optional :: replace

    !> Pointer to the child node (with the provided name)
    type(fnode), pointer, optional :: child

    !> If quotes around the string should be omitted
    logical, intent(in), optional :: omitQuotes

    !> Optional modifier for the child
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
          &variableValue='"'//variableValue//'"')
    else
      call createChild_local(node, name, .false., tReplace, child2, variableValue=variableValue)
    end if

    if (present(child)) then
      child => child2
    end if
    if (present(modifier)) then
      call setAttribute(child2, attrModifier, modifier)
    end if

  end subroutine setChVal_char


  !> Sets the value (child) of a child with given name.
  subroutine setChVal_charR1(node, name, variableValue, replace, child, modifier)

    !> The node to investigate
    type(fnode), pointer :: node

    !> Name of the child to look for
    character(len=*), intent(in) :: name

    !> Value to set
    character(len=*), intent(in) :: variableValue(:)

    !> Replace if child with same name already exists
    logical, intent(in), optional :: replace

    !> Pointer to the child node (with the provided name)
    type(fnode), pointer, optional :: child

    !> Optional modifier for the child
    character(len=*), optional, intent(in) :: modifier

    type(string) :: strBuffer
    type(fnode), pointer :: child2
    logical :: tReplace

    if (present(replace)) then
      tReplace = replace
    else
      tReplace = .false.
    end if
    call getAsString(variableValue, strBuffer)
    call createChild_local(node, name, .true., tReplace, child2, &
        &variableValue=char(strBuffer))
    if (present(child)) then
      child => child2
    end if
    if (present(modifier)) then
      call setAttribute(child2, attrModifier, modifier)
    end if

  end subroutine setChVal_charR1


  !> Writes the text representation of a node and its value to an xmlwriter.
  subroutine writeChVal_charR1(xf, name, variableValue)

    !> Xmlwriter stream
    type(xmlf_t), intent(inout) :: xf

    !> Name of the node
    character(len=*), intent(in) :: name

    !> Value of the node
    character(len=*), intent(in) :: variableValue(:)

    type(string) :: strBuffer

    call getAsString(variableValue, strBuffer)
    call writeChild_local(xf, name, char(strBuffer))

  end subroutine writeChVal_charR1


  !> Returns the text representation of the passed object
  subroutine getAsString_charR1(variableValue, strBuffer)

    !> Value to represent
    character(len=*), intent(in) :: variableValue(:)

    !> Text representation on exit
    type(string), intent(inout) :: strBuffer

    integer :: buffLen, len
    integer :: ii

    call resize_string(strBuffer, preAllocSize)
    len = 0
    do ii = 1, size(variableValue)
      buffLen = len_trim(variableValue(ii))
      len = len + buffLen
      if (len > lineLength) then
        call append_to_string(strBuffer, newline // '"'//trim(variableValue(ii))//'"')
        len = buffLen
      else
        call append_to_string(strBuffer, space // '"'//trim(variableValue(ii))//'"')
      end if
    end do

  end subroutine getAsString_charR1


  !> Sets the value (child) of a child with given name.
  subroutine setChVal_intR2RealR2(node, name, intValue, realValue, replace, child, modifier)

    !> The node to investigate
    type(fnode), pointer :: node

    !> Name of the child to look for
    character(len=*), intent(in) :: name

    !> Value for the integers
    integer, intent(in) :: intValue(:,:)

    !> Value for the reals
    real(dp), intent(in) :: realValue(:,:)

    !> Replace if child with same name already exists
    logical, intent(in), optional :: replace

    !> Pointer to the child node (with the provided name)
    type(fnode), pointer, optional :: child

    !> Optional modifier for the child
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
        &variableValue=char(strBuffer))
    if (present(child)) then
      child => child2
    end if
    if (present(modifier)) then
      call setAttribute(child2, attrModifier, modifier)
    end if

  end subroutine setChVal_intR2RealR2


  !> Writes the text representation of a node and its value to an xmlwriter.
  subroutine writeChVal_intR2RealR2(xf, name, intValue, realValue)

    !> Xmlwriter stream
    type(xmlf_t), intent(inout) :: xf

    !> Name of the node
    character(len=*), intent(in) :: name

    !> Integer value of the node
    integer, intent(in) :: intValue(:,:)

    !> real values of the node
    real(dp), intent(in) :: realValue(:,:)

    type(string) :: strBuffer

    call getAsString(intValue, realValue, strBuffer)
    call writeChild_local(xf, name, char(strBuffer))

  end subroutine writeChVal_intR2RealR2


  !> Returns the text representation of the passed object
  subroutine getAsString_intR2RealR2(intValue, realValue, strBuffer)

    !> integer value in node
    integer, intent(in) :: intValue(:,:)

    !> real value in node
    real(dp), intent(in) :: realValue(:,:)

    !> Text representation on exit
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


  !> Sets the value (child) of a child with given name.
  subroutine setChVal_charR1IntR2RealR2(node, name, charValue, intValue, realValue, replace, &
      & child, modifier)

    !> The node to investigate
    type(fnode), pointer :: node

    !> Name of the child to look for
    character(len=*), intent(in) :: name

    !> Value for the characters
    character(len=*), intent(in) :: charValue(:)

    !> Value for the integers
    integer, intent(in) :: intValue(:,:)

    !> Value for the reals
    real(dp), intent(in) :: realValue(:,:)

    !> Replace if child with same name already exists
    logical, intent(in), optional :: replace

    !> Pointer to the child node (with the provided name)
    type(fnode), pointer, optional :: child

    !> Optional modifier for the child
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
        &variableValue=char(strBuffer))
    if (present(child)) then
      child => child2
    end if
    if (present(modifier)) then
      call setAttribute(child2, attrModifier, modifier)
    end if

  end subroutine setChVal_charR1IntR2RealR2


  !> Writes the text representation of a node and its value to an xmlwriter.
  subroutine writeChVal_charR1IntR2RealR2(xf, name, charValue, intValue, realValue)

    !> Xmlwriter stream
    type(xmlf_t), intent(inout) :: xf

    !> Name of the node
    character(len=*), intent(in) :: name

    !> character part of node
    character(len=*), intent(in) :: charValue(:)

    !> integer part of node
    integer, intent(in) :: intValue(:,:)

    !> real value part of node
    real(dp), intent(in) :: realValue(:,:)

    type(string) :: strBuffer

    call getAsString(charValue, intValue, realValue, strBuffer)
    call writeChild_local(xf, name, char(strBuffer))

  end subroutine writeChVal_charR1IntR2RealR2


  !> Returns the text representation of the passed object
  subroutine getAsString_charR1IntR2RealR2(charValue, intValue, realValue, strBuffer)

    !> character part of node
    character(len=*), intent(in) :: charValue(:)

    !> integer part of node
    integer, intent(in) :: intValue(:,:)

    !> real value part of node
    real(dp), intent(in) :: realValue(:,:)

    !> Text representation on exit
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


  !> Sets the value (child) of a child with given name.
  subroutine setChVal_node(node, name, variableValue, replace, child, modifier, list)

    !> The node to investigate
    type(fnode), pointer :: node

    !> Name of the child to look for
    character(len=*), intent(in) :: name

    !> Value to set
    type(fnode), pointer :: variableValue

    !> Replace if child with same name already exists
    logical, intent(in), optional :: replace

    !> Pointer to the child node (with the provided name)
    type(fnode), pointer, optional :: child

    !> Optional modifier for the child
    character(len=*), optional, intent(in) :: modifier

    !> If created child should be marked as a list.
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
    if (associated(variableValue)) then
      dummy => appendChild(child2, variableValue)
    end if
    if (present(child)) then
      child => child2
    end if
    if (present(modifier)) then
      call setAttribute(child2, attrModifier, modifier)
    end if

  end subroutine setChVal_node


  !> Workhorse for the setChildValue routines
  !>
  !> If an empty string is provided as child name, no child is created, and the current node is
  !> replace instead. The pointer "node" becames associated with the new node, since the old
  !> instance will be destroyed.
  subroutine createChild_local(node, name, list, replace, child, variableValue)

    !> The node to investigate
    type(fnode), pointer :: node

    !> Name of the child to create
    character(len=*), intent(in) :: name

    !> True, if child should be signed as a list
    logical, intent(in) :: list

    !> Replace if child with same name already exists
    logical, intent(in) :: replace

    !> Pointer to the created child on return
    type(fnode), pointer :: child

    !> Value to set (if empty, no child is appended to the created child)
    character(len=*), intent(in), optional :: variableValue

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
    if (present(variableValue)) then
      text => createTextNode(variableValue)
      dummy => appendChild(child, text)
    end if

  end subroutine createChild_local


  !> new child in the xml
  subroutine writeChild_local(xf, name, variableValue)

    !> xmlWriter stream
    type(xmlf_t), intent(inout) :: xf

    !> node name
    character(len=*), intent(in) :: name

    !> stored variale string
    character(len=*), intent(in) :: variableValue

    call xml_NewElement(xf, name)
    call xml_AddPCData(xf, variableValue)
    call xml_EndElement(xf, name)

  end subroutine writeChild_local


  !> Creates a child with the given name
  subroutine setChild(node, name, child, replace, list, modifier)

    !> Node to append the child to
    type(fnode), pointer :: node

    !> Name of the child node to append
    character(len=*), intent(in) :: name

    !> Contains the pointer to the added child node on return
    type(fnode), pointer :: child

    !> If an already existing child with the same name should be replaced
    logical, intent(in), optional :: replace

    !> If child should be signed as a list tag
    logical, intent(in), optional :: list

    !> Optional modifier for the child
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


  !> Returns the content of the first TEXT_NODE child of a given node or empty string, if such a
  !> node does not exist.
  !>
  !> Note: the document tree is normalized, every node has only one TEXT_NODE child.
  subroutine getFirstTextChild(node, str)

    !> The node to investigate.
    type(fnode), pointer :: node

    !> String representation of the TEXT_NODE.
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


  !> Checks if error flag signals an error. If yes, raises error.
  subroutine checkError(node, iErr, msg)

    !> Node which the error flag was set for
    type(fnode), pointer :: node

    !> Content of the error flag.
    integer, intent(in) :: iErr

    !> Message to print, if error occured
    character(len=*), intent(in) :: msg

    if (iErr == TOKEN_ERROR) then
      call detailedError(node, msg)
    elseif (iErr == TOKEN_EOS) then
      call detailedError(node, "Unexpected end of data")
    end if

  end subroutine checkError


  !> Issues an error, if the string from a given position contains non-whitespace characters.
  subroutine checkNoData(node, str, start)

    !> Node which is being processed (for error message)
    type(fnode), pointer :: node

    !> String content of the child.
    character(len=*), intent(in) :: str

    !> Starting position, after which the string should not contain any whitespace characters.
    integer, intent(in) :: start

    if (complementaryScan(str(start:), whiteSpaces) > 0) then
      call detailedError(node, "Superfluous data found.")
    end if

  end subroutine checkNoData


  !> Prints detailed error, including line number and path
  subroutine detailedError(node, msg)

    !> Node where the error occured.
    type(fnode), pointer :: node

    !> Message to print
    character(len=*), intent(in) :: msg

    type(string) :: str

    str = msg
    call appendPathAndLine(node, str)
    call error(char(str) // newline)

  end subroutine detailedError


  !> Prints detailed warning, including line number and path
  subroutine detailedWarning(node, msg)

    !> Node where the error occured.
    type(fnode), pointer :: node

    !> Message to print
    character(len=*), intent(in) :: msg

    type(string) :: str

    str = msg
    call appendPathAndLine(node, str)
    call warning(char(str) // newline)

  end subroutine detailedWarning


  !> Appends path and line information to a string.
  subroutine appendPathAndLine(node, str)

    !> Node, for which path and line should be added
    type(fnode), pointer :: node

    !> String prepending the path and line information
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

end module dftbp_hsdutils
