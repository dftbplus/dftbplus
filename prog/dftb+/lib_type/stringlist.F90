!--------------------------------------------------------------------------------------------------!
!  DFTB+: general package for performing fast atomistic simulations                                !
!  Copyright (C) 2017  DFTB+ developers group                                                      !
!                                                                                                  !
!  See the LICENSE file for terms of usage and distribution.                                       !
!--------------------------------------------------------------------------------------------------!

!!* Contains a linked list for string types.
!!* @desc The linked list for strings takes character variables with arbitary
!!*   length and stores those in a linked list. The elements of the lists
!!*   can be returned as character variables again, where the affected element
!!*   is eventually chopped, if the querying variable is not long enough.
!!* @note This linked list is separated from the general linkedList module
!!*   since it depends on the string module, which not all compiler are able
!!*   to compile.
module stringlist
#include "assert.h"
#include "allocate.h"  
  use xmlf90
  implicit none

  private

  !!* One node of the list.
  type nodeString
    type(string) :: value
    type(nodeString), pointer :: pNext
  end type nodeString


  !!* The list.
  type listString
    private
    integer                   :: length
    logical                   :: tUnishaped
    type(nodeString), pointer :: pFirst
    type(nodeString), pointer :: pLast
    integer                   :: iCache
    type(nodeString), pointer :: pCache
    logical                   :: tInitialized = .false.
  end type listString
  

  !!* Generic interface for initializing lists
  interface init
    module procedure initString
  end interface

  !!* Generic interface for destroying lists
  interface destroy
    module procedure destroyString
  end interface

  !!* Generic interface for appending elements to a list
  interface append
    module procedure appendString
  end interface

  !!* Generic interface for getting the length(nr. of elements) of a list
  interface len
    module procedure lenString
  end interface

  !!* Generic interface for getting the positon of an element in the list
  interface find
    module procedure findString
  end interface

  !!* Generic interface for checking if an element is in the list
  interface hasElement
    module procedure hasElementString
  end interface
  
  !!* Generic interface for getting an element with a given index
  interface get
    module procedure getString
  end interface

  !!* Generic interface for checking if all the list members have equal shape
  interface isUnishaped
    module procedure isUnishapedString
  end interface
  
  !!* Generic interface for getting the list as an array
  interface asArray
    module procedure asArrayString
  end interface


  public :: listString
  public :: init, destroy, append, len, find, hasElement, get, isUnishaped
  public :: asArray


contains

  
  !!* Initializes a list containing characters
  !!* @param list       The list to initialize.
  subroutine initString(list)
    type(listString), intent(inout) :: list

    ASSERT(.not. list%tInitialized)

    list%length = 0
    list%tUnishaped = .true.
    INIT_P(list%pFirst)
    INIT_P(list%pLast)
    list%iCache = 0
    INIT_P(list%pCache)
    list%tInitialized = .true.

  end subroutine initString



  !!* Destroys a list containing characters
  !!* @param list The list to destroy.
  subroutine destroyString(list)
    type(listString), intent(inout) :: list

    type(nodeString), pointer :: pCur, pNext

    ASSERT(list%tInitialized)

    pCur => list%pFirst
    do while(associated(pCur))
      call unstring(pCur%value)
      pNext => pCur%pNext
      DEALLOCATE_P(pCur)
      pCur => pNext
    end do
    list%tInitialized = .false.

  end subroutine destroyString



  !!* Appends an element to the list.
  !!* @param list  The list to extend.
  !!* @param value The value to add.
  subroutine appendString(list, value)
    type(listString), intent(inout) :: list
    character(len=*), intent(in)    :: value                ! type specific

    ASSERT(list%tInitialized)

    !! List contains already elements -> append to the end otherwise as first
    if(associated(list%pLast)) then
      INITALLOCATE_P(list%pLast%pNext)
      list%pLast => list%pLast%pNext
    else
      INITALLOCATE_P(list%pFirst)
      list%pLast => list%pFirst
    end if
    list%length = list%length + 1

    !! initialize node
    INIT_P(list%pLast%pNext)
    list%pLast%value = value

  end subroutine appendString



  !!* Returns the length(nr. of elements) of the list
  !!* @param list The list to get the length of.
  !!* @return     Nr. of elements in the list.
  integer function lenString(list) result(len)
    type(listString), intent(in) :: list
    ASSERT(list%tInitialized)
    len = list%length
  end function lenString



  !!* Returns the index of an element in the list.
  !!* @param list  The list object.
  !!* @param value The value to look for.
  !!* @return      Index of the element or zero if not found
  integer function findString(list, value)
    type(listString), intent(inout) :: list
    character(len=*), intent(in)    :: value

    type(nodeString), pointer :: pCur
    integer                   :: ii

    ASSERT(list%tInitialized)

    pCur => list%pFirst
    ii = 1
    do while(associated(pCur))
      if (pCur%value == value) then
        exit
      end if
      pCur => pCur%pNext
      ii = ii + 1
    end do

    if (associated(pCur)) then
      findString = ii
      list%iCache = ii
      list%pCache => pCur
    else
      findString = 0
    endif

  end function findString



  !!* Check if given element is in the list
  !!* @param list   The list object
  !!* @param value  Element to look for
  !!* @return       True if element had been found, false otherwise
  logical function hasElementString(list, value) result(hasElement)
    type(listString), intent(inout) :: list
    character(len=*), intent(in)    :: value

    ASSERT(list%tInitialized)

    if (find(list, value) == 0) then
      hasElement = .false.
    else
      hasElement = .true.
    end if

  end function hasElementString





  !!* Fills a variable with the speciesfied element of the list
  !!* @param list  The list object.
  !!* @param value The variable to put the element in.
  !!* @param index Index of the element (0 < index < length of the list)
  subroutine getString(list, value, index)
    type(listString), intent(inout) :: list
    character(len=*), intent(out)   :: value
    integer,          intent(in)    :: index

    type(nodeString), pointer :: pCur

    ASSERT(list%tInitialized)
    ASSERT(index > 0 .and. index <= list%length)

    pCur => getNodeString(list, index)
    value = pCur%value

  end subroutine getString


  !!* Checks if list contains members with equal shaped
  !!* @param list The list object.
  !!* @return     True, if elements have equals shaped, False otherwise.
  logical function isUnishapedString(list) result(isUnishaped)
    type(listString), intent(in) :: list
    ASSERT(list%tInitialized)
    isUnishaped = list%tUnishaped
  end function isUnishapedString



  !!* Returns the list as an array of elements.
  !!* @param list The list to get the elements from.
  !!* @param val  Array which will be filled with the elements of the list.
  !!* @note
  !!*   The passed array has to have the rank of the list elements + 1.
  !!*   According to Fortran traditions, the last index of the array addresses
  !!*   the list elements, the indexes before address the elements inside
  !!*   the list elements.
  !!* @note Only unishaped lists can be converted to array!
  !!* @assert Array has the shape(:, :, :, ..., :, <length of the list>)
  !!*         and the dimensions before the last one are compatible with the
  !!*         shape of the elements in the list.
  subroutine asArrayString(list, val, optPad)
    type(listString), intent(in)           :: list
    character(len=*), intent(out)          :: val(:)     ! type specific
    character(len=*), intent(in), optional :: optPad

    type(nodeString), pointer :: pCur
    integer                   :: lenVal
    integer                   :: ii

    lenVal = size(val, dim=size(shape(val)))

    ASSERT(list%tInitialized)
    ASSERT(lenVal >= list%length)

    pCur => list%pFirst
    ii = 1
    do while(associated(pCur))
      val(ii) = pCur%value                                  ! type specific
      ii = ii + 1
      pCur => pCur%pNext
    end do
    if (ii <= lenVal) then
      if (present(optPad)) then
        val(ii:lenVal) = optPad
      else
        val(ii:lenVal) = ""
      end if
    end if
    
  end subroutine asArrayString



  !!* Returns a pointer to a node with a given index
  !!* @param list  The list object.
  !!* @param pNode Pointer to set to the wanted node.
  !!* @param index Index of the wanted node.
  function getNodeString(list, index) result(getNode)
    type(nodeString), pointer                :: getNode
    type(listString),          intent(inout) :: list
    integer,                   intent(in)    :: index

    integer :: ii, iStart

    ASSERT(index > 0 .and. index <= list%length)

    if (list%iCache == index) then
      getNode => list%pCache
      return
    end if

    if (list%iCache > 0 .and. list%iCache < index) then
      iStart = list%iCache
      getNode => list%pCache
    else
      iStart = 1
      getNode => list%pFirst
    end if

    do ii = iStart + 1, index
      getNode => getNode%pNext
    end do
    list%pCache => getNode
    list%iCache = index

  end function getNodeString


end module stringlist
