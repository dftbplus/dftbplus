!--------------------------------------------------------------------------------------------------!
!  DFTB+: general package for performing fast atomistic simulations                                !
!  Copyright (C) 2017  DFTB+ developers group                                                      !
!                                                                                                  !
!  See the LICENSE file for terms of usage and distribution.                                       !
!--------------------------------------------------------------------------------------------------!

!!! -*- f90 -*-
!!!
!!!  Method definitions for linked list of type: character (strings)
!!!


  !!* Initializes a list containing characters
  !!* @param list       The list to initialize.
  subroutine initCharMc(list)
    type(listCharMc), intent(inout) :: list

    ASSERT(.not. list%tInitialized)

    list%length = 0
    list%tUnishaped = .true.
    INIT_P(list%pFirst)
    INIT_P(list%pLast)
    list%iCache = 0
    INIT_P(list%pCache)
    list%tInitialized = .true.

  end subroutine initCharMc

  

  !!* Destroys a list containing characters
  !!* @param list The list to destroy.
  subroutine destroyCharMc(list)
    type(listCharMc), intent(inout) :: list

    type(nodeCharMc), pointer :: pCur, pNext

    ASSERT(list%tInitialized)

    pCur => list%pFirst
    do while(associated(pCur))
      pNext => pCur%pNext
      DEALLOCATE_P(pCur)
      pCur => pNext
    end do
    list%tInitialized = .false.

  end subroutine destroyCharMc


  
  !!* Appends an element to the list.
  !!* @param list  The list to extend.
  !!* @param value The value to add.
  subroutine appendCharMc(list, value)
    type(listCharMc), intent(inout) :: list
    character(mc),    intent(in)    :: value                ! type specific

    ASSERT(list%tInitialized)

    !! List contains already elements -> append to the end otherwise as first
    if(associated(list%pLast)) then
      ALLOCATE_P(list%pLast%pNext)
      list%pLast => list%pLast%pNext
    else
      ALLOCATE_P(list%pFirst)
      list%pLast => list%pFirst
    end if
    list%length = list%length + 1

    !! initialize node
    INIT_P(list%pLast%pNext)
    list%pLast%value = value
        
  end subroutine appendCharMc

  

  !!* Returns the length(nr. of elements) of the list
  !!* @param list The list to get the length of.
  !!* @return     Nr. of elements in the list.
  integer function lenCharMc(list) result(len)
    type(listCharMc), intent(in) :: list
    ASSERT(list%tInitialized)
    len = list%length
  end function lenCharMc



  !!* Returns the index of an element in the list.
  !!* @param list  The list object.
  !!* @param value The value to look for.
  !!* @return      Index of the element or zero if not found
  integer function findCharMc(list, value)
    type(listCharMc), intent(inout) :: list
    character(mc),    intent(in)    :: value

    type(nodeCharMc), pointer :: pCur
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
      findCharMc = ii
      list%iCache = ii
      list%pCache => pCur
    else
      findCharMc = 0
    endif

  end function findCharMc



  !!* Check if given element is in the list
  !!* @param list   The list object
  !!* @param value  Element to look for
  !!* @return       True if element had been found, false otherwise
  logical function hasElementCharMc(list, value) result(hasElement)
    type(listCharMc), intent(inout) :: list
    character(mc),    intent(in)    :: value

    ASSERT(list%tInitialized)

    if (find(list, value) == 0) then
      hasElement = .false.
    else
      hasElement = .true.
    end if

  end function hasElementCharMc

  
  
  !!* Fills a variable with the speciefied element of the list
  !!* @param list  The list object.
  !!* @param value The variable to put the element in.
  !!* @param index Index of the element (0 < index < length of the list)
  subroutine getCharMc(list, value, index)
    type(listCharMc), intent(inout) :: list
    character(mc),    intent(out)   :: value
    integer,          intent(in)    :: index

    type(nodeCharMc), pointer :: pCur

    ASSERT(list%tInitialized)
    ASSERT(index > 0 .and. index <= list%length)

    pCur => getNodeCharMc(list, index)
    value = pCur%value

  end subroutine getCharMc
  

  
  !!* Fills the speciefied element of the list with a variable.
  !!* @param list  The list object.
  !!* @param value The variable to put the element in.
  !!* @param index Index of the element (0 < index < length of the list)
  subroutine setCharMc(list, value, index)
    type(listCharMc), intent(inout) :: list
    character(mc),    intent(in)   :: value
    integer,          intent(in)    :: index

    type(nodeCharMc), pointer :: pCur

    ASSERT(list%tInitialized)
    ASSERT(index > 0 .and. index <= list%length)

    pCur => getNodeCharMc(list, index)
    pCur%value = value

  end subroutine setCharMc
  


  !!* Checks if list contains members with equal shaped
  !!* @param list The list object.
  !!* @return     True, if elements have equals shaped, False otherwise.
  logical function isUnishapedCharMc(list) result(isUnishaped)
    type(listCharMc), intent(in) :: list
    ASSERT(list%tInitialized)
    isUnishaped = list%tUnishaped
  end function isUnishapedCharMc

  

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
  subroutine asArrayCharMc(list, val, optPad)
    type(listCharMc), intent(in)           :: list
    character(mc),    intent(out)          :: val(:)     ! type specific
    character(mc),    intent(in), optional :: optPad

    type(nodeCharMc), pointer :: pCur
    integer                   :: lenVal
    character(mc)             :: pad
    integer                   :: ii

    lenVal = size(val, dim=size(shape(val)))
    
    ASSERT(list%tInitialized)
    ASSERT(lenVal >= list%length)

    if (present(optPad)) then
      pad = optPad
    else
      pad = ''
    end if
    
    pCur => list%pFirst
    ii = 1
    do while(associated(pCur))
      val(ii) = pCur%value                                  ! type specific
      ii = ii + 1
      pCur => pCur%pNext
    end do
    if (ii <= lenVal) then
      val(ii:lenVal) = pad
    end if
    
  end subroutine asArrayCharMc



  !!* Returns a pointer to a node with a given index
  !!* @param list  The list object.
  !!* @param pNode Pointer to set to the wanted node.
  !!* @param index Index of the wanted node.
  function getNodeCharMc(list, index) result(getNode)
    type(nodeCharMc), pointer                :: getNode
    type(listCharMc),          intent(inout) :: list
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

   end function getNodeCharMc
