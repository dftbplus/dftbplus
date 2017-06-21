!--------------------------------------------------------------------------------------------------!
!  DFTB+: general package for performing fast atomistic simulations                                !
!  Copyright (C) 2017  DFTB+ developers group                                                      !
!                                                                                                  !
!  See the LICENSE file for terms of usage and distribution.                                       !
!--------------------------------------------------------------------------------------------------!

!!! -*- f90 -*-
!!!
!!!  Method definitions for linked list of type: reals
!!!


  !!* Initializes the list.
  !!* @param list The list to initialize.
  subroutine initReal(list)
    type(listReal), intent(inout) :: list

    ASSERT(.not. list%tInitialized)

    list%length = 0
    list%tUnishaped = .true.
    nullify(list%pFirst)
    nullify(list%pLast)
    list%iCache = 0
    nullify(list%pCache)
    list%tInitialized = .true.

  end subroutine initReal

  

  !!* Destructs the list.
  !!* @param list The list to destruct.
  subroutine destructReal(list)
    type(listReal), intent(inout) :: list

    type(nodeReal), pointer :: pCur, pNext

    if (.not. list%tInitialized) then
      return
    end if

    pCur => list%pFirst
    do while(associated(pCur))
      pNext => pCur%pNext
      deallocate(pCur)
      pCur => pNext
    end do
    list%tInitialized = .false.

  end subroutine destructReal


  
  !!* Appends an element to the list.
  !!* @param list  The list to extend.
  !!* @param value The value to add.
  subroutine appendReal(list, value)
    type(listReal), intent(inout) :: list
    real(dp),       intent(in)    :: value                ! type specific

    ASSERT(list%tInitialized)

    !! List contains already elements -> append to the end otherwise as first
    if(associated(list%pLast)) then
      allocate(list%pLast%pNext)
      list%pLast => list%pLast%pNext
    else
      allocate(list%pFirst)
      list%pLast => list%pFirst
    end if
    list%length = list%length + 1

    !! initialize node
    nullify(list%pLast%pNext)
    list%pLast%value = value
        
  end subroutine appendReal

  

  !!* Returns the length(nr. of elements) of the list
  !!* @param list The list to get the length of.
  !!* @return     Nr. of elements in the list.
  integer function lenReal(list) result(len)
    type(listReal), intent(in) :: list
    ASSERT(list%tInitialized)
    len = list%length
  end function lenReal



  !!* Returns the index of an element in the list.
  !!* @param list  The list object.
  !!* @param value The value to look for.
  !!* @return      Index of the element or zero if not found
  integer function findReal(list, value)
    type(listReal), intent(inout) :: list
    real(dp),       intent(in)    :: value

    type(nodeReal), pointer :: pCur
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
      findReal = ii
      list%iCache = ii
      list%pCache => pCur
    else
      findReal = 0
    endif

  end function findReal



  !!* Check if given element is in the list
  !!* @param list   The list object
  !!* @param value  Element to look for
  !!* @return       True if element had been found, false otherwise
  logical function hasElementReal(list, value) result(hasElement)
    type(listReal), intent(inout) :: list
    real(dp),    intent(in)    :: value

    ASSERT(list%tInitialized)

    if (find(list, value) == 0) then
      hasElement = .false.
    else
      hasElement = .true.
    end if

  end function hasElementReal
  
  


  !!* Fills a variable with the speciefied element of the list
  !!* @param list  The list object.
  !!* @param value The variable to put the element in.
  !!* @param index Index of the element (0 < index < length of the list)
  subroutine getReal(list, value, index)
    type(listReal), intent(inout) :: list
    real(dp),       intent(out)   :: value
    integer,        intent(in)    :: index

    type(nodeReal), pointer :: pCur

    ASSERT(list%tInitialized)
    ASSERT(index > 0 .and. index <= list%length)

    pCur => getNodeReal(list, index)
    value = pCur%value

  end subroutine getReal
  


  !!* Fills the speciefied element of the list with a variable.
  !!* @param list  The list object.
  !!* @param value The variable to put the element in.
  !!* @param index Index of the element (0 < index < length of the list)
  subroutine setReal(list, value, index)
    type(listReal), intent(inout) :: list
    real(dp),       intent(in)   :: value
    integer,        intent(in)    :: index

    type(nodeReal), pointer :: pCur

    ASSERT(list%tInitialized)
    ASSERT(index > 0 .and. index <= list%length)

    pCur => getNodeReal(list, index)
    pCur%value = value

  end subroutine setReal
  


  !!* Checks if list contains members with equal shaped
  !!* @param list The list object.
  !!* @return     True, if elements have equals shaped, False otherwise.
  logical function isUnishapedReal(list) result(isUnishaped)
    type(listReal), intent(in) :: list
    ASSERT(list%tInitialized)
    isUnishaped = list%tUnishaped
  end function isUnishapedReal

  

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
  subroutine asArrayReal(list, val, optPad)
    type(listReal), intent(in)           :: list
    real(dp),       intent(out)          :: val(:)     ! type specific
    real(dp),       intent(in), optional :: optPad

    type(nodeReal), pointer :: pCur
    integer                :: lenVal
    real(dp)               :: pad
    integer                :: ii

    lenVal = size(val, dim=size(shape(val)))
    
    ASSERT(list%tInitialized)
    ASSERT(lenVal >= list%length)

    if (present(optPad)) then
      pad = optPad
    else
      pad = 0
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
    
  end subroutine asArrayReal



  !!* Returns a pointer to a node with a given index
  !!* @param list  The list object.
  !!* @param pNode Pointer to set to the wanted node.
  !!* @param index Index of the wanted node.
  function getNodeReal(list, index) result(getNode)
    type(nodeReal), pointer       :: getNode
    type(listReal), intent(inout) :: list
    integer,       intent(in)    :: index

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

   end function getNodeReal
