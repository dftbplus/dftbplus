!--------------------------------------------------------------------------------------------------!
!  DFTB+: general package for performing fast atomistic simulations                                !
!  Copyright (C) 2017  DFTB+ developers group                                                      !
!                                                                                                  !
!  See the LICENSE file for terms of usage and distribution.                                       !
!--------------------------------------------------------------------------------------------------!

!!! -*- f90 -*-
!!!
!!!  Method definitions for linked list of type: Rank 1, Integer array
!!!


  !!* Initializes a list containing rank 1 arrays
  !!* @param list       The list to initialize.
  subroutine initIntR1(list)
    type(listIntR1), intent(inout) :: list

    ASSERT(.not. list%tInitialized)

    list%length = 0
    list%elemShape = (/ 0 /)
    list%tUnishaped = .true.
    nullify(list%pFirst)
    nullify(list%pLast)
    list%iCache = 0
    nullify(list%pCache)
    list%tInitialized = .true.

  end subroutine initIntR1

  

  !!* Destructs a list containing rank 1 arrays
  !!* @param list The list to destruct.
  subroutine destructIntR1(list)
    type(listIntR1), intent(inout) :: list

    type(nodeIntR1), pointer :: pCur, pNext

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

  end subroutine destructIntR1


  
  !!* Appends an element to the list.
  !!* @param list  The list to extend.
  !!* @param value The value to add.
  subroutine appendIntR1(list, value)
    type(listIntR1), intent(inout) :: list
    integer,         intent(in)    :: value(:)              ! type specific

    integer :: valshape(1)                                   ! type specific
    
    ASSERT(list%tInitialized)

    valshape = shape(value)

    !! List contains already elements -> append to the end otherwise as first
    !! check or set unishapeness
    if(associated(list%pLast)) then
      allocate(list%pLast%pNext)
      list%pLast => list%pLast%pNext
      if (list%tUnishaped .and. &
          &(.not. all(valshape(:) == list%elemShape(:)))) then
        list%tUnishaped = .false.
      end if
    else
      allocate(list%pFirst)
      list%pLast => list%pFirst
      list%elemShape(:) = valshape(:)
    end if
    list%length = list%length + 1

    !! initialize node
    nullify(list%pLast%pNext)
    ALLOCATE_(list%pLast%pValue, (valshape(1)))
    list%pLast%pValue(:) = value(:)
        
  end subroutine appendIntR1

  

  !!* Returns the length(nr. of elements) of the list
  !!* @param list The list to get the length of.
  !!* @return     Nr. of elements in the list.
  integer function lenIntR1(list) result(len)
    type(listIntR1), intent(in) :: list
    ASSERT(list%tInitialized)
    len = list%length
  end function lenIntR1



  !!* Returns the index of an element in the list.
  !!* @param list  The list object.
  !!* @param value The value to look for.
  !!* @return      Index of the element or zero if not found
  integer function findIntR1(list, value)
    type(listIntR1), intent(inout) :: list
    integer,         intent(in)    :: value(:)

    type(nodeIntR1), pointer :: pCur
    integer                  :: ii
    integer                  :: valshape(1)

    ASSERT(list%tInitialized)

    valshape = shape(value)
    pCur => list%pFirst
    ii = 1
    do while(associated(pCur))
      if (all(shape(pCur%pValue) == valshape)) then
        if (all(pCur%pValue(:) == value(:))) then
          exit
        end if
      end if
      pCur => pCur%pNext
      ii = ii + 1
    end do

    if (associated(pCur)) then
      findIntR1 = ii
      list%iCache = ii
      list%pCache => pCur
    else
      findIntR1 = 0
    endif

  end function findIntR1



  !!* Check if given element is in the list
  !!* @param list   The list object
  !!* @param value  Element to look for
  !!* @return       True if element had been found, false otherwise
  logical function hasElementIntR1(list, value) result(hasElement)
    type(listIntR1), intent(inout) :: list
    integer,         intent(in)    :: value(:)

    ASSERT(list%tInitialized)

    if (find(list, value) == 0) then
      hasElement = .false.
    else
      hasElement = .true.
    end if

  end function hasElementIntR1

  
  
  !!* Get the shape of the elements in the list
  !!* @param list     The list object
  !!* @param valshape contains the shape at return
  subroutine getElemShapeIntR1(list, valshape, index)
    type(listIntR1), intent(inout) :: list
    integer,          intent(out)   :: valshape(:)
    integer,          intent(in)    :: index

    type(nodeIntR1), pointer :: pCur

    ASSERT(list%tInitialized)
    ASSERT(index > 0 .and. index <= list%length)
    ASSERT(size(valshape, dim=1) == 1)

    if (list%tUnishaped) then
      valshape(:) = list%elemShape(:)
    else
      pCur => getNodeIntR1(list, index)
      valshape(:) = shape(pCur%pValue)
    end if

  end subroutine getElemShapeIntR1



  !!* Fills a variable with the speciefied element of the list
  !!* @param list  The list object.
  !!* @param value The variable to put the element in.
  !!* @param index Index of the element (0 < index < length of the list)
  subroutine getIntR1(list, value, index)
    type(listIntR1), intent(inout) :: list
    integer,         intent(out)   :: value(:)
    integer,         intent(in)    :: index

    type(nodeIntR1), pointer :: pCur

    ASSERT(list%tInitialized)
    ASSERT(index > 0 .and. index <= list%length)

    pCur => getNodeIntR1(list, index)
    ASSERT(all(shape(pCur%pValue) == shape(value)))
    value(:) = pCur%pValue(:)

  end subroutine getIntR1

  

  !!* Fills the speciefied element of the list with a variable
  !!* @param list  The list object.
  !!* @param value The variable to put the element in.
  !!* @param index Index of the element (0 < index < length of the list)
  subroutine setIntR1(list, value, index)
    type(listIntR1), intent(inout) :: list
    integer,         intent(in)   :: value(:)
    integer,         intent(in)    :: index

    type(nodeIntR1), pointer :: pCur

    ASSERT(list%tInitialized)
    ASSERT(index > 0 .and. index <= list%length)

    pCur => getNodeIntR1(list, index)
    ASSERT(all(shape(pCur%pValue) == shape(value)))
    pCur%pValue(:) = value(:)

  end subroutine setIntR1
  
  

  !!* Checks if list contains members with equal shaped
  !!* @param list The list object.
  !!* @return     True, if elements have equals shaped, False otherwise.
  logical function isUnishapedIntR1(list) result(isUnishaped)
    type(listIntR1), intent(in) :: list

    ASSERT(list%tInitialized)

    isUnishaped = list%tUnishaped

  end function isUnishapedIntR1

  

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
  subroutine asArrayIntR1(list, val, optPad)
    type(listIntR1), intent(in)           :: list
    integer,         intent(out)          :: val(:, :)      ! type specific
    integer,         intent(in), optional :: optPad

    type(nodeIntR1), pointer :: pCur
    integer                  :: lenVal
    integer                  :: pad
    integer                  :: ii

    ASSERT_ENV(integer       :: shapeVal(2))

    lenVal = size(val, dim=size(shape(val)))

    ASSERT(list%tInitialized)
    ASSERT(list%tUnishaped)
    ASSERT(lenVal >= list%length)
    ASSERT_ENV(shapeVal = shape(val))
    ASSERT(all(shapeVal(1:size(shapeVal)-1) == list%elemShape(:)))

    if (present(optPad)) then
      pad = optPad
    else
      pad = 0
    end if

    pCur => list%pFirst
    ii = 1
    do while(associated(pCur))
      val(:, ii) = pCur%pValue(:)                              ! type specific
      ii = ii + 1
      pCur => pCur%pNext
    end do
    if (ii <= lenVal) then
      val(:,ii:lenVal) = pad
    end if

  end subroutine asArrayIntR1



  !!* Returns the list as a vector (rank 1 array).
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
  subroutine asVectorIntR1(list, val, optPad)
    type(listIntR1), intent(in)           :: list
    integer,         intent(out)          :: val(:)      ! type specific
    integer,         intent(in), optional :: optPad

    type(nodeIntR1), pointer :: pCur
    integer                  :: lenVal, lenElem
    integer                  :: pad
    integer                  :: ii

    lenVal = size(val)
    lenElem = product(list%elemShape)

    ASSERT(list%tInitialized)
    ASSERT(list%tUnishaped)
    ASSERT(lenVal >= list%length * lenElem)

    if (present(optPad)) then
      pad = optPad
    else
      pad = 0
    end if
    
    pCur => list%pFirst
    ii = 1
    do while(associated(pCur))
      val(ii:ii+lenElem-1) = pCur%pValue(:)    ! type specific
      ii = ii + lenElem
      pCur => pCur%pNext
    end do
    if (ii <= lenVal) then
      val(ii:) = pad
    end if
    
  end subroutine asVectorIntR1

  

  !!* Fills a variable with the speciefied element of the list
  !!* @param list  The list object.
  !!* @param value The variable to put the element in.
  !!* @param nElem number of elements returned in value
  !!* @param index Index of the element (0 < index < length of the list)
  subroutine intoArrayIntR1(list, value, nElem, index)
    type(listIntR1), intent(inout) :: list
    integer, intent(inout) :: value(:)
    integer, intent(out) :: nElem
    integer, intent(in) :: index

    type(nodeIntR1), pointer :: pCur

    ASSERT(list%tInitialized)
    ASSERT(index > 0 .and. index <= list%length)

    pCur => getNodeIntR1(list, index)
    nElem = size(pCur%pValue)
    ASSERT(size(value) >= nElem)
    value(:nElem) = pCur%pValue(:)

  end subroutine intoArrayIntR1



  !!* Returns a pointer to a node with a given index
  !!* @param list  The list object.
  !!* @param pNode Pointer to set to the wanted node.
  !!* @param index Index of the wanted node.
  function getNodeIntR1(list, index) result(getNode)
    type(nodeIntR1), pointer       :: getNode
    type(listIntR1), intent(inout) :: list
    integer,         intent(in)    :: index

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
    
   end function getNodeIntR1
   
  
