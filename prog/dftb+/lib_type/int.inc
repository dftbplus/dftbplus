!--------------------------------------------------------------------------------------------------!
!  DFTB+: general package for performing fast atomistic simulations                                !
!  Copyright (C) 2017  DFTB+ developers group                                                      !
!                                                                                                  !
!  See the LICENSE file for terms of usage and distribution.                                       !
!--------------------------------------------------------------------------------------------------!

!!! -*- f90 -*-
!!!
!!!  Method definitions for linked list of type: integers
!!!


  !!* Initializes the list.
  !!* @param list The list to initialize.
  subroutine initInt(list)
    type(listInt), intent(inout) :: list

    @:ASSERT(.not. list%tInitialized)

    list%length = 0
    list%tUnishaped = .true.
    nullify(list%pFirst)
    nullify(list%pLast)
    list%iCache = 0
    nullify(list%pCache)
    list%tInitialized = .true.

  end subroutine initInt

  

  !!* Destructs the list.
  !!* @param list The list to destruct.
  subroutine destructInt(list)
    type(listInt), intent(inout) :: list

    type(nodeInt), pointer :: pCur, pNext

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

  end subroutine destructInt


  
  !!* Appends an element to the list.
  !!* @param list  The list to extend.
  !!* @param value The value to add.
  subroutine appendInt(list, value)
    type(listInt), intent(inout) :: list
    integer,       intent(in)    :: value                ! type specific

    @:ASSERT(list%tInitialized)

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
        
  end subroutine appendInt

  

  !!* Returns the length(nr. of elements) of the list
  !!* @param list The list to get the length of.
  !!* @return     Nr. of elements in the list.
  integer function lenInt(list) result(len)
    type(listInt), intent(in) :: list
    @:ASSERT(list%tInitialized)
    len = list%length
  end function lenInt



  !!* Returns the index of an element in the list.
  !!* @param list  The list object.
  !!* @param value The value to look for.
  !!* @return      Index of the element or zero if not found
  integer function findInt(list, value)
    type(listInt), intent(inout) :: list
    integer,       intent(in)    :: value

    type(nodeInt), pointer :: pCur
    integer                   :: ii

    @:ASSERT(list%tInitialized)

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
      findInt = ii
      list%iCache = ii
      list%pCache => pCur
    else
      findInt = 0
    endif

  end function findInt



  !!* Check if given element is in the list
  !!* @param list   The list object
  !!* @param value  Element to look for
  !!* @return       True if element had been found, false otherwise
  logical function hasElementInt(list, value) result(hasElement)
    type(listInt), intent(inout) :: list
    integer,    intent(in)    :: value

    @:ASSERT(list%tInitialized)

    if (find(list, value) == 0) then
      hasElement = .false.
    else
      hasElement = .true.
    end if

  end function hasElementInt

  
  
  !!* Fills a variable with the speciefied element of the list
  !!* @param list  The list object.
  !!* @param value The variable to put the element in.
  !!* @param index Index of the element (0 < index < length of the list)
  subroutine getInt(list, value, index)
    type(listInt), intent(inout) :: list
    integer,       intent(out)   :: value
    integer,       intent(in)    :: index

    type(nodeInt), pointer :: pCur

    @:ASSERT(list%tInitialized)
    @:ASSERT(index > 0 .and. index <= list%length)

    pCur => getNodeInt(list, index)
    value = pCur%value

  end subroutine getInt
  


  !!* Fills a speciefied element of the list with a variable
  !!* @param list  The list object.
  !!* @param value The variable to put the element in.
  !!* @param index Index of the element (0 < index < length of the list)
  subroutine setInt(list, value, index)
    type(listInt), intent(inout) :: list
    integer,       intent(in)   :: value
    integer,       intent(in)    :: index

    type(nodeInt), pointer :: pCur

    @:ASSERT(list%tInitialized)
    @:ASSERT(index > 0 .and. index <= list%length)

    pCur => getNodeInt(list, index)
    pCur%value = value

  end subroutine setInt
  
  
  
  !!* Checks if list contains members with equal shaped
  !!* @param list The list object.
  !!* @return     True, if elements have equals shaped, False otherwise.
  logical function isUnishapedInt(list) result(isUnishaped)
    type(listInt), intent(in) :: list
    @:ASSERT(list%tInitialized)
    isUnishaped = list%tUnishaped
  end function isUnishapedInt

  

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
  subroutine asArrayInt(list, val, optPad)
    type(listInt), intent(in)           :: list
    integer,       intent(out)          :: val(:)     ! type specific
    integer,       intent(in), optional :: optPad

    type(nodeInt), pointer :: pCur
    integer                :: lenVal
    integer                :: pad
    integer                :: ii

    lenVal = size(val, dim=size(shape(val)))
    
    @:ASSERT(list%tInitialized)
    @:ASSERT(lenVal >= list%length)

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
    
  end subroutine asArrayInt



  !!* Returns a pointer to a node with a given index
  !!* @param list  The list object.
  !!* @param pNode Pointer to set to the wanted node.
  !!* @param index Index of the wanted node.
  function getNodeInt(list, index) result(getNode)
    type(nodeInt), pointer       :: getNode
    type(listInt), intent(inout) :: list
    integer,       intent(in)    :: index

    integer :: ii, iStart
    
    @:ASSERT(index > 0 .and. index <= list%length)

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

   end function getNodeInt
