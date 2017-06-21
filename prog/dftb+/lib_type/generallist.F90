!--------------------------------------------------------------------------------------------------!
!  DFTB+: general package for performing fast atomistic simulations                                !
!  Copyright (C) 2017  DFTB+ developers group                                                      !
!                                                                                                  !
!  See the LICENSE file for terms of usage and distribution.                                       !
!--------------------------------------------------------------------------------------------------!

!!* Types and routines for the general purpose linked list.
module generallist
#include "assert.h"
#include "allocate.h"  
  implicit none
  private
  
  public :: OList, OListIterator, listData
  public :: init, length, append, getItem, index, isValid, next

  !!* List node (private).
  type OListNode
    private
    type(OListNode), pointer :: next => null()
    character, allocatable :: data(:)
  end type OListNode
  
  !!* General purpose linked list.
  type OList
    private
    type(OListNode), pointer :: first => null()
    type(OListNode), pointer :: last => null()
    integer :: size = 0
    integer :: dataSize
  contains
    final :: OList_destruct
  end type OList
  
  !!* Iterator for iterating through the elements of a list.
  type OListIterator
    private
    type(OListNode), pointer :: curNode
    integer :: dataSize
  end type OListIterator

  !!* Constant representing the data type stored by the list.
  character, parameter :: listData(1) = (/ '0' /)
  
  !!* Initialisation of the various derived types.
  interface init
    module procedure OList_init
    module procedure OListIterator_init
  end interface

  !!* Returns the length of the list.
  interface length
    module procedure OList_length
  end interface

  !!* Appends an element.
  interface append
    module procedure OList_append
  end interface

  !!* Returns an item.
  interface getItem
    module procedure OList_getItem
  end interface

  !!* Returns the index of an item.
  interface index
    module procedure OList_index
  end interface

  !!* Returns validity status.
  interface isValid
    module procedure OListIterator_isValid
  end interface
  
  !!* Returns the next element
  interface next
    module procedure OListIterator_next
  end interface
  

contains

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!! OList
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !!* Initialises the list.
  !!* @param sf List instance (self).
  !!* @param dataSize Data units needed to store the data for one element.
  subroutine OList_init(sf, dataSize)
    type(OList), intent(inout) :: sf
    integer, intent(in) :: dataSize

    sf%dataSize = dataSize

  end subroutine OList_init

  
  !!* Destructs the list.
  !!* @param sf List instance (self).
  subroutine OList_destruct(sf)
    type(OList), intent(inout) :: sf

    type(OListNode), pointer :: node, curNode

    node => sf%first
    do while (associated(node))
      curNode => node
      node => node%next
      deallocate(curNode)
    end do
    
  end subroutine OList_destruct


  !!* Return the length of the list.
  !!* @param sf List instance (self).
  function OList_length(sf) result(length)
    type(OList), intent(in) :: sf
    integer :: length

    length = sf%size

  end function OList_length

  
  !!* Appends an element to the list.
  !!* @param sf List instance (self).
  !!* @param data Data to append as a new element to the list.
  subroutine OList_append(sf, data)
    type(OList), intent(inout) :: sf
    character, intent(in) :: data(:)

    type(OListNode), pointer :: node

    ASSERT(size(data) == sf%dataSize)

    allocate(node)
    ALLOCATE_(node%data, (sf%dataSize))
    node%data = data
    
    if (sf%size == 0) then
      sf%first => node
      sf%last => node
    else
      sf%last%next => node
      sf%last => node
    end if
    sf%size = sf%size + 1

  end subroutine OList_append


  !!* Returns an item with a certain index.
  !!* The index must be between 1 and length of the list.
  !!* @param sf List instance (self).
  !!* @param index Element index.
  !!* @return Element with given index.
  function OList_getItem(sf, index) result(item)
    type(OList), intent(in) :: sf
    integer, intent(in) :: index
    character :: item(sf%dataSize)

    type(OListNode), pointer :: node
    integer :: ii

    ASSERT(index > 0 .or. index <= sf%size)

    node => sf%first
    do ii = 2, index
      node => node%next
    end do
    item = node%data

  end function OList_getItem


  !!* Returns the index for a given element.
  !!* @param sf List instance (self).
  !!* @param data Element to look for.
  !!* @return Index of the element or 0 if not found.
  function OList_index(sf, data) result(index)
    type(OList), intent(in) :: sf
    character :: data(:)
    integer :: index

    type(OListNode), pointer :: node
    integer :: ii

    index = 0
    node => sf%first
    do ii = 1, sf%size
      if (all(node%data == data)) then
        index = ii
        exit
      end if
      node => node%next
    end do
    
  end function OList_index

  
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!* OListIterator
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !!* Initialises the iterator.
  !!* @param sf Iterator instance (self).
  !!* @param list The list for which the iterator should be created.
  subroutine OListIterator_init(sf, list)
    type(OListIterator), intent(inout) :: sf
    type(OList), intent(in) :: list

    sf%curNode => list%first
    sf%dataSize = list%dataSize

  end subroutine OListIterator_init

  
  !!* Returns the validity of the iterator.
  !!* @param sf Iterator instance (self).
  !!* @return True, if iterator is valid (pointing to an existing element in
  !!*  the list), False otherwise.
  function OListIterator_isValid(sf) result(valid)
    type(OListIterator), intent(in) :: sf
    logical :: valid

    valid = associated(sf%curNode)

  end function OListIterator_isValid


  !!* Returns element which the iterator is pointing to and advances iterator.
  !!* @param sf Iterator instance (self).
  !!* @return Current element.
  !!* @note You should never invoke next(), unless isValid returned True.
  !!* @note Since next() also advances the iterator, isValid must be called
  !!* before every getItem call.
  function OListIterator_next(sf) result(item)
    type(OListIterator), intent(inout) :: sf
    character :: item(sf%dataSize)

    item = sf%curNode%data
    sf%curNode => sf%curNode%next

  end function OListIterator_next


end module generallist
