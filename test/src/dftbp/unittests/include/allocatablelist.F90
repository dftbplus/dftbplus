!--------------------------------------------------------------------------------------------------!
!  DFTB+: general package for performing fast atomistic simulations                                !
!  Copyright (C) 2006 - 2023  DFTB+ developers group                                               !
!                                                                                                  !
!  See the LICENSE file for terms of usage and distribution.                                       !
!--------------------------------------------------------------------------------------------------!

#:include "fytest.fypp"


#:block TEST_SUITE("allocatablelist")
  use alloclisthelper, only : TIntList
  implicit none

#:contains

  #:block TEST("uninitialized")
    type(TIntList) :: list

    @:ASSERT(list%size() == 0)
  #:endblock TEST


  #:block TEST("pushAndPopItems")

    type(TIntList) :: list
    integer, allocatable :: item
    integer, pointer :: pItem
    integer :: ii

    @:ASSERT(list%size() == 0)
    do ii = 1, 5
      item = ii
      call list%push(item)
      @:ASSERT(list%size() == ii)
      @:ASSERT(.not. allocated(item))
    end do

    do ii = 1, 5
      pItem => null()
      call list%view(ii, pItem)
      @:ASSERT(associated(pItem))
      @:ASSERT(pItem == ii)
      @:ASSERT(list%size() == 5)
    end do

    do ii = 5, 1, -1
      call list%pop(item)
      @:ASSERT(allocated(item))
      @:ASSERT(item == ii)
      @:ASSERT(list%size() == ii - 1)
    end do

  #:endblock TEST


  #:block TEST_FIXTURE("smallList")
    type(TIntList), target :: list
    integer, parameter :: initSize = 10

  #:contains

    #:block TEST_FIXTURE_INIT()
      integer, allocatable :: item
      integer :: ii

      do ii = 1, initSize
        item = ii
        call list%push(item)
      end do
    #:endblock TEST_FIXTURE_INIT


    #:block TEST("pushFront")
      integer, allocatable :: item
      integer, pointer :: pItem

      allocate(item)
      item = 0
      call list%push(item, pos=1)
      @:ASSERT(list%size() == initSize + 1)
      call list%view(2, pItem)
      @:ASSERT(pItem == 1)
      call list%view(1, pItem)
      @:ASSERT(pItem == 0)
    #:endblock TEST


    #:block TEST("pushFrontNegative")
      integer, allocatable :: item
      integer, pointer :: pItem

      item = 0
      call list%push(item, pos=-initSize)
      @:ASSERT(list%size() == initSize + 1)
      call list%view(2, pItem)
      @:ASSERT(pItem == 1)
      call list%view(1, pItem)
      @:ASSERT(pItem == 0)
    #:endblock TEST


    #:block TEST("pushMiddle")
      integer, allocatable :: item
      integer, pointer :: pItem
      integer :: pos

      pos = int(initSize / 2)
      item = 0
      call list%push(item, pos=pos)
      @:ASSERT(list%size() == initSize + 1)
      call list%view(pos + 1, pItem)
      @:ASSERT(pItem == pos)
      call list%view(pos, pItem)
      @:ASSERT(pItem == 0)
      if (pos > 1) then
        call list%view(pos - 1, pItem)
        @:ASSERT(pItem == pos - 1)
      end if
    #:endblock TEST


    #:block TEST("pushMiddleNegative")
      integer, allocatable :: item
      integer, pointer :: pItem
      integer :: pos

      pos = int(initSize / 2)
      item = 0
      call list%push(item, pos=-pos)
      @:ASSERT(list%size() == initSize + 1)
      call list%view(-pos, pItem)
      @:ASSERT(pItem == initSize + 1 - pos)
      call list%view(-pos - 1, pItem)
      @:ASSERT(pItem == 0)
    #:endblock TEST


    #:block TEST("pushBack")
      integer, allocatable :: item
      integer, pointer :: pItem
      integer :: pos

      pos = initSize + 1
      item = 0
      call list%push(item, pos=pos)
      @:ASSERT(list%size() == initSize + 1)
      call list%view(pos, pItem)
      @:ASSERT(pItem == 0)
    #:endblock TEST


    #:block TEST("pushBackZero")
      integer, allocatable :: item
      integer, pointer :: pItem
      integer :: pos

      pos = 0
      item = 0
      call list%push(item, pos=pos)
      @:ASSERT(list%size() == initSize + 1)
      call list%view(-1, pItem)
      @:ASSERT(pItem == 0)
      call list%view(initSize + 1, pItem)
      @:ASSERT(pItem == 0)
    #:endblock TEST


    #:block TEST("popFront")
      integer, allocatable :: item
      integer, pointer :: pItem

      call list%pop(item, pos=1)
      @:ASSERT(allocated(item))
      @:ASSERT(item == 1)
      @:ASSERT(list%size() == initSize - 1)
      call list%view(1, pItem)
      @:ASSERT(pItem == 2)
    #:endblock TEST


    #:block TEST("popFrontNegative")
      integer, allocatable :: item
      integer, pointer :: pItem

      call list%pop(item, pos=-initSize)
      @:ASSERT(allocated(item))
      @:ASSERT(item == 1)
      @:ASSERT(list%size() == initSize - 1)
      call list%view(1, pItem)
      @:ASSERT(pItem == 2)
    #:endblock TEST


    #:block TEST("popMiddle")
      integer, allocatable :: item
      integer, pointer :: pItem
      integer :: pos

      pos = int(initSize / 2)
      call list%pop(item, pos=pos)
      @:ASSERT(allocated(item))
      @:ASSERT(item == pos)
      @:ASSERT(list%size() == initSize - 1)
      call list%view(pos, pItem)
      @:ASSERT(pItem == pos + 1)
    #:endblock TEST


    #:block TEST("popMiddleNegative")
      integer, allocatable :: item
      integer, pointer :: pItem
      integer :: pos

      pos = int(initSize / 2)
      call list%pop(item, pos=-pos)
      @:ASSERT(allocated(item))
      @:ASSERT(item == initSize + 1 - pos)
      @:ASSERT(list%size() == initSize - 1)
      call list%view(-pos + 1, pItem)
      @:ASSERT(pItem == initSize + 1 - pos + 1)
    #:endblock TEST


    #:block TEST("popBack")
      integer, allocatable :: item
      integer, pointer :: pItem
      integer :: pos

      pos = initSize
      call list%pop(item, pos=pos)
      @:ASSERT(allocated(item))
      @:ASSERT(item == initSize)
      @:ASSERT(list%size() == initSize - 1)
      call list%view(pos - 1, pItem)
      @:ASSERT(pItem == pos - 1)
    #:endblock TEST


    #:block TEST("popBackNegative")
      integer, allocatable :: item
      integer, pointer :: pItem
      integer :: pos

      pos = -1
      call list%pop(item, pos=pos)
      @:ASSERT(allocated(item))
      @:ASSERT(item == initSize)
      @:ASSERT(list%size() == initSize - 1)
      call list%view(-1, pItem)
      @:ASSERT(pItem == initSize - 1)
      call list%view(initSize - 1, pItem)
      @:ASSERT(pItem == initSize - 1)
    #:endblock TEST

  #:endblock TEST_FIXTURE

#:endblock TEST_SUITE


@:TEST_DRIVER()
