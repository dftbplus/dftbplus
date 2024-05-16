!--------------------------------------------------------------------------------------------------!
!  DFTB+: general package for performing fast atomistic simulations                                !
!  Copyright (C) 2006 - 2024  DFTB+ developers group                                               !
!                                                                                                  !
!  See the LICENSE file for terms of usage and distribution.                                       !
!--------------------------------------------------------------------------------------------------!

#:include "fytest.fypp"

#:block TEST_SUITE("pointerlist")
  use ptrlisthelper, only : TInt, TIntList
  implicit none

#:contains


  #:block TEST("uninitialized")
    type(TIntList) :: list

    @:ASSERT(list%size() == 0)
  #:endblock TEST


  #:block TEST("pushAndPopItems")

    type(TIntList) :: list
    type(TInt), pointer :: item
    integer :: ii

    @:ASSERT(list%size() == 0)
    do ii = 1, 5
      allocate(item)
      item%value = ii
      call list%push(item)
      @:ASSERT(list%size() == ii)
      @:ASSERT(.not. associated(item))
    end do

    do ii = 1, 5
      item => null()
      call list%view(ii, item)
      @:ASSERT(associated(item))
      @:ASSERT(item%value == ii)
      @:ASSERT(list%size() == 5)
    end do

    do ii = 5, 1, -1
      item => null()
      call list%pop(item)
      @:ASSERT(associated(item))
      @:ASSERT(item%value == ii)
      @:ASSERT(list%size() == ii - 1)
    end do

  #:endblock TEST


  #:block TEST_FIXTURE("smallList")
    type(TIntList) :: list
    integer, parameter :: initSize = 10

  #:contains

    #:block TEST_FIXTURE_INIT()
      type(TInt), pointer :: item
      integer :: ii

      do ii = 1, initSize
        allocate(item)
        item%value = ii
        call list%push(item)
      end do
    #:endblock TEST_FIXTURE_INIT


    #:block TEST("pushFront")
      type(TInt), pointer :: item

      allocate(item)
      item%value = 0
      call list%push(item, pos=1)
      @:ASSERT(list%size() == initSize + 1)
      call list%view(2, item)
      @:ASSERT(item%value == 1)
      call list%view(1, item)
      @:ASSERT(item%value == 0)
    #:endblock TEST


    #:block TEST("pushFrontNegative")
      type(TInt), pointer :: item

      allocate(item)
      item%value = 0
      call list%push(item, pos=-initSize)
      @:ASSERT(list%size() == initSize + 1)
      call list%view(2, item)
      @:ASSERT(item%value == 1)
      call list%view(1, item)
      @:ASSERT(item%value == 0)
    #:endblock TEST


    #:block TEST("pushMiddle")
      type(TInt), pointer :: item
      integer :: pos

      pos = int(initSize / 2)
      allocate(item)
      item%value = 0
      call list%push(item, pos=pos)
      @:ASSERT(list%size() == initSize + 1)
      call list%view(pos + 1, item)
      @:ASSERT(item%value == pos)
      call list%view(pos, item)
      @:ASSERT(item%value == 0)
      if (pos > 1) then
        call list%view(pos - 1, item)
        @:ASSERT(item%value == pos - 1)
      end if
    #:endblock TEST


    #:block TEST("pushMiddleNegative")
      type(TInt), pointer :: item
      integer :: pos

      pos = int(initSize / 2)
      allocate(item)
      item%value = 0
      call list%push(item, pos=-pos)
      @:ASSERT(list%size() == initSize + 1)
      call list%view(-pos, item)
      @:ASSERT(item%value == initSize + 1 - pos)
      call list%view(-pos - 1, item)
      @:ASSERT(item%value == 0)
    #:endblock TEST


    #:block TEST("pushBack")
      type(TInt), pointer :: item
      integer :: pos

      pos = initSize + 1
      allocate(item)
      item%value = 0
      call list%push(item, pos=pos)
      @:ASSERT(list%size() == initSize + 1)
      call list%view(pos, item)
      @:ASSERT(item%value == 0)
    #:endblock TEST


    #:block TEST("pushBackZero")
      type(TInt), pointer :: item
      integer :: pos

      pos = 0
      allocate(item)
      item%value = 0
      call list%push(item, pos=pos)
      @:ASSERT(list%size() == initSize + 1)
      call list%view(-1, item)
      @:ASSERT(item%value == 0)
      call list%view(initSize + 1, item)
      @:ASSERT(item%value == 0)
    #:endblock TEST


    #:block TEST("popFront")
      type(TInt), pointer :: item

      item => null()
      call list%pop(item, pos=1)
      @:ASSERT(associated(item))
      @:ASSERT(item%value == 1)
      deallocate(item)
      @:ASSERT(list%size() == initSize - 1)
      call list%view(1, item)
      @:ASSERT(item%value == 2)
    #:endblock TEST


    #:block TEST("popFrontNegative")
      type(TInt), pointer :: item

      item => null()
      call list%pop(item, pos=-initSize)
      @:ASSERT(associated(item))
      @:ASSERT(item%value == 1)
      deallocate(item)
      @:ASSERT(list%size() == initSize - 1)
      call list%view(1, item)
      @:ASSERT(item%value == 2)
    #:endblock TEST


    #:block TEST("popMiddle")
      type(TInt), pointer :: item
      integer :: pos

      pos = int(initSize / 2)
      item => null()
      call list%pop(item, pos=pos)
      @:ASSERT(associated(item))
      @:ASSERT(item%value == pos)
      deallocate(item)
      @:ASSERT(list%size() == initSize - 1)
      call list%view(pos, item)
      @:ASSERT(item%value == pos + 1)
    #:endblock TEST


    #:block TEST("popMiddleNegative")
      type(TInt), pointer :: item
      integer :: pos

      pos = int(initSize / 2)
      item => null()
      call list%pop(item, pos=-pos)
      @:ASSERT(associated(item))
      @:ASSERT(item%value == initSize + 1 - pos)
      deallocate(item)
      @:ASSERT(list%size() == initSize - 1)
      call list%view(-pos + 1, item)
      @:ASSERT(item%value == initSize + 1 - pos + 1)
    #:endblock TEST


    #:block TEST("popBack")
      type(TInt), pointer :: item
      integer :: pos

      pos = initSize
      item => null()
      call list%pop(item, pos=pos)
      @:ASSERT(associated(item))
      @:ASSERT(item%value == initSize)
      deallocate(item)
      @:ASSERT(list%size() == initSize - 1)
      call list%view(pos - 1, item)
      @:ASSERT(item%value == pos - 1)
    #:endblock TEST


    #:block TEST("popBackNegative")
      type(TInt), pointer :: item
      integer :: pos

      pos = -1
      item => null()
      call list%pop(item, pos=pos)
      @:ASSERT(associated(item))
      @:ASSERT(item%value == initSize)
      deallocate(item)
      @:ASSERT(list%size() == initSize - 1)
      call list%view(-1, item)
      @:ASSERT(item%value == initSize - 1)
      call list%view(initSize - 1, item)
      @:ASSERT(item%value == initSize - 1)
    #:endblock TEST

  #:endblock TEST_FIXTURE

#:endblock TEST_SUITE


@:TEST_DRIVER()
