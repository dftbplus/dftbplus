!--------------------------------------------------------------------------------------------------!
!  DFTB+: general package for performing fast atomistic simulations                                !
!  Copyright (C) 2006 - 2025  DFTB+ developers group                                               !
!                                                                                                  !
!  See the LICENSE file for terms of usage and distribution.                                       !
!--------------------------------------------------------------------------------------------------!

#:include "fortuno_serial.fypp"

module test_include_allocatablelist
  use alloclisthelper, only : TIntList
  use fortuno_serial, only : serial_case_base, suite => serial_suite_item, test_item, test_list
  $:FORTUNO_SERIAL_IMPORTS()
  implicit none

  private
  public :: tests


  ! Fixtured test providing an initialized small list at startup
  type, extends(serial_case_base) :: small_list_case
      procedure(small_list_case_proc), pointer, nopass :: proc
  contains
    procedure :: run => small_list_case_run
  end type small_list_case

  abstract interface
    subroutine small_list_case_proc(list)
      import :: TIntList
      type(TIntList), target, intent(inout) :: list
    end subroutine small_list_case_proc
  end interface


  ! Initial size of the small list
  integer, parameter :: initSize = 10

contains


  $:TEST("uninitialized")
    type(TIntList) :: list

    @:ASSERT(list%size() == 0)
  $:END_TEST()


  $:TEST("pushAndPopItems")
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

  $:END_TEST()


  $:TEST("pushFront", args=["list"], label="smallList")
    type(TIntList), target, intent(inout) :: list

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
  $:END_TEST()


  $:TEST("pushFrontNegative", args=["list"], label="smallList")
    type(TIntList), target, intent(inout) :: list

    integer, allocatable :: item
    integer, pointer :: pItem

    item = 0
    call list%push(item, pos=-initSize)
    @:ASSERT(list%size() == initSize + 1)
    call list%view(2, pItem)
    @:ASSERT(pItem == 1)
    call list%view(1, pItem)
    @:ASSERT(pItem == 0)
  $:END_TEST()


  $:TEST("pushMiddle", args=["list"], label="smallList")
    type(TIntList), target, intent(inout) :: list

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
  $:END_TEST()


  $:TEST("pushMiddleNegative", args=["list"], label="smallList")
    type(TIntList), target, intent(inout) :: list

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
  $:END_TEST()


  $:TEST("pushBack", args=["list"], label="smallList")
    type(TIntList), target, intent(inout) :: list

    integer, allocatable :: item
    integer, pointer :: pItem
    integer :: pos

    pos = initSize + 1
    item = 0
    call list%push(item, pos=pos)
    @:ASSERT(list%size() == initSize + 1)
    call list%view(pos, pItem)
    @:ASSERT(pItem == 0)
  $:END_TEST()


  $:TEST("pushBackZero", args=["list"], label="smallList")
    type(TIntList), target, intent(inout) :: list

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
  $:END_TEST()


  $:TEST("popFront", args=["list"], label="smallList")
    type(TIntList), target, intent(inout) :: list

    integer, allocatable :: item
    integer, pointer :: pItem

    call list%pop(item, pos=1)
    @:ASSERT(allocated(item))
    @:ASSERT(item == 1)
    @:ASSERT(list%size() == initSize - 1)
    call list%view(1, pItem)
    @:ASSERT(pItem == 2)
  $:END_TEST()


  $:TEST("popFrontNegative", args=["list"], label="smallList")
    type(TIntList), target, intent(inout) :: list

    integer, allocatable :: item
    integer, pointer :: pItem

    call list%pop(item, pos=-initSize)
    @:ASSERT(allocated(item))
    @:ASSERT(item == 1)
    @:ASSERT(list%size() == initSize - 1)
    call list%view(1, pItem)
    @:ASSERT(pItem == 2)
  $:END_TEST()


  $:TEST("popMiddle", args=["list"], label="smallList")
    type(TIntList), target, intent(inout) :: list

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
  $:END_TEST()


  $:TEST("popMiddleNegative", args=["list"], label="smallList")
    type(TIntList), target, intent(inout) :: list

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
  $:END_TEST()


  $:TEST("popBack", args=["list"], label="smallList")
    type(TIntList), target, intent(inout) :: list

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
  $:END_TEST()


  $:TEST("popBackNegative", args=["list"], label="smallList")
    type(TIntList), target, intent(inout) :: list

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
  $:END_TEST()

  ! Returns the tests from this module.
  function tests()
    type(test_list) :: tests

    tests = test_list([&
        suite("allocatablelist", test_list([&
          $:TEST_ITEMS(label="", suffix=",")
          $:TEST_ITEMS(constructor="small_list_test", label="smallList")
        ]))&
    ])
    $:STOP_ON_MISSING_TEST_ITEMS()

  end function tests


  ! Convenience function returning a small_list_case instance wrapped as test_item.
  function small_list_test(name, proc) result(testitem)
    character(*), intent(in) :: name
    procedure(small_list_case_proc) :: proc
    type(test_item) :: testitem

    testitem = test_item(small_list_case(name=name, proc=proc))

  end function small_list_test


   ! Run procedure of the small_list_case class.
  subroutine small_list_case_run(this)
    class(small_list_case), intent(in) :: this

    type(TIntList), target :: list

    call initSmallList_(list)
    call this%proc(list)
    ! Workaround:gfortran:14.1 (bug 116679)
    ! GFortran fails to finalize allocatable components of derived types within a pointer array
    ! {+
    block
      call finalSmallList(list)
    end block
    ! +}

  end subroutine small_list_case_run


  subroutine initSmallList_(list)
    type(TIntList), intent(out) :: list

    integer, allocatable :: item
    integer :: ii

    do ii = 1, initSize
      item = ii
      call list%push(item)
    end do

  end subroutine initSmallList_


  ! Workaround:gfortran:14.1 (bug 116679)
  ! GFortran fails to finalize allocatable components of derived types within a pointer array
  ! {+
  subroutine finalSmallList(list)
    type(TIntList), intent(inout) :: list

    integer, allocatable :: item
    integer :: ii

    do ii = 1, list%size()
      call list%pop(item)
      deallocate(item)
    end do

  end subroutine finalSmallList
  ! +}

end module test_include_allocatablelist
