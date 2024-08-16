!--------------------------------------------------------------------------------------------------!
!  DFTB+: general package for performing fast atomistic simulations                                !
!  Copyright (C) 2006 - 2025  DFTB+ developers group                                               !
!                                                                                                  !
!  See the LICENSE file for terms of usage and distribution.                                       !
!--------------------------------------------------------------------------------------------------!

#:include "fortuno_serial.fypp"

module test_common_schedule
  use dftbp_common_schedule, only : getChunkRanges, getIndicesWithWorkload
  use fortuno_serial, only : suite => serial_suite_item, test_list
  $:FORTUNO_SERIAL_IMPORTS()
  implicit none

  private
  public :: tests

contains


  $:TEST("singleRank")
    integer :: localFirst, localLast

    call getChunkRanges(1, 0, 1, 42, localFirst, localLast)
    @:ASSERT(localFirst == 1)
    @:ASSERT(localLast == 42)
  $:END_TEST()

  $:TEST("singleRankOffset")
    integer :: localFirst, localLast

    call getChunkRanges(1, 0, 7, 42, localFirst, localLast)
    @:ASSERT(localFirst == 7)
    @:ASSERT(localLast == 42)
  $:END_TEST()


  $:TEST("multipleRanks")
    integer :: localFirst, localLast

    call getChunkRanges(3, 0, 1, 17, localFirst, localLast)
    @:ASSERT(localFirst == 1)
    @:ASSERT(localLast == 6)

    call getChunkRanges(3, 1, 1, 17, localFirst, localLast)
    @:ASSERT(localFirst == 7)
    @:ASSERT(localLast == 12)

    call getChunkRanges(3, 2, 1, 17, localFirst, localLast)
    @:ASSERT(localFirst == 13)
    @:ASSERT(localLast == 17)
  $:END_TEST()


  $:TEST("multipleRanksOffset")
    integer :: localFirst, localLast

    call getChunkRanges(3, 0, 3, 20, localFirst, localLast)
    @:ASSERT(localFirst == 3)
    @:ASSERT(localLast == 8)

    call getChunkRanges(3, 1, 3, 20, localFirst, localLast)
    @:ASSERT(localFirst == 9)
    @:ASSERT(localLast == 14)

    call getChunkRanges(3, 2, 3, 20, localFirst, localLast)
    @:ASSERT(localFirst == 15)
    @:ASSERT(localLast == 20)
  $:END_TEST()


  $:TEST("singleRankEqualWorkload")
    integer :: workload(10), ii
    integer, allocatable :: indices(:)

    workload(:) = 1
    call getIndicesWithWorkload(1, 0, 1, 10, workload, indices)
    @:ASSERT(size(indices) == 10)
    do ii = 1, 10
      @:ASSERT(indices(ii) == ii)
    end do

    workload(:) = 5
    call getIndicesWithWorkload(1, 0, 1, 10, workload, indices)
    @:ASSERT(size(indices) == 10)
    do ii = 1, 10
      @:ASSERT(indices(ii) == ii)
    end do
  $:END_TEST()


  $:TEST("singleRankZeroWorkload")
    integer :: workload(10), ii
    integer, allocatable :: indices(:)

    workload(1:5) = 0
    workload(6:10) = 1
    call getIndicesWithWorkload(1, 0, 1, 10, workload, indices)
    @:ASSERT(size(indices) == 10)
    do ii = 1, 10
      @:ASSERT(indices(ii) == ii)
    end do

    workload(:) = 0
    call getIndicesWithWorkload(1, 0, 1, 10, workload, indices)
    @:ASSERT(size(indices) == 10)
    do ii = 1, 10
      @:ASSERT(indices(ii) == ii)
    end do
  $:END_TEST()


  $:TEST("singleRankOffsetEqualWorkload")
    integer :: workload(10), ii
    integer, allocatable :: indices(:)

    workload(:) = -100
    workload(7:10) = 1
    call getIndicesWithWorkload(1, 0, 7, 10, workload, indices)
    @:ASSERT(size(indices) == 4)
    do ii = 1, 4
      @:ASSERT(indices(ii) == ii + 6)
    end do

    workload(7:10) = 5
    call getIndicesWithWorkload(1, 0, 7, 10, workload, indices)
    @:ASSERT(size(indices) == 4)
    do ii = 1, 4
      @:ASSERT(indices(ii) == ii + 6)
    end do
  $:END_TEST()


  $:TEST("twoRanksDifferentWorkload")
    integer :: workload(10)
    integer, allocatable :: indices(:)

    workload(1:2) = 1
    workload(3) = 5
    workload(4:5) = 0
    workload(6:8) = 1
    workload(9) = 10
    workload(10) = 2

    call getIndicesWithWorkload(2, 0, 1, 10, workload, indices)
    @:ASSERT(size(indices) == 3)
    @:ASSERT(all(indices == [1,3,9]))

    call getIndicesWithWorkload(2, 1, 1, 10, workload, indices)
    @:ASSERT(size(indices) == 7)
    @:ASSERT(all(indices == [2,4,5,6,7,8,10]))
  $:END_TEST()


  $:TEST("threeRanksDifferentWorkload")
    integer :: workload(10)
    integer, allocatable :: indices(:)

    workload(1) = 10
    workload(2) = 2
    workload(3) = 3
    workload(4) = 0
    workload(5) = 10
    workload(6) = 1
    workload(7) = 1
    workload(8) = 1
    workload(9) = 10
    workload(10) = 1

    call getIndicesWithWorkload(3, 0, 1, 10, workload, indices)
    @:ASSERT(size(indices) == 2)
    @:ASSERT(all(indices == [1,10]))

    call getIndicesWithWorkload(3, 1, 1, 10, workload, indices)
    @:ASSERT(size(indices) == 3)
    @:ASSERT(all(indices == [2,4,5]))

    call getIndicesWithWorkload(3, 2, 1, 10, workload, indices)
    @:ASSERT(size(indices) == 5)
    @:ASSERT(all(indices == [3,6,7,8,9]))
  $:END_TEST()


  $:TEST("threeRanksUnbalancedWorkload")
    integer :: workload(10)
    integer, allocatable :: indices(:)

    workload(:) = 1
    workload(1) = 100
    call getIndicesWithWorkload(3, 0, 1, 10, workload, indices)
    @:ASSERT(size(indices) == 1)
    @:ASSERT(all(indices == [1]))
    call getIndicesWithWorkload(3, 1, 1, 10, workload, indices)
    @:ASSERT(size(indices) == 5)
    @:ASSERT(all(indices == [2,4,6,8,10]))
    call getIndicesWithWorkload(3, 2, 1, 10, workload, indices)
    @:ASSERT(size(indices) == 4)
    @:ASSERT(all(indices == [3,5,7,9]))

    workload(:) = 1
    workload(10) = 100
    call getIndicesWithWorkload(3, 0, 1, 10, workload, indices)
    @:ASSERT(size(indices) == 4)
    @:ASSERT(all(indices == [1,4,7,10]))
    call getIndicesWithWorkload(3, 1, 1, 10, workload, indices)
    @:ASSERT(size(indices) == 3)
    @:ASSERT(all(indices == [2,5,8]))
    call getIndicesWithWorkload(3, 2, 1, 10, workload, indices)
    @:ASSERT(size(indices) == 3)
    @:ASSERT(all(indices == [3,6,9]))
  $:END_TEST()


  function tests()
    type(test_list) :: tests

    tests = test_list([&
        suite("schedule", test_list([&
            $:TEST_ITEMS()
        ]))&
    ])
    $:STOP_ON_MISSING_TEST_ITEMS()

  end function tests

end module test_common_schedule
