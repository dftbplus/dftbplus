!--------------------------------------------------------------------------------------------------!
!  DFTB+: general package for performing fast atomistic simulations                                !
!  Copyright (C) 2006 - 2023  DFTB+ developers group                                               !
!                                                                                                  !
!  See the LICENSE file for terms of usage and distribution.                                       !
!--------------------------------------------------------------------------------------------------!

#:include "fytest.fypp"

#:block TEST_SUITE("schedule")
  use dftbp_common_schedule, only : getChunkRanges, getChunkIterWithWorkload, TChunkIterator
  implicit none

#:contains

  #:block TEST_FIXTURE("getChunkRanges")

    integer :: localFirst, localLast

  #:contains

    #:block TEST("singleRank")
      call getChunkRanges(1, 0, 1, 42, localFirst, localLast)
      @:ASSERT(localFirst == 1)
      @:ASSERT(localLast == 42)
    #:endblock

    #:block TEST("singleRankOffset")
      call getChunkRanges(1, 0, 7, 42, localFirst, localLast)
      @:ASSERT(localFirst == 7)
      @:ASSERT(localLast == 42)
    #:endblock

    #:block TEST("multipleRanks")
      call getChunkRanges(3, 0, 1, 17, localFirst, localLast)
      @:ASSERT(localFirst == 1)
      @:ASSERT(localLast == 6)

      call getChunkRanges(3, 1, 1, 17, localFirst, localLast)
      @:ASSERT(localFirst == 7)
      @:ASSERT(localLast == 12)

      call getChunkRanges(3, 2, 1, 17, localFirst, localLast)
      @:ASSERT(localFirst == 13)
      @:ASSERT(localLast == 17)
    #:endblock

    #:block TEST("multipleRanksOffset")
      call getChunkRanges(3, 0, 3, 20, localFirst, localLast)
      @:ASSERT(localFirst == 3)
      @:ASSERT(localLast == 8)

      call getChunkRanges(3, 1, 3, 20, localFirst, localLast)
      @:ASSERT(localFirst == 9)
      @:ASSERT(localLast == 14)

      call getChunkRanges(3, 2, 3, 20, localFirst, localLast)
      @:ASSERT(localFirst == 15)
      @:ASSERT(localLast == 20)
    #:endblock

  #:endblock TEST_FIXTURE

  #:block TEST_FIXTURE("getChunkIterWithWorkload")

    integer :: workload(10), i
    type(TChunkIterator) :: iter

  #:contains

    #:block TEST("singleRankEqualWorkload")
      workload(:) = 1
      call getChunkIterWithWorkload(1, 0, 1, 10, workload, iter)
      @:ASSERT(iter%getNumIndices() == 10)
      do i = 1, 10
        @:ASSERT(iter%getNextIndex() == i)
      end do
      @:ASSERT(.not. iter%hasNextIndex())

      workload(:) = 5
      call getChunkIterWithWorkload(1, 0, 1, 10, workload, iter)
      @:ASSERT(iter%getNumIndices() == 10)
      do i = 1, 10
        @:ASSERT(iter%getNextIndex() == i)
      end do
      @:ASSERT(.not. iter%hasNextIndex())
    #:endblock

    #:block TEST("singleRankZeroWorkload")
      workload(1:5) = 0
      workload(6:10) = 1
      call getChunkIterWithWorkload(1, 0, 1, 10, workload, iter)
      @:ASSERT(iter%getNumIndices() == 10)
      do i = 1, 10
        @:ASSERT(iter%getNextIndex() == i)
      end do
      @:ASSERT(.not. iter%hasNextIndex())

      workload(:) = 0
      call getChunkIterWithWorkload(1, 0, 1, 10, workload, iter)
      @:ASSERT(iter%getNumIndices() == 10)
      do i = 1, 10
        @:ASSERT(iter%getNextIndex() == i)
      end do
      @:ASSERT(.not. iter%hasNextIndex())
    #:endblock

    #:block TEST("singleRankOffsetEqualWorkload")
      workload(:) = -100
      workload(7:10) = 1
      call getChunkIterWithWorkload(1, 0, 7, 10, workload, iter)
      @:ASSERT(iter%getNumIndices() == 4)
      do i = 7, 10
        @:ASSERT(iter%getNextIndex() == i)
      end do
      @:ASSERT(.not. iter%hasNextIndex())

      workload(7:10) = 5
      call getChunkIterWithWorkload(1, 0, 7, 10, workload, iter)
      @:ASSERT(iter%getNumIndices() == 4)
      do i = 7, 10
        @:ASSERT(iter%getNextIndex() == i)
      end do
      @:ASSERT(.not. iter%hasNextIndex())
    #:endblock

    #:block TEST("twoRanksDifferentWorkload")
      workload(1:2) = 1
      workload(3) = 5
      workload(4:5) = 0
      workload(6:8) = 1
      workload(9) = 10
      workload(10) = 2

      call getChunkIterWithWorkload(2, 0, 1, 10, workload, iter)
      @:ASSERT(iter%getNumIndices() == 3)
      @:ASSERT(iter%getNextIndex() == 1)
      @:ASSERT(iter%getNextIndex() == 3)
      @:ASSERT(iter%getNextIndex() == 9)
      @:ASSERT(.not. iter%hasNextIndex())

      call getChunkIterWithWorkload(2, 1, 1, 10, workload, iter)
      @:ASSERT(iter%getNumIndices() == 7)
      @:ASSERT(iter%getNextIndex() == 2)
      @:ASSERT(iter%getNextIndex() == 4)
      @:ASSERT(iter%getNextIndex() == 5)
      @:ASSERT(iter%getNextIndex() == 6)
      @:ASSERT(iter%getNextIndex() == 7)
      @:ASSERT(iter%getNextIndex() == 8)
      @:ASSERT(iter%getNextIndex() == 10)
      @:ASSERT(.not. iter%hasNextIndex())
    #:endblock

    #:block TEST("threeRanksDifferentWorkload")
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

      call getChunkIterWithWorkload(3, 0, 1, 10, workload, iter)
      @:ASSERT(iter%getNumIndices() == 2)
      @:ASSERT(iter%getNextIndex() == 1)
      @:ASSERT(iter%getNextIndex() == 10)
      @:ASSERT(.not. iter%hasNextIndex())

      call getChunkIterWithWorkload(3, 1, 1, 10, workload, iter)
      @:ASSERT(iter%getNumIndices() == 3)
      @:ASSERT(iter%getNextIndex() == 2)
      @:ASSERT(iter%getNextIndex() == 4)
      @:ASSERT(iter%getNextIndex() == 5)
      @:ASSERT(.not. iter%hasNextIndex())

      call getChunkIterWithWorkload(3, 2, 1, 10, workload, iter)
      @:ASSERT(iter%getNumIndices() == 5)
      @:ASSERT(iter%getNextIndex() == 3)
      @:ASSERT(iter%getNextIndex() == 6)
      @:ASSERT(iter%getNextIndex() == 7)
      @:ASSERT(iter%getNextIndex() == 8)
      @:ASSERT(iter%getNextIndex() == 9)
      @:ASSERT(.not. iter%hasNextIndex())
    #:endblock

    #:block TEST("threeRanksUnbalancedWorkload")
      workload(:) = 1
      workload(1) = 100
      call getChunkIterWithWorkload(3, 0, 1, 10, workload, iter)
      @:ASSERT(iter%getNumIndices() == 1)
      @:ASSERT(iter%getNextIndex() == 1)
      @:ASSERT(.not. iter%hasNextIndex())
      call getChunkIterWithWorkload(3, 1, 1, 10, workload, iter)
      @:ASSERT(iter%getNumIndices() == 5)
      @:ASSERT(iter%getNextIndex() == 2)
      @:ASSERT(iter%getNextIndex() == 4)
      @:ASSERT(iter%getNextIndex() == 6)
      @:ASSERT(iter%getNextIndex() == 8)
      @:ASSERT(iter%getNextIndex() == 10)
      @:ASSERT(.not. iter%hasNextIndex())
      call getChunkIterWithWorkload(3, 2, 1, 10, workload, iter)
      @:ASSERT(iter%getNumIndices() == 4)
      @:ASSERT(iter%getNextIndex() == 3)
      @:ASSERT(iter%getNextIndex() == 5)
      @:ASSERT(iter%getNextIndex() == 7)
      @:ASSERT(iter%getNextIndex() == 9)
      @:ASSERT(.not. iter%hasNextIndex())

      workload(:) = 1
      workload(10) = 100
      call getChunkIterWithWorkload(3, 0, 1, 10, workload, iter)
      @:ASSERT(iter%getNumIndices() == 4)
      @:ASSERT(iter%getNextIndex() == 1)
      @:ASSERT(iter%getNextIndex() == 4)
      @:ASSERT(iter%getNextIndex() == 7)
      @:ASSERT(iter%getNextIndex() == 10)
      @:ASSERT(.not. iter%hasNextIndex())
      call getChunkIterWithWorkload(3, 1, 1, 10, workload, iter)
      @:ASSERT(iter%getNumIndices() == 3)
      @:ASSERT(iter%getNextIndex() == 2)
      @:ASSERT(iter%getNextIndex() == 5)
      @:ASSERT(iter%getNextIndex() == 8)
      @:ASSERT(.not. iter%hasNextIndex())
      call getChunkIterWithWorkload(3, 2, 1, 10, workload, iter)
      @:ASSERT(iter%getNumIndices() == 3)
      @:ASSERT(iter%getNextIndex() == 3)
      @:ASSERT(iter%getNextIndex() == 6)
      @:ASSERT(iter%getNextIndex() == 9)
      @:ASSERT(.not. iter%hasNextIndex())
    #:endblock

  #:endblock TEST_FIXTURE

#:endblock TEST_SUITE


@:TEST_DRIVER()
