!--------------------------------------------------------------------------------------------------!
!  DFTB+: general package for performing fast atomistic simulations                                !
!  Copyright (C) 2006 - 2022  DFTB+ developers group                                               !
!                                                                                                  !
!  See the LICENSE file for terms of usage and distribution.                                       !
!--------------------------------------------------------------------------------------------------!

#:include "fytest.fypp"

#:block TEST_SUITE("schedule")
  use dftbp_common_schedule, only : getChunkRanges
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

#:endblock TEST_SUITE


@:TEST_DRIVER()
