!--------------------------------------------------------------------------------------------------!
!  DFTB+: general package for performing fast atomistic simulations                                !
!  Copyright (C) 2006 - 2022  DFTB+ developers group                                               !
!                                                                                                  !
!  See the LICENSE file for terms of usage and distribution.                                       !
!--------------------------------------------------------------------------------------------------!

#:include "fytest.fypp"

#:block TEST_SUITE("period")
  use dftbp_common_accuracy, only : dp
  use dftbp_dftb_periodic, only : distributeAtoms, reallocateArrays2, fillNeighbourArrays,&
      & TNeighbourList
  use dftbp_common_status, only : TStatus
  implicit none

#:contains

  #:block TEST_FIXTURE("distributeAtoms")

    integer :: startAtom, endAtom
    logical :: error

  #:contains

    #:block TEST("singleRank")
      call distributeAtoms(0, 1, 42, startAtom, endAtom, error)
      @:ASSERT(startAtom == 1)
      @:ASSERT(endAtom == 42)
      @:ASSERT(.not. error)
    #:endblock

    #:block TEST("multipleRanks")
      call distributeAtoms(0, 2, 13, startAtom, endAtom, error)
      @:ASSERT(startAtom == 1)
      @:ASSERT(endAtom == 7)
      @:ASSERT(.not. error)

      call distributeAtoms(1, 2, 13, startAtom, endAtom, error)
      @:ASSERT(startAtom == 8)
      @:ASSERT(endAtom == 13)
      @:ASSERT(.not. error)
    #:endblock

    #:block TEST("tooManyRanks")
      integer :: i
      do i = 1, 4
        call distributeAtoms(i, 4, 2, startAtom, endAtom, error)
        @:ASSERT(error)
      end do
    #:endblock

  #:endblock TEST_FIXTURE

  #:block TEST_FIXTURE("reallocateArrays2")

    integer, allocatable :: iNeighbour(:,:)
    real(dp), allocatable :: neighDist2(:,:), neighDist2Test(:,:)
    integer :: maxNeighbour, startAtom, endAtom

  #:contains

    #:block TEST("increaseSize")
      maxNeighbour = 3
      startAtom = 17
      endAtom = 27
      allocate(iNeighbour(1:maxNeighbour, startAtom:endAtom))
      allocate(neighDist2(1:maxNeighbour, startAtom:endAtom))
      allocate(neighDist2Test(1:maxNeighbour, startAtom:endAtom))

      call createTestArray(neighDist2Test)
      neighDist2(:,:) = neighDist2Test(:,:)

      maxNeighbour = 4
      call reallocateArrays2(iNeighbour, neighDist2, maxNeighbour)

      @:ASSERT(lbound(iNeighbour, dim=1) == 1)
      @:ASSERT(ubound(iNeighbour, dim=1) == maxNeighbour)
      @:ASSERT(lbound(neighDist2, dim=1) == 1)
      @:ASSERT(ubound(neighDist2, dim=1) == maxNeighbour)

      @:ASSERT(lbound(iNeighbour, dim=2) == startAtom)
      @:ASSERT(ubound(iNeighbour, dim=2) == endAtom)
      @:ASSERT(lbound(neighDist2, dim=2) == startAtom)
      @:ASSERT(ubound(neighDist2, dim=2) == endAtom)

      @:ASSERT(all(neighDist2(:3,startAtom:endAtom) == neighDist2Test(:3,startAtom:endAtom)))
      @:ASSERT(all(neighDist2(maxNeighbour:,startAtom:endAtom) == 0.0_dp))

      deallocate(iNeighbour, neighDist2, neighDist2Test)
    #:endblock

    #:block TEST("decreaseSize")
      maxNeighbour = 3
      startAtom = 17
      endAtom = 27
      allocate(iNeighbour(1:maxNeighbour, startAtom:endAtom))
      allocate(neighDist2(1:maxNeighbour, startAtom:endAtom))
      allocate(neighDist2Test(1:maxNeighbour, startAtom:endAtom))

      call createTestArray(neighDist2Test)
      neighDist2(:,:) = neighDist2Test(:,:)

      maxNeighbour = 2
      call reallocateArrays2(iNeighbour, neighDist2, maxNeighbour)

      @:ASSERT(lbound(iNeighbour, dim=1) == 1)
      @:ASSERT(ubound(iNeighbour, dim=1) == maxNeighbour)
      @:ASSERT(lbound(neighDist2, dim=1) == 1)
      @:ASSERT(ubound(neighDist2, dim=1) == maxNeighbour)

      @:ASSERT(lbound(iNeighbour, dim=2) == startAtom)
      @:ASSERT(ubound(iNeighbour, dim=2) == endAtom)
      @:ASSERT(lbound(neighDist2, dim=2) == startAtom)
      @:ASSERT(ubound(neighDist2, dim=2) == endAtom)

      @:ASSERT(all(neighDist2(:maxNeighbour,startAtom:endAtom) == neighDist2Test(:maxNeighbour,startAtom:endAtom)))

      deallocate(iNeighbour, neighDist2, neighDist2Test)
    #:endblock

  #:endblock TEST_FIXTURE

  #:block TEST_FIXTURE("fillNeighbourArrays")

    type(TNeighbourList) :: neigh
    integer, allocatable :: iNeighbour(:,:)
    real(dp), allocatable :: neighDist2(:,:)
    integer :: maxNeighbour, nAtom, startAtom, endAtom

  #:contains

    #:block TEST("notParallel")
      integer :: i

      maxNeighbour = 3
      nAtom = 4
      allocate(iNeighbour(1:maxNeighbour, 1:nAtom))
      allocate(neighDist2(1:maxNeighbour, 1:nAtom))

      call createTestArray(neighDist2)

      call fillNeighbourArrays(neigh, iNeighbour, neighDist2, startAtom, endAtom, maxNeighbour, nAtom, .false.)

      @:ASSERT(lbound(neigh%iNeighbour, dim=1) == 0)
      @:ASSERT(ubound(neigh%iNeighbour, dim=1) == maxNeighbour)
      @:ASSERT(lbound(neigh%neighDist2, dim=1) == 0)
      @:ASSERT(ubound(neigh%neighDist2, dim=1) == maxNeighbour)

      @:ASSERT(lbound(neigh%iNeighbour, dim=2) == 1)
      @:ASSERT(ubound(neigh%iNeighbour, dim=2) == nAtom)
      @:ASSERT(lbound(neigh%neighDist2, dim=2) == 1)
      @:ASSERT(ubound(neigh%neighDist2, dim=2) == nAtom)

      @:ASSERT(all(neigh%neighDist2(1:maxNeighbour,1:nAtom) == neighDist2(1:maxNeighbour,1:nAtom)))
      @:ASSERT(all(neigh%neighDist2(0,1:nAtom) == 0.0_dp))
      do i = 1, nAtom
        @:ASSERT(neigh%iNeighbour(0,i) == i)
      end do

      deallocate(iNeighbour, neighDist2)
    #:endblock

  #:endblock TEST_FIXTURE

  subroutine createTestArray(input)
    real(dp), intent(inout) :: input(:,:)
    real(dp) :: val
    integer :: i, j

    val = 0.0_dp
    do i = 1, size(input, dim=1)
      do j = lbound(input, dim=2), ubound(input, dim=2)
        input(i,j) = val
        val = val + 1.0_dp
      end do
    end do
  end subroutine

#:endblock TEST_SUITE


@:TEST_DRIVER()
