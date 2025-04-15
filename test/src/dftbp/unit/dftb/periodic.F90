!--------------------------------------------------------------------------------------------------!
!  DFTB+: general package for performing fast atomistic simulations                                !
!  Copyright (C) 2006 - 2025  DFTB+ developers group                                               !
!                                                                                                  !
!  See the LICENSE file for terms of usage and distribution.                                       !
!--------------------------------------------------------------------------------------------------!

#:include "fortuno_serial.fypp"

module test_dftb_periodic
  use dftbp_common_accuracy, only : dp, minNeighDist
  use dftbp_common_environment, only : TEnvironment, TEnvironment_init
  use dftbp_common_globalenv, only : destructGlobalEnv, initGlobalEnv, stdOut
  use dftbp_dftb_periodic, only : distributeAtoms, reallocateArrays2, allocateNeighbourArrays,&
      & fillNeighbourArrays, TNeighbourList, updateNeighbourList, TNeighbourlist_init
  use dftbp_common_status, only : TStatus
  use fortuno_serial, only : suite => serial_suite_item, test_list
  $:FORTUNO_SERIAL_IMPORTS()
  implicit none

  private
  public :: tests

  type :: test_env
    type(TEnvironment) :: env
  contains
    final :: final_test_env
  end type test_env

contains


  !> Intializes the test environment
  subroutine init_test_env(this)
    type(test_env), intent(out) :: this

    call initGlobalEnv()
    call TEnvironment_init(this%env)
    ! temporary fix
    this%env%stdOut = stdOut

  end subroutine init_test_env


  !> Finalizes the test environment
  subroutine final_test_env(this)
    type(test_env), intent(inout) :: this

    call this%env%destruct()
    call destructGlobalEnv()

  end subroutine final_test_env


  $:TEST("singleRank")
    type(test_env) :: tenv
    integer :: startAtom, endAtom
    logical :: error

    call init_test_env(tenv)

    call distributeAtoms(tenv%env%stdOut, 0, 1, 42, startAtom, endAtom, error)
    @:ASSERT(startAtom == 1)
    @:ASSERT(endAtom == 42)
    @:ASSERT(.not. error)
  $:END_TEST()


  $:TEST("multipleRanks")
    type(test_env) :: tenv
    integer :: startAtom, endAtom
    logical :: error

    call init_test_env(tenv)

    call distributeAtoms(tenv%env%stdOut, 0, 2, 13, startAtom, endAtom, error)
    @:ASSERT(startAtom == 1)
    @:ASSERT(endAtom == 7)
    @:ASSERT(.not. error)
    call distributeAtoms(tenv%env%stdOut, 1, 2, 13, startAtom, endAtom, error)
    @:ASSERT(startAtom == 8)
    @:ASSERT(endAtom == 13)
    @:ASSERT(.not. error)
  $:END_TEST()


  $:TEST("tooManyRanks")
    type(test_env) :: tenv
    integer :: startAtom, endAtom, ii
    logical :: error

    call init_test_env(tenv)

    do ii = 1, 4
      call distributeAtoms(tenv%env%stdOut, ii, 4, 2, startAtom, endAtom, error)
      @:ASSERT(error)
    end do
  $:END_TEST()


  $:TEST("increaseSizeR2")
    integer, allocatable :: iNeighbour(:,:)
    real(dp), allocatable :: neighDist2(:,:), neighDist2Test(:,:)
    integer :: maxNeighbour, startAtom, endAtom

    maxNeighbour = 3
    startAtom = 17
    endAtom = 27
    allocate(iNeighbour(1:maxNeighbour, startAtom:endAtom))
    allocate(neighDist2(1:maxNeighbour, startAtom:endAtom))
    allocate(neighDist2Test(1:maxNeighbour, startAtom:endAtom))

    call createTestArray_(neighDist2Test)
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
  $:END_TEST()


  $:TEST("decreaseSizeR2")
    integer, allocatable :: iNeighbour(:,:)
    real(dp), allocatable :: neighDist2(:,:), neighDist2Test(:,:)
    integer :: maxNeighbour, startAtom, endAtom

    maxNeighbour = 3
    startAtom = 17
    endAtom = 27
    allocate(iNeighbour(1:maxNeighbour, startAtom:endAtom))
    allocate(neighDist2(1:maxNeighbour, startAtom:endAtom))
    allocate(neighDist2Test(1:maxNeighbour, startAtom:endAtom))

    call createTestArray_(neighDist2Test)
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

    @:ASSERT(all(neighDist2(:maxNeighbour,startAtom:endAtom) ==&
        & neighDist2Test(:maxNeighbour,startAtom:endAtom)))
    deallocate(iNeighbour, neighDist2, neighDist2Test)
  $:END_TEST()


  $:TEST("fillNeighbourArrays")
    type(TNeighbourList) :: neigh
    integer, allocatable :: iNeighbour(:,:)
    real(dp), allocatable :: neighDist2(:,:)
    integer :: maxNeighbour, nAtom, startAtom, endAtom, ii

    maxNeighbour = 3
    nAtom = 4
    allocate(iNeighbour(1:maxNeighbour, 1:nAtom))
    allocate(neighDist2(1:maxNeighbour, 1:nAtom))

    call createTestArray_(neighDist2)

    call allocateNeighbourArrays(neigh, maxNeighbour, nAtom)

    call fillNeighbourArrays(neigh, iNeighbour, neighDist2, startAtom, endAtom, maxNeighbour,&
        & nAtom)

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
    do ii = 1, nAtom
      @:ASSERT(neigh%iNeighbour(0, ii) == ii)
    end do
    deallocate(iNeighbour, neighDist2)
  $:END_TEST()


  $:TEST("neighbourdistance")
    real(dp), allocatable :: coords(:,:), coords0(:,:), rCellVec(:,:)
    integer, allocatable :: img2CentCell(:), iCellVec(:)
    type(TNeighbourList) :: neigh
    real(dp) :: mCutoff
    type(TStatus) :: errStatus
    integer :: nAtom
    type(TNeighbourList) :: neighs

    @:ASSERT(.true.)
    nAtom = 2
    allocate(coords(3, nAtom))
    allocate(img2CentCell(nAtom))
    allocate(iCellVec(nAtom))
    allocate(coords0(3, nAtom))
    allocate(rCellVec(3, 1))
    rCellVec(:, :) = 0.0_dp
    iCellVec(:) = 0
    mCutoff = 1.0_dp
    coords0(:, :) = 0.0_dp
    call TNeighbourlist_init(neighs, nAtom, 2)
    coords0(1, 2) = minNeighDist
    call updateNeighbourList(coords, img2CentCell, iCellVec, neighs, nAtom, coords0, mCutoff,&
        & rCellVec, errStatus)
    @:ASSERT(.not.errStatus%hasError())
    ! tolerance is small, so this value is now smaller and corresponds to a previous code
    ! regression (but _should_ trigger an error):
    coords0(1, 2) = minNeighDist**2
    call updateNeighbourList(coords, img2CentCell, iCellVec, neighs, nAtom, coords0, mCutoff,&
        & rCellVec, errStatus)
    @:ASSERT(errStatus%hasError())

  $:END_TEST()


  function tests()
    type(test_list) :: tests

    tests = test_list([&
        suite("periodic", test_list([&
            $:TEST_ITEMS()
        ]))&
    ])
    $:STOP_ON_MISSING_TEST_ITEMS()

  end function tests


  subroutine createTestArray_(input)
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
  end subroutine createTestArray_

end module test_dftb_periodic
