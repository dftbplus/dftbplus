!--------------------------------------------------------------------------------------------------!
!  DFTB+: general package for performing fast atomistic simulations                                !
!  Copyright (C) 2017  DFTB+ developers group                                                      !
!                                                                                                  !
!  See the LICENSE file for terms of usage and distribution.                                       !
!--------------------------------------------------------------------------------------------------!

#:include 'common.fypp'

!> Contains BLACS environmental settings.
module blacsenv
  use mpienv
  use message
  use scalapackfx
  use assert
  implicit none
  private

  public :: TBlacsEnv, TBlacsEnv_init


  !> Contains various BLACS related settings
  type :: TBlacsEnv

    !> Group grid for (nOrb, nOrb) shaped square matrices (nAtom) shaped vectors
    !> Note: the grid always contains all processes in the group
    type(blacsgrid) :: orbitalGrid

    !> Group grid for (nAtom, nAtom) shaped square matrices (nAtom) shaped vectors
    !> Note: Some processes in the group may be outside of this grid!
    type(blacsgrid) :: atomGrid

    !> Row block size
    integer :: rowBlockSize

    !> Column block size
    integer :: columnBlockSize

  end type TBlacsEnv


contains


  !> Initializes BLACS grids
  subroutine TBlacsEnv_init(this, myMpiEnv, rowBlock, colBlock, nOrb, nAtom)

    !> Initialized instance at exit.
    type(TBlacsEnv), intent(out) :: this

    !> Initialised MPI environment
    type(TMpiEnv), intent(in) :: myMpiEnv

    !> Row block size
    integer, intent(in) :: rowBlock

    !> Column block size
    integer, intent(in) :: colBlock

    !> Nr. of orbitals
    integer, intent(in) :: nOrb

    !> Nr. of atoms
    integer, intent(in) :: nAtom

    integer, allocatable :: gridMap(:,:)
    integer :: nProcRow, nProcCol, maxProcRow, maxProcColMax
    character(200) :: buffer

    ! Create orbital grid for each processor group
    call getSquareGridParams(myMpiEnv%groupSize, nProcRow, nProcCol)
    ! Check whether all processes have some portions of the H and S matrices otherwise
    ! diagonalisers may return garbage
    maxProcRow = (nOrb - 1) / rowBlock + 1
    maxProcColMax = (nOrb - 1) / colblock + 1
    if (nProcRow > maxProcRow .or. nProcCol > maxProcColMax) then
      write(buffer, "(A,I0,A,I0,A,I0,A,I0,A)") "Processor grid (", nProcRow, " x ",  nProcCol,&
          & ") too big (> ", maxProcRow, " x ", maxProcColMax, ")"
      call error(buffer)
    end if
    call getGridMap(myMpiEnv%groupMembers, nProcRow, nProcCol, gridMap)
    call this%orbitalGrid%initmappedgrids(gridMap)

    ! Create atom grid for each processor group
    maxProcRow = (nAtom - 1) / rowBlock + 1
    maxProcColMax = (nAtom - 1) / colBlock + 1
    nProcRow = min(nProcRow, maxProcRow)
    nProcCol = min(nProcCol, maxProcColMax)
    call getGridMap(myMpiEnv%groupMembers, nProcRow, nProcCol, gridMap)
    call this%atomGrid%initmappedgrids(gridMap)

    this%rowBlockSize = rowBlock
    this%columnBlockSize = colBlock

  end subroutine TBlacsEnv_init


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!! Private routines
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Returns rows and columns for a 2D processor grid closest possible to a square one.
  subroutine getSquareGridParams(nProc, nRow, nCol)

    !> Total number of processors
    integer, intent(in) :: nProc

    !> Number of rows in the 2D processor grid
    integer, intent(out) :: nRow

    !> Number of columns in the 2D processor grid
    integer, intent(out) :: nCol

    do nRow = int(sqrt(real(nProc))), 1, -1
      if (mod(nProc, nRow) == 0) then
        exit
      end if
    end do
    nCol = nProc / nRow

  end subroutine getSquareGridParams


  !> Returns the gridmap which can be used to initialise a 2D BLACS grid.
  !>
  subroutine getGridMap(groupMembers, nProcRow, nProcCol, gridMap)

    !> All members of the current group (must be >= nProcRow * nProcCol) which should
    !> be mapped on a 2D BLACS-grid.
    integer, intent(in) :: groupMembers(:)

    !> Number of process rows in the BLACS-grid
    integer, intent(in) :: nProcRow

    !> Number of process columns in the BLACS-grid
    integer, intent(in) :: nProcCol

    !> Grid map, where gridMap(i,j) contains the id of the process which in the ith row and jth
    !> column in the 2D BLACS-grid. This can be used as argument to the initmappedgrids()
    !> method of blacsgrid.
    integer, allocatable, intent(out) :: gridMap(:,:)

    integer :: iProcRow, iProcCol, ind

    @:ASSERT(size(groupMembers) >= nProcRow * nProcCol)

    allocate(gridMap(nProcRow, nProcCol))
    ind = 1
    do iProcRow = 1, nProcRow
      do iProcCol = 1, nProcCol
        gridmap(iProcRow, iProcCol) = groupMembers(ind)
        ind = ind + 1
      end do
    end do

  end subroutine getGridMap


end module blacsenv
