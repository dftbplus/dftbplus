!--------------------------------------------------------------------------------------------------!
!  DFTB+: general package for performing fast atomistic simulations                                !
!  Copyright (C) 2017  DFTB+ developers group                                                      !
!                                                                                                  !
!  See the LICENSE file for terms of usage and distribution.                                       !
!--------------------------------------------------------------------------------------------------!

#:include 'common.fypp'

!> Contains BLACS environmental settings.
module blacsenv
  use globalenv, only : stdOut
  use message
  use scalapackfx
  use accuracy, only : lc
  implicit none
  private

  public :: TBlacsEnv, TBlacsEnv_init


  !> Contains various BLACS related settings
  type :: TBlacsEnv

    !> Grid containing all processes.
    type(blacsgrid) :: gridAll

    !> Group grid for processing square matrices of the size nOrb x nOrb
    type(blacsgrid) :: gridOrbSqr

    !> Group grid for processing square matrices of the size nAtom x nAtom
    type(blacsgrid) :: gridAtomSqr

    !> Nr. of processor groups
    integer :: nGroup

    !> Index of the group current process is in
    integer :: iGroup

    !> Nr. of processor within each group
    integer :: groupSize

  end type TBlacsEnv


contains


  !> Initializes BLACS grids
  subroutine TBlacsEnv_init(this, rowBlock, colBlock, nGroup, nOrb, nAtom)

    !> Initialized instance at exit.
    type(TBlacsEnv), intent(out) :: this

    !> Row block size
    integer, intent(in) :: rowBlock

    !> Column block size
    integer, intent(in) :: colBlock

    !> Nr. of processor groups
    integer, intent(in) :: nGroup

    !> Nr. of orbitals
    integer, intent(in) :: nOrb

    !> Nr. of atoms
    integer, intent(in) :: nAtom

    integer :: nProcRow, nProcCol, maxProcRow, maxProcColMax
    integer :: nProc, iProc
    character(lc) :: buffer

    call blacsfx_pinfo(iProc, nProc)
    if (nGroup < 1 .or. nGroup > nProc) then
      call error("Nr. of  groups must be between 1 and nr. of processes")
    end if
    if (mod(nProc, nGroup) /= 0) then
      call error("Nr. of groups must be a divisor of nr. of processes")
    end if

    this%nGroup = nGroup
    this%groupSize = nProc / nGroup
    this%iGroup = iProc / this%groupSize

    ! Check whether all processes have some portions of the H and S matrices
    call getSquareGridParams(this%groupSize, nProcRow, nProcCol)
    maxProcRow = (nOrb - 1) / rowBlock + 1
    maxProcColMax = (nOrb - 1) / colblock + 1
    if (nProcRow > maxProcRow .or. nProcCol > maxProcColMax) then
      write(buffer, "(A,I0,A,I0,A,I0,A,I0,A)") "Processor grid (", nProcRow, " x ",  nProcCol,&
          & ") too big (> ", maxProcRow, " x ", maxProcColMax, ")"
      call error(buffer)
    end if

    ! Subgrid for communication inside the groups
    call this%gridOrbSqr%initgrids(nGroup, nProcRow, nProcCol)
    write(stdOut, "(1X,3(A,I0))") "PGRID:ORBITAL: ", nGroup, " x ", nProcRow, " x ", nProcCol

    ! Global grid for broadcasting messages to all processes
    call getSquareGridParams(nProc, nProcRow, nProcCol)
    call this%gridAll%initgrid(nProcRow, nProcCol)
    write(stdOut, "(1X,2(A,I0))") "PGRID:ALLPROC: ", nProcRow, " x ", nProcCol

    ! Create grid for atomic quantities
    maxProcRow = (nAtom - 1) / rowBlock + 1
    maxProcColMax = (nAtom - 1) / colBlock + 1
    nProcRow = min(nProcRow, maxProcRow)
    nProcCol = min(nProcCol, maxProcColMax)
    call this%gridAtomSqr%initgrid(nProcRow, nProcCol)
    write(stdOut, "(1X,2(A,I0))") "PGRID:ATOM: ", nProcRow, " x ", nProcCol

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


end module blacsenv
