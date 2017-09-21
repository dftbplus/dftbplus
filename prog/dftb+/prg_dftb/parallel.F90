#:include 'common.fypp'

module parallel
  use accuracy
  use environment
  use message
#:if WITH_SCALAPACK
  use scalapackfx
#:endif
  implicit none
  private

#:if WITH_SCALAPACK
  public :: TBlacsGrids, TBlacsGrids_init
#:endif


#:if WITH_SCALAPACK

  !> Contains various BLACS related options
  type :: TBlacsGrids

    !> Grid containing all processes.
    type(blacsgrid) :: gridAll

    !> Group grid for processing square matrices of the size nOrb x nOrb
    type(blacsgrid) :: gridOrbSqr

    !> Group grid for processing square matrices of the size nAtom x nAtom
    type(blacsgrid) :: gridAtomSqr

    !> Nr. of processor groups
    integer :: nGroup

    !> Nr. of processor within each group
    integer :: groupSize

    !> Contains k-point (1, ii) and spin (2, ii) tuples to be processed in current group.
    integer, allocatable :: groupKS(:,:)

  end type TBlacsGrids

#:endif

  
contains

#:if WITH_SCALAPACK
  
  !> Initializes BLACS grids
  subroutine TBlacsGrids_init(this, rowBlock, colBlock, nGroup, nOrb, nAtom, nKpoint, nSpin,&
      & tPauliHS)

    !> Initialized instance at exit.
    type(TBlacsGrids), intent(out) :: this

    !> Row block size
    integer, intent(in) :: rowBlock

    !> Column block size
    integer, intent(in) :: colBlock

    !> Nr. of processor groups.
    integer, intent(in) :: nGroup

    !> Nr. of orbitals
    integer, intent(in) :: nOrb

    !> Nr. of atoms
    integer, intent(in) :: nAtom

    !> Nr. of K-points
    integer, intent(in) :: nKpoint

    !> Nr. of spin channels
    integer, intent(in) :: nSpin

    !> Whether we need a 2x2 Pauli type Hamiltonian and overlap
    logical, intent(in) :: tPauliHS

    integer :: nProcRow, nProcCol, maxDim, nProcRowmax, nProcColmax
    integer :: nProc, iProc, iGroup
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
    iGroup = iProc / this%groupSize
    call getSquareGridParams_(this%groupSize, nProcRow, nProcCol)
    if (.not. tPauliHS) then
      call getGroupKS_(nGroup, nKpoint, nSpin, iGroup, this%groupKS)
    else
      call getGroupKS_(nGroup, nKpoint, 1, iGroup, this%groupKS)
    end if
    
    ! Check whether all processes have some portions of the H and S matrices
    if (tPauliHS) then
      maxDim = 2 * nOrb
    else
      maxDim = nOrb
    end if
    nProcRowmax = (maxDim - 1) / rowBlock + 1
    nProcColmax = (maxDim - 1) / colblock + 1
    if (nProcRow > nProcRowmax .or. nProcCol > nProcColmax) then
      write(buffer, "(A,I0,A,I0,A,I0,A,I0,A)") "Processor grid (", nProcRow, " x ",  nProcCol,&
          & ") too big (> ", nProcRowmax, " x ", nProcColmax, ")"
      call error(buffer)
    end if

    ! Subgrid for communication inside the groups
    call this%gridOrbSqr%initgrids(nGroup, nProcRow, nProcCol)
    write(stdOut, "(1X,3(A,I0))") "PGRID:ORBITAL: ", nGroup, " x ", nProcRow, " x ", nProcCol
    
    ! Global grid for broadcasting messages to all processes
    call getSquareGridParams_(nProc, nProcRow, nProcCol)
    call this%gridAll%initgrid(nProcRow, nProcCol)
    write(stdOut, "(1X,2(A,I0))") "PGRID:ALLPROC: ", nProcRow, " x ", nProcCol

    ! Create grid for atomic quantities
    nProcRowmax = (nAtom - 1) / rowBlock + 1
    nProcColmax = (nAtom - 1) / colBlock + 1
    nProcRow = min(nProcRow, nProcRowmax)
    nProcCol = min(nProcCol, nProcColmax)
    call this%gridAtomSqr%initgrid(nProcRow, nProcCol)
    write(stdOut, "(1X,2(A,I0))") "PGRID:ATOM: ", nProcRow, " x ", nProcCol

  end subroutine TBlacsGrids_init

#:endif


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!! Private routines
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    !> Returns rows and columns for a 2D processor grid closest possible to a square one.
  subroutine getSquareGridParams_(nProc, nRow, nCol)

    !> Total number of processors
    integer, intent(in) :: nProc

    !> Number of rows in the 2D processor grid
    integer, intent(out) :: nRow

    !> Number of columns in the 2D processor grid
    integer, intent(out) :: nCol

    do nRow = int(sqrt(real(nProc, dp))), 1, -1
      if (mod(nProc, nRow) == 0) then
        exit
      end if
    end do
    nCol = nProc / nRow

  end subroutine getSquareGridParams_


  !> Returns the (k-point, spin) tuples to be processed by current processor grid.
  subroutine getGroupKS_(nGroup, nKpoint, nSpin, iGroup, groupKS)

    !> Number of processor groups
    integer, intent(in) :: nGroup

    !> Number of k-points in calculation.
    integer, intent(in) :: nKpoint

    !> Number of spin channels in calculation
    integer, intent(in) :: nSpin

    !> Group index of current group
    integer, intent(in) :: iGroup

    !> Array of (k-point, spin) tuples (groupKS(:, ii) = [iK, iS])
    integer, intent(out), allocatable :: groupKS(:,:)

    integer :: nHam, nHamAll, res
    integer :: ind, iHam, iS, iK
    
    nHamAll = nKpoint * nSpin
    nHam = nHamAll / nGroup
    res = nHamAll - nHam * nGroup
    if (iGroup < res) then
      nHam = nHam + 1
    end if
    allocate(groupKS(2, nHam))
    ind = 0
    iHam = 1
    do iS = 1, nSpin
      do iK = 1, nKpoint
        if (mod(ind, nGroup) == iGroup) then
          groupKS(1, iHam) = iK
          groupKS(2, iHam) = iS
          iHam = iHam + 1
        end if
        ind = ind + 1
      end do
    end do

  end subroutine getGroupKS_


end module parallel
