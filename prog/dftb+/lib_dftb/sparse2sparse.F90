!--------------------------------------------------------------------------------------------------!
!  DFTB+: general package for performing fast atomistic simulations                                !
!  Copyright (C) 2018  DFTB+ developers group                                                      !
!                                                                                                  !
!  See the LICENSE file for terms of usage and distribution.                                       !
!--------------------------------------------------------------------------------------------------!

#:include 'common.fypp'

!> Contains subroutines for packing/unpacking Hamiltonian-like matrices between the DFTB+ compressed
!> block sparse representation and distributed CSR/CSC matrices
module sparse2sparse
  use assert
  use accuracy
  use commontypes
  use environment
  use commontypes
  use solvers
  use message
#:if WITH_MPI
  use mpifx
#:endif
#:if WITH_ELSI
  use elsi
#:endif
  implicit none

  private
#:if WITH_ELSI
  public :: calcdensity_parallel_elsi
#:endif

contains


#:if WITH_ELSI
  !> Calculates density matrix using the elsi routine.
  subroutine calcdensity_parallel_elsi(env, parallelKS, electronicSolver, ham, over, iNeighbour,&
      & nNeighbourSK, iAtomStart, iSparseStart, img2CentCell, orb, rho, Eband)

    !> Environment settings
    type(TEnvironment), intent(in) :: env

    !> Contains (iK, iS) tuples to be processed in parallel by various processor groups
    type(TParallelKS), intent(in) :: parallelKS

    !> Electronic solver information
    type(TElectronicSolver), intent(inout) :: electronicSolver

    !> hamiltonian in sparse storage
    real(dp), intent(in) :: ham(:,:)

    !> overlap matrix in sparse storage
    real(dp), intent(in) :: over(:)

    !> Neighbour list for the atoms (First index from 0!)
    integer, intent(in) :: iNeighbour(0:,:)

    !> Nr. of neighbours for the atoms.
    integer, intent(in) :: nNeighbourSK(:)

    !> Atom offset for the squared matrix
    integer, intent(in) :: iAtomStart(:)

    !> indexing array for the sparse Hamiltonian
    integer, intent(in) :: iSparseStart(0:,:)

    !> Mapping between image atoms and corresponding atom in the central cell.
    integer, intent(in) :: img2CentCell(:)

    !> data structure with atomic orbital information
    type(TOrbitals), intent(in) :: orb

    !> Density matrix in DFTB+ sparse format
    real(dp), intent(out) :: rho(:,:)

    !> Band energy
    real(dp), intent(out) :: Eband(:)

    integer :: nn
    integer :: nKS, iKS, iS
    integer :: nnz_global, nnzLocal, numColLocal
    integer, allocatable :: colptrLocal(:), rowindLocal(:)
    real(dp), allocatable :: HnzvalLocal(:), SnzvalLocal(:)
    real(dp), allocatable :: DMnzvalLocal(:)
    integer :: colStartLocal, colEndLocal, numCol, numRow, numColDefault
    integer, allocatable :: blockrow(:)

    nKS = size(parallelKS%localKS, dim=2)

    numRow = electronicSolver%ELSI_n_basis
    numCol = electronicSolver%ELSI_n_basis
    numColDefault = int(numCol / env%mpi%groupComm%size)

    colStartLocal = env%mpi%groupComm%rank * numColDefault + 1
    if (env%mpi%groupComm%rank /= env%mpi%groupComm%size - 1) then
      numColLocal = numColDefault
    else
      numColLocal = numColDefault + modulo(numCol, env%mpi%groupComm%size)
    end if
    colEndLocal = colStartLocal + numColLocal - 1

    ! Could be stored in a derived type and reused between SCC iterations
    ALLOCATE(colptrLocal(numColLocal + 1))

    call pack2colptr_parallel(iNeighbour, nNeighbourSK, iAtomStart, iSparseStart, &
        & img2CentCell, colStartLocal, colEndLocal, nnzLocal, colptrLocal)

    nnz_global = 0
    call mpifx_allreduce(env%mpi%groupComm, nnzLocal, nnz_global, MPI_SUM)

    ALLOCATE(rowindLocal(nnzLocal))
    ALLOCATE(HnzvalLocal(nnzLocal))
    ALLOCATE(SnzvalLocal(nnzLocal))
    ALLOCATE(DMnzvalLocal(nnzLocal))

    nn = size(ham, 1)
    ALLOCATE(blockRow(nn))

    ! Could be stored in a derived type end reused between SCC iterations
    call pack2elsi_parallel(over, iNeighbour, nNeighbourSK, iAtomStart, iSparseStart,&
        & img2CentCell, colStartLocal, colEndLocal, colptrLocal, rowindLocal, blockRow, SnzvalLocal)

    call elsi_set_csc(electronicSolver%elsiHandle, nnz_global, nnzLocal,&
        & colEndLocal-colStartLocal+1, rowindLocal, colptrLocal)

    do iKS = 1, nKS
      iS = parallelKS%localKS(2, iKS)
      call pack2elsi_parallel(ham(:,iS), iNeighbour, nNeighbourSK, iAtomStart, iSparseStart,&
          & img2CentCell, colStartLocal, colEndLocal, colptrLocal, rowindLocal, blockRow,&
          & HnzvalLocal)

      ! Load the matrix into the internal data structure
      call elsi_dm_real_sparse(electronicSolver%elsiHandle, HnzvalLocal, SnzvalLocal, DMnzvalLocal,&
          & Eband(iS))

      ! get results back
      call elsi2pack_parallel(colStartLocal, colEndLocal, iNeighbour, nNeighbourSK, orb%mOrb,&
          & iAtomStart, iSparseStart, img2CentCell, colptrLocal, DMnzvalLocal, blockRow, rho(:,iS))

    end do

    write(*,*)'Eband',EBand

  end subroutine calcdensity_parallel_elsi


  !> Creating colptr and nnz for CSC matrix format from packed format
  subroutine pack2colptr_parallel(iNeighbour, nNeighbourSK, iSquare, iSparseStart, img2CentCell,&
      & colStartLocal, colEndLocal, nnzLocal, colptrLocal)

    !> Neighbour list for each atom (first index from 0!).
    integer, intent(in) :: iNeighbour(0:,:)

    !> Nr. of neighbours for each atom (incl. itself).
    integer, intent(in) :: nNeighbourSK(:)

    !> Atom offset for the squared Hamiltonian
    integer, intent(in) :: iSquare(:)

    !> Indexing array for the sparse Hamiltonian
    integer, intent(in) :: iSparseStart(0:,:)

    !> Atomic mapping indexes
    integer, intent(in) :: img2CentCell(:)

    !> Column of global matrix where local matrix starts
    integer, intent(in) :: colStartLocal

    !> Column of global matrix where local matrix ends
    integer, intent(in) :: colEndLocal

    !> Local number of  non-zero elements
    integer, intent(out) :: nnzLocal

    !> Local column pointer
    integer, intent(out) :: colptrLocal(:)

    integer :: nAtom
    integer :: iOrig, ii, jj, nOrb1, nOrb2
    integer :: iNeigh
    integer :: iAtom1, iAtom2, iAtom2f, kk

    logical, allocatable :: blockList(:)

    nAtom = size(iNeighbour, dim=2)

    colptrLocal(:) = 0
    nnzLocal = 0

    ALLOCATE(blockList(nAtom))

    ! Loop over all atom blocks
    do iAtom1 = 1, nAtom
      ii = iSquare(iAtom1)
      nOrb1 = iSquare(iAtom1+1) - ii
      blockList(:) = .false.
      do iNeigh = 0, nNeighbourSK(iAtom1)
        iOrig = iSparseStart(iNeigh,iAtom1) + 1
        iAtom2 = iNeighbour(iNeigh, iAtom1)
        iAtom2f = img2CentCell(iAtom2)
        jj = iSquare(iAtom2f)
        nOrb2 = iSquare(iAtom2f+1) - jj

        if (blockList(iAtom2f)) then
          cycle
        else
          blockList(iAtom2f) = .true.
        end if

        if (isBlockInLocal(jj, jj+nOrb2-1, ii, ii+nOrb1-1, colStartLocal, colEndLocal)) then
          do kk = ii, ii + nOrb1 - 1
            if (isColumnInLocal(kk, colStartLocal, colEndLocal)) then
              nnzLocal = nnzLocal + nOrb2
              colptrLocal(kk - colStartLocal + 2) = colptrLocal(kk - colStartLocal + 2) + nOrb2
            end if
          end do
        end if

        if (ii == jj) then
          cycle
        end if

        do kk = jj, jj + nOrb2 - 1
          if (isColumnInLocal(kk, colStartLocal, colEndLocal)) then
            nnzLocal = nnzLocal + nOrb1
            colptrLocal(kk -colStartLocal + 2) = colptrLocal(kk - colStartLocal + 2) + nOrb1
          end if
        end do
      end do
    end do
    colptrLocal(1) = 1
    do ii = 2, size(colptrLocal)
      colptrLocal(ii) = colptrLocal(ii-1) + colptrLocal(ii)
    end do
    colptrLocal(size(colptrLocal)) = nnzLocal + 1

  end subroutine pack2colptr_parallel


  !> Convert sparse DFTB+ matrix to distributed CSC matrix format for ELSI calculations
  !>
  !> NOTE: ELSI needs the full matrix (both triangles)
  subroutine pack2elsi_parallel(orig, iNeighbour, nNeighbourSK, iSquare, iSparseStart,&
      & img2CentCell, colStartLocal, colEndLocal, colptrLocal, rowindLocal, blockRow, nzvalLocal)

    !> Sparse Hamiltonian
    real(dp), intent(in) :: orig(:)

    !> Neighbour list for each atom (first index from 0!).
    integer, intent(in) :: iNeighbour(0:,:)

    !> Nr. of neighbours for each atom (incl. itself).
    integer, intent(in) :: nNeighbourSK(:)

    !> Atom offset for the squared Hamiltonian
    integer, intent(in) :: iSquare(:)

    !> Indexing array for the sparse Hamiltonian
    integer, intent(in) :: iSparseStart(0:,:)

    !> Atomic mapping indexes.
    integer, intent(in) :: img2CentCell(:)

    !> Column of global matrix where local matrix starts
    integer, intent(in) :: colStartLocal

    !> Column of global matrix where local matrix ends
    integer, intent(in) :: colEndLocal

    !> Local column pointer
    integer, intent(in) :: colptrLocal(:)

    !> Local row index pointer
    integer, intent(out) :: rowindLocal(:)

    !> Saves starting row of blocks in nzvalLocal
    integer, intent(inout) :: blockRow(:)

    !> Local non-zero elements
    real(dp), intent(out) :: nzvalLocal(:)

    integer :: nAtom
    integer :: iOrig, ii, jj, nOrb1, nOrb2
    integer :: iNeigh
    integer :: iAtom1, iAtom2, iAtom2f

    integer :: iNext
    integer, allocatable :: blockList(:,:)
    logical, allocatable :: tRowTrans(:)

    nAtom = size(iNeighbour, dim=2)

    @:ASSERT(nAtom > 0)
    @:ASSERT(size(nNeighbourSK) == nAtom)
    @:ASSERT(size(iSquare) == nAtom + 1)
    @:ASSERT(size(orig) == size(blockRow))

    ALLOCATE(blockList(nAtom,2))
    ALLOCATE(tRowTrans(nAtom))
    rowindLocal(:) = 0
    nzvalLocal(:) = 0.0_dp

    blockRow(:) = 0

    ! Offset in column belonging to transposed triangle
    blockList(:,2) = 1
    do iAtom1 = 1, nAtom
      ii = iSquare(iAtom1)
      nOrb1 = iSquare(iAtom1+1) - ii
      ! Offset in current column
      blockList(:,1) = 0
      tRowTrans = .false.
      ! Starting index, because probaly column already contains blocks comming from transposing when
      ! processing previous columns.
      iNext = blockList(iAtom1, 2)
      do iNeigh = 0, nNeighbourSK(iAtom1)
        iOrig = iSparseStart(iNeigh,iAtom1) + 1
        iAtom2 = iNeighbour(iNeigh, iAtom1)
        iAtom2f = img2CentCell(iAtom2)

        jj = iSquare(iAtom2f)
        @:ASSERT(jj >= ii)
        nOrb2 = iSquare(iAtom2f+1) - jj

        if (blockList(iAtom2f,1) == 0) then
          blockList(iAtom2f,1) = iNext
          iNext = iNext + nOrb2
        end if

        blockRow(iOrig) = blockList(iAtom2f,1)

        if ( .not. isBlockInLocal(jj,jj+nOrb2-1, ii, &
            & ii+nOrb1-1, colStartLocal, colEndLocal)) then
          cycle
        end if
        call addBlock2Elsi(reshape(orig(iOrig:iOrig+nOrb1*nOrb2-1), [ nOrb2, nOrb1 ]),&
            & colStartLocal, colEndLocal, jj, ii, colptrLocal, blockList(iAtom2f,1)-1, rowindLocal,&
            & nzvalLocal)

        if (ii == jj) then
          cycle
        end if

        ! Because of folding of periodic images, it can happen, that transposed block has already
        ! been registered.
        if (.not. tRowTrans(iAtom2f)) then
          blockList(iAtom2f,2) = blockList(iAtom2f,2) + nOrb1
          tRowTrans(iAtom2f) = .true.
        end if

        call addBlock2Elsi(transpose(reshape(orig(iOrig:iOrig+nOrb1*nOrb2-1), [ nOrb2, nOrb1 ])),&
            & colStartLocal, colEndLocal, ii, jj, colptrLocal, blockList(iAtom2f,2)-nOrb1-1,&
            & rowindLocal, nzvalLocal)
      end do
    end do

  end subroutine pack2elsi_parallel


  !> Checks if atom block is part of local matrix (Elsi)
  pure logical function isBlockInLocal(rowStartBlock, rowEndBlock, colStartBlock, colEndBlock,&
      & colStartLocal, colEndLocal)

    !> Row of global matrix where block starts
    integer, intent(in) :: rowStartBlock

    !> Row of global matrix where block ends
    integer, intent(in) :: rowEndBlock

    !> Column of global matrix where block starts
    integer, intent(in) :: colStartBlock

    !> Column of global matrix where block ends
    integer, intent(in) :: colEndBlock

    !> Column of global matrix where local part starts
    integer, intent(in) :: colStartLocal

    !> Column of global matrix where local part ends
    integer, intent(in) :: colEndLocal

    isBlockInLocal = .false.
    if ( (colStartLocal <= colStartBlock) .and. (colStartBlock <= colEndLocal) ) then
      isBlockInLocal = .true.
    else if ( (colStartLocal <= colEndBlock) .and. (colEndBlock <= colEndLocal) ) then
      isBlockInLocal = .true.
    else if ( (colStartBlock <= colStartLocal) .and. (colStartLocal <= colEndBlock)) then
      isBlockInLocal = .true.
    else if ( (colStartBlock <= colEndLocal) .and. (colEndLocal <= colEndBlock)) then
      isBlockInLocal = .true.
    else if ( (colStartLocal <= rowStartBlock) .and. (rowStartBlock <= colEndLocal) ) then
      isBlockInLocal = .true.
    else if ( (colStartLocal <= rowEndBlock) .and. (rowEndBlock <= colEndLocal) ) then
      isBlockInLocal = .true.
    else if ( (rowStartBlock <= colStartLocal) .and. (colStartLocal <= rowEndBlock)) then
      isBlockInLocal = .true.
    else if ( (rowStartBlock <= colEndLocal) .and. (colEndLocal <= rowEndBlock)) then
      isBlockInLocal = .true.
    end if

  end function isBlockInLocal


  !> Checks whether global column is in local matrix
  pure logical function isColumnInLocal(col, colStartLocal, colEndLocal)

    !> global column
    integer, intent(in) :: col

    !> Column of global matrix where local matrix starts
    integer, intent(in) :: colStartLocal

    !> Column of global matrix where local matrix ends
    integer, intent(in) :: colEndLocal

    if ( (colStartLocal <= col) .and. (col <= colEndLocal) ) then
      isColumnInLocal = .true.
    else
      isColumnInLocal = .false.
    end if

  end function isColumnInLocal


  !> Add the content of a local matrix block to ELSI CSC format
  subroutine addBlock2Elsi(loc, colStart, colEnd, ii, jj, colptr, rowOffset, rowind, val)

    !> Local matrix.
    real(dp), intent(in) :: loc(:,:)

    !> Column of global matrix where local matrix starts
    integer, intent(in) :: colStart

    !> Column of global matrix where local matrix ends
    integer, intent(in) :: colEnd

    !> Starting row in the global matrix.
    integer, intent(in) :: ii

    !> Starting column in the global matrix
    integer, intent(in) :: jj

    !> column pointer
    integer, intent(in) :: colptr(:)

    !> index of the next element per column
    integer, intent(in) :: rowOffset

    !> row index pointer
    integer, intent(inout) :: rowind(:)

    !> values in CSC format
    real(dp), intent(inout) :: val(:)

    integer :: j2, iloc, jloc
    integer :: jStart, jEnd

    jStart = max(jj, colStart) - colStart + 1
    jEnd = min(jj+size(loc,dim=2)-1, colEnd) - colStart + 1
    jloc = max(jj, colStart) - jj + 1

    do j2 = jStart, jEnd
      do iloc = 1, size(loc,1)
        val(colptr(j2)+rowOffset+iloc-1) = val(colptr(j2) + rowOffset+iloc-1) + loc(iloc,jloc)
        rowind(colptr(j2) + rowOffset + iloc-1) = iloc + ii - 1
      end do
      jloc = jloc + 1
    end do

  end subroutine addBlock2Elsi


  !> Creating colptr and nnz for CSC matrix format from packed format
  !>
  !> Note: primitive will not be set to zero on startup, and values are added to enable addition of
  !> spin components. Make sure, you set it to zero before invoking this routine the first time.
  subroutine elsi2pack_parallel(colStart, colEnd, iNeighbour, nNeighbourSK, mOrb, iSquare,&
      & iSparseStart, img2CentCell, colptr, nzval, blockRow, primitive)

    !> Column of global matrix where local matrix starts
    integer, intent(in) :: colStart

    !> Column of global matrix where local matrix ends
    integer, intent(in) :: colEnd

    !> Neighbour list for each atom (first index from 0!).
    integer, intent(in) :: iNeighbour(0:,:)

    !> Nr. of neighbours for each atom (incl. itself).
    integer, intent(in) :: nNeighbourSK(:)

    !> Maximal number of orbitals on an atom.
    integer, intent(in) :: mOrb

    !> Atom offset for the squared Hamiltonian
    integer, intent(in) :: iSquare(:)

    !> Indexing array for the sparse Hamiltonian
    integer, intent(in) :: iSparseStart(0:,:)

    !> Atomic mapping indexes.
    integer, intent(in) :: img2CentCell(:)

    !> Local column pointer
    integer, intent(in) :: colptr(:)

    !> Local number of  non-zero elements
    real(dp), intent(in) :: nzval(:)

    !> Saves starting row of blocks in nzvalLocal
    integer, intent(in) :: blockRow(:)

    !> Sparse Hamiltonian
    real(dp), intent(inout) :: primitive(:)

    integer :: nAtom
    integer :: iOrig, ii, jj, kk
    integer :: iNeigh
    integer :: iAtom1, iAtom2, iAtom2f
    integer :: nOrb1, nOrb2
    real(dp) :: tmpSqr(mOrb,mOrb)

    nAtom = size(iNeighbour, dim=2)

    tmpSqr(:,:) = 0.0_dp

    do iAtom1 = 1, nAtom
      ii = iSquare(iAtom1)
      nOrb1 = iSquare(iAtom1+1) - ii
      do iNeigh = 0, nNeighbourSK(iAtom1)
        iOrig = iSparseStart(iNeigh,iAtom1) + 1
        iAtom2 = iNeighbour(iNeigh, iAtom1)
        iAtom2f = img2CentCell(iAtom2)
        jj = iSquare(iAtom2f)
        nOrb2 = iSquare(iAtom2f+1) - jj

        if ( .not. isBlockInLocal(jj,jj+nOrb2-1, ii, ii+nOrb1-1, colStart, colEnd)) then
          cycle
        end if

        call cpElsi2Block(colStart, colEnd, colptr, nzval, ii, blockRow(iOrig),&
            & tmpSqr(1:nOrb2,1:nOrb1))

        ! Symmetrize the on-site block before packing, just in case
        if (iAtom1 == iAtom2f) then
          do kk = 1, nOrb2
            tmpSqr(kk, kk+1:nOrb1) = tmpSqr(kk+1:nOrb1, kk)
          end do
        end if

        primitive(iOrig : iOrig + nOrb1 * nOrb2 - 1) = primitive(iOrig : iOrig + nOrb1 * nOrb2 - 1)&
            & + reshape(tmpSqr(1:nOrb2,1:nOrb1), [ nOrb1*nOrb2 ])

      end do
    end do

  end subroutine elsi2pack_parallel


  !> Copies the content from the ELSI structure to block
  subroutine cpElsi2Block(colStart, colEnd, colptr, nzval, jj, blockRow, loc)

    !> Column of global matrix where local matrix starts
    integer, intent(in) :: colStart

    !> Column of global matrix where local matrix ends
    integer, intent(in) :: colEnd

    !> column pointer
    integer, intent(in) :: colptr(:)

    !> non-zero values
    real(dp), intent(in) :: nzval(:)

    !> Starting column in the global matrix
    integer, intent(in) :: jj

    !> Saves starting row of blocks in nzvalLocal
    integer, intent(in) :: blockRow

    !> Local matrix.
    real(dp), intent(out) :: loc(:,:)

    integer :: j2, iloc

    loc(:,:) = 0.0_dp

    do j2 = 1, size(loc, dim=2)
      if ( isColumnInLocal(jj + j2 - 1, colStart, colEnd )) then
        iloc = blockRow + colptr(jj-colStart+j2) - 1
        loc(:, j2) = nzval(iloc:iloc+size(loc,dim=1)-1)
      end if
    end do

  end subroutine cpElsi2Block
#:endif


end module sparse2sparse
