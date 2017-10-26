!!* Contains routines for folding internal sparse matrixes into Compressed
!!* Sparse Row Format and vica-versa.
!!* @caveat All folding and unfoulding routines assume, that the CSR matrices
!!*   have the same sparsity structure as the internal sparse matrices. The
!!*   CSR matrices must be created using the <code>createEquivCSR</code>
!!*   routines.
#:include "common.fypp"

module mat_conv
  use assert
  use accuracy
  use libnegf, only: r_CSR, z_CSR, r_DNS, z_DNS, create, destroy 
  use Constants, only : pi
  use CommonTypes, only : TOrbitals
  implicit none
  private

  public :: init,   destruct
  public :: foldToCSR, unfoldFromCSR
  public :: TDftbMat_r, TDftbMat_z

  type TDftbMat_r
    real(dp), pointer    :: sparse(:)       !The sparse matrix to convert
    integer, pointer     :: nOrb(:)         !Orbitals per atom
    integer, pointer     :: iSquare(:)      !Starting positions of the atoms
    integer, pointer     :: iPair(:,:)      !Starting position of atom-neighbor
    integer, pointer     :: iNeighbor(:,:)  !Index of neighbors
    integer, pointer     :: nNeighbor(:)    !Number of neighbors
    integer, pointer     :: img2CentCell(:) !Image in the central cell
    integer, pointer     :: mAngAtom(:)     !Angular momentum of each atom
    integer              :: mmAng           !Maximal angular momentum 
  end type TDftbMat_r

  type TDftbMat_z
    complex(dp), pointer :: sparse(:)       ! The sparse matrix to convert
    integer, pointer     :: nOrb(:)         ! Orbitals per atom
    integer, pointer     :: iSquare(:)      ! Starting positions of the atoms
    integer, pointer     :: iPair(:,:)      ! Starting position of atom-neighbor
    integer, pointer     :: iNeighbor(:,:)  ! Index of neighbors
    integer, pointer     :: nNeighbor(:)    ! Number of neighbors
    integer, pointer     :: img2CentCell(:) ! Image in the central cell
    type(TOrbitals), pointer :: orb         ! Orbitals
  end type TDftbMat_z

  !!* Creates a Compressed Sparse Row matrix with equivalent sparsity structure
  !!* as the internal sparse matrixes.
  interface init
    module procedure createEmptyCSR_real
    module procedure createEmptyCSR_cplx
    module procedure createEquivCSR_real
    module procedure createEquivCSR_cplx
    module procedure cloneSparsityMap_real_real
    module procedure cloneSparsityMap_cplx_cplx
    module procedure cloneSparsityMap_real_cplx
    module procedure cloneSparsityMap_cplx_real
  end interface

  !!* Destructor for the different structures.
  interface destruct
    module procedure rdestroy_CSR, zdestroy_CSR
    module procedure rdestroy_DNS, zdestroy_DNS
  end interface

  !!* Folds internal sparse storage to Compressed Sparse Row format.
  interface foldToCSR
    module procedure foldToCSR_real
    module procedure foldToCSR_cplx
  end interface

  !!* Unfolds from Compressed Sparse Row format into internal sparse storage.
  interface unfoldFromCSR
    module procedure unfoldFromCSR_real
    module procedure unfoldFromCSR_cplx
  end interface

contains

  !!* Creates a Compressed Sparse Row matrix with a sparsity structure in
  !!* accordance with current geometry.
  !!* @param csr  Compressed sparse row matrix on exit
  !!* @param iAtomStart  Starting position of the atoms in the square H/S.
  !!* @param iNeighbor  Neighborlist for each atom.
  !!* @param nNeighbor  Number of neighbors for each atom.
  !!* @param img2CentCell  Folded image in the central cell for each atom.
  !!* @param orb  Orbitals in the system.
  !!* @note The subroutine only creates the structure of the CSR matrix, but
  !!*   does not fill up the nonzero values with anything usefull.
  !!* @note The resulting CSR matrix has storage for both triangles of the
  !!*   square matrix.
  subroutine createEquivCSR_real(csr, iAtomStart, iNeighbor, nNeighbor, &
      &img2CentCell, orb)
    type(r_CSR), intent(out) :: csr
    integer, intent(in) :: iAtomStart(:)
    integer, intent(in) :: iNeighbor(0:,:)
    integer, intent(in) :: nNeighbor(:)
    integer, intent(in) :: img2CentCell(:)
    type(TOrbitals), intent(in) :: orb

    integer, allocatable :: nColAtom(:), nCols(:), iNonZero(:)
    integer, allocatable :: rowpnt(:)
    logical, allocatable :: zero(:)
    integer :: nAtom, nNonZero
    integer :: iOrb1, iOrb2, nOrb1, nOrb2, iAt1, iAt2f, iNeigh
    integer :: iRow, iCol, ind, nrow
    integer :: ii, jj

    nAtom = size(orb%nOrbAtom)

    !! Count nr. of nonzero columns in the square (folded) form for each atom
    allocate(nColAtom(nAtom))  ! Nr. of nonzero columns for each atom
    allocate(zero(nAtom))
    nColAtom(:) = 0
    do iAt1 = 1, nAtom
      zero(:) = .true.
      nOrb1 = orb%nOrbAtom(iAt1)
      do iNeigh = 0, nNeighbor(iAt1)
        iAt2f = img2CentCell(iNeighbor(iNeigh,iAt1))
        if (zero(iAt2f)) then
          nColAtom(iAt1) = nColAtom(iAt1) + orb%nOrbAtom(iAt2f)
          zero(iAt2f) = .false.
          if (iAt1 /= iAt2f) then
            nColAtom(iAt2f) = nColAtom(iAt2f) + nOrb1
          end if
        end if
      end do
    end do

    !! Calculate CSR row pointers:
    !! A row for a certain orbital of a certain atom is as long, as the 
    !! nr. of columns determined previously for the atom.
    nrow = iAtomStart(nAtom+1) - 1
    allocate(rowpnt(nrow+1))
    rowpnt(1) = 1
    ind = 1
    do iAt1 = 1, nAtom
      nOrb1 = orb%nOrbAtom(iAt1)
      do iOrb1 = 1, nOrb1
        rowpnt(ind+iOrb1) = rowpnt(ind) + iOrb1 * nColAtom(iAt1)
      end do
      ind = ind + nOrb1
    end do
    deallocate(nColAtom)

    call create(csr, nrow, nrow, rowpnt(nrow + 1) - 1)
    csr%rowpnt = rowpnt
    deallocate(rowpnt)

    !! Nr. of CSR columns already filled
    allocate(nCols(csr%nRow))     
    nCols(:) = 0
    !! Index of the nonzero blocks
    allocate(iNonZero(nAtom))

    !! Loop over all atoms (over all atomic columns in the rectangular picture)
    lpAt1: do iAt1 = 1, nAtom
      iCol = iAtomStart(iAt1)
      nOrb1 = orb%nOrbAtom(iAt1)           ! Width of the atomic column

      !! Mark nonzero blocks in the folded matrix for the current atom
      nNonZero = 0
      zero(:) = .true.
      do iNeigh = 0, nNeighbor(iAt1)
        iAt2f = img2CentCell(iNeighbor(iNeigh,iAt1))
        if (zero(iAt2f)) then
          zero(iAt2f) = .false.
          nNonZero = nNonZero + 1
          iNonzero(nNonzero) = iAt2f
        end if
      end do

      !! Initialise CSR internal pointers according the non-zero blocks
      lpNonZero: do ii = 1, nNonZero
        iAt2f = iNonZero(ii)
        nOrb2 = orb%nOrbAtom(iAt2f)
        iRow = iAtomStart(iAt2f)

        !! Correspond rows of the current atomic block column to the appropriate
        !! partial rows in the lower triangle of the CSR matrix.
        do iOrb2 = 0, nOrb2 - 1
          jj = iRow + iOrb2
          ind = csr%rowpnt(jj) + nCols(jj)
          do iOrb1 = 0, nOrb1 - 1
            csr%colind(ind+iOrb1) = iCol + iOrb1
          end do
          nCols(jj) = nCols(jj) + nOrb1
        end do

        !! Correspond the columns of the current block to appropriate 
        !! partial rows in the upper triangle of the CSR matrix.
        if (iAt2f /= iAt1) then
          do iOrb1 = 0, nOrb1 - 1
            jj = iCol + iOrb1
            ind = csr%rowpnt(jj) + nCols(jj)
            do iOrb2 = 0, nOrb2 - 1
              csr%colind(ind+iOrb2) = iRow + iOrb2
            end do
            nCols(jj) = nCols(jj) + nOrb2
          end do
        end if
      end do lpNonZero
    end do lpAt1

    deallocate(zero)
    deallocate(nCols)
    deallocate(iNonZero)
    
  end subroutine createEquivCSR_real

  

  !!* Folds the internal sparse formatted matrix to the compressed sparse row
  !!* format (real version).
  !!* @param csr Resulting CSR matrix
  !!* @param sparse The internal sparse matrix to fold.
  !!* @param iAtomStart Starting positions of the atoms in the square matrix
  !!* @param iPair Starting position of atom-neighbor interaction in the sparse
  !!*   matrix.
  !!* @param iNeighbor Index of neighbors
  !!* @param nNeighbor Number of neighbors
  !!* @param img2CentCell Image of the atoms in the central cell.
  !!* @param orb  Orbitals in the system.
  !!* @note In the resulting CSR format both triangles of the matrix are set.
  !!* @caveat The routine should only applied on CSR matrixes created by the
  !!*   createEquiv_real subroutine, since it exploits the ordering of the
  !!*   column indexes.
  subroutine foldToCSR_real(csr, sparse, iAtomStart, iPair, iNeighbor, &
      &nNeighbor, img2CentCell, orb)
    type(r_CSR), intent(inout) :: csr
    real(dp), intent(in) :: sparse(:)
    integer, intent(in) :: iAtomStart(:)
    integer, intent(in) :: iPair(0:,:)
    integer, intent(in) :: iNeighbor(0:,:)
    integer, intent(in) :: nNeighbor(:)
    integer, intent(in) :: img2CentCell(:)
    type(TOrbitals), intent(in) :: orb

    integer, allocatable :: nCols(:)
    real(dp), allocatable :: tmpCol(:,:)
    integer :: nAtom
    integer :: iOrb1, nOrb1, nOrb2, iAt1, iAt2f, iNeigh
    integer :: iRow, iCol, iVal, ind
    integer :: ii, jj, kk

    nAtom = size(orb%nOrbAtom)

    csr%sorted = .false.

    allocate(tmpCol(csr%nRow, orb%mOrb))  ! One atomic block column.
    lpAt1: do iAt1 = 1, nAtom
      !! Unfold the columns for the current atom from the sparse matrix.
      nOrb1 = orb%nOrbAtom(iAt1)
      tmpCol(:,:) = 0.0_dp
      do iNeigh = 0, nNeighbor(iAt1)
        ind = iPair(iNeigh,iAt1) + 1
        iAt2f = img2CentCell(iNeighbor(iNeigh,iAt1))
        nOrb2 = orb%nOrbAtom(iAt2f)
        iRow = iAtomStart(iAt2f)
        tmpCol(iRow:iRow+nOrb2-1, 1:nOrb1) = &
            &tmpCol(iRow:iRow+nOrb2-1, 1:nOrb1) &
            &+ reshape(sparse(ind:ind+nOrb2*nOrb1-1), (/ nOrb2, nOrb1 /))
      end do

      !! Copy every column into the appropriate row in the upper triangle of
      !! the CSR matrix (the copy is masked by the sparsity structure stored
      !! in the CSR matrix)
      nOrb1 = orb%nOrbAtom(iAt1)
      iCol = iAtomStart(iAt1)
      do iOrb1 = 1, nOrb1
        ii = csr%rowpnt(iCol + iOrb1 - 1)
        jj = csr%rowpnt(iCol + iOrb1) - 1
        do kk=ii,jj
           csr%nzval(kk) = tmpCol(csr%colind(kk), iOrb1)
        enddo
      end do
    end do lpAt1
    deallocate(tmpCol)

    !! Fill the lower triangle of the CSR matrix
    allocate(nCols(csr%nRow))
    nCols(:) = 0
    do iRow = 1, csr%nRow
      !! Starting from the tail of the matrix row
      iVal = csr%rowpnt(iRow + 1) - 1
      iCol = csr%colind(iVal)
      do while (iVal >= csr%rowpnt(iRow) .and. iCol > iRow)
        csr%nzval(csr%rowpnt(iCol)+nCols(iCol)) = csr%nzval(iVal)
        nCols(iCol) = nCols(iCol) + 1
        iVal = iVal - 1
        iCol = csr%colind(iVal)
      end do
    end do
    deallocate(nCols)

  end subroutine foldToCSR_real



  !!* Unfolds a matrix from the CSR form into the internal sparse
  !!* representation (real version).
  !!* @param sparse Updated sparse matrix.
  !!* @param csr CSR matrix
  !!* @param iAtomStart Starting positions of the atoms in the square matrix
  !!* @param iPair Starting position of atom-neighbor interaction in the sparse
  !!*   matrix.
  !!* @param iNeighbor Index of neighbors
  !!* @param nNeighbor Number of neighbors
  !!* @param img2CentCell Image of the atoms in the central cell.
  !!* @param orb  Orbitals in the system.
  !!* @note The CSR matrix must be symmetric.
  !!* @note The unfolded matrix is added to the passed sparse matrix.
  subroutine unfoldFromCSR_real(sparse, csr, iAtomStart, iPair, iNeighbor, &
      &nNeighbor, img2CentCell, orb)
    real(dp), intent(inout) :: sparse(:)
    type(r_CSR), intent(in) :: csr
    integer, intent(in) :: iAtomStart(:)
    integer, intent(in) :: iPair(0:,:)
    integer, intent(in) :: iNeighbor(0:,:)
    integer, intent(in) :: nNeighbor(:)
    integer, intent(in) :: img2CentCell(:)
    type(TOrbitals), intent(in) :: orb

    real(dp), allocatable :: tmpCol(:,:)
    integer :: nAtom
    integer :: iOrb1, nOrb1, nOrb2, iAt1, iAt2, iAt2f, iNeigh
    integer :: iRow, iCol, ind
    integer :: ii

    nAtom = size(orb%nOrbAtom)

    @:ASSERT(csr%nRow == iAtomStart(nAtom+1) - 1)

    allocate(tmpCol(csr%nRow, orb%mOrb))

    do iAt1 = 1, nAtom
      !! Put the rows belonging to a certain atom into the appropriate column
      !! of the block column. (Matrix is assumed to be symmetric!)
      tmpCol(:,:) = 0.0_dp
      nOrb1 = orb%nOrbAtom(iAt1)
      iCol = iAtomStart(iAt1)
      do iOrb1 = 1, nOrb1
        ii = iCol + iOrb1 - 1
        do ind = csr%rowpnt(ii), csr%rowpnt(ii+1) - 1
          tmpCol(csr%colind(ind), iOrb1) = real(csr%nzval(ind))
        end do
      end do

      !! Unfold the block column into the internal sparse format
      do iNeigh = 0, nNeighbor(iAt1)
        ind = iPair(iNeigh,iAt1) + 1
        iAt2 = iNeighbor(iNeigh, iAt1)
        iAt2f = img2CentCell(iAt2)
        nOrb2 = orb%nOrbAtom(iAt2f)
        iRow = iAtomStart(iAt2f)
        sparse(ind:ind+nOrb2*nOrb1-1) = &
            &sparse(ind:ind+nOrb2*nOrb1-1) &
            &+ reshape(tmpCol(iRow:iRow+nOrb2-1, 1:nOrb1), (/nOrb2*nOrb1/))
      end do
    end do

    deallocate(tmpCol)

  end subroutine unfoldFromCSR_real



  !!* Creates a Compressed Sparse Row matrix with a sparsity structure in
  !!* accordance with current geometry.
  !!* @param csr  Compressed sparse row matrix on exit
  !!* @param iAtomStart  Starting position of the atoms in the square H/S.
  !!* @param iNeighbor  Neighborlist for each atom.
  !!* @param nNeighbor  Number of neighbors for each atom.
  !!* @param img2CentCell  Folded image in the central cell for each atom.
  !!* @param orb  Orbitals in the system.
  !!* @note The subroutine only creates the structure of the CSR matrix, but
  !!*   does not fill up the nonzero values with anything usefull.
  !!* @note The resulting CSR matrix has storage for both triangles of the
  !!*   square matrix.
  subroutine createEquivCSR_cplx(csr, iAtomStart, iNeighbor, nNeighbor, &
      &img2CentCell, orb)
    type(z_CSR), intent(out) :: csr
    integer, intent(in) :: iAtomStart(:)
    integer, intent(in) :: iNeighbor(0:,:)
    integer, intent(in) :: nNeighbor(:)
    integer, intent(in) :: img2CentCell(:)
    type(TOrbitals), intent(in) :: orb

    integer, allocatable :: nColAtom(:), nCols(:), iNonZero(:)
    integer, allocatable :: rowpnt(:)
    logical, allocatable :: zero(:)
    integer :: nAtom, nNonZero
    integer :: iOrb1, iOrb2, nOrb1, nOrb2, iAt1, iAt2f, iNeigh
    integer :: iRow, iCol, ind, nrow
    integer :: ii, jj

    nAtom = size(orb%nOrbAtom)

    !! Count nr. of nonzero columns in the square (folded) form for each atom
    allocate(nColAtom(nAtom))  ! Nr. of nonzero columns for each atom
    allocate(zero(nAtom))
    nColAtom(:) = 0
    do iAt1 = 1, nAtom
      zero(:) = .true.
      nOrb1 = orb%nOrbAtom(iAt1)
      do iNeigh = 0, nNeighbor(iAt1)
        iAt2f = img2CentCell(iNeighbor(iNeigh,iAt1))
        if (zero(iAt2f)) then
          nColAtom(iAt1) = nColAtom(iAt1) + orb%nOrbAtom(iAt2f)
          zero(iAt2f) = .false.
          if (iAt1 /= iAt2f) then
            nColAtom(iAt2f) = nColAtom(iAt2f) + nOrb1
          end if
        end if
      end do
    end do

    nrow = iAtomStart(nAtom+1) - 1
    allocate(rowpnt(nrow+1) ) 
    !! Calculate CSR row pointers:
    !! A row for a certain orbital of a certain atom is as long, as the 
    !! nr. of columns determined previously for the atom.
    rowpnt(1) = 1
    ind = 1
    do iAt1 = 1, nAtom
      nOrb1 = orb%nOrbAtom(iAt1)
      do iOrb1 = 1, nOrb1
        rowpnt(ind+iOrb1) = rowpnt(ind) + iOrb1 * nColAtom(iAt1)
      end do
      ind = ind + nOrb1
    end do
    deallocate(nColAtom)

    call create(csr, nrow, nrow, rowpnt(nrow + 1) - 1)
    csr%rowpnt = rowpnt
    deallocate(rowpnt)

    !! Nr. of CSR columns already filled
    allocate(nCols(csr%nRow))     
    nCols(:) = 0
    !! Index of the nonzero blocks
    allocate(iNonZero(nAtom))

    !! Loop over all atoms (over all atomic columns in the rectangular picture)
    lpAt1: do iAt1 = 1, nAtom
      iCol = iAtomStart(iAt1)
      nOrb1 = orb%nOrbAtom(iAt1)           ! Width of the atomic column

      !! Mark nonzero blocks in the folded matrix for the current atom
      nNonZero = 0
      zero(:) = .true.
      do iNeigh = 0, nNeighbor(iAt1)
        iAt2f = img2CentCell(iNeighbor(iNeigh,iAt1))
        if (zero(iAt2f)) then
          zero(iAt2f) = .false.
          nNonZero = nNonZero + 1
          iNonzero(nNonzero) = iAt2f
        end if
      end do

      !! Initialise CSR internal pointers according the non-zero blocks
      lpNonZero: do ii = 1, nNonZero
        iAt2f = iNonZero(ii)
        nOrb2 = orb%nOrbAtom(iAt2f)
        iRow = iAtomStart(iAt2f)

        !! Correspond rows of the current atomic block column to the appropriate
        !! partial rows in the lower triangle of the CSR matrix.
        do iOrb2 = 0, nOrb2 - 1
          jj = iRow + iOrb2
          ind = csr%rowpnt(jj) + nCols(jj)
          do iOrb1 = 0, nOrb1 - 1
            csr%colind(ind+iOrb1) = iCol + iOrb1
          end do
          nCols(jj) = nCols(jj) + nOrb1
        end do

        !! Correspond the columns of the current block to appropriate 
        !! partial rows in the upper triangle of the CSR matrix.
        if (iAt2f /= iAt1) then
          do iOrb1 = 0, nOrb1 - 1
            jj = iCol + iOrb1
            ind = csr%rowpnt(jj) + nCols(jj)
            do iOrb2 = 0, nOrb2 - 1
              csr%colind(ind+iOrb2) = iRow + iOrb2
            end do
            nCols(jj) = nCols(jj) + nOrb2
          end do
        end if
      end do lpNonZero
    end do lpAt1

    deallocate(zero)
    deallocate(nCols)
    deallocate(iNonZero)
    
  end subroutine createEquivCSR_cplx

  

  !!* Folds the internal sparse formatted matrix to the compressed sparse row
  !!* format (complex version).
  !!* @param csr Resulting CSR matrix.
  !!* @param sparse The internal sparse matrix to fold.
  !!* @param iAtomStart Starting positions of the atoms in the square matrix
  !!* @param iPair Starting position of atom-neighbor interaction in the sparse
  !!*   matrix.
  !!* @param iNeighbor Index of neighbors.
  !!* @param nNeighbor Number of neighbors.
  !!* @param img2CentCell Image of the atoms in the central cell.
  !!* @param orb  Orbitals in the system.
  !!* @note In the resulting CSR format both triangles of the matrix are set.
  !!* @caveat The routine should only applied on CSR matrixes created by the
  !!*   createEquiv_cplx subroutine, since it exploits the ordering of the
  !!*   column indexes.
  subroutine foldToCSR_cplx(csr, sparse, kPoint, iAtomStart, iPair, iNeighbor, &
      &nNeighbor, img2CentCell, iCellVec, cellVec, orb)
    type(z_CSR), intent(inout) :: csr
    real(dp), intent(in) :: sparse(:)
    real(dp), intent(in) :: kPoint(:)
    integer, intent(in) :: iAtomStart(:)
    integer, intent(in) :: iPair(0:,:)
    integer, intent(in) :: iNeighbor(0:,:)
    integer, intent(in) :: nNeighbor(:)
    integer, intent(in) :: img2CentCell(:)
    integer, intent(in) :: iCellVec(:)
    real(dp), intent(in) :: cellVec(:,:)
    type(TOrbitals), intent(in) :: orb

    integer, allocatable :: nCols(:)
    complex(dp), allocatable :: tmpCol(:,:), phases(:)
    integer :: nAtom, nCellVec
    integer :: iOrb1, nOrb1, nOrb2, iAt1, iAt2, iAt2f, iNeigh
    integer :: iRow, iCol, iVal, ind
    integer :: ii, jj, kk

    nAtom = size(orb%nOrbAtom)
    nCellVec = size(cellVec, dim=2)

    !! Create necessary phase factors
    allocate(phases(nCellVec))
    do ii = 1, nCellVec
      phases(ii) = exp(2.0_dp * pi * (0.0_dp, 1.0_dp) &
          &* dot_product(kPoint, cellVec(:,ii)))
    end do

    allocate(tmpCol(csr%nRow, orb%mOrb))   ! One atomic block column.
    lpAt1: do iAt1 = 1, nAtom
      !! Unfold the columns for the current atom from the sparse matrix.
      nOrb1 = orb%nOrbAtom(iAt1)
      tmpCol(:,:) = 0.0_dp
      do iNeigh = 0, nNeighbor(iAt1)
        ind = iPair(iNeigh,iAt1) + 1
        iAt2 = iNeighbor(iNeigh,iAt1)
        iAt2f = img2CentCell(iAt2)
        nOrb2 = orb%nOrbAtom(iAt2f)
        iRow = iAtomStart(iAt2f)
        tmpCol(iRow:iRow+nOrb2-1, 1:nOrb1) = &
            &tmpCol(iRow:iRow+nOrb2-1, 1:nOrb1) &
            &+ phases(iCellVec(iAt2)) &
            &* reshape(sparse(ind:ind+nOrb2*nOrb1-1), (/ nOrb2, nOrb1 /))
      end do

      !! Copy every column into the appropriate row in the upper triangle of
      !! the CSR matrix (the copy is masked by the sparsity structure stored
      !! in the CSR matrix)
      nOrb1 = orb%nOrbAtom(iAt1)
      iCol = iAtomStart(iAt1)
      do iOrb1 = 1, nOrb1
        ii = csr%rowpnt(iCol + iOrb1 - 1)
        jj = csr%rowpnt(iCol + iOrb1) - 1
        do kk=ii,jj
           csr%nzval(kk) = conjg(tmpCol(csr%colind(kk), iOrb1))
        enddo
      end do
    end do lpAt1
    deallocate(tmpCol)
    deallocate(phases)

    !! Fill the lower triangle of the CSR matrix
    allocate(nCols(csr%nRow))
    nCols(:) = 0
    do iRow = 1, csr%nRow
      !! Starting from the tail of the matrix row
      iVal = csr%rowpnt(iRow + 1) - 1
      iCol = csr%colind(iVal)
      do while (iVal >= csr%rowpnt(iRow) .and. iCol > iRow)
        csr%nzval(csr%rowpnt(iCol)+nCols(iCol)) = conjg(csr%nzval(iVal))
        nCols(iCol) = nCols(iCol) + 1
        iVal = iVal - 1
        iCol = csr%colind(iVal)
      end do
    end do
    deallocate(nCols)

  end subroutine foldToCSR_cplx



  !!* Unfolds a matrix from the CSR form into the internal sparse
  !!* representation (real version).
  !!* @param sparse Updated sparse matrix.
  !!* @param csr CSR matrix
  !!* @param kPoint K-point for the unfolding.
  !!* @param kWeight Weight of the K-point for the unfolding.
  !!* @param iAtomStart Starting positions of the atoms in the square matrix
  !!* @param iPair Starting position of atom-neighbor interaction in the sparse
  !!*   matrix.
  !!* @param iNeighbor Index of neighbors
  !!* @param nNeighbor Number of neighbors
  !!* @param img2CentCell Image of the atoms in the central cell.
  !!* @param iCellVec  Index of the cell shifting vector for each atom.
  !!* @param cellVec  Cell shifting vectors (relative coordinates)
  !!* @param orb  Orbitals in the system.
  !!* @note The CSR matrix must be hermitian.
  !!* @note The unfolded matrix is added to the passed sparse matrix.
  subroutine unfoldFromCSR_cplx(sparse, csr, kPoint, kWeight, iAtomStart, &
      &iPair, iNeighbor, nNeighbor, img2CentCell, iCellVec, cellVec, orb)
    real(dp), intent(inout) :: sparse(:)
    type(z_CSR), intent(in) :: csr
    real(dp), intent(in) :: kPoint(:)
    real(dp), intent(in) :: kWeight
    integer, intent(in) :: iAtomStart(:)
    integer, intent(in) :: iPair(0:,:)
    integer, intent(in) :: iNeighbor(0:,:)
    integer, intent(in) :: nNeighbor(:)
    integer, intent(in) :: img2CentCell(:)
    integer, intent(in) :: iCellVec(:)
    real(dp), intent(in) :: cellVec(:,:)
    type(TOrbitals), intent(in) :: orb

    complex(dp), allocatable :: tmpCol(:,:), phases(:)
    integer :: nAtom, nCellVec
    integer :: iOrb1, nOrb1, nOrb2, iAt1, iAt2, iAt2f, iNeigh
    integer :: iRow, iCol, ind
    integer :: ii

    nAtom = size(orb%nOrbAtom)

    @:ASSERT(csr%nRow == iAtomStart(nAtom+1) - 1)
    @:ASSERT(size(kPoint) == 3)

    allocate(tmpCol(csr%nRow, orb%mOrb))

    !! Create phase factors
    nCellVec = size(cellVec, dim=2)
    allocate(phases(nCellVec))
    do ii = 1, nCellVec
      phases(ii) = exp(-2.0_dp * pi * (0.0_dp, 1.0_dp) &
          &* dot_product(kPoint, cellVec(:,ii)))
    end do

    do iAt1 = 1, nAtom
      !! Put the rows belonging to a certain atom into the appropriate column
      !! of the block column. (Matrix is assumed to be hermitian!)
      tmpCol(:,:) = (0.0_dp, 0.0_dp)
      nOrb1 = orb%nOrbAtom(iAt1)
      iCol = iAtomStart(iAt1)
      do iOrb1 = 1, nOrb1
        ii = iCol + iOrb1 - 1
        do ind = csr%rowpnt(ii), csr%rowpnt(ii+1) - 1
          tmpCol(csr%colind(ind), iOrb1) = conjg(csr%nzval(ind))
          ! NOTE: why conjg ?? 
        end do
      end do

      !! Unfold the block column into the internal sparse format
      do iNeigh = 0, nNeighbor(iAt1)
        ind = iPair(iNeigh,iAt1) + 1
        iAt2 = iNeighbor(iNeigh, iAt1)
        iAt2f = img2CentCell(iAt2)
        nOrb2 = orb%nOrbAtom(iAt2f)
        iRow = iAtomStart(iAt2f)
        sparse(ind:ind+nOrb2*nOrb1-1) = &
            &sparse(ind:ind+nOrb2*nOrb1-1) &
            &+ kWeight&
            &* real(  phases(iCellVec(iAt2)) &
            &* reshape( tmpCol(iRow:iRow+nOrb2-1, 1:nOrb1),(/nOrb2*nOrb1/) ), dp ) 
            
      end do
    end do

    deallocate(tmpCol)
    deallocate(phases)

  end subroutine unfoldFromCSR_cplx

  

  !!* Creates a CSR matrix by cloning the sparsity structure of an other CSR
  !!* matrix (real version).
  !!* @param csrOut  New CSR matrix on exit.
  !!* @param csrIn  CSR matrix with the sparsity map to be cloned.
  !!* @note The elements of the original matrix are not copied.
  subroutine cloneSparsityMap_real_real(csrOut, csrIn)
    type(r_CSR), intent(inout) :: csrOut
    type(r_CSR), intent(in) :: csrIn

    call create(csrOut, csrIn%nRow, csrIn%nCol, csrIn%nnz)
    csrOut%colind = csrIn%colind
    csrOut%rowpnt = csrIn%rowpnt
    csrOut%sorted = csrIn%sorted

  end subroutine cloneSparsityMap_real_real
  

  
  !!* Creates a CSR matrix by cloning the sparsity structure of an other CSR
  !!* matrix (complex version).
  !!* @param csrOut  New CSR matrix on exit.
  !!* @param csrIn  CSR matrix with the sparsity map to be cloned.
  !!* @note The elements of the original matrix are not copied.
  subroutine cloneSparsityMap_cplx_cplx(csrOut, csrIn)
    type(z_CSR), intent(inout) :: csrOut
    type(z_CSR), intent(in) :: csrIn

    call create(csrOut, csrIn%nRow, csrIn%nCol, csrIn%nnz)
    csrOut%colind = csrIn%colind
    csrOut%rowpnt = csrIn%rowpnt
    csrOut%sorted = csrIn%sorted

  end subroutine cloneSparsityMap_cplx_cplx



  !!* Creates a CSR matrix by cloning the sparsity structure of an other CSR
  !!* matrix (complex version).
  !!* @param csrOut  New CSR matrix on exit.
  !!* @param csrIn  CSR matrix with the sparsity map to be cloned.
  !!* @note The elements of the original matrix are not copied.
  subroutine cloneSparsityMap_real_cplx(csrOut, csrIn)
    type(r_CSR), intent(inout) :: csrOut
    type(z_CSR), intent(in) :: csrIn

    call create(csrOut, csrIn%nRow, csrIn%nCol, csrIn%nnz)
    csrOut%colind = csrIn%colind
    csrOut%rowpnt = csrIn%rowpnt
    csrOut%sorted = csrIn%sorted
    
  end subroutine cloneSparsityMap_real_cplx


  
  !!* Creates a CSR matrix by cloning the sparsity structure of an other CSR
  !!* matrix (complex version).
  !!* @param csrOut  New CSR matrix on exit.
  !!* @param csrIn  CSR matrix with the sparsity map to be cloned.
  !!* @note The elements of the original matrix are not copied.
  subroutine cloneSparsityMap_cplx_real(csrOut, csrIn)
    type(z_CSR), intent(inout) :: csrOut
    type(r_CSR), intent(in) :: csrIn

    call create(csrOut, csrIn%nRow, csrIn%nCol, csrIn%nnz)
    csrOut%colind = csrIn%colind
    csrOut%rowpnt = csrIn%rowpnt
    csrOut%sorted = csrIn%sorted
    
  end subroutine cloneSparsityMap_cplx_real



  !!* Creates an empty CSR matrix containing zero elements (real version).
  !!* @param csr  Empty CSR matrix on exit.
  subroutine createEmptyCSR_real(csr)
    type(r_CSR), intent(out) :: csr

    call create(csr, 0, 0, 0)
    csr%nnz = 0
    csr%nRow = 0
    csr%nCol = 0
    csr%sorted = .false. 

  end subroutine createEmptyCSR_real

  

  !!* Creates an empty CSR matrix containing zero elements (complex version).
  !!* @param csr  Empty CSR matrix on exit.
  subroutine createEmptyCSR_cplx(csr)
    type(z_CSR), intent(out) :: csr

    call create(csr, 0, 0, 0)
    csr%nnz = 0
    csr%nRow = 0
    csr%nCol = 0
    csr%sorted = .false. 

  end subroutine createEmptyCSR_cplx

  subroutine rdestroy_CSR(mat)
    type(r_CSR) :: mat
    !print*,'destruct location',%LOC(mat%nzval)
    if (allocated(mat%nzval)) call destroy(mat)
  end subroutine rdestroy_CSR

  subroutine rdestroy_DNS(mat)
    type(r_DNS) :: mat
    if (allocated(mat%val)) call destroy(mat)
  end subroutine rdestroy_DNS

  subroutine zdestroy_CSR(mat)
    type(z_CSR) :: mat
    if (allocated(mat%nzval)) call destroy(mat)
  end subroutine zdestroy_CSR

  subroutine zdestroy_DNS(mat)
    type(z_DNS) :: mat
    if (allocated(mat%val)) call destroy(mat)
  end subroutine zdestroy_DNS


end module mat_conv
