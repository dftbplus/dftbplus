!--------------------------------------------------------------------------------------------------!
!  DFTB+: general package for performing fast atomistic simulations                                !
!  Copyright (C) 2006 - 2023  DFTB+ developers group                                               !
!                                                                                                  !
!  See the LICENSE file for terms of usage and distribution.                                       !
!--------------------------------------------------------------------------------------------------!


#:include "common.fypp"


!> Contains routines for folding internal sparse matrixes into Compressed
!> Sparse Row Format and vice-versa.
!>
!> Caveat All folding and unfolding routines assume, that the CSR matrices have the same sparsity
!> structure as the internal sparse matrices. The CSR matrices must be created using the
!> createEquivCSR routines.
module dftbp_transport_matconv
  use dftbp_common_accuracy, only : dp
  use dftbp_common_constants, only : pi
  use dftbp_type_commontypes, only : TOrbitals
  use libnegf, only: r_CSR, z_CSR, r_DNS, z_DNS, create, destroy
  implicit none

  private
  public :: init,   destruct
  public :: foldToCSR, unfoldFromCSR


  !> Creates a Compressed Sparse Row matrix with equivalent sparsity pattern to the internal sparse
  !> matrices.
  interface init
    module procedure createEmptyCSR_real
    module procedure createEmptyCSR_cplx
    module procedure createEquivCSR_real
    module procedure createEquivCSR_cplx
    module procedure createEquivCSR_pauli
    module procedure cloneSparsityMap_real_real
    module procedure cloneSparsityMap_cplx_cplx
    module procedure cloneSparsityMap_real_cplx
    module procedure cloneSparsityMap_cplx_real
  end interface init

  !> Destructor for the different structures.
  interface destruct
    module procedure rdestroy_CSR, zdestroy_CSR
    module procedure rdestroy_DNS, zdestroy_DNS
  end interface destruct

  !> Folds internal sparse storage to Compressed Sparse Row format.
  interface foldToCSR
    module procedure foldToCSR_real
    module procedure foldToCSR_cplx
    module procedure foldToCSR_pauli
  end interface foldToCSR

  !> Unfolds from Compressed Sparse Row format into internal sparse storage.
  interface unfoldFromCSR
    module procedure unfoldFromCSR_real
    module procedure unfoldFromCSR_cplx
    module procedure unfoldFromCSR_pauli
  end interface unfoldFromCSR

contains

  !> Creates a Compressed Sparse Row matrix with a sparsity structure in
  !> accordance with current geometry.
  !>
  !> Note: The subroutine only creates the structure of the CSR matrix, but
  !> does not fill up the nonzero values with anything useful.
  !> The resulting CSR matrix has storage for both triangles of the square matrix.
  subroutine createEquivCSR_real(csr, iAtomStart, iNeighbor, nNeighbor, img2CentCell, orb)

    !> Compressed sparse row matrix on exit
    type(r_CSR), intent(out) :: csr

    !> Starting position of the atoms in the square H/S.
    integer, intent(in) :: iAtomStart(:)

    !> Neighborlist for each atom.
    integer, intent(in) :: iNeighbor(0:,:)

    !> Number of neighbors for each atom.
    integer, intent(in) :: nNeighbor(:)

    !> Folded image in the central cell for each atom.
    integer, intent(in) :: img2CentCell(:)

    !> orb  Orbitals in the system.
    type(TOrbitals), intent(in) :: orb

    integer, allocatable :: nColAtom(:), nCols(:), iNonZero(:)
    integer, allocatable :: rowpnt(:)
    logical, allocatable :: zero(:)
    integer :: nAtom, nNonZero
    integer :: iOrb1, iOrb2, nOrb1, nOrb2, iAt1, iAt2f, iNeigh
    integer :: iRow, iCol, ind, nrow
    integer :: ii, jj

    nAtom = size(orb%nOrbAtom)

    ! Count nr. of nonzero columns in the square (folded) form for each atom.  Holds nr. of nonzero
    ! columns for each atom.
    allocate(nColAtom(nAtom))
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

    ! Calculate CSR row pointers:
    ! A row for a certain orbital of a certain atom is as long, as the
    ! nr. of columns determined previously for the atom.
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

    ! Nr. of CSR columns already filled
    allocate(nCols(csr%nRow))
    nCols(:) = 0
    ! Index of the nonzero blocks
    allocate(iNonZero(nAtom))

    ! Loop over all atoms (over all atomic columns in the rectangular picture)
    lpAt1: do iAt1 = 1, nAtom
      iCol = iAtomStart(iAt1)
      ! Width of the atomic column
      nOrb1 = orb%nOrbAtom(iAt1)

      ! Mark nonzero blocks in the folded matrix for the current atom
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

      ! Initialise CSR internal pointers according the non-zero blocks
      lpNonZero: do ii = 1, nNonZero
        iAt2f = iNonZero(ii)
        nOrb2 = orb%nOrbAtom(iAt2f)
        iRow = iAtomStart(iAt2f)

        ! Correspond rows of the current atomic block column to the appropriate
        ! partial rows in the lower triangle of the CSR matrix.
        do iOrb2 = 0, nOrb2 - 1
          jj = iRow + iOrb2
          ind = csr%rowpnt(jj) + nCols(jj)
          do iOrb1 = 0, nOrb1 - 1
            csr%colind(ind+iOrb1) = iCol + iOrb1
          end do
          nCols(jj) = nCols(jj) + nOrb1
        end do

        ! Correspond the columns of the current block to appropriate
        ! partial rows in the upper triangle of the CSR matrix.
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



  !> Folds the internal sparse formatted matrix to the compressed sparse row
  !> format (real version).
  !>
  !> Note: In the resulting CSR format both triangles of the matrix are set.
  !> Caveat: The routine should only applied on CSR matrixes created by the
  !> createEquiv_real subroutine, since it exploits the ordering of the
  !> column indexes.
  subroutine foldToCSR_real(csr, sparse, iAtomStart, iPair, iNeighbor, nNeighbor, img2CentCell, orb)

    !> csr Resulting CSR matrix
    type(r_CSR), intent(inout) :: csr

    !> sparse The internal sparse matrix to fold.
    real(dp), intent(in) :: sparse(:)

    !> iAtomStart Starting positions of the atoms in the square matrix
    integer, intent(in) :: iAtomStart(:)

    !> iPair Starting position of atom-neighbor interaction in the sparse matrix.
    integer, intent(in) :: iPair(0:,:)

    !> iNeighbor Index of neighbors
    integer, intent(in) :: iNeighbor(0:,:)

    !> nNeighbor Number of neighbors
    integer, intent(in) :: nNeighbor(:)

    !> img2CentCell Image of the atoms in the central cell.
    integer, intent(in) :: img2CentCell(:)

    !> orb  Orbitals in the system.
    type(TOrbitals), intent(in) :: orb

    integer, allocatable :: nCols(:)
    real(dp), allocatable :: tmpCol(:,:)
    integer :: nAtom
    integer :: iOrb1, nOrb1, nOrb2, iAt1, iAt2f, iNeigh
    integer :: iRow, iCol, iVal, ind
    integer :: ii, jj, kk

    nAtom = size(orb%nOrbAtom)

    csr%sorted = .false.

    ! One atomic block column.
    allocate(tmpCol(csr%nRow, orb%mOrb))
    lpAt1: do iAt1 = 1, nAtom
      ! Unfold the columns for the current atom from the sparse matrix.
      nOrb1 = orb%nOrbAtom(iAt1)
      tmpCol(:,:) = 0.0_dp
      do iNeigh = 0, nNeighbor(iAt1)
        ind = iPair(iNeigh,iAt1) + 1
        iAt2f = img2CentCell(iNeighbor(iNeigh,iAt1))
        nOrb2 = orb%nOrbAtom(iAt2f)
        iRow = iAtomStart(iAt2f)
        tmpCol(iRow:iRow+nOrb2-1, 1:nOrb1) = tmpCol(iRow:iRow+nOrb2-1, 1:nOrb1)&
            & + reshape(sparse(ind:ind+nOrb2*nOrb1-1), (/ nOrb2, nOrb1 /))
      end do

      ! Copy every column into the appropriate row in the upper triangle of
      ! the CSR matrix (the copy is masked by the sparsity structure stored
      ! in the CSR matrix)
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

    ! Fill the lower triangle of the CSR matrix
    allocate(nCols(csr%nRow))
    nCols(:) = 0
    do iRow = 1, csr%nRow
      ! Starting from the tail of the matrix row
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



  !> Unfolds a matrix from the CSR form into the internal DFTB+ sparse representation (real
  !> version).
  !>
  !> Note: The CSR matrix must be symmetric. The unfolded matrix is added to the passed sparse
  !> matrix.
  subroutine unfoldFromCSR_real(sparse, csr, iAtomStart, iPair, iNeighbor, nNeighbor, img2CentCell,&
      & orb)

    !> sparse Updated sparse matrix.
    real(dp), intent(inout) :: sparse(:)

    !> CSR matrix
    type(r_CSR), intent(in) :: csr

    !> Starting positions of the atoms in the square matrix
    integer, intent(in) :: iAtomStart(:)

    !> Starting position of atom-neighbor interaction in the sparse matrix.
    integer, intent(in) :: iPair(0:,:)

    !> Index of neighbors
    integer, intent(in) :: iNeighbor(0:,:)

    !> Number of neighbors
    integer, intent(in) :: nNeighbor(:)

    !> Image of the atoms in the central cell.
    integer, intent(in) :: img2CentCell(:)

    !> Orbitals in the system.
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
      ! Put the rows belonging to a certain atom into the appropriate column
      ! of the block column. (Matrix is assumed to be symmetric!)
      tmpCol(:,:) = 0.0_dp
      nOrb1 = orb%nOrbAtom(iAt1)
      iCol = iAtomStart(iAt1)
      do iOrb1 = 1, nOrb1
        ii = iCol + iOrb1 - 1
        do ind = csr%rowpnt(ii), csr%rowpnt(ii+1) - 1
          tmpCol(csr%colind(ind), iOrb1) = real(csr%nzval(ind))
        end do
      end do

      ! Unfold the block column into the internal sparse format
      do iNeigh = 0, nNeighbor(iAt1)
        ind = iPair(iNeigh,iAt1) + 1
        iAt2 = iNeighbor(iNeigh, iAt1)
        iAt2f = img2CentCell(iAt2)
        nOrb2 = orb%nOrbAtom(iAt2f)
        iRow = iAtomStart(iAt2f)
        sparse(ind:ind+nOrb2*nOrb1-1) = sparse(ind:ind+nOrb2*nOrb1-1)&
            & + reshape(tmpCol(iRow:iRow+nOrb2-1, 1:nOrb1), (/nOrb2*nOrb1/))
      end do
    end do

    deallocate(tmpCol)

  end subroutine unfoldFromCSR_real



  !> Creates a Compressed Sparse Row matrix with a sparsity structure in accordance with current
  !> geometry.
  !>
  !> Note: The subroutine only creates the structure of the CSR matrix, but does not fill up the
  !> nonzero values with anything useful. The resulting CSR matrix has storage for both triangles
  !> of the square matrix.
  subroutine createEquivCSR_cplx(csr, iAtomStart, iNeighbor, nNeighbor, img2CentCell, orb)

    !> Compressed sparse row matrix on exit
    type(z_CSR), intent(out) :: csr

    !> Starting position of the atoms in the square H/S.
    integer, intent(in) :: iAtomStart(:)

    !> Neighborlist for each atom.
    integer, intent(in) :: iNeighbor(0:,:)

    !> Number of neighbors for each atom.
    integer, intent(in) :: nNeighbor(:)

    !> Folded image in the central cell for each atom.
    integer, intent(in) :: img2CentCell(:)

    !> Orbitals in the system.
    type(TOrbitals), intent(in) :: orb

    integer, allocatable :: nColAtom(:), nCols(:), iNonZero(:)
    integer, allocatable :: rowpnt(:)
    logical, allocatable :: zero(:)
    integer :: nAtom, nNonZero
    integer :: iOrb1, iOrb2, nOrb1, nOrb2, iAt1, iAt2f, iNeigh
    integer :: iRow, iCol, ind, nrow
    integer :: ii, jj

    nAtom = size(orb%nOrbAtom)

    ! Count nr. of nonzero columns in the square (folded) form for each atom
    ! Nr. of nonzero columns for each atom
    allocate(nColAtom(nAtom))
    allocate(zero(nAtom))
    nColAtom(:) = 0
    do iAt1 = 1, nAtom
      zero(:) = .true.
      nOrb1 = orb%nOrbAtom(iAt1)
      do iNeigh = 0, nNeighbor(iAt1)
        iAt2f = img2CentCell(iNeighbor(iNeigh,iAt1))
        ! Do not count twice atoms in the central cell
        if (zero(iAt2f)) then
          nColAtom(iAt1) = nColAtom(iAt1) + orb%nOrbAtom(iAt2f)
          zero(iAt2f) = .false.
          ! We need this count because the nMap is not symmetric
          if (iAt2f /= iAt1) then
            nColAtom(iAt2f) = nColAtom(iAt2f) + nOrb1
          end if
        end if
      end do
    end do

    nrow = iAtomStart(nAtom+1) - 1
    allocate(rowpnt(nrow+1) )
    ! Calculate CSR row pointers:
    ! A row for a certain orbital of a certain atom is as long, as the
    ! nr. of columns determined previously for the atom.
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

    ! Nr. of CSR columns already filled
    allocate(nCols(csr%nRow))
    nCols(:) = 0
    ! Index of the nonzero blocks
    allocate(iNonZero(nAtom))

    ! Loop over all atoms (over all atomic columns in the rectangular picture)
    lpAt1: do iAt1 = 1, nAtom
      iCol = iAtomStart(iAt1)
      nOrb1 = orb%nOrbAtom(iAt1)           ! Width of the atomic column

      ! Mark nonzero blocks in the folded matrix for the current atom
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

      ! Initialise CSR internal pointers according the non-zero blocks
      lpNonZero: do ii = 1, nNonZero
        iAt2f = iNonZero(ii)
        nOrb2 = orb%nOrbAtom(iAt2f)
        iRow = iAtomStart(iAt2f)

        ! Correspond rows of the current atomic block column to the appropriate
        ! partial rows in the lower triangle of the CSR matrix.
        do iOrb2 = 0, nOrb2 - 1
          jj = iRow + iOrb2
          ind = csr%rowpnt(jj) + nCols(jj)
          do iOrb1 = 0, nOrb1 - 1
            csr%colind(ind+iOrb1) = iCol + iOrb1
          end do
          nCols(jj) = nCols(jj) + nOrb1
        end do

        ! Correspond the columns of the current block to appropriate
        ! partial rows in the upper triangle of the CSR matrix.
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


  !> Creates a Compressed Sparse Row matrix with a sparsity structure in accordance with current
  !> geometry.
  !>
  !> Note: The subroutine only creates the structure of the CSR matrix, but does not fill up the
  !> nonzero values with anything usefull. The resulting CSR matrix has storage for both triangles
  !> of the square matrix.
  subroutine createEquivCSR_pauli(csr, iAtomStart, iNeighbor, nNeighbor, &
      &img2CentCell, orb, tDualSpinOrbit)

    !> Resulting CSR matrix.
    type(z_CSR), intent(out) :: csr

    !> Starting position of the atoms in the square H/S.
    integer, intent(in) :: iAtomStart(:)

    !> Neighborlist for each atom.
    integer, intent(in) :: iNeighbor(0:,:)

    !> Number of neighbors for each atom.
    integer, intent(in) :: nNeighbor(:)

    !> Folded image in the central cell for each atom.
    integer, intent(in) :: img2CentCell(:)

    !> orbital definition structure
    type(TOrbitals), intent(in) :: orb

    !> enables dual spin-orbit coupling
    logical, intent(in) :: tDualSpinOrbit


    integer, allocatable :: nColAtom(:), nCols(:), iNonZero(:)
    integer, allocatable :: iAtomStart2(:), rowpnt(:)
    logical, allocatable :: zero(:)
    integer :: nAtom, nNonZero
    integer :: iOrb1, iOrb2, nOrb1, nOrb2, iAt1, iAt2f, iNeigh
    integer :: iRow, iCol, ind, nrow
    integer :: ii, jj


    nAtom = size(orb%nOrbAtom)

    ! Count nr. of nonzero columns in the square (folded) form for each atom
    ! Nr. of nonzero columns for each atom.
    allocate(nColAtom(nAtom))
    allocate(zero(nAtom))
    allocate(iAtomStart2(nAtom+1))

    nColAtom(:) = 0
    do iAt1 = 1, nAtom
      zero(:) = .true.
      ! Compute the atom indexes of the 2x2 Hamiltonian
      iAtomStart2(iAt1) = 2*iAtomStart(iAt1) - 1
      nOrb1 = orb%nOrbAtom(iAt1)
      do iNeigh = 0, nNeighbor(iAt1)
        iAt2f = img2CentCell(iNeighbor(iNeigh,iAt1))
        ! Do not count twice atoms in the central cell
        if (zero(iAt2f)) then
          nColAtom(iAt1) = nColAtom(iAt1) + 2*orb%nOrbAtom(iAt2f)
          zero(iAt2f) = .false.
          ! We need this count because the nMap is not symmetric
          if (iAt2f /= iAt1) then
            nColAtom(iAt2f) = nColAtom(iAt2f) + 2*nOrb1
          end if
        end if
      end do
    end do
    iAtomStart2(nAtom+1) = 2*iAtomStart(nAtom+1) - 1

    nrow = iAtomStart2(nAtom+1) - 1
    allocate(rowpnt(nrow+1) )
    ! Calculate CSR row pointers:
    ! A row for a certain orbital of a certain atom is as long, as the
    ! nr. of columns determined previously for the atom.
    rowpnt(1) = 1
    ind = 1
    do iAt1 = 1, nAtom
      nOrb1 = 2*orb%nOrbAtom(iAt1)
      do iOrb1 = 1, nOrb1
        rowpnt(ind+iOrb1) = rowpnt(ind) + iOrb1 * nColAtom(iAt1)
      end do
      ind = ind + nOrb1
    end do
    deallocate(nColAtom)

    call create(csr, nrow, nrow, rowpnt(nrow + 1) - 1)
    csr%rowpnt = rowpnt
    deallocate(rowpnt)

    ! Nr. of CSR columns already filled
    allocate(nCols(csr%nRow))
    nCols(:) = 0
    ! Index of the nonzero blocks
    allocate(iNonZero(nAtom))

    ! Loop over all atoms (over all atomic columns in the rectangular picture)
    lpAt1: do iAt1 = 1, nAtom
      iCol = iAtomStart2(iAt1)
      nOrb1 = 2*orb%nOrbAtom(iAt1)         ! Width of the atomic column

      ! Mark nonzero blocks in the folded matrix for the current atom
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

      ! Initialise CSR internal pointers according the non-zero blocks
      lpNonZero: do ii = 1, nNonZero
        iAt2f = iNonZero(ii)
        nOrb2 = 2*orb%nOrbAtom(iAt2f)
        iRow = iAtomStart2(iAt2f)

        ! Correspond rows of the current atomic block column to the appropriate
        ! partial rows in the lower triangle of the CSR matrix.
        do iOrb2 = 0, nOrb2 - 1
          jj = iRow + iOrb2
          ind = csr%rowpnt(jj) + nCols(jj)
          do iOrb1 = 0, nOrb1 - 1
            csr%colind(ind+iOrb1) = iCol + iOrb1
          end do
          nCols(jj) = nCols(jj) + nOrb1
        end do

        ! Correspond the columns of the current block to appropriate
        ! partial rows in the upper triangle of the CSR matrix.
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

  end subroutine createEquivCSR_pauli


  !> Folds the internal sparse formatted matrix to the compressed sparse row format (complex
  !> version).
  !> It allows for 1x1 or 2x2 spin matrices, resulting useful to convert the overlap matrix
  !>
  !> Note: In the resulting CSR format both triangles of the matrix are set.
  !> Caveat The routine should only applied on CSR matrixes created by the createEquiv_cplx
  !> subroutine, since it exploits the ordering of the column indexes.
  subroutine foldToCSR_cplx(csr, sparse, kPoint, iAtomStart, iPair, iNeighbor, nNeighbor,&
      & img2CentCell, iCellVec, cellVec, orb, nComp)

    !> Resulting CSR matrix.
    type(z_CSR), intent(inout) :: csr

    !> The internal sparse matrix to fold.
    real(dp), intent(in) :: sparse(:)

    !> Current k-point
    real(dp), intent(in) :: kPoint(:)

    !> Starting positions of the atoms in the square matrix
    integer, intent(in) :: iAtomStart(:)

    !> iPair Starting position of atom-neighbor interaction in the sparse matrix.
    integer, intent(in) :: iPair(0:,:)

    !> iNeighbor Index of neighbors.
    integer, intent(in) :: iNeighbor(0:,:)

    !> nNeighbor Number of neighbors.
    integer, intent(in) :: nNeighbor(:)

    !> img2CentCell Image of the atoms in the central cell.
    integer, intent(in) :: img2CentCell(:)

    !> Index for cells atomic images are sited in
    integer, intent(in) :: iCellVec(:)

    !> vectors to crystal unit cells
    real(dp), intent(in) :: cellVec(:,:)

    !> orb  Orbitals in the system.
    type(TOrbitals), intent(in) :: orb

    !> number of Spin components (1 or 2)
    integer, intent(in) :: nComp


    integer, allocatable :: nCols(:), iAtomStart2(:)
    complex(dp), allocatable :: tmpCol(:,:), phases(:), tmpMat(:,:)
    integer :: nAtom, nCellVec, shift
    integer :: iOrb1, nOrb1, nOrb2, iAt1, iAt2, iAt2f, iNeigh
    integer :: iRow, iCol, iVal, ind
    integer :: ii, jj, kk

    if (nComp == 2) then
      shift = 1
    else
      shift = 0
    end if

    nAtom = size(orb%nOrbAtom)
    nCellVec = size(cellVec, dim=2)

    ! Create necessary phase factors
    allocate(phases(nCellVec))
    do ii = 1, nCellVec
      phases(ii) = exp(2.0_dp * pi * (0.0_dp, 1.0_dp) * dot_product(kPoint, cellVec(:,ii)))
    end do

    iAtomStart2 = nComp*iAtomStart - shift

    allocate(tmpCol(csr%nRow, nComp*orb%mOrb))   ! One atomic block column
    allocate(tmpMat(nComp*orb%mOrb, nComp*orb%mOrb))
    csr%nzval = (0.0_dp, 0.0_dp)

    lpAt1: do iAt1 = 1, nAtom
      ! Unfold the columns for the current atom from the sparse matrix.
      nOrb1 = nComp * orb%nOrbAtom(iAt1)
      tmpCol(:,:) = 0.0_dp
      do iNeigh = 0, nNeighbor(iAt1)
        ind = iPair(iNeigh,iAt1) + 1
        iAt2 = iNeighbor(iNeigh,iAt1)
        iAt2f = img2CentCell(iAt2)
        nOrb2 = nComp * orb%nOrbAtom(iAt2f)
        iRow = iAtomStart2(iAt2f)
        call getSparseBlock(sparse, ind, tmpMat, nOrb2/nComp, nOrb1/nComp, nComp)
        tmpCol(iRow:iRow+nOrb2-1, 1:nOrb1) = tmpCol(iRow:iRow+nOrb2-1, 1:nOrb1)&
            & + phases(iCellVec(iAt2)) * tmpMat(1:nOrb2, 1:nOrb1)
      end do

      ! Copy every column into the appropriate row in the upper triangle of
      ! the CSR matrix (the copy is masked by the sparsity structure stored
      ! in the CSR matrix)
      iCol = iAtomStart2(iAt1)
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

    ! Fill the lower triangle of the CSR matrix
    allocate(nCols(csr%nRow))
    nCols(:) = 0
    do iRow = 1, csr%nRow
      ! Starting from the tail of the matrix row
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


  !> Folds the internal sparse matrix to the compressed sparse row format
  !> (Spin up-down representation from Pauli matrices).
  !>
  !> Note: In the resulting CSR format both triangles of the matrix are set.
  !> Caveat The routine should only applied on CSR matrixes created by the createEquiv_cplx
  !> subroutine, since it exploits the ordering of the column indexes.
  subroutine foldToCSR_Pauli(csr, sparse, kPoint, iAtomStart, iPair, iNeighbor, nNeighbor,&
      & img2CentCell, iCellVec, cellVec, orb, tSpinOrbit, iHam)

    !> Resulting CSR matrix.
    type(z_CSR), intent(inout) :: csr

    !> The internal sparse matrix to fold (4 blocks)
    real(dp), intent(in) :: sparse(:,:)

    !> Current k-point
    real(dp), intent(in) :: kPoint(:)

    !> Starting positions of the atoms in the square matrix
    integer, intent(in) :: iAtomStart(:)

    !> iPair Starting position of atom-neighbor interaction in the sparse matrix.
    integer, intent(in) :: iPair(0:,:)

    !> iNeighbor Index of neighbors.
    integer, intent(in) :: iNeighbor(0:,:)

    !> nNeighbor Number of neighbors.
    integer, intent(in) :: nNeighbor(:)

    !> img2CentCell Image of the atoms in the central cell.
    integer, intent(in) :: img2CentCell(:)

    !> Index for cells atomic images are sited in
    integer, intent(in) :: iCellVec(:)

    !> vectors to crystal unit cells
    real(dp), intent(in) :: cellVec(:,:)

    !> orb  Orbitals in the system.
    type(TOrbitals), intent(in) :: orb

    !> whether SOC is active
    logical, intent(in) :: tSpinOrbit

    !> imaginary part of sparse hamiltonian
    real(dp), intent(in), allocatable :: iHam(:, :)


    integer, allocatable :: nCols(:), iAtomStart2(:)
    complex(dp), allocatable :: tmpCol(:,:), tmpMat(:,:), phases(:)
    integer :: nAtom, nCellVec
    integer :: iOrb1, nOrb1, nOrb2, iAt1, iAt2, iAt2f, iNeigh
    integer :: iRow, iCol, iVal, ind
    integer :: ii, jj, kk

    nAtom = size(orb%nOrbAtom)
    nCellVec = size(cellVec, dim=2)

    ! Create necessary phase factors
    allocate(phases(nCellVec))
    do ii = 1, nCellVec
      phases(ii) = exp(2.0_dp * pi * (0.0_dp, 1.0_dp) * dot_product(kPoint, cellVec(:,ii)))
    end do

    iAtomStart2 = 2*iAtomStart - 1

    allocate(tmpCol(csr%nRow, 2*orb%mOrb))   ! One atomic block column.
    allocate(tmpMat(2*orb%mOrb, 2*orb%mOrb))

    lpAt1: do iAt1 = 1, nAtom
      ! Unfold the columns for the current atom from the sparse matrix.
      nOrb1 = 2*orb%nOrbAtom(iAt1)
      tmpCol(:,:) = 0.0_dp
      do iNeigh = 0, nNeighbor(iAt1)
        ind = iPair(iNeigh,iAt1) + 1
        iAt2 = iNeighbor(iNeigh,iAt1)
        iAt2f = img2CentCell(iAt2)
        nOrb2 = 2*orb%nOrbAtom(iAt2f)
        iRow = iAtomStart2(iAt2f)
        call getPauliBlock(sparse, ind, tmpMat, nOrb2/2, nOrb1/2, tSpinOrbit, iHam)
        tmpCol(iRow:iRow+nOrb2-1, 1:nOrb1) = tmpCol(iRow:iRow+nOrb2-1, 1:nOrb1)&
            & + phases(iCellVec(iAt2)) * tmpMat(1:nOrb2, 1:nOrb1)
      end do

      ! Copy every column into the appropriate row in the upper triangle of
      ! the CSR matrix (the copy is masked by the sparsity structure stored
      ! in the CSR matrix)
      iCol = iAtomStart2(iAt1)
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

    ! Fill the lower triangle of the CSR matrix
    allocate(nCols(csr%nRow))
    nCols(:) = 0
    do iRow = 1, csr%nRow
      ! Starting from the tail of the matrix row
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

  end subroutine foldToCSR_Pauli

  subroutine getSparseBlock(sparse, ind, tmpMat, nOrb2, nOrb1, nComp)

    !> sparse hamiltonian
    real(dp), intent(in) :: sparse(:)

    !> start index of atom block in sparse matrix
    integer, intent(in) :: ind

    !> dense atom block
    complex(dp), intent(out) :: tmpMat(:, :)

    !> number of orbitals along rows
    integer, intent(in) :: nOrb2

    !> number of orbitals along cols
    integer, intent(in) :: nOrb1

    !> number of Spin components (1 or 2)
    integer, intent(in) :: nComp


    real(dp), allocatable :: work(:, :)

    tmpMat = (0.0_dp, 0.0_dp)
    work = reshape(sparse(ind:ind+nOrb2*nOrb1-1), [nOrb2, nOrb1])
    tmpMat(1:nOrb2, 1:nOrb1) = work(1:nOrb2, 1:nOrb1)

    ! 2x2 spin matrix
    if (nComp == 2) then
      tmpMat(nOrb2+1:2*nOrb2, nOrb1+1:2*nOrb1) = work(1:nOrb2, 1:nOrb1)
    end if

  end subroutine getSparseBlock



  subroutine getPauliBlock(ham, ind, tmpMat, nOrb2, nOrb1, tSpinOrbit, iHam)

    !> sparse hamiltonian
    real(dp), intent(in) :: ham(:, :)

    !> start index of atom block in sparse matrix
    integer, intent(in) :: ind

    !> dense atom block
    complex(dp), intent(out) :: tmpMat(:, :)

    !> number of orbitals along rows
    integer, intent(in) :: nOrb2

    !> number of orbitals along cols
    integer, intent(in) :: nOrb1

    !> true if SOC is active
    logical, intent(in) :: tSpinOrbit

    !> imaginary part of sparse hamiltonian
    real(dp), intent(in), allocatable :: iHam(:, :)


    real(dp), allocatable :: work(:, :)
    integer :: n2, n1p1, n1t2, n2p1, n2t2

    n2 = nOrb2*nOrb1
    n1p1 = nOrb1+1
    n2p1 = nOrb2+1
    n1t2 = 2*nOrb1
    n2t2 = 2*nOrb2

    tmpMat = (0.0_dp, 0.0_dp)
    ! 1 0 charge part
    ! 0 1
    work = reshape(ham(ind:ind+n2-1, 1), [nOrb2, nOrb1])
    tmpMat(1:nOrb2, 1:nOrb1) = 0.5_dp*work(1:nOrb2, 1:nOrb1)
    tmpMat(n2p1:n2t2, n1p1:n1t2) = 0.5_dp*work(1:nOrb2, 1:nOrb1)

    if (allocated(iHam)) then
      work = reshape(iHam(ind:ind+n2-1, 1), [nOrb2, nOrb1])
      tmpMat(1:nOrb2, 1:nOrb1) = tmpMat(1:nOrb2,1:nOrb1) &
          &+ 0.5_dp*cmplx(0,1,dp) * work(1:nOrb2, 1:nOrb1)
      tmpMat(n2p1:n2t2, n1p1:n1t2) = tmpMat(n2p1:n2t2, n1p1:n1t2) &
          &+ 0.5_dp*cmplx(0,1,dp) * work(1:nOrb2, 1:nOrb1)
    end if

    if (.not.tSpinOrbit) return

    ! 0 1 x part
    ! 1 0
    work = reshape(ham(ind:ind+n2-1, 2), [nOrb2, nOrb1])
    tmpMat(n2p1:n2t2, 1:nOrb1) = tmpMat(n2p1:n2t2, 1:nOrb1) &
         & + 0.5_dp * work(1:nOrb2, 1:nOrb1)
    tmpMat(1:nOrb2, n1p1:n1t2) = tmpMat(1:nOrb2, n1p1:n1t2) &
         & + 0.5_dp * work(1:nOrb2, 1:nOrb1)

    ! 0 -i y part
    ! i  0
    work = reshape(ham(ind:ind+nOrb2*nOrb1-1, 3), [nOrb2, nOrb1])
    tmpMat(n2p1:n2t2, 1:nOrb1) = tmpMat(n2p1:n2t2, 1:nOrb1) &
         & + 0.5_dp * cmplx(0,1,dp) * work(1:nOrb2, 1:nOrb1)
    tmpMat(1:nOrb2, n1p1:n1t2) = tmpMat(1:nOrb2, n1p1:n1t2) &
         & - 0.5_dp * cmplx(0,1,dp) * work(1:nOrb2, 1:nOrb1)

    ! 1  0 z part
    ! 0 -1
    work = reshape(ham(ind:ind+nOrb2*nOrb1-1, 4), [nOrb2, nOrb1])
    tmpMat(1:nOrb2, 1:nOrb1) = tmpMat(1:nOrb2, 1:nOrb1) &
         & + 0.5_dp * work(1:nOrb2, 1:nOrb1)
    tmpMat(n2p1:n2t2, 1:nOrb1) = tmpMat(n2p1:n2t2, 1:nOrb1) &
         & - 0.5_dp * work(1:nOrb2, 1:nOrb1)

    ! iHam(:,2:4) = (Im{L+}, -Re{L+}, Im{Lz})
    ! iHam(:,2:4) = (Im{L-}, +Re{L-}, Im{Lz})
    ! => L+ = i* iHam(:,2) - iHam(:,3)
    ! => L- = i* iHam(:,2) + iHam(:,3)
    ! (Lz L-)
    ! (L+ Lz)
    if (allocated(iHam)) then
      work = reshape(iHam(ind:ind+n2-1, 2), [nOrb2, nOrb1])
      tmpMat(n2p1:n2t2, 1:nOrb1) = tmpMat(n2p1:n2t2, 1:nOrb1) &
          & + 0.5_dp*cmplx(0,1,dp) * work(1:nOrb2, 1:nOrb1)
      tmpMat(1:nOrb2, n1p1:n1t2) = tmpMat(1:nOrb2, n1p1:n1t2) &
          & + 0.5_dp*cmplx(0,1,dp) * work(1:nOrb2, 1:nOrb1)

      work = reshape(iHam(ind:ind+n2-1, 3), [nOrb2, nOrb1])
      tmpMat(n2p1:n2t2, 1:nOrb1) = tmpMat(n2p1:n2t2, 1:nOrb1) &
          & - 0.5_dp * work(1:nOrb2, 1:nOrb1)
      tmpMat(1:nOrb2, n1p1:n1t2) = tmpMat(1:nOrb2, n1p1:n1t2) &
          & + 0.5_dp * work(1:nOrb2, 1:nOrb1)

      work = reshape(iHam(ind:ind+n2-1, 4), [nOrb2, nOrb1])
      tmpMat(1:nOrb2, 1:nOrb1) = tmpMat(1:nOrb2,1:nOrb1) &
          &+ 0.5_dp*cmplx(0,1,dp) * work(1:nOrb2, 1:nOrb1)
      tmpMat(n2p1:n2t2, n1p1:n1t2) = tmpMat(n2p1:n2t2, n1p1:n1t2) &
          &- 0.5_dp*cmplx(0,1,dp) * work(1:nOrb2, 1:nOrb1)
    end if

  end subroutine getPauliBlock


  !> Unfolds a matrix from the CSR form into the internal sparse representation (real version).
  !> Allows for 1 or 2 spin components, resulting useful to conver the ERho matrix
  !>
  !> Note: The CSR matrix must be hermitian. The unfolded matrix is added to the passed sparse
  !> matrix.
  subroutine unfoldFromCSR_cplx(sparse, csr, kPoint, kWeight, iAtomStart, iPair, iNeighbor,&
      & nNeighbor, img2CentCell, iCellVec, cellVec, orb, nComp)

    !> sparse Updated sparse matrix.
    real(dp), intent(inout) :: sparse(:)

    !> csr CSR matrix
    type(z_CSR), intent(in) :: csr

    !> kPoint K-point for the unfolding.
    real(dp), intent(in) :: kPoint(:)

    !> kWeight Weight of the K-point for the unfolding.
    real(dp), intent(in) :: kWeight

    !> iAtomStart Starting positions of the atoms in the square matrix
    integer, intent(in) :: iAtomStart(:)

    !> iPair Starting position of atom-neighbor interaction in the sparse matrix.
    integer, intent(in) :: iPair(0:,:)

    !> iNeighbor Index of neighbors
    integer, intent(in) :: iNeighbor(0:,:)

    !> nNeighbor Number of neighbors
    integer, intent(in) :: nNeighbor(:)

    !> img2CentCell Image of the atoms in the central cell.
    integer, intent(in) :: img2CentCell(:)

    !> iCellVec  Index of the cell shifting vector for each atom.
    integer, intent(in) :: iCellVec(:)

    !> cellVec  Cell shifting vectors (relative coordinates)
    real(dp), intent(in) :: cellVec(:,:)

    !> orb  Orbitals in the system.
    type(TOrbitals), intent(in) :: orb

    !> number of Spin components (1 or 2)
    integer, intent(in) :: nComp

    complex(dp), allocatable :: tmpCol(:,:), phases(:)
    integer, allocatable :: iAtomStart2(:)
    integer :: nAtom, nCellVec, shift
    integer :: iOrb1, nOrb1, nOrb2, iAt1, iAt2, iAt2f, iNeigh
    integer :: iRow, iCol, ind
    integer :: ii, n12

    if (nComp == 2) then
      shift = 1
    else
      shift = 0
    end if

    nAtom = size(orb%nOrbAtom)
    iAtomStart2 = nComp*iAtomStart - shift

    @:ASSERT(csr%nRow == iAtomStart2(nAtom+1) - 1)
    @:ASSERT(size(kPoint) == 3)

    allocate(tmpCol(csr%nRow, nComp*orb%mOrb))   ! One atomic block column.

    ! Create phase factors
    nCellVec = size(cellVec, dim=2)
    allocate(phases(nCellVec))
    do ii = 1, nCellVec
      phases(ii) = exp(-2.0_dp * pi * (0.0_dp, 1.0_dp) * dot_product(kPoint, cellVec(:,ii)))
    end do

    do iAt1 = 1, nAtom
      ! Put the rows belonging to a certain atom into the appropriate column
      ! of the block column. (Matrix is assumed to be hermitian!)
      tmpCol(:,:) = (0.0_dp, 0.0_dp)
      nOrb1 = orb%nOrbAtom(iAt1)
      iCol = iAtomStart2(iAt1)
      do iOrb1 = 1, nComp*nOrb1
        ii = iCol + iOrb1 - 1
        do ind = csr%rowpnt(ii), csr%rowpnt(ii+1) - 1
          tmpCol(csr%colind(ind), iOrb1) = conjg(csr%nzval(ind))
          ! conjg is needed to transform rows into columns
        end do
      end do

      ! Unfold the block column into the internal sparse format
      do iNeigh = 0, nNeighbor(iAt1)
        ind = iPair(iNeigh,iAt1) + 1
        iAt2 = iNeighbor(iNeigh, iAt1)
        iAt2f = img2CentCell(iAt2)
        nOrb2 = orb%nOrbAtom(iAt2f)
        iRow = iAtomStart2(iAt2f)
        n12 = nOrb1 * nOrb2
        sparse(ind:ind+n12-1) = sparse(ind:ind+n12-1)&
            & + kWeight * real(  phases(iCellVec(iAt2))&
            & * reshape(tmpCol(iRow:iRow+nOrb2-1, 1:nOrb1), (/n12/) ), dp)
        if (nComp==2) then
          sparse(ind:ind+n12-1) = sparse(ind:ind+n12-1)&
             & + kWeight * real(  phases(iCellVec(iAt2))&
             & * reshape(tmpCol(iRow+nOrb2:iRow+nComp*nOrb2-1, nOrb1+1:nComp*nOrb1), (/n12/) ), dp)
        end if
      end do
    end do

    deallocate(tmpCol)
    deallocate(phases)

  end subroutine unfoldFromCSR_cplx


  !> Unfolds a matrix from the CSR form into the internal sparse representation (real version).
  !>
  !> Note: The CSR matrix must be hermitian. The unfolded matrix is added to the passed sparse
  !> matrix.
  subroutine unfoldFromCSR_Pauli(sparse, csr, kPoint, kWeight, iAtomStart, iPair, iNeighbor,&
      & nNeighbor, img2CentCell, iCellVec, cellVec, orb, isparse)

    !> sparse Updated sparse matrix.
    real(dp), intent(inout) :: sparse(:,:)

    !> csr CSR matrix
    type(z_CSR), intent(in) :: csr

    !> kPoint K-point for the unfolding.
    real(dp), intent(in) :: kPoint(:)

    !> kWeight Weight of the K-point for the unfolding.
    real(dp), intent(in) :: kWeight

    !> iAtomStart Starting positions of the atoms in the square matrix
    integer, intent(in) :: iAtomStart(:)

    !> iPair Starting position of atom-neighbor interaction in the sparse matrix.
    integer, intent(in) :: iPair(0:,:)

    !> iNeighbor Index of neighbors
    integer, intent(in) :: iNeighbor(0:,:)

    !> nNeighbor Number of neighbors
    integer, intent(in) :: nNeighbor(:)

    !> img2CentCell Image of the atoms in the central cell.
    integer, intent(in) :: img2CentCell(:)

    !> iCellVec  Index of the cell shifting vector for each atom.
    integer, intent(in) :: iCellVec(:)

    !> cellVec  Cell shifting vectors (relative coordinates)
    real(dp), intent(in) :: cellVec(:,:)

    !> orb  Orbitals in the system.
    type(TOrbitals), intent(in) :: orb

    !> Imaginary part of sparse matrix.
    real(dp), intent(inout), allocatable :: isparse(:,:)


    complex(dp), allocatable :: tmpCol(:,:), tmpMat(:,:), phases(:)
    integer, allocatable :: iAtomStart2(:)
    integer :: nAtom, nCellVec, shift
    integer :: iOrb1, nOrb1, nOrb2, iAt1, iAt2, iAt2f, iNeigh
    integer :: ii, iRow, iCol, ind
    integer :: n1, n2, n1p1, n1t2, n2p1, n2t2, n12

    nAtom = size(orb%nOrbAtom)
    iAtomStart2 = 2*iAtomStart - 1

    @:ASSERT(size(kPoint) == 3)
    @:ASSERT(csr%nRow == iAtomStart2(nAtom+1) - 1)

    allocate(tmpCol(csr%nRow, 2*orb%mOrb))   ! One atomic block column.
    allocate(tmpMat(2*orb%mOrb, 2*orb%mOrb))

    ! Create phase factors
    nCellVec = size(cellVec, dim=2)
    allocate(phases(nCellVec))
    do ii = 1, nCellVec
      phases(ii) = exp(-2.0_dp * pi * (0.0_dp, 1.0_dp) * dot_product(kPoint, cellVec(:,ii)))
    end do

    do iAt1 = 1, nAtom
      ! Put the rows belonging to a certain atom into the appropriate column
      ! of the block column. (Matrix is assumed to be hermitian!)
      tmpCol(:,:) = (0.0_dp, 0.0_dp)
      nOrb1 = 2*orb%nOrbAtom(iAt1)
      iCol = iAtomStart2(iAt1)
      do iOrb1 = 1, nOrb1
        ii = iCol + iOrb1 - 1
        do ind = csr%rowpnt(ii), csr%rowpnt(ii+1) - 1
          tmpCol(csr%colind(ind), iOrb1) = conjg(csr%nzval(ind))
          ! conjg is needed to transform rows into columns
        end do
      end do

      ! Unfold the block column into the internal sparse format
      do iNeigh = 0, nNeighbor(iAt1)
        ind = iPair(iNeigh,iAt1) + 1
        iAt2 = iNeighbor(iNeigh, iAt1)
        iAt2f = img2CentCell(iAt2)
        nOrb2 = 2*orb%nOrbAtom(iAt2f)
        iRow = iAtomStart2(iAt2f)
        tmpMat(1:nOrb2,1:nOrb1) = kWeight * 0.5_dp * phases(iCellVec(iAt2)) &
              & * tmpCol(iRow:iRow+nOrb2-1, 1:nOrb1)
        n12 = nOrb1*nOrb2/4
        n1=nOrb1/2; n2=nOrb2/2
        n1p1=n1+1; n2p1=n2+1
        n1t2=n1*2; n2t2=n2*2
        !rho_q = Re{up-up + down-down}
        sparse(ind:ind+n12-1, 1) = sparse(ind:ind+n12-1, 1) &
            & + reshape(real(tmpMat(1:n2, 1:n1)), (/n12/)) &
            & + reshape(real(tmpMat(n2p1:n2t2, n1p1:n1t2)), (/n12/))

        !rho_x = Re{down-up}
        sparse(ind:ind+n12-1, 2) = sparse(ind:ind+n12-1, 2) &
            & + reshape(real(tmpMat(n2p1:n2t2, 1:n1)), (/n12/)) &
            & + reshape(real(tmpMat(1:n2, n1p1:n1t2)), (/n12/))

        !rho_y = Im{down-up}
        sparse(ind:ind+n12-1, 3) = sparse(ind:ind+n12-1, 3) &
            & + reshape(aimag(tmpMat(n2p1:n2t2, 1:n1)), (/n12/)) &
            & - reshape(aimag(tmpMat(1:n2, n1p1:n1t2)), (/n12/))

        !rho_z = Re{up-up - down-down}
        sparse(ind:ind+n12-1, 4) = sparse(ind:ind+n12-1, 4) &
            & + reshape(real(tmpMat(1:n2, 1:n1)), (/n12/)) &
            & - reshape(real(tmpMat(n2p1:n2t2, n1p1:n1t2)), (/n12/))

        if (allocated(isparse)) then
          !rho_q = up-up + down-down
          isparse(ind:ind+n12-1, 1) = isparse(ind:ind+n12-1, 1) &
              & + reshape(aimag(tmpMat(1:n2, 1:n1)), (/n12/)) &
              & + reshape(aimag(tmpMat(n2p1:n2t2, n1p1:n1t2)), (/n12/))

          !rho_x = Im{down-up}
          isparse(ind:ind+n12-1, 2) = isparse(ind:ind+n12-1, 2) &
              & + reshape(aimag(tmpMat(n2p1:n2t2, 1:n1)), (/n12/)) &
              & + reshape(aimag(tmpMat(1:n2, n1p1:n1t2)), (/n12/))

          !rho_y = -Re{down-up}
          isparse(ind:ind+n12-1, 3) = isparse(ind:ind+n12-1, 3) &
              & - reshape(real(tmpMat(n2p1:n2t2, 1:n1)), (/n12/)) &
              & + reshape(real(tmpMat(1:n2, n1p1:n1t2)), (/n12/))

          !rho_z = Im{up-up - down-down}
          isparse(ind:ind+n12-1, 4) = isparse(ind:ind+n12-1, 4) &
              & + reshape(aimag(tmpMat(1:n2, 1:n1)), (/n12/)) &
              & - reshape(aimag(tmpMat(n2p1:n2t2, n1p1:n1t2)), (/n12/))
        end if
      end do
    end do

    deallocate(tmpCol)
    deallocate(tmpMat)
    deallocate(phases)


  end subroutine unfoldFromCSR_Pauli


  !> Creates a CSR matrix by cloning the sparsity structure of another CSR matrix (real version).
  !>
  !> Note: The elements of the original matrix are not copied.
  subroutine cloneSparsityMap_real_real(csrOut, csrIn)

    !> New CSR matrix on exit.
    type(r_CSR), intent(inout) :: csrOut

    !> CSR matrix with the sparsity map to be cloned.
    type(r_CSR), intent(in) :: csrIn

    call create(csrOut, csrIn%nRow, csrIn%nCol, csrIn%nnz)
    csrOut%colind = csrIn%colind
    csrOut%rowpnt = csrIn%rowpnt
    csrOut%sorted = csrIn%sorted

  end subroutine cloneSparsityMap_real_real



  !> Creates a CSR matrix by cloning the sparsity structure of another CSR matrix (complex version).
  !>
  !> Note: The elements of the original matrix are not copied.
  subroutine cloneSparsityMap_cplx_cplx(csrOut, csrIn)

    !> New CSR matrix on exit.
    type(z_CSR), intent(inout) :: csrOut

    !> CSR matrix with the sparsity map to be cloned.
    type(z_CSR), intent(in) :: csrIn

    call create(csrOut, csrIn%nRow, csrIn%nCol, csrIn%nnz)
    csrOut%colind = csrIn%colind
    csrOut%rowpnt = csrIn%rowpnt
    csrOut%sorted = csrIn%sorted

  end subroutine cloneSparsityMap_cplx_cplx



  !> Creates a CSR matrix by cloning the sparsity structure of an other CSR matrix (complex
  !> version).
  !>
  !> Note: The elements of the original matrix are not copied.
  subroutine cloneSparsityMap_real_cplx(csrOut, csrIn)

    !> New CSR matrix on exit.
    type(r_CSR), intent(inout) :: csrOut

    !> CSR matrix with the sparsity map to be cloned.
    type(z_CSR), intent(in) :: csrIn

    call create(csrOut, csrIn%nRow, csrIn%nCol, csrIn%nnz)
    csrOut%colind = csrIn%colind
    csrOut%rowpnt = csrIn%rowpnt
    csrOut%sorted = csrIn%sorted

  end subroutine cloneSparsityMap_real_cplx



  !> Creates a CSR matrix by cloning the sparsity structure of another CSR matrix (complex version).
  !>
  !> Note: The elements of the original matrix are not actually copied (only the sparsity pattern).
  subroutine cloneSparsityMap_cplx_real(csrOut, csrIn)

    !> New CSR matrix on exit.
    type(z_CSR), intent(inout) :: csrOut

    !> CSR matrix with the sparsity map to be cloned.
    type(r_CSR), intent(in) :: csrIn

    call create(csrOut, csrIn%nRow, csrIn%nCol, csrIn%nnz)
    csrOut%colind = csrIn%colind
    csrOut%rowpnt = csrIn%rowpnt
    csrOut%sorted = csrIn%sorted

  end subroutine cloneSparsityMap_cplx_real



  !> Creates an empty CSR matrix containing zero elements (real version).
  subroutine createEmptyCSR_real(csr)

    !> Empty CSR matrix on exit.
    type(r_CSR), intent(out) :: csr

    call create(csr, 0, 0, 0)
    csr%nnz = 0
    csr%nRow = 0
    csr%nCol = 0
    csr%sorted = .false.

  end subroutine createEmptyCSR_real



  !> Creates an empty CSR matrix containing zero elements (complex version).
  subroutine createEmptyCSR_cplx(csr)

    !> Empty CSR matrix on exit.
    type(z_CSR), intent(out) :: csr

    call create(csr, 0, 0, 0)
    csr%nnz = 0
    csr%nRow = 0
    csr%nCol = 0
    csr%sorted = .false.

  end subroutine createEmptyCSR_cplx

  !> Destroy a real CSR matrix
  subroutine rdestroy_CSR(mat)

    !> matrix
    type(r_CSR) :: mat

    if (allocated(mat%nzval)) then
      call destroy(mat)
    end if

  end subroutine rdestroy_CSR

  !> Destroy a real DNS matrix
  subroutine rdestroy_DNS(mat)

    !> matrix
    type(r_DNS) :: mat

    if (allocated(mat%val)) then
      call destroy(mat)
    end if

  end subroutine rdestroy_DNS

  !> Destroy a complex CSR matrix
  subroutine zdestroy_CSR(mat)

    !> matrix
    type(z_CSR) :: mat

    if (allocated(mat%nzval)) then
      call destroy(mat)
    end if

  end subroutine zdestroy_CSR

  !> Destroy a complex DNS matrix
  subroutine zdestroy_DNS(mat)

    !> matrix
    type(z_DNS) :: mat

    if (allocated(mat%val)) then
      call destroy(mat)
    end if

  end subroutine zdestroy_DNS


end module dftbp_transport_matconv
