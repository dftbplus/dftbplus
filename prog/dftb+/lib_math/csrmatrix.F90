!--------------------------------------------------------------------------------------------------!
!  DFTB+: general package for performing fast atomistic simulations                                !
!  Copyright (C) 2017  DFTB+ developers group                                                      !
!                                                                                                  !
!  See the LICENSE file for terms of usage and distribution.                                       !
!--------------------------------------------------------------------------------------------------!

#:include 'common.fypp'

!!* Contains functions for transforming between the compressed sparse format
!!* and the internal sparse matrix format.
module csrmatrix
  use assert
  use accuracy
  implicit none
  private

  public :: r_CSR, foldToCSR, unfoldFromCSR


  type r_CSR
    integer :: nnz
    integer :: nrow
    integer :: ncol
    real(dp), allocatable :: nzval(:)
    integer, allocatable :: colind(:)
    integer, allocatable :: rowpnt(:)
  end type r_CSR


  interface foldToCSR
    module procedure foldToCSR_real
  end interface


  interface unfoldFromCSR
    module procedure unfoldFromCSR_real
  end interface


contains


  !!* Folds the internal sparse formatted matrix and converts it to
  !!* the compressed sparse row (csr) format (real version).
  !!* @param sparse The sparse matrix to convert
  !!* @param iAtomStart Starting positions of the atoms in the square matrix
  !!* @param iPair Starting position of atom-neighbor interaction in the sparse
  !!*   matrix.
  !!* @param iNeighbor Index of neighbors
  !!* @param nNeighbor Number of neighbors
  !!* @param img2CentCell Image of the atoms in the central cell.
  !!* @param mAngAtom Angular momentum of each atom
  !!* @param mmAng Maximal angular momentum in the system.
  !!* @param csr Resulting CSR matrix
  !!* @note In the resulting CSR format both triangles of the matrix are
  !!*   set.
  subroutine foldToCSR_real(sparse, iAtomStart, iPair, &
      &iNeighbor, nNeighbor, img2CentCell, mAngAtom, mmAng, &
      &csr)
    real(dp), intent(in) :: sparse(:)
    integer, intent(in) :: iAtomStart(:)
    integer, intent(in) :: iPair(0:,:)
    integer, intent(in) :: iNeighbor(0:,:)
    integer, intent(in) :: nNeighbor(:)
    integer, intent(in) :: img2CentCell(:)
    integer, intent(in) :: mAngAtom(:)
    integer, intent(in) :: mmAng
    type(r_CSR), intent(out) :: csr

    integer :: nOrb(size(mAngAtom))
    integer, allocatable :: nColAtom(:), nCols(:), iNonZero(:)
    real(dp), allocatable :: tmpCol(:,:)
    logical, allocatable :: zero(:)
    integer :: nAtom, nNonZero
    integer :: iOrb1, iOrb2, nOrb1, nOrb2, iAt1, iAt2, iAt2f, iNeigh
    integer :: iRow, iCol, ind
    integer :: ii, jj


    nOrb(:) = (mAngAtom(:)+1)**2
    nAtom = size(mAngAtom)

    !! Count nr. of nonzero columns in the square (folded) form for each atom
    allocate(nColAtom(nAtom))
    allocate(zero(nAtom))
    nColAtom(:) = 0
    do iAt1 = 1, nAtom
      zero(:) = .true.
      nOrb1 = nOrb(iAt1)
      do iNeigh = 0, nNeighbor(iAt1)
        iAt2f = img2CentCell(iNeighbor(iNeigh,iAt1))
        if (zero(iAt2f)) then
          nColAtom(iAt1) = nColAtom(iAt1) + nOrb(iAt2f)
          zero(iAt2f) = .false.
          if (iAt1 /= iAt2f) then
            nColAtom(iAt2f) = nColAtom(iAt2f) + nOrb1
          end if
        end if
      end do
    end do

    csr%nRow = iAtomStart(nAtom) + nOrb(nAtom) - 1
    csr%nCol = csr%nRow
    allocate(csr%rowpnt(csr%nRow + 1))

    !! Calculate CSR row pointers
    csr%rowpnt(1) = 1
    ind = 2
    do iAt1 = 1, nAtom
      ii = nColAtom(iAt1)
      do iOrb1 = 1, nOrb(iAt1)
        csr%rowpnt(ind) = csr%rowpnt(ind-1) + ii
        ind = ind + 1
      end do
    end do

    csr%nnz = csr%rowpnt(csr%nRow + 1) - 1
    allocate(csr%nzval(csr%nnz))
    allocate(csr%colind(csr%nnz))

    !! Initialize auxiliary arrays

    ! Nr. of CSR columns already filled
    allocate(nCols(csr%nRow))
    nCols(:) = 0
    ! One block column of the mtx
    allocate(tmpCol(csr%nRow, (mmAng+1)**2))
    ! Index of the nonzero blocks
    allocate(iNonZero(nAtom))

    !! Loop over all atoms (over all block columns in the rectangular picture)
    lpAt1: do iAt1 = 1, nAtom
      iCol = iAtomStart(iAt1)
      nOrb1 = nOrb(iAt1)
      tmpCol(:,:) = 0.0_dp
      nNonZero = 0                         ! Nr. of the nonzero blocks
      zero(:) = .true.

      !! Fold back contributions from the periodic images
      do iNeigh = 0, nNeighbor(iAt1)
        ind = iPair(iNeigh,iAt1) + 1       ! pos. in the sparse mtx
        iAt2 = iNeighbor(iNeigh,iAt1)
        iAt2f = img2CentCell(iAt2)
        nOrb2 = nOrb(iAt2f)
        iRow = iAtomStart(iAt2f)
        tmpCol(iRow:iRow+nOrb2-1, 1:nOrb1) = &
            &tmpCol(iRow:iRow+nOrb2-1, 1:nOrb1) &
            &+ reshape(sparse(ind:ind+nOrb2*nOrb1-1), (/ nOrb2, nOrb1 /))
        if (zero(iAt2f)) then
          zero(iAt2f) = .false.
          nNonZero = nNonZero + 1
          iNonzero(nNonzero) = iAt2f
        end if
      end do

      !! Convert the nonzero blocks of the folded column into CSR format
      lpNonZero: do ii = 1, nNonZero
        iAt2f = iNonZero(ii)
        nOrb2 = nOrb(iAt2f)
        iRow = iAtomStart(iAt2f)

        !! Symmetrize onsite blocks (just in case)
        if (iAt2f == iAt1) then
          do iOrb2 = 1, nOrb2 - 1
            tmpCol(iRow+iOrb2-1,iOrb2+1:nOrb2) = &
                &tmpCol(iRow+iOrb2:iRow+nOrb2-1,iOrb2)
          end do
        end if

        !! Put the rows of current the block into the lower triangle of the CSR
        !! matrix
        do iOrb2 = 1, nOrb2
          jj = iRow + iOrb2 - 1
          ind = csr%rowpnt(jj) + nCols(jj)
          csr%nzval(ind:ind+nOrb1-1) = tmpCol(iRow+iOrb2-1,1:nOrb1)
          nCols(jj) = nCols(jj) + nOrb1
          do iOrb1 = 0, nOrb1 - 1
            csr%colind(ind+iOrb1) = iCol + iOrb1
          end do
        end do

        !! Put the columns of the current block as rows into upper triangle
        !! of the CSR matrix
        if (iAt2f /= iAt1) then
          do iOrb1 = 1, nOrb1
            jj = iCol+iOrb1-1
            ind = csr%rowpnt(jj) + nCols(jj)
            csr%nzval(ind:ind+nOrb2-1) = tmpCol(iRow:iRow+nOrb2-1,iOrb1)
            nCols(jj) = nCols(jj) + nOrb2
            do iOrb2 = 0, nOrb2 - 1
              csr%colind(ind+iOrb2) = iRow + iOrb2
            end do
          end do
        end if
      end do lpNonZero
    end do lpAt1

  end subroutine foldToCSR_real



  !!* Unfolds a matrix from the CSR form into the internal sparse
  !!* representation (real version).
  !!* @param csr CSR matrix
  !!* @param iAtomStart Starting positions of the atoms in the square matrix
  !!* @param iPair Starting position of atom-neighbor interaction in the sparse
  !!*   matrix.
  !!* @param iNeighbor Index of neighbors
  !!* @param nNeighbor Number of neighbors
  !!* @param img2CentCell Image of the atoms in the central cell.
  !!* @param mAngAtom Angular momentum of each atom
  !!* @param mmAng Maximal angular momentum in the system.
  !!* @param sparse The sparse matrix to convert
  !!* @note The CSR matrix must be symmetric.
  subroutine unfoldFromCSR_real(csr, iAtomStart, iPair, iNeighbor, nNeighbor, &
      &img2CentCell, mAngAtom, mmAng, sparse)
    type(r_CSR), intent(in) :: csr
    integer, intent(in) :: iAtomStart(:)
    integer, intent(in) :: iPair(0:,:)
    integer, intent(in) :: iNeighbor(0:,:)
    integer, intent(in) :: nNeighbor(:)
    integer, intent(in) :: img2CentCell(:)
    integer, intent(in) :: mAngAtom(:)
    integer, intent(in) :: mmAng
    real(dp), intent(inout) :: sparse(:)

    integer :: nOrb(size(mAngAtom))
    real(dp), allocatable :: tmpCol(:,:)
    integer :: nAtom
    integer :: iOrb1, nOrb1, nOrb2, iAt1, iAt2, iAt2f, iNeigh
    integer :: iRow, iCol, ind
    integer :: ii

    nOrb(:) = (mAngAtom(:)+1)**2
    nAtom = size(mAngAtom)

    @:ASSERT(csr%nRow == iAtomStart(nAtom) + nOrb(nAtom) - 1)
    allocate(tmpCol(csr%nRow, (mmAng+1)**2))

    do iAt1 = 1, nAtom
      !! Put the rows belonging to a certain atom into the appropriate column
      !! of the block column. (Matrix must be symmetric!)
      tmpCol(:,:) = 0.0_dp
      nOrb1 = nOrb(iAt1)
      iCol = iAtomStart(iAt1)
      do iOrb1 = 1, nOrb(iAt1)
        ii = iCol + iOrb1 - 1
        do ind = csr%rowpnt(ii), csr%rowpnt(ii+1) - 1
          tmpCol(csr%colind(ind), iOrb1) = csr%nzval(ind)
        end do
      end do

      !! Unfold the block column into the internal sparse format
      do iNeigh = 0, nNeighbor(iAt1)
        ind = iPair(iNeigh,iAt1) + 1
        iAt2 = iNeighbor(iNeigh, iAt1)
        iAt2f = img2CentCell(iAt2)
        nOrb2 = nOrb(iAt2f)
        iRow = iAtomStart(iAt2f)
        sparse(ind:ind+nOrb2*nOrb1-1) = &
            &sparse(ind:ind+nOrb2*nOrb1-1) &
            &+ reshape(tmpCol(iRow:iRow+nOrb2-1, 1:nOrb1), (/nOrb2*nOrb1/))
      end do
    end do

  end subroutine unfoldFromCSR_real



end module csrmatrix
