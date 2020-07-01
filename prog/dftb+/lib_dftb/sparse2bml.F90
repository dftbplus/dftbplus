!--------------------------------------------------------------------------------------------------!
!  DFTB+: general package for performing fast atomistic simulations                                !
!  Copyright (C) 2018  DFTB+ developers group                                                      !
!                                                                                                  !
!  See the LICENSE file for terms of usage and distribution.                                       !
!--------------------------------------------------------------------------------------------------!

#:include 'common.fypp'

!> Contains routines to convert between DFTB+ sparse matrices and BML matrices.
module sparse2bml
  use, intrinsic :: iso_c_binding, only : C_INT, C_DOUBLE_COMPLEX
  use accuracy
  use commontypes
  use constants
  use dftbp_bml
  implicit none
  private

  public :: foldToRealBml, foldToComplexBml
  public :: unfoldFromRealBml, unfoldFromComplexBml

contains

  !> Folds a sparse matrix to ELL format as used in the BML library (real).
  !>
  subroutine foldToRealBml(sparse, iNeighbor, nNeighbor, orb, iSquare, iPair, img2CentCell,&
      & bmlMatrix)

    !> Sparse matrix
    real(dp), intent(in) :: sparse(:)

    !> Neighbor list for each atom.
    integer, intent(in) :: iNeighbor(0:,:)

    !> Nr. of neighbors for each atom (incl. itself)
    integer, intent(in) :: nNeighbor(:)

    !> Orbital (basis) information.
    type(TOrbitals), intent(in) :: orb

    !> Atom offset for the squared Hamiltonian.
    integer, intent(in) :: iSquare(:)

    !> Indexing array for the sparse Hamiltonian.
    integer, intent(in) :: iPair(0:,:)

    !> Atomic mapping indexes
    integer, intent(in) :: img2CentCell(:)

    !> Sparse matrix in BML format
    type(bml_matrix_t), intent(inout) :: bmlMatrix

    real(dp), allocatable :: buffer(:,:)
    integer :: iAt1, iAt2f, iNeigh, iOrb
    integer :: nAtom, nOrb1, nOrb2, globCol, globRow, indSparse

    nAtom = size(iNeighbor, dim=2)
    allocate(buffer(orb%nOrb, orb%mOrb))
    do iAt1 = 1, size(iNeighbor, dim=2)
      globCol = iSquare(iAt1)
      buffer(:,:) = 0.0_dp
      nOrb1 = orb%nOrbAtom(iAt1)
      do iNeigh = 0, nNeighbor(iAt1)
        indSparse = iPair(iNeigh, iAt1) + 1
        iAt2f = img2CentCell(iNeighbor(iNeigh, iAt1))
        globRow = iSquare(iAt2f)
        nOrb2 = orb%nOrbAtom(iAt2f)
        buffer(globRow : globRow + nOrb2 - 1, 1:nOrb1) =&
            & buffer(globRow : globRow + nOrb2 - 1, 1:nOrb1) &
            & + reshape(sparse(indSparse : indSparse + nOrb2 * nOrb1 - 1), [nOrb2, nOrb1])
      end do
      do iOrb = 1, nOrb1
        call bml_set_row(bmlMatrix, globCol + iOrb - 1, buffer(:, iOrb))
      end do
    end do

    call bml_adjungate_triangle(bmlMatrix, "u")

  end subroutine foldToRealBml


  !> Folds sparse matrix to ELL format as used in the BML library (complex).
  !>
  subroutine foldToComplexBml(sparse, kPoint, iNeighbor, nNeighbor, iCellVec, cellVec, orb,&
      & iSquare, iPair, img2CentCell, bmlMatrix)

    !> Sparse matrix
    real(dp), intent(in) :: sparse(:)

    !> K-point which the matrix belongs to
    real(dp), intent(in) :: kPoint(:)

    !> Neighbor list for each atom.
    integer, intent(in) :: iNeighbor(0:,:)

    !> Nr. of neighbors for each atom (incl. itself).
    integer, intent(in) :: nNeighbor(:)

    !> Index of the cell translation vector for each atom.
    integer, intent(in) :: iCellVec(:)

    !> Relative coordinates of the cell translation vectors.
    real(dp), intent(in) :: cellVec(:,:)

    !> Orbital (basis) information.
    type(TOrbitals), intent(in) :: orb

    !> Atom offset for the squared Hamiltonian.
    integer, intent(in) :: iSquare(:)

    !> Indexing array for the sparse Hamiltonian.
    integer, intent(in) :: iPair(0:,:)

    !> Atomic mapping indexes
    integer, intent(in) :: img2CentCell(:)

    !> Sparse matrix in BML format.
    type(bml_matrix_t), intent(inout) :: bmlMatrix

    complex(dp), allocatable :: buffer(:,:)
    real(dp) :: kPoint2p(3)
    complex(dp) :: phase
    integer :: iAt1, iAt2, iAt2f, iNeigh, iOrb, iVec, iOldVec
    integer :: nAtom, nOrb1, nOrb2, globCol, globRow, indSparse

    nAtom = size(iNeighbor, dim=2)
    allocate(buffer(orb%nOrb, orb%mOrb))
    kPoint2p(:) = 2.0_dp * pi * kPoint
    iOldVec = 0
    phase = (1.0_dp, 0.0_dp)
    do iAt1 = 1, size(iNeighbor, dim=2)
      globCol = iSquare(iAt1)
      buffer(:,:) = (0.0_dp, 0.0_dp)
      nOrb1 = orb%nOrbAtom(iAt1)
      do iNeigh = 0, nNeighbor(iAt1)
        indSparse = iPair(iNeigh,iAt1) + 1
        iAt2 = iNeighbor(iNeigh, iAt1)
        iAt2f = img2CentCell(iAt2)
        globRow = iSquare(iAt2f)
        nOrb2 = orb%nOrbAtom(iAt2f)
        iVec = iCellVec(iAt2)
        if (iVec /= iOldVec) then
          phase = exp((0.0_dp, 1.0_dp) * dot_product(kPoint2p, cellVec(:, iVec)))
          iOldVec = iVec
        end if
        buffer(globRow : globRow + nOrb2 - 1, 1:nOrb1) = &
            & buffer(globRow : globRow + nOrb2 - 1, 1:nOrb1) &
            & + phase * reshape(sparse(indSparse : indSparse + nOrb1 * nOrb2 - 1), [nOrb2, nOrb1])
      end do
      do iOrb = 1, nOrb1
        ! We set a column of sparse as a row of bml -> conjugate
        call bml_set_row(bmlMatrix, iOrb + globCol - 1, conjg(buffer(:, iOrb)))
      end do
    end do
    call bml_adjungate_triangle(bmlMatrix, "u")

  end subroutine foldToComplexBml


  !> Unfolds the contribution of a BML matrix to the sparse form (real).
  !>
  !> The unfolded matrix is *added* to the sparse one.
  !>
  subroutine unfoldFromRealBml(bmlMatrix, iNeighbor, nNeighbor, orb, iSquare, iPair, img2CentCell,&
      & sparse)

    !> BML matrix
    type(bml_matrix_t), intent(in) :: bmlMatrix

    !> Neighbor list for the atoms.
    integer, intent(in) :: iNeighbor(0:,:)

    !> Nr. of neighbors for the atoms.
    integer, intent(in) :: nNeighbor(:)

    !> Orbital information.
    type(TOrbitals), intent(in) :: orb

    !> Atom offset for the squared matrix.
    integer, intent(in) :: iSquare(:)

    !> Indexing array for the sparse Hamiltonian.
    integer, intent(in) :: iPair(0:,:)

    !> Mapping between image atoms and correspondent atom in the central cell.
    integer, intent(in) :: img2CentCell(:)

    !> Updated sparse matrix.
    real(dp), intent(inout) :: sparse(:)

    integer :: iAt1, iAt2f, iNeigh, iOrb, iOrig
    integer :: nAtom, nOrb1, nOrb2, globCol, globRow
    real(dp), allocatable :: buffer(:,:), rowBuffer(:)

    nAtom = size(iNeighbor, dim=2)
    allocate(buffer(orb%nOrb, orb%mOrb))
    do iAt1 = 1, nAtom
      nOrb1 = orb%nOrbAtom(iAt1)
      globCol = iSquare(iAt1)
      do iOrb = 1, nOrb1
        call bml_get_row(bmlMatrix, globCol + iOrb - 1, rowBuffer)
        buffer(:,iOrb) = rowBuffer
      end do
      do iNeigh = 0, nNeighbor(iAt1)
        iOrig = iPair(iNeigh, iAt1) + 1
        iAt2f = img2CentCell(iNeighbor(iNeigh, iAt1))
        globRow = iSquare(iAt2f)
        nOrb2 = orb%nOrbAtom(iAt2f)
        sparse(iOrig : iOrig + nOrb1 * nOrb2 - 1) = &
            & sparse(iOrig : iOrig + nOrb1 * nOrb2 - 1) &
            & + reshape(buffer(globRow : globRow + nOrb2 - 1, 1:nOrb1), [nOrb1 * nOrb2])
      end do
    end do

  end subroutine unfoldFromRealBml


  !> Unfolds the contribution of a BML matrix to the sparse form (complex).
  !>
  !> The unfolded matrix is *added* to the sparse one weighted by the passed k-point weight.
  !>
  subroutine unfoldFromComplexBml(bmlMatrix, kPoint, kWeight, iNeighbor, nNeighbor, orb, iCellVec,&
      & cellVec, iSquare, iPair, img2CentCell, sparse)

    !> BML matrix.
    type(bml_matrix_t), intent(in) :: bmlMatrix

    !> K-point to which the BML matrix belongs to.
    real(dp), intent(in) :: kPoint(:)

    !> Weight of the k-point.
    real(dp), intent(in) :: kWeight

    !> Neighbor list for the atoms.
    integer, intent(in) :: iNeighbor(0:,:)

    !> Nr. of neighbors for the atoms.
    integer, intent(in) :: nNeighbor(:)

    !> Orbital information.
    type(TOrbitals), intent(in) :: orb

    !> Index of the cell translation vector for each atom.
    integer, intent(in) :: iCellVec(:)

    !> Relative coordinates of the cell translation vectors.
    real(dp), intent(in) :: cellVec(:,:)

    !> Atom offset for the squared matrix
    integer, intent(in) :: iSquare(:)

    !> Indexing array for the sparse Hamiltonian.
    integer, intent(in) :: iPair(0:,:)

    !> Mapping between image atoms and correspondent atom in the central cell.
    integer, intent(in) :: img2CentCell(:)

    !> Updated sparse matrix.
    real(dp), intent(inout) :: sparse(:)

    real(dp) :: kPoint2p(3)
    complex(dp) :: phase
    complex(C_DOUBLE_COMPLEX), allocatable :: buffer(:,:), rowBuffer(:)
    integer :: nAtom
    integer :: iAt1, iAt2, iAt2f, iNeigh, iOrb, iVec, iOldVec, iOrig
    integer :: nOrb1, nOrb2, globCol, globRow

    nAtom = size(iNeighbor, dim=2)

    allocate(buffer(orb%nOrb, orb%mOrb))
    iOldVec = 0
    kPoint2p = 2.0_dp * pi * kPoint
    phase = 1.0_dp
    do iAt1 = 1, nAtom
      nOrb1 = orb%nOrbAtom(iAt1)
      globCol = iSquare(iAt1)
      do iOrb = 1, nOrb1
        call bml_get_row(bmlMatrix, int(globCol + iOrb - 1, kind=C_INT), rowBuffer)
        buffer(:,iOrb) = rowBuffer
      end do
      ! Taking a row of BML as column in sparse -> conjugate
      buffer(:,:) = conjg(buffer)
      do iNeigh = 0, nNeighbor(iAt1)
        iOrig = iPair(iNeigh, iAt1) + 1
        iAt2 = iNeighbor(iNeigh, iAt1)
        iAt2f = img2CentCell(iAt2)
        globRow = iSquare(iAt2f)
        nOrb2 = orb%nOrbAtom(iAt2f)
        iVec = iCellVec(iAt2)
        if (iVec /= iOldVec) then
          phase = exp(cmplx(0, -1, dp) &
              & * dot_product(kPoint2p, cellVec(:, iVec)))
          iOldVec = iVec
        end if

        sparse(iOrig : iOrig + nOrb1 * nOrb2 - 1) = sparse(iOrig : iOrig + nOrb1 * nOrb2 - 1) &
            & + kWeight * real(phase * reshape(buffer(globRow : globRow + nOrb2 - 1, 1:nOrb1), &
            & [nOrb1 * nOrb2]), dp)
      end do
    end do

  end subroutine unfoldFromComplexBml


end module sparse2bml
