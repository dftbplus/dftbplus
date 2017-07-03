!--------------------------------------------------------------------------------------------------!
!  DFTB+: general package for performing fast atomistic simulations                                !
!  Copyright (C) 2017  DFTB+ developers group                                                      !
!                                                                                                  !
!  See the LICENSE file for terms of usage and distribution.                                       !
!--------------------------------------------------------------------------------------------------!

#:include 'common.fypp'

!!* Contains subroutines for packing/unpacking Hamiltonian-like matrices
!!* between the square and 1-dimensional representations
module sparse2dense
  use assert
  use accuracy
  use constants, only : pi
  use commontypes
  use memman
  use periodic, only : TNeighborList
  implicit none

  private

  public :: unpackHS, packHS, iPackHS, packErho
  public :: blockSymmetrizeHS,  blockHermitianHS, blockAntiSymmetrizeHS

  !!* Unpack sparse matrix (Hamiltonian, overlap, etc.) to square form
  interface unpackHS
    module procedure unpackHS_real
    module procedure unpackHS_real_kpts
    module procedure unpackHS_cmplx
  end interface

  !!* Pack square matrix to sparse form.
  interface packHS
    module procedure packHS_real
    module procedure packHS_cmplx
    module procedure packHSPauli
    module procedure packHSPauli_kpts
  end interface

  !!* Pack square matrix to sparse form.
  interface iPackHS
    module procedure packHSPauliImag
    module procedure packHSPauliImag_kpts
  end interface


  !!* Pack energy weighted Pauli idenity square matrix to sparse form.
  interface packErho
    module procedure packHSPauliERho
    module procedure packHSPauliERho_kpts
  end interface

  !!* Symmetrize the square matrix except the on-site blocks
  interface blockSymmetrizeHS
    module procedure blockSymmetrizeHS_real
    module procedure blockSymmetrizeHS_cmplx
  end interface

  !!* Hermitian the square matrix except the on-site blocks
  interface blockHermitianHS
    module procedure blockSymmetrizeHS_real
    module procedure blockHermitianHS_cmplx
  end interface

   !!* Symmetrize the square matrix except the on-site blocks
  interface blockAntiSymmetrizeHS
    module procedure blockAntiSymmetrizeHS_real
  end interface


contains

  !!* Unpacks sparse matrix to square form (complex version)
  !!* @param square       Square form matrix on exit.
  !!* @param orig         Sparse matrix
  !!* @param kPoint       <b>Relative</b> coordinates of the K-point where the
  !!*   sparse  matrix    should be unfolded.
  !!* @param iNeighbor    Neighbor list for each atom (First index from 0!)
  !!* @param nNeighbor    Nr. of neighbors for each atom (incl. itself).
  !!* @param iCellVec     Index of the cell translation vector for each atom.
  !!* @param cellVec      Relative coordinates of the cell translation vectors.
  !!* @param iAtomStart   Atom offset for the squared Hamiltonian
  !!* @param iPair indexing array for the sparse Hamiltonian
  !!* @param img2CentCell Atomic mapping indexes.
  !!* @note The non on-site blocks are only filled in the lower triangle part
  !!*   of the matrix. To fill the matrix completely, apply the
  !!*   blockSymmetrizeHS subroutine.
  subroutine unpackHS_cmplx(square, orig, kPoint, iNeighbor, nNeighbor, &
      &iCellVec, cellVec, iAtomStart, iPair, img2CentCell)
    complex(dp), intent(out) :: square(:,:)
    real(dp),    intent(in)  :: orig(:)
    real(dp),    intent(in)  :: kPoint(:)
    integer,     intent(in)  :: iNeighbor(0:,:)
    integer,     intent(in)  :: nNeighbor(:)
    integer,     intent(in)  :: iCellVec(:)
    real(dp),    intent(in)  :: cellVec(:,:)
    integer,     intent(in)  :: iAtomStart(:)
    integer,     intent(in)  :: iPair(0:,:)
    integer,     intent(in)  :: img2CentCell(:)

    complex(dp) :: phase
    integer     :: nAtom
    integer     :: iOrig, ii, jj
    integer     :: iNeigh
    integer     :: iOldVec, iVec
    integer     :: iAtom1, iAtom2, iAtom2f
    integer     :: nOrb1, nOrb2
    real(dp)    :: kPoint2p(3)

    nAtom = size(iNeighbor, dim=2)

    @:ASSERT(nAtom > 0)
    @:ASSERT(size(square, dim=1) == size(square, dim=2))
    @:ASSERT(size(square, dim=1) == iAtomStart(nAtom+1) - 1)
    @:ASSERT(all(shape(kPoint) == (/ 3 /)))
    @:ASSERT(all(shape(nNeighbor) == (/ nAtom /)))
    @:ASSERT(size(iAtomStart) == nAtom + 1)

    square(:,:) = cmplx(0, 0, dp)
    kPoint2p(:) = 2.0_dp * pi * kPoint(:)
    iOldVec = 0
    phase = 1.0_dp
    do iAtom1 = 1, nAtom
      ii = iAtomStart(iAtom1)
      nOrb1 = iAtomStart(iAtom1+1) - ii
      do iNeigh = 0, nNeighbor(iAtom1)
        iOrig = iPair(iNeigh,iAtom1) + 1
        iAtom2 = iNeighbor(iNeigh, iAtom1)
        iAtom2f = img2CentCell(iAtom2)
        jj = iAtomStart(iAtom2f)
        @:ASSERT(jj >= ii)
        nOrb2 = iAtomStart(iAtom2f+1) - jj
        iVec = iCellVec(iAtom2)
        if (iVec /= iOldVec) then
          phase = exp((0.0_dp, 1.0_dp) &
              &* dot_product(kPoint2p, cellVec(:, iVec)))
          iOldVec = iVec
        end if
        square(jj:jj+nOrb2-1, ii:ii+nOrb1-1) = &
            & square(jj:jj+nOrb2-1, ii:ii+nOrb1-1) + phase &
            & * reshape(orig(iOrig:iOrig+nOrb1*nOrb2-1), (/nOrb2, nOrb1/))
      end do
    end do

  end subroutine unpackHS_cmplx

  !!* Unpacks sparse matrix to square form (only real part of the phase factor
  !!* is considered)
  !!* @param square       Square form matrix on exit.
  !!* @param orig         Sparse matrix
  !!* @param kPoint       <b>Relative</b> coordinates of the K-point where the
  !!*   sparse  matrix    should be unfolded.
  !!* @param iNeighbor    Neighbor list for each atom (First index from 0!)
  !!* @param nNeighbor    Nr. of neighbors for each atom (incl. itself).
  !!* @param iCellVec     Index of the cell translation vector for each atom.
  !!* @param cellVec      Relative coordinates of the cell translation vectors.
  !!* @param iAtomStart   Atom offset for the squared Hamiltonian
  !!* @param iPair indexing array for the sparse Hamiltonian
  !!* @param img2CentCell Atomic mapping indexes.
  !!* @note The non on-site blocks are only filled in the lower triangle part
  !!*   of the matrix. To fill the matrix completely, apply the
  !!*   blockSymmetrizeHS subroutine.
  subroutine unpackHS_real_kpts(square, orig, kPoint, iNeighbor, nNeighbor, &
      &iCellVec, cellVec, iAtomStart, iPair, img2CentCell)
    real(dp), intent(out) :: square(:,:)
    real(dp), intent(in)  :: orig(:)
    real(dp), intent(in)  :: kPoint(:)
    integer,  intent(in)  :: iNeighbor(0:,:)
    integer,  intent(in)  :: nNeighbor(:)
    integer,  intent(in)  :: iCellVec(:)
    real(dp), intent(in)  :: cellVec(:,:)
    integer,  intent(in)  :: iAtomStart(:)
    integer,     intent(in)  :: iPair(0:,:)
    integer,  intent(in)  :: img2CentCell(:)

    real(dp) :: phase
    integer  :: nAtom
    integer  :: iOrig, ii, jj
    integer  :: iNeigh
    integer  :: iOldVec, iVec
    integer  :: iAtom1, iAtom2, iAtom2f
    integer  :: nOrb1, nOrb2
    real(dp) :: kPoint2p(3)


    nAtom = size(iNeighbor, dim=2)

    @:ASSERT(nAtom > 0)
    @:ASSERT(size(square, dim=1) == size(square, dim=2))
    @:ASSERT(size(square, dim=1) == iAtomStart(nAtom+1) -1)
    @:ASSERT(all(shape(kPoint) == (/ 3 /)))
    @:ASSERT(all(shape(nNeighbor) == (/ nAtom /)))
    @:ASSERT(size(iAtomStart) == nAtom + 1)

    square(:,:) = 0.0_dp
    kPoint2p(:) = 2.0_dp * pi * kPoint(:)
    iOldVec = 0
    phase = 1.0_dp
    do iAtom1 = 1, nAtom
      ii = iAtomStart(iAtom1)
      nOrb1 = iAtomStart(iAtom1+1) - ii
      do iNeigh = 0, nNeighbor(iAtom1)
        iOrig = iPair(iNeigh,iAtom1) + 1
        iAtom2 = iNeighbor(iNeigh, iAtom1)
        iAtom2f = img2CentCell(iAtom2)
        jj = iAtomStart(iAtom2f)
        @:ASSERT(jj >= ii)
        nOrb2 = iAtomStart(iAtom2f+1) - jj
        iVec = iCellVec(iAtom2)
        if (iVec /= iOldVec) then
          phase = cos(dot_product(kPoint2p, cellVec(:, iVec)))
          iOldVec = iVec
        end if
        square(jj:jj+nOrb2-1,ii:ii+nOrb1-1) = &
            & square(jj:jj+nOrb2-1,ii:ii+nOrb1-1) + phase &
            &* reshape(orig(iOrig:iOrig+nOrb1*nOrb2-1),(/nOrb2,nOrb1/))
      end do
    end do
  end subroutine unpackHS_real_kpts



  !!* Unpacks sparse matrix to square form (real version for Gamma point)
  !!* @param square       Square form matrix on exit.
  !!* @param orig         Sparse matrix
  !!* @param iNeighbor    Neighbor list for each atom (First index from 0!)
  !!* @param nNeighbor    Nr. of neighbors for each atom (incl. itself).
  !!* @param iAtomStart   Atom offset for the squared Hamiltonian
  !!* @param iPair indexing array for the sparse Hamiltonian
  !!* @param img2CentCell Atomic mapping indexes.
  !!* @note The non on-site blocks are only filled in the lower triangle part
  !!*   of the matrix. To fill the matrix completely, apply the
  !!*   blockSymmetrizeHS subroutine.
  subroutine unpackHS_real(square, orig, iNeighbor, nNeighbor, iAtomStart, &
      &iPair, img2CentCell)
    real(dp), intent(out) :: square(:,:)
    real(dp), intent(in)  :: orig(:)
    integer,  intent(in)  :: iNeighbor(0:,:)
    integer,  intent(in)  :: nNeighbor(:)
    integer,  intent(in)  :: iAtomStart(:)
    integer,     intent(in)  :: iPair(0:,:)
    integer,  intent(in)  :: img2CentCell(:)

    integer     :: nAtom
    integer     :: iOrig, ii, jj
    integer     :: iNeigh
    integer     :: iAtom1, iAtom2, iAtom2f
    integer     :: nOrb1, nOrb2

    nAtom = size(iNeighbor, dim=2)

    @:ASSERT(nAtom > 0)
    @:ASSERT(size(square, dim=1) == size(square, dim=2))
    @:ASSERT(size(square, dim=1) == iAtomStart(nAtom+1) - 1)
    @:ASSERT(all(shape(nNeighbor) == (/ nAtom /)))
    @:ASSERT(size(iAtomStart) == nAtom + 1)

    square(:,:) = 0.0_dp

    do iAtom1 = 1, nAtom
      ii = iAtomStart(iAtom1)
      nOrb1 = iAtomStart(iAtom1+1) - ii
      do iNeigh = 0, nNeighbor(iAtom1)
        iOrig = iPair(iNeigh,iAtom1) + 1
        iAtom2 = iNeighbor(iNeigh, iAtom1)
        iAtom2f = img2CentCell(iAtom2)
        jj = iAtomStart(iAtom2f)
        @:ASSERT(jj >= ii)
        nOrb2 = iAtomStart(iAtom2f+1) - jj
        square(jj:jj+nOrb2-1, ii:ii+nOrb1-1) = &
            & square(jj:jj+nOrb2-1, ii:ii+nOrb1-1) &
            & + reshape(orig(iOrig:iOrig+nOrb1*nOrb2-1), (/nOrb2,nOrb1/))
      end do
    end do

  end subroutine unpackHS_real

  !!* Pack squared matrix in the sparse form (complex version).
  !!* @param primitive    Sparse matrix
  !!* @param square       Squared form matrix
  !!* @param kPoint       Relative coordinates of the K-point
  !!* @param kWeight      Weight of the K-point
  !!* @param iNeighbor    Neighbor list for the atoms (First index from 0!)
  !!* @param nNeighbor    Nr. of neighbors for the atoms.
  !!* @param mOrb         Maximal number of orbitals on an atom.
  !!* @param iCellVec     Index of the cell translation vector for each atom.
  !!* @param cellVec      Relative coordinates of the cell translation vectors.
  !!* @param iAtomStart   Atom offset for the squared matrix
  !!* @param iPair indexing array for the sparse Hamiltonian
  !!* @param img2CentCell Mapping between image atoms and correspondent atom in
  !!*   the central cell.
  subroutine packHS_cmplx(primitive, square, kPoint, kWeight, iNeighbor, &
      & nNeighbor, mOrb, iCellVec, cellVec, iAtomStart, iPair, img2CentCell)
    real(dp),    intent(inout) :: primitive(:)
    complex(dp), intent(in)    :: square(:,:)
    real(dp),    intent(in)    :: kPoint(:)
    real(dp),    intent(in)    :: kweight
    integer,     intent(in)    :: iNeighbor(0:,:)
    integer,     intent(in)    :: nNeighbor(:)
    integer,     intent(in)    :: mOrb
    integer,     intent(in)    :: iCellVec(:)
    real(dp),    intent(in)    :: cellVec(:,:)
    integer,     intent(in)    :: iAtomStart(:)
    integer,     intent(in)    :: iPair(0:,:)
    integer,     intent(in)    :: img2CentCell(:)

    complex(dp) :: phase
    integer     :: nAtom
    integer     :: iOrig, ii, jj, kk
    integer     :: iNeigh
    integer     :: iOldVec, iVec
    integer     :: iAtom1, iAtom2, iAtom2f
    integer     :: nOrb1, nOrb2
    real(dp)    :: kPoint2p(3)
    complex(dp) :: tmpSqr(mOrb,mOrb)
  #:call ASSERT_CODE
    integer :: sizePrim
  #:endcall ASSERT_CODE

    nAtom = size(iNeighbor, dim=2)
  #:call ASSERT_CODE
    sizePrim = size(primitive)
  #:endcall ASSERT_CODE

    @:ASSERT(nAtom > 0)
    @:ASSERT(size(square, dim=1) == size(square, dim=2))
    @:ASSERT(size(square, dim=1) == iAtomStart(nAtom+1) - 1)
    @:ASSERT(all(shape(kPoint) == (/ 3 /)))
    @:ASSERT(all(shape(nNeighbor) == (/ nAtom /)))
    @:ASSERT(kWeight > 0.0_dp)
    @:ASSERT(size(iAtomStart) == nAtom + 1)

    kPoint2p(:) = 2.0_dp * pi * kPoint(:)
    iOldVec = 0
    phase = 1.0_dp
    do iAtom1 = 1, nAtom
      ii = iAtomStart(iAtom1)
      nOrb1 = iAtomStart(iAtom1+1) - ii
      do iNeigh = 0, nNeighbor(iAtom1)
        iOrig = iPair(iNeigh,iAtom1) + 1
        iAtom2 = iNeighbor(iNeigh, iAtom1)
        iAtom2f = img2CentCell(iAtom2)
        jj = iAtomStart(iAtom2f)
        @:ASSERT(jj >= ii)
        nOrb2 = iAtomStart(iAtom2f+1) - jj
        iVec = iCellVec(iAtom2)
        if (iVec /= iOldVec) then
          phase = exp(cmplx(0,-1,dp) &
              &* dot_product(kPoint2p(:), cellVec(:, iVec)))
          iOldVec = iVec
        end if
        tmpSqr(1:nOrb2, 1:nOrb1) = square(jj:jj+nOrb2-1, ii:ii+nOrb1-1)

        !! Hermitian the on-site block before packing, just in case
        if (iAtom1 == iAtom2f) then
          do kk = 1, nOrb2
            tmpSqr(kk, kk+1:nOrb1) = conjg(tmpSqr(kk+1:nOrb1, kk))
          end do
        end if

        @:ASSERT(sizePrim >= iOrig + nOrb1*nOrb2 - 1)
        primitive(iOrig : iOrig + nOrb1*nOrb2 - 1) = &
            & primitive(iOrig : iOrig + nOrb1*nOrb2 - 1) &
            &+ kWeight &
            &* real(phase &
            &* reshape(tmpSqr(1:nOrb2, 1:nOrb1), (/nOrb1*nOrb2/)), dp)
      end do
    end do
  end subroutine packHS_cmplx



  !!* Pack squared matrix in the sparse form (real version).
  !!* @param primitive    Sparse matrix
  !!* @param square       Squared form matrix
  !!* @param iNeighbor    Neighbor list for the atoms (First index from 0!)
  !!* @param nNeighbor    Nr. of neighbors for the atoms.
  !!* @param mOrb         Maximal number of orbitals on an atom.
  !!* @param iCellVec     Index of the cell translation vector for each atom.
  !!* @param cellVec      Relative coordinates of the cell translation vectors.
  !!* @param iAtomStart   Atom offset for the squared matrix
  !!* @param iPair indexing array for the sparse Hamiltonian
  !!* @param img2CentCell Mapping between image atoms and correspondent atom in
  !!*   the central cell.
  subroutine packHS_real(primitive, square, iNeighbor, nNeighbor, mOrb, &
      &iAtomStart, iPair, img2CentCell)
    real(dp), intent(inout) :: primitive(:)
    real(dp), intent(in)    :: square(:,:)
    integer,  intent(in)    :: iNeighbor(0:,:)
    integer,  intent(in)    :: nNeighbor(:)
    integer,  intent(in)    :: mOrb
    integer,  intent(in)    :: iAtomStart(:)
    integer,  intent(in)    :: iPair(0:,:)
    integer,  intent(in)    :: img2CentCell(:)

    integer     :: nAtom
    integer     :: iOrig, ii, jj, kk
    integer     :: iNeigh
    integer     :: iAtom1, iAtom2, iAtom2f
    integer     :: nOrb1, nOrb2
    real(dp)    :: tmpSqr(mOrb,mOrb)
  #:call ASSERT_CODE
    integer :: sizePrim
  #:endcall ASSERT_CODE

    nAtom = size(iNeighbor, dim=2)
  #:call ASSERT_CODE
    sizePrim = size(primitive)
  #:endcall ASSERT_CODE

    @:ASSERT(nAtom > 0)
    @:ASSERT(size(square, dim=1) == size(square, dim=2))
    @:ASSERT(size(square, dim=1) == iAtomStart(nAtom+1) - 1)
    @:ASSERT(all(shape(nNeighbor) == (/ nAtom /)))

    do iAtom1 = 1, nAtom
      ii = iAtomStart(iAtom1)
      nOrb1 = iAtomStart(iAtom1+1) - ii
      do iNeigh = 0, nNeighbor(iAtom1)
        iOrig = iPair(iNeigh,iAtom1) + 1
        iAtom2 = iNeighbor(iNeigh, iAtom1)
        iAtom2f = img2CentCell(iAtom2)
        jj = iAtomStart(iAtom2f)
        @:ASSERT(jj >= ii)
        nOrb2 = iAtomStart(iAtom2f+1) - jj
        tmpSqr(1:nOrb2, 1:nOrb1) = square(jj:jj+nOrb2-1, ii:ii+nOrb1-1)

        ! Symmetrize the on-site block before packing, just in case
        if (iAtom1 == iAtom2f) then
          do kk = 1, nOrb2
            tmpSqr(kk, kk+1:nOrb1) = tmpSqr(kk+1:nOrb1, kk)
          end do
        end if

        @:ASSERT(sizePrim >= iOrig + nOrb1*nOrb2 - 1)
        primitive(iOrig : iOrig + nOrb1*nOrb2 - 1) = &
            &primitive(iOrig : iOrig + nOrb1*nOrb2 - 1) &
            &+ reshape(tmpSqr(1:nOrb2,1:nOrb1), (/nOrb1*nOrb2/))
      end do
    end do
  end subroutine packHS_real


  !!* Pack squared matrix in the sparse form (real Pauli version).
  !!* @param primitive    Sparse matrix
  !!* @param square       Squared form matrix
  !!* @param iNeighbor    Neighbor list for the atoms (First index from 0!)
  !!* @param nNeighbor    Nr. of neighbors for the atoms.
  !!* @param mOrb         Maximal number of orbitals on an atom.
  !!* @param iCellVec     Index of the cell translation vector for each atom.
  !!* @param cellVec      Relative coordinates of the cell translation vectors.
  !!* @param iAtomStart   Atom offset for the squared matrix
  !!* @param iPair indexing array for the sparse Hamiltonian
  !!* @param img2CentCell Mapping between image atoms and correspondent atom in
  !!*   the central cell.
  subroutine packHSPauli(primitive, square, iNeighbor, nNeighbor, mOrb, &
      &iAtomStart, iPair, img2CentCell)
    real(dp), intent(inout)   :: primitive(:,:)
    complex(dp), intent(in) :: square(:,:)
    integer,  intent(in)    :: iNeighbor(0:,:)
    integer,  intent(in)    :: nNeighbor(:)
    integer,  intent(in)    :: mOrb
    integer,  intent(in)    :: iAtomStart(:)
    integer,  intent(in)    :: iPair(0:,:)
    integer,  intent(in)    :: img2CentCell(:)

    integer     :: nAtom
    integer     :: iOrig, ii, jj, kk
    integer     :: iNeigh, iBlock
    integer     :: iAtom1, iAtom2, iAtom2f
    integer     :: nOrb1, nOrb2, nOrb
    complex(dp) :: tmpSqr(mOrb,mOrb)
  #:call ASSERT_CODE
    integer :: sizePrim
  #:endcall ASSERT_CODE

    nAtom = size(iNeighbor, dim=2)
    nOrb = (iAtomStart(nAtom+1) - 1) ! number of orbitals in a regular
                                     ! spin block

  #:call ASSERT_CODE
    sizePrim = size(primitive,dim=1)
  #:endcall ASSERT_CODE

    @:ASSERT(nAtom > 0)
    @:ASSERT(size(square, dim=1) == size(square, dim=2))
    @:ASSERT(size(square, dim=1) == 2 * nOrb )
    @:ASSERT(all(shape(nNeighbor) == (/ nAtom /)))
    @:ASSERT(size(iAtomStart) == nAtom + 1)
    @:ASSERT(size(primitive,dim=2)==4)


    do iBlock = 0, 1
      do iAtom1 = 1, nAtom
        ii = iAtomStart(iAtom1)
        nOrb1 = iAtomStart(iAtom1+1) - ii
        do iNeigh = 0, nNeighbor(iAtom1)
          iOrig = iPair(iNeigh,iAtom1) + 1
          iAtom2 = iNeighbor(iNeigh, iAtom1)
          iAtom2f = img2CentCell(iAtom2)
          jj = iAtomStart(iAtom2f)
          @:ASSERT(jj >= ii)
          nOrb2 = iAtomStart(iAtom2f+1) - jj
          tmpSqr(1:nOrb2, 1:nOrb1) = &
              & square(jj+iBlock*nOrb:jj+nOrb2-1+iBlock*nOrb, &
              & ii+iBlock*nOrb:ii+nOrb1-1+iBlock*nOrb)
          ! Hermitian the on-site block before packing, as only one triangle
          ! usually supplied
          if (iAtom1 == iAtom2f) then
            do kk = 1, nOrb2
              tmpSqr(kk, kk+1:nOrb1) = conjg(tmpSqr(kk+1:nOrb1, kk))
            end do
          end if
          @:ASSERT(sizePrim >= iOrig + nOrb1*nOrb2 - 1)
          primitive(iOrig : iOrig + nOrb1*nOrb2 - 1,1) = &
              &primitive(iOrig : iOrig + nOrb1*nOrb2 - 1,1) &
              &+ 0.5_dp*reshape(real(tmpSqr(1:nOrb2,1:nOrb1)), (/nOrb1*nOrb2/))
          primitive(iOrig : iOrig + nOrb1*nOrb2 - 1,4) = &
              &primitive(iOrig : iOrig + nOrb1*nOrb2 - 1,4) &
              & + real(1-2*iBlock,dp) * &
              & 0.5_dp*reshape(real(tmpSqr(1:nOrb2,1:nOrb1)), (/nOrb1*nOrb2/))
        end do
      end do
    end do

    do iAtom1 = 1, nAtom
      ii = iAtomStart(iAtom1)
      nOrb1 = iAtomStart(iAtom1+1) - ii
      do iNeigh = 0, nNeighbor(iAtom1)
        iOrig = iPair(iNeigh,iAtom1) + 1
        iAtom2 = iNeighbor(iNeigh, iAtom1)
        iAtom2f = img2CentCell(iAtom2)
        jj = iAtomStart(iAtom2f)
        @:ASSERT(jj >= ii)
        nOrb2 = iAtomStart(iAtom2f+1) - jj
        ! take the Hermitian part of the block
        tmpSqr(1:nOrb2, 1:nOrb1) = &
            & 0.5_dp * (square(jj+nOrb:jj+nOrb2-1+nOrb, &
            & ii:ii+nOrb1-1) &
            & + transpose(square(nOrb+ii:nOrb+ii+nOrb1-1, &
            & jj:jj+nOrb2-1)))
        @:ASSERT(sizePrim >= iOrig + nOrb1*nOrb2 - 1)
        primitive(iOrig : iOrig + nOrb1*nOrb2 - 1,2) = &
            &primitive(iOrig : iOrig + nOrb1*nOrb2 - 1,2) &
            &+ reshape(real(tmpSqr(1:nOrb2,1:nOrb1)), (/nOrb1*nOrb2/))
        primitive(iOrig : iOrig + nOrb1*nOrb2 - 1,3) = &
            &primitive(iOrig : iOrig + nOrb1*nOrb2 - 1,3) &
            & +reshape(aimag(tmpSqr(1:nOrb2,1:nOrb1)), (/nOrb1*nOrb2/))
      end do
    end do

  end subroutine packHSPauli

  !!* Pack squared matrix into the sparse form (complex Pauli version).
  !!* @param primitive    Sparse matrix
  !!* @param square       Squared form matrix
  !!* @param iNeighbor    Neighbor list for the atoms (First index from 0!)
  !!* @param nNeighbor    Nr. of neighbors for the atoms.
  !!* @param mOrb         Maximal number of orbitals on an atom.
  !!* @param iCellVec     Index of the cell translation vector for each atom.
  !!* @param cellVec      Relative coordinates of the cell translation vectors.
  !!* @param iAtomStart   Atom offset for the squared matrix
  !!* @param iPair indexing array for the sparse Hamiltonian
  !!* @param img2CentCell Mapping between image atoms and correspondent atom in
  !!*   the central cell.
  subroutine packHSPauli_kpts(primitive, square, kPoint, kWeight, iNeighbor, &
      & nNeighbor, mOrb, iCellVec, cellVec, iAtomStart, iPair, img2CentCell)
    real(dp), intent(inout)   :: primitive(:,:)
    complex(dp), intent(in) :: square(:,:)
    real(dp),    intent(in) :: kPoint(:)
    real(dp),    intent(in) :: kweight
    integer,  intent(in)    :: iNeighbor(0:,:)
    integer,  intent(in)    :: nNeighbor(:)
    integer,  intent(in)    :: mOrb
    integer,     intent(in) :: iCellVec(:)
    real(dp),    intent(in) :: cellVec(:,:)
    integer,  intent(in)    :: iAtomStart(:)
    integer,  intent(in)    :: iPair(0:,:)
    integer,  intent(in)    :: img2CentCell(:)

    complex(dp) :: phase
    integer     :: nAtom
    integer     :: iOrig, ii, jj, kk
    integer     :: iNeigh, iBlock
    integer     :: iOldVec, iVec
    integer     :: iAtom1, iAtom2, iAtom2f
    integer     :: nOrb1, nOrb2, nOrb
    real(dp)    :: kPoint2p(3)
    complex(dp) :: tmpSqr(mOrb,mOrb)
  #:call ASSERT_CODE
    integer :: sizePrim
  #:endcall ASSERT_CODE

    nAtom = size(iNeighbor, dim=2)
    nOrb = (iAtomStart(nAtom+1) - 1) ! number of orbitals in a regular
                                     ! spin block

  #:call ASSERT_CODE
    sizePrim = size(primitive,dim=1)
  #:endcall ASSERT_CODE

    @:ASSERT(nAtom > 0)
    @:ASSERT(size(square, dim=1) == size(square, dim=2))
    @:ASSERT(size(square, dim=1) == 2 * nOrb )
    @:ASSERT(all(shape(kPoint) == (/ 3 /)))
    @:ASSERT(all(shape(nNeighbor) == (/ nAtom /)))
    @:ASSERT(size(iAtomStart) == nAtom + 1)
    @:ASSERT(kWeight > 0.0_dp)
    @:ASSERT(size(primitive,dim=2)==4)

    kPoint2p(:) = 2.0_dp * pi * kPoint(:)

    ! sigma_I and sigma_z blocks
    do iBlock = 0, 1
      iOldVec = 0
      phase = 1.0_dp
      do iAtom1 = 1, nAtom
        ii = iAtomStart(iAtom1)
        nOrb1 = iAtomStart(iAtom1+1) - ii
        do iNeigh = 0, nNeighbor(iAtom1)
          iOrig = iPair(iNeigh,iAtom1) + 1
          iAtom2 = iNeighbor(iNeigh, iAtom1)
          iAtom2f = img2CentCell(iAtom2)
          jj = iAtomStart(iAtom2f)
          @:ASSERT(jj >= ii)
          nOrb2 = iAtomStart(iAtom2f+1) - jj
          iVec = iCellVec(iAtom2)
          if (iVec /= iOldVec) then
            phase = exp(cmplx(0,-1,dp) &
                &* dot_product(kPoint2p(:), cellVec(:, iVec)))
            iOldVec = iVec
          end if
          tmpSqr(1:nOrb2, 1:nOrb1) = &
              & square(jj+iBlock*nOrb:jj+nOrb2-1+iBlock*nOrb, &
              & ii+iBlock*nOrb:ii+nOrb1-1+iBlock*nOrb)
          ! Hermitian the on-site block before packing, as only one triangle
          ! usually supplied
          if (iAtom1 == iAtom2f) then
            do kk = 1, nOrb2
              tmpSqr(kk, kk+1:nOrb1) = conjg(tmpSqr(kk+1:nOrb1, kk))
            end do
          end if
          @:ASSERT(sizePrim >= iOrig + nOrb1*nOrb2 - 1)
          primitive(iOrig : iOrig + nOrb1*nOrb2 - 1,1) = &
              & primitive(iOrig : iOrig + nOrb1*nOrb2 - 1,1) &
              &+ kWeight * 0.5_dp*reshape(real(phase* &
              & tmpSqr(1:nOrb2,1:nOrb1)), (/nOrb1*nOrb2/))
          primitive(iOrig : iOrig + nOrb1*nOrb2 - 1,4) = &
              & primitive(iOrig : iOrig + nOrb1*nOrb2 - 1,4) &
              & + real(1-2*iBlock,dp) * &
              & kWeight * 0.5_dp*reshape( &
              & real(phase*tmpSqr(1:nOrb2,1:nOrb1)), (/nOrb1*nOrb2/))
        end do
      end do
    end do

    ! sigma_x and sigma_y blocks
    iOldVec = 0
    phase = 1.0_dp
    do iAtom1 = 1, nAtom
      ii = iAtomStart(iAtom1)
      nOrb1 = iAtomStart(iAtom1+1) - ii
      do iNeigh = 0, nNeighbor(iAtom1)
        iOrig = iPair(iNeigh,iAtom1) + 1
        iAtom2 = iNeighbor(iNeigh, iAtom1)
        iAtom2f = img2CentCell(iAtom2)
        jj = iAtomStart(iAtom2f)
        @:ASSERT(jj >= ii)
        nOrb2 = iAtomStart(iAtom2f+1) - jj
        iVec = iCellVec(iAtom2)
        if (iVec /= iOldVec) then
          phase = exp(cmplx(0,-1,dp) &
              &* dot_product(kPoint2p(:), cellVec(:, iVec)))
          iOldVec = iVec
        end if
        @:ASSERT(sizePrim >= iOrig + nOrb1*nOrb2 - 1)
        ! take the Pauli part of the block

        tmpSqr(1:nOrb2, 1:nOrb1) = &
            & 0.5_dp * (square(jj+nOrb:jj+nOrb2-1+nOrb, &
            & ii:ii+nOrb1-1) &
            & + conjg(transpose(square(nOrb+ii:nOrb+ii+nOrb1-1, &
            & jj:jj+nOrb2-1))))
        if (iAtom1 == iAtom2f) then !make up the other side of the on-site block
          do kk = 1, nOrb2
            tmpSqr(kk, kk+1:nOrb1) = tmpSqr(kk+1:nOrb1, kk)
          end do
        end if
        primitive(iOrig : iOrig + nOrb1*nOrb2 - 1,2) = &
            & primitive(iOrig : iOrig + nOrb1*nOrb2 - 1,2) &
            & + kWeight*reshape( &
            & real(phase*tmpSqr(1:nOrb2,1:nOrb1)), (/nOrb1*nOrb2/))
        tmpSqr(1:nOrb2, 1:nOrb1) = &
            & 0.5_dp * (square(jj+nOrb:jj+nOrb2-1+nOrb, &
            & ii:ii+nOrb1-1) &
            & - conjg(transpose(square(nOrb+ii:nOrb+ii+nOrb1-1, &
            & jj:jj+nOrb2-1))))
        if (iAtom1 == iAtom2f) then !make up the other side of the on-site block
          do kk = 1, nOrb2
            tmpSqr(kk, kk+1:nOrb1) = tmpSqr(kk+1:nOrb1, kk)
          end do
        end if
        primitive(iOrig : iOrig + nOrb1*nOrb2 - 1,3) = &
            &primitive(iOrig : iOrig + nOrb1*nOrb2 - 1,3) &
            & + kWeight*reshape( &
            & aimag(phase*tmpSqr(1:nOrb2,1:nOrb1)), (/nOrb1*nOrb2/))
      end do
    end do

  end subroutine packHSPauli_kpts

  !!* Pack imaginary coefficient part of Pauli square matrix into the sparse
  !!* form.
  !!* @param primitive    Sparse matrix
  !!* @param square       Squared form matrix
  !!* @param iNeighbor    Neighbor list for the atoms (First index from 0!)
  !!* @param nNeighbor    Nr. of neighbors for the atoms.
  !!* @param mOrb         Maximal number of orbitals on an atom.
  !!* @param iAtomStart   Atom offset for the squared matrix
  !!* @param iPair indexing array for the sparse Hamiltonian
  !!* @param img2CentCell Mapping between image atoms and correspondent atom in
  !!*   the central cell.
  subroutine packHSPauliImag(primitive, square, iNeighbor, nNeighbor, mOrb, &
      &iAtomStart, iPair, img2CentCell)
    real(dp), intent(inout)   :: primitive(:,:)
    complex(dp), intent(in) :: square(:,:)
    integer,  intent(in)    :: iNeighbor(0:,:)
    integer,  intent(in)    :: nNeighbor(:)
    integer,  intent(in)    :: mOrb
    integer,  intent(in)    :: iAtomStart(:)
    integer,  intent(in)    :: iPair(0:,:)
    integer,  intent(in)    :: img2CentCell(:)

    integer     :: nAtom
    integer     :: iOrig, ii, jj, kk
    integer     :: iNeigh, iBlock
    integer     :: iAtom1, iAtom2, iAtom2f
    integer     :: nOrb1, nOrb2, nOrb
    complex(dp) :: tmpSqr(mOrb,mOrb)
  #:call ASSERT_CODE
    integer :: sizePrim
  #:endcall ASSERT_CODE

    nAtom = size(iNeighbor, dim=2)
    nOrb = (iAtomStart(nAtom+1) - 1) ! number of orbitals in a regular
                                     ! spin block

  #:call ASSERT_CODE
    sizePrim = size(primitive,dim=1)
  #:endcall ASSERT_CODE

    @:ASSERT(nAtom > 0)
    @:ASSERT(size(square, dim=1) == size(square, dim=2))
    @:ASSERT(size(square, dim=1) == 2 * nOrb )
    @:ASSERT(all(shape(nNeighbor) == (/ nAtom /)))
    @:ASSERT(size(iAtomStart) == nAtom + 1)
    @:ASSERT(size(primitive,dim=2)==4)


    do iBlock = 0, 1
      do iAtom1 = 1, nAtom
        ii = iAtomStart(iAtom1)
        nOrb1 = iAtomStart(iAtom1+1) - ii
        do iNeigh = 0, nNeighbor(iAtom1)
          iOrig = iPair(iNeigh,iAtom1) + 1
          iAtom2 = iNeighbor(iNeigh, iAtom1)
          iAtom2f = img2CentCell(iAtom2)
          jj = iAtomStart(iAtom2f)
          @:ASSERT(jj >= ii)
          nOrb2 = iAtomStart(iAtom2f+1) - jj
          tmpSqr(1:nOrb2, 1:nOrb1) = &
              & square(jj+iBlock*nOrb:jj+nOrb2-1+iBlock*nOrb, &
              & ii+iBlock*nOrb:ii+nOrb1-1+iBlock*nOrb)
          ! Hermitian the on-site block before packing, as only one triangle
          ! usually supplied
          if (iAtom1 == iAtom2f) then
            do kk = 1, nOrb2
              tmpSqr(kk, kk+1:nOrb1) = conjg(tmpSqr(kk+1:nOrb1, kk))
            end do
          end if
          @:ASSERT(sizePrim >= iOrig + nOrb1*nOrb2 - 1)
          primitive(iOrig : iOrig + nOrb1*nOrb2 - 1,1) = &
              &primitive(iOrig : iOrig + nOrb1*nOrb2 - 1,1) &
              &+ 0.5_dp*reshape(aimag(tmpSqr(1:nOrb2,1:nOrb1)), (/nOrb1*nOrb2/))
          primitive(iOrig : iOrig + nOrb1*nOrb2 - 1,4) = &
              &primitive(iOrig : iOrig + nOrb1*nOrb2 - 1,4) &
              & + real(1-2*iBlock,dp) * &
              & 0.5_dp*reshape(aimag(tmpSqr(1:nOrb2,1:nOrb1)), (/nOrb1*nOrb2/))
        end do
      end do
    end do

    do iAtom1 = 1, nAtom
      ii = iAtomStart(iAtom1)
      nOrb1 = iAtomStart(iAtom1+1) - ii
      do iNeigh = 0, nNeighbor(iAtom1)
        iOrig = iPair(iNeigh,iAtom1) + 1
        iAtom2 = iNeighbor(iNeigh, iAtom1)
        iAtom2f = img2CentCell(iAtom2)
        jj = iAtomStart(iAtom2f)
        @:ASSERT(jj >= ii)
        nOrb2 = iAtomStart(iAtom2f+1) - jj
        ! take the anti-Hermitian part of the block
        tmpSqr(1:nOrb2, 1:nOrb1) = 0.5_dp * ( &
            & square(nOrb+jj:nOrb+jj+nOrb2-1, ii:ii+nOrb1-1) &
            & - transpose(square(nOrb+ii:nOrb+ii+nOrb1-1,jj:jj+nOrb2-1)))
        @:ASSERT(sizePrim >= iOrig + nOrb1*nOrb2 - 1)
        ! sigma_x : i * 1 = i use + imaginary part of block
        primitive(iOrig : iOrig + nOrb1*nOrb2 - 1,2) = &
            & primitive(iOrig : iOrig + nOrb1*nOrb2 - 1,2) &
            & + reshape(aimag(tmpSqr(1:nOrb2,1:nOrb1)), (/nOrb1*nOrb2/))
        ! sigma_y : i * i = -1 use - real part of block
        primitive(iOrig : iOrig + nOrb1*nOrb2 - 1,3) = &
            & primitive(iOrig : iOrig + nOrb1*nOrb2 - 1,3) &
            & -reshape(real(tmpSqr(1:nOrb2,1:nOrb1)), (/nOrb1*nOrb2/))
      end do
    end do

  end subroutine packHSPauliImag

  !!* Pack imaginary coefficient part of Pauli square matrix into the sparse
  !!* form (complex version).
  !!* @param primitive    Sparse matrix
  !!* @param square       Squared form matrix
  !!* @param kPoint       Relative coordinates of the K-point
  !!* @param kWeight      Weight of the K-point
  !!* @param iNeighbor    Neighbor list for the atoms (First index from 0!)
  !!* @param nNeighbor    Nr. of neighbors for the atoms.
  !!* @param mOrb         Maximal number of orbitals on an atom.
  !!* @param iCellVec     Index of the cell translation vector for each atom.
  !!* @param cellVec      Relative coordinates of the cell translation vectors.
  !!* @param iAtomStart   Atom offset for the squared matrix
  !!* @param iPair indexing array for the sparse Hamiltonian
  !!* @param img2CentCell Mapping between image atoms and correspondent atom in
  !!*   the central cell.
  subroutine packHSPauliImag_kpts(primitive, square, kPoint, kWeight, &
      & iNeighbor, nNeighbor, mOrb, iCellVec, cellVec, iAtomStart, iPair, &
      & img2CentCell)
    real(dp), intent(inout) :: primitive(:,:)
    complex(dp), intent(in) :: square(:,:)
    real(dp),    intent(in) :: kPoint(:)
    real(dp),    intent(in) :: kweight
    integer,  intent(in)    :: iNeighbor(0:,:)
    integer,  intent(in)    :: nNeighbor(:)
    integer,  intent(in)    :: mOrb
    integer,     intent(in) :: iCellVec(:)
    real(dp),    intent(in) :: cellVec(:,:)
    integer,  intent(in)    :: iAtomStart(:)
    integer,  intent(in)    :: iPair(0:,:)
    integer,  intent(in)    :: img2CentCell(:)

    complex(dp) :: phase
    integer     :: nAtom
    integer     :: iOrig, ii, jj, kk
    integer     :: iNeigh, iBlock
    integer     :: iOldVec, iVec
    integer     :: iAtom1, iAtom2, iAtom2f
    integer     :: nOrb1, nOrb2, nOrb
    real(dp)    :: kPoint2p(3)
    complex(dp) :: tmpSqr(mOrb,mOrb)
  #:call ASSERT_CODE
    integer :: sizePrim
  #:endcall ASSERT_CODE

    nAtom = size(iNeighbor, dim=2)
    nOrb = (iAtomStart(nAtom+1) - 1) ! number of orbitals in a regular
                                     ! spin block

  #:call ASSERT_CODE
    sizePrim = size(primitive,dim=1)
  #:endcall ASSERT_CODE

    @:ASSERT(nAtom > 0)
    @:ASSERT(size(square, dim=1) == size(square, dim=2))
    @:ASSERT(size(square, dim=1) == 2 * nOrb )
    @:ASSERT(all(shape(kPoint) == (/ 3 /)))
    @:ASSERT(all(shape(nNeighbor) == (/ nAtom /)))
    @:ASSERT(size(iAtomStart) == nAtom + 1)
    @:ASSERT(kWeight > 0.0_dp)
    @:ASSERT(size(primitive,dim=2)==4)

    kPoint2p(:) = 2.0_dp * pi * kPoint(:)

    ! sigma_I and sigma_z blocks
    do iBlock = 0, 1
      iOldVec = 0
      phase = 1.0_dp
      do iAtom1 = 1, nAtom
        ii = iAtomStart(iAtom1)
        nOrb1 = iAtomStart(iAtom1+1) - ii
        do iNeigh = 0, nNeighbor(iAtom1)
          iOrig = iPair(iNeigh,iAtom1) + 1
          iAtom2 = iNeighbor(iNeigh, iAtom1)
          iAtom2f = img2CentCell(iAtom2)
          jj = iAtomStart(iAtom2f)
          @:ASSERT(jj >= ii)
          nOrb2 = iAtomStart(iAtom2f+1) - jj
          iVec = iCellVec(iAtom2)
          if (iVec /= iOldVec) then
            phase = exp(cmplx(0,-1,dp) &
                &* dot_product(kPoint2p(:), cellVec(:, iVec)))
            iOldVec = iVec
          end if
          tmpSqr(1:nOrb2, 1:nOrb1) = &
              & square(jj+iBlock*nOrb:jj+nOrb2-1+iBlock*nOrb, &
              & ii+iBlock*nOrb:ii+nOrb1-1+iBlock*nOrb)
          ! Hermitian the on-site block before packing, as only one triangle
          ! usually supplied
          if (iAtom1 == iAtom2f) then
            do kk = 1, nOrb2
              tmpSqr(kk, kk+1:nOrb1) = conjg(tmpSqr(kk+1:nOrb1, kk))
            end do
          end if
          tmpSqr = phase*kWeight*tmpSqr
          @:ASSERT(sizePrim >= iOrig + nOrb1*nOrb2 - 1)
          primitive(iOrig : iOrig + nOrb1*nOrb2 - 1,1) = &
              &primitive(iOrig : iOrig + nOrb1*nOrb2 - 1,1) &
              &+ 0.5_dp*reshape(aimag(tmpSqr(1:nOrb2,1:nOrb1)), (/nOrb1*nOrb2/))
          primitive(iOrig : iOrig + nOrb1*nOrb2 - 1,4) = &
              &primitive(iOrig : iOrig + nOrb1*nOrb2 - 1,4) &
              & + real(1-2*iBlock,dp) * &
              & 0.5_dp*reshape(aimag(tmpSqr(1:nOrb2,1:nOrb1)), (/nOrb1*nOrb2/))
        end do
      end do
    end do

    ! sigma_x and sigma_y blocks
    iOldVec = 0
    phase = 1.0_dp
    do iAtom1 = 1, nAtom
      ii = iAtomStart(iAtom1)
      nOrb1 = iAtomStart(iAtom1+1) - ii
      do iNeigh = 0, nNeighbor(iAtom1)
        iOrig = iPair(iNeigh,iAtom1) + 1
        iAtom2 = iNeighbor(iNeigh, iAtom1)
        iAtom2f = img2CentCell(iAtom2)
        jj = iAtomStart(iAtom2f)
        @:ASSERT(jj >= ii)
        nOrb2 = iAtomStart(iAtom2f+1) - jj
        iVec = iCellVec(iAtom2)
        if (iVec /= iOldVec) then
          phase = exp(cmplx(0,-1,dp)*dot_product(kPoint2p(:), cellVec(:, iVec)))
          iOldVec = iVec
        end if
        @:ASSERT(sizePrim >= iOrig + nOrb1*nOrb2 - 1)

        ! take the anti-Hermitian part of the block
        tmpSqr(1:nOrb2, 1:nOrb1) = 0.5_dp * ( &
            & square(nOrb+jj:nOrb+jj+nOrb2-1, ii:ii+nOrb1-1) &
            & + conjg(transpose(square(nOrb+ii:nOrb+ii+nOrb1-1,jj:jj+nOrb2-1))))
        tmpSqr = phase*kWeight*tmpSqr
        ! sigma_x : i * 1 = i use + imaginary part of block
        primitive(iOrig : iOrig + nOrb1*nOrb2 - 1,2) = &
            & primitive(iOrig : iOrig + nOrb1*nOrb2 - 1,2) &
            & +reshape(aimag(tmpSqr(1:nOrb2,1:nOrb1)), (/nOrb1*nOrb2/))

        ! take the anti-Hermitian part of the block
        tmpSqr(1:nOrb2, 1:nOrb1) = 0.5_dp * ( &
            & square(nOrb+jj:nOrb+jj+nOrb2-1, ii:ii+nOrb1-1) &
            & - conjg(transpose(square(nOrb+ii:nOrb+ii+nOrb1-1,jj:jj+nOrb2-1))))
        tmpSqr = phase*kWeight*tmpSqr
        ! sigma_y : i * i = -1 use - real part of block
        primitive(iOrig : iOrig + nOrb1*nOrb2 - 1,3) = &
            & primitive(iOrig : iOrig + nOrb1*nOrb2 - 1,3) &
            & -reshape(real(tmpSqr(1:nOrb2,1:nOrb1)), (/nOrb1*nOrb2/))



      end do
    end do

  end subroutine packHSPauliImag_kpts

  !!* Pack only the charge (spin channel 1) part of a 2 component matrix
  !!* @param primitive    Sparse matrix
  !!* @param square       Squared form matrix
  !!* @param iNeighbor    Neighbor list for the atoms (First index from 0!)
  !!* @param nNeighbor    Nr. of neighbors for the atoms.
  !!* @param mOrb         Maximal number of orbitals on an atom.
  !!* @param iCellVec     Index of the cell translation vector for each atom.
  !!* @param cellVec      Relative coordinates of the cell translation vectors.
  !!* @param iAtomStart   Atom offset for the squared matrix
  !!* @param iPair indexing array for the sparse Hamiltonian
  !!* @param img2CentCell Mapping between image atoms and correspondent atom in
  !!*   the central cell.
  subroutine packHSPauliERho(primitive, square, iNeighbor, nNeighbor, mOrb, &
      &iAtomStart, iPair, img2CentCell)
    real(dp), intent(inout)   :: primitive(:)
    complex(dp), intent(in) :: square(:,:)
    integer,  intent(in)    :: iNeighbor(0:,:)
    integer,  intent(in)    :: nNeighbor(:)
    integer,  intent(in)    :: mOrb
    integer,  intent(in)    :: iAtomStart(:)
    integer,  intent(in)    :: iPair(0:,:)
    integer,  intent(in)    :: img2CentCell(:)

    integer     :: nAtom
    integer     :: iOrig, ii, jj, kk
    integer     :: iNeigh, iBlock
    integer     :: iAtom1, iAtom2, iAtom2f
    integer     :: nOrb1, nOrb2, nOrb
    complex(dp) :: tmpSqr(mOrb,mOrb)
  #:call ASSERT_CODE
    integer :: sizePrim
  #:endcall ASSERT_CODE

    nAtom = size(iNeighbor, dim=2)
    nOrb = (iAtomStart(nAtom+1) - 1) ! number of orbitals in a regular spin block

  #:call ASSERT_CODE
    sizePrim = size(primitive,dim=1)
  #:endcall ASSERT_CODE

    @:ASSERT(nAtom > 0)
    @:ASSERT(size(square, dim=1) == size(square, dim=2))
    @:ASSERT(size(square, dim=1) == 2 * nOrb )
    @:ASSERT(all(shape(nNeighbor) == (/ nAtom /)))
    @:ASSERT(size(iAtomStart) == nAtom + 1)

    do iBlock = 0, 1
      do iAtom1 = 1, nAtom
        ii = iAtomStart(iAtom1)
        nOrb1 = iAtomStart(iAtom1+1) - ii
        do iNeigh = 0, nNeighbor(iAtom1)
          iOrig = iPair(iNeigh,iAtom1) + 1
          iAtom2 = iNeighbor(iNeigh, iAtom1)
          iAtom2f = img2CentCell(iAtom2)
          jj = iAtomStart(iAtom2f)
          @:ASSERT(jj >= ii)
          nOrb2 = iAtomStart(iAtom2f+1) - jj
          tmpSqr(1:nOrb2, 1:nOrb1) = &
              & square(jj+iBlock*nOrb:jj+nOrb2-1+iBlock*nOrb, &
              & ii+iBlock*nOrb:ii+nOrb1-1+iBlock*nOrb)
          ! Symmetrize the on-site block before packing, as only one
          ! triangle supplied
          if (iAtom1 == iAtom2f) then
            do kk = 1, nOrb2
              tmpSqr(kk, kk+1:nOrb1) = conjg(tmpSqr(kk+1:nOrb1, kk))
            end do
          end if
          @:ASSERT(sizePrim >= iOrig + nOrb1*nOrb2 - 1)
          primitive(iOrig : iOrig + nOrb1*nOrb2 - 1) = &
              &primitive(iOrig : iOrig + nOrb1*nOrb2 - 1) &
              & + reshape(real(tmpSqr(1:nOrb2,1:nOrb1)), (/nOrb1*nOrb2/))
        end do
      end do
    end do

  end subroutine packHSPauliERho

  !!* Pack squared matrix in the sparse form (real version).
  !!* @param primitive    Sparse matrix
  !!* @param square       Squared form matrix
  !!* @param iNeighbor    Neighbor list for the atoms (First index from 0!)
  !!* @param nNeighbor    Nr. of neighbors for the atoms.
  !!* @param mOrb         Maximal number of orbitals on an atom.
  !!* @param iCellVec     Index of the cell translation vector for each atom.
  !!* @param cellVec      Relative coordinates of the cell translation vectors.
  !!* @param iAtomStart   Atom offset for the squared matrix
  !!* @param iPair indexing array for the sparse Hamiltonian
  !!* @param img2CentCell Mapping between image atoms and correspondent atom in
  !!*   the central cell.
  subroutine packHSPauliERho_kpts(primitive, square, kPoint, kWeight, &
      & iNeighbor, nNeighbor, mOrb, iCellVec, cellVec, iAtomStart, iPair, &
      & img2CentCell)
    real(dp), intent(inout)   :: primitive(:)
    complex(dp), intent(in) :: square(:,:)
    real(dp),    intent(in) :: kPoint(:)
    real(dp),    intent(in) :: kweight
    integer,  intent(in)    :: iNeighbor(0:,:)
    integer,  intent(in)    :: nNeighbor(:)
    integer,  intent(in)    :: mOrb
    integer,     intent(in) :: iCellVec(:)
    real(dp),    intent(in) :: cellVec(:,:)
    integer,  intent(in)    :: iAtomStart(:)
    integer,  intent(in)    :: iPair(0:,:)
    integer,  intent(in)    :: img2CentCell(:)

    complex(dp) :: phase
    integer     :: nAtom
    integer     :: iOrig, ii, jj, kk
    integer     :: iNeigh, iBlock
    integer     :: iOldVec, iVec
    integer     :: iAtom1, iAtom2, iAtom2f
    integer     :: nOrb1, nOrb2, nOrb
    real(dp)    :: kPoint2p(3)
    complex(dp) :: tmpSqr(mOrb,mOrb)
  #:call ASSERT_CODE
    integer :: sizePrim
  #:endcall ASSERT_CODE

    nAtom = size(iNeighbor, dim=2)
    nOrb = (iAtomStart(nAtom+1) - 1) ! number of orbitals in a regular
                                     ! spin block
  #:call ASSERT_CODE
    sizePrim = size(primitive,dim=1)
  #:endcall ASSERT_CODE

    @:ASSERT(nAtom > 0)
    @:ASSERT(size(square, dim=1) == size(square, dim=2))
    @:ASSERT(size(square, dim=1) == 2 * nOrb )
    @:ASSERT(all(shape(kPoint) == (/ 3 /)))
    @:ASSERT(all(shape(nNeighbor) == (/ nAtom /)))
    @:ASSERT(size(iAtomStart) == nAtom + 1)
    @:ASSERT(kWeight > 0.0_dp)

    kPoint2p(:) = 2.0_dp * pi * kPoint(:)

    do iBlock = 0, 1
      iOldVec = 0
      phase = 1.0_dp
      do iAtom1 = 1, nAtom
        ii = iAtomStart(iAtom1)
        nOrb1 = iAtomStart(iAtom1+1) - ii
        do iNeigh = 0, nNeighbor(iAtom1)
          iOrig = iPair(iNeigh,iAtom1) + 1
          iAtom2 = iNeighbor(iNeigh, iAtom1)
          iAtom2f = img2CentCell(iAtom2)
          jj = iAtomStart(iAtom2f)
          @:ASSERT(jj >= ii)
          nOrb2 = iAtomStart(iAtom2f+1) - jj
          iVec = iCellVec(iAtom2)
          if (iVec /= iOldVec) then
            phase = exp(cmplx(0,-1,dp) &
                &* dot_product(kPoint2p(:), cellVec(:, iVec)))
            iOldVec = iVec
          end if
          tmpSqr(1:nOrb2, 1:nOrb1) = &
              & square(jj+iBlock*nOrb:jj+nOrb2-1+iBlock*nOrb, &
              & ii+iBlock*nOrb:ii+nOrb1-1+iBlock*nOrb)
          ! Symmetrize the on-site block before packing, as only one
          ! triangle supplied
          if (iAtom1 == iAtom2f) then
            do kk = 1, nOrb2
              tmpSqr(kk, kk+1:nOrb1) = conjg(tmpSqr(kk+1:nOrb1, kk))
            end do
          end if
          @:ASSERT(sizePrim >= iOrig + nOrb1*nOrb2 - 1)
          primitive(iOrig : iOrig + nOrb1*nOrb2 - 1) = &
              & primitive(iOrig : iOrig + nOrb1*nOrb2 - 1) + kWeight * &
              & reshape(real(phase*tmpSqr(1:nOrb2,1:nOrb1)), (/nOrb1*nOrb2/))
        end do
      end do
    end do

  end subroutine packHSPauliERho_kpts

  !!* Symmetrize a squared matrix leaving the on-site atomic blocks alone.
  !!* (Complex version)
  !!* @param square   Square form matrix.
  !!* @param iAtomStart  Returns the offset array for each atom.
  subroutine blockSymmetrizeHS_cmplx(square, iAtomStart)
    complex(dp), intent(inout) :: square(:,:)
    integer, intent(in) :: iAtomStart(:)

    integer     :: nAtom, iAtom, iStart, iEnd, mOrb

    nAtom = size(iAtomStart, dim=1) - 1
    mOrb = iAtomStart(nAtom+1) - 1

    @:ASSERT(nAtom > 0)
    @:ASSERT(size(square, dim=1) == size(square, dim=2))
    @:ASSERT((size(square, dim=1) == 2*mOrb) .or. (size(square, dim=1) == mOrb))
    @:ASSERT(size(square, dim=1) == mOrb)

    do iAtom = 1, nAtom
      iStart = iAtomStart(iAtom)
      iEnd = iAtomStart(iAtom+1) - 1
      square(iStart:iEnd, iEnd+1:mOrb) = &
          & transpose(square(iEnd+1:mOrb,iStart:iEnd))
    end do

  end subroutine blockSymmetrizeHS_cmplx


  !!* Symmetrize a squared matrix leaving the on-site atomic blocks alone.
  !!* (Complex version)
  !!* @param square   Square form matrix.
  !!* @param iAtomStart  Returns the offset array for each atom.
  subroutine blockHermitianHS_cmplx(square, iAtomStart)
    complex(dp), intent(inout) :: square(:,:)
    integer, intent(in) :: iAtomStart(:)

    integer     :: nAtom, iAtom, iStart, iEnd, mOrb

    nAtom = size(iAtomStart, dim=1) - 1
    mOrb = iAtomStart(nAtom+1) - 1

    @:ASSERT(nAtom > 0)
    @:ASSERT(size(square, dim=1) == size(square, dim=2))
    @:ASSERT((size(square, dim=1) == 2*mOrb) .or. (size(square, dim=1) == mOrb))
    @:ASSERT(size(square, dim=1) == mOrb)

    do iAtom = 1, nAtom
      iStart = iAtomStart(iAtom)
      iEnd = iAtomStart(iAtom+1) - 1
      square(iStart:iEnd, iEnd+1:mOrb) = &
          & transpose(conjg(square(iEnd+1:mOrb,iStart:iEnd)))
    end do

  end subroutine blockHermitianHS_cmplx

  !!* Symmetrize a squared matrix leaving the on-site atomic blocks alone.
  !!* (Real version)
  !!* @param square   Square form matrix.
  !!* @param iAtomStart  Returns the offset array for each atom.
  subroutine blockSymmetrizeHS_real(square, iAtomStart)
    real(dp), intent(inout) :: square(:,:)
    integer,  intent(in)    :: iAtomStart(:)

    integer     :: nAtom, iAtom, iStart, iEnd, mOrb

    nAtom = size(iAtomStart) - 1
    mOrb = iAtomStart(nAtom+1) - 1

    @:ASSERT(nAtom > 0)
    @:ASSERT(size(square, dim=1) == size(square, dim=2))
    @:ASSERT(size(square, dim=1) == mOrb)

    do iAtom = 1, nAtom
      iStart = iAtomStart(iAtom)
      iEnd = iAtomStart(iAtom+1) - 1
      square(iStart:iEnd, iEnd+1:mOrb) = &
          & transpose(square(iEnd+1:mOrb,iStart:iEnd))
    end do

  end subroutine blockSymmetrizeHS_real


  !!* Anti-symmetrize a squared matrix leaving the on-site atomic blocks alone.
  !!* (Real version)
  !!* @param square   Square form matrix.
  !!* @param iAtomStart  Contains the offset in the array for each atom.
  subroutine blockAntiSymmetrizeHS_real(square, iAtomStart)
    real(dp), intent(inout) :: square(:,:)
    integer,  intent(in)    :: iAtomStart(:)

    integer     :: nAtom, iAtom, iStart, iEnd, mOrb

    nAtom = size(iAtomStart) - 1
    mOrb = iAtomStart(nAtom+1) - 1

    @:ASSERT(nAtom > 0)
    @:ASSERT(size(square, dim=1) == size(square, dim=2))
    @:ASSERT(size(square, dim=1) == mOrb)

    do iAtom = 1, nAtom
      iStart = iAtomStart(iAtom)
      iEnd = iAtomStart(iAtom+1) - 1
      square(iStart:iEnd, iEnd+1:mOrb) = &
          & -transpose(square(iEnd+1:mOrb,iStart:iEnd))
    end do

  end subroutine blockAntiSymmetrizeHS_real

end module sparse2dense
