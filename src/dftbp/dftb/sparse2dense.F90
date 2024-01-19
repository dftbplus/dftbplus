!--------------------------------------------------------------------------------------------------!
!  DFTB+: general package for performing fast atomistic simulations                                !
!  Copyright (C) 2006 - 2023  DFTB+ developers group                                               !
!                                                                                                  !
!  See the LICENSE file for terms of usage and distribution.                                       !
!--------------------------------------------------------------------------------------------------!

#:include 'common.fypp'

!> Contains subroutines for packing/unpacking Hamiltonian-like matrices between the square and
!> 1-dimensional representations
!>
module dftbp_dftb_sparse2dense
  use dftbp_common_accuracy, only : dp
  use dftbp_common_constants, only : pi, imag
  use dftbp_dftb_periodic, only : TNeighbourList
  use dftbp_math_angmomentum, only : rotateZ
  use dftbp_math_matrixops, only : adjointLowerTriangle
  use dftbp_type_commontypes, only : TOrbitals
  use dftbp_type_densedescr, only : TDenseDescr
#:if WITH_SCALAPACK
  use dftbp_common_blacsenv, only : TBlacsEnv
  use dftbp_extlibs_scalapackfx, only : scalafx_cpg2l, scalafx_addl2g
#:endif
  implicit none

  private
  public :: unpackHS, packHS, iPackHS, packErho, unpackDQ
  public :: packHSPauli, packHSPauliImag, unpackHPauli, unpackSPauli
  public :: unpackHelicalHS, packHelicalHS
  public :: getSparseDescriptor

#:if WITH_SCALAPACK
  public :: unpackHSRealBlacs, unpackHSCplxBlacs, unpackHPauliBlacs, unpackSPauliBlacs
  public :: packRhoRealBlacs, packRhoCplxBlacs, packRhoPauliBlacs, packERhoPauliBlacs
  public :: unpackHSHelicalRealBlacs, unpackHSHelicalCplxBlacs
  public :: packRhoHelicalRealBlacs, packRhoHelicalCplxBlacs
#:endif


  !> Unpack sparse matrix (Hamiltonian, overlap, etc.) to square form
  interface unpackHS
    module procedure unpackHS_real
    module procedure unpackHS_cmplx_kpts
  end interface unpackHS


  !> Unpack sparse matrix (dipole, quadrupole integrals) to square form
  interface unpackDQ
    module procedure unpackDQ_real
    module procedure unpackDQ_cmplx_kpts
  end interface unpackDQ


  !> Unpack sparse matrix (Hamiltonian, overlap, etc.) to square form for helical geometries
  interface unpackHelicalHS
    module procedure unpackHSHelical_real
    module procedure unpackHSHelical_cmplx
  end interface unpackHelicalHS


  !> Pack square matrix to sparse form
  interface packHS
    module procedure packHS_real
    module procedure packHS_cmplx_kpts
    module procedure packHSPauli
    module procedure packHSPauli_kpts
    module procedure packhs_cmplx
  end interface packHS


  !> Pack square matrix to sparse form for helical geometries
  interface packHelicalHS
    module procedure packHShelical_real
    module procedure packHShelical_cmplx
  end interface packHelicalHS


  !> Pack square matrix to sparse form.
  interface iPackHS
    module procedure packHSPauliImag
    module procedure packHSPauliImag_kpts
  end interface iPackHS


  !> Pack energy weighted Pauli idenity square matrix to sparse form.
  interface packErho
    module procedure packHSPauliERho
    module procedure packHSPauliERho_kpts
  end interface packErho


contains

  !> Unpacks sparse matrix to square form (complex version) Note the non on-site blocks are only
  !> filled in the lower triangle part of the matrix.
  subroutine unpackHS_cmplx_kpts(square, orig, kPoint, iNeighbour, nNeighbourSK, iCellVec, cellVec,&
      & iAtomStart, iSparseStart, img2CentCell)

    !> Square form matrix on exit.
    complex(dp), intent(out) :: square(:, :)

    !> Sparse matrix
    real(dp), intent(in) :: orig(:)

    !> Relative coordinates of the K-point where the sparse matrix should be unfolded.
    real(dp), intent(in) :: kPoint(:)

    !> Neighbour list for each atom (First index from 0!)
    integer, intent(in) :: iNeighbour(0:, :)

    !> Nr. of neighbours for each atom (incl. itself).
    integer, intent(in) :: nNeighbourSK(:)

    !> Index of the cell translation vector for each atom.
    integer, intent(in) :: iCellVec(:)

    !> Relative coordinates of the cell translation vectors.
    real(dp), intent(in) :: cellVec(:, :)

    !> Atom offset for the square Hamiltonian
    integer, intent(in) :: iAtomStart(:)

    !> indexing array for the sparse Hamiltonian
    integer, intent(in) :: iSparseStart(0:, :)

    !> Map from images of atoms to central cell atoms
    integer, intent(in) :: img2CentCell(:)

    complex(dp) :: phase
    integer :: nAtom
    integer :: iOrig, ii, jj
    integer :: iNeigh
    integer :: iOldVec, iVec
    integer :: iAtom1, iAtom2, iAtom2f
    integer :: nOrb1, nOrb2
    real(dp) :: kPoint2p(3)

    nAtom = size(iNeighbour, dim=2)

    @:ASSERT(nAtom > 0)
    @:ASSERT(size(square, dim=1) == size(square, dim=2))
    @:ASSERT(size(square, dim=1) == iAtomStart(nAtom+1) - 1)
    @:ASSERT(all(shape(kPoint) == [3]))
    @:ASSERT(all(shape(nNeighbourSK) == [nAtom]))
    @:ASSERT(size(iAtomStart) == nAtom + 1)

    square(:, :) = cmplx(0, 0, dp)
    kPoint2p(:) = 2.0_dp * pi * kPoint
    iOldVec = 0
    phase = 1.0_dp
    do iAtom1 = 1, nAtom
      ii = iAtomStart(iAtom1)
      nOrb1 = iAtomStart(iAtom1 + 1) - ii
      do iNeigh = 0, nNeighbourSK(iAtom1)
        iOrig = iSparseStart(iNeigh, iAtom1) + 1
        iAtom2 = iNeighbour(iNeigh, iAtom1)
        iAtom2f = img2CentCell(iAtom2)
        jj = iAtomStart(iAtom2f)
        @:ASSERT(jj >= ii)
        nOrb2 = iAtomStart(iAtom2f + 1) - jj
        iVec = iCellVec(iAtom2)
        if (iVec /= iOldVec) then
          phase = exp((0.0_dp, 1.0_dp) * dot_product(kPoint2p, cellVec(:, iVec)))
          iOldVec = iVec
        end if
        square(jj:jj+nOrb2-1, ii:ii+nOrb1-1) = square(jj:jj+nOrb2-1, ii:ii+nOrb1-1)&
            & + phase * reshape(orig(iOrig:iOrig+nOrb1*nOrb2-1), [nOrb2, nOrb1])
      end do
    end do

  end subroutine unpackHS_cmplx_kpts


  !> Unpacks sparse matrix to square form (real version for Gamma point)
  !!
  !! Note: The non on-site blocks are only filled in the lower triangle part of the matrix.
  subroutine unpackHS_real(square, orig, iNeighbour, nNeighbourSK, iAtomStart, iSparseStart,&
      & img2CentCell)

    !> Square form matrix on exit.
    real(dp), intent(out) :: square(:,:)

    !> Sparse matrix
    real(dp), intent(in) :: orig(:)

    !> Neighbour list for each atom (First index from 0!)
    integer, intent(in) :: iNeighbour(0:,:)

    !> Nr. of neighbours for each atom (incl. itself).
    integer, intent(in) :: nNeighbourSK(:)

    !> Atom offset for the square Hamiltonian
    integer, intent(in) :: iAtomStart(:)

    !> indexing array for the sparse Hamiltonian
    integer, intent(in) :: iSparseStart(0:,:)

    !> Map from images of atoms to central cell atoms
    integer, intent(in) :: img2CentCell(:)

    integer :: nAtom
    integer :: iOrig, ii, jj
    integer :: iNeigh
    integer :: iAtom1, iAtom2, iAtom2f
    integer :: nOrb1, nOrb2

    nAtom = size(iNeighbour, dim=2)

    @:ASSERT(nAtom > 0)
    @:ASSERT(size(square, dim=1) == size(square, dim=2))
    @:ASSERT(size(square, dim=1) == iAtomStart(nAtom+1) - 1)
    @:ASSERT(all(shape(nNeighbourSK) == [nAtom]))
    @:ASSERT(size(iAtomStart) == nAtom + 1)

    square(:,:) = 0.0_dp

    do iAtom1 = 1, nAtom
      ii = iAtomStart(iAtom1)
      nOrb1 = iAtomStart(iAtom1 + 1) - ii
      do iNeigh = 0, nNeighbourSK(iAtom1)
        iOrig = iSparseStart(iNeigh, iAtom1) + 1
        iAtom2 = iNeighbour(iNeigh, iAtom1)
        iAtom2f = img2CentCell(iAtom2)
        jj = iAtomStart(iAtom2f)
        @:ASSERT(jj >= ii)
        nOrb2 = iAtomStart(iAtom2f + 1) - jj
        square(jj:jj+nOrb2-1, ii:ii+nOrb1-1) = square(jj:jj+nOrb2-1, ii:ii+nOrb1-1)&
            & + reshape(orig(iOrig:iOrig+nOrb1*nOrb2-1), [nOrb2, nOrb1])
      end do
    end do

  end subroutine unpackHS_real


  !> Unpacks sparse matrix to square form (complex version) for helical geometries. Note the non
  !> on-site blocks are only filled in the lower triangle part of the matrix.
  subroutine unpackHSHelical_cmplx(square, orig, kPoint, iNeighbour, nNeighbourSK, iCellVec,&
      & cellVec, iAtomStart, iSparseStart, img2CentCell, orb, species, coord)

    !> Square form matrix on exit.
    complex(dp), intent(out) :: square(:,:)

    !> Sparse matrix
    real(dp), intent(in) :: orig(:)

    !> Relative coordinates of the K-point where the sparse matrix should be unfolded.
    real(dp), intent(in) :: kPoint(:)

    !> Neighbour list for each atom (First index from 0!)
    integer, intent(in) :: iNeighbour(0:, :)

    !> Nr. of neighbours for each atom (incl. itself).
    integer, intent(in) :: nNeighbourSK(:)

    !> Index of the cell translation vector for each atom.
    integer, intent(in) :: iCellVec(:)

    !> Relative coordinates of the cell translation vectors.
    real(dp), intent(in) :: cellVec(:,:)

    !> Atom offset for the square Hamiltonian
    integer, intent(in) :: iAtomStart(:)

    !> indexing array for the sparse Hamiltonian
    integer, intent(in) :: iSparseStart(0:, :)

    !> Map from images of atoms to central cell atoms
    integer, intent(in) :: img2CentCell(:)

    !> data type for atomic orbital information
    type(TOrbitals), intent(in) :: orb

    !> Species of each atom
    integer :: species(:)

    !> Coordinates of all atoms
    real(dp), intent(in) :: coord(:,:)

    complex(dp) :: phase
    integer :: nAtom, iOrig, ii, jj, iNeigh, iOldVec, iVec, iAtom1, iAtom2, iAtom2f
    integer :: nOrb1, nOrb2, iSh, iSp
    real(dp) :: kPoint2p(2), rotZ(orb%mOrb,orb%mOrb), theta, tmpSqr(orb%mOrb,orb%mOrb)
    integer :: lShellVals(orb%mShell)

    nAtom = size(iNeighbour, dim=2)
    square(:, :) = cmplx(0, 0, dp)
    kPoint2p(:) = 2.0_dp * pi * kPoint
    iOldVec = 0
    phase = 1.0_dp
    lShellVals(:) = 0
    rotZ(:,:) = 0.0_dp
    do iAtom1 = 1, nAtom
      ii = iAtomStart(iAtom1)
      nOrb1 = iAtomStart(iAtom1+1) - ii
      do iNeigh = 0, nNeighbourSK(iAtom1)
        iOrig = iSparseStart(iNeigh,iAtom1) + 1
        iAtom2 = iNeighbour(iNeigh, iAtom1)
        iAtom2f = img2CentCell(iAtom2)
        jj = iAtomStart(iAtom2f)
        nOrb2 = iAtomStart(iAtom2f+1) - jj
        iVec = iCellVec(iAtom2)
        if (iVec /= iOldVec) then
          phase = exp((0.0_dp, 1.0_dp) * dot_product(kPoint2p(:2), cellVec(:2, iVec)) )
          iOldVec = iVec
        end if
        tmpSqr(:nOrb2,:nOrb1) = reshape(orig(iOrig:iOrig+nOrb1*nOrb2-1), (/nOrb2, nOrb1/))
        iSp = species(iAtom2f)
        iSh = orb%nShell(iSp)
        lShellVals(:iSh) = orb%angShell(:iSh,iSp)
        theta = -atan2(coord(2,iAtom2),coord(1,iAtom2))&
            & + atan2(coord(2,iAtom2f),coord(1,iAtom2f))
        theta = mod(theta,2.0_dp*pi)
        call rotateZ(rotZ, lShellVals(:iSh), theta)
        tmpSqr(:nOrb2,:nOrb1) = matmul(rotZ(:nOrb2,:nOrb2),tmpSqr(:nOrb2,:nOrb1))
        square(jj:jj+nOrb2-1, ii:ii+nOrb1-1) = square(jj:jj+nOrb2-1, ii:ii+nOrb1-1)&
            & + phase * tmpSqr(:nOrb2,:nOrb1)
      end do
    end do

  end subroutine unpackHSHelical_cmplx


  !> Unpacks sparse matrix to square form (real version for Gamma point) for helical geometry
  !>
  !> Note: The non on-site blocks are only filled in the lower triangle part of the matrix.
  subroutine unpackHSHelical_real(square, orig, iNeighbour, nNeighbourSK, iAtomStart, iSparseStart,&
      & img2CentCell, orb, species, coord)

    !> Square form matrix on exit.
    real(dp), intent(out) :: square(:, :)

    !> Sparse matrix
    real(dp), intent(in) :: orig(:)

    !> Neighbour list for each atom (First index from 0!)
    integer, intent(in) :: iNeighbour(0:, :)

    !> Nr. of neighbours for each atom (incl. itself).
    integer, intent(in) :: nNeighbourSK(:)

    !> Atom offset for the square Hamiltonian
    integer, intent(in) :: iAtomStart(:)

    !> indexing array for the sparse Hamiltonian
    integer, intent(in) :: iSparseStart(0:, :)

    !> Map from images of atoms to central cell atoms
    integer, intent(in) :: img2CentCell(:)

    !> data type for atomic orbital information
    type(TOrbitals), intent(in) :: orb

    !> Species of each atom
    integer :: species(:)

    !> Coordinates of all atoms
    real(dp), intent(in) :: coord(:,:)

    integer :: nAtom, iOrig, ii, jj, iNeigh, iAtom1, iAtom2, iAtom2f, nOrb1, nOrb2
    real(dp) :: rotZ(orb%mOrb,orb%mOrb), theta, tmpSqr(orb%mOrb,orb%mOrb)
    integer :: lShellVals(orb%mShell), iSh, iSp

    nAtom = size(iNeighbour, dim=2)
    square(:, :) = 0.0_dp

    lShellVals(:) = 0
    rotZ(:,:) = 0.0_dp
    do iAtom1 = 1, nAtom
      ii = iAtomStart(iAtom1)
      nOrb1 = iAtomStart(iAtom1+1) - ii
      do iNeigh = 0, nNeighbourSK(iAtom1)
        iOrig = iSparseStart(iNeigh,iAtom1) + 1
        iAtom2 = iNeighbour(iNeigh, iAtom1)
        iAtom2f = img2CentCell(iAtom2)
        jj = iAtomStart(iAtom2f)
        nOrb2 = iAtomStart(iAtom2f+1) - jj
        tmpSqr(:nOrb2,:nOrb1) = reshape(orig(iOrig:iOrig+nOrb1*nOrb2-1), (/nOrb2,nOrb1/))
        iSp = species(iAtom2f)
        iSh = orb%nShell(iSp)
        lShellVals(:iSh) = orb%angShell(:iSh,iSp)
        theta = -atan2(coord(2,iAtom2),coord(1,iAtom2))&
            & + atan2(coord(2,iAtom2f),coord(1,iAtom2f))
        theta = mod(theta,2.0_dp*pi)
        call rotateZ(rotZ,lShellVals(:iSh), theta)
        tmpSqr(:nOrb2,:nOrb1) = matmul(rotZ(:nOrb2,:nOrb2),tmpSqr(:nOrb2,:nOrb1))
        square(jj:jj+nOrb2-1, ii:ii+nOrb1-1) = square(jj:jj+nOrb2-1, ii:ii+nOrb1-1)&
            & + tmpSqr(:nOrb2,:nOrb1)
      end do
    end do

  end subroutine unpackHSHelical_real


  !> Unpacks sparse matrix to square form (complex version)
  subroutine unpackDQ_cmplx_kpts(square, origBra, origKet, kPoint, iNeighbour, nNeighbourSK,&
      & iCellVec, cellVec, iAtomStart, iSparseStart, img2CentCell)

    !> Square form matrix on exit.
    complex(dp), intent(out) :: square(:, :, :)

    !> Sparse matrix
    real(dp), intent(in) :: origBra(:, :)

    !> Sparse matrix
    real(dp), intent(in) :: origKet(:, :)

    !> Relative coordinates of the K-point where the sparse matrix should be unfolded.
    real(dp), intent(in) :: kPoint(:)

    !> Neighbour list for each atom (First index from 0!)
    integer, intent(in) :: iNeighbour(0:, :)

    !> Nr. of neighbours for each atom (incl. itself).
    integer, intent(in) :: nNeighbourSK(:)

    !> Index of the cell translation vector for each atom.
    integer, intent(in) :: iCellVec(:)

    !> Relative coordinates of the cell translation vectors.
    real(dp), intent(in) :: cellVec(:, :)

    !> Atom offset for the square Hamiltonian
    integer, intent(in) :: iAtomStart(:)

    !> indexing array for the sparse Hamiltonian
    integer, intent(in) :: iSparseStart(0:, :)

    !> Map from images of atoms to central cell atoms
    integer, intent(in) :: img2CentCell(:)

    complex(dp) :: phase
    integer :: nAtom
    integer :: iOrig, ii, jj
    integer :: iNeigh
    integer :: iOldVec, iVec
    integer :: iAtom1, iAtom2, iAtom2f
    integer :: nOrb1, nOrb2, nmp
    real(dp) :: kPoint2p(3)

    nAtom = size(iNeighbour, dim=2)
    nmp = size(origBra, dim=1)

    @:ASSERT(nAtom > 0)
    @:ASSERT(size(origKet, dim=1) == size(origBra, dim=1))
    @:ASSERT(size(square, dim=1) == size(origBra, dim=1))
    @:ASSERT(size(square, dim=2) == size(square, dim=3))
    @:ASSERT(size(square, dim=2) == iAtomStart(nAtom+1) - 1)
    @:ASSERT(all(shape(kPoint) == [3]))
    @:ASSERT(all(shape(nNeighbourSK) == [nAtom]))
    @:ASSERT(size(iAtomStart) == nAtom + 1)

    square(:, :, :) = cmplx(0, 0, dp)
    kPoint2p(:) = 2.0_dp * pi * kPoint
    iOldVec = 0
    phase = 1.0_dp
    do iAtom1 = 1, nAtom
      ii = iAtomStart(iAtom1)
      nOrb1 = iAtomStart(iAtom1 + 1) - ii
      do iNeigh = 0, nNeighbourSK(iAtom1)
        iOrig = iSparseStart(iNeigh, iAtom1) + 1
        iAtom2 = iNeighbour(iNeigh, iAtom1)
        iAtom2f = img2CentCell(iAtom2)
        jj = iAtomStart(iAtom2f)
        @:ASSERT(jj >= ii)
        nOrb2 = iAtomStart(iAtom2f + 1) - jj
        iVec = iCellVec(iAtom2)
        if (iVec /= iOldVec) then
          phase = exp((0.0_dp, 1.0_dp) * dot_product(kPoint2p, cellVec(:, iVec)))
          iOldVec = iVec
        end if
        square(:, jj:jj+nOrb2-1, ii:ii+nOrb1-1) = square(:, jj:jj+nOrb2-1, ii:ii+nOrb1-1)&
            & + phase * reshape(origKet(:, iOrig:iOrig+nOrb1*nOrb2-1), [nmp, nOrb2, nOrb1])
        if (iAtom2f == iAtom1) cycle
        square(:, ii:ii+nOrb1-1, jj:jj+nOrb2-1) = square(:, ii:ii+nOrb1-1, jj:jj+nOrb2-1)&
            & + phase * reshape(origBra(:, iOrig:iOrig+nOrb1*nOrb2-1), [nmp, nOrb1, nOrb2], &
            &                   order=[1, 3, 2])
      end do
    end do

  end subroutine unpackDQ_cmplx_kpts


  !> Unpacks sparse matrix to square form (real version for Gamma point)
  subroutine unpackDQ_real(square, origBra, origKet, iNeighbour, nNeighbourSK, iAtomStart,&
      & iSparseStart, img2CentCell)

    !> Square form matrix on exit.
    real(dp), intent(out) :: square(:, :, :)

    !> Sparse matrix
    real(dp), intent(in) :: origBra(:, :)

    !> Sparse matrix
    real(dp), intent(in) :: origKet(:, :)

    !> Neighbour list for each atom (First index from 0!)
    integer, intent(in) :: iNeighbour(0:, :)

    !> Nr. of neighbours for each atom (incl. itself).
    integer, intent(in) :: nNeighbourSK(:)

    !> Atom offset for the square Hamiltonian
    integer, intent(in) :: iAtomStart(:)

    !> indexing array for the sparse Hamiltonian
    integer, intent(in) :: iSparseStart(0:, :)

    !> Map from images of atoms to central cell atoms
    integer, intent(in) :: img2CentCell(:)

    integer :: nAtom
    integer :: iOrig, ii, jj
    integer :: iNeigh
    integer :: iAtom1, iAtom2, iAtom2f
    integer :: nOrb1, nOrb2, nmp

    nAtom = size(iNeighbour, dim=2)
    nmp = size(origBra, dim=1)

    @:ASSERT(nAtom > 0)
    @:ASSERT(size(origKet, dim=1) == size(origBra, dim=1))
    @:ASSERT(size(square, dim=1) == size(origBra, dim=1))
    @:ASSERT(size(square, dim=2) == size(square, dim=3))
    @:ASSERT(size(square, dim=2) == iAtomStart(nAtom+1) - 1)
    @:ASSERT(all(shape(nNeighbourSK) == [nAtom]))
    @:ASSERT(size(iAtomStart) == nAtom + 1)

    square(:, :, :) = 0.0_dp

    do iAtom1 = 1, nAtom
      ii = iAtomStart(iAtom1)
      nOrb1 = iAtomStart(iAtom1 + 1) - ii
      do iNeigh = 0, nNeighbourSK(iAtom1)
        iOrig = iSparseStart(iNeigh, iAtom1) + 1
        iAtom2 = iNeighbour(iNeigh, iAtom1)
        iAtom2f = img2CentCell(iAtom2)
        jj = iAtomStart(iAtom2f)
        @:ASSERT(jj >= ii)
        nOrb2 = iAtomStart(iAtom2f + 1) - jj
        square(:, jj:jj+nOrb2-1, ii:ii+nOrb1-1) = square(:, jj:jj+nOrb2-1, ii:ii+nOrb1-1)&
            & + reshape(origKet(:, iOrig:iOrig+nOrb1*nOrb2-1), [nmp, nOrb2, nOrb1])
        if (iAtom2f == iAtom1) cycle
        square(:, ii:ii+nOrb1-1, jj:jj+nOrb2-1) = square(:, ii:ii+nOrb1-1, jj:jj+nOrb2-1)&
            & + reshape(origBra(:, iOrig:iOrig+nOrb1*nOrb2-1), [nmp, nOrb1, nOrb2], &
            &           order=[1, 3, 2])
      end do
    end do

  end subroutine unpackDQ_real


  !> Unpacks sparse matrices to square form (2 component version for k-points)
  !>
  !> Note: The non on-site blocks are only filled in the lower triangle part of the matrix.
  subroutine unpackHPauli(ham, kPoint, iNeighbour, nNeighbourSK, iSparseStart, iAtomStart,&
      & img2CentCell, iCellVec, cellVec, HSqrCplx, iHam)

    !> sparse hamiltonian
    real(dp), intent(in) :: ham(:, :)

    !> The k-point at which to unpack
    real(dp), intent(in) :: kPoint(:)

    !> Neighbour list for each atom (First index from 0!)
    integer, intent(in) :: iNeighbour(0:, :)

    !> Nr. of neighbours for each atom (incl. itself).
    integer, intent(in) :: nNeighbourSK(:)

    !> indexing array for the sparse Hamiltonian
    integer, intent(in) :: iSparseStart(:, :)

    !> Atom offset for the square Hamiltonian
    integer, intent(in) :: iAtomStart(:)

    !> Map from images of atoms to central cell atoms
    integer, intent(in) :: img2CentCell(:)

    !> index to vector to unit cell containing specified atom
    integer, intent(in) :: iCellVec(:)

    !> vectors to periodic unit cells
    real(dp), intent(in) :: cellVec(:, :)

    !> dense hamiltonian matrix
    complex(dp), intent(out) :: HSqrCplx(:, :)

    !> imaginary part of sparse hamiltonian
    real(dp), intent(in), optional :: iHam(:, :)

    complex(dp), allocatable :: work(:, :)
    integer :: nOrb
    integer :: ii

    nOrb = size(HSqrCplx, dim=1) / 2

    ! for the moment, but will use S as workspace in the future
    allocate(work(nOrb, nOrb))
    HSqrCplx(:, :) = 0.0_dp

    ! 1 0 charge part
    ! 0 1
    call unpackHS(work, ham(:, 1), kPoint, iNeighbour, nNeighbourSK, iCellVec,&
        & cellVec, iAtomStart, iSparseStart, img2CentCell)
    HSqrCplx(1:nOrb, 1:nOrb) = 0.5_dp*work(1:nOrb, 1:nOrb)
    HSqrCplx(nOrb+1:2*nOrb, nOrb+1:2*nOrb) = 0.5_dp*work(1:nOrb, 1:nOrb)
    if (present(iHam)) then
      call unpackHS(work, iHam(:, 1), kPoint, iNeighbour, nNeighbourSK, iCellVec,&
          & cellVec, iAtomStart, iSparseStart, img2CentCell)
      HSqrCplx(1:nOrb, 1:nOrb) = HSqrCplx(1:nOrb, 1:nOrb)&
          & + 0.5_dp*cmplx(0, 1, dp)*work(1:nOrb, 1:nOrb)
      HSqrCplx(nOrb+1:2*nOrb, nOrb+1:2*nOrb) =&
          & HSqrCplx(nOrb+1:2*nOrb, nOrb+1:2*nOrb)&
          & + 0.5_dp*cmplx(0, 1, dp)*work(1:nOrb, 1:nOrb)
    end if

    ! 0 1 x part
    ! 1 0
    call unpackHS(work, ham(:, 2), kPoint, iNeighbour, nNeighbourSK, iCellVec, cellVec, iAtomStart,&
        & iSparseStart, img2CentCell)
    do ii = 1, nOrb
      work(ii, ii+1:) = conjg(work(ii+1:, ii))
    end do

    HSqrCplx(nOrb+1:2*nOrb, 1:nOrb) = HSqrCplx(nOrb+1:2*nOrb, 1:nOrb)&
        & + 0.5_dp * work(1:nOrb, 1:nOrb)
    if (present(iHam)) then
      call unpackHS(work, iHam(:, 2), kPoint, iNeighbour, nNeighbourSK, iCellVec, cellVec,&
          & iAtomStart, iSparseStart, img2CentCell)
      do ii = 1, nOrb
        work(ii, ii+1:) = -conjg(work(ii+1:, ii))
      end do
      HSqrCplx(nOrb+1:2*nOrb, 1:nOrb) = HSqrCplx(nOrb+1:2*nOrb, 1:nOrb)&
          & + 0.5_dp * cmplx(0, 1, dp) * work(1:nOrb, 1:nOrb)
    end if

    ! 0 -i y part
    ! i  0
    call unpackHS(work, ham(:, 3), kPoint, iNeighbour, nNeighbourSK, iCellVec,&
        & cellVec, iAtomStart, iSparseStart, img2CentCell)
    do ii = 1, nOrb
      work(ii, ii+1:) = conjg(work(ii+1:, ii))
    end do

    HSqrCplx(nOrb+1:2*nOrb, 1:nOrb) = HSqrCplx(nOrb+1:2*nOrb, 1:nOrb)&
        & + cmplx(0.0, 0.5, dp) * work(1:nOrb, 1:nOrb)
    if (present(iHam)) then
      call unpackHS(work, iHam(:, 3), kPoint, iNeighbour, nNeighbourSK, iCellVec, cellVec,&
          & iAtomStart, iSparseStart, img2CentCell)

      ! Apply hermitian symmetry just in case
      do ii = 1, nOrb
        work(ii, ii+1:) = -conjg(work(ii+1:, ii))
      end do

      HSqrCplx(nOrb+1:2*nOrb, 1:nOrb) = HSqrCplx(nOrb+1:2*nOrb, 1:nOrb)&
          & - 0.5_dp * work(1:nOrb, 1:nOrb)
    end if

    ! 1  0 z part
    ! 0 -1
    call unpackHS(work, ham(:, 4), kPoint, iNeighbour, nNeighbourSK, iCellVec,&
        & cellVec, iAtomStart, iSparseStart, img2CentCell)
    HSqrCplx(1:nOrb, 1:nOrb) = HSqrCplx(1:nOrb, 1:nOrb)&
        & + 0.5_dp * work(1:nOrb, 1:nOrb)
    HSqrCplx(nOrb+1:2*nOrb, nOrb+1:2*nOrb) = HSqrCplx(nOrb+1:2*nOrb, nOrb+1:2*nOrb)&
        & - 0.5_dp * work(1:nOrb, 1:nOrb)
    if (present(iHam)) then
      call unpackHS(work, iHam(:, 4), kPoint, iNeighbour, nNeighbourSK, iCellVec,&
          & cellVec, iAtomStart, iSparseStart, img2CentCell)
      HSqrCplx(1:nOrb, 1:nOrb) = HSqrCplx(1:nOrb, 1:nOrb)&
          & + 0.5_dp * cmplx(0, 1, dp) * work(1:nOrb, 1:nOrb)
      HSqrCplx(nOrb+1:2*nOrb, nOrb+1:2*nOrb) = HSqrCplx(nOrb+1:2*nOrb, nOrb+1:2*nOrb)&
          & - 0.5_dp * cmplx(0, 1, dp) * work(1:nOrb, 1:nOrb)
    end if

  end subroutine unpackHPauli


  !> Unpacks sparse overlap matrices to square form (2 component version for k-points)
  !>
  !> Note: The non on-site blocks are only filled in the lower triangle part of the matrix.
  subroutine unpackSPauli(over, kPoint, iNeighbour, nNeighbourSK, iAtomStart, iSparseStart,&
      & img2CentCell, iCellVec, cellVec, SSqrCplx)

    !> sparse overlap matrix
    real(dp), intent(in) :: over(:)

    !> The k-point at which to unpack
    real(dp), intent(in) :: kPoint(:)

    !> Neighbour list for each atom (First index from 0!)
    integer, intent(in) :: iNeighbour(0:, :)

    !> Nr. of neighbours for each atom (incl. itself).
    integer, intent(in) :: nNeighbourSK(:)

    !> Atom offset for the square Hamiltonian
    integer, intent(in) :: iAtomStart(:)

    !> indexing array for the sparse Hamiltonian
    integer, intent(in) :: iSparseStart(:, :)

    !> Map from images of atoms to central cell atoms
    integer, intent(in) :: img2CentCell(:)

    !> index to vector to unit cell containing specified atom
    integer, intent(in) :: iCellVec(:)

    !> vectors to periodic unit cells
    real(dp), intent(in) :: cellVec(:, :)

    !> dense overlap matrix
    complex(dp), intent(out) :: SSqrCplx(:, :)

    complex(dp), allocatable :: work(:, :)
    integer :: nOrb

    nOrb = size(SSqrCplx, dim=1) / 2
    allocate(work(nOrb, nOrb))
    SSqrCplx(:, :) = 0.0_dp
    call unpackHS(work, over, kPoint, iNeighbour, nNeighbourSK, iCellVec, cellVec, iAtomStart,&
        & iSparseStart, img2CentCell)
    SSqrCplx(1:nOrb, 1:nOrb) = work(1:nOrb, 1:nOrb)
    SSqrCplx(nOrb + 1 : 2 * nOrb, nOrb + 1 : 2 * nOrb) = work(1:nOrb, 1:nOrb)

  end subroutine unpackSPauli


  !> Pack squared matrix in the sparse form (complex version).
  subroutine packHS_cmplx_kpts(primitive, square, kPoint, kWeight, iNeighbour, nNeighbourSK, mOrb,&
      & iCellVec, cellVec, iAtomStart, iSparseStart, img2CentCell)

    !> Sparse matrix
    real(dp), intent(inout) :: primitive(:)

    !> Square form matrix
    complex(dp), intent(in) :: square(:,:)

    !> Relative coordinates of the K-point
    real(dp), intent(in) :: kPoint(:)

    !> Weight of the K-point
    real(dp), intent(in) :: kweight

    !> Neighbour list for the atoms (First index from 0!)
    integer, intent(in) :: iNeighbour(0:,:)

    !> Nr. of neighbours for the atoms.
    integer, intent(in) :: nNeighbourSK(:)

    !> Maximal number of orbitals on an atom.
    integer, intent(in) :: mOrb

    !> Index of the cell translation vector for each atom.
    integer, intent(in) :: iCellVec(:)

    !> Relative coordinates of the cell translation vectors.
    real(dp), intent(in) :: cellVec(:,:)

    !> Atom offset for the square matrix
    integer, intent(in) :: iAtomStart(:)

    !> indexing array for the sparse Hamiltonian
    integer, intent(in) :: iSparseStart(0:,:)

    !> Mapping between image atoms and corresponding atom in the central cell.
    integer, intent(in) :: img2CentCell(:)

    complex(dp) :: phase
    integer :: nAtom
    integer :: iOrig, ii, jj, kk
    integer :: iNeigh
    integer :: iOldVec, iVec
    integer :: iAtom1, iAtom2, iAtom2f
    integer :: nOrb1, nOrb2
    real(dp) :: kPoint2p(3)
    complex(dp) :: tmpSqr(mOrb, mOrb)
  #:block DEBUG_CODE
    integer :: sizePrim
  #:endblock DEBUG_CODE

    nAtom = size(iNeighbour, dim=2)
  #:block DEBUG_CODE
    sizePrim = size(primitive)
  #:endblock DEBUG_CODE

    @:ASSERT(nAtom > 0)
    @:ASSERT(size(square, dim=1) == size(square, dim=2))
    @:ASSERT(size(square, dim=1) == iAtomStart(nAtom+1) - 1)
    @:ASSERT(all(shape(kPoint) == [3]))
    @:ASSERT(all(shape(nNeighbourSK) == [nAtom]))
    @:ASSERT(kWeight > 0.0_dp)
    @:ASSERT(size(iAtomStart) == nAtom + 1)

    kPoint2p(:) = 2.0_dp * pi * kPoint
    iOldVec = 0
    phase = 1.0_dp
    do iAtom1 = 1, nAtom
      ii = iAtomStart(iAtom1)
      nOrb1 = iAtomStart(iAtom1 + 1) - ii
      do iNeigh = 0, nNeighbourSK(iAtom1)
        iOrig = iSparseStart(iNeigh, iAtom1) + 1
        iAtom2 = iNeighbour(iNeigh, iAtom1)
        iAtom2f = img2CentCell(iAtom2)
        jj = iAtomStart(iAtom2f)
        @:ASSERT(jj >= ii)
        nOrb2 = iAtomStart(iAtom2f + 1) - jj
        iVec = iCellVec(iAtom2)
        if (iVec /= iOldVec) then
          phase = exp(cmplx(0, -1, dp) * dot_product(kPoint2p(:), cellVec(:, iVec)))
          iOldVec = iVec
        end if
        tmpSqr(1:nOrb2, 1:nOrb1) = square(jj:jj+nOrb2-1, ii:ii+nOrb1-1)

        ! Hermitian the on-site block before packing, just in case
        if (iAtom1 == iAtom2f) then
          do kk = 1, nOrb2
            tmpSqr(kk, kk+1:nOrb1) = conjg(tmpSqr(kk+1:nOrb1, kk))
          end do
        end if

        @:ASSERT(sizePrim >= iOrig + nOrb1*nOrb2 - 1)
        primitive(iOrig : iOrig + nOrb1 * nOrb2 - 1) = primitive(iOrig : iOrig + nOrb1*nOrb2 - 1)&
            & + kWeight * real(phase * reshape(tmpSqr(1:nOrb2, 1:nOrb1), [nOrb1*nOrb2]), dp)
      end do
    end do

  end subroutine packHS_cmplx_kpts


  !> Pack square matrix in the sparse form (real version).
  subroutine packHS_real(primitive, square, iNeighbour, nNeighbourSK, mOrb, iAtomStart,&
      & iSparseStart, img2CentCell)

    !> Sparse matrix
    real(dp), intent(inout) :: primitive(:)

    !> Square form matrix
    real(dp), intent(in) :: square(:, :)

    !> Neighbour list for the atoms (First index from 0!)
    integer, intent(in) :: iNeighbour(0:, :)

    !> Nr. of neighbours for the atoms.
    integer, intent(in) :: nNeighbourSK(:)

    !> Maximal number of orbitals on an atom.
    integer, intent(in) :: mOrb

    !> Atom offset for the square matrix
    integer, intent(in) :: iAtomStart(:)

    !> indexing array for the sparse Hamiltonian
    integer, intent(in) :: iSparseStart(0:, :)

    !> Mapping between image atoms and corresponding atom in the central cell.
    integer, intent(in) :: img2CentCell(:)

    integer :: nAtom
    integer :: iOrig, ii, jj, kk
    integer :: iNeigh
    integer :: iAtom1, iAtom2, iAtom2f
    integer :: nOrb1, nOrb2
    real(dp) :: tmpSqr(mOrb, mOrb)
  #:block DEBUG_CODE
    integer :: sizePrim
  #:endblock DEBUG_CODE

    nAtom = size(iNeighbour, dim=2)
  #:block DEBUG_CODE
    sizePrim = size(primitive)
  #:endblock DEBUG_CODE

    @:ASSERT(nAtom > 0)
    @:ASSERT(size(square, dim=1) == size(square, dim=2))
    @:ASSERT(size(square, dim=1) == iAtomStart(nAtom+1) - 1)
    @:ASSERT(all(shape(nNeighbourSK) == [nAtom]))

    do iAtom1 = 1, nAtom
      ii = iAtomStart(iAtom1)
      nOrb1 = iAtomStart(iAtom1 + 1) - ii
      do iNeigh = 0, nNeighbourSK(iAtom1)
        iOrig = iSparseStart(iNeigh, iAtom1) + 1
        iAtom2 = iNeighbour(iNeigh, iAtom1)
        iAtom2f = img2CentCell(iAtom2)
        jj = iAtomStart(iAtom2f)
        @:ASSERT(jj >= ii)
        nOrb2 = iAtomStart(iAtom2f + 1) - jj
        tmpSqr(1:nOrb2, 1:nOrb1) = square(jj:jj+nOrb2-1, ii:ii+nOrb1-1)

        ! Symmetrize the on-site block before packing, just in case
        if (iAtom1 == iAtom2f) then
          do kk = 1, nOrb2
            tmpSqr(kk, kk+1:nOrb1) = tmpSqr(kk+1:nOrb1, kk)
          end do
        end if

        @:ASSERT(sizePrim >= iOrig + nOrb1*nOrb2 - 1)
        primitive(iOrig : iOrig + nOrb1*nOrb2 - 1) = primitive(iOrig : iOrig + nOrb1 * nOrb2 - 1)&
            & + reshape(tmpSqr(1:nOrb2, 1:nOrb1), [nOrb1 * nOrb2])
      end do
    end do

  end subroutine packHS_real


  !> Pack square matrix in the sparse form (complex version) for helical boundary conditions.
  subroutine packHShelical_cmplx(primitive, square, kPoint, kWeight, iNeighbour, nNeighbourSK,&
      & mOrb, iCellVec, cellVec, iAtomStart, iSparseStart, img2CentCell, orb, species, coord)

    !> Sparse matrix
    real(dp), intent(inout) :: primitive(:)

    !> Square form matrix
    complex(dp), intent(in) :: square(:, :)

    !> Relative coordinates of the K-point
    real(dp), intent(in) :: kPoint(:)

    !> Weight of the K-point
    real(dp), intent(in) :: kweight

    !> Neighbour list for the atoms (First index from 0!)
    integer, intent(in) :: iNeighbour(0:, :)

    !> Nr. of neighbours for the atoms.
    integer, intent(in) :: nNeighbourSK(:)

    !> Maximal number of orbitals on an atom.
    integer, intent(in) :: mOrb

    !> Index of the cell translation vector for each atom.
    integer, intent(in) :: iCellVec(:)

    !> Relative coordinates of the cell translation vectors.
    real(dp), intent(in) :: cellVec(:, :)

    !> Atom offset for the square matrix
    integer, intent(in) :: iAtomStart(:)

    !> indexing array for the sparse Hamiltonian
    integer, intent(in) :: iSparseStart(0:, :)

    !> Mapping between image atoms and corresponding atom in the central cell.
    integer, intent(in) :: img2CentCell(:)

    !> data type for atomic orbital information
    type(TOrbitals), intent(in) :: orb

    !> Species of each atom
    integer :: species(:)

    !> Coordinates of all atoms
    real(dp), intent(in) :: coord(:,:)

    complex(dp) :: phase
    integer :: nAtom, iOrig, ii, jj, kk, iNeigh, iOldVec, iVec, iAtom1, iAtom2, iAtom2f
    integer :: nOrb1, nOrb2, lShellVals(orb%mShell), iSh, iSp
    real(dp) :: kPoint2p(2), tmpSqrR(mOrb, mOrb), rotZ(orb%mOrb,orb%mOrb), theta
    complex(dp) :: tmpSqr(mOrb, mOrb)

    nAtom = size(iNeighbour, dim=2)
    kPoint2p(:) = 2.0_dp * pi * kPoint
    iOldVec = 0
    phase = 1.0_dp
    lShellVals(:) = 0
    do iAtom1 = 1, nAtom
      ii = iAtomStart(iAtom1)
      nOrb1 = iAtomStart(iAtom1 + 1) - ii
      do iNeigh = 0, nNeighbourSK(iAtom1)
        iOrig = iSparseStart(iNeigh, iAtom1) + 1
        iAtom2 = iNeighbour(iNeigh, iAtom1)
        iAtom2f = img2CentCell(iAtom2)
        jj = iAtomStart(iAtom2f)
        nOrb2 = iAtomStart(iAtom2f + 1) - jj
        iVec = iCellVec(iAtom2)
        iSp = species(iAtom2f)
        iSh = orb%nShell(iSp)
        lShellVals(:iSh) = orb%angShell(:iSh,iSp)
        theta = atan2(coord(2,iAtom2),coord(1,iAtom2)) &
            & - atan2(coord(2,iAtom2f),coord(1,iAtom2f))
        theta = mod(theta,2.0_dp*pi)
        call rotateZ(rotZ,lShellVals(:iSh), theta)
        if (iVec /= iOldVec) then
          phase = exp(cmplx(0, -1, dp) * dot_product(kPoint2p(:), cellVec(:, iVec)))
          iOldVec = iVec
        end if
        tmpSqr(1:nOrb2, 1:nOrb1) = square(jj:jj+nOrb2-1, ii:ii+nOrb1-1)
        ! Hermitian the on-site block before packing, just in case
        if (iAtom1 == iAtom2f) then
          do kk = 1, nOrb2
            tmpSqr(kk, kk+1:nOrb1) = conjg(tmpSqr(kk+1:nOrb1, kk))
          end do
        end if
        ! rotate
        tmpSqrR(:nOrb2,:nOrb1) = real(phase*tmpSqr(:nOrb2,:nOrb1))
        tmpSqrR(:nOrb2,:nOrb1) = matmul(rotZ(:nOrb2,:nOrb2),tmpSqrR(:nOrb2,:nOrb1))
        primitive(iOrig : iOrig + nOrb1 * nOrb2 - 1) = primitive(iOrig : iOrig + nOrb1*nOrb2 - 1)&
            & + kWeight * reshape(tmpSqrR(1:nOrb2, 1:nOrb1), [nOrb1*nOrb2])
      end do
    end do

  end subroutine packHShelical_cmplx


  !> Pack square matrix in the sparse form (real version).
  subroutine packHShelical_real(primitive, square, iNeighbour, nNeighbourSK, iAtomStart,&
      & iSparseStart, img2CentCell, orb, species, coord)

    !> Sparse matrix
    real(dp), intent(inout) :: primitive(:)

    !> Square form matrix
    real(dp), intent(in) :: square(:, :)

    !> Neighbour list for the atoms (First index from 0!)
    integer, intent(in) :: iNeighbour(0:, :)

    !> Nr. of neighbours for the atoms.
    integer, intent(in) :: nNeighbourSK(:)

    !> Atom offset for the square matrix
    integer, intent(in) :: iAtomStart(:)

    !> indexing array for the sparse Hamiltonian
    integer, intent(in) :: iSparseStart(0:, :)

    !> Mapping between image atoms and corresponding atom in the central cell.
    integer, intent(in) :: img2CentCell(:)

    !> data type for atomic orbital information
    type(TOrbitals), intent(in) :: orb

    !> Species of each atom
    integer :: species(:)

    !> Coordinates of all atoms
    real(dp), intent(in) :: coord(:,:)

    integer :: nAtom, iOrig, ii, jj, kk, iNeigh, iAtom1, iAtom2, iAtom2f
    integer :: nOrb1, nOrb2, iSp, iSh, lShellVals(orb%mShell)
    real(dp) :: tmpSqr(orb%mOrb, orb%mOrb), rotZ(orb%mOrb,orb%mOrb), theta

    nAtom = size(iNeighbour, dim=2)
    lShellVals(:) = 0
    do iAtom1 = 1, nAtom
      ii = iAtomStart(iAtom1)
      nOrb1 = iAtomStart(iAtom1 + 1) - ii
      do iNeigh = 0, nNeighbourSK(iAtom1)
        iOrig = iSparseStart(iNeigh, iAtom1) + 1
        iAtom2 = iNeighbour(iNeigh, iAtom1)
        iAtom2f = img2CentCell(iAtom2)
        jj = iAtomStart(iAtom2f)
        nOrb2 = iAtomStart(iAtom2f + 1) - jj
        tmpSqr(1:nOrb2, 1:nOrb1) = square(jj:jj+nOrb2-1, ii:ii+nOrb1-1)
        ! Symmetrize the on-site block before packing, just in case
        if (iAtom1 == iAtom2f) then
          do kk = 1, nOrb2
            tmpSqr(kk, kk+1:nOrb1) = tmpSqr(kk+1:nOrb1, kk)
          end do
        end if
        iSp = species(iAtom2f)
        iSh = orb%nShell(iSp)
        lShellVals(:iSh) = orb%angShell(:iSh,iSp)
        theta = atan2(coord(2,iAtom2),coord(1,iAtom2)) &
            & - atan2(coord(2,iAtom2f),coord(1,iAtom2f))
        theta = mod(theta,2.0_dp*pi)
        call rotateZ(rotZ, lShellVals(:iSh), theta)
        tmpSqr(:nOrb2,:nOrb1) =  matmul(rotZ(:nOrb2,:nOrb2),tmpSqr(:nOrb2,:nOrb1))

        primitive(iOrig : iOrig + nOrb1*nOrb2 - 1) = primitive(iOrig : iOrig + nOrb1 * nOrb2 - 1)&
            & + reshape(tmpSqr(1:nOrb2, 1:nOrb1), [nOrb1 * nOrb2])
      end do
    end do

  end subroutine packHShelical_real


  !> Pack square matrix in the sparse form (real Pauli version).
  subroutine packHSPauli(primitive, square, iNeighbour, nNeighbourSK, mOrb, iAtomStart,&
      & iSparseStart, img2CentCell)

    !> Sparse matrix
    real(dp), intent(inout) :: primitive(:, :)

    !> Square form matrix
    complex(dp), intent(in) :: square(:, :)

    !> Neighbour list for the atoms (First index from 0!)
    integer, intent(in) :: iNeighbour(0:, :)

    !> Nr. of neighbours for the atoms.
    integer, intent(in) :: nNeighbourSK(:)

    !> Maximal number of orbitals on an atom.
    integer, intent(in) :: mOrb

    !> Atom offset for the square matrix
    integer, intent(in) :: iAtomStart(:)

    !> indexing array for the sparse Hamiltonian
    integer, intent(in) :: iSparseStart(0:, :)

    !> Mapping between image atoms and corresponding atom in the central cell.
    integer, intent(in) :: img2CentCell(:)

    integer :: nAtom
    integer :: iOrig, ii, jj, kk
    integer :: iNeigh, iBlock
    integer :: iAtom1, iAtom2, iAtom2f
    integer :: nOrb1, nOrb2, nOrb
    complex(dp) :: tmpSqr(mOrb, mOrb)
  #:block DEBUG_CODE
    integer :: sizePrim
  #:endblock DEBUG_CODE

    nAtom = size(iNeighbour, dim=2)
    ! number of orbitals in a regular spin block
    nOrb = (iAtomStart(nAtom+1) - 1)

  #:block DEBUG_CODE
    sizePrim = size(primitive, dim=1)
  #:endblock DEBUG_CODE

    @:ASSERT(nAtom > 0)
    @:ASSERT(size(square, dim=1) == size(square, dim=2))
    @:ASSERT(size(square, dim=1) == 2 * nOrb )
    @:ASSERT(all(shape(nNeighbourSK) == [nAtom]))
    @:ASSERT(size(iAtomStart) == nAtom + 1)
    @:ASSERT(size(primitive, dim=2)==4)

    do iBlock = 0, 1
      do iAtom1 = 1, nAtom
        ii = iAtomStart(iAtom1)
        nOrb1 = iAtomStart(iAtom1 + 1) - ii
        do iNeigh = 0, nNeighbourSK(iAtom1)
          iOrig = iSparseStart(iNeigh, iAtom1) + 1
          iAtom2 = iNeighbour(iNeigh, iAtom1)
          iAtom2f = img2CentCell(iAtom2)
          jj = iAtomStart(iAtom2f)
          @:ASSERT(jj >= ii)
          nOrb2 = iAtomStart(iAtom2f + 1) - jj
          tmpSqr(1:nOrb2, 1:nOrb1) =&
              & square(jj+iBlock*nOrb:jj+nOrb2-1+iBlock*nOrb, ii+iBlock*nOrb:ii+nOrb1-1+iBlock*nOrb)
          ! Hermitian the on-site block before packing, as only one triangle usually supplied
          if (iAtom1 == iAtom2f) then
            do kk = 1, nOrb2
              tmpSqr(kk, kk+1:nOrb1) = conjg(tmpSqr(kk+1:nOrb1, kk))
            end do
          end if
          @:ASSERT(sizePrim >= iOrig + nOrb1*nOrb2 - 1)
          primitive(iOrig : iOrig + nOrb1*nOrb2 - 1, 1) =&
              & primitive(iOrig : iOrig + nOrb1*nOrb2 - 1, 1)&
              & + 0.5_dp*reshape(real(tmpSqr(1:nOrb2, 1:nOrb1)), [nOrb1*nOrb2])
          primitive(iOrig : iOrig + nOrb1*nOrb2 - 1, 4) =&
              &primitive(iOrig : iOrig + nOrb1*nOrb2 - 1, 4)&
              & + real(1-2*iBlock, dp) *&
              & 0.5_dp*reshape(real(tmpSqr(1:nOrb2, 1:nOrb1)), [nOrb1*nOrb2])
        end do
      end do
    end do

    do iAtom1 = 1, nAtom
      ii = iAtomStart(iAtom1)
      nOrb1 = iAtomStart(iAtom1 + 1) - ii
      do iNeigh = 0, nNeighbourSK(iAtom1)
        iOrig = iSparseStart(iNeigh, iAtom1) + 1
        iAtom2 = iNeighbour(iNeigh, iAtom1)
        iAtom2f = img2CentCell(iAtom2)
        jj = iAtomStart(iAtom2f)
        @:ASSERT(jj >= ii)
        nOrb2 = iAtomStart(iAtom2f + 1) - jj
        ! take the Hermitian part of the block
        tmpSqr(1:nOrb2, 1:nOrb1) = 0.5_dp * (square(jj+nOrb:jj+nOrb2-1+nOrb, ii:ii+nOrb1-1)&
            & + transpose(square(nOrb+ii:nOrb+ii+nOrb1-1, jj:jj+nOrb2-1)))
        @:ASSERT(sizePrim >= iOrig + nOrb1*nOrb2 - 1)
        primitive(iOrig : iOrig + nOrb1*nOrb2 - 1, 2) =&
            & primitive(iOrig : iOrig + nOrb1*nOrb2 - 1, 2)&
            & + reshape(real(tmpSqr(1:nOrb2, 1:nOrb1)), [nOrb1*nOrb2])
        primitive(iOrig : iOrig + nOrb1*nOrb2 - 1, 3) =&
            &primitive(iOrig : iOrig + nOrb1*nOrb2 - 1, 3)&
            & +reshape(aimag(tmpSqr(1:nOrb2, 1:nOrb1)), [nOrb1*nOrb2])
      end do
    end do

  end subroutine packHSPauli


  !> Pack square matrix into the sparse form (complex Pauli version).
  subroutine packHSPauli_kpts(primitive, square, kPoint, kWeight, iNeighbour, nNeighbourSK, mOrb,&
      & iCellVec, cellVec, iAtomStart, iSparseStart, img2CentCell)

    !> Sparse matrix
    real(dp), intent(inout) :: primitive(:, :)

    !> Square form matrix
    complex(dp), intent(in) :: square(:, :)

    !> location in the BZ in units of 2pi
    real(dp), intent(in) :: kPoint(:)

    !> Weight of the k-point
    real(dp), intent(in) :: kweight

    !> Neighbour list for the atoms (First index from 0!)
    integer, intent(in) :: iNeighbour(0:, :)

    !> Nr. of neighbours for the atoms.
    integer, intent(in) :: nNeighbourSK(:)

    !> Maximal number of orbitals on an atom.
    integer, intent(in) :: mOrb

    !> Index of the cell translation vector for each atom.
    integer, intent(in) :: iCellVec(:)

    !> Relative coordinates of the cell translation vectors.
    real(dp), intent(in) :: cellVec(:, :)

    !> Atom offset for the square matrix
    integer, intent(in) :: iAtomStart(:)

    !> indexing array for the sparse Hamiltonian
    integer, intent(in) :: iSparseStart(0:, :)

    !> Mapping between image atoms and corresponding atom in the central cell.
    integer, intent(in) :: img2CentCell(:)

    complex(dp) :: phase
    integer :: nAtom
    integer :: iOrig, ii, jj, kk
    integer :: iNeigh, iBlock
    integer :: iOldVec, iVec
    integer :: iAtom1, iAtom2, iAtom2f
    integer :: nOrb1, nOrb2, nOrb
    real(dp) :: kPoint2p(3)
    complex(dp) :: tmpSqr(mOrb, mOrb)
  #:block DEBUG_CODE
    integer :: sizePrim
  #:endblock DEBUG_CODE

    nAtom = size(iNeighbour, dim=2)
    ! number of orbitals in a regular spin block
    nOrb = (iAtomStart(nAtom+1) - 1)

  #:block DEBUG_CODE
    sizePrim = size(primitive, dim=1)
  #:endblock DEBUG_CODE

    @:ASSERT(nAtom > 0)
    @:ASSERT(size(square, dim=1) == size(square, dim=2))
    @:ASSERT(size(square, dim=1) == 2 * nOrb )
    @:ASSERT(all(shape(kPoint) == [3]))
    @:ASSERT(all(shape(nNeighbourSK) == [nAtom]))
    @:ASSERT(size(iAtomStart) == nAtom + 1)
    @:ASSERT(kWeight > 0.0_dp)
    @:ASSERT(size(primitive, dim=2)==4)

    kPoint2p(:) = 2.0_dp * pi * kPoint

    ! sigma_I and sigma_z blocks
    do iBlock = 0, 1
      iOldVec = 0
      phase = 1.0_dp
      do iAtom1 = 1, nAtom
        ii = iAtomStart(iAtom1)
        nOrb1 = iAtomStart(iAtom1 + 1) - ii
        do iNeigh = 0, nNeighbourSK(iAtom1)
          iOrig = iSparseStart(iNeigh, iAtom1) + 1
          iAtom2 = iNeighbour(iNeigh, iAtom1)
          iAtom2f = img2CentCell(iAtom2)
          jj = iAtomStart(iAtom2f)
          @:ASSERT(jj >= ii)
          nOrb2 = iAtomStart(iAtom2f + 1) - jj
          iVec = iCellVec(iAtom2)
          if (iVec /= iOldVec) then
            phase = exp(cmplx(0, -1, dp) * dot_product(kPoint2p(:), cellVec(:, iVec)))
            iOldVec = iVec
          end if
          tmpSqr(1:nOrb2, 1:nOrb1) = square(jj+iBlock*nOrb:jj+nOrb2-1+iBlock*nOrb,&
              & ii+iBlock*nOrb:ii+nOrb1-1+iBlock*nOrb)
          ! Hermitian the on-site block before packing, as only one triangle usually supplied
          if (iAtom1 == iAtom2f) then
            do kk = 1, nOrb2
              tmpSqr(kk, kk+1:nOrb1) = conjg(tmpSqr(kk+1:nOrb1, kk))
            end do
          end if
          @:ASSERT(sizePrim >= iOrig + nOrb1*nOrb2 - 1)
          primitive(iOrig : iOrig + nOrb1*nOrb2 - 1, 1) =&
              & primitive(iOrig : iOrig + nOrb1*nOrb2 - 1, 1)&
              & + kWeight * 0.5_dp*reshape(real(phase* tmpSqr(1:nOrb2, 1:nOrb1)), [nOrb1*nOrb2])
          primitive(iOrig : iOrig + nOrb1*nOrb2 - 1, 4) =&
              & primitive(iOrig : iOrig + nOrb1*nOrb2 - 1, 4)&
              & + real(1-2*iBlock, dp) *&
              & kWeight * 0.5_dp*reshape(real(phase*tmpSqr(1:nOrb2, 1:nOrb1)), [nOrb1*nOrb2])
        end do
      end do
    end do

    ! sigma_x and sigma_y blocks
    iOldVec = 0
    phase = 1.0_dp
    do iAtom1 = 1, nAtom
      ii = iAtomStart(iAtom1)
      nOrb1 = iAtomStart(iAtom1 + 1) - ii
      do iNeigh = 0, nNeighbourSK(iAtom1)
        iOrig = iSparseStart(iNeigh, iAtom1) + 1
        iAtom2 = iNeighbour(iNeigh, iAtom1)
        iAtom2f = img2CentCell(iAtom2)
        jj = iAtomStart(iAtom2f)
        @:ASSERT(jj >= ii)
        nOrb2 = iAtomStart(iAtom2f + 1) - jj
        iVec = iCellVec(iAtom2)
        if (iVec /= iOldVec) then
          phase = exp(cmplx(0, -1, dp) * dot_product(kPoint2p(:), cellVec(:, iVec)))
          iOldVec = iVec
        end if
        @:ASSERT(sizePrim >= iOrig + nOrb1*nOrb2 - 1)

        ! take the Pauli part of the block
        tmpSqr(1:nOrb2, 1:nOrb1) = 0.5_dp * (square(jj+nOrb:jj+nOrb2-1+nOrb, ii:ii+nOrb1-1)&
            & + conjg(transpose(square(nOrb+ii:nOrb+ii+nOrb1-1, jj:jj+nOrb2-1))))
        if (iAtom1 == iAtom2f) then
          ! make up the other side of the on-site block
          do kk = 1, nOrb2
            tmpSqr(kk, kk+1:nOrb1) = tmpSqr(kk+1:nOrb1, kk)
          end do
        end if
        primitive(iOrig : iOrig + nOrb1 * nOrb2 - 1, 2) =&
            & primitive(iOrig : iOrig + nOrb1 * nOrb2 - 1, 2)&
            & + kWeight * reshape(real(phase * tmpSqr(1:nOrb2, 1:nOrb1)), [nOrb1 * nOrb2])
        tmpSqr(1:nOrb2, 1:nOrb1) =&
            & 0.5_dp * (square(jj + nOrb : jj + nOrb2 - 1 + nOrb, ii : ii + nOrb1 - 1)&
            & - conjg(transpose(square(nOrb +ii : nOrb + ii + nOrb1 - 1, jj : jj + nOrb2 - 1))))
        if (iAtom1 == iAtom2f) then !make up the other side of the on-site block
          do kk = 1, nOrb2
            tmpSqr(kk, kk+1:nOrb1) = tmpSqr(kk+1:nOrb1, kk)
          end do
        end if
        primitive(iOrig : iOrig + nOrb1 * nOrb2 - 1, 3) =&
            & primitive(iOrig : iOrig + nOrb1 * nOrb2 - 1, 3)&
            & + kWeight*reshape(aimag(phase * tmpSqr(1:nOrb2, 1:nOrb1)), [nOrb1 * nOrb2])
      end do
    end do

  end subroutine packHSPauli_kpts


  !> Pack imaginary coefficient part of Pauli square matrix into the sparse form.
  subroutine packHSPauliImag(primitive, square, iNeighbour, nNeighbourSK, mOrb, iAtomStart,&
      & iSparseStart, img2CentCell)

    !> Sparse matrix
    real(dp), intent(inout) :: primitive(:, :)

    !> Square form matrix
    complex(dp), intent(in) :: square(:, :)

    !> Neighbour list for the atoms (First index from 0!)
    integer, intent(in) :: iNeighbour(0:, :)

    !> Nr. of neighbours for the atoms.
    integer, intent(in) :: nNeighbourSK(:)

    !> Maximal number of orbitals on an atom.
    integer, intent(in) :: mOrb

    !> Atom offset for the square matrix
    integer, intent(in) :: iAtomStart(:)

    !> indexing array for the sparse Hamiltonian
    integer, intent(in) :: iSparseStart(0:, :)

    !> Mapping between image atoms and corresponding atom in the central cell.
    integer, intent(in) :: img2CentCell(:)

    integer :: nAtom
    integer :: iOrig, ii, jj, kk
    integer :: iNeigh, iBlock
    integer :: iAtom1, iAtom2, iAtom2f
    integer :: nOrb1, nOrb2, nOrb
    complex(dp) :: tmpSqr(mOrb, mOrb)
  #:block DEBUG_CODE
    integer :: sizePrim
  #:endblock DEBUG_CODE

    nAtom = size(iNeighbour, dim=2)
    ! number of orbitals in a regular spin block
    nOrb = (iAtomStart(nAtom+1) - 1)


  #:block DEBUG_CODE
    sizePrim = size(primitive, dim=1)
  #:endblock DEBUG_CODE

    @:ASSERT(nAtom > 0)
    @:ASSERT(size(square, dim=1) == size(square, dim=2))
    @:ASSERT(size(square, dim=1) == 2 * nOrb )
    @:ASSERT(all(shape(nNeighbourSK) == [nAtom]))
    @:ASSERT(size(iAtomStart) == nAtom + 1)
    @:ASSERT(size(primitive, dim=2)==4)

    do iBlock = 0, 1
      do iAtom1 = 1, nAtom
        ii = iAtomStart(iAtom1)
        nOrb1 = iAtomStart(iAtom1 + 1) - ii
        do iNeigh = 0, nNeighbourSK(iAtom1)
          iOrig = iSparseStart(iNeigh, iAtom1) + 1
          iAtom2 = iNeighbour(iNeigh, iAtom1)
          iAtom2f = img2CentCell(iAtom2)
          jj = iAtomStart(iAtom2f)
          @:ASSERT(jj >= ii)
          nOrb2 = iAtomStart(iAtom2f + 1) - jj
          tmpSqr(1:nOrb2, 1:nOrb1) = square(jj + iBlock *nOrb :jj + nOrb2 - 1 + iBlock * nOrb,&
              & ii + iBlock * nOrb : ii + nOrb1 - 1 + iBlock * nOrb)
          ! Hermitian the on-site block before packing, as only one triangle usually supplied
          if (iAtom1 == iAtom2f) then
            do kk = 1, nOrb2
              tmpSqr(kk, kk+1:nOrb1) = conjg(tmpSqr(kk+1:nOrb1, kk))
            end do
          end if
          @:ASSERT(sizePrim >= iOrig + nOrb1*nOrb2 - 1)
          primitive(iOrig : iOrig + nOrb1*nOrb2 - 1, 1) =&
              & primitive(iOrig : iOrig + nOrb1 * nOrb2 - 1, 1)&
              & + 0.5_dp * reshape(aimag(tmpSqr(1:nOrb2, 1:nOrb1)), [nOrb1 * nOrb2])
          primitive(iOrig : iOrig + nOrb1*nOrb2 - 1, 4) =&
              & primitive(iOrig : iOrig + nOrb1*nOrb2 - 1, 4)&
              & + real(1-2*iBlock, dp)&
              & * 0.5_dp * reshape(aimag(tmpSqr(1:nOrb2, 1:nOrb1)), [nOrb1 * nOrb2])
        end do
      end do
    end do

    do iAtom1 = 1, nAtom
      ii = iAtomStart(iAtom1)
      nOrb1 = iAtomStart(iAtom1 + 1) - ii
      do iNeigh = 0, nNeighbourSK(iAtom1)
        iOrig = iSparseStart(iNeigh, iAtom1) + 1
        iAtom2 = iNeighbour(iNeigh, iAtom1)
        iAtom2f = img2CentCell(iAtom2)
        jj = iAtomStart(iAtom2f)
        @:ASSERT(jj >= ii)
        nOrb2 = iAtomStart(iAtom2f + 1) - jj
        ! take the anti-Hermitian part of the block
        tmpSqr(1:nOrb2, 1:nOrb1) = 0.5_dp * ( square(nOrb+jj:nOrb+jj+nOrb2-1, ii:ii+nOrb1-1)&
            & - transpose(square(nOrb+ii:nOrb+ii+nOrb1-1, jj:jj+nOrb2-1)))
        @:ASSERT(sizePrim >= iOrig + nOrb1*nOrb2 - 1)
        ! sigma_x : i * 1 = i use + imaginary part of block
        primitive(iOrig : iOrig + nOrb1*nOrb2 - 1, 2) =&
            & primitive(iOrig : iOrig + nOrb1 * nOrb2 - 1, 2)&
            & + reshape(aimag(tmpSqr(1:nOrb2, 1:nOrb1)), [nOrb1 * nOrb2])
        ! sigma_y : i * i = -1 use - real part of block
        primitive(iOrig : iOrig + nOrb1*nOrb2 - 1, 3) =&
            & primitive(iOrig : iOrig + nOrb1*nOrb2 - 1, 3)&
            & - reshape(real(tmpSqr(1:nOrb2, 1:nOrb1)), [nOrb1 * nOrb2])
      end do
    end do

  end subroutine packHSPauliImag


  !> Pack imaginary coefficient part of Pauli square matrix into the sparse form (complex version).
  subroutine packHSPauliImag_kpts(primitive, square, kPoint, kWeight, iNeighbour, nNeighbourSK,&
      & mOrb, iCellVec, cellVec, iAtomStart, iSparseStart, img2CentCell)

    !> Sparse matrix
    real(dp), intent(inout) :: primitive(:, :)

    !> Square form matrix
    complex(dp), intent(in) :: square(:, :)

    !> Relative coordinates of the K-point
    real(dp), intent(in) :: kPoint(:)

    !> Weight of the K-point
    real(dp), intent(in) :: kweight

    !> Neighbour list for the atoms (First index from 0!)
    integer, intent(in) :: iNeighbour(0:, :)

    !> Nr. of neighbours for the atoms.
    integer, intent(in) :: nNeighbourSK(:)

    !> Maximal number of orbitals on an atom.
    integer, intent(in) :: mOrb

    !> Index of the cell translation vector for each atom.
    integer, intent(in) :: iCellVec(:)

    !> Relative coordinates of the cell translation vectors.
    real(dp), intent(in) :: cellVec(:, :)

    !> Atom offset for the square matrix
    integer, intent(in) :: iAtomStart(:)

    !> indexing array for the sparse Hamiltonian
    integer, intent(in) :: iSparseStart(0:, :)

    !> Mapping between image atoms and corresponding atom in the central cell.
    integer, intent(in) :: img2CentCell(:)

    complex(dp) :: phase
    integer :: nAtom
    integer :: iOrig, ii, jj, kk
    integer :: iNeigh, iBlock
    integer :: iOldVec, iVec
    integer :: iAtom1, iAtom2, iAtom2f
    integer :: nOrb1, nOrb2, nOrb
    real(dp) :: kPoint2p(3)
    complex(dp) :: tmpSqr(mOrb, mOrb)
  #:block DEBUG_CODE
    integer :: sizePrim
  #:endblock DEBUG_CODE

    nAtom = size(iNeighbour, dim=2)
    ! number of orbitals in a regular spin block
    nOrb = (iAtomStart(nAtom+1) - 1)

  #:block DEBUG_CODE
    sizePrim = size(primitive, dim=1)
  #:endblock DEBUG_CODE

    @:ASSERT(nAtom > 0)
    @:ASSERT(size(square, dim=1) == size(square, dim=2))
    @:ASSERT(size(square, dim=1) == 2 * nOrb )
    @:ASSERT(all(shape(kPoint) == [3]))
    @:ASSERT(all(shape(nNeighbourSK) == [nAtom]))
    @:ASSERT(size(iAtomStart) == nAtom + 1)
    @:ASSERT(kWeight > 0.0_dp)
    @:ASSERT(size(primitive, dim=2)==4)

    kPoint2p(:) = 2.0_dp * pi * kPoint

    ! sigma_I and sigma_z blocks
    do iBlock = 0, 1
      iOldVec = 0
      phase = 1.0_dp
      do iAtom1 = 1, nAtom
        ii = iAtomStart(iAtom1)
        nOrb1 = iAtomStart(iAtom1 + 1) - ii
        do iNeigh = 0, nNeighbourSK(iAtom1)
          iOrig = iSparseStart(iNeigh, iAtom1) + 1
          iAtom2 = iNeighbour(iNeigh, iAtom1)
          iAtom2f = img2CentCell(iAtom2)
          jj = iAtomStart(iAtom2f)
          @:ASSERT(jj >= ii)
          nOrb2 = iAtomStart(iAtom2f + 1) - jj
          iVec = iCellVec(iAtom2)
          if (iVec /= iOldVec) then
            phase = exp(cmplx(0, -1, dp) * dot_product(kPoint2p(:), cellVec(:, iVec)))
            iOldVec = iVec
          end if
          tmpSqr(1:nOrb2, 1:nOrb1) = square(jj+iBlock*nOrb:jj+nOrb2-1+iBlock*nOrb,&
              & ii+iBlock*nOrb:ii+nOrb1-1+iBlock*nOrb)
          ! Hermitian the on-site block before packing, as only one triangle usually supplied
          if (iAtom1 == iAtom2f) then
            do kk = 1, nOrb2
              tmpSqr(kk, kk+1:nOrb1) = conjg(tmpSqr(kk+1:nOrb1, kk))
            end do
          end if
          tmpSqr = phase*kWeight*tmpSqr
          @:ASSERT(sizePrim >= iOrig + nOrb1*nOrb2 - 1)
          primitive(iOrig : iOrig + nOrb1*nOrb2 - 1, 1) =&
              & primitive(iOrig : iOrig + nOrb1*nOrb2 - 1, 1)&
              & + 0.5_dp*reshape(aimag(tmpSqr(1:nOrb2, 1:nOrb1)), [nOrb1*nOrb2])
          primitive(iOrig : iOrig + nOrb1*nOrb2 - 1, 4) =&
              & primitive(iOrig : iOrig + nOrb1*nOrb2 - 1, 4)&
              & + real(1-2*iBlock, dp) *&
              & 0.5_dp*reshape(aimag(tmpSqr(1:nOrb2, 1:nOrb1)), [nOrb1*nOrb2])
        end do
      end do
    end do

    ! sigma_x and sigma_y blocks
    iOldVec = 0
    phase = 1.0_dp
    do iAtom1 = 1, nAtom
      ii = iAtomStart(iAtom1)
      nOrb1 = iAtomStart(iAtom1 + 1) - ii
      do iNeigh = 0, nNeighbourSK(iAtom1)
        iOrig = iSparseStart(iNeigh, iAtom1) + 1
        iAtom2 = iNeighbour(iNeigh, iAtom1)
        iAtom2f = img2CentCell(iAtom2)
        jj = iAtomStart(iAtom2f)
        @:ASSERT(jj >= ii)
        nOrb2 = iAtomStart(iAtom2f + 1) - jj
        iVec = iCellVec(iAtom2)
        if (iVec /= iOldVec) then
          phase = exp(cmplx(0, -1, dp)*dot_product(kPoint2p(:), cellVec(:, iVec)))
          iOldVec = iVec
        end if
        @:ASSERT(sizePrim >= iOrig + nOrb1*nOrb2 - 1)

        ! take the anti-Hermitian part of the block
        tmpSqr(1:nOrb2, 1:nOrb1) = 0.5_dp * ( square(nOrb+jj:nOrb+jj+nOrb2-1, ii:ii+nOrb1-1)&
            & + conjg(transpose(square(nOrb+ii:nOrb+ii+nOrb1-1, jj:jj+nOrb2-1))))
        tmpSqr = phase*kWeight*tmpSqr
        ! sigma_x : i * 1 = i use + imaginary part of block
        primitive(iOrig : iOrig + nOrb1*nOrb2 - 1, 2) =&
            & primitive(iOrig : iOrig + nOrb1 * nOrb2 - 1, 2)&
            & + reshape(aimag(tmpSqr(1:nOrb2, 1:nOrb1)), [nOrb1 * nOrb2])

        ! take the anti-Hermitian part of the block
        tmpSqr(1:nOrb2, 1:nOrb1) = 0.5_dp * (square(nOrb+jj:nOrb+jj+nOrb2-1, ii:ii+nOrb1-1)&
            & - conjg(transpose(square(nOrb+ii:nOrb+ii+nOrb1-1, jj:jj+nOrb2-1))))
        tmpSqr = phase*kWeight*tmpSqr
        ! sigma_y : i * i = -1 use - real part of block
        primitive(iOrig : iOrig + nOrb1*nOrb2 - 1, 3) =&
            & primitive(iOrig : iOrig + nOrb1 * nOrb2 - 1, 3)&
            & -reshape(real(tmpSqr(1:nOrb2, 1:nOrb1)), [nOrb1 * nOrb2])
      end do
    end do

  end subroutine packHSPauliImag_kpts


  !> Pack squared matrix in the sparse form (complex version without k-points).
  subroutine packhs_cmplx(prim, iPrim, square, iNeighbour, nNeighbourSK, mOrb, iAtomStart,&
      & iSparseStart, img2CentCell)

    !> Sparse matrix, real part
    real(dp), intent(inout) :: prim(:)

    !> Sparse matrix, imaginary part
    real(dp), intent(inout) :: iPrim(:)

    !> Hermitian matrix
    complex(dp), intent(in) :: square(:, :)

    !> Neighbour list for the atoms (First index from 0!)
    integer, intent(in) :: iNeighbour(0:, :)

    !> Nr. of neighbours for the atoms.
    integer, intent(in) :: nNeighbourSK(:)

    !> Maximal number of orbitals on an atom.
    integer, intent(in) :: mOrb

    !> Atom offset for the squared matrix
    integer, intent(in) :: iAtomStart(:)

    !> indexing array for the sparse Hamiltonian
    integer, intent(in) :: iSparseStart(0:, :)

    !> Mapping between image atoms and corresponding atom in the central cell.
    integer, intent(in) :: img2CentCell(:)

    integer :: nAtom
    integer :: iOrig, ii, jj, kk
    integer :: iNeigh
    integer :: iAtom1, iAtom2, iAtom2f
    integer :: nOrb1, nOrb2
    complex(dp) :: tmpSqr(mOrb, mOrb)

    nAtom = size(iNeighbour, dim=2)

    do iAtom1 = 1, nAtom
      ii = iAtomStart(iAtom1)
      nOrb1 = iAtomStart(iAtom1 + 1) - ii
      do iNeigh = 0, nNeighbourSK(iAtom1)
        iOrig = iSparseStart(iNeigh, iAtom1) + 1
        iAtom2 = iNeighbour(iNeigh, iAtom1)
        iAtom2f = img2CentCell(iAtom2)
        jj = iAtomStart(iAtom2f)
        nOrb2 = iAtomStart(iAtom2f + 1) - jj
        tmpSqr(1:nOrb2, 1:nOrb1) = square(jj:jj+nOrb2-1, ii:ii+nOrb1-1)

        ! Hermitian symmetrise the on-site block before packing, just in case
        if (iAtom1 == iAtom2f) then
          do kk = 1, nOrb2
            tmpSqr(kk, kk+1:nOrb1) = conjg(tmpSqr(kk+1:nOrb1, kk))
          end do
        end if

        prim(iOrig : iOrig + nOrb1*nOrb2 - 1) = prim(iOrig : iOrig + nOrb1 * nOrb2 - 1)&
            & + reshape(real(tmpSqr(1:nOrb2, 1:nOrb1),dp), [nOrb1 * nOrb2])
        iPrim(iOrig : iOrig + nOrb1*nOrb2 - 1) = iPrim(iOrig : iOrig + nOrb1 * nOrb2 - 1)&
            & + reshape(aimag(tmpSqr(1:nOrb2, 1:nOrb1)), [nOrb1 * nOrb2])
      end do
    end do

  end subroutine packhs_cmplx


  !> Pack only the charge (spin channel 1) part of a 2 component matrix
  subroutine packHSPauliERho(primitive, square, iNeighbour, nNeighbourSK, mOrb, iAtomStart,&
      & iSparseStart, img2CentCell)

    !> Sparse matrix
    real(dp), intent(inout) :: primitive(:)

    !> Square form matrix
    complex(dp), intent(in) :: square(:, :)

    !> Neighbour list for the atoms (First index from 0!)
    integer, intent(in) :: iNeighbour(0:, :)

    !> Nr. of neighbours for the atoms.
    integer, intent(in) :: nNeighbourSK(:)

    !> Maximal number of orbitals on an atom.
    integer, intent(in) :: mOrb

    !> Atom offset for the square matrix
    integer, intent(in) :: iAtomStart(:)

    !> indexing array for the sparse Hamiltonian
    integer, intent(in) :: iSparseStart(0:, :)

    !> Mapping between image atoms and corresponding atom in the central cell.
    integer, intent(in) :: img2CentCell(:)

    integer :: nAtom
    integer :: iOrig, ii, jj, kk
    integer :: iNeigh, iBlock
    integer :: iAtom1, iAtom2, iAtom2f
    integer :: nOrb1, nOrb2, nOrb
    complex(dp) :: tmpSqr(mOrb, mOrb)
  #:block DEBUG_CODE
    integer :: sizePrim
  #:endblock DEBUG_CODE

    nAtom = size(iNeighbour, dim=2)
    ! number of orbitals in a regular spin block
    nOrb = (iAtomStart(nAtom+1) - 1)

  #:block DEBUG_CODE
    sizePrim = size(primitive, dim=1)
  #:endblock DEBUG_CODE

    @:ASSERT(nAtom > 0)
    @:ASSERT(size(square, dim=1) == size(square, dim=2))
    @:ASSERT(size(square, dim=1) == 2 * nOrb )
    @:ASSERT(all(shape(nNeighbourSK) == [nAtom]))
    @:ASSERT(size(iAtomStart) == nAtom + 1)

    do iBlock = 0, 1
      do iAtom1 = 1, nAtom
        ii = iAtomStart(iAtom1)
        nOrb1 = iAtomStart(iAtom1 + 1) - ii
        do iNeigh = 0, nNeighbourSK(iAtom1)
          iOrig = iSparseStart(iNeigh, iAtom1) + 1
          iAtom2 = iNeighbour(iNeigh, iAtom1)
          iAtom2f = img2CentCell(iAtom2)
          jj = iAtomStart(iAtom2f)
          @:ASSERT(jj >= ii)
          nOrb2 = iAtomStart(iAtom2f+1) - jj
          tmpSqr(1:nOrb2, 1:nOrb1) = square(jj+iBlock*nOrb:jj+nOrb2-1+iBlock*nOrb,&
              & ii+iBlock*nOrb:ii+nOrb1-1+iBlock*nOrb)
          ! Symmetrize the on-site block before packing, as only one triangle supplied
          if (iAtom1 == iAtom2f) then
            do kk = 1, nOrb2
              tmpSqr(kk, kk+1:nOrb1) = conjg(tmpSqr(kk+1:nOrb1, kk))
            end do
          end if
          @:ASSERT(sizePrim >= iOrig + nOrb1*nOrb2 - 1)
          primitive(iOrig : iOrig + nOrb1*nOrb2 - 1) = primitive(iOrig : iOrig + nOrb1*nOrb2 - 1)&
              & + reshape(real(tmpSqr(1:nOrb2, 1:nOrb1)), (/nOrb1*nOrb2/))
        end do
      end do
    end do

  end subroutine packHSPauliERho


  !> Pack square matrix in the sparse form (real version).
  subroutine packHSPauliERho_kpts(primitive, square, kPoint, kWeight, iNeighbour, nNeighbourSK,&
      & mOrb, iCellVec, cellVec, iAtomStart, iSparseStart, img2CentCell)

    !> Sparse matrix
    real(dp), intent(inout) :: primitive(:)

    !> Square form matrix
    complex(dp), intent(in) :: square(:, :)

    !> Relative coordinates of the K-point where the sparse matrix should be unfolded.
    real(dp), intent(in) :: kPoint(:)

    !> weight for the k-point
    real(dp), intent(in) :: kweight

    !> Neighbour list for the atoms (First index from 0!)
    integer, intent(in) :: iNeighbour(0:, :)

    !> Nr. of neighbours for the atoms.
    integer, intent(in) :: nNeighbourSK(:)

    !> Maximal number of orbitals on an atom.
    integer, intent(in) :: mOrb

    !> Index of the cell translation vector for each atom.
    integer, intent(in) :: iCellVec(:)

    !> Relative coordinates of the cell translation vectors.
    real(dp), intent(in) :: cellVec(:, :)

    !> Atom offset for the square matrix
    integer, intent(in) :: iAtomStart(:)

    !> indexing array for the sparse Hamiltonian
    integer, intent(in) :: iSparseStart(0:, :)

    !> Mapping between image atoms and corresponding atom in the central cell.
    integer, intent(in) :: img2CentCell(:)

    complex(dp) :: phase
    integer :: nAtom
    integer :: iOrig, ii, jj, kk
    integer :: iNeigh, iBlock
    integer :: iOldVec, iVec
    integer :: iAtom1, iAtom2, iAtom2f
    integer :: nOrb1, nOrb2, nOrb
    real(dp) :: kPoint2p(3)
    complex(dp) :: tmpSqr(mOrb, mOrb)
  #:block DEBUG_CODE
    integer :: sizePrim
  #:endblock DEBUG_CODE

    nAtom = size(iNeighbour, dim=2)
    ! number of orbitals in a regular spin block
    nOrb = (iAtomStart(nAtom+1) - 1)

  #:block DEBUG_CODE
    sizePrim = size(primitive, dim=1)
  #:endblock DEBUG_CODE

    @:ASSERT(nAtom > 0)
    @:ASSERT(size(square, dim=1) == size(square, dim=2))
    @:ASSERT(size(square, dim=1) == 2 * nOrb )
    @:ASSERT(all(shape(kPoint) == [3]))
    @:ASSERT(all(shape(nNeighbourSK) == [nAtom]))
    @:ASSERT(size(iAtomStart) == nAtom + 1)
    @:ASSERT(kWeight > 0.0_dp)

    kPoint2p(:) = 2.0_dp * pi * kPoint

    do iBlock = 0, 1
      iOldVec = 0
      phase = 1.0_dp
      do iAtom1 = 1, nAtom
        ii = iAtomStart(iAtom1)
        nOrb1 = iAtomStart(iAtom1 + 1) - ii
        do iNeigh = 0, nNeighbourSK(iAtom1)
          iOrig = iSparseStart(iNeigh, iAtom1) + 1
          iAtom2 = iNeighbour(iNeigh, iAtom1)
          iAtom2f = img2CentCell(iAtom2)
          jj = iAtomStart(iAtom2f)
          @:ASSERT(jj >= ii)
          nOrb2 = iAtomStart(iAtom2f + 1) - jj
          iVec = iCellVec(iAtom2)
          if (iVec /= iOldVec) then
            phase = exp(cmplx(0, -1, dp) * dot_product(kPoint2p(:), cellVec(:, iVec)))
            iOldVec = iVec
          end if
          tmpSqr(1:nOrb2, 1:nOrb1) = square(jj+iBlock*nOrb:jj+nOrb2-1+iBlock*nOrb,&
              & ii+iBlock*nOrb:ii+nOrb1-1+iBlock*nOrb)
          ! Symmetrize the on-site block before packing, as only one triangle supplied
          if (iAtom1 == iAtom2f) then
            do kk = 1, nOrb2
              tmpSqr(kk, kk+1:nOrb1) = conjg(tmpSqr(kk+1:nOrb1, kk))
            end do
          end if
          @:ASSERT(sizePrim >= iOrig + nOrb1*nOrb2 - 1)
          primitive(iOrig : iOrig + nOrb1*nOrb2 - 1) =&
              & primitive(iOrig : iOrig + nOrb1*nOrb2 - 1) + kWeight *&
              & reshape(real(phase*tmpSqr(1:nOrb2, 1:nOrb1)), [nOrb1*nOrb2])
        end do
      end do
    end do

  end subroutine packHSPauliERho_kpts


#:if WITH_SCALAPACK
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!! Scalapack routines
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


  !> Unpacks sparse H or S into dense (real, blacs).
  !>
  !> Note: In contrast to the serial routines, both triangles of the resulting matrix are filled.
  !>
  subroutine unpackHSRealBlacs(myBlacs, orig, iNeighbour, nNeighbourSK, iSparseStart, img2CentCell,&
      & desc, square)

    !> BLACS matrix descriptor
    type(TBlacsEnv), intent(in) :: myBlacs

    !> sparse matrix
    real(dp), intent(in) :: orig(:)

    !> Neighbour list for each atom (First index from 0!)
    integer, intent(in) :: iNeighbour(0:, :)

    !> Nr. of neighbours for each atom (incl. itself).
    integer, intent(in) :: nNeighbourSK(:)

    !> indexing array for the sparse Hamiltonian
    integer, intent(in) :: iSparseStart(0:, :)

    !> Map from images of atoms to central cell atoms
    integer, intent(in) :: img2CentCell(:)

    !> Dense matrix description
    type(TDenseDescr), intent(in) :: desc

    !> dense matrix part of distributed whole
    real(dp), intent(out) :: square(:, :)

    integer :: nAtom
    integer :: iOrig, ii, jj, nOrb1, nOrb2
    integer :: iNeigh
    integer :: iAtom1, iAtom2, iAtom2f

    nAtom = size(iNeighbour, dim=2)

    @:ASSERT(nAtom > 0)
    @:ASSERT(size(nNeighbourSK) == nAtom)
    @:ASSERT(size(desc%iAtomStart) == nAtom + 1)

    square(:, :) = 0.0_dp

    do iAtom1 = 1, nAtom
      ii = desc%iAtomStart(iAtom1)
      nOrb1 = desc%iAtomStart(iAtom1 + 1) - ii
      do iNeigh = 0, nNeighbourSK(iAtom1)
        iOrig = iSparseStart(iNeigh, iAtom1) + 1
        iAtom2 = iNeighbour(iNeigh, iAtom1)
        iAtom2f = img2CentCell(iAtom2)
        jj = desc%iAtomStart(iAtom2f)
        nOrb2 = desc%iAtomStart(iAtom2f + 1) - jj
        call scalafx_addl2g(myBlacs%orbitalGrid,&
            & reshape(orig(iOrig : iOrig + nOrb1 * nOrb2 - 1), [nOrb2, nOrb1]),&
            & desc%blacsOrbSqr, jj, ii, square)
        if (iAtom1 /= iAtom2f) then
          call scalafx_addl2g(myBlacs%orbitalGrid,&
              & transpose(reshape(orig(iOrig : iOrig + nOrb1 * nOrb2 - 1), [nOrb2, nOrb1])),&
              & desc%blacsOrbSqr, ii, jj, square)
        end if
      end do
    end do

  end subroutine unpackHSRealBlacs


  !> Unpacks sparse H into distributed dense matrix (complex)
  !>
  !> Note: In contrast to the serial routines, both triangles of the resulting matrix are filled.
  !>
  subroutine unpackHSCplxBlacs(myBlacs, orig, kPoint, iNeighbour, nNeighbourSK, iCellVec, cellVec,&
      & iSparseStart, img2CentCell, desc, square)

    !> BLACS matrix descriptor
    type(TBlacsEnv), intent(in) :: myBlacs

    !> sparse matrix
    real(dp), intent(in) :: orig(:)

    !> Relative coordinates of the K-point where the sparse matrix should be unfolded.
    real(dp), intent(in) :: kPoint(:)

    !> Neighbour list for each atom (First index from 0!)
    integer, intent(in) :: iNeighbour(0:, :)

    !> Nr. of neighbours for each atom (incl. itself).
    integer, intent(in) :: nNeighbourSK(:)

    !> Index of the cell translation vector for each atom.
    integer, intent(in) :: iCellVec(:)

    !> Relative coordinates of the cell translation vectors.
    real(dp), intent(in) :: cellVec(:, :)

    !> indexing array for the sparse Hamiltonian
    integer, intent(in) :: iSparseStart(0:, :)

    !> Mapping between image atoms and corresponding atom in the central cell.
    integer, intent(in) :: img2CentCell(:)

    !> Dense matrix description
    type(TDenseDescr), intent(in) :: desc

    !> dense matrix part of distributed whole
    complex(dp), intent(out) :: square(:, :)

    complex(dp) :: phase
    integer :: nAtom
    integer :: iOrig, nOrb1, nOrb2, ii, jj
    integer :: iNeigh
    integer :: iOldVec, iVec
    integer :: iAtom1, iAtom2, iAtom2f
    real(dp) :: kPoint2p(3)

    nAtom = size(iNeighbour, dim=2)

    @:ASSERT(nAtom > 0)
    @:ASSERT(size(kPoint) == 3)
    @:ASSERT(size(nNeighbourSK) == nAtom)
    @:ASSERT(size(desc%iAtomStart) == nAtom + 1)

    square(:, :) = cmplx(0, 0, dp)
    kPoint2p(:) = 2.0_dp * pi * kPoint
    iOldVec = 0
    phase = 1.0_dp
    do iAtom1 = 1, nAtom
      ii = desc%iAtomStart(iAtom1)
      nOrb1 = desc%iAtomStart(iAtom1 + 1) - ii
      do iNeigh = 0, nNeighbourSK(iAtom1)
        iOrig = iSparseStart(iNeigh, iAtom1) + 1
        iAtom2 = iNeighbour(iNeigh, iAtom1)
        iAtom2f = img2CentCell(iAtom2)
        jj = desc%iAtomStart(iAtom2f)
        nOrb2 = desc%iAtomStart(iAtom2f + 1) - jj
        iVec = iCellVec(iAtom2)
        if (iVec /= iOldVec) then
          phase = exp((0.0_dp, 1.0_dp) * dot_product(kPoint2p, cellVec(:, iVec)))
          iOldVec = iVec
        end if
        call scalafx_addl2g(myBlacs%orbitalGrid, phase&
            & * reshape(orig(iOrig : iOrig + nOrb1 * nOrb2 - 1), [nOrb2, nOrb1]),&
            & desc%blacsOrbSqr, jj, ii, square)
        if (iAtom1 /= iAtom2f) then
          call scalafx_addl2g(myBlacs%orbitalGrid, transpose(conjg(phase&
              & * reshape(orig(iOrig : iOrig + nOrb1 * nOrb2 - 1), [nOrb2, nOrb1]))),&
              & desc%blacsOrbSqr, ii, jj, square)
        end if
      end do
    end do

  end subroutine unpackHSCplxBlacs


  !> Unpacks sparse Hamiltonian to square form (Pauli-type Hamiltonian).
  !>
  !> Note: In contrast to the serial routines, both triangles of the resulting matrix are filled.
  !>
  subroutine unpackHPauliBlacs(myBlacs, orig, kPoint, iNeighbour, nNeighbourSK, iCellVec, cellVec,&
      & iSparseStart, img2CentCell, mOrb, desc, square, iorig)

    !> BLACS matrix descriptor
    type(TBlacsEnv), intent(in) :: myBlacs

    !> sparse matrix
    real(dp), intent(in) :: orig(:,:)

    !> Relative coordinates of the K-point where the sparse matrix should be unfolded.
    real(dp), intent(in) :: kPoint(:)

    !> Neighbour list for the atoms (First index from 0!)
    integer, intent(in) :: iNeighbour(0:, :)

    !> Nr. of neighbours for the atoms.
    integer, intent(in) :: nNeighbourSK(:)

    !> Index of the cell translation vector for each atom.
    integer, intent(in) :: iCellVec(:)

    !> Relative coordinates of the cell translation vectors.
    real(dp), intent(in) :: cellVec(:, :)

    !> indexing array for the sparse Hamiltonian
    integer, intent(in) :: iSparseStart(0:, :)

    !> Mapping between image atoms and corresponding atom in the central cell.
    integer, intent(in) :: img2CentCell(:)

    !> Maximal number of orbitals on an atom.
    integer, intent(in) :: mOrb

    !> Dense matrix description
    type(TDenseDescr), intent(in) :: desc

    !> local part of distributed whole dense matrix
    complex(dp), intent(out) :: square(:,:)

    !> imaginary part of sparse matrix
    real(dp), intent(in), optional :: iorig(:,:)

    square(:,:) = cmplx(0, 0, dp)
    call unpackHPauliBlacsHelper(myBlacs, orig, kPoint, iNeighbour, nNeighbourSK, iCellVec,&
        & cellVec, iSparseStart, img2CentCell, mOrb, cmplx(1, 0, dp), cmplx(1, 0, dp), desc, square)
    if (present(iorig)) then
      call unpackHPauliBlacsHelper(myBlacs, iorig, kPoint, iNeighbour, nNeighbourSK, iCellVec,&
          & cellVec, iSparseStart, img2CentCell, mOrb, cmplx(0, 1, dp), cmplx(-1, 0, dp), desc,&
          & square)
    end if

  end subroutine unpackHPauliBlacs


  !> Helper routine for unpacking into Pauli-type Hamiltonians.
  !!
  !! The routine creates both triangle of the 2x2 Pauli Hamiltonian
  !! 1*orig(:, 1) + sigma1*orig(:, 2) + sigma2*orig(:, 3) + sigma3*orig(:, 4).
  !!
  subroutine unpackHPauliBlacsHelper(myBlacs, orig, kPoint, iNeighbour, nNeighbourSK, iCellVec,&
      & cellVec, iSparseStart, img2CentCell, mOrb, imagPrefac, hermPrefac, desc, square)

    !> BLACS matrix descriptor
    type(TBlacsEnv), intent(in) :: myBlacs

    !> sparse matrix to unpack
    real(dp), intent(in) :: orig(:, :)

    !> The k-point at which to unpack
    real(dp), intent(in) :: kPoint(:)

    !> Neighbour list for the atoms (First index from 0!)
    integer, intent(in) :: iNeighbour(0:, :)

    !> Nr. of neighbours for the atoms.
    integer, intent(in) :: nNeighbourSK(:)

    !> Index of the cell translation vector for each atom.
    integer, intent(in) :: iCellVec(:)

    !> Relative coordinates of the cell translation vectors.
    real(dp), intent(in) :: cellVec(:, :)

    !> indexing array for the sparse Hamiltonian
    integer, intent(in) :: iSparseStart(0:, :)

    !> Mapping between image atoms and corresponding atom in the central cell.
    integer, intent(in) :: img2CentCell(:)

    !> Maximal number of orbitals on an atom.
    integer, intent(in) :: mOrb

    !> prefactor for the imaginary related part of a block
    complex(dp), intent(in) :: imagPrefac

    !> prefactor for the hermitian related part of a block
    complex(dp), intent(in) :: hermPrefac

    !> Dense matrix description
    type(TDenseDescr), intent(in) :: desc

    !> local part of distributed whole dense matrix to add to
    complex(dp), intent(inout) :: square(:, :)

    integer :: iAtom1, iAtom2, iAtom2f, nAtom, nOrb1, nOrb2, nOrb
    integer :: iNeigh, iOrig, ii, jj, kk, iOldVec, iVec
    complex(dp), target :: tmpSqr(mOrb, mOrb, 4)
    complex(dp), pointer :: ptmp(:, :, :)
    real(dp) :: kPoint2p(3)
    complex(dp) :: phase

    nAtom = size(nNeighbourSK)
    nOrb = desc%iAtomStart(nAtom + 1) - 1
    kPoint2p(:) = 2.0_dp * pi * kPoint
    iOldVec = 0
    phase = 1.0_dp
    do iAtom1 = 1, nAtom
      ii = desc%iAtomStart(iAtom1)
      nOrb1 = desc%iAtomStart(iAtom1 + 1) - ii
      do iNeigh = 0, nNeighbourSK(iAtom1)
        iOrig = iSparseStart(iNeigh, iAtom1) + 1
        iAtom2 = iNeighbour(iNeigh, iAtom1)
        iAtom2f = img2CentCell(iAtom2)
        jj = desc%iAtomStart(iAtom2f)
        nOrb2 = desc%iAtomStart(iAtom2f + 1) - jj
        iVec = iCellVec(iAtom2)
        if (iVec /= iOldVec) then
          phase = exp(imag * dot_product(kPoint2p, cellVec(:, iVec)))
          iOldVec = iVec
        end if
        ptmp => tmpSqr(1:nOrb2, 1:nOrb1, :)
        ptmp(:, :, :) = 0.5_dp * phase&
            & * reshape(orig(iOrig:iOrig+nOrb1*nOrb2-1, :), [nOrb2, nOrb1, 4])
        ! up-up and down-down components
        call scalafx_addl2g(myBlacs%orbitalGrid, imagPrefac * (ptmp(:, :, 1) + ptmp(:, :, 4)),&
            & desc%blacsOrbSqr, jj, ii, square)
        call scalafx_addl2g(myBlacs%orbitalGrid, imagPrefac * (ptmp(:, :, 1) - ptmp(:, :, 4)),&
            & desc%blacsOrbSqr, jj + nOrb, ii + nOrb, square)
        if (iAtom1 /= iAtom2f) then
          call scalafx_addl2g(myBlacs%orbitalGrid,&
              & transpose(conjg(imagPrefac * (ptmp(:, :, 1) + ptmp(:, :, 4)))),&
              & desc%blacsOrbSqr, ii, jj, square)
          call scalafx_addl2g(myBlacs%orbitalGrid,&
              & transpose(conjg(imagPrefac * (ptmp(:, :, 1) - ptmp(:, :, 4)))),&
              & desc%blacsOrbSqr, ii + nOrb, jj + nOrb, square)
        end if
        ! down-up component
        ! also upper triangle of the down-up component must be filled
        if (iAtom1 == iAtom2f) then
          ! symmetrize/antisymmetrize onsite block
          do kk = 1, nOrb1
            ptmp(kk, kk+1:, 2:3) = hermPrefac * conjg(ptmp(kk+1:, kk, 2:3))
          end do
          call scalafx_addl2g(myBlacs%orbitalGrid,&
              & imagPrefac * (ptmp(:, :, 2) + imag * ptmp(:, :, 3)), desc%blacsOrbSqr,&
              & jj + nOrb, ii, square)
          ! Other triangle
          call scalafx_addl2g(myBlacs%orbitalGrid,&
              & -hermPrefac*transpose(conjg(imagPrefac * (ptmp(:, :, 2) + imag * ptmp(:, :, 3)))),&
              & desc%blacsOrbSqr, ii, jj + nOrb, square)
        else
          call scalafx_addl2g(myBlacs%orbitalGrid,&
              & imagPrefac * (ptmp(:, :, 2) + imag * ptmp(:, :, 3)), desc%blacsOrbSqr,&
              & jj + nOrb, ii, square)
          call scalafx_addl2g(myBlacs%orbitalGrid,&
              & -hermPrefac*transpose(conjg(imagPrefac * (ptmp(:, :, 2) + imag * ptmp(:, :, 3)))),&
              & desc%blacsOrbSqr, ii, jj + nOrb, square)

          call scalafx_addl2g(myBlacs%orbitalGrid, imagPrefac * hermPrefac&
              & * conjg(transpose(ptmp(:, :, 2) - imag * ptmp(:, :, 3))), desc%blacsOrbSqr,&
              & ii + nOrb, jj, square)
          call scalafx_addl2g(myBlacs%orbitalGrid, conjg(imagPrefac * hermPrefac)&
              & * (ptmp(:, :, 2) - imag * ptmp(:, :, 3))&
              & , desc%blacsOrbSqr, jj, ii + nOrb, square)
        end if
      end do
    end do

  end subroutine unpackHPauliBlacsHelper


  !> Unpacking the overlap for Pauli-type matrices.
  !>
  !> Note: In contrast to the serial routines, both triangles of the resulting matrix are filled.
  !>
  subroutine unpackSPauliBlacs(myBlacs, orig, kPoint, iNeighbour, nNeighbourSK, iCellVec, cellVec,&
      & iSparseStart, img2CentCell, mOrb, desc, square)

    !> BLACS matrix descriptor
    type(TBlacsEnv), intent(in) :: myBlacs

    !> sparse matrix to unpack
    real(dp), intent(in) :: orig(:)

    !> The k-point at which to unpack
    real(dp), intent(in) :: kPoint(:)

    !> Neighbour list for the atoms (First index from 0!)
    integer, intent(in) :: iNeighbour(0:, :)

    !> Nr. of neighbours for the atoms.
    integer, intent(in) :: nNeighbourSK(:)

    !> Index of the cell translation vector for each atom.
    integer, intent(in) :: iCellVec(:)

    !> Relative coordinates of the cell translation vectors.
    real(dp), intent(in) :: cellVec(:, :)

    !> indexing array for the sparse Hamiltonian
    integer, intent(in) :: iSparseStart(0:, :)

    !> Mapping between image atoms and corresponding atom in the central cell.
    integer, intent(in) :: img2CentCell(:)

    !> Maximal number of orbitals on an atom.
    integer, intent(in) :: mOrb

    !> Dense matrix description
    type(TDenseDescr), intent(in) :: desc

    !> Local part of of distributed whole dense matrix
    complex(dp), intent(out) :: square(:, :)

    integer :: iAtom1, iAtom2, iAtom2f, nAtom, nOrb1, nOrb2, nOrb
    integer :: iNeigh, iOrig, ii, jj, iOldVec, iVec
    complex(dp), target :: tmpSqr(mOrb, mOrb)
    complex(dp), pointer :: ptmp(:, :)
    real(dp) :: kPoint2p(3)
    complex(dp) :: phase

    square(:, :) = (0.0_dp, 0.0_dp)
    nAtom = size(nNeighbourSK)
    nOrb = desc%iAtomStart(nAtom + 1) - 1
    kPoint2p(:) = 2.0_dp * pi * kPoint
    iOldVec = 0
    phase = 1.0_dp
    do iAtom1 = 1, nAtom
      ii = desc%iAtomStart(iAtom1)
      nOrb1 = desc%iAtomStart(iAtom1 + 1) - ii
      do iNeigh = 0, nNeighbourSK(iAtom1)
        iOrig = iSparseStart(iNeigh, iAtom1) + 1
        iAtom2 = iNeighbour(iNeigh, iAtom1)
        iAtom2f = img2CentCell(iAtom2)
        jj = desc%iAtomStart(iAtom2f)
        nOrb2 = desc%iAtomStart(iAtom2f + 1) - jj
        iVec = iCellVec(iAtom2)
        if (iVec /= iOldVec) then
          phase = exp(imag * dot_product(kPoint2p, cellVec(:, iVec)))
          iOldVec = iVec
        end if
        ptmp => tmpSqr(1:nOrb2, 1:nOrb1)
        ptmp(:, :) = phase * reshape(orig(iOrig:iOrig+nOrb1*nOrb2-1), [nOrb2, nOrb1])
        ! up-up component
        call scalafx_addl2g(myBlacs%orbitalGrid, ptmp, desc%blacsOrbSqr, jj, ii, square)
        ! down-down component
        call scalafx_addl2g(myBlacs%orbitalGrid, ptmp, desc%blacsOrbSqr, jj + nOrb, ii + nOrb,&
            & square)
        if (iAtom1 /= iAtom2f) then
          call scalafx_addl2g(myBlacs%orbitalGrid, transpose(conjg(ptmp)), desc%blacsOrbSqr, ii,&
              & jj, square)
          call scalafx_addl2g(myBlacs%orbitalGrid, transpose(conjg(ptmp)), desc%blacsOrbSqr,&
              & ii + nOrb, jj + nOrb, square)
        end if
      end do
    end do

  end subroutine unpackSPauliBlacs


  !> Packs distributed dense real matrix into sparse form (real).
  subroutine packRhoRealBlacs(myBlacs, desc, square, iNeighbour, nNeighbourSK, mOrb, iSparseStart,&
      & img2CentCell, primitive)

    !> BLACS matrix descriptor
    type(TBlacsEnv), intent(in) :: myBlacs

    !> Dense matrix description
    type(TDenseDescr), intent(in) :: desc

    !> distributed dense matrix to pack
    real(dp), intent(in) :: square(:, :)

    !> Neighbour list for the atoms (First index from 0!)
    integer, intent(in) :: iNeighbour(0:, :)

    !> Nr. of neighbours for the atoms.
    integer, intent(in) :: nNeighbourSK(:)

    !> Maximal number of orbitals on an atom.
    integer, intent(in) :: mOrb

    !> indexing array for the sparse Hamiltonian
    integer, intent(in) :: iSparseStart(0:, :)

    !> Mapping between image atoms and corresponding atom in the central cell.
    integer, intent(in) :: img2CentCell(:)

    !> sparse matrix to add this contribution into
    real(dp), intent(inout) :: primitive(:)

    integer :: nAtom
    integer :: iOrig, ii, jj, kk
    integer :: iNeigh
    integer :: iAtom1, iAtom2, iAtom2f
    integer :: nOrb1, nOrb2
    real(dp) :: tmpSqr(mOrb, mOrb)
  #:block DEBUG_CODE
    integer :: sizePrim
  #:endblock DEBUG_CODE

    nAtom = size(iNeighbour, dim=2)
  #:block DEBUG_CODE
    sizePrim = size(primitive)
  #:endblock DEBUG_CODE
    @:ASSERT(nAtom > 0)
    @:ASSERT(size(nNeighbourSK) == nAtom)

    do iAtom1 = 1, nAtom
      ii = desc%iAtomStart(iAtom1)
      nOrb1 = desc%iAtomStart(iAtom1 + 1) - ii
      do iNeigh = 0, nNeighbourSK(iAtom1)
        iOrig = iSparseStart(iNeigh, iAtom1) + 1
        iAtom2 = iNeighbour(iNeigh, iAtom1)
        iAtom2f = img2CentCell(iAtom2)
        jj = desc%iAtomStart(iAtom2f)
        nOrb2 = desc%iAtomStart(iAtom2f + 1) - jj
        call scalafx_cpg2l(myBlacs%orbitalGrid, desc%blacsOrbSqr, jj, ii, square,&
            & tmpSqr(1:nOrb2, 1:nOrb1))

        ! Symmetrize the on-site block before packing, just in case
        if (iAtom1 == iAtom2f) then
          do kk = 1, nOrb2
            tmpSqr(kk, kk+1:nOrb1) = tmpSqr(kk + 1 : nOrb1, kk)
          end do
        end if

        primitive(iOrig : iOrig + nOrb1 * nOrb2 - 1) =&
            & primitive(iOrig : iOrig + nOrb1 * nOrb2 - 1)&
            & + reshape(tmpSqr(1:nOrb2, 1:nOrb1), [nOrb1 * nOrb2])
      end do
    end do

  end subroutine packRhoRealBlacs


  !> Packs distributed dense matrix into sparse form (complex).
  subroutine packRhoCplxBlacs(myblacs, desc, square, kPoint, kWeight, iNeighbour, nNeighbourSK,&
      & mOrb, iCellVec, cellVec, iSparseStart, img2CentCell, primitive)

    !> BLACS matrix descriptor
    type(TBlacsEnv), intent(in) :: myBlacs

    !> Dense matrix description
    type(TDenseDescr), intent(in) :: desc

    !> Distributed dense matrix to pack
    complex(dp), intent(in) :: square(:, :)

    !> The k-point at which to pack
    real(dp), intent(in) :: kPoint(:)

    !> weight for this k-point
    real(dp), intent(in) :: kWeight

    !> Neighbour list for the atoms (First index from 0!)
    integer, intent(in) :: iNeighbour(0:, :)

    !> Nr. of neighbours for the atoms.
    integer, intent(in) :: nNeighbourSK(:)

    !> Maximal number of orbitals on an atom.
    integer, intent(in) :: mOrb

    !> Index of the cell translation vector for each atom.
    integer, intent(in) :: iCellVec(:)

    !> Relative coordinates of the cell translation vectors.
    real(dp), intent(in) :: cellVec(:, :)

    !> indexing array for the sparse Hamiltonian
    integer, intent(in) :: iSparseStart(0:, :)

    !> Mapping between image atoms and corresponding atom in the central cell.
    integer, intent(in) :: img2CentCell(:)

    !> sparse matrix to add this contribution into
    real(dp), intent(inout) :: primitive(:)

    complex(dp) :: phase
    integer :: nAtom
    integer :: iOrig, ii, jj, kk
    integer :: iNeigh
    integer :: iOldVec, iVec
    integer :: iAtom1, iAtom2, iAtom2f
    integer :: nOrb1, nOrb2
    real(dp) :: kPoint2p(3)
    complex(dp) :: tmpSqr(mOrb, mOrb)
  #:block DEBUG_CODE
    integer :: sizePrim
  #:endblock DEBUG_CODE

    nAtom = size(iNeighbour, dim=2)
  #:block DEBUG_CODE
    sizePrim = size(primitive)
  #:endblock DEBUG_CODE

    @:ASSERT(nAtom > 0)
    @:ASSERT(size(kPoint) == 3)
    @:ASSERT(size(nNeighbourSK) == nAtom)
    @:ASSERT(kWeight > 0.0_dp)
    @:ASSERT(size(desc%iAtomStart) == nAtom + 1)

    kPoint2p(:) = 2.0_dp * pi * kPoint
    iOldVec = 0
    phase = 1.0_dp
    do iAtom1 = 1, nAtom
      ii = desc%iAtomStart(iAtom1)
      nOrb1 = desc%iAtomStart(iAtom1 + 1) - ii
      do iNeigh = 0, nNeighbourSK(iAtom1)
        iOrig = iSparseStart(iNeigh, iAtom1) + 1
        iAtom2 = iNeighbour(iNeigh, iAtom1)
        iAtom2f = img2CentCell(iAtom2)
        jj = desc%iAtomStart(iAtom2f)
        nOrb2 = desc%iAtomStart(iAtom2f + 1) - jj
        iVec = iCellVec(iAtom2)
        if (iVec /= iOldVec) then
          phase = exp(cmplx(0, -1, dp) * dot_product(kPoint2p(:), cellVec(:, iVec)))
          iOldVec = iVec
        end if
        call scalafx_cpg2l(myBlacs%orbitalGrid, desc%blacsOrbSqr, jj, ii, square,&
            & tmpSqr(1:nOrb2, 1:nOrb1))

        ! Hermitian the on-site block before packing, just in case
        if (iAtom1 == iAtom2f) then
          do kk = 1, nOrb2
            tmpSqr(kk, kk + 1 : nOrb1) = conjg(tmpSqr(kk + 1 : nOrb1, kk))
          end do
        end if

        primitive(iOrig : iOrig + nOrb1*nOrb2 - 1) =&
            & primitive(iOrig : iOrig + nOrb1*nOrb2 - 1)&
            & + kWeight * real(phase * reshape(tmpSqr(1:nOrb2, 1:nOrb1), [nOrb1 * nOrb2]), dp)
      end do
    end do

  end subroutine packRhoCplxBlacs


  !> Pack square dense matrix into the sparse form (complex Pauli version).
  subroutine packRhoPauliBlacs(myBlacs, desc, square, kPoint, kWeight, iNeighbour, nNeighbourSK,&
      & mOrb, iCellVec, cellVec, iSparseStart, img2CentCell, primitive, iprimitive)

    !> BLACS matrix descriptor
    type(TBlacsEnv), intent(in) :: myBlacs

    !> Dense matrix description
    type(TDenseDescr), intent(in) :: desc

    !> dense distributed matrix to pack
    complex(dp), intent(in) :: square(:, :)

    !> The k-point at which to pack
    real(dp), intent(in) :: kPoint(:)

    !> weight for this k-point
    real(dp), intent(in) :: kWeight

    !> Neighbour list for the atoms (First index from 0!)
    integer, intent(in) :: iNeighbour(0:, :)

    !> Nr. of neighbours for the atoms.
    integer, intent(in) :: nNeighbourSK(:)

    !> Maximal number of orbitals on an atom.
    integer, intent(in) :: mOrb

    !> Index of the cell translation vector for each atom.
    integer, intent(in) :: iCellVec(:)

    !> Relative coordinates of the cell translation vectors.
    real(dp), intent(in) :: cellVec(:, :)

    !> indexing array for the sparse Hamiltonian
    integer, intent(in) :: iSparseStart(0:, :)

    !> Mapping between image atoms and corresponding atom in the central cell.
    integer, intent(in) :: img2CentCell(:)

    !> sparse matrix to pack into
    real(dp), intent(inout) :: primitive(:, :)

    !> Imaginary part of sparse matrix
    real(dp), intent(inout), optional :: iprimitive(:, :)

    call packRhoPauliBlacsHelper(myBlacs, desc, square, kPoint, kWeight, iNeighbour, nNeighbourSK,&
        & mOrb, iCellVec, cellVec, iSparseStart, img2CentCell, cmplx(1, 0, dp), .true., primitive)
    if (present(iprimitive)) then
      call packRhoPauliBlacsHelper(myBlacs, desc, square, kPoint, kWeight, iNeighbour,&
          & nNeighbourSK, mOrb, iCellVec, cellVec, iSparseStart, img2CentCell, cmplx(0, -1, dp),&
          & .false., iprimitive)
    end if

  end subroutine packRhoPauliBlacs


  !> Helper routine for the Pauli packing.
  subroutine packRhoPauliBlacsHelper(myBlacs, desc, square, kPoint, kWeight, iNeighbour,&
      & nNeighbourSK, mOrb, iCellVec, cellVec, iSparseStart, img2CentCell, imagprefac, symmetrize,&
      & primitive)

    !> BLACS matrix descriptor
    type(TBlacsEnv), intent(in) :: myBlacs

    !> Dense matrix description
    type(TDenseDescr), intent(in) :: desc

    !> distributed sparse matrix to pack
    complex(dp), intent(in) :: square(:, :)

    !> The k-point at which to unpack
    real(dp), intent(in) :: kPoint(:)

    !> weight for this k-point
    real(dp), intent(in) :: kWeight

    !> Neighbour list for the atoms (First index from 0!)
    integer, intent(in) :: iNeighbour(0:, :)

    !> Nr. of neighbours for the atoms.
    integer, intent(in) :: nNeighbourSK(:)

    !> Maximal number of orbitals on an atom.
    integer, intent(in) :: mOrb

    !> Index of the cell translation vector for each atom.
    integer, intent(in) :: iCellVec(:)

    !> Relative coordinates of the cell translation vectors.
    real(dp), intent(in) :: cellVec(:, :)

    !> indexing array for the sparse Hamiltonian
    integer, intent(in) :: iSparseStart(0:, :)

    !> Mapping between image atoms and corresponding atom in the central cell.
    integer, intent(in) :: img2CentCell(:)

    !> prefactor for the imaginary related part of a block
    complex(dp), intent(in) :: imagprefac

    !> symmetrisation of onsite blocks
    logical, intent(in) :: symmetrize

    !> sparse matrix to add this contribution to
    real(dp), intent(inout) :: primitive(:, :)

    complex(dp) :: phase
    real(dp) :: zprefac
    integer :: nAtom
    integer :: iOrig, ii, jj, kk
    integer :: iNeigh, iBlock
    integer :: iOldVec, iVec
    integer :: iAtom1, iAtom2, iAtom2f
    integer :: nOrb1, nOrb2, nOrb
    real(dp) :: kPoint2p(3)
    complex(dp), target :: tmpSqr(mOrb, mOrb)
    complex(dp), target :: tmpSqr1(mOrb, mOrb), tmpSqr2(mOrb, mOrb)
    complex(dp), pointer :: ptmp(:, :), ptmp1(:, :), ptmp2(:, :)
  #:block DEBUG_CODE
    integer :: sizePrim
  #:endblock DEBUG_CODE

    nAtom = size(iNeighbour, dim=2)
    ! number of orbitals in a regular spin block
    nOrb = desc%iAtomStart(nAtom + 1) - 1

  #:block DEBUG_CODE
    sizePrim = size(primitive, dim=1)
  #:endblock DEBUG_CODE

    @:ASSERT(nAtom > 0)
    @:ASSERT(all(shape(kPoint) == [3]))
    @:ASSERT(all(shape(nNeighbourSK) == [nAtom]))
    @:ASSERT(size(desc%iAtomStart) == nAtom + 1)
    @:ASSERT(kWeight > 0.0_dp)
    @:ASSERT(size(primitive, dim=2) == 4)

    kPoint2p(:) = 2.0_dp * pi * kPoint

    ! sigma_I and sigma_z blocks
    do iBlock = 0, 1
      zprefac = real(1 - 2 * iBlock, dp)
      iOldVec = 0
      phase = 1.0_dp
      do iAtom1 = 1, nAtom
        ii = desc%iAtomStart(iAtom1)
        nOrb1 = desc%iAtomStart(iAtom1 + 1) - ii
        do iNeigh = 0, nNeighbourSK(iAtom1)
          iOrig = iSparseStart(iNeigh, iAtom1) + 1
          iAtom2 = iNeighbour(iNeigh, iAtom1)
          iAtom2f = img2CentCell(iAtom2)
          jj = desc%iAtomStart(iAtom2f)
          nOrb2 = desc%iAtomStart(iAtom2f + 1) - jj
          iVec = iCellVec(iAtom2)
          if (iVec /= iOldVec) then
            phase = exp(cmplx(0, -1, dp) * dot_product(kPoint2p, cellVec(:, iVec)))
            iOldVec = iVec
          end if
          ptmp => tmpSqr(1:nOrb2, 1:nOrb1)
          call scalafx_cpg2l(myBlacs%orbitalGrid, desc%blacsOrbSqr, jj + iBlock * nOrb,&
              & ii + iBlock * nOrb, square, ptmp)
          ! Hermitian the on-site block before packing, as only one triangle
          ! usually supplied
          if (iAtom1 == iAtom2f) then
            do kk = 1, nOrb2
              ptmp(kk, kk+1:) = conjg(ptmp(kk+1:, kk))
            end do
          end if
          ptmp(:, :) = 0.5_dp * imagprefac * phase * kWeight * ptmp
          primitive(iOrig : iOrig + nOrb1 * nOrb2 - 1, 1) =&
              & primitive(iOrig : iOrig + nOrb1 * nOrb2 - 1, 1)&
              & + reshape(real(ptmp), [nOrb1 * nOrb2])
          primitive(iOrig : iOrig + nOrb1 * nOrb2 - 1, 4) =&
              & primitive(iOrig : iOrig + nOrb1 * nOrb2 - 1, 4)&
              & + zprefac * reshape(real(ptmp), [nOrb1 * nOrb2])
        end do
      end do
    end do

    ! sigma_x and sigma_y blocks
    iOldVec = 0
    phase = 1.0_dp
    do iAtom1 = 1, nAtom
      ii = desc%iAtomStart(iAtom1)
      nOrb1 = desc%iAtomStart(iAtom1 + 1) - ii
      do iNeigh = 0, nNeighbourSK(iAtom1)
        iOrig = iSparseStart(iNeigh, iAtom1) + 1
        iAtom2 = iNeighbour(iNeigh, iAtom1)
        iAtom2f = img2CentCell(iAtom2)
        jj = desc%iAtomStart(iAtom2f)
        nOrb2 = desc%iAtomStart(iAtom2f + 1) - jj
        iVec = iCellVec(iAtom2)
        if (iVec /= iOldVec) then
          phase = exp(cmplx(0, -1, dp) * dot_product(kPoint2p, cellVec(:, iVec)))
          iOldVec = iVec
        end if
        ptmp => tmpSqr(1:nOrb2, 1:nOrb1)
        ptmp1 => tmpSqr1(1:nOrb2, 1:nOrb1)
        ptmp2 => tmpSqr2(1:nOrb1, 1:nOrb2)
        ! take the Pauli part of the block
        call scalafx_cpg2l(myBlacs%orbitalGrid, desc%blacsOrbSqr, jj + nOrb, ii, square, ptmp1)
        call scalafx_cpg2l(myBlacs%orbitalGrid, desc%blacsOrbSqr, ii + nOrb, jj, square, ptmp2)
        ptmp(:, :) = ptmp1 + conjg(transpose(ptmp2))
        if (symmetrize .and. iAtom1 == iAtom2f) then
          do kk = 1, nOrb2
            ptmp(kk, kk+1:) = ptmp(kk+1:, kk)
          end do
        end if
        ptmp(:, :) = 0.5_dp * imagprefac * phase * kWeight * ptmp
        primitive(iOrig : iOrig + nOrb1 * nOrb2 - 1, 2) =&
            & primitive(iOrig : iOrig + nOrb1 * nOrb2 - 1, 2)&
            & + reshape(real(ptmp), [nOrb1 * nOrb2])
        ptmp(:, :) = ptmp1 - conjg(transpose(ptmp2))
        if (symmetrize .and. iAtom1 == iAtom2f) then
          do kk = 1, nOrb2
            ptmp(kk, kk+1:) = ptmp(kk+1:, kk)
          end do
        end if
        ptmp(:, :) = 0.5_dp * imagprefac * phase * kWeight * ptmp
        primitive(iOrig : iOrig + nOrb1 * nOrb2 - 1, 3) =&
            & primitive(iOrig : iOrig + nOrb1 * nOrb2 - 1, 3)&
            & + reshape(aimag(ptmp), [nOrb1 * nOrb2])
      end do
    end do

  end subroutine packRhoPauliBlacsHelper


  !> Pack only the charge (spin channel 1) part of a 2 component matrix.
  subroutine packERhoPauliBlacs(myBlacs, desc, square, kPoint, kWeight, iNeighbour, nNeighbourSK,&
      & mOrb, iCellVec, cellVec, iSparseStart, img2CentCell, primitive)

    !> BLACS matrix descriptor
    type(TBlacsEnv), intent(in) :: myBlacs

    !> Dense matrix description
    type(TDenseDescr), intent(in) :: desc

    !> local part of the distributed matrix
    complex(dp), intent(in) :: square(:, :)

    !> The k-point at which to pack
    real(dp), intent(in) :: kPoint(:)

    !> weight at this k-point
    real(dp), intent(in) :: kWeight

    !> Neighbour list for the atoms (First index from 0!)
    integer, intent(in) :: iNeighbour(0:, :)

    !> Nr. of neighbours for the atoms.
    integer, intent(in) :: nNeighbourSK(:)

    !> Maximal number of orbitals on an atom.
    integer, intent(in) :: mOrb

    !> Index of the cell translation vector for each atom.
    integer, intent(in) :: iCellVec(:)

    !> Relative coordinates of the cell translation vectors.
    real(dp), intent(in) :: cellVec(:, :)

    !> indexing array for the sparse Hamiltonian
    integer, intent(in) :: iSparseStart(0:, :)

    !> Mapping between image atoms and corresponding atom in the central cell.
    integer, intent(in) :: img2CentCell(:)

    !> sparse matrix to add this contribution to
    real(dp), intent(inout) :: primitive(:)

    complex(dp) :: phase
    integer :: nAtom, iAtom1, iAtom2, iAtom2f, iNeigh, iBlock
    integer :: iOrig, ii, jj, kk
    integer :: iOldVec, iVec
    integer :: nOrb1, nOrb2, nOrb
    real(dp) :: kPoint2p(3)
    complex(dp) :: tmpSqr(mOrb, mOrb)
  #:block DEBUG_CODE
    integer :: sizePrim
  #:endblock DEBUG_CODE

    nAtom = size(iNeighbour, dim=2)
    ! number of orbitals in a regular spin block
    nOrb = desc%iAtomStart(nAtom + 1) - 1

  #:block DEBUG_CODE
    sizePrim = size(primitive, dim=1)
  #:endblock DEBUG_CODE

    @:ASSERT(nAtom > 0)
    @:ASSERT(all(shape(kPoint) == [3]))
    @:ASSERT(all(shape(nNeighbourSK) == [nAtom]))
    @:ASSERT(size(desc%iAtomStart) == nAtom + 1)
    @:ASSERT(kWeight > 0.0_dp)

    kPoint2p(:) = 2.0_dp * pi * kPoint

    do iBlock = 0, 1
      iOldVec = 0
      phase = 1.0_dp
      do iAtom1 = 1, nAtom
        ii = desc%iAtomStart(iAtom1)
        nOrb1 = desc%iAtomStart(iAtom1 + 1) - ii
        do iNeigh = 0, nNeighbourSK(iAtom1)
          iOrig = iSparseStart(iNeigh, iAtom1) + 1
          iAtom2 = iNeighbour(iNeigh, iAtom1)
          iAtom2f = img2CentCell(iAtom2)
          jj = desc%iAtomStart(iAtom2f)
          nOrb2 = desc%iAtomStart(iAtom2f + 1) - jj
          iVec = iCellVec(iAtom2)
          if (iVec /= iOldVec) then
            phase = exp(cmplx(0, -1, dp) * dot_product(kPoint2p, cellVec(:, iVec)))
            iOldVec = iVec
          end if
          call scalafx_cpg2l(myBlacs%orbitalGrid, desc%blacsOrbSqr, jj + iBlock * nOrb,&
              & ii + iBlock * nOrb, square, tmpSqr(1:nOrb2, 1:nOrb1))
          ! Symmetrize the on-site block before packing, as only one
          ! triangle supplied
          if (iAtom1 == iAtom2f) then
            do kk = 1, nOrb2
              tmpSqr(kk, kk+1:nOrb1) = conjg(tmpSqr(kk+1:nOrb1, kk))
            end do
          end if
          primitive(iOrig : iOrig + nOrb1*nOrb2 - 1) =&
              & primitive(iOrig : iOrig + nOrb1*nOrb2 - 1) + kWeight *&
              & reshape(real(phase*tmpSqr(1:nOrb2, 1:nOrb1)), [nOrb1*nOrb2])
        end do
      end do
    end do

  end subroutine packERhoPauliBlacs


  !> Unpacks sparse H or S into dense (real, blacs) for helical geometry.
  !>
  !> Note: In contrast to the serial routines, both triangles of the resulting matrix are filled.
  !>
  subroutine unpackHSHelicalRealBlacs(myBlacs, orig, iNeighbour, nNeighbourSK, iSparseStart,&
      & img2CentCell, orb, species, coord, desc, square)

    !> BLACS matrix descriptor
    type(TBlacsEnv), intent(in) :: myBlacs

    !> sparse matrix
    real(dp), intent(in) :: orig(:)

    !> Neighbour list for each atom (First index from 0!)
    integer, intent(in) :: iNeighbour(0:, :)

    !> Nr. of neighbours for each atom (incl. itself).
    integer, intent(in) :: nNeighbourSK(:)

    !> indexing array for the sparse Hamiltonian
    integer, intent(in) :: iSparseStart(0:, :)

    !> Map from images of atoms to central cell atoms
    integer, intent(in) :: img2CentCell(:)

    !> data type for atomic orbital information
    type(TOrbitals), intent(in) :: orb

    !> Species of each atom
    integer :: species(:)

    !> Coordinates of all atoms
    real(dp), intent(in) :: coord(:,:)

    !> Dense matrix description
    type(TDenseDescr), intent(in) :: desc

    !> dense matrix part of distributed whole
    real(dp), intent(out) :: square(:, :)

    integer :: nAtom, iOrig, ii, jj, nOrb1, nOrb2, iNeigh, iAtom1, iAtom2, iAtom2f
    real(dp) :: rotZ(orb%mOrb,orb%mOrb), theta, tmpSqr(orb%mOrb,orb%mOrb)
    integer :: lShellVals(orb%mShell), iSh, iSp

    nAtom = size(iNeighbour, dim=2)
    square(:, :) = 0.0_dp
    do iAtom1 = 1, nAtom
      ii = desc%iAtomStart(iAtom1)
      nOrb1 = desc%iAtomStart(iAtom1 + 1) - ii
      do iNeigh = 0, nNeighbourSK(iAtom1)
        iOrig = iSparseStart(iNeigh, iAtom1) + 1
        iAtom2 = iNeighbour(iNeigh, iAtom1)
        iAtom2f = img2CentCell(iAtom2)
        jj = desc%iAtomStart(iAtom2f)
        nOrb2 = desc%iAtomStart(iAtom2f + 1) - jj
        iSp = species(iAtom2f)
        iSh = orb%nShell(iSp)
        lShellVals(:iSh) = orb%angShell(:iSh,iSp)
        tmpSqr(:nOrb2,:nOrb1) = reshape(orig(iOrig:iOrig+nOrb1*nOrb2-1), (/nOrb2, nOrb1/))
        theta = -atan2(coord(2,iAtom2),coord(1,iAtom2))&
            & + atan2(coord(2,iAtom2f),coord(1,iAtom2f))
        theta = mod(theta,2.0_dp*pi)
        call rotateZ(rotZ, lShellVals(:iSh), theta)
        tmpSqr(:nOrb2,:nOrb1) = matmul(rotZ(:nOrb2,:nOrb2),tmpSqr(:nOrb2,:nOrb1))
        call scalafx_addl2g(myBlacs%orbitalGrid, tmpSqr(:nOrb2,:nOrb1), desc%blacsOrbSqr,&
            & jj, ii, square)
        if (iAtom1 /= iAtom2f) then
          call scalafx_addl2g(myBlacs%orbitalGrid, transpose(tmpSqr(:nOrb2,:nOrb1)),&
              & desc%blacsOrbSqr, ii, jj, square)
        end if
      end do
    end do

  end subroutine unpackHSHelicalRealBlacs


  !> Unpacks sparse H into distributed dense matrix (complex) for helical geometry
  !>
  !> Note: In contrast to the serial routines, both triangles of the resulting matrix are filled.
  !>
  subroutine unpackHSHelicalCplxBlacs(myBlacs, orig, kPoint, iNeighbour, nNeighbourSK, iCellVec,&
      & cellVec, iSparseStart, img2CentCell, orb, species, coord, desc, square)

    !> BLACS matrix descriptor
    type(TBlacsEnv), intent(in) :: myBlacs

    !> sparse matrix
    real(dp), intent(in) :: orig(:)

    !> Relative coordinates of the K-point where the sparse matrix should be unfolded.
    real(dp), intent(in) :: kPoint(:)

    !> Neighbour list for each atom (First index from 0!)
    integer, intent(in) :: iNeighbour(0:, :)

    !> Nr. of neighbours for each atom (incl. itself).
    integer, intent(in) :: nNeighbourSK(:)

    !> Index of the cell translation vector for each atom.
    integer, intent(in) :: iCellVec(:)

    !> Relative coordinates of the cell translation vectors.
    real(dp), intent(in) :: cellVec(:, :)

    !> indexing array for the sparse Hamiltonian
    integer, intent(in) :: iSparseStart(0:, :)

    !> Mapping between image atoms and corresponding atom in the central cell.
    integer, intent(in) :: img2CentCell(:)

    !> data type for atomic orbital information
    type(TOrbitals), intent(in) :: orb

    !> Species of each atom
    integer :: species(:)

    !> Coordinates of all atoms
    real(dp), intent(in) :: coord(:,:)

    !> Dense matrix description
    type(TDenseDescr), intent(in) :: desc

    !> dense matrix part of distributed whole
    complex(dp), intent(out) :: square(:, :)

    complex(dp) :: phase, tmpSqr(orb%mOrb,orb%mOrb)
    integer :: nAtom, iOrig, nOrb1, nOrb2, ii, jj, iNeigh, iOldVec, iVec, iAtom1, iAtom2, iAtom2f
    real(dp) :: kPoint2p(2), rotZ(orb%mOrb,orb%mOrb), theta
    integer :: iSh, iSp, lShellVals(orb%mShell)

    nAtom = size(iNeighbour, dim=2)
    square(:, :) = cmplx(0, 0, dp)
    kPoint2p(:) = 2.0_dp * pi * kPoint
    iOldVec = 0
    phase = 1.0_dp
    lShellVals(:) = 0
    rotZ(:,:) = 0.0_dp
    do iAtom1 = 1, nAtom
      ii = desc%iAtomStart(iAtom1)
      nOrb1 = desc%iAtomStart(iAtom1 + 1) - ii
      do iNeigh = 0, nNeighbourSK(iAtom1)
        iOrig = iSparseStart(iNeigh, iAtom1) + 1
        iAtom2 = iNeighbour(iNeigh, iAtom1)
        iAtom2f = img2CentCell(iAtom2)
        jj = desc%iAtomStart(iAtom2f)
        nOrb2 = desc%iAtomStart(iAtom2f + 1) - jj
        iVec = iCellVec(iAtom2)
        if (iVec /= iOldVec) then
          phase = exp((0.0_dp, 1.0_dp) * dot_product(kPoint2p, cellVec(:, iVec)))
          iOldVec = iVec
        end if

        iSp = species(iAtom2f)
        iSh = orb%nShell(iSp)
        lShellVals(:iSh) = orb%angShell(:iSh,iSp)

        tmpSqr(:nOrb2,:nOrb1) = reshape(orig(iOrig:iOrig+nOrb1*nOrb2-1), (/nOrb2, nOrb1/))
        theta = -atan2(coord(2,iAtom2),coord(1,iAtom2))&
            & + atan2(coord(2,iAtom2f),coord(1,iAtom2f))
        theta = mod(theta,2.0_dp*pi)
        call rotateZ(rotZ, lShellVals(:iSh), theta)
        tmpSqr(:nOrb2,:nOrb1) = phase * matmul(rotZ(:nOrb2,:nOrb2),tmpSqr(:nOrb2,:nOrb1))

        call scalafx_addl2g(myBlacs%orbitalGrid, tmpSqr(:nOrb2,:nOrb1), desc%blacsOrbSqr, jj,&
            & ii, square)

        if (iAtom1 /= iAtom2f) then

          call scalafx_addl2g(myBlacs%orbitalGrid, transpose(conjg(tmpSqr(:nOrb2,:nOrb1))),&
              & desc%blacsOrbSqr, ii, jj, square)

        end if

      end do
    end do

  end subroutine unpackHSHelicalCplxBlacs


  !> Packs distributed dense real matrix into sparse form (real) for helical boundary conditions.
  subroutine packRhoHelicalRealBlacs(myBlacs, desc, square, iNeighbour, nNeighbourSK, iSparseStart,&
      & img2CentCell, orb, species, coord, primitive)

    !> BLACS matrix descriptor
    type(TBlacsEnv), intent(in) :: myBlacs

    !> Dense matrix description
    type(TDenseDescr), intent(in) :: desc

    !> distributed dense matrix to pack
    real(dp), intent(in) :: square(:, :)

    !> Neighbour list for the atoms (First index from 0!)
    integer, intent(in) :: iNeighbour(0:, :)

    !> Nr. of neighbours for the atoms.
    integer, intent(in) :: nNeighbourSK(:)

    !> indexing array for the sparse Hamiltonian
    integer, intent(in) :: iSparseStart(0:, :)

    !> Mapping between image atoms and corresponding atom in the central cell.
    integer, intent(in) :: img2CentCell(:)

    !> data type for atomic orbital information
    type(TOrbitals), intent(in) :: orb

    !> Species of each atom
    integer :: species(:)

    !> Coordinates of all atoms
    real(dp), intent(in) :: coord(:,:)

    !> sparse matrix to add this contribution into
    real(dp), intent(inout) :: primitive(:)

    integer :: nAtom, iOrig, ii, jj, kk, iNeigh, iAtom1, iAtom2, iAtom2f, nOrb1, nOrb2
    integer :: iSp, iSh, lShellVals(orb%mShell)
    real(dp) :: tmpSqr(orb%mOrb, orb%mOrb), rotZ(orb%mOrb,orb%mOrb), theta

    nAtom = size(iNeighbour, dim=2)
    lShellVals(:) = 0
    do iAtom1 = 1, nAtom
      ii = desc%iAtomStart(iAtom1)
      nOrb1 = desc%iAtomStart(iAtom1 + 1) - ii
      do iNeigh = 0, nNeighbourSK(iAtom1)
        iOrig = iSparseStart(iNeigh, iAtom1) + 1
        iAtom2 = iNeighbour(iNeigh, iAtom1)
        iAtom2f = img2CentCell(iAtom2)
        jj = desc%iAtomStart(iAtom2f)
        nOrb2 = desc%iAtomStart(iAtom2f + 1) - jj
        call scalafx_cpg2l(myBlacs%orbitalGrid, desc%blacsOrbSqr, jj, ii, square,&
            & tmpSqr(1:nOrb2, 1:nOrb1))
        iSp = species(iAtom2f)
        iSh = orb%nShell(iSp)
        lShellVals(:iSh) = orb%angShell(:iSh,iSp)
        theta = atan2(coord(2,iAtom2),coord(1,iAtom2)) &
            & - atan2(coord(2,iAtom2f),coord(1,iAtom2f))
        theta = mod(theta,2.0_dp*pi)
        call rotateZ(rotZ, lShellVals(:iSh), theta)
        ! Symmetrize the on-site block before packing, just in case
        if (iAtom1 == iAtom2f) then
          do kk = 1, nOrb2
            tmpSqr(kk, kk+1:nOrb1) = tmpSqr(kk + 1 : nOrb1, kk)
          end do
        end if
        tmpSqr(:nOrb2,:nOrb1) =  matmul(rotZ(:nOrb2,:nOrb2),tmpSqr(:nOrb2,:nOrb1))
        primitive(iOrig : iOrig + nOrb1 * nOrb2 - 1) =&
            & primitive(iOrig : iOrig + nOrb1 * nOrb2 - 1)&
            & + reshape(tmpSqr(1:nOrb2, 1:nOrb1), [nOrb1 * nOrb2])
      end do
    end do

  end subroutine packRhoHelicalRealBlacs


  !> Packs distributed dense real matrix into sparse form (real).
  subroutine packRhoHelicalCplxBlacs(myblacs, desc, square, kPoint, kWeight, iNeighbour,&
      & nNeighbourSK, iCellVec, cellVec, iSparseStart, img2CentCell, orb, species, coord,&
      & primitive)

    !> BLACS matrix descriptor
    type(TBlacsEnv), intent(in) :: myBlacs

    !> Dense matrix description
    type(TDenseDescr), intent(in) :: desc

    !> Distributed dense matrix to pack
    complex(dp), intent(in) :: square(:, :)

    !> The k-point at which to pack
    real(dp), intent(in) :: kPoint(:)

    !> weight for this k-point
    real(dp), intent(in) :: kWeight

    !> Neighbour list for the atoms (First index from 0!)
    integer, intent(in) :: iNeighbour(0:, :)

    !> Nr. of neighbours for the atoms.
    integer, intent(in) :: nNeighbourSK(:)

    !> Index of the cell translation vector for each atom.
    integer, intent(in) :: iCellVec(:)

    !> Relative coordinates of the cell translation vectors.
    real(dp), intent(in) :: cellVec(:, :)

    !> indexing array for the sparse Hamiltonian
    integer, intent(in) :: iSparseStart(0:, :)

    !> Mapping between image atoms and corresponding atom in the central cell.
    integer, intent(in) :: img2CentCell(:)

    !> data type for atomic orbital information
    type(TOrbitals), intent(in) :: orb

    !> Species of each atom
    integer :: species(:)

    !> Coordinates of all atoms
    real(dp), intent(in) :: coord(:,:)

    !> sparse matrix to add this contribution into
    real(dp), intent(inout) :: primitive(:)

    complex(dp) :: phase
    integer :: nAtom, iOrig, ii, jj, kk, iNeigh, iOldVec, iVec, iAtom1, iAtom2, iAtom2f
    integer :: nOrb1, nOrb2, iSp, iSh, lShellVals(orb%mShell)
    real(dp) :: kPoint2p(2), rotZ(orb%mOrb,orb%mOrb), theta
    complex(dp) :: tmpSqr(orb%mOrb, orb%mOrb)

    nAtom = size(iNeighbour, dim=2)
    kPoint2p(:) = 2.0_dp * pi * kPoint
    lShellVals(:) = 0
    iOldVec = 0
    phase = 1.0_dp
    do iAtom1 = 1, nAtom
      ii = desc%iAtomStart(iAtom1)
      nOrb1 = desc%iAtomStart(iAtom1 + 1) - ii
      do iNeigh = 0, nNeighbourSK(iAtom1)
        iOrig = iSparseStart(iNeigh, iAtom1) + 1
        iAtom2 = iNeighbour(iNeigh, iAtom1)
        iAtom2f = img2CentCell(iAtom2)
        jj = desc%iAtomStart(iAtom2f)
        nOrb2 = desc%iAtomStart(iAtom2f + 1) - jj
        iVec = iCellVec(iAtom2)
        if (iVec /= iOldVec) then
          phase = exp(cmplx(0, -1, dp) * dot_product(kPoint2p(:), cellVec(:, iVec)))
          iOldVec = iVec
        end if
        call scalafx_cpg2l(myBlacs%orbitalGrid, desc%blacsOrbSqr, jj, ii, square,&
            & tmpSqr(1:nOrb2, 1:nOrb1))
        iSp = species(iAtom2f)
        iSh = orb%nShell(iSp)
        lShellVals(:iSh) = orb%angShell(:iSh,iSp)
        theta = atan2(coord(2,iAtom2),coord(1,iAtom2)) &
            & - atan2(coord(2,iAtom2f),coord(1,iAtom2f))
        theta = mod(theta,2.0_dp*pi)
        call rotateZ(rotZ, lShellVals(:iSh), theta)
        ! Hermitian the on-site block before packing, just in case
        if (iAtom1 == iAtom2f) then
          do kk = 1, nOrb2
            tmpSqr(kk, kk + 1 : nOrb1) = conjg(tmpSqr(kk + 1 : nOrb1, kk))
          end do
        end if
        tmpSqr(:nOrb2,:nOrb1) = matmul(rotZ(:nOrb2,:nOrb2),tmpSqr(:nOrb2,:nOrb1))
        primitive(iOrig : iOrig + nOrb1*nOrb2 - 1) =&
            & primitive(iOrig : iOrig + nOrb1*nOrb2 - 1)&
            & + kWeight * real(phase * reshape(tmpSqr(1:nOrb2, 1:nOrb1), [nOrb1 * nOrb2]), dp)
      end do
    end do

  end subroutine packRhoHelicalCplxBlacs

#:endif


  !> Calculate indexing array and number of elements in sparse arrays like the real space overlap
  subroutine getSparseDescriptor(iNeighbour, nNeighbourSK, img2CentCell, orb, iPair, sparseSize)

    !> Neighbours of each atom
    integer, intent(in) :: iNeighbour(0:,:)

    !> Number of neighbours of each atom
    integer, intent(in) :: nNeighbourSK(:)

    !> Indexing for mapping image atoms to central cell
    integer, intent(in) :: img2CentCell(:)

    !> Atomic orbital information
    type(TOrbitals), intent(in) :: orb

    !> Sparse array indexing for the start of atomic blocks in data structures
    integer, allocatable, intent(inout) :: iPair(:,:)

    !> Total number of elements in a sparse structure (ignoring extra indices like spin)
    integer, intent(out) :: sparseSize

    integer :: nAtom, mNeighbour
    integer :: ind, iAt1, nOrb1, iNeigh1, nOrb2

    nAtom = size(iNeighbour, dim=2)
    mNeighbour = size(iNeighbour, dim=1)

    @:ASSERT(allocated(iPair))
    @:ASSERT(size(iPair, dim=2) == nAtom)

    if (mNeighbour > size(iPair, dim=1)) then
      deallocate(iPair)
      allocate(iPair(0 : mNeighbour - 1, nAtom))
      iPair(:,:) = 0
    end if
    ind = 0
    do iAt1 = 1, nAtom
      nOrb1 = orb%nOrbAtom(iAt1)
      do iNeigh1 = 0, nNeighbourSK(iAt1)
        iPair(iNeigh1, iAt1) = ind
        nOrb2 = orb%nOrbAtom(img2CentCell(iNeighbour(iNeigh1, iAt1)))
        ind = ind + nOrb1 * nOrb2
      end do
    end do
    sparseSize = ind

  end subroutine getSparseDescriptor

end module dftbp_dftb_sparse2dense
