!--------------------------------------------------------------------------------------------------!
!  DFTB+: general package for performing fast atomistic simulations                                !
!  Copyright (C) 2018  DFTB+ developers group                                                      !
!                                                                                                  !
!  See the LICENSE file for terms of usage and distribution.                                       !
!--------------------------------------------------------------------------------------------------!

#:include 'common.fypp'

!> Module for linear response derivative calculations using perturbation methods
module dftbp_perturbkderivs
  use dftbp_accuracy
  use dftbp_constants
  use dftbp_globalenv
  use dftbp_message
  use dftbp_commontypes
  use dftbp_potentials
  use dftbp_scc
  use dftbp_orbitalequiv
  use dftbp_populations
  use dftbp_spin
  use dftbp_thirdorder_module, only : ThirdOrder
  use dftbp_dftbplusu
  use dftbp_onsitecorrection
  use dftbp_mainio
  use dftbp_shift
  use dftbp_mixer
  use dftbp_finitethelper
  use dftbp_scalapackfx
  use dftbp_environment
  use dftbp_periodic
  use dftbp_densedescr
  use dftbp_sparse2dense
  use dftbp_degeneratePerturb
  use dftbp_taggedoutput
#:if WITH_MPI
  use dftbp_mpifx
#:endif
#:if WITH_SCALAPACK
  use dftbp_scalafxext
#:else
  use dftbp_blasroutines
#:endif

  implicit none

  private
  public :: dPsidK

contains

  subroutine dPsidK(env, parallelKS, eigvals, eigVecsCplx, ham, over, orb, nAtom, species,&
      & neighbourList, nNeighbourSK, denseDesc, iSparseStart, img2CentCell, coord, kPoint, kWeight,&
      & cellVec, iCellVec, latVec)

    !> Environment settings
    type(TEnvironment), intent(inout) :: env

    !> K-points and spins to process
    type(TParallelKS), intent(in) :: parallelKS

    !> Eigenvalue of each level, kpoint and spin channel
    real(dp), intent(in) :: eigvals(:,:,:)

    !> ground state complex eigenvectors
    complex(dp), intent(in), allocatable :: eigvecsCplx(:,:,:)

    !> Sparse Hamiltonian
    real(dp), intent(in) :: ham(:,:)

    !> sparse overlap matrix
    real(dp), intent(in) :: over(:)

    !> Atomic orbital information
    type(TOrbitals), intent(in) :: orb

    !> Number of central cell atoms
    integer, intent(in) :: nAtom

    !> chemical species
    integer, intent(in) :: species(:)

    !> list of neighbours for each atom
    type(TNeighbourList), intent(in) :: neighbourList

    !> Number of neighbours for each of the atoms
    integer, intent(in) :: nNeighbourSK(:)

    !> Dense matrix descriptor
    type(TDenseDescr), intent(in) :: denseDesc

    !> Index array for the start of atomic blocks in sparse arrays
    integer, intent(in) :: iSparseStart(:,:)

    !> map from image atoms to the original unique atom
    integer, intent(in) :: img2CentCell(:)

    !> atomic coordinates
    real(dp), intent(in) :: coord(:,:)

    !> k-points
    real(dp), intent(in) :: kPoint(:,:)

    !> Weights for k-points
    real(dp), intent(in) :: kWeight(:)

    !> Vectors (in units of the lattice constants) to cells of the lattice
    real(dp), intent(in) :: cellVec(:,:)

    !> Index for which unit cell atoms are associated with
    integer, intent(in) :: iCellVec(:)

    real(dp), intent(in) :: latVec(:,:)

    integer :: nIndepHam, nSpin, nOrbs, nKpts, iAt1, iNeigh, iAt2, iAt2f, ii, jj, kk, iCell, iSpin
    integer :: iDir, iK, iOrb, iKS, iS
    real(dp), allocatable :: dHam(:,:), dOver(:), dEi(:,:)
    real(dp) :: vec(3)
    complex(dp), allocatable :: dPsi(:,:,:,:), dHamSqr(:,:,:), dOverSqr(:,:,:), workLocal(:,:)
    complex(dp), allocatable :: hamSqr(:,:,:), overSqr(:,:,:), eCi(:,:)

  #:if WITH_SCALAPACK
    ! need distributed matrix descriptors
    integer :: desc(DLEN_), nn

    nn = denseDesc%fullSize
    call scalafx_getdescriptor(env%blacs%orbitalGrid, nn, nn, env%blacs%rowBlockSize,&
        & env%blacs%columnBlockSize, desc)
  #:endif

    if (.not.allocated(eigVecsCplx)) then
      call error("Missing complex eigenvectors")
    end if

  #:if WITH_SCALAPACK
  #:else

    nSpin = size(ham, dim=2)
    nOrbs = size(eigVecsCplx,dim=1)
    nKpts = size(kpoint, dim=2)

    allocate(dHamSqr(nOrbs, nOrbs, 3))
    allocate(dOverSqr(nOrbs, nOrbs, 3))
    allocate(eCi(nOrbs, nOrbs))
    allocate(dEi(nOrbs, nOrbs))


    do iKS = 1, parallelKS%nLocalKS
      iK = parallelKS%localKS(1, iKS)
      iS = parallelKS%localKS(2, iKS)
      do iOrb = 1, nOrbs
        eCi(:,iOrb) = eigvals(iOrb,iK,iS) * eigvecsCplx(:,iOrb,iKS)
      end do

      dOverSqr(:,:,:) = cmplx(0,0,dp)
      dHamSqr(:,:,:) = cmplx(0,0,dp)
      call unpackHSdk(dOverSqr, over, kPoint(:,iK), neighbourList%iNeighbour, nNeighbourSK,&
          & iCellVec, cellVec, denseDesc%iAtomStart, iSparseStart, img2CentCell)
      call unpackHSdk(dHamSqr, ham(:,iS), kPoint(:,iK), neighbourList%iNeighbour, nNeighbourSK,&
          & iCellVec, cellVec, denseDesc%iAtomStart, iSparseStart, img2CentCell)
      do kk = 1, 3
        do iOrb = 1, nOrbs
          dOverSqr(iOrb, iOrb+1:,kk) = conjg(dOverSqr(iOrb+1:, iOrb, kk))
          dHamSqr(iOrb, iOrb+1:,kk) = conjg(dHamSqr(iOrb+1:, iOrb, kk))
          dOverSqr(iOrb, iOrb, kk) = real(dOverSqr(iOrb, iOrb, kk))
          dHamSqr(iOrb, iOrb, kk) = real(dHamSqr(iOrb, iOrb, kk))
        end do
      end do

      do kk = 1, 3
        dEi(:,:) = real(matmul(transpose(conjg(eigvecsCplx(:,:,iKS))),&
            & matmul(dHamSqr(:,:,kk), eigvecsCplx(:,:,iKS))&
            & - matmul(dOverSqr(:,:,kk), eCi)))
        write(*,*)iK,kk,(dEi(iOrb, iOrb)*Hartree__eV, iOrb = 1, nOrbs)
      end do

    end do

  #:endif

  end subroutine dPsidK

end module dftbp_perturbkderivs
