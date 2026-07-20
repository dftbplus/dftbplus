!--------------------------------------------------------------------------------------------------!
!  DFTB+: general package for performing fast atomistic simulations                                !
!  Copyright (C) 2006 - 2024  DFTB+ developers group                                               !
!                                                                                                  !
!  See the LICENSE file for terms of usage and distribution.                                       !
!--------------------------------------------------------------------------------------------------!

#:include 'common.fypp'
#:include 'error.fypp'

!> Module for evaluating momentum matrix elements between filled and empty states, see doi:
!! 10.1103/PhysRevB.63.201101 10.1103/PhysRevB.98.115115
module dftbp_derivs_momentumelements
  use dftbp_common_accuracy, only : dp
  use dftbp_common_environment, only : TEnvironment
  use dftbp_common_status, only : TStatus
  use dftbp_dftb_densitymatrix, only : TDensityMatrix
  use dftbp_dftb_hybridxc, only : THybridXcFunc
  use dftbp_dftb_periodic, only : TAuxNeighbourList, TNeighbourList
  use dftbp_dftb_sparse2dense, only : unpackHSdk
  use dftbp_type_commontypes, only : TOrbitals, TParallelKS
  use dftbp_type_integral, only : TIntegral
  use dftbp_type_densedescr, only : TDenseDescr
  use dftbp_type_wrappedintr, only : TWrappedInt2, TWrappedReal1
  use dftbp_math_blasroutines, only : hemm

  implicit none

  private
  public :: momentumElements

contains


  !> Calculate the momentum elements between occupied and virtual states
  subroutine momentumElements(env, pElement, parallelKS, eigvals, eigVecsCplx, ints, nTrans, getIA,&
      & wij, neighbourList, nNeighbourSK, symNeighbourList, nNeighbourCamSym, orb, denseDesc,&
      & iSparseStart, img2CentCell, kPoint, kWeight, rCellVecs, cellVec, iCellVec, densityMatrix,&
      & hybridXc, dab, errStatus)

    !> Environment settings
    type(TEnvironment), intent(inout) :: env

    complex(dp), allocatable, intent(out) :: pElement(:,:,:)

    !> K-points and spins to process
    type(TParallelKS), intent(in) :: parallelKS

    !> Eigenvalue of each level, kpoint and spin channel
    real(dp), intent(in) :: eigvals(:,:,:)

    !> ground state complex eigenvectors
    complex(dp), intent(in), allocatable :: eigvecsCplx(:,:,:)

    !> Integral container
    type(TIntegral), intent(in) :: ints

    integer, intent(in) :: nTrans(:)

    type(TWrappedInt2), intent(in) :: getIA(:)

    type(TWrappedReal1), intent(in) :: wij(:)

    !> list of neighbours for each atom
    type(TNeighbourList), intent(in) :: neighbourList

    !> Number of neighbours for each of the atoms
    integer, intent(in) :: nNeighbourSK(:)

    !> List of neighbouring atoms (symmetric version)
    type(TAuxNeighbourList), intent(in), allocatable :: symNeighbourList

    !> Symmetric neighbour list version of nNeighbourCam
    integer, intent(in), allocatable :: nNeighbourCamSym(:)

    !> Atomic orbital information
    type(TOrbitals), intent(in) :: orb

    !> Dense matrix descriptor
    type(TDenseDescr), intent(in) :: denseDesc

    !> Index array for the start of atomic blocks in sparse arrays
    integer, intent(in) :: iSparseStart(:,:)

    !> map from image atoms to the original unique atom
    integer, intent(in) :: img2CentCell(:)

    !> k-points
    real(dp), intent(in) :: kPoint(:,:)

    !> Weights for k-points
    real(dp), intent(in) :: kWeight(:)

    !> Vectors to units cells in absolute units
    real(dp), intent(in) :: rCellVecs(:,:)

    !> Vectors (in units of the lattice constants) to cells of the lattice
    real(dp), intent(in) :: cellVec(:,:)

    !> Index for which unit cell atoms are associated with
    integer, intent(in) :: iCellVec(:)

    !> Holds real and complex delta density matrices and pointers
    type(TDensityMatrix), intent(in) :: densityMatrix

    !> Hybrid functional, if in use
    class(THybridXcFunc), allocatable, intent(inout) :: hybridXc

    !> On-site atomic dipole matrix elements, if available
    real(dp), allocatable, intent(in) :: dab(:,:,:,:)

    !> Status of operation
    type(TStatus), intent(out) :: errStatus

    integer :: iAt, iK, iKS, iS, iCart, nOrbs, iTrans, nLocCol, nLocRow, iFil, iEmp, iOrb, ii, jj
    complex(dp), allocatable :: work(:,:), work2(:,:), dHamHyb_dk(:,:,:,:)

  #:if WITH_SCALAPACK

    @:ASSERT(.false.)

  #:endif

    allocate(pElement(3, maxval(nTrans), parallelKS%nLocalKS), source=(0.0_dp,0.0_dp))

    nOrbs = size(eigVecsCplx, dim=1)
    nLocCol = size(eigVecsCplx, dim=1)
    nLocRow = size(eigVecsCplx, dim=2)

    allocate(work(nLocCol, nLocRow))
    allocate(work2(nLocCol, nLocRow))

    if (allocated(hybridXc)) then

      call hybridXc%getCamHamiltonian_dkpts(env, denseDesc, orb, ints, densityMatrix,&
          & neighbourList, nNeighbourSK, symNeighbourList, nNeighbourCamSym, iCellVec, cellVec,&
          & rCellVecs, iSparseStart, img2CentCell, kPoint, kWeight, dHamHyb_dk, errStatus)
      @:PROPAGATE_ERROR(errStatus)
    end if

    do iKS = 1, parallelKS%nLocalKS

      iK = parallelKS%localKS(1, iKS)
      iS = parallelKS%localKS(2, iKS)

      do iCart = 1, 3

        ! First term of (7) from 10.1103/PhysRevB.98.115115 since basis is non-orthogonal
        ! Pulay term <psi|dS/dk|psi>
        work(:,:) = cmplx(0,0,dp)
        call unpackHSdk(work, ints%overlap, kPoint(:,iK), neighbourList%iNeighbour, nNeighbourSK,&
            & iCellVec, cellVec, denseDesc%iAtomStart, iSparseStart, img2CentCell, iCart)
        call hemm(work2, 'l', work, eigVecsCplx(:,:,iKS))
        do iOrb = 1, nOrbs
          work2(:,iOrb) = -eigvals(iOrb, iK, iS) * work2(:,iOrb)
        end do

        ! <psi|dH/dk|psi>
        work(:,:) = cmplx(0,0,dp)
        call unpackHSdk(work, ints%hamiltonian(:,iS), kPoint(:,iK), neighbourList%iNeighbour,&
            & nNeighbourSK, iCellVec, cellVec, denseDesc%iAtomStart, iSparseStart, img2CentCell,&
            & iCart)
        if (allocated(hybridXc)) then
          work(:,:) = work + dHamHyb_dk(:, :, iKS, iCart)
        end if
        call hemm(work2, 'l', work, eigVecsCplx(:,:,iKS), beta=cmplx(1,0,dp))

        do iTrans = 1, nTrans(iKS)
          iFil = getIA(iKS)%data(iTrans, 1)
          iEmp = getIA(iKS)%data(iTrans, 2)
          pElement(iCart, iTrans, iKS) = dot_product(eigVecsCplx(:,iEmp,iKS), work2(:,iFil))
        end do

        if (allocated(dab)) then
          ! On-site dipole elements

          work(:,:) = cmplx(0,0,dp)
          do iAt = 1, size(nNeighbourSK)
            ii = denseDesc%iAtomStart(iAt)
            jj = denseDesc%iAtomStart(iAt+1)-1
            work(ii:jj,ii:jj) = dab(orb%nOrbAtom(iAt), orb%nOrbAtom(iAt), iAt, iCart)
          end do

          call hemm(work2, 'l', work, eigVecsCplx(:,:,iKS))

          do iTrans = 1, nTrans(iKS)
            iFil = getIA(iKS)%data(iTrans, 1)
            iEmp = getIA(iKS)%data(iTrans, 2)
            pElement(iCart, iTrans, iKS) = pElement(iCart, iTrans, iKS)&
                & + wij(iKS)%data(nTrans(iKS)) * dot_product(eigVecsCplx(:,iEmp,iKS),work2(:,iFil))
          end do

        end if

      end do

    end do

  end subroutine momentumElements


end module dftbp_derivs_momentumelements
