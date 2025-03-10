!--------------------------------------------------------------------------------------------------!
!  DFTB+: general package for performing fast atomistic simulations                                !
!  Copyright (C) 2006 - 2025  DFTB+ developers group                                               !
!                                                                                                  !
!  See the LICENSE file for terms of usage and distribution.                                       !
!--------------------------------------------------------------------------------------------------!

#:include 'common.fypp'
#:include 'error.fypp'

!> Calculate the whole single particle density matrix and energy weighted density matrix.
!! Note: Dense serial code based on suggestions from Thomas Heine.
!! Caveat: The routines create the transposed and complex conjugated of the density matrices
!! (cc* instead of the conventional c*c).
module dftbp_dftb_densitymatrix
  use dftbp_common_accuracy, only : dp, lc
  use dftbp_common_constants, only : pi, imag
  use dftbp_common_status, only : TStatus
  use dftbp_elecsolvers_dmsolvertypes, only : densityMatrixTypes
  use dftbp_math_blasroutines, only : herk
  use dftbp_math_sorting, only : unique, heap_sort
  use dftbp_type_commontypes, only : TOrbitals, TParallelKS
  use dftbp_common_globalenv, only : stdOut
#:if WITH_SCALAPACK
  use dftbp_extlibs_scalapackfx, only : blacsgrid, blocklist, size, pblasfx_pgemm, pblasfx_ptranc
#:endif
#:if WITH_MAGMA
  use magma, only : magmaf_queue_create, magmaf_dsetmatrix, magmaf_dsyrk, magmaf_dgetmatrix,&
      & magmaf_dmalloc, magmaf_free, magmaf_zsetmatrix, magmaf_zherk, magmaf_zgetmatrix,&
      & magmaf_zmalloc, magmaf_queue_destroy
  use iso_fortran_env, only : int64
#:endif

  implicit none
  private

  public :: TDensityMatrix, TDensityMatrix_init, transformDualSpaceToBvKRealSpace
#:if WITH_SCALAPACK
  public :: makeDensityMtxRealBlacs, makeDensityMtxCplxBlacs
#:endif


  !> Holds density matrix related data and control variables
  type :: TDensityMatrix

    !> Method by which the density matrix is being generated
    integer :: iDensityMatrixAlgorithm

    !! Real and complex delta (i.e. change between cycles) density matrices.

    !> Real, square dual-space deltaRho input
    real(dp), allocatable :: deltaRhoIn(:,:,:)

    !> Real, square dual-space deltaRho output
    real(dp), allocatable :: deltaRhoOut(:,:,:)

    !> Complex, square dual-space deltaRho input
    complex(dp), allocatable :: deltaRhoInCplx(:,:,:)

    !> Complex, square dual-space deltaRho output
    complex(dp), allocatable :: deltaRhoOutCplx(:,:,:)

    !> Real-space, square deltaRho input
    real(dp), allocatable :: deltaRhoInCplxHS(:,:,:,:,:,:)

    !> Real-space, square deltaRho output
    real(dp), allocatable :: deltaRhoOutCplxHS(:,:,:,:,:,:)

    !> Composite index for mapping iK/iS --> iGlobalKS for arrays present at every MPI rank
    integer, allocatable :: iKiSToiGlobalKS(:,:)

    !> The k'-points that are possibly different from the current k-points in case of a
    !! bandstructure calculation
    real(dp), allocatable :: kPointPrime(:,:)

    !> Weights of the k'-points
    real(dp), allocatable :: kWeightPrime(:)

  contains

    !> Returns the density matrix
    procedure, private :: getDensityMatrix_real
    procedure, private :: getDensityMatrix_cmplx
    generic :: getDensityMatrix => getDensityMatrix_real, getDensityMatrix_cmplx

    !> Returns the energy weighted density matrix
    procedure, private :: getEDensityMatrix_real
    procedure, private :: getEDensityMatrix_cmplx
    generic :: getEDensityMatrix => getEDensityMatrix_real, getEDensityMatrix_cmplx

  end type TDensityMatrix


  !> Used to shift occupations in density matrix build away from 0
  real(dp), parameter :: arbitraryConstant = 0.1_dp

contains


  !> Initialise the density matrix container
  subroutine TDensityMatrix_init(this, iDensityMatrixAlgorithm)

    !> Instance
    type(TDensityMatrix), intent(out) :: this

    !> Density matrix generation method
    integer, intent(in) :: iDensityMatrixAlgorithm

    this%iDensityMatrixAlgorithm = iDensityMatrixAlgorithm

  end subroutine TDensityMatrix_init


  !> Returns the density matrix
  subroutine getDensityMatrix_real(this, densityMatrix, eigenvecs, filling, errStatus)

    !> Instance
    class(TDensityMatrix), intent(in) :: this

    !> Resulting density matrix
    real(dp), intent(inout) :: densityMatrix(:,:)

    !> Eigenvectors of the system
    real(dp), intent(inout) :: eigenvecs(:,:)

    !> Occupation numbers of the orbitals
    real(dp), intent(in) :: filling(:)

    !> Error status
    type(TStatus), intent(out) :: errStatus

    select case(this%iDensityMatrixAlgorithm)

    case (densityMatrixTypes%magma_fromEigenVecs)

    #:if WITH_MAGMA

      call makeDensityMtxRealGPU(densityMatrix, eigenvecs, filling, errStatus)
      @:PROPAGATE_ERROR(errStatus)

    #:else

      @:RAISE_ERROR(errStatus, -1, "Internal error, compiled without MAGMA, but this is requested&
          & for density matrix generation method")

    #:endif

    case (densityMatrixTypes%fromEigenVecs)

      call fullDensityMatrix_real(densityMatrix, eigenvecs, filling)

    case (densityMatrixTypes%elecSolverProvided)

      return

    case (densityMatrixTypes%none)

      densityMatrix(:,:) = 0.0_dp

    case default

      @:RAISE_ERROR(errStatus, -1, "Internal error, unknown density matrix generation method")

    end select

  end subroutine getDensityMatrix_real


  !> Returns the energy weighted density matrix
  subroutine getEDensityMatrix_real(this, egyDensityMatrix, eigenvecs, filling, eigenvals,&
      & errStatus)

    !> Instance
    class(TDensityMatrix), intent(in) :: this

    !> Resulting density matrix
    real(dp), intent(inout) :: egyDensityMatrix(:,:)

    !> Eigenvectors of the system
    real(dp), intent(inout) :: eigenvecs(:,:)

    !> Occupation numbers of the orbitals
    real(dp), intent(in) :: filling(:)

    !> Eigenvalues of the system
    real(dp), intent(in) :: eigenvals(:)

    !> Error status
    type(TStatus), intent(out) :: errStatus

    select case(this%iDensityMatrixAlgorithm)

    case (densityMatrixTypes%magma_fromEigenVecs)

    #:if WITH_MAGMA

      call makeEnergyDensityMtxRealGPU(egyDensityMatrix, eigenvecs, filling, eigenvals, errStatus)
      @:PROPAGATE_ERROR(errStatus)

    #:else

      @:RAISE_ERROR(errStatus, -1, "Internal error, compiled without MAGMA, but this is requested&
          & for energy weighted density matrix generation method")

    #:endif

    case (densityMatrixTypes%fromEigenVecs)

      call fullEnergyDensityMatrix_real(egyDensityMatrix, eigenvecs, filling, eigenvals)

    case (densityMatrixTypes%elecSolverProvided)

      return

    case (densityMatrixTypes%none)

      egyDensityMatrix(:,:) = 0.0_dp

    case default

      @:RAISE_ERROR(errStatus, -1, "Internal error, unknown density matrix generation method")

    end select

  end subroutine getEDensityMatrix_real


  !> Returns the density matrix
  subroutine getDensityMatrix_cmplx(this, densityMatrix, eigenvecs, filling, errStatus)

    !> Instance
    class(TDensityMatrix), intent(in) :: this

    !> Resulting density matrix
    complex(dp), intent(inout) :: densityMatrix(:,:)

    !> Eigenvectors of the system
    complex(dp), intent(inout) :: eigenvecs(:,:)

    !> Occupation numbers of the orbitals
    real(dp), intent(in) :: filling(:)

    !> Error status
    type(TStatus), intent(out) :: errStatus

    select case(this%iDensityMatrixAlgorithm)

    case (densityMatrixTypes%magma_fromEigenVecs)

    #:if WITH_MAGMA

      call makeDensityMtxCmplxGPU(densityMatrix, eigenvecs, filling, errStatus)
      @:PROPAGATE_ERROR(errStatus)

    #:else

      @:RAISE_ERROR(errStatus, -1, "Internal error, compiled without MAGMA, but this is requested&
          & for density matrix generation method")

    #:endif

    case (densityMatrixTypes%fromEigenVecs)

      call fullDensityMatrix_cmplx(densityMatrix, eigenvecs, filling)

    case (densityMatrixTypes%elecSolverProvided)

      return

    case (densityMatrixTypes%none)

      densityMatrix(:,:) = 0.0_dp

    case default

      @:RAISE_ERROR(errStatus, -1, "Internal error, unknown density matrix generation method")

    end select

  end subroutine getDensityMatrix_cmplx


  !> Returns the energy-weighted density matrix
  subroutine getEDensityMatrix_cmplx(this, egyDensityMatrix, eigenvecs, filling, eigenvals,&
      & errStatus)

    !> Instance
    class(TDensityMatrix), intent(in) :: this

    !> Resulting density matrix
    complex(dp), intent(inout) :: egyDensityMatrix(:,:)

    !> Eigenvectors of the system
    complex(dp), intent(inout) :: eigenvecs(:,:)

    !> Occupation numbers of the orbitals
    real(dp), intent(in) :: filling(:)

    !> Eigenvalues of the system
    real(dp), intent(in) :: eigenvals(:)

    !> Error status
    type(TStatus), intent(out) :: errStatus

    select case(this%iDensityMatrixAlgorithm)

    case (densityMatrixTypes%magma_fromEigenVecs)

    #:if WITH_MAGMA

      call makeEnergyDensityMtxCmplxGPU(egyDensityMatrix, eigenvecs, filling, eigenvals, errStatus)
      @:PROPAGATE_ERROR(errStatus)

    #:else

      @:RAISE_ERROR(errStatus, -1, "Internal error, compiled without MAGMA, but this is requested&
          & for density matrix generation method")

    #:endif

    case (densityMatrixTypes%fromEigenVecs)

      call fullEnergyDensityMatrix_cmplx(egyDensityMatrix, eigenvecs, filling, eigenvals)

    case (densityMatrixTypes%elecSolverProvided)

      return

    case (densityMatrixTypes%none)

      egyDensityMatrix(:,:) = cmplx(0,0,dp)

    case default

      @:RAISE_ERROR(errStatus, -1, "Internal error, unknown density matrix generation method")

    end select

  end subroutine getEDensityMatrix_cmplx


  !> Transforms dense, square density matrix for all spins/k-points to real-space (BvK cell).
  subroutine transformDualSpaceToBvKRealSpace(rhoDual, parallelKS, kPoint, kWeight, bvKShifts,&
      & rhoBvK)

    !> Complex, dense, square dual-space rho of all spins/k-points
    complex(dp), intent(in) :: rhoDual(:,:,:)

    !> The k-points and spins to process
    type(TParallelKS), intent(in) :: parallelKS

    !> The k-points in relative units
    real(dp), intent(in) :: kPoint(:,:)

    !> Weights of k-points
    real(dp), intent(in) :: kWeight(:)

    !> The k-point compatible BvK real-space shifts in relative coordinates (units of latVecs)
    real(dp), intent(in) :: bvKShifts(:,:)

    !> Real-space, dense, square rho for BvK cell
    real(dp), intent(inout) :: rhoBvK(:,:,:,:,:,:)

    !! K-point-spin composite index and k-point/spin index
    integer :: iKS, iK, iSpin

    !! Iterates over all BvK real-space vectors
    integer :: iG

    !! Phase factor
    complex(dp) :: phase

    !! Integer BvK real-space shift translated to density matrix indices
    integer :: bvKIndex(3)

    do iG = 1, size(bvKShifts, dim=2)
      bvKIndex(:) = nint(bvKShifts(:, iG)) + 1
      do iKS = 1, parallelKS%nLocalKS
        iK = parallelKS%localKS(1, iKS)
        iSpin = parallelKS%localKS(2, iKS)
        phase = exp(-imag * dot_product(2.0_dp * pi * kPoint(:, iK), bvKShifts(:, iG)))
        rhoBvK(:,:, bvKIndex(1), bvKIndex(2), bvKIndex(3), iSpin)&
            & = rhoBvK(:,:, bvKIndex(1), bvKIndex(2), bvKIndex(3), iSpin)&
            & + kWeight(iK) * real(rhoDual(:,:, iKS) * phase, dp)
      end do
    end do

  end subroutine transformDualSpaceToBvKRealSpace


  !> Make a regular density matrix for the real wave-function case
  !> Note: In order to save memory, the eigenvectors (which should be intent(in) parameters) are
  !> overwritten and then restored again
  subroutine fullDensityMatrix_real(dm, eigenvecs, filling)

    !> The resulting nOrb*nOrb density matrix
    real(dp), intent(out) :: dm(:,:)

    !> The eigenvectors of the system
    real(dp), intent(inout) :: eigenvecs(:,:)

    !> The occupation numbers of the orbitals
    real(dp), intent(in) :: filling(:)

    integer :: ii, nLevels
    real(dp) :: shift

    @:ASSERT(all(shape(eigenvecs) == shape(dm)))
    @:ASSERT(size(eigenvecs,dim=1) == size(eigenvecs,dim=2))
    @:ASSERT(size(eigenvecs,dim=1) == size(filling))

    dm(:,:) = 0.0_dp
    do ii =  size(filling), 1, -1
      nLevels = ii
      if (abs(filling(ii)) >= epsilon(1.0_dp)) then
        exit
      end if
    end do
    shift = minval(filling(1:nLevels))
    if (shift > epsilon(1.0_dp)) then
      ! all fillings are definitely positive

      !$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(ii) SCHEDULE(RUNTIME)
      do ii = 1, nLevels
        eigenvecs(:,ii) = sqrt(filling(ii)) * eigenvecs(:,ii)
      end do
      !$OMP  END PARALLEL DO

      call herk(dm, eigenvecs(:,1:nLevels))
      !$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(ii) SCHEDULE(RUNTIME)
      do ii = 1, nLevels
        eigenvecs(:,ii) = eigenvecs(:,ii) / sqrt(filling(ii))
      end do
      !$OMP  END PARALLEL DO

    else

      ! shift matrix so that filling operations are positive
      call herk(dm, eigenvecs(:,1:nLevels))
      shift = shift - arbitraryConstant
      !$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(ii) SCHEDULE(RUNTIME)
      do ii = 1, nLevels
        eigenvecs(:,ii) = sqrt(filling(ii)-shift) * eigenvecs(:,ii)
      end do
      !$OMP  END PARALLEL DO

      call herk(dm, eigenvecs(:,1:nLevels), beta=shift)
      !$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(ii) SCHEDULE(RUNTIME)
      do ii = 1, nLevels
        eigenvecs(:,ii) = eigenvecs(:,ii) / sqrt(filling(ii)-shift)
      end do
      !$OMP  END PARALLEL DO

    end if
  end subroutine fullDensityMatrix_real


  !> Make a regular density matrix for the complex wave-function case
  !> Note: in order to save memory, the eigenvectors (which should be intent(in) parameters) are
  !> overwritten and then restored again
  subroutine fullDensityMatrix_cmplx(dm, eigenvecs, filling)

    !> The resulting nOrb*nOrb density matrix
    complex(dp), intent(out) :: dm(:,:)

    !> The eigenvectors of the system
    complex(dp), intent(inout) :: eigenvecs(:,:)

    !> The occupation numbers of the orbitals
    real(dp), intent(in) :: filling(:)

    integer :: ii, nLevels
    real(dp) :: shift

    @:ASSERT(all(shape(eigenvecs) == shape(dm)))
    @:ASSERT(size(eigenvecs,dim=1) == size(eigenvecs,dim=2))
    @:ASSERT(size(eigenvecs,dim=1) == size(filling))

    dm(:,:) = cmplx(0.0_dp,0.0_dp,dp)

    do ii =  size(filling), 1, -1
      nLevels = ii
      if (abs(filling(ii)) >= epsilon(1.0_dp)) then
        exit
      end if
    end do
    shift = minval(filling(1:nLevels))
    if (shift > epsilon(1.0_dp)) then
      ! all fillings are definitely positive
      !$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(ii) SCHEDULE(RUNTIME)
      do ii = 1, nLevels
        eigenvecs(:,ii) = sqrt(filling(ii))*eigenvecs(:,ii)
      end do
      !$OMP  END PARALLEL DO

      call herk(dm, eigenvecs(:,1:nLevels))
      !$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(ii) SCHEDULE(RUNTIME)
      do ii = 1, nLevels
        eigenvecs(:,ii) = eigenvecs(:,ii) / sqrt(filling(ii))
      end do
      !$OMP  END PARALLEL DO

    else
      ! shift matrix so that filling operations are positive
      call herk(dm, eigenvecs(:,1:nLevels))
      shift = shift - arbitraryConstant
      !$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(ii) SCHEDULE(RUNTIME)
      do ii = 1, nLevels
        eigenvecs(:,ii) = sqrt(filling(ii)-shift)*eigenvecs(:,ii)
      end do
      !$OMP  END PARALLEL DO

      call herk(dm, eigenvecs(:,1:nLevels), beta=shift)
      !$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(ii) SCHEDULE(RUNTIME)
      do ii = 1, nLevels
        eigenvecs(:,ii) = eigenvecs(:,ii) / sqrt(filling(ii)-shift)
      end do
      !$OMP  END PARALLEL DO

    end if

  end subroutine fullDensityMatrix_cmplx


  !> Make an energy weighted density matrix for the real wave-function case
  !> Note: in order to save memory, the eigenvectors (which should be intent(in) parameters) are
  !> overwritten and then restored again
  subroutine fullEnergyDensityMatrix_real(dm, eigenvecs, filling, eigen)

    !> The resulting nOrb*nOrb density matrix
    real(dp), intent(out) :: dm(:,:)

    !> The eigenvectors of the system
    real(dp), intent(inout) :: eigenvecs(:,:)

    !> The occupation numbers of the orbitals
    real(dp), intent(in) :: filling(:)

    !> Eigenvalues of the system
    real(dp), intent(in) :: eigen(:)

    integer :: ii, nLevels
    real(dp) :: shift
    real(dp) :: fillProduct(size(filling))

    @:ASSERT(all(shape(eigenvecs) == shape(dm)))
    @:ASSERT(size(eigenvecs,dim=1) == size(eigenvecs,dim=2))
    @:ASSERT(size(eigenvecs,dim=1) == size(filling))
    @:ASSERT(size(eigen) == size(filling))

    dm(:,:) = 0.0_dp
    do ii =  size(filling), 1, -1
      nLevels = ii
      if (abs(filling(ii)) >= epsilon(1.0_dp)) then
        exit
      end if
    end do
    fillProduct(1:nLevels) = filling(1:nLevels) * eigen(1:nLevels)
    if ((minval(fillProduct(1:nLevels)) < 0.0_dp&
        & .eqv. maxval(fillProduct(1:nLevels)) < 0.0_dp)&
        & .and. abs(minval(fillProduct(1:nLevels))) > epsilon(1.0_dp)&
        & .and. abs(maxval(fillProduct(1:nLevels))) > epsilon(1.0_dp)) then
      ! all fillings the same sign, and fairly large

      !$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(ii) SCHEDULE(RUNTIME)
      do ii = 1, nLevels
        eigenvecs(:,ii) = sqrt(abs(fillProduct(ii)))*eigenvecs(:,ii)
      end do
      !$OMP  END PARALLEL DO

      call herk(dm, eigenvecs(:,1:nLevels), alpha=sign(1.0_dp, maxval(fillProduct(1:nLevels))))
      !$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(ii) SCHEDULE(RUNTIME)
      do ii = 1, nLevels
        eigenvecs(:,ii) = eigenvecs(:,ii) / sqrt(abs(fillProduct(ii)))
      end do
      !$OMP  END PARALLEL DO

    else

      ! shift matrix so that filling operations are positive
      call herk(dm, eigenvecs(:,1:nLevels))
      shift = minval(fillProduct(1:nLevels)) - arbitraryConstant
      !$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(ii) SCHEDULE(RUNTIME)
      do ii = 1, nLevels
        eigenvecs(:,ii) = sqrt(fillProduct(ii)-shift) * eigenvecs(:,ii)
      end do
      !$OMP  END PARALLEL DO

      call herk(dm, eigenvecs(:,1:nLevels), beta=shift)
      !$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(ii) SCHEDULE(RUNTIME)
      do ii = 1, nLevels
        eigenvecs(:,ii) = eigenvecs(:,ii) / sqrt(fillProduct(ii)-shift)
      end do
      !$OMP  END PARALLEL DO

    end if
  end subroutine fullEnergyDensityMatrix_real


  !> Make an energy weighted density matrix for the complex wave-function case
  !> Note: in order to save memory, the eigenvectors (which should be intent(in) parameters) are
  !> overwritten and then restored again
  subroutine fullEnergyDensityMatrix_cmplx(dm, eigenvecs, filling, eigen)

    !> The resulting nOrb*nOrb density matrix
    complex(dp), intent(out) :: dm(:,:)

    !> The eigenvectors of the system
    complex(dp), intent(inout) :: eigenvecs(:,:)

    !> The occupation numbers of the orbitals
    real(dp), intent(in) :: filling(:)

    !> Eigenvalues of the system
    real(dp), intent(in) :: eigen(:)

    integer :: ii, nLevels
    real(dp) :: shift
    real(dp) :: fillProduct(size(filling))

    @:ASSERT(all(shape(eigenvecs) == shape(dm)))
    @:ASSERT(size(eigenvecs,dim=1) == size(eigenvecs,dim=2))
    @:ASSERT(size(eigenvecs,dim=1) == size(filling))
    @:ASSERT(size(eigen) == size(filling))

    dm(:,:) = cmplx(0.0_dp,0.0_dp,dp)

    do ii =  size(filling), 1, -1
      nLevels = ii
      if (abs(filling(ii)) >= epsilon(1.0_dp)) then
        exit
      end if
    end do

    fillProduct(1:nLevels) = filling(1:nLevels) * eigen(1:nLevels)
    if ((minval(fillProduct(1:nLevels)) < 0.0_dp&
        & .eqv. maxval(fillProduct(1:nLevels)) < 0.0_dp)&
        & .and. abs(minval(fillProduct(1:nLevels))) > epsilon(1.0_dp)&
        & .and. abs(maxval(fillProduct(1:nLevels))) > epsilon(1.0_dp)) then
      ! all fillings the same sign, and fairly large
      !$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(ii) SCHEDULE(RUNTIME)
      do ii = 1, nLevels
        eigenvecs(:,ii) = sqrt(abs(fillProduct(ii))) * eigenvecs(:,ii)
      end do
      !$OMP  END PARALLEL DO

      call herk(dm, eigenvecs(:,1:nLevels), alpha=sign(1.0_dp, maxval(fillProduct(1:nLevels))))
      !$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(ii) SCHEDULE(RUNTIME)
      do ii = 1, nLevels
        eigenvecs(:,ii) = eigenvecs(:,ii) / sqrt(abs(fillProduct(ii)))
      end do
      !$OMP  END PARALLEL DO

    else

      ! shift matrix so that filling operations are positive
      call herk(dm, eigenvecs(:,1:nLevels))
      shift = minval(fillProduct(1:nLevels)) - arbitraryConstant
      !$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(ii) SCHEDULE(RUNTIME)
      do ii = 1, nLevels
        eigenvecs(:,ii) = sqrt(fillProduct(ii)-shift)*eigenvecs(:,ii)
      end do
      !$OMP  END PARALLEL DO

      call herk(dm, eigenvecs(:,1:nLevels), beta=shift)
      !$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(ii) SCHEDULE(RUNTIME)
      do ii = 1, nLevels
        eigenvecs(:,ii) = eigenvecs(:,ii) / sqrt(fillProduct(ii)-shift)
      end do
      !$OMP  END PARALLEL DO

    end if
  end subroutine fullEnergyDensityMatrix_cmplx


#:if WITH_SCALAPACK
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!! Scalapack routines
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Create density or energy weighted density matrix (real) for both triangles.
  subroutine makeDensityMtxRealBlacs(myBlacs, desc, filling, eigenVecs, densityMtx, eigenVals)

    !> BLACS grid information
    type(blacsgrid), intent(in) :: myBlacs

    !> Matrix descriptor
    integer, intent(in) :: desc(:)

    !> Occupations of levels
    real(dp), target, intent(in) :: filling(:)

    !> Eigenvectors of system
    real(dp), intent(inout) :: eigenVecs(:,:)

    !> Resulting (symmetric) density matrix
    real(dp), intent(out) :: densityMtx(:,:)

    !> Eigenvalues, if energy weighted density matrix required
    real(dp), intent(in), optional :: eigenVals(:)

    integer  :: ii, jj, iGlob, iLoc, blockSize
    type(blocklist) :: blocks
    real(dp), allocatable :: work(:,:)

    densityMtx(:, :) = 0.0_dp
    work = densityMtx

    ! Scale a copy of the eigenvectors
    call blocks%init(myBlacs, desc, "c")
    if (present(eigenVals)) then
      do ii = 1, size(blocks)
        call blocks%getblock(ii, iGlob, iLoc, blockSize)
        do jj = 0, blockSize - 1
          work(:, iLoc + jj) = eigenVecs(:, iLoc + jj) * eigenVals(iGlob + jj) * filling(iGlob + jj)
        end do
      end do
    else
      do ii = 1, size(blocks)
        call blocks%getblock(ii, iGlob, iLoc, blockSize)
        do jj = 0, blockSize - 1
          work(:, iLoc + jj) = eigenVecs(:, iLoc + jj) * filling(iGlob + jj)
        end do
      end do
    end if

    ! Create matrix
    call pblasfx_pgemm(eigenVecs, desc, work, desc, densityMtx, desc, transb="T")

  end subroutine makeDensityMtxRealBlacs


  !> Create density or energy weighted density matrix (complex) for both triangles.
  subroutine makeDensityMtxCplxBlacs(myBlacs, desc, filling, eigenVecs, densityMtx, eigenVals)

    !> BLACS grid information
    type(blacsgrid), intent(in) :: myBlacs

    !> Matrix descriptor
    integer, intent(in) :: desc(:)

    !> Occupations of levels
    real(dp), target, intent(in) :: filling(:)

    !> Eigenvectors of system
    complex(dp), intent(inout) :: eigenVecs(:,:)

    !> Resulting (symmetric) density matrix
    complex(dp), intent(out) :: densityMtx(:,:)

    !> Eigenvalues, if energy weighted density matrix required
    real(dp), intent(in), optional :: eigenVals(:)

    integer  :: ii, jj, iGlob, iLoc, blockSize
    type(blocklist) :: blocks
    complex(dp), allocatable :: work(:,:)

    densityMtx(:, :) = cmplx(0,0,dp)
    work = densityMtx

    ! Scale a copy of the eigenvectors
    call blocks%init(myBlacs, desc, "c")
    if (present(eigenVals)) then
      do ii = 1, size(blocks)
        call blocks%getblock(ii, iGlob, iLoc, blockSize)
        do jj = 0, blockSize - 1
          work(:, iLoc + jj) = eigenVecs(:, iLoc + jj) * eigenVals(iGlob + jj) * filling(iGlob + jj)
        end do
      end do
    else
      do ii = 1, size(blocks)
        call blocks%getblock(ii, iGlob, iLoc, blockSize)
        do jj = 0, blockSize - 1
          work(:, iLoc + jj) = eigenVecs(:, iLoc + jj) * filling(iGlob + jj)
        end do
      end do
    end if

    ! Create matrix
    call pblasfx_pgemm(eigenVecs, desc, work, desc, densityMtx, desc, transb="C")
    ! hermitian symmetrize
    work(:,:) = densityMtx
    call pblasfx_ptranc(work, desc, densityMtx, desc, alpha=(0.5_dp,0.0_dp), beta=(0.5_dp,0.0_dp))

  end subroutine makeDensityMtxCplxBlacs

#:endif

#:if WITH_MAGMA
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!! MAGMA GPU routines
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Make a regular density matrix for the real wave-function case on a single GPU
  !! Note: In order to save CPU memory, the eigenvectors (which should be intent(in) parameters) are
  !! overwritten with scaled copies and then restored again
  subroutine makeDensityMtxRealGPU(dm_CPU, eigenvecs, filling, errStatus)

    !> Resulting nOrb*nOrb density matrix
    real(dp), intent(out) :: dm_CPU(:,:)

    !> Eigenvectors of the system
    real(dp), intent(inout) :: eigenvecs(:,:)

    !> Occupation numbers of the orbitals
    real(dp), intent(in) :: filling(:)

    !> Error status
    type(TStatus), intent(out) :: errStatus

    integer :: ii, nLevels
    real(dp) :: shift

    !! BLAS specific variables
    character, parameter :: uplo = 'L'
    character, parameter :: trans = 'N'
    integer :: n, k, ldda, lddc, info

    ! MAGA queue, density matrix and eigenvectors
    ! Note: this variable type should be equivalent to the magma_devptr_t type on most systems, and
    ! will be caught at compile time otherwise.
    integer(int64) :: queue, dm_GPU, evec_GPU

    character(lc) :: strErr

    @:ASSERT(all(shape(eigenvecs) == shape(dm_CPU)))
    @:ASSERT(size(eigenvecs, dim=1) == size(eigenvecs, dim=2))
    @:ASSERT(size(eigenvecs, dim=1) == size(filling))

    dm_CPU(:,:) = 0.0_dp

    n = size(dm_CPU, dim=2)
    k = size(eigenvecs, dim=2)
    lddc = size(dm_CPU, dim=1)
    ldda = size(eigenvecs, dim=1)

    ! Allocate GPU memory
    info = magmaf_dmalloc(evec_GPU, ldda*k)
    if (info /= 0) then
      @:RAISE_FORMATTED_ERROR(errStatus, -1, "(A,I0)",&
          & 'Error: magmaf_dmalloc( evec_GPU ) failed. Info=', info)
    endif
    info = magmaf_dmalloc(dm_GPU, lddc*n)
    if (info /= 0) then
      @:RAISE_FORMATTED_ERROR(errStatus, -1, "(A,I0)",&
          & 'Error: magmaf_dmalloc( dm_GPU ) failed. Info=', info)
    endif

    call magmaf_queue_create(0, queue)

    do ii = size(filling), 1, -1
      nLevels = ii
      if (abs(filling(ii)) >= epsilon(1.0_dp)) then
        exit
      end if
    end do

    shift = minval(filling(:nLevels))

    if (shift > epsilon(1.0_dp)) then
      ! All fillings are positive definite

      !$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(ii) SCHEDULE(RUNTIME)
      do ii = 1, nLevels
        eigenvecs(:,ii) = sqrt(filling(ii)) * eigenvecs(:,ii)
      end do
      !$OMP  END PARALLEL DO

      ! Load eigenvector matrix to GPU. DM is already allocated and will be over-written by dsyrk
      call magmaf_dsetmatrix( n, nLevels, eigenvecs(:, :nLevels), ldda, evec_GPU, ldda, queue)

      ! Rank-k update
      call magmaf_dsyrk(uplo, trans, n, nLevels, 1.0_dp, evec_GPU, ldda, 0.0_dp, dm_GPU, lddc,&
          & queue)

      ! Undo eigenvector changes
      !$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(ii) SCHEDULE(RUNTIME)
      do ii = 1, nLevels
        eigenvecs(:,ii) = eigenvecs(:,ii) / sqrt(filling(ii))
      end do
      !$OMP  END PARALLEL DO

    else

      ! shift matrix so that filling occupations are positive, then un-shift to correct values

      ! Load eigenvector matrix to GPU. DM is already allocated and will be over-written by dsyrk
      call magmaf_dsetmatrix(n, nLevels, eigenvecs(:, :nLevels), ldda, evec_GPU, ldda, queue)

      ! Rank-k update for identity over the rank of the occupied orbitals (dm=1)_occ
      call magmaf_dsyrk(uplo, trans, n, nLevels, 1.0_dp, evec_GPU, ldda, 0.0_dp, dm_GPU, lddc,&
          & queue)

      shift = shift - arbitraryConstant
      !$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(ii) SCHEDULE(RUNTIME)
      do ii = 1, nLevels
        eigenvecs(:,ii) = sqrt(filling(ii)-shift) * eigenvecs(:,ii)
      end do
      !$OMP  END PARALLEL DO

      ! Load eigenvectors with modified scaling into GPU
      call magmaf_dsetmatrix(n, nLevels, eigenvecs(:, :nLevels), ldda, evec_GPU, ldda, queue)

      ! Rank-k to add to existing occupied identity with shifted eigenvectors
      ! dm = |psi_i> (occ_i-shift) <psi_i| + shift * 1_i
      call magmaf_dsyrk( uplo, trans, n, nLevels, 1.0_dp, evec_GPU, ldda, shift, dm_GPU, lddc,&
          & queue)

      ! undo eigenvector change
      !$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(ii) SCHEDULE(RUNTIME)
      do ii = 1, nLevels
        eigenvecs(:,ii) = eigenvecs(:,ii) / sqrt(filling(ii)-shift)
      end do
      !$OMP  END PARALLEL DO

    end if

    ! Get dm_GPU back to CPU
    call magmaf_dgetmatrix( n, n, dm_GPU, lddc, dm_CPU, lddc, queue)

    ! Free GPU memory
    info = magmaf_free( evec_GPU )
    if (info /= 0) then
      @:RAISE_FORMATTED_ERROR(errStatus, -1, "(A,I0)",&
          & 'Error: magmaf_free( evec_GPU ) failed. Info=', info)
    endif
    info = magmaf_free( dm_GPU )
    if (info /= 0) then
      @:RAISE_FORMATTED_ERROR(errStatus, -1, "(A,I0)",&
          & 'Error: magmaf_free( dm_GPU ) failed. Info=', info)
    endif

    call magmaf_queue_destroy(queue)

  end subroutine makeDensityMtxRealGPU


  !> Make an energy weighted density matrix for the real wave-function case on a single GPU
  !! Note: In order to save CPU memory, the eigenvectors (which should be intent(in) parameters) are
  !! overwritten with scaled copies and then restored again
  subroutine makeEnergyDensityMtxRealGPU(edm_CPU, eigenvecs, filling, eigen, errStatus)

    !> Resulting nOrb*nOrb energy weighted density matrix
    real(dp), intent(out) :: edm_CPU(:,:)

    !> Eigenvectors of the system
    real(dp), intent(inout) :: eigenvecs(:,:)

    !> Occupation numbers of the orbitals
    real(dp), intent(in) :: filling(:)

    !> Eigenvalues of the system
    real(dp), intent(in) :: eigen(:)

    !> Error status
    type(TStatus), intent(out) :: errStatus

    integer :: ii, nLevels
    real(dp) :: shift, fillProduct(size(filling))

    ! BLAS specific variables
    character, parameter :: uplo = 'L'
    character, parameter :: trans = 'N'
    integer :: n, k, ldda, lddc, info

    ! MAGA queue, density matrix and eigenvectors
    ! Note: this variable type should be equivalent to the magma_devptr_t type on most systems, and
    ! will be caught at compile time otherwise.
    integer(int64) :: queue, edm_GPU, evec_GPU

    character(lc) :: strErr

    @:ASSERT(all(shape(eigenvecs) == shape(edm_CPU)))
    @:ASSERT(size(eigenvecs, dim=1) == size(eigenvecs, dim=2))
    @:ASSERT(size(eigenvecs, dim=1) == size(filling))
    @:ASSERT(size(eigen) == size(filling))

    edm_CPU(:,:) = 0.0_dp

    n = size(edm_CPU, dim=2)
    k = size(eigenvecs, dim=2)
    lddc = size(edm_CPU, dim=1)
    ldda = size(eigenvecs, dim=1)

    ! Allocate GPU memory
    info = magmaf_dmalloc(evec_GPU, ldda*k)
    if (info /= 0) then
      @:RAISE_FORMATTED_ERROR(errStatus, -1, "(A,I0)",&
          & 'Error: magmaf_dmalloc( evec_GPU ) failed. Info=', info)
    endif
    info = magmaf_dmalloc(edm_GPU, lddc*n)
    if (info /= 0) then
      @:RAISE_FORMATTED_ERROR(errStatus, -1, "(A,I0)",&
          & 'Error: magmaf_dmalloc( edm_GPU ) failed. Info=', info)
    endif

    call magmaf_queue_create(0, queue)

    do ii = size(filling), 1, -1
      nLevels = ii
      if (abs(filling(ii)) >= epsilon(1.0_dp)) then
        exit
      end if
    end do

    fillProduct(:nLevels) = filling(:nLevels) * eigen(:nLevels)
    if ((minval(fillProduct(:nLevels)) < 0.0_dp&
        & .eqv. maxval(fillProduct(:nLevels)) < 0.0_dp)&
        & .and. abs(minval(fillProduct(:nLevels))) > epsilon(1.0_dp)&
        & .and. abs(maxval(fillProduct(:nLevels))) > epsilon(1.0_dp)) then
      ! all fillings the same sign, and fairly large

      !$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(ii) SCHEDULE(RUNTIME)
      do ii = 1, nLevels
        eigenvecs(:,ii) = sqrt(abs(fillProduct(ii))) * eigenvecs(:,ii)
      end do
      !$OMP  END PARALLEL DO

      ! Load eigenvector matrix to GPU. DM is already allocated and will be over-written by dsyrk
      call magmaf_dsetmatrix( n, nLevels, eigenvecs(:, :nLevels), ldda, evec_GPU, ldda, queue)

      ! Rank-k update
      call magmaf_dsyrk(uplo, trans, n, nLevels,&
          & sign(1.0_dp, maxval(fillProduct(1:nLevels))), evec_GPU, ldda, 0.0_dp, edm_GPU, lddc,&
          & queue)

      ! Undo eigenvector changes
      !$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(ii) SCHEDULE(RUNTIME)
      do ii = 1, nLevels
        eigenvecs(:,ii) = eigenvecs(:,ii) / sqrt(abs(fillProduct(ii)))
      end do
      !$OMP  END PARALLEL DO

    else

      ! shift matrix so that filling occupations are positive, then un-shift to correct values

      ! Load eigenvector matrix to GPU. EDM is already allocated and will be over-written by dsyrk
      call magmaf_dsetmatrix(n, nLevels, eigenvecs(:, :nLevels), ldda, evec_GPU, ldda, queue)

      ! Rank-k update for identity over the rank of the occupied orbitals (edm=1)_occ
      call magmaf_dsyrk(uplo, trans, n, nLevels, 1.0_dp, evec_GPU, ldda, 0.0_dp, edm_GPU, lddc,&
          & queue)

      shift = minval(fillProduct(:nLevels)) - arbitraryConstant

      !$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(ii) SCHEDULE(RUNTIME)
      do ii = 1, nLevels
        eigenvecs(:,ii) = sqrt(fillProduct(ii)-shift) * eigenvecs(:,ii)
      end do
      !$OMP  END PARALLEL DO

      ! Load eigenvectors with modified scaling into GPU
      call magmaf_dsetmatrix(n, nLevels, eigenvecs(:, :nLevels), ldda, evec_GPU, ldda, queue)

      ! Rank-k to add to existing occupied identity with shifted eigenvectors, giving
      ! edm = |psi_i> (e_i*occ_i-shift) <psi_i| + shift * 1_i
      call magmaf_dsyrk( uplo, trans, n, nLevels, 1.0_dp, evec_GPU, ldda, shift, edm_GPU, lddc,&
          & queue)

      ! undo eigenvector change
      !$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(ii) SCHEDULE(RUNTIME)
      do ii = 1, nLevels
        eigenvecs(:,ii) = eigenvecs(:,ii) / sqrt(fillProduct(ii)-shift)
      end do
      !$OMP  END PARALLEL DO

    end if

    ! Get edm_GPU back to CPU
    call magmaf_dgetmatrix( n, n, edm_GPU, lddc, edm_CPU, lddc, queue)

    ! Free GPU memory
    info = magmaf_free( evec_GPU )
    if (info /= 0) then
      @:RAISE_FORMATTED_ERROR(errStatus, -1, "(A,I0)",&
          & 'Error: magmaf_free( evec_GPU ) failed. Info=', info)
    endif
    info = magmaf_free( edm_GPU )
    if (info /= 0) then
      @:RAISE_FORMATTED_ERROR(errStatus, -1, "(A,I0)",&
          & 'Error: magmaf_free( edm_GPU ) failed. Info=', info)
    endif

    call magmaf_queue_destroy(queue)

  end subroutine makeEnergyDensityMtxRealGPU


  !> Make a regular density matrix for the complex wave-function case on the GPU
  !! Note: In order to save memory, the eigenvectors (which should actually be intent(in)
  !! parameters) are overwritten and then restored again
  subroutine makeDensityMtxCmplxGPU(dm_CPU, eigenvecs, filling, errStatus)

    !> Resulting nOrb*nOrb density matrix
    complex(dp), intent(out) :: dm_CPU(:,:)

    !> Eigenvectors of the system
    complex(dp), intent(inout) :: eigenvecs(:,:)

    !> Occupation numbers of the orbitals
    real(dp), intent(in) :: filling(:)

    !> Error status
    type(TStatus), intent(out) :: errStatus

    integer :: ii, nLevels
    real(dp) :: shift

    !! BLAS specific variables
    character, parameter :: uplo = 'L'
    character, parameter :: trans = 'N'
    integer :: n, k, ldda, lddc, info

    ! MAGA queue, density matrix and eigenvectors
    ! Note: this variable type should be equivalent to the magma_devptr_t type on most systems, and
    ! will be caught at compile time otherwise.
    integer(int64) :: queue, dm_GPU, evec_GPU

    character(lc) :: strErr

    @:ASSERT(all(shape(eigenvecs) == shape(dm_CPU)))
    @:ASSERT(size(eigenvecs, dim=1) == size(eigenvecs, dim=2))
    @:ASSERT(size(eigenvecs, dim=1) == size(filling))

    dm_CPU(:,:) = cmplx(0, 0, dp)

    n = size(dm_CPU, dim=2)
    k = size(eigenvecs, dim=2)
    lddc = size(dm_CPU, dim=1)
    ldda = size(eigenvecs, dim=1)

    ! Allocate GPU memory
    info = magmaf_zmalloc(evec_GPU, ldda*k)
    if (info /= 0) then
      @:RAISE_FORMATTED_ERROR(errStatus, -1, "(A,I0)",&
          & 'Error: magmaf_zmalloc( evec_GPU ) failed. Info=', info)
    endif
    info = magmaf_zmalloc(dm_GPU, lddc*n)
    if (info /= 0) then
      @:RAISE_FORMATTED_ERROR(errStatus, -1, "(A,I0)",&
          & 'Error: magmaf_zmalloc( dm_GPU ) failed. Info=', info)
    endif

    call magmaf_queue_create(0, queue)

    do ii = size(filling), 1, -1
      nLevels = ii
      if (abs(filling(ii)) >= epsilon(1.0_dp)) then
        exit
      end if
    end do

    shift = minval(filling(:nLevels))

    if (shift > epsilon(1.0_dp)) then
      ! All fillings are definitely positve

      !$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(ii) SCHEDULE(RUNTIME)
      do ii = 1, nLevels
        eigenvecs(:,ii) = sqrt(filling(ii)) * eigenvecs(:,ii)
      end do
      !$OMP  END PARALLEL DO

      ! Load eigenvector matrix to GPU. DM is already allocated and will be over-written by dsyrk
      call magmaf_zsetmatrix( n, nLevels, eigenvecs(:, :nLevels), ldda, evec_GPU, ldda, queue)

      ! Rank-k update
      call magmaf_zherk(uplo, trans, n, nLevels, 1.0_dp, evec_GPU, ldda, 0.0_dp, dm_GPU, lddc,&
          & queue)

      ! Undo eigenvector changes
      !$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(ii) SCHEDULE(RUNTIME)
      do ii = 1, nLevels
        eigenvecs(:,ii) = eigenvecs(:,ii) / sqrt(filling(ii))
      end do
      !$OMP  END PARALLEL DO

    else

      ! shift matrix so that filling occupations are positive, then un-shift to correct values

      ! Load eigenvector matrix to GPU. DM is already allocated and will be over-written by dsyrk
      call magmaf_zsetmatrix(n, nLevels, eigenvecs(:, :nLevels), ldda, evec_GPU, ldda, queue)

      ! Rank-k update for identity over the rank of the occupied orbitals (dm=1)_occ
      call magmaf_zherk(uplo, trans, n, nLevels, 1.0_dp, evec_GPU, ldda, 0.0_dp, dm_GPU, lddc,&
          & queue)

      shift = shift - arbitraryConstant
      !$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(ii) SCHEDULE(RUNTIME)
      do ii = 1, nLevels
        eigenvecs(:,ii) = sqrt(filling(ii)-shift) * eigenvecs(:,ii)
      end do
      !$OMP  END PARALLEL DO

      ! Load eigenvectors with modified scaling into GPU
      call magmaf_zsetmatrix(n, nLevels, eigenvecs(:, :nLevels), ldda, evec_GPU, ldda, queue)

      ! Rank-k to add to existing occupied identity with shifted eigenvectors
      ! dm = |psi_i> (occ_i-shift) <psi_i| + shift * 1_i
      call magmaf_zherk( uplo, trans, n, nLevels, 1.0_dp, evec_GPU, ldda, shift, dm_GPU, lddc,&
          & queue)

      ! undo eigenvector change
      !$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(ii) SCHEDULE(RUNTIME)
      do ii = 1, nLevels
        eigenvecs(:,ii) = eigenvecs(:,ii) / sqrt(filling(ii)-shift)
      end do
      !$OMP  END PARALLEL DO

    end if

    ! Get dm_GPU back to CPU
    call magmaf_zgetmatrix( n, n, dm_GPU, lddc, dm_CPU, lddc, queue)

    ! Free GPU memory
    info = magmaf_free( evec_GPU )
    if (info /= 0) then
      @:RAISE_FORMATTED_ERROR(errStatus, -1, "(A,I0)",&
          & 'Error: magmaf_free( evec_GPU ) failed. Info=', info)
    endif
    info = magmaf_free( dm_GPU )
    if (info /= 0) then
      @:RAISE_FORMATTED_ERROR(errStatus, -1, "(A,I0)",&
          & 'Error: magmaf_free( dm_GPU ) failed. Info=', info)
    endif

    call magmaf_queue_destroy(queue)

  end subroutine makeDensityMtxCmplxGPU


  !> Make an energy weighted density matrix for the complex wave-function case on the GPU
  !! Note: In order to save memory, the eigenvectors (which should actually be intent(in)
  !! parameters) are overwritten and then restored again
  subroutine makeEnergyDensityMtxCmplxGPU(edm_CPU, eigenvecs, filling, eigenvals, errStatus)

    !> Resulting nOrb*nOrb energy weighted density matrix
    complex(dp), intent(out) :: edm_CPU(:,:)

    !> Eigenvectors of the system
    complex(dp), intent(inout) :: eigenvecs(:,:)

    !> Occupation numbers of the orbitals
    real(dp), intent(in) :: filling(:)

    !> Eigenvalues of the system
    real(dp), intent(in) :: eigenvals(:)

    !> Error status
    type(TStatus), intent(out) :: errStatus

    integer :: ii, nLevels
    real(dp) :: shift, fillProduct(size(filling))

    ! BLAS specific variables
    character, parameter :: uplo = 'L'
    character, parameter :: trans = 'N'
    integer :: n, k, ldda, lddc, info

    ! MAGA queue, density matrix and eigenvectors
    ! Note: this variable type should be equivalent to the magma_devptr_t type on most systems, and
    ! will be caught at compile time otherwise.
    integer(int64) :: queue, edm_GPU, evec_GPU

    character(lc) :: strErr

    @:ASSERT(all(shape(eigenvecs) == shape(edm_CPU)))
    @:ASSERT(size(eigenvecs, dim=1) == size(eigenvecs, dim=2))
    @:ASSERT(size(eigenvecs, dim=1) == size(filling))
    @:ASSERT(size(eigenvals) == size(filling))

    edm_CPU(:,:) = cmplx(0, 0, dp)

    n = size(edm_CPU, dim=2)
    k = size(eigenvecs, dim=2)
    lddc = size(edm_CPU, dim=1)
    ldda = size(eigenvecs, dim=1)

    ! Allocate GPU memory
    info = magmaf_zmalloc(evec_GPU, ldda*k)
    if (info /= 0) then
      @:RAISE_FORMATTED_ERROR(errStatus, -1, "(A,I0)",&
          & 'Error: magmaf_zmalloc( evec_GPU ) failed. Info=', info)
    endif
    info = magmaf_zmalloc(edm_GPU, lddc*n)
    if (info /= 0) then
      @:RAISE_FORMATTED_ERROR(errStatus, -1, "(A,I0)",&
          & 'Error: magmaf_zmalloc( edm_GPU ) failed. Info=', info)
    endif

    call magmaf_queue_create(0, queue)

    do ii = size(filling), 1, -1
      nLevels = ii
      if (abs(filling(ii)) >= epsilon(1.0_dp)) then
        exit
      end if
    end do

    fillProduct(:nLevels) = filling(:nLevels) * eigenvals(:nLevels)
    if ((minval(fillProduct(:nLevels)) < 0.0_dp&
        & .eqv. maxval(fillProduct(:nLevels)) < 0.0_dp)&
        & .and. abs(minval(fillProduct(:nLevels))) > epsilon(1.0_dp)&
        & .and. abs(maxval(fillProduct(:nLevels))) > epsilon(1.0_dp)) then
      ! all fillings the same sign, and fairly large

      !$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(ii) SCHEDULE(RUNTIME)
      do ii = 1, nLevels
        eigenvecs(:,ii) = sqrt(abs(fillProduct(ii))) * eigenvecs(:,ii)
      end do
      !$OMP  END PARALLEL DO

      ! Load eigenvector matrix to GPU. DM is already allocated and will be over-written by dsyrk
      call magmaf_zsetmatrix( n, nLevels, eigenvecs(:, :nLevels), ldda, evec_GPU, ldda, queue)

      ! Rank-k update
      call magmaf_zherk(uplo, trans, n, nLevels,&
          & sign(1.0_dp, maxval(fillProduct(1:nLevels))), evec_GPU, ldda, 0.0_dp, edm_GPU, lddc,&
          & queue)

      ! Undo eigenvector changes
      !$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(ii) SCHEDULE(RUNTIME)
      do ii = 1, nLevels
        eigenvecs(:,ii) = eigenvecs(:,ii) / sqrt(abs(fillProduct(ii)))
      end do
      !$OMP  END PARALLEL DO

    else

      ! shift matrix so that filling occupations are positive, then un-shift to correct values

      ! Load eigenvector matrix to GPU. EDM is already allocated and will be over-written by dsyrk
      call magmaf_zsetmatrix(n, nLevels, eigenvecs(:, :nLevels), ldda, evec_GPU, ldda, queue)

      ! Rank-k update for identity over the rank of the occupied orbitals (edm=1)_occ
      call magmaf_zherk(uplo, trans, n, nLevels, 1.0_dp, evec_GPU, ldda, 0.0_dp, edm_GPU, lddc,&
          & queue)

      shift = minval(fillProduct(:nLevels)) - arbitraryConstant

      !$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(ii) SCHEDULE(RUNTIME)
      do ii = 1, nLevels
        eigenvecs(:,ii) = sqrt(fillProduct(ii)-shift) * eigenvecs(:,ii)
      end do
      !$OMP  END PARALLEL DO

      ! Load eigenvectors with modified scaling into GPU
      call magmaf_zsetmatrix(n, nLevels, eigenvecs(:, :nLevels), ldda, evec_GPU, ldda, queue)

      ! Rank-k to add to existing occupied identity with shifted eigenvectors, giving
      ! edm = |psi_i> (e_i*occ_i-shift) <psi_i| + shift * 1_i
      call magmaf_zherk( uplo, trans, n, nLevels, 1.0_dp, evec_GPU, ldda, shift, edm_GPU, lddc,&
          & queue)

      ! undo eigenvector change
      !$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(ii) SCHEDULE(RUNTIME)
      do ii = 1, nLevels
        eigenvecs(:,ii) = eigenvecs(:,ii) / sqrt(fillProduct(ii)-shift)
      end do
      !$OMP  END PARALLEL DO

    end if

    ! Get edm_GPU back to CPU
    call magmaf_zgetmatrix( n, n, edm_GPU, lddc, edm_CPU, lddc, queue)

    ! Free GPU memory
    info = magmaf_free( evec_GPU )
    if (info /= 0) then
      @:RAISE_FORMATTED_ERROR(errStatus, -1, "(A,I0)",&
          & 'Error: magmaf_free( evec_GPU ) failed. Info=', info)
    endif
    info = magmaf_free( edm_GPU )
    if (info /= 0) then
      @:RAISE_FORMATTED_ERROR(errStatus, -1, "(A,I0)",&
          & 'Error: magmaf_free( edm_GPU ) failed. Info=', info)
    endif

    call magmaf_queue_destroy(queue)

  end subroutine makeEnergyDensityMtxCmplxGPU

#:endif

end module dftbp_dftb_densitymatrix
