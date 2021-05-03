!--------------------------------------------------------------------------------------------------!
!  DFTB+: general package for performing fast atomistic simulations                                !
!  Copyright (C) 2006 - 2020  DFTB+ developers group                                               !
!                                                                                                  !
!  See the LICENSE file for terms of usage and distribution.                                       !
!--------------------------------------------------------------------------------------------------!

#:include 'common.fypp'

!> Module for linear response calculations
module dftbp_linearresponse
  use dftbp_accuracy, only : dp
  use dftbp_commontypes, only : TOrbitals
  use dftbp_rangeseparated, only : TRangeSepFunc
  use dftbp_fermihelper, only : theta, deltamn, invDiff
  use dftbp_environment, only : TEnvironment
  use dftbp_periodic, only : TNeighbourList
  use dftbp_densedescr, only : TDenseDescr
  use dftbp_rotatedegen, only : TRotateDegen, TRotateDegen_init
  use dftbp_parallelks, only : TParallelKS
#:if WITH_SCALAPACK
  use dftbp_sparse2dense, only : unpackHSRealBlacs, packRhoRealBlacs
  use dftbp_sparse2dense, only : unpackHPauliBlacs, packRhoPauliBlacs
  use dftbp_scalapackfx, only : CSRC_, DLEN_, MB_, NB_, RSRC_, pblasfx_pgemm, pblasfx_ptranc
  use dftbp_scalapackfx, only : pblasfx_ptran, pblasfx_phemm, pblasfx_psymm, scalafx_indxl2g
#:else
  use dftbp_sparse2dense, only : unpackHS, packHS, packHSPauli, unpackHPauli, packHSPauliImag
  use dftbp_blasroutines, only : hemm, symm
#:endif

  implicit none

  private
  public :: dRhoStaticReal, dRhoFermiChangeStaticReal, dRhoStaticPauli, dRhoFermiChangeStaticPauli

  !> small complex value for frequency dependent cases
  complex(dp), parameter :: eta = (0.0_dp,1.0E-8_dp)

contains

  !> Calculate the derivative of density matrix from derivative of hamiltonian in static case at
  !> q=0, k=0
  subroutine dRhoStaticReal(env, dHam, neighbourList, nNeighbourSK, iSparseStart, img2CentCell,&
      & denseDesc, iKS, parallelKS, nFilled, nEmpty, eigVecsReal, eigVals, Ef, tempElec, orb,&
      & dRhoSparse, dRhoSqr, rangeSep, over, nNeighbourLC, transform,&
    #:if WITH_SCALAPACK
      & desc,&
    #:endif
      & dEi, dPsi)

    !> Environment settings
    type(TEnvironment), intent(inout) :: env

    !> Derivative of the hamiltonian
    real(dp), intent(in) :: dHam(:,:)

    !> list of neighbours for each atom
    type(TNeighbourList), intent(in) :: neighbourList

    !> Number of neighbours for each of the atoms
    integer, intent(in) :: nNeighbourSK(:)

    !> Index array for the start of atomic blocks in sparse arrays
    integer, intent(in) :: iSparseStart(:,:)

    !> map from image atoms to the original unique atom
    integer, intent(in) :: img2CentCell(:)

    !> Dense matrix descriptor
    type(TDenseDescr), intent(in) :: denseDesc

    !> Particular spin/k-point
    integer, intent(in) :: iKS

    !> K-points and spins to process
    type(TParallelKS), intent(in) :: parallelKS

    !> ground state eigenvectors
    real(dp), intent(in) :: eigVecsReal(:,:,:)

    !> Eigenvalue of each level, kpoint and spin channel
    real(dp), intent(in) :: eigvals(:,:,:)

    !> Fermi level(s)
    real(dp), intent(in) :: Ef(:)

    !> Last (partly) filled level in each spin channel
    integer, intent(in) :: nFilled(:)

    !> First (partly) empty level in each spin channel
    integer, intent(in) :: nEmpty(:)

    !> Electron temperature
    real(dp), intent(in) :: tempElec

    !> Atomic orbital information
    type(TOrbitals), intent(in) :: orb

    !> returning dRhoSparse on exit
    real(dp), intent(out) :: dRhoSparse(:)

    !> Derivative of rho as a square matrix, if needed
    real(dp), pointer :: dRhoSqr(:,:,:)

    !> Data for range-separated calculation
    type(TRangeSepFunc), allocatable, intent(inout) :: rangeSep

    !> sparse overlap matrix
    real(dp), intent(in) :: over(:)

    !> Number of neighbours for each of the atoms for the exchange contributions in the long range
    !> functional
    integer, intent(inout), allocatable :: nNeighbourLC(:)

    !> Transformation structure for degenerate orbitals
    type(TRotateDegen), intent(inout) :: transform

  #:if WITH_SCALAPACK
    !> BLACS matrix descriptor
    integer, intent(in) :: desc(DLEN_)
  #:endif

    !> Derivative of single particle eigenvalues
    real(dp), allocatable, intent(inout) :: dEi(:,:,:)

    !> Optional derivatives of single particle wavefunctions
    real(dp), allocatable, intent(inout) :: dPsi(:,:,:)

  #:if WITH_SCALAPACK
    integer :: iGlob, jGlob
  #:else
    integer :: iFilled, iEmpty, nOrb
  #:endif
    integer :: ii, jj, iS, iK

    real(dp) :: workLocal(size(eigVecsReal,dim=1), size(eigVecsReal,dim=2))
    real(dp), allocatable :: dRho(:,:)
    real(dp), allocatable :: eigVecsTransformed(:,:)
    logical :: isTransformed

    iK = parallelKS%localKS(1, iKS)
    iS = parallelKS%localKS(2, iKS)

    if (allocated(dEi)) then
      dEi(:, iK, iS) = 0.0_dp
    end if
    if (allocated(dPsi)) then
      dPsi(:, :, iS) = 0.0_dp
    end if

    workLocal(:,:) = 0.0_dp
    allocate(dRho(size(eigVecsReal,dim=1), size(eigVecsReal,dim=2)))
    dRho(:,:) = 0.0_dp

  #:if WITH_SCALAPACK

    ! dH in square form
    call unpackHSRealBlacs(env%blacs, dHam(:,iS), neighbourList%iNeighbour, nNeighbourSK,&
        & iSparseStart, img2CentCell, denseDesc, workLocal)

    ! dH times c_i
    call pblasfx_psymm(workLocal, denseDesc%blacsOrbSqr, eigVecsReal(:,:,iKS),&
        & denseDesc%blacsOrbSqr, dRho, denseDesc%blacsOrbSqr)

    ! c_i times dH times c_i
    call pblasfx_pgemm(eigVecsReal(:,:,iKS), denseDesc%blacsOrbSqr, dRho, denseDesc%blacsOrbSqr,&
        & workLocal, denseDesc%blacsOrbSqr, transa="T")

    eigvecsTransformed = eigVecsReal(:,:,iKS)
    call transform%generateUnitary(env, worklocal, eigvals(:,1,iS), eigVecsTransformed, denseDesc,&
        & isTransformed)
    ! now have states orthogonalised agains the operator in degenerate cases, |c~>

    if (isTransformed) then
      ! re-form |c~> H' <c~| with the transformed vectors

      dRho(:,:) = 0.0_dp
      workLocal(:,:) = 0.0_dp

      ! dH in square form
      call unpackHSRealBlacs(env%blacs, dHam(:,iS), neighbourList%iNeighbour, nNeighbourSK,&
          & iSparseStart, img2CentCell, denseDesc, workLocal)

      ! dH times c_i
      call pblasfx_psymm(workLocal, denseDesc%blacsOrbSqr, eigVecsTransformed,&
          & denseDesc%blacsOrbSqr, dRho, denseDesc%blacsOrbSqr)

      ! c_i times dH times c_i
      call pblasfx_pgemm(eigVecsTransformed, denseDesc%blacsOrbSqr, dRho, denseDesc%blacsOrbSqr,&
          & workLocal, denseDesc%blacsOrbSqr, transa="T")

    end if

    ! derivative of eigenvalues stored in diagonal of matrix workLocal, from <c|h'|c>
    if (allocated(dEi)) then
      do jj = 1, size(workLocal,dim=2)
        jGlob = scalafx_indxl2g(jj, desc(NB_), env%blacs%orbitalGrid%mycol, desc(CSRC_),&
            & env%blacs%orbitalGrid%ncol)
        do ii = 1, size(workLocal,dim=1)
          iGlob = scalafx_indxl2g(ii, desc(MB_), env%blacs%orbitalGrid%myrow, desc(RSRC_),&
              & env%blacs%orbitalGrid%nrow)
          if (iGlob == jGlob) then
            !if (iGlob == jGlob) then workLocal(ii,jj) contains a derivative of an eigenvalue
            dEi(iGlob, iK, iS) = workLocal(ii,jj)
          end if
        end do
      end do
    end if

    ! weight matrix with inverse of energy differences
    do jj = 1, size(workLocal,dim=2)
      jGlob = scalafx_indxl2g(jj, desc(NB_), env%blacs%orbitalGrid%mycol, desc(CSRC_),&
          & env%blacs%orbitalGrid%ncol)
      do ii = 1, size(workLocal,dim=1)
        iGlob = scalafx_indxl2g(ii, desc(MB_), env%blacs%orbitalGrid%myrow, desc(RSRC_),&
            & env%blacs%orbitalGrid%nrow)
        if (abs(eigvals(jGlob,iK,iS) - eigvals(iGlob,iK,iS)) < epsilon(0.0_dp)&
            & .and. iGlob /= jGlob) then
          ! degenerate, so no contribution
          workLocal(ii,jj) = 0.0_dp
        else
          workLocal(ii,jj) = workLocal(ii,jj) * &
              & invDiff(eigvals(jGlob,1,iS),eigvals(iGlob,1,iS),Ef(iS),tempElec)&
              & * theta(eigvals(jGlob,1,iS),eigvals(iGlob,1,iS),tempElec)
        end if
      end do
    end do

    ! Derivatives of states
    call pblasfx_pgemm(eigvecsTransformed, denseDesc%blacsOrbSqr, workLocal, denseDesc%blacsOrbSqr,&
        & dRho, denseDesc%blacsOrbSqr)

    if (allocated(dPsi)) then
      dPsi(:, :, iS) = workLocal
    end if

    ! Form derivative of occupied density matrix
    call pblasfx_pgemm(dRho, denseDesc%blacsOrbSqr, eigvecsTransformed, denseDesc%blacsOrbSqr,&
        & workLocal, denseDesc%blacsOrbSqr, transb="T", kk=nFilled(iS))

    dRho(:,:) = workLocal
    ! and symmetrize
    call pblasfx_ptran(workLocal, denseDesc%blacsOrbSqr, dRho, denseDesc%blacsOrbSqr, beta=1.0_dp)

  #:else

    ! serial case
    nOrb = size(dRho, dim = 1)

    ! dH matrix in square form
    call unpackHS(dRho, dHam(:,iS), neighbourList%iNeighbour, nNeighbourSK,&
        & denseDesc%iAtomStart, iSparseStart, img2CentCell)

    if (allocated(rangeSep)) then
      call unpackHS(workLocal, over, neighbourList%iNeighbour, nNeighbourSK,&
          & denseDesc%iAtomStart, iSparseStart, img2CentCell)
      call rangeSep%addLRHamiltonian(env, dRhoSqr(:,:,iS), over, neighbourList%iNeighbour,&
          & nNeighbourLC, denseDesc%iAtomStart, iSparseStart, orb, dRho, workLocal)
    end if

    ! form |c> H' <c|
    call symm(workLocal, 'l', dRho, eigVecsReal(:,:,iS))
    workLocal(:,:) = matmul(transpose(eigVecsReal(:,:,iS)), workLocal)

    ! orthogonalise degenerate states against perturbation, producing |c~> H' <c~|
    call transform%generateUnitary(workLocal, eigvals(:,iK,iS))
    call transform%degenerateTransform(workLocal)

    ! diagonal elements of workLocal are now derivatives of eigenvalues if needed
    if (allocated(dEi)) then
      do ii = 1, nOrb
        dEi(ii, iK, iS) = workLocal(ii,ii)
      end do
    end if

    ! Form actual perturbation U matrix for eigenvectors (potentially at finite T) by weighting
    ! the elements
    do iFilled = 1, nFilled(iS)
      do iEmpty = nEmpty(iS), nOrb
        if (.not.transform%degenerate(iFilled,iEmpty) .or. iEmpty == iFilled) then
          workLocal(iEmpty, iFilled) = workLocal(iEmpty, iFilled) * &
              & invDiff(eigvals(iFilled, iK, iS), eigvals(iEmpty, iK, iS), Ef(iS), tempElec)&
              & *theta(eigvals(iFilled, iK, iS), eigvals(iEmpty, iK, iS), tempElec)
        else
          ! rotation should already have set these elements to zero
          workLocal(iEmpty, iFilled) = 0.0_dp
        end if
      end do
    end do

    eigvecsTransformed = eigVecsReal(:,:,iS)
    call transform%applyUnitary(eigvecsTransformed)

    ! calculate the derivatives of the eigenvectors
    workLocal(:, :nFilled(iS)) =&
        & matmul(eigvecsTransformed(:, nEmpty(iS):), workLocal(nEmpty(iS):, :nFilled(iS)))

    if (allocated(dPsi)) then
      dPsi(:, :, iS) = workLocal
    end if

    ! zero the uncalculated virtual states
    workLocal(:, nFilled(iS)+1:) = 0.0_dp

    ! form the derivative of the density matrix
    dRho(:,:) = matmul(workLocal(:, :nFilled(iS)), transpose(eigvecsTransformed(:, :nFilled(iS))))&
        & + matmul(eigvecsTransformed(:, :nFilled(iS)), transpose(workLocal(:, :nFilled(iS))))

  #:endif

    dRhoSparse(:) = 0.0_dp
  #:if WITH_SCALAPACK
    call packRhoRealBlacs(env%blacs, denseDesc, dRho, neighbourList%iNeighbour, nNeighbourSK,&
        & orb%mOrb, iSparseStart, img2CentCell, dRhoSparse)
  #:else
    call packHS(dRhoSparse, dRho, neighbourList%iNeighbour, nNeighbourSK, orb%mOrb,&
        & denseDesc%iAtomStart, iSparseStart, img2CentCell)
  #:endif

    if (associated(dRhoSqr)) then
      dRhoSqr(:,:,iS) = dRho
    end if

  end subroutine dRhoStaticReal


  !> Calculate the change in the density matrix due to shift in the Fermi energy
  subroutine dRhoFermiChangeStaticReal(dRhoExtra, env, parallelKS, iKS, neighbourList,&
      & nNEighbourSK, img2CentCell, iSparseStart, dEfdE, Ef, nFilled, nEmpty, eigVecsReal, orb,&
      & denseDesc, tempElec, eigVals, dRhoSqr&
    #:if WITH_SCALAPACK
      &, desc&
    #:endif
      &)

    !> Additional contribution to the density matrix to cancel effect of Fermi energy change
    real(dp), intent(inout) :: dRhoExtra(:)

    !> Environment settings
    type(TEnvironment), intent(inout) :: env

    !> K-points and spins to process
    type(TParallelKS), intent(in) :: parallelKS

    !> spin/kpoint channel
    integer, intent(in) :: iKS

    !> list of neighbours for each atom
    type(TNeighbourList), intent(in) :: neighbourList

    !> Number of neighbours for each of the atoms
    integer, intent(in) :: nNeighbourSK(:)

    !> map from image atoms to the original unique atom
    integer, intent(in) :: img2CentCell(:)

    !> Index array for the start of atomic blocks in sparse arrays
    integer, intent(in) :: iSparseStart(:,:)

    !> Fermi level derivative
    real(dp), intent(in) :: dEfdE(:)

    !> Fermi level
    real(dp), intent(in) :: Ef(:)

    !> Last (partly) filled level in each spin channel
    integer, intent(in) :: nFilled(:)

    !> First (partly) empty level in each spin channel
    integer, intent(in) :: nEmpty(:)

    !> ground state eigenvectors
    real(dp), intent(in) :: eigVecsReal(:,:,:)

    !> Atomic orbital information
    type(TOrbitals), intent(in) :: orb

    !> Dense matrix descriptor
    type(TDenseDescr), intent(in) :: denseDesc

    !> Electron temperature
    real(dp), intent(in) :: tempElec

    !> Eigenvalue of each level, kpoint and spin channel
    real(dp), intent(in) :: eigvals(:,:,:)

    !> Derivative of rho as a square matrix, if needed
    real(dp), pointer :: dRhoSqr(:,:,:)

  #:if WITH_SCALAPACK
    !> BLACS matrix descriptor
    integer, intent(in) :: desc(DLEN_)
  #:endif

  #:if WITH_SCALAPACK
    integer :: jj, jGlob
  #:else
    integer :: iFilled
  #:endif

    integer :: nSpin, iS
    real(dp) :: workReal(size(eigVecsReal, dim=1), size(eigVecsReal, dim=2))

    iS = parallelKS%localKS(2, iKS)
    nSpin = size(eigVecsReal, dim=3)

    workReal(:,:) = 0.0_dp

  #:if WITH_SCALAPACK

    do jj = 1, size(workReal,dim=2)
      jGlob = scalafx_indxl2g(jj,desc(NB_), env%blacs%orbitalGrid%mycol, desc(CSRC_),&
          & env%blacs%orbitalGrid%ncol)
      if (jGlob >= nEmpty(iS) .and. jGlob <= nFilled(iS)) then
        workReal(:,jj) = eigVecsReal(:,jj,iKS) * &
            & deltamn(eigVals(jGlob, 1, iKS), Ef(iS), tempElec) * dEfdE(iS)
      end if
    end do
    call pblasfx_pgemm(workReal, denseDesc%blacsOrbSqr,eigVecsReal(:,:,iKS),&
        & denseDesc%blacsOrbSqr, workReal, denseDesc%blacsOrbSqr, transb="T",&
        & alpha=real(3-nSpin,dp))

  #:else

    do iFilled = nEmpty(iS), nFilled(iS)
      workReal(:, iFilled) = eigVecsReal(:, iFilled, iS) * &
          & deltamn(eigvals(iFilled, 1, iS), Ef(iS), tempElec) * dEfdE(iS)
    end do
    workReal(:, :) = real(3-nSpin,dp)&
        & * matmul(workReal(:, nEmpty(iS):nFilled(iS)),&
        & transpose(eigVecsReal(:, nEmpty(iS):nFilled(iS), iS)))

  #:endif

    ! pack extra term into density matrix
  #:if WITH_SCALAPACK
    call packRhoRealBlacs(env%blacs, denseDesc, workReal, neighbourList%iNeighbour, nNeighbourSK,&
        & orb%mOrb, iSparseStart, img2CentCell, drhoExtra)
  #:else
    call packHS(drhoExtra, workReal, neighbourList%iNeighbour, nNeighbourSK, orb%mOrb,&
        & denseDesc%iAtomStart, iSparseStart, img2CentCell)
  #:endif

    if  (associated(dRhoSqr)) then
      dRhoSqr(:,:,iS) = dRhoSqr(:,:,iS) + workReal
    end if

  end subroutine dRhoFermiChangeStaticReal


  !> Calculate the derivative of density matrix from derivative of hamiltonian in static case at q=0
  subroutine dRhoStaticPauli(env, dHam, idHam, neighbourList, nNeighbourSK, iSparseStart,&
      & img2CentCell, denseDesc, parallelKS, nFilled, nEmpty, eigVecsCplx, eigVals, Ef, tempElec,&
      & orb, dRhoSparse, idRhoSparse, kPoint, kWeight, iCellVec, cellVec, iKS, transform,&
    #:if WITH_SCALAPACK
      & desc,&
    #:endif
      & dEi, dPsi)

    !> Environment settings
    type(TEnvironment), intent(in) :: env

    !> Derivative of the hamiltonian
    real(dp), intent(in) :: dHam(:,:)

    !> Derivative of the imaginary part of the hamiltonian
    real(dp), intent(inout), allocatable :: idHam(:,:)

    !> list of neighbours for each atom
    type(TNeighbourList), intent(in) :: neighbourList

    !> Number of neighbours for each of the atoms
    integer, intent(in) :: nNeighbourSK(:)

    !> Index array for the start of atomic blocks in sparse arrays
    integer, intent(in) :: iSparseStart(:,:)

    !> map from image atoms to the original unique atom
    integer, intent(in) :: img2CentCell(:)

    !> Dense matrix descriptor
    type(TDenseDescr), intent(in) :: denseDesc

    !> K-points and spins to process
    type(TParallelKS), intent(in) :: parallelKS

    !> ground state eigenvectors
    complex(dp), intent(in) :: eigVecsCplx(:,:,:)

    !> Eigenvalue of each level, kpoint and spin channel
    real(dp), intent(in) :: eigvals(:,:,:)

    !> Fermi level(s)
    real(dp), intent(in) :: Ef(:)

    !> Last (partly) filled level in each spin channel
    integer, intent(in) :: nFilled(:,:)

    !> First (partly) empty level in each spin channel
    integer, intent(in) :: nEmpty(:,:)

    !> Electron temperature
    real(dp), intent(in) :: tempElec

    !> Atomic orbital information
    type(TOrbitals), intent(in) :: orb

    !> returning dRhoSparse on exit
    real(dp), intent(out) :: dRhoSparse(:,:)

    !> returning imaginary part of dRhoSparse on exit
    real(dp), intent(inout), allocatable :: idRhoSparse(:,:)

    !> k-points
    real(dp), intent(in) :: kPoint(:,:)

    !> Weights for k-points
    real(dp), intent(in) :: kWeight(:)

    !> Vectors (in units of the lattice constants) to cells of the lattice
    real(dp), intent(in) :: cellVec(:,:)

    !> Index for which unit cell atoms are associated with
    integer, intent(in) :: iCellVec(:)

    !> spin/kpoint channel
    integer, intent(in) :: iKS

    !> Transformation structure for degenerate orbitals
    type(TRotateDegen), intent(inout) :: transform

  #:if WITH_SCALAPACK
    !> BLACS matrix descriptor
    integer, intent(in) :: desc(DLEN_)
  #:endif

    !> Derivative of single particle eigenvalues
    real(dp), allocatable, intent(inout) :: dEi(:,:,:)

    !> Optional derivatives of single particle wavefunctions
    complex(dp), allocatable, intent(inout) :: dPsi(:,:,:,:)

  #:if WITH_SCALAPACK
    integer :: jj, iGlob, jGlob
  #:else
    integer :: iFilled, iEmpty, nOrb
  #:endif

    integer :: ii, iK, iS
    complex(dp) :: workLocal(size(eigVecsCplx,dim=1), size(eigVecsCplx,dim=2))
    complex(dp) :: dRho(size(eigVecsCplx,dim=1), size(eigVecsCplx,dim=2))
    complex(dp), allocatable :: eigVecsTransformed(:,:)
    logical :: isTransformed

    iK = parallelKS%localKS(1, iKS)
    iS = parallelKS%localKS(2, iKS)

    if (allocated(dEi)) then
      dEi(:, iK, iS) = 0.0_dp
    end if
    if (allocated(dPsi)) then
      dPsi(:, :, iK, iS) = cmplx(0,0,dp)
    end if

    workLocal(:,:) = cmplx(0,0,dp)
    dRho(:,:) = cmplx(0,0,dp)

  #:if WITH_SCALAPACK

    ! dH in square form
    if (allocated(idHam)) then
      call unpackHPauliBlacs(env%blacs, dHam, kPoint(:,iK), neighbourList%iNeighbour,&
          & nNeighbourSK, iCellVec, cellVec, iSparseStart, img2CentCell, orb%mOrb, denseDesc,&
          & workLocal, iorig=idHam)
    else
      call unpackHPauliBlacs(env%blacs, dHam, kPoint(:,iK), neighbourList%iNeighbour,&
          & nNeighbourSK, iCellVec, cellVec, iSparseStart, img2CentCell, orb%mOrb, denseDesc,&
          & workLocal)
    end if

    ! dH times c_i
    call pblasfx_phemm(workLocal, denseDesc%blacsOrbSqr, eigVecsCplx(:,:,iKS),&
        & denseDesc%blacsOrbSqr, dRho, denseDesc%blacsOrbSqr)

    ! c_i times dH times c_i
    call pblasfx_pgemm(eigVecsCplx(:,:,iKS), denseDesc%blacsOrbSqr, dRho,&
        & denseDesc%blacsOrbSqr, workLocal, denseDesc%blacsOrbSqr, transa="C")

    eigvecsTransformed = eigVecsCplx(:,:,iKS)
    call transform%generateUnitary(env, worklocal, eigvals(:,1,iS), eigVecsTransformed, denseDesc,&
        & isTransformed)
    ! now have states orthogonalised agains the operator in degenerate cases, |c~>

    if (isTransformed) then
      ! re-form |c~> H' <c~| with the transformed vectors

      dRho(:,:) = 0.0_dp
      workLocal(:,:) = 0.0_dp

      ! dH in square form
      if (allocated(idHam)) then
        call unpackHPauliBlacs(env%blacs, dHam, kPoint(:,iK), neighbourList%iNeighbour,&
            & nNeighbourSK, iCellVec, cellVec, iSparseStart, img2CentCell, orb%mOrb, denseDesc,&
            & workLocal, iorig=idHam)
      else
        call unpackHPauliBlacs(env%blacs, dHam, kPoint(:,iK), neighbourList%iNeighbour,&
            & nNeighbourSK, iCellVec, cellVec, iSparseStart, img2CentCell, orb%mOrb, denseDesc,&
            & workLocal)
      end if

      ! dH times c_i
      call pblasfx_phemm(workLocal, denseDesc%blacsOrbSqr, eigVecsTransformed,&
          & denseDesc%blacsOrbSqr, dRho, denseDesc%blacsOrbSqr)

      ! c_i times dH times c_i
      call pblasfx_pgemm(eigVecsTransformed, denseDesc%blacsOrbSqr, dRho,&
          & denseDesc%blacsOrbSqr, workLocal, denseDesc%blacsOrbSqr, transa="C")

    end if

    ! derivative of eigenvalues stored in diagonal of matrix workLocal, from <c|h'|c>
    if (allocated(dEi)) then
      do jj = 1, size(workLocal,dim=2)
        jGlob = scalafx_indxl2g(jj, desc(NB_), env%blacs%orbitalGrid%mycol, desc(CSRC_),&
            & env%blacs%orbitalGrid%ncol)
        do ii = 1, size(workLocal,dim=1)
          iGlob = scalafx_indxl2g(ii, desc(MB_), env%blacs%orbitalGrid%myrow, desc(RSRC_),&
              & env%blacs%orbitalGrid%nrow)
          if (iGlob == jGlob) then
            dEi(iGlob, iK, iS) = real(workLocal(ii,jj),dp)
          end if
        end do
      end do
    end if

    ! weight matrix with inverse of energy differences
    do jj = 1, size(workLocal,dim=2)
      jGlob = scalafx_indxl2g(jj, desc(NB_), env%blacs%orbitalGrid%mycol, desc(CSRC_),&
          & env%blacs%orbitalGrid%ncol)
      if (jGlob > nFilled(iS, iK)) then
        workLocal(:, jj) = 0.0_dp
        cycle
      end if
      do ii = 1, size(workLocal,dim=1)
        iGlob = scalafx_indxl2g(ii, desc(MB_), env%blacs%orbitalGrid%myrow, desc(RSRC_),&
            & env%blacs%orbitalGrid%nrow)
        if (iGlob < nEmpty(iS, iK)) then
          workLocal(ii, :) = 0.0_dp
          cycle
        end if
        if (abs(eigvals(jGlob,iK,iS) - eigvals(iGlob,iK,iS)) < epsilon(0.0_dp)&
            & .and. iGlob /= jGlob) then
          ! degenerate, so no contribution
          workLocal(ii,jj) = 0.0_dp
        else
          workLocal(ii, jj) = workLocal(ii, jj) * &
              & invDiff(eigvals(jGlob,iK,iS),eigvals(iGlob,iK,iS),Ef(iS),tempElec)&
              & * theta(eigvals(jGlob,iK,iS),eigvals(iGlob,iK,iS),tempElec)
        end if
      end do
    end do

    ! Derivatives of states
    call pblasfx_pgemm(eigvecsTransformed, denseDesc%blacsOrbSqr, workLocal, denseDesc%blacsOrbSqr,&
        & dRho, denseDesc%blacsOrbSqr)

    if (allocated(dPsi)) then
      dPsi(:, :, iK, iS) = workLocal
    end if

    ! Form derivative of occupied density matrix
    call pblasfx_pgemm(dRho, denseDesc%blacsOrbSqr,eigvecsTransformed, denseDesc%blacsOrbSqr,&
        & workLocal, denseDesc%blacsOrbSqr, transb="C", kk=nFilled(iS, iK))
    dRho(:,:) = workLocal
    ! and hermitize
    call pblasfx_ptranc(workLocal, denseDesc%blacsOrbSqr, dRho, denseDesc%blacsOrbSqr,&
        & beta=(1.0_dp,0.0_dp))

  #:else

    ! serial case
    nOrb = size(dRho, dim = 1)

    call unpackHPauli(dHam, kPoint(:,iK), neighbourList%iNeighbour, nNeighbourSK, iSparseStart,&
        & denseDesc%iAtomStart, img2CentCell, iCellVec, cellVec, dRho, idHam)

    ! form |c> H' <c|
    call hemm(workLocal, 'l', dRho, eigVecsCplx(:,:,iKS))
    workLocal(:,:) = matmul(transpose(conjg(eigVecsCplx(:,:,iKS))), workLocal)

    ! orthogonalise degenerate states against perturbation
    call transform%generateUnitary(workLocal, eigvals(:,iK,iS))
    call transform%degenerateTransform(workLocal)

    ! diagonal elements of workLocal are now derivatives of eigenvalues if needed
    if (allocated(dEi)) then
      do ii = 1, nOrb
        dEi(ii, iK, iS) = real(workLocal(ii,ii),dp)
      end do
    end if

    ! static case

    ! Form actual perturbation U matrix for eigenvectors (potentially at finite T) by
    ! weighting the elements
    do iFilled = 1, nFilled(iS, iK)
      do iEmpty = nEmpty(iS, iK), nOrb
        if (.not.transform%degenerate(iFilled,iEmpty) .or. iEmpty == iFilled) then
          workLocal(iEmpty, iFilled) = workLocal(iEmpty, iFilled) * &
              & invDiff(eigvals(iFilled, iK, 1), eigvals(iEmpty, iK, 1), Ef(1), tempElec)&
              & *theta(eigvals(iFilled, iK, 1), eigvals(iEmpty, iK, 1), tempElec)
        else
          ! rotation should already have set these elements to zero
          workLocal(iEmpty, iFilled) = 0.0_dp
        end if
      end do
    end do

    eigvecsTransformed = eigVecsCplx(:,:,iS)
    call transform%applyUnitary(eigvecsTransformed)

    ! calculate the derivatives of the eigenvectors
    workLocal(:, :nFilled(iS, iK)) =&
        & matmul(eigvecsTransformed(:,nEmpty(iS, iK):), workLocal(nEmpty(iS, iK):,:nFilled(iS, iK)))

    if (allocated(dPsi)) then
      dPsi(:, :, iK, iS) = workLocal
    end if

    ! zero the uncalculated virtual states
    workLocal(:, nFilled(iS, iK)+1:) = 0.0_dp

    ! form the derivative of the density matrix
    dRho(:,:) = matmul(workLocal(:, :nFilled(iS, iK)),&
        & transpose(conjg(eigvecsTransformed(:, :nFilled(iS, iK)))) )&
        & + matmul(eigvecsTransformed(:, :nFilled(iS, iK)),&
        & transpose(conjg(workLocal(:, :nFilled(iS, iK)))) )

  #:endif

    dRhoSparse(:,:) = 0.0_dp
    if (allocated(idRhoSparse)) then
      idRhoSparse(:,:) = 0.0_dp
    end if

  #:if WITH_SCALAPACK
    if (allocated(idRhoSparse)) then
      call packRhoPauliBlacs(env%blacs, denseDesc, dRho, kPoint(:,iK), kWeight(iK),&
          & neighbourList%iNeighbour, nNeighbourSK, orb%mOrb, iCellVec, cellVec, iSparseStart,&
          & img2CentCell, dRhoSparse, idRhoSparse)
      else
        call packRhoPauliBlacs(env%blacs, denseDesc, dRho, kPoint(:,iK), kWeight(iK),&
            & neighbourList%iNeighbour, nNeighbourSK, orb%mOrb, iCellVec, cellVec, iSparseStart,&
            & img2CentCell, dRhoSparse)
      end if
  #:else
      call packHSPauli(dRhoSparse, dRho, neighbourlist%iNeighbour, nNeighbourSK, orb%mOrb,&
          & denseDesc%iAtomStart, iSparseStart, img2CentCell)
      if (allocated(idRhoSparse)) then
        call packHSPauliImag(idRhoSparse, dRho, neighbourlist%iNeighbour, nNeighbourSK,&
            & orb%mOrb, denseDesc%iAtomStart, iSparseStart, img2CentCell)
      end if
  #:endif

      ! adjustment from Pauli to charge/spin
      dRhoSparse(:,:) = 2.0_dp * dRhoSparse

    end subroutine dRhoStaticPauli


  !> Calculate the change in the density matrix due to shift in the Fermi energy
  subroutine dRhoFermiChangeStaticPauli(dRhoExtra, idRhoExtra, env, parallelKS, iKS, kPoint,&
      & kWeight, iCellVec, cellVec, neighbourList, nNEighbourSK, img2CentCell, iSparseStart, dEfdE,&
      & Ef, nFilled, nEmpty, eigVecsCplx, orb, denseDesc, tempElec, eigVals&
    #:if WITH_SCALAPACK
      &, desc&
    #:endif
      &)

    !> Additional contribution to the density matrix to cancel effect of Fermi energy change
    real(dp), intent(inout) :: dRhoExtra(:,:)

    !> Imaginary part of additional contribution to the density matrix to cancel effect of Fermi
    !> energy change
    real(dp), intent(inout), allocatable :: idRhoExtra(:,:)

    !> Environment settings
    type(TEnvironment), intent(inout) :: env

    !> K-points and spins to process
    type(TParallelKS), intent(in) :: parallelKS

    !> spin/kpoint channel
    integer, intent(in) :: iKS

    !> k-points
    real(dp), intent(in) :: kPoint(:,:)

    !> Weights for k-points
    real(dp), intent(in) :: kWeight(:)

    !> Vectors (in units of the lattice constants) to cells of the lattice
    real(dp), intent(in) :: cellVec(:,:)

    !> Index for which unit cell atoms are associated with
    integer, intent(in) :: iCellVec(:)

    !> list of neighbours for each atom
    type(TNeighbourList), intent(in) :: neighbourList

    !> Number of neighbours for each of the atoms
    integer, intent(in) :: nNeighbourSK(:)

    !> map from image atoms to the original unique atom
    integer, intent(in) :: img2CentCell(:)

    !> Index array for the start of atomic blocks in sparse arrays
    integer, intent(in) :: iSparseStart(:,:)

    !> Fermi level derivative
    real(dp), intent(in) :: dEfdE(:)

    !> Fermi level
    real(dp), intent(in) :: Ef(:)

    !> Last (partly) filled level in each spin channel
    integer, intent(in) :: nFilled(:,:)

    !> First (partly) empty level in each spin channel
    integer, intent(in) :: nEmpty(:,:)

    !> ground state eigenvectors
    complex(dp), intent(in) :: eigVecsCplx(:,:,:)

    !> Atomic orbital information
    type(TOrbitals), intent(in) :: orb

    !> Dense matrix descriptor
    type(TDenseDescr), intent(in) :: denseDesc

    !> Electron temperature
    real(dp), intent(in) :: tempElec

    !> Eigenvalue of each level, kpoint and spin channel
    real(dp), intent(in) :: eigvals(:,:,:)

  #:if WITH_SCALAPACK
    !> BLACS matrix descriptor
    integer, intent(in) :: desc(DLEN_)
  #:endif

  #:if WITH_SCALAPACK
    integer :: jj, jGlob
  #:else
    integer :: iFilled
  #:endif
    integer :: iK, iS
    complex(dp) :: workLocal(size(eigVecsCplx, dim=1), size(eigVecsCplx, dim=2))
  #:if WITH_SCALAPACK
    complex(dp) :: workLocal2(size(eigVecsCplx, dim=1), size(eigVecsCplx, dim=2))
    complex(dp) :: workLocal3(size(eigVecsCplx, dim=1), size(eigVecsCplx, dim=2))
  #:endif

    workLocal(:,:) = cmplx(0,0,dp)

    iK = parallelKS%localKS(1, iKS)
    iS = parallelKS%localKS(2, iKS)

  #:if WITH_SCALAPACK

    do jj = 1, size(workLocal,dim=2)
      jGlob = scalafx_indxl2g(jj, desc(NB_), env%blacs%orbitalGrid%mycol, desc(CSRC_),&
          & env%blacs%orbitalGrid%ncol)
      if (jGlob >= nEmpty(iS,iK) .and. jGlob <= nFilled(iS,iK)) then
        workLocal(:, jj) = eigVecsCplx(:, jj, iKS) * &
            & deltamn(eigVals(jGlob, iK, iS), Ef(iS), tempElec) * dEfdE(iS)
      end if
    end do

    workLocal3(:,:) = eigVecsCplx(:,:,iKS)
    call pblasfx_pgemm(workLocal, denseDesc%blacsOrbSqr, workLocal3,&
        & denseDesc%blacsOrbSqr, workLocal2, denseDesc%blacsOrbSqr, transb="C")
    workLocal(:,:) = workLocal2
    ! Hermitian symmetry
    call pblasfx_ptranc(workLocal2, denseDesc%blacsOrbSqr, workLocal, denseDesc%blacsOrbSqr,&
        & alpha=(0.5_dp,0.0_dp), beta=(0.5_dp,0.0_dp))
  #:else

    do iFilled = nEmpty(iS,iK), nFilled(iS,iK)
      workLocal(:, iFilled) = eigVecsCplx(:, iFilled, 1) * &
          & deltamn(eigvals(iFilled, 1, 1), Ef(1), tempElec) * dEfdE(1)
    end do

    workLocal(:, :) = matmul(workLocal(:, nEmpty(iS,iK):nFilled(iS,iK)),&
        & transpose(conjg(eigVecsCplx(:, nEmpty(iS,iK):nFilled(iS,iK), 1))))

  #:endif

    ! pack extra term into density matrix
  #:if WITH_SCALAPACK
    if (allocated(idRhoExtra)) then
      call packRhoPauliBlacs(env%blacs, denseDesc, workLocal, kPoint(:,iK), kWeight(iK),&
          & neighbourList%iNeighbour, nNeighbourSK, orb%mOrb, iCellVec, cellVec, iSparseStart,&
          & img2CentCell, dRhoExtra, idRhoExtra)
    else
      call packRhoPauliBlacs(env%blacs, denseDesc, workLocal, kPoint(:,iK), kWeight(iK),&
          & neighbourList%iNeighbour, nNeighbourSK, orb%mOrb, iCellVec, cellVec, iSparseStart,&
          & img2CentCell, dRhoExtra)
    end if
  #:else
    call packHSPauli(dRhoExtra, workLocal, neighbourlist%iNeighbour, nNeighbourSK, orb%mOrb,&
        & denseDesc%iAtomStart, iSparseStart, img2CentCell)
    if (allocated(idRhoExtra)) then
      call packHSPauliImag(idRhoExtra, workLocal, neighbourlist%iNeighbour, nNeighbourSK,&
          & orb%mOrb, denseDesc%iAtomStart, iSparseStart, img2CentCell)
    end if

  #:endif

  end subroutine dRhoFermiChangeStaticPauli

end module dftbp_linearresponse
