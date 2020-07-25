!--------------------------------------------------------------------------------------------------!
!  DFTB+: general package for performing fast atomistic simulations                                !
!  Copyright (C) 2006 - 2020  DFTB+ developers group                                               !
!                                                                                                  !
!  See the LICENSE file for terms of usage and distribution.                                       !
!--------------------------------------------------------------------------------------------------!

#:include 'common.fypp'

!> REKS and SI-SA-REKS formulation in DFTB as developed by Lee et al.
!>
!> The functionality of the module has some limitation:
!> * Third order does not work.
!> * Periodic system do not work yet appart from Gamma point.
!> * Orbital potentials or spin-orbit or external E-field does not work yet.
!> * Only for closed shell system.
!> * Onsite corrections are not included in this version
module dftbp_reksgrad

#:if WITH_OMP
  use omp_lib
#:endif
  use dftbp_accuracy
  use dftbp_blasroutines, only : gemm, gemv
  use dftbp_coulomb, only : addInvRPrime
  use dftbp_densedescr
  use dftbp_environment
  use dftbp_globalenv
  use dftbp_lapackroutines, only : getrf, getri
  use dftbp_message
  use dftbp_nonscc
  use dftbp_orbitals
  use dftbp_periodic
  use dftbp_rangeseparated
  use dftbp_scc
  use dftbp_schedule
  use dftbp_slakocont
  use dftbp_sparse2dense
  use dftbp_rekscommon
  use dftbp_reksvar, only : reksTypes

  implicit none

  private

  public :: getEnergyWeightedDensityL
  public :: derivative_blockL, weightGradient

  public :: getSccSpinLrPars, getHxcKernel, getG1ILOmegaRab, getSuperAMatrix
  public :: buildSaReksVectors, buildInteractionVectors, buildLstateVector, solveZT
  public :: getRmat, getRdel, getZmat, getQ1mat, getQ1del, getQ2mat

  public :: SaToSsrXT, SaToSsrWeight, SaToSsrGradient, addSItoRQ
  public :: SSRshift, SIshift, Lshift, RTshift
  public :: getOtherSAgrad, getReksNAC

  public :: getExtChrgGradients

  contains

  !> Calculate energy weighted density matrix for each microstate
  subroutine getEnergyWeightedDensityL(env, denseDesc, neighbourList, &
      & nNeighbourSK, iSparseStart, img2CentCell, orb, hamSqrL, hamSpL, &
      & fillingL, eigenvecs, Lpaired, Efunc, isRangeSep, edmSpL)

    !> Environment settings
    type(TEnvironment), intent(inout) :: env

    !> Dense matrix descriptor
    type(TDenseDescr), intent(in) :: denseDesc

    !> neighbours to atoms
    type(TNeighbourList), intent(in) :: neighbourList

    !> Number of atomic neighbours
    integer, intent(in) :: nNeighbourSK(:)

    !> Index for atomic blocks in sparse data
    integer, intent(in) :: iSparseStart(:,:)

    !> map from image atom to real atoms
    integer, intent(in) :: img2CentCell(:)

    !> Atomic orbital information
    type(TOrbitals), intent(in) :: orb

    !> Dense Hamiltonian matrix for each microstate
    real(dp), allocatable, intent(inout) :: hamSqrL(:,:,:,:)

    !> Sparse Hamiltonian matrix for each microstate
    real(dp), allocatable, intent(in) :: hamSpL(:,:,:)

    !> Filling for each microstate
    real(dp), intent(in) :: fillingL(:,:,:)

    !> Eigenvectors on eixt
    real(dp), intent(in) :: eigenvecs(:,:)

    !> Number of spin-paired microstates
    integer, intent(in) :: Lpaired

    !> Minimized energy functional
    integer, intent(in) :: Efunc

    !> Whether to run a range separated calculation
    logical, intent(in) :: isRangeSep

    !> sparse energy-weighted density matrix for each microstate
    real(dp), intent(out) :: edmSpL(:,:)

    real(dp), allocatable :: tmpEps(:,:)
    real(dp), allocatable :: tmpHam(:,:)

    integer :: tmpL, tmpchk, iL, Lmax, i, j, nOrb

    Lmax = size(fillingL,dim=3)
    nOrb = size(fillingL,dim=1)

    allocate(tmpEps(nOrb,nOrb))
    if (.not. isRangeSep) then
      allocate(tmpHam(nOrb,nOrb))
    end if

    edmSpL(:,:) = 0.0_dp
    do iL = 1, Lmax

      if (isRangeSep) then
        if (Efunc == 1) then
          ! For single-state REKS, current hamSqrL is still in AO basis
          ! since the secular equation routine is not used
          ! convert the hamiltonians from AO basis to MO basis
          call matAO2MO(hamSqrL(:,:,1,iL), eigenvecs)
        end if
      else
        tmpHam(:,:) = 0.0_dp
        ! hamSpL has (my_ud) component
        call env%globalTimer%startTimer(globalTimers%sparseToDense)
        call unpackHS(tmpHam, hamSpL(:,1,iL), &
            & neighbourList%iNeighbour, nNeighbourSK, &
            & denseDesc%iAtomStart, iSparseStart, img2CentCell)
        call env%globalTimer%stopTimer(globalTimers%sparseToDense)
        call blockSymmetrizeHS(tmpHam, denseDesc%iAtomStart)
        ! convert the hamiltonians from AO basis to MO basis
        call matAO2MO(tmpHam, eigenvecs)
      end if

      ! calculate the lagrange multipliers
      tmpEps(:,:) = 0.0_dp
      do i = 1, nOrb
        do j = 1, nOrb
          if (isRangeSep) then
            tmpEps(i,j) = (fillingL(i,1,iL) + fillingL(j,1,iL)) &
                & * hamSqrL(i,j,1,iL) * 0.5_dp
          else
            tmpEps(i,j) = (fillingL(i,1,iL) + fillingL(j,1,iL)) &
                & * tmpHam(i,j) * 0.5_dp
          end if
        end do
      end do

      if (iL <= Lpaired) then
        tmpL = iL
        tmpchk = 0
        tmpEps(:,:) = 2.0_dp * tmpEps
      else
        if (mod(iL,2) == 1) then
          tmpL = iL
          tmpchk = 0
        else
          tmpL = iL - 1
          tmpchk = 1
        end if
      end if

      ! convert the multipliers from MO basis to AO basis
      call matMO2AO(tmpEps, eigenvecs)

      call env%globalTimer%startTimer(globalTimers%denseToSparse)
      call packHS(edmSpL(:,tmpL), tmpEps, neighbourlist%iNeighbour, &
          & nNeighbourSK, orb%mOrb, denseDesc%iAtomStart, iSparseStart, &
          & img2CentCell)
      call env%globalTimer%stopTimer(globalTimers%denseToSparse)

      if (tmpchk == 1) then
        edmSpL(:,iL) = edmSpL(:,tmpL)
      end if

    end do

  end subroutine getEnergyWeightedDensityL


  !> The SCC and spin electronic force contribution for all atoms from the matrix derivatives, self
  !> consistent potential and the density and energy-density matrices
  subroutine derivative_blockL(env, deriv, derivator, DM, EDM, skHamCont, &
      & skOverCont, coords, species, iNeighbour, nNeighbourSK, img2CentCell, &
      & iPair, iSquare, orb, shift, Lpaired, gradL, Hderiv, Sderiv)

    !> Computational environment settings
    type(TEnvironment), intent(in) :: env

    !> x,y,z derivatives for each real atom in the system
    real(dp), intent(out) :: deriv(:,:)

    !> Differentiatior for the non-scc components
    class(TNonSccDiff), intent(in) :: derivator

    !> density matrix in unpacked format
    real(dp), intent(in) :: DM(:,:,:,:)

    !> energy-weighted density matrix in packed format
    real(dp), intent(in) :: EDM(:,:)

    !> Container for SK Hamiltonian integrals
    type(TSlakoCont) :: skHamCont

    !> Container for SK overlap integrals
    type(TSlakoCont) :: skOverCont

    !> list of all atomic coordinates
    real(dp), intent(in) :: coords(:,:)

    !> list of all atomic species
    integer, intent(in) :: species(:)

    !> neighbour list for atoms
    integer, intent(in) :: iNeighbour(0:,:)

    !> number of neighbours of each atom
    integer, intent(in) :: nNeighbourSK(:)

    !> indexing array for periodic image atoms
    integer, intent(in) :: img2CentCell(:)

    !> indexing array for the Hamiltonian
    integer, intent(in) :: iPair(0:,:)

    !> Position of each atom in the rows/columns of the square matrices. Shape: (nAtom)
    integer, intent(in) :: iSquare(:)

    !> Information about the shells and orbitals in the system.
    type(TOrbitals), intent(in) :: orb

    !> block shift from the potential
    real(dp), intent(in) :: shift(:,:,:,:,:)

    !> Number of spin-paired microstates
    integer, intent(in) :: Lpaired

    !> gradients for each microstate except orbital derivative terms
    real(dp), intent(out) :: gradL(:,:,:)

    !> Dense non-scc Hamiltonian derivative in AO basis
    real(dp), intent(out) :: Hderiv(:,:,:)

    !> Dense overlap derivative in AO basis
    real(dp), intent(out) :: Sderiv(:,:,:)

    integer :: iOrig, iSpin, ii, nSpin, nAtom
    integer :: iNeigh, iAtom1, iAtom2, iAtom2f, iSp1, iSp2
    integer :: nOrb1, nOrb2

    real(dp) :: sqrDMTmp(orb%mOrb,orb%mOrb), sqrEDMTmp(orb%mOrb,orb%mOrb)
    real(dp) :: shiftSprime(orb%mOrb,orb%mOrb)
    real(dp) :: hPrimeTmp(orb%mOrb,orb%mOrb,3), sPrimeTmp(orb%mOrb,orb%mOrb,3)
    real(dp) :: derivTmp(3)
    real(dp) :: fac
    integer :: iAtFirst, iAtLast
    integer :: tmpL, tmpLs, iL, Lmax

    nAtom = size(orb%nOrbAtom)
    nSpin = size(shift,dim=4)
    Lmax = size(DM,dim=4)

    do iL = 1, Lmax

      deriv(:,:) = 0.0_dp
      if (iL == 1) then
        Hderiv(:,:,:) = 0.0_dp
        Sderiv(:,:,:) = 0.0_dp
      end if

      if (iL <= Lpaired) then
        tmpL = iL
      else
        if (mod(iL,2) == 1) then
          tmpL = iL
          tmpLs = iL + 1
          fac = 1.0_dp
        else
          tmpL = iL - 1
          tmpLs = iL
          fac = -1.0_dp
        end if
      end if

      call distributeRangeInChunks(env, 1, nAtom, iAtFirst, iAtLast)

      !$OMP PARALLEL DO PRIVATE(iAtom1,iSp1,nOrb1,iNeigh,iAtom2,iAtom2f,iSp2,nOrb2,iOrig,sqrDMTmp, &
      !$OMP& sqrEDMTmp,hPrimeTmp,sPrimeTmp,derivTmp,shiftSprime,iSpin,ii) DEFAULT(SHARED) &
      !$OMP& SCHEDULE(RUNTIME) REDUCTION(+:deriv)
      do iAtom1 = iAtFirst, iAtLast
        iSp1 = species(iAtom1)
        nOrb1 = orb%nOrbSpecies(iSp1)
        do iNeigh = 1, nNeighbourSK(iAtom1)
          iAtom2 = iNeighbour(iNeigh, iAtom1)
          iAtom2f = img2CentCell(iAtom2)
          iSp2 = species(iAtom2f)
          if (iAtom1 /= iAtom2f) then
            nOrb2 = orb%nOrbSpecies(iSp2)
            iOrig = iPair(iNeigh,iAtom1) + 1
            sqrDMTmp(1:nOrb2,1:nOrb1) = &
                & DM(iSquare(iAtom2f):iSquare(iAtom2f+1)-1, &
                & iSquare(iAtom1):iSquare(iAtom1+1)-1,1,tmpL)
            sqrEDMTmp(1:nOrb2,1:nOrb1) = &
                & reshape(EDM(iOrig:iOrig+nOrb1*nOrb2-1,iL),(/nOrb2,nOrb1/))
            call derivator%getFirstDeriv(hPrimeTmp, skHamCont, &
                & coords, species, iAtom1, iAtom2, orb)
            call derivator%getFirstDeriv(sPrimeTmp, skOverCont, &
                & coords, species, iAtom1, iAtom2, orb)

            if (iL == 1) then
              ! H0 & S derivative with respect to AO basis
              do ii = 1, 3
                Hderiv(iSquare(iAtom2f):iSquare(iAtom2f+1)-1, &
                    & iSquare(iAtom1):iSquare(iAtom1+1)-1,ii) &
                    & = hPrimeTmp(1:nOrb2,1:nOrb1,ii)
                Sderiv(iSquare(iAtom2f):iSquare(iAtom2f+1)-1, &
                    & iSquare(iAtom1):iSquare(iAtom1+1)-1,ii) &
                    & = sPrimeTmp(1:nOrb2,1:nOrb1,ii)
              end do
            end if
            
            derivTmp(:) = 0.0_dp
            ! note factor of 2 for implicit summation over lower triangle of density matrix:
            do ii = 1, 3
              derivTmp(ii) = 2.0_dp * (&
                  & sum(sqrDMTmp(1:nOrb2,1:nOrb1)*hPrimeTmp(1:nOrb2,1:nOrb1,ii)) &
                  & - sum(sqrEDMTmp(1:nOrb2,1:nOrb1)*sPrimeTmp(1:nOrb2,1:nOrb1,ii)))
            end do

            do iSpin = 1, nSpin
              do ii = 1, 3
                shiftSprime(1:nOrb2,1:nOrb1) = 0.5_dp * (&
                    & matmul(sPrimeTmp(1:nOrb2,1:nOrb1,ii), &
                    & shift(1:nOrb1,1:nOrb1,iAtom1,iSpin,iL) )&
                    & + matmul(shift(1:nOrb2,1:nOrb2,iAtom2f,iSpin,iL), &
                    & sPrimeTmp(1:nOrb2,1:nOrb1,ii)) )
                ! again factor of 2 from lower triangle, cf published force expressions for SCC:
                if (iSpin == 1) then
                  derivTmp(ii) = derivTmp(ii) + 2.0_dp * ( sum( shiftSprime(1:nOrb2,1:nOrb1) *&
                      & sqrDMTmp(1:nOrb2,1:nOrb1) ) )
                else
                  if (iL > Lpaired) then
                    derivTmp(ii) = derivTmp(ii)&
                        & + fac * 2.0_dp * ( sum( shiftSprime(1:nOrb2,1:nOrb1) *&
                        & DM(iSquare(iAtom2f):iSquare(iAtom2f+1)-1, &
                        & iSquare(iAtom1):iSquare(iAtom1+1)-1,1,tmpLs) ) )
                  end if
                end if
              end do
            end do

            ! forces from atom 1 on atom 2f and 2f onto 1
            deriv(:,iAtom1) = deriv(:,iAtom1) + derivTmp(:)
            deriv(:,iAtom2f) = deriv(:,iAtom2f) - derivTmp(:)

          end if
        enddo
      enddo
      !$OMP END PARALLEL DO

      if (iL == 1) then
        do ii = 1, 3
          call blockSymmetrizeHS(Hderiv(:,:,ii), iSquare)
          call blockSymmetrizeHS(Sderiv(:,:,ii), iSquare)
        end do
      end if

      call assembleChunks(env, deriv)

      gradL(:,:,iL) = deriv

    end do

  end subroutine derivative_blockL


  !> Calculate averaged gradient with weighting factors of each microstate
  subroutine weightGradient(gradL, weight, grad)

    !> gradients for each microstate except orbital derivative terms
    real(dp), intent(in) :: gradL(:,:,:)

    !> Weight of each microstate for state to be optimized; weight = weightL * SAweight
    real(dp), intent(in) :: weight(:)

    !> averaged gradient with weighting factors
    real(dp), intent(out) :: grad(:,:)

    integer :: iL, Lmax

    Lmax = size(weight,dim=1)

    grad(:,:) = 0.0_dp
    do iL = 1, Lmax
      grad(:,:) = grad + weight(iL) * gradL(:,:,iL)
    end do

  end subroutine weightGradient


  !> Calculate SCC, spin, LC parameters with matrix form
  subroutine getSccSpinLrPars(env, sccCalc, rangeSep, coords, species, &
      & iNeighbour, img2CentCell, iSquare, spinW, getAtomIndex, isRangeSep, &
      & GammaAO, GammaDeriv, SpinAO, LrGammaAO, LrGammaDeriv)

    !> Computational environment settings
    type(TEnvironment), intent(in) :: env

    !> SCC module internal variables
    type(TScc), allocatable, intent(inout) :: sccCalc

    !> Range separation contributions
    type(TRangeSepFunc), allocatable, intent(inout) :: rangeSep

    !> list of all atomic coordinates
    real(dp), intent(in) :: coords(:,:)

    !> list of all atomic species
    integer, intent(in) :: species(:)

    !> neighbour list for atoms
    integer, intent(in) :: iNeighbour(0:,:)

    !> indexing array for periodic image atoms
    integer, intent(in) :: img2CentCell(:)

    !> Position of each atom in the rows/columns of the square matrices. Shape: (nAtom)
    integer, intent(in) :: iSquare(:)

    !> spin constants
    real(dp), intent(in) :: spinW(:,:,:)

    !> get atom index from AO index
    integer, intent(in) :: getAtomIndex(:)

    !> Whether to run a range separated calculation
    logical, intent(in) :: isRangeSep

    !> scc gamma integrals in AO basis
    real(dp), intent(out) :: GammaAO(:,:)

    !> scc gamma derivative integrals
    real(dp), intent(out) :: GammaDeriv(:,:,:)

    !> spin W in AO basis
    real(dp), intent(out) :: SpinAO(:,:)

    !> long-range gamma integrals in AO basis
    real(dp), allocatable, intent(inout) :: LrGammaAO(:,:)

    !> long-range gamma derivative integrals
    real(dp), allocatable, intent(inout) :: LrGammaDeriv(:,:,:)

    real(dp), allocatable :: tmpGamma(:,:)
    real(dp), allocatable :: tmpLrGamma(:,:)

    real(dp) :: fac1, fac2
    integer :: nOrb, nAtom, mu, nu, G1, G2, ii, iAt1, iAt2, iSp1, iSp2

    nOrb = size(GammaAO,dim=1)
    nAtom = size(iSquare,dim=1) - 1

    allocate(tmpGamma(nAtom,nAtom))
    if (isRangeSep) then
      allocate(tmpLrGamma(nAtom,nAtom))
    end if

    ! get total gamma (gamma = 1/R - S)
    tmpGamma(:,:) = 0.0_dp
    call sccCalc%getAtomicGammaMatrix(tmpGamma, iNeighbour, img2CentCell)
    call symmetrizeHS(tmpGamma)
    ! convert from atom to AO
    GammaAO(:,:) = 0.0_dp
    do iAt1 = 1, nAtom
      do iAt2 = 1, nAtom
        GammaAO(iSquare(iAt2):iSquare(iAt2+1)-1,iSquare(iAt1):iSquare(iAt1+1)-1) =&
            & tmpGamma(iAt2,iAt1)
      end do
    end do

    ! get total gamma derivative (gamma = 1/R - S)
    GammaDeriv(:,:,:) = 0.0_dp
    call sccCalc%getGammaDeriv(env, species, iNeighbour, img2CentCell, GammaDeriv)
    do ii = 1, 3
      call symmetrizeHS(GammaDeriv(:,:,ii))
    end do

    ! get spinW with respect to AO
    SpinAO(:,:) = 0.0_dp
    do mu = 1, nOrb
      do nu = mu, nOrb
        ! find proper atom index
        G1 = getAtomIndex(mu)
        G2 = getAtomIndex(nu)
        if (G1 == G2) then
          ! set the orbital shell for mu index
          call findShellOfAO(mu, mu, getAtomIndex, iSquare, iSp1, fac1, fac2)
          ! set the orbital shell for nu index
          call findShellOfAO(nu, nu, getAtomIndex, iSquare, iSp2, fac1, fac2)
          SpinAO(nu,mu) = spinW(iSp2,iSp1,species(G1))
          if (nu /= mu) then
            SpinAO(mu,nu) = SpinAO(nu,mu)
          end if
        end if
      end do
    end do

    if (isRangeSep) then

      ! get total long-range gamma
      tmpLrGamma(:,:) = 0.0_dp
      call rangeSep%getLrGamma(tmpLrGamma)
      ! convert from atom to AO
      LrGammaAO(:,:) = 0.0_dp
      do iAt1 = 1, nAtom
        do iAt2 = 1, nAtom
          LrGammaAO(iSquare(iAt2):iSquare(iAt2+1)-1,iSquare(iAt1):iSquare(iAt1+1)-1) =&
              & tmpLrGamma(iAt2,iAt1)
        end do
      end do

      ! get long-range gamma derivative
      LrGammaDeriv(:,:,:) = 0.0_dp
      call rangeSep%getLrGammaDeriv(coords, species, LrGammaDeriv)
      do ii = 1, 3
        call symmetrizeHS(LrGammaDeriv(:,:,ii))
      end do

    end if

  end subroutine getSccSpinLrPars


  !> Interface routine to calculate H-XC kernel in REKS
  subroutine getHxcKernel(getDenseAO, over, overSqr, GammaAO, SpinAO, LrGammaAO, Glevel, tSaveMem,&
      & isRangeSep, HxcSpS, HxcSpD, HxcHalfS, HxcHalfD, HxcSqrS, HxcSqrD)

    !> get dense AO index from sparse AO array
    integer, intent(in) :: getDenseAO(:,:)

    !> sparse overlap matrix
    real(dp), intent(in) :: over(:)

    !> Dense overlap matrix
    real(dp), intent(in) :: overSqr(:,:)

    !> scc gamma integrals in AO basis
    real(dp), intent(in) :: GammaAO(:,:)

    !> spin W in AO basis
    real(dp), intent(in) :: SpinAO(:,:)

    !> long-range gamma integrals in AO basis
    real(dp), allocatable, intent(in) :: LrGammaAO(:,:)

    !> Algorithms to calculate analytic gradients
    integer, intent(in) :: Glevel

    !> Save 'A' and 'Hxc' to memory in gradient calculation
    logical, intent(in) :: tSaveMem

    !> Whether to run a range separated calculation
    logical, intent(in) :: isRangeSep

    !> Hartree-XC kernel with sparse form with same spin part
    real(dp), allocatable, intent(inout) :: HxcSpS(:,:)

    !> Hartree-XC kernel with sparse form with different spin part
    real(dp), allocatable, intent(inout) :: HxcSpD(:,:)

    !> Hartree-XC kernel with half dense form with same spin part
    real(dp), allocatable, intent(inout) :: HxcHalfS(:,:)

    !> Hartree-XC kernel with half dense form with different spin part
    real(dp), allocatable, intent(inout) :: HxcHalfD(:,:)

    !> Hartree-XC kernel with dense form with same spin part
    real(dp), allocatable, intent(inout) :: HxcSqrS(:,:,:,:)

    !> Hartree-XC kernel with dense form with different spin part
    real(dp), allocatable, intent(inout) :: HxcSqrD(:,:,:,:)

    if (Glevel == 1 .or. Glevel == 2) then

      if (tSaveMem) then

        if (isRangeSep) then

          ! get Hxc kernel for DFTB with respect to AO basis
          ! for LC case, we use half dense form.
          call HxcKernelHalf_(getDenseAO, overSqr, GammaAO, SpinAO, LrGammaAO, isRangeSep,&
              & HxcHalfS, HxcHalfD)

        else

          ! get Hxc kernel for DFTB with respect to AO basis with sparse form
          call HxcKernelSparse_(getDenseAO, over, GammaAO, SpinAO, &
              & HxcSpS, HxcSpD)

        end if

      end if

    else if (Glevel == 3) then

      ! get Hxc kernel for DFTB with respect to AO basis
      call HxcKernelDense_(overSqr, GammaAO, SpinAO, LrGammaAO, isRangeSep, HxcSqrS, HxcSqrD)

    end if

  end subroutine getHxcKernel


  !> Calculate G1, weightIL, omega, Rab variables
  subroutine getG1ILOmegaRab(env, denseDesc, neighbourList, &
      & nNeighbourSK, iSparseStart, img2CentCell, eigenvecs, hamSqrL, &
      & hamSpL, fockFa, fillingL, FONs, SAweight, enLtot, hess, &
      & Nc, Na, reksAlg, tSSR, isRangeSep, G1, weightIL, omega, Rab)

    !> Environment settings
    type(TEnvironment), intent(inout) :: env

    !> Dense matrix descriptor
    type(TDenseDescr), intent(in) :: denseDesc

    !> neighbours to atoms
    type(TNeighbourList), intent(in) :: neighbourList

    !> Number of atomic neighbours
    integer, intent(in) :: nNeighbourSK(:)

    !> Index for atomic blocks in sparse data
    integer, intent(in) :: iSparseStart(:,:)

    !> map from image atom to real atoms
    integer, intent(in) :: img2CentCell(:)

    !> Eigenvectors on eixt
    real(dp), intent(in) :: eigenvecs(:,:,:)

    !> Dense Hamiltonian matrix for each microstate
    real(dp), allocatable, intent(in) :: hamSqrL(:,:,:,:)

    !> Sparse Hamiltonian matrix for each microstate
    real(dp), allocatable, intent(in) :: hamSpL(:,:,:)

    !> dense fock matrix for active orbitals
    real(dp), intent(in) :: fockFa(:,:,:)

    !> Filling for each microstate
    real(dp), intent(in) :: fillingL(:,:,:)

    !> Fractional occupation numbers of active orbitals
    real(dp), intent(in) :: FONs(:,:)

    !> Weights used in state-averaging
    real(dp), intent(in) :: SAweight(:)

    !> total energy for each microstate
    real(dp), intent(in) :: enLtot(:)

    !> Hessian of FONs
    real(dp), intent(in) :: hess

    !> Number of core orbitals
    integer, intent(in) :: Nc

    !> Number of active orbitals
    integer, intent(in) :: Na

    !> Type of REKS calculations
    integer, intent(in) :: reksAlg

    !> Calculate SSR state with inclusion of SI, otherwise calculate SA-REKS state
    logical, intent(in) :: tSSR

    !> Whether to run a range separated calculation
    logical, intent(in) :: isRangeSep

    !> constant calculated from hessian and energy of microstates
    real(dp), intent(out) :: G1

    !> modified weight of each microstate
    real(dp), intent(out) :: weightIL(:)

    !> anti-symmetric matrices originated from Hamiltonians
    real(dp), intent(out) :: omega(:)

    !> state-interaction term used in SSR gradients
    real(dp), allocatable, intent(inout) :: Rab(:,:)

    real(dp), allocatable :: tmpHam(:,:)

    real(dp) :: tmpdE
    integer :: nOrb, Lmax, superN
    integer :: iL, ij, i, j, Nv

    nOrb = size(fillingL,dim=1)
    superN = size(omega,dim=1)
    Lmax = size(fillingL,dim=3)
    Nv = nOrb - Nc - Na

    if (.not. isRangeSep) then
      allocate(tmpHam(nOrb,nOrb))
    end if

    ! set G1 value
    select case (reksAlg)
    case (reksTypes%noReks)
    case (reksTypes%ssr22)
      tmpdE = enLtot(5) + enLtot(6) &
           & - enLtot(4) - enLtot(3)
      G1 = 2.0_dp / hess / tmpdE
    case (reksTypes%ssr44)
      call error("SSR(4,4) is not implemented yet")
    end select

    ! set weightIL value
    select case (reksAlg)
    case (reksTypes%noReks)
    case (reksTypes%ssr22)
      weightIL(1) =  1.0_dp
      weightIL(2) = -1.0_dp
      weightIL(3) = (enLtot(1) - enLtot(2)) / tmpdE
      weightIL(4) =  weightIL(3)
      weightIL(5) = -weightIL(3)
      weightIL(6) = -weightIL(3)
    case (reksTypes%ssr44)
      call error("SSR(4,4) is not implemented yet")
    end select

    if (tSSR) then
      Rab(:,:) = 0.0_dp
    end if
    omega(:) = 0.0_dp
    do iL = 1, Lmax

      if (isRangeSep) then

        ! set omega value
        do ij = 1, superN
          ! assign index i and j from ij
          call assignIndex(Nc, Na, Nv, reksAlg, ij, i, j)
          omega(ij) = omega(ij) + 2.0_dp * weightIL(iL) * hamSqrL(i,j,1,iL) &
              & * ( fillingL(i,1,iL) - fillingL(j,1,iL) )
        end do

        ! get Rab value
        if (tSSR) then
          select case (reksAlg)
          case (reksTypes%noReks)
          case (reksTypes%ssr22)
            call getRab22_1st_(hamSqrL(:,:,1,iL), fillingL, weightIL, &
                & SAweight, FONs, Nc, iL, Rab)
          case (reksTypes%ssr44)
            call error("SSR(4,4) is not implemented yet")
          end select
        end if

      else

        ! convert from sparse to dense for hamSpL in MO basis
        tmpHam(:,:) = 0.0_dp
        ! hamSpL has (my_ud) component
        call env%globalTimer%startTimer(globalTimers%sparseToDense)
        call unpackHS(tmpHam, hamSpL(:,1,iL), &
            & neighbourList%iNeighbour, nNeighbourSK, &
            & denseDesc%iAtomStart, iSparseStart, img2CentCell)
        call env%globalTimer%stopTimer(globalTimers%sparseToDense)
        call blockSymmetrizeHS(tmpHam, denseDesc%iAtomStart)
        ! convert the multipliers from MO basis to AO basis
        call matAO2MO(tmpHam, eigenvecs(:,:,1))

        ! set omega value
        do ij = 1, superN
          ! assign index i and j from ij
          call assignIndex(Nc, Na, Nv, reksAlg, ij, i, j)
          omega(ij) = omega(ij) + 2.0_dp * weightIL(iL) * tmpHam(i,j) &
              & * ( fillingL(i,1,iL) - fillingL(j,1,iL) )
        end do

        ! get Rab value
        if (tSSR) then
          select case (reksAlg)
          case (reksTypes%noReks)
          case (reksTypes%ssr22)
            call getRab22_1st_(tmpHam, fillingL, weightIL, &
                & SAweight, FONs, Nc, iL, Rab)
          case (reksTypes%ssr44)
            call error("SSR(4,4) is not implemented yet")
          end select
        end if

      end if

    end do

    ! get Rab value
    if (tSSR) then
      select case (reksAlg)
      case (reksTypes%noReks)
      case (reksTypes%ssr22)
        call getRab22_2nd_(fockFa, FONs, SAweight, Nc, Rab)
      case (reksTypes%ssr44)
        call error("SSR(4,4) is not implemented yet")
      end select
    end if

  end subroutine getG1ILOmegaRab


  !> Calculate super A hessian matrix with and without H-XC kernel
  subroutine getSuperAMatrix(eigenvecs, HxcSqrS, HxcSqrD, Fc, Fa, &
      & omega, fillingL, weight, SAweight, FONs, G1, Lpaired, Nc, &
      & Na, Glevel, reksAlg, tSaveMem, A1e, A1ePre, Aall)

    !> Eigenvectors on eixt
    real(dp), intent(in) :: eigenvecs(:,:,:)

    !> Hartree-XC kernel with dense form with same spin part
    real(dp), allocatable, intent(in) :: HxcSqrS(:,:,:,:)

    !> Hartree-XC kernel with dense form with different spin part
    real(dp), allocatable, intent(in) :: HxcSqrD(:,:,:,:)

    !> dense fock matrix for core orbitals
    real(dp), intent(in) :: Fc(:,:)

    !> dense fock matrix for active orbitals
    real(dp), intent(in) :: Fa(:,:,:)

    !> anti-symmetric matrices originated from Hamiltonians
    real(dp), intent(in) :: omega(:)

    !> Filling for each microstate
    real(dp), intent(in) :: fillingL(:,:,:)

    !> Weight of each microstate for state to be optimized; weight = weightL * SAweight
    real(dp), intent(in) :: weight(:)

    !> Weights used in state-averaging
    real(dp), intent(in) :: SAweight(:)

    !> Fractional occupation numbers of active orbitals
    real(dp), intent(in) :: FONs(:,:)

    !> constant calculated from hessian and energy of microstates
    real(dp), intent(in) :: G1

    !> Number of spin-paired microstates
    integer, intent(in) :: Lpaired

    !> Number of core orbitals
    integer, intent(in) :: Nc

    !> Number of active orbitals
    integer, intent(in) :: Na

    !> Algorithms to calculate analytic gradients
    integer, intent(in) :: Glevel

    !> Type of REKS calculations
    integer, intent(in) :: reksAlg

    !> Save 'A' and 'Hxc' to memory in gradient calculation
    logical, intent(in) :: tSaveMem

    !> super A hessian matrix with one-electron term in front of orbital derivatives
    real(dp), allocatable, intent(inout) :: A1e(:,:)

    !> preconditioner of super A hessian matrix with one-electron term in front of orbital
    !> derivatives
    real(dp), allocatable, intent(inout) :: A1ePre(:,:)

    !> super A hessian matrix in front of orbital derivatives
    real(dp), allocatable, intent(inout) :: Aall(:,:)

    if (Glevel == 1 .or. Glevel == 2) then

      if (tSaveMem) then

        ! build super A matrix except H-XC kernel
        call buildA1e_(Fc, Fa, omega, SAweight, FONs, G1, Nc, Na, &
            & Glevel, reksAlg, A1e, A1ePre)

      end if

    else if (Glevel == 3) then

      ! build super A matrix including H-XC kernel
      call buildAall_(eigenvecs, HxcSqrS, HxcSqrD, Fc, Fa, omega, &
          & fillingL, weight, SAweight, FONs, G1, Lpaired, Nc, &
          & Na, reksAlg, Aall)

    end if

  end subroutine getSuperAMatrix


  !> Calculate X^T vectors for state X = PPS, OSS, etc
  subroutine buildSaReksVectors(env, denseDesc, neighbourList, nNeighbourSK, &
      & iSparseStart, img2CentCell, eigenvecs, hamSqrL, hamSpL, fillingL, &
      & weightL, Nc, Na, rstate, reksAlg, tSSR, isRangeSep, XT)

    !> Environment settings
    type(TEnvironment), intent(inout) :: env

    !> Dense matrix descriptor
    type(TDenseDescr), intent(in) :: denseDesc

    !> neighbours to atoms
    type(TNeighbourList), intent(in) :: neighbourList

    !> Number of atomic neighbours
    integer, intent(in) :: nNeighbourSK(:)

    !> Index for atomic blocks in sparse data
    integer, intent(in) :: iSparseStart(:,:)

    !> map from image atom to real atoms
    integer, intent(in) :: img2CentCell(:)

    !> Eigenvectors on eixt
    real(dp), intent(in) :: eigenvecs(:,:,:)

    !> Dense Hamiltonian matrix for each microstate
    real(dp), allocatable, intent(in) :: hamSqrL(:,:,:,:)

    !> Sparse Hamiltonian matrix for each microstate
    real(dp), allocatable, intent(in) :: hamSpL(:,:,:)

    !> Filling for each microstate
    real(dp), intent(in) :: fillingL(:,:,:)

    !> Weight for each microstate per state
    real(dp), intent(in) :: weightL(:,:)

    !> Number of core orbitals
    integer, intent(in) :: Nc

    !> Number of active orbitals
    integer, intent(in) :: Na

    !> Target SSR state
    integer, intent(in) :: rstate

    !> Type of REKS calculations
    integer, intent(in) :: reksAlg

    !> Calculate SSR state with inclusion of SI, otherwise calculate SA-REKS state
    logical, intent(in) :: tSSR

    !> Whether to run a range separated calculation
    logical, intent(in) :: isRangeSep

    !> SA-REKS state vector
    real(dp), intent(out) :: XT(:,:)

    real(dp), allocatable :: tmpHam(:,:)

    real(dp) :: tmp1, tmp2
    integer :: Lmax, nOrb, Nv, superN, nstates
    integer :: iL, ij, i, j, ist

    nOrb = size(eigenvecs,dim=1)
    Lmax = size(fillingL,dim=3)
    nstates = size(weightL,dim=1)
    superN = size(XT,dim=1)
    Nv = nOrb - Nc - Na

    if (.not. isRangeSep) then
      allocate(tmpHam(nOrb,nOrb))
    end if

    XT(:,:) = 0.0_dp
    do iL = 1, Lmax

      if (isRangeSep) then

        if (tSSR) then
          do ist = 1, nstates
            do ij = 1, superN
              ! assign index i and j from ij
              call assignIndex(Nc, Na, Nv, reksAlg, ij, i, j)
              tmp1 = weightL(ist,iL)*fillingL(i,1,iL)*hamSqrL(i,j,1,iL)
              tmp2 = weightL(ist,iL)*fillingL(j,1,iL)*hamSqrL(i,j,1,iL)
              XT(ij,ist) = XT(ij,ist) + 2.0_dp * (tmp1 - tmp2)
            end do
          end do
        else
          do ij = 1, superN
            ! assign index i and j from ij
            call assignIndex(Nc, Na, Nv, reksAlg, ij, i, j)
            tmp1 = weightL(rstate,iL)*fillingL(i,1,iL)*hamSqrL(i,j,1,iL)
            tmp2 = weightL(rstate,iL)*fillingL(j,1,iL)*hamSqrL(i,j,1,iL)
            XT(ij,1) = XT(ij,1) + 2.0_dp * (tmp1 - tmp2)
          end do
        end if

      else

        ! convert from sparse to dense for hamSpL in MO basis
        tmpHam(:,:) = 0.0_dp
        ! hamSpL has (my_ud) component
        call env%globalTimer%startTimer(globalTimers%sparseToDense)
        call unpackHS(tmpHam, hamSpL(:,1,iL), &
            & neighbourList%iNeighbour, nNeighbourSK, &
            & denseDesc%iAtomStart, iSparseStart, img2CentCell)
        call env%globalTimer%stopTimer(globalTimers%sparseToDense)
        call blockSymmetrizeHS(tmpHam, denseDesc%iAtomStart)
        ! convert the multipliers from MO basis to AO basis
        call matAO2MO(tmpHam, eigenvecs(:,:,1))

        if (tSSR) then
          do ist = 1, nstates
            do ij = 1, superN
              ! assign index i and j from ij
              call assignIndex(Nc, Na, Nv, reksAlg, ij, i, j)
              tmp1 = weightL(ist,iL)*fillingL(i,1,iL)*tmpHam(i,j)
              tmp2 = weightL(ist,iL)*fillingL(j,1,iL)*tmpHam(i,j)
              XT(ij,ist) = XT(ij,ist) + 2.0_dp * (tmp1 - tmp2)
            end do
          end do
        else
          do ij = 1, superN
            ! assign index i and j from ij
            call assignIndex(Nc, Na, Nv, reksAlg, ij, i, j)
            tmp1 = weightL(rstate,iL)*fillingL(i,1,iL)*tmpHam(i,j)
            tmp2 = weightL(rstate,iL)*fillingL(j,1,iL)*tmpHam(i,j)
            XT(ij,1) = XT(ij,1) + 2.0_dp * (tmp1 - tmp2)
          end do
        end if

      end if

    end do

  end subroutine buildSaReksVectors


  !> Calculate X^T_del vectors for state-interaction
  subroutine buildInteractionVectors(eigenvecs, ZdelL, Fc, Fa, FONs, &
      & fillingL, weight, SAweight, omega, Rab, G1, Nc, Na, &
      & ia, ib, reksAlg, XTdel)

    !> Eigenvectors on eixt
    real(dp), intent(in) :: eigenvecs(:,:,:)

    !> auxiliary matrix in AO basis related to state-interaction term
    real(dp), intent(in) :: ZdelL(:,:,:)

    !> dense fock matrix for core orbitals
    real(dp), intent(in) :: Fc(:,:)

    !> dense fock matrix for active orbitals
    real(dp), intent(in) :: Fa(:,:,:)

    !> Fractional occupation numbers of active orbitals
    real(dp), intent(in) :: FONs(:,:)

    !> Filling for each microstate
    real(dp), intent(in) :: fillingL(:,:,:)

    !> Weight of each microstate for state to be optimized; weight = weightL * SAweight
    real(dp), intent(in) :: weight(:)

    !> Weights used in state-averaging
    real(dp), intent(in) :: SAweight(:)

    !> anti-symmetric matrices originated from Hamiltonians
    real(dp), intent(in) :: omega(:)

    !> state-interaction term used in SSR gradients
    real(dp), intent(in) :: Rab(:,:)

    !> constant calculated from hessian and energy of microstates
    real(dp), intent(in) :: G1

    !> Number of core orbitals
    integer, intent(in) :: Nc

    !> Number of active orbitals
    integer, intent(in) :: Na

    !> Two indices related to state-interaction term
    integer, intent(in) :: ia, ib

    !> Type of REKS calculations
    integer, intent(in) :: reksAlg

    !> state-interaction term vector
    real(dp), intent(inout) :: XTdel(:)

    real(dp), allocatable :: tmpMat(:,:)
    real(dp), allocatable :: Zmo(:,:)

    real(dp) :: fac
    integer :: Lmax, nOrb, Nv, superN
    integer :: iL, ij, i, j

    nOrb = size(fillingL,dim=1)
    Lmax = size(fillingL,dim=3)
    Nv = nOrb - Nc - Na
    superN = size(XTdel,dim=1)

    allocate(tmpMat(nOrb,nOrb))
    allocate(Zmo(nOrb,nOrb))

    ! 2-electron part
    XTdel(:) = 0.0_dp
    do iL = 1, Lmax
      ! convert ZdelL from AO basis to MO basis
      tmpMat(:,:) = 0.0_dp
      call gemm(tmpMat, ZdelL(:,:,iL), eigenvecs(:,:,1))
      Zmo(:,:) = 0.0_dp
      call gemm(Zmo, eigenvecs(:,:,1), tmpMat, transA='T')
      do ij = 1, superN
        ! assign index p and q from pq
        call assignIndex(Nc, Na, Nv, reksAlg, ij, i, j)
        fac = 2.0_dp * weight(iL) * Zmo(i,j) &
           & * (fillingL(i,1,iL) - fillingL(j,1,iL))
        XTdel(ij) = XTdel(ij) - fac
      end do
    end do

    ! 1-electron part
    select case (reksAlg)
    case (reksTypes%noReks)
    case (reksTypes%ssr22)
      call getInteraction1e22_(Fc, Fa, FONs, SAweight, omega, Rab, G1, &
          & Nc, Na, ia, ib, reksAlg, XTdel)
    case (reksTypes%ssr44)
      call error("SSR(4,4) is not implemented yet")
    end select

  end subroutine buildInteractionVectors


  !> Calculate X^T_L vector for L-th microstate
  subroutine buildLstateVector(env, denseDesc, neighbourList, nNeighbourSK, &
      & iSparseStart, img2CentCell, eigenvecs, hamSqrL, hamSpL, fillingL, &
      & Nc, Na, Lstate, Lpaired, reksAlg, isRangeSep, XTL)

    !> Environment settings
    type(TEnvironment), intent(inout) :: env

    !> Dense matrix descriptor
    type(TDenseDescr), intent(in) :: denseDesc

    !> neighbours to atoms
    type(TNeighbourList), intent(in) :: neighbourList

    !> Number of atomic neighbours
    integer, intent(in) :: nNeighbourSK(:)

    !> Index for atomic blocks in sparse data
    integer, intent(in) :: iSparseStart(:,:)

    !> map from image atom to real atoms
    integer, intent(in) :: img2CentCell(:)

    !> Eigenvectors on eixt
    real(dp), intent(in) :: eigenvecs(:,:,:)

    !> Dense Hamiltonian matrix for each microstate
    real(dp), allocatable, intent(in) :: hamSqrL(:,:,:,:)

    !> Sparse Hamiltonian matrix for each microstate
    real(dp), allocatable, intent(in) :: hamSpL(:,:,:)

    !> Filling for each microstate
    real(dp), intent(in) :: fillingL(:,:,:)

    !> Number of core orbitals
    integer, intent(in) :: Nc

    !> Number of active orbitals
    integer, intent(in) :: Na

    !> Target microstate
    integer, intent(in) :: Lstate

    !> Number of spin-paired microstates
    integer, intent(in) :: Lpaired

    !> Type of REKS calculations
    integer, intent(in) :: reksAlg

    !> Whether to run a range separated calculation
    logical, intent(in) :: isRangeSep

    !> L-th microstate vector
    real(dp), intent(out) :: XTL(:)

    real(dp), allocatable :: tmpHam(:,:)

    real(dp) :: tmp1, tmp2
    integer :: nSpin, nOrb, Nv, superN
    integer :: iS, iL, tmpL, ij, i, j

    nOrb = size(eigenvecs,dim=1)
    nSpin = size(fillingL,dim=2)
    superN = size(XTL,dim=1)
    Nv = nOrb - Nc - Na

    if (.not. isRangeSep) then
      allocate(tmpHam(nOrb,nOrb))
    end if

    if (Lstate <= Lpaired) then
      tmpL = Lstate
    else
      if (mod(Lstate,2) == 1) then
        tmpL = Lstate + 1
      else
        tmpL = Lstate - 1
      end if
    end if

    XTL(:) = 0.0_dp
    do iS = 1, nSpin

      if (iS == 1) then
        iL = Lstate
      else
        iL = tmpL
      end if

      if (isRangeSep) then

        do ij = 1, superN
          ! assign index i and j from ij
          call assignIndex(Nc, Na, Nv, reksAlg, ij, i, j)
          tmp1 = fillingL(i,1,iL)*hamSqrL(i,j,1,iL)
          tmp2 = fillingL(j,1,iL)*hamSqrL(i,j,1,iL)
          XTL(ij) = XTL(ij) + (tmp1 - tmp2)
        end do

      else

        ! convert from sparse to dense for hamSpL in MO basis
        tmpHam(:,:) = 0.0_dp
        ! hamSpL has (my_ud) component
        call env%globalTimer%startTimer(globalTimers%sparseToDense)
        call unpackHS(tmpHam, hamSpL(:,1,iL), &
            & neighbourList%iNeighbour, nNeighbourSK, &
            & denseDesc%iAtomStart, iSparseStart, img2CentCell)
        call env%globalTimer%stopTimer(globalTimers%sparseToDense)
        call blockSymmetrizeHS(tmpHam, denseDesc%iAtomStart)
        ! convert the multipliers from MO basis to AO basis
        call matAO2MO(tmpHam, eigenvecs(:,:,1))

        do ij = 1, superN
          ! assign index i and j from ij
          call assignIndex(Nc, Na, Nv, reksAlg, ij, i, j)
          tmp1 = fillingL(i,1,iL)*tmpHam(i,j)
          tmp2 = fillingL(j,1,iL)*tmpHam(i,j)
          XTL(ij) = XTL(ij) + (tmp1 - tmp2)
        end do

      end if

    end do

  end subroutine buildLstateVector


  !> solve A * Z = X equation using direct matrix inversion
  subroutine solveZT(Aall, XT, ZT)

    !> super A hessian matrix in front of orbital derivatives
    real(dp), intent(in) :: Aall(:,:)

    !> SA-REKS state vector
    real(dp), intent(in) :: XT(:)

    !> solution of A * Z = X equation with X is XT
    real(dp), intent(out) :: ZT(:)

    real(dp), allocatable :: Alu(:,:)
    integer, allocatable :: ipiv(:)
    integer :: superN

    superN = size(XT,dim=1)

    allocate(Alu(superN,superN))
    allocate(ipiv(superN))

    Alu(:,:) = Aall
    ! Alu = LU decompose
    call getrf(Alu, ipiv)
    ! Alu = A inverse
    call getri(Alu, ipiv)
    ! solve ZT = A^(-1) * XT
    call gemv(ZT, Alu, XT)

  end subroutine solveZT


  !> Calculate RmatL used in CP-REKS equations
  subroutine getRmat(eigenvecs, ZT, fillingL, Nc, Na, reksAlg, RmatL)

    !> Eigenvectors on eixt
    real(dp), intent(in) :: eigenvecs(:,:,:)

    !> solution of A * Z = X equation with X is XT
    real(dp), intent(in) :: ZT(:)

    !> Filling for each microstate
    real(dp), intent(in) :: fillingL(:,:,:)

    !> Number of core orbitals
    integer, intent(in) :: Nc

    !> Number of active orbitals
    integer, intent(in) :: Na

    !> Type of REKS calculations
    integer, intent(in) :: reksAlg

    !> auxiliary matrix in AO basis related to SA-REKS term
    real(dp), intent(out) :: RmatL(:,:,:)

    select case (reksAlg)
    case (reksTypes%noReks)
    case (reksTypes%ssr22)
      call getRmat22_(eigenvecs, ZT, fillingL, Nc, Na, &
          & reksAlg, RmatL)
    case (reksTypes%ssr44)
      call error("SSR(4,4) is not implemented yet")
    end select

  end subroutine getRmat


  !> Calculate RdelL used in CP-REKS equations
  subroutine getRdel(eigenvecs, fillingL, FONs, Nc, &
      & nstates, reksAlg, RdelL)

    !> Eigenvectors on eixt
    real(dp), intent(in) :: eigenvecs(:,:,:)

    !> Filling for each microstate
    real(dp), intent(in) :: fillingL(:,:,:)

    !> Fractional occupation numbers of active orbitals
    real(dp), intent(in) :: FONs(:,:)

    !> Number of core orbitals
    integer, intent(in) :: Nc

    !> Number of states
    integer, intent(in) :: nstates

    !> Type of REKS calculations
    integer, intent(in) :: reksAlg

    !> auxiliary matrix in AO basis related to state-interaction term
    real(dp), intent(out) :: RdelL(:,:,:,:)

    select case (reksAlg)
    case (reksTypes%noReks)
    case (reksTypes%ssr22)
      call getRdel22_(eigenvecs, fillingL, FONs, Nc, nstates, RdelL)
    case (reksTypes%ssr44)
      call error("SSR(4,4) is not implemented yet")
    end select

  end subroutine getRdel


  !> Calculate ZmatL or ZdelL used in CP-REKS equations
  subroutine getZmat(env, denseDesc, neighbourList, nNeighbourSK, &
      & iSparseStart, img2CentCell, orb, RmatL, HxcSqrS, HxcSqrD, HxcHalfS, &
      & HxcHalfD, HxcSpS, HxcSpD, overSqr, over, GammaAO, SpinAO, LrGammaAO, &
      & orderRmatL, getDenseAO, Lpaired, Glevel, tSaveMem, isRangeSep, ZmatL)

    !> Environment settings
    type(TEnvironment), intent(inout) :: env

    !> Dense matrix descriptor
    type(TDenseDescr), intent(in) :: denseDesc

    !> neighbours to atoms
    type(TNeighbourList), intent(in) :: neighbourList

    !> Number of atomic neighbours
    integer, intent(in) :: nNeighbourSK(:)

    !> Index for atomic blocks in sparse data
    integer, intent(in) :: iSparseStart(:,:)

    !> map from image atom to real atoms
    integer, intent(in) :: img2CentCell(:)

    !> atomic orbital information
    type(TOrbitals), intent(in) :: orb

    !> auxiliary matrix in AO basis related to SA-REKS term
    real(dp), intent(in) :: RmatL(:,:,:)

    !> Hartree-XC kernel with dense form with same spin part
    real(dp), allocatable, intent(in) :: HxcSqrS(:,:,:,:)

    !> Hartree-XC kernel with dense form with different spin part
    real(dp), allocatable, intent(in) :: HxcSqrD(:,:,:,:)

    !> Hartree-XC kernel with half dense form with same spin part
    real(dp), allocatable, intent(in) :: HxcHalfS(:,:)

    !> Hartree-XC kernel with half dense form with different spin part
    real(dp), allocatable, intent(in) :: HxcHalfD(:,:)

    !> Hartree-XC kernel with sparse form with same spin part
    real(dp), allocatable, intent(in) :: HxcSpS(:,:)

    !> Hartree-XC kernel with sparse form with different spin part
    real(dp), allocatable, intent(in) :: HxcSpD(:,:)

    !> Dense overlap matrix
    real(dp), intent(in) :: overSqr(:,:)

    !> sparse overlap matrix
    real(dp), intent(in) :: over(:)

    !> scc gamma integrals in AO basis
    real(dp), intent(in) :: GammaAO(:,:)

    !> spin W in AO basis
    real(dp), intent(in) :: SpinAO(:,:)

    !> long-range gamma integrals in AO basis
    real(dp), allocatable, intent(in) :: LrGammaAO(:,:)

    !> Ordering between RmatL and fillingL
    integer, intent(in) :: orderRmatL(:)

    !> get dense AO index from sparse AO array
    integer, intent(in) :: getDenseAO(:,:)

    !> Number of spin-paired microstates
    integer, intent(in) :: Lpaired

    !> Algorithms to calculate analytic gradients
    integer, intent(in) :: Glevel

    !> Save 'A' and 'Hxc' to memory in gradient calculation
    logical, intent(in) :: tSaveMem

    !> Whether to run a range separated calculation
    logical, intent(in) :: isRangeSep

    !> auxiliary matrix in AO basis related to SA-REKS term
    real(dp), intent(out) :: ZmatL(:,:,:)

    ! calculate ZmatL or ZdelL using same routine
    if (Glevel == 1 .or. Glevel == 2) then

      if (tSaveMem) then

        if (isRangeSep) then

          call getZmatHalf_(HxcHalfS, HxcHalfD, orderRmatL, Lpaired, RmatL, ZmatL)

        else

          call getZmatSparse_(env, denseDesc, neighbourList, nNeighbourSK, &
              & iSparseStart, img2CentCell, orb, over, HxcSpS, HxcSpD, &
              & orderRmatL, getDenseAO, Lpaired, RmatL, ZmatL)

        end if

      else

        call getZmatNoHxc_(env, denseDesc, neighbourList, nNeighbourSK, &
            & iSparseStart, img2CentCell, orb, getDenseAO, GammaAO, SpinAO, &
            & LrGammaAO, overSqr, RmatL, orderRmatL, Lpaired, isRangeSep, ZmatL)

      end if

    else if (Glevel == 3) then

      call getZmatDense_(HxcSqrS, HxcSqrD, orderRmatL, Lpaired, RmatL, ZmatL)

    end if

  end subroutine getZmat


  !> Calculate Q1mat used in CP-REKS equations
  subroutine getQ1mat(ZT, Fc, Fa, SAweight, FONs, &
      & Nc, Na, reksAlg, Q1mat)

    !> solution of A * Z = X equation with X is XT
    real(dp), intent(in) :: ZT(:)

    !> dense fock matrix for core orbitals
    real(dp), intent(in) :: Fc(:,:)

    !> dense fock matrix for active orbitals
    real(dp), intent(in) :: Fa(:,:,:)

    !> Weights used in state-averaging
    real(dp), intent(in) :: SAweight(:)

    !> Fractional occupation numbers of active orbitals
    real(dp), intent(in) :: FONs(:,:)

    !> Number of core orbitals
    integer, intent(in) :: Nc

    !> Number of active orbitals
    integer, intent(in) :: Na

    !> Type of REKS calculations
    integer, intent(in) :: reksAlg

    !> auxiliary matrix in MO basis related to SA-REKS term
    real(dp), intent(out) :: Q1mat(:,:)

    real(dp) :: e1, e2
    integer :: nOrb, p, q, t, pq, Nv, superN

    superN = size(ZT,dim=1)
    nOrb = size(Fc,dim=1)
    Nv = nOrb - Nc - Na

    Q1mat(:,:) = 0.0_dp
!$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(p,q,e1,e2) SCHEDULE(RUNTIME)
    do t = 1, nOrb
      do pq = 1, superN

        call assignIndex(Nc, Na, Nv, reksAlg, pq, p, q)

        call assignEpsilon(Fc, Fa, SAweight, FONs, Nc, p, q, t, &
            & 1, reksAlg, e1, e2)
        Q1mat(p,t) = Q1mat(p,t) + 0.5_dp*ZT(pq)*(e1-e2)
        call assignEpsilon(Fc, Fa, SAweight, FONs, Nc, q, p, t, &
            & 1, reksAlg, e1, e2)
        Q1mat(q,t) = Q1mat(q,t) - 0.5_dp*ZT(pq)*(e1-e2)

      end do
    end do
!$OMP END PARALLEL DO

  end subroutine getQ1mat


  !> Calculate Q1del used in CP-REKS equations
  subroutine getQ1del(eigenvecs, Fc, Fa, FONs, SAweight, &
      & Nc, nstates, reksAlg, Q1del)

    !> Eigenvectors on eixt
    real(dp), intent(in) :: eigenvecs(:,:,:)

    !> dense fock matrix for core orbitals
    real(dp), intent(in) :: Fc(:,:)

    !> dense fock matrix for active orbitals
    real(dp), intent(in) :: Fa(:,:,:)

    !> Fractional occupation numbers of active orbitals
    real(dp), intent(in) :: FONs(:,:)

    !> Weights used in state-averaging
    real(dp), intent(in) :: SAweight(:)

    !> Number of active orbitals
    integer, intent(in) :: Nc

    !> Number of states
    integer, intent(in) :: nstates

    !> Type of REKS calculations
    integer, intent(in) :: reksAlg

    !> auxiliary matrix in AO basis related to state-interaction term
    real(dp), intent(out) :: Q1del(:,:,:)

    select case (reksAlg)
    case (reksTypes%noReks)
    case (reksTypes%ssr22)
      call getQ1del22_(eigenvecs, Fc, Fa, FONs, SAweight, &
          & Nc, nstates, reksAlg, Q1del)
    case (reksTypes%ssr44)
      call error("SSR(4,4) is not implemented yet")
    end select

  end subroutine getQ1del


  !> Calculate Q2mat or Q2del used in CP-REKS equations
  subroutine getQ2mat(eigenvecs, fillingL, weight, ZmatL, Q2mat)

    !> Eigenvectors on eixt
    real(dp), intent(in) :: eigenvecs(:,:,:)

    !> Filling for each microstate
    real(dp), intent(in) :: fillingL(:,:,:)

    !> Weight of each microstate for state to be optimized; weight = weightL * SAweight
    real(dp), intent(in) :: weight(:)

    !> auxiliary matrix in AO basis related to SA-REKS term
    real(dp), intent(in) :: ZmatL(:,:,:)

    !> auxiliary matrix in AO basis related to SA-REKS term
    real(dp), intent(out) :: Q2mat(:,:)

    real(dp), allocatable :: tmpMat(:,:)
    real(dp), allocatable :: Zmo(:,:)
    real(dp) :: fac
    integer :: nOrb, Lmax, iL, p, q

    nOrb = size(fillingL,dim=1)
    Lmax = size(fillingL,dim=3)

    allocate(tmpMat(nOrb,nOrb))
    allocate(Zmo(nOrb,nOrb))

    Q2mat(:,:) = 0.0_dp
    do iL = 1, Lmax
      tmpMat(:,:) = 0.0_dp
      Zmo(:,:) = 0.0_dp
      call gemm(tmpMat, ZmatL(:,:,iL), eigenvecs(:,:,1))
      call gemm(Zmo, eigenvecs(:,:,1), tmpMat, transA='T')
      do p = 1, nOrb
        do q = p, nOrb
          fac = 2.0_dp * weight(iL) * Zmo(p,q) * (fillingL(p,1,iL) + fillingL(q,1,iL))
          Q2mat(q,p) = Q2mat(q,p) + 0.5_dp * fac
          if (p /= q) then
            Q2mat(p,q) = Q2mat(q,p)
          end if
        end do
      end do
    end do

  end subroutine getQ2mat


  !> Convert the state vectors from SA-REKS to target SSR state
  subroutine SaToSsrXT(XTdel, eigvecsSSR, rstate, XT)

    !> state-interaction term vector
    real(dp), intent(in) :: XTdel(:,:)

    !> eigenvectors from SA-REKS state
    real(dp), intent(in) :: eigvecsSSR(:,:)

    !> Target SSR state
    integer, intent(in) :: rstate

    !> SA-REKS state vector
    real(dp), intent(inout) :: XT(:,:)

    real(dp), allocatable :: tmpXT(:)
    integer :: nstates, superN, ist, jst, kst

    nstates = size(XT,dim=2)
    superN = size(XT,dim=1)

    allocate(tmpXT(superN))

    kst = 0
    tmpXT(:) = 0.0_dp
    do ist = 1, nstates
      do jst = ist, nstates
        if (ist == jst) then
          tmpXT(:) = tmpXT + eigvecsSSR(ist,rstate)**2 * XT(:,ist)
        else
          kst = kst + 1
          tmpXT(:) = tmpXT - 2.0_dp * eigvecsSSR(ist,rstate) * &
              & eigvecsSSR(jst,rstate) * XTdel(:,kst)
        end if
      end do
    end do
    XT(:,rstate) = tmpXT

  end subroutine SaToSsrXT


  !> Convert the weighting factors from SA-REKS to target SSR state
  subroutine SaToSsrWeight(Rab, weightIL, G1, eigvecsSSR, rstate, weightL)

    !> state-interaction term used in SSR gradients
    real(dp), intent(in) :: Rab(:,:)

    !> modified weight of each microstate
    real(dp), intent(in) :: weightIL(:)

    !> constant calculated from hessian and energy of microstates
    real(dp), intent(in) :: G1

    !> eigenvectors from SA-REKS state
    real(dp), intent(in) :: eigvecsSSR(:,:)

    !> Target SSR state
    integer, intent(in) :: rstate

    !> Weight for each microstate per state
    real(dp), intent(inout) :: weightL(:,:)

    real(dp), allocatable :: tmpCL(:)
    integer :: nstates, Lmax, ist, jst, iL

    nstates = size(weightL,dim=1)
    Lmax = size(weightL,dim=2)

    allocate(tmpCL(Lmax))

    tmpCL(:) = 0.0_dp
    do iL = 1, Lmax
      do ist = 1, nstates
        do jst = ist, nstates
          if (ist == jst) then
            tmpCL(iL) = tmpCL(iL) + eigvecsSSR(ist,rstate)**2 * weightL(ist,iL)
          else
            tmpCL(iL) = tmpCL(iL) - 2.0_dp * eigvecsSSR(ist,rstate) * &
                & eigvecsSSR(jst,rstate) * G1*weightIL(iL)*Rab(ist,jst)
          end if
        end do
      end do
    end do
    weightL(rstate,:) = tmpCL

  end subroutine SaToSsrWeight


  !> Convert the gradients from SA-REKS to target SSR state
  subroutine SaToSsrGradient(SAgrad, SIgrad, eigvecsSSR, SSRgrad)

    !> gradient of SA-REKS state
    real(dp), intent(in) :: SAgrad(:,:,:)

    !> gradient of state-interaction term
    real(dp), intent(in) :: SIgrad(:,:,:)

    !> eigenvectors from SA-REKS state
    real(dp), intent(in) :: eigvecsSSR(:,:)

    !> gradient of SSR state
    real(dp), intent(out) :: SSRgrad(:,:,:)

    integer :: nstates, ist, jst, rstate, kst

    nstates = size(SAgrad,dim=3)

    SSRgrad(:,:,:) = 0.0_dp
    do rstate = 1, nstates
      kst = 0
      do ist = 1, nstates
        do jst = ist, nstates
          if (ist == jst) then
            SSRgrad(:,:,rstate) = SSRgrad(:,:,rstate) &
                & + eigvecsSSR(ist,rstate)**2 * SAgrad(:,:,ist)
          else
            kst = kst + 1
            SSRgrad(:,:,rstate) = SSRgrad(:,:,rstate) &
                & + 2.0_dp * eigvecsSSR(ist,rstate) * &
                & eigvecsSSR(jst,rstate) * SIgrad(:,:,kst)
          end if
        end do
      end do
    end do

  end subroutine SaToSsrGradient


  !> Add SI contribution to R and Q matricex
  subroutine addSItoRQ(eigenvecs, RdelL, Q1del, Q2del, &
      & eigvecsSSR, rstate, RmatL, Q1mat, Q2mat, Qmat)

    !> Eigenvectors on eixt
    real(dp), intent(in) :: eigenvecs(:,:,:)

    !> auxiliary matrix in AO basis related to state-interaction term
    real(dp), intent(in) :: RdelL(:,:,:,:)

    !> auxiliary matrix in AO basis related to state-interaction term
    real(dp), intent(in) :: Q1del(:,:,:)

    !> auxiliary matrix in MO basis related to state-interaction term
    real(dp), intent(in) :: Q2del(:,:,:)

    !> eigenvectors from SA-REKS state
    real(dp), intent(in) :: eigvecsSSR(:,:)

    !> Target SSR state
    integer, intent(in) :: rstate

    !> auxiliary matrix in AO basis related to SA-REKS term
    real(dp), intent(inout) :: RmatL(:,:,:)

    !> auxiliary matrix in MO basis related to SA-REKS term
    real(dp), intent(inout) :: Q1mat(:,:)

    !> auxiliary matrix in MO basis related to SA-REKS term
    real(dp), intent(inout) :: Q2mat(:,:)

    !> auxiliary matrix in AO basis related to SA-REKS term
    real(dp), intent(out) :: Qmat(:,:)

    integer :: nstates, ist, jst, kst, iL, LmaxR

    nstates = size(eigvecsSSR,dim=1)
    LmaxR = size(RmatL,dim=3)

    do iL = 1, LmaxR
      kst = 0
      do ist = 1, nstates
        do jst = ist + 1, nstates
          kst = kst + 1
          RmatL(:,:,iL) = RmatL(:,:,iL) - 2.0_dp * eigvecsSSR(ist,rstate) * &
              & eigvecsSSR(jst,rstate) * RdelL(:,:,iL,kst)
        end do
      end do
    end do

    ! convert Q1mat from MO basis to AO basis
    call matMO2AO(Q1mat, eigenvecs(:,:,1))

    kst = 0
    do ist = 1, nstates
      do jst = ist + 1, nstates
        kst = kst + 1
        Q1mat(:,:) = Q1mat - 2.0_dp * eigvecsSSR(ist,rstate) * &
            & eigvecsSSR(jst,rstate) * Q1del(:,:,kst)
      end do
    end do

    kst = 0
    do ist = 1, nstates
      do jst = ist + 1, nstates
        kst = kst + 1
        Q2mat(:,:) = Q2mat - 2.0_dp * eigvecsSSR(ist,rstate) * &
            & eigvecsSSR(jst,rstate) * Q2del(:,:,kst)
      end do
    end do

    ! convert Q2mat + Q2del from MO basis to AO basis
    call matMO2AO(Q2mat, eigenvecs(:,:,1))

    Qmat(:,:) = Q1mat + Q2mat

  end subroutine addSItoRQ


  !> Calculate SSR gradient with 1st and 3rd terms
  subroutine SSRshift(eigenvecs, gradL, Qmat, Sderiv, ZT, SAweight, &
      & weightL, omega, weightIL, G1, iSquare, mOrb, grad, option)

    !> Eigenvectors on eixt
    real(dp), intent(in) :: eigenvecs(:,:,:)

    !> gradients for each microstate except orbital derivative terms
    real(dp), intent(in) :: gradL(:,:,:)

    !> auxiliary matrix in AO basis related to SA-REKS term
    real(dp), intent(inout) :: Qmat(:,:)

    !> Dense overlap derivative in AO basis
    real(dp), intent(in) :: Sderiv(:,:,:)

    !> solution of A * Z = X equation with X is XT
    real(dp), intent(in) :: ZT(:)

    !> Weights used in state-averaging
    real(dp), intent(in) :: SAweight(:)

    !> Weight for each microstate per state
    real(dp), intent(in) :: weightL(:)

    !> anti-symmetric matrices originated from Hamiltonians
    real(dp), intent(in) :: omega(:)

    !> modified weight of each microstate
    real(dp), intent(in) :: weightIL(:)

    !> constant calculated from hessian and energy of microstates
    real(dp), intent(in) :: G1

    !> Position of each atom in the rows/columns of the square matrices. Shape: (nAtom)
    integer, intent(in) :: iSquare(:)

    !> Max. nr. of orbitals for any species
    integer, intent(in) :: mOrb

    !> gradient choice, 1: SA-REKS, 2: SSR
    integer, intent(in) :: option

    !> gradient from 1st(H-F term) and 3rd(Q*S) terms
    real(dp), intent(out) :: grad(:,:)

    real(dp) :: tmpValue
    integer :: iL, Lmax

    Lmax = size(weightIL,dim=1)

    if (option == 1) then
      call matMO2AO(Qmat, eigenvecs(:,:,1))
    end if

    tmpValue = sum(ZT(:)*omega(:))

    grad(:,:) = 0.0_dp
    do iL = 1, Lmax
      grad(:,:) = grad + gradL(:,:,iL) * &
          & (weightL(iL) + SAweight(1)*G1*weightIL(iL)*tmpValue)
    end do

    call shiftQSgrad_(Qmat + transpose(Qmat), Sderiv, &
        & 1.0_dp, iSquare, mOrb, grad)

  end subroutine SSRshift


  !> Calculate SI gradient with 1st and 3rd terms
  subroutine SIshift(eigenvecs, gradL, Q1del, Q2del, tmpQ1, tmpQ2, &
      & Qmat, Sderiv, ZTdel, SAweight, omega, weightIL, Rab, G1, &
      & iSquare, mOrb, ist, nstates, grad)

    !> Eigenvectors on eixt
    real(dp), intent(in) :: eigenvecs(:,:,:)

    !> gradients for each microstate except orbital derivative terms
    real(dp), intent(in) :: gradL(:,:,:)

    !> auxiliary matrix in AO basis related to state-interaction term
    real(dp), intent(inout) :: Q1del(:,:)

    !> auxiliary matrix in MO basis related to state-interaction term
    real(dp), intent(inout) :: Q2del(:,:)

    !> auxiliary matrix in MO basis related to state-interaction term
    real(dp), intent(inout) :: tmpQ1(:,:)

    !> auxiliary matrix in MO basis related to state-interaction term
    real(dp), intent(inout) :: tmpQ2(:,:)

    !> auxiliary matrix in AO basis related to state-interaction term
    real(dp), intent(inout) :: Qmat(:,:)

    !> Dense overlap derivative in AO basis
    real(dp), intent(in) :: Sderiv(:,:,:)

    !> solution of A * Z = X equation with X is XTdel
    real(dp), intent(in) :: ZTdel(:)

    !> Weights used in state-averaging
    real(dp), intent(in) :: SAweight(:)

    !> anti-symmetric matrices originated from Hamiltonians
    real(dp), intent(in) :: omega(:)

    !> modified weight of each microstate
    real(dp), intent(in) :: weightIL(:)

    !> state-interaction term used in SSR gradients
    real(dp), intent(in) :: Rab(:,:)

    !> constant calculated from hessian and energy of microstates
    real(dp), intent(in) :: G1

    !> Position of each atom in the rows/columns of the square matrices. Shape: (nAtom)
    integer, intent(in) :: iSquare(:)

    !> Max. nr. of orbitals for any species
    integer, intent(in) :: mOrb

    !> Current index for loop half states
    integer, intent(in) :: ist

    !> Number of states
    integer, intent(in) :: nstates

    !> gradient from 1st(H-F term) and 3rd(Q*S) terms
    real(dp), intent(out) :: grad(:,:)

    real(dp) :: tmpValue
    integer :: iL, Lmax, ia, ib

    Lmax = size(weightIL,dim=1)

    call getTwoIndices(nstates, ist, ia, ib, 1)

    ! tmp_Q1, tmp_Q2, Q2_del : MO index, Q1_del : AO index
    call matMO2AO(tmpQ1, eigenvecs(:,:,1))
    call matMO2AO(tmpQ2, eigenvecs(:,:,1))
    call matMO2AO(Q2del, eigenvecs(:,:,1))
    Qmat(:,:) = tmpQ1 + tmpQ2 + Q1del + Q2del

    tmpValue = sum(ZTdel(:)*omega(:))

    grad(:,:) = 0.0_dp
    do iL = 1, Lmax
      grad(:,:) = grad - gradL(:,:,iL) * &
          & (G1*Rab(ia,ib)*weightIL(iL) + SAweight(1)*G1*weightIL(iL)*tmpValue)
    end do

    call shiftQSgrad_(Qmat + transpose(Qmat), Sderiv, &
        & -1.0_dp, iSquare, mOrb, grad)

  end subroutine SIshift


  !> Calculate L-th microstate gradient with 1st and 3rd terms
  subroutine Lshift(eigenvecs, gradL, Qmat, Sderiv, ZT, SAweight, &
      & omega, weightIL, G1, iSquare, mOrb, Lstate, grad)

    !> Eigenvectors on eixt
    real(dp), intent(in) :: eigenvecs(:,:,:)

    !> gradients for each microstate except orbital derivative terms
    real(dp), intent(in) :: gradL(:,:,:)

    !> auxiliary matrix in AO basis related to L state
    real(dp), intent(inout) :: Qmat(:,:)

    !> Dense overlap derivative in AO basis
    real(dp), intent(in) :: Sderiv(:,:,:)

    !> solution of A * Z = X equation with X is XTL
    real(dp), intent(in) :: ZT(:)

    !> Weights used in state-averaging
    real(dp), intent(in) :: SAweight(:)

    !> anti-symmetric matrices originated from Hamiltonians
    real(dp), intent(in) :: omega(:)

    !> modified weight of each microstate
    real(dp), intent(in) :: weightIL(:)

    !> constant calculated from hessian and energy of microstates
    real(dp), intent(in) :: G1

    !> Position of each atom in the rows/columns of the square matrices. Shape: (nAtom)
    integer, intent(in) :: iSquare(:)

    !> Max. nr. of orbitals for any species
    integer, intent(in) :: mOrb

    !> Target microstate
    integer, intent(in) :: Lstate

    !> gradient from 1st(H-F term) and 3rd(Q*S) terms
    real(dp), intent(out) :: grad(:,:)

    real(dp) :: tmpValue
    integer :: iL, Lmax

    Lmax = size(weightIL,dim=1)

    call matMO2AO(Qmat, eigenvecs(:,:,1))

    tmpValue = sum(ZT(:)*omega(:))

    grad(:,:) = 0.0_dp
    do iL = 1, Lmax
      grad(:,:) = grad + gradL(:,:,iL) * &
          & SAweight(1)*G1*weightIL(iL)*tmpValue
    end do
    grad(:,:) = grad + gradL(:,:,Lstate)

    call shiftQSgrad_(Qmat + transpose(Qmat), Sderiv, &
        & 1.0_dp, iSquare, mOrb, grad)

  end subroutine Lshift


  !> Calculate R*T contribution of gradient (2nd term)
  subroutine RTshift(env, sccCalc, denseDesc, neighbourList, nNeighbourSK, &
      & iSparseStart, img2CentCell, orb, coord0, Hderiv, Sderiv, rhoSqrL, overSqr, &
      & deltaRhoSqrL, qOutputL, q0, GammaAO, GammaDeriv, SpinAO, LrGammaAO, &
      & LrGammaDeriv, RmatL, RdelL, tmpRL, weight, extCharges, blurWidths, &
      & rVec, gVec, alpha, vol, getDenseAO, getDenseAtom, getAtomIndex, &
      & orderRmatL, Lpaired, SAstates, tNAC, isRangeSep, tExtChrg, tPeriodic, &
      & tBlur, SAgrad, SIgrad, SSRgrad)

    !> Environment settings
    type(TEnvironment), intent(inout) :: env

    !> SCC module internal variables
    type(TScc), allocatable, intent(inout) :: sccCalc

    !> Dense matrix descriptor
    type(TDenseDescr), intent(in) :: denseDesc

    !> neighbours to atoms
    type(TNeighbourList), intent(in) :: neighbourList

    !> Number of atomic neighbours
    integer, intent(in) :: nNeighbourSK(:)

    !> Index for atomic blocks in sparse data
    integer, intent(in) :: iSparseStart(:,:)

    !> map from image atom to real atoms
    integer, intent(in) :: img2CentCell(:)

    !> atomic orbital information
    type(TOrbitals), intent(in) :: orb

    !> central cell coordinates of atoms
    real(dp), intent(inout) :: coord0(:,:)


    !> Dense non-scc Hamiltonian derivative in AO basis
    real(dp), intent(in) :: Hderiv(:,:,:)

    !> Dense overlap derivative in AO basis
    real(dp), intent(in) :: Sderiv(:,:,:)

    !> Dense density matrix for each microstate
    real(dp), intent(inout) :: rhoSqrL(:,:,:,:)

    !> Dense overlap matrix
    real(dp), intent(in) :: overSqr(:,:)

    !> Dense delta density matrix for each microstate
    real(dp), allocatable, intent(in) :: deltaRhoSqrL(:,:,:,:)

    !> Mulliken population for each microstate
    real(dp), intent(in) :: qOutputL(:,:,:,:)

    !> reference atomic occupations
    real(dp), intent(in) :: q0(:,:,:)


    !> scc gamma integrals in AO basis
    real(dp), intent(in) :: GammaAO(:,:)

    !> scc gamma derivative integrals
    real(dp), intent(in) :: GammaDeriv(:,:,:)

    !> spin W in AO basis
    real(dp), intent(in) :: SpinAO(:,:)

    !> long-range gamma integrals in AO basis
    real(dp), allocatable, intent(in) :: LrGammaAO(:,:)

    !> long-range gamma derivative integrals
    real(dp), allocatable, intent(in) :: LrGammaDeriv(:,:,:)


    !> auxiliary matrix in AO basis related to SA-REKS term
    real(dp), intent(in) :: RmatL(:,:,:,:)

    !> auxiliary matrix in AO basis related to state-interaction term
    real(dp), allocatable, intent(in) :: RdelL(:,:,:,:)

    !> auxiliary matrix in AO basis related to state-interaction term
    real(dp), allocatable, intent(in) :: tmpRL(:,:,:,:)

    !> Weight of each microstate for state to be optimized; weight = weightL * SAweight
    real(dp), intent(in) :: weight(:)


    !> coordinates and charges of external point charges
    real(dp), allocatable, intent(in) :: extCharges(:,:)

    !> Width of the Gaussians if the charges are blurred
    real(dp), allocatable, intent(in) :: blurWidths(:)

    !> real lattice points for Ewald-sum
    real(dp), allocatable, intent(in) :: rVec(:,:)

    !> lattice points for reciprocal Ewald
    real(dp), allocatable, intent(in) :: gVec(:,:)

    !> parameter for Ewald
    real(dp), intent(in) :: alpha

    !> parameter for cell volume
    real(dp), intent(in) :: vol


    !> get dense AO index from sparse AO array
    integer, intent(in) :: getDenseAO(:,:)

    !> get dense atom index from sparse atom array
    integer, intent(in) :: getDenseAtom(:,:)

    !> get atom index from AO index
    integer, intent(in) :: getAtomIndex(:)

    !> Ordering between RmatL and fillingL
    integer, intent(in) :: orderRmatL(:)

    !> Number of spin-paired microstates
    integer, intent(in) :: Lpaired

    !> Number of states used in state-averaging
    integer, intent(in) :: SAstates


    !> Calculate nonadiabatic coupling vectors
    logical, intent(in) :: tNAC

    !> Whether to run a range separated calculation
    logical, intent(in) :: isRangeSep

    !> If external charges must be considered
    logical, intent(in) :: tExtChrg

    !> if calculation is periodic
    logical, intent(in) :: tPeriodic

    !> If charges should be blured
    logical, intent(in) :: tBlur


    !> gradient of SA-REKS state
    real(dp), allocatable, intent(inout) :: SAgrad(:,:,:)

    !> gradient of state-interaction term
    real(dp), allocatable, intent(inout) :: SIgrad(:,:,:)

    !> gradient of SSR state
    real(dp), intent(inout) :: SSRgrad(:,:,:)


    ! for common term
    real(dp), allocatable :: deriv1(:,:,:)
    real(dp), allocatable :: deriv2(:,:,:)

    ! for scc & pc (sparseSize), LC (nOrbHalf) term
    real(dp), allocatable :: tmpRmatL(:,:,:)
    real(dp), allocatable :: tmpRdelL(:,:,:)

    ! for LC term
    real(dp), allocatable :: SP(:,:,:)
    real(dp), allocatable :: SPS(:,:,:)

    integer :: nAtom, ist, nstates, nstHalf
    integer :: nOrb, nOrbhalf, sparseSize, LmaxR, Lmax

    nAtom = size(GammaDeriv,dim=1)
    nstates = size(RmatL,dim=4)
    nstHalf = nstates * (nstates - 1) / 2

    nOrb = size(GammaAO,dim=1)
    nOrbHalf = nOrb * (nOrb + 1) / 2
    sparseSize = size(getDenseAO,dim=1)

    LmaxR = size(RmatL,dim=3)
    Lmax = size(qOutputL,dim=4)

    allocate(deriv1(3,nAtom,nstates))
    allocate(deriv2(3,nAtom,nstHalf))
    if (tNAC) then
      allocate(tmpRmatL(sparseSize,LmaxR,nstates))
      allocate(tmpRdelL(sparseSize,LmaxR,nstHalf))
    else
      allocate(tmpRmatL(sparseSize,LmaxR,1))
    end if

    ! pack R matrix obtained from CP-REKS equation
    if (tNAC) then
      tmpRmatL(:,:,:) = 0.0_dp
      ! SA-REKS state
      do ist = 1, nstates
        if (ist /= SAstates) then
          call getRmatSp(env, denseDesc, neighbourList, nNeighbourSK, &
              & iSparseStart, img2CentCell, orb, RmatL(:,:,:,ist), &
              & tmpRmatL(:,:,ist))
        end if
      end do
      ! state-interaction term
      tmpRdelL(:,:,:) = 0.0_dp
      do ist = 1, nstHalf
        call getRmatSp(env, denseDesc, neighbourList, nNeighbourSK, &
            & iSparseStart, img2CentCell, orb, RdelL(:,:,:,ist) + &
            & tmpRL(:,:,:,ist), tmpRdelL(:,:,ist))
      end do
    else
      ! SA-REKS or SSR or L state
      tmpRmatL(:,:,:) = 0.0_dp
      call getRmatSp(env, denseDesc, neighbourList, nNeighbourSK, &
          & iSparseStart, img2CentCell, orb, RmatL(:,:,:,1), &
          & tmpRmatL(:,:,1))
    end if

    deriv1(:,:,:) = 0.0_dp
    deriv2(:,:,:) = 0.0_dp

    ! scc, spin, and pc term with sparse R and T variables
    call getSccPcTerms_(sccCalc, Hderiv, Sderiv, rhoSqrL, overSqr, &
        & qOutputL, q0, GammaAO, GammaDeriv, SpinAO, tmpRmatL, tmpRdelL, &
        & weight, getDenseAO, getAtomIndex, denseDesc%iAtomStart, &
        & orderRmatL, Lpaired, SAstates, tNAC, tExtChrg, deriv1, deriv2)

    ! point charge term with sparse R and T variables
    if (tExtChrg) then
      call getPc2ndTerms_(env, coord0, overSqr, tmpRmatL, tmpRdelL, weight, &
          & extCharges, blurWidths, rVec, gVec, alpha, vol, getDenseAO, getAtomIndex, &
          & orderRmatL, SAstates, tNAC, tPeriodic, tBlur, deriv1, deriv2)
    end if

    if (isRangeSep) then

      deallocate(tmpRmatL)
      if (tNAC) then
        deallocate(tmpRdelL)
      end if

      allocate(SP(nOrb,nOrb,Lmax))
      allocate(SPS(nOrb,nOrb,Lmax))

      call getSPmatrices(deltaRhoSqrL, overSqr, SP, SPS)

      ! LC term (gamma derivative) with dense R and T variables
      call getLr2ndTerms_(deltaRhoSqrL, overSqr, LrGammaDeriv, SP, SPS, &
          & RmatL, RdelL, tmpRL, weight, denseDesc%iAtomStart, orderRmatL, &
          & SAstates, orb%mOrb, tNAC, deriv1, deriv2)

      deallocate(SPS)

      if (tNAC) then
        allocate(tmpRmatL(nOrbHalf,LmaxR,nstates))
        allocate(tmpRdelL(nOrbHalf,LmaxR,nstHalf))
      else
        allocate(tmpRmatL(nOrbHalf,LmaxR,1))
      end if

      ! pack R matrix obtained from CP-REKS equation
      if (tNAC) then
        tmpRmatL(:,:,:) = 0.0_dp
        ! SA-REKS state
        do ist = 1, nstates
          if (ist /= SAstates) then
            call getRmatHalf(RmatL(:,:,:,ist), tmpRmatL(:,:,ist))
          end if
        end do
        ! state-interaction term
        tmpRdelL(:,:,:) = 0.0_dp
        do ist = 1, nstHalf
          call getRmatHalf(RdelL(:,:,:,ist) + tmpRL(:,:,:,ist), &
              & tmpRdelL(:,:,ist))
        end do
      else
        ! SA-REKS or SSR or L state
        tmpRmatL(:,:,:) = 0.0_dp
        call getRmatHalf(RmatL(:,:,:,1), tmpRmatL(:,:,1))
      end if

      ! LC term (overlap derivative) with half dense R and T variables
      call getLr1stTerms_(Sderiv, deltaRhoSqrL, overSqr, LrGammaAO, SP, &
          & tmpRmatL, tmpRdelL, weight, getDenseAtom, denseDesc%iAtomStart, &
          & orderRmatL, SAstates, orb%mOrb, tNAC, deriv1, deriv2)

    end if

    ! calculate the final gradient
    if (tNAC) then
      SAgrad(:,:,:) = SAgrad - deriv1
      SIgrad(:,:,:) = SIgrad + deriv2
    else
      SSRgrad(:,:,1) = SSRgrad(:,:,1) - deriv1(:,:,1)
    end if

    contains

      !> Convert RmatL from dense to sparse form
      subroutine getRmatSp(env, denseDesc, neighbourList, nNeighbourSK, &
          & iSparseStart, img2CentCell, orb, RmatL, RmatSpL)

        !> Environment settings
        type(TEnvironment), intent(inout) :: env

        !> Dense matrix descriptor
        type(TDenseDescr), intent(in) :: denseDesc

        !> neighbours to atoms
        type(TNeighbourList), intent(in) :: neighbourList

        !> Number of atomic neighbours
        integer, intent(in) :: nNeighbourSK(:)

        !> Index for atomic blocks in sparse data
        integer, intent(in) :: iSparseStart(:,:)

        !> map from image atom to real atoms
        integer, intent(in) :: img2CentCell(:)

        !> atomic orbital information
        type(TOrbitals), intent(in) :: orb

        !> auxiliary matrix in AO basis related to SA-REKS term
        real(dp), intent(in) :: RmatL(:,:,:)

        !> auxiliary matrix in AO basis related to SA-REKS term with sparse form
        real(dp), intent(inout) :: RmatSpL(:,:)

        real(dp), allocatable :: tmpMat(:,:)
        integer :: nOrb, LmaxR, mu, iL

        nOrb = size(RmatL,dim=1)
        LmaxR = size(RmatL,dim=3)

        allocate(tmpMat(nOrb,nOrb))

        do iL = 1, LmaxR

          ! In this case onsite elements remain except diagonal
          ! elements, but it does not affect the gradient
          ! since T derivative for scc, spin, pc shows only 
          ! diagonal elements in onsite block. Thus, tr(R*T)
          ! can show correct values in sparse case
          tmpMat(:,:) = RmatL(:,:,iL) + transpose(RmatL(:,:,iL))
          do mu = 1, nOrb
            tmpMat(mu,mu) = RmatL(mu,mu,iL)
          end do

          ! convert from dense to sparse for RmatL in AO basis
          call env%globalTimer%startTimer(globalTimers%denseToSparse)
          call packHS(RmatSpL(:,iL), tmpMat, &
              & neighbourlist%iNeighbour, nNeighbourSK, orb%mOrb, &
              & denseDesc%iAtomStart, iSparseStart, img2CentCell)
          call env%globalTimer%stopTimer(globalTimers%denseToSparse)

        end do

      end subroutine getRmatSp

      !> Convert RmatL from dense to half dense form
      subroutine getRmathalf(RmatL, RmatHalfL)

        !> auxiliary matrix in AO basis related to SA-REKS term
        real(dp), intent(in) :: RmatL(:,:,:)

        !> auxiliary matrix in AO basis related to SA-REKS term with half dense form
        real(dp), intent(out) :: RmatHalfL(:,:)

        real(dp), allocatable :: tmpMat(:,:)
        integer :: iL, LmaxR, mu, nu, nOrb, k, nOrbHalf

        nOrb = size(RmatL,dim=1)
        nOrbHalf = nOrb * (nOrb + 1) / 2
        LmaxR = size(RmatL,dim=3)

        allocate(tmpMat(nOrb,nOrb))

        do iL = 1, LmaxR

          tmpMat(:,:) = RmatL(:,:,iL) + transpose(RmatL(:,:,iL))
          do mu = 1, nOrb
            tmpMat(mu,mu) = RmatL(mu,mu,iL)
          end do

          do k = 1, nOrbHalf
            call getTwoIndices(nOrb, k, mu, nu, 2)
            RmatHalfL(k,iL) = tmpMat(mu,nu)
          end do

        end do

      end subroutine getRmathalf

      !> Calculate matrix product of overlap and density matrix
      subroutine getSPmatrices(deltaRhoSqrL, overSqr, SP, SPS)

        !> Dense delta density matrix for each microstate
        real(dp), intent(in) :: deltaRhoSqrL(:,:,:,:)

        !> Dense overlap matrix
        real(dp), intent(in) :: overSqr(:,:)

        !> Dense overlap * delta density in AO basis
        real(dp), intent(out) :: SP(:,:,:)

        !> Dense overlap * delta density * overlap in AO basis
        real(dp), intent(out) :: SPS(:,:,:)

        real(dp), allocatable :: tmpMat(:,:)
        integer :: nOrb, Lmax, iL

        nOrb = size(deltaRhoSqrL,dim=1)
        Lmax = size(deltaRhoSqrL,dim=4)

        allocate(tmpMat(nOrb,nOrb))

        SP(:,:,:) = 0.0_dp
        SPS(:,:,:) = 0.0_dp
        do iL = 1, Lmax
          tmpMat(:,:) = 0.0_dp
          call gemm(tmpMat, overSqr, deltaRhoSqrL(:,:,1,iL))
          ! S * Delta P calculation
          SP(:,:,iL) = tmpMat
          ! S * Delta P * S calculation
          call gemm(SPS(:,:,iL), tmpMat, overSqr)
        end do

      end subroutine getSPmatrices

  end subroutine RTshift


  !> compute the gradient for remaining SA-REKS state
  subroutine getOtherSAgrad(avgGrad, reksAlg, SAgrad)

    !> gradient of averaged state
    real(dp), intent(in) :: avgGrad(:,:)

    !> Type of REKS calculations
    integer, intent(in) :: reksAlg

    !> gradient of SA-REKS state
    real(dp), intent(inout) :: SAgrad(:,:,:)

    select case (reksAlg)
    case (reksTypes%noReks)
    case (reksTypes%ssr22)
      call getOtherSAgrad22_(avgGrad, SAgrad)
    case (reksTypes%ssr44)
      call error("SSR(4,4) is not implemented yet")
    end select

  end subroutine getOtherSAgrad


  !> In this routine, I do not use the equation in REKS document
  !> (calculated from g and h vectors), I directly calculate 
  !> G vector from SSR gradients and H vector from SA and SI gradients
  !> compute the NAC vectors (G, H)
  subroutine getReksNAC(SAgrad, SIgrad, SSRgrad, eigvecsSSR, energy, nacG, nacH)

    !> gradient of SA-REKS state
    real(dp), intent(in) :: SAgrad(:,:,:)

    !> gradient of state-interaction term
    real(dp), intent(in) :: SIgrad(:,:,:)

    !> gradient of SSR state
    real(dp), intent(in) :: SSRgrad(:,:,:)

    !> eigenvectors from SA-REKS state
    real(dp), intent(in) :: eigvecsSSR(:,:)

    !> energy of states
    real(dp), intent(in) :: energy(:)

    !> difference gradient vector, G
    real(dp), intent(out) :: nacG(:,:,:)

    !> nonadiabatic coupling vector, H
    real(dp), intent(out) :: nacH(:,:,:)

    integer :: ast, bst, cst, ist, jst, kst, nstates

    nstates = size(SAgrad,dim=3)

    nacG(:,:,:) = 0.0_dp
    nacH(:,:,:) = 0.0_dp
    cst = 0
    do ast = 1, nstates
      do bst = ast + 1, nstates

        cst = cst + 1

        ! calculate G vector from SSR gradient
        nacG(:,:,cst) = 0.5_dp * (SSRgrad(:,:,ast) - SSRgrad(:,:,bst))

        ! calculate H vector, non-adiabatic coupling
        kst = 0
        do ist = 1, nstates
          do jst = ist, nstates
            if (ist == jst) then
              nacH(:,:,cst) = nacH(:,:,cst) + SAgrad(:,:,ist) * &
                  & eigvecsSSR(ist,ast) * eigvecsSSR(ist,bst)
            else
              kst = kst + 1
              nacH(:,:,cst) = nacH(:,:,cst) + SIgrad(:,:,kst) * &
                  & ( eigvecsSSR(ist,ast) * eigvecsSSR(jst,bst) &
                  & + eigvecsSSR(jst,ast) * eigvecsSSR(ist,bst) )
            end if
          end do
        end do
        nacH(:,:,cst) = nacH(:,:,cst) / (energy(bst) - energy(ast))

      end do
    end do

  end subroutine getReksNAC


  !> Calculate external charge gradients for target state
  subroutine getExtChrgGradients(env, qmCoords, pcCoords, qOutput, q0, &
      & pcCharges, blurWidths, rVec, gVec, alpha, vol, tPeriodic, tBlur, chrgForces)

    !> Environment settings
    type(TEnvironment), intent(inout) :: env

    !> atomic coordinates
    real(dp), intent(in) :: qmCoords(:,:)

    !> coordinates of external point charges
    real(dp), intent(in) :: pcCoords(:,:)

    !> Output electrons
    real(dp), intent(in) :: qOutput(:,:,:)

    !> reference atomic charges
    real(dp), intent(in) :: q0(:,:,:)

    !> charges of external point charges
    real(dp), intent(in) :: pcCharges(:)

    !> Width of the Gaussians if the charges are blurred
    real(dp), allocatable, intent(in) :: blurWidths(:)

    !> real lattice points for Ewald-sum
    real(dp), allocatable, intent(in) :: rVec(:,:)

    !> lattice points for reciprocal Ewald
    real(dp), allocatable, intent(in) :: gVec(:,:)

    !> parameter for Ewald
    real(dp), intent(in) :: alpha

    !> parameter for cell volume
    real(dp), intent(in) :: vol

    !> if calculation is periodic
    logical, intent(in) :: tPeriodic

    !> If charges should be blured
    logical, intent(in) :: tBlur

    !> forces on external charges
    real(dp), intent(inout) :: chrgForces(:,:)

    real(dp), allocatable :: qmCharges(:)
    real(dp), allocatable :: deriv(:,:)

    integer :: iAt, nAtom, nAtomPc

    nAtom = size(qmCoords,dim=2)
    nAtomPc = size(pcCoords,dim=2)

    allocate(qmCharges(nAtom))
    allocate(deriv(3,nAtom))

    ! get the charge of QM system.
    ! actually, charge should be q0 - qOutput, but the defined charge
    ! in dftb+ is qOutput - q0 for simplicity of coding.
    qmCharges(:) = 0.0_dp
    do iAt = 1, nAtom
      qmCharges(iAt) = sum(qOutput(:,iAt,1) - q0(:,iAt,1),dim=1)
    end do

    ! chrgForces is gradient, not force
    ! currently, qmCoords, extCharges : au, not A
    deriv(:,:) = 0.0_dp
    chrgForces(:,:) = 0.0_dp
    if (tPeriodic) then
      if (tBlur) then
        call addInvRPrime(env, nAtom, nAtomPc, qmCoords, pcCoords, qmCharges, &
            & pcCharges, rVec, gVec, alpha, vol, deriv, chrgForces, tHamDeriv=.false., &
            & blurWidths1=blurWidths)
      else
        call addInvRPrime(env, nAtom, nAtomPc, qmCoords, pcCoords, qmCharges, &
            & pcCharges, rVec, gVec, alpha, vol, deriv, chrgForces, tHamDeriv=.false.)
      end if
    else
      if (tBlur) then
        call addInvRPrime(env, nAtom, nAtomPc, qmCoords, pcCoords, qmCharges, &
            & pcCharges, deriv, chrgForces, tHamDeriv=.false., blurWidths1=blurWidths)
      else
        call addInvRPrime(env, nAtom, nAtomPc, qmCoords, pcCoords, qmCharges, &
            & pcCharges, deriv, chrgForces, tHamDeriv=.false.)
      end if
    end if

  end subroutine getExtChrgGradients


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!! Private routines
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Calculate H-XC kernel for DFTB in AO basis with dense form
  subroutine HxcKernelDense_(overSqr, GammaAO, SpinAO, LrGammaAO, isRangeSep, HxcSqrS, HxcSqrD)

    !> Dense overlap matrix
    real(dp), intent(in) :: overSqr(:,:)

    !> scc gamma integrals in AO basis
    real(dp), intent(in) :: GammaAO(:,:)

    !> spin W in AO basis
    real(dp), intent(in) :: SpinAO(:,:)

    !> long-range gamma integrals in AO basis
    real(dp), allocatable, intent(in) :: LrGammaAO(:,:)

    !> Whether to run a range separated calculation
    logical, intent(in) :: isRangeSep

    !> Hartree-XC kernel with dense form with same spin part
    real(dp), allocatable, intent(inout) :: HxcSqrS(:,:,:,:)

    !> Hartree-XC kernel with dense form with different spin part
    real(dp), allocatable, intent(inout) :: HxcSqrD(:,:,:,:)

    ! common variables
    integer :: nOrb, mu, nu, tau, gam

    ! scc/spin variables
    real(dp) :: tmpG1, tmpG2, tmpG3, tmpG4
    real(dp) :: tmpS1, tmpS2, tmpS3, tmpS4

    ! LC variables
    real(dp) :: tmpL1, tmpL2, tmpL3, tmpL4

    nOrb = size(overSqr,dim=1)

    ! zeroing for dense H-XC kernel
    HxcSqrS(:,:,:,:) = 0.0_dp
    HxcSqrD(:,:,:,:) = 0.0_dp

!$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(tmpG1,tmpG2,tmpG3,tmpG4, &
!$OMP& tmpS1,tmpS2,tmpS3,tmpS4,tmpL1,tmpL2,tmpL3,tmpL4) SCHEDULE(RUNTIME)
    do mu = 1, nOrb
    do nu = 1, nOrb

      do tau = 1, nOrb
      do gam = 1, nOrb

        ! s = same spin (s1=s2), d = different spin (s1/=s2)
        tmpG1 = GammaAO(mu,tau)
        tmpG2 = GammaAO(nu,tau)
        tmpG3 = GammaAO(mu,gam)
        tmpG4 = GammaAO(nu,gam)

        tmpS1 = SpinAO(mu,tau)
        tmpS2 = SpinAO(nu,tau)
        tmpS3 = SpinAO(mu,gam)
        tmpS4 = SpinAO(nu,gam)

        ! compute Hxc kernel with respect to AO basis
        HxcSqrS(mu,nu,tau,gam) = 0.25_dp * overSqr(nu,mu) * overSqr(gam,tau) * &
            & ( (tmpG1+tmpG2+tmpG3+tmpG4) + (tmpS1+tmpS2+tmpS3+tmpS4) )
        HxcSqrD(mu,nu,tau,gam) = 0.25_dp * overSqr(nu,mu) * overSqr(gam,tau) * &
            & ( (tmpG1+tmpG2+tmpG3+tmpG4) - (tmpS1+tmpS2+tmpS3+tmpS4) )

        if (isRangeSep) then

          tmpL1 = LRgammaAO(mu,gam)
          tmpL2 = LRgammaAO(mu,nu)
          tmpL3 = LRgammaAO(tau,gam)
          tmpL4 = LRgammaAO(tau,nu)

          HxcSqrS(mu,nu,tau,gam) = HxcSqrS(mu,nu,tau,gam) - 0.25_dp &
              & * overSqr(mu,tau) * overSqr(nu,gam) * (tmpL1+tmpL2+tmpL3+tmpL4)

        end if

      end do
      end do

    end do
    end do
!$OMP END PARALLEL DO

  end subroutine HxcKernelDense_


  !> Calculate H-XC kernel for DFTB in AO basis with half dense form
  subroutine HxcKernelHalf_(getDenseAO, overSqr, GammaAO, SpinAO, LrGammaAO, isRangeSep, HxcHalfS,&
      & HxcHalfD)

    !> get dense AO index from sparse AO array
    integer, intent(in) :: getDenseAO(:,:)

    !> Dense overlap matrix
    real(dp), intent(in) :: overSqr(:,:)

    !> scc gamma integrals in AO basis
    real(dp), intent(in) :: GammaAO(:,:)

    !> spin W in AO basis
    real(dp), intent(in) :: SpinAO(:,:)

    !> long-range gamma integrals in AO basis
    real(dp), intent(in) :: LrGammaAO(:,:)

    !> Whether to run a range separated calculation
    logical, intent(in) :: isRangeSep

    !> Hartree-XC kernel with half dense form with same spin part
    real(dp), allocatable, intent(inout) :: HxcHalfS(:,:)

    !> Hartree-XC kernel with half dense form with different spin part
    real(dp), allocatable, intent(inout) :: HxcHalfD(:,:)

    ! common variables
    real(dp) :: tmp22
    integer :: ii, jj, sparseSize, k, l, nOrbHalf
    integer :: mu, nu, tau, gam, nOrb

    ! scc/spin variables
    real(dp) :: tmpG1, tmpG2, tmpG3, tmpG4
    real(dp) :: tmpS1, tmpS2, tmpS3, tmpS4

    ! LC variables
    real(dp) :: tmpL1, tmpL2, tmpL3, tmpL4
    real(dp) :: tmpvalue1, tmpvalue2

    nOrb = size(overSqr,dim=1)
    sparseSize = size(getDenseAO,dim=1)
    nOrbHalf = nOrb * (nOrb + 1) / 2

    ! zeroing for half dense H-XC kernel
    HxcHalfS(:,:) = 0.0_dp
    HxcHalfD(:,:) = 0.0_dp

    ! scc, spin part
!$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(mu,nu,k,tau,gam,l, &
!$OMP& tmpG1,tmpG2,tmpG3,tmpG4,tmpS1,tmpS2,tmpS3,tmpS4) SCHEDULE(RUNTIME)
    do ii = 1, sparseSize

      ! set the AO indices with respect to sparsity
      mu = getDenseAO(ii,1)
      nu = getDenseAO(ii,2)

      if (mu <= nu .and. abs(overSqr(mu,nu)) >= epsilon(1.0_dp)) then

        ! calculate the index in terms of half dense form
        k = (mu-1)*nOrb - mu*(mu-1)/2 + nu

        do jj = 1, sparseSize

          ! set the AO indices with respect to sparsity
          tau = getDenseAO(jj,1)
          gam = getDenseAO(jj,2)

          if (tau <= gam .and. abs(overSqr(tau,gam)) /= epsilon(1.0_dp)) then

            ! calculate the index in terms of half dense form
            l = (tau-1)*nOrb - tau*(tau-1)/2 + gam

            ! scc part
            tmpG1 = GammaAO(mu,tau)
            tmpG2 = GammaAO(nu,tau)
            tmpG3 = GammaAO(mu,gam)
            tmpG4 = GammaAO(nu,gam)

            ! spin part
            tmpS1 = SpinAO(mu,tau)
            tmpS2 = SpinAO(nu,tau)
            tmpS3 = SpinAO(mu,gam)
            tmpS4 = SpinAO(nu,gam)

            ! s = same spin (s1=s2), d = different spin (s1/=s2)
            HxcHalfS(k,l) = 0.25_dp * overSqr(mu,nu) * overSqr(tau,gam) * &
                & ( (tmpG1+tmpG2+tmpG3+tmpG4) + (tmpS1+tmpS2+tmpS3+tmpS4) )
            HxcHalfD(k,l) = 0.25_dp * overSqr(mu,nu) * overSqr(tau,gam) * &
                & ( (tmpG1+tmpG2+tmpG3+tmpG4) - (tmpS1+tmpS2+tmpS3+tmpS4) )

          end if

        end do

      end if

    end do
!$OMP END PARALLEL DO

    ! LC terms
    if (isRangeSep) then

!$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(mu,nu,tau,gam,tmp22, &
!$OMP& tmpL1,tmpL2,tmpL3,tmpL4,tmpvalue1,tmpvalue2) SCHEDULE(RUNTIME)
      do k = 1, nOrbHalf

        call getTwoIndices(nOrb, k, mu, nu, 2)

        do l = 1, nOrbHalf

          ! TODO : use 'getTwoIndices' -> more times
          !        use tmp22 part -> relatively less times
          !        Why show different cost?
          !call getTwoIndices(nOrb, l, tau, gam, 2)
          tmp22 = ( real(2.0_dp*nOrb+3.0_dp, dp) - sqrt( (2.0_dp*nOrb+ &
              & 3.0_dp)**2.0_dp - 8.0_dp*(nOrb+l) ) )/2.0_dp
          tau = int( real(tmp22, dp) )
          gam = tau**2/2 - tau/2 - nOrb*tau + nOrb + l

          ! LC terms
          if (isRangeSep) then

            ! (mu,nu,tau,gam)
            tmpvalue1 = 0.0_dp
            if (abs(overSqr(mu,tau)) >= epsilon(1.0_dp) .and. &
                & abs(overSqr(nu,gam)) >= epsilon(1.0_dp)) then
              tmpL1 = LrGammaAO(mu,gam)
              tmpL2 = LrGammaAO(mu,nu)
              tmpL3 = LrGammaAO(tau,gam)
              tmpL4 = LrGammaAO(tau,nu)
              tmpvalue1 = -0.125_dp * overSqr(mu,tau) * &
                  & overSqr(nu,gam) * (tmpL1+tmpL2+tmpL3+tmpL4)
            end if

            ! (mu,nu,gam,tau)
            tmpvalue2 = 0.0_dp
            if (abs(overSqr(mu,gam)) >= epsilon(1.0_dp) .and. &
                & abs(overSqr(nu,tau)) >= epsilon(1.0_dp)) then
              tmpL1 = LrGammaAO(mu,tau)
              tmpL2 = LrGammaAO(mu,nu)
              tmpL3 = LrGammaAO(gam,tau)
              tmpL4 = LrGammaAO(gam,nu)
              tmpvalue2 = -0.125_dp * overSqr(mu,gam) * &
                  & overSqr(nu,tau) * (tmpL1+tmpL2+tmpL3+tmpL4)
            end if

            HxcHalfS(k,l) = HxcHalfS(k,l) + tmpvalue1 + tmpvalue2

          end if

        end do

      end do
!$OMP END PARALLEL DO

    end if

  end subroutine HxcKernelHalf_


  !> Calculate H-XC kernel for DFTB in AO basis with sparse form
  subroutine HxcKernelSparse_(getDenseAO, over, GammaAO, SpinAO, &
      & HxcSpS, HxcSpD)

    !> get dense AO index from sparse AO array
    integer, intent(in) :: getDenseAO(:,:)

    !> sparse overlap matrix
    real(dp), intent(in) :: over(:)

    !> scc gamma integrals in AO basis
    real(dp), intent(in) :: GammaAO(:,:)

    !> spin W in AO basis
    real(dp), intent(in) :: SpinAO(:,:)

    !> Hartree-XC kernel with sparse form with same spin part
    real(dp), allocatable, intent(inout) :: HxcSpS(:,:)

    !> Hartree-XC kernel with sparse form with different spin part
    real(dp), allocatable, intent(inout) :: HxcSpD(:,:)

    ! common variables
    integer :: sparseSize, ii, jj, mu, nu, tau, gam

    ! scc/spin variables
    real(dp) :: tmpG1, tmpG2, tmpG3, tmpG4
    real(dp) :: tmpS1, tmpS2, tmpS3, tmpS4

    sparseSize = size(over,dim=1)

    ! zeroing for sparse H-XC kernel
    HxcSpS(:,:) = 0.0_dp
    HxcSpD(:,:) = 0.0_dp

!$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(mu,nu,tau,gam,tmpG1,tmpG2, &
!$OMP& tmpG3,tmpG4,tmpS1,tmpS2,tmpS3,tmpS4) SCHEDULE(RUNTIME)
    do ii = 1, sparseSize
      if (abs(over(ii)) >= epsilon(1.0_dp)) then

        ! set the AO indices with respect to sparsity
        mu = getDenseAO(ii,1)
        nu = getDenseAO(ii,2)

        if (mu <= nu) then

          do jj = 1, sparseSize
            if (abs(over(jj)) >= epsilon(1.0_dp)) then

              ! set the AO indices with respect to sparsity
              tau = getDenseAO(jj,1)
              gam = getDenseAO(jj,2)

              if (tau <= gam) then

                tmpG1 = GammaAO(mu,tau)
                tmpG2 = GammaAO(nu,tau)
                tmpG3 = GammaAO(mu,gam)
                tmpG4 = GammaAO(nu,gam)

                tmpS1 = SpinAO(mu,tau)
                tmpS2 = SpinAO(nu,tau)
                tmpS3 = SpinAO(mu,gam)
                tmpS4 = SpinAO(nu,gam)

                ! s = same spin (s1=s2), d = different spin (s1/=s2)
                HxcSpS(jj,ii) = 0.25_dp * over(ii) * over(jj) * &
                    & ( (tmpG1+tmpG2+tmpG3+tmpG4) + (tmpS1+tmpS2+tmpS3+tmpS4) )
                HxcSpD(jj,ii) = 0.25_dp * over(ii) * over(jj) * &
                    & ( (tmpG1+tmpG2+tmpG3+tmpG4) - (tmpS1+tmpS2+tmpS3+tmpS4) )

              end if

            end if
          end do

        end if

      end if
    end do
!$OMP END PARALLEL DO

  end subroutine HxcKernelSparse_


  !> Calculate 1st contribution of Rab in REKS(2,2)
  subroutine getRab22_1st_(hamSqr, fillingL, weightIL, SAweight, FONs, Nc, iL, Rab)

    !> Dense Hamiltonian matrix for each microstate
    real(dp), intent(in) :: hamSqr(:,:)

    !> Filling for each microstate
    real(dp), intent(in) :: fillingL(:,:,:)

    !> modified weight of each microstate
    real(dp), intent(in) :: weightIL(:)

    !> Weights used in state-averaging
    real(dp), intent(in) :: SAweight(:)

    !> Fractional occupation numbers of active orbitals
    real(dp), intent(in) :: FONs(:,:)

    !> Number of core orbitals
    integer, intent(in) :: Nc

    !> Current index for loop L
    integer, intent(in) :: iL

    !> state-interaction term used in SSR gradients
    real(dp), intent(inout) :: Rab(:,:)

    real(dp) :: n_a, n_b
    integer :: a, b, nstates

    n_a = FONs(1,1)
    n_b = FONs(2,1)

    nstates = size(Rab,dim=1)
    a = Nc + 1
    b = Nc + 2

    Rab(1,2) = Rab(1,2) + 2.0_dp * SAweight(1) * weightIL(iL) * &
        & ( sqrt(n_a)*fillingL(a,1,iL)*hamSqr(b,a) &
        & - sqrt(n_b)*fillingL(b,1,iL)*hamSqr(a,b) )
    if (nstates == 3) then
      Rab(2,3) = Rab(2,3) + 2.0_dp * SAweight(1) * weightIL(iL) * &
          & ( sqrt(n_a)*fillingL(a,1,iL)*hamSqr(b,a) &
          & + sqrt(n_b)*fillingL(b,1,iL)*hamSqr(a,b) )
    end if

  end subroutine getRab22_1st_


  !> Calculate 2nd contribution of Rab in REKS(2,2)
  subroutine getRab22_2nd_(fockFa, FONs, SAweight, Nc, Rab)

    !> dense fock matrix for active orbitals
    real(dp), intent(in) :: fockFa(:,:,:)

    !> Fractional occupation numbers of active orbitals
    real(dp), intent(in) :: FONs(:,:)

    !> Weights used in state-averaging
    real(dp), intent(in) :: SAweight(:)

    !> Number of core orbitals
    integer, intent(in) :: Nc

    !> state-interaction term used in SSR gradients
    real(dp), intent(inout) :: Rab(:,:)

    real(dp) :: n_a, n_b, e1
    integer :: a, b, nstates

    n_a = FONs(1,1)
    n_b = FONs(2,1)

    nstates = size(Rab,dim=1)
    a = Nc + 1
    b = Nc + 2

    e1 = fockFa(b,a,1) * (SAweight(1)*n_a + SAweight(2))
    Rab(1,2) = Rab(1,2) + (1.0_dp/sqrt(n_a) + 1.0_dp/sqrt(n_b)) * e1
    Rab(2,1) = Rab(1,2)
    if (nstates == 3) then
      Rab(2,3) = Rab(2,3) + (1.0_dp/sqrt(n_a) - 1.0_dp/sqrt(n_b)) * e1
      Rab(3,2) = Rab(2,3)
    end if

  end subroutine getRab22_2nd_


  !> Calculate super A hessian matrix without H-XC kernel
  subroutine buildA1e_(Fc, Fa, omega, SAweight, FONs, G1, Nc, Na, &
      & Glevel, reksAlg, A1e, A1ePre)

    !> dense fock matrix for core orbitals
    real(dp), intent(in) :: Fc(:,:)

    !> dense fock matrix for active orbitals
    real(dp), intent(in) :: Fa(:,:,:)

    !> anti-symmetric matrices originated from Hamiltonians
    real(dp), intent(in) :: omega(:)

    !> Weights used in state-averaging
    real(dp), intent(in) :: SAweight(:)

    !> Fractional occupation numbers of active orbitals
    real(dp), intent(in) :: FONs(:,:)

    !> constant calculated from hessian and energy of microstates
    real(dp), intent(in) :: G1

    !> Number of core orbitals
    integer, intent(in) :: Nc

    !> Number of active orbitals
    integer, intent(in) :: Na

    !> Algorithms to calculate analytic gradients
    integer, intent(in) :: Glevel

    !> Type of REKS calculations
    integer, intent(in) :: reksAlg

    !> super A hessian matrix with one-electron term in front of orbital derivatives
    real(dp), allocatable, intent(inout) :: A1e(:,:)

    !> preconditioner of super A hessian matrix with one-electron term in front of orbital
    !> derivatives
    real(dp), allocatable, intent(inout) :: A1ePre(:,:)

    real(dp) :: e1, e2
    integer :: nOrb, superN, Nv, superNhalf
    integer :: ij, pq, i, j, p, q, ijpq

    nOrb = size(Fc,dim=1)
    Nv = nOrb - Nc - Na
    superN = size(A1e,dim=1)
    superNhalf = superN * (superN + 1) / 2

    A1e(:,:) = 0.0_dp
!$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(ij,i,j,pq,p,q,e1,e2) SCHEDULE(RUNTIME)
    do ijpq = 1, superNhalf

      call getTwoIndices(superN, ijpq, ij, pq, 2)

      ! assign index i and j from ij
      call assignIndex(Nc, Na, Nv, reksAlg, ij, i, j)

      ! assign index p and q from pq
      call assignIndex(Nc, Na, Nv, reksAlg, pq, p, q)

      ! get lagrange multipliers with delta function
      if (i == p) then
        e1 = 0.0_dp; e2 = 0.0_dp;
        call assignEpsilon(Fc, Fa, SAweight, FONs, Nc, i, j, q, &
            & 1, reksAlg, e1, e2)
        A1e(ij,pq) = A1e(ij,pq) + 0.5_dp*(e1 - e2)
      end if

      if (j == q) then
        e1 = 0.0_dp; e2 = 0.0_dp;
        call assignEpsilon(Fc, Fa, SAweight, FONs, Nc, i, j, p, &
            & 2, reksAlg, e1, e2)
        A1e(ij,pq) = A1e(ij,pq) - 0.5_dp*(e1 - e2)
      end if

      if (i == q) then
        e1 = 0.0_dp; e2 = 0.0_dp;
        call assignEpsilon(Fc, Fa, SAweight, FONs, Nc, i, j, p, &
            & 1, reksAlg, e1, e2)
        A1e(ij,pq) = A1e(ij,pq) - 0.5_dp*(e1 - e2)
      end if

      if (j == p) then
        e1 = 0.0_dp; e2 = 0.0_dp;
        call assignEpsilon(Fc, Fa, SAweight, FONs, Nc, i, j, q, &
            & 2, reksAlg, e1, e2)
        A1e(ij,pq) = A1e(ij,pq) + 0.5_dp*(e1 - e2)
      end if

      ! SAweight(1) is equal to W0
      A1e(ij,pq) = A1e(ij,pq) - SAweight(1) * G1 * omega(ij) * omega(pq)

      ! remaining part of super A matrix
      if (ij /= pq) then
        A1e(pq,ij) = A1e(ij,pq)
      end if

    end do
!$OMP END PARALLEL DO

    A1ePre(:,:) = 0.0_dp
    do ij = 1, superN

      ! assign index i and j from ij
      call assignIndex(Nc, Na, Nv, reksAlg, ij, i, j)

      if (Glevel == 1) then

        ! get lagrange multipliers with delta function
        e1 = 0.0_dp; e2 = 0.0_dp;
        call assignEpsilon(Fc, Fa, SAweight, FONs, Nc, i, j, j, &
            & 1, reksAlg, e1, e2)
        A1ePre(ij,ij) = A1ePre(ij,ij) + 0.5_dp*(e1 - e2)

        e1 = 0.0_dp; e2 = 0.0_dp;
        call assignEpsilon(Fc, Fa, SAweight, FONs, Nc, i, j, i, &
            & 2, reksAlg, e1, e2)
        A1ePre(ij,ij) = A1ePre(ij,ij) - 0.5_dp*(e1 - e2)

        ! SAweight(1) is equal to W0
        A1ePre(ij,ij) = A1ePre(ij,ij) - SAweight(1) * G1 * omega(ij) * omega(ij)

      else if (Glevel == 2) then

        A1ePre(ij,ij) = 1.0_dp

      end if

      ! check singularity for preconditioner
      if (abs(A1ePre(ij,ij)) <= epsilon(1.0_dp)) then
        write(stdOut,'(A,f15.8)') " Current preconditioner value = ", A1ePre(ij,ij)
        call error("A singularity exists in preconditioner for PCG, set Preconditioner = No")
      end if

      ! preconditioner part for CG
      A1ePre(ij,ij) = 1.0_dp / A1ePre(ij,ij)

    end do

  end subroutine buildA1e_


  !> Calculate super A hessian matrix with H-XC kernel
  subroutine buildAall_(eigenvecs, HxcSqrS, HxcSqrD, Fc, Fa, omega, &
      & fillingL, weight, SAweight, FONs, G1, Lpaired, Nc, Na, &
      & reksAlg, Aall)

    !> Eigenvectors on eixt
    real(dp), intent(in) :: eigenvecs(:,:,:)

    !> Hartree-XC kernel with dense form with same spin part
    real(dp), intent(in) :: HxcSqrS(:,:,:,:)

    !> Hartree-XC kernel with dense form with different spin part
    real(dp), intent(in) :: HxcSqrD(:,:,:,:)

    !> dense fock matrix for core orbitals
    real(dp), intent(in) :: Fc(:,:)

    !> dense fock matrix for active orbitals
    real(dp), intent(in) :: Fa(:,:,:)

    !> anti-symmetric matrices originated from Hamiltonians
    real(dp), intent(in) :: omega(:)

    !> Filling for each microstate
    real(dp), intent(in) :: fillingL(:,:,:)

    !> Weight of each microstate for state to be optimized; weight = weightL * SAweight
    real(dp), intent(in) :: weight(:)

    !> Weights used in state-averaging
    real(dp), intent(in) :: SAweight(:)

    !> Fractional occupation numbers of active orbitals
    real(dp), intent(in) :: FONs(:,:)

    !> constant calculated from hessian and energy of microstates
    real(dp), intent(in) :: G1

    !> Number of spin-paired microstates
    integer, intent(in) :: Lpaired

    !> Number of core orbitals
    integer, intent(in) :: Nc

    !> Number of active orbitals
    integer, intent(in) :: Na

    !> Type of REKS calculations
    integer, intent(in) :: reksAlg

    !> super A hessian matrix in front of orbital derivatives
    real(dp), allocatable, intent(inout) :: Aall(:,:)

    real(dp), allocatable :: HxcTot(:,:,:,:)
    real(dp) :: e1, e2
    integer :: nOrb, superN, Nv
    integer :: ij, pq, i, j, p, q

    nOrb = size(Fc,dim=1)
    Nv = nOrb - Nc - Na
    superN = size(Aall,dim=1)

    allocate(HxcTot(nOrb,nOrb,nOrb,nOrb))

    HxcTot(:,:,:,:) = 0.0_dp
    call getHxcMo_(eigenvecs(:,:,1), HxcSqrS, HxcSqrD, fillingL, weight, &
        & superN, Lpaired, Nc, Na, reksAlg, HxcTot)

    Aall(:,:) = 0.0_dp
!$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(ij,i,j,pq,p,q,e1,e2) SCHEDULE(RUNTIME)
    do ij = 1, superN

      ! assign index i and j from ij
      call assignIndex(Nc, Na, Nv, reksAlg, ij, i, j)

      do pq = ij, superN

        ! assign index p and q from pq
        call assignIndex(Nc, Na, Nv, reksAlg, pq, p, q)

        ! get lagrange multipliers with delta function
        if (i == p) then
          e1 = 0.0_dp; e2 = 0.0_dp;
          call assignEpsilon(Fc, Fa, SAweight, FONs, Nc, i, j, q, &
              & 1, reksAlg, e1, e2)
          Aall(ij,pq) = Aall(ij,pq) + 0.5_dp*(e1 - e2)
        end if

        if (j == q) then
          e1 = 0.0_dp; e2 = 0.0_dp;
          call assignEpsilon(Fc, Fa, SAweight, FONs, Nc, i, j, p, &
              & 2, reksAlg, e1, e2)
          Aall(ij,pq) = Aall(ij,pq) - 0.5_dp*(e1 - e2)
        end if

        if (i == q) then
          e1 = 0.0_dp; e2 = 0.0_dp;
          call assignEpsilon(Fc, Fa, SAweight, FONs, Nc, i, j, p, &
              & 1, reksAlg, e1, e2)
          Aall(ij,pq) = Aall(ij,pq) - 0.5_dp*(e1 - e2)
        end if

        if (j == p) then
          e1 = 0.0_dp; e2 = 0.0_dp;
          call assignEpsilon(Fc, Fa, SAweight, FONs, Nc, i, j, q, &
              & 2, reksAlg, e1, e2)
          Aall(ij,pq) = Aall(ij,pq) + 0.5_dp*(e1 - e2)
        end if

        ! SAweight(1) is equal to W0
        Aall(ij,pq) = Aall(ij,pq) - SAweight(1) * G1 * omega(ij) * omega(pq)

        ! get HxcTot with respect to MO basis
        Aall(ij,pq) = Aall(ij,pq) + HxcTot(i,j,q,p)

        ! remaining part of super A matrix
        if (ij /= pq) then
          Aall(pq,ij) = Aall(ij,pq)
        end if

      end do
    end do
!$OMP END PARALLEL DO

  end subroutine buildAall_


  !> Calculate H-XC kernel in MO basis, it shows very high computational cost
  subroutine getHxcMo_(eigenvecs, HxcSqrS, HxcSqrD, fillingL, weight, &
      & superN, Lpaired, Nc, Na, reksAlg, HxcTot)

    !> Eigenvectors on eixt
    real(dp), intent(in) :: eigenvecs(:,:)

    !> Hartree-XC kernel with dense form with same spin part
    real(dp), intent(in) :: HxcSqrS(:,:,:,:)

    !> Hartree-XC kernel with dense form with different spin part
    real(dp), intent(in) :: HxcSqrD(:,:,:,:)

    !> Filling for each microstate
    real(dp), intent(in) :: fillingL(:,:,:)

    !> Weight of each microstate for state to be optimized; weight = weightL * SAweight
    real(dp), intent(in) :: weight(:)

    !> Size of super matrix
    integer, intent(in) :: superN

    !> Number of spin-paired microstates
    integer, intent(in) :: Lpaired

    !> Number of core orbitals
    integer, intent(in) :: Nc

    !> Number of active orbitals
    integer, intent(in) :: Na

    !> Type of REKS calculations
    integer, intent(in) :: reksAlg

    !> total averaged H-XC kernel in MO basis
    real(dp), intent(inout) :: HxcTot(:,:,:,:)

    real(dp), allocatable :: tmpMat(:,:)
    real(dp), allocatable :: HxcSmo(:,:,:,:)
    real(dp), allocatable :: HxcDmo(:,:,:,:)

    real(dp) :: tmpHxcS, tmpHxcD
    integer :: nOrb, i, j, p, q, ij, pq
    integer :: Lmax, iL, Nv

    nOrb = size(eigenvecs,dim=1)
    Lmax = size(fillingL,dim=3)
    Nv = nOrb - Nc - Na

    allocate(tmpMat(nOrb,nOrb))
    allocate(HxcSmo(nOrb,nOrb,nOrb,nOrb))
    allocate(HxcDmo(nOrb,nOrb,nOrb,nOrb))

    HxcSmo(:,:,:,:) = 0.0_dp
    HxcDmo(:,:,:,:) = 0.0_dp
!$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(i,j,tmpMat) SCHEDULE(RUNTIME)
    do i = 1, nOrb
    do j = 1, nOrb
      tmpMat(:,:) = 0.0_dp
      call gemm(tmpMat,HxcSqrS(:,:,i,j),eigenvecs,1.0_dp,0.0_dp,'N','N')
      call gemm(HxcSmo(:,:,i,j),eigenvecs,tmpMat,1.0_dp,0.0_dp,'T','N')
      tmpMat(:,:) = 0.0_dp
      call gemm(tmpMat,HxcSqrD(:,:,i,j),eigenvecs,1.0_dp,0.0_dp,'N','N')
      call gemm(HxcDmo(:,:,i,j),eigenvecs,tmpMat,1.0_dp,0.0_dp,'T','N')
    end do
    end do
!$OMP END PARALLEL DO
!$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(i,j,tmpMat) SCHEDULE(RUNTIME)
    do i = 1, nOrb
    do j = 1, nOrb
      tmpMat(:,:) = 0.0_dp
      call gemm(tmpMat,HxcSmo(i,j,:,:),eigenvecs,1.0_dp,0.0_dp,'N','N')
      call gemm(HxcSmo(i,j,:,:),eigenvecs,tmpMat,1.0_dp,0.0_dp,'T','N')
      tmpMat(:,:) = 0.0_dp
      call gemm(tmpMat,HxcDmo(i,j,:,:),eigenvecs,1.0_dp,0.0_dp,'N','N')
      call gemm(HxcDmo(i,j,:,:),eigenvecs,tmpMat,1.0_dp,0.0_dp,'T','N')
    end do
    end do
!$OMP END PARALLEL DO

!$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(ij,i,j,pq,p,q, &
!$OMP& iL,tmpHxcS,tmpHxcD) SCHEDULE(RUNTIME)
    do ij = 1, superN

      ! assign index i and j from ij
      call assignIndex(Nc, Na, Nv, reksAlg, ij, i, j)

      do pq = 1, superN

        ! assign index p and q from pq
        call assignIndex(Nc, Na, Nv, reksAlg, pq, p, q)

        tmpHxcS = 0.5_dp * (HxcSmo(i,j,p,q) + HxcSmo(i,j,q,p))
        tmpHxcD = 0.5_dp * (HxcDmo(i,j,p,q) + HxcDmo(i,j,q,p))

        do iL = 1, Lmax

          if (iL <= Lpaired) then
            HxcTot(i,j,q,p) = HxcTot(i,j,q,p) + 2.0_dp * weight(iL) * &
                 & ( fillingL(i,1,iL) - fillingL(j,1,iL) ) &
                 & * ( tmpHxcS * ( fillingL(p,1,iL) - fillingL(q,1,iL) ) &
                 & + tmpHxcD * ( fillingL(p,1,iL) - fillingL(q,1,iL) ) )
          else
            if (mod(iL,2) == 1) then
              HxcTot(i,j,q,p) = HxcTot(i,j,q,p) + 2.0_dp * weight(iL) * &
                   & ( fillingL(i,1,iL) - fillingL(j,1,iL) ) &
                   & * ( tmpHxcS * ( fillingL(p,1,iL) - fillingL(q,1,iL) ) &
                   & + tmpHxcD * ( fillingL(p,1,iL+1) - fillingL(q,1,iL+1) ) )
            else
              HxcTot(i,j,q,p) = HxcTot(i,j,q,p) + 2.0_dp * weight(iL) * &
                   & ( fillingL(i,1,iL) - fillingL(j,1,iL) ) &
                   & * ( tmpHxcS * ( fillingL(p,1,iL) - fillingL(q,1,iL) ) &
                   & + tmpHxcD * ( fillingL(p,1,iL-1) - fillingL(q,1,iL-1) ) )
            end if
          end if

        end do

      end do
    end do
!$OMP END PARALLEL DO

  end subroutine getHxcMo_


  !> Calculate X^T_del vectors for 1e contribution in REKS(2,2)
  subroutine getInteraction1e22_(Fc, Fa, FONs, SAweight, omega, Rab, G1, &
      & Nc, Na, ia, ib, reksAlg, XTdel)

    !> dense fock matrix for core orbitals
    real(dp), intent(in) :: Fc(:,:)

    !> dense fock matrix for active orbitals
    real(dp), intent(in) :: Fa(:,:,:)

    !> Fractional occupation numbers of active orbitals
    real(dp), intent(in) :: FONs(:,:)

    !> Weights used in state-averaging
    real(dp), intent(in) :: SAweight(:)

    !> anti-symmetric matrices originated from Hamiltonians
    real(dp), intent(in) :: omega(:)

    !> state-interaction term used in SSR gradients
    real(dp), intent(in) :: Rab(:,:)

    !> constant calculated from hessian and energy of microstates
    real(dp), intent(in) :: G1

    !> Number of core orbitals
    integer, intent(in) :: Nc

    !> Number of active orbitals
    integer, intent(in) :: Na

    !> Two indices related to state-interaction term
    integer, intent(in) :: ia, ib

    !> Type of REKS calculations
    integer, intent(in) :: reksAlg

    !> state-interaction term vector
    real(dp), intent(inout) :: XTdel(:)

    real(dp) :: e1, e2, n_a, n_b
    integer :: a, b, ij, i, j, Nv, nOrb, superN

    superN = size(XTdel,dim=1)
    nOrb = size(Fc,dim=1)
    Nv = nOrb - Nc - Na
    a = Nc + 1
    b = Nc + 2

    n_a = FONs(1,1)
    n_b = FONs(2,1)

    ! 1-electron part
    do ij = 1, superN

      ! assign index i and j from ij
      call assignIndex(Nc, Na, Nv, reksAlg, ij, i, j)

      if (j == b) then
        call assignEpsilon(Fc, Fa, SAweight, FONs, Nc, a, b, i, &
            & 3, reksAlg, e1, e2)
        if (ia == 1 .and. ib == 2) then
          XTdel(ij) = XTdel(ij) + 0.5_dp*( sqrt(n_a)*e1 - sqrt(n_b)*e2 )
        else if (ia == 2 .and. ib == 3) then
          XTdel(ij) = XTdel(ij) + 0.5_dp*( sqrt(n_a)*e1 + sqrt(n_b)*e2 )
        end if
      end if
      if (j == a) then
        call assignEpsilon(Fc, Fa, SAweight, FONs, Nc, b, a, i, &
            & 3, reksAlg, e1, e2)
        if (ia == 1 .and. ib == 2) then
          XTdel(ij) = XTdel(ij) + 0.5_dp*( sqrt(n_a)*e2 - sqrt(n_b)*e1 )
        else if (ia == 2 .and. ib == 3) then
          XTdel(ij) = XTdel(ij) + 0.5_dp*( sqrt(n_a)*e2 + sqrt(n_b)*e1 )
        end if
      end if

      ! the following parts are not written in REKS document (summary ver.)
      ! it is described in REKS full document (derivation ver.)
      if (i == a) then
        call assignEpsilon(Fc, Fa, SAweight, FONs, Nc, b, a, j, &
            & 3, reksAlg, e1, e2)
        if (ia == 1 .and. ib == 2) then
          XTdel(ij) = XTdel(ij) - 0.5_dp*( sqrt(n_a)*e2 )
        else if (ia == 2 .and. ib == 3) then
          XTdel(ij) = XTdel(ij) - 0.5_dp*( sqrt(n_a)*e2 )
        end if
      end if
      if (i == b) then
        call assignEpsilon(Fc, Fa, SAweight, FONs, Nc, a, b, j, &
            & 3, reksAlg, e1, e2)
        if (ia == 1 .and. ib == 2) then
          XTdel(ij) = XTdel(ij) + 0.5_dp*( sqrt(n_b)*e2 )
        else if (ia == 2 .and. ib == 3) then
          XTdel(ij) = XTdel(ij) - 0.5_dp*( sqrt(n_b)*e2 )
        end if
      end if
      if (i == a .and. j == b) then
        call assignEpsilon(Fc, Fa, SAweight, FONs, Nc, b, b, b, &
            & 3, reksAlg, e1, e2)
        if (ia == 1 .and. ib == 2) then
          XTdel(ij) = XTdel(ij) + 0.5_dp*( sqrt(n_b)*e1 )
        else if (ia == 2 .and. ib == 3) then
          XTdel(ij) = XTdel(ij) - 0.5_dp*( sqrt(n_b)*e1 )
        end if
      end if

      XTdel(ij) = XTdel(ij) + G1 * Rab(ia,ib) * omega(ij)

    end do

  end subroutine getInteraction1e22_


  !> Calculate RmatL used in CP-REKS equations in REKS(2,2)
  subroutine getRmat22_(eigenvecs, ZT, fillingL, Nc, Na, &
      & reksAlg, RmatL)

    !> Eigenvectors on eixt
    real(dp), intent(in) :: eigenvecs(:,:,:)

    !> solution of A * Z = X equation with X is XT
    real(dp), intent(in) :: ZT(:)

    !> Filling for each microstate
    real(dp), intent(in) :: fillingL(:,:,:)

    !> Number of core orbitals
    integer, intent(in) :: Nc

    !> Number of active orbitals
    integer, intent(in) :: Na

    !> Type of REKS calculations
    integer, intent(in) :: reksAlg

    !> auxiliary matrix in AO basis related to SA-REKS term
    real(dp), intent(out) :: RmatL(:,:,:)

    real(dp), allocatable :: tmpZT(:,:)
    real(dp), allocatable :: tmpMat(:,:)

    integer :: superN, LmaxR, nOrb, Nv, ij, i, j, iL

    superN = size(ZT,dim=1)
    LmaxR = size(RmatL,dim=3)
    nOrb = size(eigenvecs,dim=1)
    Nv = nOrb - Nc - Na

    allocate(tmpZT(nOrb,nOrb))
    allocate(tmpMat(nOrb,nOrb))

    RmatL(:,:,:) = 0.0_dp
    do iL = 1, LmaxR
      tmpZT(:,:) = 0.0_dp
      do ij = 1, superN
        ! assign index i and j from ij
        call assignIndex(Nc, Na, Nv, reksAlg, ij, i, j)
        if (iL <= 2) then
          tmpZT(i,j) = ZT(ij) * (fillingL(i,1,iL) - fillingL(j,1,iL))
        else
          tmpZT(i,j) = ZT(ij) * (fillingL(i,1,iL+2) - fillingL(j,1,iL+2))
        end if
      end do
      tmpMat(:,:) = 0.0_dp
      call gemm(tmpMat, tmpZT, eigenvecs(:,:,1), transB='T')
      call gemm(RmatL(:,:,iL), eigenvecs(:,:,1), tmpMat)
    end do

  end subroutine getRmat22_


  !> Calculate RdelL used in CP-REKS equations in REKS(2,2)
  subroutine getRdel22_(eigenvecs, fillingL, FONs, Nc, nstates, RdelL)

    !> Eigenvectors on eixt
    real(dp), intent(in) :: eigenvecs(:,:,:)

    !> Filling for each microstate
    real(dp), intent(in) :: fillingL(:,:,:)

    !> Fractional occupation numbers of active orbitals
    real(dp), intent(in) :: FONs(:,:)

    !> Number of core orbitals
    integer, intent(in) :: Nc

    !> Number of states
    integer, intent(in) :: nstates

    !> auxiliary matrix in AO basis related to state-interaction term
    real(dp), intent(out) :: RdelL(:,:,:,:)

    real(dp) :: n_a, n_b
    integer :: mu, nu, nOrb, a, b

    nOrb = size(eigenvecs,dim=1)
    a = Nc + 1
    b = Nc + 2

    n_a = FONs(1,1)
    n_b = FONs(2,1)

    RdelL(:,:,:,:) = 0.0_dp
    do mu = 1, nOrb
      do nu = 1, nOrb

        ! 1st microstate (up) = 1st down = 3rd up = 4th down
        RdelL(nu,mu,1,1) = eigenvecs(nu,b,1) * eigenvecs(mu,a,1) * &
          & ( sqrt(n_a)*fillingL(a,1,1) - sqrt(n_b)*fillingL(b,1,1) )
        ! 2nd microstate (up) = 2nd down = 3rd down = 4th up
        RdelL(nu,mu,2,1) = eigenvecs(nu,b,1) * eigenvecs(mu,a,1) * &
          & ( sqrt(n_a)*fillingL(a,1,2) - sqrt(n_b)*fillingL(b,1,2) )
        ! 5th microstate (up) = 6th down
        RdelL(nu,mu,3,1) = eigenvecs(nu,b,1) * eigenvecs(mu,a,1) * &
          & ( sqrt(n_a)*fillingL(a,1,5) - sqrt(n_b)*fillingL(b,1,5) )
        ! 6th microstate (up) = 5th down
        RdelL(nu,mu,4,1) = eigenvecs(nu,b,1) * eigenvecs(mu,a,1) * &
          & ( sqrt(n_a)*fillingL(a,1,6) - sqrt(n_b)*fillingL(b,1,6) )

        if (nstates == 3) then
          ! 1st microstate (up) = 1st down = 3rd up = 4th down
          RdelL(nu,mu,1,3) = eigenvecs(nu,b,1) * eigenvecs(mu,a,1) * &
            & ( sqrt(n_a)*fillingL(a,1,1) + sqrt(n_b)*fillingL(b,1,1) )
          ! 2nd microstate (up) = 2nd down = 3rd down = 4th up
          RdelL(nu,mu,2,3) = eigenvecs(nu,b,1) * eigenvecs(mu,a,1) * &
            & ( sqrt(n_a)*fillingL(a,1,2) + sqrt(n_b)*fillingL(b,1,2) )
          ! 5th microstate (up) = 6th down
          RdelL(nu,mu,3,3) = eigenvecs(nu,b,1) * eigenvecs(mu,a,1) * &
            & ( sqrt(n_a)*fillingL(a,1,5) + sqrt(n_b)*fillingL(b,1,5) )
          ! 6th microstate (up) = 5th down
          RdelL(nu,mu,4,3) = eigenvecs(nu,b,1) * eigenvecs(mu,a,1) * &
            & ( sqrt(n_a)*fillingL(a,1,6) + sqrt(n_b)*fillingL(b,1,6) )
        end if

      end do
    end do

  end subroutine getRdel22_


  !> Calculate ZmatL with dense form used in CP-REKS equations in REKS(2,2)
  subroutine getZmatDense_(HxcSqrS, HxcSqrD, orderRmatL, Lpaired, RmatL, ZmatL)

    !> Hartree-XC kernel with dense form with same spin part
    real(dp), intent(in) :: HxcSqrS(:,:,:,:)

    !> Hartree-XC kernel with dense form with different spin part
    real(dp), intent(in) :: HxcSqrD(:,:,:,:)

    !> Ordering between RmatL and fillingL
    integer, intent(in) :: orderRmatL(:)

    !> Number of spin-paired microstates
    integer, intent(in) :: Lpaired

    !> auxiliary matrix in AO basis related to SA-REKS term
    real(dp), intent(in) :: RmatL(:,:,:)

    !> auxiliary matrix in AO basis related to SA-REKS term
    real(dp), intent(out) :: ZmatL(:,:,:)

    real(dp), allocatable :: tmpHxcS(:,:)
    real(dp), allocatable :: tmpHxcD(:,:)
    real(dp), allocatable :: tmpZM(:,:)

    integer :: iL, Lmax, LmaxR, nOrb, nOrbHalf, gam, tau, ii

    nOrb = size(RmatL,dim=1)
    nOrbHalf = nOrb * (nOrb + 1) / 2
    Lmax = size(ZmatL,dim=3)
    LmaxR = size(RmatL,dim=3)

    allocate(tmpHxcS(nOrb,nOrb))
    allocate(tmpHxcD(nOrb,nOrb))
    allocate(tmpZM(2,LmaxR))

    ZmatL(:,:,:) = 0.0_dp
!$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(gam,tau, &
!$OMP& tmpHxcS,tmpHxcD,tmpZM) SCHEDULE(RUNTIME)
    do ii = 1, nOrbHalf

      call getTwoIndices(nOrb, ii, gam, tau, 2)

      tmpHxcS(:,:) = HxcSqrS(:,:,gam,tau) + HxcSqrS(:,:,tau,gam)
      tmpHxcD(:,:) = HxcSqrD(:,:,gam,tau) + HxcSqrD(:,:,tau,gam)

      ! calculate the ZmatL
      do iL = 1, LmaxR
        tmpZM(1,iL) = sum(RmatL(:,:,iL)*tmpHxcS)
        tmpZM(2,iL) = sum(RmatL(:,:,iL)*tmpHxcD)
      end do
      tmpZM(:,:) = 0.5_dp * tmpZM
      call getZmatLoop_(tmpZM, orderRmatL, Lpaired, gam, tau, ZmatL)

    end do
!$OMP END PARALLEL DO

    do iL = 1, Lmax
      call symmetrizeHS(ZmatL(:,:,iL))
    end do

  end subroutine getZmatDense_


  !> Calculate ZmatL with half dense form used in CP-REKS equations in REKS(2,2)
  subroutine getZmatHalf_(HxcHalfS, HxcHalfD, orderRmatL, Lpaired, RmatL, ZmatL)

    !> Hartree-XC kernel with half dense form with same spin part
    real(dp), intent(in) :: HxcHalfS(:,:)

    !> Hartree-XC kernel with half dense form with different spin part
    real(dp), intent(in) :: HxcHalfD(:,:)

    !> Ordering between RmatL and fillingL
    integer, intent(in) :: orderRmatL(:)

    !> Number of spin-paired microstates
    integer, intent(in) :: Lpaired

    !> auxiliary matrix in AO basis related to SA-REKS term
    real(dp), intent(in) :: RmatL(:,:,:)

    !> auxiliary matrix in AO basis related to SA-REKS term
    real(dp), intent(out) :: ZmatL(:,:,:)

    real(dp), allocatable :: tmpHxcS(:)
    real(dp), allocatable :: tmpHxcD(:)
    real(dp), allocatable :: RmatHalfL(:,:)
    real(dp), allocatable :: tmpMat(:,:)
    real(dp), allocatable :: tmpZM(:,:)

    integer :: Lmax, LmaxR, nOrb, nOrbHalf, iL, gam, tau, ii

    nOrb = size(RmatL,dim=1)
    LmaxR = size(RmatL,dim=3)
    Lmax = size(ZmatL,dim=3)
    nOrbHalf = nOrb * (nOrb + 1) / 2

    allocate(tmpHxcS(nOrbHalf))
    allocate(tmpHxcD(nOrbHalf))
    allocate(RmatHalfL(nOrbHalf,LmaxR))
    allocate(tmpMat(nOrb,nOrb))
    allocate(tmpZM(2,LmaxR))

    RmatHalfL(:,:) = 0.0_dp
    do iL = 1, LmaxR
      tmpMat(:,:) = RmatL(:,:,iL) + transpose(RmatL(:,:,iL))
      do gam = 1, nOrb
        tmpMat(gam,gam) = RmatL(gam,gam,iL)
      end do
      do ii = 1, nOrbHalf
        call getTwoIndices(nOrb, ii, gam, tau, 2)
        RmatHalfL(ii,iL) = tmpMat(gam,tau)
      end do
    end do

    ZmatL(:,:,:) = 0.0_dp
!$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(gam,tau, &
!$OMP& tmpHxcS,tmpHxcD,tmpZM) SCHEDULE(RUNTIME)
    do ii = 1, nOrbHalf

      call getTwoIndices(nOrb, ii, gam, tau, 2)

      tmpHxcS(:) = HxcHalfS(:,ii)
      tmpHxcD(:) = HxcHalfD(:,ii)

      ! calculate the ZmatL
      do iL = 1, LmaxR
        tmpZM(1,iL) = sum(RmatHalfL(:,iL)*tmpHxcS)
        tmpZM(2,iL) = sum(RmatHalfL(:,iL)*tmpHxcD)
      end do
      call getZmatLoop_(tmpZM, orderRmatL, Lpaired, gam, tau, ZmatL)

    end do
!$OMP END PARALLEL DO

    do iL = 1, Lmax
      call symmetrizeHS(ZmatL(:,:,iL))
    end do

  end subroutine getZmatHalf_


  !> Calculate ZmatL with sparse form used in CP-REKS equations in REKS(2,2)
  subroutine getZmatSparse_(env, denseDesc, neighbourList, &
      & nNeighbourSK, iSparseStart, img2CentCell, orb, over, &
      & HxcSpS, HxcSpD, orderRmatL, getDenseAO, Lpaired, RmatL, ZmatL)

    !> Environment settings
    type(TEnvironment), intent(inout) :: env

    !> Dense matrix descriptor
    type(TDenseDescr), intent(in) :: denseDesc

    !> neighbours to atoms
    type(TNeighbourList), intent(in) :: neighbourList

    !> Number of atomic neighbours
    integer, intent(in) :: nNeighbourSK(:)

    !> Index for atomic blocks in sparse data
    integer, intent(in) :: iSparseStart(:,:)

    !> map from image atom to real atoms
    integer, intent(in) :: img2CentCell(:)

    !> atomic orbital information
    type(TOrbitals), intent(in) :: orb

    !> sparse overlap matrix
    real(dp), intent(in) :: over(:)

    !> Hartree-XC kernel with sparse form with same spin part
    real(dp), intent(in) :: HxcSpS(:,:)

    !> Hartree-XC kernel with sparse form with different spin part
    real(dp), intent(in) :: HxcSpD(:,:)

    !> Ordering between RmatL and fillingL
    integer, intent(in) :: orderRmatL(:)

    !> get dense AO index from sparse AO array
    integer, intent(in) :: getDenseAO(:,:)

    !> Number of spin-paired microstates
    integer, intent(in) :: Lpaired

    !> auxiliary matrix in AO basis related to SA-REKS term
    real(dp), intent(in) :: RmatL(:,:,:)

    !> auxiliary matrix in AO basis related to SA-REKS term
    real(dp), intent(out) :: ZmatL(:,:,:)

    real(dp), allocatable :: tmpHxcS(:)
    real(dp), allocatable :: tmpHxcD(:)
    real(dp), allocatable :: RmatSpL(:,:)
    real(dp), allocatable :: tmpMat(:,:)
    real(dp), allocatable :: tmpZM(:,:)

    integer :: ii, tau, gam, nOrb, sparseSize, iL, LmaxR, Lmax

    nOrb = size(RmatL,dim=1)
    LmaxR = size(RmatL,dim=3)
    Lmax = size(ZmatL,dim=3)
    sparseSize = size(HxcSpS,dim=1)

    allocate(tmpHxcS(sparseSize))
    allocate(tmpHxcD(sparseSize))
    allocate(RmatSpL(sparseSize,LmaxR))
    allocate(tmpMat(nOrb,nOrb))
    allocate(tmpZM(2,LmaxR))

    ! pack R matrix
    RmatSpL(:,:) = 0.0_dp
    do iL = 1, LmaxR

      ! In this case onsite elements remain except diagonal
      ! elements, but it does not affect the gradient
      ! since H-XC kernel for scc, spin shows only 
      ! diagonal elements in onsite block. Thus, sum(R*Hxc)
      ! can show correct values in sparse case
      tmpMat(:,:) = RmatL(:,:,iL) + transpose(RmatL(:,:,iL))
      do gam = 1, nOrb
        tmpMat(gam,gam) = RmatL(gam,gam,iL)
      end do

      ! convert from dense to sparse for RmatL in AO basis
      call env%globalTimer%startTimer(globalTimers%denseToSparse)
      call packHS(RmatSpL(:,iL), tmpMat, neighbourlist%iNeighbour, &
          & nNeighbourSK, orb%mOrb, denseDesc%iAtomStart, iSparseStart, &
          & img2CentCell)
      call env%globalTimer%stopTimer(globalTimers%denseToSparse)

    end do

    ZmatL(:,:,:) = 0.0_dp
!$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(tau,gam,tmpHxcS,tmpHxcD,tmpZM) SCHEDULE(RUNTIME)
    do ii = 1, size(over,dim=1) ! tau and gam

      if (abs(over(ii)) >= epsilon(1.0_dp)) then

        ! set the AO indices with respect to sparsity
        tau = getDenseAO(ii,1)
        gam = getDenseAO(ii,2)

        tmpHxcS(:) = HxcSpS(:,ii)
        tmpHxcD(:) = HxcSpD(:,ii)

        ! calculate the ZmatL
        do iL = 1, LmaxR
          tmpZM(1,iL) = sum(RmatSpL(:,iL)*tmpHxcS)
          tmpZM(2,iL) = sum(RmatSpL(:,iL)*tmpHxcD)
        end do
        call getZmatLoop_(tmpZM, orderRmatL, Lpaired, tau, gam, ZmatL)

      end if

    end do
!$OMP END PARALLEL DO

    do iL = 1, Lmax
      call symmetrizeHS(ZmatL(:,:,iL))
    end do

  end subroutine getZmatSparse_


  !> Calculate ZmatL without saving H-XC kernel used in CP-REKS equations in REKS(2,2)
  subroutine getZmatNoHxc_(env, denseDesc, neighbourList, nNeighbourSK, &
      & iSparseStart, img2CentCell, orb, getDenseAO, GammaAO, SpinAO, &
      & LrGammaAO, overSqr, RmatL, orderRmatL, Lpaired, isRangeSep, ZmatL)

    !> Environment settings
    type(TEnvironment), intent(inout) :: env

    !> Dense matrix descriptor
    type(TDenseDescr), intent(in) :: denseDesc

    !> neighbours to atoms
    type(TNeighbourList), intent(in) :: neighbourList

    !> Number of atomic neighbours
    integer, intent(in) :: nNeighbourSK(:)

    !> Index for atomic blocks in sparse data
    integer, intent(in) :: iSparseStart(:,:)

    !> map from image atom to real atoms
    integer, intent(in) :: img2CentCell(:)

    !> atomic orbital information
    type(TOrbitals), intent(in) :: orb

    !> get dense AO index from sparse AO array
    integer, intent(in) :: getDenseAO(:,:)

    !> scc gamma integrals in AO basis
    real(dp), intent(in) :: GammaAO(:,:)

    !> spin W in AO basis
    real(dp), intent(in) :: SpinAO(:,:)

    !> long-range gamma integrals in AO basis
    real(dp), allocatable, intent(in) :: LrGammaAO(:,:)

    !> Dense overlap matrix
    real(dp), intent(in) :: overSqr(:,:)

    !> auxiliary matrix in AO basis related to SA-REKS term
    real(dp), intent(in) :: RmatL(:,:,:)

    !> Ordering between RmatL and fillingL
    integer, intent(in) :: orderRmatL(:)

    !> Number of spin-paired microstates
    integer, intent(in) :: Lpaired

    !> Whether to run a range separated calculation
    logical, intent(in) :: isRangeSep

    !> auxiliary matrix in AO basis related to SA-REKS term
    real(dp), intent(out) :: ZmatL(:,:,:)

    real(dp), allocatable :: tmpRmatL(:,:)
    real(dp), allocatable :: tmpHxcS(:)
    real(dp), allocatable :: tmpHxcD(:)
    real(dp), allocatable :: tmpMat(:,:)
    real(dp), allocatable :: tmpZM(:,:)

    ! common variables
    real(dp) :: tmp22!, tmp11
!    real(dp) :: tmpZ1, tmpZ2, tmp1, tmp2
    integer :: mu, nu, tau, gam, nOrb, iL, Lmax, LmaxR
    integer :: ii, jj, sparseSize, nOrbHalf

    ! scc/spin variables
    real(dp) :: tmpG1, tmpG2, tmpG3, tmpG4
    real(dp) :: tmpS1, tmpS2, tmpS3, tmpS4

    ! LC variables
    real(dp) :: tmpL1, tmpL2, tmpL3, tmpL4
    real(dp) :: tmpvalue1, tmpvalue2

    nOrb = size(RmatL,dim=1)
    LmaxR = size(RmatL,dim=3)
    Lmax = size(ZmatL,dim=3)
    sparseSize = size(getDenseAO,dim=1)
    nOrbHalf = nOrb * (nOrb + 1) / 2

    allocate(tmpRmatL(sparseSize,LmaxR))
    allocate(tmpHxcS(sparseSize))
    allocate(tmpHxcD(sparseSize))
    allocate(tmpMat(nOrb,nOrb))
    allocate(tmpZM(2,LmaxR))

    ! pack R matrix
    tmpRmatL(:,:) = 0.0_dp
    do iL = 1, LmaxR

      ! In this case onsite elements remain except diagonal
      ! elements, but it does not affect the gradient
      ! since H-XC kernel for scc, spin shows only 
      ! diagonal elements in onsite block. Thus, sum(R*Hxc)
      ! can show correct values in sparse case
      tmpMat(:,:) = RmatL(:,:,iL) + transpose(RmatL(:,:,iL))
      do mu = 1, nOrb
        tmpMat(mu,mu) = RmatL(mu,mu,iL)
      end do

      ! convert from dense to sparse for RmatL in AO basis
      call env%globalTimer%startTimer(globalTimers%denseToSparse)
      call packHS(tmpRmatL(:,iL), tmpMat, &
          & neighbourlist%iNeighbour, nNeighbourSK, orb%mOrb, &
          & denseDesc%iAtomStart, iSparseStart, img2CentCell)
      call env%globalTimer%stopTimer(globalTimers%denseToSparse)

    end do

    ! zeroing for ZmatL
    ZmatL(:,:,:) = 0.0_dp

    ! calculate ZmatL for scc and spin term

!$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(mu,nu,tau,gam, &
!$OMP& tmpG1,tmpG2,tmpG3,tmpG4,tmpS1,tmpS2,tmpS3,tmpS4, &
!$OMP& tmpHxcS,tmpHxcD,tmpZM) SCHEDULE(RUNTIME)
    do ii = 1, sparseSize

      ! set the AO indices with respect to sparsity
      tau = getDenseAO(ii,1)
      gam = getDenseAO(ii,2)

      if (tau <= gam .and. abs(overSqr(tau,gam)) >= epsilon(1.0_dp)) then

        tmpHxcS(:) = 0.0_dp
        tmpHxcD(:) = 0.0_dp
        do jj = 1, sparseSize

          ! set the AO indices with respect to sparsity
          mu = getDenseAO(jj,1)
          nu = getDenseAO(jj,2)

          ! calculate the H-XC kernel for scc, spin term
          if (mu <= nu .and. abs(overSqr(mu,nu)) >= epsilon(1.0_dp)) then

            tmpG1 = GammaAO(mu,tau)
            tmpG2 = GammaAO(nu,tau)
            tmpG3 = GammaAO(mu,gam)
            tmpG4 = GammaAO(nu,gam)

            tmpS1 = SpinAO(mu,tau)
            tmpS2 = SpinAO(nu,tau)
            tmpS3 = SpinAO(mu,gam)
            tmpS4 = SpinAO(nu,gam)

            ! s = same spin (s1=s2), d = different spin (s1/=s2)
            tmpHxcS(jj) = 0.25_dp * overSqr(mu,nu) * overSqr(tau,gam) * &
                & ( (tmpG1+tmpG2+tmpG3+tmpG4) + (tmpS1+tmpS2+tmpS3+tmpS4) )
            tmpHxcD(jj) = 0.25_dp * overSqr(mu,nu) * overSqr(tau,gam) * &
                & ( (tmpG1+tmpG2+tmpG3+tmpG4) - (tmpS1+tmpS2+tmpS3+tmpS4) )

          end if

        end do
        ! end of loop mu, nu

        ! calculate the ZmatL
        do iL = 1, LmaxR
          tmpZM(1,iL) = sum(tmpRmatL(:,iL)*tmpHxcS)
          tmpZM(2,iL) = sum(tmpRmatL(:,iL)*tmpHxcD)
        end do
        call getZmatLoop_(tmpZM, orderRmatL, Lpaired, tau, gam, ZmatL)

      end if

    end do
    ! end of loop tau, gam
!$OMP END PARALLEL DO

    ! calculate the ZmatL for LC term

    if (isRangeSep) then

      deallocate(tmpRmatL)
      deallocate(tmpHxcS)

      allocate(tmpRmatL(nOrbHalf,LmaxR))
      allocate(tmpHxcS(nOrbHalf))

      ! pack R matrix
      tmpRmatL(:,:) = 0.0_dp
      do iL = 1, LmaxR

        tmpMat(:,:) = RmatL(:,:,iL) + transpose(RmatL(:,:,iL))
        do mu = 1, nOrb
          tmpMat(mu,mu) = RmatL(mu,mu,iL)
        end do

        do ii = 1, nOrbHalf
          call getTwoIndices(nOrb, ii, mu, nu, 2)
          tmpRmatL(ii,iL) = tmpMat(mu,nu)
        end do

      end do

!$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(tau,gam, &
!$OMP& tmpHxcS,mu,nu,tmpvalue1,tmpvalue2,tmpL1, &
!$OMP& tmpL2,tmpL3,tmpL4,tmpZM,tmp22) SCHEDULE(RUNTIME)
!!$OMP& tmpL2,tmpL3,tmpL4,tmpZM) SCHEDULE(RUNTIME)
      do ii = 1, nOrbHalf

        call getTwoIndices(nOrb, ii, tau, gam, 2)
        !tmp11 = ( real(2.0_dp*nOrb+3.0_dp, dp) - sqrt( (2.0_dp*nOrb+ &
        !    & 3.0_dp)**2.0_dp - 8.0_dp*(nOrb+ii) ) )/2.0_dp
        !tau = int( real(tmp11, dp) )
        !gam = tau**2/2 - tau/2 - nOrb*tau + nOrb + ii

        tmpHxcS(:) = 0.0_dp
        do jj = 1, nOrbHalf

          ! TODO : use 'getTwoIndices' -> 5.1 sec per 1 CG-loop
          !        use tmp22 part -> 3.9 sec per 1 CG-loop
          !        Why show different cost?
          !call getTwoIndices(nOrb, jj, mu, nu, 2)
          tmp22 = ( real(2.0_dp*nOrb+3.0_dp, dp) - sqrt( (2.0_dp*nOrb+ &
              & 3.0_dp)**2.0_dp - 8.0_dp*(nOrb+jj) ) )/2.0_dp
          mu = int( real(tmp22, dp) )
          nu = mu**2/2 - mu/2 - nOrb*mu + nOrb + jj

          ! calculate the H-XC kernel for LC term
          if (isRangeSep) then

            ! (mu,nu,tau,gam)
            tmpL1 = LrGammaAO(mu,gam)
            tmpL2 = LrGammaAO(mu,nu)
            tmpL3 = LrGammaAO(tau,gam)
            tmpL4 = LrGammaAO(tau,nu)
            tmpvalue1 = -0.125_dp * overSqr(mu,tau) * &
                & overSqr(nu,gam) * (tmpL1+tmpL2+tmpL3+tmpL4)

            ! (mu,nu,gam,tau)
            tmpL1 = LrGammaAO(mu,tau)
            tmpL2 = LrGammaAO(mu,nu)
            tmpL3 = LrGammaAO(gam,tau)
            tmpL4 = LrGammaAO(gam,nu)
            tmpvalue2 = -0.125_dp * overSqr(mu,gam) * &
                & overSqr(nu,tau) * (tmpL1+tmpL2+tmpL3+tmpL4)

            tmpHxcS(jj) = tmpHxcS(jj) + tmpvalue1 + tmpvalue2

          end if

        end do
        ! end of loop mu, nu

        ! calculate the ZmatL
        do iL = 1, LmaxR
          tmpZM(1,iL) = sum(tmpRmatL(:,iL)*tmpHxcS)
          tmpZM(2,iL) = 0.0_dp
        end do
        call getZmatLoop_(tmpZM, orderRmatL, Lpaired, tau, gam, ZmatL)

      end do
      ! end of loop tau, gam
!$OMP END PARALLEL DO

    end if

    do iL = 1, Lmax
      call symmetrizeHS(ZmatL(:,:,iL))
    end do

  end subroutine getZmatNoHxc_


  !> Calculate ZmatL used in CP-REKS equations in REKS
  subroutine getZmatLoop_(tmpZM, orderRmatL, Lpaired, gam, tau, ZmatL)

    !> temporary matrix obtained from summation of RmatL and H-XC kernel
    real(dp), intent(in) :: tmpZM(:,:)

    !> Ordering between RmatL and fillingL
    integer, intent(in) :: orderRmatL(:)

    !> Number of spin-paired microstates
    integer, intent(in) :: Lpaired

    !> current indices for ZmatL
    integer, intent(in) :: gam, tau

    !> auxiliary matrix in AO basis related to SA-REKS term
    real(dp), intent(inout) :: ZmatL(:,:,:)

    integer :: iL, Lmax, tmpL1, tmpL2

    Lmax = size(ZmatL,dim=3)

    do iL = 1, Lmax
      if (iL <= Lpaired) then
        tmpL1 = iL
        tmpL2 = iL
      else
        if (mod(iL,2) == 1) then
          tmpL1 = orderRmatL(iL)
          tmpL2 = orderRmatL(iL) + 1
        else
          tmpL1 = orderRmatL(iL)
          tmpL2 = orderRmatL(iL) - 1
        end if
      end if
      ZmatL(tau,gam,iL) = ZmatL(tau,gam,iL) + tmpZM(1,tmpL1) + tmpZM(2,tmpL2)
    end do

  end subroutine getZmatLoop_


  !> Calculate Q1del used in CP-REKS equations in REKS(2,2)
  subroutine getQ1del22_(eigenvecs, Fc, Fa, FONs, SAweight, &
      & Nc, nstates, reksAlg, Q1del)

    !> Eigenvectors on eixt
    real(dp), intent(in) :: eigenvecs(:,:,:)

    !> dense fock matrix for core orbitals
    real(dp), intent(in) :: Fc(:,:)

    !> dense fock matrix for active orbitals
    real(dp), intent(in) :: Fa(:,:,:)

    !> Fractional occupation numbers of active orbitals
    real(dp), intent(in) :: FONs(:,:)

    !> Weights used in state-averaging
    real(dp), intent(in) :: SAweight(:)

    !> Number of core orbitals
    integer, intent(in) :: Nc

    !> Number of states
    integer, intent(in) :: nstates

    !> Type of REKS calculations
    integer, intent(in) :: reksAlg

    !> auxiliary matrix in AO basis related to state-interaction term
    real(dp), intent(out) :: Q1del(:,:,:)

    real(dp) :: n_a, n_b, e1, e2
    integer :: nOrb, mu, nu, t, a, b

    nOrb = size(Fc,dim=1)
    a = Nc + 1
    b = Nc + 2

    n_a = FONs(1,1)
    n_b = FONs(2,1)

    Q1del(:,:,:) = 0.0_dp
!$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(e1,e2) SCHEDULE(RUNTIME)
    do mu = 1, nOrb
      do nu = 1, nOrb

        do t = 1, nOrb

          e1 = 0.0_dp; e2 = 0.0_dp
          call assignEpsilon(Fc, Fa, SAweight, FONs, Nc, a, b, t, &
              & 2, reksAlg, e1, e2)
          Q1del(nu,mu,1) = Q1del(nu,mu,1) + ( sqrt(n_a)*e1 - sqrt(n_b)*e2 ) * &
            & eigenvecs(mu,t,1) * eigenvecs(nu,b,1) * 0.5_dp
          if (nstates == 3) then
            Q1del(nu,mu,3) = Q1del(nu,mu,3) + ( sqrt(n_a)*e1 + sqrt(n_b)*e2 ) * &
              & eigenvecs(mu,t,1) * eigenvecs(nu,b,1) * 0.5_dp
          end if

          e1 = 0.0_dp; e2 = 0.0_dp
          call assignEpsilon(Fc, Fa, SAweight, FONs, Nc, a, b, t, &
              & 1, reksAlg, e1, e2)
          Q1del(nu,mu,1) = Q1del(nu,mu,1) + ( sqrt(n_a)*e1 - sqrt(n_b)*e2 ) * &
            & eigenvecs(mu,t,1) * eigenvecs(nu,a,1) * 0.5_dp
          if (nstates == 3) then
            Q1del(nu,mu,3) = Q1del(nu,mu,3) + ( sqrt(n_a)*e1 + sqrt(n_b)*e2 ) * &
              & eigenvecs(mu,t,1) * eigenvecs(nu,a,1) * 0.5_dp
          end if

        end do

      end do
    end do
!$OMP END PARALLEL DO

  end subroutine getQ1del22_


  !> Calculate R*T contribution of gradient from SCC, spin, pc terms
  subroutine getSccPcTerms_(sccCalc, Hderiv, Sderiv, rhoSqrL, overSqr, &
      & qOutputL, q0, GammaAO, GammaDeriv, SpinAO, RmatSpL, RdelSpL, &
      & weight, getDenseAO, getAtomIndex, iSquare, orderRmatL, Lpaired, &
      & SAstates, tNAC, tExtChrg, deriv1, deriv2)

    !> SCC module internal variables
    type(TScc), allocatable, intent(inout) :: sccCalc

    !> Dense non-scc Hamiltonian derivative in AO basis
    real(dp), intent(in) :: Hderiv(:,:,:)

    !> Dense overlap derivative in AO basis
    real(dp), intent(in) :: Sderiv(:,:,:)

    !> Dense density matrix for each microstate
    real(dp), intent(inout) :: rhoSqrL(:,:,:,:)

    !> Dense overlap matrix
    real(dp), intent(in) :: overSqr(:,:)

    !> Mulliken population for each microstate
    real(dp), intent(in) :: qOutputL(:,:,:,:)

    !> reference atomic occupations
    real(dp), intent(in) :: q0(:,:,:)

    !> scc gamma integrals in AO basis
    real(dp), intent(in) :: GammaAO(:,:)

    !> scc gamma derivative integrals
    real(dp), intent(in) :: GammaDeriv(:,:,:)

    !> spin W in AO basis
    real(dp), intent(in) :: SpinAO(:,:)

    !> auxiliary matrix in AO basis related to SA-REKS term with sparse form
    real(dp), intent(in) :: RmatSpL(:,:,:)

    !> auxiliary matrix in AO basis related to state-interaction term with sparse form
    real(dp), allocatable, intent(in) :: RdelSpL(:,:,:)

    !> Weight of each microstate for state to be optimized; weight = weightL * SAweight
    real(dp), intent(in) :: weight(:)

    !> get dense AO index from sparse AO array
    integer, intent(in) :: getDenseAO(:,:)

    !> get atom index from AO index
    integer, intent(in) :: getAtomIndex(:)

    !> Position of each atom in the rows/columns of the square matrices. Shape: (nAtom)
    integer, intent(in) :: iSquare(:)

    !> Ordering between RmatL and fillingL
    integer, intent(in) :: orderRmatL(:)

    !> Number of spin-paired microstates
    integer, intent(in) :: Lpaired

    !> Number of states used in state-averaging
    integer, intent(in) :: SAstates

    !> Calculate nonadiabatic coupling vectors
    logical, intent(in) :: tNAC

    !> If external charges must be considered
    logical, intent(in) :: tExtChrg

    !> computed tr(R*T) gradient for SA-REKS, SSR, or L state
    real(dp), intent(inout) :: deriv1(:,:,:)

    !> computed tr(R*T) gradient for state-interaction term
    real(dp), intent(inout) :: deriv2(:,:,:)

    real(dp), allocatable :: q0AO(:)                 ! nOrb
    real(dp), allocatable :: qOutputAO(:,:)          ! nOrb, Lmax; up+down
    real(dp), allocatable :: qSpinAO(:,:)            ! nOrb, Lmax; up-down
    real(dp), allocatable :: GammaQ(:,:)             ! nOrb, Lmax
    real(dp), allocatable :: SpinQ(:,:)              ! nOrb, Lmax
    real(dp), allocatable :: GammaDerivQ(:,:,:,:)    ! nAtom, nAtom, 3, Lmax

    real(dp), allocatable :: QinvR(:)                ! nAtom

    real(dp), allocatable :: tmpS(:,:,:)             ! mOrb, mOrb, 3
    real(dp), allocatable :: tmpD(:,:,:,:)           ! mOrb, mOrb, nSpin, Lmax
    real(dp), allocatable :: tmpGM(:,:,:)            ! nOrb, mOrb, 2
    real(dp), allocatable :: tmpPM(:,:,:)            ! nOrb, mOrb, 2
    real(dp), allocatable :: tmpQ1(:,:,:)            ! mOrb, 1, nSpin
    real(dp), allocatable :: tmpQ2(:,:,:)            ! mOrb, 1, nSpin

    real(dp), allocatable :: GammaQderiv(:,:,:,:)    ! nOrb, 1, 3, Lmax
    real(dp), allocatable :: SpinQderiv(:,:,:,:)     ! nOrb, 1, 3, Lmax
    real(dp), allocatable :: TderivL(:,:,:,:)        ! sparseSize, 3, Lmax, Ncpu

    real(dp) :: tmpCoulomb, tmpG1, tmpS1, tmpG2(3), tmpV1(3)
    integer :: iAtom1, iAtom2, iAtom3, iAtom4, nAtom, k, nAtomPair
    integer :: ist, nstates, nstHalf, mOrb, mu, nu, nOrb, l, sparseSize
    integer :: iS, nSpin, iL, Lmax, id, Ncpu
    integer :: ii, G1, G2, G3, G4

    nAtom = size(GammaDeriv,dim=1)
    nAtomPair = nAtom * (nAtom - 1) / 2
    nstates = size(RmatSpL,dim=3)
    nstHalf = nstates * (nstates - 1) / 2

    mOrb = size(qOutputL,dim=1)
    nOrb = size(GammaAO,dim=1)
    sparseSize = size(RmatSpL,dim=1)

    nSpin = size(qOutputL,dim=3)
    Lmax = size(qOutputL,dim=4)
  #:if WITH_OMP
!$OMP PARALLEL
    Ncpu = OMP_GET_NUM_THREADS()
!$OMP END PARALLEL
  #:else
    Ncpu = 1
  #:endif

    allocate(q0AO(nOrb))
    allocate(qOutputAO(nOrb,Lmax))
    allocate(qSpinAO(nOrb,Lmax))
    allocate(GammaQ(nOrb,Lmax))
    allocate(SpinQ(nOrb,Lmax))
    allocate(GammaDerivQ(nAtom,nAtom,3,Lmax))

    if (tExtChrg) then
      allocate(QinvR(nAtom))
    end if

    allocate(tmpS(mOrb,mOrb,3))
    allocate(tmpD(mOrb,mOrb,nSpin,Lmax))
    allocate(tmpGM(nOrb,mOrb,2))
    allocate(tmpPM(nOrb,mOrb,2))
    allocate(tmpQ1(mOrb,1,nSpin))
    allocate(tmpQ2(mOrb,1,nSpin))

    allocate(GammaQderiv(nOrb,1,3,Lmax))
    allocate(SpinQderiv(nOrb,1,3,Lmax))
    allocate(TderivL(sparseSize,3,Lmax,Ncpu))

    ! get q0/qOutput/qSpin with respect to AO
    call getq0qOutputqSpin(qOutputL, q0, iSquare, qOutputAO, qSpinAO, q0AO)
    ! rearrange gamma*Q, spin*Q
    call getGammaQSpinQ(qOutputAO, qSpinAO, q0AO, GammaAO, SpinAO, GammaQ, SpinQ)
    ! rearrange (gamma derivtive) * Q
    call getGammaDerivQ(qOutputL, q0, GammaDeriv, GammaDerivQ)

    ! contributions related to gradient of external charges
    ! get Q(PC) * 1/R(between QM and PC)
    if (tExtChrg) then
      call sccCalc%getShiftOfPC(QinvR)
    end if

    ! compute R*T shift with only up-spin part of TderivL due to symmetry

!$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(id,iAtom1,iAtom2,G1, &
!$OMP& G2,G3,G4,tmpS,tmpD,tmpGM,tmpPM,tmpQ1,tmpQ2,GammaQderiv, &
!$OMP& SpinQderiv,mu,nu,iAtom3,iAtom4,tmpV1,tmpG2,tmpG1,tmpS1, &
!$OMP& tmpCoulomb,ist) REDUCTION(+:deriv1,deriv2) SCHEDULE(RUNTIME)
    do k = 1, nAtomPair

    #:if WITH_OMP
      id = OMP_GET_THREAD_NUM() + 1
    #:else
      id = 1
    #:endif

      call getTwoIndices(nAtom, k, iAtom1, iAtom2, 1)

      G1 = iSquare(iAtom1)
      G2 = iSquare(iAtom1+1)-1
      G3 = iSquare(iAtom2)
      G4 = iSquare(iAtom2+1)-1

      ! zeroing for temporary S_deriv
      tmpS(:,:,:) = 0.0_dp
      do ii = 1, 3
        tmpS(1:G2-G1+1,1:G4-G3+1,ii) = Sderiv(G1:G2,G3:G4,ii)
      end do

      ! zeroing for temporary density matrix
      ! rhoSqrL has (my_qm) component
      ! tmpD has (my_qm) component
      tmpD(:,:,:,:) = 0.0_dp
      do iL = 1, Lmax
        tmpD(1:G2-G1+1,1:G4-G3+1,1,iL) = rhoSqrL(G1:G2,G3:G4,1,iL)
      end do
      ! tmpD has (qm) component
      call qmExpandL(tmpD, Lpaired)

      ! zeroing for temporary gamma & spin parts
      tmpGM(:,:,:) = 0.0_dp
      tmpGM(:,1:G2-G1+1,1) = GammaAO(:,G1:G2)
      tmpGM(:,1:G4-G3+1,2) = GammaAO(:,G3:G4)
      tmpPM(:,:,:) = 0.0_dp
      tmpPM(:,1:G2-G1+1,1) = SpinAO(:,G1:G2)
      tmpPM(:,1:G4-G3+1,2) = SpinAO(:,G3:G4)

      ! calculate the charge derivative part
      do iL = 1, Lmax
        do ii = 1, 3

          ! zeroing for tmpQ1 & tmpQ2 with each k pair
          ! tmpQ1 & tmpQ2 has (qm) component
          tmpQ1(:,:,:) = 0.0_dp; tmpQ2(:,:,:) = 0.0_dp
          do iS = 1, nSpin
            ! G3 ~ G4
            tmpQ1(:,1,iS) = sum(tmpS(:,:,ii)*tmpD(:,:,iS,iL),dim=1)
            ! G1 ~ G2
            tmpQ2(:,1,iS) = sum(tmpS(:,:,ii)*tmpD(:,:,iS,iL),dim=2)
          end do

          ! rearrange gamma * (Q derivative), spin * (Q derivative)
          GammaQderiv(:,:,ii,iL) = matmul(tmpGM(:,:,1),tmpQ2(:,:,1)) &
              & + matmul(tmpGM(:,:,2),tmpQ1(:,:,1))
          SpinQderiv(:,:,ii,iL) = matmul(tmpPM(:,:,1),tmpQ2(:,:,2)) &
              & + matmul(tmpPM(:,:,2),tmpQ1(:,:,2))

        end do
      end do

      ! zeroing for temporary TderivL in each k atom pair
      ! calculate sparse TderivL in AO basis
      TderivL(:,:,:,id) = 0.0_dp
      do l = 1, sparseSize

        ! set the AO indices with respect to sparsity
        mu = getDenseAO(l,1)
        nu = getDenseAO(l,2)
        ! find proper atom index
        iAtom3 = getAtomIndex(mu)
        iAtom4 = getAtomIndex(nu)

        if (mu <= nu) then
          do iL = 1, Lmax

            ! calculate charge derivative part (overlap derivative in charge)
            ! this is non-zero for all (mu,nu) regardless of R_(alpha,beta) pair
            if (abs(overSqr(mu,nu)) >= epsilon(1.0_dp)) then
              tmpV1(:) = 0.0_dp
              tmpV1(:) = tmpV1 + (GammaQderiv(mu,1,:,iL) + GammaQderiv(nu,1,:,iL))
              tmpV1(:) = tmpV1 + (SpinQderiv(mu,1,:,iL) + SpinQderiv(nu,1,:,iL))
              TderivL(l,:,iL,id) = TderivL(l,:,iL,id) + 0.5_dp*overSqr(mu,nu)*tmpV1
            end if

            ! calculate (gamma derivative) part
            if (abs(overSqr(mu,nu)) >= epsilon(1.0_dp)) then
              tmpG2(:) = 0.0_dp
              ! mu in alpha
              if (iAtom3 == iAtom1) then
                tmpG2(:) = tmpG2 + GammaDerivQ(iAtom2,iAtom1,:,iL)
              end if
              ! mu in beta
              if (iAtom3 == iAtom2) then
                tmpG2(:) = tmpG2 + GammaDerivQ(iAtom1,iAtom2,:,iL)
              end if
              ! nu in alpha
              if (iAtom4 == iAtom1) then
                tmpG2(:) = tmpG2 + GammaDerivQ(iAtom2,iAtom1,:,iL)
              end if
              ! nu in beta
              if (iAtom4 == iAtom2) then
                tmpG2(:) = tmpG2 + GammaDerivQ(iAtom1,iAtom2,:,iL)
              end if
              TderivL(l,:,iL,id) = TderivL(l,:,iL,id) + 0.5_dp*overSqr(mu,nu)*tmpG2
            end if

            ! calculate Hderiv, Sderiv for scc/spin part
            ! mu in alpha and nu in beta
            if (iAtom3 == iAtom1 .and. iAtom4 == iAtom2) then
              tmpG1 = GammaQ(mu,iL) + GammaQ(nu,iL)
              tmpS1 = SpinQ(mu,iL) + SpinQ(nu,iL)
              TderivL(l,:,iL,id) = TderivL(l,:,iL,id) + Hderiv(mu,nu,:) &
                  & + 0.5_dp*Sderiv(mu,nu,:)*tmpG1 + 0.5_dp*Sderiv(mu,nu,:)*tmpS1
              if (tExtChrg) then
                tmpCoulomb = QinvR(iAtom1) + QinvR(iAtom2)
                TderivL(l,:,iL,id) = TderivL(l,:,iL,id) &
                    & + 0.5_dp*Sderiv(mu,nu,:)*tmpCoulomb
              end if
            end if

          end do
          ! end of loop iL
        end if

      end do
      ! end of loop l

      if (tNAC) then
        do ist = 1, nstates
          if (ist /= SAstates) then
            call shiftRTgradSparse_(deriv1(:,:,ist), RmatSpL(:,:,ist), &
                & TderivL(:,:,:,id), weight, orderRmatL, iAtom1, iAtom2)
          end if
        end do
        do ist = 1, nstHalf
          call shiftRTgradSparse_(deriv2(:,:,ist), RdelSpL(:,:,ist), &
              & TderivL(:,:,:,id), weight, orderRmatL, iAtom1, iAtom2)
        end do
      else
        call shiftRTgradSparse_(deriv1(:,:,1), RmatSpL(:,:,1), &
            & TderivL(:,:,:,id), weight, orderRmatL, iAtom1, iAtom2)
      end if

    end do
    ! end of loop k
!$OMP END PARALLEL DO

    contains

      !> get q0/qOutput/qSpin with respect to AO
      subroutine getq0qOutputqSpin(qOutputL, q0, iSquare, &
          & qOutputAO, qSpinAO, q0AO)

        !> Mulliken population for each microstate
        real(dp), intent(in) :: qOutputL(:,:,:,:)

        !> reference atomic occupations
        real(dp), intent(in) :: q0(:,:,:)

        !> Position of each atom in the rows/columns of the square matrices. Shape: (nAtom)
        integer, intent(in) :: iSquare(:)

        !> mulliken population with up+down in AO basis
        real(dp), intent(out) :: qOutputAO(:,:)

        !> mulliken population with up-down in AO basis
        real(dp), intent(out) :: qSpinAO(:,:)

        !> mulliken population with reference in AO basis
        real(dp), intent(out) :: q0AO(:)

        integer :: nAtom, Lmax, iAt, i, iL

        nAtom = size(qOutputL,dim=2)
        Lmax = size(qOutputL,dim=4)

        q0AO(:) = 0.0_dp
        qOutputAO(:,:) = 0.0_dp
        qSpinAO(:,:) = 0.0_dp
        do iAt = 1, nAtom
          i = iSquare(iAt+1) - iSquare(iAt)
          q0AO(iSquare(iAt):iSquare(iAt+1)-1) = q0(1:i,iAt,1)
          do iL = 1, Lmax
            qOutputAO(iSquare(iAt):iSquare(iAt+1)-1,iL) = qOutputL(1:i,iAt,1,iL)
            qSpinAO(iSquare(iAt):iSquare(iAt+1)-1,iL) = qOutputL(1:i,iAt,2,iL)
          end do
        end do

      end subroutine getq0qOutputqSpin

      !> rearrange gamma*Q, spin*Q
      subroutine getGammaQSpinQ(qOutputAO, qSpinAO, q0AO, &
          & GammaAO, SpinAO, GammaQ, SpinQ)

        !> mulliken population with up+down in AO basis
        real(dp), intent(in) :: qOutputAO(:,:)

        !> mulliken population with up-down in AO basis
        real(dp), intent(in) :: qSpinAO(:,:)

        !> mulliken population with reference in AO basis
        real(dp), intent(in) :: q0AO(:)

        !> scc gamma integrals in AO basis
        real(dp), intent(in) :: GammaAO(:,:)

        !> spin W in AO basis
        real(dp), intent(in) :: SpinAO(:,:)

        !> (gamma) * (mulliken population) in AO basis
        real(dp), intent(out) :: GammaQ(:,:)

        !> (spin constants) * (spin population) in AO basis
        real(dp), intent(out) :: SpinQ(:,:)

        integer :: nOrb, Lmax, iL, mu

        nOrb = size(GammaAO,dim=1)
        Lmax = size(qOutputAO,dim=2)

        GammaQ(:,:) = 0.0_dp
        SpinQ(:,:) = 0.0_dp
        do iL = 1, Lmax
          do mu = 1, nOrb
            GammaQ(mu,iL) = sum(GammaAO(:,mu)*(qOutputAO(:,iL)-q0AO(:)),dim=1)
            SpinQ(mu,iL) = sum(SpinAO(:,mu)*qSpinAO(:,iL),dim=1)
          end do
        end do

      end subroutine getGammaQSpinQ

      !> rearrange (gamma derivtive) * Q
      subroutine getGammaDerivQ(qOutputL, q0, GammaDeriv, GammaDerivQ)

        !> Mulliken population for each microstate
        real(dp), intent(in) :: qOutputL(:,:,:,:)

        !> reference atomic occupations
        real(dp), intent(in) :: q0(:,:,:)

        !> scc gamma derivative integrals
        real(dp), intent(in) :: GammaDeriv(:,:,:)

        !> (gamma deriv) * (mulliken population) in AO basis
        real(dp), intent(out) :: GammaDerivQ(:,:,:,:)

        integer :: nAtom, Lmax, iL, iAt1, iAt2

        Lmax = size(qOutputL,dim=4)
        nAtom = size(qOutputL,dim=2)

        GammaDerivQ(:,:,:,:) = 0.0_dp
        do iL = 1, Lmax
          do iAt1 = 1, nAtom
            do iAt2 = 1, nAtom
              GammaDerivQ(iAt1,iAt2,:,iL) = GammaDeriv(iAt1,iAt2,:) * &
                   & sum(qOutputL(:,iAt1,1,iL) - q0(:,iAt1,1),dim=1)
            end do
          end do
        end do

      end subroutine getGammaDerivQ

  end subroutine getSccPcTerms_


  !> Calculate R*T contribution of gradient from LC term with overlap derivative
  subroutine getLr1stTerms_(Sderiv, deltaRhoSqrL, overSqr, LrGammaAO, SP, &
      & RmatHalfL, RdelHalfL, weight, getDenseAtom, iSquare, orderRmatL, &
      & SAstates, mOrb, tNAC, deriv1, deriv2)

    !> Dense overlap derivative in AO basis
    real(dp), intent(in) :: Sderiv(:,:,:)

    !> Dense delta density matrix for each microstate
    real(dp), intent(in) :: deltaRhoSqrL(:,:,:,:)

    !> Dense overlap matrix
    real(dp), intent(in) :: overSqr(:,:)

    !> long-range gamma integrals in AO basis
    real(dp), intent(in) :: LrGammaAO(:,:)

    !> Dense overlap * delta density in AO basis
    real(dp), intent(in) :: SP(:,:,:)

    !> auxiliary matrix in AO basis related to SA-REKS term with half dense form
    real(dp), intent(in) :: RmatHalfL(:,:,:)

    !> auxiliary matrix in AO basis related to state-interaction term with half dense form
    real(dp), allocatable, intent(in) :: RdelHalfL(:,:,:)

    !> Weight of each microstate for state to be optimized; weight = weightL * SAweight
    real(dp), intent(in) :: weight(:)

    !> get dense atom index from sparse atom array
    integer, intent(in) :: getDenseAtom(:,:)

    !> Position of each atom in the rows/columns of the square matrices. Shape: (nAtom)
    integer, intent(in) :: iSquare(:)

    !> Ordering between RmatL and fillingL
    integer, intent(in) :: orderRmatL(:)

    !> Number of states used in state-averaging
    integer, intent(in) :: SAstates

    !> Max. nr. of orbitals for any species
    integer, intent(in) :: mOrb

    !> Calculate nonadiabatic coupling vectors
    logical, intent(in) :: tNAC

    !> computed tr(R*T) gradient for SA-REKS, SSR, or L state
    real(dp), intent(inout) :: deriv1(:,:,:)

    !> computed tr(R*T) gradient for state-interaction term
    real(dp), intent(inout) :: deriv2(:,:,:)

    real(dp), allocatable :: tmpS(:,:,:)            ! mOrb, mOrb, 3

    real(dp), allocatable :: shiftPP1(:,:)          ! mOrb, nOrb
    real(dp), allocatable :: shiftIM1(:,:)          ! mOrb, nOrb
    real(dp), allocatable :: shiftFM1(:,:,:,:)      ! mOrb, nOrb, 3, Lmax

    real(dp), allocatable :: shiftPP2(:,:)          ! nOrb, mOrb
    real(dp), allocatable :: shiftIM2(:,:)          ! nOrb, mOrb
    real(dp), allocatable :: shiftFM2(:,:,:,:)      ! nOrb, mOrb, 3, Lmax

    real(dp), allocatable :: TderivL(:,:,:,:)       ! nOrbHalf, 3, Lmax, Ncpu

    integer :: iAtom1, iAtom2, k, nAtomSparse, ist, nstates, nstHalf
    integer :: mu, nu, al, be, nOrb, l, nOrbHalf
    integer :: iL, Lmax, id, Ncpu
    integer :: ii, G1, G2, G3, G4

    nAtomSparse = size(getDenseAtom,dim=1)
    nstates = size(RmatHalfL,dim=3)
    nstHalf = nstates * (nstates - 1) /2

    nOrb = size(deltaRhoSqrL,dim=1)
    nOrbHalf = size(RmatHalfL,dim=1)

    Lmax = size(deltaRhoSqrL,dim=4)
  #:if WITH_OMP
!$OMP PARALLEL
    Ncpu = OMP_GET_NUM_THREADS()
!$OMP END PARALLEL
  #:else
    Ncpu = 1
  #:endif

    allocate(tmpS(mOrb,mOrb,3))

    allocate(shiftPP1(mOrb,nOrb))
    allocate(shiftIM1(mOrb,nOrb))
    allocate(shiftFM1(mOrb,nOrb,3,Lmax))

    allocate(shiftPP2(nOrb,mOrb))
    allocate(shiftIM2(nOrb,mOrb))
    allocate(shiftFM2(nOrb,mOrb,3,Lmax))

    allocate(TderivL(nOrbHalf,3,Lmax,Ncpu))

    ! compute R*T shift with only up-spin part of TderivL due to symmetry

!$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(id,iAtom1,iAtom2, &
!$OMP& G1,G2,G3,G4,tmpS,shiftPP1,shiftPP2,shiftIM1,shiftIM2, &
!$OMP& shiftFM1,shiftFM2,iL,ii,mu,nu,al,be,l,ist) &
!$OMP& REDUCTION(+:deriv1,deriv2) SCHEDULE(RUNTIME)
    loopKK: do k = 1, nAtomSparse

    #:if WITH_OMP
      id = OMP_GET_THREAD_NUM() + 1
    #:else
      id = 1
    #:endif

      ! set the atom indices with respect to sparsity
      iAtom1 = getDenseAtom(k,1)
      iAtom2 = getDenseAtom(k,2)

      if (iAtom1 /= iAtom2) then

        G1 = iSquare(iAtom1)
        G2 = iSquare(iAtom1+1)-1
        G3 = iSquare(iAtom2)
        G4 = iSquare(iAtom2+1)-1

        ! zeroing for temporary Sderiv
        tmpS(:,:,:) = 0.0_dp
        do ii = 1, 3
          tmpS(1:G2-G1+1,1:G4-G3+1,ii) = Sderiv(G1:G2,G3:G4,ii)
        end do

        ! zeroing for temporary shift
        shiftFM1(:,:,:,:) = 0.0_dp
        shiftFM2(:,:,:,:) = 0.0_dp

        ! a-1
        do iL = 1, Lmax
          shiftPP1(:,:) = 0.0_dp
          shiftPP1(1:G4-G3+1,:) = deltaRhoSqrL(G3:G4,:,1,iL)
          do ii = 1, 3
            shiftIM1(:,:) = 0.0_dp
            call gemm(shiftIM1, tmpS(:,:,ii), shiftPP1)
            do mu = G1, G2
              do be = 1, nOrb
                shiftIM1(mu-G1+1,be) = shiftIM1(mu-G1+1,be) * LrGammaAO(mu,be)
              end do
            end do
            call gemm(shiftFM1(:,:,ii,iL), shiftIM1, overSqr)
          end do
        end do
        ! b-1
        do iL = 1, Lmax
          shiftPP1(:,:) = 0.0_dp
          shiftPP1(1:G4-G3+1,:) = transpose(SP(:,G3:G4,iL))
          do ii = 1, 3
            shiftIM1(:,:) = 0.0_dp
            call gemm(shiftIM1, tmpS(:,:,ii), shiftPP1)
            do mu = G1, G2
              do nu = 1, nOrb
                shiftFM1(mu-G1+1,nu,ii,iL) = shiftFM1(mu-G1+1,nu,ii,iL) &
                    & + shiftIM1(mu-G1+1,nu) * LrGammaAO(mu,nu)
              end do
            end do
          end do
        end do
        ! c-1
        do iL = 1, Lmax
          shiftPP1(:,:) = 0.0_dp
          do al = G3, G4
            do be = 1, nOrb
              shiftPP1(al-G3+1,be) = deltaRhoSqrL(al,be,1,iL) * LrGammaAO(al,be)
            end do
          end do
          shiftIM1(:,:) = 0.0_dp
          call gemm(shiftIM1, shiftPP1, overSqr)
          do ii = 1, 3
            call gemm(shiftFM1(:,:,ii,iL), tmpS(:,:,ii), shiftIM1, alpha=1.0_dp, beta=1.0_dp)
          end do
        end do
        ! d-1
        do iL = 1, Lmax
          shiftIM1(:,:) = 0.0_dp
          do al = G3, G4
            do nu = 1, nOrb
              shiftIM1(al-G3+1,nu) = SP(nu,al,iL) * LrGammaAO(al,nu)
            end do
          end do
          do ii = 1, 3
            call gemm(shiftFM1(:,:,ii,iL), tmpS(:,:,ii), shiftIM1, alpha=1.0_dp, beta=1.0_dp)
          end do
        end do

        ! a-2
        do iL = 1, Lmax
          shiftIM2(:,:) = 0.0_dp
          do mu = 1, nOrb
            do be = G1, G2
              shiftIM2(mu,be-G1+1) = SP(mu,be,iL) * LrGammaAO(mu,be)
            end do
          end do
          do ii = 1, 3
            call gemm(shiftFM2(:,:,ii,iL), shiftIM2, tmpS(:,:,ii))
          end do
        end do
        ! b-2
        do iL = 1, Lmax
          shiftPP2(:,:) = 0.0_dp
          shiftPP2(:,1:G2-G1+1) = SP(:,G1:G2,iL)
          do ii = 1, 3
            shiftIM2(:,:) = 0.0_dp
            call gemm(shiftIM2, shiftPP2, tmpS(:,:,ii))
            do mu = 1, nOrb
              do nu = G3, G4
                shiftFM2(mu,nu-G3+1,ii,iL) = shiftFM2(mu,nu-G3+1,ii,iL) &
                    & + shiftIM2(mu,nu-G3+1) * LrGammaAO(mu,nu)
              end do
            end do
          end do
        end do
        ! c-2
        do iL = 1, Lmax
          shiftPP2(:,:) = 0.0_dp
          do al = 1, nOrb
            do be = G1, G2
              shiftPP2(al,be-G1+1) = deltaRhoSqrL(al,be,1,iL) * LrGammaAO(al,be)
            end do
          end do
          shiftIM2(:,:) = 0.0_dp
          call gemm(shiftIM2, overSqr, shiftPP2)
          do ii = 1, 3
            call gemm(shiftFM2(:,:,ii,iL), shiftIM2, tmpS(:,:,ii), alpha=1.0_dp, beta=1.0_dp)
          end do
        end do
        ! d-2
        do iL = 1, Lmax
          shiftPP2(:,:) = 0.0_dp
          shiftPP2(:,1:G2-G1+1) = deltaRhoSqrL(:,G1:G2,1,iL)
          do ii = 1, 3
            shiftIM2(:,:) = 0.0_dp
            call gemm(shiftIM2, shiftPP2, tmpS(:,:,ii))
            do al = 1, nOrb
              do nu = G3, G4
                shiftIM2(al,nu-G3+1) = shiftIM2(al,nu-G3+1) * LrGammaAO(al,nu)
              end do
            end do
            call gemm(shiftFM2(:,:,ii,iL), overSqr, shiftIM2, alpha=1.0_dp, beta=1.0_dp)
          end do
        end do

        ! zeroing for temporary TderivL in each k atom pair
        ! calculate half dense TderivL in AO basis
        TderivL(:,:,:,id) = 0.0_dp
        loopLL: do l = 1, nOrbHalf

          call getTwoIndices(nOrb, l, mu, nu, 2)

          if (mu >= G1 .and. mu <= G2) then
            TderivL(l,:,:,id) = TderivL(l,:,:,id) + shiftFM1(mu-G1+1,nu,:,:)
          end if
          if (nu >= G1 .and. nu <= G2) then
            TderivL(l,:,:,id) = TderivL(l,:,:,id) + shiftFM1(nu-G1+1,mu,:,:)
          end if
          if (mu >= G3 .and. mu <= G4) then
            TderivL(l,:,:,id) = TderivL(l,:,:,id) + shiftFM2(nu,mu-G3+1,:,:)
          end if
          if (nu >= G3 .and. nu <= G4) then
            TderivL(l,:,:,id) = TderivL(l,:,:,id) + shiftFM2(mu,nu-G3+1,:,:)
          end if

        end do loopLL
        TderivL(:,:,:,id) = TderivL(:,:,:,id) * (-0.25_dp)

        if (tNAC) then
          do ist = 1, nstates
            if (ist /= SAstates) then
              call shiftRTgradSparse_(deriv1(:,:,ist), RmatHalfL(:,:,ist), &
                  & TderivL(:,:,:,id), weight, orderRmatL, iAtom1, iAtom2)
            end if
          end do
          do ist = 1, nstHalf
            call shiftRTgradSparse_(deriv2(:,:,ist), RdelHalfL(:,:,ist), &
                & TderivL(:,:,:,id), weight, orderRmatL, iAtom1, iAtom2)
          end do
        else
          call shiftRTgradSparse_(deriv1(:,:,1), RmatHalfL(:,:,1), &
              & TderivL(:,:,:,id), weight, orderRmatL, iAtom1, iAtom2)
        end if

      end if

    end do loopKK
!$OMP END PARALLEL DO

  end subroutine getLr1stTerms_


  !> Calculate R*T contribution of gradient from LC term with long-range gamma derivative
  subroutine getLr2ndTerms_(deltaRhoSqrL, overSqr, LrGammaDeriv, SP, SPS, &
      & RmatL, RdelL, tmpRL, weight, iSquare, orderRmatL, SAstates, mOrb, &
      & tNAC, deriv1, deriv2)

    !> Dense delta density matrix for each microstate
    real(dp), intent(in) :: deltaRhoSqrL(:,:,:,:)

    !> Dense overlap matrix
    real(dp), intent(in) :: overSqr(:,:)

    !> long-range gamma derivative integrals
    real(dp), intent(in) :: LrGammaDeriv(:,:,:)

    !> Dense overlap * delta density in AO basis
    real(dp), intent(in) :: SP(:,:,:)

    !> Dense overlap * delta density * overlap in AO basis
    real(dp), intent(in) :: SPS(:,:,:)

    !> auxiliary matrix in AO basis related to SA-REKS term
    real(dp), intent(in) :: RmatL(:,:,:,:)

    !> auxiliary matrix in AO basis related to state-interaction term
    real(dp), allocatable, intent(in) :: RdelL(:,:,:,:)

    !> auxiliary matrix in AO basis related to state-interaction term
    real(dp), allocatable, intent(in) :: tmpRL(:,:,:,:)

    !> Weight of each microstate for state to be optimized; weight = weightL * SAweight
    real(dp), intent(in) :: weight(:)

    !> Position of each atom in the rows/columns of the square matrices. Shape: (nAtom)
    integer, intent(in) :: iSquare(:)

    !> Ordering between RmatL and fillingL
    integer, intent(in) :: orderRmatL(:)

    !> Number of states used in state-averaging
    integer, intent(in) :: SAstates

    !> Max. nr. of orbitals for any species
    integer, intent(in) :: mOrb

    !> Calculate nonadiabatic coupling vectors
    logical, intent(in) :: tNAC

    !> computed tr(R*T) gradient for SA-REKS, SSR, or L state
    real(dp), intent(inout) :: deriv1(:,:,:)

    !> computed tr(R*T) gradient for state-interaction term
    real(dp), intent(inout) :: deriv2(:,:,:)

    real(dp), allocatable :: SR(:,:,:)
    real(dp), allocatable :: RS(:,:,:)
    real(dp), allocatable :: SRS(:,:,:)

    integer :: iAtom1, iAtom2, nAtom, k, nAtomPair
    integer :: ist, nstates, nstHalf, nOrb, LmaxR

    nAtom = size(LrGammaDeriv,dim=1)
    nAtomPair = nAtom * (nAtom - 1) / 2
    nstates = size(RmatL,dim=4)
    nstHalf = nstates * (nstates - 1) / 2
    nOrb = size(overSqr,dim=1)
    LmaxR = size(RmatL,dim=3)

    allocate(SR(nOrb,nOrb,LmaxR))
    allocate(RS(nOrb,nOrb,LmaxR))
    allocate(SRS(nOrb,nOrb,LmaxR))

    if (tNAC) then

      do ist = 1, nstates
        if (ist /= SAstates) then

          call getSRmatrices(RmatL(:,:,:,ist), overSqr, SR, RS, SRS)

!$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(iAtom1,iAtom2) &
!$OMP& REDUCTION(+:deriv1) SCHEDULE(RUNTIME)
          do k = 1, nAtomPair

            call getTwoIndices(nAtom, k, iAtom1, iAtom2, 1)

            call shiftRTgradLr1st_(deriv1(:,:,ist), RmatL(:,:,:,ist), SPS, LrGammaDeriv, &
                & weight, iSquare, orderRmatL, mOrb, iAtom1, iAtom2)
            call shiftRTgradLr1st_(deriv1(:,:,ist), SRS, deltaRhoSqrL(:,:,1,:), LrGammaDeriv, &
                & weight, iSquare, orderRmatL, mOrb, iAtom1, iAtom2)
            call shiftRTgradLr2nd_(deriv1(:,:,ist), RS, SP, LrGammaDeriv, &
                & weight, iSquare, orderRmatL, mOrb, iAtom1, iAtom2, 1)
            call shiftRTgradLr2nd_(deriv1(:,:,ist), SR, SP, LrGammaDeriv, &
                & weight, iSquare, orderRmatL, mOrb, iAtom1, iAtom2, 2)

          end do
!$OMP END PARALLEL DO

        end if
      end do

      do ist = 1, nstHalf

        call getSRmatrices(RdelL(:,:,:,ist) + tmpRL(:,:,:,ist), &
            & overSqr, SR, RS, SRS)

!$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(iAtom1,iAtom2) &
!$OMP& REDUCTION(+:deriv2) SCHEDULE(RUNTIME)
        do k = 1, nAtomPair

          call getTwoIndices(nAtom, k, iAtom1, iAtom2, 1)

          call shiftRTgradLr1st_(deriv2(:,:,ist), RdelL(:,:,:,ist) + &
              & tmpRL(:,:,:,ist), SPS, LrGammaDeriv, weight, iSquare, &
              & orderRmatL, mOrb, iAtom1, iAtom2)
          call shiftRTgradLr1st_(deriv2(:,:,ist), SRS, deltaRhoSqrL(:,:,1,:), LrGammaDeriv, &
              & weight, iSquare, orderRmatL, mOrb, iAtom1, iAtom2)
          call shiftRTgradLr2nd_(deriv2(:,:,ist), RS, SP, LrGammaDeriv, &
              & weight, iSquare, orderRmatL, mOrb, iAtom1, iAtom2, 1)
          call shiftRTgradLr2nd_(deriv2(:,:,ist), SR, SP, LrGammaDeriv, &
              & weight, iSquare, orderRmatL, mOrb, iAtom1, iAtom2, 2)

        end do
!$OMP END PARALLEL DO
      end do

    else

      call getSRmatrices(RmatL(:,:,:,1), overSqr, SR, RS, SRS)

!$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(iAtom1,iAtom2) &
!$OMP& REDUCTION(+:deriv1) SCHEDULE(RUNTIME)
      do k = 1, nAtomPair
  
        call getTwoIndices(nAtom, k, iAtom1, iAtom2, 1)
  
        call shiftRTgradLr1st_(deriv1(:,:,1), RmatL(:,:,:,1), SPS, LrGammaDeriv, &
            & weight, iSquare, orderRmatL, mOrb, iAtom1, iAtom2)
        call shiftRTgradLr1st_(deriv1(:,:,1), SRS, deltaRhoSqrL(:,:,1,:), LrGammaDeriv, &
            & weight, iSquare, orderRmatL, mOrb, iAtom1, iAtom2)
        call shiftRTgradLr2nd_(deriv1(:,:,1), RS, SP, LrGammaDeriv, &
            & weight, iSquare, orderRmatL, mOrb, iAtom1, iAtom2, 1)
        call shiftRTgradLr2nd_(deriv1(:,:,1), SR, SP, LrGammaDeriv, &
            & weight, iSquare, orderRmatL, mOrb, iAtom1, iAtom2, 2)
  
      end do
!$OMP END PARALLEL DO

    end if

    contains

      !> Calculate matrix product of overlap and R matrix
      subroutine getSRmatrices(RmatL, overSqr, SR, RS, SRS)

        !> auxiliary matrix in AO basis related to SA-REKS term
        real(dp), intent(in) :: RmatL(:,:,:)

        !> Dense overlap matrix
        real(dp), intent(in) :: overSqr(:,:)

        !> Dense overlap * R matrix in AO basis
        real(dp), intent(out) :: SR(:,:,:)

        !> Dense R matrix * overlap in AO basis
        real(dp), intent(out) :: RS(:,:,:)

        !> Dense overlap * R matrix * overlap in AO basis
        real(dp), intent(out) :: SRS(:,:,:)

        real(dp), allocatable :: tmpMat(:,:)
        integer :: iL, LmaxR, nOrb

        nOrb = size(overSqr,dim=1)
        LmaxR = size(RmatL,dim=3)

        allocate(tmpMat(nOrb,nOrb))

        SR(:,:,:) = 0.0_dp
        do iL = 1, LmaxR
          ! S * R calculation
          call gemm(SR(:,:,iL), overSqr, RmatL(:,:,iL))
        end do
        RS(:,:,:) = 0.0_dp
        SRS(:,:,:) = 0.0_dp
        do iL = 1, LmaxR
          tmpMat(:,:) = 0.0_dp
          call gemm(tmpMat, RmatL(:,:,iL), overSqr)
          ! R * S calculation
          RS(:,:,iL) = tmpMat
          ! S * R * S calculation
          call gemm(SRS(:,:,iL), overSqr, tmpMat)
        end do

      end subroutine getSRmatrices

  end subroutine getLr2ndTerms_


  !> Calculate R*T contribution of gradient from pc terms
  subroutine getPc2ndTerms_(env, coord0, overSqr, RmatSpL, RdelSpL, &
      & weight, extCharges, blurWidths, rVec, gVec, alpha, vol, getDenseAO, &
      & getAtomIndex, orderRmatL, SAstates, tNAC, tPeriodic, tBlur, deriv1, deriv2)

    !> Environment settings
    type(TEnvironment), intent(inout) :: env

    !> central cell coordinates of atoms
    real(dp), intent(inout) :: coord0(:,:)

    !> Dense overlap matrix
    real(dp), intent(in) :: overSqr(:,:)

    !> auxiliary matrix in AO basis related to SA-REKS term with sparse form
    real(dp), intent(in) :: RmatSpL(:,:,:)

    !> auxiliary matrix in AO basis related to state-interaction term with sparse form
    real(dp), allocatable, intent(in) :: RdelSpL(:,:,:)

    !> Weight of each microstate for state to be optimized; weight = weightL * SAweight
    real(dp), intent(in) :: weight(:)

    !> coordinates and charges of external point charges
    real(dp), intent(in) :: extCharges(:,:)

    !> Width of the Gaussians if the charges are blurred
    real(dp), allocatable, intent(in) :: blurWidths(:)

    !> real lattice points for Ewald-sum
    real(dp), allocatable, intent(in) :: rVec(:,:)

    !> lattice points for reciprocal Ewald
    real(dp), allocatable, intent(in) :: gVec(:,:)

    !> parameter for Ewald
    real(dp), intent(in) :: alpha

    !> parameter for cell volume
    real(dp), intent(in) :: vol

    !> get dense AO index from sparse AO array
    integer, intent(in) :: getDenseAO(:,:)

    !> get atom index from AO index
    integer, intent(in) :: getAtomIndex(:)

    !> Ordering between RmatL and fillingL
    integer, intent(in) :: orderRmatL(:)

    !> Number of states used in state-averaging
    integer, intent(in) :: SAstates

    !> Calculate nonadiabatic coupling vectors
    logical, intent(in) :: tNAC

    !> if calculation is periodic
    logical, intent(in) :: tPeriodic

    !> If charges should be blured
    logical, intent(in) :: tBlur

    !> computed tr(R*T) gradient for SA-REKS, SSR, or L state
    real(dp), intent(inout) :: deriv1(:,:,:)

    !> computed tr(R*T) gradient for state-interaction term
    real(dp), intent(inout) :: deriv2(:,:,:)

    real(dp), allocatable :: QinvRderiv(:,:)      ! ... nAtom, 3
    real(dp), allocatable :: Tderiv(:,:,:)        ! ... sparseSize, 3, Ncpu

    real(dp) :: tmpCoulombDeriv(3)
    integer :: iAtom1, iAtom3, iAtom4, nAtom
    integer :: ist, nstates, nstHalf, mu, nu, l, sparseSize, id, Ncpu

    nAtom = size(coord0,dim=2)
    nstates = size(RmatSpL,dim=3)
    nstHalf = nstates * (nstates - 1) / 2

    sparseSize = size(getDenseAO,dim=1)

  #:if WITH_OMP
!$OMP PARALLEL
    Ncpu = OMP_GET_NUM_THREADS()
!$OMP END PARALLEL
  #:else
    Ncpu = 1
  #:endif

    allocate(QinvRderiv(3,nAtom))
    allocate(Tderiv(sparseSize,3,Ncpu))

    ! contributions related to gradient of external charges
    ! Q_{pc} * (-1/R**2) between QM and PC
    call getQinvRderiv(env, coord0, extCharges(1:3,:), extCharges(4,:), &
        & blurWidths, rVec, gVec, alpha, vol, tPeriodic, tBlur, QinvRderiv)

    ! compute R*T shift with only up-spin part of Tderiv due to symmetry
    ! contribution for derivative of 1/r coulomb potential w.r.t. R_atom, not R_atompair
!$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(id,mu,nu,iAtom3,iAtom4, &
!$OMP& tmpCoulombDeriv) REDUCTION(+:deriv1,deriv2) SCHEDULE(RUNTIME)
    do iAtom1 = 1, nAtom
      
    #:if WITH_OMP
      id = OMP_GET_THREAD_NUM() + 1
    #:else
      id = 1
    #:endif

      ! zeroing for temporary Tderiv in each k atom pair
      ! calculate sparse Tderiv in AO basis
      Tderiv(:,:,id) = 0.0_dp
      do l = 1, sparseSize

        ! set the AO indices with respect to sparsity
        mu = getDenseAO(l,1)
        nu = getDenseAO(l,2)
        ! find proper atom index
        iAtom3 = getAtomIndex(mu)
        iAtom4 = getAtomIndex(nu)

        if (mu <= nu) then
          if (abs(overSqr(mu,nu)) >= epsilon(1.0_dp)) then

            tmpCoulombDeriv(:) = 0.0_dp
            ! mu in alpha
            if (iAtom3 == iAtom1) then
              tmpCoulombDeriv(:) = tmpCoulombDeriv + QinvRderiv(:,iAtom1)
            end if
            ! nu in alpha
            if (iAtom4 == iAtom1) then
              tmpCoulombDeriv(:) = tmpCoulombDeriv + QinvRderiv(:,iAtom1)
            end if

            ! calculate the contribution of external charges to T derivative
            ! originally, the hamiltonian has minus sign, but
            ! we already multiply -1 in externalcharges.F90 routine.
            ! thus, we have to use + sign in the process of gradient
            Tderiv(l,:,id) = Tderiv(l,:,id) &
                & + 0.5_dp*overSqr(mu,nu)*tmpCoulombDeriv(:)

          end if
        end if

      end do
      ! end of loop l

      if (tNAC) then
        do ist = 1, nstates
          if (ist /= SAstates) then
            call shiftRTgradPc_(deriv1(:,:,ist), RmatSpL(:,:,ist), &
                & Tderiv(:,:,id), weight, orderRmatL, iAtom1)
          end if
        end do
        do ist = 1, nstHalf
          call shiftRTgradPc_(deriv2(:,:,ist), RdelSpL(:,:,ist), &
              & Tderiv(:,:,id), weight, orderRmatL, iAtom1)
        end do
      else
        call shiftRTgradPc_(deriv1(:,:,1), RmatSpL(:,:,1), &
            & Tderiv(:,:,id), weight, orderRmatL, iAtom1)
      end if

    end do
!$OMP END PARALLEL DO

    contains

      !> contributions related to gradient of external charges
      !> Q_{pc} * (-1/R**2) between QM and PC
      subroutine getQinvRderiv(env, qmCoords, pcCoords, pcCharges, &
          & blurWidths, rVec, gVec, alpha, vol, tPeriodic, tBlur, QinvRderiv)

        !> Environment settings
        type(TEnvironment), intent(inout) :: env

        !> central cell coordinates of atoms
        real(dp), intent(in) :: qmCoords(:,:)

        !> coordinates of external point charges
        real(dp), intent(in) :: pcCoords(:,:)

        !> charges of external point charges
        real(dp), intent(in) :: pcCharges(:)

        !> Width of the Gaussians if the charges are blurred
        real(dp), allocatable, intent(in) :: blurWidths(:)

        !> real lattice points for Ewald-sum
        real(dp), allocatable, intent(in) :: rVec(:,:)

        !> lattice points for reciprocal Ewald
        real(dp), allocatable, intent(in) :: gVec(:,:)

        !> parameter for Ewald
        real(dp), intent(in) :: alpha

        !> parameter for cell volume
        real(dp), intent(in) :: vol

        !> if calculation is periodic
        logical, intent(in) :: tPeriodic

        !> If charges should be blured
        logical, intent(in) :: tBlur

        !> Q_{pc} * (-1/R**2) between QM and PC
        real(dp), intent(out) :: QinvRderiv(:,:)

        real(dp), allocatable :: tmpCharges(:)
        real(dp), allocatable :: tmpDeriv(:,:)
        integer :: nAtom, nAtomPc

        nAtom = size(qmCoords,dim=2)
        nAtomPc = size(pcCoords,dim=2)

        allocate(tmpCharges(nAtom))
        allocate(tmpDeriv(3, nAtomPc))

        QinvRderiv(:,:) = 0.0_dp
        if (tPeriodic) then
          if (tBlur) then
            call addInvRPrime(env, nAtom, nAtomPc, qmCoords, pcCoords, tmpCharges, pcCharges, &
                & rVec, gVec, alpha, vol, QinvRderiv, tmpDeriv, tHamDeriv=.true., &
                & blurWidths1=blurWidths)
          else
            call addInvRPrime(env, nAtom, nAtomPc, qmCoords, pcCoords, tmpCharges, pcCharges, &
                & rVec, gVec, alpha, vol, QinvRderiv, tmpDeriv, tHamDeriv=.true.)
          end if
        else
          if (tBlur) then
            call addInvRPrime(env, nAtom, nAtomPc, qmCoords, pcCoords, tmpCharges, &
                & pcCharges, QinvRderiv, tmpDeriv, tHamDeriv=.true., blurWidths1=blurWidths)
          else
            call addInvRPrime(env, nAtom, nAtomPc, qmCoords, pcCoords, tmpCharges, &
                & pcCharges, QinvRderiv, tmpDeriv, tHamDeriv=.true.)
          end if
        end if

      end subroutine getQinvRderiv

  end subroutine getPc2ndTerms_


  !> Calculate Q*S gradient contribution
  subroutine shiftQSgrad_(Qmat, Sderiv, fac, iSquare, mOrb, deriv)

    !> auxiliary matrix in AO basis related to SA-REKS term
    real(dp), intent(in) :: Qmat(:,:)

    !> Dense overlap derivative in AO basis
    real(dp), intent(in) :: Sderiv(:,:,:)

    !> factor of Q*S term, for SI term it becomes -1 otherwise 1
    real(dp), intent(in) :: fac

    !> Position of each atom in the rows/columns of the square matrices. Shape: (nAtom)
    integer, intent(in) :: iSquare(:)

    !> Max. nr. of orbitals for any species
    integer, intent(in) :: mOrb

    !> gradient contribution from tr(Q*S)
    real(dp), intent(inout) :: deriv(:,:)

    real(dp), allocatable :: sqrQtmp(:,:)
    real(dp), allocatable :: sqrStmp(:,:,:)

    real(dp) :: derivTmp(3)
    integer :: iAtom1, nOrb1, iAtom2, nOrb2, ii
    integer :: nAtom, nAO1, nAO2, iOrb1, iOrb2

    nAtom = size(deriv,dim=2)

    allocate(sqrQtmp(mOrb,mOrb))
    allocate(sqrStmp(mOrb,mOrb,3))

    do iAtom1 = 1, nAtom
      nOrb1 = iSquare(iAtom1+1) - iSquare(iAtom1)
      nAO1 = iSquare(iAtom1) - 1
      do iAtom2 = iAtom1, nAtom
        nOrb2 = iSquare(iAtom2+1) - iSquare(iAtom2)
        nAO2 = iSquare(iAtom2) - 1
        if (iAtom1 /= iAtom2) then

          sqrQtmp(:,:) = 0.0_dp
          sqrStmp(:,:,:) = 0.0_dp
          do iOrb1 = 1, nOrb1
            do iOrb2 = 1, nOrb2
              sqrQtmp(iOrb2,iOrb1) = Qmat(nAO2+iOrb2,nAO1+iOrb1)
              do ii = 1, 3
                sqrStmp(iOrb2,iOrb1,ii) = Sderiv(nAO2+iOrb2,nAO1+iOrb1,ii)
              end do
            end do
          end do

          derivTmp(:) = 0.0_dp
          do ii = 1, 3
            derivTmp(ii) = sum(sqrQtmp(1:nOrb2,1:nOrb1)*sqrStmp(1:nOrb2,1:nOrb1,ii))
          end do

          ! forces from atom 1 on atom 2f and 2f onto 1
          deriv(:,iAtom1) = deriv(:,iAtom1) + fac*derivTmp(:)
          deriv(:,iAtom2) = deriv(:,iAtom2) - fac*derivTmp(:)

        end if
      end do
    end do

  end subroutine shiftQSgrad_


  !> Calculate R*T gradient contribution with sparse or half dense form
  subroutine shiftRTgradSparse_(deriv, RmatSpL, TderivL, &
      & weight, orderRmatL, iAtom1, iAtom2)

    !> auxiliary matrix in AO basis related to SA-REKS term with sparse form
    real(dp), intent(in) :: RmatSpL(:,:)

    !> T derivative for each microstate with sparse form
    real(dp), intent(in) :: TderivL(:,:,:)

    !> Weight of each microstate for state to be optimized; weight = weightL * SAweight
    real(dp), intent(in) :: weight(:)

    !> Ordering between RmatL and fillingL
    integer, intent(in) :: orderRmatL(:)

    !> index related to gradient
    integer, intent(in) :: iAtom1, iAtom2

    !> temporary derivatives
    real(dp), intent(inout) :: deriv(:,:)

    real(dp) :: derivTmp(3)
    integer :: ii, iL, Lmax, tmpL

    Lmax = size(weight,dim=1)

    do iL = 1, Lmax
      tmpL = orderRmatL(iL)
      derivTmp(:) = 0.0_dp
      do ii = 1, 3
        derivTmp(ii) = sum(RmatSpL(:,tmpL)*TderivL(:,ii,iL))
      end do
      ! forces from atom 1 on atom 2f and 2f onto 1
      deriv(:,iAtom1) = deriv(:,iAtom1) + 2.0_dp*weight(iL)*derivTmp(:)
      deriv(:,iAtom2) = deriv(:,iAtom2) - 2.0_dp*weight(iL)*derivTmp(:)
    end do

  end subroutine shiftRTgradSparse_


  !> Calculate R*T gradient contribution with LC term
  subroutine shiftRTgradLr1st_(deriv, RmatL, PmatL, LrGammaDeriv, &
      & weight, iSquare, orderRmatL, mOrb, iAtom1, iAtom2)

    !> auxiliary matrix in AO basis related to SA-REKS term
    real(dp), intent(in) :: RmatL(:,:,:)

    !> product of overlap and density matrix
    real(dp), intent(in) :: PmatL(:,:,:)

    !> long-range gamma derivative integrals
    real(dp), intent(in) :: LrGammaDeriv(:,:,:)

    !> Weight of each microstate for state to be optimized; weight = weightL * SAweight
    real(dp), intent(in) :: weight(:)

    !> Position of each atom in the rows/columns of the square matrices. Shape: (nAtom)
    integer, intent(in) :: iSquare(:)

    !> Ordering between RmatL and fillingL
    integer, intent(in) :: orderRmatL(:)

    !> Max. nr. of orbitals for any species
    integer, intent(in) :: mOrb

    !> index related to gradient
    integer, intent(in) :: iAtom1, iAtom2

    !> temporary derivatives
    real(dp), intent(inout) :: deriv(:,:)

    real(dp), allocatable :: sqrPtmp(:,:,:)
    real(dp), allocatable :: sqrRtmp(:,:,:)

    real(dp) :: sqrLRtmp(3), derivTmp(3)
    real(dp) :: fac, tmpval
    integer :: nAO1, nOrb1, nAO2, nOrb2
    integer :: mu, nu, tmpL, iL, LmaxR, Lmax, ii

    Lmax = size(PmatL,dim=3)
    LmaxR = size(RmatL,dim=3)
    fac = -0.25_dp

    allocate(sqrPtmp(mOrb,mOrb,Lmax))
    allocate(sqrRtmp(mOrb,mOrb,LmaxR))

    nOrb1 = iSquare(iAtom1+1) - iSquare(iAtom1)
    nAO1 = iSquare(iAtom1) - 1
    nOrb2 = iSquare(iAtom2+1) - iSquare(iAtom2)
    nAO2 = iSquare(iAtom2) - 1

    sqrPtmp(:,:,:) = 0.0_dp
    sqrRtmp(:,:,:) = 0.0_dp
    sqrLRtmp(:) = 0.0_dp
    do mu = 1, nOrb1
      do nu = 1, nOrb2
        do iL = 1, Lmax
          sqrPtmp(nu,mu,iL) = PmatL(nAO2+nu,nAO1+mu,iL)
          if (iL <= LmaxR) then
            sqrRtmp(nu,mu,iL) = RmatL(nAO2+nu,nAO1+mu,iL) &
                & + RmatL(nAO1+mu,nAO2+nu,iL)
          end if
        end do
      end do
    end do
    do ii = 1, 3
      sqrLRtmp(ii) = LrGammaDeriv(iAtom2,iAtom1,ii)
    end do

    do iL = 1, Lmax
      derivTmp(:) = 0.0_dp
      tmpL = orderRmatL(iL)
      tmpval = sum(sqrRtmp(:,:,tmpL)*sqrPtmp(:,:,iL))
      do ii = 1, 3
        derivTmp(ii) = fac * tmpval * sqrLRtmp(ii)
      end do
      ! forces from atom 1 on atom 2f and 2f onto 1
      deriv(:,iAtom1) = deriv(:,iAtom1) + 2.0_dp*weight(iL)*derivTmp(:)
      deriv(:,iAtom2) = deriv(:,iAtom2) - 2.0_dp*weight(iL)*derivTmp(:)
    end do

  end subroutine shiftRTgradLr1st_


  !> Calculate R*T gradient contribution with LC term
  subroutine shiftRTgradLr2nd_(deriv, RmatL, PmatL, LrGammaDeriv, &
      & weight, iSquare, orderRmatL, mOrb, iAtom1, iAtom2, option)

    !> auxiliary matrix in AO basis related to SA-REKS term
    real(dp), intent(in) :: RmatL(:,:,:)

    !> product of overlap and density matrix
    real(dp), intent(in) :: PmatL(:,:,:)

    !> long-range gamma derivative integrals
    real(dp), intent(in) :: LrGammaDeriv(:,:,:)

    !> Weight of each microstate for state to be optimized; weight = weightL * SAweight
    real(dp), intent(in) :: weight(:)

    !> Position of each atom in the rows/columns of the square matrices. Shape: (nAtom)
    integer, intent(in) :: iSquare(:)

    !> Ordering between RmatL and fillingL
    integer, intent(in) :: orderRmatL(:)

    !> Max. nr. of orbitals for any species
    integer, intent(in) :: mOrb

    !> index related to gradient
    integer, intent(in) :: iAtom1, iAtom2

    !> option for S*P product
    integer, intent(in) :: option

    !> temporary derivatives
    real(dp), intent(inout) :: deriv(:,:)

    real(dp), allocatable :: sqrPtmp1(:,:,:)
    real(dp), allocatable :: sqrPtmp2(:,:,:)
    real(dp), allocatable :: sqrRtmp1(:,:,:)
    real(dp), allocatable :: sqrRtmp2(:,:,:)

    real(dp) :: sqrLRtmp(3), derivTmp(3)
    real(dp) :: fac, tmpval1, tmpval2
    integer :: nAO1, nOrb1, nAO2, nOrb2
    integer :: mu, nu, tmpL, iL, LmaxR, Lmax, ii

    Lmax = size(PmatL,dim=3)
    LmaxR = size(RmatL,dim=3)
    fac = -0.25_dp

    allocate(sqrPtmp1(mOrb,mOrb,Lmax))
    allocate(sqrRtmp1(mOrb,mOrb,LmaxR))
    allocate(sqrPtmp2(mOrb,mOrb,Lmax))
    allocate(sqrRtmp2(mOrb,mOrb,LmaxR))

    nOrb1 = iSquare(iAtom1+1) - iSquare(iAtom1)
    nAO1 = iSquare(iAtom1) - 1
    nOrb2 = iSquare(iAtom2+1) - iSquare(iAtom2)
    nAO2 = iSquare(iAtom2) - 1

    sqrPtmp1(:,:,:) = 0.0_dp
    sqrPtmp2(:,:,:) = 0.0_dp
    sqrRtmp1(:,:,:) = 0.0_dp
    sqrRtmp2(:,:,:) = 0.0_dp
    sqrLRtmp(:) = 0.0_dp
    do mu = 1, nOrb1
      do nu = 1, nOrb2
        do iL = 1, Lmax
          if (option == 1) then
            sqrPtmp1(nu,mu,iL) = PmatL(nAO2+nu,nAO1+mu,iL)
            sqrPtmp2(nu,mu,iL) = PmatL(nAO1+mu,nAO2+nu,iL)
          else if (option == 2) then
            sqrPtmp2(nu,mu,iL) = PmatL(nAO2+nu,nAO1+mu,iL)
            sqrPtmp1(nu,mu,iL) = PmatL(nAO1+mu,nAO2+nu,iL)
          end if
          if (iL <= LmaxR) then
            sqrRtmp1(nu,mu,iL) = RmatL(nAO2+nu,nAO1+mu,iL)
            sqrRtmp2(nu,mu,iL) = RmatL(nAO1+mu,nAO2+nu,iL)
          end if
        end do
      end do
    end do
    do ii = 1, 3
      sqrLRtmp(ii) = LrGammaDeriv(iAtom2,iAtom1,ii)
    end do

    do iL = 1, Lmax
      derivTmp(:) = 0.0_dp
      tmpL = orderRmatL(iL)
      tmpval1 = sum(sqrRtmp1(:,:,tmpL)*sqrPtmp1(:,:,iL))
      tmpval2 = sum(sqrRtmp2(:,:,tmpL)*sqrPtmp2(:,:,iL))
      do ii = 1, 3
        derivTmp(ii) = fac * (tmpval1 + tmpval2) * sqrLRtmp(ii)
      end do
      ! forces from atom 1 on atom 2f and 2f onto 1
      deriv(:,iAtom1) = deriv(:,iAtom1) + 2.0_dp*weight(iL)*derivTmp(:)
      deriv(:,iAtom2) = deriv(:,iAtom2) - 2.0_dp*weight(iL)*derivTmp(:)
    end do

  end subroutine shiftRTgradLr2nd_


  !> Calculate R*T gradient contribution with sparse form in point charges
  subroutine shiftRTgradPc_(deriv, RmatSpL, Tderiv, weight, &
      & orderRmatL, iAtom1)

    !> auxiliary matrix in AO basis related to SA-REKS term with sparse form
    real(dp), intent(in) :: RmatSpL(:,:)

    !> T derivative with sparse form
    real(dp), intent(in) :: Tderiv(:,:)

    !> Weight of each microstate for state to be optimized; weight = weightL * SAweight
    real(dp), intent(in) :: weight(:)

    !> Ordering between RmatL and fillingL
    integer, intent(in) :: orderRmatL(:)

    !> index related to gradient
    integer, intent(in) :: iAtom1

    !> temporary derivatives
    real(dp), intent(inout) :: deriv(:,:)

    real(dp) :: derivTmp(3)
    integer :: ii, iL, Lmax, tmpL

    Lmax = size(weight,dim=1)

    do iL = 1, Lmax
      tmpL = orderRmatL(iL)
      derivTmp(:) = 0.0_dp
      do ii = 1, 3
        derivTmp(ii) = sum(RmatSpL(:,tmpL)*Tderiv(:,ii))
      end do
      ! forces from atom 1 on atom 2f and 2f onto 1
      deriv(:,iAtom1) = deriv(:,iAtom1) + 2.0_dp*weight(iL)*derivTmp(:)
    end do

  end subroutine shiftRTgradPc_


  !> compute the gradient for remaining SA-REKS(2,2) state
  subroutine getOtherSAgrad22_(avgGrad, SAgrad)

    !> gradient of averaged state
    real(dp), intent(in) :: avgGrad(:,:)

    !> gradient of SA-REKS state
    real(dp), intent(inout) :: SAgrad(:,:,:)

    SAgrad(:,:,2) = 2.0_dp*avgGrad(:,:) - SAgrad(:,:,1)

  end subroutine getOtherSAgrad22_


end module dftbp_reksgrad
