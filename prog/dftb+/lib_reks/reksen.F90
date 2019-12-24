!--------------------------------------------------------------------------------------------------!
!  DFTB+: general package for performing fast atomistic simulations                                !
!  Copyright (C) 2006 - 2019  DFTB+ developers group                                               !
!                                                                                                  !
!  See the LICENSE file for terms of usage and distribution.                                       !
!--------------------------------------------------------------------------------------------------!

!> REKS and SI-SA-REKS formulation in DFTB as developed by Lee et al.
!>
!> The functionality of the module has some limitation:
!> * Third order does not work.
!> * Periodic system do not work yet appart from Gamma point.
!> * Orbital potentials or spin-orbit or external E-field does not work yet.
!> * Only for closed shell system.
!> * Onsite corrections are not included in this version
module dftbp_reksen

  use dftbp_accuracy
  use dftbp_blasroutines, only : gemm
  use dftbp_densedescr
  use dftbp_eigenvects
  use dftbp_elecsolvers
  use dftbp_energies
  use dftbp_environment
  use dftbp_globalenv
  use dftbp_mainio
  use dftbp_message
  use dftbp_periodic
  use dftbp_sparse2dense
  use dftbp_rekscommon

  implicit none

  private

  public :: constructMicrostates, calcWeights
  public :: activeOrbSwap, getFilling, calcSaReksEnergy
  public :: getFockandDiag, guessNewEigvecs
  public :: adjustEigenval, solveSecularEqn

  contains

  !> Construct L, spin dependent microstates from identical KS orbitals
  subroutine constructMicrostates(Nc, tSSR22, tSSR44, fillingL)

    !> Number of core orbitals
    integer, intent(in) :: Nc

    !> Calculate DFTB/SSR(2,2) formalism
    logical, intent(in) :: tSSR22

    !> Calculate DFTB/SSR(4,4) formalism
    logical, intent(in) :: tSSR44

    !> Filling for each microstate
    real(dp), intent(out) :: fillingL(:,:,:)

    if (tSSR22) then
      call getFillingL22_(Nc, fillingL)
    else if (tSSR44) then
      call error("SSR(4,4) is not implemented yet")
    end if

  end subroutine constructMicrostates


  !> Calculate the weight of each microstate for current cycle, C_L
  subroutine calcWeights(FONs, delta, SAweight, tSSR22, tSSR44, weightL, weight)

    !> Fractional occupation numbers of active orbitals
    real(dp), intent(in) :: FONs(:,:)

    !> Smoothing factor used in FON optimization
    real(dp), intent(in) :: delta

    !> Weights used in state-averaging
    real(dp), intent(in) :: SAweight(:)

    !> Calculate DFTB/SSR(2,2) formalism
    logical, intent(in) :: tSSR22

    !> Calculate DFTB/SSR(4,4) formalism
    logical, intent(in) :: tSSR44

    !> Weight for each microstate per state
    real(dp), intent(out) :: weightL(:,:)

    !> Weight of each microstate for state to be optimized; weight = weightL * SAweight
    real(dp), intent(out) :: weight(:)

    if (tSSR22) then
      call getWeightL22_(FONs, delta, SAweight, weightL, weight)
    else if (tSSR44) then
      call error("SSR(4,4) is not implemented yet")
    end if

  end subroutine calcWeights


  !> Swap the active orbitals for feasible occupation in REKS
  subroutine activeOrbSwap(eigenvecs, SAweight, FONs, Efunction, &
      & Nc, tSSR22, tSSR44)

    !> eigenvectors
    real(dp), intent(inout) :: eigenvecs(:,:)

    !> Weights used in state-averaging
    real(dp), intent(in) :: SAweight(:)

    !> Fractional occupation numbers of active orbitals
    real(dp), intent(in) :: FONs(:,:)

    !> Minimized energy functional
    integer, intent(in) :: Efunction

    !> Number of core orbitals
    integer, intent(in) :: Nc

    !> Calculate DFTB/SSR(2,2) formalism
    logical, intent(in) :: tSSR22

    !> Calculate DFTB/SSR(4,4) formalism
    logical, intent(in) :: tSSR44

    if (tSSR22) then
      call MOswap22_(eigenvecs, SAweight, FONs, Efunction, Nc)
    else if (tSSR44) then
      call error("SSR(4,4) is not implemented yet")
    end if

  end subroutine activeOrbSwap


  !> Calculate filling for minimzed state with optimized FONs
  subroutine getFilling(filling, SAweight, FONs, Efunction, &
      & Nc, tSSR22, tSSR44)

    !> occupations (level)
    real(dp), intent(out) :: filling(:)

    !> Weights used in state-averaging
    real(dp), intent(in) :: SAweight(:)

    !> Fractional occupation numbers of active orbitals
    real(dp), intent(in) :: FONs(:,:)

    !> Minimized energy functional
    integer, intent(in) :: Efunction

    !> Number of core orbitals
    integer, intent(in) :: Nc

    !> Calculate DFTB/SSR(2,2) formalism
    logical, intent(in) :: tSSR22

    !> Calculate DFTB/SSR(4,4) formalism
    logical, intent(in) :: tSSR44

    if (tSSR22) then
      call getFilling22_(filling, SAweight, FONs, Efunction, Nc)
    else if (tSSR44) then
      call error("SSR(4,4) is not implemented yet")
    end if

  end subroutine getFilling


  !> Calculate the energy of SA-REKS states and averaged state
  subroutine calcSaReksEnergy(SAweight, weightL, enLnonSCC, enLscc, &
      & enLspin, enL3rd, enLfock, enLtot, rstate, t3rd, tRangeSep, &
      & energy, EnonSCC, Escc, Espin, e3rd, Efock, Eelec, Etotal)

    !> Weights used in state-averaging
    real(dp), intent(in) :: SAweight(:)

    !> Weight for each microstate per state
    real(dp), intent(in) :: weightL(:,:)

    !> non SCC energy for each microstate
    real(dp), intent(in) :: enLnonSCC(:)

    !> SCC energy for each microstate
    real(dp), intent(in) :: enLscc(:)

    !> spin-polarized energy for each microstate
    real(dp), intent(in) :: enLspin(:)

    !> 3rd order SCC energy for each microstate
    real(dp), intent(in) :: enL3rd(:)

    !> Long-range corrected energy for each microstate
    real(dp), intent(in) :: enLfock(:)

    !> total energy for each microstate
    real(dp), intent(in) :: enLtot(:)

    !> Target SSR state
    integer, intent(in) :: rstate

    !> Third order DFTB
    logical, intent(in) :: t3rd

    !> Whether to run a range separated calculation
    logical, intent(in) :: tRangeSep

    !> energy of states
    real(dp), intent(out) :: energy(:)

    !> Non-SCC energy
    real(dp), intent(inout) :: EnonSCC

    !> SCC energy
    real(dp), intent(inout) :: Escc

    !> spin energy
    real(dp), intent(inout) :: Espin

    !> Total 3rd order
    real(dp), intent(inout) :: e3rd

    !> range-separation energy
    real(dp), intent(inout) :: Efock

    !> total electronic energy
    real(dp), intent(inout) :: Eelec

    !> energy of averaged state
    real(dp), intent(out) :: Etotal

    integer :: ist, iL, SAstates, nstates, Lmax

    Lmax = size(weightL,dim=2)
    SAstates = size(SAweight,dim=1)
    nstates = size(weightL,dim=1)

    ! Compute the energy contributions of target SA-REKS state
    ! energy = nonSCC + scc + spin + 3rd + fock
    EnonSCC = 0.0_dp
    Escc = 0.0_dp
    Espin = 0.0_dp
    if (t3rd) then
      e3rd = 0.0_dp
    end if
    if (tRangeSep) then
      Efock = 0.0_dp
    end if
    do iL = 1, Lmax
      EnonSCC = EnonSCC + weightL(rstate,iL) * enLnonSCC(iL)
      Escc = Escc + weightL(rstate,iL) * enLSCC(iL)
      Espin = Espin + weightL(rstate,iL) * enLspin(iL)
      if (t3rd) then
        e3rd = e3rd + weightL(rstate,iL) * enL3rd(iL)
      end if
      if (tRangeSep) then
        Efock = Efock + weightL(rstate,iL) * enLfock(iL)
      end if
    end do
    Eelec = EnonSCC + Escc + Espin
    if (t3rd) then
      Eelec = Eelec + e3rd
    end if
    if (tRangeSep) then
      Eelec = Eelec + Efock
    end if

    ! Compute the energy of SA-REKS states
    energy(:) = 0.0_dp
    do ist = 1, nstates
      do iL = 1, Lmax
        energy(ist) = energy(ist) + weightL(ist,iL) * enLtot(iL)
      end do
    end do

    ! In this step Etotal is energy of averaged state, not individual states
    ! From this energy we can check the variational principle
    Etotal = 0.0_dp
    do ist = 1, SAstates
      Etotal = Etotal + SAweight(ist) * energy(ist)
    end do

  end subroutine calcSaReksEnergy


  !> Make pseudo-fock operator with Hamiltonian of each microstate
  !> and diagonalize the fock matrix
  subroutine getFockandDiag(env, denseDesc, neighbourList, &
      & nNeighbourSK, iSparseStart, img2CentCell, eigenvecs, &
      & electronicSolver, eigen, hamSqrL, hamSpL, weight, &
      & fillingL, shift, Nc, Na, Lpaired, tRangeSep, fockFc, &
      & fockFa, fock, eigvecsFock)

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

    !> Electronic solver information
    type(TElectronicSolver), intent(inout) :: electronicSolver

    !> eigenvalues
    real(dp), intent(out) :: eigen(:,:,:)

    !> Dense Hamiltonian matrix for each microstate
    real(dp), intent(inout) :: hamSqrL(:,:,:,:)

    !> Sparse Hamiltonian matrix for each microstate
    real(dp), intent(inout) :: hamSpL(:,:,:)

    !> Weight of each microstate for state to be optimized; weight = weightL * SAweight
    real(dp), intent(in) :: weight(:)

    !> Filling for each microstate
    real(dp), intent(in) :: fillingL(:,:,:)

    !> Shift value in SCC cycle
    real(dp), intent(in) :: shift

    !> Number of core orbitals
    integer, intent(in) :: Nc

    !> Number of active orbitals
    integer, intent(in) :: Na

    !> Number of spin-paired microstates
    integer, intent(in) :: Lpaired

    !> Whether to run a range separated calculation
    logical, intent(in) :: tRangeSep

    !> dense fock matrix for core orbitals
    real(dp), intent(inout) :: fockFc(:,:)

    !> dense fock matrix for active orbitals
    real(dp), intent(inout) :: fockFa(:,:,:)

    !> dense pseudo-fock matrix
    real(dp), intent(inout) :: fock(:,:)

    !> eigenvectors from pesudo-fock matrix
    real(dp), intent(inout) :: eigvecsFock(:,:)

    real(dp), allocatable :: orbFON(:)
    real(dp), allocatable :: tmpOver(:,:)
    real(dp), allocatable :: tmpMat(:,:)

    integer :: ii, nOrb

    nOrb = size(fockFc,dim=1)

    allocate(orbFON(nOrb))
    allocate(tmpOver(nOrb,nOrb))
    allocate(tmpMat(nOrb,nOrb))

    call getFockFcFa_(env, denseDesc, neighbourList, nNeighbourSK, &
        & iSparseStart, img2CentCell, hamSqrL, hamSpL, weight, &
        & fillingL, Nc, Na, Lpaired, tRangeSep, orbFON, fockFc, fockFa)

    call matAO2MO(fockFc, eigenvecs(:,:,1))
    do ii = 1, Na
      call matAO2MO(fockFa(:,:,ii), eigenvecs(:,:,1))
    end do

    call getPseudoFock_(fockFc, fockFa, orbFON, Nc, Na, fock)

    call levelShifting_(fock, shift, Nc, Na)

    ! Diagonalize the pesudo-Fock matrix
    tmpOver(:,:) = 0.0_dp
    do ii = 1, nOrb
      tmpOver(ii,ii) = 1.0_dp
    end do
    tmpMat(:,:) = fock

    eigen(:,1,1) = 0.0_dp
    call diagDenseMtx(electronicSolver, 'V', tmpMat, tmpOver, eigen(:,1,1))
    eigvecsFock(:,:) = tmpMat

  end subroutine getFockandDiag


  !> guess new eigenvectors from Fock eigenvectors
  subroutine guessNewEigvecs(eigenvecs, eigvecsFock)

    !> Eigenvectors on eixt
    real(dp), intent(inout) :: eigenvecs(:,:)

    !> eigenvectors from pesudo-fock matrix
    real(dp), intent(in) :: eigvecsFock(:,:)

    real(dp), allocatable :: tmpVec(:,:)
    integer :: nOrb

    nOrb = size(eigvecsFock,dim=1)

    allocate(tmpVec(nOrb,nOrb))

    tmpVec(:,:) = 0.0_dp
    call gemm(tmpVec, eigenvecs, eigvecsFock)
    eigenvecs(:,:) = tmpVec

  end subroutine guessNewEigvecs


  !> adjust the eigenvalues (eliminate shift values)
  subroutine adjustEigenval(eigen, shift, Nc, Na)

    !> eigenvalues
    real(dp), intent(inout) :: eigen(:,:,:)

    !> Shift value in SCC cycle
    real(dp), intent(in) :: shift

    !> Number of core orbitals
    integer, intent(in) :: Nc

    !> Number of active orbitals
    integer, intent(in) :: Na

    integer :: nOrb, ind, ii

    nOrb = size(eigen,dim=1)

    do ii = Nc + 1, Nc + Na
      ind = ii - Nc
      eigen(ii,1,1) = eigen(ii,1,1) - dble(ind) * shift
    end do

    do ii = Nc + Na + 1, nOrb
      ind = Na + 1
      eigen(ii,1,1) = eigen(ii,1,1) - dble(ind) * shift
    end do

  end subroutine adjustEigenval


  !> Solve secular equation with coupling element between SA-REKS states
  subroutine solveSecularEqn(env, denseDesc, neighbourList, &
      & nNeighbourSK, iSparseStart, img2CentCell, electronicSolver, &
      & eigenvecs, hamSqrL, hamSpL, weight, FONs, fillingL, &
      & Elevel, useSSR, Lpaired, Nc, Na, tRangeSep, tSSR22, &
      & tSSR44, energy, eigvecsSSR)

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

    !> Electronic solver information
    type(TElectronicSolver), intent(inout) :: electronicSolver

    !> Eigenvectors on eixt
    real(dp), intent(in) :: eigenvecs(:,:,:)

    !> Dense Hamiltonian matrix for each microstate
    real(dp), intent(inout) :: hamSqrL(:,:,:,:)

    !> Sparse Hamiltonian matrix for each microstate
    real(dp), intent(inout) :: hamSpL(:,:,:)

    !> Weight of each microstate for state to be optimized; weight = weightL * SAweight
    real(dp), intent(in) :: weight(:)

    !> Fractional occupation numbers of active orbitals
    real(dp), intent(in) :: FONs(:,:)

    !> Filling for each microstate
    real(dp), intent(in) :: fillingL(:,:,:)

    !> Calculated energy states in SA-REKS
    integer, intent(in) :: Elevel

    !> Calculate SSR state (SI term is included)
    integer, intent(in) :: useSSR

    !> Number of spin-paired microstates
    integer, intent(in) :: Lpaired

    !> Number of core orbitals
    integer, intent(in) :: Nc

    !> Number of active orbitals
    integer, intent(in) :: Na

    !> Whether to run a range separated calculation
    logical, intent(in) :: tRangeSep

    !> Calculate DFTB/SSR(2,2) formalism
    logical, intent(in) :: tSSR22

    !> Calculate DFTB/SSR(4,4) formalism
    logical, intent(in) :: tSSR44

    !> energy of states
    real(dp), intent(inout) :: energy(:)

    !> eigenvectors from SA-REKS state
    real(dp), intent(inout) :: eigvecsSSR(:,:)

    real(dp), allocatable :: Wab(:,:)
    real(dp), allocatable :: StateCoup(:,:)
    real(dp), allocatable :: tmpOver(:,:)
    real(dp), allocatable :: tmpState(:,:)
    real(dp), allocatable :: tmpEigen(:)
    real(dp), allocatable :: tmpEn(:)

    integer :: ist, jst, nstates, nActPair

    nActPair = Na * (Na - 1) / 2
    nstates = size(energy,dim=1)

    allocate(Wab(nActPair,2))
    allocate(StateCoup(nstates,nstates))
    allocate(tmpOver(nstates,nstates))
    allocate(tmpState(nstates,nstates))
    allocate(tmpEigen(nstates))
    allocate(tmpEn(nstates))

    call getLagrangians_(env, denseDesc, neighbourList, nNeighbourSK, &
        & iSparseStart, img2CentCell, eigenvecs(:,:,1), hamSqrL, &
        & hamSpL, weight, fillingL, Nc, Na, Lpaired, tRangeSep, Wab)

    if (tSSR22) then
      call getStateCoup22_(Wab, FONs, StateCoup)
    else if (tSSR44) then
      call error("SSR(4,4) is not implemented yet")
    end if

    ! diagonalize the state energies
    ! obtain SSR energies & state-interaction term
    tmpOver(:,:) = 0.0_dp
    do ist = 1, nstates
      tmpOver(ist,ist) = 1.0_dp
    end do
    tmpEigen(:) = 0.0_dp

    tmpState(:,:) = 0.0_dp
    do ist = 1, nstates
      do jst = 1, nstates
        if (ist == jst) then
          tmpState(ist,jst) = energy(ist)
        else
          tmpState(ist,jst) = StateCoup(ist,jst)
        end if
      end do
    end do

    ! save state energies to print information
    tmpEn(:) = energy
    if (useSSR == 1) then
      call diagDenseMtx(electronicSolver, 'V', tmpState, tmpOver, tmpEigen)
      eigvecsSSR(:,:) = tmpState
      energy(:) = tmpEigen
    end if

    ! print state energies and couplings
    call printReksSSRInfo(Wab, tmpEn, StateCoup, energy, eigvecsSSR, &
        & Elevel, useSSR, Na, tSSR22, tSSR44)

  end subroutine solveSecularEqn


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!! Private routines
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Calculate filling of each microstate in REKS(2,2)
  subroutine getFillingL22_(Nc, fillingL)
    
    !> Number of core orbitals
    integer, intent(in) :: Nc

    !> Filling for each microstate
    real(dp), intent(out) :: fillingL(:,:,:)

    integer :: iL, iSpin, ii, nSpin, Lmax

    nSpin = size(fillingL,dim=2)
    Lmax = size(fillingL,dim=3)

    ! Filling of core orbitals
    do iL = 1, Lmax
      do iSpin = 1, nSpin
        do ii = 1, Nc
          fillingL(ii,iSpin,iL) = 1.0_dp
        end do
      end do
    end do

    ! Filling of active orbitals for REKS(2,2) case
    ! 1 = a: up + down
    fillingL(Nc+1,1,1) = 1.0_dp
    fillingL(Nc+1,2,1) = 1.0_dp
    ! 2 = b: up + down
    fillingL(Nc+2,1,2) = 1.0_dp
    fillingL(Nc+2,2,2) = 1.0_dp
    ! 3 = a: up, b: down
    fillingL(Nc+1,1,3) = 1.0_dp
    fillingL(Nc+2,2,3) = 1.0_dp
    ! 4 = a: down, b: up
    fillingL(Nc+2,1,4) = 1.0_dp
    fillingL(Nc+1,2,4) = 1.0_dp
    ! 5 = a: up, b: up
    fillingL(Nc+1,1,5) = 1.0_dp
    fillingL(Nc+2,1,5) = 1.0_dp
    ! 6 = a: down, b: down
    fillingL(Nc+1,2,6) = 1.0_dp
    fillingL(Nc+2,2,6) = 1.0_dp

  end subroutine getFillingL22_


  !> Make (2e,2o) weights, C_L used in SA-REKS
  subroutine getWeightL22_(FONs, delta, SAweight, weightL, weight)

    !> Fractional occupation numbers of active orbitals
    real(dp), intent(in) :: FONs(:,:)

    !> Smoothing factor used in FON optimization
    real(dp), intent(in) :: delta

    !> Weights used in state-averaging
    real(dp), intent(in) :: SAweight(:)

    !> Weight for each microstate per state
    real(dp), intent(out) :: weightL(:,:)

    !> Weight of each microstate for state to be optimized; weight = weightL * SAweight
    real(dp), intent(out) :: weight(:)

    integer :: iL, Lmax, ist, SAstates, nstates
    real(dp) :: n_a, n_b, fac

    Lmax = size(weightL,dim=2)
    SAstates = size(SAweight,dim=1)
    nstates = size(weightL,dim=1)

    n_a = FONs(1,1)
    n_b = FONs(2,1)

    fac = getFactor(n_a, n_b, delta)

    weightL(1,1) = 0.5_dp*n_a
    weightL(1,2) = 0.5_dp*n_b
    weightL(1,3) = fac
    weightL(1,4) = fac
    weightL(1,5) = -fac
    weightL(1,6) = -fac

    if(nstates >= 2) then
      weightL(2,1) = 0.0_dp
      weightL(2,2) = 0.0_dp
      weightL(2,3) = 1.0_dp
      weightL(2,4) = 1.0_dp
      weightL(2,5) = -0.5_dp
      weightL(2,6) = -0.5_dp
    end if

    if(nstates >= 3) then
      weightL(3,1) = 0.5_dp*n_b
      weightL(3,2) = 0.5_dp*n_a
      weightL(3,3) = -fac
      weightL(3,4) = -fac
      weightL(3,5) = fac
      weightL(3,6) = fac
    end if

    ! Decide which state will be optimized
    ! SAstates = 1 -> PPS state is optimized
    ! SAstates = 2 -> (PPS+OSS)/2 state (averaged state) is optimized
    weight(:) = 0.0_dp
    do iL = 1, Lmax
      do ist = 1, SAstates
        weight(iL) = weight(iL) + SAweight(ist)*weightL(ist,iL)
      end do
    end do

  end subroutine getWeightL22_


  !> Swap active orbitals when fa < fb in REKS(2,2) case
  subroutine MOswap22_(eigenvecs, SAweight, FONs, Efunction, Nc)

    !> eigenvectors
    real(dp), intent(inout) :: eigenvecs(:,:)

    !> Weights used in state-averaging
    real(dp), intent(in) :: SAweight(:)

    !> Fractional occupation numbers of active orbitals
    real(dp), intent(in) :: FONs(:,:)

    !> Minimized energy functional
    integer, intent(in) :: Efunction

    !> Number of core orbitals
    integer, intent(in) :: Nc

    real(dp), allocatable :: tmpMO(:)

    real(dp) :: n_a, n_b, fa, fb
    integer :: nOrb

    nOrb = size(eigenvecs,dim=1)

    n_a = FONs(1,1)
    n_b = FONs(2,1)

    allocate(tmpMO(nOrb))

    if (Efunction == 1) then
      ! REKS charge
      fa = n_a * 0.5_dp
      fb = n_b * 0.5_dp
    else if (Efunction == 2) then
      ! SA-REKS charge
      fa = (n_a*SAweight(1) + SAweight(2)) * 0.5_dp
      fb = (n_b*SAweight(1) + SAweight(2)) * 0.5_dp
    end if

    if (fa < fb) then
      write(stdOut,'(A6,F9.6,A20,I4,A8,I4,A8)') " fa = ", fa, &
          & ", MO swap between a(", Nc+1, ") and b(", Nc+2, ") occurs"
      tmpMO(:) = eigenvecs(:,Nc+1)
      eigenvecs(:,Nc+1) = eigenvecs(:,Nc+2)
      eigenvecs(:,Nc+2) = tmpMO
    end if

  end subroutine MOswap22_


  !> Calculate filling for minimzed state with optimized FONs in REKS(2,2)
  subroutine getFilling22_(filling, SAweight, FONs, Efunction, Nc)

    !> occupations (level)
    real(dp), intent(out) :: filling(:)

    !> Weights used in state-averaging
    real(dp), intent(in) :: SAweight(:)

    !> Fractional occupation numbers of active orbitals
    real(dp), intent(in) :: FONs(:,:)

    !> Minimized energy functional
    integer, intent(in) :: Efunction

    !> Number of core orbitals
    integer, intent(in) :: Nc

    real(dp) :: n_a, n_b
    integer :: ii

    n_a = FONs(1,1)
    n_b = FONs(2,1)

    filling(:) = 0.0_dp
    do ii = 1, Nc
      filling(ii) = 2.0_dp
    end do
    if (Efunction == 1) then
      ! REKS charge
      filling(Nc+1) = n_a
      filling(Nc+2) = n_b
    else if (Efunction == 2) then
      ! SA-REKS charge
      filling(Nc+1) = n_a*SAweight(1) + SAweight(2)
      filling(Nc+2) = n_b*SAweight(1) + SAweight(2)
    end if

  end subroutine getFilling22_


  !> Calculate Fc and Fa from Hamiltonian of each microstate
  subroutine getFockFcFa_(env, denseDesc, neighbourList, nNeighbourSK, &
      & iSparseStart, img2CentCell, hamSqrL, hamSpL, weight, fillingL, &
      & Nc, Na, Lpaired, tRangeSep, orbFON, Fc, Fa)

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

    !> state-averaged occupation numbers
    real(dp), intent(inout) :: orbFON(:)

    !> dense fock matrix for core orbitals
    real(dp), intent(out) :: Fc(:,:)

    !> dense fock matrix for active orbitals
    real(dp), intent(out) :: Fa(:,:,:)

    !> Dense Hamiltonian matrix for each microstate
    real(dp), intent(inout) :: hamSqrL(:,:,:,:)

    !> Sparse Hamiltonian matrix for each microstate
    real(dp), intent(in) :: hamSpL(:,:,:)

    !> Weight of each microstate for state to be optimized; weight = weightL * SAweight
    real(dp), intent(in) :: weight(:)

    !> Filling for each microstate
    real(dp), intent(in) :: fillingL(:,:,:)

    !> Number of core orbitals
    integer, intent(in) :: Nc

    !> Number of active orbitals
    integer, intent(in) :: Na

    !> Number of spin-paired microstates
    integer, intent(in) :: Lpaired

    !> Whether to run a range separated calculation
    logical, intent(in) :: tRangeSep

    real(dp), allocatable :: tmpHam(:,:)

    integer :: iL, Lmax, nOrb

    nOrb = size(Fc,dim=1)
    Lmax = size(weight,dim=1)

    if (.not. tRangeSep) then
      allocate(tmpHam(nOrb,nOrb))
    end if

    call fockFON_(fillingL, weight, orbFON)

    Fc(:,:) = 0.0_dp
    Fa(:,:,:) = 0.0_dp
    do iL = 1, Lmax

      if (.not. tRangeSep) then
        tmpHam(:,:) = 0.0_dp
        ! convert from sparse to dense for hamSpL in AO basis
        ! hamSpL has (my_ud) component
        call env%globalTimer%startTimer(globalTimers%sparseToDense)
        call unpackHS(tmpHam, hamSpL(:,1,iL), neighbourList%iNeighbour, nNeighbourSK, &
            & denseDesc%iAtomStart, iSparseStart, img2CentCell)
        call env%globalTimer%stopTimer(globalTimers%sparseToDense)
        call blockSymmetrizeHS(tmpHam, denseDesc%iAtomStart)
      end if

      ! compute the Fock operator with core, r, s orbitals in AO basis
      if (tRangeSep) then
        call fockFcAO_(hamSqrL(:,:,1,iL), weight, Lpaired, iL, Fc)
        call fockFaAO_(hamSqrL(:,:,1,iL), weight, fillingL, orbFON, &
            & Nc, Na, Lpaired, iL, Fa)
      else
        call fockFcAO_(tmpHam, weight, Lpaired, iL, Fc)
        call fockFaAO_(tmpHam, weight, fillingL, orbFON, &
            & Nc, Na, Lpaired, iL, Fa)
      end if

    end do

  end subroutine getFockFcFa_


  !> Calculate pseudo-fock matrix from Fc and Fa
  subroutine getPseudoFock_(Fc, Fa, orbFON, Nc, Na, fock)

    !> dense pseudo-fock matrix
    real(dp), intent(out) :: fock(:,:)

    !> dense fock matrix for core orbitals
    real(dp), intent(in) :: Fc(:,:)

    !> dense fock matrix for active orbitals
    real(dp), intent(in) :: Fa(:,:,:)

    !> state-averaged occupation numbers
    real(dp), intent(in) :: orbFON(:)

    !> Number of core orbitals
    integer, intent(in) :: Nc

    !> Number of active orbitals
    integer, intent(in) :: Na

    real(dp) :: res
    integer :: ii, jj, ind1, ind2, nOrb

    nOrb = size(fock,dim=1)

    fock(:,:) = 0.0_dp
    do ii = 1, Nc
      do jj = ii, Nc
        fock(jj,ii) = Fc(ii,jj)
      end do
      do jj = Nc + 1, Nc + Na
        ind1 = jj - Nc
        call fockFijMO_(res, Fc(ii,jj), Fa(ii,jj,ind1), &
            & orbFON(ii), orbFON(jj))
        fock(jj,ii) = res
      end do
      do jj = Nc + Na + 1, nOrb
        fock(jj,ii) = Fc(ii,jj)
      end do
    end do

    do jj = Nc + Na + 1, nOrb
      do ii = Nc + 1, Nc + Na
        ind1 = ii - Nc
        call fockFijMO_(res, Fc(jj,ii), Fa(jj,ii,ind1), &
            & orbFON(jj), orbFON(ii))
        fock(jj,ii) = res
      end do
      do ii = jj, nOrb
        fock(ii,jj) = Fc(jj,ii)
      end do
    end do

    do ii = Nc + 1, Nc + Na
      ind1 = ii - Nc
      do jj = Nc + 1, Nc + Na
        ind2 = jj - Nc
        if (ii == jj) then
          fock(jj,ii) = Fa(ii,jj,ind1)
        else
          call fockFijMO_(res, Fa(ii,jj,ind1), Fa(ii,jj,ind2), &
              & orbFON(ii), orbFON(jj))
          fock(jj,ii) = res
        end if
      end do
    end do

    call symmetrizeHS(fock)

  end subroutine getPseudoFock_


  !> Avoid changing the order of MOs
  !> Required number of cycles increases as the number of shift increases
  subroutine levelShifting_(fock, shift, Nc, Na)

    !> dense pseudo-fock matrix
    real(dp), intent(inout) :: fock(:,:)

    !> Shift value in SCC cycle
    real(dp), intent(in) :: shift

    !> Number of core orbitals
    integer, intent(in) :: Nc

    !> Number of active orbitals
    integer, intent(in) :: Na

    integer :: nOrb, ind, ii

    nOrb = size(fock,dim=1)

    do ii = Nc + 1, Nc + Na
      ind = ii - Nc
      fock(ii,ii) = fock(ii,ii) + dble(ind) * shift
    end do

    do ii = Nc + Na + 1, nOrb
      ind = Na + 1
      fock(ii,ii) = fock(ii,ii) + dble(ind) * shift
    end do

  end subroutine levelShifting_


  !> Calculate state-averaged FONs
  subroutine fockFON_(fillingL, weight, orbFON)

    !> state-averaged occupation numbers
    real(dp), intent(out) :: orbFON(:)

    !> Filling for each microstate
    real(dp), intent(in) :: fillingL(:,:,:)

    !> Weight of each microstate for state to be optimized; weight = weightL * SAweight
    real(dp), intent(in) :: weight(:)

    integer :: Lmax, iL

    Lmax = size(weight,dim=1)

    orbFON(:) = 0.0_dp
    do iL = 1, Lmax
      orbFON(:) = orbFON(:) + 0.5_dp * weight(iL) * &
          & ( fillingL(:,1,iL) + fillingL(:,2,iL) )
    end do

  end subroutine fockFON_


  !> Calculate fock matrix for core orbitals in AO basis
  subroutine fockFcAO_(hamSqr, weight, Lpaired, iL, Fc)

    !> dense fock matrix for core orbitals
    real(dp), intent(out) :: Fc(:,:)

    !> Dense Hamiltonian matrix for each microstate
    real(dp), intent(in) :: hamSqr(:,:)

    !> Weight of each microstate for state to be optimized; weight = weightL * SAweight
    real(dp), intent(in) :: weight(:)

    !> Number of spin-paired microstates
    integer, intent(in) :: Lpaired

    !> current index in loop L
    integer, intent(in) :: iL

    if (iL <= Lpaired) then
      Fc(:,:) = Fc + 0.5_dp * hamSqr * &
          & ( weight(iL) + weight(iL) )
    else
      if (mod(iL,2) == 1) then
        Fc(:,:) = Fc + 0.5_dp * hamSqr * &
            & ( weight(iL) + weight(iL+1) )
      else
        Fc(:,:) = Fc + 0.5_dp * hamSqr * &
            & ( weight(iL) + weight(iL-1) )
      end if
    end if

  end subroutine fockFcAO_


  !> Calculate fock matrix for active orbitals in AO basis
  subroutine fockFaAO_(hamSqr, weight, fillingL, orbFON, Nc, Na, &
      & Lpaired, iL, Fa)

    !> dense fock matrix for active orbitals
    real(dp), intent(out) :: Fa(:,:,:)

    !> Dense Hamiltonian matrix for each microstate
    real(dp), intent(in) :: hamSqr(:,:)

    !> Weight of each microstate for state to be optimized; weight = weightL * SAweight
    real(dp), intent(in) :: weight(:)

    !> Filling for each microstate
    real(dp), intent(in) :: fillingL(:,:,:)

    !> state-averaged occupation numbers
    real(dp), intent(in) :: orbFON(:)

    !> Number of core orbitals
    integer, intent(in) :: Nc

    !> Number of active orbitals
    integer, intent(in) :: Na

    !> Number of spin-paired microstates
    integer, intent(in) :: Lpaired

    !> current index in loop L
    integer, intent(in) :: iL

    integer :: ind, ind_a

    do ind = 1, Na
      ind_a = Nc + ind
      if (iL <= Lpaired) then
        Fa(:,:,ind) = Fa(:,:,ind) + 0.5_dp * fillingL(ind_a,1,iL) * &
            & ( weight(iL) + weight(iL) ) * hamSqr / orbFON(ind_a)
      else
        if (mod(iL,2) == 1) then
          Fa(:,:,ind) = Fa(:,:,ind) + 0.5_dp * fillingL(ind_a,1,iL) * &
              & ( weight(iL) + weight(iL+1) ) * hamSqr / orbFON(ind_a)
        else
          Fa(:,:,ind) = Fa(:,:,ind) + 0.5_dp * fillingL(ind_a,1,iL) * &
              & ( weight(iL) + weight(iL-1) ) * hamSqr / orbFON(ind_a)
        end if
      end if
    end do

  end subroutine fockFaAO_


  !> Calculate pseudo-fock off-diagonal element in MO basis
  subroutine fockFijMO_(res, fock_i, fock_j, f_i, f_j)

    !> temporary pseudo-fock value
    real(dp), intent(out) :: res

    !> temporary Fc or Fa values
    real(dp), intent(in) :: fock_i, fock_j

    !> temporary orbFON values
    real(dp), intent(in) :: f_i, f_j

    real(dp) :: eps = 1.0E-3_dp

    res = 0.0_dp
    if (abs(f_j-f_i) .LT. eps) then
      if (f_j >= f_i) then
        res = ( f_j*fock_j - f_i*fock_i )
      else
        res = -( f_j*fock_j - f_i*fock_i )
      end if
    else
      res = ( f_j*fock_j - f_i*fock_i ) / (f_j - f_i)
    end if

  end subroutine fockFijMO_


  !> Calculate converged Lagrangian values
  subroutine getLagrangians_(env, denseDesc, neighbourList, nNeighbourSK, &
      & iSparseStart, img2CentCell, eigenvecs, hamSqrL, hamSpL, weight, &
      & fillingL, Nc, Na, Lpaired, tRangeSep, Wab)

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
    real(dp), intent(in) :: eigenvecs(:,:)

    !> Dense Hamiltonian matrix for each microstate
    real(dp), intent(inout) :: hamSqrL(:,:,:,:)

    !> Sparse Hamiltonian matrix for each microstate
    real(dp), intent(inout) :: hamSpL(:,:,:)

    !> Weight of each microstate for state to be optimized; weight = weightL * SAweight
    real(dp), intent(in) :: weight(:)

    !> Filling for each microstate
    real(dp), intent(in) :: fillingL(:,:,:)

    !> Number of core orbitals
    integer, intent(in) :: Nc

    !> Number of active orbitals
    integer, intent(in) :: Na

    !> Number of spin-paired microstates
    integer, intent(in) :: Lpaired

    !> Whether to run a range separated calculation
    logical, intent(in) :: tRangeSep

    !> converged Lagrangian values within active space
    real(dp), intent(out) :: Wab(:,:)

    real(dp), allocatable :: tmpHam(:,:)
    real(dp), allocatable :: tmpHamL(:,:,:)

    integer :: nOrb, iL, Lmax
    integer :: ia, ib, ist, nActPair

    nOrb = size(eigenvecs,dim=1)
    Lmax = size(fillingL,dim=3)
    nActPair = Na * (Na - 1) / 2

    if (.not. tRangeSep) then
      allocate(tmpHam(nOrb,nOrb))
    end if
    allocate(tmpHamL(nActPair,1,Lmax))

    tmpHamL(:,:,:) = 0.0_dp
    do ist = 1, nActPair

      call getTwoIndices(Na, ist, ia, ib, 1)

      do iL = 1, Lmax

        if (tRangeSep) then
          ! convert hamSqrL from AO basis to MO basis
          ! hamSqrL has (my_ud) component
          if (ist == 1) then
            call matAO2MO(hamSqrL(:,:,1,iL), eigenvecs)
          end if
          tmpHamL(ist,1,iL) = hamSqrL(Nc+ia,Nc+ib,1,iL)
        else
          tmpHam(:,:) = 0.0_dp
          ! convert from sparse to dense for hamSpL in AO basis
          ! hamSpL has (my_ud) component
          call env%globalTimer%startTimer(globalTimers%sparseToDense)
          call unpackHS(tmpHam, hamSpL(:,1,iL), &
              & neighbourList%iNeighbour, nNeighbourSK, &
              & denseDesc%iAtomStart, iSparseStart, img2CentCell)
          call env%globalTimer%stopTimer(globalTimers%sparseToDense)
          call blockSymmetrizeHS(tmpHam, denseDesc%iAtomStart)
          ! convert tmpHam from AO basis to MO basis
          call matAO2MO(tmpHam, eigenvecs)
          ! save F_{L,ab}^{\sigma} in MO basis
          tmpHamL(ist,1,iL) = tmpHam(Nc+ia,Nc+ib)
        end if

      end do

      ! calculate the Lagrangian eps_{ab} and state-interaction term
      Wab(ist,1) = 0.0_dp
      Wab(ist,2) = 0.0_dp
      do iL = 1, Lmax
        if (iL <= Lpaired) then
          Wab(ist,1) = Wab(ist,1) + fillingL(Nc+ia,1,iL) * &
              & tmpHamL(ist,1,iL) * ( weight(iL) + weight(iL) )
          Wab(ist,2) = Wab(ist,2) + fillingL(Nc+ib,1,iL) * &
              & tmpHamL(ist,1,iL) * ( weight(iL) + weight(iL) )
        else
          if (mod(iL,2) == 1) then
            Wab(ist,1) = Wab(ist,1) + fillingL(Nc+ia,1,iL) * &
                & tmpHamL(ist,1,iL) * ( weight(iL) + weight(iL+1) )
            Wab(ist,2) = Wab(ist,2) + fillingL(Nc+ib,1,iL) * &
                & tmpHamL(ist,1,iL) * ( weight(iL) + weight(iL+1) )
          else
            Wab(ist,1) = Wab(ist,1) + fillingL(Nc+ia,1,iL) * &
                & tmpHamL(ist,1,iL) * ( weight(iL) + weight(iL-1) )
            Wab(ist,2) = Wab(ist,2) + fillingL(Nc+ib,1,iL) * &
                & tmpHamL(ist,1,iL) * ( weight(iL) + weight(iL-1) )
          end if
        end if
      end do

    end do

  end subroutine getLagrangians_


  !> calculate state-interaction terms between SA-REKS states in (2,2) case
  subroutine getStateCoup22_(Wab, FONs, StateCoup)

    !> converged Lagrangian values within active space
    real(dp), intent(in) :: Wab(:,:)

    !> Fractional occupation numbers of active orbitals
    real(dp), intent(in) :: FONs(:,:)

    !> state-interaction term between SA-REKS states
    real(dp), intent(out) :: StateCoup(:,:)

    real(dp) :: n_a, n_b
    integer :: nstates

    n_a = FONs(1,1)
    n_b = FONs(2,1)
    nstates = size(StateCoup,dim=1)

    StateCoup(:,:) = 0.0_dp
    StateCoup(1,2) = sqrt(n_a) * Wab(1,1) - sqrt(n_b) * Wab(1,1)
    StateCoup(2,1) = StateCoup(1,2)
    if (nstates == 3) then
      StateCoup(2,3) = sqrt(n_a) * Wab(1,1) + sqrt(n_b) * Wab(1,1)
      StateCoup(3,2) = StateCoup(2,3)
    end if

  end subroutine getStateCoup22_


  !> Calculate factor from n_a, n_b, and delta for certain active orbital set
  function getFactor(n_a, n_b, delta) result(factor)

    !> Fractional occupation numbers of active orbitals
    real(dp), intent(in) :: n_a, n_b

    !> Smoothing factor used in FON optimization
    real(dp), intent(in) :: delta

    !> factor of n_a and n_b
    real(dp) :: factor

    factor = -0.5_dp*(n_a*n_b)**&
        & (1.0_dp-0.5_dp*(n_a*n_b+delta)/(1.0_dp+delta))

  end function getFactor


end module dftbp_reksen
