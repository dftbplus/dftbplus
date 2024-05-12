!--------------------------------------------------------------------------------------------------!
!  DFTB+: general package for performing fast atomistic simulations                                !
!  Copyright (C) 2006 - 2023  DFTB+ developers group                                               !
!                                                                                                  !
!  See the LICENSE file for terms of usage and distribution.                                       !
!--------------------------------------------------------------------------------------------------!

#:include 'common.fypp'

!> REKS and SI-SA-REKS formulation in DFTB as developed by Lee et al.
!>
!> The functionality of the module has some limitation:
!> * Third order does not work.
!> * Periodic system do not work yet apart from Gamma point.
!> * Orbital potentials or spin-orbit or external E-field does not work yet.
!> * Only for closed shell system.
!> * Onsite corrections are not included in this version
module dftbp_reks_reksen
  use dftbp_common_accuracy, only : dp
  use dftbp_common_environment, only : globalTimers, TEnvironment
  use dftbp_common_globalenv, only : stdOut
  use dftbp_dftb_energytypes, only : TEnergies
  use dftbp_dftb_periodic, only : TNeighbourList
  use dftbp_dftb_sparse2dense, only : unpackHS
  use dftbp_elecsolvers_elecsolvers, only: TElectronicSolver
  use dftbp_io_message, only : error
  use dftbp_math_blasroutines, only : gemm
  use dftbp_math_eigensolver, only : heev
  use dftbp_math_matrixops, only : adjointLowerTriangle
  use dftbp_reks_rekscommon, only : getTwoIndices, matAO2MO
  use dftbp_reks_reksio, only : printReksSSRInfo
  use dftbp_reks_reksvar, only : TReksCalc, reksTypes
  use dftbp_type_densedescr, only : TDenseDescr
#:if WITH_MPI
  use dftbp_extlibs_mpifx, only : MPI_SUM, mpifx_allreduceip
#:endif
#:if WITH_SCALAPACK
  use dftbp_dftb_sparse2dense, only : unpackHSRealBlacs
  use dftbp_reks_rekscommon, only : matAO2MOBlacs
  use dftbp_extlibs_scalapackfx, only : CSRC_, RSRC_, MB_, NB_, scalafx_indxl2g,&
      & scalafx_psyev, pblasfx_pgemm, blocklist, size
#:endif

  implicit none

  private
  public :: constructMicrostates, calcWeights
  public :: activeOrbSwap, getFilling, calcSaReksEnergy
  public :: getFockandDiag, guessNewEigvecs
  public :: adjustEigenval, solveSecularEqn
  public :: setReksTargetEnergy

  contains

  !> Construct L, spin dependent microstates from identical KS orbitals
  subroutine constructMicrostates(this)

    !> data type for REKS
    type(TReksCalc), intent(inout) :: this

    select case (this%reksAlg)
    case (reksTypes%noReks)
    case (reksTypes%ssr22)
      call getFillingL22_(this%Nc, this%fillingL)
    case (reksTypes%ssr44)
      call error("SSR(4,4) is not implemented yet")
    end select

  end subroutine constructMicrostates


  !> Calculate the weight of each microstate for current cycle, C_L
  subroutine calcWeights(this)

    !> data type for REKS
    type(TReksCalc), intent(inout) :: this

    select case (this%reksAlg)
    case (reksTypes%noReks)
    case (reksTypes%ssr22)
      call getWeightL22_(this%FONs, this%delta, this%SAweight, this%weightL, this%weight)
    case (reksTypes%ssr44)
      call error("SSR(4,4) is not implemented yet")
    end select

  end subroutine calcWeights


  !> Swap the active orbitals for feasible occupation in REKS
  subroutine activeOrbSwap(env, denseDesc, this, eigenvecs)

    !> Environment settings
    type(TEnvironment), intent(inout) :: env

    !> Dense matrix descriptor
    type(TDenseDescr), intent(in) :: denseDesc

    !> data type for REKS
    type(TReksCalc), intent(inout) :: this

    !> eigenvectors
    real(dp), intent(inout) :: eigenvecs(:,:)

    select case (this%reksAlg)
    case (reksTypes%noReks)
    case (reksTypes%ssr22)
      call MOswap22_(env, denseDesc, eigenvecs, this%SAweight, this%FONs, this%Efunction,&
          & this%Nc)
    case (reksTypes%ssr44)
      call error("SSR(4,4) is not implemented yet")
    end select

  end subroutine activeOrbSwap


  !> Calculate filling for minimzed state with optimized FONs
  subroutine getFilling(this, filling)

    !> data type for REKS
    type(TReksCalc), intent(inout) :: this

    !> occupations (level)
    real(dp), intent(out) :: filling(:)

    select case (this%reksAlg)
    case (reksTypes%noReks)
    case (reksTypes%ssr22)
      call getFilling22_(filling, this%SAweight, this%FONs, this%Efunction, this%Nc)
    case (reksTypes%ssr44)
      call error("SSR(4,4) is not implemented yet")
    end select

  end subroutine getFilling


  !> Calculate the energy of SA-REKS states and averaged state
  subroutine calcSaReksEnergy(this, energy)

    !> data type for REKS
    type(TReksCalc), intent(inout) :: this

    !> Energy terms in the system
    type(TEnergies), intent(inout) :: energy

    integer :: ist

    ! Compute the energy contributions for target SA-REKS state
    ! electronic energy = nonSCC + scc + spin + 3rd + fock
    energy%EnonSCC = sum(this%weightL(this%rstate,:)*this%enLnonSCC(:))
    energy%Escc = sum(this%weightL(this%rstate,:)*this%enLscc(:))
    energy%Espin = sum(this%weightL(this%rstate,:)*this%enLspin(:))
    if (this%t3rd) then
      energy%e3rd = sum(this%weightL(this%rstate,:)*this%enL3rd(:))
    end if
    if (this%isHybridXc) then
      energy%Efock = sum(this%weightL(this%rstate,:)*this%enLfock(:))
    end if
    if (this%isDispersion) then
      energy%Edisp = sum(this%weightL(this%rstate,:)*this%enLdisp(:))
    end if

    energy%Eelec = energy%EnonSCC + energy%Escc + energy%Espin + &
        & energy%e3rd + energy%Efock
    energy%Etotal = energy%Eelec + energy%Erep + energy%Edisp

    ! Compute the total energy for SA-REKS states
    do ist = 1, this%nstates
      this%energy(ist) = sum(this%weightL(ist,:)*this%enLtot(:))
    end do

!    if (abs(energy%Etotal - this%energy(this%rstate)) >= epsilon(1.0_dp)) then
    if (abs(energy%Etotal - this%energy(this%rstate)) >= 1.0e-8_dp) then
      call error("Wrong energy contribution for target SA-REKS state")
    end if

    ! In this step Eavg becomes the energy of averaged state
    ! From this energy we can check the variational principle
    energy%Eavg = 0.0_dp
    do ist = 1, this%SAstates
      energy%Eavg = energy%Eavg + this%SAweight(ist) * this%energy(ist)
    end do

  end subroutine calcSaReksEnergy


  !> Make pseudo-fock operator with Hamiltonian of each microstate
  !> and diagonalize the fock matrix
  subroutine getFockandDiag(env, denseDesc, neighbourList, &
      & nNeighbourSK, iSparseStart, img2CentCell, eigenvecs, &
      & electronicSolver, eigen, this)

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

    !> data type for REKS
    type(TReksCalc), intent(inout) :: this

    real(dp), allocatable :: orbFON(:)
    real(dp), allocatable :: tmpMat(:,:)

    integer :: ii, nOrb

    nOrb = denseDesc%fullSize

    allocate(orbFON(nOrb))
  #:if not WITH_SCALAPACK
    allocate(tmpMat(nOrb,nOrb))
  #:endif

    call getFockFcFa_(env, denseDesc, neighbourList, nNeighbourSK, &
        & iSparseStart, img2CentCell, this%hamSqrL, this%hamSpL, this%weight, &
        & this%fillingL, this%Nc, this%Na, this%Lpaired, this%isHybridXc, &
        & orbFON, this%fockFc, this%fockFa)

  #:if WITH_SCALAPACK
    call matAO2MOBlacs(this%fockFc, denseDesc%blacsOrbSqr, eigenvecs(:,:,1))
    do ii = 1, this%Na
      call matAO2MOBlacs(this%fockFa(:,:,ii), denseDesc%blacsOrbSqr, eigenvecs(:,:,1))
    end do
  #:else
    call matAO2MO(this%fockFc, eigenvecs(:,:,1))
    do ii = 1, this%Na
      call matAO2MO(this%fockFa(:,:,ii), eigenvecs(:,:,1))
    end do
  #:endif

  #:if WITH_SCALAPACK
    call getPseudoFockBlacs_(env, denseDesc, this%fockFc, this%fockFa, orbFON, this%Nc,&
        & this%Na, this%fock)
  #:else
    call getPseudoFock_(this%fockFc, this%fockFa, orbFON, this%Nc, this%Na, this%fock)
  #:endif

  #:if WITH_SCALAPACK
    call levelShiftingBlacs_(env, denseDesc, this%fock, this%shift, this%Nc, this%Na)
  #:else
    call levelShifting_(this%fock, this%shift, this%Nc, this%Na)
  #:endif

    ! Diagonalize the pesudo-Fock matrix
    eigen(:,1,1) = 0.0_dp
    call env%globalTimer%startTimer(globalTimers%diagonalization)
  #:if WITH_SCALAPACK
    call scalafx_psyev(this%fock, denseDesc%blacsOrbSqr, eigen(:,1,1), this%eigvecsFock,&
        & denseDesc%blacsOrbSqr, uplo="U", jobz="V")
  #:else
    tmpMat(:,:) = this%fock
    call heev(tmpMat, eigen(:,1,1), 'U', 'V')
    this%eigvecsFock(:,:) = tmpMat
  #:endif
    call env%globalTimer%stopTimer(globalTimers%diagonalization)

  end subroutine getFockandDiag


  !> guess new eigenvectors from Fock eigenvectors
  subroutine guessNewEigvecs(eigenvecs, desc, eigvecsFock)

    !> Eigenvectors on eixt
    real(dp), intent(inout) :: eigenvecs(:,:)

    !> BLACS matrix descriptor
    integer, intent(in) :: desc(:)

    !> eigenvectors from pesudo-fock matrix
    real(dp), intent(in) :: eigvecsFock(:,:)

    real(dp), allocatable :: tmpVec(:,:)
    integer :: nLocalRows, nLocalCols

    nLocalRows = size(eigvecsFock,dim=1)
    nLocalCols = size(eigvecsFock,dim=2)

    allocate(tmpVec(nLocalRows,nLocalCols))

    tmpVec(:,:) = 0.0_dp
  #:if WITH_SCALAPACK
    call pblasfx_pgemm(eigenvecs, desc, eigvecsFock, desc, tmpVec, desc)
  #:else
    call gemm(tmpVec, eigenvecs, eigvecsFock)
  #:endif
    eigenvecs(:,:) = tmpVec

  end subroutine guessNewEigvecs


  !> adjust the eigenvalues (eliminate shift values)
  subroutine adjustEigenval(this, eigen)

    !> data type for REKS
    type(TReksCalc), intent(inout) :: this

    !> eigenvalues
    real(dp), intent(inout) :: eigen(:,:,:)

    integer :: nOrb, ind, ii

    nOrb = size(eigen,dim=1)

    do ii = this%Nc + 1, this%Nc + this%Na
      ind = ii - this%Nc
      eigen(ii,1,1) = eigen(ii,1,1) - real(ind, dp) * this%shift
    end do

    do ii = this%Nc + this%Na + 1, nOrb
      ind = this%Na + 1
      eigen(ii,1,1) = eigen(ii,1,1) - real(ind, dp) * this%shift
    end do

  end subroutine adjustEigenval


  !> Solve secular equation with coupling element between SA-REKS states
  subroutine solveSecularEqn(env, denseDesc, neighbourList, &
      & nNeighbourSK, iSparseStart, img2CentCell, electronicSolver, &
      & eigenvecs, this)

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

    !> data type for REKS
    type(TReksCalc), intent(inout) :: this

    real(dp), allocatable :: Wab(:,:)
    real(dp), allocatable :: StateCoup(:,:)
    real(dp), allocatable :: tmpState(:,:)
    real(dp), allocatable :: tmpEigen(:)
    real(dp), allocatable :: tmpEn(:)

    integer :: ist, jst, nActPair

    nActPair = this%Na * (this%Na - 1) / 2

    allocate(Wab(nActPair,2))
    allocate(StateCoup(this%nstates,this%nstates))
    allocate(tmpState(this%nstates,this%nstates))
    allocate(tmpEigen(this%nstates))
    allocate(tmpEn(this%nstates))

    call getLagrangians_(env, denseDesc, neighbourList, nNeighbourSK, &
        & iSparseStart, img2CentCell, eigenvecs(:,:,1), this%hamSqrL, &
        & this%hamSpL, this%weight, this%fillingL, this%Nc, this%Na, &
        & this%Lpaired, this%isHybridXc, Wab)

    select case (this%reksAlg)
    case (reksTypes%noReks)
    case (reksTypes%ssr22)
      call getStateCoup22_(Wab, this%FONs, StateCoup)
    case (reksTypes%ssr44)
      call error("SSR(4,4) is not implemented yet")
    end select

    ! diagonalize the state energies
    ! obtain SSR energies & state-interaction term
    tmpEigen(:) = 0.0_dp

    tmpState(:,:) = 0.0_dp
    do ist = 1, this%nstates
      do jst = 1, this%nstates
        if (ist == jst) then
          tmpState(ist,jst) = this%energy(ist)
        else
          tmpState(ist,jst) = StateCoup(ist,jst)
        end if
      end do
    end do

    ! save state energies to print information
    tmpEn(:) = this%energy
    if (this%tSSR) then
      call env%globalTimer%startTimer(globalTimers%diagonalization)
      call heev(tmpState, tmpEigen, 'U', 'V')
      this%eigvecsSSR(:,:) = tmpState
      call env%globalTimer%stopTimer(globalTimers%diagonalization)
      this%energy(:) = tmpEigen
    end if

    ! print state energies and couplings
    call printReksSSRInfo(this, Wab, tmpEn, StateCoup)

  end subroutine solveSecularEqn


  !> Set correct final energy values for target state or microstate
  subroutine setReksTargetEnergy(this, energy, cellVol, pressure)

    !> data type for REKS
    type(TReksCalc), intent(in) :: this

    !> Energy terms in the system
    type(TEnergies), intent(inout) :: energy

    !> Unit cell volume
    real(dp), intent(in) :: cellVol

    !> External pressure
    real(dp), intent(in) :: pressure

    ! get correct energy values
    if (this%Lstate == 0) then

      ! get energy contributions for target state
      energy%Etotal = this%energy(this%rstate)
      if (this%nstates > 1) then
        energy%Eexcited = this%energy(this%rstate) - this%energy(1)
      else
        energy%Eexcited = 0.0_dp
      end if

    else

      ! get energy contributions for target microstate
      energy%EnonSCC = this%enLnonSCC(this%Lstate)
      energy%ESCC = this%enLSCC(this%Lstate)
      energy%Espin = this%enLspin(this%Lstate)
      if (this%t3rd) then
        energy%e3rd = this%enL3rd(this%Lstate)
      end if
      if (this%isHybridXc) then
        energy%Efock = this%enLfock(this%Lstate)
      end if
      if (this%isDispersion) then
        energy%Edisp = this%enLdisp(this%Lstate)
      end if

      energy%Eelec = energy%EnonSCC + energy%Escc + energy%Espin + &
          & energy%e3rd + energy%Efock
      energy%Etotal = energy%Eelec + energy%Erep + energy%Edisp
      energy%Eexcited = 0.0_dp

!      if (abs(energy%Etotal - this%enLtot(this%Lstate)) > epsilon(1.0_dp)) then
      if (abs(energy%Etotal - this%enLtot(this%Lstate)) > 1.0e-8_dp) then
        call error("Wrong energy contribution for target microstate")
      end if

    end if

    ! REKS is not affected by filling, so TS becomes 0
    energy%EMermin = energy%Etotal
    ! extrapolated to 0 K
    energy%Ezero = energy%Etotal
    energy%EGibbs = energy%EMermin + cellVol * pressure
    energy%EForceRelated = energy%EGibbs

  end subroutine setReksTargetEnergy


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

    fillingL(:,:,:) = 0.0_dp

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
  subroutine MOswap22_(env, denseDesc, eigenvecs, SAweight, FONs, Efunction, Nc)

    !> Environment settings
    type(TEnvironment), intent(inout) :: env

    !> Dense matrix descriptor
    type(TDenseDescr), intent(in) :: denseDesc

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

    real(dp), allocatable :: tmpMO(:,:)
    real(dp), allocatable :: tmpMOa(:), tmpMOb(:)

    real(dp) :: n_a, n_b, fa, fb
    integer :: nOrb, nLocalRows, nLocalCols, iOrb, iGlob, jOrb, jGlob

    nOrb = denseDesc%fullSize
    nLocalRows = size(eigenvecs,dim=1)
    nLocalCols = size(eigenvecs,dim=2)

    n_a = FONs(1,1)
    n_b = FONs(2,1)

  #:if WITH_SCALAPACK
    allocate(tmpMO(nOrb,2))
    allocate(tmpMOa(nLocalRows))
    allocate(tmpMOb(nLocalRows))
  #:else
    allocate(tmpMO(nOrb,1))
  #:endif

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
    #:if WITH_SCALAPACK
      tmpMO(:,:) = 0.0_dp
      tmpMOa(:) = 0.0_dp
      tmpMOb(:) = 0.0_dp
      do jOrb = 1, nLocalCols
        jGlob = scalafx_indxl2g(jOrb, denseDesc%blacsOrbSqr(NB_), env%blacs%orbitalGrid%mycol,&
            & denseDesc%blacsOrbSqr(CSRC_), env%blacs%orbitalGrid%ncol)
        if (jGlob == Nc + 1) then
          tmpMOa(:) = eigenvecs(:,jOrb)
        else if (jGlob == Nc + 2) then
          tmpMOb(:) = eigenvecs(:,jOrb)
        end if
      end do
      do iOrb = 1, nLocalRows
        iGlob = scalafx_indxl2g(iOrb, denseDesc%blacsOrbSqr(MB_), env%blacs%orbitalGrid%myrow,&
            & denseDesc%blacsOrbSqr(RSRC_), env%blacs%orbitalGrid%nrow)
        tmpMO(iGlob,1) = tmpMOa(iOrb)
        tmpMO(iGlob,2) = tmpMOb(iOrb)
      end do
      call mpifx_allreduceip(env%mpi%globalComm, tmpMO, MPI_SUM)
      tmpMOa(:) = 0.0_dp
      tmpMOb(:) = 0.0_dp
      do iOrb = 1, nLocalRows
        iGlob = scalafx_indxl2g(iOrb, denseDesc%blacsOrbSqr(MB_), env%blacs%orbitalGrid%myrow,&
            & denseDesc%blacsOrbSqr(RSRC_), env%blacs%orbitalGrid%nrow)
        tmpMOb(iOrb) = tmpMO(iGlob,1)
        tmpMOa(iOrb) = tmpMO(iGlob,2)
      end do
      do jOrb = 1, nLocalCols
        jGlob = scalafx_indxl2g(jOrb, denseDesc%blacsOrbSqr(NB_), env%blacs%orbitalGrid%mycol,&
            & denseDesc%blacsOrbSqr(CSRC_), env%blacs%orbitalGrid%ncol)
        if (jGlob == Nc + 1) then
          eigenvecs(:,jOrb) = tmpMOa
        else if (jGlob == Nc + 2) then
          eigenvecs(:,jOrb) = tmpMOb
        end if
      end do
    #:else
      tmpMO(:,1) = eigenvecs(:,Nc+1)
      eigenvecs(:,Nc+1) = eigenvecs(:,Nc+2)
      eigenvecs(:,Nc+2) = tmpMO(:,1)
    #:endif
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
      & Nc, Na, Lpaired, isHybridXc, orbFON, Fc, Fa)

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
    real(dp), allocatable, intent(inout) :: hamSqrL(:,:,:,:)

    !> Sparse Hamiltonian matrix for each microstate
    real(dp), allocatable, intent(in) :: hamSpL(:,:,:)

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
    logical, intent(in) :: isHybridXc

    real(dp), allocatable :: tmpHam(:,:)

    integer :: iL, Lmax, nLocalRows, nLocalCols

    nLocalRows = size(Fc,dim=1)
    nLocalCols = size(Fc,dim=2)
    Lmax = size(weight,dim=1)

    if (.not. isHybridXc) then
      allocate(tmpHam(nLocalRows,nLocalCols))
    end if

    call fockFON_(fillingL, weight, orbFON)

    Fc(:,:) = 0.0_dp
    Fa(:,:,:) = 0.0_dp
    do iL = 1, Lmax

      if (.not. isHybridXc) then
        tmpHam(:,:) = 0.0_dp
        ! convert from sparse to dense for hamSpL in AO basis
        ! hamSpL has (my_ud) component
        call env%globalTimer%startTimer(globalTimers%sparseToDense)
      #:if WITH_SCALAPACK
        call unpackHSRealBlacs(env%blacs, hamSpL(:,1,iL), neighbourList%iNeighbour,&
            & nNeighbourSK, iSparseStart, img2CentCell, denseDesc, tmpHam)
      #:else
        call unpackHS(tmpHam, hamSpL(:,1,iL), neighbourList%iNeighbour, nNeighbourSK, &
            & denseDesc%iAtomStart, iSparseStart, img2CentCell)
        call adjointLowerTriangle(tmpHam)
      #:endif
        call env%globalTimer%stopTimer(globalTimers%sparseToDense)
      end if

      ! compute the Fock operator with core, a, b orbitals in AO basis
      if (isHybridXc) then
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


#:if WITH_SCALAPACK
  !> Calculate pseudo-fock matrix from Fc and Fa
  subroutine getPseudoFockBlacs_(env, denseDesc, Fc, Fa, orbFON, Nc, Na, fock)

    !> Environment settings
    type(TEnvironment), intent(inout) :: env

    !> Dense matrix descriptor
    type(TDenseDescr), intent(in) :: denseDesc

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

    !> dense pseudo-fock matrix
    real(dp), intent(out) :: fock(:,:)

    type(blocklist) :: blocksRow, blocksCol
    real(dp) :: res
    integer :: ii, iGlob, iLoc, iSize, blockSize1, indGlobRow, indLocRow
    integer :: jj, jGlob, jLoc, jSize, blockSize2, indGlobCol, indLocCol
    integer :: ind1, ind2

    fock(:,:) = 0.0_dp

    call blocksRow%init(env%blacs%orbitalGrid, denseDesc%blacsOrbSqr, "r")
    call blocksCol%init(env%blacs%orbitalGrid, denseDesc%blacsOrbSqr, "c")
    do ii = 1, size(blocksRow)
      call blocksRow%getblock(ii, iGlob, iLoc, blockSize1)
      do jj = 1, size(blocksCol)
        call blocksCol%getblock(jj, jGlob, jLoc, blockSize2)

        do iSize = 0, blockSize1 - 1
          indGlobRow = iGlob + iSize
          indLocRow = iLoc + iSize

          if (indGlobRow <= Nc) then

            do jSize = 0, blockSize2 - 1
              indGlobCol = jGlob + jSize
              indLocCol = jLoc + jSize

              if (indGlobCol <= Nc) then
                fock(indLocRow,indLocCol) = Fc(indLocRow,indLocCol)
              else if (indGlobCol > Nc .and. indGlobCol <= Nc + Na) then
                ind1 = indGlobCol - Nc
                call fockFijMO_(res, Fc(indLocRow,indLocCol), Fa(indLocRow,indLocCol,ind1),&
                    & orbFON(indGlobRow), orbFON(indGlobCol))
                fock(indLocRow,indLocCol) = res
              else
                fock(indLocRow,indLocCol) = Fc(indLocRow,indLocCol)
              end if

            end do

          else if (indGlobRow > Nc .and. indGlobRow <= Nc + Na) then

            do jSize = 0, blockSize2 - 1
              indGlobCol = jGlob + jSize
              indLocCol = jLoc + jSize

              if (indGlobCol <= Nc) then
                ind1 = indGlobRow - Nc
                call fockFijMO_(res, Fa(indLocRow,indLocCol,ind1), Fc(indLocRow,indLocCol),&
                    & orbFON(indGlobRow), orbFON(indGlobCol))
                fock(indLocRow,indLocCol) = res
              else if (indGlobCol > Nc .and. indGlobCol <= Nc + Na) then
                ind1 = indGlobRow - Nc
                ind2 = indGlobCol - Nc
                if (indGlobRow == indGlobCol) then
                  fock(indLocRow,indLocCol) = Fa(indLocRow,indLocCol,ind1)
                else
                  call fockFijMO_(res, Fa(indLocRow,indLocCol,ind1), Fa(indLocRow,indLocCol,ind2),&
                      & orbFON(indGlobRow), orbFON(indGlobCol))
                  fock(indLocRow,indLocCol) = res
                end if
              else
                ind1 = indGlobRow - Nc
                call fockFijMO_(res, Fa(indLocRow,indLocCol,ind1), Fc(indLocRow,indLocCol),&
                    & orbFON(indGlobRow), orbFON(indGlobCol))
                fock(indLocRow,indLocCol) = res
              end if

            end do

          else

            do jSize = 0, blockSize2 - 1
              indGlobCol = jGlob + jSize
              indLocCol = jLoc + jSize

              if (indGlobCol <= Nc) then
                fock(indLocRow,indLocCol) = Fc(indLocRow,indLocCol)
              else if (indGlobCol > Nc .and. indGlobCol <= Nc + Na) then
                ind1 = indGlobCol - Nc
                call fockFijMO_(res, Fc(indLocRow,indLocCol), Fa(indLocRow,indLocCol,ind1),&
                    & orbFON(indGlobRow), orbFON(indGlobCol))
                fock(indLocRow,indLocCol) = res
              else
                fock(indLocRow,indLocCol) = Fc(indLocRow,indLocCol)
              end if

            end do

          end if

        end do

      end do
    end do

  end subroutine getPseudoFockBlacs_
#:else
  !> Calculate pseudo-fock matrix from Fc and Fa
  subroutine getPseudoFock_(Fc, Fa, orbFON, Nc, Na, fock)

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

    !> dense pseudo-fock matrix
    real(dp), intent(out) :: fock(:,:)

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

    call adjointLowerTriangle(fock)

  end subroutine getPseudoFock_
#:endif


#:if WITH_SCALAPACK
  !> Avoid changing the order of MOs
  !> Required number of cycles increases as the number of shift increases
  subroutine levelShiftingBlacs_(env, denseDesc, fock, shift, Nc, Na)

    !> Environment settings
    type(TEnvironment), intent(inout) :: env

    !> Dense matrix descriptor
    type(TDenseDescr), intent(in) :: denseDesc

    !> dense pseudo-fock matrix
    real(dp), intent(inout) :: fock(:,:)

    !> Shift value in SCC cycle
    real(dp), intent(in) :: shift

    !> Number of core orbitals
    integer, intent(in) :: Nc

    !> Number of active orbitals
    integer, intent(in) :: Na

    real(dp) :: res
    integer :: ind, iOrb, iGlob, jOrb, jGlob

    do jOrb = 1, size(fock, dim=2)
      jGlob = scalafx_indxl2g(jOrb, denseDesc%blacsOrbSqr(NB_), env%blacs%orbitalGrid%mycol,&
          & denseDesc%blacsOrbSqr(CSRC_), env%blacs%orbitalGrid%ncol)
      do iOrb = 1, size(fock, dim=1)
        iGlob = scalafx_indxl2g(iOrb, denseDesc%blacsOrbSqr(MB_), env%blacs%orbitalGrid%myrow,&
            & denseDesc%blacsOrbSqr(RSRC_), env%blacs%orbitalGrid%nrow)

        if (iGlob == jGlob .and. iGlob > Nc) then
          ind = iGlob - Nc
          if (ind >= Na + 1) ind = Na + 1
          res = fock(iOrb,jOrb) + real(ind, dp) * shift
          fock(iOrb,jOrb) = res
        end if

      end do
    end do

  end subroutine levelShiftingBlacs_
#:else
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
      fock(ii,ii) = fock(ii,ii) + real(ind, dp) * shift
    end do

    do ii = Nc + Na + 1, nOrb
      ind = Na + 1
      fock(ii,ii) = fock(ii,ii) + real(ind, dp) * shift
    end do

  end subroutine levelShifting_
#:endif


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
    real(dp), intent(inout) :: Fc(:,:)

    !> Dense Hamiltonian matrix for each microstate
    real(dp), intent(in) :: hamSqr(:,:)

    !> Weight of each microstate for state to be optimized; weight = weightL * SAweight
    real(dp), intent(in) :: weight(:)

    !> Number of spin-paired microstates
    integer, intent(in) :: Lpaired

    !> current index in loop L
    integer, intent(in) :: iL

    if (iL <= Lpaired) then
      Fc(:,:) = Fc + 0.5_dp * hamSqr * (weight(iL) + weight(iL))
    else
      if (mod(iL,2) == 1) then
        Fc(:,:) = Fc + 0.5_dp * hamSqr * (weight(iL) + weight(iL+1))
      else
        Fc(:,:) = Fc + 0.5_dp * hamSqr * (weight(iL) + weight(iL-1))
      end if
    end if

  end subroutine fockFcAO_


  !> Calculate fock matrix for active orbitals in AO basis
  subroutine fockFaAO_(hamSqr, weight, fillingL, orbFON, Nc, Na, &
      & Lpaired, iL, Fa)

    !> dense fock matrix for active orbitals
    real(dp), intent(inout) :: Fa(:,:,:)

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
      & fillingL, Nc, Na, Lpaired, isHybridXc, Wab)

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
    real(dp), allocatable, intent(inout) :: hamSqrL(:,:,:,:)

    !> Sparse Hamiltonian matrix for each microstate
    real(dp), allocatable, intent(in) :: hamSpL(:,:,:)

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
    logical, intent(in) :: isHybridXc

    !> converged Lagrangian values within active space
    real(dp), intent(out) :: Wab(:,:)

    real(dp), allocatable :: tmpHam(:,:)
    real(dp), allocatable :: tmpHamL(:,:,:)

    integer :: nLocalRows, nLocalCols, iL, Lmax
    integer :: ia, ib, ist, nActPair, aa, bb

    nLocalRows = size(eigenvecs,dim=1)
    nLocalCols = size(eigenvecs,dim=2)
    Lmax = size(fillingL,dim=3)
    nActPair = Na * (Na - 1) / 2

    if (.not. isHybridXc) then
      allocate(tmpHam(nLocalRows,nLocalCols))
    end if
    allocate(tmpHamL(nActPair,1,Lmax))

    tmpHamL(:,:,:) = 0.0_dp
    do ist = 1, nActPair

      call getTwoIndices(Na, ist, ia, ib, 1)
      aa = Nc + ia
      bb = Nc + ib

      do iL = 1, Lmax

        if (isHybridXc) then
          ! convert hamSqrL from AO basis to MO basis
          ! save F_{L,ab}^{\sigma} in MO basis
          ! hamSqrL has (my_ud) component
        #:if WITH_SCALAPACK
          if (ist == 1) then
            call matAO2MOBlacs(hamSqrL(:,:,1,iL), denseDesc%blacsOrbSqr, eigenvecs)
          end if
          call findActiveElements_(env, denseDesc, hamSqrL(:,:,1,iL), aa, bb, tmpHamL(ist,1,iL))
        #:else
          if (ist == 1) then
            call matAO2MO(hamSqrL(:,:,1,iL), eigenvecs)
          end if
          tmpHamL(ist,1,iL) = hamSqrL(aa,bb,1,iL)
        #:endif
        else
          tmpHam(:,:) = 0.0_dp
          ! convert from sparse to dense for hamSpL in AO basis
          ! hamSpL has (my_ud) component
          call env%globalTimer%startTimer(globalTimers%sparseToDense)
        #:if WITH_SCALAPACK
          call unpackHSRealBlacs(env%blacs, hamSpL(:,1,iL), neighbourList%iNeighbour,&
              & nNeighbourSK, iSparseStart, img2CentCell, denseDesc, tmpHam)
        #:else
          call unpackHS(tmpHam, hamSpL(:,1,iL), neighbourList%iNeighbour,&
              & nNeighbourSK, denseDesc%iAtomStart, iSparseStart, img2CentCell)
          call adjointLowerTriangle(tmpHam)
        #:endif
          call env%globalTimer%stopTimer(globalTimers%sparseToDense)
          ! convert tmpHam from AO basis to MO basis
          ! save F_{L,ab}^{\sigma} in MO basis
        #:if WITH_SCALAPACK
          call matAO2MOBlacs(tmpHam, denseDesc%blacsOrbSqr, eigenvecs)
          call findActiveElements_(env, denseDesc, tmpHam, aa, bb, tmpHamL(ist,1,iL))
        #:else
          call matAO2MO(tmpHam, eigenvecs)
          tmpHamL(ist,1,iL) = tmpHam(aa,bb)
        #:endif
        end if

      end do
    #:if WITH_SCALAPACK
      call mpifx_allreduceip(env%mpi%globalComm, tmpHamL(ist,1,:), MPI_SUM)
    #:endif

      ! calculate the Lagrangian eps_{ab} and state-interaction term
      Wab(ist,1) = 0.0_dp
      Wab(ist,2) = 0.0_dp
      do iL = 1, Lmax
        if (iL <= Lpaired) then
          Wab(ist,1) = Wab(ist,1) + fillingL(aa,1,iL) * &
              & tmpHamL(ist,1,iL) * ( weight(iL) + weight(iL) )
          Wab(ist,2) = Wab(ist,2) + fillingL(bb,1,iL) * &
              & tmpHamL(ist,1,iL) * ( weight(iL) + weight(iL) )
        else
          if (mod(iL,2) == 1) then
            Wab(ist,1) = Wab(ist,1) + fillingL(aa,1,iL) * &
                & tmpHamL(ist,1,iL) * ( weight(iL) + weight(iL+1) )
            Wab(ist,2) = Wab(ist,2) + fillingL(bb,1,iL) * &
                & tmpHamL(ist,1,iL) * ( weight(iL) + weight(iL+1) )
          else
            Wab(ist,1) = Wab(ist,1) + fillingL(aa,1,iL) * &
                & tmpHamL(ist,1,iL) * ( weight(iL) + weight(iL-1) )
            Wab(ist,2) = Wab(ist,2) + fillingL(bb,1,iL) * &
                & tmpHamL(ist,1,iL) * ( weight(iL) + weight(iL-1) )
          end if
        end if
      end do

    end do

  end subroutine getLagrangians_


#:if WITH_SCALAPACK
  !> Find elements of active space to calculate Lagrangian values
  subroutine findActiveElements_(env, denseDesc, hamSqr, aa, bb, hamActive)

    !> Environment settings
    type(TEnvironment), intent(inout) :: env

    !> Dense matrix descriptor
    type(TDenseDescr), intent(in) :: denseDesc

    !> Dense Hamiltonian matrix
    real(dp), intent(in) :: hamSqr(:,:)

    !> Index for first active orbital
    integer, intent(in) :: aa

    !> Index for second active orbital
    integer, intent(in) :: bb

    !> Hamiltonian matrix element for active space
    real(dp), intent(inout) :: hamActive

    integer :: ind_a, ind_b, iOrb, iGlob, jOrb, jGlob

    ind_a = 0
    ind_b = 0
    do jOrb = 1, size(hamSqr, dim=2)
      jGlob = scalafx_indxl2g(jOrb, denseDesc%blacsOrbSqr(NB_), env%blacs%orbitalGrid%mycol,&
          & denseDesc%blacsOrbSqr(CSRC_), env%blacs%orbitalGrid%ncol)
      do iOrb = 1, size(hamSqr, dim=1)
        iGlob = scalafx_indxl2g(iOrb, denseDesc%blacsOrbSqr(MB_), env%blacs%orbitalGrid%myrow,&
            & denseDesc%blacsOrbSqr(RSRC_), env%blacs%orbitalGrid%nrow)

        if (iGlob == aa .and. jGlob == bb) then
          ind_a = iOrb
          ind_b = jOrb
          exit
        end if

      end do
    end do
    if (ind_a /= 0 .and. ind_b /= 0) hamActive = hamSqr(ind_a,ind_b)

  end subroutine findActiveElements_
#:endif


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


end module dftbp_reks_reksen
