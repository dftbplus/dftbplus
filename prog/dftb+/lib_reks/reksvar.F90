!--------------------------------------------------------------------------------------------------!
!  DFTB+: general package for performing fast atomistic simulations                                !
!  Copyright (C) 2006 - 2020  DFTB+ developers group                                               !
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
module dftbp_reksvar

  use dftbp_accuracy
  use dftbp_message
  use dftbp_orbitals

  implicit none

  private

  public :: TReksInp, TReksCalc, REKS_init, reksTypes

  type :: TReksTypesEnum

    !> Do not use REKS calculations
    integer :: noReks = 0

    !> SSR(2,2) calculations
    integer :: ssr22 = 1

    !> SSR(4,4) calculations
    integer :: ssr44 = 2

  end type TReksTypesEnum

  !> Container for enumerated REKS types
  type(TReksTypesEnum), parameter :: reksTypes = TReksTypesEnum()

  !> Data type for initial values for REKS calculations
  type :: TReksInp

    !> Type of REKS calculations
    integer :: reksAlg

    !> REKS: energy variables

    !> Minimized energy functional
    !> In SSR(2,2), 1: PPS, 2: (PPS+OSS)/2
    integer :: Efunction

    !> Decide the energy states in SA-REKS
    !> If true, it includes all possible states in current active space
    !> If false, it includes the states used in minimized energy functional
    logical :: tAllStates

    !> Calculate SSR state with inclusion of SI, otherwise calculate SA-REKS state
    logical :: tSSR


    !> Target SSR state
    integer :: rstate

    !> Target microstate
    integer :: Lstate

    !> Read initial guess for eigenvectors in REKS
    !> If true, initial eigenvectors are obtained from 'eigenvec.bin'
    !> If false, initial eigenvectors are obtained from diagonalization of H0
    logical :: tReadMO

    !> Maximum iteration used in FON optimization
    integer :: FonMaxIter

    !> Shift value in SCC cycle
    real(dp) :: shift

    !> Scaling constant for atomic spin constants
    real(dp), allocatable :: Tuning(:)

    !> Calculate transition dipole moments
    logical :: tTDP


    !> REKS: gradient variables

    !> Algorithms to calculate analytic gradients
    !> 1: preconditioned conjugate gradient (PCG)
    !> 2: conjugate gradient (CG)
    !> 3: direct inverse-matrix multiplication
    integer :: Glevel

    !> Maximum iteration used in calculation of gradient with PCG and CG
    integer :: CGmaxIter

    !> Tolerance used in calculation of gradient with PCG and CG
    real(dp) :: Glimit

    !> Use preconditioner for conjugate gradient algorithm
    logical :: tPrecond

    !> Save 'A' and 'Hxc' to memory in gradient calculation
    !> Usually, it shows fast speed but requires large memory ~O(N^4) orders
    logical :: tSaveMem


    !> Calculate relaxed density of SSR or SA-REKS state
    logical :: tRD

    !> Calculate nonadiabatic coupling vectors
    logical :: tNAC


    !> REKS: system variables

    !> Print level in standard output file
    !> 0: energy information
    !> 1: 0 + gradient information
    !> 2: 1 + detailed energy contribution of each microstate
    integer :: Plevel

  end type TReksInp

  !> Data type for REKS internal settings
  type :: TReksCalc

    !> Type of REKS calculations
    integer :: reksAlg

    !> REKS: energy variables (input)

    !> Minimized energy functional
    integer :: Efunction

    !> Decide the energy states in SA-REKS
    logical :: tAllStates

    !> Calculate SSR state with inclusion of SI, otherwise calculate SA-REKS state
    logical :: tSSR


    !> Target SSR state
    integer :: rstate

    !> Target microstate
    integer :: Lstate

    !> Read initial guess for eigenvectors in REKS
    logical :: tReadMO

    !> Maximum iteration used in FON optimization
    integer :: FonMaxIter

    !> Shift value in SCC cycle
    real(dp) :: shift

    !> Scaling constant for atomic spin constants
    real(dp), allocatable :: Tuning(:)

    !> Calculate transition dipole moments
    logical :: tTDP


    !> REKS: gradient variables (input)

    !> Algorithms to calculate analytic gradients
    integer :: Glevel

    !> Maximum iteration used in calculation of gradient with PCG and CG
    integer :: CGmaxIter

    !> Tolerance used in calculation of gradient with PCG and CG
    real(dp) :: Glimit

    !> Save 'A' and 'Hxc' to memory in gradient calculation
    logical :: tSaveMem


    !> Calculate relaxed density of SSR or SA-REKS state
    logical :: tRD

    !> Calculate nonadiabatic coupling vectors
    logical :: tNAC


    !> REKS: system variables (input)

    !> Print level in standard output file
    integer :: Plevel


    !> REKS: energy variables

    !> Number of microstates
    integer :: Lmax

    !> Number of independent R matrix used in gradient calculations
    integer :: LmaxR

    !> Number of spin-paired microstates
    integer :: Lpaired

    !> Number of core orbitals
    integer :: Nc

    !> Number of active orbitals
    integer :: Na

    !> Number of states
    integer :: nstates

    !> Number of states used in state-averaging
    integer :: SAstates

    !> Smoothing factor used in FON optimization
    real(dp) :: delta = 0.4_dp

    !> Hessian of FONs
    real(dp) :: hess

    !> Third order DFTB
    logical :: t3rd

    !> Whether to run a range separated calculation
    logical :: isRangeSep

    !> Do we need forces?
    logical :: tForces

    !> if calculation is periodic
    logical :: tPeriodic

    !> Can stress be calculated?
    logical :: tStress

    !> If external charges must be considered
    logical :: tExtChrg

    !> If charges should be blured
    logical :: tBlur


    !> get atom index from AO index
    integer, allocatable :: getAtomIndex(:)

    !> get dense AO index from sparse AO array
    integer, allocatable :: getDenseAO(:,:)

    !> get dense atom index from sparse atom array
    integer, allocatable :: getDenseAtom(:,:)


    !> Dense overlap matrix
    real(dp), allocatable :: overSqr(:,:)

    !> Filling for each microstate
    real(dp), allocatable :: fillingL(:,:,:)

    !> Dense density matrix for each microstate
    real(dp), allocatable :: rhoSqrL(:,:,:,:)

    !> Sparse density matrix for each microstate
    real(dp), allocatable :: rhoSpL(:,:,:)

    !> Dense delta density matrix for each microstate
    real(dp), allocatable :: deltaRhoSqrL(:,:,:,:)

    !> Mulliken population for each microstate
    real(dp), allocatable :: qOutputL(:,:,:,:)

    !> charge per atomic shell for each microstate
    real(dp), allocatable :: chargePerShellL(:,:,:,:)


    !> internal shell and spin resolved potential for each microstate
    real(dp), allocatable :: intShellL(:,:,:,:)

    !> internal block and spin resolved potential for each microstate
    real(dp), allocatable :: intBlockL(:,:,:,:,:)


    !> Dense Hamiltonian matrix for each microstate
    real(dp), allocatable :: hamSqrL(:,:,:,:)

    !> Sparse Hamiltonian matrix for each microstate
    real(dp), allocatable :: hamSpL(:,:,:)


    !> Weight for each microstate per state
    real(dp), allocatable :: weightL(:,:)

    !> Weights used in state-averaging
    real(dp), allocatable :: SAweight(:)

    !> Weight of each microstate for state to be optimized; weight = weightL * SAweight
    real(dp), allocatable :: weight(:)


    !> Fractional occupation numbers of active orbitals
    real(dp), allocatable :: FONs(:,:)

    !> energy of states
    real(dp), allocatable :: energy(:)


    !> non SCC energy for each microstate
    real(dp), allocatable :: enLnonSCC(:)

    !> SCC energy for each microstate
    real(dp), allocatable :: enLscc(:)

    !> spin-polarized energy for each microstate
    real(dp), allocatable :: enLspin(:)

    !> 3rd order SCC energy for each microstate
    real(dp), allocatable :: enL3rd(:)

    !> Long-range corrected energy for each microstate
    real(dp), allocatable :: enLfock(:)

    !> total energy for each microstate
    real(dp), allocatable :: enLtot(:)


    !> dense fock matrix for core orbitals
    real(dp), allocatable :: fockFc(:,:)

    !> dense fock matrix for active orbitals
    real(dp), allocatable :: fockFa(:,:,:)

    !> dense pseudo-fock matrix
    real(dp), allocatable :: fock(:,:)


    !> eigenvectors from pesudo-fock matrix
    real(dp), allocatable :: eigvecsFock(:,:)

    !> eigenvectors from SA-REKS state
    real(dp), allocatable :: eigvecsSSR(:,:)


    !> REKS: gradient variables

    !> constant calculated from hessian and energy of microstates
    real(dp) :: G1

    !> Ordering between RmatL and fillingL
    integer, allocatable :: orderRmatL(:)

    !> Dense non-scc Hamiltonian derivative in AO basis
    real(dp), allocatable :: Hderiv(:,:,:)

    !> Dense overlap derivative in AO basis
    real(dp), allocatable :: Sderiv(:,:,:)

    !> sparse energy-weighted density matrix for each microstate
    real(dp), allocatable :: edmSpL(:,:)

    !> gradients for each microstate except orbital derivative terms
    real(dp), allocatable :: gradL(:,:,:)


    !> scc gamma integrals in AO basis
    real(dp), allocatable :: GammaAO(:,:)

    !> scc gamma derivative integrals
    real(dp), allocatable :: GammaDeriv(:,:,:)

    !> spin W in AO basis
    real(dp), allocatable :: SpinAO(:,:)

    !> long-range gamma integrals in AO basis
    real(dp), allocatable :: LrGammaAO(:,:)

    !> long-range gamma derivative integrals
    real(dp), allocatable :: LrGammaDeriv(:,:,:)


    !> Hartree-XC kernel with sparse form with same spin part
    real(dp), allocatable :: HxcSpS(:,:)

    !> Hartree-XC kernel with sparse form with different spin part
    real(dp), allocatable :: HxcSpD(:,:)

    !> Hartree-XC kernel with half dense form with same spin part
    real(dp), allocatable :: HxcHalfS(:,:)

    !> Hartree-XC kernel with half dense form with different spin part
    real(dp), allocatable :: HxcHalfD(:,:)

    !> Hartree-XC kernel with dense form with same spin part
    real(dp), allocatable :: HxcSqrS(:,:,:,:)

    !> Hartree-XC kernel with dense form with different spin part
    real(dp), allocatable :: HxcSqrD(:,:,:,:)


    !> modified weight of each microstate
    real(dp), allocatable :: weightIL(:)

    !> anti-symmetric matrices originated from Hamiltonians
    real(dp), allocatable :: omega(:)

    !> state-interaction term used in SSR gradients
    real(dp), allocatable :: Rab(:,:)

    !> super A hessian matrix in front of orbital derivatives
    real(dp), allocatable :: Aall(:,:)

    !> super A hessian matrix with one-electron term in front of orbital derivatives
    real(dp), allocatable :: A1e(:,:)

    !> preconditioner of super A hessian matrix with one-electron term in front of orbital
    !> derivatives
    real(dp), allocatable :: A1ePre(:,:)


    !> SA-REKS state vector
    real(dp), allocatable :: XT(:,:)

    !> state-interaction term vector
    real(dp), allocatable :: XTdel(:,:)

    !> auxiliary matrix in AO basis related to state-interaction term
    real(dp), allocatable :: RdelL(:,:,:,:)

    !> auxiliary matrix in AO basis related to state-interaction term
    real(dp), allocatable :: ZdelL(:,:,:)

    !> auxiliary matrix in AO basis related to state-interaction term
    real(dp), allocatable :: Q1del(:,:,:)

    !> auxiliary matrix in MO basis related to state-interaction term
    real(dp), allocatable :: Q2del(:,:,:)


    !> solution of A * Z = X equation with X is XT
    real(dp), allocatable :: ZT(:,:)

    !> auxiliary matrix in AO basis related to SA-REKS term
    real(dp), allocatable :: RmatL(:,:,:,:)

    !> auxiliary matrix in AO basis related to SA-REKS term
    real(dp), allocatable :: ZmatL(:,:,:)

    !> auxiliary matrix in MO basis related to SA-REKS term
    real(dp), allocatable :: Q1mat(:,:)

    !> auxiliary matrix in MO basis related to SA-REKS term
    real(dp), allocatable :: Q2mat(:,:)


    !> solution of A * Z = X equation with X is XTdel
    real(dp), allocatable :: ZTdel(:,:)

    !> auxiliary matrix in AO basis related to state-interaction term
    real(dp), allocatable :: tmpRL(:,:,:,:)


    !> gradient of SSR state
    real(dp), allocatable :: SSRgrad(:,:,:)

    !> gradient of SA-REKS state
    real(dp), allocatable :: SAgrad(:,:,:)

    !> gradient of state-interaction term
    real(dp), allocatable :: SIgrad(:,:,:)

    !> gradient of averaged state
    real(dp), allocatable :: avgGrad(:,:)


    !> difference gradient vector, G
    real(dp), allocatable :: nacG(:,:,:)

    !> nonadiabatic coupling vector, H
    real(dp), allocatable :: nacH(:,:,:)


    !> REKS: relaxed density & transition dipole variables

    !> unrelaxed density matrix for target SSR or SA-REKS state
    real(dp), allocatable :: unrelRhoSqr(:,:)

    !> unrelaxed transition density matrix between SSR or SA-REKS states
    real(dp), allocatable :: unrelTdm(:,:,:)

    !> relaxed density matrix for target SSR or SA-REKS state
    real(dp), allocatable :: relRhoSqr(:,:)

    !> transition dipole moment between states
    real(dp), allocatable :: tdp(:,:)


    !> REKS: point charges (QM/MM) variables

    !> coordinates and charges of external point charges
    real(dp), allocatable :: extCharges(:,:)

    !> Width of the Gaussians if the charges are blurred
    real(dp), allocatable :: blurWidths(:)


    !> REKS: periodic variables

    !> parameter for Ewald
    real(dp) :: alpha

    !> parameter for cell volume
    real(dp) :: volume

    !> real lattice points for Ewald-sum
    real(dp), allocatable :: rVec(:,:)

    !> lattice points for reciprocal Ewald
    real(dp), allocatable :: gVec(:,:)


    !> REKS: stress variables

    !> electronic stress part for lattice optimization
    real(dp), allocatable :: elecStressL(:,:,:)

  contains

    !> Reallocate sparse arrays used in REKS
    procedure :: reallocate

  end type TReksCalc

  contains

  !> Initialize REKS data from REKS input
  subroutine REKS_init(this, inp, orb, spinW, nSpin, nEl, nChrgs, extChrg, &
      & blurWidths, t3rd, isRangeSep, tForces, tPeriodic, tStress, tDipole)
    
    !> data type for REKS
    type(TReksCalc), intent(out) :: this

    !> data type for REKS input
    type(TReksInp), intent(inout) :: inp

    !> Atomic orbital information
    type(TOrbitals), intent(in) :: orb

    !> Spin W values
    real(dp), intent(inout) :: spinW(:,:,:)

    !> Number of spin components, 1 is unpolarised, 2 is polarised, 4 is noncolinear / spin-orbit
    integer, intent(inout) :: nSpin

    !> nr. of electrons
    real(dp), intent(in) :: nEl

    !> Nr. of external charges
    integer, intent(in) :: nChrgs

    !> coordinates and charges of external point charges
    real(dp), allocatable, intent(in) :: extChrg(:,:)

    !> Width of the Gaussians if the charges are blurred
    real(dp), allocatable, intent(in) :: blurWidths(:)

    !> Third order DFTB
    logical, intent(in) :: t3rd

    !> Whether to run a range separated calculation
    logical, intent(in) :: isRangeSep

    !> Do we need forces?
    logical, intent(in) :: tForces

    !> if calculation is periodic
    logical, intent(in) :: tPeriodic

    !> Can stress be calculated?
    logical, intent(in) :: tStress

    !> calculate an electric dipole?
    logical, intent(inout) :: tDipole

    integer :: nOrb, mOrb, mShell, nOrbHalf
    integer :: nstates, SAstates, nstHalf
    integer :: Nc, Na, Nv, superN, LmaxR, Lmax
    integer :: iAt, nAtom, nType

    ! Set REKS input variables

    this%reksAlg = inp%reksAlg

    this%Efunction = inp%Efunction
    this%tAllStates = inp%tAllStates
    this%tSSR = inp%tSSR

    this%rstate = inp%rstate
    this%Lstate = inp%Lstate
    this%tReadMO = inp%tReadMO
    this%FonMaxIter = inp%FonMaxIter
    this%shift = inp%shift

    this%tTDP = inp%tTDP

    this%Glevel = inp%Glevel
    if (this%Glevel == 1 .or. this%Glevel == 2) then
      this%CGmaxIter = inp%CGmaxIter
      this%Glimit = inp%Glimit
      this%tSaveMem = inp%tSaveMem
    end if

    this%Plevel = inp%Plevel

    this%tRD = inp%tRD
    this%tNAC = inp%tNAC

    ! Set REKS variables

    select case (this%reksAlg)
    case (reksTypes%noReks)
    case (reksTypes%ssr22)
      call setSSR22conditions(this, nEl)
    case (reksTypes%ssr44)
      call error("SSR(4,4) is not implemented yet")
    end select

    nOrb = orb%nOrb
    mOrb = orb%mOrb
    mShell = orb%mShell
    nstates = this%nstates
    SAstates = this%SAstates

    ! Here, nSpin changes to two for calculation of microstates
    nSpin = 2
    Lmax = this%Lmax
    LmaxR = this%LmaxR
    Nc = this%Nc
    Na = this%Na
    Nv = nOrb - Nc - Na
    nOrbHalf = nOrb * (nOrb + 1) / 2
    nstHalf = nstates * (nstates - 1) / 2
    superN = Nc*Nv + Na*(Nc+Nv) + Na*(Na-1)/2

    nAtom = size(orb%nOrbAtom,dim=1)
    nType = size(inp%Tuning,dim=1)

    this%t3rd = t3rd
    this%isRangeSep = isRangeSep

    this%tForces = tForces
    this%tPeriodic = tPeriodic
    this%tStress = tStress
    this%tExtChrg = (nChrgs > 0)
    ! The standard 1.0e-7 is given in setExternalCharges routine of externalcharges.F90
    if (allocated(blurWidths)) then
      this%tBlur = any(blurWidths > 1.0e-7_dp)
    else
      this%tBlur = .false.
    end if
    ! Set tDipole is true when we compute the relaxed density for REKS
    tDipole = (tDipole .and. this%tRD)

    ! Check REKS requirements

    call checkReksRequirements(this)

    select case (this%reksAlg)
    case (reksTypes%noReks)
    case (reksTypes%ssr22)
      call checkSSR22Requirements(this)
    case (reksTypes%ssr44)
      call error("SSR(4,4) is not implemented yet")
    end select

    ! Allocate REKS variables

    ! REKS: energy variables

    allocate(this%getAtomIndex(nOrb))
    allocate(this%getDenseAO(0,2))
    allocate(this%getDenseAtom(0,2))

    allocate(this%overSqr(nOrb,nOrb))
    allocate(this%fillingL(nOrb,nSpin,Lmax))

    if (this%tForces) then
      allocate(this%rhoSqrL(nOrb,nOrb,1,Lmax))
    else
      allocate(this%rhoSpL(0,1,Lmax))
    end if

    if (this%isRangeSep) then
      allocate(this%deltaRhoSqrL(nOrb,nOrb,1,Lmax))
    end if

    allocate(this%qOutputL(mOrb,nAtom,nSpin,Lmax))
    allocate(this%chargePerShellL(mShell,nAtom,nSpin,Lmax))

    allocate(this%intShellL(mShell,nAtom,nSpin,Lmax))
    allocate(this%intBlockL(mOrb,mOrb,nAtom,nSpin,Lmax))

    if (this%isRangeSep) then
      allocate(this%hamSqrL(nOrb,nOrb,1,Lmax))
    else
      allocate(this%hamSpL(0,1,Lmax))
    end if

    allocate(this%weightL(nstates,Lmax))
    allocate(this%SAweight(SAstates))
    allocate(this%weight(Lmax))

    allocate(this%energy(nstates))

    allocate(this%enLnonSCC(Lmax))
    allocate(this%enLSCC(Lmax))
    allocate(this%enLspin(Lmax))

    if (this%t3rd) then
      allocate(this%enL3rd(Lmax))
    end if

    if (this%isRangeSep) then
      allocate(this%enLfock(Lmax))
    end if

    allocate(this%enLtot(Lmax))

    allocate(this%fockFc(nOrb,nOrb))
    allocate(this%fockFa(nOrb,nOrb,Na))
    allocate(this%fock(nOrb,nOrb))

    allocate(this%eigvecsFock(nOrb,nOrb))
    allocate(this%eigvecsSSR(nstates,nstates))

    ! REKS: gradient variables

    if (this%tForces) then

      allocate(this%Hderiv(nOrb,nOrb,3))
      allocate(this%Sderiv(nOrb,nOrb,3))
      allocate(this%edmSpL(0,Lmax))
      allocate(this%gradL(3,nAtom,Lmax))

      if (this%Efunction > 1) then

        allocate(this%orderRmatL(Lmax))

        allocate(this%GammaAO(nOrb,nOrb))
        allocate(this%GammaDeriv(nAtom,nAtom,3))
        allocate(this%SpinAO(nOrb,nOrb))

        if (this%isRangeSep) then
          allocate(this%LrGammaAO(nOrb,nOrb))
          allocate(this%LrGammaDeriv(nAtom,nAtom,3))
        end if

        allocate(this%weightIL(Lmax))
        allocate(this%omega(superN))
        if (this%tSSR) then
          allocate(this%Rab(nstates,nstates))
        end if

        if (this%Glevel == 1 .or. this%Glevel == 2) then
          if (this%tSaveMem) then
            if (this%isRangeSep) then
              allocate(this%HxcHalfS(nOrbHalf,nOrbHalf))
              allocate(this%HxcHalfD(nOrbHalf,nOrbHalf))
            else
              allocate(this%HxcSpS(0,0))
              allocate(this%HxcSpD(0,0))
            end if
            allocate(this%A1e(superN,superN))
            allocate(this%A1ePre(superN,superN))
          end if
        else if (this%Glevel == 3) then
          allocate(this%HxcSqrS(nOrb,nOrb,nOrb,nOrb))
          allocate(this%HxcSqrD(nOrb,nOrb,nOrb,nOrb))
          allocate(this%Aall(superN,superN))
        end if

        if (this%tSSR) then
          allocate(this%XT(superN,nstates))
          allocate(this%XTdel(superN,nstHalf))
          allocate(this%RdelL(nOrb,nOrb,LmaxR,nstHalf))
          allocate(this%ZdelL(nOrb,nOrb,Lmax))
          allocate(this%Q1del(nOrb,nOrb,nstHalf))
          allocate(this%Q2del(nOrb,nOrb,nstHalf))
        else
          allocate(this%XT(superN,1))
        end if

        if (this%tNAC) then
          allocate(this%ZT(superN,nstates))
          allocate(this%RmatL(nOrb,nOrb,LmaxR,nstates))
        else
          allocate(this%ZT(superN,1))
          allocate(this%RmatL(nOrb,nOrb,LmaxR,1))
        end if
        allocate(this%ZmatL(nOrb,nOrb,Lmax))
        allocate(this%Q1mat(nOrb,nOrb))
        allocate(this%Q2mat(nOrb,nOrb))

        if (this%tNAC) then
          allocate(this%ZTdel(superN,nstHalf))
          allocate(this%tmpRL(nOrb,nOrb,LmaxR,nstHalf))
        end if

        if (this%tNAC) then
          allocate(this%SSRgrad(3,nAtom,nstates))
          allocate(this%SAgrad(3,nAtom,nstates))
          allocate(this%SIgrad(3,nAtom,nstHalf))
          allocate(this%avgGrad(3,nAtom))
        else
          allocate(this%SSRgrad(3,nAtom,1))
        end if

        if (this%tNAC) then
          allocate(this%nacG(3,nAtom,nstHalf))
          allocate(this%nacH(3,nAtom,nstHalf))
        end if

      end if

    end if


    ! REKS: relaxed density & transition dipole variables

    allocate(this%unrelRhoSqr(nOrb,nOrb))

    if (this%tTDP .and. this%Lstate == 0) then
      allocate(this%unrelTdm(nOrb,nOrb,nstHalf))
    end if

    if (this%tRD) then
      allocate(this%relRhoSqr(nOrb,nOrb))
    end if

    if (this%tTDP .and. this%Lstate == 0) then
      allocate(this%tdp(3,nstHalf))
    end if

    if (this%tForces) then

      ! REKS: point charges (QM/MM) variables

      if (this%tExtChrg) then
        allocate(this%extCharges(4,nChrgs))
        if (this%tBlur) then
          allocate(this%blurWidths(nChrgs))
        end if
      end if

      ! REKS: stress variables

      if (this%tStress) then
        allocate(this%elecStressL(3,3,Lmax))
      end if

    end if


    ! REKS: energy variables

    this%getAtomIndex(:) = 0
    this%getDenseAO(:,:) = 0
    this%getDenseAtom(:,:) = 0

    this%overSqr(:,:) = 0.0_dp
    this%fillingL(:,:,:) = 0.0_dp

    if (this%tForces) then
      this%rhoSqrL(:,:,:,:) = 0.0_dp
    else
      this%rhoSpL(:,:,:) = 0.0_dp
    end if

    if (this%isRangeSep) then
      this%deltaRhoSqrL(:,:,:,:) = 0.0_dp
    end if

    this%qOutputL(:,:,:,:) = 0.0_dp
    this%chargePerShellL(:,:,:,:) = 0.0_dp

    this%intShellL(:,:,:,:) = 0.0_dp
    this%intBlockL(:,:,:,:,:) = 0.0_dp

    if (this%isRangeSep) then
      this%hamSqrL(:,:,:,:) = 0.0_dp
    else
      this%hamSpL(:,:,:) = 0.0_dp
    end if

    this%weightL(:,:) = 0.0_dp
    this%weight(:) = 0.0_dp

    this%FONs(:,:) = 0.0_dp
    this%energy(:) = 0.0_dp

    this%enLnonSCC(:) = 0.0_dp
    this%enLSCC(:) = 0.0_dp
    this%enLspin(:) = 0.0_dp

    if (this%t3rd) then
      this%enL3rd(:) = 0.0_dp
    end if

    if (this%isRangeSep) then
      this%enLfock(:) = 0.0_dp
    end if

    this%enLtot(:) = 0.0_dp

    this%fockFc(:,:) = 0.0_dp
    this%fockFa(:,:,:) = 0.0_dp
    this%fock(:,:) = 0.0_dp

    this%eigvecsFock(:,:) = 0.0_dp
    this%eigvecsSSR(:,:) = 0.0_dp

    ! REKS: gradient variables

    if (this%tForces) then

      this%Hderiv(:,:,:) = 0.0_dp
      this%Sderiv(:,:,:) = 0.0_dp
      this%edmSpL(:,:) = 0.0_dp
      this%gradL(:,:,:) = 0.0_dp

      if (this%Efunction > 1) then

        this%GammaAO(:,:) = 0.0_dp
        this%GammaDeriv(:,:,:) = 0.0_dp
        this%SpinAO(:,:) = 0.0_dp

        if (this%isRangeSep) then
          this%LrGammaAO(:,:) = 0.0_dp
          this%LrGammaDeriv(:,:,:) = 0.0_dp
        end if

        this%weightIL(:) = 0.0_dp
        this%omega(:) = 0.0_dp
        if (this%tSSR) then
          this%Rab(:,:) = 0.0_dp
        end if

        if (this%Glevel == 1 .or. this%Glevel == 2) then
          if (this%tSaveMem) then
            if (this%isRangeSep) then
              this%HxcHalfS(:,:) = 0.0_dp
              this%HxcHalfD(:,:) = 0.0_dp
            else
              this%HxcSpS(:,:) = 0.0_dp
              this%HxcSpD(:,:) = 0.0_dp
            end if
            this%A1e(:,:) = 0.0_dp
            this%A1ePre(:,:) = 0.0_dp
          end if
        else if (this%Glevel == 3) then
          this%HxcSqrS(:,:,:,:) = 0.0_dp
          this%HxcSqrD(:,:,:,:) = 0.0_dp
          this%Aall(:,:) = 0.0_dp
        end if

        this%XT = 0.0_dp
        if (this%tSSR) then
          this%XTdel(:,:) = 0.0_dp
          this%RdelL(:,:,:,:) = 0.0_dp
          this%ZdelL(:,:,:) = 0.0_dp
          this%Q1del(:,:,:) = 0.0_dp
          this%Q2del(:,:,:) = 0.0_dp
        end if

        this%ZT(:,:) = 0.0_dp
        this%RmatL(:,:,:,:) = 0.0_dp
        this%ZmatL(:,:,:) = 0.0_dp
        this%Q1mat(:,:) = 0.0_dp
        this%Q2mat(:,:) = 0.0_dp

        if (this%tNAC) then
          this%ZTdel(:,:) = 0.0_dp
          this%tmpRL(:,:,:,:) = 0.0_dp
        end if

        this%SSRgrad(:,:,:) = 0.0_dp
        if (this%tNAC) then
          this%SAgrad(:,:,:) = 0.0_dp
          this%SIgrad(:,:,:) = 0.0_dp
          this%avgGrad(:,:) = 0.0_dp
        end if

        if (this%tNAC) then
          this%nacG(:,:,:) = 0.0_dp
          this%nacH(:,:,:) = 0.0_dp
        end if

      end if

    end if


    ! REKS: relaxed density & transition dipole variables

    this%unrelRhoSqr(:,:) = 0.0_dp

    if (this%tTDP .and. this%Lstate == 0) then
      this%unrelTdm(:,:,:) = 0.0_dp
    end if

    if (this%tRD) then
      this%relRhoSqr(:,:) = 0.0_dp
    end if

    if (this%tTDP .and. this%Lstate == 0) then
      this%tdp(:,:) = 0.0_dp
    end if

    if (this%tForces) then

      ! REKS: point charges (QM/MM) variables

      if (this%tExtChrg) then
        this%extCharges(1:3,:) = extChrg(1:3,:)
        ! Adapt to internal sign convention for potentials (electron has positive charge)
        this%extCharges(4,:) = -extChrg(4,:)
        if (this%tBlur) then
          this%blurWidths(:) = blurWidths
        end if
      end if

      ! REKS: stress variables

      if (this%tStress) then
        this%elecStressL(:,:,:) = 0.0_dp
      end if

    end if


    ! REKS: initialize variables

    this%SAweight(:) = 1.0_dp / real(this%SAstates, dp)

    call move_alloc(inp%Tuning, this%Tuning)

    ! Scale up or down the atomic spin constants w.r.t. the systems
    ! iAt : loop for atomic species
    ! order of iAt = order in .gen file
    do iAt = 1, nType
      spinW(:,:,iAt) = this%Tuning(iAt) * spinW(:,:,iAt)
    end do

    ! Set the ordering information betwen R_mat_L and filling_L
    if (this%Efunction > 1 .and. this%tForces) then
      select case (this%reksAlg)
      case (reksTypes%noReks)
      case (reksTypes%ssr22)
        ! R_mat_L has 4 elements and filling_L has 6 elements in (2,2) case
        this%orderRmatL(:) = [1, 2, 1, 2, 3, 4]
      case (reksTypes%ssr44)
        call error("SSR(4,4) is not implemented yet")
      end select
    end if

  contains

    !> Check REKS common requirements
    subroutine setSSR22conditions(this, nEl)

      !> data type for REKS
      type(TReksCalc), intent(inout) :: this

      !> nr. of electrons
      real(dp), intent(in) :: nEl

      this%Nc = int(real(nEl, dp)) / 2 - 1
      this%Na = 2
      this%Lmax = 6
      this%LmaxR = 4
      this%Lpaired = 2
      if (this%Efunction == 1) then
        ! Only PPS state is minimized; single-state REKS
        this%SAstates = 1
        if (.not. this%tAllStates) then
          ! PPS state
          this%nstates = 1
        else
          call error("IncludeAllStates should be set to No in single-state REKS")
        end if
      else if (this%Efunction == 2) then
        ! (PPS+OSS)/2 state is minimized; SA-REKS
        this%SAstates = 2
        if (.not. this%tAllStates) then
          ! PPS, OSS state
          this%nstates = 2
        else
          ! PPS, OSS, DES state
          this%nstates = 3
        end if
      end if
      ! Fractional occupation numbers, n_a and n_b
      allocate(this%FONs(this%Na,1))

    end subroutine setSSR22conditions

    !> Check REKS common requirements
    subroutine checkReksRequirements(this)

      !> data type for REKS
      type(TReksCalc), intent(in) :: this

      ! REKS energy requirements

      if (this%tTDP) then
        if (this%Lstate > 0) then
          call error("Transition dipole is not compatible with L-th microstate")
        else if (this%Efunction == 1) then
          call error("Transition dipole is not compatible with single-state REKS")
        end if
      end if

      if (this%tSSR .and. this%Efunction == 1) then
        call error("Single-state REKS is not SSR state")
      end if

      if (this%shift < -epsilon(1.0_dp)) then
        call error("Too small shift value in REKS")
      end if

      ! REKS gradient requirements

      if (this%tForces) then

        if (this%t3rd) then
          call error("3rd order scc is not compatible with force calculation in REKS")
        end if

        if (this%Lstate > 0) then
          if (this%Efunction == 1) then
            call error("gradient of microstate is not compatible with single-state REKS")
          else if (this%tSSR) then
            call error("For gradient of microstate, please set StateInteractions = No")
          end if
        end if

        if (this%tNAC) then
          if (this%Lstate > 0) then
            call error("Nonadiabatic coupling is not compatible with gradient of microstate")
          else if (.not. this%tSSR .or. this%Efunction == 1) then
            call error("Nonadiabatic coupling is not compatible with &
                & single-state REKS or SA-REKS calculation")
          end if
        end if

      else

        if (this%tNAC) then
          call error("Nonadiabatic coupling requires force evaluation")
        else if (this%tRD) then
          call error("Relaxed density requires force evaluation")
        end if

      end if

      ! REKS stress requirements

      if (this%Efunction /= 1 .and. this%tStress) then
        ! tStress = tForces * tPeriodic * (not tExtChrg)
        call error("Only single-state REKS can evaluate stress")
      end if

      ! REKS system requirements

      if (this%Plevel > 2 .or. this%Plevel < 0) then
        call error("Wrong VerbosityLevel given, please write 0 to 2")
      end if

    end subroutine checkReksRequirements

    !> Check REKS(2,2) requirements
    subroutine checkSSR22Requirements(this)

      !> data type for REKS
      type(TReksCalc), intent(in) :: this

      ! REKS energy requirements

      if (this%rstate > 3 .or. this%rstate < 1) then
        call error("Wrong TargetState given, please write 1 to 3")
      else if (this%Lstate > 6 .or. this%Lstate < 0) then
        call error("Wrong TargetMicrostate given, please write 0 to 6")
      end if

    end subroutine checkSSR22Requirements

  end subroutine REKS_init


  !> Reallocate sparse arrays used in REKS
  subroutine reallocate(this, sparseSize)

    !> data type for REKS
    class(TReksCalc), intent(inout) :: this

    !> Total size of orbitals in the sparse data structures, where the decay of the overlap sets the
    !> sparsity pattern
    integer, intent(in) :: sparseSize

    deallocate(this%getDenseAO)
    if (.not. this%tForces) then
      deallocate(this%rhoSpL)
    end if
    if (.not. this%isRangeSep) then
      deallocate(this%hamSpL)
    end if

    if (this%tForces) then
      deallocate(this%edmSpL)
      if (this%Efunction > 1) then
        if (this%Glevel == 1 .or. this%Glevel == 2) then
          if (this%tSaveMem) then
            if (.not. this%isRangeSep) then
              deallocate(this%HxcSpS)
              deallocate(this%HxcSpD)
            end if
          end if
        end if
      end if
    end if

    allocate(this%getDenseAO(sparseSize,2))
    if (.not. this%tForces) then
      allocate(this%rhoSpL(sparseSize,1,this%Lmax))
    end if
    if (.not. this%isRangeSep) then
      allocate(this%hamSpL(sparseSize,1,this%Lmax))
    end if

    if (this%tForces) then
      allocate(this%edmSpL(sparseSize,this%Lmax))
      if (this%Efunction > 1) then
        if (this%Glevel == 1 .or. this%Glevel == 2) then
          if (this%tSaveMem) then
            if (.not. this%isRangeSep) then
              allocate(this%HxcSpS(sparseSize,sparseSize))
              allocate(this%HxcSpD(sparseSize,sparseSize))
            end if
          end if
        end if
      end if
    end if

    this%getDenseAO = 0
    if (.not. this%tForces) then
      this%rhoSpL = 0.0_dp
    end if
    if (.not. this%isRangeSep) then
      this%hamSpL = 0.0_dp
    end if

    if (this%tForces) then
      this%edmSpL = 0.0_dp
      if (this%Efunction > 1) then
        if (this%Glevel == 1 .or. this%Glevel == 2) then
          if (this%tSaveMem) then
            if (.not. this%isRangeSep) then
              this%HxcSpS = 0.0_dp
              this%HxcSpD = 0.0_dp
            end if
          end if
        end if
      end if
    end if

  end subroutine reallocate


end module dftbp_reksvar
