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


    !> internal atom and spin resolved potential
    real(dp), allocatable :: intAtom(:,:)

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
  subroutine REKS_init(self, inp, orb, spinW, nSpin, nEl, nChrgs, extChrg, &
      & blurWidths, t3rd, isRangeSep, tForces, tPeriodic, tStress, tDipole)
    
    !> data type for REKS
    type(TReksCalc), intent(out) :: self

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

    self%reksAlg = inp%reksAlg

    self%Efunction = inp%Efunction
    self%tAllStates = inp%tAllStates
    self%tSSR = inp%tSSR

    self%rstate = inp%rstate
    self%Lstate = inp%Lstate
    self%tReadMO = inp%tReadMO
    self%FonMaxIter = inp%FonMaxIter
    self%shift = inp%shift

    self%tTDP = inp%tTDP

    self%Glevel = inp%Glevel
    if (self%Glevel == 1 .or. self%Glevel == 2) then
      self%CGmaxIter = inp%CGmaxIter
      self%Glimit = inp%Glimit
      self%tSaveMem = inp%tSaveMem
    end if

    self%Plevel = inp%Plevel

    self%tRD = inp%tRD
    self%tNAC = inp%tNAC

    ! Set REKS variables

    select case (self%reksAlg)
    case (reksTypes%noReks)
    case (reksTypes%ssr22)
      call setSSR22conditions(self, nEl)
    case (reksTypes%ssr44)
      call error("SSR(4,4) is not implemented yet")
    end select

    nOrb = orb%nOrb
    mOrb = orb%mOrb
    mShell = orb%mShell
    nstates = self%nstates
    SAstates = self%SAstates

    ! Here, nSpin changes to two for calculation of microstates
    nSpin = 2
    Lmax = self%Lmax
    LmaxR = self%LmaxR
    Nc = self%Nc
    Na = self%Na
    Nv = nOrb - Nc - Na
    nOrbHalf = nOrb * (nOrb + 1) / 2
    nstHalf = nstates * (nstates - 1) / 2
    superN = Nc*Nv + Na*(Nc+Nv) + Na*(Na-1)/2

    nAtom = size(orb%nOrbAtom,dim=1)
    nType = size(inp%Tuning,dim=1)

    self%t3rd = t3rd
    self%isRangeSep = isRangeSep

    self%tForces = tForces
    self%tPeriodic = tPeriodic
    self%tStress = tStress
    self%tExtChrg = (nChrgs > 0)
    ! The standard 1.0e-7 is given in setExternalCharges routine of externalcharges.F90
    if (allocated(blurWidths)) then
      self%tBlur = any(blurWidths > 1.0e-7_dp)
    else
      self%tBlur = .false.
    end if
    ! Set tDipole is true when we compute the relaxed density for REKS
    tDipole = (tDipole .and. self%tRD)

    ! Check REKS requirements

    call checkReksRequirements(self)

    select case (self%reksAlg)
    case (reksTypes%noReks)
    case (reksTypes%ssr22)
      call checkSSR22Requirements(self)
    case (reksTypes%ssr44)
      call error("SSR(4,4) is not implemented yet")
    end select

    ! Allocate REKS variables

    ! REKS: energy variables

    allocate(self%getAtomIndex(nOrb))
    allocate(self%getDenseAO(0,2))
    allocate(self%getDenseAtom(0,2))

    allocate(self%overSqr(nOrb,nOrb))
    allocate(self%fillingL(nOrb,nSpin,Lmax))

    if (self%tForces) then
      allocate(self%rhoSqrL(nOrb,nOrb,1,Lmax))
    else
      allocate(self%rhoSpL(0,1,Lmax))
    end if

    if (self%isRangeSep) then
      allocate(self%deltaRhoSqrL(nOrb,nOrb,1,Lmax))
    end if

    allocate(self%qOutputL(mOrb,nAtom,nSpin,Lmax))
    allocate(self%chargePerShellL(mShell,nAtom,nSpin,Lmax))

    allocate(self%intAtom(nAtom,nSpin))
    allocate(self%intShellL(mShell,nAtom,nSpin,Lmax))
    allocate(self%intBlockL(mOrb,mOrb,nAtom,nSpin,Lmax))

    if (self%isRangeSep) then
      allocate(self%hamSqrL(nOrb,nOrb,1,Lmax))
    else
      allocate(self%hamSpL(0,1,Lmax))
    end if

    allocate(self%weightL(nstates,Lmax))
    allocate(self%SAweight(SAstates))
    allocate(self%weight(Lmax))

    allocate(self%energy(nstates))

    allocate(self%enLnonSCC(Lmax))
    allocate(self%enLSCC(Lmax))
    allocate(self%enLspin(Lmax))

    if (self%t3rd) then
      allocate(self%enL3rd(Lmax))
    end if

    if (self%isRangeSep) then
      allocate(self%enLfock(Lmax))
    end if

    allocate(self%enLtot(Lmax))

    allocate(self%fockFc(nOrb,nOrb))
    allocate(self%fockFa(nOrb,nOrb,Na))
    allocate(self%fock(nOrb,nOrb))

    allocate(self%eigvecsFock(nOrb,nOrb))
    allocate(self%eigvecsSSR(nstates,nstates))

    ! REKS: gradient variables

    if (self%tForces) then

      allocate(self%Hderiv(nOrb,nOrb,3))
      allocate(self%Sderiv(nOrb,nOrb,3))
      allocate(self%edmSpL(0,Lmax))
      allocate(self%gradL(3,nAtom,Lmax))

      if (self%Efunction > 1) then

        allocate(self%orderRmatL(Lmax))

        allocate(self%GammaAO(nOrb,nOrb))
        allocate(self%GammaDeriv(nAtom,nAtom,3))
        allocate(self%SpinAO(nOrb,nOrb))

        if (self%isRangeSep) then
          allocate(self%LrGammaAO(nOrb,nOrb))
          allocate(self%LrGammaDeriv(nAtom,nAtom,3))
        end if

        allocate(self%weightIL(Lmax))
        allocate(self%omega(superN))
        if (self%tSSR) then
          allocate(self%Rab(nstates,nstates))
        end if

        if (self%Glevel == 1 .or. self%Glevel == 2) then
          if (self%tSaveMem) then
            if (self%isRangeSep) then
              allocate(self%HxcHalfS(nOrbHalf,nOrbHalf))
              allocate(self%HxcHalfD(nOrbHalf,nOrbHalf))
            else
              allocate(self%HxcSpS(0,0))
              allocate(self%HxcSpD(0,0))
            end if
            allocate(self%A1e(superN,superN))
            allocate(self%A1ePre(superN,superN))
          end if
        else if (self%Glevel == 3) then
          allocate(self%HxcSqrS(nOrb,nOrb,nOrb,nOrb))
          allocate(self%HxcSqrD(nOrb,nOrb,nOrb,nOrb))
          allocate(self%Aall(superN,superN))
        end if

        if (self%tSSR) then
          allocate(self%XT(superN,nstates))
          allocate(self%XTdel(superN,nstHalf))
          allocate(self%RdelL(nOrb,nOrb,LmaxR,nstHalf))
          allocate(self%ZdelL(nOrb,nOrb,Lmax))
          allocate(self%Q1del(nOrb,nOrb,nstHalf))
          allocate(self%Q2del(nOrb,nOrb,nstHalf))
        else
          allocate(self%XT(superN,1))
        end if

        if (self%tNAC) then
          allocate(self%ZT(superN,nstates))
          allocate(self%RmatL(nOrb,nOrb,LmaxR,nstates))
        else
          allocate(self%ZT(superN,1))
          allocate(self%RmatL(nOrb,nOrb,LmaxR,1))
        end if
        allocate(self%ZmatL(nOrb,nOrb,Lmax))
        allocate(self%Q1mat(nOrb,nOrb))
        allocate(self%Q2mat(nOrb,nOrb))

        if (self%tNAC) then
          allocate(self%ZTdel(superN,nstHalf))
          allocate(self%tmpRL(nOrb,nOrb,LmaxR,nstHalf))
        end if

        if (self%tNAC) then
          allocate(self%SSRgrad(3,nAtom,nstates))
          allocate(self%SAgrad(3,nAtom,nstates))
          allocate(self%SIgrad(3,nAtom,nstHalf))
          allocate(self%avgGrad(3,nAtom))
        else
          allocate(self%SSRgrad(3,nAtom,1))
        end if

        if (self%tNAC) then
          allocate(self%nacG(3,nAtom,nstHalf))
          allocate(self%nacH(3,nAtom,nstHalf))
        end if

      end if

    end if


    ! REKS: relaxed density & transition dipole variables

    allocate(self%unrelRhoSqr(nOrb,nOrb))

    if (self%tTDP .and. self%Lstate == 0) then
      allocate(self%unrelTdm(nOrb,nOrb,nstHalf))
    end if

    if (self%tRD) then
      allocate(self%relRhoSqr(nOrb,nOrb))
    end if

    if (self%tTDP .and. self%Lstate == 0) then
      allocate(self%tdp(3,nstHalf))
    end if

    if (self%tForces) then

      ! REKS: point charges (QM/MM) variables

      if (self%tExtChrg) then
        allocate(self%extCharges(4,nChrgs))
        if (self%tBlur) then
          allocate(self%blurWidths(nChrgs))
        end if
      end if

      ! REKS: stress variables

      if (self%tStress) then
        allocate(self%elecStressL(3,3,Lmax))
      end if

    end if


    ! REKS: energy variables

    self%getAtomIndex(:) = 0
    self%getDenseAO(:,:) = 0
    self%getDenseAtom(:,:) = 0

    self%overSqr(:,:) = 0.0_dp
    self%fillingL(:,:,:) = 0.0_dp

    if (self%tForces) then
      self%rhoSqrL(:,:,:,:) = 0.0_dp
    else
      self%rhoSpL(:,:,:) = 0.0_dp
    end if

    if (self%isRangeSep) then
      self%deltaRhoSqrL(:,:,:,:) = 0.0_dp
    end if

    self%qOutputL(:,:,:,:) = 0.0_dp
    self%chargePerShellL(:,:,:,:) = 0.0_dp

    self%intAtom(:,:) = 0.0_dp
    self%intShellL(:,:,:,:) = 0.0_dp
    self%intBlockL(:,:,:,:,:) = 0.0_dp

    if (self%isRangeSep) then
      self%hamSqrL(:,:,:,:) = 0.0_dp
    else
      self%hamSpL(:,:,:) = 0.0_dp
    end if

    self%weightL(:,:) = 0.0_dp
    self%weight(:) = 0.0_dp

    self%FONs(:,:) = 0.0_dp
    self%energy(:) = 0.0_dp

    self%enLnonSCC(:) = 0.0_dp
    self%enLSCC(:) = 0.0_dp
    self%enLspin(:) = 0.0_dp

    if (self%t3rd) then
      self%enL3rd(:) = 0.0_dp
    end if

    if (self%isRangeSep) then
      self%enLfock(:) = 0.0_dp
    end if

    self%enLtot(:) = 0.0_dp

    self%fockFc(:,:) = 0.0_dp
    self%fockFa(:,:,:) = 0.0_dp
    self%fock(:,:) = 0.0_dp

    self%eigvecsFock(:,:) = 0.0_dp
    self%eigvecsSSR(:,:) = 0.0_dp

    ! REKS: gradient variables

    if (self%tForces) then

      self%Hderiv(:,:,:) = 0.0_dp
      self%Sderiv(:,:,:) = 0.0_dp
      self%edmSpL(:,:) = 0.0_dp
      self%gradL(:,:,:) = 0.0_dp

      if (self%Efunction > 1) then

        self%GammaAO(:,:) = 0.0_dp
        self%GammaDeriv(:,:,:) = 0.0_dp
        self%SpinAO(:,:) = 0.0_dp

        if (self%isRangeSep) then
          self%LrGammaAO(:,:) = 0.0_dp
          self%LrGammaDeriv(:,:,:) = 0.0_dp
        end if

        self%weightIL(:) = 0.0_dp
        self%omega(:) = 0.0_dp
        if (self%tSSR) then
          self%Rab(:,:) = 0.0_dp
        end if

        if (self%Glevel == 1 .or. self%Glevel == 2) then
          if (self%tSaveMem) then
            if (self%isRangeSep) then
              self%HxcHalfS(:,:) = 0.0_dp
              self%HxcHalfD(:,:) = 0.0_dp
            else
              self%HxcSpS(:,:) = 0.0_dp
              self%HxcSpD(:,:) = 0.0_dp
            end if
            self%A1e(:,:) = 0.0_dp
            self%A1ePre(:,:) = 0.0_dp
          end if
        else if (self%Glevel == 3) then
          self%HxcSqrS(:,:,:,:) = 0.0_dp
          self%HxcSqrD(:,:,:,:) = 0.0_dp
          self%Aall(:,:) = 0.0_dp
        end if

        self%XT = 0.0_dp
        if (self%tSSR) then
          self%XTdel(:,:) = 0.0_dp
          self%RdelL(:,:,:,:) = 0.0_dp
          self%ZdelL(:,:,:) = 0.0_dp
          self%Q1del(:,:,:) = 0.0_dp
          self%Q2del(:,:,:) = 0.0_dp
        end if

        self%ZT(:,:) = 0.0_dp
        self%RmatL(:,:,:,:) = 0.0_dp
        self%ZmatL(:,:,:) = 0.0_dp
        self%Q1mat(:,:) = 0.0_dp
        self%Q2mat(:,:) = 0.0_dp

        if (self%tNAC) then
          self%ZTdel(:,:) = 0.0_dp
          self%tmpRL(:,:,:,:) = 0.0_dp
        end if

        self%SSRgrad(:,:,:) = 0.0_dp
        if (self%tNAC) then
          self%SAgrad(:,:,:) = 0.0_dp
          self%SIgrad(:,:,:) = 0.0_dp
          self%avgGrad(:,:) = 0.0_dp
        end if

        if (self%tNAC) then
          self%nacG(:,:,:) = 0.0_dp
          self%nacH(:,:,:) = 0.0_dp
        end if

      end if

    end if


    ! REKS: relaxed density & transition dipole variables

    self%unrelRhoSqr(:,:) = 0.0_dp

    if (self%tTDP .and. self%Lstate == 0) then
      self%unrelTdm(:,:,:) = 0.0_dp
    end if

    if (self%tRD) then
      self%relRhoSqr(:,:) = 0.0_dp
    end if

    if (self%tTDP .and. self%Lstate == 0) then
      self%tdp(:,:) = 0.0_dp
    end if

    if (self%tForces) then

      ! REKS: point charges (QM/MM) variables

      if (self%tExtChrg) then
        self%extCharges(1:3,:) = extChrg(1:3,:)
        ! Adapt to internal sign convention for potentials (electron has positive charge)
        self%extCharges(4,:) = -extChrg(4,:)
        if (self%tBlur) then
          self%blurWidths(:) = blurWidths
        end if
      end if

      ! REKS: stress variables

      if (self%tStress) then
        self%elecStressL(:,:,:) = 0.0_dp
      end if

    end if


    ! REKS: initialize variables

    self%SAweight(:) = 1.0_dp / real(self%SAstates, dp)

    call move_alloc(inp%Tuning, self%Tuning)

    ! Scale up or down the atomic spin constants w.r.t. the systems
    ! iAt : loop for atomic species
    ! order of iAt = order in .gen file
    do iAt = 1, nType
      spinW(:,:,iAt) = self%Tuning(iAt) * spinW(:,:,iAt)
    end do

    ! Set the ordering information betwen R_mat_L and filling_L
    if (self%Efunction > 1 .and. self%tForces) then
      select case (self%reksAlg)
      case (reksTypes%noReks)
      case (reksTypes%ssr22)
        ! R_mat_L has 4 elements and filling_L has 6 elements in (2,2) case
        self%orderRmatL(:) = [1, 2, 1, 2, 3, 4]
      case (reksTypes%ssr44)
        call error("SSR(4,4) is not implemented yet")
      end select
    end if

  contains

    !> Check REKS common requirements
    subroutine setSSR22conditions(self, nEl)

      !> data type for REKS
      type(TReksCalc), intent(inout) :: self

      !> nr. of electrons
      real(dp), intent(in) :: nEl

      self%Nc = int(real(nEl, dp)) / 2 - 1
      self%Na = 2
      self%Lmax = 6
      self%LmaxR = 4
      self%Lpaired = 2
      if (self%Efunction == 1) then
        ! Only PPS state is minimized; single-state REKS
        self%SAstates = 1
        if (.not. self%tAllStates) then
          ! PPS state
          self%nstates = 1
        else
          call error("IncludeAllStates should be set to No in single-state REKS")
        end if
      else if (self%Efunction == 2) then
        ! (PPS+OSS)/2 state is minimized; SA-REKS
        self%SAstates = 2
        if (.not. self%tAllStates) then
          ! PPS, OSS state
          self%nstates = 2
        else
          ! PPS, OSS, DES state
          self%nstates = 3
        end if
      end if
      ! Fractional occupation numbers, n_a and n_b
      allocate(self%FONs(self%Na,1))

    end subroutine setSSR22conditions

    !> Check REKS common requirements
    subroutine checkReksRequirements(self)

      !> data type for REKS
      type(TReksCalc), intent(in) :: self

      ! REKS energy requirements

      if (self%tTDP) then
        if (self%Lstate > 0) then
          call error("Transition dipole is not compatible with L-th microstate")
        else if (self%Efunction == 1) then
          call error("Transition dipole is not compatible with single-state REKS")
        end if
      end if

      if (self%tSSR .and. self%Efunction == 1) then
        call error("Single-state REKS is not SSR state")
      end if

      if (self%shift < -epsilon(1.0_dp)) then
        call error("Too small shift value in REKS")
      end if

      ! REKS gradient requirements

      if (self%tForces) then

        if (self%t3rd) then
          call error("3rd order scc is not compatible with force calculation in REKS")
        end if

        if (self%Lstate > 0) then
          if (self%Efunction == 1) then
            call error("gradient of microstate is not compatible with single-state REKS")
          else if (self%tSSR) then
            call error("For gradient of microstate, please set StateInteractions = No")
          end if
        end if

        if (self%tNAC) then
          if (self%Lstate > 0) then
            call error("Nonadiabatic coupling is not compatible with gradient of microstate")
          else if (.not. self%tSSR .or. self%Efunction == 1) then
            call error("Nonadiabatic coupling is not compatible with &
                & single-state REKS or SA-REKS calculation")
          end if
        end if

      else

        if (self%tNAC) then
          call error("Nonadiabatic coupling requires force evaluation")
        else if (self%tRD) then
          call error("Relaxed density requires force evaluation")
        end if

      end if

      ! REKS stress requirements

      if (self%Efunction /= 1 .and. self%tStress) then
        ! tStress = tForces * tPeriodic * (not tExtChrg)
        call error("Only single-state REKS can evaluate stress")
      end if

      ! REKS system requirements

      if (self%Plevel > 2 .or. self%Plevel < 0) then
        call error("Wrong VerbosityLevel given, please write 0 to 2")
      end if

    end subroutine checkReksRequirements

    !> Check REKS(2,2) requirements
    subroutine checkSSR22Requirements(self)

      !> data type for REKS
      type(TReksCalc), intent(in) :: self

      ! REKS energy requirements

      if (self%rstate > 3 .or. self%rstate < 1) then
        call error("Wrong TargetState given, please write 1 to 3")
      else if (self%Lstate > 6 .or. self%Lstate < 0) then
        call error("Wrong TargetMicrostate given, please write 0 to 6")
      end if

    end subroutine checkSSR22Requirements

  end subroutine REKS_init


  !> Reallocate sparse arrays used in REKS
  subroutine reallocate(self, sparseSize)

    !> data type for REKS
    class(TReksCalc), intent(inout) :: self

    !> Total size of orbitals in the sparse data structures, where the decay of the overlap sets the
    !> sparsity pattern
    integer, intent(in) :: sparseSize

    deallocate(self%getDenseAO)
    if (.not. self%tForces) then
      deallocate(self%rhoSpL)
    end if
    if (.not. self%isRangeSep) then
      deallocate(self%hamSpL)
    end if

    if (self%tForces) then
      deallocate(self%edmSpL)
      if (self%Efunction > 1) then
        if (self%Glevel == 1 .or. self%Glevel == 2) then
          if (self%tSaveMem) then
            if (.not. self%isRangeSep) then
              deallocate(self%HxcSpS)
              deallocate(self%HxcSpD)
            end if
          end if
        end if
      end if
    end if

    allocate(self%getDenseAO(sparseSize,2))
    if (.not. self%tForces) then
      allocate(self%rhoSpL(sparseSize,1,self%Lmax))
    end if
    if (.not. self%isRangeSep) then
      allocate(self%hamSpL(sparseSize,1,self%Lmax))
    end if

    if (self%tForces) then
      allocate(self%edmSpL(sparseSize,self%Lmax))
      if (self%Efunction > 1) then
        if (self%Glevel == 1 .or. self%Glevel == 2) then
          if (self%tSaveMem) then
            if (.not. self%isRangeSep) then
              allocate(self%HxcSpS(sparseSize,sparseSize))
              allocate(self%HxcSpD(sparseSize,sparseSize))
            end if
          end if
        end if
      end if
    end if

    self%getDenseAO = 0
    if (.not. self%tForces) then
      self%rhoSpL = 0.0_dp
    end if
    if (.not. self%isRangeSep) then
      self%hamSpL = 0.0_dp
    end if

    if (self%tForces) then
      self%edmSpL = 0.0_dp
      if (self%Efunction > 1) then
        if (self%Glevel == 1 .or. self%Glevel == 2) then
          if (self%tSaveMem) then
            if (.not. self%isRangeSep) then
              self%HxcSpS = 0.0_dp
              self%HxcSpD = 0.0_dp
            end if
          end if
        end if
      end if
    end if

  end subroutine reallocate


end module dftbp_reksvar
