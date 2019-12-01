!--------------------------------------------------------------------------------------------------!
!  DFTB+: general package for performing fast atomistic simulations                                !
!  Copyright (C) 2006 - 2019  DFTB+ developers group                                               !
!                                                                                                  !
!  See the LICENSE file for terms of usage and distribution.                                       !
!--------------------------------------------------------------------------------------------------!

! TODO
!!!!#:include 'common.fypp'

!> REKS and SI-SA-REKS formulation in DFTB as developed by Lee et al.
!>
!> The functionality of the module has some limitation:
!> * Third order does not work.
!> * Periodic system do not work yet appart from Gamma point.
!> * Orbital potentials or spin-orbit or external E-field does not work yet.
!> * Only for closed shell system.
!> * Onsite corrections are not included in this version
! TODO
!> * Dispersion would be combined with REKS
module dftbp_reksvar

  use dftbp_accuracy
  use dftbp_message
  use dftbp_orbitals
!  use dftbp_sccinit, only : initQFromShellChrg

  implicit none

  private

  public :: TReksIni, TReksCalc
  public :: REKS_init, REKS_reallocate, REKS_destroy

  !> Data type for initial values for REKS calculations
  type :: TReksIni

    !> Is this DFTB/SSR formalism
    logical :: tREKS = .false.

    !> Calculate DFTB/SSR(2,2) formalism
    logical :: tSSR22 = .false.

    !> Calculate DFTB/SSR(4,4) formalism
    logical :: tSSR44 = .false.

    !> REKS: energy variables

    !> Minimized energy functional
    !> In SSR(2,2), 1: PPS, 2: (PPS+OSS)/2
    integer :: Efunction

    !> Calculated energy states in SA-REKS
    !> 1: calculate states only included in 'EnergyFuncitonal'
    !> 2: calculate all possible states
    integer :: Elevel

    !> Calculate SSR state (SI term is included)
    !> 1: calculate SSR state, 0: calculate SA-REKS state
    integer :: useSSR

    !> Target SSR state
    integer :: rstate

    !> Target microstate
    integer :: Lstate

    !> Initial guess for eigenvectors in REKS
    !> 1: diagonalize H0, 2: read external file, 'eigenvec.bin'
    integer :: guess

    !> Maximum iteration used in FON optimization
    integer :: FONmaxIter

    !> Shift value in SCC cycle
    real(dp) :: shift

    !> Scaling constant for atomic spin constants
    real(dp), allocatable :: Tuning(:)

    !> Calculate transition dipole moments
    logical :: tTDP

    !> Swap initial eigenvectors once compared with external file in 1st SCC cycle, 'eigenvec.bin'
    logical :: tAdjustMO

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

    !> Memory level used in calculation of gradient
    !> Usually, 1 is more faster than 2
    !> 1: save 'A_1e' and 'Hxc' variables in memory (requires large memory about 4 orders)
    !> 2: do not save 'A_1e' and 'Hxc' variables in memory
    integer :: Mlevel

  end type TReksIni

  !> Data type for REKS internal settings
  type :: TReksCalc

    !> Is this DFTB/SSR formalism
    logical :: tREKS = .false.

    !> Calculate DFTB/SSR(2,2) formalism
    logical :: tSSR22 = .false.

    !> Calculate DFTB/SSR(4,4) formalism
    logical :: tSSR44 = .false.

    !> REKS: energy variables (input)

    !> Minimized energy functional
    integer :: Efunction

    !> Calculated energy states in SA-REKS
    integer :: Elevel

    !> Calculate SSR state (SI term is included)
    integer :: useSSR

    !> Target SSR state
    integer :: rstate

    !> Target microstate
    integer :: Lstate

    !> Initial guess for eigenvectors in REKS
    integer :: guess

    !> Maximum iteration used in FON optimization
    integer :: FONmaxIter

    !> Shift value in SCC cycle
    real(dp) :: shift

    !> Scaling constant for atomic spin constants
    real(dp), allocatable :: Tuning(:)

    !> Calculate transition dipole moments
    logical :: tTDP

    !> Swap initial eigenvectors once compared with external file in 1st SCC cycle, 'eigenvec.bin'
    logical :: tAdjustMO

    !> REKS: gradient variables (input)

    !> Algorithms to calculate analytic gradients
    integer :: Glevel

    !> Maximum iteration used in calculation of gradient with PCG and CG
    integer :: CGmaxIter

    !> Tolerance used in calculation of gradient with PCG and CG
    real(dp) :: Glimit

    !> Calculate relaxed density of SSR or SA-REKS state
    logical :: tRD

    !> Calculate nonadiabatic coupling vectors
    logical :: tNAC

    !> REKS: system variables (input)

    !> Print level in standard output file
    integer :: Plevel

    !> Memory level used in calculation of gradient
    integer :: Mlevel


    !> REKS: energy variables

    !> Number of microstates
    integer :: Lmax

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
    logical :: tRangeSep

    !> Do we need forces?
    logical :: tForces

    !> if calculation is periodic
    logical :: tPeriodic

    !> Can stress be calculated?
    logical :: tStress

    !> If external charges must be considered
    logical :: tExtChrg

    !> get atom index from AO index
    integer, allocatable :: getAtomIndex(:)

    !> get dense AO index from sparse AO array
    integer, allocatable :: getDenseAO(:,:)

    !> get dense atom index from sparse atom array
    integer, allocatable :: getDenseAtom(:,:)

    !> Dense overlap matrix
    real(dp), allocatable :: over(:,:)

    !> Filling for each microstate
    real(dp), allocatable :: filling_L(:,:,:)

    !> Dense density matrix for each microstate
    real(dp), allocatable :: dm_L(:,:,:,:)

    !> Sparse density matrix for each microstate
    real(dp), allocatable :: dm_sp_L(:,:,:)

    !> Dense delta density matrix for each microstate
    real(dp), allocatable :: Deltadm_L(:,:,:,:)

    !> Mulliken population for each microstate
    real(dp), allocatable :: qOutput_L(:,:,:,:)

    ! TODO : it is not needed
!    !> reference neutral atomic occupations
!    real(dp), allocatable :: q0(:,:,:)

    !> charge per atomic shell for each microstate
    real(dp), allocatable :: chargePerShell_L(:,:,:,:)

    !> internal atom and spin resolved potential
    real(dp), allocatable :: intAtom(:,:)

    !> internal shell and spin resolved potential for each microstate
    real(dp), allocatable :: intShell_L(:,:,:,:)

    !> internal block and spin resolved potential for each microstate
    real(dp), allocatable :: intBlock_L(:,:,:,:,:)

    !> Dense Hamiltonian matrix for each microstate
    real(dp), allocatable :: ham_L(:,:,:,:)

    !> Sparse Hamiltonian matrix for each microstate
    real(dp), allocatable :: ham_sp_L(:,:,:)

    !> Weight for each microstate per state
    real(dp), allocatable :: weight_L(:,:)

    !> Weights used in state-averaging
    real(dp), allocatable :: SAweight(:)

    !> Weight of each microstate for state to be optimized; weight = weight_L * SAweight
    real(dp), allocatable :: weight(:)

    !> Fractional occupation numbers of active orbitals
    real(dp), allocatable :: FONs(:,:)

    !> energy of states
    real(dp), allocatable :: energy(:)

    !> non SCC energy for each microstate
    real(dp), allocatable :: en_L_nonscc(:)

    !> SCC energy for each microstate
    real(dp), allocatable :: en_L_scc(:)

    !> spin-polarized energy for each microstate
    real(dp), allocatable :: en_L_spin(:)

    !> 3rd order SCC energy for each microstate
    real(dp), allocatable :: en_L_3rd(:)

    !> Long-range corrected energy for each microstate
    real(dp), allocatable :: en_L_fock(:)

    !> total energy for each microstate
    real(dp), allocatable :: en_L_tot(:)

    !> dense fock matrix for core orbitals
    real(dp), allocatable :: fock_Fc(:,:)

    !> dense fock matrix for active orbitals
    real(dp), allocatable :: fock_Fa(:,:,:)

    !> dense pseudo-fock matrix
    real(dp), allocatable :: fock(:,:)

    !> eigenvectors from pesudo-fock matrix
    real(dp), allocatable :: eigvecsFock(:,:)

    !> eigenvectors from SA-REKS state
    real(dp), allocatable :: eigvecsSSR(:,:)

    !> Eigenvectors for previous step
    ! TODO : this variable should be removed
!    real(dp), allocatable :: eigvecs_be(:,:)

    !> REKS: gradient variables

    !> Ordering between R_mat_L and filling_L
    integer, allocatable :: orderRmat_L(:)

    !> REKS: relaxed density & transition dipole variables

    !> unrelaxed density matrix for target SSR or SA-REKS state
    real(dp), allocatable :: P_m_o(:,:)

    !> unrelaxed transition density matrix between SSR or SA-REKS states
    real(dp), allocatable :: P_m_del_o(:,:,:)

    !> transition dipole moment between states
    real(dp), allocatable :: tdp(:,:)

  end type TReksCalc

  contains

!  subroutine REKS_init(self, ini, orb, referenceN0, species0, spinW, nSpin,&
!      & nEl, nChrgs, t3rd, tOnSite, tRangeSep, tForces, tPeriodic, tStress)
  subroutine REKS_init(self, ini, orb, spinW, nSpin, nEl,&
      & nChrgs, t3rd, tRangeSep, tForces, tPeriodic, tStress)
    
    !> data type for REKS
    type(TReksCalc), intent(out) :: self

    !> data type for REKS input
    type(TReksIni), intent(inout) :: ini

    !> Atomic orbital information
    type(TOrbitals), intent(in) :: orb

!    !> reference n_0 charges for each atom
!    real(dp), intent(in) :: referenceN0(:,:)
!
!    !> type of the atoms (nAtom)
!    integer, intent(in) :: species0(:)

    !> Spin W values
    real(dp), intent(inout) :: spinW(:,:,:)

    !> Number of spin components, 1 is unpolarised, 2 is polarised, 4 is noncolinear / spin-orbit
    integer, intent(inout) :: nSpin

    !> nr. of electrons
    real(dp), intent(in) :: nEl

    !> Nr. of external charges
    integer, intent(in) :: nChrgs

    !> Third order DFTB
    logical, intent(in) :: t3rd

    !> Whether to run a range separated calculation
    logical, intent(in) :: tRangeSep

    !> Do we need forces?
    logical, intent(in) :: tForces

    !> if calculation is periodic
    logical, intent(in) :: tPeriodic

    !> Can stress be calculated?
    logical, intent(in) :: tStress

    integer :: nOrb, mOrb, mShell, nOrbHalf
    integer :: nstates, SAstates, nstHalf
    integer :: Nc, Na, Nv, superN, Lmax
    integer :: iAt, nAtom, nType

    ! Set REKS input variables

    self%tREKS      = ini%tREKS
    self%tSSR22     = ini%tSSR22
    self%tSSR44     = ini%tSSR44

    self%Efunction  = ini%Efunction
    self%Elevel     = ini%Elevel
    self%useSSR     = ini%useSSR
    self%rstate     = ini%rstate
    self%Lstate     = ini%Lstate
    self%guess      = ini%guess
    self%FONmaxIter = ini%FONmaxIter
    self%shift      = ini%shift

    self%Glevel     = ini%Glevel
    self%CGmaxIter  = ini%CGmaxIter
    self%Glimit     = ini%Glimit

    self%Plevel     = ini%Plevel
    self%Mlevel     = ini%Mlevel

    self%tTDP       = ini%tTDP
    self%tAdjustMO  = ini%tAdjustMO

    self%tRD        = ini%tRD
    self%tNAC       = ini%tNAC

    ! Set REKS variables

    if (self%tSSR22) then

      self%Nc = int(dble(nEl)) / 2 - 1
      self%Na = 2
      self%Lmax = 6
      self%Lpaired = 2
      if (self%Efunction == 1) then
        ! Only PPS state is minimized; single-state REKS
        self%SAstates = 1
        if (self%Elevel == 1) then
          ! PPS state
          self%nstates = 1
        else
          call error("EnergyLevel should be 1 in single-state REKS")
        end if
      else if (self%Efunction == 2) then
        ! (PPS+OSS)/2 state is minimized; SA-REKS
        self%SAstates = 2
        if (self%Elevel == 1) then
          ! PPS, OSS state
          self%nstates = 2
        else if (self%Elevel == 2) then
          ! PPS, OSS, DES state
          self%nstates = 3
        else
          call error("EnergyLevel should be 1 or 2 in SI-SA-REKS")
        end if
      else
        call error("EnergyFunctional should be 1 or 2 in REKS calculation")
      end if
      ! Fractional occupation numbers, n_a and n_b
      allocate(self%FONs(self%Na,1))

    else if (self%tSSR44) then

      call error("SSR(4,4) not implemented yet")

    end if

    nOrb     = orb%nOrb
    mOrb     = orb%mOrb      ! ... s px py pz ...
    mShell   = orb%mShell    ! ... s p ...
    nstates  = self%nstates
    SAstates = self%SAstates

    ! Here, nSpin changes to two for calculation of microstates
    nSpin    = 2
    Lmax     = self%Lmax
    Nc       = self%Nc
    Na       = self%Na
    Nv       = nOrb - Nc - Na
    nOrbHalf = nOrb * (nOrb + 1) / 2 ! original nHalf
    nstHalf  = nstates * (nstates - 1) / 2 ! original nst_h
    superN   = Nc*Nv + Na*(Nc+Nv) + Na*(Na-1)/2

    nAtom  = size(orb%nOrbAtom,dim=1)
    nType  = size(ini%Tuning,dim=1)

    self%t3rd = t3rd
    self%tRangeSep = tRangeSep

    self%tForces = tForces
    self%tPeriodic = tPeriodic
    self%tStress = tStress
    self%tExtChrg = (nChrgs > 0)

    ! Check REKS requirements

    call checkReksRequirements(self)

    if (self%tSSR22) then
      call checkSSR22Requirements(self)
    else if (self%tSSR44) then
      call error("SSR(4,4) not implemented yet")
    end if

    ! Allocate REKS variables

!    allocate(self%q0(mOrb,nAtom,nSpin))
!    self%q0(:,:,:) = 0.0_dp
!    call initQFromShellChrg(self%q0, referenceN0, species0, orb)

    allocate(self%Tuning(nType))

    ! REKS: energy variables

    allocate(self%getAtomIndex(nOrb))
    allocate(self%getDenseAO(0,2))
    allocate(self%getDenseAtom(0,2))

    allocate(self%over(nOrb,nOrb))
    allocate(self%filling_L(nOrb,nSpin,Lmax))

    if (self%tForces) then
      allocate(self%dm_L(nOrb,nOrb,1,Lmax))
    else
      allocate(self%dm_sp_L(0,1,Lmax))
    end if

    if (self%tRangeSep) then
      allocate(self%Deltadm_L(nOrb,nOrb,1,Lmax))
    end if

    allocate(self%qOutput_L(mOrb,nAtom,nSpin,Lmax))
    allocate(self%chargePerShell_L(mShell,nAtom,nSpin,Lmax))

    allocate(self%intAtom(nAtom,nSpin))
    allocate(self%intShell_L(mShell,nAtom,nSpin,Lmax))
    allocate(self%intBlock_L(mOrb,mOrb,nAtom,nSpin,Lmax))

    if (self%tRangeSep) then
      allocate(self%ham_L(nOrb,nOrb,1,Lmax))
    else
      allocate(self%ham_sp_L(0,1,Lmax))
    end if

    allocate(self%weight_L(nstates,Lmax))
    allocate(self%SAweight(SAstates))
    allocate(self%weight(Lmax))

    allocate(self%energy(nstates))

    allocate(self%en_L_nonSCC(Lmax))
    allocate(self%en_L_SCC(Lmax))
    allocate(self%en_L_spin(Lmax))

    if (self%t3rd) then
      allocate(self%en_L_3rd(Lmax))
    end if

    if (self%tRangeSep) then
      allocate(self%en_L_fock(Lmax))
    end if

    allocate(self%en_L_tot(Lmax))

    allocate(self%fock_Fc(nOrb,nOrb))
    allocate(self%fock_Fa(nOrb,nOrb,Na))
    allocate(self%fock(nOrb,nOrb))

    allocate(self%eigvecsFock(nOrb,nOrb))
    allocate(self%eigvecsSSR(nstates,nstates))

!    if (self%tAdjustMO) then
!      allocate(self%eigvecs_be(nOrb,nOrb))
!    end if

    ! REKS: gradient variables

    if (self%tForces) then

      allocate(self%orderRmat_L(Lmax))

    end if

    ! REKS: relaxed density & transition dipole variables

    allocate(self%P_m_o(nOrb,nOrb))

    if (self%tTDP .and. self%Lstate == 0) then
      allocate(self%P_m_del_o(nOrb,nOrb,nstHalf))
    end if

    ! ... REKS: unrelaxed transition dipole variables
    if (self%tTDP .and. self%Lstate == 0) then
      allocate(self%tdp(3,nstHalf))
    end if


    ! REKS: energy variables

    self%getAtomIndex      = 0
    self%getDenseAO        = 0
    self%getDenseAtom      = 0

    self%over              = 0.0_dp
    self%filling_L         = 0.0_dp

    if (self%tForces) then
      self%dm_L = 0.0_dp
    else
      self%dm_sp_L = 0.0_dp
    end if

    if (self%tRangeSep) then
      self%Deltadm_L = 0.0_dp
    end if

    self%qOutput_L         = 0.0_dp
    self%chargePerShell_L  = 0.0_dp

    self%intAtom           = 0.0_dp
    self%intShell_L        = 0.0_dp
    self%intBlock_L        = 0.0_dp

    if (self%tRangeSep) then
      self%ham_L = 0.0_dp
    else
      self%ham_sp_L = 0.0_dp
    end if

    self%weight_L          = 0.0_dp
    self%weight            = 0.0_dp

    self%FONs              = 0.0_dp
    self%energy            = 0.0_dp

    self%en_L_nonSCC       = 0.0_dp
    self%en_L_SCC          = 0.0_dp
    self%en_L_spin         = 0.0_dp

    if (self%t3rd) then
      self%en_L_3rd      = 0.0_dp
    end if

    if (self%tRangeSep) then
      self%en_L_fock   = 0.0_dp
    end if

    self%en_L_tot          = 0.0_dp

    self%fock_Fc           = 0.0_dp
    self%fock_Fa           = 0.0_dp
    self%fock              = 0.0_dp

    self%eigvecsFock       = 0.0_dp
    self%eigvecsSSR        = 0.0_dp

!    if (self%tAdjustMO) then
!      self%eigvecs_be = 0.0_dp
!    end if

    ! REKS: gradient variables

    if (self%tForces) then

    end if

    ! REKS: relaxed density & transition dipole variables

    self%P_m_o     = 0.0_dp

    if (self%tTDP .and. self%Lstate == 0) then
      self%P_m_del_o = 0.0_dp
    end if

    if (self%tTDP .and. self%Lstate == 0) then
      self%tdp = 0.0_dp
    end if

    ! REKS: initialize variables

    self%SAweight(:) = 1.0_dp / dble(self%SAstates)

    self%Tuning(:) = ini%Tuning(:)
    deallocate(ini%Tuning)

    ! Scale up or down the atomic spin constants w.r.t. the systems
    ! iAt : loop for atomic species
    ! order of iAt = order in .gen file
    do iAt = 1, nType
      spinW(:,:,iAt) = self%Tuning(iAt) * spinW(:,:,iAt)
    end do

    ! Set the ordering information betwen R_mat_L and filling_L
    if (self%Efunction > 1 .and. self%tForces) then
      if (self%tSSR22) then
        ! R_mat_L has 4 elements and filling_L has 6 elements in (2,2) case
        self%orderRmat_L(:) = [1, 2, 1, 2, 3, 4]
      else if (self%tSSR44) then
        call error("SSR(4,4) not implemented yet")
      end if
    end if

  contains

    subroutine checkReksRequirements(self)

      !> data type for REKS
      type(TReksCalc), intent(in) :: self

      ! REKS energy requirements

      if (self%tTDP .and. self%Lstate > 0) then
        call error("Transition dipole is not compatible with L-th microstate")
      end if

      ! REKS gradient requirements

      if (self%tForces) then

        if (self%t3rd) then
          call error("3rd order scc is not compatible with force calculation in REKS")
        end if

        if (self%Lstate > 0) then
          if (self%Efunction == 1) then
            call error("gradient of microstate is not compatible with single-state REKS")
          else if (self%useSSR == 1) then
            call error("For gradient of microstate, please set useSSRstate = 0")
          end if
        end if

        if (self%tNAC) then
          if (self%Lstate > 0) then
            call error("Nonadiabatic coupling is not compatible with gradient of microstate")
          else if (self%useSSR == 0 .or. self%Efunction == 1) then
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

      if (self%Efunction /= 1 .and. tStress) then
        ! tStress = tForces * tPeriodic * (not tExtChrg)
        call error("Only single-state REKS can evaluate stress")
      end if

    end subroutine checkReksRequirements

    subroutine checkSSR22Requirements(self)

      !> data type for REKS
      type(TReksCalc), intent(in) :: self

      ! REKS energy requirements

      if (self%rstate > 3 .or. self%rstate < 1) then
        call error("Wrong TargetState given, please write 1 to 3")
      else if (self%Lstate > 6 .or. self%Lstate < 0) then
        call error("Wrong TargetStateL given, please write 0 to 6")
      else if (self%useSSR > 1 .or. self%useSSR < 0) then
        call error("Wrong useSSRstate given, please select 0 or 1")
      else if (self%guess > 2 .or. self%guess < 1) then
        call error("Wrong InitialGuess given, please select 1 or 2")
      end if

      ! REKS gradient requirements

      if (self%tForces) then
        if (self%Glevel > 3 .or. self%Glevel < 1) then
          call error("Wrong GradientLevel option, please write 1 to 3")
        end if
        if (self%Mlevel > 2 .or. self%Mlevel < 1) then
          call error("Wrong memory option, please select 1 or 2")
        end if
      end if

      ! REKS system requirements

      if (self%Plevel > 2 .or. self%Plevel < 0) then
        call error("Wrong printing option, please write 0 to 2")
      end if

    end subroutine checkSSR22Requirements

  end subroutine REKS_init


  subroutine REKS_reallocate(self, sparseSize)

    !> data type for REKS
    type(TReksCalc), intent(inout) :: self

    !> Total size of orbitals in the sparse data structures, where the decay of the overlap sets the
    !> sparsity pattern
    integer, intent(in) :: sparseSize

    deallocate(self%getDenseAO)
    if (.not. self%tForces) deallocate(self%dm_sp_L)
    if (.not. self%tRangeSep) deallocate(self%ham_sp_L)

    allocate(self%getDenseAO(sparseSize,2))
    if (.not. self%tForces) allocate(self%dm_sp_L(sparseSize,1,self%Lmax))
    if (.not. self%tRangeSep) allocate(self%ham_sp_L(sparseSize,1,self%Lmax))

    self%getDenseAO = 0
    if (.not. self%tForces) self%dm_sp_L = 0.0_dp
    if (.not. self%tRangeSep) self%ham_sp_L = 0.0_dp

  end subroutine REKS_reallocate


  subroutine REKS_destroy(self)

    !> data type for REKS
    type(TReksCalc), intent(inout) :: self

    ! Deallocate REKS variables

!    deallocate(self%q0)

    deallocate(self%Tuning)

    ! REKS: energy variables

    deallocate(self%getAtomIndex)
    deallocate(self%getDenseAO)
    deallocate(self%getDenseAtom)

    deallocate(self%over)
    deallocate(self%filling_L)

    if (self%tForces) then
      deallocate(self%dm_L)
    else
      deallocate(self%dm_sp_L)
    end if

    if (self%tRangeSep) then
      deallocate(self%Deltadm_L)
    end if

    deallocate(self%qOutput_L)
    deallocate(self%chargePerShell_L)

    deallocate(self%intAtom)
    deallocate(self%intShell_L)
    deallocate(self%intBlock_L)

    if (self%tRangeSep) then
      deallocate(self%ham_L)
    else
      deallocate(self%ham_sp_L)
    end if

    deallocate(self%weight_L)
    deallocate(self%SAweight)
    deallocate(self%weight)

    deallocate(self%FONs)
    deallocate(self%energy)

    deallocate(self%en_L_nonSCC)
    deallocate(self%en_L_SCC)
    deallocate(self%en_L_spin)

    if (self%t3rd) then
      deallocate(self%en_L_3rd)
    end if

    if (self%tRangeSep) then
      deallocate(self%en_L_fock)
    end if

    deallocate(self%en_L_tot)

    deallocate(self%fock_Fc)
    deallocate(self%fock_Fa)
    deallocate(self%fock)

    deallocate(self%eigvecsFock)
    deallocate(self%eigvecsSSR)

!    if (self%tAdjustMO) then
!      deallocate(self%eigvecs_be)
!    end if

    ! REKS: gradient variables

    if (self%tForces) then

      deallocate(self%orderRmat_L)

    end if

    ! REKS: relaxed density & transition dipole variables

    deallocate(self%P_m_o)

    if (self%tTDP .and. self%Lstate == 0) then
      deallocate(self%P_m_del_o)
    end if

    if (self%tTDP .and. self%Lstate == 0) then
      deallocate(self%tdp)
    end if

  end subroutine REKS_destroy


end module dftbp_reksvar
