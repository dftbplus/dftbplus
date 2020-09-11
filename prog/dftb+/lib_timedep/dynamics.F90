!--------------------------------------------------------------------------------------------------!
!  DFTB+: general package for performing fast atomistic simulations                                !
!  Copyright (C) 2006 - 2020  DFTB+ developers group                                               !
!                                                                                                  !
!  See the LICENSE file for terms of usage and distribution.                                       !
!--------------------------------------------------------------------------------------------------!

#:include 'common.fypp'

!> Implements real-time Ehrenfest time-dependent DFTB by numerically propagating the electronic
!> one-electron density matrix of the system and the nuclei in the presence of an external
!> perturbation (kick or laser) Bonafé, F. P., Aradi, B., Hourahine, B., Medrano, C. R., Hernandez,
!> F. J., Frauenheim, T., & Sánchez, C. G.  Journal of Chemical Theory and Computation (2020)
!> https://doi.org/10.1021/acs.jctc.9b01217

module dftbp_timeprop
  use dftbp_globalenv
  use dftbp_commontypes
  use dftbp_potentials
  use dftbp_scc
  use dftbp_shift
  use dftbp_accuracy
  use dftbp_constants
  use dftbp_sparse2dense
  use dftbp_densitymatrix
  use dftbp_blasroutines
  use dftbp_lapackroutines
  use dftbp_populations
  use dftbp_bondpops
  use dftbp_blas
  use dftbp_lapack
  use dftbp_spin
  use dftbp_forces
  use dftbp_repulsive
  use dftbp_slakocont
  use dftbp_repcont
  use dftbp_thermostat
  use dftbp_mdintegrator
  use dftbp_dummytherm
  use dftbp_mdcommon
  use dftbp_ranlux
  use dftbp_periodic
  use dftbp_velocityverlet
  use dftbp_nonscc
  use dftbp_energytypes, only : TEnergies, TEnergies_init
  use dftbp_getenergies, only : calcEnergies, calcRepulsiveEnergy, calcDispersionEnergy, sumEnergies
  use dftbp_thirdorder, only : TThirdOrder
  use dftbp_solvation, only : TSolvation
  use dftbp_populations
  use dftbp_eigenvects
  use dftbp_sk
  use dftbp_dispiface
  use dftbp_dispersions
  use dftbp_environment
  use dftbp_repcont
  use dftbp_timer
  use dftbp_taggedoutput
  use dftbp_hamiltonian
  use dftbp_elstattypes
  use dftbp_onsitecorrection
  use dftbp_message
  use dftbp_elecsolvers, only : TElectronicSolver
  use dftbp_simplealgebra
  use dftbp_RangeSeparated, only : TRangeSepFunc
  use dftbp_qdepextpotproxy, only : TQDepExtPotProxy
  use dftbp_reks, only : TReksCalc
  implicit none
  private

  public :: runDynamics, TElecDynamics_init
  public :: TElecDynamicsInp, TElecDynamics
  public :: pertTypes, envTypes, tdSpinTypes

  !> Data type to  initialize electronic dynamics variables from parser
  type TElecDynamicsInp
    !> External field peak intensity
    real(dp) :: tdfield

    !> Timestep for propagation
    real(dp) :: dt

    !> Electric field angular frequency
    real(dp) :: omega

    !> Real part of the polarization direction for laser fields
    real(dp) :: reFieldPolVec(3)

    !> Imaginary part of the polarization direction for laser fields
    real(dp) :: imFieldPolVec(3)

    !> Initial time of pulse
    real(dp) :: time0

    !> Final time of pulse
    real(dp) :: time1

    !> Phase applied of the laser field
    real(dp) :: phase

    !> Polarization direction of kicks
    integer :: polDir

    !> Number of steps to propagate
    integer :: steps

    !> Frequency of output writing for charges and populations
    integer :: writeFreq

    !> Should additional data files be printed, covering atom dynamics, charges and energies
    logical :: tdWriteExtras

    !> Type of external perturbation applied
    integer :: pertType

    !> Envelope shape for laser perturbation
    integer :: envType

    !> Spin type of absorption spectrum for spin polarised calculations
    integer :: spType

    !> Frequency of variable dumping to restart file
    integer :: restartFreq

    !> If time-dependent populations should be calculated
    logical :: tPopulations

    !> If calculation should be restarted from dump file
    logical :: tReadRestart

    !> If dump file should be written during the dynamics
    logical :: tWriteRestart

    !> If a dump file is read, should it be ascii (T) or binary (F)
    logical :: tReadRestartAscii = .false.

    !> If a dump file is read, should it be ascii (T) or binary (F)
    logical :: tWriteRestartAscii = .false.

    !> Index of the moved atoms
    integer, allocatable :: indMovedAtom(:)

    !> Index of the excited atoms
    integer, allocatable :: indExcitedAtom(:)

    !> Number of moved atoms
    integer :: nMovedAtom

    !> Number of steps every which an Euler step is applied (to kill numerical noise)
    integer :: eulerFreq

    !> Number of total system snapshots to be saved during pump simulation
    integer :: tdPPFrames

    !> Number of excited atoms
    integer :: nExcitedAtom

    !> Initial and final times to save the snapshots during pump simulation
    real(dp) :: tdPpRange(2)

    !> if forces should be calculated during propagation
    logical :: tForces

    !> if nuclei should be moved during propagation
    logical :: tIons

    !> if velocities are supplied from input
    logical :: tReadMDVelocities

    !> if Euler steps should be done during simulation
    logical :: tEulers

    !> if pairwise bond energy should be calculated and written
    logical :: tBondE

    !> if pairwise bond population should be calculated and written
    logical :: tBondP

    !> if this is a pump trajectory (for a pump-probe simulation)
    logical :: tPump

    !> if this is a probe trajectory (for a pump-probe simulation)
    logical :: tProbe

    !> atomic (initial) kinetic temperature
    real(dp) :: tempAtom

    !> field strength for KickAndLaser perturbations
    real(dp) :: tdLaserField = 0.0_dp

    !> intial atomic velocities if supplied
    real(dp), allocatable :: initialVelocities(:,:)

    !> if initial fillings are provided in an external file
    logical :: tFillingsFromFile

  end type TElecDynamicsInp

  !> Data type for electronic dynamics internal settings
  type TElecDynamics
    private
    real(dp) :: field, dt, omega, time0, time1
    complex(dp) :: fieldDir(3)
    real(dp), allocatable :: tdFunction(:, :), phase
    integer :: nSteps, writeFreq, pertType, envType, spType
    integer :: nAtom, nOrbs, nSpin=1, currPolDir=1, restartFreq
    logical :: tdWriteExtras
    integer, allocatable :: species(:), polDirs(:), speciesAll(:)
    character(mc), allocatable :: speciesName(:)
    logical :: tPopulations, tSpinPol=.false.
    logical :: tReadRestart, tWriteRestart, tRestartAscii, tWriteRestartAscii, tWriteAutotest
    logical :: tLaser = .false., tKick = .false., tKickAndLaser = .false., tEnvFromFile = .false.
    type(TScc), allocatable :: sccCalc
    character(mc) :: autotestTag

    real(dp), allocatable :: initialVelocities(:,:), movedVelo(:,:), movedMass(:,:)
    real(dp) :: mCutoff, skCutoff, laserField
    real(dp), allocatable :: rCellVec(:,:), cellVec(:,:), kPoint(:,:), KWeight(:)
    real(dp), allocatable :: atomEigVal(:,:)
    integer :: nExcitedAtom, nMovedAtom, nSparse, eulerFreq, PpFreq, PpIni, PpEnd
    integer, allocatable :: iCellVec(:), indMovedAtom(:), indExcitedAtom(:)
    logical :: tIons, tForces, ReadMDVelocities, tPump, tProbe, tRealHS
    logical :: isRangeSep
    logical :: FirstIonStep = .true., tEulers = .false., tBondE = .false., tBondP = .false.
    logical :: tPeriodic = .false., tFillingsFromFile = .false.
    logical :: tNetCharges = .false.
    type(TThermostat), allocatable :: pThermostat
    type(TMDIntegrator), allocatable :: pMDIntegrator
    class(TDispersionIface), allocatable :: dispersion
    type(TNonSccDiff), allocatable :: derivator
    type(TParallelKS), allocatable :: parallelKS
    real(dp), allocatable :: latVec(:,:), invLatVec(:,:)
    real(dp), allocatable :: initCoord(:,:)

    !> count of the number of times dynamics has been initialised
    integer :: nDynamicsInit = 0

  end type TElecDynamics

  type :: TDPertTypeEnum
    !> Dirac delta kick to DM
    integer :: kick = 1

    !> Sinusoidal external field
    integer :: laser = 2

    !> No external perturbation, free dynamics
    integer :: noTDPert = 3

    !> Simultaneous kick (at t=0) and external monochromatic field
    integer :: kickAndLaser = 4

  end type TDPertTypeEnum

  !> Container for enumerated available types of perturbation
  type(TDPertTypeEnum), parameter :: pertTypes = TDPertTypeEnum()

  type :: TDEnvelopeFunctionEnum

    !> Constant envelope
    integer :: constant = 1

    !> Gaussian envelope
    integer :: gaussian = 2

    !> Sin^2 envelope
    integer :: sin2 = 3

    !> Read field from file
    integer :: fromFile = 4

  end type TDEnvelopeFunctionEnum

  !> Container for enumerated available types of envelope function
  type(TDEnvelopeFunctionEnum), parameter :: envTypes = TDEnvelopeFunctionEnum()

  type :: TDSpinTypesEnum

    ! only singlet excitations (no change of total spin)
    integer :: singlet = 1

    ! only triplet excitations (with change of total spin = 1)
    integer :: triplet = 2

  end type TDSpinTypesEnum

  !> Container for enumerated types of spin polarized spectra
  type(TDSpinTypesEnum), parameter :: tdSpinTypes = TDSpinTypesEnum()

  !> version number for restart format, please increment if you change the file format (and consider
  !> adding backward compatibility)
  integer, parameter :: tdDumpFormat = 1

  !> Prefix for dump files for restart
  character(*), parameter :: restartFileName = 'tddump'

  !> Used to return runtime diagnostics
  character(len=120) :: error_string

contains

  !> Initialisation of input variables
  subroutine TElecDynamics_init(this, inp, species, speciesName, tWriteAutotest, autotestTag,&
      & randomThermostat, mass, nAtom, skCutoff, mCutoff, atomEigVal, dispersion, nonSccDeriv,&
      & tPeriodic, parallelKS, tRealHS, kPoint, kWeight, isRangeSep)

    !> ElecDynamics instance
    type(TElecDynamics), intent(out) :: this

    !> ElecDynamicsInp instance
    type(TElecDynamicsInp), intent(in) :: inp

    !> label for each atomic chemical species
    character(mc), allocatable, intent(in) :: speciesName(:)

    !> produce tagged output?
    logical, intent(in) :: tWriteAutotest

    !> Tagged output files (machine readable)
    character(*), intent(in) :: autotestTag

    !> self energy (orbital, atom)
    real(dp), intent(in), allocatable :: atomEigVal(:,:)

    !> nr. of atoms
    integer, intent(in) :: nAtom

    ! thermostat object
    type(TRanlux), allocatable, intent(inout) :: randomThermostat

    !> longest range of interactions for which neighbours are required
    real(dp), intent(in) :: mCutoff

    !> Cut off distance for Slater-Koster interactions
    real(dp) :: skCutoff

    !> list of atomic masses
    real(dp) :: mass(:)

    !> dispersion data and calculations
    class(TDispersionIface), allocatable, intent(inout) :: dispersion

    !> Differentiation method for (H^0,S)
    type(TNonSccDiff), intent(in) :: nonSccDeriv

    !> types of the atoms (nAtom)
    integer, intent(in) :: species(:)

    !> if calculation is periodic
    logical, intent(in) :: tPeriodic

    !> dummy thermostat objetct
    type(TDummyThermostat), allocatable :: pDummyTherm

    !> MD Framework
    type(TMDCommon), allocatable :: pMDFrame

    !> Contains (iK, iS) tuples to be processed in parallel by various processor groups
    type(TParallelKS), intent(in) :: parallelKS

    !> H and S are real
    logical, intent(in) :: tRealHS

    !> K-points
    real(dp) :: kPoint(:,:)

    !> weight of the K-Points
    real(dp) :: KWeight(:)

    !> LC correction
    logical, intent(in) :: isRangeSep

    real(dp) :: norm, tempAtom
    logical :: tMDstill
    integer :: iAtom

    this%field = inp%tdField
    this%dt = inp%dt
    this%nSteps = inp%steps
    this%pertType = inp%pertType
    this%envType = inp%envType
    this%spType = inp%spType
    this%tPopulations = inp%tPopulations
    this%tReadRestart = inp%tReadRestart
    this%tWriteRestart = inp%tWriteRestart
    this%tRestartAscii = inp%tReadRestartAscii
    this%tWriteRestartAscii = inp%tWriteRestartAscii
    this%phase = inp%phase
    this%writeFreq = inp%writeFreq
    this%tdWriteExtras = inp%tdWriteExtras
    this%restartFreq = inp%restartFreq
    this%speciesName = speciesName
    this%tFillingsFromFile = inp%tFillingsFromFile
    this%tRealHS = tRealHS
    this%kPoint = kPoint
    this%KWeight = KWeight
    allocate(this%parallelKS, source=parallelKS)

    if (inp%envType /= envTypes%constant) then
      this%time0 = inp%time0
      this%time1 = inp%time1
    end if

    select case(inp%pertType)
    case(pertTypes%laser)
      this%tLaser = .true.
    case(pertTypes%kick)
      this%tKick = .true.
    case(pertTypes%kickAndLaser)
      this%tKick = .true.
      this%tLaser = .true.
      this%laserField = inp%tdLaserField
      this%tKickAndLaser = .true.
    case(pertTypes%noTDPert)
      this%tKick = .false.
      this%tLaser = .false.
    case default
      call error("Wrong type of perturbation.")
    end select

    if (this%tLaser) then
      if (tPeriodic) then
        call warning('Polarization components of the laser in a periodic direction do not work. &
            & Make sure you are polarizing the field in non-periodic directions.')
        if (any(inp%imFieldPolVec > epsilon(1.0_dp))) then
          call warning('Using circular or elliptical polarization with periodic structures might&
              & not work.')
        end if
      end if
      this%omega = inp%omega
      this%fieldDir = inp%reFieldPolVec + imag * inp%imFieldPolVec
      norm = sqrt(dot_product(real(this%fieldDir, dp),real(this%fieldDir, dp)))
      this%fieldDir = this%fieldDir / norm
      allocate(this%tdFunction(3, 0:this%nSteps))
      this%tEnvFromFile = (this%envType == envTypes%fromFile)
      this%indExcitedAtom = inp%indExcitedAtom
      this%nExcitedAtom = inp%nExcitedAtom
    end if

    if (this%tKick) then
      if (inp%polDir == 4) then
        this%polDirs = [1, 2, 3]
      else
        this%polDirs = [inp%polDir]
      end if
    end if

    this%tWriteAutotest = tWriteAutotest
    this%autotestTag = autotestTag

    this%tIons = inp%tIons
    this%tForces = inp%tForces
    this%tEulers = inp%tEulers
    this%eulerFreq = inp%eulerFreq
    this%tBondE = inp%tBondE
    if (this%tBondE .and. .not. this%tRealHS) then
      call error("Real hamiltonian required for bond energies")
    end if
    this%tBondP = inp%tBondP
    this%species = species
    this%tPeriodic = tPeriodic
    this%isRangeSep = isRangeSep

    if (this%tIons) then
      if (.not. this%tRealHS) then
        call error("Ion dynamics is not implemented yet for imaginary Hamiltonians.")
      elseif (isRangeSep) then
        call error("Ion dynamics is not implemented yet for range separated calculations.")
      end if
      this%tForces = .true.
      this%indMovedAtom = inp%indMovedAtom
      this%nMovedAtom = inp%nMovedAtom
      tempAtom = inp%tempAtom
      tMDstill = .false.
      allocate(this%initialVelocities(3, this%nMovedAtom))

      this%ReadMDVelocities = allocated(inp%initialVelocities)
      if (this%ReadMDVelocities) then
        this%initialVelocities(:,:) = inp%initialVelocities
      end if

      allocate(this%movedVelo(3, this%nMovedAtom))
      allocate(this%movedMass(3, this%nMovedAtom))
      this%movedMass(:,:) = spread(mass(this%indMovedAtom), 1, 3)
      allocate(this%pThermostat)
      allocate(pMDFrame)
      call init(pMDFrame, this%nMovedAtom, nAtom, tMDstill)
      allocate(pDummyTherm)
      call init(pDummyTherm, tempAtom, mass(this%indMovedAtom), randomThermostat, pMDFrame)
      call init(this%pThermostat, pDummyTherm)
      allocate(this%derivator, source=nonSccDeriv)
    else
      if (this%tForces) then
        this%nMovedAtom = nAtom
        allocate(this%movedMass(3, nAtom))
        this%movedMass(:,:) = spread(mass, 1, 3)
        allocate(this%derivator, source=nonSccDeriv)
        this%indMovedAtom = [(iAtom, iAtom = 1, nAtom)]
      else
        this%nMovedAtom = 0
      end if
      allocate(this%movedVelo(3, nAtom))
      this%movedVelo(:,:) = 0.0_dp
    end if

    if (this%tIons .or. this%tForces) then
      if (this%nExcitedAtom /= nAtom) then
        if (this%tLaser) then
          call error("Ion dynamics and forces are not implemented for excitation of a subgroup of&
              & atoms")
        else
          this%nExcitedAtom = nAtom
        end if
      end if
    end if

    this%tNetCharges = .false.
    if (allocated(dispersion)) then
      allocate(this%dispersion, source=dispersion)
    #:if WITH_MBD
      select type (dispersion)
      type is (TDispMbd)
        this%tNetCharges = .true.
      end select
    #:endif
    end if

    this%skCutoff = skCutoff
    this%mCutoff = mCutoff
    allocate(this%atomEigVal, source=atomEigVal)

    this%tPump = inp%tPump
    if (inp%tPump) then
      this%PpIni = int(inp%tdPpRange(1)/this%dt)
      this%PpEnd = int(inp%tdPpRange(2)/this%dt)
      this%PpFreq = int((this%PpEnd - this%PpIni) / inp%tdPPFrames)
    end if
    this%tProbe = inp%tProbe
    if (this%tProbe) then
      this%tReadRestart = .true.
      this%tWriteRestart = .false.
    end if
  end subroutine TElecDynamics_init


  !> Driver of time dependent propagation to calculate with either spectrum or laser
  subroutine runDynamics(this, eigvecs, ham, H0, speciesAll, q0, referenceN0, over, filling,&
      & neighbourList, nNeighbourSK, nNeighbourLC, iSquare, iSparseStart, img2CentCell, orb, coord,&
      & spinW, pRepCont, sccCalc, env, tDualSpinOrbit, xi, thirdOrd, solvation, rangeSep,&
      & qDepExtPot, nDftbUFunc, UJ, nUJ, iUJ, niUJ, iAtInCentralRegion, tFixEf, Ef, coordAll,&
      & onSiteElements, skHamCont, skOverCont, latVec, invLatVec, iCellVec, rCellVec, cellVec,&
      & electronicSolver, eigvecsCplx, taggedWriter, refExtPot)

    !> ElecDynamics instance
    type(TElecDynamics) :: this

    !> Real Eigenvectors
    real(dp), intent(inout), allocatable :: eigvecs(:,:,:)

    !> Complex Eigevenctors
    complex(dp), intent(inout), allocatable :: eigvecsCplx(:,:,:)

    !> Sparse non-SCC hamiltonian
    real(dp), intent(in) :: H0(:)

    !> species of all atoms in the system
    integer, intent(in) :: speciesAll(:)

    !> reference atomic occupations
    real(dp), intent(inout) :: q0(:,:,:)

    !> Reference charges from the Slater-Koster file
    real(dp), intent(in) :: referenceN0(:,:)

    !> resulting hamiltonian (sparse)
    real(dp), allocatable, intent(inout) :: ham(:,:)

    !> overlap (sparse)
    real(dp), allocatable, intent(inout) :: over(:)

    !> central atomic coordinates
    real(dp), allocatable, intent(inout) :: coord(:,:)

    !> all atomic coordinates
    real(dp), allocatable, intent(inout) :: coordAll(:,:)

    !> spin constants
    real(dp), allocatable, intent(in) :: spinW(:,:,:)

    !> occupations
    real(dp), intent(inout) :: filling(:,:,:)

    !> Number of neighbours for each of the atoms
    integer, intent(inout) :: nNeighbourSK(:)

    !> Number of neighbours for each of the atoms with the range separated hybrid
    integer, intent(inout), allocatable :: nNeighbourLC(:)

    !> index array for location of atomic blocks in large sparse arrays
    integer, allocatable, intent(inout) :: iSparseStart(:,:)

    !> image atoms to their equivalent in the central cell
    integer, allocatable, intent(inout) :: img2CentCell(:)

    !> Index array for start of atomic block in dense matrices
    integer, intent(in) :: iSquare(:)

    !> list of neighbours for each atom
    type(TNeighbourList), intent(inout) :: neighbourList

    !> repulsive information
    type(TRepCont), intent(in) :: pRepCont

    !> Atomic orbital information
    type(TOrbitals), intent(in) :: orb

    !> SCC module internal variables
    type(TScc), intent(in), allocatable :: sccCalc

    !> Environment settings
    type(TEnvironment), intent(inout) :: env

    type(TSlakoCont), intent(in) :: skHamCont, skOverCont

    !> Is dual spin orbit being used (block potentials)
    logical, intent(in) :: tDualSpinOrbit

    !> Spin orbit constants if required
    real(dp), allocatable, intent(in) :: xi(:,:)

    !> 3rd order settings
    type(TThirdOrder), intent(inout), allocatable :: thirdOrd

    !> Solvation model
    class(TSolvation), allocatable, intent(inout) :: solvation

    !> Range separation contributions
    type(TRangeSepFunc), allocatable, intent(inout) :: rangeSep

    !> Proxy for querying Q-dependant external potentials
    type(TQDepExtPotProxy), intent(inout), allocatable :: qDepExtPot

    !> which DFTB+U functional (if used)
    integer, intent(in), optional :: nDftbUFunc

    !> U-J prefactors in DFTB+U
    real(dp), intent(in), allocatable :: UJ(:,:)

    !> Number DFTB+U blocks of shells for each atom type
    integer, intent(in), allocatable :: nUJ(:)

    !> which shells are in each DFTB+U block
    integer, intent(in), allocatable :: iUJ(:,:,:)

    !> Number of shells in each DFTB+U block
    integer, intent(in), allocatable :: niUJ(:,:)

    !> Atoms over which to sum the total energies
    integer, intent(in) :: iAtInCentralRegion(:)

    !> Whether fixed Fermi level(s) should be used. (No charge conservation!)
    logical, intent(in) :: tFixEf

    !> If tFixEf is .true. contains reservoir chemical potential, otherwise the Fermi levels found
    !> from the given number of electrons
    real(dp), intent(inout) :: Ef(:)

    !> Corrections terms for on-site elements
    real(dp), intent(in), allocatable :: onSiteElements(:,:,:,:)

    !> Lattice vectors if periodic
    real(dp), intent(in) :: latVec(:,:)

    !> Inverse of the lattice vectors
    real(dp), intent(in) :: invLatVec(:,:)

    !> cell vectors in absolute units
    real(dp), intent(in) :: rCellVec(:,:)

    !> Vectors (in units of the lattice constants) to cells of the lattice
    real(dp), intent(in) :: cellVec(:,:)

    !> index of cell in cellVec and rCellVec for each atom
    integer, allocatable, intent(in) :: iCellVec(:)

    !> Electronic solver information
    type(TElectronicSolver), intent(inout) :: electronicSolver

    !> Tagged writer object
    type(TTaggedWriter), intent(inout) :: taggedWriter

    !> Reference external potential (usual provided via API)
    type(TRefExtPot) :: refExtPot

    integer :: iPol, iCall
    logical :: tWriteAutotest

    if (.not. allocated(sccCalc)) then
      call error("SCC calculations are currently required for dynamics")
    end if
    this%sccCalc = sccCalc
    this%speciesAll = speciesAll
    this%nSpin = size(ham(:,:), dim=2)
    if (this%nSpin > 1) then
      call qm2ud(q0)
    end if

    if (this%tRealHS .and. .not. this%isRangeSep) then
      this%nOrbs = size(eigvecs, dim=1)
    else
      this%nOrbs = size(eigvecsCplx, dim=1)
    end if

    this%nAtom = size(coord, dim=2)
    this%latVec = latVec
    this%invLatVec = invLatVec
    this%iCellVec = iCellVec
    this%rCellVec = rCellVec
    this%cellVec = cellVec
    tWriteAutotest = this%tWriteAutotest
    iCall = 1
    if (allocated(this%polDirs)) then
      if (size(this%polDirs) > 1) then
        this%initCoord = coord
      end if
    end if
    if (this%tKick) then
      do iPol = 1, size(this%polDirs)
        this%currPolDir = this%polDirs(iPol)
        ! Make sure only last component enters autotest
        tWriteAutotest = tWriteAutotest .and. (iPol == size(this%polDirs))
        call doDynamics(this, eigvecs, ham, H0, q0, referenceN0, over, filling, neighbourList,&
            & nNeighbourSK, nNeighbourLC, iSquare, iSparseStart, img2CentCell, orb, coord, spinW,&
            & pRepCont, env, tDualSpinOrbit, xi, thirdOrd, solvation, rangeSep, qDepExtPot,&
            & nDftbUFunc, UJ, nUJ, iUJ, niUJ, iAtInCentralRegion, tFixEf, Ef, tWriteAutotest,&
            & coordAll, onSiteElements, skHamCont, skOverCont, iCall, electronicSolver,&
            & eigvecsCplx, taggedWriter, refExtPot)
        iCall = iCall + 1
      end do
    else
      call doDynamics(this, eigvecs, ham, H0, q0, referenceN0, over, filling, neighbourList,&
          & nNeighbourSK, nNeighbourLC, iSquare, iSparseStart, img2CentCell, orb, coord, spinW,&
          & pRepCont, env, tDualSpinOrbit, xi, thirdOrd, solvation, rangeSep, qDepExtPot,&
          & nDftbUFunc, UJ, nUJ, iUJ, niUJ, iAtInCentralRegion, tFixEf, Ef, tWriteAutotest,&
          & coordAll, onSiteElements, skHamCont, skOverCont, iCall, electronicSolver, eigvecsCplx,&
          & taggedWriter, refExtPot)
    end if

  end subroutine runDynamics


  !> Runs the electronic dynamics of the system
  subroutine doDynamics(this, eigvecsReal, ham, H0, q0, referenceN0, over, filling, neighbourList,&
      & nNeighbourSK, nNeighbourLC, iSquare, iSparseStart, img2CentCell, orb, coord, spinW,&
      & pRepCont, env, tDualSpinOrbit, xi, thirdOrd, solvation, rangeSep, qDepExtPot, nDftbUFunc,&
      & UJ, nUJ, iUJ, niUJ, iAtInCentralRegion, tFixEf, Ef, tWriteAutotest, coordAll,&
      & onSiteElements, skHamCont, skOverCont, iCall, electronicSolver, eigvecsCplx, taggedWriter,&
      & refExtPot)

    !> ElecDynamics instance
    type(TElecDynamics) :: this

    !> Real Eigenvectors
    real(dp), intent(inout), allocatable :: eigvecsReal(:,:,:)

    !> Complex Eigevenctors
    complex(dp), intent(inout), allocatable :: eigvecsCplx(:,:,:)

    !> Sparse storage for non-SCC hamiltonian
    real(dp), intent(in) :: H0(:)

    !> reference atomic occupations
    real(dp), intent(inout) :: q0(:,:,:)

    !> Reference charges from the Slater-Koster file
    real(dp), intent(in) :: referenceN0(:,:)

    !> resulting hamiltonian (sparse)
    real(dp), allocatable, intent(inout) :: ham(:,:)

    !> overlap (sparse)
    real(dp), allocatable, intent(inout) :: over(:)

    !> atomic coordinates
    real(dp), allocatable, intent(inout) :: coord(:,:)

    !> all atomic coordinates
    real(dp), allocatable, intent(inout) :: coordAll(:,:)

    !> spin constants
    real(dp), allocatable, intent(in) :: spinW(:,:,:)

    !> occupations
    real(dp), intent(inout) :: filling(:,:,:)

    !> Number of neighbours for each of the atoms
    integer, intent(inout) :: nNeighbourSK(:)

    !> Number of neighbours for each of the atoms with the range separated hybrid
    integer, intent(inout), allocatable :: nNeighbourLC(:)

    !> index array for location of atomic blocks in large sparse arrays
    integer, allocatable, intent(inout) :: iSparseStart(:,:)

    !> image atoms to their equivalent in the central cell
    integer, allocatable, intent(inout) :: img2CentCell(:)

    !> Index array for start of atomic block in dense matrices
    integer, intent(in) :: iSquare(:)

    !> list of neighbours for each atom
    type(TNeighbourList), intent(inout) :: neighbourList

    !> repulsive information
    type(TRepCont), intent(in) :: pRepCont

    !> Atomic orbital information
    type(TOrbitals), intent(in) :: orb

    !> Environment settings
    type(TEnvironment), intent(inout) :: env

    !> Raw H^0 hamiltonian data
    type(TSlakoCont), intent(in) :: skHamCont

    !> Raw overlap data
    type(TSlakoCont), intent(in) :: skOverCont

    !> Is dual spin orbit being used (block potentials)
    logical, intent(in) :: tDualSpinOrbit

    !> which DFTB+U functional (if used)
    integer, intent(in), optional :: nDftbUFunc

    !> U-J prefactors in DFTB+U
    real(dp), intent(in), allocatable :: UJ(:,:)

    !> Number DFTB+U blocks of shells for each atom type
    integer, intent(in), allocatable :: nUJ(:)

    !> which shells are in each DFTB+U block
    integer, intent(in), allocatable :: iUJ(:,:,:)

    !> Number of shells in each DFTB+U block
    integer, intent(in), allocatable :: niUJ(:,:)

    !> Spin orbit constants if required
    real(dp), allocatable, intent(in) :: xi(:,:)

    !> 3rd order settings
    type(TThirdOrder), intent(inout), allocatable :: thirdOrd

    !> Solvation model
    class(TSolvation), allocatable, intent(inout) :: solvation

    !> Range separation contributions
    type(TRangeSepFunc), allocatable, intent(inout) :: rangeSep

    !> Proxy for querying Q-dependant external potentials
    type(TQDepExtPotProxy), intent(inout), allocatable :: qDepExtPot

    !> Atoms over which to sum the total energies
    integer, intent(in) :: iAtInCentralRegion(:)

    !> Whether fixed Fermi level(s) should be used. (No charge conservation!)
    logical, intent(in) :: tFixEf

    !> If tFixEf is .true. contains reservoir chemical potential, otherwise the Fermi levels found
    !> from the given number of electrons
    real(dp), intent(inout) :: Ef(:)

    !> Should autotest data be written?
    logical, intent(in) :: tWriteAutotest

    !> Corrections terms for on-site elements
    real(dp), intent(in), allocatable :: onSiteElements(:,:,:,:)

    !> Number of times this has been called
    integer, intent(in) :: iCall

    !> Electronic solver information
    type(TElectronicSolver), intent(inout) :: electronicSolver

    !> Tagged writer object
    type(TTaggedWriter), intent(inout) :: taggedWriter

    !> Reference external potential (usual provided via API)
    type(TRefExtPot) :: refExtPot

    complex(dp), allocatable :: Ssqr(:,:,:)
    complex(dp), allocatable :: Sinv(:,:,:)
    complex(dp), allocatable, target :: trho(:,:,:)
    complex(dp), allocatable, target :: trhoOld(:,:,:)
    complex(dp), pointer :: rho(:,:,:), rhoOld(:,:,:)
    complex(dp), allocatable :: H1(:,:,:)
    complex(dp), allocatable :: Eiginv(:,:,:), EiginvAdj(:,:,:)
    real(dp) :: qq(orb%mOrb, this%nAtom, this%nSpin), deltaQ(this%nAtom,this%nSpin)
    real(dp) :: dipole(3,this%nSpin), chargePerShell(orb%mShell,this%nAtom,this%nSpin)
    real(dp), allocatable :: rhoPrim(:,:), ham0(:), ErhoPrim(:)
    complex(dp), allocatable :: H1LC(:,:), deltaRho(:,:,:)
    real(dp) :: time, startTime, timeElec
    integer :: dipoleDat, qDat, energyDat, populDat(this%parallelKS%nLocalKS)
    integer :: forceDat, coorDat
    integer :: fdBondPopul, fdBondEnergy
    integer :: iStep, iKS
    type(TPotentials) :: potential
    type(TEnergies) :: energy
    type(TTimer) :: loopTime
    real(dp), allocatable :: qBlock(:,:,:,:), qNetAtom(:)


    complex(dp), allocatable :: RdotSprime(:,:)
    real(dp) :: coordNew(3, this%nAtom), totalForce(3, this%nAtom)
    real(dp), allocatable :: bondWork(:, :)
    real(dp) :: movedAccel(3, this%nMovedAtom), energyKin, new3Coord(3, this%nMovedAtom)
    real(dp) :: occ(this%nOrbs)
    character(sc) :: dumpIdx
    logical :: tProbeFrameWrite
    real(dp) :: lastBondPopul
    real(dp), allocatable :: velInternal(:,:)

    allocate(Ssqr(this%nOrbs,this%nOrbs,this%parallelKS%nLocalKS))
    allocate(Sinv(this%nOrbs,this%nOrbs,this%parallelKS%nLocalKS))
    allocate(trho(this%nOrbs,this%nOrbs,this%parallelKS%nLocalKS))
    allocate(trhoOld(this%nOrbs,this%nOrbs,this%parallelKS%nLocalKS))
    allocate(H1(this%nOrbs,this%nOrbs,this%parallelKS%nLocalKS))
    allocate(RdotSprime(this%nOrbs,this%nOrbs))

    call env%globalTimer%startTimer(globalTimers%elecDynInit)

    iStep = 0
    startTime = 0.0_dp
    timeElec = 0.0_dp
    RdotSprime(:,:) = 0.0_dp

    if (this%tReadRestart) then
      call readRestartFile(trho, trhoOld, coord, this%movedVelo, startTime, this%dt,&
          & restartFileName, this%tRestartAscii)
      call updateH0S(this, Ssqr, Sinv, coord, orb, neighbourList, nNeighbourSK, iSquare,&
          & iSparseStart, img2CentCell, skHamCont, skOverCont, ham, ham0, over, env, rhoPrim,&
          & ErhoPrim, coordAll)
      if (this%tIons) then

        this%initialVelocities(:,:) = this%movedVelo

        this%ReadMDVelocities = .true.
      end if
    else if (iCall > 1 .and. this%tIons) then
      coord(:,:) = this%initCoord
      call updateH0S(this, Ssqr, Sinv, coord, orb, neighbourList, nNeighbourSK, iSquare,&
          & iSparseStart, img2CentCell, skHamCont, skOverCont, ham, ham0, over, env, rhoPrim,&
          & ErhoPrim, coordAll)
      this%initialVelocities(:,:) = this%movedVelo
      this%ReadMDVelocities = .true.
    end if
    if (this%tLaser) then
      call getTDFunction(this, startTime)
    end if

    call initializeTDVariables(this, trho, H1, Ssqr, Sinv, H0, ham0, over, ham, eigvecsReal,&
        & filling, orb, rhoPrim, potential, neighbourList%iNeighbour, nNeighbourSK, iSquare,&
        & iSparseStart, img2CentCell, Eiginv, EiginvAdj, energy, ErhoPrim, skOverCont, qBlock,&
        & qNetAtom,UJ, onSiteElements, eigvecsCplx, H1LC, bondWork, fdBondEnergy, fdBondPopul,&
        & lastBondPopul, time)

    if (this%tPeriodic) then
      call initLatticeVectors(this)
    end if

    call this%sccCalc%updateCoords(env, coordAll, this%speciesAll, neighbourList)
    if (allocated(this%dispersion)) then
      call this%dispersion%updateCoords(env, neighbourList, img2CentCell, coordAll,&
          & this%speciesAll)
      this%mCutOff = max(this%mCutOff, this%dispersion%getRCutOff())
    end if

    call initTDOutput(this, dipoleDat, qDat, energyDat, populDat, forceDat, coorDat)

    call getChargeDipole(this, deltaQ, qq, dipole, q0, trho, Ssqr, coord, iSquare, qBlock, qNetAtom)
    if (allocated(this%dispersion)) then
      call this%dispersion%updateOnsiteCharges(qNetAtom, orb, referenceN0,&
          & this%speciesAll(:this%nAtom), .true.)
    end if

    call updateH(this, H1, ham, over, ham0, this%speciesAll, qq, q0, coord, orb, potential,&
        & neighbourList, nNeighbourSK, iSquare, iSparseStart, img2CentCell, 0, chargePerShell,&
        & spinW, env, tDualSpinOrbit, xi, thirdOrd, qBlock, nDftbUFunc, UJ, nUJ, iUJ, niUJ,&
        & onSiteElements, refExtPot, deltaRho, H1LC, Ssqr, solvation, rangeSep, this%dispersion,&
        & trho)

    if (this%tForces) then
      totalForce(:,:) = 0.0_dp
      call getForces(this, movedAccel, totalForce, trho, H1, Sinv, neighbourList, nNeighbourSK, &
          & img2CentCell, iSparseStart, iSquare, potential, orb, skHamCont, skOverCont, qq, q0, &
          & pRepCont, coordAll, rhoPrim, ErhoPrim, 0, env, rangeSep, deltaRho)
      if (this%tIons) then
        call initIonDynamics(this, coordNew, coord, movedAccel)
      end if
    end if

    ! Apply kick to rho if necessary (in restart case, check it starttime is 0 or not)
    if (this%tKick .and. startTime < this%dt / 10.0_dp) then
      call kickDM(this, trho, Ssqr, Sinv, iSquare, coord)
    end if

    ! had to add the "or tKick" option to override rhoOld if tReadRestart = yes, otherwise it will
    ! be badly initialised
    if (.not.this%tReadRestart .or. (this%tKick .and. startTime < this%dt / 10.0_dp)) then
      ! Initialize electron dynamics
      ! rhoOld is now the GS DM, rho will be the DM at time=dt
      trhoOld(:,:,:) = trho
      call initializePropagator(this, this%dt, trhoOld, trho, H1, Sinv, coordAll, skOverCont, orb,&
          & neighbourList, nNeighbourSK, img2CentCell, iSquare, rangeSep)
    end if

    call getPositionDependentEnergy(this, energy, coordAll, img2CentCell, nNeighbourSK,&
        & neighbourList, pRepCont, iAtInCentralRegion)

    call getTDEnergy(this, energy, rhoPrim, trhoOld, neighbourList, nNeighbourSK, orb,&
        & iSquare, iSparseStart, img2CentCell, ham0, qq, q0, potential, chargePerShell,&
        & energyKin, tDualSpinOrbit, thirdOrd, solvation, rangeSep, qDepExtPot, qBlock,&
        & nDftbUFunc, UJ, nUJ, iUJ, niUJ, xi, iAtInCentralRegion, tFixEf, Ef, onSiteElements)

    ! after calculating the TD function, set initial time to zero for probe simulations
    ! this is to properly calculate the dipole fourier transform after the simulation
    if (this%tProbe) then
      startTime = 0.0_dp
    end if

    call env%globalTimer%stopTimer(globalTimers%elecDynInit)

    ! Main loop
    call env%globalTimer%startTimer(globalTimers%elecDynLoop)
    call loopTime%start()

    write(stdOut, "(A)") 'Starting dynamics'

    rho => trho
    rhoOld => trhoOld

    do iStep = 0, this%nSteps

      time = iStep * this%dt + startTime

      if (this%tWriteRestart .and. iStep > 0 .and. mod(iStep, max(this%restartFreq,1)) == 0) then
        allocate(velInternal(3,size(this%movedVelo, dim=2)))
        if (this%tIons) then
          call state(this%pMDIntegrator, velocities=velInternal)
        else
          velInternal(:,:) = 0.0_dp
        end if
        call writeRestartFile(rho, rhoOld, coord, velInternal, time, this%dt, restartFileName,&
            & this%tWriteRestartAscii)
        deallocate(velInternal)
      end if

      ! WORKAROUND for gfort9 using max() as mod(0,0) can be reached and fails with this compiler
      tProbeFrameWrite = this%tPump .and. (iStep >= this%PpIni) .and. (iStep <= this%PpEnd)&
          & .and. (mod(iStep-this%PpIni, max(this%PpFreq,1)) == 0)
      if (tProbeFrameWrite) then
        write(dumpIdx,'(I0)')int((iStep-this%PpIni)/this%PpFreq)
        allocate(velInternal(3,size(this%movedVelo, dim=2)))
        if (this%tIons) then
          call state(this%pMDIntegrator, velocities=velInternal)
        else
          velInternal(:,:) = 0.0_dp
        end if
        call writeRestartFile(rho, rhoOld, coord, velInternal, time, this%dt,&
            & trim(dumpIdx) // 'ppdump', this%tWriteRestartAscii)
        deallocate(velInternal)
      end if

      if (.not. this%tReadRestart .or. (iStep > 0) .or. this%tProbe) then
        call writeTDOutputs(this, dipoleDat, qDat, energyDat, forceDat, coorDat, fdBondPopul,&
            & fdBondEnergy, time, energy, energyKin, dipole, deltaQ, coord, totalForce, iStep)
      end if

      if (this%tIons) then
        coord(:,:) = coordNew
        call updateH0S(this, Ssqr, Sinv, coord, orb, neighbourList, nNeighbourSK, iSquare,&
            & iSparseStart, img2CentCell, skHamCont, skOverCont, ham, ham0, over, env, rhoPrim,&
            & ErhoPrim, coordAll)
      end if

      call getChargeDipole(this, deltaQ, qq, dipole, q0, rho, Ssqr, coord, iSquare, qBlock,&
          & qNetAtom)
      if (allocated(this%dispersion)) then
        call this%dispersion%updateOnsiteCharges(qNetAtom, orb, referenceN0,&
            & this%speciesAll(:this%nAtom), .true.)
      end if

      call updateH(this, H1, ham, over, ham0, this%speciesAll, qq, q0, coord, orb, potential,&
          & neighbourList, nNeighbourSK, iSquare, iSparseStart, img2CentCell, iStep,&
          & chargePerShell, spinW, env, tDualSpinOrbit, xi, thirdOrd, qBlock, nDftbUFunc, UJ, nUJ,&
          & iUJ, niUJ, onSiteElements, refExtPot, deltaRho, H1LC, Ssqr, solvation, rangeSep,&
          & this%dispersion,rho)

      if (this%tForces) then
        call getForces(this, movedAccel, totalForce, rho, H1, Sinv, neighbourList,&  !F_1
            & nNeighbourSK, img2CentCell, iSparseStart, iSquare, potential, orb, skHamCont, &
            & skOverCont, qq, q0, pRepCont, coordAll, rhoPrim, ErhoPrim, iStep, env, rangeSep,&
            & deltaRho)
      end if

      if (this%tIons) then
        new3Coord(:,:) = coordNew(:, this%indMovedAtom)
        call next(this%pMDIntegrator, movedAccel, new3Coord, this%movedVelo)
        coordNew(:, this%indMovedAtom) = new3Coord
        call getRdotSprime(this, RdotSprime, coordAll, skOverCont, orb, img2CentCell, &
            &neighbourList, nNeighbourSK, iSquare)
        if ((this%tPopulations) .and. (mod(iStep, this%writeFreq) == 0)) then
          call updateBasisMatrices(this, electronicSolver, Eiginv, EiginvAdj, H1, Ssqr)
        end if

        call getPositionDependentEnergy(this, energy, coordAll, img2CentCell, nNeighbourSK,&
            & neighbourList, pRepCont, iAtInCentralRegion)
      end if

      call getTDEnergy(this, energy, rhoPrim, rho, neighbourList, nNeighbourSK, orb, iSquare,&
          & iSparseStart, img2CentCell, ham0, qq, q0, potential, chargePerShell, energyKin,&
          & tDualSpinOrbit, thirdOrd, solvation, rangeSep, qDepExtPot, qBlock, nDftbUFunc, UJ, nUJ,&
          & iUJ, niUJ, xi, iAtInCentralRegion, tFixEf, Ef, onSiteElements)

      if ((mod(iStep, this%writeFreq) == 0)) then
        call getBondPopulAndEnergy(this, bondWork, lastBondPopul, rhoPrim, ham0, over,&
            & neighbourList%iNeighbour, nNeighbourSK, iSparseStart, img2CentCell, iSquare,&
            & fdBondEnergy, fdBondPopul, tWriteAutotest, time)
      end if

      do iKS = 1, this%parallelKS%nLocalKS
        if (this%tIons .or. (.not. this%tRealHS) .or. this%isRangeSep) then
          H1(:,:,iKS) = RdotSprime + imag * H1(:,:,iKS)

          if (this%tEulers .and. (iStep > 0) .and. (mod(iStep, max(this%eulerFreq,1)) == 0)) then
            call zcopy(this%nOrbs*this%nOrbs, rho(:,:,iKS), 1, rhoOld(:,:,iKS), 1)
            call propagateRho(this, rhoOld(:,:,iKS), rho(:,:,iKS), H1(:,:,iKS), Sinv(:,:,iKS),&
                & this%dt)
          else
            call propagateRho(this, rhoOld(:,:,iKS), rho(:,:,iKS), H1(:,:,iKS), Sinv(:,:,iKS),&
                & 2.0_dp * this%dt)
          end if
        else
          call propagateRhoRealH(this, rhoOld(:,:,iKS), rho(:,:,iKS), H1(:,:,iKS), Sinv(:,:,iKS),&
              & 2.0_dp * this%dt)
        end if
      end do

      if (mod(iStep, 2) == 1) then
        rho => trho
        rhoOld => trhoOld
      else
        rho => trhoOld
        rhoOld => trho
      end if

      if ((this%tPopulations) .and. (mod(iStep, this%writeFreq) == 0)) then
        do iKS = 1, this%parallelKS%nLocalKS
          call getTDPopulations(this, occ, rho, Eiginv, EiginvAdj, populDat, time, iKS)
        end do
      end if

      if (mod(iStep, max(this%nSteps / 10, 1)) == 0) then
        call loopTime%stop()
        timeElec  = loopTime%getWallClockTime()
        write(stdOut, "(A,2x,I6,2(2x,A,F10.6))") 'Step ', iStep, 'elapsed loop time: ',&
            & timeElec, 'average time per loop ', timeElec / (iStep + 1)
      end if
    end do

    write(stdOut, "(A)") 'Dynamics finished OK!'
    call env%globalTimer%stopTimer(globalTimers%elecDynLoop)

    if (tWriteAutotest) then
      call writeTDAutotest(this, dipole, energy, deltaQ, coord, totalForce, occ, lastBondPopul,&
          & taggedWriter)
    end if

    call closeTDOutputs(this, dipoleDat, qDat, energyDat, populDat, forceDat, coorDat, fdBondPopul,&
        & fdBondEnergy)

    if (this%tIons) then
      if (allocated(this%polDirs)) then
        if (size(this%polDirs) <  (this%nDynamicsInit + 1)) then
          deallocate(this%pMDIntegrator)
        end if
      end if
    end if

  end subroutine doDynamics


  !> Updates the hamiltonian with SCC and external TD field (if any) contributions
  subroutine updateH(this, H1, ham, over, H0, speciesAll, qq, q0, coord, orb, potential,&
      & neighbourList, nNeighbourSK, iSquare, iSparseStart, img2CentCell, iStep, chargePerShell,&
      & spinW, env, tDualSpinOrbit, xi, thirdOrd, qBlock, nDftbUFunc, UJ, nUJ, iUJ, niUJ,&
      & onSiteElements, refExtPot, deltaRho, H1LC, Ssqr, solvation, rangeSep, dispersion, rho)

    !> ElecDynamics instance
    type(TElecDynamics) :: this

    !> Square hamiltonian at each spin and k-point
    complex(dp), intent(inout) :: H1(:,:,:)

    !> resulting hamiltonian (sparse)
    real(dp), intent(inout) :: ham(:,:)

    !> overlap (sparse)
    real(dp), intent(inout) :: over(:)

    !> Sparse storage for non-SCC hamiltonian
    real(dp), intent(in) :: H0(:)

    !> species of all atoms in the system
    integer, intent(in) :: speciesAll(:)

    !> atomic ocupations
    real(dp), intent(inout) :: qq(:,:,:)

    !> reference atomic occupations
    real(dp), intent(inout) :: q0(:,:,:)

    !> atomic coordinates
    real(dp), allocatable, intent(inout) :: coord(:,:)

    !> Atomic orbital information
    type(TOrbitals), intent(in) :: orb

    !> potential acting on the system
    type(TPotentials), intent(inout) :: potential

    !> list of neighbours for each atom
    type(TNeighbourList), intent(inout) :: neighbourList

    !> Number of neighbours for each of the atoms
    integer, intent(in) :: nNeighbourSK(:)

    !> Index array for start of atomic block in dense matrices
    integer, intent(in) :: iSquare(:)

    !> index array for location of atomic blocks in large sparse arrays
    integer, intent(in) :: iSparseStart(0:,:)

    !> image atoms to their equivalent in the central cell
    integer, intent(in) :: img2CentCell(:)

    !> current step of the propagation
    integer, intent(in) :: iStep

    !> electrons in each atomi shell
    real(dp), intent(inout) :: chargePerShell(:,:,:)

    !> spin constants
    real(dp), allocatable, intent(in) :: spinW(:,:,:)

    !> Environment settings
    type(TEnvironment), intent(inout) :: env

    !> Is dual spin orbit being used (block potentials)
    logical, intent(in) :: tDualSpinOrbit

    !> Spin orbit constants if required
    real(dp), allocatable, intent(in) :: xi(:,:)

    !> 3rd order settings
    type(TThirdOrder), intent(inout), allocatable :: thirdOrd

    !> block (dual) atomic populations
    real(dp), intent(inout), allocatable :: qBlock(:,:,:,:)

    !> which DFTB+U functional (if used)
    integer, intent(in), optional :: nDftbUFunc

    !> U-J prefactors in DFTB+U
    real(dp), intent(in), allocatable :: UJ(:,:)

    !> Number DFTB+U blocks of shells for each atom type
    integer, intent(in), allocatable :: nUJ(:)

    !> which shells are in each DFTB+U block
    integer, intent(in), allocatable :: iUJ(:,:,:)

    !> Number of shells in each DFTB+U block
    integer, intent(in), allocatable :: niUJ(:,:)

    !> Corrections terms for on-site elements
    real(dp), intent(in), allocatable :: onSiteElements(:,:,:,:)

    !> Reference external potential (usual provided via API)
    type(TRefExtPot) :: refExtPot

    !> Change in density matrix
    complex(dp), allocatable, intent(inout) :: deltaRho(:,:,:)

    !> LC contribution to hamiltonian
    complex(dp), allocatable, intent(inout) :: H1LC(:,:)

    !> Square overlap
    complex(dp), intent(in) :: Ssqr(:,:,:)

    !> Solvation model
    class(TSolvation), allocatable, intent(inout) :: solvation

    !> Range separation contributions
    type(TRangeSepFunc), allocatable, intent(inout) :: rangeSep

    !> dispersion data and calculations
    class(TDispersionIface), allocatable, intent(in) :: dispersion

    !> Density matrix
    complex(dp), intent(in) :: rho(:,:,:)

    real(dp), allocatable :: qiBlock(:,:,:,:) ! not allocated since no imaginary ham
    real(dp), allocatable :: iHam(:,:) ! not allocated since no imaginary ham
    real(dp), allocatable :: T2(:,:)
    integer :: iAtom, iEatom, iSpin, iKS, iK
    logical :: tDFTBU, tImHam

    ! left over from Poisson shift upload from transport being messy
    real(dp), allocatable :: dummy(:,:)

    allocate(T2(this%nOrbs,this%nOrbs))

    ham(:,:) = 0.0_dp

    if (this%nSpin == 2) then
      call ud2qm(qq)
      call ud2qm(q0)
    end if

    tDFTBU = .false.
    if (size(UJ) > 0) then
      tDFTBU = .true.
    end if
    tImHam = .false. ! for the moment

    call resetExternalPotentials(refExtPot, potential)
    call resetInternalPotentials(tDualSpinOrbit, xi, orb, speciesAll, potential)

    call getChargePerShell(qq, orb, speciesAll, chargePerShell)
    call addChargePotentials(env, this%sccCalc, qq, q0, chargePerShell, orb, speciesAll,&
        & neighbourList, img2CentCell, spinW, solvation, thirdOrd, potential,&
        & elstatTypes%gammaFunc, .false., .false., dummy, dispersion)

    if ((size(UJ) /= 0) .or. allocated(onSiteElements)) then
      ! convert to qm representation
      call ud2qm(qBlock)
    end if

    if (size(UJ) /= 0) then
      call addBlockChargePotentials(qBlock, qiBlock, tDftbU, .false., speciesAll, orb, nDftbUFunc,&
          & UJ, nUJ, iUJ, niUJ, potential)
    end if
    if (allocated(onSiteElements)) then
      call addOnsShift(potential%intBlock, potential%iOrbitalBlock, qBlock, qiBlock, q0,&
          & onSiteElements, speciesAll, orb)
    end if

    ! Add time dependent field if necessary
    if (this%tLaser) then
      do iAtom = 1, this%nExcitedAtom
        iEatom = this%indExcitedAtom(iAtom)
        potential%extAtom(iEatom, 1) = dot_product(coord(:,iEatom), this%tdFunction(:, iStep))
      end do
      call total_shift(potential%extShell, potential%extAtom, orb, speciesAll)
      call total_shift(potential%extBlock, potential%extShell, orb, speciesAll)
    end if

    potential%intBlock = potential%intBlock + potential%extBlock
    potential%intShell = potential%intShell + potential%extShell

    call getSccHamiltonian(H0, over, nNeighbourSK, neighbourList, speciesAll, orb, iSparseStart,&
        & img2CentCell, potential, .false., ham, iHam)

    ! Hack due to not using Pauli-type structure outside of this part of the routine
    if (this%nSpin == 2) then
      ham = 2.0_dp * ham
      call qm2ud(ham)
      call qm2ud(q0)
      call qm2ud(qq)
    end if

    do iKS = 1, this%parallelKS%nLocalKS
      iK = this%parallelKS%localKS(1, iKS)
      iSpin = this%parallelKS%localKS(2, iKS)
      if (this%tRealHS) then
        call unpackHS(T2, ham(:,iSpin), neighbourList%iNeighbour, nNeighbourSK, iSquare,&
            & iSparseStart, img2CentCell)
        call blockSymmetrizeHS(T2, iSquare)
        H1(:,:,iSpin) = cmplx(T2, 0.0_dp, dp)
      else
        call unpackHS(H1(:,:,iKS), ham(:,iSpin), this%kPoint(:,iK), neighbourList%iNeighbour,&
            & nNeighbourSK, this%iCellVec, this%cellVec, iSquare, iSparseStart, img2CentCell)
        call blockHermitianHS(H1(:,:,iKS), iSquare)
      end if
    end do

    ! add LC correction
    if (this%isRangeSep) then
      deltaRho = rho
      select case(this%nSpin)
      case(2)
        do iSpin = 1, this%nSpin
          call denseSubtractDensityOfAtoms(q0, iSquare, deltaRho, iSpin)
        end do
      case(1)
        call denseSubtractDensityOfAtoms(q0, iSquare, deltaRho)
      case default
        call error("Range separation not implemented for non-colinear spin")
      end select
      do iSpin = 1, this%nSpin
        H1LC(:,:) = (0.0_dp, 0.0_dp)
        call rangeSep%addLrHamiltonianMatrixCmplx(iSquare, sSqr(:,:,iSpin), deltaRho(:,:,iSpin),&
            & H1LC)
        call blockHermitianHS(H1LC, iSquare)
        H1(:,:,iSpin) = H1(:,:,iSpin) + H1LC
      end do
    end if

  end subroutine updateH


  !> Kick the density matrix for spectrum calculations
  subroutine kickDM(this, rho, Ssqr, Sinv, iSquare, coord)

    !> ElecDynamics instance
    type(TElecDynamics), intent(in) :: this

    !> Square overlap
    complex(dp), intent(in) :: Ssqr(:,:,:)

    !> Square overlap inverse
    complex(dp), intent(in) :: Sinv(:,:,:)

    !> Density matrix
    complex(dp), intent(inout) :: rho(:,:,:)

    !> atomic coordinates
    real(dp), allocatable, intent(in) :: coord(:,:)

    !> Index array for start of atomic block in dense matrices
    integer, intent(in) :: iSquare(:)

    complex(dp), allocatable :: T1(:, :, :), T2(:, :), T3(:, :, :), T4(:, :)
    integer :: iAt, iStart, iEnd, iKS, iSpin, iOrb
    real(dp) :: pkick(this%nSpin)

    character(1), parameter :: localDir(3) = ['x', 'y', 'z']

    allocate(T1(this%nOrbs, this%nOrbs, this%parallelKS%nLocalKS))
    allocate(T2(this%nOrbs, this%nOrbs))
    allocate(T3(this%nOrbs, this%nOrbs, this%parallelKS%nLocalKS))
    allocate(T4(this%nOrbs, this%nOrbs))

    T1(:,:,:) = cmplx(0,0,dp)
    T2(:,:) = cmplx(0,0,dp)
    T3(:,:,:) = cmplx(0,0,dp)
    T4(:,:) = cmplx(0,0,dp)

    pkick(1) = this%field

    if (this%nSpin == 2) then
      select case(this%spType)
      case (tdSpinTypes%singlet)
        pkick(2) = pkick(1)
      case(tdSpinTypes%triplet)
        pkick(2) = -pkick(1)
      end select
    end if

    do iKS = 1, this%parallelKS%nLocalKS
      iSpin = this%parallelKS%localKS(2, iKS)
      do iAt = 1, this%nAtom
        iStart = iSquare(iAt)
        iEnd = iSquare(iAt + 1) - 1
        do iOrb = iStart, iEnd
          T1(iOrb, iOrb, iKS) = exp(cmplx(0, -pkick(iSpin) * coord(this%currPolDir, iAt), dp))
          T3(iOrb, iOrb, iKS) = exp(cmplx(0,  pkick(iSpin) * coord(this%currPolDir, iAt), dp))
        end do
      end do
    end do

    do iKS = 1, this%parallelKS%nLocalKS
      call gemm(T2, T1(:,:,iKS), rho(:,:,iKS))
      call gemm(T4, T2, Ssqr(:,:,iKS), cmplx(1, 0, dp))
      call gemm(T2, T4, T3(:,:,iKS))
      call gemm(rho(:,:,iKS), T2, Sinv(:,:,iKS), cmplx(0.5, 0, dp))
      call gemm(rho(:,:,iKS), Sinv(:,:,iKS), T2, cmplx(0.5, 0, dp), cmplx(1, 0, dp), 'N', 'C')
    end do

    write(stdout,"(A)")'Density kicked along ' // localDir(this%currPolDir) //'!'

  end subroutine kickDM


  !> Creates array for external TD field
  subroutine getTDFunction(this, startTime)

    !> ElecDynamics instance
    type(TElecDynamics), intent(inout) :: this

    !> starting time of the simulation, if relevant
    real(dp), intent(in) :: startTime

    real(dp) :: midPulse, deltaT, angFreq, E0, time, envelope
    real(dp) :: tdfun(3)
    integer :: iStep, laserDat

    midPulse = (this%time0 + this%time1)/2.0_dp
    deltaT = this%time1 - this%time0
    angFreq = this%omega
    E0 = this%field
    if (this%tKickAndLaser) then
      E0 = this%laserField
    end if
    this%tdFunction(:,:) = 0.0_dp
    if (this%tEnvFromFile) then
      E0 = 0.0_dp !this is to make sure we never sum the current field with the read from file
    end if

    call openFile(this, laserDat, 'laser.dat')
    if (.not. this%tEnvFromFile) then
      write(laserDat, "(A)") "#     time (fs)  |  E_x (eV/ang)  | E_y (eV/ang) | E_z (eV/ang)"
    end if

    do iStep = 0,this%nSteps
      time = iStep * this%dt + startTime

      if (this%envType == envTypes%constant) then
        envelope = 1.0_dp
      else if (this%envType == envTypes%gaussian) then
        envelope = exp(-4.0_dp*pi*(time-midPulse)**2 / deltaT**2)
      else if (this%envType == envTypes%sin2 .and. time >= this%time0 .and. time <= this%time1) then
        envelope = sin(pi*(time-this%time0)/deltaT)**2
      else
        envelope = 0.0_dp
      end if

      if (this%tEnvFromFile) then
        read(laserDat, *)time, tdfun(1), tdfun(2), tdfun(3)
        this%tdFunction(:, iStep) = tdfun * (Bohr__AA / Hartree__eV)
      else
        this%tdFunction(:, iStep) = E0 * envelope * aimag(exp(imag*(time*angFreq + this%phase))&
            & * this%fieldDir)
        write(laserDat, "(5F15.8)") time * au__fs,&
            & this%tdFunction(:, iStep) * (Hartree__eV / Bohr__AA)
      end if

    end do

    close(laserDat)

  end subroutine getTDFunction


  !> Calculate charges, dipole moments
  subroutine getChargeDipole(this, deltaQ, qq, dipole, q0, rho, Ssqr, coord, iSquare, qBlock,&
      & qNetAtom)

    !> ElecDynamics instance
    type(TElecDynamics), intent(in) :: this

    !> Negative gross charge
    real(dp), intent(out) :: deltaQ(:,:)

    !> Dipole moment
    real(dp), intent(out) :: dipole(:,:)

    !> atomic ocupations
    real(dp), intent(out) :: qq(:,:,:)

    !> reference atomic occupations
    real(dp), intent(in) :: q0(:,:,:)

    !> atomic coordinates
    real(dp), intent(in) :: coord(:,:)

    !> Density matrix
    complex(dp), intent(in) :: rho(:,:,:)

    !> Square overlap matrix
    complex(dp), intent(in) :: Ssqr(:,:,:)

    !> Index array for start of atomic block in dense matrices
    integer, intent(in) :: iSquare(:)

    !> Mulliken block charges
    real(dp), allocatable, intent(inout) :: qBlock(:,:,:,:)

    !> Net (on-site only) atomic charge
    real(dp), allocatable, intent(inout) :: qNetAtom(:)

    integer :: iAt, iSpin, iOrb1, iOrb2, nOrb, iKS, iK, ii

    qq(:,:,:) = 0.0_dp
    if (this%tRealHS) then

      do iSpin = 1, this%nSpin
        do iAt = 1, this%nAtom
          iOrb1 = iSquare(iAt)
          iOrb2 = iSquare(iAt+1)-1
          ! hermitian transpose used as real part only is needed
          qq(:iOrb2-iOrb1+1,iAt,iSpin) = real(sum(&
              & rho(:,iOrb1:iOrb2,iSpin)*Ssqr(:,iOrb1:iOrb2,iSpin), dim=1), dp)
        end do
      end do

    else

      do iKS = 1, this%parallelKS%nLocalKS
        iK = this%parallelKS%localKS(1, iKS)
        iSpin = this%parallelKS%localKS(2, iKS)

        do iAt = 1, this%nAtom
          iOrb1 = iSquare(iAt)
          iOrb2 = iSquare(iAt+1)-1
          ! only real part is needed
          qq(:iOrb2-iOrb1+1,iAt,iSpin) = qq(:iOrb2-iOrb1+1,iAt,iSpin) + this%kWeight(iK) * &
              & real(sum(rho(:,iOrb1:iOrb2,iKS) * conjg(Ssqr(:,iOrb1:iOrb2,iKS)), dim=1), dp)
        end do
      end do

    end if

    deltaQ(:,:) = sum((qq - q0), dim=1)
    dipole(:,:) = -matmul(coord, deltaQ)

    if (allocated(qBlock)) then
      if (.not. this%tRealHS) then
        call error("Not implemented yet")
      end if
      qBlock(:,:,:,:) = 0.0_dp
      do iKS = 1, this%parallelKS%nLocalKS
        iK = this%parallelKS%localKS(1, iKS)
        iSpin = this%parallelKS%localKS(2, iKS)
        do iAt = 1, this%nAtom
          iOrb1 = iSquare(iAt)
          iOrb2 = iSquare(iAt+1)
          nOrb = iOrb2 - iOrb1
          qBlock(:nOrb,:nOrb,iAt,iSpin) = qBlock(:nOrb,:nOrb,iAt,iSpin) + this%kWeight(iK) *&
              & real(matmul(Ssqr(iOrb1:iOrb2-1,:,iSpin), rho(:,iOrb1:iOrb2-1,iSpin)), dp)
        end do
      end do
      do iAt = 1, this%nAtom
        iOrb1 = iSquare(iAt)
        iOrb2 = iSquare(iAt+1)
        nOrb = iOrb2 - iOrb1
        qBlock(:nOrb,:nOrb,iAt,iSpin) = 0.5_dp * (qBlock(:nOrb,:nOrb,iAt,iSpin)&
            & + transpose(qBlock(:nOrb,:nOrb,iAt,iSpin)) )
      end do
    end if

    if (allocated(qNetAtom)) then
      qNetAtom(:) = 0.0_dp
      do iKS = 1, this%parallelKS%nLocalKS
        iK = this%parallelKS%localKS(1, iKS)
        do iAt = 1, this%nAtom
          iOrb1 = iSquare(iAt)
          iOrb2 = iSquare(iAt+1)-1
          do ii = iOrb1, iOrb2
            ! only real part is needed
            qNetAtom(iAt) = qNetAtom(iAt) + this%kWeight(iK) * real(rho(ii,ii,iKS))
          end do
        end do
      end do
    end if

  end subroutine getChargeDipole


  !> Calculate energy - modify to include new way to calculate energy
  !> Repulsive energy and dispersion energies must be calculated before calling this subroutine
  subroutine getTDEnergy(this, energy, rhoPrim, rho, neighbourList, nNeighbourSK, orb, iSquare,&
      & iSparseStart, img2CentCell, ham0, qq, q0, potential, chargePerShell, energyKin,&
      & tDualSpinOrbit, thirdOrd, solvation, rangeSep, qDepExtPot, qBlock, nDftbUFunc, UJ, nUJ,&
      & iUJ, niUJ, xi, iAtInCentralRegion, tFixEf, Ef, onSiteElements)

    !> ElecDynamics instance
    type(TElecDynamics), intent(inout) :: this

    !> data type for energy components and total
    type(TEnergies), intent(inout) :: energy

    !> sparse density matrix
    real(dp), allocatable, intent(inout) :: rhoPrim(:,:)

    !> Density matrix
    complex(dp), intent(in) :: rho(:,:,:)

    !> Sparse storage for non-SCC hamiltonian
    real(dp), intent(in) :: ham0(:)

    !> atomic ocupations
    real(dp), intent(inout) :: qq(:,:,:)

    !> reference atomic occupations
    real(dp), intent(in) :: q0(:,:,:)

    !> Atomic orbital information
    type(TOrbitals), intent(in) :: orb

    !> neighbour list
    type(TNeighbourList), intent(in) :: neighbourList

    !> Number of neighbours for each of the atoms
    integer, intent(in) :: nNeighbourSK(:)

    !> Index array for start of atomic block in dense matrices
    integer, intent(in) :: iSquare(:)

    !> index array for location of atomic blocks in large sparse arrays
    integer, intent(in) :: iSparseStart(0:,:)

    !> image atoms to their equivalent in the central cell
    integer, intent(in) :: img2CentCell(:)

    !> electrons in each atomi shell
    real(dp), intent(in) :: chargePerShell(:,:,:)

    !> potential acting on the system
    type(TPotentials), intent(in) :: potential

    !> kinetic energy
    real(dp), intent(out) :: energyKin

    !> Is dual spin orbit being used
    logical, intent(in) :: tDualSpinOrbit

    !> 3rd order settings
    type(TThirdOrder), intent(inout), allocatable :: thirdOrd

    !> Solvation model
    class(TSolvation), allocatable, intent(inout) :: solvation

    !> Range separation contributions
    type(TRangeSepFunc), allocatable, intent(inout) :: rangeSep

    !> Proxy for querying Q-dependant external potentials
    type(TQDepExtPotProxy), intent(inout), allocatable :: qDepExtPot

    !> block (dual) atomic populations
    real(dp), intent(in), allocatable :: qBlock(:,:,:,:)

    !> which DFTB+U functional (if used)
    integer, intent(in), optional :: nDftbUFunc

    !> U-J prefactors in DFTB+U
    real(dp), intent(in), allocatable :: UJ(:,:)

    !> Number DFTB+U blocks of shells for each atom type
    integer, intent(in), allocatable :: nUJ(:)

    !> which shells are in each DFTB+U block
    integer, intent(in), allocatable :: iUJ(:,:,:)

    !> Number of shells in each DFTB+U block
    integer, intent(in), allocatable :: niUJ(:,:)

    !> Spin orbit constants
    real(dp), intent(in), allocatable :: xi(:,:)

    !> Atoms over which to sum the total energies
    integer, intent(in) :: iAtInCentralRegion(:)

    !> Whether fixed Fermi level(s) should be used. (No charge conservation!)
    logical, intent(in) :: tFixEf

    !> If tFixEf is .true. contains reservoir chemical potential, otherwise the Fermi levels found
    !> from the given number of electrons
    real(dp), intent(inout) :: Ef(:)

    !> Corrections terms for on-site elements
    real(dp), intent(in), allocatable :: onSiteElements(:,:,:,:)

    real(dp), allocatable :: qiBlock(:,:,:,:) ! never allocated
    integer :: iKS, iK, iSpin
    real(dp) :: TS(this%nSpin)
    logical :: tDFTBU
    type(TReksCalc), allocatable :: reks ! never allocated

    ! if Forces are calculated, rhoPrim has already been calculated
    ! check allways that calcEnergy is called AFTER getForces
    if (.not. this%tForces) then
      rhoPrim(:,:) = 0.0_dp
      do iKS = 1, this%parallelKS%nLocalKS
        iSpin = this%parallelKS%localKS(2, iKS)
        if (this%tRealHS) then
          call packHS(rhoPrim(:,iSpin), real(rho(:,:,iSpin), dp), neighbourlist%iNeighbour,&
              & nNeighbourSK, orb%mOrb, iSquare, iSparseStart, img2CentCell)
        else
          iK = this%parallelKS%localKS(1, iKS)
          call packHS(rhoPrim(:,iSpin), rho(:,:,iKS), this%kPoint(:,iK), this%kWeight(iK),&
              & neighbourList%iNeighbour, nNeighbourSK, orb%mOrb, this%iCellVec, this%cellVec,&
              & iSquare, iSparseStart, img2CentCell)
        end if
      end do
    end if
    call ud2qm(rhoPrim)

    tDFTBU = .false.
    if (size(UJ) > 0) then
      tDFTBU = .true.
    end if

    TS = 0.0_dp
    call calcEnergies(this%sccCalc, qq, q0, chargePerShell, this%speciesAll, this%tLaser, .false.,&
        & tDFTBU, tDualSpinOrbit, rhoPrim, ham0, orb, neighbourList, nNeighbourSK, img2CentCell,&
        & iSparseStart, 0.0_dp, 0.0_dp, TS, potential, energy, thirdOrd, solvation, rangeSep, reks,&
        & qDepExtPot, qBlock, qiBlock, nDftbUFunc, UJ, nUJ, iUJ, niUJ, xi, iAtInCentralRegion,&
        & tFixEf, Ef, onSiteElements)
    call sumEnergies(energy)
    ! calcEnergies then sumEnergies returns the total energy Etotal including repulsive and
    ! dispersions energies

    ! Calculate nuclear kinetic energy
    energyKin = 0.0_dp
    if (this%tIons) then
      energyKin = 0.5_dp * sum(this%movedMass(:,:) * this%movedVelo(:,:)**2)
      energy%Etotal = energy%Etotal + energyKin
    end if

  end subroutine getTDEnergy


  !> Create all necessary matrices and instances for dynamics
  subroutine initializeTDVariables(this, rho, H1, Ssqr, Sinv, H0, ham0, over, ham, eigvecsReal,&
      & filling, orb, rhoPrim, potential, iNeighbour, nNeighbourSK, iSquare, iSparseStart,&
      & img2CentCell, Eiginv, EiginvAdj, energy, ErhoPrim, skOverCont, qBlock, qNetAtom, UJ,&
      & onSiteElements, eigvecsCplx, H1LC, bondWork, fdBondEnergy, fdBondPopul, lastBondPopul, time)

    !> ElecDynamics instance
    type(TElecDynamics), intent(inout) :: this

    !> Real Eigenvectors
    real(dp), intent(inout), allocatable :: eigvecsReal(:,:,:)

    !> Complex Eigevenctors
    complex(dp), intent(inout), allocatable :: eigvecsCplx(:,:,:)

    !> overlap (sparse)
    real(dp), allocatable, intent(in) :: over(:)

    !> resulting hamiltonian (sparse)
    real(dp), allocatable, intent(in) :: ham(:,:)

    !> occupations
    real(dp), intent(inout) :: filling(:,:,:)

    !> Atomic neighbour data
    integer, intent(in) :: iNeighbour(0:,:)

    !> Number of neighbours for each of the atoms
    integer, intent(in) :: nNeighbourSK(:)

    !> Index array for start of atomic block in dense matrices
    integer, intent(in) :: iSquare(:)

    !> index array for location of atomic blocks in large sparse arrays
    integer, intent(in) :: iSparseStart(0:,:)

    !> image atoms to their equivalent in the central cell
    integer, intent(in) :: img2CentCell(:)

    !> Atomic orbital information
    type(TOrbitals), intent(in) :: orb

    !> potential acting on the system
    type(TPotentials), intent(out) :: potential

    !> data type for energy components and total
    type(TEnergies), intent(out) :: energy

    !> sparse density matrix
    real(dp), allocatable, intent(out) :: rhoPrim(:,:)

    !> Square overlap matrix
    complex(dp), intent(inout) :: Ssqr(:,:,:)

    !> Square overlap inverse
    complex(dp), intent(inout) :: Sinv(:,:,:)

    !> Square hamiltonian
    complex(dp), intent(out) :: H1(:,:,:)

    !> Density matrix
    complex(dp), intent(inout) :: rho(:,:,:)

    !> Inverse of eigenvectors matrix (for populations)
    complex(dp), allocatable, intent(out) :: Eiginv(:,:,:)

    !> Adjoint of the inverse of eigenvectors matrix (for populations)
    complex(dp), allocatable, intent(out) :: EiginvAdj(:,:,:)

    !> Non-SCC hamiltonian
    real(dp), intent(in) :: H0(:)

    !> Local sparse storage for non-SCC hamiltonian
    real(dp), allocatable, intent(out) :: ham0(:)

    !> Raw overlap data
    type(TSlakoCont), intent(in) :: skOverCont

    !> Energy weighted density matrix
    real(dp), allocatable, intent(out) :: ErhoPrim(:)

    !> block (dual) atomic populations
    real(dp), intent(inout), allocatable :: qBlock(:,:,:,:)

    !> net (onsite only) atomic charges
    real(dp), intent(inout), allocatable :: qNetAtom(:)

    !> U-J prefactors in DFTB+U
    real(dp), intent(in), allocatable :: UJ(:,:)

    !> Corrections terms for on-site elements
    real(dp), intent(in), allocatable :: onSiteElements(:,:,:,:)

    !> LC contribution to hamiltonian
    complex(dp), allocatable, intent(inout) :: H1LC(:,:)

    !> Container for either bond populations or bond energy
    real(dp), allocatable, intent(inout) :: bondWork(:, :)

    !> Last calculated bond population (for tagged output)
    real(dp), intent(inout) :: lastBondPopul

    !> Pairwise bond population output file ID
    integer, intent(inout) :: fdBondPopul

    !> Pairwise bond energy output file ID
    integer, intent(inout) :: fdBondEnergy

    !> simulation time (in atomic units)
    real(dp), intent(in) :: time

    real(dp), allocatable :: T2(:,:), T3(:,:)
    complex(dp), allocatable :: T4(:,:)
    integer :: iSpin, iOrb, iOrb2, fillingsIn, iKS, iK

    allocate(rhoPrim(size(ham, dim=1), this%nSpin))
    allocate(ErhoPrim(size(ham, dim=1)))
    this%nSparse = size(H0)
    allocate(ham0(size(H0)))
    ham0(:) = H0

    if (this%tRealHS) then
      allocate(T2(this%nOrbs,this%nOrbs))
      allocate(T3(this%nOrbs, this%nOrbs))
    else
      allocate(T4(this%nOrbs,this%nOrbs))
    end if

    if (.not. this%tReadRestart) then
      Ssqr(:,:,:) = 0.0_dp
      Sinv(:,:,:) = 0.0_dp
      do iKS = 1, this%parallelKS%nLocalKS
        if (this%tRealHS) then
          call unpackHS(T2, over, iNeighbour, nNeighbourSK, iSquare, iSparseStart, img2CentCell)
          call blockSymmetrizeHS(T2, iSquare)
          Ssqr(:,:,iKS) = cmplx(T2, 0, dp)
          T3(:,:) = 0.0_dp
          do iOrb = 1, this%nOrbs
            T3(iOrb, iOrb) = 1.0_dp
          end do
          call gesv(T2, T3)
          Sinv(:,:,iKS) = cmplx(T3, 0, dp)
        else
          iK = this%parallelKS%localKS(1, iKS)
          iSpin = this%parallelKS%localKS(2, iKS)
          T4(:,:) = cmplx(0,0,dp)
          call unpackHS(T4, over, this%kPoint(:,iK), iNeighbour, nNeighbourSK, this%iCellVec,&
              & this%cellVec, iSquare, iSparseStart, img2CentCell)
          call blockHermitianHS(T4, iSquare)
          Ssqr(:,:,iKS) = T4
          Sinv(:,:,iKS) = cmplx(0,0,dp)
          do iOrb = 1, this%nOrbs
            Sinv(iOrb, iOrb, iKS) = 1.0_dp
          end do
          call gesv(T4, Sinv(:,:,iKS))
        end if
      end do
      write(stdOut,"(A)")'S inverted'

      do iKS = 1, this%parallelKS%nLocalKS
        iK = this%parallelKS%localKS(1, iKS)
        iSpin = this%parallelKS%localKS(2, iKS)
        if (this%tRealHS) then
          call unpackHS(T3, ham(:,iSpin), iNeighbour, nNeighbourSK, iSquare, iSparseStart,&
              & img2CentCell)
          call blockSymmetrizeHS(T3, iSquare)
          H1(:,:,iKS) = cmplx(T3, 0, dp)
        else
          call unpackHS(H1(:,:,iKS), ham(:,iSpin), this%kPoint(:,iK), iNeighbour, nNeighbourSK,&
              & this%iCellVec, this%cellVec, iSquare, iSparseStart, img2CentCell)
          call blockHermitianHS(H1(:,:,iKS), iSquare)
        end if
      end do
    end if

    if (this%tPopulations) then
      allocate(Eiginv(this%nOrbs, this%nOrbs, this%parallelKS%nLocalKS))
      allocate(EiginvAdj(this%nOrbs, this%nOrbs, this%parallelKS%nLocalKS))
      do iKS = 1, this%parallelKS%nLocalKS
        if (this%tRealHS) then
          call tdPopulInit(this, Eiginv(:,:,iKS), EiginvAdj(:,:,iKS), eigvecsReal(:,:,iKS))
        else
          call tdPopulInit(this, Eiginv(:,:,iKS), EiginvAdj(:,:,iKS), &
              & eigvecsCplx=eigvecsCplx(:,:,iKS))
        end if
      end do
    end if

    if (this%tFillingsFromFile) then
      filling(:,:,:) = 0.0_dp
      open(newunit=fillingsIn, file='fillings.in')
      do iSpin=1,this%nSpin
        do iOrb=1,this%nOrbs
          read(fillingsIn, *) filling(iOrb,1,iSpin)
        end do
      end do
    end if

    if (.not.this%tReadRestart) then
      rho(:,:,:) = 0.0_dp
      do iKS = 1, this%parallelKS%nLocalKS
        iK = this%parallelKS%localKS(1, iKS)
        iSpin = this%parallelKS%localKS(2, iKS)
        if (this%tRealHS) then
          T2(:,:) = 0.0_dp
          call makeDensityMatrix(T2, eigvecsReal(:,:,iKS), filling(:,1,iSpin))
          rho(:,:,iKS) = cmplx(T2, 0, dp)
        else
          call makeDensityMatrix(rho(:,:,iKS), eigvecsCplx(:,:,iKS), filling(:,iK,iSpin))
        end if
        do iOrb = 1, this%nOrbs-1
          do iOrb2 = iOrb+1, this%nOrbs
            rho(iOrb, iOrb2, iKS) = conjg(rho(iOrb2, iOrb, iKS))
          end do
        end do
      end do
    end if

    call init(potential, orb, this%nAtom, this%nSpin)
    call TEnergies_init(energy, this%nAtom)

    if ((size(UJ) /= 0) .or. allocated(onSiteElements)) then
      allocate(qBlock(orb%mOrb, orb%mOrb, this%nAtom, this%nSpin))
    end if

    if (this%tNetCharges) then
      allocate(qNetAtom(this%nAtom))
    end if

    if (this%isRangeSep) then
      allocate(H1LC(this%nOrbs, this%nOrbs))
    end if

    if (this%tBondE .or. this%tBondP) then
      allocate(bondWork(this%nAtom, this%nAtom))
    end if
    if (this%tBondE) then
      call openFile(this, fdBondEnergy, 'bondenergy.bin', isBinary = .true.)
    end if
    if (this%tBondP) then
      call openFile(this, fdBondPopul, 'bondpop.bin', isBinary = .true.)
    end if
    call getBondPopulAndEnergy(this, bondWork, lastBondPopul, rhoPrim, ham0, over, iNeighbour,&
        & nNeighbourSK, iSparseStart, img2CentCell, iSquare, fdBondEnergy, fdBondPopul, .false.,&
        & time)

  end subroutine initializeTDVariables


  !> Perfoms a step backwards to boot the dynamics using the Euler algorithm.
  !> Output is rho(deltaT) called rhoNew, input is rho(t=0) (ground state) called rho
  subroutine initializePropagator(this, step, rho, rhoNew, H1, Sinv, coordAll, skOverCont,&
      & orb, neighbourList, nNeighbourSK, img2CentCell, iSquare, rangeSep)

    !> ElecDynamics instance
    type(TElecDynamics), intent(inout) :: this

    !> Density matrix at next step
    complex(dp), intent(inout) :: rhoNew(:,:,:)

    !> Square overlap inverse
    complex(dp), intent(in) :: Sinv(:,:,:)

    !> Square hamiltonian
    complex(dp), intent(inout) :: H1(:,:,:)

    !> Density matrix
    complex(dp), intent(in) :: rho(:,:,:)

    !> Time step in atomic units (with sign, to perform step backwards or forwards)
    real(dp), intent(in) :: step

    !> Coords of the atoms (3, nAllAtom)
    real(dp), intent(in) :: coordAll(:,:)

    !> Raw overlap data
    type(TSlakoCont), intent(in) :: skOverCont

    !> data type for atomic orbital information
    type(TOrbitals), intent(in) :: orb

    !> ADT for neighbour parameters
    type(TNeighbourList), intent(in) :: neighbourList

    !> nr. of neighbours for atoms out to max interaction distance (excluding Ewald terms)
    integer, intent(in) :: nNeighbourSK(:)

    !> Index array for start of atomic block in dense matrices
    integer, intent(in) :: iSquare(:)

    !> image atoms to their equivalent in the central cell
    integer, intent(in) :: img2CentCell(:)

    !> Range separation contributions
    type(TRangeSepFunc), allocatable, intent(inout) :: rangeSep

    integer :: iKS
    complex(dp), allocatable :: RdotSprime(:,:)

    allocate(RdotSprime(this%nOrbs,this%nOrbs))

    if (this%tIons) then
      call getRdotSprime(this, RdotSprime, coordAll, skOverCont, orb, img2CentCell, &
          &neighbourList, nNeighbourSK, iSquare)
    else
      RdotSprime(:,:) = 0.0_dp
    end if

    do iKS = 1, this%parallelKS%nLocalKS
      if (this%tIons .or. (.not. this%tRealHS) .or. allocated(rangeSep)) then
        H1(:,:,iKS) = RdotSprime + imag * H1(:,:,iKS)
        call propagateRho(this, rhoNew(:,:,iKS), rho(:,:,iKS), H1(:,:,iKS), Sinv(:,:,iKS), step)
      else
        ! The following line is commented to make the fast propagate work since it needs a real H
        !H1(:,:,iKS) = imag * H1(:,:,iKS)
        call propagateRhoRealH(this, rhoNew(:,:,iKS), rho(:,:,iKS), H1(:,:,iKS), Sinv(:,:,iKS),&
            & step)
      end if
    end do

  end subroutine initializePropagator


  !> Propagate rho, notice that H = iH (coeficients are real)
  subroutine propagateRho(this, rhoOld, rho, H1, Sinv, step)

    !> ElecDynamics instance
    type(TElecDynamics), intent(inout) :: this

    !> Density matrix at previous step
    complex(dp), intent(inout) :: rhoOld(:,:)

    !> Density matrix
    complex(dp), intent(in) :: rho(:,:)

    !> Square imaginary hamiltonian plus non-adiabatic contribution
    complex(dp), intent(in) :: H1(:,:)

    !> Square overlap inverse
    complex(dp), intent(in) :: Sinv(:,:)

    !> Time step in atomic units
    real(dp), intent(in) :: step

    complex(dp), allocatable :: T1(:,:)
    allocate(T1(this%nOrbs,this%nOrbs))

    T1(:,:) = 0.0_dp
    call gemm(T1, Sinv, H1)
    call gemm(rhoOld, T1, rho, cmplx(-step, 0, dp), cmplx(1, 0, dp))
    call gemm(rhoOld, rho, T1, cmplx(-step, 0, dp), cmplx(1, 0, dp), 'N', 'C')

  end subroutine propagateRho


  !> Propagate rho for real Hamiltonian (used for frozen nuclei dynamics and gamma point periodic)
  subroutine propagateRhoRealH(this, rhoOld, rho, H1, Sinv, step)

    !> ElecDynamics instance
    type(TElecDynamics), intent(inout) :: this

    !> Density matrix at previous step
    complex(dp), intent(inout) :: rhoOld(:,:)

    !> Density matrix
    complex(dp), intent(in) :: rho(:,:)

    !> Square hamiltonian
    complex(dp), intent(in) :: H1(:,:)

    !> Square overlap inverse
    complex(dp), intent(in) :: Sinv(:,:)

    !> Time step in atomic units
    real(dp), intent(in) :: step

    real(dp), allocatable :: T1R(:,:),T2R(:,:)
    real(dp), allocatable :: T3R(:,:),T4R(:,:),T5R(:,:)
    integer :: i,j

    allocate(T1R(this%nOrbs,this%nOrbs))
    allocate(T2R(this%nOrbs,this%nOrbs))
    allocate(T3R(this%nOrbs,this%nOrbs))
    allocate(T4R(this%nOrbs,this%nOrbs))
    allocate(T5R(this%nOrbs,this%nOrbs))


    ! The code below takes into account that Sinv and H1 are real, this is twice as fast as the
    ! original above (propageteRho)

    ! get the real part of Sinv and H1
    T1R(:,:) = real(H1)
    T2R(:,:) = real(Sinv)
    call gemm(T3R,T2R,T1R)

    ! calculate the first term products for the real and imaginary parts independently
    T1R(:,:) = real(rho)
    T2R(:,:) = aimag(rho)
    call gemm(T4R,T3R,T1R)
    call gemm(T5R,T3R,T2R)

    ! build the commutator combining the real and imaginary parts of the previous result
    !$omp parallel do private(i,j)
    do i=1,this%nOrbs
      do j=1,this%nOrbs
        rhoOld(i,j) = rhoOld(i,j) + cmplx(0, -step, dp) * (T4R(i,j) + imag * T5R(i,j)) &
            + cmplx(0, step, dp) * conjg(T4R(j,i) + imag * T5R(j,i))
      enddo
    enddo
    !$omp end parallel do

  end subroutine propagateRhoRealH


  !> Initialize output files
  subroutine initTDOutput(this, dipoleDat, qDat, energyDat, populDat, forceDat, coorDat)
    !> ElecDynamics instance
    type(TElecDynamics), intent(in) :: this

    !> Dipole output file ID
    integer, intent(out) :: dipoleDat

    !> Charge output file ID
    integer, intent(out) :: qDat

    !> Energy output file ID
    integer, intent(out) :: energyDat

    !> Populations  output file ID
    integer, intent(out) :: populDat(:)

    !> Forces output file ID
    integer, intent(out) :: forceDat

    !> Coords  output file ID
    integer, intent(out) :: coorDat

    character(20) :: dipoleFileName
    character(1) :: strSpin
    character(3) :: strK
    integer :: iSpin, iKS, iK

    if (this%tKick) then
      if (this%currPolDir == 1) then
        dipoleFileName = 'mux.dat'
      else if (this%currPolDir == 2) then
        dipoleFileName = 'muy.dat'
      else if (this%currPolDir == 3) then
        dipoleFileName = 'muz.dat'
      end if
    else
      dipoleFileName = 'mu.dat'
    end if
    call openFile(this, dipoleDat, dipoleFileName)

    write(dipoleDat, "(A)", advance = "NO")"#           time (fs)    |"
    select case(this%nSpin)
    case(1)
      write(dipoleDat, "(A)", advance = "NO")"     mu_x (e.angstrom)   |"
      write(dipoleDat, "(A)", advance = "NO")"     mu_y (e.angstrom)   |"
      write(dipoleDat, "(A)", advance = "NO")"     mu_z (e.angstrom)   |"
    case(2)
      write(dipoleDat, "(A)", advance = "NO")"  mu_x (up) (e.angstrom) |"
      write(dipoleDat, "(A)", advance = "NO")"  mu_y (up) (e.angstrom) |"
      write(dipoleDat, "(A)", advance = "NO")"  mu_z (up) (e.angstrom) |"
      write(dipoleDat, "(A)", advance = "NO")" mu_x (down) (e.angstrom)|"
      write(dipoleDat, "(A)", advance = "NO")" mu_y (down) (e.angstrom)|"
      write(dipoleDat, "(A)", advance = "NO")" mu_z (down) (e.angstrom)|"
    end select
    write(dipoleDat, "(A)")

    if (this%tdWriteExtras) then
      call openFile(this, qDat, 'qsvst.dat')
      write(qDat, "(A)", advance = "NO")"#             time (fs)      |"
      write(qDat, "(A)", advance = "NO")"   total net charge (e)  |"
      write(qDat, "(A)", advance = "NO")"   charge (atom_1) (e)   |"
      write(qDat, "(A)", advance = "NO")"   charge (atom_2) (e)   |"
      write(qDat, "(A)", advance = "NO")"        ...        |"
      write(qDat, "(A)", advance = "NO")"   charge (atom_N) (e)   |"
      write(qDat, "(A)")

      call openFile(this, energyDat, 'energyvst.dat')
      write(energyDat, "(A)", advance = "NO")"#                  time (fs)         |"
      write(energyDat, "(A)", advance = "NO")"        E total (H)         |"
      write(energyDat, "(A)", advance = "NO")"        E non-SCC (H)       |"
      write(energyDat, "(A)", advance = "NO")"            E SCC (H)       |"
      write(energyDat, "(A)", advance = "NO")"           E spin (H)       |"
      write(energyDat, "(A)", advance = "NO")"       E external (H)       |"
      write(energyDat, "(A)", advance = "NO")"            E rep (H)       |"
      write(energyDat, "(A)", advance = "NO")"E kinetic nuclear (H)       |"
      write(energyDat, "(A)", advance = "NO")"     E dispersion (H)       |"
      write(energyDat, "(A)")

      if (this%tForces) then
        call openFile(this, forceDat, 'forcesvst.dat')
        write(forceDat, "(A)", advance = "NO")"#           time (fs)       |"
        write(forceDat, "(A)", advance = "NO")" force (atom_1) (H/b)   |  force (atom_2) (H/b)  |"
        write(forceDat, "(A)", advance = "NO")"           ...          |  force (atom_N) (H/b)  |"
        write(forceDat, "(A)")
      end if

      if (this%tIons) then
        call openFile(this, coorDat, 'tdcoords.xyz')
      end if
    end if

    if (this%tPopulations) then
      do iKS = 1, this%parallelKS%nLocalKS
        iSpin = this%parallelKS%localKS(2, iKS)
        write(strSpin,'(i1)')iSpin
        if (this%tRealHS) then
          call openFile(this, populDat(iKS), 'molpopul' // trim(strSpin) // '.dat')
          write(populDat(iKS), "(A,A)") "#  GS molecular orbital populations, spin channel : ",&
              & trim(strSpin)
        else
          iK = this%parallelKS%localKS(1, iKS)
          write(strK,'(i0.3)')iK
          call openFile(this, populDat(iKS), 'molpopul' // trim(strSpin) // '-' // trim(strK) //&
              & '.dat')
          write(populDat(iKS), "(A,A,A,A)") "#  GS molecular orbital populations, spin channel : ",&
              & trim(strSpin), ", k-point number: ", trim(strK)
        end if
        write(populDat(iKS), "(A)", advance = "NO")"#          time (fs)            |"
        write(populDat(iKS), "(A)", advance = "NO")"   population (orb 1)       |"
        write(populDat(iKS), "(A)", advance = "NO")"    population (orb 2)      |"
        write(populDat(iKS), "(A)", advance = "NO")"           ...              |"
        write(populDat(iKS), "(A)", advance = "NO")"    population (orb N)      |"
        write(populDat(iKS), "(A)")
      end do
    end if

  end subroutine initTDOutput


  !> Close output files
  subroutine closeTDOutputs(this, dipoleDat, qDat, energyDat, populDat, forceDat, coorDat,&
      & fdBondPopul, fdBondEnergy)

    !> ElecDynamics instance
    type(TElecDynamics), intent(in) :: this

    !> Dipole output file ID
    integer, intent(in) :: dipoleDat

    !> Charge output file ID
    integer, intent(in) :: qDat

    !> Energy output file ID
    integer, intent(in) :: energyDat

    !> Populations output file ID
    integer, intent(in) :: populDat(:)

    !> Forces output file ID
    integer, intent(in) :: forceDat

    !> Coords output file ID
    integer, intent(in) :: coorDat

    !> Pairwise bond population output file ID
    integer, intent(in) :: fdBondPopul

    !> Pairwise bond energy output file ID
    integer, intent(in) :: fdBondEnergy

    integer :: iKS

    close(dipoleDat)

    if (this%tdWriteExtras) then
      close(qDat)
      close(energyDat)

      if (this%tPopulations) then
        do iKS = 1, this%parallelKS%nLocalKS
          close(populDat(iKS))
        end do
      end if

      if (this%tIons) then
        close(coorDat)
      end if

      if (this%tForces) then
        close(forceDat)
      end if

    end if

    if (this%tBondP) then
      close(fdBondPopul)
    end if

    if (this%tBondE) then
      close(fdBondEnergy)
    end if

  end subroutine closeTDOutputs


  !> Open files in different ways depending on their previous existance
  subroutine openFile(this, unitName, fileName, isBinary)

    !> ElecDynamics instance
    type(TElecDynamics), intent(in) :: this

    !> File ID
    integer, intent(out) :: unitName

    !> Name of the file to open
    character(*), intent(in) :: fileName

    !> should this be a binary file?
    logical, intent(in), optional :: isBinary

    character(lc) :: newName

    character(lc) :: strCount

    logical :: exist=.false.

    integer :: iCount

    logical :: isBinary_

    if (present(isBinary)) then
      isBinary_ = isBinary
    else
      isBinary_ = .false.
    end if

    newName = fileName
    ! changed the append by this block to rename the restarted output
    if (this%tReadRestart) then
      inquire(file=fileName, exist=exist)
      iCount = 1
      do while (exist)
        write(strCount,'(I0)') iCount
        newName = "rest" // trim(strCount) // "_" // fileName
        inquire(file=newName, exist=exist)
        iCount = iCount + 1
      end do
    end if

    if (isBinary_) then
      open(newunit=unitName, file=newName, form='unformatted', access='stream')
    else
      open(newunit=unitName, file=newName, action="write")
    end if

  end subroutine openFile


  !> Write to restart file
  subroutine writeRestartFile(rho, rhoOld, coord, veloc, time, dt, fileName, tAsciiFile)

    !> Density matrix
    complex(dp), intent(in) :: rho(:,:,:)

    !> Density matrix at previous time step
    complex(dp), intent(in) :: rhoOld(:,:,:)

    !> atomic coordinates
    real(dp), intent(in) :: coord(:,:)

    !> atomic velocities
    real(dp), intent(in) :: veloc(:,:)

    !> simulation time (in atomic units)
    real(dp), intent(in) :: time

    !> time step being used (in atomic units)
    real(dp), intent(in) :: dt

    !> name of the dump file
    character(len=*), intent(in) :: fileName

    !> Should restart data be written as ascii (cross platform, but potentially lower
    !> reproducibility) or binary files
    logical, intent(in) :: tAsciiFile

    integer :: fd, ii, jj, kk, iErr

    if (tAsciiFile) then
      open(newunit=fd, file=trim(fileName) // '.dat', position="rewind", status="replace",&
          & iostat=iErr)
    else
      open(newunit=fd, file=trim(fileName) // '.bin', form='unformatted', access='stream',&
          & action='write', iostat=iErr)
    end if

    if (iErr /= 0) then
      if (tAsciiFile) then
        write(error_string, "(A,A,A)") "Failure to open external restart file ",trim(fileName),&
            & ".dat for writing"
      else
        write(error_string, "(A,A,A)") "Failure to open external restart file ",trim(fileName),&
            & ".bin for writing"
      end if
      call error(error_string)
    end if

    if (tAsciiFile) then

      write(fd, *)tdDumpFormat
      write(fd, *)size(rho, dim=1), size(rho, dim=3), size(coord, dim=2), time, dt
      do ii = 1, size(rho, dim=3)
        do jj = 1, size(rho, dim=2)
          do kk = 1, size(rho, dim=1)
            write(fd, *)rho(kk,jj,ii)
          end do
        end do
      end do
      do ii = 1, size(rhoOld, dim=3)
        do jj = 1, size(rhoOld, dim=2)
          do kk = 1, size(rhoOld, dim=1)
            write(fd, *)rhoOld(kk,jj,ii)
          end do
        end do
      end do
      do ii = 1, size(coord, dim=2)
        write(fd, *)coord(:,ii)
      end do
      do ii = 1, size(veloc, dim=2)
        write(fd, *)veloc(:,ii)
      end do

    else

      write(fd)tdDumpFormat
      write(fd)size(rho, dim=1), size(rho, dim=3), size(coord, dim=2), time, dt
      write(fd) rho, rhoOld, coord, veloc

    end if

    close(fd)

  end subroutine writeRestartFile


  !> read a restart file containing density matrix, overlap, coordinates and time step
  subroutine readRestartFile(rho, rhoOld, coord, veloc, time, dt, fileName, tAsciiFile)

    !> Density Matrix
    complex(dp), intent(out) :: rho(:,:,:)

    !> Previous density Matrix
    complex(dp), intent(out) :: rhoOld(:,:,:)

    !> atomic coordinates
    real(dp), intent(out) :: coord(:,:)

    !> Previous simulation elapsed time until restart file writing
    real(dp), intent(out) :: time

    !> time step being currently used (in atomic units) for checking compatibility
    real(dp), intent(in) :: dt

    !> Name of the file to open
    character(*), intent(in) :: fileName

    !> atomic velocities
    real(dp), intent(out) :: veloc(:,:)

    !> Should restart data be read as ascii (cross platform, but potentially lower reproducibility)
    !> or binary files
    logical, intent(in) :: tAsciiFile

    integer :: fd, ii, jj, kk, nOrb, nSpin, nAtom, version, iErr
    real(dp) :: deltaT
    logical :: tExist

    if (tAsciiFile) then
      inquire(file=trim(fileName)//'.dat', exist=tExist)
      if (.not. tExist) then
        call error("TD restart file " // trim(fileName)//'.dat' // " is missing")
      end if
    else
      inquire(file=trim(fileName)//'.bin', exist=tExist)
      if (.not. tExist) then
        call error("TD restart file " // trim(fileName)//'.bin' // " is missing")
      end if
    end if

    if (tAsciiFile) then
      open(newunit=fd, file=trim(fileName)//'.dat', status='old', action='READ', iostat=iErr)
    else
      open(newunit=fd, file=trim(fileName)//'.bin', form='unformatted', access='stream',&
          & action='read', iostat=iErr)
    end if

    if (iErr /= 0) then
      if (tAsciiFile) then
        write(error_string, "(A,A,A)") "Failure to open external tddump file",trim(fileName), ".dat"
      else
        write(error_string, "(A,A,A)") "Failure to open external tddump file",trim(fileName), ".bin"
      end if
      call error(error_string)
    end if
    rewind(fd)

    if (tAsciiFile) then
      read(fd, *)version
      if (version /= tdDumpFormat) then
        call error("Unknown TD format")
      end if
      read(fd, *) nOrb, nSpin, nAtom, time, deltaT
      if (nOrb /= size(rho, dim=1)) then
        write(error_string, "(A,I0,A,I0)")"Incorrect number of orbitals, ",nOrb,&
            & " in tddump file, should be ",size(rho, dim=1)
        call error(error_string)
      end if
      if (nSpin /= size(rho, dim=3)) then
        write(error_string, "(A,I1,A,I1)")"Incorrect number of spin channels, ",nSpin,&
            & " in tddump file, should be ",size(rho, dim=3)
        call error(error_string)
      end if
      if (nAtom /= size(coord, dim=2)) then
        write(error_string, "(A,I0,A,I0)")"Incorrect number of atoms, ",nAtom,&
            & " in tddump file, should be ", size(coord, dim=2)
        call error(error_string)
      end if
      if (abs(deltaT - dt) > epsilon(0.0_dp)) then
        write(error_string, "(A,E14.8,A,E14.8)")"Restart file generated for time step",&
            & deltaT, " instead of current timestep of", dt
      end if
      do ii = 1, size(rho, dim=3)
        do jj = 1, size(rho, dim=2)
          do kk = 1, size(rho, dim=1)
            read(fd, *)rho(kk,jj,ii)
          end do
        end do
      end do
      do ii = 1, size(rhoOld, dim=3)
        do jj = 1, size(rhoOld, dim=2)
          do kk = 1, size(rhoOld, dim=1)
            read(fd, *)rhoOld(kk,jj,ii)
          end do
        end do
      end do
      do ii = 1, size(coord, dim=2)
        read(fd, *)coord(:,ii)
      end do
      do ii = 1, size(veloc, dim=2)
        read(fd, *)veloc(:,ii)
      end do
    else
      read(fd)version
      if (version /= tdDumpFormat) then
        call error("Unknown TD format")
      end if
      read(fd) nOrb, nSpin, nAtom, time, deltaT
      if (nOrb /= size(rho, dim=1)) then
        write(error_string, "(A,I0,A,I0)")"Incorrect number of orbitals, ",nOrb,&
            & " in tddump file, should be ",size(rho, dim=1)
        call error(error_string)
      end if
      if (nSpin /= size(rho, dim=3)) then
        write(error_string, "(A,I1,A,I1)")"Incorrect number of spin channels, ",nSpin,&
            & " in tddump file, should be ",size(rho, dim=3)
        call error(error_string)
      end if
      if (nAtom /= size(coord, dim=2)) then
        write(error_string, "(A,I0,A,I0)")"Incorrect number of atoms, ",nAtom,&
            & " in tddump file, should be ", size(coord, dim=2)
        call error(error_string)
      end if
      if (abs(deltaT - dt) > epsilon(0.0_dp)) then
        write(error_string, "(A,E14.8,A,E14.8)")"Restart file generated for time step",&
            & deltaT, " instead of current timestep of", dt
      end if
      read(fd) rho, rhoOld, coord, veloc
    end if
    close(fd)

  end subroutine readRestartFile


  !> Write results to file
  subroutine writeTDOutputs(this, dipoleDat, qDat, energyDat, forceDat, coorDat, fdBondPopul,&
      & fdBondEnergy, time, energy, energyKin, dipole, deltaQ, coord, totalForce, iStep)

    !> ElecDynamics instance
    type(TElecDynamics), intent(in) :: this

    !> data type for energy components and total
    type(TEnergies), intent(in) :: energy

    !> Dipole output file ID
    integer, intent(in) :: dipoleDat

    !> Charge output file ID
    integer, intent(in) :: qDat

    !> Energy output file ID
    integer, intent(in) :: energyDat

    !> Elapsed simulation time
    real(dp), intent(in) :: time

    !> Dipole moment
    real(dp), intent(in) :: dipole(:,:)

    !> Negative gross charge
    real(dp), intent(in) :: deltaQ(:,:)

    !> current step of the propagation
    integer, intent(in) :: iStep

    !> Forces output file ID
    integer, intent(in) :: forceDat

    !> Coords output file ID
    integer, intent(in) :: coorDat

    !> Pairwise bond population output file ID
    integer, intent(in) :: fdBondPopul

    !> Pairwise bond energy output file ID
    integer, intent(in) :: fdBondEnergy

    !> atomic coordinates
    real(dp), intent(in) :: coord(:,:)

    !> Kinetic energy
    real(dp), intent(in) :: energyKin

    !> forces (3, nAtom)
    real(dp), intent(in) :: totalForce(:,:)

    real(dp) :: auxVeloc(3, this%nAtom)
    integer :: iAtom, iSpin, iDir

    write(dipoleDat, '(7F25.15)') time * au__fs, ((dipole(iDir, iSpin) * Bohr__AA, iDir=1, 3),&
        & iSpin=1, this%nSpin)

    if (this%tdWriteExtras) then
      write(energydat, '(9F30.15)') time * au__fs, energy%Etotal, energy%EnonSCC, energy%eSCC,&
          & energy%Espin, energy%Eext, energy%Erep, energyKin, energy%eDisp
    end if

    if (mod(iStep, this%writeFreq) == 0) then
      if (this%tdWriteExtras) then
        write(qDat, "(2X,2F25.15)", advance="no") time * au__fs, -sum(deltaQ)
        do iAtom = 1, this%nAtom
          write(qDat, "(F25.15)", advance="no")-sum(deltaQ(iAtom,:))
        end do
        write(qDat,*)
      end if
    end if

    if (this%tIons .and. (mod(iStep,this%writeFreq) == 0)) then
      if (this%tdWriteExtras) then
        auxVeloc = 0.0_dp
        auxVeloc(:, this%indMovedAtom) = this%movedVelo
        write(coorDat,'(I5)')this%nAtom
        write(coorDat,*) 'MD step:', iStep, 'time', time * au__fs
        do iAtom=1,this%nAtom
          write(coorDat, '(A2, 6F16.8)') trim(this%speciesName(this%species(iAtom))), &
              &coord(:, iAtom) * Bohr__AA, auxVeloc(:, iAtom) * Bohr__AA / au__fs
        end do
      endif
    end if

    if (this%tForces .and. (mod(iStep,this%writeFreq) == 0)) then
      if (this%tdWriteExtras) then
        write(forceDat, "(F25.15)", advance="no") time * au__fs
        do iAtom = 1, this%nAtom
          write(forceDat, "(3F25.15)", advance="no") totalForce(:,iAtom)
        end do
        write(forceDat,*)
      end if
    end if

    ! Flush output every 5% of the simulation
    if (mod(iStep, max(this%nSteps / 20, 1)) == 0 .and. iStep > this%writeFreq) then
      if (this%tdWriteExtras) then
        flush(qDat)
        flush(energyDat)
        if (this%tIons) then
          flush(coorDat)
        end if
        if (this%tForces) then
          flush(forceDat)
        end if
        if (this%tBondP) then
          flush(fdBondPopul)
        end if
        if (this%tBondE) then
          flush(fdBondEnergy)
        end if
      end if
    end if

  end subroutine writeTDOutputs


  !> Initialize matrices for populations
  !> Note, this will need to get generalised for complex eigenvectors
  subroutine tdPopulInit(this, Eiginv, EiginvAdj, eigvecsReal, eigvecsCplx)

    !> ElecDynamics instance
    type(TElecDynamics), intent(in) :: this

    !> Inverse of eigenvectors matrix (for populations)
    complex(dp), intent(out) :: Eiginv(:,:)

    !> Adjoint of the inverse of eigenvectors matrix (for populations)
    complex(dp), intent(out) :: EiginvAdj(:,:)

    !> Eigenvectors
    real(dp), intent(in), optional :: eigvecsReal(:,:)

    !> Complex Eigevenctors
    complex(dp), intent(in), optional :: eigvecsCplx(:,:)

    complex(dp), allocatable :: T2(:,:), T3(:,:)
    integer :: iOrb

    allocate(T2(this%nOrbs, this%nOrbs), T3(this%nOrbs, this%nOrbs))

    if (this%tRealHS) then
      T2 = cmplx(eigvecsReal, 0, dp)
    else
      T2 = eigvecsCplx
    end if

    T3 = 0.0_dp
    do iOrb = 1, this%nOrbs
      T3(iOrb, iOrb) = 1.0_dp
    end do
    call gesv(T2,T3)
    Eiginv(:,:) = T3

    if (this%tRealHS) then
      T2 = cmplx(transpose(eigvecsReal), 0, dp)
    else
      T2 = conjg(transpose(eigvecsCplx))
    end if

    T3 = 0.0_dp
    do iOrb = 1, this%nOrbs
      T3(iOrb, iOrb) = 1.0_dp
    end do
    call gesv(T2,T3)
    EiginvAdj(:,:) = T3

    deallocate(T2, T3)

  end subroutine tdPopulInit


  ! updates Eiginv and EiginvAdj if nuclear dynamics is done
  ! important to call after H1 has been updated with new charges and before D is included in H1
  subroutine updateBasisMatrices(this, electronicSolver, Eiginv, EiginvAdj, H1, Ssqr)
    !> ElecDynamics instance
    type(TElecDynamics), intent(in) :: this

    !> Electronic solver information
    type(TElectronicSolver), intent(inout) :: electronicSolver

    !> Inverse of eigenvectors matrix (for populations)
    complex(dp), intent(inout), allocatable :: Eiginv(:,:,:)

    !> Adjoint of the inverse of eigenvectors matrix (for populations)
    complex(dp), intent(inout), allocatable :: EiginvAdj(:,:,:)

    !> Square hamiltonian
    complex(dp), intent(in) :: H1(:,:,:)

    !> Square overlap matrix
    complex(dp), intent(inout) :: Ssqr(:,:,:)

    !> Auxiliary matrix
    complex(dp), allocatable :: T1(:,:)

    !> Auxiliary matrix
    real(dp), allocatable :: T2(:,:)

    !> K-Spin mixed index
    integer :: iKS

    real(dp) :: eigen(this%nOrbs)

    allocate(T1(this%nOrbs,this%nOrbs))
    allocate(T2(this%nOrbs,this%nOrbs))
    do iKS = 1, this%parallelKS%nLocalKS
      !check if this works with both complex and real
      T1(:,:) = H1(:,:,iKS)
      call diagDenseMtx(electronicSolver, 'V', T1, Ssqr(:,:,iKS), eigen)
      if (this%tRealHS) then
        T2(:,:) = real(T1, dp)
        call tdPopulInit(this, Eiginv(:,:,iKS), EiginvAdj(:,:,iKS), T2)
      else
        call tdPopulInit(this, Eiginv(:,:,iKS), EiginvAdj(:,:,iKS), eigvecsCplx=T1)
      end if
    end do
    deallocate(T1, T2)

  end subroutine updateBasisMatrices


  !> Calculate populations at each time step
  subroutine getTDPopulations(this, occ, rho, Eiginv, EiginvAdj, populDat, time, iKS)

    !> ElecDynamics instance
    type(TElecDynamics), intent(in) :: this

    !> Density Matrix
    complex(dp), intent(in) :: rho(:,:,:)

    !> Inverse of eigenvectors matrix (for populations)
    complex(dp), intent(inout), allocatable :: Eiginv(:,:,:)

    !> Adjoint of the inverse of eigenvectors matrix (for populations)
    complex(dp), intent(inout), allocatable :: EiginvAdj(:,:,:)

    !> Elapsed simulation time
    real(dp), intent(in) :: time

    !> Populations output file ID
    integer, intent(in) :: populDat(:)

    !> K-Spin mixed index
    integer, intent(in) :: iKS

    !> Molecular orbital occupations
    real(dp), intent(inout) :: occ(:)

    !> Auxiliary matrix
    complex(dp) :: T1(this%nOrbs,this%nOrbs)

    integer :: ii

    call gemm(T1, rho(:,:,iKS), EiginvAdj(:,:,iKS))
    T1 = transpose(Eiginv(:,:,iKS)) * T1

    occ = real(sum(T1,dim=1), dp)
    write(populDat(iKS),'(*(2x,F25.15))', advance='no') time * au__fs
    do ii = 1, size(occ)
      write(populDat(iKS),'(*(2x,F25.15))', advance='no')occ(ii)
    end do
    write(populDat(iKS),*)

  end subroutine getTDPopulations


  !> Write time-dependent tagged information to autotestTag file
  subroutine writeTDAutotest(this, dipole, energy, deltaQ, coord, totalForce, occ, lastBondPopul,&
      & taggedWriter)

    !> ElecDynamics instance
    type(TElecDynamics), intent(in) :: this

    !> Dipole moment
    real(dp), intent(in) :: dipole(:,:)

    !> data type for energy components and total
    type(TEnergies), intent(in) :: energy

    !> Negative gross charge
    real(dp), intent(in) :: deltaQ(:,:)

    !> atomic coordinates
    real(dp), intent(in) :: coord(:,:)

    !> forces (3, nAtom)
    real(dp), intent(in) :: totalForce(:,:)

    !> molecular orbital projected populations
    real(dp), intent(in) :: occ(:)

    !> Last bond population in the run
    real(dp), intent(in) :: lastBondPopul

    !> Tagged writer object
    type(TTaggedWriter), intent(inout) :: taggedWriter

    integer :: fdAutotest

    open(newunit=fdAutotest, file=trim(this%autotestTag), position="append")

    call taggedWriter%write(fdAutotest, tagLabels%tdenergy, energy%eSCC)
    call taggedWriter%write(fdAutotest, tagLabels%tddipole, dipole)
    call taggedWriter%write(fdAutotest, tagLabels%tdcharges, deltaQ)
    if (this%tIons) then
      call taggedWriter%write(fdAutotest, tagLabels%ehrencoords, coord)
      call taggedWriter%write(fdAutotest, tagLabels%ehrenvelos, this%movedVelo)
    end if
    if (this%tForces) then
      call taggedWriter%write(fdAutotest, tagLabels%ehrenforces, totalForce)
    end if
    if (this%tPopulations) then
      call taggedWriter%write(fdAutotest, tagLabels%tdprojocc, occ)
    end if
    if (this%tBondP) then
      call taggedWriter%write(fdAutotest, tagLabels%sumBondPopul, lastBondPopul)
    end if

    close(fdAutotest)

  end subroutine writeTDAutotest


  !> Initialize ion dynamics
  subroutine initIonDynamics(this, coordNew, coord, movedAccel)

    !> ElecDynamics instance
    type(TElecDynamics), intent(inout) :: this

    !> coordinates of next step
    real(dp), intent(out) :: coordNew(:,:)

    !> atomic coordinates
    real(dp), intent(in) :: coord(:,:)

    !> acceleration on moved atoms (3, nMovedAtom)
    real(dp), intent(in) :: movedAccel(:,:)

    ! Data for the velocity verlet integrator
    type(TVelocityVerlet), allocatable :: pVelocityVerlet

    if (this%nDynamicsInit == 0) then
      allocate(pVelocityVerlet)
    end if

    if (this%ReadMDVelocities) then
      this%movedVelo(:, :) = this%initialVelocities
    else
      this%movedVelo(:, :) = 0.0_dp
    end if

    if (this%nDynamicsInit == 0) then
      if (this%tReadRestart) then
        call init(pVelocityVerlet, this%dt, coord(:, this%indMovedAtom), this%pThermostat,&
            & this%movedVelo, this%ReadMDVelocities, tHalfVelocities=.true.)
      else
        call init(pVelocityVerlet, this%dt, coord(:, this%indMovedAtom), this%pThermostat,&
            & this%movedVelo, this%ReadMDVelocities, tHalfVelocities=.true.)
      end if
      this%initialVelocities(:, this%indMovedAtom) = this%movedVelo
    else
      call reset(this%pMDIntegrator, coordNew(:, this%indMovedAtom), this%initialVelocities,&
          & tHalfVelocities=.true.)
    end if

    ! Euler step from 1st VV step
    ! Ensures good initialization and puts velocity and coords on common time step
    this%movedVelo(:,:) = this%movedVelo - 0.5_dp * movedAccel * this%dt
    coordNew(:,:) = coord
    coordNew(:,this%indMovedAtom) = coord(:,this%indMovedAtom) &
        & + this%movedVelo(:,:) * this%dt + 0.5_dp * movedAccel(:,:) * this%dt**2
    ! This re-initializes the VVerlet propagator with coordNew
    this%movedVelo(:,:) = this%movedVelo + 0.5_dp * movedAccel * this%dt

    if (this%nDynamicsInit == 0) then
      call reset(pVelocityVerlet, coordNew(:, this%indMovedAtom), this%movedVelo,&
          & tHalfVelocities=.true.)
      allocate(this%pMDIntegrator)
      call init(this%pMDIntegrator, pVelocityVerlet)
    else
      call reset(this%pMDIntegrator, coordNew(:, this%indMovedAtom), this%movedVelo,&
          & tHalfVelocities=.true.)
    end if

    this%nDynamicsInit = this%nDynamicsInit + 1

  end subroutine initIonDynamics


  !> Calculates nonscc hamiltonian and overlap for new geometry and reallocates sparse arrays
  subroutine updateH0S(this, Ssqr, Sinv, coord, orb, neighbourList, nNeighbourSK, iSquare,&
      & iSparseStart, img2CentCell, skHamCont, skOverCont, ham, ham0, over, env, rhoPrim,&
      & ErhoPrim, coordAll)

    !> ElecDynamics instance
    type(TElecDynamics), intent(inout), target :: this

    !> Square overlap inverse
    complex(dp), intent(inout) :: Sinv(:,:,:)

    !> Square overlap matrix
    complex(dp), intent(inout), allocatable :: Ssqr(:,:,:)

    !> Local sparse storage for non-SCC hamiltonian
    real(dp), allocatable, intent(inout) :: ham0(:)

    !> scc hamiltonian (sparse)
    real(dp), allocatable, intent(inout) :: ham(:,:)

    !> overlap (sparse)
    real(dp), allocatable, intent(inout) :: over(:)

    !> atomic coordinates
    real(dp), allocatable, intent(inout) :: coord(:,:)

    !> Coords of the atoms (3, nAllAtom)
    real(dp), allocatable, intent(inout) :: coordAll(:,:)

    !> ADT for neighbour parameters
    type(TNeighbourList), intent(inout) :: neighbourList

    !> nr. of neighbours for atoms out to max interaction distance (excluding Ewald terms)
    integer, intent(inout) :: nNeighbourSK(:)

    !> index array for location of atomic blocks in large sparse arrays
    integer, allocatable, intent(inout) :: iSparseStart(:,:)

    !> image atoms to their equivalent in the central cell
    integer, allocatable, intent(inout) :: img2CentCell(:)

    !> Index array for start of atomic block in dense matrices
    integer, intent(in) :: iSquare(:)

    !> data type for atomic orbital information
    type(TOrbitals), intent(in) :: orb

    !> Raw H^0 hamiltonian data
    type(TSlakoCont), intent(in) :: skHamCont

    !> Raw overlap data
    type(TSlakoCont), intent(in) :: skOverCont

    !> Environment settings
    type(TEnvironment), intent(inout) :: env

    !> sparse density matrix
    real(dp), allocatable, intent(inout) :: rhoPrim(:,:)

    !> Energy weighted density matrix
    real(dp), allocatable, intent(inout) :: ErhoPrim(:)

    real(dp), allocatable :: Sreal(:,:), SinvReal(:,:)
    complex(dp), allocatable :: T4(:,:)
    real(dp) :: coord0Fold(3,this%nAtom)
    integer :: nAllAtom, iSpin, sparseSize, iOrb, iKS, iK

    coord0Fold(:,:) = coord
    if (this%tPeriodic) then
      call foldCoordToUnitCell(coord0Fold, this%latVec, this%invLatVec)
    end if

    call updateNeighbourListAndSpecies(coordAll, this%speciesAll, img2CentCell, this%iCellVec, &
        &neighbourList, nAllAtom, coord0Fold, this%species, this%mCutoff, this%rCellVec)
    call getNrOfNeighboursForAll(nNeighbourSK, neighbourList, this%skCutoff)
    call getSparseDescriptor(neighbourList%iNeighbour, nNeighbourSK, img2CentCell, orb,&
        & iSparseStart, sparseSize)

    this%nSparse = sparseSize
    if (.not. allocated(ham0)) then
      allocate(ham0(this%nSparse))
    end if
    if (.not. allocated(rhoPrim)) then
      allocate(rhoPrim(this%nSparse, this%nSpin))
    end if
    call reallocateTDSparseArrays(this, ham, over, ham0, rhoPrim, ErhoPrim)

    call this%sccCalc%updateCoords(env, coordAll, this%speciesAll, neighbourList)

    if (allocated(this%dispersion)) then
      call this%dispersion%updateCoords(env, neighbourList, img2CentCell, coordAll,&
          & this%speciesAll)
    end if

    call buildH0(env, ham0, skHamCont, this%atomEigVal, coordAll, nNeighbourSK, &
        & neighbourList%iNeighbour, this%speciesAll, iSparseStart, orb)
    call buildS(env, over, skOverCont, coordAll, nNeighbourSK, neighbourList%iNeighbour,&
        & this%speciesAll, iSparseStart, orb)

    if (this%tRealHS) then
      allocate(Sreal(this%nOrbs,this%nOrbs))
      allocate(SinvReal(this%nOrbs,this%nOrbs))
      Sreal = 0.0_dp
      call unpackHS(Sreal, over, neighbourList%iNeighbour, nNeighbourSK, iSquare, iSparseStart,&
          & img2CentCell)
      call blockSymmetrizeHS(Sreal, iSquare)
      do iKS = 1, this%parallelKS%nLocalKS
        Ssqr(:,:,iKS) = cmplx(Sreal, 0, dp)
      end do

      SinvReal = 0.0_dp
      do iOrb = 1, this%nOrbs
        SinvReal(iOrb, iOrb) = 1.0_dp
      end do
      call gesv(Sreal, SinvReal)

      do iKS = 1, this%parallelKS%nLocalKS
        Sinv(:,:,iKS) = cmplx(SinvReal, 0, dp)
      end do

    else

      allocate(T4(this%nOrbs,this%nOrbs))
      Ssqr(:,:,:) = cmplx(0,0,dp)
      do iKS = 1, this%parallelKS%nLocalKS
        iK = this%parallelKS%localKS(1, iKS)
        iSpin = this%parallelKS%localKS(2, iKS)
        T4(:,:) = cmplx(0,0,dp)
        call unpackHS(T4, over, this%kPoint(:,iK), neighbourList%iNeighbour, nNeighbourSK,&
            & this%iCellVec, this%cellVec, iSquare, iSparseStart, img2CentCell)
        call blockHermitianHS(T4, iSquare)
        Ssqr(:,:,iKS) = T4
        Sinv(:,:,iKS) = cmplx(0,0,dp)
        do iOrb = 1, this%nOrbs
          Sinv(iOrb, iOrb, iKS) = cmplx(1,0,dp)
        end do
        call gesv(T4, Sinv(:,:,iKS))
      end do
      deallocate(T4)

    end if

  end subroutine updateH0S


  !> Calculates force
  subroutine getForces(this, movedAccel, totalForce, rho, H1, Sinv, neighbourList, nNeighbourSK,&
      & img2CentCell, iSparseStart, iSquare, potential, orb, skHamCont, skOverCont, qq, q0,&
      & pRepCont, coordAll, rhoPrim, ErhoPrim, iStep, env, rangeSep, deltaRho)

    !> ElecDynamics instance
    type(TElecDynamics), intent(inout) :: this

    !> Density Matrix
    complex(dp), intent(in) :: rho(:,:,:)

    !> Square hamiltonian
    complex(dp), intent(in) :: H1(:,:,:)

    !> Square inverse overlap
    complex(dp), intent(in) :: Sinv(:,:,:)

    !> ADT for neighbour parameters
    type(TNeighbourList), intent(inout) :: neighbourList

    !> nr. of neighbours for atoms out to max interaction distance (excluding Ewald terms)
    integer, intent(in) :: nNeighbourSK(:)

    !> index array for location of atomic blocks in large sparse arrays
    integer, intent(in) :: iSparseStart(:,:)

    !> image atoms to their equivalent in the central cell
    integer, intent(in) :: img2CentCell(:)

    !> Index array for start of atomic block in dense matrices
    integer, intent(in) :: iSquare(:)

    !> data type for atomic orbital information
    type(TOrbitals), intent(in) :: orb

    !> Raw H^0 hamiltonian data
    type(TSlakoCont), intent(in) :: skHamCont

    !> Raw overlap data
    type(TSlakoCont), intent(in) :: skOverCont

    !> sparse density matrix
    real(dp), intent(inout) :: rhoPrim(:,:)

    !> Energy weighted density matrix
    real(dp), intent(inout) :: ErhoPrim(:)

    !> forces (3, nAtom)
    real(dp), intent(out) :: totalForce(:,:)

    !> acceleration on moved atoms (3, nMovedAtom)
    real(dp), intent(out) :: movedAccel(:,:)

    !> potential acting on the system
    type(TPotentials), intent(in) :: potential

    !> atomic ocupations
    real(dp), intent(inout) :: qq(:,:,:)

    !> reference atomic occupations
    real(dp), intent(inout) :: q0(:,:,:)

    !> repulsive information
    type(TRepCont), intent(in) :: pRepCont

    !> Coords of the atoms (3, nAllAtom)
    real(dp), intent(in) :: coordAll(:,:)

    !> current step of the propagation
    integer, intent(in) :: iStep

    !> Environment settings
    type(TEnvironment), intent(inout) :: env

    !> Range separation contributions
    type(TRangeSepFunc), allocatable, intent(inout) :: rangeSep

    !> Real part of density matrix, adjusted by reference charges
    complex(dp), allocatable, intent(inout) :: deltaRho(:,:,:)

    real(dp), allocatable :: T1R(:,:), T2R(:,:)
    complex(dp), allocatable :: T1C(:,:), T2C(:,:)
    real(dp) :: derivs(3,this%nAtom), repulsiveDerivs(3,this%nAtom), totalDeriv(3, this%nAtom)
    integer :: iSpin, iDir, iKS, iK

    ErhoPrim(:) = 0.0_dp
    rhoPrim(:,:) = 0.0_dp

    !    do iSpin = 1, this%nSpin
    !       call gemm(T1, real(rho(:,:,iSpin), dp), real(H1(:,:,iSpin), dp))
    !       call her2k(T2, real(Sinv(:,:,iSpin), dp), T1, 0.5_dp)
    !       call packHS(rhoPrim(:,iSpin), real(rho(:,:,iSpin), dp), neighbourList%iNeighbour,&
    !           & nNeighbourSK, orb%mOrb, iSquare, iSparseStart, img2CentCell)
    !       call packHS(ErhoPrim, T2, neighbourList%iNeighbour, nNeighbourSK, orb%mOrb, iSquare,&
    !           & iSparseStart, img2CentCell)
    !    end do

    if (this%tRealHS) then
      allocate(T1R(this%nOrbs,this%nOrbs))
      allocate(T2R(this%nOrbs,this%nOrbs))
      do iKS = 1, this%parallelKS%nLocalKS
        iK = this%parallelKS%localKS(1, iKS)
        iSpin = this%parallelKS%localKS(2, iKS)
        call packHS(rhoPrim(:,iSpin), real(rho(:,:,iKS), dp), neighbourList%iNeighbour,&
            & nNeighbourSK, orb%mOrb, iSquare, iSparseStart, img2CentCell)
        call gemm(T1R, real(rho(:,:,iKS), dp), real(H1(:,:,iKS), dp))
        call her2k(T2R, real(Sinv(:,:,iKS), dp), T1R, 0.5_dp)
        call packHS(ErhoPrim, T2R, neighbourList%iNeighbour, nNeighbourSK, orb%mOrb, iSquare,&
            & iSparseStart, img2CentCell)
      end do
    else
      allocate(T1C(this%nOrbs,this%nOrbs))
      allocate(T2C(this%nOrbs,this%nOrbs))
      do iKS = 1, this%parallelKS%nLocalKS
        iK = this%parallelKS%localKS(1, iKS)
        iSpin = this%parallelKS%localKS(2, iKS)
        call packHS(rhoPrim(:,iSpin), rho(:,:,iKS), this%kPoint(:,iK), this%kWeight(iK),&
            & neighbourList%iNeighbour, nNeighbourSK, orb%mOrb, this%iCellVec, this%cellVec,&
            & iSquare, iSparseStart, img2CentCell)
        call gemm(T1C, rho(:,:,iKS), H1(:,:,iKS))
        call her2k(T2C, Sinv(:,:,iKS), T1C, (0.5_dp,0.0_dp))
        call packHS(ErhoPrim, T2C, this%kPoint(:,iK), this%kWeight(iK), neighbourList%iNeighbour,&
            & nNeighbourSK, orb%mOrb, this%iCellVec, this%cellVec, iSquare, iSparseStart,&
            & img2CentCell)
      end do
    end if

    call ud2qm(qq)
    call ud2qm(q0)
    call ud2qm(rhoPrim)

    derivs(:,:) = 0.0_dp
    repulsiveDerivs(:,:) = 0.0_dp

    call derivative_shift(env, derivs, this%derivator, rhoPrim, ErhoPrim, skHamCont,&
        & skOverCont, coordAll, this%speciesAll, neighbourList%iNeighbour, nNeighbourSK, &
        & img2CentCell, iSparseStart, orb, potential%intBlock)
    call this%sccCalc%updateCharges(env, qq, q0, orb, this%speciesAll)
    call this%sccCalc%addForceDc(env, derivs, this%speciesAll, neighbourList%iNeighbour, &
        & img2CentCell)
    call getERepDeriv(repulsiveDerivs, coordAll, nNeighbourSK, neighbourList%iNeighbour,&
        & this%speciesAll, pRepCont, img2CentCell)

    if (this%isRangeSep) then
      call error("Ehrenfest forces not implemented yet with range separated calculations.")
      !call rangeSep%addLRGradients(derivs, this%derivator, deltaRho, skHamCont, skOverCont,&
      ! & coordAll, this%speciesAll, orb, iSquare, sSqr, neighbourList%iNeighbour, nNeighbourSK)
    end if

    if (this%tLaser) then
      do iDir = 1, 3
        derivs(iDir,:) = derivs(iDir,:)&
            & - sum(q0(:,:,1) - qq(:,:,1), dim=1) * this%TDFunction(iDir, iStep)
      end do
    end if

    totalDeriv(:,:) = repulsiveDerivs + derivs
    if (allocated(this%dispersion)) then
      call this%dispersion%addGradients(totalDeriv)
    end if

    totalForce(:,:) = - totalDeriv
    movedAccel(:,:) = totalForce(:, this%indMovedAtom) / this%movedMass

    call qm2ud(qq)
    call qm2ud(q0)
    call qm2ud(rhoPrim)

  end subroutine getForces


  !> Calculates nonadiabatic matrix: overlap gradient (Sprime) times velocities (Rdot)
  subroutine getRdotSprime(this, RdotSprime, coordAll, skOverCont, orb, img2CentCell, &
      &neighbourList, nNeighbourSK, iSquare)

    !> ElecDynamics instance
    type(TElecDynamics), intent(in), target :: this

    !> Raw overlap data
    type(TSlakoCont), intent(in) :: skOverCont

    ! nonadiabatic coupling matrix elements
    complex(dp), intent(out) :: RdotSprime(:,:)

    !> data type for atomic orbital information
    type(TOrbitals), intent(in) :: orb

    !> Coords of the atoms (3, nAllAtom)
    real(dp), intent(in) :: coordAll(:,:)

    !> ADT for neighbour parameters
    type(TNeighbourList), intent(in) :: neighbourList

    !> nr. of neighbours for atoms out to max interaction distance (excluding Ewald terms)
    integer, intent(in) :: nNeighbourSK(:)

    !> Index array for start of atomic block in dense matrices
    integer, intent(in) :: iSquare(:)

    !> image atoms to their equivalent in the central cell
    integer, intent(in) :: img2CentCell(:)

    real(dp) :: sPrimeTmp(orb%mOrb,orb%mOrb,3)
    real(dp) :: sPrimeTmp2(orb%mOrb,orb%mOrb), dcoord(3,this%nAtom)
    integer :: iAtom1, iStart1, iEnd1, iSp1, nOrb1, iDir
    integer :: iNeigh, iStart2, iEnd2, iAtom2, iAtom2f, iSp2, nOrb2

    dcoord(:,:) = 0.0_dp
    dcoord(:, this%indMovedAtom) = this%movedVelo(:,1:this%nMovedAtom)
    sPrimeTmp(:,:,:) = 0.0_dp
    RdotSprime(:,:) = 0.0_dp

    !$OMP PARALLEL DO PRIVATE(iAtom1,iStart1,iEnd1,iSp1,nOrb1,sPrimeTmp2,iNeigh,iAtom2, &
    !$OMP& iAtom2f,iStart2,iEnd2,iSp2,nOrb2,sPrimeTmp,iDir) DEFAULT(SHARED) &
    !$OMP& SCHEDULE(RUNTIME)
    do iAtom1 = 1, this%nAtom
      iStart1 = iSquare(iAtom1)
      iEnd1 = iSquare(iAtom1+1)-1
      iSp1 = this%species(iAtom1)
      nOrb1 = orb%nOrbAtom(iAtom1)

      do iNeigh = 1, nNeighbourSK(iAtom1)
        iAtom2 = neighbourList%iNeighbour(iNeigh, iAtom1)
        iAtom2f = img2CentCell(iAtom2)
        iStart2 = iSquare(iAtom2f)
        iEnd2 = iSquare(iAtom2f+1)-1
        iSp2 = this%species(iAtom2f)
        nOrb2 = orb%nOrbAtom(iAtom2f)
        if (iAtom2f /= iAtom1) then
          call this%derivator%getFirstDeriv(sPrimeTmp, skOverCont, coordAll, this%speciesAll,&
              & iAtom1, iAtom2, orb)

          sPrimeTmp2(:,:) = 0.0_dp
          do iDir=1,3
            sPrimeTmp2(:,:) = sPrimeTmp2 + sPrimeTmp(:,:,iDir) * dcoord(iDir,iAtom1)
          end do
          RdotSprime(iStart2:iEnd2,iStart1:iEnd1) = RdotSprime(iStart2:iEnd2,iStart1:iEnd1)&
              & + cmplx(sPrimeTmp2(1:nOrb2,1:nOrb1), 0, dp)

          sPrimeTmp2(:,:) = 0.0_dp
          do iDir=1,3
            sPrimeTmp2(:,:) = sPrimeTmp2 - sPrimeTmp(:,:,iDir) * dcoord(iDir,iAtom2f)
          end do
          RdotSprime(iStart1:iEnd1,iStart2:iEnd2) = RdotSprime(iStart1:iEnd1,iStart2:iEnd2)&
              & + cmplx(transpose(sPrimeTmp2(1:nOrb2,1:nOrb1)), 0, dp)
        end if
      end do
    end do
    !$OMP END PARALLEL DO

  end subroutine getRdotSprime


  !> Reallocates sparse arrays after change of coordinates
  subroutine reallocateTDSparseArrays(this, ham, over, ham0, rhoPrim, ErhoPrim)

    !> ElecDynamics instance
    type(TElecDynamics), intent(in), target :: this

    !> scc hamiltonian (sparse)
    real(dp), allocatable, intent(inout) :: ham(:,:)

    !> overlap (sparse)
    real(dp), allocatable, intent(inout) :: over(:)

    !> Local sparse storage for non-SCC hamiltonian
    real(dp), allocatable, intent(inout) :: ham0(:)

    !> sparse density matrix
    real(dp), allocatable, intent(inout) :: rhoPrim(:,:)

    !> Energy weighted density matrix
    real(dp), allocatable, intent(inout) :: ErhoPrim(:)

    deallocate(ham)
    deallocate(over)
    deallocate(ham0)
    deallocate(rhoPrim)
    allocate(ham(this%nSparse, this%nSpin))
    allocate(over(this%nSparse))
    allocate(ham0(this%nSparse))
    allocate(rhoPrim(this%nSparse, this%nSpin))

    if (allocated(ErhoPrim)) then
      deallocate(ErhoPrim)
      allocate(ErhoPrim(this%nSparse))
    end if

  end subroutine reallocateTDSparseArrays


  !updates SCC module with lattice vectors
  subroutine initLatticeVectors(this)

    !> ElecDynamics instance
    type(TElecDynamics), intent(inout), target :: this

    real(dp) :: cellVol, recVecs(3,3), recVecs2p(3,3)

    cellVol = abs(determinant33(this%latVec))
    recVecs2p(:,:) = this%latVec
    call matinv(recVecs2p)
    recVecs2p = transpose(recVecs2p)
    recVecs = 2.0_dp * pi * recVecs2p
    call this%sccCalc%updateLatVecs(this%latVec, recVecs, cellVol)
    this%mCutOff = max(this%mCutOff, this%sccCalc%getCutOff())

    if (allocated(this%dispersion)) then
      call this%dispersion%updateLatVecs(this%latVec)
      this%mCutOff = max(this%mCutOff, this%dispersion%getRCutOff())
    end if
  end subroutine initLatticeVectors


  !> Calculates repulsive and dispersion energies
  subroutine  getPositionDependentEnergy(this, energy, coordAll, img2CentCell, nNeighbourSK,&
      & neighbourList, pRepCont, iAtInCentralRegion)

    !> ElecDynamics instance
    type(TElecDynamics), intent(inout), target :: this

    !> data type for energy components and total
    type(TEnergies), intent(inout) :: energy

    !> All atomic coordinates
    real(dp), intent(in) :: coordAll(:,:)

    !> Image atom indices to central cell atoms
    integer, intent(in) :: img2CentCell(:)

    !> Number of neighbours for each of the atoms
    integer, intent(in) :: nNeighbourSK(:)

    !> List of neighbours for each atom
    type(TNeighbourList), intent(in) :: neighbourList

    !> Repulsive interaction data
    type(TRepCont), intent(in) :: pRepCont

    !> atoms in the central cell
    integer, intent(in) :: iAtInCentralRegion(:)

    ! using here nNeighbourSK instead of nNeighbourRep
    call calcRepulsiveEnergy(coordAll, this%speciesAll, img2CentCell, nNeighbourSK, neighbourList,&
        & pRepCont, energy%atomRep, energy%Erep, iAtInCentralRegion)

    if (allocated(this%dispersion)) then
      call calcDispersionEnergy(this%dispersion, energy%atomDisp, energy%eDisp, iAtInCentralRegion)
    else
      energy%atomDisp(:) = 0.0_dp
      energy%eDisp = 0.0_dp
    end if

  end subroutine getPositionDependentEnergy


  !> Calculates bond populations and bond energies if requested
  subroutine getBondPopulAndEnergy(this, bondWork, lastBondPopul, rhoPrim, ham0, over, iNeighbour,&
      & nNeighbourSK, iSparseStart, img2CentCell, iSquare,  fdBondEnergy, fdBondPopul,&
      & tWriteAutotest, time)

    !> ElecDynamics instance
    type(TElecDynamics), intent(inout) :: this

    !> Container for either bond populations or bond energy
    real(dp), allocatable, intent(inout) :: bondWork(:, :)

    !> Last calculated bond population (for tagged output)
    real(dp), intent(inout) :: lastBondPopul

    !> sparse density matrix
    real(dp), intent(in) :: rhoPrim(:,:)

    !> overlap (sparse)
    real(dp), intent(in) :: over(:)

    !> Local sparse storage for non-SCC hamiltonian
    real(dp), intent(in) :: ham0(:)

    !> Atomic neighbour data
    integer, intent(in) :: iNeighbour(0:,:)

    !> nr. of neighbours for atoms out to max interaction distance (excluding Ewald terms)
    integer, intent(in) :: nNeighbourSK(:)

    !> index array for location of atomic blocks in large sparse arrays
    integer, intent(in) :: iSparseStart(:,:)

    !> image atoms to their equivalent in the central cell
    integer, intent(in) :: img2CentCell(:)

    !> Index array for start of atomic block in dense matrices
    integer, intent(in) :: iSquare(:)

    !> File descriptor for bond energy
    integer, intent(in) :: fdBondEnergy

    !> File descriptor for bond populations
    integer, intent(in) :: fdBondPopul

    !> produce tagged output?
    logical, intent(in) :: tWriteAutotest

    !> Elapsed simulation time
    real(dp), intent(in) :: time

    integer :: iSpin

    if (this%tBondE) then
      bondWork(:,:) = 0.0_dp
      do iSpin = 1, this%nSpin
        call addPairWiseBondInfo(bondWork, rhoPrim(:,iSpin), ham0, iSquare,&
            & iNeighbour, nNeighbourSK, img2CentCell, iSparseStart)
      end do
      write(fdBondEnergy) time * au__fs, sum(bondWork), bondWork
    end if
    if (this%tBondP) then
      bondWork(:,:) = 0.0_dp
      do iSpin = 1, this%nSpin
        call addPairWiseBondInfo(bondWork, rhoPrim(:,1), over, iSquare,&
            & iNeighbour, nNeighbourSK, img2CentCell, iSparseStart)
      end do
      write(fdBondPopul) time * au__fs, sum(bondWork), bondWork
      if (tWriteAutotest) then
        lastBondPopul = sum(bondWork)
      end if
    end if

  end subroutine getBondPopulAndEnergy

end module dftbp_timeprop
