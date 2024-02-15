!--------------------------------------------------------------------------------------------------!
!  DFTB+: general package for performing fast atomistic simulations                                !
!  Copyright (C) 2006 - 2023  DFTB+ developers group                                               !
!                                                                                                  !
!  See the LICENSE file for terms of usage and distribution.                                       !
!--------------------------------------------------------------------------------------------------!

#:include 'common.fypp'
#:include 'error.fypp'

!> Implements real-time Ehrenfest time-dependent DFTB by numerically propagating the electronic
!! one-electron density matrix of the system and the nuclei in the presence of an external
!! perturbation (kick or laser) Bonafé, F. P., Aradi, B., Hourahine, B., Medrano, C. R., Hernandez,
!! F. J., Frauenheim, T., & Sánchez, C. G.  Journal of Chemical Theory and Computation (2020)
!! https://doi.org/10.1021/acs.jctc.9b01217
module dftbp_timedep_timeprop
  use dftbp_common_accuracy, only : dp, sc, lc, mc
  use dftbp_common_constants, only : au__fs, pi, Bohr__AA, imag, Hartree__eV
  use dftbp_common_environment, only : TEnvironment, globalTimers
  use dftbp_common_file, only : TFileDescr, TOpenOptions, openFile, closeFile
  use dftbp_common_globalenv, only : stdOut
  use dftbp_common_hamiltoniantypes, only : hamiltonianTypes
  use dftbp_common_status, only : TStatus
  use dftbp_common_timer, only : TTimer
  use dftbp_dftb_bondpopulations, only : addPairWiseBondInfo
  use dftbp_dftb_boundarycond, only : TBoundaryConditions
  use dftbp_dftb_densitymatrix, only : makeDensityMatrix
  use dftbp_dftb_dftbplusu, only : TDftbU
  use dftbp_dftb_dispersions, only : TDispersionIface
  use dftbp_dftb_energytypes, only : TEnergies, TEnergies_init
  use dftbp_dftb_forces, only : derivative_shift
  use dftbp_dftb_getenergies, only : calcEnergies, calcDispersionEnergy, sumEnergies
  use dftbp_dftb_hamiltonian, only : TRefExtPot, resetExternalPotentials, resetInternalPotentials,&
      & addBlockChargePotentials, addChargePotentials, getSccHamiltonian
  use dftbp_dftb_nonscc, only : TNonSccDiff, buildH0, buildS
  use dftbp_dftb_onsitecorrection, only : addOnsShift
  use dftbp_dftb_periodic, only : TNeighbourList, updateNeighbourListAndSpecies,&
      & getNrOfNeighboursForAll
  use dftbp_dftb_populations, only :  getChargePerShell, denseSubtractDensityOfAtoms,&
      & denseSubtractDensityOfAtoms_nospin_real_nonperiodic_reks,&
      & denseSubtractDensityOfAtoms_spin_real_nonperiodic_reks
  use dftbp_dftb_potentials, only : TPotentials, TPotentials_init
  use dftbp_dftb_hybridxc, only : THybridXcFunc
  use dftbp_dftb_repulsive_repulsive, only : TRepulsive
  use dftbp_dftb_scc, only : TScc
  use dftbp_dftb_shift, only : totalShift
  use dftbp_dftb_slakocont, only : TSlakoCont
  use dftbp_dftb_sparse2dense, only : packHS, unpackHS, unpackDQ, getSparseDescriptor
  use dftbp_dftb_spin, only : ud2qm, qm2ud
  use dftbp_dftb_thirdorder, only : TThirdOrder
  use dftbp_dftbplus_eigenvects, only : diagDenseMtx
  use dftbp_dftbplus_qdepextpotproxy, only : TQDepExtPotProxy
  use dftbp_elecsolvers_elecsolvers, only : TElectronicSolver
  use dftbp_extlibs_tblite, only : TTBLite
  use dftbp_io_message, only : warning
  use dftbp_io_taggedoutput, only : TTaggedWriter, tagLabels
  use dftbp_math_blasroutines, only : gemm, her2k
  use dftbp_math_lapackroutines, only : matinv, gesv
  use dftbp_math_matrixops, only : adjointLowerTriangle
  use dftbp_math_ranlux, only : TRanlux
  use dftbp_math_simplealgebra, only : determinant33
  use dftbp_md_dummytherm, only : TDummyThermostat
  use dftbp_md_mdcommon, only : TMDCommon
  use dftbp_md_mdintegrator, only : TMDIntegrator, reset, init, state, next
  use dftbp_md_thermostat, only : TThermostat
  use dftbp_md_velocityverlet, only : TVelocityVerlet
  use dftbp_reks_reks, only : TReksCalc
  use dftbp_solvation_solvation, only : TSolvation
  use dftbp_solvation_fieldscaling, only : TScaleExtEField
  use dftbp_timedep_dynamicsrestart, only : writeRestartFile, readRestartFile
  use dftbp_type_commontypes, only : TParallelKS, TOrbitals
  use dftbp_type_densedescr, only: TDenseDescr
  use dftbp_type_integral, only : TIntegral
  use dftbp_type_multipole, only : TMultipole, TMultipole_init
  use dftbp_io_message, only : error, warning
#:if WITH_SCALAPACK
  use dftbp_dftb_densitymatrix, only : makeDensityMtxRealBlacs
  use dftbp_dftb_sparse2dense, only : unpackHSRealBlacs, packRhoRealBlacs
  use dftbp_dftb_populations, only : mulliken
  use dftbp_extlibs_scalapackfx, only : pblasfx_psymm, pblasfx_ptran
  use dftbp_extlibs_mpifx, only : MPI_SUM, mpifx_allreduceip
  use dftbp_math_scalafxext, only : psymmatinv, phermatinv
#:endif
#:if WITH_MBD
  use dftbp_dftb_dispmbd, only : TDispMbd
#:endif
  implicit none

  private
  public :: runDynamics, TElecDynamics_init
  public :: initializeDynamics, finalizeDynamics, doTdStep
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

    !> Number of steps at which an Euler step is applied (to kill numerical noise)
    integer :: eulerFreq

    !> Number of total system snapshots to be saved during pump simulation
    integer :: tdPPFrames

    !> Number of excited atoms
    integer :: nExcitedAtom

    !> Initial and final times to save the snapshots during pump simulation
    real(dp) :: tdPpRange(2)

    !> If forces should be calculated during propagation
    logical :: tForces

    !> If nuclei should be moved during propagation
    logical :: tIons

    !> If velocities are supplied from input
    logical :: tReadMDVelocities

    !> If Euler steps should be done during simulation
    logical :: tEulers

    !> If pairwise bond energy should be calculated and written
    logical :: tBondE

    !> If pairwise bond population should be calculated and written
    logical :: tBondP

    !> If atom-resolved energies should be written
    logical :: tWriteAtomEnergies

    !> If this is a pump trajectory (for a pump-probe simulation)
    logical :: tPump

    !> If this is a probe trajectory (for a pump-probe simulation)
    logical :: tProbe

    !> Atomic (initial) kinetic temperature
    real(dp) :: tempAtom

    !> Field strength for KickAndLaser perturbations
    real(dp) :: tdLaserField = 0.0_dp

    !> Initial atomic velocities if supplied
    real(dp), allocatable :: initialVelocities(:,:)

    !> If initial fillings are provided in an external file
    logical :: tFillingsFromFile

  end type TElecDynamicsInp


  !> Data type for electronic dynamics internal settings
  type TElecDynamics
    private

    !> External field peak intensity
    real(dp) :: field

    !> Electric field angular frequency
    real(dp) :: omega

    !> Initial time of pulse
    real(dp) :: time0

    !> Final time of pulse
    real(dp) :: time1

    !> Phase applied of the laser field
    real(dp) :: phase

    !> Cartesian components of time dependent electric field at each time step
    real(dp), allocatable :: tdFunction(:, :)

    !> Complex electric field direction
    complex(dp) :: fieldDir(3)

    !> Frequency of output writing for charges and populations
    integer :: writeFreq

    !> Type of external perturbation applied
    integer :: pertType

    !> Envelope shape for laser perturbation
    integer :: envType

    !> Spin type of absorption spectrum for spin polarised calculations (singlet/triplet)
    integer :: spType

    !> Number of atoms in the central cell
    integer :: nAtom

    !> Number of basis functions
    integer :: nOrbs

    !> Number of spin channels
    integer :: nSpin = 1

    !> Electric field polarisation direction (x,y,z)
    integer :: currPolDir = 1

    !> How often to save state to disc
    integer :: restartFreq

    !> Number of atomic dipoles
    integer :: nDipole = 0

    !> Number of atomic quadrupoles
    integer :: nQuadrupole = 0

    !> Should additional data files be printed, covering atom dynamics, charges and energies
    logical :: tdWriteExtras

    !> Atomic species
    integer, allocatable :: species(:)

    !> Electric field polarisation directions to evaluate
    integer, allocatable :: polDirs(:)

    !> Species of all atoms in the system, including image atoms
    integer, allocatable :: speciesAll(:)

    !> Labels for atomic species
    character(mc), allocatable :: speciesName(:)

    !> If time-dependent populations should be calculated
    logical :: tPopulations

    !> If calculation should be restarted from dump file
    logical :: tReadRestart

    !> If dump file should be written during the dynamics
    logical :: tWriteRestart

    !> If a dump file is read, should it be ascii (T) or binary (F)
    logical :: tRestartAscii

    !> If a dump file is read, should it be ascii (T) or binary (F)
    logical :: tWriteRestartAscii

    !> Produce tagged output
    logical :: tWriteAutotest

    !> Is a laser field applied
    logical :: tLaser = .false.

    !> Is an initial kick applied
    logical :: tKick = .false.

    !> Is both a laser and kick applied
    logical :: tKickAndLaser = .false.

    !> Read laser field data from a file
    logical :: tEnvFromFile = .false.

    !> Hamiltonian model
    integer :: hamiltonianType

    !> Is this a calculation involving electrostatics
    type(TScc), allocatable :: sccCalc

    !> TBLITE object
    type(TTBLite), allocatable :: tblite

    !> Electrostatic multipoles
    type(TMultipole) :: multipole

    !> Tagged output files (machine readable)
    character(mc) :: autotestTag

    !> Initial atomic velocities, if supplied
    real(dp), allocatable :: initialVelocities(:,:)

    !> Masses of moving atoms
    real(dp), allocatable :: movedMass(:,:)

    !> Longest cutoff distance for neighbours
    real(dp) :: mCutoff

    !> Slater-Koster cutoff distance
    real(dp) :: skCutoff

    !> Electric field magnitude for laster
    real(dp) :: laserField

    !> Real space cell vectors for periodic structures
    real(dp), allocatable :: rCellVec(:,:)

    !> Cell vectors (in units of the lattice constants) to cells of the lattice
    real(dp), allocatable :: cellVec(:,:)

    !> The k-points of the system
    real(dp), allocatable :: kPoint(:,:)

    !> Weights for the k-points
    real(dp), allocatable :: KWeight(:)

    !> Atomic (on-site) orbital energies
    real(dp), allocatable :: atomEigVal(:,:)

    !> Number of atoms receiving initial excitation kicks
    integer :: nExcitedAtom

    !> Number of atoms moving
    integer :: nMovedAtom

    !> Number of elements in the sparse  hamiltonian
    integer :: nSparse

    !> Number of steps at which an Euler step is applied (to kill numerical noise)
    integer :: eulerFreq

    !> Frequency of writing pump-probe data
    integer :: PpFreq

    !> Initial time to start pump-probe writing
    integer :: PpIni

    !> Final time for pump-probe writing
    integer :: PpEnd

    !> Index for the cell vector in which a periodic image of an atom is sited
    integer, allocatable :: iCellVec(:)

    !> List of atoms being excited in a kick
    integer, allocatable :: indExcitedAtom(:)

    !> Evaluate and store forces
    logical :: tForces

    !> Read initial molecular dynamics velocities
    logical :: ReadMDVelocities

    !> Is a pump pulse/beam used
    logical :: tPump

    !> Is a probe pulse/beam used
    logical :: tProbe

    !> Is this a real hamiltonian
    logical :: tRealHS

    !> Whether this is a hybrid xc-functional calculation
    logical :: isHybridXc

    !> If Euler steps should be done during simulation
    logical :: tEulers = .false.

    !> Should bond energy be evaluated?
    logical :: tBondE = .false.

    !> Should bond populations be evaluated?
    logical :: tBondP = .false.

    !> Is this a periodic geometry
    logical :: tPeriodic = .false.

    !> Read electronic fillings from a file
    logical :: tFillingsFromFile = .false.

    !> Are net populations (onsite only) evaluated
    logical :: tNetCharges = .false.

    !> Write atom resolved energies
    logical :: tWriteAtomEnergies = .false.

    !> Thermostat
    type(TThermostat), allocatable :: pThermostat

    !> Molecular dynamics integrator
    type(TMDIntegrator), allocatable :: pMDIntegrator

    !> Dispersion interaction container
    class(TDispersionIface), allocatable :: dispersion

    !> Differentiator for properties like overlap matrix
    type(TNonSccDiff), allocatable :: derivator

    !> Container for parallel distribution pattern of k-points and spins
    type(TParallelKS), allocatable :: parallelKS

    !> Lattice vectors
    real(dp), allocatable :: latVec(:,:)

    !> Inverse of the lattice vectors
    real(dp), allocatable :: invLatVec(:,:)

    !> Initial coordinates
    real(dp), allocatable :: initCoord(:,:)

    !> Square matrix for the overlap
    complex(dp), allocatable :: Ssqr(:,:,:)

    !> Square matrix for the inverse of the overlap
    complex(dp), allocatable :: Sinv(:,:,:)

    !> Dipole matrix
    complex(dp), allocatable :: Dsqr(:,:,:,:)

    !> Quadrupole matrix
    complex(dp), allocatable :: Qsqr(:,:,:,:)

    !> Square hamiltonian at each spin and k-point
    complex(dp), allocatable :: H1(:,:,:)

    !> Non-adiabatic coupling matrix elements
    complex(dp), allocatable :: RdotSprime(:,:)

    !> Density matrix pointer
    complex(dp), pointer :: rho(:,:,:)

    !> Pointer to density matrix before update
    complex(dp), pointer :: rhoOld(:,:,:)

    !> Density matrix target
    complex(dp), allocatable :: trho(:,:,:)

    !> Target to density matrix before update
    complex(dp), allocatable :: trhoOld(:,:,:)

    !> Atomic populations
    real(dp), allocatable :: qq(:,:,:)

    !> Density matrix in sparse storage
    real(dp), allocatable :: rhoPrim(:,:)

    !> Non-self-consistent hamiltonian components
    real(dp), allocatable :: ham0(:)

    !> Energy-weighted density matrix
    real(dp), allocatable :: ErhoPrim(:)

    !> Shell-resolved Mulliken charges
    real(dp), allocatable :: chargePerShell(:,:,:)

    !> Range separated hamiltonian
    complex(dp), allocatable :: H1LC(:,:)

    !> Density matrix, adjusted by reference charges
    complex(dp), allocatable :: deltaRho(:,:,:)

    !> Acceleration of atoms that can move
    real(dp), allocatable :: movedAccel(:,:)

    !> Block populations
    real(dp), allocatable :: qBlock(:,:,:,:)

    !> Net (on-site only) atomic charge
    real(dp), allocatable :: qNetAtom(:)

    !> Inverse of eigenvectors matrix (for populations)
    complex(dp), allocatable :: Eiginv(:,:,:)

    !> Adjoint of inverse of eigenvectors matrix
    complex(dp), allocatable :: EiginvAdj(:,:,:)

    !> Workspace array for bond data
    real(dp), allocatable :: bondWork(:, :)

    !> Simulation time (in atomic units)
    real(dp) :: time

    !> Start time within simulation
    real(dp) :: startTime

    !> Kinetic energy
    real(dp) :: energyKin

    !> Most recently calculated bond population (retained for final write out)
    real(dp) :: lastBondPopul

    !> IDs for Population output files
    type(TFileDescr), allocatable :: populDat(:)

    !> Dipole output file ID
    type(TFileDescr) :: dipoleDat

    !> Charge output file ID
    type(TFileDescr) :: qDat

    !> Energy output file ID
    type(TFileDescr) :: energyDat

    !> Atom-resolved energy output file ID
    type(TFileDescr) :: atomEnergyDat

    !> Forces output file ID
    type(TFileDescr) :: forceDat

    !> Coords  output file ID
    type(TFileDescr) :: coorDat

    !> Pairwise bond population output file ID
    type(TFileDescr) :: fdBondPopul

    !> Pairwise bond energy output file ID
    type(TFileDescr) :: fdBondEnergy

    !> Potential acting on the system
    type(TPotentials) :: potential

    !> Count of the number of times dynamics has been initialised
    integer :: nDynamicsInit = 0

    !> Number of times doDynamics() has been called, for example with different kick directions
    integer, public :: iCall

    !> Has the propagator been initialised? (Needed for the API)
    logical, public :: tPropagatorsInitialized = .false.

    !> Has the field been initialised? (Needed for the API)
    logical, public :: tdFieldIsSet = .false.

    !> Is the field set in the API
    logical, public :: tdFieldThroughAPI = .false.

    !> Are coordinates and velocities set
    logical, public :: tdCoordsAndVelosAreSet = .false.

    !> Were these set in the API
    logical, public :: tdCoordsAndVelosThroughAPI = .false.

    !> Are ion dynamics required
    logical, public :: tIons

    !> Updated coordinates after dynamics step
    real(dp), allocatable, public :: coordNew(:,:)

    !> Velocities of moving atoms
    real(dp), allocatable, public :: movedVelo(:,:)

    !> Index for the atoms moving
    integer, allocatable, public :: indMovedAtom(:)

    !> Data structure for the energy components
    type(TEnergies), public :: energy

    !> Dipole moment
    real(dp), allocatable, public :: dipole(:,:)

    !> Forces on atoms
    real(dp), allocatable, public :: totalForce(:,:)

    !> Occupations in canonical orbitals
    real(dp), allocatable, public :: occ(:)

    !> Negative gross atomic charges
    real(dp), allocatable, public :: deltaQ(:,:)

    !> Electric field at current step
    real(dp), public :: presentField(3)

    !> Time step
    real(dp), public :: dt

    !> Number of dynamics steps to perform
    integer, public :: nSteps
    !> Dense matrix descriptor
    type(TDenseDescr) :: denseDesc

  end type TElecDynamics


  !> Nature of perturbing field applied to the initial state
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


  !> Envelope function for pulses and other fields
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


  !> Spin type for calculation
  type :: TDSpinTypesEnum

    !> Only singlet excitations (no change of total spin)
    integer :: singlet = 1

    !> Only triplet excitations (with change of total spin = 1)
    integer :: triplet = 2

  end type TDSpinTypesEnum


  !> Container for enumerated types of spin polarized spectra
  type(TDSpinTypesEnum), parameter :: tdSpinTypes = TDSpinTypesEnum()


  !> Prefix for dump files for restart
  character(*), parameter :: restartFileName = 'tddump'


  !> Prefix for dump files for pump-probe
  character(*), parameter :: pumpFilesDir = 'pump_frames'


contains


  !> Initialisation of input variables
  subroutine TElecDynamics_init(this, inp, species, speciesName, tWriteAutotest, autotestTag,&
      & randomThermostat, mass, nAtom, skCutoff, mCutoff, atomEigVal, dispersion, nonSccDeriv,&
      & tPeriodic, parallelKS, tRealHS, kPoint, kWeight, isHybridXc, sccCalc, tblite,&
      & eFieldScaling, hamiltonianType, errStatus, denseDesc)

    !> ElecDynamics instance
    type(TElecDynamics), intent(out) :: this

    !> ElecDynamicsInp instance
    type(TElecDynamicsInp), intent(in) :: inp

    !> Label for each atomic chemical species
    character(mc), allocatable, intent(in) :: speciesName(:)

    !> Produce tagged output?
    logical, intent(in) :: tWriteAutotest

    !> Tagged output files (machine readable)
    character(*), intent(in) :: autotestTag

    !> Self energy (orbital, atom)
    real(dp), intent(in), allocatable :: atomEigVal(:,:)

    !> Nr. of atoms
    integer, intent(in) :: nAtom

    ! thermostat object
    type(TRanlux), allocatable, intent(inout) :: randomThermostat

    !> Longest range of interactions for which neighbours are required
    real(dp), intent(in) :: mCutoff

    !> Cut off distance for Slater-Koster interactions
    real(dp) :: skCutoff

    !> List of atomic masses
    real(dp) :: mass(:)

    !> Dispersion data and calculations
    class(TDispersionIface), allocatable, intent(inout) :: dispersion

    !> Differentiation method for (H^0,S)
    type(TNonSccDiff), intent(in) :: nonSccDeriv

    !> Types of the atoms (nAtom)
    integer, intent(in) :: species(:)

    !> If calculation is periodic
    logical, intent(in) :: tPeriodic

    !> Dummy thermostat object
    type(TDummyThermostat), allocatable :: pDummyTherm

    !> MD Framework
    type(TMDCommon), allocatable :: pMDFrame

    !> Contains (iK, iS) tuples to be processed in parallel by various processor groups
    type(TParallelKS), intent(in) :: parallelKS

    !> H and S are real
    logical, intent(in) :: tRealHS

    !> The k-points
    real(dp) :: kPoint(:,:)

    !> Weight of the K-Points
    real(dp) :: KWeight(:)

    !> Whether this is a hybrid xc-functional calculation
    logical, intent(in) :: isHybridXc

    !> SCC module internal variables
    type(TScc), intent(in), allocatable :: sccCalc

    !> Library interface handler
    type(TTBLite), intent(in), allocatable :: tblite

    !> Any dielectric environment scaling
    class(TScaleExtEField), intent(in) :: eFieldScaling

    !> Type of Hamiltonian used
    integer, intent(in) :: hamiltonianType

    !> Error status
    type(TStatus), intent(out) :: errStatus

    !> Dense matrix descriptor
    type(TDenseDescr), intent(in) :: denseDesc

    real(dp) :: norm, tempAtom
    logical :: tMDstill
    integer :: iAtom

    this%field = eFieldScaling%scaledExtEField(inp%tdField)
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
    this%hamiltonianType = hamiltonianType
    allocate(this%parallelKS, source=parallelKS)
    allocate(this%populDat(this%parallelKS%nLocalKS))
    if (.not.any([allocated(sccCalc), allocated(tblite)])) then
      @:RAISE_ERROR(errStatus, -1, "SCC calculations are currently required for dynamics")
    end if
    if (allocated(sccCalc)) then
      this%sccCalc = sccCalc
    end if
    if (allocated(tblite)) then
      this%tblite = tblite
    end if
    ! we assume denseDesc is allocated
    this%denseDesc = denseDesc

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
      @:RAISE_ERROR(errStatus, -1, "Wrong type of perturbation.")
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
      this%fieldDir(:) = inp%reFieldPolVec + imag * inp%imFieldPolVec
      norm = sqrt(dot_product(real(this%fieldDir, dp),real(this%fieldDir, dp)))
      this%fieldDir(:) = this%fieldDir / norm
      this%tEnvFromFile = (this%envType == envTypes%fromFile)
      this%indExcitedAtom = inp%indExcitedAtom
      this%nExcitedAtom = inp%nExcitedAtom
    end if

    if (this%tKick) then
      if (inp%polDir == 4) then
        this%polDirs = [1, 2, 3]
      else
        this%polDirs = [inp%polDir]
        this%currPolDir = inp%polDir
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
      @:RAISE_ERROR(errStatus, -1, "Real hamiltonian required for bond energies")
    end if
    this%tBondP = inp%tBondP
    this%species = species
    this%tPeriodic = tPeriodic
    this%isHybridXc = isHybridXc
    this%tWriteAtomEnergies = inp%tWriteAtomEnergies

    if (this%tIons) then
      if (.not. this%tRealHS) then
        @:RAISE_ERROR(errStatus, -1, "Ion dynamics is not implemented yet for complex&
            & hamiltonians (k-points, spin-orbit, ...).")
      end if
      this%tForces = .true.
      this%indMovedAtom = inp%indMovedAtom
      this%nMovedAtom = inp%nMovedAtom
      tempAtom = inp%tempAtom
      tMDstill = .false.

      allocate(this%movedVelo(3, this%nMovedAtom))
      allocate(this%movedMass(3, this%nMovedAtom))
      this%movedMass(:,:) = spread(mass(this%indMovedAtom), 1, 3)

      allocate(this%initialVelocities(3, this%nMovedAtom))
      this%ReadMDVelocities = allocated(inp%initialVelocities)
      if (this%ReadMDVelocities) then
        this%initialVelocities(:,:) = inp%initialVelocities
        this%movedVelo(:, :) = this%initialVelocities
      else
        this%movedVelo(:, :) = 0.0_dp
      end if

      allocate(this%coordNew(3, nAtom))
      allocate(this%movedAccel(3, this%nMovedAtom))

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
        this%nExcitedAtom = nAtom
      end if
      if (allocated(tblite)) then
        @:RAISE_ERROR(errStatus, -1, "Ion dynamics and forces not available for xTB Hamiltonian")
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
    if (allocated(atomEigVal)) then
      allocate(this%atomEigVal, source=atomEigVal)
    end if

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


  !> Driver for time dependent propagation to calculate either a spectrum or laser driving
  subroutine runDynamics(this, boundaryCond, eigvecs, H0, speciesAll, q0, referenceN0, ints,&
      & filling, neighbourList, nNeighbourSK, iSquare, iSparseStart, img2CentCell, orb, coord,&
      & spinW, repulsive, env, tDualSpinOrbit, xi, thirdOrd, solvation, eFieldScaling, hybridXc,&
      & qDepExtPot, dftbU, iAtInCentralRegion, tFixEf, Ef, coordAll, onSiteElements, skHamCont,&
      & skOverCont, latVec, invLatVec, iCellVec, rCellVec, cellVec, electronicSolver, eigvecsCplx,&
      & taggedWriter, refExtPot, errStatus)

    !> ElecDynamics instance
    type(TElecDynamics) :: this

    !> Boundary conditions on the calculation
    type(TBoundaryConditions), intent(in) :: boundaryCond

    !> Real Eigenvectors
    real(dp), intent(inout), allocatable :: eigvecs(:,:,:)

    !> Complex Eigevenctors
    complex(dp), intent(inout), allocatable :: eigvecsCplx(:,:,:)

    !> Sparse non-SCC hamiltonian
    real(dp), intent(in) :: H0(:)

    !> Species of all atoms in the system
    integer, intent(in) :: speciesAll(:)

    !> Reference atomic occupations
    real(dp), intent(inout) :: q0(:,:,:)

    !> Reference charges from the Slater-Koster file
    real(dp), intent(in) :: referenceN0(:,:)

    !> Integral container
    type(TIntegral), intent(inout) :: ints

    !> Central atomic coordinates
    real(dp), allocatable, intent(inout) :: coord(:,:)

    !> All atomic coordinates
    real(dp), allocatable, intent(inout) :: coordAll(:,:)

    !> Spin constants
    real(dp), allocatable, intent(in) :: spinW(:,:,:)

    !> Occupations
    real(dp), intent(inout) :: filling(:,:,:)

    !> Number of neighbours for each of the atoms
    integer, intent(inout) :: nNeighbourSK(:)

    !> Index array for location of atomic blocks in large sparse arrays
    integer, allocatable, intent(inout) :: iSparseStart(:,:)

    !> Image atoms to their equivalent in the central cell
    integer, allocatable, intent(inout) :: img2CentCell(:)

    !> Index array for start of atomic block in dense matrices
    integer, intent(in) :: iSquare(:)

    !> List of neighbours for each atom
    type(TNeighbourList), intent(inout) :: neighbourList

    !> Repulsive information
    class(TRepulsive), allocatable, intent(inout) :: repulsive

    !> Atomic orbital information
    type(TOrbitals), intent(in) :: orb

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

    !> Any dielectric environment scaling
    type(TScaleExtEField), intent(in) :: eFieldScaling

    !> Range separation contributions
    class(THybridXcFunc), allocatable, intent(inout) :: hybridXc

    !> Proxy for querying Q-dependant external potentials
    type(TQDepExtPotProxy), intent(inout), allocatable :: qDepExtPot

    !> DFTB+U functional (if used)
    type(TDftbU), intent(in), allocatable :: dftbU

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

    !> Cell vectors in absolute units
    real(dp), intent(in) :: rCellVec(:,:)

    !> Cell vectors (in units of the lattice constants) to cells of the lattice
    real(dp), intent(in) :: cellVec(:,:)

    !> Index of cell in cellVec and rCellVec for each atom
    integer, allocatable, intent(in) :: iCellVec(:)

    !> Electronic solver information
    type(TElectronicSolver), intent(inout) :: electronicSolver

    !> Tagged writer object
    type(TTaggedWriter), intent(inout) :: taggedWriter

    !> Reference external potential (usual provided via API)
    type(TRefExtPot) :: refExtPot

    !> Error status
    type(TStatus), intent(inout) :: errStatus

    integer :: iPol
    logical :: tWriteAutotest

    tWriteAutotest = this%tWriteAutotest
    this%iCall = 1
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
        call doDynamics(this, boundaryCond, eigvecs, H0, q0, referenceN0, ints, filling,&
            & neighbourList, nNeighbourSK, iSquare, iSparseStart, img2CentCell, orb, coord, spinW,&
            & repulsive, env, tDualSpinOrbit, xi, thirdOrd, solvation, eFieldScaling, hybridXc,&
            & qDepExtPot, dftbU, iAtInCentralRegion, tFixEf, Ef, tWriteAutotest, coordAll,&
            & onSiteElements, skHamCont, skOverCont, electronicSolver, speciesAll, eigvecsCplx,&
            & taggedWriter, refExtPot, latVec, invLatVec, iCellVec, rCellVec, cellVec, errStatus)
        this%iCall = this%iCall + 1
      end do
    else
      call doDynamics(this, boundaryCond, eigvecs, H0, q0, referenceN0, ints, filling,&
          & neighbourList, nNeighbourSK, iSquare, iSparseStart, img2CentCell, orb, coord, spinW,&
          & repulsive, env, tDualSpinOrbit, xi, thirdOrd, solvation, eFieldScaling, hybridXc,&
          & qDepExtPot, dftbU, iAtInCentralRegion, tFixEf, Ef, tWriteAutotest, coordAll,&
          & onSiteElements, skHamCont, skOverCont, electronicSolver, speciesAll, eigvecsCplx,&
          & taggedWriter, refExtPot, latVec, invLatVec, iCellVec, rCellVec, cellVec, errStatus)
    end if

  end subroutine runDynamics


  !> Runs the electronic dynamics of the system
  subroutine doDynamics(this, boundaryCond, eigvecsReal, H0, q0, referenceN0, ints, filling,&
      & neighbourList, nNeighbourSK, iSquare, iSparseStart, img2CentCell, orb, coord, spinW,&
      & repulsive, env, tDualSpinOrbit, xi, thirdOrd, solvation, eFieldScaling, hybridXc,&
      & qDepExtPot, dftbU, iAtInCentralRegion, tFixEf, Ef, tWriteAutotest, coordAll,&
      & onSiteElements, skHamCont, skOverCont, electronicSolver, speciesAll, eigvecsCplx,&
      & taggedWriter, refExtPot, latVec, invLatVec, iCellVec, rCellVec, cellVec, errStatus)

    !> ElecDynamics instance
    type(TElecDynamics) :: this

    !> Boundary conditions on the calculation
    type(TBoundaryConditions), intent(in) :: boundaryCond

    !> Real Eigenvectors
    real(dp), intent(inout), allocatable :: eigvecsReal(:,:,:)

    !> Complex Eigevenctors
    complex(dp), intent(inout), allocatable :: eigvecsCplx(:,:,:)

    !> Sparse storage for non-SCC hamiltonian
    real(dp), intent(in) :: H0(:)

    !> Reference atomic occupations
    real(dp), intent(inout) :: q0(:,:,:)

    !> Reference charges from the Slater-Koster file
    real(dp), intent(in) :: referenceN0(:,:)

    !> Integral container
    type(TIntegral), intent(inout) :: ints

    !> Atomic coordinates
    real(dp), allocatable, intent(inout) :: coord(:,:)

    !> All atomic coordinates
    real(dp), allocatable, intent(inout) :: coordAll(:,:)

    !> Spin constants
    real(dp), allocatable, intent(in) :: spinW(:,:,:)

    !> Occupations
    real(dp), intent(inout) :: filling(:,:,:)

    !> Number of neighbours for each of the atoms
    integer, intent(inout) :: nNeighbourSK(:)

    !> Index array for location of atomic blocks in large sparse arrays
    integer, allocatable, intent(inout) :: iSparseStart(:,:)

    !> Image atoms to their equivalent in the central cell
    integer, allocatable, intent(inout) :: img2CentCell(:)

    !> Index array for start of atomic block in dense matrices
    integer, intent(in) :: iSquare(:)

    !> List of neighbours for each atom
    type(TNeighbourList), intent(inout) :: neighbourList

    !> Repulsive information
    class(TRepulsive), allocatable, intent(inout) :: repulsive

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

    !> DFTB+U functional (if used)
    type(TDftbU), intent(in), allocatable :: dftbU

    !> Spin orbit constants if required
    real(dp), allocatable, intent(in) :: xi(:,:)

    !> 3rd order settings
    type(TThirdOrder), intent(inout), allocatable :: thirdOrd

    !> Solvation model
    class(TSolvation), allocatable, intent(inout) :: solvation

    !> Any dielectric environment scaling
    class(TScaleExtEField), intent(in) :: eFieldScaling

    !> Range separation contributions
    class(THybridXcFunc), allocatable, intent(inout) :: hybridXc

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

    !> Electronic solver information
    type(TElectronicSolver), intent(inout) :: electronicSolver

    !> Tagged writer object
    type(TTaggedWriter), intent(inout) :: taggedWriter

    !> Reference external potential (usual provided via API)
    type(TRefExtPot) :: refExtPot

    !> Lattice vectors if periodic
    real(dp), intent(in) :: latVec(:,:)

    !> Inverse of the lattice vectors
    real(dp), intent(in) :: invLatVec(:,:)

    !> Cell vectors in absolute units
    real(dp), intent(in) :: rCellVec(:,:)

    !> Vectors (in units of the lattice constants) to cells of the lattice
    real(dp), intent(in) :: cellVec(:,:)

    !> Index of cell in cellVec and rCellVec for each atom
    integer, allocatable, intent(in) :: iCellVec(:)

    !> Species of all atoms in the system
    integer, intent(in) :: speciesAll(:)

    !> Error status
    type(TStatus), intent(inout) :: errStatus

    type(TTimer) :: loopTime
    integer :: iStep, i, j
    real(dp) :: timeElec

    call env%globalTimer%startTimer(globalTimers%elecDynInit)

    call initializeDynamics(this, boundaryCond, coord, orb, neighbourList, nNeighbourSK,&
       & iSquare, iSparseStart, img2CentCell, skHamCont, skOverCont, ints, env, coordAll,&
       & H0, spinW, tDualSpinOrbit, xi, thirdOrd, dftbU, onSiteElements,&
       & refExtPot, solvation, eFieldScaling, hybridXc, referenceN0, q0, repulsive,&
       & iAtInCentralRegion, eigvecsReal, eigvecsCplx, filling, qDepExtPot, tFixEf, Ef, latVec,&
       & invLatVec, iCellVec, rCellVec, cellVec, speciesAll, errStatus)
    @:PROPAGATE_ERROR(errStatus)

    call env%globalTimer%stopTimer(globalTimers%elecDynInit)

    call env%globalTimer%startTimer(globalTimers%elecDynLoop)
    call loopTime%start()

    write(stdOut, "(A)")
    write(stdOut, "(A)") 'Starting electronic dynamics...'
    write(stdOut, "(A80)") repeat("-", 80)

    ! Main loop
    do iStep = 1, this%nSteps

      call doTdStep(this, boundaryCond, iStep, coord, orb, neighbourList, nNeighbourSK,&
       & iSquare, iSparseStart, img2CentCell, skHamCont, skOverCont, ints, env,&
       & coordAll, q0, referenceN0, spinW, tDualSpinOrbit, xi, thirdOrd, dftbU,&
       & onSiteElements, refExtPot, solvation, eFieldScaling, hybridXc, repulsive,&
       & iAtInCentralRegion, tFixEf, Ef, electronicSolver, qDepExtPot, errStatus)
      @:PROPAGATE_ERROR(errStatus)

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
      call writeTDAutotest(this, this%dipole, this%energy, this%deltaQ, coord, this%totalForce,&
          & this%occ, this%lastBondPopul, taggedWriter)
    end if

  !  call finalizeDynamics(this)

  end subroutine doDynamics


  !> Updates the hamiltonian with SCC and external TD field (if any) contributions
  subroutine updateH(this, H1, ints, H0, speciesAll, qq, q0, coord, orb, potential,&
      & neighbourList, nNeighbourSK, iSquare, iSparseStart, img2CentCell, iStep, chargePerShell,&
      & spinW, env, tDualSpinOrbit, xi, thirdOrd, qBlock, dftbU, onSiteElements, refExtPot,&
      & deltaRho, H1LC, Ssqr, solvation, hybridXc, dispersion, rho, errStatus)

    !> ElecDynamics instance
    type(TElecDynamics) :: this

    !> Square hamiltonian at each spin and k-point
    complex(dp), intent(inout) :: H1(:,:,:)

    !> Integral container
    type(TIntegral), intent(inout) :: ints

    !> Sparse storage for non-SCC hamiltonian
    real(dp), intent(in) :: H0(:)

    !> Species of all atoms in the system
    integer, intent(in) :: speciesAll(:)

    !> Atomic populations
    real(dp), intent(inout) :: qq(:,:,:)

    !> Reference atomic occupations
    real(dp), intent(inout) :: q0(:,:,:)

    !> Atomic coordinates
    real(dp), allocatable, intent(inout) :: coord(:,:)

    !> Atomic orbital information
    type(TOrbitals), intent(in) :: orb

    !> Potential acting on the system
    type(TPotentials), intent(inout) :: potential

    !> List of neighbours for each atom
    type(TNeighbourList), intent(inout) :: neighbourList

    !> Number of neighbours for each of the atoms
    integer, intent(in) :: nNeighbourSK(:)

    !> Index array for start of atomic block in dense matrices
    integer, intent(in) :: iSquare(:)

    !> Index array for location of atomic blocks in large sparse arrays
    integer, intent(in) :: iSparseStart(0:,:)

    !> Image atoms to their equivalent in the central cell
    integer, intent(in) :: img2CentCell(:)

    !> Current step of the propagation
    integer, intent(in) :: iStep

    !> Electrons in each atomic shell
    real(dp), intent(inout) :: chargePerShell(:,:,:)

    !> Spin constants
    real(dp), allocatable, intent(in) :: spinW(:,:,:)

    !> Environment settings
    type(TEnvironment), intent(inout) :: env

    !> Is dual spin orbit being used (block potentials)
    logical, intent(in) :: tDualSpinOrbit

    !> Spin orbit constants if required
    real(dp), allocatable, intent(in) :: xi(:,:)

    !> 3rd order settings
    type(TThirdOrder), intent(inout), allocatable :: thirdOrd

    !> Block (dual) atomic populations
    real(dp), intent(inout), allocatable :: qBlock(:,:,:,:)

    !> DFTB+U functional (if used)
    type(TDftbU), intent(in), allocatable :: dftbU

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
    class(THybridXcFunc), allocatable, intent(inout) :: hybridXc

    !> Dispersion data and calculations
    class(TDispersionIface), allocatable, intent(inout) :: dispersion

    !> Density matrix
    complex(dp), intent(in) :: rho(:,:,:)

    !> Error status
    type(TStatus), intent(out) :: errStatus

    real(dp), allocatable :: qiBlock(:,:,:,:) ! not allocated since no imaginary ham
    real(dp), allocatable :: iHam(:,:) ! not allocated since no imaginary ham
    real(dp), allocatable :: T2(:,:)
    integer :: iAtom, iEatom, iSpin, iKS, iK
    logical :: tImHam
    integer :: nLocalCols, nLocalRows

#:if WITH_SCALAPACK
    nLocalRows = size(H1, dim=1)
    nLocalCols = size(H1, dim=2)
#:else
    nLocalRows = this%nOrbs
    nLocalCols = this%nOrbs
#:endif

    if (this%tRealHS) then
      allocate(T2(nLocalRows,nLocalCols))
    end if

    ints%hamiltonian(:,:) = 0.0_dp

    if (this%nSpin == 2) then
      call ud2qm(qq)
      call ud2qm(q0)
    end if

    tImHam = .false. ! for the moment

    call resetExternalPotentials(refExtPot, potential)
    call resetInternalPotentials(tDualSpinOrbit, xi, orb, speciesAll, potential)

    call getChargePerShell(qq, orb, speciesAll, chargePerShell)
    call addChargePotentials(env, this%sccCalc, this%tblite, .true., qq, q0, chargePerShell,&
        & orb, this%multipole, speciesAll, neighbourList, img2CentCell, spinW, solvation,&
        & thirdOrd, dispersion, potential)

    if (allocated(dftbU) .or. allocated(onSiteElements)) then
      ! convert to qm representation
      call ud2qm(qBlock)
    end if

    if (allocated(dftbU)) then
      call addBlockChargePotentials(qBlock, qiBlock, dftbU, .false., speciesAll, orb, potential)
    end if
    if (allocated(onSiteElements)) then
      call addOnsShift(potential%intBlock, potential%iOrbitalBlock, qBlock, qiBlock,&
          & onSiteElements, speciesAll, orb, q0)
    end if

    ! Add time dependent field if necessary
    if (this%tLaser) then
      call setPresentField(this, iStep, errStatus)
      @:PROPAGATE_ERROR(errStatus)
      do iAtom = 1, this%nExcitedAtom
        iEatom = this%indExcitedAtom(iAtom)
        potential%extAtom(iEatom, 1) = dot_product(coord(:,iEatom), this%presentField)
      end do
      call totalShift(potential%extShell, potential%extAtom, orb, speciesAll)
      call totalShift(potential%extBlock, potential%extShell, orb, speciesAll)
    end if

    potential%intBlock = potential%intBlock + potential%extBlock
    potential%intShell = potential%intShell + potential%extShell

    call getSccHamiltonian(env, H0, ints, nNeighbourSK, neighbourList, speciesAll, orb,&
        & iSparseStart, img2CentCell, potential, .false., ints%hamiltonian, iHam)

    ! Hack due to not using Pauli-type structure outside of this part of the routine
    if (this%nSpin == 2) then
      ints%hamiltonian = 2.0_dp * ints%hamiltonian
      call qm2ud(ints%hamiltonian)
      call qm2ud(q0)
      call qm2ud(qq)
    end if

#:if WITH_SCALAPACK
    do iKS = 1, this%parallelKS%nLocalKS
      iK = this%parallelKS%localKS(1, iKS)
      iSpin = this%parallelKS%localKS(2, iKS)
      if (this%tRealHS) then
        call  unpackHSRealBlacs(env%blacs, ints%hamiltonian(:,iSpin), neighbourList%iNeighbour,&
               & nNeighbourSK, iSparseStart, img2CentCell, this%denseDesc, T2)
        H1(:,:,iSpin) = cmplx(T2, 0.0_dp, dp)
      ! TODO: add here the unpacking of H1 for kpoints
      end if
    end do
#:else
    do iKS = 1, this%parallelKS%nLocalKS
      iK = this%parallelKS%localKS(1, iKS)
      iSpin = this%parallelKS%localKS(2, iKS)
      if (this%tRealHS) then
        call unpackHS(T2, ints%hamiltonian(:,iSpin), neighbourList%iNeighbour, nNeighbourSK,&
            & iSquare, iSparseStart, img2CentCell)
        call adjointLowerTriangle(T2)
        H1(:,:,iSpin) = cmplx(T2, 0.0_dp, dp)
      else
        call unpackHS(H1(:,:,iKS), ints%hamiltonian(:,iSpin), this%kPoint(:,iK),&
            & neighbourList%iNeighbour, nNeighbourSK, this%iCellVec, this%cellVec, iSquare,&
            & iSparseStart, img2CentCell)
        call adjointLowerTriangle(H1(:,:,iKS))
      end if
    end do
#:endif

    ! add hybrid xc-functional contribution
    if (this%isHybridXc) then
      deltaRho = rho
      if (this%nSpin > 2) then
        @:RAISE_ERROR(errStatus, -1, "HybridXc: Not implemented for non-colinear spin.")
      end if
      ! denseSubtractDensityOfAtoms_cmplx_nonperiodic
      call denseSubtractDensityOfAtoms(q0, iSquare, deltaRho)
      do iSpin = 1, this%nSpin
        H1LC(:,:) = (0.0_dp, 0.0_dp)
        call hybridXc%addCamHamiltonianMatrix_cluster_cmplx(iSquare, sSqr(:,:, iSpin),&
            & deltaRho(:,:, iSpin), H1LC)
        call adjointLowerTriangle(H1LC)
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

    !> Atomic coordinates
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


  !> Creates array for an external TD field
  subroutine getTDFunction(this, startTime)

    !> ElecDynamics instance
    type(TElecDynamics), intent(inout) :: this

    !> Starting time of the simulation, if relevant
    real(dp), intent(in) :: startTime

    real(dp) :: midPulse, deltaT, angFreq, E0, time, envelope
    real(dp) :: tdfun(3)
    integer :: iStep
    type(TFileDescr) :: laserDat

    allocate(this%tdFunction(3, 0:this%nSteps))
    this%tdFunction(:,:) = 0.0_dp

    midPulse = (this%time0 + this%time1)/2.0_dp
    deltaT = this%time1 - this%time0
    angFreq = this%omega
    E0 = this%field
    if (this%tKickAndLaser) then
      E0 = this%laserField
    end if
    if (this%tEnvFromFile) then
      E0 = 0.0_dp !this is to make sure we never sum the current field with that read from file
    end if

    if (this%tEnvFromFile) then
      call openFile(laserDat, "laser.dat", mode="r")
    else
      call openOutputFile(this, laserDat, 'laser.dat')
    end if

    if (.not. this%tEnvFromFile) then
      write(laserDat%unit, "(A)") "#     time (fs)  |  E_x (eV/ang)  | E_y (eV/ang) | E_z (eV/ang)"
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
        read(laserDat%unit, *)time, tdfun(1), tdfun(2), tdfun(3)
        this%tdFunction(:, iStep) = tdfun * (Bohr__AA / Hartree__eV)
      else
        this%tdFunction(:, iStep) = E0 * envelope * aimag(exp(imag*(time*angFreq + this%phase))&
            & * this%fieldDir)
        write(laserDat%unit, "(5F15.8)") time * au__fs,&
            & this%tdFunction(:, iStep) * (Hartree__eV / Bohr__AA)
      end if

    end do

    call closeFile(laserDat)

  end subroutine getTDFunction


  !> Calculate charges, dipole moments
  subroutine getChargeDipole(this, deltaQ, qq, multipole, dipole, q0, rho, Ssqr, Dsqr, Qsqr,&
      & coord, iSquare, eFieldScaling, qBlock, qNetAtom, errStatus, &
      & iNeighbour, nNeighbourSK, orb, iSparseStart, img2CentCell, env, ints)

    !> ElecDynamics instance
    type(TElecDynamics), intent(inout) :: this

    !> Negative gross charge
    real(dp), intent(out) :: deltaQ(:,:)

    !> Dipole moment
    real(dp), intent(out) :: dipole(:,:)

    !> Atomic occupations
    real(dp), intent(out) :: qq(:,:,:)

    !> Reference atomic occupations
    real(dp), intent(in) :: q0(:,:,:)

    !> Multipole moments
    type(TMultipole), intent(inout) :: multipole

    !> Atomic coordinates
    real(dp), intent(in) :: coord(:,:)

    !> Density matrix
    complex(dp), intent(in) :: rho(:,:,:)

    !> Square overlap matrix
    complex(dp), intent(in) :: Ssqr(:,:,:)

    !> Square dipole matrix
    complex(dp), intent(in), optional :: Dsqr(:,:,:,:)

    !> Square quadrupole matrix
    complex(dp), intent(in), optional :: Qsqr(:,:,:,:)

    !> Index array for start of atomic block in dense matrices
    integer, intent(in) :: iSquare(:)

    !> Any dielectric environment scaling
    class(TScaleExtEField), intent(in) :: eFieldScaling

    !> Mulliken block charges
    real(dp), allocatable, intent(inout) :: qBlock(:,:,:,:)

    !> Net (on-site only) atomic charge
    real(dp), allocatable, intent(inout) :: qNetAtom(:)

    !> Error status
    type(TStatus), intent(out) :: errStatus

    !> Atomic neighbour data
    integer, intent(in) :: iNeighbour(0:,:)

    !> Number of neighbours for each of the atoms
    integer, intent(in) :: nNeighbourSK(:)

    !> Atomic orbital information
    type(TOrbitals), intent(in) :: orb

    !> index array for location of atomic blocks in large sparse arrays
    integer, intent(in) :: iSparseStart(0:,:)

    !> image atoms to their equivalent in the central cell
    integer, intent(in) :: img2CentCell(:)

    !> Environment settings
    type(TEnvironment), intent(inout) :: env

    !> Integral container
    type(TIntegral), intent(inout) :: ints

    integer :: iAt, iSpin, iOrb1, iOrb2, nOrb, iKS, iK, ii
    real(dp), allocatable :: tmp(:,:)

    qq(:,:,:) = 0.0_dp

    if (this%tRealHS) then

#:if WITH_SCALAPACK
      this%rhoPrim(:,:) = 0.0_dp
      allocate(tmp (size(rho,dim=1),size(rho,dim=2)))
      do iSpin = 1, this%nSpin
        tmp = real(rho(:,:,iSpin), dp)
        call packRhoRealBlacs(env%blacs, this%denseDesc, tmp, iNeighbour, nNeighbourSK,&
        & orb%mOrb, iSparseStart, img2CentCell, this%rhoPrim(:,iSpin))
      end do
      deallocate(tmp)
      call mpifx_allreduceip(env%mpi%globalComm, this%rhoPrim, MPI_SUM)

      do iSpin = 1, this%nSpin
       call mulliken(env, qq(:,:,iSpin), ints%overlap, this%rhoPrim(:,iSpin), orb, iNeighbour,&
        & nNeighbourSK, img2CentCell, iSparseStart)
      end do
#:else
      do iSpin = 1, this%nSpin
        do iAt = 1, this%nAtom
          iOrb1 = iSquare(iAt)
          iOrb2 = iSquare(iAt+1)-1
          ! hermitian transpose used as real part only is needed
          qq(:iOrb2-iOrb1+1,iAt,iSpin) = real(sum(&
              & rho(:,iOrb1:iOrb2,iSpin)*Ssqr(:,iOrb1:iOrb2,iSpin), dim=1), dp)
        end do
      end do
#:endif

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

    if (present(Dsqr)) then
      multipole%dipoleAtom(:, :, :) = 0.0_dp
      if (this%tRealHS) then
        do iSpin = 1, this%nSpin
          do iAt = 1, this%nAtom
            iOrb1 = iSquare(iAt)
            iOrb2 = iSquare(iAt+1)-1
            do ii = iOrb1, iOrb2
              multipole%dipoleAtom(:,iAt,iSpin) = multipole%dipoleAtom(:,iAt,iSpin) &
                  & + real(matmul(Dsqr(:,:,ii,iSpin), rho(:,ii,iSpin)), dp)
            end do
          end do
        end do
      else
        do iKS = 1, this%parallelKS%nLocalKS
          iK = this%parallelKS%localKS(1, iKS)
          iSpin = this%parallelKS%localKS(2, iKS)

          do iAt = 1, this%nAtom
            iOrb1 = iSquare(iAt)
            iOrb2 = iSquare(iAt+1)-1
            do ii = iOrb1, iOrb2
              multipole%dipoleAtom(:,iAt,iSpin) = multipole%dipoleAtom(:,iAt,iSpin) &
                & + this%kWeight(ik)*real(matmul(conjg(Dsqr(:,:,ii,iKS)), rho(:,ii,iKS)), dp)
            end do
          end do
        end do
      end if
    end if

    if (present(Qsqr)) then
      multipole%quadrupoleAtom(:, :, :) = 0.0_dp
      if (this%tRealHS) then
        do iSpin = 1, this%nSpin
          do iAt = 1, this%nAtom
            iOrb1 = iSquare(iAt)
            iOrb2 = iSquare(iAt+1)-1
            do ii = iOrb1, iOrb2
              multipole%quadrupoleAtom(:,iAt,iSpin) = multipole%quadrupoleAtom(:,iAt,iSpin) &
                  & + real(matmul(Qsqr(:,:,ii,iSpin), rho(:,ii,iSpin)), dp)
            end do
          end do
        end do
      else
        do iKS = 1, this%parallelKS%nLocalKS
          iK = this%parallelKS%localKS(1, iKS)
          iSpin = this%parallelKS%localKS(2, iKS)

          do iAt = 1, this%nAtom
            iOrb1 = iSquare(iAt)
            iOrb2 = iSquare(iAt+1)-1
            do ii = iOrb1, iOrb2
              multipole%quadrupoleAtom(:,iAt,iSpin) = multipole%quadrupoleAtom(:,iAt,iSpin) &
                & + this%kWeight(ik)*real(matmul(conjg(Qsqr(:,:,ii,iKS)), rho(:,ii,iKS)), dp)
            end do
          end do
        end do
      end if
    end if

    deltaQ(:,:) = sum((qq - q0), dim=1)
    dipole(:,:) = -matmul(coord, deltaQ)
    if (allocated(multipole%dipoleAtom)) then
      dipole(:,:) = dipole - sum(multipole%dipoleAtom(:, :, :), dim=2)
    end if

    dipole(:,1) = eFieldScaling%scaledSoluteDipole(dipole(:,1))

    if (allocated(qBlock)) then
      if (.not. this%tRealHS) then
        @:RAISE_ERROR(errStatus, -1, "Block populations not implemented yet")
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
      do iSpin = 1, this%nSpin
        do iAt = 1, this%nAtom
          iOrb1 = iSquare(iAt)
          iOrb2 = iSquare(iAt+1)
          nOrb = iOrb2 - iOrb1
          qBlock(:nOrb,:nOrb,iAt,iSpin) = 0.5_dp * (qBlock(:nOrb,:nOrb,iAt,iSpin)&
              & + transpose(qBlock(:nOrb,:nOrb,iAt,iSpin)) )
        end do
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
  subroutine getTDEnergy(this, env, energy, rhoPrim, rho, neighbourList, nNeighbourSK, orb,&
      & iSquare, iSparseStart, img2CentCell, ham0, qq, q0, potential, chargePerShell, energyKin,&
      & tDualSpinOrbit, thirdOrd, solvation, hybridXc, qDepExtPot, qBlock, dftbU, xi,&
      & iAtInCentralRegion, tFixEf, Ef, onSiteElements, errStatus)

    !> ElecDynamics instance
    type(TElecDynamics), intent(inout) :: this

    !> Environment settings
    type(TEnvironment), intent(in) :: env

    !> Data type for energy components and total
    type(TEnergies), intent(inout) :: energy

    !> Sparse density matrix
    real(dp), allocatable, intent(inout) :: rhoPrim(:,:)

    !> Density matrix
    complex(dp), intent(in) :: rho(:,:,:)

    !> Sparse storage for non-SCC hamiltonian
    real(dp), intent(in) :: ham0(:)

    !> Atomic occupations
    real(dp), intent(inout) :: qq(:,:,:)

    !> Reference atomic occupations
    real(dp), intent(in) :: q0(:,:,:)

    !> Atomic orbital information
    type(TOrbitals), intent(in) :: orb

    !> Neighbour list
    type(TNeighbourList), intent(in) :: neighbourList

    !> Number of neighbours for each of the atoms
    integer, intent(in) :: nNeighbourSK(:)

    !> Index array for start of atomic block in dense matrices
    integer, intent(in) :: iSquare(:)

    !> Index array for location of atomic blocks in large sparse arrays
    integer, intent(in) :: iSparseStart(0:,:)

    !> Image atoms to their equivalent in the central cell
    integer, intent(in) :: img2CentCell(:)

    !> Electrons in each atomic shell
    real(dp), intent(in) :: chargePerShell(:,:,:)

    !> Potential acting on the system
    type(TPotentials), intent(in) :: potential

    !> Kinetic energy
    real(dp), intent(out) :: energyKin

    !> Is dual spin orbit being used
    logical, intent(in) :: tDualSpinOrbit

    !> 3rd order settings
    type(TThirdOrder), intent(inout), allocatable :: thirdOrd

    !> Solvation model
    class(TSolvation), allocatable, intent(inout) :: solvation

    !> Range separation contributions
    class(THybridXcFunc), allocatable, intent(inout) :: hybridXc

    !> Proxy for querying Q-dependant external potentials
    type(TQDepExtPotProxy), intent(inout), allocatable :: qDepExtPot

    !> Block (dual) atomic populations
    real(dp), intent(in), allocatable :: qBlock(:,:,:,:)

    !> DFTB+U functional (if used)
    type(TDftbU), intent(in), allocatable :: dftbU

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

    !> Error status
    type(TStatus), intent(inout) :: errStatus

    real(dp), allocatable :: qiBlock(:,:,:,:), tmp(:,:)
    integer :: iKS, iK, iSpin
    real(dp) :: TS(this%nSpin)
    type(TReksCalc), allocatable :: reks ! never allocated

    ! if Forces are calculated, rhoPrim has already been calculated
    ! check allways that calcEnergy is called AFTER getForces
    if (.not. this%tForces) then
      rhoPrim(:,:) = 0.0_dp

#:if WITH_SCALAPACK
      do iKS = 1, this%parallelKS%nLocalKS
        iSpin = this%parallelKS%localKS(2, iKS)
        if (this%tRealHS) then
          allocate(tmp (size(rho,dim=1),size(rho,dim=2)))
            tmp = real(rho(:,:,iSpin), dp)
            call packRhoRealBlacs(env%blacs, this%denseDesc, tmp, neighbourlist%iNeighbour, nNeighbourSK,&
            & orb%mOrb, iSparseStart, img2CentCell, rhoPrim(:,iSpin))
          deallocate(tmp)
        ! TODO: add here the case for complex Hamiltonian
        end if
      end do
      call mpifx_allreduceip(env%mpi%globalComm, rhoPrim, MPI_SUM)

#:else
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
#:endif

    end if
    call ud2qm(rhoPrim)

    TS = 0.0_dp
    call calcEnergies(env, this%sccCalc, this%tblite, qq, q0, chargePerShell, this%multipole,&
        & this%speciesAll, this%tLaser, .false., dftbU, tDualSpinOrbit, rhoPrim, ham0, orb,&
        & neighbourList, nNeighbourSK, img2CentCell, iSparseStart, 0.0_dp, 0.0_dp, TS,&
        & potential, energy, thirdOrd, solvation, hybridXc, reks, qDepExtPot, qBlock,&
        & qiBlock, xi, iAtInCentralRegion, tFixEf, Ef, .true., onSiteElements, errStatus)
    @:PROPAGATE_ERROR(errStatus)
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
  subroutine initializeTDVariables(this, rho, H1, Ssqr, Sinv, H0, ham0, Dsqr, Qsqr, ints,&
      & eigvecsReal, filling, orb, rhoPrim, potential, iNeighbour, nNeighbourSK, iSquare,&
      & iSparseStart, img2CentCell, Eiginv, EiginvAdj, energy, ErhoPrim, qBlock, qNetAtom, isDftbU,&
      & onSiteElements, eigvecsCplx, H1LC, bondWork, fdBondEnergy, fdBondPopul, lastBondPopul, time, &
      & env, errStatus)

    !> ElecDynamics instance
    type(TElecDynamics), intent(inout) :: this

    !> Real Eigenvectors
    real(dp), intent(inout), allocatable :: eigvecsReal(:,:,:)

    !> Complex Eigevenctors
    complex(dp), intent(inout), allocatable :: eigvecsCplx(:,:,:)

    !> Integral container
    type(TIntegral), intent(inout) :: ints

    !> Occupations
    real(dp), intent(inout) :: filling(:,:,:)

    !> Atomic neighbour data
    integer, intent(in) :: iNeighbour(0:,:)

    !> Number of neighbours for each of the atoms
    integer, intent(in) :: nNeighbourSK(:)

    !> Index array for start of atomic block in dense matrices
    integer, intent(in) :: iSquare(:)

    !> Index array for location of atomic blocks in large sparse arrays
    integer, intent(in) :: iSparseStart(0:,:)

    !> Image atoms to their equivalent in the central cell
    integer, intent(in) :: img2CentCell(:)

    !> Atomic orbital information
    type(TOrbitals), intent(in) :: orb

    !> Potential acting on the system
    type(TPotentials), intent(out) :: potential

    !> Data type for energy components and total
    type(TEnergies), intent(out) :: energy

    !> Sparse density matrix
    real(dp), allocatable, intent(out) :: rhoPrim(:,:)

    !> Square overlap matrix
    complex(dp), intent(inout) :: Ssqr(:,:,:)

    !> Square overlap inverse
    complex(dp), intent(inout) :: Sinv(:,:,:)

    !> Square dipole matrix
    complex(dp), intent(inout), optional :: Dsqr(:,:,:,:)

    !> Square quadrupole matrix
    complex(dp), intent(inout), optional :: Qsqr(:,:,:,:)

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

    !> Energy weighted density matrix
    real(dp), allocatable, intent(out) :: ErhoPrim(:)

    !> Block (dual) atomic populations
    real(dp), intent(inout), allocatable :: qBlock(:,:,:,:)

    !> Net (onsite only) atomic charges
    real(dp), intent(inout), allocatable :: qNetAtom(:)

    !> Is this a DFTB+U calculation
    logical, intent(in) :: isDftbU

    !> Corrections terms for on-site elements
    real(dp), intent(in), allocatable :: onSiteElements(:,:,:,:)

    !> LC contribution to hamiltonian
    complex(dp), allocatable, intent(inout) :: H1LC(:,:)

    !> Container for either bond populations or bond energy
    real(dp), allocatable, intent(inout) :: bondWork(:, :)

    !> Last calculated bond population (for tagged output)
    real(dp), intent(inout) :: lastBondPopul

    !> Pairwise bond population output file ID
    type(TFileDescr), intent(out) :: fdBondPopul

    !> Pairwise bond energy output file ID
    type(TFileDescr), intent(out) :: fdBondEnergy

    !> Simulation time (in atomic units)
    real(dp), intent(in) :: time

    !> Environment settings
    type(TEnvironment), intent(in) :: env

    !> Error status
    type(TStatus), intent(out) :: errStatus

    real(dp), allocatable :: T2(:,:), T3(:,:)
    complex(dp), allocatable :: T4(:,:)
    integer :: iSpin, iOrb, iOrb2, iKS, iK
    type(TFileDescr) :: fillingsIn
    integer :: nLocalCols, nLocalRows, i, j

    allocate(rhoPrim(size(ints%hamiltonian, dim=1), this%nSpin))
    allocate(ErhoPrim(size(ints%hamiltonian, dim=1)))
    this%nSparse = size(H0)
    allocate(ham0(size(H0)))
    ham0(:) = H0

#:if WITH_SCALAPACK
    nLocalRows = size(eigvecsReal, dim=1)
    nLocalCols = size(eigvecsReal, dim=2)
#:else
    nLocalRows = this%denseDesc%fullSize
    nLocalCols = this%denseDesc%fullSize
#:endif

    if (this%tRealHS) then
      allocate(T2(nLocalRows,nLocalCols))
      allocate(T3(nLocalRows,nLocalCols))
    else
      allocate(T4(nLocalRows,nLocalCols))
    end if

    if (.not. this%tReadRestart) then
      Ssqr(:,:,:) = 0.0_dp
      Sinv(:,:,:) = 0.0_dp

#:if WITH_SCALAPACK
      do iKS = 1, this%parallelKS%nLocalKS
        if (this%tRealHS) then
          call  unpackHSRealBlacs(env%blacs, ints%overlap, iNeighbour,&
           & nNeighbourSK, iSparseStart, img2CentCell, this%denseDesc, T2)
          Ssqr(:,:,iKS) = cmplx(T2, 0, dp)

          call psymmatinv(this%denseDesc%blacsOrbSqr, T2, errStatus)
          Sinv(:,:,iKS) = cmplx(T2, 0, dp)
        ! TODO: add here complex overlap matrix with blacs
        end if
      end do
#:else
      do iKS = 1, this%parallelKS%nLocalKS
        if (this%tRealHS) then
          call unpackHS(T2, ints%overlap, iNeighbour, nNeighbourSK, iSquare, iSparseStart,&
              & img2CentCell)
          call adjointLowerTriangle(T2)
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
          call unpackHS(T4, ints%overlap, this%kPoint(:,iK), iNeighbour, nNeighbourSK,&
              & this%iCellVec, this%cellVec, iSquare, iSparseStart, img2CentCell)
          call adjointLowerTriangle(T4)
          Ssqr(:,:,iKS) = T4
          Sinv(:,:,iKS) = cmplx(0,0,dp)
          do iOrb = 1, this%nOrbs
            Sinv(iOrb, iOrb, iKS) = 1.0_dp
          end do
          call gesv(T4, Sinv(:,:,iKS))
        end if
      end do
#:endif

      write(stdOut,"(A)")'S inverted'

#:if WITH_SCALAPACK
      do iKS = 1, this%parallelKS%nLocalKS
        iK = this%parallelKS%localKS(1, iKS)
        iSpin = this%parallelKS%localKS(2, iKS)
        if (this%tRealHS) then
          call  unpackHSRealBlacs(env%blacs, ints%hamiltonian(:,iSpin), iNeighbour,&
               & nNeighbourSK, iSparseStart, img2CentCell, this%denseDesc, T3)
          H1(:,:,iKS) = cmplx(T3, 0, dp)
        ! TODO: add here complex Hamiltomian matrix with blacs
        end if
      end do
#:else
      do iKS = 1, this%parallelKS%nLocalKS
        iK = this%parallelKS%localKS(1, iKS)
        iSpin = this%parallelKS%localKS(2, iKS)
        if (this%tRealHS) then
          call unpackHS(T3, ints%hamiltonian(:,iSpin), iNeighbour, nNeighbourSK, iSquare,&
              & iSparseStart, img2CentCell)
          call adjointLowerTriangle(T3)
          H1(:,:,iKS) = cmplx(T3, 0, dp)
        else
          call unpackHS(H1(:,:,iKS), ints%hamiltonian(:,iSpin), this%kPoint(:,iK), iNeighbour,&
              & nNeighbourSK, this%iCellVec, this%cellVec, iSquare, iSparseStart, img2CentCell)
          call adjointLowerTriangle(H1(:,:,iKS))
        end if
      end do
#:endif

      call updateDQ(this, ints, iNeighbour, nNeighbourSK, img2CentCell, iSquare,&
          & iSparseStart, Dsqr, Qsqr)

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
      call openFile(fillingsIn, "fillings.in", mode="r")
      do iSpin=1,this%nSpin
        do iOrb=1,this%nOrbs
          read(fillingsIn%unit, *) filling(iOrb,1,iSpin)
        end do
      end do
      call closeFile(fillingsIn)
    end if

    if (.not.this%tReadRestart) then
      rho(:,:,:) = 0.0_dp
#:if WITH_SCALAPACK
      do iKS = 1, this%parallelKS%nLocalKS
        iK = this%parallelKS%localKS(1, iKS)
        iSpin = this%parallelKS%localKS(2, iKS)
        if (this%tRealHS) then
          call makeDensityMtxRealBlacs(env%blacs%orbitalGrid, this%denseDesc%blacsOrbSqr, filling(:,1,iSpin),&
          & eigvecsReal(:,:,iKS), T2)
          rho(:,:,iKS) = cmplx(T2, 0, dp)
        ! TODO: add here density matrix for complex HS with MPI
        end if
      end do
#:else
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
#:endif
    end if

    call TPotentials_init(potential, orb, this%nAtom, this%nSpin, &
        & this%nDipole, this%nQuadrupole)
    call TEnergies_init(energy, this%nAtom, this%nSpin)

    if (isDftbU .or. allocated(onSiteElements)) then
      allocate(qBlock(orb%mOrb, orb%mOrb, this%nAtom, this%nSpin))
    end if

    if (this%tNetCharges) then
      allocate(qNetAtom(this%nAtom))
    end if

    if (this%isHybridXc) then
      allocate(H1LC(this%nOrbs, this%nOrbs))
    end if

    if (this%tBondE .or. this%tBondP) then
      allocate(bondWork(this%nAtom, this%nAtom))
    end if
    if (this%tBondE) then
      call openOutputFile(this, fdBondEnergy, 'bondenergy.bin', isBinary = .true.)
    end if
    if (this%tBondP) then
      call openOutputFile(this, fdBondPopul, 'bondpop.bin', isBinary = .true.)
    end if
    call getBondPopulAndEnergy(this, bondWork, lastBondPopul, rhoPrim, ham0, ints, iNeighbour,&
        & nNeighbourSK, iSparseStart, img2CentCell, iSquare, fdBondEnergy, fdBondPopul, time)

    if (this%tRealHS) then
      deallocate(T2)
      deallocate(T3)
    else
      deallocate(T4)
    end if

  end subroutine initializeTDVariables


  !> Performs a step backwards to boot the dynamics using the Euler algorithm.
  !> Output is rho(deltaT) called rhoNew, input is rho(t=0) (ground state) called rho
  subroutine initializePropagator(this, env, step, rho, rhoNew, H1, Sinv, coordAll, skOverCont,&
      & orb, neighbourList, nNeighbourSK, img2CentCell, iSquare)

    !> ElecDynamics instance
    type(TElecDynamics), intent(inout) :: this

    !> Environment settings
    type(TEnvironment), intent(inout) :: env

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

    !> Data type for atomic orbital information
    type(TOrbitals), intent(in) :: orb

    !> ADT for neighbour parameters
    type(TNeighbourList), intent(in) :: neighbourList

    !> Nr. of neighbours for atoms out to max interaction distance (excluding Ewald terms)
    integer, intent(in) :: nNeighbourSK(:)

    !> Index array for start of atomic block in dense matrices
    integer, intent(in) :: iSquare(:)

    !> Image atoms to their equivalent in the central cell
    integer, intent(in) :: img2CentCell(:)

    integer :: iKS
    complex(dp), allocatable :: RdotSprime(:,:)

    allocate(RdotSprime(this%nOrbs,this%nOrbs))

    if (this%tIons) then
      if (allocated(this%tblite)) then
        call this%tblite%buildRdotSprime(env, RdotSprime, coordAll, this%movedVelo, &
            & this%species, nNeighbourSK, neighbourList%iNeighbour, img2CentCell, iSquare, orb)
      else
        call getRdotSprime(this, RdotSprime, coordAll, skOverCont, orb, img2CentCell, &
            &neighbourList, nNeighbourSK, iSquare)
      end if
    else
      RdotSprime(:,:) = 0.0_dp
    end if

    do iKS = 1, this%parallelKS%nLocalKS
      if (this%tIons .or. (.not. this%tRealHS)) then
        H1(:,:,iKS) = RdotSprime + imag * H1(:,:,iKS)
        call propagateRho(this, rhoNew(:,:,iKS), rho(:,:,iKS), H1(:,:,iKS), Sinv(:,:,iKS), step)
      else
        ! The following line is commented to make the fast propagate work since it needs a real H
        !H1(:,:,iKS) = imag * H1(:,:,iKS)
#:if WITH_SCALAPACK
        call propagateRhoRealHBlacs(this, rhoNew(:,:,iKS), rho(:,:,iKS), H1(:,:,iKS), Sinv(:,:,iKS), &
            & step)
#:else
        call propagateRhoRealH(this, rhoNew(:,:,iKS), rho(:,:,iKS), H1(:,:,iKS), Sinv(:,:,iKS),&
            & step)
#:endif
      end if
    end do

  end subroutine initializePropagator


  !> Propagate rho, notice that H = iH (coefficients are real)
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

#:if WITH_SCALAPACK
  !> Propagate rho with Scalapack for real Hamiltonian (used for frozen nuclei dynamics and gamma point periodic)
  subroutine propagateRhoRealHBlacs(this, rhoOld, rho, H1, Sinv, step)

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

    real(dp), allocatable :: T1R(:,:), T2R(:,:), T3R(:,:),T4R(:,:)

    integer :: nLocalCols, nLocalRows, i, j

    nLocalRows = size(rho, dim=1)
    nLocalCols = size(rho, dim=2)
    allocate(T1R(nLocalRows,nLocalCols))
    allocate(T2R(nLocalRows,nLocalCols))
    allocate(T3R(nLocalRows,nLocalCols))
    allocate(T4R(nLocalRows,nLocalCols))

    T2R(:,:) = real(H1)
    T1R(:,:) = real(rho)

    call pblasfx_psymm(T2R, this%denseDesc%blacsOrbSqr, T1R, this%denseDesc%blacsOrbSqr,&
         & T3R, this%denseDesc%blacsOrbSqr, side="R")
    ! T3R= Re[Rho]*Re[H]

    T2R(:,:) = T3R
    T1R(:,:) = real(Sinv)

    call pblasfx_psymm(T1R, this%denseDesc%blacsOrbSqr, T2R, this%denseDesc%blacsOrbSqr,&
       & T3R, this%denseDesc%blacsOrbSqr, side="R")
    ! T3R= Re[Rho]*Re[H]*Re[Sinv]

    T1R(:,:) = aimag(rho)
    T2R(:,:) = real(H1)

    call pblasfx_psymm(T2R, this%denseDesc%blacsOrbSqr, T1R, this%denseDesc%blacsOrbSqr,&
    & T4R, this%denseDesc%blacsOrbSqr, side="R")
    ! T4R= Im[Rho]*Im[H]

    T2R(:,:) = T4R
    T1R(:,:) = real(Sinv)

    call pblasfx_psymm(T1R, this%denseDesc%blacsOrbSqr, T2R, this%denseDesc%blacsOrbSqr,&
       & T4R, this%denseDesc%blacsOrbSqr, side="R")
    ! T4R= Im[Rho]*Im[H]*Im[Sinv]

    call pblasfx_ptran(T3R, this%denseDesc%blacsOrbSqr, T1R, this%denseDesc%blacsOrbSqr)
    call pblasfx_ptran(T4R, this%denseDesc%blacsOrbSqr, T2R, this%denseDesc%blacsOrbSqr)
    ! T1R = Transpose(T3R)
    ! T4R = Transpose(T2R)
    ! build the commutator combining the real and imaginary parts of the previous result

    !$OMP WORKSHARE
    rhoOld(:,:) = rhoOld + cmplx(0, step, dp) * (T3R + imag * T4R)&
        & + cmplx(0, -step, dp) * (T1R - imag * T2R)
    !$OMP END WORKSHARE

    deallocate(T1R, T2R, T3R, T4R)

  end subroutine propagateRhoRealHBlacs
#:endif

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

    real(dp), allocatable :: T1R(:,:), T2R(:,:), T3R(:,:),T4R(:,:)

    integer :: i, j

    allocate(T1R(this%nOrbs,this%nOrbs))
    allocate(T2R(this%nOrbs,this%nOrbs))
    allocate(T3R(this%nOrbs,this%nOrbs))
    allocate(T4R(this%nOrbs,this%nOrbs))

    ! The code below takes into account that Sinv and H1 are real, this is twice as fast as the
    ! original above (propageteRho)

    ! get the real part of Sinv and H1
    T1R(:,:) = real(H1)
    T2R(:,:) = real(Sinv)

    call gemm(T3R,T2R,T1R)
    T2R(:,:) = T3R

    ! calculate the first term products for the real and imaginary parts independently
    T1R(:,:) = real(rho)
    call gemm(T3R,T2R,T1R)

    T1R(:,:) = aimag(rho)
    call gemm(T4R,T2R,T1R)

    ! build the commutator combining the real and imaginary parts of the previous result
    !$OMP WORKSHARE
    rhoOld(:,:) = rhoOld + cmplx(0, -step, dp) * (T3R + imag * T4R)&
        & + cmplx(0, step, dp) * transpose(T3R - imag * T4R)
    !$OMP END WORKSHARE

  end subroutine propagateRhoRealH


  !> Initialize output files
  subroutine initTDOutput(this, dipoleDat, qDat, energyDat, populDat, forceDat, coorDat,&
      & atomEnergyDat)

    !> ElecDynamics instance
    type(TElecDynamics), intent(in) :: this

    !> Dipole output file ID
    type(TFileDescr), intent(out) :: dipoleDat

    !> Charge output file ID
    type(TFileDescr), intent(out) :: qDat

    !> Energy output file ID
    type(TFileDescr), intent(out) :: energyDat

    !> IDs for Population output files
    type(TFileDescr), intent(out) :: populDat(:)

    !> Forces output file ID
    type(TFileDescr), intent(out) :: forceDat

    !> Coords  output file ID
    type(TFileDescr), intent(out) :: coorDat

    !> Atom-resolved energy output file ID
    type(TFileDescr), intent(out) :: atomEnergyDat

    character(20) :: dipoleFileName
    character(1) :: strSpin
    character(3) :: strK
    integer :: iSpin, iKS, iK, iErr

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
    call openOutputFile(this, dipoleDat, dipoleFileName)

    write(dipoleDat%unit, "(A)", advance = "NO")"#           time (fs)    |"
    select case(this%nSpin)
    case(1)
      write(dipoleDat%unit, "(A)", advance = "NO")"     mu_x (e.angstrom)   |"
      write(dipoleDat%unit, "(A)", advance = "NO")"     mu_y (e.angstrom)   |"
      write(dipoleDat%unit, "(A)", advance = "NO")"     mu_z (e.angstrom)   |"
    case(2)
      write(dipoleDat%unit, "(A)", advance = "NO")"  mu_x (up) (e.angstrom) |"
      write(dipoleDat%unit, "(A)", advance = "NO")"  mu_y (up) (e.angstrom) |"
      write(dipoleDat%unit, "(A)", advance = "NO")"  mu_z (up) (e.angstrom) |"
      write(dipoleDat%unit, "(A)", advance = "NO")" mu_x (down) (e.angstrom)|"
      write(dipoleDat%unit, "(A)", advance = "NO")" mu_y (down) (e.angstrom)|"
      write(dipoleDat%unit, "(A)", advance = "NO")" mu_z (down) (e.angstrom)|"
    end select
    write(dipoleDat%unit, "(A)")

    if (this%tdWriteExtras) then
      call openOutputFile(this, qDat, 'qsvst.dat')
      write(qDat%unit, "(A)", advance = "NO")"#             time (fs)      |"
      write(qDat%unit, "(A)", advance = "NO")"   total net charge (e)  |"
      write(qDat%unit, "(A)", advance = "NO")"   charge (atom_1) (e)   |"
      write(qDat%unit, "(A)", advance = "NO")"   charge (atom_2) (e)   |"
      write(qDat%unit, "(A)", advance = "NO")"        ...        |"
      write(qDat%unit, "(A)", advance = "NO")"   charge (atom_N) (e)   |"
      write(qDat%unit, "(A)")

      call openOutputFile(this, energyDat, 'energyvst.dat')
      write(energyDat%unit, "(A)", advance = "NO")"#                  time (fs)         |"
      write(energyDat%unit, "(A)", advance = "NO")"        E total (H)         |"
      write(energyDat%unit, "(A)", advance = "NO")"        E non-SCC (H)       |"
      write(energyDat%unit, "(A)", advance = "NO")"            E SCC (H)       |"
      write(energyDat%unit, "(A)", advance = "NO")"           E spin (H)       |"
      write(energyDat%unit, "(A)", advance = "NO")"       E external (H)       |"
      write(energyDat%unit, "(A)", advance = "NO")"            E rep (H)       |"
      write(energyDat%unit, "(A)", advance = "NO")"E kinetic nuclear (H)       |"
      write(energyDat%unit, "(A)", advance = "NO")"     E dispersion (H)       |"
      write(energyDat%unit, "(A)")

      if (this%tForces) then
        call openOutputFile(this, forceDat, 'forcesvst.dat')
        write(forceDat%unit, "(A)", advance = "NO")"#           time (fs)       |"
        write(forceDat%unit, "(A)", advance = "NO")&
            & " force (atom_1) (H/b)   |  force (atom_2) (H/b)  |"
        write(forceDat%unit, "(A)", advance = "NO")&
            & "           ...          |  force (atom_N) (H/b)  |"
        write(forceDat%unit, "(A)")
      end if

      if (this%tIons) then
        call openOutputFile(this, coorDat, 'tdcoords.xyz')
      end if
    end if

    if (this%tPopulations) then
      do iKS = 1, this%parallelKS%nLocalKS
        iSpin = this%parallelKS%localKS(2, iKS)
        write(strSpin,'(i1)')iSpin
        if (this%tRealHS) then
          call openOutputFile(this, populDat(iKS), 'molpopul' // trim(strSpin) // '.dat')
          write(populDat(iKS)%unit, "(A,A)")&
              & "#  GS molecular orbital populations, spin channel : ", trim(strSpin)
        else
          iK = this%parallelKS%localKS(1, iKS)
          write(strK,'(i0.3)')iK
          call openOutputFile(this, populDat(iKS),&
              & 'molpopul' // trim(strSpin) // '-' // trim(strK) // '.dat')
          write(populDat(iKS)%unit, "(A,A,A,A)")&
              & "#  GS molecular orbital populations, spin channel : ", trim(strSpin), ",&
              & k-point number: ", trim(strK)
        end if
        write(populDat(iKS)%unit, "(A)", advance = "NO")"#          time (fs)            |"
        write(populDat(iKS)%unit, "(A)", advance = "NO")"   population (orb 1)       |"
        write(populDat(iKS)%unit, "(A)", advance = "NO")"    population (orb 2)      |"
        write(populDat(iKS)%unit, "(A)", advance = "NO")"           ...              |"
        write(populDat(iKS)%unit, "(A)", advance = "NO")"    population (orb N)      |"
        write(populDat(iKS)%unit, "(A)")
      end do
    end if

    iErr = -999
    if (this%tPump) then
      call execute_command_line("mkdir "//trim(pumpFilesDir), exitstat=iErr)
      if (iErr /= 0) then
        write (stdOut,*) 'cannot create '//trim(pumpFilesDir)//', error status of mkdir: ', iErr
      end if
    end if

    if (this%tWriteAtomEnergies) then
      call openOutputFile(this, atomEnergyDat, 'atomenergies.dat')
      write(atomEnergyDat%unit, "(A)", advance = "NO")"#             time (fs)      |"
      write(atomEnergyDat%unit, "(A)", advance = "NO")"   E total (H)  |"
      write(atomEnergyDat%unit, "(A)", advance = "NO")"   E (atom_1) (H)   |"
      write(atomEnergyDat%unit, "(A)", advance = "NO")"   E (atom_2) (H)   |"
      write(atomEnergyDat%unit, "(A)", advance = "NO")"        ...        |"
      write(atomEnergyDat%unit, "(A)", advance = "NO")"   E (atom_N) (H)   |"
      write(atomEnergyDat%unit, "(A)")
    end if

  end subroutine initTDOutput


  !> Close output files
  subroutine closeTDOutputs(dipoleDat, qDat, energyDat, populDat, forceDat, coorDat,&
      & fdBondPopul, fdBondEnergy, atomEnergyDat)

    !> Dipole output file ID
    type(TFileDescr), intent(inout) :: dipoleDat

    !> Charge output file ID
    type(TFileDescr), intent(inout) :: qDat

    !> Energy output file ID
    type(TFileDescr), intent(inout) :: energyDat

    !> Populations output file ID
    type(TFileDescr), intent(inout) :: populDat(:)

    !> Forces output file ID
    type(TFileDescr), intent(inout) :: forceDat

    !> Coords output file ID
    type(TFileDescr), intent(inout) :: coorDat

    !> Pairwise bond population output file ID
    type(TFileDescr), intent(inout) :: fdBondPopul

    !> Pairwise bond energy output file ID
    type(TFileDescr), intent(inout) :: fdBondEnergy

    !> Atom-resolved energy output file ID
    type(TFileDescr), intent(inout) :: atomEnergyDat

    call closeFile(dipoleDat)
    call closeFile(qDat)
    call closeFile(energyDat)
    call closeFile(forceDat)
    call closeFile(coorDat)
    call closeFile(populDat)
    call closeFile(fdBondPopul)
    call closeFile(fdBondEnergy)
    call closeFile(atomEnergyDat)

  end subroutine closeTDOutputs


  !> Open files in different ways depending on their previous existence
  subroutine openOutputFile(this, fileDescr, fileName, isBinary)

    !> ElecDynamics instance
    type(TElecDynamics), intent(in) :: this

    !> File descriptor
    type(TFileDescr), intent(out) :: fileDescr

    !> Name of the file to open
    character(*), intent(in) :: fileName

    !> Should this be a binary file?
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
      call openFile(fileDescr, newName,&
          & options=TOpenOptions(form="unformatted", access="stream", action="write"))
    else
      call openFile(fileDescr, newName, mode="w")
    end if

  end subroutine openOutputFile


  !> Write results to file
  subroutine writeTDOutputs(this, dipoleDat, qDat, energyDat, forceDat, coorDat, fdBondPopul,&
      & fdBondEnergy, atomEnergyDat, time, energy, energyKin, dipole, deltaQ, coord, totalForce,&
      & iStep)

    !> ElecDynamics instance
    type(TElecDynamics), intent(in) :: this

    !> Data type for energy components and total
    type(TEnergies), intent(in) :: energy

    !> Dipole output file ID
    type(TFileDescr), intent(in) :: dipoleDat

    !> Charge output file ID
    type(TFileDescr), intent(in) :: qDat

    !> Energy output file ID
    type(TFileDescr), intent(in) :: energyDat

    !> Elapsed simulation time
    real(dp), intent(in) :: time

    !> Dipole moment
    real(dp), intent(in) :: dipole(:,:)

    !> Negative gross charge
    real(dp), intent(in) :: deltaQ(:,:)

    !> Current step of the propagation
    integer, intent(in) :: iStep

    !> Forces output file ID
    type(TFileDescr), intent(in) :: forceDat

    !> Coords output file ID
    type(TFileDescr), intent(in) :: coorDat

    !> Pairwise bond population output file ID
    type(TFileDescr), intent(in) :: fdBondPopul

    !> Pairwise bond energy output file ID
    type(TFileDescr), intent(in) :: fdBondEnergy

    !> Atomic coordinates
    real(dp), intent(in) :: coord(:,:)

    !> Kinetic energy
    real(dp), intent(in) :: energyKin

    !> Forces (3, nAtom)
    real(dp), intent(in) :: totalForce(:,:)

    !> Atom-resolved energy output file ID
    type(TFileDescr), intent(in) :: atomEnergyDat

    real(dp) :: auxVeloc(3, this%nAtom)
    integer :: iAtom, iSpin, iDir

    write(dipoleDat%unit, '(7F25.15)') time * au__fs, ((dipole(iDir, iSpin) * Bohr__AA, iDir=1, 3),&
        & iSpin=1, this%nSpin)

    if (this%tdWriteExtras) then
      write(energydat%unit, '(9F30.15)') time * au__fs, energy%Etotal, energy%EnonSCC, energy%eSCC,&
          & energy%Espin, energy%Eext, energy%Erep, energyKin, energy%eDisp
    end if

    if (mod(iStep, this%writeFreq) == 0) then
      if (this%tdWriteExtras) then
        write(qDat%unit, "(2X,2F25.15)", advance="no") time * au__fs, -sum(deltaQ)
        do iAtom = 1, this%nAtom
          write(qDat%unit, "(F25.15)", advance="no")-sum(deltaQ(iAtom,:))
        end do
        write(qDat%unit,*)
      end if
    end if

    if (this%tIons .and. (mod(iStep,this%writeFreq) == 0)) then
      if (this%tdWriteExtras) then
        auxVeloc = 0.0_dp
        auxVeloc(:, this%indMovedAtom) = this%movedVelo
        write(coorDat%unit,'(I5)')this%nAtom
        write(coorDat%unit,*) 'MD step:', iStep, 'time', time * au__fs
        do iAtom=1,this%nAtom
          write(coorDat%unit, '(A2, 6F16.8)') trim(this%speciesName(this%species(iAtom))), &
              &coord(:, iAtom) * Bohr__AA, auxVeloc(:, iAtom) * Bohr__AA / au__fs
        end do
      endif
    end if

    if (this%tForces .and. (mod(iStep,this%writeFreq) == 0)) then
      if (this%tdWriteExtras) then
        write(forceDat%unit, "(F25.15)", advance="no") time * au__fs
        do iAtom = 1, this%nAtom
          write(forceDat%unit, "(3F25.15)", advance="no") totalForce(:,iAtom)
        end do
        write(forceDat%unit,*)
      end if
    end if

    if (this%tWriteAtomEnergies) then
        write(atomEnergyDat%unit, "(2X,2F25.15)", advance="no") time * au__fs, energy%Etotal
        do iAtom = 1, this%nAtom
          write(atomEnergyDat%unit, "(F25.15)", advance="no")this%energy%atomTotal(iAtom)
        end do
      write(atomEnergyDat%unit, *)
    end if

    ! Flush output every 5% of the simulation
    if (mod(iStep, max(this%nSteps / 20, 1)) == 0 .and. iStep > this%writeFreq) then
      if (this%tdWriteExtras) then
        flush(qDat%unit)
        flush(energyDat%unit)
        if (this%tIons) then
          flush(coorDat%unit)
        end if
        if (this%tForces) then
          flush(forceDat%unit)
        end if
        if (this%tBondP) then
          flush(fdBondPopul%unit)
        end if
        if (this%tBondE) then
          flush(fdBondEnergy%unit)
        end if
        if (this%tWriteAtomEnergies) then
          flush(atomEnergyDat%unit)
        end if
      end if
    end if

  end subroutine writeTDOutputs


  !> Initialize matrices for populations
  !! Note, this will need to get generalised for complex eigenvectors
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


  !> Updates Eiginv and EiginvAdj if nuclear dynamics is done
  !! important to call after H1 has been updated with new charges and before D is included in H1
  subroutine updateBasisMatrices(this, env, electronicSolver, Eiginv, EiginvAdj, H1, Ssqr,&
      & errStatus)

    !> ElecDynamics instance
    type(TElecDynamics), intent(in) :: this

    !> Environment
    type(TEnvironment), intent(in) :: env

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

    !> Error status
    type(TStatus), intent(out) :: errStatus

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
      call diagDenseMtx(env, electronicSolver, 'V', T1, Ssqr(:,:,iKS), eigen, errStatus)
      @:PROPAGATE_ERROR(errStatus)
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
    type(TFileDescr), intent(in) :: populDat(:)

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
    write(populDat(iKS)%unit,'(*(2x,F25.15))', advance='no') time * au__fs
    do ii = 1, size(occ)
      write(populDat(iKS)%unit,'(*(2x,F25.15))', advance='no')occ(ii)
    end do
    write(populDat(iKS)%unit,*)

  end subroutine getTDPopulations


  !> Write time-dependent tagged information to autotestTag file
  subroutine writeTDAutotest(this, dipole, energy, deltaQ, coord, totalForce, occ, lastBondPopul,&
      & taggedWriter)

    !> ElecDynamics instance
    type(TElecDynamics), intent(in) :: this

    !> Dipole moment
    real(dp), intent(in) :: dipole(:,:)

    !> Data type for energy components and total
    type(TEnergies), intent(in) :: energy

    !> Negative gross charge
    real(dp), intent(in) :: deltaQ(:,:)

    !> Atomic coordinates
    real(dp), intent(in) :: coord(:,:)

    !> Forces (3, nAtom)
    real(dp), intent(in) :: totalForce(:,:)

    !> Molecular orbital projected populations
    real(dp), intent(in) :: occ(:)

    !> Last bond population in the run
    real(dp), intent(in) :: lastBondPopul

    !> Tagged writer object
    type(TTaggedWriter), intent(inout) :: taggedWriter

    type(TFileDescr) :: fdAutotest

    call openFile(fdAutotest, trim(this%autotestTag), mode="a")
    call taggedWriter%write(fdAutotest%unit, tagLabels%tdenergy, energy%eSCC)
    call taggedWriter%write(fdAutotest%unit, tagLabels%tddipole, dipole)
    call taggedWriter%write(fdAutotest%unit, tagLabels%tdcharges, deltaQ)
    if (this%tIons) then
      call taggedWriter%write(fdAutotest%unit, tagLabels%ehrencoords, coord)
      call taggedWriter%write(fdAutotest%unit, tagLabels%ehrenvelos, this%movedVelo)
    end if
    if (this%tForces) then
      call taggedWriter%write(fdAutotest%unit, tagLabels%ehrenforces, totalForce)
    end if
    if (this%tPopulations) then
      call taggedWriter%write(fdAutotest%unit, tagLabels%tdprojocc, occ)
    end if
    if (this%tBondP) then
      call taggedWriter%write(fdAutotest%unit, tagLabels%sumBondPopul, lastBondPopul)
    end if
    if (this%tWriteAtomEnergies) then
      call taggedWriter%write(fdAutotest%unit, tagLabels%atomenergies, energy%atomTotal)
    end if
    call closeFile(fdAutotest)

  end subroutine writeTDAutotest


  !> Initialize ion dynamics
  subroutine initIonDynamics(this, coordNew, coord)

    !> ElecDynamics instance
    type(TElecDynamics), intent(inout) :: this

    !> Coordinates of next step
    real(dp), intent(out) :: coordNew(:,:)

    !> Atomic coordinates
    real(dp), intent(in) :: coord(:,:)

    ! Data for the velocity verlet integrator
    type(TVelocityVerlet), allocatable :: pVelocityVerlet

    if (this%nDynamicsInit == 0) then
      allocate(pVelocityVerlet)
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
    ! Velocities should actually be v(t+0.5*dt), not v(t),
    ! like this: this%movedVelo(:,:) = this%movedVelo + 0.5_dp * movedAccel * this%dt
    coordNew(:,:) = coord
    coordNew(:,this%indMovedAtom) = coordNew(:,this%indMovedAtom) &
        & + this%movedVelo(:,:) * this%dt

    ! This re-initializes the velocity Verlet propagator with coordNew
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


  !> Calculates non-SCC hamiltonian and overlap for new geometry and reallocates sparse arrays
  subroutine updateH0S(this, boundaryCond, Ssqr, Sinv, coord, orb, neighbourList, nNeighbourSK,&
      & iSquare, iSparseStart, img2CentCell, skHamCont, skOverCont, ham0, ints, env, rhoPrim,&
      & ErhoPrim, coordAll, errStatus, Dsqr, Qsqr)

    !> ElecDynamics instance
    type(TElecDynamics), intent(inout), target :: this

    !> Boundary conditions on the calculation
    type(TBoundaryConditions), intent(in) :: boundaryCond

    !> Square overlap inverse
    complex(dp), intent(inout) :: Sinv(:,:,:)

    !> Square overlap matrix
    complex(dp), intent(inout), allocatable :: Ssqr(:,:,:)

    !> Square dipole matrix
    complex(dp), intent(inout), optional :: Dsqr(:,:,:,:)

    !> Square quadrupole matrix
    complex(dp), intent(inout), optional :: Qsqr(:,:,:,:)

    !> Local sparse storage for non-SCC hamiltonian
    real(dp), allocatable, intent(inout) :: ham0(:)

    !> Integral container
    type(TIntegral), intent(inout) :: ints

    !> Atomic coordinates
    real(dp), allocatable, intent(inout) :: coord(:,:)

    !> Coords of the atoms (3, nAllAtom)
    real(dp), allocatable, intent(inout) :: coordAll(:,:)

    !> ADT for neighbour parameters
    type(TNeighbourList), intent(inout) :: neighbourList

    !> Nr. of neighbours for atoms out to max interaction distance (excluding Ewald terms)
    integer, intent(inout) :: nNeighbourSK(:)

    !> Index array for location of atomic blocks in large sparse arrays
    integer, allocatable, intent(inout) :: iSparseStart(:,:)

    !> Image atoms to their equivalent in the central cell
    integer, allocatable, intent(inout) :: img2CentCell(:)

    !> Index array for start of atomic block in dense matrices
    integer, intent(in) :: iSquare(:)

    !> Data type for atomic orbital information
    type(TOrbitals), intent(in) :: orb

    !> Raw H^0 hamiltonian data
    type(TSlakoCont), intent(in) :: skHamCont

    !> Raw overlap data
    type(TSlakoCont), intent(in) :: skOverCont

    !> Environment settings
    type(TEnvironment), intent(inout) :: env

    !> Sparse density matrix
    real(dp), allocatable, intent(inout) :: rhoPrim(:,:)

    !> Energy weighted density matrix
    real(dp), allocatable, intent(inout) :: ErhoPrim(:)

    !> Error status
    type(TStatus), intent(out) :: errStatus

    real(dp), allocatable :: Sreal(:,:), SinvReal(:,:)
    complex(dp), allocatable :: T4(:,:)
    real(dp) :: coord0Fold(3,this%nAtom)
    integer :: nAllAtom, iSpin, sparseSize, iOrb, iKS, iK

    coord0Fold(:,:) = coord
    call boundaryCond%foldCoordsToCell(coord0Fold, this%latVec)

    call updateNeighbourListAndSpecies(env, coordAll, this%speciesAll, img2CentCell, this%iCellVec,&
        & neighbourList, nAllAtom, coord0Fold, this%species, this%mCutoff, this%rCellVec, errStatus)
    @:PROPAGATE_ERROR(errStatus)
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
    call reallocateTDSparseArrays(this, ints, ham0, rhoPrim, ErhoPrim)

    if (this%tPeriodic) then
      call initLatticeVectors(this, boundaryCond)
    end if
    if (allocated(this%sccCalc)) then
      call this%sccCalc%updateCoords(env, coord, coordAll, this%speciesAll, neighbourList)
    end if
    if (allocated(this%tblite)) then
      call this%tblite%updateCoords(env, neighbourList, img2CentCell, coordAll,&
          & this%speciesAll)
    end if

    if (allocated(this%dispersion)) then
      call this%dispersion%updateCoords(env, neighbourList, img2CentCell, coordAll,&
          & this%speciesAll, errStatus)
      @:PROPAGATE_ERROR(errStatus)
    end if

    select case(this%hamiltonianType)
    case default
      @:ASSERT(.false.)
    case(hamiltonianTypes%dftb)
      call buildH0(env, ham0, skHamCont, this%atomEigVal, coordAll, nNeighbourSK, &
          & neighbourList%iNeighbour, this%speciesAll, iSparseStart, orb)
      call buildS(env, ints%overlap, skOverCont, coordAll, nNeighbourSK, neighbourList%iNeighbour,&
          & this%speciesAll, iSparseStart, orb)
    case(hamiltonianTypes%xtb)
      @:ASSERT(allocated(this%tblite))
      call this%tblite%buildSH0(env, this%speciesAll, coordAll, nNeighbourSk, &
          & neighbourList%iNeighbour, img2CentCell, iSparseStart, orb, ham0,&
          & ints%overlap, ints%dipoleBra, ints%dipoleKet, &
          & ints%quadrupoleBra, ints%quadrupoleKet)
    end select

    ! TODO: modify this routine to enable restart and ion dynamics with MPI
    if (this%tRealHS) then
      allocate(Sreal(this%nOrbs,this%nOrbs))
      allocate(SinvReal(this%nOrbs,this%nOrbs))
      Sreal = 0.0_dp

      call unpackHS(Sreal, ints%overlap, neighbourList%iNeighbour, nNeighbourSK, iSquare,&
          & iSparseStart, img2CentCell)
      call adjointLowerTriangle(Sreal)
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
        call unpackHS(T4, ints%overlap, this%kPoint(:,iK), neighbourList%iNeighbour, nNeighbourSK,&
            & this%iCellVec, this%cellVec, iSquare, iSparseStart, img2CentCell)
        call adjointLowerTriangle(T4)
        Ssqr(:,:,iKS) = T4
        Sinv(:,:,iKS) = cmplx(0,0,dp)
        do iOrb = 1, this%nOrbs
          Sinv(iOrb, iOrb, iKS) = cmplx(1,0,dp)
        end do
        call gesv(T4, Sinv(:,:,iKS))
      end do
      deallocate(T4)

    end if

    call updateDQ(this, ints, neighbourList%iNeighbour, nNeighbourSK, img2CentCell, iSquare,&
        & iSparseStart, Dsqr, Qsqr)

  end subroutine updateH0S


  !> Updates the dipole and quadrupole matrix elements
  subroutine updateDQ(this, ints, iNeighbour, nNeighbourSK, img2CentCell, iSquare,&
      & iSparseStart, Dsqr, Qsqr)

    !> ElecDynamics instance
    type(TElecDynamics), intent(inout) :: this

    !> Integral container
    type(TIntegral), intent(inout) :: ints

    !> Atomic neighbour data
    integer, intent(in) :: iNeighbour(0:,:)

    !> Number of neighbours for each of the atoms
    integer, intent(in) :: nNeighbourSK(:)

    !> Index array for start of atomic block in dense matrices
    integer, intent(in) :: iSquare(:)

    !> Index array for location of atomic blocks in large sparse arrays
    integer, intent(in) :: iSparseStart(0:,:)

    !> Image atoms to their equivalent in the central cell
    integer, intent(in) :: img2CentCell(:)

    !> Square dipole matrix
    complex(dp), intent(inout), optional :: Dsqr(:,:,:,:)

    !> Square quadrupole matrix
    complex(dp), intent(inout), optional :: Qsqr(:,:,:,:)

    real(dp), allocatable :: M3(:, :, :)
    integer :: iSpin, iKS, iK

    if (present(Dsqr)) then
      Dsqr(:,:,:,:) = 0.0_dp
      if (this%tRealHS) then
        allocate(M3(this%nDipole, this%nOrbs, this%nOrbs))
        do iKS = 1, this%parallelKS%nLocalKS
          iK = this%parallelKS%localKS(1, iKS)
          iSpin = this%parallelKS%localKS(2, iKS)
          call unpackDQ(M3, ints%dipoleBra, ints%dipoleKet, iNeighbour,&
              & nNeighbourSK, iSquare, iSparseStart, img2CentCell)
          Dsqr(:,:,:,iKS) = cmplx(M3, 0, dp)
        end do
      else
        do iKS = 1, this%parallelKS%nLocalKS
          iK = this%parallelKS%localKS(1, iKS)
          iSpin = this%parallelKS%localKS(2, iKS)
          call unpackDQ(Dsqr(:,:,:,iKS), ints%dipoleBra, ints%dipoleKet,&
              & this%kPoint(:,iK), iNeighbour, nNeighbourSK, this%iCellVec, this%cellVec,&
              & iSquare, iSparseStart, img2CentCell)
        end do
      end if
    end if

    if (present(Qsqr)) then
      Qsqr(:,:,:,:) = 0.0_dp
      if (this%tRealHS) then
        if (allocated(M3)) deallocate(M3)
        allocate(M3(this%nQuadrupole, this%nOrbs, this%nOrbs))
        do iKS = 1, this%parallelKS%nLocalKS
          iK = this%parallelKS%localKS(1, iKS)
          iSpin = this%parallelKS%localKS(2, iKS)
          call unpackDQ(M3, ints%quadrupoleBra, ints%quadrupoleKet, iNeighbour,&
              & nNeighbourSK, iSquare, iSparseStart, img2CentCell)
          Qsqr(:,:,:,iKS) = cmplx(M3, 0, dp)
        end do
      else
        do iKS = 1, this%parallelKS%nLocalKS
          iK = this%parallelKS%localKS(1, iKS)
          iSpin = this%parallelKS%localKS(2, iKS)
          call unpackDQ(Qsqr(:,:,:,iKS), ints%quadrupoleBra, ints%quadrupoleKet,&
              & this%kPoint(:,iK), iNeighbour, nNeighbourSK, this%iCellVec, this%cellVec,&
              & iSquare, iSparseStart, img2CentCell)
        end do
      end if
    end if

  end subroutine updateDQ


  !> Calculates force
  subroutine getForces(this, movedAccel, totalForce, rho, H1, Sinv, neighbourList, nNeighbourSK,&
      & img2CentCell, iSparseStart, iSquare, potential, orb, skHamCont, skOverCont, qq, q0,&
      & repulsive, coordAll, rhoPrim, ErhoPrim, iStep, env, hybridXc, deltaRho, sSqr, errStatus)

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

    !> Nr. of neighbours for atoms out to max interaction distance (excluding Ewald terms)
    integer, intent(in) :: nNeighbourSK(:)

    !> Index array for location of atomic blocks in large sparse arrays
    integer, intent(in) :: iSparseStart(:,:)

    !> Image atoms to their equivalent in the central cell
    integer, intent(in) :: img2CentCell(:)

    !> Index array for start of atomic block in dense matrices
    integer, intent(in) :: iSquare(:)

    !> Data type for atomic orbital information
    type(TOrbitals), intent(in) :: orb

    !> Raw H^0 hamiltonian data
    type(TSlakoCont), intent(in) :: skHamCont

    !> Raw overlap data
    type(TSlakoCont), intent(in) :: skOverCont

    !> Sparse density matrix
    real(dp), intent(inout) :: rhoPrim(:,:)

    !> Energy weighted density matrix
    real(dp), intent(inout) :: ErhoPrim(:)

    !> Forces (3, nAtom)
    real(dp), intent(out) :: totalForce(:,:)

    !> Acceleration on moved atoms (3, nMovedAtom)
    real(dp), intent(out) :: movedAccel(:,:)

    !> Potential acting on the system
    type(TPotentials), intent(in) :: potential

    !> Atomic occupations
    real(dp), intent(inout) :: qq(:,:,:)

    !> Reference atomic occupations
    real(dp), intent(inout) :: q0(:,:,:)

    !> Repulsive information
    class(TRepulsive), allocatable, intent(in) :: repulsive

    !> Coords of the atoms (3, nAllAtom)
    real(dp), intent(in) :: coordAll(:,:)

    !> Current step of the propagation
    integer, intent(in) :: iStep

    !> Environment settings
    type(TEnvironment), intent(inout) :: env

    !> Range separation contributions
    class(THybridXcFunc), allocatable, intent(inout) :: hybridXc

    !> Real part of density matrix, adjusted by reference charges
    complex(dp), allocatable, intent(inout) :: deltaRho(:,:,:)

    !> Square overlap matrix
    complex(dp), intent(in) :: Ssqr(:,:,:)

    !> Error status
    type(TStatus), intent(out) :: errStatus

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


    if (allocated(this%tblite)) then
      call this%tblite%updateCharges(env, this%speciesAll, neighbourList, qq, q0,&
          & this%multipole%dipoleAtom, this%multipole%quadrupoleAtom, img2CentCell, orb)
      call this%tblite%buildDerivativeShift(env, rhoPrim, ERhoPrim, coordAll, this%speciesAll,&
          & nNeighbourSK, neighbourList%iNeighbour, img2CentCell, iSparseStart, orb,&
          & potential%intBlock, potential%dipoleAtom, potential%quadrupoleAtom)
      call this%tblite%addGradients(env, neighbourList, this%speciesAll, coordAll,&
          & img2CentCell, derivs)
    else
      call derivative_shift(env, derivs, this%derivator, rhoPrim, ErhoPrim, skHamCont,&
          & skOverCont, coordAll, this%speciesAll, neighbourList%iNeighbour, nNeighbourSK, &
          & img2CentCell, iSparseStart, orb, potential%intBlock)
    end if
    if (allocated(this%sccCalc)) then
      call this%sccCalc%updateCharges(env, qq, orb, this%speciesAll, q0)
      call this%sccCalc%addForceDc(env, derivs, this%speciesAll, neighbourList%iNeighbour, &
          & img2CentCell)
    end if
    if (allocated(repulsive)) then
      call repulsive%getGradients(coordAll, this%speciesAll, img2CentCell, neighbourList,&
          & repulsiveDerivs)
    else
      repulsiveDerivs(:,:) = 0.0_dp
    end if

    if (this%isHybridXc) then
    #:if WITH_SCALAPACK
      @:RAISE_ERROR(errStatus, -1, "MPI-parallel hybrid-DFTB matrix-based force evaluation not&
          & implemented for rTD-DFTB.")
    #:else
      call hybridXc%addCamGradients_real(real(deltaRho), real(sSqr(:,:,1)), skOverCont, orb,&
          & iSquare, neighbourList%iNeighbour, nNeighbourSK, this%derivator, img2CentCell,&
          & this%speciesAll, coordAll, .false., derivs, errStatus)
      @:PROPAGATE_ERROR(errStatus)
    #:endif
    end if

    if (this%tLaser) then
      call setPresentField(this, iStep, errStatus)
      @:PROPAGATE_ERROR(errStatus)
      do iDir = 1, 3
        derivs(iDir,:) = derivs(iDir,:)&
            & - sum(q0(:,:,1) - qq(:,:,1), dim=1) * this%presentField(iDir)
      end do
    end if

    totalDeriv(:,:) = repulsiveDerivs + derivs
    if (allocated(this%dispersion)) then
      call this%dispersion%addGradients(env, neighbourList, this%speciesAll, coordAll,  &
          & img2CentCell, totalDeriv)
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

    !> Non-adiabatic coupling matrix elements
    complex(dp), intent(out) :: RdotSprime(:,:)

    !> Data type for atomic orbital information
    type(TOrbitals), intent(in) :: orb

    !> Coords of the atoms (3, nAllAtom)
    real(dp), intent(in) :: coordAll(:,:)

    !> ADT for neighbour parameters
    type(TNeighbourList), intent(in) :: neighbourList

    !> Nr. of neighbours for atoms out to max interaction distance (excluding Ewald terms)
    integer, intent(in) :: nNeighbourSK(:)

    !> Index array for start of atomic block in dense matrices
    integer, intent(in) :: iSquare(:)

    !> Image atoms to their equivalent in the central cell
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
  subroutine reallocateTDSparseArrays(this, ints, ham0, rhoPrim, ErhoPrim)

    !> ElecDynamics instance
    type(TElecDynamics), intent(in), target :: this

    !> Integral container
    type(TIntegral), intent(inout) :: ints

    !> Local sparse storage for non-SCC hamiltonian
    real(dp), allocatable, intent(inout) :: ham0(:)

    !> Sparse density matrix
    real(dp), allocatable, intent(inout) :: rhoPrim(:,:)

    !> Energy weighted density matrix
    real(dp), allocatable, intent(inout) :: ErhoPrim(:)

    deallocate(ints%hamiltonian)
    deallocate(ints%overlap)
    deallocate(ham0)
    deallocate(rhoPrim)
    allocate(ints%hamiltonian(this%nSparse, this%nSpin))
    allocate(ints%overlap(this%nSparse))
    allocate(ham0(this%nSparse))
    allocate(rhoPrim(this%nSparse, this%nSpin))

    if (allocated(ErhoPrim)) then
      deallocate(ErhoPrim)
      allocate(ErhoPrim(this%nSparse))
    end if

  end subroutine reallocateTDSparseArrays


  !> Updates SCC module with lattice vectors
  subroutine initLatticeVectors(this, boundarycond)

    !> ElecDynamics instance
    type(TElecDynamics), intent(inout), target :: this

    !> Boundary conditions on the calculation
    type(TBoundaryConditions), intent(in) :: boundaryCond

    real(dp) :: cellVol, recVecs(3,3), recVecs2p(3,3)

    cellVol = abs(determinant33(this%latVec))
    recVecs2p(:,:) = this%latVec
    call matinv(recVecs2p)
    recVecs2p = transpose(recVecs2p)
    recVecs = 2.0_dp * pi * recVecs2p
    if (allocated(this%sccCalc)) then
      call this%sccCalc%updateLatVecs(this%latVec, recVecs, boundaryCond, cellVol)
      this%mCutOff = max(this%mCutOff, this%sccCalc%getCutOff())
    end if
    if (allocated(this%tblite)) then
      call this%tblite%updateLatVecs(this%latVec)
      this%mCutOff = max(this%mCutOff, this%tblite%getRCutOff())
    end if

    if (allocated(this%dispersion)) then
      call this%dispersion%updateLatVecs(this%latVec)
      this%mCutOff = max(this%mCutOff, this%dispersion%getRCutOff())
    end if
  end subroutine initLatticeVectors


  !> Calculates repulsive and dispersion energies
  subroutine getPositionDependentEnergy(this, env, energy, coordAll, img2CentCell, neighbourList,&
      & repulsive, iAtInCentralRegion, hybridXc)

    !> ElecDynamics instance
    type(TElecDynamics), intent(inout), target :: this

    !> Environment settings
    type(TEnvironment), intent(inout) :: env

    !> Data type for energy components and total
    type(TEnergies), intent(inout) :: energy

    !> All atomic coordinates
    real(dp), intent(in) :: coordAll(:,:)

    !> Image atom indices to central cell atoms
    integer, intent(in) :: img2CentCell(:)

    !> List of neighbours for each atom
    type(TNeighbourList), intent(in) :: neighbourList

    !> Repulsive interaction data
    class(TRepulsive), allocatable, intent(inout) :: repulsive

    !> Atoms in the central cell
    integer, intent(in) :: iAtInCentralRegion(:)

    !> Range separation contributions
    class(THybridXcFunc), allocatable, intent(inout) :: hybridXc

    if (allocated(repulsive)) then
      call repulsive%updateCoords(coordAll, this%speciesAll, img2CentCell, neighbourList)
      call repulsive%getEnergy(coordAll, this%speciesAll, img2CentCell, neighbourList,&
          & energy%atomRep, energy%Erep, iAtInCentralRegion=iAtInCentralRegion)
    else
      energy%atomRep(:) = 0.0_dp
      energy%Erep = 0.0_dp
    end if
    if (allocated(this%dispersion)) then
      call calcDispersionEnergy(this%dispersion, energy%atomDisp, energy%eDisp, iAtInCentralRegion)
    else
      energy%atomDisp(:) = 0.0_dp
      energy%eDisp = 0.0_dp
    end if
    if (allocated(hybridXc)) then
      call hybridXc%addCamEnergy_real(env, energy%Efock)
    else
      energy%Efock = 0.0_dp
    end if

  end subroutine getPositionDependentEnergy


  !> Calculates bond populations and bond energies if requested
  subroutine getBondPopulAndEnergy(this, bondWork, lastBondPopul, rhoPrim, ham0, ints, iNeighbour,&
      & nNeighbourSK, iSparseStart, img2CentCell, iSquare,  fdBondEnergy, fdBondPopul, time)

    !> ElecDynamics instance
    type(TElecDynamics), intent(inout) :: this

    !> Container for either bond populations or bond energy
    real(dp), allocatable, intent(inout) :: bondWork(:, :)

    !> Last calculated bond population (for tagged output)
    real(dp), intent(inout) :: lastBondPopul

    !> Sparse density matrix
    real(dp), intent(in) :: rhoPrim(:,:)

    !> Integral container
    type(TIntegral), intent(inout) :: ints

    !> Local sparse storage for non-SCC hamiltonian
    real(dp), intent(in) :: ham0(:)

    !> Atomic neighbour data
    integer, intent(in) :: iNeighbour(0:,:)

    !> Nr. of neighbours for atoms out to max interaction distance (excluding Ewald terms)
    integer, intent(in) :: nNeighbourSK(:)

    !> Index array for location of atomic blocks in large sparse arrays
    integer, intent(in) :: iSparseStart(:,:)

    !> Image atoms to their equivalent in the central cell
    integer, intent(in) :: img2CentCell(:)

    !> Index array for start of atomic block in dense matrices
    integer, intent(in) :: iSquare(:)

    !> File descriptor for bond energy
    type(TFileDescr), intent(in) :: fdBondEnergy

    !> File descriptor for bond populations
    type(TFileDescr), intent(in) :: fdBondPopul

    !> Elapsed simulation time
    real(dp), intent(in) :: time

    integer :: iSpin

    if (this%tBondE) then
      bondWork(:,:) = 0.0_dp
      do iSpin = 1, this%nSpin
        call addPairWiseBondInfo(bondWork, rhoPrim(:,iSpin), ham0, iSquare,&
            & iNeighbour, nNeighbourSK, img2CentCell, iSparseStart)
      end do
      write(fdBondEnergy%unit) time * au__fs, sum(bondWork), bondWork
    end if
    if (this%tBondP) then
      bondWork(:,:) = 0.0_dp
      do iSpin = 1, this%nSpin
        call addPairWiseBondInfo(bondWork, rhoPrim(:,1), ints%overlap, iSquare,&
            & iNeighbour, nNeighbourSK, img2CentCell, iSparseStart)
      end do
      write(fdBondPopul%unit) time * au__fs, sum(bondWork), bondWork
      if (this%tWriteAutotest) then
        lastBondPopul = sum(bondWork)
      end if
    end if

  end subroutine getBondPopulAndEnergy


  !> Sets electric field at present timestep
  subroutine setPresentField(this, iStep, errStatus)

    !> ElecDynamics instance
    type(TElecDynamics), intent(inout) :: this

    !> Current step of the propagation
    integer, intent(in) :: iStep

    !> Error status
    type(TStatus), intent(out) :: errStatus

    if (.not. this%tdFieldThroughAPI) then
      this%presentField(:) = this%tdFunction(:, iStep)
    elseif (.not. this%tdFieldIsSet) then
      if (iStep == 0) then
        this%presentField(:) = 0.0_dp
      else
        @:RAISE_ERROR(errStatus, -1, "External field has not been set.")
      end if
    end if

  end subroutine setPresentField


  !> Handles the initializations of the variables needed for the time propagation
  subroutine initializeDynamics(this, boundaryCond, coord, orb, neighbourList, nNeighbourSK,&
      & iSquare, iSparseStart, img2CentCell, skHamCont, skOverCont, ints, env, coordAll, H0, spinW,&
      & tDualSpinOrbit, xi, thirdOrd, dftbU, onSiteElements, refExtPot, solvation, eFieldScaling,&
      & hybridXc, referenceN0, q0, repulsive, iAtInCentralRegion, eigvecsReal, eigvecsCplx,&
      & filling, qDepExtPot, tFixEf, Ef, latVec, invLatVec, iCellVec, rCellVec, cellVec,&
      & speciesAll, errStatus)

    !> ElecDynamics instance
    type(TElecDynamics), intent(inout), target :: this

    !> Boundary conditions on the calculation
    type(TBoundaryConditions), intent(in) :: boundaryCond

    !> Real Eigenvectors
    real(dp), intent(inout), allocatable :: eigvecsReal(:,:,:)

    !> Complex Eigevenctors
    complex(dp), intent(inout), allocatable :: eigvecsCplx(:,:,:)

    !> Sparse storage for non-SCC hamiltonian
    real(dp), intent(in) :: H0(:)

    !> Reference atomic occupations
    real(dp), intent(inout) :: q0(:,:,:)

    !> Reference charges from the Slater-Koster file
    real(dp), intent(in) :: referenceN0(:,:)

    !> Integral container
    type(TIntegral), intent(inout) :: ints

    !> Atomic coordinates
    real(dp), allocatable, intent(inout) :: coord(:,:)

    !> All atomic coordinates
    real(dp), allocatable, intent(inout) :: coordAll(:,:)

    !> Spin constants
    real(dp), allocatable, intent(in) :: spinW(:,:,:)

    !> Occupations
    real(dp), intent(inout) :: filling(:,:,:)

    !> Number of neighbours for each of the atoms
    integer, intent(inout) :: nNeighbourSK(:)

    !> Index array for location of atomic blocks in large sparse arrays
    integer, allocatable, intent(inout) :: iSparseStart(:,:)

    !> Image atoms to their equivalent in the central cell
    integer, allocatable, intent(inout) :: img2CentCell(:)

    !> Index array for start of atomic block in dense matrices
    integer, intent(in) :: iSquare(:)

    !> List of neighbours for each atom
    type(TNeighbourList), intent(inout) :: neighbourList

    !> Repulsive information
    class(TRepulsive), allocatable, intent(inout) :: repulsive

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

    !> DFTB+U functional (if used)
    type(TDftbU), intent(in), allocatable :: dftbU

    !> Spin orbit constants if required
    real(dp), allocatable, intent(in) :: xi(:,:)

    !> 3rd order settings
    type(TThirdOrder), intent(inout), allocatable :: thirdOrd

    !> Solvation model
    class(TSolvation), allocatable, intent(inout) :: solvation

    !> Any dielectric environment scaling
    class(TScaleExtEField), intent(in) :: eFieldScaling

    !> Range separation contributions
    class(THybridXcFunc), allocatable, intent(inout) :: hybridXc

    !> Proxy for querying Q-dependant external potentials
    type(TQDepExtPotProxy), intent(inout), allocatable :: qDepExtPot

    !> Atoms over which to sum the total energies
    integer, intent(in) :: iAtInCentralRegion(:)

    !> Whether fixed Fermi level(s) should be used. (No charge conservation!)
    logical, intent(in) :: tFixEf

    !> If tFixEf is .true. contains reservoir chemical potential, otherwise the Fermi levels found
    !> from the given number of electrons
    real(dp), intent(inout) :: Ef(:)

    !> Corrections terms for on-site elements
    real(dp), intent(in), allocatable :: onSiteElements(:,:,:,:)

    !> Reference external potential (usual provided via API)
    type(TRefExtPot) :: refExtPot

    !> Lattice vectors if periodic
    real(dp), intent(in) :: latVec(:,:)

    !> Inverse of the lattice vectors
    real(dp), intent(in) :: invLatVec(:,:)

    !> Cell vectors in absolute units
    real(dp), intent(in) :: rCellVec(:,:)

    !> Vectors (in units of the lattice constants) to cells of the lattice
    real(dp), intent(in) :: cellVec(:,:)

    !> Index of cell in cellVec and rCellVec for each atom
    integer, allocatable, intent(in) :: iCellVec(:)

    !> Species of all atoms in the system
    integer, intent(in) :: speciesAll(:)

    !> Error status
    type(TStatus), intent(inout) :: errStatus

    integer :: nLocalCols, nLocalRows
    real(dp), allocatable :: velInternal(:,:)

    this%startTime = 0.0_dp

    this%speciesAll = speciesAll
    this%nSpin = size(ints%hamiltonian(:,:), dim=2)
    if (this%nSpin > 1 .and. this%iCall == 1) then
      call qm2ud(q0)
    end if

    if (this%tRealHS) then
      this%nOrbs = size(eigvecsReal, dim=1)
    else
      this%nOrbs = size(eigvecsCplx, dim=1)
    end if

    this%nAtom = size(coord, dim=2)
    this%latVec = latVec
    this%invLatVec = invLatVec
    this%iCellVec = iCellVec
    this%rCellVec = rCellVec
    this%cellVec = cellVec

    if (allocated(this%tblite)) then
      call this%tblite%getMultipoleInfo(this%nDipole, this%nQuadrupole)
    end if
    call TMultipole_init(this%multipole, this%nAtom, this%nDipole, this%nQuadrupole, &
        & this%nSpin)

#:if WITH_SCALAPACK
    nLocalRows = size(eigvecsReal, dim=1)
    nLocalCols = size(eigvecsReal, dim=2)
#:else
    nLocalRows = this%denseDesc%fullSize
    nLocalCols = this%denseDesc%fullSize
#:endif

    allocate(this%trho(nLocalRows,nLocalCols,this%parallelKS%nLocalKS))
    allocate(this%trhoOld(nLocalRows,nLocalCols,this%parallelKS%nLocalKS))
    allocate(this%Ssqr(nLocalRows,nLocalCols,this%parallelKS%nLocalKS))
    allocate(this%Sinv(nLocalRows,nLocalCols,this%parallelKS%nLocalKS))
    allocate(this%H1(nLocalRows,nLocalCols,this%parallelKS%nLocalKS))

    if (this%nDipole > 0) then
      allocate(this%Dsqr(this%nDipole,this%nOrbs,this%nOrbs,this%parallelKS%nLocalKS))
    end if
    if (this%nQuadrupole > 0) then
      allocate(this%Qsqr(this%nQuadrupole,this%nOrbs,this%nOrbs,this%parallelKS%nLocalKS))
    end if
    allocate(this%qq(orb%mOrb, this%nAtom, this%nSpin))
    allocate(this%deltaQ(this%nAtom,this%nSpin))
    allocate(this%dipole(3,this%nSpin))
    allocate(this%chargePerShell(orb%mShell,this%nAtom,this%nSpin))

    allocate(this%occ(this%nOrbs))
    allocate(this%RdotSprime(this%nOrbs,this%nOrbs))
    allocate(this%totalForce(3, this%nAtom))
    this%RdotSprime(:,:) = 0.0_dp
    this%totalForce(:,:) = 0.0_dp
    this%occ(:) = 0.0_dp

    if (this%tReadRestart) then
      call readRestartFile(this%trho, this%trhoOld, coord, this%movedVelo, this%startTime, this%dt,&
          & restartFileName, this%tRestartAscii, errStatus)
      @:PROPAGATE_ERROR(errStatus)
      call updateH0S(this, boundaryCond, this%Ssqr, this%Sinv, coord, orb, neighbourList,&
          & nNeighbourSK, iSquare, iSparseStart, img2CentCell, skHamCont, skOverCont, this%ham0,&
          & ints, env, this%rhoPrim, this%ErhoPrim, coordAll, errStatus, this%Dsqr, this%Qsqr)
      @:PROPAGATE_ERROR(errStatus)
      if (this%tIons) then

        this%initialVelocities(:,:) = this%movedVelo

        this%ReadMDVelocities = .true.
      end if
    else if (this%iCall > 1 .and. this%tIons) then
      coord(:,:) = this%initCoord
      call updateH0S(this, boundaryCond, this%Ssqr, this%Sinv, coord, orb, neighbourList,&
          & nNeighbourSK, iSquare, iSparseStart, img2CentCell, skHamCont, skOverCont, this%ham0,&
          & ints, env, this%rhoPrim, this%ErhoPrim, coordAll, errStatus, this%Dsqr, this%Qsqr)
      @:PROPAGATE_ERROR(errStatus)
      this%initialVelocities(:,:) = this%movedVelo
      this%ReadMDVelocities = .true.
    end if
    if (this%tLaser .and. .not. this%tdFieldThroughAPI .and. this%iCall == 1) then
      call getTDFunction(this, this%startTime)
    end if

    call initializeTDVariables(this, this%trho, this%H1, this%Ssqr, this%Sinv, H0, this%ham0, &
        & this%Dsqr, this%Qsqr, ints, eigvecsReal, filling, orb, this%rhoPrim, this%potential, &
        & neighbourList%iNeighbour, nNeighbourSK, iSquare, iSparseStart, img2CentCell, this%Eiginv,&
        & this%EiginvAdj, this%energy, this%ErhoPrim, this%qBlock, this%qNetAtom, allocated(dftbU),&
        & onSiteElements, eigvecsCplx, this%H1LC, this%bondWork, this%fdBondEnergy,&
        & this%fdBondPopul, this%lastBondPopul, this%time, env, errStatus)

    if (this%tPeriodic) then
      call initLatticeVectors(this, boundaryCond)
    end if

    call initTDOutput(this, this%dipoleDat, this%qDat, this%energyDat,&
        & this%populDat, this%forceDat, this%coorDat, this%atomEnergyDat)

    ! Write density at t=0
    if (this%tPump .and. .not. this%tReadRestart) then
      allocate(velInternal(3,size(this%movedVelo, dim=2)))
        velInternal(:,:) = 0.0_dp
      call writeRestartFile(this%trho, this%trho, coord, velInternal, this%startTime, this%dt,&
          & trim(pumpFilesDir) // '/0ppdump', this%tWriteRestartAscii, errStatus)
      @:PROPAGATE_ERROR(errStatus)
      deallocate(velInternal)
    end if

    if (allocated(this%sccCalc)) then
      call this%sccCalc%updateCoords(env, coord, coordAll, this%speciesAll, neighbourList)
    end if
    if (allocated(this%tblite)) then
      call this%tblite%updateCoords(env, neighbourList, img2CentCell, coordAll,&
          & this%speciesAll)
    end if

    if (allocated(this%dispersion)) then
      call this%dispersion%updateCoords(env, neighbourList, img2CentCell, coordAll,&
          & this%speciesAll, errStatus)
      @:PROPAGATE_ERROR(errStatus)
      this%mCutOff = max(this%mCutOff, this%dispersion%getRCutOff())
    end if

    call getChargeDipole(this, this%deltaQ, this%qq, this%multipole, this%dipole, q0,&
        & this%trho, this%Ssqr, this%Dsqr, this%Qsqr, coord, iSquare, eFieldScaling, this%qBlock,&
        & this%qNetAtom, errStatus, neighbourList%iNeighbour, nNeighbourSK, orb, iSparseStart, &
        & img2CentCell,  env, ints)
    @:PROPAGATE_ERROR(errStatus)
    if (allocated(this%dispersion)) then
      call this%dispersion%updateOnsiteCharges(this%qNetAtom, orb, referenceN0,&
          & this%speciesAll(:this%nAtom), .true.)
    end if

    call updateH(this, this%H1, ints, this%ham0, this%speciesAll, this%qq, q0, coord, orb,&
        & this%potential, neighbourList, nNeighbourSK, iSquare, iSparseStart, img2CentCell, 0,&
        & this%chargePerShell, spinW, env, tDualSpinOrbit, xi, thirdOrd, this%qBlock, dftbU,&
        & onSiteElements, refExtPot, this%deltaRho, this%H1LC, this%Ssqr, solvation, hybridXc,&
        & this%dispersion, this%trho, errStatus)
    @:PROPAGATE_ERROR(errStatus)

    if (this%tForces) then
      this%totalForce(:,:) = 0.0_dp
      call getForces(this, this%movedAccel, this%totalForce, this%trho, this%H1, this%Sinv,&
          & neighbourList, nNeighbourSK, img2CentCell, iSparseStart, iSquare, this%potential, orb,&
          & skHamCont, skOverCont, this%qq, q0, repulsive, coordAll, this%rhoPrim, this%ErhoPrim,&
          & 0, env, hybridXc, this%deltaRho, this%Ssqr, errStatus)
      @:PROPAGATE_ERROR(errStatus)
    end if

    ! the ion dynamics init must be done here, as it needs the DM and outputs the velocities
    ! needed to initialise the electronic dynamics
    ! coordNew stores the coordinates at t=dt
    if (this%tIons) then
      call initIonDynamics(this, this%coordNew, coord)
    end if

    ! after calculating the TD function, set initial time to zero for probe simulations
    ! this is to properly calculate the dipole fourier transform after the simulation
    if (this%tProbe) then
      this%startTime = 0.0_dp
    end if

    ! Apply kick to rho if necessary (in restart case, check it starttime is 0 or not)
    ! TODO: Kick for MPI
    if (this%tKick .and. this%startTime < this%dt / 10.0_dp) then
      call kickDM(this, this%trho, this%Ssqr, this%Sinv, iSquare, coord)
    end if

    call getPositionDependentEnergy(this, env, this%energy, coordAll, img2CentCell, neighbourList,&
        & repulsive, iAtInCentralRegion, hybridXc)

    call getTDEnergy(this, env, this%energy, this%rhoPrim, this%trho, neighbourList, nNeighbourSK,&
        & orb, iSquare, iSparseStart, img2CentCell, this%ham0, this%qq, q0, this%potential,&
        & this%chargePerShell, this%energyKin, tDualSpinOrbit, thirdOrd, solvation, hybridXc,&
        & qDepExtPot, this%qBlock, dftbu, xi, iAtInCentralRegion, tFixEf, Ef, onSiteElements,&
        & errStatus)
    @:PROPAGATE_ERROR(errStatus)

    if (.not. this%tReadRestart .or. this%tProbe) then
      ! output ground state data
      call writeTDOutputs(this, this%dipoleDat, this%qDat, this%energyDat, &
          & this%forceDat, this%coorDat, this%fdBondPopul, this%fdBondEnergy, this%atomEnergyDat,&
          & 0.0_dp, this%energy, this%energyKin, this%dipole, this%deltaQ, coord, this%totalForce,&
          & 0)
    end if

    ! now first step of dynamics is computed (init of leapfrog and first step of nuclei)

    ! had to add the "or tKick" option to override rhoOld if tReadRestart = yes, otherwise it will
    ! be badly initialised
    if (.not.this%tReadRestart .or. (this%tKick .and. this%startTime < this%dt / 10.0_dp)) then
      ! Initialize electron dynamics
      ! rhoOld is now the GS DM, rho will be the DM at time=dt
      this%trhoOld(:,:,:) = this%trho
      call initializePropagator(this, env, this%dt, this%trhoOld, this%trho, this%H1, this%Sinv,&
          & coordAll, skOverCont, orb, neighbourList, nNeighbourSK, img2CentCell, iSquare)
    end if

    this%rho => this%trho
    this%rhoOld => this%trhoOld

    if (this%tIons) then
      coord(:,:) = this%coordNew
      call updateH0S(this, boundaryCond, this%Ssqr, this%Sinv, coord, orb, neighbourList,&
          & nNeighbourSK, iSquare, iSparseStart, img2CentCell, skHamCont, skOverCont, this%ham0,&
          & ints, env, this%rhoPrim, this%ErhoPrim, coordAll, errStatus, this%Dsqr, this%Qsqr)
      @:PROPAGATE_ERROR(errStatus)
    end if

    call getChargeDipole(this, this%deltaQ, this%qq, this%multipole, this%dipole, q0,&
        & this%rho, this%Ssqr, this%Dsqr, this%Qsqr, coord, iSquare, eFieldScaling, this%qBlock,&
        & this%qNetAtom, errStatus, &
        & neighbourList%iNeighbour, nNeighbourSK, orb, iSparseStart, img2CentCell,  env, ints)
    @:PROPAGATE_ERROR(errStatus)
    if (allocated(this%dispersion)) then
      call this%dispersion%updateOnsiteCharges(this%qNetAtom, orb, referenceN0,&
          & this%speciesAll(:this%nAtom), .true.)
    end if

    call updateH(this, this%H1, ints, this%ham0, this%speciesAll, this%qq, q0, coord, orb,&
        & this%potential, neighbourList, nNeighbourSK, iSquare, iSparseStart, img2CentCell, 0,&
        & this%chargePerShell, spinW, env, tDualSpinOrbit, xi, thirdOrd, this%qBlock, dftbU,&
        & onSiteElements, refExtPot, this%deltaRho, this%H1LC, this%Ssqr, solvation, hybridXc,&
        & this%dispersion,this%rho, errStatus)
    @:PROPAGATE_ERROR(errStatus)

    if (this%tForces) then
      call getForces(this, this%movedAccel, this%totalForce, this%rho, this%H1, this%Sinv,&
          & neighbourList, nNeighbourSK, img2CentCell, iSparseStart, iSquare, this%potential, orb,&
          & skHamCont,  skOverCont, this%qq, q0, repulsive, coordAll, this%rhoPrim, this%ErhoPrim,&
          & 0, env, hybridXc, this%deltaRho, this%Ssqr, errStatus)
      @:PROPAGATE_ERROR(errStatus)
    end if

    this%tPropagatorsInitialized = .true.

  end subroutine initializeDynamics


  !> Do one TD step, propagating electrons and nuclei (if IonDynamics is enabled)
  subroutine doTdStep(this, boundaryCond, iStep, coord, orb, neighbourList, nNeighbourSK,&
      & iSquare, iSparseStart, img2CentCell, skHamCont, skOverCont, ints, env, coordAll, q0,&
      & referenceN0, spinW, tDualSpinOrbit, xi, thirdOrd, dftbU, onSiteElements, refExtPot,&
      & solvation, eFieldScaling, hybridXc, repulsive, iAtInCentralRegion, tFixEf, Ef,&
      & electronicSolver, qDepExtPot, errStatus)

    !> ElecDynamics instance
    type(TElecDynamics), intent(inout), target :: this

    !> Boundary conditions on the calculation
    type(TBoundaryConditions), intent(in) :: boundaryCond

    !> Current step of the propagation
    integer, intent(in) :: iStep

    !> Reference atomic occupations
    real(dp), intent(inout) :: q0(:,:,:)

    !> Reference charges from the Slater-Koster file
    real(dp), intent(in) :: referenceN0(:,:)

    !> Integral container
    type(TIntegral), intent(inout) :: ints

    !> Atomic coordinates
    real(dp), allocatable, intent(inout) :: coord(:,:)

    !> All atomic coordinates
    real(dp), allocatable, intent(inout) :: coordAll(:,:)

    !> Spin constants
    real(dp), allocatable, intent(in) :: spinW(:,:,:)

    !> Number of neighbours for each of the atoms
    integer, intent(inout) :: nNeighbourSK(:)

    !> Index array for location of atomic blocks in large sparse arrays
    integer, allocatable, intent(inout) :: iSparseStart(:,:)

    !> Image atoms to their equivalent in the central cell
    integer, allocatable, intent(inout) :: img2CentCell(:)

    !> Index array for start of atomic block in dense matrices
    integer, intent(in) :: iSquare(:)

    !> List of neighbours for each atom
    type(TNeighbourList), intent(inout) :: neighbourList

    !> Repulsive information
    class(TRepulsive), allocatable, intent(inout) :: repulsive

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

    !> DFTB+U functional (if used)
    type(TDftbU), intent(in), allocatable :: dftbU

    !> Spin orbit constants if required
    real(dp), allocatable, intent(in) :: xi(:,:)

    !> 3rd order settings
    type(TThirdOrder), intent(inout), allocatable :: thirdOrd

    !> Solvation model
    class(TSolvation), allocatable, intent(inout) :: solvation

    !> Any dielectric environment scaling
    class(TScaleExtEField), intent(in) :: eFieldScaling

    !> Range separation contributions
    class(THybridXcFunc), allocatable, intent(inout) :: hybridXc

    !> Proxy for querying Q-dependant external potentials
    type(TQDepExtPotProxy), intent(inout), allocatable :: qDepExtPot

    !> Atoms over which to sum the total energies
    integer, intent(in) :: iAtInCentralRegion(:)

    !> Whether fixed Fermi level(s) should be used. (No charge conservation!)
    logical, intent(in) :: tFixEf

    !> If tFixEf is .true. contains reservoir chemical potential, otherwise the Fermi levels found
    !> from the given number of electrons
    real(dp), intent(inout) :: Ef(:)

    !> Corrections terms for on-site elements
    real(dp), intent(in), allocatable :: onSiteElements(:,:,:,:)

    !> Reference external potential (usual provided via API)
    type(TRefExtPot) :: refExtPot

    !> Electronic solver information
    type(TElectronicSolver), intent(inout) :: electronicSolver

    !> Error status
    type(TStatus), intent(inout) :: errStatus

    real(dp), allocatable :: velInternal(:,:)
    real(dp) :: new3Coord(3, this%nMovedAtom)
    character(sc) :: dumpIdx
    logical :: tProbeFrameWrite
    integer :: iKS

    this%time = iStep * this%dt + this%startTime

    if (this%tIons) then

      if (.not. this%tdCoordsAndVelosThroughAPI) then
        ! update coordNew (saved for later), get velocities for current step
        new3Coord(:,:) = this%coordNew(:, this%indMovedAtom)
        call next(this%pMDIntegrator, this%movedAccel, new3Coord, this%movedVelo)
        this%coordNew(:, this%indMovedAtom) = new3Coord

      elseif (.not. this%tdCoordsAndVelosAreSet) then
        @:RAISE_ERROR(errStatus, -1, "Coordinates and velocities were not set externally.")
      end if

      select case(this%hamiltonianType)
      case default
        @:ASSERT(.false.)
      case(hamiltonianTypes%dftb)
        call getRdotSprime(this, this%RdotSprime, coordAll, skOverCont, orb, img2CentCell, &
            &neighbourList, nNeighbourSK, iSquare)
      case(hamiltonianTypes%xtb)
        @:RAISE_ERROR(errStatus, -1, "Nuclei dynamic not implemented for xTB Hamiltonian yet")
      end select
      if ((this%tPopulations) .and. (mod(iStep, this%writeFreq) == 0)) then
        call updateBasisMatrices(this, env, electronicSolver, this%Eiginv, this%EiginvAdj, this%H1,&
            & this%Ssqr, errStatus)
        @:PROPAGATE_ERROR(errStatus)
      end if

      call getPositionDependentEnergy(this, env, this%energy, coordAll, img2CentCell,&
          & neighbourList, repulsive, iAtInCentralRegion, hybridXc)
    end if

    call getTDEnergy(this, env, this%energy, this%rhoPrim, this%rho, neighbourList, nNeighbourSK,&
        & orb, iSquare, iSparseStart, img2CentCell, this%ham0, this%qq, q0, this%potential,&
        & this%chargePerShell, this%energyKin, tDualSpinOrbit, thirdOrd, solvation, hybridXc,&
        & qDepExtPot, this%qBlock, dftbU, xi, iAtInCentralRegion, tFixEf, Ef, onSiteElements,&
        & errStatus)
    @:PROPAGATE_ERROR(errStatus)

    if ((mod(iStep, this%writeFreq) == 0)) then
      call getBondPopulAndEnergy(this, this%bondWork, this%lastBondPopul, this%rhoPrim, this%ham0,&
          & ints, neighbourList%iNeighbour, nNeighbourSK, iSparseStart, img2CentCell, iSquare,&
          & this%fdBondEnergy, this%fdBondPopul, this%time)
    end if

    do iKS = 1, this%parallelKS%nLocalKS
      if (this%tIons .or. (.not. this%tRealHS) .or. this%isHybridXc) then
        this%H1(:,:,iKS) = this%RdotSprime + imag * this%H1(:,:,iKS)

        if (this%tEulers .and. (iStep > 0) .and. (mod(iStep, max(this%eulerFreq,1)) == 0)) then
          call zcopy(this%nOrbs*this%nOrbs, this%rho(:,:,iKS), 1, this%rhoOld(:,:,iKS), 1)
          call propagateRho(this, this%rhoOld(:,:,iKS), this%rho(:,:,iKS),&
              & this%H1(:,:,iKS), this%Sinv(:,:,iKS), this%dt)
        else
          call propagateRho(this, this%rhoOld(:,:,iKS), this%rho(:,:,iKS),&
              & this%H1(:,:,iKS), this%Sinv(:,:,iKS), 2.0_dp * this%dt)
        end if
      else
#:if WITH_SCALAPACK
        call propagateRhoRealHBlacs(this, this%rhoOld(:,:,iKS), this%rho(:,:,iKS),&
      & this%H1(:,:,iKS), this%Sinv(:,:,iKS), 2.0_dp * this%dt)
#:else
        call propagateRhoRealH(this, this%rhoOld(:,:,iKS), this%rho(:,:,iKS),&
            & this%H1(:,:,iKS), this%Sinv(:,:,iKS), 2.0_dp * this%dt)
#:endif
      end if
    end do

    if (mod(iStep, 2) == 0) then
      this%rho => this%trho
      this%rhoOld => this%trhoOld
    else
      this%rho => this%trhoOld
      this%rhoOld => this%trho
    end if

    if ((this%tPopulations) .and. (mod(iStep, this%writeFreq) == 0)) then
      do iKS = 1, this%parallelKS%nLocalKS
        ! time-dt is due to the fact that populations were always written one step later than the
        ! rest of the quantities but with the same time label.
        ! TODO: fix tests values for populations so that it becomes exactly syncronized with the
        ! other outputs
        call getTDPopulations(this, this%occ, this%rho, this%Eiginv, this%EiginvAdj, this%populDat,&
            & this%time-this%dt, iKS)
      end do
    end if

    if (.not. this%tReadRestart .or. (iStep > 0) .or. this%tProbe) then
      call writeTDOutputs(this, this%dipoleDat, this%qDat, this%energyDat, &
          & this%forceDat, this%coorDat, this%fdBondPopul, this%fdBondEnergy, this%atomEnergyDat,&
          & this%time, this%energy, this%energyKin, this%dipole, this%deltaQ, coord,&
          & this%totalForce, iStep)
    end if

    if (this%tWriteRestart .and. iStep > 0 .and. mod(iStep, max(this%restartFreq,1)) == 0) then
      allocate(velInternal(3,size(this%movedVelo, dim=2)))
      if (this%tIons) then
        call state(this%pMDIntegrator, velocities=velInternal)
      else
        velInternal(:,:) = 0.0_dp
      end if
      call writeRestartFile(this%rho, this%rhoOld, coord, velInternal, this%time, this%dt, &
          &restartFileName, this%tWriteRestartAscii, errStatus)
      @:PROPAGATE_ERROR(errStatus)
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
      call writeRestartFile(this%rho, this%rhoOld, coord, velInternal, this%time, this%dt,&
          & trim(pumpFilesDir) // '/' // trim(dumpIdx) // 'ppdump', this%tWriteRestartAscii,&
          & errStatus)
      @:PROPAGATE_ERROR(errStatus)
      deallocate(velInternal)
    end if

    if (this%tIons) then
      coord(:,:) = this%coordNew
      call updateH0S(this, boundaryCond, this%Ssqr, this%Sinv, coord, orb, neighbourList,&
          & nNeighbourSK, iSquare, iSparseStart, img2CentCell, skHamCont, skOverCont, this%ham0,&
          & ints, env, this%rhoPrim, this%ErhoPrim, coordAll, errStatus, this%Dsqr, this%Qsqr)
      @:PROPAGATE_ERROR(errStatus)
    end if

    call getChargeDipole(this, this%deltaQ, this%qq, this%multipole, this%dipole, q0,&
        & this%rho, this%Ssqr, this%Dsqr, this%Qsqr, coord, iSquare, eFieldScaling, this%qBlock,&
        & this%qNetAtom, errStatus, neighbourList%iNeighbour, nNeighbourSK, orb, iSparseStart,&
        & img2CentCell, env, ints)
    @:PROPAGATE_ERROR(errStatus)
    if (allocated(this%dispersion)) then
      call this%dispersion%updateOnsiteCharges(this%qNetAtom, orb, referenceN0,&
          & this%speciesAll(:this%nAtom), .true.)
    end if

    call updateH(this, this%H1, ints, this%ham0, this%speciesAll, this%qq, q0, coord, orb,&
        & this%potential, neighbourList, nNeighbourSK, iSquare, iSparseStart, img2CentCell, iStep,&
        & this%chargePerShell, spinW, env, tDualSpinOrbit, xi, thirdOrd, this%qBlock, dftbU,&
        & onSiteElements, refExtPot, this%deltaRho, this%H1LC, this%Ssqr, solvation, hybridXc,&
        & this%dispersion,this%rho, errStatus)
    @:PROPAGATE_ERROR(errStatus)

    if (this%tForces) then
      call getForces(this, this%movedAccel, this%totalForce, this%rho, this%H1, this%Sinv,&
          & neighbourList, nNeighbourSK, img2CentCell, iSparseStart, iSquare, this%potential, orb,&
          & skHamCont, skOverCont, this%qq, q0, repulsive, coordAll, this%rhoPrim, this%ErhoPrim,&
          & iStep, env, hybridXc, this%deltaRho, this%Ssqr, errStatus)
      @:PROPAGATE_ERROR(errStatus)
    end if

    ! unset coordinates and velocities at the end of the step
    if (this%tdCoordsAndVelosThroughAPI) then
      this%tdCoordsAndVelosAreSet = .false.
    end if

  end subroutine doTdStep


  !> Handles deallocation, closing outputs and autotest writing
  subroutine finalizeDynamics(this)

    !> ElecDynamics instance
    type(TElecDynamics), intent(inout) :: this

    call closeTDOutputs(this%dipoleDat, this%qDat, this%energyDat, this%populDat,&
        & this%forceDat, this%coorDat, this%fdBondPopul, this%fdBondEnergy, this%atomEnergyDat)

    deallocate(this%Ssqr)
    deallocate(this%Sinv)
    deallocate(this%H1)
    deallocate(this%trho)
    deallocate(this%trhoOld)
    deallocate(this%RdotSprime)
    deallocate(this%qq)
    deallocate(this%deltaQ)
    deallocate(this%dipole)
    deallocate(this%chargePerShell)
    deallocate(this%occ)
    deallocate(this%totalForce)
    if (allocated(this%Dsqr)) then
      deallocate(this%Dsqr)
    end if
    if (allocated(this%Qsqr)) then
      deallocate(this%Qsqr)
    end if

    deallocate(this%rhoPrim)
    deallocate(this%ErhoPrim)
    if (allocated(this%H1LC)) then
      deallocate(this%H1LC)
    end if
    if (allocated(this%deltaRho)) then
      deallocate(this%deltaRho)
    end if
    if (allocated(this%qBlock)) then
      deallocate(this%qBlock)
    end if
    if (allocated(this%qNetAtom)) then
      deallocate(this%qNetAtom)
    end if
    if (this%tPopulations) then
      deallocate(this%Eiginv)
      deallocate(this%EiginvAdj)
    end if
    if (this%tBondE .or. this%tBondP) then
      deallocate(this%bondWork)
    end if

    if (this%tIons) then
      if (allocated(this%polDirs)) then
        if (size(this%polDirs) < (this%nDynamicsInit + 1)) then
          deallocate(this%pMDIntegrator)
        end if
      end if
    end if

  end subroutine finalizeDynamics

end module dftbp_timedep_timeprop
