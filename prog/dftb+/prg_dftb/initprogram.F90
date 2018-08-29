!--------------------------------------------------------------------------------------------------!
!  DFTB+: general package for performing fast atomistic simulations                                !
!  Copyright (C) 2018  DFTB+ developers group                                                      !
!                                                                                                  !
!  See the LICENSE file for terms of usage and distribution.                                       !
!--------------------------------------------------------------------------------------------------!

#:include 'common.fypp'

!> Global variables and initialization for the main program.
module initprogram
  use omp_lib
  use mainio, only : initOutputFile
  use assert
  use globalenv
  use environment
  use scalapackfx
  use inputdata_module
  use densedescr
  use constants
  use periodic
  use accuracy
  use intrinsicpr
  use shortgamma
  use coulomb
  use message
  use mixer
  use simplemixer
  use andersonmixer
  use broydenmixer
  use diismixer

  use geoopt
  use conjgrad
  use steepdesc
  use gdiis
  use lbfgs

  use randomgenpool
  use ranlux
  use mdcommon
  use mdintegrator
  use velocityverlet
  use thermostat
  use dummytherm
  use andersentherm
  use berendsentherm
  use nhctherm
  use tempprofile
  use numderivs2
  use lapackroutines
  use simplealgebra
  use nonscc
  use scc
  use h5correction
  use sccinit
  use slakocont
  use repcont
  use spin, only: Spin_getOrbitalEquiv, ud2qm, qm2ud
  use dftbplusu
  use dispersions
  use thirdorder_module
  use linresp_module
  use stress
  use orbitalequiv
  use commontypes
  use sorting, only : heap_sort
  use linkedlist
  use xlbomd_module
  use etemp, only : Fermi
#:if WITH_SOCKETS
  use mainio, only : receiveGeometryFromSocket
  use ipisocket
#:endif
  use elstatpot
  use pmlocalisation
  use energies
  use potentials
  use taggedoutput
  use formatout
  implicit none

  !> Tagged output files (machine readable)
  character(*), parameter :: autotestTag = "autotest.tag"

  !> Detailed user output
  character(*), parameter :: userOut = "detailed.out"

  !> band structure and filling information
  character(*), parameter :: bandOut = "band.out"

  !> File accumulating data during an MD run
  character(*), parameter :: mdOut = "md.out"

  !> Machine readable tagged output
  character(*), parameter :: resultsTag = "results.tag"

  !> Second derivative of the energy with respect to atomic positions
  character(*), parameter :: hessianOut = "hessian.out"

  !> file name prefix for charge data
  character(*), parameter :: fCharges = "charges"

  !> file to stop code during geometry driver
  character(*), parameter :: fStopDriver = "stop_driver"

  !> file to stop code during scc cycle
  character(*), parameter :: fStopSCC = "stop_scc"


  !> Is the calculation SCC?
  logical :: tSccCalc

  !> SCC module internal variables
  type(TScc), allocatable :: sccCalc

  !> nr. of atoms
  integer :: nAtom

  !> nr. of all (boundary condition images and original) atoms
  integer :: nAllAtom

  !> nr. of original atom in central cell
  integer, allocatable :: Img2CentCell(:)

  !> nr of different types (nAtom)
  integer :: nType


  !> data type for atomic orbital information
  type(TOrbitals), target :: orb

  !> nr. of orbitals in the system
  integer :: nOrb

  !> nr. of orbitals for all atoms
  integer :: nAllOrb

  !> types of the atoms (nAllAtom)
  integer, allocatable :: species(:)

  !> type of the atoms (nAtom)
  integer, allocatable :: species0(:)

  !> Coords of the atoms (3, nAllAtom)
  real(dp), allocatable :: coord(:,:)

  !> Coords in central cell (3, nAtom)
  real(dp), allocatable, target :: coord0(:,:)

  !> temporary coordinates
  real(dp), allocatable :: tmpCoords(:)

  !> temporary weights
  real(dp), allocatable :: tmpWeight(:)

  !> temporary array of coords (3,:)
  real(dp), allocatable :: tmp3Coords(:,:)


  !> if calculation is periodic
  logical :: tPeriodic

  !> Should central cell coordinates be output?
  logical :: tShowFoldedCoord


  !> How to calculate forces
  integer :: forceType

  !> are atomic coordinates fractional?
  logical :: tFracCoord

  !> Tolerance for SCC cycle
  real(dp) :: sccTol


  !> lattice vectors as columns
  real(dp), allocatable, target :: latVec(:,:)

  !> reciprocal lattice vectors as columns
  real(dp), allocatable, target :: recVec(:,:)


  !> original lattice vectors used for optimizing
  real(dp) :: origLatVec(3,3)

  !> normalized vectors in those directions
  real(dp) :: normOrigLatVec(3,3)


  !> reciprocal vectors in 2 pi units
  real(dp), allocatable :: invLatVec(:,:)

  !> cell volume
  real(dp) :: CellVol

  !> reciprocal cell volume
  real(dp) :: recCellVol

  !> translation vecs for interacting image cells (3, nImgCell + 1)
  real(dp), allocatable :: cellVec(:,:)

  !> cell vectors in absolute units
  real(dp), allocatable :: rCellVec(:,:)

  !> index in cellVec for each atom
  integer, allocatable :: iCellVec(:)


  !> ADT for neighbour parameters
  type(TNeighbourList), allocatable, save :: neighbourList

  !> nr. of neighbours for atoms out to max interaction distance (excluding Ewald terms)
  integer, allocatable :: nNeighbourSK(:)

  !> nr. of neighbours for atoms within Erep interaction distance (usually short)
  integer, allocatable :: nNeighbourRep(:)

  !> H/S sparse matrices indexing array for atomic blocks
  integer, allocatable :: iSparseStart(:,:)

  !> Hubbard Us (orbital, atom)
  real(dp), allocatable, target :: hubbU(:,:)

  !> self energy (orbital, atom)
  real(dp), allocatable :: atomEigVal(:,:)

  !> reference n_0 charges for each atom
  real(dp), allocatable :: referenceN0(:,:)

  !> list of atomic masses
  real(dp), allocatable :: mass(:)

  !> list of atomic masses for each species
  real(dp), allocatable :: speciesMass(:)

  !> Raw H^0 hamiltonian data
  type(OSlakoCont) :: skHamCont

  !> Raw overlap hamiltonian data
  type(OSlakoCont) :: skOverCont

  !> Repulsive interaction raw data
  type(ORepCont) :: pRepCont

  !> Cut off distance for Slater-Koster interactions
  real(dp) :: skCutOff

  !> Cut off distance for repulsive interactions
  real(dp) :: repCutOff

  !> longest range of interactions for which neighbours are required
  real(dp) :: mCutOff

  !> Sparse hamiltonian matrix
  real(dp), allocatable :: ham(:,:)

  !> imaginary part of the Hamiltonian
  real(dp), allocatable :: iHam(:,:)

  !> Charge per atomic shell (shell, atom, spin channel)
  real(dp), allocatable :: chargePerShell(:,:,:)

  !> sparse overlap
  real(dp), allocatable :: over(:)


  !> nr. of K-points
  integer :: nKPoint

  !> K-points
  real(dp), allocatable :: kPoint(:,:)

  !> weight of the K-Points
  real(dp), allocatable :: KWeight(:)


  !> external pressure if periodic
  real(dp) :: extPressure

  !> Barostat used if MD and periodic
  logical :: tBarostat

  !> Barostat coupling strength
  real(dp) :: BarostatStrength


  !> H and S are real
  logical :: tRealHS


  !> nr. of electrons
  real(dp), allocatable :: nEl(:)

  !> Nr. of all electrons if neutral
  real(dp) :: nEl0


  !> Spin W values
  real(dp), allocatable :: spinW(:,:,:)

  !> Spin orbit constants
  real(dp), allocatable :: xi(:,:)


  !> is this a DFTB+U calculation?
  logical :: tDFTBU

  !> Choice of orbital functional
  integer :: nDFTBUfunc

  !> list of U-J for species
  real(dp), allocatable :: UJ(:,:)

  !> How many U-J for each species
  integer, allocatable :: nUJ(:)

  !> number of l-values of U-J for each block
  integer, allocatable :: niUJ(:,:)

  !> l-values of U-J for each block
  integer, allocatable :: iUJ(:,:,:)


  !> electron temperature
  real(dp) :: tempElec

  !> If K points should filled separately
  logical :: tFillKSep

  !> Fix Fermi energy at specified value
  logical :: tFixEf

  !> Fermi energy
  real(dp) :: Ef(2)


  !> Filling temp updated by MD.
  logical :: tSetFillingTemp

  !> Choice of electron distribution function, defaults to Fermi
  integer :: iDistribFn = 0

  !> atomic kinetic temperature
  real(dp) :: tempAtom

  !> MD stepsize
  real(dp) :: deltaT


  !> eigensolver
  integer :: solver

  !> maximal number of SCC iterations
  integer :: maxSccIter


  !> is this a spin polarized calculation?
  logical :: tSpin

  !> Number of spin components, 1 is unpolarised, 2 is polarised, 4 is noncolinear / spin-orbit
  integer :: nSpin

  !> is there spin-orbit coupling
  logical :: tSpinOrbit

  !> Use block like dual representation for spin orbit
  logical :: tDualSpinOrbit

  !> Is there a complex hamiltonian contribution in real space
  logical :: tImHam

  !> is this a two component calculation (spin orbit or non-collinear spin)
  logical :: t2Component

  !> Common Fermi level accross spin channels
  logical :: tSpinSharedEf


  !> Geometry optimization needed?
  logical :: tGeoOpt

  !> optimize coordinates inside unit cell (periodic)?
  logical :: tCoordOpt

  !> optimize lattice constants?
  logical :: tLatOpt

  !> Fix angles between lattice vectors when optimizing?
  logical :: tLatOptFixAng

  !> Fix length of specified lattice vectors when optimizing?
  logical :: tLatOptFixLen(3)

  !> Optimise lattice isotropically
  logical :: tLatOptIsotropic

  !> Is this a MD calculation?
  logical :: tMD

  !> Is this a derivatives calc?
  logical :: tDerivs

  !> Do we need Mulliken charges?
  logical :: tMulliken

  !> Electrostatic potentials if requested
  type(TElStatPotentials), allocatable :: esp

  !> Calculate localised orbitals?
  logical :: tLocalise

  !> Do we need to show Mulliken charges?
  logical :: tPrintMulliken

  !> calculate an electric dipole?
  logical :: tDipole

  !> Do we need atom resolved E?
  logical :: tAtomicEnergy

  !> Print out eigenvectors?
  logical :: tPrintEigVecs

  !> Store eigenvectors as a text file
  logical :: tPrintEigVecsTxt

  !> Print eigenvector projections?
  logical :: tProjEigenvecs

  !> Do we need forces?
  logical :: tForces

  !> are forces being returned
  logical :: tPrintForces

  !> Number of moved atoms
  integer :: nMovedAtom

  !> Index of the moved atoms
  integer, allocatable :: indMovedAtom(:)

  !> Nr. of moved coordinates
  integer :: nMovedCoord

  !> Nr. of geo movements to do
  integer :: nGeoSteps

  !> Index of constrained atoms
  integer, allocatable :: conAtom(:)

  !> Constraint vectors
  real(dp), allocatable :: conVec(:,:)

  !> Pipek-Mezey localisation calculator
  type(TPipekMezey), allocatable :: pipekMezey

  !> use commands from socket communication to control the run
  logical :: tSocket

  !> socket details
#:if WITH_SOCKETS
  type(IpiSocketComm), allocatable :: socket
#:endif

  !> File containing output geometry
  character(lc) :: geoOutFile


  !> Append geometries in the output?
  logical :: tAppendGeo


  !> Only use converged forces if SCC
  logical :: tUseConvergedForces


  !> labels of atomic species
  character(mc), allocatable :: speciesName(:)

  !> General geometry optimizer
  type(OGeoOpt), allocatable :: pGeoCoordOpt

  !> Geometry optimizer for lattice consts
  type(OGeoOpt), allocatable :: pGeoLatOpt


  !> Charge mixer
  type(OMixer), allocatable :: pChrgMixer


  !> MD Framework
  type(OMDCommon), allocatable :: pMDFrame

  !> MD integrator
  type(OMDIntegrator), allocatable :: pMDIntegrator

  !> Temperature profile driver in MD
  type(OTempProfile), allocatable, target :: temperatureProfile

  !> geometry optimiser
  type(OnumDerivs), allocatable, target :: derivDriver

  !> reference neutral atomic occupations
  real(dp), allocatable :: q0(:, :, :)

  !> shell resolved neutral reference
  real(dp), allocatable :: qShell0(:,:)

  !> input charges (for potentials)
  real(dp), allocatable :: qInput(:, :, :)

  !> output charges
  real(dp), allocatable :: qOutput(:, :, :)

  !> input Mulliken block charges (diagonal part == Mulliken charges)
  real(dp), allocatable :: qBlockIn(:, :, :, :)

  !> Output Mulliken block charges
  real(dp), allocatable :: qBlockOut(:, :, :, :)

  !> Imaginary part of input Mulliken block charges
  real(dp), allocatable :: qiBlockIn(:, :, :, :)

  !> Imaginary part of output Mulliken block charges
  real(dp), allocatable :: qiBlockOut(:, :, :, :)

  !> input charges packed into unique equivalence elements
  real(dp), allocatable :: qInpRed(:)

  !> output charges packed into unique equivalence elements
  real(dp), allocatable :: qOutRed(:)

  !> charge differences packed into unique equivalence elements
  real(dp), allocatable :: qDiffRed(:)

  !> Orbital equivalence relations
  integer, allocatable :: iEqOrbitals(:,:,:)

  !> nr. of inequivalent orbitals
  integer :: nIneqOrb

  !> nr. of elements to go through the mixer - may contain reduced orbitals and also orbital blocks
  !> (if tDFTBU)
  integer :: nMixElements

  !> Orbital equivalency for orbital blocks
  integer, allocatable :: iEqBlockDFTBU(:,:,:,:)

  !> Orbital equivalency for orbital blocks with spin-orbit
  integer, allocatable :: iEqBlockDFTBULS(:,:,:,:)

  ! External charges

  !> If external charges must be considered
  logical :: tExtChrg

  !> Nr. of external charges
  integer :: nExtChrg


  !> external electric field
  logical :: tEField = .false.

  !> field strength
  real(dp) :: EFieldStrength = 0.0_dp

  !> field direction
  real(dp) :: EfieldVector(3) = 0.0_dp

  !> time dependent
  logical :: tTDEfield = .false.

  !> angular frequency
  real(dp) :: EfieldOmega = 0.0_dp

  !> phase of field at step 0
  integer :: EfieldPhase = 0


  !> Partial density of states (PDOS) projection regions
  type(listIntR1), save :: iOrbRegion

  !> PDOS region labels
  type(listCharLc), save :: regionLabels

  !> Third order DFTB
  logical :: t3rd

  !> Full 3rd order or only atomic site
  logical :: t3rdFull

  !> data structure for 3rd order
  type(ThirdOrder), allocatable :: thirdOrd


  !> Calculate Casida linear response excitations
  logical :: tLinResp

  !> calculate Z vector for excited properties
  logical :: tLinRespZVect

  !> Print eigenvectors
  logical :: tPrintExcitedEigVecs = .false.

  !> data type for linear response
  type(linresp), save :: lresp

  !> If initial charges/dens mtx. from external file.
  logical :: tReadChrg

  !> should charges written to disc be in ascii or binary format?
  logical :: tWriteChrgAscii

  !> produce tagged output?
  logical :: tWriteAutotest

  !> Produce detailed.xml
  logical :: tWriteDetailedXML

  !> Produce detailed.tag
  logical :: tWriteResultsTag

  !> Produce detailed.out
  logical :: tWriteDetailedOut

  !> Produce band.dat
  logical :: tWriteBandDat

  !> Should HS (square) be printed?
  logical :: tWriteHS

  !> Should HS (sparse) be printed?
  logical :: tWriteRealHS

  !> Program run id
  integer :: runId

  !> Frequency for saving restart info
  integer :: restartFreq

  !> If dispersion should be calculated
  logical :: tDispersion

  !> dispersion data and calculations
  class(DispersionIface), allocatable :: dispersion

  !> Can stress be calculated? - start by assuming it can
  logical :: tStress = .true.

  !> should XLBOMD be used in MD
  logical :: tXlbomd

  !> XLBOMD related parameters
  type(Xlbomd), allocatable :: xlbomdIntegrator

  !> Differentiation method for (H^0,S)
  type(NonSccDiff), save :: nonSccDeriv

  !> Whether lattice has changed since last geometry iteration
  logical :: tLatticeChanged

  !> Whether atomic coordindates have changed since last geometry iteration
  logical :: tCoordsChanged

  !> Dense matrix descriptor for H and S
  type(TDenseDescr) :: denseDesc

  !> MD velocities
  real(dp), allocatable :: velocities(:,:)

  !> MD velocities for moved atoms
  real(dp), allocatable :: movedVelo(:,:)

  !> MD acceleration for moved atoms
  real(dp), allocatable :: movedAccel(:,:)

  !> Mass of the moved atoms
  real(dp), allocatable :: movedMass(:,:)

  !> Sparse storage of density matrix
  real(dp), allocatable :: rhoPrim(:,:)

  !> Imaginary part of density matrix in sparse storage
  real(dp), allocatable :: iRhoPrim(:,:)

  !> Energy weighted density matrix
  real(dp), allocatable :: ERhoPrim(:)

  !> Non-SCC part of the hamiltonian in sparse storage
  real(dp), allocatable :: h0(:)

  !> electronic filling
  real(dp), allocatable :: filling(:,:,:)

  !> band structure energy
  real(dp), allocatable :: Eband(:)

  !> entropy of electrons at temperature T
  real(dp), allocatable :: TS(:)

  !> zero temperature electronic energy
  real(dp), allocatable :: E0(:)

  !> Square dense hamiltonian storage for cases with k-points
  complex(dp), allocatable :: HSqrCplx(:,:)

  !> Square dense overlap storage for cases with k-points
  complex(dp), allocatable :: SSqrCplx(:,:)

  !> Complex eigenvectors
  complex(dp), allocatable :: eigvecsCplx(:,:,:)

  !> Square dense hamiltonian storage
  real(dp), allocatable :: HSqrReal(:,:)

  !> Square dense overlap storage
  real(dp), allocatable :: SSqrReal(:,:)

  !> Real eigenvectors
  real(dp), allocatable :: eigvecsReal(:,:,:)

  !> Eigenvalues
  real(dp), allocatable :: eigen(:,:,:)

  !> density matrix
  real(dp), allocatable :: rhoSqrReal(:,:,:)

  !> Total energy components
  type(TEnergies) :: energy

  !> Potentials for orbitals
  type(TPotentials) :: potential

  !> Energy derivative with respect to atomic positions
  real(dp), allocatable :: derivs(:,:)

  !> Forces on any external charges
  real(dp), allocatable :: chrgForces(:,:)

  !> excited state force addition
  real(dp), allocatable :: excitedDerivs(:,:)

  !> dipole moments when available
  real(dp), allocatable :: dipoleMoment(:)

  !> Coordinates to print out
  real(dp), pointer :: pCoord0Out(:,:)

  !> Folded coords (3, nAtom)
  real(dp), allocatable, target :: coord0Fold(:,:)

  !> New coordinates returned by the MD routines
  real(dp), allocatable :: newCoords(:,:)

  !> Orbital angular momentum
  real(dp), allocatable :: orbitalL(:,:,:)

  !> Natural orbitals for excited state density matrix, if requested
  real(dp), allocatable, target :: occNatural(:)

  !> Dynamical (Hessian) matrix
  real(dp), pointer :: pDynMatrix(:,:)

  !> File descriptor for the human readable output
  integer :: fdDetailedOut

  !> File descriptor for extra MD output
  integer :: fdMD

  !> Contains (iK, iS) tuples to be processed in parallel by various processor groups
  type(TParallelKS) :: parallelKS

  private :: createRandomGenerators

contains


  !> Initializes the variables in the module based on the parsed input
  subroutine initProgramVariables(input, env)

    !> Holds the parsed input data.
    type(inputData), intent(inout), target :: input

    !> Environment settings
    type(TEnvironment), intent(inout) :: env

    ! Mixer related local variables
    integer :: nGeneration
    real(dp) :: mixParam

    !> mixer number
    integer :: iMixer

    !> simple mixer (if used)
    type(OSimpleMixer), allocatable :: pSimpleMixer

    !> Anderson mixer (if used)
    type(OAndersonMixer), allocatable :: pAndersonMixer

    !> Broyden mixer (if used)
    type(OBroydenMixer), allocatable :: pBroydenMixer

    !> DIIS mixer (if used)
    type(ODIISMixer), allocatable :: pDIISMixer

    ! Geometry optimizer related local variables

    !> Conjugate gradient driver
    type(OConjGrad), allocatable :: pConjGrad

    !> Steepest descent driver
    type(OSteepDesc), allocatable :: pSteepDesc

    !> Conjugate gradient driver
    type(OConjGrad), allocatable :: pConjGradLat

    !> Steepest descent driver
    type(OSteepDesc), allocatable :: pSteepDescLat

    !> gradient DIIS driver
    type(ODIIS), allocatable :: pDIIS

    !> lBFGS driver for geometry  optimisation
    type(TLbfgs), allocatable :: pLbfgs

    !> lBFGS driver for lattice optimisation
    type(TLbfgs), allocatable :: pLbfgsLat

    ! MD related local variables
    type(OThermostat), allocatable :: pThermostat
    type(ODummyThermostat), allocatable :: pDummyTherm
    type(OAndersenThermostat), allocatable :: pAndersenTherm
    type(OBerendsenThermostat), allocatable :: pBerendsenTherm
    type(ONHCThermostat), allocatable :: pNHCTherm

    type(OVelocityVerlet), allocatable :: pVelocityVerlet
    type(OTempProfile), pointer :: pTempProfile

    type(ORanlux), allocatable :: randomInit, randomThermostat
    integer :: iSeed

    integer :: ind, ii, jj, kk, iS, iAt, iSp, iSh, iOrb
    integer :: iStart, iEnd

    ! Dispersion
    type(DispSlaKirk), allocatable :: slaKirk
    type(DispUFF), allocatable :: uff
  #:if WITH_DFTD3
    type(DispDftD3), allocatable :: dftd3
  #:endif

    ! H5 correction
    type(H5Corr), allocatable :: pH5Correction
    logical :: tHHRepulsion

    character(lc) :: tmpStr
    integer, allocatable :: tmpir1(:)

    character(lc) :: strTmp, strTmp2

    !> flag to check for first cycle through a loop
    logical :: tFirst

    !> Nr. of Hamiltonians to diagonalise independently
    integer :: nIndepHam

    real(dp) :: rTmp

    !> Flag if some files do exist or not
    logical :: tExist

    ! Orbital equivalency for SCC and Spin
    integer, allocatable :: iEqOrbSCC(:,:,:), iEqOrbSpin(:,:,:)
    ! Orbital equivalency for orbital potentials
    integer, allocatable :: iEqOrbDFTBU(:,:,:)

    ! Damped interactions
    logical, allocatable, target :: tDampedShort(:)
    type(ThirdOrderInp) :: thirdInp

    ! PDOS stuff
    integer :: iReg, nAtomRegion, nOrbRegion, iTmp
    integer, allocatable :: iAtomRegion(:)
    integer :: valshape(1)

    !> Is SCC cycle initialised
    type(TSccInp), allocatable :: sccInp

    !> Used for indexing linear response
    integer :: homoLoc(1)

    !> Whether seed was randomly created
    logical :: tRandomSeed

    !> First guess for nr. of neighbours.
    integer, parameter :: nInitNeighbour = 40


    @:ASSERT(input%tInitialized)

    write(stdOut, "(/, A)") "Starting initialization..."
    write(stdOut, "(A80)") repeat("-", 80)

    call env%initGlobalTimer(input%ctrl%timingLevel, "DFTB+ running times", stdOut)
    call env%globalTimer%startTimer(globalTimers%globalInit)

    ! Basic variables
    tSccCalc = input%ctrl%tScc
    tDFTBU = input%ctrl%tDFTBU
    tSpin = input%ctrl%tSpin
    if (tSpin) then
      nSpin = 2
    else
      nSpin = 1
    end if
    tSpinSharedEf = input%ctrl%tSpinSharedEf
    tSpinOrbit = input%ctrl%tSpinOrbit
    tDualSpinOrbit = input%ctrl%tDualSpinOrbit
    t2Component = input%ctrl%t2Component

    if (t2Component) then
      nSpin = 4
      nIndepHam = 1
    else
      nIndepHam = nSpin
    end if
    if (nSpin /= 2 .and. tSpinSharedEf) then
      call error("Colinear spin polarization required for shared Ef over spin channels")
    end if

    nAtom = input%geom%nAtom
    nType = input%geom%nSpecies
    orb = input%slako%orb
    nOrb = orb%nOrb
    tPeriodic = input%geom%tPeriodic

    ! Brillouin zone sampling
    if (tPeriodic) then
      nKPoint = input%ctrl%nKPoint
      allocate(kPoint(3, nKPoint))
      allocate(kWeight(nKPoint))
      @:ASSERT(all(shape(kPoint) == shape(input%ctrl%KPoint)))
      @:ASSERT(all(shape(kWeight) == shape(input%ctrl%kWeight)))
      kPoint(:,:) = input%ctrl%KPoint(:,:)
      if (sum(input%ctrl%kWeight(:)) < epsilon(1.0_dp)) then
        call error("Sum of k-point weights should be greater than zero!")
      end if
      kWeight(:) = input%ctrl%kWeight / sum(input%ctrl%kWeight)
    else
      nKPoint = 1
      allocate(kPoint(3, nKPoint))
      allocate(kWeight(nKpoint))
      kPoint(:,1) = 0.0_dp
      kWeight(1) = 1.0_dp
    end if

    if ((.not. tPeriodic) .or. (nKPoint == 1 .and. all(kPoint(:, 1) == [0.0_dp, 0.0_dp, 0.0_dp])))&
        & then
      tRealHS = .true.
    else
      tRealHS = .false.
    end if

  #:if WITH_MPI
    call env%initMpi(input%ctrl%parallelOpts%nGroup)
  #:endif
  #:if WITH_SCALAPACK
    call initScalapack(input%ctrl%parallelOpts%blacsOpts, nAtom, nOrb, t2Component, env)
  #:endif
    call TParallelKS_init(parallelKS, env, nKPoint, nIndepHam)

    sccTol = input%ctrl%sccTol
    tShowFoldedCoord = input%ctrl%tShowFoldedCoord
    if (tShowFoldedCoord .and. .not. tPeriodic) then
      call error("Folding coordinates back into the central cell is meaningless for molecular&
          & boundary conditions!")
    end if
    tFracCoord = input%geom%tFracCoord
    solver = input%ctrl%iSolver
    if (tSccCalc) then
      maxSccIter = input%ctrl%maxIter
    else
      maxSccIter = 1
    end if

    if (tPeriodic) then
      tLatticeChanged = .true.
      allocate(latVec(3, 3))
      @:ASSERT(all(shape(input%geom%latVecs) == shape(latVec)))
      latVec(:,:) = input%geom%latVecs(:,:)
      allocate(recVec(3, 3))
      allocate(invLatVec(3, 3))
      invLatVec = latVec(:,:)
      call matinv(invLatVec)
      invLatVec = reshape(invLatVec, (/3, 3/), order=(/2, 1/))
      recVec = 2.0_dp * pi * invLatVec
      CellVol = abs(determinant33(latVec))
      recCellVol = abs(determinant33(recVec))
    else
      allocate(latVec(0, 0))
      allocate(recVec(0, 0))
      allocate(invLatVec(0, 0))
      CellVol = 0.0_dp
      recCellVol = 0.0_dp
      tLatticeChanged = .false.
    end if

    ! Slater-Koster tables
    skHamCont = input%slako%skHamCont
    skOverCont = input%slako%skOverCont
    pRepCont = input%slako%repCont

    allocate(atomEigVal(orb%mShell, nType))
    @:ASSERT(size(input%slako%skSelf, dim=1) == orb%mShell)
    @:ASSERT(size(input%slako%skSelf, dim=2) == size(atomEigVal, dim=2))
    atomEigVal(:,:) = input%slako%skSelf(1:orb%mShell, :)

    @:ASSERT(size(input%slako%skOcc, dim=1) >= orb%mShell)
    allocate(referenceN0(orb%mShell, nType))
    referenceN0(:,:) = input%slako%skOcc(1:orb%mShell, :)
    @:ASSERT(size(input%slako%mass) == nType)
    allocate(speciesMass(nType))
    speciesMass(:) = input%slako%mass(:)

    ! Spin W's !'
    if (allocated(input%ctrl%spinW)) then
      allocate(spinW(orb%mShell, orb%mShell, nType))
      spinW(:,:,:) = 0.0_dp
      do iSp = 1, nType
        do jj = 1, orb%nShell(iSp)
          do kk = 1, orb%nShell(iSp)
            spinW(jj, kk, iSp) = input%ctrl%spinW(jj, kk, iSp)
          end do
        end do
      end do
    end if

    if (tSpinOrbit) then
      allocate(xi(orb%mShell,nType))
      xi(:,:) = 0.0_dp
      do iSp=1,nType
        do jj=1, orb%nShell(iSp)
          xi(jj,iSp)=input%ctrl%xi(jj,iSp)
        end do
      end do
    end if

    ! DFTB+U parameters
    if (tDFTBU) then
      nDFTBUfunc = input%ctrl%DFTBUfunc
      allocate(UJ(size(input%ctrl%UJ,dim=1),size(input%ctrl%UJ,dim=2)))
      allocate(nUJ(size(input%ctrl%nUJ)))
      allocate(niUJ(size(input%ctrl%niUJ,dim=1),size(input%ctrl%niUJ,dim=2)))
      allocate(iUJ(size(input%ctrl%iUJ,dim=1), size(input%ctrl%iUJ,dim=2),&
          & size(input%ctrl%iUJ,dim=3)))

      UJ(:,:) = input%ctrl%UJ(:,:)
      nUJ(:) = input%ctrl%nUJ(:)
      niUJ(:,:) = input%ctrl%niUJ(:,:)
      iUJ(:,:,:) = input%ctrl%iUJ(:,:,:)
      do iSp = 1, nType
        do jj = 1, nUJ(iSp)
          if (niUJ(jj,iSp)>1) then
            call heap_sort(iUJ(1:niUJ(jj,iSp),jj,iSp))
          end if
        end do
      end do
    else
      allocate(UJ(0,0))
      allocate(nUJ(0))
      allocate(niUJ(0,0))
      allocate(iUJ(0,0,0))
    end if

    ! Cut-offs from SlaKo and repulsive
    skCutOff = max(getCutOff(skHamCont), getCutOff(skOverCont))
    repCutOff = getCutOff(pRepCont)
    mCutOff = max(skCutOff, repCutOff)

    ! Get species names and output file
    geoOutFile = input%ctrl%outFile
    allocate(speciesName(size(input%geom%speciesNames)))
    speciesName(:) = input%geom%speciesNames(:)

    do iSp = 1, nType
      do jj = iSp+1, nType
        if (speciesName(iSp) == speciesName(jj)) then
          write(tmpStr,"('Duplicate identical species labels in the geometry: ',A)")speciesName(iSp)
          call error(tmpStr)
        end if
      end do
    end do

    ! Initialise the SCC module (the two copies of the Hubbard Us are rather
    ! artifical, since the copy for the main program is only used for dumping
    ! into the tagged format for autotest)
    allocate(hubbU(orb%mShell, nType))
    @:ASSERT(size(input%slako%skHubbU, dim=1) >= orb%mShell)
    @:ASSERT(size(input%slako%skHubbU, dim=2) == nType)
    hubbU(:,:) = input%slako%skHubbU(1:orb%mShell, :)
    if (allocated(input%ctrl%hubbU)) then
      where (input%ctrl%hubbU > 0.0_dp)
        hubbU = input%ctrl%hubbU
      end where
    end if
    if (tSccCalc) then
      allocate(sccInp)
      allocate(sccCalc)
      sccInp%orb => orb
      if (tPeriodic) then
        sccInp%latVecs = latVec
        sccInp%recVecs = recVec
        sccInp%volume = CellVol
      end if
      sccInp%hubbU = hubbU
      allocate(tDampedShort(nType))
      if (input%ctrl%tDampH) then
        tDampedShort = (speciesMass < 3.5_dp * amu__au)
        !tDampedShort(:) = (speciesName == "H" .or. speciesName == "h")
      else
        tDampedShort(:) = .false.
      end if
      sccInp%tDampedShort = tDampedShort
      sccInp%dampExp = input%ctrl%dampExp

      ! H5 correction
      if (input%ctrl%h5SwitchedOn) then
        if (.not. any(speciesMass < 3.5_dp * amu__au)) then
          call error("H5 correction used without H atoms present")
        end if
        if (any(tDampedShort)) then
          call error("H5 correction is not compatible with X-H damping")
        end if
        allocate(pH5Correction)
        call H5Corr_init(pH5Correction, speciesName, input%ctrl%h5RScale, input%ctrl%h5WScale,&
            & input%ctrl%h5ElementPara)
        sccInp%h5Correction = pH5Correction
      end if

      nExtChrg = input%ctrl%nExtChrg
      tExtChrg = (nExtChrg > 0)
      if (tExtChrg) then
        if (.not.tSccCalc) then
          call error("External charges can only be used in an SCC calculation")
        end if
        tStress = .false. ! Stress calculations not allowed
        @:ASSERT(size(input%ctrl%extChrg, dim=1) == 4)
        @:ASSERT(size(input%ctrl%extChrg, dim=2) == nExtChrg)
        sccInp%extCharges = input%ctrl%extChrg
        if (allocated(input%ctrl%extChrgBlurWidth)) then
          sccInp%blurWidths = input%ctrl%extChrgblurWidth
          if (any(sccInp%blurWidths < 0.0_dp)) then
            call error("Gaussian blur widths for charges may not be negative")
          end if
        end if
      end if
      if (allocated(input%ctrl%chrgConstr)) then
        @:ASSERT(all(shape(input%ctrl%chrgConstr) == (/ nAtom, 2 /)))
        if (any(abs(input%ctrl%chrgConstr(:,2)) > epsilon(1.0_dp))) then
          sccInp%chrgConstraints = input%ctrl%chrgConstr
        end if
      end if

      if (allocated(input%ctrl%thirdOrderOn)) then
        @:ASSERT(tSccCalc)
        @:ASSERT(all(shape(input%ctrl%thirdOrderOn) == (/ nAtom, 2 /)))
        sccInp%thirdOrderOn = input%ctrl%thirdOrderOn
      end if

      sccInp%ewaldAlpha = input%ctrl%ewaldAlpha
      sccInp%tolEwald = input%ctrl%tolEwald
      call initialize(sccCalc, env, sccInp)
      deallocate(sccInp)

      ! Longest cut-off including the softening part of gamma
      mCutOff = max(mCutOff, sccCalc%getCutOff())

      if (input%ctrl%t3rd .and. input%ctrl%tOrbResolved) then
        call error("Onsite third order DFTB only compatible with orbital non resolved SCC")
      end if

      ! Initialize full 3rd order module
      t3rd = input%ctrl%t3rd
      t3rdFull = input%ctrl%t3rdFull
      if (t3rdFull) then
        @:ASSERT(tSccCalc)
        thirdInp%orb => orb
        thirdInp%hubbUs = hubbU
        thirdInp%hubbUDerivs = input%ctrl%hubDerivs
        allocate(thirdInp%damped(nType))
        thirdInp%damped(:) = tDampedShort
        thirdInp%dampExp = input%ctrl%dampExp
        thirdInp%shellResolved = input%ctrl%tOrbResolved
        allocate(thirdOrd)
        call ThirdOrder_init(thirdOrd, thirdInp)
        mCutOff = max(mCutOff, thirdOrd%getCutOff())
      end if
    end if

    ! Initial coordinates
    allocate(coord0(3, nAtom))
    @:ASSERT(all(shape(coord0) == shape(input%geom%coords)))
    coord0(:,:) = input%geom%coords(:,:)
    tCoordsChanged = .true.

    allocate(species0(nAtom))
    @:ASSERT(all(shape(species0) == shape(input%geom%species)))
    species0(:) = input%geom%species(:)

    allocate(mass(nAtom))
    mass = speciesMass(species0)
    if (allocated(input%ctrl%masses)) then
      @:ASSERT(size(input%ctrl%masses) == nAtom)
      where (input%ctrl%masses >= 0.0_dp)
        mass = input%ctrl%masses
      end where
    end if

    if (tPeriodic) then
      ! Make some guess for the nr. of all interacting atoms
      nAllAtom = int((real(nAtom, dp)**(1.0_dp/3.0_dp) + 3.0_dp)**3)
    else
      nAllAtom = nAtom
    end if
    allocate(coord(3, nAllAtom))
    allocate(species(nAllAtom))
    allocate(img2CentCell(nAllAtom))
    allocate(iCellVec(nAllAtom))

    ! Intialize Hamilton and overlap
    tImHam = tDualSpinOrbit .or. (tSpinOrbit .and. tDFTBU) ! .or. tBField
    if (tSccCalc) then
      allocate(chargePerShell(orb%mShell,nAtom,nSpin))
    else
       allocate(chargePerShell(0,0,0))
    end if
    allocate(ham(0, nSpin))
    if (tImHam) then
      allocate(iHam(0, nSpin))
    end if
    allocate(over(0))
    allocate(iSparseStart(0, nAtom))

    if (nSpin == 4) then
      allocate(nEl(1))
    else
      allocate(nEl(nSpin))
    end if

    nEl0 = 0.0_dp
    do ii = 1, nAtom
      nEl0 = nEl0 + sum(input%slako%skOcc(1:orb%nShell(species0(ii)), species0(ii)))
    end do
    nEl(:) = 0.0_dp
    if (nSpin == 1 .or. nSpin == 4) then
      nEl(1) = nEl0 - input%ctrl%nrChrg
      if(ceiling(nEl(1)) > 2.0_dp*nOrb) then
        call error("More electrons than basis functions!")
      end if
    else
      nEl(1) = 0.5_dp * (nEl0 - input%ctrl%nrChrg + input%ctrl%nrSpinPol)
      nEl(2) = 0.5_dp * (nEl0 - input%ctrl%nrChrg - input%ctrl%nrSpinPol)
      if (any(ceiling(nEl(:)) > nOrb)) then
        call error("More electrons than basis functions!")
      end if
    end if

    if (.not.all(nEl(:) >= 0.0_dp)) then
      call error("Less than 0 electrons!")
    end if

    iDistribFn = input%ctrl%iDistribFn
    tempElec = input%ctrl%tempElec
    tFixEf = input%ctrl%tFixEf
    if (tFixEf) then
      Ef = input%ctrl%Ef
    else
      Ef = 0.0_dp
    end if
    tSetFillingTemp = input%ctrl%tSetFillingTemp
    tFillKSep = input%ctrl%tFillKSep
    tempAtom = input%ctrl%tempAtom
    deltaT = input%ctrl%deltaT


    ! Create equivalency relations
    if (tSccCalc) then
      allocate(iEqOrbitals(orb%mOrb, nAtom, nSpin))
      allocate(iEqOrbSCC(orb%mOrb, nAtom, nSpin))
      call sccCalc%getOrbitalEquiv(orb, species0, iEqOrbSCC)
      if (nSpin == 1) then
        iEqOrbitals(:,:,:) = iEqOrbSCC(:,:,:)
      else
        allocate(iEqOrbSpin(orb%mOrb, nAtom, nSpin))
        call Spin_getOrbitalEquiv(orb, species0, iEqOrbSpin)
        call OrbitalEquiv_merge(iEqOrbSCC, iEqOrbSpin, orb, iEqOrbitals)
        deallocate(iEqOrbSpin)
      end if
      deallocate(iEqOrbSCC)
      nIneqOrb = maxval(iEqOrbitals)
      nMixElements = nIneqOrb
      if (tDFTBU) then
        allocate(iEqOrbSpin(orb%mOrb, nAtom, nSpin))
        allocate(iEqOrbDFTBU(orb%mOrb, nAtom, nSpin))
        call DFTBplsU_getOrbitalEquiv(iEqOrbDFTBU,orb, species0, nUJ, niUJ, iUJ)
        call OrbitalEquiv_merge(iEqOrbitals, iEqOrbDFTBU, orb, iEqOrbSpin)
        iEqOrbitals(:,:,:) = iEqOrbSpin(:,:,:)
        nIneqOrb = maxval(iEqOrbitals)
        deallocate(iEqOrbSpin)
        deallocate(iEqOrbDFTBU)
        allocate(iEqBlockDFTBU(orb%mOrb, orb%mOrb, nAtom, nSpin))
        call DFTBU_blockIndx(iEqBlockDFTBU, nIneqOrb, orb, species0, nUJ, niUJ, iUJ)
        nMixElements = max(nMixElements,maxval(iEqBlockDFTBU)) ! as
        !  iEqBlockDFTBU does not include diagonal elements, so in the case of
        !  a purely s-block DFTB+U calculation, maxval(iEqBlockDFTBU) would
        !  return 0
        if (tImHam) then
          allocate(iEqBlockDFTBULS(orb%mOrb, orb%mOrb, nAtom, nSpin))
          call DFTBU_blockIndx(iEqBlockDFTBULS, nMixElements, orb, species0, nUJ, niUJ, iUJ)
          nMixElements = max(nMixElements,maxval(iEqBlockDFTBULS))
        end if
      end if
    else
      nIneqOrb = nOrb
      nMixElements = 0
    end if

    ! Initialize mixer
    ! (at the moment, the mixer does not need to know about the size of the
    ! vector to mix.)
    if (tSccCalc) then
      allocate(pChrgMixer)
      iMixer = input%ctrl%iMixSwitch
      nGeneration = input%ctrl%iGenerations
      mixParam = input%ctrl%almix
      select case (iMixer)
      case (1)
        allocate(pSimplemixer)
        call init(pSimpleMixer, mixParam)
        call init(pChrgMixer, pSimpleMixer)
      case (2)
        allocate(pAndersonMixer)
        if (input%ctrl%andersonNrDynMix > 0) then
          call init(pAndersonMixer, nGeneration, mixParam, input%ctrl%andersonInitMixing,&
              & input%ctrl%andersonDynMixParams, input%ctrl%andersonOmega0)
        else
          call init(pAndersonMixer, nGeneration, mixParam, input%ctrl%andersonInitMixing,&
              & omega0=input%ctrl%andersonOmega0)
        end if
        call init(pChrgMixer, pAndersonMixer)
      case (3)
        allocate(pBroydenMixer)
        call init(pBroydenMixer, maxSccIter, mixParam, input%ctrl%broydenOmega0,&
            & input%ctrl%broydenMinWeight, input%ctrl%broydenMaxWeight, input%ctrl%broydenWeightFac)
        call init(pChrgMixer, pBroydenMixer)
      case(4)
        allocate(pDIISMixer)
        call init(pDIISMixer,nGeneration, mixParam, input%ctrl%tFromStart)
        call init(pChrgMixer, pDIISMixer)
      case default
        call error("Unknown charge mixer type.")
      end select
    end if

    ! initialise in cases where atoms move
    tGeoOpt = input%ctrl%tGeoOpt
    tCoordOpt = input%ctrl%tCoordOpt
    tLatOpt = (input%ctrl%tLatOpt .and. tPeriodic)
    if (tLatOpt) then
      if (tExtChrg) then
        ! Stop as not sure, what to do with the coordinates of the
        ! external charges, when the lattice changes.
        call error("External charges and lattice optimisation can not be used together.")
      end if
    end if
    if (tLatOpt) then
      tLatOptFixAng = input%ctrl%tLatOptFixAng
      tLatOptFixLen = input%ctrl%tLatOptFixLen
      tLatOptIsotropic = input%ctrl%tLatOptIsotropic
      if (tLatOptFixAng .or. any(tLatOptFixLen) .or. tLatOptIsotropic) then
        origLatVec(:,:) = latVec(:,:)
        do ii = 1, 3
           normOrigLatVec(:,ii) = origLatVec(:,ii) / sqrt(sum(origLatVec(:,ii)**2))
        end do
      end if
    end if
    extPressure = input%ctrl%pressure
    tBarostat = input%ctrl%tBarostat
    BarostatStrength = input%ctrl%BarostatStrength

#:if WITH_SOCKETS
    tSocket = allocated(input%ctrl%socketInput)
    if (tSocket) then
      input%ctrl%socketInput%nAtom = nAtom
      call initSocket(env, input%ctrl%socketInput, tPeriodic, coord0, latVec, socket,&
          & tCoordsChanged, tLatticeChanged)
      tForces = .true.
      tGeoOpt = .false.
      tMD = .false.
    end if
#:else
    tSocket = .false.
#:endif

    tAppendGeo = input%ctrl%tAppendGeo
    tUseConvergedForces = (input%ctrl%tConvrgForces .and. tSccCalc) ! no point if not SCC
    tMD = input%ctrl%tMD
    tDerivs = input%ctrl%tDerivs
    tPrintMulliken = input%ctrl%tPrintMulliken
    tEField = input%ctrl%tEfield ! external electric field
    tMulliken = input%ctrl%tMulliken .or. tPrintMulliken .or. tEField
    tAtomicEnergy = input%ctrl%tAtomicEnergy
    tPrintEigVecs = input%ctrl%tPrintEigVecs
    tPrintEigVecsTxt = input%ctrl%tPrintEigVecsTxt

    tPrintForces = input%ctrl%tPrintForces
    tForces = input%ctrl%tForces .or. tPrintForces
    if (tSccCalc) then
      forceType = input%ctrl%forceType
    else
      if (input%ctrl%forceType /= 0) then
        call error("Invalid force evaluation method for non-SCC calculations.")
      end if
    end if
    if (forceType == 2 .and. tempElec > minTemp) then
       call error("This ForceEvaluation method requires the electron temperature to be zero")
    end if
    if (tForces) then
      select case(input%ctrl%iDerivMethod)
      case (1)
        ! set step size from input
        if (input%ctrl%deriv1stDelta < epsilon(1.0_dp)) then
          write(tmpStr, "(A,E12.4)") 'Too small value for finite difference step :',&
              & input%ctrl%deriv1stDelta
          call error(tmpStr)
        end if
        call NonSccDiff_init(nonSccDeriv, diffTypes%finiteDiff, input%ctrl%deriv1stDelta)
      case (2)
        call NonSccDiff_init(nonSccDeriv, diffTypes%richardson)
      end select
    end if

    ! requires stress to already be possible and it being a periodic calculation
    ! with forces
    tStress = ((tPeriodic .and. tForces).and.tStress)

    nMovedAtom = input%ctrl%nrMoved
    nMovedCoord = 3 * nMovedAtom

    if (input%ctrl%maxRun == -1) then
      nGeoSteps = huge(1) - 1
      ! Workaround:PGI 17.10 -> do i = 0, huge(1) executes 0 times
      ! nGeoSteps = huge(1)
    else
      nGeoSteps = input%ctrl%maxRun
    end if

    if (nMovedAtom > 0) then
      allocate(indMovedAtom(size(input%ctrl%indMovedAtom)))
      indMovedAtom(:) = input%ctrl%indMovedAtom(:)
    else
      allocate(indMovedAtom(0))
    end if

    allocate(pGeoCoordOpt)
    if (tCoordOpt) then
      allocate(tmpCoords(nMovedCoord))
      tmpCoords(1:nMovedCoord) = reshape(coord0(:, indMovedAtom), (/ nMovedCoord /))
      select case (input%ctrl%iGeoOpt)
      case(1)
        allocate(tmpWeight(nMovedCoord))
        tmpWeight(1:nMovedCoord) = 0.5_dp * deltaT**2 / reshape(spread(mass(indMovedAtom), 1, 3),&
            & (/nMovedCoord/))
        allocate(pSteepDesc)
        call init(pSteepDesc, size(tmpCoords), input%ctrl%maxForce,&
            & input%ctrl%maxAtomDisp,tmpWeight)
        deallocate(tmpWeight)
        call init(pGeoCoordOpt, pSteepDesc)
      case (2)
        allocate(pConjGrad)
        call init(pConjGrad, size(tmpCoords), input%ctrl%maxForce, input%ctrl%maxAtomDisp)
        call init(pGeoCoordOpt, pConjGrad)
      case (3)
        allocate(pDIIS)
        call init(pDIIS, size(tmpCoords), input%ctrl%maxForce, input%ctrl%deltaGeoOpt,&
            & input%ctrl%iGenGeoOpt)
        call init(pGeoCoordOpt, pDIIS)
      case (4)
        allocate(pLbfgs)
        call TLbfgs_init(pLbfgs, size(tmpCoords), input%ctrl%maxForce, tolSameDist,&
            & input%ctrl%maxAtomDisp, input%ctrl%lbfgsInp%memory)
        call init(pGeoCoordOpt, pLbfgs)
      end select
      call reset(pGeoCoordOpt, tmpCoords)
    end if

    allocate(pGeoLatOpt)
    if (tLatOpt) then
      select case (input%ctrl%iGeoOpt)
      case(1)
        allocate(tmpWeight(9))
        tmpWeight = 1.0_dp
        allocate(pSteepDescLat)
        call init(pSteepDescLat, 9, input%ctrl%maxForce, input%ctrl%maxLatDisp, tmpWeight )
        deallocate(tmpWeight)
        call init(pGeoLatOpt, pSteepDescLat)
      case(2,3) ! use CG lattice for both DIIS and CG
        allocate(pConjGradLat)
        call init(pConjGradLat, 9, input%ctrl%maxForce, input%ctrl%maxLatDisp)
        call init(pGeoLatOpt, pConjGradLat)
      case (4)
        allocate(pLbfgsLat)
        call TLbfgs_init(pLbfgsLat, 9, input%ctrl%maxForce, tolSameDist, input%ctrl%maxLatDisp,&
            & input%ctrl%lbfgsInp%memory)
        call init(pGeoLatOpt, pLbfgsLat)
      end select
      ! optimization uses scaling factor of unit cell
      if (tLatOptIsotropic ) then
        call reset(pGeoLatOpt, (/1.0_dp,0.0_dp,0.0_dp,0.0_dp,0.0_dp,0.0_dp,0.0_dp,0.0_dp,0.0_dp/))
      else if (tLatOptFixAng) then
        ! optimization uses scaling factor of lattice vectors
        call reset( pGeoLatOpt, (/1.0_dp,1.0_dp,1.0_dp,0.0_dp,0.0_dp,0.0_dp,0.0_dp,0.0_dp,0.0_dp/))
      else
        call reset( pGeoLatOpt, reshape(latVec, (/ 9 /)) )
      end if
    end if

    if (.not.(tGeoOpt.or.tMD.or.tSocket)) then
      nGeoSteps = 0
    end if

    ! Initialize constraints
    if (input%ctrl%nrConstr > 0) then
      allocate(conAtom(input%ctrl%nrConstr))
      allocate(conVec(3, input%ctrl%nrConstr))
      conAtom(:) = input%ctrl%conAtom
      conVec(:,:) = input%ctrl%conVec
      do ii = 1, input%ctrl%nrConstr
        conVec(:,ii) = conVec(:,ii) / sqrt(sum(conVec(:,ii)**2))
      end do
    end if

    ! Dispersion
    tHHRepulsion = .false.
    tDispersion = allocated(input%ctrl%dispInp)
    if (tDispersion) then
      if (allocated(input%ctrl%dispInp%slakirk)) then
        tStress = .false.
        if (tLatOpt) then
          call error("Sorry, lattice optimisation and Slater-Kirkwood type dispersion can not be&
              & used together")
        end if
        if (tBarostat) then
          call error("Sorry, barostatic MD and Slater-Kirkwood type dispersion can not be used&
              & together")
        end if
        allocate(slaKirk)
        if (tPeriodic) then
          call DispSlaKirk_init(slaKirk, input%ctrl%dispInp%slakirk, latVec)
        else
          call DispSlaKirk_init(slaKirk, input%ctrl%dispInp%slakirk)
        end if
        call move_alloc(slaKirk, dispersion)

      elseif (allocated(input%ctrl%dispInp%uff)) then
        allocate(uff)
        if (tPeriodic) then
          call DispUff_init(uff, input%ctrl%dispInp%uff, nAtom, species0, latVec)
        else
          call DispUff_init(uff, input%ctrl%dispInp%uff, nAtom)
        end if
        call move_alloc(uff, dispersion)

    #:if WITH_DFTD3
      elseif (allocated(input%ctrl%dispInp%dftd3)) then
        allocate(dftd3)
        tHHRepulsion = input%ctrl%dispInp%dftd3%hhrepulsion
        if (tHHRepulsion .and. .not. any(speciesMass < 3.5_dp * amu__au)) then
          call error("H-H repulsion correction used without H atoms present")
        end if
        if (tPeriodic) then
          call DispDftD3_init(dftd3, input%ctrl%dispInp%dftd3, nAtom, species0, speciesName, latVec)
        else
          call DispDftD3_init(dftd3, input%ctrl%dispInp%dftd3, nAtom, species0, speciesName)
        end if
        call move_alloc(dftd3, dispersion)
    #:endif
      end if
      mCutOff = max(mCutOff, dispersion%getRCutOff())

    end if

    if (input%ctrl%nrChrg == 0.0_dp .and. (.not.tPeriodic) .and. tMulliken) then
      tDipole = .true.
    else
      tDipole = .false.
    end if

    if (allocated(input%ctrl%elStatPotentialsInp)) then
      if (.not.tSccCalc) then
        call error("Electrostatic potentials only available for SCC calculations")
      end if
      allocate(esp)
      call TElStatPotentials_init(esp, input%ctrl%elStatPotentialsInp, tEField .or. tExtChrg)
    end if

    tLocalise = input%ctrl%tLocalise
    if (tLocalise .and. (nSpin > 2 .or. t2Component)) then
      call error("Localisation of electronic states currently unsupported for non-collinear and&
          & spin orbit calculations")
    end if
    if (allocated(input%ctrl%pipekMezeyInp)) then
      allocate(pipekMezey)
      call initialise(pipekMezey, input%ctrl%pipekMezeyInp)
    end if
    tLocalise = allocated(pipekMezey)

    ! Linear response
    tLinResp = input%ctrl%lrespini%tInit

    if (tLinResp) then

      ! input sanity checking
    #:if not WITH_ARPACK
      call error("This binary has been compiled without support for linear response calculations.")
    #:endif
      if (.not. tSccCalc) then
        call error("Linear response excitation requires SCC=Yes")
      end if
      if (nspin > 2) then
        call error("Linear reponse does not work with non-colinear spin polarization yet")
      elseif (tSpin .and. tForces) then
        call error("excited state relaxation is not implemented yet for spin-polarized systems")
      elseif (tPeriodic .and. tForces) then
        call error("excited state relaxation is not implemented yet periodic systems")
      elseif (tPeriodic .and. .not.tRealHS) then
        call error("Linear response only works with non-periodic or gamma-point molecular crystals")
      elseif (tSpinOrbit) then
        call error("Linear response does not support spin orbit coupling at the moment.")
      elseif (tDFTBU) then
        call error("Linear response does not support LDA+U yet")
      elseif (input%ctrl%tOrbResolved) then
        call error("Linear response does not support orbital resolved scc yet")
      end if
      if (tempElec > 0.0_dp .and. tForces) then
        write(tmpStr, "(A,E12.4,A)")"Excited state forces are not implemented yet for fractional&
            & occupations, kT=", tempElec/Boltzmann,"K"
        call warning(tmpStr)
      end if

      if (input%ctrl%lrespini%nstat == 0) then
        if (tForces) then
          call error("Excited forces only available for StateOfInterest non zero.")
        end if
        if (input%ctrl%lrespini%tPrintEigVecs .or. input%ctrl%lrespini%tCoeffs) then
          call error("Excited eigenvectors only available for StateOfInterest non zero.")
        end if
      end if
      if (input%ctrl%lrespini%energyWindow < 0.0_dp) then
        call error("Negative energy window for excitations")
      end if

      ! Hubbard U and spin constants for excitations (W only needed for triplet/spin polarised)
      allocate(input%ctrl%lrespini%HubbardU(nType))
      allocate(input%ctrl%lrespini%spinW(nType))
      input%ctrl%lrespini%HubbardU = 0.0_dp
      input%ctrl%lrespini%spinW = 0.0_dp

      ! calculate linear response Gamma values from HOAO shell Hubbard U (non
      ! shell resolved)
      do iSp = 1, nType
        homoLoc = maxloc(atomEigVal(:orb%nShell(iSp), iSp),&
            & mask=input%slako%skOcc(:orb%nShell(iSp), iSp) > 0.0_dp)
        input%ctrl%lrespini%HubbardU(iSp) = hubbU(homoLoc(1), iSp)
      end do

      ! and atomic HOAO spin W value if needed
      input%ctrl%lrespini%spinW(:) = 0.0_dp
      select case(input%ctrl%lrespini%sym)
      case("S")
        ! Singlet case, no need for spin constants
      case("T","B"," ")
        ! triplet or spin-polarised
        do iSp = 1, nType
          homoLoc = maxloc(atomEigVal(:orb%nShell(iSp), iSp),&
              & mask=input%slako%skOcc(:orb%nShell(iSp), iSp) > 0.0_dp)
          input%ctrl%lrespini%spinW(iSp) = spinW(homoLoc(1), homoLoc(1), iSp)
        end do
      case default
        call error("Unknown excitation type requested")
      end select

      tPrintExcitedEigVecs = input%ctrl%lrespini%tPrintEigVecs
      tLinRespZVect = (input%ctrl%lrespini%tMulliken .or. tForces .or. input%ctrl%lrespini%tCoeffs&
          & .or. tPrintExcitedEigVecs)

      call init(lresp, input%ctrl%lrespini, nAtom, nEl(1), orb)

    end if

    iSeed = input%ctrl%iSeed
    tRandomSeed = (iSeed < 1)
    ! Note: This routine may not be called multiple times. If you need further random generators,
    ! extend the routine and create them within this call.
    call createRandomGenerators(env, iSeed, randomInit, randomThermostat)

    call getRandom(randomInit, rTmp)
    runId = int(real(huge(runId) - 1, dp) * rTmp) + 1

    ! MD stuff
    if (tMD) then
      ! Create MD framework.
      allocate(pMDFrame)
      call init(pMDFrame, nMovedAtom, nAtom, input%ctrl%tMDstill)

      ! Create temperature profile, if thermostat is not the dummy one
      if (input%ctrl%iThermostat /= 0) then
        allocate(temperatureProfile)
        call init(temperatureProfile, input%ctrl%tempMethods, input%ctrl%tempSteps,&
            & input%ctrl%tempValues)
        pTempProfile => temperatureProfile
      else
        nullify(pTempProfile)
      end if

      ! Create thermostat
      allocate(pThermostat)
      select case (input%ctrl%iThermostat)
      case (0) ! No thermostat
        allocate(pDummyTherm)
        call init(pDummyTherm, tempAtom, mass(indMovedAtom), randomThermostat, pMDFrame)
        call init(pThermostat, pDummyTherm)
      case (1) ! Andersen thermostat
        allocate(pAndersenTherm)
        call init(pAndersenTherm, randomThermostat, mass(indMovedAtom), pTempProfile,&
            & input%ctrl%tRescale, input%ctrl%wvScale, pMDFrame)
        call init(pThermostat, pAndersenTherm)
      case (2) ! Berendsen thermostat
        allocate(pBerendsenTherm)
        call init(pBerendsenTherm, randomThermostat, mass(indMovedAtom), pTempProfile,&
            & input%ctrl%wvScale, pMDFrame)
        call init(pThermostat, pBerendsenTherm)
      case (3) ! Nose-Hoover-Chain thermostat
        allocate(pNHCTherm)
        if (input%ctrl%tInitNHC) then
          call init(pNHCTherm, randomThermostat, mass(indMovedAtom), pTempProfile,&
              & input%ctrl%wvScale, pMDFrame, input%ctrl%deltaT, input%ctrl%nh_npart,&
              & input%ctrl%nh_nys, input%ctrl%nh_nc, input%ctrl%xnose, input%ctrl%vnose,&
              & input%ctrl%gnose)
        else
          call init(pNHCTherm, randomThermostat, mass(indMovedAtom), pTempProfile,&
              & input%ctrl%wvScale, pMDFrame, input%ctrl%deltaT, input%ctrl%nh_npart,&
              & input%ctrl%nh_nys, input%ctrl%nh_nc)
        end if
        call init(pThermostat, pNHCTherm)
      end select

      ! Create MD integrator
      allocate(pVelocityVerlet)
      if (input%ctrl%tReadMDVelocities) then
        if (tBarostat) then
          call init(pVelocityVerlet, deltaT, coord0(:,indMovedAtom), pThermostat,&
              & input%ctrl%initialVelocities, BarostatStrength, extPressure, input%ctrl%tIsotropic)
        else
          call init(pVelocityVerlet, deltaT, coord0(:,indMovedAtom), pThermostat,&
              & input%ctrl%initialVelocities)
        end if
      else
        if (tBarostat) then
          call init(pVelocityVerlet, deltaT, coord0(:,indMovedAtom), pThermostat, BarostatStrength,&
              & extPressure, input%ctrl%tIsotropic)
        else
          call init(pVelocityVerlet, deltaT, coord0(:,indMovedAtom), pThermostat)
        end if
      end if
      allocate(pMDIntegrator)
      call init(pMDIntegrator, pVelocityVerlet)
    end if

    ! Check for extended Born-Oppenheimer MD
    tXlbomd = allocated(input%ctrl%xlbomd)
    if (tXlbomd) then
      if (input%ctrl%iThermostat /= 0) then
        call error("XLBOMD does not work with thermostats yet")
      elseif (tBarostat) then
        call error("XLBOMD does not work with barostats yet")
      elseif (nSpin /= 1 .or. tDFTBU) then
        call error("XLBOMD does not work for spin or DFTB+U yet")
      elseif (forceType /= 2 .and. forceType /= 3) then
        call error("Force evaluation method incompatible with XLBOMD")
      elseif (iDistribFn /= Fermi) then
        call error("Filling function incompatible with XLBOMD")
      end if
      allocate(xlbomdIntegrator)
      call Xlbomd_init(xlbomdIntegrator, input%ctrl%xlbomd, nIneqOrb)
    end if

    if (tDerivs) then
      allocate(tmp3Coords(3,nMovedAtom))
      tmp3Coords = coord0(:,indMovedAtom)
      call create(derivDriver, tmp3Coords, input%ctrl%deriv2ndDelta)
      coord0(:,indMovedAtom) = tmp3Coords
      deallocate(tmp3Coords)
      nGeoSteps = 2 * 3 * nMovedAtom - 1
    end if

    if (tEField) then
      EFieldStrength = input%ctrl%EFieldStrength
      EfieldVector(:) = input%ctrl%EfieldVector(:)
      tTDEfield = input%ctrl%tTDEfield
      EfieldOmega = input%ctrl%EfieldOmega
      EfieldPhase = input%ctrl%EfieldPhase
      if (tTDEfield .and. .not. tMD) then
        call error ("Time dependent electric fields only possible for MD!")
      end if
      ! parser should catch all of these:
      @:ASSERT(.not.tTDEfield .or. tMD)
    else
      tEField = .false.
      EFieldStrength = 0.0_dp
      EfieldVector(:) = 0.0_dp
      tTDEfield = .false.
      EfieldOmega = 0.0_dp
      EfieldPhase = 0
    end if

    ! Allocate charge arrays
    if (tMulliken) then ! automatically true if tSccCalc
      allocate(q0(orb%mOrb, nAtom, nSpin))
      q0(:,:,:) = 0.0_dp

      allocate(qShell0(orb%mShell, nAtom))
      qShell0(:,:) = 0.0_dp
    else
      allocate(q0(0,0,0))
      allocate(qShell0(0,0))
    end if

    allocate(qInput(orb%mOrb, nAtom, nSpin))
    allocate(qOutput(orb%mOrb, nAtom, nSpin))
    qInput(:,:,:) = 0.0_dp
    qOutput(:,:,:) = 0.0_dp

    if (tDFTBU) then
      allocate(qBlockIn(orb%mOrb, orb%mOrb, nAtom, nSpin))
      allocate(qBlockOut(orb%mOrb, orb%mOrb, nAtom, nSpin))
      qBlockIn(:,:,:,:) = 0.0_dp
      qBlockOut(:,:,:,:) = 0.0_dp
      if (tImHam) then
        allocate(qiBlockIn(orb%mOrb, orb%mOrb, nAtom, nSpin))
        qiBlockIn(:,:,:,:) = 0.0_dp
      end if
    end if

    if (tImHam) then
      allocate(qiBlockOut(orb%mOrb, orb%mOrb, nAtom, nSpin))
      qiBlockOut(:,:,:,:) = 0.0_dp
    end if

    if (tSccCalc) then
      allocate(qDiffRed(nMixElements))
      allocate(qInpRed(nMixElements))
      allocate(qOutRed(nMixElements))
      qDiffRed = 0.0_dp
      qInpRed = 0.0_dp
      qOutRed = 0.0_dp
    end if

    ! Initialize Mulliken charges
    if (tMulliken .or. tLinResp) then
      call initQFromShellChrg(q0, referenceN0, species0, orb)
    end if

    tReadChrg = input%ctrl%tReadChrg
    tWriteChrgAscii = input%ctrl%tWriteChrgAscii
    if (tSccCalc) then
      do iAt = 1, nAtom
        iSp = species0(iAt)
        do iSh = 1, orb%nShell(iSp)
          qShell0 (iSh,iAt) = sum(q0(orb%posShell(iSh,iSp):orb%posShell(iSh+1,iSp)-1,iAt,1))
        end do
      end do
      if (tReadChrg) then
        if (tDFTBU) then
          if (nSpin == 2) then
            if (tFixEf .or. input%ctrl%tSkipChrgChecksum) then
              ! do not check charge or magnetisation from file
              call initQFromFile(qInput, fCharges, input%ctrl%tReadChrgAscii, orb, qBlock=qBlockIn)
            else
              call initQFromFile(qInput, fCharges, input%ctrl%tReadChrgAscii, orb, nEl = sum(nEl),&
                  & magnetisation=nEl(1)-nEl(2), qBlock=qBlockIn)
            end if
          else
            if (tImHam) then
              if (tFixEf .or. input%ctrl%tSkipChrgChecksum) then
                ! do not check charge or magnetisation from file
                call initQFromFile(qInput, fCharges, input%ctrl%tReadChrgAscii, orb,&
                    & qBlock=qBlockIn,qiBlock=qiBlockIn)
              else
                call initQFromFile(qInput, fCharges, input%ctrl%tReadChrgAscii, orb, nEl = nEl(1),&
                    & qBlock=qBlockIn,qiBlock=qiBlockIn)
              end if
            else
              if (tFixEf .or. input%ctrl%tSkipChrgChecksum) then
                ! do not check charge or magnetisation from file
                call initQFromFile(qInput, fCharges, input%ctrl%tReadChrgAscii, orb,&
                    & qBlock=qBlockIn)
              else
                call initQFromFile(qInput, fCharges, input%ctrl%tReadChrgAscii, orb, nEl = nEl(1),&
                    & qBlock=qBlockIn)
              end if
            end if
          end if
        else
          ! hack again caused by going from up/down to q and M
          if (nSpin == 2) then
            if (tFixEf .or. input%ctrl%tSkipChrgChecksum) then
              ! do not check charge or magnetisation from file
              call initQFromFile(qInput, fCharges, input%ctrl%tReadChrgAscii, orb)
            else
              call initQFromFile(qInput, fCharges, input%ctrl%tReadChrgAscii, orb, nEl = sum(nEl),&
                  & magnetisation=nEl(1)-nEl(2))
            end if
          else
            if (tFixEf .or. input%ctrl%tSkipChrgChecksum) then
              ! do not check charge or magnetisation from file
              call initQFromFile(qInput, fCharges, input%ctrl%tReadChrgAscii, orb)
            else
              call initQFromFile(qInput, fCharges, input%ctrl%tReadChrgAscii, orb, nEl = nEl(1))
            end if
          end if
        end if
      else
        if (allocated(input%ctrl%initialCharges)) then
          if (abs(sum(input%ctrl%initialCharges) - input%ctrl%nrChrg) > 1e-4_dp) then
            write(strTmp, "(A,G13.6,A,G13.6,A,A)") "Sum of initial charges does not match specified&
                & total charge. (", sum(input%ctrl%initialCharges), " vs. ", input%ctrl%nrChrg,&
                & ") ", "Your initial charge distribution will be rescaled."
            call warning(strTmp)
          end if
          call initQFromAtomChrg(qInput, input%ctrl%initialCharges, referenceN0, species0,&
              & speciesName, orb)
        else
          qInput(:,:,:) = q0
        end if
        if (.not. input%ctrl%tSkipChrgChecksum) then
          ! Rescaling to ensure correct number of electrons in the system
          qInput(:,:,1) = qInput(:,:,1) *  sum(nEl) / sum(qInput(:,:,1))
        end if

        select case (nSpin)
        case (1)
          ! nothing to do
        case (2)
          if (allocated(input%ctrl%initialSpins)) then
            do ii = 1, nAtom
              ! does not actually matter if additional spin polarization pushes
              ! charges to <0 as the initial charges are not mixed in to later
              ! iterations
              qInput(1:orb%nOrbAtom(ii),ii,2) = qInput(1:orb%nOrbAtom(ii),ii,1)&
                  & * input%ctrl%initialSpins(1,ii) / sum(qInput(1:orb%nOrbAtom(ii),ii,1))
            end do
          else
            if (.not. input%ctrl%tSkipChrgChecksum) then
              do ii = 1, nAtom
                qInput(1:orb%nOrbAtom(ii),ii,2) = qInput(1:orb%nOrbAtom(ii),ii,1)&
                    & * (nEl(1)-nEl(2))/sum(qInput(:,:,1))
              end do
            end if
          end if
        case (4)
          if (tSpin) then
            if (.not. allocated(input%ctrl%initialSpins)) then
              call error("Missing initial spins!")
            end if
            if (any(shape(input%ctrl%initialSpins)/=(/3,nAtom/))) then
              call error("Incorrect shape initialSpins array!")
            end if
            if (.not. input%ctrl%tSkipChrgChecksum) then
              do ii = 1, nAtom
                do jj = 1, 3
                  qInput(1:orb%nOrbAtom(ii),ii,jj+1) = qInput(1:orb%nOrbAtom(ii),ii,1)&
                      & * input%ctrl%initialSpins(jj,ii) / sum(qInput(1:orb%nOrbAtom(ii),ii,1))
                end do
              end do
            end if
          end if
        end select
        if (tDFTBU) then
          qBlockIn = 0.0_dp
          do iS = 1, nSpin
            do iAt = 1, nAtom
              iSp = species0(iAt)
              do iSh = 1, orb%nShell(iSp)
                iStart = orb%posShell(iSh,iSp)
                iEnd = orb%posShell(iSh+1,iSp)-1
                rTmp = sum(qInput(iStart:iEnd,iAt,iS))
                rTmp = rTmp / real(iEnd+1-iStart,dp)
                do ii = iStart, iEnd
                  qBlockIn(ii,ii,iAt,iS) = rTmp
                end do
              end do
            end do
          end do
          if (tImHam) then
            qiBlockIn = 0.0_dp
          end if
        end if
      end if

      qInpRed = 0.0_dp
      if (nSpin == 2) then
        call qm2ud(qInput)
        if (tDFTBU) then
          call qm2ud(qBlockIn)
        end if
      end if

      call OrbitalEquiv_reduce(qInput, iEqOrbitals, orb, qInpRed(1:nIneqOrb))
      if (tDFTBU) then
        call AppendBlock_reduce(qBlockIn, iEqBlockDFTBU, orb, qInpRed )
        if (tImHam) then
          call AppendBlock_reduce(qiBlockIn, iEqBlockDFTBULS, orb, qInpRed, skew=.true. )
        end if
      end if

      if (nSpin == 2) then
        call ud2qm(qInput)
        if (tDFTBU) then
          call ud2qm(qBlockIn)
        end if
      end if
    end if

    ! Initialise images (translations)
    if (tPeriodic) then
      call getCellTranslations(cellVec, rCellVec, latVec, invLatVec, mCutOff)
    else
      allocate(cellVec(3, 1))
      allocate(rCellVec(3, 1))
      cellVec(:,1) = [0.0_dp, 0.0_dp, 0.0_dp]
      rCellVec(:,1) = [0.0_dp, 0.0_dp, 0.0_dp]
    end if

    ! Initialize neighbourlist.
    allocate(neighbourList)
    call init(neighbourList, nAtom, nInitNeighbour)
    allocate(nNeighbourSK(nAtom))
    allocate(nNeighbourRep(nAtom))

    ! Set various options
    tWriteAutotest = env%tGlobalMaster .and. input%ctrl%tWriteTagged
    tWriteDetailedXML = env%tGlobalMaster .and. input%ctrl%tWriteDetailedXML
    tWriteResultsTag = env%tGlobalMaster .and. input%ctrl%tWriteResultsTag
    tWriteDetailedOut = env%tGlobalMaster .and. input%ctrl%tWriteDetailedOut
    tWriteBandDat = env%tGlobalMaster .and. input%ctrl%tWriteBandDat
    tWriteHS = input%ctrl%tWriteHS
    tWriteRealHS = input%ctrl%tWriteRealHS

    ! Check if stopfiles already exist and quit if yes
    inquire(file=fStopSCC, exist=tExist)
    if (tExist) then
      call error("Stop file '" // fStopSCC // "' already present at startup")
    end if
    inquire(file=fStopDriver, exist=tExist)
    if (tExist) then
      call error("Stop file '" // fStopDriver // "' already present at startup")
    end if

    restartFreq = input%ctrl%restartFreq

    if (env%tGlobalMaster) then
      call initOutputFiles(env, tWriteAutotest, tWriteResultsTag, tWriteBandDat, tDerivs,&
          & tWriteDetailedOut, tMd, tGeoOpt, geoOutFile, fdDetailedOut, fdMd, esp)
    end if

    call getDenseDescCommon(orb, nAtom, t2Component, denseDesc)

  #:if WITH_SCALAPACK
    associate (blacsOpts => input%ctrl%parallelOpts%blacsOpts)
      call getDenseDescBlacs(env, blacsOpts%blockSize, blacsOpts%blockSize, denseDesc)
    end associate
  #:endif

    call initArrays(env, tForces, tExtChrg, tLinResp, tLinRespZVect, tMd, tMulliken, tSpinOrbit,&
        & tImHam, tWriteRealHS, tWriteHS, t2Component, tRealHS, tPrintExcitedEigvecs, tDipole, orb,&
        & nAtom, nMovedAtom, nKPoint, nSpin, nExtChrg, indMovedAtom, mass, denseDesc, rhoPrim, h0,&
        & iRhoPrim, excitedDerivs, ERhoPrim, derivs, chrgForces, energy, potential, TS, E0, Eband,&
        & eigen, filling, coord0Fold, newCoords, orbitalL, HSqrCplx, SSqrCplx, eigvecsCplx,&
        & HSqrReal, SSqrReal, eigvecsReal, rhoSqrReal, chargePerShell, occNatural, velocities,&
        & movedVelo, movedAccel, movedMass, dipoleMoment)

    if (tShowFoldedCoord) then
      pCoord0Out => coord0Fold
    else
      pCoord0Out => coord0
    end if


    ! Projection of eigenstates onto specific regions of the system
    tProjEigenvecs = input%ctrl%tProjEigenvecs
    if (tProjEigenvecs) then
      call init(iOrbRegion)
      call init(regionLabels)
      do iReg = 1, size(input%ctrl%tShellResInRegion)
        call elemShape(input%ctrl%iAtInRegion, valshape, iReg)
        nAtomRegion = valshape(1)
        allocate(iAtomRegion(nAtomRegion))
        call intoArray(input%ctrl%iAtInRegion, iAtomRegion, iTmp, iReg)
        if (input%ctrl%tOrbResInRegion(iReg) .or. input%ctrl%tShellResInRegion(iReg)) then

          if (input%ctrl%tOrbResInRegion(iReg)) then
            iSp = species0(iAtomRegion(1)) ! all atoms the same in the region
            @:ASSERT(all(species0(iAtomRegion) == iSp))
            nOrbRegion = nAtomRegion
            ! Create orbital index.
            allocate(tmpir1(nOrbRegion))
            do iOrb = 1, orb%nOrbSpecies(iSp)
              tmpir1 = 0
              ind = 1
              do iAt = 1, nAtomRegion
                tmpir1(ind) = denseDesc%iAtomStart(iAt) + iOrb - 1
                ind = ind + 1
              end do
              call append(iOrbRegion, tmpir1)
              write(tmpStr, "(A,'.',I0,'.',I0,'.out')")trim(input%ctrl%RegionLabel(iReg)),&
                  & orb%iShellOrb(iOrb,iSp), iOrb-orb%posShell(orb%iShellOrb(iOrb,iSp),iSp)&
                  & -orb%angShell(orb%iShellOrb(iOrb,iSp),iSp)
              call append(regionLabels, tmpStr)
            end do
            deallocate(tmpir1)
          end if

          if (input%ctrl%tShellResInRegion(iReg)) then
            iSp = species0(iAtomRegion(1)) ! all atoms the same in the region
            @:ASSERT(all(species0(iAtomRegion) == iSp))
            ! Create a separate region for each shell. It will contain
            ! the orbitals of that given shell for each atom in the region.
            do iSh = 1, orb%nShell(iSp)
              nOrbRegion = nAtomRegion * (orb%posShell(iSh + 1, iSp) - orb%posShell(iSh, iSp))
              ind = 1
              ! Create orbital index.
              allocate(tmpir1(nOrbRegion))
              do ii = 1, nAtomRegion
                iAt = iAtomRegion(ii)
                do jj = orb%posShell(iSh, iSp), orb%posShell(iSh + 1, iSp) - 1
                  tmpir1(ind) = denseDesc%iAtomStart(iAt) + jj - 1
                  ind = ind + 1
                end do
              end do
              call append(iOrbRegion, tmpir1)
              deallocate(tmpir1)
              write(tmpStr, "(A,'.',I0,'.out')")trim(input%ctrl%RegionLabel(iReg)), iSh
              call append(regionLabels, tmpStr)
            end do
          end if

        else
          ! We take all orbitals from all atoms.
          nOrbRegion = 0
          do ii = 1, nAtomRegion
            nOrbRegion = nOrbRegion + orb%nOrbAtom(iAtomRegion(ii))
          end do
          ind = 1
          allocate(tmpir1(nOrbRegion))
          ! Create an index of the orbitals
          do ii = 1, nAtomRegion
            iAt = iAtomRegion(ii)
            do jj = 1, orb%nOrbAtom(iAt)
              tmpir1(ind) = denseDesc%iAtomStart(iAt) + jj - 1
              ind = ind + 1
            end do
          end do
          call append(iOrbRegion, tmpir1)
          deallocate(tmpir1)
          write(tmpStr, "(A,'.out')") trim(input%ctrl%RegionLabel(iReg))
          call append(regionLabels, tmpStr)
        end if
        deallocate(iAtomRegion)
      end do
    end if

  #:if WITH_MPI
    if (env%mpi%nGroup > 1) then
      write(stdOut, "('MPI processes: ',T30,I0,' (split into ',I0,' groups)')")&
          & env%mpi%globalComm%size, env%mpi%nGroup
    else
      write(stdOut, "('MPI processes:',T30,I0)") env%mpi%globalComm%size
    end if
  #:endif

    write(stdOut, "('OpenMP threads: ', T30, I0)") omp_get_max_threads()

  #:if WITH_MPI
    if (omp_get_max_threads() > 1 .and. .not. input%ctrl%parallelOpts%tOmpThreads) then
      write(stdOut, *)
      call error("You must explicitely enable OpenMP threads (UseOmpThreads = Yes) if you wish to&
          & run an MPI-parallelised binary with OpenMP threads. If not, make sure that the&
          & environment variable OMP_NUM_THREADS is set to 1.")
      write(stdOut, *)
    end if
  #:endif

  #:if WITH_SCALAPACK
    write(stdOut, "('BLACS orbital grid size:', T30, I0, ' x ', I0)")env%blacs%orbitalGrid%nRow,&
        & env%blacs%orbitalGrid%nCol
    write(stdOut, "('BLACS atom grid size:', T30, I0, ' x ', I0)")env%blacs%atomGrid%nRow,&
        & env%blacs%atomGrid%nCol
  #:endif

    if (tRandomSeed) then
      write(stdOut, "(A,':',T30,I14)") "Chosen random seed", iSeed
    else
      write(stdOut, "(A,':',T30,I14)") "Specified random seed", iSeed
    end if

    if (input%ctrl%tMD) then
      select case(input%ctrl%iThermostat)
      case (0)
        if (tBarostat) then
          write(stdOut, "('Mode:',T30,A,/,T30,A)") 'MD without scaling of velocities',&
              & '(a.k.a. "NPE" ensemble)'
        else
          write(stdOut, "('Mode:',T30,A,/,T30,A)") 'MD without scaling of velocities',&
              & '(a.k.a. NVE ensemble)'
        end if
      case (1)
        if (tBarostat) then
          write(stdOut, "('Mode:',T30,A,/,T30,A)")&
              & "MD with re-selection of velocities according to temperature",&
              & "(a.k.a. NPT ensemble using Andersen thermostating + Berensen barostat)"
        else
          write(stdOut, "('Mode:',T30,A,/,T30,A)")&
              & "MD with re-selection of velocities according to temperature",&
              & "(a.k.a. NVT ensemble using Andersen thermostating)"
        end if
      case(2)
        if (tBarostat) then
          write(stdOut, "('Mode:',T30,A,/,T30,A)")&
              & "MD with scaling of velocities according to temperature",&
              & "(a.k.a. 'not' NVP ensemble using Berendsen thermostating and barostat)"
        else
          write(stdOut, "('Mode:',T30,A,/,T30,A)")&
              & "MD with scaling of velocities according to temperature",&
              & "(a.k.a. 'not' NVT ensemble using Berendsen thermostating)"
        end if
      case(3)
        if (tBarostat) then
          write(stdOut, "('Mode:',T30,A,/,T30,A)")"MD with scaling of velocities according to",&
              & "Nose-Hoover-Chain thermostat + Berensen barostat"
        else
          write(stdOut, "('Mode:',T30,A,/,T30,A)")"MD with scaling of velocities according to",&
              & "Nose-Hoover-Chain thermostat"
        end if

      case default
        call error("Unknown thermostat mode")
      end select
    elseif (tGeoOpt) then
      if (allocated(conAtom)) then
        strTmp = "with constraints"
      else
        strTmp = ""
      end if
      select case (input%ctrl%iGeoOpt)
      case (1)
        write(stdOut, "('Mode:',T30,A)")'Steepest descent' // trim(strTmp)
      case (2)
        write(stdOut, "('Mode:',T30,A)") 'Conjugate gradient relaxation' // trim(strTmp)
      case (3)
        write(stdOut, "('Mode:',T30,A)") 'Modified gDIIS relaxation' // trim(strTmp)
      case (4)
        write(stdout, "('Mode:',T30,A)") 'LBFGS relaxation' // trim(strTmp)
      case default
        call error("Unknown optimisation mode")
      end select
    elseif (tDerivs) then
      write(stdOut, "('Mode:',T30,A)") "2nd derivatives calculation"
      write(stdOut, "('Mode:',T30,A)") "Calculated for atoms:"
      write(stdOut, *) indMovedAtom
    elseif (tSocket) then
      write(stdOut, "('Mode:',T30,A)") "Socket controlled calculation"
    else
      write(stdOut, "('Mode:',T30,A)") "Static calculation"
    end if

    if (tSccCalc) then
      write(stdOut, "(A,':',T30,A)")"Self consistent charges", "Yes"
      write(stdOut, "(A,':',T30,E14.6)")"SCC-tolerance", sccTol
      write(stdOut, "(A,':',T30,I14)")"Max. scc iterations", maxSccIter
      if (tDFTBU) then
        write(stdOut, "(A,':',T35,A)")"Orbitally dependant functional", "Yes"
        write(stdOut, "(A,':',T30,A)")"Orbital functional", trim(plusUFunctionals%names(nDFTBUfunc))
      end if
    else
      write(stdOut, "(A,':',T30,A)") "Self consistent charges", "No"
    end if

    select case (nSpin)
    case(1)
      write(stdOut, "(A,':',T30,A)") "Spin polarisation", "No"
      write(stdOut, "(A,':',T30,F12.6,/,A,':',T30,F12.6)") "Nr. of up electrons", 0.5_dp*nEl(1),&
          & "Nr. of down electrons", 0.5_dp*nEl(1)
    case(2)
      write(stdOut, "(A,':',T30,A)") "Spin polarisation", "Yes"
      write(stdOut, "(A,':',T30,F12.6,/,A,':',T30,F12.6)") "Nr. of up electrons", nEl(1),&
          & "Nr. of down electrons", nEl(2)
    case(4)
      write(stdOut, "(A,':',T30,A)") "Non-collinear calculation", "Yes"
      write(stdOut, "(A,':',T30,F12.6)") "Nr. of electrons", nEl(1)
    end select

    if (tPeriodic) then
      write(stdOut, "(A,':',T30,A)") "Periodic boundaries", "Yes"
      if (tLatOpt) then
        write(stdOut, "(A,':',T30,A)") "Lattice optimisation", "Yes"
        write(stdOut, "(A,':',T30,f12.6)") "Pressure", extPressure
      end if
    else
      write(stdOut, "(A,':',T30,A)") "Periodic boundaries", "No"
    end if

    select case (solver)
    case(1)
      write (strTmp, "(A)") "Standard"
    case(2)
      write (strTmp, "(A)") "Divide and Conquer"
    case(3)
      write (strTmp, "(A)") "Relatively robust (version 1)"
    case(4)
      write (strTmp, "(A)") "Relativel robust (version 2)"
    case default
      call error("Unknown eigensolver!")
    end select
    write(stdOut, "(A,':',T30,A)") "Diagonalizer", trim(strTmp)

    if (tSccCalc) then
      select case (iMixer)
      case(1)
        write (strTmp, "(A)") "Simple"
      case(2)
        write (strTmp, "(A)") "Anderson"
      case(3)
        write (strTmp, "(A)") "Broyden"
      case(4)
        write (strTmp, "(A)") "DIIS"
      end select
      write(stdOut, "(A,':',T30,A,' ',A)") "Mixer", trim(strTmp), "mixer"
      write(stdOut, "(A,':',T30,F14.6)") "Mixing parameter", mixParam
      write(stdOut, "(A,':',T30,I14)") "Maximal SCC-cycles", maxSccIter
      select case (iMixer)
      case(2)
        write(stdOut, "(A,':',T30,I14)") "Nr. of chrg. vectors to mix", nGeneration
      case(3)
        write(stdOut, "(A,':',T30,I14)") "Nr. of chrg. vec. in memory", nGeneration
      case(4)
        write(stdOut, "(A,':',T30,I14)") "Nr. of chrg. vectors to mix", nGeneration
      end select
    end if

    if (tCoordOpt) then
      write(stdOut, "(A,':',T30,I14)") "Nr. of moved atoms", nMovedAtom
    end if
    if (tGeoOpt) then
      write(stdOut, "(A,':',T30,I14)") "Max. nr. of geometry steps", nGeoSteps
      write(stdOut, "(A,':',T30,E14.6)") "Force tolerance", input%ctrl%maxForce
      if (input%ctrl%iGeoOpt == 1) then
        write(stdOut, "(A,':',T30,E14.6)") "Step size", deltaT
      end if
    end if

    if (tForces) then
      select case (forceType)
      case(0)
        strTmp = "Traditional"
      case(2)
        strTmp = "Dynamics, zero electronic temp."
      case(3)
        strTmp = "Dynamics, finite electronic temp."
      end select
      write(stdOut, "(A,':',T30,A)") "Force evaluation method", strTmp
    end if

    tFirst = .true.
    if (allocated(conAtom)) then
      do ii = 1, nAtom
        do jj = 1, size(conAtom)
          if (conAtom(jj) == ii) then
            if (tFirst) then
              write(strTmp, "(A,':')") "Geometry constraints"
              tFirst = .false.
            else
              write(strTmp, "(A)") ""
            end if
            write(stdOut, "(A,T30,'At',I4,': ',3F10.6)") trim(strTmp), ii, (conVec(kk,jj), kk=1,3)
          end if
        end do
      end do
    end if

    if (.not.input%ctrl%tSetFillingTemp) then
      write(stdOut, "(A,':',T30,E14.6)") "Electronic temperature", tempElec
    end if
    if (tMD) then
      write(stdOut, "(A,':',T30,E14.6)") "Time step", deltaT
      if (input%ctrl%iThermostat == 0 .and. .not.input%ctrl%tReadMDVelocities) then
        write(stdOut, "(A,':',T30,E14.6)") "Temperature", tempAtom
      end if
      if (input%ctrl%tRescale) then
        write(stdOut, "(A,':',T30,E14.6)") "Rescaling probability", input%ctrl%wvScale
      end if
    end if

    if (tSccCalc) then
      if (input%ctrl%tReadChrg) then
        write (strTmp, "(A,A,A)") "Read in from '", trim(fCharges), "'"
      else
        write (strTmp, "(A,E11.3,A)") "Set automatically (system chrg: ", input%ctrl%nrChrg, ")"
      end if
      write(stdOut, "(A,':',T30,A)") "Initial charges", trim(strTmp)
    end if

    do iSp = 1, nType
      if (iSp == 1) then
        write (strTmp, "(A,':')") "Included shells"
      else
        write (strTmp, "(A)") ""
      end if
      do jj = 1, orb%nShell(iSp)
        if (jj == 1) then
          strTmp2 = trim(orbitalNames(orb%angShell(jj, iSp) + 1))
        else
          strTmp2 = trim(strTmp2) // ", " // trim(orbitalNames(orb%angShell(jj, iSp) + 1))
        end if
      end do
      write(stdOut, "(A,T29,A2,':  ',A)") trim(strTmp), trim(speciesName(iSp)), trim(strTmp2)
    end do

    if (tPeriodic) then
      do ii = 1, nKPoint
        if (ii == 1) then
          write(strTmp, "(A,':')") "K-points and weights"
        else
          write(strTmp, "(A)") ""
        end if
        write(stdOut, "(A,T28,I6,':',3F10.6,3X,F10.6)") trim(strTmp), ii,&
            & (kPoint(jj, ii), jj=1, 3), kWeight(ii)
      end do
      write(stdout,*)
      do ii = 1, nKPoint
        if (ii == 1) then
          write(strTmp, "(A,':')") "K-points in absolute space"
        else
          write(strTmp, "(A)") ""
        end if
        write(stdout, "(A,T28,I6,':',3F10.6)") trim(strTmp), ii, matmul(invLatVec,kPoint(:,ii))
      end do
      write(stdout, *)
    end if

    if (tDispersion) then
      select type (dispersion)
      type is (DispSlaKirk)
        write(stdOut, "(A)") "Using Slater-Kirkwood dispersion corrections"
      type is (DispUff)
        write(stdOut, "(A)") "Using Lennard-Jones dispersion corrections"
    #:if WITH_DFTD3
      type is (DispDftD3)
        write(stdOut, "(A)") "Using DFT-D3 dispersion corrections"
    #:endif
      class default
        call error("Unknown dispersion model - this should not happen!")
      end select
    end if

    if (tSccCalc) then
      ! Have the SK values of U been replaced?
      if (allocated(input%ctrl%hubbU)) then
        strTmp = ""
        tFirst = .true.
        do iSp = 1, nType
          if (all(input%ctrl%hubbU(1:orb%nShell(iSp), iSp) == 0.0_dp)) then
            cycle
          end if
          do jj = 1, orb%nShell(iSp)
            if (tFirst) then
              write(strTmp, "(A,':')") "Non-default Hubbard U"
              tFirst = .false.
            else
              write(strTmp, "(A)") ""
            end if
            write(stdOut, "(A,T30,A2,2X,I1,'(',A1,'): ',E14.6)") trim(strTmp), speciesName(iSp),&
                & jj, orbitalNames(orb%angShell(jj, iSp)+1), hubbU(jj, iSp)
          end do
        end do
      end if
    end if

    tFirst = .true.
    if (tSpin) then
      do iSp = 1, nType
        do jj = 1, orb%nShell(iSp)
          do kk = 1, orb%nShell(iSp)
            if (tFirst) then
              write(strTmp, "(A,':')") "Spin coupling constants"
              tFirst = .false.
            else
              write(strTmp, "(A)") ""
            end if
            write(stdOut, "(A,T30,A2,2X,I1,'(',A1,')-',I1,'(',A1,'): ',E14.6)")trim(strTmp),&
                & speciesName(iSp), jj, orbitalNames(orb%angShell(jj, iSp)+1), kk,&
                & orbitalNames(orb%angShell(kk, iSp)+1), spinW(kk, jj, iSp)
          end do
        end do
      end do
    end if

    tFirst = .true.
    if (tSpinOrbit) then
      if (tDualSpinOrbit) then
        write(stdOut, "(A)")"Dual representation spin orbit"
      end if
      do iSp = 1, nType
        do jj = 1, orb%nShell(iSp)
          if (tFirst) then
            write(strTmp, "(A,':')") "Spin orbit constants"
            tFirst = .false.
          else
            write(strTmp, "(A)") ""
          end if
          write(stdOut, "(A,T30,A2,2X,I1,'(',A1,'): ',E14.6)")trim(strTmp), speciesName(iSp),&
                & jj, orbitalNames(orb%angShell(jj, iSp)+1), xi(jj, iSp)
          if (xi(jj, iSp) /= 0.0_dp .and. orb%angShell(jj, iSp) == 0) then
            call error("Program halt due to non-zero s-orbital spin-orbit coupling constant!")
          end if
        end do
      end do
    end if

    if (tSccCalc) then
      if (t3rdFull) then
        write(stdOut, "(A,T30,A)") "Full 3rd order correction", "Yes"
        if (input%ctrl%tOrbResolved) then
          write(stdOut, "(A,T30,A)") "Orbital-resolved 3rd order", "Yes"
          write(stdOut, "(A30)") "Shell-resolved Hubbard derivs:"
          write(stdOut, "(A)") "        s-shell   p-shell   d-shell   f-shell"
          do iSp = 1, nType
            write(stdOut, "(A3,A3,4F10.4)") "  ", trim(speciesName(iSp)),&
                & input%ctrl%hubDerivs(:orb%nShell(iSp),iSp)
          end do
        end if
      end if

      if (any(tDampedShort)) then
        write(stdOut, "(A,T30,A)") "Damped SCC", "Yes"
        ii = count(tDampedShort)
        write(strTmp, "(A,I0,A)") "(A,T30,", ii, "(A,1X))"
        write(stdOut, strTmp) "Damped species(s):", pack(speciesName, tDampedShort)
        deallocate(tDampedShort)
      end if

      if (input%ctrl%h5SwitchedOn) then
        write(stdOut, "(A,T30,A)") "H-bond correction:", "H5"
      end if
      if (tHHRepulsion) then
        write(stdOut, "(A,T30,A)") "H-H repulsion correction:", "H5"
      end if
    end if

    write(stdOut, "(A,':')") "Extra options"
    if (tPrintMulliken) then
      write(stdOut, "(T30,A)") "Mulliken analysis"
    end if
    if (tPrintForces .and. .not. (tMD .or. tGeoOpt .or. tDerivs)) then
      write(stdOut, "(T30,A)") "Force calculation"
    end if
    if (tPrintEigVecs) then
      write(stdOut, "(T30,A)") "Eigenvector printing"
    end if
    if (tExtChrg) then
      write(stdOut, "(T30,A)") "External charges specified"
    end if

    if (tEField) then
      if (tTDEfield) then
        write(stdOut, "(T30,A)") "External electric field specified"
        write(stdOut, "(A,':',T30,E14.6)") "Angular frequency", EfieldOmega
      else
        write(stdOut, "(T30,A)") "External static electric field specified"
      end if
      write(stdOut, "(A,':',T30,E14.6)") "Field strength", EFieldStrength
      write(stdOut, "(A,':',T30,3F9.6)") "Direction", EfieldVector
      if (tPeriodic) then
        call warning("Saw tooth potential used for periodic geometry - make sure there is a vacuum&
            & region!")
      end if
    end if

    if (tDFTBU) then
      do iSp = 1, nType
        if (nUJ(iSp)>0) then
          write(strTmp, "(A,':')") "U-J coupling constants"
          write(stdOut, "(A,T25,A2)")trim(strTmp), speciesName(iSp)
          do jj = 1, nUJ(iSp)
            write(strTmp, "(A,I1,A)")'(A,',niUJ(jj,iSp),'I2,T25,A,F6.4)'
            write(stdOut, trim(strTmp))'Shells:',iUJ(1:niUJ(jj,iSp),jj,iSp),'UJ:', UJ(jj,iSp)
          end do
        end if
      end do

    end if

    select case (forceType)
    case(0)
      write(stdOut, "(A,T30,A)") "Force type", "original"
    case(1)
      write(stdOut, "(A,T30,A)") "Force type", "erho with re-diagonalized eigenvalues"
    case(2)
      write(stdOut, "(A,T30,A)") "Force type", "erho with DHD-product (T_elec = 0K)"
    case(3)
      write(stdOut, "(A,T30,A)") "Force type", "erho with S^-1 H D (Te <> 0K)"
    end select

    if (tSpinOrbit .and. tDFTBU .and. .not. tDualSpinOrbit)  then
      call error("Only dual spin orbit currently supported for orbital potentials")
    end if

    if (.not.tStress) then
      if (tBarostat) then
        call error("Sorry, MD with a barostat requires stress evaluation")
      end if
      if (tLatOpt) then
        call error("Sorry, lattice optimization requires stress tensor evaluation")
      end if
    end if

    if (tSpinOrbit .and. (tWriteHS .or.(tWriteRealHS.and..not.tDualSpinOrbit))) then
      call error("Writing of Hamiltonian currently not possible with spin orbit coupling enabled.")
    end if

    if (tLinResp) then
      if (tDFTBU) then
        call error("Linear response is not compatible with Orbitally dependant functionals yet")
      end if

      if (tForces .and. nSpin > 1) then
        call error("Linear response is not available for spin polarised forces yet")
      end if

      if (t2Component) then
        call error("Linear response is not compatibile with this spinor Hamiltonian")
      end if

      if (tStress) then
        call error("Excited state stresses not implemented")
      end if

      if (.not.tRealHS) then
        call error("Linear response does not support k-points")
      end if

    end if

    call env%globalTimer%stopTimer(globalTimers%globalInit)

  end subroutine initProgramVariables


  !> Clean up things that do not automatically get removed on going out of scope
  subroutine destructProgramVariables()

    if (tProjEigenvecs) then
      call destruct(iOrbRegion)
      call destruct(RegionLabels)
    end if

  end subroutine destructProgramVariables


  !> Creates all random generators needed in the code.
  !!
  subroutine createRandomGenerators(env, seed, randomInit, randomThermostat)

    !> Environment settings
    type(TEnvironment), intent(in) :: env

    !> Global seed used for initialisation of the random generator pool. If less than one, random
    !! initialisation of the seed will occur.
    integer, intent(inout) :: seed

    !> Random generator for initprogram.
    type(ORanlux), allocatable, intent(out) :: randomInit

    !> Random generator for the actual thermostat.
    type(ORanlux), allocatable, intent(out) :: randomThermostat

    type(ORandomGenPool) :: randGenPool

    call init(randGenPool, env, seed, oldCompat=.true.)

    ! DO NOT CHANGE the ORDER of calls below, as this will destroy backwards compatibility and
    ! reproduciblity of random number sequences in the code. If further random generators are needed
    ! *append* similar getGenerator() calls. All random generators used within the code must be
    ! generated here.
    call randGenPool%getGenerator(env, randomThermostat)
    call randGenPool%getGenerator(env, randomInit)

  end subroutine createRandomGenerators

#:if WITH_SOCKETS
  !> Initializes the socket and recieves and broadcasts initial geometry.
  subroutine initSocket(env, socketInput, tPeriodic, coord0, latVec, socket, tCoordsChanged,&
      & tLatticeChanged)

    !> Environment settings.
    type(TEnvironment), intent(in) :: env

    !> Input data for the socket.
    type(IpiSocketCommInp), intent(inout) :: socketInput

    !> Is the system periodic?
    logical, intent(in) :: tPeriodic

    !> Received atom coordinates in the unit cell.
    real(dp), intent(inout) :: coord0(:,:)

    !> Received lattice vectors
    real(dp), intent(inout) :: latVec(:,:)

    !> Initialised socket.
    type(IpiSocketComm), allocatable, intent(out) :: socket

    !> Whether coordinates has been changed
    logical, intent(out) :: tCoordsChanged

    !> Whether lattice vectors has been changed
    logical, intent(out) :: tLatticeChanged

    logical :: tDummy

    if (env%tGlobalMaster) then
      write(stdOut, "(A,1X,A)") "Initialising for socket communication to host",&
          & trim(socketInput%host)
      socket = IpiSocketComm(socketInput)
    end if
    call receiveGeometryFromSocket(env, socket, tPeriodic, coord0, latVec, tCoordsChanged,&
        & tLatticeChanged, tDummy)

  end subroutine initSocket
#:endif

  !> Initialises (clears) output files.
  subroutine initOutputFiles(env, tWriteAutotest, tWriteResultsTag, tWriteBandDat, tDerivs,&
      & tWriteDetailedOut, tMd, tGeoOpt, geoOutFile, fdDetailedOut, fdMd, esp)

    !> Environment
    type(TEnvironment), intent(inout) :: env

    !> Should tagged regression test data be printed
    logical, intent(in) :: tWriteAutotest

    !> Write tagged output for machine readable results
    logical, intent(in) :: tWriteResultsTag

    !> Band structure and fillings
    logical, intent(in) :: tWriteBandDat

    !> Finite different second derivatives
    logical, intent(in) :: tDerivs

    !> Write detail on the calculation to file
    logical, intent(in) :: tWriteDetailedOut

    !> Is this a molecular dynamics calculation
    logical, intent(in) :: tMd

    !> Are atomic coodinates being optimised
    logical, intent(in) :: tGeoOpt

    !> Filename for geometry output
    character(*), intent(in) :: geoOutFile

    !> File unit for detailed.out
    integer, intent(out) :: fdDetailedOut

    !> File unit for information during molecular dynamics
    integer, intent(out) :: fdMd

    !> Electrostatic potentials if requested
    type(TElStatPotentials), allocatable, intent(inout) :: esp

    call initTaggedWriter()
    if (tWriteAutotest) then
      call initOutputFile(autotestTag)
    end if
    if (tWriteResultsTag) then
      call initOutputFile(resultsTag)
    end if
    if (tWriteBandDat) then
      call initOutputFile(bandOut)
    end if
    if (tDerivs) then
      call initOutputFile(hessianOut)
    end if
    if (tWriteDetailedOut) then
      call initOutputFile(userOut, fdDetailedOut)
      call env%fileFinalizer%register(fdDetailedOut)
    end if
    if (tMD) then
      call initOutputFile(mdOut, fdMD)
      call env%fileFinalizer%register(fdMd)
    end if
    if (tGeoOpt .or. tMD) then
      call clearFile(trim(geoOutFile) // ".gen")
      call clearFile(trim(geoOutFile) // ".xyz")
    end if
    if (allocated(esp)) then
      call initOutputFile(esp%espOutFile)
    end if

  end subroutine initOutputFiles


  !> Allocates most of the large arrays needed during the DFTB run.
  subroutine initArrays(env, tForces, tExtChrg, tLinResp, tLinRespZVect, tMd, tMulliken,&
      & tSpinOrbit, tImHam, tWriteRealHS, tWriteHS, t2Component, tRealHS, tPrintExcitedEigvecs,&
      & tDipole, orb, nAtom, nMovedAtom, nKPoint, nSpin, nExtChrg, indMovedAtom, mass, denseDesc,&
      & rhoPrim, h0, iRhoPrim, excitedDerivs, ERhoPrim, derivs, chrgForces, energy, potential, TS,&
      & E0, Eband, eigen, filling, coord0Fold, newCoords, orbitalL, HSqrCplx, SSqrCplx,&
      & eigvecsCplx, HSqrReal, SSqrReal, eigvecsReal, rhoSqrReal, chargePerShell, occNatural,&
      & velocities, movedVelo, movedAccel, movedMass, dipoleMoment)

    !> Current environment
    type(TEnvironment), intent(in) :: env

    !> Are forces required
    logical, intent(in) :: tForces

    !> Are the external charges
    logical, intent(in) :: tExtChrg

    !> Are excitation energies being calculated
    logical, intent(in) :: tLinResp

    !> Are excited state properties being calculated
    logical, intent(in) :: tLinRespZVect

    !> Is this a molecular dynamics calculation
    logical, intent(in) :: tMd

    !> Is Mulliken analysis being performed
    logical, intent(in) :: tMulliken

    !> Are there spin orbit interactions
    logical, intent(in) :: tSpinOrbit

    !> Are there imaginary parts to the hamiltonian
    logical, intent(in) :: tImHam

    !> Should the sparse hamiltonian and overlap be writen to disc?
    logical, intent(in) :: tWriteRealHS

    !> Should the hamiltonian and overlap be written out as dense matrices
    logical, intent(in) :: tWriteHS

    !> Is the hamiltonian for a two component (Pauli) system
    logical, intent(in) :: t2Component

    !> Is the hamiltonian real?
    logical, intent(in) :: tRealHS

    !> Print the excited state eigenvectors
    logical, intent(in) :: tPrintExcitedEigvecs

    !> Print the dipole moment
    logical, intent(in) :: tDipole

    !> data structure with atomic orbital information
    type(TOrbitals), intent(in) :: orb

    !> Number of atoms
    integer, intent(in) :: nAtom

    !> Number of atoms moved about during the calculation
    integer, intent(in) :: nMovedAtom

    !> Number of k-points
    integer, intent(in) :: nKPoint

    !> Number of spin channels
    integer, intent(in) :: nSpin

    !> Number of external charges
    integer, intent(in) :: nExtChrg

    !> Indices for any moving atoms
    integer, intent(in) :: indMovedAtom(:)

    !> Masses of each species of atom
    real(dp), intent(in) :: mass(:)

    !> Dense matrix descriptor for H and S
    type(TDenseDescr), intent(in) :: denseDesc

    !> Sparse storage density matrix
    real(dp), intent(out), allocatable :: rhoPrim(:,:)

    !> Non-scc part of the hamiltonian
    real(dp), intent(out), allocatable :: h0(:)

    !> Imaginary part of the sparse density matrix
    real(dp), intent(out), allocatable :: iRhoPrim(:,:)

    !> Excitation energy derivatives with respect to atomic coordinates
    real(dp), intent(out), allocatable :: excitedDerivs(:,:)

    !> Energy weighted density matrix
    real(dp), intent(out), allocatable :: ERhoPrim(:)

    !> Derivatives of total energy with respect to atomic coordinates
    real(dp), intent(out), allocatable :: derivs(:,:)

    !> Forces on (any) external charges
    real(dp), intent(out), allocatable :: chrgForces(:,:)

    !> Energy terms
    type(TEnergies), intent(out) :: energy

    !> Potentials acting on the system
    type(TPotentials), intent(out) :: potential

    !> Electron entropy contribution at T
    real(dp), intent(out), allocatable :: TS(:)

    !> zero temperature extrapolated electronic energy
    real(dp), intent(out), allocatable :: E0(:)

    !> band  energy
    real(dp), intent(out), allocatable :: Eband(:)

    !> single particle energies (band structure)
    real(dp), intent(out), allocatable :: eigen(:,:,:)

    !> Occupations of single particle states
    real(dp), intent(out), allocatable :: filling(:,:,:)

    !> Coordinates in central cell
    real(dp), intent(out), allocatable :: coord0Fold(:,:)

    !> Updated coordinates
    real(dp), intent(out), allocatable :: newCoords(:,:)

    !> Orbital angular momentum
    real(dp), intent(out), allocatable :: orbitalL(:,:,:)

    !> Complex dense hamiltonian
    complex(dp), intent(out), allocatable :: HSqrCplx(:,:)

    !> overlap matrix dense storage
    complex(dp), intent(out), allocatable :: SSqrCplx(:,:)

    !> Complex eigenvectors
    complex(dp), intent(out), allocatable :: eigvecsCplx(:,:,:)

    !> real dense hamiltonian
    real(dp), intent(out), allocatable :: HSqrReal(:,:)

    !> overlap matrix dense storage
    real(dp), intent(out), allocatable :: SSqrReal(:,:)

    !> Real eigenvectors
    real(dp), intent(out), allocatable :: eigvecsReal(:,:,:)

    !> density matrix dense storage
    real(dp), intent(out), allocatable :: rhoSqrReal(:,:,:)

    !> Number of electron in each atomic shell
    real(dp), intent(out), allocatable :: chargePerShell(:,:,:)

    !> Occupations for natural orbitals
    real(dp), intent(out), allocatable :: occNatural(:)

    !> Atomic velocities
    real(dp), intent(out), allocatable :: velocities(:,:)

    !> Array for moving atom velocities
    real(dp), intent(out), allocatable :: movedVelo(:,:)

    !> moving atom accelerations
    real(dp), intent(out), allocatable :: movedAccel(:,:)

    !> moving atoms masses
    real(dp), intent(out), allocatable :: movedMass(:,:)

    !> system dipole moment
    real(dp), intent(out), allocatable :: dipoleMoment(:)


    integer :: nSpinHams, sqrHamSize

    allocate(rhoPrim(0, nSpin))
    allocate(h0(0))
    if (tImHam) then
      allocate(iRhoPrim(0, nSpin))
    end if

    allocate(excitedDerivs(0,0))
    if (tForces) then
      allocate(ERhoPrim(0))
      allocate(derivs(3, nAtom))
      if (tExtChrg) then
        allocate(chrgForces(3, nExtChrg))
      end if
      if (tLinRespZVect) then
        deallocate(excitedDerivs)
        allocate(excitedDerivs(3, nAtom))
      end if
    end if

    call init(energy, nAtom)
    call init(potential, orb, nAtom, nSpin)

    ! Nr. of independent spin Hamiltonians
    select case (nSpin)
    case (1)
      nSpinHams = 1
    case (2)
      nSpinHams = 2
    case (4)
      nSpinHams = 1
    end select

    sqrHamSize = denseDesc%fullSize
    allocate(TS(nSpinHams))
    allocate(E0(nSpinHams))
    allocate(Eband(nSpinHams))
    allocate(eigen(sqrHamSize, nKPoint, nSpinHams))
    allocate(filling(sqrHamSize, nKpoint, nSpinHams))

    allocate(coord0Fold(3, nAtom))

    if (tMD) then
      allocate(newCoords(3, nAtom))
    end if

    if ((tMulliken .and. tSpinOrbit) .or. tImHam) then
      allocate(orbitalL(3, orb%mShell, nAtom))
    end if

    ! If only H/S should be printed, no allocation for square HS is needed
    if (.not. (tWriteRealHS .or. tWriteHS)) then
      call allocateDenseMatrices(env, denseDesc, parallelKS%localKS, t2Component, tRealHS,&
          & HSqrCplx, SSqrCplx, eigVecsCplx, HSqrReal, SSqrReal, eigvecsReal)
    end if

    if (tLinResp) then
      if (withMpi) then
        call error("Linear response calc. does not work with MPI yet")
      end if
      if (tLinRespZVect) then
        allocate(rhoSqrReal(sqrHamSize, sqrHamSize, nSpin))
      end if
    end if
    allocate(chargePerShell(orb%mShell, nAtom, nSpin))

    if (tLinResp .and. tPrintExcitedEigVecs) then
      allocate(occNatural(orb%nOrb))
    end if

    if (tMD) then
      allocate(velocities(3, nAtom))
      allocate(movedVelo(3, nMovedAtom))
      allocate(movedAccel(3, nMovedAtom))
      allocate(movedMass(3, nMovedAtom))
      movedMass(:,:) = spread(mass(indMovedAtom),1,3)
    end if

    if (tDipole) then
      allocate(dipoleMoment(3))
    end if

  end subroutine initArrays


  !> Set up storage for dense matrices, either on a single processor, or as BLACS matrices
  subroutine allocateDenseMatrices(env, denseDesc, localKS, t2Component, tRealHS, HSqrCplx,&
      & SSqrCplx, eigvecsCplx, HSqrReal, SSqrReal, eigvecsReal)

    !> Computing environment
    type(TEnvironment), intent(in) :: env

    !> Descriptor of the large square matrices in the program
    type(TDenseDescr), intent(in) :: denseDesc

    !> Index array for spin and k-point index
    integer, intent(in) :: localKS(:,:)

    !> Is this a two component calculation
    logical, intent(in) :: t2Component

    !> Is this a real hamiltonian
    logical, intent(in) :: tRealHS

    !> Square H matrix
    complex(dp), allocatable, intent(out) :: HSqrCplx(:,:)

    !> Square S matrix
    complex(dp), allocatable, intent(out) :: SSqrCplx(:,:)

    !> Eigenvectors for complex eigenproblem
    complex(dp), allocatable, intent(out) :: eigvecsCplx(:,:,:)

    !> Square H matrix
    real(dp), allocatable, intent(out) :: HSqrReal(:,:)

    !> Square S matrix
    real(dp), allocatable, intent(out) :: SSqrReal(:,:)

    !> Eigenvectors for real eigenproblem
    real(dp), allocatable, intent(out) :: eigvecsReal(:,:,:)

    integer :: nLocalCols, nLocalRows, nLocalKS

    nLocalKS = size(localKS, dim=2)
  #:if WITH_SCALAPACK
    call scalafx_getlocalshape(env%blacs%orbitalGrid, denseDesc%blacsOrbSqr, nLocalRows, nLocalCols)
  #:else
    nLocalRows = denseDesc%fullSize
    nLocalCols = denseDesc%fullSize
  #:endif

    if (t2Component .or. .not. tRealHS) then
      allocate(HSqrCplx(nLocalRows, nLocalCols))
      allocate(SSqrCplx(nLocalRows, nLocalCols))
      allocate(eigVecsCplx(nLocalRows, nLocalCols, nLocalKS))
    else
      allocate(HSqrReal(nLocalRows, nLocalCols))
      allocate(SSqrReal(nLocalRows, nLocalCols))
      allocate(eigVecsReal(nLocalRows, nLocalCols, nLocalKS))
    end if

  end subroutine allocateDenseMatrices


#:if WITH_SCALAPACK
  #!
  #! SCALAPACK related routines
  #!

  !> Initialise parallel large matrix decomposition methods
  subroutine initScalapack(blacsOpts, nAtom, nOrb, t2Component, env)

    !> BLACS settings
    type(TBlacsOpts), intent(in) :: blacsOpts

    !> Number of atoms
    integer, intent(in) :: nAtom

    !> Number of orbitals
    integer, intent(in) :: nOrb

    !> Is this a two component calculation
    logical, intent(in) :: t2Component

    !> Computing enviroment data
    type(TEnvironment), intent(inout) :: env

    integer :: sizeHS

    if (t2Component) then
      sizeHS = 2 * nOrb
    else
      sizeHS = nOrb
    end if
    call env%initBlacs(blacsOpts%blockSize, blacsOpts%blockSize, sizeHS, nAtom)

  end subroutine initScalapack


  !> Generate descriptions of large dense matrices in BLACS decomposition
  !>
  !> Note: It must be called after getDenseDescCommon() has been called.
  !>
  subroutine getDenseDescBlacs(env, rowBlock, colBlock, denseDesc)

    !> parallel environment
    type(TEnvironment), intent(in) :: env

    !> Number of matrix rows
    integer, intent(in) :: rowBlock

    !> Number of matrix columns
    integer, intent(in) :: colBlock

    !> Descriptor of the dense matrix
    type(TDenseDescr), intent(inout) :: denseDesc

    integer :: nn

    nn = denseDesc%fullSize
    call scalafx_getdescriptor(env%blacs%orbitalGrid, nn, nn, rowBlock, colBlock,&
        & denseDesc%blacsOrbSqr)

  end subroutine getDenseDescBlacs

#:endif


  !> Generate description of the total large square matrices, on the basis of atomic orbital
  !> orderings
  !>
  subroutine getDenseDescCommon(orb, nAtom, t2Component, denseDesc)

    !> Orbital information for species
    type(TOrbitals), intent(in) :: orb

    !> Number of atoms in the system
    integer, intent(in) :: nAtom

    !> Is this a two component calculation
    logical, intent(in) :: t2Component

    !> Resulting descriptor
    type(TDenseDescr), intent(out) :: denseDesc

    integer :: nOrb

    allocate(denseDesc%iAtomStart(nAtom + 1))
    call buildSquaredAtomIndex(denseDesc%iAtomStart, orb)
    nOrb = denseDesc%iAtomStart(nAtom + 1) - 1
    denseDesc%t2Component = t2Component
    denseDesc%nOrb = nOrb
    if (t2Component) then
      denseDesc%fullSize = 2 * nOrb
    else
      denseDesc%fullSize = nOrb
    end if

  end subroutine getDenseDescCommon



end module initprogram
