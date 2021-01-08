!--------------------------------------------------------------------------------------------------!
!  DFTB+: general package for performing fast atomistic simulations                                !
!  Copyright (C) 2006 - 2020  DFTB+ developers group                                               !
!                                                                                                  !
!  See the LICENSE file for terms of usage and distribution.                                       !
!--------------------------------------------------------------------------------------------------!

#:include 'common.fypp'


!> Global variables and initialization for the main program.
module dftbp_initprogram
#:if WITH_OMP
  use omp_lib
#:endif
  use dftbp_mainio, only : initOutputFile
  use dftbp_assert
  use dftbp_globalenv
  use dftbp_coherence
  use dftbp_environment
  use dftbp_scalapackfx
  use dftbp_inputdata
  use dftbp_densedescr
  use dftbp_constants
  use dftbp_elecsolvers
  use dftbp_elsisolver, only : TElsiSolver_init, TElsiSolver_final
  use dftbp_elsiiface
  use dftbp_arpack, only : withArpack
  use dftbp_periodic, only : TNeighbourList, TNeighbourlist_init, buildSquaredAtomIndex
  use dftbp_periodic, only : getCellTranslations
  use dftbp_accuracy
  use dftbp_intrinsicpr
  use dftbp_shortgamma
  use dftbp_message
  use dftbp_mixer
  use dftbp_simplemixer
  use dftbp_andersonmixer
  use dftbp_broydenmixer
  use dftbp_diismixer

  use dftbp_geoopt
  use dftbp_conjgrad
  use dftbp_steepdesc
  use dftbp_gdiis
  use dftbp_lbfgs
  use dftbp_fire

  use dftbp_hamiltoniantypes

  use dftbp_randomgenpool
  use dftbp_ranlux
  use dftbp_mdcommon
  use dftbp_mdintegrator
  use dftbp_velocityverlet
  use dftbp_thermostat
  use dftbp_dummytherm
  use dftbp_andersentherm
  use dftbp_berendsentherm
  use dftbp_nhctherm
  use dftbp_tempprofile, only : TTempProfile, TempProfile_init
  use dftbp_numderivs2
  use dftbp_lapackroutines
  use dftbp_simplealgebra
  use dftbp_nonscc
  use dftbp_scc
  use dftbp_sccinit
  use dftbp_onsitecorrection
  use dftbp_hamiltonian, only : TRefExtPot
  use dftbp_h5correction
  use dftbp_halogenx
  use dftbp_slakocont
  use dftbp_repcont
  use dftbp_fileid
  use dftbp_spin, only: Spin_getOrbitalEquiv, ud2qm, qm2ud
  use dftbp_dftbplusu, only : TDftbU, TDftbU_init
  use dftbp_blockpothelper, only : appendBlockReduced
  use dftbp_dispersions
  use dftbp_thirdorder
  use dftbp_linresp
  use dftbp_linresptypes
  use dftbp_pprpa, only : TppRPAcal
  use dftbp_RangeSeparated
  use dftbp_stress
  use dftbp_orbitalequiv
  use dftbp_orbitals
  use dftbp_commontypes
  use dftbp_linkedlist
  use dftbp_wrappedintr
  use dftbp_timeprop
  use dftbp_xlbomd
  use dftbp_etemp, only : fillingTypes
#:if WITH_SOCKETS
  use dftbp_mainio, only : receiveGeometryFromSocket
  use dftbp_ipisocket
#:endif
  use dftbp_elstatpot
  use dftbp_pmlocalisation
  use dftbp_energytypes
  use dftbp_potentials
  use dftbp_taggedoutput
  use dftbp_formatout
  use dftbp_qdepextpotproxy, only : TQDepExtPotProxy
  use dftbp_forcetypes, only : forceTypes
  use dftbp_elstattypes, only : elstatTypes
  use dftbp_reks
  use dftbp_plumed, only : withPlumed, TPlumedCalc, TPlumedCalc_init
  use dftbp_cm5, only : TChargeModel5, TChargeModel5_init
  use dftbp_solvation, only : TSolvation
  use dftbp_solvinput, only : createSolvationModel, writeSolvationInfo

#:if WITH_TRANSPORT
  use libnegf_vars
  use negf_int
#:endif
  use poisson_init
  use dftbp_transportio
  use dftbp_determinants
  implicit none

#:if not WITH_TRANSPORT

  !> Dummy type for negf interface
  type :: TNegfInt
  end type TNegfInt

#:endif

! private

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

  !> file name for shift data
  character(*), parameter :: fShifts = "shifts.dat"


  !> Interaction cutoff distances
  type, public :: TCutoffs

    real(dp) :: skCutOff
    real(dp) :: repCutOff
    real(dp) :: lcCutOff
    real(dp) :: mCutOff

  end type TCutoffs


  type, public :: TDftbPlusMain

    !> Is this calculation using a restarted input that does not require self consistency before
    !> moving to the post SCC loop part (i.e. Ehrenfest)
    logical :: tRestartNoSC = .false.

    !> Is the calculation SCC?
    logical :: tSccCalc

    !> SCC module internal variables
    type(TScc), allocatable :: sccCalc

    !> nr. of atoms
    integer :: nAtom

    !> nr. of all (boundary condition images and original) atoms
    integer :: nAllAtom

    !> nr. of original atom in central cell
    integer, allocatable :: img2CentCell(:)

    !> nr of different types (nAtom)
    integer :: nType

    !> data type for atomic orbital information
    type(TOrbitals) :: orb

    !> nr. of orbitals in the system
    integer :: nOrb

    !> types of the atoms (nAllAtom)
    integer, allocatable :: species(:)

    !> type of the atoms (nAtom)
    integer, allocatable :: species0(:)

    !> Coords of the atoms (3, nAllAtom)
    real(dp), allocatable :: coord(:,:)

    !> Coords in central cell (3, nAtom)
    real(dp), allocatable :: coord0(:,:)

    !> if calculation is periodic
    logical :: tPeriodic

    !> If the calculation is helical geometry
    logical :: tHelical

    !> Should central cell coordinates be output?
    logical :: tShowFoldedCoord

    !> How to calculate forces
    integer :: forceType

    !> are atomic coordinates fractional?
    logical :: tFracCoord

    !> Tolerance for SCC cycle
    real(dp) :: sccTol

    !> lattice vectors as columns
    real(dp), allocatable :: latVec(:,:)

    !> Origin of coordinate system for periodic systems
    real(dp), allocatable :: origin(:)

    !> reciprocal lattice vectors as columns
    real(dp), allocatable :: recVec(:,:)

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
    type(TNeighbourList), allocatable :: neighbourList

    !> nr. of neighbours for atoms out to max interaction distance (excluding Ewald terms)
    integer, allocatable :: nNeighbourSK(:)

    !> nr. of neighbours for atoms within Erep interaction distance (usually short)
    integer, allocatable :: nNeighbourRep(:)

    !> Number of neighbours for each of the atoms for the exchange contributions in the long range
    !> functional
    integer, allocatable :: nNeighbourLC(:)

    !> H/S sparse matrices indexing array for atomic blocks
    integer, allocatable :: iSparseStart(:,:)

    !> Hubbard Us (orbital, atom)
    real(dp), allocatable :: hubbU(:,:)

    !> self energy (orbital, atom)
    real(dp), allocatable :: atomEigVal(:,:)

    !> reference n_0 charges for each atom
    real(dp), allocatable :: referenceN0(:,:)

    !> list of atomic masses
    real(dp), allocatable :: mass(:)

    !> list of atomic masses for each species
    real(dp), allocatable :: speciesMass(:)

    !> Hamiltonian type
    integer :: hamiltonianType

    !> Raw H^0 hamiltonian data
    type(TSlakoCont) :: skHamCont

    !> Raw overlap hamiltonian data
    type(TSlakoCont) :: skOverCont

    !> Repulsive interaction raw data
    type(TRepCont) :: pRepCont

    !> Cut off distances for various types of interaction
    type(TCutoffs) :: cutOff

    !> Cut off distance for repulsive interactions
    real(dp) :: repCutOff

    !> Sparse hamiltonian matrix
    real(dp), allocatable :: ham(:,:)

    !> imaginary part of the Hamiltonian
    real(dp), allocatable :: iHam(:,:)

    !> Charge per atomic shell (shell, atom, spin channel)
    real(dp), allocatable :: chargePerShell(:,:,:)

    !> Charge par atom (atom, spin channel)
    real(dp), allocatable :: chargePerAtom(:,:)

    !> sparse overlap
    real(dp), allocatable :: over(:)


    !> nr. of K-points
    integer :: nKPoint

    !> K-points
    real(dp), allocatable :: kPoint(:,:)

    !> weight of the K-Points
    real(dp), allocatable :: kWeight(:)


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

    !> DFTB+U calculation, if relevant
    type(TDftbU), allocatable :: dftbU

    !> electron temperature
    real(dp) :: tempElec

    !> If K points should filled separately
    logical :: tFillKSep

    !> Fix Fermi energy at specified value
    logical :: tFixEf

    !> Fermi energy per spin
    real(dp), allocatable :: Ef(:)

    !> Filling temp updated by MD.
    logical :: tSetFillingTemp

    !> Choice of electron distribution function, defaults to Fermi
    integer :: iDistribFn = fillingTypes%Fermi

    !> atomic kinetic temperature
    real(dp) :: tempAtom

    !> MD stepsize
    real(dp) :: deltaT

    !> maximal number of SCC iterations
    integer :: maxSccIter

    !> Minimal number of SCC iterations
    integer :: minSccIter

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
    logical :: isGeoOpt

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

    !> Do we need to show net atomic charges?
    logical :: tNetAtomCharges

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

    !> Is the contribution from an excited state needed for the forces
    logical :: tCasidaForces

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
    type(ipiSocketComm), allocatable :: socket
  #:endif

    !> File containing output geometry
    character(lc) :: geoOutFile

    !> Append geometries in the output?
    logical :: tAppendGeo

    !> Use converged SCC charges for properties like forces or charge dependent dispersion
    logical :: isSccConvRequired

    !> labels of atomic species
    character(mc), allocatable :: speciesName(:)

    !> General geometry optimizer
    type(TGeoOpt), allocatable :: pGeoCoordOpt

    !> Geometry optimizer for lattice consts
    type(TGeoOpt), allocatable :: pGeoLatOpt

    !> Charge mixer
    type(TMixer), allocatable :: pChrgMixer

    !> MD Framework
    type(TMDCommon), allocatable :: pMDFrame

    !> MD integrator
    type(TMDIntegrator), allocatable :: pMDIntegrator

    !> Temperature profile driver in MD
    type(TTempProfile), allocatable :: temperatureProfile

    !> geometry optimiser
    type(TNumDerivs), allocatable :: derivDriver

    !> Total charge
    real(dp) :: nrChrg

    !> Spin polarisation
    real(dp) :: nrSpinPol

    !> Is the check-sum for charges read externally to be used?
    logical :: tSkipChrgChecksum

    !> reference neutral atomic occupations
    real(dp), allocatable :: q0(:, :, :)

    !> shell resolved neutral reference
    real(dp), allocatable :: qShell0(:,:)

    !> input charges (for potentials)
    real(dp), allocatable :: qInput(:, :, :)

    !> output charges
    real(dp), allocatable :: qOutput(:, :, :)

    !> charge differences between input and output charges
    real(dp), allocatable :: qDiff(:, :, :)

    !> net (on-site only contributions) charge per atom
    real(dp), allocatable :: qNetAtom(:)

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

    !> nr. of elements to go through the mixer - may contain reduced orbitals and also orbital
    !> blocks (if tDFTBU or onsite corrections)
    integer :: nMixElements

    !> Orbital equivalency for orbital blocks
    integer, allocatable :: iEqBlockDFTBU(:,:,:,:)

    !> Equivalences for onsite block corrections if needed
    integer, allocatable :: iEqBlockOnSite(:,:,:,:)

    !> Orbital equivalency for orbital blocks with spin-orbit
    integer, allocatable :: iEqBlockDFTBULS(:,:,:,:)

    !> Equivalences for onsite block corrections if needed with spin orbit
    integer, allocatable :: iEqBlockOnSiteLS(:,:,:,:)


    ! External charges

    !> If external charges must be considered
    logical :: tExtChrg

    !> Nr. of external charges
    integer :: nExtChrg

    !> external electric field
    logical :: tEField

    !> Arbitrary external field (including electric)
    logical :: tExtField

    !> field strength
    real(dp) :: EFieldStrength

    !> field direction
    real(dp) :: EfieldVector(3)

    !> time dependent
    logical :: tTDEfield

    !> angular frequency
    real(dp) :: EfieldOmega

    !> phase of field at step 0
    integer :: EfieldPhase


    !> Partial density of states (PDOS) projection regions
    type(TListIntR1) :: iOrbRegion

    !> PDOS region labels
    type(TListCharLc) :: regionLabels

    !> Third order DFTB
    logical :: t3rd

    !> Full 3rd order or only atomic site
    logical :: t3rdFull

    !> data structure for 3rd order
    type(TThirdOrder), allocatable :: thirdOrd

    !> Correction to energy from on-site matrix elements
    real(dp), allocatable :: onSiteElements(:,:,:,:)

    !> Correction to dipole momements on-site matrix elements
    real(dp), allocatable :: onSiteDipole(:,:)

    !> Should block charges be mixed as well as charges
    logical :: tMixBlockCharges

    !> Calculate Casida linear response excitations
    logical :: isLinResp

    !> calculate Z vector for excited properties
    logical :: tLinRespZVect

    !> data type for pp-RPA
    type(TppRPAcal), allocatable :: ppRPA

    !> Print eigenvectors
    logical :: tPrintExcitedEigVecs

    !> data type for linear response
    type(TLinResp), allocatable :: linearResponse

    !> Whether to run a range separated calculation
    logical :: isRangeSep

    !> Range Separation data
    type(TRangeSepFunc), allocatable :: rangeSep

    !> DeltaRho input for calculation of range separated Hamiltonian
    real(dp), allocatable :: deltaRhoIn(:)

    !> DeltaRho output from calculation of range separated Hamiltonian
    real(dp), allocatable :: deltaRhoOut(:)

    !> Holds change in deltaRho between SCC steps for range separation
    real(dp), allocatable :: deltaRhoDiff(:)

    !> DeltaRho input for range separation in matrix form
    real(dp), pointer :: deltaRhoInSqr(:,:,:) => null()

    !> DeltaRho output from range separation in matrix form
    real(dp), pointer :: deltaRhoOutSqr(:,:,:) => null()

    !> Linear response calculation with range-separated functional
    logical :: isRS_LinResp

    !> If initial charges/dens mtx. from external file.
    logical :: tReadChrg

    !> Whether potential shifts are read from file
    logical :: tReadShifts

    !> Should charges be read in ascii format?
    logical :: tReadChrgAscii

    !> Whether potential shifts are read from file
    logical :: tWriteShifts

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

    !> dispersion data and calculations
    class(TDispersionIface), allocatable :: dispersion

    !> Solvation data and calculations
    class(TSolvation), allocatable :: solvation

    !> Charge Model 5 for printout
    type(TChargeModel5), allocatable :: cm5Cont

    !> Can stress be calculated?
    logical :: tStress

    !> should XLBOMD be used in MD
    logical :: isXlbomd

    !> XLBOMD related parameters
    type(TXLBOMD), allocatable :: xlbomdIntegrator

    !> Differentiation method for (H^0,S)
    type(TNonSccDiff) :: nonSccDeriv

    !> Whether lattice has changed since last geometry iteration
    logical :: tLatticeChanged

    !> Whether atomic coordindates have changed since last geometry iteration
    logical :: tCoordsChanged

    !> Plumed calculator
    type(TPlumedCalc), allocatable :: plumedCalc

    !> Dense matrix descriptor for H and S
    type(TDenseDescr), allocatable :: denseDesc

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

    !> Total energy components (potentially for multiple determinants)
    type(TEnergies), allocatable :: dftbEnergy(:)

    !> Potentials for orbitals
    type(TPotentials) :: potential

    !> Reference external potential (usual provided via API)
    type(TRefExtPot) :: refExtPot

    !> Proxy for querying population dependant external potenials
    type(TQDepExtPotProxy), allocatable :: qDepExtPot

    !> Energy derivative with respect to atomic positions
    real(dp), allocatable :: derivs(:,:)

    !> Energy derivative for ground state determinant
    real(dp), allocatable :: groundDerivs(:,:)

    !> Energy derivative for triplet determinant (TI-DFTB excited states)
    real(dp), allocatable :: tripletDerivs(:,:)

    !> Energy derivative for mixed determinant (TI-DFTB excited states)
    real(dp), allocatable :: mixedDerivs(:,:)

    !> Forces on any external charges
    real(dp), allocatable :: chrgForces(:,:)

    !> excited state force addition
    real(dp), allocatable :: excitedDerivs(:,:)

    !> dipole moments, when available, for whichever determinants are present
    real(dp), allocatable :: dipoleMoment(:, :)

    !> Coordinates to print out
    real(dp), pointer :: pCoord0Out(:,:)

    !> Folded coords (3, nAtom)
    real(dp), allocatable :: coord0Fold(:,:)

    !> New coordinates returned by the MD routines
    real(dp), allocatable :: newCoords(:,:)

    !> Orbital angular momentum
    real(dp), allocatable :: orbitalL(:,:,:)

    !> Natural orbitals for excited state density matrix, if requested
    real(dp), allocatable :: occNatural(:)

    !> Dynamical (Hessian) matrix
    real(dp), pointer :: pDynMatrix(:,:)

    !> File descriptor for the human readable output
    integer :: fdDetailedOut

    !> File descriptor for extra MD output
    integer :: fdMD

    !> Contains (iK, iS) tuples to be processed in parallel by various processor groups
    type(TParallelKS) :: parallelKS

    !> Electron dynamics
    type(TElecDynamics), allocatable :: electronDynamics

    !> external electric field
    real(dp) :: Efield(3), absEfield

    !> Electronic structure solver
    type(TElectronicSolver) :: electronicSolver

    !> Are large dense matrices required?
    logical :: tLargeDenseMatrices

    !> derivative of cell volume wrt to lattice vectors, needed for pV term
    real(dp) :: extLatDerivs(3,3)

    !> internal pressure within the cell
    real(dp) :: intPressure

    !> Derivative of total energy with respect to lattice vectors
    !> Sign convention: This is in the uphill energy direction for the lattice vectors (each row
    !> pertaining to a separate lattice vector), i.e. opposite to the force.
    !>
    !> The component of a derivative vector that is orthogonal to the plane containing the other two
    !> lattice vectors will expand (contract) the supercell if it is on the opposite (same) same
    !> side of the plane as its associated lattice vector.
    !>
    !> In the special case of cartesian axis aligned orthorhombic lattice vectors, negative diagonal
    !> elements expand the supercell.
    real(dp) :: totalLatDeriv(3,3)

    !> Stress tensors for various contribution in periodic calculations
    !> Sign convention: Positive diagonal elements expand the supercell
    real(dp) :: totalStress(3,3)

    !> Stress tensors for determinants if using TI-DFTB
    real(dp), allocatable :: mixedStress(:,:), tripletStress(:,:)

    ! Tagged writer
    type(TTaggedWriter) :: taggedWriter

    !> Container for the atomistic structure for poisson
    type(TPoissonStructure) :: poissStr

  #:if WITH_TRANSPORT
    !> Transport variables
    type(TTransPar) :: transpar
    type(TNEGFInfo) :: ginfo
  #:endif

    !> Transport interface (may be dummy placeholder, if built without transport)
    type(TNegfInt) :: negfInt

    !> Whether contact Hamiltonians are uploaded
    !> Synonym for G.F. calculation of density
    logical :: tUpload

    !> Whether contact Hamiltonians are computed
    logical :: tContCalc

    !> Whether Poisson solver is invoked
    logical :: tPoisson

    !> Whether recompute Poisson after every SCC
    logical :: tPoissonTwice

    !> Calculate terminal tunneling and current
    logical :: tTunn

    !> True if we use any part of Negf (green solver, landauer etc.)
    logical :: tNegf

    !> Whether local currents are computed
    logical :: tLocalCurrents

    !> True if LDOS is stored on separate files for k-points
    logical :: tWriteLDOS

    !> Labels for LDOS regions, if needed
    character(lc), allocatable :: regionLabelLDOS(:)

    !> True if Tunneling is stored on separate files
    logical :: writeTunn

    !> Holds spin-dependent electrochemical potentials of contacts
    !> This is because libNEGF is not spin-aware
    real(dp), allocatable :: mu(:,:)

    !> Variables for Transport NEGF/Poisson solver
    !> Tunneling, local DOS and current
    real(dp), allocatable :: tunneling(:,:), ldos(:,:), current(:,:)
    real(dp), allocatable :: leadCurrents(:)
    !> Array storing local (bond) currents
    real(dp), allocatable :: lCurrArray(:,:)

    !> Shell-resolved Potential shifts uploaded from contacts
    real(dp), allocatable :: shiftPerLUp(:,:)

    !> Orbital-resolved charges uploaded from contacts
    real(dp), allocatable :: chargeUp(:,:,:)

    !> Shell-resolved block potential shifts uploaded from contacts
    real(dp), allocatable :: shiftBlockUp(:,:,:,:)

    !> Block populations uploaded from contacts
    real(dp), allocatable :: blockUp(:,:,:,:)

    !> Details of energy interval for tunneling used in output
    real(dp) :: Emin, Emax, Estep

    !> Electrostatics type (either gammafunctional or poisson)
    integer :: electrostatics

    !> list of atoms in the central cell (or device region if transport)
    integer, allocatable :: iAtInCentralRegion(:)

    !> Correction for {O,N}-X bonds
    type(THalogenX), allocatable :: halogenXCorrection

    !> All of the excited energies actuall solved by Casida routines (if used)
    real(dp), allocatable :: energiesCasida(:)

    !> Type for determinant control in DFTB (Delta DFTB)
    type(TDftbDeterminants) :: deltaDftb

    !> Number of determinants in use in the calculation
    integer :: nDets

    !> Final SCC charges if multiple determinants being used
    real(dp), allocatable :: qDets(:,:,:,:)

    !> Final SCC block charges if multiple determinants being used
    real(dp), allocatable :: qBlockDets(:,:,:,:,:)

    !> Final density matrices if multiple determinants being used
    real(dp), allocatable :: deltaRhoDets(:,:)

    !> data type for REKS
    type(TReksCalc), allocatable :: reks

    !> atomic charge contribution in excited state
    real(dp), allocatable :: dQAtomEx(:)

  contains

    procedure :: initProgramVariables
    procedure :: setEquivalencyRelations
    procedure :: initializeCharges
    procedure :: initializeReferenceCharges
    procedure :: setNElectrons
  #:if WITH_TRANSPORT
    procedure :: initTransport
    procedure :: initTransportArrays
  #:endif
    procedure :: destructProgramVariables
  #:if WITH_SOCKETS
    procedure :: initSocket
  #:endif
    procedure :: initOutputFiles
    procedure :: initArrays
    procedure :: initDetArrays
    procedure :: allocateDenseMatrices
  #:if WITH_SCALAPACK
    procedure :: initScalapack
  #:endif
    procedure :: getDenseDescCommon
    procedure :: ensureRangeSeparatedReqs
    procedure :: initRangeSeparated
    procedure :: initPlumed

  end type TDftbPlusMain

contains


  !> Initializes the variables in the module based on the parsed input
  subroutine initProgramVariables(this, input, env)

    !> Instance
    class(TDftbPlusMain), intent(inout), target :: this

    !> Holds the parsed input data.
    type(TInputData), intent(inout), target :: input

    !> Environment settings
    type(TEnvironment), intent(inout) :: env

    ! Mixer related local variables
    integer :: nGeneration
    real(dp) :: mixParam

    !> mixer number
    integer :: iMixer

    !> simple mixer (if used)
    type(TSimpleMixer), allocatable :: pSimpleMixer

    !> Anderson mixer (if used)
    type(TAndersonMixer), allocatable :: pAndersonMixer

    !> Broyden mixer (if used)
    type(TBroydenMixer), allocatable :: pBroydenMixer

    !> DIIS mixer (if used)
    type(TDIISMixer), allocatable :: pDIISMixer

    ! Geometry optimizer related local variables

    !> Conjugate gradient driver
    type(TConjGrad), allocatable :: pConjGrad

    !> Steepest descent driver
    type(TSteepDesc), allocatable :: pSteepDesc

    !> Conjugate gradient driver
    type(TConjGrad), allocatable :: pConjGradLat

    !> Steepest descent driver
    type(TSteepDesc), allocatable :: pSteepDescLat

    !> gradient DIIS driver
    type(TDIIS), allocatable :: pDIIS

    !> lBFGS driver for geometry optimisation
    type(TLbfgs), allocatable :: pLbfgs

    !> lBFGS driver for lattice optimisation
    type(TLbfgs), allocatable :: pLbfgsLat

    !> FIRE driver for geometry optimisation
    type(TFire), allocatable :: pFire

    !> FIRE driver for lattice optimisation
    type(TFire), allocatable :: pFireLat

    ! MD related local variables
    type(TThermostat), allocatable :: pThermostat
    type(TDummyThermostat), allocatable :: pDummyTherm
    type(TAndersenThermostat), allocatable :: pAndersenTherm
    type(TBerendsenThermostat), allocatable :: pBerendsenTherm
    type(TNHCThermostat), allocatable :: pNHCTherm

    type(TVelocityVerlet), allocatable :: pVelocityVerlet
    type(TTempProfile), pointer :: pTempProfile

    real(dp), allocatable :: tmpCoords(:), tmpWeight(:), tmp3Coords(:,:)

    type(TRanlux), allocatable :: randomInit, randomThermostat
    integer :: iSeed

    integer :: ind, ii, jj, kk, iAt, iSp, iSh, iOrb

    ! Dispersion
    type(TDispSlaKirk), allocatable :: slaKirk
    type(TDispUFF), allocatable :: uff
  #:if WITH_DFTD3
    type(TDispDftD3), allocatable :: dftd3
  #:endif
    type(TSimpleDftD3), allocatable :: sdftd3
    type(TDispDftD4), allocatable :: dftd4
  #:if WITH_MBD
    type(TDispMbd), allocatable :: mbd
  #:endif

    ! H5 correction
    type(TH5Corr), allocatable :: pH5Correction
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

    ! Damped interactions
    logical, allocatable, target :: tDampedShort(:)
    type(TThirdOrderInp) :: thirdInp

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
    integer :: nInitNeighbour = 40

    !> Spin loop index
    integer :: iSpin

    !> Nr. of buffered Cholesky-decompositions
    integer :: nBufferedCholesky

    character(sc), allocatable :: shellNamesTmp(:)
    logical :: tRequireDerivator

    logical :: tInitialized

    !> Format for two using exponential notation values with units
    character(len=*), parameter :: format2Ue = "(A, ':', T30, E14.6, 1X, A, T50, E14.6, 1X, A)"

    @:ASSERT(input%tInitialized)

    write(stdOut, "(/, A)") "Starting initialization..."
    write(stdOut, "(A80)") repeat("-", 80)

    call env%initGlobalTimer(input%ctrl%timingLevel, "DFTB+ running times", stdOut)
    call env%globalTimer%startTimer(globalTimers%globalInit)

    ! Basic variables
    this%hamiltonianType = input%ctrl%hamiltonian
    this%tSccCalc = input%ctrl%tScc
    if (allocated(input%ctrl%dftbUInp)) then
      allocate(this%dftbU)
      call TDftbU_init(this%dftbU, input%ctrl%dftbUInp)
    end if
    this%tSpin = input%ctrl%tSpin
    this%nSpin = 1
    if (input%ctrl%reksInp%reksAlg /= reksTypes%noReks) then
      ! REKS follows spin-restricted open-shell scheme so nSpin should be two in the main code, but
      ! some variables such as this%qOutput should be treated in a restricted scheme. Here nSpin is
      ! set to one and changes to two later in the initialization.
      allocate(this%reks)
    else if (this%tSpin) then
      ! unrestricted spin polarisation
      this%nSpin = 2
    end if
    nIndepHam = this%nSpin

    this%tSpinSharedEf = input%ctrl%tSpinSharedEf
    this%tSpinOrbit = input%ctrl%tSpinOrbit
    this%tDualSpinOrbit = input%ctrl%tDualSpinOrbit
    this%t2Component = input%ctrl%t2Component
    this%isRangeSep = allocated(input%ctrl%rangeSepInp)

    if (this%t2Component) then
      this%nSpin = 4
      nIndepHam = 1
    end if

    if (this%nSpin /= 2 .and. this%tSpinSharedEf) then
      call error("Colinear spin polarization required for shared Ef over spin channels")
    end if

    this%nAtom = input%geom%nAtom
    this%nType = input%geom%nSpecies
    select case(this%hamiltonianType)
    case default
      call error("Invalid Hamiltonian")
    case(hamiltonianTypes%dftb)
      this%orb = input%slako%orb
    case(hamiltonianTypes%xtb)
      ! TODO
      call error("xTB calculation currently not supported")
    end select
    this%nOrb = this%orb%nOrb
    this%tPeriodic = input%geom%tPeriodic
    this%tHelical = input%geom%tHelical

    ! start by assuming stress can be calculated if periodic
    this%tStress = this%tPeriodic

    ! Brillouin zone sampling
    if (this%tPeriodic .or. this%tHelical) then
      this%nKPoint = input%ctrl%nKPoint
      allocate(this%kPoint(size(input%ctrl%KPoint,dim=1), this%nKPoint))
      allocate(this%kWeight(this%nKPoint))
      this%kPoint(:,:) = input%ctrl%KPoint
      if (sum(input%ctrl%kWeight(:)) < epsilon(1.0_dp)) then
        call error("Sum of k-point weights should be greater than zero!")
      end if
      this%kWeight(:) = input%ctrl%kWeight / sum(input%ctrl%kWeight)
    else
      this%nKPoint = 1
      allocate(this%kPoint(3, this%nKPoint))
      allocate(this%kWeight(this%nKPoint))
      this%kPoint(:,1) = 0.0_dp
      this%kWeight(1) = 1.0_dp
    end if

    this%tRealHS = .true.
    if (this%tPeriodic .or. this%tHelical) then
      if ( size(this%kPoint,dim=2) == 1 .and. all(this%kPoint(:, 1) == 0.0_dp)) then
        this%tRealHS = .true.
      else
        this%tRealHS = .false.
      end if
    end if

  #:if WITH_MPI

    if (input%ctrl%parallelOpts%nGroup > nIndepHam * this%nKPoint) then
      write(stdOut, *)"Parallel groups only relevant for tasks split over sufficent spins and/or&
          & k-points"
      write(tmpStr,"('Nr. groups:',I4,', Nr. indepdendent spins times k-points:',I4)")&
          & input%ctrl%parallelOpts%nGroup, nIndepHam * this%nKPoint
      call error(trim(tmpStr))
    end if

    call env%initMpi(input%ctrl%parallelOpts%nGroup)
  #:endif


  #:if WITH_SCALAPACK
    call this%initScalapack(input%ctrl%parallelOpts%blacsOpts, this%nAtom, this%nOrb,&
        & this%t2Component, env)
  #:endif
    call TParallelKS_init(this%parallelKS, env, this%nKPoint, nIndepHam)

    this%sccTol = input%ctrl%sccTol
    this%tShowFoldedCoord = input%ctrl%tShowFoldedCoord
    if (this%tShowFoldedCoord .and. .not. (this%tPeriodic .or. this%tHelical)) then
      call error("Folding coordinates back into the central cell is meaningless for molecular&
          & boundary conditions!")
    end if
    this%tFracCoord = input%geom%tFracCoord

    ! no point if not SCC
    this%isSccConvRequired = (input%ctrl%isSccConvRequired .and. this%tSccCalc)

    if (this%tSccCalc) then
      this%maxSccIter = input%ctrl%maxIter
    else
      this%maxSccIter = 1
    end if
    if (this%maxSccIter < 1) then
      call error("SCC iterations must be larger than 0")
    end if
    if (this%tSccCalc) then
      if (allocated(input%ctrl%elecDynInp)) then
        if (input%ctrl%elecDynInp%tReadRestart .and. .not.input%ctrl%elecDynInp%tPopulations) then
          this%maxSccIter = 0
          this%isSccConvRequired = .false.
          this%tRestartNoSC = .true.
        end if
      end if
    end if

    this%tWriteHS = input%ctrl%tWriteHS
    this%tWriteRealHS = input%ctrl%tWriteRealHS

    if (this%tPeriodic) then
      this%tLatticeChanged = .true.
      allocate(this%latVec(3, 3))
      @:ASSERT(all(shape(input%geom%latVecs) == shape(this%latVec)))
      this%latVec(:,:) = input%geom%latVecs(:,:)
      allocate(this%origin(3))
      this%origin(:) = input%geom%origin
      allocate(this%recVec(3, 3))
      allocate(this%invLatVec(3, 3))
      this%invLatVec = this%latVec(:,:)
      call matinv(this%invLatVec)
      this%invLatVec = reshape(this%invLatVec, (/3, 3/), order=(/2, 1/))
      this%recVec = 2.0_dp * pi * this%invLatVec
      this%CellVol = abs(determinant33(this%latVec))
      this%recCellVol = abs(determinant33(this%recVec))
    else if (this%tHelical) then
      this%origin = input%geom%origin
      this%latVec = input%geom%latVecs(:,:)
      allocate(this%recVec(1, 1))
      this%recVec = 1.0_dp / this%latVec(1,1)
      allocate(this%invLatVec(0, 0))
    else
      allocate(this%latVec(0, 0))
      allocate(this%origin(0))
      allocate(this%recVec(0, 0))
      allocate(this%invLatVec(0, 0))
      this%CellVol = 0.0_dp
      this%recCellVol = 0.0_dp
      this%tLatticeChanged = .false.
    end if

    select case(this%hamiltonianType)
    case default
      call error("Invalid Hamiltonian")
    case(hamiltonianTypes%dftb)
      ! Slater-Koster tables
      this%skHamCont = input%slako%skHamCont
      this%skOverCont = input%slako%skOverCont
      this%pRepCont = input%slako%repCont

      allocate(this%atomEigVal(this%orb%mShell, this%nType))
      @:ASSERT(size(input%slako%skSelf, dim=1) == this%orb%mShell)
      @:ASSERT(size(input%slako%skSelf, dim=2) == size(this%atomEigVal, dim=2))
      this%atomEigVal(:,:) = input%slako%skSelf(1:this%orb%mShell, :)

      @:ASSERT(size(input%slako%skOcc, dim=1) >= this%orb%mShell)
      @:ASSERT(size(input%slako%mass) == this%nType)
      allocate(this%speciesMass(this%nType))
      this%speciesMass(:) = input%slako%mass(:)
    case(hamiltonianTypes%xtb)
      ! TODO
      call error("xTB calculation currently not supported")
    end select

    ! Spin W's !'
    if (allocated(input%ctrl%spinW)) then
      allocate(this%spinW(this%orb%mShell, this%orb%mShell, this%nType))
      this%spinW(:,:,:) = 0.0_dp
      do iSp = 1, this%nType
        do jj = 1, this%orb%nShell(iSp)
          do kk = 1, this%orb%nShell(iSp)
            this%spinW(jj, kk, iSp) = input%ctrl%spinW(jj, kk, iSp)
          end do
        end do
      end do
    end if

    if (this%tSpinOrbit) then
      allocate(this%xi(this%orb%mShell,this%nType))
      this%xi(:,:) = 0.0_dp
      do iSp=1,this%nType
        do jj=1, this%orb%nShell(iSp)
          this%xi(jj,iSp)=input%ctrl%xi(jj,iSp)
        end do
      end do
    end if

    ! on-site corrections
    if (allocated(input%ctrl%onSiteElements)) then
      allocate(this%onSiteElements(this%orb%mShell, this%orb%mShell, 2, this%nType))
      this%onSiteElements(:,:,:,:) = input%ctrl%onSiteElements(:,:,:,:)
    end if

    this%tMixBlockCharges = allocated(this%dftbU) .or. allocated(this%onSiteElements)

    select case(this%hamiltonianType)
    case default
      call error("Invalid Hamiltonian")
    case(hamiltonianTypes%dftb)
      ! Cut-offs for SlaKo, repulsive
      this%cutOff%skCutOff = max(getCutOff(this%skHamCont), getCutOff(this%skOverCont))
      this%cutOff%repCutOff = getCutOff(this%pRepCont)
      this%cutOff%mCutOff = maxval([this%cutOff%skCutOff, this%cutOff%repCutOff])
    case(hamiltonianTypes%xtb)
      ! TODO
      call error("xTB calculation currently not supported")
    end select

    ! Get species names and output file
    this%geoOutFile = input%ctrl%outFile
    allocate(this%speciesName(size(input%geom%speciesNames)))
    this%speciesName(:) = input%geom%speciesNames(:)

    do iSp = 1, this%nType
      do jj = iSp+1, this%nType
        if (this%speciesName(iSp) == this%speciesName(jj)) then
          write (tmpStr,"('Duplicate identical species labels in the geometry: ',A)")&
              & this%speciesName(iSp)
          call error(tmpStr)
        end if
      end do
    end do

    ! Initialise the SCC module (the two copies of the Hubbard Us are rather
    ! artifical, since the copy for the main program is only used for dumping
    ! into the tagged format for autotest)
    allocate(this%hubbU(this%orb%mShell, this%nType))

    select case(this%hamiltonianType)
    case default
      call error("Invalid Hamiltonian")
    case(hamiltonianTypes%dftb)
      @:ASSERT(size(input%slako%skHubbU, dim=1) >= this%orb%mShell)
      @:ASSERT(size(input%slako%skHubbU, dim=2) == this%nType)
      this%hubbU(:,:) = input%slako%skHubbU(1:this%orb%mShell, :)
    case(hamiltonianTypes%xtb)
      ! TODO
      call error("xTB calculation currently not supported")
    end select

    if (allocated(input%ctrl%hubbU)) then
      where (input%ctrl%hubbU > 0.0_dp)
        this%hubbU = input%ctrl%hubbU
      end where
    end if

    this%tPoisson = input%ctrl%tPoisson

    if (this%tSccCalc) then
      allocate(sccInp)
      allocate(this%sccCalc)
      sccInp%orb => this%orb

      sccInp%hasExternalShifts = this%tPoisson

      if (this%tPeriodic) then
        sccInp%latVecs = this%latVec
        sccInp%recVecs = this%recVec
        sccInp%volume = this%CellVol
      else if (this%tHelical) then
        call error("Scc calculations not currently supported for helical boundary conditions")
      end if
      sccInp%hubbU = this%hubbU
      allocate(tDampedShort(this%nType))
      if (input%ctrl%tDampH) then
        tDampedShort = (this%speciesMass < 3.5_dp * amu__au)
        !tDampedShort(:) = (this%speciesName == "H" .or. this%speciesName == "h")
      else
        tDampedShort(:) = .false.
      end if
      sccInp%tDampedShort = tDampedShort
      sccInp%dampExp = input%ctrl%dampExp

      ! H5 correction
      if (input%ctrl%h5SwitchedOn) then
        if (.not. any(this%speciesMass < 3.5_dp * amu__au)) then
          call error("H5 correction used without H atoms present")
        end if
        if (any(tDampedShort)) then
          call error("H5 correction is not compatible with X-H damping")
        end if
        allocate(pH5Correction)
        call H5Corr_init(pH5Correction, this%speciesName, input%ctrl%h5RScale, input%ctrl%h5WScale,&
            & input%ctrl%h5ElementPara)
        sccInp%h5Correction = pH5Correction
      end if

      this%nExtChrg = input%ctrl%nExtChrg
      this%tExtChrg = (this%nExtChrg > 0)
      if (this%tExtChrg) then
        if (.not.this%tSccCalc) then
          call error("External charges can only be used in an SCC calculation")
        end if
        this%tStress = .false.
        @:ASSERT(size(input%ctrl%extChrg, dim=1) == 4)
        @:ASSERT(size(input%ctrl%extChrg, dim=2) == this%nExtChrg)
        sccInp%extCharges = input%ctrl%extChrg
        if (allocated(input%ctrl%extChrgBlurWidth)) then
          sccInp%blurWidths = input%ctrl%extChrgblurWidth
          if (any(sccInp%blurWidths < 0.0_dp)) then
            call error("Gaussian blur widths for charges may not be negative")
          end if
        end if
      end if
      if (allocated(input%ctrl%chrgConstr)) then
        @:ASSERT(all(shape(input%ctrl%chrgConstr) == (/ this%nAtom, 2 /)))
        if (any(abs(input%ctrl%chrgConstr(:,2)) > epsilon(1.0_dp))) then
          sccInp%chrgConstraints = input%ctrl%chrgConstr
        end if
      end if

      if (allocated(input%ctrl%thirdOrderOn)) then
        @:ASSERT(this%tSccCalc)
        @:ASSERT(all(shape(input%ctrl%thirdOrderOn) == (/ this%nAtom, 2 /)))
        sccInp%thirdOrderOn = input%ctrl%thirdOrderOn
      end if

      sccInp%coulombInput%ewaldAlpha = input%ctrl%ewaldAlpha
      sccInp%coulombInput%tolEwald = input%ctrl%tolEwald
      call initialize(this%sccCalc, env, sccInp)
      deallocate(sccInp)

      ! Longest cut-off including the softening part of gamma
      this%cutOff%mCutOff = max(this%cutOff%mCutOff, this%sccCalc%getCutOff())

      if (input%ctrl%t3rd .and. input%ctrl%tShellResolved) then
        call error("Onsite third order DFTB only compatible with shell non-resolved SCC")
      end if

      ! Initialize full 3rd order module
      this%t3rd = input%ctrl%t3rd
      this%t3rdFull = input%ctrl%t3rdFull
      if (this%t3rdFull) then
        @:ASSERT(this%tSccCalc)
        thirdInp%orb => this%orb
        thirdInp%hubbUs = this%hubbU
        thirdInp%hubbUDerivs = input%ctrl%hubDerivs
        allocate(thirdInp%damped(this%nType))
        thirdInp%damped(:) = tDampedShort
        thirdInp%dampExp = input%ctrl%dampExp
        thirdInp%shellResolved = input%ctrl%tShellResolved
        allocate(this%thirdOrd)
        call ThirdOrder_init(this%thirdOrd, thirdInp)
        this%cutOff%mCutOff = max(this%cutOff%mCutOff, this%thirdOrd%getCutOff())
      end if
    end if

    ! Initial coordinates
    allocate(this%coord0(3, this%nAtom))
    @:ASSERT(all(shape(this%coord0) == shape(input%geom%coords)))
    this%coord0(:,:) = input%geom%coords(:,:)

    this%tCoordsChanged = .true.

    allocate(this%species0(this%nAtom))
    @:ASSERT(all(shape(this%species0) == shape(input%geom%species)))
    this%species0(:) = input%geom%species(:)

  #:block DEBUG_CODE
    call inputCoherenceCheck(env, this%hamiltonianType, this%nSpin, this%nAtom, this%coord0,&
       & this%species0, this%speciesName, this%tSccCalc, this%tPeriodic, this%tFracCoord,&
       & this%latVec, this%origin)
  #:endblock DEBUG_CODE

    if (input%ctrl%tHalogenX) then
      if (.not. (this%t3rd .or. this%t3rdFull)) then
        call error("Halogen correction only fitted for 3rd order models")
      end if
      if (this%tPeriodic) then
        call error("Halogen correction was not fitted in periodic systems in original paper")
      end if
      allocate(this%halogenXCorrection)
      call THalogenX_init(this%halogenXCorrection, this%species0, this%speciesName)
    end if

    allocate(this%referenceN0(this%orb%mShell, this%nType))
    allocate(this%mass(this%nAtom))
    this%mass = this%speciesMass(this%species0)
    if (allocated(input%ctrl%masses)) then
      @:ASSERT(size(input%ctrl%masses) == this%nAtom)
      where (input%ctrl%masses >= 0.0_dp)
        this%mass = input%ctrl%masses
      end where
    end if

    if (this%tPeriodic) then
      ! Make some guess for the nr. of all interacting atoms
      this%nAllAtom = int((real(this%nAtom, dp)**(1.0_dp/3.0_dp) + 3.0_dp)**3)
    else if (this%tHelical) then
      ! 1D system, so much lower number of initial interactions
      this%nAllAtom = this%nAtom + 3
      if (size(this%latVec,dim=1)==3) then
        this%nAllAtom = this%nAllAtom * nint(this%latVec(3,1))
      end if
    else
      this%nAllAtom = this%nAtom
    end if
    allocate(this%coord(3, this%nAllAtom))
    allocate(this%species(this%nAllAtom))
    allocate(this%img2CentCell(this%nAllAtom))
    allocate(this%iCellVec(this%nAllAtom))

    ! Intialize Hamilton and overlap
    this%tImHam = this%tDualSpinOrbit .or. (this%tSpinOrbit .and. allocated(this%dftbU))
    if (this%tSccCalc .and. .not.allocated(this%reks)) then
      allocate(this%chargePerShell(this%orb%mShell,this%nAtom,this%nSpin))
    else
      allocate(this%chargePerShell(0,0,0))
    end if
    if (.not.allocated(this%reks)) then
      allocate(this%ham(0, this%nSpin))
    end if
    if (this%tImHam) then
      allocate(this%iHam(0, this%nSpin))
    end if
    allocate(this%over(0))
    allocate(this%iSparseStart(0, this%nAtom))

    if (this%nSpin == 4) then
      allocate(this%nEl(1))
      allocate(this%Ef(1))
    else
      allocate(this%nEl(this%nSpin))
      allocate(this%Ef(this%nSpin))
    end if

    this%iDistribFn = input%ctrl%iDistribFn
    this%tempElec = input%ctrl%tempElec

    this%tFixEf = input%ctrl%tFixEf
    if (allocated(input%ctrl%Ef)) then
      this%Ef(:) = input%ctrl%Ef
    else
      this%Ef(:) = 0.0_dp
    end if
    this%tSetFillingTemp = input%ctrl%tSetFillingTemp
    this%tFillKSep = input%ctrl%tFillKSep
    this%tempAtom = input%ctrl%tempAtom
    this%deltaT = input%ctrl%deltaT

    ! Orbital equivalency relations
    call this%setEquivalencyRelations()

    ! Initialize mixer
    ! (at the moment, the mixer does not need to know about the size of the vector to mix.)
    if (this%tSccCalc .and. .not. allocated(this%reks) .and. .not. this%tRestartNoSC) then
      allocate(this%pChrgMixer)
      iMixer = input%ctrl%iMixSwitch
      nGeneration = input%ctrl%iGenerations
      mixParam = input%ctrl%almix
      select case (iMixer)
      case (mixerTypes%simple)
        allocate(pSimplemixer)
        call init(pSimpleMixer, mixParam)
        call init(this%pChrgMixer, pSimpleMixer)
      case (mixerTypes%anderson)
        allocate(pAndersonMixer)
        if (input%ctrl%andersonNrDynMix > 0) then
          call init(pAndersonMixer, nGeneration, mixParam, input%ctrl%andersonInitMixing,&
              & input%ctrl%andersonDynMixParams, input%ctrl%andersonOmega0)
        else
          call init(pAndersonMixer, nGeneration, mixParam, input%ctrl%andersonInitMixing,&
              & omega0=input%ctrl%andersonOmega0)
        end if
        call init(this%pChrgMixer, pAndersonMixer)
      case (mixerTypes%broyden)
        allocate(pBroydenMixer)
        call init(pBroydenMixer, this%maxSccIter, mixParam, input%ctrl%broydenOmega0,&
            & input%ctrl%broydenMinWeight, input%ctrl%broydenMaxWeight, input%ctrl%broydenWeightFac)
        call init(this%pChrgMixer, pBroydenMixer)
      case(mixerTypes%diis)
        allocate(pDIISMixer)
        call init(pDIISMixer,nGeneration, mixParam, input%ctrl%tFromStart)
        call init(this%pChrgMixer, pDIISMixer)
      case default
        call error("Unknown charge mixer type.")
      end select
    end if

    ! initialise in cases where atoms move
    this%isGeoOpt = input%ctrl%isGeoOpt
    this%tCoordOpt = input%ctrl%tCoordOpt
    this%tLatOpt = (input%ctrl%tLatOpt .and. this%tPeriodic)
    if (this%tLatOpt) then
      if (this%tExtChrg) then
        ! Stop as not sure, what to do with the coordinates of the
        ! external charges, when the lattice changes.
        call error("External charges and lattice optimisation can not be used together.")
      end if
    end if
    if (this%tLatOpt) then
      this%tLatOptFixAng = input%ctrl%tLatOptFixAng
      this%tLatOptFixLen = input%ctrl%tLatOptFixLen
      this%tLatOptIsotropic = input%ctrl%tLatOptIsotropic
      if (this%tLatOptFixAng .or. any(this%tLatOptFixLen) .or. this%tLatOptIsotropic) then
        this%origLatVec(:,:) = this%latVec(:,:)
        do ii = 1, 3
           this%normOrigLatVec(:,ii) = this%origLatVec(:,ii) / sqrt(sum(this%origLatVec(:,ii)**2))
        end do
      end if
    end if
    this%extPressure = input%ctrl%pressure
    this%tBarostat = input%ctrl%tBarostat
    this%BarostatStrength = input%ctrl%BarostatStrength

  #:if WITH_SOCKETS
    this%tSocket = allocated(input%ctrl%socketInput)
    if (this%tSocket) then
      input%ctrl%socketInput%nAtom = this%nAtom
      call this%initSocket(env, input%ctrl%socketInput)
      this%tForces = .true.
      this%isGeoOpt = .false.
      this%tMD = .false.
    end if
  #:else
    this%tSocket = .false.
  #:endif

    this%tAppendGeo = input%ctrl%tAppendGeo
    this%isSccConvRequired = input%ctrl%isSccConvRequired
    this%tMD = input%ctrl%tMD
    this%tDerivs = input%ctrl%tDerivs
    this%tPrintMulliken = input%ctrl%tPrintMulliken
    this%tEField = input%ctrl%tEfield
    this%tExtField = this%tEField
    this%tMulliken = input%ctrl%tMulliken .or. this%tPrintMulliken .or. this%tExtField .or.&
        & this%tFixEf .or. this%isRangeSep
    this%tAtomicEnergy = input%ctrl%tAtomicEnergy
    this%tPrintEigVecs = input%ctrl%tPrintEigVecs
    this%tPrintEigVecsTxt = input%ctrl%tPrintEigVecsTxt

    this%tPrintForces = input%ctrl%tPrintForces
    this%tForces = input%ctrl%tForces .or. this%tPrintForces
    this%isLinResp = input%ctrl%lrespini%tInit
    if (this%isLinResp) then
      allocate(this%linearResponse)
    end if

    select case(this%hamiltonianType)
    case default
      call error("Invalid Hamiltonian")
    case(hamiltonianTypes%dftb)
      this%referenceN0(:,:) = input%slako%skOcc(1:this%orb%mShell, :)
    case(hamiltonianTypes%xtb)
      ! TODO
      call error("xTB calculation currently not supported")
    end select

    this%nrChrg = input%ctrl%nrChrg
    this%nrSpinPol = input%ctrl%nrSpinPol

    if (this%isLinResp) then
      allocate(this%dQAtomEx(this%nAtom))
      this%dQAtomEx(:) = 0.0_dp
    end if

    call this%initializeReferenceCharges(input%ctrl%customOccAtoms, input%ctrl%customOccFillings)
    call this%setNElectrons()

    ! DFTB related variables if multiple determinants are used
    call TDftbDeterminants_init(this%deltaDftb, input%ctrl%isNonAufbau, input%ctrl%isSpinPurify,&
        & input%ctrl%isGroundGuess, this%nEl, this%dftbEnergy)

    if (this%tForces) then
      this%tCasidaForces = input%ctrl%tCasidaForces
    else
      this%tCasidaForces = .false.
    end if
    if (this%tSccCalc) then
      this%forceType = input%ctrl%forceType
    else
      if (input%ctrl%forceType /= forceTypes%orig) then
        call error("Invalid force evaluation method for non-SCC calculations.")
      end if
    end if
    if (this%forceType == forceTypes%dynamicT0 .and. this%tempElec > minTemp) then
       call error("This ForceEvaluation method requires the electron temperature to be zero")
     end if

     tRequireDerivator = this%tForces
     if (.not. tRequireDerivator .and. allocated(input%ctrl%elecDynInp)) then
       tRequireDerivator = input%ctrl%elecDynInp%tIons
     end if
     if (tRequireDerivator) then
      select case(input%ctrl%iDerivMethod)
      case (1)
        ! set step size from input
        if (input%ctrl%deriv1stDelta < epsilon(1.0_dp)) then
          write(tmpStr, "(A,E12.4)") 'Too small value for finite difference step :',&
              & input%ctrl%deriv1stDelta
          call error(tmpStr)
        end if
        call NonSccDiff_init(this%nonSccDeriv, diffTypes%finiteDiff, input%ctrl%deriv1stDelta)
      case (2)
        call NonSccDiff_init(this%nonSccDeriv, diffTypes%richardson)
      end select
    end if

    call this%getDenseDescCommon() !this%orb, this%nAtom, this%t2Component, this%denseDesc)

    call ensureSolverCompatibility(input%ctrl%solver%iSolver, this%tSpin, this%kPoint,&
        & input%ctrl%parallelOpts, nIndepHam, this%tempElec)

    if (this%tRealHS) then
      nBufferedCholesky = 1
    else
      nBufferedCholesky = this%parallelKS%nLocalKS
    end if
    call TElectronicSolver_init(this%electronicSolver, input%ctrl%solver%iSolver, nBufferedCholesky)

    if (this%electronicSolver%isElsiSolver) then
      @:ASSERT(this%parallelKS%nLocalKS == 1)

      if (input%ctrl%parallelOpts%nGroup /= nIndepHam * this%nKPoint) then
        if (this%nSpin == 2) then
          write(tmpStr, "(A,I0,A,I0,A)")"ELSI solvers require as many groups as spin and k-point&
              & combinations. There are ", nIndepHam * this%nKPoint, " spin times k-point&
              & combinations and ", input%ctrl%parallelOpts%nGroup, " groups"
        else
          write(tmpStr, "(A,I0,A,I0,A)")"ELSI solvers require as many groups as k-points. There&
              & are ", nIndepHam * this%nKPoint, " k-points and ", input%ctrl%parallelOpts%nGroup,&
              & " groups"
        end if
        call error(tmpStr)
      end if

      #:if WITH_OMP
        if (omp_get_max_threads() > 1) then
          call error("The ELSI-solvers should not be run with multiple threads. Set the&
              & environment variable OMP_NUM_THREADS to 1 in order to disable multi-threading.")
        end if
      #:endif

      ! Would be using the ELSI matrix writing mechanism, so set this as always false
        this%tWriteHS = .false.

      call TElsiSolver_init(this%electronicSolver%elsi, input%ctrl%solver%elsi, env,&
          & this%denseDesc%fullSize, this%nEl, this%iDistribFn, this%nSpin,&
          & this%parallelKS%localKS(2, 1), this%nKPoint, this%parallelKS%localKS(1, 1),&
          & this%kWeight(this%parallelKS%localKS(1, 1)), input%ctrl%tWriteHS,&
          & this%electronicSolver%providesElectronEntropy)

    end if

    if (this%deltaDftb%isNonAufbau .and. .not.this%electronicSolver%providesEigenvals) then
      call error("Eigensolver that calculates eigenvalues is required for Delta DFTB")
    end if

    if (allocated(this%reks)) then
      this%electronicSolver%providesElectronEntropy = .false.
    end if

    if (this%forceType /= forceTypes%orig .and. .not. this%electronicSolver%providesEigenvals) then
      call error("Alternative force evaluation methods are not supported by this electronic&
          & solver.")
    end if

  #:if WITH_TRANSPORT
    ! whether tunneling is computed
    this%tTunn = input%ginfo%tundos%defined
    ! whether local currents are computed
    this%tLocalCurrents = input%ginfo%greendens%doLocalCurr

    ! Do we use any part of negf (solver, tunnelling etc.)?
    this%tNegf = (this%electronicSolver%iSolver == electronicSolverTypes%GF) .or. this%tTunn .or.&
        & this%tLocalCurrents

  #:else

    this%tTunn = .false.
    this%tNegf = .false.

  #:endif

    ! temporary disables for various issues with NEGF
    if (this%tNegf) then
      if (this%nSpin > 2) then
        call error("Non-collinear spin polarization disabled for transport calculations at the&
            & moment.")
      end if
      if (this%tExtChrg) then
        call error("External charges temporarily disabled for transport calculations&
            & (electrostatic gates are available).")
      end if
    #:if WITH_TRANSPORT
      if (this%isRangeSep .and. this%transpar%nCont > 0) then
        call error("Range separated calculations do not work with transport calculations yet")
      end if
    #:endif
    end if


    ! requires stress to already be possible and it being a periodic calculation
    ! with forces
    this%tStress = (this%tPeriodic .and. this%tForces .and. .not.this%tNegf .and. this%tStress)

    this%nMovedAtom = input%ctrl%nrMoved
    this%nMovedCoord = 3 * this%nMovedAtom

    if (input%ctrl%maxRun == -1) then
      this%nGeoSteps = huge(1) - 1
      ! Workaround:PGI 17.10 -> do i = 0, huge(1) executes 0 times
      ! this%nGeoSteps = huge(1)
    else
      this%nGeoSteps = input%ctrl%maxRun
    end if

    if (this%nMovedAtom > 0) then
      allocate(this%indMovedAtom(size(input%ctrl%indMovedAtom)))
      this%indMovedAtom(:) = input%ctrl%indMovedAtom(:)
    else
      allocate(this%indMovedAtom(0))
    end if

    allocate(this%pGeoCoordOpt)
    if (this%tCoordOpt) then
      allocate(tmpCoords(this%nMovedCoord))
      tmpCoords(1:this%nMovedCoord) = reshape(this%coord0(:, this%indMovedAtom),&
          & (/ this%nMovedCoord /))
      select case (input%ctrl%iGeoOpt)
      case(geoOptTypes%steepestDesc)
        allocate(tmpWeight(this%nMovedCoord))
        tmpWeight(1:this%nMovedCoord) = 0.5_dp * this%deltaT**2 /&
            & reshape(spread(this%mass(this%indMovedAtom), 1, 3), (/this%nMovedCoord/))
        allocate(pSteepDesc)
        call init(pSteepDesc, size(tmpCoords), input%ctrl%maxForce, input%ctrl%maxAtomDisp,&
            & tmpWeight )
        deallocate(tmpWeight)
        call init(this%pGeoCoordOpt, pSteepDesc)
      case (geoOptTypes%conjugateGrad)
        allocate(pConjGrad)
        call init(pConjGrad, size(tmpCoords), input%ctrl%maxForce, input%ctrl%maxAtomDisp)
        call init(this%pGeoCoordOpt, pConjGrad)
      case (geoOptTypes%diis)
        allocate(pDIIS)
        call init(pDIIS, size(tmpCoords), input%ctrl%maxForce, input%ctrl%deltaGeoOpt,&
            & input%ctrl%iGenGeoOpt)
        call init(this%pGeoCoordOpt, pDIIS)
      case (geoOptTypes%lbfgs)
        allocate(pLbfgs)
        call TLbfgs_init(pLbfgs, size(tmpCoords), input%ctrl%maxForce, tolSameDist,&
            & input%ctrl%maxAtomDisp, input%ctrl%lbfgsInp%memory, input%ctrl%lbfgsInp%isLineSearch,&
            & input%ctrl%lbfgsInp%isOldLS, input%ctrl%lbfgsInp%MaxQNStep)
        call init(this%pGeoCoordOpt, pLbfgs)
      case (geoOptTypes%fire)
        allocate(pFire)
        call TFire_init(pFire, size(tmpCoords), input%ctrl%maxForce, input%ctrl%deltaT)
        call init(this%pGeoCoordOpt, pFire)
      end select
      call reset(this%pGeoCoordOpt, tmpCoords)
    end if

    allocate(this%pGeoLatOpt)
    if (this%tLatOpt) then
      select case (input%ctrl%iGeoOpt)
      case(geoOptTypes%steepestDesc)
        allocate(tmpWeight(9))
        tmpWeight = 1.0_dp
        allocate(pSteepDescLat)
        call init(pSteepDescLat, 9, input%ctrl%maxForce, input%ctrl%maxLatDisp, tmpWeight)
        deallocate(tmpWeight)
        call init(this%pGeoLatOpt, pSteepDescLat)
      case(geoOptTypes%conjugateGrad, geoOptTypes%diis) ! use CG lattice for both DIIS and CG
        allocate(pConjGradLat)
        call init(pConjGradLat, 9, input%ctrl%maxForce, input%ctrl%maxLatDisp)
        call init(this%pGeoLatOpt, pConjGradLat)
      case (geoOptTypes%LBFGS)
        allocate(pLbfgsLat)
        call TLbfgs_init(pLbfgsLat, 9, input%ctrl%maxForce, tolSameDist, input%ctrl%maxLatDisp,&
            & input%ctrl%lbfgsInp%memory, input%ctrl%lbfgsInp%isLineSearch,&
            & input%ctrl%lbfgsInp%isOldLS, input%ctrl%lbfgsInp%MaxQNStep)
        call init(this%pGeoLatOpt, pLbfgsLat)
      case (geoOptTypes%FIRE)
        allocate(pFireLat)
        call TFire_init(pFireLat, 9, input%ctrl%maxForce, input%ctrl%deltaT)
        call init(this%pGeoLatOpt, pFireLat)
      end select
      if (this%tLatOptIsotropic ) then
        ! optimization uses scaling factor of unit cell
        call reset(this%pGeoLatOpt,&
            & (/1.0_dp,0.0_dp,0.0_dp,0.0_dp,0.0_dp,0.0_dp,0.0_dp,0.0_dp,0.0_dp/))
      else if (this%tLatOptFixAng) then
        ! optimization uses scaling factor of lattice vectors
        call reset(this%pGeoLatOpt,&
            & (/1.0_dp,1.0_dp,1.0_dp,0.0_dp,0.0_dp,0.0_dp,0.0_dp,0.0_dp,0.0_dp/))
      else
        call reset(this%pGeoLatOpt, reshape(this%latVec, (/ 9 /)) )
      end if
    end if

    if (.not.(this%isGeoOpt.or.this%tMD.or.this%tSocket)) then
      this%nGeoSteps = 0
    end if

    if (this%tSocket .and. this%tHelical) then
      call error("The socket protocol does not understand helical geometries")
    end if

    ! Initialize constraints
    if (input%ctrl%nrConstr > 0) then
      allocate(this%conAtom(input%ctrl%nrConstr))
      allocate(this%conVec(3, input%ctrl%nrConstr))
      this%conAtom(:) = input%ctrl%conAtom
      this%conVec(:,:) = input%ctrl%conVec
      do ii = 1, input%ctrl%nrConstr
        this%conVec(:,ii) = this%conVec(:,ii) / sqrt(sum(this%conVec(:,ii)**2))
      end do
    end if

    ! Dispersion
    tHHRepulsion = .false.
    if (allocated(input%ctrl%dispInp)) then
      if (this%tHelical) then
        call error("Dispersion not currently supported for helical boundary conditions")
      end if
      if (allocated(input%ctrl%dispInp%slakirk)) then
        allocate(slaKirk)
        if (this%tPeriodic) then
          call DispSlaKirk_init(slaKirk, input%ctrl%dispInp%slakirk, this%latVec)
        else if (this%tHelical) then
          call error("Slater-Kirkwood incompatible with helical boundary conditions")
        else
          call DispSlaKirk_init(slaKirk, input%ctrl%dispInp%slakirk)
        end if
        call move_alloc(slaKirk, this%dispersion)

      elseif (allocated(input%ctrl%dispInp%uff)) then
        allocate(uff)
        if (this%tPeriodic) then
          call DispUff_init(uff, input%ctrl%dispInp%uff, this%nAtom, this%species0, this%latVec)
        else
          call DispUff_init(uff, input%ctrl%dispInp%uff, this%nAtom)
        end if
        call move_alloc(uff, this%dispersion)

    #:if WITH_DFTD3
      elseif (allocated(input%ctrl%dispInp%dftd3)) then
        allocate(dftd3)
        tHHRepulsion = input%ctrl%dispInp%dftd3%hhrepulsion
        if (tHHRepulsion .and. .not. any(this%speciesMass < 3.5_dp * amu__au)) then
          call error("H-H repulsion correction used without H atoms present")
        end if
        if (this%tPeriodic) then
          call DispDftD3_init(dftd3, input%ctrl%dispInp%dftd3, this%nAtom, this%species0,&
              & this%speciesName, this%latVec)
        else
          call DispDftD3_init(dftd3, input%ctrl%dispInp%dftd3, this%nAtom, this%species0,&
              & this%speciesName)
        end if
        call move_alloc(dftd3, this%dispersion)
    #:endif
      else if (allocated(input%ctrl%dispInp%sdftd3)) then
        allocate(sdftd3)
        if (this%tPeriodic) then
          call init(sdftd3, input%ctrl%dispInp%sdftd3, this%nAtom, this%species0, this%speciesName,&
              & this%latVec)
        else
          call init(sdftd3, input%ctrl%dispInp%sdftd3, this%nAtom, this%species0, this%speciesName)
        end if
        call move_alloc(sdftd3, this%dispersion)
      else if (allocated(input%ctrl%dispInp%dftd4)) then
        allocate(dftd4)
        if (this%tPeriodic) then
          call init(dftd4, input%ctrl%dispInp%dftd4, this%nAtom, this%speciesName, this%latVec)
        else
          call init(dftd4, input%ctrl%dispInp%dftd4, this%nAtom, this%speciesName)
        end if
        call move_alloc(dftd4, this%dispersion)
    #:if WITH_MBD
      else if (allocated(input%ctrl%dispInp%mbd)) then
        if (this%isLinResp) then
          call error("MBD model not currently supported for Casida linear response")
        end if
        allocate (mbd)
        associate (inp => input%ctrl%dispInp%mbd)
          inp%calculate_forces = this%tForces
          inp%atom_types = this%speciesName(this%species0)
          inp%coords = this%coord0
          if (this%tPeriodic) then
            inp%lattice_vectors = this%latVec
          end if
          call TDispMbd_init(mbd, inp, input%geom, isPostHoc=.true.)
        end associate
        call mbd%checkError()
        call move_alloc(mbd, this%dispersion)
        if (input%ctrl%dispInp%mbd%method == 'ts' .and. this%tForces) then
          call warning("Forces for the TS-dispersion model are calculated by finite differences&
              & which may result in long gradient calculation times for large systems")
        end if
    #:endif
      end if
      this%cutOff%mCutOff = max(this%cutOff%mCutOff, this%dispersion%getRCutOff())
    end if

    if (allocated(input%ctrl%solvInp)) then
      if (allocated(input%ctrl%solvInp%GBInp)) then
        if (this%tPeriodic) then
          call createSolvationModel(this%solvation, input%ctrl%solvInp%GBInp, &
              & this%nAtom, this%species0, this%speciesName, this%latVec)
        else
          call createSolvationModel(this%solvation, input%ctrl%solvInp%GBInp, &
              & this%nAtom, this%species0, this%speciesName)
        end if
      else if (allocated(input%ctrl%solvInp%SASAInp)) then
        if (this%tPeriodic) then
          call createSolvationModel(this%solvation, input%ctrl%solvInp%SASAInp, &
              & this%nAtom, this%species0, this%speciesName, this%latVec)
        else
          call createSolvationModel(this%solvation, input%ctrl%solvInp%SASAInp, &
              & this%nAtom, this%species0, this%speciesName)
        end if
      end if
      if (.not.allocated(this%solvation)) then
        call error("Could not initialize solvation model!")
      end if
      this%cutOff%mCutOff = max(this%cutOff%mCutOff, this%solvation%getRCutOff())
    end if

    if (allocated(this%halogenXCorrection)) then
      this%cutOff%mCutOff = max(this%cutOff%mCutOff, this%halogenXCorrection%getRCutOff())
    end if

    if (input%ctrl%nrChrg == 0.0_dp .and. .not.(this%tPeriodic.or.this%tHelical) .and.&
        & this%tMulliken) then
      this%tDipole = .true.
    else
      this%tDipole = .false.
    end if

    if (this%tMulliken) then
      this%tNetAtomCharges = input%ctrl%tNetAtomCharges
      if (allocated(input%ctrl%cm5Input)) then
        allocate(this%cm5Cont)
        if (this%tPeriodic) then
          call TChargeModel5_init(this%cm5Cont, input%ctrl%cm5Input, input%geom%nAtom, &
              & input%geom%speciesNames, .false., input%geom%latVecs)
        else
          call TChargeModel5_init(this%cm5Cont, input%ctrl%cm5Input, input%geom%nAtom, &
              & input%geom%speciesNames, .false.)
        end if
        this%cutOff%mCutOff = max(this%cutOff%mCutOff, this%cm5Cont%getRCutOff())
      end if
    else
      this%tNetAtomCharges = .false.
    end if

    if (allocated(input%ctrl%elStatPotentialsInp)) then
      if (.not.this%tSccCalc) then
        call error("Electrostatic potentials only available for SCC calculations")
      end if
      allocate(this%esp)
      call TElStatPotentials_init(this%esp, input%ctrl%elStatPotentialsInp, this%tEField .or.&
          & this%tExtChrg)
    end if

    if (allocated(input%ctrl%pipekMezeyInp)) then
      allocate(this%pipekMezey)
      call initialise(this%pipekMezey, input%ctrl%pipekMezeyInp)
    end if
    this%tLocalise = allocated(this%pipekMezey)
    if (this%tLocalise .and. (this%nSpin > 2 .or. this%t2Component)) then
      call error("Localisation of electronic states currently unsupported for non-collinear and&
          & spin orbit calculations")
    end if

    if (this%isLinResp) then

      ! input checking for linear response
      if (.not. withArpack) then
        call error("This binary has been compiled without support for linear response&
            & calculations.")
      end if
      call ensureLinRespConditions(this%tSccCalc, this%t3rd .or. this%t3rdFull, this%tRealHS,&
          & this%tPeriodic, this%tCasidaForces, this%solvation, this%isRS_LinResp, this%nSpin,&
          & this%tSpin, this%tHelical, this%tSpinOrbit, allocated(this%dftbU), this%tempElec, input)

      ! Hubbard U and spin constants for excitations (W only needed for triplet/spin polarised)
      allocate(input%ctrl%lrespini%HubbardU(this%nType))
      allocate(input%ctrl%lrespini%spinW(this%nType))
      input%ctrl%lrespini%HubbardU = 0.0_dp
      input%ctrl%lrespini%spinW = 0.0_dp

      ! calculate linear response Gamma values from HOAO shell Hubbard U (non
      ! shell resolved)
      do iSp = 1, this%nType
        homoLoc = maxloc(this%atomEigVal(:this%orb%nShell(iSp), iSp),&
            & mask=this%referenceN0(:this%orb%nShell(iSp), iSp) > 0.0_dp)
        input%ctrl%lrespini%HubbardU(iSp) = this%hubbU(homoLoc(1), iSp)
      end do

      ! and atomic HOAO spin W value if needed
      input%ctrl%lrespini%spinW(:) = 0.0_dp
      select case(input%ctrl%lrespini%sym)
      case("S")
        ! Singlet case, no need for spin constants
      case("T","B"," ")
        ! triplet or spin-polarised
        do iSp = 1, this%nType
          homoLoc = maxloc(this%atomEigVal(:this%orb%nShell(iSp), iSp),&
              & mask=this%referenceN0(:this%orb%nShell(iSp), iSp) > 0.0_dp)
          input%ctrl%lrespini%spinW(iSp) = this%spinW(homoLoc(1), homoLoc(1), iSp)
        end do
      case default
        call error("Unknown excitation type requested")
      end select

      this%tPrintExcitedEigVecs = input%ctrl%lrespini%tPrintEigVecs
      this%tLinRespZVect = (input%ctrl%lrespini%tMulliken .or. this%tCasidaForces .or.&
          & input%ctrl%lrespini%tCoeffs .or. this%tPrintExcitedEigVecs .or.&
          & input%ctrl%lrespini%tWriteDensityMatrix)

      if (allocated(this%onSiteElements) .and. this%tLinRespZVect) then
        call error("Excited state property evaluation currently incompatible with onsite&
            & corrections")
      end if

      call LinResp_init(this%linearResponse, input%ctrl%lrespini, this%nAtom, this%nEl(1),&
          & this%onSiteElements)

    end if

    ! turn on if LinResp and RangSep turned on, no extra input required for now
    this%isRS_LinResp = this%isLinResp .and. this%isRangeSep

    ! ppRPA stuff
    if (allocated(input%ctrl%ppRPA)) then

      if (abs(input%ctrl%nrChrg - 2.0_dp) > elecTolMax) then
        call warning("Particle-particle RPA should be for a reference system with a charge of +2.")
      end if

    #:for VAR, ERR in [("this%tSpinOrbit","spin orbit coupling"), &
      & ("this%tSpin","spin polarised ground state"), ("this%t3rd","third order"),&
      & ("any(this%kPoint /= 0.0_dp)","non-gamma k-points"),&
      & ("this%tFixEf", "a fixed Fermi level"), ("this%tPoisson", "use of the Poisson solver")]
      if (${VAR}$) then
        call error("PP-RPA does not support ${ERR}$")
      end if
    #:endfor
    #:for VAR, ERR in [("tShellResolved","shell resolved hamiltonians"),&
      & ("h5SwitchedOn","H5"), ("tDampH","H damping")]
      if (input%ctrl%${VAR}$) then
        call error("PP-RPA does not support ${ERR}$")
      end if
    #:endfor
    #:for VAR, ERR in [("this%solvation","solvation"), ("this%dftbU","DFTB+U/pSIC"),&
      & ("this%onSiteElements","onsite corrections"), ("this%reks","REKS")]
      if (allocated(${VAR}$)) then
        call error("PP-RPA does not support ${ERR}$")
      end if
    #:endfor
    #:if WITH_TRANSPORT
      if (input%transpar%defined) then
        call error("PP-RPA does not support transport calculations")
      end if
    #:endif

      if (this%isGeoOpt .or. this%tMD .or. this%tSocket) then
        call warning ("Geometry optimisation with ppRPA is probably not what you want - forces in&
            & the (N-2) electron ground state system do not match the targeted system for the&
            & excited states")
      end if

      call move_alloc(input%ctrl%ppRPA, this%ppRPA)

    end if

    iSeed = input%ctrl%iSeed
    tRandomSeed = (iSeed < 1)
    ! Note: This routine may not be called multiple times. If you need further random generators,
    ! extend the routine and create them within this call.
    call createRandomGenerators(env, iSeed, randomInit, randomThermostat)

    call getRandom(randomInit, rTmp)
    this%runId = int(real(huge(this%runId) - 1, dp) * rTmp) + 1


    ! MD stuff
    if (this%tMD) then
      ! Create MD framework.
      allocate(this%pMDFrame)
      call init(this%pMDFrame, this%nMovedAtom, this%nAtom, input%ctrl%tMDstill)

      ! Create temperature profile, if thermostat is not the dummy one
      if (input%ctrl%iThermostat /= 0) then
        allocate(this%temperatureProfile)
        call TempProfile_init(this%temperatureProfile, input%ctrl%tempMethods,&
            & input%ctrl%tempSteps, input%ctrl%tempValues)
        pTempProfile => this%temperatureProfile
      else
        nullify(pTempProfile)
      end if

      ! Create thermostat
      allocate(pThermostat)
      select case (input%ctrl%iThermostat)
      case (0) ! No thermostat
        allocate(pDummyTherm)
        call init(pDummyTherm, this%tempAtom, this%mass(this%indMovedAtom), randomThermostat,&
            & this%pMDFrame)
        call init(pThermostat, pDummyTherm)
      case (1) ! Andersen thermostat
        allocate(pAndersenTherm)
        call init(pAndersenTherm, randomThermostat, this%mass(this%indMovedAtom), pTempProfile,&
            & input%ctrl%tRescale, input%ctrl%wvScale, this%pMDFrame)
        call init(pThermostat, pAndersenTherm)
      case (2) ! Berendsen thermostat
        allocate(pBerendsenTherm)
        call init(pBerendsenTherm, randomThermostat, this%mass(this%indMovedAtom), pTempProfile,&
            & input%ctrl%wvScale, this%pMDFrame)
        call init(pThermostat, pBerendsenTherm)
      case (3) ! Nose-Hoover-Chain thermostat
        allocate(pNHCTherm)
        if (input%ctrl%tInitNHC) then
          call init(pNHCTherm, randomThermostat, this%mass(this%indMovedAtom), pTempProfile,&
              & input%ctrl%wvScale, this%pMDFrame, input%ctrl%deltaT, input%ctrl%nh_npart,&
              & input%ctrl%nh_nys, input%ctrl%nh_nc, input%ctrl%xnose, input%ctrl%vnose,&
              & input%ctrl%gnose)
        else
          call init(pNHCTherm, randomThermostat, this%mass(this%indMovedAtom), pTempProfile,&
              & input%ctrl%wvScale, this%pMDFrame, input%ctrl%deltaT, input%ctrl%nh_npart,&
              & input%ctrl%nh_nys, input%ctrl%nh_nc)
        end if
        call init(pThermostat, pNHCTherm)
      end select

      ! Create MD integrator
      allocate(pVelocityVerlet)
      if (input%ctrl%tReadMDVelocities) then
        if (this%tBarostat) then
          call init(pVelocityVerlet, this%deltaT, this%coord0(:,this%indMovedAtom), pThermostat,&
              & input%ctrl%initialVelocities, this%BarostatStrength, this%extPressure,&
              & input%ctrl%tIsotropic)
        else
          call init(pVelocityVerlet, this%deltaT, this%coord0(:,this%indMovedAtom), pThermostat,&
              & input%ctrl%initialVelocities, .true., .false.)
        end if
      else
        if (this%tBarostat) then
          call init(pVelocityVerlet, this%deltaT, this%coord0(:,this%indMovedAtom), pThermostat,&
              & this%BarostatStrength, this%extPressure, input%ctrl%tIsotropic)
        else
          call init(pVelocityVerlet, this%deltaT, this%coord0(:,this%indMovedAtom), pThermostat,&
              & input%ctrl%initialVelocities, .false., .false.)
        end if
      end if
      allocate(this%pMDIntegrator)
      call init(this%pMDIntegrator, pVelocityVerlet)
    end if

    call this%initPlumed(env, input%ctrl%tPlumed, this%tMD, this%plumedCalc)

    ! Check for extended Born-Oppenheimer MD
    this%isXlbomd = allocated(input%ctrl%xlbomd)
    if (this%isXlbomd) then
      if (input%ctrl%iThermostat /= 0) then
        call error("XLBOMD does not work with thermostats yet")
      elseif (this%tBarostat) then
        call error("XLBOMD does not work with barostats yet")
      elseif (this%nSpin /= 1 .or. allocated(this%dftbU) .or. allocated(this%onSiteElements)) then
        call error("XLBOMD does not work for spin, DFTB+U or onsites yet")
      elseif (this%forceType /= forceTypes%dynamicT0 .and. this%forceType /=&
          & forceTypes%dynamicTFinite) then
        call error("Force evaluation method incompatible with XLBOMD")
      elseif (this%iDistribFn /= fillingTypes%Fermi) then
        call error("Filling function incompatible with XLBOMD")
      end if
      allocate(this%xlbomdIntegrator)
      call Xlbomd_init(this%xlbomdIntegrator, input%ctrl%xlbomd, this%nIneqOrb)
    end if

    this%minSccIter = getMinSccIters(this%tSccCalc, allocated(this%dftbU), this%nSpin)
    this%minSccIter = min(this%minSccIter, this%maxSccIter)
    if (this%isXlbomd) then
      call this%xlbomdIntegrator%setDefaultSCCParameters(this%minSccIter, this%maxSccIter,&
          & this%sccTol)
    end if

    if (this%tDerivs) then
      allocate(tmp3Coords(3,this%nMovedAtom))
      tmp3Coords = this%coord0(:,this%indMovedAtom)
      call create(this%derivDriver, tmp3Coords, input%ctrl%deriv2ndDelta)
      this%coord0(:,this%indMovedAtom) = tmp3Coords
      deallocate(tmp3Coords)
      this%nGeoSteps = 2 * 3 * this%nMovedAtom - 1
    end if

    if (this%tEField) then
      this%EFieldStrength = input%ctrl%EFieldStrength
      this%EfieldVector(:) = input%ctrl%EfieldVector(:)
      this%tTDEfield = input%ctrl%tTDEfield
      this%EfieldOmega = input%ctrl%EfieldOmega
      this%EfieldPhase = input%ctrl%EfieldPhase
      if (this%tTDEfield .and. .not. this%tMD) then
        call error ("Time dependent electric fields only possible for MD!")
      end if
      ! parser should catch all of these:
      @:ASSERT(.not.this%tTDEfield .or. this%tMD)
    else
      this%EFieldStrength = 0.0_dp
      this%EfieldVector(:) = 0.0_dp
      this%tTDEfield = .false.
      this%EfieldOmega = 0.0_dp
      this%EfieldPhase = 0
    end if

    this%tReadChrg = input%ctrl%tReadChrg
    if (this%tReadChrg .and. this%deltaDftb%isNonAufbau) then
      call error("Charge restart not currently supported for Delta DFTB")
    end if

    if (this%isRangeSep) then
      call this%ensureRangeSeparatedReqs(input%ctrl%tShellResolved, input%ctrl%rangeSepInp)
      call getRangeSeparatedCutoff(input%ctrl%rangeSepInp%cutoffRed, this%cutOff)
      call this%initRangeSeparated(this%nAtom, this%species0, this%hubbU, input%ctrl%rangeSepInp,&
          & this%tSpin, allocated(this%reks), this%rangeSep, this%deltaRhoIn, this%deltaRhoOut,&
          & this%deltaRhoDiff, this%deltaRhoInSqr, this%deltaRhoOutSqr, this%nMixElements)
    end if

    this%tReadShifts = input%ctrl%tReadShifts
    this%tWriteShifts = input%ctrl%tWriteShifts

    ! Both temporarily removed until debugged:
    @:ASSERT(.not. this%tReadShifts)
    @:ASSERT(.not. this%tWriteShifts)

    this%tReadChrgAscii  = input%ctrl%tReadChrgAscii
    this%tWriteChrgAscii = input%ctrl%tWriteChrgAscii
    this%tSkipChrgChecksum = input%ctrl%tSkipChrgChecksum .or. this%tNegf

    call this%initializeCharges(input%ctrl%initialSpins, input%ctrl%initialCharges)

    ! Initialise images (translations)
    if (this%tPeriodic .or. this%tHelical) then
      call getCellTranslations(this%cellVec, this%rCellVec, this%latVec, this%invLatVec,&
          & this%cutOff%mCutOff)
    else
      allocate(this%cellVec(3, 1))
      allocate(this%rCellVec(3, 1))
      this%cellVec(:,1) = [0.0_dp, 0.0_dp, 0.0_dp]
      this%rCellVec(:,1) = [0.0_dp, 0.0_dp, 0.0_dp]
    end if

    ! Initialize neighbourlist.
    allocate(this%neighbourList)
    call TNeighbourlist_init(this%neighbourList, this%nAtom, nInitNeighbour)
    allocate(this%nNeighbourSK(this%nAtom))
    allocate(this%nNeighbourRep(this%nAtom))
    if (this%isRangeSep) then
      allocate(this%nNeighbourLC(this%nAtom))
    end if

    ! Set various options
    this%tWriteAutotest = env%tGlobalLead .and. input%ctrl%tWriteTagged
    this%tWriteDetailedXML = env%tGlobalLead .and. input%ctrl%tWriteDetailedXML
    this%tWriteResultsTag = env%tGlobalLead .and. input%ctrl%tWriteResultsTag
    this%tWriteDetailedOut = env%tGlobalLead .and. input%ctrl%tWriteDetailedOut .and.&
        & .not. this%tRestartNoSC
    this%tWriteBandDat = input%ctrl%tWriteBandDat .and. env%tGlobalLead&
        & .and. this%electronicSolver%providesEigenvals

    ! Check if stopfiles already exist and quit if yes
    inquire(file=fStopSCC, exist=tExist)
    if (tExist) then
      call error("Stop file '" // fStopSCC // "' already present at startup")
    end if
    inquire(file=fStopDriver, exist=tExist)
    if (tExist) then
      call error("Stop file '" // fStopDriver // "' already present at startup")
    end if

    this%restartFreq = input%ctrl%restartFreq

  #:if WITH_TRANSPORT
    if (this%tLatOpt .and. this%tNegf) then
      call error("Lattice optimisation currently incompatible with transport calculations")
    end if

    this%tUpload = input%transpar%taskUpload
    ! NOTE: originally EITHER 'contact calculations' OR 'upload' was possible
    !       introducing 'TransportOnly' option the logic is bit more
    !       involved: Contacts are not uploded in case of non-scc calculations
    if (this%electronicSolver%iSolver == electronicSolverTypes%OnlyTransport .and.&
        & .not.this%tSccCalc) then
      this%tUpload = .false.
    end if

    call this%initTransportArrays(input%transpar)
    call this%initTransport(env, input)

  #:else
    this%tNegf = .false.
    this%tUpload = .false.
  #:endif

    if (this%tPoisson) then
      this%poissStr%nAtom = this%nAtom
      this%poissStr%nSpecies = this%nType
      this%poissStr%specie0 => this%species0
      this%poissStr%x0 => this%coord0
      this%poissStr%nel = this%nEl0
      this%poissStr%isPeriodic = this%tPeriodic
      if (this%tPeriodic) then
        this%poissStr%latVecs(:,:) = this%latVec(:,:)
      else
        this%poissStr%latVecs(:,:) = 0.0_dp
      end if
      this%poissStr%tempElec = this%tempElec

      call poiss_init(env, this%poissStr, this%orb, this%hubbU, input%poisson,&
        #:if WITH_TRANSPORT
          & input%transpar,&
        #:endif
          & tInitialized)

      if (.not. tInitialized) then
        call error("Poisson solver not initialized")
      end if
    end if

    this%tPoissonTwice = input%poisson%solveTwice

    if (this%tNegf) then
      if (allocated(this%dispersion)) then
        call error("Dispersion not currently avalable with transport calculations")
      end if
      if (this%isLinResp) then
        call error("Linear response is not compatible with transport calculations")
      end if
      if (this%nSpin > 2) then
        call error("Non-colinear spin not currently compatible with transport calculations")
      end if
      if (allocated(this%solvation)) then
        call error("Solvation is currently not available with transport calculations")
      end if
    end if

    if (env%tGlobalLead) then
      call this%initOutputFiles(env)
    end if

    if (this%tPoisson) then
      this%electrostatics = elstatTypes%poisson
    else
      this%electrostatics = elstatTypes%gammaFunc
    end if

  #:if WITH_SCALAPACK
    associate (blacsOpts => input%ctrl%parallelOpts%blacsOpts)
      call getDenseDescBlacs(env, blacsOpts%blockSize, blacsOpts%blockSize, this%denseDesc)
    end associate
  #:endif

    if (allocated(this%reks)) then
      call checkReksConsistency(input%ctrl%reksInp, this%solvation, this%onSiteElements,&
          & this%kPoint, this%nEl, this%nKPoint, this%tSccCalc, this%tSpin, this%tSpinOrbit,&
          & allocated(this%dftbU), this%tEField, this%isLinResp, this%tPeriodic, this%tLatOpt,&
          & this%tReadChrg, this%tPoisson, input%ctrl%tShellResolved)
      ! here, this%nSpin changes to 2 for REKS
      call TReksCalc_init(this%reks, input%ctrl%reksInp, this%electronicSolver, this%orb,&
          & this%spinW, this%nEl, input%ctrl%extChrg, input%ctrl%extChrgBlurWidth,&
          & this%hamiltonianType, this%nSpin, this%nExtChrg, this%t3rd.or.this%t3rdFull,&
          & this%isRangeSep, this%tForces, this%tPeriodic, this%tStress, this%tDipole)
    end if

    call this%initDetArrays()
    call this%initArrays(env)

  #:if WITH_TRANSPORT
    ! note, this has the side effect of setting up module variable transpar as copy of
    ! input%transpar

    if (this%tUpload) then
      ! check geometry details are consistent with transport with contacts
      call checkTransportRanges(this%nAtom, input%transpar)
    end if

    if (this%tContCalc) then
      ! geometry is reduced to contacts only
      allocate(this%iAtInCentralRegion(this%nAtom))
    else
      allocate(this%iAtInCentralRegion(this%transpar%idxdevice(2)))
    end if

    if (this%transpar%tPeriodic1D) then
      if ( any(abs(this%kPoint(2:3, :)) > 0.0_dp) ) then
        call error("For transport in wire-like cases, only k-points in the first index should be&
            & non-zero")
      end if
    end if

    if (this%transpar%taskUpload .and. this%transpar%ncont > 0) then
      if (this%tPeriodic .and. .not. this%transpar%tPeriodic1D) then
        do ii = 1, this%transpar%ncont
          do jj = 1, 3
            if (abs(dot_product(this%transpar%contacts(ii)%lattice, this%latVec(:,jj)))&
                & > epsilon(0.0) .and. any(abs(this%kPoint(jj,:)) > 0.0_dp)) then
              call error("The k-points along transport direction(s) should zero in that direction")
            end if
          end do
        end do
      end if
    end if

  #:else
    allocate(this%iAtInCentralRegion(this%nAtom))
  #:endif
    ! atoms in central cell/device region/all atoms depending on boundary conditions
    do iAt = 1, size(this%iAtInCentralRegion)
      this%iAtInCentralRegion(iAt) = iAt
    end do

    if (this%tShowFoldedCoord) then
      this%pCoord0Out => this%coord0Fold
    else
      this%pCoord0Out => this%coord0
    end if


    ! Projection of eigenstates onto specific regions of the system
    this%tProjEigenvecs = input%ctrl%tProjEigenvecs
    if (this%tProjEigenvecs) then
      call init(this%iOrbRegion)
      call init(this%regionLabels)
      do iReg = 1, size(input%ctrl%tShellResInRegion)
        call elemShape(input%ctrl%iAtInRegion, valshape, iReg)
        nAtomRegion = valshape(1)
        allocate(iAtomRegion(nAtomRegion))
        call intoArray(input%ctrl%iAtInRegion, iAtomRegion, iTmp, iReg)
        if (input%ctrl%tOrbResInRegion(iReg) .or. input%ctrl%tShellResInRegion(iReg)) then

          if (input%ctrl%tOrbResInRegion(iReg)) then
            iSp = this%species0(iAtomRegion(1)) ! all atoms the same in the region
            @:ASSERT(all(this%species0(iAtomRegion) == iSp))
            nOrbRegion = nAtomRegion
            ! Create orbital index.
            allocate(tmpir1(nOrbRegion))
            do iOrb = 1, this%orb%nOrbSpecies(iSp)
              tmpir1 = 0
              ind = 1
              do iAt = 1, nAtomRegion
                tmpir1(ind) = this%denseDesc%iAtomStart(iAtomRegion(iAt)) + iOrb - 1
                ind = ind + 1
              end do
              call append(this%iOrbRegion, tmpir1)
              write(tmpStr, "(A,'.',I0,'.',I0,'.out')")trim(input%ctrl%RegionLabel(iReg)),&
                  & this%orb%iShellOrb(iOrb,iSp), iOrb&
                  & - this%orb%posShell(this%orb%iShellOrb(iOrb,iSp),iSp)&
                  & - this%orb%angShell(this%orb%iShellOrb(iOrb,iSp),iSp)
              call append(this%regionLabels, tmpStr)
            end do
            deallocate(tmpir1)
          end if

          if (input%ctrl%tShellResInRegion(iReg)) then
            iSp = this%species0(iAtomRegion(1)) ! all atoms the same in the region
            @:ASSERT(all(this%species0(iAtomRegion) == iSp))
            ! Create a separate region for each shell. It will contain
            ! the orbitals of that given shell for each atom in the region.
            do iSh = 1, this%orb%nShell(iSp)
              nOrbRegion = nAtomRegion * (this%orb%posShell(iSh + 1, iSp) -&
                  & this%orb%posShell(iSh, iSp))
              ind = 1
              ! Create orbital index.
              allocate(tmpir1(nOrbRegion))
              do ii = 1, nAtomRegion
                iAt = iAtomRegion(ii)
                do jj = this%orb%posShell(iSh, iSp), this%orb%posShell(iSh + 1, iSp) - 1
                  tmpir1(ind) = this%denseDesc%iAtomStart(iAt) + jj - 1
                  ind = ind + 1
                end do
              end do
              call append(this%iOrbRegion, tmpir1)
              deallocate(tmpir1)
              write(tmpStr, "(A,'.',I0,'.out')")trim(input%ctrl%RegionLabel(iReg)), iSh
              call append(this%regionLabels, tmpStr)
            end do
          end if

        else
          ! We take all orbitals from all atoms.
          nOrbRegion = 0
          do ii = 1, nAtomRegion
            nOrbRegion = nOrbRegion + this%orb%nOrbAtom(iAtomRegion(ii))
          end do
          ind = 1
          allocate(tmpir1(nOrbRegion))
          ! Create an index of the orbitals
          do ii = 1, nAtomRegion
            iAt = iAtomRegion(ii)
            do jj = 1, this%orb%nOrbAtom(iAt)
              tmpir1(ind) = this%denseDesc%iAtomStart(iAt) + jj - 1
              ind = ind + 1
            end do
          end do
          call append(this%iOrbRegion, tmpir1)
          deallocate(tmpir1)
          write(tmpStr, "(A,'.out')") trim(input%ctrl%RegionLabel(iReg))
          call append(this%regionLabels, tmpStr)
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

  #:if WITH_OMP
    write(stdOut, "('OpenMP threads: ', T30, I0)") omp_get_max_threads()
  #:endif

  #:if WITH_MPI and WITH_OMP
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
      write(stdOut, "(A,':',T30,I0)") "Chosen random seed", iSeed
    else
      write(stdOut, "(A,':',T30,I0)") "Specified random seed", iSeed
    end if

    if (input%ctrl%tMD) then
      select case(input%ctrl%iThermostat)
      case (0)
        if (this%tBarostat) then
          write(stdOut, "('Mode:',T30,A,/,T30,A)") 'MD without scaling of velocities',&
              & '(a.k.a. "NPE" ensemble)'
        else
          write(stdOut, "('Mode:',T30,A,/,T30,A)") 'MD without scaling of velocities',&
              & '(a.k.a. NVE ensemble)'
        end if
      case (1)
        if (this%tBarostat) then
          write(stdOut, "('Mode:',T30,A,/,T30,A)")&
              & "MD with re-selection of velocities according to temperature",&
              & "(a.k.a. NPT ensemble using Andersen thermostating + Berensen barostat)"
        else
          write(stdOut, "('Mode:',T30,A,/,T30,A)")&
              & "MD with re-selection of velocities according to temperature",&
              & "(a.k.a. NVT ensemble using Andersen thermostating)"
        end if
      case(2)
        if (this%tBarostat) then
          write(stdOut, "('Mode:',T30,A,/,T30,A)")&
              & "MD with scaling of velocities according to temperature",&
              & "(a.k.a. 'not' NVP ensemble using Berendsen thermostating and barostat)"
        else
          write(stdOut, "('Mode:',T30,A,/,T30,A)")&
              & "MD with scaling of velocities according to temperature",&
              & "(a.k.a. 'not' NVT ensemble using Berendsen thermostating)"
        end if
      case(3)
        if (this%tBarostat) then
          write(stdOut, "('Mode:',T30,A,/,T30,A)")"MD with scaling of velocities according to",&
              & "Nose-Hoover-Chain thermostat + Berensen barostat"
        else
          write(stdOut, "('Mode:',T30,A,/,T30,A)")"MD with scaling of velocities according to",&
              & "Nose-Hoover-Chain thermostat"
        end if

      case default
        call error("Unknown thermostat mode")
      end select
    elseif (this%isGeoOpt) then
      if (allocated(this%conAtom)) then
        strTmp = "with constraints"
      else
        strTmp = ""
      end if
      select case (input%ctrl%iGeoOpt)
      case (geoOptTypes%steepestDesc)
        write(stdOut, "('Mode:',T30,A)")'Steepest descent' // trim(strTmp)
      case (geoOptTypes%conjugateGrad)
        write(stdOut, "('Mode:',T30,A)") 'Conjugate gradient relaxation' // trim(strTmp)
      case (geoOptTypes%diis)
        write(stdOut, "('Mode:',T30,A)") 'Modified gDIIS relaxation' // trim(strTmp)
      case (geoOptTypes%lbfgs)
        write(stdout, "('Mode:',T30,A)") 'LBFGS relaxation' // trim(strTmp)
      case (geoOptTypes%fire)
        write(stdout, "('Mode:',T30,A)") 'FIRE relaxation' // trim(strTmp)
      case default
        call error("Unknown optimisation mode")
      end select
    elseif (this%tDerivs) then
      write(stdOut, "('Mode:',T30,A)") "2nd derivatives calculation"
      write(stdOut, "('Mode:',T30,A)") "Calculated for atoms:"
      write(stdOut, *) this%indMovedAtom
    elseif (this%tSocket) then
      write(stdOut, "('Mode:',T30,A)") "Socket controlled calculation"
    else
      write(stdOut, "('Mode:',T30,A)") "Static calculation"
    end if

    if (this%tSccCalc) then
      if (.not.this%tRestartNoSC) then
        write(stdOut, "(A,':',T30,A)") "Self consistent charges", "Yes"
        write(stdOut, "(A,':',T30,E14.6)") "SCC-tolerance", this%sccTol
        write(stdOut, "(A,':',T30,I14)") "Max. scc iterations", this%maxSccIter
      end if
      if (this%tPeriodic) then
        write(stdout, "(A,':',T30,E14.6)") "Ewald alpha parameter", this%sccCalc%getEwaldPar()
      end if
      if (input%ctrl%tShellResolved) then
         write(stdOut, "(A,':',T30,A)") "Shell resolved Hubbard", "Yes"
      else
         write(stdOut, "(A,':',T30,A)") "Shell resolved Hubbard", "No"
      end if
      if (allocated(this%dftbU)) then
        write(stdOut, "(A,':',T35,A)")"Orbitally dependant functional", "Yes"
        write(stdOut, "(A,':',T30,A)")"Orbital functional", this%dftbU%funcName()
      end if
      if (allocated(this%onSiteElements)) then
        write(stdOut, "(A,':',T35,A)")"On-site corrections", "Yes"
      end if
    else
      write(stdOut, "(A,':',T30,A)") "Self consistent charges", "No"
    end if

    if (allocated(this%reks)) then
      write(stdOut, "(A,':',T30,A)") "Spin polarisation", "No"
      write(stdOut, "(A,':',T30,F12.6,/,A,':',T30,F12.6)") "Nr. of up electrons",&
          & 0.5_dp*this%nEl(1), "Nr. of down electrons", 0.5_dp*this%nEl(1)
    else
      select case (this%nSpin)
      case(1)
        write(stdOut, "(A,':',T30,A)") "Spin polarisation", "No"
        write(stdOut, "(A,':',T30,F12.6,/,A,':',T30,F12.6)") "Nr. of up electrons",&
            & 0.5_dp*this%nEl(1), "Nr. of down electrons", 0.5_dp*this%nEl(1)
      case(2)
        write(stdOut, "(A,':',T30,A)") "Spin polarisation", "Yes"
        write(stdOut, "(A,':',T30,F12.6,/,A,':',T30,F12.6)") "Nr. of up electrons", this%nEl(1),&
            & "Nr. of down electrons", this%nEl(2)
      case(4)
        write(stdOut, "(A,':',T30,A)") "Non-collinear calculation", "Yes"
        write(stdOut, "(A,':',T30,F12.6)") "Nr. of electrons", this%nEl(1)
      end select
    end if

    if (this%tPeriodic) then
      write(stdOut, "(A,':',T30,A)") "Periodic boundaries", "Yes"
      if (this%tLatOpt) then
        write(stdOut, "(A,':',T30,A)") "Lattice optimisation", "Yes"
        write(stdOut, "(A,':',T30,f12.6)") "Pressure", this%extPressure
      end if
    else if (this%tHelical) then
      write (stdOut, "(A,':',T30,A)") "Helical boundaries", "Yes"
      if (this%tLatOpt) then
        write (stdOut, "(A,':',T30,A)") "Lattice optimisation", "Yes"
      end if
    else
      write(stdOut, "(A,':',T30,A)") "Periodic boundaries", "No"
    end if

    if (.not.this%tRestartNoSC) then
      write(stdOut, "(A,':',T30,A)") "Electronic solver", this%electronicSolver%getSolverName()
    end if

    if (this%electronicSolver%iSolver == electronicSolverTypes%magma_gvd) then
      #:if WITH_GPU
        call env%initGpu()
      #:else
        call error("Magma-solver selected, but program was compiled without MAGMA")
      #:endif
    endif

    if (this%tSccCalc .and. .not.this%tRestartNoSC) then
      if (.not. allocated(this%reks)) then
        select case (iMixer)
        case(mixerTypes%simple)
          write (strTmp, "(A)") "Simple"
        case(mixerTypes%anderson)
          write (strTmp, "(A)") "Anderson"
        case(mixerTypes%broyden)
          write (strTmp, "(A)") "Broyden"
        case(mixerTypes%diis)
          write (strTmp, "(A)") "DIIS"
        end select
        write(stdOut, "(A,':',T30,A,' ',A)") "Mixer", trim(strTmp), "mixer"
        write(stdOut, "(A,':',T30,F14.6)") "Mixing parameter", mixParam
        write(stdOut, "(A,':',T30,I14)") "Maximal SCC-cycles", this%maxSccIter
        select case (iMixer)
        case(mixerTypes%anderson)
          write(stdOut, "(A,':',T30,I14)") "Nr. of chrg. vectors to mix", nGeneration
        case(mixerTypes%broyden)
          write(stdOut, "(A,':',T30,I14)") "Nr. of chrg. vec. in memory", nGeneration
        case(mixerTypes%diis)
          write(stdOut, "(A,':',T30,I14)") "Nr. of chrg. vectors to mix", nGeneration
        end select
      else
        write(stdOut, "(A,':',T30,I14)") "Maximal SCC-cycles", this%maxSccIter
      end if
    end if

    if (this%tCoordOpt) then
      write(stdOut, "(A,':',T30,I14)") "Nr. of moved atoms", this%nMovedAtom
    end if
    if (this%isGeoOpt) then
      write(stdOut, "(A,':',T30,I14)") "Max. nr. of geometry steps", this%nGeoSteps
      write(stdOut, "(A,':',T30,E14.6)") "Force tolerance", input%ctrl%maxForce
      if (input%ctrl%iGeoOpt == geoOptTypes%steepestDesc) then
        write(stdOut, "(A,':',T30,E14.6)") "Step size", this%deltaT
      end if
    end if

    tFirst = .true.
    if (allocated(this%conAtom)) then
      do ii = 1, this%nAtom
        do jj = 1, size(this%conAtom)
          if (this%conAtom(jj) == ii) then
            if (tFirst) then
              write(strTmp, "(A,':')") "Geometry constraints"
              tFirst = .false.
            else
              write(strTmp, "(A)") ""
            end if
            write(stdOut, "(A,T30,'At',I4,': ',3F10.6)") trim(strTmp), ii, (this%conVec(kk,jj),&
                & kk=1,3)
          end if
        end do
      end do
    end if

    if (.not. allocated(this%reks) .and. .not.this%tRestartNoSC) then
      if (.not.input%ctrl%tSetFillingTemp) then
        write(stdOut, format2Ue) "Electronic temperature", this%tempElec, 'H',&
            & Hartree__eV * this%tempElec, 'eV'
      end if
    end if
    if (this%tMD) then
      write(stdOut, "(A,':',T30,E14.6)") "Time step", this%deltaT
      if (input%ctrl%iThermostat == 0 .and. .not.input%ctrl%tReadMDVelocities) then
        write(stdOut, "(A,':',T30,E14.6)") "Temperature", this%tempAtom
      end if
      if (input%ctrl%tRescale) then
        write(stdOut, "(A,':',T30,E14.6)") "Rescaling probability", input%ctrl%wvScale
      end if
    end if

    if (this%tSccCalc .and. .not.this%tRestartNoSC) then
      if (this%tReadChrg) then
        write (strTmp, "(A,A,A)") "Read in from '", trim(fCharges), "'"
      else
        write (strTmp, "(A,E11.3,A)") "Set automatically (system chrg: ", input%ctrl%nrChrg, ")"
      end if
      write(stdOut, "(A,':',T30,A)") "Initial charges", trim(strTmp)
    end if

    do iSp = 1, this%nType
      call getShellNames(iSp, this%orb, shellNamesTmp)
      if (iSp == 1) then
        write (strTmp, "(A,':')") "Included shells"
      else
        write (strTmp, "(A)") ""
      end if
      do jj = 1, this%orb%nShell(iSp)
        if (jj == 1) then
          strTmp2 = trim(shellNamesTmp(jj))
        else
          strTmp2 = trim(strTmp2) // ", " // trim(shellNamesTmp(jj))
        end if
      end do
      write(stdOut, "(A,T29,A2,':  ',A)") trim(strTmp), trim(this%speciesName(iSp)), trim(strTmp2)
      deallocate(shellNamesTmp)
    end do

    if (this%tMulliken) then
      if (allocated(input%ctrl%customOccAtoms)) then
        call printCustomReferenceOccupations(this%orb, input%geom%species, &
            & input%ctrl%customOccAtoms, input%ctrl%customOccFillings)
      end if
    end if

    if (this%tPeriodic) then
      do ii = 1, this%nKPoint
        if (ii == 1) then
          write(strTmp, "(A,':')") "K-points and weights"
        else
          write(strTmp, "(A)") ""
        end if
        write(stdOut, "(A,T28,I6,':',3F10.6,3X,F10.6)") trim(strTmp), ii,&
            & (this%kPoint(jj, ii), jj=1, 3), this%kWeight(ii)
      end do
      write(stdout,*)
      do ii = 1, this%nKPoint
        if (ii == 1) then
          write(strTmp, "(A,':')") "K-points in absolute space"
        else
          write(strTmp, "(A)") ""
        end if
        write(stdout, "(A,T28,I6,':',3F10.6)") trim(strTmp), ii,&
            & matmul(this%invLatVec,this%kPoint(:,ii))
      end do
      write(stdout, *)
    end if

    if (this%tHelical) then
      do ii = 1, this%nKPoint
        if (ii == 1) then
          write(strTmp, "(A,':')") "K-points and weights"
        else
          write(strTmp, "(A)") ""
        end if
        write(stdOut,"(A,T28,I6,':',2F10.6,3X,F10.6)") trim(strTmp), ii, this%kPoint(:, ii),&
            & this%kWeight(ii)
      end do
    end if

    if (allocated(this%dispersion)) then
      select type (o=>this%dispersion)
      type is (TDispSlaKirk)
        write(stdOut, "(A)") "Using Slater-Kirkwood dispersion corrections"
      type is (TDispUff)
        write(stdOut, "(A)") "Using Lennard-Jones dispersion corrections"
    #:if WITH_DFTD3
      type is (TDispDftD3)
        write(stdOut, "(A)") "Using DFT-D3 dispersion corrections"
    #:endif
      type is (TSimpleDftD3)
        write(stdOut, "(A)") "Using simple DFT-D3 dispersion corrections"
      type is (TDispDftD4)
        write(stdOut, "(A)") "Using DFT-D4 dispersion corrections"
    #:if WITH_MBD
      type is (TDispMbd)
        call writeMbdInfo(input%ctrl%dispInp%mbd)
    #:endif
      class default
        call error("Unknown dispersion model - this should not happen!")
      end select
    end if

    if (allocated(this%solvation)) then
      call writeSolvationInfo(stdOut, this%solvation)
    end if

    if (this%tSccCalc) then
      ! Have the SK values of U been replaced?
      if (allocated(input%ctrl%hubbU)) then
        strTmp = ""
        tFirst = .true.
        do iSp = 1, this%nType
          if (all(input%ctrl%hubbU(1:this%orb%nShell(iSp), iSp) == 0.0_dp)) then
            cycle
          end if
          do jj = 1, this%orb%nShell(iSp)
            if (tFirst) then
              write(strTmp, "(A,':')") "Non-default Hubbard U"
              tFirst = .false.
            else
              write(strTmp, "(A)") ""
            end if
            write(stdOut, "(A,T30,A2,2X,I1,'(',A1,'): ',E14.6)") trim(strTmp),&
                & this%speciesName(iSp), jj, shellNames(this%orb%angShell(jj, iSp)+1),&
                & this%hubbU(jj, iSp)
          end do
        end do
      end if
    end if

    tFirst = .true.
    if (this%tSpin .or. allocated(this%reks)) then
      do iSp = 1, this%nType
        do jj = 1, this%orb%nShell(iSp)
          do kk = 1, this%orb%nShell(iSp)
            if (tFirst) then
              write(strTmp, "(A,':')") "Spin coupling constants"
              tFirst = .false.
            else
              write(strTmp, "(A)") ""
            end if
            if (allocated(this%reks)) then
              write(stdOut, "(A,T30,A2,2X,I1,'(',A1,')-',I1,'(',A1,'): ',E14.6)")trim(strTmp),&
                  & this%speciesName(iSp), jj, shellNames(this%orb%angShell(jj, iSp)+1), kk,&
                  & shellNames(this%orb%angShell(kk, iSp)+1), this%spinW(kk, jj, iSp) /&
                  & this%reks%Tuning(iSp)
            else
              write(stdOut, "(A,T30,A2,2X,I1,'(',A1,')-',I1,'(',A1,'): ',E14.6)")trim(strTmp),&
                  & this%speciesName(iSp), jj, shellNames(this%orb%angShell(jj, iSp)+1), kk,&
                  & shellNames(this%orb%angShell(kk, iSp)+1), this%spinW(kk, jj, iSp)
            end if
          end do
        end do
      end do
    end if

    tFirst = .true.
    if (this%tSpinOrbit) then
      if (this%tDualSpinOrbit) then
        write(stdOut, "(A)")"Dual representation spin orbit"
      end if
      do iSp = 1, this%nType
        do jj = 1, this%orb%nShell(iSp)
          if (tFirst) then
            write(strTmp, "(A,':')") "Spin orbit constants"
            tFirst = .false.
          else
            write(strTmp, "(A)") ""
          end if
          write(stdOut, "(A,T30,A2,2X,I1,'(',A1,'): ',E14.6)")trim(strTmp), this%speciesName(iSp),&
                & jj, shellNames(this%orb%angShell(jj, iSp)+1), this%xi(jj, iSp)
          if (this%xi(jj, iSp) /= 0.0_dp .and. this%orb%angShell(jj, iSp) == 0) then
            call error("Program halt due to non-zero s-orbital spin-orbit coupling constant!")
          end if
        end do
      end do
    end if

    if (this%tSccCalc) then
      if (this%t3rdFull) then
        write(stdOut, "(A,T30,A)") "Full 3rd order correction", "Yes"
        if (input%ctrl%tShellResolved) then
          write(stdOut, "(A,T30,A)") "Shell-resolved 3rd order", "Yes"
          write(stdOut, "(A30)") "Shell-resolved Hubbard derivs:"
          write(stdOut, "(A)") "        s-shell   p-shell   d-shell   f-shell"
          do iSp = 1, this%nType
            write(stdOut, "(A3,A3,4F10.4)") "  ", trim(this%speciesName(iSp)),&
                & input%ctrl%hubDerivs(:this%orb%nShell(iSp),iSp)
          end do
        end if
      end if

      if (any(tDampedShort)) then
        write(stdOut, "(A,T30,A)") "Damped SCC", "Yes"
        ii = count(tDampedShort)
        write(strTmp, "(A,I0,A)") "(A,T30,", ii, "(A,1X))"
        write(stdOut, strTmp) "Damped species(s):", pack(this%speciesName, tDampedShort)
        deallocate(tDampedShort)
      end if

      if (input%ctrl%h5SwitchedOn) then
        write(stdOut, "(A,T30,A)") "H-bond correction:", "H5"
      end if
      if (tHHRepulsion) then
        write(stdOut, "(A,T30,A)") "H-H repulsion correction:", "H5"
      end if
    end if

    if (this%isRangeSep) then
      write(stdOut, "(A,':',T30,A)") "Range separated hybrid", "Yes"
      write(stdOut, "(2X,A,':',T30,E14.6)") "Screening parameter omega",&
          & input%ctrl%rangeSepInp%omega

      select case(input%ctrl%rangeSepInp%rangeSepAlg)
      case (rangeSepTypes%neighbour)
        write(stdOut, "(2X,A,':',T30,2X,A)") "Range separated algorithm", "NeighbourBased"
        write(stdOut, "(2X,A,':',T30,E14.6,A)") "Spatially cutoff at",&
            & input%ctrl%rangeSepInp%cutoffRed * Bohr__AA," A"
      case (rangeSepTypes%threshold)
        write(stdOut, "(2X,A,':',T30,2X,A)") "Range separated algorithm", "Thresholded"
        write(stdOut, "(2X,A,':',T30,E14.6)") "Thresholded to",&
            & input%ctrl%rangeSepInp%screeningThreshold
      case (rangeSepTypes%matrixBased)
        write(stdOut, "(2X,A,':',T30,2X,A)") "Range separated algorithm", "MatrixBased"
      case default
        call error("Unknown range separated hybrid method")
      end select
    end if


    write(stdOut, "(A,':')") "Extra options"
    if (this%tPrintMulliken) then
      write(stdOut, "(T30,A)") "Mulliken analysis"
    end if
    if (this%tPrintForces .and. .not. (this%tMD .or. this%isGeoOpt .or. this%tDerivs)) then
      write(stdOut, "(T30,A)") "Force calculation"
    end if
    if (this%tForces) then
      select case (this%forceType)
      case(forceTypes%orig)
        write(stdOut, "(A,T30,A)") "Force type", "original"
      case(forceTypes%dynamicT0)
        write(stdOut, "(A,T30,A)") "Force type", "erho with re-diagonalized eigenvalues"
        write(stdOut, "(A,T30,A)") "Force type", "erho with DHD-product (T_elec = 0K)"
      case(forceTypes%dynamicTFinite)
        write(stdOut, "(A,T30,A)") "Force type", "erho with S^-1 H D (Te <> 0K)"
      end select
    end if
    if (this%tPrintEigVecs) then
      write(stdOut, "(T30,A)") "Eigenvector printing"
    end if
    if (this%tExtChrg) then
      write(stdOut, "(T30,A)") "External charges specified"
    end if

    if (this%tEField) then
      if (this%tTDEfield) then
        write(stdOut, "(T30,A)") "External electric field specified"
        write(stdOut, "(A,':',T30,E14.6)") "Angular frequency", this%EfieldOmega
      else
        write(stdOut, "(T30,A)") "External static electric field specified"
      end if
      write(stdOut, "(A,':',T30,E14.6)") "Field strength", this%EFieldStrength
      write(stdOut, "(A,':',T30,3F9.6)") "Direction", this%EfieldVector
      if (this%tPeriodic) then
        call warning("Saw tooth potential used for periodic geometry - make sure there is a vacuum&
            & region!")
      end if
    end if

    if (allocated(this%dftbU)) then
      do iSp = 1, this%nType
        if (this%dftbU%nUJ(iSp)>0) then
          write(strTmp, "(A,':')") "U-J coupling constants"
          write(stdOut, "(A,T25,A2)")trim(strTmp), this%speciesName(iSp)
          do jj = 1, this%dftbU%nUJ(iSp)
            write(strTmp, "(A,I1,A)")'(A,',this%dftbU%niUJ(jj,iSp),'I2,T25,A,F6.4)'
            write(stdOut, trim(strTmp)) 'Shells:',&
                & this%dftbU%iUJ(1:this%dftbU%niUJ(jj,iSp),jj,iSp), 'UJ:', this%dftbU%UJ(jj,iSp)
          end do
        end if
      end do
    end if

    tFirst = .true.
    if (allocated(this%onSiteElements)) then
      do iSp = 1, this%nType
        do iSpin = 1, 2
          if (iSpin == 1) then
            write(strTmp2, "(A,':')") "uu"
          else
            write(strTmp2, "(A,':')") "ud"
          end if
          do jj = 1, this%orb%nShell(iSp)
            do kk = 1, this%orb%nShell(iSp)
              if (tFirst) then
                write(strTmp, "(A,':')") "On-site coupling constants"
                tFirst = .false.
              else
                write(strTmp, "(A)") ""
              end if
              write(stdOut, "(A,T30,A5,2X,I1,'(',A1,')-',I1,'(',A1,'): ',E14.6)")trim(strTmp),&
                  & trim(this%speciesName(iSp))//trim(strTmp2), jj,&
                  & shellNames(this%orb%angShell(jj, iSp)+1), kk,&
                  & shellNames(this%orb%angShell(kk, iSp)+1),&
                  & this%onSiteElements(kk, jj, iSpin, iSp)
            end do
          end do
        end do
      end do
    end if

    if (this%deltaDftb%isNonAufbau) then
      if (this%nSpin /= 2) then
        call error("Internal error, Delta DFTB requires two spin channels")
      end if
      if (this%nEl(1) /= this%nEl(2)) then
        call error("Internal error, Delta DFTB requires a spin free reference")
      end if
      if (abs(this%nEl(1) - nint(this%nEl(1))) > epsilon(0.0)) then
        call error("Delta DFTB requires an integer number of electrons in the reference state")
      end if
      if (mod(sum(nint(this%nEl)),2) /= 0) then
        call error("Delta DFTB requires an even number of electrons in reference state")
      end if
      if (sum(this%nEl) >= 2*this%nOrb) then
        call error("Delta DFTB requires at least one empty orbita in the system")
      end if
      if (sum(this%nEl) < 2) then
        call error("Delta DFTB requires at least one full orbital in the system")
      end if
    end if
    if (this%deltaDftb%isNonAufbau .and. .not.this%tSccCalc) then
      call error("Delta DFTB must use SCC = Yes")
    end if
    if (this%deltaDftb%isNonAufbau .and. this%isLinResp) then
      call error("Delta DFTB incompatible with linear response")
    end if
    if (this%deltaDftb%isNonAufbau .and. allocated(this%ppRPA)) then
      call error("Delta DFTB incompatible with ppRPA")
    end if
    if (this%deltaDftb%isNonAufbau .and. allocated(input%ctrl%elecDynInp)) then
      call error("Delta DFTB incompatible with electron dynamics")
    end if
    if (this%deltaDftb%isNonAufbau .and. this%tFixEf) then
      call error("Delta DFTB incompatible with fixed Fermi energy")
    end if
    if (this%deltaDftb%isNonAufbau .and. this%tSpinSharedEf) then
      call error("Delta DFTB incompatible with shared Fermi energy")
    end if
    if (this%deltaDftb%isNonAufbau .and. allocated(this%reks)) then
      call error("Delta DFTB incompatible with REKS")
    end if
    if (this%deltaDftb%isNonAufbau .and. this%tNegf) then
      call error("Delta DFTB incompatible with transport")
    end if
    if (this%deltaDftb%isNonAufbau .and. this%tLocalise) then
      call error("Delta DFTB incompatible with localisation")
    end if

    if (this%tSpinOrbit .and. allocated(this%dftbU) .and. .not. this%tDualSpinOrbit)  then
      call error("Only dual spin orbit currently supported for orbital potentials")
    end if

    if (this%tHelical .and. this%tSpinOrbit) then
      call error("L.S coupling not yet supported for helical boundary conditions.")
    end if

    if (this%tHelical .and. this%nSpin > 2) then
      call error("Non-collinear not yet supported for helical boundary conditions.")
    end if

    if (.not.this%tStress) then
      if (this%tBarostat) then
        call error("Sorry, MD with a barostat requires stress evaluation")
      end if
      if (this%tLatOpt) then
        call error("Sorry, lattice optimization requires stress tensor evaluation")
      end if
    end if

    if (this%tSpinOrbit .and. (this%tWriteHS .or. (this%tWriteRealHS .and. .not.&
        & this%tDualSpinOrbit))) then
      call error("Writing of Hamiltonian currently not possible with spin orbit coupling enabled.")
    end if

    if (this%isLinResp) then
      if (allocated(this%dftbU)) then
        call error("Linear response is not compatible with Orbitally dependant functionals yet")
      end if

      if (this%t2Component) then
        call error("Linear response is not compatibile with this spinor Hamiltonian")
      end if

      if (this%tStress) then
        call error("Excited state stresses not implemented")
      end if

      if (.not.this%tRealHS) then
        call error("Linear response does not support k-points")
      end if

      if (this%t3rd .or. this%t3rdFull) then
        call error ("Third order DFTB is not currently compatible with linear response excitations")
      end if

    end if

    ! Electron dynamics stuff
    if (allocated(input%ctrl%elecDynInp)) then

      if (this%t2Component) then
        call error("Electron dynamics is not compatibile with this spinor Hamiltonian")
      end if

      if (withMpi) then
        call error("Electron dynamics does not work with MPI yet")
      end if

      if (this%tFixEf) then
        call error("Electron dynamics does not work with fixed Fermi levels yet")
      end if

      if (this%tSpinSharedEf) then
        call error("Electron dynamics does not work with spin shared Fermi levels yet")
      end if

      if (this%tMD) then
        call error("Electron dynamics does not work with MD")
      end if

      if (this%isRangeSep) then
        call error("Electron dynamics does not work with range separated calculations yet.")
      end if

      if (.not. this%tRealHS .and. input%ctrl%elecDynInp%tBondE) then
        call error("Bond energies during electron dynamics currently requires a real hamiltonian.")
      end if

      allocate(this%electronDynamics)

      call TElecDynamics_init(this%electronDynamics, input%ctrl%elecDynInp, this%species0,&
          & this%speciesName, this%tWriteAutotest, autotestTag, randomThermostat, this%mass,&
          & this%nAtom, this%cutOff%skCutoff, this%cutOff%mCutoff, this%atomEigVal,&
          & this%dispersion, this%nonSccDeriv, this%tPeriodic, this%parallelKS, this%tRealHS,&
          & this%kPoint, this%kWeight, this%isRangeSep)

    end if

    if (allocated(this%reks)) then
      call printReksInitInfo(this%reks, this%orb, this%speciesName, this%nType)
    end if

    call env%globalTimer%stopTimer(globalTimers%globalInit)

  end subroutine initProgramVariables


  !> Check coherence across processes for various key variables (relevant if running in MPI,
  !> particularly for external driving via API)
  subroutine inputCoherenceCheck(env, hamiltonianType, nSpin, nAtom, coord0, species0, speciesName,&
       & tSccCalc, tPeriodic, tFracCoord, latVec, origin)

    !> Environment settings
    type(TEnvironment), intent(in) :: env

    !> Hamiltonian type
    integer, intent(in) :: hamiltonianType

    !> Number of spin components
    integer, intent(in) :: nSpin

    !> Atoms in the system
    integer, intent(in) :: nAtom

    ! Atom coordinates (in the central unit cell, if relevant).
    real(dp), intent(in) :: coord0(:,:)

    !> Species of atoms in the central cell
    integer, intent(in) :: species0(:)

    !> Names of chemical species
    character(*), intent(in) :: speciesName(:)

    !> Is the calculation SCC?
    logical, intent(in) :: tSccCalc

    !> Is the calculation periodic?
    logical, intent(in) :: tPeriodic

    !> If periodic, are the atomic positions in fractional coordinates?
    logical, intent(in) :: tFracCoord

    !> lattice vectors, stored columnwise
    real(dp), intent(in) :: latVec(:,:)

    !> Origin of coordinate system for periodic systems
    real(dp), intent(in) :: origin(:)

    integer :: iSp

    call checkExactCoherence(env, hamiltonianType, "hamiltonianType in initProgramVariables")
    call checkExactCoherence(env, nSpin, "spin integer in initProgramVariables")
    call checkExactCoherence(env, nAtom, "the number of atoms in initProgramVariables")
    call checkToleranceCoherence(env, coord0, "coord0 in initProgramVariables", tol=1.e-10_dp)
    call checkExactCoherence(env, species0, "atomic species in initProgramVariables")
    call checkExactCoherence(env, tSccCalc, &
         & "the type of calculation, SCC, in initProgramVariables")
    do iSp = 1, size(speciesName)
       call checkExactCoherence(env, speciesName(iSp), "species names in initProgramVariables")
    enddo

    if (tPeriodic) then
       call checkExactCoherence(env, tFracCoord, "tFracCoord in initProgramVariables")
       call checkToleranceCoherence(env, latVec, &
            & "lattice vectors in initProgramVariables", tol=1.e-10_dp)
       call checkToleranceCoherence(env, origin, &
            & "coordinate origin in initProgramVariables", tol=1.e-10_dp)
    endif

  end subroutine inputCoherenceCheck


  !> Create equivalency relations
  ! Note, this routine should not be called
  subroutine setEquivalencyRelations(this)

    !> Instance
    class(TDftbPlusMain), intent(inout) :: this

    !> Orbital equivalency for SCC and Spin
    integer, allocatable :: iEqOrbSCC(:,:,:), iEqOrbSpin(:,:,:)
    !> Orbital equivalency for orbital potentials
    integer, allocatable :: iEqOrbDFTBU(:,:,:)

    if (this%tSccCalc) then
      if(.not. allocated(this%iEqOrbitals)) then
        allocate(this%iEqOrbitals(this%orb%mOrb, this%nAtom, this%nSpin))
      endif
      allocate(iEqOrbSCC(this%orb%mOrb, this%nAtom, this%nSpin))
      @:ASSERT(allocated(this%sccCalc))
      call this%sccCalc%getOrbitalEquiv(this%orb, this%species0, iEqOrbSCC)
      if (this%nSpin == 1) then
        this%iEqOrbitals(:,:,:) = iEqOrbSCC(:,:,:)
      else
        allocate(iEqOrbSpin(this%orb%mOrb, this%nAtom, this%nSpin))
        call Spin_getOrbitalEquiv(this%orb, this%species0, iEqOrbSpin)
        call OrbitalEquiv_merge(iEqOrbSCC, iEqOrbSpin, this%orb, this%iEqOrbitals)
        deallocate(iEqOrbSpin)
      end if
      deallocate(iEqOrbSCC)
      this%nIneqOrb = maxval(this%iEqOrbitals)
      this%nMixElements = this%nIneqOrb

      if (allocated(this%dftbU)) then
        allocate(iEqOrbSpin(this%orb%mOrb, this%nAtom, this%nSpin))
        allocate(iEqOrbDFTBU(this%orb%mOrb, this%nAtom, this%nSpin))
        call this%dftbU%getOrbitalEquiv(iEqOrbDFTBU, this%orb, this%species0)
        call OrbitalEquiv_merge(this%iEqOrbitals, iEqOrbDFTBU, this%orb, iEqOrbSpin)
        this%iEqOrbitals(:,:,:) = iEqOrbSpin(:,:,:)
        this%nIneqOrb = maxval(this%iEqOrbitals)
        deallocate(iEqOrbSpin)
        deallocate(iEqOrbDFTBU)
      end if

      if (allocated(this%onSiteElements)) then
        allocate(iEqOrbSpin(this%orb%mOrb, this%nAtom, this%nSpin))
        iEqOrbSpin(:,:,:) = 0.0_dp
        allocate(iEqOrbDFTBU(this%orb%mOrb, this%nAtom, this%nSpin))
        iEqOrbDFTBU(:,:,:) = 0.0_dp
        call Ons_getOrbitalEquiv(iEqOrbDFTBU, this%orb, this%species0)
        call OrbitalEquiv_merge(this%iEqOrbitals, iEqOrbDFTBU, this%orb, iEqOrbSpin)
        this%iEqOrbitals(:,:,:) = iEqOrbSpin(:,:,:)
        this%nIneqOrb = maxval(this%iEqOrbitals)
        deallocate(iEqOrbSpin)
        deallocate(iEqOrbDFTBU)
      end if

      if (allocated(this%onSiteElements)) then
        ! all onsite blocks are full of unique elements
        if(.not. allocated(this%iEqBlockOnSite)) then
          allocate(this%iEqBlockOnSite(this%orb%mOrb, this%orb%mOrb, this%nAtom, this%nSpin))
        endif
        if (this%tImHam) then
          if(.not. allocated(this%iEqBlockOnSiteLS))then
            allocate(this%iEqBlockOnSiteLS(this%orb%mOrb, this%orb%mOrb, this%nAtom, this%nSpin))
          endif
        end if
        call Ons_blockIndx(this%iEqBlockOnSite, this%iEqBlockOnSiteLS, this%nIneqOrb, this%orb)
        this%nMixElements = max(this%nMixElements, maxval(this%iEqBlockOnSite))
        if (allocated(this%iEqBlockOnSiteLS)) then
          this%nMixElements = max(this%nMixElements, maxval(this%iEqBlockOnSiteLS))
        end if
      else if (allocated(this%dftbU)) then
        ! only a sub-set of onsite blocks are reduced/expanded
        if(.not. allocated(this%iEqBlockDFTBU))then
          allocate(this%iEqBlockDFTBU(this%orb%mOrb, this%orb%mOrb, this%nAtom, this%nSpin))
        endif
        call this%dftbU%blockIndx(this%iEqBlockDFTBU, this%nIneqOrb, this%orb, this%species0)
        this%nMixElements = max(this%nMixElements, maxval(this%iEqBlockDFTBU)) ! as
        !  iEqBlockDFTBU does not include diagonal elements, so in the case of
        !  a purely s-block DFTB+U calculation, maxval(iEqBlockDFTBU) would
        !  return 0
        if (this%tImHam) then
          if(.not. allocated(this%iEqBlockDFTBULS))then
            allocate(this%iEqBlockDFTBULS(this%orb%mOrb, this%orb%mOrb, this%nAtom, this%nSpin))
          endif
          call this%dftbU%blockIndx(this%iEqBlockDFTBULS, this%nMixElements, this%orb,&
              & this%species0)
          this%nMixElements = max(this%nMixElements, maxval(this%iEqBlockDFTBULS))
        end if
      end if

    !Non-SCC
    else
      this%nIneqOrb = this%nOrb
      this%nMixElements = 0
    end if

  end subroutine setEquivalencyRelations


  !> Initialise partial charges
  !>
  subroutine initializeCharges(this, initialSpins, initialCharges)

    !> Instance
    class(TDftbPlusMain), intent(inout) :: this

    !> Initial spins
    real(dp), allocatable, intent(in) :: initialSpins(:,:)

    !> Set of atom-resolved atomic charges
    real(dp), allocatable, intent(in) :: initialCharges(:)

    !> Tolerance in difference between total charge and sum of initial charges
    real(dp), parameter :: deltaChargeTol = 1.e-4_dp

    integer :: iAt, iSp, iSh, ii, jj, iStart, iEnd, iS
    real(dp) :: rTmp
    character(lc) :: message
    logical :: tAllocate

    ! Charge arrays may have already been initialised
    @:ASSERT(size(this%species0) == this%nAtom)

    if (.not. allocated(this%qInput)) then
      allocate(this%qInput(this%orb%mOrb, this%nAtom, this%nSpin))
    endif
    this%qInput(:,:,:) = 0.0_dp

    if (.not. allocated(this%qOutput)) then
      allocate(this%qOutput(this%orb%mOrb, this%nAtom, this%nSpin))
    endif
    this%qOutput(:,:,:) = 0.0_dp

    if (allocated(this%reks)) then
      if (.not. allocated(this%qDiff)) then
        allocate(this%qDiff(this%orb%mOrb, this%nAtom, this%nSpin))
      endif
      this%qDiff(:,:,:) = 0.0_dp
    endif

    tAllocate = .false.
  #:if WITH_MBD
    if (allocated(this%dispersion)) then
      select type (o=>this%dispersion)
      type is (TDispMbd)
        tAllocate = .true.
      end select
    end if
  #:endif
    tAllocate = tAllocate .or. this%tNetAtomCharges
    if (tAllocate) then
      if (.not. allocated(this%qNetAtom)) then
        allocate(this%qNetAtom(this%nAtom))
      endif
      this%qNetAtom(:) = 0.0_dp
    end if

    if (this%tMixBlockCharges) then
      if (.not. allocated(this%reks)) then
        if (.not. allocated(this%qBlockIn)) then
          allocate(this%qBlockIn(this%orb%mOrb, this%orb%mOrb, this%nAtom, this%nSpin))
        endif
        this%qBlockIn(:,:,:,:) = 0.0_dp
      endif
      if (.not. allocated(this%qBlockOut)) then
        allocate(this%qBlockOut(this%orb%mOrb, this%orb%mOrb, this%nAtom, this%nSpin))
      endif
      this%qBlockOut(:,:,:,:) = 0.0_dp
      if (this%tImHam) then
        if(.not. allocated(this%qiBlockIn)) then
          allocate(this%qiBlockIn(this%orb%mOrb, this%orb%mOrb, this%nAtom, this%nSpin))
        endif
        this%qiBlockIn(:,:,:,:) = 0.0_dp
      end if
    end if

    if (this%tImHam) then
      if(.not. allocated(this%qiBlockOut))then
        allocate(this%qiBlockOut(this%orb%mOrb, this%orb%mOrb, this%nAtom, this%nSpin))
      endif
      this%qiBlockOut(:,:,:,:) = 0.0_dp
    end if

    if (.not. this%tSccCalc) return

    ! Charges read from file
    if (this%tReadChrg) then
      if (this%tFixEf .or. this%tSkipChrgChecksum) then
        ! do not check charge or magnetisation from file
        call initQFromFile(this%qInput, fCharges, this%tReadChrgAscii, this%orb, this%qBlockIn,&
            & this%qiBlockIn, this%deltaRhoIn)
      else
        ! check number of electrons in file
        if (this%nSpin /= 2) then
          call initQFromFile(this%qInput, fCharges, this%tReadChrgAscii, this%orb, this%qBlockIn,&
              & this%qiBlockIn, this%deltaRhoIn, nEl = sum(this%nEl))
        else
          ! check magnetisation in addition
          call initQFromFile(this%qInput, fCharges, this%tReadChrgAscii, this%orb, this%qBlockIn,&
              & this%qiBlockIn, this%deltaRhoIn, nEl = sum(this%nEl),&
              & magnetisation=this%nEl(1)-this%nEl(2))
        end if
      end if
    endif

    if (.not. allocated(this%reks)) then
      !Input charges packed into unique equivalence elements
      #:for NAME in [('this%qDiffRed'),('this%qInpRed'),('this%qOutRed')]
        if (.not. allocated(${NAME}$)) then
          allocate(${NAME}$(this%nMixElements))
        end if
        ${NAME}$(:) = 0.0_dp
      #:endfor
    end if


    !TODO(Alex) Could definitely split the code here
    if (allocated(this%reks)) return

    ! Charges not read from file
    notChrgRead: if (.not. this%tReadChrg) then

      if (allocated(initialCharges)) then
        if (abs(sum(initialCharges) - this%nrChrg) > deltaChargeTol) then
          write(message, "(A,G13.6,A,G13.6,A,A)") "Sum of initial charges does not match&
              & specified total charge. (", sum(initialCharges), " vs. ", this%nrChrg, ") ",&
              & "Your initial charge distribution will be rescaled."
          call warning(message)
        end if
        call initQFromAtomChrg(this%qInput, initialCharges, this%referenceN0, this%species0,&
            & this%speciesName, this%orb)
      else
        this%qInput(:,:,:) = this%q0
      end if

      if (.not. this%tSkipChrgChecksum) then
        ! Rescaling to ensure correct number of electrons in the system
        this%qInput(:,:,1) = this%qInput(:,:,1) *  sum(this%nEl) / sum(this%qInput(:,:,1))
      end if


      select case (this%nSpin)
      case (1)
        continue
      case (2)
        if (allocated(initialSpins)) then
          do ii = 1, this%nAtom
            ! does not actually matter if additional spin polarization pushes
            ! charges to <0 as the initial charges are not mixed in to later
            ! iterations
            this%qInput(1:this%orb%nOrbAtom(ii),ii,2) = this%qInput(1:this%orb%nOrbAtom(ii),ii,1)&
                & * initialSpins(1,ii) / sum(this%qInput(1:this%orb%nOrbAtom(ii),ii,1))
          end do
        else
          if (.not. this%tSkipChrgChecksum) then
            do ii = 1, this%nAtom
              this%qInput(1:this%orb%nOrbAtom(ii),ii,2) = this%qInput(1:this%orb%nOrbAtom(ii),ii,1)&
                  & * (this%nEl(1)-this%nEl(2)) / sum(this%qInput(:,:,1))
            end do
          end if
        end if
      case (4)
        if (this%tSpin) then
          if (.not. allocated(initialSpins)) then
            call error("Missing initial spins!")
          end if
          if (any(shape(initialSpins)/=(/3,this%nAtom/))) then
            call error("Incorrect shape initialSpins array!")
          end if
          ! Rescaling to ensure correct number of electrons in the system
          if (.not. this%tSkipChrgChecksum) then
            do ii = 1, this%nAtom
              do jj = 1, 3
                this%qInput(1:this%orb%nOrbAtom(ii),ii,jj+1) =&
                    & this%qInput(1:this%orb%nOrbAtom(ii),ii,1) *&
                    & initialSpins(jj,ii) / sum(this%qInput(1:this%orb%nOrbAtom(ii),ii,1))
              end do
            end do
          end if
        end if
      end select

      if (this%tMixBlockCharges) then
        this%qBlockIn = 0.0_dp
        do iS = 1, this%nSpin
          do iAt = 1, this%nAtom
            iSp = this%species0(iAt)
            do iSh = 1, this%orb%nShell(iSp)
              iStart = this%orb%posShell(iSh,iSp)
              iEnd = this%orb%posShell(iSh+1,iSp)-1
              rTmp = sum(this%qInput(iStart:iEnd,iAt,iS))
              rTmp = rTmp / real(iEnd+1-iStart,dp)
              do ii = iStart, iEnd
                this%qBlockIn(ii,ii,iAt,iS) = rTmp
              end do
            end do
          end do
        end do
        if (this%tImHam) then
          this%qiBlockIn = 0.0_dp
        end if
      end if

    endif notChrgRead

    !Swap from charge/magnetisation to up/down
    if (this%nSpin == 2) then
      call qm2ud(this%qInput)
      if (this%tMixBlockCharges) then
        call qm2ud(this%qBlockIn)
      end if
    end if

    call OrbitalEquiv_reduce(this%qInput, this%iEqOrbitals, this%orb, this%qInpRed(1:this%nIneqOrb))

    if (allocated(this%onSiteElements)) then
      call appendBlockReduced(this%qBlockIn, this%iEqBlockOnSite, this%orb, this%qInpRed)
      if (this%tImHam) then
        call appendBlockReduced(this%qiBlockIn, this%iEqBlockOnSiteLS, this%orb, this%qInpRed,&
            & isSkew=.true.)
      end if
    else if (allocated(this%dftbU)) then
      call appendBlockReduced(this%qBlockIn, this%iEqBlockDFTBU, this%orb, this%qInpRed)
      if (this%tImHam) then
        call appendBlockReduced(this%qiBlockIn, this%iEqBlockDFTBULS, this%orb, this%qInpRed,&
            & isSkew=.true.)
      end if
    end if

    !Convert up/down set back to charge/magnetization
    if (this%nSpin == 2) then
      call ud2qm(this%qInput)
      if (this%tMixBlockCharges) call ud2qm(this%qBlockIn)
    end if

  end subroutine initializeCharges


  !> Assign reference charge arrays, q0 and qShell0
  !
  subroutine initializeReferenceCharges(this, customOccAtoms, customOccFillings)

    !> Instance
    class(TDftbPlusMain), intent(inout) :: this

    !> Atom indices corresponding to user defined reference atomic charges
    !  Array of occupation arrays, one for each atom
    type(TWrappedInt1), allocatable, intent(in) :: customOccAtoms(:)

    !> User-defined reference atomic shell charges
    real(dp), allocatable, intent(in) :: customOccFillings(:,:)

    integer :: iAt, iSp, iSh

    if (.not. allocated(this%q0)) then
      allocate(this%q0(this%orb%mOrb, this%nAtom, this%nSpin))
    endif
    this%q0(:,:,:) = 0.0_dp

    if (allocated(customOccAtoms)) then
      if (this%isLinResp) then
        call error("Custom occupation not compatible with linear response")
      end if
      call applyCustomReferenceOccupations(customOccAtoms, customOccFillings, this%species0,&
          & this%orb, this%referenceN0, this%q0)
    else
      call initQFromShellChrg(this%q0, this%referenceN0, this%species0, this%orb)
    end if

    if (.not. allocated(this%qShell0)) then
      allocate(this%qShell0(this%orb%mShell, this%nAtom))
    endif

    do iAt = 1, this%nAtom
      iSp = this%species0(iAt)
      do iSh = 1, this%orb%nShell(iSp)
        this%qShell0(iSh,iAt) = sum(this%q0(this%orb%posShell(iSh,iSp):&
            & this%orb%posShell(iSh+1,iSp)-1, iAt, 1))
      end do
    end do

  end subroutine initializeReferenceCharges


  !> Set number of electrons
  !
  subroutine setNElectrons(this)

    !> Instance
    class(TDftbPlusMain), intent(inout) :: this

    @:ASSERT(allocated(this%q0))
    @:ASSERT(allocated(this%nEl))

    this%nEl0 = sum(this%q0(:,:,1))
    if (abs(this%nEl0 - nint(this%nEl0)) < elecTolMax) then
      this%nEl0 = nint(this%nEl0)
    end if
    this%nEl(:) = 0.0_dp
    if (this%nSpin == 1 .or. this%nSpin == 4) then
      this%nEl(1) = this%nEl0 - this%nrChrg
      if(ceiling(this%nEl(1)) > 2.0_dp * this%nOrb) then
        call error("More electrons than basis functions!")
      end if
    else
      this%nEl(1) = 0.5_dp * (this%nEl0 - this%nrChrg + this%nrSpinPol)
      this%nEl(2) = 0.5_dp * (this%nEl0 - this%nrChrg - this%nrSpinPol)
      if (any(ceiling(this%nEl(:)) > this%nOrb)) then
        call error("More electrons than basis functions!")
      end if
    end if

    if (.not.all(this%nEl >= 0.0_dp)) then
      call error("Less than 0 electrons!")
    end if

  end subroutine setNElectrons


#:if WITH_TRANSPORT
  !> Check for inconsistencies in transport atom region definitions
  subroutine checkTransportRanges(nAtom, transpar)

    !> Count of all atoms in the system
    integer :: nAtom

    !> Transport parameters
    type(TTransPar), intent(in) :: transpar

    character(lc) :: strTmp
    integer :: ii, jj
    logical :: tFailCheck
    logical, allocatable :: notInRegion(:)

    ! check for atoms occurring inside both the device and a contact
    do ii = 1, transpar%ncont
      if (transpar%contacts(ii)%idxrange(1)<=transpar%idxdevice(2)) then
        write(strTmp,"(A,I0,A,A,A,I0,A,I0)") "The device and contact overlap in their atom index&
            & ranges, the device ends at ", transpar%idxdevice(2), ', contact "',&
            & trim(transpar%contacts(ii)%name), '" is between ', transpar%contacts(ii)%idxrange(1),&
            & ' and ',transpar%contacts(ii)%idxrange(2)
        call error(trim(strTmp))
      end if
    end do

    ! Check for atom(s) occuring in multiple contacts
    do ii = 1, transpar%ncont
      do jj = 1, transpar%ncont
        if (ii == jj) then
          cycle
        end if
        tFailCheck = .false.
        if (transpar%contacts(ii)%idxrange(1) <= transpar%contacts(jj)%idxrange(1)) then
          if (transpar%contacts(ii)%idxrange(2) >= transpar%contacts(jj)%idxrange(1)) then
            tFailCheck = .true.
          end if
        else
          if (transpar%contacts(jj)%idxrange(2) >= transpar%contacts(ii)%idxrange(1)) then
            tFailCheck = .true.
          end if
        end if
        if (tFailCheck) then
          write(strTmp,"(A,A,A,A,A)")"Contact '",trim(transpar%contacts(ii)%name),"' and '",&
              & trim(transpar%contacts(jj)%name),"' share atoms"
          call error(trim(strTmp))
        end if
      end do
    end do

    ! check for additional atoms outside of the device and all contacts
    if (maxval(transpar%contacts(:)%idxrange(2)) < nAtom) then
      call error("Atoms present that are not in the device or any contact region")
    end if

    ! Check for gaps in atom ranges between regions
    allocate(notInRegion(nAtom))
    notInRegion = .true.
    notInRegion(transpar%idxdevice(1):transpar%idxdevice(2)) = .false.
    do ii = 1, transpar%ncont
      notInRegion(transpar%contacts(ii)%idxrange(1):transpar%contacts(ii)%idxrange(2)) = .false.
    end do
    if (any(notInRegion)) then
      call error("Atom(s) present that are not in any region of the device or contacts")
    end if

  end subroutine checkTransportRanges
#:endif


  !> Clean up things that did not automatically get removed by going out of scope
  subroutine destructProgramVariables(this)

    !> Instance
    class(TDftbPlusMain), intent(inout) :: this

    if (this%electronicSolver%isElsiSolver) then
      call TElsiSolver_final(this%electronicSolver%elsi)
    end if

    if (this%tProjEigenvecs) then
      call destruct(this%iOrbRegion)
      call destruct(this%regionLabels)
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
    type(TRanlux), allocatable, intent(out) :: randomInit

    !> Random generator for the actual thermostat.
    type(TRanlux), allocatable, intent(out) :: randomThermostat

    type(TRandomGenPool) :: randGenPool

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
  subroutine initSocket(this, env, socketInput)

    !> Instance
    class(TDftbPlusMain), intent(inout) :: this

    !> Environment settings.
    type(TEnvironment), intent(in) :: env

    !> Input data for the socket.
    type(ipiSocketCommInp), intent(inout) :: socketInput


    logical :: tDummy

    if (env%tGlobalLead) then
      write(stdOut, "(A,1X,A)") "Initialising for socket communication to host",&
          & trim(socketInput%host)
      this%socket = IpiSocketComm(socketInput)
    end if
    call receiveGeometryFromSocket(env, this%socket, this%tPeriodic, this%coord0, this%latVec,&
        & this%tCoordsChanged, this%tLatticeChanged, tDummy)

  end subroutine initSocket
#:endif

#:if WITH_TRANSPORT
  subroutine initTransport(this, env, input)

    !> Instance
    class(TDftbPlusMain), intent(inout) :: this

    !> Computational environment
    type(TEnvironment), intent(inout) :: env

    !> Input data
    type(TInputData), intent(in) :: input

    !> Electronic structure solver and its capabilities
    type(TElectronicSolver) :: electronicSolver

    logical :: tAtomsOutside
    integer :: iSpin
    integer :: nSpinChannels, iCont, jCont
    real(dp) :: mu1, mu2

    electronicSolver = this%electronicSolver

    ! contact calculation in case some contact is computed
    this%tContCalc = (input%transpar%taskContInd /= 0)

    if (this%nSpin <=2) then
      nSpinChannels = this%nSpin
    else
      nSpinChannels = 1
    endif

    associate(transpar=>input%transpar, greendens=>input%ginfo%greendens)
      ! Non colinear spin not yet supported
      ! Include the built-in potential as in negf init, but the whole
      ! scc only works for
      ! calculation without spin (poisson does not support spin dependent
      ! built in potentials)
      if (transpar%ncont > 0) then
        allocate(this%mu(transpar%ncont, nSpinChannels))
        this%mu = 0.0_dp
        do iSpin = 1, nSpinChannels
          this%mu(1:transpar%ncont, iSpin) =&
               & minval(transpar%contacts(1:transpar%ncont)%eFermi(iSpin))&
               & - transpar%contacts(1:transpar%ncont)%potential
        end do
        ! Test if this is a non-equilibrium system
        lpConts: do iCont = 1, transpar%ncont
          do jCont = iCont + 1, transpar%ncont
            do iSpin = 1, nSpinChannels
              mu1 = transpar%contacts(iCont)%eFermi(iSpin) - transpar%contacts(iCont)%potential
              mu2 = transpar%contacts(jCont)%eFermi(iSpin) - transpar%contacts(jCont)%potential
              if (abs(mu1 - mu2) > tolEfEquiv) then
                ! there is no global chemical potential so associated free energy currently
                ! undefined.
                electronicSolver%elecChemPotAvailable = .false.
                exit lpConts
              end if
            end do
          end do
        end do lpConts

      else
        allocate(this%mu(1, nSpinChannels))
        this%mu(1,1:nSpinChannels) = greendens%oneFermi(1:nSpinChannels)
      end if

    end associate

    if (this%tNegf) then
      write(stdOut,*) 'init negf'
      if (size(this%denseDesc%iAtomStart) /= this%nAtom+1) then
        call error('Internal error: DenseDesc not created')
      end if

      ! Some checks and initialization of GDFTB/NEGF
      call TNegfInt_init(this%negfInt, input%transpar, env, input%ginfo%greendens,&
          & input%ginfo%tundos, this%tempElec, electronicSolver%iSolver)

      this%ginfo = input%ginfo

      if (allocated(input%ctrl%indMovedAtom)) then
        tAtomsOutside = any(input%ctrl%indMovedAtom < input%transpar%idxdevice(1))&
            & .or. any(input%ctrl%indMovedAtom > input%transpar%idxdevice(2))
        if (tAtomsOutside) then
          call error("There are moving atoms specified outside of the device region")
        end if
      end if

      if (input%ctrl%tLatOpt) then
        call error("Lattice optimization is not currently possible with transport")
      end if

    end if

    this%transpar = input%transpar

    !Write Dos and tunneling on separate files?
    this%writeTunn = this%ginfo%tundos%writeTunn
    this%tWriteLDOS = this%ginfo%tundos%writeLDOS

    if (this%tWriteLDOS) then
      call move_alloc(this%ginfo%tundos%dosLabels, this%regionLabelLDOS)
    end if

  end subroutine initTransport

#:endif

  !> Initialises (clears) output files.
  subroutine initOutputFiles(this, env)

    !> Instance
    class(TDftbPlusMain), intent(inout) :: this

    !> Environment
    type(TEnvironment), intent(inout) :: env


    call TTaggedWriter_init(this%taggedWriter)

    if (this%tWriteAutotest) then
      call initOutputFile(autotestTag)
    end if
    if (this%tWriteResultsTag) then
      call initOutputFile(resultsTag)
    end if
    if (this%tWriteBandDat) then
      call initOutputFile(bandOut)
    end if
    if (this%tDerivs) then
      call initOutputFile(hessianOut)
    end if
    if (this%tWriteDetailedOut) then
      call initOutputFile(userOut, this%fdDetailedOut)
      call env%fileFinalizer%register(this%fdDetailedOut)
    end if
    if (this%tMD) then
      call initOutputFile(mdOut, this%fdMD)
      call env%fileFinalizer%register(this%fdMD)
    end if
    if (this%isGeoOpt .or. this%tMD) then
      call clearFile(trim(this%geoOutFile) // ".gen")
      call clearFile(trim(this%geoOutFile) // ".xyz")
    end if
    if (allocated(this%esp)) then
      call initOutputFile(this%esp%espOutFile)
    end if

  end subroutine initOutputFiles


  !> Allocates most of the large arrays needed during the DFTB run.
  subroutine initArrays(this, env)

    !> Instance
    class(TDftbPlusMain), intent(inout) :: this

    !> Current environment
    type(TEnvironment), intent(in) :: env


    logical :: isREKS
    integer :: nSpinHams, sqrHamSize, iDet

    isREKS = allocated(this%reks)

    if (isREKS) then
      allocate(this%rhoPrim(0, 1))
    else
      allocate(this%rhoPrim(0, this%nSpin))
    end if
    allocate(this%h0(0))
    if (this%tImHam) then
      allocate(this%iRhoPrim(0, this%nSpin))
    end if

    allocate(this%excitedDerivs(0,0))
    if (this%tForces) then
      if (.not.isREKS) then
        allocate(this%ERhoPrim(0))
      end if
      allocate(this%derivs(3, this%nAtom))
      if (this%deltaDftb%isNonAufbau) then
        allocate(this%groundDerivs(3, this%nAtom))
        allocate(this%tripletDerivs(3, this%nAtom))
        allocate(this%mixedDerivs(3, this%nAtom))
        if (this%tStress) then
          allocate(this%mixedStress(3,3))
          allocate(this%tripletStress(3,3))
        end if
      end if
      if (this%tExtChrg) then
        allocate(this%chrgForces(3, this%nExtChrg))
      end if
      if (this%tLinRespZVect) then
        deallocate(this%excitedDerivs)
        allocate(this%excitedDerivs(3, this%nAtom))
      end if
    end if

    call init(this%potential, this%orb, this%nAtom, this%nSpin)

    ! Nr. of independent spin Hamiltonians
    select case (this%nSpin)
    case (1)
      nSpinHams = 1
    case (2)
      nSpinHams = 2
    case (4)
      nSpinHams = 1
    end select
    if (isREKS) then
      nSpinHams = 1
    end if

    do iDet = 1, size(this%dftbEnergy)
      call TEnergies_init(this%dftbEnergy(iDet), this%nAtom, nSpinHams)
    end do

    sqrHamSize = this%denseDesc%fullSize

    if (this%electronicSolver%providesEigenvals) then
      allocate(this%eigen(sqrHamSize, this%nKPoint, nSpinHams))
      allocate(this%filling(sqrHamSize, this%nKPoint, nSpinHams))
    else
      ! due to use of the shape elsewhere in determining kpoints and spin channels:
      allocate(this%eigen(0, this%nKPoint, nSpinHams))
      allocate(this%filling(0, this%nKPoint, nSpinHams))
    end if
    this%eigen(:,:,:) = 0.0_dp
    this%filling(:,:,:) = 0.0_dp

    allocate(this%coord0Fold(3, this%nAtom))

    if (this%tMD) then
      allocate(this%newCoords(3, this%nAtom))
    end if

    if ((this%tMulliken .and. this%tSpinOrbit) .or. this%tImHam) then
      allocate(this%orbitalL(3, this%orb%mShell, this%nAtom))
    end if

    ! Decides whether large dense matricese should be allocated
    ! Currently needed by dense eigensolvers, hence not needed if
    ! 1. only H/S should be printed
    ! 2. Solver == GreensFunctions
    ! 3. Solver == TransportOnly
    ! 4. Solver == ELSI using a sparse solver
    this%tLargeDenseMatrices = .not. (this%tWriteRealHS .or. this%tWriteHS .or. &
          &   (this%electronicSolver%iSolver == electronicSolverTypes%GF) .or. &
          &   (this%electronicSolver%iSolver == electronicSolverTypes%OnlyTransport) )
    if (this%electronicSolver%isElsiSolver) then
      this%tLargeDenseMatrices = this%tLargeDenseMatrices&
          & .and. .not. this%electronicSolver%elsi%isSparse
    end if
    if (this%tLargeDenseMatrices) then
      call this%allocateDenseMatrices(env)
    end if

    if (this%isLinResp) then
      if (withMpi) then
        call error("Linear response calc. does not work with MPI yet")
      end if
      if (this%tLinRespZVect) then
        allocate(this%rhoSqrReal(sqrHamSize, sqrHamSize, this%nSpin))
      end if
    end if

    if (this%isLinResp .and. this%tPrintExcitedEigVecs) then
      allocate(this%occNatural(this%orb%nOrb))
    end if

    if (this%tMD) then
      allocate(this%velocities(3, this%nAtom))
      allocate(this%movedVelo(3, this%nMovedAtom))
      allocate(this%movedAccel(3, this%nMovedAtom))
      allocate(this%movedMass(3, this%nMovedAtom))
      this%movedMass(:,:) = spread(this%mass(this%indMovedAtom),1,3)
    end if

    if (this%tDipole) then
      if (this%deltaDftb%isSpinPurify) then
        allocate(this%dipoleMoment(3, this%deltaDftb%nDeterminant()+1))
      else
        allocate(this%dipoleMoment(3, this%deltaDftb%nDeterminant()))
      end if
    end if

  end subroutine initArrays


  !> Initialize storage for multi-determinantal calculations
  subroutine initDetArrays(this)

    !> Instance
    class(TDftbPlusMain), intent(inout) :: this


    this%nDets = this%deltaDftb%nDeterminant()
    if (this%nDets > 1) then
      ! must be SCC and also need storage for final charges
      allocate(this%qDets(this%orb%mOrb, this%nAtom, this%nSpin, this%nDets))
      this%qDets(:,:,:,:) = 0.0_dp
      ! When block charges are needed
      if (allocated(this%dftbU) .or. allocated(this%onSiteElements)) then
        allocate(this%qBlockDets(this%orb%mOrb, this%orb%mOrb, this%nAtom, this%nSpin, this%nDets))
        this%qBlockDets(:,:,:,:,:) = 0.0_dp
      end if
      if (this%isRangeSep) then
        allocate(this%deltaRhoDets(this%nOrb * this%nOrb * this%nSpin, this%nDets))
        this%deltaRhoDets(:,:) = 0.0_dp
      end if
    end if

  end subroutine initDetArrays


#:if WITH_TRANSPORT

  !> initialize arrays for tranpsport
  subroutine initTransportArrays(this, transpar)

    !> Instance
    class(TDftbPlusMain), intent(inout) :: this

    !> Transport parameters
    type(TTransPar), intent(inout) :: transpar


    !> Format for two values with units
    character(len=*), parameter :: format2U = "(1X,A, ':', T32, F18.10, T51, A, T54, F16.4, T71, A)"

    if (this%tUpload) then
      allocate(this%shiftPerLUp(this%orb%mShell, this%nAtom))
      allocate(this%chargeUp(this%orb%mOrb, this%nAtom, this%nSpin))
      !> When block charges and potentials are present
      if (allocated(this%qBlockIn)) then
        allocate(this%shiftBlockUp(this%orb%mOrb, this%orb%mOrb, this%nAtom, this%nSpin))
        allocate(this%blockUp(this%orb%mOrb, this%orb%mOrb, this%nAtom, this%nSpin))
      end if
      call readContactShifts(this%shiftPerLUp, this%chargeUp, transpar, this%orb,&
          & this%shiftBlockUp, this%blockUp)
    end if


  end subroutine initTransportArrays

#:endif


  !> Set up storage for dense matrices, either on a single processor, or as BLACS matrices
  subroutine allocateDenseMatrices(this, env)

    !> Instance
    class(TDftbPlusMain), intent(inout) :: this

    !> Computing environment
    type(TEnvironment), intent(in) :: env


    integer :: nLocalCols, nLocalRows, nLocalKS

    nLocalKS = size(this%parallelKS%localKS, dim=2)
  #:if WITH_SCALAPACK
    call scalafx_getlocalshape(env%blacs%orbitalGrid, this%denseDesc%blacsOrbSqr, nLocalRows,&
        & nLocalCols)
  #:else
    nLocalRows = this%denseDesc%fullSize
    nLocalCols = this%denseDesc%fullSize
  #:endif

    if (this%t2Component .or. .not. this%tRealHS) then
      allocate(this%HSqrCplx(nLocalRows, nLocalCols))
      allocate(this%SSqrCplx(nLocalRows, nLocalCols))
      allocate(this%eigVecsCplx(nLocalRows, nLocalCols, nLocalKS))
    else
      allocate(this%HSqrReal(nLocalRows, nLocalCols))
      allocate(this%SSqrReal(nLocalRows, nLocalCols))
      allocate(this%eigVecsReal(nLocalRows, nLocalCols, nLocalKS))
    end if

  end subroutine allocateDenseMatrices


#:if WITH_SCALAPACK
  #!
  #! SCALAPACK related routines
  #!

  !> Initialise parallel large matrix decomposition methods
  subroutine initScalapack(this, blacsOpts, nAtom, nOrb, t2Component, env)

    !> Instance
    class(TDftbPlusMain), intent(inout) :: this

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
  subroutine getDenseDescCommon(this)

    !> Instance
    class(TDftbPlusMain), intent(inout) :: this

    integer :: nOrb

    if (allocated(this%denseDesc)) then
      deallocate(this%denseDesc)
    end if
    allocate(this%denseDesc)
    allocate(this%denseDesc%iAtomStart(this%nAtom + 1))
    call buildSquaredAtomIndex(this%denseDesc%iAtomStart, this%orb)
    nOrb = this%denseDesc%iAtomStart(this%nAtom + 1) - 1
    this%denseDesc%t2Component = this%t2Component
    this%denseDesc%nOrb = nOrb
    if (this%t2Component) then
      this%denseDesc%fullSize = 2 * nOrb
    else
      this%denseDesc%fullSize = nOrb
    end if

  end subroutine getDenseDescCommon


#:if WITH_MBD
  !> Writes MBD-related info
  subroutine writeMbdInfo(input)

    !> MBD input parameters
    type(TDispMbdInp), intent(in) :: input

    character(lc) :: tmpStr

    write(stdOut, "(A)") ''
    select case (input%method)
    case ('ts')
      write(stdOut, "(A)") "Using TS dispersion corrections [Phys. Rev. B 80, 205414 (2009)]"
      write(stdOut, "(A)") "PLEASE CITE: J. Chem. Phys. 144, 151101 (2016)"
    case ('mbd-rsscs')
      write(stdOut,"(A)") "Using MBD dispersion corrections [Phys. Rev. Lett. 108, 236402 (2012)]"
      write(stdOut,"(A)") "PLEASE CITE: J. Chem. Phys. 144, 151101 (2016)"
    end select
    select case (trim(input%vdw_params_kind))
    case ('tssurf')
      write(tmpStr, "(A)") "vdw-surf [Phys. Rev. Lett. 108, 146103 (2012)] "
    case ('ts')
      write(tmpStr, "(A)") "vdw(TS) [Phys. Rev. Lett. 102, 073005 (2009)] "
    end select
    write(stdOut, "(A,T30,A)") "  Parameters", tmpStr
    select case (input%method)
    case ('ts')
      write(tmpStr,"(E18.6)") input%ts_f_acc
      write(stdOut, "(A,T30,A)") '  TSForceAccuracy', trim(adjustl(tmpStr))
      write(tmpStr, "(E18.6)") input%ts_ene_acc
      write(stdOut, "(A,T30,A)") '  TSEnergyAccuracy', trim(adjustl(tmpStr))
    case ('mbd-rsscs')
      write(tmpStr, "(3(I3,1X))") input%k_grid
      write(stdOut,"(A,T30,A)") "  MBD k-Grid", trim(adjustl(tmpStr))
      write(tmpStr, "(3(F4.3,1X))") input%k_grid_shift
      write(stdOut,"(A,T30,A)") "  MBD k-Grid shift", trim(adjustl(tmpStr))
      write(tmpStr, "(I3)") input%n_omega_grid
      write(stdOut, "(A,T30,A)") "  Gridsize (frequencies)", trim(adjustl(tmpStr))
    end select
    write(stdOut,"(A)") ""

  end subroutine writeMbdInfo
#:endif


  !> Check for compatibility between requested electronic solver and features of the calculation
  subroutine ensureSolverCompatibility(iSolver, tSpin, kPoints, parallelOpts, nIndepHam, tempElec)

    !> Solver number (see dftbp_elecsolvertypes)
    integer, intent(in) :: iSolver

    !> Is this a spin polarised calculation
    logical, intent(in) :: tSpin

    !> Set of k-points used in calculation (or [0,0,0] if molecular)
    real(dp), intent(in) :: kPoints(:,:)

    !> Options for a parallel calculation, if needed
    type(TParallelOpts), intent(in), allocatable :: parallelOpts

    !> Number of indepdent hamiltonian matrices at a given k-value
    integer, intent(in) :: nIndepHam

    !> Temperature of the electrons
    real(dp), intent(in) :: tempElec

    logical :: tElsiSolver
    integer :: nKPoint

    ! Temporary error test for PEXSI bug (July 2019)
    if (iSolver == electronicSolverTypes%pexsi .and. any(kPoints /= 0.0_dp)) then
      call warning("A temporary PEXSI bug may prevent correct evaluation at general k-points.&
          & This should be fixed soon.")
    end if

    tElsiSolver = any(iSolver ==&
        & [electronicSolverTypes%elpa, electronicSolverTypes%omm, electronicSolverTypes%pexsi,&
        & electronicSolverTypes%ntpoly, electronicSolverTypes%elpadm])
    if (.not. withELSI .and. tElsiSolver) then
      call error("This binary was not compiled with ELSI support enabled")
    end if

    nKPoint = size(kPoints, dim=2)
    if (withMpi) then
      if (tElsiSolver .and. parallelOpts%nGroup /= nIndepHam * nKPoint) then
        call error("This solver requires as many parallel processor groups as there are independent&
            & spin and k-point combinations")
      end if
    end if

    if (iSolver == electronicSolverTypes%pexsi .and. tempElec < epsilon(0.0)) then
      call error("This solver requires a finite electron broadening")
    end if

  end subroutine ensureSolverCompatibility


  !> Modify the reference atomic shell charges for the neutral atom
  subroutine applyCustomReferenceOccupations(customOccAtoms, customOccFillings, species, orb,&
      & referenceN0, q0)

    !> Array of occupation arrays, one for each atom
    type(TWrappedInt1), allocatable, intent(in) :: customOccAtoms(:)

    !> Reference fillings for atomic shells
    real(dp), intent(in) :: customOccFillings(:,:)

    !> Species of atoms
    integer, intent(in) :: species(:)

    !> Atomic orbital data
    type(TOrbitals), intent(in) :: orb

    !> Reference charges from the Slater-Koster file
    real(dp), intent(in) :: referenceN0(:,:)

    !> Charges required for atomic neutrality in reference state
    real(dp), intent(inout) :: q0(:,:,:)

    integer :: nCustomBlock, iCustomBlock, iCustomAtom, nAtom, iAt, iSp
    real(dp), allocatable :: refOcc(:,:)

    nAtom = size(species)
    ! note that all arrays, referenceN0, customOccAtoms, refOcc
    ! are allocated to orb%mShell so assignments vecA(:,) = vecB(:,) work
    allocate(refOcc(orb%mShell, nAtom))
    ! initialize to referenceN0
    do iAt = 1, nAtom
      iSp = species(iAt)
      refOcc(:, iAt) = referenceN0(:, iSp)
    end do

    ! override to customOccupation
    if (allocated(customOccAtoms)) then
      nCustomBlock = size(customOccAtoms)
      do iCustomBlock = 1, nCustomBlock
        do iCustomAtom = 1, size(customOccAtoms(iCustomBlock)%data)
          iAt =  customOccAtoms(iCustomBlock)%data(iCustomAtom)
          refOcc(:, iAt) = customOccFillings(:,iCustomBlock)
        end do
      end do
    end if

    ! initialize q0 with right orbital order
    call initQFromUsrChrg(q0, refOcc, species, orb)

  end subroutine applyCustomReferenceOccupations


  !> Print out the reference occupations for atoms
  subroutine printCustomReferenceOccupations(orb, species, customOccAtoms, customOccFillings)

    !> Atomic orbital information
    type(TOrbitals), intent(in) :: orb

    !> Chemical species of atoms
    integer, intent(in) :: species(:)

    !> Array of occupation arrays, one for each atom
    type(TWrappedInt1), intent(in) :: customOccAtoms(:)

    !> Fillings for each atomic shell
    real(dp), intent(in) :: customOccFillings(:,:)

    character(lc) :: formstr, outStr
    integer :: nCustomBlock, iCustomBlock, iSp, nShell, nAtom, iSh
    character(sc), allocatable :: shellnames(:)

    nCustomBlock = size(customOccFillings)
    if (nCustomBlock == 0) then
      return
    end if
    write(stdout, "(A)") "Custom defined reference occupations:"
    do iCustomBlock = 1, size(customOccAtoms)
      nAtom = size(customOccAtoms(iCustomBlock)%data)
      if (nAtom == 1) then
        write(outStr, "(A)") "Atom:"
      else
        write(outStr, "(A)") "Atoms:"
      end if
      write(formstr, "(I0,A)") nAtom, "(I0,1X))"
      write(stdout, "(A,T30,"//trim(formstr)//")") trim(outStr), customOccAtoms(iCustomBlock)%data
      iSp = species(customOccAtoms(iCustomBlock)%data(1))
      nShell = orb%nShell(iSp)
      call getShellNames(iSp, orb, shellnames)
      outStr = ""
      do iSh = 1, nShell
        if (iSh > 1) then
          write(outStr,"(A,',')")trim(outStr)
        end if
        write(outStr,"(A,1X,A,F8.4)")trim(outStr), trim(shellnames(iSh)),&
            & customOccFillings(iSh, iCustomBlock)
      end do
      write(stdout,"(A,T29,A)")"Fillings:",trim(outStr)
      deallocate(shellnames)
    end do
  end subroutine printCustomReferenceOccupations


  !> Initialises SCC related parameters before geometry loop starts
  function getMinSccIters(tSccCalc, isDftbU, nSpin) result(minSccIter)

    !> Is this a self consistent calculation
    logical, intent(in) :: tSccCalc

    !> Are there orbital potentials present
    logical, intent(in) :: isDftbU

    !> Number of spin channels
    integer, intent(in) :: nSpin

    !> Minimum possible number of self consistent iterations
    integer :: minSccIter

    if (tSccCalc) then
      if (isDftbU) then
        minSccIter = 2
      else
        if (nSpin == 1) then
          minSccIter = 1
        else
          minSccIter = 2
        end if
      end if
    else
      minSccIter = 1
    end if

  end function getMinSccIters


  !> Stop if any range separated incompatible setting is found
  subroutine ensureRangeSeparatedReqs(this, tShellResolved, rangeSepInp)

    !> Instance
    class(TDftbPlusMain), intent(inout) :: this

    !> Is this a shell resolved calculation
    logical, intent(in) :: tShellResolved

    !> Parameters for the range separated calculation
    type(TRangeSepInp), intent(in) :: rangeSepInp

    if (withMpi) then
      call error("Range separated calculations do not work with MPI yet")
    end if

    if (this%tPeriodic) then
      call error("Range separated functionality only works with non-periodic structures at the&
          & moment")
    end if

    if (this%tHelical) then
      call error("Range separated functionality only works with non-helical structures at the&
          & moment")
    end if

    if (this%tReadChrg .and. rangeSepInp%rangeSepAlg == rangeSepTypes%threshold) then
      call error("Restart on thresholded range separation not working correctly")
    end if

    if (tShellResolved) then
      call error("Range separated functionality currently does not yet support shell-resolved scc")
    end if

    if (this%tAtomicEnergy) then
      call error("Atomic resolved energies cannot be calculated with the range-separated&
          & hybrid functional at the moment")
    end if

    if (this%nSpin > 2) then
      call error("Range separated calculations not implemented for non-colinear calculations")
    end if

    if (this%tSpinOrbit) then
      call error("Range separated calculations not currently implemented for spin orbit")
    end if

    if (this%isXlbomd) then
      call error("Range separated calculations not currently implemented for XLBOMD")
    end if

    if (this%t3rd) then
      call error("Range separated calculations not currently implemented for 3rd order DFTB")
    end if

    if (allocated(this%dftbU)) then
      call error("Range separated calculations not currently implemented for DFTB+U")
    end if

    if (this%isRS_LinResp) then

      if (allocated(this%onSiteElements)) then
        call error("Excited state range separated calculations not implemented for onsite&
            & corrections")
      end if

      if (this%nSpin > 1) then
        call error("Excited state range separated calculations not implemented for spin polarized&
            & calculations")
      end if

      if (this%linearResponse%symmetry /= "S") then
        call error("Excited state range separated calculations currently only implemented for&
            & singlet excitaions")
      end if

    end if

  end subroutine ensureRangeSeparatedReqs


  !> Stop if linear response module can not be invoked due to unimplemented combinations of
  !> features.
  subroutine ensureLinRespConditions(tSccCalc, t3rd, tRealHS, tPeriodic, tCasidaForces, solvation,&
      & isRS_LinResp, nSpin, tSpin, tHelical, tSpinOrbit, isDftbU, tempElec, input)

    !> Is the calculation SCC?
    logical, intent(in) :: tSccCalc

    !> 3rd order hamiltonian contributions included
    logical, intent(in) :: t3rd

    !> a real hamiltonian
    logical, intent(in) :: tRealHs

    !> periodic boundary conditions
    logical, intent(in) :: tPeriodic

    !> forces being evaluated in the excited state
    logical, intent(in) :: tCasidaForces

    !> Solvation data and calculations
    class(TSolvation), allocatable :: solvation

    !> Is this an excited state calculation with range separation
    logical, intent(in) :: isRS_LinResp

    !> Number of spin components, 1 is unpolarised, 2 is polarised, 4 is noncolinear / spin-orbit
    integer, intent(in) :: nSpin

    !> Is this a spin polarised calculation
    logical, intent(in) :: tSpin

    !> If the calculation is helical geometry
    logical :: tHelical

    !> Is there spin-orbit coupling
    logical, intent(in) :: tSpinOrbit

    !> Is this a DFTB+U calculation?
    logical, intent(in) :: isDftbU

    !> Temperature of the electrons
    real(dp), intent(in) :: tempElec

    !> Holds the parsed input data.
    type(TInputData), intent(in), target :: input

    character(lc) :: tmpStr

    if (withMpi) then
      call error("Linear response calc. does not work with MPI yet")
    end if

    if (.not. tSccCalc) then
      call error("Linear response excitation requires SCC=Yes")
    end if

    if (t3rd) then
      call error("Third order currently incompatible with excited state")
    end if
    if (.not. tRealHS) then
      call error("Only real systems are supported for excited state calculations")
    end if
    if ((tPeriodic .or. tHelical) .and. .not.tRealHS) then
      call error("Linear response only works with non-periodic or gamma-point molecular crystals")
    end if
    if (tPeriodic .and. tCasidaForces) then
      call error("Forces in the excited state for periodic geometries are currently unavailable")
    end if

    if (allocated(solvation)) then
      call error("Solvation models do not work with linear response yet.")
    end if

    if (nspin > 2) then
      call error("Linear reponse does not work with non-colinear spin polarization yet")
    end if

    if (tSpinOrbit) then
      call error("Linear response does not support spin orbit coupling at the moment.")
    end if
    if (isDftbU) then
      call error("Linear response does not support LDA+U yet")
    end if
    if (input%ctrl%tShellResolved) then
      call error("Linear response does not support shell resolved scc yet")
    end if

    if (tempElec > minTemp .and. tCasidaForces) then
      write(tmpStr, "(A,E12.4,A)")"Excited state forces are not implemented yet for fractional&
          & occupations, kT=", tempElec/Boltzmann,"K"
      call warning(tmpStr)
    end if

    if (input%ctrl%lrespini%nstat == 0) then
      if (tCasidaForces) then
        call error("Excited forces only available for StateOfInterest non zero.")
      end if
      if (input%ctrl%lrespini%tPrintEigVecs .or. input%ctrl%lrespini%tCoeffs) then
        call error("Excited eigenvectors only available for StateOfInterest non zero.")
      end if
    end if

    if (isRS_LinResp) then
      if (tPeriodic) then
        call error("Range separated excited states for periodic geometries are currently&
            & unavailable")
      end if
      if (nSpin > 1) then
        call error("Range separated excited states for spin polarized calculations are currently&
            & unavailable")
      end if
    else
      if (input%ctrl%lrespini%energyWindow < 0.0_dp) then
        call error("Negative energy window for excitations")
      end if
    end if

  end subroutine ensureLinRespConditions


  !> Determine range separated cut-off and also update maximal cutoff
  subroutine getRangeSeparatedCutOff(cutoffRed, cutOff)

    !> Reduction in cut-off
    real(dp), intent(in) :: cutoffRed

    !> Resulting cut-off
    type(TCutoffs), intent(inout) :: cutOff

    cutOff%lcCutOff = 0.0_dp
    if (cutoffRed < 0.0_dp) then
      call error("Cutoff reduction for range-separated neighbours should be zero or positive.")
    end if
    cutOff%lcCutOff = cutOff%skCutOff - cutoffRed
    if (cutOff%lcCutOff < 0.0_dp) then
      call error("Screening cutoff for range-separated neighbours too short.")
    end if
    cutOff%mCutoff = max(cutOff%mCutOff, cutoff%lcCutOff)

  end subroutine getRangeSeparatedCutOff


  !> Initialise range separated extension.
  subroutine initRangeSeparated(this, nAtom, species0, hubbU, rangeSepInp, tSpin, isREKS, rangeSep,&
      & deltaRhoIn, deltaRhoOut, deltaRhoDiff, deltaRhoInSqr, deltaRhoOutSqr, nMixElements)

    !> Instance
    class(TDftbPlusMain), intent(inout) :: this

    !> Number of atoms in the system
    integer, intent(in) :: nAtom

    !> species of atoms
    integer, intent(in) :: species0(:)

    !> Hubbard values for species
    real(dp), intent(in) :: hubbU(:,:)

    !> input for range separated calculation
    type(TRangeSepInp), intent(in) :: rangeSepInp

    !> Is this spin restricted (F) or unrestricted (T)
    logical, intent(in) :: tSpin

    !> Is this DFTB/SSR formalism
    logical, intent(in) :: isREKS

    !> Resulting settings for range separation
    type(TRangeSepFunc), allocatable, intent(out) :: rangeSep

    !> Change in input density matrix flattened to 1D array
    real(dp), allocatable, target, intent(out) :: deltaRhoIn(:)

    !> Change in output density matrix flattened to 1D array
    real(dp), allocatable, target, intent(out) :: deltaRhoOut(:)

    !> Change in density matrix between in and out
    real(dp), allocatable, intent(out) :: deltaRhoDiff(:)

    !> Change in input density matrix
    real(dp), pointer, intent(out) :: deltaRhoInSqr(:,:,:)

    !> Change in output density matrix
    real(dp), pointer, intent(out) :: deltaRhoOutSqr(:,:,:)

    !> Number of mixer elements
    integer, intent(out) :: nMixElements

    allocate(rangeSep)
    call RangeSepFunc_init(rangeSep, nAtom, species0, hubbU(1,:), rangeSepInp%screeningThreshold,&
        & rangeSepInp%omega, tSpin, isREKS, rangeSepInp%rangeSepAlg)
    allocate(deltaRhoIn(this%nOrb * this%nOrb * this%nSpin))
    allocate(deltaRhoOut(this%nOrb * this%nOrb * this%nSpin))
    allocate(deltaRhoDiff(this%nOrb * this%nOrb * this%nSpin))
    deltaRhoInSqr(1:this%nOrb, 1:this%nOrb, 1:this%nSpin) =>&
        & deltaRhoIn(1 : this%nOrb * this%nOrb * this%nSpin)
    deltaRhoOutSqr(1:this%nOrb, 1:this%nOrb, 1:this%nSpin) =>&
        & deltaRhoOut(1 : this%nOrb * this%nOrb * this%nSpin)
    nMixElements = this%nOrb * this%nOrb * this%nSpin
    deltaRhoInSqr(:,:,:) = 0.0_dp

  end subroutine initRangeSeparated


  !> Initializes PLUMED calculator.
  subroutine initPlumed(this, env, tPlumed, tMD, plumedCalc)

    !> Instance
    class(TDftbPlusMain), intent(inout) :: this

    !> Environment settings
    type(TEnvironment), intent(in) :: env

    !> Whether plumed should be used
    logical, intent(in) :: tPlumed

    !> Whether this is an MD-run
    logical, intent(in) :: tMD

    !> Plumed calculator (allocated only on demand)
    type(TPlumedCalc), allocatable, intent(out) :: plumedCalc


    ! Minimal plumed API version (as in Plumed 2.5.3)
    ! Earlier versions may work but were not tested
    integer, parameter :: minApiVersion = 6

    integer :: apiVersion
    character(300) :: strTmp

    if (.not. tPlumed) then
      return
    end if
    if (.not. withPlumed) then
      call error("Code was compiled without PLUMED support")
    end if
    if (.not. tMD) then
      call error("Metadynamics via PLUMED is only possible in MD-simulations")
    end if
    allocate(plumedCalc)
    call TPlumedCalc_init(plumedCalc)
    call plumedCalc%sendCmdPtr("getApiVersion", apiVersion)
    if (apiVersion < minApiVersion) then
      write(strTmp, "(A,I0,A)") "PLUMED interface has not been tested with PLUMED API version < ",&
          & minApiVersion, ". Your PLUMED library provides API version ", apiVersion, ". Check your&
          & results carefully and consider to use a more recent PLUMED library if in doubt!"
      call warning(strTmp)
    end if
    call plumedCalc%sendCmdVal("setNatoms", this%nAtom)
    call plumedCalc%sendCmdVal("setPlumedDat", "plumed.dat")
    call plumedCalc%sendCmdVal("setNoVirial", 0)
    call plumedCalc%sendCmdVal("setTimestep", this%deltaT)
    call plumedCalc%sendCmdVal("setMDEnergyUnits", Hartree__kJ_mol)
    call plumedCalc%sendCmdVal("setMDLengthUnits", Bohr__nm)
    call plumedCalc%sendCmdVal("setMDTimeUnits", au__ps)
    #:if WITH_MPI
      call plumedCalc%sendCmdVal("setMPIFComm", env%mpi%globalComm%id)
    #:endif
    call plumedCalc%sendCmdVal("init", 0)

  end subroutine initPlumed


  subroutine checkReksConsistency(reksInp, solvation, onSiteElements, kPoint, nEl, nKPoint,&
      & tSccCalc, tSpin, tSpinOrbit, isDftbU, tEField, isLinResp, tPeriodic, tLatOpt, tReadChrg,&
      & tPoisson, isShellResolved)

    !> data type for REKS input
    type(TReksInp), intent(in) :: reksInp

    !> Solvation data and calculations
    class(TSolvation), allocatable, intent(in) :: solvation

    !> Correction to energy from on-site matrix elements
    real(dp), allocatable, intent(in) :: onSiteElements(:,:,:,:)

    !> K-points
    real(dp), intent(in) :: kPoint(:,:)

    !> nr. of electrons
    real(dp), intent(in) :: nEl(:)

    !> nr. of K-points
    integer, intent(in) :: nKPoint

    !> Is the calculation SCC?
    logical, intent(in) :: tSccCalc

    !> is this a spin polarized calculation?
    logical, intent(in) :: tSpin

    !> is there spin-orbit coupling
    logical, intent(in) :: tSpinOrbit

    !> is this a DFTB+U calculation?
    logical, intent(in) :: isDftbU

    !> external electric field
    logical, intent(in) :: tEField

    !> Calculate Casida linear response excitations
    logical, intent(in) :: isLinResp

    !> if calculation is periodic
    logical, intent(in) :: tPeriodic

    !> optimize lattice constants?
    logical, intent(in) :: tLatOpt

    !> If initial charges/dens mtx. from external file.
    logical, intent(in) :: tReadChrg

    !> Whether Poisson solver is invoked
    logical, intent(in) :: tPoisson

    !> l-shell resolved SCC
    logical, intent(in) :: isShellResolved

    if (.not. tSccCalc) then
      call error("REKS requires SCC=Yes")
    end if
    if (tSpin) then
      call error("REKS is not compatible with standard DFTB spin polarization, only the&
          & SpinConstants block is required")
    end if

    if (tSpinOrbit) then
      call error("REKS is not compatible with spin-orbit (LS-coupling) calculation")
    else if (isDftbU) then
      call error("REKS is not compatible with DFTB+U calculation")
    else if (tEField) then
      call error("REKS is not compatible with external electric field, only point charge&
          & embedding is implemented")
    else if (isLinResp) then
      call error("REKS is not compatible with standard linear response excitation")
    else if (allocated(onSiteElements)) then
      call error("REKS is not compatible with onsite corrections")
    end if

    if (allocated(solvation)) then
      call error("REKS is currently not available with solvation")
    else if (tPoisson) then
      call error("Poisson solver is not compatible with REKS")
    elseif (isShellResolved) then
      call error("REKS does not support shell resolved scc yet")
    end if

    if (tPeriodic) then
      if ( .not. (nKPoint == 1 .and. all(kPoint(:, 1) == [0.0_dp, 0.0_dp, 0.0_dp])) ) then
        call error("REKS can compute only gamma-point in periodic case")
      end if
    end if

    if (reksInp%Efunction /= 1 .and. tLatOpt) then
      call error("Lattice optimization is only possible&
          & with single-state REKS, not SA-REKS or SI-SA-REKS")
    end if

    if (tReadChrg) then
      call error("Reading of initial charges is currently incompatible with REKS calculations")
    end if

    ! REKS can treat only closed shell systems.
    if (mod(nint(nEl(1)),2) /= 0) then
      call error("Current system is not a closed shell system, please check charge if using REKS")
    end if
    if (abs(nint(nEl(1)) - nEl(1)) >= elecTolMax) then
      call error("Current system is fractionally charged, please check charge if using REKS")
    end if

  end subroutine checkReksConsistency


  subroutine TReksCalc_init(reks, reksInp, electronicSolver, orb, spinW, nEl, extChrg, blurWidths,&
      & hamiltonianType, nSpin, nExtChrg, is3rd, isRangeSep, tForces, tPeriodic, tStress, tDipole)

    !> data type for REKS
    type(TReksCalc), intent(out) :: reks

    !> data type for REKS input
    type(TReksInp), intent(inout) :: reksInp

    !> electronic solver for the system
    type(TElectronicSolver), intent(in) :: electronicSolver

    !> Atomic orbital information
    type(TOrbitals), intent(in) :: orb

    !> Spin W values
    real(dp), intent(inout) :: spinW(:,:,:)

    !> nr. of electrons
    real(dp), intent(in) :: nEl(:)

    !> coordinates and charges of external point charges
    real(dp), allocatable, intent(in) :: extChrg(:,:)

    !> Width of the Gaussians if the charges are blurred
    real(dp), allocatable, intent(in) :: blurWidths(:)

    !> Hamiltonian type
    integer, intent(in) :: hamiltonianType

    !> Number of spin components, 1 is unpolarised, 2 is polarised, 4 is noncolinear / spin-orbit
    integer, intent(inout) :: nSpin

    !> Nr. of external charges
    integer, intent(in) :: nExtChrg

    !> Third order DFTB
    logical, intent(in) :: is3rd

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

    ! Condition for Hamiltonian types
    select case(hamiltonianType)
    case default
      call error("Invalid Hamiltonian")
    case(hamiltonianTypes%dftb)

      ! Condition for electronicSolver
      select case (electronicSolver%iSolver)
      case (electronicSolverTypes%GF)
        call error("REKS is not compatible with Green's function solver")
      case (electronicSolverTypes%onlyTransport)
        call error("REKS is not compatible with OnlyTransport-solver")
      case(electronicSolverTypes%qr, electronicSolverTypes%divideandconquer,&
          & electronicSolverTypes%relativelyrobust)
        call REKS_init(reks, reksInp, orb, spinW, nSpin, nEl(1), nExtChrg, extChrg,&
            & blurWidths, is3rd, isRangeSep, tForces, tPeriodic, tStress, tDipole)
      case(electronicSolverTypes%magma_gvd)
        call error("REKS is not compatible with MAGMA GPU solver")
      case(electronicSolverTypes%omm, electronicSolverTypes%pexsi, electronicSolverTypes%ntpoly,&
          & electronicSolverTypes%elpadm, electronicSolverTypes%elpa)
        call error("REKS is not compatible with ELSI-solvers")
      end select

    case(hamiltonianTypes%xtb)
      call error("xTB calculation currently not supported for REKS")
    end select

  end subroutine TReksCalc_init


  subroutine printReksInitInfo(reks, orb, speciesName, nType)

    !> data type for REKS
    type(TReksCalc), intent(in) :: reks

    !> data structure with atomic orbital information
    type(TOrbitals), intent(in) :: orb

    !> labels of atomic species
    character(mc), intent(in) :: speciesName(:)

    !> nr of different types (nAtom)
    integer, intent(in) :: nType

    integer :: ii, iType
    character(lc) :: strTmp

    write (stdOut,*)
    write (stdOut,*)
    write (stdOut, "(A,':',T30,A)") "REKS Calcuation", "Yes"

    select case (reks%reksAlg)
    case (reksTypes%noReks)
    case (reksTypes%ssr22)
      write (stdOut, "(A,':',T30,A)") "SSR(2,2) Calcuation", "Yes"
      if (reks%Efunction == 1) then
        write (stdOut, "(A,':',T30,A)") "Energy Functional", "PPS"
      else if (reks%Efunction == 2) then
        write (stdOut, "(A,':',T30,A)") "Energy Functional", "(PPS+OSS)/2"
      end if
    case (reksTypes%ssr44)
      call error("SSR(4,4) is not implemented yet")
    end select

    write (stdOut, "(A,':',T30,I14)") "Number of Core Orbitals", reks%Nc
    write (stdOut, "(A,':',T30,I14)") "Number of Active Orbitals", reks%Na
    write (stdOut, "(A,':',T30,I14)") "Number of Basis", orb%nOrb
    write (stdOut, "(A,':',T30,I14)") "Number of States", reks%nstates
    do ii = 1, reks%SAstates
      if (ii == 1) then
        write (strTmp, "(A,':')") "State-Averaging Weight"
      else
        write (strTmp, "(A)") ""
      end if
      write (stdOut, "(A,T30,F12.6)") trim(strTmp), reks%SAweight(ii)
    end do
    write (stdOut, "(A,':',T30,I14)") "State of Interest", reks%rstate

    if (reks%tReadMO) then
      write (stdOut, "(A,':',T30,A)") "Initial Guess", "Read Eigenvec.bin file"
    else
      write (stdOut, "(A,':',T30,A)") "Initial Guess", "Diagonalize H0 matrix"
    end if

    write (stdOut, "(A,':',T30,A)") "Newton-Raphson for FON opt", "Yes"
    write (stdOut, "(A,':',T30,I14)") "NR max. Iterations", reks%FonMaxIter
    if (reks%shift > epsilon(1.0_dp)) then
      write (stdOut, "(A,':',T30,A)") "Level Shifting", "Yes"
    else
      write (stdOut, "(A,':',T30,A)") "Level Shifting", "No"
    end if
    write (stdOut, "(A,':',T30,F12.6)") "Shift Value", reks%shift

    do iType = 1, nType
      if (iType == 1) then
        write (strTmp, "(A,':')") "W Scale Factor"
      else
        write (strTmp, "(A)") ""
      end if
      write (stdOut, "(A,T30,A3,'=',F12.6)") trim(strTmp), &
          & speciesName(iType), reks%Tuning(iType)
    end do

    if (reks%tTDP) then
      write (stdOut, "(A,':',T30,A)") "Transition Dipole", "Yes"
    else
      write (stdOut, "(A,':',T30,A)") "Transition Dipole", "No"
    end if

    if (reks%tForces) then

      if (reks%Lstate > 0) then
        write (stdOut, "(A,':',T30,A)") "Gradient of Microstate", "Yes"
        write (stdOut, "(A,':',T30,I14)") "Index of Interest", reks%Lstate
      else
        write (stdOut, "(A,':',T30,A)") "Gradient of Microstate", "No"
      end if

      if (reks%Efunction /= 1) then
        if (reks%Glevel == 1) then
          write (stdOut, "(A,':',T30,A)") "CP-REKS Solver", "Preconditioned Conjugate-Gradient"
          write (stdOut, "(A,':',T30,I14)") "CG max. Iterations", reks%CGmaxIter
          write (stdOut, "(A,':',T30,E14.6)") "CG Tolerance", reks%Glimit
          if (reks%tSaveMem) then
            write (stdOut, "(A,':',T30,A)") "Memory for A and Hxc", "Save in Cache Memory"
          else
            write (stdOut, "(A,':',T30,A)") "Memory for A and Hxc", "Direct Updating Without Saving"
          end if
        else if (reks%Glevel == 2) then
          write (stdOut, "(A,':',T30,A)") "CP-REKS Solver", "Conjugate-Gradient"
          write (stdOut, "(A,':',T30,I14)") "CG max. Iterations", reks%CGmaxIter
          write (stdOut, "(A,':',T30,E14.6)") "CG Tolerance", reks%Glimit
          if (reks%tSaveMem) then
            write (stdOut, "(A,':',T30,A)") "Memory for A and Hxc", "Save in Cache Memory"
          else
            write (stdOut, "(A,':',T30,A)") "Memory for A and Hxc", "Direct Updating Without Saving"
          end if
        else if (reks%Glevel == 3) then
          write (stdOut, "(A,':',T30,A)") "CP-REKS Solver", "Direct Matrix Multiplication"
        end if
        if (reks%tNAC) then
          write (stdOut, "(A,':',T30,A)") "Non-Adiabatic Coupling", "Yes"
        end if
      end if

      if (reks%tRD) then
        write (stdOut, "(A,':',T30,A)") "Relaxed Density for QM/MM", "Yes"
      end if

    end if

  end subroutine printReksInitInfo


end module dftbp_initprogram
