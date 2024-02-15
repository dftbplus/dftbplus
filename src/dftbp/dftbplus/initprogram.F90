!--------------------------------------------------------------------------------------------------!
!  DFTB+: general package for performing fast atomistic simulations                                !
!  Copyright (C) 2006 - 2023  DFTB+ developers group                                               !
!                                                                                                  !
!  See the LICENSE file for terms of usage and distribution.                                       !
!--------------------------------------------------------------------------------------------------!

#:include 'common.fypp'
#:include 'error.fypp'

!> Global variables and initialization for the main program.
module dftbp_dftbplus_initprogram
  use dftbp_common_accuracy, only : dp, lc, mc, sc, elecTolMax, minTemp, tolSameDist, tolEfEquiv
  use dftbp_common_atomicmass, only : getAtomicMass
  use dftbp_common_coherence, only : checkToleranceCoherence, checkExactCoherence
  use dftbp_common_constants, only : shellNames, Hartree__eV, Bohr__AA, amu__au, pi, au__ps,&
      & Bohr__nm, Hartree__kJ_mol, Boltzmann
  use dftbp_common_envcheck, only : checkStackSize
  use dftbp_common_environment, only : TEnvironment, globalTimers
  use dftbp_common_file, only : TFileDescr, setDefaultBinaryAccess, clearFile
  use dftbp_common_globalenv, only : stdOut, withMpi
  use dftbp_common_hamiltoniantypes, only : hamiltonianTypes
  use dftbp_common_status, only : TStatus
  use dftbp_dftb_densitymatrix, only : TDensityMatrix
  use dftbp_derivs_numderivs2, only : TNumDerivs, create
  use dftbp_derivs_perturb, only : TResponse, TResponse_init, responseSolverTypes
  use dftbp_dftb_blockpothelper, only : appendBlockReduced
  use dftbp_dftb_boundarycond, only : boundaryConditions, TBoundaryConditions,&
      & TBoundaryConditions_init
  use dftbp_dftb_coulomb, only : TCoulombInput
  use dftbp_dftb_dense, only : buildSquaredAtomIndex
  use dftbp_dftb_determinants, only : TDftbDeterminants, TDftbDeterminants_init
  use dftbp_dftb_dftbplusu, only : TDftbU, TDftbU_init
  use dftbp_dftb_dispdftd4, only : writeDftD4Info
  use dftbp_dftb_dispersions, only : TDispersionIface, TDispSlaKirk, TDispUFF, TSimpleDftD3,&
      & TDispDftD4, init, DispSlaKirk_init, DispUff_init
  use dftbp_dftb_elstatpot, only : TElStatPotentials, TElStatPotentials_init
  use dftbp_dftb_energytypes, only : TEnergies, TEnergies_init
  use dftbp_dftb_etemp, only : fillingTypes
  use dftbp_dftb_extfields, only : TEField
  use dftbp_dftb_h5correction, only : TH5CorrectionInput
  use dftbp_dftb_halogenx, only : THalogenX, THalogenX_init
  use dftbp_dftb_hamiltonian, only : TRefExtPot
  use dftbp_dftb_nonscc, only : TNonSccDiff, NonSccDiff_init, diffTypes
  use dftbp_dftb_onsitecorrection, only : Ons_getOrbitalEquiv, Ons_blockIndx
  use dftbp_dftb_orbitalequiv, only : OrbitalEquiv_merge, OrbitalEquiv_reduce
  use dftbp_dftb_periodic, only : TNeighbourList, TNeighbourlist_init, TSymNeighbourList,&
      & getCellTranslations
  use dftbp_dftb_pmlocalisation, only : TPipekMezey, initialise
  use dftbp_dftb_potentials, only : TPotentials, TPotentials_init
  use dftbp_dftb_hybridxc, only : THybridXcFunc, hybridXcAlgo, hybridXcFunc, hybridXcGammaTypes,&
      & THybridXcFunc_init, checkSupercellFoldingMatrix
  use dftbp_dftb_repulsive_chimesrep, only : TChimesRepInp, TChimesRep, TChimesRep_init
  use dftbp_dftb_repulsive_pairrepulsive, only : TPairRepulsiveItem
  use dftbp_dftb_repulsive_repulsive, only : TRepulsive
  use dftbp_dftb_repulsive_repulsivecont, only : TRepulsiveCont, TRepulsiveCont_init
  use dftbp_dftb_repulsive_repulsivelist, only : TRepulsiveList
  use dftbp_dftb_repulsive_twobodyrep, only : TTwoBodyRepInp, TTwoBodyRep, TTwoBodyRep_init
  use dftbp_dftb_rshgamma, only : getCoulombTruncationCutoff
  use dftbp_dftb_scc, only : TSccInput, TScc, TScc_init
  use dftbp_dftb_sccinit, only : initQFromFile, initQFromUsrChrg, initQFromAtomChrg,&
      & initQFromShellChrg
  use dftbp_dftb_shortgamma, only : TShortGammaInput, TShortGammaDamp
  use dftbp_dftb_slakocont, only : TSlakoCont, getCutOff
  use dftbp_dftb_spin, only: Spin_getOrbitalEquiv, ud2qm, qm2ud
  use dftbp_dftb_thirdorder, only : TThirdOrderInp, TThirdOrder, ThirdOrder_init
  use dftbp_dftb_uniquehubbard, only : TUniqueHubbard, TUniqueHubbard_init
  use dftbp_dftb_elecconstraints, only : TElecConstraint, TElecConstraint_init, TElecConstraintInput
  use dftbp_dftbplus_elstattypes, only : elstatTypes
  use dftbp_dftbplus_forcetypes, only : forceTypes
  use dftbp_dftbplus_inputdata, only : TParallelOpts, TInputData, THybridXcInp, TControl, TBlacsOpts
  use dftbp_dftbplus_outputfiles, only : autotestTag, bandOut, derivEBandOut, hessianOut, mdOut,&
      & resultsTag, userOut, fCharges, fStopDriver, fStopSCC
  use dftbp_dftbplus_qdepextpotproxy, only : TQDepExtPotProxy
  use dftbp_dftbplus_transportio, only : readContactShifts
  use dftbp_elecsolvers_elecsolvers, only : TElectronicSolver, electronicSolverTypes,&
      & TElectronicSolver_init
  use dftbp_elecsolvers_elsisolver, only : TElsiSolver_init, TElsiSolver_final
  use dftbp_extlibs_arpack, only : withArpack
  use dftbp_extlibs_elsiiface, only : withELSI
  use dftbp_extlibs_plumed, only : withPlumed, TPlumedCalc, TPlumedCalc_init
  use dftbp_extlibs_poisson, only : TPoissonInput
  use dftbp_extlibs_sdftd3, only : TSDFTD3, TSDFTD3_init, writeSDFTD3Info
  use dftbp_extlibs_tblite, only : TTBLite, TTBLite_init, writeTBLiteInfo
  use dftbp_geoopt_conjgrad, only : TConjGrad
  use dftbp_geoopt_filter, only : TFilter, TFilter_init
  use dftbp_geoopt_fire, only : TFire, TFire_init
  use dftbp_geoopt_gdiis, only : TDIIS
  use dftbp_geoopt_geoopt, only : TGeoOpt, geoOptTypes, reset, init
  use dftbp_geoopt_lbfgs, only : TLbfgs, TLbfgs_init
  use dftbp_geoopt_package, only : TOptimizer, createOptimizer, TOptTolerance
  use dftbp_geoopt_deprecated_steepdesc, only : TSteepDescDepr
  use dftbp_io_commonformats, only : format2Ue
  use dftbp_io_message, only : error, warning
  use dftbp_io_taggedoutput, only : TTaggedWriter, TTaggedWriter_init
  use dftbp_math_duplicate, only : isRepeated
  use dftbp_math_lapackroutines, only : matinv
  use dftbp_math_randomgenpool, only : TRandomGenPool, init
  use dftbp_math_ranlux, only : TRanlux, getRandom
  use dftbp_math_simplealgebra, only : determinant33
  use dftbp_md_andersentherm, only : TAndersenThermostat, init
  use dftbp_md_berendsentherm, only :TBerendsenThermostat, init
  use dftbp_md_dummytherm, only : TDummyThermostat, init
  use dftbp_md_mdcommon, only : TMDCommon, init
  use dftbp_md_mdintegrator, only : TMDIntegrator, init
  use dftbp_md_nhctherm, only : TNHCThermostat, init
  use dftbp_md_tempprofile, only : TTempProfile, TempProfile_init
  use dftbp_md_thermostat, only : TThermostat, init
  use dftbp_md_velocityverlet, only : TVelocityVerlet, init
  use dftbp_md_xlbomd, only : TXLBOMD, Xlbomd_init
  use dftbp_mixer_andersonmixer, only : TAndersonMixer, init
  use dftbp_mixer_broydenmixer, only : TBroydenMixer, init
  use dftbp_mixer_diismixer, only : TDIISMixer, init
  use dftbp_mixer_mixer, only : TMixer, mixerTypes, init
  use dftbp_mixer_simplemixer, only : TSimpleMixer, init
  use dftbp_reks_reks, only : TReksInp, TReksCalc, reksTypes, REKS_init
  use dftbp_solvation_cm5, only : TChargeModel5, TChargeModel5_init
  use dftbp_solvation_fieldscaling, only : TScaleExtEField, init_TScaleExtEField
  use dftbp_solvation_solvation, only : TSolvation
  use dftbp_solvation_solvinput, only : createSolvationModel, writeSolvationInfo
  use dftbp_timedep_linresp, only : LinResp_init
  use dftbp_timedep_linresptypes, only : TLinResp, linRespSolverTypes
  use dftbp_timedep_pprpa, only : TppRPAcal
  use dftbp_timedep_timeprop, only : TElecDynamics, TElecDynamics_init, tdSpinTypes
  use dftbp_type_commontypes, only : TOrbitals, TParallelKS, TParallelKS_init
  use dftbp_type_densedescr, only : TDenseDescr
  use dftbp_type_integral, only : TIntegral, TIntegral_init
  use dftbp_type_linkedlist, only : TListIntR1, TListCharLc, init, destruct, elemShape, intoArray,&
      & append
  use dftbp_type_multipole, only : TMultipole, TMultipole_init
  use dftbp_type_orbitals, only : getShellNames
  use dftbp_type_wrappedintr, only : TWrappedInt1
#:if WITH_MBD
  use dftbp_dftb_dispmbd, only :TDispMbd, TDispMbdInp, TDispMbd_init
#:endif
#:if WITH_OMP
  use omp_lib, only : omp_get_max_threads
#:endif
#:if WITH_SCALAPACK
  use dftbp_extlibs_scalapackfx, only : scalafx_getdescriptor, scalafx_getlocalshape
#:endif
#:if WITH_SOCKETS
  use dftbp_dftbplus_mainio, only : receiveGeometryFromSocket
  use dftbp_io_ipisocket, only : ipiSocketCommInp, ipiSocketComm
#:endif
#:if WITH_TRANSPORT
  use dftbp_dftbplus_inputdata, only : TNEGFInfo
  use dftbp_transport_negfint, only : TNegfInt, TNegfInt_init
  use dftbp_transport_negfvars, only : TTransPar
#:endif
  implicit none

  private
  public :: TDftbPlusMain, TCutoffs, TNegfInt
  public :: initReferenceCharges, updateReferenceShellCharges, initElectronNumber
#:if WITH_TRANSPORT
  public :: overrideContactCharges
#:endif
#:if WITH_SCALAPACK
  public :: getDenseDescBlacs
#:endif

  !> Interaction cutoff distances
  type :: TCutoffs

    !> Cutoff for overlap and Hamiltonian according to SK-files
    real(dp) :: skCutOff

    !> Cutoff for overlap and Hamiltonian according to SK-files minus possible cutoff reduction
    real(dp) :: camCutOff

    !> Cutoff for real-space g-summation in CAM Hartree-Fock contributions
    real(dp), allocatable :: gSummationCutoff

    !> Number of unit cells along each supercell folding direction to substract from MIC
    !! Wigner-Seitz cell construction
    integer, allocatable :: wignerSeitzReduction

    !> Cutoff for truncated long-range Gamma integral
    real(dp), allocatable :: gammaCutoff

    !> Max. of all cutoffs
    real(dp) :: mCutOff

  end type TCutoffs


#:if not WITH_TRANSPORT

  !> Dummy type for negf interface
  type :: TNegfInt
  end type TNegfInt

#:endif


  type :: TDftbPlusMain

    !> Is this calculation using a restarted input that does not require self consistency before
    !> moving to the post SCC loop part (i.e. Ehrenfest)
    logical :: tRestartNoSC = .false.

    !> Is the calculation SCC?
    logical :: tSccCalc

    !> SCC module internal variables
    type(TScc), allocatable :: scc

    !> Nr. of atoms
    integer :: nAtom

    !> Nr. of all (boundary condition images and original) atoms
    integer :: nAllAtom

    !> Nr. of original atom in central cell
    integer, allocatable :: img2CentCell(:)

    !> Nr. of different types (nAtom)
    integer :: nType

    !> Data type for atomic orbital information
    type(TOrbitals), allocatable :: orb

    !> Nr. of orbitals in the system
    integer :: nOrb

    !> Types of the atoms (nAllAtom)
    integer, allocatable :: species(:)

    !> Type of the atoms (nAtom)
    integer, allocatable :: species0(:)

    !> Coords of the atoms (3, nAllAtom)
    real(dp), allocatable :: coord(:,:)

    !> Coords in central cell (3, nAtom)
    real(dp), allocatable :: coord0(:,:)

    !> If calculation is periodic
    logical :: tPeriodic

    !> If the calculation is helical geometry
    logical :: tHelical

    !> Should central cell coordinates be output?
    logical :: tShowFoldedCoord

    !> How to calculate forces
    integer :: forceType

    !> Are atomic coordinates fractional?
    logical :: tFracCoord

    !> Tolerance for SCC cycle
    real(dp) :: sccTol

    !> Lattice vectors as columns
    real(dp), allocatable :: latVec(:,:)

    !> Origin of coordinate system for periodic systems
    real(dp), allocatable :: origin(:)

    !> Reciprocal lattice vectors as columns
    real(dp), allocatable :: recVec(:,:)

    !> Original lattice vectors used for optimising
    real(dp) :: origLatVec(3,3)

    !> Normalized vectors in those directions
    real(dp) :: normOrigLatVec(3,3)

    !> Reciprocal vectors in 2 pi units
    real(dp), allocatable :: invLatVec(:,:)

    !> Cell volume
    real(dp) :: cellVol

    !> Reciprocal cell volume
    real(dp) :: recCellVol

    !> Translation vecs for interacting image cells (3, nImgCell + 1)
    real(dp), allocatable :: cellVec(:,:)

    !> Cell vectors in absolute units
    real(dp), allocatable :: rCellVec(:,:)

    !> Index in cellVec for each atom
    integer, allocatable :: iCellVec(:)

    !> ADT for neighbour parameters
    type(TNeighbourList), allocatable :: neighbourList

    !> ADT for neighbour parameters, symmetric version for CAM calculations
    type(TSymNeighbourList), allocatable :: symNeighbourList

    !> Nr. of neighbours for atoms out to max interaction distance (excluding Ewald terms)
    integer, allocatable :: nNeighbourSK(:)

    !> Number of neighbours for each of the atoms for the exchange contributions of CAM functionals
    integer, allocatable :: nNeighbourCam(:)

    !> Symmetric neighbour list version of nNeighbourCam
    integer, allocatable :: nNeighbourCamSym(:)

    !> H/S sparse matrices indexing array for atomic blocks
    integer, allocatable :: iSparseStart(:,:)

    !> Self energy (orbital, atom)
    real(dp), allocatable :: atomEigVal(:,:)

    !> Reference n_0 charges for each atom
    real(dp), allocatable :: referenceN0(:,:)

    !> List of atomic masses
    real(dp), allocatable :: mass(:)

    !> List of atomic masses for each species
    real(dp), allocatable :: speciesMass(:)

    !> Hamiltonian type
    integer :: hamiltonianType

    !> Raw H^0 hamiltonian data
    type(TSlakoCont) :: skHamCont

    !> Raw overlap hamiltonian data
    type(TSlakoCont) :: skOverCont

    !> Repulsive (force-field like) interactions
    class(TRepulsive), allocatable :: repulsive

    !> Cut off distances for various types of interaction
    type(TCutoffs) :: cutOff

    !> Charge per atomic shell (shell, atom, spin channel)
    real(dp), allocatable :: chargePerShell(:,:,:)

    !> Charge par atom (atom, spin channel)
    real(dp), allocatable :: chargePerAtom(:,:)

    !> Integral container
    type(TIntegral) :: ints

    !> Nr. of K-points
    integer :: nKPoint

    !> The k-points
    real(dp), allocatable :: kPoint(:,:)

    !> Weight of the K-Points
    real(dp), allocatable :: kWeight(:)

    !> Coefficients of the lattice vectors in the linear combination for the super lattice vectors
    !! (should be integer values) and shift of the grid along the three small reciprocal lattice
    !! vectors (between 0.0 and 1.0)
    real(dp), allocatable :: supercellFoldingMatrix(:,:)

    !> Three diagonal elements of supercell folding coefficient matrix
    integer, allocatable :: supercellFoldingDiag(:)

    !> External pressure if periodic
    real(dp) :: extPressure

    !> Barostat used if MD and periodic
    logical :: tBarostat

    !> Barostat coupling strength
    real(dp) :: BarostatStrength

    !> H and S are real
    logical :: tRealHS

    !> Nr. of electrons
    real(dp), allocatable :: nEl(:)

    !> Nr. of all electrons if neutral
    real(dp) :: nEl0

    !> Spin W values
    real(dp), allocatable :: spinW(:,:,:)

    !> Spin orbit constants
    real(dp), allocatable :: xi(:,:)

    !> DFTB+U calculation, if relevant
    type(TDftbU), allocatable :: dftbU

    !> Electron temperature
    real(dp) :: tempElec

    !> If K points should filled separately
    logical :: tFillKSep

    !> Fix Fermi energy at specified value
    logical :: tFixEf

    !> Fermi energy for each spin
    real(dp), allocatable :: Ef(:)

    !> Filling temp, as updated by MD.
    logical :: tSetFillingTemp

    !> Choice of electron distribution function, defaults to Fermi
    integer :: iDistribFn = fillingTypes%Fermi

    !> Atomic kinetic temperature
    real(dp) :: tempAtom

    !> MD stepsize
    real(dp) :: deltaT

    !> Maximal number of SCC iterations
    integer :: maxSccIter

    !> Minimal number of SCC iterations
    integer :: minSccIter

    !> Is this a spin polarized calculation?
    logical :: tSpin

    !> Number of spin components, 1 is unpolarised, 2 is polarised, 4 is noncolinear / spin-orbit
    integer :: nSpin

    !> Nr. of independent spin cases to diagonalise/solve independently (nSpin {1,4}:1, 2:2)
    integer :: nIndepSpin

    !> Is there spin-orbit coupling
    logical :: tSpinOrbit

    !> Use block like dual representation for spin orbit
    logical :: tDualSpinOrbit

    !> Number of atomic dipole moment components
    integer :: nDipole = 0

    !> Number of atomic quadrupole moment components
    integer :: nQuadrupole = 0

    !> Is there a complex hamiltonian contribution in real space
    logical :: tImHam

    !> Is this a two component calculation (spin orbit or non-collinear spin)
    logical :: t2Component

    !> Common Fermi level across spin channels
    logical :: tSpinSharedEf

    !> Geometry optimisation needed?
    logical :: isGeoOpt

    !> Optimise coordinates inside unit cell (periodic)?
    logical :: tCoordOpt

    !> Optimise lattice constants?
    logical :: tLatOpt

    !> Fix angles between lattice vectors when optimising?
    logical :: tLatOptFixAng

    !> Fix length of specified lattice vectors when optimising?
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
    type(TElStatPotentials), allocatable :: electrostatPot

    !> Calculate localised orbitals?
    logical :: tLocalise

    !> Do we need to show Mulliken charges?
    logical :: tPrintMulliken

    !> Logical to determine whether to calculate net charge per atom (qNetAtom)
    logical :: isQNetAllocated

    !> Do need net atomic charges?
    logical :: tNetAtomCharges

    !> Should the net charges be printed
    logical :: tPrintNetAtomCharges

    !> Calculate an electric dipole?
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

    !> Optimization of conical intersections
    logical :: isCIopt = .false.

    !> Are forces being returned
    logical :: tPrintForces

    !> Number of moved atoms
    integer :: nMovedAtom

    !> Index of the moved atoms
    integer, allocatable :: indMovedAtom(:)

    !> Index of the atoms for which second derivatives are computed
    integer, allocatable :: indDerivAtom(:)

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

    !> Density functional tight binding perturbation theory
    logical :: doPerturbation = .false.

    !> Self-consistency tolerance for perturbation (if SCC)
    real(dp) :: perturbSccTol = 0.0_dp

    !> Maximum iteration for self-consistency in the perturbation routines
    integer :: maxPerturbIter = 0

    !> Require converged perturbation (if true, terminate on failure, otherwise return NaN for
    !> non-converged)
    logical :: isPerturbConvRequired = .true.

    !> Density functional tight binding perturbation for each geometry step
    logical :: doPerturbEachGeom = .false.

    !> Response property calculations
    type(TResponse), allocatable :: response

    !> Static polarisability
    logical :: isEResp = .false.

    !> Dynamic polarisability at finite frequencies
    real(dp), allocatable :: dynRespEFreq(:)

    !> Is the response kernel (and frontier eigenvalue derivatives) calculated by perturbation
    logical :: isKernelResp

    !> Should the response Kernel use RPA (non-SCC) or self-consistent
    logical :: isRespKernelRPA

    !> Dynamic polarisability at finite frequencies
    real(dp), allocatable :: dynKernelFreq(:)

    !> Electric static polarisability
    real(dp), allocatable :: polarisability(:,:,:)

    !> Number of electrons  at the Fermi energy
    real(dp), allocatable :: neFermi(:)

    !> Derivatives of the Fermi energy [spin, perturbation]
    real(dp), allocatable :: dEfdE(:,:)

    !> Derivatives of the DFTB eigenvalues with respect to perturbation
    real(dp), allocatable :: dEidE(:,:,:,:)

    !> Derivatives of Mulliken charges with respect to perturbation
    real(dp), allocatable :: dqOut(:,:,:,:)

    !> Use commands from socket communication to control the run
    logical :: tSocket

    !> Socket details
  #:if WITH_SOCKETS
    type(ipiSocketComm), allocatable :: socket
  #:endif

    !> File containing output geometry
    character(lc) :: geoOutFile

    !> Append geometries in the output?
    logical :: tAppendGeo

    !> Use converged SCC charges for properties like forces or charge dependent dispersion
    logical :: isSccConvRequired

    !> Labels of atomic species
    character(mc), allocatable :: speciesName(:)

    !> General geometry optimiser
    type(TGeoOpt), allocatable :: pGeoCoordOpt

    !> Geometry optimiser for lattice consts
    type(TGeoOpt), allocatable :: pGeoLatOpt

    !> Coordinate transformation filter
    type(TFilter), allocatable :: filter

    !> General geometry optimiser
    class(TOptimizer), allocatable :: geoOpt

    !> Convergence thresholds for geometry optimiser
    type(TOptTolerance) :: optTol

    real(dp) :: elast
    real(dp), allocatable :: gcurr(:), glast(:), displ(:)

    !> Charge mixer
    type(TMixer), allocatable :: pChrgMixer

    !> MD Framework
    type(TMDCommon), allocatable :: pMDFrame

    !> MD integrator
    type(TMDIntegrator), allocatable :: pMDIntegrator

    !> Temperature profile driver in MD
    type(TTempProfile), allocatable :: temperatureProfile

    !> Geometry optimiser
    type(TNumDerivs), allocatable :: derivDriver

    !> Total charge
    real(dp) :: nrChrg

    !> Spin polarisation
    real(dp) :: nrSpinPol

    !> Is the check-sum for charges read externally to be used?
    logical :: tSkipChrgChecksum

    !> Reference neutral atomic occupations
    real(dp), allocatable :: q0(:, :, :)

    !> Shell resolved neutral reference
    real(dp), allocatable :: qShell0(:,:)

    !> Input charges (for potentials)
    real(dp), allocatable :: qInput(:, :, :)

    !> Output charges
    real(dp), allocatable :: qOutput(:, :, :)

    !> Charge differences between input and output charges
    real(dp), allocatable :: qDiff(:, :, :)

    !> Net (on-site only contributions) charge per atom
    real(dp), allocatable :: qNetAtom(:)

    !> Input Mulliken block charges (diagonal part == Mulliken charges)
    real(dp), allocatable :: qBlockIn(:, :, :, :)

    !> Output Mulliken block charges
    real(dp), allocatable :: qBlockOut(:, :, :, :)

    !> Imaginary part of input Mulliken block charges
    real(dp), allocatable :: qiBlockIn(:, :, :, :)

    !> Imaginary part of output Mulliken block charges
    real(dp), allocatable :: qiBlockOut(:, :, :, :)

    !> Input charges packed into unique equivalence elements
    real(dp), allocatable :: qInpRed(:)

    !> Output charges packed into unique equivalence elements
    real(dp), allocatable :: qOutRed(:)

    !> Charge differences packed into unique equivalence elements
    real(dp), allocatable :: qDiffRed(:)

    !> Orbital equivalence relations
    integer, allocatable :: iEqOrbitals(:,:,:)

    !> Dipolar equivalence relations
    integer, allocatable :: iEqDipole(:,:)

    !> Quadrupolar equivalence relations
    integer, allocatable :: iEqQuadrupole(:,:)

    !> Nr. of inequivalent orbitals
    integer :: nIneqOrb

    !> Nr. of inequivalent dipoles
    integer :: nIneqDip

    !> Nr. of inequivalent quadrupoles
    integer :: nIneqQuad

    !> Multipole moments for the input charges
    type(TMultipole) :: multipoleInp

    !> Multipole moments for the output charges
    type(TMultipole) :: multipoleOut

    !> Nr. of elements to go through the mixer - may contain reduced orbitals and also orbital
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

    !> Unique Hubbard U (needed for being able to calculate equivalency relations)
    type(TUniqueHubbard), allocatable :: uniqHubbU

    ! External charges

    !> If external charges must be considered
    logical :: tExtChrg

    !> Nr. of external charges
    integer :: nExtChrg

    !> Electric field
    type(TEfield), allocatable :: eField

    !> Is an arbitrary external field (including electric) present
    logical :: isExtField

    !> Partial density of states (PDOS) projection regions
    type(TListIntR1) :: iOrbRegion

    !> PDOS region labels
    type(TListCharLc) :: regionLabels

    !> Third order DFTB
    logical :: t3rd = .false.

    !> Full 3rd order or only atomic site
    logical :: t3rdFull = .false.

    !> Data structure for 3rd order
    type(TThirdOrder), allocatable :: thirdOrd

    !> Correction to energy from on-site matrix elements
    real(dp), allocatable :: onSiteElements(:,:,:,:)

    !> Correction to dipole momements on-site matrix elements
    real(dp), allocatable :: onSiteDipole(:,:)

    !> Should block charges be mixed as well as charges
    logical :: tMixBlockCharges

    !> Calculate Casida linear response excitations
    logical :: isLinResp

    !> Calculate Z vector for excited properties
    logical :: tLinRespZVect = .false.

    !> Data type for pp-RPA
    type(TppRPAcal), allocatable :: ppRPA

    !> Print eigenvectors
    logical :: tPrintExcitedEigVecs

    !> Data type for linear response
    type(TLinResp), allocatable :: linearResponse

    !> Whether to run a hybrid xc-functional calculation
    logical :: isHybridXc

    !> Whether constraints are imposed on electronic ground state
    logical :: isElecConstr

    !> Hybrid xc-functional data
    class(THybridXcFunc), allocatable :: hybridXc

    !> Holds real and complex delta density matrices and pointers
    type(TDensityMatrix) :: densityMatrix

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

    !> Produce charges.bin
    logical :: tWriteCharges

    !> Should charges written to disc be in ascii or binary format?
    logical :: tWriteChrgAscii

    !> Produce tagged output?
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

    !> Dispersion data and calculations
    class(TDispersionIface), allocatable :: dispersion

    !> Solvation data and calculations
    class(TSolvation), allocatable :: solvation

    !> Dielectric scaling of electric fields (if relevant)
    type(TScaleExtEField) :: eFieldScaling

    !> Charge Model 5 for printout
    type(TChargeModel5), allocatable :: cm5Cont

    !> Write cavity information as cosmo file
    logical :: tWriteCosmoFile

    !> Structure holding electronic constraints
    type(TElecConstraint), allocatable :: elecConstraint

    !> Library interface handler
    type(TTBLite), allocatable :: tblite

    !> Can stress be calculated?
    logical :: tStress

    !> Should XLBOMD be used in MD
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

    !> Electronic filling
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

    !> Density matrix
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

    !> Excited state force addition (xyz,atom,state)
    real(dp), allocatable :: excitedDerivs(:,:,:)

    !> Nonadiabatic coupling vectors
    real(dp), allocatable :: naCouplings(:,:,:)

    !> Dipole moments, when available, for whichever determinants are present
    real(dp), allocatable :: dipoleMoment(:, :)

    !> Additional dipole moment related message to write out
    character(lc) :: dipoleMessage

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

    !> File descriptor for the human readable output
    type(TFileDescr) :: fdDetailedOut

    !> File descriptor for extra MD output
    type(TFileDescr) :: fdMd

    !> Contains (iK, iS) tuples to be processed in parallel by various processor groups
    type(TParallelKS) :: parallelKS

    !> True, if electron dynamics input block is present
    logical :: isElecDyn

    !> Electron dynamics
    type(TElecDynamics), allocatable :: electronDynamics

    !> Electronic structure solver
    type(TElectronicSolver) :: electronicSolver

    !> Are large dense matrices required?
    logical :: tLargeDenseMatrices

    !> Derivative of cell volume wrt to lattice vectors, needed for pV term
    real(dp) :: extLatDerivs(3,3)

    !> Internal pressure within the cell
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
    real(dp) :: totalStress(3,3) = 0.0_dp

    !> Stress tensors for determinants if using TI-DFTB
    real(dp), allocatable :: mixedStress(:,:), tripletStress(:,:)

    ! Tagged writer
    type(TTaggedWriter) :: taggedWriter

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

    !> Whether a contact Hamiltonian is being computed and stored
    logical :: isAContactCalc

    !> Whether Poisson solver is invoked
    logical :: tPoisson

    !> Whether the scc (2nd-order) potentials should be updated after the diagonalisation
    logical :: updateSccAfterDiag

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

    !> Shell-resolved electrostatic potential shifts uploaded from contacts
    real(dp), allocatable :: shiftPerLUp(:,:)

    !> Orbital-resolved charges uploaded from contacts
    real(dp), allocatable :: chargeUp(:,:,:)

    !> Block populations uploaded from contacts
    real(dp), allocatable :: blockUp(:,:,:,:)

    !> Details of energy interval for tunneling used in output
    real(dp) :: Emin, Emax, Estep

    !> Electrostatics type (either gammafunctional or poisson)
    integer :: electrostatics

    !> List of atoms in the central cell (or device region if transport)
    integer, allocatable :: iAtInCentralRegion(:)

    !> Correction for {O,N}-X bonds
    type(THalogenX), allocatable :: halogenXCorrection

    !> All of the excited energies actually solved by Casida routines (if used)
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
    real(dp), allocatable :: deltaRhoDets(:,:,:,:)

    !> Data type for REKS
    type(TReksCalc), allocatable :: reks

    !> Atomic charge contribution in excited state
    real(dp), allocatable :: dQAtomEx(:)

    !> Boundary condition
    type(TBoundaryConditions) :: boundaryCond

    !> Whether the order of the atoms matter. Typically the case, when properties were specified
    !> based on atom numbers (e.g. custom occupations). In that case setting a different order
    !> of the atoms via the API is forbidden.
    logical :: atomOrderMatters = .false.

  #:if WITH_SCALAPACK

    !> Should the dense matrices be re-ordered for sparse operations
    logical :: isSparseReorderRequired = .false.

    !> Re-orded real eigenvectors (if required)
    real(dp), allocatable :: eigVecsRealReordered(:,:,:)

    !> Re-orded complex eigenvectors (if required)
    complex(dp), allocatable :: eigVecsCplxReordered(:,:,:)

  #:endif

  contains

    procedure :: initProgramVariables
    procedure :: setEquivalencyRelations
    procedure :: initializeCharges
    procedure :: destructProgramVariables
  #:if WITH_SOCKETS
    procedure :: initSocket
  #:endif
    procedure :: initOutputFiles
    procedure :: initArrays
    procedure :: initDetArrays
    procedure :: allocateDenseMatrices
    procedure :: getDenseDescCommon
    procedure :: ensureConstrainedDftbReqs
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

    !> Mixer number
    integer :: iMixer

    !> Simple mixer (if used)
    type(TSimpleMixer), allocatable :: pSimpleMixer

    !> Anderson mixer (if used)
    type(TAndersonMixer), allocatable :: pAndersonMixer

    !> Broyden mixer (if used)
    type(TBroydenMixer), allocatable :: pBroydenMixer

    !> DIIS mixer (if used)
    type(TDIISMixer), allocatable :: pDIISMixer

    ! Geometry optimiser related local variables

    !> Conjugate gradient driver
    type(TConjGrad), allocatable :: pConjGrad

    !> Steepest descent driver
    type(TSteepDescDepr), allocatable :: pSteepDesc

    !> Conjugate gradient driver
    type(TConjGrad), allocatable :: pConjGradLat

    !> Steepest descent driver
    type(TSteepDescDepr), allocatable :: pSteepDescLat

    !> Gradient DIIS driver
    type(TDIIS), allocatable :: pDIIS

    !> LBFGS driver for geometry optimisation
    type(TLbfgs), allocatable :: pLbfgs

    !> LBFGS driver for lattice optimisation
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
    type(TSimpleDftD3), allocatable :: sdftd3
    type(TDispDftD4), allocatable :: dftd4
  #:if WITH_MBD
    type(TDispMbd), allocatable :: mbd
  #:endif

    logical :: tHHRepulsion

    character(lc) :: tmpStr
    integer, allocatable :: tmpir1(:)

    character(lc) :: strTmp, strTmp2

    !> Flag to check for first cycle through a loop
    logical :: tFirst

    real(dp) :: rTmp

    !> Flag if some files do exist or not
    logical :: tExist

    ! Damped interactions
    type(TShortGammaDamp) :: shortGammaDamp
    type(TThirdOrderInp) :: thirdInp

    ! PDOS stuff
    integer :: iReg, nAtomRegion, nOrbRegion, iTmp
    integer, allocatable :: iAtomRegion(:)
    integer :: valshape(1)

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

    real(dp), allocatable :: hubbU(:,:)
    type(TShortGammaInput), allocatable :: shortGammaInput
    type(TCoulombInput), allocatable :: coulombInput
    type(TPoissonInput), allocatable :: poissonInput

    logical :: tGeoOptRequiresEgy, isOnsiteCorrected
    type(TStatus) :: errStatus
    integer :: nLocalRows, nLocalCols

    @:ASSERT(input%tInitialized)
    write(stdOut, "(/, A)") "Starting initialization..."
    write(stdOut, "(A80)") repeat("-", 80)

    call env%initGlobalTimer(input%ctrl%timingLevel, "DFTB+ running times", stdOut)
    call env%globalTimer%startTimer(globalTimers%globalInit)

    ! Set the same access for readwrite as for write (we do not open any files in readwrite mode)
    call setDefaultBinaryAccess(input%ctrl%binaryAccessTypes(1), input%ctrl%binaryAccessTypes(2),&
        & input%ctrl%binaryAccessTypes(2))

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
    this%nIndepSpin = this%nSpin

    this%tSpinSharedEf = input%ctrl%tSpinSharedEf
    this%tSpinOrbit = input%ctrl%tSpinOrbit
    this%tDualSpinOrbit = input%ctrl%tDualSpinOrbit
    this%t2Component = input%ctrl%t2Component
    this%isXlbomd = allocated(input%ctrl%xlbomd)
    this%isElecConstr = allocated(input%ctrl%elecConstraintInp)
    this%isElecDyn = allocated(input%ctrl%elecDynInp)
    this%isHybridXc = allocated(input%ctrl%hybridXcInp)
    if (this%isHybridXc) then
      allocate(this%symNeighbourList)
    end if

    if (this%t2Component) then
      this%nSpin = 4
      this%nIndepSpin = 1
    end if

    if (this%nSpin /= 2 .and. this%tSpinSharedEf) then
      call error("Colinear spin polarization required for shared Ef over spin channels")
    end if

    call initGeometry_(input, this%nAtom, this%nType, this%tPeriodic, this%tHelical,&
        & this%boundaryCond, this%coord0, this%species0, this%tCoordsChanged, this%tLatticeChanged,&
        & this%latVec, this%origin, this%recVec, this%invLatVec, this%cellVol, this%recCellVol,&
        & errStatus)
    if (errStatus%hasError()) then
      call error(errStatus%message)
    end if

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

    select case(this%hamiltonianType)
    case default
      call error("Invalid Hamiltonian")
    case(hamiltonianTypes%dftb)
      this%orb = input%slako%orb
    case(hamiltonianTypes%xtb)
      if (.not.allocated(input%ctrl%tbliteInp)) then
        call error("xTB calculation supported only with tblite library")
      end if
      allocate(this%tblite)
      if (this%tPeriodic) then
        call TTBlite_init(this%tblite, input%ctrl%tbliteInp, this%nAtom, this%species0, &
            & this%speciesName, this%coord0, this%latVec)
      else
        call TTBlite_init(this%tblite, input%ctrl%tbliteInp, this%nAtom, this%species0, &
            & this%speciesName, this%coord0)
      end if
      allocate(input%slako%orb)
      call this%tblite%getOrbitalInfo(this%species0, input%slako%orb)
      allocate(input%slako%skOcc(input%slako%orb%mShell, input%geom%nSpecies))
      call this%tblite%getReferenceN0(this%species0, input%slako%skOcc)
      ! Workaround: ifort 2021.7
      ! Assignment of derived type instances with allocatable components seems to be broken,
      ! resulting in a strange run-time error message. Turning it into move_alloc seems to avoid it.
      !this%orb = input%slako%orb
      call move_alloc(input%slako%orb, this%orb)
    end select
    this%nOrb = this%orb%nOrb

    ! start by assuming stress can be calculated if periodic
    this%tStress = this%tPeriodic

    ! Brillouin zone sampling
    if (this%tPeriodic .or. this%tHelical) then
      if (input%ctrl%poorKSampling) then
        call warning("The suplied k-points are probably not accurate for properties requiring&
            & integration over the Brillouin zone")
      end if
      this%nKPoint = input%ctrl%nKPoint
      allocate(this%kPoint(size(input%ctrl%KPoint,dim=1), this%nKPoint))
      allocate(this%kWeight(this%nKPoint))
      this%kPoint(:,:) = input%ctrl%KPoint
      if (sum(input%ctrl%kWeight) < epsilon(1.0_dp)) then
        call error("Sum of k-point weights should be greater than zero!")
      end if
      this%kWeight(:) = input%ctrl%kWeight / sum(input%ctrl%kWeight)
      if (this%tHelical) then
        if (any(abs(this%kPoint(2,:) * nint(this%latVec(3,1)) - nint(this%kPoint(2,:) *&
            & nint(this%latVec(3,1)))) > input%ctrl%helicalSymTol)) then
          call warning("Specified k-value(s) incommensurate with C_n symmetry operation.")
        end if
      end if
    else
      this%nKPoint = 1
      allocate(this%kPoint(3, this%nKPoint))
      allocate(this%kWeight(this%nKPoint))
      this%kPoint(:,1) = 0.0_dp
      this%kWeight(1) = 1.0_dp
    end if

    this%tRealHS = .true.
    if (this%tPeriodic .or. this%tHelical) then
      if (size(this%kPoint, dim=2) == 1 .and. all(this%kPoint(:, 1) == 0.0_dp)) then
        this%tRealHS = .true.
      else
        this%tRealHS = .false.
      end if
    end if

  #:if WITH_MPI
    if (input%ctrl%parallelOpts%nGroup > this%nIndepSpin * this%nKPoint&
        & .and. (.not. (this%isHybridXc .and. (.not. this%tRealHS)))) then
      write(stdOut, *) "Parallel groups only relevant for tasks split over sufficient spins and/or&
          & k-points"
      write(tmpStr,"('Nr. groups:',I4,', Nr. indepdendent spins times k-points:',I4)")&
          & input%ctrl%parallelOpts%nGroup, this%nIndepSpin * this%nKPoint
      call error(trim(tmpStr))
    end if

    call env%initMpi(input%ctrl%parallelOpts%nGroup)

    if (this%isHybridXc .and. (.not. this%tRealHS)&
        & .and. (input%ctrl%parallelOpts%nGroup /= env%mpi%globalComm%size)) then
      ! General k-point case (parallelized over k-points)
      write(tmpStr, "(A,I0,A,I0,A)") 'For range-separated calculations beyond the Gamma point, the&
          & number of MPI' // NEW_LINE('A') // '   groups must match the total number of MPI&
          & processes.' // NEW_LINE('A') // '   Obtained (', input%ctrl%parallelOpts%nGroup, ')&
          & groups for (', env%mpi%globalComm%size, ') total MPI processes!'
      call error(trim(tmpStr))
    end if
  #:endif

  #:if WITH_SCALAPACK
    call initBlacs(input%ctrl%parallelOpts%blacsOpts, this%nAtom, this%nOrb,&
        & this%t2Component, env, errStatus)
    if (errStatus%hasError()) then
      if (errStatus%code == -1) then
        call warning("Insufficient atoms for this number of MPI processors")
      end if
      call error(errStatus%message)
    end if
  #:endif
    call TParallelKS_init(this%parallelKS, env, this%nKPoint, this%nIndepSpin)

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
      this%maxSccIter = input%ctrl%maxSccIter
    else
      this%maxSccIter = 1
    end if
    if (this%maxSccIter < 1) then
      call error("SCC iterations must be larger than 0")
    end if
    if (this%tSccCalc) then
      if (this%isElecDyn) then
        if (input%ctrl%elecDynInp%tReadRestart .and. .not.input%ctrl%elecDynInp%tPopulations) then
          this%maxSccIter = 0
          this%isSccConvRequired = .false.
          this%tRestartNoSC = .true.
        end if
      end if
    end if

    this%tWriteHS = input%ctrl%tWriteHS
    this%tWriteRealHS = input%ctrl%tWriteRealHS

    select case(this%hamiltonianType)
    case default
      call error("Invalid Hamiltonian")
    case(hamiltonianTypes%dftb)
      ! Slater-Koster tables
      this%skHamCont = input%slako%skHamCont
      this%skOverCont = input%slako%skOverCont
      call initRepulsive_(this%nAtom, this%tPeriodic, this%tHelical, input%slako%pairRepulsives,&
          & input%ctrl%chimesRepInput, this%speciesName, this%species0, this%repulsive)
      allocate(this%atomEigVal(this%orb%mShell, this%nType))
      @:ASSERT(size(input%slako%skSelf, dim=1) == this%orb%mShell)
      @:ASSERT(size(input%slako%skSelf, dim=2) == size(this%atomEigVal, dim=2))
      this%atomEigVal(:,:) = input%slako%skSelf(1:this%orb%mShell, :)

      @:ASSERT(size(input%slako%skOcc, dim=1) >= this%orb%mShell)
      @:ASSERT(size(input%slako%mass) == this%nType)
      allocate(this%speciesMass(this%nType))
      this%speciesMass(:) = input%slako%mass(:)
    case(hamiltonianTypes%xtb)
      allocate(this%speciesMass(this%nType))
      this%speciesMass(:) = getAtomicMass(this%speciesName)
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
      this%cutOff%mCutoff = this%cutOff%skCutOff
      if (allocated(this%repulsive)) then
        this%cutOff%mCutOff = max(this%cutOff%mCutOff, this%repulsive%getRCutOff())
      end if
    case(hamiltonianTypes%xtb)
      this%cutOff%skCutoff = this%tblite%getRCutoff()
      this%cutOff%mCutoff = this%cutOff%skCutoff
    end select

    call initHubbardUs_(input, this%orb, this%hamiltonianType, hubbU)
    if (this%tSccCalc) then
      allocate(this%uniqHubbU)
      call TUniqueHubbard_init(this%uniqHubbU, hubbU, this%orb)
    end if

    call initReferencePopulation_(input, this%orb, this%hamiltonianType, this%referenceN0)

    this%atomOrderMatters = this%atomOrderMatters .or. allocated(input%ctrl%customOccAtoms)
    call initReferenceCharges(this%species0, this%orb, this%referenceN0, this%nSpin, this%q0,&
        & this%qShell0, input%ctrl%customOccAtoms, input%ctrl%customOccFillings)

    this%nrChrg = input%ctrl%nrChrg
    this%nrSpinPol = input%ctrl%nrSpinPol
    call initElectronNumber(this%q0, this%nrChrg, this%nrSpinPol, this%nSpin, this%orb,&
        & this%nEl0, this%nEl)
    call initElectronFilling_(input, this%nSpin, this%Ef, this%iDistribFn, this%tempElec,&
        & this%tFixEf, this%tSetFillingTemp, this%tFillKSep)

    call ensureSolverCompatibility(input%ctrl%solver%iSolver, this%kPoint, input%ctrl%parallelOpts,&
        & this%nIndepSpin, this%tempElec)
    call getBufferedCholesky_(this%tRealHS, this%parallelKS%nLocalKS, nBufferedCholesky)
    call TElectronicSolver_init(this%electronicSolver, input%ctrl%solver%iSolver, nBufferedCholesky)

  #:if WITH_TRANSPORT
    this%tTunn = input%ginfo%tundos%defined
    this%tLocalCurrents = input%ginfo%greendens%doLocalCurr
    this%tNegf = (this%electronicSolver%iSolver == electronicSolverTypes%GF) .or. this%tTunn&
        & .or. this%tLocalCurrents
    this%tUpload = input%transpar%taskUpload&
        & .and. .not. (this%electronicSolver%iSolver == electronicSolverTypes%OnlyTransport&
        & .and. .not. this%tSccCalc)
    if (this%tUpload) then
      call initUploadArrays_(input%transpar, this%orb, this%nSpin, this%tMixBlockCharges,&
          & this%shiftPerLUp, this%chargeUp, this%blockUp)
    end if
    call initTransport_(this, env, input, this%electronicSolver, this%nSpin, this%tempElec,&
        & this%tNegf, this%isAContactCalc, this%mu, this%negfInt, this%ginfo, this%transpar,&
        & this%writeTunn, this%tWriteLDOS, this%regionLabelLDOS, errStatus)
    if (errStatus%hasError()) then
      call error(errStatus%message)
    end if
  #:else
    this%tTunn = .false.
    this%tLocalCurrents = .false.
    this%tNegf = .false.
    this%tUpload = .false.
    this%isAContactCalc = .false.
  #:endif

    this%tPoisson = input%ctrl%tPoisson .and. this%tSccCalc
    this%updateSccAfterDiag = input%ctrl%updateSccAfterDiag

    this%tExtChrg = .false.
    if (this%tSccCalc .and. .not.allocated(this%tblite)) then
      call initShortGammaDamping_(input%ctrl, this%speciesMass, shortGammaDamp)
      if (this%tPoisson) then
        #:block REQUIRES_COMPONENT('Poisson-solver', WITH_POISSON)
          call initPoissonInput_(input, this%nAtom, this%nType, this%species0, this%coord0,&
              & this%tPeriodic, this%latVec, this%orb, hubbU, poissonInput, this%shiftPerLUp)
        #:endblock
      else
        call initShortGammaInput_(input%ctrl, this%speciesMass, this%uniqHubbU, shortGammaDamp,&
            & shortGammaInput)
        call initCoulombInput_(env, input%ctrl%ewaldAlpha, input%ctrl%tolEwald,&
            & this%boundaryCond%iBoundaryCondition, coulombInput)
      end if
      call initSccCalculator_(env, this%orb, input%ctrl, this%boundaryCond%iBoundaryCondition,&
          & coulombInput, shortGammaInput, poissonInput, this%scc)

      ! Stress calculation does not work if external charges are involved
      this%nExtChrg = input%ctrl%nExtChrg
      this%tExtChrg = this%nExtChrg > 0
      this%tStress = this%tStress .and. .not. this%tExtChrg

      if (this%tExtChrg .and. this%hamiltonianType == hamiltonianTypes%xtb) then
        call error("External charges not currently supported for xTB hamiltonians")
      end if

      ! Longest cut-off including the softening part of gamma
      this%cutOff%mCutOff = max(this%cutOff%mCutOff, this%scc%getCutOff())

      if (input%ctrl%t3rd .and. input%ctrl%tShellResolved) then
        call error("Onsite third order DFTB only compatible with shell non-resolved SCC")
      end if

      ! Initialize full 3rd order module
      this%t3rd = input%ctrl%t3rd
      this%t3rdFull = input%ctrl%t3rdFull
      if (this%t3rdFull) then
        @:ASSERT(this%tSccCalc)
        thirdInp%orb => this%orb
        thirdInp%hubbUs = hubbU
        thirdInp%hubbUDerivs = input%ctrl%hubDerivs
        allocate(thirdInp%damped(this%nType))
        thirdInp%damped(:) = shortGammaDamp%isDamped
        thirdInp%dampExp = shortGammaDamp%exponent
        thirdInp%shellResolved = input%ctrl%tShellResolved
        allocate(this%thirdOrd)
        call ThirdOrder_init(this%thirdOrd, thirdInp)
        this%cutOff%mCutOff = max(this%cutOff%mCutOff, this%thirdOrd%getCutOff())
      end if
    end if

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

    if (allocated(this%symNeighbourList)) then
      allocate(this%symNeighbourList%coord(3, this%nAllAtom))
      allocate(this%symNeighbourList%species(this%nAllAtom))
      allocate(this%symNeighbourList%img2CentCell(this%nAllAtom))
      allocate(this%symNeighbourList%iCellVec(this%nAllAtom))
      allocate(this%symNeighbourList%iPair(0, this%nAtom))
    end if

    ! Check if multipolar contributions are required
    if (allocated(this%tblite)) then
      call this%tblite%getMultipoleInfo(this%nDipole, this%nQuadrupole)
    end if
    call TMultipole_init(this%multipoleOut, this%nAtom, this%nDipole, this%nQuadrupole,&
        & this%nSpin)
    this%multipoleInp = this%multipoleOut

    ! Initialize Hamilton and overlap
    this%tImHam = this%tDualSpinOrbit .or. (this%tSpinOrbit .and. allocated(this%dftbU))
    if (this%tSccCalc .and. .not.allocated(this%reks)) then
      allocate(this%chargePerShell(this%orb%mShell,this%nAtom,this%nSpin))
    else
      allocate(this%chargePerShell(0,0,0))
    end if
    call TIntegral_init(this%ints, this%nSpin, .not. allocated(this%reks), this%tImHam,&
        & this%nDipole, this%nQuadrupole)
    allocate(this%iSparseStart(0, this%nAtom))

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
    this%tWriteCosmoFile = input%ctrl%tWriteCosmoFile

    if (allocated(input%ctrl%electricField) .or. allocated(input%ctrl%atomicExtPotential)) then
      allocate(this%eField)
    end if
    this%isExtField = allocated(this%eField)

    if (allocated(input%ctrl%electricField)) then
      allocate(this%eField%EFieldStrength)
      this%eField%EFieldStrength = input%ctrl%electricField%EFieldStrength
      this%eField%EfieldVector(:) = input%ctrl%electricField%EfieldVector(:)
      this%eField%isTDEfield = input%ctrl%electricField%isTDEfield
      this%eField%EfieldOmega = input%ctrl%electricField%EfieldOmega
      this%eField%EfieldPhase = input%ctrl%electricField%EfieldPhase
      if (this%eField%isTDEfield .and. .not. this%tMD) then
        call error ("Time dependent electric fields only possible for MD!")
      end if
      ! parser should catch all of these:
      @:ASSERT(.not.this%eField%isTDEfield .or. this%tMD)
    end if

    this%tMulliken = input%ctrl%tMulliken .or. this%tPrintMulliken .or. this%isExtField .or.&
        & this%tFixEf .or. this%tSpinSharedEf .or. this%isHybridXc .or.&
        & this%electronicSolver%iSolver == electronicSolverTypes%GF
    this%tAtomicEnergy = input%ctrl%tAtomicEnergy
    this%tPrintEigVecs = input%ctrl%tPrintEigVecs
    this%tPrintEigVecsTxt = input%ctrl%tPrintEigVecsTxt

    this%tPrintForces = input%ctrl%tPrintForces
    this%tForces = input%ctrl%tForces .or. this%tPrintForces
    this%isLinResp = allocated(input%ctrl%lrespini)
    if (this%isLinResp) then
      allocate(this%linearResponse)
    end if

    if (this%isLinResp) then
      allocate(this%dQAtomEx(this%nAtom))
      this%dQAtomEx(:) = 0.0_dp
    end if

    if (allocated(input%ctrl%customOccAtoms) .and. this%isLinResp) then
      call error("Custom occupation not compatible with linear response")
    end if

    ! DFTB related variables if multiple determinants are used
    call TDftbDeterminants_init(this%deltaDftb, input%ctrl%isNonAufbau, input%ctrl%isSpinPurify,&
        & input%ctrl%isGroundGuess, this%nEl, this%dftbEnergy)

    if (this%tForces) then
      this%tCasidaForces = input%ctrl%tCasidaForces
    else
      this%tCasidaForces = .false.
    end if
    this%forceType = input%ctrl%forceType
    if (.not. this%tSccCalc .and. input%ctrl%forceType /= forceTypes%orig) then
      call error("Invalid force evaluation method for non-SCC calculations.")
    end if
    if (this%forceType == forceTypes%dynamicT0 .and. this%tempElec > minTemp) then
       call error("This ForceEvaluation method requires the electron temperature to be zero")
     end if
     if (this%isLinResp) then
       tRequireDerivator = (this%tForces .or. input%ctrl%lrespini%tNaCoupling)
     else
       tRequireDerivator = this%tForces
     end if
     if (.not. tRequireDerivator .and. this%isElecDyn) then
       tRequireDerivator = input%ctrl%elecDynInp%tIons
     end if
     if (tRequireDerivator) then
      select case(input%ctrl%iDerivMethod)
      case (diffTypes%finiteDiff)
        ! set step size from input
        if (input%ctrl%deriv1stDelta < epsilon(1.0_dp)) then
          write(tmpStr, "(A,E12.4)") 'Too small value for finite difference step :',&
              & input%ctrl%deriv1stDelta
          call error(tmpStr)
        end if
        call NonSccDiff_init(this%nonSccDeriv, diffTypes%finiteDiff, input%ctrl%deriv1stDelta)
      case (diffTypes%richardson)
        call NonSccDiff_init(this%nonSccDeriv, diffTypes%richardson)
      end select
    end if

    call this%getDenseDescCommon()

    if (this%electronicSolver%isElsiSolver) then
      @:ASSERT(this%parallelKS%nLocalKS == 1)

      if (input%ctrl%parallelOpts%nGroup /= this%nIndepSpin * this%nKPoint) then
        if (this%nSpin == 2) then
          write(tmpStr, "(A,I0,A,I0,A)")"ELSI solvers require as many groups as spin and k-point&
              & combinations. There are ", this%nIndepSpin * this%nKPoint, " spin times k-point&
              & combinations and ", input%ctrl%parallelOpts%nGroup, " groups"
        else
          write(tmpStr, "(A,I0,A,I0,A)")"ELSI solvers require as many groups as k-points. There&
              & are ", this%nIndepSpin * this%nKPoint, " k-points and ",&
              & input%ctrl%parallelOpts%nGroup, " groups"
        end if
        call error(tmpStr)
      end if

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
    ! Check for incompatible options if this is a transport calculation
    if (this%transpar%nCont > 0 .or. this%isAContactCalc) then
      if (allocated(this%dispersion)) then
        call error ("Dispersion interactions are not currently available for transport&
            & calculations")
      end if
      if (this%nSpin > 2) then
        call error("Non-collinear spin polarization disabled for transport calculations at the&
            & moment.")
      end if
      if (this%tExtChrg) then
        call error("External charges temporarily disabled for transport calculations&
            & (electrostatic gates are available).")
      end if
      if (this%t3rdFull .or. this%t3rd) then
        call error ("Third order DFTB is not currently available for transport calculations")
      end if
      if (this%isHybridXc) then
        call error("Range separated calculations do not yet work with transport calculations")
      end if
      if (this%isElecConstr) then
        call error("Constrained DFTB calculations do not yet support electron transport.")
      end if
    end if
  #:endif

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
      this%indMovedAtom = input%ctrl%indMovedAtom
      if (allocated(input%ctrl%indDerivAtom)) then
        this%indDerivAtom = input%ctrl%indDerivAtom
      end if
    else
      allocate(this%indMovedAtom(0))
      allocate(this%indDerivAtom(0))
    end if

    if (allocated(input%ctrl%geoOpt)) then
      if (this%tHelical) then
        call error("GeometryOptimisation driver currently does not support helical geometries")
      end if

      allocate(this%filter)
      call TFilter_init(this%filter, input%ctrl%geoOpt%filter, this%coord0, this%latVec)
      call createOptimizer(input%ctrl%geoOpt%optimiser, this%filter%getDimension(), this%geoOpt)
      this%optTol = input%ctrl%geoOpt%tolerance
      allocate(this%gcurr(this%filter%getDimension()))
      allocate(this%displ(this%filter%getDimension()))
      this%gcurr(:) = 0.0_dp
      this%displ(:) = 0.0_dp
      this%elast = 0.0_dp
      this%nGeoSteps = input%ctrl%geoOpt%nGeoSteps
      this%geoOutFile = input%ctrl%geoOpt%outFile
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
        ! optimisation uses scaling factor of unit cell
        call reset(this%pGeoLatOpt,&
            & (/1.0_dp,0.0_dp,0.0_dp,0.0_dp,0.0_dp,0.0_dp,0.0_dp,0.0_dp,0.0_dp/))
      else if (this%tLatOptFixAng) then
        ! optimisation uses scaling factor of lattice vectors
        call reset(this%pGeoLatOpt,&
            & (/1.0_dp,1.0_dp,1.0_dp,0.0_dp,0.0_dp,0.0_dp,0.0_dp,0.0_dp,0.0_dp/))
      else
        call reset(this%pGeoLatOpt, reshape(this%latVec, (/ 9 /)) )
      end if
    end if

    if (.not.(this%isGeoOpt.or.this%tMD.or.this%tSocket.or.allocated(this%geoOpt))) then
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

      else if (allocated(input%ctrl%dispInp%dftd3)) then
        block
          type(TSDFTD3), allocatable :: dftd3
          allocate(dftd3)
          if (this%tPeriodic) then
            call TSDFTD3_init(dftd3, input%ctrl%dispInp%dftd3, this%nAtom, this%species0, &
                & this%speciesName, this%coord0, this%latVec)
          else
            call TSDFTD3_init(dftd3, input%ctrl%dispInp%dftd3, this%nAtom, this%species0, &
                & this%speciesName, this%coord0)
          end if
          call move_alloc(dftd3, this%dispersion)
        end block

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
        if (allocated(this%reks)) then
          if (input%ctrl%dispInp%dftd4%selfConsistent .and. this%tForces) then
            call error("Calculation of self-consistent dftd4 is not currently compatible with&
                & force calculation in REKS")
          end if
        end if
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
        else if (allocated(this%reks)) then
          call error("Selfconsistent MBD/TS dispersion is blocked from REKS")
          if (this%tForces) then
            call error("Calculation of self-consistent MBD/TS is not currently compatible with&
                & force calculation in REKS")
          end if
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
      else if (allocated(input%ctrl%solvInp%CosmoInp)) then
        if (this%tPeriodic) then
          call createSolvationModel(this%solvation, input%ctrl%solvInp%CosmoInp, &
              & this%nAtom, this%species0, this%speciesName, this%latVec)
        else
          call createSolvationModel(this%solvation, input%ctrl%solvInp%CosmoInp, &
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

      call init_TScaleExtEField(this%eFieldScaling, this%solvation,&
          & input%ctrl%isSolvatedFieldRescaled)

      if (allocated(this%eField)) then
        if (allocated(this%eField%EFieldStrength)) then
          this%eField%EFieldStrength =&
              & this%eFieldScaling%scaledExtEField(this%eField%EFieldStrength)
        end if
      end if

    end if

    if (allocated(this%halogenXCorrection)) then
      this%cutOff%mCutOff = max(this%cutOff%mCutOff, this%halogenXCorrection%getRCutOff())
    end if

    if (allocated(input%ctrl%elecConstraintInp)) then
      call this%ensureConstrainedDftbReqs(input%ctrl%elecConstraintInp)
      allocate(this%elecConstraint)
      call TElecConstraint_init(this%elecConstraint, input%ctrl%elecConstraintInp, this%orb)
    end if

    this%tDipole = this%tMulliken
    if (this%tDipole) then
      block
        logical :: isDipoleDefined
        isDipoleDefined = .true.
        if (abs(input%ctrl%nrChrg) > epsilon(0.0_dp)) then
          call warning("Dipole printed for a charged system : origin dependent quantity")
          isDipoleDefined = .false.
        end if
        if (this%tPeriodic .or. this%tHelical) then
          call warning("Dipole printed for extended system : value printed is not well defined")
          isDipoleDefined = .false.
        end if
        if (.not. isDipoleDefined) then
          write(this%dipoleMessage, "(A)") "Warning: dipole moment is not defined absolutely!"
        else
          write(this%dipoleMessage, "(A)") ""
        end if
      end block
    end if

    if (this%tMulliken) then
      this%tPrintNetAtomCharges = input%ctrl%tPrintNetAtomCharges
      this%tNetAtomCharges = input%ctrl%tNetAtomCharges .or. input%ctrl%tPrintNetAtomCharges
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
      if (allocated(this%solvation)) then
        if (this%solvation%isEFieldModified()) then
          call error("Electrostatic potentials not currently available in the presence of a solvent&
              & which modifies the electrostatics")
        end if
      end if
      allocate(this%electrostatPot)
      call TElStatPotentials_init(this%electrostatPot, input%ctrl%elStatPotentialsInp,&
          & allocated(this%eField) .or. this%tExtChrg, this%hamiltonianType, errStatus)
      if (errStatus%hasError()) then
        call error(errStatus%message)
      end if
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

    this%doPerturbation = allocated(input%ctrl%perturbInp)
    this%doPerturbEachGeom = this%tDerivs .and. this%doPerturbation ! needs work

    if (this%doPerturbation .or. this%doPerturbEachGeom) then

      this%perturbSccTol = input%ctrl%perturbInp%perturbSccTol
      this%maxPerturbIter = input%ctrl%perturbInp%maxPerturbIter
      this%isPerturbConvRequired = input%ctrl%perturbInp%isPerturbConvRequired

      allocate(this%response)
      call TResponse_init(this%response, responseSolverTypes%spectralSum, this%tFixEf,&
          & input%ctrl%perturbInp%tolDegenDFTBPT, input%ctrl%perturbInp%etaFreq)

      this%isEResp = allocated(input%ctrl%perturbInp%dynEFreq)
      if (this%isEResp) then
        call move_alloc(input%ctrl%perturbInp%dynEFreq, this%dynRespEFreq)
        if (this%isHybridXc .and. any(this%dynRespEFreq /= 0.0_dp)) then
          call error("Finite frequency hybrid xc-functional calculation not currently supported.")
        end if
      end if

      this%isKernelResp = allocated(input%ctrl%perturbInp%dynKernelFreq)
      if (this%isKernelResp) then
        call move_alloc(input%ctrl%perturbInp%dynKernelFreq, this%dynKernelFreq)
        if (this%isHybridXc .and. any(this%dynKernelFreq /= 0.0_dp)) then
          call error("Finite frequency hybrid xc-functional calculation not currently supported.")
        end if
        this%isRespKernelRPA = input%ctrl%perturbInp%isRespKernelRPA
        if (.not.this%isRespKernelRPA .and. .not.this%tSccCalc) then
          call error("RPA option only relevant for SCC calculations of response kernel")
        end if
      end if

      if (this%iDistribFn /= fillingTypes%Fermi) then
        call error("Choice of filling function currently incompatible with perturbation&
            & calculations")
      end if

      if (this%tNegf) then
        call error("Currently the perturbation expressions for NEGF are not implemented")
      end if
      if (.not. this%electronicSolver%providesEigenvals) then
        call error("Perturbation expressions currently require eigenvalues and eigenvectors")
      end if

      if (this%t3rd) then
        call error("Only full 3rd order currently supported for perturbation")
      end if
      if (allocated(this%reks)) then
        call error("REKS not currently supported for perturbation")
      end if
      if (this%deltaDftb%isNonAufbau) then
        call error("Delta-DFTB not currently supported for perturbation")
      end if
      if (allocated(this%solvation)) then
        call error("Solvation not currently implemented for perturbation")
      end if
      if (allocated(this%dispersion)) then
        call error("Dispersion (particularly charge dependent) not currently implemented for&
            & perturbation")
      end if

      if (this%isEResp) then
        allocate(this%polarisability(3, 3, size(this%dynRespEFreq)))
        this%polarisability(:,:,:) = 0.0_dp
        if (input%ctrl%tWriteBandDat) then
          ! only one frequency at the moment if dynamic!
          allocate(this%dEidE(this%denseDesc%fullSize, this%nKpoint, this%nIndepSpin, 3))
          this%dEidE(:,:,:,:) = 0.0_dp
        end if
        ! only one frequency at the moment if dynamic!
        allocate(this%dqOut(this%orb%mOrb, this%nAtom, this%nSpin, 3))
        this%dqOut(:,:,:,:) = 0.0_dp
      end if

    else
      this%isEResp = .false.
      this%isKernelResp = .false.
    end if

    ! turn on if LinResp and RangSep turned on, no extra input required for now
    this%isRS_LinResp = this%isLinResp .and. this%isHybridXc

    if (this%isLinResp) then

      ! input checking for linear response
      if (input%ctrl%lrespini%iLinRespSolver==linrespSolverTypes%Arpack .and. .not.withArpack) then
        call error("This binary has been compiled without support for linear response&
            & calculations using the Arpack solver.")
      end if
      isOnsiteCorrected = allocated(this%onSiteElements)
      call ensureLinRespConditions(this%tSccCalc, this%t3rd .or. this%t3rdFull, this%tRealHS,&
          & this%tPeriodic, this%tCasidaForces, this%solvation, this%isRS_LinResp, this%nSpin,&
          & this%tHelical, this%tSpinOrbit, allocated(this%dftbU), this%tempElec,&
          & isOnsiteCorrected, input)

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
        input%ctrl%lrespini%HubbardU(iSp) = hubbU(homoLoc(1), iSp)
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
          & input%ctrl%lrespini%tWriteDensityMatrix .or. input%ctrl%lrespini%tNaCoupling)

      if (allocated(this%onSiteElements) .and. this%tLinRespZVect) then
        call error("Excited state property evaluation currently incompatible with onsite&
            & corrections")
      end if

      call LinResp_init(this%linearResponse, input%ctrl%lrespini, this%nAtom, this%nEl(1),&
          & this%nSpin, this%onSiteElements)

    end if

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
    #:for VAR, ERR in [("tShellResolved", "shell resolved hamiltonians"),&
      & ("tDampH", "H damping")]
      if (input%ctrl%${VAR}$) then
        call error("PP-RPA does not support ${ERR}$")
      end if
    #:endfor
      if (allocated(input%ctrl%h5Input)) then
        call error("PP-RPA does not support H5 correction")
      end if
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
        call error("Choice of filling function incompatible with XLBOMD")
      end if
      if (this%tExtChrg .or. this%isExtField) then
        call error("External fields currently disabled for XLBOMD calculations")
      end if
      if (this%hamiltonianType /= hamiltonianTypes%dftb) then
        call error("XLBOMD calculations currently only supported for the DFTB hamiltonian")
      end if
      if (allocated(this%solvation)) then
        call error("XLBOMD does not work with solvation models yet!")
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
      call create(this%derivDriver, tmp3Coords, size(this%indDerivAtom), input%ctrl%deriv2ndDelta,&
          &this%tDipole, this%doPerturbEachGeom)
      this%coord0(:,this%indMovedAtom) = tmp3Coords
      deallocate(tmp3Coords)
      this%nGeoSteps = 2 * 3 * this%nMovedAtom - 1
    end if

    this%tReadChrg = input%ctrl%tReadChrg
    if (this%tReadChrg .and. this%deltaDftb%isNonAufbau) then
      call error("Charge restart not currently supported for Delta DFTB")
    end if

    this%tReadShifts = input%ctrl%tReadShifts
    this%tWriteShifts = input%ctrl%tWriteShifts

    ! Both temporarily removed until debugged:
    @:ASSERT(.not. this%tReadShifts)
    @:ASSERT(.not. this%tWriteShifts)

    this%tReadChrgAscii = input%ctrl%tReadChrgAscii
    this%tWriteChrgAscii = input%ctrl%tWriteChrgAscii
    this%tSkipChrgChecksum = input%ctrl%tSkipChrgChecksum .or. this%tNegf

  #:if WITH_SCALAPACK
      associate (blacsOpts => input%ctrl%parallelOpts%blacsOpts)
        call getDenseDescBlacs(env, blacsOpts%blockSize, blacsOpts%blockSize, this%denseDesc,&
            & this%isSparseReorderRequired)
      end associate
  #:endif

    if (allocated(this%symNeighbourList)) then
      if ((.not. this%tReadChrg) .and. (this%tPeriodic) .and. this%tPeriodic) then
        this%supercellFoldingMatrix = input%ctrl%supercellFoldingMatrix
        this%supercellFoldingDiag = input%ctrl%supercellFoldingDiag
      end if
      if (this%isHybridXc) then
        call ensureHybridXcReqs(this, input%ctrl%tShellResolved, input%ctrl%hybridXcInp)
        if (.not. this%tReadChrg) then
          if (this%tPeriodic .and. this%tRealHS) then
            ! Periodic system (Gamma-point only), dense Hamiltonian and overlap are real-valued
            call getHybridXcCutOff_gamma(this%cutOff, input%geom%latVecs,&
                & input%ctrl%hybridXcInp%cutoffRed,&
                & gSummationCutoff=input%ctrl%hybridXcInp%gSummationCutoff,&
                & gammaCutoff=input%ctrl%hybridXcInp%gammaCutoff)
          elseif (.not. this%tRealHS) then
            ! Dense Hamiltonian and overlap are complex-valued (general k-point case)
            call getHybridXcCutOff_kpts(this%cutOff, input%geom%latVecs,&
                & input%ctrl%hybridXcInp%cutoffRed, this%supercellFoldingDiag,&
                & gammaCutoff=input%ctrl%hybridXcInp%gammaCutoff,&
                & wignerSeitzReduction=input%ctrl%hybridXcInp%wignerSeitzReduction,&
                & gSummationCutoff=input%ctrl%hybridXcInp%gSummationCutoff)
          end if
        end if

        ! Non-periodic system (cluster)
        if (.not. this%tPeriodic) then
          call getHybridXcCutOff_cluster(this%cutOff, input%ctrl%hybridXcInp%cutoffRed)
        end if
      end if

    #:if WITH_SCALAPACK
      call scalafx_getlocalshape(env%blacs%orbitalGrid, this%denseDesc%blacsOrbSqr, nLocalRows,&
          & nLocalCols)
    #:else
      nLocalRows = this%denseDesc%fullSize
      nLocalCols = this%denseDesc%fullSize
    #:endif

      if (this%isHybridXc) then
        ! allocation is necessary to hint "initializeCharges" what information to extract
        call reallocateHybridXc(this, nLocalRows, nLocalCols,&
            & size(this%parallelKS%localKS, dim=2))
      end if

    end if

    call this%initializeCharges(errStatus, input%ctrl%initialSpins, input%ctrl%initialCharges)
    if (errStatus%hasError()) then
      call error(errStatus%message)
    end if

    ! When restarting and reading charges from charges.bin, the supercell folding matrix of the
    ! initial run is only known after invoking this%initializeCharges(). Inferring the Coulomb
    ! truncation cutoff, therefore calling getHybridXcCutoff(), needs this information.
    if (this%isHybridXc .and. this%tReadChrg) then
      ! First, check if supercell folding matrix is identical to previous run, if specified in input
      if (allocated(input%ctrl%supercellFoldingMatrix)) then
        if (any(abs(input%ctrl%supercellFoldingMatrix&
            & - this%supercellFoldingMatrix) > 1e-06_dp)) then
          write(tmpStr, "(A,3I5,A,3I5,A,3I5,A,3F10.6)")&
              & 'Error while processing k-point sampling for range-separated run.'&
              & // NEW_LINE('A')&
              & // '   When restarting, only identical k-point samplings to previous run are'&
              & // NEW_LINE('A') // '   supported. In this case this would correspond to the&
              & following supercell' // NEW_LINE('A') // '   folding matrix:'&
              & // NEW_LINE('A'),&
              & nint(this%supercellFoldingMatrix(1, 1:3)), NEW_LINE('A'),&
              & nint(this%supercellFoldingMatrix(2, 1:3)), NEW_LINE('A'),&
              & nint(this%supercellFoldingMatrix(3, 1:3)), NEW_LINE('A'),&
              & this%supercellFoldingMatrix(:, 4)
          call error(trim(tmpStr))
        end if
      end if
      if (this%tPeriodic .and. this%tRealHS) then
        ! Periodic system (Gamma-point only), dense Hamiltonian and overlap are real-valued
        call getHybridXcCutOff_gamma(this%cutOff, input%geom%latVecs,&
            & input%ctrl%hybridXcInp%cutoffRed,&
            & gammaCutoff=input%ctrl%hybridXcInp%gammaCutoff)
      elseif (.not. this%tRealHS) then
        ! Dense Hamiltonian and overlap are complex-valued (general k-point case)
        call getHybridXcCutOff_kpts(this%cutOff, input%geom%latVecs,&
            & input%ctrl%hybridXcInp%cutoffRed, this%supercellFoldingDiag,&
            & gammaCutoff=input%ctrl%hybridXcInp%gammaCutoff,&
            & wignerSeitzReduction=input%ctrl%hybridXcInp%wignerSeitzReduction,&
            & gSummationCutoff=input%ctrl%hybridXcInp%gSummationCutoff)
      end if
    end if

    if (this%isHybridXc) then
      call THybridXcFunc_init(this%hybridXc, this%nAtom, this%species0, hubbU(1, :),&
          & input%ctrl%hybridXcInp%screeningThreshold, input%ctrl%hybridXcInp%omega,&
          & input%ctrl%hybridXcInp%camAlpha, input%ctrl%hybridXcInp%camBeta,&
          & this%tSpin, allocated(this%reks), input%ctrl%hybridXcInp%hybridXcAlg,&
          & input%ctrl%hybridXcInp%hybridXcType, input%ctrl%hybridXcInp%gammaType, this%tPeriodic,&
          & this%tRealHS, errStatus, coeffsDiag=this%supercellFoldingDiag,&
          & gammaCutoff=this%cutOff%gammaCutoff,&
          & gSummationCutoff=this%cutOff%gSummationCutoff,&
          & wignerSeitzReduction=this%cutOff%wignerSeitzReduction,&
          & latVecs=input%geom%latVecs)
      if (errStatus%hasError()) then
        call error(errStatus%message)
      end if
      ! now all information is present to properly allocate density matrices and associate pointers
      call reallocateHybridXc(this, nLocalRows, nLocalCols, size(this%parallelKS%localKS, dim=2))
      ! reset number of mixer elements, so that there is enough space for density matrices
      if (this%tRealHS) then
        this%nMixElements = size(this%densityMatrix%deltaRhoIn)
      else
        this%nMixElements = size(this%densityMatrix%deltaRhoInCplxHS)
      end if
    end if

    ! Initialise images (translations)
    if (this%tPeriodic .or. this%tHelical) then
      call getCellTranslations(this%cellVec, this%rCellVec, this%latVec, this%invLatVec,&
          & this%cutOff%mCutOff)
    else
      allocate(this%cellVec(3, 1))
      allocate(this%rCellVec(3, 1))
      this%cellVec(:, 1) = [0.0_dp, 0.0_dp, 0.0_dp]
      this%rCellVec(:, 1) = [0.0_dp, 0.0_dp, 0.0_dp]
    end if

    ! Initialize neighbourlist(s)
    allocate(this%neighbourList)
    call TNeighbourlist_init(this%neighbourList, this%nAtom, nInitNeighbour)
    allocate(this%nNeighbourSK(this%nAtom))
    if (allocated(this%symNeighbourList)) then
      allocate(this%symNeighbourList%neighbourList)
      call TNeighbourlist_init(this%symNeighbourList%neighbourList, this%nAtom, nInitNeighbour)
      if (this%isHybridXc) then
        allocate(this%nNeighbourCam(this%nAtom))
        allocate(this%nNeighbourCamSym(this%nAtom))
      end if
    end if

    ! Set various options
    this%tWriteAutotest = env%tGlobalLead .and. input%ctrl%tWriteTagged
    this%tWriteDetailedXML = env%tGlobalLead .and. input%ctrl%tWriteDetailedXML
    this%tWriteResultsTag = env%tGlobalLead .and. input%ctrl%tWriteResultsTag
    this%tWriteCharges = env%tGlobalLead .and. input%ctrl%tWriteCharges
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
  #:endif

    if (this%tNegf) then
      if (allocated(this%dispersion)) then
        call error("Dispersion not currently available with transport calculations")
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

    if (allocated(this%reks)) then
      call checkReksConsistency(input%ctrl%reksInp, this%solvation, this%onSiteElements,&
          & this%kPoint, this%nEl, this%nKPoint, this%tSccCalc, this%tSpin, this%tSpinOrbit,&
          & allocated(this%dftbU), this%isExtField, this%isLinResp, this%tPeriodic, this%tLatOpt,&
          & this%tReadChrg, this%tPoisson, input%ctrl%tShellResolved)
      ! here, this%nSpin changes to 2 for REKS
      call TReksCalc_init(this%reks, input%ctrl%reksInp, this%electronicSolver, this%orb,&
          & this%spinW, this%nEl, input%ctrl%extChrg, input%ctrl%extChrgBlurWidth,&
          & this%hamiltonianType, this%nSpin, this%nExtChrg, this%t3rd.or.this%t3rdFull,&
          & this%isHybridXc, allocated(this%dispersion), this%isQNetAllocated, this%tForces,&
          & this%tPeriodic, this%tStress, this%tDipole)
    end if

    call this%initDetArrays(nLocalRows, nLocalCols)
    call this%initArrays(env, input)

  #:if WITH_TRANSPORT
    if (this%tUpload) then
      ! check geometry details are consistent with transport with contacts
      call checkTransportRanges(this%nAtom, input%transpar)
    end if

    if (this%isAContactCalc) then
      ! geometry is reduced to contacts only
      allocate(this%iAtInCentralRegion(this%nAtom))
      ! for storage of the electrostatic potential in the contact
      allocate(this%potential%coulombShell(this%orb%mShell,this%nAtom,1))
      this%potential%coulombShell(:,:,:) = 0.0_dp
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

    if (allocated(this%eField)) then
      if (allocated(input%ctrl%atomicExtPotential)) then
        ii = 0
        jj = 0
        if (allocated(input%ctrl%atomicExtPotential%iAtOnSite)) then
          ii = size(input%ctrl%atomicExtPotential%iAtOnSite)
          if (ii > 1) then
            if (any(input%ctrl%atomicExtPotential%iAtOnSite < 1) .or.&
                & any(input%ctrl%atomicExtPotential%iAtOnSite > size(this%iAtInCentralRegion))) then
              call error("Net potential atom(s) outside of range of real atoms")
            end if
            if (isRepeated(input%ctrl%atomicExtPotential%iAtOnSite)) then
              call error("Repeating atom(s) for net potentials")
            end if
          end if
        end if
        if (allocated(input%ctrl%atomicExtPotential%iAt)) then
          jj = size(input%ctrl%atomicExtPotential%iAt)
          if (jj > 1) then
            if (any(input%ctrl%atomicExtPotential%iAt < 1) .or.&
                & any(input%ctrl%atomicExtPotential%iAt > size(this%iAtInCentralRegion))) then
              call error("Gross potential atom(s) outside of range of real atoms")
            end if
            if (isRepeated(input%ctrl%atomicExtPotential%iAt)) then
              call error("Repeating atom(s) for gross potentials")
            end if
          end if
        end if
        allocate(this%eField%atomicSites(ii+jj))
        allocate(this%eField%atomicPotential(ii+jj))
        allocate(this%eField%atomicOnSite(ii+jj))
        if (allocated(input%ctrl%atomicExtPotential%iAtOnSite)) then
          this%eField%atomicSites(:ii) = input%ctrl%atomicExtPotential%iAtOnSite
          this%eField%atomicPotential(:ii) = input%ctrl%atomicExtPotential%VextOnSite
          this%eField%atomicOnSite(:ii) = .true.
        end if
        if (allocated(input%ctrl%atomicExtPotential%iAt)) then
          this%eField%atomicSites(ii+1:) = input%ctrl%atomicExtPotential%iAt
          this%eField%atomicPotential(ii+1:) = input%ctrl%atomicExtPotential%Vext
          this%eField%atomicOnSite(ii+1:) = .false.
        end if
      end if
    end if

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
    if (.not. (this%isHybridXc .and. this%tRealHS .and. this%tPeriodic)) then
      write(stdOut, "('BLACS orbital grid size:', T30, I0, ' x ', I0)")env%blacs%orbitalGrid%nRow,&
          & env%blacs%orbitalGrid%nCol
      write(stdOut, "('BLACS atom grid size:', T30, I0, ' x ', I0)")env%blacs%atomGrid%nRow,&
          & env%blacs%atomGrid%nCol
    end if
  #:endif

    if (tRandomSeed) then
      write(stdOut, "(A,':',T30,I0)") "Chosen random seed", iSeed
    else
      write(stdOut, "(A,':',T30,I0)") "Specified random seed", iSeed
    end if

    call checkStackSize(env)

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
      tGeoOptRequiresEgy = .true.
      select case (input%ctrl%iGeoOpt)
      case (geoOptTypes%steepestDesc)
        write(stdOut, "('Mode:',T30,A)")'Steepest descent' // trim(strTmp)
      case (geoOptTypes%conjugateGrad)
        write(stdOut, "('Mode:',T30,A)") 'Conjugate gradient relaxation' // trim(strTmp)
      case (geoOptTypes%diis)
        write(stdOut, "('Mode:',T30,A)") 'Modified gDIIS relaxation' // trim(strTmp)
        tGeoOptRequiresEgy = .false.
      case (geoOptTypes%lbfgs)
        write(stdout, "('Mode:',T30,A)") 'LBFGS relaxation' // trim(strTmp)
      case (geoOptTypes%fire)
        write(stdout, "('Mode:',T30,A)") 'FIRE relaxation' // trim(strTmp)
        tGeoOptRequiresEgy = .false.
      case default
        call error("Unknown optimisation mode")
      end select
      if (tGeoOptRequiresEgy .neqv. this%electronicSolver%providesFreeEnergy) then
        call warning("This geometry optimisation method requires force related energies for&
            & accurate minimisation.")
      end if
    elseif (this%tDerivs) then
      write(stdOut, "('Mode:',T30,A)") "2nd derivatives calculation"
      write(stdOut, "('Mode:',T30,A)") "Calculated for atoms:"
      write(stdOut, *) this%indDerivAtom
      if (size(this%indDerivAtom) > size(this%indMovedAtom)) then
        write(stdOut, "('Mode:',T30,A)") "Moved atoms:"
        write(stdOut, *) this%indMovedAtom
      end if
    elseif (this%tSocket) then
      write(stdOut, "('Mode:',T30,A)") "Socket controlled calculation"
    else
      write(stdOut, "('Mode:',T30,A)") "Static calculation"
    end if

    if (this%tSccCalc) then
      if (.not. this%tRestartNoSC) then
        write(stdOut, "(A,':',T30,A)") "Self consistent charges", "Yes"
        write(stdOut, "(A,':',T30,E14.6)") "SCC-tolerance", this%sccTol
        write(stdOut, "(A,':',T30,I14)") "Max. scc iterations", this%maxSccIter
      end if
      !if (this%tPeriodic) then
      !  write(stdout, "(A,':',T30,E14.6)") "Ewald alpha parameter", this%scc%getEwaldPar()
      !end if
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
      case(2)
        write(stdOut, "(A,':',T30,A)") "Spin polarisation", "Yes"
      case(4)
        write(stdOut, "(A,':',T30,A)") "Non-collinear calculation", "Yes"
      end select
      if (any(this%electronicSolver%iSolver ==&
          & [electronicSolverTypes%GF,electronicSolverTypes%onlyTransport]) .or. this%tFixEf) then
        ! number of electrons not determined at this stage
      else
        select case (this%nSpin)
        case(1)
          write(stdOut, "(A,':',T30,F12.6,/,A,':',T30,F12.6)") "Nr. of up electrons",&
              & 0.5_dp*this%nEl(1), "Nr. of down electrons", 0.5_dp*this%nEl(1)
        case(2)
          if (this%tSpinSharedEf) then
            write(stdOut, "(A,':',T30,F12.6,/,A,':',T30,F12.6)") "Initial nr. of up electrons",&
                & this%nEl(1), "Initial Nr. of down electrons", this%nEl(2)
          else
            write(stdOut, "(A,':',T30,F12.6,/,A,':',T30,F12.6)") "Nr. of up electrons",&
                & this%nEl(1), "Nr. of down electrons", this%nEl(2)
          end if
        case(4)
          write(stdOut, "(A,':',T30,F12.6)") "Nr. of electrons", this%nEl(1)
        end select
      end if
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
      #:if WITH_MAGMA
        call env%initGpu()
      #:else
        call error("Magma-solver selected, but program was compiled without MAGMA")
      #:endif
    endif

    if(this%isLinResp) then
      select case(input%ctrl%lrespini%iLinRespSolver)
      case (linrespSolverTypes%None)
        call error("Casida solver has not been selected")
      case (linrespSolverTypes%Arpack)
        write(stdOut, "(A,':',T30,A)") "Casida solver", "Arpack"
      case (linrespSolverTypes%Stratmann)
        write(stdOut, "(A,':',T30,A,i4)") "Casida solver", "Stratmann, SubSpace: ",&
            & input%ctrl%lrespini%subSpaceFactorStratmann
      case default
        call error("Unknown Casida solver")
      end select
    end if

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
          write(stdOut, "(A,':',T30,I14)") "Nr. of chrg. vec. in memory", this%maxSccIter
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

    if (allocated(this%tblite)) then
      call writeTBLiteInfo(stdOut, this%tblite)
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

    if (this%tSccCalc .and. .not. allocated(this%reks) .and. .not.this%tRestartNoSC) then
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
      type is (TSDFTD3)
        call writeSDFTD3Info(stdout, o)
      type is (TSimpleDftD3)
        write(stdOut, "(A)") "Using simple DFT-D3 dispersion corrections"
      type is (TDispDftD4)
        call writeDftD4Info(stdOut, o)
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
      if (this%eFieldScaling%isRescaled) then
        write(stdOut, "(A,':',T30,A)")"Solvated fields rescaled", "Yes"
      else
        write(stdOut, "(A,':',T30,A)")"Solvated fields rescaled", "No"
      end if
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
                & hubbU(jj, iSp)
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
      if (allocated(this%scc)) then
        if (any(shortGammaDamp%isDamped)) then
          write(stdOut, "(A,T30,A)") "Damped SCC", "Yes"
          ii = count(shortGammaDamp%isDamped)
          write(strTmp, "(A,I0,A)") "(A,T30,", ii, "(A,1X))"
          write(stdOut, strTmp) "Damped species(s):", pack(this%speciesName,&
              & shortGammaDamp%isDamped)
        end if
      end if

      if (allocated(input%ctrl%h5Input)) then
        write(stdOut, "(A,T30,A)") "H-bond correction:", "H5"
      end if
      if (tHHRepulsion) then
        write(stdOut, "(A,T30,A)") "H-H repulsion correction:", "H5"
      end if
    end if

    if (this%isHybridXc) then
      if (input%ctrl%hybridXcInp%hybridXcType == hybridXcFunc%hyb) then
        write(stdOut, "(A,':',T30,A)") "Global hybrid", "Yes"
        write(stdOut, "(2X,A,':',T30,E14.6)") "Fraction of exchange",&
            & input%ctrl%hybridXcInp%camAlpha
      elseif (input%ctrl%hybridXcInp%hybridXcType == hybridXcFunc%lc) then
        write(stdOut, "(A,':',T30,A)") "Long-range corrected hybrid", "Yes"
        write(stdOut, "(2X,A,':',T30,E14.6)") "Screening parameter omega",&
            & input%ctrl%hybridXcInp%omega
      elseif (input%ctrl%hybridXcInp%hybridXcType == hybridXcFunc%cam) then
        write(stdOut, "(A,':',T30,A)") "CAM range-separated hybrid", "Yes"
        write(stdOut, "(2X,A,':',T30,E14.6)") "Screening parameter omega",&
            & input%ctrl%hybridXcInp%omega
        write(stdOut, "(2X,A,':',T30,E14.6,E14.6)") "CAM parameters alpha/beta",&
            & input%ctrl%hybridXcInp%camAlpha, input%ctrl%hybridXcInp%camBeta
      end if
      if (this%tPeriodic) then
        if (input%ctrl%hybridXcInp%gammaType == hybridXcGammaTypes%full) then
          write(stdOut, "(2X,A,':',T30,2X,A)") "Gamma function", "full"
        elseif (input%ctrl%hybridXcInp%gammaType == hybridXcGammaTypes%mic) then
          write(stdOut, "(2X,A,':',T30,2X,A)") "Gamma function", "minimum image convention"
          write(stdOut, "(2X,A,':',T30,2X,I0,A)") "Wigner-Seitz cell reduction",&
              & this%cutOff%wignerSeitzReduction, " primitive cell(s)"
        elseif (input%ctrl%hybridXcInp%gammaType == hybridXcGammaTypes%truncated) then
          write(stdOut, "(2X,A,':',T30,2X,A)") "Gamma function", "truncated"
        elseif (input%ctrl%hybridXcInp%gammaType == hybridXcGammaTypes%truncatedAndDamped) then
          write(stdOut, "(2X,A,':',T30,2X,A)") "Gamma function", "truncated+poly5zero"
        end if
        if (input%ctrl%hybridXcInp%gammaType /= hybridXcGammaTypes%mic) then
          write(stdOut, "(2X,A,':',T30,E14.6,A)") "G-Summation Cutoff",&
              & this%cutOff%gSummationCutoff, " Bohr"
        end if
        if (input%ctrl%hybridXcInp%gammaType == hybridXcGammaTypes%truncated&
            & .or. input%ctrl%hybridXcInp%gammaType == hybridXcGammaTypes%truncatedAndDamped) then
          write(stdOut, "(2X,A,':',T30,E14.6,A)") "Coulomb Truncation",&
              & this%cutOff%gammaCutoff, " Bohr"
        end if
      end if

      select case(input%ctrl%hybridXcInp%hybridXcAlg)
      case (hybridXcAlgo%neighbourBased)
        write(stdOut, "(2X,A,':',T30,2X,A)") "Screening algorithm", "NeighbourBased"
        write(stdOut, "(2X,A,':',T30,E14.6,A)") "Reduce neighlist cutoff by",&
            & input%ctrl%hybridXcInp%cutoffRed * Bohr__AA, " A"
        if (this%tPeriodic) then
          write(stdOut, "(2X,A,':',T30,E14.6)") "Thresholded to",&
              & input%ctrl%hybridXcInp%screeningThreshold
        end if
      case (hybridXcAlgo%thresholdBased)
        write(stdOut, "(2X,A,':',T30,2X,A)") "Screening algorithm", "Thresholded"
        write(stdOut, "(2X,A,':',T30,E14.6)") "Thresholded to",&
            & input%ctrl%hybridXcInp%screeningThreshold
      case (hybridXcAlgo%matrixBased)
        write(stdOut, "(2X,A,':',T30,2X,A)") "Screening algorithm", "MatrixBased"
      case default
        call error("Unknown hybrid xc-functional screening algorithm")
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
        write(stdOut, "(A,T30,A)") "Force type", "erho with re-diagonalised eigenvalues"
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

    if (this%isExtField) then

      if (allocated(this%eField%EFieldStrength)) then
        if (this%eField%isTDEfield) then
          write(stdOut, "(T30,A)") "External electric field specified"
          write(stdOut, "(A,':',T30,E14.6)") "Angular frequency", this%eField%EfieldOmega
        else
          write(stdOut, "(T30,A)") "External static electric field specified"
        end if
        write(stdOut, "(A,':',T30,E14.6)") "Field strength", this%eField%EFieldStrength
        write(stdOut, "(A,':',T30,3F9.6)") "Direction", this%eField%EfieldVector
        if (this%tPeriodic) then
          call warning("Saw tooth potential used for periodic geometry - make sure there is a&
              & vacuum region!")
        end if
      end if

      if (allocated(input%ctrl%atomicExtPotential)) then
        if (allocated(input%ctrl%atomicExtPotential%iAtOnSite)) then
          write(stdOut, "(A)")'Net on-site potentials at atoms (/ H)'
          do ii = 1, size(input%ctrl%atomicExtPotential%iAtOnSite)
            write(stdOut,"(1X,I6,' : ',E14.6)")input%ctrl%atomicExtPotential%iAtOnSite(ii),&
                & input%ctrl%atomicExtPotential%VextOnSite(ii)
          end do
        end if
        if (allocated(input%ctrl%atomicExtPotential%iAt)) then
          write(stdOut, "(A)")'Gross on-site potentials at atoms (/ H)'
          do ii = 1, size(input%ctrl%atomicExtPotential%iAt)
            write(stdOut,"(1X,I6,' : ',E14.6)")input%ctrl%atomicExtPotential%iAt(ii),&
                & input%ctrl%atomicExtPotential%Vext(ii)
          end do
        end if
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
    if (this%deltaDftb%isNonAufbau .and. this%isElecDyn) then
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
        call error("Sorry, lattice optimisation requires stress tensor evaluation")
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
    if (this%isElecDyn) then

      if (this%t2Component) then
        call error("Electron dynamics is not compatibile with this spinor Hamiltonian")
      end if

      if (input%ctrl%elecDynInp%tIons .and. withMPI) then
        call error("Ion dynamics time propagation does not work with MPI yet")
      end if

      if (.not. this%tRealHS .and. withMpi) then
        call error("Electron dynamics of periodic systems does not work with MPI yet")
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

      if (.not. this%tRealHS .and. input%ctrl%elecDynInp%tBondE) then
        call error("Bond energies during electron dynamics currently requires a real hamiltonian.")
      end if

      if (this%isHybridXc) then
        if (input%ctrl%elecDynInp%spType == tdSpinTypes%triplet) then
          call error("Triplet perturbations currently disabled for electron dynamics with&
              & range-separated functionals")
        end if
        if (input%ctrl%elecDynInp%tForces) then
          call error("Forces for time propagation currently disabled for range-separated")
        end if
        if (input%ctrl%elecDynInp%tIons) then
          call error("Ion dynamics time propagation currently disabled for range-separated")
        end if
      end if

      allocate(this%electronDynamics)

      call TElecDynamics_init(this%electronDynamics, input%ctrl%elecDynInp, this%species0,&
          & this%speciesName, this%tWriteAutotest, autotestTag, randomThermostat, this%mass,&
          & this%nAtom, this%cutOff%skCutoff, this%cutOff%mCutoff, this%atomEigVal,&
          & this%dispersion, this%nonSccDeriv, this%tPeriodic, this%parallelKS, this%tRealHS,&
          & this%kPoint, this%kWeight, this%isHybridXc, this%scc, this%tblite, this%eFieldScaling,&
          & this%hamiltonianType, errStatus, this%denseDesc)
      if (errStatus%hasError()) then
        call error(errStatus%message)
      end if

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

    !> Lattice vectors, stored columnwise
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
  subroutine setEquivalencyRelations(this)

    !> Instance
    class(TDftbPlusMain), intent(inout) :: this

    ! Orbital equivalency for SCC and Spin
    integer, allocatable :: iEqOrbSCC(:,:,:), iEqOrbSpin(:,:,:)
    ! Orbital equivalency for orbital potentials
    integer, allocatable :: iEqOrbDFTBU(:,:,:)
    ! Equivalency for multipolar contributions
    integer, allocatable :: iEqDip(:,:), iEqQuad(:,:)

    this%nIneqOrb = 0
    this%nIneqDip = 0
    this%nIneqQuad = 0
    this%nMixElements = 0

    if (this%tSccCalc) then
      if(.not. allocated(this%iEqOrbitals)) then
        allocate(this%iEqOrbitals(this%orb%mOrb, this%nAtom, this%nSpin))
      endif
      if (.not.allocated(this%iEqDipole)) then
        allocate(this%iEqDipole(this%nDipole, this%nAtom))
      end if
      if (.not.allocated(this%iEqQuadrupole)) then
        allocate(this%iEqQuadrupole(this%nQuadrupole, this%nAtom))
      end if
      allocate(iEqOrbSCC(this%orb%mOrb, this%nAtom, this%nSpin))
      @:ASSERT(allocated(this%uniqHubbU))
      call this%uniqHubbU%getOrbitalEquiv(this%orb, this%species0, iEqOrbSCC)
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

      if (allocated(this%tblite)) then
        allocate(iEqOrbSpin(this%orb%mOrb, this%nAtom, this%nSpin), source=0)
        allocate(iEqOrbDFTBU(this%orb%mOrb, this%nAtom, this%nSpin), source=0)
        allocate(iEqDip(this%nDipole, this%nAtom), source=0)
        allocate(iEqQuad(this%nQuadrupole, this%nAtom), source=0)
        call this%tblite%getOrbitalEquiv(iEqOrbDFTBU, iEqDip, iEqQuad, this%orb, this%species0)
        call OrbitalEquiv_merge(this%iEqOrbitals, iEqOrbDFTBU, this%orb, iEqOrbSpin)
        this%iEqOrbitals(:,:,:) = iEqOrbSpin
        this%iEqDipole(:,:) = iEqDip
        this%iEqQuadrupole(:,:) = iEqQuad
        this%nIneqOrb = maxval(this%iEqOrbitals)
        this%nIneqDip = maxval(this%iEqDipole)
        this%nIneqQuad = maxval(this%iEqQuadrupole)
        this%nMixElements = this%nIneqOrb + this%nIneqDip + this%nIneqQuad
        deallocate(iEqOrbSpin)
        deallocate(iEqOrbDFTBU)
        deallocate(iEqDip)
        deallocate(iEqQuad)
      end if

      if (allocated(this%onSiteElements)) then
        allocate(iEqOrbSpin(this%orb%mOrb, this%nAtom, this%nSpin))
        iEqOrbSpin(:,:,:) = 0
        allocate(iEqOrbDFTBU(this%orb%mOrb, this%nAtom, this%nSpin))
        iEqOrbDFTBU(:,:,:) = 0
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
    else
      !Non-SCC
      this%nIneqOrb = this%nOrb
      this%nMixElements = 0
    end if

  end subroutine setEquivalencyRelations


  !> Initialise partial charges
  !>
  subroutine initializeCharges(this, errStatus, initialSpins, initialCharges)

    !> Instance
    class(TDftbPlusMain), intent(inout) :: this

    !> Error status
    type(TStatus), intent(inout) :: errStatus

    !> Initial spins
    real(dp), optional, intent(in) :: initialSpins(:,:)

    !> Set of atom-resolved atomic charges
    real(dp), optional, intent(in) :: initialCharges(:)

    !> Tolerance in difference between total charge and sum of initial charges
    real(dp), parameter :: deltaChargeTol = 1.e-4_dp

    integer :: iAt, iSp, iSh, ii, jj, iStart, iEnd, iS
    real(dp) :: rTmp
    character(lc) :: message

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

    this%isQNetAllocated = .false.
  #:if WITH_MBD
    if (allocated(this%dispersion)) then
      select type (o=>this%dispersion)
      type is (TDispMbd)
        this%isQNetAllocated = .true.
      end select
    end if
  #:endif
    this%isQNetAllocated = this%isQNetAllocated .or. this%tNetAtomCharges
    if (this%isQNetAllocated) then
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

    if (this%isHybridXc .and. this%tReadChrg) then
      allocate(this%supercellFoldingMatrix(3, 4))
    end if

    ! Charges read from file
    if (this%tReadChrg) then
      if (this%tFixEf .or. this%tSkipChrgChecksum) then
        ! do not check charge or magnetisation from file
        call initQFromFile(this%qInput, fCharges, this%tReadChrgAscii, this%orb, this%qBlockIn,&
            & this%qiBlockIn, this%densityMatrix, this%tRealHS, errStatus,&
            & multipoles=this%multipoleInp, coeffsAndShifts=this%supercellFoldingMatrix)
        @:PROPAGATE_ERROR(errStatus)
      else
        ! check number of electrons in file
        if (this%nSpin /= 2) then
          call initQFromFile(this%qInput, fCharges, this%tReadChrgAscii, this%orb, this%qBlockIn,&
              & this%qiBlockIn, this%densityMatrix, this%tRealHS, errStatus, nEl=sum(this%nEl),&
              & multipoles=this%multipoleInp, coeffsAndShifts=this%supercellFoldingMatrix)
          @:PROPAGATE_ERROR(errStatus)
        else
          ! check magnetisation in addition
          call initQFromFile(this%qInput, fCharges, this%tReadChrgAscii, this%orb, this%qBlockIn,&
              & this%qiBlockIn, this%densityMatrix, this%tRealHS, errStatus, nEl=sum(this%nEl),&
              & magnetisation=this%nEl(1)-this%nEl(2), multipoles=this%multipoleInp,&
              & coeffsAndShifts=this%supercellFoldingMatrix)
          @:PROPAGATE_ERROR(errStatus)
        end if
      end if

      ! Check if obtained supercell folding matrix meets current requirements
      if (this%isHybridXc .and. this%tPeriodic) then
        allocate(this%supercellFoldingDiag(3))
        call checkSupercellFoldingMatrix(this%supercellFoldingMatrix, errStatus,&
            & supercellFoldingDiagOut=this%supercellFoldingDiag)
        @:PROPAGATE_ERROR(errStatus)
      end if

    #:if WITH_TRANSPORT
      if (this%tUpload) then
        call overrideContactCharges(this%qInput, this%chargeUp, this%transpar, this%qBlockIn,&
            & this%blockUp)
      end if
    #:endif

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

    if (allocated(this%reks)) return

    ! Charges not read from file
    notChrgRead: if (.not. this%tReadChrg) then

      if (present(initialCharges)) then
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
        this%qInput(:,:,1) = this%qInput(:,:,1) * sum(this%nEl) / sum(this%qInput(:,:,1))
      end if

      select case (this%nSpin)
      case (1)
        continue
      case (2)
        if (present(initialSpins)) then
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
          if (.not. present(initialSpins)) then
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

    if (allocated(this%multipoleInp%dipoleAtom)) then
      ! FIXME: Assumes we always mix all dipole moments
      this%qInpRed(this%nIneqOrb+1:this%nIneqOrb+this%nIneqDip) =&
          & reshape(this%multipoleInp%dipoleAtom, [this%nIneqDip])
    end if
    if (allocated(this%multipoleInp%quadrupoleAtom)) then
      ! FIXME: Assumes we always mix all quadrupole moments
      this%qInpRed(this%nIneqOrb+this%nIneqDip+1:this%nIneqOrb+this%nIneqDip+this%nIneqQuad) =&
          & reshape(this%multipoleInp%quadrupoleAtom, [this%nIneqQuad])
    end if

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


  ! Assign reference charge arrays, q0 and qShell0
  subroutine initReferenceCharges(species0, orb, referenceN0, nSpin, q0, qShell0, customOccAtoms,&
      & customOccFillings)

    !> Species of central cell atoms
    integer, intent(in) :: species0(:)

    !> Atomic orbital data
    type(TOrbitals), intent(in) :: orb

    !> Reference neutral atom charges for each species
    real(dp), intent(in) :: referenceN0(:,:)

    !> Number of spin channels
    integer, intent(in) :: nSpin

    !> Reference 'neutral' charges for each atom orbital
    real(dp), allocatable, intent(out) :: q0(:,:,:)

    !> Shell resolved neutral reference
    real(dp), allocatable, intent(out) :: qShell0(:,:)

    !> Array of lists of atoms where the 'neutral' shell occupation is modified
    type(TWrappedInt1), optional, intent(in) :: customOccAtoms(:)

    !> Modified occupations for shells of the groups atoms in customOccAtoms
    real(dp), optional, intent(in) :: customOccFillings(:,:)

    integer :: nAtom

    @:ASSERT(present(customOccAtoms) .eqv. present(customOccFillings))

    nAtom = size(orb%nOrbAtom)
    allocate(q0(orb%mOrb, nAtom, nSpin))
    q0(:,:,:) = 0.0_dp
    if (present(customOccAtoms)) then
      call applyCustomReferenceOccupations(customOccAtoms, customOccFillings, species0, orb,&
          & referenceN0, q0)
    else
      call initQFromShellChrg(q0, referenceN0, species0, orb)
    end if

    allocate(qShell0(orb%mShell, nAtom))
    call updateReferenceShellCharges(qShell0, q0, orb, nAtom, species0)

  end subroutine initReferenceCharges


  !> Updates the reference shell charges
  subroutine updateReferenceShellCharges(qShell0, q0, orb, nAtom, species0)

    !> Shell resolved neutral reference
    real(dp), intent(out) :: qShell0(:,:)

    !> Reference 'neutral' charges for each atom's orbitals
    real(dp), intent(in) :: q0(:,:,:)

    !> Atomic orbital data
    type(TOrbitals), intent(in) :: orb

    !> Atoms in the system
    integer, intent(in) :: nAtom

    !> Species of central cell atoms
    integer, intent(in) :: species0(:)

    integer :: iAt, iSp, iSh

    do iAt = 1, nAtom
      iSp = species0(iAt)
      do iSh = 1, orb%nShell(iSp)
        qShell0(iSh, iAt) = sum(q0(orb%posShell(iSh, iSp) : orb%posShell(iSh + 1, iSp) - 1, iAt, 1))
      end do
    end do

  end subroutine updateReferenceShellCharges


  !> Set number of electrons
  subroutine initElectronNumber(q0, nrChrg, nSpinPol, nSpin, orb, nEl0, nEl)

    !> Reference 'neutral' charges for each atom's orbitals
    real(dp), intent(in) :: q0(:,:,:)

    !> Total charge
    real(dp), intent(in) :: nrChrg

    !> Spin polarisation
    real(dp), intent(in) :: nSpinPol

    !> Number of spin components, 1 is unpolarised, 2 is polarised, 4 is noncolinear / spin-orbit
    integer, intent(in) :: nSpin

    !> Atomic orbital data
    type(TOrbitals), intent(in) :: orb

    !> Nr. of all electrons if neutral
    real(dp), intent(out) :: nEl0

    !> Nr. of electrons
    real(dp), allocatable, intent(out) :: nEl(:)

    nEl0 = sum(q0(:,:,1))
    if (abs(nEl0 - nint(nEl0)) < elecTolMax) then
      nEl0 = nint(nEl0)
    end if

    if (nSpin == 4) then
      allocate(nEl(1))
    else
      allocate(nEl(nSpin))
    end if
    nEl(:) = 0.0_dp
    if (nSpin == 1 .or. nSpin == 4) then
      nEl(1) = nEl0 - nrChrg
      if(ceiling(nEl(1)) > 2.0_dp * orb%nOrb) then
        call error("More electrons than basis functions!")
      end if
    else
      nEl(1) = 0.5_dp * (nEl0 - nrChrg + nSpinPol)
      nEl(2) = 0.5_dp * (nEl0 - nrChrg - nSpinPol)
      if (any(ceiling(nEl) > orb%nOrb)) then
        call error("More electrons than basis functions!")
      end if
    end if

    if (.not. all(nEl >= 0.0_dp)) then
      call error("Less than 0 electrons!")
    end if

  end subroutine initElectronNumber


  ! Set up reference population
  subroutine initReferencePopulation_(input, orb, hamiltonianType, referenceN0)
    type(TInputData), intent(in) :: input
    type(TOrbitals), intent(in) :: orb
    integer, intent(in) :: hamiltonianType
    real(dp), allocatable, intent(out) :: referenceN0(:,:)

    integer :: nSpecies

    nSpecies = size(orb%nOrbSpecies)
    allocate(referenceN0(orb%mShell, nSpecies))
    select case (hamiltonianType)
    case default
      call error("Invalid Hamiltonian")
    case(hamiltonianTypes%dftb)
      referenceN0(:,:) = input%slako%skOcc(1:orb%mShell, :)
    case(hamiltonianTypes%xtb)
      ! It's an xTB Hamiltonian masquerading as a DFTB one (for now)
      referenceN0(:,:) = input%slako%skOcc(1:orb%mShell, :)
    end select

  end subroutine initReferencePopulation_


  !> Initializes electron filling related variables
  subroutine initElectronFilling_(input, nSpin, Ef, iDistribFn, tempElec, tFixEf, tSetFillingTemp,&
      & tFillKSep)

    !> Holds the parsed input data.
    type(TInputData), intent(in) :: input

    !> Number of spin components, 1 is unpolarised, 2 is polarised, 4 is noncolinear / spin-orbit
    integer, intent(in) :: nSpin

    !> Fermi energy for each spin
    real(dp), allocatable, intent(out) :: Ef(:)

    !> Choice of electron distribution function, defaults to Fermi
    integer, intent(out) :: iDistribFn

    !> Electron temperature
    real(dp), intent(out) :: tempElec

    !> Fix Fermi energy at specified value
    logical, intent(out) :: tFixEf

    !> Filling temp, as updated by MD.
    logical, intent(out) :: tSetFillingTemp

    !> If K points should filled separately
    logical, intent(out) :: tFillKSep

    if (nSpin == 4) then
      allocate(Ef(1))
    else
      allocate(Ef(nSpin))
    end if

    iDistribFn = input%ctrl%iDistribFn
    tempElec = input%ctrl%tempElec

    tFixEf = input%ctrl%tFixEf
    if (allocated(input%ctrl%Ef)) then
      Ef(:) = input%ctrl%Ef
    else
      Ef(:) = 0.0_dp
    end if
    tSetFillingTemp = input%ctrl%tSetFillingTemp
    tFillKSep = input%ctrl%tFillKSep

  end subroutine initElectronFilling_


  ! Set up Hubbard U values
  subroutine initHubbardUs_(input, orb, hamiltonianType, hubbU)
    type(TInputData), intent(in) :: input
    type(TOrbitals), intent(in) :: orb
    integer, intent(in) :: hamiltonianType
    real(dp), allocatable, intent(out) :: hubbU(:,:)

    integer :: nSpecies

    nSpecies = size(orb%nOrbSpecies)

    allocate(hubbU(orb%mShell, nSpecies), source=0.0_dp)

    select case(hamiltonianType)
    case default
      call error("Invalid Hamiltonian")
    case(hamiltonianTypes%dftb)
      @:ASSERT(size(input%slako%skHubbU, dim=1) >= orb%mShell)
      @:ASSERT(size(input%slako%skHubbU, dim=2) == nSpecies)
      hubbU(:,:) = input%slako%skHubbU(1:orb%mShell, :)
    case(hamiltonianTypes%xtb)
      ! handled elsewhere
    end select

    if (allocated(input%ctrl%hubbU)) then
      where (input%ctrl%hubbU > 0.0_dp)
        hubbU = input%ctrl%hubbU
      end where
    end if

  end subroutine initHubbardUs_


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

    ! Check for atom(s) occurring in multiple contacts
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
  !> Initializes the socket and receives and broadcasts initial geometry.
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

  !> Initialise a transport calculation
  subroutine initTransport_(this, env, input, electronicSolver, nSpin, tempElec, tNegf,&
      & isAContactCalc, mu, negfInt, ginfo, transpar, writeTunn, tWriteLDOS, regionLabelLDOS,&
      & errStatus)

    !> Instance
    class(TDftbPlusMain), intent(in) :: this

    !> Environment settings
    type(TEnvironment), intent(inout) :: env

    !> Holds the parsed input data.
    type(TInputData), intent(in) :: input

    !> Electronic solver for the system
    type(TElectronicSolver), intent(inout) :: electronicSolver

    !> Number of spin components, 1 is unpolarised, 2 is polarised, 4 is noncolinear / spin-orbit
    integer, intent(in) :: nSpin

    !> Electron temperature
    real(dp), intent(in) :: tempElec

    !> Transport interface
    logical, intent(in) :: tNegf

    !> Whether a contact Hamiltonian is being computed and stored
    logical, intent(out) :: isAContactCalc

    !> Holds spin-dependent electrochemical potentials of contacts
    !> This is because libNEGF is not spin-aware
    real(dp), allocatable, intent(out) :: mu(:,:)

    !> Transport interface
    type(TNegfInt), intent(out) :: negfInt

    !> Container for data needed by libNEGF
    type(TNegfInfo), intent(out) :: ginfo

    !> Transport calculation parameters
    type(TTransPar), intent(out) :: transpar

    !> Should tunnelling be written out
    logical, intent(out) :: writeTunn

    !> Should local density of states be written out
    logical, intent(out) :: tWriteLDOS

    !> Labels for different regions for DOS output
    character(lc), allocatable, intent(out) :: regionLabelLDOS(:)

    !> Operation status, if an error needs to be returned
    type(TStatus), intent(inout) :: errStatus

    logical :: tAtomsOutside
    integer :: iSpin
    integer :: nSpinChannels, iCont, jCont
    real(dp) :: mu1, mu2

    ! It is a contact calculation in the case that some contact is computed
    isAContactCalc = (input%transpar%taskContInd /= 0)

    if (nSpin <= 2) then
      nSpinChannels = nSpin
    else
      nSpinChannels = 1
    endif

    associate(transpar => input%transpar, greendens => input%ginfo%greendens)
      ! Non colinear spin not yet supported
      ! Include the built-in potential as in negf init, but the whole
      ! scc only works for
      ! calculation without spin (poisson does not support spin dependent
      ! built in potentials)
      if (transpar%ncont > 0) then
        call testTransportAtoms_(transpar, this%nAtom, errStatus)
        @:PROPAGATE_ERROR(errStatus)
        allocate(mu(transpar%ncont, nSpinChannels))
        mu(:,:) = 0.0_dp
        do iSpin = 1, nSpinChannels
          mu(1:transpar%ncont, iSpin) =&
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
        allocate(mu(1, nSpinChannels))
        mu(1, 1:nSpinChannels) = greendens%oneFermi(1:nSpinChannels)
      end if

    end associate

    if (tNegf) then
      write(stdOut,*) 'init negf'

      ! Some checks and initialization of GDFTB/NEGF
      call TNegfInt_init(negfInt, input%transpar, env, input%ginfo%greendens,&
          & input%ginfo%tundos, tempElec, this%coord0, this%cutOff%skCutoff, this%tPeriodic)

      ginfo = input%ginfo

      if (allocated(input%ctrl%indMovedAtom)) then
        tAtomsOutside = any(input%ctrl%indMovedAtom < input%transpar%idxdevice(1))&
            & .or. any(input%ctrl%indMovedAtom > input%transpar%idxdevice(2))
        if (tAtomsOutside) then
          call error("There are moving atoms specified outside of the device region")
        end if
      end if

      if (input%ctrl%tLatOpt) then
        call error("Lattice optimisation is not currently possible with transport")
      end if

    end if

    transpar = input%transpar

    !Write Dos and tunneling on separate files?
    writeTunn = ginfo%tundos%writeTunn
    tWriteLDOS = ginfo%tundos%writeLDOS

    if (tWriteLDOS) then
      call move_alloc(ginfo%tundos%dosLabels, regionLabelLDOS)
    end if

  end subroutine initTransport_


  !> Checks for atoms in two regions or not in any region of a transport calculation
  subroutine testTransportAtoms_(transpar, nAtom, errStatus)

    !> Transport calculation parameters
    type(TTransPar), intent(in) :: transpar

    !> Number of atoms in the provided device + contacts
    integer, intent(in) :: nAtom

    !> Operation status, if an error needs to be returned
    type(TStatus), intent(inout) :: errStatus

    logical, allocatable :: atomInRegion(:)
    integer :: ii

    allocate(atomInRegion(nAtom), source=.false.)

    ! check for atoms in multiple contact ranges/device
    atomInRegion(transpar%idxdevice(1):transpar%idxdevice(2)) = .true.
    do ii = 1, transpar%nCont
      if (any(atomInRegion(transpar%contacts(ii)%idxrange(1):transpar%contacts(ii)%idxrange(2))))&
          & then
        @:RAISE_FORMATTED_ERROR(errStatus, -1, "('Contact ''', A,''' contains atom ', I0,&
            & ' which is also in another contact or device region')",&
            & trim(transpar%contacts(ii)%name),&
            & findloc(atomInRegion(transpar%contacts(ii)%idxrange(1):&
            & transpar%contacts(ii)%idxrange(2)), .true.) + transpar%contacts(ii)%idxrange(1))
      end if
      atomInRegion(transpar%contacts(ii)%idxrange(1):transpar%contacts(ii)%idxrange(2)) = .true.
    end do

    ! Check for atoms missing from any of these regions
    if (any(.not.atomInRegion)) then
      @:RAISE_FORMATTED_ERROR(errStatus, -1, "('Atom ',I0, ' not in any contact or the device')",&
          & findloc(atomInRegion, .false.))
    end if

  end subroutine testTransportAtoms_

#:endif

  !> Initialises (clears) output files.
  subroutine initOutputFiles(this, env)

    !> Instance
    class(TDftbPlusMain), intent(inout) :: this

    !> Environment
    type(TEnvironment), intent(inout) :: env

    call TTaggedWriter_init(this%taggedWriter)

    if (this%tWriteAutotest) then
      call clearFile(autotestTag)
    end if
    if (this%tWriteResultsTag) then
      call clearFile(resultsTag)
    end if
    if (this%tWriteBandDat) then
      call clearFile(bandOut)
      if (this%doPerturbation .and. this%isEResp) then
        call clearFile(derivEBandOut)
      end if
    end if
    if (this%tDerivs) then
      call clearFile(hessianOut)
    end if
    if (this%tWriteDetailedOut) then
      call clearFile(userOut)
    end if
    if (this%tMD) then
      call clearFile(mdOut)
    end if
    if (this%isGeoOpt .or. this%tMD) then
      call clearFile(trim(this%geoOutFile) // ".gen")
      call clearFile(trim(this%geoOutFile) // ".xyz")
    end if
    if (allocated(this%electrostatPot)) then
      call clearFile(this%electrostatPot%espOutFile)
    end if

  end subroutine initOutputFiles


  !> Allocates most of the large arrays needed during the DFTB run.
  subroutine initArrays(this, env, input)

    !> Instance
    class(TDftbPlusMain), intent(inout) :: this

    !> Current environment
    type(TEnvironment), intent(in) :: env

    !> Input data structure
    type(TInputData), intent(in) :: input

    logical :: isREKS
    ! dLev is the number of states to include for non-adiabatic couplings
    integer :: nSpinHams, sqrHamSize, iDet, dLev

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

      if (this%isLinResp) then
        ! For CI optimization store gradient for several states,
        ! otherwise store excited state gradient for state of interest only
        if(this%linearResponse%isCIopt) then
          if (.not. this%linearResponse%tNaCoupling) then
            call error("Optimization of CI requires StateCouplings keyword.")
          end if
          dLev = this%linearResponse%indNACouplings(2) - this%linearResponse%indNACouplings(1) + 1
          if (this%linearResponse%indNACouplings(1) == 0) then
            allocate(this%excitedDerivs(3, this%nAtom, dLev-1))
          else
            allocate(this%excitedDerivs(3, this%nAtom, dLev))
          end if
          else  if (this%tLinRespZVect .and. this%tCasidaForces) then
            allocate(this%excitedDerivs(3, this%nAtom, 1))
        end if
        this%isCIopt = this%linearResponse%isCIopt
      end if
    end if

    if (this%isLinResp) then
      if(this%linearResponse%tNaCoupling) then
        dLev = this%linearResponse%indNACouplings(2) - this%linearResponse%indNACouplings(1) + 1
        allocate(this%naCouplings(3, this%nAtom, dLev*(dLev-1)/2))
      end if
    end if

    call TPotentials_init(this%potential, this%orb, this%nAtom, this%nSpin, &
        & this%nDipole, this%nQuadrupole, input%ctrl%atomicExtPotential)

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

    if (this%tPrintEigVecs .and. .not. this%electronicSolver%providesEigenvals) then
      call error("Eigenvectors are not available with this solver, so cannot be written to disc")
    end if

    if (this%electronicSolver%iSolver == electronicSolverTypes%OnlyTransport) then
      if (this%tForces) then
        call error("TransportOnly calculations cannot evaluate forces")
      end if
      if (this%tMulliken) then
        call error("TransportOnly calculations cannot evaluate charges")
      end if
      if (this%tAtomicEnergy) then
        call error("TransportOnly calculations cannot evaluate atom energies")
      end if
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
  subroutine initDetArrays(this, nLocalRows, nLocalCols)

    !> Instance
    class(TDftbPlusMain), intent(inout) :: this

    !> Size descriptors for MPI parallel execution
    integer, intent(in) :: nLocalRows, nLocalCols

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
      if (this%isHybridXc) then
        allocate(this%deltaRhoDets(nLocalRows, nLocalCols, this%nIndepSpin, this%nDets))
        this%deltaRhoDets(:,:,:,:) = 0.0_dp
      end if
    end if

  end subroutine initDetArrays


#:if WITH_TRANSPORT

  ! initialize arrays for tranpsport
  subroutine initUploadArrays_(transpar, orb, nSpin, hasBlockCharges, shiftPerLUp, chargeUp,&
      & blockUp)

    !> Transport calculation parameters
    type(TTransPar), intent(inout) :: transpar

    !> Data structure for atomic orbitals
    type(TOrbitals), intent(in) :: orb

    !> Number of spin channels
    integer, intent(in) :: nSpin

    !> Are block charges expected for the contact?
    logical, intent(in) :: hasBlockCharges

    !> Electrostatic shifts from contact(s)
    real(dp), allocatable, intent(out) :: shiftPerLUp(:,:)

    !> Uploaded charges from contact(s)
    real(dp), allocatable, intent(out) :: chargeUp(:,:,:)

    !> Block charges from contact(s), if present
    real(dp), allocatable, intent(out) :: blockUp(:,:,:,:)

    integer :: nAtom

    nAtom = size(orb%nOrbAtom)
    allocate(shiftPerLUp(orb%mShell, nAtom))
    allocate(chargeUp(orb%mOrb, nAtom, nSpin))
    if (hasBlockCharges) then
      allocate(blockUp(orb%mOrb, orb%mOrb, nAtom, nSpin))
    end if
    call readContactShifts(shiftPerLUp, chargeUp, transpar, orb, blockUp)

  end subroutine initUploadArrays_

#:endif


  !> Set up storage for dense matrices, either on a single processor, or as BLACS matrices
  subroutine allocateDenseMatrices(this, env)

    !> Instance
    class(TDftbPlusMain), intent(inout) :: this

    !> Computing environment
    type(TEnvironment), intent(in) :: env

    !! True, if hybrid xc-functional calculation requested and MPI-ready algorithm
    !! has been selected
    logical :: hybridXcAlgoNonDistributed

    integer :: nLocalCols, nLocalRows, nLocalKS

    if (this%isHybridXc) then
      hybridXcAlgoNonDistributed = .not. (this%tRealHS&
          & .and. (this%hybridXc%hybridXcAlg == hybridXcAlgo%matrixBased))
    else
      hybridXcAlgoNonDistributed = .false.
    end if

    nLocalKS = size(this%parallelKS%localKS, dim=2)
  #:if WITH_SCALAPACK
    if (hybridXcAlgoNonDistributed) then
      nLocalRows = this%denseDesc%fullSize
      nLocalCols = this%denseDesc%fullSize
    else
      call scalafx_getlocalshape(env%blacs%orbitalGrid, this%denseDesc%blacsOrbSqr, nLocalRows,&
          & nLocalCols)
    end if
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

  #:if WITH_SCALAPACK

    if (this%isSparseReorderRequired .and. (.not. hybridXcAlgoNonDistributed)) then
      call scalafx_getlocalshape(env%blacs%rowOrbitalGrid, this%denseDesc%blacsColumnSqr,&
          & nLocalRows, nLocalCols)
      if (this%t2Component .or. .not. this%tRealHS) then
        allocate(this%eigVecsCplxReordered(nLocalRows, nLocalCols, nLocalKS))
      else
        allocate(this%eigVecsRealReordered(nLocalRows, nLocalCols, nLocalKS))
      end if
    end if

  #:endif

  end subroutine allocateDenseMatrices


#:if WITH_SCALAPACK
  #!
  #! SCALAPACK related routines
  #!

  !> Initialise parallel large matrix decomposition methods
  subroutine initBlacs(blacsOpts, nAtom, nOrb, t2Component, env, errStatus)

    !> BLACS settings
    type(TBlacsOpts), intent(in) :: blacsOpts

    !> Number of atoms
    integer, intent(in) :: nAtom

    !> Number of orbitals
    integer, intent(in) :: nOrb

    !> Is this a two component calculation
    logical, intent(in) :: t2Component

    !> Computing environment data
    type(TEnvironment), intent(inout) :: env

    !> Operation status, if an error needs to be returned
    type(TStatus), intent(inout) :: errStatus

    integer :: sizeHS

    if (t2Component) then
      sizeHS = 2 * nOrb
    else
      sizeHS = nOrb
    end if
    call env%initBlacs(blacsOpts%blockSize, blacsOpts%blockSize, sizeHS, nAtom, errStatus)
    @:PROPAGATE_ERROR(errStatus)

  end subroutine initBlacs


  !> Generate descriptions of large dense matrices in BLACS decomposition
  !>
  !> Note: It must be called after getDenseDescCommon() has been called.
  !>
  subroutine getDenseDescBlacs(env, rowBlock, colBlock, denseDesc, isSparseReorderRequired)

    !> Parallel environment
    type(TEnvironment), intent(in) :: env

    !> Size of matrix row blocks
    integer, intent(in) :: rowBlock

    !> Size of matrix column blocks
    integer, intent(in) :: colBlock

    !> Descriptor of the dense matrix
    type(TDenseDescr), intent(inout) :: denseDesc

    !> Is a data distribution for sparse with dense operations needed?
    logical, intent(in) :: isSparseReorderRequired

    integer :: nn

    nn = denseDesc%fullSize
    call scalafx_getdescriptor(env%blacs%orbitalGrid, nn, nn, rowBlock, colBlock,&
        & denseDesc%blacsOrbSqr)

    if (isSparseReorderRequired) then
      ! Distribution to put entire rows on each processor
      call scalafx_getdescriptor(env%blacs%rowOrbitalGrid, nn, nn, nn, 1, denseDesc%blacsColumnSqr)
    end if

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
    select case (input%method)
    case ('mbd-rsscs')
      write(stdOut, "(A,T30,A)") "  Parameters", tmpStr
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
  subroutine ensureSolverCompatibility(iSolver, kPoints, parallelOpts, nIndepSpin, tempElec)

    !> Solver number (see dftbp_elecsolvers_elecsolvertypes)
    integer, intent(in) :: iSolver

    !> Set of k-points used in calculation (or [0,0,0] if molecular)
    real(dp), intent(in) :: kPoints(:,:)

    !> Options for a parallel calculation, if needed
    type(TParallelOpts), intent(in), allocatable :: parallelOpts

    !> Number of indepdent hamiltonian matrices at a given k-value
    integer, intent(in) :: nIndepSpin

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
      if (tElsiSolver .and. parallelOpts%nGroup /= nIndepSpin * nKPoint) then
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
    type(TWrappedInt1), intent(in) :: customOccAtoms(:)

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

    do iAt = 1, nAtom
      iSp = species(iAt)
      refOcc(:, iAt) = referenceN0(:, iSp)
    end do

    nCustomBlock = size(customOccAtoms)
    do iCustomBlock = 1, nCustomBlock
      do iCustomAtom = 1, size(customOccAtoms(iCustomBlock)%data)
        iAt =  customOccAtoms(iCustomBlock)%data(iCustomAtom)
        refOcc(:, iAt) = customOccFillings(:,iCustomBlock)
      end do
    end do

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


  !> Stop if any setting incompatible with the constrained DFTB formalism is found.
  subroutine ensureConstrainedDftbReqs(this, elecConstraintInp)

    !> Instance
    class(TDftbPlusMain), intent(inout) :: this

    !> Input parameters for electronic constraints
    type(TElecConstraintInput), intent(in) :: elecConstraintInp

    if (.not. this%tSccCalc) then
      call error("Electronically constrained calculations do not yet support non-SCC calculations.")
    end if

    if (this%isXlbomd) then
      call error("Electronically constrained calculations do not yet support XLBOMD.")
    end if

    if (allocated(this%reks)) then
      call error("Electronically constrained calculations do not yet support REKS.")
    end if

    if (this%deltaDftb%isNonAufbau) then
      call error("Electronically constrained calculations do not yet support delta-DFTB.")
    end if

    if (this%isElecDyn) then
      call error("Electronically constrained calculations do not yet support electron dynamics.")
    end if

  end subroutine ensureConstrainedDftbReqs


  !> Stop if any setting incompatible with using hybrid xc-functionals are found.
  subroutine ensureHybridXcReqs(this, tShellResolved, hybridXcInp)

    !> Instance
    type(TDftbPlusMain), intent(inout) :: this

    !> True, if this is a shell resolved calculation
    logical, intent(in) :: tShellResolved

    !> Parameters for the range separated calculation
    type(THybridXcInp), intent(in) :: hybridXcInp

    if (withMpi .and. (.not. this%tPeriodic) .and. (hybridXcInp%hybridXcAlg&
        & /= hybridXcAlgo%matrixBased)) then
      call error("MPI parallelized hybrid calculations of non-periodic systems only available for&
          & the matrix-based algorithm.")
    end if

    if (withMpi .and. this%tRealHS .and. this%tReadChrg) then
      call error("Hybrid functionality for non-periodic or Gamma-point only calculations does not&
          & yet support restarts for MPI-enabled builds.")
    end if

    if (this%tPeriodic) then
      if ((.not. this%tRealHS) .and. (hybridXcInp%hybridXcAlg /= hybridXcAlgo%neighbourBased)) then
        call error("Hybrid functionality for periodic systems beyond the Gamma-point currently only&
            & working for the neighbour list based algorithm.")
      end if
      if (this%tRealHS .and. hybridXcInp%hybridXcAlg == hybridXcAlgo%thresholdBased) then
        call error("Hybrid functionality at the Gamma-point not implemented for the threshold&
            & algorithm.")
      end if
    end if

    if ((.not. this%tRealHS) .and. this%tForces&
        & .and. (hybridXcInp%gammaType /= hybridXcGammaTypes%truncated)) then
      call error("Hybrid functionals don't yet support gradient calculations for periodic systems&
          & beyond the Gamma point for CoulombMatrix settings other than 'Truncated'.")
    end if

    if (this%tPeriodic .and. this%tRealHS&
        & .and. hybridXcInp%gammaType /= hybridXcGammaTypes%truncated&
        & .and. hybridXcInp%gammaType /= hybridXcGammaTypes%truncatedAndDamped) then
      call error("Hybrid functionals don't yet support Gamma point calculations for CoulombMatrix&
          & types other than 'Truncated' or 'Truncated+Damping'.")
    end if

    if (this%tHelical) then
      call error("Hybrid functionality only works with non-helical structures at the moment.")
    end if

    if (this%tReadChrg .and. hybridXcInp%hybridXcAlg == hybridXcAlgo%thresholdBased) then
      call error("Restart on thresholded range separation not working correctly.")
    end if

    if (tShellResolved) then
      call error("Hybrid functionality currently does not yet support shell-resolved SCF.")
    end if

    if (this%tAtomicEnergy) then
      call error("Atomic resolved energies cannot be calculated with hybrid functionals at the&
          & moment.")
    end if

    if (this%nSpin > 2) then
      call error("Hybrid calculations not implemented for non-colinear calculations.")
    end if

    if ((.not. this%tRealHS) .and. this%nSpin == 2&
        & .and. hybridXcInp%gammaType /= hybridXcGammaTypes%truncated) then
      call error("Hybrid functionality does not yet support spin-polarized calculations of periodic&
          & systems beyond the Gamma-point for CoulombMatrix settings other than 'Truncated'.")
    end if

    if ((.not. this%tRealHS) .and. this%nSpin == 2 .and. this%tForces) then
      call error("Hybrid functionality currently does not yet support spin-polarized gradient&
          & evaluation for periodic systems beyond the Gamma-point.")
    end if

    if (this%tSpinOrbit) then
      call error("Hybrid calculations not currently implemented for spin orbit coupling.")
    end if

    if (this%isXlbomd) then
      call error("Hybrid calculations not currently implemented for XLBOMD.")
    end if

    if (this%t3rd) then
      call error("Hybrid calculations not currently implemented for 3rd-order DFTB.")
    end if

    if (this%isRS_LinResp .and. hybridXcInp%hybridXcType == hybridXcFunc%cam) then
      call error("General CAM functionals not currently implemented for linear response.")
    end if

    if (this%isRS_LinResp .and. hybridXcInp%hybridXcType == hybridXcFunc%hyb) then
      call error("Global hybrid functionals not currently implemented for linear response.")
    end if

  end subroutine ensureHybridXcReqs


  !> Stop if linear response module can not be invoked due to unimplemented combinations of
  !! features.
  subroutine ensureLinRespConditions(tSccCalc, t3rd, tRealHS, tPeriodic, tCasidaForces, solvation,&
      & isRS_LinResp, nSpin, tHelical, tSpinOrbit, isDftbU, tempElec, isOnsiteCorrected, input)

    !> Is the calculation SCC?
    logical, intent(in) :: tSccCalc

    !> 3rd order hamiltonian contributions included
    logical, intent(in) :: t3rd

    !> A real hamiltonian
    logical, intent(in) :: tRealHs

    !> Periodic boundary conditions
    logical, intent(in) :: tPeriodic

    !> Forces being evaluated in the excited state
    logical, intent(in) :: tCasidaForces

    !> Solvation data and calculations
    class(TSolvation), allocatable :: solvation

    !> Is this an excited state calculation with range separation
    logical, intent(in) :: isRS_LinResp

    !> Number of spin components, 1 is unpolarised, 2 is polarised, 4 is noncolinear / spin-orbit
    integer, intent(in) :: nSpin

    !> If the calculation is helical geometry
    logical :: tHelical

    !> Is there spin-orbit coupling
    logical, intent(in) :: tSpinOrbit

    !> Is this a DFTB+U calculation?
    logical, intent(in) :: isDftbU

    !> Temperature of the electrons
    real(dp), intent(in) :: tempElec

    !> Is this a onsite corrected TD-DFTB calculation?
    logical, intent(in) :: isOnsiteCorrected

    !> Holds the parsed input data.
    type(TInputData), intent(in), target :: input

    character(lc) :: tmpStr

    @:ASSERT(allocated(input%ctrl%lrespini))

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
      call error("Linear response does not work with non-colinear spin polarization yet")
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

    if (input%ctrl%lrespini%nstat == 0 .and. (.not. input%ctrl%lrespini%isCIopt)) then
      if (tCasidaForces) then
        call error("Excited forces only available for StateOfInterest non zero.")
      end if
      if (input%ctrl%lrespini%tPrintEigVecs .or. input%ctrl%lrespini%tCoeffs) then
        call error("Excited eigenvectors only available for StateOfInterest non zero.")
      end if
    end if

    if (isOnsiteCorrected .and. input%ctrl%lrespini%iLinRespSolver == linRespSolverTypes%Stratmann)&
        & then
      call error("Onsite corrections not implemented for Stratmann diagonaliser.")
    end if

    if (isRS_LinResp) then
      if (input%ctrl%lrespini%iLinRespSolver /= linRespSolverTypes%Stratmann) then
        call error("TD-LC-DFTB implemented only for Stratmann diagonaliser.")
      end if
      if (tPeriodic) then
        call error("Range separated excited states for periodic geometries are currently&
            & unavailable")
      end if
      if (input%ctrl%lrespini%tEnergyWindow .or. input%ctrl%lrespini%tOscillatorWindow) then
        call error("Range separated excited states not available for window options.")
      end if
      if (input%ctrl%lrespini%sym == 'B' .or. input%ctrl%lrespini%sym == 'T') then
        call warning("Range separated excited states not well tested for triplet excited states!")
      end if
      if (input%ctrl%tSpin) then
        call warning("Range separated excited states not well tested for spin-polarized systems!")
      end if
    else
      if (input%ctrl%lrespini%energyWindow < 0.0_dp) then
        call error("Negative energy window for excitations")
      end if
    end if


  end subroutine ensureLinRespConditions


  !> Determine range separated cut-off and also update maximal cutoff
  subroutine getHybridXcCutOff_cluster(cutOff, cutoffRed)

    !> Resulting cutoff
    type(TCutoffs), intent(inout) :: cutOff

    !> CAM-neighbour list cutoff reduction
    real(dp), intent(in) :: cutoffRed

    if (cutoffRed < 0.0_dp) then
      call error("Cutoff reduction for range-separated neighbours should be zero or positive.")
    end if

    cutOff%camCutOff = cutOff%skCutOff - cutoffRed

    if (cutOff%camCutOff < 0.0_dp) then
      call error("Screening cutoff for range-separated neighbours too short.")
    end if

    cutOff%mCutoff = max(cutOff%mCutOff, cutoff%camCutOff)

    ! dummy cutoff values
    cutOff%gammaCutoff = 1.0_dp
    cutOff%gSummationCutoff = 1.0_dp

  end subroutine getHybridXcCutOff_cluster


  !> Determine range separated cut-off and also update maximal cutoff
  subroutine getHybridXcCutOff_gamma(cutOff, latVecs, cutoffRed, gSummationCutoff,&
      & gammaCutoff)

    !> Resulting cutoff
    type(TCutoffs), intent(inout) :: cutOff

    !> Lattice vectors
    real(dp), intent(in) :: latVecs(:,:)

    !> CAM-neighbour list cutoff reduction
    real(dp), intent(in) :: cutoffRed

    !> Cutoff for real-space g-summation
    real(dp), intent(in), optional :: gSummationCutoff

    !> Coulomb truncation cutoff for Gamma electrostatics
    real(dp), intent(in), optional :: gammaCutoff

    !! Error status
    type(TStatus) :: errStatus

    if (cutoffRed < 0.0_dp) then
      call error("Cutoff reduction for range-separated neighbours should be zero or positive.")
    end if

    cutOff%camCutOff = cutOff%skCutOff - cutoffRed

    if (cutOff%camCutOff < 0.0_dp) then
      call error("Screening cutoff for range-separated neighbours too short.")
    end if

    cutOff%mCutoff = max(cutOff%mCutOff, cutoff%camCutOff)

    if (present(gammaCutoff)) then
      cutOff%gammaCutoff = gammaCutoff
    else
      allocate(cutOff%gammaCutoff)
      call getCoulombTruncationCutoff(latVecs, cutOff%gammaCutoff, errStatus)
      @:PROPAGATE_ERROR(errStatus)
      if (errStatus%hasError()) call error(errStatus%message)
    end if

    if (present(gSummationCutoff)) then
      cutOff%gSummationCutoff = gSummationCutoff
    else
      ! This would correspond to "the safest option"
      cutOff%gSummationCutoff = 2.0_dp * cutOff%gammaCutoff
    end if

  end subroutine getHybridXcCutOff_gamma


  !> Determine range separated cut-off and also update maximal cutoff
  subroutine getHybridXcCutOff_kpts(cutOff, latVecs, cutoffRed, supercellFoldingDiag,&
      & gSummationCutoff, wignerSeitzReduction, gammaCutoff)

    !> Resulting cutoff
    type(TCutoffs), intent(inout) :: cutOff

    !> Lattice vectors
    real(dp), intent(in) :: latVecs(:,:)

    !> CAM-neighbour list cutoff reduction
    real(dp), intent(in) :: cutoffRed

    !> Supercell folding coefficients and shifts
    integer, intent(in) :: supercellFoldingDiag(:)

    !> Cutoff for real-space g-summation
    real(dp), intent(in), optional :: gSummationCutoff

    !> Number of unit cells along each supercell folding direction to substract from MIC
    !! Wigner-Seitz cell construction
    integer, intent(in), optional :: wignerSeitzReduction

    !> Coulomb truncation cutoff for Gamma electrostatics
    real(dp), intent(in), optional :: gammaCutoff

    !! Error status
    type(TStatus) :: errStatus

    if (cutoffRed < 0.0_dp) then
      call error("Cutoff reduction for range-separated neighbours should be zero or positive.")
    end if

    cutOff%camCutOff = cutOff%skCutOff - cutoffRed

    if (cutOff%camCutOff < 0.0_dp) then
      call error("Screening cutoff for range-separated neighbours too short.")
    end if

    cutOff%mCutoff = max(cutOff%mCutOff, cutoff%camCutOff)

    if (present(gammaCutoff)) then
      cutOff%gammaCutoff = gammaCutoff
    else
      allocate(cutOff%gammaCutoff)
      call getCoulombTruncationCutoff(latVecs, cutOff%gammaCutoff, errStatus,&
          & nK=supercellFoldingDiag)
      @:PROPAGATE_ERROR(errStatus)
      if (errStatus%hasError()) call error(errStatus%message)
    end if

    if (present(gSummationCutoff)) then
      cutOff%gSummationCutoff = gSummationCutoff
    else
      ! This would correspond to "the safest option"
      cutOff%gSummationCutoff = 2.0_dp * cutOff%gammaCutoff
    end if

    if (present(wignerSeitzReduction)) then
      cutOff%wignerSeitzReduction = wignerSeitzReduction
    else
      cutOff%wignerSeitzReduction = 1
    end if

  end subroutine getHybridXcCutOff_kpts


  !> Pre-allocate density matrix for range-separation.
  subroutine reallocateHybridXc(this, nLocalRows, nLocalCols, nLocalKS)

    !> Instance
    type(TDftbPlusMain), intent(inout) :: this

    !> Size descriptors for MPI parallel execution
    integer, intent(in) :: nLocalRows, nLocalCols, nLocalKS

    !! Global iKS composite index
    integer :: iGlobalKS

    !! Spin and k-point indices
    integer :: iSpin, iK

    ! deallocate arrays, if already allocated
    if (allocated(this%densityMatrix%deltaRhoOut)) deallocate(this%densityMatrix%deltaRhoOut)
    if (allocated(this%densityMatrix%deltaRhoOutCplx))&
        & deallocate(this%densityMatrix%deltaRhoOutCplx)
    if (allocated(this%densityMatrix%iKiSToiGlobalKS))&
        & deallocate(this%densityMatrix%iKiSToiGlobalKS)
    if (allocated(this%densityMatrix%deltaRhoOutCplxHS))&
        & deallocate(this%densityMatrix%deltaRhoOutCplxHS)

    if (this%tRealHS) then
      ! Prevent for deleting charges read in from file
      if (.not. allocated(this%densityMatrix%deltaRhoIn)) then
        allocate(this%densityMatrix%deltaRhoIn(nLocalRows, nLocalCols, nLocalKS))
        this%densityMatrix%deltaRhoIn(:,:,:) = 0.0_dp
      end if
      allocate(this%densityMatrix%deltaRhoOut(nLocalRows, nLocalCols, nLocalKS))
      this%densityMatrix%deltaRhoOut(:,:,:) = 0.0_dp
    elseif (this%tReadChrg .and. (.not. allocated(this%supercellFoldingDiag))) then
      ! in case of k-points and restart from file, we have to wait until charges.bin was read
      if (.not. allocated(this%densityMatrix%deltaRhoInCplxHS)) then
        allocate(this%densityMatrix%deltaRhoInCplxHS(0, 0, 0, 0, 0, 0))
        this%densityMatrix%deltaRhoInCplxHS(:,:,:,:,:,:) = 0.0_dp
      end if
      allocate(this%densityMatrix%deltaRhoOutCplxHS(0, 0, 0, 0, 0, 0))
      this%densityMatrix%deltaRhoOutCplxHS(:,:,:,:,:,:) = 0.0_dp
    else
      ! normal k-point case, without restart from file
      if (.not. allocated(this%densityMatrix%deltaRhoInCplxHS)) then
        allocate(this%densityMatrix%deltaRhoInCplxHS(this%nOrb, this%nOrb,&
            & this%supercellFoldingDiag(1), this%supercellFoldingDiag(2),&
            & this%supercellFoldingDiag(3), this%nIndepSpin))
        this%densityMatrix%deltaRhoInCplxHS(:,:,:,:,:,:) = 0.0_dp
      end if
        allocate(this%densityMatrix%deltaRhoOutCplxHS(this%nOrb, this%nOrb,&
            & this%supercellFoldingDiag(1), this%supercellFoldingDiag(2),&
            & this%supercellFoldingDiag(3), this%nIndepSpin))
        this%densityMatrix%deltaRhoOutCplxHS(:,:,:,:,:,:) = 0.0_dp
    end if

    if (.not. this%tRealHS) then

      allocate(this%densityMatrix%deltaRhoOutCplx(this%nOrb, this%nOrb,&
          & size(this%parallelKS%localKS, dim=2)))
      this%densityMatrix%deltaRhoOutCplx(:,:,:) = 0.0_dp

      ! Build spin/k-point composite index for all spins and k-points (global)
      iGlobalKS = 1
      allocate(this%densityMatrix%iKiSToiGlobalKS(this%nKPoint, this%nSpin))
      do iSpin = 1, this%nSpin
        do iK = 1, this%nKPoint
          this%densityMatrix%iKiSToiGlobalKS(iK, iSpin) = iGlobalKS
          iGlobalKS = iGlobalKS + 1
        end do
      end do
    end if

  end subroutine reallocateHybridXc


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
      & tSccCalc, tSpin, tSpinOrbit, isDftbU, isExtField, isLinResp, tPeriodic, tLatOpt, tReadChrg,&
      & tPoisson, isShellResolved)

    !> Data type for REKS input
    type(TReksInp), intent(in) :: reksInp

    !> Solvation data and calculations
    class(TSolvation), allocatable, intent(in) :: solvation

    !> Correction to energy from on-site matrix elements
    real(dp), allocatable, intent(in) :: onSiteElements(:,:,:,:)

    !> The k-points
    real(dp), intent(in) :: kPoint(:,:)

    !> Nr. of electrons
    real(dp), intent(in) :: nEl(:)

    !> Nr. of K-points
    integer, intent(in) :: nKPoint

    !> Is the calculation SCC?
    logical, intent(in) :: tSccCalc

    !> Is this a spin polarized calculation?
    logical, intent(in) :: tSpin

    !> Is there spin-orbit coupling
    logical, intent(in) :: tSpinOrbit

    !> Is this a DFTB+U calculation?
    logical, intent(in) :: isDftbU

    !> External electric field
    logical, intent(in) :: isExtField

    !> Calculate Casida linear response excitations
    logical, intent(in) :: isLinResp

    !> If calculation is periodic
    logical, intent(in) :: tPeriodic

    !> Optimise lattice constants?
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
    else if (isExtField) then
      call error("REKS is not compatible with external field and potentials, only point charge&
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
      call error("Lattice optimisation is only possible&
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
      & hamiltonianType, nSpin, nExtChrg, is3rd, isHybridXc, isDispersion, isQNetAllocated,&
      & tForces, tPeriodic, tStress, tDipole)

    !> Data type for REKS
    type(TReksCalc), intent(out) :: reks

    !> Data type for REKS input
    type(TReksInp), intent(inout) :: reksInp

    !> Electronic solver for the system
    type(TElectronicSolver), intent(in) :: electronicSolver

    !> Atomic orbital information
    type(TOrbitals), intent(in) :: orb

    !> Spin W values
    real(dp), intent(inout) :: spinW(:,:,:)

    !> Nr. of electrons
    real(dp), intent(in) :: nEl(:)

    !> Coordinates and charges of external point charges
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
    logical, intent(in) :: isHybridXc

    !> Whether to run a dispersion calculation
    logical, intent(in) :: isDispersion

    !> Logical to determine whether to calculate net charge per atom (qNetAtom)
    logical, intent(in) :: isQNetAllocated

    !> Do we need forces?
    logical, intent(in) :: tForces

    !> If calculation is periodic
    logical, intent(in) :: tPeriodic

    !> Can stress be calculated?
    logical, intent(in) :: tStress

    !> Calculate an electric dipole?
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
            & blurWidths, is3rd, isHybridXc, isDispersion, isQNetAllocated, tForces,&
            & tPeriodic, tStress, tDipole)
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

    !> Data type for REKS
    type(TReksCalc), intent(in) :: reks

    !> Data structure with atomic orbital information
    type(TOrbitals), intent(in) :: orb

    !> Labels of atomic species
    character(mc), intent(in) :: speciesName(:)

    !> Nr of different types (nAtom)
    integer, intent(in) :: nType

    integer :: ii, iType
    character(lc) :: strTmp

    write(stdOut,*)
    write(stdOut,*)
    write(stdOut, "(A,':',T30,A)") "REKS Calculation", "Yes"

    select case (reks%reksAlg)
    case (reksTypes%noReks)
    case (reksTypes%ssr22)
      write(stdOut, "(A,':',T30,A)") "SSR(2,2) Calculation", "Yes"
      if (reks%Efunction == 1) then
        write(stdOut, "(A,':',T30,A)") "Energy Functional", "PPS"
      else if (reks%Efunction == 2) then
        write(stdOut, "(A,':',T30,A)") "Energy Functional", "(PPS+OSS)/2"
      end if
    case (reksTypes%ssr44)
      call error("SSR(4,4) is not implemented yet")
    end select

    write(stdOut, "(A,':',T30,I14)") "Number of Core Orbitals", reks%Nc
    write(stdOut, "(A,':',T30,I14)") "Number of Active Orbitals", reks%Na
    write(stdOut, "(A,':',T30,I14)") "Number of Basis", orb%nOrb
    write(stdOut, "(A,':',T30,I14)") "Number of States", reks%nstates
    do ii = 1, reks%SAstates
      if (ii == 1) then
        write(strTmp, "(A,':')") "State-Averaging Weight"
      else
        write(strTmp, "(A)") ""
      end if
      write(stdOut, "(A,T30,F12.6)") trim(strTmp), reks%SAweight(ii)
    end do
    write(stdOut, "(A,':',T30,I14)") "State of Interest", reks%rstate

    if (reks%tReadMO) then
      write(stdOut, "(A,':',T30,A)") "Initial Guess", "Read Eigenvec.bin file"
    else
      write(stdOut, "(A,':',T30,A)") "Initial Guess", "Diagonalise H0 matrix"
    end if

    write(stdOut, "(A,':',T30,A)") "Newton-Raphson for FON opt", "Yes"
    write(stdOut, "(A,':',T30,I14)") "NR max. Iterations", reks%FonMaxIter
    if (reks%shift > epsilon(1.0_dp)) then
      write(stdOut, "(A,':',T30,A)") "Level Shifting", "Yes"
    else
      write(stdOut, "(A,':',T30,A)") "Level Shifting", "No"
    end if
    write(stdOut, "(A,':',T30,F12.6)") "Shift Value", reks%shift

    do iType = 1, nType
      if (iType == 1) then
        write(strTmp, "(A,':')") "W Scale Factor"
      else
        write(strTmp, "(A)") ""
      end if
      write(stdOut, "(A,T30,A3,'=',F12.6)") trim(strTmp), speciesName(iType), reks%Tuning(iType)
    end do

    if (reks%tTDP) then
      write(stdOut, "(A,':',T30,A)") "Transition Dipole", "Yes"
    else
      write(stdOut, "(A,':',T30,A)") "Transition Dipole", "No"
    end if

    if (reks%tForces) then

      if (reks%Lstate > 0) then
        write(stdOut, "(A,':',T30,A)") "Gradient of Microstate", "Yes"
        write(stdOut, "(A,':',T30,I14)") "Index of Interest", reks%Lstate
      else
        write(stdOut, "(A,':',T30,A)") "Gradient of Microstate", "No"
      end if

      if (reks%Efunction /= 1) then
        if (reks%Glevel == 1) then
          write(stdOut, "(A,':',T30,A)") "CP-REKS Solver", "Preconditioned Conjugate-Gradient"
          write(stdOut, "(A,':',T30,I14)") "CG max. Iterations", reks%CGmaxIter
          write(stdOut, "(A,':',T30,E14.6)") "CG Tolerance", reks%Glimit
          if (reks%tSaveMem) then
            write(stdOut, "(A,':',T30,A)") "Memory for A and Hxc", "Save in Cache Memory"
          else
            write(stdOut, "(A,':',T30,A)") "Memory for A and Hxc", "Direct Updating Without Saving"
          end if
        elseif (reks%Glevel == 2) then
          write(stdOut, "(A,':',T30,A)") "CP-REKS Solver", "Conjugate-Gradient"
          write(stdOut, "(A,':',T30,I14)") "CG max. Iterations", reks%CGmaxIter
          write(stdOut, "(A,':',T30,E14.6)") "CG Tolerance", reks%Glimit
          if (reks%tSaveMem) then
            write(stdOut, "(A,':',T30,A)") "Memory for A and Hxc", "Save in Cache Memory"
          else
            write(stdOut, "(A,':',T30,A)") "Memory for A and Hxc", "Direct Updating Without Saving"
          end if
        elseif (reks%Glevel == 3) then
          write(stdOut, "(A,':',T30,A)") "CP-REKS Solver", "Direct Matrix Multiplication"
        end if
        if (reks%tNAC) then
          write(stdOut, "(A,':',T30,A)") "Non-Adiabatic Coupling", "Yes"
        end if
      end if

      if (reks%tRD) then
        write(stdOut, "(A,':',T30,A)") "Relaxed Density for QM/MM", "Yes"
      end if

    end if

  end subroutine printReksInitInfo


  ! Initializes the variables directly related to the user specified geometry.
  subroutine initGeometry_(input, nAtom, nType, tPeriodic, tHelical, boundaryCond, coord0,&
      & species0, tCoordsChanged, tLatticeChanged, latVec, origin, recVec, invLatVec, cellVol,&
      & recCellVol, errStatus)
    type(TInputData), intent(in) :: input
    integer, intent(out) :: nAtom, nType
    logical, intent(out) :: tPeriodic, tHelical
    type(TBoundaryConditions), intent(out) :: boundaryCond
    real(dp), allocatable, intent(out) :: coord0(:,:)
    integer, allocatable, intent(out) :: species0(:)
    logical, intent(out) :: tCoordsChanged, tLatticeChanged
    real(dp), allocatable, intent(out) :: latVec(:,:), origin(:), recVec(:,:), invLatVec(:,:)
    real(dp), intent(out) :: cellVol, recCellVol

    !> Operation status, if an error needs to be returned
    type(TStatus), intent(inout) :: errStatus

    nAtom = input%geom%nAtom
    nType = input%geom%nSpecies

    tPeriodic = input%geom%tPeriodic
    tHelical = input%geom%tHelical


    if (tPeriodic) then
      call TBoundaryConditions_init(boundaryCond, boundaryConditions%pbc3d, errStatus)
    else if (tHelical) then
      call TBoundaryConditions_init(boundaryCond, boundaryConditions%helical, errStatus)
    else
      call TBoundaryConditions_init(boundaryCond, boundaryConditions%cluster, errStatus)
    end if
    @:PROPAGATE_ERROR(errStatus)

    coord0 = input%geom%coords
    species0 = input%geom%species
    tCoordsChanged = .true.

    cellVol = 0.0_dp
    recCellVol = 0.0_dp
    if (tPeriodic) then
      tLatticeChanged = .true.
      latVec = input%geom%latVecs
      origin = input%geom%origin
      allocate(recVec(3, 3))
      allocate(invLatVec(3, 3))
      invLatVec(:,:) = latVec
      call matinv(invLatVec)
      invLatVec = reshape(invLatVec, [3, 3], order=[2, 1])
      recVec = 2.0_dp * pi * invLatVec
      cellVol = abs(determinant33(latVec))
      recCellVol = abs(determinant33(recVec))
    elseif (tHelical) then
      origin = input%geom%origin
      latVec = input%geom%latVecs
      allocate(recVec(1, 1))
      recVec = 1.0_dp / latVec(1,1)
      allocate(invLatVec(0, 0))
    else
      allocate(latVec(0, 0))
      allocate(origin(0))
      allocate(recVec(0, 0))
      allocate(invLatVec(0, 0))
      cellVol = 0.0_dp
      recCellVol = 0.0_dp
      tLatticeChanged = .false.
    end if

  end subroutine initGeometry_


  ! Initializes short gamma damping
  subroutine initShortGammaDamping_(ctrl, speciesMass, damping)
    type(TControl), intent(in) :: ctrl
    real(dp), intent(in) :: speciesMass(:)
    type(TShortGammaDamp), intent(out) :: damping

    integer :: nSpecies

    nSpecies = size(speciesMass)

    if (ctrl%tDampH) then
      ! Damping species: only H or any isotopes of it, but not He or anything heavier
      damping%isDamped = (speciesMass < 3.5_dp * amu__au)
      damping%exponent = ctrl%dampExp
    else
      allocate(damping%isDamped(nSpecies))
      damping%isDamped(:) = .false.
      damping%exponent = 0.0_dp
    end if

  end subroutine initShortGammaDamping_


  ! Initializes short gamma calculator
  subroutine initShortGammaInput_(ctrl, speciesMass, uniqHubbU, shortGammaDamp,&
      & shortGammaInp)
    type(TControl), intent(in) :: ctrl
    real(dp), intent(in) :: speciesMass(:)
    type(TUniqueHubbard), intent(in) :: uniqHubbU
    type(TShortGammaDamp), intent(in) :: shortGammaDamp
    type(TShortGammaInput), allocatable, intent(out) :: shortGammaInp

    allocate(shortGammaInp)
    shortGammaInp%damping = shortGammaDamp
    if (allocated(ctrl%h5Input)) then
      if (.not. any(speciesMass < 3.5_dp * amu__au)) then
        call error("H5 correction used without H atoms present")
      end if
      if (any(shortGammaInp%damping%isDamped)) then
        call error("H5 correction is not compatible with X-H damping")
      end if
      shortGammaInp%h5CorrectionInp = ctrl%h5Input
    end if
    shortGammaInp%hubbU = uniqHubbU

  end subroutine initShortGammaInput_


  ! Initializes a Coulomb-calculator
  subroutine initCoulombInput_(env, ewaldAlpha, tolEwald, boundaryCond, coulombInput)
    type(TEnvironment), intent(in) :: env
    real(dp), intent(in) :: ewaldAlpha
    real(dp), intent(in) :: tolEwald
    integer, intent(in) :: boundaryCond
    type(TCoulombInput), allocatable, intent(out) :: coulombInput

    allocate(coulombInput)
    coulombInput%ewaldAlpha = ewaldAlpha
    coulombInput%tolEwald = tolEwald
    coulombInput%boundaryCond = boundaryCond

  end subroutine initCoulombInput_


#:if WITH_POISSON

  ! Initializes a Poisson solver
  subroutine initPoissonInput_(input, nAtom, nType, species0, coord0, tPeriodic, latVec, orb,&
      & hubbU, poissonInput, shiftPerLUp)
    type(TInputData), intent(in) :: input
    integer, intent(in) :: nAtom, nType
    integer, target, intent(in) :: species0(:)
    real(dp), target, intent(in) :: coord0(:,:)
    logical, intent(in) :: tPeriodic
    real(dp), intent(in) :: latVec(:,:)
    type(TOrbitals), intent(in) :: orb
    real(dp), intent(in) :: hubbU(:,:)
    type(TPoissonInput), allocatable, intent(out) :: poissonInput
    real(dp), intent(in), optional :: shiftPerLUp(:,:)

    allocate(poissonInput)
    associate(poissStr => poissonInput%poissonStruct)
      poissStr%nAtom = nAtom
      poissStr%nSpecies = nType
      poissStr%specie0 = species0
      poissStr%x0 = coord0
      poissStr%isPeriodic = tPeriodic
      if (tPeriodic) then
        poissStr%latVecs = latVec
      else
        allocate(poissStr%latVecs(3, 3))
        poissStr%latVecs(:,:) = 0.0_dp
      end if
    end associate

    poissonInput%poissonInfo = input%poisson
    poissonInput%hubbU = hubbU
    if (present(shiftPerLUp)) then
      poissonInput%shellPotUpload = shiftPerLUp
    end if

  #:if WITH_TRANSPORT

    poissonInput%transPar = input%transPar

  #:endif

  end subroutine initPoissonInput_

#:endif


  ! Initializes the scc calculator
  subroutine initSccCalculator_(env, orb, ctrl, boundaryCond, coulombInput, shortGammaInput,&
      & poissonInput, sccCalc)
    type(TEnvironment), intent(inout) :: env
    type(TOrbitals), intent(in) :: orb
    type(TControl), intent(in) :: ctrl
    integer, intent(in) :: boundaryCond
    type(TCoulombInput), allocatable, intent(inout) :: coulombInput
    type(TShortGammaInput), allocatable, intent(inout) :: shortGammaInput
    type(TPoissonInput), allocatable, intent(inout) :: poissonInput
    type(TScc), allocatable, intent(out) :: sccCalc

    type(TSccInput) :: sccInput

    call move_alloc(coulombInput, sccInput%coulombInput)
    call move_alloc(shortGammaInput, sccInput%shortGammaInput)
    call move_alloc(poissonInput, sccInput%poissonInput)

    sccInput%boundaryCond = boundaryCond
    if (boundaryCond == boundaryConditions%helical) then
      call error("Scc calculations not currently supported for helical boundary conditions")
    end if

    if (ctrl%nExtChrg > 0) then
      sccInput%extCharges = ctrl%extChrg
      if (allocated(ctrl%extChrgBlurWidth)) then
        sccInput%blurWidths = ctrl%extChrgblurWidth
        if (any(sccInput%blurWidths < 0.0_dp)) then
          call error("Gaussian blur widths for charges may not be negative")
        end if
      end if
    end if

    if (allocated(ctrl%chrgConstr)) then
      if (any(abs(ctrl%chrgConstr(:,2)) > epsilon(1.0_dp))) then
        sccInput%chrgConstraints = ctrl%chrgConstr
      end if
    end if

    if (allocated(ctrl%thirdOrderOn)) then
      sccInput%thirdOrderOn = ctrl%thirdOrderOn
    end if

    allocate(sccCalc)
    call TScc_init(sccCalc, sccInput, env, orb)

  end subroutine initSccCalculator_


  subroutine initRepulsive_(nAtom, isPeriodic, isHelical, pairRepulsives, chimesInp, speciesNames,&
        & species0, repulsive)
    integer, intent(in) :: nAtom
    logical, intent(in) :: isPeriodic, isHelical
    type(TPairRepulsiveItem), allocatable, intent(inout) :: pairRepulsives(:,:)
    type(TChimesRepInp), allocatable, intent(in) :: chimesInp
    character(*), intent(in) :: speciesNames(:)
    integer, intent(in) :: species0(:)
    class(TRepulsive), allocatable, intent(out) :: repulsive

    type(TRepulsiveList), allocatable :: repulsiveList
    type(TTwoBodyRepInp) :: twoBodyInp

    allocate(repulsiveList)

    twoBodyInp%nAtom = nAtom
    twoBodyInp%isHelical = isHelical
    call move_alloc(pairRepulsives, twoBodyInp%pairRepulsives)
    @:CREATE_CLASS(repulsive, TTwoBodyRep, TTwoBodyRep_init, twoBodyInp)
    call repulsiveList%push(repulsive)

    #:if WITH_CHIMES
      if (allocated(chimesInp)) then
        if (.not. isPeriodic) then
          call error("ChIMES repulsives currently require periodic boundary conditions")
        end if
        if (isHelical) then
          call error("ChIMES repulsive is not compatible with helical boundary conditions")
        end if
        @:CREATE_CLASS(repulsive, TChimesRep, TChimesRep_init, chimesInp, speciesNames, species0)
        call repulsiveList%push(repulsive)
      end if
    #:endif

    ! If multiple repulsives, wrap via container, otherwise use the one directly
    if (repulsiveList%size() > 1) then
        @:CREATE_CLASS(repulsive, TRepulsiveCont, TRepulsiveCont_init, repulsiveList)
    else
      call repulsiveList%pop(repulsive)
    end if

  end subroutine initRepulsive_


  ! Decides how many Cholesky-decompositions should be buffered
  subroutine getBufferedCholesky_(tRealHS, nLocalKS, nBufferedCholesky)
    logical, intent(in) :: tRealHS
    integer, intent(in) :: nLocalKS
    integer, intent(out) :: nBufferedCholesky

    if (tRealHS) then
      nBufferedCholesky = 1
    else
      nBufferedCholesky = nLocalKS
    end if

  end subroutine getBufferedCholesky_


  #:if WITH_TRANSPORT

  !> Replace charges with those from the stored contact values
  subroutine overrideContactCharges(qOrb, qOrbUp, transpar, qBlock, qBlockUp)

    !> Input charges
    real(dp), intent(inout) :: qOrb(:,:,:)

    !> Uploaded charges
    real(dp), intent(in) :: qOrbUp(:,:,:)

    !> Transport parameters
    type(TTransPar), intent(in) :: transpar

    !> Block charges, for example from DFTB+U
    real(dp), allocatable, intent(inout) :: qBlock(:,:,:,:)

    !> Uploaded block charges
    real(dp), allocatable, intent(in) :: qBlockUp(:,:,:,:)

    integer :: ii, iStart, iEnd

    do ii = 1, transpar%ncont
      iStart = transpar%contacts(ii)%idxrange(1)
      iEnd = transpar%contacts(ii)%idxrange(2)
      qOrb(:,iStart:iEnd,:) = qOrbUp(:,iStart:iEnd,:)
    end do

    @:ASSERT(allocated(qBlock) .eqv. allocated(qBlockUp))
    if (allocated(qBlock)) then
      do ii = 1, transpar%ncont
        iStart = transpar%contacts(ii)%idxrange(1)
        iEnd = transpar%contacts(ii)%idxrange(2)
        qBlock(:,:,iStart:iEnd,:) = qBlockUp(:,:,iStart:iEnd,:)
      end do
    end if

  end subroutine overrideContactCharges

#:endif


end module dftbp_dftbplus_initprogram
