!--------------------------------------------------------------------------------------------------!
!  DFTB+: general package for performing fast atomistic simulations                                !
!  Copyright (C) 2006 - 2022  DFTB+ developers group                                               !
!                                                                                                  !
!  See the LICENSE file for terms of usage and distribution.                                       !
!--------------------------------------------------------------------------------------------------!

#:include 'common.fypp'

!> Contains data type representing the input data for DFTB
module dftbp_dftbplus_inputdata
  use dftbp_common_accuracy, only : dp, lc
  use dftbp_common_hamiltoniantypes, only : hamiltonianTypes
  use dftbp_dftb_dftbplusu, only : TDftbUInp
  use dftbp_dftb_dispersions, only : TDispersionInp
  use dftbp_dftb_elstatpot, only : TElStatPotentialsInp
  use dftbp_dftb_etemp, only : fillingTypes
  use dftbp_dftb_extfields, only : TElecFieldInput
  use dftbp_dftb_h5correction, only : TH5CorrectionInput
  use dftbp_dftb_repulsive_chimesrep, only : TChimesRepInp
  use dftbp_dftb_repulsive_pairrepulsive, only : TPairRepulsiveItem
  use dftbp_dftb_pmlocalisation, only : TPipekMezeyInp
  use dftbp_dftb_potentials, only : TAtomExtPotInput
  use dftbp_dftb_slakocont, only : TSlakoCont
  use dftbp_dftbplus_input_geoopt, only : TGeoOptInput
  use dftbp_elecsolvers_elecsolvers, only : TElectronicSolverInp
  use dftbp_extlibs_poisson, only : TPoissonInfo
  use dftbp_extlibs_tblite, only : TTBLiteInput
  use dftbp_io_message, only : error, warning
  use dftbp_md_xlbomd, only : TXLBOMDInp
  use dftbp_reks_reks, only : TReksInp
  use dftbp_solvation_cm5, only : TCM5Input
  use dftbp_solvation_solvinput, only : TSolvationInp
  use dftbp_timedep_linresp, only : TLinrespini
  use dftbp_timedep_pprpa, only : TppRPAcal
  use dftbp_timedep_timeprop, only : TElecDynamicsInp
  use dftbp_type_commontypes, only : TOrbitals
  use dftbp_type_linkedlist, only : TListIntR1, destruct
  use dftbp_type_typegeometry, only : TGeometry
  use dftbp_type_wrappedintr, only : TWrappedInt1
#:if WITH_SOCKETS
  use dftbp_io_ipisocket, only : IpiSocketCommInp
#:endif
#:if WITH_TRANSPORT
  use dftbp_transport_negfvars, only : TNEGFTunDos, TNEGFGreenDensInfo, TTransPar
#:endif
  implicit none

  private
  public :: TControl, TSlater, TInputData, TParallelOpts
  public :: TBlacsOpts
  public :: TRangeSepInp
  public :: init, destruct
  public :: TNEGFInfo


  !> Contains Blacs specific options.
  type :: TBlacsOpts

    !> Block size for matrix rows and columns.
    integer :: blockSize

  end type TBlacsOpts


  !> Contains the parallel options
  type :: TParallelOpts

    !> Number of processor groups
    integer :: nGroup

    !> Blacs options
    type(TBlacsOpts) :: blacsOpts

    !> Whether hybrid parallelisation is enable
    logical :: tOmpThreads

  end type TParallelOpts


  !> LBFGS input settings
  type TLbfgsInput

    !> Number of stored steps
    integer :: memory

    !> Is a line search followed along quasi-Newton directions
    logical :: isLineSearch

    !> Should the quasi-Newton step be limited?
    logical :: MaxQNStep

    !> If performing line search, should the original implementation be used
    logical :: isOldLS

  end type TLbfgsInput


  !> Range separation input
  type TRangeSepInp

    !> Threshold for integral screening
    real(dp) :: screeningThreshold

    !> Reduction of cutoff in spatial screening
    real(dp) :: cutoffRed

    !> Separation parameter
    real(dp) :: omega

    !> Choice of range separation method
    integer :: rangeSepAlg

  end type TRangeSepInp


  !> Main control data for program as extracted by the parser
  type TControl

    !> Choice of electronic hamiltonian
    integer :: hamiltonian = hamiltonianTypes%none

    !> random number generator seed
    integer :: iSeed       = 0

    !> maximum force for geometry convergence
    real(dp) :: maxForce    = 0.0_dp

    !> SCC calculation?
    logical :: tScc        = .false.

    !> l-shell resolved SCC
    logical :: tShellResolved = .false.

    !> SCC tolerance
    real(dp) :: sccTol = 0.0_dp

    !> Read starting charges from disc
    logical :: tReadChrg = .false.

    logical :: tSkipChrgChecksum = .false.

    !> Disc charges are stored as ascii or binary files
    logical :: tReadChrgAscii = .true.

    !> Disc charges should be written as ascii or binary files
    logical :: tWriteChrgAscii = .true.

    !> should probably be packaged
    logical :: isGeoOpt = .false.

    !> coordinate optimisation
    logical :: tCoordOpt   = .false.

    !> maximum line search step for atoms
    real(dp) :: maxAtomDisp = 0.2_dp

    !> should probably be packaged
    logical :: tLatOpt     = .false.

    !> Fix angles during lattice optimisation
    logical :: tLatOptFixAng = .false.

    !> Fix lengths of specified vectors
    logical :: tLatOptFixLen(3) = .false.

    !> Isotropically scale instead
    logical :: tLatOptIsotropic = .false.

    !> maximum possible linesearch step
    real(dp) :: maxLatDisp = 0.2_dp

    !> add new geometries at the end of files
    logical :: tAppendGeo  = .false.

    !> use converged SCC charges for properties like forces or charge dependent dispersion
    logical :: isSccConvRequired = .true.

    !> geometry step
    integer :: iGeoOpt     = 0

    !> used for gDIIS
    real(dp) :: deltaGeoOpt = 0.0_dp

    !> used for gDIIS
    integer :: iGenGeoOpt = 0

    !> internal variable for requirement of Mulliken analysis
    logical :: tMulliken = .false.

    !> printout of Mulliken
    logical :: tPrintMulliken = .false.

    !> Net atomic charges (i.e. on-site only part of Mulliken charges)
    logical :: tNetAtomCharges = .false.

    !> Should net atomic charges be printed
    logical :: tPrintNetAtomCharges = .false.

    !> Input for CM5 corrected Mulliken charges
    type(TCM5Input), allocatable :: cm5Input

    !> electrostatic potential evaluation and printing
    type(TElStatPotentialsInp), allocatable :: elStatPotentialsInp

    !> Localise electronic states
    logical :: tLocalise   = .false.

    !> Input data for Pipek-Mezey localisation
    type(TPipekMezeyInp), allocatable :: pipekMezeyInp

    !> Is a perturbation expression in use
    logical :: isDFTBPT = .false.

    !> Tolerance for idenfifying need for degenerate perturbation theory
    real(dp) :: tolDegenDFTBPT = 128.0_dp

    !> Is this is a static electric field perturbation calculation
    logical :: isEPerturb = .false.

    !> Frequencies for perturbation (0 being static case)
    real(dp), allocatable :: dynEFreq(:)

    !> Frequency dependent perturbation eta
    real(dp), allocatable :: etaFreq

    !> Is the response kernel (and frontier eigenvalue derivatives) calculated by perturbation
    logical :: isRespKernelPert = .false.

    !> Is the response kernel evaluated at the RPA level, or (if SCC) self-consistent
    logical :: isRespKernelRPA

    !> Frequencies for perturbation (0 being static case)
    real(dp), allocatable :: dynKernelFreq(:)

    !> printing of atom resolved energies
    logical :: tAtomicEnergy = .false.

    !> print eigenvectors to disc
    logical :: tPrintEigVecs  = .false.

    !> text file of eigenvectors?
    logical :: tPrintEigVecsTxt = .false.

    !> project eigenvectors spatially
    logical :: tProjEigenvecs = .false.

    !> Evaluate forces
    logical :: tForces = .false.

    !> Evaluate force contributions from the excited state if required and (tForces)
    logical :: tCasidaForces = .false.

    !> force evaluation method
    integer :: forceType

    !> Output forces
    logical :: tPrintForces = .false.

    !> method for calculating derivatives
    integer :: iDerivMethod = 0

    !> 1st derivative finite difference step
    real(dp) :: deriv1stDelta = 0.0_dp


    !> Molecular dynamics
    logical :: tMD         = .false.

    !> Use Plumed
    logical :: tPlumed = .false.

    !> Finite difference derivatives calculation?
    logical :: tDerivs     = .false.

    !> Should central cell coordinates be output?
    logical :: tShowFoldedCoord

    real(dp) :: nrChrg        = 0.0_dp
    real(dp) :: nrSpinPol     = 0.0_dp
    logical :: tSpin         = .false.
    logical :: tSpinSharedEf = .false.
    logical :: tSpinOrbit    = .false.
    logical :: tDualSpinOrbit = .false.
    logical :: t2Component   = .false.

    !> initial spin pattern
    real(dp), allocatable :: initialSpins(:,:)

    !> initial charges
    real(dp), allocatable :: initialCharges(:)

    !> Electronic/eigenvalue solver options
    type(TElectronicSolverInp) :: solver

    integer :: iMixSwitch    = 0
    integer :: maxIter       = 0
    real(dp) :: almix         = 0.0_dp
    integer :: iGenerations  = 0
    logical :: tFromStart    = .true.
    real(dp) :: broydenOmega0 = 0.01_dp
    real(dp) :: broydenMinWeight = 1.0_dp
    real(dp) :: broydenMaxWeight = 1.0e5_dp
    real(dp) :: broydenWeightFac = 1.0e-2_dp
    real(dp) :: andersonInitMixing = 0.01_dp
    integer :: andersonNrDynMix = 0
    real(dp), allocatable :: andersonDynMixParams(:,:)
    real(dp) :: andersonOmega0 = 1.0e-2_dp
    integer :: nrMoved       = 0
    integer, allocatable :: indMovedAtom(:)
    integer, allocatable :: indDerivAtom(:)
    integer :: nrConstr      = 0
    integer, allocatable :: conAtom(:)
    real(dp), allocatable :: conVec(:,:)
    character(lc) :: outFile       = ''

    !> do we have MD velocities
    logical :: tReadMDVelocities = .false.

    !> initial MD velocities
    real(dp), allocatable :: initialVelocities(:,:)
    real(dp) :: deltaT        = 0.0_dp

    real(dp) :: tempAtom      = 0.0_dp
    integer :: iThermostat   = 0

    !> whether to initialize internal state of the Nose-Hoover thermostat from input
    logical :: tInitNHC = .false.
    real(dp), allocatable :: xnose(:)
    real(dp), allocatable :: vnose(:)
    real(dp), allocatable :: gnose(:)


    !> whether to shift to a co-moving frame for MD
    logical :: tMDstill
    logical :: tRescale = .false.
    integer, allocatable :: tempMethods(:)
    integer, allocatable :: tempSteps(:)
    real(dp), allocatable :: tempValues(:)
    logical :: tSetFillingTemp = .false.

    real(dp) :: tempElec      = 0.0_dp
    logical :: tFixEf        = .false.
    real(dp), allocatable :: Ef(:)
    logical :: tFillKSep     = .false.
    integer :: iDistribFn    = fillingTypes%Fermi
    real(dp) :: wvScale       = 0.0_dp

    !> default chain length for Nose-Hoover
    integer :: nh_npart      = 3

    !> default order of NH integration
    integer :: nh_nys        = 3

    !> default multiple time steps for N-H propagation
    integer :: nh_nc         = 1

    integer :: maxRun        = -2


    !> second derivative finite difference step
    real(dp) :: deriv2ndDelta    = 0.0_dp

    integer :: nKPoint       = 0
    real(dp), allocatable :: kPoint(:,:)
    real(dp), allocatable :: kWeight(:)


    !> cell pressure if periodic
    real(dp) :: pressure       = 0.0_dp
    logical :: tBarostat = .false.

    !> use isotropic scaling if barostatting
    logical :: tIsotropic = .true.
    real(dp) :: BarostatStrength = 0.0_dp


    !> read atomic masses from the input not the SK data
    real(dp), allocatable :: masses(:)


    !> spin constants
    real(dp), allocatable :: spinW(:,:,:)

    !> customised Hubbard U values
    real(dp), allocatable :: hubbU(:,:)

    !> spin-orbit constants
    real(dp), allocatable :: xi(:,:)

    !> DFTB+U input, if present
    type(TDftbUInp), allocatable :: dftbUInp

    !> Correction to energy from on-site matrix elements
    real(dp), allocatable :: onSiteElements(:,:,:,:)

    !> Correction to dipole momements on-site matrix elements
    real(dp), allocatable :: onSiteDipole(:,:)

    !> Number of external charges
    integer :: nExtChrg = 0

    !> external charge values and locations
    real(dp), allocatable :: extChrg(:,:)

    !> finite charge width if needed
    real(dp), allocatable :: extChrgBlurWidth(:)

    !> Homogeneous external electric field
    type(TElecFieldInput), allocatable :: electricField

    !> Potential(s) at atomic sites
    type(TAtomExtPotInput), allocatable :: atomicExtPotential

    !> Projection of eigenvectors
    type(TListIntR1) :: iAtInRegion
    logical, allocatable :: tShellResInRegion(:)
    logical, allocatable :: tOrbResInRegion(:)
    character(lc), allocatable :: RegionLabel(:)


    !> H short range damping
    logical :: tDampH = .false.
    real(dp) :: dampExp = 0.0_dp

    type(TH5CorrectionInput), allocatable :: h5Input

    !> Halogen X correction
    logical :: tHalogenX = .false.

    !> Old repulsive
    logical :: useBuggyRepSum


    !> Old kinetic energy stress contribution in MD
    logical :: useBuggyKEStress = .false.


    !> Ewald alpha
    real(dp) :: ewaldAlpha = 0.0_dp

    !> Ewald tolerance
    real(dp) :: tolEwald = 1.0E-9_dp

    !> Various options
    logical :: tWriteTagged = .false.

    !> Nr. of SCC iterations without restart info
    integer :: restartFreq  = 20
    logical :: tWriteDetailedXML = .false.
    logical :: tWriteResultsTag = .false.
    logical :: tWriteDetailedOut = .true.
    logical :: tWriteBandDat = .true.
    logical :: oldSKInter = .false.
    logical :: tWriteHS = .false.
    logical :: tWriteRealHS = .false.
    logical :: tMinMemory = .false.

    !> potential shifts are read from file
    logical :: tReadShifts = .false.
    !> potential shifts are written on file
    logical :: tWriteShifts = .false.

    !> use Poisson solver for electrostatics
    logical :: tPoisson = .false.


    !> Dispersion related stuff
    type(TDispersionInp), allocatable :: dispInp

    !> Solvation
    class(TSolvationInp), allocatable :: solvInp

    !> Rescaling of electric fields (applied or dipole) if the system is solvated
    logical :: isSolvatedFieldRescaled = .false.

    !> Input for tblite library
    type(TTBLiteInput), allocatable :: tbliteInp


    !> Local potentials
    real(dp), allocatable :: chrgConstr(:,:)
    real(dp), allocatable :: thirdOrderOn(:,:)


    !> 3rd order
    real(dp), allocatable :: hubDerivs(:,:)
    logical :: t3rd = .false.
    logical :: t3rdFull = .false.


    !> XLBOMD
    type(TXLBOMDInp), allocatable :: xlbomd

    !> TD Linear response input
    type(TLinrespini), allocatable :: lrespini

    !> ElectronDynamics
    type(TElecDynamicsInp), allocatable :: elecDynInp

    !> input for particle-particle RPA
    type(TppRPAcal), allocatable :: ppRPA

    !> LBFGS input
    type(TLbfgsInput), allocatable :: lbfgsInp

    !> Geometry optimizer input
    type(TGeoOptInput), allocatable :: geoOpt

    !> Range separated input
    type(TRangeSepInp), allocatable :: rangeSepInp

  #:if WITH_SOCKETS
    !> socket communication
    type(ipiSocketCommInp), allocatable :: socketInput
  #:endif

    type(TParallelOpts), allocatable :: parallelOpts

    !> Maximal timing level to show in output
    integer :: timingLevel

    ! Custom occupations
    type(TWrappedInt1), allocatable :: customOccAtoms(:)
    real(dp), allocatable :: customOccFillings(:,:)

    ! TI-DFTB variables

    !> Non-Aufbau filling
    logical :: isNonAufbau = .false.

    !> SpinPurify
    logical :: isSpinPurify = .false.

    !> GroundGuess
    logical :: isGroundGuess = .false.

    !> REKS input
    type(TReksInp) :: reksInp

    !> Whether Scc should be updated with the output charges (obtained after diagonalization)
    !> Could be set to .false. to prevent costly recalculations (e.g. when using Poisson-solver)
    logical :: updateSccAfterDiag = .true.

    !> Write cavity information as COSMO file
    logical :: tWriteCosmoFile = .false.

    !> Whether ChIMES correction for repulsives should be applied.
    type(TChimesRepInp), allocatable :: chimesRepInput

  end type TControl


  !> Slater-Koster data
  type TSlater
    real(dp), allocatable :: skSelf(:, :)
    real(dp), allocatable :: skHubbU(:, :)
    real(dp), allocatable :: skOcc(:, :)
    real(dp), allocatable :: mass(:)

    type(TSlakoCont), allocatable :: skHamCont
    type(TSlakoCont), allocatable :: skOverCont
    type(TOrbitals), allocatable :: orb
    type(TPairRepulsiveItem), allocatable :: pairRepulsives(:,:)

  end type TSlater

#:if WITH_TRANSPORT

  !> container for data needed by libNEGF
  type TNEGFInfo
    type(TNEGFTunDos) :: tundos  !Transport section informations
    type(TNEGFGreenDensInfo) :: greendens  !NEGF solver section informations
  end type TNEGFInfo

#:else

  !> Dummy type replacement
  type TNegfInfo
  end type TNegfInfo

#:endif


  !> container for input data constituents
  type TInputData
    logical :: tInitialized = .false.
    type(TControl) :: ctrl
    type(TGeometry) :: geom
    type(TSlater) :: slako
  #:if WITH_TRANSPORT
    type(TTransPar) :: transpar
    type(TNEGFInfo) :: ginfo
  #:endif
    type(TPoissonInfo) :: poisson
  end type TInputData


  !> Initialise the input data
  interface init
    module procedure InputData_init
  end interface init


  !> destroy input data for variables that do not go out of scope
  interface destruct
    module procedure InputData_destruct
  end interface destruct

contains


  !> Mark data structure as initialised
  subroutine InputData_init(this)

    !> Instance
    type(TInputData), intent(out) :: this

    this%tInitialized = .true.

  end subroutine InputData_init


  !> destructor for parts that are not cleaned up when going out of scope
  subroutine InputData_destruct(this)

    !> Instance
    type(TInputData), intent(inout) :: this

    call Control_destruct(this%ctrl)

  end subroutine InputData_destruct


  !> destructor for parts that are not cleaned up when going out of scope
  subroutine Control_destruct(this)

    !> Instance
    type(TControl), intent(inout) :: this

    if (allocated(this%tShellResInRegion)) then
      call destruct(this%iAtInRegion)
    end if

  end subroutine Control_destruct

end module dftbp_dftbplus_inputdata
