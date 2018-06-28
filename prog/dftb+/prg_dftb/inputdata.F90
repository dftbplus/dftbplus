!--------------------------------------------------------------------------------------------------!
!  DFTB+: general package for performing fast atomistic simulations                                !
!  Copyright (C) 2018  DFTB+ developers group                                                      !
!                                                                                                  !
!  See the LICENSE file for terms of usage and distribution.                                       !
!--------------------------------------------------------------------------------------------------!

#:include 'common.fypp'

!> Contains data type representing the input data for DFTB
module inputdata_module
  use assert
  use accuracy
  use typegeometry
  use message
  use dispersions, only : DispersionInp
  use linresp_module, only : linrespini
  use slakocont
  use commontypes
  use repcont
  use linkedlist
  use xlbomd_module
#:if WITH_SOCKETS
  use ipisocket, only : IpiSocketCommInp
#:endif
  use pmlocalisation, only : TPipekMezeyInp
  use elstatpot, only : TElStatPotentialsInp
  implicit none
  private
  save

  public :: control, TGeometry, slater, inputData, XLBOMDInp, TParallelOpts
  public :: TBlacsOpts
  public :: init, destruct


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
    integer :: memory
  end type TLbfgsInput


  !> Main control data for program as extracted by the parser
  type control

    !> random number generator seed
    integer :: iSeed       = 0

    !> maximum force for geometry convergence
    real(dp) :: maxForce    = 0.0_dp

    !> SCC calculation?
    logical :: tScc        = .false.

    !> l-shell resolved SCC
    logical :: tOrbResolved = .false.

    !> SCC tolerance
    real(dp) :: sccTol      = 0.0_dp

    !> Read starting charges from disc
    logical :: tReadChrg = .false.

    logical :: tSkipChrgChecksum = .false.

    !> Disc charges are stored as ascii or binary files
    logical :: tReadChrgAscii = .true.

    !> Disc charges should be written as ascii or binary files
    logical :: tWriteChrgAscii = .true.

    !> should probably be packaged
    logical :: tGeoOpt     = .false.

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

    !> use converged SCC forces only
    logical :: tConvrgForces = .true.

    !> geometry step
    integer :: iGeoOpt     = 0

    !> used for gDIIS
    real(dp) :: deltaGeoOpt = 0.0_dp

    !> used for gDIIS
    integer :: iGenGeoOpt = 0

    !> internal variable for requirement of Mulliken analysis
    logical :: tMulliken   = .false.

    !> printout of Mulliken
    logical :: tPrintMulliken   = .false.

    !> electrostatic potential evaluation and printing
    type(TElStatPotentialsInp), allocatable :: elStatPotentialsInp

    !> Localise electronic states
    logical :: tLocalise   = .false.

    type(TPipekMezeyInp), allocatable :: pipekMezeyInp

    !> printing of atom resolved energies
    logical :: tAtomicEnergy = .false.

    !> print eigenvectors to disc
    logical :: tPrintEigVecs  = .false.

    !> text file of eigenvectors?
    logical :: tPrintEigVecsTxt = .false.

    !> project eigenvectors spatially
    logical :: tProjEigenvecs = .false.

    !> Evaluate forces
    logical :: tForces     = .false.

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
    logical :: tDFTBU        = .false.
    integer :: iSolver       = 0
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
    real(dp) :: Ef(2)         = 0.0_dp
    logical :: tFillKSep     = .false.
    integer :: iDistribFn    = 0
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


    !> cell presure if periodic
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


    !> choice of the DFTB+U functional
    integer :: DFTBUfunc     = 0

    !> list of U-J for species
    real(dp), allocatable :: UJ(:,:)

    !> How many U-J for each species
    integer, allocatable :: nUJ(:)

    !> number of l-values of U-J for each block
    integer, allocatable :: niUJ(:,:)

    !> l-values of U-J for each block
    integer, allocatable :: iUJ(:,:,:)


    !> Number of external charges
    integer :: nExtChrg = 0

    !> external charge values and locations
    real(dp), allocatable :: extChrg(:,:)

    !> finite charge width if needed
    real(dp), allocatable :: extChrgBlurWidth(:)


    !> External homogeneous electric field
    logical :: tEField = .false.

    !> time dependent field in MD
    logical :: tTDEfield = .false.

    !> strength
    real(dp) :: EFieldStrength = 0.0_dp

    !> direction
    real(dp) :: EfieldVector(3)

    !> frequency of time dependent field
    real(dp) :: EfieldOmega

    !> relative phase of field
    integer :: EfieldPhase = 0


    !> Projection of eigenvectors
    type(listIntR1) :: iAtInRegion
    logical, allocatable :: tShellResInRegion(:)
    logical, allocatable :: tOrbResInRegion(:)
    character(lc), allocatable :: RegionLabel(:)


    !> H short range damping
    logical :: tDampH = .false.
    real(dp) :: dampExp = 0.0_dp


    ! H5 correction
    !> H5 correction On/Off(default) flag
    logical ::h5SwitchedOn = .false.
    !> Global parameters - set to -1 to identify they were not initialized
    real(dp) :: h5RScale = -1.0_dp
    real(dp) :: h5WScale = -1.0_dp
    real(dp), allocatable :: h5ElementPara(:)
    ! H5 correction end


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


    !> Dispersion related stuff
    type(DispersionInp), allocatable :: dispInp


    !> Local potentials
    real(dp), allocatable :: chrgConstr(:,:)
    real(dp), allocatable :: thirdOrderOn(:,:)


    !> 3rd order
    real(dp), allocatable :: hubDerivs(:,:)
    logical :: t3rd, t3rdFull


    !> XLBOMD
    type(XLBOMDInp), allocatable :: xlbomd

    type(linrespini) :: lrespini

    !> LBFGS input
    type(TLbfgsInput), allocatable :: lbfgsInp


  #:if WITH_SOCKETS
    !> socket communication
    type(IpiSocketCommInp), allocatable :: socketInput
  #:endif

    type(TParallelOpts), allocatable :: parallelOpts

    !> Maximal timing level to show in output
    integer :: timingLevel

  end type control


  !> Slater-Koster data
  type slater
    real(dp), allocatable :: skSelf(:, :)
    real(dp), allocatable :: skHubbU(:, :)
    real(dp), allocatable :: skOcc(:, :)
    real(dp), allocatable :: mass(:)

    type(OSlakoCont), allocatable :: skHamCont
    type(OSlakoCont), allocatable :: skOverCont
    type(ORepCont), allocatable :: repCont
    type(TOrbitals), allocatable :: orb
  end type slater


  !> container for input data constituents
  type inputData
    type(control) :: ctrl
    type(TGeometry) :: geom
    type(slater) :: slako
    logical :: tInitialized = .false.
  end type inputData


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
  subroutine InputData_init(self)
    type(inputData), intent(out) :: self

    self%tInitialized = .true.

  end subroutine InputData_init


  !> destructor for parts that are not cleaned up when going out of scope
  subroutine InputData_destruct(self)
    type(inputData), intent(inout) :: self

    call Control_destruct(self%ctrl)

  end subroutine InputData_destruct


  !> destructor for parts that are not cleaned up when going out of scope
  subroutine Control_destruct(self)
    type(control), intent(inout) :: self

    if (allocated(self%tShellResInRegion)) then
      call destruct(self%iAtInRegion)
    end if

  end subroutine Control_destruct

end module inputdata_module
