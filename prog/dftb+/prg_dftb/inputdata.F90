!--------------------------------------------------------------------------------------------------!
!  DFTB+: general package for performing fast atomistic simulations                                !
!  Copyright (C) 2017  DFTB+ developers group                                                      !
!                                                                                                  !
!  See the LICENSE file for terms of usage and distribution.                                       !
!--------------------------------------------------------------------------------------------------!

#:include 'common.fypp'

!!* Contains data type representing the input data for DFTB
!!* @todo
!!*   Has to be replaced by a more general approach, which should be
!!*   a kind of database (quite similar to Thomas' global library, but !'
!!*   with a dynamical approach).
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
  use ipisocket, only : IpiSocketCommInp
  implicit none
  private
  save

  type control
    integer       :: iSeed       = 0
    real(dp)      :: maxForce    = 0.0_dp
    logical       :: tScc        = .false.
    logical       :: tOrbResolved = .false. ! l-shell resolved SCC
    real(dp)      :: sccTol      = 0.0_dp
    logical       :: tReadChrg   = .false.
    logical       :: tGeoOpt     = .false. ! should probably be packaged
    logical       :: tCoordOpt   = .false.
    real(dp)      :: maxAtomDisp = 0.2_dp ! maximum line search step for atoms
    logical       :: tLatOpt     = .false. ! should probably be packaged
    logical       :: tLatOptFixAng = .false.
    logical       :: tLatOptFixLen(3) = .false.
    logical       :: tLatOptIsotropic = .false.
    real(dp)      :: maxLatDisp = 0.2_dp ! maximum possible linesearch step
    logical       :: tAppendGeo  = .false.
    logical       :: tConvrgForces = .true.
    integer       :: iGeoOpt     = 0
    real(dp)      :: deltaGeoOpt = 0.0_dp ! used for gDIIS
    integer       :: iGenGeoOpt = 0 ! used for gDIIS
    logical       :: tMulliken   = .false. ! internal requirement for Mulliken
    logical       :: tPrintMulliken   = .false. ! printout of Mulliken
    logical       :: tLocalise   = .false.
    logical       :: tPipekMezey = .false.
    logical       :: tPipekDense = .false.
    real(dp), allocatable :: sparsePipekTols(:)
    real(dp)      :: PipekTol
    integer       :: PipekMaxIter = 100 ! cycles to localise charges
    logical       :: tAtomicEnergy = .false. !
    logical       :: tPrintEigVecs  = .false. ! print eigenvectors to disc
    logical       :: tPrintEigVecsTxt = .false. ! text file of eigenvectors?
    logical       :: tProjEigenvecs = .false. !project eigenvectors spatially

    logical       :: tForces     = .false. ! Evaluate forces
    integer :: forceType
    logical       :: tPrintForces = .false. ! Output forces
    integer       :: iDerivMethod = 0
    real(dp)      :: deriv1stDelta = 0.0_dp ! 1st derivative finite
                                            ! difference step

    logical       :: tMD         = .false. ! Molecular dynamics
    logical       :: tDerivs     = .false. ! Finite difference
                                           ! derivatives calculation?
    logical       :: tShowFoldedCoord !* Should central cell coordinates
    !* be output?

    real(dp)               :: nrChrg        = 0.0_dp
    real(dp)               :: nrSpinPol     = 0.0_dp
    logical                :: tSpin         = .false.
    logical                :: tSpinSharedEf = .false.
    logical                :: tSpinOrbit    = .false.
    logical                :: tDualSpinOrbit = .false.
    logical                :: t2Component   = .false.
    real(dp), allocatable  :: initialSpins(:,:)   !initial spin pattern
    real(dp), allocatable  :: initialCharges(:)
    logical                :: tDFTBU        = .false.
    integer                :: iSolver       = 0
    integer                :: iMixSwitch    = 0
    integer                :: maxIter       = 0
    real(dp)               :: almix         = 0.0_dp
    integer                :: iGenerations  = 0
    logical                :: tFromStart    = .true.
    real(dp)               :: broydenOmega0 = 0.01_dp
    real(dp)               :: broydenMinWeight = 1.0_dp
    real(dp)               :: broydenMaxWeight = 1.0e5_dp
    real(dp)               :: broydenWeightFac = 1.0e-2_dp
    real(dp)               :: andersonInitMixing = 0.01_dp
    integer                :: andersonNrDynMix = 0
    real(dp), allocatable  :: andersonDynMixParams(:,:)
    real(dp)               :: andersonOmega0 = 1.0e-2_dp
    integer                :: nrMoved       = 0
    integer, allocatable   :: indMovedAtom(:)
    integer                :: nrConstr      = 0
    integer, allocatable   :: conAtom(:)
    real(dp), allocatable  :: conVec(:,:)
    character(lc)          :: outFile       = ''
    logical                :: tReadMDVelocities = .false. ! do we have MD velocities
    real(dp), allocatable  :: initialVelocities(:,:) ! initial MD velocities
    real(dp)               :: deltaT        = 0.0_dp

    real(dp)               :: tempAtom      = 0.0_dp
    integer                :: iThermostat   = 0
    logical                :: tInitNHC = .false. ! whether to initialize
    !  internal state of the Nose-Hoover thermostat from input
    real(dp), allocatable :: xnose(:)
    real(dp), allocatable :: vnose(:)
    real(dp), allocatable :: gnose(:)

    logical                :: tMDstill ! whether to shift to a co-moving
    ! frame for MD
    logical                :: tRescale = .false.
    integer, allocatable   :: tempMethods(:)
    integer, allocatable   :: tempSteps(:)
    real(dp), allocatable  :: tempValues(:)
    logical                :: tSetFillingTemp = .false.

    real(dp)               :: tempElec      = 0.0_dp
    logical                :: tFixEf        = .false.
    real(dp)               :: Ef(2)         = 0.0_dp
    logical                :: tFillKSep     = .false.
    integer                :: iDistribFn    = 0
    real(dp)               :: wvScale       = 0.0_dp
    integer                :: nh_npart      = 3 ! chain length for Nose-Hoover
    integer                :: nh_nys        = 3 ! order of integration
    integer                :: nh_nc         = 1 ! multiple time steps for N-H
    !  propagation

    integer                :: maxRun        = -2

    real(dp)               :: deriv2ndDelta    = 0.0_dp ! second
    ! derivative finite difference step

    integer                :: nKPoint       = 0
    real(dp), allocatable :: kPoint(:,:)
    real(dp), allocatable :: kWeight(:)

    real(dp)               :: pressure       = 0.0_dp ! cell presure if periodic
    logical                :: tBarostat = .false.
    logical                :: tIsotropic = .true. ! use isotropic scaling if
    !  barostatting
    real(dp)               :: BarostatStrength = 0.0_dp

    ! read atomic masses from the input not the SK data
    real(dp), allocatable :: masses(:)

    real(dp), allocatable :: spinW(:,:,:)  ! spin constants
    real(dp), allocatable :: hubbU(:,:)    ! customised Hubbard U values
    real(dp), allocatable :: xi(:,:)       ! spin-orbit constants

    integer                :: DFTBUfunc     = 0 ! choice of the DFTB+U
    ! functional
    real(dp), allocatable :: UJ(:,:)     ! list of U-J for species
    integer,  allocatable :: nUJ(:)      ! How many U-J for each
    ! species
    integer,  allocatable :: niUJ(:,:)   ! number of l-values of
    ! U-J for each block
    integer,  allocatable :: iUJ(:,:,:)  ! l-values of U-J for each
    ! block

    !! External charges
    integer :: nExtChrg = 0
    real(dp), allocatable :: extChrg(:,:)
    real(dp), allocatable :: extChrgBlurWidth(:)

    !! External electric field
    logical  :: tEField = .false.
    logical  :: tTDEfield = .false.
    real(dp) :: EFieldStrength = 0.0_dp
    real(dp) :: EfieldVector(3)
    real(dp) :: EfieldOmega
    integer  :: EfieldPhase = 0

    !! Projection of eigenvectors
    type(listIntR1) :: iAtInRegion
    logical, allocatable :: tShellResInRegion(:)
    logical, allocatable :: tOrbResInRegion(:)
    character(lc), allocatable :: RegionLabel(:)

    !! H short range damping
    logical :: tDampH = .false.
    real(dp) :: dampExp = 0.0_dp

    !! Old repulsive
    logical :: useBuggyRepSum

    !! Old kinetic energy stress contribution in MD
    logical :: useBuggyKEStress = .false.

    !! Ewald alpha
    real(dp) :: ewaldAlpha = 0.0_dp

    !! Various options
    logical :: tWriteTagged = .false.
    integer :: restartFreq  = 20 ! Nr. of SCC iterations without restart info
    logical :: tWriteDetailedXML = .false.
    logical :: tWriteResultsTag = .false.
    logical :: tWriteDetailedOut = .true.
    logical :: tWriteBandDat = .true.
    logical :: oldSKInter = .false.
    logical :: tWriteHS = .false.
    logical :: tWriteRealHS = .false.
    logical :: tMinMemory = .false.

    !! Dispersion related stuff
    type(DispersionInp), allocatable :: dispInp

    !! Local potentials
    real(dp), allocatable :: chrgConstr(:,:)
    real(dp), allocatable :: thirdOrderOn(:,:)

    !! 3rd order
    real(dp), allocatable :: hubDerivs(:,:)
    logical :: t3rd, t3rdFull

    !! XLBOMD
    type(XLBOMDInp), allocatable :: xlbomd

    type(linrespini) :: lrespini

    !! socket communication
    type(IpiSocketCommInp), allocatable :: socketInput
  end type control


  type geometry
    integer           :: nrAtoms         = 0
    logical           :: tPeriodic       = .false.
    logical           :: tFracCoord      = .false.
    integer,  allocatable :: types(:)
    real(dp), allocatable :: coords(:, :)
    integer           :: nrTypes         = 0
    real(dp), allocatable :: origo(:)
    real(dp), allocatable :: latVecs(:, :)
    character(mc), allocatable :: speciesName(:)
  end type geometry


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


  type inputData
    type(control)  :: ctrl
    type(TGeometry) :: geom
    type(slater)   :: slako
    logical        :: tInitialized = .false.
  end type inputData



  interface init
    module procedure InputData_init
  end interface init

  interface destruct
    module procedure InputData_destruct
  end interface destruct


  public :: control, TGeometry, slater, inputData, XLBOMDInp
  public :: init, destruct


contains

  subroutine InputData_init(self)
    type(inputData), intent(out) :: self

    self%tInitialized = .true.

  end subroutine InputData_init


  subroutine InputData_destruct(self)
    type(inputData), intent(inout) :: self

    call Control_destruct(self%ctrl)

  end subroutine InputData_destruct


  subroutine Control_destruct(self)
    type(control), intent(inout) :: self

    if (allocated(self%tShellResInRegion)) then
      call destruct(self%iAtInRegion)
    end if

  end subroutine Control_destruct



end module inputdata_module
