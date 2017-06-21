!--------------------------------------------------------------------------------------------------!
!  DFTB+: general package for performing fast atomistic simulations                                !
!  Copyright (C) 2017  DFTB+ developers group                                                      !
!                                                                                                  !
!  See the LICENSE file for terms of usage and distribution.                                       !
!--------------------------------------------------------------------------------------------------!

!!* Contains data type representing the input data for DFTB
!!* @todo
!!*   Has to be replaced by a more general approach, which should be
!!*   a kind of database (quite similar to Thomas' global library, but !'
!!*   with a dynamical approach).
module inputdata_module
#include "assert.h"
#include "allocate.h"  
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
    real(dp), pointer :: sparsePipekTols(:) => null()
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
    real(dp), pointer      :: initialSpins(:,:) => null() !initial spin pattern
    real(dp), pointer      :: initialCharges(:) => null() 
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
    real(dp), pointer      :: andersonDynMixParams(:,:) => null()
    real(dp)               :: andersonOmega0 = 1.0e-2_dp
    integer                :: nrMoved       = 0
    integer, pointer       :: indMovedAtom(:)  => null()
    integer                :: nrConstr      = 0
    integer, pointer       :: conAtom(:)    => null()
    real(dp), pointer      :: conVec(:,:)   => null()
    character(lc)          :: outFile       = ''
    logical                :: tReadMDVelocities = .false. ! do we have MD velocities
    real(dp), pointer      :: initialVelocities(:,:) => null() ! initial MD velocities
    real(dp)               :: deltaT        = 0.0_dp

    real(dp)               :: tempAtom      = 0.0_dp
    integer                :: iThermostat   = 0
    logical                :: tInitNHC = .false. ! whether to initialize
    !  internal state of the Nose-Hoover thermostat from input
    real(dp), pointer      :: xnose(:) => null()
    real(dp), pointer      :: vnose(:) => null()
    real(dp), pointer      :: gnose(:) => null()
    
    logical                :: tMDstill ! whether to shift to a co-moving
    ! frame for MD
    logical                :: tRescale = .false.
    integer, pointer       :: tempMethods(:) => null()
    integer, pointer       :: tempSteps(:)  => null()
    real(dp), pointer      :: tempValues(:) => null()
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
    real(dp), pointer      :: UJ(:,:)       => null() ! list of U-J for species
    integer,  pointer      :: nUJ(:)        => null() ! How many U-J for each
    ! species
    integer,  pointer      :: niUJ(:,:)     => null() ! number of l-values of
    ! U-J for each block    
    integer,  pointer      :: iUJ(:,:,:)    => null() ! l-values of U-J for each
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
    logical, pointer :: tShellResInRegion(:) => null()
    logical, pointer :: tOrbResInRegion(:) => null()
    character(lc), pointer :: RegionLabel(:) => null()
    
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
    real(dp), pointer :: chrgConstr(:,:) => null()
    real(dp), pointer :: thirdOrderOn(:,:) => null()

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
    integer,  pointer :: types(:)        => null()
    real(dp), pointer :: coords(:, :)    => null()
    integer           :: nrTypes         = 0
    real(dp), pointer :: origo(:)        => null()
    real(dp), pointer :: latVecs(:, :)   => null()
    character(mc), pointer :: speciesName(:) => null()
  end type geometry

  
  type slater
    real(dp), pointer :: skSelf(:, :)          => null()
    real(dp), pointer :: skHubbU(:, :)         => null()
    real(dp), pointer :: skOcc(:, :)           => null()
    real(dp), pointer :: mass(:)               => null()

    type(OSlakoCont), pointer :: skHamCont
    type(OSlakoCont), pointer :: skOverCont
    type(ORepCont), pointer :: repCont
    type(TOrbitals), pointer :: orb
  end type slater


  type inputData
    type(control)  :: ctrl
    type(TGeometry) :: geom
    type(slater)   :: slako
    logical        :: tInitialized = .false.
  contains
    final :: inputData_destruct
  end type inputData


  
  interface init
    module procedure InputData_init
  end interface
  
  
  public :: control, TGeometry, slater, inputData, XLBOMDInp
  public :: init

  
contains

  subroutine InputData_init(self)
    type(inputData), intent(out) :: self

    call initControl(self%ctrl)
    call initSlater(self%slako)
    self%tInitialized = .true.

  end subroutine InputData_init



  subroutine InputData_destruct(self)
    type(inputData), intent(inout) :: self

    if (.not. self%tInitialized) then
      return
    end if
    
    call destroyControl(self%ctrl)
    call destroySlater(self%slako)
    self%tInitialized = .false.

  end subroutine InputData_destruct



  !!* Initializes control data.
  !!* @param ctrl Holds control data.
  subroutine initControl(ctrl)
    type(control), intent(inout) :: ctrl
    
    INIT_PARR(ctrl%conAtom)
    INIT_PARR(ctrl%conVec)
    INIT_PARR(ctrl%initialVelocities)
    INIT_PARR(ctrl%indMovedAtom)
    
  end subroutine initControl
  
  
  
  !!* Destroys control data.
  !!* @param ctrl Holds control data.
  subroutine destroyControl(ctrl)
    type(control), intent(inout) :: ctrl

    DEALLOCATE_PARR(ctrl%initialSpins)
    DEALLOCATE_PARR(ctrl%initialCharges)
    DEALLOCATE_PARR(ctrl%andersonDynMixParams)
    DEALLOCATE_PARR(ctrl%indMovedAtom)
    DEALLOCATE_PARR(ctrl%conAtom)
    DEALLOCATE_PARR(ctrl%conVec)
    DEALLOCATE_PARR(ctrl%initialVelocities)
    DEALLOCATE_PARR(ctrl%tempMethods)
    DEALLOCATE_PARR(ctrl%tempSteps)
    DEALLOCATE_PARR(ctrl%tempValues)
    DEALLOCATE_PARR(ctrl%UJ)
    DEALLOCATE_PARR(ctrl%iUJ)
    DEALLOCATE_PARR(ctrl%niUJ)
    DEALLOCATE_PARR(ctrl%nUJ)
    
    if (associated(ctrl%xnose)) then
      DEALLOCATE_PARR(ctrl%xnose)
      DEALLOCATE_PARR(ctrl%vnose)
      DEALLOCATE_PARR(ctrl%gnose)
    end if

    if (associated(ctrl%tShellResInRegion)) then
      DEALLOCATE_PARR(ctrl%tShellResInRegion)
      DEALLOCATE_PARR(ctrl%RegionLabel)
    end if
    
    DEALLOCATE_PARR(ctrl%chrgConstr)
    DEALLOCATE_PARR(ctrl%thirdOrderOn)

    
    
  end subroutine destroyControl



  !!* Initalizes Slater-Koster data.
  !!* @param slako Holds Slater-Koster data.
  subroutine initSlater(slako)
    type(slater), intent(inout) :: slako

    INIT_PARR(slako%skSelf)
    INIT_PARR(slako%skHubbU)
    INIT_PARR(slako%skOcc)
    INIT_PARR(slako%mass)


  end subroutine initSlater



  !!* Initalizes Slater-Koster data.
  !!* @param slako Holds Slater-Koster data.
  subroutine destroySlater(slako)
    type(slater), intent(inout) :: slako

    DEALLOCATE_PARR(slako%skSelf)
    DEALLOCATE_PARR(slako%skHubbU)
    DEALLOCATE_PARR(slako%skOcc)
    DEALLOCATE_PARR(slako%mass)

  end subroutine destroySlater

end module inputdata_module
