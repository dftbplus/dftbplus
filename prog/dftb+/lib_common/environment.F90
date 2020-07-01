!--------------------------------------------------------------------------------------------------!
!  DFTB+: general package for performing fast atomistic simulations                                !
!  Copyright (C) 2018  DFTB+ developers group                                                      !
!                                                                                                  !
!  See the LICENSE file for terms of usage and distribution.                                       !
!--------------------------------------------------------------------------------------------------!

#:include 'common.fypp'

!> Contains computer environment settings
module environment
  use globalenv, only : shutdown, stdOut
  use timerarray
  use fileregistry
#:if WITH_MPI
  use mpienv
#:endif
#:if WITH_SCALAPACK
  use blacsenv
#:endif
  implicit none
  private

  public :: TEnvironment, TEnvironment_init
  public :: globalTimers


  !> Contains environment settings.
  type :: TEnvironment
    private

    !> Whether this process is the master?
    logical, public :: tGlobalMaster = .true.

    !> Nr. of groups in the system
    integer, public :: nGroup = 1

    !> Id of current group (starts with 0)
    integer, public :: myGroup = 0

    !> Global timers
    type(TTimerArray), public :: globalTimer

    !> Registry of files, which may be open and must be closed when environment is shut down
    type(TFileRegistry), public :: fileFinalizer

  #:if WITH_MPI
    !> Global mpi settings
    type(TMpiEnv), public :: mpi
  #:endif
  #:if WITH_SCALAPACK
    !> Global scalapack settings
    type(TBlacsEnv), public :: blacs
  #:endif

  contains
    procedure :: destruct => TEnvironment_destruct
    procedure :: shutdown => TEnvironment_shutdown
    procedure :: initGlobalTimer => TEnvironment_initGlobalTimer
  #:if WITH_MPI
    procedure :: initMpi => TEnvironment_initMpi
  #:endif
  #:if WITH_SCALAPACK
    procedure :: initBlacs => TEnvironment_initBlacs
  #:endif

  end type TEnvironment

  type(TTimerItem), parameter :: globalTimerItems(14) = [&
      & TTimerItem("Global initialisation", 1),&
      & TTimerItem("Pre-SCC initialisation", 1),&
      & TTimerItem("Sparse H0 and S build", 4),&
      & TTimerItem("SCC", 1),&
      & TTimerItem("Diagonalisation", 2),&
      & TTimerItem("Sparse to dense", 4),&
      & TTimerItem("Dense to sparse", 4),&
      & TTimerItem("Density matrix creation", 2),&
      & TTimerItem("Post-SCC processing", 1),&
      & TTimerItem("Eigenvector writing", 2),&
      & TTimerItem("Energy-density matrix creation", 2),&
      & TTimerItem("Force calculation", 2),&
      & TTimerItem("Stress calculation", 2),&
      & TTimerItem("Post-geometry optimisation", 1)&
      & ]

  type :: TGlobalTimersHelper
    integer :: globalInit = 1
    integer :: preSccInit = 2
    integer :: sparseH0S = 3
    integer :: scc = 4
    integer :: diagonalization = 5
    integer :: sparseToDense = 6
    integer :: denseToSparse = 7
    integer :: densityMatrix = 8
    integer :: postScc = 9
    integer :: eigvecWriting = 10
    integer :: energyDensityMatrix = 11
    integer :: forceCalc = 12
    integer :: stressCalc = 13
    integer :: postGeoOpt = 14
  end type TGlobalTimersHelper

  type(TGlobalTimersHelper), parameter :: globalTimers = TGlobalTimersHelper()


contains

  !> Returns an initialized instance.
  subroutine TEnvironment_init(this)

    !> Instance
    type(TEnvironment), intent(out) :: this

    call TFileRegistry_init(this%fileFinalizer)

  end subroutine TEnvironment_init


  !> Finalizes the environment.
  subroutine TEnvironment_destruct(this)

    !> Instance
    class(TEnvironment), intent(inout) :: this

    call this%globalTimer%writeTimings()
    call this%fileFinalizer%closeAll()
    flush(stdOut)

  end subroutine TEnvironment_destruct
    
  
  !> Gracefully cleans up and shuts down.
  !>
  !> Note: This routine must be collectively called by all processes.
  !>
  subroutine TEnvironment_shutdown(this)

    !> Instance
    class(TEnvironment), intent(inout) :: this

    call this%destruct()
    call shutdown()

  end subroutine TEnvironment_shutdown


  !> Initlaizes the global timer of the environment.
  subroutine TEnvironment_initGlobalTimer(this, timingLevel, header, unit)

    !> Instance
    class(TEnvironment), intent(inout) :: this

    !> Timing level up to which timings should be printed (-1 for all)
    integer, intent(in) :: timingLevel

    !> Header of the timing table
    character(*), intent(in) :: header

    !> File unit into which the table should be written
    integer, intent(in) :: unit

    call TTimerArray_init(this%globalTimer, globalTimerItems, maxLevel=timingLevel, header=header,&
        & unit=unit)

  end subroutine TEnvironment_initGlobalTimer



#:if WITH_MPI

  !> Initializes MPI environment.
  subroutine TEnvironment_initMpi(this, nGroup)

    !> Instance
    class(TEnvironment), intent(inout) :: this

    !> Number of process groups to create
    integer, intent(in) :: nGroup

    ! MPI settings
    call TMpiEnv_init(this%mpi, nGroup)
    this%tGlobalMaster = this%mpi%tGlobalMaster
    this%nGroup = this%mpi%nGroup
    this%myGroup = this%mpi%myGroup

  end subroutine TEnvironment_initMpi

#:endif


#:if WITH_SCALAPACK

  !> Initializes BLACS environment
  subroutine TEnvironment_initBlacs(this, rowBlock, colBlock, nOrb, nAtom)

    !> Instance
    class(TEnvironment), intent(inout) :: this

    !> Row block size
    integer, intent(in) :: rowBlock

    !> Column block size
    integer, intent(in) :: colBlock

    !> Nr. of orbitals
    integer, intent(in) :: nOrb

    !> Nr. of atoms
    integer, intent(in) :: nAtom

    call TBlacsEnv_init(this%blacs, this%mpi, rowBlock, colBlock, nOrb, nAtom)

  end subroutine TEnvironment_initBlacs

#:endif


end module environment
