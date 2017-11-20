!--------------------------------------------------------------------------------------------------!
!  DFTB+: general package for performing fast atomistic simulations                                !
!  Copyright (C) 2017  DFTB+ developers group                                                      !
!                                                                                                  !
!  See the LICENSE file for terms of usage and distribution.                                       !
!--------------------------------------------------------------------------------------------------!

#:include 'common.fypp'

!> Contains computer environment settings
module environment
  use timerarray
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

    !> Whether this process is the master?
    logical, public :: tGlobalMaster = .true.

    !> Nr. of groups in the system
    integer, public :: nGroup = 1

    !> Id of current group (starts with 0)
    integer, public :: myGroup = 0

    !> Global timers
    type(TTimerArray) :: globalTimer

  #:if WITH_MPI
    !> Global mpi settings
    type(TMpiEnv), public :: mpi
  #:endif
  #:if WITH_SCALAPACK
    !> Global scalapack settings
    type(TBlacsEnv), public :: blacs
  #:endif

  contains
  #:if WITH_MPI
    procedure :: initMpi
  #:endif
  #:if WITH_SCALAPACK
    procedure :: initBlacs
  #:endif

  end type TEnvironment

  type(TTimerItem), parameter :: globalTimerItems(*) = [&
      & TTimerItem("Global initialisation", 1),&
      & TTimerItem("Pre-SCC initialisation", 1),&
      & TTimerItem("SCC", 1),&
      & TTimerItem("Diagonalisation", 2),&
      & TTimerItem("Density matrix creation", 2),&
      & TTimerItem("Post-SCC processing", 1),&
      & TTimerItem("Eigenvector writing", 2),&
      & TTimerItem("Force calculation", 2),&
      & TTimerItem("Post-geometry optimisation", 1)]

  type :: TGlobalTimersHelper
    integer :: globalInit = 1
    integer :: preSccInit = 2
    integer :: scc = 3
    integer :: diagonalization = 4
    integer :: densityMatrix = 5
    integer :: postScc = 6
    integer :: eigvecWriting = 7
    integer :: forceCalc = 8
    integer :: postGeoOpt = 9
  end type TGlobalTimersHelper

  type(TGlobalTimersHelper), parameter :: globalTimers = TGlobalTimersHelper()


contains

  !> Returns an initialized instance.
  subroutine TEnvironment_init(this)

    !> Instance
    type(TEnvironment), intent(out) :: this

    call TTimerArray_init(this%globalTimer, globalTimerItems)

  end subroutine TEnvironment_init


#:if WITH_MPI

  !> Initializes MPI environment.
  subroutine initMpi(this, nGroup)

    !> Instance
    class(TEnvironment), intent(inout) :: this

    !> Number of process groups to create
    integer, intent(in) :: nGroup

    ! MPI settings
    call TMpiEnv_init(this%mpi, nGroup)
    this%tGlobalMaster = this%mpi%tGlobalMaster
    this%nGroup = this%mpi%nGroup
    this%myGroup = this%mpi%myGroup

  end subroutine initMpi

#:endif


#:if WITH_SCALAPACK

  !> Initializes BLACS environment
  subroutine initBlacs(this, rowBlock, colBlock, nOrb, nAtom)

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

  end subroutine initBlacs

#:endif


end module environment
