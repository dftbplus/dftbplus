#:include 'common.fypp'

module environment
  use, intrinsic :: iso_fortran_env
#:if WITH_MPI
  use mpifx
  use scalapackfx
#:endif
  implicit none
  private

  public :: initGlobalEnv, destructGlobalEnv
  public :: TEnvironment, init
  public :: abort
  public :: stdOut, stdOut0
  public :: includesAllProcesses

  !> Standard out file handler
  integer, protected :: stdOut

  !> Unredirected standard out
  integer, protected :: stdOut0

#:if WITH_MPI
  !> Global MPI communicator (used for aborts)
  type(mpifx_comm) :: globalMpiComm
#:endif

  !> Contains environment settings.
  type :: TEnvironment
    private

    !> Whether this process is the master?
    logical, public :: tMaster

    !> Whether this process is supposed to do I/O
    logical, public :: tIoProc

    !> Standard output to use for messages (redirected to /dev/null for non-IO processes) 
    integer, public :: stdOut

  #:if WITH_MPI

    !> Rank of the process designated for I/O.
    integer, public :: ioProcId

    !> MPI communicator
    type(mpifx_comm), public :: mpiComm

    !> BLACS communicator
    type(blacsgrid), public :: blacsComm

  #:endif

  contains
    procedure :: includesAllProcesses
  end type TEnvironment


  !> Initializes environment instance
  interface init
    module procedure Environment_init
  end interface init


contains

  !> Initializes global environment (must be the first statement of a program)
  subroutine initGlobalEnv()

  #:if WITH_MPI
    call mpifx_init()
    call globalMpiComm%init()
    stdOut0 = output_unit
    if (globalMpiComm%master) then
      stdOut = stdOut0
    else
      stdOut = 1
      open(stdOut, file="/dev/null", action="write")
    end if
  #:endif

  end subroutine initGlobalEnv


  !> Finalizes global environment (must be the last statement of a program)
  subroutine destructGlobalEnv()

  #:if WITH_MPI
    call mpifx_finalize()
  #:endif
    
  end subroutine destructGlobalEnv


  !> Initializes environment instance.
  subroutine Environment_init(this)

    !> Instance
    type(TEnvironment), intent(out) :: this

  #:if WITH_MPI
  #! MPI settings
    call this%mpiComm%init()
    this%ioProcId = this%mpiComm%masterrank
    this%tMaster = this%mpiComm%master
    this%tIoProc = this%tMaster
    this%stdOut = stdOut
  #:else
  #! NON-MPI settings
    this%tMaster = .true.
    this%tIoProc = .true.
    this%stdOut = stdOut0
  #:endif

  end subroutine Environment_init


  !> Aborts program execution.
  subroutine abort(errorCode)

    !> Error code to emit (default: 1)
    integer, intent(in), optional :: errorCode

    integer :: error, errorCode0

    if (.not. present(errorCode)) then
      errorCode0 = 1
    else
      errorCode0 = errorCode
    end if

  #:if WITH_MPI
    call mpifx_abort(globalMpiComm, errorCode0, error)
    if (error /= 0) then
      write(stdOut0, "(A,I0,A)") "Process ", globalMpiComm%rank, " could not be aborted."
    end if
  #:endif
    error stop

  end subroutine abort


  function includesAllProcesses(this) result(tIncludesAll)
    class(TEnvironment), intent(in) :: this
    logical :: tIncludesAll

  #:if WITH_MPI
    tIncludesAll = (this%mpiComm%id == globalMpiComm%id)
  #:else
    tIncludesAll = .true.
  #:endif
    
  end function includesAllProcesses

end module environment
