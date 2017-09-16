#:include 'common.fypp'

module environment
  use, intrinsic :: iso_fortran_env
#:if WITH_MPI
  use mpifx
  use scalapackfx
#:endif
  implicit none
  private

  public :: initializeGlobalEnv, finalizeGlobalEnv
  public :: TEnvironment, initialize
  public :: abort
  public :: stdOut, stdOut0

  !> Standard out file handler
  integer, protected :: stdOut

  !> Unredirected standard out
  integer, protected :: stdOut0

#:if WITH_MPI
  !> Global MPI communicator (used for aborts)
  type(mpifx_comm) :: mpiCommGlobal
#:endif

  !> Contains environment settings.
  type :: TEnvironment
    private

    !> Whether this process is the master?
    logical, public :: tMaster

    !> Whether this process is supposed to do I/O
    logical, public :: tIO

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

  end type TEnvironment


  !> Initializes environment instance
  interface initialize
    module procedure Environment_initialize
  end interface initialize


contains

  !> Initializes global environment (must be the first statement of a program)
  subroutine initializeGlobalEnv()

  #:if WITH_MPI
    call mpifx_init()
    call mpiCommGlobal%init()
    stdOut0 = output_unit
    if (mpiCommGlobal%master) then
      stdOut = stdOut0
    else
      stdOut = 1
      open(stdOut, file="/dev/null", action="write")
    end if
  #:endif

  end subroutine initializeGlobalEnv


  !> Finalizes global environment (must be the last statement of a program)
  subroutine finalizeGlobalEnv()

  #:if WITH_MPI
    call mpifx_finalize()
  #:endif
    
  end subroutine finalizeGlobalEnv


  !> Initializes environment instance.
  subroutine Environment_initialize(this)

    !> Instance
    type(TEnvironment), intent(out) :: this

  #:if WITH_MPI
  #! MPI settings
    call this%mpiComm%init()
    this%ioProcId = this%mpiComm%masterrank
    this%tMaster = this%mpiComm%master
    this%tIO = this%tMaster
    this%stdOut = stdOut
  #:else
  #! NON-MPI settings
    this%tMaster = .true.
    this%tIO = .true.
    this%stdOut = stdOut0
  #:endif

  end subroutine Environment_initialize


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
    call mpifx_abort(mpiCommGlobal, errorCode0, error)
    if (error /= 0) then
      write(stdOut0, "(A,I0,A)") "Process ", mpiCommGlobal%rank, " could not be aborted."
    end if
  #:endif
    error stop

  end subroutine abort

end module environment
