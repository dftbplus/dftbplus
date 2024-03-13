!--------------------------------------------------------------------------------------------------!
!  DFTB+: general package for performing fast atomistic simulations                                !
!  Copyright (C) 2006 - 2024  DFTB+ developers group                                               !
!                                                                                                  !
!  See the LICENSE file for terms of usage and distribution.                                       !
!--------------------------------------------------------------------------------------------------!

#:include 'common.fypp'

!> Contains fundamental global computing environment settings.
!>
!> It contains global settings and routines, which can be used already before the parsing of the
!> input has taken place and the details of the user settings for the running-environment are
!> known. Also, it can be used by routines which are not MPI-aware but wish to make I/O or abort the
!> code.
!>
module dftbp_common_globalenv
  use, intrinsic :: iso_fortran_env, only : output_unit, error_unit
#:if WITH_MPI
  use mpi, only : MPI_COMM_WORLD
  use dftbp_extlibs_mpifx, only : mpifx_comm, MPI_THREAD_FUNNELED, mpifx_init_thread, mpifx_abort,&
      & mpifx_barrier, mpifx_finalize
#:endif
  implicit none

  private
  public :: initGlobalEnv, destructGlobalEnv
  public :: abortProgram, shutdown, synchronizeAll
  public :: stdOut, stdOut0, stdErr, stdErr0, tIoProc
  public :: withScalapack, withMpi
  public :: instanceSafeBuild
  #:if WITH_MPI
    public :: globalMpiComm
  #:endif


  !> Unredirected standard out
  integer, parameter :: stdOut0 = output_unit

  !> Unredirected standard error
  integer, parameter :: stdErr0 = error_unit

  !> Standard out file handler
  integer, protected :: stdOut = stdOut0

  !> Standard error file handler
  integer, protected :: stdErr = stdErr0

  !> Whether current process is the global lead process
  logical, protected :: tIoProc = .true.

#:if WITH_MPI
  !> Global MPI communicator (used for aborts)
  type(mpifx_comm), protected :: globalMpiComm
#:endif

  !> Whether code was compiled with MPI support
  logical, parameter :: withMpi = ${FORTRAN_LOGICAL(WITH_MPI)}$

  !> Whether code was compiled with Scalapack
  logical, parameter :: withScalapack = ${FORTRAN_LOGICAL(WITH_SCALAPACK)}$

#:if WITH_MPI
  !> Whether MPI finalization should be performed at the end
  logical :: doMpiFinalization = .true.
#:endif

  !> Whether code was compiled with many-body dispersion support
  logical, parameter :: withMbd = ${FORTRAN_LOGICAL(WITH_MBD)}$

  !> Whether the code had been built with instance safe components only
  logical, parameter :: instanceSafeBuild = ${FORTRAN_LOGICAL(INSTANCE_SAFE_BUILD)}$



contains

  !> Initializes global environment (must be the first statement of a program)
  subroutine initGlobalEnv(outputUnit, mpiComm, errorUnit, devNull)

    !> Customised global standard output
    integer, intent(in), optional :: outputUnit

    !> Customised global MPI communicator
    integer, intent(in), optional :: mpiComm

    !> Customised global standard error
    integer, intent(in), optional :: errorUnit

    !> Unit of the null device (needed for follow processes to suppress their output)
    integer, intent(in), optional :: devNull


    integer :: outputUnit0, errorUnit0, devNull0

  #:if WITH_MPI
    integer :: mpiComm0
  #:endif

    if (present(outputUnit)) then
      outputUnit0 = outputUnit
    else
      outputUnit0 = stdOut0
    end if

    if (present(errorUnit)) then
      errorUnit0 = errorUnit
    else
      errorUnit0 = stdErr0
    end if

  #:if WITH_MPI
    if (present(mpiComm)) then
      mpiComm0 = mpiComm
      doMpiFinalization = .false.
    else
      mpiComm0 = MPI_COMM_WORLD
      call mpifx_init_thread(requiredThreading=MPI_THREAD_FUNNELED)
    end if

    call globalMpiComm%init(commid=mpiComm0)
    if (globalMpiComm%lead) then
      stdOut = outputUnit0
      stdErr = errorUnit0
    else
      if (present(devNull)) then
        devNull0 = devNull
      else
        open(newunit=devNull0, file="/dev/null", action="write")
      end if
      stdOut = devNull0
      stdErr = devNull0
    end if
    tIoProc = globalMpiComm%lead
  #:else
    stdOut = outputUnit0
    stdErr = errorUnit0
  #:endif

  end subroutine initGlobalEnv


  !> Finalizes global environment (must be the last statement of a program)
  subroutine destructGlobalEnv()

  #:if WITH_MPI
    if (doMpiFinalization) then
      call mpifx_finalize()
    end if
  #:endif

  end subroutine destructGlobalEnv


  !> Gracefully shuts down environment and stops execution.
  !>
  !> Note: this routine must be called collectively by all processes.
  !>
  subroutine shutdown()

    call synchronizeAll()
    call destructGlobalEnv()
    stop

  end subroutine shutdown


  !> Aborts program execution immediately.
  !>
  !> Note: if this routine is called by any the processes, execution immediately stops
  !> without waiting for any other processes.
  !>
  subroutine abortProgram(errorCode)

    !> Error code to emit (default: 1)
    integer, intent(in), optional :: errorCode

    integer :: errorCode0
  #:if WITH_MPI
    integer :: error
  #:endif

    if (.not. present(errorCode)) then
      errorCode0 = 1
    else
      errorCode0 = errorCode
    end if

  #:if WITH_MPI
    call mpifx_abort(globalMpiComm, errorCode0, error)
    if (error /= 0) then
      write(stdErr0, "(A,I0,A)") "Process ", globalMpiComm%rank, " could not be aborted."
    end if
  #:endif
    error stop

  end subroutine abortProgram


  !> Waits until all processes reach this point
  subroutine synchronizeAll()

  #:if WITH_MPI
    call mpifx_barrier(globalMpiComm)
  #:endif

  end subroutine synchronizeAll


end module dftbp_common_globalenv
