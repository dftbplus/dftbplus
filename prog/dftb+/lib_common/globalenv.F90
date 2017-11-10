!--------------------------------------------------------------------------------------------------!
!  DFTB+: general package for performing fast atomistic simulations                                !
!  Copyright (C) 2017  DFTB+ developers group                                                      !
!                                                                                                  !
!  See the LICENSE file for terms of usage and distribution.                                       !
!--------------------------------------------------------------------------------------------------!

#:include 'common.fypp'

!> Contains fundamental global environment settings
module globalenv
  use, intrinsic :: iso_fortran_env, only : output_unit
#:if WITH_MPI
  use mpifx
#:endif
  implicit none
  private

  public :: initGlobalEnv, destructGlobalEnv
  public :: abort, synchronizeAll
  public :: stdOut

  !> Standard out file handler
  integer, protected :: stdOut

#:if WITH_MPI
  !> Global MPI communicator (used for aborts)
  type(mpifx_comm) :: globalMpiComm
#:endif

  !> Unredirected standard out
  integer, parameter :: stdOut0 = output_unit


contains

  !> Initializes global environment (must be the first statement of a program)
  subroutine initGlobalEnv()

  #:if WITH_MPI
    call mpifx_init()
    call globalMpiComm%init()
    if (globalMpiComm%master) then
      stdOut = stdOut0
    else
      stdOut = 1
      open(stdOut, file="/dev/null", action="write")
    end if
  #:else
    stdOut = stdOut0
  #:endif

  end subroutine initGlobalEnv


  !> Finalizes global environment (must be the last statement of a program)
  subroutine destructGlobalEnv()

  #:if WITH_MPI
    call mpifx_finalize()
  #:endif

  end subroutine destructGlobalEnv


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


  !> Waits until all processes reach this point
  subroutine synchronizeAll()

  #:if WITH_MPI
    call mpifx_barrier(globalMpiComm)
  #:endif

  end subroutine synchronizeAll


end module globalenv
