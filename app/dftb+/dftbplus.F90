!--------------------------------------------------------------------------------------------------!
!  DFTB+: general package for performing fast atomistic simulations                                !
!  Copyright (C) 2006 - 2025  DFTB+ developers group                                               !
!                                                                                                  !
!  See the LICENSE file for terms of usage and distribution.                                       !
!--------------------------------------------------------------------------------------------------!

#:include 'common.fypp'

program dftbplus
  use dftbp_common_environment, only : TEnvironment, TEnvironment_init
  use dftbp_common_exception, only : TEarlyExit, TException, TGeneralError, TInputError,&
      & TInternalError, TModelError, TSystemError
  use dftbp_common_globalenv, only : destructGlobalEnv, initGlobalEnv, stdOut0
  use dftbp_common_release, only : releaseName, releaseYear
  use dftbp_dftbplus_hsdhelpers, only : parseHsdInput
  use dftbp_dftbplus_initprogram, only : TDftbPlusMain
  use dftbp_dftbplus_inputdata, only : TInputData
  use dftbp_dftbplus_main, only : runDftbPlus
  use dftbp_io_formatout, only : printDftbHeader
  implicit none

  type(TEnvironment) :: env
  type(TException), allocatable :: exc
  type(TInputData), allocatable :: input
  type(TDftbPlusMain), allocatable, target :: main
  integer :: exitCode

  call initGlobalEnv()
  call TEnvironment_init(env)

  try: block
    call printDftbHeader(releaseName, releaseYear)

    allocate(input)
    call parseHsdInput(exc, input)
    if (allocated(exc)) exit try

    allocate(main)
    call main%initProgramVariables(exc, env, input)
    deallocate(input)
    if (allocated(exc)) exit try

    call runDftbPlus(main, exc, env)
    call main%destructProgramVariables()
    deallocate(main)
  #:if WITH_MAGMA
    call magmaf_finalize()
  #:endif
  end block try

  call handleException(exc, exitCode)
  call env%destruct(writeTimings=(exitCode == 0))
  call destructGlobalEnv()
  if (exitCode /= 0) error stop exitCode, quiet=.true.

contains

  !> Handle eventual exception
  subroutine handleException(exc, exitCode)

    !> Exception, might be unallocated, if no exception occured.
    type(TException), allocatable, intent(inout) :: exc

    !> Exit code the program should generate
    integer, intent(out) :: exitCode

    logical :: withPropagationPath

    withPropagationPath = .false.
    exitCode = 0
    if (allocated(exc)) then
      select type (cause => exc%cause)
      class is (TInputError)
        write(stdOut0, "(/, a, /)") "*** Input error occured! ***"
        exitCode = 2
      class is (TModelError)
        write(stdOut0, "(/, a, /)") "*** Model error occured! ***"
        exitCode = 3
      class is (TSystemError)
        write(stdOut0, "(/, a, /)") "*** System error occured! ***"
        exitCode = 4
      class is (TInternalError)
        write(stdOut0, "(/, a, /)") "*** Internal error occured (please file a bug report)! ***"
        exitCode = 100
        withPropagationPath = .true.
      ! TGeneralError is superclass of all other errors, must be checked after specific ones
      class is (TGeneralError)
        write(stdOut0, "(/, a, /)") "*** Error occured ***"
        exitCode = 1
      class is (TEarlyExit)
        exitCode = 0
      class default
        write(stdOut0, "(/, a, /)") "*** Unknown error occured (please file a bug report)!"
        exitCode = 200
        withPropagationPath = .true.
      end select
      call exc%writeTo(stdOut0, withPropagationPath=withPropagationPath)
      call exc%deactivate()
    end if

  end subroutine handleException

end program dftbplus
