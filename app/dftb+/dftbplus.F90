!--------------------------------------------------------------------------------------------------!
!  DFTB+: general package for performing fast atomistic simulations                                !
!  Copyright (C) 2006 - 2025  DFTB+ developers group                                               !
!                                                                                                  !
!  See the LICENSE file for terms of usage and distribution.                                       !
!--------------------------------------------------------------------------------------------------!

#:include 'common.fypp'

program dftbplus
  use dftbp_common_environment, only : TEnvironment, TEnvironment_init
  use dftbp_common_exception, only : TEarlyExit, TException, TException_destroy, TFatalError,&
      & TInternalError
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
  logical :: hasError

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

  ! Handle eventual exception
  hasError = allocated(exc)
  if (hasError) then
    select type (event => exc%event)
    class is (TFatalError)
      write(stdOut0, "(/, a, /)") "*** Fatal error occured! ***"
    class is (TInternalError)
      write(stdOut0, "(/, a, /)") "*** Internal error occured (please file bug report)! ***"
    class is (TEarlyExit)
      hasError = .false.
    class default
      write(stdOut0, "(/, a, /)") "*** Unknown error occured! ***"
    end select
    call exc%writeTo(stdOut0, withPropagationPath=hasError)
    call TException_destroy(exc)
  end if

  call env%destruct(writeTimings=(.not. hasError))
  call destructGlobalEnv()

  if (hasError) error stop

end program dftbplus
