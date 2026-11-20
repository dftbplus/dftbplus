!--------------------------------------------------------------------------------------------------!
!  DFTB+: general package for performing fast atomistic simulations                                !
!  Copyright (C) 2006 - 2025  DFTB+ developers group                                               !
!                                                                                                  !
!  See the LICENSE file for terms of usage and distribution.                                       !
!--------------------------------------------------------------------------------------------------!

#:include 'common.fypp'

program setupgeom
  use dftbp_common_environment, only : TEnvironment, TEnvironment_init
  use transporttools_inputdata, only : TInputData
  use transporttools_parser, only : parseHsdInput
  use dftbp_common_globalenv, only : destructGlobalEnv, ensureSerialRun, initGlobalEnv
  use dftbp_common_release, only : releaseYear
  use dftbp_io_formatout, only : printDftbHeader
  implicit none (type, external)

  type(TEnvironment) :: env
  type(TInputData), allocatable :: input

  call initGlobalEnv()
  call ensureSerialRun()
  call TEnvironment_init(env)

  call printDftbHeader(env%stdOut, '(setupgeom)', releaseYear)
  allocate(input)
  call parseHsdInput(env, input)
  deallocate(input)
  call env%destruct()
  call destructGlobalEnv()

end program setupgeom
