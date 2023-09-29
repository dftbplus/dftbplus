!--------------------------------------------------------------------------------------------------!
!  DFTB+: general package for performing fast atomistic simulations                                !
!  Copyright (C) 2006 - 2025  DFTB+ developers group                                               !
!                                                                                                  !
!  See the LICENSE file for terms of usage and distribution.                                       !
!--------------------------------------------------------------------------------------------------!

#:include 'common.fypp'

program dftbplus
  use dftbp_common_environment, only : TEnvironment, TEnvironment_init
  use dftbp_common_globalenv, only : initGlobalEnv, destructGlobalEnv, stdOut
  use dftbp_common_release, only : releaseName, releaseYear
  use dftbp_dftbplus_hsdhelpers, only : parseHsdInput
  use dftbp_dftbplus_initprogram, only : TDftbPlusMain
  use dftbp_dftbplus_inputdata, only : TInputData
  use dftbp_dftbplus_main, only : runDftbPlus
  use dftbp_io_formatout, only : printDftbHeader
  implicit none

  type(TEnvironment) :: env
  type(TInputData), allocatable :: input
  type(TDftbPlusMain), allocatable, target :: main

  call initGlobalEnv()
  allocate(input)
  call TEnvironment_init(env)
  ! temporary fix
  env%stdOut = stdOut
  call printDftbHeader(env%stdOut, releaseName, releaseYear)
  call parseHsdInput(env, input)
  allocate(main)
  call main%initProgramVariables(input, env)
  deallocate(input)
  call runDftbPlus(main, env)
  call main%destructProgramVariables()
  deallocate(main)
#:if WITH_MAGMA
  call magmaf_finalize()
#:endif
  call env%destruct()
  call destructGlobalEnv()

end program dftbplus
