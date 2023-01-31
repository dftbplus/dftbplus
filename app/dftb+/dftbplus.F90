!--------------------------------------------------------------------------------------------------!
!  DFTB+: general package for performing fast atomistic simulations                                !
!  Copyright (C) 2006 - 2023  DFTB+ developers group                                               !
!                                                                                                  !
!  See the LICENSE file for terms of usage and distribution.                                       !
!--------------------------------------------------------------------------------------------------!

#:include 'common.fypp'

program dftbplus
  use dftbp_common_environment, only : TEnvironment, TEnvironment_init
  use dftbp_common_globalenv, only : initGlobalEnv, destructGlobalEnv
  use dftbp_dftbplus_hsdhelpers, only : parseHsdInput
  use dftbp_dftbplus_initprogram, only : TDftbPlusMain
  use dftbp_dftbplus_inputdata, only : TInputData
  use dftbp_dftbplus_main, only : runDftbPlus
  use dftbp_io_formatout, only : printDftbHeader
  implicit none

  character(len=*), parameter :: releaseName = '${RELEASE}$'
  integer, parameter :: releaseYear = 2022

  type(TEnvironment) :: env
  type(TInputData), allocatable :: input
  type(TDftbPlusMain), allocatable, target :: main

  call initGlobalEnv()
  call printDftbHeader(releaseName, releaseYear)
  allocate(input)
  call parseHsdInput(input)
  call TEnvironment_init(env)
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
