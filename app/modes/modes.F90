!--------------------------------------------------------------------------------------------------!
!  DFTB+: general package for performing fast atomistic simulations                                !
!  Copyright (C) 2006 - 2025  DFTB+ developers group                                               !
!                                                                                                  !
!  See the LICENSE file for terms of usage and distribution.                                       !
!--------------------------------------------------------------------------------------------------!

#:include 'common.fypp'

program modes
  use dftbp_common_environment, only : TEnvironment, TEnvironment_init
  use dftbp_common_globalenv, only : initGlobalEnv, destructGlobalEnv
  use dftbp_common_release, only : releaseYear
  use dftbp_io_formatout, only : printDftbHeader
  use modes_hsdhelpers, only : parseHsdInput
  use modes_initmodes, only : TModesMain
  use modes_inputdata, only : TInputData
  use modes_main, only : runModes
  implicit none

  !> Program version
  character(len=*), parameter :: version = "0.03"

  type(TEnvironment) :: env
  type(TInputData), allocatable :: input
  type(TModesMain), allocatable, target :: main

  call initGlobalEnv()
  call printDftbHeader('(MODES '// version //')', releaseYear)
  allocate(input)
  call parseHsdInput(input)
  call TEnvironment_init(env)
  allocate(main)
  call main%initProgramVariables(input, env)
  deallocate(input)
  call runModes(main, env)
  ! call main%destructProgramVariables()
  deallocate(main)
#:if WITH_MAGMA
  call magmaf_finalize()
#:endif
  call env%destruct()
  call destructGlobalEnv()

end program modes
