!--------------------------------------------------------------------------------------------------!
!  DFTB+: general package for performing fast atomistic simulations                                !
!  Copyright (C) 2018  DFTB+ developers group                                                      !
!                                                                                                  !
!  See the LICENSE file for terms of usage and distribution.                                       !
!--------------------------------------------------------------------------------------------------!

#:include 'common.fypp'

program dftbplus
  use globalenv
  use environment
  use main, only : runDftbPlus
  use inputdata_module, only : inputData
  use formatout, only : printDftbHeader
  use parser, only : parseHsdInput
  use initprogram, only : initProgramVariables
  implicit none

  character(len=*), parameter :: releaseName = '${RELEASE}$'
  integer, parameter :: releaseYear = 2018

  type(TEnvironment) :: env
  type(inputData), allocatable :: input

  call initGlobalEnv()
  call printDftbHeader(releaseName, releaseYear)
  allocate(input)
  call parseHsdInput(input)
  call TEnvironment_init(env)
  call initProgramVariables(input, env)
  deallocate(input)
  call runDftbPlus(env)
  call env%destruct()
  call destructGlobalEnv()

end program dftbplus
