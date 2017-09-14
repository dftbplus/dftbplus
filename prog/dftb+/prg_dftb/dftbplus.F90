!--------------------------------------------------------------------------------------------------!
!  DFTB+: general package for performing fast atomistic simulations                                !
!  Copyright (C) 2017  DFTB+ developers group                                                      !
!                                                                                                  !
!  See the LICENSE file for terms of usage and distribution.                                       !
!--------------------------------------------------------------------------------------------------!

#:include 'common.fypp'

program dftbplus
  use environment
  use main, only : runDftbPlus
  use inputdata_module, only : inputData
  use formatout, only : printDftbHeader
  use parser, only : parseHsdInput
  use initprogram, only : initProgramVariables
  implicit none

  character(len=*), parameter :: RELEASE_VERSION = '17.1'
  integer, parameter :: RELEASE_YEAR = 2017

  type(TEnvironment) :: env
  type(inputData), allocatable :: input

  call initialize(env)
  call printDftbHeader(env, RELEASE_VERSION, RELEASE_YEAR)
#:if WITH_MPI
  ! Following parts not MPI-aware -> stop here
  call env%abort(1)
#:endif
  allocate(input)
  call parseHsdInput(input)
  call initProgramVariables(input)
  deallocate(input)
  call runDftbPlus()
  call finalize(env)

end program dftbplus
