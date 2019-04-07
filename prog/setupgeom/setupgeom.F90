!--------------------------------------------------------------------------------------------------!
!  DFTB+: general package for performing fast atomistic simulations                                !
!  Copyright (C) 2018  DFTB+ developers group                                                      !
!                                                                                                  !
!  See the LICENSE file for terms of usage and distribution.                                       !
!--------------------------------------------------------------------------------------------------!

#:include 'common.fypp'

program setupGeometry 
  use globalenv
  use environment
  use inputdata_module, only : inputData
  use formatout, only : printDftbHeader
  use parser_setup, only : parseHsdInput
  implicit none

  character(len=*), parameter :: releaseName = ''
  integer, parameter :: releaseYear = 2018

  type(TEnvironment) :: env
  type(inputData), allocatable :: input

  call initGlobalEnv()
  call printDftbHeader(releaseName, releaseYear)
  allocate(input)
  call parseHsdInput(input)
  deallocate(input)
  call env%destruct()
  call destructGlobalEnv()

end program setupGeometry
