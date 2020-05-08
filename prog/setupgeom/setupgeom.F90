!--------------------------------------------------------------------------------------------------!
!  DFTB+: general package for performing fast atomistic simulations                                !
!  Copyright (C) 2006 - 2020  DFTB+ developers group                                               !
!                                                                                                  !
!  See the LICENSE file for terms of usage and distribution.                                       !
!--------------------------------------------------------------------------------------------------!

#:include 'common.fypp'

program setupGeometry 
  use dftbp_globalenv
  use dftbp_formatout, only : printDftbHeader
  use dftbp_inputsetup, only : TInputData
  use dftbp_parsersetup, only : parseHsdInput
  implicit none

  character(len=*), parameter :: releaseName = ''
  integer, parameter :: releaseYear = 2020

  type(TInputData), allocatable :: input

  call initGlobalEnv()
  call printDftbHeader(releaseName, releaseYear)
  allocate(input)
  call parseHsdInput(input)
  deallocate(input)
  call destructGlobalEnv()

end program setupGeometry
