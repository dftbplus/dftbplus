!--------------------------------------------------------------------------------------------------!
!  DFTB+: general package for performing fast atomistic simulations                                !
!  Copyright (C) 2018  DFTB+ developers group                                                      !
!                                                                                                  !
!  See the LICENSE file for terms of usage and distribution.                                       !
!--------------------------------------------------------------------------------------------------!

module io
  use, intrinsic :: iso_fortran_env
  implicit none
  private

  public :: stdin, stdout, stderr


  !> File id of standard input
  integer, parameter :: stdIn = input_unit

  !> File id of standard output
  integer, parameter :: stdOut = output_unit

  !> File id of standard error
  integer, parameter :: stdErr = error_unit

end module io
