!--------------------------------------------------------------------------------------------------!
!  DFTB+: general package for performing fast atomistic simulations                                !
!  Copyright (C) 2018  DFTB+ developers group                                                      !
!                                                                                                  !
!  See the LICENSE file for terms of usage and distribution.                                       !
!--------------------------------------------------------------------------------------------------!

#:include 'common.fypp'

!> Auxiliary subroutines for the ASSERT command
module assert
  use io
  implicit none
  private

#:call ASSERT_CODE
  public :: assertError
#:endcall ASSERT_CODE

contains

#:call ASSERT_CODE


  !> Prints assertion error and abort program execution.
  subroutine assertError(fileName, lineNr)

    !> Name of the file in which the error occured.
    character(*), intent(in) :: fileName

    !> Nr. of the line at which the error occured.
    integer, intent(in) :: lineNr

    write(stdout, '(A)') "!!! UNFULLFILLED ASSERTION"
    write(stdout, '(A,A)') "!!! FILE:      ", fileName
    write(stdout, '(A,I0)') "!!! LINE NR.:  ", lineNr
    stop

  end subroutine assertError

#:endcall ASSERT_CODE

end module assert
