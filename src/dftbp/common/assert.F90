!--------------------------------------------------------------------------------------------------!
!  DFTB+: general package for performing fast atomistic simulations                                !
!  Copyright (C) 2006 - 2025  DFTB+ developers group                                               !
!                                                                                                  !
!  See the LICENSE file for terms of usage and distribution.                                       !
!--------------------------------------------------------------------------------------------------!

#:include 'common.fypp'

!> Auxiliary subroutines for the ASSERT command
module dftbp_common_assert
  use dftbp_common_globalenv, only : stdOut0
  use dftbp_common_globalenv, only : abortProgram
  implicit none

  private
#:block DEBUG_CODE
  public :: assertError
#:endblock DEBUG_CODE

contains

#:block DEBUG_CODE


  !> Prints assertion error and abort program execution.
  subroutine assertError(fileName, lineNr, message)

    !> Name of the file in which the error occurred.
    character(*), intent(in) :: fileName

    !> Nr. of the line at which the error occurred.
    integer, intent(in) :: lineNr

    !> Additional message for error
    character(*), intent(in), optional :: message

    write(stdOut0, '(A)') "!!! UNFULLFILLED ASSERTION"
    write(stdOut0, '(A,A)') "!!! FILE:      ", fileName
    write(stdOut0, '(A,I0)') "!!! LINE NR.:  ", lineNr
    if (present(message)) then
      write(stdOut0, '(A,A,A)') '!!! MESSAGE:  "', trim(message), '"'
    end if
    call abortProgram()

  end subroutine assertError

#:endblock DEBUG_CODE

end module dftbp_common_assert
