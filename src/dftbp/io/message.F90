!--------------------------------------------------------------------------------------------------!
!  DFTB+: general package for performing fast atomistic simulations                                !
!  Copyright (C) 2006 - 2023  DFTB+ developers group                                               !
!                                                                                                  !
!  See the LICENSE file for terms of usage and distribution.                                       !
!--------------------------------------------------------------------------------------------------!

!> Contains rudimentary functions for warnings and error functions for the code to report problems
!> during run time.
!> Provides routines to call with a string or array of strings if problems occur of a fatal (error)
!> or recoverable (warning) nature.
module dftbp_io_message
  use dftbp_common_globalenv, only : stdOut, synchronizeAll, abortProgram
  implicit none

  private


  !> recoverable error warnings
  interface warning
    module procedure warning_single
    module procedure warning_array
  end interface


  !> fatal error warnings, terminating the code
  interface error
    module procedure error_single
    module procedure error_array
  end interface


  !> terminating the code with a message
  interface cleanShutdown
    module procedure shutdown_single
    module procedure shutdown_array
  end interface


  public :: warning, error, cleanShutdown, getJointMessage

contains


  !> Gives a warning message.
  subroutine warning_single(message)

    !> Warning message to print to standard out.
    character (len=*), intent(in) :: message

    write(stdOut, '(1a)') 'WARNING!'
    write(stdOut, '(2a)') '-> ', trim(message)

  end subroutine warning_single


  !> Gives a warning message.
  subroutine warning_array(messages)

    !> Lines of the error message to print to standard out.
    character(len=*), intent(in) :: messages(:)

    integer :: ii

    write(stdOut, '(1a)') 'WARNING!'
    do ii = 1, size(messages)
      write(stdOut, '(2a)') '-> ', trim(messages(ii))
    end do

  end subroutine warning_array


  !> Gives an error message and stops the code.
  subroutine error_single(message)

    !> Error message to print to standard out.
    character (len=*), intent(in) :: message

    write(stdOut, '(1a)') 'ERROR!'
    write(stdOut, '(2a)') '-> ', trim(message)
    flush(stdOut)
    call abortProgram()

  end subroutine error_single


  !> Gives an error messages and stops the code.
  subroutine error_array(messages)

    !> Lines of the error message to print to standard out.
    character(len=*), intent(in) :: messages(:)

    integer :: ii

    write(stdOut, '(1a)') 'ERROR!'
    do ii = 1, size(messages)
      write(stdOut, '(2a)') '-> ', trim(messages(ii))
    end do
    flush(stdOut)
    call abortProgram()

  end subroutine error_array


  !> Joins an array of error messages.
  subroutine getJointMessage(messages, message)

    !> Lines of the error message to print to standard out
    character(len=*), intent(in) :: messages(:)

    !> Joint message
    character(len=:), intent(out), allocatable :: message

    integer, parameter :: lenNewline = 1
    integer :: totLength, iMessage, iAppend

    if (size(messages, dim=1) >= 1) then

      totLength = 0

      do iMessage = 1, size(messages, dim=1)
        totLength = totLength + len(messages(iMessage))
      end do
      totLength = totLength + (size(messages, dim=1) - 1) * lenNewline

      allocate(character(len=totLength) :: message)
      message(:len(messages(1))) = messages(1)

      iAppend = len(messages(1)) + 1
      do iMessage = 2, size(messages, dim=1)
        message(iAppend:len(messages(iMessage)) + 1) = messages(iMessage) // NEW_LINE("A")
        iAppend = len(messages(iMessage)) + 2
      end do

    end if

  end subroutine getJointMessage


  !> Prints a message and stops the code cleanly.
  subroutine shutdown_single(message)

    !> Shutdown message to print to standard out.
    character (len=*), intent(in) :: message

    write(stdOut, '(A)') trim(message)
    flush(stdOut)
    call synchronizeAll()
    call abortProgram()

  end subroutine shutdown_single


  !> Prints messages and stops the code cleanly.
  subroutine shutdown_array(messages)

    !> Lines of the shutdown message to print to standard out.
    character(len=*), intent(in) :: messages(:)

    integer :: ii

    do ii = 1, size(messages)
      write(stdOut, '(A)') trim(messages(ii))
    end do
    flush(stdOut)
    call synchronizeAll()
    call abortProgram()

  end subroutine shutdown_array

end module dftbp_io_message
