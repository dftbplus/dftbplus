!--------------------------------------------------------------------------------------------------!
!  DFTB+: general package for performing fast atomistic simulations                                !
!  Copyright (C) 2018  DFTB+ developers group                                                      !
!                                                                                                  !
!  See the LICENSE file for terms of usage and distribution.                                       !
!--------------------------------------------------------------------------------------------------!

!> Contains rudimentary functions for warnings and error functions for the code to report problems
!> during run time.
!> Provides routines to call with a string or array of strings if problems occur of a fatal (error)
!> or recoverable (warning) nature.
module message
  use globalenv
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

  public :: warning, error

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
    call synchronizeAll()
    call abort()

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
    call synchronizeAll()
    call abort()

  end subroutine error_array

end module message
