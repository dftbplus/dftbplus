!--------------------------------------------------------------------------------------------------!
!  DFTB+: general package for performing fast atomistic simulations                                !
!  Copyright (C) 2017  DFTB+ developers group                                                      !
!                                                                                                  !
!  See the LICENSE file for terms of usage and distribution.                                       !
!--------------------------------------------------------------------------------------------------!

!!* Contains rudimentary warning and error functions for the code to report
!!* problems during run time.
!!* @desc Provides routines to call with a string or array of strings if
!!*   problems occur of a fatal (error) or recoverable (warning) nature.
module message
  implicit none
  private

  interface warning
    module procedure warning_single
    module procedure warning_array
  end interface

  interface error
    module procedure error_single
    module procedure error_array
  end interface

  public :: warning, error

contains

  !!* Gives a warning message.
  !!* @param message Warning message to print to standard out.
  subroutine warning_single(message)
    character (len=*), intent(in) :: message
    write(*, '(1a)') 'WARNING!'
    write(*, '(2a)') '-> ', trim(message)
  end subroutine warning_single



  !!* Gives a warning message.
  !!* @param messages Lines of the error message to print to standard out.
  subroutine warning_array(messages)
    character(len=*), intent(in) :: messages(:)

    integer :: ii

    write(*,'(1a)') 'WARNING!'
    do ii = 1, size(messages)
      write(*,'(2a)') '-> ', trim(messages(ii))
    end do
    stop

  end subroutine warning_array



  !!* Gives an error message and stops the code.
  !!* @param message Error message to print to standard out.
  subroutine error_single(message)
    character (len=*), intent(in) :: message
    write(*,'(1a)') 'ERROR!'
    write(*,'(2a)') '-> ', trim(message)
    stop
  end subroutine error_single



  !!* Gives an error messages and stops the code.
  !!* @param messages Lines of the error message to print to standard out.
  subroutine error_array(messages)
    character(len=*), intent(in) :: messages(:)

    integer :: ii

    write(*,'(1a)') 'ERROR!'
    do ii = 1, size(messages)
      write(*,'(2a)') '-> ', trim(messages(ii))
    end do
    stop

  end subroutine error_array

end module message


