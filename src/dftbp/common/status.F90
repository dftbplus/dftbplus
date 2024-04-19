!--------------------------------------------------------------------------------------------------!
!  DFTB+: general package for performing fast atomistic simulations                                !
!  Copyright (C) 2006 - 2023  DFTB+ developers group                                               !
!                                                                                                  !
!  See the LICENSE file for terms of usage and distribution.                                       !
!--------------------------------------------------------------------------------------------------!

#:include 'common.fypp'

!> Implements a status object to indicate possible errors.
module dftbp_common_status
  implicit none

  private
  public :: TStatus


  !> Represents status of an operation.
  type :: TStatus

    !> Status code. Zero signals no error, every other value is interpreted as an error.
    integer :: code = 0

    !> Error message (in case error code is non-zero)
    character(:), allocatable :: message

    !> Source file where error was generated (in case error code is non-zero)
    character(:), allocatable :: file

    !> Number of the line where error was generated (in case error code is non-zero)
    integer :: line = 0

  contains

    procedure :: setError
    procedure :: hasError
    procedure :: isOk

  end type TStatus


contains

  !> Sets the status to error.
  subroutine setError(this, code, message, file, line)

    !> Instance.
    class(TStatus), intent(inout) :: this

    !> Error code.
    integer, intent(in) :: code

    !> Error message.
    character(*), intent(in) :: message

    !> File where error occurred.
    character(*), intent(in) :: file

    !> Line where error occurred.
    integer, intent(in) :: line

    this%code = code
    this%message = message
    this%file = file
    this%line = line

  end subroutine setError


  !> Whether status contains an error.
  function hasError(this)

    !> Instance
    class(TStatus), intent(in) :: this

    !> True if an error occurred
    logical :: hasError

    hasError = this%code /= 0

  end function hasError


  !> Whether status is OK (contains no error).
  function isOk(this)

    !> Instance
    class(TStatus), intent(in) :: this

    !> Whether status is error-free
    logical :: isOk

    isOk = .not. this%hasError()

  end function isOk


end module dftbp_common_status
