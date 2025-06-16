!--------------------------------------------------------------------------------------------------!
!  DFTB+: general package for performing fast atomistic simulations                                !
!  Copyright (C) 2006 - 2025  DFTB+ developers group                                               !
!                                                                                                  !
!  See the LICENSE file for terms of usage and distribution.                                       !
!--------------------------------------------------------------------------------------------------!

#:include 'common.fypp'

!> Implements a status object to indicate possible errors.
module dftbp_common_exception
  use, intrinsic :: iso_fortran_env, only : stdErr => error_unit
  implicit none

  private
  public :: TException, TException_create
  public :: TExceptionCause
  public :: TEarlyExit, TGeneralError, TInputError, TModelError, TSystemError, TInternalError


  !> Represents an exception propagation path item
  type :: TPathItem
    character(:), allocatable :: file
    integer :: line = 0
    type(TPathItem), allocatable :: next
  end type TPathItem


  !> Represents the propagation path obtained during exception propagation
  type :: TPropagationPath
    type(TPathItem), allocatable :: lastItem
  contains
    procedure :: addItem => TPropagationPath_addItem
    procedure :: writeTo => TPropagationPath_writeTo
  end type TPropagationPath


  !> Represents the details about why the exception was triggered
  type, abstract :: TExceptionCause
  contains
    procedure(TExceptionCause_asChar), deferred :: asChar
  end type TExceptionCause


  abstract interface
    !> Character representation of the exception details
    function TExceptionCause_asChar(this) result(str)
      import :: TExceptionCause
      implicit none

      !> Instance
      class(TExceptionCause), intent(in) :: this

      !> Character representation of the cause
      character(:), allocatable :: str

    end function TExceptionCause_asChar
  end interface


  !> Represents an exception
  type :: TException
    type(TPropagationPath) :: propagationPath
    class(TExceptionCause), allocatable :: cause
    logical :: active = .true.
  contains
    procedure, non_overridable :: addLocation => TException_addLocation
    procedure, non_overridable :: writeTo => TException_writeTo
    procedure, non_overridable :: deactivate => TException_deactivate
    final :: TException_final
  end type TException


  !> General error, base class for more specific error classes.
  !!
  !! Avoid if possible and use more specific sub-classes instead.
  !!
  type, extends(TExceptionCause) :: TGeneralError
    character(:), allocatable :: message
  contains
    procedure :: asChar => TGeneralError_asChar
  end type TGeneralError


  !> Error raised due to invalid input (either syntactically or semantically incorrect)
  type, extends(TGeneralError) :: TInputError
  end type TInputError


  !> The model set up is not solvable (e.g. singular, non-convergent, etc.)
  type, extends(TGeneralError) :: TModelError
  end type TModelError


  !> Unrecoverable system error (e.g. I/O or memory problems)
  type, extends(TGeneralError) :: TSystemError
  end type TSystemError


  !> Internal error (which in optimal case should never occur)
  type, extends(TGeneralError) :: TInternalError
  end type TInternalError


  !> Normal but early termination (can be used to propagate up early termination as an exception)
  type, extends(TExceptionCause) :: TEarlyExit
  contains
    procedure :: asChar => TEarlyExit_asChar
  end type TEarlyExit


  !> Whether program should stop, if an exception was active when going out of scope
  logical, parameter :: stopOnActiveException = ${FORTRAN_LOGICAL(DEBUG > 0)}$


contains

  !> Creates (allocates and initializes) an exception
  subroutine TException_create(this, cause, file, line)

    !> Instance
    type(TException), allocatable, intent(out) :: this

    !> Cause of the exception
    class(TExceptionCause), intent(in) :: cause

    !> File where the exception had been triggered
    character(*), optional, intent(in) :: file

    !> Line where the exception had been triggered
    integer, optional, intent(in) :: line

    allocate(this)
    this%cause = cause
    call this%propagationPath%addItem(file, line)

  end subroutine TException_create


  !> Finalizer of the exception, might stop the code, if the exception is still active
  subroutine TException_final(this)

    !> Instance
    type(TException), intent(inout) :: this

    if (stopOnActiveException .and. this%active) then
      write(stdErr, "(a)") "*** Unhandled exception went out of scope ***"
      call this%writeTo(stdErr, withPropagationPath=.true.)
      error stop 255
    end if

  end subroutine TException_final


  !> Deactivates an exception, so that no error stop occurs, if it goes out of scope
  subroutine TException_deactivate(this)
    class(TException), intent(inout) :: this

    this%active = .false.

  end subroutine TException_deactivate


  !> Appends an additional location to the propagation path of the exception
  subroutine TException_addLocation(this, file, line)

    !> Instance
    class(TException), intent(inout) :: this

    !> File name
    character(*), optional, intent(in) :: file

    !> Line number
    integer, optional, intent(in) :: line

    call this%propagationPath%addItem(file=file, line=line)

  end subroutine TException_addLocation


  !> Writes the exception to a file unit
  subroutine TException_writeTo(this, unit, withPropagationPath)

    !> Instance
    class(TException), intent(in) :: this

    !> File unit
    integer, intent(in) :: unit

    !> Whether the exception propagation path should also be written
    logical, optional :: withPropagationPath

    if (present(withPropagationPath)) then
      if (withPropagationPath) call this%propagationPath%writeTo(unit)
    end if
    write(unit, "(a)") this%cause%asChar()

  end subroutine TException_writeTo


  !> Adds a new path item to the propagation path
  subroutine TPropagationPath_addItem(this, file, line)

    !> Instance
    class(TPropagationPath), intent(inout) :: this

    !> File where the exception had been propagated
    character(*), optional, intent(in) :: file

    !> Line where the exception had been propagated
    integer, optional, intent(in) :: line

    type(TPathItem), allocatable :: item

    item = TPathItem(file=file, line=line)
    if (allocated(this%lastItem)) call move_alloc(this%lastItem, item%next)
    call move_alloc(item, this%lastItem)

  end subroutine TPropagationPath_addItem


  !> Writes the propagation path to a file
  subroutine TPropagationPath_writeTo(this, unit)

    !> Instance
    class(TPropagationPath), target, intent(in) :: this

    !> File unit
    integer, intent(in) :: unit

    type(TPathItem), pointer :: pItem

    if (.not. allocated(this%lastItem)) return
    pItem => this%lastItem
    do
      write(unit, "(a, 1x, a)") "File:", pItem%file
      write(unit, "(a, 1x, i0)") "Line:", pItem%line
      if (.not. allocated(pItem%next)) exit
      write(unit, "(a)") "|"
      pItem => pItem%next
    end do

  end subroutine TPropagationPath_writeTo


  !> Returns the character representation of the exception cause
  function TGeneralError_asChar(this) result(str)

    !> Instance
    class(TGeneralError), intent(in) :: this

    !> Character representation of the cause
    character(:), allocatable :: str

    str = this%message

  end function TGeneralError_asChar


  !> Returns the character representation of the exception cause
  function TEarlyExit_asChar(this) result(msg)

    !> Instance
    class(TEarlyExit), intent(in) :: this

    !> !> Character representation of the cause
    character(:), allocatable :: msg

    msg = "Terminating program."

  end function TEarlyExit_asChar

end module dftbp_common_exception
