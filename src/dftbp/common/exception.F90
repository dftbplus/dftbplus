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
  public :: TException, TException_create, TException_destroy
  public :: TExceptionEvent, TEarlyExit, TFatalError, TInternalError


  !> Represents an expecption propagation path item
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


  !> Represents the event which triggered the exception
  type, abstract :: TExceptionEvent
  contains
    procedure(TExceptionEvent_asChar), deferred :: asChar
  end type TExceptionEvent


  abstract interface
    !> Character representation of the event
    function TExceptionEvent_asChar(this) result(res)
      import :: TExceptionEvent
      implicit none

      !> Instance
      class(TExceptionEvent), intent(in) :: this

      !> Character representation of the event
      character(:), allocatable :: res

    end function TExceptionEvent_asChar
  end interface


  !> Represents an exception
  type :: TException
    type(TPropagationPath) :: propagationPath
    class(TExceptionEvent), allocatable :: event
    logical :: active = .true.
  contains
    procedure, non_overridable :: propagate => TException_propagate
    procedure, non_overridable :: writeTo => TException_writeTo
    procedure, non_overridable :: deactivate => TException_deactivate
    final :: TException_final
  end type TException


  !> Unrecoverable fatal error
  type, extends(TExceptionEvent) :: TFatalError
    character(:), allocatable :: message
  contains
      procedure :: asChar => TFatalError_asChar
  end type TFatalError


  !> Internal error (which in optimal case should never occur)
  type, extends(TExceptionEvent) :: TInternalError
    character(:), allocatable :: message
  contains
      procedure :: asChar => TInternalError_asChar
  end type TInternalError


  !> Normal but early termination (can be used to propagate up early termination as an exception)
  type, extends(TExceptionEvent) :: TEarlyExit
  contains
      procedure :: asChar => TEarlyExit_asChar
  end type TEarlyExit


  !> Whether program should stop, if an exception was active when going out of scope
  logical, parameter :: stopOnActiveException = ${FORTRAN_LOGICAL(DEBUG > 0)}$


contains

  !> Creates (allocates and initializes) an exception
  subroutine TException_create(this, event, file, line)

    !> Instance
    type(TException), allocatable, intent(out) :: this

    !> Event triggering the exception
    class(TExceptionEvent), intent(in) :: event

    !> File where the exception had been triggered
    character(*), optional, intent(in) :: file

    !> Line where the exception had been triggered
    integer, optional, intent(in) :: line

    allocate(this)
    this%event = event
    call this%propagationPath%addItem(file, line)

  end subroutine TException_create


  !> Safely destroys (deactivates and deallocates) an exception
  subroutine TException_destroy(this)

    !> Instance
    type(TException), allocatable, intent(inout) :: this

    call this%deactivate()
    deallocate(this)

  end subroutine TException_destroy


  !> Finalizer of the exception, might stop the code, if the exception is still active
  subroutine TException_final(this)

    !> Instance
    type(TException), intent(inout) :: this

    if (stopOnActiveException .and. this%active) then
      write(stdErr, "(a)") "*** Unhandled exception went out of scope ***"
      call this%writeTo(stdErr, withPropagationPath=.true.)
      error stop 1
    end if

  end subroutine TException_final


  !> Deactivates an exception, so that no error stop occurs, if it goes out of scope
  subroutine TException_deactivate(this)
    class(TException), intent(inout) :: this

    this%active = .false.

  end subroutine TException_deactivate


  !> Adds propagation information to an exception
  subroutine TException_propagate(this, file, line)

    !> Instance
    class(TException), intent(inout) :: this

    !> File, where the exception had been propagated
    character(*), optional, intent(in) :: file

    !> Line, where the exception had been propagated
    integer, optional, intent(in) :: line

    call this%propagationPath%addItem(file=file, line=line)

  end subroutine TException_propagate


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
    write(unit, "(a)") this%event%asChar()

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


  !> Returns the character representation of the event
  function TFatalError_asChar(this) result(msg)

    !> Instance
    class(TFatalError), intent(in) :: this

    !> Character representation of the event
    character(:), allocatable :: msg

    msg = "FatalError: " // this%message

  end function TFatalError_asChar


  !> Returns the character representation of the event
  function TInternalError_asChar(this) result(msg)

    !> Instance
    class(TInternalError), intent(in) :: this

    !> Character representation of the event
    character(:), allocatable :: msg

    msg = "InternalError: " // this%message

  end function TInternalError_asChar


  !> Returns the character representation of the error
  function TEarlyExit_asChar(this) result(msg)

    !> Instance
    class(TEarlyExit), intent(in) :: this

    !> Character representation of the event
    character(:), allocatable :: msg

    msg = "Program terminated"

  end function TEarlyExit_asChar

end module dftbp_common_exception
