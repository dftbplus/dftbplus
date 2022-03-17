!--------------------------------------------------------------------------------------------------!
!  DFTB+: general package for performing fast atomistic simulations                                !
!  Copyright (C) 2006 - 2021  DFTB+ developers group                                               !
!                                                                                                  !
!  See the LICENSE file for terms of usage and distribution.                                       !
!--------------------------------------------------------------------------------------------------!

#:include 'common.fypp'

!> Implements a simple exception-like mechanism.
module dftbp_common_exception
  use iso_fortran_env, only : stdErr => error_unit
  use dftbp_common_globalenv, only : abortProgram
  implicit none

  private
  public :: TException, TException_create, TException_destroy


  type :: TLinkedLocation
    character(:), allocatable :: file
    integer :: line
    type(TLinkedLocation), allocatable :: previous
  contains
    procedure :: writeLocations => TLinkedLocation_writeLocations
  end type TLinkedLocation


  !> Implements an exception type
  !>
  !> The exception should be created using TException_create() and destroyed with
  !> TException_destruct(). If the exception goes out of scope without having been explicitely
  !> distructed, it will stop the program (in DEBUG mode).
  !>
  type :: TException
    private

    !> Error code
    integer, public :: code

    !> Error message
    character(:), allocatable, public :: message

    !> Path the error had been propagated along (may be unallocated, if not recorded)
    type(TLinkedLocation), allocatable :: propagationPath

    logical :: active = .false.

  contains

    !> Adds propagation info
    procedure :: activate => TException_activate
    procedure :: deactivate => TException_deactivate
    procedure :: addPropagationInfo => TException_addPropagationInfo
    #:block DEBUG_CODE
      final :: TException_final
    #:endblock DEBUG_CODE
  end type TException


contains

  !> Creates an exception.
  subroutine TException_create(this, code, message)

    !> Allocated and initialised instance on exit.
    type(TException), allocatable, intent(out) :: this

    !> Error code
    integer, intent(in) :: code

    !> Error message
    character(*), intent(in) :: message

    allocate(this)
    call this%activate(code, message)

  end subroutine TException_create


  !> Safely destroys an exception (after deactivating it).
  subroutine TException_destroy(this)

    !> Instance.
    class(TException), allocatable, intent(inout) :: this

    call this%deactivate()
    deallocate(this)

  end subroutine TException_destroy


  !> Activates an exception.
  !>
  !> Note, if compiled in DEBUG mode, an activated exception simply stops the code, if it goes
  !> out of scope. Make sure, you deactivate it via the deactivate() method or destroy
  !> it via TException_destroy() before this happens.
  !>
  subroutine TException_activate(this, code, message)

    !> Initialised instance on exit.
    class(TException), intent(out) :: this

    !> Error code
    integer, intent(in) :: code

    !> Error message
    character(*), intent(in) :: message

    this%code = code
    this%message = message
    this%active = .true.

  end subroutine TException_activate


  !> Deactivates the error (going out of scope will not trigger a stop in DEBUG mode).
  subroutine TException_deactivate(this)

    !> Instance
    class(TException), intent(inout) :: this

    this%active = .false.

  end subroutine TException_deactivate


  !> Adds propagation info to the error.
  subroutine TException_addPropagationInfo(this, file, line)

    !> Instance.
    class(TException), intent(inout) :: this

    !> File, in which error was created/propagated.
    character(*), intent(in) :: file

    !> Line where error was created/propagated.
    integer, intent(in) :: line

    type(TLinkedLocation), allocatable :: currentLoc

    allocate(currentLoc)
    currentLoc%file = file
    currentLoc%line = line
    call move_alloc(this%propagationPath, currentLoc%previous)
    call move_alloc(currentLoc, this%propagationPath)

  end subroutine TException_addPropagationInfo


#:block DEBUG_CODE

  !> Finalizes the exception.
  !>
  !>  Note, that this stops the code, if the exception is active.
  !>
  subroutine TException_final(this)
    type(TException), intent(inout) :: this

    if (this%active) then
      write(stdErr, "('Unhandled error went out of scope!')")
      write(stdErr, "('Error message: ', a)") this%message
      write(stdErr, "('Error code: ', i0)") this%code
      if (allocated(this%propagationPath)) then
        write(stdErr, "('Propagation path:')")
        call this%propagationPath%writeLocations(stdErr)
      end if
      call abortProgram()
    end if

  end subroutine TException_final

#:endblock DEBUG_CODE


  !> Writes the locations contained in a recursive TLinkedLocation structure.
  recursive subroutine TLinkedLocation_writeLocations(this, unit)

    !> Instance.
    class(TLinkedLocation), intent(in) :: this

    !> Unit to write the information to.
    integer, intent(in) :: unit

    if (allocated(this%previous)) then
      call this%previous%writeLocations(unit)
    end if
    write(unit, "(a, ':', i0)") this%file, this%line

  end subroutine TLinkedLocation_writeLocations


end module  dftbp_common_exception
