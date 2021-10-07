!--------------------------------------------------------------------------------------------------!
!  DFTB+: general package for performing fast atomistic simulations                                !
!  Copyright (C) 2006 - 2021  DFTB+ developers group                                               !
!                                                                                                  !
!  See the LICENSE file for terms of usage and distribution.                                       !
!--------------------------------------------------------------------------------------------------!

#:include 'common.fypp'

!> Implements a simple error mechanism.
module dftbp_common_error
  use iso_fortran_env, stdErr => error_unit
  implicit none


  private
  public :: TError, TError_create, TError_destroy


  type :: TLinkedLocation
    character(:), allocatable :: file
    integer :: line
    type(TLinkedLocation), allocatable :: previous
  contains
    procedure :: writeLocations => TLinkedLocation_writeLocations
  end type TLinkedLocation


  !> Implements an error type
  !>
  !> The error should be created using TError_create() and destroyed with TError_destruct().
  !> If the error goes out of scope without having been explicitely distructed, it will stop
  !> the program (in DEBUG mode).
  !>
  type :: TError
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
    procedure :: activate => TError_activate
    procedure :: deactivate => TError_deactivate
    procedure :: addPropagationInfo => TError_addPropagationInfo
    #:block DEBUG_CODE
      final :: TError_final
    #:endblock DEBUG_CODE
  end type TError


contains

  !> Creates an error
  subroutine TError_create(this, code, message)

    !> Allocated and initialised instance on exit.
    type(TError), allocatable, intent(out) :: this

    !> Error code
    integer, intent(in) :: code

    !> Error message
    character(*), intent(in) :: message

    allocate(this)
    call this%activate(code, message)

  end subroutine TError_create


  !> Safely destroys an error (after deactivating it).
  subroutine TError_destroy(this)

    !> Instance.
    class(TError), allocatable, intent(inout) :: this

    call this%deactivate()
    deallocate(this)

  end subroutine TError_destroy


  !> Activates an error.
  !>
  !> Note, if compiled in DEBUG mode, an activated error simply stops the code, if it goes
  !> out of scope. Make sure, you deactivate it via the deactivate() method or destroy
  !> it via TError_destroy() before this happens.
  !>
  subroutine TError_activate(this, code, message)

    !> Initialised instance on exit.
    class(TError), intent(out) :: this

    !> Error code
    integer, intent(in) :: code

    !> Error message
    character(*), intent(in) :: message

    this%code = code
    this%message = message
    this%active = .true.

  end subroutine TError_activate


  !> Deactivates the error (going out of scope will not trigger a stop in DEBUG mode).
  subroutine TError_deactivate(this)

    !> Instance
    class(TError), intent(inout) :: this

    this%active = .false.

  end subroutine TError_deactivate


  !> Adds propagation info to the error.
  subroutine TError_addPropagationInfo(this, file, line)

    !> Instance.
    class(TError), intent(inout) :: this

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

  end subroutine TError_addPropagationInfo


#:block DEBUG_CODE

  !> Finalizes the error.
  !>
  !>  Note, that this stops the code, if the error is active.
  !>
  subroutine TError_final(this)
    type(TError), intent(inout) :: this

    if (this%active) then
      write(stdErr, "('Unhandled error went out of scope!')")
      write(stdErr, "('Error message: ', a)") this%message
      write(stdErr, "('Error code: ', i0)") this%code
      if (allocated(this%propagationPath)) then
        write(stdErr, "('Propagation path:')")
        call this%propagationPath%writeLocations(stdErr)
      end if
      error stop
    end if

  end subroutine TError_final

#:endblock DEBUG_CODE


  !> Writes the locations contained in a recursive TLinkedLocation structure.
  subroutine TLinkedLocation_writeLocations(this, unit)

    !> Instance.
    class(TLinkedLocation), intent(in) :: this

    !> Unit to write the information to.
    integer, intent(in) :: unit

    if (allocated(this%previous)) then
      call this%previous%writeLocations(unit)
    end if
    write(unit, "(a, ':', i0)") this%file, this%line

  end subroutine TLinkedLocation_writeLocations


end module  dftbp_common_error