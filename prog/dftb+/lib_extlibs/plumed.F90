!--------------------------------------------------------------------------------------------------!
!  DFTB+: general package for performing fast atomistic simulations                                !
!  Copyright (C) 2006 - 2020  DFTB+ developers group                                               !
!                                                                                                  !
!  See the LICENSE file for terms of usage and distribution.                                       !
!--------------------------------------------------------------------------------------------------!

#:include "common.fypp"

#:set PLUMED_CMD_RANKS = range(0, 3)
#:set PLUMED_CMD_TYPES = [('integer', 'Int'), ('real(dp)', 'Real')]


!> Exporting the functionality from the plumed library
module dftbp_plumed
  use, intrinsic :: iso_c_binding, only : c_int, c_char, c_ptr, c_loc
  use dftbp_accuracy, only : dp
  use dftbp_message, only : error
  implicit none
  private


  public :: TPlumedCalc, TPlumedCalc_init, TPlumedCalc_final
  public :: withPlumed


  type, bind(c) :: c_plumed
    type(c_ptr) :: cptr
  end type c_plumed


  type :: TPlumedCalc
    private
    type(c_plumed) :: desc
  contains
    private
    #:for _, SUFFIX in PLUMED_CMD_TYPES
      #:for RANK in PLUMED_CMD_RANKS
        procedure :: sendCmdVal${SUFFIX}$${RANK}$
        generic, public :: sendCmdVal => sendCmdVal${SUFFIX}$${RANK}$
        procedure :: sendCmdPtr${SUFFIX}$${RANK}$
        generic, public :: sendCmdPtr => sendCmdPtr${SUFFIX}$${RANK}$
      #:endfor
    #:endfor
    procedure :: sendCmdValChar
    generic, public :: sendCmdVal => sendCmdValChar
  end type TPlumedCalc


  ! Explicit interfaces for the C-API provided by PLUMED
  interface

    function plumed_create() result(instance) bind(C, name='plumed_create')
      import :: c_plumed
      type(c_plumed) :: instance
    end function plumed_create

    subroutine plumed_cmd(instance, key, val) bind(C, name='plumed_cmd')
      import :: c_plumed, c_char, c_ptr
      type(c_plumed), value :: instance
      character(kind=c_char), intent(in) :: key(*)
      type(c_ptr), value :: val
    end subroutine plumed_cmd

    subroutine plumed_finalize(instance) bind(C, name='plumed_finalize')
      import :: c_plumed
      type(c_plumed), value :: instance
    end subroutine plumed_finalize

    function plumed_installed(instance) result(installed) bind(C, name='plumed_installed')
      import :: c_plumed, c_int
      type(c_plumed), value :: instance
      integer(c_int) :: installed
    end function plumed_installed

  end interface


  !> Whether package was build with PLUMED support
  logical, parameter :: withPlumed = ${FORTRAN_LOGICAL(WITH_PLUMED)}$


contains

  !> Initializes plumed.
  !>
  subroutine TPlumedCalc_init(this)
    type(TPlumedCalc), intent(out) :: this

    #:if WITH_PLUMED
      this%desc = plumed_create()
    #:else
      call stubError("TPlumedCalc_init")
    #:endif

  end subroutine TPlumedCalc_init


  !> Destroys the PLUMED instance.
  subroutine TPlumedCalc_final(this)
    type(TPlumedCalc), intent(inout) :: this

    #:if WITH_PLUMED
      call plumed_finalize(this%desc)
    #:else
      call stubError("TPlumedCalc_final")
    #:endif

  end subroutine TPlumedCalc_final


  #:for TYPE, SUFFIX in PLUMED_CMD_TYPES
    #:for RANK in PLUMED_CMD_RANKS
      #:set CONTIGUOUS = ", contiguous" if RANK > 0 else ""

      !> Wrapper for passing a value to a PLUMED instance.
      !>
      !> NOTE: This wrapper should only be used to pass values to PLUMED which are
      !> are immediately COPIED in PLUMED before returning, as the argument may contain
      !> temporary expression
      !>
      subroutine sendCmdVal${SUFFIX}$${RANK}$(this, key, val)

        !> Instance
        class(TPlumedCalc), intent(inout) :: this

        !> Key (will be automatically extended with the necessary termination character).
        character(len=*, kind=c_char), intent(in) :: key

        !> Value to pass.
        ${TYPE}$, target, intent(in) ${CONTIGUOUS}$ :: val${FORTRAN_ARG_DIM_SUFFIX(RANK)}$

        #:if WITH_PLUMED
          call plumed_cmd(this%desc, key // char(0), c_loc(val))
        #:else
          call stubError("sendCmdVal${SUFFIX}$${RANK}$")
        #:endif

      end subroutine sendCmdVal${SUFFIX}$${RANK}$

    #:endfor
  #:endfor


  !> Wrapper for passing a value to a PLUMED instance (character version).
  !>
  !> NOTE: This wrapper should only be used to pass values to PLUMED which are
  !> are immediately COPIED in PLUMED before returning, as the argument may contain
  !> temporary expression.
  !>
  subroutine sendCmdValChar(this, key, val)

    !> Instance
    class(TPlumedCalc), intent(inout) :: this
    
    !> Key (will be automatically extended with the necessary termination character)
    character(len=*, kind=c_char), intent(in) :: key

    !> Value to pass (will be automatically extended with the necessary termination character)
    character(len=*, kind=c_char), intent(in) :: val

    character(len=(len(val) + 1), kind=c_char), target :: buffer

    #:if WITH_PLUMED
      buffer = val // char(0)
      call plumed_cmd(this%desc, key // char(0), c_loc(buffer))
    #:else
      call stubError("sendCmdValChar")
    #:endif

  end subroutine sendCmdValChar


  #:for TYPE, SUFFIX in PLUMED_CMD_TYPES
    #:for RANK in PLUMED_CMD_RANKS
      #:set CONTIGUOUS = ", contiguous" if RANK > 0 else ""

      !> Wrapper for passing a reference to a PLUMED instance.
      !>
      !> NOTE: This wrapper passes the adress of the value object. Make sure, the object
      !> exists long enough, that PLUMED can access and eventually modify it when necessary.
      !>
      subroutine sendCmdPtr${SUFFIX}$${RANK}$(this, key, val)

        !> Instance
        class(TPlumedCalc), intent(inout) :: this

        !> Key (will be automatically extended with the necessary termination character)
        character(len=*, kind=c_char), intent(in) :: key

        !> Object which should be passed as a reference.
        !> Contains workaround for bug in Intel 19 compiler (pointer => target)
        ${TYPE}$, target, intent(in) ${CONTIGUOUS}$ :: val${FORTRAN_ARG_DIM_SUFFIX(RANK)}$

        #:if WITH_PLUMED
          call plumed_cmd(this%desc, key // char(0), c_loc(val))
        #:else
          call stubError("sendCmdPtr${SUFFIX}$${RANK}$")
        #:endif

      end subroutine sendCmdPtr${SUFFIX}$${RANK}$

    #:endfor
  #:endfor


  #:if not WITH_PLUMED

    !> Raises an error signalizing the call of a stub-function.
    subroutine stubError(name)

      !> Name of the stub function which was called.
      character(*), intent(in) :: name

      call error("Internal error: Function '" // name // "' called but code was compiled without&
          & PLUMED support")

    end subroutine stubError

  #:endif


end module dftbp_plumed
