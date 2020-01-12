!--------------------------------------------------------------------------------------------------!
!  DFTB+: general package for performing fast atomistic simulations                                !
!  Copyright (C) 2006 - 2020  DFTB+ developers group                                               !
!                                                                                                  !
!  See the LICENSE file for terms of usage and distribution.                                       !
!--------------------------------------------------------------------------------------------------!

#:include "common.fypp"

#:set GCMD_RANKS = range(0, 3)
#:set GCMD_TYPES = [('integer', 'Int'), ('real(dp)', 'Real')]

!> Exporting the functionality from the plumed library
module dftbp_plumed
  use, intrinsic :: iso_c_binding, only : c_char, c_ptr, c_loc
  use dftbp_accuracy, only : dp
  use dftbp_message, only : error
  implicit none
  private

  public :: plumedInit, plumedFinal, plumedGlobalCmdVal, plumedGlobalCmdPtr
  public :: withPlumed


  interface plumedGlobalCmdVal
    module procedure plumedGlobalCmdValChar
    #:for _, SUFFIX in GCMD_TYPES
      #:for RANK in GCMD_RANKS
        module procedure :: plumedGlobalCmdVal${SUFFIX}$${RANK}$
      #:endfor
    #:endfor
  end interface plumedGlobalCmdVal


  interface plumedGlobalCmdPtr
    #:for _, SUFFIX in GCMD_TYPES
      #:for RANK in GCMD_RANKS
        module procedure :: plumedGlobalCmdPtr${SUFFIX}$${RANK}$
      #:endfor
    #:endfor
  end interface plumedGlobalCmdPtr


  !> Whether package was build with PLUMED support
  logical, parameter :: withPlumed = ${FORTRAN_LOGICAL(WITH_PLUMED)}$


  ! Explicit interfaces for the API provided by PLUMED
  interface
    subroutine plumed_f_gcreate() bind(C, name='plumed_f_gcreate')
    end subroutine plumed_f_gcreate

    subroutine plumed_f_gfinalize() bind(C, name='plumed_f_gfinalize')
    end subroutine plumed_f_gfinalize

    subroutine plumed_f_gcmd(key, val) bind(C, name='plumed_f_gcmd')
      import :: c_char, c_ptr
      implicit none
      character(kind=c_char), intent(in) :: key(*)
      type(c_ptr), value :: val
    end subroutine plumed_f_gcmd
  end interface

  ! Number of global instances (should be 0 or 1)
  integer :: instances = 0

contains

  !> Initializes plumed.
  !>
  !> Note: Currently the global PLUMED interfaces is used, which only allows
  !> one PLUMED instance at a time. Once plumedInit() had been called, it should not be called
  !> again until the instance had been destroyed by calling plumedFinal().
  !>
  subroutine plumedInit()
    #:if WITH_PLUMED
      if (instances == 0) then
        call plumed_f_gcreate()
        instances = instances + 1
      else
        call error("Internal error: Plumed interface only allows one instance at a time")
      end if
    #:else
      call stubError("plumedInit")
    #:endif
  end subroutine plumedInit


  !> Destroys the PLUMED instance.
  subroutine plumedFinal()
    #:if WITH_PLUMED
    if (instances == 1) then
      call plumed_f_gfinalize()
      instances = instances - 1
    else
      call error("Internal error: plumedFinal() called before plumedInit()")
    end if
    #:else
      call stubError("plumedFinal")
    #:endif

  end subroutine plumedFinal


  #:for TYPE, SUFFIX in GCMD_TYPES
    #:for RANK in GCMD_RANKS
      #:set CONTIGUOUS = ", contiguous" if RANK > 0 else ""

      !> Wrapper for passing a value to a global PLUMED instance.
      !>
      !> NOTE: This wrapper should only be used to pass values to PLUMED which are
      !> are immediately COPIED in PLUMED before returning, as the argument may contain
      !> temporary expression
      !>
      subroutine plumedGlobalCmdVal${SUFFIX}$${RANK}$(key, val)

        !> Key (will be automatically extended with the necessary termination character).
        character(len=*, kind=c_char), intent(in) :: key

        !> Value to pass.
        ${TYPE}$, target, intent(in) ${CONTIGUOUS}$ :: val${FORTRAN_ARG_DIM_SUFFIX(RANK)}$

        #:if WITH_PLUMED
          call plumed_f_gcmd(key // char(0), c_loc(val))
        #:else
          call stubError("plumedGlobalCmdVal${SUFFIX}$${RANK}$")
        #:endif

      end subroutine plumedGlobalCmdVal${SUFFIX}$${RANK}$

    #:endfor
  #:endfor


  !> Wrapper for passing a value to a global PLUMED instance (character version).
  !>
  !> NOTE: This wrapper should only be used to pass values to PLUMED which are
  !> are immediately COPIED in PLUMED before returning, as the argument may contain
  !> temporary expression.
  !>
  subroutine plumedGlobalCmdValChar(key, val)

    !> Key (will be automatically extended with the necessary termination character)
    character(len=*, kind=c_char), intent(in) :: key

    !> Value to pass (will be automatically extended with the necessary termination character)
    character(len=*, kind=c_char), intent(in) :: val

    character(len=(len(val) + 1), kind=c_char), target :: buffer

    #:if WITH_PLUMED
      buffer = val // char(0)
      call plumed_f_gcmd(key // char(0), c_loc(buffer))
    #:else
      call stubError("plumedGlobalCmdValChar")
    #:endif

  end subroutine plumedGlobalCmdValChar


  #:for TYPE, SUFFIX in GCMD_TYPES
    #:for RANK in GCMD_RANKS
      #:set CONTIGUOUS = ", contiguous" if RANK > 0 else ""

      !> Wrapper for passing a reference to a global PLUMED instance.
      !>
      !> NOTE: This wrapper passes the adress of the value object. Make sure, the object
      !> exists long enough, that PLUMED can access and eventually modify it when necessary.
      !>
      subroutine plumedGlobalCmdPtr${SUFFIX}$${RANK}$(key, val)

        !> Key (will be automatically extended with the necessary termination character)
        character(len=*, kind=c_char), intent(in) :: key

        !> Object which should be passed as a reference.
        !> Contains workaround for bug in Intel 19 compiler (pointer => target)
        ${TYPE}$, target, intent(in) ${CONTIGUOUS}$ :: val${FORTRAN_ARG_DIM_SUFFIX(RANK)}$

        #:if WITH_PLUMED
          call plumed_f_gcmd(key // char(0), c_loc(val))
        #:else
          call stubError("plumedGlobalCmdPtr${SUFFIX}$${RANK}$")
        #:endif

      end subroutine plumedGlobalCmdPtr${SUFFIX}$${RANK}$

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
