!--------------------------------------------------------------------------------------------------!
!  DFTB+: general package for performing fast atomistic simulations                                !
!  Copyright (C) 2006 - 2019  DFTB+ developers group                                               !
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


  logical, parameter :: withPlumed = ${FORTRAN_LOGICAL(WITH_PLUMED)}$


  interface
    subroutine plumed_f_gcreate() bind(C, name='plumed_f_gcreate')
    end subroutine plumed_f_gcreate

    subroutine plumed_f_gfinalize() bind(C, name='plumed_f_gfinalize')
    end subroutine plumed_f_gfinalize

    subroutine plumed_f_gcmd(key, val) bind(C, name='plumed_f_gcmd')
      import :: c_char, c_ptr
      character(kind=c_char), intent(in) :: key(*)
      type(c_ptr), value :: val
    end subroutine plumed_f_gcmd
  end interface


  integer :: instances = 0

contains

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

      subroutine plumedGlobalCmdVal${SUFFIX}$${RANK}$(key, val)
        character(len=*, kind=c_char), intent(in) :: key
        ${TYPE}$, target, intent(in) ${CONTIGUOUS}$ :: val${FORTRAN_ARG_DIM_SUFFIX(RANK)}$

        #:if WITH_PLUMED
          call plumed_f_gcmd(key // char(0), c_loc(val))
        #:else
          call stubError("plumedGlobalCmdVal${SUFFIX}$${RANK}$")
        #:endif

      end subroutine plumedGlobalCmdVal${SUFFIX}$${RANK}$

    #:endfor
  #:endfor


  subroutine plumedGlobalCmdValChar(key, val)
    character(len=*, kind=c_char), intent(in) :: key
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

      subroutine plumedGlobalCmdPtr${SUFFIX}$${RANK}$(key, val)
        character(len=*, kind=c_char), intent(in) :: key
        ${TYPE}$, pointer, intent(in) ${CONTIGUOUS}$ :: val${FORTRAN_ARG_DIM_SUFFIX(RANK)}$

        #:if WITH_PLUMED
          call plumed_f_gcmd(key // char(0), c_loc(val))
        #:else
          call stubError("plumedGlobalCmdPtr${SUFFIX}$${RANK}$")
        #:endif

      end subroutine plumedGlobalCmdPtr${SUFFIX}$${RANK}$

    #:endfor
  #:endfor


  #:if not WITH_PLUMED
   subroutine stubError(name)
     character(*), intent(in) :: name

     call error("Internal error: Function '" // name // "' called but code was compiled without&
         & PLUMED support")

   end subroutine stubError
 #:endif


end module dftbp_plumed
