!--------------------------------------------------------------------------------------------------!
!  DFTB+: general package for performing fast atomistic simulations                                !
!  Copyright (C) 2006 - 2020  DFTB+ developers group                                               !
!                                                                                                  !
!  See the LICENSE file for terms of usage and distribution.                                       !
!--------------------------------------------------------------------------------------------------!

#:include 'common.fypp'


!> Interface wrapper for the ELSI_RCI library
!>
module dftbp_elsirciiface
  use dftbp_accuracy, only : dp
#:if WITH_ELSI_RCI
  use elsi_rci
#:else
  use iso_c_binding, only : r8 => c_double, i4 => c_int32_t
  use dftbp_message, only : error
#:endif
  implicit none
  private

  public :: withElsiRCI


  !> Whether code was built with ELSI support
  logical, parameter :: withElsiRCI = #{if WITH_ELSI_RCI}# .true. #{else}# .false. #{endif}#

#:if not WITH_ELSI_RCI

  ! Placeholder types when compiled without ELSI_RCI support

  type :: rci_handle
  end type rci_handle

  type :: rci_instr
  end type rci_instr

#:endif

contains

#:if not WITH_ELSI_RCI

  !> Generates error message, if a stub was called
  subroutine stubError(routineName)
    character(*), intent(in) :: routineName

    call error("Internal error: " // trim(routineName) // "() called in a build without ELSI_RCI&
        & support")

  end subroutine stubError


  !
  ! Placeholder routines when compiled without ELSI support
  !


#:endif

end module dftbp_elsirciiface
