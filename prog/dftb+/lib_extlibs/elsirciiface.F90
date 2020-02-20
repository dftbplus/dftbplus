!--------------------------------------------------------------------------------------------------!
!  DFTB+: general package for performing fast atomistic simulations                                !
!  Copyright (C) 2006 - 2020  DFTB+ developers group                                               !
!                                                                                                  !
!  See the LICENSE file for terms of usage and distribution.                                       !
!--------------------------------------------------------------------------------------------------!

#:include 'common.fypp'


!> Interface wrapper for the ELSI_RCI reverse communication interface eigensolver library
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

  public :: withElsiRCI, rci_init, rci_solve_allocate, rci_solve, rci_solve_deallocate

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
  ! Placeholder routines when compiled without ELSI_RCI support
  !

  subroutine rci_init(r_h, solver, n_basis, n_state, tol_iter, max_iter, verbose)
    type(rci_handle), intent(out) :: r_h
    integer(i4), intent(in) :: solver
    integer(i4), intent(in) :: n_basis
    integer(i4), intent(in) :: n_state
    real(r8), intent(in) :: tol_iter
    integer(i4), intent(in) :: max_iter
    integer(i4), intent(in) :: verbose
    call stubError("rci_init")
  end subroutine rci_init

  subroutine rci_solve_allocate(r_h, ijob, iS, task)
    type(rci_handle), intent(in) :: r_h
    integer, intent(out) :: task
    integer, intent(inout) :: ijob
    type(rci_instr), intent(inout) :: iS
    call stubError("rci_solve_allocate")
  end subroutine rci_solve_allocate

  subroutine rci_solve(r_h, ijob, iS, task, resvec)
    type(rci_handle), intent(inout) :: r_h
    integer, intent(inout) :: ijob
    type(rci_instr), intent(inout) :: iS
    integer, intent(out) :: task
    real(r8), intent(inout) :: resvec(:)
    call stubError("rci_solve")
  end subroutine rci_solve

  subroutine rci_solve_deallocate(r_h, ijob, iS, task)
    type(rci_handle), intent(in) :: r_h
    integer, intent(out) :: task
    integer, intent(inout) :: ijob
    type(rci_instr), intent(inout) :: iS
    call stubError("rci_solve_deallocate")
  end subroutine rci_solve_deallocate

#:endif

end module dftbp_elsirciiface
