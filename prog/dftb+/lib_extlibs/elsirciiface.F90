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
  public :: rci_handle, rci_instr
  !> Solver choices
  public :: RCI_SOLVER, RCI_SOLVER_DAVIDSON, RCI_SOLVER_OMM, RCI_SOLVER_PPCG,&
      & RCI_SOLVER_CHEBFILTER
  !> Constant for ijob
  public :: RCI_INIT_IJOB
  !> Tasks
  public :: RCI_NULL, RCI_STOP, RCI_CONVERGE, RCI_ALLOCATE, RCI_DEALLOCATE, RCI_H_MULTI
  public :: RCI_S_MULTI, RCI_P_MULTI, RCI_COPY, RCI_SUBCOPY, RCI_SUBCOL, RCI_SUBROW
  public :: RCI_TRACE, RCI_DOT, RCI_SCALE, RCI_COLSCALE, RCI_ROWSCALE, RCI_AXPY
  public :: RCI_COL_NORM, RCI_GEMM, RCI_HEEV, RCI_HEGV, RCI_POTRF, RCI_TRSM

  !> Whether code was built with ELSI support
  logical, parameter :: withElsiRCI = #{if WITH_ELSI_RCI}# .true. #{else}# .false. #{endif}#

#:if not WITH_ELSI_RCI

  ! Placeholders for when compiled without ELSI_RCI support

  !> Dummy derived type
  type :: rci_handle
    ! Systems
    integer(i4) :: n_basis, n_state, max_n
    logical :: ovlp_is_unit
    ! Iteration
    integer(i4) :: solver, max_iter
    integer(i4) :: n_res = 0
    real(r8) :: tol_iter
    ! Output
    integer(i4) :: total_iter
    real(r8) :: total_energy
    ! OMM
    real(r8) :: omm_est_ub
    ! PPCG
    integer(i4) :: ppcg_sbsize, ppcg_rrstep
    real(r8) :: ppcg_tol_lock
    ! ChebFilter
    real(r8) :: cheb_est_lb, cheb_est_ub
    integer(i4) :: cheb_max_inneriter
  end type rci_handle

  !> Dummy derived type
  type :: rci_instr
    character :: jobz, uplo, side, trH, trS, trP, trA, trB
    integer(i4) :: m, n, k, lda, ldb, ldc, rAoff, cAoff, rBoff, cBoff, Aidx, Bidx, Cidx
    real(r8) :: alpha, beta
  end type rci_instr

  ! Dummy constants from the interface
  ! Constant of solver
  integer(i4), parameter :: RCI_SOLVER = -1
  integer(i4), parameter :: RCI_SOLVER_DAVIDSON = 1
  integer(i4), parameter :: RCI_SOLVER_OMM = 2
  integer(i4), parameter :: RCI_SOLVER_PPCG = 3
  integer(i4), parameter :: RCI_SOLVER_CHEBFILTER = 4

  ! Constant of ijob
  integer(i4), parameter :: RCI_INIT_IJOB = -1

  ! Constant of task
  integer(i4), parameter :: RCI_NULL = 0
  integer(i4), parameter :: RCI_STOP = 1
  integer(i4), parameter :: RCI_CONVERGE = 2
  integer(i4), parameter :: RCI_ALLOCATE = 3
  integer(i4), parameter :: RCI_DEALLOCATE = 4
  integer(i4), parameter :: RCI_H_MULTI = 5
  integer(i4), parameter :: RCI_S_MULTI = 6
  integer(i4), parameter :: RCI_P_MULTI = 7
  integer(i4), parameter :: RCI_COPY = 8
  integer(i4), parameter :: RCI_SUBCOPY = 9
  integer(i4), parameter :: RCI_SUBCOL = 10
  integer(i4), parameter :: RCI_SUBROW = 11
  integer(i4), parameter :: RCI_TRACE = 12
  integer(i4), parameter :: RCI_DOT = 13
  integer(i4), parameter :: RCI_SCALE = 14
  integer(i4), parameter :: RCI_COLSCALE = 15
  integer(i4), parameter :: RCI_ROWSCALE = 16
  integer(i4), parameter :: RCI_AXPY = 17
  integer(i4), parameter :: RCI_COL_NORM = 18
  integer(i4), parameter :: RCI_GEMM = 19
  integer(i4), parameter :: RCI_HEEV = 20
  integer(i4), parameter :: RCI_HEGV = 21
  integer(i4), parameter :: RCI_POTRF = 22
  integer(i4), parameter :: RCI_TRSM = 23
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
    integer(i4), intent(out) :: task
    integer(i4), intent(inout) :: ijob
    type(rci_instr), intent(inout) :: iS
    call stubError("rci_solve_allocate")
  end subroutine rci_solve_allocate

  subroutine rci_solve(r_h, ijob, iS, task, resvec)
    type(rci_handle), intent(inout) :: r_h
    integer(i4), intent(inout) :: ijob
    type(rci_instr), intent(inout) :: iS
    integer(i4), intent(out) :: task
    real(r8), intent(inout) :: resvec(:)
    call stubError("rci_solve")
  end subroutine rci_solve

  subroutine rci_solve_deallocate(r_h, ijob, iS, task)
    type(rci_handle), intent(in) :: r_h
    integer(i4), intent(out) :: task
    integer(i4), intent(inout) :: ijob
    type(rci_instr), intent(inout) :: iS
    call stubError("rci_solve_deallocate")
  end subroutine rci_solve_deallocate

#:endif

end module dftbp_elsirciiface
