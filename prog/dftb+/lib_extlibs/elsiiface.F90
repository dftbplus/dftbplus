!--------------------------------------------------------------------------------------------------!
!  DFTB+: general package for performing fast atomistic simulations                                !
!  Copyright (C) 2006 - 2020  DFTB+ developers group                                               !
!                                                                                                  !
!  See the LICENSE file for terms of usage and distribution.                                       !
!--------------------------------------------------------------------------------------------------!

#:include 'common.fypp'


!> Interface wrapper for the ELSI library
!>
module dftbp_elsiiface
  use dftbp_accuracy, only : dp
#:if WITH_ELSI
  use elsi
#:else
  use iso_c_binding, only : r8 => c_double, i4 => c_int32_t
  use dftbp_message, only : error
#:endif
  implicit none
  private

  public :: withElsi, withPexsi
  public :: elsi_handle, elsi_rw_handle
  public :: elsi_init, elsi_reinit, elsi_finalize
  public :: elsi_init_rw, elsi_finalize_rw
  public :: elsi_set_csc, elsi_set_csc_blk, elsi_set_sparsity_mask
  public :: elsi_set_zero_def, elsi_set_rw_zero_def
  public :: elsi_dm_real, elsi_dm_complex
  public :: elsi_dm_real_sparse, elsi_dm_complex_sparse
  public :: elsi_get_edm_real, elsi_get_edm_complex
  public :: elsi_get_edm_real_sparse, elsi_get_edm_complex_sparse
  public :: elsi_write_mat_real_sparse, elsi_write_mat_complex_sparse
  public :: elsi_write_mat_real, elsi_write_mat_complex
  public :: elsi_get_pexsi_mu_min, elsi_get_pexsi_mu_max
  public :: elsi_set_pexsi_mu_min, elsi_set_pexsi_mu_max
  public :: elsi_get_mu, elsi_get_entropy
  public :: elsi_set_mu_mp_order, elsi_set_mu_broaden_width, elsi_set_mu_broaden_scheme
  public :: elsi_set_elpa_solver
  public :: elsi_set_omm_flavor, elsi_set_omm_n_elpa, elsi_set_omm_tol
  public :: elsi_set_pexsi_method
  public :: elsi_set_pexsi_np_per_pole, elsi_set_pexsi_temp, elsi_set_pexsi_n_pole
  public :: elsi_set_pexsi_n_mu, elsi_set_pexsi_np_symbo, elsi_set_pexsi_delta_e
  public :: elsi_set_ntpoly_method, elsi_set_ntpoly_filter, elsi_set_ntpoly_tol
  public :: elsi_set_spin, elsi_set_kpoint
  public :: elsi_set_output, elsi_set_output_log
  public :: elsi_ev_complex, elsi_ev_real
  public :: elsi_set_rw_csc
  public :: elsi_set_rw_mpi, elsi_set_rw_blacs
  public :: elsi_set_mpi, elsi_set_mpi_global
  public :: elsi_set_blacs
  public :: elsi_set_sing_check
  public :: elsi_get_version, elsi_get_datestamp


  !> Whether code was built with ELSI support
  logical, parameter :: withElsi = #{if WITH_ELSI}# .true. #{else}# .false. #{endif}#

  !> Whether code was built with PEXSI support
  logical, parameter :: withPexsi = #{if WITH_PEXSI}# .true. #{else}# .false. #{endif}#


#:if not WITH_ELSI

  ! Placeholder types when compiled without ELSI support

  type :: elsi_handle
  end type elsi_handle

  type :: elsi_rw_handle
  end type elsi_rw_handle

#:endif


contains

#:if not WITH_ELSI

  !> Generates error message, if a stub was called
  subroutine stubError(routineName)
    character(*), intent(in) :: routineName

    call error("Internal error: " // trim(routineName) // "() called in a build without ELSI&
        & support")

  end subroutine stubError


  !
  ! Placeholder routines when compiled without ELSI support
  !

  subroutine elsi_init(eh, solver, parallel_mode, matrix_format, n_basis, n_electron, n_state)
    type(elsi_handle), intent(out) :: eh
    integer(i4), intent(in) :: solver
    integer(i4), intent(in) :: parallel_mode
    integer(i4), intent(in) :: matrix_format
    integer(i4), intent(in) :: n_basis
    real(r8), intent(in) :: n_electron
    integer(i4), intent(in) :: n_state
    call stubError("elsi_init")
  end subroutine elsi_init

  subroutine elsi_reinit(eh)
    type(elsi_handle), intent(inout) :: eh
    call stubError("elsi_reinit")
  end subroutine elsi_reinit

  subroutine elsi_finalize(eh)
    type(elsi_handle), intent(inout) :: eh
    call stubError("elsi_finalize")
  end subroutine elsi_finalize


  subroutine elsi_set_csc(eh, nnz_g, nnz_l, n_lcol, row_ind, col_ptr)
    type(elsi_handle), intent(inout) :: eh
    integer(i4), intent(in) :: nnz_g
    integer(i4), intent(in) :: nnz_l
    integer(i4), intent(in) :: n_lcol
    integer(i4), intent(in) :: row_ind(:)
    integer(i4), intent(in) :: col_ptr(:)
    call stubError("elsi_set_csc")
  end subroutine elsi_set_csc


  subroutine elsi_set_csc_blk(eh, blk)
    type(elsi_handle), intent(inout) :: eh
    integer(i4), intent(in) :: blk
    call stubError("elsi_set_csc_blk")
  end subroutine elsi_set_csc_blk

  subroutine elsi_set_sparsity_mask(eh, msk)
    type(elsi_handle), intent(inout) :: eh
    integer(i4), intent(in) :: msk
    call stubError("elsi_set_sparsity_mask")
  end subroutine elsi_set_sparsity_mask

  subroutine elsi_set_zero_def(eh, zero)
    type(elsi_handle), intent(inout) :: eh
    real(r8), intent(in) :: zero
  end subroutine elsi_set_zero_def

  subroutine elsi_set_rw_zero_def(eh, zero)
    type(elsi_rw_handle), intent(inout) :: eh
    real(r8), intent(in) :: zero
  end subroutine elsi_set_rw_zero_def


  subroutine elsi_dm_real(eh, ham, ovlp, dm, ebs)
    type(elsi_handle), intent(inout) :: eh
    real(r8), intent(inout) :: ham(:,:)
    real(r8), intent(inout) :: ovlp(:,:)
    real(r8), intent(inout) :: dm(:,:)
    real(r8), intent(inout) :: ebs
    call stubError("elsi_dm_real")
  end subroutine elsi_dm_real


  subroutine elsi_dm_complex(eh, ham, ovlp, dm, ebs)
    type(elsi_handle), intent(inout) :: eh
    complex(r8), intent(inout) :: ham(:,:)
    complex(r8), intent(inout) :: ovlp(:,:)
    complex(r8), intent(inout) :: dm(:,:)
    real(r8), intent(inout) :: ebs
    call stubError("elsi_dm_complex")
  end subroutine elsi_dm_complex


  subroutine elsi_dm_real_sparse(eh, ham, ovlp, dm, ebs)
    type(elsi_handle), intent(inout) :: eh
    real(r8), intent(inout) :: ham(:)
    real(r8), intent(inout) :: ovlp(:)
    real(r8), intent(inout) :: dm(:)
    real(r8), intent(inout) :: ebs
    call stubError("elsi_dm_real_sparse")
  end subroutine elsi_dm_real_sparse


  subroutine elsi_dm_complex_sparse(eh, ham, ovlp, dm, ebs)
    type(elsi_handle), intent(inout) :: eh
    complex(r8), intent(inout) :: ham(:)
    complex(r8), intent(inout) :: ovlp(:)
    complex(r8), intent(inout) :: dm(:)
    real(r8), intent(inout) :: ebs
    call stubError("elsi_dm_complex_sparse")
  end subroutine elsi_dm_complex_sparse


  subroutine elsi_get_edm_real(eh, edm)
    type(elsi_handle), intent(inout) :: eh
    real(r8), intent(out) :: edm(:,:)
    call stubError("elsi_get_edm_real")
  end subroutine elsi_get_edm_real


  subroutine elsi_get_edm_complex(eh, edm)
    type(elsi_handle), intent(inout) :: eh
    complex(r8), intent(out) :: edm(:,:)
    call stubError("elsi_get_edm_complex")
  end subroutine elsi_get_edm_complex


  subroutine elsi_get_edm_real_sparse(eh, edm)
    type(elsi_handle), intent(inout) :: eh
    real(r8), intent(out) :: edm(:)
    call stubError("elsi_get_edm_real_sparse")
  end subroutine elsi_get_edm_real_sparse


  subroutine elsi_get_edm_complex_sparse(eh, edm)
    type(elsi_handle), intent(inout) :: eh
    complex(r8), intent(out) :: edm(:)
    call stubError("elsi_get_edm_complex_sparse")
  end subroutine elsi_get_edm_complex_sparse


  subroutine elsi_get_mu(eh, mu)
    type(elsi_handle), intent(inout) :: eh
    real(r8), intent(out) :: mu
    call stubError("elsi_get_mu")
  end subroutine elsi_get_mu


  subroutine elsi_get_entropy(eh, entropy)
    type(elsi_handle), intent(inout) :: eh
    real(r8), intent(out) :: entropy
    call stubError("elsi_get_entropy")
  end subroutine elsi_get_entropy


  subroutine elsi_ev_real(eh, ham, ovlp, eval, evec)
    type(elsi_handle), intent(inout) :: eh
    real(r8), intent(inout) :: ham(:,:)
    real(r8), intent(inout) :: ovlp(:,:)
    real(r8), intent(inout) :: eval(:)
    real(r8), intent(inout) :: evec(:,:)
    call stubError("elsi_ev_real")
  end subroutine elsi_ev_real


  subroutine elsi_ev_complex(eh, ham, ovlp, eval, evec)
    type(elsi_handle), intent(inout) :: eh
    complex(r8), intent(inout) :: ham(:,:)
    complex(r8), intent(inout) :: ovlp(:,:)
    real(r8), intent(inout) :: eval(:)
    complex(r8), intent(inout) :: evec(:,:)
    call stubError("elsi_ev_complex")
  end subroutine elsi_ev_complex


  subroutine elsi_set_blacs(eh, blacs_ctxt, block_size)
    type(elsi_handle), intent(inout) :: eh
    integer(i4), intent(in) :: blacs_ctxt
    integer(i4), intent(in) :: block_size
    call stubError("elsi_set_blacs")
  end subroutine elsi_set_blacs


  subroutine elsi_set_elpa_solver(eh, solver)
    type(elsi_handle), intent(inout) :: eh
    integer(i4), intent(in) :: solver
    call stubError("elsi_set_set_elpa_solver")
  end subroutine elsi_set_elpa_solver


  subroutine elsi_set_kpoint(eh, n_kpt, i_kpt, i_wt)
    type(elsi_handle), intent(inout) :: eh
    integer(i4), intent(in) :: n_kpt
    integer(i4), intent(in) :: i_kpt
    real(r8), intent(in) :: i_wt
    call stubError("elsi_set_kpoint")
  end subroutine elsi_set_kpoint


  subroutine elsi_set_spin(eh, n_spin, i_spin)
    type(elsi_handle), intent(inout) :: eh
    integer(i4), intent(in) :: n_spin
    integer(i4), intent(in) :: i_spin
    call stubError("elsi_set_spin")
  end subroutine elsi_set_spin


  subroutine elsi_set_mpi(eh, comm)
    type(elsi_handle), intent(inout) :: eh
    integer(i4), intent(in) :: comm
    call stubError("elsi_set_mpi")
  end subroutine elsi_set_mpi


  subroutine elsi_set_mpi_global(eh, comm_all)
    type(elsi_handle), intent(inout) :: eh
    integer(i4), intent(in) :: comm_all
    call stubError("elsi_set_mpi_global")
  end subroutine elsi_set_mpi_global


  subroutine elsi_set_sing_check(eh, illcond_check)
    type(elsi_handle), intent(inout) :: eh
    integer(i4), intent(in) :: illcond_check
    call stubError("elsi_set_sing_check")
  end subroutine elsi_set_sing_check


  subroutine elsi_set_mu_broaden_scheme(eh, broaden_scheme)
    type(elsi_handle), intent(inout) :: eh
    integer(i4), intent(in) :: broaden_scheme
    call stubError("elsi_set_mu_broaden_scheme")
  end subroutine elsi_set_mu_broaden_scheme


  subroutine elsi_set_mu_broaden_width(eh, broaden_width)
    type(elsi_handle), intent(inout) :: eh
    real(r8), intent(in) :: broaden_width
    call stubError("elsi_set_mu_broaden_width")
  end subroutine elsi_set_mu_broaden_width


  subroutine elsi_set_mu_mp_order(eh, mp_order)
    type(elsi_handle), intent(inout) :: eh
    integer(i4), intent(in) :: mp_order
    call stubError("elsi_set_mu_broaden_width")
  end subroutine elsi_set_mu_mp_order


  subroutine elsi_set_ntpoly_method(eh, method)
    type(elsi_handle), intent(inout) :: eh
    integer(i4), intent(in) :: method
    call stubError("elsi_set_ntpoly_method")
  end subroutine elsi_set_ntpoly_method


  subroutine elsi_set_ntpoly_filter(eh, filter)
    type(elsi_handle), intent(inout) :: eh
    real(r8), intent(in) :: filter
    call stubError("elsi_set_ntpoly_filter")
  end subroutine elsi_set_ntpoly_filter


  subroutine elsi_set_ntpoly_tol(eh, tol)
    type(elsi_handle), intent(inout) :: eh
    real(r8), intent(in) :: tol
    call stubError("elsi_set_ntpoly_tol")
  end subroutine elsi_set_ntpoly_tol


  subroutine elsi_set_omm_flavor(eh, flavor)
    type(elsi_handle), intent(inout) :: eh
    integer(i4), intent(in) :: flavor
    call stubError("elsi_set_omm_flavor")
  end subroutine elsi_set_omm_flavor


  subroutine elsi_set_omm_n_elpa(eh, n_elpa)
    type(elsi_handle), intent(inout) :: eh
    integer(i4), intent(in) :: n_elpa
    call stubError("elsi_set_omm_n_elpa")
  end subroutine elsi_set_omm_n_elpa


  subroutine elsi_set_omm_tol(eh, tol)
    type(elsi_handle), intent(inout) :: eh
    real(r8), intent(in) :: tol
    call stubError("elsi_set_omm_tol")
  end subroutine elsi_set_omm_tol


  subroutine elsi_set_pexsi_method(eh, method)
    type(elsi_handle), intent(inout) :: eh
    integer(i4), intent(in) :: method
    call stubError("elsi_set_pexsi_method")
  end subroutine elsi_set_pexsi_method


  subroutine elsi_set_pexsi_mu_min(eh, mu_min)
    type(elsi_handle), intent(inout) :: eh
    real(r8), intent(in) :: mu_min
    call stubError("elsi_set_pexsi_mu_min")
  end subroutine elsi_set_pexsi_mu_min


  subroutine elsi_set_pexsi_mu_max(eh, mu_max)
    type(elsi_handle), intent(inout) :: eh
    real(r8), intent(in) :: mu_max
    call stubError("elsi_set_pexsi_mu_max")
  end subroutine elsi_set_pexsi_mu_max


  subroutine elsi_get_pexsi_mu_min(eh, mu_min)
    type(elsi_handle), intent(inout) :: eh
    real(r8), intent(out) :: mu_min
    call stubError("elsi_get_pexsi_mu_min")
  end subroutine elsi_get_pexsi_mu_min


  subroutine elsi_get_pexsi_mu_max(eh, mu_max)
    type(elsi_handle), intent(inout) :: eh
    real(r8), intent(out) :: mu_max
    call stubError("elsi_get_pexsi_mu_max")
  end subroutine elsi_get_pexsi_mu_max


  subroutine elsi_set_pexsi_delta_e(eh, delta_e)
    type(elsi_handle), intent(inout) :: eh
    real(r8), intent(in) :: delta_e
    call stubError("elsi_get_pexsi_delta_e")
  end subroutine elsi_set_pexsi_delta_e


  subroutine elsi_set_pexsi_temp(eh, temp)
    type(elsi_handle), intent(inout) :: eh
    real(r8), intent(in) :: temp
    call stubError("elsi_set_pexsi_temp")
  end subroutine elsi_set_pexsi_temp


  subroutine elsi_set_pexsi_n_mu(eh, n_mu)
    type(elsi_handle), intent(inout) :: eh
    integer(i4), intent(in) :: n_mu
    call stubError("elsi_set_pexsi_n_mu")
  end subroutine elsi_set_pexsi_n_mu


  subroutine elsi_set_pexsi_np_symbo(eh, np_symbo)
    type(elsi_handle), intent(inout) :: eh
    integer(i4), intent(in) :: np_symbo
    call stubError("elsi_set_pexsi_np_symbo")
  end subroutine elsi_set_pexsi_np_symbo


  subroutine elsi_set_pexsi_n_pole(eh, n_pole)
    type(elsi_handle), intent(inout) :: eh
    integer(i4), intent(in) :: n_pole
    call stubError("elsi_set_pexsi_n_pole")
  end subroutine elsi_set_pexsi_n_pole


  subroutine elsi_set_pexsi_np_per_pole(eh, np_per_pole)
    type(elsi_handle), intent(inout) :: eh
    integer(i4), intent(in) :: np_per_pole
    call stubError("elsi_set_pexsi_np_per_pole")
  end subroutine elsi_set_pexsi_np_per_pole


  subroutine elsi_set_output(eh, output)
    type(elsi_handle), intent(inout) :: eh
    integer(i4), intent(in) :: output
    call stubError("elsi_set_output")
  end subroutine elsi_set_output


  subroutine elsi_set_output_log(eh, output_log)
    type(elsi_handle), intent(inout) :: eh
    integer(i4), intent(in) :: output_log
    call stubError("elsi_set_output_log")
  end subroutine elsi_set_output_log


  subroutine elsi_init_rw(rwh, task, parallel_mode, n_basis, n_electron)
    type(elsi_rw_handle), intent(out) :: rwh
    integer(i4), intent(in) :: task
    integer(i4), intent(in) :: parallel_mode
    integer(i4), intent(in) :: n_basis
    real(r8), intent(in) :: n_electron
    call stubError("elsi_init_rw")
  end subroutine elsi_init_rw


  subroutine elsi_finalize_rw(rwh)
    type(elsi_rw_handle), intent(inout) :: rwh
    call stubError("elsi_finalize_rw")
  end subroutine elsi_finalize_rw


  subroutine elsi_set_rw_csc(rwh, nnz_g, nnz_l_sp, n_lcol_sp)
    type(elsi_rw_handle), intent(inout) :: rwh
    integer(i4), intent(in) :: nnz_g
    integer(i4), intent(in) :: nnz_l_sp
    integer(i4), intent(in) :: n_lcol_sp
    call stubError("elsi_set_rw_csc")
  end subroutine elsi_set_rw_csc


  subroutine elsi_set_rw_mpi(rwh, mpi_comm)
    type(elsi_rw_handle), intent(inout) :: rwh
    integer(i4), intent(in) :: mpi_comm
    call stubError("elsi_set_rw_mpi")
  end subroutine elsi_set_rw_mpi



  subroutine elsi_set_rw_blacs(rwh, blacs_ctxt, block_size)
    type(elsi_rw_handle), intent(inout) :: rwh
    integer(i4), intent(in) :: blacs_ctxt
    integer(i4), intent(in) :: block_size
    call stubError("elsi_set_rw_blacs")
  end subroutine elsi_set_rw_blacs


  subroutine elsi_write_mat_real_sparse(rwh, f_name, row_ind, col_ptr, mat)
    type(elsi_rw_handle), intent(in) :: rwh
    character(*), intent(in) :: f_name
    integer(i4), intent(in) :: row_ind(:)
    integer(i4), intent(in) :: col_ptr(:)
    real(r8), intent(in) :: mat(:)
    call stubError("elsi_write_mat_real_sparse")
  end subroutine elsi_write_mat_real_sparse


  subroutine elsi_write_mat_complex_sparse(rwh, f_name, row_ind, col_ptr, mat)
    type(elsi_rw_handle), intent(in) :: rwh
    character(*), intent(in) :: f_name
    integer(i4), intent(in) :: row_ind(:)
    integer(i4), intent(in) :: col_ptr(:)
    complex(r8), intent(in) :: mat(:)
    call stubError("elsi_write_mat_complex_sparse")
  end subroutine elsi_write_mat_complex_sparse


  subroutine elsi_write_mat_real(rwh, f_name, mat)
    type(elsi_rw_handle), intent(in) :: rwh
    character(*), intent(in) :: f_name
    real(r8), intent(in) :: mat(:,:)
    call stubError("elsi_write_mat_real")
  end subroutine elsi_write_mat_real


  subroutine elsi_write_mat_complex(rwh, f_name, mat)
    type(elsi_rw_handle), intent(in) :: rwh
    character(*), intent(in) :: f_name
    complex(r8), intent(in) :: mat(:,:)
    call stubError("elsi_write_mat_complex")
  end subroutine elsi_write_mat_complex


  subroutine elsi_get_version(major, minor, patch)
    integer, intent(out) :: major, minor, patch
  end subroutine elsi_get_version

  subroutine elsi_get_datestamp(datestamp)
    integer, intent(out) :: datestamp
  end subroutine elsi_get_datestamp

#:endif

end module dftbp_elsiiface
