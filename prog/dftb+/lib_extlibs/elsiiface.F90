!--------------------------------------------------------------------------------------------------!
!  DFTB+: general package for performing fast atomistic simulations                                !
!  Copyright (C) 2017  DFTB+ developers group                                                      !
!                                                                                                  !
!  See the LICENSE file for terms of usage and distribution.                                       !
!--------------------------------------------------------------------------------------------------!

#:include 'common.fypp'


!> Interface wrapper for the ELSI library
!>
module elsiiface
  use accuracy, only : dp
#:if WITH_ELSI
  use elsi
#:endif
  implicit none
  private

  public :: withElsi, withPexsi
  public :: elsi_handle, elsi_rw_handle
  public :: elsi_init, elsi_finalize
  public :: elsi_init_rw, elsi_finalize_rw
  public :: elsi_set_csc, elsi_set_csc_blk
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

  !
  ! Placeholder routines when compiled without ELSI support
  !

  $:STUB_SUBROUTINE('elsi_init', ['type(elsi_handle)', 'integer', 'integer', 'integer', 'integer',&
      & 'real(dp)', 'integer'])

  $:STUB_SUBROUTINE('elsi_finalize', ['type(elsi_handle)'])

  $:STUB_SUBROUTINE('elsi_set_csc', ['type(elsi_handle)', 'integer', 'integer', 'integer',&
      & 'integer, dimension(:)', 'integer, dimension(:)'])

  $:STUB_SUBROUTINE('elsi_set_csc_blk', ['type(elsi_handle)', 'integer'])

  $:STUB_SUBROUTINE('elsi_dm_real', ['type(elsi_handle)', 'real(dp), dimension(:,:)',&
      & 'real(dp), dimension(:,:)', 'real(dp), dimension(:,:)', 'real(dp)'])

  $:STUB_SUBROUTINE('elsi_dm_complex', ['type(elsi_handle)', 'complex(dp), dimension(:,:)',&
      & 'complex(dp), dimension(:,:)', 'complex(dp), dimension(:,:)', 'real(dp)'])

  $:STUB_SUBROUTINE('elsi_dm_real_sparse', ['type(elsi_handle)', 'real(dp), dimension(:)',&
      & 'real(dp), dimension(:)', 'real(dp), dimension(:)', 'real(dp)'])

  $:STUB_SUBROUTINE('elsi_dm_complex_sparse', ['type(elsi_handle)', 'complex(dp), dimension(:)',&
      & 'complex(dp), dimension(:)', 'complex(dp), dimension(:)', 'real(dp)'])

  $:STUB_SUBROUTINE('elsi_get_edm_real', ['type(elsi_handle)', 'real(dp), dimension(:,:)'])

  $:STUB_SUBROUTINE('elsi_get_edm_complex', ['type(elsi_handle)', 'complex(dp), dimension(:,:)'])

  $:STUB_SUBROUTINE('elsi_get_edm_real_sparse', ['type(elsi_handle)', 'real(dp), dimension(:)'])

  $:STUB_SUBROUTINE('elsi_get_edm_complex_sparse', ['type(elsi_handle)',&
      & 'complex(dp), dimension(:)'])

  $:STUB_SUBROUTINE('elsi_get_mu', ['type(elsi_handle)', 'real(dp)'])

  $:STUB_SUBROUTINE('elsi_get_entropy', ['type(elsi_handle)', 'real(dp)'])

  $:STUB_SUBROUTINE('elsi_ev_real', ['type(elsi_handle)', 'real(dp), dimension(:,:)',&
      & 'real(dp), dimension(:,:)', 'real(dp), dimension(:)', 'real(dp), dimension(:,:)'])

  $:STUB_SUBROUTINE('elsi_ev_complex', ['type(elsi_handle)', 'complex(dp), dimension(:,:)',&
      & 'complex(dp), dimension(:,:)', 'real(dp), dimension(:)', 'complex(dp), dimension(:,:)'])

  $:STUB_SUBROUTINE('elsi_set_blacs', ['type(elsi_handle)', 'integer', 'integer'])

  $:STUB_SUBROUTINE('elsi_set_elpa_solver', ['type(elsi_handle)', 'integer'])

  $:STUB_SUBROUTINE('elsi_set_kpoint', ['type(elsi_handle)', 'integer', 'integer', 'real(dp)'])

  $:STUB_SUBROUTINE('elsi_set_spin', ['type(elsi_handle)', 'integer', 'integer'])

  $:STUB_SUBROUTINE('elsi_set_mpi', ['type(elsi_handle)', 'integer'])

  $:STUB_SUBROUTINE('elsi_set_mpi_global', ['type(elsi_handle)', 'integer'])

  $:STUB_SUBROUTINE('elsi_set_sing_check', ['type(elsi_handle)', 'integer'])

  $:STUB_SUBROUTINE('elsi_set_mu_broaden_scheme', ['type(elsi_handle)', 'integer'])

  $:STUB_SUBROUTINE('elsi_set_mu_broaden_width', ['type(elsi_handle)', 'real(dp)'])

  $:STUB_SUBROUTINE('elsi_set_mu_mp_order', ['type(elsi_handle)', 'integer'])

  $:STUB_SUBROUTINE('elsi_set_ntpoly_method', ['type(elsi_handle)', 'integer'])

  $:STUB_SUBROUTINE('elsi_set_ntpoly_filter', ['type(elsi_handle)', 'real(dp)'])

  $:STUB_SUBROUTINE('elsi_set_ntpoly_tol', ['type(elsi_handle)', 'real(dp)'])

  $:STUB_SUBROUTINE('elsi_set_omm_flavor', ['type(elsi_handle)', 'integer'])

  $:STUB_SUBROUTINE('elsi_set_omm_n_elpa', ['type(elsi_handle)', 'integer'])

  $:STUB_SUBROUTINE('elsi_set_omm_tol', ['type(elsi_handle)', 'real(dp)'])

  $:STUB_SUBROUTINE('elsi_get_pexsi_mu_min', ['type(elsi_handle)', 'real(dp)'])

  $:STUB_SUBROUTINE('elsi_get_pexsi_mu_max', ['type(elsi_handle)', 'real(dp)'])

  $:STUB_SUBROUTINE('elsi_set_pexsi_mu_min', ['type(elsi_handle)', 'real(dp)'])

  $:STUB_SUBROUTINE('elsi_set_pexsi_mu_max', ['type(elsi_handle)', 'real(dp)'])

  $:STUB_SUBROUTINE('elsi_set_pexsi_delta_e', ['type(elsi_handle)', 'real(dp)'])

  $:STUB_SUBROUTINE('elsi_set_pexsi_temp', ['type(elsi_handle)', 'real(dp)'])

  $:STUB_SUBROUTINE('elsi_set_pexsi_n_mu', ['type(elsi_handle)', 'integer'])

  $:STUB_SUBROUTINE('elsi_set_pexsi_np_symbo', ['type(elsi_handle)', 'integer'])

  $:STUB_SUBROUTINE('elsi_set_pexsi_n_pole', ['type(elsi_handle)', 'integer'])

  $:STUB_SUBROUTINE('elsi_set_pexsi_np_per_pole', ['type(elsi_handle)', 'integer'])

  $:STUB_SUBROUTINE('elsi_set_output', ['type(elsi_handle)', 'integer'])

  $:STUB_SUBROUTINE('elsi_set_output_log', ['type(elsi_handle)', 'integer'])

  $:STUB_SUBROUTINE('elsi_init_rw', ['type(elsi_rw_handle)', 'integer', 'integer', 'integer',&
      & 'real(dp)'])

  $:STUB_SUBROUTINE('elsi_finalize_rw', ['type(elsi_rw_handle)'])

  $:STUB_SUBROUTINE('elsi_set_rw_csc', ['type(elsi_rw_handle)', 'integer', 'integer', 'integer'])

  $:STUB_SUBROUTINE('elsi_set_rw_mpi', ['type(elsi_rw_handle)', 'integer'])

  $:STUB_SUBROUTINE('elsi_set_rw_blacs', ['type(elsi_rw_handle)', 'integer', 'integer'])

  $:STUB_SUBROUTINE('elsi_write_mat_real_sparse', ['type(elsi_rw_handle)', 'character(*)',&
      & 'integer, dimension(:)', 'integer, dimension(:)', 'real(dp), dimension(:)'])

  $:STUB_SUBROUTINE('elsi_write_mat_complex_sparse', ['type(elsi_rw_handle)', 'character(*)',&
      & 'integer, dimension(:)', 'integer, dimension(:)', 'complex(dp), dimension(:)'])

  $:STUB_SUBROUTINE('elsi_write_mat_real', ['type(elsi_rw_handle)', 'character(*)',&
      & 'real(dp), dimension(:,:)'])

  $:STUB_SUBROUTINE('elsi_write_mat_complex', ['type(elsi_rw_handle)', 'character(*)',&
      & 'complex(dp), dimension(:,:)'])

#:endif

end module elsiiface
