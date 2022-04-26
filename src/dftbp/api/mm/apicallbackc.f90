!--------------------------------------------------------------------------------------------------!
!  DFTB+: general package for performing fast atomistic simulations                                !
!  Copyright (C) 2006 - 2022  DFTB+ developers group                                               !
!                                                                                                  !
!  See the LICENSE file for terms of usage and distribution.                                       !
!--------------------------------------------------------------------------------------------------!

!> Contains types for the callback interface that exports density matrix, overlap, and hamiltonians
module dftbp_apicallbackc
  use iso_c_binding
  use dftbp_dftbplus_apicallback, only: dmhs_callback_t

  implicit none

  private
  public :: dmhs_callback_c_wrapper_ptr, TCAuxWrapper

  type :: TCAuxWrapper
    type(c_ptr) :: aux_ptr
    type(c_funptr) :: callback
  end type TCAuxWrapper

  !> C-style wrapper for dmhs_callback_t from apicallback.f90 (see for details).
  abstract interface
    subroutine dmhs_callback_c_t(aux_ptr, i_kpoint, i_spin, blacs_descr, data_ptr) bind(c)
      use iso_c_binding
      !> Pointer to auxilary data that is set when callback is registered. Can be NULL.
      type(c_ptr), value :: aux_ptr
      !> 1-based indices of k-point and spin chanel of the matrix
      integer(c_int), value :: i_kpoint, i_spin
      !> BLACS descriptor of the matrix. Can be NULL if DFTB+ is built without SCALAPACK support
      type(c_ptr), value :: blacs_descr
      !> Pointer to the matrix elements, that can be real or complex
      type(c_ptr), value :: data_ptr
    end subroutine dmhs_callback_c_t
  end interface
  
  procedure(dmhs_callback_t), pointer:: dmhs_callback_c_wrapper_ptr => dmhs_callback_c_wrapper

contains

  !> Register callback to be invoked on each density matrix evaluation. See apicallback.f90 for details.
  subroutine  dmhs_callback_c_wrapper(aux_obj, i_kpoint, i_spin, blacs_descr, data_buf_real, data_buf_cplx)
    use dftbp_common_accuracy, only : dp
    use iso_c_binding

    !> Pointer to auxilary data that is set when callback is registered. Can be NULL.
    class(*), intent(inout) :: aux_obj
    !> 1-based indices of k-point and spin chanel of the matrix
    integer, value :: i_kpoint, i_spin
    !> BLACS descriptor of the matrix. Can be NULL if DFTB+ is built without SCALAPACK support
    integer, intent(in), target, optional :: blacs_descr(:)
    !> Matrix, that can be either real or complex
    real(dp),    intent(in), target, optional :: data_buf_real(:,:)
    complex(dp), intent(in), target, optional :: data_buf_cplx(:,:)
    
    procedure(dmhs_callback_c_t), pointer :: callback_proc
    type(c_ptr) :: blacs_descr_ptr
    type(c_ptr) :: data_ptr

    if (present(blacs_descr)) then
      blacs_descr_ptr = c_loc(blacs_descr(1))
    else
      blacs_descr_ptr = c_null_ptr
    endif
    
    if (present(data_buf_real)) then
      data_ptr = c_loc(data_buf_real(1,1))
    else
      data_ptr = c_loc(data_buf_cplx(1,1))
    endif
    select type(aux_obj)
    type is (TCAuxWrapper)
      call c_f_procpointer(aux_obj%callback, callback_proc)
      call callback_proc(aux_obj%aux_ptr, i_kpoint, i_spin, blacs_descr_ptr, data_ptr)
    end select
  end subroutine dmhs_callback_c_wrapper
  

end module dftbp_apicallbackc
