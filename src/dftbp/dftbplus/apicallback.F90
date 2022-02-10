!--------------------------------------------------------------------------------------------------!
!  DFTB+: general package for performing fast atomistic simulations                                !
!  Copyright (C) 2006 - 2022  DFTB+ developers group                                               !
!                                                                                                  !
!  See the LICENSE file for terms of usage and distribution.                                       !
!--------------------------------------------------------------------------------------------------!

!> TODO
module dftbp_dftbplus_apicallback
  use dftbp_common_accuracy, only : dp
  use iso_c_binding
  implicit none

  private
  public :: TAPICallback
  public :: null_apicallback


  abstract interface
    subroutine dm_callback_t(aux_ptr, i_kpoint, i_spin, blacs_descr, data_ptr) bind(c)
      use iso_c_binding
      type(c_ptr), value :: aux_ptr
      integer(c_int), value :: i_kpoint, i_spin
      type(c_ptr), value :: blacs_descr
      type(c_ptr), value :: data_ptr
    end subroutine dm_callback_t
  end interface

  abstract interface
    subroutine hs_callback_t(aux_ptr,blacs_descr, data_ptr) bind(c)
      use iso_c_binding
      type(c_ptr), value :: aux_ptr
      type(c_ptr), value :: blacs_descr
      type(c_ptr), value :: data_ptr
    end subroutine hs_callback_t
  end interface

  
  !> TODO
  type :: TAPICallback
    private
    !> callback
    type(c_funptr) :: dm_callback
    !> callback context
    type(c_ptr) :: dm_aux_ptr

    !> callback
    type(c_funptr) :: s_callback
    !> callback context
    type(c_ptr) :: s_aux_ptr

    !> callback
    type(c_funptr) :: h_callback
    !> callback context
    type(c_ptr) :: h_aux_ptr

  contains
    !> TODO
    procedure :: registerDM => TAPICallback_registerDM
    procedure :: invokeDM_real => TAPICallback_invokeDM_real
    procedure :: invokeDM_cplx => TAPICallback_invokeDM_cplx
    generic :: invokeDM => invokeDM_real, invokeDM_cplx

    procedure :: registerS => TAPICallback_registerS
    procedure :: invokeS_real => TAPICallback_invokeS_real
    procedure :: invokeS_cplx => TAPICallback_invokeS_cplx
    generic :: invokeS => invokeS_real, invokeS_cplx

    procedure :: registerH => TAPICallback_registerH
    procedure :: invokeH_real => TAPICallback_invokeH_real
    procedure :: invokeH_cplx => TAPICallback_invokeH_cplx
    generic :: invokeH => invokeH_real, invokeH_cplx
  end type TAPICallback

  !> Empty TAPICallback. All invokeXX calls do nothing.
  type(TAPICallback) :: null_apicallback


contains

  subroutine TAPICallback_registerDM(this, callback, aux_ptr)
    class(TAPICallback) :: this
    type(c_funptr), value :: callback
    type(c_ptr), value :: aux_ptr
    
    this%dm_callback = callback
    this%dm_aux_ptr = aux_ptr
  end subroutine TAPICallback_registerDM
  

  subroutine TAPICallback_invokeDM_real(this, i_kpoint, i_spin, data_buf, blacs_descr)
    class(TAPICallback) :: this
    integer(c_int), value :: i_kpoint, i_spin
    real(dp), intent(in), target :: data_buf(:,:)
    integer, intent(in), target, optional :: blacs_descr(:)
    
    procedure(dm_callback_t), pointer :: callback_proc
    type(c_ptr) :: blacs_descr_ptr
    type(c_ptr) :: data_ptr
    
    if (.not. c_associated(this%dm_callback)) then
      return
    endif
    
    call c_f_procpointer(this%dm_callback, callback_proc)
    blacs_descr_ptr = merge(c_loc(blacs_descr(1)), c_null_ptr, present(blacs_descr))
    data_ptr = c_loc(data_buf(1,1))

    call callback_proc(this%dm_aux_ptr, i_kpoint, i_spin, blacs_descr_ptr, data_ptr)
  
  end subroutine TAPICallback_invokeDM_real

  subroutine TAPICallback_invokeDM_cplx(this, i_kpoint, i_spin, data_buf, blacs_descr)
    class(TAPICallback) :: this
    integer(c_int), value :: i_kpoint, i_spin
    complex(dp), intent(in), target :: data_buf(:,:)
    integer, intent(in), target, optional :: blacs_descr(:)

    procedure(dm_callback_t), pointer :: callback_proc
    type(c_ptr) :: blacs_descr_ptr
    type(c_ptr) :: data_ptr

    if (.not. c_associated(this%dm_callback)) then
      return
    endif
    
    call c_f_procpointer(this%dm_callback, callback_proc)
    blacs_descr_ptr = merge(c_loc(blacs_descr(1)), c_null_ptr, present(blacs_descr))
    data_ptr = c_loc(data_buf(1,1))

    call callback_proc(this%dm_aux_ptr, i_kpoint, i_spin, blacs_descr_ptr, data_ptr)

  end subroutine TAPICallback_invokeDM_cplx
  
  

  subroutine TAPICallback_registerS(this, callback, aux_ptr)
    class(TAPICallback) :: this
    type(c_funptr), value :: callback
    type(c_ptr), value :: aux_ptr
    
    this%s_callback = callback
    this%s_aux_ptr = aux_ptr
  end subroutine TAPICallback_registerS
  

  subroutine TAPICallback_invokeS_real(this, data_buf, blacs_descr)
    class(TAPICallback) :: this
    real(dp), intent(in), target :: data_buf(:,:)
    integer, intent(in), target, optional :: blacs_descr(:)
    
    procedure(hs_callback_t), pointer :: callback_proc
    type(c_ptr) :: blacs_descr_ptr
    type(c_ptr) :: data_ptr
    
    if (.not. c_associated(this%s_callback)) then
      return
    endif
    
    call c_f_procpointer(this%s_callback, callback_proc)
    blacs_descr_ptr = merge(c_loc(blacs_descr(1)), c_null_ptr, present(blacs_descr))
    data_ptr = c_loc(data_buf(1,1))

    call callback_proc(this%s_aux_ptr, blacs_descr_ptr, data_ptr)
  
  end subroutine TAPICallback_invokeS_real

  subroutine TAPICallback_invokeS_cplx(this, data_buf, blacs_descr)
    class(TAPICallback) :: this
    complex(dp), intent(in), target :: data_buf(:,:)
    integer, intent(in), target, optional :: blacs_descr(:)

    procedure(hs_callback_t), pointer :: callback_proc
    type(c_ptr) :: blacs_descr_ptr
    type(c_ptr) :: data_ptr
    
    if (.not. c_associated(this%s_callback)) then
      return
    endif
    
    call c_f_procpointer(this%s_callback, callback_proc)
    blacs_descr_ptr = merge(c_loc(blacs_descr(1)), c_null_ptr, present(blacs_descr))
    data_ptr = c_loc(data_buf(1,1))

    call callback_proc(this%s_aux_ptr, blacs_descr_ptr, data_ptr)

  end subroutine TAPICallback_invokeS_cplx

  subroutine TAPICallback_registerH(this, callback, aux_ptr)
    class(TAPICallback) :: this
    type(c_funptr), value :: callback
    type(c_ptr), value :: aux_ptr
    
    this%h_callback = callback
    this%h_aux_ptr = aux_ptr
  end subroutine TAPICallback_registerH
  

  subroutine TAPICallback_invokeH_real(this, data_buf, blacs_descr)
    class(TAPICallback) :: this
    real(dp), intent(in), target :: data_buf(:,:)
    integer, intent(in), target, optional :: blacs_descr(:)

    procedure(hs_callback_t), pointer :: callback_proc
    type(c_ptr) :: blacs_descr_ptr
    type(c_ptr) :: data_ptr
    
    if (.not. c_associated(this%h_callback)) then
      return
    endif
    
    call c_f_procpointer(this%h_callback, callback_proc)
    blacs_descr_ptr = merge(c_loc(blacs_descr(1)), c_null_ptr, present(blacs_descr))
    data_ptr = c_loc(data_buf(1,1))

    call callback_proc(this%h_aux_ptr, blacs_descr_ptr, data_ptr)
  
  end subroutine TAPICallback_invokeH_real

  subroutine TAPICallback_invokeH_cplx(this, data_buf, blacs_descr)
    class(TAPICallback) :: this
    complex(dp), intent(in), target :: data_buf(:,:)
    integer, intent(in), target, optional :: blacs_descr(:)

    procedure(hs_callback_t), pointer :: callback_proc
    type(c_ptr) :: blacs_descr_ptr
    type(c_ptr) :: data_ptr
    
    if (.not. c_associated(this%h_callback)) then
      return
    endif
    
    call c_f_procpointer(this%h_callback, callback_proc)
    blacs_descr_ptr = merge(c_loc(blacs_descr(1)), c_null_ptr, present(blacs_descr))
    data_ptr = c_loc(data_buf(1,1))

    call callback_proc(this%h_aux_ptr, blacs_descr_ptr, data_ptr)

  end subroutine TAPICallback_invokeH_cplx

end module dftbp_dftbplus_apicallback
