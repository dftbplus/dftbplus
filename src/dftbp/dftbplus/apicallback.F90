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
    
  end type TAPICallback



contains

  subroutine TAPICallback_registerDM(this, callback, aux_ptr)
    class(TAPICallback) :: this
    type(c_funptr), value :: callback
    type(c_ptr), value :: aux_ptr
    
    this%dm_callback = callback
    this%dm_aux_ptr = aux_ptr
  end subroutine TAPICallback_registerDM
  

  subroutine TAPICallback_invokeDM_real(this, i_kpoint, i_spin, blacs_descr, data_buf)
    class(TAPICallback) :: this
    integer(c_int), value :: i_kpoint, i_spin
    integer, intent(in), target :: blacs_descr(:)
    real(dp), intent(in), target :: data_buf(:,:)
    
    procedure(dm_callback_t), pointer :: callback_proc
    type(c_ptr) :: blacs_descr_ptr
    type(c_ptr) :: data_ptr
    
    if (.not. c_associated(this%dm_callback)) then
      return
    endif
    
    call c_f_procpointer(this%dm_callback, callback_proc)
    blacs_descr_ptr = c_loc(blacs_descr(1))
    data_ptr = c_loc(data_buf(1,1))

    call callback_proc(this%dm_aux_ptr, i_kpoint, i_spin, blacs_descr_ptr, data_ptr)
  
  end subroutine TAPICallback_invokeDM_real

  subroutine TAPICallback_invokeDM_cplx(this, i_kpoint, i_spin, blacs_descr, data_buf)
    class(TAPICallback) :: this
    integer(c_int), value :: i_kpoint, i_spin
    integer, intent(in), target :: blacs_descr(:)
    complex(dp), intent(in), target :: data_buf(:,:)

    procedure(dm_callback_t), pointer :: callback_proc
    type(c_ptr) :: blacs_descr_ptr
    type(c_ptr) :: data_ptr

    if (.not. c_associated(this%dm_callback)) then
      return
    endif
    
    call c_f_procpointer(this%dm_callback, callback_proc)
    blacs_descr_ptr = c_loc(blacs_descr(1))
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
  

  subroutine TAPICallback_invokeS_real(this, blacs_descr, data_buf)
    class(TAPICallback) :: this
    integer, intent(in), target :: blacs_descr(:)
    real(dp), intent(in), target :: data_buf(:,:)
    
    procedure(hs_callback_t), pointer :: callback_proc
    type(c_ptr) :: blacs_descr_ptr
    type(c_ptr) :: data_ptr
    
    if (.not. c_associated(this%s_callback)) then
      return
    endif
    
    call c_f_procpointer(this%s_callback, callback_proc)
    blacs_descr_ptr = c_loc(blacs_descr(1))
    data_ptr = c_loc(data_buf(1,1))

    call callback_proc(this%s_aux_ptr, blacs_descr_ptr, data_ptr)
  
  end subroutine TAPICallback_invokeS_real

  subroutine TAPICallback_invokeS_cplx(this, blacs_descr, data_buf)
    class(TAPICallback) :: this
    integer, intent(in), target :: blacs_descr(:)
    complex(dp), intent(in), target :: data_buf(:,:)
    
    procedure(hs_callback_t), pointer :: callback_proc
    type(c_ptr) :: blacs_descr_ptr
    type(c_ptr) :: data_ptr
    
    if (.not. c_associated(this%s_callback)) then
      return
    endif
    
    call c_f_procpointer(this%s_callback, callback_proc)
    blacs_descr_ptr = c_loc(blacs_descr(1))
    data_ptr = c_loc(data_buf(1,1))

    call callback_proc(this%s_aux_ptr, blacs_descr_ptr, data_ptr)

  end subroutine TAPICallback_invokeS_cplx

  
end module dftbp_dftbplus_apicallback
