!--------------------------------------------------------------------------------------------------!
!  DFTB+: general package for performing fast atomistic simulations                                !
!  Copyright (C) 2006 - 2022  DFTB+ developers group                                               !
!                                                                                                  !
!  See the LICENSE file for terms of usage and distribution.                                       !
!--------------------------------------------------------------------------------------------------!

!> Contains types for the callback interface that exports density matrix, overlap, and hamiltonians
module dftbp_dftbplus_apicallback
  use dftbp_common_accuracy, only : dp
  use iso_c_binding
  implicit none

  private
  public :: TAPICallback
  public :: null_apicallback


  !> Callback function signature for density matrix export in square dense format
  !> DFTB+ would call it after *each* density matrix evaluation. The density matrix is in BLACS dense 
  !> format, with zero lower triangular part, due to it's symmetry. Type of the density matrix 
  !> elements is either double or complex double, depending on the task (number of k-points) 
  !> Total matrix size is NxN, where N - number of basis functions.
  abstract interface
    subroutine dm_callback_t(aux_ptr, i_kpoint, i_spin, blacs_descr, data_ptr) bind(c)
      use iso_c_binding
      !> Pointer to auxilary data that is set when callback is registered. Can be NULL.
      type(c_ptr), value :: aux_ptr
      !> 1-based indices of k-point and spin chanel of the current density matrix
      integer(c_int), value :: i_kpoint, i_spin
      !> BLACS descriptor of the matrix. Can be NULL if DFTB+ is built without SCALAPACK support
      type(c_ptr), value :: blacs_descr
      !> Pointer to the matrix elements, that can be real or complex
      type(c_ptr), value :: data_ptr
    end subroutine dm_callback_t
  end interface

  !> Callback function signature for the overlap or hamiltonian matrices export in square dense format.
  !> DFTB+ would call it after the *first* overlap or hamiltonian evaluation. The matrix is exported in 
  !> BLACS format. Type of the matrix elements is either double or complex double, depending on the 
  !> task (number of k-points) 
  !> Total matrix size is NxN, where N - number of basis functions.
  abstract interface
    subroutine hs_callback_t(aux_ptr,blacs_descr, data_ptr) bind(c)
      use iso_c_binding
      !> Pointer to auxilary data that is set when callback is registered. Can be NULL.
      type(c_ptr), value :: aux_ptr
      !> BLACS descriptor of the matrix. Can be NULL if DFTB+ is built without SCALAPACK support
      type(c_ptr), value :: blacs_descr
      !> Pointer to the matrix elements, that can be real or complex
      type(c_ptr), value :: data_ptr
    end subroutine hs_callback_t
  end interface

  
  !> This type incapsulates registering and invocation of callbacks for export of the density,
  !> overlap, and hamiltonian matrices.
  type :: TAPICallback
    private
    !> Callback for density matrix export
    type(c_funptr) :: dm_callback = c_null_funptr
    !> Pointer to auxilary data that is set when the density matrix callback is registered. Can be NULL.
    type(c_ptr) :: dm_aux_ptr

    !> Callback for the overlap matrix export
    type(c_funptr) :: s_callback = c_null_funptr
    !> Pointer to auxilary data that is set when the overlap matrix callback is registered. Can be NULL.
    type(c_ptr) :: s_aux_ptr

    !> Callback for the hamiltonian matrix export
    type(c_funptr) :: h_callback = c_null_funptr
    !> Pointer to auxilary data that is set when the hamiltonian matrix callback is registered. Can be NULL.
    type(c_ptr) :: h_aux_ptr

  contains
    !> Register callback to be invoked on each density matrix evaluation
    procedure :: registerDM => TAPICallback_registerDM
    !> This function must be invoked on each density matrix evaluation
    procedure :: invokeDM_real => TAPICallback_invokeDM_real
    procedure :: invokeDM_cplx => TAPICallback_invokeDM_cplx
    generic :: invokeDM => invokeDM_real, invokeDM_cplx

    !> Register callback to be invocked on the first overlap matrix evaluation
    procedure :: registerS => TAPICallback_registerS
    !> This function must be invoked on the first overlap matrix evaluation
    procedure :: invokeS_real => TAPICallback_invokeS_real
    procedure :: invokeS_cplx => TAPICallback_invokeS_cplx
    generic :: invokeS => invokeS_real, invokeS_cplx

    !> Register callback to be invocked on the first hamiltonian matrix evaluation
    procedure :: registerH => TAPICallback_registerH
   !> This function must be invoked on the first hamiltonian matrix evaluation
    procedure :: invokeH_real => TAPICallback_invokeH_real
    procedure :: invokeH_cplx => TAPICallback_invokeH_cplx
    generic :: invokeH => invokeH_real, invokeH_cplx
  end type TAPICallback

  !> Empty TAPICallback. Value for unregistered callbacks
  type(TAPICallback) :: null_apicallback


contains
  !> Register callback to be invoked on each density matrix evaluation
  subroutine TAPICallback_registerDM(this, callback, aux_ptr)
    class(TAPICallback) :: this
    type(c_funptr), value :: callback
    type(c_ptr), value :: aux_ptr
    
    this%dm_callback = callback
    this%dm_aux_ptr = aux_ptr
  end subroutine TAPICallback_registerDM
  

  subroutine TAPICallback_invokeDM_real(this, i_kpoint, i_spin, data_buf, blacs_descr)
    class(TAPICallback) :: this
    !> Indices if k-point and spin chanel
    integer(c_int), value :: i_kpoint, i_spin
    !> Density matrix in dense format
    real(dp), intent(in), target :: data_buf(:,:)
    !> Optional BLACS descriptor for the matrix in data_buf. Not present if SCALAPACK is not supported
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
    !> Indices if k-point and spin chanel
    integer(c_int), value :: i_kpoint, i_spin
    !> Density matrix in dense format
    complex(dp), intent(in), target :: data_buf(:,:)
    !> Optional BLACS descriptor for the matrix in data_buf. Not present if SCALAPACK is not supported
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
    !> Overlap matrix in dense format
    real(dp), intent(in), target :: data_buf(:,:)
    !> Optional BLACS descriptor for the matrix in data_buf. Not present if SCALAPACK is not supported
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
    !> Overlap matrix in dense format
    complex(dp), intent(in), target :: data_buf(:,:)
    !> Optional BLACS descriptor for the matrix in data_buf. Not present if SCALAPACK is not supported
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
    !> Hamiltonian matrix in dense format
    real(dp), intent(in), target :: data_buf(:,:)
    !> Optional BLACS descriptor for the matrix in data_buf. Not present if SCALAPACK is not supported
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
    !> Hamiltonian matrix in dense format
    complex(dp), intent(in), target :: data_buf(:,:)
    !> Optional BLACS descriptor for the matrix in data_buf. Not present if SCALAPACK is not supported
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
