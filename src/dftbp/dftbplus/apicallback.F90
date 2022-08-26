!--------------------------------------------------------------------------------------------------!
!  DFTB+: general package for performing fast atomistic simulations                                !
!  Copyright (C) 2006 - 2022  DFTB+ developers group                                               !
!                                                                                                  !
!  See the LICENSE file for terms of usage and distribution.                                       !
!--------------------------------------------------------------------------------------------------!

!> Contains types for the callback interface that exports density matrix, overlap, and hamiltonians
module dftbp_dftbplus_apicallback
  use dftbp_common_accuracy, only : dp
  implicit none

  private
  public :: TAPICallback
  public :: null_apicallback, dmhs_callback_t


  !> Callback function signature for overlap, or hamiltonian, or density matrix export in square 
  !> dense BLACS format.Type of the matrix elements is either double or complex double, depending on 
  !> the task (number of k-points) 
  !> Total matrix size is NxN, where N - number of basis functions.
  abstract interface
    subroutine dmhs_callback_t(aux_obj, i_kpoint, i_spin, blacs_descr, data_buf_real, data_buf_cplx)
      use dftbp_common_accuracy, only : dp
      !> Pointer to auxilary data that is set when callback is registered. Can be NULL.
      class(*), intent(inout) :: aux_obj
      !> 1-based indices of k-point and spin chanel of the matrix
      integer, value :: i_kpoint, i_spin
      !> BLACS descriptor of the matrix. Can be NULL if DFTB+ is built without SCALAPACK support
      integer, intent(in), target, optional :: blacs_descr(:)
      !> Matrix, that can be either real or complex
      real(dp),    intent(inout), target, optional :: data_buf_real(:,:)
      complex(dp), intent(inout), target, optional :: data_buf_cplx(:,:)
    end subroutine dmhs_callback_t
  end interface

  
  !> This type incapsulates registering and invocation of callbacks for export of the density,
  !> overlap, and hamiltonian matrices.
  type :: TAPICallback

    !> Flag that signals that the density matrix callback is associated with a function
    logical :: dm_callback_associated = .false.
    !> Callback for density matrix export
    procedure(dmhs_callback_t), nopass, pointer :: dm_callback
    !> Pointer to auxilary data that is set when the density matrix callback is registered. Can be NULL.
    class(*), pointer :: dm_aux_ptr

    !> Flag that signals that the overlap matrix callback is associated with a function
    logical :: s_callback_associated = .false.
    !> Callback for the overlap matrix export
    procedure(dmhs_callback_t), pointer, nopass :: s_callback
    !> Pointer to auxilary data that is set when the overlap matrix callback is registered. Can be NULL.
    class(*), pointer :: s_aux_ptr

    !> Flag that signals that the hamiltonian matrix callback is associated with a function
    logical :: h_callback_associated = .false.
    !> Callback for the hamiltonian matrix export
    procedure(dmhs_callback_t), pointer, nopass :: h_callback
    !> Pointer to auxilary data that is set when the hamiltonian matrix callback is registered. Can be NULL.
    class(*), pointer :: h_aux_ptr

    !> Number (index) of the current self-consistent charge iteration. Meant to be used by the calling code
    integer :: iSCCIter

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
    !> Instance
    class(TAPICallback) :: this

    procedure(dmhs_callback_t), pointer, intent(in) :: callback
    class(*), pointer :: aux_ptr
    
    this%dm_callback => callback
    this%dm_aux_ptr => aux_ptr
    this%dm_callback_associated = .true.
  end subroutine TAPICallback_registerDM
  

  subroutine TAPICallback_invokeDM_real(this, i_kpoint, i_spin, data_buf, blacs_descr)
    class(TAPICallback) :: this
    !> Indices of k-point and spin chanel
    integer, value :: i_kpoint, i_spin
    !> Density matrix in dense format
    real(dp), intent(inout), target :: data_buf(:,:)
    !> Optional BLACS descriptor for the matrix in data_buf. Not present if SCALAPACK is not supported
    integer, intent(in), target, optional :: blacs_descr(:)

    if (.not. this%dm_callback_associated) then
      return
    endif

    if (present(blacs_descr)) then
      call this%dm_callback(this%dm_aux_ptr, i_kpoint, i_spin, blacs_descr, data_buf_real=data_buf)
    else
      call this%dm_callback(this%dm_aux_ptr, i_kpoint, i_spin, data_buf_real=data_buf)
    endif

  end subroutine TAPICallback_invokeDM_real

  subroutine TAPICallback_invokeDM_cplx(this, i_kpoint, i_spin, data_buf, blacs_descr)
    class(TAPICallback) :: this
    !> Indices of k-point and spin chanel
    integer, value :: i_kpoint, i_spin
    !> Density matrix in dense format
    complex(dp), intent(inout), target :: data_buf(:,:)
    !> Optional BLACS descriptor for the matrix in data_buf. Not present if SCALAPACK is not supported
    integer, intent(in), target, optional :: blacs_descr(:)

    if (.not. this%dm_callback_associated) then
      return
    endif
    
    if (present(blacs_descr)) then
      call this%dm_callback(this%dm_aux_ptr, i_kpoint, i_spin, blacs_descr, data_buf_cplx=data_buf)
    else
      call this%dm_callback(this%dm_aux_ptr, i_kpoint, i_spin, data_buf_cplx=data_buf)
    endif

  end subroutine TAPICallback_invokeDM_cplx
  

  subroutine TAPICallback_registerS(this, callback, aux_ptr)
    class(TAPICallback) :: this
    procedure(dmhs_callback_t), pointer :: callback
    class(*), pointer :: aux_ptr
    
    this%s_callback => callback
    this%s_aux_ptr => aux_ptr
    this%s_callback_associated = .true.
  end subroutine TAPICallback_registerS
  

  subroutine TAPICallback_invokeS_real(this, i_kpoint, i_spin, data_buf, blacs_descr)
    class(TAPICallback) :: this
    !> Indices of k-point and spin chanel
    integer, value :: i_kpoint, i_spin
    !> Overlap matrix in dense format
    real(dp), intent(inout), target :: data_buf(:,:)
    !> Optional BLACS descriptor for the matrix in data_buf. Not present if SCALAPACK is not supported
    integer, intent(in), target, optional :: blacs_descr(:)
    
    if (.not. this%s_callback_associated) then
      return
    endif

    if (present(blacs_descr)) then
      call this%s_callback(this%s_aux_ptr, i_kpoint, i_spin, blacs_descr, data_buf_real=data_buf)
    else
      call this%s_callback(this%s_aux_ptr, i_kpoint, i_spin, data_buf_real=data_buf)
    endif
  
  end subroutine TAPICallback_invokeS_real

  subroutine TAPICallback_invokeS_cplx(this, i_kpoint, i_spin, data_buf, blacs_descr)
    class(TAPICallback) :: this
    !> Indices of k-point and spin chanel
    integer, value :: i_kpoint, i_spin
    !> Overlap matrix in dense format
    complex(dp), intent(inout), target :: data_buf(:,:)
    !> Optional BLACS descriptor for the matrix in data_buf. Not present if SCALAPACK is not supported
    integer, intent(in), target, optional :: blacs_descr(:)

    if (.not. this%s_callback_associated) then
      return
    endif

    if (present(blacs_descr)) then
      call this%s_callback(this%s_aux_ptr, i_kpoint, i_spin, blacs_descr, data_buf_cplx=data_buf)
    else
      call this%s_callback(this%s_aux_ptr, i_kpoint, i_spin, data_buf_cplx=data_buf)
    endif

  end subroutine TAPICallback_invokeS_cplx

  subroutine TAPICallback_registerH(this, callback, aux_ptr)
    class(TAPICallback) :: this
    procedure(dmhs_callback_t), pointer :: callback
    class(*), pointer :: aux_ptr
    
    this%h_callback => callback
    this%h_aux_ptr => aux_ptr
    this%h_callback_associated = .true.
  end subroutine TAPICallback_registerH
  

  subroutine TAPICallback_invokeH_real(this, i_kpoint, i_spin, data_buf, blacs_descr)
    class(TAPICallback) :: this
    !> Indices of k-point and spin chanel
    integer, value :: i_kpoint, i_spin
    !> Hamiltonian matrix in dense format
    real(dp), intent(inout), target :: data_buf(:,:)
    !> Optional BLACS descriptor for the matrix in data_buf. Not present if SCALAPACK is not supported
    integer, intent(in), target, optional :: blacs_descr(:)

    if (.not. this%h_callback_associated) then
      return
    endif

    if (present(blacs_descr)) then
      call this%h_callback(this%h_aux_ptr, i_kpoint, i_spin, blacs_descr, data_buf_real=data_buf)
    else
      call this%h_callback(this%h_aux_ptr, i_kpoint, i_spin, data_buf_real=data_buf)
    endif
  end subroutine TAPICallback_invokeH_real

  subroutine TAPICallback_invokeH_cplx(this, i_kpoint, i_spin, data_buf, blacs_descr)
    class(TAPICallback) :: this
    !> Indices of k-point and spin chanel
    integer, value :: i_kpoint, i_spin
    !> Hamiltonian matrix in dense format
    complex(dp), intent(inout), target :: data_buf(:,:)
    !> Optional BLACS descriptor for the matrix in data_buf. Not present if SCALAPACK is not supported
    integer, intent(in), target, optional :: blacs_descr(:)

    if (.not. this%h_callback_associated) then
      return
    endif
    
    if (present(blacs_descr)) then
      call this%h_callback(this%h_aux_ptr, i_kpoint, i_spin, blacs_descr, data_buf_cplx=data_buf)
    else
      call this%h_callback(this%h_aux_ptr, i_kpoint, i_spin, data_buf_cplx=data_buf)
    endif
  end subroutine TAPICallback_invokeH_cplx

end module dftbp_dftbplus_apicallback
