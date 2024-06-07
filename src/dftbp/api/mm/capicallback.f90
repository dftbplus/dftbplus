!--------------------------------------------------------------------------------------------------!
!  DFTB+: general package for performing fast atomistic simulations                                !
!  Copyright (C) 2006 - 2024  DFTB+ developers group                                               !
!                                                                                                  !
!  See the LICENSE file for terms of usage and distribution.                                       !
!--------------------------------------------------------------------------------------------------!

#:include 'common.fypp'

!> Contains types for the callback interface that exports dense matrix formatted density matrix,
!! overlap and hamiltonian data
module dftbp_capicallback
  use iso_c_binding, only: c_ptr, c_int, c_funptr, c_loc, c_null_ptr, c_f_procpointer
  use dftbp_common_accuracy, only : dp

  implicit none

  private
  public :: dmhs_callback_c_wrapper, set_dmhs_callback_c_wrapper, TCAuxWrapper, TCMatrixDescr

  !> Encapsulates C pointers to a callback and auxilary parameter for it
  type :: TCAuxWrapper
    !> Auxilary parameter for callback
    type(c_ptr) :: auxPtr
    !> Pointer to a callback function
    type(c_funptr) :: callback
  end type TCAuxWrapper

  !> Descriptor for matrix structure and storage types
  !> See ASI_matrix_descr_t in dftbplus.h
  type :: TCMatrixDescr
    !> Type of matrix: 0 - generic; 1 - hermitian or symmetric
    integer(c_int) :: matrix_type
    !> Type of storage: 0 - dense full; 1 - only lower triangle (upper for row-major); 2 - only upper triangle (lower for row-major).
    integer(c_int) :: storage_type
  end type TCMatrixDescr


  !> C-style wrapper for TDMHSCallbackFunc
  abstract interface
    subroutine dmhs_callback_c_t(auxPtr, iKpoint, iSpin, blacsDescr, dataPtr, matrixDescr) bind(c)
      use iso_c_binding, only: c_int, c_ptr
      !> Pointer to auxilary data that is set when callback is registered. Can be NULL.
      type(c_ptr), value :: auxPtr
      !> 1-based index of k-points of the matrix
      integer(c_int), value :: iKpoint
      !> 1-based index of spin chanel of the matrix
      integer(c_int), value :: iSpin
      !> BLACS descriptor of the matrix. Can be NULL if DFTB+ is built without SCALAPACK support
      type(c_ptr), value :: blacsDescr
      !> Pointer to the matrix elements, that can be real or complex
      type(c_ptr), value :: dataPtr
      !> Pointer to the descriptor of matrix structure and storage types
      type(c_ptr), value :: matrixDescr
    end subroutine dmhs_callback_c_t
  end interface

  !> C-style wrapper for TSetDMHSCallbackFunc
  abstract interface
    integer function set_dmhs_callback_c_t(auxPtr, iKpoint, iSpin, blacsDescr, dataPtr, matrixDescr) bind(c)
      use iso_c_binding, only: c_int, c_ptr
      !> Pointer to auxilary data that is set when callback is registered. Can be NULL.
      type(c_ptr), value :: auxPtr
      !> 1-based index of k-points of the matrix
      integer(c_int), value :: iKpoint
      !> 1-based index of spin chanel of the matrix
      integer(c_int), value :: iSpin
      !> BLACS descriptor of the matrix. Can be NULL if DFTB+ is built without SCALAPACK support
      type(c_ptr), value :: blacsDescr
      !> Pointer to the matrix elements, that can be real or complex
      type(c_ptr), value :: dataPtr
      !> Pointer to the descriptor of matrix structure and storage types
      type(c_ptr), value :: matrixDescr
    end function set_dmhs_callback_c_t
  end interface

contains

  !> Wrapper for C callbacks to be invoked for matrix export.
  subroutine  dmhs_callback_c_wrapper(auxObj, iKpoint, iSpin, blacsDescr, dataBufReal, &
      & dataBufCplx)

    !> Pointer to auxilary data that is set when callback is registered. Can be NULL.
    class(*), intent(inout) :: auxObj

    !> 1-based index of k-points of the matrix
    integer, value :: iKpoint

    !> 1-based index of spin chanel of the matrix
    integer, value :: iSpin

    !> BLACS descriptor of the matrix. Can be NULL if DFTB+ is built without SCALAPACK support
    integer, intent(in), target, optional :: blacsDescr(:)

    !> Matrix data buffer, for real case
    real(dp), intent(inout), target, optional, contiguous :: dataBufReal(:,:)

    !> Matrix data buffer, for complex case
    complex(dp), intent(inout), target, optional, contiguous :: dataBufCplx(:,:)

    procedure(dmhs_callback_c_t), pointer :: callbackProc
    type(c_ptr) :: blacsDescrPtr
    type(c_ptr) :: dataPtr
    type(TCMatrixDescr), target :: mdescr
    
    mdescr%matrix_type = 1  ! symmetric / hermitian
    #:if WITH_SCALAPACK
      mdescr%storage_type = 0 ! dense full
    #:else
      mdescr%storage_type = 1 ! lower triangle, unpacked
    #:endif

    if (present(blacsDescr)) then
      blacsDescrPtr = c_loc(blacsDescr(1))
    else
      blacsDescrPtr = c_null_ptr
    endif

    if (present(dataBufReal)) then
      dataPtr = c_loc(dataBufReal(1,1))
    else
      dataPtr = c_loc(dataBufCplx(1,1))
    endif
    select type(auxObj)
    type is (TCAuxWrapper)
      call c_f_procpointer(auxObj%callback, callbackProc)
      call callbackProc(auxObj%auxPtr, iKpoint, iSpin, blacsDescrPtr, dataPtr, c_loc(mdescr))
    end select

  end subroutine dmhs_callback_c_wrapper

  !> Wrapper for C callbacks to be invoked for matrix import.
  integer function set_dmhs_callback_c_wrapper(auxObj, iKpoint, iSpin, blacsDescr, dataBufReal, &
      & dataBufCplx)

    !> Pointer to auxilary data that is set when callback is registered. Can be NULL.
    class(*), intent(inout) :: auxObj

    !> 1-based index of k-points of the matrix
    integer, value :: iKpoint

    !> 1-based index of spin chanel of the matrix
    integer, value :: iSpin

    !> BLACS descriptor of the matrix. Can be NULL if DFTB+ is built without SCALAPACK support
    integer, intent(in), target, optional :: blacsDescr(:)

    !> Matrix data buffer, for real case
    real(dp), intent(inout), target, optional, contiguous :: dataBufReal(:,:)

    !> Matrix data buffer, for complex case
    complex(dp), intent(inout), target, optional, contiguous :: dataBufCplx(:,:)

    procedure(set_dmhs_callback_c_t), pointer :: callbackProc
    type(c_ptr) :: blacsDescrPtr
    type(c_ptr) :: dataPtr
    type(TCMatrixDescr), target :: mdescr
    
    #:if WITH_SCALAPACK
      mdescr%storage_type = 0 ! dense full
    #:else
      mdescr%storage_type = 1 ! lower triangle, unpacked
    #:endif

    if (present(blacsDescr)) then
      blacsDescrPtr = c_loc(blacsDescr(1))
    else
      blacsDescrPtr = c_null_ptr
    endif

    if (present(dataBufReal)) then
      dataPtr = c_loc(dataBufReal(1,1))
    else
      dataPtr = c_loc(dataBufCplx(1,1))
    endif
    select type(auxObj)
    type is (TCAuxWrapper)
      call c_f_procpointer(auxObj%callback, callbackProc)
      set_dmhs_callback_c_wrapper = callbackProc(auxObj%auxPtr, iKpoint, iSpin, blacsDescrPtr, dataPtr, c_loc(mdescr))
    end select

  end function set_dmhs_callback_c_wrapper

end module dftbp_capicallback
