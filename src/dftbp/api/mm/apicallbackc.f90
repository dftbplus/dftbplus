!--------------------------------------------------------------------------------------------------!
!  DFTB+: general package for performing fast atomistic simulations                                !
!  Copyright (C) 2006 - 2023  DFTB+ developers group                                               !
!                                                                                                  !
!  See the LICENSE file for terms of usage and distribution.                                       !
!--------------------------------------------------------------------------------------------------!

!> Contains types for the callback interface that exports dense matrix formatted density matrix,
!> overlap and hamiltonian data
module dftbp_apicallbackc
  use iso_c_binding, only: c_ptr, c_funptr, c_loc, c_null_ptr, c_f_procpointer
  use dftbp_dftbplus_apicallback, only: TDMHSCallbackFunc
  use dftbp_common_accuracy, only : dp

  implicit none

  private
  public :: dmhs_callback_c_wrapper, TCAuxWrapper

  !> Encapsulates C pointers to a callback and auxilary parameter for it
  type :: TCAuxWrapper
    !> Auxilary parameter for callback
    type(c_ptr) :: auxPtr
    !> Pointer to a callback function
    type(c_funptr) :: callback
  end type TCAuxWrapper

  !> C-style wrapper for TDMHSCallbackFunc
  abstract interface
    subroutine dmhs_callback_c_t(auxPtr, iKpoint, iSpin, blacsDescr, dataPtr) bind(c)
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
    end subroutine dmhs_callback_c_t
  end interface

contains

  !> Register callback to be invoked on each density matrix evaluation.
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
      call callbackProc(auxObj%auxPtr, iKpoint, iSpin, blacsDescrPtr, dataPtr)
    end select

  end subroutine dmhs_callback_c_wrapper
  

end module dftbp_apicallbackc
