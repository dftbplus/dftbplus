!--------------------------------------------------------------------------------------------------!
!  DFTB+: general package for performing fast atomistic simulations                                !
!  Copyright (C) 2006 - 2023  DFTB+ developers group                                               !
!                                                                                                  !
!  See the LICENSE file for terms of usage and distribution.                                       !
!--------------------------------------------------------------------------------------------------!

!> Contains types for the callback interface that exports density, overlap, and hamiltonian matrices
module dftbp_dftbplus_apicallback
  use dftbp_common_accuracy, only : dp
  implicit none

  private
  public :: TAPICallback
  public :: TDMHSCallbackFunc


  !> Callback function signature for overlap, or hamiltonian, or density matrix export in square
  !> dense BLACS format.Type of the matrix elements is either double or complex double, depending on
  !> the task (number of k-points)
  !> Total matrix size is NxN, where N - number of basis functions.
  abstract interface
    subroutine TDMHSCallbackFunc(auxObj, iKpoint, iSpin, blacsDescr, dataBufReal, dataBufCplx)
      use dftbp_common_accuracy, only : dp
      !> Pointer to auxilary data that is set when callback is registered. Can be NULL.
      class(*), intent(inout) :: auxObj
      !> 1-based indices of k-point of the matrix
      integer, value :: iKpoint
      !> 1-based indices of spin chanel of the matrix
      integer, value :: iSpin
      !> BLACS descriptor of the matrix. Can be NULL if DFTB+ is built without SCALAPACK support
      integer, intent(in), target, optional :: blacsDescr(:)
      !> Matrix, that can be either real or complex, buffer for real case
      real(dp), intent(inout), target, optional, contiguous :: dataBufReal(:,:)
      !> Matrix, that can be either real or complex, buffer for complex case
      complex(dp), intent(inout), target, optional, contiguous :: dataBufCplx(:,:)
    end subroutine TDMHSCallbackFunc
  end interface


  !> This type encapsulates registering and invocation of callbacks for export of the density,
  !> overlap, and hamiltonian matrices.
  type :: TAPICallback

    !> Callback for density matrix export
    procedure(TDMHSCallbackFunc), nopass, pointer :: dm_callback => null()
    !> Pointer to auxilary data that is set when the density matrix callback is registered. Can be
    !> NULL.
    class(*), pointer :: dmAuxPtr => null()

    !> Flag that signals that the overlap matrix callback is associated with a function
    !> Callback for the overlap matrix export
    procedure(TDMHSCallbackFunc), pointer, nopass :: s_callback => null()
    !> Pointer to auxilary data that is set when the overlap matrix callback is registered. Can be
    !> NULL.
    class(*), pointer :: sAuxPtr => null()

    !> Callback for the hamiltonian matrix export
    procedure(TDMHSCallbackFunc), pointer, nopass :: h_callback => null()
    !> Pointer to auxilary data that is set when the hamiltonian matrix callback is registered. Can
    !> be NULL.
    class(*), pointer :: hAuxPtr => null()

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

contains


  !> Register callback to be invoked on each density matrix evaluation
  subroutine TAPICallback_registerDM(this, callback, auxPtr)

    !> Instance
    class(TAPICallback) :: this

    !> Callback procedure to invoke
    procedure(TDMHSCallbackFunc) :: callback

    !> Auxillary data set by the callback
    class(*), pointer :: auxPtr

    this%dm_callback => callback
    this%dmAuxPtr => auxPtr

  end subroutine TAPICallback_registerDM


  !> Callback for real density matrix
  subroutine TAPICallback_invokeDM_real(this, iKpoint, iSpin, dataBuf, blacsDescr)

    !> Instance of the callback
    class(TAPICallback) :: this

    !> Index of k-point
    integer, value :: iKpoint

    !> Index of spin chanel
    integer, value :: iSpin

    !> Density matrix in dense format
    real(dp), intent(inout), target :: dataBuf(:,:)

    !> Optional BLACS descriptor for the matrix in dataBuf. Not present if SCALAPACK is not
    !> supported
    integer, intent(in), target, optional :: blacsDescr(:)

    if (.not. associated(this%dm_callback)) then
      return
    endif

    call this%dm_callback(this%dmAuxPtr, iKpoint, iSpin, blacsDescr=blacsDescr, dataBufReal=dataBuf)

  end subroutine TAPICallback_invokeDM_real


  !> Callback for complex density matrix
  subroutine TAPICallback_invokeDM_cplx(this, iKpoint, iSpin, dataBuf, blacsDescr)

    !> Instance
    class(TAPICallback) :: this

    !> Index of k-point
    integer, value :: iKpoint

    !> Index of spin chanel
    integer, value :: iSpin

    !> Density matrix in dense format
    complex(dp), intent(inout), target :: dataBuf(:,:)

    !> Optional BLACS descriptor for the matrix in dataBuf. Not present if SCALAPACK is not
    !> supported
    integer, intent(in), target, optional :: blacsDescr(:)

    if (.not. associated(this%dm_callback)) then
      return
    endif

    call this%dm_callback(this%dmAuxPtr, iKpoint, iSpin, blacsDescr=blacsDescr, dataBufCplx=dataBuf)

  end subroutine TAPICallback_invokeDM_cplx


  !> Register an overlap matrix callback
  subroutine TAPICallback_registerS(this, callback, auxPtr)

    !> Instance
    class(TAPICallback) :: this

    !> Callback procedure to invoke
    procedure(TDMHSCallbackFunc) :: callback

    !> Auxillary data set by the callback
    class(*), pointer :: auxPtr

    this%s_callback => callback
    this%sAuxPtr => auxPtr
  end subroutine TAPICallback_registerS


  !> Callback for real overlap matrix
  subroutine TAPICallback_invokeS_real(this, iKpoint, iSpin, dataBuf, blacsDescr)

    !> Instance
    class(TAPICallback) :: this

    !> Index of k-point
    integer, value :: iKpoint

    !> Index of spin chanel
    integer, value :: iSpin

    !> Overlap matrix in dense format
    real(dp), intent(inout), target :: dataBuf(:,:)

    !> Optional BLACS descriptor for the matrix in dataBuf. Not present if SCALAPACK is not
    !> supported
    integer, intent(in), target, optional :: blacsDescr(:)

    if (.not. associated(this%s_callback)) then
      return
    endif

    call this%s_callback(this%sAuxPtr, iKpoint, iSpin, blacsDescr=blacsDescr, dataBufReal=dataBuf)

  end subroutine TAPICallback_invokeS_real


  !> Callback for complex overlap matrix
  subroutine TAPICallback_invokeS_cplx(this, iKpoint, iSpin, dataBuf, blacsDescr)

    !> Instance
    class(TAPICallback) :: this

    !> Indices of k-point and spin chanel
    integer, value :: iKpoint, iSpin

    !> Overlap matrix in dense format
    complex(dp), intent(inout), target :: dataBuf(:,:)

    !> Optional BLACS descriptor for the matrix in dataBuf. Not present if SCALAPACK is not
    !> supported
    integer, intent(in), target, optional :: blacsDescr(:)

    if (.not. associated(this%s_callback)) then
      return
    endif

    call this%s_callback(this%sAuxPtr, iKpoint, iSpin, blacsDescr=blacsDescr, dataBufCplx=dataBuf)

  end subroutine TAPICallback_invokeS_cplx


  !> Register an hamiltonian matrix callback
  subroutine TAPICallback_registerH(this, callback, auxPtr)

    !> Instance
    class(TAPICallback) :: this

    !> Callback procedure to invoke
    procedure(TDMHSCallbackFunc) :: callback

    !> Auxillary data set by the callback
    class(*), pointer :: auxPtr

    this%h_callback => callback
    this%hAuxPtr => auxPtr

  end subroutine TAPICallback_registerH


  !> Callback for real hamiltonian matrix
  subroutine TAPICallback_invokeH_real(this, iKpoint, iSpin, dataBuf, blacsDescr)

    !> Instance
    class(TAPICallback) :: this

    !> Index of k-point
    integer, value :: iKpoint

    !> Index of spin chanel
    integer, value :: iSpin

    !> Hamiltonian matrix in dense format
    real(dp), intent(inout), target :: dataBuf(:,:)

    !> Optional BLACS descriptor for the matrix in dataBuf. Not present if SCALAPACK is not
    !> supported
    integer, intent(in), target, optional :: blacsDescr(:)

    if (.not. associated(this%h_callback)) then
      return
    endif

    call this%h_callback(this%hAuxPtr, iKpoint, iSpin, blacsDescr=blacsDescr, dataBufReal=dataBuf)

  end subroutine TAPICallback_invokeH_real


  !> Hamiltonian complex matrix callback invocation
  subroutine TAPICallback_invokeH_cplx(this, iKpoint, iSpin, dataBuf, blacsDescr)

    !> Instance
    class(TAPICallback) :: this

    !> Index of k-point
    integer, value :: iKpoint

    !> Index of spin chanel
    integer, value :: iSpin

    !> Hamiltonian matrix in dense format
    complex(dp), intent(inout), target :: dataBuf(:,:)

    !> Optional BLACS descriptor for the matrix in dataBuf. Not present if SCALAPACK is not
    !> supported
    integer, intent(in), target, optional :: blacsDescr(:)

    if (.not. associated(this%h_callback)) then
      return
    endif

    call this%h_callback(this%hAuxPtr, iKpoint, iSpin, blacsDescr=blacsDescr, dataBufCplx=dataBuf)

  end subroutine TAPICallback_invokeH_cplx

end module dftbp_dftbplus_apicallback
