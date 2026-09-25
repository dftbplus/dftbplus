!--------------------------------------------------------------------------------------------------!
!  DFTB+: general package for performing fast atomistic simulations                                !
!  Copyright (C) 2006 - 2025  DFTB+ developers group                                               !
!                                                                                                  !
!  See the LICENSE file for terms of usage and distribution.                                       !
!--------------------------------------------------------------------------------------------------!

#:include 'common.fypp'

!> Contains types for the callback interface that exports density, overlap, and hamiltonian matrices
module dftbp_dftbplus_apicallback
  use dftbp_common_accuracy, only : dp
  use dftbp_type_densedescr, only : TDenseDescr
  implicit none

  private
  public :: TAPICallback, TDMHSCallbackFunc, TSetDMHSCallbackFunc


  !> Callback function signature for overlap, or hamiltonian, or density matrix export in square
  !! dense BLACS format. Type of the matrix elements is either double or complex double, depending
  !! on the task (number of k-points)
  !! Total matrix size is NxN, where N - number of basis functions.
  abstract interface

    !> Interface example for dense BLACS matrices containing various information
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


  !> Callback function signature for overlap, or hamiltonian, or density matrix export in square
  !! dense BLACS format. Type of the matrix elements is either double or complex double, depending
  !! on the task (number of k-points)
  !! Total matrix size is NxN, where N - number of basis functions.
  abstract interface

    !> Interface example for dense BLACS matrices containing various information, returning a status
    !! value
    integer function TSetDMHSCallbackFunc(auxObj, iKpoint, iSpin, blacsDescr, dataBufReal,&
        & dataBufCplx)
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
    end function TSetDMHSCallbackFunc

  end interface


  !> This type encapsulates registering and invocation of callbacks for export of the density,
  !! overlap, and hamiltonian matrices.
  type :: TAPICallback

    !> Callback for density matrix export
    procedure(TDMHSCallbackFunc), nopass, pointer :: dm_callback => null()
    !> Pointer to auxilary data that is set when the density matrix callback is registered. Can be
    !! NULL.
    class(*), pointer :: dmAuxPtr => null()

    !> Flag that signals that the overlap matrix callback is associated with a function
    !! Callback for the overlap matrix export
    procedure(TDMHSCallbackFunc), pointer, nopass :: s_callback => null()
    !> Pointer to auxilary data that is set when the overlap matrix callback is registered. Can be
    !> NULL.
    class(*), pointer :: sAuxPtr => null()

    !> Flag that signals that the overlap matrix callback is associated with a function
    !! Callback for the overlap matrix export
    procedure(TSetDMHSCallbackFunc), pointer, nopass :: set_s_callback => null()
    !> Pointer to auxilary data that is set when the overlap matrix callback is registered. Can be
    !> NULL.
    class(*), pointer :: set_sAuxPtr => null()

    !> Callback for the hamiltonian matrix export
    procedure(TDMHSCallbackFunc), pointer, nopass :: h_callback => null()
    !> Pointer to auxilary data that is set when the hamiltonian matrix callback is registered. Can
    !! be NULL.
    class(*), pointer :: hAuxPtr => null()

    !> Callback for the hamiltonian matrix import
    procedure(TSetDMHSCallbackFunc), pointer, nopass :: set_h_callback => null()
    !> Pointer to auxilary data that is set when the hamiltonian matrix callback is registered. Can
    !! be NULL.
    class(*), pointer :: set_hAuxPtr => null()

  contains

    !> Register callback to be invoked on each density matrix evaluation
    procedure :: registerDM => TAPICallback_registerDM
    !> This function must be invoked on each density matrix evaluation
    procedure :: invokeDM_real => TAPICallback_invokeDM_real
    procedure :: invokeDM_cplx => TAPICallback_invokeDM_cplx
    generic :: invokeDM => invokeDM_real, invokeDM_cplx

    !> Register callback to be invoked on the first overlap matrix evaluation
    procedure :: registerS => TAPICallback_registerS
    !> Register callback to be invoked to import overlap matrix
    procedure :: registerSetS => TAPICallback_registerSetS
    !> This function must be invoked on the first overlap matrix evaluation
    procedure :: invokeS_real => TAPICallback_invokeS_real
    procedure :: invokeS_cplx => TAPICallback_invokeS_cplx
    generic :: invokeS => invokeS_real, invokeS_cplx

    !> Register callback to be invoked on the first hamiltonian matrix evaluation
    procedure :: registerH => TAPICallback_registerH
    !> Register callback to be invoked to import hamiltonian matrix
    procedure :: registerSetH => TAPICallback_registerSetH
    !> This function must be invoked on the first hamiltonian matrix evaluation
    procedure :: invokeH_real => TAPICallback_invokeH_real
    procedure :: invokeH_cplx => TAPICallback_invokeH_cplx
    generic :: invokeH => invokeH_real, invokeH_cplx

    procedure :: invokeHS_real => TAPIinvokeHS_real
    procedure :: invokeHS_cplx => TAPIinvokeHS_cplx
    generic :: invokeHSCallBack => invokeHS_real, invokeHS_cplx

    !> Are the H/S callbacks expected to modify the model?
    procedure :: canAsiChangeTheModel => TAPICallback_canAsiChangeTheModel

  end type TAPICallback


contains

  !> Are the H/S callbacks expected to modify the model?
  function TAPICallback_canAsiChangeTheModel(this) result(canIt)

    !> Instance
    class(TAPICallback) :: this

    !> Can the ASI binding change the state?
    logical :: canIt

    canIt = associated(this%set_s_callback) .or. associated(this%set_h_callback)

  end function TAPICallback_canAsiChangeTheModel


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
    !! supported
    integer, intent(in), target, optional :: blacsDescr(:)

    if (associated(this%dm_callback)) then
      ! Expose the DFTB+ density matrix to the external code (assumes they do not modify it!)
      call this%dm_callback(this%dmAuxPtr, iKpoint, iSpin, blacsDescr=blacsDescr,&
          & dataBufReal=dataBuf)
    end if

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
    !! supported
    integer, intent(in), target, optional :: blacsDescr(:)

    if (associated(this%dm_callback)) then
      ! Expose the DFTB+ density matrix to the external code (assumes they do not modify it!)
      call this%dm_callback(this%dmAuxPtr, iKpoint, iSpin, blacsDescr=blacsDescr,&
          & dataBufCplx=dataBuf)
    end if

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


  !> Register an overlap matrix setting callback
  subroutine TAPICallback_registerSetS(this, callback, auxPtr)

    !> Instance
    class(TAPICallback) :: this

    !> Callback procedure to invoke
    procedure(TSetDMHSCallbackFunc) :: callback

    !> Auxillary data set by the callback
    class(*), pointer :: auxPtr

    this%set_s_callback => callback
    this%set_sAuxPtr => auxPtr

  end subroutine TAPICallback_registerSetS


  !> Callback for real overlap matrix
  subroutine TAPICallback_invokeS_real(this, iKpoint, iSpin, dataBuf, status, blacsDescr)

    !> Instance
    class(TAPICallback) :: this

    !> Index of k-point
    integer, value :: iKpoint

    !> Index of spin chanel
    integer, value :: iSpin

    !> Overlap matrix in dense format
    real(dp), intent(inout), target :: dataBuf(:,:)

    !> Return value from the ASI matrix setting callback, set -1 if not registered, 0 if unchanged
    integer, intent(out) :: status

    !> Optional BLACS descriptor for the matrix in dataBuf. Not present if SCALAPACK is not
    !! supported
    integer, intent(in), target, optional :: blacsDescr(:)

    if (associated(this%s_callback)) then
      ! Expose the DFTB+ S matrix to the external code (assumes they do not modify it!)
      call this%s_callback(this%sAuxPtr, iKpoint, iSpin, blacsDescr=blacsDescr, dataBufReal=dataBuf)
    endif

    status = -1
    if (associated(this%set_s_callback)) then
      ! Allow the external code to (potentially) modify the S matrix, status on return
      status = this%set_s_callback(this%set_sAuxPtr, iKpoint, iSpin, blacsDescr=blacsDescr,&
          & dataBufReal=dataBuf)
    endif

  end subroutine TAPICallback_invokeS_real


  !> Callback for complex overlap matrix
  subroutine TAPICallback_invokeS_cplx(this, iKpoint, iSpin, dataBuf, status, blacsDescr)

    !> Instance
    class(TAPICallback) :: this

    !> Indices of k-point and spin chanel
    integer, value :: iKpoint, iSpin

    !> Overlap matrix in dense format
    complex(dp), intent(inout), target :: dataBuf(:,:)

    !> Return value from the ASI matrix setting callback, set -1 if not registered, 0 if unchanged
    integer, intent(out) :: status

    !> Optional BLACS descriptor for the matrix in dataBuf. Not present if SCALAPACK is not
    !! supported
    integer, intent(in), target, optional :: blacsDescr(:)

    if (associated(this%s_callback)) then
      ! Expose the DFTB+ S matrix to the external code (assumes they do not modify it!)
      call this%s_callback(this%sAuxPtr, iKpoint, iSpin, blacsDescr=blacsDescr, dataBufCplx=dataBuf)
    endif

    status = -1
    if (associated(this%set_s_callback)) then
      ! Expose the DFTB+ S matrix to the external code, possibly with modification happening
      status = this%set_s_callback(this%set_sAuxPtr, iKpoint, iSpin, blacsDescr=blacsDescr,&
          & dataBufCplx=dataBuf)
    endif

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


  !> Register a hamiltonian matrix setting callback
  subroutine TAPICallback_registerSetH(this, callback, auxPtr)

    !> Instance
    class(TAPICallback) :: this

    !> Callback procedure to invoke
    procedure(TSetDMHSCallbackFunc) :: callback

    !> Auxillary data set by the callback
    class(*), pointer :: auxPtr

    this%set_h_callback => callback
    this%set_hAuxPtr => auxPtr

  end subroutine TAPICallback_registerSetH


  !> Callback for real hamiltonian matrix
  subroutine TAPICallback_invokeH_real(this, iKpoint, iSpin, dataBuf, status, blacsDescr)

    !> Instance
    class(TAPICallback) :: this

    !> Index of k-point
    integer, value :: iKpoint

    !> Index of spin chanel
    integer, value :: iSpin

    !> Hamiltonian matrix in dense format
    real(dp), intent(inout), target :: dataBuf(:,:)

    !> Return value from the ASI matrix setting callback, set -1 if not registered, 0 if unchanged
    integer, intent(out) :: status

    !> Optional BLACS descriptor for the matrix in dataBuf. Not present if SCALAPACK is not
    !! supported
    integer, intent(in), target, optional :: blacsDescr(:)

    if (associated(this%h_callback)) then
      ! Expose the DFTB+ H matrix to the external code (assumes they do not modify it!)
      call this%h_callback(this%hAuxPtr, iKpoint, iSpin, blacsDescr=blacsDescr, dataBufReal=dataBuf)
    endif

    status = -1
    if (associated(this%set_h_callback)) then
      ! Expose the DFTB+ H matrix to the external code, possibly with modification happening
      status = this%set_h_callback(this%set_hAuxPtr, iKpoint, iSpin, blacsDescr=blacsDescr,&
          & dataBufReal=dataBuf)
    endif

  end subroutine TAPICallback_invokeH_real


  !> Hamiltonian complex matrix callback invocation
  subroutine TAPICallback_invokeH_cplx(this, iKpoint, iSpin, dataBuf, status, blacsDescr)

    !> Instance
    class(TAPICallback) :: this

    !> Index of k-point
    integer, value :: iKpoint

    !> Index of spin chanel
    integer, value :: iSpin

    !> Hamiltonian matrix in dense format
    complex(dp), intent(inout), target :: dataBuf(:,:)

    !> Return value from the ASI matrix setting callback, set -1 if not registered, 0 if unchanged
    integer, intent(out) :: status

    !> Optional BLACS descriptor for the matrix in dataBuf. Not present if SCALAPACK is not
    !! supported
    integer, intent(in), target, optional :: blacsDescr(:)

    if (associated(this%h_callback)) then
      ! Expose the DFTB+ H matrix to the external code (assumes they do not modify it!)
      call this%h_callback(this%hAuxPtr, iKpoint, iSpin, blacsDescr=blacsDescr, dataBufCplx=dataBuf)
    endif

    status = -1
    if (associated(this%set_h_callback)) then
      ! Expose the DFTB+ H matrix to the external code, possibly with modification happening
      status = this%set_h_callback(this%set_hAuxPtr, iKpoint, iSpin, blacsDescr=blacsDescr,&
          & dataBufCplx=dataBuf)
    endif

  end subroutine TAPICallback_invokeH_cplx


#:if WITH_SCALAPACK

  !> Callback invocation
  subroutine TAPIinvokeHS_real(this, iKS, isCholesky, S, H, isSChanged, isHChanged, denseDesc)

    !> Instance
    class(TAPICallback) :: this

    !> The combined spin and k-point index
    integer, intent(in) :: iKS(:)

    !> Is overlap already Cholesky factorized?
    logical, intent(in) :: isCholesky

    !> Dense overlap matrix
    real(dp), intent(inout) :: S(:,:)

    !> Dense hamiltonian matrix
    real(dp), intent(out) :: H(:,:)

    !> Is the overlap modified?
    logical, intent(out) :: isSChanged

    !> Is the overlap modified?
    logical, intent(out) :: isHChanged

    !> Dense matrix descriptor
    type(TDenseDescr), intent(in) :: denseDesc

    integer :: iK, iS, status

    iK = iKS(1)
    iS = iKS(2)
    status = 0
    if (.not.isCholesky) then
      call this%invokeS(iK, iS, S, status, denseDesc%blacsOrbSqr)
    end if
    isSChanged = status > 0
    call this%invokeH(iK, iS, H, status, denseDesc%blacsOrbSqr)
    isHChanged = status > 0

  end subroutine TAPIinvokeHS_real


  !> Callback invocation
  subroutine TAPIinvokeHS_cplx(this, iKS, isCholesky, S, H, isSChanged, isHChanged, denseDesc)

    !> Instance
    class(TAPICallback) :: this

    !> The combined spin and k-point index
    integer, intent(in) :: iKS(:)

    !> Is overlap already Cholesky factorized?
    logical, intent(in) :: isCholesky

    !> Dense overlap matrix
    complex(dp), intent(inout) :: S(:,:)

    !> Dense hamiltonian matrix
    complex(dp), intent(out) :: H(:,:)

    !> Dense matrix descriptor
    type(TDenseDescr), intent(in) :: denseDesc

    !> Is the overlap modified?
    logical, intent(out) :: isSChanged

    !> Is the overlap modified?
    logical, intent(out) :: isHChanged

    integer :: iK, iS, status

    iK = iKS(1)
    iS = iKS(2)
    status = 0
    if (.not.isCholesky) then
      call this%invokeS(iK, iS, S, status, denseDesc%blacsOrbSqr)
    end if
    isSChanged = status > 0
    call this%invokeH(iK, iS, H, status, denseDesc%blacsOrbSqr)
    isHChanged = status > 0

  end subroutine TAPIinvokeHS_cplx


#:else

  !> Callback invocation
  subroutine TAPIinvokeHS_real(this, iKS, isCholesky, S, H, isSChanged, isHChanged)

    !> Instance
    class(TAPICallback) :: this

    !> The combined spin and k-point index
    integer, intent(in) :: iKS(:)

    !> Is overlap already Cholesky factorized?
    logical, intent(in) :: isCholesky

    !> Dense overlap matrix
    real(dp), intent(inout) :: S(:,:)

    !> Dense hamiltonian matrix
    real(dp), intent(out) :: H(:,:)

    !> Is the overlap modified?
    logical, intent(out) :: isSChanged

    !> Is the overlap modified?
    logical, intent(out) :: isHChanged

    integer :: iK, iS, status

    iK = iKS(1)
    iS = iKS(2)
    status = 0
    if (.not.isCholesky) then
      call this%invokeS(iK, iS, S, status)
    end if
    isSChanged = status > 0
    call this%invokeH(iK, iS, H, status)
    isHChanged = status > 0

  end subroutine TAPIinvokeHS_real


  !> Callback invocation
  subroutine TAPIinvokeHS_cplx(this, iKS, isCholesky, S, H, isSChanged, isHChanged)

    !> Instance
    class(TAPICallback) :: this

    !> The combined spin and k-point index
    integer, intent(in) :: iKS(:)

    !> Is overlap already Cholesky factorized?
    logical, intent(in) :: isCholesky

    !> Dense overlap matrix
    complex(dp), intent(inout) :: S(:,:)

    !> Dense hamiltonian matrix
    complex(dp), intent(out) :: H(:,:)

    !> Is the overlap modified?
    logical, intent(out) :: isSChanged

    !> Is the overlap modified?
    logical, intent(out) :: isHChanged

    integer :: iK, iS, status

    iK = iKS(1)
    iS = iKS(2)
    status = 0
    if (.not.isCholesky) then
      call this%invokeS(iK, iS, S, status)
    end if
    isSChanged = status > 0
    call this%invokeH(iK, iS, H, status)
    isHChanged = status > 0

  end subroutine TAPIinvokeHS_cplx

#:endif

end module dftbp_dftbplus_apicallback
