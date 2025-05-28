!--------------------------------------------------------------------------------------------------!
!  DFTB+: general package for performing fast atomistic simulations                                !
!  Copyright (C) 2006 - 2025  DFTB+ developers group                                               !
!                                                                                                  !
!  See the LICENSE file for terms of usage and distribution.                                       !
!--------------------------------------------------------------------------------------------------!

#:include 'common.fypp'
#:include 'error.fypp'

!> Module to wrap around the process of converting from a Hamiltonian and overlap in sparse form
!> into eigenvectors
module dftbp_dftbplus_eigenvects
  use dftbp_common_accuracy, only : dp
  use dftbp_common_environment, only : TEnvironment
  use dftbp_common_status, only : TStatus
  use dftbp_elecsolvers_elecsolvers, only : electronicSolverTypes, TElectronicSolver
  use dftbp_extlibs_elsiiface, only : elsi_ev_complex, elsi_ev_real, elsi_finalize_rw,&
      & elsi_write_mat_complex, elsi_write_mat_real
  use dftbp_io_message, only : cleanShutdown
  use dftbp_math_eigensolver, only : hegv, hegvd, hegvr
#:if WITH_MAGMA
  use dftbp_math_eigensolver, only : magmaHegvd
#:endif
#:if WITH_SCALAPACK
  use dftbp_extlibs_scalapackfx, only : CSRC_, DLEN_, MB_, NB_, RSRC_, scalafx_indxl2g,&
      & scalafx_phegv, scalafx_phegvd, scalafx_phegvr, scalafx_psygv, scalafx_psygvd,&
      & scalafx_psygvr
#:endif
  implicit none

  private
  public :: diagDenseMtx
#:if WITH_SCALAPACK
  public :: diagDenseMtxBlacs
#:endif


  interface diagDenseMtx
    module procedure diagDenseRealMtx
    module procedure diagDenseComplexMtx
  end interface diagDenseMtx

#:if WITH_SCALAPACK

  interface diagDenseMtxBlacs
    module procedure diagDenseRealMtxBlacs
    module procedure diagDenseCplxMtxBlacs
  end interface diagDenseMtxBlacs

#:endif

contains

  !> Diagonalizes a sparse represented Hamiltonian and overlap to give the eigenvectors and values,
  !> as well as often the Cholesky factorized overlap matrix (due to a side effect of lapack)
  subroutine diagDenseRealMtx(env, electronicSolver, jobz, HSqrReal, SSqrReal, eigen, errStatus)

    !> Environment
    type(TEnvironment), intent(in) :: env

    !> Electronic solver information
    type(TElectronicSolver), intent(inout) :: electronicSolver

    !> type of eigen-problem, either 'V'/'v' with vectors or 'N'/'n' eigenvalues only
    character, intent(in) :: jobz

    !> Large square matrix for the resulting eigenvectors
    real(dp), intent(inout) :: HSqrReal(:,:)

    !> Large square matrix for the overlap workspace, often overwritten with the Cholesky factorized
    !> form.
    real(dp), intent(inout) :: SSqrReal(:,:)

    !> Eigen values.
    real(dp), intent(out) :: eigen(:)

    !> Status of routine
    type(TStatus), intent(out) :: errStatus

    @:ASSERT(size(HSqrReal, dim=1) == size(HSqrReal, dim=2))
    @:ASSERT(all(shape(HSqrReal) == shape(SSqrReal)))
    @:ASSERT(size(HSqrReal, dim=1) == size(eigen))
    @:ASSERT(jobz == 'n' .or. jobz == 'N' .or. jobz == 'v' .or. jobz == 'V')
    select case(electronicSolver%iSolver)
    case(electronicSolverTypes%QR)
      call hegv(HSqrReal,SSqrReal,eigen,'L',jobz)
    case(electronicSolverTypes%divideandconquer)
      call hegvd(HSqrReal,SSqrReal,eigen,'L',jobz)
    case(electronicSolverTypes%relativelyrobust)
      call hegvr(HSqrReal,SSqrReal,eigen,'L',jobz)
    case(electronicSolverTypes%magmaGvd)
  #:if WITH_MAGMA
      call magmaHegvd(env%gpu%nGpu, HSqrReal, SSqrReal, eigen, 'L', jobz)
  #:else
      @:RAISE_ERROR(errStatus, -1, "This binary is compiled without GPU support")
  #:endif
    case default
      @:RAISE_ERROR(errStatus, -1, "Unknown eigensolver")
    end select

  end subroutine diagDenseRealMtx


  !> Diagonalizes a sparse represented Hamiltonian and overlap with k-points to give the
  !> eigenvectors and values, as well as often the Cholesky factorized overlap matrix (due to a side
  !> effect of lapack)
  subroutine diagDenseComplexMtx(env, electronicSolver, jobz, HSqrCplx, SSqrCplx, eigen, errStatus)

    !> Environment
    type(TEnvironment), intent(in) :: env

    !> Electronic solver information
    type(TElectronicSolver), intent(inout) :: electronicSolver

    !> type of eigen-problem, either 'V'/'v' vectors or 'N'/'n' eigenvalues only
    character, intent(in) :: jobz

    !> Large square matrix for the resulting eigenvectors
    complex(dp), intent(inout) :: HSqrCplx(:,:)

    !> Large square matrix for the overlap workspace, overwritten with the Cholesky factorized form.
    complex(dp), intent(inout) :: SSqrCplx(:,:)

    !> The eigenvalues of the matrices
    real(dp), intent(out) :: eigen(:)

    !> Status of routine
    type(TStatus), intent(out) :: errStatus

    @:ASSERT(size(HSqrCplx, dim=1) == size(HSqrCplx, dim=2))
    @:ASSERT(all(shape(HSqrCplx) == shape(SSqrCplx)))
    @:ASSERT(size(HSqrCplx, dim=1) == size(eigen))
    @:ASSERT(jobz == 'n' .or. jobz == 'N' .or. jobz == 'v' .or. jobz == 'V')

  #:block DEBUG_CODE
    call checkDiagonal(hSqrCplx, "hamiltonian", errStatus)
    @:PROPAGATE_ERROR(errStatus)
    call checkDiagonal(sSqrCplx, "overlap", errStatus)
    @:PROPAGATE_ERROR(errStatus)
  #:endblock DEBUG_CODE

    select case(electronicSolver%iSolver)
    case(electronicSolverTypes%QR)
      call hegv(HSqrCplx,SSqrCplx,eigen,'L',jobz)
    case(electronicSolverTypes%divideandconquer)
      call hegvd(HSqrCplx,SSqrCplx,eigen,'L',jobz)
    case(electronicSolverTypes%relativelyrobust)
      call hegvr(HSqrCplx,SSqrCplx,eigen,'L',jobz)
    case(electronicSolverTypes%magmaGvd)
  #:if WITH_MAGMA
      call magmaHegvd(env%gpu%nGpu, HSqrCplx, SSqrCplx, eigen, 'L', jobz)
  #:else
      @:RAISE_ERROR(errStatus, -1, "This binary is compiled without GPU support")
  #:endif
    case default
      @:RAISE_ERROR(errStatus, -1, "Unknown eigensolver")
    end select

  end subroutine diagDenseComplexMtx


#:block DEBUG_CODE
  !> Check diagonal elements of a matrix are real
  subroutine checkDiagonal(matrix, label, errStatus)

    !> Matrix to test
    complex(dp), intent(in) :: matrix(:,:)

    !> label string for error
    character(*), intent(in) :: label

    !> Status of routine
    type(TStatus), intent(out) :: errStatus

    integer :: iOrb

    @:ASSERT(size(matrix, dim=1) == size(matrix, dim=2))

    do iOrb = 1, size(matrix, dim=2)
      if (abs(aimag(matrix(iOrb, iOrb))) > 1024.0_dp*epsilon(0.0_dp)) then
        @:RAISE_FORMATTED_ERROR(errStatus, -1,&
            & "('Diagonal element ', I0, ' of ',A,' is complex:',2E20.12)", iOrb, trim(label),&
            & matrix(iOrb, iOrb))
      end if
    end do

  end subroutine checkDiagonal
#:endblock DEBUG_CODE


#:if WITH_SCALAPACK

  !> Diagonalizes a sparse represented Hamiltonian and overlap to give the eigenValsvectors and
  !> values, as well as often the Cholesky factorized overlap matrix (due to a side effect of
  !> lapack)
  subroutine diagDenseRealMtxBlacs(electronicSolver, iCholesky, jobz, desc, HSqr, SSqr, eigenVals,&
      & eigenVecs, errStatus)

    !> Electronic solver information
    type(TElectronicSolver), intent(inout) :: electronicSolver

    !> Index of the overlap matrix where Cholesky should be stored.
    integer, intent(in) :: iCholesky

    !> type of eigenVals-problem, either 'V'/'v' with vectors or 'N'/'n' eigenValsvalues only
    character, intent(in) :: jobz

    !> Dense descriptor
    integer, intent(in) :: desc(DLEN_)

    !> Large square matrix with hamiltonian
    real(dp), intent(inout) :: HSqr(:,:)

    !> Large square matrix for the overlap workspace, often overwritten with the Cholesky factorized
    !> form.
    real(dp), intent(inout) :: SSqr(:,:)

    !> Eigenvalues.
    real(dp), intent(out) :: eigenVals(:)

    !> Eigenvectors
    real(dp), intent(out) :: eigenVecs(:,:)

    !> Status of routine
    type(TStatus), intent(out) :: errStatus

    @:ASSERT(jobz == 'n' .or. jobz == 'N' .or. jobz == 'v' .or. jobz == 'V')

    if (electronicSolver%hasCholesky(iCholesky)) then
      call electronicSolver%getCholesky(iCholesky, SSqr)
    end if

    select case(electronicSolver%iSolver)
    case(electronicSolverTypes%QR)
      call scalafx_psygv(HSqr, desc, SSqr, desc, eigenVals, eigenVecs, desc, uplo="L", jobz=jobz,&
          & skipchol=electronicSolver%hasCholesky(iCholesky))
    case(electronicSolverTypes%divideandconquer)
      call scalafx_psygvd(HSqr, desc, SSqr, desc, eigenVals, eigenVecs, desc, uplo="L", jobz="V",&
          & allocfix=.true., skipchol=electronicSolver%hasCholesky(iCholesky))
    case(electronicSolverTypes%relativelyrobust)
      call scalafx_psygvr(HSqr, desc, SSqr, desc, eigenVals, eigenVecs, desc, uplo="L", jobz="V",&
          & skipchol=electronicSolver%hasCholesky(iCholesky))
    case(electronicSolverTypes%elpa)
      if (electronicSolver%elsi%tWriteHS) then
        call elsi_write_mat_real(electronicSolver%elsi%rwHandle, "ELSI_Hreal.bin", HSqr)
        call elsi_write_mat_real(electronicSolver%elsi%rwHandle, "ELSI_Sreal.bin", SSqr)
        call elsi_finalize_rw(electronicSolver%elsi%rwHandle)
        call cleanShutdown("Finished dense matrix write")
      end if
      ! ELPA solver, returns eigenstates
      ! note, this only factorises overlap on first call - no skipchol equivalent
      call elsi_ev_real(electronicSolver%elsi%handle, HSqr, SSqr, eigenVals, eigenVecs)

    case default
      @:RAISE_ERROR(errStatus, -1, "Unknown eigensolver")
    end select

    if (.not. electronicSolver%hasCholesky(1)) then
      call electronicSolver%storeCholesky(1, SSqr)
    end if

  end subroutine diagDenseRealMtxBlacs


  !> Diagonalizes a sparse represented Hamiltonian and overlap to give the eigenValsvectors and
  !> values, as well as often the Cholesky factorized overlap matrix (due to a side effect of
  !> lapack)
  subroutine diagDenseCplxMtxBlacs(env, electronicSolver, iCholesky, jobz, desc, HSqr, SSqr,&
      & eigenVals, eigenVecs, errStatus)

    !> Environment
    type(TEnvironment), intent(in) :: env

    !> Electronic solver information
    type(TElectronicSolver), intent(inout) :: electronicSolver

    !> Index of the overlap matrix where Cholesky should be stored.
    integer, intent(in) :: iCholesky

    !> type of eigenVals-problem, either 'V'/'v' with vectors or 'N'/'n' eigenValsvalues only
    character, intent(in) :: jobz

    !> Dense descriptor
    integer, intent(in) :: desc(DLEN_)

    !> Large square matrix for the resulting eigenValsvectors
    complex(dp), intent(inout) :: HSqr(:,:)

    !> Large square matrix for the overlap workspace, often overwritten with the Cholesky factorized
    !> form.
    complex(dp), intent(inout) :: SSqr(:,:)

    !> Eigenvalues.
    real(dp), intent(out) :: eigenVals(:)

    !> Eigenvectors
    complex(dp), intent(out) :: eigenVecs(:,:)

    !> Status of routine
    type(TStatus), intent(out) :: errStatus

    @:ASSERT(jobz == 'n' .or. jobz == 'N' .or. jobz == 'v' .or. jobz == 'V')

  #:block DEBUG_CODE
    call checkDiagonalBlacs(env, hSqr, "hamiltonian", desc, errStatus)
    @:PROPAGATE_ERROR(errStatus)
    call checkDiagonalBlacs(env, sSqr, "overlap", desc, errStatus)
    @:PROPAGATE_ERROR(errStatus)
  #:endblock DEBUG_CODE

    if (electronicSolver%hasCholesky(iCholesky)) then
      call electronicSolver%getCholesky(iCholesky, SSqr)
    end if

    select case(electronicSolver%iSolver)

    case(electronicSolverTypes%QR)
      call scalafx_phegv(HSqr, desc, SSqr, desc, eigenVals, eigenVecs, desc, uplo="L", jobz=jobz, &
          & skipchol=electronicSolver%hasCholesky(iCholesky))

    case(electronicSolverTypes%divideandconquer)
      call scalafx_phegvd(HSqr, desc, SSqr, desc, eigenVals, eigenVecs, desc, uplo="L", jobz=jobz,&
          & allocfix=.true., skipchol=electronicSolver%hasCholesky(iCholesky))

    case(electronicSolverTypes%relativelyrobust)
      call scalafx_phegvr(HSqr, desc, SSqr, desc, eigenVals, eigenVecs, desc, uplo="L", jobz=jobz,&
          & skipchol=electronicSolver%hasCholesky(iCholesky))

    case(electronicSolverTypes%elpa)
      if (electronicSolver%elsi%tWriteHS) then
        call elsi_write_mat_complex(electronicSolver%elsi%rwHandle, "ELSI_Hcmplx.bin", HSqr)
        call elsi_write_mat_complex(electronicSolver%elsi%rwHandle, "ELSI_Scmplx.bin", SSqr)
        call elsi_finalize_rw(electronicSolver%elsi%rwHandle)
        call cleanShutdown("Finished dense matrix write")
      end if
      ! ELPA solver, returns eigenstates
      ! note, this only factorises overlap on first call - no skipchol equivalent
      call elsi_ev_complex(electronicSolver%elsi%handle, HSqr, SSqr, eigenVals, eigenVecs)

    case default
      @:RAISE_ERROR(errStatus, -1, "Unknown eigensolver")
    end select

    if (.not. electronicSolver%hasCholesky(iCholesky)) then
      call electronicSolver%storeCholesky(iCholesky, SSqr)
    end if

  end subroutine diagDenseCplxMtxBlacs


#:block DEBUG_CODE
  !> Check for imaginary part on diagonals of (global) BLACS matrix
  subroutine checkDiagonalBlacs(env, matrix, label, desc, errStatus)

    !> Environment
    type(TEnvironment), intent(in) :: env

    !> Matrix to test the (global) diagonal elements
    complex(dp), intent(in) :: matrix(:,:)

    !> label string for error
    character(*), intent(in) :: label

    !> Dense descriptor
    integer, intent(in) :: desc(DLEN_)

    !> Status of routine
    type(TStatus), intent(out) :: errStatus

    integer :: iOrb, jOrb, iGlob, jGlob

    do jOrb = 1, size(matrix, dim=2)
      jGlob = scalafx_indxl2g(jOrb, desc(NB_), env%blacs%orbitalGrid%mycol, desc(CSRC_),&
          & env%blacs%orbitalGrid%ncol)
      do iOrb = 1, size(matrix, dim=1)
        iGlob = scalafx_indxl2g(iOrb, desc(MB_), env%blacs%orbitalGrid%myrow, desc(RSRC_),&
            & env%blacs%orbitalGrid%nrow)
        if (iGlob == jGlob) then
          if (abs(aimag(matrix(iOrb, jOrb))) > 1024.0_dp*epsilon(0.0_dp)) then
            @:RAISE_FORMATTED_ERROR(errStatus, -1,&
                & "('Diagonal element ', I0, ' of ',A,' is complex:',2E20.12)", iGlob, trim(label),&
                & matrix(iOrb, jOrb))
          end if
        end if
      end do
    end do

  end subroutine checkDiagonalBlacs
#:endblock DEBUG_CODE

#:endif

end module dftbp_dftbplus_eigenvects
