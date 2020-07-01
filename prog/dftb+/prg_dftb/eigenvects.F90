!--------------------------------------------------------------------------------------------------!
!  DFTB+: general package for performing fast atomistic simulations                                !
!  Copyright (C) 2018  DFTB+ developers group                                                      !
!                                                                                                  !
!  See the LICENSE file for terms of usage and distribution.                                       !
!--------------------------------------------------------------------------------------------------!

#:include 'common.fypp'

!> Contains routines to obtain eigenvalues and eigenvectors.
module eigenvects
  use assert
  use accuracy
  use eigensolver
  use dftbp_solvers
  use message
#:if WITH_SCALAPACK
  use scalapackfx
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
  subroutine diagDenseRealMtx(solver, jobz, HSqrReal, SSqrReal, eigen)

    !> Choice of eigensolver, currently only lapack dense solvers supported
    type(TSolver), intent(in) :: solver

    !> type of eigen-problem, either 'V'/'v' with vectors or 'N'/'n' eigenvalues only
    character, intent(in) :: jobz

    !> Large square matrix for the resulting eigenvectors
    real(dp), intent(inout) :: HSqrReal(:,:)

    !> Large square matrix for the overlap workspace, often overwritten with the Cholesky factorized
    !> form.
    real(dp), intent(inout) :: SSqrReal(:,:)

    !> Eigen values.
    real(dp), intent(out) :: eigen(:)


    @:ASSERT(size(HSqrReal, dim=1) == size(HSqrReal, dim=2))
    @:ASSERT(all(shape(HSqrReal) == shape(SSqrReal)))
    @:ASSERT(size(HSqrReal, dim=1) == size(eigen))
    @:ASSERT(jobz == 'n' .or. jobz == 'N' .or. jobz == 'v' .or. jobz == 'V')

    select case (solver%solverType)
    case (solverTypes%lapackQr)
      call hegv(HSqrReal,SSqrReal,eigen,'L',jobz)
    case (solverTypes%lapackDivAndConq)
      call hegvd(HSqrReal,SSqrReal,eigen,'L',jobz)
    case (solverTypes%lapackRelRobust)
      call gvr(HSqrReal,SSqrReal,eigen,'L',jobz)
    case default
      call error('Internal error: Unknown eigensolver in diagDenseRealMtx')
    end select

  end subroutine diagDenseRealMtx


  !> Diagonalizes a sparse represented Hamiltonian and overlap with k-points to give the
  !> eigenvectors and values, as well as often the Cholesky factorized overlap matrix (due to a side
  !> effect of lapack)
  subroutine diagDenseComplexMtx(solver, jobz, HSqrCplx, SSqrCplx, eigen)

    !> Choice of eigensolver, currently only lapack dense solvers supported
    type(TSolver), intent(in) :: solver

    !> type of eigen-problem, either 'V'/'v' vectors or 'N'/'n' eigenvalues only
    character, intent(in) :: jobz

    !> Large square matrix for the resulting eigenvectors
    complex(dp), intent(inout) :: HSqrCplx(:,:)

    !> Large square matrix for the overlap workspace, overwritten with the Cholesky factorized form.
    complex(dp), intent(inout) :: SSqrCplx(:,:)

    !> The eigenvalues of the matrices
    real(dp), intent(out) :: eigen(:)


    @:ASSERT(size(HSqrCplx, dim=1) == size(HSqrCplx, dim=2))
    @:ASSERT(all(shape(HSqrCplx) == shape(SSqrCplx)))
    @:ASSERT(size(HSqrCplx, dim=1) == size(eigen))
    @:ASSERT(jobz == 'n' .or. jobz == 'N' .or. jobz == 'v' .or. jobz == 'V')

    select case(solver%solverType)
    case (solverTypes%lapackQr)
      call hegv(HSqrCplx,SSqrCplx,eigen,'L',jobz)
    case (solverTypes%lapackDivAndConq)
      call hegvd(HSqrCplx,SSqrCplx,eigen,'L',jobz)
    case (solverTypes%lapackRelRobust)
      call gvr(HSqrCplx,SSqrCplx,eigen,'L',jobz)
    case default
      call error('Internal error: Unknown eigensolver in diagDenseComplexMtx')
    end select

  end subroutine diagDenseComplexMtx


#:if WITH_SCALAPACK

  !> Diagonalizes a sparse represented Hamiltonian and overlap to give the eigenValsvectors and
  !> values, as well as often the Cholesky factorized overlap matrix (due to a side effect of
  !> lapack)
  subroutine diagDenseRealMtxBlacs(solver, jobz, desc, HSqr, SSqr, eigenVals, eigenVecs)

    !> Choice of eigensolver, currently only lapack dense solvers supported
    type(TSolver), intent(in) :: solver

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


    @:ASSERT(jobz == 'n' .or. jobz == 'N' .or. jobz == 'v' .or. jobz == 'V')

    select case (solver%solverType)
    case (solverType%lapackQr)
      call scalafx_psygv(HSqr, desc, SSqr, desc, eigenVals, eigenVecs, desc, uplo="L", jobz=jobz)
    case (solverType%lapackDivAndConq)
      call scalafx_psygvd(HSqr, desc, SSqr, desc, eigenVals, eigenVecs, desc, uplo="L", jobz="V",&
          & allocfix=.true.)
    case (solverType%lapackRelRobust)
      call scalafx_psygvr(HSqr, desc, SSqr, desc, eigenVals, eigenVecs, desc, uplo="L", jobz="V")
    case default
      call error('Internal error: Unknown eigensolver in diagDenseRealMtxBlacs')
    end select

  end subroutine diagDenseRealMtxBlacs


  !> Diagonalizes a sparse represented Hamiltonian and overlap to give the eigenValsvectors and
  !> values, as well as often the Cholesky factorized overlap matrix (due to a side effect of
  !> lapack)
  subroutine diagDenseCplxMtxBlacs(solver, jobz, desc, HSqr, SSqr, eigenVals, eigenVecs)

    !> Choice of eigensolver, currently only lapack dense solvers supported
    type(TSolver), intent(in) :: solver

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


    @:ASSERT(jobz == 'n' .or. jobz == 'N' .or. jobz == 'v' .or. jobz == 'V')

    select case (solver%solverType)
    case (solverTypes%lapackQr)
      call scalafx_phegv(HSqr, desc, SSqr, desc, eigenVals, eigenVecs, desc, uplo="L", jobz=jobz)
    case (solverTypes%lapackDivAndConq)
      call scalafx_phegvd(HSqr, desc, SSqr, desc, eigenVals, eigenVecs, desc, uplo="L", jobz=jobz,&
          & allocfix=.true.)
    case (solverTypes%lapackRelRobust)
      call scalafx_phegvr(HSqr, desc, SSqr, desc, eigenVals, eigenVecs, desc, uplo="L", jobz=jobz)
    case default
      call error('Internal error: Unknown eigensolver in diagDenseCplxMtxBlacs')
    end select

  end subroutine diagDenseCplxMtxBlacs

#:endif

end module eigenvects
