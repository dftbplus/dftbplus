!--------------------------------------------------------------------------------------------------!
!  DFTB+: general package for performing fast atomistic simulations                                !
!  Copyright (C) 2018  DFTB+ developers group                                                      !
!                                                                                                  !
!  See the LICENSE file for terms of usage and distribution.                                       !
!--------------------------------------------------------------------------------------------------!

#:include 'common.fypp'

!> Module to wrap around the process of converting from a Hamiltonian and overlap in sparse form
!> into eigenvectors
module eigenvects
  use assert
  use accuracy
  use eigensolver
  use message
  use solvertypes
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
  subroutine diagDenseRealMtx(iSolver, jobz, HSqrReal, SSqrReal, eigen)

    !> Choice of eigensolver, 4 different lapack dense solvers currently supported
    integer, intent(in) :: iSolver

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

    select case(iSolver)
    case(solverQR)
      call hegv(HSqrReal,SSqrReal,eigen,'L',jobz)
    case(solverDAC)
      call hegvd(HSqrReal,SSqrReal,eigen,'L',jobz)
    case(solverRR)
      call gvr(HSqrReal,SSqrReal,eigen,'L',jobz)
    case default
      call error('Unknown eigensolver')
    end select

  end subroutine diagDenseRealMtx


  !> Diagonalizes a sparse represented Hamiltonian and overlap with k-points to give the
  !> eigenvectors and values, as well as often the Cholesky factorized overlap matrix (due to a side
  !> effect of lapack)
  subroutine diagDenseComplexMtx(iSolver, jobz, HSqrCplx, SSqrCplx, eigen)

    !> Choice of eigensolver, 4 different lapack dense solvers currently supported
    integer, intent(in) :: iSolver

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

    select case(iSolver)
    case(solverQR)
      call hegv(HSqrCplx,SSqrCplx,eigen,'L',jobz)
    case(solverDAC)
      call hegvd(HSqrCplx,SSqrCplx,eigen,'L',jobz)
    case(solverRR)
      call gvr(HSqrCplx,SSqrCplx,eigen,'L',jobz)
    case default
      call error('Unknown eigensolver')
    end select

  end subroutine diagDenseComplexMtx


#:if WITH_SCALAPACK

  !> Diagonalizes a sparse represented Hamiltonian and overlap to give the eigenValsvectors and
  !> values, as well as often the Cholesky factorized overlap matrix (due to a side effect of
  !> lapack)
  subroutine diagDenseRealMtxBlacs(iSolver, jobz, desc, HSqr, SSqr, eigenVals, eigenVecs)

    !> Choice of eigensolver, 4 different lapack dense solvers currently supported
    integer, intent(in) :: iSolver

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

    select case(iSolver)
    case(solverQR)
      call scalafx_psygv(HSqr, desc, SSqr, desc, eigenVals, eigenVecs, desc, uplo="L", jobz=jobz)
    case(solverDAC)
      call scalafx_psygvd(HSqr, desc, SSqr, desc, eigenVals, eigenVecs, desc, uplo="L", jobz="V",&
          & allocfix=.true.)
    case(solverRR)
      call scalafx_psygvr(HSqr, desc, SSqr, desc, eigenVals, eigenVecs, desc, uplo="L", jobz="V")
    case default
      call error('Unknown eigensolver')
    end select

  end subroutine diagDenseRealMtxBlacs



  !> Diagonalizes a sparse represented Hamiltonian and overlap to give the eigenValsvectors and
  !> values, as well as often the Cholesky factorized overlap matrix (due to a side effect of
  !> lapack)
  subroutine diagDenseCplxMtxBlacs(iSolver, jobz, desc, HSqr, SSqr, eigenVals, eigenVecs)

    !> Choice of eigensolver, 4 different lapack dense solvers currently supported
    integer, intent(in) :: iSolver

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

    select case(iSolver)
    case(solverQR)
      call scalafx_phegv(HSqr, desc, SSqr, desc, eigenVals, eigenVecs, desc, uplo="L", jobz=jobz)
    case(solverDAC)
      call scalafx_phegvd(HSqr, desc, SSqr, desc, eigenVals, eigenVecs, desc, uplo="L", jobz=jobz,&
          & allocfix=.true.)
    case(solverRR)
      call scalafx_phegvr(HSqr, desc, SSqr, desc, eigenVals, eigenVecs, desc, uplo="L", jobz=jobz)
    case default
      call error('Unknown eigensolver')
    end select

  end subroutine diagDenseCplxMtxBlacs

#:endif

end module eigenvects
