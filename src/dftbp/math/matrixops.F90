!--------------------------------------------------------------------------------------------------!
!  DFTB+: general package for performing fast atomistic simulations                                !
!  Copyright (C) 2006 - 2023  DFTB+ developers group                                               !
!                                                                                                  !
!  See the LICENSE file for terms of usage and distribution.                                       !
!--------------------------------------------------------------------------------------------------!

#:include 'common.fypp'
#:include 'error.fypp'

!> Simple matrix operations for which LAPACK does not have a direct call
module dftbp_math_matrixops
  use dftbp_common_accuracy, only : dp
  use dftbp_common_environment, only : TEnvironment
  use dftbp_common_schedule, only : assembleChunks
  use dftbp_common_status, only : TStatus
  use dftbp_math_determinant, only : det
  use dftbp_math_lapackroutines, only : matinv, gesvd
#:if WITH_SCALAPACK
  use dftbp_extlibs_mpifx, only : MPI_SUM, mpifx_allreduceip
  use dftbp_extlibs_scalapackfx, only : DLEN_, CSRC_, RSRC_, MB_, NB_, pblasfx_ptranc,&
      & pblasfx_ptran, scalafx_indxl2g
#:endif

  implicit none

  private
  public :: adjointLowerTriangle, orthonormalizeVectors, adjugate
#:if WITH_SCALAPACK
  public :: adjointLowerTriangle_BLACS

  !> Copy lower triangle of distributed matrix into the upper triangle, obeying hermitian symmetry
  !! if appropriate
  interface adjointLowerTriangle_BLACS
    module procedure symmetrize_BLACS
    module procedure hermitian_BLACS
  end interface adjointLowerTriangle_BLACS

#:endif

  !> Copy lower triangle into the upper triangle of a square matrix, obeying hermitian symmetry if
  !! appropriate
  interface adjointLowerTriangle
    module procedure symmetrizeSquareMatrix
    module procedure hermitianSquareMatrix
  end interface adjointLowerTriangle


  !> Adjugate of a well behaved square matrix
  interface adjugate
    module procedure adjugate_stable
    module procedure adjugate_simple
  end interface adjugate


contains


  !> Copy lower triangle to upper for a square matrix.
  subroutine symmetrizeSquareMatrix(matrix)

    !> matrix to symmetrize
    real(dp), intent(inout) :: matrix(:,:)
    integer :: ii, matSize

    matSize = size(matrix, dim = 1)
    do ii = 1, matSize -1
      matrix(ii, ii+1:) = matrix(ii+1:, ii)
    end do

  end subroutine symmetrizeSquareMatrix


  !> Copy lower triangle to upper for a square matrix with Hermitian symmetry
  subroutine hermitianSquareMatrix(matrix)

    !> matrix to symmetrize
    complex(dp), intent(inout) :: matrix(:,:)
    integer :: ii, matSize

    matSize = size(matrix, dim = 1)
    do ii = 1, matSize -1
      matrix(ii, ii+1:) = conjg(matrix(ii+1:, ii))
    end do

  end subroutine hermitianSquareMatrix


#:if WITH_SCALAPACK

  !> Copy upper triangle into lower triangle of distributed matrix
  subroutine symmetrize_BLACS(desc, myCol, myRow, nCol, nRow, matrix)

    !> BLACS matrix descriptor
    integer, intent(in) :: desc(DLEN_)

    !> Column of the current process in the BLACS grid
    integer, intent(in) :: myCol

    !> Row of the current process in the BLACS grid
    integer, intent(in) :: myRow

    !> Number. of process columns in the grid
    integer, intent(in) :: nCol

    !> Number. of process rows in the grid
    integer, intent(in) :: nRow

    !> Matrix to symmetrize
    real(dp), intent(inout) :: matrix(:,:)

    real(dp), allocatable :: work(:,:)
    integer :: ii, jj, iGlob, jGlob

    allocate(work(size(matrix, dim=1), size(matrix, dim=2)), source = 0.0_dp)
    call pblasfx_ptran(matrix, desc, work, desc)
    do jj = 1, size(matrix,dim=2)
      jGlob = scalafx_indxl2g(jj, desc(NB_), mycol, desc(CSRC_), ncol)
      do ii = 1, size(matrix,dim=1)
        iGlob = scalafx_indxl2g(ii, desc(MB_), myrow, desc(RSRC_), nrow)
        if (iGlob < jGlob) then
          matrix(ii, jj) = work(ii, jj)
        end if
      end do
    end do

  end subroutine symmetrize_BLACS


  !> Copy upper triangle into lower triangle of distributed hermitian matrix
  subroutine hermitian_BLACS(desc, myCol, myRow, nCol, nRow, matrix)

    !> BLACS matrix descriptor
    integer, intent(in) :: desc(DLEN_)

    !> Column of the current process in the BLACS grid
    integer, intent(in) :: myCol

    !> Row of the current process in the BLACS grid
    integer, intent(in) :: myRow

    !> Number. of process columns in the grid
    integer, intent(in) :: nCol

    !> Number. of process rows in the grid
    integer, intent(in) :: nRow

    !> Matrix to hermitian symmetrize
    complex(dp), intent(inout) :: matrix(:,:)

    complex(dp), allocatable :: work(:,:)
    integer :: ii, jj, iGlob, jGlob

    allocate(work(size(matrix, dim=1), size(matrix, dim=2)), source = (0.0_dp, 0.0_dp))
    call pblasfx_ptranc(matrix, desc, work, desc)
    do jj = 1, size(matrix,dim=2)
      jGlob = scalafx_indxl2g(jj, desc(NB_), mycol, desc(CSRC_), ncol)
      do ii = 1, size(matrix,dim=1)
        iGlob = scalafx_indxl2g(ii, desc(MB_), myrow, desc(RSRC_), nrow)
        if (iGlob < jGlob) then
          matrix(ii, jj) = work(ii, jj)
        elseif (iGlob == jGlob) then
          matrix(ii, jj) = real(matrix(ii, jj),dp)
        end if
      end do
    end do

  end subroutine hermitian_BLACS

#:endif


  !> Perform modified Gram-Schmidt orthonormalization of the vectors in the columns of
  !! vecs(:,start:end), while also keeping them orthogonal to vecs(:,:start-1) (which are assumed to
  !! already be orthogonal)
  subroutine orthonormalizeVectors(env, iStart, iEnd, vecs)

    !> Environment settings
    type(TEnvironment), intent(in) :: env

    !> Starting place in vectors to work from
    integer, intent(in) :: iStart

    !> Ending place in vectors
    integer, intent(in) :: iEnd

    !> Vectors to be orthogonalized against 1:end vectors
    real(dp), intent(inout) :: vecs(:,:)

    integer :: ii, jj
    real(dp) :: dummyReal

    ! Obviously, not optimal in terms of communication, can be optimized if necessary. Assumes vecs
    ! are column block distributed (not block cyclic column and row) if parallel.
    do ii = iStart, iEnd
      do jj = 1, ii - 1
        dummyReal = dot_product(vecs(:,ii), vecs(:,jj))
        call assembleChunks(env, dummyReal)
        vecs(:,ii) = vecs(:,ii) - dummyReal * vecs(:,jj)
      end do
      dummyReal = dot_product(vecs(:,ii), vecs(:,ii))
      call assembleChunks(env, dummyReal)
      vecs(:,ii) = vecs(:,ii) / sqrt(dummyReal)
    end do

  end subroutine orthonormalizeVectors


  !> Evaluate adjugate of matrix in stable way
  subroutine adjugate_stable(A)

    !> Matrix for which to evaluate adjugate, over-written on exit
    real(dp), intent(inout) :: A(:,:)

    real(dp), allocatable :: aTmp(:,:), U(:,:), Sigma(:), Vt(:,:), Xi(:)
    real(dp) :: eta
    integer :: ii, jj, kk, n

    n = size(A, dim=2)
    allocate(U(n,n))
    allocate(Sigma(n))
    allocate(Vt(n,n))
    allocate(Xi(n))

    call gesvd(A, U, Sigma, Vt)

    do ii = 1, n
      Xi(ii) = product_real(Sigma, ii)
    end do

    A(:,:) = 0.0_dp
    do ii = 1, n
      do jj = 1, n
        do kk = 1, n
          A(kk,jj) = A(kk,jj) + Vt(ii,kk) * Xi(ii) * U(jj,ii)
        end do
      end do
    end do

    U(:,:) = matmul(U, Vt)
    eta = det(U)
    A(:,:) = eta * A

  end subroutine adjugate_stable


  !> Evaluate adjugate of matrix from det(A) A^-1, note can be unstable even when resulting quantity
  !! is well defined
  subroutine adjugate_simple(A, status)

    !> Matrix for which to evaluate adjugate, over-written on exit
    real(dp), intent(inout) :: A(:,:)

    !> Status of operation
    type(TStatus), intent(out) :: status

    real(dp), allocatable :: aTmp(:,:)
    integer :: iErr

    aTmp = A
    call matinv(A, iError=iErr)
    if (iErr /= 0) then
      @:RAISE_FORMATTED_ERROR(status, -1, "('Matrix inversion failure in adjugate, info: ',I0)",&
          & iErr)
    end if
    A(:,:) = det(aTmp) * A

  end subroutine adjugate_simple


  !> Numerically stable product of an array of numbers, excluding a specified value
  function product_real(Sigma, ii) result(res)

    !> Numbers to product
    real(dp), intent(in) :: Sigma(:)

    !> Index of value to omit
    integer, intent(in) :: ii

    real(dp), parameter :: base = 10.0_dp
    real(dp), parameter :: invBase = 1.0_dp / base
    real(dp) :: res
    integer :: jj, exponent

    res = 1.0_dp
    exponent = 0
    do jj = 1, size(sigma)
      if (ii == jj) cycle
      res = Sigma(jj) * res
      if (res < epsilon(0.0_dp)) then
        res = 0.0_dp
        return
      end if
      do while (abs(res) > base)
        res = res * invBase
        exponent = exponent + 1
      end do
      do while (abs(res) < invBase)
        res = res * base
        exponent = exponent - 1
      end do
    end do
    res = res * base ** exponent

  end function product_real

end module dftbp_math_matrixops
