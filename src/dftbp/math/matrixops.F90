!--------------------------------------------------------------------------------------------------!
!  DFTB+: general package for performing fast atomistic simulations                                !
!  Copyright (C) 2006 - 2025  DFTB+ developers group                                               !
!                                                                                                  !
!  See the LICENSE file for terms of usage and distribution.                                       !
!--------------------------------------------------------------------------------------------------!

#:include 'common.fypp'
#:include 'error.fypp'

#:set FLAVOURS = [('Real', 'real', 'real'), ('Cmplx', 'complex', 'cmplx')]

!> Matrix operations for which LAPACK/ScaLAPACK does not have a direct call
module dftbp_math_matrixops
  use dftbp_common_accuracy, only : dp
  use dftbp_common_environment, only : TEnvironment
  use dftbp_common_schedule, only : assembleChunks
  use dftbp_common_status, only : TStatus
  use dftbp_io_message, only : error
  use dftbp_math_blasroutines, only : gemm
  use dftbp_math_eigensolver, only : heev
  use dftbp_math_lapackroutines, only : gesvd, getrf, getri, hetrf, hetri, sytrf, sytri
#:if WITH_SCALAPACK
  use dftbp_extlibs_mpifx, only : MPI_SUM, mpifx_allreduceip, mpifx_comm
  use dftbp_extlibs_scalapackfx, only : blacsgrid, CSRC_, DLEN_, M_, MB_, N_, NB_, pblasfx_ptran,&
      & pblasfx_ptranc, RSRC_, scalafx_indxl2g, scalafx_islocal, scalafx_pgetrf
#:endif

  implicit none

  private
  public :: adjointLowerTriangle, orthonormalizeVectors, adjugate, makeSimilarityTrans
  public :: matinv, symmatinv, hermatinv, det, calcMatrixSqrt, pseudoInv
#:if WITH_SCALAPACK
  public :: adjointLowerTriangle_BLACS

  !> Copy lower triangle of distributed matrix into the upper triangle, obeying hermitian symmetry
  !! if appropriate
  interface adjointLowerTriangle_BLACS
    module procedure symmetrize_BLACS
    module procedure hermitian_BLACS
  end interface adjointLowerTriangle_BLACS

#:endif

  !> Determinant of a matrix
  interface det
#:for SUFFIX, _, _ in FLAVOURS
    module procedure det${SUFFIX}$
#:if WITH_SCALAPACK
    module procedure detScaLAPACK${SUFFIX}$
#:endif
#:endfor
  end interface det

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


  !> perform a similarity (or unitary) transformation of a matrix
  interface makeSimilarityTrans
    module procedure makeSimilarityTrans_real
    module procedure makeSimilarityTrans_cmplx
  end interface makeSimilarityTrans


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

    real(dp), allocatable :: U(:,:), Sigma(:), Vt(:,:), Xi(:)
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


  !> unitary (or projective) transformation of a matrix X' = U X U^T or X' = U^T X U
  subroutine makeSimilarityTrans_real(xx, uu, side)

    !> matrix in original basis, U X U^T* on return.
    real(dp), intent(inout) :: xx(:,:)

    !> unitary matrix
    real(dp), intent(in) :: uu(:,:)

    !> which transform order to use, i.e. to which side the original unitary is applied
    character(1), intent(in), optional :: side

    real(dp) :: work(size(xx,dim=1), size(xx,dim=2))
    character(1) :: iSide

    if (present(side)) then
      iSide(:) = side
    else
      iSide(:) = 'L'
    end if

    @:ASSERT(all(shape(xx) == shape(uu)))
    @:ASSERT(size(xx, dim=1) == size(xx, dim=2))

    ! should blasify:
    select case(iSide)
    case ('L', 'l')
      work(:,:) = matmul(xx, transpose(uu))
      xx(:,:) = matmul(uu, work)
    case ('R', 'r')
      work(:,:) = matmul(xx, uu)
      xx(:,:) = matmul(transpose(uu), work)
    case default
      call error("Unknown unitary transform request")
    end select

  end subroutine makeSimilarityTrans_real


  !> unitary (or projective) transformation of a matrix X' = U X U^T* or X' = U^T* X U
  subroutine makeSimilarityTrans_cmplx(xx, uu, side)

    !> matrix in original basis, U X U^T* on return.
    complex(dp), intent(inout) :: xx(:,:)

    !> unitary matrix
    complex(dp), intent(in) :: uu(:,:)

    !> which transform order to use, i.e. which to side the original unitary is applied
    character(1), intent(in), optional :: side

    complex(dp) :: work(size(xx,dim=1), size(xx,dim=2))
    character(1) :: iSide

    if (present(side)) then
      iSide(:) = side
    else
      iSide(:) = 'L'
    end if

    @:ASSERT(all(shape(xx) == shape(uu)))
    @:ASSERT(size(xx, dim=1) == size(xx, dim=2))

    ! should blasify:
    select case(iSide)
    case ('L', 'l')
      work(:,:) = matmul(xx, transpose(conjg(uu)))
      xx(:,:) = matmul(uu, work)
    case ('R', 'r')
      work(:,:) = matmul(xx, uu)
      xx(:,:) = matmul(transpose(conjg(uu)), work)
    case default
      call error("Unknown unitary transform request")
    end select

  end subroutine makeSimilarityTrans_cmplx


  !> Inverts a general matrix.
  subroutine matinv(aa, nRow, iError)

    !> Matrix to invert on entry, inverted matrix on exit
    real(dp), intent(inout) :: aa(:,:)

    !> Nr. of rows of the matrix (if different from size(aa, dim=1)
    integer, intent(in), optional :: nRow

    !> iError Error flag. Returns 0 on successful operation. If this variable is not specified, any
    !> occurring error (e.g. singular matrix) stops the program.
    integer, intent(out), optional :: iError

    integer :: nn, info
    integer, allocatable :: ipiv(:)
    character(len=100) :: error_string

    nn = size(aa, dim=1)
    if (present(nRow)) then
      @:ASSERT(nRow >= 1 .and. nRow <= nn)
      nn = nRow
    end if
    @:ASSERT(size(aa, dim=2) >= nn)

    allocate(ipiv(nn))
    call getrf(aa, ipiv, nRow=nn, nColumn=nn, iError=info)
    if (info == 0) then
      call getri(aa, ipiv, nRow=nn, iError=info)
    end if

    if (present(iError)) then
      iError = info
    elseif (info /= 0) then
99120 format ('Matrix inversion failed because of error in getrf or getri.', &
          & ' Info flag: ',i10)
      write (error_string, 99120) info
      call error(error_string)
    end if

  end subroutine matinv


#:for SUFFIX, TYPE, KIND, NAME in [('sym', 'sy', 'real', 'symmetric'),&
  & ('her', 'he', 'complex', 'hermitian')]

  !> Inverts a ${NAME}$ matrix.
  subroutine ${SUFFIX}$matinv(aa, status, uplo)

    !> Symmetric matrix to invert on entry, inverted matrix on exit.
    ${KIND}$(dp), intent(inout) :: aa(:,:)

    !> Status of operation
    type(TStatus), intent(out) :: status

    !> Upper ('U') or lower ('L') matrix. Default: 'L'.
    character, intent(in), optional :: uplo

    integer :: nn
    integer, allocatable :: ipiv(:)

    nn = size(aa, dim=1)
    allocate(ipiv(nn))

    call ${TYPE}$trf(aa, ipiv, status, uplo=uplo)
    @:PROPAGATE_ERROR(status)

    call ${TYPE}$tri(aa, ipiv, status, uplo=uplo)
    @:PROPAGATE_ERROR(status)

  end subroutine ${SUFFIX}$matinv

#:endfor


#:for SUFFIX, TYPE, CONVERT in FLAVOURS

  !> Determinant of a matrix, matrix destroyed in process
  function det${SUFFIX}$(A) result(det)

    !> The matrix
    ${TYPE}$(dp), intent(inout) :: A(:,:)

    !> resulting determinant
    ${TYPE}$(dp) :: det

    integer, allocatable  :: ipiv(:)
    integer :: ii, n, exponent

    n = minval(shape(A))
    allocate(ipiv(n))

    call getrf(A,ipiv)

    det = ${CONVERT}$(1, kind=dp)
    exponent = 0
    do ii = 1, n
      if (ipiv(ii) /= ii) then
        det = -det * A(ii,ii)
      else
        det = det * A(ii,ii)
      end if
      if (det == 0.0_dp) then
        return
      end if
      do while (abs(det) > 2.0_dp)
        det = det / 2.0_dp
        exponent = exponent + 1
      end do
      do while (abs(det) < 0.5_dp)
        det = det * 2.0_dp
        exponent = exponent - 1
      end do
    end do
    det = det * 2.0_dp ** exponent

  end function det${SUFFIX}$

#:if WITH_SCALAPACK

  !> Determinant of a matrix, matrix destroyed in process
  function detScaLAPACK${SUFFIX}$(A, descA, grid, myComm) result(det)

    !> The matrix
    ${TYPE}$(dp), intent(inout) :: A(:,:)

    !> Dense descriptor
    integer, intent(in) :: descA(DLEN_)

    !> BLACS grid involved in calculation
    type(blacsgrid), intent(in) :: grid

    !> Communicator for the region involved in the BLACS grid
    type(mpifx_comm), intent(in) :: myComm

    !> resulting determinant
    ${TYPE}$(dp) :: det

    integer, allocatable  :: ipiv(:)
    integer :: ii, jj, iLoc, jLoc, mm, nn
    logical :: tDiagBlock, tAnyDiag
    ${TYPE}$(dp) :: detLocal
    ${TYPE}$(dp), allocatable :: detBuffer(:)
    integer :: expLocal
    integer, allocatable :: expBuffer(:)

    @:ASSERT(grid%nProc == myComm%size)

    if (grid%iproc /= -1) then
      mm = descA(M_)
      nn = descA(N_)

      allocate(detBuffer(grid%nProc))
      allocate(expBuffer(grid%nProc))
      detBuffer = 0.0_dp
      expbuffer = 0.0_dp

      allocate(ipiv(min(mm,nn)))
      ipiv = 0
      call scalafx_pgetrf(A,descA,ipiv)

      ! note, this includes under-/over-flow protection similar to LINPACK routine dgedi.f
      detLocal = 1.0_dp
      expLocal = 0
      tAnyDiag = .false.
      lpLocal: do ii = 1, size(A,dim=2)

        ! Look for diagonal blocks
        jj = scalafx_indxl2g(ii, descA(NB_), grid%mycol, descA(CSRC_), grid%ncol)
        call scalafx_islocal(grid, descA, jj, jj, tDiagBlock, iLoc, jLoc)

        tAnyDiag = tAnyDiag .or. tDiagBlock

        if (tDiagBlock) then
          if (jj /= ipiv(ii)) then
            detLocal = -detLocal * A(iLoc,jLoc)
          else
            detLocal = detLocal * A(iLoc,jLoc)
          end if

          if (detLocal == 0.0_dp) then
            exit lpLocal
          end if

          do while (abs(detLocal) > 2)
            detLocal = detLocal / 2.0_dp
            expLocal = expLocal + 1
          end do
          do while (abs(detLocal) < 0.5_dp)
            detLocal = detLocal * 2.0_dp
            expLocal = expLocal - 1
          end do
        end if

      end do lpLocal

      if (tAnyDiag) then
        detBuffer(grid%iProc+1) = detLocal
        expBuffer(grid%iProc+1) = expLocal
      else
        ! node did not have any diagonal elements, so does not contribute to det
        detBuffer(grid%iProc+1) = 1.0_dp
        expBuffer(grid%iProc+1) = 0.0_dp
      end if

      ! now product the full det from the sub-processes
      call mpifx_allreduceip(myComm, detBuffer, MPI_SUM)
      call mpifx_allreduceip(myComm, expBuffer, MPI_SUM)

      detLocal = ${CONVERT}$(1, kind=dp)
      expLocal = 0
      lpTotal: do ii = 1, grid%nProc
        detLocal = detLocal * detBuffer(ii)
        expLocal = expLocal + expBuffer(ii)
        if (detLocal == 0.0_dp) then
          exit lpTotal
        end if
        do while (abs(detLocal) > 2)
          detLocal = detLocal / 2.0_dp
          expLocal = expLocal + 1
        end do
        do while (abs(detLocal) < 0.5_dp)
          detLocal = detLocal * 2.0_dp
          expLocal = expLocal - 1
        end do
      end do lpTotal

      det = detLocal * 2.0_dp ** expLocal

    end if

  end function detScaLAPACK${SUFFIX}$

#:endif

#:endfor


  !> Calculate square root and inverse of sqrt of a real, symmetric positive definite matrix.
  !! Note, faster and more stable algorithm would be to Cholesky factor the matrix as A = R* R, then
  !! Use the polar decomposition of R^-1 H U, where U is unitary and H is hermitian positive
  !! definite, as performed with the polar-Newton method in Chapter 6 of Functions of Matrices by
  !! N. Higham. Then A^.5 = UR and A^-.5 = R^-1 U*
  subroutine calcMatrixSqrt(matIn, matSqrt, matSqrtInv)

    !> Matrix to operate on
    real(dp), intent(in) :: matIn(:,:)

    !> Matrix square root
    real(dp), intent(out) :: matSqrt(:,:)

    !> Inverse of matrix square root
    real(dp), intent(out) :: matSqrtInv(:,:)

    real(dp), allocatable :: dummyEV(:), dummyM(:, :), dummyM2(:, :)
    integer :: ii, spaceDim

    spaceDim = size(matIn, dim=2)
    allocate(dummyEV(spaceDim))
    allocate(dummyM(spaceDim, spaceDim))
    allocate(dummyM2(spaceDim, spaceDim))

    dummyM(:,:) = matIn

    call heev(dummyM, dummyEV, 'U', 'V')

    ! Calculate matrix sqrt
    do ii = 1, spaceDim
      dummyM2(:,ii) = sqrt(dummyEV(ii)) * dummyM(:,ii)
    end do

    call gemm(matSqrt, dummyM2, dummyM, transB='T')

    ! Calculate invverse of matrix sqrt
    do ii = 1, spaceDim
      dummyM2(:,ii) = dummyM(:,ii) / sqrt(dummyEV(ii))
    end do

    call gemm(matSqrtInv, dummyM2, dummyM, transB='T')

  end subroutine calcMatrixSqrt


  !> Moore-Penrose pseudo-inverse of general rectangular matrix
  subroutine pseudoInv(A, Ainv)

    !> Matrix A(m,n), overwritten on output
    real(dp), intent(inout) :: A(:,:)

    !> Pseudoinverse, returned as transpose Ainv(m,n)
    real(dp), intent(out) :: Ainv(:,:)

    real(dp), allocatable :: U(:,:), sigma(:), Vt(:,:)
    integer :: m, n, mn, ii

    m = size(A, dim=1)
    n = size(A, dim=2)
    mn = min(m, n)

    @:ASSERT(all(shape(Ainv) == [m,n]))

    allocate(U(m,mn))
    allocate(Vt(mn,n))
    allocate(sigma(mn))

    call gesvd(A, U, sigma, Vt)

    where(sigma > epsilon(0.0_dp))
      sigma = 1.0_dp / sigma
    elsewhere
      sigma = 0.0_dp
    end where

    do ii = 1, mn
      U(:, ii) = U(:, ii) * sigma(ii)
    end do

    Ainv(:,:) = matmul(U,Vt)

  end subroutine pseudoInv

end module dftbp_math_matrixops
