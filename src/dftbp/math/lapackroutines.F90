!--------------------------------------------------------------------------------------------------!
!  DFTB+: general package for performing fast atomistic simulations                                !
!  Copyright (C) 2006 - 2023  DFTB+ developers group                                               !
!                                                                                                  !
!  See the LICENSE file for terms of usage and distribution.                                       !
!--------------------------------------------------------------------------------------------------!

#:include 'common.fypp'
#:include 'error.fypp'

#! suffix and kinds for real types
#:set REAL_KIND_PARAMS = [('real', 's'), ('dble', 'd')]

!> Contains F90 wrapper functions for some commonly used lapack calls needed in the code. The
!> interface of all LAPACK calls must be defined in the module lapack.
module dftbp_math_lapackroutines
  use dftbp_common_accuracy, only : dp, rdp, rsp
  use dftbp_common_status, only : TStatus
  use dftbp_io_message, only : error, warning
  ! use dftbp_extlibs_lapack
  implicit none

  private
  public :: gesv, getri, getrf, sytri, sytrf, matinv, symmatinv, larnv
  public :: hermatinv, hetri, hetrf, gesvd, potrf, trsm, getrs

  !> Computes the solution to a real system of linear equations A * X = B, where A is an N-by-N
  !> matrix and X and B are N-by-NRHS matrices
  interface gesv
    module procedure gesv_real
    module procedure gesv_dble
    module procedure gesv_dcomplex
  end interface gesv


  !> Computes the LU decomposition of a general rectangular matrix using partial pivoting with row
  !> interchanges.
  !> The decomposition has the form: A = P*L*U, where P is a permutation matrix, L is a lower
  !> triangular matrix with unit diagonal elements and U is an upper triangular matrix.
  interface getrf
    module procedure getrf_real
    module procedure getrf_dble
    module procedure getrf_complex
    module procedure getrf_dcomplex
  end interface getrf


  !> Bunch-Kaufman factorization of a symmetric matrix.
  interface sytrf
    module procedure sytrf_real
    module procedure sytrf_dreal
  end interface sytrf


  !> Bunch-Kaufman factorization of a Hermitian matrix.
  interface hetrf
    module procedure hetrf_cmplx
    module procedure hetrf_dcmplx
  end interface hetrf


  !> Inverts a symmetric matrix.
  interface sytri
    module procedure sytri_real
    module procedure sytri_dreal
  end interface sytri


  !> Inverts a Hermitian matrix.
  interface hetri
    module procedure hetri_cmplx
    module procedure hetri_dcmplx
  end interface hetri


  !> Computes the inverse of a matrix using LU factorization computed by getrf.
  interface getri
    module procedure getri_real
    module procedure getri_dble
  end interface getri


  !> returns a vector of random numers, either from a uniform or normal distribution
  interface larnv
    module procedure larnv_real
    module procedure larnv_dble
    module procedure larnv_cplx
    module procedure larnv_dblecplx
  end interface larnv

  !> svd decomposition of matrix A into left and right vectors and singular values U S V^dag
  interface gesvd
    module procedure sgesvd_real
    module procedure dgesvd_dble
    module procedure cgesvd_cplx
    module procedure zgesvd_dblecplx
  end interface gesvd

  interface potrf
  #:for suffix, kind in REAL_KIND_PARAMS
    module procedure ${kind}$potrf_${suffix}$
  #:endfor
  end interface potrf

  interface trsm
  #:for suffix, kind in REAL_KIND_PARAMS
    module procedure ${kind}$trsm_${suffix}$
  #:endfor
  end interface trsm


  !> Solves a system of linear equations
  !>  A * X = B  or  A**T * X = B
  !> with a general N-by-N matrix A using the LU factorization computed by getrf.
  interface getrs
    module procedure :: getrs_dble
    module procedure :: getrs1_dble
    module procedure :: getrs_real
    module procedure :: getrs1_real
  end interface getrs




contains


  !> Single precision version of gesv
  subroutine gesv_real(aa, bb, nEquation, nSolution, iError)

    !> Contains the coefficients on entry, the LU factorisation on exit.
    real(rsp), intent(inout) :: aa(:,:)

    !> Right hand side(s) of the linear equation on entry, solution(s) on exit.
    real(rsp), intent(inout) :: bb(:,:)

    !> The size of the problem (nr. of variables and equations). Must be only specified if different
    !> from size(aa, dim=1).
    integer, intent(in), optional :: nEquation

    !> Nr. of right hand sides (nr. of solutions). Must be only specified if different from size(b,
    !> dim=2).
    integer, intent(in), optional :: nSolution

    !> Error flag. If present, Lapack error flags are reported and noncritical errors (iError > 0)
    !> will not abort the program.
    integer, intent(out), optional :: iError
    integer :: info
    integer :: nn, nrhs, lda, ldb
    integer, allocatable :: ipiv(:)
    character(len=100) :: error_string


    lda = size(aa, dim=1)
    if (present(nEquation)) then
      @:ASSERT(nEquation >= 1 .and. nEquation <= lda)
      nn = nEquation
    else
      nn = lda
    end if
    @:ASSERT(size(aa, dim=2) >= nn)

    ldb = size(bb, dim=1)
    @:ASSERT(ldb >= nn)
    nrhs = size(bb, dim=2)
    if (present(nSolution)) then
      @:ASSERT(nSolution <= nrhs)
      nrhs = nSolution
    end if

    info = 0
    allocate(ipiv(nn))
    call sgesv(nn, nrhs, aa, lda, ipiv, bb, ldb, info)

    if (info < 0) then
99000 format ('Failure in linear equation solver sgesv,', &
          & ' illegal argument at position ',i10)
      write (error_string, 99000) info
      call error(error_string)
    else
      if (present(iError)) then
        iError = info
      elseif (info > 0) then
99010   format ('Linear dependent system in linear equation solver sgetrf,', &
            & ' info flag: ',i10)
        write (error_string, 99010) info
        call error(error_string)
      end if
    end if

  end subroutine gesv_real


  !> Double precision version of gesv
  subroutine gesv_dble(aa, bb, nEquation, nSolution, iError)

    !> Contains the coefficients on entry, the LU factorisation on exit.
    real(rdp), intent(inout) :: aa(:,:)

    !> Right hand side(s) of the linear equation on entry, solution(s) on exit.
    real(rdp), intent(inout) :: bb(:,:)

    !> The size of the problem (nr. of variables and equations). Must be only specified if different
    !> from size(aa, dim=1).
    integer, intent(in), optional :: nEquation

    !> Nr. of right hand sides (nr. of solutions). Must be only specified if different from size(b,
    !> dim=2).
    integer, intent(in), optional :: nSolution

    !> Error flag. If present, Lapack error flags are reported and noncritical errors (iError > 0)
    !> will not abort the program.
    integer, intent(out), optional :: iError

    integer :: info
    integer :: nn, nrhs, lda, ldb
    integer, allocatable :: ipiv(:)
    character(len=100) :: error_string

    lda = size(aa, dim=1)
    if (present(nEquation)) then
      @:ASSERT(nEquation >= 1 .and. nEquation <= lda)
      nn = nEquation
    else
      nn = lda
    end if
    @:ASSERT(size(aa, dim=2) >= nn)

    ldb = size(bb, dim=1)
    @:ASSERT(ldb >= nn)
    nrhs = size(bb, dim=2)
    if (present(nSolution)) then
      @:ASSERT(nSolution <= nrhs)
      nrhs = nSolution
    end if

    info = 0
    allocate(ipiv(nn))
    call dgesv(nn, nrhs, aa, lda, ipiv, bb, ldb, info)

    if (info < 0) then
99020 format ('Failure in linear equation solver dgesv,', &
          & ' illegal argument at position ',i10)
      write (error_string, 99020) info
      call error(error_string)
    else
      if (present(iError)) then
        iError = info
      elseif (info > 0) then
99030   format ('Linear dependent system in linear equation solver dgesv,', &
            & ' info flag: ',i10)
        write (error_string, 99030) info
        call error(error_string)
      end if
    end if

  end subroutine gesv_dble


  !> Double precision version of gesv
  subroutine gesv_dcomplex(aa, bb, nEquation, nSolution, iError)

    !> Contains the coefficients on entry, the LU factorisation on exit.
    complex(rdp), intent(inout) :: aa(:,:)

    !> Right hand side(s) of the linear equation on entry, solution(s) on exit.
    complex(rdp), intent(inout) :: bb(:,:)

    !> The size of the problem (nr. of variables and equations). Must be only specified if different
    !> from size(aa, dim=1).
    integer, intent(in), optional :: nEquation

    !> Nr. of right hand sides (nr. of solutions). Must be only specified if different from size(b,
    !> dim=2).
    integer, intent(in), optional :: nSolution

    !> Error flag. If present, Lapack error flags are reported and noncritical errors (iError > 0)
    !> will not abort the program.
    integer, intent(out), optional :: iError

    integer :: info
    integer :: nn, nrhs, lda, ldb
    integer, allocatable :: ipiv(:)
    character(len=100) :: error_string

    lda = size(aa, dim=1)
    if (present(nEquation)) then
      @:ASSERT(nEquation >= 1 .and. nEquation <= lda)
      nn = nEquation
    else
      nn = lda
    end if
    @:ASSERT(size(aa, dim=2) >= nn)

    ldb = size(bb, dim=1)
    @:ASSERT(ldb >= nn)
    nrhs = size(bb, dim=2)
    if (present(nSolution)) then
      @:ASSERT(nSolution <= nrhs)
      nrhs = nSolution
    end if

    info = 0
    allocate(ipiv(nn))
    call zgesv(nn, nrhs, aa, lda, ipiv, bb, ldb, info)

    if (info < 0) then
99020 format ('Failure in linear equation solver dgesv,', &
          & ' illegal argument at position ',i10)
      write (error_string, 99020) info
      call error(error_string)
    else
      if (present(iError)) then
        iError = info
      elseif (info > 0) then
99030   format ('Linear dependent system in linear equation solver dgesv,', &
            & ' info flag: ',i10)
        write (error_string, 99030) info
        call error(error_string)
      end if
    end if

  end subroutine gesv_dcomplex


  !> Single precision version of getrf.
  subroutine getrf_real(aa, ipiv, nRow, nColumn, iError)

    !> Matrix to decompose on entry, L and U on exit. Unit diagonal elements of L are not stored.
    real(rsp), intent(inout) :: aa(:,:)

    !> Pivot indices, row i of the matrix was interchanged with row ipiv(i).
    integer, intent(out) :: ipiv(:)

    !> Number of rows of the matrix to decomposea. (Necessary if different from the number of rows
    !> of the passed matrix)
    integer, optional, intent(in) :: nRow

    !> Number of rows of the matrix to decompose. (Necessary if different from the number of columns
    !> of the passed matrix)
    integer, optional, intent(in) :: nColumn

    !> Error flag. Zero on successful exit. If not present, any lapack error causes program
    !> termination. If passed only fatal lapack errors with error flag < 0 cause abort.
    integer, optional, intent(out) :: iError

    integer :: mm, nn, lda, info
    character(len=100) :: error_string

    lda = size(aa, dim=1)
    nn = size(aa, dim=2)
    if (present(nRow)) then
      @:ASSERT(nRow >= 1 .and. nRow <= lda)
      mm = nRow
    else
      mm = lda
    end if
    if (present(nColumn)) then
      @:ASSERT(nColumn >= 1 .and. nColumn <= nn)
      nn = nColumn
    end if
    @:ASSERT(size(ipiv) == min(mm, nn))

    call sgetrf(mm, nn, aa, lda, ipiv, info)

    if (info < 0) then
99040 format ('Failure in LU factorisation sgetrf,', &
          & ' illegal argument at position ',i10)
      write (error_string, 99040) info
      call error(error_string)
    else
      if (present(iError)) then
        iError = info
      elseif (info > 0) then
99050   format ('Factor U is exactly zero in sgetrf,', &
            & ' info flag is ',i10)
        write (error_string, 99050) info
        call error(error_string)
      end if
    end if

  end subroutine getrf_real


  !> Double precision version of getrf.
  subroutine getrf_dble(aa, ipiv, nRow, nColumn, iError)

    !> Matrix to decompose on entry, L and U on exit. Unit diagonal elements of L are not stored.
    real(rdp), intent(inout) :: aa(:,:)

    !> Pivot indices, row i of the matrix was interchanged with row ipiv(i).
    integer, intent(out) :: ipiv(:)

    !> Number of rows of the matrix to decomposea. (Necessary if different from the number of rows
    !> of the passed matrix)
    integer, optional, intent(in) :: nRow

    !> Number of rows of the matrix to decompose. (Necessary if different from the number of columns
    !> of the passed matrix)
    integer, optional, intent(in) :: nColumn

    !> Error flag. Zero on successful exit. If not present, any lapack error causes program
    !> termination. If passed only fatal lapack errors with error flag < 0 cause abort.
    integer, optional, intent(out) :: iError

    integer :: mm, nn, lda, info
    character(len=100) :: error_string

    lda = size(aa, dim=1)
    nn = size(aa, dim=2)
    if (present(nRow)) then
      @:ASSERT(nRow >= 1 .and. nRow <= lda)
      mm = nRow
    else
      mm = lda
    end if
    if (present(nColumn)) then
      @:ASSERT(nColumn >= 1 .and. nColumn <= nn)
      nn = nColumn
    end if
    @:ASSERT(size(ipiv) == min(mm, nn))

    call dgetrf(mm, nn, aa, lda, ipiv, info)

    if (info < 0) then
99060 format ('Failure in LU factorisation dgetrf,', &
          & ' illegal argument at position ',i10)
      write (error_string, 99060) info
      call error(error_string)
    else
      if (present(iError)) then
        iError = info
      elseif (info > 0) then
99070   format ('Factor U is exactly zero in dgetrf,', &
            & ' info flag is ',i10)
        write (error_string, 99070) info
        call error(error_string)
      end if
    end if

  end subroutine getrf_dble


  !> Complex precision version of getrf.
  subroutine getrf_complex(aa, ipiv, nRow, nColumn, iError)

    !> Matrix to decompose on entry, L and U on exit. Unit diagonal elements of L are not stored.
    complex(rsp), intent(inout) :: aa(:,:)

    !> Pivot indices, row i of the matrix was interchanged with row ipiv(i).
    integer, intent(out) :: ipiv(:)

    !> Number of rows of the matrix to decomposea. (Necessary if different from the number of rows
    !> of the passed matrix)
    integer, optional, intent(in) :: nRow

    !> Number of rows of the matrix to decompose. (Necessary if different from the number of columns
    !> of the passed matrix)
    integer, optional, intent(in) :: nColumn

    !> Error flag. Zero on successful exit. If not present, any lapack error causes program
    !> termination. If passed only fatal lapack errors with error flag < 0 cause abort.
    integer, optional, intent(out) :: iError

    integer :: mm, nn, lda, info
    character(len=100) :: error_string

    lda = size(aa, dim=1)
    nn = size(aa, dim=2)
    if (present(nRow)) then
      @:ASSERT(nRow >= 1 .and. nRow <= lda)
      mm = nRow
    else
      mm = lda
    end if
    if (present(nColumn)) then
      @:ASSERT(nColumn >= 1 .and. nColumn <= nn)
      nn = nColumn
    end if
    @:ASSERT(size(ipiv) == min(mm, nn))

    call cgetrf(mm, nn, aa, lda, ipiv, info)

    if (info < 0) then
99045 format ('Failure in LU factorisation sgetrf,', &
          & ' illegal argument at position ',i10)
      write (error_string, 99045) info
      call error(error_string)
    else
      if (present(iError)) then
        iError = info
      elseif (info > 0) then
99055   format ('Factor U is exactly zero in sgetrf,', &
            & ' info flag is ',i10)
        write (error_string, 99055) info
        call error(error_string)
      end if
    end if

  end subroutine getrf_complex


  !> Double precision version of getrf.
  subroutine getrf_dcomplex(aa, ipiv, nRow, nColumn, iError)

    !> Matrix to decompose on entry, L and U on exit. Unit diagonal elements of L are not stored.
    complex(rdp), intent(inout) :: aa(:,:)

    !> Pivot indices, row i of the matrix was interchanged with row ipiv(i).
    integer, intent(out) :: ipiv(:)

    !> Number of rows of the matrix to decomposea. (Necessary if different from the number of rows
    !> of the passed matrix)
    integer, optional, intent(in) :: nRow

    !> Number of rows of the matrix to decompose. (Necessary if different from the number of columns
    !> of the passed matrix)
    integer, optional, intent(in) :: nColumn

    !> Error flag. Zero on successful exit. If not present, any lapack error causes program
    !> termination. If passed only fatal lapack errors with error flag < 0 cause abort.
    integer, optional, intent(out) :: iError

    integer :: mm, nn, lda, info
    character(len=100) :: error_string

    lda = size(aa, dim=1)
    nn = size(aa, dim=2)
    if (present(nRow)) then
      @:ASSERT(nRow >= 1 .and. nRow <= lda)
      mm = nRow
    else
      mm = lda
    end if
    if (present(nColumn)) then
      @:ASSERT(nColumn >= 1 .and. nColumn <= nn)
      nn = nColumn
    end if
    @:ASSERT(size(ipiv) == min(mm, nn))

    call zgetrf(mm, nn, aa, lda, ipiv, info)

    if (info < 0) then
99065 format ('Failure in LU factorisation dgetrf,', &
          & ' illegal argument at position ',i10)
      write (error_string, 99065) info
      call error(error_string)
    else
      if (present(iError)) then
        iError = info
      elseif (info > 0) then
99075   format ('Factor U is exactly zero in dgetrf,', &
            & ' info flag is ',i10)
        write (error_string, 99075) info
        call error(error_string)
      end if
    end if

  end subroutine getrf_dcomplex


  !> Single precision version of getri.
  subroutine getri_real(aa, ipiv, nRow, iError)

    !> Matrix to decompose on entry, L and U on exit. Unit diagonal elements of L are not stored.
    real(rsp), intent(inout) :: aa(:,:)

    !> Pivot indices, as calculated by getri
    integer, intent(in) :: ipiv(:)

    !> Number of rows of the matrix to decompose. (Necessary if different from the number of rows of
    !> the passed matrix)
    integer, intent(in), optional :: nRow

    !> iError Error flag. Zero on successful exit. If not present, any lapack error causes program
    !> termination. If present, only fatal lapack errors with error flag < 0 cause abort.
    integer, intent(out), optional :: iError

    integer :: nn, lda, info, lwork
    real(rsp), allocatable :: work(:)
    real(rsp) :: work2(1)
    character(len=100) :: error_string

    lda = size(aa, dim=1)
    if (present(nRow)) then
      @:ASSERT(nRow >= 1 .and. nRow <= lda)
      nn = nRow
    else
      nn = lda
    end if
    @:ASSERT(size(aa, dim=2) >= nn)
    @:ASSERT(size(ipiv) == nn)

    lwork = -1
    call sgetri(nn, aa, lda, ipiv, work2, lwork, info)
    lwork = int(work2(1))

    allocate(work(lwork))
    call sgetri(nn, aa, lda, ipiv, work, lwork, info)

    if (info < 0) then
99080 format ('Failure in LU factorisation (sgetri),', &
          & ' illegal argument at position ',i10)
      write (error_string, 99080) info
      call error(error_string)
    else
      if (present(iError)) then
        iError = info
      elseif (info > 0) then
99090   format ('Factor U is exactly zero in sgetri,', &
            & ' info flag is ',i10)
        write (error_string, 99090) info
        call error(error_string)
      end if
    end if

  end subroutine getri_real


  !> Double precision version of getri.
  subroutine getri_dble(aa, ipiv, nRow, iError)

    !> Matrix to decompose on entry, L and U on exit. Unit diagonal elements of L are not stored.
    real(rdp), intent(inout) :: aa(:,:)

    !> Pivot indices, as calculated by getri
    integer, intent(in) :: ipiv(:)

    !> Number of rows of the matrix to decompose. (Necessary if different from the number of rows of
    !> the passed matrix)
    integer, intent(in), optional :: nRow

    !> iError Error flag. Zero on successful exit. If not present, any lapack error causes program
    !> termination. If present, only fatal lapack errors with error flag < 0 cause abort.
    integer, intent(out), optional :: iError

    integer :: nn, lda, info, lwork
    real(rdp), allocatable :: work(:)
    real(rdp) :: work2(1)
    character(len=100) :: error_string

    lda = size(aa, dim=1)
    if (present(nRow)) then
      @:ASSERT(nRow >= 1 .and. nRow <= lda)
      nn = nRow
    else
      nn = lda
    end if
    @:ASSERT(size(aa, dim=2) >= nn)
    @:ASSERT(size(ipiv) == nn)

    lwork = -1
    call dgetri(nn, aa, lda, ipiv, work2, lwork, info)
    lwork = int(work2(1))

    allocate(work(lwork))
    call dgetri(nn, aa, lda, ipiv, work, lwork, info)

    if (info < 0) then
99100 format ('Failure in LU factorisation (dgetri), illegal argument at&
          & position ', i10)
      write (error_string, 99100) info
      call error(error_string)
    else
      if (present(iError)) then
        iError = info
      elseif (info > 0) then
99110   format ('Factor U is exactly zero in dgetri,', &
            & ' info flag is ',i10)
        write (error_string, 99110) info
        call error(error_string)
      end if
    end if

  end subroutine getri_dble


  !> Inverts a matrix.
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


#:for TYPE, KIND, NAME, PRC, LABEL, PRF in [('sy', 'real', 'symmetric', 'rsp', 'real', 'ssytrf'),&
  & ('sy', 'real', 'symmetric', 'rdp', 'dreal', 'dsytrf'),&
  & ('he', 'complex', 'hermitian', 'rsp', 'cmplx', 'chetrf'),&
  & ('he', 'complex', 'hermitian', 'rdp', 'dcmplx', 'zhetrf')]

  !> Computes the Bunch-Kaufman factorization of a ${NAME}$ matrix.
  subroutine ${TYPE}$trf_${LABEL}$(aa, ipiv, status, uplo)

    !> ${NAME}$ matrix
    ${KIND}$(${PRC}$), intent(inout) :: aa(:,:)

    !> Interchanges of blocks on exit.
    integer, intent(out) :: ipiv(:)

    !> Status of operation
    type(TStatus), intent(out) :: status

    !> Signals whether upper (U) or lower (L) triangle should be used (default: lower).
    character, intent(in), optional :: uplo

    integer :: nn, info, lwork
    ${KIND}$(${PRC}$), allocatable :: work(:)
    ${KIND}$(${PRC}$) :: tmpwork(1)
    character :: uplo0

    uplo0 = uploHelper(uplo)
    nn = size(aa, dim=2)
    lwork = -1
    call ${PRF}$(uplo0, nn, aa, nn, ipiv, tmpwork, lwork, info)
    if (info /= 0) then
      @:RAISE_FORMATTED_ERROR(status, -1, "('Failure in ${PRF}$ memory check, info: ',I0)", info)
    end if
    lwork = int(tmpwork(1))
    allocate(work(lwork), stat=info)
    if (info /= 0) then
      @:RAISE_ERROR(status, -1, "Out of memory in ${PRF}$")
    end if
    call ${PRF}$(uplo0, nn, aa, nn, ipiv, work, lwork, info)
    if (info /= 0) then
      @:RAISE_FORMATTED_ERROR(status, -1, "('Failure in ${PRF}$, info: ',I0)", info)
    end if

  end subroutine ${TYPE}$trf_${LABEL}$

#:endfor

#:for TYPE, KIND, NAME, PRC, LABEL, PRF in [('sy', 'real', 'symmetric', 'rsp', 'real', 'ssytri'),&
  & ('sy', 'real', 'symmetric', 'rdp', 'dreal', 'dsytri'),&
  & ('he', 'complex', 'hermitian', 'rsp', 'cmplx', 'chetri'),&
  & ('he', 'complex', 'hermitian', 'rdp', 'dcmplx', 'zhetri')]

  !> Computes the inverse of a ${NAME}$ matrix.
  subroutine ${TYPE}$tri_${LABEL}$(aa, ipiv, status, uplo)

    !> ${NAME}$ matrix to be inverted.
    ${KIND}$(${PRC}$), intent(in) :: aa(:,:)

    !> Block interchanges as created by the sytrf() routine.
    integer, intent(in) :: ipiv(:)

    !> Status of operation
    type(TStatus), intent(out) :: status

    !> Upper ('U') or lower ('L') matrix (default: 'L')
    character, intent(in), optional :: uplo

    integer :: info, nn
    character :: uplo0
    ${KIND}$(${PRC}$), allocatable :: work(:)

    uplo0 = uploHelper(uplo)
    nn = size(aa, dim=1)
    allocate(work(max(1, 2 * nn)), stat=info)
    if (info /= 0) then
      @:RAISE_ERROR(status, -1, "Out of memory in ${PRF}$ wrapper")
    end if
    call ${PRF}$(uplo0, nn, aa, nn, ipiv, work, info)
    if (info /= 0) then
      @:RAISE_FORMATTED_ERROR(status, -1, "('Failure in ${PRF}$, info: ',I0)", info)
    end if

  end subroutine ${TYPE}$tri_${LABEL}$

#:endfor


  !> single precision version of larnv
  subroutine larnv_real(iDist,iSeed,x)

    !> choice of distribution (1: uniform (0,1), 2: uniform (-1,1), 3: normal (0,1)
    integer, intent(in) :: iDist

    !> On entry, the seed of the random number generator; the array elements must be between 0 and
    !> 4095, and ISEED(4) must be odd. On exit, the seed is updated.
    integer, intent(inout) :: iSeed(4)

    !> On exit, vector of random numbers
    real(rsp), intent(out) :: x(:)

    integer :: n

    @:ASSERT(iDist > 0)
    @:ASSERT(iDist < 4)
    @:ASSERT(all(iSeed(:) >= 0))
    @:ASSERT(all(iSeed(:) <= 4095))
    @:ASSERT(mod(iSeed(4),2) == 1)
    @:ASSERT(size(x) > 0)
    n = size(x)
    x(:) = 0.0
    call SLARNV( iDist, iSeed, n, x )
  end subroutine larnv_real


  !> double precision version of larnv
  subroutine larnv_dble(iDist,iSeed,x)

    !> choice of distribution (1: uniform (0,1), 2: uniform (-1,1), 3: normal (0,1)
    integer, intent(in) :: iDist

    !> On entry, the seed of the random number generator; the array elements must be between 0 and
    !> 4095, and ISEED(4) must be odd. On exit, the seed is updated.
    integer, intent(inout) :: iSeed(4)

    !> On exit, vector of random numbers
    real(rdp), intent(out) :: x(:)

    integer :: n

    @:ASSERT(iDist > 0)
    @:ASSERT(iDist < 4)
    @:ASSERT(all(iSeed(:) >= 0))
    @:ASSERT(all(iSeed(:) <= 4095))
    @:ASSERT(mod(iSeed(4),2) == 1)
    @:ASSERT(size(x) > 0)
    n = size(x)
    x(:) = 0.0d0
    call DLARNV( iDist, iSeed, n, x )
  end subroutine larnv_dble


  !> complex version of larnv
  subroutine larnv_cplx(iDist,iSeed,x)

    !> choice of distribution (1: uniform (0,1), 2: uniform (-1,1), 3: normal (0,1)
    integer, intent(in) :: iDist

    !> On entry, the seed of the random number generator; the array elements must be between 0 and
    !> 4095, and ISEED(4) must be odd. On exit, the seed is updated.
    integer, intent(inout) :: iSeed(4)

    !> On exit, vector of random numbers
    complex(rsp), intent(out) :: x(:)

    integer :: n

    @:ASSERT(iDist > 0)
    @:ASSERT(iDist < 4)
    @:ASSERT(all(iSeed(:) >= 0))
    @:ASSERT(all(iSeed(:) <= 4095))
    @:ASSERT(mod(iSeed(4),2) == 1)
    @:ASSERT(size(x) > 0)
    n = size(x)
    x(:) = 0.0
    call CLARNV( iDist, iSeed, n, x )
  end subroutine larnv_cplx


  !> double complex precision version of larnv
  subroutine larnv_dblecplx(iDist,iSeed,x)

    !> choice of distribution (1: uniform (0,1), 2: uniform (-1,1), 3: normal (0,1)
    integer, intent(in) :: iDist

    !> INTEGER array, dimension (4) On entry, the seed of the random number generator; the array
    !> elements must be between 0 and 4095, and ISEED(4) must be odd. On exit, the seed is updated.
    integer, intent(inout) :: iSeed(4)

    !> On exit, vector of random numbers
    complex(rdp), intent(out) :: x(:)

    integer :: n

    @:ASSERT(iDist > 0)
    @:ASSERT(iDist < 4)
    @:ASSERT(all(iSeed(:) >= 0))
    @:ASSERT(all(iSeed(:) <= 4095))
    @:ASSERT(mod(iSeed(4),2) == 1)
    @:ASSERT(size(x) > 0)
    n = size(x)
    x(:) = 0.0d0
    call ZLARNV( iDist, iSeed, n, x )
  end subroutine larnv_dblecplx


  !> real svd decomposition of matrix A into left and right vectors and singular values
  subroutine sgesvd_real(A,u,sigma,vt)

    !> matrix to decompose, warning the matrix is over-written by the routine
    real(rsp), intent(inout) :: A(:,:)

    !> first min(m,n) columns of u hold the left singular vector on return
    real(rsp), intent(out) :: u(:,:)

    !> holds the singular values on return
    real(rsp), intent(out) :: sigma(:)

    !> first min(m,n) columns of vt hold the right singular vector on return - warning this matrix
    !> is returned transpose(conjugated()) i.e. A = u.s.vt and all non-returned singular vectors are
    !> zero!
    real(rsp), intent(out) :: vt(:,:)

    integer :: n, m, mn, lda, lwork, ldu, ldvt, info
    real(rsp), allocatable :: work(:)
    character(len=100) :: error_string

    m = size(A,dim=1)
    n = size(A,dim=2)
    mn = min(m,n)
    lda = size(A,dim=1)
    ldu = size(U,dim=1)
    ldvt = size(Vt,dim=1)
    @:ASSERT(all(shape(u) == (/m,mn/)))
    @:ASSERT(all(shape(vt) == (/mn,n/)))
    @:ASSERT(size(sigma) == mn)

    lwork = max(1,3*min(m,n)+max(m,n),5*min(m,n))

    allocate(work(lwork))

    ! get only the minimum(m,n) singular vectors
    call sgesvd('S', 'S', m, n, A, lda, sigma, u, ldu, vt, ldvt, work, lwork, info)

    if (info /= 0) then
      write(error_string, "(A,I10)") "SVD failed. Info: ", info
      call error(error_string)
    end if

    deallocate(work)

  end subroutine sgesvd_real


  !> double precision svd decomposition of matrix A into left and right vectors and singular values
  subroutine dgesvd_dble(A,u,sigma,vt)

    !> matrix to decompose, warning the matrix is over-written by the routine
    real(rdp), intent(inout) :: A(:,:)

    !> first min(m,n) columns of u hold the left singular vector on return
    real(rdp), intent(out) :: u(:,:)

    !> holds the singular values on return
    real(rdp), intent(out) :: sigma(:)

    !> first min(m,n) columns of vt hold the right singular vector on return - warning this matrix
    !> is returned transpose(conjugated()) i.e. A = u.s.vt and all non-returned singular vectors are
    !> zero!
    real(rdp), intent(out) :: vt(:,:)

    integer :: n, m, mn, lda, lwork, ldu, ldvt, info
    real(rdp), allocatable :: work(:)
    character(len=100) :: error_string

    m = size(A,dim=1)
    n = size(A,dim=2)
    mn = min(m,n)
    lda = size(A,dim=1)
    ldu = size(U,dim=1)
    ldvt = size(Vt,dim=1)
    @:ASSERT(all(shape(u) == (/m,mn/)))
    @:ASSERT(all(shape(vt) == (/mn,n/)))
    @:ASSERT(size(sigma) == mn)

    lwork = max(1,3*min(m,n)+max(m,n),5*min(m,n))

    allocate(work(lwork))

    ! get only the minimum(m,n) singular vectors
    call dgesvd('S', 'S', m, n, A, lda, sigma, u, ldu, vt, ldvt, work, lwork, info)

    if (info /= 0) then
      write(error_string, "(A,I10)") "SVD failed. Info: ", info
      call error(error_string)
    end if

    deallocate(work)

  end subroutine dgesvd_dble

  !> complex svd decomposition of matrix A into left and right vectors and singular values
  subroutine cgesvd_cplx(A,u,sigma,vt)

    !> matrix to decompose, warning the matrix is over-written by the routine
    complex(rsp), intent(inout) :: A(:,:)

    !> first min(m,n) columns of u hold the left singular vector on return
    complex(rsp), intent(out) :: u(:,:)

    !> holds the singular values on return
    real(rsp), intent(out) :: sigma(:)

    !> first min(m,n) columns of vt hold the right singular vector on return - warning this matrix
    !> is returned transpose(conjugated()) i.e. A = u.s.vt and all non-returned singular vectors are
    !> zero!
    complex(rsp), intent(out) :: vt(:,:)

    integer :: n, m, mn, lda, lwork, ldu, ldvt, info
    real(rsp), allocatable :: rwork(:)
    complex(rsp), allocatable :: work(:)
    character(len=100) :: error_string

    m = size(A,dim=1)
    n = size(A,dim=2)
    mn = min(m,n)
    lda = size(A,dim=1)
    ldu = size(U,dim=1)
    ldvt = size(Vt,dim=1)
    @:ASSERT(all(shape(u) == (/m,mn/)))
    @:ASSERT(all(shape(vt) == (/mn,n/)))
    @:ASSERT(size(sigma) == mn)

    lwork = 2*min(m,n)+max(m,n)

    allocate(rwork(5*mn))
    allocate(work(lwork))

    ! get only the minimum(m,n) singular vectors
    call cgesvd('S', 'S', m, n, A, lda, sigma, u, ldu, vt, ldvt, work, lwork, rwork, info)

    if (info /= 0) then
      write(error_string, "(A,I10)") "SVD failed. Info: ", info
      call error(error_string)
    end if

    deallocate(rwork)
    deallocate(work)

  end subroutine cgesvd_cplx


  !> double complex svd decomposition of matrix A into left and right vectors and singular values
  subroutine zgesvd_dblecplx(A,u,sigma,vt)

    !> matrix to decompose, warning the matrix is over-written by the routine
    complex(rdp), intent(inout) :: A(:,:)

    !> first min(m,n) columns of u hold the left singular vector on return
    complex(rdp), intent(out) :: u(:,:)

    !> holds the singular values on return
    real(rdp), intent(out) :: sigma(:)

    !> first min(m,n) columns of vt hold the right singular vector on return - warning this matrix
    !> is returned transpose(conjugated()) i.e. A = u.s.vt and all non-returned singular vectors are
    !> zero!
    complex(rdp), intent(out) :: vt(:,:)

    integer :: n, m, mn, lda, lwork, ldu, ldvt, info
    real(rdp), allocatable :: rwork(:)
    complex(rdp), allocatable :: work(:)
    character(len=100) :: error_string

    m = size(A,dim=1)
    n = size(A,dim=2)
    mn = min(m,n)
    lda = size(A,dim=1)
    ldu = size(U,dim=1)
    ldvt = size(Vt,dim=1)
    @:ASSERT(all(shape(u) == (/m,mn/)))
    @:ASSERT(all(shape(vt) == (/mn,n/)))
    @:ASSERT(size(sigma) == mn)

    lwork = 2*min(m,n)+max(m,n)

    allocate(rwork(5*mn))
    allocate(work(lwork))

    ! get only the minimum(m,n) singular vectors
    call zgesvd('S', 'S', m, n, A, lda, sigma, u, ldu, vt, ldvt, work, lwork, rwork, info)

    if (info /= 0) then
      write(error_string, "(A,I10)") "SVD failed. Info: ", info
      call error(error_string)
    end if

    deallocate(rwork)
    deallocate(work)

  end subroutine zgesvd_dblecplx


#:for suffix, kind in REAL_KIND_PARAMS

  !> Choleskii factorization of a matrix
  subroutine ${kind}$potrf_${suffix}$(b, uplo, info)

    !> Matrix to be factorised, over-written on return
    real(r${kind}$p), intent(inout) :: b(:,:)

    !> upper or lower triangle of the matrix, defaults to lower
    character, intent(in), optional :: uplo

    !> Info flag. If not present and an error occurs, the subroutine stops.
    integer, intent(out), optional :: info

    integer :: info0, n, ldb
    character :: uplo0
    character(len=100) :: error_string

    uplo0 = uploHelper(uplo)
    n = size(b, dim=2)
    ldb = size(b, dim=1)
  @:ASSERT(ldb >= n)

    call ${kind}$potrf(uplo0, n, b, ldb, info0)
    if (present(info)) then
      info = info0
    elseif (info0 /= 0) then
      write(error_string, "(A,I10)") "Routine ${kind}$potrf failed. Info: ", info0
      call error(error_string)
    end if

  end subroutine ${kind}$potrf_${suffix}$

#:endfor


#:for suffix, kind in REAL_KIND_PARAMS

  !> solve one of the matrix equations op( A )*X = alpha*B, or X*op( A ) = alpha*B
  subroutine ${kind}$trsm_${suffix}$(side, A, B, m, n, diag, alpha, transa, uplo)

    !> matrix A on 'l'eft or 'r'ight of X
    character, intent(in) :: side

    real(r${kind}$p), intent(inout) :: A(:,:)

    real(r${kind}$p), intent(inout) :: B(:,:)

    integer, intent(in) :: m

    integer, intent(in) :: n

    !> 'U'nit triangular or 'N'ot
    character, intent(in) :: diag

    real(r${kind}$p), intent(in) :: alpha

    !> optional transpose of A matrix (defaults to 'n'), allowed choices are 'n', 'N', 't', 'T', 'c'
    !> and 'C'
    character, intent(in), optional :: transA

    !> upper or lower triangle of the matrix, defaults to lower
    character, intent(in), optional :: uplo

    integer :: lda, ldb
    character :: uplo0, iTransA

    uplo0 = uploHelper(uplo)
    lda = size(A, dim=1)
    ldb = size(B, dim=1)
  @:ASSERT(m > 0)
  @:ASSERT(n > 0)
  @:ASSERT(((side == 'r' .or. side == 'R') .and. lda > n) .or.&
      & ((side == 'l' .or. side == 'L') .and. lda > m))
  @:ASSERT(ldb >= m)
  @:ASSERT(size(B,dim=2) >= n)
  @:ASSERT(side == 'r' .or. side == 'R' .or. side == 'l' .or. side == 'L')
  @:ASSERT(diag == 'u' .or. diag == 'U' .or. diag == 'n' .or. diag == 'N')
  if (present(transa)) then
    iTransA = transA
  else
    iTransA = 'n'
  end if
  @:ASSERT(iTransA == 'n' .or. iTransA == 'N' .or. iTransA == 't'&
      & .or. iTransA == 'T' .or. iTransA == 'c' .or. iTransA == 'C')

    call ${kind}$trsm ( side, uplo, iTransa, diag, m, n, alpha, a, lda, b, ldb )

  end subroutine ${kind}$trsm_${suffix}$

#:endfor


  !> Helper function for matrix triangle options to choose optional triangle
  pure function uploHelper(uplo)

    !> upper or lower triangle of the matrix, defaults to lower if not present
    character, intent(in), optional :: uplo

    !> Resulting triangle to use
    character :: uploHelper

    if (present(uplo)) then
      uploHelper = uplo
    else
      uploHelper = "L"
    end if

  end function uploHelper


  !> Solves a system of linear equations with multiple right hand sides
  subroutine getrs_dble(amat, ipiv, bmat, trans, iError)

    !> Matrix of the linear system
    real(rdp), intent(in) :: amat(:, :)

    !> Pivot indices, row i of the matrix was interchanged with row ipiv(i).
    integer, intent(in) :: ipiv(:)

    !> Matrix of the right hand side vectors
    real(rdp), intent(inout) :: bmat(:, :)

    !> Optional transpose (defaults to 'n')
    character(len=1), intent(in), optional :: trans

    !> Error flag, zero on successful exit
    integer, intent(out), optional :: iError

    character(len=1) :: atr
    integer :: info, nn, nrhs, lda, ldb

    @:ASSERT(size(amat, 1) == size(amat, dim=2))
    @:ASSERT(size(amat, 1) == size(bmat, dim=1))
    if(present(trans)) then
      @:ASSERT(any(trans == ['n', 'N', 't', 'T', 'c', 'C']))
      atr = trans
    else
      atr = 'n'
    endif
    lda = max(1, size(amat, 1))
    ldb = max(1, size(bmat, 1))
    nn = size(amat, 2)
    nrhs = size(bmat, 2)
    call dgetrs(atr, nn, nrhs, amat, lda, ipiv, bmat, ldb, info)
    if(present(iError)) then
      iError = info
    else
      if (info /= 0) then
        call error("Failed to solve linear system by diagonal pivoting")
      end if
    endif

  end subroutine getrs_dble


  !> Solves a system of linear equations with one right hand sides
  subroutine getrs1_dble(amat, ipiv, bvec, trans, iError)

    !> Matrix of the linear system
    real(rdp), intent(in) :: amat(:, :)

    !> Pivot indices, row i of the matrix was interchanged with row ipiv(i).
    integer, intent(in) :: ipiv(:)

    !> Right hand side vector
    real(rdp), intent(inout), target :: bvec(:)

    !> optional transpose (defaults to 'n')
    character(len=1), intent(in), optional :: trans

    !> Error flag, zero on successful exit
    integer, intent(out), optional :: iError

    real(rdp), pointer :: bptr(:, :)

    bptr(1:size(bvec, 1), 1:1) => bvec(1:size(bvec, 1))
    call getrs(amat, ipiv, bptr, trans, iError)

  end subroutine getrs1_dble


  !> Solves a system of linear equations with multiple right hand sides
  subroutine getrs_real(amat, ipiv, bmat, trans, iError)

    !> Matrix of the linear system
    real(rsp), intent(in) :: amat(:, :)

    !> Pivot indices, row i of the matrix was interchanged with row ipiv(i).
    integer, intent(in) :: ipiv(:)

    !> Matrix of the right hand side vectors
    real(rsp), intent(inout) :: bmat(:, :)

    !> Optional transpose (defaults to 'n')
    character(len=1), intent(in), optional :: trans

    !> Error flag, zero on successful exit
    integer, intent(out), optional :: iError

    character(len=1) :: atr
    integer :: info, nn, nrhs, lda, ldb

    @:ASSERT(size(amat, 1) == size(amat, dim=2))
    @:ASSERT(size(amat, 1) == size(bmat, dim=1))
    if(present(trans)) then
      @:ASSERT(any(trans == ['n', 'N', 't', 'T', 'c', 'C']))
      atr = trans
    else
      atr = 'n'
    endif
    lda = max(1, size(amat, 1))
    ldb = max(1, size(bmat, 1))
    nn = size(amat, 2)
    nrhs = size(bmat, 2)
    call sgetrs(atr, nn, nrhs, amat, lda, ipiv, bmat, ldb, info)
    if(present(iError)) then
      iError = info
    else
      if (info /= 0) then
        call error("Failed to solve linear system by diagonal pivoting")
      end if
    endif

  end subroutine getrs_real


  !> Solves a system of linear equations with one right hand sides
  subroutine getrs1_real(amat, ipiv, bvec, trans, iError)

    !> Matrix of the linear system
    real(rsp), intent(in) :: amat(:, :)

    !> Pivot indices, row i of the matrix was interchanged with row ipiv(i).
    integer, intent(in) :: ipiv(:)

    !> Right hand side vector
    real(rsp), intent(inout), target :: bvec(:)

    !> Optional transpose (defaults to 'n')
    character(len=1), intent(in), optional :: trans

    !> Error flag, zero on successful exit
    integer, intent(out), optional :: iError

    real(rsp), pointer :: bptr(:, :)

    bptr(1:size(bvec, 1), 1:1) => bvec(1:size(bvec, 1))
    call getrs(amat, ipiv, bptr, trans, iError)

  end subroutine getrs1_real


end module dftbp_math_lapackroutines
