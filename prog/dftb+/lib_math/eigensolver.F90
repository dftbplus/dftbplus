!--------------------------------------------------------------------------------------------------!
!  DFTB+: general package for performing fast atomistic simulations                                !
!  Copyright (C) 2006 - 2020  DFTB+ developers group                                               !
!                                                                                                  !
!  See the LICENSE file for terms of usage and distribution.                                       !
!--------------------------------------------------------------------------------------------------!

#:include 'common.fypp'
#:include "error.fypp"

!> Contains F90 wrapper functions for some commonly used lapack calls needed in the code.
!> Contains some fixes for lapack 3.0 bugs, if this gets corrected in lapack 4.x they should be
!> removed.
module dftbp_eigensolver
  use dftbp_assert
  use dftbp_message
  use dftbp_accuracy, only : rsp, rdp
  use dftbp_blas
  use dftbp_lapack
#:if WITH_GPU
  use dftbp_magma,  only : magmaf_ssygvd_m, magmaf_dsygvd_m, magmaf_chegvd_m, magmaf_zhegvd_m
#:endif
  implicit none
  private

  public :: heev, hegv, hegvd, gvr, bgv, geev
#:if WITH_GPU
  public :: gpu_gvd
#:endif


  !> Simple eigensolver for a symmetric/Hermitian matrix
  !> Caveat: the matrix a is overwritten
  interface heev
    module procedure real_ssyev
    module procedure dble_dsyev
    module procedure cmplx_cheev
    module procedure dblecmplx_zheev
  end interface heev


  !> Simple eigensolver for a symmetric/Hermitian generalized matrix problem
  !> caveat: the matrix a is overwritten
  !> caveat: the matrix b is overwritten with Cholesky factorization
  interface hegv
    module procedure real_ssygv
    module procedure dble_dsygv
    module procedure cmplx_chegv
    module procedure dblecmplx_zhegv
  end interface hegv


  !> Simple eigensolver for a symmetric/Hermitian generalized matrix problem using divide and
  !> conquer eigensolver
  !> caveat: the matrix a is overwritten
  !> caveat: the matrix b is overwritten with Cholesky factorization
  interface hegvd
    module procedure real_ssygvd
    module procedure dble_dsygvd
    module procedure cmplx_chegvd
    module procedure dblecmplx_zhegvd
  end interface hegvd


  !> Simple eigensolver for a symmetric/Hermitian generalized matrix problem using the lapack
  !> relatively robust representation solver, based on the SYGV source. If the requested number of
  !> eigenvalues is lower than the size of H/S suspace mode is used (optionally the range can be set
  !> using il and ul) to return the lowest eigenvalues/vectors of number size(w)
  interface gvr
    module procedure real_ssygvr
    module procedure dble_dsygvr
    module procedure cmplx_chegvr
    module procedure dblecmplx_zhegvr
  end interface


  !> Eigensolver for a symmetric/Hermitian banded generalized matrix
  !> problem of the form A*x=(lambda)*B*x
  interface bgv
    module procedure real_ssbgv
    module procedure dble_dsbgv
    module procedure cmplx_chbgv
    module procedure dblecmplx_zhbgv
  end interface

#:if WITH_GPU
  !> Divide and conquer MAGMA GPU eigensolver
  interface gpu_gvd
    module procedure real_magma_ssygvd
    module procedure dble_magma_dsygvd
    module procedure cmplx_magma_chegvd
    module procedure dblecmplx_magma_zhegvd
  end interface
#:endif

  !> Simple eigensolver for a general matrix
  interface geev
    module procedure real_sgeev
    module procedure dble_dgeev
  end interface geev


contains

  !> Real eigensolver for a symmetric matrix
  subroutine real_ssyev(a,w,uplo,jobz)

    !> contains the matrix for the solver, returns as eigenvectors if requested (matrix always
    !> overwritten on return anyway)
    real(rsp), intent(inout) :: a(:,:)

    !> eigenvalues
    real(rsp), intent(out) :: w(:)

    !> upper or lower triangle of the matrix
    character, intent(in) :: uplo

    !> compute eigenvalues 'N' or eigenvalues and eigenvectors 'V'
    character, intent(in) :: jobz

    real(rsp), allocatable :: work(:)
    integer n, info
    integer :: int_idealwork
    real(rsp) :: idealwork(1)
    character(len=100) :: error_string

    @:ASSERT(uplo == 'u' .or. uplo == 'U' .or. uplo == 'l' .or. uplo == 'L')
    @:ASSERT(jobz == 'n' .or. jobz == 'N' .or. jobz == 'v' .or. jobz == 'V')
    @:ASSERT(all(shape(a)==size(w,dim=1)))
    n=size(a,dim=1)
    @:ASSERT(n>0)
    call ssyev(jobz, uplo, n, a, n, w, idealwork, -1, info)
    if (info/=0) then
      call error("Failure in SSYEV to determine optimum workspace")
    endif
    int_idealwork=floor(idealwork(1))
    allocate(work(int_idealwork))
    call ssyev(jobz, uplo, n, a, n, w, work, int_idealwork, info)
    if (info/=0) then
      if (info<0) then
99000 format ('Failure in diagonalisation routine ssyev,', &
          & ' illegal argument at position ',i6)
        write(error_string, 99000) info
        call error(error_string)
      else
99010 format ('Failure in diagonalisation routine ssyev,', &
          & ' diagonal element ',i6,' did not converge to zero.')
        write(error_string, 99010) info
        call error(error_string)
      endif
    endif

  end subroutine real_ssyev


  !> Double precision eigensolver for a symmetric matrix
  subroutine dble_dsyev(a,w,uplo,jobz)

    !> contains the matrix for the solver, returns as eigenvectors if requested (matrix always
    !> overwritten on return anyway)
    real(rdp), intent(inout) :: a(:,:)

    !> eigenvalues
    real(rdp), intent(out) :: w(:)

    !> upper or lower triangle of the matrix
    character, intent(in) :: uplo

    !> compute eigenvalues 'N' or eigenvalues and eigenvectors 'V'
    character, intent(in) :: jobz

    real(rdp), allocatable :: work(:)
    integer n, info
    integer :: int_idealwork
    real(rdp) :: idealwork(1)
    character(len=100) :: error_string

    @:ASSERT(uplo == 'u' .or. uplo == 'U' .or. uplo == 'l' .or. uplo == 'L')
    @:ASSERT(jobz == 'n' .or. jobz == 'N' .or. jobz == 'v' .or. jobz == 'V')
    @:ASSERT(all(shape(a)==size(w,dim=1)))
    n=size(a,dim=1)
    @:ASSERT(n>0)
    call dsyev(jobz, uplo, n, a, n, w, idealwork, -1, info)
    if (info/=0) then
      call error("Failure in DSYEV to determine optimum workspace")
    endif
    int_idealwork=floor(idealwork(1))
    allocate(work(int_idealwork))
    call dsyev(jobz, uplo, n, a, n, w, work, int_idealwork, info)
    if (info/=0) then
      if (info<0) then
99020 format ('Failure in diagonalisation routine dsyev,', &
          & ' illegal argument at position ',i6)
        write(error_string, 99020) info
        call error(error_string)
      else
99030 format ('Failure in diagonalisation routine dsyev,', &
          & ' diagonal element ',i6,' did not converge to zero.')
        write(error_string, 99030) info
        call error(error_string)
      endif
    endif

  end subroutine dble_dsyev


  !> Complex eigensolver for a Hermitian matrix
  subroutine cmplx_cheev(a,w,uplo,jobz)

    !> contains the matrix for the solver, returns as eigenvectors if requested (matrix always
    !> overwritten on return anyway)
    complex(rsp), intent(inout) :: a(:,:)

    !> eigenvalues
    real(rsp), intent(out) :: w(:)

    !> upper or lower triangle of the matrix
    character, intent(in) :: uplo

    !> compute eigenvalues 'N' or eigenvalues and eigenvectors 'V'
    character, intent(in) :: jobz

    real(rsp), allocatable :: rwork(:)
    complex(rsp), allocatable :: work(:)
    integer n, info
    integer :: int_idealwork
    complex(rsp) :: idealwork(1)
    character(len=100) :: error_string

    @:ASSERT(uplo == 'u' .or. uplo == 'U' .or. uplo == 'l' .or. uplo == 'L')
    @:ASSERT(jobz == 'n' .or. jobz == 'N' .or. jobz == 'v' .or. jobz == 'V')
    @:ASSERT(all(shape(a)==size(w,dim=1)))
    n=size(a,dim=1)
    @:ASSERT(n>0)
    allocate(rwork(3*n-2))
    call cheev(jobz, uplo, n, a, n, w, idealwork, -1, rwork, info)
    if (info/=0) then
      call error("Failure in CHEEV to determine optimum workspace")
    endif
    int_idealwork=floor(real(idealwork(1)))
    allocate(work(int_idealwork))
    call cheev(jobz, uplo, n, a, n, w, work, int_idealwork, rwork, info)
    if (info/=0) then
      if (info<0) then
99040 format ('Failure in diagonalisation routine cheev,', &
          & ' illegal argument at position ',i6)
        write(error_string, 99040) info
        call error(error_string)
      else
99050 format ('Failure in diagonalisation routine cheev,', &
          & ' diagonal element ',i6,' did not converge to zero.')
        write(error_string, 99050) info
        call error(error_string)
      endif
    endif

  end subroutine cmplx_cheev


  !> Double complex eigensolver for a Hermitian matrix
  subroutine dblecmplx_zheev(a,w,uplo,jobz)

    !> contains the matrix for the solver, returns as eigenvectors if requested (matrix always
    !> overwritten on return anyway)
    complex(rdp), intent(inout) :: a(:,:)

    !> eigenvalues
    real(rdp), intent(out) :: w(:)

    !> upper or lower triangle of the matrix
    character, intent(in) :: uplo

    !> compute eigenvalues 'N' or eigenvalues and eigenvectors 'V'
    character, intent(in) :: jobz

    real(rdp), allocatable :: rwork(:)
    complex(rdp), allocatable :: work(:)
    integer n, info
    integer :: int_idealwork
    complex(rdp) :: idealwork(1)
    character(len=100) :: error_string

    @:ASSERT(uplo == 'u' .or. uplo == 'U' .or. uplo == 'l' .or. uplo == 'L')
    @:ASSERT(jobz == 'n' .or. jobz == 'N' .or. jobz == 'v' .or. jobz == 'V')
    @:ASSERT(all(shape(a)==size(w,dim=1)))
    n=size(a,dim=1)
    @:ASSERT(n>0)
    allocate(rwork(3*n-2))
    call zheev(jobz, uplo, n, a, n, w, idealwork, -1, rwork, info)
    if (info/=0) then
      call error("Failure in ZHEEV to determine optimum workspace")
    endif
    int_idealwork=floor(real(idealwork(1)))
    allocate(work(int_idealwork))
    call zheev(jobz, uplo, n, a, n, w, work, int_idealwork, rwork, info)
    if (info/=0) then
      if (info<0) then
99060 format ('Failure in diagonalisation routine zheev,', &
          & ' illegal argument at position ',i6)
        write(error_string, 99060) info
        call error(error_string)
      else
99070 format ('Failure in diagonalisation routine zheev,', &
          & ' diagonal element ',i6,' did not converge to zero.')
        write(error_string, 99070) info
        call error(error_string)
      endif
    endif

  end subroutine dblecmplx_zheev


#:for NAME, VTYPE, RP in [("Double", "d", "dble"),("Single", "s", "real")]
  !> ${NAME}$ precision eigensolver for generalized symmetric matrix problem
  subroutine ${RP}$_${VTYPE}$sygv(a, b, w, uplo, jobz, itype)

    !> contains the matrix for the solver, returns eigenvectors if requested (matrix always
    !> overwritten on return anyway)
    real(r${VTYPE}$p), intent(inout) :: a(:,:)

    !> contains the second matrix for the solver (overwritten by Cholesky factorization)
    real(r${VTYPE}$p), intent(inout) :: b(:,:)

    !> eigenvalues
    real(r${VTYPE}$p), intent(out) :: w(:)

    !> upper or lower triangle of both matrices
    character, intent(in) :: uplo

    !> compute eigenvalues 'N' or eigenvalues and eigenvectors 'V'
    character, intent(in) :: jobz

    !> specifies the problem type to be solved 1:A*x=(lambda)*B*x, 2:A*B*x=(lambda)*x,
    !> 3:B*A*x=(lambda)*x default is 1
    integer, optional, intent(in) :: itype

    real(r${VTYPE}$p), allocatable :: work(:)
    integer n, lda, info, iitype, ldb
    integer :: int_idealwork
    real(r${VTYPE}$p) :: idealwork(1)
    character(len=100) :: error_string

  @:ASSERT(uplo == 'u' .or. uplo == 'U' .or. uplo == 'l' .or. uplo == 'L')
  @:ASSERT(jobz == 'n' .or. jobz == 'N' .or. jobz == 'v' .or. jobz == 'V')
    n = size(a,dim=2)
  @:ASSERT(n>0)
  @:ASSERT(all(shape(b) >= n))
  @:ASSERT(size(w) >= n)
    lda = size(a,dim=1)
  @:ASSERT(lda >= n)
    ldb = size(b,dim=1)
    if (present(itype)) then
      iitype = itype
    else
      iitype = 1
    end if
  @:ASSERT(iitype >= 1 .and. iitype <= 3 )
    call ${VTYPE}$sygv(iitype, jobz, uplo, n, a, lda, b, ldb, w, idealwork, -1, info)
    if (info/=0) then
       call error("Failure in ${VTYPE}$sygv to determine optimum workspace")
    endif
    int_idealwork=nint(idealwork(1))
    allocate(work(int_idealwork))
    call ${VTYPE}$sygv(iitype, jobz, uplo, n, a, lda, b, ldb, w, work, int_idealwork, info)
    if (info/=0) then
      if (info<0) then
        write(error_string, "('Failure in diagonalisation routine ${VTYPE}$sygv, illegal ',&
            & 'argument at position ',i6)") info
        call error(error_string)
      else if (info <= n) then
        write(error_string, "('Failure in diagonalisation routine ${VTYPE}$sygv, diagonal ',&
            & 'element ', i6, ' did not converge to zero.')") info
      else
        write(error_string, "('Failure in diagonalisation routine ${VTYPE}$sygv,', &
            & ' non-positive definite overlap! Minor ',i6,' responsible.')") info - n
        call error(error_string)
      endif
    endif

  end subroutine ${RP}$_${VTYPE}$sygv
#:endfor


  !> Complex eigensolver for generalized Hermitian matrix problem
  subroutine cmplx_chegv(a,b,w,uplo,jobz,itype)

    !> contains the matrix for the solver, returns eigenvectors if requested (matrix always
    !> overwritten on return anyway)
    complex(rsp), intent(inout) :: a(:,:)

    !> contains the second matrix for the solver (overwritten by Cholesky factorization)
    complex(rsp), intent(inout) :: b(:,:)

    !> eigenvalues
    real(rsp), intent(out) :: w(:)

    !> upper or lower triangle of both matrices
    character, intent(in) :: uplo

    !> compute eigenvalues 'N' or eigenvalues and eigenvectors 'V'
    character, intent(in) :: jobz

    !> specifies the problem type to be solved 1:A*x=(lambda)*B*x, 2:A*B*x=(lambda)*x,
    !> 3:B*A*x=(lambda)*x default is 1
    integer, optional, intent(in) :: itype

    complex(rsp), allocatable :: work(:)
    real(rsp), allocatable :: rwork(:)
    integer n, info, iitype, int_idealwork
    complex(rsp) :: idealwork(1)
    character(len=100) :: error_string

    @:ASSERT(uplo == 'u' .or. uplo == 'U' .or. uplo == 'l' .or. uplo == 'L')
    @:ASSERT(jobz == 'n' .or. jobz == 'N' .or. jobz == 'v' .or. jobz == 'V')
    @:ASSERT(all(shape(a)==shape(b)))
    @:ASSERT(all(shape(a)==size(w,dim=1)))
    n=size(a,dim=1)
    @:ASSERT(n>0)
    if (present(itype)) then
      iitype = itype
    else
      iitype = 1
    end if
    @:ASSERT(iitype >= 1 .and. iitype <= 3 )
    allocate(rwork(3*n-2))
    call CHEGV(iitype, jobz, uplo, n, a, n, b, n, w, idealwork, -1, rwork, info)
    if (info/=0) then
       call error("Failure in CHEGV to determine optimum workspace")
    endif
    int_idealwork=floor(real(idealwork(1)))
    allocate(work(int_idealwork))
    ! A*x = (lambda)*B*x upper triangles to be used
    call chegv(iitype, 'V', 'L', n, a, n, b, n, w, work, int_idealwork, rwork, info)
    if (info/=0) then
       if (info<0) then
99220 format ('Failure in diagonalisation routine chegv,', &
          & ' illegal argument at position ',i6)
          write(error_string, 99220) info
          call error(error_string)
       else if (info <= n) then
99230 format ('Failure in diagonalisation routine chegv,', &
          & ' diagonal element ',i6,' did not converge to zero.')
          write(error_string, 99230) info
          call error(error_string)
       else
99240 format ('Failure in diagonalisation routine chegv,', &
          & ' non-positive definite overlap! Minor ',i6,' responsible.')
          write(error_string, 99240) info - n
          call error(error_string)
       endif
    endif

  end subroutine cmplx_chegv


  !> Double complex eigensolver for generalized Hermitian matrix problem
  subroutine dblecmplx_zhegv(a,b,w,uplo,jobz,itype)

    !> contains the matrix for the solver, returns eigenvectors if requested (matrix always
    !> overwritten on return anyway)
    complex(rdp), intent(inout) :: a(:,:)

    !> contains the second matrix for the solver (overwritten by Cholesky factorization)
    complex(rdp), intent(inout) :: b(:,:)

    !> eigenvalues
    real(rdp), intent(out) :: w(:)

    !> upper or lower triangle of both matrices
    character, intent(in) :: uplo

    !> compute eigenvalues 'N' or eigenvalues and eigenvectors 'V'
    character, intent(in) :: jobz

    !> specifies the problem type to be solved 1:A*x=(lambda)*B*x, 2:A*B*x=(lambda)*x,
    !> 3:B*A*x=(lambda)*x default is 1
    integer, optional, intent(in) :: itype

    complex(rdp), allocatable :: work(:)
    real(rdp), allocatable :: rwork(:)
    integer n, info, iitype, int_idealwork
    complex(rdp) :: idealwork(1)
    character(len=100) :: error_string

    @:ASSERT(uplo == 'u' .or. uplo == 'U' .or. uplo == 'l' .or. uplo == 'L')
    @:ASSERT(jobz == 'n' .or. jobz == 'N' .or. jobz == 'v' .or. jobz == 'V')
    @:ASSERT(all(shape(a)==shape(b)))
    @:ASSERT(all(shape(a)==size(w,dim=1)))
    n=size(a,dim=1)
    @:ASSERT(n>0)
    if (present(itype)) then
      iitype = itype
    else
      iitype = 1
    end if
    @:ASSERT(iitype >= 1 .and. iitype <= 3 )
    allocate(rwork(3*n-2))
    call ZHEGV(iitype, jobz, uplo, n, a, n, b, n, w, idealwork, -1, rwork, info)
    if (info/=0) then
       call error("Failure in CHEGV to determine optimum workspace")
    endif
    int_idealwork=floor(real(idealwork(1)))
    allocate(work(int_idealwork))
    ! A*x = (lambda)*B*x upper triangles to be used
    call zhegv(iitype, 'V', 'L', n, a, n, b, n, w, work, int_idealwork, rwork, info)
    if (info/=0) then
       if (info<0) then
99250 format ('Failure in diagonalisation routine zhegv,', &
          & ' illegal argument at position ',i6)
          write(error_string, 99250) info
          call error(error_string)
       else if (info <= n) then
99260 format ('Failure in diagonalisation routine zhegv,', &
          & ' diagonal element ',i6,' did not converge to zero.')
          write(error_string, 99260) info
          call error(error_string)
       else
99270 format ('Failure in diagonalisation routine zhegv,', &
          & ' non-positive definite overlap! Minor ',i6,' responsible.')
          write(error_string, 99270) info - n
          call error(error_string)
       endif
    endif

  end subroutine dblecmplx_zhegv


  !> Real eigensolver for generalized symmetric matrix problem - divide and conquer
  subroutine real_ssygvd(a,b,w,uplo,jobz,itype)

    !> contains the matrix for the solver, returns eigenvectors if requested (matrix always
    !> overwritten on return anyway)
    real(rsp), intent(inout) :: a(:,:)

    !> contains the second matrix for the solver (overwritten by Cholesky factorization)
    real(rsp), intent(inout) :: b(:,:)

    !> eigenvalues
    real(rsp), intent(out) :: w(:)

    !> upper or lower triangle of the matrix
    character, intent(in) :: uplo

    !> compute eigenvalues 'N' or eigenvalues and eigenvectors 'V'
    character, intent(in) :: jobz

    !> optional specifies the problem type to be solved 1:A*x=(lambda)*B*x, 2:A*B*x=(lambda)*x,
    !> 3:B*A*x=(lambda)*x default is 1
    integer, optional, intent(in) :: itype

    real(rsp), allocatable :: work(:)
    integer n, info, iitype
    integer :: int_idealwork, iidealwork(1)
    integer, allocatable :: iwork(:)
    real(rsp) :: idealwork(1)
    character(len=100) :: error_string

    @:ASSERT(uplo == 'u' .or. uplo == 'U' .or. uplo == 'l' .or. uplo == 'L')
    @:ASSERT(jobz == 'n' .or. jobz == 'N' .or. jobz == 'v' .or. jobz == 'V')
    @:ASSERT(all(shape(a)==shape(b)))
    @:ASSERT(all(shape(a)==size(w,dim=1)))
    n=size(a,dim=1)
    @:ASSERT(n>0)
    if (present(itype)) then
      iitype = itype
    else
      iitype = 1
    end if
    @:ASSERT(iitype >= 1 .and. iitype <= 3 )
    call ssygvd(iitype, jobz, uplo, n, a, n, b, n, w, idealwork, -1, &
         & iidealwork, -1, info)
    if (info/=0) then
       call error("Failure in SSYGVD to determine optimum workspace")
    endif
    int_idealwork=floor(idealwork(1))
    allocate(work(int_idealwork))
    allocate(iwork(iidealwork(1)))
    call ssygvd(iitype, jobz, uplo, n, a, n, b, n, w, work, int_idealwork, &
         & iwork, iidealwork(1), info)
    if (info/=0) then
       if (info<0) then
99280 format ('Failure in diagonalisation routine ssygvd,', &
          & ' illegal argument at position ',i6)
          write(error_string, 99280) info
          call error(error_string)
       else if (info <= n) then
99290 format ('Failure in diagonalisation routine ssygvd,', &
          & ' diagonal element ',i6,' did not converge to zero.')
          write(error_string, 99290) info
          call error(error_string)
       else
99300 format ('Failure in diagonalisation routine ssygvd,', &
          & ' non-positive definite overlap! Minor ',i6,' responsible.')
          write(error_string, 99300) info - n
          call error(error_string)
       endif
    endif

  end subroutine real_ssygvd


  !> Double precision eigensolver for generalized symmetric matrix problem divide and conquer
  subroutine dble_dsygvd(a,b,w,uplo,jobz,itype)

    !> contains the matrix for the solver, returns eigenvectors if requested (matrix always
    !> overwritten on return anyway)
    real(rdp), intent(inout) :: a(:,:)

    !> contains the second matrix for the solver (overwritten by Cholesky factorization)
    real(rdp), intent(inout) :: b(:,:)

    !> eigenvalues
    real(rdp), intent(out) :: w(:)

    !> upper or lower triangle of the matrix
    character, intent(in) :: uplo

    !> compute eigenvalues 'N' or eigenvalues and eigenvectors 'V'
    character, intent(in) :: jobz

    !> optional specifies the problem type to be solved 1:A*x=(lambda)*B*x, 2:A*B*x=(lambda)*x,
    !> 3:B*A*x=(lambda)*x default is 1
    integer, optional, intent(in) :: itype

    real(rdp), allocatable :: work(:)
    integer n, info, iitype
    integer :: int_idealwork, iidealwork(1)
    integer, allocatable :: iwork(:)
    real(rdp) :: idealwork(1)
    character(len=100) :: error_string

    @:ASSERT(uplo == 'u' .or. uplo == 'U' .or. uplo == 'l' .or. uplo == 'L')
    @:ASSERT(jobz == 'n' .or. jobz == 'N' .or. jobz == 'v' .or. jobz == 'V')
    @:ASSERT(all(shape(a)==shape(b)))
    @:ASSERT(all(shape(a)==size(w,dim=1)))
    n=size(a,dim=1)
    @:ASSERT(n>0)
    if (present(itype)) then
      iitype = itype
    else
      iitype = 1
    end if
    @:ASSERT(iitype >= 1 .and. iitype <= 3 )
    call dsygvd(iitype, jobz, uplo, n, a, n, b, n, w, idealwork, -1, &
        & iidealwork, -1, info)
    if (info/=0) then
      call error("Failure in DSYGVD to determine optimum workspace")
    endif
    int_idealwork=floor(idealwork(1))
    allocate(work(int_idealwork))
    allocate(iwork(iidealwork(1)))
    call dsygvd(iitype, jobz, uplo, n, a, n, b, n, w, work, int_idealwork, &
        & iwork, iidealwork(1), info)
    if (info/=0) then
      if (info<0) then
99310 format ('Failure in diagonalisation routine dsygvd,', &
          & ' illegal argument at position ',i6)
        write(error_string, 99310) info
        call error(error_string)
      else if (info <= n) then
99320 format ('Failure in diagonalisation routine dsygvd,', &
          & ' diagonal element ',i6,' did not converge to zero.')
        write(error_string, 99320) info
        call error(error_string)
      else
99330 format ('Failure in diagonalisation routine dsygvd,', &
          & ' non-positive definite overlap! Minor ',i6,' responsible.')
        write(error_string, 99330) info - n
        call error(error_string)
      endif
    endif

  end subroutine dble_dsygvd


  !> Complex eigensolver for generalized Hermitian matrix problem divide and conquer
  subroutine cmplx_chegvd(a,b,w,uplo,jobz,itype)

    !> contains the matrix for the solver, returns eigenvectors if requested (matrix always
    !> overwritten on return anyway)
    complex(rsp), intent(inout) :: a(:,:)

    !> contains the second matrix for the solver (overwritten by Cholesky factorization)
    complex(rsp), intent(inout) :: b(:,:)

    !> eigenvalues
    real(rsp), intent(out) :: w(:)

    !> upper or lower triangle of the matrix
    character, intent(in) :: uplo

    !> compute eigenvalues 'N' or eigenvalues and eigenvectors 'V'
    character, intent(in) :: jobz

    !> optional specifies the problem type to be solved 1:A*x=(lambda)*B*x, 2:A*B*x=(lambda)*x,
    !> 3:B*A*x=(lambda)*x default is 1
    integer, optional, intent(in) :: itype

    complex(rsp), allocatable :: work(:)
    real(rsp), allocatable :: rwork(:)
    integer n, info, iitype
    integer :: int_idealwork, int_ridealwork, iidealwork(1)
    integer, allocatable :: iwork(:)
    complex(rsp) :: idealwork(1)
    real(rsp) :: ridealwork(1)
    character(len=100) :: error_string

    @:ASSERT(uplo == 'u' .or. uplo == 'U' .or. uplo == 'l' .or. uplo == 'L')
    @:ASSERT(jobz == 'n' .or. jobz == 'N' .or. jobz == 'v' .or. jobz == 'V')
    @:ASSERT(all(shape(a)==shape(b)))
    @:ASSERT(all(shape(a)==size(w,dim=1)))
    n=size(a,dim=1)
    @:ASSERT(n>0)
    if (present(itype)) then
      iitype = itype
    else
      iitype = 1
    end if
    @:ASSERT(iitype >= 1 .and. iitype <= 3 )
    call chegvd(iitype, jobz, uplo, n, a, n, b, n, w, idealwork, -1, &
        & ridealwork, -1, iidealwork, -1, info)
    if (info/=0) then
      call error("Failure in CHEGVD to determine optimum workspace")
    endif
    int_idealwork=floor(real(idealwork(1)))
    int_ridealwork=floor(ridealwork(1))
    allocate(work(int_idealwork))
    allocate(rwork(int_ridealwork))
    allocate(iwork(iidealwork(1)))
    call chegvd(iitype, jobz, uplo, n, a, n, b, n, w, work, int_idealwork,&
        & rwork, int_ridealwork, iwork, iidealwork(1), info)
    if (info/=0) then
      if (info<0) then
99340 format ('Failure in diagonalisation routine chegvd,', &
          & ' illegal argument at position ',i6)
        write(error_string, 99340) info
        call error(error_string)
      else if (info <= n) then
99350 format ('Failure in diagonalisation routine chegvd,', &
          & ' diagonal element ',i6,' did not converge to zero.')
        write(error_string, 99350) info
        call error(error_string)
      else
99360 format ('Failure in diagonalisation routine chegvd,', &
          & ' non-positive definite overlap! Minor ',i6,' responsible.')
        write(error_string, 99360) info - n
        call error(error_string)
      endif
    endif

  end subroutine cmplx_chegvd


  !> Double complex eigensolver for generalized Hermitian matrix problem divide and conquer
  subroutine dblecmplx_zhegvd(a,b,w,uplo,jobz,itype)

    !> contains the matrix for the solver, returns eigenvectors if requested (matrix always
    !> overwritten on return anyway)
    complex(rdp), intent(inout) :: a(:,:)

    !> contains the second matrix for the solver (overwritten by Cholesky factorization)
    complex(rdp), intent(inout) :: b(:,:)

    !> eigenvalues
    real(rdp), intent(out) :: w(:)

    !> upper or lower triangle of the matrix
    character, intent(in) :: uplo

    !> compute eigenvalues 'N' or eigenvalues and eigenvectors 'V'
    character, intent(in) :: jobz

    !> optional specifies the problem type to be solved 1:A*x=(lambda)*B*x, 2:A*B*x=(lambda)*x,
    !> 3:B*A*x=(lambda)*x default is 1
    integer, optional, intent(in) :: itype

    complex(rdp), allocatable :: work(:)
    real(rdp), allocatable :: rwork(:)
    integer n, info, iitype
    integer :: int_idealwork, int_ridealwork, iidealwork(1)
    integer, allocatable :: iwork(:)
    complex(rdp) :: idealwork(1)
    real(rdp) :: ridealwork(1)
    character(len=100) :: error_string

    @:ASSERT(uplo == 'u' .or. uplo == 'U' .or. uplo == 'l' .or. uplo == 'L')
    @:ASSERT(jobz == 'n' .or. jobz == 'N' .or. jobz == 'v' .or. jobz == 'V')
    @:ASSERT(all(shape(a)==shape(b)))
    @:ASSERT(all(shape(a)==size(w,dim=1)))
    n=size(a,dim=1)
    @:ASSERT(n>0)
    if (present(itype)) then
      iitype = itype
    else
      iitype = 1
    end if
    @:ASSERT(iitype >= 1 .and. iitype <= 3 )
    call zhegvd(iitype, jobz, uplo, n, a, n, b, n, w, idealwork, -1, &
        & ridealwork, -1, iidealwork, -1, info)
    if (info/=0) then
      call error("Failure in ZHEGVD to determine optimum workspace")
    endif
    int_idealwork=floor(real(idealwork(1)))
    int_ridealwork=floor(ridealwork(1))
    allocate(work(int_idealwork))
    allocate(rwork(int_ridealwork))
    allocate(iwork(iidealwork(1)))
    call zhegvd(iitype, jobz, uplo, n, a, n, b, n, w, work, int_idealwork, &
        & rwork, int_ridealwork, iwork, iidealwork(1), info)
    if (info/=0) then
      if (info<0) then
99370 format ('Failure in diagonalisation routine zhegvd,', &
          & ' illegal argument at position ',i6)
        write(error_string, 99370) info
        call error(error_string)
      else if (info <= n) then
99380 format ('Failure in diagonalisation routine zhegvd,', &
          & ' diagonal element ',i6,' did not converge to zero.')
        write(error_string, 99380) info
        call error(error_string)
      else
99390 format ('Failure in diagonalisation routine zhegvd,', &
          & ' non-positive definite overlap! Minor ',i6,' responsible.')
        write(error_string, 99390) info - n
        call error(error_string)
      endif
    endif

  end subroutine dblecmplx_zhegvd


  !> Real eigensolver for generalized symmetric matrix problem - Relatively Robust
  !> Representation, optionally use the subspace form if w is smaller than the size of a and b, then
  !> only the first n eigenvalues/eigenvectors are found.
  !> This version re-uses a triangle of a matrix (saving an additional allocation that was in the
  !> previous version).
  !> Based in part on deMon routine from T. Heine
  subroutine real_ssygvr(a,b,w,uplo,jobz,itype,ilIn,iuIn)

    !> contains the matrix for the solver, returns eigenvectors if requested (matrix always
    !> overwritten on return anyway)
    real(rsp), intent(inout) :: a(:,:)

    !> contains the second matrix for the solver (overwritten by Cholesky factorization)
    real(rsp), intent(inout) :: b(:,:)

    !> eigenvalues
    real(rsp), intent(out) :: w(:)

    !> upper or lower triangle of both matrices
    character, intent(in) :: uplo

    !> compute eigenvalues 'N' or eigenvalues and eigenvectors 'V'
    character, intent(in) :: jobz

    !> specifies the problem type to be solved 1:A*x=(lambda)*B*x, 2:A*B*x=(lambda)*x,
    !> 3:B*A*x=(lambda)*x default is 1
    integer, optional, intent(in) :: itype

    !> lower range of eigenstates
    integer, optional, intent(in) :: ilIn

    !> upper range of eigenstates
    integer, optional, intent(in) :: iuIn

    real(rsp), allocatable :: work(:)
    real(rsp) :: tmpWork(1)
    real(rsp), allocatable :: tmpChole(:)
    integer :: lwork
    integer, allocatable :: iwork(:)
    integer :: tmpIWork(1)
    integer :: liwork
    integer, allocatable :: isuppz(:)
    integer :: n, info, iitype
    integer :: m, neig
    real(rsp) :: abstol
    logical :: wantz,upper
    character :: trans
    real(rsp) :: vl, vu
    integer :: il, iu
    logical :: subspace
    integer :: ii, jj
    character :: uplo_new
    character(len=100) :: error_string

    n = size(a, dim=1)

    @:ASSERT(n > 0)
    @:ASSERT(uplo == 'u' .or. uplo == 'U' .or. uplo == 'l' .or. uplo == 'L')
    @:ASSERT(jobz == 'n' .or. jobz == 'N' .or. jobz == 'v' .or. jobz == 'V')
    @:ASSERT(all(shape(a)==shape(b)))
    @:ASSERT(present(ilIn) .eqv. present(iuIn))

    subspace = (size(w) < n)
    if (subspace) then
      if (present(ilIn)) then
        @:ASSERT(ilIn <= iuIn)
        @:ASSERT(ilIn > 0)
        @:ASSERT(n >= iuIn)
        @:ASSERT(size(w,dim=1) == (iuIn - ilIn + 1))
        il = ilIn
        iu = iuIn
      else
        il = 1
        iu = size(w,dim=1)
      end if
    else
      if (present(ilIn)) then
        @:ASSERT(ilIn == 1 .and. iuIn == n)
      end if
      il = 1
      iu = n
    end if

    if (present(itype)) then
      iitype = itype
    else
      iitype = 1
    end if
    @:ASSERT(iitype >= 1 .and. iitype <= 3 )

    allocate(isuppz(2*n))

    allocate(tmpChole(size(a,dim=1)))

    wantz = ( jobz == 'V' .or. jobz == 'v' )
    upper = ( uplo == 'U' .or. uplo == 'u' )
    abstol = slamch( 'Safe minimum' )

    info = 0

    if (subspace) then
      call ssyevr(jobz,'I',uplo,n,a,size(a,dim=1),vl,vu,il,iu,abstol,m,&
          & w,b,size(b,dim=1),isuppz,tmpWork,-1,tmpIWork,-1,info)
    else
      call ssyevr(jobz,'A',uplo,n,a,size(a,dim=1),vl,vu,il,iu,abstol,m,&
          & w,b,size(b,dim=1),isuppz,tmpWork,-1,tmpIWork,-1,info)
    end if

    if (info/=0) then
      call error("Failure in SSYGVR to determine optimum workspace")
    endif
    lwork = floor(tmpWork(1))
    liwork = floor(real(tmpIWork(1)))
    allocate(work(lwork))
    allocate(iwork(liwork))

    ! Form a Cholesky factorization of B.
    call spotrf( uplo, n, b, n, info )
    if( info /= 0 ) then
      info = n + info
99400 format ('Failure in diagonalisation routine SSYGVR,', &
          & ' unable to complete Cholesky factorization of B ',i6)
      write(error_string, 99400) info
      call error(error_string)
    end if
    ! Transform problem to standard eigenvalue problem and solve.
    call ssygst( iitype, uplo, n, a, n, b, n, info )

    if( info /= 0 ) then
      write(error_string, *)'Failure in SSYGVR to transform to standard',info
      call error(error_string)
    end if

    if ( wantz ) then
      ! Save Cholesky factor in the other triangle of H and tmpChole
      do ii = 1, n
        tmpChole(ii) = b(ii,ii)
      end do
      if ( upper ) then
        do jj = 1, n
          do ii = jj+1, n
            a(ii,jj) = b(jj,ii)
          end do
        end do
      else
        do jj = 1, n
          do ii = 1, jj-1
            a(ii,jj) = b(jj,ii)
          end do
        end do
      end if
    end if

    if (subspace) then
      call ssyevr( jobz, 'I', uplo, n, a, n, vl, vu, il, iu, &
          & abstol, m, w, b, n, isuppz, work, lwork, iwork, &
          & liwork, info )
    else
      call ssyevr( jobz, 'A', uplo, n, a, n, vl, vu, il, iu, &
          & abstol, m, w, b, n, isuppz, work, lwork, iwork, &
          & liwork, info )
    end if
    if( info /= 0 ) then
99410 format ('Failure in SSYEVR ',I6)
      write(error_string, 99410) info
      call error(error_string)
    end if

    if ( wantz ) then
      ! Backtransform eigenvectors to the original problem.

      do ii = 1, n
        a(ii,ii) = tmpChole(ii)
      end do

      if ( upper ) then
        uplo_new = 'L'
        upper = .false.
      else
        uplo_new = 'U'
        upper = .true.
      end if

      neig = n
      if( info > 0 ) then
        neig = info - 1
      end if
      if( iitype == 1 .or. iitype == 2 ) then
        ! For A*x=(lambda)*B*x and A*B*x=(lambda)*x;
        ! backtransform eigenvectors: x = inv(L)'*y or inv(U)*y          !'
        if( upper ) then
          trans = 'N'
        else
          trans = 'T'
        end if
        call strsm('Left',uplo_new,trans,'Non-unit',n,neig,1.0,A,n, B,n)
      else if( iitype == 3 ) then
        ! For B*A*x=(lambda)*x;
        ! backtransform eigenvectors: x = L*y or U'*y     !'
        if( upper ) then
          trans = 'T'
        else
          trans = 'N'
        end if
        call strmm('Left',uplo_new,trans,'Non-unit',n,neig,1.0,a,n, b,n)
      end if
      do ii = 1,m
        a( 1:n, ii) = b( 1:n, ii )
      end do
      a( 1:n, m+1:n ) = 0.0
    end if

  end subroutine real_ssygvr


  !> Double precision eigensolver for generalized symmetric matrix problem - Relatively Robust
  !> Representation, optionally use the subspace form if w is smaller than the size of a and b, then
  !> only the first n eigenvalues/eigenvectors are found.
  !> This version re-uses a triangle of a matrix (saving an additional allocation that was in the
  !> previous version).
  !> Based in part on deMon routine from T. Heine
  subroutine dble_dsygvr(a,b,w,uplo,jobz,itype,ilIn,iuIn)

    !> contains the matrix for the solver, returns eigenvectors if requested (matrix always
    !> overwritten on return anyway)
    real(rdp), intent(inout) :: a(:,:)

    !> contains the second matrix for the solver (overwritten by Cholesky factorization)
    real(rdp), intent(inout) :: b(:,:)

    !> eigenvalues
    real(rdp), intent(out) :: w(:)

    !> upper or lower triangle of both matrices
    character, intent(in) :: uplo

    !> compute eigenvalues 'N' or eigenvalues and eigenvectors 'V'
    character, intent(in) :: jobz

    !> specifies the problem type to be solved 1:A*x=(lambda)*B*x, 2:A*B*x=(lambda)*x,
    !> 3:B*A*x=(lambda)*x default is 1
    integer, optional, intent(in) :: itype

    !> lower range of eigenstates
    integer, optional, intent(in) :: ilIn

    !> upper range of eigenstates
    integer, optional, intent(in) :: iuIn

    real(rdp), allocatable :: work(:)
    real(rdp) :: tmpWork(1)
    real(rdp), allocatable :: tmpChole(:)
    integer :: lwork
    integer, allocatable :: iwork(:)
    integer :: tmpIWork(1)
    integer :: liwork
    integer, allocatable :: isuppz(:)
    integer :: n, info, iitype
    integer :: m, neig
    real(rdp) :: abstol
    logical :: wantz,upper
    character :: trans
    real(rdp) :: vl, vu
    integer :: il, iu
    logical :: subspace
    integer :: ii, jj
    character :: uplo_new
    character(len=100) :: error_string

    n = size(a, dim=1)

    @:ASSERT(n > 0)
    @:ASSERT(uplo == 'u' .or. uplo == 'U' .or. uplo == 'l' .or. uplo == 'L')
    @:ASSERT(jobz == 'n' .or. jobz == 'N' .or. jobz == 'v' .or. jobz == 'V')
    @:ASSERT(all(shape(a)==shape(b)))
    @:ASSERT(present(ilIn) .eqv. present(iuIn))

    subspace = (size(w) < n)
    if (subspace) then
      if (present(ilIn)) then
        @:ASSERT(ilIn <= iuIn)
        @:ASSERT(ilIn > 0)
        @:ASSERT(n >= iuIn)
        @:ASSERT(size(w,dim=1) == (iuIn - ilIn + 1))
        il = ilIn
        iu = iuIn
      else
        il = 1
        iu = size(w,dim=1)
      end if
    else
      if (present(ilIn)) then
        @:ASSERT(ilIn == 1 .and. iuIn == n)
      end if
      il = 1
      iu = n
    end if

    if (present(itype)) then
      iitype = itype
    else
      iitype = 1
    end if
    @:ASSERT(iitype >= 1 .and. iitype <= 3 )

    allocate(isuppz(2*n))

    allocate(tmpChole(size(a,dim=1)))

    wantz = ( jobz == 'V' .or. jobz == 'v' )
    upper = ( uplo == 'U' .or. uplo == 'u' )
    abstol = dlamch( 'Safe minimum' )

    info = 0

    if (subspace) then
      call dsyevr(jobz,'I',uplo,n,a,size(a,dim=1),vl,vu,il,iu,abstol,m,&
          & w,b,size(b,dim=1),isuppz,tmpWork,-1,tmpIWork,-1,info)
    else
      call dsyevr(jobz,'A',uplo,n,a,size(a,dim=1),vl,vu,il,iu,abstol,m,&
          & w,b,size(b,dim=1),isuppz,tmpWork,-1,tmpIWork,-1,info)
    end if

    if (info/=0) then
      call error("Failure in DSYGVR to determine optimum workspace")
    endif
    lwork = floor(tmpWork(1))
    liwork = floor(real(tmpIWork(1)))
    allocate(work(lwork))
    allocate(iwork(liwork))

    ! Form a Cholesky factorization of B.
    call dpotrf( uplo, n, b, n, info )
    if( info /= 0 ) then
      info = n + info
99400 format ('Failure in diagonalisation routine DSYGVR,', &
          & ' unable to complete Cholesky factorization of B ',i6)
      write(error_string, 99400) info
      call error(error_string)
    end if
    ! Transform problem to standard eigenvalue problem and solve.
    call dsygst( iitype, uplo, n, a, n, b, n, info )

    if( info /= 0 ) then
      write(error_string, *)'Failure in DSYGVR to transform to standard',info
      call error(error_string)
    end if

    if ( wantz ) then
      ! Save Cholesky factor in the other triangle of H and tmpChole
      do ii = 1, n
        tmpChole(ii) = b(ii,ii)
      end do
      if ( upper ) then
        do jj = 1, n
          do ii = jj + 1, n
            a(ii,jj) = b(jj,ii)
          end do
        end do
      else
        do jj = 1, n
          do ii = 1, jj - 1
            a( ii, jj ) = b( jj, ii )
          end do
        end do
      end if
    end if

    if (subspace) then
      call dsyevr( jobz, 'I', uplo, n, a, n, vl, vu, il, iu, &
          & abstol, m, w, b, n, isuppz, work, lwork, iwork, &
          & liwork, info )
    else
      call dsyevr( jobz, 'A', uplo, n, a, n, vl, vu, il, iu, &
          & abstol, m, w, b, n, isuppz, work, lwork, iwork, &
          & liwork, info )
    end if
    if( info /= 0 ) then
99410 format ('Failure in DSYEVR ',I6)
      write(error_string, 99410) info
      call error(error_string)
    end if

    if ( wantz ) then
      ! Backtransform eigenvectors to the original problem.

      do ii = 1, n
        a(ii,ii) = tmpChole(ii)
      end do

      if ( upper ) then
        uplo_new = 'L'
        upper = .false.
      else
        uplo_new = 'U'
        upper = .true.
      end if

      neig = n
      if( info > 0 ) then
        neig = info - 1
      end if
      if( iitype == 1 .or. iitype == 2 ) then
        ! For A*x=(lambda)*B*x and A*B*x=(lambda)*x;
        ! backtransform eigenvectors: x = inv(L)'*y or inv(U)*y          !'
        if( upper ) then
          trans = 'N'
        else
          trans = 'T'
        end if
        call dtrsm('Left',uplo_new,trans,'Non-unit',n,neig,1.0d0,A,n, B,n)
      else if( iitype == 3 ) then
        ! For B*A*x=(lambda)*x;
        ! backtransform eigenvectors: x = L*y or U'*y               !'
        if( upper ) then
          trans = 'T'
        else
          trans = 'N'
        end if
        call dtrmm('Left',uplo_new,trans,'Non-unit',n,neig,1.0d0,a,n, b,n)
      end if
      do ii = 1,m
        a( 1:n, ii) = b( 1:n, ii )
      end do
      a( 1:n, m+1:n ) = 0.0d0
    end if

  end subroutine dble_dsygvr


  !> Complex precision eigensolver for generalized symmetric matrix problem - Relatively Robust
  !> Representation, optionally use the subspace form if w is smaller than the size of a and b, then
  !> only the first n eigenvalues/eigenvectors are found.
  !> This version re-uses a triangle of a matrix (saving an additional allocation that was in the
  !> previous version).
  !> Based in part on deMon routine from T. Heine
  subroutine cmplx_chegvr(a,b,w,uplo,jobz,itype,ilIn,iuIn)

    !> contains the matrix for the solver, returns eigenvectors if requested (matrix always
    !> overwritten on return anyway)
    complex(rsp), intent(inout) :: a(:,:)

    !> contains the second matrix for the solver (overwritten by Cholesky factorization)
    complex(rsp), intent(inout) :: b(:,:)

    !> eigenvalues
    real(rsp), intent(out) :: w(:)

    !> upper or lower triangle of both matrices
    character, intent(in) :: uplo

    !> compute eigenvalues 'N' or eigenvalues and eigenvectors 'V'
    character, intent(in) :: jobz

    !> specifies the problem type to be solved 1:A*x=(lambda)*B*x, 2:A*B*x=(lambda)*x,
    !> 3:B*A*x=(lambda)*x default is 1
    integer, optional, intent(in) :: itype

    !> lower range of eigenstates
    integer, optional, intent(in) :: ilIn

    !> upper range of eigenstates
    integer, optional, intent(in) :: iuIn

    complex(rsp), allocatable :: work(:)
    complex(rsp) :: tmpWork(1)
    real(rsp), allocatable :: rWork(:)
    real(rsp) :: tmpRWork(1)
    complex(rsp), allocatable :: tmpChole(:)
    integer :: lwork
    integer, allocatable :: iwork(:)
    integer :: tmpIWork(1)
    integer :: liwork
    integer :: lrwork
    integer, allocatable :: isuppz(:)
    integer :: n, info, iitype
    integer :: m, neig
    real(rsp) :: abstol
    logical :: wantz,upper
    character :: trans
    real(rsp) :: vl, vu
    integer :: il, iu
    logical :: subspace
    integer :: ii, jj
    character :: uplo_new
    character(len=100) :: error_string

    n = size(a, dim=1)
    @:ASSERT(uplo == 'u' .or. uplo == 'U' .or. uplo == 'l' .or. uplo == 'L')
    @:ASSERT(jobz == 'n' .or. jobz == 'N' .or. jobz == 'v' .or. jobz == 'V')
    @:ASSERT(all(shape(a)==shape(b)))

    @:ASSERT(present(ilIn) .eqv. present(iuIn))

    subspace = (size(w) < n)
    if (subspace) then
      if (present(ilIn)) then
        @:ASSERT(ilIn <= iuIn)
        @:ASSERT(ilIn > 0)
        @:ASSERT(n >= iuIn)
        @:ASSERT(size(w,dim=1) == (iuIn - ilIn + 1))
        il = ilIn
        iu = iuIn
      else
        il = 1
        iu = size(w,dim=1)
      end if
    else
      if (present(ilIn)) then
        @:ASSERT(ilIn == 1 .and. iuIn == n)
      end if
      il = 1
      iu = n
    end if

    if (present(itype)) then
      iitype = itype
    else
      iitype = 1
    end if
    @:ASSERT(iitype >= 1 .and. iitype <= 3 )

    allocate(isuppz(2*n))

    allocate(tmpChole(size(a,dim=1)))

    wantz = ( jobz == 'V' .or. jobz == 'v' )
    upper = ( uplo == 'U' .or. uplo == 'u' )
    abstol = SLAMCH( 'Safe minimum' )

    info = 0

    if (subspace) then
      call cheevr(jobz,'I',uplo,n,a,size(a,dim=1),vl,vu,il,iu,abstol,m,&
          & w,b,size(b,dim=1),isuppz,tmpWork,-1,tmpRwork,-1,tmpIWork,-1,info)
    else
      call cheevr(jobz,'A',uplo,n,a,size(a,dim=1),vl,vu,il,iu,abstol,m,&
          & w,b,size(b,dim=1),isuppz,tmpWork,-1,tmpRwork,-1,tmpIWork,-1,info)
    end if

    if (info/=0) then
      call error("Failure in CHEEVR to determine optimum workspace")
    endif
    lwork = floor(real(tmpWork(1)))
    liwork = floor(real(tmpIWork(1)))
    lrwork = floor(tmpRwork(1))
    allocate(work(lwork))
    allocate(iwork(liwork))
    allocate(rwork(lrwork))

    ! Form a Cholesky factorization of B.
    call cpotrf( uplo, n, b, n, info )
    if( info /= 0 ) then
      info = n + info
99400 format ('Failure in diagonalisation routine CHEGVR,', &
          & ' unable to complete Cholesky factorization of B ',i6)
      write(error_string, 99400) info
      call error(error_string)
    end if
    ! Transform problem to standard eigenvalue problem and solve.
    call chegst( iitype, uplo, n, a, n, b, n, info )

    if( info /= 0 ) then
      write(error_string, *)'Failure in CHEGST to transform to standard',info
      call error(error_string)
    end if

    if ( wantz ) then
      ! Save Cholesky factor in the other triangle of H and tmpChole
      do ii = 1, n
        tmpChole(ii) = b(ii,ii)
      end do
      if ( upper ) then
        do jj = 1, n
          do ii = jj+1, n
            a(ii,jj) = conjg(b(jj,ii))
          end do
        end do
      else
        do jj = 1, n
          do ii = 1, jj-1
            a(ii,jj) = conjg(b(jj,ii))
          end do
        end do
      end if
    end if

    if (subspace) then
      call cheevr( jobz, 'I', uplo, n, a, n, vl, vu, il, iu, &
          & abstol, m, w, b, n, isuppz, work, lwork, rwork, lrwork, iwork, &
          & liwork, info )
    else
      call cheevr( jobz, 'A', uplo, n, a, n, vl, vu, il, iu, &
          & abstol, m, w, b, n, isuppz, work, lwork, rwork, lrwork, iwork, &
          & liwork, info )
    end if
    if( info /= 0 ) then
99410 format ('Failure in CHEEVR ',I6)
      write(error_string, 99410) info
      call error(error_string)
    end if

    if ( wantz ) then
      ! Backtransform eigenvectors to the original problem.

      do ii = 1, n
        a(ii,ii) = tmpChole(ii)
      end do

      if ( upper ) then
        uplo_new = 'L'
        upper = .false.
      else
        uplo_new = 'U'
        upper = .true.
      end if

      neig = n
      if( info > 0 ) then
        neig = info - 1
      end if
      if( iitype == 1 .or. iitype == 2 ) then
        ! For A*x=(lambda)*B*x and A*B*x=(lambda)*x;
        ! backtransform eigenvectors: x = inv(L)'*y or inv(U)*y          !'
        if( upper ) then
          trans = 'N'
        else
          trans = 'C'
        end if
        call ctrsm('Left',uplo_new,trans,'Non-unit',n,neig, &
            & cmplx(1.0,0.0,rsp),A,n, B,n)
      else if( iitype == 3 ) then
        ! For B*A*x=(lambda)*x;
        ! backtransform eigenvectors: x = L*y or U'*y               !'
        if( upper ) then
          trans = 'C'
        else
          trans = 'N'
        end if
        call ctrmm('Left',uplo_new,trans,'Non-unit',n,neig, &
            & cmplx(1.0,0.0,rsp),a,n, b,n)
      end if
      do ii = 1,m
        a( 1:n, ii) = b( 1:n, ii )
      end do
      a( 1:n, m+1:n ) = cmplx(0.0,0.0,rsp)
    end if

  end subroutine cmplx_chegvr


  !> Double complex precision eigensolver for generalized symmetric matrix problem - Relatively
  !> Robust Representation, optionally use the subspace form if w is smaller than the size of a and
  !> b, then only the first n eigenvalues/eigenvectors are found.
  !> This version re-uses a triangle of a matrix (saving an additional allocation that was in the
  !> previous version).
  !> Based in part on deMon routine from T. Heine
  subroutine dblecmplx_zhegvr(a,b,w,uplo,jobz,itype,ilIn,iuIn)

    !> contains the matrix for the solver, returns eigenvectors if requested (matrix always
    !> overwritten on return anyway)
    complex(rdp), intent(inout) :: a(:,:)

    !> contains the second matrix for the solver (overwritten by Cholesky factorization)
    complex(rdp), intent(inout) :: b(:,:)

    !> eigenvalues
    real(rdp), intent(out) :: w(:)

    !> upper or lower triangle of both matrices
    character, intent(in) :: uplo

    !> compute eigenvalues 'N' or eigenvalues and eigenvectors 'V'
    character, intent(in) :: jobz

    !> specifies the problem type to be solved 1:A*x=(lambda)*B*x, 2:A*B*x=(lambda)*x,
    !> 3:B*A*x=(lambda)*x default is 1
    integer, optional, intent(in) :: itype

    !> lower range of eigenstates
    integer, optional, intent(in) :: ilIn

    !> upper range of eigenstates
    integer, optional, intent(in) :: iuIn

    complex(rdp), allocatable :: work(:)
    complex(rdp) :: tmpWork(1)
    real(rdp), allocatable :: rWork(:)
    real(rdp) :: tmpRWork(1)
    complex(rdp), allocatable :: tmpChole(:)
    integer :: lwork
    integer, allocatable :: iwork(:)
    integer :: tmpIWork(1)
    integer :: liwork
    integer :: lrwork
    integer, allocatable :: isuppz(:)
    integer :: n, info, iitype
    integer :: m, neig
    real(rdp) :: abstol
    logical :: wantz,upper
    character :: trans
    real(rdp) :: vl, vu
    integer :: il, iu
    logical :: subspace
    integer :: ii, jj
    character :: uplo_new
    character(len=100) :: error_string

    n = size(a, dim=1)

    @:ASSERT(uplo == 'u' .or. uplo == 'U' .or. uplo == 'l' .or. uplo == 'L')
    @:ASSERT(jobz == 'n' .or. jobz == 'N' .or. jobz == 'v' .or. jobz == 'V')
    @:ASSERT(all(shape(a)==shape(b)))
    @:ASSERT(present(ilIn) .eqv. present(iuIn))

    subspace = (size(w) < n)
    if (subspace) then
      if (present(ilIn)) then
        @:ASSERT(ilIn <= iuIn)
        @:ASSERT(ilIn > 0)
        @:ASSERT(n >= iuIn)
        @:ASSERT(size(w,dim=1) == (iuIn - ilIn + 1))
        il = ilIn
        iu = iuIn
      else
        il = 1
        iu = size(w,dim=1)
      end if
    else
      if (present(ilIn)) then
        @:ASSERT(ilIn == 1 .and. iuIn == n)
      end if
      il = 1
      iu = n
    end if

    if (present(itype)) then
      iitype = itype
    else
      iitype = 1
    end if
    @:ASSERT(iitype >= 1 .and. iitype <= 3 )

    allocate(isuppz(2*n))

    allocate(tmpChole(size(a,dim=1)))

    wantz = ( jobz == 'V' .or. jobz == 'v' )
    upper = ( uplo == 'U' .or. uplo == 'u' )
    abstol = DLAMCH( 'Safe minimum' )

    info = 0

    if (subspace) then
      call zheevr(jobz,'I',uplo,n,a,size(a,dim=1),vl,vu,il,iu,abstol,m,&
          & w,b,size(b,dim=1),isuppz,tmpWork,-1,tmpRwork,-1,tmpIWork,-1,info)
    else
      call zheevr(jobz,'A',uplo,n,a,size(a,dim=1),vl,vu,il,iu,abstol,m,&
          & w,b,size(b,dim=1),isuppz,tmpWork,-1,tmpRwork,-1,tmpIWork,-1,info)
    end if

    if (info/=0) then
      call error("Failure in ZHEEVR to determine optimum workspace")
    endif
    lwork = floor(real(tmpWork(1)))
    liwork = floor(real(tmpIWork(1)))
    lrwork = floor(tmpRwork(1))
    allocate(work(lwork))
    allocate(iwork(liwork))
    allocate(rwork(lrwork))

    ! Form a Cholesky factorization of B.
    call zpotrf( uplo, n, b, n, info )
    if( info /= 0 ) then
      info = n + info
99400 format ('Failure in diagonalisation routine ZHEEVR,', &
          & ' unable to complete Cholesky factorization of B ',i6)
      write(error_string, 99400) info
      call error(error_string)
    end if
    ! Transform problem to standard eigenvalue problem and solve.
    call zhegst( iitype, uplo, n, a, n, b, n, info )

    if( info /= 0 ) then
      write(error_string, *)'Failure in ZHEGST to transform to standard',info
      call error(error_string)
    end if

    if ( wantz ) then
      ! Save Cholesky factor in the other triangle of H and tmpChole
      do ii = 1, n
        tmpChole(ii) = b(ii,ii)
      end do
      if ( upper ) then
        do jj = 1, n
          do ii = jj+1, n
            a(ii,jj) = conjg(b(jj,ii))
          end do
        end do
      else
        do jj = 1, n
          do ii = 1, jj-1
            a(ii,jj) = conjg(b(jj,ii))
          end do
        end do
      end if
    end if

    if (subspace) then
      call zheevr( jobz, 'I', uplo, n, a, n, vl, vu, il, iu, &
          & abstol, m, w, b, n, isuppz, work, lwork, rwork, lrwork, iwork, &
          & liwork, info )
    else
      call zheevr( jobz, 'A', uplo, n, a, n, vl, vu, il, iu, &
          & abstol, m, w, b, n, isuppz, work, lwork, rwork, lrwork, iwork, &
          & liwork, info )
    end if
    if( info /= 0 ) then
99410 format ('Failure in ZHEEVR ',I6)
      write(error_string, 99410) info
      call error(error_string)
    end if

    if ( wantz ) then
      ! Backtransform eigenvectors to the original problem.

      do ii = 1, n
        a(ii,ii) = tmpChole(ii)
      end do

      if ( upper ) then
        uplo_new = 'L'
        upper = .false.
      else
        uplo_new = 'U'
        upper = .true.
      end if

      neig = n
      if( info > 0 ) then
        neig = info - 1
      end if
      if( iitype == 1 .or. iitype == 2 ) then
        ! For A*x=(lambda)*B*x and A*B*x=(lambda)*x;
        ! backtransform eigenvectors: x = inv(L)'*y or inv(U)*y          !'
        if( upper ) then
          trans = 'N'
        else
          trans = 'C'
        end if
        call ztrsm('Left',uplo_new,trans,'Non-unit',n,neig, &
            & cmplx(1.0,0.0,rdp),A,n, B,n)
      else if( iitype == 3 ) then
        ! For B*A*x=(lambda)*x;
        ! backtransform eigenvectors: x = L*y or U'*y               !'
        if( upper ) then
          trans = 'C'
        else
          trans = 'N'
        end if
        call ztrmm('Left',uplo_new,trans,'Non-unit',n,neig, &
            & cmplx(1.0,0.0,rdp),a,n, b,n)
      end if
      do ii = 1,m
        a( 1:n, ii) = b( 1:n, ii )
      end do
      a( 1:n, m+1:n ) = cmplx(0.0,0.0,rdp)
    end if

  end subroutine dblecmplx_zhegvr


  !> Single precision banded symmetric generalised matrix eigensolver
  subroutine real_ssbgv(ab, bb, w, uplo, z)

    !> contains the matrix for the solver (overwritten before exit)
    real(rsp), intent(inout) :: ab(:,:)

    !> contains the second matrix for the solver (overwritten by split Cholesky factorization)
    real(rsp), intent(inout) :: bb(:,:)

    !> eigenvalues
    real(rsp), intent(out) :: w(:)

    !> upper or lower triangle of both matrices
    character, intent(in) :: uplo

    !> returns calculated eigenvectors if present
    real(rsp), optional, intent(out) :: z(:,:)

    real(rsp), allocatable :: work(:)
    integer :: n, ka, kb, ldab, ldbb, ldz, info
    character :: jobz
    real(rsp) :: zTmp(1,1)
    character(len=100) :: error_string

    if (present(z)) then
      jobz = 'v'
    else
      jobz = 'n'
    end if

    @:ASSERT(uplo == 'u' .or. uplo == 'U' .or. uplo == 'l' .or. uplo == 'L')
    @:ASSERT(jobz == 'n' .or. jobz == 'N' .or. jobz == 'v' .or. jobz == 'V')
    @:ASSERT(all(shape(ab)==shape(bb)))

    n = size(ab,dim=2)

    ldab = size(ab,dim=1)
    ldbb = size(bb,dim=1)
    ka = ldab - 1
    kb = ldbb - 1
    @:ASSERT(ka >= 0)
    @:ASSERT(kb >= 0)

    if (present(z)) then
      ldz = n
    else
      ldz = 1
    end if

    allocate(work(3*n))

    if (present(z)) then
      call ssbgv( jobz,uplo,n,ka,kb,ab,ldab,bb,ldbb,w,z,ldz,work,info )
    else
      call ssbgv( jobz,uplo,n,ka,kb,ab,ldab,bb,ldbb,w,zTmp,ldz,work,info )
    end if

    if (info/=0) then
      if (info<0) then
99440 format ('Failure in diagonalisation routine ssbgv,', &
          & ' illegal argument at position ',i6)
          write(error_string, 99440) info
          call error(error_string)
       else if (info <= n) then
99450 format ('Failure in diagonalisation routine ssbgv,', &
          & ' tri diagonal element ',i6,' did not converge to zero.')
          write(error_string, 99450) info
          call error(error_string)
       else
99460 format ('Failure in diagonalisation routine ssbgv,', &
          & ' non-positive definite overlap! ',i6)
          write(error_string, 99460) info
          call error(error_string)
       endif
    endif

  end subroutine real_ssbgv


  !> Double precision banded symmetric generalised matrix eigensolver
  subroutine dble_dsbgv(ab, bb, w, uplo, z)

    !> contains the matrix for the solver (overwritten before exit)
    real(rdp), intent(inout) :: ab(:,:)

    !> contains the second matrix for the solver (overwritten by split Cholesky factorization)
    real(rdp), intent(inout) :: bb(:,:)

    !> eigenvalues
    real(rdp), intent(out) :: w(:)

    !> upper or lower triangle of both matrices
    character, intent(in) :: uplo

    !> returns calculated eigenvectors if present
    real(rdp), optional, intent(out) :: z(:,:)

    real(rdp), allocatable :: work(:)
    integer :: n, ka, kb, ldab, ldbb, ldz, info
    character :: jobz
    real(rdp) :: zTmp(1,1)
    character(len=100) :: error_string

    if (present(z)) then
      jobz = 'v'
    else
      jobz = 'n'
    end if

    @:ASSERT(uplo == 'u' .or. uplo == 'U' .or. uplo == 'l' .or. uplo == 'L')
    @:ASSERT(jobz == 'n' .or. jobz == 'N' .or. jobz == 'v' .or. jobz == 'V')
    @:ASSERT(all(shape(ab)==shape(bb)))

    n = size(ab,dim=2)

    ldab = size(ab,dim=1)
    ldbb = size(bb,dim=1)
    ka = ldab - 1
    kb = ldbb - 1
    @:ASSERT(ka >= 0)
    @:ASSERT(kb >= 0)

    if (present(z)) then
      ldz = n
    else
      ldz = 1
    end if

    allocate(work(3*n))

    if (present(z)) then
      call dsbgv( jobz,uplo,n,ka,kb,ab,ldab,bb,ldbb,w,z,ldz,work,info )
    else
      call dsbgv( jobz,uplo,n,ka,kb,ab,ldab,bb,ldbb,w,zTmp,ldz,work,info )
    end if

    if (info/=0) then
      if (info<0) then
99470 format ('Failure in diagonalisation routine dsbgv,', &
          & ' illegal argument at position ',i6)
          write(error_string, 99470) info
          call error(error_string)
       else if (info <= n) then
99480 format ('Failure in diagonalisation routine dsbgv,', &
          & ' tri diagonal element ',i6,' did not converge to zero.')
          write(error_string, 99480) info
          call error(error_string)
       else
99490 format ('Failure in diagonalisation routine dsbgv,', &
          & ' non-positive definite overlap! ',i6)
          write(error_string, 99490) info
          call error(error_string)
       endif
    endif

  end subroutine dble_dsbgv


  !> Complex banded symmetric generalised matrix eigensolver
  subroutine cmplx_chbgv(ab, bb, w, uplo, z)

    !> contains the matrix for the solver (overwritten before exit)
    complex(rsp), intent(inout) :: ab(:,:)

    !> contains the second matrix for the solver (overwritten by split Cholesky factorization)
    complex(rsp), intent(inout) :: bb(:,:)

    !> eigenvalues
    real(rsp), intent(out) :: w(:)

    !> upper or lower triangle of both matrices
    character, intent(in) :: uplo

    !> returns calculated eigenvectors if present
    complex(rsp), optional, intent(out) :: z(:,:)

    complex(rsp), allocatable :: work(:)
    real(rsp), allocatable :: rwork(:)
    integer :: n, ka, kb, ldab, ldbb, ldz, info
    character :: jobz
    complex(rsp) :: zTmp(1,1)
    character(len=100) :: error_string

    if (present(z)) then
      jobz = 'v'
    else
      jobz = 'n'
    end if

    @:ASSERT(uplo == 'u' .or. uplo == 'U' .or. uplo == 'l' .or. uplo == 'L')
    @:ASSERT(jobz == 'n' .or. jobz == 'N' .or. jobz == 'v' .or. jobz == 'V')
    @:ASSERT(all(shape(ab)==shape(bb)))

    n = size(ab,dim=2)

    ldab = size(ab,dim=1)
    ldbb = size(bb,dim=1)
    ka = ldab - 1
    kb = ldbb - 1
    @:ASSERT(ka >= 0)
    @:ASSERT(kb >= 0)

    if (present(z)) then
      ldz = n
    else
      ldz = 1
    end if

    allocate(work(n))
    allocate(rwork(3*n))

    if (present(z)) then
      call chbgv( jobz,uplo,n,ka,kb,ab,ldab,bb,ldbb,w,z,ldz,work,rwork,info )
    else
      call chbgv( jobz,uplo,n,ka,kb,ab,ldab,bb,ldbb,w,zTmp,ldz,work,rwork, &
          & info )
    end if

    if (info/=0) then
      if (info<0) then
99500 format ('Failure in diagonalisation routine chbgv,', &
          & ' illegal argument at position ',i6)
          write(error_string, 99500) info
          call error(error_string)
       else if (info <= n) then
99510 format ('Failure in diagonalisation routine chbgv,', &
          & ' tri diagonal element ',i6,' did not converge to zero.')
          write(error_string, 99510) info
          call error(error_string)
       else
99520 format ('Failure in diagonalisation routine chbgv,', &
          & ' non-positive definite overlap! ',i6)
          write(error_string, 99520) info
          call error(error_string)
       endif
    endif

  end subroutine cmplx_chbgv


  !> Complex double precision banded symmetric generalised matrix eigensolver
  subroutine dblecmplx_zhbgv(ab, bb, w, uplo, z)

    !> contains the matrix for the solver (overwritten before exit)
    complex(rdp), intent(inout) :: ab(:,:)

    !> contains the second matrix for the solver (overwritten by split Cholesky factorization)
    complex(rdp), intent(inout) :: bb(:,:)

    !> eigenvalues
    real(rdp), intent(out) :: w(:)

    !> upper or lower triangle of both matrices
    character, intent(in) :: uplo

    !> returns calculated eigenvectors if present
    complex(rdp), optional, intent(out) :: z(:,:)

    complex(rdp), allocatable :: work(:)
    real(rdp), allocatable :: rwork(:)
    integer :: n, ka, kb, ldab, ldbb, ldz, info
    character :: jobz
    complex(rdp) :: zTmp(1,1)
    character(len=100) :: error_string

    if (present(z)) then
      jobz = 'v'
    else
      jobz = 'n'
    end if

    @:ASSERT(uplo == 'u' .or. uplo == 'U' .or. uplo == 'l' .or. uplo == 'L')
    @:ASSERT(jobz == 'n' .or. jobz == 'N' .or. jobz == 'v' .or. jobz == 'V')
    @:ASSERT(all(shape(ab)==shape(bb)))

    n = size(ab,dim=2)

    ldab = size(ab,dim=1)
    ldbb = size(bb,dim=1)
    ka = ldab - 1
    kb = ldbb - 1
    @:ASSERT(ka >= 0)
    @:ASSERT(kb >= 0)

    if (present(z)) then
      ldz = n
    else
      ldz = 1
    end if

    allocate(work(n))
    allocate(rwork(3*n))

    if (present(z)) then
      call zhbgv( jobz,uplo,n,ka,kb,ab,ldab,bb,ldbb,w,z,ldz,work,rwork,info )
    else
      call zhbgv( jobz,uplo,n,ka,kb,ab,ldab,bb,ldbb,w,zTmp,ldz,work,rwork, &
          & info )
    end if

    if (info/=0) then
      if (info<0) then
99530 format ('Failure in diagonalisation routine zhbgv,', &
          & ' illegal argument at position ',i6)
          write(error_string, 99530) info
          call error(error_string)
       else if (info <= n) then
99540 format ('Failure in diagonalisation routine zhbgv,', &
          & ' tri diagonal element ',i6,' did not converge to zero.')
          write(error_string, 99540) info
          call error(error_string)
       else
99550 format ('Failure in diagonalisation routine zhbgv,', &
          & ' non-positive definite overlap! ',i6)
          write(error_string, 99550) info
          call error(error_string)
       endif
    endif

  end subroutine dblecmplx_zhbgv


#:if WITH_GPU

#:for DTYPE, VPREC, VTYPE, NAME in [('real', 's', 'real', 'ssygvd'),&
  & ('dble', 'd', 'real', 'dsygvd'), ('cmplx', 's', 'complex', 'chegvd'),&
  & ('dblecmplx', 'd', 'complex', 'zhegvd')]
  !> Generalised eigensolution for symmetric/hermitian matrices on GPU(s)
  subroutine ${DTYPE}$_magma_${NAME}$(ngpus, a, b, w, uplo, jobz, itype)

    !> Number of GPUs to use
    integer, intent(in) :: ngpus

    !> contains the matrix for the solver, returns eigenvectors if requested (matrix always
    !> overwritten on return anyway)
    ${VTYPE}$(r${VPREC}$p), intent(inout) :: a(:,:)

    !> contains the second matrix for the solver (overwritten by Cholesky factorization)
    ${VTYPE}$(r${VPREC}$p), intent(inout) :: b(:,:)

    !> eigenvalues
    real(r${VPREC}$p), intent(out) :: w(:)

    !> upper or lower triangle of the matrix
    character, intent(in) :: uplo

    !> compute eigenvalues 'N' or eigenvalues and eigenvectors 'V'
    character, intent(in) :: jobz

    !> optional specifies the problem type to be solved 1:A*x=(lambda)*B*x, 2:A*B*x=(lambda)*x,
    !> 3:B*A*x=(lambda)*x default is 1
    integer, optional, intent(in) :: itype

    ${VTYPE}$(r${VPREC}$p), allocatable :: work(:)

  #:if VTYPE == 'complex'
    real(r${VPREC}$p), allocatable :: rwork(:)
    integer :: lrwork
  #:endif

    integer, allocatable :: iwork(:)
    integer :: lwork, liwork, n, info, iitype
    character(len=100) :: error_string

    @:ASSERT(uplo == 'u' .or. uplo == 'U' .or. uplo == 'l' .or. uplo == 'L')
    @:ASSERT(jobz == 'n' .or. jobz == 'N' .or. jobz == 'v' .or. jobz == 'V')
    @:ASSERT(all(shape(a)==shape(b)))
    @:ASSERT(all(shape(a)==size(w,dim=1)))
    n=size(a,dim=1)
    @:ASSERT(n>0)
    if (present(itype)) then
      iitype = itype
    else
      iitype = 1
    end if
    @:ASSERT(iitype >= 1 .and. iitype <= 3 )

    ! Workspace query
    allocate(work(1))
    allocate(iwork(1))
  #:if VTYPE == 'complex'
    allocate(rwork(1))
  #:endif

    call magmaf_${NAME}$_m(ngpus, iitype, jobz, uplo, n, a, n, b, n, w, work, -1,&
      #:if VTYPE == 'complex'
        & rwork, -1,&
      #:endif
        & iwork, -1, info)

    if (info /= 0) then
      call error("Failure in MAGMA_${NAME}$ to determine optimum workspace")
    endif

 #:if VTYPE == 'complex'
    lwork = floor(real(work(1)))
    lrwork = floor(rwork(1))
    liwork = int(iwork(1))
    deallocate(work) ;  allocate(work(lwork))
    deallocate(rwork) ;  allocate(rwork(lrwork))
    deallocate(iwork) ; allocate(iwork(liwork))
  #:endif
  #:if VTYPE == 'real'
    lwork = floor(work(1))
    liwork = int(iwork(1))
    deallocate(work) ;  allocate(work(lwork))
    deallocate(iwork) ; allocate(iwork(liwork))
   #:endif

    ! MAGMA Diagonalization
    call magmaf_${NAME}$_m(ngpus, iitype, jobz, uplo, n, a, n, b, n, w, work, lwork,&
      #:if VTYPE == 'complex'
        & rwork, lrwork,&
      #:endif
        & iwork, liwork, info)

    ! test for errors
    if (info /= 0) then
      if (info < 0) then
        write(error_string, "('Failure in diagonalisation routine magmaf_${NAME}$_m,&
            & illegal argument at position ',i6)") info
        call error(error_string)
      else if (info <= n) then
        write(error_string, "('Failure in diagonalisation routine magmaf_${NAME}$_m,&
            & diagonal element ',i6,' did not converge to zero.')") info
        call error(error_string)
      else
        write(error_string, "('Failure in diagonalisation routine magmaf_${NAME}$_m,&
            & non-positive definite overlap! ',i6)") info - n
        call error(error_string)
      endif
    endif

  end subroutine ${DTYPE}$_magma_${NAME}$

#:endfor

#:endif


#:for DTYPE, VPREC, VTYPE, NAME in [('real', 's', 'real', 'sgeev'), ('dble', 'd', 'real', 'dgeev')]

  !> Simple general matrix eigensolver
  subroutine ${DTYPE}$_${NAME}$(a, wr, wi, vl, vr, err)

    !> Matrix, overwritten on exit
    real(r${VPREC}$p), intent(inout) :: a(:,:)

    !> Real part of eigenvalues
    real(r${VPREC}$p), intent(out) :: wr(:)

    !> Imaginary part of eigenvalues
    real(r${VPREC}$p), intent(out) :: wi(:)

    !> Left eigenvectors
    real(r${VPREC}$p), intent(out), optional :: vl(:,:)

    !> Right eigenvectors
    real(r${VPREC}$p), intent(out), optional :: vr(:,:)

    !> Error code return, 0 if no problems
    integer, intent(out), optional :: err

    real(r${VPREC}$p), allocatable :: work(:)
    integer :: n, lda, info, int_idealwork, ldvl, ldvr
    real(r${VPREC}$p) :: idealwork(1)
    character :: jobvl, jobvr
    character(len=100) :: error_string

    ! If no eigenvectors requested, need a dummy array for lapack call
    real(r${VPREC}$p) :: dummyvl(1,1), dummyvr(1,1)

    if (present(err)) then
      err = 0
    end if

    lda = size(a, dim=1)
    n = size(a, dim=2)

    @:ASSERT(n>0)
    @:ASSERT(size(wr) >= n)
    @:ASSERT(size(wi) >= n)

    if (present(vl)) then
      jobvl = 'V'
      ldvl = size(vl, dim=1)
      @:ASSERT(all(shape(vl)>=[n,n]))
    else
      jobvl = 'N'
      ldvl = 1
    end if
    if (present(vr)) then
      jobvr = 'V'
      ldvr = size(vr, dim=1)
      @:ASSERT(all(shape(vr)>=[n,n]))
    else
      jobvr = 'N'
      ldvr = 1
    end if

    if (jobvl == 'V' .and. jobvr == 'V') then
      call ${VPREC}$geev(jobvl, jobvr, n, a, lda, wr, wi, vl, ldvl, vr, ldvr, idealwork, -1, info)
    else if (jobvl == 'V' .and. jobvr == 'N') then
      call ${VPREC}$geev(jobvl, jobvr, n, a, lda, wr, wi, vl, ldvl, dummyvr, ldvr, idealwork, -1,&
          & info)
    else if (jobvl == 'N' .and. jobvr == 'V') then
      call ${VPREC}$geev(jobvl, jobvr, n, a, lda, wr, wi, dummyvl, ldvl, vr, ldvr, idealwork, -1,&
          & info)
    else if (jobvl == 'N' .and. jobvr == 'N') then
      call ${VPREC}$geev(jobvl, jobvr, n, a, lda, wr, wi, dummyvl, ldvl, dummyvr, ldvr, idealwork,&
          & -1, info)
    end if
    if (info/=0) then
      @:ERROR_HANDLING(err, -1, "Failue in ${VPREC}$geev to determine optimum workspace")
    endif
    int_idealwork=nint(idealwork(1))
    allocate(work(int_idealwork))

    if (jobvl == 'V' .and. jobvr == 'V') then
      call ${VPREC}$geev(jobvl, jobvr, n, a, lda, wr, wi, vl, ldvl, vr, ldvr, work, int_idealwork,&
          & info)
    else if (jobvl == 'V' .and. jobvr == 'N') then
      call ${VPREC}$geev(jobvl, jobvr, n, a, lda, wr, wi, vl, ldvl, dummyvr, ldvr, work,&
          & int_idealwork, info)
    else if (jobvl == 'N' .and. jobvr == 'V') then
      call ${VPREC}$geev(jobvl, jobvr, n, a, lda, wr, wi, dummyvl, ldvl, vr, ldvr, work,&
          & int_idealwork, info)
    else if (jobvl == 'N' .and. jobvr == 'N') then
      call ${VPREC}$geev(jobvl, jobvr, n, a, lda, wr, wi, dummyvl, ldvl, dummyvr, ldvr, work,&
          & int_idealwork, info)
    end if

    if (info/=0) then
      if (info<0) then
        @:FORMATTED_ERROR_HANDLING(err, info, "(A,I0)", 'Failure in diagonalisation routine&
            & ${VPREC}$geev, illegal argument at position ', info)
      else
        @:FORMATTED_ERROR_HANDLING(err, info, "(A,I0,A)", 'Failure in diagonalisation routine&
            & ${VPREC}$geev, diagonal element ', info, ' did not converge to zero.')
      endif
    endif

  end subroutine ${DTYPE}$_${NAME}$

#:endfor

end module dftbp_eigensolver
