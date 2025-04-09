!--------------------------------------------------------------------------------------------------!
!  DFTB+: general package for performing fast atomistic simulations                                !
!  Copyright (C) 2006 - 2025  DFTB+ developers group                                               !
!                                                                                                  !
!  See the LICENSE file for terms of usage and distribution.                                       !
!--------------------------------------------------------------------------------------------------!

!> Interface wrappers for LAPACK routines. See the <a href="http://www.netlib.org/lapack/">lapack
!! project documentation</a> for documentation of the routines.
module dftbp_extlibs_lapack
  use dftbp_common_accuracy, only : rsp, rdp
  implicit none

  interface
  #:for LAPACK_ROUTINE, KIND in [("ssyev", "rsp"), ("dsyev", "rdp")]

    !> Real symmetric eigensolver
    subroutine ${LAPACK_ROUTINE}$(jobz, uplo, nn, aa, lda, ww, work, lwork, info)
      import :: ${KIND}$
      character, intent(in) :: jobz
      character, intent(in) :: uplo
      integer, intent(in) :: nn
      integer, intent(in) :: lda
      real(${KIND}$), intent(inout) :: aa(lda, *)
      real(${KIND}$), intent(out) :: ww(*)
      real(${KIND}$), intent(inout) :: work(*)
      integer, intent(in) :: lwork
      integer, intent(out) :: info
    end subroutine ${LAPACK_ROUTINE}$

  #:endfor
  end interface


  interface
  #:for LAPACK_ROUTINE, KIND in [("cheev", "rsp"), ("zheev", "rdp")]

    !> Complex symmetric eigensolver
    subroutine ${LAPACK_ROUTINE}$(jobz, uplo, nn, aa, lda, ww, work, lwork, rwork, info)
      import :: ${KIND}$
      character, intent(in) :: jobz
      character, intent(in) :: uplo
      integer, intent(in) :: nn
      integer, intent(in) :: lda
      complex(${KIND}$), intent(inout) :: aa(lda, *)
      real(${KIND}$), intent(out) :: ww(*)
      complex(${KIND}$), intent(inout) :: work(*)
      integer, intent(in) :: lwork
      real(${KIND}$), intent(inout) :: rwork(*)
      integer, intent(out) :: info
    end subroutine ${LAPACK_ROUTINE}$

  #:endfor
  end interface


  interface
  #:for LAPACK_ROUTINE, KIND in [("ssyevd", "rsp"), ("dsyevd", "rdp")]

    !> Real symmetric eigensolver (relatively robust)
    subroutine ${LAPACK_ROUTINE}$(jobz, uplo, nn, aa, lda, ww, work, lwork, iwork, liwork, info)
      import :: ${KIND}$
      character, intent(in) :: jobz
      character, intent(in) :: uplo
      integer, intent(in) :: nn
      integer, intent(in) :: lda
      real(${KIND}$), intent(inout) :: aa(lda, *)
      real(${KIND}$), intent(out) :: ww(*)
      real(${KIND}$), intent(inout) :: work(*)
      integer, intent(in) :: lwork
      integer, intent(inout) :: iwork(*)
      integer, intent(in) :: liwork
      integer, intent(out) :: info
    end subroutine ${LAPACK_ROUTINE}$

  #:endfor
  end interface


  interface
  #:for LAPACK_ROUTINE, KIND in [("cheevd", "rsp"), ("zheevd", "rdp")]

    !> Complex hermitian eigensolver (relatively robust)
    subroutine ${LAPACK_ROUTINE}$(jobz, uplo, nn, aa, lda, ww, work, lwork, rwork, lrwork, iwork,&
        & liwork, info)
      import :: ${KIND}$
      character, intent(in) :: jobz
      character, intent(in) :: uplo
      integer, intent(in) :: nn
      integer, intent(in) :: lda
      complex(${KIND}$), intent(inout) :: aa(lda, *)
      real(${KIND}$), intent(out) :: ww(*)
      complex(${KIND}$), intent(inout) :: work(*)
      integer, intent(in) :: lwork
      real(${KIND}$), intent(inout) :: rwork(*)
      integer, intent(in) :: lrwork
      integer, intent(inout) :: iwork(*)
      integer, intent(in) :: liwork
      integer, intent(out) :: info
    end subroutine ${LAPACK_ROUTINE}$

  #:endfor
  end interface


  interface
  #:for LAPACK_ROUTINE, KIND in [("ssyevr", "rsp"), ("dsyevr", "rdp")]

    !> Real symmetric eigensolver (relatively robust)
    subroutine ${LAPACK_ROUTINE}$(jobz, range, uplo, nn, aa, lda, vl, vu, il, iu, abstol, mm, ww,&
        & zz, ldz, isuppz, work, lwork, iwork, liwork, info)
      import :: ${KIND}$
      character, intent(in) :: jobz
      character, intent(in) :: range
      character, intent(in) :: uplo
      integer, intent(in) :: nn
      integer, intent(in) :: lda
      real(${KIND}$), intent(inout) :: aa(lda, *)
      real(${KIND}$), intent(in) :: vl
      real(${KIND}$), intent(in) :: vu
      integer, intent(in) :: il
      integer, intent(in) :: iu
      real(${KIND}$), intent(in) :: abstol
      integer, intent(out) :: mm
      real(${KIND}$), intent(out) :: ww(*)
      integer, intent(in) :: ldz
      real(${KIND}$), intent(out) :: zz(ldz, *)
      integer, intent(out) :: isuppz(*)
      real(${KIND}$), intent(inout) :: work(*)
      integer, intent(in) :: lwork
      integer, intent(inout) :: iwork(*)
      integer, intent(in) :: liwork
      integer, intent(out) :: info
    end subroutine ${LAPACK_ROUTINE}$

  #:endfor
  end interface


  interface
  #:for LAPACK_ROUTINE, KIND in [("cheevr", "rsp"), ("zheevr", "rdp")]

    !> Complex hermitian eigensolver (relatively robust)
    subroutine ${LAPACK_ROUTINE}$(jobz, range, uplo, nn, aa, lda, vl, vu, il, iu, abstol,&
        & mm, ww, zz, ldz, isuppz, work, lwork, rwork, lrwork, iwork, liwork, info)
      import :: ${KIND}$
      character, intent(in) :: jobz
      character, intent(in) :: range
      character, intent(in) :: uplo
      integer, intent(in) :: nn
      integer, intent(in) :: lda
      complex(${KIND}$), intent(inout) :: aa(lda, *)
      real(${KIND}$), intent(in) :: vl
      real(${KIND}$), intent(in) :: vu
      integer, intent(in) :: il
      integer, intent(in) :: iu
      real(${KIND}$), intent(in) :: abstol
      integer, intent(out) :: mm
      real(${KIND}$), intent(out) :: ww(*)
      integer, intent(in) :: ldz
      complex(${KIND}$), intent(out) :: zz(ldz, *)
      integer, intent(out) :: isuppz(*)
      complex(${KIND}$), intent(inout) :: work(*)
      integer, intent(in) :: lwork
      real(${KIND}$), intent(inout) :: rwork(*)
      integer, intent(in) :: lrwork
      integer, intent(inout) :: iwork(*)
      integer, intent(in) :: liwork
      integer, intent(out) :: info
    end subroutine ${LAPACK_ROUTINE}$

  #:endfor
  end interface


  interface
  #:for LAPACK_ROUTINE, KIND in [("ssygv", "rsp"), ("dsygv", "rdp")]

    !> Real symmetric generalized eigensolver
    subroutine ${LAPACK_ROUTINE}$(itype, jobz, uplo, nn, aa, lda, bb, ldb, ww, work, lwork, info)
      import :: ${KIND}$
      integer, intent(in) :: itype
      character, intent(in) :: jobz
      character, intent(in) :: uplo
      integer, intent(in) :: nn
      integer, intent(in) :: lda
      real(${KIND}$), intent(inout) :: aa(lda, *)
      integer, intent(in) :: ldb
      real(${KIND}$), intent(inout) :: bb(ldb, *)
      real(${KIND}$), intent(out) :: ww(*)
      real(${KIND}$), intent(inout) :: work(*)
      integer, intent(in) :: lwork
      integer, intent(out) :: info
    end subroutine ${LAPACK_ROUTINE}$

  #:endfor
  end interface


  interface
  #:for LAPACK_ROUTINE, KIND in [("chegv", "rsp"), ("zhegv", "rdp")]

    !> Complex hermitian generalized eigensolver
    subroutine ${LAPACK_ROUTINE}$(itype, jobz, uplo, nn, aa, lda, bb, ldb, ww, work, lwork, rwork,&
          & info)
      import :: ${KIND}$
      integer, intent(in) :: itype
      character, intent(in) :: jobz
      character, intent(in) :: uplo
      integer, intent(in) :: nn
      integer, intent(in) :: lda
      complex(${KIND}$), intent(inout) :: aa(lda, *)
      integer, intent(in) :: ldb
      complex(${KIND}$), intent(inout) :: bb(ldb, *)
      real(${KIND}$), intent(out) :: ww(*)
      complex(${KIND}$), intent(inout) :: work(*)
      integer, intent(in) :: lwork
      real(${KIND}$), intent(inout) :: rwork(*)
      integer, intent(out) :: info
    end subroutine ${LAPACK_ROUTINE}$

  #:endfor
  end interface


  interface
    #:for LAPACK_ROUTINE, KIND in [("ssygvd", "rsp"), ("dsygvd", "rdp")]

    !> Real symmetric generalized eigensolver (divide and conquer)
    subroutine ${LAPACK_ROUTINE}$(itype, jobz, uplo, nn, aa, lda, bb, ldb, ww, work, lwork, iwork,&
          & liwork, info)
      import :: ${KIND}$
      integer, intent(in) :: itype
      character, intent(in) :: jobz
      character, intent(in) :: uplo
      integer, intent(in) :: nn
      integer, intent(in) :: lda
      real(${KIND}$), intent(inout) :: aa(lda, *)
      integer, intent(in) :: ldb
      real(${KIND}$), intent(inout) :: bb(ldb, *)
      real(${KIND}$), intent(out) :: ww(*)
      real(${KIND}$), intent(inout) :: work(*)
      integer, intent(in) :: lwork
      integer, intent(inout) :: iwork(*)
      integer, intent(in) :: liwork
      integer, intent(out) :: info
    end subroutine ${LAPACK_ROUTINE}$
  #:endfor
  end interface


  interface
  #:for LAPACK_ROUTINE, KIND in [("chegvd", "rsp"), ("zhegvd", "rdp")]

    !> Complex hermitian generalized eigensolver (divide and conquer)
    subroutine ${LAPACK_ROUTINE}$(itype, jobz, uplo, nn, aa, lda, bb, ldb, ww, work, lwork, rwork,&
          & lrwork, iwork, liwork, info)
      import :: ${KIND}$
      integer, intent(in) :: itype
      character, intent(in) :: jobz
      character, intent(in) :: uplo
      integer, intent(in) :: nn
      integer, intent(in) :: lda
      complex(${KIND}$), intent(inout) :: aa(lda, *)
      integer, intent(in) :: ldb
      complex(${KIND}$), intent(inout) :: bb(ldb, *)
      real(${KIND}$), intent(out) :: ww(*)
      complex(${KIND}$), intent(inout) :: work(*)
      integer, intent(in) :: lwork
      real(${KIND}$), intent(inout) :: rwork(*)
      integer, intent(in) :: lrwork
      integer, intent(inout) :: iwork(*)
      integer, intent(in) :: liwork
      integer, intent(out) :: info
    end subroutine ${LAPACK_ROUTINE}$

  #:endfor
  end interface




  interface
  #:for LAPACK_ROUTINE, TYPE, KIND in &
      & [("spotrf", "real", "rsp"), ("dpotrf", "real", "rdp"),&
      & ("cpotrf", "complex", "rsp"), ("zpotrf", "complex", "rdp")]

    !> Cholesky factorization of symmetric / hermitian positive definite matrix
    subroutine ${LAPACK_ROUTINE}$(uplo, nn, aa, lda, info)
      import :: ${KIND}$
      character, intent(in) :: uplo
      integer, intent(in) :: nn
      integer, intent(in) :: lda
      ${TYPE}$(${KIND}$), intent(inout) :: aa(lda, *)
      integer, intent(out) :: info
    end subroutine ${LAPACK_ROUTINE}$

  #:endfor
  end interface


  interface
  #:for LAPACK_ROUTINE, TYPE, KIND in &
      & [("ssygst", "real", "rsp"), ("dsygst", "real", "rdp"),&
      & ("chegst", "complex", "rsp"), ("zhegst", "complex", "rdp")]

    !> Reduces generalized eigenproblem to standard form
    subroutine ${LAPACK_ROUTINE}$(itype, uplo, nn, aa, lda, bb, ldb, info)
      import :: ${KIND}$
      integer, intent(in) :: itype
      character, intent(in) :: uplo
      integer, intent(in) :: nn
      integer, intent(in) :: lda
      ${TYPE}$(${KIND}$), intent(inout) :: aa(lda, *)
      integer, intent(in) :: ldb
      ${TYPE}$(${KIND}$), intent(in) :: bb(ldb, *)
      integer, intent(out) :: info
    end subroutine ${LAPACK_ROUTINE}$

  #:endfor
  end interface


  interface
  #:for LAPACK_ROUTINE, KIND in [("sgesv", "rsp"), ("dgesv", "rdp")]

    !> Solve overdetermined or underdetermined real linear systems
    subroutine ${LAPACK_ROUTINE}$(nn, nrhs, aa, lda, ipiv, bb, ldb, info)
      import :: ${KIND}$
      integer, intent(in) :: nn
      integer, intent(in) :: nrhs
      integer, intent(in) :: lda
      real(${KIND}$), intent(inout) :: aa(lda, *)
      integer, intent(out) :: ipiv(*)
      integer, intent(in) :: ldb
      real(${KIND}$), intent(inout) :: bb(ldb, *)
      integer, intent(out) :: info
    end subroutine ${LAPACK_ROUTINE}$
  #:endfor
  end interface


  interface
  #:for LAPACK_ROUTINE, TYPE, KIND in &
      & [("sgetrf", "real", "rsp"), ("dgetrf", "real", "rdp"),&
      &  ("cgetrf", "complex", "rsp"), ("zgetrf", "complex", "rdp")]

    !> Computes LU factorization of real matrix
    subroutine ${LAPACK_ROUTINE}$(mm, nn, aa, lda, ipiv, info)
      import :: ${KIND}$
      integer, intent(in) :: mm
      integer, intent(in) :: nn
      integer, intent(in) :: lda
      ${TYPE}$(${KIND}$), intent(inout) :: aa(lda, *)
      integer, intent(out) :: ipiv(*)
      integer, intent(out) :: info
    end subroutine ${LAPACK_ROUTINE}$
  #:endfor
  end interface


  interface
  #:for LAPACK_ROUTINE, KIND in [("sgetri", "rsp"), ("dgetri", "rdp")]

    !> Computes inverse of a real matrix using LU factorisation
    subroutine ${LAPACK_ROUTINE}$(nn, aa, lda, ipiv, work, lwork, info)
      import :: ${KIND}$
      integer, intent(in) :: nn
      integer, intent(in) :: lda
      real(${KIND}$), intent(inout) :: aa(lda, *)
      integer, intent(in) :: ipiv(*)
      real(${KIND}$), intent(inout) :: work(*)
      integer, intent(in) :: lwork
      integer, intent(out) :: info
    end subroutine ${LAPACK_ROUTINE}$

  #:endfor
  end interface


  interface
  #:for LAPACK_ROUTINE, TYPE_KIND in [("ssytrf", "rsp"), ("dsytrf", "rdp")]

    !> Factorise a real symmetric matrix as A = U*D*U**T or A = L*D*L**T
    subroutine ${LAPACK_ROUTINE}$(uplo, nn, aa, lda, ipiv, work, lwork, info)
      import :: ${KIND}$
      character, intent(in) :: uplo
      integer, intent(in) :: nn
      integer, intent(in) :: lda
      real(${KIND}$), intent(inout) :: aa(lda, *)
      integer, intent(out) :: ipiv(*)
      real(${KIND}$), intent(inout) :: work(*)
      integer, intent(in) :: lwork
      integer, intent(out) :: info
    end subroutine ${LAPACK_ROUTINE}$

  #:endfor
  end interface


  interface
  #:for LAPACK_ROUTINE, TYPE, KIND in &
      & [("slarnv", "real", "rsp"), ("dlarnv", "real", "rdp"),&
      &  ("clarnv", "complex", "rsp"), ("zlarnv", "complex", "rdp")]

    !> Returns a vector of real random numbers from a uniform or normal distribution
    subroutine ${LAPACK_ROUTINE}$(idist, iseed, nn, xx)
      import :: ${KIND}$
      integer, intent(in) :: idist
      integer, intent(inout) :: iseed(4)
      integer, intent(in) :: nn
      ${TYPE}$(${KIND}$), intent(out) :: xx(*)
    end subroutine ${LAPACK_ROUTINE}$

  #:endfor
  end interface


  interface
    !> Provides problem-dependent LAPACK routine parameters for the local environment
    function ilaenv(ispec, name, opts, n1, n2, n3, n4)
      integer, intent(in) :: ispec
      character, intent(in) :: name
      character, intent(in) :: opts
      integer, intent(in) :: n1
      integer, intent(in) :: n2
      integer, intent(in) :: n3
      integer, intent(in) :: n4
      integer :: ilaenv
    end function ilaenv
  end interface


  interface
  #:for LAPACK_ROUTINE, KIND in [("slamch", "rsp"), ("dlamch", "rdp")]

    !> Machine parameters
    function ${LAPACK_ROUTINE}$(cmach)
      import :: ${KIND}$
      character, intent(in) :: cmach
      real(${KIND}$) :: ${LAPACK_ROUTINE}$
    end function ${LAPACK_ROUTINE}$

  #:endfor
  end interface


  interface
    !> Error handler for the LAPACK routines
    subroutine xerbla(srname, info)
      character(6), intent(in) :: srname
      integer, intent(in) :: info
    end subroutine xerbla
  end interface


  interface
  #:for LAPACK_ROUTINE, KIND in [("rgesvd", "rsp"), ("dgesvd", "rdp")]

  !> Real singular value decomposition
    subroutine ${LAPACK_ROUTINE}$(jobu, jobvt, m, n, a, lda, s, u, ldu, vt, ldvt, work, lwork, info)
      import :: ${KIND}$
      character, intent(in) :: jobvt
      character, intent(in) :: jobu
      integer, intent(in) :: m
      integer, intent(in) :: n
      integer, intent(in) :: lda
      integer, intent(in) :: ldu
      integer, intent(in) :: ldvt
      real(${KIND}$), intent(inout) :: a(lda,*)
      real(${KIND}$), intent(out) :: s(*)
      real(${KIND}$), intent(out) :: u(ldu,*)
      real(${KIND}$), intent(out) :: vt(ldvt,*)
      real(${KIND}$), intent(out) :: work(*)
      integer, intent(in) :: lwork
      integer, intent(in) :: info
    end subroutine ${LAPACK_ROUTINE}$

  #:endfor
  end interface


  interface
  #:for LAPACK_ROUTINE, KIND in [("cgesvd", "rsp"), ("zgesvd", "rdp")]

    !> Complex singular value decomposition
    subroutine ${LAPACK_ROUTINE}$(jobu, jobvt, m, n, a, lda, s, u, ldu, vt, ldvt, work, lwork, rwork, info)
      import :: ${KIND}$
      character, intent(in) :: jobvt
      character, intent(in) :: jobu
      integer, intent(in) :: m
      integer, intent(in) :: n
      integer, intent(in) :: lda
      integer, intent(in) :: ldu
      integer, intent(in) :: ldvt
      complex(${KIND}$), intent(inout) :: a(lda,*)
      real(${KIND}$), intent(out) :: s(*)
      real(${KIND}$), intent(out) :: rwork(*)
      complex(${KIND}$), intent(out) :: u(ldu,*)
      complex(${KIND}$), intent(out) :: vt(ldvt,*)
      complex(${KIND}$), intent(out) :: work(*)
      integer, intent(in) :: lwork
      integer, intent(in) :: info
    end subroutine ${LAPACK_ROUTINE}$

  #:endfor
  end interface


  interface
  #:for LAPACK_ROUTINE, KIND in [("sgeev", "rsp"), ("dgeev", "rdp")]

    !> general matrix eigensolver
    subroutine ${LAPACK_ROUTINE}$(jobvl, jobvr, n, a, lda, wr, wi, vl, ldvl, vr, ldvr, work,&
          & lwork, info)
      import :: ${KIND}$
      character, intent(in) :: jobvl
      character, intent(in) :: jobvr
      integer, intent(in) :: n
      real(${KIND}$), intent(inout) :: a(lda,*)
      integer, intent(in) :: lda
      real(${KIND}$), intent(out) :: wr(*)
      real(${KIND}$), intent(out) :: wi(*)
      real(${KIND}$), intent(out) :: vl(ldvl, *)
      integer, intent(in) :: ldvl
      real(${KIND}$), intent(out) :: vr(ldvr, *)
      integer, intent(in) :: ldvr
      real(${KIND}$), intent(out) :: work(*)
      integer, intent(in) :: lwork
      integer, intent(in) :: info
    end subroutine ${LAPACK_ROUTINE}$

  #:endfor
  end interface


  interface
  #:for LAPACK_ROUTINE, TYPE, KIND in &
      & [("strsm", "real", "rsp"), ("dtrsm", "real", "rdp"),&
      &  ("ctrsm", "complex", "rsp"), ("ztrsm", "complex", "rdp")]

    !> Triangular solve with multiple right-hand sides
    subroutine ${LAPACK_ROUTINE}$(side, uplo, transa, diag, m, n, alpha, a, lda, b, ldb)
      import :: ${KIND}$
      character, intent(in) :: side
      character, intent(in) :: uplo
      character, intent(in) :: transa
      character, intent(in) :: diag
      integer, intent(in) :: m
      integer, intent(in) :: n
      ${TYPE}$(${KIND}$), intent(in) :: alpha
      integer, intent(in) :: lda
      ${TYPE}$(${KIND}$), intent(in) :: a(lda, *)
      integer, intent(in) :: ldb
      ${TYPE}$(${KIND}$), intent(inout) :: b(ldb, *)
    end subroutine ${LAPACK_ROUTINE}$

  #:endfor
  end interface


  interface
  #:for LAPACK_ROUTINE, TYPE, KIND in &
      & [("strmm", "real", "rsp"), ("dtrmm", "real", "rdp"),&
      &  ("ctrmm", "complex", "rsp"), ("ztrmm", "complex", "rdp")]

    !> Triangular matrix-matrix multiplication
    subroutine ${LAPACK_ROUTINE}$(side, uplo, transa, diag, m, n, alpha, a, lda, b, ldb)
      import :: ${KIND}$
      character, intent(in) :: side
      character, intent(in) :: uplo
      character, intent(in) :: transa
      character, intent(in) :: diag
      integer, intent(in) :: m
      integer, intent(in) :: n
      ${TYPE}$(${KIND}$), intent(in) :: alpha
      integer, intent(in) :: lda
      ${TYPE}$(${KIND}$), intent(in) :: a(lda, *)
      integer, intent(in) :: ldb
      ${TYPE}$(${KIND}$), intent(inout) :: b(ldb, *)
    end subroutine ${LAPACK_ROUTINE}$

  #:endfor
  end interface

end module dftbp_extlibs_lapack
