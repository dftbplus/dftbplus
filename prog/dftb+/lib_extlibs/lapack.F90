!--------------------------------------------------------------------------------------------------!
!  DFTB+: general package for performing fast atomistic simulations                                !
!  Copyright (C) 2017  DFTB+ developers group                                                      !
!                                                                                                  !
!  See the LICENSE file for terms of usage and distribution.                                       !
!--------------------------------------------------------------------------------------------------!

!> Interface wrapper for the lapack routines.
module lapack
  use accuracy, only : rsp, rdp
  implicit none
  public

  interface

    subroutine ssyev(jobz, uplo, nn, aa, lda, ww, work, lwork, info)
      import rsp
      character, intent(in) :: jobz, uplo
      integer, intent(in) :: nn
      integer, intent(in) :: lda
      real(rsp), intent(inout) :: aa(lda, *)
      real(rsp), intent(out) :: ww(*)
      real(rsp), intent(inout) :: work(*)
      integer, intent(in) :: lwork
      integer, intent(out) :: info
    end subroutine ssyev

    subroutine dsyev(jobz, uplo, nn, aa, lda, ww, work, lwork, info)
      import rdp
      character, intent(in) :: jobz, uplo
      integer, intent(in) :: nn
      integer, intent(in) :: lda
      real(rdp), intent(inout) :: aa(lda, *)
      real(rdp), intent(out) :: ww(*)
      real(rdp), intent(inout) :: work(*)
      integer, intent(in) :: lwork
      integer, intent(out) :: info
    end subroutine dsyev

    subroutine cheev(jobz, uplo, nn, aa, lda, ww, work, lwork, rwork, info)
      import rsp
      character, intent(in) :: jobz, uplo
      integer, intent(in) :: nn
      integer, intent(in) :: lda
      complex(rsp), intent(inout) :: aa(lda, *)
      real(rsp), intent(out) :: ww(*)
      complex(rsp), intent(inout) :: work(*)
      integer, intent(in) :: lwork
      real(rsp), intent(inout) :: rwork(*)
      integer, intent(out) :: info
    end subroutine cheev

    subroutine zheev(jobz, uplo, nn, aa, lda, ww, work, lwork, rwork, info)
      import rdp
      character, intent(in) :: jobz, uplo
      integer, intent(in) :: nn
      integer, intent(in) :: lda
      complex(rdp), intent(inout) :: aa(lda, *)
      real(rdp), intent(out) :: ww(*)
      complex(rdp), intent(inout) :: work(*)
      integer, intent(in) :: lwork
      real(rdp), intent(inout) :: rwork(*)
      integer, intent(out) :: info
    end subroutine zheev

    subroutine ssygv(itype, jobz, uplo, nn, aa, lda, bb, ldb, ww, work, lwork,&
        & info)
      import rsp
      integer, intent(in) :: itype
      character, intent(in) :: jobz, uplo
      integer, intent(in) :: nn
      integer, intent(in) :: lda
      real(rsp), intent(inout) :: aa(lda, *)
      integer, intent(in) :: ldb
      real(rsp), intent(inout) :: bb(ldb, *)
      real(rsp), intent(out) :: ww(*)
      real(rsp), intent(inout) :: work(*)
      integer, intent(in) :: lwork
      integer, intent(out) :: info
    end subroutine ssygv

    subroutine dsygv(itype, jobz, uplo, nn, aa, lda, bb, ldb, ww, work, lwork,&
        & info)
      import rdp
      integer, intent(in) :: itype
      character, intent(in) :: jobz, uplo
      integer, intent(in) :: nn
      integer, intent(in) :: lda
      real(rdp), intent(inout) :: aa(lda, *)
      integer, intent(in) :: ldb
      real(rdp), intent(inout) :: bb(ldb, *)
      real(rdp), intent(out) :: ww(*)
      real(rdp), intent(inout) :: work(*)
      integer, intent(in) :: lwork
      integer, intent(out) :: info
    end subroutine dsygv

    subroutine chegv(itype, jobz, uplo, nn, aa, lda, bb, ldb, ww, work, lwork,&
        & rwork, info)
      import rsp
      integer, intent(in) :: itype
      character, intent(in) :: jobz, uplo
      integer, intent(in) :: nn
      integer, intent(in) :: lda
      complex(rsp), intent(inout) :: aa(lda, *)
      integer, intent(in) :: ldb
      complex(rsp), intent(inout) :: bb(ldb, *)
      real(rsp), intent(out) :: ww(*)
      complex(rsp), intent(inout) :: work(*)
      integer, intent(in) :: lwork
      real(rsp), intent(inout) :: rwork(*)
      integer, intent(out) :: info
    end subroutine chegv

    subroutine zhegv(itype, jobz, uplo, nn, aa, lda, bb, ldb, ww, work, lwork,&
        & rwork, info)
      import rdp
      integer, intent(in) :: itype
      character, intent(in) :: jobz, uplo
      integer, intent(in) :: nn
      integer, intent(in) :: lda
      complex(rdp), intent(inout) :: aa(lda, *)
      integer, intent(in) :: ldb
      complex(rdp), intent(inout) :: bb(ldb, *)
      real(rdp), intent(out) :: ww(*)
      complex(rdp), intent(inout) :: work(*)
      integer, intent(in) :: lwork
      real(rdp), intent(inout) :: rwork(*)
      integer, intent(out) :: info
    end subroutine zhegv

    subroutine ssygvd(itype, jobz, uplo, nn, aa, lda, bb, ldb, ww, work,&
        & lwork, iwork, liwork, info)
      import rsp
      integer, intent(in) :: itype
      character, intent(in) :: jobz, uplo
      integer, intent(in) :: nn
      integer, intent(in) :: lda
      real(rsp), intent(inout) :: aa(lda, *)
      integer, intent(in) :: ldb
      real(rsp), intent(inout) :: bb(ldb, *)
      real(rsp), intent(out) :: ww(*)
      real(rsp), intent(inout) :: work(*)
      integer, intent(in) :: lwork
      integer, intent(inout) :: iwork(*)
      integer, intent(in) :: liwork
      integer, intent(out) :: info
    end subroutine ssygvd

    subroutine dsygvd(itype, jobz, uplo, nn, aa, lda, bb, ldb, ww, work,&
        & lwork, iwork, liwork, info)
      import rdp
      integer, intent(in) :: itype
      character, intent(in) :: jobz, uplo
      integer, intent(in) :: nn
      integer, intent(in) :: lda
      real(rdp), intent(inout) :: aa(lda, *)
      integer, intent(in) :: ldb
      real(rdp), intent(inout) :: bb(ldb, *)
      real(rdp), intent(out) :: ww(*)
      real(rdp), intent(inout) :: work(*)
      integer, intent(in) :: lwork
      integer, intent(inout) :: iwork(*)
      integer, intent(in) :: liwork
      integer, intent(out) :: info
    end subroutine dsygvd

    subroutine chegvd(itype, jobz, uplo, nn, aa, lda, bb, ldb, ww, work,&
        & lwork, rwork, lrwork, iwork, liwork, info)
      import rsp
      integer, intent(in) :: itype
      character, intent(in) :: jobz, uplo
      integer, intent(in) :: nn
      integer, intent(in) :: lda
      complex(rsp), intent(inout) :: aa(lda, *)
      integer, intent(in) :: ldb
      complex(rsp), intent(inout) :: bb(ldb, *)
      real(rsp), intent(out) :: ww(*)
      complex(rsp), intent(inout) :: work(*)
      integer, intent(in) :: lwork
      real(rsp), intent(inout) :: rwork(*)
      integer, intent(in) :: lrwork
      integer, intent(inout) :: iwork(*)
      integer, intent(in) :: liwork
      integer, intent(out) :: info
    end subroutine chegvd

    subroutine zhegvd(itype, jobz, uplo, nn, aa, lda, bb, ldb, ww, work,&
        & lwork, rwork, lrwork, iwork, liwork, info)
      import rdp
      integer, intent(in) :: itype
      character, intent(in) :: jobz, uplo
      integer, intent(in) :: nn
      integer, intent(in) :: lda
      complex(rdp), intent(inout) :: aa(lda, *)
      integer, intent(in) :: ldb
      complex(rdp), intent(inout) :: bb(ldb, *)
      real(rdp), intent(out) :: ww(*)
      complex(rdp), intent(inout) :: work(*)
      integer, intent(in) :: lwork
      real(rdp), intent(inout) :: rwork(*)
      integer, intent(in) :: lrwork
      integer, intent(inout) :: iwork(*)
      integer, intent(in) :: liwork
      integer, intent(out) :: info
    end subroutine zhegvd

    subroutine ssyevr(jobz, range, uplo, nn, aa, lda, vl, vu, il, iu, abstol,&
        & mm, ww, zz, ldz, isuppz, work, lwork, iwork, liwork, info)
      import rsp
      character, intent(in) :: jobz, range, uplo
      integer, intent(in) :: nn
      integer, intent(in) :: lda
      real(rsp), intent(inout) :: aa(lda, *)
      real(rsp), intent(in) :: vl, vu
      integer, intent(in) :: il, iu
      real(rsp), intent(in) :: abstol
      integer, intent(out) :: mm
      real(rsp), intent(out) :: ww(*)
      integer, intent(in) :: ldz
      real(rsp), intent(out) :: zz(ldz, *)
      integer, intent(out) :: isuppz(*)
      real(rsp), intent(inout) :: work(*)
      integer, intent(in) :: lwork
      integer, intent(inout) :: iwork(*)
      integer, intent(in) :: liwork
      integer, intent(out) :: info
    end subroutine ssyevr

    subroutine dsyevr(jobz, range, uplo, nn, aa, lda, vl, vu, il, iu, abstol,&
        & mm, ww, zz, ldz, isuppz, work, lwork, iwork, liwork, info)
      import rdp
      character, intent(in) :: jobz, range, uplo
      integer, intent(in) :: nn
      integer, intent(in) :: lda
      real(rdp), intent(inout) :: aa(lda, *)
      real(rdp), intent(in) :: vl, vu
      integer, intent(in) :: il, iu
      real(rdp), intent(in) :: abstol
      integer, intent(out) :: mm
      real(rdp), intent(out) :: ww(*)
      integer, intent(in) :: ldz
      real(rdp), intent(out) :: zz(ldz, *)
      integer, intent(out) :: isuppz(*)
      real(rdp), intent(inout) :: work(*)
      integer, intent(in) :: lwork
      integer, intent(inout) :: iwork(*)
      integer, intent(in) :: liwork
      integer, intent(out) :: info
    end subroutine dsyevr

    subroutine cheevr(jobz, range, uplo, nn, aa, lda, vl, vu, il, iu, abstol,&
        & mm, ww, zz, ldz, isuppz, work, lwork, rwork, lrwork, iwork, liwork,&
        & info)
      import rsp
      character, intent(in) :: jobz, range, uplo
      integer, intent(in) :: nn
      integer, intent(in) :: lda
      complex(rsp), intent(inout) :: aa(lda, *)
      real(rsp), intent(in) :: vl, vu
      integer, intent(in) :: il, iu
      real(rsp), intent(in) :: abstol
      integer, intent(out) :: mm
      real(rsp), intent(out) :: ww(*)
      integer, intent(in) :: ldz
      complex(rsp), intent(out) :: zz(ldz, *)
      integer, intent(out) :: isuppz(*)
      complex(rsp), intent(inout) :: work(*)
      integer, intent(in) :: lwork
      real(rsp), intent(inout) :: rwork(*)
      integer, intent(in) :: lrwork
      integer, intent(inout) :: iwork(*)
      integer, intent(in) :: liwork
      integer, intent(out) :: info
    end subroutine cheevr

    subroutine zheevr(jobz, range, uplo, nn, aa, lda, vl, vu, il, iu, abstol,&
        & mm, ww, zz, ldz, isuppz, work, lwork, rwork, lrwork, iwork, liwork,&
        & info)
      import rdp
      character, intent(in) :: jobz, range, uplo
      integer, intent(in) :: nn
      integer, intent(in) :: lda
      complex(rdp), intent(inout) :: aa(lda, *)
      real(rdp), intent(in) :: vl, vu
      integer, intent(in) :: il, iu
      real(rdp), intent(in) :: abstol
      integer, intent(out) :: mm
      real(rdp), intent(out) :: ww(*)
      integer, intent(in) :: ldz
      complex(rdp), intent(out) :: zz(ldz, *)
      integer, intent(out) :: isuppz(*)
      complex(rdp), intent(inout) :: work(*)
      integer, intent(in) :: lwork
      real(rdp), intent(inout) :: rwork(*)
      integer, intent(in) :: lrwork
      integer, intent(inout) :: iwork(*)
      integer, intent(in) :: liwork
      integer, intent(out) :: info
    end subroutine zheevr

    subroutine spotrf(uplo, nn, aa, lda, info)
      import rsp
      character, intent(in) :: uplo
      integer, intent(in) :: nn
      integer, intent(in) :: lda
      real(rsp), intent(inout) :: aa(lda, *)
      integer, intent(out) :: info
    end subroutine spotrf

    subroutine dpotrf(uplo, nn, aa, lda, info)
      import rdp
      character, intent(in) :: uplo
      integer, intent(in) :: nn
      integer, intent(in) :: lda
      real(rdp), intent(inout) :: aa(lda, *)
      integer, intent(out) :: info
    end subroutine dpotrf

    subroutine cpotrf(uplo, nn, aa, lda, info)
      import rsp
      character, intent(in) :: uplo
      integer, intent(in) :: nn
      integer, intent(in) :: lda
      complex(rsp), intent(inout) :: aa(lda, *)
      integer, intent(out) :: info
    end subroutine cpotrf

    subroutine zpotrf(uplo, nn, aa, lda, info)
      import rdp
      character, intent(in) :: uplo
      integer, intent(in) :: nn
      integer, intent(in) :: lda
      complex(rdp), intent(inout) :: aa(lda, *)
      integer, intent(out) :: info
    end subroutine zpotrf

    subroutine ssygst(itype, uplo, nn, aa, lda, bb, ldb, info)
      import rsp
      integer, intent(in) :: itype
      character, intent(in) :: uplo
      integer, intent(in) :: nn
      integer, intent(in) :: lda
      real(rsp), intent(inout) :: aa(lda, *)
      integer, intent(in) :: ldb
      real(rsp), intent(in) :: bb(ldb, *)
      integer, intent(out) :: info
    end subroutine ssygst

    subroutine dsygst(itype, uplo, nn, aa, lda, bb, ldb, info)
      import rdp
      integer, intent(in) :: itype
      character, intent(in) :: uplo
      integer, intent(in) :: nn
      integer, intent(in) :: lda
      real(rdp), intent(inout) :: aa(lda, *)
      integer, intent(in) :: ldb
      real(rdp), intent(in) :: bb(ldb, *)
      integer, intent(out) :: info
    end subroutine dsygst

    subroutine chegst(itype, uplo, nn, aa, lda, bb, ldb, info)
      import rsp
      integer, intent(in) :: itype
      character, intent(in) :: uplo
      integer, intent(in) :: nn
      integer, intent(in) :: lda
      complex(rsp), intent(inout) :: aa(lda, *)
      integer, intent(in) :: ldb
      complex(rsp), intent(in) :: bb(ldb, *)
      integer, intent(out) :: info
    end subroutine chegst

    subroutine zhegst(itype, uplo, nn, aa, lda, bb, ldb, info)
      import rdp
      integer, intent(in) :: itype
      character, intent(in) :: uplo
      integer, intent(in) :: nn
      integer, intent(in) :: lda
      complex(rdp), intent(inout) :: aa(lda, *)
      integer, intent(in) :: ldb
      complex(rdp), intent(in) :: bb(ldb, *)
      integer, intent(out) :: info
    end subroutine zhegst

    subroutine ssbgv(jobz, uplo, nn, ka, kb, ab, ldab, bb, ldbb, ww, zz, ldz,&
        & work, info)
      import rsp
      character, intent(in) :: jobz, uplo
      integer, intent(in) :: nn, ka, kb
      integer, intent(in) :: ldab
      real(rsp), intent(inout) :: ab(ldab, *)
      integer, intent(in) :: ldbb
      real(rsp), intent(inout) :: bb(ldbb, *)
      real(rsp), intent(out) :: ww(*)
      integer, intent(in) :: ldz
      real(rsp), intent(out) :: zz(ldz, *)
      real(rsp), intent(inout) :: work(*)
      integer, intent(out) :: info
    end subroutine ssbgv

    subroutine dsbgv(jobz, uplo, nn, ka, kb, ab, ldab, bb, ldbb, ww, zz, ldz,&
        & work, info)
      import rdp
      character, intent(in) :: jobz, uplo
      integer, intent(in) :: nn, ka, kb
      integer, intent(in) :: ldab
      real(rdp), intent(inout) :: ab(ldab, *)
      integer, intent(in) :: ldbb
      real(rdp), intent(inout) :: bb(ldbb, *)
      real(rdp), intent(out) :: ww(*)
      integer, intent(in) :: ldz
      real(rdp), intent(out) :: zz(ldz, *)
      real(rdp), intent(inout) :: work(*)
      integer, intent(out) :: info
    end subroutine dsbgv

    subroutine chbgv(jobz, uplo, nn, ka, kb, ab, ldab, bb, ldbb, ww, zz, ldz,&
        & work, rwork, info)
      import rsp
      character, intent(in) :: jobz, uplo
      integer, intent(in) :: nn, ka, kb
      integer, intent(in) :: ldab
      complex(rsp), intent(inout) :: ab(ldab, *)
      integer, intent(in) :: ldbb
      complex(rsp), intent(inout) :: bb(ldbb, *)
      real(rsp), intent(out) :: ww(*)
      integer, intent(in) :: ldz
      complex(rsp), intent(out) :: zz(ldz, *)
      complex(rsp), intent(inout) :: work(*)
      real(rsp), intent(inout) :: rwork(*)
      integer, intent(out) :: info
    end subroutine chbgv

    subroutine zhbgv(jobz, uplo, nn, ka, kb, ab, ldab, bb, ldbb, ww, zz, ldz,&
        & work, rwork, info)
      import rdp
      character, intent(in) :: jobz, uplo
      integer, intent(in) :: nn, ka, kb
      integer, intent(in) :: ldab
      complex(rdp), intent(inout) :: ab(ldab, *)
      integer, intent(in) :: ldbb
      complex(rdp), intent(inout) :: bb(ldbb, *)
      real(rdp), intent(out) :: ww(*)
      integer, intent(in) :: ldz
      complex(rdp), intent(out) :: zz(ldz, *)
      complex(rdp), intent(inout) :: work(*)
      real(rdp), intent(inout) :: rwork(*)
      integer, intent(out) :: info
    end subroutine zhbgv

    subroutine sgesv(nn, nrhs, aa, lda, ipiv, bb, ldb, info)
      import rsp
      integer, intent(in) :: nn, nrhs
      integer, intent(in) :: lda
      real(rsp), intent(inout) :: aa(lda, *)
      integer, intent(out) :: ipiv(*)
      integer, intent(in) :: ldb
      real(rsp), intent(inout) :: bb(ldb, *)
      integer, intent(out) :: info
    end subroutine sgesv

    subroutine dgesv(nn, nrhs, aa, lda, ipiv, bb, ldb, info)
      import rdp
      integer, intent(in) :: nn, nrhs
      integer, intent(in) :: lda
      real(rdp), intent(inout) :: aa(lda, *)
      integer, intent(out) :: ipiv(*)
      integer, intent(in) :: ldb
      real(rdp), intent(inout) :: bb(ldb, *)
      integer, intent(out) :: info
    end subroutine dgesv

    subroutine sgetrf(mm, nn, aa, lda, ipiv, info)
      import rsp
      integer, intent(in) :: mm, nn
      integer, intent(in) :: lda
      real(rsp), intent(inout) :: aa(lda, *)
      integer, intent(out) :: ipiv(*)
      integer, intent(out) :: info
    end subroutine sgetrf

    subroutine dgetrf(mm, nn, aa, lda, ipiv, info)
      import rdp
      integer, intent(in) :: mm, nn
      integer, intent(in) :: lda
      real(rdp), intent(inout) :: aa(lda, *)
      integer, intent(out) :: ipiv(*)
      integer, intent(out) :: info
    end subroutine dgetrf

    subroutine sgetri(nn, aa, lda, ipiv, work, lwork, info)
      import rsp
      integer, intent(in) :: nn
      integer, intent(in) :: lda
      real(rsp), intent(inout) :: aa(lda, *)
      integer, intent(in) :: ipiv(*)
      real(rsp), intent(inout) :: work(*)
      integer, intent(in) :: lwork
      integer, intent(out) :: info
    end subroutine sgetri

    subroutine dgetri(nn, aa, lda, ipiv, work, lwork, info)
      import rdp
      integer, intent(in) :: nn
      integer, intent(in) :: lda
      real(rdp), intent(inout) :: aa(lda, *)
      integer, intent(in) :: ipiv(*)
      real(rdp), intent(inout) :: work(*)
      integer, intent(in) :: lwork
      integer, intent(out) :: info
    end subroutine dgetri

    subroutine ssytrf(uplo, nn, aa, lda, ipiv, work, lwork, info)
      import rsp
      character, intent(in) :: uplo
      integer, intent(in) :: nn
      integer, intent(in) :: lda
      real(rsp), intent(inout) :: aa(lda, *)
      integer, intent(out) :: ipiv(*)
      real(rsp), intent(inout) :: work(*)
      integer, intent(in) :: lwork
      integer, intent(out) :: info
    end subroutine ssytrf

    subroutine dsytrf(uplo, nn, aa, lda, ipiv, work, lwork, info)
      import rdp
      character, intent(in) :: uplo
      integer, intent(in) :: nn
      integer, intent(in) :: lda
      real(rdp), intent(inout) :: aa(lda, *)
      integer, intent(out) :: ipiv(*)
      real(rdp), intent(inout) :: work(*)
      integer, intent(in) :: lwork
      integer, intent(out) :: info
    end subroutine dsytrf

    subroutine ssytrs(uplo, nn, nrhs, aa, lda, ipiv, bb, ldb, info)
      import rsp
      character, intent(in) :: uplo
      integer, intent(in) :: nn, nrhs
      integer, intent(in) :: lda
      real(rsp), intent(in) :: aa(lda, *)
      integer, intent(in) :: ipiv(*)
      integer, intent(in) :: ldb
      real(rsp), intent(inout) :: bb(ldb, *)
      integer, intent(out) :: info
    end subroutine ssytrs

    subroutine dsytrs(uplo, nn, nrhs, aa, lda, ipiv, bb, ldb, info)
      import rdp
      character, intent(in) :: uplo
      integer, intent(in) :: nn, nrhs
      integer, intent(in) :: lda
      real(rdp), intent(in) :: aa(lda, *)
      integer, intent(in) :: ipiv(*)
      integer, intent(in) :: ldb
      real(rdp), intent(inout) :: bb(ldb, *)
      integer, intent(out) :: info
    end subroutine dsytrs

    subroutine slarnv(idist, iseed, nn, xx)
      import rsp
      integer, intent(in) :: idist
      integer, intent(inout) :: iseed(4)
      integer, intent(in) :: nn
      real(rsp), intent(out) :: xx(*)
    end subroutine slarnv

    subroutine dlarnv(idist, iseed, nn, xx)
      import rdp
      integer, intent(in) :: idist
      integer, intent(inout) :: iseed(4)
      integer, intent(in) :: nn
      real(rdp), intent(out) :: xx(*)
    end subroutine dlarnv

    subroutine clarnv(idist, iseed, nn, xx)
      import rsp
      integer, intent(in) :: idist
      integer, intent(inout) :: iseed(4)
      integer, intent(in) :: nn
      complex(rsp), intent(out) :: xx(*)
    end subroutine clarnv

    subroutine zlarnv(idist, iseed, nn, xx)
      import rdp
      integer, intent(in) :: idist
      integer, intent(inout) :: iseed(4)
      integer, intent(in) :: nn
      complex(rdp), intent(out) :: xx(*)
    end subroutine zlarnv

    function ilaenv(ispec, name, opts, n1, n2, n3, n4)
      character, intent(in) :: name, opts
      integer, intent(in) :: ispec, n1, n2, n3, n4
      integer :: ilaenv
    end function ilaenv

    function slamch(cmach)
      import rsp
      character, intent(in) :: cmach
      real(rsp) :: slamch
    end function slamch

    function dlamch(cmach)
      import rdp
      character, intent(in) :: cmach
      real(rdp) :: dlamch
    end function dlamch

    subroutine xerbla(srname, info)
      character(6), intent(in) :: srname
      integer, intent(in) :: info
    end subroutine xerbla


  end interface

end module lapack

