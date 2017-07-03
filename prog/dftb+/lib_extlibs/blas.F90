!--------------------------------------------------------------------------------------------------!
!  DFTB+: general package for performing fast atomistic simulations                                !
!  Copyright (C) 2017  DFTB+ developers group                                                      !
!                                                                                                  !
!  See the LICENSE file for terms of usage and distribution.                                       !
!--------------------------------------------------------------------------------------------------!

!> Interface wrapper for the blas routines.
!!
!! All BLAS routines must be included here, which are called from the main code.
module blas
  use accuracy, only : rsp, rdp
  public

  interface

    subroutine ssyr(uplo, nn, alpha, xx, incx, aa, lda)
      import rsp
      character, intent(in) :: uplo
      integer, intent(in) :: nn
      real(rsp), intent(in) :: alpha
      real(rsp), intent(in) :: xx(*)
      integer, intent(in) :: incx
      integer, intent(in) :: lda
      real(rsp), intent(inout) :: aa(lda, *)
    end subroutine ssyr

    subroutine dsyr(uplo, nn, alpha, xx, incx, aa, lda)
      import rdp
      character, intent(in) :: uplo
      integer, intent(in) :: nn
      real(rdp), intent(in) :: alpha
      real(rdp), intent(in) :: xx(*)
      integer, intent(in) :: incx
      integer, intent(in) :: lda
      real(rdp), intent(inout) :: aa(lda, *)
    end subroutine dsyr

    subroutine cher(uplo, nn, alpha, xx, incx, aa, lda)
      import rsp
      character, intent(in) :: uplo
      integer, intent(in) :: nn
      real(rsp), intent(in) :: alpha
      complex(rsp), intent(in) :: xx(*)
      integer, intent(in) :: incx
      integer, intent(in) :: lda
      complex(rsp), intent(inout) :: aa(lda, *)
    end subroutine cher

    subroutine zher(uplo, nn, alpha, xx, incx, aa, lda)
      import rdp
      character, intent(in) :: uplo
      integer, intent(in) :: nn
      real(rdp), intent(in) :: alpha
      complex(rdp), intent(in) :: xx(*)
      integer, intent(in) :: incx
      integer, intent(in) :: lda
      complex(rdp), intent(inout) :: aa(lda, *)
    end subroutine zher

    subroutine sger(mm, nn, alpha, xx, incx, yy, incy, aa, lda)
      import rsp
      integer, intent(in) :: mm, nn
      real(rsp), intent(in) :: alpha
      real(rsp), intent(in) :: xx(*)
      integer, intent(in) :: incx
      real(rsp), intent(in) :: yy(*)
      integer, intent(in) :: incy
      integer, intent(in) :: lda
      real(rsp), intent(inout) :: aa(lda, *)
    end subroutine sger

    subroutine dger(mm, nn, alpha, xx, incx, yy, incy, aa, lda)
      import rdp
      integer, intent(in) :: mm, nn
      real(rdp), intent(in) :: alpha
      real(rdp), intent(in) :: xx(*)
      integer, intent(in) :: incx
      real(rdp), intent(in) :: yy(*)
      integer, intent(in) :: incy
      integer, intent(in) :: lda
      real(rdp), intent(inout) :: aa(lda, *)
    end subroutine dger

    subroutine cgerc(mm, nn, alpha, xx, incx, yy, incy, aa, lda)
      import rsp
      integer, intent(in) :: mm, nn
      complex(rsp), intent(in) :: alpha
      complex(rsp), intent(in) :: xx(*)
      integer, intent(in) :: incx
      complex(rsp), intent(in) :: yy(*)
      integer, intent(in) :: incy
      integer, intent(in) :: lda
      complex(rsp), intent(inout) :: aa(lda, *)
    end subroutine cgerc

    subroutine zgerc(mm, nn, alpha, xx, incx, yy, incy, aa, lda)
      import rdp
      integer, intent(in) :: mm, nn
      complex(rdp), intent(in) :: alpha
      complex(rdp), intent(in) :: xx(*)
      integer, intent(in) :: incx
      complex(rdp), intent(in) :: yy(*)
      integer, intent(in) :: incy
      integer, intent(in) :: lda
      complex(rdp), intent(inout) :: aa(lda, *)
    end subroutine zgerc

    subroutine ssymv(uplo, nn, alpha, aa, lda, xx, incx, beta, yy, incy)
      import rsp
      character, intent(in) :: uplo
      integer, intent(in) :: nn
      real(rsp), intent(in) :: alpha
      integer, intent(in) :: lda
      real(rsp), intent(in) :: aa(lda, *)
      real(rsp), intent(in) :: xx(*)
      integer, intent(in) :: incx
      real(rsp), intent(in) :: beta
      real(rsp), intent(inout) :: yy(*)
      integer, intent(in) :: incy
    end subroutine ssymv

    subroutine dsymv(uplo, nn, alpha, aa, lda, xx, incx, beta, yy, incy)
      import rdp
      character, intent(in) :: uplo
      integer, intent(in) :: nn
      real(rdp), intent(in) :: alpha
      integer, intent(in) :: lda
      real(rdp), intent(in) :: aa(lda, *)
      real(rdp), intent(in) :: xx(*)
      integer, intent(in) :: incx
      real(rdp), intent(in) :: beta
      real(rdp), intent(inout) :: yy(*)
      integer, intent(in) :: incy
    end subroutine dsymv

    subroutine chemv(uplo, nn, alpha, aa, lda, xx, incx, beta, yy, incy)
      import rsp
      character, intent(in) :: uplo
      integer, intent(in) :: nn
      complex(rsp), intent(in) :: alpha
      integer, intent(in) :: lda
      complex(rsp), intent(in) :: aa(lda, *)
      complex(rsp), intent(in) :: xx(*)
      integer, intent(in) :: incx
      complex(rsp), intent(in) :: beta
      complex(rsp), intent(inout) :: yy(*)
      integer, intent(in) :: incy
    end subroutine chemv

    subroutine zhemv(uplo, nn, alpha, aa, lda, xx, incx, beta, yy, incy)
      import rdp
      character, intent(in) :: uplo
      integer, intent(in) :: nn
      complex(rdp), intent(in) :: alpha
      integer, intent(in) :: lda
      complex(rdp), intent(in) :: aa(lda, *)
      complex(rdp), intent(in) :: xx(*)
      integer, intent(in) :: incx
      complex(rdp), intent(in) :: beta
      complex(rdp), intent(inout) :: yy(*)
      integer, intent(in) :: incy
    end subroutine zhemv

    subroutine sgemv(trans, mm, nn, alpha, aa, lda, xx, incx, beta, yy, incy)
      import rsp
      character, intent(in) :: trans
      integer, intent(in) :: mm, nn
      real(rsp), intent(in) :: alpha
      integer, intent(in) :: lda
      real(rsp), intent(in) :: aa(lda, *)
      real(rsp), intent(in) :: xx(*)
      integer, intent(in) :: incx
      real(rsp), intent(in) :: beta
      real(rsp), intent(inout) :: yy(*)
      integer, intent(in) :: incy
    end subroutine sgemv

    subroutine dgemv(trans, mm, nn, alpha, aa, lda, xx, incx, beta, yy, incy)
      import rdp
      character, intent(in) :: trans
      integer, intent(in) :: mm, nn
      real(rdp), intent(in) :: alpha
      integer, intent(in) :: lda
      real(rdp), intent(in) :: aa(lda, *)
      real(rdp), intent(in) :: xx(*)
      integer, intent(in) :: incx
      real(rdp), intent(in) :: beta
      real(rdp), intent(inout) :: yy(*)
      integer, intent(in) :: incy
    end subroutine dgemv

    subroutine ssbmv(uplo, nn, kk, alpha, aa, lda, xx, incx, beta, yy, incy)
      import rsp
      character, intent(in) :: uplo
      integer, intent(in) :: nn, kk
      real(rsp), intent(in) :: alpha
      integer, intent(in) :: lda
      real(rsp), intent(in) :: aa(lda, *)
      real(rsp), intent(in) :: xx(*)
      integer, intent(in) :: incx
      real(rsp), intent(in) :: beta
      real(rsp), intent(inout) :: yy(*)
      integer, intent(in) :: incy
    end subroutine ssbmv

    subroutine dsbmv(uplo, nn, kk, alpha, aa, lda, xx, incx, beta, yy, incy)
      import rdp
      character, intent(in) :: uplo
      integer, intent(in) :: nn, kk
      real(rdp), intent(in) :: alpha
      integer, intent(in) :: lda
      real(rdp), intent(in) :: aa(lda, *)
      real(rdp), intent(in) :: xx(*)
      integer, intent(in) :: incx
      real(rdp), intent(in) :: beta
      real(rdp), intent(inout) :: yy(*)
      integer, intent(in) :: incy
    end subroutine dsbmv

    subroutine ssymm(side, uplo, mm, nn, alpha, aa, lda, bb, ldb, beta, cc, ldc)
      import rsp
      character, intent(in) :: side, uplo
      integer, intent(in) :: mm, nn
      real(rsp), intent(in) :: alpha
      integer, intent(in) :: lda
      real(rsp), intent(in) :: aa(lda, *)
      integer, intent(in) :: ldb
      real(rsp), intent(in) :: bb(ldb, *)
      real(rsp), intent(in) :: beta
      integer, intent(in) :: ldc
      real(rsp), intent(inout) :: cc(ldc, *)
    end subroutine ssymm

    subroutine dsymm(side, uplo, mm, nn, alpha, aa, lda, bb, ldb, beta, cc, ldc)
      import rdp
      character, intent(in) :: side, uplo
      integer, intent(in) :: mm, nn
      real(rdp), intent(in) :: alpha
      integer, intent(in) :: lda
      real(rdp), intent(in) :: aa(lda, *)
      integer, intent(in) :: ldb
      real(rdp), intent(in) :: bb(ldb, *)
      real(rdp), intent(in) :: beta
      integer, intent(in) :: ldc
      real(rdp), intent(inout) :: cc(ldc, *)
    end subroutine dsymm

    subroutine chemm(side, uplo, mm, nn, alpha, aa, lda, bb, ldb, beta, cc, ldc)
      import rsp
      character, intent(in) :: side, uplo
      integer, intent(in) :: mm, nn
      complex(rsp), intent(in) :: alpha
      integer, intent(in) :: lda
      complex(rsp), intent(in) :: aa(lda, *)
      integer, intent(in) :: ldb
      complex(rsp), intent(in) :: bb(ldb, *)
      complex(rsp), intent(in) :: beta
      integer, intent(in) :: ldc
      complex(rsp), intent(inout) :: cc(ldc, *)
    end subroutine chemm

    subroutine zhemm(side, uplo, mm, nn, alpha, aa, lda, bb, ldb, beta, cc, ldc)
      import rdp
      character, intent(in) :: side, uplo
      integer, intent(in) :: mm, nn
      complex(rdp), intent(in) :: alpha
      integer, intent(in) :: lda
      complex(rdp), intent(in) :: aa(lda, *)
      integer, intent(in) :: ldb
      complex(rdp), intent(in) :: bb(ldb, *)
      complex(rdp), intent(in) :: beta
      integer, intent(in) :: ldc
      complex(rdp), intent(inout) :: cc(ldc, *)
    end subroutine zhemm

    subroutine sgemm(transa, transb, mm, nn, kk, alpha, aa, lda, bb, ldb, beta,&
        & cc, ldc)
      import rsp
      character, intent(in) :: transa, transb
      integer, intent(in) :: mm, nn, kk
      real(rsp), intent(in) :: alpha
      integer, intent(in) :: lda
      real(rsp), intent(in) :: aa(lda, *)
      integer, intent(in) :: ldb
      real(rsp), intent(in) :: bb(ldb, *)
      real(rsp), intent(in) :: beta
      integer, intent(in) :: ldc
      real(rsp), intent(inout) :: cc(ldc, *)
    end subroutine sgemm

    subroutine dgemm(transa, transb, mm, nn, kk, alpha, aa, lda, bb, ldb, beta,&
        & cc, ldc)
      import rdp
      character, intent(in) :: transa, transb
      integer, intent(in) :: mm, nn, kk
      real(rdp), intent(in) :: alpha
      integer, intent(in) :: lda
      real(rdp), intent(in) :: aa(lda, *)
      integer, intent(in) :: ldb
      real(rdp), intent(in) :: bb(ldb, *)
      real(rdp), intent(in) :: beta
      integer, intent(in) :: ldc
      real(rdp), intent(inout) :: cc(ldc, *)
    end subroutine dgemm

    subroutine cgemm(transa, transb, mm, nn, kk, alpha, aa, lda, bb, ldb, beta,&
        & cc, ldc)
      import rsp
      character, intent(in) :: transa, transb
      integer, intent(in) :: mm, nn, kk
      complex(rsp), intent(in) :: alpha
      integer, intent(in) :: lda
      complex(rsp), intent(in) :: aa(lda, *)
      integer, intent(in) :: ldb
      complex(rsp), intent(in) :: bb(ldb, *)
      complex(rsp), intent(in) :: beta
      integer, intent(in) :: ldc
      complex(rsp), intent(inout) :: cc(ldc, *)
    end subroutine cgemm

    subroutine zgemm(transa, transb, mm, nn, kk, alpha, aa, lda, bb, ldb, beta,&
        & cc, ldc)
      import rdp
      character, intent(in) :: transa, transb
      integer, intent(in) :: mm, nn, kk
      complex(rdp), intent(in) :: alpha
      integer, intent(in) :: lda
      complex(rdp), intent(in) :: aa(lda, *)
      integer, intent(in) :: ldb
      complex(rdp), intent(in) :: bb(ldb, *)
      complex(rdp), intent(in) :: beta
      integer, intent(in) :: ldc
      complex(rdp), intent(inout) :: cc(ldc, *)
    end subroutine zgemm

    subroutine ssyrk(uplo, trans, nn, kk, alpha, aa, lda, beta, cc, ldc)
      import rsp
      character, intent(in) :: uplo, trans
      integer, intent(in) :: nn, kk
      real(rsp), intent(in) :: alpha
      integer, intent(in) :: lda
      real(rsp), intent(in) :: aa(lda, *)
      real(rsp), intent(in) :: beta
      integer, intent(in) :: ldc
      real(rsp), intent(inout) :: cc(ldc, *)
    end subroutine ssyrk

    subroutine dsyrk(uplo, trans, nn, kk, alpha, aa, lda, beta, cc, ldc)
      import rdp
      character, intent(in) :: uplo, trans
      integer, intent(in) :: nn, kk
      real(rdp), intent(in) :: alpha
      integer, intent(in) :: lda
      real(rdp), intent(in) :: aa(lda, *)
      real(rdp), intent(in) :: beta
      integer, intent(in) :: ldc
      real(rdp), intent(inout) :: cc(ldc, *)
    end subroutine dsyrk

    subroutine cherk(uplo, trans, nn, kk, alpha, aa, lda, beta, cc, ldc)
      import rsp
      character, intent(in) :: uplo, trans
      integer, intent(in) :: nn, kk
      real(rsp), intent(in) :: alpha
      integer, intent(in) :: lda
      complex(rsp), intent(in) :: aa(lda, *)
      real(rsp), intent(in) :: beta
      integer, intent(in) :: ldc
      complex(rsp), intent(inout) :: cc(ldc, *)
    end subroutine cherk

    subroutine zherk(uplo, trans, nn, kk, alpha, aa, lda, beta, cc, ldc)
      import rdp
      character, intent(in) :: uplo, trans
      integer, intent(in) :: nn, kk
      real(rdp), intent(in) :: alpha
      integer, intent(in) :: lda
      complex(rdp), intent(in) :: aa(lda, *)
      real(rdp), intent(in) :: beta
      integer, intent(in) :: ldc
      complex(rdp), intent(inout) :: cc(ldc, *)
    end subroutine zherk

    subroutine strsm(side, uplo, transa, diag, mm, nn, alpha, aa, lda, bb, ldb)
      import rsp
      character, intent(in) :: side, uplo, transa, diag
      integer, intent(in) :: mm, nn
      real(rsp), intent(in) :: alpha
      integer, intent(in) :: lda
      real(rsp), intent(in) :: aa(lda, *)
      integer, intent(in) :: ldb
      real(rsp), intent(inout) :: bb(ldb, *)
    end subroutine strsm

    subroutine dtrsm(side, uplo, transa, diag, mm, nn, alpha, aa, lda, bb, ldb)
      import rdp
      character, intent(in) :: side, uplo, transa, diag
      integer, intent(in) :: mm, nn
      real(rdp), intent(in) :: alpha
      integer, intent(in) :: lda
      real(rdp), intent(in) :: aa(lda, *)
      integer, intent(in) :: ldb
      real(rdp), intent(inout) :: bb(ldb, *)
    end subroutine dtrsm

    subroutine ctrsm(side, uplo, transa, diag, mm, nn, alpha, aa, lda, bb, ldb)
      import rsp
      character, intent(in) :: side, uplo, transa, diag
      integer, intent(in) :: mm, nn
      complex(rsp), intent(in) :: alpha
      integer, intent(in) :: lda
      complex(rsp), intent(in) :: aa(lda, *)
      integer, intent(in) :: ldb
      complex(rsp), intent(inout) :: bb(ldb, *)
    end subroutine ctrsm

    subroutine ztrsm(side, uplo, transa, diag, mm, nn, alpha, aa, lda, bb, ldb)
      import rdp
      character, intent(in) :: side, uplo, transa, diag
      integer, intent(in) :: mm, nn
      complex(rdp), intent(in) :: alpha
      integer, intent(in) :: lda
      complex(rdp), intent(in) :: aa(lda, *)
      integer, intent(in) :: ldb
      complex(rdp), intent(inout) :: bb(ldb, *)
    end subroutine ztrsm

    subroutine strmm(side, uplo, transa, diag, mm, nn, alpha, aa, lda, bb, ldb)
      import rsp
      character, intent(in) :: side, uplo, transa, diag
      integer, intent(in) :: mm, nn
      real(rsp), intent(in) :: alpha
      integer, intent(in) :: lda
      real(rsp), intent(in) :: aa(lda, *)
      integer, intent(in) :: ldb
      real(rsp), intent(inout) :: bb(ldb, *)
    end subroutine strmm

    subroutine dtrmm(side, uplo, transa, diag, mm, nn, alpha, aa, lda, bb, ldb)
      import rdp
      character, intent(in) :: side, uplo, transa, diag
      integer, intent(in) :: mm, nn
      real(rdp), intent(in) :: alpha
      integer, intent(in) :: lda
      real(rdp), intent(in) :: aa(lda, *)
      integer, intent(in) :: ldb
      real(rdp), intent(inout) :: bb(ldb, *)
    end subroutine dtrmm

    subroutine ctrmm(side, uplo, transa, diag, mm, nn, alpha, aa, lda, bb, ldb)
      import rsp
      character, intent(in) :: side, uplo, transa, diag
      integer, intent(in) :: mm, nn
      complex(rsp), intent(in) :: alpha
      integer, intent(in) :: lda
      complex(rsp), intent(in) :: aa(lda, *)
      integer, intent(in) :: ldb
      complex(rsp), intent(inout) :: bb(ldb, *)
    end subroutine ctrmm

    subroutine ztrmm(side, uplo, transa, diag, mm, nn, alpha, aa, lda, bb, ldb)
      import rdp
      character, intent(in) :: side, uplo, transa, diag
      integer, intent(in) :: mm, nn
      complex(rdp), intent(in) :: alpha
      integer, intent(in) :: lda
      complex(rdp), intent(in) :: aa(lda, *)
      integer, intent(in) :: ldb
      complex(rdp), intent(inout) :: bb(ldb, *)
    end subroutine ztrmm

  end interface

end module blas
