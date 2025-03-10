!--------------------------------------------------------------------------------------------------!
!  DFTB+: general package for performing fast atomistic simulations                                !
!  Copyright (C) 2006 - 2025  DFTB+ developers group                                               !
!                                                                                                  !
!  See the LICENSE file for terms of usage and distribution.                                       !
!--------------------------------------------------------------------------------------------------!

!> Interface wrapper for the lapack routines. See the <a href="http://www.netlib.org/lapack/">lapack
!> project documentation</a> for more details
module dftbp_extlibs_lapack
  use dftbp_common_accuracy, only : rsp, rdp
  implicit none
  public


  !> Real symmetric eigensolver
  interface ssyev

    !> Real symmetric eigensolver
    subroutine ssyev(jobz, uplo, nn, aa, lda, ww, work, lwork, info)
      import rsp

      !> Job type
      character, intent(in) :: jobz

      !> Upper 'U' or lower 'L' triangle
      character, intent(in) :: uplo

      !> Matrix dimension
      integer, intent(in) :: nn

      !> Leading dimension of A
      integer, intent(in) :: lda

      !> Matrix A
      real(rsp), intent(inout) :: aa(lda, *)

      !> Eigenvalues
      real(rsp), intent(out) :: ww(*)

      !> Workspace
      real(rsp), intent(inout) :: work(*)

      !> Workspace sizing
      integer, intent(in) :: lwork

      !> State of routine on return
      integer, intent(out) :: info
    end subroutine ssyev
  end interface ssyev


  !> Double precision symmetric eigensolver
  interface dsyev

    !> Double precision symmetric eigensolver
    subroutine dsyev(jobz, uplo, nn, aa, lda, ww, work, lwork, info)
      import rdp

      !> Job type
      character, intent(in) :: jobz

      !> Upper 'U' or lower 'L' triangle
      character, intent(in) :: uplo

      !> Matrix dimension
      integer, intent(in) :: nn

      !> Leading dimension of A
      integer, intent(in) :: lda

      !> Matrix A
      real(rdp), intent(inout) :: aa(lda, *)

      !> Eigenvalues
      real(rdp), intent(out) :: ww(*)

      !> Workspace
      real(rdp), intent(inout) :: work(*)

      !> Workspace sizing
      integer, intent(in) :: lwork

      !> State of routine on return
      integer, intent(out) :: info
    end subroutine dsyev
  end interface dsyev


  !> Complex hermitian eigensolver
  interface cheev

    !> Complex hermitian eigensolver
    subroutine cheev(jobz, uplo, nn, aa, lda, ww, work, lwork, rwork, info)
      import rsp

      !> Job type
      character, intent(in) :: jobz

      !> Upper 'U' or lower 'L' triangle
      character, intent(in) :: uplo

      !> Matrix dimension
      integer, intent(in) :: nn

      !> Leading dimension of A
      integer, intent(in) :: lda

      !> Matrix A
      complex(rsp), intent(inout) :: aa(lda, *)

      !> Eigenvalues
      real(rsp), intent(out) :: ww(*)

      !> Workspace
      complex(rsp), intent(inout) :: work(*)

      !> Workspace sizing
      integer, intent(in) :: lwork

      !> Real workspace
      real(rsp), intent(inout) :: rwork(*)

      !> State of routine on return
      integer, intent(out) :: info
    end subroutine cheev
  end interface cheev


  !> Double complex hermitian eigensolver
  interface zheev

    !> Double complex hermitian eigensolver
    subroutine zheev(jobz, uplo, nn, aa, lda, ww, work, lwork, rwork, info)
      import rdp

      !> Job type
      character, intent(in) :: jobz

      !> Upper 'U' or lower 'L' triangle
      character, intent(in) :: uplo

      !> Matrix dimension
      integer, intent(in) :: nn

      !> Leading dimension of A
      integer, intent(in) :: lda

      !> Matrix A
      complex(rdp), intent(inout) :: aa(lda, *)

      !> Eigenvalues
      real(rdp), intent(out) :: ww(*)

      !> Workspace
      complex(rdp), intent(inout) :: work(*)

      !> Workspace sizing
      integer, intent(in) :: lwork

      !> Real workspace
      real(rdp), intent(inout) :: rwork(*)

      !> State of routine on return
      integer, intent(out) :: info
    end subroutine zheev
  end interface zheev


  !> Real symmetric eigensolver (divide and conquer)
  interface ssyevd

    !> Real symmetric eigensolver, divide and conquer
    subroutine ssyevd(jobz, uplo, nn, aa, lda, ww, work, lwork, iwork, liwork, info)
      import rsp

      !> Job type
      character, intent(in) :: jobz

      !> Upper 'U' or lower 'L' triangle
      character, intent(in) :: uplo

      !> Matrix dimension
      integer, intent(in) :: nn

      !> Leading dimension of A
      integer, intent(in) :: lda

      !> Matrix A
      real(rsp), intent(inout) :: aa(lda, *)

      !> Eigenvalues
      real(rsp), intent(out) :: ww(*)

      !> Workspace
      real(rsp), intent(inout) :: work(*)

      !> Workspace sizing
      integer, intent(in) :: lwork

      !> Integer workspace
      integer, intent(inout) :: iwork(*)

      !> Integer workspace sizing
      integer, intent(in) :: liwork

      !> State of routine on return
      integer, intent(out) :: info
    end subroutine ssyevd
  end interface ssyevd


  !> Double precision symmetric eigensolver (divide and conquer)
  interface dsyevd

    !> Real symmetric eigensolver, divide and conquer
    subroutine dsyevd(jobz, uplo, nn, aa, lda, ww, work, lwork, iwork, liwork, info)
      import rdp

      !> Job type
      character, intent(in) :: jobz

      !> Upper 'U' or lower 'L' triangle
      character, intent(in) :: uplo

      !> Matrix dimension
      integer, intent(in) :: nn

      !> Leading dimension of A
      integer, intent(in) :: lda

      !> Matrix A
      real(rdp), intent(inout) :: aa(lda, *)

      !> Eigenvalues
      real(rdp), intent(out) :: ww(*)

      !> Workspace
      real(rdp), intent(inout) :: work(*)

      !> Workspace sizing
      integer, intent(in) :: lwork

      !> Integer workspace
      integer, intent(inout) :: iwork(*)

      !> Integer workspace sizing
      integer, intent(in) :: liwork

      !> State of routine on return
      integer, intent(out) :: info
    end subroutine dsyevd
  end interface dsyevd


  !> Complex hermitian eigensolver (divide and conquer)
  interface cheevd

    !> Complex hermitian eigensolver
    subroutine cheevd(jobz, uplo, nn, aa, lda, ww, work, lwork, rwork, lrwork, iwork, liwork, info)
      import rsp

      !> Job type
      character, intent(in) :: jobz

      !> Upper 'U' or lower 'L' triangle
      character, intent(in) :: uplo

      !> Matrix dimension
      integer, intent(in) :: nn

      !> Leading dimension of A
      integer, intent(in) :: lda

      !> Matrix A
      complex(rsp), intent(inout) :: aa(lda, *)

      !> Eigenvalues
      real(rsp), intent(out) :: ww(*)

      !> Workspace
      complex(rsp), intent(inout) :: work(*)

      !> Workspace sizing
      integer, intent(in) :: lwork

      !> Real workspace
      real(rsp), intent(inout) :: rwork(*)

      !> Size of rwork
      integer, intent(in) :: lrwork

      !> Integer workspace
      integer, intent(inout) :: iwork(*)

      !> Integer workspace sizing
      integer, intent(in) :: liwork

      !> State of routine on return
      integer, intent(out) :: info
    end subroutine cheevd
  end interface cheevd


  !> Double complex hermitian eigensolver (divide and conquer)
  interface zheevd

    !> Double complex hermitian eigensolver
    subroutine zheevd(jobz, uplo, nn, aa, lda, ww, work, lwork, rwork, lrwork, iwork, liwork, info)
      import rdp

      !> Job type
      character, intent(in) :: jobz

      !> Upper 'U' or lower 'L' triangle
      character, intent(in) :: uplo

      !> Matrix dimension
      integer, intent(in) :: nn

      !> Leading dimension of A
      integer, intent(in) :: lda

      !> Matrix A
      complex(rdp), intent(inout) :: aa(lda, *)

      !> Eigenvalues
      real(rdp), intent(out) :: ww(*)

      !> Workspace
      complex(rdp), intent(inout) :: work(*)

      !> Workspace sizing
      integer, intent(in) :: lwork

      !> Real workspace
      real(rdp), intent(inout) :: rwork(*)

      !> Size of rwork
      integer, intent(in) :: lrwork

      !> Integer workspace
      integer, intent(inout) :: iwork(*)

      !> Integer workspace sizing
      integer, intent(in) :: liwork

      !> State of routine on return
      integer, intent(out) :: info
    end subroutine zheevd
  end interface zheevd


  !> Real symmetric generalised eigensolver
  interface ssygv

    !> Real symmetric generalised eigensolver
    subroutine ssygv(itype, jobz, uplo, nn, aa, lda, bb, ldb, ww, work, lwork,&
        & info)
      import rsp

      !> Specifies the problem type to be solved
      integer, intent(in) :: itype

      !> Job type
      character, intent(in) :: jobz

      !> Upper 'U' or lower 'L' triangle
      character, intent(in) :: uplo

      !> Matrix dimension
      integer, intent(in) :: nn

      !> Leading dimension of A
      integer, intent(in) :: lda

      !> Matrix A
      real(rsp), intent(inout) :: aa(lda, *)

      !> Leading dimension of B
      integer, intent(in) :: ldb

      !> Matrix B
      real(rsp), intent(inout) :: bb(ldb, *)

      !> Eigenvalues
      real(rsp), intent(out) :: ww(*)

      !> Workspace
      real(rsp), intent(inout) :: work(*)

      !> Workspace sizing
      integer, intent(in) :: lwork

      !> State of routine on return
      integer, intent(out) :: info
    end subroutine ssygv
  end interface ssygv


  !> Double precision generalised symmetric eigensolver
  interface dsygv

    !> Double precision generalised symmetric eigensolver
    subroutine dsygv(itype, jobz, uplo, nn, aa, lda, bb, ldb, ww, work, lwork,&
        & info)
      import rdp

      !> Specifies the problem type to be solved
      integer, intent(in) :: itype

      !> Job type
      character, intent(in) :: jobz

      !> Upper 'U' or lower 'L' triangle
      character, intent(in) :: uplo

      !> Matrix dimension
      integer, intent(in) :: nn

      !> Leading dimension of A
      integer, intent(in) :: lda

      !> Matrix A
      real(rdp), intent(inout) :: aa(lda, *)

      !> Leading dimension of B
      integer, intent(in) :: ldb

      !> Matrix B
      real(rdp), intent(inout) :: bb(ldb, *)

      !> Eigenvalues
      real(rdp), intent(out) :: ww(*)

      !> Workspace
      real(rdp), intent(inout) :: work(*)

      !> Workspace sizing
      integer, intent(in) :: lwork

      !> State of routine on return
      integer, intent(out) :: info
    end subroutine dsygv
  end interface dsygv


  !> Complex generalised hermitian eigensolver
  interface chegv

    !> Complex generalised hermitian eigensolver
    subroutine chegv(itype, jobz, uplo, nn, aa, lda, bb, ldb, ww, work, lwork,&
        & rwork, info)
      import rsp

      !> Specifies the problem type to be solved
      integer, intent(in) :: itype

      !> Job type
      character, intent(in) :: jobz

      !> Upper 'U' or lower 'L' triangle
      character, intent(in) :: uplo

      !> Matrix dimension
      integer, intent(in) :: nn

      !> Leading dimension of A
      integer, intent(in) :: lda

      !> Matrix A
      complex(rsp), intent(inout) :: aa(lda, *)

      !> Leading dimension of B
      integer, intent(in) :: ldb

      !> Matrix B
      complex(rsp), intent(inout) :: bb(ldb, *)

      !> Eigenvalues
      real(rsp), intent(out) :: ww(*)

      !> Workspace
      complex(rsp), intent(inout) :: work(*)

      !> Workspace sizing
      integer, intent(in) :: lwork

      !> Real workspace
      real(rsp), intent(inout) :: rwork(*)

      !> State of routine on return
      integer, intent(out) :: info
    end subroutine chegv
  end interface chegv


  !> Double complex generalised hermitian eigensolver
  interface zhegv

    !> Double complex generalised hermitian eigensolver
    subroutine zhegv(itype, jobz, uplo, nn, aa, lda, bb, ldb, ww, work, lwork,&
        & rwork, info)
      import rdp

      !> Specifies the problem type to be solved
      integer, intent(in) :: itype

      !> Job type
      character, intent(in) :: jobz

      !> Upper 'U' or lower 'L' triangle
      character, intent(in) :: uplo

      !> Matrix dimension
      integer, intent(in) :: nn

      !> Leading dimension of A
      integer, intent(in) :: lda

      !> Matrix A
      complex(rdp), intent(inout) :: aa(lda, *)

      !> Leading dimension of B
      integer, intent(in) :: ldb

      !> Matrix B
      complex(rdp), intent(inout) :: bb(ldb, *)

      !> Eigenvalues
      real(rdp), intent(out) :: ww(*)

      !> Workspace
      complex(rdp), intent(inout) :: work(*)

      !> Workspace sizing
      integer, intent(in) :: lwork

      !> Real workspace
      real(rdp), intent(inout) :: rwork(*)

      !> State of routine on return
      integer, intent(out) :: info
    end subroutine zhegv
  end interface zhegv


  !> Real symmetric generalised eigensolver, divide and conquer
  interface ssygvd

    !> Real symmetric generalised eigensolver, divide and conquer
    subroutine ssygvd(itype, jobz, uplo, nn, aa, lda, bb, ldb, ww, work,&
        & lwork, iwork, liwork, info)
      import rsp

      !> Specifies the problem type to be solved
      integer, intent(in) :: itype

      !> Job type
      character, intent(in) :: jobz

      !> Upper 'U' or lower 'L' triangle
      character, intent(in) :: uplo

      !> Matrix dimension
      integer, intent(in) :: nn

      !> Leading dimension of A
      integer, intent(in) :: lda

      !> Matrix A
      real(rsp), intent(inout) :: aa(lda, *)

      !> Leading dimension of B
      integer, intent(in) :: ldb

      !> Matrix B
      real(rsp), intent(inout) :: bb(ldb, *)

      !> Eigenvalues
      real(rsp), intent(out) :: ww(*)

      !> Workspace
      real(rsp), intent(inout) :: work(*)

      !> Workspace sizing
      integer, intent(in) :: lwork

      !> Integer workspace
      integer, intent(inout) :: iwork(*)

      !> Size of integer workspace
      integer, intent(in) :: liwork

      !> State of routine on return
      integer, intent(out) :: info
    end subroutine ssygvd
  end interface ssygvd


  !> Double precision generalised symmetric eigensolver, divide and conquer
  interface dsygvd

    !> Double precision generalised symmetric eigensolver, divide and conquer
    subroutine dsygvd(itype, jobz, uplo, nn, aa, lda, bb, ldb, ww, work,&
        & lwork, iwork, liwork, info)
      import rdp

      !> Specifies the problem type to be solved
      integer, intent(in) :: itype

      !> Job type
      character, intent(in) :: jobz

      !> Upper 'U' or lower 'L' triangle
      character, intent(in) :: uplo

      !> Matrix dimension
      integer, intent(in) :: nn

      !> Leading dimension of A
      integer, intent(in) :: lda

      !> Matrix A
      real(rdp), intent(inout) :: aa(lda, *)

      !> Leading dimension of B
      integer, intent(in) :: ldb

      !> Matrix B
      real(rdp), intent(inout) :: bb(ldb, *)

      !> Eigenvalues
      real(rdp), intent(out) :: ww(*)

      !> Workspace
      real(rdp), intent(inout) :: work(*)

      !> Workspace sizing
      integer, intent(in) :: lwork

      !> Integer workspace
      integer, intent(inout) :: iwork(*)

      !> Size of integer workspace
      integer, intent(in) :: liwork

      !> State of routine on return
      integer, intent(out) :: info
    end subroutine dsygvd
  end interface dsygvd


  !> Complex generalised hermitian eigensolver, divide and conquer
  interface chegvd

    !> Complex generalised hermitian eigensolver, divide and conquer
    subroutine chegvd(itype, jobz, uplo, nn, aa, lda, bb, ldb, ww, work,&
        & lwork, rwork, lrwork, iwork, liwork, info)
      import rsp

      !> Specifies the problem type to be solved
      integer, intent(in) :: itype

      !> Job type
      character, intent(in) :: jobz

      !> Upper 'U' or lower 'L' triangle
      character, intent(in) :: uplo

      !> Matrix dimension
      integer, intent(in) :: nn

      !> Leading dimension of A
      integer, intent(in) :: lda

      !> Matrix A
      complex(rsp), intent(inout) :: aa(lda, *)

      !> Leading dimension of B
      integer, intent(in) :: ldb

      !> Matrix B
      complex(rsp), intent(inout) :: bb(ldb, *)

      !> Eigenvalues
      real(rsp), intent(out) :: ww(*)

      !> Workspace
      complex(rsp), intent(inout) :: work(*)

      !> Workspace sizing
      integer, intent(in) :: lwork

      !> Real workspace
      real(rsp), intent(inout) :: rwork(*)

      !> Size of rwork
      integer, intent(in) :: lrwork

      !> Integer workspace
      integer, intent(inout) :: iwork(*)

      !> Size of integer workspace
      integer, intent(in) :: liwork

      !> State of routine on return
      integer, intent(out) :: info
    end subroutine chegvd
  end interface chegvd


  !> Double complex generalised hermitian eigensolver, divide and conquer
  interface zhegvd

    !> Double complex generalised hermitian eigensolver, divide and conquer
    subroutine zhegvd(itype, jobz, uplo, nn, aa, lda, bb, ldb, ww, work,&
        & lwork, rwork, lrwork, iwork, liwork, info)
      import rdp

      !> Specifies the problem type to be solved
      integer, intent(in) :: itype

      !> Job type
      character, intent(in) :: jobz

      !> Upper 'U' or lower 'L' triangle
      character, intent(in) :: uplo

      !> Matrix dimension
      integer, intent(in) :: nn

      !> Leading dimension of A
      integer, intent(in) :: lda

      !> Matrix A
      complex(rdp), intent(inout) :: aa(lda, *)

      !> Leading dimension of B
      integer, intent(in) :: ldb

      !> Matrix B
      complex(rdp), intent(inout) :: bb(ldb, *)

      !> Eigenvalues
      real(rdp), intent(out) :: ww(*)

      !> Workspace
      complex(rdp), intent(inout) :: work(*)

      !> Workspace sizing
      integer, intent(in) :: lwork

      !> Real workspace
      real(rdp), intent(inout) :: rwork(*)

      !> Workspace size for rwork
      integer, intent(in) :: lrwork

      !> Integer workspace
      integer, intent(inout) :: iwork(*)

      !> Size of integer workspace
      integer, intent(in) :: liwork

      !> State of routine on return
      integer, intent(out) :: info
    end subroutine zhegvd
  end interface zhegvd


  !> Real symmetric generalised eigensolver, relatively robust
  interface ssyevr

    !> Real symmetric generalised eigensolver, relatively robust
    subroutine ssyevr(jobz, range, uplo, nn, aa, lda, vl, vu, il, iu, abstol,&
        & mm, ww, zz, ldz, isuppz, work, lwork, iwork, liwork, info)
      import rsp

      !> Job type
      character, intent(in) :: jobz

      !> Choice for range of eigenstates, 'A'll, 'V' half range (VL,VU], 'I' IL-th through IU-th
      !> eigenvalues
      character, intent(in) :: range

      !> Upper 'U' or lower 'L' triangle
      character, intent(in) :: uplo

      !> Matrix dimension
      integer, intent(in) :: nn

      !> Leading dimension of A
      integer, intent(in) :: lda

      !> Matrix A
      real(rsp), intent(inout) :: aa(lda, *)

      !> Lower range if in mode range =  V
      real(rsp), intent(in) :: vl

      !> Upper range
      real(rsp), intent(in) :: vu

      !> Lower number if in mode range =  I
      integer, intent(in) :: il

      !> Upper number if in mode range =  I
      integer, intent(in) :: iu

      !> Absolute error tolerance for the eigenvalues
      real(rsp), intent(in) :: abstol

      !> Total number of eigenvalues found
      integer, intent(out) :: mm

      !> Eigenvalues
      real(rsp), intent(out) :: ww(*)

      !> Leading dimension of Z
      integer, intent(in) :: ldz

      !> Matrix Z
      real(rsp), intent(out) :: zz(ldz, *)

      !> Support of the eigenvectors in Z
      integer, intent(out) :: isuppz(*)

      !> Workspace
      real(rsp), intent(inout) :: work(*)

      !> Workspace sizing
      integer, intent(in) :: lwork

      !> Integer workspace
      integer, intent(inout) :: iwork(*)

      !> Size of integer workspace
      integer, intent(in) :: liwork

      !> State of routine on return
      integer, intent(out) :: info
    end subroutine ssyevr
  end interface ssyevr


  !> Double precision generalised symmetric eigensolver, relatively robust
  interface dsyevr

    !> Double precision generalised symmetric eigensolver, relatively robust
    subroutine dsyevr(jobz, range, uplo, nn, aa, lda, vl, vu, il, iu, abstol,&
        & mm, ww, zz, ldz, isuppz, work, lwork, iwork, liwork, info)
      import rdp

      !> Job type
      character, intent(in) :: jobz

      !> Choice for range of eigenstates, 'A'll, 'V' half range (VL,VU], 'I' IL-th through IU-th
      !> eigenvalues
      character, intent(in) :: range

      !> Upper 'U' or lower 'L' triangle
      character, intent(in) :: uplo

      !> Matrix dimension
      integer, intent(in) :: nn

      !> Leading dimension of A
      integer, intent(in) :: lda

      !> Matrix A
      real(rdp), intent(inout) :: aa(lda, *)

      !> Lower range if in mode range =  V
      real(rdp), intent(in) :: vl

      !> Upper range
      real(rdp), intent(in) :: vu

      !> Lower number if in mode range =  I
      integer, intent(in) :: il

      !> Upper number if in mode range =  I
      integer, intent(in) :: iu

      !> Absolute error tolerance for the eigenvalues
      real(rdp), intent(in) :: abstol

      !> Total number of eigenvalues found
      integer, intent(out) :: mm

      !> Eigenvalues
      real(rdp), intent(out) :: ww(*)

      !> Leading dimension of Z
      integer, intent(in) :: ldz

      !> Matrix Z
      real(rdp), intent(out) :: zz(ldz, *)

      !> Support of the eigenvectors in Z
      integer, intent(out) :: isuppz(*)

      !> Workspace
      real(rdp), intent(inout) :: work(*)

      !> Workspace sizing
      integer, intent(in) :: lwork

      !> Integer workspace
      integer, intent(inout) :: iwork(*)

      !> Size of integer workspace
      integer, intent(in) :: liwork

      !> State of routine on return
      integer, intent(out) :: info
    end subroutine dsyevr
  end interface dsyevr


  !> Complex generalised hermitian eigensolver, relatively robust
  interface cheevr

    !> Complex generalised hermitian eigensolver, relatively robust
    subroutine cheevr(jobz, range, uplo, nn, aa, lda, vl, vu, il, iu, abstol,&
        & mm, ww, zz, ldz, isuppz, work, lwork, rwork, lrwork, iwork, liwork,&
        & info)
      import rsp

      !> Job type
      character, intent(in) :: jobz

      !> Choice for range of eigenstates, 'A'll, 'V' half range (VL,VU], 'I' IL-th through IU-th
      !> eigenvalues
      character, intent(in) :: range

      !> Upper 'U' or lower 'L' triangle
      character, intent(in) :: uplo

      !> Matrix dimension
      integer, intent(in) :: nn

      !> Leading dimension of A
      integer, intent(in) :: lda

      !> Matrix A
      complex(rsp), intent(inout) :: aa(lda, *)

      !> Lower range if in mode range =  V
      real(rsp), intent(in) :: vl

      !> Upper range
      real(rsp), intent(in) :: vu

      !> Lower number if in mode range =  I
      integer, intent(in) :: il

      !> Upper number if in mode range =  I
      integer, intent(in) :: iu

      !> Absolute error tolerance for the eigenvalues
      real(rsp), intent(in) :: abstol

      !> Total number of eigenvalues found
      integer, intent(out) :: mm

      !> Eigenvalues
      real(rsp), intent(out) :: ww(*)

      !> Leading dimension of Z
      integer, intent(in) :: ldz

      !> Matrix Z
      complex(rsp), intent(out) :: zz(ldz, *)

      !> Support of the eigenvectors in Z
      integer, intent(out) :: isuppz(*)

      !> Workspace
      complex(rsp), intent(inout) :: work(*)

      !> Workspace sizing
      integer, intent(in) :: lwork

      !> Real workspace
      real(rsp), intent(inout) :: rwork(*)

      !> Size of rwork
      integer, intent(in) :: lrwork

      !> Integer workspace
      integer, intent(inout) :: iwork(*)

      !> Size of integer workspace
      integer, intent(in) :: liwork

      !> State of routine on return
      integer, intent(out) :: info
    end subroutine cheevr
  end interface cheevr


  !> Complex generalised hermitian eigensolver, relatively robust
  interface zheevr

    !> Complex generalised hermitian eigensolver, relatively robust
    subroutine zheevr(jobz, range, uplo, nn, aa, lda, vl, vu, il, iu, abstol,&
        & mm, ww, zz, ldz, isuppz, work, lwork, rwork, lrwork, iwork, liwork,&
        & info)
      import rdp

      !> Job type
      character, intent(in) :: jobz

      !> Choice for range of eigenstates, 'A'll, 'V' half range (VL,VU], 'I' IL-th through IU-th
      !> eigenvalues
      character, intent(in) :: range

      !> Upper 'U' or lower 'L' triangle
      character, intent(in) :: uplo

      !> Matrix dimension
      integer, intent(in) :: nn

      !> Leading dimension of A
      integer, intent(in) :: lda

      !> Matrix A
      complex(rdp), intent(inout) :: aa(lda, *)

      !> Lower range if in mode range =  V
      real(rdp), intent(in) :: vl

      !> Upper range
      real(rdp), intent(in) :: vu

      !> Lower number if in mode range =  I
      integer, intent(in) :: il

      !> Upper number if in mode range =  I
      integer, intent(in) :: iu

      !> Absolute error tolerance for the eigenvalues
      real(rdp), intent(in) :: abstol

      !> Total number of eigenvalues found
      integer, intent(out) :: mm

      !> Eigenvalues
      real(rdp), intent(out) :: ww(*)

      !> Leading dimension of Z
      integer, intent(in) :: ldz

      !> Matrix Z
      complex(rdp), intent(out) :: zz(ldz, *)

      !> Support of the eigenvectors in Z
      integer, intent(out) :: isuppz(*)

      !> Workspace
      complex(rdp), intent(inout) :: work(*)

      !> Workspace sizing
      integer, intent(in) :: lwork

      !> Real workspace
      real(rdp), intent(inout) :: rwork(*)

      !> Size of rwork
      integer, intent(in) :: lrwork

      !> Integer workspace
      integer, intent(inout) :: iwork(*)

      !> Size of integer workspace
      integer, intent(in) :: liwork

      !> State of routine on return
      integer, intent(out) :: info
    end subroutine zheevr
  end interface zheevr


  !> Cholesky factorization of real symmetric positive definite matrix
  interface spotrf

    !> Cholesky factorization of real symmetric positive definite matrix
    subroutine spotrf(uplo, nn, aa, lda, info)
      import rsp

      !> Upper 'U' or lower 'L' triangle
      character, intent(in) :: uplo

      !> Matrix dimension
      integer, intent(in) :: nn

      !> Leading dimension of A
      integer, intent(in) :: lda

      !> Matrix A
      real(rsp), intent(inout) :: aa(lda, *)

      !> State of routine on return
      integer, intent(out) :: info
    end subroutine spotrf
  end interface spotrf


  !> Cholesky factorization of double precision symmetric positive definite matrix
  interface dpotrf

    !> Cholesky factorization of double precision symmetric positive definite matrix
    subroutine dpotrf(uplo, nn, aa, lda, info)
      import rdp

      !> Upper 'U' or lower 'L' triangle
      character, intent(in) :: uplo

      !> Matrix dimension
      integer, intent(in) :: nn

      !> Leading dimension of A
      integer, intent(in) :: lda

      !> Matrix A
      real(rdp), intent(inout) :: aa(lda, *)

      !> State of routine on return
      integer, intent(out) :: info
    end subroutine dpotrf
  end interface dpotrf


  !> Cholesky factorization of complex hermitian positive definite matrix
  interface cpotrf

    !> Cholesky factorization of complex hermitian positive definite matrix
    subroutine cpotrf(uplo, nn, aa, lda, info)
      import rsp

      !> Upper 'U' or lower 'L' triangle
      character, intent(in) :: uplo

      !> Matrix dimension
      integer, intent(in) :: nn

      !> Leading dimension of A
      integer, intent(in) :: lda

      !> Matrix A
      complex(rsp), intent(inout) :: aa(lda, *)

      !> State of routine on return
      integer, intent(out) :: info
    end subroutine cpotrf
  end interface cpotrf


  !> Cholesky factorization of double complex hermitian positive definite matrix
  interface zpotrf

    !> Cholesky factorization of double complex hermitian positive definite matrix
    subroutine zpotrf(uplo, nn, aa, lda, info)
      import rdp

      !> Upper 'U' or lower 'L' triangle
      character, intent(in) :: uplo

      !> Matrix dimension
      integer, intent(in) :: nn

      !> Leading dimension of A
      integer, intent(in) :: lda

      !> Matrix A
      complex(rdp), intent(inout) :: aa(lda, *)

      !> State of routine on return
      integer, intent(out) :: info
    end subroutine zpotrf
  end interface zpotrf


  !> Reduce real symmetric-definite generalized eigenproblem to standard form
  interface ssygst

    !> Reduce real symmetric-definite generalized eigenproblem to standard form
    subroutine ssygst(itype, uplo, nn, aa, lda, bb, ldb, info)
      import rsp

      !> Specifies the problem type to be solved
      integer, intent(in) :: itype

      !> Upper 'U' or lower 'L' triangle
      character, intent(in) :: uplo

      !> Matrix dimension
      integer, intent(in) :: nn

      !> Leading dimension of A
      integer, intent(in) :: lda

      !> Matrix A
      real(rsp), intent(inout) :: aa(lda, *)

      !> Leading dimension of B
      integer, intent(in) :: ldb

      !> Matrix B
      real(rsp), intent(in) :: bb(ldb, *)

      !> State of routine on return
      integer, intent(out) :: info
    end subroutine ssygst
  end interface ssygst


  !> Reduce double precision symmetric-definite generalized eigenproblem to standard form
  interface dsygst

    !> Reduce double precision symmetric-definite generalized eigenproblem to standard form
    subroutine dsygst(itype, uplo, nn, aa, lda, bb, ldb, info)
      import rdp

      !> Specifies the problem type to be solved
      integer, intent(in) :: itype

      !> Upper 'U' or lower 'L' triangle
      character, intent(in) :: uplo

      !> Matrix dimension
      integer, intent(in) :: nn

      !> Leading dimension of A
      integer, intent(in) :: lda

      !> Matrix A
      real(rdp), intent(inout) :: aa(lda, *)

      !> Leading dimension of B
      integer, intent(in) :: ldb

      !> Matrix B
      real(rdp), intent(in) :: bb(ldb, *)

      !> State of routine on return
      integer, intent(out) :: info
    end subroutine dsygst
  end interface dsygst


  !> Reduce complex hermitian-definite generalized eigenproblem to standard form
  interface chegst

    !> Reduce complex hermitian-definite generalized eigenproblem to standard form
    subroutine chegst(itype, uplo, nn, aa, lda, bb, ldb, info)
      import rsp

      !> Specifies the problem type to be solved
      integer, intent(in) :: itype

      !> Upper 'U' or lower 'L' triangle
      character, intent(in) :: uplo

      !> Matrix dimension
      integer, intent(in) :: nn

      !> Leading dimension of A
      integer, intent(in) :: lda

      !> Matrix A
      complex(rsp), intent(inout) :: aa(lda, *)

      !> Leading dimension of B
      integer, intent(in) :: ldb

      !> Matrix B
      complex(rsp), intent(in) :: bb(ldb, *)

      !> State of routine on return
      integer, intent(out) :: info
    end subroutine chegst
  end interface chegst


  !> Reduce double complex hermitian-definite generalized eigenproblem to standard form
  interface zhegst

    !> Reduce double complex hermitian-definite generalized eigenproblem to standard form
    subroutine zhegst(itype, uplo, nn, aa, lda, bb, ldb, info)
      import rdp

      !> Specifies the problem type to be solved
      integer, intent(in) :: itype

      !> Upper 'U' or lower 'L' triangle
      character, intent(in) :: uplo

      !> Matrix dimension
      integer, intent(in) :: nn

      !> Leading dimension of A
      integer, intent(in) :: lda

      !> Matrix A
      complex(rdp), intent(inout) :: aa(lda, *)

      !> Leading dimension of B
      integer, intent(in) :: ldb

      !> Matrix B
      complex(rdp), intent(in) :: bb(ldb, *)

      !> State of routine on return
      integer, intent(out) :: info
    end subroutine zhegst
  end interface zhegst


  !> Solve overdetermined or underdetermined real linear systems
  interface sgesv

    !> Solve overdetermined or underdetermined real linear systems
    subroutine sgesv(nn, nrhs, aa, lda, ipiv, bb, ldb, info)
      import rsp

      !> Matrix dimension
      integer, intent(in) :: nn

      !> Number of right hand side equations
      integer, intent(in) :: nrhs

      !> Leading dimension of A
      integer, intent(in) :: lda

      !> Matrix A
      real(rsp), intent(inout) :: aa(lda, *)

      !> Pivot array
      integer, intent(out) :: ipiv(*)

      !> Leading dimension of B
      integer, intent(in) :: ldb

      !> Matrix B
      real(rsp), intent(inout) :: bb(ldb, *)

      !> State of routine on return
      integer, intent(out) :: info
    end subroutine sgesv
  end interface sgesv


  !> Solve overdetermined or underdetermined double precision linear systems
  interface dgesv

    !> Solve overdetermined or underdetermined double precision linear systems
    subroutine dgesv(nn, nrhs, aa, lda, ipiv, bb, ldb, info)
      import rdp

      !> Matrix dimension
      integer, intent(in) :: nn

      !> Number of right hand side equations
      integer, intent(in) :: nrhs

      !> Leading dimension of A
      integer, intent(in) :: lda

      !> Matrix A
      real(rdp), intent(inout) :: aa(lda, *)

      !> Pivot array
      integer, intent(out) :: ipiv(*)

      !> Leading dimension of B
      integer, intent(in) :: ldb

      !> Matrix B
      real(rdp), intent(inout) :: bb(ldb, *)

      !> State of routine on return
      integer, intent(out) :: info
    end subroutine dgesv
  end interface dgesv


  !> Computes LU factorization of real matrix
  interface sgetrf

    !> Computes LU factorization of real matrix
    subroutine sgetrf(mm, nn, aa, lda, ipiv, info)
      import rsp

      !> Number of rows of the matrix
      integer, intent(in) :: mm

      !> Matrix dimension
      integer, intent(in) :: nn

      !> Leading dimension of A
      integer, intent(in) :: lda

      !> Matrix A
      real(rsp), intent(inout) :: aa(lda, *)

      !> Pivot array
      integer, intent(out) :: ipiv(*)

      !> State of routine on return
      integer, intent(out) :: info
    end subroutine sgetrf
  end interface sgetrf


  !> Computes LU factorization of double precision matrix
  interface dgetrf

    !> Computes LU factorization of double precision matrix
    subroutine dgetrf(mm, nn, aa, lda, ipiv, info)
      import rdp

      !> Number of rows of the matrix
      integer, intent(in) :: mm

      !> Matrix dimension
      integer, intent(in) :: nn

      !> Leading dimension of A
      integer, intent(in) :: lda

      !> Matrix A
      real(rdp), intent(inout) :: aa(lda, *)

      !> Pivot array
      integer, intent(out) :: ipiv(*)

      !> State of routine on return
      integer, intent(out) :: info
    end subroutine dgetrf
  end interface dgetrf


  !> Computes LU factorization of complex matrix
  interface cgetrf

    !> Computes LU factorization of real matrix
    subroutine cgetrf(mm, nn, aa, lda, ipiv, info)
      import rsp

      !> Number of rows of the matrix
      integer, intent(in) :: mm

      !> Matrix dimension
      integer, intent(in) :: nn

      !> Leading dimension of A
      integer, intent(in) :: lda

      !> Matrix A
      complex(rsp), intent(inout) :: aa(lda, *)

      !> Pivot array
      integer, intent(out) :: ipiv(*)

      !> State of routine on return
      integer, intent(out) :: info
    end subroutine cgetrf
  end interface cgetrf


  !> Computes LU factorization of double precision complex matrix
  interface zgetrf

    !> Computes LU factorization of double precision matrix
    subroutine zgetrf(mm, nn, aa, lda, ipiv, info)
      import rdp

      !> Number of rows of the matrix
      integer, intent(in) :: mm

      !> Matrix dimension
      integer, intent(in) :: nn

      !> Leading dimension of A
      integer, intent(in) :: lda

      !> Matrix A
      complex(rdp), intent(inout) :: aa(lda, *)

      !> Pivot array
      integer, intent(out) :: ipiv(*)

      !> State of routine on return
      integer, intent(out) :: info
    end subroutine zgetrf
  end interface zgetrf


  !> Computes inverse of a real matrix using LU factorisation
  interface sgetri

    !> Computes inverse of a real matrix using LU factorisation
    subroutine sgetri(nn, aa, lda, ipiv, work, lwork, info)
      import rsp

      !> Matrix dimension
      integer, intent(in) :: nn

      !> Leading dimension of A
      integer, intent(in) :: lda

      !> Matrix A
      real(rsp), intent(inout) :: aa(lda, *)

      !> Pivot array
      integer, intent(in) :: ipiv(*)

      !> Workspace
      real(rsp), intent(inout) :: work(*)

      !> Workspace sizing
      integer, intent(in) :: lwork

      !> State of routine on return
      integer, intent(out) :: info
    end subroutine sgetri
  end interface sgetri


  !> Computes inverse of a double precision matrix using LU factorisation
  interface dgetri

    !> Computes inverse of a double precision matrix using LU factorisation
    subroutine dgetri(nn, aa, lda, ipiv, work, lwork, info)
      import rdp

      !> Matrix dimension
      integer, intent(in) :: nn

      !> Leading dimension of A
      integer, intent(in) :: lda

      !> Matrix A
      real(rdp), intent(inout) :: aa(lda, *)

      !> Pivot array
      integer, intent(in) :: ipiv(*)

      !> Workspace
      real(rdp), intent(inout) :: work(*)

      !> Workspace sizing
      integer, intent(in) :: lwork

      !> State of routine on return
      integer, intent(out) :: info
    end subroutine dgetri
  end interface dgetri


  !> Factorise a real symmetric matrix as A = U*D*U**T or A = L*D*L**T
  interface ssytrf

    !> Factorise a real symmetric matrix as A = U*D*U**T or A = L*D*L**T
    subroutine ssytrf(uplo, nn, aa, lda, ipiv, work, lwork, info)
      import rsp

      !> Upper 'U' or lower 'L' triangle
      character, intent(in) :: uplo

      !> Matrix dimension
      integer, intent(in) :: nn

      !> Leading dimension of A
      integer, intent(in) :: lda

      !> Matrix A
      real(rsp), intent(inout) :: aa(lda, *)

      !> Pivot array
      integer, intent(out) :: ipiv(*)

      !> Workspace
      real(rsp), intent(inout) :: work(*)

      !> Workspace sizing
      integer, intent(in) :: lwork

      !> State of routine on return
      integer, intent(out) :: info
    end subroutine ssytrf
  end interface ssytrf


  !> Factorise a double precision symmetric matrix as A = U*D*U**T or A = L*D*L**T
  interface dsytrf

    !> Factorise a double precision symmetric matrix as A = U*D*U**T or A = L*D*L**T
    subroutine dsytrf(uplo, nn, aa, lda, ipiv, work, lwork, info)
      import rdp

      !> Upper 'U' or lower 'L' triangle
      character, intent(in) :: uplo

      !> Matrix dimension
      integer, intent(in) :: nn

      !> Leading dimension of A
      integer, intent(in) :: lda

      !> Matrix A
      real(rdp), intent(inout) :: aa(lda, *)

      !> Pivot array
      integer, intent(out) :: ipiv(*)

      !> Workspace
      real(rdp), intent(inout) :: work(*)

      !> Workspace sizing
      integer, intent(in) :: lwork

      !> State of routine on return
      integer, intent(out) :: info
    end subroutine dsytrf
  end interface dsytrf


  !> Returns a vector of real random numbers from a uniform or normal distribution
  interface slarnv

    !> Returns a vector of real random numbers from a uniform or normal distribution
    subroutine slarnv(idist, iseed, nn, xx)
      import rsp

      !> Distribution choice
      integer, intent(in) :: idist

      !> Generator seed
      integer, intent(inout) :: iseed(4)

      !> Vector dimension
      integer, intent(in) :: nn

      !> Random values on exit
      real(rsp), intent(out) :: xx(*)
    end subroutine slarnv
  end interface slarnv


  !> Returns a vector of double precision random numbers from a uniform or normal distribution
  interface dlarnv

    !> Returns a vector of double precision random numbers from a uniform or normal distribution
    subroutine dlarnv(idist, iseed, nn, xx)
      import rdp

      !> Distribution choice
      integer, intent(in) :: idist

      !> Generator seed
      integer, intent(inout) :: iseed(4)

      !> Vector dimension
      integer, intent(in) :: nn

      !> Random values on exit
      real(rdp), intent(out) :: xx(*)
    end subroutine dlarnv
  end interface dlarnv


  !> Returns a vector of complex random numbers from a uniform or normal distribution
  interface clarnv

    !> Returns a vector of complex random numbers from a uniform or normal distribution
    subroutine clarnv(idist, iseed, nn, xx)
      import rsp

      !> Distribution choice
      integer, intent(in) :: idist

      !> Generator seed
      integer, intent(inout) :: iseed(4)

      !> Vector dimension
      integer, intent(in) :: nn

      !> Random values on exit
      complex(rsp), intent(out) :: xx(*)
    end subroutine clarnv
  end interface clarnv


  !> Returns a vector of double complex random numbers from a uniform or normal distribution
  interface zlarnv

    !> Returns a vector of double complex random numbers from a uniform or normal distribution
    subroutine zlarnv(idist, iseed, nn, xx)
      import rdp

      !> Distribution choice
      integer, intent(in) :: idist

      !> Generator seed
      integer, intent(inout) :: iseed(4)

      !> Vector dimension
      integer, intent(in) :: nn

      !> Random values on exit
      complex(rdp), intent(out) :: xx(*)
    end subroutine zlarnv
  end interface zlarnv


  !> Provides problem-dependent LAPACK routine parameters for the local environment
  interface ilaenv

    !> Provides problem-dependent LAPACK routine parameters for the local environment
    function ilaenv(ispec, name, opts, n1, n2, n3, n4)

      !> Specifies the parameter to be returned
      integer, intent(in) :: ispec

      !> Name of alling subroutine
      character, intent(in) :: name

      !> The character options to the subroutine NAME, concatenated together
      character, intent(in) :: opts

      !> Problem dimensions for the subroutine NAME
      integer, intent(in) :: n1

      !> Problem dimensions for the subroutine NAME
      integer, intent(in) :: n2

      !> Problem dimensions for the subroutine NAME
      integer, intent(in) :: n3

      !> Problem dimensions for the subroutine NAME
      integer, intent(in) :: n4

      !> Returned parameter
      integer :: ilaenv
    end function ilaenv
  end interface ilaenv


  !> Single precision machine parameters
  interface slamch

    !> Single precision machine parameters
    function slamch(cmach)
      import rsp

      !> Name of parameter to return
      character, intent(in) :: cmach

      !> Parameter value
      real(rsp) :: slamch
    end function slamch
  end interface slamch


  !> Double precision machine parameters
  interface dlamch

    !> Double precision machine parameters
    function dlamch(cmach)
      import rdp

      !> Name of parameter to return
      character, intent(in) :: cmach

      !> Parameter value
      real(rdp) :: dlamch
    end function dlamch
  end interface dlamch


  !> Error handler for the LAPACK routines
  interface xerbla

    !> Error handler for the LAPACK routines
    subroutine xerbla(srname, info)

      !> Calling subroutine name
      character(6), intent(in) :: srname

      !> Info state of the routine
      integer, intent(in) :: info
    end subroutine xerbla
  end interface xerbla

  !> Real singular value decomposition
  interface rgesvd

    !> Real singular value decomposition
    subroutine rgesvd(jobu, jobvt, m, n, a, lda, s, u, ldu, vt, ldvt, work, lwork, info)
      import rsp

      !> Job type for vt
      character, intent(in) :: jobvt

      !> Job type for u
      character, intent(in) :: jobu

      !> First matrix dimension for A
      integer, intent(in) :: m

      !> Second matrix dimension for A
      integer, intent(in) :: n

      !> Leading dimension of A
      integer, intent(in) :: lda

      !> Leading dimension of U
      integer, intent(in) :: ldu

      !> Leading dimension of Vt
      integer, intent(in) :: ldvt

      !> Matrix to decompose
      real(rsp), intent(inout) :: a(lda,*)

      !> Singular values on return min(m,n)
      real(rsp), intent(out) :: s(*)

      !> Left singular vectors
      real(rsp), intent(out) :: u(ldu,*)

      !> Right singular vectors
      real(rsp), intent(out) :: vt(ldvt,*)

      !> Work space
      real(rsp), intent(out) :: work(*)

      !> Size of real work space
      integer, intent(in) :: lwork

      !> State of routine on return
      integer, intent(in) :: info

    end subroutine rgesvd

  end interface rgesvd

  !> Double real singular value decomposition
  interface dgesvd

    !> Double real singular value decomposition
    subroutine dgesvd(jobu, jobvt, m, n, a, lda, s, u, ldu, vt, ldvt, work, lwork, info)
      import rdp

      !> Job type for vt
      character, intent(in) :: jobvt

      !> Job type for u
      character, intent(in) :: jobu

      !> First matrix dimension for A
      integer, intent(in) :: m

      !> Second matrix dimension for A
      integer, intent(in) :: n

      !> Leading dimension of A
      integer, intent(in) :: lda

      !> Leading dimension of U
      integer, intent(in) :: ldu

      !> Leading dimension of Vt
      integer, intent(in) :: ldvt

      !> Matrix to decompose
      real(rdp), intent(inout) :: a(lda,*)

      !> Singular values on return min(m,n)
      real(rdp), intent(out) :: s(*)

      !> Left singular vectors
      real(rdp), intent(out) :: u(ldu,*)

      !> Right singular vectors
      real(rdp), intent(out) :: vt(ldvt,*)

      !> Work space
      real(rdp), intent(out) :: work(*)

      !> Size of work space
      integer, intent(in) :: lwork

      !> State of routine on return
      integer, intent(in) :: info

    end subroutine dgesvd

  end interface dgesvd

  !> Complex singular value decomposition
  interface cgesvd

    !> Complex singular value decomposition
    subroutine cgesvd(jobu, jobvt, m, n, a, lda, s, u, ldu, vt, ldvt, work, lwork, rwork, info)
      import rsp

      !> Job type for vt
      character, intent(in) :: jobvt

      !> Job type for u
      character, intent(in) :: jobu

      !> First matrix dimension for A
      integer, intent(in) :: m

      !> Second matrix dimension for A
      integer, intent(in) :: n

      !> Leading dimension of A
      integer, intent(in) :: lda

      !> Leading dimension of U
      integer, intent(in) :: ldu

      !> Leading dimension of Vt
      integer, intent(in) :: ldvt

      !> Matrix to decompose
      complex(rsp), intent(inout) :: a(lda,*)

      !> Singular values on return min(m,n)
      real(rsp), intent(out) :: s(*)

      !> Real workspace
      real(rsp), intent(out) :: rwork(*)

      !> Left singular vectors
      complex(rsp), intent(out) :: u(ldu,*)

      !> Right singular vectors
      complex(rsp), intent(out) :: vt(ldvt,*)

      !> Complex work space
      complex(rsp), intent(out) :: work(*)

      !> Size of complex work space
      integer, intent(in) :: lwork

      !> State of routine on return
      integer, intent(in) :: info

    end subroutine cgesvd

  end interface cgesvd

  !> Double complex singular value decomposition
  interface zgesvd

    !> Double complex singular value decomposition
    subroutine zgesvd(jobu, jobvt, m, n, a, lda, s, u, ldu, vt, ldvt, work, lwork, rwork, info)
      import rdp

      !> Job type for vt
      character, intent(in) :: jobvt

      !> Job type for u
      character, intent(in) :: jobu

      !> First matrix dimension for A
      integer, intent(in) :: m

      !> Second matrix dimension for A
      integer, intent(in) :: n

      !> Leading dimension of A
      integer, intent(in) :: lda

      !> Leading dimension of U
      integer, intent(in) :: ldu

      !> Leading dimension of Vt
      integer, intent(in) :: ldvt

      !> Matrix to decompose
      complex(rdp), intent(inout) :: a(lda,*)

      !> Singular values on return min(m,n)
      real(rdp), intent(out) :: s(*)

      !> Real workspace
      real(rdp), intent(out) :: rwork(*)

      !> Left singular vectors
      complex(rdp), intent(out) :: u(ldu,*)

      !> Right singular vectors
      complex(rdp), intent(out) :: vt(ldvt,*)

      !> Complex work space
      complex(rdp), intent(out) :: work(*)

      !> Size of complex work space
      integer, intent(in) :: lwork

      !> State of routine on return
      integer, intent(in) :: info

    end subroutine zgesvd

  end interface zgesvd


#:for VPREC in [('s'), ('d')]

  interface ${VPREC}$geev

    !> General matrix eigensolver
    subroutine ${VPREC}$geev(jobvl, jobvr, n, a, lda, wr, wi, vl, ldvl, vr, ldvr, work, lwork, info)

      import r${VPREC}$p

      !> Job type for left vector
      character, intent(in) :: jobvl

      !> Job type for right vector
      character, intent(in) :: jobvr

      !> Second matrix dimension for A
      integer, intent(in) :: n

      !> Matrix to decompose
      real(r${VPREC}$p), intent(inout) :: a(lda,*)

      !> Leading dimension of A
      integer, intent(in) :: lda

      !> Real part of eigenvalues
      real(r${VPREC}$p), intent(out) :: wr(*)

      !> Imaginary part of eigenvalues
      real(r${VPREC}$p), intent(out) :: wi(*)

      !> Left eigenvectors
      real(r${VPREC}$p), intent(out) :: vl(ldvl, *)

      !> Leading dimension of vl
      integer, intent(in) :: ldvl

      !> Right eigenvectors
      real(r${VPREC}$p), intent(out) :: vr(ldvr, *)

      !> Leading dimension of vr
      integer, intent(in) :: ldvr

      !> Work array
      real(r${VPREC}$p), intent(out) :: work(*)

      !> Work array size
      integer, intent(in) :: lwork

      !> State of routine on return
      integer, intent(in) :: info

    end subroutine ${VPREC}$geev

  end interface ${VPREC}$geev

#:endfor

end module dftbp_extlibs_lapack
