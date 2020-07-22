!--------------------------------------------------------------------------------------------------!
!  DFTB+: general package for performing fast atomistic simulations                                !
!  Copyright (C) 2006 - 2020  DFTB+ developers group                                               !
!                                                                                                  !
!  See the LICENSE file for terms of usage and distribution.                                       !
!--------------------------------------------------------------------------------------------------!

#! Note: This module contains preprocessor variable substitutions in subroutine names (${NAME}$)
#! which may break the documentation system. Make sure you preprocess this file before passing it
#! to a source code documentation tool.

#:include 'common.fypp'

!> Interface wrapper for the blas routines.
!>
!> ALL BLAS routines which are called from the main code must be included here.
module dftbp_blas
  use dftbp_accuracy, only : rsp, rdp
  public

  interface


    !> performs the symmetric rank 1 operation
    !> A := alpha*x*x**T + A
    subroutine ssyr(uplo, nn, alpha, xx, incx, aa, lda)
      import rsp

      !> Upper 'U' or lower 'L' triangle
      character, intent(in) :: uplo

      !> matrix  size
      integer, intent(in) :: nn

      !> scaling factor
      real(rsp), intent(in) :: alpha

      !> vector
      real(rsp), intent(in) :: xx(*)

      !> stride
      integer, intent(in) :: incx

      !> leading matrix dimension
      integer, intent(in) :: lda

      !> matrix A
      real(rsp), intent(inout) :: aa(lda, *)
    end subroutine ssyr


    !> performs the symmetric rank 1 operation
    !> A := alpha*x*x**T + A
    subroutine dsyr(uplo, nn, alpha, xx, incx, aa, lda)
      import rdp

      !> Upper 'U' or lower 'L' triangle
      character, intent(in) :: uplo

      !> matrix  size
      integer, intent(in) :: nn

      !> scale factor
      real(rdp), intent(in) :: alpha

      !> vector
      real(rdp), intent(in) :: xx(*)

      !> stride
      integer, intent(in) :: incx

      !> leading matrix dimension
      integer, intent(in) :: lda

      !> matrix A
      real(rdp), intent(inout) :: aa(lda, *)
    end subroutine dsyr


    !> performs the hermitian rank 1 operation
    !> A := alpha*x*x**H + A,
    subroutine cher(uplo, nn, alpha, xx, incx, aa, lda)
      import rsp

      !> Upper 'U' or lower 'L' triangle
      character, intent(in) :: uplo

      !> matrix  size
      integer, intent(in) :: nn

      !> scaling factor
      real(rsp), intent(in) :: alpha

      !> vector
      complex(rsp), intent(in) :: xx(*)

      !> stride
      integer, intent(in) :: incx

      !> leading matrix dimension
      integer, intent(in) :: lda

      !> matrix A
      complex(rsp), intent(inout) :: aa(lda, *)
    end subroutine cher


    !> performs the hermitian rank 1 operation
    !> A := alpha*x*x**H + A,
    subroutine zher(uplo, nn, alpha, xx, incx, aa, lda)
      import rdp

      !> Upper 'U' or lower 'L' triangle
      character, intent(in) :: uplo

      !> matrix  size
      integer, intent(in) :: nn

      !> scale factor
      real(rdp), intent(in) :: alpha

      !> vector
      complex(rdp), intent(in) :: xx(*)

      !> stride
      integer, intent(in) :: incx

      !> leading matrix dimension
      integer, intent(in) :: lda

      !> matrix A
      complex(rdp), intent(inout) :: aa(lda, *)
    end subroutine zher


    !> performs the rank 1 operation
    !> A := alpha*x*y**T + A,
    subroutine sger(mm, nn, alpha, xx, incx, yy, incy, aa, lda)
      import rsp

      !> matrix sizing
      integer, intent(in) :: mm

      !> matrix  size
      integer, intent(in) :: nn

      !> scaling factor
      real(rsp), intent(in) :: alpha

      !> vector
      real(rsp), intent(in) :: xx(*)

      !> stride
      integer, intent(in) :: incx

      !> vector
      real(rsp), intent(in) :: yy(*)

      !> stride
      integer, intent(in) :: incy

      !> leading matrix dimension
      integer, intent(in) :: lda

      !> matrix A
      real(rsp), intent(inout) :: aa(lda, *)
    end subroutine sger


    !> performs the rank 1 operation
    !> A := alpha*x*y**T + A,
    subroutine dger(mm, nn, alpha, xx, incx, yy, incy, aa, lda)
      import rdp

      !> matrix sizing
      integer, intent(in) :: mm

      !> matrix  size
      integer, intent(in) :: nn

      !> scale factor
      real(rdp), intent(in) :: alpha

      !> vector
      real(rdp), intent(in) :: xx(*)

      !> stride
      integer, intent(in) :: incx

      !> vector
      real(rdp), intent(in) :: yy(*)

      !> stride
      integer, intent(in) :: incy

      !> leading matrix dimension
      integer, intent(in) :: lda

      !> matrix A
      real(rdp), intent(inout) :: aa(lda, *)
    end subroutine dger


    !> performs the rank 1 operation
    !> A := alpha*x*y**H + A,
    subroutine cgerc(mm, nn, alpha, xx, incx, yy, incy, aa, lda)
      import rsp

      !> matrix sizing
      integer, intent(in) :: mm

      !> matrix  size
      integer, intent(in) :: nn

      !> scale factor
      complex(rsp), intent(in) :: alpha

      !> vector
      complex(rsp), intent(in) :: xx(*)

      !> stride
      integer, intent(in) :: incx

      !> vector
      complex(rsp), intent(in) :: yy(*)

      !> stride
      integer, intent(in) :: incy

      !> leading matrix dimension
      integer, intent(in) :: lda

      !> matrix A
      complex(rsp), intent(inout) :: aa(lda, *)
    end subroutine cgerc


    !> performs the rank 1 operation
    !> A := alpha*x*y**H + A,
    subroutine zgerc(mm, nn, alpha, xx, incx, yy, incy, aa, lda)
      import rdp

      !> matrix sizing
      integer, intent(in) :: mm

      !> matrix  size
      integer, intent(in) :: nn

      !> scale factor
      complex(rdp), intent(in) :: alpha

      !> vector
      complex(rdp), intent(in) :: xx(*)

      !> stride
      integer, intent(in) :: incx

      !> vector
      complex(rdp), intent(in) :: yy(*)

      !> stride
      integer, intent(in) :: incy

      !> leading matrix dimension
      integer, intent(in) :: lda

      !> matrix A
      complex(rdp), intent(inout) :: aa(lda, *)
    end subroutine zgerc


    !> performs the matrix-vector operation
    !> y := alpha*A*x + beta*y
    subroutine ssymv(uplo, nn, alpha, aa, lda, xx, incx, beta, yy, incy)
      import rsp

      !> Upper 'U' or lower 'L' triangle
      character, intent(in) :: uplo

      !> matrix  size
      integer, intent(in) :: nn

      !> scaling factor
      real(rsp), intent(in) :: alpha

      !> leading matrix dimension
      integer, intent(in) :: lda

      !> matrix A
      real(rsp), intent(in) :: aa(lda, *)

      !> vector
      real(rsp), intent(in) :: xx(*)

      !> stride
      integer, intent(in) :: incx

      !> scale factor
      real(rsp), intent(in) :: beta

      !> vector
      real(rsp), intent(inout) :: yy(*)

      !> stride
      integer, intent(in) :: incy
    end subroutine ssymv


    !> performs the matrix-vector operation
    !> y := alpha*A*x + beta*y
    subroutine dsymv(uplo, nn, alpha, aa, lda, xx, incx, beta, yy, incy)
      import rdp

      !> Upper 'U' or lower 'L' triangle
      character, intent(in) :: uplo

      !> matrix  size
      integer, intent(in) :: nn

      !> scale factor
      real(rdp), intent(in) :: alpha

      !> leading matrix dimension
      integer, intent(in) :: lda

      !> matrix A
      real(rdp), intent(in) :: aa(lda, *)

      !> vector
      real(rdp), intent(in) :: xx(*)

      !> stride
      integer, intent(in) :: incx

      !> scaling factor
      real(rdp), intent(in) :: beta

      !> vector
      real(rdp), intent(inout) :: yy(*)

      !> stride
      integer, intent(in) :: incy
    end subroutine dsymv


    !> performs the matrix-vector operation
    !> y := alpha*A*x + beta*y
    subroutine chemv(uplo, nn, alpha, aa, lda, xx, incx, beta, yy, incy)
      import rsp

      !> Upper 'U' or lower 'L' triangle
      character, intent(in) :: uplo

      !> matrix  size
      integer, intent(in) :: nn

      !> scale factor
      complex(rsp), intent(in) :: alpha

      !> leading matrix dimension
      integer, intent(in) :: lda

      !> matrix A
      complex(rsp), intent(in) :: aa(lda, *)

      !> vector
      complex(rsp), intent(in) :: xx(*)

      !> stride
      integer, intent(in) :: incx

      !> scaling factor
      complex(rsp), intent(in) :: beta

      !> vector
      complex(rsp), intent(inout) :: yy(*)

      !> stride
      integer, intent(in) :: incy
    end subroutine chemv


    !> performs the matrix-vector operation
    !> y := alpha*A*x + beta*y
    subroutine zhemv(uplo, nn, alpha, aa, lda, xx, incx, beta, yy, incy)
      import rdp

      !> Upper 'U' or lower 'L' triangle
      character, intent(in) :: uplo

      !> matrix  size
      integer, intent(in) :: nn

      !> scale factor
      complex(rdp), intent(in) :: alpha

      !> leading matrix dimension
      integer, intent(in) :: lda

      !> matrix A
      complex(rdp), intent(in) :: aa(lda, *)

      !> vector
      complex(rdp), intent(in) :: xx(*)

      !> stride
      integer, intent(in) :: incx

      !> scale factor
      complex(rdp), intent(in) :: beta

      !> vector
      complex(rdp), intent(inout) :: yy(*)

      !> stride
      integer, intent(in) :: incy
    end subroutine zhemv


    !> performs one of the matrix-vector operations
    !> y := alpha*A*x + beta*y, or y := alpha*A**T*x + beta*y,
    subroutine sgemv(trans, mm, nn, alpha, aa, lda, xx, incx, beta, yy, incy)
      import rsp

      !> should transposition be used
      character, intent(in) :: trans

      !> matrix sizing
      integer, intent(in) :: mm

      !> matrix  size
      integer, intent(in) :: nn

      !> scaling factor
      real(rsp), intent(in) :: alpha

      !> leading matrix dimension
      integer, intent(in) :: lda

      !> matrix A
      real(rsp), intent(in) :: aa(lda, *)

      !> vector
      real(rsp), intent(in) :: xx(*)

      !> stride
      integer, intent(in) :: incx

      !> scale factor
      real(rsp), intent(in) :: beta

      !> vector
      real(rsp), intent(inout) :: yy(*)

      !> stride
      integer, intent(in) :: incy
    end subroutine sgemv


    !> performs one of the matrix-vector operations
    !> y := alpha*A*x + beta*y, or y := alpha*A**T*x + beta*y,
    subroutine dgemv(trans, mm, nn, alpha, aa, lda, xx, incx, beta, yy, incy)
      import rdp

      !> should transposition be used
      character, intent(in) :: trans

      !> matrix sizing
      integer, intent(in) :: mm

      !> matrix  size
      integer, intent(in) :: nn

      !> scale factor
      real(rdp), intent(in) :: alpha

      !> leading matrix dimension
      integer, intent(in) :: lda

      !> matrix A
      real(rdp), intent(in) :: aa(lda, *)

      !> vector
      real(rdp), intent(in) :: xx(*)

      !> stride
      integer, intent(in) :: incx

      !> scaling factor
      real(rdp), intent(in) :: beta

      !> vector
      real(rdp), intent(inout) :: yy(*)

      !> stride
      integer, intent(in) :: incy
    end subroutine dgemv


    !> performs the matrix-vector  operation
    !> y := alpha*A*x + beta*y,
    subroutine ssbmv(uplo, nn, kk, alpha, aa, lda, xx, incx, beta, yy, incy)
      import rsp

      !> Upper 'U' or lower 'L' triangle
      character, intent(in) :: uplo

      !> matrix  size
      integer, intent(in) :: nn

      !> number of superdiagonals
      integer, intent(in) :: kk

      !> scaling factor
      real(rsp), intent(in) :: alpha

      !> leading matrix dimension
      integer, intent(in) :: lda

      !> matrix A
      real(rsp), intent(in) :: aa(lda, *)

      !> vector
      real(rsp), intent(in) :: xx(*)

      !> stride
      integer, intent(in) :: incx

      !> scale factor
      real(rsp), intent(in) :: beta

      !> vector
      real(rsp), intent(inout) :: yy(*)

      !> stride
      integer, intent(in) :: incy
    end subroutine ssbmv


    !> performs the matrix-vector  operation
    !> y := alpha*A*x + beta*y,
    subroutine dsbmv(uplo, nn, kk, alpha, aa, lda, xx, incx, beta, yy, incy)
      import rdp

      !> Upper 'U' or lower 'L' triangle
      character, intent(in) :: uplo

      !> matrix  size
      integer, intent(in) :: nn

      !> number of superdiagonals
      integer, intent(in) :: kk

      !> scale factor
      real(rdp), intent(in) :: alpha

      !> leading matrix dimension
      integer, intent(in) :: lda

      !> matrix A
      real(rdp), intent(in) :: aa(lda, *)

      !> vector
      real(rdp), intent(in) :: xx(*)

      !> stride
      integer, intent(in) :: incx

      !> scaling factor
      real(rdp), intent(in) :: beta

      !> vector
      real(rdp), intent(inout) :: yy(*)

      !> stride
      integer, intent(in) :: incy
    end subroutine dsbmv


    !> performs one of the matrix-matrix operations
    !> C := alpha*A*B + beta*C, or C := alpha*B*A + beta*C,
    subroutine ssymm(side, uplo, mm, nn, alpha, aa, lda, bb, ldb, beta, cc, ldc)
      import rsp

      !> side of the product
      character, intent(in) :: side

      !> upper or lower matrix
      character, intent(in) :: uplo

      !> matrix sizing
      integer, intent(in) :: mm

      !> matrix  size
      integer, intent(in) :: nn

      !> scaling factor
      real(rsp), intent(in) :: alpha

      !> leading matrix dimension
      integer, intent(in) :: lda

      !> matrix A
      real(rsp), intent(in) :: aa(lda, *)

      !> leading matrix dimension
      integer, intent(in) :: ldb

      !> matrix B
      real(rsp), intent(in) :: bb(ldb, *)

      !> scale factor
      real(rsp), intent(in) :: beta

      !> leading matrix dimension
      integer, intent(in) :: ldc

      !> matrix C
      real(rsp), intent(inout) :: cc(ldc, *)
    end subroutine ssymm


    !> performs one of the matrix-matrix operations
    !> C := alpha*A*B + beta*C, or C := alpha*B*A + beta*C,
    subroutine dsymm(side, uplo, mm, nn, alpha, aa, lda, bb, ldb, beta, cc, ldc)
      import rdp

      !> side of the product
      character, intent(in) :: side

      !> upper or lower matrix
      character, intent(in) :: uplo

      !> matrix sizing
      integer, intent(in) :: mm

      !> matrix  size
      integer, intent(in) :: nn

      !> scale factor
      real(rdp), intent(in) :: alpha

      !> leading matrix dimension
      integer, intent(in) :: lda

      !> matrix A
      real(rdp), intent(in) :: aa(lda, *)

      !> leading matrix dimension
      integer, intent(in) :: ldb

      !> matrix B
      real(rdp), intent(in) :: bb(ldb, *)

      !> scaling factor
      real(rdp), intent(in) :: beta

      !> leading matrix dimension
      integer, intent(in) :: ldc

      !> matrix C
      real(rdp), intent(inout) :: cc(ldc, *)
    end subroutine dsymm


    !> performs one of the matrix-matrix operations
    !> C := alpha*A*B + beta*C, or C := alpha*B*A + beta*C,
    subroutine chemm(side, uplo, mm, nn, alpha, aa, lda, bb, ldb, beta, cc, ldc)
      import rsp

      !> side of the product
      character, intent(in) :: side

      !> upper or lower matrix
      character, intent(in) :: uplo

      !> matrix sizing
      integer, intent(in) :: mm

      !> matrix  size
      integer, intent(in) :: nn

      !> scale factor
      complex(rsp), intent(in) :: alpha

      !> leading matrix dimension
      integer, intent(in) :: lda

      !> matrix A
      complex(rsp), intent(in) :: aa(lda, *)

      !> leading matrix dimension
      integer, intent(in) :: ldb

      !> matrix B
      complex(rsp), intent(in) :: bb(ldb, *)

      !> scaling factor
      complex(rsp), intent(in) :: beta

      !> leading matrix dimension
      integer, intent(in) :: ldc

      !> matrix C
      complex(rsp), intent(inout) :: cc(ldc, *)
    end subroutine chemm


    !> performs one of the matrix-matrix operations
    !> C := alpha*A*B + beta*C, or C := alpha*B*A + beta*C,
    subroutine zhemm(side, uplo, mm, nn, alpha, aa, lda, bb, ldb, beta, cc, ldc)
      import rdp

      !> side of the product
      character, intent(in) :: side

      !> upper or lower matrix
      character, intent(in) :: uplo

      !> matrix sizing
      integer, intent(in) :: mm

      !> matrix  size
      integer, intent(in) :: nn

      !> scale factor
      complex(rdp), intent(in) :: alpha

      !> leading matrix dimension
      integer, intent(in) :: lda

      !> matrix A
      complex(rdp), intent(in) :: aa(lda, *)

      !> leading matrix dimension
      integer, intent(in) :: ldb

      !> matrix B
      complex(rdp), intent(in) :: bb(ldb, *)

      !> scale factor
      complex(rdp), intent(in) :: beta

      !> leading matrix dimension
      integer, intent(in) :: ldc

      !> matrix C
      complex(rdp), intent(inout) :: cc(ldc, *)
    end subroutine zhemm


    !> performs one of the matrix-matrix operations
    !> C := alpha*op( A )*op( B ) + beta*C, where  op( X ) is one of
    !> op( X ) = X, or op( X ) = X**T
    subroutine sgemm(transa, transb, mm, nn, kk, alpha, aa, lda, bb, ldb, beta,&
        & cc, ldc)
      import rsp

      !> On entry, TRANSA specifies the form of op( A ) to be used
      character, intent(in) :: transa

      !> on entry specifies op(B)
      !> should transposition be used
      character, intent(in) :: transb

      !> matrix sizing
      integer, intent(in) :: mm

      !> matrix  size
      integer, intent(in) :: nn

      !> shared index size
      integer, intent(in) :: kk

      !> scaling factor
      real(rsp), intent(in) :: alpha

      !> leading matrix dimension
      integer, intent(in) :: lda

      !> matrix A
      real(rsp), intent(in) :: aa(lda, *)

      !> leading matrix dimension
      integer, intent(in) :: ldb

      !> matrix B
      real(rsp), intent(in) :: bb(ldb, *)

      !> scale factor
      real(rsp), intent(in) :: beta

      !> leading matrix dimension
      integer, intent(in) :: ldc

      !> matrix C
      real(rsp), intent(inout) :: cc(ldc, *)
    end subroutine sgemm


    !> performs one of the matrix-matrix operations
    !> C := alpha*op( A )*op( B ) + beta*C, where  op( X ) is one of
    !> op( X ) = X, or op( X ) = X**T
    subroutine dgemm(transa, transb, mm, nn, kk, alpha, aa, lda, bb, ldb, beta,&
        & cc, ldc)
      import rdp

      !> On entry, TRANSA specifies the form of op( A ) to be used
      character, intent(in) :: transa

      !> on entry specifies op(B)
      !> should transposition be used
      character, intent(in) :: transb

      !> matrix sizing
      integer, intent(in) :: mm

      !> matrix  size
      integer, intent(in) :: nn

      !> shared index size
      integer, intent(in) :: kk

      !> scale factor
      real(rdp), intent(in) :: alpha

      !> leading matrix dimension
      integer, intent(in) :: lda

      !> matrix A
      real(rdp), intent(in) :: aa(lda, *)

      !> leading matrix dimension
      integer, intent(in) :: ldb

      !> matrix B
      real(rdp), intent(in) :: bb(ldb, *)

      !> scaling factor
      real(rdp), intent(in) :: beta

      !> leading matrix dimension
      integer, intent(in) :: ldc

      !> matrix C
      real(rdp), intent(inout) :: cc(ldc, *)
    end subroutine dgemm


    !> performs one of the matrix-matrix operations
    !> C := alpha*op( A )*op( B ) + beta*C, where  op( X ) is one of
    !> op( X ) = X, op( X ) = X**T, or op( X ) = X**H
    subroutine cgemm(transa, transb, mm, nn, kk, alpha, aa, lda, bb, ldb, beta,&
        & cc, ldc)
      import rsp

      !> On entry, TRANSA specifies the form of op( A ) to be used
      character, intent(in) :: transa

      !> on entry specifies op(B)
      !> should transposition be used
      character, intent(in) :: transb

      !> matrix sizing
      integer, intent(in) :: mm

      !> matrix  size
      integer, intent(in) :: nn

      !> shared index size
      integer, intent(in) :: kk

      !> scale factor
      complex(rsp), intent(in) :: alpha

      !> leading matrix dimension
      integer, intent(in) :: lda

      !> matrix A
      complex(rsp), intent(in) :: aa(lda, *)

      !> leading matrix dimension
      integer, intent(in) :: ldb

      !> matrix B
      complex(rsp), intent(in) :: bb(ldb, *)

      !> scaling factor
      complex(rsp), intent(in) :: beta

      !> leading matrix dimension
      integer, intent(in) :: ldc

      !> matrix C
      complex(rsp), intent(inout) :: cc(ldc, *)
    end subroutine cgemm


    !> performs one of the matrix-matrix operations
    !> C := alpha*op( A )*op( B ) + beta*C, where  op( X ) is one of
    !> op( X ) = X, op( X ) = X**T, or op( X ) = X**H
    subroutine zgemm(transa, transb, mm, nn, kk, alpha, aa, lda, bb, ldb, beta,&
        & cc, ldc)
      import rdp

      !> On entry, TRANSA specifies the form of op( A ) to be used
      character, intent(in) :: transa

      !> on entry specifies op(B)
      !> should transposition be used
      character, intent(in) :: transb

      !> matrix sizing
      integer, intent(in) :: mm

      !> matrix  size
      integer, intent(in) :: nn

      !> shared index size
      integer, intent(in) :: kk

      !> scale factor
      complex(rdp), intent(in) :: alpha

      !> leading matrix dimension
      integer, intent(in) :: lda

      !> matrix A
      complex(rdp), intent(in) :: aa(lda, *)

      !> leading matrix dimension
      integer, intent(in) :: ldb

      !> matrix B
      complex(rdp), intent(in) :: bb(ldb, *)

      !> scale factor
      complex(rdp), intent(in) :: beta

      !> leading matrix dimension
      integer, intent(in) :: ldc

      !> matrix C
      complex(rdp), intent(inout) :: cc(ldc, *)
    end subroutine zgemm


  #:for VPREC, VTYPE, NAME in [('rsp', 'real', 'ssyrk'), ('rdp', 'real', 'dsyrk'),&
    & ('rsp', 'complex', 'cherk'), ('rdp', 'complex', 'zherk')]
    !> performs one of the symmetric rank k operations
    !> C := alpha*A*A**T + beta*C, or C := alpha*A**T*A + beta*C
    subroutine ${NAME}$(uplo, trans, nn, kk, alpha, aa, lda, beta, cc, ldc)
      import ${VPREC}$

      !> Upper 'U' or lower 'L' triangle
      character, intent(in) :: uplo

      !> should transposition be used
      character, intent(in) :: trans

      !> matrix  size
      integer, intent(in) :: nn

      !> rank
      integer, intent(in) :: kk

      !> scaling factor
      real(${VPREC}$), intent(in) :: alpha

      !> leading matrix dimension
      integer, intent(in) :: lda

      !> matrix A
      ${VTYPE}$(${VPREC}$), intent(in) :: aa(lda, *)

      !> scale factor
      real(${VPREC}$), intent(in) :: beta

      !> leading matrix dimension
      integer, intent(in) :: ldc

      !> matrix C
      ${VTYPE}$(${VPREC}$), intent(inout) :: cc(ldc, *)

    end subroutine ${NAME}$
  #:endfor


  #:for VPREC, VTYPE, NAME in [('rsp', 'real', 'ssyr2k'), ('rdp', 'real', 'dsyr2k'),&
    & ('rsp', 'complex', 'cher2k'), ('rdp', 'complex', 'zher2k')]
    !> performs one of the symmetric rank 2k operations
    !> C := 0.5*alpha*(A*B**H + B**H*A)+ beta*C, or C := 0.5*alpha*(B*A**H + A**H*B)+ beta*C
    subroutine ${NAME}$(uplo, trans, nn, kk, alpha, aa, lda, bb, ldb, beta, cc, ldc)
      import ${VPREC}$

      !> Upper 'U' or lower 'L' triangle
      character, intent(in) :: uplo

      !> should transposition be used
      character, intent(in) :: trans

      !> matrix  size
      integer, intent(in) :: nn

      !> rank
      integer, intent(in) :: kk

      !> scaling factor
      ${VTYPE}$(${VPREC}$), intent(in) :: alpha

      !> leading matrix dimension for a
      integer, intent(in) :: lda

      !> matrix A
      ${VTYPE}$(${VPREC}$), intent(in) :: aa(lda, *)

      !> leading matrix dimension for b
      integer, intent(in) :: ldb

      !> matrix B
      ${VTYPE}$(${VPREC}$), intent(in) :: bb(lda, *)

      !> scale factor
      ${VTYPE}$(${VPREC}$), intent(in) :: beta

      !> leading matrix dimension
      integer, intent(in) :: ldc

      !> matrix C
      ${VTYPE}$(${VPREC}$), intent(inout) :: cc(ldc, *)

    end subroutine ${NAME}$
  #:endfor


    !> solves one of the matrix equations
    !> op( A )*X = alpha*B, or X*op( A ) = alpha*B
    subroutine strsm(side, uplo, transa, diag, mm, nn, alpha, aa, lda, bb, ldb)
      import rsp

      !> side of the product
      character, intent(in) :: side

      !> Upper 'U' or lower 'L' triangle
      character, intent(in) :: uplo

      !> On entry, TRANSA specifies the form of op( A ) to be used
      character, intent(in) :: transa

      !> On entry, DIAG specifies whether or not A is unit triangular
      character, intent(in) :: diag

      !> matrix sizing
      integer, intent(in) :: mm

      !> matrix  size
      integer, intent(in) :: nn

      !> scaling factor
      real(rsp), intent(in) :: alpha

      !> leading matrix dimension
      integer, intent(in) :: lda

      !> matrix A
      real(rsp), intent(in) :: aa(lda, *)

      !> leading matrix dimension
      integer, intent(in) :: ldb

      !> matrix B
      real(rsp), intent(inout) :: bb(ldb, *)
    end subroutine strsm


    !> solves one of the matrix equations
    !> op( A )*X = alpha*B, or X*op( A ) = alpha*B
    subroutine dtrsm(side, uplo, transa, diag, mm, nn, alpha, aa, lda, bb, ldb)
      import rdp

      !> side of the product
      character, intent(in) :: side

      !> Upper 'U' or lower 'L' triangle
      character, intent(in) :: uplo

      !> On entry, TRANSA specifies the form of op( A ) to be used
      character, intent(in) :: transa

      !> On entry, DIAG specifies whether or not A is unit triangular
      character, intent(in) :: diag

      !> matrix sizing
      integer, intent(in) :: mm

      !> matrix  size
      integer, intent(in) :: nn

      !> scale factor
      real(rdp), intent(in) :: alpha

      !> leading matrix dimension
      integer, intent(in) :: lda

      !> matrix A
      real(rdp), intent(in) :: aa(lda, *)

      !> leading matrix dimension
      integer, intent(in) :: ldb

      !> matrix B
      real(rdp), intent(inout) :: bb(ldb, *)
    end subroutine dtrsm


    !> solves one of the matrix equations
    !> op( A )*X = alpha*B, or X*op( A ) = alpha*B
    subroutine ctrsm(side, uplo, transa, diag, mm, nn, alpha, aa, lda, bb, ldb)
      import rsp

      !> side of the product
      character, intent(in) :: side

      !> Upper 'U' or lower 'L' triangle
      character, intent(in) :: uplo

      !> On entry, TRANSA specifies the form of op( A ) to be used
      character, intent(in) :: transa

      !> On entry, DIAG specifies whether or not A is unit triangular
      character, intent(in) :: diag

      !> matrix sizing
      integer, intent(in) :: mm

      !> matrix  size
      integer, intent(in) :: nn

      !> scale factor
      complex(rsp), intent(in) :: alpha

      !> leading matrix dimension
      integer, intent(in) :: lda

      !> matrix A
      complex(rsp), intent(in) :: aa(lda, *)

      !> leading matrix dimension
      integer, intent(in) :: ldb

      !> matrix B
      complex(rsp), intent(inout) :: bb(ldb, *)
    end subroutine ctrsm


    !> solves one of the matrix equations
    !> op( A )*X = alpha*B, or X*op( A ) = alpha*B
    subroutine ztrsm(side, uplo, transa, diag, mm, nn, alpha, aa, lda, bb, ldb)
      import rdp

      !> side of the product
      character, intent(in) :: side

      !> Upper 'U' or lower 'L' triangle
      character, intent(in) :: uplo

      !> On entry, TRANSA specifies the form of op( A ) to be used
      character, intent(in) :: transa

      !> On entry, DIAG specifies whether or not A is unit triangular
      character, intent(in) :: diag

      !> matrix sizing
      integer, intent(in) :: mm

      !> matrix  size
      integer, intent(in) :: nn

      !> scale factor
      complex(rdp), intent(in) :: alpha

      !> leading matrix dimension
      integer, intent(in) :: lda

      !> matrix A
      complex(rdp), intent(in) :: aa(lda, *)

      !> leading matrix dimension
      integer, intent(in) :: ldb

      !> matrix B
      complex(rdp), intent(inout) :: bb(ldb, *)
    end subroutine ztrsm


    !> performs one of the matrix-matrix operations
    !> B := alpha*op( A )*B, or B := alpha*B*op( A ), where op( A ) = A or op( A ) = A**T.
    subroutine strmm(side, uplo, transa, diag, mm, nn, alpha, aa, lda, bb, ldb)
      import rsp

      !> side of the product
      character, intent(in) :: side

      !> Upper 'U' or lower 'L' triangle
      character, intent(in) :: uplo

      !> On entry, TRANSA specifies the form of op( A ) to be used
      character, intent(in) :: transa

      !> On entry, DIAG specifies whether or not A is unit triangular
      character, intent(in) :: diag

      !> matrix sizing
      integer, intent(in) :: mm

      !> matrix  size
      integer, intent(in) :: nn

      !> scaling factor
      real(rsp), intent(in) :: alpha

      !> leading matrix dimension
      integer, intent(in) :: lda

      !> matrix A
      real(rsp), intent(in) :: aa(lda, *)

      !> leading matrix dimension
      integer, intent(in) :: ldb

      !> matrix B
      real(rsp), intent(inout) :: bb(ldb, *)
    end subroutine strmm


    !> performs one of the matrix-matrix operations
    !> B := alpha*op( A )*B, or B := alpha*B*op( A ), where op( A ) = A or op( A ) = A**T.
    subroutine dtrmm(side, uplo, transa, diag, mm, nn, alpha, aa, lda, bb, ldb)
      import rdp

      !> side of the product
      character, intent(in) :: side

      !> Upper 'U' or lower 'L' triangle
      character, intent(in) :: uplo

      !> On entry, TRANSA specifies the form of op( A ) to be used
      character, intent(in) :: transa

      !> On entry, DIAG specifies whether or not A is unit triangular
      character, intent(in) :: diag

      !> matrix sizing
      integer, intent(in) :: mm

      !> matrix  size
      integer, intent(in) :: nn

      !> scale factor
      real(rdp), intent(in) :: alpha

      !> leading matrix dimension
      integer, intent(in) :: lda

      !> matrix A
      real(rdp), intent(in) :: aa(lda, *)

      !> leading matrix dimension
      integer, intent(in) :: ldb

      !> matrix B
      real(rdp), intent(inout) :: bb(ldb, *)
    end subroutine dtrmm


    !> performs one of the matrix-matrix operations
    !> B := alpha*op( A )*B, or B := alpha*B*op( A ),
    !> where op( A ) = A, op( A ) = A**T or op( A ) = A**H.
    subroutine ctrmm(side, uplo, transa, diag, mm, nn, alpha, aa, lda, bb, ldb)
      import rsp

      !> side of the product
      character, intent(in) :: side

      !> Upper 'U' or lower 'L' triangle
      character, intent(in) :: uplo

      !> On entry, TRANSA specifies the form of op( A ) to be used
      character, intent(in) :: transa

      !> On entry, DIAG specifies whether or not A is unit triangular
      character, intent(in) :: diag

      !> matrix sizing
      integer, intent(in) :: mm

      !> matrix  size
      integer, intent(in) :: nn

      !> scale factor
      complex(rsp), intent(in) :: alpha

      !> leading matrix dimension
      integer, intent(in) :: lda

      !> matrix A
      complex(rsp), intent(in) :: aa(lda, *)

      !> leading matrix dimension
      integer, intent(in) :: ldb

      !> matrix B
      complex(rsp), intent(inout) :: bb(ldb, *)
    end subroutine ctrmm


    !> performs one of the matrix-matrix operations
    !> B := alpha*op( A )*B, or B := alpha*B*op( A ),
    !> where op( A ) = A, op( A ) = A**T or op( A ) = A**H.
    subroutine ztrmm(side, uplo, transa, diag, mm, nn, alpha, aa, lda, bb, ldb)
      import rdp

      !> side of the product
      character, intent(in) :: side

      !> Upper 'U' or lower 'L' triangle
      character, intent(in) :: uplo

      !> On entry, TRANSA specifies the form of op( A ) to be used
      character, intent(in) :: transa

      !> On entry, DIAG specifies whether or not A is unit triangular
      character, intent(in) :: diag

      !> matrix sizing
      integer, intent(in) :: mm

      !> matrix  size
      integer, intent(in) :: nn

      !> scale factor
      complex(rdp), intent(in) :: alpha

      !> leading matrix dimension
      integer, intent(in) :: lda

      !> matrix A
      complex(rdp), intent(in) :: aa(lda, *)

      !> leading matrix dimension
      integer, intent(in) :: ldb

      !> matrix B
      complex(rdp), intent(inout) :: bb(ldb, *)
    end subroutine ztrmm

  end interface

end module dftbp_blas
