!--------------------------------------------------------------------------------------------------!
!  DFTB+: general package for performing fast atomistic simulations                                !
!  Copyright (C) 2017  DFTB+ developers group                                                      !
!                                                                                                  !
!  See the LICENSE file for terms of usage and distribution.                                       !
!--------------------------------------------------------------------------------------------------!

#:include 'common.fypp'

!!* Contains F90 wrapper functions for some commonly used lapack calls needed
!!* in the code
!!* @caveat contains some fixes for lapack 3.0 bugs, if this gets corrected in
!!* lapack 4.x they should be removed
module eigensolver
  use assert
  use message
  use accuracy, only : rsp, rdp
  use blas
  use lapack
  implicit none
  private

  public :: heev, hegv, hegvd, gvr, bgv

  character(len=100) :: error_string !* Used to return runtime diagnostics

  !!* Simple eigensolver for a symmetric/Hermitian matrix
  !!* @param a contains the matrix for the solver, returns eigenvalues if
  !!* requested
  !!* @param w eigenvalues
  !!* @param uplo upper or lower triangle of the matrix
  !!* @param jobz compute eigenvalues 'N' or eigenvalues and eigenvectors 'V'
  !!* @caveat the matrix a is overwritten
  interface heev
    module procedure real_ssyev
    module procedure dble_dsyev
    module procedure cmplx_cheev
    module procedure dblecmplx_zheev
  end interface heev

!  !!* Simple eigensolver for a general matrix
!  !!* @param a contains the matrix for the solver, returns eigenvalues if
!  !!* requested
!  !!* @param w eigenvalues
!  !!* @param vr
!  !!* @param vl
!  !!* @param jobz compute eigenvalues 'N' or eigenvalues and eigenvectors 'V'
!  !!* @caveat the matrix a is overwritten
!  interface geev
!    module procedure real_sgeev
!    module procedure dble_dgeev
!    module procedure cmplx_cgeev
!    module procedure dblecmplx_zgeev
!  end interface geev

  !!* Simple eigensolver for a symmetric/Hermitian generalized matrix problem
  !!* @param a contains the matrix for the solver, returns eigenvalues if
  !!* requested
  !!* @param b contains the second matrix for the solver
  !!* @param w eigenvalues
  !!* @param uplo upper or lower triangle of the matrix
  !!* @param jobz compute eigenvalues 'N' or eigenvalues and eigenvectors 'V'
  !!* @param itype optional specifies the problem type to be solved
  !!* 1:A*x=(lambda)*B*x, 2:A*B*x=(lambda)*x, 3:B*A*x=(lambda)*x
  !!* default is 1
  !!* @caveat the matrix a is overwritten
  !!* @caveat the matrix b is overwritted with Cholesky factorization if
  !!* eigenvalues are computed
  interface hegv
    module procedure real_ssygv
    module procedure dble_dsygv
    module procedure cmplx_chegv
    module procedure dblecmplx_zhegv
  end interface hegv

  !!* Simple eigensolver for a symmetric/Hermitian generalized matrix problem
  !!* using divide and conquer
  !!* @param a contains the matrix for the solver, returns eigenvalues if
  !!* requested
  !!* @param b contains the second matrix for the solver
  !!* @param w eigenvalues
  !!* @param uplo upper or lower triangle of the matrix
  !!* @param jobz compute eigenvalues 'N' or eigenvalues and eigenvectors 'V'
  !!* @param itype optional specifies the problem type to be solved
  !!* 1:A*x=(lambda)*B*x, 2:A*B*x=(lambda)*x, 3:B*A*x=(lambda)*x
  !!* default is 1
  !!* @caveat the matrix a is overwritten
  !!* @caveat the matrix b is overwritted with Cholesky factorization if
  !!* eigenvalues are computed
  interface hegvd
    module procedure real_ssygvd
    module procedure dble_dsygvd
    module procedure cmplx_chegvd
    module procedure dblecmplx_zhegvd
  end interface hegvd

  !!* Simple eigensolver for a symmetric/Hermitian generalized matrix problem
  !!* using the lapack relatively robust representation solver, based on the
  !!* SYGV source. If the requested number of eigenvalues is lower than
  !!* the size of H/S suspace mode is used (optionally the range can be
  !!* set using il and ul) to return the lowest eigenvalues/vectors of number
  !!* size(w)
  !!* @param a contains the matrix for the solver, returns eigenvalues if
  !!* requested
  !!* @param b contains the second matrix for the solver
  !!* @param w eigenvalues
  !!* @param uplo upper or lower triangle of the matrix
  !!* @param jobz compute eigenvalues 'N' or eigenvalues and eigenvectors 'V'
  !!* @param itype optional specifies the problem type to be solved
  !!* 1:A*x=(lambda)*B*x, 2:A*B*x=(lambda)*x, 3:B*A*x=(lambda)*x
  !!* default is 1
  !!* @param il optional lower range
  !!* @param ul optional upper range
  !!* @caveat the matrix a is overwritten
  !!* @caveat the matrix b is overwritten
  interface gvr
    module procedure real_ssygvr
    module procedure dble_dsygvr
    module procedure cmplx_chegvr
    module procedure dblecmplx_zhegvr
  end interface

  !!* Simple eigensolver for a symmetric/Hermitian banded generalized matrix
  !!* problem of the form A*x=(lambda)*B*x
  !!* @param ab contains the matrix for the solver, returns eigenvalues if
  !!* requested
  !!* @param bb contains the second matrix for the solver
  !!* @param w eigenvalues
  !!* @param uplo upper or lower triangle of the matrix
  !!* @param itype optional specifies the problem type to be solved
  !!* @caveat the matrix ab is overwritten
  !!* @caveat the matrix bb is overwritted with a split Cholesky factorization
  interface bgv
    module procedure real_ssbgv
    module procedure dble_dsbgv
    module procedure cmplx_chbgv
    module procedure dblecmplx_zhbgv
  end interface

contains

  !!* Real eigensolver for a symmetric matrix
  subroutine real_ssyev(a,w,uplo,jobz)
    real(rsp), intent(inout) :: a(:,:)
    real(rsp), intent(out) :: w(:)
    real(rsp), allocatable :: work(:)
    character, intent(in) :: uplo
    character, intent(in) :: jobz
    integer n, info
    integer :: int_idealwork
    real(rsp) :: idealwork(1)
    @:ASSERT(uplo == 'u' .or. uplo == 'U' .or. uplo == 'l' .or. uplo == 'L')
    @:ASSERT(jobz == 'n' .or. jobz == 'N' .or. jobz == 'v' .or. jobz == 'V')
    @:ASSERT(all(shape(a)==size(w,dim=1)))
    n=size(a,dim=1)
    @:ASSERT(n>0)
    call ssyev(jobz, uplo, n, a, n, w, idealwork, -1, info)
    if (info/=0) then
      call error("Failue in SSYEV to determine optimum workspace")
    endif
    int_idealwork=floor(idealwork(1))
    allocate(work(int_idealwork))
    call SSYEV(jobz, uplo, n, a, n, w, work, int_idealwork, info)
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

  End Subroutine real_ssyev

  !!* Double precision eigensolver for a symmetric matrix
  Subroutine dble_dsyev(a,w,uplo,jobz)
    real(rdp), intent(inout) :: a(:,:)
    real(rdp), intent(out) :: w(:)
    real(rdp), allocatable :: work(:)
    character, intent(in) :: uplo
    character, intent(in) :: jobz
    integer n, info
    integer :: int_idealwork
    real(rdp) :: idealwork(1)
    @:ASSERT(uplo == 'u' .or. uplo == 'U' .or. uplo == 'l' .or. uplo == 'L')
    @:ASSERT(jobz == 'n' .or. jobz == 'N' .or. jobz == 'v' .or. jobz == 'V')
    @:ASSERT(all(shape(a)==size(w,dim=1)))
    n=size(a,dim=1)
    @:ASSERT(n>0)
    call dsyev(jobz, uplo, n, a, n, w, idealwork, -1, info)
    if (info/=0) then
      call error("Failue in DSYEV to determine optimum workspace")
    endif
    int_idealwork=floor(idealwork(1))
    allocate(work(int_idealwork))
    call DSYEV(jobz, uplo, n, a, n, w, work, int_idealwork, info)
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

  End Subroutine dble_dsyev

  !!* Complex eigensolver for a Hermitian matrix
  Subroutine cmplx_cheev(a,w,uplo,jobz)
    complex(rsp), intent(inout) :: a(:,:)
    real(rsp), intent(out) :: w(:)
    character, intent(in) :: uplo
    character, intent(in) :: jobz
    real(rsp), allocatable :: rwork(:)
    complex(rsp), allocatable :: work(:)
    integer n, info
    integer :: int_idealwork
    complex(rsp) :: idealwork(1)
    @:ASSERT(uplo == 'u' .or. uplo == 'U' .or. uplo == 'l' .or. uplo == 'L')
    @:ASSERT(jobz == 'n' .or. jobz == 'N' .or. jobz == 'v' .or. jobz == 'V')
    @:ASSERT(all(shape(a)==size(w,dim=1)))
    n=size(a,dim=1)
    @:ASSERT(n>0)
    allocate(rwork(3*n-2))
    call CHEEV(jobz, uplo, n, a, n, w, idealwork, -1, rwork, info)
    if (info/=0) then
      call error("Failue in CHEEV to determine optimum workspace")
    endif
    int_idealwork=floor(real(idealwork(1)))
    allocate(work(int_idealwork))
    call CHEEV(jobz, uplo, n, a, n, w, work, int_idealwork, rwork, info)
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

  End Subroutine cmplx_cheev

  !!* Double complex eigensolver for a Hermitian matrix
  Subroutine dblecmplx_zheev(a,w,uplo,jobz)
    complex(rdp), intent(inout) :: a(:,:)
    real(rdp), intent(out) :: w(:)
    character, intent(in) :: uplo
    character, intent(in) :: jobz
    real(rdp), allocatable :: rwork(:)
    complex(rdp), allocatable :: work(:)
    integer n, info
    integer :: int_idealwork
    complex(rdp) :: idealwork(1)
    @:ASSERT(uplo == 'u' .or. uplo == 'U' .or. uplo == 'l' .or. uplo == 'L')
    @:ASSERT(jobz == 'n' .or. jobz == 'N' .or. jobz == 'v' .or. jobz == 'V')
    @:ASSERT(all(shape(a)==size(w,dim=1)))
    n=size(a,dim=1)
    @:ASSERT(n>0)
    allocate(rwork(3*n-2))
    call ZHEEV(jobz, uplo, n, a, n, w, idealwork, -1, rwork, info)
    if (info/=0) then
      call error("Failue in ZHEEV to determine optimum workspace")
    endif
    int_idealwork=floor(real(idealwork(1)))
    allocate(work(int_idealwork))
    call ZHEEV(jobz, uplo, n, a, n, w, work, int_idealwork, rwork, info)
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

  End Subroutine dblecmplx_zheev


  !!* Real eigensolver for generalized symmetric matrix problem
  Subroutine real_ssygv(a,b,w,uplo,jobz,itype)
    real(rsp), intent(inout) :: a(:,:)
    real(rsp), intent(inout) :: b(:,:)
    real(rsp), intent(out) :: w(:)
    character, intent(in) :: uplo
    character, intent(in) :: jobz
    integer, optional, intent(in) :: itype
    real(rsp), allocatable :: work(:)
    integer n, info, iitype
    integer :: int_idealwork
    real(rsp) :: idealwork(1)
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
    call SSYGV(iitype, jobz, uplo, n, a, n, b, n, w, idealwork, -1, info)
    if (info/=0) then
       call error("Failue in SSYGV to determine optimum workspace")
    endif
    int_idealwork=floor(idealwork(1))
    allocate(work(int_idealwork))
    call SSYGV(iitype, jobz, uplo, n, a, n, b, n, w, work, int_idealwork, info)
    if (info/=0) then
       if (info<0) then
99160 format ('Failure in diagonalisation routine ssygv,', &
          & ' illegal argument at position ',i6)
          write(error_string, 99160) info
          call error(error_string)
       else if (info <= n) then
99170 format ('Failure in diagonalisation routine ssygv,', &
          & ' diagonal element ',i6,' did not converge to zero.')
          write(error_string, 99170) info
          call error(error_string)
       else
99180 format ('Failure in diagonalisation routine ssygv,', &
          & ' non-positive definite overlap! Minor ',i6,' responsible.')
          write(error_string, 99180) info - n
          call error(error_string)
       endif
    endif

  End Subroutine real_ssygv

  !!* Double precision eigensolver for generalized symmetric matrix problem
  Subroutine dble_dsygv(a,b,w,uplo,jobz,itype)
    real(rdp), intent(inout) :: a(:,:)
    real(rdp), intent(inout) :: b(:,:)
    real(rdp), intent(out) :: w(:)
    character, intent(in) :: uplo
    character, intent(in) :: jobz
    integer, optional, intent(in) :: itype
    real(rdp), allocatable :: work(:)
    integer n, info, iitype
    integer :: int_idealwork
    real(rdp) :: idealwork(1)
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
    call DSYGV(iitype, jobz, uplo, n, a, n, b, n, w, idealwork, -1, info)
    if (info/=0) then
       call error("Failue in DSYGV to determine optimum workspace")
    endif
    int_idealwork=floor(idealwork(1))
    allocate(work(int_idealwork))
    call DSYGV(iitype, jobz, uplo, n, a, n, b, n, w, work, int_idealwork, info)
    if (info/=0) then
       if (info<0) then
99190 format ('Failure in diagonalisation routine dsygv,', &
          & ' illegal argument at position ',i6)
          write(error_string, 99190) info
          call error(error_string)
       else if (info <= n) then
99200 format ('Failure in diagonalisation routine dsygv,', &
          & ' diagonal element ',i6,' did not converge to zero.')
          write(error_string, 99200) info
          call error(error_string)
       else
99210 format ('Failure in diagonalisation routine dsygv,', &
          & ' non-positive definite overlap! Minor ',i6,' responsible.')
          write(error_string, 99210) info - n
          call error(error_string)
       endif
    endif

  End Subroutine dble_dsygv

  !!* Complex eigensolver for generalized Hermitian matrix problem
  Subroutine cmplx_chegv(a,b,w,uplo,jobz,itype)
    complex(rsp), intent(inout) :: a(:,:)
    complex(rsp), intent(inout) :: b(:,:)
    real(rsp), intent(out) :: w(:)
    character, intent(in) :: uplo
    character, intent(in) :: jobz
    integer, optional, intent(in) :: itype
    complex(rsp), allocatable :: work(:)
    real(rsp), allocatable :: rwork(:)
    integer n, info, iitype
    integer ::  NB, LWKOPT

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
!   Bugfix for fault in allocation query in zhegv. Should be removed if
!   lapack gets fixed
!
!   In zhegv/chegv the reference lapack reads :
!
!         IF( INFO.NE.0 ) THEN
!            CALL XERBLA( 'ZHEGV ', -INFO )
!            RETURN
!         END IF
!
!   But in dsygv/ssygv instead the test reads :
!
!         IF( INFO.NE.0 ) THEN
!            CALL XERBLA( 'DSYGV ', -INFO )
!            RETURN
!         ELSE IF( LQUERY ) THEN
!            RETURN
!         END IF
!
!
!   Hence the complex routines attempt to solve the eigenproblem even when
!   called as a workspace query.
    NB = ILAENV( 1, 'CHETRD', UPLO, N, -1, -1, -1 )
    LWKOPT = ( NB+1 )*N
!   end bug fix
    allocate(work(LWKOPT))
    ! A*x = (lambda)*B*x upper triangles to be used
    call CHEGV(iitype, 'V', 'L', n, a, n, b, n, w, work, LWKOPT, rwork, info)
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

  End Subroutine cmplx_chegv

  !!* Double complex eigensolver for generalized Hermitian matrix problem
  Subroutine dblecmplx_zhegv(a,b,w,uplo,jobz,itype)
    complex(rdp), intent(inout) :: a(:,:)
    complex(rdp), intent(inout) :: b(:,:)
    real(rdp), intent(out) :: w(:)
    character, intent(in) :: uplo
    character, intent(in) :: jobz
    integer, optional, intent(in) :: itype
    complex(rdp), allocatable :: work(:)
    real(rdp), allocatable :: rwork(:)
    integer n, info, iitype
    integer ::  NB, LWKOPT

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
!   Bugfix for fault in allocation query in zhegv. Should be removed if
!   lapack gets fixed
!
!   In zhegv/zhegv the reference lapack reads :
!
!         IF( INFO.NE.0 ) THEN
!            CALL XERBLA( 'ZHEGV ', -INFO )
!            RETURN
!         END IF
!
!   But in dsygv/ssygv instead the test reads :
!
!         IF( INFO.NE.0 ) THEN
!            CALL XERBLA( 'DSYGV ', -INFO )
!            RETURN
!         ELSE IF( LQUERY ) THEN
!            RETURN
!         END IF
!
!
!   Hence the complex routines attempt to solve the eigenproblem even when
!   called as a workspace query.
    NB = ILAENV( 1, 'CHETRD', UPLO, N, -1, -1, -1 )
    LWKOPT = ( NB+1 )*N
!   end bug fix
    allocate(work(LWKOPT))
    ! A*x = (lambda)*B*x upper triangles to be used
    call ZHEGV(iitype, 'V', 'L', n, a, n, b, n, w, work, LWKOPT, rwork, info)
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

  End Subroutine dblecmplx_zhegv

  !!* Real eigensolver for generalized symmetric matrix problem - divide and
  !!* conquer
  Subroutine real_ssygvd(a,b,w,uplo,jobz,itype)
    real(rsp), intent(inout) :: a(:,:)
    real(rsp), intent(inout) :: b(:,:)
    real(rsp), intent(out) :: w(:)
    character, intent(in) :: uplo
    character, intent(in) :: jobz
    integer, optional, intent(in) :: itype
    real(rsp), allocatable :: work(:)
    integer n, info, iitype
    integer :: int_idealwork, iidealwork(1)
    integer, allocatable :: iwork(:)
    real(rsp) :: idealwork(1)
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
    call SSYGVD(iitype, jobz, uplo, n, a, n, b, n, w, idealwork, -1, &
         & iidealwork, -1, info)
    if (info/=0) then
       call error("Failue in SSYGVD to determine optimum workspace")
    endif
    int_idealwork=floor(idealwork(1))
    allocate(work(int_idealwork))
    allocate(iwork(iidealwork(1)))
    call SSYGVD(iitype, jobz, uplo, n, a, n, b, n, w, work, int_idealwork, &
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

  End Subroutine real_ssygvd

  !!* Double precision eigensolver for generalized symmetric matrix problem
  !!* divide and conquer
  Subroutine dble_dsygvd(a,b,w,uplo,jobz,itype)
    real(rdp), intent(inout) :: a(:,:)
    real(rdp), intent(inout) ::  b(:,:)
    real(rdp), intent(out) :: w(:)
    character, intent(in) :: uplo
    character, intent(in) :: jobz
    integer, optional, intent(in) :: itype
    real(rdp), allocatable :: work(:)
    integer n, info, iitype
    integer :: int_idealwork, iidealwork(1)
    integer, allocatable :: iwork(:)
    real(rdp) :: idealwork(1)
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
    call DSYGVD(iitype, jobz, uplo, n, a, n, b, n, w, idealwork, -1, &
        & iidealwork, -1, info)
    if (info/=0) then
      call error("Failue in DSYGVD to determine optimum workspace")
    endif
    int_idealwork=floor(idealwork(1))
    allocate(work(int_idealwork))
    allocate(iwork(iidealwork(1)))
    call DSYGVD(iitype, jobz, uplo, n, a, n, b, n, w, work, int_idealwork, &
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

  End Subroutine dble_dsygvd

  !!* Complex eigensolver for generalized Hermitian matrix problem divide and
  !!* conquer
  Subroutine cmplx_chegvd(a,b,w,uplo,jobz,itype)
    complex(rsp), intent(inout) :: a(:,:)
    complex(rsp), intent(inout) :: b(:,:)
    real(rsp), intent(out) :: w(:)
    character, intent(in) :: uplo
    character, intent(in) :: jobz
    integer, optional, intent(in) :: itype
    complex(rsp), allocatable :: work(:)
    real(rsp), allocatable :: rwork(:)
    integer n, info, iitype
    integer :: int_idealwork, int_ridealwork, iidealwork(1)
    integer, allocatable :: iwork(:)
    complex(rsp) :: idealwork(1)
    real(rsp) :: ridealwork(1)
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
    call CHEGVD(iitype, jobz, uplo, n, a, n, b, n, w, idealwork, -1, &
        & ridealwork, -1, iidealwork, -1, info)
    if (info/=0) then
      call error("Failue in CHEGVD to determine optimum workspace")
    endif
    int_idealwork=floor(real(idealwork(1)))
    int_ridealwork=floor(ridealwork(1))
    allocate(work(int_idealwork))
    allocate(rwork(int_ridealwork))
    allocate(iwork(iidealwork(1)))
    call CHEGVD(iitype, jobz, uplo, n, a, n, b, n, w, work, int_idealwork,&
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

  End Subroutine cmplx_chegvd

  !!* Double complex eigensolver for generalized Hermitian matrix problem
  !!* divide and conquer
  Subroutine dblecmplx_zhegvd(a,b,w,uplo,jobz,itype)
    complex(rdp), intent(inout) :: a(:,:)
    complex(rdp), intent(inout) :: b(:,:)
    real(rdp), intent(out) :: w(:)
    character, intent(in) :: uplo
    character, intent(in) :: jobz
    integer, optional, intent(in) :: itype
    complex(rdp), allocatable :: work(:)
    real(rdp), allocatable :: rwork(:)
    integer n, info, iitype
    integer :: int_idealwork, int_ridealwork, iidealwork(1)
    integer, allocatable :: iwork(:)
    complex(rdp) :: idealwork(1)
    real(rdp) :: ridealwork(1)
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
    call ZHEGVD(iitype, jobz, uplo, n, a, n, b, n, w, idealwork, -1, &
        & ridealwork, -1, iidealwork, -1, info)
    if (info/=0) then
      call error("Failue in ZHEGVD to determine optimum workspace")
    endif
    int_idealwork=floor(real(idealwork(1)))
    int_ridealwork=floor(ridealwork(1))
    allocate(work(int_idealwork))
    allocate(rwork(int_ridealwork))
    allocate(iwork(iidealwork(1)))
    call ZHEGVD(iitype, jobz, uplo, n, a, n, b, n, w, work, int_idealwork, &
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

  End Subroutine dblecmplx_zhegvd

  !!* Real eigensolver for generalized symmetric matrix problem -
  !!* Relatively Robust Representation, optionally use the subspace form if w
  !!* is smaller than the  size of a and b, then only the first n
  !!* eigenvalues/eigenvectors are found
  !!* This version re-uses a triangle of a matrix (saving an additional
  !!* allocation that was in the previous version)
  !!* @author B. Hourahine, based in part on deMon routine from T. Heine
  subroutine real_ssygvr(a,b,w,uplo,jobz,itype,ilIn,iuIn)
    real(rsp), intent(inout) :: a(:,:)
    real(rsp), intent(inout) :: b(:,:)
    real(rsp), intent(out) :: w(:)
    character, intent(in) :: uplo
    character, intent(in) :: jobz
    integer, optional, intent(in) :: itype
    integer, optional, intent(in) :: ilIn
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
      call SSYEVR(jobz,'I',uplo,n,a,size(a,dim=1),vl,vu,il,iu,abstol,m,&
          & w,b,size(b,dim=1),isuppz,tmpWork,-1,tmpIWork,-1,info)
    else
      call SSYEVR(jobz,'A',uplo,n,a,size(a,dim=1),vl,vu,il,iu,abstol,m,&
          & w,b,size(b,dim=1),isuppz,tmpWork,-1,tmpIWork,-1,info)
    end if

    if (info/=0) then
      call error("Failue in SSYGVR to determine optimum workspace")
    endif
    lwork = floor(tmpWork(1))
    liwork = floor(real(tmpIWork(1)))
    allocate(work(lwork))
    allocate(iwork(liwork))

    ! Form a Cholesky factorization of B.
    call SPOTRF( uplo, n, b, n, info )
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

  !!* Double precision  eigensolver for generalized symmetric matrix problem -
  !!* Relatively Robust Representation, optionally use the subspace form if w
  !!* is smaller than the  size of a and b, then only the first n
  !!* eigenvalues/eigenvectors are found
  !!* This version re-uses a triangle of a matrix (saving an additional
  !!* allocation that was in the previous version)
  !!* @author B. Hourahine, based in part on deMon routine from T. Heine
  subroutine dble_dsygvr(a,b,w,uplo,jobz,itype,ilIn,iuIn)
    real(rdp), intent(inout) :: a(:,:)
    real(rdp), intent(inout) :: b(:,:)
    real(rdp), intent(out) :: w(:)
    character, intent(in) :: uplo
    character, intent(in) :: jobz
    integer, optional, intent(in) :: itype
    integer, optional, intent(in) :: ilIn
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
      call DSYEVR(jobz,'I',uplo,n,a,size(a,dim=1),vl,vu,il,iu,abstol,m,&
          & w,b,size(b,dim=1),isuppz,tmpWork,-1,tmpIWork,-1,info)
    else
      call DSYEVR(jobz,'A',uplo,n,a,size(a,dim=1),vl,vu,il,iu,abstol,m,&
          & w,b,size(b,dim=1),isuppz,tmpWork,-1,tmpIWork,-1,info)
    end if

    if (info/=0) then
      call error("Failue in DSYGVR to determine optimum workspace")
    endif
    lwork = floor(tmpWork(1))
    liwork = floor(real(tmpIWork(1)))
    allocate(work(lwork))
    allocate(iwork(liwork))

    ! Form a Cholesky factorization of B.
    call DPOTRF( uplo, n, b, n, info )
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

  !!* Complex eigensolver for generalized symmetric matrix problem -
  !!* Relatively Robust Representation, optionally use the subspace form if w
  !!* is smaller than the  size of a and b, then only the first n
  !!* eigenvalues/eigenvectors are found
  !!* @author B. Hourahine, based in part on deMon routine from T. Heine
  subroutine cmplx_chegvr(a,b,w,uplo,jobz,itype,ilIn,iuIn)
    complex(rsp), intent(inout) :: a(:,:)
    complex(rsp), intent(inout) :: b(:,:)
    real(rsp), intent(out) :: w(:)
    character, intent(in) :: uplo
    character, intent(in) :: jobz
    integer, optional, intent(in) :: itype
    integer, optional, intent(in) :: ilIn
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
      call CHEEVR(jobz,'I',uplo,n,a,size(a,dim=1),vl,vu,il,iu,abstol,m,&
          & w,b,size(b,dim=1),isuppz,tmpWork,-1,tmpRwork,-1,tmpIWork,-1,info)
    else
      call CHEEVR(jobz,'A',uplo,n,a,size(a,dim=1),vl,vu,il,iu,abstol,m,&
          & w,b,size(b,dim=1),isuppz,tmpWork,-1,tmpRwork,-1,tmpIWork,-1,info)
    end if

    if (info/=0) then
      call error("Failue in CHEEVR to determine optimum workspace")
    endif
    lwork = floor(real(tmpWork(1)))
    liwork = floor(real(tmpIWork(1)))
    lrwork = floor(tmpRwork(1))
    allocate(work(lwork))
    allocate(iwork(liwork))
    allocate(rwork(lrwork))

    ! Form a Cholesky factorization of B.
    call CPOTRF( uplo, n, b, n, info )
    if( info /= 0 ) then
      info = n + info
99400 format ('Failure in diagonalisation routine CHEGVR,', &
          & ' unable to complete Cholesky factorization of B ',i6)
      write(error_string, 99400) info
      call error(error_string)
    end if
    ! Transform problem to standard eigenvalue problem and solve.
    call CHEGST( iitype, uplo, n, a, n, b, n, info )

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
      call CHEEVR( jobz, 'I', uplo, n, a, n, vl, vu, il, iu, &
          & abstol, m, w, b, n, isuppz, work, lwork, rwork, lrwork, iwork, &
          & liwork, info )
    else
      call CHEEVR( jobz, 'A', uplo, n, a, n, vl, vu, il, iu, &
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
        call CTRSM('Left',uplo_new,trans,'Non-unit',n,neig, &
            & cmplx(1.0,0.0,rsp),A,n, B,n)
      else if( iitype == 3 ) then
        ! For B*A*x=(lambda)*x;
        ! backtransform eigenvectors: x = L*y or U'*y               !'
        if( upper ) then
          trans = 'C'
        else
          trans = 'N'
        end if
        call CTRMM('Left',uplo_new,trans,'Non-unit',n,neig, &
            & cmplx(1.0,0.0,rsp),a,n, b,n)
      end if
      do ii = 1,m
        a( 1:n, ii) = b( 1:n, ii )
      end do
      a( 1:n, m+1:n ) = cmplx(0.0,0.0,rsp)
    end if

  end subroutine cmplx_chegvr

  !!* Double complex eigensolver for generalized symmetric matrix problem -
  !!* Relatively Robust Representation, optionally use the subspace form if w
  !!* is smaller than the  size of a and b, then only the first n
  !!* eigenvalues/eigenvectors are found
  !!* @author B. Hourahine, based in part on deMon routine from T. Heine
  subroutine dblecmplx_zhegvr(a,b,w,uplo,jobz,itype,ilIn,iuIn)
    complex(rdp), intent(inout) :: a(:,:)
    complex(rdp), intent(inout) :: b(:,:)
    real(rdp), intent(out) :: w(:)
    character, intent(in) :: uplo
    character, intent(in) :: jobz
    integer, optional, intent(in) :: itype
    integer, optional, intent(in) :: ilIn
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
      call ZHEEVR(jobz,'I',uplo,n,a,size(a,dim=1),vl,vu,il,iu,abstol,m,&
          & w,b,size(b,dim=1),isuppz,tmpWork,-1,tmpRwork,-1,tmpIWork,-1,info)
    else
      call ZHEEVR(jobz,'A',uplo,n,a,size(a,dim=1),vl,vu,il,iu,abstol,m,&
          & w,b,size(b,dim=1),isuppz,tmpWork,-1,tmpRwork,-1,tmpIWork,-1,info)
    end if

    if (info/=0) then
      call error("Failue in ZHEEVR to determine optimum workspace")
    endif
    lwork = floor(real(tmpWork(1)))
    liwork = floor(real(tmpIWork(1)))
    lrwork = floor(tmpRwork(1))
    allocate(work(lwork))
    allocate(iwork(liwork))
    allocate(rwork(lrwork))

    ! Form a Cholesky factorization of B.
    call ZPOTRF( uplo, n, b, n, info )
    if( info /= 0 ) then
      info = n + info
99400 format ('Failure in diagonalisation routine ZHEEVR,', &
          & ' unable to complete Cholesky factorization of B ',i6)
      write(error_string, 99400) info
      call error(error_string)
    end if
    ! Transform problem to standard eigenvalue problem and solve.
    call ZHEGST( iitype, uplo, n, a, n, b, n, info )

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
      call ZHEEVR( jobz, 'I', uplo, n, a, n, vl, vu, il, iu, &
          & abstol, m, w, b, n, isuppz, work, lwork, rwork, lrwork, iwork, &
          & liwork, info )
    else
      call ZHEEVR( jobz, 'A', uplo, n, a, n, vl, vu, il, iu, &
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
        call ZTRSM('Left',uplo_new,trans,'Non-unit',n,neig, &
            & cmplx(1.0,0.0,rdp),A,n, B,n)
      else if( iitype == 3 ) then
        ! For B*A*x=(lambda)*x;
        ! backtransform eigenvectors: x = L*y or U'*y               !'
        if( upper ) then
          trans = 'C'
        else
          trans = 'N'
        end if
        call ZTRMM('Left',uplo_new,trans,'Non-unit',n,neig, &
            & cmplx(1.0,0.0,rdp),a,n, b,n)
      end if
      do ii = 1,m
        a( 1:n, ii) = b( 1:n, ii )
      end do
      a( 1:n, m+1:n ) = cmplx(0.0,0.0,rdp)
    end if

  end subroutine dblecmplx_zhegvr



    !!* simple single precision banded matrix eigensolver
  subroutine real_ssbgv(ab, bb, w, uplo, z)
    real(rsp), intent(inout) :: ab(:,:)
    real(rsp), intent(inout) :: bb(:,:)
    real(rsp), intent(out) :: w(:)
    character, intent(in) :: uplo
    real(rsp), optional, intent(out) :: z(:,:)

    real(rsp), allocatable :: work(:)
    integer :: n, ka, kb, ldab, ldbb, ldz, info
    character :: jobz
    real(rsp) :: zTmp(1,1)

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
      call SSBGV( jobz,uplo,n,ka,kb,ab,ldab,bb,ldbb,w,z,ldz,work,info )
    else
      call SSBGV( jobz,uplo,n,ka,kb,ab,ldab,bb,ldbb,w,zTmp,ldz,work,info )
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

  !!* Simple double precision banded matrix eigen solver
  subroutine dble_dsbgv(ab, bb, w, uplo, z)
    real(rdp), intent(inout) :: ab(:,:)
    real(rdp), intent(inout) :: bb(:,:)
    real(rdp), intent(out) :: w(:)
    character, intent(in) :: uplo
    real(rdp), optional, intent(out) :: z(:,:)

    real(rdp), allocatable :: work(:)
    integer :: n, ka, kb, ldab, ldbb, ldz, info
    character :: jobz
    real(rdp) :: zTmp(1,1)

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
      call DSBGV( jobz,uplo,n,ka,kb,ab,ldab,bb,ldbb,w,z,ldz,work,info )
    else
      call DSBGV( jobz,uplo,n,ka,kb,ab,ldab,bb,ldbb,w,zTmp,ldz,work,info )
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

  !!* Simple complex precision banded matrix eigen solver
  subroutine cmplx_chbgv(ab, bb, w, uplo, z)
    complex(rsp), intent(inout) :: ab(:,:)
    complex(rsp), intent(inout) :: bb(:,:)
    real(rsp), intent(out) :: w(:)
    character, intent(in) :: uplo
    complex(rsp), optional, intent(out) :: z(:,:)

    complex(rsp), allocatable :: work(:)
    real(rsp), allocatable :: rwork(:)
    integer :: n, ka, kb, ldab, ldbb, ldz, info
    character :: jobz
    complex(rsp) :: zTmp(1,1)

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
      call CHBGV( jobz,uplo,n,ka,kb,ab,ldab,bb,ldbb,w,z,ldz,work,rwork,info )
    else
      call CHBGV( jobz,uplo,n,ka,kb,ab,ldab,bb,ldbb,w,zTmp,ldz,work,rwork, &
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

  !!* Simple double complex precision banded matrix eigen solver
  subroutine dblecmplx_zhbgv(ab, bb, w, uplo, z)
    complex(rdp), intent(inout) :: ab(:,:)
    complex(rdp), intent(inout) :: bb(:,:)
    real(rdp), intent(out) :: w(:)
    character, intent(in) :: uplo
    complex(rdp), optional, intent(out) :: z(:,:)

    complex(rdp), allocatable :: work(:)
    real(rdp), allocatable :: rwork(:)
    integer :: n, ka, kb, ldab, ldbb, ldz, info
    character :: jobz
    complex(rdp) :: zTmp(1,1)

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
      call ZHBGV( jobz,uplo,n,ka,kb,ab,ldab,bb,ldbb,w,z,ldz,work,rwork,info )
    else
      call ZHBGV( jobz,uplo,n,ka,kb,ab,ldab,bb,ldbb,w,zTmp,ldz,work,rwork, &
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

end module eigensolver
