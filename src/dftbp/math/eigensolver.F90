!--------------------------------------------------------------------------------------------------!
!  DFTB+: general package for performing fast atomistic simulations                                !
!  Copyright (C) 2006 - 2025  DFTB+ developers group                                               !
!                                                                                                  !
!  See the LICENSE file for terms of usage and distribution.                                       !
!--------------------------------------------------------------------------------------------------!

#:include 'common.fypp'
#:include "error.fypp"

!> Contains F90 wrapper functions for some commonly used lapack calls needed in the code.
!> Contains some fixes for lapack 3.0 bugs, if this gets corrected in lapack 4.x they should be
!> removed.
module dftbp_math_eigensolver
  use dftbp_common_accuracy, only : rdp, rsp
  use dftbp_extlibs_lapack, only : cheev, cheevd, cheevr, chegst, chegv, chegvd, cpotrf, ctrmm,&
      & ctrsm, dgeev, dlamch, dpotrf, dsyev, dsyevd, dsyevr, dsygst, dsygv, dsygvd, dtrmm, dtrsm,&
      & sgeev, slamch, spotrf, ssyev, ssyevd, ssyevr, ssygst, ssygv, ssygvd, strmm, strsm, zheev,&
      & zheevd, zheevr, zhegst, zhegv, zhegvd, zpotrf, ztrmm, ztrsm
  use dftbp_io_message, only : error
#:if WITH_MAGMA
  use dftbp_extlibs_magma, only : magmaf_cheevd_m, magmaf_chegvd_m, magmaf_dsyevd_m,&
      & magmaf_dsygvd_m, magmaf_ssyevd_m, magmaf_ssygvd_m, magmaf_zheevd_m, magmaf_zhegvd_m
#:endif
  implicit none (type, external)

  private
  public :: heev, heevd, heevr, hegv, hegvd, hegvr, geev
#:if WITH_MAGMA
  public :: magmaHeevd, magmaHegvd
#:endif


  !> QR eigensolver for a symmetric/hermitian matrix
  !! Caveat: the matrix a is overwritten
  interface heev
  #:for SUFFIX in ["real", "dreal", "complex", "dcomplex"]
    module procedure heev_${SUFFIX}$
  #:endfor
  end interface heev


  !> Divide and conquer eigensolver for a symmetric/hermitian matrix
  !! Caveat: the matrix a is overwritten
  interface heevd
  #:for SUFFIX in ["real", "dreal", "complex", "dcomplex"]
    module procedure heevd_${SUFFIX}$
  #:endfor
  end interface heevd


  !> Relatively robust eigensolver for a symmetric/hermitian matrix
  !! Caveat: the matrix a is overwritten
  interface heevr
  #:for SUFFIX in ["real", "dreal", "complex", "dcomplex"]
    module procedure heevr_${SUFFIX}$
  #:endfor
  end interface heevr


  !> QR eigensolver for a symmetric/hermitian generalized matrix problem
  !! caveat: the matrix a is overwritten
  !! caveat: the matrix b is overwritten with Cholesky factorization
  interface hegv
  #:for SUFFIX in ["real", "dreal", "complex", "dcomplex"]
    module procedure hegv_${SUFFIX}$
  #:endfor
  end interface hegv


  !> Divide and conquer eigensolver for a symmetric/hermitian generalized matrix
  !! conquer eigensolver
  !! caveat: the matrix a is overwritten
  !! caveat: the matrix b is overwritten with Cholesky factorization
  interface hegvd
  #:for SUFFIX in ["real", "dreal", "complex", "dcomplex"]
    module procedure hegvd_${SUFFIX}$
  #:endfor
  end interface hegvd


  !> Relatively robust eigensolver for a symmetric/hermitian generalized matrix problem.
  !!
  !! Using the lapack ! relatively robust representation solver, based on the SYGV source. If the
  !requested number of ! eigenvalues is lower than the size of H/S suspace mode is used (optionally
  !the range can be set ! using il and ul) to return the lowest eigenvalues/vectors of number
  !! size(w)
  interface hegvr
  #:for SUFFIX in ["real", "dreal", "complex", "dcomplex"]
    module procedure hegvr_${SUFFIX}$
  #:endfor
  end interface


  !> QR eigensolver for a general matrix
  interface geev
  #:for SUFFIX in ["real", "dreal"]
    module procedure geev_${SUFFIX}$
  #:endfor
  end interface geev


#:if WITH_MAGMA

  !> Divide and conquer MAGMA GPU eigensolver
  interface magmaHeevd
  #:for SUFFIX in ["real", "dreal", "complex", "dcomplex"]
    module procedure magmaHeevd_${SUFFIX}$
  #:endfor
  end interface magmaHeevd

  !> Divide and conquer MAGMA GPU generalised eigensolver
  interface magmaHegvd
  #:for SUFFIX in ["real", "dreal", "complex", "dcomplex"]
    module procedure magmaHegvd_${SUFFIX}$
  #:endfor
  end interface magmaHegvd

#:endif

contains


#:for TYPE, KIND, SUFFIX, LAPACK_ROUTINE in &
    & [("real", "rsp", "real", "ssyev"),&
    & ("real", "rdp", "dreal", "dsyev"),&
    & ("complex", "rsp", "complex", "cheev"),&
    & ("complex", "rdp", "dcomplex", "zheev")]

  !> Real eigensolver for a symmetric matrix
  subroutine heev_${SUFFIX}$(a, w, uplo, jobz, info)

    !> contains the matrix for the solver, returns as eigenvectors if requested (matrix always
    !> overwritten on return anyway)
    ${TYPE}$(${KIND}$), intent(inout) :: a(:,:)

    !> eigenvalues
    real(${KIND}$), intent(out) :: w(:)

    !> upper or lower triangle of the matrix
    character, intent(in) :: uplo

    !> compute eigenvalues 'N' or eigenvalues and eigenvectors 'V'
    character, intent(in) :: jobz

    !> if present and info /= 0 job is to be terminated by the calling routine
    integer, optional, intent(out) :: info

    ${TYPE}$(${KIND}$), allocatable :: work(:)
    ${TYPE}$(${KIND}$) :: workDummy(1)
    integer :: workSize
  #:if TYPE == "complex"
    real(${KIND}$), allocatable :: rwork(:)
  #:endif
    integer n, info_, iStep
    character(len=100) :: errorMsg

    @:ASSERT(uplo == 'u' .or. uplo == 'U' .or. uplo == 'l' .or. uplo == 'L')
    @:ASSERT(jobz == 'n' .or. jobz == 'N' .or. jobz == 'v' .or. jobz == 'V')
    @:ASSERT(all(shape(a) == size(w, dim=1)))
    n = size(a, dim=1)
    @:ASSERT(n > 0)

  #:if TYPE == "complex"
    allocate(rwork(3 * n - 2))
  #:endif

    errorGuard: block
      iStep = 1
    #:if TYPE == "real"
      call ${LAPACK_ROUTINE}$(jobz, uplo, n, a, n, w, workDummy, -1, info_)
    #:else
      call ${LAPACK_ROUTINE}$(jobz, uplo, n, a, n, w, workDummy, -1, rwork, info_)
    #:endif
      if (info_ /= 0) exit errorGuard

      iStep = 2
      workSize = nint(real(workDummy(1), kind=${KIND}$))
      allocate(work(workSize))
    #:if TYPE == "real"
      call ${LAPACK_ROUTINE}$(jobz, uplo, n, a, n, w, work, workSize, info_)
    #:else
      call ${LAPACK_ROUTINE}$(jobz, uplo, n, a, n, w, work, workSize, rwork, info_)
    #:endif
    end block errorGuard

    if (present(info)) info = info_
    if (info_ == 0 .or. present(info)) return

    select case (iStep)
    case (1)
      call error("Failure in LAPACK routine ${LAPACK_ROUTINE}$ to determine optimum workspace")
    case (2)
      if (info_ < 0) then
        write(errorMsg, "(a, i0)") "Failure LAPACK routine ${LAPACK_ROUTINE}$, illegal argument at&
            & position ", -info_
      else
        write(errorMsg, "(a, i0, a)") "Failure in LAPACK routine ${LAPACK_ROUTINE}$, ", info_,&
            & " off-diagonal elements did not converge to zero."
      end if
    end select
    call error(errorMsg)

  end subroutine heev_${SUFFIX}$

#:endfor


#:for TYPE, KIND, SUFFIX, LAPACK_ROUTINE in &
    & [("real", "rsp", "real", "ssyevd"),&
    & ("real", "rdp", "dreal", "dsyevd"),&
    & ("complex", "rsp", "complex", "cheevd"),&
    & ("complex", "rdp", "dcomplex", "zheevd")]

  !> Real eigensolver for a symmetric matrix
  subroutine heevd_${SUFFIX}$(a, w, uplo, jobz, info)

    !> contains the matrix for the solver, returns as eigenvectors if requested (matrix always
    !> overwritten on return anyway)
    ${TYPE}$(${KIND}$), intent(inout) :: a(:,:)

    !> eigenvalues
    real(${KIND}$), intent(out) :: w(:)

    !> upper or lower triangle of the matrix
    character, intent(in) :: uplo

    !> compute eigenvalues 'N' or eigenvalues and eigenvectors 'V'
    character, intent(in) :: jobz

    !> if present and info /= 0 job is to be terminated by the calling routine
    integer, optional, intent(out) :: info

    ${TYPE}$(${KIND}$), allocatable :: work(:)
    ${TYPE}$(${KIND}$) :: workDummy(1)
    integer :: workSize
    integer, allocatable :: iwork(:)
    integer :: iworkDummy(1)
    integer :: iworkSize
  #:if TYPE == "complex"
    real(${KIND}$), allocatable :: rwork(:)
    real(${KIND}$) :: rworkDummy(1)
    integer :: rworkSize
  #:endif
    integer n, info_, iStep
    character(len=100) :: errorMsg

    @:ASSERT(uplo == 'u' .or. uplo == 'U' .or. uplo == 'l' .or. uplo == 'L')
    @:ASSERT(jobz == 'n' .or. jobz == 'N' .or. jobz == 'v' .or. jobz == 'V')
    @:ASSERT(all(shape(a) == size(w, dim=1)))
    n = size(a, dim=1)
    @:ASSERT(n > 0)

    errorGuard: block
      iStep = 1
    #:if TYPE == "real"
      call ${LAPACK_ROUTINE}$(jobz, uplo, n, a, n, w, workDummy, -1, iworkDummy, -1, info_)
    #:else
      call ${LAPACK_ROUTINE}$(jobz, uplo, n, a, n, w, workDummy, -1, rworkDummy, -1, iworkDummy,&
          & -1, info_)
    #:endif
      if (info_ /= 0) exit errorGuard

      iStep = 2
      workSize = nint(real(workDummy(1), kind=${KIND}$))
      allocate(work(workSize))
      iworkSize = iworkDummy(1)
      allocate(iwork(iworkSize))
    #:if TYPE == "real"
      call ${LAPACK_ROUTINE}$(jobz, uplo, n, a, n, w, work, workSize, iwork, iworkSize, info_)
    #:else
      rworkSize = nint(rworkDummy(1))
      allocate(rwork(rworkSize))
      call ${LAPACK_ROUTINE}$(jobz, uplo, n, a, n, w, work, workSize, rwork, rworkSize, iwork,&
          & iworkSize, info_)
    #:endif
    end block errorGuard

    if (present(info)) info = info_
    if (info_ == 0 .or. present(info)) return

    select case (iStep)
    case (1)
      call error("Failure in LAPACK routine ${LAPACK_ROUTINE}$ to determine optimum workspace")
    case (2)
      if (info_ < 0) then
        write(errorMsg, "(a, i0)") "Failure in LAPACK routine ${LAPACK_ROUTINE}$, illegal argument&
            & at position ", -info_
      else
        write(errorMsg, "(a, i0, a)") "Failure in LAPACK routine ${LAPACK_ROUTINE}$, ", info_,&
            & " off-diagonal elements did not converge to zero."
        end if
    end select
    call error(errorMsg)

  end subroutine heevd_${SUFFIX}$

#:endfor


#:for TYPE, KIND, SUFFIX, LAPACK_HEEVR, LAPACK_LAMCH in &
    & [("real", "rsp", "real", "ssyevr", "slamch"),&
    & ("real", "rdp", "dreal", "dsyevr", "dlamch"),&
    & ("complex", "rsp", "complex", "cheevr", "slamch"),&
    & ("complex", "rdp", "dcomplex", "zheevr", "dlamch")]

  !> Relatively robust eigensolver for symmetric/hermitian matrix problem
  !> Representation, optionally use the subspace form if w is smaller than the size of a and b, then
  !> only the first n eigenvalues/eigenvectors are found.
  subroutine heevr_${SUFFIX}$(a, w, uplo, jobz, ilIn, iuIn, info)

    !> Contains the matrix for the solver, returns eigenvectors if requested (matrix always
    !! overwritten on return anyway)
    ${TYPE}$(${KIND}$), intent(inout) :: a(:,:)

    !> Eigenvalues on return
    real(${KIND}$), intent(out) :: w(:)

    !> Upper or lower triangle of both matrices
    character, intent(in) :: uplo

    !> Compute eigenvalues 'N' or eigenvalues and eigenvectors 'V'
    character, intent(in) :: jobz

    !> Lower range of eigenstates
    integer, optional, intent(in) :: ilIn

    !> Upper range of eigenstates
    integer, optional, intent(in) :: iuIn

    !> if present and info /= 0 job is to be terminated by the calling routine
    integer, optional, intent(out) :: info

    ${TYPE}$(${KIND}$), allocatable :: work(:)
    ${TYPE}$(${KIND}$) :: workDummy(1)
    integer :: workSize
    integer, allocatable :: iwork(:)
    integer :: iworkDummy(1)
    integer :: iworkSize
  #:if TYPE == "complex"
    real(${KIND}$), allocatable :: rwork(:)
    real(${KIND}$) :: rworkDummy(1)
    integer :: rworkSize
  #:endif
    ${TYPE}$(${KIND}$), allocatable :: z(:,:)
    real(${KIND}$) :: abstol, vl, vu
    integer :: n, info_, m, il, iu, ldz, iStep
    integer, allocatable :: isuppz(:)
    logical :: subspace
    character :: range
    character(len=100) :: errorMsg

    n = size(a, dim=1)

    @:ASSERT(n > 0)
    @:ASSERT(uplo == 'u' .or. uplo == 'U' .or. uplo == 'l' .or. uplo == 'L')
    @:ASSERT(jobz == 'n' .or. jobz == 'N' .or. jobz == 'v' .or. jobz == 'V')
    @:ASSERT(present(ilIn) .eqv. present(iuIn))

    subspace = (size(w) < n)
    if (subspace) then
      range = "I"
      if (present(ilIn)) then
        @:ASSERT(ilIn <= iuIn)
        @:ASSERT(ilIn > 0)
        @:ASSERT(n >= iuIn)
        @:ASSERT(size(w) == (iuIn - ilIn + 1))
        il = ilIn
        iu = iuIn
      else
        il = 1
        iu = size(w)
      end if
    else
      range = "A"
      if (present(ilIn)) then
        @:ASSERT(ilIn == 1 .and. iuIn == n)
      end if
      il = 1
      iu = n
    end if

    if (jobz == 'v' .or. jobz == 'V') then
      allocate(z(n, iu - il + 1))
      ldz = n
    else
      allocate(z(1, 1))
      ldz = 1
    end if

    allocate(isuppz(2 * n))
    abstol = ${LAPACK_LAMCH}$( 'Safe minimum' )

    errorGuard: block
      iStep = 1
      #:if TYPE == "real"
        call ${LAPACK_HEEVR}$(jobz, range, uplo, n, a, n, vl, vu, il, iu, abstol, m, w, z, ldz,&
            & isuppz, workDummy, -1, iworkDummy, -1, info_)
      #:else
        call ${LAPACK_HEEVR}$(jobz, range, uplo, n, a, n, vl, vu, il, iu, abstol, m, w, z, ldz,&
            & isuppz, workDummy, -1, rworkDummy, -1, iworkDummy, -1, info_)
      #:endif
      if (info_ /= 0) exit errorGuard

      iStep = 2
      workSize = nint(real(workDummy(1), kind=${KIND}$))
      allocate(work(workSize))
      iworkSize = iworkDummy(1)
      allocate(iwork(iworkSize))
    #:if TYPE == "real"
      call ${LAPACK_HEEVR}$(jobz, range, uplo, n, a, n, vl, vu, il, iu, abstol, m, w, z, ldz,&
          & isuppz, work, workSize, iwork, iworkSize, info_)
    #:else
      rworkSize = nint(rworkDummy(1))
      allocate(rwork(rworkSize))
      call ${LAPACK_HEEVR}$(jobz, range, uplo, n, a, n, vl, vu, il, iu, abstol, m, w, z, ldz,&
          & isuppz, work, workSize, rwork, rworkSize, iwork, iworkSize, info_)
    #:endif
      if (info_ /= 0) exit errorGuard

      iStep = 3
      a(:,:) = 0.0_${KIND}$
      if (jobz == 'v' .or. jobz == 'V') then
        a(:n, :iu-il+1) = z
      end if
    end block errorGuard

    if (present(info)) info = info_
    if (info_ == 0 .or. present(info)) return

    select case (iStep)
    case (1)
      call error("Failure in LAPACK routine ${LAPACK_HEEVR}$ to determine optimum workspace")
    case (2)
      if (info_ < 0) then
        write(errorMsg, "(a, i0)") "Failure in LAPACK routine ${LAPACK_HEEVR}$, illegal argument at&
            & position ", -info_
      else
        write(errorMsg, "(a, i0, a)") "Failure in LAPACK routine ${LAPACK_HEEVR}$, internal error"
      end if
    end select
    call error(errorMsg)

  end subroutine heevr_${SUFFIX}$
#:endfor


#:for TYPE, KIND, SUFFIX, LAPACK_ROUTINE in &
    & [("real", "rsp", "real", "ssygv"),&
    & ("real", "rdp", "dreal", "dsygv"),&
    & ("complex", "rsp", "complex", "chegv"),&
    & ("complex", "rdp", "dcomplex", "zhegv")]

  !> QR eigensolver for generalized symmetric/hermitian matrix problem
  subroutine hegv_${SUFFIX}$(a, b, w, uplo, jobz, itype, info)

    !> contains the matrix for the solver, returns eigenvectors if requested (matrix always
    !> overwritten on return anyway)
    ${TYPE}$(${KIND}$), intent(inout) :: a(:,:)

    !> contains the second matrix for the solver (overwritten by Cholesky factorization)
    ${TYPE}$(${KIND}$), intent(inout) :: b(:,:)

    !> eigenvalues
    real(${KIND}$), intent(out) :: w(:)

    !> upper or lower triangle of both matrices
    character, intent(in) :: uplo

    !> compute eigenvalues 'N' or eigenvalues and eigenvectors 'V'
    character, intent(in) :: jobz

    !> specifies the problem type to be solved 1:A*x=(lambda)*B*x, 2:A*B*x=(lambda)*x,
    !> 3:B*A*x=(lambda)*x default is 1
    integer, optional, intent(in) :: itype

    !> if present and info/=0 job is to be terminated by the calling routine
    integer, optional, intent(out) :: info

    ${TYPE}$(${KIND}$), allocatable :: work(:)
    ${TYPE}$(${KIND}$) :: workDummy(1)
    integer :: workSize
  #:if TYPE == "complex"
    real(${KIND}$), allocatable :: rwork(:)
    real(${KIND}$) :: rworkDummy(1)
  #:endif
    integer n, lda, info_, iitype, ldb, iStep
    character(len=100) :: errorMsg

    @:ASSERT(uplo == 'u' .or. uplo == 'U' .or. uplo == 'l' .or. uplo == 'L')
    @:ASSERT(jobz == 'n' .or. jobz == 'N' .or. jobz == 'v' .or. jobz == 'V')
    n = size(a, dim=2)
    @:ASSERT(n > 0)
    @:ASSERT(all(shape(b) >= n))
    @:ASSERT(size(w) >= n)
    lda = size(a, dim=1)
    @:ASSERT(lda >= n)
    ldb = size(b, dim=1)
    if (present(itype)) then
      iitype = itype
    else
      iitype = 1
    end if
    @:ASSERT(iitype >= 1 .and. iitype <= 3)

    errorGuard: block
      iStep = 1
    #:if TYPE == "real"
      call ${LAPACK_ROUTINE}$(iitype, jobz, uplo, n, a, lda, b, ldb, w, workDummy, -1, info_)
    #:else
      call ${LAPACK_ROUTINE}$(iitype, jobz, uplo, n, a, lda, b, ldb, w, workDummy, -1,&
          & rworkDummy, info_)
    #:endif
      if (info_/=0) exit errorGuard

      iStep = 2
      workSize = nint(real(workDummy(1), kind=${KIND}$))
      allocate(work(workSize))
    #:if TYPE == "real"
      call ${LAPACK_ROUTINE}$(iitype, jobz, uplo, n, a, lda, b, ldb, w, work, workSize, info_)
    #:else
      allocate(rwork(3 * n - 2))
      call ${LAPACK_ROUTINE}$(iitype, jobz, uplo, n, a, lda, b, ldb, w, work, workSize, rwork,&
          & info_)
    #:endif
    end block errorGuard

    if (present(info)) info = info_
    if (info_ == 0 .or. present(info)) return

    select case (iStep)
    case (1)
      call error("Failure in LAPACK routine ${LAPACK_ROUTINE}$ to determine optimum workspace")
    case (2)
      if (info_ < 0) then
        write(errorMsg, "(a, i0)") "Failure in LAPACK routine ${LAPACK_ROUTINE}$, illegal argument&
            & at position ", -info_
      else if (info_ <= n) then
        write(errorMsg, "(a, i0, a)") "Failure in LAPACK routine ${LAPACK_ROUTINE}$, ", info_,&
            & " off diagonal elements did not converge to zero."
      else
        write(errorMsg, "(a, i0, a)") "Failure in LAPACK routine ${LAPACK_ROUTINE}$, non-positive&
            & definite overlap, minor ", info_ - n, " responsible."
      end if
    end select
    call error(errorMsg)

  end subroutine hegv_${SUFFIX}$

#:endfor



#:for TYPE, KIND, SUFFIX, LAPACK_ROUTINE in&
    & [("real", "rsp", "real", "ssygvd"),&
    & ("real", "rdp", "dreal", "dsygvd"),&
    & ("complex", "rsp", "complex", "chegvd"),&
    & ("complex", "rdp", "dcomplex", "zhegvd")]

  !> Real eigensolver for generalized symmetric matrix problem - divide and conquer
  subroutine hegvd_${SUFFIX}$(a, b, w, uplo, jobz, itype, info)

    !> contains the matrix for the solver, returns eigenvectors if requested (matrix always
    !> overwritten on return anyway)
    ${TYPE}$(${KIND}$), intent(inout) :: a(:,:)

    !> contains the second matrix for the solver (overwritten by Cholesky factorization)
    ${TYPE}$(${KIND}$), intent(inout) :: b(:,:)

    !> eigenvalues
    real(${KIND}$), intent(out) :: w(:)

    !> upper or lower triangle of the matrix
    character, intent(in) :: uplo

    !> compute eigenvalues 'N' or eigenvalues and eigenvectors 'V'
    character, intent(in) :: jobz

    !> optional specifies the problem type to be solved 1:A*x=(lambda)*B*x, 2:A*B*x=(lambda)*x,
    !> 3:B*A*x=(lambda)*x default is 1
    integer, optional, intent(in) :: itype

    !> if present and info/=0 job is to be terminated by the calling routine
    integer, optional, intent(out) :: info

    ${TYPE}$(${KIND}$), allocatable :: work(:)
    ${TYPE}$(${KIND}$) :: workDummy(1)
    integer :: workSize
    integer, allocatable :: iwork(:)
    integer :: iworkDummy(1)
    integer :: iworkSize
  #:if TYPE == "complex"
    real(${KIND}$), allocatable :: rwork(:)
    real(${KIND}$) :: rworkDummy(1)
    integer :: rworkSize
  #:endif
    integer n, info_, iitype, iStep
    character(len=100) :: errorMsg

    @:ASSERT(uplo == 'u' .or. uplo == 'U' .or. uplo == 'l' .or. uplo == 'L')
    @:ASSERT(jobz == 'n' .or. jobz == 'N' .or. jobz == 'v' .or. jobz == 'V')
    @:ASSERT(all(shape(a) == shape(b)))
    @:ASSERT(all(shape(a) == size(w, dim=1)))
    n = size(a, dim=1)
    @:ASSERT(n > 0)
    if (present(itype)) then
      iitype = itype
    else
      iitype = 1
    end if
    @:ASSERT(iitype >= 1 .and. iitype <= 3)

    errorGuard: block
      iStep = 1
    #:if TYPE == "real"
      call ${LAPACK_ROUTINE}$(iitype, jobz, uplo, n, a, n, b, n, w, workDummy, -1, &
         & iworkDummy, -1, info_)
    #:else
      call ${LAPACK_ROUTINE}$(iitype, jobz, uplo, n, a, n, b, n, w, workDummy, -1,&
         & rworkDummy, -1, iworkDummy, -1, info_)
    #:endif
      if (info_/=0) exit errorGuard

      iStep = 2
      worksize = nint(real(workDummy(1), kind=${KIND}$))
      allocate(work(workSize))
      iworkSize = iworkDummy(1)
      allocate(iwork(iworkSize))
    #:if TYPE == "real"
      call ${LAPACK_ROUTINE}$(iitype, jobz, uplo, n, a, n, b, n, w, work, workSize, iwork,&
          & iworkSize, info_)
    #:else
      rworkSize = nint(rworkDummy(1))
      allocate(rwork(rworkSize))
      call ${LAPACK_ROUTINE}$(iitype, jobz, uplo, n, a, n, b, n, w, work, workSize, rwork,&
          & rworkSize, iwork, iworkSize, info_)
    #:endif
    end block errorGuard

    if (present(info)) info = info_
    if (info_ == 0 .or. present(info)) return

    select case (iStep)
    case (1)
      call error("Failure in LAPACK routine ${LAPACK_ROUTINE}$ to determine optimum workspace")
    case (2)
      if (info_ < 0) then
        write(errorMsg, "(a, i0)") "Failure in LAPACK routine ${LAPACK_ROUTINE}$, illegal argument&
            & at position ", -info_
      else if (info_ <= n) then
        write(errorMsg, "(a, i0, a)") "Failure in LAPACK routine ${LAPACK_ROUTINE}$, ", info_,&
            & " off diagonal elements did not converge to zero."
      else
        write(errorMsg, "(a, i0, a)") "Failure in LAPACK routine ${LAPACK_ROUTINE}$, non-positive&
            & definite overlap, minor ", info_ - n, " responsible."
      end if
    end select
    call error(errorMsg)

  end subroutine hegvd_${SUFFIX}$

#:endfor


#:for TYPE, KIND, SUFFIX, CALLS in &
    & [("real", "rsp", "real", ("ssyevr", "spotrf", "ssygst", "strsm", "strmm", "slamch")),&
    &  ("real", "rdp", "dreal", ("dsyevr", "dpotrf", "dsygst", "dtrsm", "dtrmm", "dlamch")),&
    &  ("complex", "rsp", "complex", ("cheevr", "cpotrf", "chegst", "ctrsm", "ctrmm", "slamch")),&
    &  ("complex", "rdp", "dcomplex", ("zheevr", "zpotrf", "zhegst", "ztrsm", "ztrmm", "dlamch"))]

    #:set LAPACK_HEEVR, LAPACK_POTRF, LAPACK_HEGST, LAPACK_TRSM, LAPACK_TRMM, LAPACK_LAMCH = CALLS

  !> Relatively robust eigensolver for generalized symmetric matrix problem.
  !!
  !! Representation, optionally use the subspace form if w is smaller than the size of a and b, then
  !! only the first n eigenvalues/eigenvectors are found.
  !! This version re-uses a triangle of a matrix (saving an additional allocation that was in the
  !! previous version).
  !! Based in part on deMon routine from T. Heine
  subroutine hegvr_${SUFFIX}$(a, b, w, uplo, jobz, itype, ilIn, iuIn, info)

    !> contains the matrix for the solver, returns eigenvectors if requested (matrix always
    !! overwritten on return anyway)
    ${TYPE}$(${KIND}$), intent(inout) :: a(:,:)

    !> contains the second matrix for the solver (overwritten by Cholesky factorization)
    ${TYPE}$(${KIND}$), intent(inout) :: b(:,:)

    !> eigenvalues
    real(${KIND}$), intent(out) :: w(:)

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

    !> if present and info/=0 job is to be terminated by the calling routine
    integer, optional, intent(out) :: info

    ${TYPE}$(${KIND}$), allocatable :: work(:)
    ${TYPE}$(${KIND}$) :: workDummy(1)
    integer :: workSize
    integer, allocatable :: iwork(:)
    integer :: iworkDummy(1)
    integer :: iworkSize
  #:if TYPE == "complex"
    real(${KIND}$), allocatable :: rwork(:)
    real(${KIND}$) :: rworkDummy(1)
    integer :: rworkSize
  #:endif
    ${TYPE}$(${KIND}$), allocatable :: tmpChole(:)
    integer, allocatable :: isuppz(:)
    integer :: n, info_, iitype, iStep
    integer :: m, neig
    real(${KIND}$) :: abstol
    logical :: wantz, upper
    character :: trans, range
    real(${KIND}$) :: vl, vu
    integer :: il, iu
    logical :: subspace
    integer :: ii, jj
    character :: uploNew
    character(len=100) :: errorMsg

    character, parameter :: transNormal = "N"
  #:if TYPE == "real"
    character, parameter :: transAdjoint = "T"
    real(${KIND}$), parameter :: one = 1.0_${KIND}$
  #:else
    character, parameter :: transAdjoint = "C"
    complex(${KIND}$), parameter :: one =  (1.0_${KIND}$, 0.0_${KIND}$)
  #:endif

    n = size(a, dim=1)

    @:ASSERT(n > 0)
    @:ASSERT(uplo == 'u' .or. uplo == 'U' .or. uplo == 'l' .or. uplo == 'L')
    @:ASSERT(jobz == 'n' .or. jobz == 'N' .or. jobz == 'v' .or. jobz == 'V')
    @:ASSERT(all(shape(a) == shape(b)))
    @:ASSERT(present(ilIn) .eqv. present(iuIn))

    subspace = (size(w) < n)
    if (subspace) then
      if (present(ilIn)) then
        @:ASSERT(ilIn <= iuIn)
        @:ASSERT(ilIn > 0)
        @:ASSERT(n >= iuIn)
        @:ASSERT(size(w, dim=1) == (iuIn - ilIn + 1))
        il = ilIn
        iu = iuIn
      else
        il = 1
        iu = size(w, dim=1)
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
    @:ASSERT(iitype >= 1 .and. iitype <= 3)

    allocate(isuppz(2 * n))
    allocate(tmpChole(size(a, dim=1)))

    wantz = (jobz == 'V' .or. jobz == 'v')
    upper = (uplo == 'U' .or. uplo == 'u')
    abstol = ${LAPACK_LAMCH}$('Safe minimum')
    if (subspace) then
      range = 'I'
    else
      range = 'A'
    end if

    errorGuard: block
      ! Determine the optimal workspace sizes
      iStep = 1
    #:if TYPE == "real"
      call ${LAPACK_HEEVR}$(jobz, range, uplo, n, a, size(a, dim=1), vl, vu, il, iu, abstol, m,&
          & w, b, size(b, dim=1), isuppz, workDummy, -1, iworkDummy, -1, info_)
    #:else
      call ${LAPACK_HEEVR}$(jobz, range, uplo, n, a, size(a, dim=1), vl, vu, il, iu, abstol, m,&
          & w, b, size(b, dim=1), isuppz, workDummy, -1, rworkDummy, -1, iworkDummy, -1, info_)
    #:endif
      if (info_ /= 0) exit errorGuard

      ! Form a Cholesky factorization of B.
      iStep = 2
      workSize = nint(real(workDummy(1), kind=${KIND}$))
      allocate(work(workSize))
      iworkSize = iworkDummy(1)
      allocate(iwork(iworkSize))
    #:if TYPE == "complex"
      rworkSize = nint(rworkDummy(1))
      allocate(rwork(rworkSize))
    #:endif
      call ${LAPACK_POTRF}$(uplo, n, b, n, info_)
      if (info_ /= 0) then
        info_ = n + info_
        exit errorGuard
      end if

      ! Transform problem to standard eigenvalue
      iStep = 3
      call ${LAPACK_HEGST}$(iitype, uplo, n, a, n, b, n, info_)
      if (info_ /= 0) exit errorGuard

      ! Solve standard eigenvalue problem
      iStep = 4
      if (wantz) then
        ! Save Cholesky factor in the other triangle of H and tmpChole
        do ii = 1, n
          tmpChole(ii) = b(ii, ii)
        end do
        if (upper) then
          do jj = 1, n
            do ii = jj+1, n
              a(ii, jj) = #{if TYPE == "real"}# b(jj, ii) #{else}# conjg(b(jj, ii)) #{endif}#
            end do
          end do
        else
          do jj = 1, n
            do ii = 1, jj - 1
              a(ii, jj) = #{if TYPE == "real"}# b(jj, ii) #{else}# conjg(b(jj, ii)) #{endif}#
            end do
          end do
        end if
      end if
    #:if TYPE == "real"
      call ${LAPACK_HEEVR}$(jobz, range, uplo, n, a, n, vl, vu, il, iu, abstol, m, w, b, n,&
          & isuppz, work, workSize, iwork, iworkSize, info_)
    #:else
      call ${LAPACK_HEEVR}$(jobz, range, uplo, n, a, n, vl, vu, il, iu, abstol, m, w, b, n,&
          & isuppz, work, workSize, rwork, rworkSize, iwork, iworkSize, info_)
    #:endif
      if (info_ /= 0) exit errorGuard

      ! Backtransform eigenvectors to the original problem.
      iStep = 5
      if (wantz) then

        do ii = 1, n
          a(ii,ii) = tmpChole(ii)
        end do

        if (upper) then
          uploNew = 'L'
          upper = .false.
        else
          uploNew = 'U'
          upper = .true.
        end if

        neig = n
        if (iitype == 1 .or. iitype == 2) then
          ! For A*x=(lambda)*B*x and A*B*x=(lambda)*x;
          ! backtransform eigenvectors: x = inv(L)'*y or inv(U)*y          !'
          if (upper) then
            trans = transNormal
          else
            trans = transAdjoint
          end if
          call ${LAPACK_TRSM}$('Left', uploNew, trans, 'Non-unit', n, neig, one, A, n, B, n)
        else if (iitype == 3) then
          ! For B*A*x=(lambda)*x;
          ! backtransform eigenvectors: x = L*y or U'*y     !'
          if (upper) then
            trans = transAdjoint
          else
            trans = transNormal
          end if
          call ${LAPACK_TRMM}$('Left', uploNew, trans, 'Non-unit', n, neig, one, a, n, b, n)
        end if
        do ii = 1,m
          a(1:n, ii) = b(1:n, ii)
        end do
        a(1:n, m+1:n) = 0.0
      end if
    end block errorGuard

    if (present(info)) info = info_
    if (info_ == 0 .or. present(info)) return

    select case (iStep)
    case (1)
      errorMsg = "Failure in LAPACK routine ${LAPACK_HEEVR}$ to determine optimum workspace"
    case (2)
      write(errorMsg, "(a, i0)") "Failure in LAPACK routine ${LAPACK_POTRF}$, unable to complete&
          & Cholesky factorization of B, info: ", info_
    case (3)
      write(errorMsg, "(a, i0)") "Failure in ${LAPACK_HEGST}$ to transform to standard form,&
          & info: ", info_
    case (4)
      write(errorMsg, "(a, i0)") "Failure in LAPACK routine ${LAPACK_HEEVR}$ to solve the&
          & eigenvalue, problem, info: ", info_
    end select
    call error(errorMsg)

  end subroutine hegvr_${SUFFIX}$

#:endfor


#:for KIND, SUFFIX, LAPACK_ROUTINE in [("rsp", "real", "sgeev"), ("rdp", "dreal", "dgeev")]

  !> QR general matrix eigensolver
  subroutine geev_${SUFFIX}$(a, wr, wi, vl, vr, info)

    !> Matrix, overwritten on exit
    real(${KIND}$), intent(inout) :: a(:,:)

    !> Real part of eigenvalues
    real(${KIND}$), intent(out) :: wr(:)

    !> Imaginary part of eigenvalues
    real(${KIND}$), intent(out) :: wi(:)

    !> Left eigenvectors
    real(${KIND}$), intent(out), optional :: vl(:,:)

    !> Right eigenvectors
    real(${KIND}$), intent(out), optional :: vr(:,:)

    !> if present and info/=0 job is to be terminated by the calling routine
    integer, optional, intent(out) :: info

    real(${KIND}$), allocatable :: work(:)
    real(${KIND}$) :: workDummy(1)
    integer :: workSize
    integer :: n, lda, info_, ldvl, ldvr, iStep
    character :: jobvl, jobvr
    character(len=100) :: errorMsg

    ! If no eigenvectors requested, need a dummy array for lapack call
    real(${KIND}$) :: dummyvl(1,1), dummyvr(1,1)

    lda = size(a, dim=1)
    n = size(a, dim=2)

    @:ASSERT(n > 0)
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

    errorGuard: block
      iStep = 1
      if (jobvl == 'V' .and. jobvr == 'V') then
        call ${LAPACK_ROUTINE}$(jobvl, jobvr, n, a, lda, wr, wi, vl, ldvl, vr, ldvr,&
            & workDummy, -1, info_)
      else if (jobvl == 'V' .and. jobvr == 'N') then
        call ${LAPACK_ROUTINE}$(jobvl, jobvr, n, a, lda, wr, wi, vl, ldvl, dummyvr, ldvr,&
            & workDummy, -1, info_)
      else if (jobvl == 'N' .and. jobvr == 'V') then
        call ${LAPACK_ROUTINE}$(jobvl, jobvr, n, a, lda, wr, wi, dummyvl, ldvl, vr, ldvr,&
            & workDummy, -1, info_)
      else if (jobvl == 'N' .and. jobvr == 'N') then
        call ${LAPACK_ROUTINE}$(jobvl, jobvr, n, a, lda, wr, wi, dummyvl, ldvl, dummyvr, ldvr,&
            & workDummy, -1, info_)
      end if
      if (info_ /= 0) exit errorGuard

      iStep = 2
      workSize = nint(workDummy(1))
      allocate(work(workSize))
      if (jobvl == 'V' .and. jobvr == 'V') then
        call ${LAPACK_ROUTINE}$(jobvl, jobvr, n, a, lda, wr, wi, vl, ldvl, vr, ldvr, work,&
            & workSize, info_)
      else if (jobvl == 'V' .and. jobvr == 'N') then
        call ${LAPACK_ROUTINE}$(jobvl, jobvr, n, a, lda, wr, wi, vl, ldvl, dummyvr, ldvr, work,&
            & workSize, info_)
      else if (jobvl == 'N' .and. jobvr == 'V') then
        call ${LAPACK_ROUTINE}$(jobvl, jobvr, n, a, lda, wr, wi, dummyvl, ldvl, vr, ldvr, work,&
            & workSize, info_)
      else if (jobvl == 'N' .and. jobvr == 'N') then
        call ${LAPACK_ROUTINE}$(jobvl, jobvr, n, a, lda, wr, wi, dummyvl, ldvl, dummyvr, ldvr,&
            & work, workSize, info_)
      end if
    end block errorGuard

    if (present(info)) info = info_
    if (info_ == 0 .or. present(info)) return

    select case (iStep)
    case (1)
      errorMsg = "Failure in LAPACK routine ${LAPACK_ROUTINE}$ to determine optimum workspace"
    case (2)
      if (info_ < 0) then
        write(errorMsg, "(a, i0)") "Failure in diagonalisation routine ${LAPACK_ROUTINE}$, illegal&
            & argument at position ", -info_
      else
        write(errorMsg, "(a, i0, a)") "Failure in diagonalisation routine ${LAPACK_ROUTINE}$,&
            & QR algorithm failed to calculate ", info_, " eigenvalues."
      end if
    end select
    call error(errorMsg)

  end subroutine geev_${SUFFIX}$

#:endfor


#:if WITH_MAGMA

#:for TYPE, KIND, SUFFIX, MAGMA_ROUTINE in &
    & [("real", "rsp", "real", "magmaf_ssyevd_m"),&
    &  ("real", "rdp", "dreal", "magmaf_dsyevd_m"),&
    &  ("complex", "rsp", "complex", "magmaf_cheevd_m"),&
    &  ("complex", "rdp", "dcomplex", "magmaf_zheevd_m")]


  !> Eigensolution for symmetric/hermitian matrices on GPU(s)
  subroutine magmaHeevd_${SUFFIX}$(ngpus, a, w, uplo, jobz, info)

    !> Number of GPUs to use
    integer, intent(in) :: ngpus

    !> contains the matrix for the solver, returns eigenvectors if requested (matrix always
    !! overwritten on return anyway)
    ${TYPE}$(${KIND}$), intent(inout) :: a(:,:)

    !> eigenvalues
    real(${KIND}$), intent(out) :: w(:)

    !> upper or lower triangle of the matrix
    character, intent(in) :: uplo

    !> compute eigenvalues 'N' or eigenvalues and eigenvectors 'V'
    character, intent(in) :: jobz

    !> if present and info/=0 job is to be terminated by the calling routine
    integer, optional, intent(out) :: info

   ${TYPE}$(${KIND}$), allocatable :: work(:)
   ${TYPE}$(${KIND}$) :: workDummy(1)
   integer :: workSize
   integer, allocatable :: iwork(:)
   integer :: iworkDummy(1)
   integer :: iworkSize
 #:if TYPE == 'complex'
    real(${KIND}$), allocatable :: rwork(:)
    real(${KIND}$) :: rworkDummy(1)
    integer :: rworkSize
  #:endif
    integer :: n, info_, iStep
    character(len=100) :: errorMsg

    @:ASSERT(uplo == 'u' .or. uplo == 'U' .or. uplo == 'l' .or. uplo == 'L')
    @:ASSERT(jobz == 'n' .or. jobz == 'N' .or. jobz == 'v' .or. jobz == 'V')
    @:ASSERT(all(shape(a) == size(w)))
    n = size(a, dim=1)
    @:ASSERT(n > 0)

    errorGuard: block
      iStep = 1
    #:if TYPE == 'real'
      call ${MAGMA_ROUTINE}$(ngpus, jobz, uplo, n, a, n, w, workDummy, -1, iworkDummy, -1, info_)
    #:else
      call ${MAGMA_ROUTINE}$(ngpus, jobz, uplo, n, a, n, w, workDummy, -1, rworkDummy, -1,&
          & iworkDummy, -1, info)
    #:endif
      if (info_ /= 0) exit errorGuard

      iStep = 2
      workSize = nint(real(workDummy(1), kind=${KIND}$))
      allocate(work(workSize))
      iworkSize = iworkDummy(1)
      allocate(iwork(iworkSize))
    #:if TYPE == "real"
      call ${MAGMA_ROUTINE}$(ngpus, jobz, uplo, n, a, n, w, work, workSize, iwork, iworkSize, info_)
    #:else
      rworkSize = nint(rwork(1))
      allocate(rwork(rworkSize))
      call ${MAGMA_ROUTINE}$(ngpus, jobz, uplo, n, a, n, w, work, workSize, rwork, rworkSize,&
          & iwork, iworkSize, info_)
    #:endif
    end block errorGuard

    if (present(info)) info = info_
    if (info_ == 0 .or. present(info)) return

    select case (iStep)
    case (1)
      errorMsg = "Failure in MAGMA routine ${MAGMA_ROUTINE}$ to determine optimum workspace"
    case (2)
      if (info_ < 0) then
        write(errorMsg, "(a, i0)") "Failure in MAGMA routine ${MAGMA_ROUTINE}$, illegal argument&
            & at position ", -info_
      else if (info_ <= n) then
        write(errorMsg, "(a, i0, a)") "Failure in MAGMA routine ${MAGMA_ROUTINE}$ during&
            & diagonalization, info: ", info_
      end if
    end select
    call error(errorMsg)

  end subroutine magmaHeevd_${SUFFIX}$

#:endfor


#:for TYPE, KIND, SUFFIX, MAGMA_ROUTINE in &
    & [("real", "rsp", "real", "magmaf_ssygvd_m"),&
    &  ("real", "rdp", "dreal", "magmaf_dsygvd_m"),&
    &  ("complex", "rsp", "complex", "magmaf_chegvd_m"),&
    &  ("complex", "rdp", "dcomplex", "magmaf_zhegvd_m")]

  !> Generalised eigensolution for symmetric/hermitian matrices on GPU(s)
  subroutine magmaHegvd_${SUFFIX}$(ngpus, a, b, w, uplo, jobz, itype, info)

    !> Number of GPUs to use
    integer, intent(in) :: ngpus

    !> contains the matrix for the solver, returns eigenvectors if requested (matrix always
    !> overwritten on return anyway)
    ${TYPE}$(${KIND}$), intent(inout) :: a(:,:)

    !> contains the second matrix for the solver (overwritten by Cholesky factorization)
    ${TYPE}$(${KIND}$), intent(inout) :: b(:,:)

    !> eigenvalues
    real(${KIND}$), intent(out) :: w(:)

    !> upper or lower triangle of the matrix
    character, intent(in) :: uplo

    !> compute eigenvalues 'N' or eigenvalues and eigenvectors 'V'
    character, intent(in) :: jobz

    !> optional specifies the problem type to be solved 1:A*x=(lambda)*B*x, 2:A*B*x=(lambda)*x,
    !> 3:B*A*x=(lambda)*x default is 1
    integer, optional, intent(in) :: itype

    !> if present and info/=0 job is to be terminated by the calling routine
    integer, optional, intent(out) :: info

    ${TYPE}$(${KIND}$), allocatable :: work(:)
    ${TYPE}$(${KIND}$) :: workDummy(1)
    integer :: workSize
    integer, allocatable :: iwork(:)
    integer :: iworkDummy(1)
    integer :: iworkSize
  #:if TYPE == "complex"
    real(${KIND}$), allocatable :: rwork(:)
    real(${KIND}$) :: rworkDummy(1)
    integer :: rworkSize
  #:endif
    integer :: n, info_, iitype, iStep
    character(len=100) :: errorMsg

    @:ASSERT(uplo == 'u' .or. uplo == 'U' .or. uplo == 'l' .or. uplo == 'L')
    @:ASSERT(jobz == 'n' .or. jobz == 'N' .or. jobz == 'v' .or. jobz == 'V')
    @:ASSERT(all(shape(a) == shape(b)))
    @:ASSERT(all(shape(a) == size(w, dim=1)))
    n = size(a, dim=1)
    @:ASSERT(n > 0)
    if (present(itype)) then
      iitype = itype
    else
      iitype = 1
    end if
    @:ASSERT(iitype >= 1 .and. iitype <= 3)

    errorGuard: block
      iStep = 1
    #:if TYPE == "real"
      call ${MAGMA_ROUTINE}$(ngpus, iitype, jobz, uplo, n, a, n, b, n, w, workDummy, -1,&
          & iworkDummy, -1, info_)
    #:else
      call ${MAGMA_ROUTINE}$(ngpus, iitype, jobz, uplo, n, a, n, b, n, w, workDummy, -1,&
          & rworkDummy, -1, iworkDummy, -1, info_)
    #:endif
      if (info_ /= 0) exit errorGuard

      iStep = 2
      workSize = nint(real(workDummy(1), kind=${KIND}$))
      allocate(work(workSize))
      iworkSize = iworkDummy(1)
      allocate(iwork(iworkSize))
    #:if TYPE == "real"
      call ${MAGMA_ROUTINE}$(ngpus, iitype, jobz, uplo, n, a, n, b, n, w, work, workSize, iwork,&
          & iworkSize, info_)
    #:else
      rworkSize = nint(rworkDummy(1))
      allocate(rwork(rworkSize))
      call ${MAGMA_ROUTINE}$(ngpus, iitype, jobz, uplo, n, a, n, b, n, w, work, workSize, rwork,&
          & rworkSize, iwork, iworkSize, info_)
    #:endif
    end block errorGuard

    if (present(info)) info = info_
    if (info_ == 0 .or. present(info)) return

    select case (iStep)
    case (1)
      errorMsg = "Failure in routine ${MAGMA_ROUTINE}$ to determine optimum workspace"
    case (2)
      if (info_ < 0) then
        write(errorMsg, "(i0, a)") "Failure in diagonalisation routine ${MAGMA_ROUTINE}$, illegal&
          & argument at position ", -info_
      else if (info_ <= n) then
        write(errorMsg, "(a, i0, a)") "Failure in diagonalisation routine ${MAGMA_ROUTINE}$, ",&
          & info_, " off-diagonal elements did not converge to zero."
      else
        write(errorMsg, "(a, i0, a)") "Failure in diagonalisation routine ${MAGMA_ROUTINE}$,&
          & non-positive definite overlap, minor ", info_ - n, " responsible."
      end if
    end select
    call error(errorMsg)

  end subroutine magmaHegvd_${SUFFIX}$

#:endfor

#:endif

end module dftbp_math_eigensolver
