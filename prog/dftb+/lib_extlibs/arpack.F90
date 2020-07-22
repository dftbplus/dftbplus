!--------------------------------------------------------------------------------------------------!
!  DFTB+: general package for performing fast atomistic simulations                                !
!  Copyright (C) 2006 - 2020  DFTB+ developers group                                               !
!                                                                                                  !
!  See the LICENSE file for terms of usage and distribution.                                       !
!--------------------------------------------------------------------------------------------------!

#:include 'common.fypp'

!> Interfaces for the ARPACK routines needed in DFTB+ (currently for the linear response excited
!> state calculations).
module dftbp_arpack
  use dftbp_accuracy, only : rsp, rdp
#:if not WITH_ARPACK
  use dftbp_message
#:endif
  implicit none
  private

  public :: withArpack, saupd, seupd

#:if WITH_ARPACK

  !> Whether code was built with ARPACK support
  logical, parameter :: withArpack = .true.

#:else

  !> Whether code was built with ARPACK support
  logical, parameter :: withArpack = .false.

  !> Dummy routines, as ARPACK library is not compiled in
  interface saupd
  #:for PREC in [("s"),("d")]
    module procedure ${PREC}$saupd
  #:endfor
  end interface saupd

  !> Dummy routines, as ARPACK library is not compiled in
  interface seupd
  #:for PREC in [("s"),("d")]
    module procedure ${PREC}$seupd
  #:endfor
  end interface

contains

  !> Generates error message, if a stub was called
  subroutine stubError(routineName)
    character(*), intent(in) :: routineName

    call error("Internal error: " // trim(routineName) // "() called in a build without ARPACK&
        & support")

  end subroutine stubError

#:endif

#:if WITH_ARPACK
  !> Wrapper around ARPACK routines ssaupd/dsaupd.
  interface saupd
#:endif
  #:for PREC, LABEL, VTYPE in [("s","single","rsp"),("d","double","rdp")]
  #:if not WITH_ARPACK
    !> Dummy ARPACK routine
  #:endif
    !> ${LABEL}$ precision Arnoldi solver call
    subroutine ${PREC}$saupd(ido, bmat, n, which, nev, tol, resid, ncv, v, ldv, iparam, ipntr,&
        & workd, workl, lworkl, info)
    #:if WITH_ARPACK
      import :: ${VTYPE}$
    #:endif
      integer, intent(inout) :: ido
      character, intent(in) :: bmat
      integer, intent(in) :: n
      character(2), intent(in) :: which
      integer, intent(in) :: nev
      real(${VTYPE}$), intent(in) :: tol
      real(${VTYPE}$), intent(inout) :: resid(n)
      integer, intent(in) :: ncv
      integer, intent(in) :: ldv
      real(${VTYPE}$), intent(out) :: v(ldv, ncv)
      integer, intent(inout) :: iparam(11)
      integer, intent(out) :: ipntr(11)
      real(${VTYPE}$), intent(inout) :: workd(3 * n)
      integer, intent(in) :: lworkl
      real(${VTYPE}$), intent(inout) :: workl(lworkl)
      integer, intent(inout) :: info
     #:if not WITH_ARPACK
      call stubError("${PREC}$saupd")
     #:endif
    end subroutine ${PREC}$saupd
  #:endfor
#:if WITH_ARPACK
  end interface saupd
#:endif

#:if WITH_ARPACK
  !> Wrapper around ARPACK routines sseupd/dseupd.
  interface seupd
#:endif
  #:for PREC, LABEL, VTYPE in [("s","single","rsp"),("d","double","rdp")]
  #:if not WITH_ARPACK
    !> Dummy ARPACK routine
  #:endif
    !> ${LABEL}$ precision return from the results of the solver
    subroutine ${PREC}$seupd(rvec, howmny, sel, d, z, ldz, sigma, bmat, n, which, nev, tol, resid,&
        & ncv, v, ldv, iparam, ipntr, workd, workl, lworkl, info)
    #:if WITH_ARPACK
      import :: ${VTYPE}$
    #:endif
      logical, intent(in) :: rvec
      character, intent(in) :: howmny
      integer, intent(in) :: ncv
      logical, intent(in) :: sel(ncv)
      integer, intent(in) :: nev
      real(${VTYPE}$), intent(out) :: d(nev)
      integer, intent(in) :: ldz
      real(${VTYPE}$), intent(out) :: z(ldz, nev)
      real(${VTYPE}$), intent(in) :: sigma
      character, intent(in) :: bmat
      integer, intent(in) :: n
      character(2), intent(in) :: which
      real(${VTYPE}$), intent(in) :: tol
      real(${VTYPE}$), intent(in) :: resid(n)
      integer, intent(in) :: ldv
      real(${VTYPE}$), intent(inout) :: v(ldv, ncv)
      integer, intent(in) :: iparam(7)
      integer, intent(inout) :: ipntr(11)
      real(${VTYPE}$), intent(inout) :: workd(2 * n)
      integer, intent(in) :: lworkl
      real(${VTYPE}$), intent(inout) :: workl(lworkl)
      integer, intent(inout) :: info
     #:if not WITH_ARPACK
      call stubError("${PREC}$seupd")
     #:endif
    end subroutine ${PREC}$seupd
  #:endfor
#:if WITH_ARPACK
  end interface seupd
#:endif

end module dftbp_arpack
