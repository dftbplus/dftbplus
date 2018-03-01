!--------------------------------------------------------------------------------------------------!
!  DFTB+: general package for performing fast atomistic simulations                                !
!  Copyright (C) 2018  DFTB+ developers group                                                      !
!                                                                                                  !
!  See the LICENSE file for terms of usage and distribution.                                       !
!--------------------------------------------------------------------------------------------------!

#:include 'common.fypp'

!> Interfaces for the ARPACK routines needed in DFTB+ (currently for the linear response excited
!> state calculations).
module arpack
  use accuracy, only : rsp, rdp
  implicit none
  private

  public :: withArpack

#:if WITH_ARPACK

  public :: saupd, seupd


  !> Whether code was built with Arpack support
  logical, parameter :: withArpack = .true.


  !> Wrapper around ARPACK routines ssaupd/dsaupd.
  interface saupd

    !> single precision Arnoldi solver call
    subroutine ssaupd(ido, bmat, n, which, nev, tol, resid, ncv, v, ldv,&
        & iparam, ipntr, workd, workl, lworkl, info)
      import :: rsp
      implicit none
      integer, intent(inout) :: ido
      character, intent(in) :: bmat
      integer, intent(in) :: n
      character(2), intent(in) :: which
      integer, intent(in) :: nev
      real(rsp), intent(in) :: tol
      real(rsp), intent(inout) :: resid(n)
      integer, intent(in) :: ncv
      integer, intent(in) :: ldv
      real(rsp), intent(out) :: v(ldv, ncv)
      integer, intent(inout) :: iparam(11)
      integer, intent(out) :: ipntr(11)
      real(rsp), intent(inout) :: workd(3 * n)
      integer, intent(in) :: lworkl
      real(rsp), intent(inout) :: workl(lworkl)
      integer, intent(inout) :: info
    end subroutine ssaupd


    !> double precision Arnoldi solver call
    subroutine dsaupd(ido, bmat, n, which, nev, tol, resid, ncv, v, ldv,&
        & iparam, ipntr, workd, workl, lworkl, info)
      import :: rdp
      implicit none
      integer, intent(inout) :: ido
      character, intent(in) :: bmat
      integer, intent(in) :: n
      character(2), intent(in) :: which
      integer, intent(in) :: nev
      real(rdp), intent(in) :: tol
      real(rdp), intent(inout) :: resid(n)
      integer, intent(in) :: ncv
      integer, intent(in) :: ldv
      real(rdp), intent(out) :: v(ldv, ncv)
      integer, intent(inout) :: iparam(11)
      integer, intent(out) :: ipntr(11)
      real(rdp), intent(inout) :: workd(3 * n)
      integer, intent(in) :: lworkl
      real(rdp), intent(inout) :: workl(lworkl)
      integer, intent(inout) :: info
    end subroutine dsaupd
  end interface saupd


  !> Wrapper around ARPACK routines sseupd/dseupd.
  interface seupd


    !> single precision return from the results of the solver
    subroutine sseupd(rvec, howmny, sel, d, z, ldz, sigma, bmat, n, which, nev,&
        & tol, resid, ncv, v, ldv, iparam, ipntr, workd, workl, lworkl, info)
      import :: rsp
      logical, intent(in) :: rvec
      character, intent(in) :: howmny
      integer, intent(in) :: ncv
      logical, intent(in) :: sel(ncv)
      integer, intent(in) :: nev
      real(rsp), intent(out) :: d(nev)
      integer, intent(in) :: ldz
      real(rsp), intent(out) :: z(ldz, nev)
      real(rsp), intent(in) :: sigma
      character, intent(in) :: bmat
      integer, intent(in) :: n
      character(2), intent(in) :: which
      real(rsp), intent(in) :: tol
      real(rsp), intent(in) :: resid(n)
      integer, intent(in) :: ldv
      real(rsp), intent(inout) :: v(ldv, ncv)
      integer, intent(in) :: iparam(7)
      integer, intent(inout) :: ipntr(11)
      real(rsp), intent(inout) :: workd(2 * n)
      integer, intent(in) :: lworkl
      real(rsp), intent(inout) :: workl(lworkl)
      integer, intent(inout) :: info
    end subroutine sseupd


    !> double precision return from the results of the solver
    subroutine dseupd(rvec, howmny, sel, d, z, ldz, sigma, bmat, n, which, nev,&
        & tol, resid, ncv, v, ldv, iparam, ipntr, workd, workl, lworkl, info)
      import :: rdp
      logical, intent(in) :: rvec
      character, intent(in) :: howmny
      integer, intent(in) :: ncv
      logical, intent(in) :: sel(ncv)
      integer, intent(in) :: nev
      real(rdp), intent(out) :: d(nev)
      integer, intent(in) :: ldz
      real(rdp), intent(out) :: z(ldz, nev)
      real(rdp), intent(in) :: sigma
      character, intent(in) :: bmat
      integer, intent(in) :: n
      character(2), intent(in) :: which
      real(rdp), intent(in) :: tol
      real(rdp), intent(in) :: resid(n)
      integer, intent(in) :: ldv
      real(rdp), intent(inout) :: v(ldv, ncv)
      integer, intent(in) :: iparam(7)
      integer, intent(inout) :: ipntr(11)
      real(rdp), intent(inout) :: workd(2 * n)
      integer, intent(in) :: lworkl
      real(rdp), intent(inout) :: workl(lworkl)
      integer, intent(inout) :: info
    end subroutine dseupd
  end interface seupd

#:else


  !> Whether code was built with ARPACK support
  logical, parameter :: withArpack = .false.

#:endif

end module arpack
