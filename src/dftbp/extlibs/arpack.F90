!--------------------------------------------------------------------------------------------------!
!  DFTB+: general package for performing fast atomistic simulations                                !
!  Copyright (C) 2006 - 2024  DFTB+ developers group                                               !
!                                                                                                  !
!  See the LICENSE file for terms of usage and distribution.                                       !
!--------------------------------------------------------------------------------------------------!

#:include 'common.fypp'

#:assert not (WITH_ARPACK and INSTANCE_SAFE_BUILD)

#:set ARPACK_PROC_PREFIX = "" if WITH_ARPACK else "module"
#:set PARPACK_PROC_PREFIX = "" if WITH_PARPACK else "module"
#:set PREFIXES_AND_KINDS =  [("s", "rsp"), ("d", "rdp")]


!> Interfaces for the ARPACK routines needed in DFTB+ (currently for the linear response excited
!> state calculations).
module dftbp_extlibs_arpack
  use dftbp_common_accuracy, only : rdp, rsp
#:if not WITH_ARPACK
  use dftbp_io_message, only : error
#:endif
  implicit none

  private
  public :: withArpack, saupd, seupd
  public :: withParpack, psaupd, pseupd


  !> Whether code was built with ARPACK support
  logical, parameter :: withArpack = ${FORTRAN_LOGICAL(WITH_ARPACK)}$

  !> Whether code was built with PARPACK support
  logical, parameter :: withParpack = ${FORTRAN_LOGICAL(WITH_PARPACK)}$


  ! Function overloading to be used within DFTB+
  interface saupd
  #:for PREFIX, _ in PREFIXES_AND_KINDS
    ${ARPACK_PROC_PREFIX}$ procedure ${PREFIX}$saupd
  #:endfor
  end interface

  ! Function overloading to be used within DFTB+
  interface seupd
  #:for PREFIX, _ in PREFIXES_AND_KINDS
    ${ARPACK_PROC_PREFIX}$ procedure ${PREFIX}$seupd
  #:endfor
  end interface

  ! Function overloading to be used within DFTB+
  interface psaupd
  #:for PREFIX, _ in PREFIXES_AND_KINDS
    ${PARPACK_PROC_PREFIX}$ procedure p${PREFIX}$saupd
  #:endfor
  end interface

  ! Function overloading to be used within DFTB+
  interface pseupd
  #:for PREFIX, _ in PREFIXES_AND_KINDS
    ${PARPACK_PROC_PREFIX}$ procedure p${PREFIX}$seupd
  #:endfor
  end interface

  ! Interface definition of the routines
  interface
  #:for PREFIX, KIND in PREFIXES_AND_KINDS

    !> Arnoldi solver call
    ${ARPACK_PROC_PREFIX}$ subroutine ${PREFIX}$saupd(ido, bmat, n, which, nev, tol, resid, ncv,&
        & v, ldv, iparam, ipntr, workd, workl, lworkl, info)
    #:if WITH_ARPACK
      import :: ${KIND}$
    #:endif
      integer, intent(inout) :: ido
      character, intent(in) :: bmat
      integer, intent(in) :: n
      character(2), intent(in) :: which
      integer, intent(in) :: nev
      real(${KIND}$), intent(in) :: tol
      real(${KIND}$), intent(inout) :: resid(n)
      integer, intent(in) :: ncv
      integer, intent(in) :: ldv
      real(${KIND}$), intent(out) :: v(ldv, ncv)
      integer, intent(inout) :: iparam(11)
      integer, intent(out) :: ipntr(11)
      real(${KIND}$), intent(inout) :: workd(3 * n)
      integer, intent(in) :: lworkl
      real(${KIND}$), intent(inout) :: workl(lworkl)
      integer, intent(inout) :: info
    end subroutine ${PREFIX}$saupd

    ${ARPACK_PROC_PREFIX}$ subroutine ${PREFIX}$seupd(rvec, howmny, sel, d, z, ldz, sigma, bmat,&
        & n, which, nev, tol, resid, ncv, v, ldv, iparam, ipntr, workd, workl, lworkl, info)
    #:if WITH_ARPACK
      import :: ${KIND}$
    #:endif
      logical, intent(in) :: rvec
      character, intent(in) :: howmny
      integer, intent(in) :: ncv
      logical, intent(in) :: sel(ncv)
      integer, intent(in) :: nev
      real(${KIND}$), intent(out) :: d(nev)
      integer, intent(in) :: ldz
      real(${KIND}$), intent(out) :: z(ldz, nev)
      real(${KIND}$), intent(in) :: sigma
      character, intent(in) :: bmat
      integer, intent(in) :: n
      character(2), intent(in) :: which
      real(${KIND}$), intent(in) :: tol
      real(${KIND}$), intent(in) :: resid(n)
      integer, intent(in) :: ldv
      real(${KIND}$), intent(inout) :: v(ldv, ncv)
      integer, intent(in) :: iparam(7)
      integer, intent(inout) :: ipntr(11)
      real(${KIND}$), intent(inout) :: workd(2 * n)
      integer, intent(in) :: lworkl
      real(${KIND}$), intent(inout) :: workl(lworkl)
      integer, intent(inout) :: info
    end subroutine ${PREFIX}$seupd

    ${PARPACK_PROC_PREFIX}$ subroutine p${PREFIX}$saupd(comm, ido, bmat, n, which, nev, tol, resid,&
        & ncv, v, ldv, iparam, ipntr, workd, workl, lworkl, info)
    #:if WITH_PARPACK
      import :: ${KIND}$
    #:endif
      integer, intent(in) :: comm
      integer, intent(inout) :: ido
      character, intent(in) :: bmat
      integer, intent(in) :: n
      character(2), intent(in) :: which
      integer, intent(in) :: nev
      real(${KIND}$), intent(in) :: tol
      real(${KIND}$), intent(inout) :: resid(n)
      integer, intent(in) :: ncv
      integer, intent(in) :: ldv
      real(${KIND}$), intent(out) :: v(ldv, ncv)
      integer, intent(inout) :: iparam(11)
      integer, intent(out) :: ipntr(11)
      real(${KIND}$), intent(inout) :: workd(3 * n)
      integer, intent(in) :: lworkl
      real(${KIND}$), intent(inout) :: workl(lworkl)
      integer, intent(inout) :: info
    end subroutine p${PREFIX}$saupd

    ${PARPACK_PROC_PREFIX}$ subroutine p${PREFIX}$seupd(comm, rvec, howmny, sel, d, z, ldz, sigma,&
        & bmat, n, which, nev, tol, resid, ncv, v, ldv, iparam, ipntr, workd, workl, lworkl, info)
    #:if WITH_PARPACK
      import :: ${KIND}$
    #:endif
      integer, intent(in) :: comm
      logical, intent(in) :: rvec
      character, intent(in) :: howmny
      integer, intent(in) :: ncv
      logical, intent(inout) :: sel(ncv)
      integer, intent(in) :: nev
      real(${KIND}$), intent(out) :: d(nev)
      integer, intent(in) :: ldz
      real(${KIND}$), intent(out) :: z(ldz, nev)
      real(${KIND}$), intent(in) :: sigma
      character, intent(in) :: bmat
      integer, intent(in) :: n
      character(2), intent(in) :: which
      real(${KIND}$), intent(in) :: tol
      real(${KIND}$), intent(in) :: resid(n)
      integer, intent(in) :: ldv
      real(${KIND}$), intent(inout) :: v(ldv, ncv)
      integer, intent(in) :: iparam(7)
      integer, intent(inout) :: ipntr(11)
      real(${KIND}$), intent(inout) :: workd(2 * n)
      integer, intent(in) :: lworkl
      real(${KIND}$), intent(inout) :: workl(lworkl)
      integer, intent(inout) :: info
    end subroutine p${PREFIX}$seupd

  #:endfor
  end interface

end module dftbp_extlibs_arpack


#:if (not WITH_ARPACK) or (not WITH_PARPACK)

!> Defines stubs for ARPACK/PARPACK routines, in case the libraries are not present.
submodule (dftbp_extlibs_arpack) dftbp_extlibs_arpack_stubs
  use dftbp_io_message
  implicit none

contains

#:for PREFIX, _ in PREFIXES_AND_KINDS

#:if not WITH_ARPACK
  module procedure ${PREFIX}$saupd
    call stubError_("${PREFIX}$saupd", "ARPACK")
  end procedure

  module procedure ${PREFIX}$seupd
    call stubError_("${PREFIX}$seupd", "ARPACK")
  end procedure
#:endif

#:if not WITH_PARPACK
  module procedure p${PREFIX}$saupd
    call stubError_("p${PREFIX}$saupd", "PARPACK")
  end procedure

  module procedure p${PREFIX}$seupd
    call stubError_("p${PREFIX}$seupd", "PARPACK")
  end procedure
#:endif

#:endfor


  !! Generates error message, if a stub was called
  subroutine stubError_(routineName, libraryName)
    character(*), intent(in) :: routineName, libraryName

    call error("Internal error: " // routineName // "() called in a build without " // libraryName&
        & // " suppport")

  end subroutine stubError_

end submodule dftbp_extlibs_arpack_stubs

#:endif
