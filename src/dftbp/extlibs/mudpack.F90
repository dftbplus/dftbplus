!--------------------------------------------------------------------------------------------------!
!  DFTB+: general package for performing fast atomistic simulations                                !
!  Copyright (C) 2006 - 2025  DFTB+ developers group                                               !
!                                                                                                  !
!  See the LICENSE file for terms of usage and distribution.                                       !
!--------------------------------------------------------------------------------------------------!

!> Interfaces for the (double precision) MUDPACK multigrid solver routines.
module dftbp_extlibs_mudpack
  use dftbp_common_accuracy, only : rdp
  implicit none (type, external)

  private
  public :: mud3, mud3sp

  ! Interfaces of the callback routines used by mudpack.
  ! The dummy arguments of the call-back routines carry no intent in order to match the
  ! characteristics of the actual call-back routines passed to the solver.
  abstract interface

    !> Coefficients of a single direction of a separable elliptic PDE (for mud3sp)
    subroutine mudpackCoef1d(x, cxx, cx, cex)
      import :: rdp
      implicit none
      real(rdp) :: x, cxx, cx, cex
    end subroutine mudpackCoef1d

    !> Coefficients of the general separable elliptic PDE in 3D (for mud3)
    subroutine mudpackCoef3d(x, y, z, cxx, cyy, czz, cx, cy, cz, ce)
      import :: rdp
      implicit none
      real(rdp) :: x, y, z, cxx, cyy, czz, cx, cy, cz, ce
    end subroutine mudpackCoef3d

    !> Mixed boundary conditions
    subroutine mudpackBndyc(kbdy, xory, yorz, alfa, gbdy)
      import :: rdp
      implicit none
      integer :: kbdy
      real(rdp) :: xory, yorz, alfa, gbdy
    end subroutine mudpackBndyc

  end interface

  ! Interfaces of the mudpack routines
  interface

    !> Multigrid solver for a general (non-separable) 3D elliptic PDE
    subroutine mud3(iparm, fparm, work, coef, bndyc, rhs, phi, mgopt, ierror)
      import :: rdp, mudpackCoef3d, mudpackBndyc
      implicit none
      integer :: iparm(23)
      real(rdp) :: fparm(8)
      real(rdp) :: work(*)
      procedure(mudpackCoef3d) :: coef
      procedure(mudpackBndyc) :: bndyc
      real(rdp) :: rhs(*)
      real(rdp) :: phi(*)
      integer :: mgopt(4)
      integer :: ierror
    end subroutine mud3

    !> Multigrid solver for a separable 3D elliptic PDE
    subroutine mud3sp(iparm, fparm, work, cfx, cfy, cfz, bndyc, rhs, phi, mgopt, ierror)
      import :: rdp, mudpackCoef1d, mudpackBndyc
      implicit none
      integer :: iparm(22)
      real(rdp) :: fparm(8)
      real(rdp) :: work(*)
      procedure(mudpackCoef1d) :: cfx
      procedure(mudpackCoef1d) :: cfy
      procedure(mudpackCoef1d) :: cfz
      procedure(mudpackBndyc) :: bndyc
      real(rdp) :: rhs(*)
      real(rdp) :: phi(*)
      integer :: mgopt(4)
      integer :: ierror
    end subroutine mud3sp

  end interface

end module dftbp_extlibs_mudpack
