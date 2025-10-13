!--------------------------------------------------------------------------------------------------!
!  DFTB+: general package for performing fast atomistic simulations                                !
!  Copyright (C) 2006 - 2025  DFTB+ developers group                                               !
!                                                                                                  !
!  See the LICENSE file for terms of usage and distribution.                                       !
!--------------------------------------------------------------------------------------------------!
!> This module defines TOrbital, which, aside from holding the angular momentum,
!! may be subclassed by one of three concrete radial function implementations:
!! Slater Type Orbitals (STO), Gaussian Type Orbitals (GTO) or a radial lookup table (LUT),
!! allowing arbitrary radial functions to be used. Any orbital may be resampled onto a LUT
!! with a new chosen resolution and cutoff.

#:include 'common.fypp'
#:set pure = "" if defined('WITH_ASSERT') else "pure"

module dftbp_wavegrid_basis_orbital
  use dftbp_common_accuracy, only : dp
  implicit none


  public :: TOrbital
  
  type, abstract :: TOrbital
    !> Angular momentum (l)
    integer :: angMom = -1
    !> Square of the Cutoff, after which the orbital is assumed to be zero
    real(dp) :: cutoffSq
  contains
    procedure(IGetRadial), deferred :: getRadial
    procedure(IAssign), deferred, pass(lhs) :: assign
    generic, public :: assignment(=) => assign
  end type TOrbital

  abstract interface
    ${pure}$ function IGetRadial(this, r) result(val)
      import :: TOrbital, dp
      class(TOrbital), intent(in) :: this
      real(dp), intent(in) :: r
      real(dp) :: val
    end function IGetRadial

    subroutine IAssign(lhs, rhs)
      import :: TOrbital
      class(TOrbital), intent(out) :: lhs
      class(TOrbital), intent(in) :: rhs
    end subroutine IAssign
  end interface

end module dftbp_wavegrid_basis_orbital

