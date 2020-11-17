!--------------------------------------------------------------------------------------------------!
!  DFTB+: general package for performing fast atomistic simulations                                !
!  Copyright (C) 2006 - 2020  DFTB+ developers group                                               !
!                                                                                                  !
!  See the LICENSE file for terms of usage and distribution.                                       !
!--------------------------------------------------------------------------------------------------!

#:include 'common.fypp'

!> Wrappers for the functions erf(x) and erfc(x).
!>
!> Based on the preprocessor settings, the error function is wrapped differently:
!>
!> a) no special definitions: the intrinsic error function is used (officially first available in
!>    the Fortran 2008 standard, but most F95/2003 compilers already implements this).
!>
!> b) INTERNAL_ERFC is defined: erf(x) and erfc(x) are internally calculated by the code.
!>
module dftbp_errorfunction
  use dftbp_accuracy
#:if INTERNAL_ERFC
  use dftbp_erfcalc, only : erf, erfc
#:endif
  implicit none
  private

  public :: erfwrap, erfcwrap

contains


  !> Calculates the value of the error function.
  elemental function erfwrap(xx) result(res)

    !> Function argument.
    real(dp), intent(in) :: xx

    !> erf(x)
    real(dp) :: res

    res = erf(xx)

  end function erfwrap


  !> Calculates the value of the complementary error function.
  elemental function erfcwrap(xx) result(res)

    !> Function argument.
    real(dp), intent(in) :: xx

    !> erf(x)
    real(dp) :: res

    res = erfc(xx)

  end function erfcwrap

end module dftbp_errorfunction
