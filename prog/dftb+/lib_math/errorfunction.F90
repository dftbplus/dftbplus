!--------------------------------------------------------------------------------------------------!
!  DFTB+: general package for performing fast atomistic simulations                                !
!  Copyright (C) 2017  DFTB+ developers group                                                      !
!                                                                                                  !
!  See the LICENSE file for terms of usage and distribution.                                       !
!--------------------------------------------------------------------------------------------------!

#:include 'common.fypp'

!> Wrappers for the functions erf(x) and erfc(x).
!>
!> Based on the preprocessor settings, the error function is wrapped differently:
!> <ol>
!>   <li> no special definitions: the intrinsic error function is used (officially
!>   first available in the Fortran 2008 standard, but most F95/2003 compilers
!>   already implements this). </li>
!>   <li> EXTERNALERFC is defined: single precision and double precision external
!>   routines are expected (erf(x), erfc(x), derf(x), derfc(x)).</li>
!>   <li>INTERNALERFC is defined: erf(x) and erfc(x) are internally calculated
!>   by the code.</li>
!> </ol>
module errorfunction
  use accuracy
#:if INTERNAL_ERFC
  use erfcalc, only: erf, erfc
#:endif
  implicit none
  private

  public :: erfwrap, erfcwrap


contains

  !> Calculates the value of the error function.
  function erfwrap(xx) result(res)
    !> Function argument.
    real(dp), intent(in) :: xx
    !> erf(x)
    real(dp) :: res

    res = erf(xx)

  end function erfwrap


  !> Calculates the value of the complementary error function.
  function erfcwrap(xx) result(res)
    !> Function argument.
    real(dp), intent(in) :: xx
    !> erf(x)
    real(dp) :: res

    res = erfc(xx)

  end function erfcwrap

end module errorfunction
