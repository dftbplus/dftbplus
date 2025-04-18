!--------------------------------------------------------------------------------------------------!
!  DFTB+: general package for performing fast atomistic simulations                                !
!  Copyright (C) 2006 - 2025  DFTB+ developers group                                               !
!                                                                                                  !
!  See the LICENSE file for terms of usage and distribution.                                       !
!--------------------------------------------------------------------------------------------------!

#:include "fortuno_serial.fypp"

module test_math_interpolation
  use fortuno_serial, only : all_close, suite => serial_suite_item, test_list
  use dftbp_common_accuracy, only : dp
  use dftbp_extlibs_dnaoad
  use dftbp_math_interpolation, only : polyInterUniform
  $:FORTUNO_SERIAL_IMPORTS()
  implicit none

  private
  public :: tests

contains

  $:TEST("interpolate_scalar")
    integer, parameter :: nPts = 10, order = 3
    real(dp) :: xStart = 1.0_dp, xEnd = 3.0_dp
    real(dp) :: xp(nPts), yp(nPts)
    type(dual_real64) :: x, y
    integer :: ii
    do ii = 0, nPts-1
      xp(ii+1) = xStart + (xEnd - xStart) * real(ii, dp) / real(nPts-1, dp)
    end do
    yp(:) = testFn(xp)
    write(*,*)'Data points'
    write(*,*)xp
    write(*,*)yp
    write(*,*)
    call initialize_dual(x, order)
    x = (xEnd + xStart) / 2.0_dp
    x%f(1) = 1.0_dp
    write(*,*)'Interpolation at'
    write(*,*)x%f
    y = polyInterUniform(xp, yp, x)
    write(*,*)'Interpolated value and derivatives'
    write(*,*)y%f
    write(*,*)'Exact values'
    x%f(:) = [testFn(x%f(0)), testdFn(x%f(0)), testddFn(x%f(0)), testdddFn(x%f(0))]
    write(*,*)x%f
    @:ASSERT(all(abs(x%f - y%f) < 10000.0_dp * epsilon(0.0_dp)))
  $:END_TEST()

  elemental function testFn(x) result(y)

    real(dp), intent(in) :: x
    real(dp) :: y

    y = 1.0_dp - 2.0_dp * x + 3.0_dp * x * x - 4.0_dp * x * x * x

  end function testFn

  elemental function testdFn(x) result(y)

    real(dp), intent(in) :: x
    real(dp) :: y

    y = -2.0_dp + 6.0_dp * x - 12.0_dp * x * x

  end function testdFn

  elemental function testddFn(x) result(y)

    real(dp), intent(in) :: x
    real(dp) :: y

    y = 6.0_dp - 24.0_dp * x

  end function testddFn

  elemental function testdddFn(x) result(y)

    real(dp), intent(in) :: x
    real(dp) :: y

    y = -24.0_dp

  end function testdddFn

  function tests()
    type(test_list) :: tests

    tests = test_list([&
        suite("interpolation", test_list([&
            $:TEST_ITEMS()
        ]))&
    ])
    $:STOP_ON_MISSING_TEST_ITEMS()

  end function tests

end module test_math_interpolation
