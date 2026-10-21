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
    type(dual_real64) :: x, y, ref
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
    x = [(xEnd + xStart) / 2.0_dp, 1.0_dp]
    write(*,*)'Interpolation at'
    write(*,*)x%get_derivatives()
    y = polyInterUniform(xp, yp, x)
    write(*,*)'Interpolated value and derivatives'
    write(*,*)y%get_derivatives()
    write(*,*)'Exact values'
    call initialize_dual(ref, order)
    ref = [testFn(x%get_derivative(0)), testdFn(x%get_derivative(0)),&
        & testddFn(x%get_derivative(0)), testdddFn(x%get_derivative(0))]
    write(*,*)ref%get_derivatives()
    write(*,*)'Fraction of eps difference'
    write(*,*)abs(ref%get_derivatives() - y%get_derivatives()) / epsilon(0.0_dp)
    @:ASSERT(all(abs(ref%get_derivatives() - y%get_derivatives()) < 10000.0_dp * epsilon(0.0_dp)))
  $:END_TEST()

  $:TEST("interpolate_vector")
    integer, parameter :: nPts = 10, order = 3
    real(dp) :: xStart = 1.0_dp, xEnd = 3.0_dp
    real(dp) :: xp(nPts), yp(2,nPts)
    type(dual_real64) :: x, y(2), ref(2)
    integer :: ii
    do ii = 0, nPts-1
      xp(ii+1) = xStart + (xEnd - xStart) * real(ii, dp) / real(nPts-1, dp)
    end do
    yp(1, :) = testFn(xp)
    yp(2, :) = 2.0_dp * yp(1, :)
    write(*,*)'Data points'
    write(*,*)xp
    write(*,*)yp(1,:)
    write(*,*)yp(2,:)
    write(*,*)
    call initialize_dual(x, order)
    x = [(xEnd + xStart) / 2.0_dp, 1.0_dp]
    write(*,*)'Interpolation at'
    write(*,*)x%get_derivatives()
    y = polyInterUniform(xp, yp, x)
    write(*,*)'Interpolated value and derivatives'
    write(*,*)y(1)%get_derivatives()
    write(*,*)y(2)%get_derivatives()
    write(*,*)'Exact values'
    call initialize_dual(ref, order)
    ref(1) = [testFn(x%get_derivative(0)), testdFn(x%get_derivative(0)),&
        & testddFn(x%get_derivative(0)), testdddFn(x%get_derivative(0))]
    ref(2) = 2.0_dp * ref(1)
    write(*,*)ref(1)%get_derivatives()
    write(*,*)ref(2)%get_derivatives()
    write(*,*)'Fraction of eps difference'
    write(*,*)abs(ref(1)%get_derivatives() - y(1)%get_derivatives()) / epsilon(0.0_dp)
    write(*,*)abs(ref(2)%get_derivatives() - y(2)%get_derivatives()) / epsilon(0.0_dp)
    @:ASSERT(all(abs(ref(1)%get_derivatives() - y(1)%get_derivatives())&
        & < 100000.0_dp * epsilon(0.0_dp)))
    @:ASSERT(all(abs(ref(2)%get_derivatives() - y(2)%get_derivatives())&
        & < 100000.0_dp * epsilon(0.0_dp)))
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
