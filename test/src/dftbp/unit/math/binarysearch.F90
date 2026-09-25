!--------------------------------------------------------------------------------------------------!
!  DFTB+: general package for performing fast atomistic simulations                                !
!  Copyright (C) 2006 - 2025  DFTB+ developers group                                               !
!                                                                                                  !
!  See the LICENSE file for terms of usage and distribution.                                       !
!--------------------------------------------------------------------------------------------------!

#:include "fortuno_serial.fypp"
module test_math_binarysearch
  use dftbp_common_accuracy, only : dp, rsp
  use dftbp_math_binarysearch, only : search_asc_real_geq, search_asc_real_gt
  use dftbp_math_binarysearch, only : search_des_real_geq, search_des_real_gt
  use dftbp_math_binarysearch, only : search_int
  use fortuno_serial, only : suite => serial_suite_item, test_list, all_close
  $:FORTUNO_SERIAL_IMPORTS()
  implicit none

  private
  public :: tests

contains

  $:TEST("binarysearch_int")
    integer, parameter :: data(8) = [1, 3, 3, 4, 5, 5, 6 ,7]
    integer, parameter :: x(*) = [0, 8, 5, 3, 7, 2, 1]
    integer, parameter :: nn = size(x)
    integer, parameter :: expect(nn) = [0, 8, 6, 3, 8, 1, 1]
    integer :: ii, res(nn)

    write(*,*)
    write(*,*)'Integer'
    do ii = 1, nn
      call search_int(res(ii), data, x(ii))
    end do
    write(*,*)
    write(*,'(A,*(I2))')'#',(ii, ii=1, size(data))
    write(*,"(1X,*(I2))")data
    write(*,"(1X,A,T10,*(I2))")'Test', x
    write(*,"(1X,A,T10,*(I2))")'Expect', expect
    write(*,"(1X,A,T10,*(I2))")'Got', res
    @:ASSERT(all(res == expect))

  $:END_TEST()


  $:TEST("binarysearch_real_simple")
    real(dp), parameter :: tol = real(epsilon(0.0_rsp),dp)
    real(dp), parameter :: data(*) = [&
        & 0.1_dp, 0.2_dp, 0.3_dp, 0.4_dp, 0.4_dp, 0.5_dp, 0.6_dp, 0.8_dp ]
    real(dp), parameter :: x(*) = [-1.0_dp, 0.1_dp + tol/2.0_dp, 10.0_dp, 0.8_dp -tol/2.0_dp,&
        & 0.4_dp]
    integer, parameter :: nn = size(x), nData = size(data)
    integer, parameter :: expect_geq(nn) = [0,1,8,8,5]
    integer, parameter :: expect_gt(nn) = [0,0,8,7,3]
    integer, parameter :: expect_ageq(nn) = [8,8,0,1,5]
    integer, parameter :: expect_agt(nn) = [8,7,0,0,3]
    integer :: ii, res(size(expect_geq))

    write(*,*)
    write(*,*)'Real ascending'
    do ii = 1, nn
      call search_asc_real_geq(res(ii), data, x(ii), tol)
    end do
    write(*,*)'max jj : x(jj) <= x - tol'
    write(*,'(A,I2,*(I4))')'#',(ii, ii=1, nData)
    write(*,"(1X,*(F4.1))")data
    write(*,"(1X,A,T10,*(F6.1))")'Test', x
    write(*,"(1X,A,T10,*(E9.2))")'+/-',tol
    write(*,"(1X,A,T10,*(I2))")'Expect', expect_geq
    write(*,"(1X,A,T10,*(I2))")'Got', res
    @:ASSERT(all(res == expect_geq))
    do ii = 1, nn
      call search_asc_real_gt(res(ii), data, x(ii), tol)
    end do
    write(*,*)'max jj : x(jj) < x - tol'
    write(*,'(A,I2,*(I4))')'#',(ii, ii=1, nData)
    write(*,"(1X,*(F4.1))")data
    write(*,"(1X,A,T10,*(F6.1))")'Test', x
    write(*,"(1X,A,T10,*(E9.2))")'+/-',tol
    write(*,"(1X,A,T10,*(I2))")'Expect', expect_gt
    write(*,"(1X,A,T10,*(I2))")'Got', res
    @:ASSERT(all(res == expect_gt))
    write(*,*)'Real descending'
    do ii = 1, nn
      call search_des_real_geq(res(ii), data(nData:1:-1), x(ii), tol)
    end do
    write(*,*)'max jj : x(jj) >= x + tol'
    write(*,'(A,I2,*(I4))')'#',(ii, ii=1, nData)
    write(*,"(1X,*(F4.1))")data(nData:1:-1)
    write(*,"(1X,A,T10,*(F6.1))")'Test', x
    write(*,"(1X,A,T10,*(E9.2))")'+/-',tol
    write(*,"(1X,A,T10,*(I2))")'Expect', expect_ageq
    write(*,"(1X,A,T10,*(I2))")'Got', res
    @:ASSERT(all(res == expect_ageq))
    do ii = 1, nn
      call search_des_real_gt(res(ii), data(nData:1:-1), x(ii), tol)
    end do
    write(*,*)'max jj : x(jj) > x + tol'
    write(*,'(A,I2,*(I4))')'#',(ii, ii=1, nData)
    write(*,"(1X,*(F4.1))")data(nData:1:-1)
    write(*,"(1X,A,T10,*(F6.1))")'Test', x
    write(*,"(1X,A,T10,*(E9.2))")'+/-',tol
    write(*,"(1X,A,T10,*(I2))")'Expect', expect_agt
    write(*,"(1X,A,T10,*(I2))")'Got', res
    @:ASSERT(all(res == expect_agt))

  $:END_TEST()

  function tests()
    type(test_list) :: tests

    tests = test_list([&
        suite("binarysearch", test_list([&
            $:TEST_ITEMS()
        ]))&
    ])
    $:STOP_ON_MISSING_TEST_ITEMS()

  end function tests

end module test_math_binarysearch
