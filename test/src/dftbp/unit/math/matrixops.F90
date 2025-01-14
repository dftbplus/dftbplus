!--------------------------------------------------------------------------------------------------!
!  DFTB+: general package for performing fast atomistic simulations                                !
!  Copyright (C) 2006 - 2025  DFTB+ developers group                                               !
!                                                                                                  !
!  See the LICENSE file for terms of usage and distribution.                                       !
!--------------------------------------------------------------------------------------------------!

#:include "fortuno_serial.fypp"

module test_math_matrixops
  use dftbp_common_accuracy, only : dp
  use dftbp_common_environment, only : TEnvironment
  use dftbp_math_matrixops, only : adjointLowerTriangle, orthonormalizeVectors
  use fortuno_serial, only : suite => serial_suite_item, test_list, all_close
  $:FORTUNO_SERIAL_IMPORTS()
  implicit none

  private
  public :: tests

contains


  $:TEST("adjointLowerTriangle_real")
    integer, parameter :: nn = 3
    real(dp) :: matReal(nn, nn)
    integer :: ii, jj, kk

    matReal(:,:) = 0.0_dp
    kk = 0
    do ii = 1, size(matReal, dim=1)
      do jj = ii, size(matReal, dim=2)
        kk = kk + 1
        matReal(jj, ii) = real(kk, dp)
      end do
    end do
    call adjointLowerTriangle(matReal)
    @:ASSERT(all(matReal > epsilon(0.0_dp)))
    @:ASSERT(all(abs(matReal - transpose(matReal)) < epsilon(0.0_dp)))
  $:END_TEST()


  $:TEST("adjointLowerTriangle_complex")
    integer, parameter :: nn = 3
    complex(dp) :: matCplx(nn, nn)
    integer :: ii, jj, kk

    matCplx(:,:) = cmplx(0, 0, dp)
    kk = 0
    do ii = 1, size(matCplx, dim=1)
      do jj = ii, size(matCplx, dim=2)
        kk = kk + 1
        matCplx(jj, ii) = real(kk, dp)
      end do
    end do
    kk = 0
    do ii = 1, size(matCplx, dim=1)
      do jj = ii + 1, size(matCplx, dim=2)
        kk = kk + 1
        matCplx(jj, ii) = matCplx(jj, ii) + cmplx(0, kk, dp)
      end do
    end do
    call adjointLowerTriangle(matCplx)
    @:ASSERT(all(abs(matCplx) > epsilon(0.0_dp)))
    @:ASSERT(all(abs(matCplx - transpose(conjg(matCplx))) < epsilon(0.0_dp)))
  $:END_TEST()


  $:TEST("orthonormalizeVectors")
    integer, parameter :: n = 50
    real(dp), parameter :: atol = 1000.0_dp * epsilon(0.0_dp)
    type(TEnvironment) :: env
    integer :: ii, jj, kk
    real(dp) :: matReal(n,n), eye(n,n)

    matReal(:,:) = 0.0_dp
    eye(:,:) = 0.0_dp
    kk = 0
    do ii = 1, n
      eye(ii,ii) = 1.0_dp
      do jj = ii, n
        kk = kk + 1
          matReal(jj, ii) = real(kk,dp)
      end do
    end do
    call orthonormalizeVectors(env, 1, n, matReal)
    matReal(:,:) = matmul(transpose(matReal), matReal)
    @:ASSERT(all_close(matReal, eye, atol=atol))

  $:END_TEST()


  function tests()
    type(test_list) :: tests

    tests = test_list([&
        suite("matrixops", test_list([&
            $:TEST_ITEMS()
        ]))&
    ])
    $:STOP_ON_MISSING_TEST_ITEMS()

  end function tests

end module test_math_matrixops
