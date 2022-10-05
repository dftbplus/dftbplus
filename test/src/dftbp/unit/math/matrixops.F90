!--------------------------------------------------------------------------------------------------!
!  DFTB+: general package for performing fast atomistic simulations                                !
!  Copyright (C) 2006 - 2025  DFTB+ developers group                                               !
!                                                                                                  !
!  See the LICENSE file for terms of usage and distribution.                                       !
!--------------------------------------------------------------------------------------------------!

#:include "fortuno_serial.fypp"

module test_math_matrixops
  use fortuno_serial, only : all_close, suite => serial_suite_item, test_list
  use dftbp_common_accuracy, only : dp
  use dftbp_common_environment, only : TEnvironment
  use dftbp_common_status, only : TStatus
  use dftbp_math_matrixops, only : adjointLowerTriangle, adjugate, orthonormalizeVectors, pseudoInv
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


  $:TEST("adjugatematrix")

    real(dp), allocatable :: matReal(:,:), CT(:,:)
    type(TStatus) :: errStatus
    integer :: ii, jj, n

    ! very small case
    matReal = reshape([1,2,3,4], [2,2])
    CT = reshape([4,-2,-3,1], [2,2])
    call adjugate(matReal)
    @:ASSERT(all(abs(matReal - CT) < 128_dp * epsilon(0.0_dp)))
    deallocate(matReal)
    deallocate(CT)

    ! slightly larger
    matReal = reshape([-3,2,-5,-1,0,-2,3,-4,1], [3,3])
    CT = reshape([-8,18,-4,-5,12,-1,4,-6,2], [3,3])
    call adjugate(matReal)
    @:ASSERT(all(abs(matReal - CT) < 1024_dp * epsilon(0.0_dp)))
    deallocate(matReal)
    deallocate(CT)

    ! Numerical test between implementations
    n = 10
    allocate(matReal(n,n))
    do ii = 1, n
      do jj = 1, n
        matReal(jj, ii) = sign(ii-jj, ii+jj)
      end do
      matReal(ii, ii) = 0.1_dp * ii
    end do
    CT = matReal
    call adjugate(matReal)
    ! use simple routine (requires matrix to be invertable):
    call adjugate(CT, errStatus)
    @:ASSERT(.not.errStatus%hasError())
    @:ASSERT(all(abs(matReal - CT) < 1024_dp * epsilon(0.0_dp)))

  $:END_TEST()


  $:TEST("pseudoInverse")

    real(dp), allocatable :: A(:,:), Ainv(:,:), work(:,:)
    integer :: ii, m, n
    character(10) :: formatString

    ! small 1D test
    m = 4
    n = 1
    allocate(A(m, n))
    allocate(Ainv(m, n))

    A(:,:) = real(reshape([1,2,3,4], [m,n]), dp)
    work = A
    call pseudoInv(work, Ainv)

    ! check A = A.A.A^~
    @:ASSERT(all(abs(A - matmul(A, matmul(transpose(A),Ainv))) < 256.0_dp * epsilon(1.0_dp)))
    deallocate(A)
    deallocate(Ainv)
    deallocate(work)

    ! small 2D test
    m = 3
    n = 2
    allocate(A(m, n))
    allocate(Ainv(m, n))

    A(:,:) = real(reshape([1,2,3,4,5,6], [m,n]), dp)
    work = A
    call pseudoInv(work, Ainv)

    ! check A = A.A.A^~
    @:ASSERT(all(abs(A - matmul(A, matmul(transpose(A),Ainv))) < 256.0_dp * epsilon(1.0_dp)))
    deallocate(A)
    deallocate(Ainv)
    deallocate(work)

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
