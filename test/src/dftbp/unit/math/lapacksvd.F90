!--------------------------------------------------------------------------------------------------!
!  DFTB+: general package for performing fast atomistic simulations                                !
!  Copyright (C) 2006 - 2025  DFTB+ developers group                                               !
!                                                                                                  !
!  See the LICENSE file for terms of usage and distribution.                                       !
!--------------------------------------------------------------------------------------------------!

#:include "fortuno_serial.fypp"

module test_math_lapackroutines
  use dftbp_common_accuracy, only : dp
  use dftbp_math_lapackroutines, only : gesvd
  use fortuno_serial, only : suite => serial_suite_item, test_list, all_close
  $:FORTUNO_SERIAL_IMPORTS()
  implicit none

  private
  public :: tests

contains

  $:TEST("realAA")
    real(dp) :: A(4, 5), U(4,4), Vt(5,5), sigma(4)
    integer :: ii
    call testFullMatrix(A)
    call gesvd(A, u, sigma, vt, jobu="A", jobvt="A")
    call testFullMatrix(A)
    do ii = 1, 4
      U(:,ii) =  U(:,ii) * sigma(ii)
    end do
    A(:,:) = A -matmul(U(:,:4), Vt(:4,:))
    write(*,*)
    write(*,"(4E12.4)")A
    write(*,*)sigma(:4) - singularValues()
    @:ASSERT(all(abs(sigma(:4) - singularValues()) < epsilon(0.0_dp)))
    @:ASSERT(all(abs(A) < 10.0_dp * epsilon(0.0_dp)))
  $:END_TEST()

  $:TEST("realSS")
    real(dp) :: A(10, 10), U(6,6), Vt(6,7), sigma(6)
    integer, parameter :: m = 4
    integer, parameter :: n = 5
    integer :: ii, mn
    mn = min(m,n)
    call testFullMatrix(A)
    call gesvd(A, u, sigma, vt, jobu="S", jobvt="S", m=m, n=n)
    call testFullMatrix(A)
    do ii = 1, mn
      U(:m,ii) =  U(:m,ii) * sigma(ii)
    end do
    A(:m,:n) = A(:m,:n) - matmul(U(:m,:mn), Vt(:mn,:n))
    write(*,*)
    write(*,"(4E12.4)")A(:m,:n)
    write(*,*)sigma(:4) - singularValues()
    @:ASSERT(all(abs(sigma(:mn) - singularValues()) < epsilon(0.0_dp)))
    @:ASSERT(all(abs(A) < 10.0_dp * epsilon(0.0_dp)))
  $:END_TEST()

  !> Test matrix for decomposition. Has
  subroutine testFullMatrix(A)

    real(dp), intent(out) :: A(:,:)

    ! Matrix from https://en.wikipedia.org/wiki/Singular_value_decomposition 2025/11/24
    A(:,:) = 0.0_dp
    A(1,1) = 1.0_dp
    A(1,5) = 2.0_dp
    A(2,3) = 3.0_dp
    A(4,2) = 2.0_dp

  end subroutine testFullMatrix

  !> Test matrix singular values
  function singularValues()

    real(dp) :: singularValues(4)

    singularValues(:) = [3.0_dp, sqrt(5.0_dp), 2.0_dp, 0.0_dp]

  end function singularValues


  function tests()
    type(test_list) :: tests

    tests = test_list([&
        suite("lapack_svd", test_list([&
            $:TEST_ITEMS()
        ]))&
    ])
    $:STOP_ON_MISSING_TEST_ITEMS()

  end function tests

end module test_math_lapackroutines
