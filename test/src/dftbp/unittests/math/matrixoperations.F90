!--------------------------------------------------------------------------------------------------!
!  DFTB+: general package for performing fast atomistic simulations                                !
!  Copyright (C) 2006 - 2023  DFTB+ developers group                                               !
!                                                                                                  !
!  See the LICENSE file for terms of usage and distribution.                                       !
!--------------------------------------------------------------------------------------------------!

#:include "fytest.fypp"

#:block TEST_SUITE("matrixoperations")
  use dftbp_common_accuracy, only : dp
  use dftbp_math_matrixoperations, only : triangleCopySquareMatrix
  implicit none

#:contains

  #:block TEST_FIXTURE("denseSymmetry")

     integer :: ii, jj, kk
     integer, parameter :: n = 3
     real(dp) :: matReal(n,n)
     complex(dp) :: matCplx(n,n)

  #:contains

    #:block TEST("real")
       @:ASSERT(size(matReal, dim=1) == size(matReal, dim=2))
       matReal(:,:) = 0.0_dp
       kk = 0
       do ii = 1, size(matReal, dim=1)
         do jj = ii, size(matReal, dim=2)
           kk = kk + 1
           matReal(jj, ii) = real(kk,dp)
         end do
       end do
       call triangleCopySquareMatrix(matReal)
       @:ASSERT(all(matReal > epsilon(0.0_dp)))
       @:ASSERT(all(abs(matReal - transpose(matReal)) < epsilon(0.0_dp)))
    #:endblock

    #:block TEST("complex")
       @:ASSERT(size(matCplx, dim=1) == size(matCplx, dim=2))
       matCplx(:,:) = cmplx(0,0,dp)
       kk = 0
       do ii = 1, size(matCplx, dim=1)
         do jj = ii, size(matCplx, dim=2)
           kk = kk + 1
           matCplx(jj, ii) = real(kk,dp)
         end do
       end do
       kk = 0
       do ii = 1, size(matCplx, dim=1)
         do jj = ii+1, size(matCplx, dim=2)
           kk = kk + 1
           matCplx(jj, ii) = matCplx(jj, ii) + cmplx(0,kk,dp)
         end do
       end do
       call triangleCopySquareMatrix(matCplx)
       @:ASSERT(all(abs(matCplx) > epsilon(0.0_dp)))
       @:ASSERT(all(abs(matCplx - transpose(conjg(matCplx))) < epsilon(0.0_dp)))
    #:endblock

  #:endblock TEST_FIXTURE

#:endblock TEST_SUITE


@:TEST_DRIVER()
