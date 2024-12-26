!--------------------------------------------------------------------------------------------------!
!  DFTB+: general package for performing fast atomistic simulations                                !
!  Copyright (C) 2006 - 2023  DFTB+ developers group                                               !
!                                                                                                  !
!  See the LICENSE file for terms of usage and distribution.                                       !
!--------------------------------------------------------------------------------------------------!

#:include "fytest.fypp"

#:block TEST_SUITE("matrixops")
  use dftbp_common_accuracy, only : dp
  use dftbp_common_environment, only : TEnvironment
  use dftbp_math_matrixops, only : adjointLowerTriangle, orthonormalizeVectors
  implicit none

#:contains

  #:block TEST_FIXTURE("denseSymmetry")

     integer :: ii, jj, kk
     integer, parameter :: n = 3
     real(dp) :: matReal(n,n)
     complex(dp) :: matCmplx(n,n)

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
      call adjointLowerTriangle(matReal)
      @:ASSERT(all(matReal > epsilon(0.0_dp)))
      @:ASSERT(all(abs(matReal - transpose(matReal)) < epsilon(0.0_dp)))
    #:endblock

    #:block TEST("complex")
      @:ASSERT(size(matCmplx, dim=1) == size(matCmplx, dim=2))
      matCmplx(:,:) = cmplx(0,0,dp)
      kk = 0
      do ii = 1, size(matCmplx, dim=1)
        do jj = ii, size(matCmplx, dim=2)
          kk = kk + 1
          matCmplx(jj, ii) = real(kk,dp)
        end do
      end do
      kk = 0
      do ii = 1, size(matCmplx, dim=1)
        do jj = ii+1, size(matCmplx, dim=2)
          kk = kk + 1
          matCmplx(jj, ii) = matCmplx(jj, ii) + cmplx(0,kk,dp)
        end do
      end do
      call adjointLowerTriangle(matCmplx)
      @:ASSERT(all(abs(matCmplx) > epsilon(0.0_dp)))
      @:ASSERT(all(abs(matCmplx - transpose(conjg(matCmplx))) < epsilon(0.0_dp)))
    #:endblock

  #:endblock TEST_FIXTURE


  #:block TEST_FIXTURE("orthogonalise")

     type(TEnvironment) :: env
     integer :: ii, jj, kk
     integer, parameter :: n = 50
     real(dp) :: matReal(n,n), overReal(n,n), eye(n,n)
     complex(dp) :: matCmplx(n,n), overCmplx(n,n)

  #:contains


    #:block TEST("real")
      matReal(:,:) = 0.0_dp
      kk = 0
      eye(:,:) = 0.0_dp
      do ii = 1, n
        eye(ii,ii) = 1.0_dp
        do jj = ii, n
          kk = kk + 1
          matReal(jj, ii) = real(kk,dp)
        end do
      end do
      call orthonormalizeVectors(env, 1, n, matReal)
      matReal(:,:) = matmul(transpose(matReal),matReal)
      @:ASSERT(all(abs(matReal - eye) < 1000.0_dp * epsilon(0.0_dp)))
    #:endblock

  #:endblock TEST_FIXTURE

#:endblock TEST_SUITE


@:TEST_DRIVER()
