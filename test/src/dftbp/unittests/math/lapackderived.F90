!--------------------------------------------------------------------------------------------------!
!  DFTB+: general package for performing fast atomistic simulations                                !
!  Copyright (C) 2006 - 2023  DFTB+ developers group                                               !
!                                                                                                  !
!  See the LICENSE file for terms of usage and distribution.                                       !
!--------------------------------------------------------------------------------------------------!

#:include "fytest.fypp"

#:block TEST_SUITE("lapackderived")
  use dftbp_common_accuracy, only : dp
  use dftbp_math_lapackderived, only : pseudoInv
  use dftbp_common_status, only : TStatus
  implicit none

#:contains

  #:block TEST_FIXTURE("pseudoInverse")

    real(dp), allocatable :: A(:,:), Ainv(:,:), work(:,:)
    integer :: ii, m, n
    character(10) :: formatString

  #:contains

    #:block TEST("1D")
      m = 4
      n = 1
      allocate(A(m, n))
      allocate(Ainv(m, n))

      A(:,:) = real(reshape([1,2,3,4], [m,n]), dp)
      work = A
      call pseudoInv(work, Ainv)

      ! check A = A.A.A^~
      @:ASSERT(all(abs(A - matmul(A, matmul(transpose(A),Ainv))) < epsilon(1.0_dp) * 2**8))
      deallocate(A)
      deallocate(Ainv)
      deallocate(work)
    #:endblock

    #:block TEST("2D")
      m = 3
      n = 2
      allocate(A(m, n))
      allocate(Ainv(m, n))

      A(:,:) = real(reshape([1,2,3,4,5,6], [m,n]), dp)
      work = A
      call pseudoInv(work, Ainv)

      ! check A = A.A.A^~
      @:ASSERT(all(abs(A - matmul(A, matmul(transpose(A),Ainv))) < epsilon(1.0_dp) * 2**8))
      deallocate(A)
      deallocate(Ainv)
      deallocate(work)
    #:endblock

  #:endblock TEST_FIXTURE

#:endblock TEST_SUITE


@:TEST_DRIVER()
