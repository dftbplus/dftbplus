!--------------------------------------------------------------------------------------------------!
!  DFTB+: general package for performing fast atomistic simulations                                !
!  Copyright (C) 2006 - 2023  DFTB+ developers group                                               !
!                                                                                                  !
!  See the LICENSE file for terms of usage and distribution.                                       !
!--------------------------------------------------------------------------------------------------!

#:include 'common.fypp'

!> Contains functions for some common operations built from lapack calls
module dftbp_math_lapackderived
  use dftbp_common_accuracy, only : dp
  use dftbp_common_status, only : TStatus
  use dftbp_math_lapackroutines, only : gesvd
  implicit none

  private
  public :: pseudoInv


contains

  !> Moore-Penrose pseudo-inverse of general rectangular matrix
  subroutine pseudoInv(A, Ainv)

    !> Matrix A(m,n), overwritten on output
    real(dp), intent(inout) :: A(:,:)

    !> Pseudoinverse, returned as transpose Ainv(m,n)
    real(dp), intent(out) :: Ainv(:,:)

    real(dp), allocatable :: U(:,:), sigma(:), Vt(:,:)
    integer :: m, n, mn, ii

    m = size(A, dim=1)
    n = size(A, dim=2)
    mn = min(m, n)

    @:ASSERT(all(shape(Ainv) == [m,n]))

    allocate(U(m,mn))
    allocate(Vt(mn,n))
    allocate(sigma(mn))

    call gesvd(A, U, sigma, Vt)

    where(sigma > epsilon(0.0_dp))
      sigma = 1.0_dp / sigma
    elsewhere
      sigma = 0.0_dp
    end where

    do ii = 1, mn
      U(:, ii) = U(:, ii) * sigma(ii)
    end do

    Ainv(:,:) = matmul(U,Vt)

  end subroutine pseudoInv

end module dftbp_math_lapackderived
