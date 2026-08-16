!--------------------------------------------------------------------------------------------------!
!  DFTB+: general package for performing fast atomistic simulations                                !
!  Copyright (C) 2006 - 2025  DFTB+ developers group                                               !
!                                                                                                  !
!  See the LICENSE file for terms of usage and distribution.                                       !
!--------------------------------------------------------------------------------------------------!

!> Contains routines to sum values accurately
module dftbp_math_summation
  use dftbp_common_accuracy, only : dp
  implicit none

  private
  public :: kahan, kahanSum

contains

  !> Kahan–Babušk compensated summation with Neumaier improvement
  !! Note, compiler should not be allowed to optimise this out
  subroutine kahan(summation, input, compensator)

    !> Value being accumulated over sum
    real(dp), intent(inout) :: summation

    !> Contribution to add
    real(dp), intent(in) :: input

    !> Compensator for underflow
    real(dp), intent(inout) :: compensator

    real(dp), volatile :: tmp

    tmp = summation + input
    ! Neumaier improvement
    if (abs(summation) >= abs(input)) then
      compensator = compensator + (summation - tmp) + input
    else
      compensator = compensator + (input - tmp) + summation
    end if
    summation = tmp

  end subroutine kahan


  !> Sum function for a 1D array
  function kahanSum(input)

    !> Array to total up
    real(dp), intent(in) :: input(:)

    !> Function result
    real(dp) :: kahanSum

    real(dp) :: compensator
    integer :: ii

    compensator = 0.0_dp
    kahanSum = 0.0_dp
    do ii = 1, size(input)
      call kahan(kahanSum, input(ii), compensator)
    end do
    kahanSum = kahanSum + compensator

  end function kahanSum

end module dftbp_math_summation
