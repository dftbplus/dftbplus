!--------------------------------------------------------------------------------------------------!
!  DFTB+: general package for performing fast atomistic simulations                                !
!  Copyright (C) 2006 - 2025  DFTB+ developers group                                               !
!                                                                                                  !
!  See the LICENSE file for terms of usage and distribution.                                       !
!--------------------------------------------------------------------------------------------------!

#:include 'error.fypp'

!> Contains routines to locate a value in a sorted array using binary search
module dftbp_math_binarysearch
  use dftbp_common_accuracy, only : dp
  use dftbp_common_status, only : TStatus
  use dftbp_io_message, only : error
  implicit none

  private
  public :: search_int, search_asc_real_geq, search_asc_real_gt, search_des_real_geq,&
      & search_des_real_gt

contains


  !> Integer case for binary search of sorted values to find the jj such that
  !! xVal < xx(jj+1), i.e., the last occurance of the value in an array, or if not present in the
  !! array, the last element smaller than xVal. If xVal < xx(1), jj = 0
  subroutine search_int(jj, xx, xVal)

    !> Located element
    integer, intent(out) :: jj

    !> Array of values in ascending order to search through
    integer, intent(in) :: xx(:)

    !> Value to locate jj for
    integer, intent(in) :: xVal

    integer :: jlower, jupper, jcurr

    jlower = 0
    jupper = size(xx)
    do while (jlower < jupper)
      jcurr = jlower + (jupper - jlower + 1) / 2
      if (xx(jcurr) <= xVal) then
        jlower = jcurr
      else
        jupper = jcurr - 1
      end if
    end do
    jj = jlower

  end subroutine search_int


  !======================================
  ! Ascending ordered real array versions
  !======================================


  !> Real case for binary search of ascending sorted values to find the jj such that
  !! xx(jj) + tol <= xVal, i.e., the last element equal to or smaller than xVal. If xVal < xx(1) -
  !! tol, jj = 0
  subroutine search_asc_real_geq(jj, xx, xVal, tol)

    !> Located element
    integer, intent(out) :: jj

    !> Array of values in ascending order to search through
    real(dp), intent(in) :: xx(:)

    !> Value to locate jj for
    real(dp), intent(in) :: xVal

    !> Tolerance for equality comparison
    real(dp), intent(in), optional :: tol

    integer :: jlower, jupper, jcurr
    real(dp) :: tol_

    tol_ = epsilon(0.0_dp)
    if (present(tol)) tol_ = tol
    jlower = 0
    jupper = size(xx)
    do while (jlower < jupper)
      jcurr = jlower + (jupper - jlower + 1) / 2
      if (xVal - xx(jcurr) >= -tol_) then
        jlower = jcurr
      else
        jupper = jcurr - 1
      end if
    end do
    jj = jlower

  end subroutine search_asc_real_geq


  !> Real case for binary search of ascending sorted values to find the jj such that
  !! xVal < xx(jj+1) - tol, i.e. the last occurance of a value smaller than xVal. If xVal < xx(1) -
  !! tol, jj = 0
  subroutine search_asc_real_gt(jj, xx, xVal, tol)

    !> Located element
    integer, intent(out) :: jj

    !> Array of values in ascending order to search through
    real(dp), intent(in) :: xx(:)

    !> Value to locate jj for
    real(dp), intent(in) :: xVal

    !> Tolerance for equality comparison
    real(dp), intent(in), optional :: tol

    integer :: jlower, jupper, jcurr
    real(dp) :: tol_

    tol_ = epsilon(0.0_dp)
    if (present(tol)) tol_ = tol
    jlower = 0
    jupper = size(xx)
    do while (jlower < jupper)
      jcurr = jlower + (jupper - jlower + 1) / 2
      if (xVal - xx(jcurr) > tol_) then
        jlower = jcurr
      else
        jupper = jcurr - 1
      end if
    end do
    jj = jlower

  end subroutine search_asc_real_gt


  !=======================================
  ! Descending ordered real array versions
  !=======================================


  !> Real case for binary search of decending sorted values to find the first jj such that
  !! xx(jj) + tol >= xVal
  subroutine search_des_real_geq(jj, xx, xVal, tol)

    !> Located element
    integer, intent(out) :: jj

    !> Array of values in ascending order to search through
    real(dp), intent(in) :: xx(:)

    !> Value to locate jj for
    real(dp), intent(in) :: xVal

    !> Tolerance for equality comparison
    real(dp), intent(in), optional :: tol

    integer :: jlower, jupper, jcurr
    real(dp) :: tol_

    tol_ = epsilon(0.0_dp)
    if (present(tol)) tol_ = tol
    jlower = 0
    jupper = size(xx)
    do while (jlower < jupper)
      jcurr = jlower + (jupper - jlower + 1) / 2
      if ((xx(jcurr) - xVal) >= -tol_) then
        jlower = jcurr
      else
        jupper = jcurr - 1
      end if
    end do
    jj = jlower

  end subroutine search_des_real_geq


  !> Real case for binary search of descending sorted values to find first element jj such that
  !! xx(jj) - tol > xVal
  subroutine search_des_real_gt(jj, xx, xVal, tol)

    !> Located element
    integer, intent(out) :: jj

    !> Array of values in ascending order to search through
    real(dp), intent(in) :: xx(:)

    !> Value to locate jj for
    real(dp), intent(in) :: xVal

    !> Tolerance for equality comparison
    real(dp), intent(in), optional :: tol

    integer :: jlower, jupper, jcurr
    real(dp) :: tol_

    tol_ = epsilon(0.0_dp)
    if (present(tol)) tol_ = tol
    jlower = 0
    jupper = size(xx)
    do while (jlower < jupper)
      jcurr = jlower + (jupper - jlower + 1) / 2
      if ((xx(jcurr) - xVal) > tol_) then
        jlower = jcurr
      else
        jupper = jcurr - 1
      end if
    end do
    jj = jlower

  end subroutine search_des_real_gt

end module dftbp_math_binarysearch
