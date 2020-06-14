!--------------------------------------------------------------------------------------------------!
!  DFTB+: general package for performing fast atomistic simulations                                !
!  Copyright (C) 2006 - 2020  DFTB+ developers group                                               !
!                                                                                                  !
!  See the LICENSE file for terms of usage and distribution.                                       !
!--------------------------------------------------------------------------------------------------!

!> Contains routines to locate a value in an ascending array using bisection
module dftbp_bisect
  use dftbp_accuracy, only : dp
  implicit none


  !> Bisection driver to find a point in an array xx(:) between xx(1) and xx(size(xx)) such that
  !> element indexed j is less than the value x queried
  interface bisection
    module procedure bisection_real
    module procedure bisection_int
  end interface bisection

contains


  !> real case for bisection search
  subroutine bisection_real(j,xx,x, tol)

    !> located element such that xx(j) < x < xx(j+1)
    integer, intent(out) :: j

    !> array of values in monotonic order to search through
    real(dp), intent(in) :: xx(:)

    !> value to locate j for
    real(dp), intent(in) :: x

    !> Tolerance for equality comparision
    real(dp), intent(in), optional :: tol

    integer :: n
    integer :: jlower,jupper,jcurr
    real(dp) :: rTol      !! real tolerance
    logical :: ascending

    n = size(xx)
    if (n == 0) then
      j = 0
      return
    end if

    if (present(tol)) then
      rTol = tol
    else
      rTol = epsilon(0.0_dp)
    end if

    if (x  < xx(1) - rTol) then
      j = 0
    else if (abs(x - xx(1)) <= rTol) then
      j = 1
    else if (abs(x - xx(n)) <= rTol) then
      j = n - 1
    else if (x > xx(n) + rTol) then
      j = n
    else
      ascending = (xx(n) >=  xx(1))
      jlower=0
      jcurr=n+1
      do while ((jcurr-jlower) > 1)
        jupper=(jcurr+jlower)/2
        if (ascending .eqv. (x >= xx(jupper) + rTol)) then
          jlower=jupper
        else
          jcurr=jupper
        end if
      end do
      j=jlower
    end if
  end subroutine bisection_real


  !> integer case for bisection search
  subroutine bisection_int(j,xx,x)

    !> located element such that xx(j) < x < xx(j+1)
    integer, intent(out) :: j

    !> array of values in monotonic order to search through
    integer, intent(in) :: xx(:)

    !> value to locate j for
    integer, intent(in) :: x

    integer :: n
    integer :: jlower,jupper,jcurr

    n = size(xx)
    if (n == 0) then
      j = 0
      return
    end if

    if (x < xx(1)) then
      j = 0
    else if (x == xx(1)) then
      j = 1
    else if(x == xx(n)) then
      j = n -1
    else if(x > xx(n)) then
      j = n
    else
      jlower=0
      jcurr=n+1
      do while ((jcurr-jlower) > 1)
        jupper=(jcurr+jlower)/2
        if((xx(n).ge.xx(1)).eqv.(x.ge.xx(jupper)))then
          jlower=jupper
        else
          jcurr=jupper
        end if
      end do
      j=jlower
    end if
  end subroutine bisection_int

end module dftbp_bisect
