!--------------------------------------------------------------------------------------------------!
!  DFTB+: general package for performing fast atomistic simulations                                !
!  Copyright (C) 2006 - 2019  DFTB+ developers group                                               !
!                                                                                                  !
!  See the LICENSE file for terms of usage and distribution.                                       !
!--------------------------------------------------------------------------------------------------!

#:include 'common.fypp'

!> Implements a repulsive potential between two atoms represented by cubic splines.
module dftbp_repspline
  use dftbp_assert
  use dftbp_accuracy
  use dftbp_bisect
  implicit none
  private

  public :: TRepSplineIn, ORepSpline, init
  public :: getCutoff, getEnergy, getEnergyDeriv


  !> Initialisation type for ORepSpline
  type TRepSplineIn

    !> Starting pos. for each spline
    real(dp), allocatable :: xStart(:)

    !> Spline coeffs (3, nSpline-1)
    real(dp), allocatable :: spCoeffs(:,:)

    !> Coeffs. of the last spline
    real(dp) :: spLastCoeffs(6)

    !> Coeffs for the exponential head
    real(dp) :: expCoeffs(3)

    !> Cutoff for the last spline
    real(dp) :: cutoff
  end type TRepSplineIn


  !> Contains the spline representation of a repulsive.
  type ORepSpline
    private

    !> Nr. of splines.
    integer :: nSpline

    !> Starting point for each spline
    real(dp), allocatable :: xStart(:)

    !> Spline coeffs (3, nSpline-1)
    real(dp), allocatable :: spCoeffs(:,:)

    !> Coeffs of the last spline
    real(dp) :: spLastCoeffs(6)

    !> Exponential head
    real(dp) :: expCoeffs(3)

    !> Cutoff of the last spline
    real(dp) :: cutoff

    !> Initialisation status
    logical :: tInit = .false.
  end type ORepSpline


  !> Initialises spline repulsive.
  interface init
    module procedure RepSpline_init
  end interface


  !> Returns cutoff of the repulsive.
  interface getCutoff
    module procedure RepSpline_getCutoff
  end interface


  !> Returns energy of the repulsive for a given distance.
  interface getEnergy
    module procedure RepSpline_getEnergy
  end interface


  !> Returns gradient of the repulsive for a given distance.
  interface getEnergyDeriv
    module procedure RepSpline_getEnergyDeriv
  end interface

contains


  !> Initialises spline repulsive.
  subroutine RepSpline_init(self, inp)

    !> Spline repulsive.
    type(ORepSpline), intent(out) :: self

    !> Input parameters for the spline repulsive.
    type(TRepSplineIn), intent(in) :: inp

    @:ASSERT(.not. self%tInit)
    @:ASSERT(size(inp%xStart) > 0)
    @:ASSERT(size(inp%spCoeffs, dim=1) == 4)
    @:ASSERT(size(inp%spCoeffs, dim=2) == size(inp%xStart) - 1)
    @:ASSERT(inp%cutoff >= 0.0_dp)

    self%nSpline = size(inp%xStart)
    allocate(self%xStart(self%nSpline))
    allocate(self%spCoeffs(4, self%nSpline - 1))
    self%xStart(:) = inp%xStart(:)
    self%spCoeffs(:,:) = inp%spCoeffs(:,:)
    self%spLastCoeffs(:) = inp%spLastCoeffs(:)
    self%expCoeffs(:) = inp%expCoeffs
    self%cutoff = inp%cutoff
    self%tInit = .true.

  end subroutine RepSpline_init


  !> Returns cutoff of the repulsive.
  function RepSpline_getCutoff(self) result(cutoff)

    !> Spline repulsive.
    type(ORepSpline), intent(in) :: self

    !> Cutoff.
    real(dp) :: cutoff

    cutoff = self%cutoff

  end function RepSpline_getCutoff


  !> Returns energy of the repulsive for a given distance.
  subroutine RepSpline_getEnergy(self, res, rr)

    !> Spline repulsive.
    type(ORepSpline), intent(in) :: self

    !> repulsive contribution
    real(dp), intent(out) :: res

    !> Distance between interacting atoms.
    real(dp), intent(in) :: rr

    integer :: imatch, ii
    real(dp) :: xh, xv

    ! below this distance, contributions are meaningless, as SK parameters and repulsive paramers
    ! are not defined
    if (rr < minNeighDist) then
      res = 0.0_dp
      ! hard repulsive
    elseif (rr < self%xStart(1)) then
      res = exp(-self%expCoeffs(1)*rr + self%expCoeffs(2)) + self%expCoeffs(3)
      ! beyond end of spline
    elseif (rr > self%cutoff) then
      res = 0.0_dp
    else
      ! find the point in the table to use
      call bisection(imatch, self%xStart, rr)

      xv = rr - self%xStart(imatch)
      xh = xv
      if (imatch < self%nSpline) then
        res = self%spCoeffs(1, imatch)
        do ii = 2, 4
          res = res + self%spCoeffs(ii, imatch) * xh
          xh = xh * xv
        end do
      else
        res = self%spLastCoeffs(1)
        do ii = 2, 6
          res = res + self%spLastCoeffs(ii) * xh
          xh = xh * xv
        end do
      end if
    end if

  end subroutine RepSpline_getEnergy


  !> Returns gradient of the repulsive for a given distance.
  subroutine RepSpline_getEnergyDeriv(self, xx, grad, grad2)

    !> Spline repulsive.
    type(ORepSpline), intent(in) :: self

    !> Actual vector between atoms
    real(dp), intent(in) :: xx(3)

    !> Resulting contribution
    real(dp), intent(out) :: grad(3)

    !> Second derivatives
    real(dp), intent(out), optional :: grad2(3,3)

    integer :: imatch, ii, jj
    real(dp) :: r2, rr, xh, xv, d1, d2

    r2 = sum(xx**2)
    rr = sqrt(r2)

    if (rr < minNeighDist .or. rr > self%cutoff) then
      d1 = 0.0_dp
      if (present(grad2)) then
        d2 = 0.0_dp
      end if
    elseif (rr < self%xStart(1)) then
      d1 = -self%expCoeffs(1) * exp(-self%expCoeffs(1) * rr + self%expCoeffs(2))
      if (present(grad2)) then
        d2 = self%expCoeffs(1)**2 * exp(-self%expCoeffs(1) * rr + self%expCoeffs(2))
      end if
    else
      call bisection(imatch, self%xStart, rr)
      xv = rr - self%xStart(imatch)
      d1 = 0.0_dp
      if (imatch < self%nSpline) then
        xh = 1.0_dp
        do ii = 2, 4
          d1 = d1 + real(ii-1, dp) * self%spCoeffs(ii, imatch) * xh
          xh = xh * xv
        end do
        if (present(grad2)) then
          xh = 1.0_dp
          d2 = 0.0_dp
          do ii = 3, 4
            d2 = d2 + real(ii-2, dp) * real(ii-1, dp) * self%spCoeffs(ii, imatch) * xh
            xh = xh * xv
          end do
        end if
      else
        xh = 1.0_dp
        do ii = 2, 6
          d1 = d1 + real(ii-1, dp) * self%spLastCoeffs(ii) * xh
          xh = xh * xv
        end do
        if (present(grad2)) then
          xh = 1.0_dp
          d2 = 0.0_dp
          do ii = 3, 6
            d2 = d2 + real(ii-2, dp) * real(ii-1, dp) * xh * self%spLastCoeffs(ii)
            xh = xh * xv
          end do
        end if
      end if
    end if
    grad(:) = d1 * xx(:) / rr

    if (present(grad2)) then
      grad2(:,:) = 0.0_dp
      do ii = 1, 3
        do jj = 1, 3
          grad2(ii,jj) = xx(ii)*xx(jj)
        end do
      end do
      grad2(:,:) = grad2 * ( d2 / r2 - d1 / (rr*r2) )
      do ii = 1, 3
        grad2(ii,ii) = grad2(ii,ii) + d1 / rr
      end do
    end if

  end subroutine RepSpline_getEnergyDeriv

end module dftbp_repspline
