!--------------------------------------------------------------------------------------------------!
!  DFTB+: general package for performing fast atomistic simulations                                !
!  Copyright (C) 2006 - 2021  DFTB+ developers group                                               !
!                                                                                                  !
!  See the LICENSE file for terms of usage and distribution.                                       !
!--------------------------------------------------------------------------------------------------!

#:include 'common.fypp'

!> Module that provides interpolation routines for Grids/Waveplot.
module waveplot_interp

  use dftbp_common_accuracy, only : dp

  implicit none
  private

  public :: linearInterpolation, polynomialInterpolation, trilinearInterpolation, locateByBisection

contains

  !> One dimensional linar interpolation.
  pure function linearInterpolation(lCc, uCc, ldata, udata, xx) result(interpVal)

    !> cartesian lower and upper bounds
    real(dp), intent(in) :: lCc, uCc

    !> corresponding data at outer bounds
    real(dp), intent(in) :: ldata, udata

    !> cartesian coordinate to calculate the interpolation value for
    real(dp), intent(in) :: xx

    !> weight factor for linear interpolation
    real(dp) :: xd

    !> obtained value by linear interpolation
    real(dp) :: interpVal

    xd = (xx - lCc) / (uCc - lCc)

    interpVal = ldata * (1.0_dp - xd) + udata * xd

  end function linearInterpolation


  !> Three dimensional linear interpolation. It is also used for interpolation on a non-cubic grid.
  pure function trilinearInterpolation(lCubeCorns, uCubeCorns, data, cartcoord) result(interpVal)

    !> cartesian points that build up the subcube, shape: [3]
    real(dp), intent(in) :: lCubeCorns(:), uCubeCorns(:)

    !> data corresponding to the corners, shape: [2, 2, 2]
    real(dp), intent(in) :: data(:,:,:)

    !> cartesian coordinate to calculate the interpolation value for, shape: [3]
    real(dp), intent(in) :: cartcoord(:)

    !> weight factors of cartesian directions
    real(dp) :: xd, yd, zd

    !> interpolation along edges
    real(dp) :: c00, c01, c10, c11

    !> interpolation between two faces
    real(dp) :: c0, c1

    !> obtained value by trilinear interpolation
    real(dp) :: interpVal

    xd = (cartcoord(1) - lCubeCorns(1)) / (uCubeCorns(1) - lCubeCorns(1))
    yd = (cartcoord(2) - lCubeCorns(2)) / (uCubeCorns(2) - lCubeCorns(2))
    zd = (cartcoord(3) - lCubeCorns(3)) / (uCubeCorns(3) - lCubeCorns(3))

    c00 = data(1, 1, 1) * (1.0_dp - xd) + data(2, 1, 1) * xd
    c01 = data(1, 1, 2) * (1.0_dp - xd) + data(2, 1, 2) * xd
    c10 = data(1, 2, 1) * (1.0_dp - xd) + data(2, 2, 1) * xd
    c11 = data(1, 2, 2) * (1.0_dp - xd) + data(2, 2, 2) * xd

    c0 = c00 * (1.0_dp - yd) + c10 * yd
    c1 = c01 * (1.0_dp - yd) + c11 * yd

    interpVal = c0 * (1.0_dp - zd) + c1 * zd

  end function trilinearInterpolation


  !> One dimensional polynomial interpolation
  pure function polynomialInterpolation(xx, yy, secDerivs, xExact) result(interpVal)

    !> tabulated x- and y-values with x-values in increasing or decreasing order
    real(dp), intent(in) :: xx(:), yy(:)

    !> second derivatives of y-values
    real(dp), intent(in) :: secDerivs(:)

    !> value to calculate interpolation for
    real(dp), intent(in) :: xExact

    !> cubic-polynomial interpolated value
    real(dp) :: interpVal

    !> bracket the input value xExact
    integer :: lower, upper, nVals

    !> auxiliary variables
    real(dp) :: coeff1, coeff2, delta

    nVals = size(xx)

    ! get table index, corresponding to the next value below given x-value xExact
    lower = max(min(locateByBisection(xx, xExact), nVals - 1), 1)
    upper = lower + 1
    delta = xx(upper) - xx(lower)

    coeff1 = (xx(upper) - xExact) / delta
    coeff2 = (xExact - xx(lower)) / delta

    interpVal = coeff1 * yy(lower) + coeff2 * yy(upper) + ((coeff1**3 - coeff1) *&
        & secDerivs(lower) + (coeff2**3 - coeff2) * secDerivs(upper)) * (delta**2) / 6.0_dp

  end function polynomialInterpolation


  !> Determines the elements of an array which surround a third value
  pure function locateByBisection(xx, yy) result(idx)

    !> tabulated values in increasing or decreasing order
    real(dp), intent(in) :: xx(:)

    !> value that will be bracket between xx(idx) and xx(idx + 1)
    real(dp), intent(in) :: yy

    !> index that will bracket yy between xx(idx) and xx(idx + 1)
    integer :: idx

    !> bracket the input value yExact
    integer :: nVals, lower, middle, upper

    !> True, if xx is in increasing order, otherwise False
    logical :: tAscnd

    nVals = size(xx)

    ! check if xx is in increasing or decreasing order
    tAscnd = (xx(nVals) >= xx(1))

    lower = 0
    upper = nVals + 1

    do

      if (upper - lower <= 1) exit
      middle = (upper + lower) / 2

      if (tAscnd .eqv. (yy >= xx(middle))) then
        lower = middle
      else
        upper = middle
      end if

    end do

    if (yy == xx(1)) then
      idx = 1
    else if (yy == xx(nVals)) then
      idx = nVals - 1
    else
      idx = lower
    end if

  end function locateByBisection


  pure function iMinLoc(arr)

    !> array to investigate
    real(dp), intent(in) :: arr(:)

    !> auxiliary variables
    integer :: iMin(1), iMinLoc

    iMin = minloc(arr(:))
    iMinLoc = imin(1)

  end function iMinLoc


end module waveplot_interp
