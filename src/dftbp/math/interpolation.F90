!--------------------------------------------------------------------------------------------------!
!  DFTB+: general package for performing fast atomistic simulations                                !
!  Copyright (C) 2006 - 2025  DFTB+ developers group                                               !
!                                                                                                  !
!  See the LICENSE file for terms of usage and distribution.                                       !
!--------------------------------------------------------------------------------------------------!

#:include 'common.fypp'

!> Contains routines for interpolation and extrapolation
module dftbp_math_interpolation
  use dftbp_common_accuracy, only : dp
  use dftbp_extlibs_dnaoad
  use dftbp_io_message, only : error
  implicit none

  private
  public :: freeCubicSpline, poly5ToZero, polyInterUniform

  !> Returns the value of a 1D free spline at a specified point
  interface freeCubicSpline
    module procedure freeCubicSpline_real
    module procedure freeCubicSpline_dual
  end interface freeCubicSpline

  !> Returns the value of a polynomial of 5th degree with specified boundary conditions
  interface poly5ToZero
    module procedure poly5ToZero_real
    module procedure poly5ToZero_dual
  end interface poly5ToZero

  !> Uniform grid polynomial interpolation
  interface polyInterUniform
    module procedure polyInterU1
    module procedure polyInterU2
    module procedure polyInterU1_dual
  end interface polyInterUniform

contains


  !> Returns the value of a polynomial of 5th degree at x.
  !! The polynomial is created with the following boundary conditions: Its value, 1st and 2nd
  !! derivatives are zero at x = 0 and agree with the provided values at x = dx.
  pure function poly5ToZero_real(y0, y0p, y0pp, xx, dx, invdx) result(yy)

    !> Value of the polynomial at x = dx.
    real(dp), intent(in) :: y0

    !> Value of the 1st derivative at x = dx.
    real(dp), intent(in) :: y0p

    !> Value of the 2nd derivative at x = dx.
    real(dp), intent(in) :: y0pp

    !> Reciprocal of dx.
    real(dp), intent(in) :: invdx

    !> The point where the polynomial should be evaluated.
    real(dp), intent(in) :: xx

    !> The point, where the polynomial value and its first two derivatives should take the provided
    !! values.
    real(dp), intent(in) :: dx

    !> Value of the polynomial at xx.
    real(dp) :: yy

    real(dp) :: dx1, dx2, dd, ee, ff, xr

    dx1 = y0p * dx
    dx2 = y0pp * dx * dx
    dd =  10.0_dp * y0 - 4.0_dp * dx1 + 0.5_dp * dx2
    ee = -15.0_dp * y0 + 7.0_dp * dx1 - 1.0_dp * dx2
    ff =   6.0_dp * y0 - 3.0_dp * dx1 + 0.5_dp * dx2
    xr = xx * invdx
    yy = ((ff * xr + ee) * xr + dd) * xr * xr * xr

  end function poly5ToZero_real


  !> Returns the value of a polynomial of 5th degree at x.
  !! The polynomial is created with the following boundary conditions: Its value, 1st and 2nd
  !! derivatives are zero at x = 0 and agree with the provided values at x = dx.
  pure function poly5ToZero_dual(y0, y0p, y0pp, xx, dx, invdx) result(yy)

    !> Value of the polynomial at x = dx.
    real(dp), intent(in) :: y0

    !> Value of the 1st derivative at x = dx.
    real(dp), intent(in) :: y0p

    !> Value of the 2nd derivative at x = dx.
    real(dp), intent(in) :: y0pp

    !> Reciprocal of dx.
    real(dp), intent(in) :: invdx

    !> The point where the polynomial should be evaluated.
    type(dual_real64), intent(in) :: xx

    !> The point where the polynomial value and its first two derivatives should take the provided
    !! values.
    real(dp), intent(in) :: dx

    !> Value of the polynomial at xx.
    type(dual_real64) :: yy

    real(dp) :: dx1, dx2, dd, ee, ff
    type(dual_real64) :: xr

    dx1 = y0p * dx
    dx2 = y0pp * dx * dx
    dd =  10.0_dp * y0 - 4.0_dp * dx1 + 0.5_dp * dx2
    ee = -15.0_dp * y0 + 7.0_dp * dx1 - 1.0_dp * dx2
    ff =   6.0_dp * y0 - 3.0_dp * dx1 + 0.5_dp * dx2
    xr = xx * invdx
    yy = ((ff * xr + ee) * xr + dd) * xr * xr * xr

  end function poly5ToZero_dual


  !> Returns the value of a free spline at a specified point.
  !! The spline is created with the following boundary conditions: Its value, 1st and 2nd
  !! derivatives agree with the provided values at x = 0 and its value agrees with the provided
  !! value at x = dx.
  !! Note: If you want the value for a derivative, you have to query both the first and second
  !! derivatives.
  subroutine freeCubicSpline_real(y0, y0p, y0pp, dx, ydx, xx, yy, yp, ypp)

    !> Function value at x = 0.
    real(dp), intent(in) :: y0

    !> First derivative at x = 0.
    real(dp), intent(in) :: y0p

    !> Second derivative at x = 0.
    real(dp), intent(in) :: y0pp

    !> Second fitting point.
    real(dp), intent(in) :: ydx

    !> Function value at dx.
    real(dp), intent(in) :: dx

    !> Point to interpolate.
    real(dp), intent(in) :: xx

    !> Value of the 3rd order polynomial at xx.
    real(dp), intent(out), optional :: yy

    !> First derivative at xx.
    real(dp), intent(out), optional :: yp

    !> Second derivative at xx.
    real(dp), intent(out), optional :: ypp

    real(dp) :: aa, bb, cc, dd, dx1

    @:ASSERT(present(yp) .eqv. present(ypp))

    aa = y0
    bb = y0p
    cc = 0.5_dp * y0pp
    dx1 = 1.0_dp / dx
    dd = (((ydx - y0)*dx1 - y0p)*dx1 - 0.5_dp*y0pp)*dx1
    if (present(yy)) then
      yy = ((dd*xx + cc)*xx + bb)*xx + aa
    end if
    if (present(yp)) then
      yp = (3.0_dp*dd*xx + 2.0_dp*cc)*xx + bb
      ypp = 6.0_dp * dd * xx + 2.0_dp * cc
    end if

  end subroutine freeCubicSpline_real


  !> Returns the value of a free spline at a specified point.
  !! The spline is created with the following boundary conditions: Its value, 1st and 2nd
  !! derivatives agree with the provided values at x = 0 and its value agrees with the provided
  !! value at x = dx.
  subroutine freeCubicSpline_dual(y0, y0p, y0pp, dx, ydx, xx, yy)

    !> Function value at x = 0.
    real(dp), intent(in) :: y0

    !> First derivative at x = 0.
    real(dp), intent(in) :: y0p

    !> Second derivative at x = 0.
    real(dp), intent(in) :: y0pp

    !> Second fitting point.
    real(dp), intent(in) :: ydx

    !> Function value at dx.
    real(dp), intent(in) :: dx

    !> Point at which to interpolate.
    type(dual_real64), intent(in) :: xx

    !> 3rd order polynomial at point xx.
    type(dual_real64), intent(out) :: yy

    real(dp) :: aa, bb, cc, dd, dx1

    aa = y0
    bb = y0p
    cc = 0.5_dp * y0pp
    dx1 = 1.0_dp / dx
    dd = (((ydx - y0) * dx1 - y0p) * dx1 - 0.5_dp * y0pp) * dx1
    yy = ((dd * xx + cc) * xx + bb) * xx + aa

  end subroutine freeCubicSpline_dual


  !> Polynomial interpolation of scalar data through given points. The algorithm is based on the one
  !! in Numerical recipes, but assumes a uniform grid spacing.
  function polyInterU1(xp, yp, xx, dy) result(yy)

    !> x-coordinates of the fit points
    real(dp), intent(in) :: xp(:)

    !> y-coordinates of the fit points
    real(dp), intent(in) :: yp(:)

    !> The point where the polynomial should be calculated
    real(dp), intent(in) :: xx

    !> Optional error estimate on calculated value
    real(dp), intent(out), optional :: dy

    !> The value of the polynomial
    real(dp) :: yy

    integer :: iCl, ii, mm, nn
    real(dp) :: cc(size(xp)), dd(size(xp)), dyy, rTmp

    nn = size(xp)

    @:ASSERT(nn > 1)
    @:ASSERT(size(yp) == nn)

    cc(:) = yp
    dd(:) = yp
    iCl = ceiling((xx-xp(1))/abs(xp(2)-xp(1)))
    yy = yp(iCl)
    iCl = iCl - 1
    do mm = 1, nn - 1
      do ii = 1, nn - mm
        rTmp = xp(ii) - xp(ii+mm)
        rTmp = (cc(ii+1) - dd(ii)) / rTmp
        cc(ii) = (xp(ii) - xx) * rTmp
        dd(ii) = (xp(ii+mm) - xx) * rTmp
      end do
      if (2 * iCl < nn - mm) then
        dyy = cc(iCl + 1)
      else
        dyy = dd(iCl)
        iCl = iCl - 1
      end if
      yy = yy + dyy
    end do

    if (present(dy)) then
      dy = dyy
    end if

  end function polyInterU1


  !> Polynomial interpolation of vector data through given points. The algorithm is based on the one
  !! in Numerical recipes, but assumes a uniform grid spacing and interpolates a vector of values.
  function polyInterU2(xp, yp, xx, dy) result(yy)

    !> x-coordinates of the fit points
    real(dp), intent(in) :: xp(:)

    !> y-coordinates of the fit points
    real(dp), intent(in) :: yp(:,:)

    !> The point where the polynomial should be calculated
    real(dp), intent(in) :: xx

    !> Optional error estimate on calculated value
    real(dp), intent(out), optional :: dy(:)

    !> The value of the polynomial
    real(dp) :: yy(size(yp,dim=1))

    integer :: iCl, ii, mm, nn
    real(dp) :: cc(size(yp,dim=1),size(xp)), dd(size(yp,dim=1),size(xp))
    real(dp) :: dyy(size(yp,dim=1)), r2Tmp(size(yp,dim=1)), delta(size(xp)-1)

    nn = size(xp)

    @:ASSERT(nn > 1)
    @:ASSERT(size(yp,dim=2) == nn)

    delta(:) = 0.0_dp
    do ii = 1, size(xp)-1
      delta(ii) = 1.0_dp / (xp(1+ii) - xp(1))
    end do
    cc(:,:) = yp
    dd(:,:) = yp
    iCl = ceiling((xx-xp(1))*delta(1))
    yy(:) = yp(:,iCl)
    iCl = iCl - 1
    do mm = 1, nn - 1
      do ii = 1, nn - mm
        r2Tmp(:) = (dd(:,ii) - cc(:,ii+1)) * delta(mm)
        cc(:,ii) = (xp(ii) - xx) * r2Tmp
        dd(:,ii) = (xp(ii+mm) - xx) * r2Tmp
      end do
      if (2 * iCl < nn - mm) then
        dyy(:) = cc(:,iCl + 1)
      else
        dyy(:) = dd(:,iCl)
        iCl = iCl - 1
      end if
      yy(:) = yy + dyy
    end do

    if (present(dy)) then
      dy(:) = dyy
    end if

  end function polyInterU2


  !> Polynomial interpolation of scalar data through given points. The algorithm is based on the one
  !! in Numerical recipes, but assumes a uniform grid spacing.
  function polyInterU1_dual(xp, yp, xx, dy) result(yy)

    !> x-coordinates of the fit points
    real(dp), intent(in) :: xp(:)

    !> y-coordinates of the fit points
    real(dp), intent(in) :: yp(:)

    !> The point where the polynomial should be calculated
    type(dual_real64), intent(in) :: xx

    !> Optional error estimate on calculated value
    real(dp), intent(out), optional :: dy

    !> The value of the polynomial
    type(dual_real64) :: yy

    integer :: iCl, ii, mm, nn, order
    type(dual_real64) :: cc(size(xp)), dd(size(xp)), dyy, rTmp

    order = size(xx%f) - 1
    call initialize_dual(yy, order)
    call initialize_dual(cc, order)
    call initialize_dual(dd, order)
    call initialize_dual(dyy, order)
    call initialize_dual(rTmp, order)

    nn = size(xp)

    @:ASSERT(nn > 1)
    @:ASSERT(size(yp) == nn)

    cc(:) = yp
    dd(:) = yp
    iCl = ceiling((xx%f(0)-xp(1))/abs(xp(2)-xp(1)))
    yy = yp(iCl)
    iCl = iCl - 1
    do mm = 1, nn - 1
      do ii = 1, nn - mm
        rTmp = xp(ii) - xp(ii+mm)
        rTmp = (cc(ii+1) - dd(ii)) / rTmp
        cc(ii) = (xp(ii) - xx) * rTmp
        dd(ii) = (xp(ii+mm) - xx) * rTmp
      end do
      if (2 * iCl < nn - mm) then
        dyy = cc(iCl + 1)
      else
        dyy = dd(iCl)
        iCl = iCl - 1
      end if
      yy = yy + dyy
    end do

    if (present(dy)) then
      dy = dyy%f(0)
    end if

  end function polyInterU1_dual

end module dftbp_math_interpolation
