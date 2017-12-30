!--------------------------------------------------------------------------------------------------!
!  DFTB+: general package for performing fast atomistic simulations                                !
!  Copyright (C) 2017  DFTB+ developers group                                                      !
!                                                                                                  !
!  See the LICENSE file for terms of usage and distribution.                                       !
!--------------------------------------------------------------------------------------------------!

#:include 'common.fypp'

!> Contains routines for interpolation and extrapolation
module interpolation
  use assert
  use accuracy
  use message
  implicit none
  private

  public :: poly5ToZero, freeCubicSpline, polyInter, polyInterUniform

  !> Uniform grid polynomial interpolation
  interface polyInterUniform
    module procedure polyInterU1
    module procedure polyInterU2
  end interface polyInterUniform

contains


  !> Returns the value of a polynomial of 5th degree at x.
  !> The polynomial is created with the following boundary conditions: Its value, and its 1st and
  !> 2nd derivatives are zero at x = 0 and agree with the provided values at x = dx.
  function poly5ToZero(y0, y0p, y0pp, xx, dx) result(yy)

    !> Value of the polynom at x = dx.
    real(dp), intent(in) :: y0

    !> Value of the 1st derivative at x = dx.
    real(dp), intent(in) :: y0p

    !> Value of the 2nd derivative at x = dx.
    real(dp), intent(in) :: y0pp

    !> The point where the polynomial should be calculated
    real(dp), intent(in) :: xx

    !> The point, where the polynomials value and first two derivatives should take the provided
    !> values.
    real(dp), intent(in) :: dx


    !> Value of the polynomial at xx.
    real(dp) :: yy

    real(dp) :: dx1, dx2, dd, ee, ff, xr

    dx1 = y0p * dx
    dx2 = y0pp * dx * dx
    dd =  10.0_dp * y0 - 4.0_dp * dx1 + 0.5_dp * dx2
    ee = -15.0_dp * y0 + 7.0_dp * dx1 - 1.0_dp * dx2
    ff =   6.0_dp * y0 - 3.0_dp * dx1 + 0.5_dp * dx2
    xr = xx / dx
    yy = ((ff*xr + ee)*xr + dd)*xr*xr*xr

  end function poly5ToZero


  !> Returns the value of a free spline at a certain point.
  !> The spline is created with the following boundary conditions: Its value, 1st and 2nd
  !> derivatives agree with the provided values at x = 0 and its value agrees with the provided
  !> value at x = dx.  Note: If you want the value for a derivative, you have to query them both.
  subroutine freeCubicSpline(y0, y0p, y0pp, dx, ydx, xx, yy, yp, ypp)

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

  end subroutine freeCubicSpline


  !> Polynomial interpolation through given points The algorithm is based on the one in Numerical
  !> recipes, but assumes a uniform grid spacing.
  function polyInterU1(xp, yp, xx, dy) result(yy)

    !> x-coordinates of the fit points
    real(dp), intent(in) :: xp(:)

    !> y-coordinates of the fit points
    real(dp), intent(in) :: yp(:)

    !> The point, where the polynomial should be calculated
    real(dp), intent(in) :: xx

    !> Optional error estimate on calculated value
    real(dp), intent(out), optional :: dy

    !> The value of the polynomial
    real(dp) :: yy

    integer :: nn
    integer :: iCl, ii, mm

    real(dp) :: cc(size(xp)), dd(size(xp))
    real(dp) :: dyy, rTmp

    nn = size(xp)

    @:ASSERT(nn > 1)
    @:ASSERT(size(yp) == nn)

    cc(:) = yp(:)
    dd(:) = yp(:)
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

  !> Polynomial interpolation through given points The algorithm is based on the one in Numerical
  !> recipes, but assumes a uniform grid spacing and interpolates a vector of values.
  function polyInterU2(xp, yp, xx, dy) result(yy)

    !> x-coordinates of the fit points
    real(dp), intent(in) :: xp(:)

    !> y-coordinates of the fit points
    real(dp), intent(in) :: yp(:,:)

    !> The point, where the polynomial should be calculated
    real(dp), intent(in) :: xx

    !> Optional error estimate on calculated value
    real(dp), intent(out), optional :: dy(:)

    !> The value of the polynomial
    real(dp) :: yy(size(yp,dim=1))

    integer :: nn
    integer :: iCl, ii, mm

    real(dp) :: cc(size(yp,dim=1),size(xp)), dd(size(yp,dim=1),size(xp))
    real(dp) :: dyy(size(yp,dim=1)), r2Tmp(size(yp,dim=1))
    real(dp) :: delta(size(xp)-1)

    nn = size(xp)

    @:ASSERT(nn > 1)
    @:ASSERT(size(yp,dim=2) == nn)

    delta = 0.0_dp
    do ii = 1, size(xp)-1
      delta(ii) = 1.0_dp / (xp(1+ii) - xp(1))
    end do
    cc = yp
    dd = yp
    iCl = ceiling((xx-xp(1))*delta(1))
    yy = yp(:,iCl)
    iCl = iCl - 1
    do mm = 1, nn - 1
      do ii = 1, nn - mm
        r2Tmp = (dd(:,ii) - cc(:,ii+1)) * delta(mm)
        cc(:,ii) = (xp(ii) - xx) * r2Tmp
        dd(:,ii) = (xp(ii+mm) - xx) * r2Tmp
      end do
      if (2 * iCl < nn - mm) then
        dyy = cc(:,iCl + 1)
      else
        dyy = dd(:,iCl)
        iCl = iCl - 1
      end if
      yy = yy + dyy
    end do

    if (present(dy)) then
      dy = dyy
    end if

  end function polyInterU2


  !> Polynomial interpolation through given points
  !> The algorithm is based on the one in Numerical recipes.
  function polyInter(xp, yp, xx, dy) result(yy)

    !> x-coordinates of the fit points
    real(dp), intent(in) :: xp(:)

    !> y-coordinates of the fit points
    real(dp), intent(in) :: yp(:)

    !> The point, where the polynomial should be calculated
    real(dp), intent(in) :: xx

    !> Optional error estimate on calculated value
    real(dp), intent(out), optional :: dy

    !> The value of the polynomial
    real(dp) :: yy

    integer :: nn
    integer :: iCl, ii, mm

    real(dp) :: cc(size(xp)), dd(size(xp))
    real(dp) :: dx, dxNew, dyy, rTmp

    nn = size(xp)

    @:ASSERT(nn > 1)
    @:ASSERT(size(yp) == nn)

    cc(:) = yp(:)
    dd(:) = yp(:)
    iCl = 1
    dx = abs(xx - xp(iCl))
    do ii = 2, nn
      dxNew = abs(xx - xp(ii))
      if (dxNew < dx) then
        iCl = ii
        dx = dxNew
      end if
    end do
    yy = yp(iCl)
    iCl = iCl - 1
    do mm = 1, nn - 1
      do ii = 1, nn - mm
        rTmp = xp(ii) - xp(ii+mm)
        if (abs(rTmp) < epsilon(1.0_dp)) then
          call error("Polint failed")
        end if
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

  end function polyInter

end module interpolation
