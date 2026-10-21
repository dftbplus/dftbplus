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
    module procedure polyInterU2_dual
  end interface polyInterUniform

contains


  !> Returns the value of a polynomial of 5th degree at xx.
  !! The polynomial is created with the following boundary conditions: Its value, 1st and 2nd
  !! derivatives are zero at x = 0 and agree with the provided values at x = dx.
  pure function poly5ToZero_real(y0, y0p, y0pp, xx, dx, invdx) result(yy)

    !> Value of the polynomial at x = dx.
    real(dp), intent(in) :: y0

    !> Value of the 1st derivative at x = dx.
    real(dp), intent(in) :: y0p

    !> Value of the 2nd derivative at x = dx.
    real(dp), intent(in) :: y0pp

    !> The point where the polynomial should be evaluated.
    real(dp), intent(in) :: xx

    !> The point, where the polynomial value and its first two derivatives should take the provided
    !! values.
    real(dp), intent(in) :: dx

    !> Reciprocal of dx.
    real(dp), intent(in) :: invdx

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


  !> Returns the value of a polynomial of 5th degree at point xx.
  !! The polynomial is created with the following boundary conditions: Its value, 1st and 2nd
  !! derivatives are zero at x = 0 and agree with the provided values at x = dx.
  !! Note, must be called with at least 2nd order derivatives of y available
  pure function poly5ToZero_dual(y, xx, dx, invdx) result(yy)

    !> Value of the polynomial (with derivatives up to at least 2nd order) at x = dx.
    type(dual_real64), intent(in) :: y

    !> The point where the polynomial should be evaluated.
    type(dual_real64), intent(in) :: xx

    !> The point where the polynomial value and its first two derivatives should take the provided
    !! values.
    real(dp), intent(in) :: dx

    !> Reciprocal of dx for efficiency.
    real(dp), intent(in) :: invdx

    !> Value of the polynomial at xx.
    type(dual_real64) :: yy

    real(dp) :: dx1, dx2, dd, ee, ff
    type(dual_real64) :: xr
    integer :: order

    order = xx%order()
    call initialize_dual(xr, order)
    call initialize_dual(yy, order)
    dx1 = y%get_derivative(1) * dx
    dx2 = y%get_derivative(2) * dx * dx
    dd =  10.0_dp * y%get_derivative(0) - 4.0_dp * dx1 + 0.5_dp * dx2
    ee = -15.0_dp * y%get_derivative(0) + 7.0_dp * dx1 - 1.0_dp * dx2
    ff =   6.0_dp * y%get_derivative(0) - 3.0_dp * dx1 + 0.5_dp * dx2
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


#:for name, vartype, xval in [("", "real(dp)", ""), ("_dual", "type(dual_real64)",&
  & "%get_derivative(0)")]

  !> Polynomial interpolation of scalar data through given points. The algorithm is based on the one
  !! in Numerical recipes, but assumes a uniform grid spacing.
  function polyInterU1${name}$(xp, yp, xx, dy) result(yy)

    !> x-coordinates of the fit points
    real(dp), intent(in) :: xp(:)

    !> y-coordinates of the fit points
    real(dp), intent(in) :: yp(:)

    !> The point where the polynomial should be calculated
    ${vartype}$, intent(in) :: xx

    !> Optional error estimate on calculated value
    ${vartype}$, intent(out), optional :: dy

    !> The value of the polynomial
    ${vartype}$ :: yy

    integer :: iCl, ii, mm, nPts
    ${vartype}$ :: cc(size(xp)), dd(size(xp)), dyy, rTmp
  #:if name == "_dual"
    integer :: order

    order = xx%order()
    call initialize_dual(yy, order)
    call initialize_dual(cc, order)
    call initialize_dual(dd, order)
    call initialize_dual(dyy, order)
    call initialize_dual(rTmp, order)
    if (present(dy)) call initialize_dual(dy, order)
  #:endif

    nPts = size(xp)

    @:ASSERT(nPts > 1)
    @:ASSERT(size(yp) == nPts)

    cc(:) = yp
    dd(:) = yp
    iCl = ceiling((xx${xval}$-xp(1)) / abs(xp(2)-xp(1)))

    yy = yp(iCl)
    iCl = iCl - 1
    do mm = 1, nPts - 1
      do ii = 1, nPts - mm
        rTmp = (dd(ii) - cc(ii+1)) / (xp(ii+mm) - xp(ii))
        cc(ii) = (xp(ii) - xx) * rTmp
        dd(ii) = (xp(ii+mm) - xx) * rTmp
      end do
      if (2 * iCl < nPts - mm) then
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

  end function polyInterU1${name}$


  !> Polynomial interpolation of vector data through given points. The algorithm is based on the one
  !! in Numerical recipes, but assumes a uniform grid spacing and interpolates a vector of values.
  function polyInterU2${name}$(xp, yp, xx, dy) result(yy)

    !> x-coordinates of the fit points
    real(dp), intent(in) :: xp(:)

    !> y-coordinates of the fit points, sized (:, nPts)
    real(dp), intent(in) :: yp(:,:)

    !> The point where the polynomial should be calculated
    ${vartype}$, intent(in) :: xx

    !> Optional error estimate on calculated value
    ${vartype}$, intent(out), optional :: dy(:)

    !> The value of the polynomial
    ${vartype}$ :: yy(size(yp,dim=1))

    integer :: iCl, ii, mm, nPts
    ${vartype}$ :: cc(size(yp,dim=1),size(xp)), dd(size(yp,dim=1),size(xp))
    ${vartype}$ :: dyy(size(yp,dim=1)), r2Tmp(size(yp,dim=1))
  #:if name == "_dual"
    integer :: order

    order = xx%order()
    call initialize_dual(yy, order)
    call initialize_dual(cc, order)
    call initialize_dual(dd, order)
    call initialize_dual(dyy, order)
    call initialize_dual(r2Tmp, order)
    if (present(dy)) call initialize_dual(dy, order)
  #:endif

    nPts = size(xp)

    @:ASSERT(nPts > 1)
    @:ASSERT(size(yp,dim=2) == nPts)

    cc(:,:) = yp
    dd(:,:) = yp
    iCl = ceiling((xx${xval}$-xp(1)) / abs(xp(2)-xp(1)))
    yy(:) = yp(:,iCl)
    iCl = iCl - 1
    do mm = 1, nPts - 1
      do ii = 1, nPts - mm
        r2Tmp(:) = (dd(:,ii) - cc(:,ii+1)) / (xp(ii+mm) - xp(ii))
        cc(:,ii) = (xp(ii) - xx) * r2Tmp
        dd(:,ii) = (xp(ii+mm) - xx) * r2Tmp
      end do
      if (2 * iCl < nPts - mm) then
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

  end function polyInterU2${name}$

#:endfor

end module dftbp_math_interpolation
