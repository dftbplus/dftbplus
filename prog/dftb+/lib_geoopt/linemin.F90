!--------------------------------------------------------------------------------------------------!
!  DFTB+: general package for performing fast atomistic simulations                                !
!  Copyright (C) 2006 - 2020  DFTB+ developers group                                               !
!                                                                                                  !
!  See the LICENSE file for terms of usage and distribution.                                       !
!--------------------------------------------------------------------------------------------------!

#:include 'common.fypp'

!> Line minimization iterator for arbitary function using its gradient.
!>
!> The line minimization algorithm works in the following way:
!> 1. One step along the unity vector along the gradient in the starting point (lenght of the step
!>    specified externally)
!> 2. Intervals between previous and current point enlarged as long as the derivative does not
!>    change sign. (Derivative in the first point is negative.)
!> 3. One point between the brackets is estimated by the secant method or by bisection (see below).
!> 4. If the derivative in the new point is not below the tolerance: A new point is searched between
!>    the two bracketing points and the intermediate point between the brackets with the following
!>    methods in fallback order
!>    a. Quadratic fit on the derivatives in the three points.
!>       The calculated root (using Muller's method) must lie between the
!>       intermediate point and the left (right) bracket if the derivative
!>       in the intermediate point is greater (less) than zero.
!>    b. Linear interpolation (secant method) between the left (right)
!>       bracket and the intermediate point if the derivative in latter is
!>       greater (less) than zero.
!>    c. Bisection between the intermediate point and the left (right)
!>       bracket and the indermediate point (depending on the sign of
!>       derivative there).
!>
!> Step 4. is repeated as long as the projected derivative of the function on the line is less than
!> the given tolerance.
module dftbp_linemin
  use dftbp_assert
  use dftbp_accuracy
  use dftbp_constants, only : goldenMeanP1
  implicit none

  private


  !> Holds data for the line minimalizer
  type TLineMin
    private

    !> Number of vector elements
    integer :: nElem

    !> If initialized
    logical :: tInitialized

    !> State of the object
    integer :: state

    !> Max. nr. of iterations
    integer :: mIter

    !> Nr. of performed steps
    integer :: iIter

    !> Starting point
    real(dp), allocatable :: x0(:)

    !> Direction of the line
    real(dp), allocatable :: d0(:)

    !> Coordinate of left and right brackets
    real(dp) :: xx(2)

    !> Derivatives in the left and right brackets
    real(dp) :: dx(2)

    !> Current position along the line
    real(dp) :: xCur

    !> Tolerance for the line derivative
    real(dp) :: tolerance

    !> Maximal displacement in any coordinate
    real(dp) :: maxDisp

    !> Maximal displacement along the line
    real(dp) :: maxX

    !> Step length of the first step
    real(dp) :: firstStep

    !> If converged
    logical :: tConverged

  end type TLineMin


  !> Creates a line minimizer
  interface init
    module procedure LineMin_init
  end interface


  !> Resets a line minimizer
  interface reset
    module procedure LineMin_reset
  end interface


  !> Gets the next point for the line minimization
  interface next
    module procedure LineMin_next
  end interface


  !> Returns the coordinate of the minimum
  interface getMinX
    module procedure LineMin_getMinX
  end interface


  !> Returns function value in the minimum
  interface getMinY
    module procedure LineMin_getMinY
  end interface


  !> Returns gradient in the minimum
  interface getMinGrad
    module procedure LineMin_getMinGrad
  end interface


  !> Returns one dimensional coordinate of the minimum
  interface getMinLambda
    module procedure LineMin_getMinLambda
  end interface

  public :: TLineMin
  public :: init, reset, next, getMinX, getMinY, getMinGrad
  public :: getMinLambda


  !> Internal state of the line minimiser algorithm
  integer, parameter :: st_1 = 1, st_2 = 2, st_3 = 3

contains


  !> Creates a new line minimizer
  subroutine LineMin_init(self, nElem, mIter, tolerance, maxDisp)

    !> Valid line minimizer instance on exit
    type(TLineMin), intent(out) :: self

    !> Nr. of elements in the coordinate/gradient vectors
    integer, intent(in) :: nElem

    !> Nr. of maximal iterations to perform (>3)
    integer, intent(in) :: mIter

    !> Convergence criteria for the projected derivative
    real(dp), intent(in) :: tolerance

    !> Maximal movement in one coordinate in one step
    real(dp), intent(in) :: maxDisp

    @:ASSERT(nElem > 0)
    @:ASSERT(mIter > 3)
    @:ASSERT(tolerance > 0.0_dp)
    @:ASSERT(maxDisp > 0.0_dp)

    self%nElem = nElem
    allocate(self%x0(nElem))
    allocate(self%d0(nElem))
    self%mIter = mIter
    self%tolerance = tolerance
    self%maxDisp = maxDisp
    self%tInitialized = .false.
  end subroutine LineMin_init


  !> Resets the line minimizer
  subroutine LineMin_reset(self, x0, d0, firstStep)

    !> Line minimizer instance
    type(TLineMin), intent(inout) :: self

    !> New starting point
    real(dp), intent(in) :: x0(:)

    !> New direction
    real(dp), intent(in) :: d0(:)

    !> Length of the first step along the line.
    real(dp), intent(in) :: firstStep

    real(dp) :: tmp

    @:ASSERT(size(x0) == self%nElem)
    @:ASSERT(size(d0) == self%nElem)

    self%state = st_1
    self%iIter = 0
    self%xx(:) = 0.0_dp
    self%dx(:) = 0.0_dp
    self%x0(:) = x0(:)
    tmp = sqrt(sum(d0**2))
    self%d0(:) = d0(:) / tmp
    self%xCur = 0.0_dp
    self%firstStep = firstStep
    self%maxX = self%maxDisp / maxval(abs(self%d0))
    self%tConverged = .false.
    self%tInitialized = .true.

  end subroutine LineMin_reset


  !> Passes the function value and the derivative of the last point to line minimizer and gives a
  !> new coordinate back.
  !> Getting back tConverged = .true. can also mean that the line minimization did not converge in
  !> the maximal nr. of steps.
  !> When calling this subroutine the first time, function value and gradient for the starting point
  !> of the minimization should be passed.
  subroutine LineMin_next(self, fx, dx, xNew, tConverged)

    !> Line minimizer instance
    type(TLineMin), intent(inout) :: self

    !> Function value for the last returned point
    real(dp), intent(in) :: fx

    !> Gradient for the last returned point
    real(dp), intent(in) :: dx(:)

    !> New point to calculate
    real(dp), intent(out) :: xNew(:)

    !> True, line minimization converged. The last passed point is the one with the lowest projected
    !> derivative.
    logical,  intent(out) :: tConverged

    @:ASSERT(self%tInitialized)
    @:ASSERT(size(xNew) == self%nElem)
    @:ASSERT(size(dx) == self%nElem)

    if (.not. self%tConverged) then
      call next_local(self%state, self%mIter, self%iIter, self%xCur, &
          &self%x0, self%d0, self%xx, self%dx, self%tConverged, self%tolerance,&
          &self%maxX, self%firstStep, fx, dx, xNew)
    else
      call getMinX(self, xNew)
    end if
    tConverged = self%tConverged

  end subroutine LineMin_next


  !> Invisible workhorse for LineMin_next.
  subroutine next_local(state, mIter, iIter, xCur, x0, d0, xx, dx, &
      &tConverged, tolerance, maxX, firstStep, fu, du, uu)

    !> State of the minimizer
    integer, intent(inout) :: state

    !> Nr. of maximal iterations
    integer, intent(inout) :: mIter

    !> Nr. of iterations so far
    integer, intent(inout) :: iIter

    !> Coordinate of the current point
    real(dp), intent(inout) :: xCur

    !> Starting point for the line minimization
    real(dp), intent(inout) :: x0(:)

    !> Direction of the line minimization
    real(dp), intent(inout) :: d0(:)

    !> Coordinates of the left and right brackets
    real(dp), intent(inout) :: xx(:)

    !> Derivatives in the left and right brackets
    real(dp), intent(inout) :: dx(:)

    !> If method converged
    logical,  intent(inout) :: tConverged

    !> Tolerance criteria for the projected derivative
    real(dp), intent(in) :: tolerance

    !> Maximal movement along one component in one step
    real(dp), intent(in) :: maxX

    !> Length of the first step along the line.
    real(dp), intent(in) :: firstStep

    !> Function value of the current point
    real(dp), intent(in) :: fu

    !> Gradient in the current point
    real(dp), intent(in) :: du(:)

    !> Suggested coordinate of the next point on exit
    real(dp), intent(out) :: uu(:)

    real(dp) :: dCur, xNew
    real(dp) :: tmp, qq, aa, bb, cc
    logical :: tDone
    integer :: nextState

    @:ASSERT(size(uu) == size(x0))
    @:ASSERT(size(uu) == size(d0))
    @:ASSERT(size(uu) == size(du))
    @:ASSERT(size(xx) == 2)
    @:ASSERT(size(dx) == 2)

    iIter = iIter + 1
    !! Projected derivative
    dCur = dot_product(du, d0)

    if ((abs(dCur) < tolerance .and. state == st_3) .or. iIter > mIter) then
      !! Got a point with low derivative -> store final variables
      x0(:) = x0(:) + xCur * d0(:)
      d0(:) = du(:)
      uu(:) = x0(:)
      xx(1) = fu
      xx(2) = xCur
      tConverged = .true.
      return
    end if

    if (state == st_1) then
      !! Make sure line direction points along decreasing gradient.
      if (dCur > 0.0_dp) then
        d0(:) = -d0(:)
        dCur = -dCur
      end if
      xx(1) = xCur
      dx(1) = dCur
      xNew = firstStep
      nextState = st_2
    end if

    if (state == st_2) then
      if (dCur < 0.0_dp) then
        !! Current derivative still negative -> enhance interval & go further
        xNew = xCur + goldenMeanP1 * (xCur - xx(1))
        xx(1) = xCur
        dx(1) = dCur
        nextState = st_2
      else
        !! Minimum bracketed -> Secant interpolation for root of derivative
        state = st_3
      endif
    end if

    if (state == st_3) then
      tDone = .false.
      nextState = st_3

      !! Do quadratic interpolation if there are enough points
      if (iIter > 2 .and. abs(xx(1)-xx(2)) > epsilon(0.0_dp)) then
        qq = (xCur - xx(1)) / (xx(1) - xx(2))
        tmp = dCur - dx(1)
        cc = qq**2 * (dx(2) - dx(1))
        aa = qq * tmp + cc
        bb = 2.0_dp * qq * tmp + cc + tmp
        cc = (qq + 1.0_dp) * dCur
        tmp = bb**2 - 4 * aa * cc
        if (tmp >= 0.0_dp) then
          tmp = sqrt(tmp)
          if (bb > 0.0_dp) then
            aa = bb + tmp
          else
            aa = bb - tmp
          end if
          if (abs(aa) > epsilon(0.0_dp)) then
            xNew = xCur - (xCur - xx(1)) * 2.0_dp * cc / aa
            if ((xx(1) < xNew .eqv. xNew < xx(2))&
                &.and. (dCur < 0.0_dp .eqv. xNew > xCur)) then
              if (xNew > xCur) then
                xx(1) = xCur
                dx(1) = dCur
              else
                xx(2) = xCur
                dx(2) = dCur
              end if
              tDone = .true.
            end if
          end if
        end if
      end if

      !! If quadratic interpolation failed or not available, do linear interp.
      if (.not. tDone) then
        if (dCur < 0.0_dp) then
          xx(1) = xCur
          dx(1) = dCur
        else
          xx(2) = xCur
          dx(2) = dCur
        end if
        tmp = dx(1) - dx(2)
        if (abs(tmp) > epsilon(1.0_dp)) then
          xNew = (xx(2) * dx(1) - xx(1) * dx(2)) / tmp
          if ((xx(1) < xNew) .eqv. (xNew < xx(2))) then
            tDone = .true.
          end if
        end if
      end if

      !! If even secant method failed, do just bisection
      if (.not. tDone) then
        xNew = 0.5_dp * (xx(1) + xx(2))
      end if
    end if

    !! If proposed displacement is too large: scale it down.
    tmp = xNew - xCur
    if (abs(tmp) <= maxX) then
      xCur = xNew
    else
      xCur = xCur + sign(maxX, tmp)
    end if
    uu = x0(:) + xCur * d0(:)
    state = nextState

  end subroutine next_local


  !> Gives the coordinate of the minimal point back
  !> The value passed back is meaningless if the subroutine is called before the line minimizer
  !> signals convergence.
  subroutine LineMin_getMinX(self, minX)

    !> Line minimizer
    type(TLineMin), intent(in) :: self

    !> Coordinate of the minimal point
    real(dp), intent(out) :: minX(:)

    minX(:) = self%x0(:)

  end subroutine LineMin_getMinX


  !> Returns the function at the minimal point
  !> The value passed back is meaningless if the subroutine is called before the line minimizer
  !> signals convergence.
  subroutine LineMin_getMinY(self, minY)

    !> Line minimizer
    type(TLineMin), intent(in) :: self

    !> Function value in the minimal point
    real(dp), intent(out) :: minY

    minY = self%xx(1)

  end subroutine LineMin_getMinY


  !> Gives the gradient in the minimal point back
  !> The value passed back is meaningless if the subroutine is called before the line minimizer
  !> signals convergence.
  subroutine LineMin_getMinGrad(self, minGrad)

    !> Line minimizer
    type(TLineMin), intent(in) :: self

    !> Gradient in the minimal point
    real(dp), intent(out) :: minGrad(:)

    minGrad(:) = self%d0(:)

  end subroutine LineMin_getMinGrad


  !> Returns the displacement to the minimum along the line
  !> The value passed back is meaningless if the subroutine is called before the line minimizer
  !> signals convergence.
  subroutine LineMin_getMinLambda(self, minLambda)

    !> Line minimizer
    type(TLineMin), intent(in) :: self

    !> Displacement along the line to the minimum
    real(dp), intent(out) :: minLambda

    minLambda = self%xx(2)

  end subroutine LineMin_getMinLambda

end module dftbp_linemin
