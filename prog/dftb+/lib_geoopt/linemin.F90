!--------------------------------------------------------------------------------------------------!
!  DFTB+: general package for performing fast atomistic simulations                                !
!  Copyright (C) 2017  DFTB+ developers group                                                      !
!                                                                                                  !
!  See the LICENSE file for terms of usage and distribution.                                       !
!--------------------------------------------------------------------------------------------------!

#:include 'common.fypp'

!!* Line minimization iterator for arbitary function using its gradient.
!!* @desc
!!*   The line minimization algorithm is working in the following way:
!!*   <ol>
!!*     <li>One step along the unity vector along the gradient in the starting
!!*       point (lenght of the step specified externally)</li>
!!*     <li>Intervals between previous and current point enlarged as long as
!!*       the derivative does not change sign. (Derivative in the first point
!!*       is negative.)</li>
!!*     <li>One point between the brackets is estimated by the secant method
!!*       or by bisection (see below).</li>
!!*     <li>If the derivative in the new point is not below the tolerance:
!!*       A new point is searched between the two bracketing
!!*       points and the intermediate point between the brackets with the
!!*       following methods in a fallback way.
!!*       <ol>
!!*         <li>Quadratic fit on the derivatives in the three points.
!!*           The calculated root (using Muller's method) must lie between the
!!*           intermediate point and the left (right) bracket if the derivative
!!*           in the intermediate point is greater (less) than zero.</li>
!!*         <li>Linear interpolation (secant method) between the left (right)
!!*           bracket and the intermediate point if the derivative in latter is
!!*           greater (less) than zero.</li>
!!*         <li>Bisection between the intermediate point and the left (right)
!!*           bracket and the indermediate point (depending on the sign of
!!*           derivative there).</li>
!!*       </ol>
!!*       This step is repeated as long as the projected derivative
!!*       of the function on the line is less than the given tolerance.
!!*     </li>
!!*   </ol>
module linemin
  use assert
  use accuracy
  use constants, only : goldenMeanP1
  implicit none

  private

  !!* Holds data for the line minimalizer
  type OLineMin
    private
    integer  :: nElem              !* Number of vector elements
    logical  :: tInitialized       !* If initialized
    integer  :: state              !* State of the object
    integer  :: mIter              !* Max. nr. of iterations
    integer  :: iIter              !* Nr. of performed steps
    real(dp), allocatable :: x0(:)     !* Starting point
    real(dp), allocatable :: d0(:)     !* Direction of the line
    real(dp) :: xx(2)              !* Coordinate of left and right brackets
    real(dp) :: dx(2)              !* Derivatives in the left and right brackets
    real(dp) :: xCur               !* Current position along the line
    real(dp) :: tolerance          !* Tolerance for the line derivative
    real(dp) :: maxDisp            !* Maximal displacement in any coordinate
    real(dp) :: maxX               !* Maximal displacement along the line
    real(dp) :: firstStep          !* Step length of the first step
    logical  :: tConverged         !* If converged
  end type OLineMin


  !!* Creates a line minimizer
  interface init
    module procedure LineMin_init
  end interface

  !!* Resets a line minimizer
  interface reset
    module procedure LineMin_reset
  end interface

  !!* Gets the next point for the line minimization
  interface next
    module procedure LineMin_next
  end interface

  !!* Returns the coordinate of the minimum
  interface getMinX
    module procedure LineMin_getMinX
  end interface

  !!* Returns function value in the minimum
  interface getMinY
    module procedure LineMin_getMinY
  end interface

  !!* Returns gradient in the minimum
  interface getMinGrad
    module procedure LineMin_getMinGrad
  end interface

  !!* Returns one dimensional coordinate of the minimum
  interface getMinLambda
    module procedure LineMin_getMinLambda
  end interface

  public :: OLineMin
  public :: init, reset, next, getMinX, getMinY, getMinGrad
  public :: getMinLambda

  integer, parameter :: st_1 = 1, st_2 = 2, st_3 = 3


contains


  !!* Creates a new line minimizer
  !!* @param self      Valid line minimizer instance on exit
  !!* @param nElem     Nr. of elements in the coordinate/gradient vectors
  !!* @param mIter     Nr. of maximal iterations to perform (>3)
  !!* @param tolerance Convergence criteria for the projected derivative
  !!* @param maxDisp   Maximal movement in one coordinate in one step
  subroutine LineMin_init(self, nElem, mIter, tolerance, maxDisp)
    type(OLineMin), intent(out) :: self
    integer,  intent(in) :: nElem
    integer,  intent(in) :: mIter
    real(dp), intent(in) :: tolerance
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



  !!* Resets the line minimizer
  !!* @param self Line minimizer instance
  !!* @param x0   New starting point
  !!* @param d0   New direction
  subroutine LineMin_reset(self, x0, d0, firstStep)
    type(OLineMin), intent(inout) :: self
    real(dp), intent(in) :: x0(:)
    real(dp), intent(in) :: d0(:)
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



  !!* Passes the function value and the derivative of the last point to
  !!* line minimizer and gives a new coordinate back.
  !!* @param self       Line minimizer instance
  !!* @param fx         Function value for the last returned point
  !!* @param dx         Gradient for the last returned point
  !!* @param xNew       New point to calculate
  !!* @param tConverged True, line minimization converged. The last passed
  !!*   point is the one with the lowest projected derivative.
  !!* @note Getting back tConverged = .true. can also mean, that the
  !!*   line minimization did not converge in the maximal nr. of steps.
  !!* @note When calling this subroutine the first time, function value and
  !!*   gradient for the starting point of the minimization should be passed.
  subroutine LineMin_next(self, fx, dx, xNew, tConverged)
    type(OLineMin), intent(inout) :: self
    real(dp), intent(in)  :: fx
    real(dp), intent(in)  :: dx(:)
    real(dp), intent(out) :: xNew(:)
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



  !!* Invisible working horse for LineMin_next.
  !!* @param state       State of the minimizer
  !!* @param mIter       Nr. of maximal iterations
  !!* @param iIter       Nr. of iterations so far
  !!* @param xCur        Coordinate of the current point
  !!* @param x0          Starting point for the line minimization
  !!* @param d0          Direction of the line minimization
  !!* @param xx          Coordinates of the left and right brackets
  !!* @param dx          Derivatives in the left and right brackets
  !!* @param tConverged  If method converged
  !!* @param tolerance   Tolerance criteria for the projected derivative
  !!* @param maxX        Maximal movement along one component in one step
  !!* @param firstStep   Length of the first step along the line.
  !!* @param fu          Function value of the current point
  !!* @param du          Gradient in the current point
  !!* @param uu          Suggested coordinate of the next point on exit
  subroutine next_local(state, mIter, iIter, xCur, x0, d0, xx, dx, &
      &tConverged, tolerance, maxX, firstStep, fu, du, uu)
    integer,  intent(inout) :: state
    integer,  intent(inout) :: mIter
    integer,  intent(inout) :: iIter
    real(dp), intent(inout) :: xCur
    real(dp), intent(inout) :: x0(:)
    real(dp), intent(inout) :: d0(:)
    real(dp), intent(inout) :: xx(:), dx(:)
    logical,  intent(inout) :: tConverged
    real(dp), intent(in)    :: tolerance
    real(dp), intent(in)    :: maxX
    real(dp), intent(in)    :: firstStep
    real(dp), intent(in)    :: fu
    real(dp), intent(in)    :: du(:)
    real(dp), intent(out)   :: uu(:)

    real(dp) :: dCur, xNew
    real(dp) :: tmp, qq, aa, bb, cc
    logical  :: tDone
    integer  :: nextState

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



  !!* Gives the coordinate of the minimal point back
  !!* @param self Line minimizer
  !!* @param minX Coordinate of the minimal point
  !!* @note The passed back value is meaningless if the subroutine is called
  !!*   before the line minimizer signalizes convergence.
  subroutine LineMin_getMinX(self, minX)
    type(OLineMin), intent(in) :: self
    real(dp), intent(out) :: minX(:)

    minX(:) = self%x0(:)

  end subroutine LineMin_getMinX



  !!* Gives the function value in the minimal point back
  !!* @param self Line minimizer
  !!* @param minY Function value in the minimal point
  !!* @note The passed back value is meaningless if the subroutine is called
  !!*   before the line minimizer signalizes convergence.
  subroutine LineMin_getMinY(self, minY)
    type(OLineMin), intent(in) :: self
    real(dp), intent(out) :: minY

    minY = self%xx(1)

  end subroutine LineMin_getMinY



  !!* Gives the gradient in the minimal point back
  !!* @param self    Line minimizer
  !!* @param minGrad Gradient in the minimal point
  !!* @note The passed back value is meaningless if the subroutine is called
  !!*   before the line minimizer signalizes convergence.
  subroutine LineMin_getMinGrad(self, minGrad)
    type(OLineMin), intent(in) :: self
    real(dp), intent(out) :: minGrad(:)

    minGrad(:) = self%d0(:)

  end subroutine LineMin_getMinGrad


  !!* Returns the displacement to the minimum along the line
  !!* @param self      Line minimizer
  !!* @param minLambda Displacement along the line to the minimum
  !!* @note The passed back value is meaningless if the subroutine is called
  !!*   before the line minimizer signalizes convergence.
  subroutine LineMin_getMinLambda(self, minLambda)
    type(OLineMin), intent(in) :: self
    real(dp), intent(out) :: minLambda

    minLambda = self%xx(2)

  end subroutine LineMin_getMinLambda


end module linemin
