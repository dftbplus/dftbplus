!--------------------------------------------------------------------------------------------------!
!  DFTB+: general package for performing fast atomistic simulations                                !
!  Copyright (C) 2017  DFTB+ developers group                                                      !
!                                                                                                  !
!  See the LICENSE file for terms of usage and distribution.                                       !
!--------------------------------------------------------------------------------------------------!

#:include 'common.fypp'

!> Function minimization with lbfgs algorithm
!> References:
!> Nocedal - Updating Quasi-Newton Matrices with Limited Storage (1980), Mathematics of Computation
!> 35, pp. 773-782.
!>
!> Nocedal, Wright - Numerical Optimization
module lbfgs
  use accuracy
  use assert
  use ieee_arithmetic
  use message
  use linemin
  implicit none
  private

  public :: Tlbfgs, init

  !> Holds data for the line minimalizer used by this implementation
  type OLineSearch
    private

    !> Elements of vector x
    integer :: nElem

    !> Number of current iteration
    integer :: iter

    !> maximum number of iterations
    integer :: mIter

    !> Search direction
    real(dp), allocatable :: d0(:)

    !> initial point
    real(dp), allocatable :: x0(:)

    !> New point
    real(dp), allocatable :: xNew(:)

    !> Gradient at initial point
    real(dp), allocatable :: dx_0(:)

    !> Gradient of fx at low value
    real(dp), allocatable :: dx_lo(:)

    !> Gradient of fx at high value
    real(dp), allocatable :: dx_hi(:)

    !> line search position on the high side
    real(dp) :: alpha_hi

    !> line search position on the low side of the gradient
    real(dp) :: alpha_lo

    !> new line search position as fraction of gradient
    real(dp) :: alpha_New

    !> phi(alpha) = f(x + alpha*d0)
    real(dp) :: phi_0

    !> Function on the high side
    real(dp) :: phi_hi

    !> Function on the low side
    real(dp) :: phi_lo

    !> Gradient of phi at initial point
    real(dp) :: dphi_0

    !> Gradient of phi on the high side
    real(dp) :: dphi_hi

    !> gradient of phi on the low side
    real(dp) :: dphi_lo

    !> maximum scale of step along gradient
    real(dp) :: maxAlpha

    !> Maximal displacement in one step
    real(dp) :: maxDisp

    !> True if program is in Zoom state
    logical :: isZoom

    !> True if converged
    logical :: tConverged

    !> Holding variable for phi between steps
    real(dp) :: phi_temp

    !> holding variable for alpha between steps
    real(dp) :: alpha_temp

    !> tolerance for line search gradient termination
    real(dp) :: tol

  end type OLineSearch

  !> Calculates search direction via lbfgs-Algorithm
  type :: Tlbfgs
    private

    !> line minimizer
    type(OLineSearch) :: lineSearch

    !> Number of elements
    integer :: nElem

    !> current point
    real(dp), allocatable :: xx(:)

    !> line search position
    real(dp) :: alpha

    !> new direction to search along
    logical :: needNewDirection

    !> has the linesearch converged
    logical :: tLineConverged

    !> Number of past iterations to save
    integer :: mem

    !> Current iteration
    integer :: iter

    !> Point
    real(dp), allocatable :: xOld(:)

    !> Gradient at point
    real(dp), allocatable :: gg(:)

    !> s_i = x_i+1 - x_i
    real(dp), allocatable :: ss(:,:)

    !> y_i = g_i+1 - g_i
    real(dp), allocatable :: yy(:,:)

    !> rho_i = 1 / (y_i * s_i)
    real(dp), allocatable :: rho(:)

    !> Search direction
    real(dp), allocatable :: dir(:)

    !> tolerance for gradient
    real(dp) :: tol

  contains

    !> Next structure from LBFGS
    procedure :: next_lbfgs

    !> reset the search
    procedure :: reset_lbfgs

  end type Tlbfgs

  !> Initialise the LBFGS minimiser
  interface init
    module procedure lbfgs_init
  end interface init

  !> Check for cubic interpolation
  real(dp), parameter :: delta1 = 0.2_dp

  !> Check for quadratic interpolation
  real(dp), parameter :: delta2 = 0.1_dp

  !> Lower Wolfe condition parameter of LBFGS
  real(dp), parameter :: WOLFE1 = 1e-4_dp

  !> Upper Wolfe condition parameter for LBFGS
  real(dp), parameter :: WOLFE2 = 0.9_dp


contains

  !> Initialize lbfgs instance
  subroutine lbfgs_init(self, nElem, tol, maxDisp, mem)

    !> Dimension of vectors
    integer, intent(in) :: nElem

    !> Tolerance for gradient
    real(dp), intent(in) :: tol

    !> Maximal displacement in one step
    real(dp), intent(in) :: maxDisp

    !> Number of past iterations which will be saved
    integer, intent(in) :: mem

    !> lbfgs instance on exit
    type(Tlbfgs), intent(out) :: self

    @:ASSERT(nElem > 0)
    @:ASSERT(tol > 0.0_dp)
    @:ASSERT(maxDisp > 0.0_dp)
    @:ASSERT(mem > 0)

    self%nElem = nElem
    self%mem = mem
    self%tol = tol
    allocate(self%xx(nElem))
    allocate(self%xOld(nElem))
    allocate(self%gg(nElem))
    allocate(self%ss(nElem, mem))
    allocate(self%yy(nElem, mem))
    allocate(self%rho(mem))
    allocate(self%dir(nElem))

    call LineMin_init(self%lineSearch, nElem, 10, tol, maxDisp)

  end subroutine lbfgs_init

  !> Passes calculated function value and gradient to the minimizer and
  !> gives a new coordinate back.
  subroutine next_lbfgs(self, fx, dx, xNew, tConverged)

    !> LBFGS Instance.
    class(Tlbfgs), intent(inout) :: self

    !> Function value in the last point
    real(dp), intent(in)  :: fx

    !> Gradient in the last point
    real(dp), intent(in)  :: dx(:)

    !> New proposed point
    real(dp), intent(out) :: xNew(:)

    !> True, if gradient got below the specified tolerance.
    logical,  intent(out) :: tConverged

    real(dp) :: dx_temp(self%nElem)

    @:ASSERT(size(xNew) == self%nElem)
    @:ASSERT(size(dx) == self%nElem)

    if (maxval(abs(dx)) < self%tol) then
      tConverged = .true.
      xNew = self%xx
      return
    else
      tConverged = .false.
    end if
    dx_temp = dx

    ! Line Search
    if (.not. self%needNewDirection) then
      call LineMin_next(self%lineSearch, fx, dx, self%xx, self%tLineConverged)
      if (self%tLineConverged) then
        self%needNewDirection = .true.
        xNew = self%xx
        call LineMin_getMinGrad(self%lineSearch, dx_temp)
      else
        xNew = self%xx
        return
      end if
    end if
    ! Calculate new search direction
    call calcDirection_lbfgs(self, self%xx, dx_temp)
    self%alpha = 1.0_dp

    call LineMin_reset(self%lineSearch, self%xx, self%dir, self%alpha)
    self%tLineConverged = .false.
    self%needNewDirection = .false.
    call LineMin_next(self%lineSearch, fx, dx, self%xx, self%tLineConverged)
    xNew = self%xx

  end subroutine next_lbfgs


  !> Resets lbfgs minimizer
  subroutine reset_lbfgs(self, xx)

    !> lbfgs minimizer
    class(Tlbfgs), intent(inout) :: self

    !> Start point
    real(dp), intent(in) :: xx(self%nElem)

    @:ASSERT(size(xx) == self%nElem)

    self%needNewDirection = .true.
    self%tLineConverged = .false.
    self%alpha = 1.0_dp
    self%xx = xx
    self%xOld = xx
    self%gg = 0.0_dp
    self%ss = 0.0_dp
    self%yy = 0.0_dp
    self%rho = 0.0_dp
    self%dir = 0.0_dp
    self%iter = 0
  end subroutine reset_lbfgs



  !> Calculates next search direction
  subroutine calcDirection_lbfgs(self, xx, gg, diag_in)

    !> lbfgs instance
    type(Tlbfgs) :: self

    !> Point
    real(dp), intent(in) :: xx(self%nElem)

    !> Gradient at point
    real(dp), intent(in) :: gg(self%nElem)

    !> Initial diagonal hessian matrix (default identity)
    real(dp), intent(in), optional :: diag_in(self%nElem)

    integer :: ii

    ! Element where new data of this iteration will be stored
    integer :: newElem
    ! Oldest Element
    integer :: oldestElem
    ! Initial hesse matrix
    real(dp) :: diag(self%nElem)

    ! Workspace variables
    real(dp) :: qq(self%nElem), alpha(self%mem), beta

    ! First Run
    if (self%iter == 0) then
      self%gg = gg
      self%dir = -gg
      self%iter = self%iter + 1
      return
    end if

    ! Checks whether new x is different to last iteration
    if ( maxval(abs(xx - self%xOld)) < epsilon(1.0_dp) ) then
      call error("Error: x need to be different in each iteration. (LBFGS Minimizer)")
    end if

    ! Checks whether new gradient is different to last iteration
    if ( maxval(abs(gg - self%gg)) < epsilon(1.0_dp) ) then
      call error("Error: g need to be different in each iteration. (LBFGS Minimizer)")
    end if

    ! Save new ss, yy and rho
    newElem = modulo(self%iter - 1, self%mem) + 1
    self%ss(:,newElem) = xx - self%xOld
    self%yy(:,newElem) = gg - self%gg

    beta = dot_product(self%yy(:,newElem), self%ss(:,newElem))

    ! y*s needs to be positive to get a positive definite hermitian matrix
    if (beta <= epsilon(1.0_dp) ) then
      call warning("Bad step. y*s <= 0 (LBFGS Minimizer)")
    end if
    self%rho(newElem) = 1.0_dp / beta

    ! Save new x and g
    self%xOld = xx
    self%gg = gg

    ! Compute guessed hessian matrix if not set by the user
    if (.not. present(diag_in)) then
      if (self%iter == 1) then
        diag = 1.0_dp
      else
        diag = dot_product(self%ss(:, newElem), self%yy(:, newElem))&
            & / dot_product(self%yy(:, newElem), self%yy(:, newElem))
      end if

      ! If hessian set by the user, do input check
    else if (any(diag_in < 0.0_dp)) then
      call error("Error: Bad Input (diag). (LBFGS Minimizer)")
    else
      diag = diag_in
    end if

    ! Compute d = -H*g using the formula given by Nocedal. Result saved in self%dir
    qq = gg
    ii = newElem
    do
      alpha(ii) = self%rho(ii)*dot_product(self%ss(:,ii), qq)
      qq = qq - alpha(ii) * self%yy(:,ii)

      ii = ii - 1
      if (ii == 0) then
        ii = self%mem
      end if

      if (ii .eq. newElem) then
        exit
      end if

    end do

    self%dir = diag*qq

    if (newElem + 1 > self%mem) then
      oldestElem = 1
    else
      oldestElem = newElem + 1
    end if
    ii = oldestElem
    do
      beta = self%rho(ii) * dot_product(self%yy(:,ii), self%dir)
      self%dir = self%dir + self%ss(:,ii) * (alpha(ii) - beta)

      ii = ii + 1
      if (ii > self%mem) then
        ii = 1
      end if

      if (ii .eq. oldestElem) then
        exit
      end if

    end do

    self%iter = self%iter + 1

    ! Change to correct sign
    self%dir = -self%dir

  end subroutine calcDirection_lbfgs

  !> Creates a new line minimizer
  subroutine LineMin_init(self, nElem, mIter, tolerance, maxDisp)

    !> Valid line minimizer instance on exit
    type(OLineSearch), intent(out) :: self

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
    @:ASSERT(maxDisp > 0.0_dp)

    self%nElem = nElem
    self%maxDisp = maxDisp
    self%mIter = mIter

    self%tol = tolerance

    allocate(self%d0(nElem))
    allocate(self%dx_lo(nElem))
    allocate(self%dx_hi(nElem))
    allocate(self%xNew(nElem))
    allocate(self%x0(nElem))
    allocate(self%dx_0(nElem))

  end subroutine LineMin_init


  !> Resets the line minimizer
  subroutine LineMin_reset(self, x0, d0, firstStep)

    !> Line minimizer instance
    type(OLineSearch), intent(inout) :: self

    !> New starting point
    real(dp), intent(in) :: x0(:)

    !> New direction
    real(dp), intent(in) :: d0(:)

    real(dp), intent(in) :: firstStep

    @:ASSERT(size(x0) == self%nElem)
    @:ASSERT(size(d0) == self%nElem)
    @:ASSERT(all(abs(d0) > 0.0_dp))

    self%d0 = d0
    self%phi_0 = 0.0_dp
    self%dphi_0 = 0.0_dp

    self%iter = 0
    self%isZoom = .false.
    self%tConverged = .false.

    self%alpha_lo = 0.0_dp
    self%alpha_hi = 0.0_dp

    self%phi_lo = 0.0_dp
    self%phi_hi = 0.0_dp
    self%dphi_lo = 0.0_dp
    self%dphi_hi = 0.0_dp

    self%dx_0 = 0.0_dp
    self%dx_lo = 0.0_dp
    self%dx_hi = 0.0_dp

    self%x0 = x0
    self%xNew = 0.0_dp

    self%alpha_new = firstStep
    self%maxAlpha = self%maxDisp / maxval(abs(d0))

    self%alpha_temp = 0.0_dp
    self%phi_temp = 0.0_dp

  end subroutine LineMin_reset


  !> Passes the function value and the derivative of the last point to line minimizer and gives a
  !> new coordinate back.
  !> Getting back tConverged = .true. can also mean that the line minimization did not converge in
  !> the maximal nr. of steps.
  !> When calling this subroutine the first time, function value and gradient for the starting point
  !> of the minimization should be passed.
  subroutine LineMin_next(self, fx, dx, xNew, tConverged)

    !> Line minimizer instance
    type(OLineSearch), intent(inout) :: self

    !> Function value for the last returned point
    real(dp), intent(in) :: fx

    !> Gradient for the last returned point
    real(dp), intent(in) :: dx(:)

    !> New point to calculate
    real(dp), intent(out) :: xNew(:)

    !> True, line minimization converged. The last passed point is the one with the lowest projected
    !> derivative.
    logical,  intent(out) :: tConverged

    real(dp) :: phi, dphi, c_check, q_check, dalpha

    @:ASSERT(size(xNew) == self%nElem)
    @:ASSERT(size(dx) == self%nElem)

    tConverged = .false.

    phi = fx
    dphi = dot_product(dx, self%d0)

    if (.not. self%isZoom) then
      self%iter = self%iter + 1

      if (self%iter == 1) then
        self%phi_0 = fx
        self%dx_0 = dx
        self%dphi_0 = dphi

        self%phi_lo = fx
        self%dx_lo = dx
        self%dphi_lo = dphi

        if (self%alpha_new > self%maxAlpha) then
          self%alpha_new = self%maxAlpha
        end if

        self%xNew = self%x0 + self%alpha_new * self%d0
        xNew = self%xNew
        return
      end if

      if ((phi > self%phi_0 + WOLFE1 * self%alpha_new * self%dphi_0)&
          & .or.  ( (phi >= self%phi_lo) .and. (self%iter > 2) )) then
        self%alpha_hi = self%alpha_new
        self%phi_hi = phi
        self%dphi_hi = dphi
        self%dx_hi = dx

        self%isZoom = .true.
        self%iter = 0

      else if (abs(dphi) <= -WOLFE2*self%dphi_0) then
        tConverged = .true.
        self%alpha_lo = self%alpha_new
        self%phi_lo = phi
        self%dphi_lo = dphi
        self%dx_lo = dx
        xNew = self%xNew
        return

      else if (dphi >= 0.0_dp) then
        self%alpha_hi = self%alpha_lo
        self%phi_hi = self%phi_lo
        self%dphi_hi = self%dphi_lo
        self%dx_hi = self%dx_lo

        self%alpha_lo = self%alpha_new
        self%phi_lo = phi
        self%dphi_lo = dphi
        self%dx_lo = dx

        self%isZoom = .true.
        self%iter = 0
        self%xNew = self%x0 + self%alpha_new * self%d0

      else

        self%alpha_lo = self%alpha_new
        self%phi_lo = phi
        self%dphi_lo = dphi
        self%dx_lo = dx

        if (abs(self%alpha_new - self%maxAlpha) < epsilon(1.0_dp)) then
          call warning("Line Search hit maximum border.")
          tConverged = .true.
          xNew = self%xNew
          return
        else if ((self%alpha_new > self%maxAlpha) .and. self%iter > 1) then
          xNew = self%x0 + self%maxAlpha * self%d0
          return
        end if

        self%alpha_new = 2.0_dp*self%alpha_new

        self%xNew = self%x0 + self%alpha_new * self%d0
        xNew = self%xNew
        return
      end if
    end if

    ! Zoom Phase

    if (self%iter > self%mIter) then
      if (self%alpha_lo*maxval(self%d0) < self%tol) then
        if (self%alpha_new*maxval(self%d0) < self%tol) then
          self%alpha_new = 0.5_dp*self%maxAlpha
          self%xNew = self%x0 + self%alpha_new*self%d0
          xNew = self%xNew
          return
        else
          xNew = self%xNew
          self%alpha_lo = self%alpha_new
          self%phi_lo = phi
          self%dphi_lo = dphi
          self%dx_lo = dx
        end if
      end if
      xNew = self%x0 + self%alpha_lo * self%d0

      call warning("LineSearch did not Converged. (zoom)")
      tConverged = .true.
      return
    end if

    ! Check new alpha
    if (self%iter > 0)then
      ! Skip in First Zoom run
      if ((phi > self%phi_0 + WOLFE1*self%alpha_new*self%dphi_0) .or. (phi >= self%phi_lo))&
          & then
        self%phi_temp = self%phi_hi
        self%alpha_temp = self%alpha_hi
        self%alpha_hi = self%alpha_new
        self%phi_hi = phi
        self%dphi_hi = dphi
        self%dx_hi = dx
      else
        if (abs(dphi) <= -WOLFE2*self%dphi_0) then
          tConverged = .true.
          self%alpha_lo = self%alpha_new
          self%phi_lo = phi
          self%dphi_lo = dphi
          self%dx_lo = dx
          xNew = self%xNew
          return
        else if (dphi*(self%alpha_hi - self%alpha_lo) >= 0.0_dp) then
          self%phi_temp = self%phi_hi
          self%alpha_temp = self%alpha_hi
          self%alpha_hi = self%alpha_lo
          self%phi_hi = self%phi_lo
          self%dphi_hi = self%dphi_lo
          self%dx_hi = self%dx_lo
        else
          self%phi_temp = self%phi_lo
          self%alpha_temp = self%alpha_lo
        end if
        self%alpha_lo = self%alpha_new
        self%phi_lo = phi
        self%dphi_lo = dphi
        self%dx_lo = dx
      end if

    end if

    ! Zoom phase
    dalpha = self%alpha_hi - self%alpha_lo

    if (self%iter > 0) then
      c_check = delta1 * dalpha
      self%alpha_new = cubicmin(self%alpha_lo, self%phi_lo, self%dphi_lo, self%alpha_hi,&
          & self%phi_hi, self%alpha_temp, self%phi_temp)
    end if

    if (ieee_is_nan(self%alpha_new) .or. (self%iter == 0) .or.&
        & (self%alpha_new > self%alpha_hi - c_check)&
        & .or. (self%alpha_new < self%alpha_lo + c_check) ) then
      q_check = delta2 * dalpha
      self%alpha_new = quadmin(self%alpha_lo, self%phi_lo, self%dphi_lo,&
          & self%alpha_hi, self%phi_hi)
      if (ieee_is_nan(self%alpha_new) .or. (self%alpha_new > self%alpha_hi - q_check)&
          & .or. (self%alpha_new < self%alpha_lo + q_check)) then
        ! Try Bisection
        self%alpha_new = self%alpha_lo + 0.5_dp*dalpha
      end if
    end if

    self%iter = self%iter + 1

    if (abs(self%alpha_new - self%alpha_lo) < self%tol) then
      self%alpha_new = self%alpha_lo
      tConverged = .true.
    end if

    self%xNew = self%x0 + self%alpha_new*self%d0
    xNew = self%xNew

  end subroutine LineMin_next

  !> Finds the minimizer for a cubic polynomial that goes through the
  function cubicmin(a, fa, fpa, b, fb, c, fc) result(xmin)

    !> Point a
    real(dp), intent(in) :: a

    !> function value at a
    real(dp), intent(in) :: fa

    !> Derviative at a
    real(dp), intent(in) :: fpa

    !> point b
    real(dp), intent(in) :: b

    !> value at b
    real(dp), intent(in) :: fb

    !> point c
    real(dp), intent(in) :: c

    !> value at c
    real(dp), intent(in) :: fc

    !> minimum for the quadratic
    real(dp) :: xmin

    real(dp) :: db, dc, denom, d1(2,2), radical, ta, tb, tc, temp(2)

    tc = fpa
    db = b - a
    dc = c - a
    denom = (db * dc)**2 * (db - dc)
    d1 = reshape([ dc**2, -dc**3 , db**2, db**3], [2,2])
    temp = matmul(d1, [fb-fa-tc*db, fc-fa-tc*dc])

    ta = temp(1) / denom
    tb = temp(2) / denom
    radical = tb*tb - 3*ta*tc
    xmin = a + (-tb + sqrt(radical)) / (3*ta)

  end function cubicmin

  !> Finds the minimizer for a quadratic polynomial that goes through the
  function quadmin(a, fa, fpa, b, fb) result(xmin)

    !> Point a
    real(dp), intent(in) :: a

    !> function value at a
    real(dp), intent(in) :: fa

    !> Derviative at a
    real(dp), intent(in) :: fpa

    !> point b
    real(dp), intent(in) :: b

    !> value at b
    real(dp), intent(in) :: fb

    !> minimum for the quadratic
    real(dp) :: xmin

    real(dp) :: td, tc, db, tb

    td = fa
    tc = fpa
    db = b - a
    tb = (fb - td - tc*db) / (db*db)
    xmin = a - tc / (2.0_dp*tb)

  end function quadmin

  !> Gives the coordinate of the minimal point back
  !> The value passed back is meaningless if the subroutine is called before the line minimizer
  !> signals convergence.
  subroutine LineMin_getMinX(self, minX)

    !> Line minimizer
    type(OLineSearch), intent(in) :: self

    !> Coordinate of the minimal point
    real(dp), intent(out) :: minX(:)

    minX(:) = self%x0(:) + self%alpha_lo * self%d0

  end subroutine LineMin_getMinX


  !> Returns the function at the minimal point
  !> The value passed back is meaningless if the subroutine is called before the line minimizer
  !> signals convergence.
  subroutine LineMin_getMinY(self, minY)

    !> Line minimizer
    type(OLineSearch), intent(in) :: self

    !> Function value in the minimal point
    real(dp), intent(out) :: minY

    minY = self%phi_lo

  end subroutine LineMin_getMinY


  !> Gives the gradient in the minimal point back
  !> The value passed back is meaningless if the subroutine is called before the line minimizer
  !> signals convergence.
  subroutine LineMin_getMinGrad(self, minGrad)

    !> Line minimizer
    type(OLineSearch), intent(in) :: self

    !> Gradient in the minimal point
    real(dp), intent(out) :: minGrad(:)

    minGrad(:) = self%dx_lo

  end subroutine LineMin_getMinGrad

  !> Returns the displacement to the minimum along the line
  !> The value passed back is meaningless if the subroutine is called before the line minimizer
  !> signals convergence.
  subroutine LineMin_getMinLambda(self, minLambda)

    !> Line minimizer
    type(OLineSearch), intent(in) :: self

    !> Displacement along the line to the minimum
    real(dp), intent(out) :: minLambda

    minLambda = self%alpha_lo

  end subroutine LineMin_getMinLambda

end module lbfgs
