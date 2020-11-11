!--------------------------------------------------------------------------------------------------!
!  DFTB+: general package for performing fast atomistic simulations                                !
!  Copyright (C) 2006 - 2020  DFTB+ developers group                                               !
!                                                                                                  !
!  See the LICENSE file for terms of usage and distribution.                                       !
!--------------------------------------------------------------------------------------------------!

#:include 'common.fypp'

!> Function minimization with lbfgs algorithm
!>
!> References:
!> Nocedal - Updating Quasi-Newton Matrices with Limited Storage (1980), Mathematics of Computation
!> 35, pp. 773-782.
!>
!> Nocedal, Wright - Numerical Optimization, Springer
!>
module dftbp_lbfgs
  use, intrinsic :: ieee_arithmetic
  use dftbp_accuracy
  use dftbp_assert
  use dftbp_message
  use dftbp_linemin, only : TLineMin, TLineMin_init
  implicit none
  private

  public :: TLbfgs, TLbfgs_init


  !> Holds data for the original line minimalizer used in this implementation
  type TLineSearch
    private

    !> Elements of vector x
    integer :: nElem

    !> Number of current iteration
    integer :: iter

    !> Maximum number of iterations
    integer :: mIter

    !> Search direction
    real(dp), allocatable :: d0(:)

    !> initial point
    real(dp), allocatable :: x0(:)

    !> New point
    real(dp), allocatable :: xNew(:)

    !> Gradient at initial point
    real(dp), allocatable :: dx0(:)

    !> Gradient of fx at low value
    real(dp), allocatable :: dxLo(:)

    !> Gradient of fx at high value
    real(dp), allocatable :: dxHi(:)

    !> new line search position as fraction of gradient
    real(dp) :: alphaNew

    !> line search position on the high side
    real(dp) :: alphaHi

    !>> line search position on the low side of the gradient
    real(dp) :: alphaLo

    !> phi(alpha) = f(x + alpha*d0)
    real(dp) :: phi0

    !> Function on the high side
    real(dp) :: phiHi

    !> Function on the low side
    real(dp) :: phiLo

    !> Gradient of phi at initial point
    real(dp) :: dphi0

    !> Gradient of phi on the high side
    real(dp) :: dphiHi

    !> gradient of phi on the low side
    real(dp) :: dphiLo

    !> minimum scale of step along gradient
    real(dp) :: minAlpha

    !> maximum scale of step along gradient
    real(dp) :: maxAlpha

    !> Minimal displacement in at least one component in one step
    real(dp) :: minDisp

    !> Maximal displacement in one step
    real(dp) :: maxDisp

    !> True if program is in Zoom state
    logical :: tZoom

    !> True if converged
    logical :: tConverged

    !> Holding variable for phi between steps
    real(dp) :: phiTemp

    !> holding variable for alpha between steps
    real(dp) :: alphaTemp

  contains

    !> Resets lbfgs minimizer
    procedure :: reset => TLineSearch_reset

    !> Passes the function value and the derivative of the last point to line minimizer and gives a
    !> new coordinate back.
    procedure :: next => TLineSearch_next

    !> Returns the function at the minimal point
    procedure :: getMinY => TLineSearch_getMinY

    !> Gives the gradient in the minimal point back.
    procedure :: getMinGrad => TLineSearch_getMinGrad

    !> Returns the displacement to the minimum along the line.
    procedure :: getMinLambda => TLineSearch_getMinLambda

    !> Gives the coordinate of the minimal point back.
    procedure :: getMinX => TLineSearch_getMinX

  end type TLineSearch


  !> Calculates search direction via lbfgs-Algorithm
  type :: TLbfgs
    private

    !> Is a line search used
    logical :: isLineSearch

    !> Is the original line search used
    logical :: isOldLSUsed

    !> line minimizer
    type(TLineMin), allocatable :: lineMin

    !> original internal line search
    type(TLineSearch), allocatable :: lineSearch

    !> Number of elements
    integer :: nElem

    !> current point
    real(dp), allocatable :: xx(:)

    !> line search position
    real(dp) :: alpha

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

    !> Maximum step size to take for the variables
    real(dp) :: maxDisp

    !> tolerance for gradient
    real(dp) :: tol

  contains

    !> Next structure from LBFGS
    procedure :: next => TLbfgs_next

    !> reset the search
    procedure :: reset => TLbfgs_reset

  end type TLbfgs


  !> Lower Wolfe condition parameter of LBFGS
  real(dp), parameter :: wolfe1 = 1e-4_dp

  !> Upper Wolfe condition parameter for LBFGS
  real(dp), parameter :: wolfe2 = 0.9_dp


contains

  !> Initialize lbfgs instance
  subroutine TLbfgs_init(this, nElem, tol, minDisp, maxDisp, mem, isLineSearch, isOldLSUsed,&
      & isQNDisp)

    !> lbfgs instance on exit
    type(TLbfgs), intent(out) :: this

    !> Dimension of vectors
    integer, intent(in) :: nElem

    !> Tolerance for gradient
    real(dp), intent(in) :: tol

    !> Minimal displacement in at leat one component in one step
    real(dp), intent(in) :: minDisp

    !> Maximal displacement in one step
    real(dp), intent(in) :: maxDisp

    !> Number of past iterations which will be saved
    integer, intent(in) :: mem

    !> Is a line search used along the quasi-Newton direction
    logical, intent(in) :: isLineSearch

    !> Is the old line search routine used, instead of the line minimizer
    logical, intent(in) :: isOldLSUsed

    !> Is the maximum step size considered for the QN
    logical, intent(in) :: isQNDisp

    @:ASSERT(nElem > 0)
    @:ASSERT(tol > 0.0_dp)
    @:ASSERT(maxDisp > 0.0_dp)
    @:ASSERT(mem > 0)

    this%nElem = nElem
    this%mem = mem
    this%tol = tol
    allocate(this%xx(nElem))
    allocate(this%xOld(nElem))
    allocate(this%gg(nElem))
    allocate(this%ss(nElem, mem))
    allocate(this%yy(nElem, mem))
    allocate(this%rho(mem))
    allocate(this%dir(nElem))

    this%isLineSearch = isLineSearch

    this%maxDisp = -1.0_dp
    if (this%isLineSearch) then
      this%isOldLSUsed = isOldLSUsed
      if (this%isOldLSUsed) then
        allocate(this%lineSearch)
        call TLineSearch_init(this%lineSearch, nElem, 10, minDisp, maxDisp)
      else
        allocate(this%lineMin)
        call TLineMin_init(this%lineMin, nElem, 10, tol, maxDisp)
      end if
    else
      if (isQNDisp) then
        this%maxDisp = maxDisp
      end if
    end if

  end subroutine TLbfgs_init


  !> Passes calculated function value and gradient to the minimizer and
  !> gives a new coordinate back.
  subroutine TLbfgs_next(this, fx, dx, xNew, tConverged)

    !> LBFGS Instance.
    class(TLbfgs), intent(inout) :: this

    !> Function value in the last point
    real(dp), intent(in)  :: fx

    !> Gradient in the last point
    real(dp), intent(in)  :: dx(:)

    !> New proposed point
    real(dp), intent(out) :: xNew(:)

    !> True, if gradient got below the specified tolerance.
    logical,  intent(out) :: tConverged

    real(dp) :: dxTemp(this%nElem), dxMax
    logical :: tLineConverged

    @:ASSERT(size(xNew) == this%nElem)
    @:ASSERT(size(dx) == this%nElem)

    if (maxval(abs(dx)) < this%tol) then
      tConverged = .true.
      xNew(:) = this%xx
      return
    else
      tConverged = .false.
    end if
    dxTemp(:) = dx

    if (this%isLineSearch) then
      if (this%iter > 0) then
        if (this%isOldLSUsed) then
          call this%lineSearch%next(fx, dx, this%xx, tLineConverged)
        else
          call this%lineMin%next(fx, dx, this%xx, tLineConverged)
        end if
        if (tLineConverged) then
          if (this%isOldLSUsed) then
            call this%lineSearch%getMinGrad(dxTemp)
          else
            call this%lineMin%getMinGrad(dxTemp)
          end if
        else
          xNew(:) = this%xx
          return
        end if
      end if
    end if

    ! Calculate new search direction
    call TLbfgs_calcDirection(this, this%xx, dxTemp)

    if (this%isLineSearch) then
      this%alpha = 1.0_dp
      if (this%isOldLSUsed) then
        call this%lineSearch%reset(this%xx, this%dir, this%alpha)
        call this%lineSearch%next(fx, dx, this%xx, tLineConverged)
      else
        call this%lineMin%reset(this%xx, this%dir, this%alpha)
        call this%lineMin%next(fx, dx, this%xx, tLineConverged)
      end if
    else
      if (this%maxDisp > 0.0_dp) then
        dxMax = maxval(abs(this%dir))
        if (dxMax > this%maxDisp) then
          this%dir(:) = this%dir * this%maxDisp / dxMax
        end if
      end if
      this%xx(:) = this%xx + this%dir
    end if

    xNew(:) = this%xx

  end subroutine TLbfgs_next


  !> Resets lbfgs minimizer
  subroutine TLbfgs_reset(this, xx)

    !> lbfgs minimizer
    class(TLbfgs), intent(inout) :: this

    !> Start point
    real(dp), intent(in) :: xx(:)

    @:ASSERT(size(xx) == this%nElem)

    this%alpha = 1.0_dp
    this%xx(:) = xx
    this%xOld(:) = xx
    this%gg(:) = 0.0_dp
    this%ss(:,:) = 0.0_dp
    this%yy(:,:) = 0.0_dp
    this%rho(:) = 0.0_dp
    this%dir(:) = 0.0_dp
    this%iter = 0

  end subroutine TLbfgs_reset


  !> Calculates next search direction
  subroutine TLbfgs_calcDirection(this, xx, gg, diagIn)

    !> lbfgs instance
    class(TLbfgs), intent(inout) :: this

    !> Point
    real(dp), intent(in) :: xx(:)

    !> Gradient at point
    real(dp), intent(in) :: gg(:)

    !> Initial diagonal hessian matrix (default identity)
    real(dp), intent(in), optional :: diagIn(:)

    ! Element where new data of this iteration will be stored
    integer :: newElem

    ! Oldest Element
    integer :: oldestElem

    ! Initial hesse matrix
    real(dp) :: diag(this%nElem)

    ! Workspace variables
    real(dp) :: qq(this%nElem), alpha(this%mem), beta

    integer :: ii, mm

    ! First Run
    if (this%iter == 0) then
      this%gg(:) = gg
      this%dir(:) = -gg
      this%iter = this%iter + 1
      return
    end if

    ! Checks whether new x is different to last iteration
    if (maxval(abs(xx - this%xOld)) < epsilon(1.0_dp)) then
      call warning("Error: x need to be different in each iteration. (LBFGS Minimizer)")
    end if

    ! Checks whether new gradient is different to last iteration
    if (maxval(abs(gg - this%gg)) < epsilon(1.0_dp)) then
      call warning("Error: g need to be different in each iteration. (LBFGS Minimizer)")
    end if

    ! Save new ss, yy and rho
    newElem = modulo(this%iter - 1, this%mem) + 1
    this%ss(:,newElem) = xx - this%xOld
    this%yy(:,newElem) = gg - this%gg

    beta = dot_product(this%yy(:,newElem), this%ss(:,newElem))

    ! y*s needs to be positive to get a positive definite hermitian matrix
    !if (beta <= epsilon(1.0_dp) ) then
    !  call warning("Bad step. y*s <= 0 (LBFGS Minimizer)")
    !end if
    this%rho(newElem) = 1.0_dp / max(abs(beta), epsilon(1.0_dp))

    ! Save new x and g
    this%xOld(:) = xx
    this%gg(:) = gg

    if (.not. present(diagIn)) then
      ! Compute guessed hessian matrix if not set by the user
      if (this%iter == 1) then
        diag(:) = 1.0_dp
      else
        diag(:) = dot_product(this%ss(:, newElem), this%yy(:, newElem))&
            & / max(dot_product(this%yy(:, newElem), this%yy(:, newElem)), epsilon(1.0_dp))
      end if

    else if (any(diagIn < 0.0_dp)) then
      ! If hessian set by the user, do input check
      call error("Error: Bad Input (diag). (LBFGS Minimizer)")
    else
      diag(:) = diagIn
    end if

    ! Compute d = -H*g using the formula given by Nocedal. Result saved in this%dir
    qq(:) = gg
    do mm = 0, this%mem - 1
      ii = modulo(newElem - mm - 1, this%mem) + 1
      alpha(ii) = this%rho(ii) * dot_product(this%ss(:,ii), qq)
      qq(:) = qq - alpha(ii) * this%yy(:,ii)
    end do

    this%dir(:) = diag * qq

    oldestElem = modulo(newElem, this%mem) + 1
    do mm = 0, this%mem - 1
      ii = modulo(oldestElem + mm - 1, this%mem) + 1
      beta = this%rho(ii) * dot_product(this%yy(:,ii), this%dir)
      this%dir(:) = this%dir + this%ss(:,ii) * (alpha(ii) - beta)
    end do

    this%iter = this%iter + 1

    ! Change to correct sign
    this%dir(:) = -this%dir

  end subroutine TLbfgs_calcDirection


! Internal linesearch routines

  !> Creates a new line minimizer
  subroutine TLineSearch_init(this, nElem, mIter, minDisp, maxDisp)

    !> Valid line minimizer instance on exit
    type(TLineSearch), intent(out) :: this

    !> Nr. of elements in the coordinate/gradient vectors
    integer, intent(in) :: nElem

    !> Nr. of maximal iterations to perform (>3)
    integer, intent(in) :: mIter

    !> Minimal displacement in at leat one component in one step
    real(dp), intent(in) :: minDisp

    !> Maximal movement in one coordinate in one step
    real(dp), intent(in) :: maxDisp

    @:ASSERT(nElem > 0)
    @:ASSERT(mIter > 3)
    @:ASSERT(maxDisp > 0.0_dp)

    this%nElem = nElem
    this%minDisp = minDisp
    this%maxDisp = maxDisp
    this%mIter = mIter

    allocate(this%d0(nElem))
    allocate(this%dxLo(nElem))
    allocate(this%dxHi(nElem))
    allocate(this%xNew(nElem))
    allocate(this%x0(nElem))
    allocate(this%dx0(nElem))

  end subroutine TLineSearch_init


  !> Resets the line minimizer
  subroutine TLineSearch_reset(this, x0, d0, firstStep)

    !> Line minimizer instance
    class(TLineSearch), intent(inout) :: this

    !> New starting point
    real(dp), intent(in) :: x0(:)

    !> New direction
    real(dp), intent(in) :: d0(:)

    !> First step size to take
    real(dp), intent(in) :: firstStep

    @:ASSERT(size(x0) == this%nElem)
    @:ASSERT(size(d0) == this%nElem)
    #!@:ASSERT(all(abs(d0) > 0.0_dp))

    this%d0(:) = d0
    this%phi0 = 0.0_dp
    this%dphi0 = 0.0_dp

    this%iter = 0
    this%tZoom = .false.
    this%tConverged = .false.

    this%alphaLo = 0.0_dp
    this%alphaHi = 0.0_dp

    this%phiLo = 0.0_dp
    this%phiHi = 0.0_dp
    this%dphiLo = 0.0_dp
    this%dphiHi = 0.0_dp

    this%dx0(:) = 0.0_dp
    this%dxLo(:) = 0.0_dp
    this%dxHi(:) = 0.0_dp

    this%x0(:) = x0
    this%xNew(:) = 0.0_dp

    this%alphaNew = firstStep
    this%minAlpha = this%minDisp / maxval(abs(d0))
    this%maxAlpha = this%maxDisp / maxval(abs(d0))

    this%alphaTemp = 0.0_dp
    this%phiTemp = 0.0_dp

  end subroutine TLineSearch_reset


  !> Passes the function value and the derivative of the last point to line minimizer and gives a
  !> new coordinate back.
  !>
  !> Getting back tConverged = .true. can also mean that the line minimization did not converge in
  !> the maximal nr. of steps.
  !> When calling this subroutine the first time, function value and gradient for the starting point
  !> of the minimization should be passed.
  !>
  subroutine TLineSearch_next(this, fx, dx, xNew, tConverged)

    !> Line minimizer instance
    class(TLineSearch), intent(inout) :: this

    !> Function value for the last returned point
    real(dp), intent(in) :: fx

    !> Gradient for the last returned point
    real(dp), intent(in) :: dx(:)

    !> New point to calculate
    real(dp), intent(out) :: xNew(:)

    !> True, if line minimization converged. The last passed point is the one with the lowest
    !> projected derivative.
    logical,  intent(out) :: tConverged

    real(dp) :: phi, dphi

    @:ASSERT(size(xNew) == this%nElem)
    @:ASSERT(size(dx) == this%nElem)

    tConverged = .false.

    phi = fx
    dphi = dot_product(dx, this%d0)

    if (.not. this%tZoom) then
      call TLineSearch_getOptimalStep(this, phi, dphi, dx, xNew, tConverged)
    end if
    if (this%tZoom) then
      call TLineSearch_zoom(this, phi, dphi, dx, xNew, tConverged)
    end if

    this%iter = this%iter + 1

  end subroutine TLineSearch_next


  !> Returns the function at the minimal point
  !>
  !> The value passed back is meaningless if the subroutine is called before the line minimizer
  !> signals convergence.
  !>
  subroutine TLineSearch_getMinY(this, minY)

    !> Line minimizer
    class(TLineSearch), intent(in) :: this

    !> Function value in the minimal point
    real(dp), intent(out) :: minY

    minY = this%phiLo

  end subroutine TLineSearch_getMinY


  !> Gives the gradient in the minimal point back.
  !>
  !> The value passed back is meaningless if the subroutine is called before the line minimizer
  !> signals convergence.
  !>
  subroutine TLineSearch_getMinGrad(this, minGrad)

    !> Line minimizer
    class(TLineSearch), intent(in) :: this

    !> Gradient in the minimal point
    real(dp), intent(out) :: minGrad(:)

    minGrad(:) = this%dxLo

  end subroutine TLineSearch_getMinGrad


  !> Returns the displacement to the minimum along the line.
  !>
  !> The value passed back is meaningless if the subroutine is called before the line minimizer
  !> signals convergence.
  !>
  subroutine TLineSearch_getMinLambda(this, minLambda)

    !> Line minimizer
    class(TLineSearch), intent(in) :: this

    !> Displacement along the line to the minimum
    real(dp), intent(out) :: minLambda

    minLambda = this%alphaLo

  end subroutine TLineSearch_getMinLambda


  !> Gives the coordinate of the minimal point back.
  !>
  !> The value passed back is meaningless if the subroutine is called before the line minimizer
  !> signals convergence.
  !>
  subroutine TLineSearch_getMinX(this, minX)

    !> Line minimizer
    class(TLineSearch), intent(in) :: this

    !> Coordinate of the minimal point
    real(dp), intent(out) :: minX(:)

    minX(:) = this%x0(:) + this%alphaLo * this%d0

  end subroutine TLineSearch_getMinX


  !> Tries to find the optimal or a bracket containing the optimal step.
  subroutine TLineSearch_getOptimalStep(this, phi, dphi, dx, xNew, tConverged)

    !> Instance.
    type(TLineSearch), intent(inout) :: this

    !> Function value at current point
    real(dp), intent(in) :: phi

    !> Derivative along search direction at current point
    real(dp), intent(in) :: dphi

    !> Derivative at current point
    real(dp), intent(in) :: dx(:)

    !> Next point to calculate
    real(dp), intent(out) :: xNew(:)

    !> Whether line minimisation converged
    logical, intent(out) :: tConverged

    if (this%iter == 0) then
      this%phi0 = phi
      this%dx0(:) = dx
      this%dphi0 = dphi

      this%phiLo = phi
      this%dxLo(:) = dx
      this%dphiLo = dphi

      this%alphaNew = min(this%alphaNew, this%maxAlpha)
      this%xNew(:) = this%x0 + this%alphaNew * this%d0
      xNew(:) = this%xNew
      return
    end if

    if ((phi > this%phi0 + wolfe1 * this%alphaNew * this%dphi0)&
        & .or.  ((phi >= this%phiLo) .and. (this%iter > 1))) then
      this%alphaHi = this%alphaNew
      this%phiHi = phi
      this%dphiHi = dphi
      this%dxHi(:) = dx

      this%tZoom = .true.
      this%iter = 0

    else if (abs(dphi) <= -wolfe2 * this%dphi0) then
      tConverged = .true.
      this%alphaLo = this%alphaNew
      this%phiLo = phi
      this%dphiLo = dphi
      this%dxLo = dx
      xNew(:) = this%xNew

    else if (dphi >= 0.0_dp) then
      this%alphaHi = this%alphaLo
      this%phiHi = this%phiLo
      this%dphiHi = this%dphiLo
      this%dxHi(:) = this%dxLo

      this%alphaLo = this%alphaNew
      this%phiLo = phi
      this%dphiLo = dphi
      this%dxLo(:) = dx

      this%tZoom = .true.
      this%iter = 0
      this%xNew(:) = this%x0 + this%alphaNew * this%d0

    else
      this%alphaLo = this%alphaNew
      this%phiLo = phi
      this%dphiLo = dphi
      this%dxLo(:) = dx

      if (abs(this%alphaNew - this%maxAlpha) < epsilon(1.0_dp)) then
        call warning("Line Search hit maximum border.")
        tConverged = .true.
        xNew(:) = this%xNew
      else if ((this%alphaNew > this%maxAlpha) .and. this%iter > 0) then
        xNew(:) = this%x0 + this%maxAlpha * this%d0
      else
        this%alphaNew = 2.0_dp * this%alphaNew
        this%xNew(:) = this%x0 + this%alphaNew * this%d0
        xNew(:) = this%xNew
      end if

    end if

  end subroutine TLineSearch_getOptimalStep


  !> Finds the optimal step by zooming into the brackets.
  subroutine TLineSearch_zoom(this, phi, dphi, dx, xNew, tConverged)

    !> Instance.
    type(TLineSearch), intent(inout) :: this

    !> Function value at current point
    real(dp), intent(in) :: phi

    !> Derivative along search direction at current point
    real(dp), intent(in) :: dphi

    !> Derivative at current point
    real(dp), intent(in) :: dx(:)

    !> Next point to calculate
    real(dp), intent(out) :: xNew(:)

    !> Whether line minimisation converged
    logical, intent(out) :: tConverged

    !> Check for cubic interpolation
    real(dp), parameter :: delta1 = 0.2_dp

    !> Check for quadratic interpolation
    real(dp), parameter :: delta2 = 0.1_dp

    real(dp) :: cCheck, qCheck, dAlpha

    if (this%iter > this%mIter) then
      if (abs(this%alphaLo * maxval(this%d0)) < this%minAlpha) then
        if (abs(this%alphaNew * maxval(this%d0)) < this%minAlpha) then
          this%alphaNew = 0.5_dp * this%maxAlpha
          this%xNew(:) = this%x0 + this%alphaNew * this%d0
          xNew(:) = this%xNew
          return
        else
          xNew(:) = this%xNew
          this%alphaLo = this%alphaNew
          this%phiLo = phi
          this%dphiLo = dphi
          this%dxLo(:) = dx
        end if
      end if
      xNew(:) = this%x0 + this%alphaLo * this%d0

      call warning("LineSearch did not Converged. (zoom)")
      tConverged = .true.
      return
    end if

    ! Check new alpha
    if (this%iter > 0)then
      ! Skip in First Zoom run
      if ((phi > this%phi0 + wolfe1 * this%alphaNew * this%dphi0) .or. (phi >= this%phiLo)) then
        this%phiTemp = this%phiHi
        this%alphaTemp = this%alphaHi
        this%alphaHi = this%alphaNew
        this%phiHi = phi
        this%dphiHi = dphi
        this%dxHi(:) = dx
      else
        if (abs(dphi) <= -wolfe2 * this%dphi0) then
          tConverged = .true.
          this%alphaLo = this%alphaNew
          this%phiLo = phi
          this%dphiLo = dphi
          this%dxLo(:) = dx
          xNew(:) = this%xNew
          return
        else if (dphi * (this%alphaHi - this%alphaLo) >= 0.0_dp) then
          this%phiTemp = this%phiHi
          this%alphaTemp = this%alphaHi
          this%alphaHi = this%alphaLo
          this%phiHi = this%phiLo
          this%dphiHi = this%dphiLo
          this%dxHi(:) = this%dxLo
        else
          this%phiTemp = this%phiLo
          this%alphaTemp = this%alphaLo
        end if
        this%alphaLo = this%alphaNew
        this%phiLo = phi
        this%dphiLo = dphi
        this%dxLo(:) = dx
      end if
    end if

    ! Zoom phase
    dAlpha = this%alphaHi - this%alphaLo

    if (this%iter > 0) then
      cCheck = delta1 * dAlpha
      this%alphaNew = cubicMin(this%alphaLo, this%phiLo, this%dphiLo, this%alphaHi, this%phiHi,&
          & this%alphaTemp, this%phiTemp)
    end if

    if (ieee_is_nan(this%alphaNew) .or. (this%iter == 0)&
        & .or. (this%alphaNew > this%alphaHi - cCheck)&
        & .or. (this%alphaNew < this%alphaLo + cCheck)) then
      qCheck = delta2 * dAlpha
      this%alphaNew = quadMin(this%alphaLo, this%phiLo, this%dphiLo, this%alphaHi, this%phiHi)
      if (ieee_is_nan(this%alphaNew) .or. (this%alphaNew > this%alphaHi - qCheck)&
          & .or. (this%alphaNew < this%alphaLo + qCheck)) then
        ! Try Bisection
        this%alphaNew = this%alphaLo + 0.5_dp * dAlpha
      end if
    end if

    if (abs(this%alphaNew - this%alphaLo) < this%minAlpha) then
      this%alphaNew = this%alphaLo
      tConverged = .true.
    end if

    this%xNew = this%x0 + this%alphaNew * this%d0
    xNew = this%xNew

  end subroutine TLineSearch_zoom


  !> Finds the minimizer for a cubic polynomial that goes through three points and has a given
  !> derivative in the first point.
  function cubicMin(aa, fa, fpa, bb, fb, cc, fc) result(xmin)

    !> Point a
    real(dp), intent(in) :: aa

    !> function value at a
    real(dp), intent(in) :: fa

    !> Derviative at a
    real(dp), intent(in) :: fpa

    !> point b
    real(dp), intent(in) :: bb

    !> value at b
    real(dp), intent(in) :: fb

    !> point c
    real(dp), intent(in) :: cc

    !> value at c
    real(dp), intent(in) :: fc

    !> minimum for the cubic
    real(dp) :: xMin

    real(dp) :: db, dc, denom, d1(2,2), radical, ta, tb, temp(2)

    db = bb - aa
    dc = cc - aa
    denom = (db * dc)**2 * (db - dc)
    d1(:,:) = reshape([ dc**2, -dc**3, -db**2, db**3], [2, 2])
    temp(:) = matmul(d1, [fb - fa - fpa * db, fc - fa - fpa * dc])

    ta = temp(1) / denom
    tb = temp(2) / denom
    radical = tb * tb - 3 * ta * fpa
    xMin = aa + (-tb + sqrt(radical)) / (3 * ta)

  end function cubicMin


  !> Finds the minimizer for a quadratic polynomial that goes through two points and has a given
  !> derivative in the first point.
  function quadMin(aa, fa, fpa, bb, fb) result(xMin)

    !> Point a
    real(dp), intent(in) :: aa

    !> function value at a
    real(dp), intent(in) :: fa

    !> Derviative at a
    real(dp), intent(in) :: fpa

    !> point b
    real(dp), intent(in) :: bb

    !> value at b
    real(dp), intent(in) :: fb

    !> minimum for the quadratic
    real(dp) :: xmin

    real(dp) :: db

    db = bb - aa
    xMin = aa - fpa * db * db / (2.0_dp * (fb - fa - fpa * db))

  end function quadMin


end module dftbp_lbfgs
