!--------------------------------------------------------------------------------------------------!
!  DFTB+: general package for performing fast atomistic simulations                                !
!  Copyright (C) 2006 - 2020  DFTB+ developers group                                               !
!                                                                                                  !
!  See the LICENSE file for terms of usage and distribution.                                       !
!--------------------------------------------------------------------------------------------------!

#:include 'common.fypp'

!> Function minimization with the standard conjugate gradient technique.  See Numerical Recipes for
!> details.
module dftbp_conjgrad
  use dftbp_assert
  use dftbp_accuracy
  use dftbp_linemin
  implicit none

  private


  !> contains data for the conjugate gradient minimizer
  type TConjGrad
    private

    !> State of the minimizer
    integer :: state

    !> Nr. of variables
    integer :: nElem

    !> Gradient in previous cycle
    real(dp), allocatable :: gg(:)

    !> conjugate gradient
    real(dp), allocatable :: hh(:)

    !> Last calculated point
    real(dp), allocatable :: uu(:)

    !> Tolerance criteria for convergence
    real(dp) :: tolerance

    !> Maximal displacement along one coordinate in one step
    real(dp) :: maxDisp

    !> If CG converged
    logical :: tConverged

    !> If object is initialized
    logical :: tInitialized

    !> Line minimizer
    type(TLineMin) :: pLinMin

  end type TConjGrad


  !> Initialises CG instance
  interface init
    module procedure conjGrad_init
  end interface init


  !> Resets CG
  interface reset
    module procedure conjGrad_reset
  end interface reset


  !> Passes calculated function value and gradient to the minimizer and gives
  !> back a new coordinate.
  interface next
    module procedure conjGrad_next
  end interface next


  !> Coordinate of the minimum.
  interface getMinX
    module procedure conjGrad_getMinX
  end interface getMinX


  !> Function value in the minimum.
  interface getMinY
    module procedure conjGrad_getMinY
  end interface getMinY


  !> Gradient in the minimum.
  interface getMinGrad
    module procedure conjGrad_getMinGrad
  end interface getMinGrad

  public :: TConjGrad
  public :: init, reset, next, getMinX, getMinY, getMinGrad


  !> State of the conjugate gradient cycle
  integer, parameter :: st_1 = 1, st_2 = 2

contains


  !> Creates a conjugate gradient instance
  subroutine conjGrad_init(this, nElem, tol, maxDisp)

    !> conjugate gradient instance on exit
    type(TConjGrad), intent(out) :: this

    !> Nr. of elements in the vectors
    integer, intent(in) :: nElem

    !> Termination tolerance for the gradient
    real(dp), intent(in) :: tol

    !> Maximal displacement in one element in one step
    real(dp), intent(in) :: maxDisp

    @:ASSERT(nElem > 0)
    @:ASSERT(tol > 0.0_dp)
    @:ASSERT(maxDisp > 0.0_dp)

    this%nElem = nElem
    this%tolerance = tol
    this%maxDisp = maxDisp
    allocate(this%gg(nElem))
    allocate(this%hh(nElem))
    allocate(this%uu(nElem))
    !! Line minimizer is created with an extrem big tolerance: it just brackets
    !! the minimum and returns an approximative minimum between them. Seems
    !! to give in most cases better results as making many line min. steps.
    call TLineMin_init(this%pLinMin, nElem, 10, 10000.0_dp, this%maxDisp)
    this%tInitialized = .false.

  end subroutine conjGrad_init


  !> Resets CG minimizer
  subroutine conjGrad_reset(this, x0)

    !> CG minimizer
    type(TConjGrad), intent(inout) :: this

    !> Point to start from
    real(dp), intent(in) :: x0(:)

    @:ASSERT(size(x0) == this%nElem)

    this%uu(:) = x0(:)
    this%state = st_1
    this%tConverged = .false.
    this%tInitialized = .true.

  end subroutine conjGrad_reset


  !> Passes calculated function value and gradient to the minimizare and gives a new coordinate
  !> back.  When calling the first time, funciton value and gradient for the starting point of the
  !> minimization should be passed.
  subroutine conjGrad_next(this, fx, dx, xNew, tConverged)

    !> CG minimizer
    type(TConjGrad), intent(inout) :: this

    !> Function value for last point returned by this routine
    real(dp), intent(in) :: fx

    !> Gradient in the last point
    real(dp), intent(in) :: dx(:)

    !> New proposed point
    real(dp), intent(out) :: xNew(:)

    !> True, if gradient goes below the specified tolerance
    logical,  intent(out) :: tConverged

    @:ASSERT(this%tInitialized)
    @:ASSERT(size(xNew) == this%nElem)
    @:ASSERT(size(dx) == this%nElem)

    if (.not. this%tConverged) then
      call next_local(this%state, this%gg, this%hh, this%uu, this%tConverged, this%tolerance,&
          & this%pLinMin, fx, dx)
    end if
    xNew(:) = this%uu(:)
    tConverged = this%tConverged

  end subroutine conjGrad_next


  !> Workhorse for the CG minimizer.
  subroutine next_local(state, gg, hh, uu, tConverged, tolerance, pLinMin, fu, du)

    !> state of the conjugate gradient minimizer
    integer, intent(inout) :: state

    !> Gradient in previous cycle
    real(dp), intent(inout) :: gg(:)

    !> conjugate gradient
    real(dp), intent(inout) :: hh(:)

    !> New calculated point
    real(dp), intent(inout) :: uu(:)

    !> Did convergence occur?
    logical,  intent(inout) :: tConverged

    !> cut out tolerance for optimisation
    real(dp), intent(in) :: tolerance

    !> Line minimizer
    type(TLineMin), intent(inout) :: pLinMin

    !> Function value in the last point
    real(dp), intent(in) :: fu

    !> Gradient in the last point
    real(dp), intent(in) :: du(:)

    real(dp) :: ggAbs, dgAbs, rTmp
    logical :: tConvLine
    real(dp), allocatable :: xi(:)

    if (state == st_1) then
      !! If first gradient converged: Reuse internal variables to store results
      if (maxval(abs(du)) < tolerance) then
        tConverged = .true.
        gg(:) = du(:)
        hh(:) = 0.0_dp
        hh(1) = fu
      else
        gg(:) = -du(:)
        hh(:) = gg(:)
      end if
      !! First step x = F/k, where F is the acting force and
      !! k is a spring constant in the magnitude of a C-C vibration mode
      call pLinMin%reset(uu, hh, 5.0_dp * sqrt(sum(gg**2)))
      state = st_2
    end if

    tConvLine = .true.
    do while ((.not. tConverged) .and. tConvLine)
      call pLinMin%next(fu, du, uu, tConvLine)
      if (tConvLine) then
        allocate(xi(size(gg)))
        call pLinMin%getMinGrad(xi)
        if (maxval(abs(xi)) < tolerance) then
          tConverged = .true.
        else
          ggAbs = dot_product(gg, gg)
          dgAbs = dot_product(xi+gg, xi)
          if (ggAbs < epsilon(0.0_dp)) then
            tConverged = .true.
          else
            gg(:) = -xi(:)
            hh(:) = gg(:) + (dgAbs / ggAbs) * hh(:)
            call pLinMin%getMinX(uu)
            call pLinMin%getMinLambda(rTmp)
            call pLinMin%reset(uu, hh, rTmp)
          end if
        end if
      end if
    end do

    !! If converged, reuse internal variables to store the result
    if (tConverged .and. state /= st_1) then
      call pLinMin%getMinX(uu)
      call pLinMin%getMinGrad(gg)
      call pLinMin%getMinY(hh(1))
    end if

  end subroutine next_local


  !> Gives the coordinate of the minimal point back
  !> The returned value is meaningless if the subroutine is called
  !>   before the CG minimizer signals convergence.
  subroutine conjGrad_getMinX(this, minX)

    !> CG minimizer
    type(TConjGrad), intent(in) :: this

    !> Coordinate of the minimal point
    real(dp), intent(out) :: minX(:)

    @:ASSERT(this%tInitialized .and. this%tConverged)
    @:ASSERT(size(minX) == this%nElem)
    call this%pLinMin%getMinX(minX)

  end subroutine conjGrad_getMinX


  !> Gives the function value in the minimal point back.
  !> The returned value is meaningless if the subroutine is called before the CG minimizer
  !>   signals convergence.
  subroutine conjGrad_getMinY(this, minY)

    !> CG minimizer
    type(TConjGrad), intent(in) :: this

    !> Coordinate of the minimal point
    real(dp), intent(out) :: minY

    @:ASSERT(this%tInitialized .and. this%tConverged)
    call this%pLinMin%getMinY(minY)

  end subroutine conjGrad_getMinY


  !> Gives the gradient in the minimal point back
  !> The returned value is meaningless if the subroutine is called before the CG minimizer signals
  !> convergence.
  subroutine conjGrad_getMinGrad(this, minGrad)

    !> CG minimizer
    type(TConjGrad), intent(in) :: this

    !> Coordinate of the minimal point
    real(dp), intent(out) :: minGrad(:)

    @:ASSERT(this%tInitialized .and. this%tConverged)
    @:ASSERT(size(minGrad) == this%nElem)
    call this%pLinMin%getMinGrad(minGrad)

  end subroutine conjGrad_getMinGrad

end module dftbp_conjgrad
