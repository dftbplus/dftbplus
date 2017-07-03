!--------------------------------------------------------------------------------------------------!
!  DFTB+: general package for performing fast atomistic simulations                                !
!  Copyright (C) 2017  DFTB+ developers group                                                      !
!                                                                                                  !
!  See the LICENSE file for terms of usage and distribution.                                       !
!--------------------------------------------------------------------------------------------------!

#:include 'common.fypp'

!!* Function minimization with steepest descent algorithm
module steepdesc
  use assert
  use accuracy
  implicit none

  private


  !!* Contains data for the steepest descent minimizer
  type OSteepDesc
    private
    integer  :: nElem                  !* Dimensionality of the space
    real(dp), allocatable :: xOld(:)
    real(dp), allocatable :: weight(:)
    real(dp) :: tolerance              !* Tolerance criteria for convergence
    real(dp) :: maxDisp                !* Maximal displacement along one
                                       !* coordinate in one step
    logical :: tConverged              !* If CG converged
    logical :: tInitialized            !* If object is initialized
  end type OSteepDesc


  !!* Creates SD instance
  interface init
    module procedure SteepDesc_init
  end interface

  !!* Resets the SD instance
  interface reset
    module procedure SteepDesc_reset
  end interface

  !!* Passes calculated gradient to the minimizer and returns a new point
  interface next
    module procedure SteepDesc_next
  end interface


  public :: OSteepDesc
  public :: init, reset, next


contains

  !!* Initialises a steepest descent instance
  !!* @param self    Steepest descent instance on exit
  !!* @param nElem   Nr. of elements in the vectors
  !!* @param tol     Tolerance for the gradient
  !!* @param maxDisp Maximal displacement in one element in one step
  !!* @param weight  The weights of the gradient components
  subroutine SteepDesc_init(self, nElem, tol, maxDisp, weight)
    type(OSteepDesc), intent(out) :: self
    integer, intent(in) :: nElem
    real(dp), intent(in) :: tol
    real(dp), intent(in) :: maxDisp
    real(dp), intent(in) :: weight(:)

    @:ASSERT(nElem > 0)
    @:ASSERT(tol > 0.0_dp)
    @:ASSERT(maxDisp > 0.0_dp)
    @:ASSERT(size(weight) == nElem)

    self%nElem = nElem
    self%tolerance = tol
    self%maxDisp = maxDisp
    allocate(self%weight(nElem))
    self%weight(:) = weight(:)
    allocate(self%xOld(nElem))
    self%tInitialized = .false.

  end subroutine SteepDesc_init


  !!* Resets CG minimizer
  !!* @param self CG minimizer
  !!* @param x0   Point to start from
  subroutine SteepDesc_reset(self, x0)
    type(OSteepDesc), intent(inout) :: self
    real(dp), intent(in) :: x0(:)

    @:ASSERT(size(x0) == self%nElem)

    self%xOld(:) = x0(:)
    self%tConverged = .false.
    self%tInitialized = .true.

  end subroutine SteepDesc_reset



  !!* Passes calculated function value and gradient to the minimizare and
  !!* gives a new coordinate back.
  !!* @param self       CG minimizer
  !!* @param fx         Function value for last point returned by this routine
  !!* @param dx         Gradient in the last point
  !!* @param xNew       New proposed point
  !!* @param tConverged True, if gradient got below the specified tolerance.
  !!* @note When calling the first time, gradient for the starting point of the
  !!*    minimization should be passed.
  subroutine SteepDesc_next(self, dx, xNew, tConverged)
    type(OSteepDesc), intent(inout) :: self
    real(dp), intent(in)  :: dx(:)
    real(dp), intent(out) :: xNew(:)
    logical,  intent(out) :: tConverged

    @:ASSERT(self%tInitialized)
    @:ASSERT(size(xNew) == self%nElem)
    @:ASSERT(size(dx) == self%nElem)

    if (.not. self%tConverged) then
      call next_local(xNew, self%xOld, dx, self%weight, self%maxDisp, &
          &self%tolerance, self%tConverged)
    else
      xNew(:) = self%xOld(:)
    end if
    tConverged = self%tConverged

  end subroutine SteepDesc_next



  !!* Working horse for the SD minimizer
  !!* @param xNew    Coordinates of the new point on exit
  !!* @param xOld    Coordinates of the previous point
  !!* @param grad    Gradient in xOld
  !!* @param weight  Weighting factors for the gradient components
  !!* @param maxDisp Maximal displacement along one component
  subroutine next_local(xNew, xOld, grad, weight, maxDisp, tolerance, &
      &tConverged)
    real(dp), intent(out) :: xNew(:)
    real(dp), intent(inout) :: xOld(:)
    real(dp), intent(in) :: grad(:)
    real(dp), intent(in) :: weight(:)
    real(dp), intent(in) :: maxDisp
    real(dp), intent(in) :: tolerance
    logical,  intent(out) :: tConverged

    real(dp) :: maxX

    @:ASSERT(size(xNew) == size(xOld))
    @:ASSERT(size(xNew) == size(xOld))
    @:ASSERT(size(xNew) == size(xOld))
    @:ASSERT(maxDisp > 0.0_dp)

    if (maxval(abs(grad)) < tolerance) then
      tConverged = .true.
      xNew(:) = xOld(:)
      return
    else
      tConverged = .false.
    end if

    xNew(:) = -1.0_dp * weight(:) * grad(:)

    maxX = maxval(abs(xNew))
    if (maxX <= maxDisp) then
      xNew(:) = xOld(:) + xNew(:)
    else
      xNew(:) = xOld(:) + (maxDisp / maxX) * xNew(:)
    end if
    xOld(:) = xNew(:)

  end subroutine next_local


end module steepdesc
