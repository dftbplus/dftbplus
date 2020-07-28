!--------------------------------------------------------------------------------------------------!
!  DFTB+: general package for performing fast atomistic simulations                                !
!  Copyright (C) 2006 - 2020  DFTB+ developers group                                               !
!                                                                                                  !
!  See the LICENSE file for terms of usage and distribution.                                       !
!--------------------------------------------------------------------------------------------------!

#:include 'common.fypp'

!> Function minimization with steepest descent algorithm
module dftbp_steepdesc
  use dftbp_assert
  use dftbp_accuracy
  implicit none

  private


  !> Contains data for the steepest descent minimizer
  type TSteepDesc
    private

    !> Dimensionality of the space
    integer :: nElem

    !> Previous coordinate
    real(dp), allocatable :: xOld(:)

    !> Weights for the coordinates
    real(dp), allocatable :: weight(:)

    !> Tolerance criteria for convergence
    real(dp) :: tolerance

    !> Maximal displacement along one coordinate in one step
    real(dp) :: maxDisp

    !> If CG converged
    logical :: tConverged

    !> If object is initialized
    logical :: tInitialized
  end type TSteepDesc


  !> Creates SD instance
  interface init
    module procedure SteepDesc_init
  end interface


  !> Resets the SD instance
  interface reset
    module procedure SteepDesc_reset
  end interface


  !> Passes calculated gradient to the minimizer and returns a new point
  interface next
    module procedure SteepDesc_next
  end interface

  public :: TSteepDesc
  public :: init, reset, next

contains


  !> Initialises a steepest descent instance
  subroutine SteepDesc_init(this, nElem, tol, maxDisp, weight)

    !> Steepest descent instance on exit
    type(TSteepDesc), intent(out) :: this

    !> Nr. of elements in the vectors
    integer, intent(in) :: nElem

    !> Tolerance for the gradient
    real(dp), intent(in) :: tol

    !> Maximal displacement in one element in one step
    real(dp), intent(in) :: maxDisp

    !> The weights of the gradient components
    real(dp), intent(in) :: weight(:)

    @:ASSERT(nElem > 0)
    @:ASSERT(tol > 0.0_dp)
    @:ASSERT(maxDisp > 0.0_dp)
    @:ASSERT(size(weight) == nElem)

    this%nElem = nElem
    this%tolerance = tol
    this%maxDisp = maxDisp
    allocate(this%weight(nElem))
    this%weight(:) = weight(:)
    allocate(this%xOld(nElem))
    this%tInitialized = .false.

  end subroutine SteepDesc_init


  !> Resets CG minimizer
  subroutine SteepDesc_reset(this, x0)

    !> minimizer object
    type(TSteepDesc), intent(inout) :: this

    !> Point to start from
    real(dp), intent(in) :: x0(:)

    @:ASSERT(size(x0) == this%nElem)

    this%xOld(:) = x0(:)
    this%tConverged = .false.
    this%tInitialized = .true.

  end subroutine SteepDesc_reset


  !> Passes calculated function value and gradient to the minimizare and gives a new coordinate
  !> back.
  !> When calling the first time, gradient for the starting point of the minimization should be
  !> passed.
  subroutine SteepDesc_next(this, dx, xNew, tConverged)

    !> CG minimizer
    type(TSteepDesc), intent(inout) :: this

    !> Gradient in the last point
    real(dp), intent(in) :: dx(:)

    !> New proposed point
    real(dp), intent(out) :: xNew(:)

    !> True, if gradient got below the specified tolerance.
    logical,  intent(out) :: tConverged

    @:ASSERT(this%tInitialized)
    @:ASSERT(size(xNew) == this%nElem)
    @:ASSERT(size(dx) == this%nElem)

    if (.not. this%tConverged) then
      call next_local(xNew, this%xOld, dx, this%weight, this%maxDisp, &
          &this%tolerance, this%tConverged)
    else
      xNew(:) = this%xOld(:)
    end if
    tConverged = this%tConverged

  end subroutine SteepDesc_next


  !> Workhorse for the SD minimizer
  subroutine next_local(xNew, xOld, grad, weight, maxDisp, tolerance, &
      &tConverged)

    !> Coordinates of the new point on exit
    real(dp), intent(out) :: xNew(:)

    !> Coordinates of the previous point
    real(dp), intent(inout) :: xOld(:)

    !> Gradient at xOld
    real(dp), intent(in) :: grad(:)

    !> Weighting factors for the gradient components
    real(dp), intent(in) :: weight(:)

    !> Maximal displacement along one component
    real(dp), intent(in) :: maxDisp

    !> Termination tolerance
    real(dp), intent(in) :: tolerance

    !> Completion of the minimiser
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

end module dftbp_steepdesc
