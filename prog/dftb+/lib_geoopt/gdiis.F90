!--------------------------------------------------------------------------------------------------!
!  DFTB+: general package for performing fast atomistic simulations                                !
!  Copyright (C) 2006 - 2020  DFTB+ developers group                                               !
!                                                                                                  !
!  See the LICENSE file for terms of usage and distribution.                                       !
!--------------------------------------------------------------------------------------------------!

#:include 'common.fypp'

!> Contains a geometry DIIS optimizer interface.
module dftbp_gdiis
  use dftbp_assert
  use dftbp_accuracy
  use dftbp_diismixer
  implicit none

  private


  !> Contains data for the DIIS mimimizer
  type TDIIS
    private

    !> DIIS object itself
    type(TDIISMixer) :: pDIIS

    !> Vector of current coordinate
    real(dp), allocatable :: x(:)

    !> Number of elements in vector
    integer :: nElem

    !> Tolerance criteria for convergence
    real(dp) :: tolerance

    !> If object is initialized
    logical :: tInitialized
  end type TDIIS


  !> Creates gDIIS instance
  interface init
    module procedure gDIIS_init
  end interface


  !> Resets the gDIIS instance
  interface reset
    module procedure gDIIS_reset
  end interface


  !> Passes calculated gradient to the minimizer and returns a new point
  interface next
    module procedure gDIIS_next
  end interface

  public :: TDIIS
  public :: init, reset, next

contains


  !> Creates a DIIS geometry optimiser instance
  subroutine gDIIS_init(self, nElem, tol, alpha, nGens)

    !> DIIS instance on exit
    type(TDIIS), intent(out) :: self

    !> Nr. of elements in the vectors
    integer, intent(in) :: nElem

    !> Termination tolerance for the gradient
    real(dp), intent(in) :: tol

    !> initial value for mixing in gradient information to DIIS space
    real(dp), intent(in) :: alpha

    !> Number of vectors to use in building DIIS space
    integer, intent(in) :: nGens

    self%nElem = nElem
    self%tolerance = tol
    allocate(self%x(self%nElem))
    call init(self%pDIIS,nGens,alpha,.true.,alpha)
    self%tInitialized = .true.

  end subroutine gDIIS_init


  !> Resets optimiser
  subroutine gDIIS_reset(self,x)

    !> Minimiser
    type(TDIIS), intent(inout) :: self

    !> Point to start from
    real(dp) :: x(:)

    call reset(self%pDIIS, self%nElem)
    self%x(:) = x(:)

  end subroutine gDIIS_reset


  !> Passes calculated function value and gradient to the minimizare and gives a new coordinate
  !> back.  When calling the first time, funciton value and gradient for the starting point of the
  !> minimization should be passed.
  subroutine gDIIS_next(self,dx, xNew, tConverged)

    !> minimiser
    type(TDIIS), intent(inout) :: self

    !> Gradient in the last point
    real(dp), intent(in) :: dx(:)

    !> New proposed point
    real(dp), intent(out) :: xNew(:)

    !> True, if gradient goes below the specified tolerance
    logical,  intent(out) :: tConverged

    @:ASSERT(self%tInitialized)
    @:ASSERT(size(xNew) == self%nElem)
    @:ASSERT(size(dx) == self%nElem)

    xNew = self%x
    call mix(self%pDIIS,self%x,dx)
    if (maxval(abs(xNew-self%x)) < self%tolerance &
        & .or. (maxval(abs(dx)) < self%tolerance)) then
      tConverged = .true.
    else
      tConverged = .false.
    end if
    xNew = self%x

  end subroutine gDIIS_next

end module dftbp_gdiis
