!--------------------------------------------------------------------------------------------------!
!  DFTB+: general package for performing fast atomistic simulations                                !
!  Copyright (C) 2006 - 2025  DFTB+ developers group                                               !
!                                                                                                  !
!  See the LICENSE file for terms of usage and distribution.                                       !
!--------------------------------------------------------------------------------------------------!

#:include 'common.fypp'

!> Contains a geometry DIIS optimizer interface.
module dftbp_geoopt_gdiis
  use dftbp_common_accuracy, only : dp
  use dftbp_mixer_diismixer, only : TDiisMixerInp, TDiisMixerReal, TDiisMixerReal_init
  implicit none

  private

  public :: TDiis
  public :: init, reset, next


  !> Contains data for the DIIS mimimizer
  type TDiis
    private

    !> DIIS object itself
    type(TDiisMixerReal) :: diis

    !> Vector of current coordinate
    real(dp), allocatable :: x(:)

    !> Number of elements in vector
    integer :: nElem

    !> Tolerance criteria for convergence
    real(dp) :: tolerance

    !> If object is initialized
    logical :: tInitialized
  end type TDiis


  !> Creates gDIIS instance
  interface init
    module procedure gDiis_init
  end interface


  !> Resets the gDIIS instance
  interface reset
    module procedure gDiis_reset
  end interface


  !> Passes calculated gradient to the minimizer and returns a new point
  interface next
    module procedure gDiis_next
  end interface

contains


  !> Creates a DIIS geometry optimiser instance.
  subroutine gDiis_init(this, nElem, tol, alpha, nGens)

    !> DIIS instance on exit
    type(TDiis), intent(out) :: this

    !> Nr. of elements in the vectors
    integer, intent(in) :: nElem

    !> Termination tolerance for the gradient
    real(dp), intent(in) :: tol

    !> Initial value for mixing in gradient information to DIIS space
    real(dp), intent(in) :: alpha

    !> Number of vectors to use in building DIIS space
    integer, intent(in) :: nGens

    type(TDiisMixerInp) :: mixerInp

    this%nElem = nElem
    this%tolerance = tol
    allocate(this%x(this%nElem))

    mixerInp = TDiisMixerInp(nGens, alpha, .true., alpha)
    call TDiisMixerReal_init(this%diis, mixerInp)
    this%tInitialized = .true.

  end subroutine gDiis_init


  !> Resets optimiser.
  subroutine gDiis_reset(this, x)

    !> Minimiser
    type(TDiis), intent(inout) :: this

    !> Point to start from
    real(dp) :: x(:)

    call this%diis%reset(this%nElem)
    this%x(:) = x

  end subroutine gDiis_reset


  !> Passes calculated function value and gradient to the minimizare and gives a new coordinate
  !! back.  When calling the first time, function value and gradient for the starting point of the
  !! minimization should be passed.
  subroutine gDiis_next(this, dx, xNew, tConverged)

    !> Minimiser
    type(TDiis), intent(inout) :: this

    !> Gradient in the last point
    real(dp), intent(in) :: dx(:)

    !> New proposed point
    real(dp), intent(out) :: xNew(:)

    !> True, if gradient goes below the specified tolerance
    logical,  intent(out) :: tConverged

    @:ASSERT(this%tInitialized)
    @:ASSERT(size(xNew) == this%nElem)
    @:ASSERT(size(dx) == this%nElem)

    xNew(:) = this%x
    call this%diis%mix1d(this%x, dx)
    tConverged = maxval(abs(xNew - this%x)) < this%tolerance .or. (maxval(abs(dx)) < this%tolerance)
    xNew(:) = this%x

  end subroutine gDiis_next

end module dftbp_geoopt_gdiis
