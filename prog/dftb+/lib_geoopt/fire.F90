!--------------------------------------------------------------------------------------------------!
!  DFTB+: general package for performing fast atomistic simulations                                !
!  Copyright (C) 2006 - 2020  DFTB+ developers group                                               !
!                                                                                                  !
!  See the LICENSE file for terms of usage and distribution.                                       !
!--------------------------------------------------------------------------------------------------!

!> FIRE optimiser
module dftbp_fire
  use dftbp_accuracy, only : dp
  implicit none
  private
  public :: TFire, TFire_init

  type TFire

    !> Iteration count
    integer :: iter

    !> Last reset in values
    integer :: resetStep

    !> Minimum steps before step size changes are allowed
    integer :: nMin

    !> Dynamical scaling factor
    real(dp) :: a

    !> Initial scaling choice
    real(dp) :: a_start

    !> Decrement scaling factor for forces
    real(dp) :: f_inc

    !> Increment scaling factor for forces
    real(dp) :: f_dec

    !> Scaling factor for dynamical variable
    real(dp) :: f_alpha

    !> Starting time step
    real(dp) :: dt_init

    !> current time step
    real(dp) :: dt

    !> Maximum step
    real(dp) :: dt_max

    !> force tolerance for halting
    real(dp) :: tol

    !> Position
    real(dp), allocatable :: x(:)

    !> Velocity
    real(dp), allocatable :: velocity(:)

  contains

    procedure :: reset
    procedure :: next

  end type TFire

contains

  !> Initialise type
  subroutine TFire_init(this, nElem, tol, maxStep)

    !> instance
    type(TFire), intent(out) :: this

    !> Number of elements to optimize
    integer, intent(in) :: nElem

    !> Halting tolerance in forces
    real(dp), intent(in) :: tol

    !> Maximum step size
    real(dp), intent(in) :: maxStep

    this%tol = tol

    if (allocated(this%velocity)) then
      if (size(this%velocity) /= nElem) then
        deallocate(this%velocity)
        allocate(this%velocity(nElem))
      end if
    else
      allocate(this%velocity(nElem))
    end if
    this%velocity(:) = 0.0_dp
    if (allocated(this%x)) then
      if (size(this%x) /= nElem) then
        deallocate(this%x)
        allocate(this%x(nElem))
      end if
    else
      allocate(this%x(nElem))
    end if
    this%x(:) = 0.0_dp

    ! default values from the paper
    this%nMin = 5
    this%a_start = 0.1_dp
    this%f_inc = 1.1_dp
    this%f_dec = 0.5_dp
    this%f_alpha = 0.99_dp
    this%dt_max = maxStep
    this%dt_init = 0.1_dp * this%dt_max

  end subroutine TFire_init

  !> Reset the integrator state
  subroutine reset(this, x0)

    !> Instance
    class(TFire), intent(inout) :: this

    !> Coordinates
    real(dp), intent(in) :: x0(:)

    this%iter = 0
    this%resetStep = 0
    this%a = this%a_start
    this%dt = this%dt_init
    this%velocity(:) = 0.0_dp
    this%x(:) = x0
    
  end subroutine reset


  !> Leap frog Verlet integrator with fire modification
  subroutine next(this, dx, xNew, isConverged)

    !> Instance
    class(TFire), intent(inout) :: this

    !> New point
    real(dp), intent(inout) :: xNew(:)

    !> Function derivative
    real(dp), intent(in) :: dx(:)

    !> Has the optimization completed
    logical, intent(out) :: isConverged

    if( maxval(abs(dx)) < this%tol ) then
      isConverged = .true.
    else
      call fireModifyVelocity(this, dx)
      ! all masses 1, so f = a, note that f = - dx
      this%velocity(:) = this%velocity - this%dt * dx
      this%x(:) = this%x + this%velocity * this%dt
      isConverged = .false.
    end if

    xNew(:) = this%x

    this%iter = this%iter + 1

  end subroutine next


  !> https://doi.org/10.1103/PhysRevLett.97.170201
  subroutine fireModifyVelocity(this, dx)

    !> Instance
    class(TFire), intent(inout) :: this

    !> Derivative of function, note dx = - f, so sign flip in expressions
    real(dp), intent(in) :: dx(:)

    real(dp) :: p

    p = -1.0_dp * dot_product(this%velocity,dx)
    this%velocity(:) = (1.0_dp-this%a)*this%velocity - this%a*dx*mag(this%velocity)/mag(dx)
    if( p < 0.0_dp ) then
      this%velocity(:) = 0.0_dp
      this%resetStep = this%iter
      this%dt = this%dt * this%f_dec
      this%a = this%a_start
    else if (this%iter - this%resetStep > this%nMin) then
      this%dt = min(this%dt * this%f_inc, this%dt_max)
      this%a = this%a * this%f_alpha
    end if

  end subroutine fireModifyVelocity


  !> magnitude of a vector
  pure function mag(v)

    !> Vector
    real(dp), intent(in) :: v(:)

    !> Its magnitude
    real(dp) :: mag

    mag = sqrt(sum(v**2))

  end function mag

end module dftbp_fire
