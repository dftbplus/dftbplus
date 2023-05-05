!--------------------------------------------------------------------------------------------------!
!  DFTB+: general package for performing fast atomistic simulations                                !
!  Copyright (C) 2006 - 2023  DFTB+ developers group                                               !
!                                                                                                  !
!  See the LICENSE file for terms of usage and distribution.                                       !
!--------------------------------------------------------------------------------------------------!

!> FIRE optimiser
module dftbp_geoopt_fire
  use dftbp_common_accuracy, only : dp
  use dftbp_io_message, only : error
  use dftbp_geoopt_optimizer, only : TOptimizer, TOptimizerInput
  implicit none

  private
  public :: TFireInput, TFire, TFire_init


  type, extends(TOptimizerInput) :: TFireInput

    !> Minimum steps before step size changes are allowed
    integer :: nMin

    !> Initial scaling choice
    real(dp) :: a_start

    !> Decrement scaling factor for forces
    real(dp) :: f_inc

    !> Increment scaling factor for forces
    real(dp) :: f_dec

    !> Scaling factor for dynamical variable
    real(dp) :: f_alpha

    !> Maximum step
    real(dp) :: dt_max

  end type TFireInput


  !> Fast inertial relaxation engine, FIRE.
  type, extends(TOptimizer) :: TFire

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
    procedure :: reset_old
    procedure :: next
    procedure :: step

  end type TFire

  interface TFire_init
    module procedure :: TFire_init
    module procedure :: TFire_init_old
  end interface TFire_init

contains

  !> Initialise type
  subroutine TFire_init(this, input, nVar)

    !> instance
    type(TFire), intent(out) :: this

    !> Input for the optimizer
    type(TFireInput), intent(in) :: input

    !> Number of variables to optimize
    integer, intent(in) :: nVar

    allocate(this%velocity(nVar))
    this%velocity(:) = 0.0_dp

    this%nMin = input%nMin
    this%a_start = input%a_start
    this%f_inc = input%f_inc
    this%f_dec = input%f_dec
    this%f_alpha = input%f_alpha
    this%dt_max = input%dt_max
    this%dt_init = 0.1_dp * input%dt_max

    this%iter = 0
    this%resetStep = 0
    this%a = this%a_start
    this%dt = this%dt_init

  end subroutine TFire_init


  !> Initialise type
  subroutine TFire_init_old(this, nElem, tol, maxStep)

    !> instance
    type(TFire), intent(out) :: this

    !> Number of elements to optimize
    integer, intent(in) :: nElem

    !> Halting tolerance in forces
    real(dp), intent(in) :: tol

    !> Maximum step size
    real(dp), intent(in) :: maxStep

    type(TFireInput) :: input

    ! default values from the paper
    input = TFireInput( &
        & nMin = 5, &
        & a_start = 0.1_dp, &
        & f_inc = 1.1_dp, &
        & f_dec = 0.5_dp, &
        & f_alpha = 0.99_dp, &
        & dt_max = maxStep &
        & )

    call TFire_init(this, input, nElem)

    this%tol = tol

    allocate(this%x(nElem))
    this%x(:) = 0.0_dp

  end subroutine TFire_init_old


  !> Reset the integrator state
  subroutine reset_old(this, xx)

    !> Instance
    class(TFire), intent(inout) :: this

    !> Coordinates
    real(dp), intent(in) :: xx(:)

    if (.not. allocated(this%x)) then
      call error("Trying to reset new FIRE optimizer with old routine.")
    end if

    this%iter = 0
    this%resetStep = 0
    this%a = this%a_start
    this%dt = this%dt_init
    this%velocity(:) = 0.0_dp
    this%x(:) = xx

  end subroutine reset_old


  !> Reset the integrator state
  subroutine reset(this)

    !> Instance
    class(TFire), intent(inout) :: this

    this%iter = 0
    this%resetStep = 0
    this%a = this%a_start
    this%dt = this%dt_init
    this%velocity(:) = 0.0_dp

  end subroutine reset


  !> Leap frog Verlet integrator with fire modification
  subroutine step(this, val, grad, displ)

    !> Instance
    class(TFire), intent(inout) :: this

    !> Current function value
    real(dp), intent(in) :: val

    !> Current gradient
    real(dp), intent(in) :: grad(:)

    !> Next displacement step
    real(dp), intent(out) :: displ(:)

    call fireModifyVelocity(this, grad)
    ! all masses 1, so f = a, note that f = - grad
    this%velocity(:) = this%velocity - this%dt * grad
    displ(:) = this%velocity * this%dt

    this%iter = this%iter + 1

  end subroutine step


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
      xNew(:) = 0.0_dp
      isConverged = .true.
    else
      call this%step(0.0_dp, dx, xNew)
      isConverged = .false.
    end if
    this%x(:) = this%x + xNew
    xNew(:) = this%x

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

end module dftbp_geoopt_fire
