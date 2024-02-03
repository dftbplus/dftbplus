!--------------------------------------------------------------------------------------------------!
!  DFTB+: general package for performing fast atomistic simulations                                !
!  Copyright (C) 2006 - 2024  DFTB+ developers group                                               !
!                                                                                                  !
!  See the LICENSE file for terms of usage and distribution.                                       !
!--------------------------------------------------------------------------------------------------!

#:include 'common.fypp'

module dftbp_geoopt_lbfgs2
  use dftbp_common_accuracy, only : dp
  use dftbp_geoopt_optimizer, only : TOptimizer, TOptimizerInput
  implicit none
  private

  public :: TLbfgsInput, TLbfgs, TLbfgs_init


  !> Input for the LBFGS optimizer
  type, extends(TOptimizerInput) :: TLbfgsInput

    !> Memory limit for the LBFGS update
    integer :: memory = 20

  end type TLbfgsInput


  !> Limited-memory BFGS optimizer
  type, extends(TOptimizer) :: TLbfgs

    !> Current iteration step
    integer :: iter

    !> Number of variables to optimize
    integer :: nvar

    !> Memory limit for the LBFGS update
    integer :: memory

    !> Last gradient
    real(dp), allocatable :: gLast(:)

    !> *Inverse* Hessian in diagonal form
    real(dp), allocatable :: hdiag(:)

    !> LBFGS scratch array of former displacements
    real(dp), allocatable :: s(:,:)

    !> LBFGS scratch array of former gradient changes
    real(dp), allocatable :: y(:,:)

    !> LBFGS scratch array of dot products between s and y
    real(dp), allocatable :: rho(:)

  contains

    !> Calculate displacement from gradient
    procedure :: step

    !> Reset optimizer
    procedure :: reset

  end type TLbfgs


contains


  !> Create new limited memory BFGS optimization driver
  subroutine TLbfgs_init(this, input, nVar)

    !> Instance of the optimizer
    type(TLbfgs), intent(out) :: this

    !> Input for the LBFGS optimizer
    type(TLbfgsInput), intent(in) :: input

    !> Number of variables to optimize
    integer, intent(in) :: nVar

    this%iter = 0
    this%nvar = nVar
    this%memory = input%memory
    allocate(this%gLast(this%nVar), source=0.0_dp)
    allocate(this%s(this%nvar, input%memory), source=0.0_dp)
    allocate(this%y(this%nvar, input%memory), source=0.0_dp)
    allocate(this%rho(input%memory), source=0.0_dp)
    allocate(this%hdiag(this%nvar), source=1.0_dp)

  end subroutine TLbfgs_init


  !> Calculate displacement from gradient
  subroutine step(this, val, grad, displ)

    !> Instance of geometry optimization driver
    class(TLbfgs), intent(inout) :: this

    !> Current function value
    real(dp), intent(in) :: val

    !> Current gradient
    real(dp), intent(in) :: grad(:)

    !> Next displacement step
    real(dp), intent(out) :: displ(:)

    this%iter = this%iter + 1

    call lbfgs_step(this%iter, this%memory, this%nvar, grad, this%glast, displ, this%hdiag, &
        & this%s, this%y, this%rho)
    this%gLast(:) = grad

  end subroutine step


  !> Reset optimizer
  subroutine reset(this)

    !> Instance of geometry optimization driver
    class(TLbfgs), intent(inout) :: this

    this%iter = 0
    this%gLast(:) = 0.0_dp
    this%hdiag(:) = 1.0_dp
    this%s(:,:) = 0.0_dp
    this%y(:,:) = 0.0_dp
    this%rho(:) = 0.0_dp

  end subroutine reset


  !> Updates displacement using the formula given by Nocedal, generally known
  !> as limited memory BFGS algorithm
  subroutine lbfgs_step(iter, memory, nvar, gradient, glast, displacement, hdiag, s, y, rho)
    !> Current iteration step
    integer, intent(in) :: iter
    !> Memory limit for the LBFGS update
    integer, intent(in) :: memory
    !> 3*natoms
    integer, intent(in) :: nvar
    !> current gradient
    real(dp), intent(in) :: gradient(:)
    !> gradient on the last point
    real(dp), intent(in) :: glast(:)
    !> on input displacement from the last step, on exit new displacement vector
    real(dp), intent(inout) :: displacement(:)
    !> *inverse* Hessian in diagonal form
    real(dp), intent(in) :: hdiag(:)

    !> LBFGS scratch array of former displacements
    real(dp), intent(inout) :: s(:, :)
    !> LBFGS scratch array of former gradient changes
    real(dp), intent(inout) :: y(:, :)
    !> LBFGS scratch array of dot products between s and y
    real(dp), intent(inout) :: rho(:)

    real(dp), allocatable :: d(:), q(:), a(:)
    real(dp) :: b

    integer :: thisiter, lastiter, mem
    integer :: i

    allocate(q(nvar), d(nvar), a(memory), source = 0.0_dp)

    thisiter = mod(iter-1, memory) + 1
    s(:, thisiter) = displacement
    y(:, thisiter) = gradient - glast

    b = dot_product(s(:, thisiter), y(:, thisiter))

    rho(thisiter) = 1.0_dp / max(abs(b), epsilon(1.0_dp))

    q = gradient
    do mem = iter, max(1, iter-memory), -1
      i = mod(mem-1, memory)+1
      a(i) = rho(i) * dot_product(s(:, i), q)
      q = q - a(i) * y(:, i)
    end do

    d = Hdiag * q

    do mem = max(1, iter-memory), iter
      i = mod(mem-1, memory)+1
      b = rho(i) * dot_product(y(:, i), d)
      d = d + s(:, i) * (a(i) - b)
    end do

    displacement = -d
  end subroutine lbfgs_step


end module dftbp_geoopt_lbfgs2
