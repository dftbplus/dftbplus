!--------------------------------------------------------------------------------------------------!
!  DFTB+: general package for performing fast atomistic simulations                                !
!  Copyright (C) 2006 - 2023  DFTB+ developers group                                               !
!                                                                                                  !
!  See the LICENSE file for terms of usage and distribution.                                       !
!--------------------------------------------------------------------------------------------------!

#:include 'common.fypp'

!> Implementation of a steepest descent optimization procedure.
module dftbp_geoopt_steepdesc
  use dftbp_common_accuracy, only : dp
  use dftbp_geoopt_optimizer, only : TOptimizer, TOptimizerInput
  implicit none

  private
  public :: TSteepdescInput, TSteepdesc, TSteepdesc_init


  !> Input for the steepest descent optimizer
  type, extends(TOptimizerInput) :: TSteepdescInput

    !> Scaling factor of gradients
    real(dp) :: scalingFactor

  end type TSteepdescInput


  !> Steepest descent optimization driver
  type, extends(TOptimizer) :: TSteepdesc

    !> Number of variables to optimize
    integer :: nVar

    !> Scaling factor of gradients
    real(dp) :: scalingFactor

  contains

    !> Calculate displacement from gradient
    procedure :: step

    !> Reset optimizer
    procedure :: reset

  end type TSteepdesc

contains


  !> Creates new steepest descent optimization driver.
  subroutine TSteepdesc_init(this, input, nVar)

    !> Instance of the optimizer
    type(TSteepdesc), intent(out) :: this

    !> Input for the steepest descent optimizer
    type(TSteepdescInput), intent(in) :: input

    !> Number of variables to optimize
    integer, intent(in) :: nVar

    @:ASSERT(nVar > 0)

    this%nVar = nVar
    this%scalingFactor = input%scalingFactor

  end subroutine TSteepdesc_init


  !> Calculates displacement from gradients.
  subroutine step(this, val, grad, displ)

    !> Instance of geometry optimization driver
    class(TSteepdesc), intent(inout) :: this

    !> Current function value
    real(dp), intent(in) :: val

    !> Current gradient
    real(dp), intent(in) :: grad(:)

    !> Next displacement step
    real(dp), intent(out) :: displ(:)

    @:ASSERT(size(grad) == this%nVar)
    @:ASSERT(size(displ) == this%nVar)

    displ(:) = -this%scalingFactor * grad

  end subroutine step


  !> Resets the optimizer.
  subroutine reset(this)

    !> Instance of geometry optimization driver
    class(TSteepdesc), intent(inout) :: this

  end subroutine reset


end module dftbp_geoopt_steepdesc
