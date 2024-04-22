!--------------------------------------------------------------------------------------------------!
!  DFTB+: general package for performing fast atomistic simulations                                !
!  Copyright (C) 2006 - 2023  DFTB+ developers group                                               !
!                                                                                                  !
!  See the LICENSE file for terms of usage and distribution.                                       !
!--------------------------------------------------------------------------------------------------!

!> Module defining abstract base class for geometry optimizers
module dftbp_geoopt_optimizer
  use dftbp_common_accuracy, only : dp
  implicit none

  private
  public :: TOptimizer, TOptimizerInput


  !> Abstract base class for optimization drivers
  type, abstract :: TOptimizer
  contains

    !> Calculate displacement from gradient
    procedure(step), deferred :: step

    !> Reset optimizer
    procedure(reset), deferred :: reset

  end type TOptimizer


  abstract interface
    !> Calculate displacement from gradient
    subroutine step(this, val, grad, displ)
      import :: TOptimizer, dp
      implicit none

      !> Instance of geometry optimization driver
      class(TOptimizer), intent(inout) :: this

      !> Current function value
      real(dp), intent(in) :: val

      !> Current gradient
      real(dp), intent(in) :: grad(:)

      !> Next displacement step
      real(dp), intent(out) :: displ(:)

    end subroutine step

    !> Reset optimizer
    subroutine reset(this)
      import :: TOptimizer
      implicit none

      !> Instance of geometry optimization driver
      class(TOptimizer), intent(inout) :: this

    end subroutine reset
  end interface


  !> Abstract base class for the optimizer input types
  type, abstract :: TOptimizerInput
  end type TOptimizerInput


end module dftbp_geoopt_optimizer
