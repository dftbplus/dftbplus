!--------------------------------------------------------------------------------------------------!
!  DFTB+: general package for performing fast atomistic simulations                                !
!  Copyright (C) 2006 - 2023  DFTB+ developers group                                               !
!                                                                                                  !
!  See the LICENSE file for terms of usage and distribution.                                       !
!--------------------------------------------------------------------------------------------------!

!> Exports the functionality of the geometry optimizers as one package
module dftbp_geoopt_package
  use dftbp_common_accuracy, only : dp
  use dftbp_geoopt_filter, only : TFilterInput, TFilter
  use dftbp_geoopt_fire, only : TFireInput, TFire, TFire_init
  use dftbp_geoopt_lbfgs2, only : TLbfgsInput, TLbfgs, TLbfgs_init
  use dftbp_geoopt_optimizer, only : TOptimizer, TOptimizerInput
  use dftbp_geoopt_rationalfunc, only : TRationalFuncInput, TRationalFunc, TRationalFunc_init
  use dftbp_geoopt_steepdesc, only : TSteepdescInput, TSteepdesc, TSteepdesc_init
  implicit none

  private
  public :: TOptimizer, TOptimizerInput
  public :: TSteepdescInput, TSteepdesc, TSteepdesc_init
  public :: TFireInput, TFire, TFire_init
  public :: TLbfgsInput, TLbfgs, TLbfgs_init
  public :: TRationalFuncInput, TRationalFunc, TRationalFunc_init
  public :: TFilterInput, TFilter
  public :: TOptTolerance
  public :: createOptimizer


  !> Tolerances for optimization
  type :: TOptTolerance

    !> Convergence threshold for energy
    real(dp) :: energy = huge(1.0_dp)

    !> Convergence threshold for gradient norm
    real(dp) :: gradNorm = huge(1.0_dp)

    !> Convergence threshold for gradient norm
    real(dp) :: gradElem = huge(1.0_dp)

    !> Convergence threshold for displacement norm
    real(dp) :: dispNorm = huge(1.0_dp)

    !> Convergence threshold for displacement norm
    real(dp) :: dispElem = huge(1.0_dp)

  end type TOptTolerance


contains

  !> Create new polymorphic optimizer object
  subroutine createOptimizer(input, nVar, optimizer)

    !> Input for defining the optimizer object
    class(TOptimizerInput), intent(inout) :: input

    !> Number of variables to optimize
    integer, intent(in) :: nVar

    !> Instance of the optimizer object
    class(TOptimizer), allocatable, intent(out) :: optimizer

    select type (input)

    type is (TLbfgsInput)
      block
        type(TLbfgs), allocatable :: tmp
        allocate(tmp)
        call TLbfgs_init(tmp, input, nVar)
        call move_alloc(tmp, optimizer)
      end block
      return

    type is (TRationalFuncInput)
      block
        type(TRationalFunc), allocatable :: tmp
        allocate(tmp)
        call TRationalFunc_init(tmp, input, nVar)
        call move_alloc(tmp, optimizer)
      end block
      return

    type is (TFireInput)
      block
        type(TFire), allocatable :: tmp
        allocate(tmp)
        call TFire_init(tmp, input, nVar)
        call move_alloc(tmp, optimizer)
      end block
      return

    type is (TSteepdescInput)
      block
        type(TSteepdesc), allocatable :: tmp
        allocate(tmp)
        call TSteepdesc_init(tmp, input, nVar)
        call move_alloc(tmp, optimizer)
      end block
      return

    end select

  end subroutine createOptimizer


end module dftbp_geoopt_package
