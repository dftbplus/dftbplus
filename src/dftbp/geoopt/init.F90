!--------------------------------------------------------------------------------------------------!
!  DFTB+: general package for performing fast atomistic simulations                                !
!  Copyright (C) 2006 - 2021  DFTB+ developers group                                               !
!                                                                                                  !
!  See the LICENSE file for terms of usage and distribution.                                       !
!--------------------------------------------------------------------------------------------------!

!> Proxy module for initialization of optimizer objects
module dftbp_geoopt_init
  use dftbp_common_accuracy, only : dp
  use dftbp_geoopt_class, only : TOptimizer
  use dftbp_geoopt_filter, only : TFilter
  use dftbp_geoopt_input, only : TGeoOptInput, TOptConv
  implicit none
  private
  public :: TOptimizer, TOptimizer_init, TOptConv


contains


  !> Create new polymorphic optimizer object
  subroutine TOptimizer_init(this, input, filter)

    !> Instance of the optimizer object
    class(TOptimizer), allocatable, intent(out) :: this

    !> Input for defining the optimizer object
    type(TGeoOptInput), intent(in) :: input

    !> Geometry transformation filter
    type(TFilter), intent(in) :: filter

    if (allocated(input%lbfgs)) then
      block
        use dftbp_geoopt_lbfgs2, only : TLBFGS, TLBFGS_init
        type(TLBFGS), allocatable :: tmp
        allocate(tmp)
        call TLBFGS_init(tmp, input%lbfgs, filter)
        call move_alloc(tmp, this)
      end block
    end if

    if (allocated(input%rf)) then
      block
        use dftbp_geoopt_rf, only : TRationalFunction, TRationalFunction_init
        type(TRationalFunction), allocatable :: tmp
        allocate(tmp)
        call TRationalFunction_init(tmp, input%rf, filter)
        call move_alloc(tmp, this)
      end block
    end if

    if (allocated(input%fire)) then
      block
        use dftbp_geoopt_fire, only : TFire, TFire_init
        type(TFire), allocatable :: tmp
        allocate(tmp)
        call TFire_init(tmp, input%fire, filter)
        call move_alloc(tmp, this)
      end block
    end if
  end subroutine TOptimizer_init


end module dftbp_geoopt_init
