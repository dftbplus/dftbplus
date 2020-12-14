!--------------------------------------------------------------------------------------------------!
!  DFTB+: general package for performing fast atomistic simulations                                !
!  Copyright (C) 2006 - 2020  DFTB+ developers group                                               !
!                                                                                                  !
!  See the LICENSE file for terms of usage and distribution.                                       !
!--------------------------------------------------------------------------------------------------!

!> General interface for the optimization algorithms
module dftbp_geoopt
  use dftbp_accuracy
  use dftbp_conjgrad
  use dftbp_steepdesc
  use dftbp_gdiis
  use dftbp_lbfgs
  use dftbp_fire
  implicit none
  private


  public :: Tgeoopt
  public :: init, reset, next
  public :: geoOptTypes

  !> Interface type for the various geometry optimization algorithms
  type Tgeoopt
    private
    integer :: iGeoOpt
    type(TConjGrad), allocatable :: pConjGrad
    type(TSteepDesc), allocatable :: pSteepDesc
    type(TDIIS), allocatable :: pDiis
    type(TLbfgs), allocatable :: pLbfgs
    type(TFire), allocatable :: pFire
  end type Tgeoopt


  !> Creates a geometry optimizer
  interface init
    module procedure GeoOpt_iniTConjGrad
    module procedure GeoOpt_iniTSteepDesc
    module procedure GeoOpt_iniTDIIS
    module procedure GeoOpt_initLbfgs
    module procedure GeoOpt_initFire
  end interface


  !> Resets the optimizer
  interface reset
    module procedure GeoOpt_reset
  end interface


  !> Delivers the next point in the minimization
  interface next
    module procedure GeoOpt_next
  end interface


  type :: TGeoOptTypesEnum
    integer :: none = 0
    integer :: steepestDesc = 1
    integer :: conjugateGrad = 2
    integer :: diis = 3
    integer :: lbfgs = 4
    integer :: fire = 5
  end type TGeoOptTypesEnum

  type(TGeoOptTypesEnum), parameter :: geoOptTypes = TGeoOptTypesEnum()


contains


  !> Creates a general geometry optimizier with a conjugate gradient instance
  subroutine GeoOpt_iniTConjGrad(this, pConjGrad)

    !> GeoOpt instance
    type(Tgeoopt), intent(out) :: this

    !> An already initialized conjugate gradient instance
    type(TConjGrad), allocatable, intent(inout) :: pConjGrad

    this%iGeoOpt = geoOptTypes%conjugateGrad
    call move_alloc(pConjGrad, this%pConjGrad)

  end subroutine GeoOpt_iniTConjGrad


  !> Creates a general geometry optimizier with a steepest descent instance
  subroutine GeoOpt_iniTSteepDesc(this, pSteepDesc)

    !> GeoOpt instance
    type(Tgeoopt), intent(out) :: this

    !> An already initialized steepest descent instance
    type(TSteepDesc), allocatable, intent(inout) :: pSteepDesc

    this%iGeoOpt = geoOptTypes%steepestDesc
    call move_alloc(pSteepDesc, this%pSteepDesc)

  end subroutine GeoOpt_iniTSteepDesc


  !> Creates a general geometry optimizier with a steepest descent instance
  subroutine GeoOpt_iniTDIIS(this, pDiis)

    !> GeoOpt instance
    type(Tgeoopt), intent(out) :: this

    !> An already initialized modified DIIS instance
    type(TDIIS), allocatable, intent(inout) :: pDiis

    this%iGeoOpt = geoOptTypes%diis
    call move_alloc(pDiis, this%pDiis)

  end subroutine GeoOpt_iniTDIIS

  !> Creates a general geometry optimizier with a limited memory BFGS driver
  subroutine GeoOpt_initLbfgs(this, pLbfgs)

    !> GeoOpt instance
    type(Tgeoopt), intent(out) :: this

    !> An already initialized modified LBFGS
    type(Tlbfgs), allocatable, intent(inout) :: pLbfgs

    this%iGeoOpt = geoOptTypes%lbfgs
    call move_alloc(pLbfgs, this%pLbfgs)

  end subroutine GeoOpt_initLbfgs


  !> Creates a general geometry optimizier with a FIRE instance
  subroutine GeoOpt_iniTFire(this, pFire)

    !> GeoOpt instance
    type(TGeoOpt), intent(out) :: this

    !> An already initialized conjugate gradient instance
    type(TFire), allocatable, intent(inout) :: pFire

    this%iGeoOpt = geoOptTypes%fire
    call move_alloc(pFire, this%pFire)

  end subroutine GeoOpt_iniTFire


  !> Resets the geometry optimizer
  subroutine GeoOpt_reset(this, x0)

    !> GeoOpt instance
    type(Tgeoopt), intent(inout) :: this

    !> Initial coordinates
    real(dp), intent(in) :: x0(:)

    select case (this%iGeoOpt)
    case(geoOptTypes%conjugateGrad)
      call reset(this%pConjGrad, x0)
    case(geoOptTypes%steepestDesc)
      call reset(this%pSteepDesc, x0)
    case(geoOptTypes%diis)
      call reset(this%pDiis, x0)
    case (geoOptTypes%lbfgs)
      call this%pLbfgs%reset(x0)
    case (geoOptTypes%fire)
      call this%pFire%reset(x0)
    end select

  end subroutine GeoOpt_reset


  !> Delivers the next point in the geometry optimization. When calling the first time, funciton
  !> value and gradient for the starting point of the minimization should be passed.
  subroutine GeoOpt_next(this, fx, dx, xNew, tConverged)

    !> Optimiser object
    type(Tgeoopt), intent(inout) :: this

    !> Function value for last point returned by this routine
    real(dp), intent(in) :: fx

    !> Gradient in the last point
    real(dp), intent(in) :: dx(:)

    !> New proposed point
    real(dp), intent(out) :: xNew(:)

    !> True, if gradient got below the specified tolerance.
    logical, intent(out) :: tConverged

    select case (this%iGeoOpt)
    case(geoOptTypes%conjugateGrad)
      call next(this%pConjGrad, fx, dx, xNew, tConverged)
    case (geoOptTypes%steepestDesc)
      call next(this%pSteepDesc, dx, xNew, tConverged)
    case (geoOptTypes%diis)
      call next(this%pDiis, dx, xNew, tConverged)
    case (geoOptTypes%lbfgs)
      call this%pLbfgs%next(fx, dx, xNew, tConverged)
    case (geoOptTypes%fire)
      call this%pFire%next(dx, xNew, tConverged)
    end select

  end subroutine GeoOpt_next

end module dftbp_geoopt
