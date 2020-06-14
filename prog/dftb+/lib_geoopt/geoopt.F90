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
  end type Tgeoopt


  !> Creates a geometry optimizer
  interface init
    module procedure GeoOpt_iniTConjGrad
    module procedure GeoOpt_iniTSteepDesc
    module procedure GeoOpt_iniTDIIS
    module procedure GeoOpt_initLbfgs
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
  end type TGeoOptTypesEnum

  type(TGeoOptTypesEnum), parameter :: geoOptTypes = TGeoOptTypesEnum()


contains


  !> Creates a general geometry optimizier with a conjugate gradient instance
  subroutine GeoOpt_iniTConjGrad(self, pConjGrad)

    !> GeoOpt instance
    type(Tgeoopt), intent(out) :: self

    !> An already initialized conjugate gradient instance
    type(TConjGrad), allocatable, intent(inout) :: pConjGrad

    self%iGeoOpt = geoOptTypes%conjugateGrad
    call move_alloc(pConjGrad, self%pConjGrad)

  end subroutine GeoOpt_iniTConjGrad


  !> Creates a general geometry optimizier with a steepest descent instance
  subroutine GeoOpt_iniTSteepDesc(self, pSteepDesc)

    !> GeoOpt instance
    type(Tgeoopt), intent(out) :: self

    !> An already initialized steepest descent instance
    type(TSteepDesc), allocatable, intent(inout) :: pSteepDesc

    self%iGeoOpt = geoOptTypes%steepestDesc
    call move_alloc(pSteepDesc, self%pSteepDesc)

  end subroutine GeoOpt_iniTSteepDesc


  !> Creates a general geometry optimizier with a steepest descent instance
  subroutine GeoOpt_iniTDIIS(self, pDiis)

    !> GeoOpt instance
    type(Tgeoopt), intent(out) :: self

    !> An already initialized modified DIIS instance
    type(TDIIS), allocatable, intent(inout) :: pDiis

    self%iGeoOpt = geoOptTypes%diis
    call move_alloc(pDiis, self%pDiis)

  end subroutine GeoOpt_iniTDIIS

  !> Creates a general geometry optimizier with a limited memory BFGS driver
  subroutine GeoOpt_initLbfgs(self, pLbfgs)

    !> GeoOpt instance
    type(Tgeoopt), intent(out) :: self

    !> An already initialized modified LBFGS
    type(Tlbfgs), allocatable, intent(inout) :: pLbfgs

    self%iGeoOpt = geoOptTypes%lbfgs
    call move_alloc(pLbfgs, self%pLbfgs)

  end subroutine GeoOpt_initLbfgs

  !> Resets the geometry optimizer
  subroutine GeoOpt_reset(self, x0)

    !> GeoOpt instance
    type(Tgeoopt), intent(inout) :: self

    !> Initial coordinates
    real(dp), intent(in) :: x0(:)

    select case (self%iGeoOpt)
    case(geoOptTypes%conjugateGrad)
      call reset(self%pConjGrad, x0)
    case(geoOptTypes%steepestDesc)
      call reset(self%pSteepDesc, x0)
    case(geoOptTypes%diis)
      call reset(self%pDiis, x0)
    case (geoOptTypes%lbfgs)
      call self%pLbfgs%reset(x0)
    end select

  end subroutine GeoOpt_reset


  !> Delivers the next point in the geometry optimization. When calling the first time, funciton
  !> value and gradient for the starting point of the minimization should be passed.
  subroutine GeoOpt_next(self, fx, dx, xNew, tConverged)

    !> Optimiser object
    type(Tgeoopt), intent(inout) :: self

    !> Function value for last point returned by this routine
    real(dp), intent(in) :: fx

    !> Gradient in the last point
    real(dp), intent(in) :: dx(:)

    !> New proposed point
    real(dp), intent(out) :: xNew(:)

    !> True, if gradient got below the specified tolerance.
    logical, intent(out) :: tConverged

    select case (self%iGeoOpt)
    case(geoOptTypes%conjugateGrad)
      call next(self%pConjGrad, fx, dx, xNew, tConverged)
    case (geoOptTypes%steepestDesc)
      call next(self%pSteepDesc, dx, xNew, tConverged)
    case (geoOptTypes%diis)
      call next(self%pDiis, dx, xNew, tConverged)
    case (geoOptTypes%lbfgs)
      call self%pLbfgs%next(fx, dx, xNew, tConverged)
    end select

  end subroutine GeoOpt_next

end module dftbp_geoopt
