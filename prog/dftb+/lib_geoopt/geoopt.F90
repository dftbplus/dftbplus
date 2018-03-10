!--------------------------------------------------------------------------------------------------!
!  DFTB+: general package for performing fast atomistic simulations                                !
!  Copyright (C) 2018  DFTB+ developers group                                                      !
!                                                                                                  !
!  See the LICENSE file for terms of usage and distribution.                                       !
!--------------------------------------------------------------------------------------------------!

!> General interface for the optimization algorithms
module geoopt
  use accuracy
  use conjgrad
  use steepdesc
  use gdiis
  use lbfgs
  implicit none

  private


  !> Interface type for the various geometry optimization algorithms
  type OGeoOpt
    private
    integer :: iGeoOpt
    type(OConjGrad), allocatable :: pConjGrad
    type(OSteepDesc), allocatable :: pSteepDesc
    type(ODiis), allocatable :: pDiis
    type(TLbfgs), allocatable :: pLbfgs
  end type OGeoOpt


  !> Creates a geometry optimizer
  interface init
    module procedure GeoOpt_initConjGrad
    module procedure GeoOpt_initSteepDesc
    module procedure GeoOpt_initDiis
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

  public :: OGeoOpt
  public :: init, reset, next

  !> Constants for the different geometry optimizers
  integer, parameter :: iConjGrad = 1
  integer, parameter :: iSteepDesc = 2
  integer, parameter :: iDiis = 3
  integer, parameter :: iLbfgs = 4

contains


  !> Creates a general geometry optimizier with a conjugate gradient instance
  subroutine GeoOpt_initConjGrad(self, pConjGrad)

    !> GeoOpt instance
    type(OGeoOpt), intent(out) :: self

    !> An already initialized conjugate gradient instance
    type(OConjGrad), allocatable, intent(inout) :: pConjGrad

    self%iGeoOpt = iConjGrad
    call move_alloc(pConjGrad, self%pConjGrad)

  end subroutine GeoOpt_initConjGrad


  !> Creates a general geometry optimizier with a steepest descent instance
  subroutine GeoOpt_initSteepDesc(self, pSteepDesc)

    !> GeoOpt instance
    type(OGeoOpt), intent(out) :: self

    !> An already initialized steepest descent instance
    type(OSteepDesc), allocatable, intent(inout) :: pSteepDesc

    self%iGeoOpt = iSteepDesc
    call move_alloc(pSteepDesc, self%pSteepDesc)

  end subroutine GeoOpt_initSteepDesc


  !> Creates a general geometry optimizier with a steepest descent instance
  subroutine GeoOpt_initDiis(self, pDiis)

    !> GeoOpt instance
    type(OGeoOpt), intent(out) :: self

    !> An already initialized modified DIIS instance
    type(ODiis), allocatable, intent(inout) :: pDiis

    self%iGeoOpt = iDiis
    call move_alloc(pDiis, self%pDiis)

  end subroutine GeoOpt_initDiis

  !> Creates a general geometry optimizier with a limited memory BFGS driver
  subroutine GeoOpt_initLbfgs(self, pLbfgs)

    !> GeoOpt instance
    type(OGeoOpt), intent(out) :: self

    !> An already initialized modified LBFGS
    type(Tlbfgs), allocatable, intent(inout) :: pLbfgs

    self%iGeoOpt = iLbfgs
    call move_alloc(pLbfgs, self%pLbfgs)

  end subroutine GeoOpt_initLbfgs

  !> Resets the geometry optimizer
  subroutine GeoOpt_reset(self, x0)

    !> GeoOpt instance
    type(OGeoOpt), intent(inout) :: self

    !> Initial coordinates
    real(dp), intent(in) :: x0(:)

    select case (self%iGeoOpt)
    case(iConjGrad)
      call reset(self%pConjGrad, x0)
    case(iSteepDesc)
      call reset(self%pSteepDesc, x0)
    case(iDiis)
      call reset(self%pDiis, x0)
    case (iLbfgs)
      call self%pLbfgs%reset(x0)
    end select

  end subroutine GeoOpt_reset


  !> Delivers the next point in the geometry optimization. When calling the first time, funciton
  !> value and gradient for the starting point of the minimization should be passed.
  subroutine GeoOpt_next(self, fx, dx, xNew, tConverged)

    !> Optimiser object
    type(OGeoOpt), intent(inout) :: self

    !> Function value for last point returned by this routine
    real(dp), intent(in) :: fx

    !> Gradient in the last point
    real(dp), intent(in) :: dx(:)

    !> New proposed point
    real(dp), intent(out) :: xNew(:)

    !> True, if gradient got below the specified tolerance.
    logical, intent(out) :: tConverged

    select case (self%iGeoOpt)
    case(iConjGrad)
      call next(self%pConjGrad, fx, dx, xNew, tConverged)
    case (iSteepDesc)
      call next(self%pSteepDesc, dx, xNew, tConverged)
    case (iDiis)
      call next(self%pDiis, dx, xNew, tConverged)
    case (iLbfgs)
      call self%pLbfgs%next(fx, dx, xNew, tConverged)
    end select

  end subroutine GeoOpt_next

end module geoopt
