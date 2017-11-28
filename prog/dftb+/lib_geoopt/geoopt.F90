!--------------------------------------------------------------------------------------------------!
!  DFTB+: general package for performing fast atomistic simulations                                !
!  Copyright (C) 2017  DFTB+ developers group                                                      !
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
    type(ODIIS), allocatable :: pDIIS
    type(TLBFGS), allocatable :: pLBFGS
  end type OGeoOpt


  !> Creates a geometry optimizer
  interface init
    module procedure GeoOpt_initConjGrad
    module procedure GeoOpt_initSteepDesc
    module procedure GeoOpt_initDIIS
    module procedure GeoOpt_initLBFGS
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
  integer, parameter :: iDIIS = 3
  integer, parameter :: iLBFGS = 4

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
  subroutine GeoOpt_initDIIS(self, pDIIS)

    !> GeoOpt instance
    type(OGeoOpt), intent(out) :: self

    !> An already initialized modified DIIS instance
    type(ODIIS), allocatable, intent(inout) :: pDIIS

    self%iGeoOpt = iDIIS
    call move_alloc(pDIIS, self%pDIIS)

  end subroutine GeoOpt_initDIIS

  !> Creates a general geometry optimizier with a limited memory BFGS driver
  subroutine GeoOpt_initLBFGS(self, pLBFGS)

    !> GeoOpt instance
    type(OGeoOpt), intent(out) :: self

    !> An already initialized modified LBFGS
    type(Tlbfgs), allocatable, intent(inout) :: pLBFGS

    self%iGeoOpt = iLBFGS
    call move_alloc(pLBFGS, self%pLBFGS)

  end subroutine GeoOpt_initLBFGS

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
    case(iDIIS)
      call reset(self%pDIIS, x0)
    case (iLBFGS)
      call self%pLBFGS%reset_lbfgs(x0)
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
    case (iDIIS)
      call next(self%pDIIS, dx, xNew, tConverged)
    case (iLBFGS)
      call self%pLBFGS%next_lbfgs(fx, dx, xNew, tConverged)
    end select

  end subroutine GeoOpt_next

end module geoopt
