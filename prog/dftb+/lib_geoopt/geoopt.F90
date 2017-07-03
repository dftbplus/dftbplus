!--------------------------------------------------------------------------------------------------!
!  DFTB+: general package for performing fast atomistic simulations                                !
!  Copyright (C) 2017  DFTB+ developers group                                                      !
!                                                                                                  !
!  See the LICENSE file for terms of usage and distribution.                                       !
!--------------------------------------------------------------------------------------------------!

!!* General interface for the optimization algorithms
module geoopt
  use accuracy
  use conjgrad
  use steepdesc
  use gdiis
  implicit none

  private

  !!* Interface type for the various geometry optimization algorithms
  type OGeoOpt
    private
    integer :: iGeoOpt
    type(OConjGrad), allocatable :: pConjGrad
    type(OSteepDesc), allocatable :: pSteepDesc
    type(ODIIS), allocatable :: pDIIS
  end type OGeoOpt

  !!* Creates a geometry optimizer
  interface init
    module procedure GeoOpt_initConjGrad
    module procedure GeoOpt_initSteepDesc
    module procedure GeoOpt_initDIIS
  end interface

  !!* Resets the optimizer
  interface reset
    module procedure GeoOpt_reset
  end interface

  !!* Delivers the next point in the minimization
  interface next
    module procedure GeoOpt_next
  end interface


  public :: OGeoOpt
  public :: init, reset, next

  !! Constanst for the different optimizers
  integer, parameter :: iConjGrad = 1
  integer, parameter :: iSteepDesc = 2
  integer, parameter :: iDIIS = 3

contains

  !!* Creates a general geometry optimizier with a conjugate gradient instance
  !!* @param self      GeoOpt instance
  !!* @param pConjGrad An already initialized conjugate gradient instance
  subroutine GeoOpt_initConjGrad(self, pConjGrad)
    type(OGeoOpt), intent(out) :: self
    type(OConjGrad), allocatable, intent(inout) :: pConjGrad

    self%iGeoOpt = iConjGrad
    call move_alloc(pConjGrad, self%pConjGrad)

  end subroutine GeoOpt_initConjGrad



  !!* Creates a general geometry optimizier with a steepest descent instance
  !!* @param self       GeoOpt instance
  !!* @param pSteepDesc An already initialized steepest descent instance
  subroutine GeoOpt_initSteepDesc(self, pSteepDesc)
    type(OGeoOpt), intent(out) :: self
    type(OSteepDesc), allocatable, intent(inout) :: pSteepDesc

    self%iGeoOpt = iSteepDesc
    call move_alloc(pSteepDesc, self%pSteepDesc)

  end subroutine GeoOpt_initSteepDesc

  !!* Creates a general geometry optimizier with a steepest descent instance
  !!* @param self       GeoOpt instance
  !!* @param pDIIS An already initialized modified DIIS instance
  subroutine GeoOpt_initDIIS(self, pDIIS)
    type(OGeoOpt), intent(out) :: self
    type(ODIIS), allocatable, intent(inout) :: pDIIS

    self%iGeoOpt = iDIIS
    call move_alloc(pDIIS, self%pDIIS)

  end subroutine GeoOpt_initDIIS


  !!* Resets the geometry optimizer
  !!* @param self GeoOpt instance
  !!* @param x0   Initial coordinates
  subroutine GeoOpt_reset(self, x0)
    type(OGeoOpt), intent(inout) :: self
    real(dp), intent(in) :: x0(:)

    select case (self%iGeoOpt)
    case(iConjGrad)
      call reset(self%pConjGrad, x0)
    case(iSteepDesc)
      call reset(self%pSteepDesc, x0)
    case(iDIIS)
      call reset(self%pDIIS, x0)
    end select

  end subroutine GeoOpt_reset


  !!* Delivers the next point in the geometry optimization
  !!* @param fx         Function value for last point returned by this routine
  !!* @param dx         Gradient in the last point
  !!* @param xNew       New proposed point
  !!* @param tConverged True, if gradient got below the specified tolerance.
  !!* @note When calling the first time, funciton value and gradient for the
  !!*   starting point of the minimization should be passed.
  subroutine GeoOpt_next(self, fx, dx, xNew, tConverged)
    type(OGeoOpt), intent(inout) :: self
    real(dp), intent(in) :: fx
    real(dp), intent(in) :: dx(:)
    real(dp), intent(out) :: xNew(:)
    logical, intent(out) :: tConverged

    select case (self%iGeoOpt)
    case(iConjGrad)
      call next(self%pConjGrad, fx, dx, xNew, tConverged)
    case (iSteepDesc)
      call next(self%pSteepDesc, dx, xNew, tConverged)
    case (iDIIS)
      call next(self%pDIIS, dx, xNew, tConverged)
    end select

  end subroutine GeoOpt_next


end module geoopt
