!--------------------------------------------------------------------------------------------------!
!  DFTB+: general package for performing fast atomistic simulations                                !
!  Copyright (C) 2017  DFTB+ developers group                                                      !
!                                                                                                  !
!  See the LICENSE file for terms of usage and distribution.                                       !
!--------------------------------------------------------------------------------------------------!

!!* General interface for the optimization algorithms
module geoopt
#include "allocate.h"
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
    type(OConjGrad), pointer :: pConjGrad
    type(OSteepDesc), pointer :: pSteepDesc
    type(ODIIS), pointer :: pDIIS
  end type OGeoOpt

  !!* Creates a geometry optimizer
  interface create
    module procedure GeoOpt_createConjGrad
    module procedure GeoOpt_createSteepDesc
    module procedure GeoOpt_createDIIS
  end interface

  !!* Destroys the optimizer
  interface destroy
    module procedure GeoOpt_destroy
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
  public :: create, destroy, reset, next

  !! Constanst for the different optimizers
  integer, parameter :: iConjGrad = 1
  integer, parameter :: iSteepDesc = 2
  integer, parameter :: iDIIS = 3
  
contains

  !!* Allocates the object and nulls the pointers in it
  !!* @param self GeoOpt instance
  subroutine GeoOpt_createCommon(self)
    type(OGeoOpt), pointer :: self

    INITALLOCATE_P(self)
    self%pConjGrad => null()
    self%pSteepDesc => null()
    self%pDIIS => null()

  end subroutine GeoOpt_createCommon

  

  !!* Creates a general geometry optimizier with a conjugate gradient instance
  !!* @param self      GeoOpt instance
  !!* @param pConjGrad An already initialized conjugate gradient instance
  subroutine GeoOpt_createConjGrad(self, pConjGrad)
    type(OGeoOpt), pointer :: self
    type(OConjGrad), pointer :: pConjGrad

    call GeoOpt_createCommon(self)
    self%iGeoOpt = iConjGrad
    self%pConjGrad => pConjGrad

  end subroutine GeoOpt_createConjGrad

  

  !!* Creates a general geometry optimizier with a steepest descent instance
  !!* @param self       GeoOpt instance
  !!* @param pSteepDesc An already initialized steepest descent instance
  subroutine GeoOpt_createSteepDesc(self, pSteepDesc)
    type(OGeoOpt), pointer :: self
    type(OSteepDesc), pointer :: pSteepDesc

    call GeoOpt_createCommon(self)
    self%iGeoOpt = iSteepDesc
    self%pSteepDesc => pSteepDesc

  end subroutine GeoOpt_createSteepDesc

  !!* Creates a general geometry optimizier with a steepest descent instance
  !!* @param self       GeoOpt instance
  !!* @param pDIIS An already initialized modified DIIS instance
  subroutine GeoOpt_createDIIS(self, pDIIS)
    type(OGeoOpt), pointer :: self
    type(ODIIS), pointer :: pDIIS

    call GeoOpt_createCommon(self)
    self%iGeoOpt = iDIIS
    self%pDIIS => pDIIS

  end subroutine GeoOpt_createDIIS

  
  

  !!* Resets the geometry optimizer
  !!* @param self GeoOpt instance
  !!* @param x0   Initial coordinates
  subroutine GeoOpt_reset(self, x0)
    type(OGeoOpt), pointer :: self
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

  

  !!* Destroys the geometry optimizer
  !!* @param self GeoOpt instance
  subroutine GeoOpt_destroy(self)
    type(OGeoOpt), pointer :: self

    if (associated(self)) then
      select case (self%iGeoOpt)
      case(iConjGrad)
        call destroy(self%pConjGrad)
      case(iSteepDesc)
        call destroy(self%pSteepDesc)
      case(iDIIS)
        call destroy(self%pDIIS)
      end select
      DEALLOCATE_P(self)
    end if
  end subroutine GeoOpt_destroy



  !!* Delivers the next point in the geometry optimization
  !!* @param fx         Function value for last point returned by this routine
  !!* @param dx         Gradient in the last point
  !!* @param xNew       New proposed point
  !!* @param tConverged True, if gradient got below the specified tolerance.
  !!* @note When calling the first time, funciton value and gradient for the
  !!*   starting point of the minimization should be passed.
  subroutine GeoOpt_next(self, fx, dx, xNew, tConverged)
    type(OGeoOpt), pointer :: self
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
