!--------------------------------------------------------------------------------------------------!
!  DFTB+: general package for performing fast atomistic simulations                                !
!  Copyright (C) 2017  DFTB+ developers group                                                      !
!                                                                                                  !
!  See the LICENSE file for terms of usage and distribution.                                       !
!--------------------------------------------------------------------------------------------------!

!!* General purpose wrapper for MD integrators.
!!* @note Currently only velocity Verlet is wrapped.
module mdintegrator
#include "assert.h"
#include "allocate.h"  
  use Accuracy
  use VelocityVerlet
  !use VelocityOmelyan
  implicit none
  private

  public :: OMDIntegrator
  public :: create, destroy, next, rescale, state

  !!* Data for the MD integrator.
  type OMDIntegrator
    private
    integer :: integrator
    type(OVelocityVerlet), pointer :: pVelocityVerlet
    !type(OVelocityOmelyan), pointer :: pVelocityOmelyan
  end type OMDIntegrator

  interface create
    module procedure MDIntegrator_create_VVerlet
    !module procedure MDIntegrator_create_VOmelyan
  end interface

  interface destroy
    module procedure MDIntegrator_destroy
  end interface

  interface next
    module procedure MDIntegrator_next
  end interface
  
  interface rescale
    module procedure MDIntegrator_rescale
  end interface

  interface state
    module procedure MDIntegrator_state
  end interface
  
  !! Type of the integrator
  integer, parameter :: velocityVerlet_ = 1
  
contains

  !!* Create integrator wrapper for velocity Verlet.
  !!* @param self Integrator wrapper instance on exit.
  !!* @param pIntegrator Velocity Verlet integrator.
  subroutine MDIntegrator_create_VVerlet(self, pIntegrator)
    type(OMDIntegrator), pointer :: self
    type(OVelocityVerlet), pointer :: pIntegrator

    INITALLOCATE_P(self)
    self%integrator = velocityVerlet_
    self%pVelocityVerlet => pIntegrator
    
  end subroutine MDIntegrator_create_VVerlet
  
  !!* Destroys MD integrator.
  !!* @param self Integrator wrapper.
  subroutine MDIntegrator_destroy(self)
    type(OMDIntegrator), pointer :: self

    if (.not. associated(self)) then
      return
    end if
    select case (self%integrator)
    case (velocityVerlet_)
      call destroy(self%pVelocityVerlet)
    !case (velocityOmelyan_)
    !  call destroy(self%velocityOmelyan)
    end select
    DEALLOCATE_P(self)
    
  end subroutine MDIntegrator_destroy

  

  !!* Delivers the next velocities
  !!* @param self Integrator wrapper instance on exit.
  !!* @param accel Accelerations.
  !!* @param newCoords Updated coordinates.
  !!* @param newVelocity Updated velocities.
  subroutine MDIntegrator_next(self, accel, newCoord, newVelocity)
    type(OMDIntegrator), pointer :: self
    real(dp), intent(in) :: accel(:,:)
    real(dp), intent(out) :: newCoord(:,:)
    real(dp), intent(out) :: newVelocity(:,:)

    ASSERT(associated(self))
    
    select case (self%integrator)
    case (velocityVerlet_)
      call next(self%pVelocityVerlet, accel, newCoord, newVelocity)
    !case (velocityOmelyan_)
    !  call next(self%pVelocityOmelyan, accel, newCoord, newVelocity)
    end select
    
  end subroutine MDIntegrator_next
  
  subroutine MDIntegrator_rescale(self,coord,latVecs,stress)
    type(OMDIntegrator), pointer :: self
    real(dp),intent(inout)       :: coord(:,:)
    real(dp),intent(inout)       :: latVecs(3,3)
    real(dp),intent(in)          :: stress(3,3)
    
    call rescale(self%pVelocityVerlet,coord,latVecs,stress)
    
  end subroutine MDIntegrator_rescale

  !!* Probe internal state of the integrator
  subroutine MDIntegrator_state(self,fd)
    type(OMDIntegrator), pointer :: self
    integer,intent(in)           :: fd
    
    call state(self%pVelocityVerlet,fd)
    
  end subroutine MDIntegrator_state
  
end module mdintegrator
