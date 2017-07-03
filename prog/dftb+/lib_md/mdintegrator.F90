!--------------------------------------------------------------------------------------------------!
!  DFTB+: general package for performing fast atomistic simulations                                !
!  Copyright (C) 2017  DFTB+ developers group                                                      !
!                                                                                                  !
!  See the LICENSE file for terms of usage and distribution.                                       !
!--------------------------------------------------------------------------------------------------!

#:include 'common.fypp'

!!* General purpose wrapper for MD integrators.
!!* @note Currently only velocity Verlet is wrapped.
module mdintegrator
  use assert
  use Accuracy
  use VelocityVerlet
  !use VelocityOmelyan
  implicit none
  private

  public :: OMDIntegrator
  public :: init, next, rescale, state

  !!* Data for the MD integrator.
  type OMDIntegrator
    private
    integer :: integrator
    type(OVelocityVerlet), allocatable :: pVelocityVerlet
    !type(OVelocityOmelyan), allocatable :: pVelocityOmelyan
  end type OMDIntegrator

  interface init
    module procedure MDIntegrator_init_VVerlet
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
  subroutine MDIntegrator_init_VVerlet(self, pIntegrator)
    type(OMDIntegrator), intent(out) :: self
    type(OVelocityVerlet), allocatable, intent(inout) :: pIntegrator

    self%integrator = velocityVerlet_
    call move_alloc(pIntegrator, self%pVelocityVerlet)

  end subroutine MDIntegrator_init_VVerlet


  !!* Delivers the next velocities
  !!* @param self Integrator wrapper instance on exit.
  !!* @param accel Accelerations.
  !!* @param newCoords Updated coordinates.
  !!* @param newVelocity Updated velocities.
  subroutine MDIntegrator_next(self, accel, newCoord, newVelocity)
    type(OMDIntegrator), intent(inout) :: self
    real(dp), intent(in) :: accel(:,:)
    real(dp), intent(out) :: newCoord(:,:)
    real(dp), intent(out) :: newVelocity(:,:)

    select case (self%integrator)
    case (velocityVerlet_)
      call next(self%pVelocityVerlet, accel, newCoord, newVelocity)
    end select

  end subroutine MDIntegrator_next


  subroutine MDIntegrator_rescale(self,coord,latVecs,stress)
    type(OMDIntegrator), intent(inout) :: self
    real(dp),intent(inout)       :: coord(:,:)
    real(dp),intent(inout)       :: latVecs(3,3)
    real(dp),intent(in)          :: stress(3,3)

    call rescale(self%pVelocityVerlet,coord,latVecs,stress)

  end subroutine MDIntegrator_rescale


  !!* Probe internal state of the integrator
  subroutine MDIntegrator_state(self,fd)
    type(OMDIntegrator), intent(in) :: self
    integer,intent(in)           :: fd

    call state(self%pVelocityVerlet,fd)

  end subroutine MDIntegrator_state

end module mdintegrator
