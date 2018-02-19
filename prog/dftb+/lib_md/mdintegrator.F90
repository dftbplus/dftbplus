!--------------------------------------------------------------------------------------------------!
!  DFTB+: general package for performing fast atomistic simulations                                !
!  Copyright (C) 2018  DFTB+ developers group                                                      !
!                                                                                                  !
!  See the LICENSE file for terms of usage and distribution.                                       !
!--------------------------------------------------------------------------------------------------!

#:include 'common.fypp'

!> General purpose wrapper for MD integrators.
!>
!> Note: Currently only velocity Verlet is wrapped.
module mdintegrator
  use assert
  use Accuracy
  use VelocityVerlet
  implicit none
  private

  public :: OMDIntegrator
  public :: init, next, rescale, state


  !> Data for the MD integrator.
  type OMDIntegrator
    private

    !> Integrator type
    integer :: integrator

    !> Verlet case
    type(OVelocityVerlet), allocatable :: pVelocityVerlet
  end type OMDIntegrator


  !> Initialise integrator
  interface init
    module procedure MDIntegrator_init_VVerlet
  end interface init


  !> Take a geometry step
  interface next
    module procedure MDIntegrator_next
  end interface next


  !> Barostat rescale if required
  interface rescale
    module procedure MDIntegrator_rescale
  end interface rescale


  !> Output state of the integrator
  interface state
    module procedure MDIntegrator_state
  end interface state


  !> Type of the integrator
  integer, parameter :: velocityVerlet_ = 1

contains


  !> Create integrator wrapper for velocity Verlet.
  subroutine MDIntegrator_init_VVerlet(self, pIntegrator)

    !> Integrator wrapper instance on exit.
    type(OMDIntegrator), intent(out) :: self

    !> Velocity Verlet integrator.
    type(OVelocityVerlet), allocatable, intent(inout) :: pIntegrator

    self%integrator = velocityVerlet_
    call move_alloc(pIntegrator, self%pVelocityVerlet)

  end subroutine MDIntegrator_init_VVerlet


  !> Delivers the next velocities
  subroutine MDIntegrator_next(self, accel, newCoord, newVelocity)

    !> Integrator wrapper instance on exit.
    type(OMDIntegrator), intent(inout) :: self

    !> Accelerations.
    real(dp), intent(in) :: accel(:,:)

    !> Updated coordinates.
    real(dp), intent(out) :: newCoord(:,:)

    real(dp), intent(out) :: newVelocity(:,:)

    select case (self%integrator)
    case (velocityVerlet_)
      call next(self%pVelocityVerlet, accel, newCoord, newVelocity)
    end select

  end subroutine MDIntegrator_next


  !> Apply Barostat type rescales if needed
  !>
  !> Note: Should be packaged in the same way as thermostats
  subroutine MDIntegrator_rescale(self,coord,latVecs,stress)

    !> Integrator instance
    type(OMDIntegrator), intent(inout) :: self

    !> coordinates of atoms
    real(dp),intent(inout) :: coord(:,:)

    !> lattice vectors
    real(dp),intent(inout) :: latVecs(3,3)

    !> stress tensor
    real(dp),intent(in) :: stress(3,3)

    call rescale(self%pVelocityVerlet,coord,latVecs,stress)

  end subroutine MDIntegrator_rescale


  !> Probe internal state of the integrator, writing this to disc
  subroutine MDIntegrator_state(self,fd)

    !> Integrator instance
    type(OMDIntegrator), intent(in) :: self

    !> file handle to write to
    integer,intent(in) :: fd

    call state(self%pVelocityVerlet,fd)

  end subroutine MDIntegrator_state

end module mdintegrator
