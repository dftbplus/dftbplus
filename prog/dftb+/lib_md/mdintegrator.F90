!--------------------------------------------------------------------------------------------------!
!  DFTB+: general package for performing fast atomistic simulations                                !
!  Copyright (C) 2006 - 2020  DFTB+ developers group                                               !
!                                                                                                  !
!  See the LICENSE file for terms of usage and distribution.                                       !
!--------------------------------------------------------------------------------------------------!

#:include 'common.fypp'

!> General purpose wrapper for MD integrators.
!>
!> Note: Currently only velocity Verlet is wrapped.
module dftbp_mdintegrator
  use dftbp_assert
  use dftbp_accuracy
  use dftbp_velocityverlet
  implicit none
  private

  public :: TMDIntegrator
  public :: init, next, rescale, state


  !> Data for the MD integrator.
  type TMDIntegrator
    private

    !> Integrator type
    integer :: integrator

    !> Verlet case
    type(TVelocityVerlet), allocatable :: pVelocityVerlet
  end type TMDIntegrator


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
    type(TMDIntegrator), intent(out) :: self

    !> Velocity Verlet integrator.
    type(TVelocityVerlet), allocatable, intent(inout) :: pIntegrator

    self%integrator = velocityVerlet_
    call move_alloc(pIntegrator, self%pVelocityVerlet)

  end subroutine MDIntegrator_init_VVerlet


  !> Delivers the next velocities
  subroutine MDIntegrator_next(self, accel, newCoord, newVelocity)

    !> Integrator wrapper instance on exit.
    type(TMDIntegrator), intent(inout) :: self

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
    type(TMDIntegrator), intent(inout) :: self

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
    type(TMDIntegrator), intent(in) :: self

    !> file handle to write to
    integer,intent(in) :: fd

    call state(self%pVelocityVerlet,fd)

  end subroutine MDIntegrator_state

end module dftbp_mdintegrator
