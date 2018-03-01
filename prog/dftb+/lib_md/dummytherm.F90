!--------------------------------------------------------------------------------------------------!
!  DFTB+: general package for performing fast atomistic simulations                                !
!  Copyright (C) 2018  DFTB+ developers group                                                      !
!                                                                                                  !
!  See the LICENSE file for terms of usage and distribution.                                       !
!--------------------------------------------------------------------------------------------------!

#:include 'common.fypp'

!> Dummy thermostat, delivers only initial velocities according to the Maxwell-Boltzmann statistics.
module dummytherm
  use assert
  use accuracy
  use mdcommon
  use ranlux
  implicit none
  private

  public :: ODummyThermostat
  public :: init, getInitVelocities, state


  !> Data for dummy thermostat
  type ODummyThermostat
    private

    !> Nr. of atoms
    integer :: nAtom

    !> Temperature
    real(dp) :: kT

    !> Mass of the atoms
    real(dp), allocatable :: mass(:)

    !> Random number generator.
    type(ORanlux), allocatable :: pRanlux

    !> MD Framwork
    type(OMDCommon) :: pMDFrame
  end type ODummyThermostat


  !> Initialise thermostat object
  interface init
    module procedure DummyThermostat_init
  end interface


  !> Velocities at start of calculation
  interface getInitVelocities
    module procedure DummyThermostat_getInitVelos
  end interface


  !> write state to disc
  interface state
    module procedure DummyThermostat_state
  end interface

contains


  !> Creates a DummyThermostat instance.
  subroutine DummyThermostat_init(self, kT, mass, pRanlux, pMDFrame)
    type(ODummyThermostat), intent(out) :: self

    !> Initialised DummyThermostat instance on return.
    real(dp), intent(in) :: kT

    !> Temperature of the thermostat
    real(dp), intent(in) :: mass(:)

    !> Random generator
    type(ORanlux), allocatable, intent(inout) :: pRanlux

    !> thermostat object
    type(OMDCommon), intent(in) :: pMDFrame

    self%kT = kT
    self%nAtom = size(mass)
    allocate(self%mass(self%nAtom))
    self%mass = mass(:)
    call move_alloc(pRanlux, self%pRanlux)
    self%pMDFrame = pMDFrame

  end subroutine DummyThermostat_init


  !> Returns the initial velocities.
  subroutine DummyThermostat_getInitVelos(self, velocities)

    !> Thermostat instance.
    type(ODummyThermostat), intent(inout) :: self

    !> Contains the velocities on return.
    real(dp), intent(out) :: velocities(:,:)

    integer :: ii

    @:ASSERT(all(shape(velocities) >= (/ 3, self%nAtom /)))

    do ii = 1, self%nAtom
      call MaxwellBoltzmann(velocities(:,ii), self%mass(ii), self%kT, &
          & self%pRanlux)
    end do
    call restFrame(self%pMDFrame, velocities(:,:), self%mass)
    call rescaleTokT(self%pMDFrame, velocities(:,:), self%mass, self%kT)

  end subroutine DummyThermostat_getInitVelos


  !> no internal state, nothing to do
  subroutine DummyThermostat_state(self, fd)

    !> thermostat object
    type(ODummyThermostat), intent(in) :: self

    !> file unit
    integer,intent(in) :: fd

  end subroutine DummyThermostat_state

end module dummytherm
