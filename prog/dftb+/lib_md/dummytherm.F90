!--------------------------------------------------------------------------------------------------!
!  DFTB+: general package for performing fast atomistic simulations                                !
!  Copyright (C) 2017  DFTB+ developers group                                                      !
!                                                                                                  !
!  See the LICENSE file for terms of usage and distribution.                                       !
!--------------------------------------------------------------------------------------------------!

#:include 'common.fypp'

!!* Dummy thermostat, delivers only initial velocities according to the
!!* Maxwell-Boltzmann statistics.
module dummytherm
  use assert
  use accuracy
  use mdcommon
  use ranlux
  implicit none
  private

  public :: ODummyThermostat
  public :: init, getInitVelocities, state

  !!* Data for dummy thermostat
  type ODummyThermostat
    private
    integer :: nAtom                  !* Nr. of atoms
    real(dp) :: kT                    !* Temperature
    real(dp), allocatable :: mass(:)      !* Mass of the atoms
    type(ORanlux), allocatable :: pRanlux !* Random number generator.
    type(OMDCommon) :: pMDFrame !* MD Framwork
  end type ODummyThermostat

  interface init
    module procedure DummyThermostat_init
  end interface

  interface getInitVelocities
    module procedure DummyThermostat_getInitVelos
  end interface

  interface state
    module procedure DummyThermostat_state
  end interface

contains

  !!* Creates a DummyThermostat instance.
  !!* @param self Initialised DummyThermostat instance on return.
  !!* @param kT Temperature of the thermostat
  !!* @param pRanlux Random generator
  subroutine DummyThermostat_init(self, kT, mass, pRanlux, pMDFrame)
    type(ODummyThermostat), intent(out) :: self
    real(dp), intent(in) :: kT
    real(dp), intent(in) :: mass(:)
    type(ORanlux), allocatable, intent(inout) :: pRanlux
    type(OMDCommon), intent(in) :: pMDFrame

    self%kT = kT
    self%nAtom = size(mass)
    allocate(self%mass(self%nAtom))
    self%mass = mass(:)
    call move_alloc(pRanlux, self%pRanlux)
    self%pMDFrame = pMDFrame

  end subroutine DummyThermostat_init


  !!* Returns the initial velocities.
  !!* @param self Thermostat instance.
  !!* @param velocities Contains the velocities on return.
  subroutine DummyThermostat_getInitVelos(self, velocities)
    type(ODummyThermostat), intent(inout) :: self
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

  subroutine DummyThermostat_state(self, fd)
    type(ODummyThermostat), intent(in) :: self
    integer,intent(in)                  :: fd

    ! no internal state, nothing to do

  end subroutine DummyThermostat_state

end module dummytherm
