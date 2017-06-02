!--------------------------------------------------------------------------------------------------!
!  DFTB+: general package for performing fast atomistic simulations                                !
!  Copyright (C) 2017  DFTB+ developers group                                                      !
!                                                                                                  !
!  See the LICENSE file for terms of usage and distribution.                                       !
!--------------------------------------------------------------------------------------------------!

!!* Dummy thermostat, delivers only initial velocities according to the
!!* Maxwell-Boltzmann statistics.
module dummytherm
#include "assert.h"
#include "allocate.h"  
  use accuracy
  use mdcommon
  use ranlux
  implicit none
  private

  public :: ODummyThermostat
  public :: create, destroy, getInitVelocities, state

  !!* Data for dummy thermostat
  type ODummyThermostat
    private
    integer :: nAtom                  !* Nr. of atoms
    real(dp) :: kT                    !* Temperature
    real(dp), pointer :: mass(:)      !* Mass of the atoms
    type(ORanlux), pointer :: pRanlux !* Random number generator.
    type(OMDCommon), pointer :: pMDFrame !* MD Framwork
  end type ODummyThermostat

  interface create
    module procedure DummyThermostat_create
  end interface

  interface destroy
    module procedure DummyThermostat_destroy
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
  subroutine DummyThermostat_create(self, kT, mass, pRanlux, pMDFrame)
    type(ODummyThermostat), pointer :: self
    real(dp), intent(in) :: kT
    real(dp), intent(in) :: mass(:)
    type(ORanlux), pointer :: pRanlux
    type(OMDCommon), pointer :: pMDFrame

    INITALLOCATE_P(self)
    self%kT = kT
    self%nAtom = size(mass)
    INITALLOCATE_PARR(self%mass, (self%nAtom))
    self%mass = mass(:)
    self%pRanlux => pRanlux
    self%pMDFrame => pMDFrame
    
  end subroutine DummyThermostat_create

  !!* Destroys a DummyThermostat instance.
  !!* @param self DummyThermostat instance.
  subroutine DummyThermostat_destroy(self)
    type(ODummyThermostat), pointer :: self

    if (.not. associated(self)) then
      return
    end if
    DEALLOCATE_PARR(self%mass)
    DEALLOCATE_P(self)
    
  end subroutine DummyThermostat_destroy

  !!* Returns the initial velocities.
  !!* @param self Thermostat instance.
  !!* @param velocities Contains the velocities on return.
  subroutine DummyThermostat_getInitVelos(self, velocities)
    type(ODummyThermostat), pointer :: self
    real(dp), intent(out) :: velocities(:,:)

    integer :: ii

    ASSERT(associated(self))
    ASSERT(all(shape(velocities) >= (/ 3, self%nAtom /)))

    do ii = 1, self%nAtom
      call MaxwellBoltzmann(velocities(:,ii), self%mass(ii), self%kT, &
          & self%pRanlux)
    end do
    call restFrame(self%pMDFrame, velocities(:,:), self%mass)
    call rescaleTokT(self%pMDFrame, velocities(:,:), self%mass, self%kT)

  end subroutine DummyThermostat_getInitVelos
  
  subroutine DummyThermostat_state(self, fd)
    type(ODummyThermostat), pointer :: self
    integer,intent(in)                  :: fd
    
    ! no internal state, nothing to do
    
  end subroutine DummyThermostat_state
  
end module dummytherm
