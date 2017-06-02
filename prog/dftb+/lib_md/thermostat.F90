!--------------------------------------------------------------------------------------------------!
!  DFTB+: general package for performing fast atomistic simulations                                !
!  Copyright (C) 2017  DFTB+ developers group                                                      !
!                                                                                                  !
!  See the LICENSE file for terms of usage and distribution.                                       !
!--------------------------------------------------------------------------------------------------!

!!* Contains wrapper for all thermostats.
module thermostat
#include "allocate.h"  
  use accuracy
  use dummytherm
  use andersentherm
  use berendsentherm
  use nhctherm
  implicit none
  private

  public :: OThermostat
  public :: create, destroy, getInitVelocities, updateVelocities, state

  !!* Data for the termostat wrapper.
  type OThermostat
    private
    integer :: thermostat
    type(ODummyThermostat), pointer :: pDummy
    type(OAndersenThermostat), pointer :: pAndersen
    type(OBerendsenThermostat), pointer :: pBerendsen
    type(ONHCThermostat), pointer :: pNHC
  end type OThermostat

  interface create
    module procedure Thermostat_create_Dummy
    module procedure Thermostat_create_Andersen
    module procedure Thermostat_create_Berendsen
    module procedure Thermostat_create_NHC
  end interface

  interface destroy
    module procedure Thermostat_destroy
  end interface

  interface getInitVelocities
    module procedure Thermostat_getInitVelocities
  end interface

  interface updateVelocities
    module procedure Thermostat_updateVelocities
  end interface
  
  interface state
    module procedure Thermostat_state
  end interface
  
  !! Thermostat types
  integer, parameter :: dummy_ = 0
  integer, parameter :: andersen_ = 1
  integer, parameter :: berendsen_ = 2
  integer, parameter :: nhc_ = 3

contains


  !!* Creates a thermostat wrapper for a DummyThermostat.
  !!* @param self Wrapper instance on exit.
  !!* @param pThermostat Pointer to a DummyThermostat.
  subroutine Thermostat_create_Dummy(self, pThermostat)
    type(OThermostat), pointer :: self
    type(ODummyThermostat), pointer :: pThermostat

    INITALLOCATE_P(self)
    self%thermostat = dummy_
    self%pDummy => pThermostat
    
  end subroutine Thermostat_create_Dummy
  

  
  !!* Creates a thermostat wrapper for an AndersenThermostat.
  !!* @param self Wrapper instance on exit.
  !!* @param pThermostat Pointer to a AndersenThermostat.
  subroutine Thermostat_create_Andersen(self, pThermostat)
    type(OThermostat), pointer :: self
    type(OAndersenThermostat), pointer :: pThermostat

    INITALLOCATE_P(self)
    self%thermostat = andersen_
    self%pAndersen => pThermostat
    
  end subroutine Thermostat_create_Andersen
    
  !!* Creates a thermostat wrapper for a BerendsenThermostat.
  !!* @param self Wrapper instance on exit.
  !!* @param pThermostat Pointer to a BerendsenThermostat.
  subroutine Thermostat_create_Berendsen(self, pThermostat)
    type(OThermostat), pointer :: self
    type(OBerendsenThermostat), pointer :: pThermostat

    INITALLOCATE_P(self)
    self%thermostat = berendsen_
    self%pBerendsen => pThermostat
    
  end subroutine Thermostat_create_Berendsen
  
  !!* Creates a thermostat wrapper for a NHCThermostat.
  !!* @param self Wrapper instance on exit.
  !!* @param pThermostat Pointer to a NHCThermostat.
  subroutine Thermostat_create_NHC(self, pThermostat)
    type(OThermostat), pointer :: self
    type(ONHCThermostat), pointer :: pThermostat

    INITALLOCATE_P(self)
    self%thermostat = nhc_
    self%pNHC => pThermostat
    
  end subroutine Thermostat_create_NHC
  
  
  !!* Destroys the thermostat wrapper.
  !!* @param self Wrapper instance.
  subroutine Thermostat_destroy(self)
    type(OThermostat), pointer :: self

    if (.not. associated(self)) then
      return
    end if
    select case (self%thermostat)
    case (dummy_)
      call destroy(self%pDummy)
    case (andersen_)
      call destroy(self%pAndersen)
    case (berendsen_)
      call destroy(self%pBerendsen) 
    case (nhc_)
      call destroy(self%pNHC) 
    end select
    DEALLOCATE_P(self)
    
  end subroutine Thermostat_destroy
  
  !!* Returns the initial velocities
  !!* @param self Wrapper instance.
  !!* @param velocities Velocities on exit.
  subroutine Thermostat_getInitVelocities(self, velocities)
    type(OThermostat), pointer :: self
    real(dp), intent(out) :: velocities(:,:)

    select case (self%thermostat)
    case (dummy_)
      call getInitVelocities(self%pDummy, velocities)
    case(andersen_)
      call getInitVelocities(self%pAndersen, velocities)
    case(berendsen_)
      call getInitVelocities(self%pBerendsen, velocities)  
    case(nhc_)
      call getInitVelocities(self%pNHC, velocities)  
    end select

  end subroutine Thermostat_getInitVelocities

  !!* Updates the velocities.
  !!* @param self Wrapper instance.
  !!* @param velocities Updated velocities on exit.
  !!* @note The DummyThermostat has no method to update the velocities,
  !!*   so the wrapper returns without touching the velocities.
  subroutine Thermostat_updateVelocities(self, velocities)
    type(OThermostat), pointer :: self
    real(dp), intent(inout) :: velocities(:,:)

    select case (self%thermostat)
    case (dummy_)
      continue
    case(andersen_)
      call updateVelocities(self%pAndersen, velocities)
    case(berendsen_)
      call updateVelocities(self%pBerendsen, velocities)
    case(nhc_)
      call updateVelocities(self%pNHC, velocities)
    end select

  end subroutine Thermostat_updateVelocities

  !!* Probe internal state of the thermostat
  !!* @param self Wrapper instance.
  !!* @param fd file handle to write state too
  subroutine Thermostat_state(self, fd)
    type(OThermostat), pointer :: self
    integer, intent(in)        :: fd

    select case (self%thermostat)
    case (dummy_)
      continue
    case(andersen_)
      call state(self%pAndersen, fd)
    case(berendsen_)
      call state(self%pBerendsen, fd)
    case(nhc_)
      call state(self%pNHC, fd)
    end select

  end subroutine Thermostat_state


  
end module thermostat
