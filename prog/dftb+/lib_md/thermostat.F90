!--------------------------------------------------------------------------------------------------!
!  DFTB+: general package for performing fast atomistic simulations                                !
!  Copyright (C) 2006 - 2020  DFTB+ developers group                                               !
!                                                                                                  !
!  See the LICENSE file for terms of usage and distribution.                                       !
!--------------------------------------------------------------------------------------------------!

!> Contains wrapper for all thermostats.
module dftbp_thermostat
  use dftbp_accuracy
  use dftbp_dummytherm
  use dftbp_andersentherm
  use dftbp_berendsentherm
  use dftbp_nhctherm
  implicit none
  private

  public :: TThermostat
  public :: init, getInitVelocities, updateVelocities, state


  !> Data for the thermostat wrapper.
  type TThermostat
    private

    !> Thermostat type
    integer :: thermostat

    !> Dummy for no temperature control
    type(TDummyThermostat), allocatable :: pDummy

    !> Anderson rescaling
    type(TAndersenThermostat), allocatable :: pAndersen

    !> Berendsen stochastic
    type(TBerendsenThermostat), allocatable :: pBerendsen

    !> Nose-Hoover and chains
    type(TNHCThermostat), allocatable :: pNHC
  end type TThermostat


  !> Initialise thermostat in use
  interface init
    module procedure Thermostat_init_Dummy
    module procedure Thermostat_init_Andersen
    module procedure Thermostat_init_Berendsen
    module procedure Thermostat_init_NHC
  end interface


  !> Initial atomic velocities
  interface getInitVelocities
    module procedure Thermostat_getInitVelocities
  end interface


  !> Update velocities, subject to thermostat in use
  interface updateVelocities
    module procedure Thermostat_updateVelocities
  end interface


  !> Return the state of the thermostat
  interface state
    module procedure Thermostat_state
  end interface


  !> Thermostat types
  integer, parameter :: dummy_ = 0
  integer, parameter :: andersen_ = 1
  integer, parameter :: berendsen_ = 2
  integer, parameter :: nhc_ = 3

contains


  !> Creates a thermostat wrapper for a DummyThermostat.
  subroutine Thermostat_init_Dummy(self, pThermostat)

    !> Wrapper instance on exit.
    type(TThermostat), intent(out) :: self

    !> A DummyThermostat.
    type(TDummyThermostat), allocatable, intent(inout) :: pThermostat

    self%thermostat = dummy_
    call move_alloc(pThermostat, self%pDummy)

  end subroutine Thermostat_init_Dummy


  !> Creates a thermostat wrapper for an AndersenThermostat.
  subroutine Thermostat_init_Andersen(self, pThermostat)

    !> Wrapper instance on exit.
    type(TThermostat), intent(out) :: self

    !> An Andersen Thermostat.
    type(TAndersenThermostat), allocatable, intent(inout) :: pThermostat

    self%thermostat = andersen_
    call move_alloc(pThermostat, self%pAndersen)

  end subroutine Thermostat_init_Andersen


  !> Creates a thermostat wrapper for a BerendsenThermostat.
  subroutine Thermostat_init_Berendsen(self, pThermostat)

    !> Wrapper instance on exit.
    type(TThermostat), intent(out) :: self

    !> A Berendsen Thermostat.
    type(TBerendsenThermostat), allocatable, intent(inout) :: pThermostat

    self%thermostat = berendsen_
    call move_alloc(pThermostat, self%pBerendsen)

  end subroutine Thermostat_init_Berendsen


  !> Creates a thermostat wrapper for a NHCThermostat.
  subroutine Thermostat_init_NHC(self, pThermostat)

    !> Wrapper instance on exit.
    type(TThermostat), intent(out) :: self

    !> A NHC Thermostat.
    type(TNHCThermostat), allocatable, intent(inout) :: pThermostat

    self%thermostat = nhc_
    call move_alloc(pThermostat, self%pNHC)

  end subroutine Thermostat_init_NHC


  !> Returns the initial velocities
  subroutine Thermostat_getInitVelocities(self, velocities)

    !> Wrapper instance.
    type(TThermostat), intent(inout) :: self

    !> Velocities on exit.
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


  !> Updates the velocities.
  !>

  !> Note: The DummyThermostat has no method to update the velocities, so the wrapper returns
  !> without touching the velocities.
  subroutine Thermostat_updateVelocities(self, velocities)

    !> Wrapper instance.
    type(TThermostat), intent(inout) :: self

    !> Updated velocities on exit.
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


  !> Probe internal state of the thermostat
  subroutine Thermostat_state(self, fd)

    !> Wrapper instance.
    type(TThermostat), intent(in) :: self

    !> file handle to write state out to
    integer, intent(in) :: fd

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

end module dftbp_thermostat
