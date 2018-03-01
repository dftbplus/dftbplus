!--------------------------------------------------------------------------------------------------!
!  DFTB+: general package for performing fast atomistic simulations                                !
!  Copyright (C) 2018  DFTB+ developers group                                                      !
!                                                                                                  !
!  See the LICENSE file for terms of usage and distribution.                                       !
!--------------------------------------------------------------------------------------------------!

#:include 'common.fypp'

!> Andersen thermostat
!> Two versions of the Andersen thermostat are implemented, either the global
!> re-select or per atom reselect of velocities from the Maxwell-Boltzmann
!> distribution
!> See Andersen J. Chem. Phys. 72. 2384 (1980)
module andersentherm
  use assert
  use accuracy
  use mdcommon
  use ranlux
  use tempprofile
  implicit none
  private

  public :: OAndersenThermostat
  public :: init, getInitVelocities, updateVelocities, state


  !> Data for the Andersen thermostat
  type OAndersenThermostat
    private

    !> Nr. of atoms
    integer :: nAtom

    !> Random number generator
    type(ORanlux), allocatable :: pRanlux

    !> Mass of the atoms
    real(dp), allocatable :: mass(:)

    !> Temperature generator
    type(OTempProfile), pointer :: pTempProfile

    !> Rescale velocities individually?
    logical :: tRescaleIndiv

    !> Rescaling probability
    real(dp) :: wvScale

    !> MD framework
    type(OMDCommon) :: pMDFramework
  end type OAndersenThermostat


  !> Initialise thermostat object
  interface init
    module procedure AndersenThermostat_init
  end interface init


  !> Velocities at start of calculation
  interface getInitVelocities
    module procedure AndersenThermostat_getInitVelos
  end interface


  !> New atomic velocities
  interface updateVelocities
    module procedure AndersenThermostat_updateVelos
  end interface


  !> write state to disc
  interface state
    module procedure AndersenThermostat_state
  end interface

contains


  !> Creates an Andersen thermostat instance.
  subroutine AndersenThermostat_init(self, pRanlux, masses, tempProfile, &
      &rescaleIndiv, wvScale, pMDFramework)

    !> Initialised instance on exit.
    type(OAndersenThermostat), intent(out) :: self

    !> Random generator
    type(ORanlux), allocatable, intent(inout) :: pRanlux

    !> Masses of the atoms.
    real(dp), intent(in) :: masses(:)

    !> Pointer to a temperature profile object.
    type(OTempProfile), pointer, intent(in) :: tempProfile

    !> If velocities should be rescaled per atom
    logical, intent(in) :: rescaleIndiv

    !> Rescaling probability.
    real(dp), intent(in) :: wvScale

    !> Molecular dynamics general specifications
    type(OMDCommon), intent(in) :: pMDFramework

    call move_alloc(pRanlux, self%pRanlux)
    self%nAtom = size(masses)
    allocate(self%mass(self%nAtom))
    self%mass(:) = masses(:)
    self%pTempProfile => tempProfile
    self%tRescaleIndiv = rescaleIndiv
    self%wvScale = wvScale
    self%pMDFramework = pMDFramework

  end subroutine AndersenThermostat_init


  !> Returns the initial velocities.
  subroutine AndersenThermostat_getInitVelos(self, velocities)

    !> AndersenThermostat instance.
    type(OAndersenThermostat), intent(inout) :: self

    !> Contains the velocities on return.
    real(dp), intent(out) :: velocities(:,:)

    real(dp) :: kT
    integer :: ii

    @:ASSERT(all(shape(velocities) <= (/ 3, self%nAtom /)))

    call getTemperature(self%pTempProfile, kT)
    do ii = 1, self%nAtom
      call MaxwellBoltzmann(velocities(:,ii), self%mass(ii), kT, self%pRanlux)
    end do
    call restFrame(self%pMDFramework, velocities, self%mass)
    call rescaleTokT(self%pMDFramework, velocities, self%mass, kT)

  end subroutine AndersenThermostat_getInitVelos


  !> Updates the provided velocities according the current temperature.
  subroutine AndersenThermostat_updateVelos(self, velocities)

    !> AndersenThermostat instance.
    type(OAndersenThermostat), intent(inout) :: self

    !> Updated velocities on exit.
    real(dp), intent(inout) :: velocities(:,:)

    real(dp) :: rescaleChance
    real(dp) :: kT
    integer :: ii

    @:ASSERT(all(shape(velocities) <= (/ 3, self%nAtom /)))

    call getTemperature(self%pTempProfile, kT)
    if (self%tRescaleIndiv) then
      do ii = 1, self%nAtom
        call getRandom(self%pRanlux, rescaleChance)
        if (rescaleChance <= self%wvScale) then
          call MaxwellBoltzmann(velocities(:,ii), self%mass(ii), kT, &
              &self%pRanlux)
        end if
      end do
      call restFrame(self%pMDFramework, velocities, self%mass)
    else
      ! all atoms re-set at random
      call getRandom(self%pRanlux, rescaleChance)
      if (rescaleChance <= self%wvScale) then
        do ii = 1, self%nAtom
          call MaxwellBoltzmann(velocities(:,ii), self%mass(ii), kT, &
              &self%pRanlux)
        end do
        call restFrame(self%pMDFramework, velocities, self%mass)
        call rescaleTokT(self%pMDFramework, velocities, self%mass, kT)
      end if
    end if

  end subroutine AndersenThermostat_updateVelos


  !> Outputs internals of thermostat
  subroutine AndersenThermostat_state(self, fd)

    !> instance of thermostat
    type(OAndersenThermostat), intent(in) :: self

    !> filehandle to write out to
    integer,intent(in) :: fd

    ! no internal state, nothing to do

  end subroutine AndersenThermostat_state

end module andersentherm
