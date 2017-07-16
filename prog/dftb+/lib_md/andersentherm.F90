!--------------------------------------------------------------------------------------------------!
!  DFTB+: general package for performing fast atomistic simulations                                !
!  Copyright (C) 2017  DFTB+ developers group                                                      !
!                                                                                                  !
!  See the LICENSE file for terms of usage and distribution.                                       !
!--------------------------------------------------------------------------------------------------!

#:include 'common.fypp'

!!* Andersen thermostat
!!* Two versions of the Andersen thermostat are implemented, either the global
!!* re-select or per atom reselect of velocities from the Maxwell-Boltzmann
!!* distribution
!!* @ref Andersen J. Chem. Phys. 72. 2384 (1980)
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

  !!* Data for the Andersen thermostat
  type OAndersenThermostat
    private
    integer :: nAtom                    !* Nr. of atoms
    type(ORanlux), allocatable :: pRanlux   !* Random number generator
    real(dp), allocatable :: mass(:)        !* Mass of the atoms
    type(OTempProfile), pointer :: pTempProfile  !* Temperature generator
    logical :: tRescaleIndiv            !* Rescale velocities individually?
    real(dp) :: wvScale                 !* Rescaling probability
    type(OMDCommon) :: pMDFramework  !* MD framework
  end type OAndersenThermostat


  interface init
    module procedure AndersenThermostat_init
  end interface init

  interface getInitVelocities
    module procedure AndersenThermostat_getInitVelos
  end interface

  interface updateVelocities
    module procedure AndersenThermostat_updateVelos
  end interface

  interface state
    module procedure AndersenThermostat_state
  end interface

contains

  !!* Creates an Andersen thermostat instance.
  !!* @param self Initialised instance on exit.
  !!* @param pRanlux Pointer to the random generator.
  !!* @param masses Masses of the atoms.
  !!* @param tempProfile Pointer to a temperature profile object.
  !!* @param rescaleIndiv If velocities should be rescaled per atom
  !!* @param wvScale Rescaling probability.
  subroutine AndersenThermostat_init(self, pRanlux, masses, tempProfile, &
      &rescaleIndiv, wvScale, pMDFramework)
    type(OAndersenThermostat), intent(out) :: self
    type(ORanlux), allocatable, intent(inout) :: pRanlux
    real(dp), intent(in) :: masses(:)
    type(OTempProfile), pointer, intent(in) :: tempProfile
    logical, intent(in) :: rescaleIndiv
    real(dp), intent(in) :: wvScale
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


  !!* Returns the initial velocities.
  !!* @param self AndersenThermostat instance.
  !!* @param velocities Contains the velocities on return.
  subroutine AndersenThermostat_getInitVelos(self, velocities)
    type(OAndersenThermostat), intent(inout) :: self
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



  !!* Updates the provided velocities according the current temperature.
  !!* @param self AndersenThermostat instance.
  !!* @param velocities Updated velocities on exit.
  subroutine AndersenThermostat_updateVelos(self, velocities)
    type(OAndersenThermostat), intent(inout) :: self
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

  subroutine AndersenThermostat_state(self, fd)
    type(OAndersenThermostat), intent(in) :: self
    integer,intent(in)                 :: fd

    ! no internal state, nothing to do

  end subroutine AndersenThermostat_state

end module andersentherm
