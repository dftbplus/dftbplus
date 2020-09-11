!--------------------------------------------------------------------------------------------------!
!  DFTB+: general package for performing fast atomistic simulations                                !
!  Copyright (C) 2006 - 2020  DFTB+ developers group                                               !
!                                                                                                  !
!  See the LICENSE file for terms of usage and distribution.                                       !
!--------------------------------------------------------------------------------------------------!

#:include 'common.fypp'

!> Andersen thermostat
!> Two versions of the Andersen thermostat are implemented, either the global
!> re-select or per atom reselect of velocities from the Maxwell-Boltzmann
!> distribution
!> See Andersen J. Chem. Phys. 72. 2384 (1980)
module dftbp_andersentherm
  use dftbp_assert
  use dftbp_accuracy
  use dftbp_mdcommon
  use dftbp_ranlux
  use dftbp_tempprofile
  use dftbp_message
  implicit none
  private

  public :: TAndersenThermostat
  public :: init, getInitVelocities, updateVelocities, state


  !> Data for the Andersen thermostat
  type TAndersenThermostat
    private

    !> Nr. of atoms
    integer :: nAtom

    !> Random number generator
    type(TRanlux), allocatable :: pRanlux

    !> Mass of the atoms
    real(dp), allocatable :: mass(:)

    !> Temperature generator
    type(TTempProfile), pointer :: pTempProfile

    !> Rescale velocities individually?
    logical :: tRescaleIndiv

    !> Rescaling probability
    real(dp) :: wvScale

    !> MD framework
    type(TMDCommon) :: pMDFramework

  end type TAndersenThermostat


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
  subroutine AndersenThermostat_init(this, pRanlux, masses, tempProfile, &
      &rescaleIndiv, wvScale, pMDFramework)

    !> Initialised instance on exit.
    type(TAndersenThermostat), intent(out) :: this

    !> Random generator
    type(TRanlux), allocatable, intent(inout) :: pRanlux

    !> Masses of the atoms.
    real(dp), intent(in) :: masses(:)

    !> Pointer to a temperature profile object.
    type(TTempProfile), pointer, intent(in) :: tempProfile

    !> If velocities should be rescaled per atom
    logical, intent(in) :: rescaleIndiv

    !> Rescaling probability.
    real(dp), intent(in) :: wvScale

    !> Molecular dynamics general specifications
    type(TMDCommon), intent(in) :: pMDFramework

    call move_alloc(pRanlux, this%pRanlux)
    this%nAtom = size(masses)
    allocate(this%mass(this%nAtom))
    this%mass(:) = masses(:)
    this%pTempProfile => tempProfile
    this%tRescaleIndiv = rescaleIndiv
    this%wvScale = wvScale
    this%pMDFramework = pMDFramework

  end subroutine AndersenThermostat_init


  !> Returns the initial velocities.
  subroutine AndersenThermostat_getInitVelos(this, velocities)

    !> AndersenThermostat instance.
    type(TAndersenThermostat), intent(inout) :: this

    !> Contains the velocities on return.
    real(dp), intent(out) :: velocities(:,:)

    real(dp) :: kT
    integer :: ii

    @:ASSERT(all(shape(velocities) <= (/ 3, this%nAtom /)))

    call this%pTempProfile%getTemperature(kT)
    if (kT < minTemp) then
      call error("Andersen thermostat not supported at zero temperature")
    end if
    do ii = 1, this%nAtom
      call MaxwellBoltzmann(velocities(:,ii), this%mass(ii), kT, this%pRanlux)
    end do
    call restFrame(this%pMDFramework, velocities, this%mass)
    call rescaleTokT(this%pMDFramework, velocities, this%mass, kT)

  end subroutine AndersenThermostat_getInitVelos


  !> Updates the provided velocities according the current temperature.
  subroutine AndersenThermostat_updateVelos(this, velocities)

    !> AndersenThermostat instance.
    type(TAndersenThermostat), intent(inout) :: this

    !> Updated velocities on exit.
    real(dp), intent(inout) :: velocities(:,:)

    real(dp) :: rescaleChance
    real(dp) :: kT
    integer :: ii

    @:ASSERT(all(shape(velocities) <= (/ 3, this%nAtom /)))

    call this%pTempProfile%getTemperature(kT)
    if (this%tRescaleIndiv) then
      do ii = 1, this%nAtom
        call getRandom(this%pRanlux, rescaleChance)
        if (rescaleChance <= this%wvScale) then
          call MaxwellBoltzmann(velocities(:,ii), this%mass(ii), kT, &
              &this%pRanlux)
        end if
      end do
      call restFrame(this%pMDFramework, velocities, this%mass)
    else
      ! all atoms re-set at random
      call getRandom(this%pRanlux, rescaleChance)
      if (rescaleChance <= this%wvScale) then
        do ii = 1, this%nAtom
          call MaxwellBoltzmann(velocities(:,ii), this%mass(ii), kT, &
              &this%pRanlux)
        end do
        call restFrame(this%pMDFramework, velocities, this%mass)
        call rescaleTokT(this%pMDFramework, velocities, this%mass, kT)
      end if
    end if

  end subroutine AndersenThermostat_updateVelos


  !> Outputs internals of thermostat
  subroutine AndersenThermostat_state(this, fd)

    !> instance of thermostat
    type(TAndersenThermostat), intent(in) :: this

    !> filehandle to write out to
    integer,intent(in) :: fd

    ! no internal state, nothing to do

  end subroutine AndersenThermostat_state

end module dftbp_andersentherm
