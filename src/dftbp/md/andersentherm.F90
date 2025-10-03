!--------------------------------------------------------------------------------------------------!
!  DFTB+: general package for performing fast atomistic simulations                                !
!  Copyright (C) 2006 - 2025  DFTB+ developers group                                               !
!                                                                                                  !
!  See the LICENSE file for terms of usage and distribution.                                       !
!--------------------------------------------------------------------------------------------------!

#:include 'common.fypp'

!> Implements the Andersen thermostat
!>
!> Two versions of the Andersen thermostat are implemented, either the global re-select or per atom
!> reselect of velocities from the Maxwell-Boltzmann distribution See Andersen J. Chem. Phys. 72.
!> 2384 (1980)
!>
module dftbp_md_andersentherm
  use dftbp_common_accuracy, only : dp, minTemp
  use dftbp_io_message, only : error
  use dftbp_math_ranlux, only : getRandom, TRanlux
  use dftbp_md_mdcommon, only : MaxwellBoltzmann, rescaleTokT, restFrame, TMDCommon
  use dftbp_md_tempprofile, only : TTempProfile
  use dftbp_md_thermostat, only : TThermostat
  implicit none

  private
  public :: TAndersenThermInput
  public :: TAndersenTherm, TAndersenTherm_init


  !> Thermostat specific input data for the Andersen thermostat
  type :: TAndersenThermInput

    !> If velocities should be rescaled per atom
    logical :: rescaleIndiv

    !> Rescaling probability.
    real(dp) :: wvScale

  end type TAndersenThermInput


  !> Andersen thermostat
  type, extends(TThermostat) :: TAndersenTherm
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

  contains

    procedure :: getInitVelocities => TAndersenTherm_getInitVelocities
    procedure :: updateVelocities => TAndersenTherm_updateVelocities
    procedure :: writeState => TAndersenTherm_writeState

  end type TAndersenTherm

contains


  !> Creates an Andersen thermostat instance.
  subroutine TAndersenTherm_init(this, input, pRanlux, masses, tempProfile, pMDFramework)

    !> Initialised instance on exit.
    type(TAndersenTherm), intent(out) :: this

    !> Thermostat specific input data
    type(TAndersenThermInput), intent(in) :: input

    !> Random generator
    type(TRanlux), allocatable, intent(inout) :: pRanlux

    !> Masses of the atoms.
    real(dp), intent(in) :: masses(:)

    !> Pointer to a temperature profile object.
    type(TTempProfile), pointer, intent(in) :: tempProfile

    !> Molecular dynamics general specifications
    type(TMDCommon), intent(in) :: pMDFramework

    call move_alloc(pRanlux, this%pRanlux)
    this%nAtom = size(masses)
    this%mass = masses
    this%pTempProfile => tempProfile
    this%tRescaleIndiv = input%rescaleIndiv
    this%wvScale = input%wvScale
    this%pMDFramework = pMDFramework

  end subroutine TAndersenTherm_init


  !> Returns the initial velocities.
  subroutine TAndersenTherm_getInitVelocities(this, velocities)

    !> Instance
    class(TAndersenTherm), intent(inout) :: this

    !> Velocities on return.
    real(dp), intent(out) :: velocities(:,:)

    real(dp) :: kT
    integer :: ii

    @:ASSERT(all(shape(velocities) <= [3, this%nAtom]))

    call this%pTempProfile%getTemperature(kT)
    if (kT < minTemp) then
      call error("Andersen thermostat not supported at zero temperature")
    end if
    do ii = 1, this%nAtom
      call MaxwellBoltzmann(velocities(:,ii), this%mass(ii), kT, this%pRanlux)
    end do
    call restFrame(this%pMDFramework, velocities, this%mass)
    call rescaleTokT(this%pMDFramework, velocities, this%mass, kT)

  end subroutine TAndersenTherm_getInitVelocities


  !> Updates the provided velocities according the current temperature.
  subroutine TAndersenTherm_updateVelocities(this, velocities)

    !> Instance
    class(TAndersenTherm), intent(inout) :: this

    !> Updated velocities on exit.
    real(dp), intent(inout) :: velocities(:,:)

    real(dp) :: rescaleChance
    real(dp) :: kT
    integer :: ii

    @:ASSERT(all(shape(velocities) <= [3, this%nAtom]))

    call this%pTempProfile%getTemperature(kT)
    if (this%tRescaleIndiv) then
      do ii = 1, this%nAtom
        call getRandom(this%pRanlux, rescaleChance)
        if (rescaleChance <= this%wvScale) then
          call MaxwellBoltzmann(velocities(:,ii), this%mass(ii), kT, this%pRanlux)
        end if
      end do
      call restFrame(this%pMDFramework, velocities, this%mass)
    else
      ! all atoms re-set at random
      call getRandom(this%pRanlux, rescaleChance)
      if (rescaleChance <= this%wvScale) then
        do ii = 1, this%nAtom
          call MaxwellBoltzmann(velocities(:,ii), this%mass(ii), kT, this%pRanlux)
        end do
        call restFrame(this%pMDFramework, velocities, this%mass)
        call rescaleTokT(this%pMDFramework, velocities, this%mass, kT)
      end if
    end if

  end subroutine TAndersenTherm_updateVelocities


  !> Writes internals of thermostat
  subroutine TAndersenTherm_writeState(this, fd)

    !> instance of thermostat
    class(TAndersenTherm), intent(in) :: this

    !> filehandle to write out to
    integer,intent(in) :: fd

    ! no internal state, nothing to do

  end subroutine TAndersenTherm_writeState

end module dftbp_md_andersentherm
