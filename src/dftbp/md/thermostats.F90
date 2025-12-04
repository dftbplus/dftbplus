!--------------------------------------------------------------------------------------------------!
!  DFTB+: general package for performing fast atomistic simulations                                !
!  Copyright (C) 2006 - 2025  DFTB+ developers group                                               !
!                                                                                                  !
!  See the LICENSE file for terms of usage and distribution.                                       !
!--------------------------------------------------------------------------------------------------!

#:include "common.fypp"

!> Provides all the thermostats to the generic thermostat interface
module dftbp_md_thermostats
  use dftbp_md_mdcommon, only : TMDCommon
  use dftbp_math_ranlux, only : TRanlux
  use dftbp_md_tempprofile, only : TTempProfile
  use dftbp_common_accuracy, only : dp
  use dftbp_md_andersentherm, only : TAndersenTherm, TAndersenTherm_init, TAndersenThermInput
  use dftbp_md_berendsentherm, only : TBerendsenTherm, TBerendsenTherm_init, TBerendsenThermInput
  use dftbp_md_dummytherm, only : TDummyTherm, TDummyTherm_init
  use dftbp_md_nhctherm, only : TNhcTherm, TNhcTherm_init, TNhcThermInput
  use dftbp_md_thermostat, only : TThermostat
  implicit none

  private
  public :: thermostatTypes
  public :: createThermostat, TThermostat, TThermostatInput
  public :: TAndersenTherm, TAndersenTherm_init, TAndersenThermInput
  public :: TBerendsenTherm, TBerendsenTherm_init, TBerendsenThermInput
  public :: TDummyTherm, TDummyTherm_init
  public :: TNhcTherm, TNhcTherm_init, TNhcThermInput


  type :: TThermostatTypes_
    integer :: dummy = 0
    integer :: andersen = 1
    integer :: berendsen = 2
    integer :: nhc = 3
  end type TThermostatTypes_

  !> Available thermostat types
  type(TThermostatTypes_), parameter :: thermostatTypes = TThermostatTypes_()


  !> Collected input data for available thermostats
  type :: TThermostatInput
    integer :: thermostatType
    type(TAndersenThermInput), allocatable :: andersen
    type(TBerendsenThermInput), allocatable :: berendsen
    type(TNhcThermInput), allocatable :: nhc
  end type TThermostatInput

contains


  !> Creates a thermostat instance based on the provided input data
  subroutine createThermostat(thermostat, input, masses, randomThermostat, pMDFrame,&
        & pTempProfile, deltaT)

    !> Created thermostat instance on exit
    class(TThermostat), allocatable, intent(out) :: thermostat

    !> Thermostat input data
    type(TThermostatInput), intent(in) :: input

    !> Masses of the thermostated atoms
    real(dp), intent(in) :: masses(:)

    !> Random number generator for the thermostat
    type(TRanlux), allocatable, intent(inout) :: randomThermostat

    !> Molecular dynamics framework
    type(TMDCommon), intent(in) :: pMdFrame

    !> Temperature profile object
    type(TTempProfile), pointer, intent(in) :: ptempProfile

    !> Time step for the MD simulation
    real(dp), intent(in) :: deltaT

    type(TDummyTherm), allocatable :: dummyTherm
    type(TAndersenTherm), allocatable :: andersenTherm
    type(TBerendsenTherm), allocatable :: berendsenTherm
    type(TNHCTherm), allocatable :: nhcTherm
    real(dp) :: tempAtom

    @:ASSERT(input%thermostatType >= 0 .and. input%thermostatType <= 3)
    @:ASSERT(all(&
        & ([thermostatTypes%andersen, thermostatTypes%berendsen, thermostatTypes%nhc]&
        & == input%thermostatType)&
        & .eqv. &
        & [allocated(input%andersen), allocated(input%berendsen), allocated(input%nhc)]&
        & ))

    select case (input%thermostatType)
    case (thermostatTypes%dummy)
      allocate(dummyTherm)
      call pTempProfile%getTemperature(tempAtom)
      call TDummyTherm_init(dummyTherm, tempAtom, masses, randomThermostat, pMdFrame)
      call move_alloc(dummyTherm, thermostat)
    case (thermostatTypes%andersen)
      allocate(andersenTherm)
      call TAndersenTherm_init(andersenTherm, input%andersen, randomThermostat, masses,&
          ptempProfile, pMdFrame)
      call move_alloc(andersenTherm, thermostat)
    case (thermostatTypes%berendsen)
      allocate(berendsenTherm)
      call TBerendsenTherm_init(berendsenTherm, input%berendsen, randomThermostat, masses,&
          ptempProfile, pMdFrame)
      call move_alloc(berendsenTherm, thermostat)
    case (thermostatTypes%nhc)
      allocate(nhcTherm)
      call TNhcTherm_init(nhcTherm, input%nhc, randomThermostat, masses, ptempProfile,&
         pMdFrame, deltaT)
      call move_alloc(nhcTherm, thermostat)
    end select

  end subroutine createThermostat

end module dftbp_md_thermostats
