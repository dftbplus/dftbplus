!--------------------------------------------------------------------------------------------------!
!  DFTB+: general package for performing fast atomistic simulations                                !
!  Copyright (C) 2006 - 2025  DFTB+ developers group                                               !
!                                                                                                  !
!  See the LICENSE file for terms of usage and distribution.                                       !
!--------------------------------------------------------------------------------------------------!

#:include 'common.fypp'

!> Berendsen thermostat - warning non-canonical distribution!
!>
!> If you do not know about the flying icecube do not use this thermostat!
!> Berendsen et al. J. Chem. Phys. 81 3684-3690 (1984).
!> Harvey, Tan and Cheatham, J. Comp. Chem. 19 726-740 (1998).
!>
module dftbp_md_berendsentherm
  use dftbp_common_accuracy, only : dp, minTemp
  use dftbp_io_message, only : error
  use dftbp_math_ranlux, only : TRanlux
  use dftbp_md_mdcommon, only : evalkT, MaxwellBoltzmann, rescaleTokT, restFrame, TMDCommon
  use dftbp_md_tempprofile, only : TTempProfile
  use dftbp_md_thermostat, only : TThermostat
  implicit none

  private
  public :: TBerendsenThermInput
  public :: TBerendsenTherm, TBerendsenTherm_init


  !> Thermostat specific input data for the Berendsen thermostat
  type :: TBerendsenThermInput

    !> Coupling parameter for the thermostat.
    real(dp) :: coupling

  end type TBerendsenThermInput


  !> Berendsen thermostat
  type, extends(TThermostat) :: TBerendsenTherm
    private

    !> Nr. of atoms
    integer :: nAtom

    !> Random number generator
    type(TRanlux), allocatable :: pRanlux

    !> Mass of the atoms
    real(dp), allocatable :: mass(:)

    !> Temperature generator
    type(tTempProfile), pointer :: pTempProfile

    !> Coupling strength to friction term
    real(dp) :: coupling

    !> MD Framework.
    type(TMDCommon) :: pMDFrame

  contains

    procedure :: getInitVelocities => TBerendsenTherm_getInitVelocities
    procedure :: updateVelocities => TBerendsenTherm_updateVelocities
    procedure :: writeState => TBerendsenTherm_writeState

  end type TBerendsenTherm

contains


  !> Initializes a Berendsen thermostat instance
  subroutine TBerendsenTherm_init(this, input, pRanlux, masses, tempProfile, pMDFrame)

    !> Initialised instance on exit.
    type(TBerendsenTherm), intent(out) :: this

    !> Thermostat specific input data
    type(TBerendsenThermInput), intent(in) :: input

    !> Pointer to the random generator.
    type(TRanlux), allocatable, intent(inout) :: pRanlux

    !> Masses of the atoms.
    real(dp), intent(in) :: masses(:)

    !> Temperature profile object.
    type(TTempProfile), pointer, intent(in) :: tempProfile

    !> Molecular dynamics generic framework
    type(TMDCommon), intent(in) :: pMDFrame

    call move_alloc(pRanlux, this%pRanlux)
    this%nAtom = size(masses)
    this%mass = masses
    this%pTempProfile => tempProfile
    this%coupling = input%coupling
    this%pMDFrame = pMDFrame

  end subroutine TBerendsenTherm_init


  !> Returns the initial velocities.
  subroutine TBerendsenTherm_getInitVelocities(this, velocities)

    !> Instance
    class(TBerendsenTherm), intent(inout) :: this

    !> Contains the velocities on return.
    real(dp), intent(out) :: velocities(:,:)

    real(dp) :: kT
    integer :: ii

    @:ASSERT(all(shape(velocities) <= [3, this%nAtom]))

    call this%pTempProfile%getTemperature(kT)
    if (kT < minTemp) then
      call error("Berendsen thermostat not supported at zero temperature")
    end if
    do ii = 1, this%nAtom
      call MaxwellBoltzmann(velocities(:,ii), this%mass(ii), kT, this%pRanlux)
    end do
    call restFrame(this%pMDFrame, velocities, this%mass)
    call rescaleTokT(this%pMDFrame, velocities, this%mass, kT)

  end subroutine TBerendsenTherm_getInitVelocities


  !> Updates the provided velocities according the current temperature.
  !> Shifts to rest frame coordinates if required - this removes some of the flying icecube effect.
  subroutine TBerendsenTherm_updateVelocities(this, velocities)

    !> Instance
    class(TBerendsenTherm), intent(inout) :: this

    !> Updated velocities on exit
    real(dp), intent(inout) :: velocities(:,:)

    real(dp) :: kTCurrent, kTTarget, scaling

    @:ASSERT(all(shape(velocities) <= [3, this%nAtom]))

    call this%pTempProfile%getTemperature(kTTarget)
    call evalkT(this%pMDFrame, kTCurrent,velocities,this%mass)
    scaling = sqrt(1.0_dp + this%coupling * (kTTarget / kTCurrent - 1.0_dp))
    velocities(:,:) = scaling * velocities
    call restFrame(this%pMDFrame, velocities, this%mass)

  end subroutine TBerendsenTherm_updateVelocities


  !> Writes internals of thermostat
  subroutine TBerendsenTherm_writeState(this, fd)

    !> thermostat object
    class(TBerendsenTherm), intent(in) :: this

    !> file unit
    integer,intent(in) :: fd

    ! no internal state, nothing to do

  end subroutine TBerendsenTherm_writeState

end module dftbp_md_berendsentherm
