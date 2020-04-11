!--------------------------------------------------------------------------------------------------!
!  DFTB+: general package for performing fast atomistic simulations                                !
!  Copyright (C) 2006 - 2020  DFTB+ developers group                                               !
!                                                                                                  !
!  See the LICENSE file for terms of usage and distribution.                                       !
!--------------------------------------------------------------------------------------------------!

#:include 'common.fypp'

!> Berendsen thermostat - warning non-canonical distribution!
!> If you do not know about the flying icecube do not use this
!> thermostat!
!> Berendsen et al. J. Chem. Phys. 81 3684-3690 (1984).
!> Harvey, Tan and Cheatham, J. Comp. Chem. 19 726-740 (1998).
module dftbp_berendsentherm
  use dftbp_assert
  use dftbp_accuracy
  use dftbp_mdcommon
  use dftbp_ranlux
  use dftbp_tempprofile
  implicit none
  private

  public :: TBerendsenThermostat
  public :: init, getInitVelocities, updateVelocities, state


  !> Data for the Berendsen thermostat
  type TBerendsenThermostat
    private

    !> Nr. of atoms
    integer :: nAtom

    !> Random number generator
    type(TRanlux), allocatable :: pRanlux

    !> Mass of the atoms
    real(dp), allocatable :: mass(:)

    !> Temperature generator
    type(tTempProfile), pointer :: pTempProfile

    !> coupling strength to friction term
    real(dp) :: couplingParameter

    !> MD Framework.
    type(TMDCommon) :: pMDFrame
  end type TBerendsenThermostat


  !> initialise object
  interface init
    module procedure Berendsen_init
  end interface


  !> initial thermal velocities if needed
  interface getInitVelocities
    module procedure Berendsen_getInitVelos
  end interface


  !> update velocites acording to the thermostat
  interface updateVelocities
    module procedure Berendsen_updateVelos
  end interface


  !> write state of the thermostat
  interface state
    module procedure Berendsen_state
  end interface

contains


  !> Creates an Berendsen thermostat instance.
  subroutine Berendsen_init(self, pRanlux, masses, tempProfile, &
      & couplingParameter, pMDFrame)

    !> Initialised instance on exit.
    type(TBerendsenThermostat), intent(out) :: self

    !> Pointer to the random generator.
    type(TRanlux), allocatable, intent(inout) :: pRanlux

    !> Masses of the atoms.
    real(dp), intent(in) :: masses(:)

    !> Temperature profile object.
    type(TTempProfile), pointer, intent(in) :: tempProfile

    !> Coupling parameter for the thermostat.
    real(dp), intent(in) :: couplingParameter

    !> Molecular dynamics generic framework
    type(TMDCommon), intent(in) :: pMDFrame

    call move_alloc(pRanlux, self%pRanlux)
    self%nAtom = size(masses)
    allocate(self%mass(self%nAtom))
    self%mass(:) = masses(:)
    self%pTempProfile => tempProfile
    self%couplingParameter = couplingParameter
    self%pMDFrame = pMDFrame

  end subroutine Berendsen_init


  !> Returns the initial velocities.
  subroutine Berendsen_getInitVelos(self, velocities)

    !> BerendsenThermostat instance.
    type(TBerendsenThermostat), intent(inout) :: self

    !> Contains the velocities on return.
    real(dp), intent(out) :: velocities(:,:)

    real(dp) :: kT
    integer :: ii

    @:ASSERT(all(shape(velocities) <= (/ 3, self%nAtom /)))

    call getTemperature(self%pTempProfile, kT)
    do ii = 1, self%nAtom
      call MaxwellBoltzmann(velocities(:,ii), self%mass(ii), kT, self%pRanlux)
    end do
    call restFrame(self%pMDFrame, velocities, self%mass)
    call rescaleTokT(self%pMDFrame, velocities, self%mass, kT)

  end subroutine Berendsen_getInitVelos


  !> Updates the provided velocities according the current temperature.
  !> Shifts to rest frame coordinates if required - this removes some of the flying icecube effect.
  subroutine Berendsen_updateVelos(self, velocities)

    !> Thermostat instance.
    type(TBerendsenThermostat), intent(inout) :: self

    !> Updated velocities on exit.
    real(dp), intent(inout) :: velocities(:,:)

    real(dp) :: kTCurrent, kTTarget, scaling

    @:ASSERT(all(shape(velocities) <= (/ 3, self%nAtom /)))

    call getTemperature(self%pTempProfile, kTTarget)
    call evalkT(self%pMDFrame, kTCurrent,velocities,self%mass)
    scaling = sqrt(1.0_dp + self%couplingParameter*(kTTarget/kTCurrent-1.0_dp))
    velocities(:,:) = scaling * velocities(:,:)
    call restFrame(self%pMDFrame, velocities, self%mass)

  end subroutine Berendsen_updateVelos


  !> Outputs internals of thermostat
  subroutine Berendsen_state(self, fd)

    !> thermostat object
    type(TBerendsenThermostat), intent(in) :: self

    !> file unit
    integer,intent(in) :: fd

    ! no internal state, nothing to do

  end subroutine Berendsen_state

end module dftbp_berendsentherm
