!--------------------------------------------------------------------------------------------------!
!  DFTB+: general package for performing fast atomistic simulations                                !
!  Copyright (C) 2017  DFTB+ developers group                                                      !
!                                                                                                  !
!  See the LICENSE file for terms of usage and distribution.                                       !
!--------------------------------------------------------------------------------------------------!

#:include 'common.fypp'

!!* Berendsen thermostat - warning non-canonical distribution!
!!* @note If you do not know about the flying icecube do not use this
!!* thermostat!
!!* @ref Berendsen <I>et al.</I> J. Chem. Phys. 81 3684-3690 (1984).
!!* @ref Harvey, Tan and Cheatham, J. Comp. Chem. 19 726-740 (1998).
module berendsentherm
  use assert
  use accuracy
  use mdcommon
  use ranlux
  use tempprofile
  implicit none
  private

  public :: OBerendsenThermostat
  public :: init, getInitVelocities, updateVelocities, state

  !!* Data for the Berendsen thermostat
  type OBerendsenThermostat
    private
    integer :: nAtom                    !* Nr. of atoms
    type(ORanlux), allocatable :: pRanlux   !* Random number generator
    real(dp), allocatable :: mass(:)        !* Mass of the atoms
    type(OTempProfile), pointer :: pTempProfile  !* Temperature generator
    real(dp) :: couplingParameter       !* coupling strength to friction term
    type(OMDCommon) :: pMDFrame !* MD Framework.
  end type OBerendsenThermostat

  interface init
    module procedure Berendsen_init
  end interface

  interface getInitVelocities
    module procedure Berendsen_getInitVelos
  end interface

  interface updateVelocities
    module procedure Berendsen_updateVelos
  end interface

  interface state
    module procedure Berendsen_state
  end interface

contains

  !!* Creates an Berendsen thermostat instance.
  !!* @param self Initialised instance on exit.
  !!* @param pRanlux Pointer to the random generator.
  !!* @param masses Masses of the atoms.
  !!* @param tempProfile Pointer to a temperature profile object.
  !!* @param couplingParameter Coupling parameter for the thermostat.
  subroutine Berendsen_init(self, pRanlux, masses, tempProfile, &
      & couplingParameter, pMDFrame)
    type(OBerendsenThermostat), intent(out) :: self
    type(ORanlux), allocatable, intent(inout) :: pRanlux
    real(dp), intent(in) :: masses(:)
    type(OTempProfile), pointer, intent(in) :: tempProfile
    real(dp), intent(in) :: couplingParameter
    type(OMDCommon) :: pMDFrame

    call move_alloc(pRanlux, self%pRanlux)
    self%nAtom = size(masses)
    allocate(self%mass(self%nAtom))
    self%mass(:) = masses(:)
    self%pTempProfile => tempProfile
    self%couplingParameter = couplingParameter
    self%pMDFrame = pMDFrame

  end subroutine Berendsen_init


  !!* Returns the initial velocities.
  !!* @param self BerendsenThermostat instance.
  !!* @param velocities Contains the velocities on return.
  subroutine Berendsen_getInitVelos(self, velocities)
    type(OBerendsenThermostat), intent(inout) :: self
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



  !!* Updates the provided velocities according the current temperature.
  !!* @param self BerendsenThermostat instance.
  !!* @param velocities Updated velocities on exit.
  !!* @note shifts to rest frame coordinates if required - this removes
  !!* some of the flying icecube effect.
  subroutine Berendsen_updateVelos(self, velocities)
    type(OBerendsenThermostat), intent(inout) :: self
    real(dp), intent(inout) :: velocities(:,:)

    real(dp) :: kTCurrent, kTTarget, scaling

    @:ASSERT(all(shape(velocities) <= (/ 3, self%nAtom /)))

    call getTemperature(self%pTempProfile, kTTarget)
    call evalkT(self%pMDFrame, kTCurrent,velocities,self%mass)
    scaling = sqrt(1.0_dp + self%couplingParameter*(kTTarget/kTCurrent-1.0_dp))
    velocities(:,:) = scaling * velocities(:,:)
    call restFrame(self%pMDFrame, velocities, self%mass)

  end subroutine Berendsen_updateVelos

  subroutine Berendsen_state(self, fd)
    type(OBerendsenThermostat), intent(in) :: self
    integer,intent(in)                  :: fd

    ! no internal state, nothing to do

  end subroutine Berendsen_state

end module berendsentherm
