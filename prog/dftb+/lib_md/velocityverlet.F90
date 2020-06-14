!--------------------------------------------------------------------------------------------------!
!  DFTB+: general package for performing fast atomistic simulations                                !
!  Copyright (C) 2006 - 2020  DFTB+ developers group                                               !
!                                                                                                  !
!  See the LICENSE file for terms of usage and distribution.                                       !
!--------------------------------------------------------------------------------------------------!

#:include 'common.fypp'

!> Velocity Verlet intergrator.
module dftbp_velocityverlet
  use dftbp_assert
  use dftbp_accuracy
  use dftbp_thermostat
  use dftbp_fileid
  use dftbp_message
  implicit none
  private

  public :: TVelocityVerlet
  public :: init, next, rescale, state


  !> Data for the integrator.
  type TVelocityVerlet
    private

    !> Nr. of atoms
    integer :: nAtom

    !> time step for the integrator
    real(dp) :: deltaT

    !> list of particle positions
    real(dp), allocatable :: positions(:,: )

    !> list of particle velocities
    real(dp), allocatable :: velocities(:,:)

    !> Thermostat
    type(TThermostat), allocatable :: pThermostat

    !> do we have the v(t-.5) internal velocity state?
    logical :: vHalfPresent = .false.

    !> do we have a barostat?
    logical :: tBarostat

    !> Strength of Berendsen coupling
    real(dp) :: BarostatStrength

    !> Pressure tensor
    real(dp) :: Pressure(3,3)

    !> is the cell scaling isotropic
    logical :: tIsotropic = .true.
  end type TVelocityVerlet


  !> initialise MD
  interface init
    module procedure VelocityVerlet_themostats
    module procedure VelocityVerlet_velocities
    module procedure VV_themostats_pressure
    module procedure VV_velocities_pressure
  end interface


  !> next geometry step
  interface next
    module procedure VelocityVerlet_next
  end interface


  !> Adjust velocities
  interface rescale
    module procedure VelocityVerlet_rescale
  end interface


  !> write state of the integrator
  interface state
    module procedure VelocityVerlet_state
  end interface

contains


  !> Creates a VelocityVerlet object from the thermostat settings
  subroutine VelocityVerlet_themostats(self, deltaT, positions, pThermostat)

    !> Initialised object on exit.
    type(TVelocityVerlet), intent(out) :: self

    !> Integration time step.
    real(dp), intent(in) :: deltaT

    !> Position of the atoms.
    real(dp), intent(in) :: positions(:,:)

    !> Thermostat if needed.
    type(TThermostat), allocatable, intent(inout) :: pThermostat

    @:ASSERT(size(positions, dim=1) == 3)

    self%nAtom = size(positions, dim=2)
    allocate(self%velocities(3, self%nAtom))
    allocate(self%positions(3, self%nAtom))

    self%deltaT = deltaT
    self%positions(:,:) = positions(:,:)
    call move_alloc(pThermostat, self%pThermostat)

    call getInitVelocities(self%pThermostat, self%velocities)

    self%vHalfPresent = .false. ! no we don't have the t-.5 velocities

    self%tBarostat = .false.

  end subroutine VelocityVerlet_themostats


  !> Creates a VelocityVerlet object from given external velocities for the t-th time step, this
  !> means later we have to reconstruct the Vel. Verlet t+.5 velocities
  subroutine VelocityVerlet_velocities(self, deltaT, positions, pThermostat, &
      & velocities)

    !> Initialised object on exit.
    type(TVelocityVerlet), intent(out) :: self

    !> Integration time step.
    real(dp), intent(in) :: deltaT

    !> Position of the atoms.
    real(dp), intent(in) :: positions(:,:)

    !> Thermostat.
    type(TThermostat), allocatable, intent(inout) :: pThermostat

    !> List of initial velocities
    real(dp), intent(in) :: velocities(:,:)

    @:ASSERT(size(positions, dim=1) == 3)

    self%nAtom = size(positions, dim=2)
    allocate(self%velocities(3, self%nAtom))
    allocate(self%positions(3, self%nAtom))

    self%deltaT = deltaT
    self%positions(:,:) = positions(:,:)
    call move_alloc(pThermostat, self%pThermostat)

    self%velocities(:,:) = velocities(:,:)

    ! assumes the V read in corresponds to the current coordinates, so we should reconstruct the
    ! t+.5 velocities when possible once forces are available for the coordinates
    self%vHalfPresent = .false.

    self%tBarostat = .false.

  end subroutine VelocityVerlet_velocities


  !> Creates a VelocityVerlet object from the thermostat settings and isotropic pressure
  subroutine VV_themostats_pressure(self, deltaT, positions, pThermostat, &
      & Barostat, Pressure, tIsotropic)

    !> Initialised object on exit.
    type(TVelocityVerlet), intent(out) :: self

    !> Integration time step.
    real(dp), intent(in) :: deltaT

    !> Position of the atoms.
    real(dp), intent(in) :: positions(:,:)

    !> Thermostat if needed.
    type(TThermostat), allocatable, intent(inout) :: pThermostat

    !> Coupling strength.
    real(dp), intent(in) :: Barostat

    !> Target isotropic pressure
    real(dp), intent(in) :: Pressure

    !> Is this an isotropic barostat, or can the cell shape change?
    logical, intent(in) :: tIsotropic

    integer :: ii

    @:ASSERT(size(positions, dim=1) == 3)

    self%nAtom = size(positions, dim=2)
    allocate(self%velocities(3, self%nAtom))
    allocate(self%positions(3, self%nAtom))

    self%deltaT = deltaT
    self%positions(:,:) = positions(:,:)
    call move_alloc(pThermostat, self%pThermostat)

    call getInitVelocities(self%pThermostat, self%velocities)

    self%vHalfPresent = .true. ! yes we have the t-.5 velocities

    self%tBarostat = .true.
    self%BarostatStrength = Barostat
    self%Pressure = 0.0_dp
    do ii = 1, 3
      self%Pressure(ii,ii) = pressure
    end do

    self%tIsotropic = tIsotropic

  end subroutine VV_themostats_pressure


  !> Creates a VelocityVerlet object from given external velocities for the t-th time step, this
  !> means later we have to reconstruct the Vel. Verlet t+.5 velocities and barostat isotropic
  !> pressure
  subroutine VV_velocities_pressure(self, deltaT, positions, pThermostat, &
      & velocities, Barostat, Pressure, tIsotropic)

    !> Initialised object on exit.
    type(TVelocityVerlet), intent(out) :: self

    !> Integration time step.
    real(dp), intent(in) :: deltaT

    !> Position of the atoms.
    real(dp), intent(in) :: positions(:,:)

    !> Thermostat.
    type(TThermostat), allocatable, intent(inout) :: pThermostat

    !> List of initial velocities
    real(dp), intent(in) :: velocities(:,:)

    !> Coupling strength
    real(dp), intent(in) :: Barostat

    !> Target isotropic pressure
    real(dp), intent(in) :: Pressure

    !> Is this an isotropic barostat, or can the cell shape change?
    logical, intent(in) :: tIsotropic

    integer :: ii

    @:ASSERT(size(positions, dim=1) == 3)

    self%nAtom = size(positions, dim=2)
    allocate(self%velocities(3, self%nAtom))
    allocate(self%positions(3, self%nAtom))

    self%deltaT = deltaT
    self%positions(:,:) = positions(:,:)
    call move_alloc(pThermostat, self%pThermostat)

    self%velocities(:,:) = velocities(:,:)

    ! assumes the V read in corresponds to the current coordinates, so we should reconstruct the
    ! t+.5 velocities when possible once forces are available for the coordinates
    self%vHalfPresent = .false.

    self%tBarostat = .true.
    self%BarostatStrength = Barostat

    self%Pressure = 0.0_dp
    do ii = 1, 3
      self%Pressure(ii,ii) = pressure
    end do
    self%tIsotropic = tIsotropic

  end subroutine VV_velocities_pressure


  !> Takes a timestep for the MD integrator, optionally with a thermostat.
  !> Due to the way the velocity Verlet is split, the returned velocity is for 1 complete MD step
  !> behind the returned positions so print positions, then call next and then print velocities to
  !> get agreement between the positions and velocities.
  subroutine VelocityVerlet_next(self, accel, newCoord, newVelocity)

    !> Integrator to propogate
    type(TVelocityVerlet), intent(inout) :: self

    !> Accelerations.
    real(dp),intent(in) :: accel(:,:)

    !> Displaced coordinates
    real(dp),intent(out) :: newCoord(:,:)

    !> Velocity of displaced coords
    real(dp),intent(out) :: newVelocity(:,:)

    newCoord(:,:) = 0.0_dp
    newVelocity(:,:) = 0.0_dp

    ! start from the usual ordering of velocity verlet method (two cycles
    ! shown):
    ! a.1 v(t+.5dt) = v(t)        + .5*a(t)*dt -- a(t) external
    ! a.2 r(t + dt) = r(t)        + v(t+.5dt)*dt
    ! a.3 v(t + dt) = v(t+.5dt)   + .5*a(t+dt)*dt -- a(t+dt) external call
    ! b.1 v(t+1.5dt) = v(t+dt)    + .5*a(t+dt)*dt -- a(t+dt) external
    ! b.2 r(t + 2dt) = r(t+dt)    + v(t+1.5dt)*dt
    ! b.3 v(t + 2dt) = v(t+1.5dt) + .5*a(t+2dt)*dt -- a(t+2dt) external call
    !
    ! and cut out a.3 b.1 b.2 as the cycle :
    ! a.3 v(t)      = v(t-.5dt)+ .5*a(t)*dt -- a(t) input
    ! b.1 v(t+.5dt) = v(t)     + .5*a(t)*dt
    ! b.2 r(t+dt)   = r(t)     + v(t+.5dt)*dt
    !
    ! so :
    ! v_out   = v_store + .5*a_input*dt
    ! v_store = v_out   + .5*a_input*dt
    ! r_out   = r_store + v_store*dt
    ! r_store = r_out
    ! where v_out is one MD step behind the positions returned.

    if (self%vHalfPresent) then
      newVelocity(:,:) = self%velocities(:,:) &
          & + 0.5_dp * accel(:,:) * self%deltaT
    else
      newVelocity(:,:) = self%velocities(:,:)
      self%vHalfPresent=.true.
    end if

    self%velocities(:,:) = newVelocity(:,:) + 0.5_dp * accel(:,:) * self%deltaT
    newCoord(:,:) = self%positions(:,:) + self%velocities(:,:) * self%deltaT
    self%positions(:,:) = newCoord(:,:)

    if (allocated(self%pThermostat)) then
      call updateVelocities(self%pThermostat, self%velocities)
    end if

  end subroutine VelocityVerlet_next


  !> Rescale the cell parameters and coordinates according to the tensorial version of the Berensen
  !> barostat (allows cell shape changes if the external pressure/stress is non isotropic)
  !> The forms of the isotropic and anisotropic Beresdsen barostats in the literature are slightly
  !> incompatible in their definitions
  subroutine VelocityVerlet_rescale(self,coord,latVecs,pressureTensor)

    !> Integrator to rescale
    type(TVelocityVerlet), intent(inout) :: self

    !> Atom coordinates to rescale
    real(dp),intent(inout) :: coord(:,:)

    !> Lattice vectors to rescale
    real(dp),intent(inout) :: latVecs(3,3)

    !> System stress tensor
    real(dp),intent(in) :: pressureTensor(3,3)

    real(dp) :: scale(3,3)
    real(dp) :: scaleIso, Pext, P
    integer :: ii

    @:ASSERT(self%tBarostat)

    ! isotropic Berendsen, not quite consistent with anisotropic but its in the literature...
    if (self%tIsotropic) then
      Pext = 0.0_dp
      P = 0.0_dp
      do ii = 1, 3
        Pext = self%Pressure(ii,ii) / 3.0_dp
        P = P + pressureTensor(ii,ii) / 3.0_dp
      end do
      scaleIso = (1.0_dp - self%BarostatStrength*(Pext - P))**(1.0_dp/3.0_dp)
      self%positions(:,:) = self%positions(:,:) * scaleIso
      coord(:,:) = coord(:,:) * scaleIso
      latVecs(:,:) = latVecs(:,:) * scaleIso
    else
      scale = 0.0_dp
      do ii = 1, 3
        scale(ii,ii) = 1.0_dp
      end do
      scale = scale - self%BarostatStrength*(self%Pressure-pressureTensor)
      do ii = 1, self%nAtom
        self%positions(:,ii) = matmul(self%positions(:,ii),scale)
        coord(:,ii) = matmul(coord(:,ii),scale)
      end do
      latVecs(:,:) = matmul(latVecs(:,:),scale)
    end if

  end subroutine VelocityVerlet_rescale


  !> Outputs internals of MD integrator
  subroutine VelocityVerlet_state(self,fd)

    !> instance of integrator
    type(TVelocityVerlet), intent(in) :: self

    !> filehandle to write out to
    integer,intent(in) :: fd

    if (allocated(self%pThermostat)) then
      call state(self%pThermostat,fd)
    end if

  end subroutine VelocityVerlet_state

end module dftbp_velocityverlet
