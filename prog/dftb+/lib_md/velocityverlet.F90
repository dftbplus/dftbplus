!--------------------------------------------------------------------------------------------------!
!  DFTB+: general package for performing fast atomistic simulations                                !
!  Copyright (C) 2017  DFTB+ developers group                                                      !
!                                                                                                  !
!  See the LICENSE file for terms of usage and distribution.                                       !
!--------------------------------------------------------------------------------------------------!

#:include 'common.fypp'

!!* Velocity Verlet intergrator.
module velocityverlet
  use assert
  use accuracy
  use thermostat
  use fileid
  use message
  implicit none
  private

  public :: OVelocityVerlet
  public :: init, next, rescale, state

  !!* Data for the integrator.
  type OVelocityVerlet
    private
    integer :: nAtom                     !* Nr. of atoms
    real(dp) :: deltaT                   !* time step for the integrator
    real(dp), allocatable :: positions(:,: ) !* list of particle positions
    real(dp), allocatable :: velocities(:,:) !* list of particle velocities
    type(OThermostat), allocatable :: pThermostat  !* Thermostat
    logical           :: vHalfPresent = .false. !* do we have the v(t-.5)
    !* internal velocity state?
    logical  :: tBarostat                !* do we have a barostat?
    real(dp) :: BarostatStrength         !* Strength of Berendsen coupling
    real(dp) :: Pressure(3,3)            !* Pressure tensor
    logical  :: tIsotropic = .true.      !* is the cell scaling isotropic
  end type OVelocityVerlet

  interface init
    module procedure VelocityVerlet_themostats
    module procedure VelocityVerlet_velocities
    module procedure VV_themostats_pressure
    module procedure VV_velocities_pressure
  end interface

  interface next
    module procedure VelocityVerlet_next
  end interface

  interface rescale
    module procedure VelocityVerlet_rescale
  end interface

  interface state
    module procedure VelocityVerlet_state
  end interface

contains

  !!* Creates a VelocityVerlet object from the thermostat settings
  !!* @param self Pointer to the initialised object on exit.
  !!* @param deltaT Integration time step.
  !!* @param positions Position of the atoms.
  !!* @param pThermostat Pointer to a thermostat if needed.
  subroutine VelocityVerlet_themostats(self, deltaT, positions, pThermostat)
    type(OVelocityVerlet), intent(out) :: self
    real(dp), intent(in)                 :: deltaT
    real(dp), intent(in)                 :: positions(:,:)
    type(OThermostat), allocatable, intent(inout) :: pThermostat

    @:ASSERT(size(positions, dim=1) == 3)

    self%nAtom = size(positions, dim=2)
    allocate(self%velocities(3, self%nAtom))
    allocate(self%positions(3, self%nAtom))

    self%deltaT = deltaT
    self%positions(:,:) = positions(:,:)
    call move_alloc(pThermostat, self%pThermostat)

    call getInitVelocities(self%pThermostat, self%velocities)

    self%vHalfPresent = .false. ! no we dont have the t-.5 velocities

    self%tBarostat = .false.

  end subroutine VelocityVerlet_themostats


  !!* Creates a VelocityVerlet object from given external velocities for the
  !!* t-th time step, this means later we have to reconstruct the Vel. Verlet
  !!* t+.5 velocities
  !!* @param self Pointer to the initialised object on exit.
  !!* @param deltaT Integration time step.
  !!* @param positions Position of the atoms.
  !!* @param pThermostat Pointer to a thermostat.
  !!* @param velocities list of initial velocities
  subroutine VelocityVerlet_velocities(self, deltaT, positions, pThermostat, &
      & velocities)
    type(OVelocityVerlet), intent(out) :: self
    real(dp), intent(in)                 :: deltaT
    real(dp), intent(in)                 :: positions(:,:)
    type(OThermostat), allocatable, intent(inout) :: pThermostat
    real(dp), intent(in)                 :: velocities(:,:)

    @:ASSERT(size(positions, dim=1) == 3)

    self%nAtom = size(positions, dim=2)
    allocate(self%velocities(3, self%nAtom))
    allocate(self%positions(3, self%nAtom))

    self%deltaT = deltaT
    self%positions(:,:) = positions(:,:)
    call move_alloc(pThermostat, self%pThermostat)

    self%velocities(:,:) = velocities(:,:)

    self%vHalfPresent = .false. ! assumes the V read in corresponds to the
    ! current coordinates, so we should reconstruct the t+.5 velocities when
    ! possible once forces are available for the coordinates

    self%tBarostat = .false.

  end subroutine VelocityVerlet_velocities

  !!* Creates a VelocityVerlet object from the thermostat settings and
  !!* isotropic pressure
  !!* @param self Pointer to the initialised object on exit.
  !!* @param deltaT Integration time step.
  !!* @param positions Position of the atoms.
  !!* @param pThermostat Pointer to a thermostat if needed.
  !!* @param Barostat coupling strength
  !!* @param Pressure target isotropic pressure
  !!* @param tIsotropic is this an isotropic barostat, or can the cell shape
  !!* change?
  subroutine VV_themostats_pressure(self, deltaT, positions, pThermostat, &
      & Barostat, Pressure, tIsotropic)
    type(OVelocityVerlet), intent(out) :: self
    real(dp), intent(in)                 :: deltaT
    real(dp), intent(in)                 :: positions(:,:)
    type(OThermostat), allocatable, intent(inout) :: pThermostat
    real(dp), intent(in)                 :: Barostat
    real(dp), intent(in)                 :: Pressure
    logical, intent(in)                  :: tIsotropic

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


  !!* Creates a VelocityVerlet object from given external velocities for the
  !!* t-th time step, this means later we have to reconstruct the Vel. Verlet
  !!* t+.5 velocities and barostat isotropic pressure
  !!* @param self Pointer to the initialised object on exit.
  !!* @param deltaT Integration time step.
  !!* @param positions Position of the atoms.
  !!* @param pThermostat Pointer to a thermostat.
  !!* @param velocities list of initial velocities
  !!* @param Barostat coupling strength
  !!* @param Pressure target target isotropic pressure
  !!* @param tIsotropic is this an isotropic barostat, or can the cell shape
  !!* change?
  subroutine VV_velocities_pressure(self, deltaT, positions, pThermostat, &
      & velocities, Barostat, Pressure, tIsotropic)
    type(OVelocityVerlet), intent(out) :: self
    real(dp), intent(in)                 :: deltaT
    real(dp), intent(in)                 :: positions(:,:)
    type(OThermostat), allocatable, intent(inout) :: pThermostat
    real(dp), intent(in)                 :: velocities(:,:)
    real(dp), intent(in)                 :: Barostat
    real(dp), intent(in)                 :: Pressure
    logical, intent(in)                  :: tIsotropic

    integer :: ii

    @:ASSERT(size(positions, dim=1) == 3)

    self%nAtom = size(positions, dim=2)
    allocate(self%velocities(3, self%nAtom))
    allocate(self%positions(3, self%nAtom))

    self%deltaT = deltaT
    self%positions(:,:) = positions(:,:)
    call move_alloc(pThermostat, self%pThermostat)

    self%velocities(:,:) = velocities(:,:)

    self%vHalfPresent = .false. ! assumes the V read in corresponds to the
    ! current coordinates, so we should reconstruct the t+.5 velocities when
    ! possible once forces are available for the coordinates

    self%tBarostat = .true.
    self%BarostatStrength = Barostat

    self%Pressure = 0.0_dp
    do ii = 1, 3
      self%Pressure(ii,ii) = pressure
    end do
    self%tIsotropic = tIsotropic

  end subroutine VV_velocities_pressure


  !!* Takes a timestep for the MD integrator, optionally with a thermostat.
  !!* @param self integrator to propogate
  !!* @param accel Accelerations.
  !!* @param newCoord displaced coordinates
  !!* @param newVelocity velocity of displaced coords
  !!* @caveat Due to the way the velocity Verlet is split, the returned
  !!* velocity is for 1 complete MD step behind the returned positions
  !!* so print positions, then call next and then print velocities
  !!* to get agreement between the positions and velocities.
  subroutine VelocityVerlet_next(self, accel, newCoord, newVelocity)
    type(OVelocityVerlet), intent(inout) :: self
    real(dp),intent(in) :: accel(:,:)
    real(dp),intent(out) :: newCoord(:,:)
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

  !!* Rescale the cell parameters and coordinates according to the tensorial
  !!* version of the Berensen barostat (allows cell shape changes if the
  !!* external pressure/stress is non isotropic)
  !!* @param self integrator to rescale
  !!* @param coord atom coordinates to rescale
  !!* @param latVecs lattice vectors to rescale
  !!* @param pressureTensor system stress tensor
  !!* @note the forms of the isotropic and anisotropic Beresdsen barostats in
  !!* the literature are slightly incompatible in their definitions
  subroutine VelocityVerlet_rescale(self,coord,latVecs,pressureTensor)
    type(OVelocityVerlet), intent(inout) :: self
    real(dp),intent(inout)         :: coord(:,:)
    real(dp),intent(inout)         :: latVecs(3,3)
    real(dp),intent(in)            :: pressureTensor(3,3)

    real(dp) :: scale(3,3)
    real(dp) :: scaleIso, Pext, P
    integer  :: ii

    @:ASSERT(self%tBarostat)

    if (self%tIsotropic) then ! isotropic Berendsen, not quite consistent
      ! with anisotropic but its in the literature...
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

  subroutine VelocityVerlet_state(self,fd)
    type(OVelocityVerlet), intent(in) :: self
    integer,intent(in)             :: fd

    if (allocated(self%pThermostat)) then
      call state(self%pThermostat,fd)
    end if

  end subroutine VelocityVerlet_state

end module velocityverlet
