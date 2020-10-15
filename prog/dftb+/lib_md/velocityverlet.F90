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
  public :: init, next, rescale, reset, state

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

    !> Is this initialised?
    logical :: tInitialised = .false.

  end type TVelocityVerlet


  !> initialise MD
  interface init
    module procedure VelocityVerlet_init
    module procedure VV_thermostats_pressure
    module procedure VV_velocities_pressure
  end interface


  !> next geometry step
  interface next
    module procedure VelocityVerlet_next
  end interface

  !> Resets the coordinates and velocities
  interface reset
    module procedure VelocityVerlet_reset
  end interface reset

  !> Adjust velocities
  interface rescale
    module procedure VelocityVerlet_rescale
  end interface


  !> write state of the integrator
  interface state
    module procedure VelocityVerlet_state
  end interface

contains


  !> Creates a VelocityVerlet object from the thermostat settings and optional starting velocity
  subroutine VelocityVerlet_init(this, deltaT, positions, pThermostat, velocities, tSetVelocities,&
      & tHalfVelocities)

    !> Initialised object on exit.
    type(TVelocityVerlet), intent(out) :: this

    !> Integration time step.
    real(dp), intent(in) :: deltaT

    !> Position of the atoms.
    real(dp), intent(in) :: positions(:,:)

    !> Thermostat if needed.
    type(TThermostat), allocatable, intent(inout) :: pThermostat

    !> On input, if tHalfVelocities these are the t=-.5 velocities, but ignored if false. On output
    !> these are the internal velocities, either at current time or t=-.5 depending on setting of
    !> tHalfVelocities if this is allocated
    real(dp), intent(inout), allocatable :: velocities(:,:)

    !> Should the velocities be read from the velocities input (T) or set from the thermostat (F)?
    logical, intent(in) :: tSetVelocities

    !> This indicates if the routine is setting the t-.5 velocities internally, otherwise they need
    !> to be regenerated later.
    logical, intent(in) :: tHalfVelocities

    @:ASSERT(.not.this%tInitialised)
    @:ASSERT(size(positions, dim=1) == 3)
    @:ASSERT(allocated(pThermostat))

    this%nAtom = size(positions, dim=2)
    allocate(this%velocities(3, this%nAtom))
    allocate(this%positions(3, this%nAtom))

    this%deltaT = deltaT
    this%positions(:,:) = positions(:,:)
    call move_alloc(pThermostat, this%pThermostat)

    if (tSetVelocities) then
      if (.not.allocated(velocities)) then
        call error("Velocities must be allocated to initialise the VV driver")
      end if
      this%velocities(:,:) = velocities
    else
      call getInitVelocities(this%pThermostat, this%velocities)
    end if

    if (allocated(velocities)) then
      velocities(:,:) = this%velocities
    end if

    this%vHalfPresent = tHalfVelocities

    this%tBarostat = .false.

    this%tInitialised = .true.

  end subroutine VelocityVerlet_Init


  !> Creates a VelocityVerlet object from the thermostat settings and isotropic pressure
  subroutine VV_thermostats_pressure(this, deltaT, positions, pThermostat, &
      & Barostat, Pressure, tIsotropic)

    !> Initialised object on exit.
    type(TVelocityVerlet), intent(out) :: this

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

    @:ASSERT(.not.this%tInitialised)
    @:ASSERT(size(positions, dim=1) == 3)

    this%nAtom = size(positions, dim=2)
    allocate(this%velocities(3, this%nAtom))
    allocate(this%positions(3, this%nAtom))

    this%deltaT = deltaT
    this%positions(:,:) = positions(:,:)
    call move_alloc(pThermostat, this%pThermostat)

    call getInitVelocities(this%pThermostat, this%velocities)

    this%vHalfPresent = .true. ! yes we have the t-.5 velocities

    this%tBarostat = .true.
    this%BarostatStrength = Barostat
    this%Pressure = 0.0_dp
    do ii = 1, 3
      this%Pressure(ii,ii) = pressure
    end do

    this%tIsotropic = tIsotropic

    this%tInitialised = .true.

  end subroutine VV_thermostats_pressure


  !> Creates a VelocityVerlet object from given external velocities for the t-th time step, this
  !> means later we have to reconstruct the Vel. Verlet t+.5 velocities and barostat isotropic
  !> pressure
  subroutine VV_velocities_pressure(this, deltaT, positions, pThermostat, &
      & velocities, Barostat, Pressure, tIsotropic)

    !> Initialised object on exit.
    type(TVelocityVerlet), intent(out) :: this

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

    @:ASSERT(.not.this%tInitialised)
    @:ASSERT(size(positions, dim=1) == 3)

    this%nAtom = size(positions, dim=2)
    allocate(this%velocities(3, this%nAtom))
    allocate(this%positions(3, this%nAtom))

    this%deltaT = deltaT
    this%positions(:,:) = positions(:,:)
    call move_alloc(pThermostat, this%pThermostat)

    this%velocities(:,:) = velocities(:,:)

    ! assumes the V read in corresponds to the current coordinates, so we should reconstruct the
    ! t+.5 velocities when possible once forces are available for the coordinates
    this%vHalfPresent = .false.

    this%tBarostat = .true.
    this%BarostatStrength = Barostat

    this%Pressure = 0.0_dp
    do ii = 1, 3
      this%Pressure(ii,ii) = pressure
    end do
    this%tIsotropic = tIsotropic

    this%tInitialised = .true.

  end subroutine VV_velocities_pressure


  !> Takes a timestep for the MD integrator, optionally with a thermostat.
  !> Due to the way the velocity Verlet is split, the returned velocity is for 1 complete MD step
  !> behind the returned positions so print positions, then call next and then print velocities to
  !> get agreement between the positions and velocities.
  subroutine VelocityVerlet_next(this, accel, newCoord, newVelocity)

    !> Integrator to propogate
    type(TVelocityVerlet), intent(inout) :: this

    !> Accelerations.
    real(dp),intent(in) :: accel(:,:)

    !> Displaced coordinates
    real(dp),intent(out) :: newCoord(:,:)

    !> Velocity of displaced coords
    real(dp),intent(out) :: newVelocity(:,:)

    @:ASSERT(this%tInitialised)

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

    if (this%vHalfPresent) then
      newVelocity(:,:) = this%velocities(:,:) &
          & + 0.5_dp * accel(:,:) * this%deltaT
    else
      newVelocity(:,:) = this%velocities(:,:)
      this%vHalfPresent=.true.
    end if

    this%velocities(:,:) = newVelocity(:,:) + 0.5_dp * accel(:,:) * this%deltaT
    newCoord(:,:) = this%positions(:,:) + this%velocities(:,:) * this%deltaT
    this%positions(:,:) = newCoord(:,:)

    if (allocated(this%pThermostat)) then
      call updateVelocities(this%pThermostat, this%velocities)
    end if

  end subroutine VelocityVerlet_next


  !> Rescale the cell parameters and coordinates according to the tensorial version of the Berensen
  !> barostat (allows cell shape changes if the external pressure/stress is non isotropic)
  !> The forms of the isotropic and anisotropic Beresdsen barostats in the literature are slightly
  !> incompatible in their definitions
  subroutine VelocityVerlet_rescale(this,coord,latVecs,pressureTensor)

    !> Integrator to rescale
    type(TVelocityVerlet), intent(inout) :: this

    !> Atom coordinates to rescale
    real(dp),intent(inout) :: coord(:,:)

    !> Lattice vectors to rescale
    real(dp),intent(inout) :: latVecs(3,3)

    !> System stress tensor
    real(dp),intent(in) :: pressureTensor(3,3)

    real(dp) :: scale(3,3)
    real(dp) :: scaleIso, Pext, P
    integer :: ii

    @:ASSERT(this%tInitialised)
    @:ASSERT(this%tBarostat)

    ! isotropic Berendsen, not quite consistent with anisotropic but its in the literature...
    if (this%tIsotropic) then
      Pext = 0.0_dp
      P = 0.0_dp
      do ii = 1, 3
        Pext = this%Pressure(ii,ii) / 3.0_dp
        P = P + pressureTensor(ii,ii) / 3.0_dp
      end do
      scaleIso = (1.0_dp - this%BarostatStrength*(Pext - P))**(1.0_dp/3.0_dp)
      this%positions(:,:) = this%positions(:,:) * scaleIso
      coord(:,:) = coord(:,:) * scaleIso
      latVecs(:,:) = latVecs(:,:) * scaleIso
    else
      scale = 0.0_dp
      do ii = 1, 3
        scale(ii,ii) = 1.0_dp
      end do
      scale = scale - this%BarostatStrength*(this%Pressure-pressureTensor)
      do ii = 1, this%nAtom
        this%positions(:,ii) = matmul(this%positions(:,ii),scale)
        coord(:,ii) = matmul(coord(:,ii),scale)
      end do
      latVecs(:,:) = matmul(latVecs(:,:),scale)
    end if

  end subroutine VelocityVerlet_rescale


  !> Outputs internals of MD integrator
  subroutine VelocityVerlet_state(this, fd, velocities)

    !> instance of integrator
    type(TVelocityVerlet), intent(in) :: this

    !> filehandle to write out to
    integer,intent(in), optional :: fd

    real(dp), intent(out), optional :: velocities(:,:)

    @:ASSERT(this%tInitialised)

    if (present(fd)) then
      if (allocated(this%pThermostat)) then
        call state(this%pThermostat,fd)
      end if
    end if

    if (present(velocities)) then
      velocities(:,:) = this%velocities(:,:)
    end if

  end subroutine VelocityVerlet_state


  !> replaces the positions and velocities in a running VV instance
  subroutine VelocityVerlet_reset(this, positions, velocities, tHalfVelocities)

    !> Instance.
    type(TVelocityVerlet), intent(inout) :: this

    !> New position of the atoms.
    real(dp), intent(in) :: positions(:,:)

    !> On input, if tHalfVelocities these are the t=-.5 velocities, but ignored if false. On output
    !> these are the internal velocities, either at current time or t=-.5 depending on setting of
    !> tHalfVelocities if this is allocated
    real(dp), intent(inout) :: velocities(:,:)

    !> This indicates if the routine is setting the t-.5 velocities internally, otherwise they need
    !> to be regenerated later.
    logical, intent(in) :: tHalfVelocities

    @:ASSERT(this%tInitialised)

    this%positions(:,:) = positions(:,:)

    this%velocities(:,:) = velocities

    this%vHalfPresent = tHalfVelocities

  end subroutine VelocityVerlet_reset

end module dftbp_velocityverlet
