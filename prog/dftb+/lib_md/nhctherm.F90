!--------------------------------------------------------------------------------------------------!
!  DFTB+: general package for performing fast atomistic simulations                                !
!  Copyright (C) 2006 - 2020  DFTB+ developers group                                               !
!                                                                                                  !
!  See the LICENSE file for terms of usage and distribution.                                       !
!--------------------------------------------------------------------------------------------------!

#:include 'common.fypp'

!> Nose-Hoover Chain thermostat
!> Based on Martyna et al. Molecular Physics 87 no. 5 1117-1157 (1996).
module dftbp_nhctherm
  use dftbp_assert
  use dftbp_accuracy
  use dftbp_mdcommon
  use dftbp_ranlux
  use dftbp_tempprofile
  use dftbp_energies
  use dftbp_message
  implicit none


  !> Short string for error returns
  character(lc) :: lcTmp

  private

  public :: TNHCThermostat
  public :: init, getInitVelocities, updateVelocities, state


  !> Data for the NHC thermostat
  type TNHCThermostat
    private

    !> Nr. of atoms
    integer :: nAtom

    !> Random number generator
    type(TRanlux), allocatable :: pRanlux

    !> Mass of the atoms
    real(dp), allocatable :: mass(:)

    !> Temperature generator
    type(TTempProfile), pointer :: pTempProfile

    !> coupling strength to friction term
    real(dp) :: couplingParameter

    !> MD Framework.
    type(TMDCommon) :: pMDFrame

    !> MD timestep
    real(dp) :: deltaT

    !> order of NHC operator - eqn 29
    integer :: nyosh

    !> weight coefficients
    real(dp), allocatable :: w(:)

    !> times steps to expand propogator of NHC part of evolution operator
    integer :: nresn

    !> number of thermostat particles in chain
    integer :: nnos

    !> internal chain positions
    real(dp), allocatable :: xnose(:)

    !> internal chain velocities
    real(dp), allocatable :: vnose(:)

    !> internal chain accelerations
    real(dp), allocatable :: gnose(:)
  end type TNHCThermostat


  !> initialise thermostat
  interface init
    module procedure NHC_init
  end interface


  !> initial thermal velocities if needed
  interface getInitVelocities
    module procedure NHC_getInitVelos
  end interface


  !> update velocites acording to the thermostat
  interface updateVelocities
    module procedure NHC_updateVelos
  end interface


  !> write state of the thermostat
  interface state
    module procedure NHC_state
  end interface

contains


  !> Creates an NHC thermostat instance.
  subroutine NHC_init(self, pRanlux, masses, tempProfile, &
      & couplingParameter, pMDFrame, deltaT, npart, nys, nc, &
      & xnose, vnose, gnose)

    !> Initialised instance on exit.
    type(TNHCThermostat), intent(out) :: self

    !> Random generator.
    type(TRanlux), allocatable, intent(inout) :: pRanlux

    !> Masses of the atoms.
    real(dp), intent(in) :: masses(:)

    !> Temperature profile object.
    type(TTempProfile), pointer, intent(in) :: tempProfile

    !> Coupling parameter for the thermostat
    real(dp), intent(in) :: couplingParameter

    !> Molecular dynamics generic framework
    type(TMDCommon), intent(in) :: pMDFrame

    !> MD time step
    real(dp), intent(in) :: deltaT
    integer, intent(in) :: npart
    integer, intent(in) :: nys
    integer, intent(in) :: nc
    real(dp), intent(in), optional :: xnose(:)
    real(dp), intent(in), optional :: vnose(:)
    real(dp), intent(in), optional :: gnose(:)

    @:ASSERT(allocated(pRanlux))
    @:ASSERT(present(xnose).eqv.present(vnose))
    @:ASSERT(present(xnose).eqv.present(gnose))
  #:call ASSERT_CODE
    if (present(xnose)) then
      @:ASSERT(size(xnose)==size(vnose))
      @:ASSERT(size(xnose)==size(gnose))
    end if
  #:endcall ASSERT_CODE

    call move_alloc(pRanlux, self%pRanlux)
    self%nAtom = size(masses)
    allocate(self%mass(self%nAtom))
    self%mass(:) = masses(:)
    self%pTempProfile => tempProfile
    self%couplingParameter = couplingParameter
    self%pMDFrame = pMDFrame
    self%deltaT = deltaT

    ! pg 1124 'For typical simulations, nc can be taken to be one.'
    self%nresn = nc
    if (self%nresn < 1) then
      call error('Nose-Hoover propogation steps must be at least 1.')
    end if

    ! particles in the chain
    self%nnos = npart
    if (self%nnos < 1) then
      call error('Nose-Hoover chains must contain at least one mass.')
    end if

    ! current choice of order
    self%nyosh = nys
    allocate(self%w(self%nyosh))
    select case (self%nyosh)
    case (3)
      self%w(1)=1.0_dp/(2.0_dp - 2.0_dp**(1.0_dp/3.0_dp))
      self%w(2)=1.0_dp-2.0_dp*self%w(1)
      self%w(3)=self%w(1)
    case (5)
      self%w(1) = 1.0_dp / (4.0_dp - 4.0_dp**(1.0_dp/3.0_dp))
      self%w(2:5) = self%w(1)
      self%w(3) = 1.0_dp - 4.0_dp*self%w(1)
    case default
      write (lcTmp, "('Order ',I0,' Nose-Hoover evolution operators are not&
          & available, only order 3 or 5.')") self%nyosh
      call error(lcTmp)
    end select
    allocate(self%xnose(self%nnos))
    allocate(self%vnose(self%nnos))
    allocate(self%gnose(self%nnos))

    ! set intial thermostat positions, velocities and forces
    if (present(xnose)) then
      self%xnose(1:self%nnos)=xnose
      self%vnose(1:self%nnos)=vnose
      self%gnose(1:self%nnos)=gnose
    else
      self%xnose(1:self%nnos)=1.0_dp
      self%vnose(1:self%nnos)=0.0_dp
      self%gnose(1:self%nnos)=0.0_dp
    end if

  end subroutine NHC_init


  !> Returns the initial velocities.
  subroutine NHC_getInitVelos(self, velocities)

    !> NHCThermostat instance.
    type(TNHCThermostat), intent(inout) :: self

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

  end subroutine NHC_getInitVelos


  !> Updates the provided velocities according the current temperature.
  !> routines based on NHCINT from reference
  subroutine NHC_updateVelos(self, velocities)

    !> NHCThermostat instance.
    type(TNHCThermostat), intent(inout) :: self

    !> Updated velocities on exit.
    real(dp), intent(inout) :: velocities(:,:)

    integer :: nnos1, iresn, iyosh, inos
    real(dp) :: qmass(self%nnos)
    real(dp) :: wdti(self%nyosh), wdti2(self%nyosh)
    real(dp) :: wdti4(self%nyosh), wdti8(self%nyosh)
    real(dp) :: scaling, gkt, gnkt, akin, aa

    @:ASSERT(all(shape(velocities) <= (/ 3, self%nAtom /)))

    nnos1=self%nnos+1

    call getTemperature(self%pTempProfile, gkt)
    gnkt=real(self%pMDFrame%Nf,dp)*gkt

    qmass(1)=gnkt/(self%couplingParameter*self%couplingParameter)
    qmass(2:self%nnos)=gkt/(self%couplingParameter*self%couplingParameter)

    ! move to init routine
    wdti(1:self%nyosh)=self%w(1:self%nyosh)*self%deltaT/real(self%nresn,dp)
    wdti2(1:self%nyosh)=wdti(1:self%nyosh)/2.0_dp
    wdti4(1:self%nyosh)=wdti(1:self%nyosh)/4.0_dp
    wdti8(1:self%nyosh)=wdti(1:self%nyosh)/8.0_dp

    ! get the total kinetic energy
    scaling=1.0_dp
    call evalKE(akin,velocities,self%mass)
    akin = 2.0_dp * akin ! Paper defines eqn 1 without 1/2 in K.E. so scale
    ! update the forces
    self%gnose(1) = (akin-gnkt)/qmass(1)
    ! start the multiple time step procedure
    do iresn=1,self%nresn
      do iyosh=1,self%nyosh
        ! update the thermostat velocities
        self%vnose(self%nnos)=self%vnose(self%nnos) &
            & + self%gnose(self%nnos)*wdti4(iyosh)
        do inos=1, self%nnos-1
          aa=exp(-wdti8(iyosh)*self%vnose(nnos1-inos))
          self%vnose(self%nnos-inos)=self%vnose(self%nnos-inos)*aa*aa &
              & +wdti4(iyosh)*self%gnose(self%nnos-inos)*aa
        end do
        ! update the particle velocities
        aa=exp(-wdti2(iyosh)*self%vnose(1))
        scaling=scaling*aa
        ! update the forces
        self%gnose(1)=(scaling*scaling*akin-gnkt)/qmass(1)
        ! update thermostat positions
        do inos= 1, self%nnos
          self%xnose(inos)=self%xnose(inos)+self%vnose(inos)*wdti2(iyosh)
        enddo

        ! update thermostat velocities
        do inos=1, self%nnos-1
          aa=exp(-wdti8(iyosh)*self%vnose(inos+1))
          self%vnose(inos)=self%vnose(inos)*aa*aa &
              & + wdti4(iyosh)*self%gnose(inos)*aa
          self%gnose(inos+1)=(qmass(inos)*self%vnose(inos) &
              & * self%vnose(inos)-gkt)/qmass(inos+1)
        enddo

        self%vnose(self%nnos)=self%vnose(self%nnos) &
            & +self%gnose(self%nnos)*wdti4(iyosh)
      end do

    end do
    ! update particle velocities
    velocities = scaling * velocities

    ! is this needed :
    call restFrame(self%pMDFrame, velocities, self%mass)

  end subroutine NHC_updateVelos


  !> Outputs internals of thermostat
  subroutine NHC_state(self, fd)

    !> instance of thermostat
    type(TNHCThermostat), intent(in) :: self

    !> filehandle to write out to
    integer,intent(in) :: fd

    write(fd,*)'Nose-Hoover chain variables'
    write(fd,*)'x:'
    write(fd,"(3E20.10)")self%xnose
    write(fd,*)'v:'
    write(fd,"(3E20.10)")self%vnose
    write(fd,*)'g:'
    write(fd,"(3E20.10)")self%gnose

  end subroutine NHC_state

end module dftbp_nhctherm
