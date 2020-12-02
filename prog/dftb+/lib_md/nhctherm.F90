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
  use dftbp_message
  implicit none

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
  subroutine NHC_init(this, pRanlux, masses, tempProfile, &
      & couplingParameter, pMDFrame, deltaT, npart, nys, nc, &
      & xnose, vnose, gnose)

    !> Initialised instance on exit.
    type(TNHCThermostat), intent(out) :: this

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

    !> Short string for error returns
    character(lc) :: lcTmp

    integer, intent(in) :: npart
    integer, intent(in) :: nys
    integer, intent(in) :: nc
    real(dp), intent(in), optional :: xnose(:)
    real(dp), intent(in), optional :: vnose(:)
    real(dp), intent(in), optional :: gnose(:)

    @:ASSERT(allocated(pRanlux))
    @:ASSERT(present(xnose).eqv.present(vnose))
    @:ASSERT(present(xnose).eqv.present(gnose))
  #:block DEBUG_CODE
    if (present(xnose)) then
      @:ASSERT(size(xnose)==size(vnose))
      @:ASSERT(size(xnose)==size(gnose))
    end if
  #:endblock DEBUG_CODE

    call move_alloc(pRanlux, this%pRanlux)
    this%nAtom = size(masses)
    allocate(this%mass(this%nAtom))
    this%mass(:) = masses(:)
    this%pTempProfile => tempProfile
    this%couplingParameter = couplingParameter
    this%pMDFrame = pMDFrame
    this%deltaT = deltaT

    ! pg 1124 'For typical simulations, nc can be taken to be one.'
    this%nresn = nc
    if (this%nresn < 1) then
      call error('Nose-Hoover propogation steps must be at least 1.')
    end if

    ! particles in the chain
    this%nnos = npart
    if (this%nnos < 1) then
      call error('Nose-Hoover chains must contain at least one mass.')
    end if

    ! current choice of order
    this%nyosh = nys
    allocate(this%w(this%nyosh))
    select case (this%nyosh)
    case (3)
      this%w(1)=1.0_dp/(2.0_dp - 2.0_dp**(1.0_dp/3.0_dp))
      this%w(2)=1.0_dp-2.0_dp*this%w(1)
      this%w(3)=this%w(1)
    case (5)
      this%w(1) = 1.0_dp / (4.0_dp - 4.0_dp**(1.0_dp/3.0_dp))
      this%w(2:5) = this%w(1)
      this%w(3) = 1.0_dp - 4.0_dp*this%w(1)
    case default
      write (lcTmp, "('Order ',I0,' Nose-Hoover evolution operators are not&
          & available, only order 3 or 5.')") this%nyosh
      call error(lcTmp)
    end select
    allocate(this%xnose(this%nnos))
    allocate(this%vnose(this%nnos))
    allocate(this%gnose(this%nnos))

    ! set intial thermostat positions, velocities and forces
    if (present(xnose)) then
      this%xnose(1:this%nnos)=xnose
      this%vnose(1:this%nnos)=vnose
      this%gnose(1:this%nnos)=gnose
    else
      this%xnose(1:this%nnos)=1.0_dp
      this%vnose(1:this%nnos)=0.0_dp
      this%gnose(1:this%nnos)=0.0_dp
    end if

  end subroutine NHC_init


  !> Returns the initial velocities.
  subroutine NHC_getInitVelos(this, velocities)

    !> NHCThermostat instance.
    type(TNHCThermostat), intent(inout) :: this

    !> Contains the velocities on return.
    real(dp), intent(out) :: velocities(:,:)

    real(dp) :: kT
    integer :: ii

    @:ASSERT(all(shape(velocities) <= (/ 3, this%nAtom /)))

    call this%pTempProfile%getTemperature(kT)
    if (kT < minTemp) then
      call error("Nose-Hover thermostat not supported at zero temperature")
    end if
    do ii = 1, this%nAtom
       call MaxwellBoltzmann(velocities(:,ii), this%mass(ii), kT, this%pRanlux)
    end do
    call restFrame(this%pMDFrame, velocities, this%mass)
    call rescaleTokT(this%pMDFrame, velocities, this%mass, kT)

  end subroutine NHC_getInitVelos


  !> Updates the provided velocities according the current temperature.
  !> routines based on NHCINT from reference
  subroutine NHC_updateVelos(this, velocities)

    !> NHCThermostat instance.
    type(TNHCThermostat), intent(inout) :: this

    !> Updated velocities on exit.
    real(dp), intent(inout) :: velocities(:,:)

    integer :: nnos1, iresn, iyosh, inos
    real(dp) :: qmass(this%nnos)
    real(dp) :: wdti(this%nyosh), wdti2(this%nyosh)
    real(dp) :: wdti4(this%nyosh), wdti8(this%nyosh)
    real(dp) :: scaling, gkt, gnkt, akin, aa

    @:ASSERT(all(shape(velocities) <= (/ 3, this%nAtom /)))

    nnos1=this%nnos+1

    call this%pTempProfile%getTemperature(gkt)
    gnkt=real(this%pMDFrame%Nf,dp)*gkt

    qmass(1)=gnkt/(this%couplingParameter*this%couplingParameter)
    qmass(2:this%nnos)=gkt/(this%couplingParameter*this%couplingParameter)

    ! move to init routine
    wdti(1:this%nyosh)=this%w(1:this%nyosh)*this%deltaT/real(this%nresn,dp)
    wdti2(1:this%nyosh)=wdti(1:this%nyosh)/2.0_dp
    wdti4(1:this%nyosh)=wdti(1:this%nyosh)/4.0_dp
    wdti8(1:this%nyosh)=wdti(1:this%nyosh)/8.0_dp

    ! get the total kinetic energy
    scaling=1.0_dp
    call evalKE(akin,velocities,this%mass)
    akin = 2.0_dp * akin ! Paper defines eqn 1 without 1/2 in K.E. so scale
    ! update the forces
    this%gnose(1) = (akin-gnkt)/qmass(1)
    ! start the multiple time step procedure
    do iresn=1,this%nresn
      do iyosh=1,this%nyosh
        ! update the thermostat velocities
        this%vnose(this%nnos)=this%vnose(this%nnos) &
            & + this%gnose(this%nnos)*wdti4(iyosh)
        do inos=1, this%nnos-1
          aa=exp(-wdti8(iyosh)*this%vnose(nnos1-inos))
          this%vnose(this%nnos-inos)=this%vnose(this%nnos-inos)*aa*aa &
              & +wdti4(iyosh)*this%gnose(this%nnos-inos)*aa
        end do
        ! update the particle velocities
        aa=exp(-wdti2(iyosh)*this%vnose(1))
        scaling=scaling*aa
        ! update the forces
        this%gnose(1)=(scaling*scaling*akin-gnkt)/qmass(1)
        ! update thermostat positions
        do inos= 1, this%nnos
          this%xnose(inos)=this%xnose(inos)+this%vnose(inos)*wdti2(iyosh)
        enddo

        ! update thermostat velocities
        do inos=1, this%nnos-1
          aa=exp(-wdti8(iyosh)*this%vnose(inos+1))
          this%vnose(inos)=this%vnose(inos)*aa*aa &
              & + wdti4(iyosh)*this%gnose(inos)*aa
          this%gnose(inos+1)=(qmass(inos)*this%vnose(inos) &
              & * this%vnose(inos)-gkt)/qmass(inos+1)
        enddo

        this%vnose(this%nnos)=this%vnose(this%nnos) &
            & +this%gnose(this%nnos)*wdti4(iyosh)
      end do

    end do
    ! update particle velocities
    velocities = scaling * velocities

    ! is this needed :
    call restFrame(this%pMDFrame, velocities, this%mass)

  end subroutine NHC_updateVelos


  !> Outputs internals of thermostat
  subroutine NHC_state(this, fd)

    !> instance of thermostat
    type(TNHCThermostat), intent(in) :: this

    !> filehandle to write out to
    integer,intent(in) :: fd

    write(fd,*)'Nose-Hoover chain variables'
    write(fd,*)'x:'
    write(fd,"(3E20.10)")this%xnose
    write(fd,*)'v:'
    write(fd,"(3E20.10)")this%vnose
    write(fd,*)'g:'
    write(fd,"(3E20.10)")this%gnose

  end subroutine NHC_state

end module dftbp_nhctherm
