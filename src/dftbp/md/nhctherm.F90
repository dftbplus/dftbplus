!--------------------------------------------------------------------------------------------------!
!  DFTB+: general package for performing fast atomistic simulations                                !
!  Copyright (C) 2006 - 2025  DFTB+ developers group                                               !
!                                                                                                  !
!  See the LICENSE file for terms of usage and distribution.                                       !
!--------------------------------------------------------------------------------------------------!

#:include 'common.fypp'

!> Nose-Hoover Chain thermostat
!>
!> Based on Martyna et al. Molecular Physics 87 no. 5 1117-1157 (1996).
!>
module dftbp_md_nhctherm
  use dftbp_common_accuracy, only : dp, lc, minTemp
  use dftbp_io_message, only : error
  use dftbp_math_ranlux, only : TRanlux
  use dftbp_md_mdcommon, only : evalKE, init, MaxwellBoltzmann, rescaleTokT, restFrame, TMDCommon
  use dftbp_md_tempprofile, only : TTempProfile
  use dftbp_md_thermostat, only : TThermostat
  implicit none

  private
  public :: TNhcThermInput
  public :: TNhcTherm, TNhcTherm_init

  !> Thermostat specific input data for the NHC thermostat
  type :: TNhcThermInput

    !> Coupling parameter for the thermostat
    real(dp) :: coupling

    !> Number of Nose-Hoover particles in the chain
    integer :: nPart

    !> Order of the Nose-Hoover operator (3 or 5)
    integer :: nYs

    !> Number of time steps to expand propagator of NHC part of evolution operator
    integer :: nC

    !> Internal chain positions (only needed when restarting the thermostat)
    real(dp), allocatable :: xnose(:)

    !> Internal chain velocities (only needed when restarting the thermostat)
    real(dp), allocatable :: vnose(:)

    !> Internal chain accelerations (only needed when restarting the thermostat)
    real(dp), allocatable :: gnose(:)

  end type TNhcThermInput


  !> NHC thermostat
  type, extends(TThermostat) :: TNhcTherm
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

    !> times steps to expand propagator of NHC part of evolution operator
    integer :: nresn

    !> number of thermostat particles in chain
    integer :: nnos

    !> internal chain positions (needed to )
    real(dp), allocatable :: xnose(:)

    !> internal chain velocities
    real(dp), allocatable :: vnose(:)

    !> internal chain accelerations
    real(dp), allocatable :: gnose(:)

  contains

    procedure :: getInitVelocities => TNhcTherm_getInitVelocities
    procedure :: updateVelocities => TNhcTherm_updateVelocities
    procedure :: writeState => TNhcTherm_writeState

  end type TNhcTherm

contains


  !> Creates an NHC thermostat instance.
  subroutine TNhcTherm_init(this, input, pRanlux, masses, tempProfile, pMDFrame, deltaT)

    !> Initialised instance on exit
    type(TNhcTherm), intent(out) :: this

    !> Thermostat specific input data
    type(TNhcThermInput), intent(in) :: input

    !> Random generator
    type(TRanlux), allocatable, intent(inout) :: pRanlux

    !> Masses of the atoms
    real(dp), intent(in) :: masses(:)

    !> Temperature profile object
    type(TTempProfile), pointer, intent(in) :: tempProfile

    !> Molecular dynamics generic framework
    type(TMDCommon), intent(in) :: pMDFrame

    !> MD time step
    real(dp), intent(in) :: deltaT

    !> Short string for error returns
    character(lc) :: lcTmp

    @:ASSERT(allocated(pRanlux))
    @:ASSERT(allocated(input%xnose) .eqv. allocated(input%vnose))
    @:ASSERT(allocated(input%xnose) .eqv. allocated(input%gnose))
  #:block DEBUG_CODE
    if (allocated(input%xnose)) then
      @:ASSERT(size(input%xnose) == size(input%vnose))
      @:ASSERT(size(input%xnose) == size(input%gnose))
    end if
  #:endblock DEBUG_CODE

    call move_alloc(pRanlux, this%pRanlux)
    this%nAtom = size(masses)
    this%mass = masses
    this%pTempProfile => tempProfile
    this%couplingParameter = input%coupling
    this%pMDFrame = pMDFrame
    this%deltaT = deltaT

    ! pg 1124 'For typical simulations, nc can be taken to be one.'
    this%nresn = input%nC
    if (this%nresn < 1) then
      call error('Nose-Hoover propagation steps must be at least 1.')
    end if

    ! particles in the chain
    this%nnos = input%nPart
    if (this%nnos < 1) then
      call error('Nose-Hoover chains must contain at least one mass.')
    end if

    ! current choice of order
    this%nyosh = input%nYs
    allocate(this%w(this%nyosh))
    select case (this%nyosh)
    case (3)
      this%w(1)= 1.0_dp / (2.0_dp - 2.0_dp**(1.0_dp / 3.0_dp))
      this%w(2)= 1.0_dp - 2.0_dp * this%w(1)
      this%w(3)=this%w(1)
    case (5)
      this%w(1) = 1.0_dp / (4.0_dp - 4.0_dp**(1.0_dp / 3.0_dp))
      this%w(2:5) = this%w(1)
      this%w(3) = 1.0_dp - 4.0_dp * this%w(1)
    case default
      write (lcTmp, "('Order ',I0,' Nose-Hoover evolution operators are not&
          & available, only order 3 or 5.')") this%nyosh
      call error(lcTmp)
    end select
    allocate(this%xnose(this%nnos))
    allocate(this%vnose(this%nnos))
    allocate(this%gnose(this%nnos))

    ! set initial thermostat positions, velocities and forces
    if (allocated(input%xnose)) then
      this%xnose(1:this%nnos) = input%xnose
      this%vnose(1:this%nnos) = input%vnose
      this%gnose(1:this%nnos) = input%gnose
    else
      this%xnose(1:this%nnos)= 1.0_dp
      this%vnose(1:this%nnos)= 0.0_dp
      this%gnose(1:this%nnos)= 0.0_dp
    end if

  end subroutine TNhcTherm_init


  !> Returns the initial velocities.
  subroutine TNhcTherm_getInitVelocities(this, velocities)

    !> NHCThermostat instance.
    class(TNhcTherm), intent(inout) :: this

    !> Contains the velocities on return.
    real(dp), intent(out) :: velocities(:,:)

    real(dp) :: kT
    integer :: ii

    @:ASSERT(all(shape(velocities) <= [3, this%nAtom]))

    call this%pTempProfile%getTemperature(kT)
    if (kT < minTemp) then
      call error("Nose-Hover thermostat not supported at zero temperature")
    end if
    do ii = 1, this%nAtom
       call MaxwellBoltzmann(velocities(:,ii), this%mass(ii), kT, this%pRanlux)
    end do
    call restFrame(this%pMDFrame, velocities, this%mass)
    call rescaleTokT(this%pMDFrame, velocities, this%mass, kT)

  end subroutine TNhcTherm_getInitVelocities


  !> Updates the provided velocities according the current temperature.
  !>
  !> routines based on NHCINT from reference
  !>
  subroutine TNhcTherm_updateVelocities(this, velocities)

    !> NHCThermostat instance.
    class(TNhcTherm), intent(inout) :: this

    !> Updated velocities on exit.
    real(dp), intent(inout) :: velocities(:,:)

    integer :: nnos1, iresn, iyosh, inos
    real(dp) :: qmass(this%nnos)
    real(dp) :: wdti(this%nyosh), wdti2(this%nyosh)
    real(dp) :: wdti4(this%nyosh), wdti8(this%nyosh)
    real(dp) :: scaling, gkt, gnkt, akin, aa

    @:ASSERT(all(shape(velocities) <= [3, this%nAtom]))

    nnos1 = this%nnos + 1

    call this%pTempProfile%getTemperature(gkt)
    gnkt = real(this%pMDFrame%Nf,dp) * gkt

    qmass(1) = gnkt / (this%couplingParameter * this%couplingParameter)
    qmass(2:this%nnos) = gkt / (this%couplingParameter * this%couplingParameter)

    ! move to init routine
    wdti(1:this%nyosh) = this%w(1:this%nyosh) * this%deltaT / real(this%nresn,dp)
    wdti2(1:this%nyosh) = wdti(1:this%nyosh) / 2.0_dp
    wdti4(1:this%nyosh) = wdti(1:this%nyosh) / 4.0_dp
    wdti8(1:this%nyosh) = wdti(1:this%nyosh) / 8.0_dp

    ! get the total kinetic energy
    scaling = 1.0_dp
    call evalKE(akin, velocities, this%mass)
    akin = 2.0_dp * akin ! Paper defines eqn 1 without 1/2 in K.E. so scale
    ! update the forces
    this%gnose(1) = (akin - gnkt) / qmass(1)
    ! start the multiple time step procedure
    do iresn = 1, this%nresn
      do iyosh = 1, this%nyosh
        ! update the thermostat velocities
        this%vnose(this%nnos) = this%vnose(this%nnos) + this%gnose(this%nnos) * wdti4(iyosh)
        do inos = 1, this%nnos - 1
          aa=exp(-wdti8(iyosh) * this%vnose(nnos1 - inos))
          this%vnose(this%nnos - inos) = this%vnose(this%nnos - inos) * aa * aa&
              & + wdti4(iyosh) * this%gnose(this%nnos - inos) * aa
        end do
        ! update the particle velocities
        aa = exp(-wdti2(iyosh) * this%vnose(1))
        scaling = scaling * aa
        ! update the forces
        this%gnose(1)=(scaling * scaling * akin - gnkt) / qmass(1)
        ! update thermostat positions
        do inos = 1, this%nnos
          this%xnose(inos) = this%xnose(inos) + this%vnose(inos) * wdti2(iyosh)
        enddo

        ! update thermostat velocities
        do inos = 1, this%nnos - 1
          aa = exp(-wdti8(iyosh) * this%vnose(inos + 1))
          this%vnose(inos) = this%vnose(inos) * aa * aa + wdti4(iyosh) * this%gnose(inos) * aa
          this%gnose(inos + 1) = (qmass(inos) * this%vnose(inos) * this%vnose(inos) - gkt)&
              & / qmass(inos+1)
        enddo

        this%vnose(this%nnos) = this%vnose(this%nnos) + this%gnose(this%nnos) * wdti4(iyosh)
      end do

    end do
    ! update particle velocities
    velocities = scaling * velocities

    ! is this needed :
    call restFrame(this%pMDFrame, velocities, this%mass)

  end subroutine TNhcTherm_updateVelocities


  !> Writes internals of thermostat
  subroutine TNhcTherm_writeState(this, fd)

    !> Instance
    class(TNhcTherm), intent(in) :: this

    !> File handle to write out to
    integer, intent(in) :: fd

    write(fd,"(a)") 'Nose-Hoover chain variables'
    write(fd,"(a)") 'x:'
    write(fd,"(3e20.10)") this%xnose
    write(fd,"(a)") 'v:'
    write(fd,"(3e20.10)") this%vnose
    write(fd,"(a)") 'g:'
    write(fd,"(3e20.10)") this%gnose

  end subroutine TNhcTherm_writeState

end module dftbp_md_nhctherm
