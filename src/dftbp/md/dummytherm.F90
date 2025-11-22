!--------------------------------------------------------------------------------------------------!
!  DFTB+: general package for performing fast atomistic simulations                                !
!  Copyright (C) 2006 - 2025  DFTB+ developers group                                               !
!                                                                                                  !
!  See the LICENSE file for terms of usage and distribution.                                       !
!--------------------------------------------------------------------------------------------------!

#:include 'common.fypp'

!> Dummy thermostat, delivers only initial velocities according to the Maxwell-Boltzmann statistics.
module dftbp_md_dummytherm
  use dftbp_common_accuracy, only : dp, minTemp
  use dftbp_math_ranlux, only : TRanlux
  use dftbp_md_mdcommon, only : MaxwellBoltzmann, rescaleTokT, restFrame, TMDCommon
  use dftbp_md_thermostat, only : TThermostat
  implicit none

  private
  public :: TDummyTherm, TDummyTherm_init


  !> Dummy thermostat
  type, extends(TThermostat) :: TDummyTherm
    private

    !> Nr. of atoms
    integer :: nAtom

    !> Temperature
    real(dp) :: kT

    !> Mass of the atoms
    real(dp), allocatable :: mass(:)

    !> Random number generator.
    type(TRanlux), allocatable :: pRanlux

    !> MD Framwork
    type(TMDCommon) :: pMDFrame

  contains

    procedure :: getInitVelocities => TDummyTherm_getInitVelocities
    procedure :: updateVelocities => TDummyTherm_updateVelocities
    procedure :: writeState => TDummyTherm_writeState

  end type TDummyTherm

contains


  !> Creates a DummyThermostat instance.
  subroutine TDummyTherm_init(this, kT, mass, pRanlux, pMDFrame)

    !> Initialized instance on exit
    type(TDummyTherm), intent(out) :: this

    !> Temperature
    real(dp), intent(in) :: kT

    !> Mass of the atoms
    real(dp), intent(in) :: mass(:)

    !> Random generator
    type(TRanlux), allocatable, intent(inout) :: pRanlux

    !> thermostat object
    type(TMDCommon), intent(in) :: pMDFrame

    this%kT = kT
    this%nAtom = size(mass)
    allocate(this%mass(this%nAtom))
    this%mass = mass(:)
    call move_alloc(pRanlux, this%pRanlux)
    this%pMDFrame = pMDFrame

  end subroutine TDummyTherm_init


  !> Returns the initial velocities.
  subroutine TDummyTherm_getInitVelocities(this, velocities)

    !> Instance
    class(TDummyTherm), intent(inout) :: this

    !> Velocities on return.
    real(dp), intent(out) :: velocities(:,:)

    integer :: ii

    @:ASSERT(all(shape(velocities) >= [3, this%nAtom]))

    if (this%kT > minTemp) then
      do ii = 1, this%nAtom
        call MaxwellBoltzmann(velocities(:,ii), this%mass(ii), this%kT, this%pRanlux)
      end do
      call restFrame(this%pMDFrame, velocities(:,:), this%mass)
      call rescaleTokT(this%pMDFrame, velocities(:,:), this%mass, this%kT)
    else
      velocities(:,:) = 0.0_dp
    end if

  end subroutine TDummyTherm_getInitVelocities


  !> Updates velocities (does nothing in this case)
  subroutine TDummyTherm_updateVelocities(this, velocities)

    !> Instance
    class(TDummyTherm), intent(inout) :: this

    !> Updated velocities on exit.
    real(dp), intent(inout) :: velocities(:,:)

  end subroutine TDummyTherm_updateVelocities


  !> Writes internals of thermostat (does nothing in this case)
  subroutine TDummyTherm_writeState(this, fd)

    !> instance of thermostat
    class(TDummyTherm), intent(in) :: this

    !> File handle to write state out to
    integer, intent(in) :: fd

  end subroutine TDummyTherm_writeState

end module dftbp_md_dummytherm
