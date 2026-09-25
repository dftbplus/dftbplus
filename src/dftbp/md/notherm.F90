!--------------------------------------------------------------------------------------------------!
!  DFTB+: general package for performing fast atomistic simulations                                !
!  Copyright (C) 2006 - 2025  DFTB+ developers group                                               !
!                                                                                                  !
!  See the LICENSE file for terms of usage and distribution.                                       !
!--------------------------------------------------------------------------------------------------!

#:include 'common.fypp'

!> No thermostat present, delivers only initial velocities according to the Maxwell-Boltzmann
!! statistics.
module dftbp_md_notherm
  use dftbp_common_accuracy, only : dp, minTemp
  use dftbp_math_ranlux, only : TRanlux
  use dftbp_md_mdcommon, only : MaxwellBoltzmann, rescaleTokT, restFrame, TMDCommon
  use dftbp_md_thermostat, only : TThermostat
  implicit none

  private
  public :: TNoTherm, TNoTherm_init


  !> No thermostat, just initial velocities
  type, extends(TThermostat) :: TNoTherm
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

    procedure :: getInitVelocities => TNoTherm_getInitVelocities
    procedure :: updateVelocities => TNoTherm_updateVelocities
    procedure :: writeState => TNoTherm_writeState

  end type TNoTherm

contains


  !> Creates a NoThermostat instance.
  subroutine TNoTherm_init(this, kT, mass, pRanlux, pMDFrame)

    !> Initialized instance on exit
    type(TNoTherm), intent(out) :: this

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

  end subroutine TNoTherm_init


  !> Returns the initial velocities.
  subroutine TNoTherm_getInitVelocities(this, velocities)

    !> Instance
    class(TNoTherm), intent(inout) :: this

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

  end subroutine TNoTherm_getInitVelocities


  !> Updates velocities (does nothing in this case)
  subroutine TNoTherm_updateVelocities(this, velocities)

    !> Instance
    class(TNoTherm), intent(inout) :: this

    !> Updated velocities on exit.
    real(dp), intent(inout) :: velocities(:,:)

  end subroutine TNoTherm_updateVelocities


  !> Writes internals of thermostat (does nothing in this case)
  subroutine TNoTherm_writeState(this, fd)

    !> Instance
    class(TNoTherm), intent(in) :: this

    !> File unit to write thermostat state out to
    integer, intent(in) :: fd

  end subroutine TNoTherm_writeState

end module dftbp_md_notherm
