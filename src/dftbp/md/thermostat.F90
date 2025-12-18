!--------------------------------------------------------------------------------------------------!
!  DFTB+: general package for performing fast atomistic simulations                                !
!  Copyright (C) 2006 - 2025  DFTB+ developers group                                               !
!                                                                                                  !
!  See the LICENSE file for terms of usage and distribution.                                       !
!--------------------------------------------------------------------------------------------------!

!> Contains the abstract interface that can provide any of the thermostats
module dftbp_md_thermostat
  use dftbp_common_accuracy, only : dp
  implicit none

  private
  public :: TThermostat

  !> Abstract interface for thermostats
  type, abstract :: TThermostat
  contains
    procedure(TThermostat_getInitVelocities), deferred :: getInitVelocities
    procedure(TThermostat_updateVelocities), deferred :: updateVelocities
    procedure(TThermostat_writeState), deferred :: writeState
  end type TThermostat


  abstract interface

    !> Returns the initial velocities
    subroutine TThermostat_getInitVelocities(this, velocities)
      import :: TThermostat, dp
      implicit none

      !> Instance
      class(TThermostat), intent(inout) :: this

      !> Velocities on exit
      real(dp), intent(out) :: velocities(:,:)

    end subroutine TThermostat_getInitVelocities


    !> Updates the velocities
    subroutine TThermostat_updateVelocities(this, velocities)
      import :: TThermostat, dp
      implicit none

      !> Instance
      class(TThermostat), intent(inout) :: this

      !> Updated velocities on exit
      real(dp), intent(inout) :: velocities(:,:)

    end subroutine TThermostat_updateVelocities


    !> Write internal state of the thermostat
    subroutine TThermostat_writeState(this, fd)
      import :: TThermostat
      implicit none

      !> Instance
      class(TThermostat), intent(in) :: this

      !> File handle to write state out to
      integer, intent(in) :: fd

    end subroutine TThermostat_writeState

  end interface

end module dftbp_md_thermostat
