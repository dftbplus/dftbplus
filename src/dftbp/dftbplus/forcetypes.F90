!--------------------------------------------------------------------------------------------------!
!  DFTB+: general package for performing fast atomistic simulations                                !
!  Copyright (C) 2006 - 2025  DFTB+ developers group                                               !
!                                                                                                  !
!  See the LICENSE file for terms of usage and distribution.                                       !
!--------------------------------------------------------------------------------------------------!

!> Contains enumeration of force calculation types
module dftbp_dftbplus_forcetypes
  implicit none

  private
  public :: forceTypes


  type :: TForceTypesEnum

    !> Uninitialised
    integer :: none = 0

    !> Conventional forces
    integer :: orig = 1

    !> convergence corrected at 0 temperature
    integer :: dynamicT0 = 2

    !> convergence corrected at finite temperature
    integer :: dynamicTFinite = 3

  end type TForceTypesEnum

  !> Container for enumerated force calculation types.
  type(TForceTypesEnum), parameter :: forceTypes = TForceTypesEnum()

end module dftbp_dftbplus_forcetypes
