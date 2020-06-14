!--------------------------------------------------------------------------------------------------!
!  DFTB+: general package for performing fast atomistic simulations                                !
!  Copyright (C) 2006 - 2020  DFTB+ developers group                                               !
!                                                                                                  !
!  See the LICENSE file for terms of usage and distribution.                                       !
!--------------------------------------------------------------------------------------------------!

!> Contains enumeration of force calculation types
module dftbp_forcetypes
  implicit none
  private

  public :: forceTypes


  type :: TForceTypesEnum

    !> Conventional forces
    integer :: orig = 0

    !> convergence corrected at 0 temperature
    integer :: dynamicT0 = 1

    !> convergence corrected at finite temperature
    integer :: dynamicTFinite = 2

  end type TForceTypesEnum

  !> Container for enumerated force calculation types.
  type(TForceTypesEnum), parameter :: forceTypes = TForceTypesEnum()

end module dftbp_forcetypes
