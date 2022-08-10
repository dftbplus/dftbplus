!--------------------------------------------------------------------------------------------------!
!  DFTB+: general package for performing fast atomistic simulations                                !
!  Copyright (C) 2006 - 2022  DFTB+ developers group                                               !
!                                                                                                  !
!  See the LICENSE file for terms of usage and distribution.                                       !
!--------------------------------------------------------------------------------------------------!

module dftbp_common_modeltypes
  implicit none

  private
  public :: modelTypes

  !> Namespace for possible models
  type :: TModelTypesEnum

    ! models

    !> Dummy none
    integer :: none = 0

    !> DFTB
    integer :: dftb = 1

    !> XTB
    integer :: xtb = 2

    !> Externally specified via API
    integer :: externalmodel = 3

  end type TModelTypesEnum

  !> Actual values for modelTypes.
  type(TModelTypesEnum), parameter :: modelTypes = TModelTypesEnum()

end module dftbp_common_modeltypes
