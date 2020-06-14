!--------------------------------------------------------------------------------------------------!
!  DFTB+: general package for performing fast atomistic simulations                                !
!  Copyright (C) 2006 - 2020  DFTB+ developers group                                               !
!                                                                                                  !
!  See the LICENSE file for terms of usage and distribution.                                       !
!--------------------------------------------------------------------------------------------------!

!> Contains enumerated electrostatic solver types.
module dftbp_elstattypes
  implicit none
  private

  public :: elstatTypes


  type :: TElstatTypesEnum

    !> Softened coulombic with gamma function
    integer :: gammaFunc = 0

    !> Poisson equation solver
    integer :: poisson = 1

  end type TElstatTypesEnum


  !> Enumerated electrostatics solver types.
  type(TElstatTypesEnum), parameter :: elstatTypes = TElstatTypesEnum()

end module dftbp_elstattypes
