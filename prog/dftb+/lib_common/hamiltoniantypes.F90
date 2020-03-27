!--------------------------------------------------------------------------------------------------!
!  DFTB+: general package for performing fast atomistic simulations                                !
!  Copyright (C) 2006 - 2020  DFTB+ developers group                                               !
!                                                                                                  !
!  See the LICENSE file for terms of usage and distribution.                                       !
!--------------------------------------------------------------------------------------------------!

module dftbp_hamiltoniantypes
  implicit none
  private

  public :: hamiltonianTypes

  !> Namespace for possible hamiltonian models
  type :: THamiltonianTypesEnum

    ! Hamiltonian models

    !> Dummy none
    integer :: none = 0

    !> DFTB
    integer :: dftb = 1

    !> XTB
    integer :: xtb = 2

  end type THamiltonianTypesEnum

  !> Actual values for hamiltonianTypes.
  type(THamiltonianTypesEnum), parameter :: hamiltonianTypes = THamiltonianTypesEnum()

end module dftbp_hamiltoniantypes
