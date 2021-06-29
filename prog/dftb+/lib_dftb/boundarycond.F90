!--------------------------------------------------------------------------------------------------!
!  DFTB+: general package for performing fast atomistic simulations                                !
!  Copyright (C) 2006 - 2021  DFTB+ developers group                                               !
!                                                                                                  !
!  See the LICENSE file for terms of usage and distribution.                                       !
!--------------------------------------------------------------------------------------------------!

#:include 'common.fypp'

!> Contains geometrical boundary condition information on the calculation
module dftbp_dftb_boundarycond
  use dftbp_common_accuracy, only : dp
  use dftbp_type_commontypes, only : TOrbitals
  implicit none

  private
  public :: zAxis
  public :: boundaryConditions


  !> z direction vector for rotation
  real(dp), parameter :: zAxis(3) = [0.0_dp, 0.0_dp, 1.0_dp]


  !> Enumerator containign possible boundary conditions
  type :: TBoundaryConditionEnum_

    !> Unknown/undefined boundary condition
    integer :: unknown = 0

    !> Finite molecular cluster
    integer :: cluster = 1

    !> Three dimensional infinite periodic boundary conditions
    integer :: pbc3d = 2

    !> Helical boundary conditions
    integer :: helical = 3

  end type TBoundaryConditionEnum_


  !> Actual instance of the boundary condition enumerator
  type(TBoundaryConditionEnum_), parameter :: boundaryConditions = TBoundaryConditionEnum_()


end module dftbp_dftb_boundarycond
