!--------------------------------------------------------------------------------------------------!
!  DFTB+: general package for performing fast atomistic simulations                                !
!  Copyright (C) 2018  DFTB+ developers group                                                      !
!                                                                                                  !
!  See the LICENSE file for terms of usage and distribution.                                       !
!--------------------------------------------------------------------------------------------------!

#:include 'common.fypp'

!> Contains geometrical boundary condition information on the calculation
module dftbp_boundarycond
  use dftbp_angmomentum, only : rotateZ
  use dftbp_quaternions, only : rotate3
  use dftbp_constants, only : pi
  use dftbp_accuracy, only : dp
  use dftbp_commontypes, only : TOrbitals
  implicit none

  private
  public :: zAxis, boundaryCondition

  !> z direction vector for rotation
  real(dp), parameter :: zAxis(3) = [0.0_dp,0.0_dp,1.0_dp]


  !> Possible boundary conditions
  type :: TBoundaryConditionEnum

    !> Finite molecular cluster
    integer :: cluster = 0

    !> One dimensional infinite periodic boundary conditions, finite in other directions
    integer :: pbd1d = 1

    !> Three dimensional infinite periodic boundary conditions
    integer :: pbc3d = 3

    !> z-axis oriented helix
    integer :: helical = 4

  end type TBoundaryConditionEnum

  !> Actual instance of the boundary condition enumerator
  type(TBoundaryConditionEnum), parameter :: boundaryCondition = TBoundaryConditionEnum()

end module dftbp_boundarycond
