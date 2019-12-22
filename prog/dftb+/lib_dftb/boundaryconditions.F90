!--------------------------------------------------------------------------------------------------!
!  DFTB+: general package for performing fast atomistic simulations                                !
!  Copyright (C) 2018  DFTB+ developers group                                                      !
!                                                                                                  !
!  See the LICENSE file for terms of usage and distribution.                                       !
!--------------------------------------------------------------------------------------------------!

#:include 'common.fypp'

!> Contains geometrical boundary condition information on the calculation
module dftbp_boundaryconditions
  use dftbp_angmomentum, only : rotateZ
  use dftbp_quaternions, only : rotate3
  use dftbp_accuracy, only : dp
  implicit none

  private
  public :: TBoundaryConditions, boundaryTypes, BoundaryConditions_init

  type :: TBoundaryTypesEnum

    !> real space cluster/molecular
    integer :: cluster = 0

    !> periodic 3D
    integer :: periodic3D = 1

    !> periodic 1D
    integer :: periodic1D = 2

    !> objective single helical operation (effectively periodic1D + twisting)
    integer :: helical = 3

    !> objective two helical operations
    integer :: helical2Op = 4

    !> Rotational symmetry in xy plane of order N
    integer :: rotational = 5

  end type TBoundaryTypesEnum


  !> Container for enumerated geometric boundary types
  type(TBoundaryTypesEnum), parameter :: boundaryTypes = TBoundaryTypesEnum()


  !> Type for geometric boundary conditions
  type :: TBoundaryConditions

    private

    !> Geometry type of the system
    integer :: boundaryType

    !> Lattice vectors for periodic cases
    real(dp), allocatable :: latVec(:,:)

    !> Objective helical boundary angle
    real(dp) :: theta

    !> Objective helical translation
    real(dp) :: T

    !> Objective helical boundary angle
    real(dp) :: theta2

    !> Objective helical translation
    real(dp) :: T2

    !> Objective operation rotational order
    integer :: Nrot

  contains

    procedure :: foldOrbsToCell => foldOrbsToCell
    procedure :: foldOrbsFromCell  => foldOrbsFromCell
    procedure :: foldCoordToCell   => foldCoordToCell
    procedure :: foldCoordFromCell => foldCoordFromCell

  end type TBoundaryConditions

contains

  !> Initialise the type of boundary condition on the geometry
  subroutine BoundaryConditions_init()


  end subroutine BoundaryConditions_init


  pure function foldCoordToCell(xIn) result xOut

    real(dp), intent(in) :: xIn(3)

    real(dp) :: xOut(3)

  end function foldCoordToCell\


  pure function foldCoordFromCell(xIn, rCellVec) result xOut

    real(dp), intent(in) :: xIn(3)

    real(dp) :: xOut(3)

  end function foldCoordFromCell


  pure function foldOrbsToCell(orbsIn) result orbsOut

    real(dp), intent(in) :: orbsIn(:,:)

    real(dp) :: orbsOut(size(orbsIn,dim=1),size(orbsIn,dim=2))

  end function foldOrbsToCell


  pure function foldOrbsFromCell(orbsIn) result orbsOut

    real(dp), intent(in) :: orbsIn(:,:)

    real(dp) :: orbsOut(size(orbsIn,dim=1),size(orbsIn,dim=2))

  end function foldOrbsFromCell

end module dftbp_boundaryconditions
