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
  public :: TBoundaryConditions, boundaryTypes, BoundaryConditions_init, zAxis

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

  contains

    procedure :: foldOrbsToCell
    procedure :: foldOrbsFromCell
    procedure :: foldCoordToCell
    procedure :: foldCoordFromCell

  end type TBoundaryConditions

  !> z direction vector for rotation
  real(dp), parameter :: zAxis(3) = [0.0_dp,0.0_dp,1.0_dp]

contains


  !> Initialise the type of boundary condition on the geometry
  subroutine BoundaryConditions_init(this)

    type(TBoundaryConditions), intent(out) :: this

  end subroutine BoundaryConditions_init


  !> Folds extended coordinates back into central unit cell
  pure subroutine foldCoordToCell(this, x)

    class(TBoundaryConditions), intent(in) :: this

    real(dp), intent(inout) :: x(:)

    select case(this%boundaryType)
    case(boundaryTypes%cluster)
      return
    case(boundaryTypes%periodic3D, boundaryTypes%periodic1D)

    case(boundaryTypes%helical)

    case(boundaryTypes%helical2Op)

    case(boundaryTypes%rotational)

    end select

  end subroutine foldCoordToCell


  !> Unfolds coordinates from central unit cell to extended structure
  pure subroutine foldCoordFromCell(this, x, cellVec)

    !> Instance
    class(TBoundaryConditions), intent(in) :: this

    !> Coordinate in central cell
    real(dp), intent(inout) :: x(:)

    !> "vector" to unfolded cell in units of the boundary conditions
    real(dp), intent(in) :: cellVec(:)

    real(dp) :: theta

    select case(this%boundaryType)
    case(boundaryTypes%cluster)
      return
    case(boundaryTypes%periodic1D,boundaryTypes%periodic3D)
      x(:) = x + matmul(this%latvec, cellVec)
      return
    case(boundaryTypes%helical)
      x(3) = x(3) + this%T * cellVec(1)
      theta = this%theta * cellVec(1)
      call rotate3(x,theta, zAxis)
    case(boundaryTypes%helical2Op)
      x(3) = x(3) + this%T * cellVec(1) + this%T2 * cellVec(2)
      theta = this%theta * cellVec(1) + this%theta2 * cellVec(2)
      call rotate3(x,theta, zAxis)
    case(boundaryTypes%rotational)
      theta = this%theta * cellVec(1)
      call rotate3(x,theta, zAxis)
    end select

  end subroutine foldCoordFromCell


#:for NAME, TO, FROM in [('foldOrbsFromCell','iAtom2','iAtom2f'),&
  & ('foldOrbsToCell','iAtom2f','iAtom2')]
  !> maps orbitals between central cell and extended sparse structure
  subroutine ${NAME}$(this, orbBlock, coords, orb, img2CentCell, species, iAtom2)

    !> Instance
    class(TBoundaryConditions), intent(in) :: this

    !> Block of orbitals to fold out
    real(dp), intent(inout) :: orbBlock(:,:)

    !> atomic coordinates
    real(dp), intent(in) :: coords(:,:)

    !> data type for atomic orbital information
    type(TOrbitals), intent(in) :: orb

    !> Map from images of atoms to central cell atoms
    integer, intent(in) :: img2CentCell(:)

    !> Species of each atom
    integer, intent(in) :: species(:)

    !> Atom being folded back
    integer, intent(in) :: iAtom2

    integer :: iAtom2f, nOrb2, lShellVals(orb%mShell), iSp, iSh
    real(dp) :: theta, rotZ(size(orbBlock,dim=1),size(orbBlock,dim=1))

    select case(this%boundaryType)
    case(boundaryTypes%cluster, boundaryTypes%periodic3D, boundaryTypes%periodic1D)
      return
    case(boundaryTypes%helical,boundaryTypes%helical2Op,boundaryTypes%rotational)
      iAtom2f = img2CentCell(iAtom2)
      iSp = species(iAtom2f)
      iSh = orb%nShell(iSp)
      lShellVals(:iSh) = orb%angShell(:iSh,iSp)
      theta = atan2(coords(2,${TO}$),coords(1,${TO}$))&
          & - atan2(coords(2,${FROM}$),coords(1,${FROM}$))
      theta = mod(theta,2.0_dp*pi)
      call rotateZ(rotZ,lShellVals(:iSh), theta)
      orbBlock(:,:) = matmul(rotZ,orbBlock)
    end select

  end subroutine ${NAME}$
#:endfor

end module dftbp_boundarycond
