!--------------------------------------------------------------------------------------------------!
!  DFTB+: general package for performing fast atomistic simulations                                !
!  Copyright (C) 2006 - 2023  DFTB+ developers group                                               !
!                                                                                                  !
!  See the LICENSE file for terms of usage and distribution.                                       !
!--------------------------------------------------------------------------------------------------!

#:include 'common.fypp'

!> data type and associated routines for specifying atomic geometry and boundary conditions
module dftbp_type_extgeometry
  use dftbp_common_accuracy, only : dp
  use dftbp_common_environment, only : TEnvironment
  use dftbp_common_status, only : TStatus
  use dftbp_dftb_boundarycond, only : TBoundaryConditions, boundaryConditions
  use dftbp_dftb_periodic, only : TNeighbourList, TNeighbourList_init, getCellTranslations,&
      & updateNeighbourListAndSpecies
  use dftbp_math_lapackroutines, only : matinv
  implicit none

  private
  public :: TExtGeometry, TExtGeometry_init


  !> Extended geometry with atoms outside of the central unit cell up to a certain cutoff distance.
  type :: TExtGeometry

    !> Number of atoms in the central cell
    integer :: nAtom = -1

    !> Number of atoms in the extended geometry
    integer :: nAllAtom = -1

    !> Asymmetric (reduced) neighbour list, containing only (A, B) but not (B, A)
    type(TNeighbourList) :: asymmNeighList

    !> Coordinates of the atoms in the extended geometry, shape: [3, nAllAtom]
    real(dp), allocatable :: coords(:,:)

    !> Species of all atoms in the extended geometry, shape: [nAllAtom]
    integer, allocatable :: species(:)

    !> Mapping of each atom to an atom in the central cell
    integer, allocatable :: img2CentCell(:)

    !> Cell shift vector index for every interacting atom
    integer, allocatable :: iCellVecs(:)

    !> Cell shift vectors in relative coordinates
    real(dp), allocatable :: cellVecs(:,:)

    !> Cell shift vectors in absolute absolute (real space) coordinates
    real(dp), allocatable :: rCellVecs(:,:)

    ! Index of the boundary condition
    integer, private :: iBoundaryCond_

    ! Vectors representing the boundary conditions in the helical case, shape [1, 3]
    real(dp), allocatable, private :: helicalBoundVecs_(:,:)

  contains

    procedure :: updateGeometry

  end type TExtGeometry

contains


  !> Initializes the external geometry
  subroutine TExtGeometry_init(this, iBoundaryCond)

    !> Instance
    type(TExtGeometry), intent(out) :: this

    !> Boundary condition
    integer, intent(in) :: iBoundaryCond

    this%iBoundaryCond_ = iBoundaryCond

  end subroutine TExtGeometry_init


  !> Updates the geometry
  subroutine updateGeometry(this, env, coords0Fold, species0, cutoff, status, boundaryVectors)

    !> Instance
    class(TExtGeometry), intent(inout) :: this

    !> Environment
    type(TEnvironment), intent(in) :: env

    !> New atom coordinates in the central cell
    real(dp), intent(in) :: coords0Fold(:,:)

    !> Species of each atom in the central cell
    integer, intent(in) :: species0(:)

    !> Real space cutoff for extended geometry
    real(dp), intent(in) :: cutoff

    !> Error status
    type(TStatus), intent(out) :: status

    !> New vectors representing the boundary conditions
    !>
    !> Boundary conditions must be passed for non-cluster geometries at the first call of
    !> updateGeometry(), and afterwards whenever they change.
    !>
    real(dp), optional, intent(in) :: boundaryVectors(:,:)

    real(dp), allocatable :: invBoundVecs(:,:)

    ! Structure is either initialized or boundary conditions must be passed for non-cluster cases
    @:ASSERT(this%nAtom /= -1&
        & .or. this%iBoundaryCond_ == boundaryConditions%cluster&
        & .or. present(boundaryVectors))

    if (this%nAtom == -1) then
      call initialize_(this, size(species0), boundaryVectors)
    end if

    if (present(boundaryVectors)) then
      if (this%iBoundaryCond_ == boundaryConditions%helical) then
        this%helicalBoundVecs_ = boundaryVectors
        ! invBoundVecs are not used, when getCellTranslation() is invoked for helical cells
        allocate(invBoundVecs(0, 0))
      else
        invBoundVecs = pbc3dInvLatVecs_(boundaryVectors)
      end if
      call getCellTranslations(this%cellVecs, this%rCellVecs, boundaryVectors, invBoundVecs, cutoff)
    end if

    call updateNeighbourListAndSpecies(env, this%coords, this%species, this%img2CentCell,&
      & this%iCellVecs, this%asymmNeighList, this%nAllAtom, coords0Fold, species0, cutoff,&
      & this%rCellVecs, status, symmetric=.false., helicalBoundConds=this%helicalBoundVecs_)
    ! Error handling missing

  end subroutine updateGeometry


  ! Carries out the first initialization.
  subroutine initialize_(this, nAtom, boundaryVectors)
    type(TExtGeometry), intent(inout) :: this
    integer, intent(in) :: nAtom
    real(dp), optional, intent(in) :: boundaryVectors(:,:)

    integer, parameter :: nInitNeighbours = 40

    @:ASSERT(present(boundaryVectors) .or. this%iBoundaryCond_ /= boundaryConditions%helical)

    this%nAtom = nAtom
    if (this%iBoundaryCond_ == boundaryConditions%helical) then
      this%helicalBoundVecs_ = boundaryVectors
    end if

    select case (this%iBoundaryCond_)
    case (boundaryConditions%pbc3d)
      ! Make some guess for the nr. of all interacting atoms
      this%nAllAtom = int((real(this%nAtom, dp)**(1.0_dp/3.0_dp) + 3.0_dp)**3)
    case (boundaryConditions%helical)
      ! 1D system, so much lower number of initial interactions
      this%nAllAtom = this%nAtom + 3
      if (size(this%helicalBoundVecs_, dim=1) == 3) then
        this%nAllAtom = this%nAllAtom * nint(this%helicalBoundVecs_(3, 1))
      end if
    case default
      this%nAllAtom = this%nAtom
    end select

    allocate(this%coords(3, this%nAllAtom))
    allocate(this%species(this%nAllAtom))
    allocate(this%img2CentCell(this%nAllAtom))
    allocate(this%iCellVecs(this%nAllAtom))
    call TNeighbourlist_init(this%asymmNeighList, this%nAtom, nInitNeighbours)

    ! For cluster geometries cellVecs and rCellVecs remain constant, for others they will be
    ! dynamically (re)allocated at every updateGeometry() call, so no need to allocate them here.
    if (this%iBoundaryCond_ == boundaryConditions%cluster) then
      allocate(this%cellVecs(3, 1), source=0.0_dp)
      allocate(this%rCellVecs(3, 1), source=0.0_dp)
    end if

  end subroutine initialize_


  ! Returns the inverse lattice vectors for 3D periodic boundary conditions.
  function pbc3dInvLatVecs_(latVecs) result(invLatVecs)
    real(dp), intent(in) :: latVecs(:,:)
    real(dp), allocatable :: invLatVecs(:,:)

    invLatVecs = latVecs
    call matinv(invLatVecs)
    invLatVecs(:,:) = reshape(invLatVecs, [3, 3], order=[2, 1])

  end function pbc3dInvLatVecs_

end module dftbp_type_extgeometry
