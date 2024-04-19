!--------------------------------------------------------------------------------------------------!
!  DFTB+: general package for performing fast atomistic simulations                                !
!  Copyright (C) 2006 - 2023  DFTB+ developers group                                               !
!                                                                                                  !
!  See the LICENSE file for terms of usage and distribution.                                       !
!--------------------------------------------------------------------------------------------------!

#:include 'common.fypp'

!> Implements interface for the repulsive (force-field like) potential
module dftbp_dftb_repulsive_repulsive
  use dftbp_common_accuracy, only : dp
  use dftbp_dftb_periodic, only : TNeighbourList
  implicit none

  private
  public :: TRepulsive


  !> Generic wrapper for force field-like (a.k.a. repulsive) contributions
  type, abstract :: TRepulsive
  contains

    procedure(getRCutOff), deferred :: getRCutOff
    procedure(updateLatVecs), deferred :: updateLatVecs
    procedure(updateCoords), deferred :: updateCoords
    procedure(getEnergy), deferred :: getEnergy
    procedure(getGradients), deferred :: getGradients
    procedure(getStress), deferred :: getStress

  end type TRepulsive


  abstract interface

    !> Returns the real space cutoff needed by the neighbour lists
    function getRCutOff(this) result(cutOff)
      import :: TRepulsive, dp
      implicit none

      !> Instance
      class(TRepulsive), intent(in) :: this

      !> Real space cutoff needed by the object
      real(dp) :: cutOff

    end function getRCutOff


    !> Transmits the updated lattice vectors to enable for pre-calculations.
    subroutine updateLatVecs(this, latVecs)
      import :: TRepulsive, dp
      implicit none

      !> Instance
      class(TRepulsive), intent(inout) :: this

      !> New lattice vectors
      real(dp), intent(in) :: latVecs(:,:)

    end subroutine updateLatVecs


    !> Transmits the updated coordinates to enable for pre-calculations.
    subroutine updateCoords(this, coords, species, img2CentCell, neighbourList)
      import :: TRepulsive, TNeighbourList, dp
      implicit none

      !> Instance
      class(TRepulsive), intent(inout) :: this

      !> New coordinates (including those in repeated cells). Shape: [3, nAllAtom]
      real(dp), intent(in) :: coords(:,:)

      !> Species of each atom. Shape: [nAllAtom]
      integer, intent(in) :: species(:)

      !> Mapping of atoms into the central cell. Shape: [nAllAtom]
      integer, intent(in) :: img2CentCell(:)

      !> Neighbour list.
      type(TNeighbourList), intent(in) :: neighbourList

    end subroutine updateCoords


    !> Returns the energy contribution of the repulsive
    subroutine getEnergy(this, coords, species, img2CentCell, neighbourList, Eatom, Etotal,&
        & iAtInCentralRegion)
      import :: TRepulsive, TNeighbourList, dp
      implicit none

      !> Instance.
      class(TRepulsive), intent(in) :: this

      !> All atomic coordinates
      real(dp), intent(in) :: coords(:,:)

      !> All atoms chemical species
      integer, intent(in) :: species(:)

      !> Image atom indices to central cell atoms
      integer, intent(in) :: img2CentCell(:)

      !> List of neighbours for each atom
      type(TNeighbourList), intent(in) :: neighbourList

      !> Energy for each atom
      real(dp), intent(out) :: Eatom(:)

      !> Total energy
      real(dp), intent(out) :: Etotal

      !> Atoms in the central cell (or device region if transport)
      integer, intent(in), optional :: iAtInCentralRegion(:)

    end subroutine getEnergy


    !> Returns the gradients of the repulsive interaction
    subroutine getGradients(this, coords, species, img2CentCell, neighbourList, grads)
      import :: TRepulsive, TNeighbourList, dp
      implicit none

      !> Instance
      class(TRepulsive), intent(in) :: this

      !> coordinates (x,y,z, all atoms including possible images)
      real(dp), intent(in) :: coords(:,:)

      !> Species of atoms in the central cell.
      integer, intent(in) :: species(:)

      !> Index of each atom in the central cell which the atom is mapped on.
      integer, intent(in) :: img2CentCell(:)

      !> Neighbour list
      type(TNeighbourList), intent(in) :: neighbourList

      !> Gradient for each atom. Shape: [3, nAtom]
      real(dp), intent(out) :: grads(:,:)

    end subroutine getGradients


    !> Returns the stress tensor contribution from the repulsive term
    subroutine getStress(this, coords, species, img2CentCell, neighbourList, cellVol, stress)
      import :: TRepulsive, TNeighbourList, dp
      implicit none

      !> Instance
      class(TRepulsive), intent(in) :: this

      !> coordinates (x,y,z, all atoms including possible images)
      real(dp), intent(in) :: coords(:,:)

      !> Species of atoms in the central cell.
      integer, intent(in) :: species(:)

      !> indexing array for periodic image atoms
      integer, intent(in) :: img2CentCell(:)

      !> Neighbour list
      type(TNeighbourList), intent(in) :: neighbourList

      !> cell volume.
      real(dp), intent(in) :: cellVol

      !> stress tensor
      real(dp), intent(out) :: stress(:,:)

    end subroutine getStress

  end interface


end module dftbp_dftb_repulsive_repulsive
