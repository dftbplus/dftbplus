!--------------------------------------------------------------------------------------------------!
!  DFTB+: general package for performing fast atomistic simulations                                !
!  Copyright (C) 2006 - 2020  DFTB+ developers group                                               !
!                                                                                                  !
!  See the LICENSE file for terms of usage and distribution.                                       !
!--------------------------------------------------------------------------------------------------!

#:include 'common.fypp'

!> Interface definition for solvation models
module dftbp_solvation
  use dftbp_accuracy, only : dp
  use dftbp_commontypes, only : TOrbitals
  use dftbp_environment, only : TEnvironment
  use dftbp_periodic, only : TNeighbourList
  implicit none
  private

  public :: TSolvation

  !> Generic wrapper of for a solvation model
  type, abstract :: TSolvation
  contains

    !> update internal copy of coordinates
    procedure(updateCoords), deferred :: updateCoords

    !> update internal copy of lattice vectors
    procedure(updateLatVecs), deferred :: updateLatVecs

    !> get real space cutoff
    procedure(getRCutoff), deferred :: getRCutoff

    !> get energy contributions
    procedure(getEnergies), deferred :: getEnergies

    !> get force contributions
    procedure(addGradients), deferred :: addGradients

    !> get stress tensor contributions
    procedure(getStress), deferred :: getStress

    !> Updates with changed charges for the instance.
    procedure(updateCharges), deferred :: updateCharges

    !> Returns shifts per atom
    procedure(getShifts), deferred :: getShifts
  end type TSolvation

  abstract interface
    !> Update internal stored coordinates
    subroutine updateCoords(this, env, neighList, img2CentCell, coords, species0)
      import :: TSolvation, TEnvironment, TNeighbourList, dp

      !> Data structure
      class(TSolvation), intent(inout) :: this

      !> Computational environment settings
      type(TEnvironment), intent(in) :: env

      !> List of neighbours to atoms
      type(TNeighbourList), intent(in) :: neighList

      !> Image to central cell atom index
      integer, intent(in) :: img2CentCell(:)

      !> Atomic coordinates
      real(dp), intent(in) :: coords(:,:)

      !> Central cell chemical species
      integer, intent(in) :: species0(:)
    end subroutine updateCoords


    !> Update internal copy of lattice vectors
    subroutine updateLatVecs(this, latVecs)
      import :: TSolvation, dp

      !> Data structure
      class(TSolvation), intent(inout) :: this

      !> Lattice vectors
      real(dp), intent(in) :: latVecs(:,:)
    end subroutine updateLatVecs


    !> Get energy contributions
    subroutine getEnergies(this, energies)
      import :: TSolvation, dp

      !> Data structure
      class(TSolvation), intent(inout) :: this

      !> Energy contributions for each atom
      real(dp), intent(out) :: energies(:)
    end subroutine getEnergies


    !> Get force contributions
    subroutine addGradients(this, env, neighList, species, coords, img2CentCell, gradients)
      import :: TSolvation, TEnvironment, TNeighbourList, dp

      !> Data structure
      class(TSolvation), intent(inout) :: this

      !> Computational environment settings
      type(TEnvironment), intent(in) :: env

      !> Neighbour list.
      type(TNeighbourList), intent(in) :: neighList

      !> Specie for each atom.
      integer, intent(in) :: species(:)

      !> Coordinate of each atom.
      real(dp), intent(in) :: coords(:,:)

      !> Mapping of atoms to cetnral cell.
      integer, intent(in) :: img2CentCell(:)

      !> Gradient contributions for each atom
      real(dp), intent(inout) :: gradients(:,:)
    end subroutine addGradients


    !> Get stress tensor contributions
    subroutine getStress(this, stress)
      import :: TSolvation, dp

      !> Data structure
      class(TSolvation), intent(inout) :: this

      !> Stress tensor contributions
      real(dp), intent(out) :: stress(:,:)
    end subroutine getStress


    !> Distance cut off for dispersion interactions
    function getRCutoff(this) result(cutoff)
      import :: TSolvation, dp

      !> Data structure
      class(TSolvation), intent(inout) :: this

      !> Resulting cutoff
      real(dp) :: cutoff
    end function getRCutoff


    !> Updates with changed charges for the instance.
    subroutine updateCharges(this, env, species, neighList, qq, q0, img2CentCell, orb)
      import :: TSolvation, TEnvironment, TNeighbourList, TOrbitals, dp

      !> Data structure
      class(TSolvation), intent(inout) :: this

      !> Computational environment settings
      type(TEnvironment), intent(in) :: env

      !> Species, shape: [nAtom]
      integer, intent(in) :: species(:)

      !> Neighbour list.
      type(TNeighbourList), intent(in) :: neighList

      !> Orbital charges.
      real(dp), intent(in) :: qq(:,:,:)

      !> Reference orbital charges.
      real(dp), intent(in) :: q0(:,:,:)

      !> Mapping on atoms in central cell.
      integer, intent(in) :: img2CentCell(:)

      !> Orbital information
      type(TOrbitals), intent(in) :: orb
    end subroutine updateCharges


    !> Returns shifts per atom
    subroutine getShifts(this, shiftPerAtom, shiftPerShell)
      import :: TSolvation, dp

      !> Data structure
      class(TSolvation), intent(inout) :: this

      !> Shift per atom
      real(dp), intent(out) :: shiftPerAtom(:)

      !> Shift per shell
      real(dp), intent(out) :: shiftPerShell(:,:)
    end subroutine getShifts
  end interface

end module dftbp_solvation
