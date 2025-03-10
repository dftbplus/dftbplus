!--------------------------------------------------------------------------------------------------!
!  DFTB+: general package for performing fast atomistic simulations                                !
!  Copyright (C) 2006 - 2025  DFTB+ developers group                                               !
!                                                                                                  !
!  See the LICENSE file for terms of usage and distribution.                                       !
!--------------------------------------------------------------------------------------------------!

#:include 'common.fypp'

!> Interface definition for solvation models
module dftbp_solvation_solvation
  use dftbp_common_accuracy, only : dp
  use dftbp_common_environment, only : TEnvironment
  use dftbp_common_status, only : TStatus
  use dftbp_dftb_periodic, only : TNeighbourList
  use dftbp_type_commontypes, only : TOrbitals
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

    !> Does the solvent model modify the electrostatics
    procedure(isEFieldModified), deferred :: isEFieldModified

    !> Relative dielectric constant in solvent region
    procedure(getEpsilon_r), deferred :: getEpsilon_r

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
    subroutine addGradients(this, env, neighList, species, coords, img2CentCell, gradients,&
        & errStatus)
      import :: dp, TEnvironment, TNeighbourList, TSolvation, TStatus

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

      !> Error status
      type(TStatus), intent(out) :: errStatus

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
    subroutine updateCharges(this, env, species, neighList, qq, q0, img2CentCell, orb, errStatus)
      import :: dp, TEnvironment, TNeighbourList, TOrbitals, TSolvation, TStatus

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

      !> Error status
      type(TStatus), intent(out) :: errStatus

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

    !> Is the electrostic field modified by this solvent model?
    pure function isEFieldModified(this) result(isChanged)
      import :: TSolvation

      !> Data structure
      class(TSolvation), intent(in) :: this

      !> Has the solvent model changed the electrostatic environment
      logical :: isChanged

    end function isEFieldModified


    !> Returns solvent region relative dielectric constant
    pure function getEpsilon_r(this) result(e_r)
      import :: TSolvation, dp

      !> Data structure
      class(TSolvation), intent(in) :: this

      !> epsilon_r for solvent
      real(dp) :: e_r

    end function getEpsilon_r

  end interface

end module dftbp_solvation_solvation
