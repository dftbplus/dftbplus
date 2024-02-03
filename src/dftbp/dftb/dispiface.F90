!--------------------------------------------------------------------------------------------------!
!  DFTB+: general package for performing fast atomistic simulations                                !
!  Copyright (C) 2006 - 2024  DFTB+ developers group                                               !
!                                                                                                  !
!  See the LICENSE file for terms of usage and distribution.                                       !
!--------------------------------------------------------------------------------------------------!

!> Common interface for all dispersion modules.
module dftbp_dftb_dispiface
  use dftbp_common_accuracy, only : dp
  use dftbp_common_environment, only : TEnvironment
  use dftbp_common_status, only : TStatus
  use dftbp_dftb_periodic, only : TNeighbourList
  use dftbp_type_commontypes, only : TOrbitals
  implicit none

  private
  public :: TDispersionIface

  !> Interface for classes providing dispersion.
  type, abstract :: TDispersionIface
  contains

    !> Update internal copy of coordinates
    procedure(updateCoordsIface), deferred :: updateCoords

    !> Update internal copy of lattice vectors
    procedure(updateLatVecsIface), deferred :: updateLatVecs

    !> Get real space cutoff
    procedure(getRCutoffIface), deferred :: getRCutoff

    !> Get energy contributions
    procedure(getEnergiesIface), deferred :: getEnergies

    !> Get force contributions
    procedure(addGradientsIface), deferred :: addGradients

    !> Get stress tensor contributions
    procedure(getStressIface), deferred :: getStress

    !> Updates with changed charges for the instance.
    procedure :: updateCharges

    !> Updates charges for dispersion models that make use of charges
    procedure :: updateOnsiteCharges

    !> Add potential contribution from  models which produce an atomic hamiltonian contribution
    procedure :: addPotential

    !> Is the dispersion energy available for use
    procedure :: energyAvailable

  end type TDispersionIface


  abstract interface

    !> Update internal stored coordinate
    subroutine updateCoordsIface(this, env, neigh, img2CentCell, coords, species0, stat)
      import :: TDispersionIface, TEnvironment, TNeighbourList, dp, TStatus

      !> Data structure
      class(TDispersionIface), intent(inout) :: this

      !> Computational environment settings
      type(TEnvironment), intent(in) :: env

      !> List of neighbours to atoms
      type(TNeighbourList), intent(in) :: neigh

      !> Image to central cell atom index
      integer, intent(in) :: img2CentCell(:)

      !> Atomic coordinates
      real(dp), intent(in) :: coords(:,:)

      !> Central cell chemical species
      integer, intent(in) :: species0(:)

      !> Status of operation
      type(TStatus), intent(out) :: stat

    end subroutine updateCoordsIface


    !> Update internal copy of lattice vectors
    subroutine updateLatVecsIface(this, latVecs)
      import :: TDispersionIface, dp

      !> Data structure
      class(TDispersionIface), intent(inout) :: this

      !> Lattice vectors
      real(dp), intent(in) :: latVecs(:,:)
    end subroutine updateLatVecsIface


    !> Get energy contributions
    subroutine getEnergiesIface(this, energies)
      import :: TDispersionIface, dp

      !> Data structure
      class(TDispersionIface), intent(inout) :: this

      !> Energy contributions for each atom
      real(dp), intent(out) :: energies(:)
    end subroutine getEnergiesIface


    !> Get force contributions
    subroutine addGradientsIface(this, env, neigh, img2CentCell, coords, species0, &
        & gradients, stat)
      import :: TDispersionIface, TEnvironment, TNeighbourList, dp

      !> Data structure
      class(TDispersionIface), intent(inout) :: this

      !> Computational environment settings
      type(TEnvironment), intent(in) :: env

      !> List of neighbours to atoms
      type(TNeighbourList), intent(in) :: neigh

      !> Image to central cell atom index
      integer, intent(in) :: img2CentCell(:)

      !> Atomic coordinates
      real(dp), intent(in) :: coords(:,:)

      !> Central cell chemical species
      integer, intent(in) :: species0(:)

      !> Gradient contributions for each atom
      real(dp), intent(inout) :: gradients(:,:)

      !> Status of operation
      integer, intent(out), optional :: stat
    end subroutine addGradientsIface


    !> Get stress tensor contributions
    subroutine getStressIface(this, stress)
      import :: TDispersionIface, dp

      !> Data structure
      class(TDispersionIface), intent(inout) :: this

      !> Stress tensor contributions
      real(dp), intent(out) :: stress(:,:)
    end subroutine getStressIface


    !> Distance cut off for dispersion interactions
    function getRCutoffIface(this) result(cutoff)
      import :: TDispersionIface, dp

      !> Data structure
      class(TDispersionIface), intent(inout) :: this

      !> Resulting cutoff
      real(dp) :: cutoff
    end function getRCutoffIface

  end interface

contains

  !> Update charges, dummy interface if not needed
  subroutine updateOnsiteCharges(this, qNetAtom, orb, referenceN0, species0, tCanUseCharges)

    !> Data structure
    class(TDispersionIface), intent(inout) :: this

    !> Net atomic populations
    real(dp), intent(in), allocatable :: qNetAtom(:)

    !> Atomic orbital information
    type(TOrbitals), intent(in) :: orb

    !> Reference neutal charge
    real(dp), intent(in) :: referenceN0(:,:)

    !> Atomic species of atoms
    integer, intent(in) :: species0(:)

    !> Are the charges from self-consistent output or otherwise suitable to use if needed
    logical, intent(in) :: tCanUseCharges

  end subroutine updateOnsiteCharges


  !> Adds the atomic potential contribution from suitable dispersion models, no effect otherwise
  subroutine addPotential(this, vDisp)

    !> Data structure
    class(TDispersionIface), intent(in) :: this

    !> Atomistic potential (dummy for most dispersion models)
    real(dp), intent(inout) :: vDisp(:)

  end subroutine addPotential


  !> Is the dispersion energy available for use in the main code after calling getEnergies
  function energyAvailable(this)

    !> Data structure
    class(TDispersionIface), intent(in) :: this

    !> Result (dummy for most dispersion models)
    logical :: energyAvailable

    energyAvailable = .true.

  end function energyAvailable


  !> Updates with changed charges for the instance.
  subroutine updateCharges(this, env, species, neigh, qq, q0, img2CentCell, orb)

    !> Data structure
    class(TDispersionIface), intent(inout) :: this

    !> Computational environment settings
    type(TEnvironment), intent(in) :: env

    !> Species, shape: [nAtom]
    integer, intent(in) :: species(:)

    !> Neighbour list.
    type(TNeighbourList), intent(in) :: neigh

    !> Orbital charges.
    real(dp), intent(in) :: qq(:,:,:)

    !> Reference orbital charges.
    real(dp), intent(in) :: q0(:,:,:)

    !> Mapping on atoms in central cell.
    integer, intent(in) :: img2CentCell(:)

    !> Orbital information
    type(TOrbitals), intent(in) :: orb

  end subroutine updateCharges


end module dftbp_dftb_dispiface
