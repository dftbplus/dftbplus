!--------------------------------------------------------------------------------------------------!
!  DFTB+: general package for performing fast atomistic simulations                                !
!  Copyright (C) 2006 - 2020  DFTB+ developers group                                               !
!                                                                                                  !
!  See the LICENSE file for terms of usage and distribution.                                       !
!--------------------------------------------------------------------------------------------------!

!> Common interface for all dispersion modules.
module dftbp_dispiface
  use dftbp_accuracy, only : dp
  use dftbp_environment, only : TEnvironment
  use dftbp_periodic, only : TNeighbourList
  use dftbp_commontypes, only : TOrbitals
  implicit none
  private
  
  public :: TDispersionIface

  !> Interface for classes providing dispersion.
  type, abstract :: TDispersionIface
  contains

    !> update internal copy of coordinates
    procedure(updateCoordsIface), deferred :: updateCoords

    !> update internal copy of lattice vectors
    procedure(updateLatVecsIface), deferred :: updateLatVecs

    !> get real space cutoff
    procedure(getRCutoffIface), deferred :: getRCutoff

    !> get energy contributions
    procedure(getEnergiesIface), deferred :: getEnergies

    !> get force contributions
    procedure(addGradientsIface), deferred :: addGradients

    !> get stress tensor contributions
    procedure(getStressIface), deferred :: getStress

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
      import :: TDispersionIface, TEnvironment, TNeighbourList, dp

      !> data structure
      class(TDispersionIface), intent(inout) :: this

      !> Computational environment settings
      type(TEnvironment), intent(in) :: env

      !> list of neighbours to atoms
      type(TNeighbourList), intent(in) :: neigh

      !> image to central cell atom index
      integer, intent(in) :: img2CentCell(:)

      !> atomic coordinates
      real(dp), intent(in) :: coords(:,:)

      !> central cell chemical species
      integer, intent(in) :: species0(:)

      !> Status of operation
      integer, intent(out), optional :: stat
    end subroutine updateCoordsIface


    !> update internal copy of lattice vectors
    subroutine updateLatVecsIface(this, latVecs)
      import :: TDispersionIface, dp

      !> data structure
      class(TDispersionIface), intent(inout) :: this

      !> lattice vectors
      real(dp), intent(in) :: latVecs(:,:)
    end subroutine updateLatVecsIface


    !> get energy contributions
    subroutine getEnergiesIface(this, energies)
      import :: TDispersionIface, dp

      !> data structure
      class(TDispersionIface), intent(inout) :: this

      !> energy contributions for each atom
      real(dp), intent(out) :: energies(:)
    end subroutine getEnergiesIface


    !> get force contributions
    subroutine addGradientsIface(this, gradients)
      import :: TDispersionIface, dp

      !> data structure
      class(TDispersionIface), intent(inout) :: this

      !> gradient contributions for each atom
      real(dp), intent(inout) :: gradients(:,:)
    end subroutine addGradientsIface


    !> get stress tensor contributions
    subroutine getStressIface(this, stress)
      import :: TDispersionIface, dp

      !> data structure
      class(TDispersionIface), intent(inout) :: this

      !> Stress tensor contributions
      real(dp), intent(out) :: stress(:,:)
    end subroutine getStressIface


    !> Distance cut off for dispersion interactions
    function getRCutoffIface(this) result(cutoff)
      import :: TDispersionIface, dp

      !> data structure
      class(TDispersionIface), intent(inout) :: this

      !> resulting cutoff
      real(dp) :: cutoff
    end function getRCutoffIface

  end interface

contains

  !> update charges, dummy interface if not needed
  subroutine updateOnsiteCharges(this, qNetAtom, orb, referenceN0, species0, tCanUseCharges)

    !> data structure
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

    !> data structure
    class(TDispersionIface), intent(in) :: this

    !> Atomistic potential (dummy for most dispersion models)
    real(dp), intent(inout) :: vDisp(:)

  end subroutine addPotential


  !> Is the dispersion energy available for use in the main code after calling getEnergies
  function energyAvailable(this)

    !> data structure
    class(TDispersionIface), intent(in) :: this

    !> result (dummy for most dispersion models)
    logical :: energyAvailable

    energyAvailable = .true.

  end function energyAvailable

end module dftbp_dispiface
