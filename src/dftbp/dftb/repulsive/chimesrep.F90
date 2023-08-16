!--------------------------------------------------------------------------------------------------!
!  DFTB+: general package for performing fast atomistic simulations                                !
!  Copyright (C) 2006 - 2025  DFTB+ developers group                                               !
!                                                                                                  !
!  See the LICENSE file for terms of usage and distribution.                                       !
!--------------------------------------------------------------------------------------------------!

#:include 'common.fypp'

!> Implements a repulsive correction using the ChIMES force field
module dftbp_dftb_repulsive_chimesrep
  use dftbp_common_environment, only : TEnvironment
  use dftbp_common_accuracy, only : dp
  use dftbp_common_constants, only : AA__Bohr, Bohr__AA, kcal_mol__Hartree, Hartree__kcal_mol
  use dftbp_dftb_periodic, only : TNeighbourList
  use dftbp_dftb_repulsive_repulsive, only : TRepulsive
  use dftbp_extlibs_chimes, only : TChimesCalc, TChimesCalc_init
  implicit none

  private
  public :: TChimesRepInp, TChimesRep, TChimesRep_init


#:if WITH_CHIMES

  !> Input for the Chimes repulsive calcualtor
  type :: TChimesRepInp

    !> Name of the Chimes input file
    character(:), allocatable :: chimesFile

  end type TChimesRepInp


  !> Contains the necessary data to query the ChiIMES force field
  type, extends(TRepulsive) :: TChimesRep
    private

    !> Chimes calculator
    type(TChimesCalc) :: chimesCalc

    !> Lattice vector buffer
    real(dp) :: latVecs(3, 3) = 0.0_dp

    !> Energy per atom buffer. Shape: [nAtom]
    real(dp), allocatable :: energyPerAtom(:)

    !> Gradient buffer. Shape: [3, nAtom]
    real(dp), allocatable :: gradients(:,:)

    !> Stress buffer
    real(dp) :: stress(3, 3) = 0.0_dp

  contains

    procedure :: getRCutOff => TChimesRep_getRCutOff
    procedure :: updateLatVecs => TChimesRep_updateLatVecs
    procedure :: updateCoords => TChimesRep_updateCoords
    procedure :: getEnergy => TChimesRep_getEnergy
    procedure :: getGradients => TChimesRep_getGradients
    procedure :: getStress => TChimesRep_getStress
    procedure :: updateSpecies => TChimesRep_updateSpecies

  end type TChimesRep

#:else

  type :: TChimesRepInp
  end type TChimesRepInp

  type :: TChimesRep
  end type TChimesRep

#:endif


contains

#:if WITH_CHIMES

  !> Initialises ChIMES repulsive
  subroutine TChimesRep_init(this, input, speciesNames, species0)

    !> Instance on exit
    type(TChimesRep), intent(out) :: this

    !> Name of the Chimes input file
    type(TChimesRepInp), intent(in) :: input

    !> Name of the species in the system. Shape: [nSpecies]
    character(*), intent(in) :: speciesNames(:)

    !> Species index of each atom in the central cell. Shape: [nAtom]
    integer, intent(in) :: species0(:)

    integer :: nAtom

    call TChimesCalc_init(this%chimesCalc, input%chimesFile, 1, 0)
    nAtom = size(species0)
    allocate(this%energyPerAtom(nAtom), source=0.0_dp)
    allocate(this%gradients(3, nAtom), source=0.0_dp)
    call this%updateSpecies(speciesNames, species0)

  end subroutine TChimesRep_init


  !> Notifies the calculator, that order of the species has been changed.
  subroutine TChimesRep_updateSpecies(this, speciesNames, species)

    !> Instance.
    class(TChimesRep), intent(inout) :: this

    !> Name of the species in the system. Shape: [nSpecies]
    character(*), intent(in) :: speciesNames(:)

    !> Species index of each atom. Shape: [nAtom]
    integer, intent(in) :: species(:)

    call this%chimesCalc%set_atom_types(speciesNames(species))

  end subroutine TChimesRep_updateSpecies


  !> Returns the real space cutoff (actually 0, as Chimes does not make use of the neighbour list)
  function TChimesRep_getRCutOff(this) result(cutOff)

    !> Instance
    class(TChimesRep), intent(in) :: this

    !> Real space cutoff needed by the object
    real(dp) :: cutOff

    cutOff = 0.0_dp

  end function TChimesRep_getRCutOff


  !> Transmits the updated lattice vectors to enable for pre-calculations.
  subroutine TChimesRep_updateLatVecs(this, latVecs)

    !> Instance
    class(TChimesRep), intent(inout) :: this

    !> New lattice vectors
    real(dp), intent(in) :: latVecs(:,:)

    this%latVecs(:,:) = latVecs

  end subroutine TChimesRep_updateLatVecs


  !> Transmits the updated coordinates to enable for pre-calculations.
  subroutine TChimesRep_updateCoords(this, env, coords, species, img2CentCell, neighbourList)

    !> Instance
    class(TChimesRep), intent(inout) :: this

    !> Environmet
    type(TEnvironment), intent(in) :: env

    !> New coordinates (including those in repeated cells). Shape: [3, nAllAtom]
    real(dp), intent(in) :: coords(:,:)

    !> Species of each atom. Shape: [nAllAtom]
    integer, intent(in) :: species(:)

    !> Mapping of atoms into the central cell. Shape: [nAllAtom]
    integer, intent(in) :: img2CentCell(:)

    !> Neighbour list.
    type(TNeighbourList), intent(in) :: neighbourList

    real(dp) :: energy
    integer :: nAtom

    energy = 0.0_dp
    this%gradients(:,:) = 0.0_dp
    this%stress(:,:) = 0.0_dp
    nAtom = size(this%energyPerAtom)
    call this%chimesCalc%calculate(coords(:, 1:nAtom) * Bohr__AA, this%latVecs * Bohr__AA,&
        & energy, this%gradients, this%stress)
    energy = energy * kcal_mol__Hartree
    ! Chimes only returns total energy, distribute it equally on each atom.
    this%energyPerAtom(:) = energy / real(size(this%energyPerAtom), dp)
    ! Chimes returns forces, not gradients
    this%gradients(:,:) = -this%gradients(:,:) *  (kcal_mol__Hartree / AA__Bohr)
    this%stress(:,:) = this%stress * (kcal_mol__Hartree / AA__Bohr**3)

  end subroutine TChimesRep_updateCoords


  !> Returns the energy of the repulsive.
  subroutine TChimesRep_getEnergy(this, coords, species, img2CentCell, neighbourList, Eatom,&
      & Etotal, iAtInCentralRegion)

    !> Instance.
    class(TChimesRep), intent(in) :: this

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

    if (present(iAtInCentralRegion)) then
      Eatom(:) = 0.0_dp
      Eatom(iAtInCentralRegion) = this%energyPerAtom(iAtInCentralRegion)
    else
      Eatom(:) = this%energyPerAtom
    end if
    Etotal = sum(Eatom)

  end subroutine TChimesRep_getEnergy


  !> Returns the gradients of the repulsive interaction
  subroutine TChimesRep_getGradients(this, coords, species, img2CentCell, neighbourList, grads)

    !> Instance
    class(TChimesRep), intent(in) :: this

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

    grads(:,:) = this%gradients

  end subroutine TChimesRep_getGradients


  !> Returns the stress tensor contribution from the repulsive term
  subroutine TChimesRep_getStress(this, coords, species, img2CentCell, neighbourList, cellVol,&
      & stress)

    !> Instance
    class(TChimesRep), intent(in) :: this

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

    stress(:,:) = this%stress

  end subroutine TChimesRep_getStress

#:else

  !> Dummy initializer in case code was compiled without ChIMES
  subroutine TChimesRep_init()
  end subroutine TChimesRep_init

#:endif


end module dftbp_dftb_repulsive_chimesrep
