!--------------------------------------------------------------------------------------------------!
!  DFTB+: general package for performing fast atomistic simulations                                !
!  Copyright (C) 2006 - 2020  DFTB+ developers group                                               !
!                                                                                                  !
!  See the LICENSE file for terms of usage and distribution.                                       !
!--------------------------------------------------------------------------------------------------!

#:include 'common.fypp'

!> Implements a repulsive correction using the ChIMES force field
module dftbp_chimesrep
  use dftbp_accuracy, only : dp
  use dftbp_chimes, only : TChimesCalc, TChimesCalc_init
  use dftbp_constants, only : AA__Bohr, Bohr__AA, kcal_mol__Hartree, Hartree__kcal_mol
  use dftbp_periodic, only : TNeighbourList
  use dftbp_repulsive, only : TRepulsive
  implicit none

  private
  public :: TChimesRepInput, TChimesRep, TChimesRep_init


  !> Input for the Chimes repulsive calcualtor
  type :: TChimesRepInput

    !> Name of the Chimes input file
    character(:), allocatable :: chimesFile

  end type TChimesRepInput


  !> Contains the necessary data to query the ChiIMES force field
  type, extends(TRepulsive) :: TChimesRep
    private

    !> Chimes calculator
    type(TChimesCalc) :: chimesCalc_

    !> Lattice vector buffer
    real(dp) :: latVecs_(3, 3) = 0.0_dp

    !> Energy per atom buffer. Shape: [nAtom]
    real(dp), allocatable :: energyPerAtom_(:)

    !> Gradient buffer. Shape: [3, nAtom]
    real(dp), allocatable :: gradients_(:,:)

    !> Stress buffer
    real(dp) :: stress_(3, 3) = 0.0_dp

  contains

    procedure :: getRCutOff => TChimesRep_getRCutOff
    procedure :: updateLatVecs => TChimesRep_updateLatVecs
    procedure :: updateCoords => TChimesRep_updateCoords
    procedure :: getEnergy => TChimesRep_getEnergy
    procedure :: getGradients => TChimesRep_getGradients
    procedure :: getStress => TChimesRep_getStress
    procedure :: updateSpecies => TChimesRep_updateSpecies

  end type TChimesRep


contains


  !> Initialises polynomial repulsive
  subroutine TChimesRep_init(this, input, speciesNames, species0)

    !> Instance on exit
    type(TChimesRep), intent(out) :: this

    !> Name of the Chimes input file
    type(TChimesRepInput), intent(in) :: input

    !> Name of the species in the system. Shape: [nSpecies]
    character(*), intent(in) :: speciesNames(:)

    !> Species index of each atom in the central cell. Shape: [nAtom]
    integer, intent(in) :: species0(:)

    integer :: nAtom

    call TChimesCalc_init(this%chimesCalc_, input%chimesFile, 0)
    nAtom = size(species0)
    allocate(this%energyPerAtom_(nAtom), source=0.0_dp)
    allocate(this%gradients_(3, nAtom), source=0.0_dp)
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

    call this%chimesCalc_%set_atom_types(speciesNames(species))

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

    this%latVecs_(:,:) = latVecs

  end subroutine TChimesRep_updateLatVecs


  !> Transmits the updated coordinates to enable for pre-calculations.
  subroutine TChimesRep_updateCoords(this, coords, species, img2CentCell, neighbourList)

    !> Instance
    class(TChimesRep), intent(inout) :: this

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
    this%gradients_(:,:) = 0.0_dp
    this%stress_(:,:) = 0.0_dp
    nAtom = size(this%energyPerAtom_)
    print *, "CHIMES Coords:"
    print "(3F12.6)", coords(:, 1:nAtom) * Bohr__AA
    print *, "CHIMES LatVecs:"
    print "(3F12.6)", this%latVecs_ * Bohr__AA

    call this%chimesCalc_%calculate(coords(:, 1:nAtom) * Bohr__AA, this%latVecs_ * Bohr__AA,&
        & energy, this%gradients_, this%stress_)
    print "(A, 2F12.6)", "CHIMES Energy: (kcal/mol)", energy
    print "(A)", "CHIMES Forces (kcal/mol / AA):"
    print "(3F12.6)", this%gradients_
    print "(A)", "CHIMES Stress (kcal/mol / AA^3)"
    print "(3F12.6)", this%stress_

    energy = energy * kcal_mol__Hartree
    ! Chimes only returns total energy, distribute it equally on each atom.
    this%energyPerAtom_(:) = energy / real(size(this%energyPerAtom_), dp)
    ! Chimes returns forces, not gradients
    this%gradients_(:,:) = -this%gradients_(:,:) *  (kcal_mol__Hartree / AA__Bohr)
    this%stress_(:,:) = this%stress_ * (kcal_mol__Hartree / AA__Bohr**3)

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
      Eatom(iAtInCentralRegion) = this%energyPerAtom_(iAtInCentralRegion)
    else
      Eatom(:) = this%energyPerAtom_
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

    grads(:,:) = this%gradients_

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

    stress(:,:) = this%stress_

  end subroutine TChimesRep_getStress


end module dftbp_chimesrep
