!--------------------------------------------------------------------------------------------------!
!  DFTB+: general package for performing fast atomistic simulations                                !
!  Copyright (C) 2006 - 2020  DFTB+ developers group                                               !
!                                                                                                  !
!  See the LICENSE file for terms of usage and distribution.                                       !
!--------------------------------------------------------------------------------------------------!

#:include 'common.fypp'

!> Non-polar solvent accessible surface area (SASA) contributions
module dftbp_sasa
  use dftbp_accuracy, only : dp
  use dftbp_blasroutines, only : gemv
  use dftbp_commontypes, only : TOrbitals
  use dftbp_environment, only : TEnvironment
  use dftbp_lebedev, only : getAngGrid, gridSize
  use dftbp_message, only : error
  use dftbp_periodic, only : TNeighbourList, getNrOfNeighboursForAll
  use dftbp_simplealgebra, only : determinant33
  use dftbp_solvation, only : TSolvation
  implicit none
  private

  public :: TSASACont, TSASAInput


  !> Global parameters for the solvent accessible surface area model
  type :: TSASAParameters
  end type TSASAParameters

  !> Input parameters to initialize the solvent accessible surface area model
  type, extends(TSASAParameters) :: TSASAInput

    !> Real space cutoff
    real(dp) :: rCutoff = 0.0_dp

  end type TSASAInput

  !> Data for the solvent accessible surface area model
  type, extends(TSolvation) :: TSASACont

    !> number of atoms
    integer :: nAtom = 0

    !> solvation free energy
    real(dp), allocatable :: energies(:)

    !> lattice vectors if periodic
    real(dp) :: latVecs(3, 3) = 0.0_dp

    !> Volume of the unit cell
    real(dp) :: volume = 0.0_dp

    !> stress tensor
    real(dp) :: stress(3, 3) = 0.0_dp

    !> is this periodic
    logical :: tPeriodic

    !> are the coordinates current?
    logical :: tCoordsUpdated = .false.

    !> are the charges current?
    logical :: tChargesUpdated = .false.

    !> Real space cutoff
    real(dp) :: rCutoff = 0.0_dp

    !> Model parameters
    type(TSASAParameters) :: param

    !> Solvent accessible surface area
    real(dp), allocatable :: sasa(:)

    !> Derivative of solvent accessible surface area w.r.t. coordinates
    real(dp), allocatable :: dsdr(:, :, :)

    !> Derivative of solvent accessible surface area w.r.t. strain deformations
    real(dp), allocatable :: dsdL(:, :, :)

  contains

    !> Initialize solvent accessible surface area model from input data
    procedure :: initialize

    !> update internal copy of coordinates
    procedure :: updateCoords

    !> update internal copy of lattice vectors
    procedure :: updateLatVecs

    !> get real space cutoff
    procedure :: getRCutoff

    !> get energy contributions
    procedure :: getEnergies

    !> get force contributions
    procedure :: addGradients

    !> get stress tensor contributions
    procedure :: getStress

    !> Updates with changed charges for the instance
    procedure :: updateCharges

    !> Returns shifts per atom
    procedure :: getShifts

  end type TSASACont


contains


  !> Initialize generalized Born model from input data
  subroutine initialize(self, input, nAtom, species0, speciesNames, latVecs)

    !> Initialised instance at return
    class(TSASACont), intent(out) :: self

    !> Specific input parameters for generalized Born
    type(TSASAInput), intent(in) :: input

    !> Nr. of atoms in the system
    integer, intent(in) :: nAtom

    !> Species of every atom in the unit cell
    integer, intent(in) :: species0(:)

    !> Symbols of the species
    character(len=*), intent(in) :: speciesNames(:)

    !> Lattice vectors, if the system is periodic
    real(dp), intent(in), optional :: latVecs(:,:)

    integer :: iAt1, iSp1

    self%tPeriodic = present(latVecs)
    if (self%tPeriodic) then
      call self%updateLatVecs(LatVecs)
    end if
    self%nAtom = nAtom

    allocate(self%energies(nAtom))
    allocate(self%sasa(nAtom))
    allocate(self%dsdr(3, nAtom, nAtom))
    allocate(self%dsdL(3, 3, nAtom))

    self%param = input%TSASAParameters
    self%rCutoff = input%rCutoff

    self%tCoordsUpdated = .false.
    self%tChargesUpdated = .false.

  end subroutine initialize

  !> Update internal stored coordinates
  subroutine updateCoords(self, env, neighList, img2CentCell, coords, species0)

    !> Data structure
    class(TSASACont), intent(inout) :: self

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

    integer, allocatable :: nNeigh(:)

    allocate(nNeigh(self%nAtom))
    call getNrOfNeighboursForAll(nNeigh, neighList, self%rCutoff)

    self%tCoordsUpdated = .true.
    self%tChargesUpdated = .false.

  end subroutine updateCoords


  !> Update internal copy of lattice vectors
  subroutine updateLatVecs(self, latVecs)

    !> Data structure
    class(TSASACont), intent(inout) :: self

    !> Lattice vectors
    real(dp), intent(in) :: latVecs(:,:)

    @:ASSERT(self%tPeriodic)
    @:ASSERT(all(shape(latvecs) == shape(self%latvecs)))

    self%volume = abs(determinant33(latVecs))
    self%latVecs(:,:) = latVecs

    self%tCoordsUpdated = .false.
    self%tChargesUpdated = .false.

  end subroutine updateLatVecs


  !> Get energy contributions
  subroutine getEnergies(self, energies)

    !> data structure
    class(TSASACont), intent(inout) :: self

    !> energy contributions for each atom
    real(dp), intent(out) :: energies(:)

    @:ASSERT(self%tCoordsUpdated)
    @:ASSERT(self%tChargesUpdated)
    @:ASSERT(size(energies) == self%nAtom)

    energies(:) = self%energies

  end subroutine getEnergies


  !> Get force contributions
  subroutine addGradients(self, env, neighList, species, coords, img2CentCell, gradients)

    !> Data structure
    class(TSASACont), intent(inout) :: self

    !> Computational environment settings
    type(TEnvironment), intent(in) :: env

    !> Neighbour list
    type(TNeighbourList), intent(in) :: neighList

    !> Specie for each atom
    integer, intent(in) :: species(:)

    !> Coordinate of each atom
    real(dp), intent(in) :: coords(:,:)

    !> Mapping of atoms to cetnral cell
    integer, intent(in) :: img2CentCell(:)

    !> Gradient contributions for each atom
    real(dp), intent(inout) :: gradients(:,:)

    @:ASSERT(self%tCoordsUpdated)
    @:ASSERT(self%tChargesUpdated)
    @:ASSERT(all(shape(gradients) == [3, self%nAtom]))

  end subroutine addGradients


  !> get stress tensor contributions
  subroutine getStress(self, stress)

    !> data structure
    class(TSASACont), intent(inout) :: self

    !> Stress tensor contributions
    real(dp), intent(out) :: stress(:,:)

    @:ASSERT(self%tCoordsUpdated)
    @:ASSERT(self%tChargesUpdated)
    @:ASSERT(all(shape(stress) == [3, 3]))
    @:ASSERT(self%tPeriodic)
    @:ASSERT(self%volume > 0.0_dp)

    stress(:,:) = self%stress / self%volume

  end subroutine getStress


  !> Distance cut off for solvent accessible surface area calculations
  function getRCutoff(self) result(cutoff)

    !> data structure
    class(TSASACont), intent(inout) :: self

    !> resulting cutoff
    real(dp) :: cutoff

    cutoff = self%rCutoff

  end function getRCutoff


  !> Updates with changed charges for the instance.
  subroutine updateCharges(self, env, species, neighList, qq, q0, img2CentCell, orb)

    !> Data structure
    class(TSASACont), intent(inout) :: self

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

    @:ASSERT(self%tCoordsUpdated)

    self%tChargesUpdated = .true.

  end subroutine updateCharges


  !> Returns shifts per atom
  subroutine getShifts(self, shiftPerAtom, shiftPerShell)

    !> Data structure
    class(TSASACont), intent(inout) :: self

    !> Shift per atom
    real(dp), intent(out) :: shiftPerAtom(:)

    !> Shift per shell
    real(dp), intent(out) :: shiftPerShell(:,:)

    @:ASSERT(self%tCoordsUpdated)
    @:ASSERT(self%tChargesUpdated)
    @:ASSERT(size(shiftPerAtom) == self%nAtoms)
    @:ASSERT(size(shiftPerShell, dim=2) == self%nAtoms)

    shiftPerAtom(:) = 0.0_dp
    shiftPerShell(:,:) = 0.0_dp

  end subroutine getShifts


end module dftbp_sasa
