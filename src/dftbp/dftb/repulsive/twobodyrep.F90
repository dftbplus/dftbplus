!--------------------------------------------------------------------------------------------------!
!  DFTB+: general package for performing fast atomistic simulations                                !
!  Copyright (C) 2006 - 2025  DFTB+ developers group                                               !
!                                                                                                  !
!  See the LICENSE file for terms of usage and distribution.                                       !
!--------------------------------------------------------------------------------------------------!

#:include 'common.fypp'

!> Implements a repulsive potential between two atoms represented by a polynomial of 9th degree
module dftbp_dftb_repulsive_twobodyrep
  use dftbp_common_accuracy, only : dp
  use dftbp_common_constants, only : pi
  use dftbp_dftb_boundarycond, only : zAxis
  use dftbp_dftb_periodic, only : TNeighbourList, getNrOfNeighboursForAll
  use dftbp_dftb_repulsive_pairrepulsive, only : TPairRepulsive, TPairRepulsiveItem
  use dftbp_dftb_repulsive_repulsive, only :TRepulsive
  use dftbp_math_quaternions, only : rotate3
  implicit none

  private
  public :: TTwoBodyRepInp
  public :: TTwoBodyRep, TTwoBodyRep_init


  !> Input for two body repulsives
  type :: TTwoBodyRepInp

    !> Array of pairwise repulsives
    type(TPairRepulsiveItem), allocatable :: pairRepulsives(:,:)

    !> Nr. of atoms in the system.
    integer :: nAtom = 0

    !> Whether the system has helical boundary conditions
    logical :: isHelical = .false.

  end type TTwoBodyRepInp


  !> Repulsive composed of two-body pair terms
  type, extends(TRepulsive) :: TTwoBodyRep
    private
    integer :: nAtom = 0
    integer :: nSpecies = 0
    logical :: isHelical = .false.
    integer, allocatable :: nNeighbours(:)
    type(TPairRepulsiveItem), allocatable :: pairRepulsives(:,:)
  contains
    procedure :: getRCutoff
    procedure :: updateLatVecs
    procedure :: updateCoords
    procedure :: getEnergy
    procedure :: getGradients
    procedure :: getStress
  end type TTwoBodyRep


contains


  !> Initializer.
  subroutine TTwoBodyRep_init(this, input)

    !> Instance
    type(TTwoBodyRep), intent(out) :: this

    !> Input data
    type(TTwoBodyRepInp), intent(inout) :: input

    @:ASSERT(allocated(input%pairRepulsives))
    @:ASSERT(size(input%pairRepulsives, dim=1) == size(input%pairRepulsives, dim=2))

    call move_alloc(input%pairRepulsives, this%pairRepulsives)
    this%nSpecies = size(this%pairRepulsives, dim=1)
    this%nAtom = input%nAtom
    this%isHelical = input%isHelical
    allocate(this%nNeighbours(this%nAtom))

  end subroutine TTwoBodyRep_init


  !> Returns the real space cutoff needed by the neighbour lists
  function getRCutOff(this) result(cutOff)

    !> Instance
    class(TTwoBodyRep), intent(in) :: this

    !> Real space cutoff needed by the object
    real(dp) :: cutOff

    integer :: iSp1, iSp2

    cutoff = 0.0_dp
    do iSp1 = 1, this%nSpecies
      do iSp2 = 1, this%nSpecies
        cutoff = max(cutoff, this%pairRepulsives(iSp2, iSp1)%item%getCutoff())
      end do
    end do

  end function getRCutOff


  !> Transmits the updated lattice vectors to enable for pre-calculations.
  subroutine updateLatVecs(this, latVecs)

    !> Instance
    class(TTwoBodyRep), intent(inout) :: this

    !> New lattice vectors
    real(dp), intent(in) :: latVecs(:,:)

    continue

  end subroutine updateLatVecs


  !> Transmits the updated coordinates to enable for pre-calculations.
  subroutine updateCoords(this, coords, species, img2CentCell, neighbourList)

    !> Instance
    class(TTwoBodyRep), intent(inout) :: this

    !> New coordinates (including those in repeated cells). Shape: [3, nAllAtom]
    real(dp), intent(in) :: coords(:,:)

    !> Species of each atom. Shape: [nAllAtom]
    integer, intent(in) :: species(:)

    !> Mapping of atoms into the central cell. Shape: [nAllAtom]
    integer, intent(in) :: img2CentCell(:)

    !> Neighbour list.
    type(TNeighbourList), intent(in) :: neighbourList

    call getNrOfNeighboursForAll(this%nNeighbours, neighbourList, this%getRCutoff())

  end subroutine updateCoords


  !> Returns the energy contribution of the repulsive
  subroutine getEnergy(this, coords, species, img2CentCell, neighbourList, Eatom, Etotal,&
      & iAtInCentralRegion)

    !> Instance.
    class(TTwoBodyRep), intent(in) :: this

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

    call getTwoBodyEnergy_(this%pairRepulsives, coords, this%nNeighbours, neighbourList%iNeighbour,&
        & species, img2CentCell, Eatom)
    if (present(iAtInCentralRegion)) then
      Etotal = sum(Eatom(iAtInCentralRegion))
    else
      Etotal = sum(Eatom)
    end if

  end subroutine getEnergy


  !> Returns the gradients of the repulsive interaction
  subroutine getGradients(this, coords, species, img2CentCell, neighbourList, grads)

    !> Instance
    class(TTwoBodyRep), intent(in) :: this

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

    call getTwoBodyGradients_(this%pairRepulsives, coords, this%nNeighbours,&
        & neighbourList%iNeighbour, species, img2CentCell, this%isHelical, grads)

  end subroutine getGradients


  !> Returns the stress tensor contribution from the repulsive term
  subroutine getStress(this, coords, species, img2CentCell, neighbourList, cellVol, stress)

    !> Instance
    class(TTwoBodyRep), intent(in) :: this

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

    call getTwoBodyStress_(this%pairRepulsives, coords, this%nNeighbours, neighbourList%iNeighbour,&
        & species, img2CentCell, cellVol, stress)

  end subroutine getStress


  !> Subroutine for repulsive energy contributions for each atom
  subroutine getTwoBodyEnergy_(pairRepulsives, coords, nNeighbours, iNeighbours, species,&
        & img2CentCell, energy)

    type(TPairRepulsiveItem), intent(in) :: pairRepulsives(:,:)
    real(dp), intent(in) :: coords(:,:)
    integer, intent(in) :: nNeighbours(:)
    integer, intent(in) :: iNeighbours(0:,:)
    integer, intent(in) :: species(:)
    integer, intent(in) :: img2CentCell(:)
    real(dp), intent(out) :: energy(:)

    integer :: iAt1, iNeigh, iAt2, iAt2f
    real(dp) :: vect(3), dist, intermed

    @:ASSERT(size(energy) == size(nNeighbours))

    energy(:) = 0.0_dp
    do iAt1 = 1, size(nNeighbours)
      do iNeigh = 1, nNeighbours(iAt1)
        iAt2 = iNeighbours(iNeigh,iAt1)
        iAt2f = img2CentCell(iAt2)
        vect(:) = coords(:,iAt1) - coords(:,iAt2)
        dist = sqrt(sum(vect**2))
        call pairRepulsives(species(iAt2f), species(iAt1))%item%getValue(dist, energy=intermed)
        energy(iAt1) = energy(iAt1) + 0.5_dp * intermed
        if (iAt2f /= iAt1) then
          energy(iAt2f) = energy(iAt2f) + 0.5_dp * intermed
        end if
      end do
    end do

  end subroutine getTwoBodyEnergy_


  !> Subroutine for force contributions of the repulsives.
  subroutine getTwoBodyGradients_(pairRepulsives, coords, nNeighbours, iNeighbours, species,&
      & img2CentCell, isHelical, gradients)

    type(TPairRepulsiveItem), intent(in) :: pairRepulsives(:,:)
    real(dp), intent(in) :: coords(:,:)
    integer, intent(in) :: nNeighbours(:)
    integer, intent(in) :: iNeighbours(0:,:)
    integer, intent(in) :: species(:)
    integer, intent(in) :: img2CentCell(:)
    logical, intent(in) :: isHelical
    real(dp), intent(out) :: gradients(:,:)

    integer :: iAt1, iNeigh, iAt2, iAt2f
    real(dp) :: vect(3), intermed(3), dist, deriv, theta
    logical :: tHelix

    tHelix = isHelical
    gradients(:,:) = 0.0_dp

    if (tHelix) then
      do iAt1 = 1, size(nNeighbours)
        lpNeighH: do iNeigh = 1, nNeighbours(iAt1)
          iAt2 = iNeighbours(iNeigh,iAt1)
          iAt2f = img2CentCell(iAt2)
          vect(:) = coords(:,iAt1) - coords(:,iAt2)
          dist = sqrt(sum(vect**2))
          call pairRepulsives(species(iAt2f), species(iAt1))%item%getValue(dist, dEnergy=deriv)
          intermed(:) = deriv * vect / dist
          gradients(:,iAt1) = gradients(:,iAt1) + intermed(:)
          theta = - atan2(coords(2,iAt2),coords(1,iAt2)) &
              & + atan2(coords(2,iAt2f),coords(1,iAt2f))
          theta = mod(theta,2.0_dp*pi)
          call rotate3(intermed, theta, zAxis)
          if (iAt2f /= iAt1) then
            gradients(:,iAt2f) = gradients(:, iAt2f) - intermed
          end if
        end do lpNeighH
      end do

    else

      do iAt1 = 1, size(nNeighbours)
        lpNeigh: do iNeigh = 1, nNeighbours(iAt1)
          iAt2 = iNeighbours(iNeigh,iAt1)
          iAt2f = img2CentCell(iAt2)
          if (iAt2f == iAt1) then
            cycle lpNeigh
          end if
          vect(:) = coords(:,iAt1) - coords(:,iAt2)
          dist = sqrt(sum(vect**2))
          call pairRepulsives(species(iAt2f), species(iAt1))%item%getValue(dist, dEnergy=deriv)
          intermed(:) = deriv * vect / dist
          gradients(:,iAt1) = gradients(:,iAt1) + intermed(:)
          gradients(:,iAt2f) = gradients(:,iAt2f) - intermed(:)
        end do lpNeigh
      end do

    end if

  end subroutine getTwoBodyGradients_


  !> The stress tensor contribution from the repulsive energy term
  subroutine getTwoBodyStress_(pairRepulsives, coords, nNeighbours, iNeighbours, species,&
      & img2CentCell, cellVol, stress)

    type(TPairRepulsiveItem), intent(in) :: pairRepulsives(:,:)
    real(dp), intent(in) :: coords(:,:)
    integer, intent(in) :: nNeighbours(:)
    integer, intent(in) :: iNeighbours(0:,:)
    integer, intent(in) :: species(:)
    integer, intent(in) :: img2CentCell(:)
    real(dp), intent(in) :: cellVol
    real(dp), intent(out) :: stress(:,:)

    integer :: iAt1, iNeigh, iAt2, iAt2f, ii, nAtom
    real(dp) :: vect(3), intermed(3), prefac, dist, deriv

    @:ASSERT(all(shape(stress) == [3, 3]))

    nAtom = size(nNeighbours)
    stress(:,:) = 0.0_dp

    do iAt1 = 1, nAtom
      do iNeigh = 1, nNeighbours(iAt1)
        iAt2 = iNeighbours(iNeigh,iAt1)
        iAt2f = img2CentCell(iAt2)
        vect(:) = coords(:,iAt1) - coords(:,iAt2)
        dist = sqrt(sum(vect**2))
        call pairRepulsives(species(iAt2f), species(iAt1))%item%getValue(dist, dEnergy=deriv)
        intermed(:) = deriv * vect / dist
        if (iAt1 == iAt2f) then
          prefac = 0.5_dp
        else
          prefac = 1.0_dp
        end if
        do ii = 1, 3
          stress(:, ii) = stress(:, ii) - prefac * intermed * vect(ii)
        end do
      end do
    end do

    stress = stress / cellVol

  end subroutine getTwoBodyStress_


end module dftbp_dftb_repulsive_twobodyrep