!--------------------------------------------------------------------------------------------------!
!  DFTB+: general package for performing fast atomistic simulations                                !
!  Copyright (C) 2006 - 2020  DFTB+ developers group                                               !
!                                                                                                  !
!  See the LICENSE file for terms of usage and distribution.                                       !
!--------------------------------------------------------------------------------------------------!

#:include 'common.fypp'

!> Contains a calculator for the repulsive (force-field like) potential
module dftbp_repulsive
  use dftbp_assert
  use dftbp_accuracy, only : dp
  use dftbp_boundarycond, only : zAxis
  use dftbp_constants, only : pi
  use dftbp_periodic, only : TNeighbourList, getNrOfNeighboursForAll
  use dftbp_quaternions, only : rotate3
  use dftbp_repcont, only : TRepCont, getCutOff, getEnergy, getEnergyDeriv
  implicit none

  private
  public :: TRepulsiveInput, TRepulsive, TRepulsive_init


  !> Input for the repulsive generator
  type :: TRepulsiveInput

    !> Number of atoms in the system
    integer :: nAtom

    !> Whether the system is helical
    logical :: isHelical = .false.

    !> Two body repulsives
    type(TRepCont), allocatable :: twoBodyCont

    !! !> Chimes calculator
    !! type(TChimes), allocatable :: chimes

  end type TRepulsiveInput


  !> Provides a calculator for repulsive (force-field type) corrections
  type :: TRepulsive
    private

    !> Container for two body repulsives
    type(TRepCont), allocatable :: twoBodyCont_

    !> Nr. of neighbours for the two body repulsive
    integer, allocatable :: nNeighTwoBody_(:)

    !> Whether repulsive must be calculated for a helical system
    logical :: isHelical_

  contains

    procedure :: getCutoff => TRepulsive_getCutoff
    procedure :: updateCoords => TRepulsive_updateCoords
    procedure :: getEnergy => TRepulsive_getEnergy
    procedure :: getGradients => TRepulsive_getGradients
    procedure :: getStress => TRepulsive_getStress

  end type TRepulsive


contains

  !> Initializes a calculator for repulsive (force-field like) contributions
  subroutine TRepulsive_init(this, input)

    !> Instance on exit.
    type(TRepulsive), intent(out) :: this

    !> Input data
    type(TRepulsiveInput), intent(inout) :: input

    call move_alloc(input%twoBodyCont, this%twoBodyCont_)
    allocate(this%nNeighTwoBody_(input%nAtom))
    this%isHelical_ = input%isHelical

  end subroutine TRepulsive_init


  !> Returns the real space cutoff needed by the neighbour lists
  function TRepulsive_getCutOff(this) result(cutOff)

    !> Instance
    class(TRepulsive), intent(in) :: this

    !> Real space cutoff needed by the object
    real(dp) :: cutOff

    cutOff = getCutOff(this%twoBodyCont_)

  end function TRepulsive_getCutOff


  !> Transmits the updated coordinates to enable for pre-calculations.
  !>
  !> Note: The same geometry must be passed, if any of the get routines are called.
  !>
  subroutine TRepulsive_updateCoords(this, neighbourList)

    !> Instance
    class(TRepulsive), intent(inout) :: this

    !> Updated neighbour list
    type(TNeighbourList), intent(in) :: neighbourList

    ! count neighbours for repulsive interactions between atoms
    call getNrOfNeighboursForAll(this%nNeighTwoBody_, neighbourList, this%getCutoff())

  end subroutine TRepulsive_updateCoords


  !> Returns the energy contribution of the repulsive
  subroutine TRepulsive_getEnergy(this, coords, species, img2CentCell, neighbourList,&
      & iAtInCentralRegion, Eatom, Etotal)

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

    !> atoms in the central cell (or device region if transport)
    integer, intent(in) :: iAtInCentralRegion(:)

    !> Energy for each atom
    real(dp), intent(out) :: Eatom(:)

    !> Total energy
    real(dp), intent(out) :: Etotal

    call getTwoBodyEnergy_(Eatom, coords, this%nNeighTwoBody_, neighbourList%iNeighbour,&
        & species, this%twoBodyCont_, img2CentCell)
    Etotal = sum(Eatom(iAtInCentralRegion))

  end subroutine TRepulsive_getEnergy


  !> Returns the gradients of the repulsive interaction
  subroutine TRepulsive_getGradients(this, coords, neighbourList, species, img2CentCell, grads)

    !> Instance
    class(TRepulsive), intent(in) :: this

    !> coordinates (x,y,z, all atoms including possible images)
    real(dp), intent(in) :: coords(:,:)

    type(TNeighbourList), intent(in) :: neighbourList

    !> Species of atoms in the central cell.
    integer, intent(in) :: species(:)

    !> Index of each atom in the central cell which the atom is mapped on.
    integer, intent(in) :: img2CentCell(:)

    !> Gradient for each atom. Shape: [3, nAtom]
    real(dp), intent(out) :: grads(:,:)

    call getTwoBodyGradients_(grads, coords, this%nNeighTwoBody_, neighbourList%iNeighbour,&
        & species, this%twoBodyCont_, img2CentCell, this%isHelical_)

  end subroutine TRepulsive_getGradients


  !> Returns the  stress tensor contribution from the repulsive term
  subroutine TRepulsive_getStress(this, coords, neighbourList, species, img2CentCell, cellVol,&
      & stress)

    !> Instance
    class(TRepulsive), intent(in) :: this

    !> coordinates (x,y,z, all atoms including possible images)
    real(dp), intent(in) :: coords(:,:)

    !> Neighbour list
    type(TNeighbourList), intent(in) :: neighbourList

    !> Species of atoms in the central cell.
    integer, intent(in) :: species(:)

    !> indexing array for periodic image atoms
    integer, intent(in) :: img2CentCell(:)

    !> cell volume.
    real(dp), intent(in) :: cellVol

    !> stress tensor
    real(dp), intent(out) :: stress(:,:)

    call getTwoBodyStress_(stress, coords, this%nNeighTwoBody_, neighbourList%iNeighbour,&
        & species, img2CentCell, this%twoBodyCont_, cellVol)

  end subroutine TRepulsive_getStress


  !> Subroutine for repulsive energy contributions for each atom
  subroutine getTwoBodyEnergy_(reslt, coords, nNeighbourRep, iNeighbours, species, repCont,&
      & img2CentCell)

    real(dp), intent(out) :: reslt(:)
    real(dp), intent(in) :: coords(:,:)
    integer, intent(in) :: nNeighbourRep(:)
    integer, intent(in) :: iNeighbours(0:,:)
    integer, intent(in) :: species(:)
    type(TRepCont), intent(in) :: repCont
    integer, intent(in) :: img2CentCell(:)

    integer :: iAt1, iNeigh, iAt2, iAt2f
    real(dp) :: vect(3), dist, intermed

    @:ASSERT(size(reslt) == size(nNeighbourRep))

    reslt(:) = 0.0_dp
    do iAt1 = 1, size(nNeighbourRep)
      do iNeigh = 1, nNeighbourRep(iAt1)
        iAt2 = iNeighbours(iNeigh,iAt1)
        iAt2f = img2CentCell(iAt2)
        vect(:) = coords(:,iAt1) - coords(:,iAt2)
        dist = sqrt(sum(vect**2))
        call getEnergy(repCont, intermed, dist, species(iAt1), species(iAt2))
        reslt(iAt1) = reslt(iAt1) + 0.5_dp * intermed
        if (iAt2f /= iAt1) then
          reslt(iAt2f) = reslt(iAt2f) + 0.5_dp * intermed
        end if
      end do
    end do

  end subroutine getTwoBodyEnergy_


  !> Subroutine for force contributions of the repulsives.
  subroutine getTwoBodyGradients_(reslt, coords, nNeighbourRep, iNeighbours, species, repCont,&
      & img2CentCell, tHelical)

    real(dp), intent(out) :: reslt(:,:)
    real(dp), intent(in) :: coords(:,:)
    integer, intent(in) :: nNeighbourRep(:)
    integer, intent(in) :: iNeighbours(0:,:)
    integer, intent(in) :: species(:)
    type(TRepCont), intent(in) :: repCont
    integer, intent(in) :: img2CentCell(:)
    logical, intent(in) :: tHelical

    integer :: iAt1, iNeigh, iAt2, iAt2f
    real(dp) :: vect(3), intermed(3), theta
    logical :: tHelix

    tHelix = tHelical
    reslt(:,:) = 0.0_dp

    if (tHelix) then
      do iAt1 = 1, size(nNeighbourRep)
        lpNeighH: do iNeigh = 1, nNeighbourRep(iAt1)
          iAt2 = iNeighbours(iNeigh,iAt1)
          iAt2f = img2CentCell(iAt2)
          vect(:) = coords(:,iAt1) - coords(:,iAt2)
          call getEnergyDeriv(repCont, intermed, vect, species(iAt1), species(iAt2))
          reslt(:,iAt1) = reslt(:,iAt1) + intermed(:)
          theta = - atan2(coords(2,iAt2),coords(1,iAt2)) &
              & + atan2(coords(2,iAt2f),coords(1,iAt2f))
          theta = mod(theta,2.0_dp*pi)
          call rotate3(intermed, theta, zAxis)
          if (iAt2f /= iAt1) then
            reslt(:,iAt2f) = reslt(:,iAt2f) - intermed(:)
          end if
        end do lpNeighH
      end do

    else

      do iAt1 = 1, size(nNeighbourRep)
        lpNeigh: do iNeigh = 1, nNeighbourRep(iAt1)
          iAt2 = iNeighbours(iNeigh,iAt1)
          iAt2f = img2CentCell(iAt2)
          if (iAt2f == iAt1) then
            cycle lpNeigh
          end if
          vect(:) = coords(:,iAt1) - coords(:,iAt2)
          call getEnergyDeriv(repCont, intermed, vect, species(iAt1), species(iAt2))
          reslt(:,iAt1) = reslt(:,iAt1) + intermed(:)
          reslt(:,iAt2f) = reslt(:,iAt2f) - intermed(:)
        end do lpNeigh
      end do

    end if

  end subroutine getTwoBodyGradients_


  !> The stress tensor contribution from the repulsive energy term
  subroutine getTwoBodyStress_(stress, coords, nNeighbourRep, iNeighbours, species, img2CentCell,&
      & repCont, cellVol)

    !> stress tensor
    real(dp), intent(out) :: stress(:,:)

    !> coordinates (x,y,z, all atoms including possible images)
    real(dp), intent(in) :: coords(:,:)

    !> Number of neighbours for atoms in the central cell
    integer, intent(in) :: nNeighbourRep(:)

    !> Index of neighbours for a given atom.
    integer, intent(in) :: iNeighbours(0:,:)

    !> Species of atoms in the central cell.
    integer, intent(in) :: species(:)

    !> indexing array for periodic image atoms
    integer, intent(in) :: img2CentCell(:)

    !> Container for repulsive potentials.
    type(TRepCont), intent(in) :: repCont

    !> cell volume.
    real(dp), intent(in) :: cellVol

    integer :: iAt1, iNeigh, iAt2, iAt2f, ii, nAtom
    real(dp) :: vect(3), intermed(3), prefac

    @:ASSERT(all(shape(stress) == [3, 3]))

    nAtom = size(nNeighbourRep)
    stress(:,:) = 0.0_dp

    do iAt1 = 1, nAtom
      do iNeigh = 1, nNeighbourRep(iAt1)
        iAt2 = iNeighbours(iNeigh,iAt1)
        iAt2f = img2CentCell(iAt2)
        vect(:) = coords(:,iAt1) - coords(:,iAt2)
        call getEnergyDeriv(repCont, intermed, vect, species(iAt1), species(iAt2))
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


end module dftbp_repulsive
