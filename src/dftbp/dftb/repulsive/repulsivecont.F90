!--------------------------------------------------------------------------------------------------!
!  DFTB+: general package for performing fast atomistic simulations                                !
!  Copyright (C) 2006 - 2023  DFTB+ developers group                                               !
!                                                                                                  !
!  See the LICENSE file for terms of usage and distribution.                                       !
!--------------------------------------------------------------------------------------------------!

#:include 'common.fypp'
#:include 'pointerlist.fypp'

!> Implements interface for the repulsive (force-field like) potential
module dftbp_dftb_repulsive_repulsivecont
  use dftbp_common_accuracy, only : dp
  use dftbp_dftb_periodic, only : TNeighbourList
  use dftbp_dftb_repulsive_repulsive, only : TRepulsive
  use dftbp_dftb_repulsive_repulsivelist, only : TRepulsiveList
  implicit none

  private
  public :: TRepulsiveCont, TRepulsiveCont_init


  !> Repulsives container, wrapping a list of repulsives with the repulsive interface.
  type, extends(TRepulsive) :: TRepulsiveCont
    private
    type(TRepulsiveList), allocatable :: repulsives
  contains
    procedure :: getRCutoff
    procedure :: updateCoords
    procedure :: updateLatVecs
    procedure :: getEnergy
    procedure :: getGradients
    procedure :: getStress
  end type TRepulsiveCont


contains

  !> Initializes a repulsive container.
  subroutine TRepulsiveCont_init(this, repulsives)

    !> Instance on exit
    type(TRepulsiveCont), intent(out) :: this

    !> List of repulsive to consider, will be transfered into the container.
    type(TRepulsiveList), allocatable, intent(inout) :: repulsives

    call move_alloc(repulsives, this%repulsives)

  end subroutine TRepulsiveCont_init


  !> Returns the real space cutoff needed by the neighbour lists
  function getRCutOff(this) result(cutOff)

    !> Instance
    class(TRepulsiveCont), intent(in) :: this

    !> Real space cutoff needed by the object
    real(dp) :: cutOff

    class(TRepulsive), pointer :: repulsive
    integer :: iRep

    cutoff = 0.0_dp
    do iRep = 1, this%repulsives%size()
      call this%repulsives%view(iRep, repulsive)
      cutoff = max(cutoff, repulsive%getRCutoff())
    end do

  end function getRCutOff


  !> Transmits the updated lattice vectors to enable for pre-calculations.
  subroutine updateLatVecs(this, latVecs)

    !> Instance
    class(TRepulsiveCont), intent(inout) :: this

    !> New lattice vectors
    real(dp), intent(in) :: latVecs(:,:)

    class(TRepulsive), pointer :: repulsive
    integer :: iRep

    do iRep = 1, this%repulsives%size()
      call this%repulsives%view(iRep, repulsive)
      call repulsive%updateLatVecs(latVecs)
    end do

  end subroutine updateLatVecs


  !> Transmits the updated coordinates to enable for pre-calculations.
  subroutine updateCoords(this, coords, species, img2CentCell, neighbourList)

    !> Instance
    class(TRepulsiveCont), intent(inout) :: this

    !> New coordinates (including those in repeated cells). Shape: [3, nAllAtom]
    real(dp), intent(in) :: coords(:,:)

    !> Species of each atom. Shape: [nAllAtom]
    integer, intent(in) :: species(:)

    !> Mapping of atoms into the central cell. Shape: [nAllAtom]
    integer, intent(in) :: img2CentCell(:)

    !> Neighbour list.
    type(TNeighbourList), intent(in) :: neighbourList

    class(TRepulsive), pointer :: repulsive
    integer :: iRep

    do iRep = 1, this%repulsives%size()
      call this%repulsives%view(iRep, repulsive)
      call repulsive%updateCoords(coords, species, img2CentCell, neighbourList)
    end do

  end subroutine updateCoords


  !> Returns the energy contribution of the repulsive
  subroutine getEnergy(this, coords, species, img2CentCell, neighbourList, Eatom, Etotal,&
      & iAtInCentralRegion)

    !> Instance.
    class(TRepulsiveCont), intent(in) :: this

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

    class(TRepulsive), pointer :: repulsive
    real(dp) :: EtotalBuffer
    real(dp), allocatable :: EatomBuffer(:)
    integer :: iRep

    Eatom(:) = 0.0_dp
    Etotal = 0.0_dp
    allocate(EatomBuffer(size(Eatom)))
    do iRep = 1, this%repulsives%size()
      call this%repulsives%view(iRep, repulsive)
      call repulsive%getEnergy(coords, species, img2CentCell, neighbourList, EatomBuffer,&
          & EtotalBuffer, iAtInCentralRegion)
      Eatom(:) = Eatom + EatomBuffer
      Etotal = Etotal + EtotalBuffer
    end do

  end subroutine getEnergy


  !> Returns the gradients of the repulsive interaction
  subroutine getGradients(this, coords, species, img2CentCell, neighbourList, grads)

    !> Instance
    class(TRepulsiveCont), intent(in) :: this

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

    class(TRepulsive), pointer :: repulsive
    integer :: iRep

    real(dp), allocatable :: gradsBuffer(:,:)

    grads(:,:) = 0.0_dp
    allocate(gradsBuffer, mold=grads)
    do iRep = 1, this%repulsives%size()
      call this%repulsives%view(iRep, repulsive)
      call repulsive%getGradients(coords, species, img2CentCell, neighbourList, gradsBuffer)
      grads(:,:) = grads + gradsBuffer
    end do

  end subroutine getGradients


  !> Returns the stress tensor contribution from the repulsive term
  subroutine getStress(this, coords, species, img2CentCell, neighbourList, cellVol, stress)

    !> Instance
    class(TRepulsiveCont), intent(in) :: this

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

    class(TRepulsive), pointer :: repulsive
    real(dp) :: stressBuffer(3, 3)
    integer :: iRep

    stress(:,:) = 0.0_dp
    do iRep = 1, this%repulsives%size()
      call this%repulsives%view(iRep, repulsive)
      call repulsive%getStress(coords, species, img2CentCell, neighbourList, cellVol, stressBuffer)
      stress(:,:) = stress + stressBuffer
    end do

  end subroutine getStress


end module dftbp_dftb_repulsive_repulsivecont
