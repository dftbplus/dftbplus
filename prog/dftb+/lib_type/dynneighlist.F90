!--------------------------------------------------------------------------------------------------!
!  DFTB+: general package for performing fast atomistic simulations                                !
!  Copyright (C) 2018  DFTB+ developers group                                                      !
!                                                                                                  !
!  See the LICENSE file for terms of usage and distribution.                                       !
!--------------------------------------------------------------------------------------------------!

#:include 'common.fypp'

!> Implements dynamic neighbor list with iterator.
!>
!> The dynamic neighbor list does not store the entire neighbor list, but creates it on the fly,
!> allowing for a low memory footprint for large neighbor lists (at the cost of speed).
!>
module dynneighlist
  use accuracy
  use assert
  use latpointiter
  use message
  implicit none
  private

  public :: TDynNeighList, TDynNeighList_init
  public :: TNeighIterator, TNeighIterator_init


  !> Dynamic neighbor list
  type :: TDynNeighList
    private

    !> Cutoff for neighbor generation
    real(dp) :: cutoff

    !> Nr. of atoms 
    integer :: nAtom

    !> Coordinates of atoms (folded into unit cell, if periodic)
    real(dp), allocatable :: coords0(:,:)

    !> Whether system is periodic
    logical :: tPeriodic

    !> Lattice vectors, if system is periodic
    real(dp), allocatable :: latVecs(:,:)

    !> Inverse lattice vectors, if system is periodic
    real(dp), allocatable :: invLatVecs(:,:)
    
  contains
    procedure :: updateLatVecs => TDynNeighList_updateLatVecs
    procedure :: updateCoords => TDynNeighList_updateCoords
  end type TDynNeighList


  !> Iterator over a dynamic neighbor list
  type :: TNeighIterator
    private

    !> Pointer to the original neighbor list
    type(TDynNeighList), pointer :: neighList

    !> Lattice point generator (if system is periodic)
    type(TLatPointIter), allocatable :: latPointGen

    !> Whether system is periodic
    logical :: tPeriodic

    !> Lattice vectors, if system is periodic
    real(dp), allocatable :: latVecs(:,:)

    !> Number of atoms
    integer :: nAtom

    !> Atom for which neighbors are returned by the iterator
    integer :: iAtom1

    !> Coordinates of the atom
    real(dp) :: coordsAtom1(3)

    !> Neighbor atom to be returned as next
    integer :: iAtom2

    !> Current cell shift vector (if system is periodic)
    real(dp) :: cellVec(3)

    !> Square of the cutoff radius
    real(dp) :: cutoff2

    !> Whether iterator has finished
    logical :: tFinished
    
  contains
    procedure :: getNextNeighbors => TNeighIterator_getNextNeighbors
  end type TNeighIterator


contains

  !> Initializes a dynamic neighbor list.
  subroutine TDynNeighList_init(this, cutoff, nAtom, tPeriodic)

    !> Initialized instance on exit.
    type(TDynNeighList), intent(out) :: this

    !> Cutoff up to which the neighbors should be generated
    real(dp), intent(in) :: cutoff

    !> Nr. of atoms in the system.
    integer, intent(in) :: nAtom

    !> Whether the system is periodic
    logical, intent(in) :: tPeriodic

    this%cutoff = cutoff
    this%nAtom = nAtom
    this%tPeriodic = tPeriodic
    allocate(this%coords0(3, this%nAtom))
    if (this%tPeriodic) then
      allocate(this%latVecs(3, 3))
      allocate(this%invLatVecs(3, 3))
    end if

  end subroutine TDynNeighList_init


  !> Updates the lattice vectors.
  !>
  !> This routine should be only called, if the neighbor list was initialized for a periodic system.
  !>
  subroutine TDynNeighList_updateLatVecs(this, latVecs, invLatVecs)

    !> Instance.
    class(TDynNeighList), intent(inout) :: this

    !> New lattice vectors.
    real(dp), intent(in) :: latVecs(:,:)

    !> Inverse lattice vectors.
    real(dp), intent(in) :: invLatVecs(:,:)

    @:ASSERT(this%tPeriodic)

    this%latVecs(:,:) = latVecs
    this%invLatVecs(:,:) = invLatVecs

  end subroutine TDynNeighList_updateLatVecs


  !> Updates the coordiantes.
  subroutine TDynNeighList_updateCoords(this, coords)

    !> Instance.
    class(TDynNeighList), intent(inout) :: this

    !> New coordinates.
    real(dp), intent(in) :: coords(:,:)

    this%coords0(:,:) = coords

  end subroutine TDynNeighList_updateCoords


  !> Initializes an iterator for the dynamic neighbors of a given atom.
  subroutine TNeighIterator_init(this, neighList, iAtom, includeSelf)

    !> Initialized instance on exit.
    type(TNeighIterator), intent(out) :: this

    !> Dynamic neighbor list containing the basic data
    type(TDynNeighList), pointer, intent(in) :: neighList

    !> Index of the atom for which the neighbors should be generated.
    integer, intent(in) :: iAtom

    !> Whether the atom itself should be also returned as first neighbor (default: false)
    logical, intent(in), optional :: includeSelf

    logical :: includeSelf0

    if (present(includeSelf)) then
      includeSelf0 = includeSelf
    else
      includeSelf0 = .false.
    end if

    this%neighList => neighList
    this%cutoff2 = this%neighList%cutoff**2
    this%nAtom = this%neighList%nAtom
    this%tPeriodic = this%neighList%tPeriodic

    this%cellVec(:) = 0.0_dp
    this%iAtom1 = iAtom
    this%coordsAtom1(:) = this%neighList%coords0(:,iAtom)
    if (includeSelf0) then
      this%iAtom2 = this%iAtom1
    else
      this%iAtom2 = this%iAtom1 + 1
    end if

    if (this%tPeriodic) then
      this%latVecs = this%neighList%latVecs
      allocate(this%latPointGen)
      call TLatPointIter_init(this%latPointGen, this%latVecs, this%neighList%invLatVecs,&
          & this%neighList%cutoff, posExtension=1, negExtension=1, excludeOrigin=.true.)
    end if

    this%tFinished = .false.

  end subroutine TNeighIterator_init


  !> Returns the next group of neighbors.
  subroutine TNeighIterator_getNextNeighbors(this, nNeighbors, coords, dists, img2CentCell)

    !> Instance.
    class(TNeighIterator), intent(inout) :: this

    !> Nr. of neighbors to return on entry, nr. of neighbors found on exit. When entry and exit
    !> values differ, the iterator has finished and can not return any more neighbors.
    integer, intent(inout) :: nNeighbors

    !> Coordinates of the neighbors. Shape: (3, nNeighbors)
    real(dp), intent(out), optional :: coords(:,:)

    !> Distances of the neighbors. Shape: (nNeighbors)
    real(dp), intent(out), optional :: dists(:)

    !> Correspoinding images of the neighbors in the central cell. Shape: (nNeighbors)
    integer, intent(out), optional :: img2CentCell(:)

    real(dp) :: neighCoords(3)
    real(dp) :: coordsTmp(3, nNeighbors), distsTmp(nNeighbors)
    integer :: img2CentCellTmp(nNeighbors)
    real(dp) :: dist2
    integer :: maxNeighs

    maxNeighs = nNeighbors
    nNeighbors = 0
    if (this%tFinished) then
      return
    end if

    do while (nNeighbors < maxNeighs)
      if (this%iAtom2 > this%nAtom) then
        if (this%tPeriodic) then
          call this%latPointGen%getNextPoint(this%cellVec, this%tFinished)
          this%cellVec = matmul(this%latVecs, this%cellVec)
        else
          this%tFinished = .true.
        end if
        if (this%tFinished) then
          exit
        end if
        this%iAtom2 = this%iAtom1
      end if

      neighCoords(:) = this%neighList%coords0(:, this%iAtom2) + this%cellVec
      dist2 = sum((this%coordsAtom1 - neighCoords)**2)
      if (dist2 <= this%cutoff2) then
        nNeighbors = nNeighbors + 1
        coordsTmp(:,nNeighbors) = neighCoords
        distsTmp(nNeighbors) = dist2
        img2CentCellTmp(nNeighbors) = this%iAtom2
      end if
      this%iAtom2 = this%iAtom2 + 1
    end do

    if (present(coords)) then
      coords(:,1:nNeighbors) = coordsTmp(:,1:nNeighbors)
    end if
    if (present(dists)) then
      dists(1:nNeighbors) = sqrt(distsTmp(1:nNeighbors))
    end if
    if (present(img2CentCell)) then
      img2CentCell(1:nNeighbors) = img2CentCellTmp(1:nNeighbors)
    end if

  end subroutine TNeighIterator_getNextNeighbors


end module dynneighlist
