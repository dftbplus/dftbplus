!--------------------------------------------------------------------------------------------------!
!  DFTB+: general package for performing fast atomistic simulations                                !
!  Copyright (C) 2006 - 2020  DFTB+ developers group                                               !
!                                                                                                  !
!  See the LICENSE file for terms of usage and distribution.                                       !
!--------------------------------------------------------------------------------------------------!

#:include 'common.fypp'

!> data type and associated routines for specifying atomic geometry and boundary conditions
module dftbp_typegeometry
  use dftbp_accuracy
  use dftbp_lapackroutines
  implicit none
  private

  public :: TGeometry, normalize
  public :: reduce, setLattice


  !> Type for containing geometrical information
  type TGeometry

    !> number of atoms
    integer :: nAtom

    !> is this periodic
    logical :: tPeriodic

    !> if periodic, is this in fractional units?
    logical :: tFracCoord

    !> atomic species
    integer, allocatable :: species(:)

    !> coordinates of atoms
    real(dp), allocatable :: coords(:,:)

    !> number of different atomic species
    integer :: nSpecies

    !> geometry origin (GEN file requirement)
    real(dp), allocatable :: origin(:)

    !> lattice vectors if periodic
    real(dp), allocatable :: latVecs(:,:)

    !> reciprocal lattice vectors in units of 2pi
    real(dp), allocatable :: recVecs2p(:,:)

    !> name(s) of the atomic species
    character(mc), allocatable :: speciesNames(:)
  end type TGeometry


  !> Interface for cleaning up a geometry against non-existent atom species
  interface normalize
    module procedure Geometry_normalize
  end interface

  !> Interface for reducing a geometry to a subset of its atoms
  interface reduce
   module procedure reduce_Geometry
  end interface

  !> Interface to set a lattice for a geometry
  interface setLattice
    module procedure setLattice_Geometry
  end interface

contains


  !> Normalises a geometry object to be safe against the absence of any atoms of a species specified
  !> in the input file
  subroutine Geometry_normalize(sf)

    !> Geometry object
    type(TGeometry), intent(inout) :: sf

    logical, allocatable :: inUse(:)
    integer, allocatable :: oldSpecies(:)
    character(mc), allocatable :: oldSpeciesNames(:)
    integer :: ind, iSp

    allocate(inUse(sf%nSpecies))
    do iSp = 1, sf%nSpecies
      inUse(iSp) = any(sf%species == iSp)
    end do
    !some of the species are redundant, so re-index
    if (.not. all(inUse)) then
      call move_alloc(sf%species, oldSpecies)
      call move_alloc(sf%speciesNames, oldSpeciesNames)
      sf%nSpecies = count(inUse)
      allocate(sf%species(size(oldSpecies)))
      allocate(sf%speciesNames(sf%nSpecies))
      ind = 1
      do iSp = 1, size(oldSpeciesNames)
        if (.not. inUse(iSp)) then
          continue
        end if
        sf%speciesNames(ind) = oldSpeciesNames(iSp)
        where (oldSpecies == iSp)
          sf%species = ind
        end where
        ind = ind + 1
      end do
    end if

  end subroutine Geometry_normalize

  !> Reduce the geometry to a subset.
  subroutine reduce_Geometry(self, iStart, iEnd, newOrigin, newLatVecs)

    !> Geometry object
    type(TGeometry), intent(inout) :: self

    !> Initial atom in the reduced geometry
    integer, intent(in) :: iStart

    !> Final atom in the reduced geometry
    integer, intent(in) :: iEnd

    !> Supercell origin - if not initially periodic, structure is converted
    real(dp), intent(in), optional :: newOrigin(:)

    !> Lattice vectors for the supercell - if not initially periodic, structure is converted
    real(dp), intent(in), optional :: newLatVecs(:,:)

    integer, allocatable :: tmpSpecies(:)
    real(dp), allocatable :: tmpCoords(:,:)


    self%nAtom = iEnd - iStart + 1
    allocate(tmpSpecies(self%nAtom))
    tmpSpecies = self%species(iStart:iEnd)
    deallocate(self%species)
    allocate(self%species(self%nAtom))
    self%species = tmpSpecies
    deallocate(tmpSpecies)

    allocate(tmpCoords(3, self%nAtom))
    tmpCoords = self%coords(:,iStart:iEnd)
    deallocate(self%coords)
    allocate(self%coords(3, self%nAtom))
    self%coords = tmpCoords
    deallocate(tmpCoords)
    if (present(newLatVecs).and.present(newOrigin)) then
      call setLattice(self, newOrigin, newLatVecs)
    end if

  end subroutine reduce_Geometry

  !> Set new lattice vectors for a geometry - if not initially periodic, structure is converted
  subroutine setLattice_Geometry(self, origin, latVecs)

    !> Geometry object
    type(TGeometry), intent(inout) :: self

    !> Supercell origin
    real(dp), intent(in) :: origin(:)

    !> Lattice vectors for the supercell
    real(dp), intent(in) :: latVecs(:,:)


    if (.not. self%tPeriodic) then
      allocate(self%origin(3))
      allocate(self%latVecs(3, 3))
      allocate(self%recVecs2p(3, 3))
      self%tPeriodic = .true.
      self%tFracCoord = .false.
    end if
    self%origin = origin
    self%latVecs = latVecs
    self%recVecs2p = self%latVecs
    call matinv(self%recVecs2p)
    self%recVecs2p = reshape(self%recVecs2p, (/3, 3/), order=(/2, 1/))

  end subroutine setLattice_Geometry


end module dftbp_typegeometry
