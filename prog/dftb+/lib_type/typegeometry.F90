!--------------------------------------------------------------------------------------------------!
!  DFTB+: general package for performing fast atomistic simulations                                !
!  Copyright (C) 2018  DFTB+ developers group                                                      !
!                                                                                                  !
!  See the LICENSE file for terms of usage and distribution.                                       !
!--------------------------------------------------------------------------------------------------!

!> data type and associated routines for specifying atomic geometry and boundary conditions
module typegeometry
  use accuracy
  implicit none
  private

  public :: TGeometry, normalize


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

end module typegeometry
