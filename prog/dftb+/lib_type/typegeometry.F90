!--------------------------------------------------------------------------------------------------!
!  DFTB+: general package for performing fast atomistic simulations                                !
!  Copyright (C) 2017  DFTB+ developers group                                                      !
!                                                                                                  !
!  See the LICENSE file for terms of usage and distribution.                                       !
!--------------------------------------------------------------------------------------------------!

module typegeometry
  use accuracy
  implicit none
  private

  public :: TGeometry, normalize

  !!* Type for containing geometrical information
  type TGeometry
    integer           :: nAtom
    logical           :: tPeriodic
    logical           :: tFracCoord
    integer,  allocatable :: species(:)
    real(dp), allocatable :: coords(:,:)
    integer           :: nSpecies
    real(dp), allocatable :: origin(:)
    real(dp), allocatable :: latVecs(:,:)
    real(dp), allocatable :: recVecs2p(:,:)
    character(mc), allocatable :: speciesNames(:)
  end type TGeometry

  !!* Interface for cleaning up a geometry against non-existent atom species
  interface normalize
    module procedure Geometry_normalize
  end interface


contains

  !!* Normalises a geometry object to be safe against the abscence of any atoms
  !!* of a species specified in the input file
  !!* @param self Geometry object
  subroutine Geometry_normalize(sf)
    type(TGeometry), intent(inout) :: sf

    logical, allocatable :: inUse(:)
    integer, allocatable :: oldSpecies(:)
    character(mc), allocatable :: oldSpeciesNames(:)
    integer :: ind, iSp

    allocate(inUse(sf%nSpecies))
    do iSp = 1, sf%nSpecies
      inUse(iSp) = any(sf%species == iSp)
    end do
    if (.not. all(inUse)) then !some of the species are redundant, so re-index
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
