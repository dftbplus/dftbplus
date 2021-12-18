!--------------------------------------------------------------------------------------------------!
!  DFTB+: general package for performing fast atomistic simulations                                !
!  Copyright (C) 2006 - 2021  DFTB+ developers group                                               !
!                                                                                                  !
!  See the LICENSE file for terms of usage and distribution.                                       !
!--------------------------------------------------------------------------------------------------!

#:include 'common.fypp'

!> Contains routines to calculate the value of one or more molecular orbitals composed from STOs on
!> an equidistant grid.
module waveplot_molorb

  use dftbp_common_accuracy, only : dp
  use dftbp_dftb_periodic, only: getCellTranslations, foldCoordToUnitCell
  use dftbp_math_simplealgebra, only : invert33
  use dftbp_type_typegeometry, only : TGeometry

  use waveplot_grids, only : TGridData, subgridsToGlobalGrid
  use waveplot_parallel, only : getStartAndEndIndices
  use waveplot_slater, only : TRadialWavefunc

#:if WITH_OMP
  use omp_lib
#:endif

  implicit none

  private
  save


  !> Data type containing information about the basis for a species
  type TSpeciesBasis

    !> Atomic number of the species
    integer :: atomicNumber

    !> Nr. of orbitals
    integer :: nOrb

    !> Angular momentum for each orb.
    integer, allocatable :: angMoms(:)

    !> Cutoff for each orbital
    real(dp), allocatable :: cutoffs(:)

    !> Radial wavefunction components
    type(TRadialWavefunc), allocatable :: rwfs(:)

    !> Occupation for each orb.
    real(dp), allocatable :: occupations(:)

  end type TSpeciesBasis


  !> Data type containing information for molecular orbital calculator
  type TMolecularOrbital

    !> Nr. of atoms
    integer :: nAtom

    !> Nr. of species
    integer :: nSpecies

    !> Species of each atom
    integer, allocatable :: species(:)

    !> Index array for STOs
    integer, allocatable :: iStos(:)

    !> All STOs sequentially
    type(TRadialWavefunc), allocatable :: rwfs(:)

    !> Cutoff for each STO
    real(dp), allocatable :: cutoffs(:)

    !> Angular mometum for each STO
    integer, allocatable :: angMoms(:)

    !> Nr. of orbitals in the system
    integer :: nOrb

    !> If sytem is periodic
    logical :: tPeriodic

    !> Lattice vectors
    real(dp), allocatable :: latVecs(:,:)

    !> Reciprocal vectors divided by 2pi
    real(dp), allocatable :: recVecs2p(:,:)

    !> Cell shift vectors
    real(dp), allocatable :: cellVec(:,:)

    !> Nr. of cell shift vectors
    integer :: nCell

    !> Coordinates in all cells
    real(dp), allocatable :: coords(:,:,:)

    !> If it is initialised
    logical :: tInitialised = .false.

  end type TMolecularOrbital


  !> Initialises a MolecularOrbital instance
  interface init
    module procedure MolecularOrbital_init
  end interface


  public :: TSpeciesBasis, TMolecularOrbital
  public :: init

contains


  !> Initialises MolecularOrbital instance.
  subroutine MolecularOrbital_init(this, geometry, basis)

    !> Molecular Orbital
    type(TMolecularOrbital), intent(out) :: this

    !> Geometrical information
    type(TGeometry), intent(in) :: geometry

    !> Basis for each species
    type(TSpeciesBasis), intent(in) :: basis(:)

    integer :: nOrb
    integer :: iCell, iSpecies, iAtom, ind, iSp
    real(dp) :: mCutoff
    real(dp), allocatable :: rCellVec(:,:)

    @:ASSERT(.not. this%tInitialised)
    @:ASSERT(geometry%nSpecies == size(basis))

    this%nAtom = geometry%nAtom
    this%nSpecies = geometry%nSpecies

    allocate(this%species(this%nAtom))
    this%species(:) = geometry%species(:)

    ! Create sequential list of STOs
    nOrb = 0
    do iSpecies = 1, this%nSpecies
      nOrb = nOrb + (basis(iSpecies)%nOrb)
    end do

    allocate(this%iStos(this%nSpecies + 1))
    allocate(this%rwfs(nOrb))
    allocate(this%cutoffs(nOrb))
    allocate(this%angMoms(nOrb))

    ind = 1

    do iSpecies = 1, this%nSpecies
      this%iStos(iSpecies) = ind
      nOrb = basis(iSpecies)%nOrb
      this%rwfs(ind:ind + nOrb - 1) = basis(iSpecies)%rwfs(1:nOrb)
      this%cutoffs(ind:ind + nOrb - 1) = basis(iSpecies)%cutoffs(1:nOrb)
      this%angMoms(ind:ind + nOrb - 1) = basis(iSpecies)%angMoms(1:nOrb)
      ind = ind + nOrb
    end do

    this%iStos(iSpecies) = ind

    ! Count all orbitals (including m-dependence)
    nOrb = 0
    do iAtom = 1, this%nAtom
      iSp = this%species(iAtom)
      nOrb = nOrb + sum(2 * this%angMoms(this%iStos(iSp):this%iStos(iSp + 1) - 1) + 1)
    end do
    this%nOrb = nOrb

    ! Get cells to look for when adding STOs from periodic images
    this%tPeriodic = geometry%tPeriodic
    if (this%tPeriodic) then

      allocate(this%latVecs(3, 3))
      allocate(this%recVecs2p(3, 3))
      this%latVecs(:,:) = geometry%latVecs(:,:)
      call invert33(this%recVecs2p, this%latVecs)
      this%recVecs2p = reshape(this%recVecs2p, [3, 3], order=[2, 1])
      mCutoff = maxval(this%cutoffs)

      ! Get cell shift vectors according to maximum cutoff
      call getCellTranslations(this%cellVec, rCellVec, this%latVecs, this%recVecs2p, mCutoff)
      this%nCell = size(this%cellVec, dim=2)

    else

      allocate(this%latVecs(3, 0))
      allocate(this%recVecs2p(3, 0))
      allocate(this%cellVec(3, 1))
      this%cellVec(:,:) = 0.0_dp
      allocate(rCellVec(3, 1))
      rCellVec(:,:) = 0.0_dp
      this%nCell = 1

    end if

    ! Create coordinates for central cell and periodic images
    allocate(this%coords(3, this%nAtom, this%nCell))
    this%coords(:,:,1) = geometry%coords
    if (this%tPeriodic) then

      call foldCoordToUnitCell(this%coords(:,:, 1), this%latVecs, this%recVecs2p)

      do iCell = 2, this%nCell
        do iAtom = 1, this%nAtom
          this%coords(:, iAtom, iCell) = this%coords(:, iAtom, 1) + rCellVec(:, iCell)
        end do
      end do

    end if

    this%tInitialised = .true.

  end subroutine MolecularOrbital_init

end module waveplot_molorb
