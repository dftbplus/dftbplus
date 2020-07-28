!--------------------------------------------------------------------------------------------------!
!  DFTB+: general package for performing fast atomistic simulations                                !
!  Copyright (C) 2006 - 2020  DFTB+ developers group                                               !
!                                                                                                  !
!  See the LICENSE file for terms of usage and distribution.                                       !
!--------------------------------------------------------------------------------------------------!

#:include 'common.fypp'

!> Contains routines to calculate the value of one or more molecular orbitals composed from STOs on
!> an equidistant grid.
module dftbp_molecularorbital
  use dftbp_assert
  use dftbp_accuracy
  use dftbp_typegeometry
  use dftbp_slater
  use dftbp_simplealgebra
  use dftbp_periodic, only: getCellTranslations, foldCoordToUnitCell
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

    !> STO for each orbital
    type(TSlaterOrbital), allocatable :: stos(:)

    !> Occupation for each orb.
    real(dp), allocatable :: occupations(:)
  end type TSpeciesBasis


  !> Data type containing information for molecular orbital calculator
  type TMolecularOrbital
    private

    !> Nr. of atoms
    integer :: nAtom

    !> Nr. of species
    integer :: nSpecies

    !> Species of each atom
    integer, allocatable :: species(:)

    !> Index array for STOs
    integer, allocatable :: iStos(:)

    !> All STOs sequentially
    type(TSlaterOrbital), allocatable :: stos(:)

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


  !> Returns the value of one or more molecular orbitals on a grid
  interface getValue
    module procedure MolecularOrbital_getValue_real
    module procedure MolecularOrbital_getValue_cmpl
  end interface

  public :: TSpeciesBasis
  public :: TMolecularOrbital, init, getValue

contains


  !> Initialises MolecularOrbital instance.
  subroutine MolecularOrbital_init(this, geometry, basis)

    !> Molecular Orbital
    type(TMolecularOrbital), intent(out) :: this

    !> Geometrical information.
    type(TGeometry), intent(in) :: geometry

    !> Basis for each species.
    type(TSpeciesBasis), intent(in) :: basis(:)

    integer :: nOrb
    integer :: ii, jj, ind, iSp
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
    do ii = 1, this%nSpecies
      nOrb = nOrb + (basis(ii)%nOrb)
    end do
    allocate(this%iStos(this%nSpecies+1))
    allocate(this%stos(nOrb))
    allocate(this%cutoffs(nOrb))
    allocate(this%angMoms(nOrb))
    ind = 1
    do ii = 1, this%nSpecies
      this%iStos(ii) = ind
      nOrb = basis(ii)%nOrb
      this%stos(ind:ind+nOrb-1) = basis(ii)%stos(1:nOrb)
      this%cutoffs(ind:ind+nOrb-1) = basis(ii)%cutoffs(1:nOrb)
      this%angMoms(ind:ind+nOrb-1) = basis(ii)%angMoms(1:nOrb)
      ind = ind + nOrb
    end do
    this%iStos(ii) = ind

    ! Count all orbitals (including m-dependence)
    nOrb = 0
    do ii = 1, this%nAtom
      iSp = this%species(ii)
      nOrb = nOrb + sum(2*this%angMoms(this%iStos(iSp):this%iStos(iSp+1)-1)+1)
    end do
    this%nOrb = nOrb

    ! Get cells to look for when adding STOs from periodic images
    this%tPeriodic = geometry%tPeriodic
    if (this%tPeriodic) then
      allocate(this%latVecs(3,3))
      allocate(this%recVecs2p(3,3))
      this%latVecs(:,:) = geometry%latVecs(:,:)
      call invert33(this%recVecs2p, this%latVecs)
      this%recVecs2p = reshape(this%recVecs2p, (/3, 3/), order=(/2, 1/))
      mCutoff = maxval(this%cutoffs)
      call getCellTranslations(this%cellVec, rCellVec, this%latVecs, &
          &this%recVecs2p, mCutoff)
      this%nCell = size(this%cellVec,dim=2)
    else
      allocate(this%latVecs(3,0))
      allocate(this%recVecs2p(3,0))
      allocate(this%cellVec(3, 1))
      this%cellVec(:,:) = 0.0_dp
      allocate(rCellVec(3, 1))
      rCellVec(:,:) = 0.0_dp
      this%nCell = 1
    end if

    ! Create coorinates for central cell and periodic images
    allocate(this%coords(3, this%nAtom, this%nCell))
    this%coords(:,:,1) = geometry%coords(:,:)
    if (this%tPeriodic) then
      call foldCoordToUnitCell(this%coords(:,:,1), this%latVecs, this%recVecs2p)
      do ii = 2, this%nCell
        do jj = 1, this%nAtom
          this%coords(:,jj, ii) = this%coords(:,jj,1) + rCellVec(:,ii)
        end do
      end do
    end if

    this%tInitialised = .true.

  end subroutine MolecularOrbital_init


  !> Returns molecular orbitals on a grid
  subroutine MolecularOrbital_getValue_real(this, origin, gridVecs, &
      &eigVecsReal, valueOnGrid, addDensities)

    !> MolecularOrbital instance
    type(TMolecularOrbital), intent(in) :: this

    !> Origin of the grid
    real(dp), intent(in) :: origin(:)

    !> Grid vectors
    real(dp), intent(in) :: gridVecs(:,:)

    !> Summation coefficients for the STOs
    real(dp), intent(in) :: eigVecsReal(:,:)

    !> Molecular orbitals on a grid
    real(dp), intent(out) :: valueOnGrid(:,:,:,:)

    !> Add densities instead of wave functions
    logical, intent(in), optional :: addDensities

    real(dp), save :: kPoints(3,0)
    integer, save :: kIndexes(0)
    complex(dp), save :: valueCmpl(0,0,0,0)
    complex(dp), save :: eigVecsCmpl(0,0)
    logical :: tAddDensities

    @:ASSERT(this%tInitialised)
    @:ASSERT(size(origin) == 3)
    @:ASSERT(all(shape(gridVecs) == (/ 3, 3 /)))
    @:ASSERT(size(eigVecsReal, dim=1) == this%nOrb)
    @:ASSERT(all(shape(valueOnGrid) > (/ 1, 1, 1, 0 /)))
    @:ASSERT(size(eigVecsReal, dim=2) == size(valueOnGrid, dim=4))

    if (present(addDensities)) then
      tAddDensities = addDensities
    else
      tAddDensities = .false.
    end if

    call local_getValue(origin, gridVecs, eigVecsReal, eigVecsCmpl, &
        &this%nAtom, this%nOrb, this%coords, this%species, this%cutoffs, &
        &this%iStos, this%angMoms, this%stos, this%tPeriodic, .true., &
        &this%latVecs, this%recVecs2p, kPoints, kIndexes, this%nCell, &
        &this%cellVec, tAddDensities, valueOnGrid, valueCmpl)

  end subroutine MolecularOrbital_getValue_real


  !> Returns molecular orbitals on a grid
  subroutine MolecularOrbital_getValue_cmpl(this, origin, gridVecs, &
      &eigVecsCmpl, kPoints, kIndexes, valueOnGrid)

    !> MolecularOrbital instance
    type(TMolecularOrbital), intent(in) :: this

    !> Origin of the grid
    real(dp), intent(in) :: origin(:)

    !> Grid vectors
    real(dp), intent(in) :: gridVecs(:,:)

    !> Summation coefficients for the STOs
    complex(dp), intent(in) :: eigVecsCmpl(:,:)

    !> Array of k-points
    real(dp), intent(in) :: kPoints(:,:)

    !> Index of the k-points in kPoints for every mol.orbital
    integer, intent(in) :: kIndexes(:)

    !> Molecular orbitals on grid on exit.
    complex(dp), intent(out) :: valueOnGrid(:,:,:,:)

    real(dp), save :: valueReal(0,0,0,0)
    real(dp), save :: eigVecsReal(0,0)
    logical, save :: tAddDensities = .false.

    @:ASSERT(this%tInitialised)
    @:ASSERT(size(origin) == 3)
    @:ASSERT(all(shape(gridVecs) == (/ 3, 3 /)))
    @:ASSERT(size(eigVecsCmpl, dim=1) == this%nOrb)
    @:ASSERT(all(shape(valueOnGrid) > (/ 0, 0, 0, 0 /)))
    @:ASSERT(size(eigVecsCmpl, dim=2) == size(valueOnGrid, dim=4))
    @:ASSERT(size(kPoints, dim=1) == 3)
    @:ASSERT(size(kPoints, dim=2) > 0)
    @:ASSERT(size(kIndexes) == size(eigVecsCmpl, dim=2))
    @:ASSERT(maxval(kIndexes) <= size(kPoints, dim=2))
    @:ASSERT(minval(kIndexes) > 0)

    call local_getValue(origin, gridVecs, eigVecsReal, eigVecsCmpl, &
        &this%nAtom, this%nOrb, this%coords, this%species, this%cutoffs, &
        &this%iStos, this%angMoms, this%stos, this%tPeriodic, .false., &
        &this%latVecs, this%recVecs2p, kPoints, kIndexes, this%nCell, &
        &this%cellVec, tAddDensities, valueReal, valueOnGrid)

  end subroutine MolecularOrbital_getValue_cmpl


  !> Returns the values of several molecular orbitals on grids.
  !> Caveat: The flag tPeriodic decides if the complex or the real version is read/written for the
  !> various parameters.
  subroutine local_getValue(origin, gridVecs, eigVecsReal, eigVecsCmpl, &
      &nAtom, nOrb, coords, species, cutoffs, iStos, angMoms, stos, tPeriodic, &
      &tReal, latVecs, recVecs2p, kPoints, kIndexes, nCell, cellVec, &
      &tAddDensities, valueReal, valueCmpl)

    !> Origin of the grid
    real(dp), intent(in) :: origin(:)

    !> Grid vectors
    real(dp), intent(in) :: gridVecs(:,:)

    !> Real eigenvectors, or null-array
    real(dp), intent(in) :: eigVecsReal(:,:)

    !> Complex eigenvectors, or null-array
    complex(dp), intent(in) :: eigVecsCmpl(:,:)

    !> Nr. of atoms
    integer, intent(in) :: nAtom

    !> Nr. of orbitals
    integer, intent(in) :: nOrb

    !> Coordinates of the atoms
    real(dp), intent(in) :: coords(:,:,:)

    !> Species for each atom
    integer, intent(in) :: species(:)

    !> Cutoff for each STO
    real(dp), intent(in) :: cutoffs(:)

    !> Starting position of the STOs for each species
    integer, intent(in) :: iStos(:)

    !> Angular moment for each STO
    integer, intent(in) :: angMoms(:)

    !> Array containing the STOs
    type(TSlaterOrbital), intent(in) :: stos(:)

    !> If the system is periodic
    logical, intent(in) :: tPeriodic

    !> If the system is real
    logical, intent(in) :: tReal

    !> Lattice vectors or null-array
    real(dp), intent(in) :: latVecs(:,:)

    !> Reciprocal vectors divided by 2pi (periodic) or a null-array (molecular)
    real(dp), intent(in) :: recVecs2p(:,:)

    !> Kpoints or null-array
    real(dp), intent(in) :: kPoints(:,:)

    !> Index of the k-points for each orbital in KPoints
    integer, intent(in) :: kIndexes(:)

    !> Nr. of cells to consider
    integer, intent(in) :: nCell

    !> Translation vector of the considered cells
    real(dp), intent(in) :: cellVec(:,:)

    !> If densities should be added instead of wave funcs
    logical, intent(in) :: tAddDensities

    !> Contains the real grid on exit
    real(dp), intent(out) :: valueReal(:,:,:,:)

    !> Contains the complex grid on exit
    complex(dp), intent(out) :: valueCmpl(:,:,:,:)

    real(dp) :: curCoords(3,3), xyz(3), diff(3), frac(3)
    real(dp) :: atomAllOrbVal(nOrb, nCell)
    logical :: nonZeroMask(nOrb), allZero
    integer :: nNonZero
    integer, target :: nonZeroIndContainer(nOrb)
    integer, pointer :: nonZeroIndices(:)
    real(dp), allocatable :: atomOrbValReal(:)
    complex(dp), allocatable :: atomOrbValCmpl(:)
    complex(dp) :: phases(nCell, size(kPoints, dim=2))
    real(dp) :: xx, val
    integer :: nPoints(4)
    integer :: ind, i1, i2, i3, iEig, iAtom, iOrb, iM, iSpecies, iL, iCell

    ! Array for the contribution of each orbital (and its periodic images)
    if (tReal) then
      allocate(atomOrbValReal(nOrb))
      nPoints = shape(valueReal)
    else
      allocate(atomOrbValCmpl(nOrb))
      nPoints = shape(valueCmpl)
    end if

    ! Phase factors for the periodic image cell. Note: This will be conjugated in the scalar product
    ! below. This is fine as, in contrast to what was published, DFTB+ uses implicitely exp(-ikr) as
    ! a phase factor, as the unpack routines assemble the lower triangular matrix with exp(ikr) as
    ! factor.
    phases(:,:) = exp((0.0_dp, 1.0_dp) * matmul(transpose(cellVec), kPoints))

    ! Loop over all grid points
    lpI3: do i3 = 1, nPoints(3)
      curCoords(:, 3) = real(i3 - 1, dp) * gridVecs(:,3)
      lpI2: do i2 = 1, nPoints(2)
        curCoords(:, 2) =  real(i2 - 1, dp) * gridVecs(:,2)
        lpI1: do i1 = 1, nPoints(1)
          curCoords(:, 1) = real(i1 - 1, dp) * gridVecs(:,1)
          xyz(:) = sum(curCoords, dim=2) + origin(:)
          if (tPeriodic) then
            frac(:) = matmul(xyz, recVecs2p)
            xyz(:) = matmul(latVecs, frac - real(floor(frac), dp))
          end if
          ! Get contribution from every atom in every cell for current point
          allZero = .true.
          lpCell: do iCell = 1, nCell
            ind = 1
            lpAtom: do iAtom = 1, nAtom
              iSpecies = species(iAtom)
              diff(:) = xyz(:) - coords(:,iAtom, iCell)
              xx = sqrt(sum(diff**2))
              lpOrb: do iOrb = iStos(iSpecies), iStos(iSpecies+1)-1
                iL = angMoms(iOrb)
                ! Calculate wave function only if atom is inside the cutoff
                if (xx <= cutoffs(iOrb)) then
                  allZero = .false.
                  call getValue(stos(iOrb), xx, val)
                  do iM = -iL, iL
                    atomAllOrbVal(ind, iCell) = val *RealTessY(iL, iM, diff, xx)
                    ind = ind + 1
                  end do
                else
                  atomAllOrbVal(ind:ind+2*iL, iCell) = 0.0_dp
                  ind = ind + 2 * iL + 1
                end if
              end do lpOrb
            end do lpAtom
          end do lpCell

          if (allZero) then
            if (tReal) then
              valueReal(i1, i2, i3, :) = 0.0_dp
            else
              valueCmpl(i1, i2, i3, :) = 0.0_dp
            end if
            cycle lpI1
          end if

          ! Establish mask and index of nonzero elements
          nonZeroMask = any(atomAllOrbVal /= 0.0_dp, dim=2)
          nNonZero = 0
          do iOrb = 1, nOrb
            if (nonZeroMask(iOrb)) then
              nNonZero = nNonZero + 1
              nonZeroIndContainer(nNonZero) = iOrb
            end if
          end do
          nonZeroIndices => nonZeroIndContainer(1:nNonZero)

          ! Sum the contribution from all cells and multiply by the provided coefficients (usually
          ! the eigenvector)
          if (tReal) then
            if (tAddDensities) then
              atomAllOrbVal = atomAllOrbVal**2
            end if
            atomOrbValReal(:) = sum(atomAllOrbVal, dim=2)
            do iEig = 1, nPoints(4)
              valueReal(i1, i2, i3, iEig) = dot_product( &
                  & atomOrbValReal(nonZeroIndices), &
                  & eigVecsReal(nonZeroIndices, iEig))
            end do
          else
            ind = 0
            do iEig = 1, nPoints(4)
              if (kIndexes(iEig) /= ind) then
                ind = kIndexes(iEig)
                atomOrbValCmpl(nonZeroIndices) = (0.0_dp, 0.0_dp)
                do iCell = 1, nCell
                  atomOrbValCmpl(nonZeroIndices) = &
                      & atomOrbValCmpl(nonZeroIndices) &
                      & + atomAllOrbVal(nonZeroIndices, iCell) &
                      & * phases(iCell, ind)
                end do
              end if
              valueCmpl(i1, i2, i3, iEig) = dot_product( &
                  & atomOrbValCmpl(nonZeroIndices), &
                  & eigVecsCmpl(nonZeroIndices, iEig))
            end do
          end if
        end do lpI1
      end do lpI2
    end do lpI3

  end subroutine local_getValue

end module dftbp_molecularorbital
