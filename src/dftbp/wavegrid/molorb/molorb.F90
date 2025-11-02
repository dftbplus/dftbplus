!--------------------------------------------------------------------------------------------------!
!  DFTB+: general package for performing fast atomistic simulations                                !
!  Copyright (C) 2006 - 2025  DFTB+ developers group                                               !
!                                                                                                  !
!  See the LICENSE file for terms of usage and distribution.                                       !
!--------------------------------------------------------------------------------------------------!

#:include 'common.fypp'

!> Provides routines to calculate the value of one or more molecular orbitals on
!! an equidistant grid.
!! Contains argument preparation and public interfaces to the private evaluateParallel routine.
module dftbp_wavegrid_molorb
  use dftbp_wavegrid_molorb_parallel, only : evaluateParallel
  use dftbp_wavegrid_molorb_types, only : TCalculationContext, TPeriodicParams, TSystemParams
  use dftbp_wavegrid_basis, only : TOrbital, TOrbitalWrapper
  use dftbp_common_accuracy, only : dp
  use dftbp_common_constants, only : imag
  use dftbp_dftb_boundarycond, only : TBoundaryConds
  use dftbp_dftb_periodic, only : getCellTranslations
  use dftbp_io_message, only : error
  use dftbp_math_simplealgebra, only : invert33
  use dftbp_type_typegeometry, only : TGeometry

  implicit none

  private
  !> Data type containing information for molecular orbital calculator.
  type TMolecularOrbital
    !> System composition and coordinates
    type(TSystemParams) :: system

    !> Periodic boundary conditions
    type(TPeriodicParams) :: periodic

    !> Basis set in AoS format
    type(TOrbitalWrapper), allocatable :: orbitals(:)

    !> Boundary conditions handler for coordinate recalculation
    type(TBoundaryConds) :: boundaryCond

    logical :: isInitialised = .false.
  contains
    private
    procedure, public :: updateCoords => TMolecularOrbital_updateCoords
    procedure :: initPeriodic => TMolecularOrbital_initPeriodic
    procedure :: initSpeciesMapping => TMolecularOrbital_initSpeciesMapping
    procedure :: flattenBasis => TMolecularOrbital_flattenBasis
  end type TMolecularOrbital

  !> Basis set for one species
  type TSpeciesBasis
    !> Array of orbitals for this species
    type(TOrbitalWrapper), allocatable :: orbitals(:)
  end type TSpeciesBasis

  !> Returns the value of one or more molecular orbitals on a grid
  interface getValue
    module procedure TMolecularOrbital_getValue_real
    module procedure TMolecularOrbital_getValue_cmpl
  end interface

  !> Returns the total charge density on a grid
  interface getTotalChrg
    module procedure TMolecularOrbital_getTotalChrg_real
    module procedure TMolecularOrbital_getTotalChrg_cmpl
  end interface

  !> Returns the atomic densities on a grid
  interface getAtomicDensities
    module procedure TMolecularOrbital_getAtomicDensities_real
  end interface

  public :: TSpeciesBasis, TMolecularOrbital, TMolecularOrbital_init
  public :: getValue, getTotalChrg, getAtomicDensities

contains


  !> Initialises MolecularOrbital instance.
  !> This prepares data structures for the molorb calculation.
  !> The coordinates may be updated later by calling updateCoords on the molorb data object.
  subroutine TMolecularOrbital_init(this, geometry, boundaryCond, basisInput, origin, gridVecs)
    !> TMolecularOrbital data object to initialise
    class(TMolecularOrbital), intent(out) :: this

    !> System geometry
    type(TGeometry), intent(in) :: geometry

    !> Boundary conditions
    type(TBoundaryConds), intent(in) :: boundaryCond

    !> Basis set for all species
    type(TSpeciesBasis), intent(in) :: basisInput(:)

    !> Output grid origin
    real(dp), intent(in) :: origin(3)

    !> Output grid vectors
    real(dp), intent(in) :: gridVecs(3,3)

    @:ASSERT(.not. this%isInitialised)
    @:ASSERT(size(origin) == 3)
    @:ASSERT(all(shape(gridVecs) == [3, 3]))
    this%boundaryCond = boundaryCond
    this%system%origin = origin
    this%system%gridVecs = gridVecs
    
    call this%initSpeciesMapping(geometry, basisInput)
    call this%flattenBasis(basisInput)
    call this%initPeriodic(geometry, 10.0_dp)!sqrt(maxval(this%orbitals%o%cutoffSq)))
    call this%updateCoords(geometry)

    this%isInitialised = .true.
  end subroutine TMolecularOrbital_init


  !> Initialises non-geometric system parameters 
  !! Counts total number of orbitals and Orbitals and creates index maps.
  subroutine TMolecularOrbital_initSpeciesMapping(this, geometry, basis)
    !> TMolecularOrbital data object to initialise
    class(TMolecularOrbital), intent(inout) :: this

    !> System geometry
    type(TGeometry), intent(in) :: geometry

    !> Basis set for all species
    type(TSpeciesBasis), intent(in) :: basis(:)

    integer :: nOrbTotal, iSpec, iAtom, ind, iOrb, angMom

    @:ASSERT(geometry%nSpecies == size(basis))

    this%system%nAtom = geometry%nAtom
    this%system%nSpecies = geometry%nSpecies
    allocate(this%system%species(this%system%nAtom))
    this%system%species(:) = geometry%species


    ! Create Orbital index map: iStos(i) points to the first Orbital of species i
    allocate(this%system%iStos(this%system%nSpecies + 1))
    ind = 1
    do iSpec = 1, this%system%nSpecies
      this%system%iStos(iSpec) = ind
      ind = ind + size(basis(iSpec)%orbitals)
    end do
    this%system%iStos(this%system%nSpecies + 1) = ind

    ! Count total number of orbitals (including m-dependence)
    nOrbTotal = 0
    do iAtom = 1, this%system%nAtom
      iSpec = this%system%species(iAtom)
      do iOrb = 1, size(basis(iSpec)%orbitals)
        angMom = basis(iSpec)%orbitals(iOrb)%o%angMom
        nOrbTotal = nOrbTotal + 1 + 2 * angMom
      end do
    end do
    this%system%nOrb = nOrbTotal

    this%system%speciesInitialised = .true.
  end subroutine TMolecularOrbital_initSpeciesMapping


  !> Flattens the basis set into an array of Orbitals. 
  subroutine TMolecularOrbital_flattenBasis(this, basisInput)
    !> TMolecularOrbital data object to initialise
    class(TMolecularOrbital), intent(inout) :: this

    !> Basis set for all species
    type(TSpeciesBasis), intent(in) :: basisInput(:)

    integer :: iSpec, iOrb, ind, nOrbitals
    
    ! Count total number of Orbitals
    nOrbitals = 0
    do iSpec = 1, this%system%nSpecies
      nOrbitals = nOrbitals + size(basisInput(iSpec)%orbitals)
    end do

    ! Allocate flat array
    allocate(this%orbitals(nOrbitals))

    ! Copy all Orbitals into array
    ind = 1
    do iSpec = 1, this%system%nSpecies
      do iOrb = 1, size(basisInput(iSpec)%orbitals)
        ! this%orbitals(ind)%o = basisInput(iSpec)%orbitals(iOrb)%o
        allocate(this%orbitals(ind)%o, source=basisInput(iSpec)%orbitals(iOrb)%o)
        ind = ind + 1
      end do
    end do
  end subroutine TMolecularOrbital_flattenBasis


  !> Initializes periodic parameters
  subroutine TMolecularOrbital_initPeriodic(this, geometry, maxCutoff)
    !> TMolecularOrbital data object to initialise
    class(TMolecularOrbital), intent(inout) :: this

    !> System geometry
    type(TGeometry), intent(in) :: geometry

    !> Maximum cutoff radius for images to include
    real(dp), intent(in) :: maxCutoff
  
    this%periodic%isPeriodic = geometry%tPeriodic
    if (this%periodic%isPeriodic) then
      allocate(this%periodic%latVecs(3, 3))
      allocate(this%periodic%recVecs2pi(3, 3))

      this%periodic%latVecs(:,:) = geometry%latVecs
      call invert33(this%periodic%recVecs2pi, this%periodic%latVecs)
      this%periodic%recVecs2pi(:,:) = transpose(this%periodic%recVecs2pi)
      call getCellTranslations(this%periodic%fCellVec, this%periodic%rCellVec, this%periodic%latVecs, &
          & this%periodic%recVecs2pi, maxCutoff)
      this%periodic%nCell = size(this%periodic%fCellVec, dim=2)
    else
      this%periodic%nCell = 1
      allocate(this%periodic%latVecs(3, 0))
      allocate(this%periodic%recVecs2pi(3, 0))
      allocate(this%periodic%fCellVec(3, 1))
      allocate(this%periodic%rCellVec(3, 1))
      this%periodic%fCellVec(:,:) = 0.0_dp
      this%periodic%rCellVec(:,:) = 0.0_dp
    end if
    this%periodic%isInitialized = .true.
  end subroutine TMolecularOrbital_initPeriodic


  !> Initializes coordinates including periodic images.
  ! Depends on initPeriodic having run first.
  subroutine TMolecularOrbital_updateCoords(this, geometry)
    !> TMolecularOrbital data object to initialise
    class(TMolecularOrbital), intent(inout) :: this

    ! System geometry
    type(TGeometry), intent(in) :: geometry

    integer :: iCell, iAtom

    @:ASSERT(this%system%speciesInitialised)
    @:ASSERT(this%periodic%isInitialized)
    allocate(this%system%coords(3, this%system%nAtom, this%periodic%nCell))
    this%system%coords(:,:,1) = geometry%coords

    if (this%periodic%isPeriodic) then
      call this%boundaryCond%foldCoordsToCell(this%system%coords(:,:,1), this%periodic%latVecs)
      do iCell = 2, this%periodic%nCell
        do iAtom = 1, this%system%nAtom
          this%system%coords(:, iAtom, iCell) = this%system%coords(:, iAtom, 1) + this%periodic%rCellVec(:, iCell)
        end do
      end do
    end if
    this%system%coordsInitialised = .true.
  end subroutine TMolecularOrbital_updateCoords

  !> Bundles input flags / flattens optional arguments into a calculation context struct.
  function bundleFlags(isRealInput, addAtomicDensities, useGpu, occupationVec) result(ctx)
    !> Is the input eigenvector in eigVecsReal and not in eigVecsCmpl?
    logical, intent(in) :: isRealInput

    !> Calculate atomic densities (squared basis contributions) instead of wave functions?
    logical, intent(in), optional :: addAtomicDensities

    !> Enable CUDA GPU offloading?
    logical, intent(in), optional :: useGpu
    
    !> If present, calculate total charge density instead of wave functions.
    real(dp), intent(in), optional :: occupationVec(:)
    
    !> Bundled output flag struc
    type(TCalculationContext) :: ctx

    ctx%isRealInput = isRealInput

    ctx%calcAtomicDensity = .false.
    if (present(addAtomicDensities)) then
      ctx%calcAtomicDensity = addAtomicDensities
    end if

    ctx%calcTotalChrg = .false.
    if (present(occupationVec)) then
      ctx%calcTotalChrg = .true.
    end if
    
    ctx%runOnGPU = .false.
    if (present(useGpu)) then
      #:if WITH_CUDA
        ctx%runOnGPU = useGpu
      #:else
      if (useGpu) then
        call error("GPU offloading requested (useGpu=.true.), but not available in this build. (missing WITH_CUDA)")
      end if
      #:endif
    end if

    ctx%isRealOutput = ctx%isRealInput .or. ctx%calcTotalChrg

  end function bundleFlags


  !> Returns molecular orbitals on a real grid. 
  subroutine TMolecularOrbital_getValue_real(this, eigVecsReal, valueOnGrid, useGpu)
    !> MolecularOrbital data instance
    type(TMolecularOrbital), intent(in) :: this

    !> Summation coefficients for the Orbitals
    real(dp), intent(in) :: eigVecsReal(:,:)

    !> Molecular orbitals on a grid
    real(dp), intent(out) :: valueOnGrid(:,:,:,:)

    !> Enable GPU offloading?
    logical, intent(in), optional :: useGpu
  
    call TMolecularOrbital_getValue_real_generic(this, eigVecsReal, valueOnGrid, useGpu)
  end subroutine TMolecularOrbital_getValue_real


  !> Returns the total charge density on a grid.
  !> This squares each state and sums them up weighted by occupationVec.
  subroutine TMolecularOrbital_getTotalChrg_real(this, eigVecsReal, valueOnGrid, occupationVec, useGpu)
    !> MolecularOrbital instance
    type(TMolecularOrbital), intent(in) :: this

    !> Summation coefficients for the Orbitals
    real(dp), intent(in) :: eigVecsReal(:,:)

    !> Molecular orbitals on a grid
    real(dp), intent(out) :: valueOnGrid(:,:,:,:)

    !> Calculate total charge. Coefficients for each squared state.
    real(dp), intent(in) :: occupationVec(:)

    !> Enable GPU offloading?
    logical, intent(in), optional :: useGpu
  
    call TMolecularOrbital_getValue_real_generic(this, eigVecsReal, valueOnGrid, useGpu, &
        & addAtomicDensities=.false., occupationVec=occupationVec)
  end subroutine TMolecularOrbital_getTotalChrg_real


  !> Calculates the atomic densities by squaring each Orbital *before* summation.
  subroutine TMolecularOrbital_getAtomicDensities_real(this, eigVecsReal, valueOnGrid, useGpu)
    !> MolecularOrbital instance
    type(TMolecularOrbital), intent(in) :: this

    !> Summation coefficients for the Orbitals
    real(dp), intent(in) :: eigVecsReal(:,:)

    !> Molecular orbitals on a grid
    real(dp), intent(out) :: valueOnGrid(:,:,:,:)

    !> Enable GPU offloading?
    logical, intent(in), optional :: useGpu

    call TMolecularOrbital_getValue_real_generic(this, eigVecsReal, valueOnGrid, useGpu, &
        & addAtomicDensities=.true.)
  end subroutine TMolecularOrbital_getAtomicDensities_real


  !> Returns molecular orbitals on a complex grid.
  subroutine TMolecularOrbital_getValue_cmpl(this, eigVecsCmpl, kPoints, kIndexes, valueOnGrid, useGpu)
    !> MolecularOrbital instance
    type(TMolecularOrbital), intent(in) :: this

    !> Summation coefficients for the Orbitals
    complex(dp), intent(in) :: eigVecsCmpl(:,:)

    !> Array of k-points
    real(dp), intent(in) :: kPoints(:,:)

    !> Index of the k-points in kPoints for every mol.orbital
    integer, intent(in) :: kIndexes(:)

    !> Molecular orbitals on grid on exit.
    complex(dp), intent(out) :: valueOnGrid(:,:,:,:)

    !> Enable GPU offloading?
    logical, intent(in), optional :: useGpu

    ! Dummy real arrays
    real(dp) :: dummyReal(0,0,0,0)

    call TMolecularOrbital_getValue_cmpl_generic(this, eigVecsCmpl, kPoints, kIndexes, dummyReal, valueOnGrid, useGpu)

  end subroutine TMolecularOrbital_getValue_cmpl


  !> Returns the total charge density on a grid.
  subroutine TMolecularOrbital_getTotalChrg_cmpl(this, eigVecsCmpl, kPoints, kIndexes, valueOnGrid, occupationVec, useGpu)
    !> MolecularOrbital instance
    type(TMolecularOrbital), intent(in) :: this

    !> Summation coefficients for the Orbitals
    complex(dp), intent(in) :: eigVecsCmpl(:,:)

    !> Array of k-points
    real(dp), intent(in) :: kPoints(:,:)

    !> Index of the k-points in kPoints for every mol.orbital
    integer, intent(in) :: kIndexes(:)

    !> Molecular orbitals on grid on exit.
    real(dp), intent(out) :: valueOnGrid(:,:,:,:)

    !> Calculate total charge. Coefficients for each squared state.
    real(dp), intent(in) :: occupationVec(:)

    !> Enable GPU offloading?
    logical, intent(in), optional :: useGpu

    ! Dummy complex arrays
    complex(dp) :: dummyCmplx(0,0,0,0)

    call TMolecularOrbital_getValue_cmpl_generic(this, eigVecsCmpl, kPoints, kIndexes, valueOnGrid, dummyCmplx, useGpu, &
        & occupationVec)

  end subroutine TMolecularOrbital_getTotalChrg_cmpl


  !> Bundles calls to addAtomicDensities, getTotalChrg and the regular molorb to allow for 
  !> Cleaner public interfaces.
  subroutine TMolecularOrbital_getValue_real_generic(this, eigVecsReal, valueOnGrid, useGpu, &
      & addAtomicDensities, occupationVec)

    !> MolecularOrbital instance
    type(TMolecularOrbital), intent(in) :: this

    !> Summation coefficients for the Orbitals
    real(dp), intent(in) :: eigVecsReal(:,:)

    !> Molecular orbitals on a grid
    real(dp), intent(out) :: valueOnGrid(:,:,:,:)

    !> Enable GPU offloading?
    logical, intent(in), optional :: useGpu

    !> Add densities instead of wave functions
    logical, intent(in), optional :: addAtomicDensities

    !> if present, calculate total charge. Coefficients for each squared state
    real(dp), intent(in), optional :: occupationVec(:)

    ! Empty complex arrays
    integer :: kIndexes(0)
    complex(dp) :: valueCmpl(0, 0, 0, 0), eigVecsCmpl(0, 0), phases(0,0)

    logical, parameter :: isRealInput = .true.
    type(TCalculationContext) :: ctx

    ctx = bundleFlags(isRealInput, addAtomicDensities, useGpu, occupationVec)

    @:ASSERT(this%isInitialised)
    @:ASSERT(all(shape(valueOnGrid) > [0, 0, 0, 0]))
    print *, this%system%nOrb, "orbitals,", size(eigVecsReal, dim=1), "states"
    @:ASSERT(size(eigVecsReal, dim=1) == this%system%nOrb)
    @:ASSERT(.not. (ctx%calcAtomicDensity .and. ctx%calcTotalChrg))

    if(ctx%calcTotalChrg) then
      @:ASSERT(size(occupationVec) == size(eigVecsReal, dim=2))
      @:ASSERT(size(valueOnGrid, dim=4) == 1)
    else
      @:ASSERT(size(eigVecsReal, dim=2) == size(valueOnGrid, dim=4))
    end if

    call evaluateParallel(this%system, this%periodic, kIndexes, phases, this%orbitals, &
        & ctx, eigVecsReal, eigVecsCmpl, valueOnGrid, valueCmpl, occupationVec)

  end subroutine TMolecularOrbital_getValue_real_generic


  !> Bundles calls to getValue_cmpl and getTotalChrg_cmpl to allow for cleaner public interfaces.
  subroutine TMolecularOrbital_getValue_cmpl_generic(this, eigVecsCmpl, kPoints, kIndexes, &
      &  valueOutReal, valueOutCmplx, useGpu, occupationVec)

    !> MolecularOrbital instance
    type(TMolecularOrbital), intent(in) :: this

    !> Summation coefficients for the Orbitals
    complex(dp), intent(in) :: eigVecsCmpl(:,:)

    !> Array of k-points
    real(dp), intent(in) :: kPoints(:,:)

    !> Index of the k-points in kPoints for every mol.orbital
    integer, intent(in) :: kIndexes(:)

    !> Density output grid (if calcTotalChrg)
    real(dp), intent(out) :: valueOutReal(:,:,:,:)

    !> Complex Molecular orbital output grid
    complex(dp), intent(out) :: valueOutCmplx(:,:,:,:)

    !> Enable GPU offloading?
    logical, intent(in), optional :: useGpu

    !> Calculate total charge. Coefficients for each squared state
    real(dp), intent(in), optional :: occupationVec(:)

    ! Dummy real arrays
    real(dp) :: eigVecsReal(0,0)

    complex(dp), allocatable :: phases(:,:)
    logical, parameter :: addAtomicDensities = .false.
    logical, parameter :: isRealInput = .false.
    type(TCalculationContext) :: ctx

    ctx = bundleFlags(isRealInput, addAtomicDensities, useGpu, occupationVec)

    @:ASSERT(this%isInitialised)
    @:ASSERT(.not. (ctx%calcAtomicDensity .and. ctx%calcTotalChrg))
    @:ASSERT(size(eigVecsCmpl, dim=1) == this%system%nOrb)
    @:ASSERT(size(kPoints, dim=1) == 3)
    @:ASSERT(size(kPoints, dim=2) > 0)
    @:ASSERT(size(kIndexes) == size(eigVecsCmpl, dim=2))
    @:ASSERT(maxval(kIndexes) <= size(kPoints, dim=2))
    @:ASSERT(minval(kIndexes) > 0)
    if(ctx%calcTotalChrg) then
      @:ASSERT(all(shape(valueOutReal) > [0, 0, 0, 0]))
      @:ASSERT(size(occupationVec) == size(eigVecsCmpl, dim=2))
      @:ASSERT(size(valueOutReal, dim=4) == 1)
    else
      @:ASSERT(all(shape(valueOutCmplx) > [0, 0, 0, 0]))
      @:ASSERT(size(eigVecsCmpl, dim=2) == size(valueOutCmplx, dim=4))
    end if

    allocate(phases(this%periodic%nCell, size(kPoints, dim =2)))
    if (this%periodic%isPeriodic) then
      phases(:,:) = exp(imag * matmul(transpose(this%periodic%fCellVec), kPoints))
    else
      phases(1,:) = (1.0_dp, 0.0_dp)
    end if

    call evaluateParallel(this%system, this%periodic, kIndexes, phases, this%orbitals, &
      & ctx, eigVecsReal, eigVecsCmpl, valueOutReal, valueOutCmplx, occupationVec)

  end subroutine TMolecularOrbital_getValue_cmpl_generic



end module dftbp_wavegrid_molorb
