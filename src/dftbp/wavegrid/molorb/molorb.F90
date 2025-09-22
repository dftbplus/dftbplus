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
  use dftbp_wavegrid_slater, only : TSlaterOrbital
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
    type(TSlaterOrbital), allocatable :: stos(:)
    !> Boundary conditions handler for coordinate recalculation
    type(TBoundaryConds) :: boundaryCond

    logical :: isInitialised = .false.
  contains
    private
    procedure, public :: init => TMolecularOrbital_init
    procedure, public :: updateCoords => TMolecularOrbital_updateCoords

    procedure :: initSpeciesMapping
    procedure :: flattenBasis
    procedure :: initPeriodic
  end type TMolecularOrbital

  !> Data type containing information about the basis for a species.
  type TSpeciesBasis
    !> Atomic number of the species
    integer :: atomicNumber
    !> Nr. of orbitals
    integer :: nOrb
    !> STO for each orbital
    type(TSlaterOrbital), allocatable :: stos(:)
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
    type(TGeometry), intent(in) :: geometry
    type(TBoundaryConds), intent(in) :: boundaryCond
    type(TSpeciesBasis), intent(in) :: basisInput(:)
    real(dp), intent(in) :: origin(3)
    real(dp), intent(in) :: gridVecs(3,3)

    @:ASSERT(.not. this%isInitialised)
    @:ASSERT(size(origin) == 3)
    @:ASSERT(all(shape(gridVecs) == [3, 3]))
    this%boundaryCond = boundaryCond
    this%system%origin = origin
    this%system%gridVecs = gridVecs
    
    call this%initSpeciesMapping(geometry, basisInput)
    call this%flattenBasis(basisInput)
    call this%initPeriodic(geometry, maxCutoff(this%stos))
    call this%updateCoords(geometry)

    this%isInitialised = .true.
  end subroutine TMolecularOrbital_init


  !> Initialises non-geometric system parameters 
  !! Counts total number of orbitals and STOs and creates index maps.
  subroutine initSpeciesMapping(this, geometry, basis)
    class(TMolecularOrbital), intent(inout) :: this
    type(TGeometry), intent(in) :: geometry
    type(TSpeciesBasis), intent(in) :: basis(:)

    integer :: nOrbTotal, iSpec, iAtom, ind, iOrb, angMom

    @:ASSERT(geometry%nSpecies == size(basis))

    this%system%nAtom = geometry%nAtom
    this%system%nSpecies = geometry%nSpecies
    allocate(this%system%species(this%system%nAtom))
    this%system%species(:) = geometry%species


    ! Create STO index map: iStos(i) points to the first STO of species i
    allocate(this%system%iStos(this%system%nSpecies + 1))
    ind = 1
    do iSpec = 1, this%system%nSpecies
      this%system%iStos(iSpec) = ind
      ind = ind + basis(iSpec)%nOrb
    end do
    this%system%iStos(this%system%nSpecies + 1) = ind

    ! Count total number of orbitals (including m-dependence)
    nOrbTotal = 0
    do iAtom = 1, this%system%nAtom
      iSpec = this%system%species(iAtom)
      do iOrb = 1, basis(iSpec)%nOrb
        angMom = basis(iSpec)%stos(iOrb)%angMom
        nOrbTotal = nOrbTotal + 1 + 2 * angMom
      end do
    end do
    this%system%nOrb = nOrbTotal

    this%system%speciesInitialised = .true.
  end subroutine initSpeciesMapping


  !> Flattens the basis set into an array of STOs. 
  subroutine flattenBasis(this, basisInput)
    class(TMolecularOrbital), intent(inout) :: this
    type(TSpeciesBasis), intent(in) :: basisInput(:)
    integer :: iSpec, ind, nStos
    nStos = sum(basisInput(:)%nOrb)

    ! Flatten the basis array
    allocate(this%stos(nStos))
    ind = 1
    do iSpec = 1, this%system%nSpecies
      this%stos(ind:ind+basisInput(iSpec)%nOrb-1) = basisInput(iSpec)%stos(:)
      ind = ind + basisInput(iSpec)%nOrb
    end do
  end subroutine flattenBasis


  !> Initializes periodic parameters
  subroutine initPeriodic(this, geometry, maxCutoff)
    class(TMolecularOrbital), intent(inout) :: this
    type(TGeometry), intent(in) :: geometry
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
  end subroutine initPeriodic


  !> Initializes coordinates including periodic images.
  ! Depends on initPeriodic having run first.
  subroutine TMolecularOrbital_updateCoords(this, geometry)
      class(TMolecularOrbital), intent(inout) :: this
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
  function bundleFlags(isRealInput, addAtomicDensities, useGPU, occupationVec) result(ctx)
    logical, intent(in) :: isRealInput
    logical, intent(in), optional :: addAtomicDensities, useGPU
    real(dp), intent(in), optional :: occupationVec(:)
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
    if (present(useGPU)) then
      #:if WITH_CUDA
        ctx%runOnGPU = useGPU
      #:else
      if (useGPU) then
        call error("GPU offloading requested (useGPU=.true.), but not available in this build. (missing WITH_CUDA)")
      end if
      #:endif
    end if

    ctx%isRealOutput = ctx%isRealInput .or. ctx%calcTotalChrg

  end function bundleFlags


  !> Returns molecular orbitals on a real grid. 
  subroutine TMolecularOrbital_getValue_real(this, eigVecsReal, valueOnGrid, useGPU)
    !> MolecularOrbital data instance
    type(TMolecularOrbital), intent(in) :: this
    !> Summation coefficients for the STOs
    real(dp), intent(in) :: eigVecsReal(:,:)
    !> Molecular orbitals on a grid
    real(dp), intent(out) :: valueOnGrid(:,:,:,:)
    !> Enable GPU offloading?
    logical, intent(in), optional :: useGPU
  
    call TMolecularOrbital_getValue_real_generic(this, eigVecsReal, valueOnGrid, useGPU)
  end subroutine TMolecularOrbital_getValue_real


  !> Returns the total charge density on a grid.
  !> This squares each state and sums them up weighted by occupationVec.
  subroutine TMolecularOrbital_getTotalChrg_real(this, eigVecsReal, valueOnGrid, occupationVec, useGPU)
    !> MolecularOrbital instance
    type(TMolecularOrbital), intent(in) :: this
    !> Summation coefficients for the STOs
    real(dp), intent(in) :: eigVecsReal(:,:)
    !> Molecular orbitals on a grid
    real(dp), intent(out) :: valueOnGrid(:,:,:,:)
    !> Calculate total charge. Coefficients for each squared state.
    real(dp), intent(in) :: occupationVec(:)
    !> Enable GPU offloading?
    logical, intent(in), optional :: useGPU
  
    call TMolecularOrbital_getValue_real_generic(this, eigVecsReal, valueOnGrid, useGPU, &
        & addAtomicDensities=.false., occupationVec=occupationVec)
  end subroutine TMolecularOrbital_getTotalChrg_real


  !> Calculates the atomic densities by squaring each STO *before* summation.
  subroutine TMolecularOrbital_getAtomicDensities_real(this, eigVecsReal, valueOnGrid, useGPU)
    !> MolecularOrbital instance
    type(TMolecularOrbital), intent(in) :: this
    !> Summation coefficients for the STOs
    real(dp), intent(in) :: eigVecsReal(:,:)
    !> Molecular orbitals on a grid
    real(dp), intent(out) :: valueOnGrid(:,:,:,:)
    !> Enable GPU offloading?
    logical, intent(in), optional :: useGPU

    call TMolecularOrbital_getValue_real_generic(this, eigVecsReal, valueOnGrid, useGPU, &
        & addAtomicDensities=.true.)
  end subroutine TMolecularOrbital_getAtomicDensities_real


  !> Returns molecular orbitals on a complex grid.
  subroutine TMolecularOrbital_getValue_cmpl(this, eigVecsCmpl, kPoints, kIndexes, valueOnGrid, useGPU)
    !> MolecularOrbital instance
    type(TMolecularOrbital), intent(in) :: this
    !> Summation coefficients for the STOs
    complex(dp), intent(in) :: eigVecsCmpl(:,:)
    !> Array of k-points
    real(dp), intent(in) :: kPoints(:,:)
    !> Index of the k-points in kPoints for every mol.orbital
    integer, intent(in) :: kIndexes(:)
    !> Molecular orbitals on grid on exit.
    complex(dp), intent(out) :: valueOnGrid(:,:,:,:)
    !> Enable GPU offloading?
    logical, intent(in), optional :: useGPU
    ! Dummy real arrays
    real(dp) :: dummyReal(0,0,0,0)

    call TMolecularOrbital_getValue_cmpl_generic(this, eigVecsCmpl, kPoints, kIndexes, dummyReal, valueOnGrid, useGPU)

  end subroutine TMolecularOrbital_getValue_cmpl


  !> Returns the total charge density on a grid.
  subroutine TMolecularOrbital_getTotalChrg_cmpl(this, eigVecsCmpl, kPoints, kIndexes, valueOnGrid, occupationVec, useGPU)
    !> MolecularOrbital instance
    type(TMolecularOrbital), intent(in) :: this
    !> Summation coefficients for the STOs
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
    logical, intent(in), optional :: useGPU
    ! Dummy complex arrays
    complex(dp) :: dummyCmplx(0,0,0,0)

    call TMolecularOrbital_getValue_cmpl_generic(this, eigVecsCmpl, kPoints, kIndexes, valueOnGrid, dummyCmplx, useGPU, &
        & occupationVec)

  end subroutine TMolecularOrbital_getTotalChrg_cmpl


  !> Bundles calls to addAtomicDensities, getTotalChrg and the regular molorb to allow for 
  !> Cleaner public interfaces.
  subroutine TMolecularOrbital_getValue_real_generic(this, eigVecsReal, valueOnGrid, useGPU, &
      & addAtomicDensities, occupationVec)

    !> MolecularOrbital instance
    type(TMolecularOrbital), intent(in) :: this
    !> Summation coefficients for the STOs
    real(dp), intent(in) :: eigVecsReal(:,:)
    !> Molecular orbitals on a grid
    real(dp), intent(out) :: valueOnGrid(:,:,:,:)
    !> Enable GPU offloading?
    logical, intent(in), optional :: useGPU
    !> Add densities instead of wave functions
    logical, intent(in), optional :: addAtomicDensities
    !> if present, calculate total charge. Coefficients for each squared state
    real(dp), intent(in), optional :: occupationVec(:)

    ! Empty complex arrays
    integer :: kIndexes(0)
    complex(dp) :: valueCmpl(0, 0, 0, 0), eigVecsCmpl(0, 0), phases(0,0)

    logical, parameter :: isRealInput = .true.
    type(TCalculationContext) :: ctx

    ctx = bundleFlags(isRealInput, addAtomicDensities, useGPU, occupationVec)

    @:ASSERT(this%isInitialised)
    @:ASSERT(all(shape(valueOnGrid) > [0, 0, 0, 0]))
    @:ASSERT(size(eigVecsReal, dim=1) == this%system%nOrb)
    @:ASSERT(.not. (ctx%calcAtomicDensity .and. ctx%calcTotalChrg))

    if(ctx%calcTotalChrg) then
      @:ASSERT(size(occupationVec) == size(eigVecsReal, dim=2))
      @:ASSERT(size(valueOnGrid, dim=4) == 1)
    else
      @:ASSERT(size(eigVecsReal, dim=2) == size(valueOnGrid, dim=4))
    end if

    call evaluateParallel(this%system, this%periodic, kIndexes, phases, this%stos, &
        & ctx, eigVecsReal, eigVecsCmpl, valueOnGrid, valueCmpl, occupationVec)

  end subroutine TMolecularOrbital_getValue_real_generic


  !> Bundles calls to getValue_cmpl and getTotalChrg_cmpl to allow for cleaner public interfaces.
  subroutine TMolecularOrbital_getValue_cmpl_generic(this, eigVecsCmpl, kPoints, kIndexes, &
      &  valueOutReal, valueOutCmplx, useGPU, occupationVec)

    !> MolecularOrbital instance
    type(TMolecularOrbital), intent(in) :: this
    !> Summation coefficients for the STOs
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
    logical, intent(in), optional :: useGPU
    !> Calculate total charge. Coefficients for each squared state
    real(dp), intent(in), optional :: occupationVec(:)
    ! Dummy real arrays
    real(dp) :: eigVecsReal(0,0)

    complex(dp), allocatable :: phases(:,:)
    logical, parameter :: addAtomicDensities = .false.
    logical, parameter :: isRealInput = .false.
    type(TCalculationContext) :: ctx

    ctx = bundleFlags(isRealInput, addAtomicDensities, useGPU, occupationVec)

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

    call evaluateParallel(this%system, this%periodic, kIndexes, phases, this%stos, &
      & ctx, eigVecsReal, eigVecsCmpl, valueOutReal, valueOutCmplx, occupationVec)

  end subroutine TMolecularOrbital_getValue_cmpl_generic


  function maxCutoff(stos) result(maxCut)
    type(TSlaterOrbital), intent(in) :: stos(:)
    real(dp) :: maxCut
    real(dp) :: maxCutSq
    integer :: iSto

    maxCutSq = 0.0_dp
    do iSto = 1, size(stos)
      if (stos(iSto)%cutoffSq > maxCutSq) then
        maxCutSq = stos(iSto)%cutoffSq
      end if
    end do
    maxCut = sqrt(maxCutSq)
  end function maxCutoff

end module dftbp_wavegrid_molorb
