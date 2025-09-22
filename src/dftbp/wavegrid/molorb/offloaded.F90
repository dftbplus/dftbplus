!--------------------------------------------------------------------------------------------------!
!  DFTB+: general package for performing fast atomistic simulations                                !
!  Copyright (C) 2006 - 2025  DFTB+ developers group                                               !
!                                                                                                  !
!  See the LICENSE file for terms of usage and distribution.                                       !
!--------------------------------------------------------------------------------------------------!

#:include 'common.fypp'

!> Contains <evaluateCuda>, an Interface to call the C bound Cuda implementation.
module dftbp_wavegrid_molorb_offloaded
  use, intrinsic :: iso_c_binding, only : c_bool, c_double, c_int, c_loc, c_ptr
  use dftbp_wavegrid_molorb_types, only : TCalculationContext, TPeriodicParams, TSystemParams
  use dftbp_wavegrid_slater, only : TSlaterOrbital
  use dftbp_common_accuracy, only : dp
  implicit none
  private

#:if WITH_CUDA
  public :: evaluateCuda
#:endif
  !> Data for the basis set in SoA format
  type TBasisParams
    integer :: nStos

    logical :: useRadialLut
    integer :: nLutPoints
    real(dp) :: invLutStep
    real(dp), allocatable :: lutGridValues(:, :)

    !! SoA
    integer :: maxNPows
    integer :: maxNAlphas
    integer, allocatable :: angMoms(:)
    real(dp), allocatable :: cutoffsSq(:)
    integer, allocatable :: nPows(:)
    integer, allocatable :: nAlphas(:)
    real(dp), allocatable :: coeffs(:,:,:)
    real(dp), allocatable :: alphas(:,:)

    logical :: isInitialized = .false.
  end type TBasisParams
  
  !> C bound structs for parameter passing.
  !> These *must* match the struct definitions in kernel.cuh.
  type, bind(c) :: TGridParamsC
    integer(c_int) :: nPointsX, nPointsY, nPointsZ
    type(c_ptr) :: origin, gridVecs
  end type

  type, bind(c) :: TSystemParamsC
    integer(c_int) :: nAtom, nCell, nSpecies, nOrb
    type(c_ptr) :: coords, species, iStos
  end type

  type, bind(c) :: TPeriodicParamsC
    logical(c_bool) :: isPeriodic
    type(c_ptr) :: latVecs, recVecs2pi, kIndexes, phases
  end type

  type, bind(c) :: TSlaterOrbitalC
    logical(c_bool) :: useRadialLut
    integer(c_int) :: nStos, nLutPoints
    real(c_double) :: inverseLutStep
    type(c_ptr) :: lutGridValues

    integer(c_int) :: maxNPows, maxNAlphas
    type(c_ptr) :: angMoms, nPows, nAlphas
    type(c_ptr) :: cutoffsSq, coeffs, alphas
  end type

  type, bind(c) :: TCalculationParamsC
    logical(c_bool) :: isRealInput, isRealOutput, calcAtomicDensity, calcTotalChrg
    integer(c_int) :: nEigIn, nEigOut
    type(c_ptr) :: eigVecsReal, eigVecsCmpl
    type(c_ptr) :: valueReal_out, valueCmpl_out
  end type


contains
#:if WITH_CUDA
  subroutine evaluateCuda(system, stos, periodic, kIndexes, phases, ctx, &
      & eigVecsReal, eigVecsCmpl, valueReal, valueCmpl)

    !> System
    type(TSystemParams), intent(in), target :: system
    !> Basis set
    type(TSlaterOrbital), intent(in), target :: stos(:)
    !> Periodic boundary conditions
    type(TPeriodicParams), intent(in), target :: periodic
    integer, intent(in), target :: kIndexes(:)
    complex(dp), intent(in), target :: phases(:, :)
    !> Calculation flags
    type(TCalculationContext), intent(in) :: ctx
    !> Eigenvectors
    real(dp), intent(in), target :: eigVecsReal(:, :)
    complex(dp), intent(in), target :: eigVecsCmpl(:, :)
    !> Output grids
    real(dp), intent(out), target :: valueReal(:, :, :, :)
    complex(dp), intent(out), target :: valueCmpl(:, :, :, :)

    interface
      subroutine evaluate_on_device_c(grid, system, periodic, basis, calc) bind(C, name ='evaluate_on_device_c')
        import
        type(TGridParamsC), intent(in) :: grid
        type(TSystemParamsC), intent(in) :: system
        type(TPeriodicParamsC), intent(in) :: periodic
        type(TSlaterOrbitalC), intent(in) :: basis
        type(TCalculationParamsC), intent(in) :: calc
      end subroutine evaluate_on_device_c
    end interface

    type(TBasisParams), target :: basis

    type(TGridParamsC) :: grid_p
    type(TSystemParamsC) :: system_p
    type(TPeriodicParamsC) :: periodic_p
    type(TSlaterOrbitalC) :: basis_p
    type(TCalculationParamsC) :: calc_p

    call prepareBasisSet(basis, stos)

    ! Output grid description
    if (ctx%isRealOutput) then
      grid_p%nPointsX = size(valueReal, dim=1)
      grid_p%nPointsY = size(valueReal, dim=2)
      grid_p%nPointsZ = size(valueReal, dim=3)
    else
      grid_p%nPointsX = size(valueCmpl, dim=1)
      grid_p%nPointsY = size(valueCmpl, dim=2)
      grid_p%nPointsZ = size(valueCmpl, dim=3)
    end if
    grid_p%origin = c_loc(system%origin)
    grid_p%gridVecs = c_loc(system%gridVecs)
  
    ! System setup
    system_p%nAtom = system%nAtom
    system_p%nCell = size(system%coords, dim=3)
    system_p%nSpecies = system%nSpecies
    system_p%nOrb = system%nOrb
    system_p%coords = c_loc(system%coords)
    system_p%species = c_loc(system%species)
    system_p%iStos = c_loc(system%iStos)
    
    ! Periodic boundary conditions
    periodic_p%isPeriodic = periodic%isPeriodic
    periodic_p%latVecs = c_loc(periodic%latVecs)
    periodic_p%recVecs2pi = c_loc(periodic%recVecs2pi)
    periodic_p%kIndexes = c_loc(kIndexes)
    periodic_p%phases = c_loc(phases)
    
    ! Basis set
    basis_p%nStos = basis%nStos
    basis_p%angMoms = c_loc(basis%angMoms)
    basis_p%cutoffsSq = c_loc(basis%cutoffsSq)
    basis_p%useRadialLut = basis%useRadialLut

    if (basis%useRadialLut) then
      basis_p%nLutPoints = basis%nLutPoints
      basis_p%inverseLutStep = basis%invLutStep
      basis_p%lutGridValues = c_loc(basis%lutGridValues)
    else
      basis_p%maxNPows = basis%maxNPows
      basis_p%maxNAlphas = basis%maxNAlphas
      basis_p%nPows = c_loc(basis%nPows)
      basis_p%nAlphas = c_loc(basis%nAlphas)
      basis_p%coeffs = c_loc(basis%coeffs)
      basis_p%alphas = c_loc(basis%alphas)
    end if

    ! Calculation params
    if (ctx%isRealInput) then
      calc_p%nEigIn = size(eigVecsReal, dim=2)
    else
      calc_p%nEigIn = size(eigVecsCmpl, dim=2)
    end if
    if (ctx%isRealOutput) then
      calc_p%nEigOut = size(valueReal, dim=4)
    else
      calc_p%nEigOut = size(valueCmpl, dim=4)
    end if
    if (ctx%calcTotalChrg) then
      @:ASSERT(calc_p%nEigOut == 1)
    end if
    calc_p%isRealInput = ctx%isRealInput
    calc_p%isRealOutput = ctx%isRealOutput
    calc_p%calcAtomicDensity = ctx%calcAtomicDensity
    calc_p%calcTotalChrg = ctx%calcTotalChrg
    calc_p%eigVecsReal = c_loc(eigVecsReal)
    calc_p%eigVecsCmpl = c_loc(eigVecsCmpl)
    calc_p%valueReal_out = c_loc(valueReal)
    calc_p%valueCmpl_out = c_loc(valueCmpl)

    call evaluate_on_device_c(grid_p, system_p, periodic_p, basis_p, calc_p)

  end subroutine evaluateCuda



 
  !> Convert the basis set to SoA format or unified Lut table
  !> Currently, mixed lut/direct calculation is not supported on the GPU.
  !> Additionally, all orbitals must use identical LUT settings.
  subroutine prepareBasisSet(this, stos)
    type(TBasisParams), intent(out) :: this
    type(TSlaterOrbital), intent(in) :: stos(:)
    integer :: iOrb

    this%nStos = size(stos)
    ! We will assert that all orbitals use identical LUT settings
    this%useRadialLut = stos(1)%useRadialLut
    this%nLutPoints = stos(1)%nGrid
    this%invLutStep = stos(1)%invLutStep
    allocate(this%angMoms(this%nStos))
    allocate(this%cutoffsSq(this%nStos))

    do iOrb = 1, this%nStos
      this%angMoms(iOrb) = stos(iOrb)%angMom
      this%cutoffsSq(iOrb) = stos(iOrb)%cutoffSq
    end do

    if (this%useRadialLut) then
      print *, "Using radial LUT for STO evaluation with ", this%nLutPoints, " points."
      allocate(this%lutGridValues(this%nLutPoints, this%nStos))
      do iOrb = 1, this%nStos
        @:ASSERT(stos(iOrb)%useRadialLut)
        @:ASSERT(stos(iOrb)%nGrid == this%nLutPoints)
        @:ASSERT(abs(stos(iOrb)%invLutStep - this%invLutStep) < 1.0e-12_dp)

        this%lutGridValues(:, iOrb) = stos(iOrb)%gridValue
      end do
    else ! Direct evaluation
      ! Allocate SoA arrays
      allocate(this%nPows(this%nStos))
      allocate(this%nAlphas(this%nStos))

      ! Populate SoA arrays
      do iOrb = 1, this%nStos
        @:ASSERT(.not. stos(iOrb)%useRadialLut)
        this%nPows(iOrb) = stos(iOrb)%nPow
        this%nAlphas(iOrb) = stos(iOrb)%nAlpha
      end do
      this%maxNPows = maxval(this%nPows)
      this%maxNAlphas = maxval(this%nAlphas)
      !print *, "nStos=", this%nStos, " maxNPows=", this%maxNPows, " maxNAlphas=", this%maxNAlphas

      ! Allocate and populate coefficient/alpha matrices
      allocate(this%coeffs(this%maxNPows, this%maxNAlphas, this%nStos))
      allocate(this%alphas(this%maxNAlphas, this%nStos))
      do iOrb = 1, this%nStos
        this%coeffs(1:stos(iOrb)%nPow, 1:stos(iOrb)%nAlpha, iOrb) = stos(iOrb)%aa
        this%alphas(1:stos(iOrb)%nAlpha, iOrb) = stos(iOrb)%alpha
      end do
    end if

    this%isInitialized = .true.
  end subroutine prepareBasisSet
#:endif

end module dftbp_wavegrid_molorb_offloaded
