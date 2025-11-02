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
  use dftbp_wavegrid_basis, only : TOrbitalWrapper, TRadialTableOrbital, TRadialTableOrbital_initFromOrbital
  use dftbp_common_accuracy, only : dp
  implicit none
  private

#:if WITH_CUDA
  public :: evaluateCuda
#:endif
  
  !> Default LUT step size if no orbital specifies one.
  real(dp), parameter :: defaultLutStep = 0.01_dp

  !> Data for the basis set in SoA format
  type TBasisParams
    integer :: nOrbitals
    
    integer :: nLutPoints
    real(dp) :: invLutStep
    real(dp), allocatable :: lutGridValues(:, :)

    integer, allocatable :: angMoms(:)
    real(dp), allocatable :: cutoffsSq(:)
  end type TBasisParams
  
  !> C bound structs for parameter passing.
  !> These *must* match the struct definitions in kernel.cuh (including order).
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

  type, bind(c) :: TOrbitalC
    integer(c_int) :: nOrbitals, nLutPoints
    real(c_double) :: inverseLutStep
    type(c_ptr) :: lutGridValues
    type(c_ptr) :: angMoms, cutoffsSq
  end type

  type, bind(c) :: TCalculationParamsC
    logical(c_bool) :: isRealInput, isRealOutput, calcAtomicDensity, calcTotalChrg
    integer(c_int) :: nEigIn, nEigOut
    type(c_ptr) :: eigVecsReal, eigVecsCmpl
    type(c_ptr) :: valueReal_out, valueCmpl_out
  end type

  !> C binding for the evaluation kernel, passing all of the above structs.
  interface
    subroutine evaluate_on_device_c(grid, system, periodic, basis, calc) bind(C, name ='evaluate_on_device_c')
      import
      type(TGridParamsC), intent(in) :: grid
      type(TSystemParamsC), intent(in) :: system
      type(TPeriodicParamsC), intent(in) :: periodic
      type(TOrbitalC), intent(in) :: basis
      type(TCalculationParamsC), intent(in) :: calc
    end subroutine evaluate_on_device_c
  end interface

contains
#:if WITH_CUDA
  subroutine evaluateCuda(system, orbitals, periodic, kIndexes, phases, ctx, &
      & eigVecsReal, eigVecsCmpl, valueReal, valueCmpl)

    !> System
    type(TSystemParams), intent(in), target :: system

    !> Basis set
    type(TOrbitalWrapper), intent(in), target :: orbitals(:)

    !> Periodic boundary conditions
    type(TPeriodicParams), intent(in), target :: periodic

    !> K-point indexes
    integer, intent(in), target :: kIndexes(:)

    !> Phase factors
    complex(dp), intent(in), target :: phases(:, :)

    !> Calculation flags
    type(TCalculationContext), intent(in) :: ctx

    !> Real Eigenvectors
    real(dp), intent(in), target :: eigVecsReal(:, :)

    !> Complex Eigenvectors (if not real)
    complex(dp), intent(in), target :: eigVecsCmpl(:, :)

    !> Real output grid
    real(dp), intent(out), target :: valueReal(:, :, :, :)

    !> Complex output grid (if not real)
    complex(dp), intent(out), target :: valueCmpl(:, :, :, :)

    type(TBasisParams), target :: basis
    type(TGridParamsC) :: grid_p
    type(TSystemParamsC) :: system_p
    type(TPeriodicParamsC) :: periodic_p
    type(TOrbitalC) :: basis_p
    type(TCalculationParamsC) :: calc_p

    call prepareBasisSet(basis, orbitals)

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
    basis_p%nOrbitals = basis%nOrbitals
    basis_p%angMoms = c_loc(basis%angMoms)
    basis_p%cutoffsSq = c_loc(basis%cutoffsSq)
    ! LUT data
    basis_p%nLutPoints = basis%nLutPoints
    basis_p%inverseLutStep = basis%invLutStep
    basis_p%lutGridValues = c_loc(basis%lutGridValues)

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


  !> Convert the basis Set to LUTs (TRadialTable).
  !> Resamples to identical resolution (highest) and cutoff (largest).
  !> Merges all Luts into a single 2D array.
  subroutine prepareBasisSet(this, orbitals)

    !> Output basis set data
    type(TBasisParams), intent(out) :: this

    !> Input orbitals to resample
    type(TOrbitalWrapper), intent(in) :: orbitals(:)

    integer :: iOrb
    real(dp) :: cutoff, resolution
    type(TRadialTableOrbital) :: lut

    ! Determine largest cutoff
    cutoff = 10.0_dp ! sqrt(maxval(orbitals%o%cutoffSq))

    ! Determine finest resolution
    resolution = defaultLutStep
    do iOrb = 1, size(orbitals)
      associate (orb=>orbitals(iOrb)%o)
      select type(orb)
        type is (TRadialTableOrbital)
          resolution = min(resolution, orb%gridDist)
        class default
          ! No resolution specified.
        end select
      end associate
    end do

    this%nOrbitals = size(orbitals)
    this%nLutPoints = floor(cutoff / resolution) + 2
    this%invLutStep = 1.0_dp / resolution

    allocate(this%angMoms(this%nOrbitals))
    allocate(this%cutoffsSq(this%nOrbitals))
    allocate(this%lutGridValues(this%nLutPoints, this%nOrbitals))


    do iOrb = 1, this%nOrbitals
      this%angMoms(iOrb) = orbitals(iOrb)%o%angMom
      this%cutoffsSq(iOrb) = orbitals(iOrb)%o%cutoffSq
      call TRadialTableOrbital_initFromOrbital(lut, orbitals(iOrb)%o, resolution, cutoff)

      @:ASSERT(this%nLutPoints == size(lut%gridValue))
      @:ASSERT(abs(lut%invLutStep - this%invLutStep) < 1.0e-12_dp)

      this%lutGridValues(:, iOrb) = lut%gridValue
      deallocate(lut%gridValue)
    end do

  end subroutine prepareBasisSet
#:endif

end module dftbp_wavegrid_molorb_offloaded
