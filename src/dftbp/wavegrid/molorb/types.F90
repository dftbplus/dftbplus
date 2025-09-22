!--------------------------------------------------------------------------------------------------!
!  DFTB+: general package for performing fast atomistic simulations                                !
!  Copyright (C) 2006 - 2025  DFTB+ developers group                                               !
!                                                                                                  !
!  See the LICENSE file for terms of usage and distribution.                                       !
!--------------------------------------------------------------------------------------------------!

#:include 'common.fypp'

!> Contains private data types for bundling input parameters to the molecular orbital calculator.
module dftbp_wavegrid_molorb_types
  use dftbp_common_accuracy, only : dp

  implicit none

  public

  !> Data for the system geometry and composition
  type TSystemParams
    !! Composition
    integer :: nAtom
    integer :: nOrb
    integer :: nSpecies
    integer, allocatable :: species(:)
    integer, allocatable :: iStos(:)
    logical :: speciesInitialised = .false.
    
    !! Geometry
    real(dp) :: origin(3)
    real(dp) :: gridVecs(3,3)
    real(dp), allocatable :: coords(:,:,:)
    logical :: coordsInitialised = .false.
  end type TSystemParams

  !> Data for periodic boundary conditions
  type TPeriodicParams
    logical :: isPeriodic
    real(dp), allocatable :: latVecs(:,:)
    real(dp), allocatable :: recVecs2pi(:,:)
    real(dp), allocatable :: rCellVec(:,:)
    real(dp), allocatable :: fCellVec(:,:)
    integer :: nCell
    logical :: isInitialized = .false.
  end type TPeriodicParams
  
  !> Control type holding calculation flags
  type :: TCalculationContext
    logical :: isRealInput
    logical :: isRealOutput
    logical :: calcAtomicDensity
    logical :: calcTotalChrg
    logical :: runOnGPU
    logical :: isInitialized = .false.
  end type

end module dftbp_wavegrid_molorb_types
