!--------------------------------------------------------------------------------------------------!
!  DFTB+: general package for performing fast atomistic simulations                                !
!  Copyright (C) 2006 - 2023  DFTB+ developers group                                               !
!                                                                                                  !
!  See the LICENSE file for terms of usage and distribution.                                       !
!--------------------------------------------------------------------------------------------------!

#:include 'common.fypp'

!> Data types for linear response routines.
module dftbp_timedep_linresptypes
  use dftbp_common_accuracy, only : dp
  implicit none

  private
  public :: linrespSolverTypes, TLinResp, TCasidaParameter, TCasidaParameter_init

  !> Types for solution of RPA equations.
  type TSolverTypesEnum
    integer :: None = 0
    integer :: arpack = 1
    integer :: stratmann = 2
  end type TSolverTypesEnum

  !> Actual values for elecSolverTypes
  type(TSolverTypesEnum), parameter :: linrespSolverTypes = TSolverTypesEnum()


  !> Data type for linear response internal settings.
  type :: TLinResp

    !> Number of excitations to be calculated
    integer :: nExc

    !> State of interest (< 0 find brightest, 0 calculate all nexc states, > 0 that specific state)
    integer :: nStat

    !> Is an energy window specified
    logical :: tEnergyWindow

    !> Energy window for transitions above energy of nexc-th single particle transtion
    real(dp) :: energyWindow

    !> Is an oscillator window specified
    logical :: tOscillatorWindow

    !> The window for transitions not included in nexc and energy window (if used)
    real(dp) :: oscillatorWindow

    !> Onsite corrections (if in use)
    real(dp), allocatable :: onSiteMatrixElements(:,:,:,:)

    !> Should occ-vir transition charges be cached or evaluated on the fly?
    logical :: tCacheChargesOccVir

    !> Same for occ-occ/vir-vir transitions
    logical :: tCacheChargesSame

    !> Number of atoms
    integer :: nAtom

    !> Number of electrons in system
    real(dp) :: nEl

    !> Symmetry required singlet ('S'), triplet ("T") or both ("B")
    character :: symmetry

    !> Is the ground state spin polarized
    logical :: tSpin

    !> Ground state spin constants for each species
    real(dp), allocatable :: spinW(:)

    !> Ground state Hubbard U values for each species
    real(dp), allocatable :: HubbardU(:)

    !> Whether X+Y data should be written
    logical :: writeXplusY = .false.

    !> Should non-adiabatic couplings be computed?
    logical :: tNaCoupling = .false.

    !> Should CI be optimized
    logical :: isCIopt

    !> Energy shift used in CI optimizer
    real(dp) :: energyShiftCI

    !> Initial and final state for non-adiabatic coupling evaluation
    integer :: indNACouplings(2)

    !> Whether coefficients for the excited states should be written to disc
    logical :: writeCoeffs = .false.

    !> Add the ground state to the excited state transition density matrix when determining the
    !! natural orbitals
    logical :: tGrndState = .true.

    !> Whether excitation information should be written
    logical :: writeExcitations = .false.

    !> Whether excited Mulliken populations should be written
    logical :: writeMulliken = .false.

    !> Whether single particle (KS) transitions should be written
    logical :: writeTrans = .false.

    !> Whether single particle (KS) transition charges should be written
    logical :: writeTransQ = .false.

    !> Whether for single particle transition dipole strengths should be written
    logical :: writeSPTrans = .false.

    !> File handle for excitation energies
    logical :: writeExc = .false.

    !> Whether transition dipole data should be written
    logical :: writeTransDip = .false.

    !> For calculations where the geometry changes, previous vectors for restarting the iterative
    !! eigensolver. Note: in the case of ARPACK this is the residual not the eigenvectors
    real(dp), allocatable :: oldEigenVectors(:,:)

    !> Should the density matrix be stored to disc?
    logical :: tWriteDensityMatrix

    ! Solver related

    !> Which solver should be used for the RPA equations?
    integer :: iLinRespSolver = linrespSolverTypes%None

    ! ARPACK related

    !> Write state of Arnoldi solver to disc
    logical :: tArnoldi

    !> Whether Arnoldi solver tests should be made (with results written to file)
    logical :: testArnoldi = .false.

    ! Stratmann related

    !> Subspace dimension factor Stratmann diagonaliser
    integer :: subSpaceFactorStratmann

    !> Whether the NACV file should be written
    logical :: writeNacv  = .false.

    ! Data structure related

    !> Is the data structure initialised?
    logical :: tInit = .false.

  end type TLinResp

  !> Data type for common parameters of the Casida routines known only at runtime 
  type :: TCasidaParameter 

    !> Occupied orbitals per spin channel
    integer, allocatable :: nocc_ud(:)

    !> Virtual orbitals per spin channel
    integer, allocatable :: nvir_ud(:)

    !> Number of occ-occ transitions per spin channel
    integer, allocatable :: nxoo_ud(:)

    !> Number of vir-vir transitions per spin channel
    integer, allocatable :: nxvv_ud(:)

    !> Number of occ-vir transitions per spin channel
    integer, allocatable :: nxov_ud(:)

    !> Number of occupied-virtual transitions (possibly reduced by windowing)
    integer :: nxov_rd

    !> array from pairs of single particles states to compound index
    integer, allocatable :: iaTrans(:,:,:)

    !> Index array for occ-vir single particle excitations
    integer, allocatable :: getIA(:,:)

    !> Index array for occ-occ single particle excitations
    integer, allocatable :: getIJ(:,:)

    !> Index array for vir-vir single particle excitations
    integer, allocatable :: getAB(:,:)

    !> Index array for single particle transitions
    integer, allocatable :: win(:)
    
    !> Single particle excitation energies
    real(dp), allocatable :: wij(:)
    
    ! Square root of occupation difference between vir and occ states
    real(dp), allocatable :: sqrOccIA(:)

    !> Is calculation range-separated?
    logical :: tHybridXc

    !> Is a Z-vector calculation required?
    logical :: tZVector
        
    end type TCasidaParameter

  contains

    !> Initialize the internal data type for the Casida calculations.
    subroutine TCasidaParameter_init(this, nocc_ud, nvir_ud, nxoo_ud, nxvv_ud, nxov_ud, nxov_rd,&
        & iaTrans, getIA, getIJ, getAB, win, wij, sqrOccIA, tHybridXc, tZVector)
    
      !> Run time parameters of the Casida routine
      type(TCasidaParameter), intent(out) :: this

      !> Occupied orbitals per spin channel
      integer, allocatable, intent(inout) :: nocc_ud(:)

      !> Virtual orbitals per spin channel
      integer, allocatable, intent(inout) :: nvir_ud(:)

      !> Number of occ-occ transitions per spin channel
      integer, allocatable, intent(inout) :: nxoo_ud(:)

      !> Number of vir-vir transitions per spin channel
      integer, allocatable, intent(inout) :: nxvv_ud(:)

      !> Number of occ-vir transitions per spin channel
      integer, allocatable, intent(inout) :: nxov_ud(:)

      !> Number of occupied-virtual transitions (possibly reduced by windowing)
      integer, intent(in) :: nxov_rd

      !> array from pairs of single particles states to compound index
      integer, allocatable, intent(inout) :: iaTrans(:,:,:)

      !> Index array for occ-vir single particle excitations
      integer, allocatable, intent(inout) :: getIA(:,:)

      !> Index array for occ-occ single particle excitations
      integer, allocatable, intent(inout) :: getIJ(:,:)

      !> Index array for vir-vir single particle excitations
      integer, allocatable, intent(inout) :: getAB(:,:)

      !> Index array for single particle transitions
      integer, allocatable, intent(inout) :: win(:)
    
      !> Single particle excitation energies
      real(dp), allocatable, intent(inout) :: wij(:)

      ! Square root of occupation difference between vir and occ states
      real(dp), allocatable, intent(inout) :: sqrOccIA(:)

      !> Is calculation range-separated?
      logical, intent(in) :: tHybridXc

      !> Is a Z-vector calculation required?
      logical, intent(in) :: tZVector     

      call move_alloc(nocc_ud, this%nocc_ud)
      call move_alloc(nvir_ud, this%nvir_ud)
      call move_alloc(nxoo_ud, this%nxoo_ud)
      call move_alloc(nxvv_ud, this%nxvv_ud)
      call move_alloc(nxov_ud, this%nxov_ud)
      this%nxov_rd = nxov_rd
      call move_alloc(iaTrans, this%iaTrans)
      call move_alloc(getIA, this%getIA)
      call move_alloc(getIJ, this%getIJ)
      call move_alloc(getAB, this%getAB)
      call move_alloc(win, this%win)
      
      ! Only windowed set of transitions is actually needed in what follows
      allocate(this%wij(nxov_rd))
      this%wij = wij(:nxov_rd)
      deallocate(wij)
      
      allocate(this%sqrOccIA(nxov_rd))
      this%sqrOccIA = sqrOccIA(:nxov_rd)
      deallocate(sqrOccIA)

      this%tHybridXc = tHybridXc
      this%tZVector = tZVector
      
    end subroutine TCasidaParameter_init

end module dftbp_timedep_linresptypes
