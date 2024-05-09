!--------------------------------------------------------------------------------------------------!
!  DFTB+: general package for performing fast atomistic simulations                                !
!  Copyright (C) 2006 - 2023  DFTB+ developers group                                               !
!                                                                                                  !
!  See the LICENSE file for terms of usage and distribution.                                       !
!--------------------------------------------------------------------------------------------------!

#:include 'common.fypp'
#:include 'error.fypp'


!> Contains hybrid xc-functional related routines.
module dftbp_dftb_hybridxc
  use dftbp_common_accuracy, only : dp
  use dftbp_common_constants, only : pi, imag
  use dftbp_common_environment, only : TEnvironment, globalTimers
  use dftbp_common_schedule, only : getStartAndEndIndex
  use dftbp_common_status, only : TStatus
  use dftbp_dftb_densitymatrix, only : TDensityMatrix
  use dftbp_dftb_nonscc, only : TNonSccDiff, buildS
  use dftbp_dftb_periodic, only : TNeighbourList, TSymNeighbourList, getCellTranslations, cart2frac
  use dftbp_dftb_rshgamma, only : getCamAnalyticalGammaValue_workhorse,&
      & getHfAnalyticalGammaValue_workhorse, getLrAnalyticalGammaValue_workhorse,&
      & getdHfAnalyticalGammaValue_workhorse, getdLrAnalyticalGammaValue_workhorse,&
      & getddLrNumericalGammaValue_workhorse
  use dftbp_dftb_slakocont, only : TSlakoCont
  use dftbp_dftb_sparse2dense, only : unpackHS
  use dftbp_math_blasroutines, only : gemm, symm, hemm
  use dftbp_math_matrixops, only : adjointLowerTriangle
  use dftbp_math_simplealgebra, only : determinant33
  use dftbp_math_sorting, only : heap_sort, index_heap_sort
  use dftbp_math_wignerseitz, only : generateWignerSeitzGrid
  use dftbp_type_commontypes, only : TOrbitals, TParallelKS
  use dftbp_type_densedescr, only : TDenseDescr
  use dftbp_type_integral, only : TIntegral
  use dftbp_type_wrappedintr, only : TWrappedInt1, TWrappedReal1, TWrappedReal2
#:if WITH_MPI
  use dftbp_extlibs_mpifx, only : mpifx_recv, mpifx_send
#:endif
#:if WITH_SCALAPACK
  use dftbp_extlibs_mpifx, only : MPI_SUM, mpifx_allreduceip, mpifx_bcast
  use dftbp_extlibs_scalapackfx, only : CSRC_, DLEN_, MB_, NB_, RSRC_, pblasfx_pgemm,&
      & pblasfx_ptranc, pblasfx_ptran, pblasfx_phemm, pblasfx_psymm, scalafx_indxl2g,&
      & scalafx_getdescriptor, scalafx_getlocalshape, scalafx_addl2g, linecomm
  use dftbp_math_bisect, only : bisection
#:endif

  implicit none
  private

  public :: THybridXcSKTag, THybridXcFunc, THybridXcFunc_init
  public :: getDirectionalCamGammaPrimeValue
  public :: hybridXcFunc, hybridXcAlgo, hybridXcGammaTypes, checkSupercellFoldingMatrix

#:if WITH_SCALAPACK
  public :: getFullFromDistributed, scatterFullToDistributed
#:endif


  !> Returns the directional derivative of long-range + full-range Hartree-Fock gamma.
  interface getDirectionalCamGammaPrimeValue
    module procedure getDirectionalCamGammaPrimeValue_cluster
    module procedure getDirectionalCamGammaPrimeValue_periodic
  end interface getDirectionalCamGammaPrimeValue


  !> Enumerator for type of hybrid functional used.
  !! (no hybrid, global hybrid, purely long-range corrected, general CAM range-separated form)
  type :: THybridXcFuncEnum

    !> (semi-)local
    integer :: none = 0

    !> Global hybrid
    integer :: hyb = 1

    !> Long-range corrected
    integer :: lc = 2

    !> General Coulomb-attenuated method range-separated hybrid
    integer :: cam = 3

  end type THybridXcFuncEnum


  !> Enumerator for algorithms to build up the Fock-type exchange contribution to the Hamiltonian.
  type :: THybridXcAlgoEnum

    !> None
    integer :: none = 0

    !> Neighbour based
    integer :: neighbourBased = 1

    !> Threshold based
    integer :: thresholdBased = 2

    !> Matrix based
    integer :: matrixBased = 3

  end type THybridXcAlgoEnum


  !> Enumerator for gamma function types.
  type :: THybridXcGammaTypesEnum

    !> None
    integer :: none = 0

    !> Full, unaltered gamma function (full)
    integer :: full = 1

    !> Truncated gamma function with hard cutoff (truncated)
    integer :: truncated = 2

    !> Truncated and poly5zero damped gamma function (truncated+damping)
    integer :: truncatedAndDamped = 3

    !> Minimum image convention, based on density matrix criteria (mic)
    integer :: mic = 4

  end type THybridXcGammaTypesEnum


  !> Container for enumerated types of hybrid functionals.
  type(THybridXcFuncEnum), parameter :: hybridXcFunc = THybridXcFuncEnum()

  !> Container for enumerated range separation algorithms
  type(THybridXcAlgoEnum), parameter :: hybridXcAlgo = THybridXcAlgoEnum()

  !> Container for enumerated range separation gamma function types
  type(THybridXcGammaTypesEnum), parameter :: hybridXcGammaTypes = THybridXcGammaTypesEnum()


  !> Slater-Koster file HybridXc tag structure.
  type :: THybridXcSKTag

    !> Range-separation parameter
    real(dp) :: omega

    !> CAM alpha parameter
    real(dp) :: camAlpha

    !> CAM beta parameter
    real(dp) :: camBeta

  end type THybridXcSKTag


  !> Abstract base type for hybrid functionality, some procedures need to be overridden.
  type, abstract :: THybridXcFunc

    !> Real-space coordinates of atoms (relative units), potentially including periodic images
    real(dp), allocatable :: coords(:,:)

    !> Real-space coordinates of atoms (absolute units), potentially including periodic images
    real(dp), allocatable :: rCoords(:,:)

    !> Evaluated (long-range + HF full-range) gamma of Atom1 and Atom2 (central cell only)
    real(dp), allocatable :: camGammaEval0(:,:), camdGammaEval0(:,:,:)

    !> Evaluated (long-range + HF full-range) gamma in the general k-point case
    type(TWrappedReal1), allocatable :: camGammaEvalG(:,:)
    type(TWrappedReal2), allocatable :: camdGammaEvalG(:,:)

    !> Number of numerically non-zero gamma's
    type(TWrappedInt1), allocatable :: nNonZeroGammaG(:,:)

    !> Range-separation parameter
    real(dp) :: omega

    !> CAM alpha parameter
    real(dp) :: camAlpha

    !> CAM beta parameter
    real(dp) :: camBeta

    !> Species-resolved Hubbard U values
    real(dp), allocatable :: hubbu(:)

    !> Full, sparse overlap as obtained by buildS() for a symmetric neighbour list
    real(dp), allocatable :: overSym(:)

    !> Previous Hamiltonian in screening by tolerance
    real(dp), allocatable :: hPrev(:,:)

    !> Previous delta density matrix in screening by tolerance
    real(dp), allocatable :: dRhoPrev(:,:)

    !> Previous Hamiltonian in screening by tolerance
    complex(dp), allocatable :: hPrevCplxHS(:,:,:)

    !> Previous delta density matrix in screening by tolerance
    real(dp), allocatable :: dRhoPrevCplxHS(:,:,:,:,:,:)

    !> Is screening initialised?
    logical :: tScreeningInited

    !> Threshold for screening by value
    real(dp) :: pScreeningThreshold

    !> Total (long-range + full-range Hartree-Fock) CAM energy
    real(dp) :: camEnergy

    !> Is this a spin restricted (false) or unrestricted (true) calculation?
    logical :: tSpin

    !> Is this a DFTB/REKS calculation (true)?
    logical :: tREKS

    !> Algorithm for RSH-Hamiltonian construction (and, if applicable, force calculation)
    integer :: hybridXcAlg = hybridXcAlgo%none

    !> Hybrid xc-functional type, as extracted from SK-file(s)
    integer :: hybridXcType = hybridXcFunc%none

    !> Species-list of atoms in the central cell
    integer, allocatable :: species0(:)

    !> Cutoff for real-space g-summation
    real(dp) :: gSummationCutoff

    !> Number of unitcells along each supercell folding direction to substract from MIC Wigner-Seitz
    !! cell construction
    integer :: wignerSeitzReduction

    !> Cutoff for truncated Gamma
    real(dp) :: gammaCutoff

    !> Damping distance for Gamma truncation
    real(dp) :: gammaDamping

    !> Value, 1st and 2nd derivative of gamma integral at damping distance
    real(dp), allocatable :: lrGammaAtDamping(:,:), lrdGammaAtDamping(:,:), lrddGammaAtDamping(:,:)

    !> Value, 1st and 2nd derivative of gamma integral at damping distance
    real(dp), allocatable :: hfGammaAtDamping(:,:), hfdGammaAtDamping(:,:), hfddGammaAtDamping(:,:)

    !> Overlap estimates
    type(TWrappedReal1), allocatable :: squareOverEst(:)

    !> Descending neighbour indices in terms of overlap estimates
    type(TWrappedInt1), allocatable :: overlapIndices(:)

    !> The k-point compatible BvK real-space shifts in relative coordinates (units of latVecs)
    real(dp), allocatable :: bvKShifts(:,:)

    !> Supercell folding coefficients (diagonal elements)
    integer, allocatable :: coeffsDiag(:)

    !> Gamma function type (mostly for periodic cases; 'full' for non-periodic systems)
    integer :: gammaType = hybridXcGammaTypes%none

    !> Wigner-Seitz grid points in units of lattice vectors
    integer, allocatable :: wsVectors(:,:)

    !> Translation vectors to lattice cells in units of lattice constants
    real(dp), allocatable :: cellVecsG(:,:)

    !> Vectors to unit cells in absolute units
    real(dp), allocatable :: rCellVecsG(:,:)

  contains

    procedure :: updateCoords_cluster, updateCoords_gamma, updateCoords_kpts

    procedure :: foldToBvK => THybridXcFunc_foldToBvK
    procedure :: foldToBvKIndex => THybridXcFunc_foldToBvKIndex

    procedure :: addCamHamiltonian_real
    procedure :: getCamHamiltonian_kpts

    procedure :: addCamHamiltonianMatrix_cluster_cmplx

    procedure :: addHybridEnergy_real
    procedure :: addHybridEnergy_kpts

    procedure :: tabulateCamdGammaEval0_cluster
    procedure :: tabulateCamdGammaEval0_gamma

    procedure :: addCamGradients_real
    procedure :: addCamGradients_kpts

    procedure :: getCentralCellSpecies
    procedure :: getCamGammaCluster
    procedure :: getCamGammaDerivCluster

    procedure(gammaFunc), deferred :: getLrGammaValue
    procedure(gammaFunc), deferred :: getLrGammaPrimeValue

    procedure(gammaFunc), deferred :: getHfGammaValue
    procedure(gammaFunc), deferred :: getHfGammaPrimeValue

  end type THybridXcFunc


  abstract interface
    !> Calculates analytical long-range or full-range HF gamma (derivative) of given type.
    function gammaFunc(this, iSp1, iSp2, dist) result(gamma)

      import :: THybridXcFunc, dp

      !> Class instance
      class(THybridXcFunc), intent(in) :: this

      !> First species
      integer, intent(in) :: iSp1

      !> Second species
      integer, intent(in) :: iSp2

      !> Distance between atoms
      real(dp), intent(in) :: dist

      !> Resulting truncated gamma
      real(dp) :: gamma

    end function gammaFunc
  end interface


  !> Base type extension for unaltered analytical gamma functions.
  type, extends(THybridXcFunc) :: THybridXcFunc_full

  contains

    procedure :: getLrGammaValue => getLrAnalyticalGammaValue
    procedure :: getLrGammaPrimeValue => getdLrAnalyticalGammaValue

    procedure :: getHfGammaValue => getHfAnalyticalGammaValue
    procedure :: getHfGammaPrimeValue => getdHfAnalyticalGammaValue

  end type THybridXcFunc_full


  !> Base type extension for truncated analytical gamma functions.
  type, extends(THybridXcFunc) :: THybridXcFunc_truncated

  contains

    procedure :: getLrGammaValue => getLrTruncatedGammaValue
    procedure :: getLrGammaPrimeValue => getdLrTruncatedGammaValue

    procedure :: getHfGammaValue => getHfTruncatedGammaValue
    procedure :: getHfGammaPrimeValue => getdHfTruncatedGammaValue

  end type THybridXcFunc_truncated


  !> Base type extension for truncated and poly5zero damped analytical gamma functions.
  type, extends(THybridXcFunc) :: THybridXcFunc_truncatedAndDamped

  contains

    procedure :: getLrGammaValue => getLrTruncatedAndDampedGammaValue
    procedure :: getLrGammaPrimeValue => getdLrTruncatedAndDampedGammaValue

    procedure :: getHfGammaValue => getHfTruncatedAndDampedGammaValue
    procedure :: getHfGammaPrimeValue => getdHfTruncatedAndDampedGammaValue

  end type THybridXcFunc_truncatedAndDamped


  !> Base type extension for minimum image convention, based on density matrix criteria.
  type, extends(THybridXcFunc) :: THybridXcFunc_mic

  contains

    procedure :: getLrGammaValue => getLrMicGammaValue
    procedure :: getLrGammaPrimeValue => getdLrMicGammaValue

    procedure :: getHfGammaValue => getHfMicGammaValue
    procedure :: getHfGammaPrimeValue => getdHfMicGammaValue

  end type THybridXcFunc_mic


contains

  !> Intitializes the range-separated hybrid DFTB module.
  subroutine THybridXcFunc_init(this, nAtom0, species0, hubbu, screeningThreshold, omega, camAlpha,&
      & camBeta, tSpin, tREKS, hybridXcAlg, hybridXcType, gammaType, tPeriodic, tRealHS, errStatus,&
      & gammaCutoff, gSummationCutoff, wignerSeitzReduction, coeffsDiag, latVecs)

    !> Class instance
    class(THybridXcFunc), intent(out), allocatable :: this

    !> Number of atoms in central cell
    integer, intent(in) :: nAtom0

    !> Species-list of atoms in the central cell
    integer, intent(in) :: species0(:)

    !> Species-resolved Hubbard U values
    real(dp), intent(in) :: hubbu(:)

    !> Screening threshold value
    real(dp), intent(in) :: screeningThreshold

    !> Range separation parameter
    real(dp), intent(in) :: omega

    !> CAM alpha parameter
    real(dp), intent(in) :: camAlpha

    !> CAM beta parameter
    real(dp), intent(in) :: camBeta

    !> Is this spin restricted (F) or unrestricted (T)
    logical, intent(in) :: tSpin

    !> Is this DFTB/SSR formalism
    logical, intent(in) :: tREKS

    !> lr-Hamiltonian construction algorithm
    integer, intent(in) :: hybridXcAlg

    !> Hybrid xc-functional type, as extracted from SK-file(s)
    integer, intent(in) :: hybridXcType

    !> Gamma function type (mostly for periodic cases)
    integer, intent(in) :: gammaType

    !> True, if system is periodic (i.e. Gamma-only or k-points)
    logical, intent(in) :: tPeriodic

    !> True, if overlap and Hamiltonian are real-valued
    logical, intent(in) :: tRealHS

    !> Error status
    type(TStatus), intent(inout) :: errStatus

    !> Cutoff for truncated Gamma
    real(dp), intent(in), optional :: gammaCutoff

    !> Cutoff for real-space g-summation
    real(dp), intent(in), optional :: gSummationCutoff

    !> Number of unitcells along each supercell folding direction to substract from MIC Wigner-Seitz
    !! cell construction
    integer, intent(in), optional :: wignerSeitzReduction

    !> Supercell folding coefficients (diagonal elements)
    integer, intent(in), optional :: coeffsDiag(:)

    !> Lattice vectors of (periodic) geometry
    real(dp), intent(in), optional :: latVecs(:,:)

    !! Species indices
    integer :: iSp1, iSp2

    !! Number of unique species in system
    integer :: nUniqueSpecies

    ! Perform basic consistency checks for optional arguments
    if (tPeriodic .and. (.not. present(gSummationCutoff))) then
      @:RAISE_ERROR(errStatus, -1, "HybridXc Module: Periodic systems require g-summation cutoff,&
          & which is not present.")
    end if
    if ((.not. tRealHS) .and. (.not. present(coeffsDiag))) then
      @:RAISE_ERROR(errStatus, -1, "HybridXc Module: General k-point case requires supercell&
          & folding coefficients, which are not present.")
    end if
    if ((.not. tRealHS) .and. gammaType == hybridXcGammaTypes%mic&
        & .and. (.not. present(wignerSeitzReduction))) then
      @:RAISE_ERROR(errStatus, -1, "HybridXc Module: General k-point case with MIC algorithm&
          & requires Wigner-Seitz reduction parameter, which is not present.")
    end if

    ! Allocate selected gamma function types
    select case(gammaType)
    case (hybridXcGammaTypes%full)
      allocate(THybridXcFunc_full:: this)
    case (hybridXcGammaTypes%mic)
      allocate(THybridXcFunc_mic:: this)
    case (hybridXcGammaTypes%truncated)
      if (.not. present(gammaCutoff)) then
        @:RAISE_ERROR(errStatus, -1, "HybridXc Module: Coulomb truncation requires cutoff, which is&
            & not present.")
      end if
      allocate(THybridXcFunc_truncated:: this)
    case (hybridXcGammaTypes%truncatedAndDamped)
      if (.not. present(gammaCutoff)) then
        @:RAISE_ERROR(errStatus, -1, "HybridXc Module: Coulomb truncation requires cutoff, which is&
            & not present.")
      end if
      allocate(THybridXcFunc_truncatedAndDamped:: this)
    case default
      @:RAISE_ERROR(errStatus, -1, "HybridXc Module: Invalid gamma function type obtained.")
    end select

    this%gammaType = gammaType

    this%tScreeningInited = .false.
    this%pScreeningThreshold = screeningThreshold

    this%omega = omega
    this%hybridXcAlg = hybridXcAlg
    this%hybridXcType = hybridXcType
    this%tSpin = tSpin
    this%tREKS = tREKS
    this%camAlpha = camAlpha
    this%camBeta = camBeta
    this%hubbu = hubbu
    this%species0 = species0

    this%camEnergy = 0.0_dp

    if (present(gSummationCutoff)) this%gSummationCutoff = gSummationCutoff
    if (present(wignerSeitzReduction)) this%wignerSeitzReduction = wignerSeitzReduction

    if (present(gammaCutoff)) then
      this%gammaCutoff = gammaCutoff

      ! Start beginning of the damping region and 95% of the gamma cutoff.
      this%gammaDamping = 0.95_dp * this%gammaCutoff

      if (this%gammaDamping <= 0.0_dp) then
        @:RAISE_ERROR(errStatus, -1, "Beginning of damped region of electrostatics must be&
            & positive.")
      end if

      ! Tabulate truncated Gamma properties for all (symmetric) combinations of species
      nUniqueSpecies = getNumberOfUniqueInt(this%species0)
      allocate(this%lrGammaAtDamping(nUniqueSpecies, nUniqueSpecies))
      allocate(this%lrdGammaAtDamping(nUniqueSpecies, nUniqueSpecies))
      allocate(this%lrddGammaAtDamping(nUniqueSpecies, nUniqueSpecies))
      allocate(this%hfGammaAtDamping(nUniqueSpecies, nUniqueSpecies))
      allocate(this%hfdGammaAtDamping(nUniqueSpecies, nUniqueSpecies))
      allocate(this%hfddGammaAtDamping(nUniqueSpecies, nUniqueSpecies))
      do iSp2 = 1, nUniqueSpecies
        do iSp1 = 1, nUniqueSpecies
          this%lrGammaAtDamping(iSp1, iSp2) = getLrAnalyticalGammaValue_workhorse(this%hubbu(iSp1),&
              & this%hubbu(iSp2), this%omega, this%gammaDamping)
          this%lrdGammaAtDamping(iSp1, iSp2)&
              & = getdLrAnalyticalGammaValue_workhorse(this%hubbu(iSp1), this%hubbu(iSp2),&
              & this%omega, this%gammaDamping)
          this%lrddGammaAtDamping(iSp1, iSp2)&
              & = getddLrNumericalGammaValue_workhorse(this%hubbu(iSp1), this%hubbu(iSp2),&
              & this%omega, this%gammaDamping, 1.0E-8_dp)
          this%hfGammaAtDamping(iSp1, iSp2) = getHfAnalyticalGammaValue_workhorse(this%hubbu(iSp1),&
              & this%hubbu(iSp2), this%gammaDamping)
          this%hfdGammaAtDamping(iSp1, iSp2)&
              & = getdHfAnalyticalGammaValue_workhorse(this%hubbu(iSp1), this%hubbu(iSp2),&
              & this%gammaDamping)
          this%hfddGammaAtDamping(iSp1, iSp2) = getddHfNumericalGammaDeriv(this, iSp1, iSp2,&
              & this%gammaDamping, 1.0E-8_dp)
        end do
      end do
    end if

    if (this%tREKS .and. this%hybridXcType == hybridXcFunc%hyb) then
      @:RAISE_ERROR(errStatus, -1, "Global hybrid functionals not currently implemented for REKS.")
    end if

    if (this%tREKS .and. this%hybridXcType == hybridXcFunc%cam) then
      @:RAISE_ERROR(errStatus, -1, "General CAM functionals not currently implemented for REKS.")
    end if

    allocate(this%coords(3, nAtom0))
    this%coords(:,:) = 0.0_dp
    allocate(this%rCoords(3, nAtom0))
    this%rCoords(:,:) = 0.0_dp

    allocate(this%camGammaEval0(nAtom0, nAtom0))
    this%camGammaEval0(:,:) = 0.0_dp

    ! Check for current restrictions
    if (this%tSpin .and. this%hybridXcAlg == hybridXcAlgo%thresholdBased) then
      @:RAISE_ERROR(errStatus, -1, "Spin-unrestricted calculation for thresholded range-separation&
          & not yet implemented!")
    end if

    if (this%tREKS .and. this%hybridXcAlg == hybridXcAlgo%thresholdBased) then
      @:RAISE_ERROR(errStatus, -1, "REKS calculation with thresholded range-separation not yet&
          & implemented!")
    end if

    if (.not. any([hybridXcAlgo%neighbourBased, hybridXcAlgo%thresholdBased,&
          & hybridXcAlgo%matrixBased] == this%hybridXcAlg)) then
      @:RAISE_ERROR(errStatus, -1, "Unknown algorithm for screening the exchange in&
          & range-separation!")
    end if

    if (tPeriodic .and. (.not. present(latVecs))) then
      @:RAISE_ERROR(errStatus, -1, "HybridXc Module: Periodic structure, but no lattice vectors&
          & handed over.")
    end if

    if (present(coeffsDiag)) then
      this%coeffsDiag = coeffsDiag
      call getBvKLatticeShifts(this%coeffsDiag, this%bvKShifts)
    end if


  contains

    !> Returns the number of unique integers of an array.
    pure function getNumberOfUniqueInt(array) result(nUnique)

      !> Array to investigate
      integer, intent(in) :: array(:)

      !> Number of unique entries
      integer :: nUnique

      !! Auxiliary variables
      integer :: ii, jj, tmp(size(array))

      tmp(:) = array
      call heap_sort(tmp)

      nUnique = 1
      jj = tmp(1)
      do ii = 2, size(tmp)
        if (tmp(ii) /= jj) then
          jj = tmp(ii)
          nUnique = nUnique + 1
        end if
      end do

    end function getNumberOfUniqueInt


    !> Returns BvK real-space shifts, compatible with k-point mesh.
    subroutine getBvKLatticeShifts(coeffsDiag, bvKShifts)

      !> Supercell folding coefficients (diagonal elements)
      integer, intent(in) :: coeffsDiag(:)

      !> The k-point compatible BvK real-space shifts in relative coordinates
      real(dp), intent(out), allocatable :: bvKShifts(:,:)

      !! Number of BvK real-space shifts
      integer :: nBvKShifts

      !! Auxiliary variables
      integer :: ii, jj, kk, ind

      nBvKShifts = product(coeffsDiag)
      allocate(bvKShifts(3, nBvKShifts))

      ind = 1
      do kk = 0, coeffsDiag(3) - 1
        do jj = 0, coeffsDiag(2) - 1
          do ii = 0, coeffsDiag(1) - 1
            bvKShifts(1, ind) = real(ii, dp)
            bvKShifts(2, ind) = real(jj, dp)
            bvKShifts(3, ind) = real(kk, dp)
            ind = ind + 1
          end do
        end do
      end do

    end subroutine getBvKLatticeShifts

  end subroutine THybridXcFunc_init


  !> Checks if obtained supercell folding matrix meets current requirements.
  subroutine checkSupercellFoldingMatrix(supercellFoldingMatrix, errStatus, supercellFoldingDiagOut)

    !> Coefficients of the lattice vectors in the linear combination for the super lattice vectors
    !! (should be integer values) and shift of the grid along the three small reciprocal lattice
    !! vectors (between 0.0 and 1.0)
    real(dp), intent(in), target :: supercellFoldingMatrix(:,:)

    !> Error status
    type(TStatus), intent(inout) :: errStatus

    !> Diagonal elements of supercell folding matrix, if present
    integer, intent(out), optional :: supercellFoldingDiagOut(:)

    !! Supercell folding coefficients and shifts
    real(dp), pointer :: coeffs(:,:), shifts(:)

    !! True, if the supercell folding does not correspond to a MP-like scheme
    logical :: tNotMonkhorstPack

    !! Auxiliary variables
    integer :: ii, jj

    if (present(supercellFoldingDiagOut)) then
      @:ASSERT(size(supercellFoldingDiagOut) == 3)
    end if

    coeffs => supercellFoldingMatrix(:, 1:3)
    shifts => supercellFoldingMatrix(:, 4)

    if (abs(determinant33(coeffs)) - 1.0_dp < -1e-06_dp) then
      @:RAISE_ERROR(errStatus, -1, "Determinant of the supercell matrix must be greater than 1.")
    end if

    if (any(abs(modulo(coeffs + 0.5_dp, 1.0_dp) - 0.5_dp) > 1e-6_dp)) then
      @:RAISE_ERROR(errStatus, -1, "The components of the supercell matrix must be integers.")
    end if

    ! Check if k-point mesh is a Monkhorst-Pack sampling with zero shift
    tNotMonkhorstPack = .false.
    lpOuter: do jj = 1, size(coeffs, dim=2)
      do ii = 1, size(coeffs, dim=1)
        if (ii == jj) cycle
        if (coeffs(ii, jj) > 1e-06_dp) then
          tNotMonkhorstPack = .true.
          exit lpOuter
        end if
      end do
    end do lpOuter
    if (tNotMonkhorstPack) then
      @:RAISE_ERROR(errStatus, -1, "Range-separated calculations with k-points require a&
          & Monkhorst-Pack-like sampling, i.e. a uniform extension of the lattice.")
    end if

    ! Check if shifts are zero
    if (any(abs(shifts) > 1e-06_dp)) then
      @:RAISE_ERROR(errStatus, -1, "Range-separated calculations with k-points require a&
          & Monkhorst-Pack-like sampling with zero shift.")
    end if

    ! All checks have passed, continue...

    ! Get diagonal elements as integers, if requested
    if (present(supercellFoldingDiagOut)) then
      do ii = 1, 3
        supercellFoldingDiagOut(ii) = nint(coeffs(ii, ii))
      end do
    end if

  end subroutine checkSupercellFoldingMatrix


  !> Folds relative real-space vector back to BvK region.
  pure function THybridXcFunc_foldToBvK(this, vector) result(bvKShift)

    !> Class instance
    class(THybridXcFunc), intent(in) :: this

    !> Vector (in relative coordinates) to fold back to BvK cell
    real(dp), intent(in) :: vector(:)

    !> Corresponding BvK vector
    integer :: bvKShift(3)

    bvKShift(:) = modulo(nint(vector), this%coeffsDiag)

  end function THybridXcFunc_foldToBvK


  !> Folds relative real-space vector back to BvK region and returns indices of density matrix.
  pure function THybridXcFunc_foldToBvKIndex(this, vector) result(bvKIndex)

    !> Class instance
    class(THybridXcFunc), intent(in) :: this

    !> Vector (in relative coordinates) to fold back to BvK cell
    real(dp), intent(in) :: vector(:)

    !> Corresponding BvK indexing
    integer :: bvKIndex(3)

    ! additionally shift by 1, so that indices start at 1 and not at 0
    bvKIndex(:) = modulo(nint(vector), this%coeffsDiag) + 1

  end function THybridXcFunc_foldToBvKIndex


  !> Updates the range-separated module on coordinate change (non-periodic version).
  subroutine updateCoords_cluster(this, rCoords)

    !> Class instance
    class(THybridXcFunc), intent(inout) :: this

    !> Atomic coordinates in absolute units
    real(dp), intent(in) :: rCoords(:,:)

    !! Indices of interacting atoms in central cell, as well as their global species index
    integer :: iAtom1, iAtom2, iSp1, iSp2

    !! Number of atoms in central cell
    integer :: nAtom0

    !! Distance between two interacting atoms
    real(dp) :: dist

    this%rCoords(:,:) = rCoords
    nAtom0 = size(this%species0)

    do iAtom1 = 1, nAtom0
      do iAtom2 = 1, iAtom1
        iSp1 = this%species0(iAtom1)
        iSp2 = this%species0(iAtom2)
        dist = norm2(this%rCoords(:, iAtom1) - this%rCoords(:, iAtom2))
        this%camGammaEval0(iAtom1, iAtom2) = getCamAnalyticalGammaValue_workhorse(this%hubbu(iSp1),&
            & this%hubbu(iSp2), this%omega, this%camAlpha, this%camBeta, dist)
        this%camGammaEval0(iAtom2, iAtom1) = this%camGammaEval0(iAtom1, iAtom2)
      end do
    end do

    if (this%tScreeningInited) then
      this%hPrev(:,:) = 0.0_dp
      this%dRhoPrev(:,:) = 0.0_dp
      this%camEnergy = 0.0_dp
    end if

  end subroutine updateCoords_cluster


  !> Updates the range-separated module on coordinate change (Gamma-only version).
  subroutine updateCoords_gamma(this, env, symNeighbourList, nNeighbourCamSym, skOverCont, orb,&
      & latVecs, recVecs2p, iSquare)

    !> Class instance
    class(THybridXcFunc), intent(inout), target :: this

    !> Environment settings
    type(TEnvironment), intent(in) :: env

    !> List of neighbours for each atom (symmetric version)
    type(TSymNeighbourList), intent(in) :: symNeighbourList

    !> Symmetric neighbour list version of nNeighbourCam
    integer, intent(in) :: nNeighbourCamSym(:)

    !> Sparse overlap part
    type(TSlakoCont), intent(in) :: skOverCont

    !> Orbital information
    type(TOrbitals), intent(in) :: orb

    !> Lattice vectors of (periodic) geometry
    real(dp), intent(in) :: latVecs(:,:)

    !> Reciprocal lattice vectors in units of 2pi
    real(dp), intent(in) :: recVecs2p(:,:)

    !> Position of each atom in the rows/columns of the square matrices. Shape: (nAtom)
    integer, intent(in) :: iSquare(:)

    !! Indices of interacting atoms in central cell, as well as their global species index
    integer :: iAtM, iAtN, iSpM, iSpN

    !! Number of atoms in central cell
    integer :: nAtom0

    ! range-separated type is initialized with nAtom0 coordinates, therefore re-allocate for
    ! periodic systems, where images beyond the central cell are accounted for
    if (allocated(this%coords)) deallocate(this%coords)
    if (allocated(this%rCoords)) deallocate(this%rCoords)
    this%rCoords = symNeighbourList%coord
    allocate(this%coords(size(this%rCoords, dim=1), size(this%rCoords, dim=2)))

    ! calculate neighbour list coordinates in relative units
    this%coords(:,:) = this%rCoords
    call cart2frac(this%coords, latVecs)

    ! re-allocate sparse overlap for symmetric neighbour list
    if (allocated(this%overSym)) deallocate(this%overSym)
    allocate(this%overSym(symNeighbourList%sparseSize))

    ! get number of atoms in central cell
    nAtom0 = size(this%species0)

    ! get all cell translations within given cutoff
    call getCellTranslations(this%cellVecsG, this%rCellVecsG, latVecs, recVecs2p,&
        & this%gSummationCutoff)

    ! build symmetric, sparse overlap
    call buildS(env, this%overSym, skOverCont, this%rCoords, nNeighbourCamSym,&
        & symNeighbourList%neighbourList%iNeighbour, symNeighbourList%species,&
        & symNeighbourList%iPair, orb)

    if (this%tScreeningInited) then
      this%hPrev(:,:) = 0.0_dp
      this%dRhoPrev(:,:) = 0.0_dp
      this%camEnergy = 0.0_dp
    end if

    ! pre-tabulate overlap estimates for neighbour-list based algorithms
    if (this%hybridXcAlg == hybridXcAlgo%neighbourBased) then
      call calculateOverlapEstimates(this, symNeighbourList, nNeighbourCamSym, iSquare)
    end if

    ! \sum\gamma(\bm{g}) pre-tabulation
    do iAtM = 1, nAtom0
      iSpM = this%species0(iAtM)
      do iAtN = 1, nAtom0
        iSpN = this%species0(iAtN)
        this%camGammaEval0(iAtM, iAtN) = getCamGammaGSum(this, iAtM, iAtN, iSpM, iSpN,&
            & this%rCellVecsG)
        this%camGammaEval0(iAtN, iAtM) = getCamGammaGSum(this, iAtN, iAtM, iSpN, iSpM,&
            & this%rCellVecsG)
      end do
    end do

  end subroutine updateCoords_gamma


  !> Updates the range-separated module on coordinate change (k-point version).
  subroutine updateCoords_kpts(this, env, symNeighbourList, nNeighbourCamSym, skOverCont, orb,&
      & latVecs, recVecs2p, iSquare)

    !> Class instance
    class(THybridXcFunc), intent(inout), target :: this

    !> Environment settings
    type(TEnvironment), intent(in) :: env

    !> List of neighbours for each atom (symmetric version)
    type(TSymNeighbourList), intent(in) :: symNeighbourList

    !> Symmetric neighbour list version of nNeighbourCam
    integer, intent(in) :: nNeighbourCamSym(:)

    !> Sparse overlap part
    type(TSlakoCont), intent(in) :: skOverCont

    !> Orbital information
    type(TOrbitals), intent(in) :: orb

    !> Lattice vectors of (periodic) geometry
    real(dp), intent(in) :: latVecs(:,:)

    !> Reciprocal lattice vectors in units of 2pi
    real(dp), intent(in) :: recVecs2p(:,:)

    !> Position of each atom in the rows/columns of the square matrices. Shape: (nAtom)
    integer, intent(in) :: iSquare(:)

    !! Dense matrix descriptor indices
    integer, parameter :: descLen = 3, iStart = 1, iEnd = 2, iNOrb = 3

    !! Index array for descending sorting of Gamma arrays
    integer, allocatable :: gammasortIdx(:)

    !! Temporary storage for g-resolved gamma values (+ directional derivatives)
    real(dp), allocatable :: gammaEvalGTmp(:), dGammaEvalGTmp(:,:)

    !! Number of non-zero entries of (sorted) array
    integer :: nNonZeroEntries

    !! Indices of interacting atoms in central cell, as well as their global species index
    integer :: iAtM, iAtN, iSpM, iSpN

    !! Number of atoms in the central cell
    integer :: nAtom0

    !! Number of g-summation summands
    integer :: nGShifts

    !! Composite index iAtM/iAtN
    integer :: ii
    integer, allocatable :: iAtMN(:,:)

    !! Dummy array with zeros
    real(dp) :: zeros(3)

    !! Iterates over g-vectors
    integer :: iG

    zeros(:) = 0.0_dp

    nAtom0 = size(this%species0)

    ! range-separated type is initialized with nAtom0 coordinates, therefore re-allocate for
    ! periodic systems, where images beyond the central cell are accounted for
    if (allocated(this%coords)) deallocate(this%coords)
    if (allocated(this%rCoords)) deallocate(this%rCoords)
    this%rCoords = symNeighbourList%coord
    allocate(this%coords(size(this%rCoords, dim=1), size(this%rCoords, dim=2)))

    ! calculate neighbour list coordinates in relative units
    this%coords(:,:) = this%rCoords
    call cart2frac(this%coords, latVecs)

    ! re-allocate sparse overlap for symmetric neighbour list
    if (allocated(this%overSym)) deallocate(this%overSym)
    allocate(this%overSym(symNeighbourList%sparseSize))

    ! build symmetric, sparse overlap
    call buildS(env, this%overSym, skOverCont, this%rCoords, nNeighbourCamSym,&
        & symNeighbourList%neighbourList%iNeighbour, symNeighbourList%species,&
        & symNeighbourList%iPair, orb)

    if (this%tScreeningInited) then
      this%hPrevCplxHS(:,:,:) = 0.0_dp
      ! Note, this is only allocated if passing through addCamHamiltonianNeighbour_kpts_mic or
      ! addCamHamiltonianNeighbour_kpts_ct :
      if (allocated(this%dRhoPrevCplxHS)) this%dRhoPrevCplxHS(:,:,:,:,:,:) = 0.0_dp
      this%camEnergy = 0.0_dp
    end if

    ! pre-tabulate overlap estimates for neighbour-list based algorithms
    ! if-branch not really necessary, since k-implementation only available for neighbour-list
    ! based algorithm anyway, but let's be on the save side
    if (this%hybridXcAlg == hybridXcAlgo%neighbourBased) then
      call calculateOverlapEstimates(this, symNeighbourList, nNeighbourCamSym, iSquare)
    end if

    ! beginning of \gamma(\bm{g}) pre-tabulation

    if (this%gammaType == hybridXcGammaTypes%mic) then
      if (allocated(this%wsVectors)) deallocate(this%wsVectors)
      ! Generate "save" Wigner-Seitz vectors for density matrix arguments
      call generateWignerSeitzGrid(max(this%coeffsDiag - this%wignerSeitzReduction, 1), latVecs,&
          & this%wsVectors)

      ! The Wigner-Seitz grid actually defines all the relevant g-vectors, therefore copy over
      this%cellVecsG = this%wsVectors
      if (allocated(this%rCellVecsG)) deallocate(this%rCellVecsG)
      allocate(this%rCellVecsG(3, size(this%cellVecsG, dim=2)))
      do iG = 1, size(this%rCellVecsG, dim=2)
        this%rCellVecsG(:, iG) = matmul(latVecs, this%cellVecsG(:, iG))
      end do

    else

      call getTwoLoopCompositeIndex(nAtom0, nAtom0, iAtMN)

      ! get all cell translations within given cutoff
      call getCellTranslations(this%cellVecsG, this%rCellVecsG, latVecs, recVecs2p,&
          & this%gSummationCutoff)

      if (allocated(this%nNonZeroGammaG)) deallocate(this%nNonZeroGammaG)
      allocate(this%nNonZeroGammaG(nAtom0, nAtom0))

      if (allocated(this%camGammaEvalG)) deallocate(this%camGammaEvalG)
      allocate(this%camGammaEvalG(nAtom0, nAtom0))
      if (allocated(this%camdGammaEvalG)) deallocate(this%camdGammaEvalG)
      allocate(this%camdGammaEvalG(nAtom0, nAtom0))

      nGShifts = size(this%cellVecsG, dim=2)
      allocate(gammaEvalGTmp(nGShifts))
      allocate(dGammaEvalGTmp(3, nGShifts))
      allocate(gammaSortIdx(nGShifts))

      ! Pre-tabulate CAM gamma integrals (+ directional derivatives)
      do ii = 1, nAtom0**2
        iAtM = iAtMN(1, ii)
        iAtN = iAtMN(2, ii)
        iSpM = this%species0(iAtM)
        iSpN = this%species0(iAtN)
        ! pre-tabulate g-resolved \gamma_{\mu\nu}(\vec{g})
        gammaEvalGTmp(:) = getCamGammaGResolved(this, iAtM, iAtN, iSpM, iSpN, this%rCellVecsG,&
            & zeros)
        call index_heap_sort(gammaSortIdx, gammaEvalGTmp)
        gammaSortIdx(:) = gammaSortIdx(size(gammaSortIdx):1:-1)
        gammaEvalGTmp(:) = gammaEvalGTmp(gammaSortIdx)
        nNonZeroEntries = getNumberOfNonZeroElements(gammaEvalGTmp)
        this%nNonZeroGammaG(iAtM, iAtN)%data = gammaSortIdx(1:nNonZeroEntries)
        this%camGammaEvalG(iAtM, iAtN)%data = gammaEvalGTmp(1:nNonZeroEntries)
        ! pre-tabulate g-resolved \partial\gamma_{\mu\nu}(\vec{g})
        dGammaEvalGTmp(:,:) = getCamGammaPrimeGResolved(this, iAtM, iAtN, iSpM, iSpN,&
            & this%rCellVecsG(:, gammaSortIdx))
        this%camdGammaEvalG(iAtM, iAtN)%data = dGammaEvalGTmp(:, 1:nNonZeroEntries)
      end do

    end if


  contains

    !> Returns the number of non-zero elements in a descending array of non-negative reals.
    function getNumberOfNonZeroElements(array) result(nNonZeroEntries)

      !> Descending, one-dimensional, real-valued array to search
      real(dp), intent(in) :: array(:)

      !> Number of non-zero entries
      integer :: nNonZeroEntries

      !! iterates over all array elements
      integer :: ii

      nNonZeroEntries = 0

      do ii = 1, size(array)
        if (array(ii) < epsilon(1.0_dp)) then
          @:ASSERT(all(array(ii:) >= 0.0_dp))
          return
        end if
        nNonZeroEntries = ii
      end do

    end function getNumberOfNonZeroElements

  end subroutine updateCoords_kpts


  !> Tabulates (descending) overlap estimates for SPS-product screening.
  subroutine calculateOverlapEstimates(this, symNeighbourList, nNeighbourCamSym, iSquare)

    !> Class instance
    class(THybridXcFunc), intent(inout), target :: this

    !> List of neighbours for each atom (symmetric version)
    type(TSymNeighbourList), intent(in) :: symNeighbourList

    !> Symmetric neighbour list version of nNeighbourCam
    integer, intent(in) :: nNeighbourCamSym(:)

    !> Position of each atom in the rows/columns of the square matrices. Shape: (nAtom)
    integer, intent(in) :: iSquare(:)

    !! Diatomic blocks from sparse, real-space overlap matrix S_{\beta\nu}
    real(dp), pointer :: pSbn(:,:)

    !! Dense matrix descriptor indices
    integer, parameter :: descLen = 3, iStart = 1, iEnd = 2, iNOrb = 3

    !! Index of atom in central cell
    integer :: iAtN

    !! Neighbour index (+corresponding atom index)
    integer :: iNeighN, iAtB

    !! Folded (to central cell) atom index
    integer :: iAtBfold

    !! Stores start/end index and number of orbitals of square matrices
    integer :: descN(descLen), descB(descLen)

    !! Number of atoms in central cell
    integer :: nAtom0

    !! Auxiliary variables for setting up 2D pointer to sparse overlap
    integer :: ind, nOrbAt, nOrbNeigh

    ! get number of atoms in central cell
    nAtom0 = size(this%species0)

    ! allocate max estimates of square overlap blocks and index array for sorting
    if (allocated(this%squareOverEst)) deallocate(this%squareOverEst)
    allocate(this%squareOverEst(nAtom0))
    if (allocated(this%overlapIndices)) deallocate(this%overlapIndices)
    allocate(this%overlapIndices(nAtom0))

    do iAtN = 1, nAtom0
      descN = getDescriptor(iAtN, iSquare)
      allocate(this%squareOverEst(iAtN)%data(nNeighbourCamSym(iAtN) + 1))
      do iNeighN = 0, nNeighbourCamSym(iAtN)
        iAtB = symNeighbourList%neighbourList%iNeighbour(iNeighN, iAtN)
        iAtBfold = symNeighbourList%img2CentCell(iAtB)
        descB = getDescriptor(iAtBfold, iSquare)
        ! get 2D pointer to Sbn overlap block
        ind = symNeighbourList%iPair(iNeighN, iAtN) + 1
        nOrbAt = descN(iNOrb)
        nOrbNeigh = descB(iNOrb)
        pSbn(1:nOrbNeigh, 1:nOrbAt) => this%overSym(ind:ind + nOrbNeigh * nOrbAt - 1)
        this%squareOverEst(iAtN)%data(iNeighN + 1) = maxval(abs(pSbn))
      end do
    end do

    ! sort max square overlap estimates (descending)
    ! this way we can exit the whole loop a.s.a. threshold has been undershot for the first time
    do iAtN = 1, nAtom0
      allocate(this%overlapIndices(iAtN)%data(nNeighbourCamSym(iAtN) + 1))
      call index_heap_sort(this%overlapIndices(iAtN)%data, this%squareOverEst(iAtN)%data)
      ! switch from ascending to descending
      this%overlapIndices(iAtN)%data(:)&
          & = this%overlapIndices(iAtN)%data(size(this%overlapIndices(iAtN)%data):1:-1)
    end do

  end subroutine calculateOverlapEstimates


  !> Builds composite index for two nested loops (e.g. over atoms in the central cell).
  subroutine getTwoLoopCompositeIndex(nInner, nOuter, iComposite)

    !> Number of iterations in inner and outer loop
    integer, intent(in) :: nInner, nOuter

    !> Composite index for two nested loops
    integer, intent(out), allocatable :: iComposite(:,:)

    !! Auxiliary variables
    integer :: ind, ii, jj

    allocate(iComposite(2, nInner * nOuter))

    ! Build up composite index for collapsing two nested loops
    ind = 1
    lpOuter: do ii = 1, nOuter
      lpInner: do jj = 1, nInner
        iComposite(1, ind) = ii
        iComposite(2, ind) = jj
        ind = ind + 1
      end do lpInner
    end do lpOuter

  end subroutine getTwoLoopCompositeIndex


  !> Builds simple composite index for two nested loops over atoms in the central cell.
  subroutine getiKSiKSPrimeCompositeIndex(env, nS, nK, nKPrime, iKSComposite)

    !> Environment settings
    type(TEnvironment), intent(in) :: env

    !> Total number of spin channels
    integer, intent(in) :: nS

    !> Total number of k-points
    integer, intent(in) :: nK, nKPrime

    !> Composite index for two nested loops over central cell
    integer, intent(out), allocatable :: iKSComposite(:,:)

    !! Indices of spins/k-points
    integer :: iS, iK, iKPrime

    !! Loop counters for global and local composite index
    integer :: indGlobal, indLocal

    !! Start and end index for MPI parallelization, if applicable
    integer :: iParallelStart, iParallelEnd

    call getStartAndEndIndex(env, nS * nK * nKPrime, iParallelStart, iParallelEnd)

    ! Composite index needs to be allocated on all ranks
    allocate(iKSComposite(3, iParallelEnd - iParallelStart + 1))

    ! Build up composite index iKSComposite for collapsing iKS-iKSPrime summations
    indGlobal = 1
    indLocal = 1
    do iS = 1, nS
      do iK = 1, nK
        do iKPrime = 1, nKPrime
          if (indGlobal >= iParallelStart .and. indGlobal <= iParallelEnd) then
            iKSComposite(1, indLocal) = iS
            iKSComposite(2, indLocal) = iK
            iKSComposite(3, indLocal) = iKPrime
            indLocal = indLocal + 1
          end if
          indGlobal = indGlobal + 1
        end do
      end do
    end do

  end subroutine getiKSiKSPrimeCompositeIndex


  !> Builds MPI rank dependent composite index for nested HFX loops.
  subroutine getFourLoopCompositeIndex(this, env, nNeighbourCamSym, pMax, compositeIndex, errStatus)

    !> Class instance
    class(THybridXcFunc), intent(in) :: this

    !> Environment settings
    type(TEnvironment), intent(in) :: env

    !> Nr. of neighbours for each atom.
    integer, intent(in) :: nNeighbourCamSym(:)

    !! Max estimate for delta-delta density matrix
    real(dp), intent(in) :: pMax

    !> Composite index for four nested loops
    integer, intent(out), allocatable :: compositeIndex(:,:)

    !> Error status
    type(TStatus), intent(inout) :: errStatus

    !! Loop counters for global and local composite index
    integer :: indGlobal, indLocal

    !! Atom indices (central cell)
    integer :: iAtM, iAtN

    !! Neighbour indices
    integer :: iNeighN, iNeighM

    !! Sorted (according to max overlap estimates) neighbour indices
    integer :: iNeighMsort, iNeighNsort

    !! Max. estimate for products of delta-delta density matrix and max overlap estimates
    real(dp) :: maxEstimate, pSbnPabMax

    !! Start and end index for MPI parallelization, if applicable
    integer :: iParallelStart, iParallelEnd

    !! Loop counter
    integer :: ind

    !! Number of atoms in central cell
    integer :: nAtom0

    if (.not. (allocated(this%squareOverEst) .and. allocated(this%overlapIndices))) then
      @:RAISE_ERROR(errStatus, -1, "HybridXc Module: 4-loop composite index construction depends on&
          & overlap estimates, that are not allocated yet.")
    end if

    nAtom0 = size(this%species0)

    ! Lead process counts composite index entries for MPI parallelization
    ! (includes integral screening)
    if (env%tGlobalLead) then
      ! First run over loops to count non-vanishing cycles
      ind = 0
      loopM1: do iAtM = 1, nAtom0
        loopN1: do iAtN = 1, nAtom0
          loopB1: do iNeighN = 0, nNeighbourCamSym(iAtN)
            iNeighNsort = this%overlapIndices(iAtN)%data(iNeighN + 1) - 1
            pSbnPabMax = pMax * this%squareOverEst(iAtN)%data(iNeighNsort + 1)
            if (pSbnPabMax < this%pScreeningThreshold) exit loopB1
            loopA1: do iNeighM = 0, nNeighbourCamSym(iAtM)
              iNeighMsort = this%overlapIndices(iAtM)%data(iNeighM + 1) - 1
              maxEstimate = pSbnPabMax * this%squareOverEst(iAtM)%data(iNeighMsort + 1)
              if (maxEstimate < this%pScreeningThreshold) exit loopA1

              ind = ind + 1

            end do loopA1
          end do loopB1
        end do loopN1
      end do loopM1
    end if

  #:if WITH_MPI
    call mpifx_bcast(env%mpi%globalComm, ind)
  #:endif

    call getStartAndEndIndex(env, ind, iParallelStart, iParallelEnd)

    ! Composite index needs to be allocated on all ranks
    allocate(compositeIndex(4, iParallelEnd - iParallelStart + 1))

    ! Flatten four nested loops into a single composite index for MPI parallelization
    indGlobal = 1
    indLocal = 1
    loopM2: do iAtM = 1, nAtom0
      loopN2: do iAtN = 1, nAtom0
        loopB2: do iNeighN = 0, nNeighbourCamSym(iAtN)
          iNeighNsort = this%overlapIndices(iAtN)%data(iNeighN + 1) - 1
          pSbnPabMax = pMax * this%squareOverEst(iAtN)%data(iNeighNsort + 1)
          if (pSbnPabMax < this%pScreeningThreshold) exit loopB2
          loopA2: do iNeighM = 0, nNeighbourCamSym(iAtM)
            iNeighMsort = this%overlapIndices(iAtM)%data(iNeighM + 1) - 1
            maxEstimate = pSbnPabMax * this%squareOverEst(iAtM)%data(iNeighMsort + 1)
            if (maxEstimate < this%pScreeningThreshold) exit loopA2
            if (indGlobal >= iParallelStart .and. indGlobal <= iParallelEnd) then
              compositeIndex(1, indLocal) = iAtM
              compositeIndex(2, indLocal) = iAtN
              compositeIndex(3, indLocal) = iNeighNsort
              compositeIndex(4, indLocal) = iNeighMsort
              indLocal = indLocal + 1
            end if
            indGlobal = indGlobal + 1
          end do loopA2
        end do loopB2
      end do loopN2
    end do loopM2

  end subroutine getFourLoopCompositeIndex


#:if WITH_SCALAPACK

  !> Interface routine for adding CAM range-separated contributions to the Hamiltonian.
  !! (non-periodic and Gamma-only version)
  subroutine addCamHamiltonian_real(this, env, denseDesc, SSqrReal, deltaRhoSqr, HSqrReal,&
      & errStatus)

    !> Class instance
    class(THybridXcFunc), intent(inout) :: this

    !> Environment settings
    type(TEnvironment), intent(inout) :: env

    !> Dense matrix descriptor
    type(TDenseDescr), intent(in) :: denseDesc

    !> Square real overlap matrix
    real(dp), intent(in) :: SSqrReal(:,:)

    !> Square (unpacked) density matrix
    real(dp), intent(in) :: deltaRhoSqr(:,:)

    !> Square (unpacked) Hamiltonian to be updated
    real(dp), intent(inout) :: HSqrReal(:,:)

    !> Error status
    type(TStatus), intent(inout) :: errStatus

    call env%globalTimer%startTimer(globalTimers%hybridXcH)

    select case(this%hybridXcAlg)
    case (hybridXcAlgo%thresholdBased)
      @:RAISE_ERROR(errStatus, -1, "HybridXc Module: Thresholded algorithm for non-periodic or&
          & Gamma-only systems does not yet support MPI parallelism (choose matrix-based algorithm&
          & instead).")
    case (hybridXcAlgo%neighbourBased)
      @:RAISE_ERROR(errStatus, -1, "HybridXc Module: Neighbour-list based algorithm for&
          & non-periodic or Gamma-only systems does not yet support MPI parallelism (choose&
          & matrix-based algorithm instead).")
    case (hybridXcAlgo%matrixBased)
      call addCamHamiltonianMatrix_real_blacs(this, env, denseDesc, SSqrReal, deltaRhoSqr, HSqrReal)
    end select

    call env%globalTimer%stopTimer(globalTimers%hybridXcH)

  end subroutine addCamHamiltonian_real

#:else

  !> Interface routine for adding CAM range-separated contributions to the Hamiltonian.
  !! (non-periodic and Gamma-only version)
  subroutine addCamHamiltonian_real(this, env, deltaRhoSqr, SSqrReal, overSparse, iNeighbour,&
      & nNeighbourCam, iSquare, iPair, orb, img2CentCell, tPeriodic, HSqrReal, errStatus)

    !> Class instance
    class(THybridXcFunc), intent(inout) :: this

    !> Environment settings
    type(TEnvironment), intent(inout) :: env

    !> Square (unpacked) density matrix
    real(dp), intent(in) :: deltaRhoSqr(:,:)

    !> Square real overlap matrix
    real(dp), intent(in) :: SSqrReal(:,:)

    !> Sparse (packed) overlap matrix
    real(dp), intent(in) :: overSparse(:)

    !> Neighbour indices
    integer, intent(in) :: iNeighbour(0:,:)

    !> Nr. of neighbours for each atom
    integer, intent(in) :: nNeighbourCam(:)

    !> Position of each atom in the rows/columns of the square matrices. Shape: (nAtom)
    integer, intent(in) :: iSquare(:)

    !> iPair Position of each (neighbour, atom) pair in the sparse matrix
    !> Shape: (0:maxNeighbour, nAtom)
    integer, intent(in) :: iPair(0:,:)

    !> Orbital information
    type(TOrbitals), intent(in) :: orb

    !> Map images of atoms to the central cell
    integer, intent(in) :: img2CentCell(:)

    !> True, if system is periodic (i.e. Gamma-only)
    logical, intent(in) :: tPeriodic

    !> Square (unpacked) Hamiltonian to be updated
    real(dp), intent(inout) :: HSqrReal(:,:)

    !> Error status
    type(TStatus), intent(inout) :: errStatus

    call env%globalTimer%startTimer(globalTimers%hybridXcH)

    select case(this%hybridXcAlg)
    case (hybridXcAlgo%thresholdBased)
      if (tPeriodic) then
        @:RAISE_ERROR(errStatus, -1, "HybridXc Module: Thresholded algorithm not implemented for&
            & Gamma-point periodic systems.")
      else
        call addCamHamiltonianThreshold_cluster(this, SSqrReal, deltaRhoSqr, iNeighbour,&
            & nNeighbourCam, iSquare, HSqrReal, orb)
      end if
    case (hybridXcAlgo%neighbourBased)
      call addCamHamiltonianNeighbour_real(this, deltaRhoSqr, overSparse, iNeighbour,&
          & nNeighbourCam, iSquare, iPair, orb, img2CentCell, HSqrReal)
    case (hybridXcAlgo%matrixBased)
      call addCamHamiltonianMatrix_real(this, iSquare, SSqrReal, deltaRhoSqr, HSqrReal)
    end select

    call env%globalTimer%stopTimer(globalTimers%hybridXcH)

  end subroutine addCamHamiltonian_real

#:endif


  !> Interface routine for adding CAM range-separated contributions to the Hamiltonian.
  !! (k-point version)
  subroutine getCamHamiltonian_kpts(this, env, denseDesc, orb, ints, densityMatrix, neighbourList,&
      & nNeighbourSK, symNeighbourList, nNeighbourCamSym, iCellVec, cellVecs, rCellVecs,&
      & iSparseStart, img2CentCell, kPoints, kWeights, HSqrCplxCam, errStatus)

    !> Class instance
    class(THybridXcFunc), intent(inout) :: this

    !> Environment settings
    type(TEnvironment), intent(inout) :: env

    !> Dense matrix descriptor
    type(TDenseDescr), intent(in) :: denseDesc

    !> Orbital information
    type(TOrbitals), intent(in) :: orb

    !> Integral container
    type(TIntegral), intent(in) :: ints

    !> Holds real and complex delta density matrices
    type(TDensityMatrix), intent(in) :: densityMatrix

    !> List of neighbours for each atom
    type(TNeighbourList), intent(in) :: neighbourList

    !> Number of neighbours for each of the atoms
    integer, intent(in) :: nNeighbourSK(:)

    !> List of neighbours for each atom (symmetric version)
    type(TSymNeighbourList), intent(in) :: symNeighbourList

    !> Nr. of neighbours for each atom.
    integer, intent(in) :: nNeighbourCamSym(:)

    !> Index for which unit cell atoms are associated with
    integer, intent(in) :: iCellVec(:)

    !> Vectors (in units of the lattice constants) to cells of the lattice
    real(dp), intent(in) :: cellVecs(:,:)

    !> Vectors to neighboring unit cells in absolute units
    real(dp), intent(in) :: rCellVecs(:,:)

    !> Index array for the start of atomic blocks in sparse arrays
    integer, intent(in) :: iSparseStart(:,:)

    !> Map from image atoms to the original unique atom
    integer, intent(in) :: img2CentCell(:)

    !> The k-points in relative coordinates to calculate delta H(k) for
    real(dp), intent(in) :: kPoints(:,:)

    !> The k-point weights
    real(dp), intent(in) :: kWeights(:)

    !> Square (unpacked) Hamiltonian for all k-point/spin composite indices to be updated
    complex(dp), intent(out), allocatable :: HSqrCplxCam(:,:,:)

    !> Error status
    type(TStatus), intent(inout) :: errStatus

    call env%globalTimer%startTimer(globalTimers%hybridXcH)

    select case(this%hybridXcAlg)
    case (hybridXcAlgo%thresholdBased)
      @:RAISE_ERROR(errStatus, -1, "Thresholded algorithm not implemented for CAM beyond the Gamma&
          & point.")
    case (hybridXcAlgo%neighbourBased)
      if (this%gammaType == hybridXcGammaTypes%mic) then
        call addCamHamiltonianNeighbour_kpts_mic(this, env, densityMatrix%deltaRhoInCplxHS,&
            & symNeighbourList, nNeighbourCamSym, rCellVecs, cellVecs, denseDesc%iAtomStart, orb,&
            & kPoints, densityMatrix%iKiSToiGlobalKS, HSqrCplxCam, errStatus)
        @:PROPAGATE_ERROR(errStatus)
    else
      call addCamHamiltonianNeighbour_kpts_ct(this, env, densityMatrix%deltaRhoInCplxHS,&
          & symNeighbourList, nNeighbourCamSym, cellVecs, denseDesc%iAtomStart, orb, kPoints,&
          & densityMatrix%iKiSToiGlobalKS, HSqrCplxCam, errStatus)
      @:PROPAGATE_ERROR(errStatus)
      end if
    case (hybridXcAlgo%matrixBased)
      call addCamHamiltonianMatrix_kpts(this, env, denseDesc, ints, densityMatrix, neighbourList,&
          & nNeighbourSK, iCellVec, cellVecs, iSparseStart, img2CentCell, kPoints, kWeights,&
          & HSqrCplxCam)
    end select

    call env%globalTimer%stopTimer(globalTimers%hybridXcH)

  end subroutine getCamHamiltonian_kpts


  !> Adds CAM range-separated contributions to Hamiltonian, using the thresholding algorithm.
  !! (non-periodic version)
  !!
  !! Eq.(31) of J. Chem. Phys. 143, 184107 (2015) (DOI: 10.1063/1.4935095)
  subroutine addCamHamiltonianThreshold_cluster(this, overlap, deltaRho, iNeighbour, nNeighbourCam,&
      & iSquare, hamiltonian, orb)

    !> Class instance
    class(THybridXcFunc), intent(inout) :: this

    !> Square real overlap matrix
    real(dp), intent(in) :: overlap(:,:)

    !> Square density matrix (deltaRho in DFTB terms)
    real(dp), intent(in) :: deltaRho(:,:)

    !> Neighbour indices
    integer, dimension(0:,:), intent(in) :: iNeighbour

    !> Nr. of neighbours for each atom.
    integer, dimension(:), intent(in) :: nNeighbourCam

    !> Mapping atom_number -> number of the first basis function of the atomic block atom_number
    integer, intent(in) :: iSquare(:)

    !> Current Hamiltonian
    real(dp), intent(inout) :: hamiltonian(:,:)

    !> Orbital information
    type(TOrbitals), intent(in) :: orb

    !! Dense matrix descriptor indices
    integer, parameter :: descLen = 3, iStart = 1, iEnd = 2, iNOrb = 3

    real(dp), allocatable :: tmpOvr(:,:), tmpDRho(:,:), testOvr(:,:), tmpDDRho(:,:), tmpDHam(:,:)
    integer, allocatable :: ovrInd(:,:)

    call allocateAndInit()
    call evaluateHamiltonian(tmpDHam)
    call adjointLowerTriangle(tmpDHam)
    this%hprev(:,:) = this%hprev + tmpDHam
    hamiltonian(:,:) = hamiltonian + this%hprev
    this%camEnergy = this%camEnergy + evaluateEnergy_real(this%hprev, tmpDRho)

  contains

    !> Allocate and initialise some necessary arrays
    subroutine allocateAndInit()

      integer :: matrixSize, nAtom0
      real(dp) :: tmp
      integer :: iAtMu, iAtNu, iNeigh

      matrixSize = size(hamiltonian, dim=1)
      nAtom0 = size(this%species0)

      allocate(tmpOvr(matrixSize, matrixSize))
      tmpOvr(:,:) = overlap
      call adjointLowerTriangle(tmpOvr)

      allocate(tmpDHam(matrixSize, matrixSize))
      tmpDHam(:,:) = 0.0_dp

      allocate(tmpDRho(matrixSize, matrixSize))
      tmpDRho(:,:) = deltaRho
      call adjointLowerTriangle(tmpDRho)
      call checkAndInitScreening(this, matrixSize, tmpDRho)

      allocate(tmpDDRho(matrixSize, matrixSize))
      tmpDDRho(:,:) = tmpDRho - this%dRhoPrev
      this%dRhoPrev(:,:) = tmpDRho
      allocate(testOvr(nAtom0, nAtom0))
      testOvr(:,:) = 0.0_dp
      allocate(ovrInd(nAtom0, nAtom0))

      do iAtMu = 1, nAtom0
        do iNeigh = 0, nNeighbourCam(iAtMu)
          iAtNu = iNeighbour(iNeigh, iAtMu)
          tmp = maxval(abs(tmpOvr(iSquare(iAtMu) : iSquare(iAtMu + 1) - 1,&
              & iSquare(iAtNu) : iSquare(iAtNu + 1) - 1)))
          testOvr(iAtMu, iAtNu) = tmp
          testOvr(iAtNu, iAtMu) = tmp
        end do
      end do
      do iAtMu = 1, nAtom0
        call index_heap_sort(ovrInd(iAtMu,:), testOvr(iAtMu,:))
      end do

    end subroutine allocateAndInit


    !> Evaluate the update to Hamiltonian due to change in the DM.
    pure subroutine evaluateHamiltonian(tmpDHam)

      !> Update for the old Hamiltonian on exit
      real(dp), intent(out) :: tmpDHam(:,:)

      integer :: nAtom0
      real(dp) :: pbound, prb
      real(dp) :: tmpvec1(orb%mOrb), tmpvec2(orb%mOrb)
      real(dp) :: tmp, tstbound, gammabatch, gammabatchtmp
      integer :: iAtMu, iAtNu, iAt1, iAt2, iSp1, iSp2, nOrb1, nOrb2
      integer :: kk, ll, jj, ii, mu, nu
      integer, dimension(descLen) :: desc1, desc2, descM, descN

      nAtom0 = size(this%species0)

      pbound = maxval(abs(tmpDDRho))
      tmpDHam = 0.0_dp
      loopMu: do iAtMu = 1, nAtom0
        descM = getDescriptor(iAtMu, iSquare)
        loopKK: do kk = 1, nAtom0
          iAt1 = ovrInd(iAtMu, nAtom0 + 1 - kk)
          desc1 = getDescriptor(iAt1, iSquare)
          iSp1 = this%species0(iAt1)
          nOrb1 = orb%nOrbSpecies(iSp1)
          prb = pbound * testOvr(iAt1, iAtMu)
          if (abs(prb) < this%pScreeningThreshold) then
            exit loopKK
          end if
          loopNu: do iAtNu = 1, iAtMu
            descN = getDescriptor(iAtNu, iSquare)
            gammabatchtmp = this%camGammaEval0(iAtMu, iAtNu) + this%camGammaEval0(iAt1, iAtNu)
            loopLL: do ll = 1, nAtom0
              iAt2 = ovrInd(iAtNu, nAtom0 + 1 - ll)
              iSp2 = this%species0(iAt2)
              nOrb2 = orb%nOrbSpecies(iSp2)
              tstbound = prb * testOvr(iAt2, iAtNu)
              if (abs(tstbound) < this%pScreeningThreshold) then
                exit loopLL
              end if
              desc2 = getDescriptor(iAt2, iSquare)
              gammabatch = (this%camGammaEval0(iAtMu, iAt2) + this%camGammaEval0(iAt1, iAt2)&
                  & + gammabatchtmp)
              gammabatch = -0.125_dp * gammabatch
              ! calculate the Q_AB
              do nu = descN(iStart), descN(iEnd)
                jj = 0
                tmpvec2(1:nOrb2) = tmpOvr(desc2(iStart):desc2(iEnd), nu)
                do ii = desc1(iStart), desc1(iEnd)
                  jj = jj + 1
                  tmpvec1(jj) = sum(tmpvec2(1:nOrb2) * tmpDDRho(ii, desc2(iStart):desc2(iEnd)))
                end do
                tmp = 0.0_dp
                do mu = descM(iStart), descM(iEnd)
                  tmp = sum(tmpOvr(desc1(iStart):desc1(iEnd), mu) * tmpvec1(1:nOrb1))
                  tmpDHam(mu, nu) = tmpDHam(mu, nu) + gammabatch * tmp
                end do
              end do
            end do loopLL
          end do loopNu
        end do loopKK
      end do loopMu

    end subroutine evaluateHamiltonian


    !> Initialise the screening matrices.
    subroutine checkAndInitScreening(this, matrixSize, tmpDRho)

      !> Class instance
      class(THybridXcFunc), intent(inout) :: this

      !> Linear dimension of matrix
      integer, intent(in) :: matrixSize

      !> Delta rho from iteration
      real(dp), intent(in), allocatable :: tmpDRho(:,:)

      if (.not. this%tScreeningInited) then
        allocate(this%hprev(matrixSize, matrixSize))
        allocate(this%dRhoPrev(matrixSize, matrixSize))
        this%hprev(:,:) = 0.0_dp
        this%dRhoPrev(:,:) = tmpDRho
        this%tScreeningInited = .true.
      end if

    end subroutine checkAndInitScreening

  end subroutine addCamHamiltonianThreshold_cluster


  !> Adds CAM range-separated contributions to Hamiltonian, using neighbour-list based algorithm.
  !! (non-periodic and Gamma-only version)
  !!
  !! Non-periodic or Gamma-only approximation of
  !! Eq.(43) of Phys. Rev. Materials 7, 063802 (DOI: 10.1103/PhysRevMaterials.7.063802)
  subroutine addCamHamiltonianNeighbour_real(this, deltaRhoSqr, overSparse, iNeighbour,&
      & nNeighbourCam, iSquare, iPair, orb, img2CentCell, HSqrReal)

    !> Instance of object
    class(THybridXcFunc), intent(inout) :: this

    !> Square (unpacked) density matrix
    real(dp), intent(in) :: deltaRhoSqr(:,:)

    !> Sparse (packed) overlap matrix
    real(dp), intent(in) :: overSparse(:)

    !> Neighbour indices
    integer, intent(in) :: iNeighbour(0:,:)

    !> Nr. of neighbours for each atom
    integer, intent(in) :: nNeighbourCam(:)

    !> Position of each atom in the rows/columns of the square matrices. Shape: (nAtom0)
    integer, intent(in) :: iSquare(:)

    !> Position of each (neighbour, atom) pair in the sparse matrix. Shape: (0:maxNeighbour, nAtom0)
    integer, intent(in) :: iPair(0:,:)

    !> Orbital information
    type(TOrbitals), intent(in) :: orb

    !> Map images of atoms to the central cell
    integer, intent(in) :: img2CentCell(:)

    !> Square (unpacked) Hamiltonian to be updated
    real(dp), intent(inout) :: HSqrReal(:,:)

    !! Dense matrix descriptor indices
    integer, parameter :: descLen = 3, iStart = 1, iEnd = 2, iNOrb = 3

    real(dp), dimension(orb%mOrb**2), target :: Sma, Sam, Snb, Sbn
    real(dp), dimension(orb%mOrb**2), target :: Pab, Pmb, Pan, Pmn
    real(dp), dimension(:,:), pointer :: pSma, pSam, pSnb, pSbn
    real(dp), dimension(:,:), pointer :: pPab, pPmb, pPan, pPmn
    real(dp) :: gamma1, gamma2, gammaTot
    integer :: nAtom0
    integer :: iAtM, iAtN, iAtA, iAtB, iNeighN, iNeighA
    integer :: iAtBfold, iAtMfold
    integer, dimension(descLen) :: descA, descB, descM, descN
    real(dp), dimension(:,:), allocatable :: tmpDRho
    real(dp), dimension(:,:), allocatable, target :: tmpHH

    call allocateAndInit(tmpHH, tmpDRho)
    call evaluateHamiltonian()

    if (this%tSpin .or. this%tREKS) then
      tmpHH(:,:) = 0.25_dp * tmpHH
    else
      tmpHH(:,:) = 0.125_dp * tmpHH
    end if

    call adjointLowerTriangle(tmpHH)

    HSqrReal(:,:) = HSqrReal + tmpHH
    this%camEnergy = this%camEnergy + evaluateEnergy_real(tmpHH, tmpDRho)

  contains

    !> Allocates storage for mapping 1D<->2D array sections.
    subroutine allocateAndInit(tmpHH, tmpDRho)

      !> Density matrix case
      real(dp), dimension(:,:), allocatable, intent(inout) :: tmpDRho

      !> Hamiltonian matrix case
      real(dp), dimension(:,:), allocatable, target, intent(inout) :: tmpHH

      allocate(tmpHH(size(HSqrReal, dim=1), size(HSqrReal, dim=2)))
      tmpHH(:,:) = 0.0_dp
      allocate(tmpDRho(size(deltaRhoSqr, dim=1), size(deltaRhoSqr, dim=1)))
      tmpDRho(:,:) = deltaRhoSqr
      call adjointLowerTriangle(tmpDRho)

    end subroutine allocateAndInit


    !> Actually evaluates the neighbour based cut-off hamiltonian.
    subroutine evaluateHamiltonian()

      nAtom0 = size(this%species0)

      loopN: do iAtN = 1, nAtom0
        descN = getDescriptor(iAtN, iSquare)
        loopB: do iNeighN = 0, nNeighbourCam(iAtN)
          iAtB = iNeighbour(iNeighN, iAtN)
          iAtBfold = img2CentCell(iAtB)
          descB = getDescriptor(iAtBfold, iSquare)
          call copyOverlapBlock(iAtN, iNeighN, descN(iNOrb), descB(iNOrb), Sbn, pSbn)
          call transposeBlock(pSbn, Snb, pSnb)
          loopA: do iAtA = 1, nAtom0
            descA = getDescriptor(iAtA, iSquare)
            call copyDensityBlock(descA, descB, Pab, pPab)
            call copyDensityBlock(descA, descN, Pan, pPan)
            gamma1 = this%camGammaEval0(iAtA, iAtN) + this%camGammaEval0(iAtA, iAtBfold)
            loopM: do iNeighA = 0, nNeighbourCam(iAtA)
              iAtM = iNeighbour(iNeighA, iAtA)
              iAtMfold = img2CentCell(iAtM)
              descM = getDescriptor(iAtMfold, iSquare)
              call copyOverlapBlock(iAtA, iNeighA, descA(iNOrb), descM(iNOrb), Sma, pSma)
              call transposeBlock(pSma, Sam, pSam)
              gamma2 = this%camGammaEval0(iAtMfold, iAtN) + this%camGammaEval0(iAtMfold, iAtBfold)
              gammaTot = gamma1 + gamma2

              if (iAtMfold >= iAtN) then
                call updateHamiltonianBlock(descM, descN, pSma, pSbn, pPab)
              end if
              if (iAtA >= iAtN .and. iAtMfold /= iAtA) then
                call copyDensityBlock(descM, descB, Pmb, pPmb)
                call updateHamiltonianBlock(descA, descN, pSam, pSbn, pPmb)
              end if
              if (iAtMfold >= iAtBfold .and. iAtN /= iAtBfold) then
                call updateHamiltonianBlock(descM, descB, pSma, pSnb, pPan)
              end if
              if (iAtA >= iAtBfold .and. iAtMfold /= iAtA .and. iAtN /= iAtBfold) then
                call copyDensityBlock(descM, descN, Pmn, pPmn)
                call updateHamiltonianBlock(descA, descB, pSam, pSnb, pPmn)
              end if
            end do loopM
          end do loopA
        end do loopB
      end do loopN

    end subroutine evaluateHamiltonian


    !> Copies atom block from sparse matrix.
    pure subroutine copyOverlapBlock(iAt, iNeigh, iNOrbAt, iNOrbNeigh, localBlock, pLocalBlock)

      !> Atom for which this is a neighbour
      integer, intent(in) :: iAt

      !> Number of neighbour for this block
      integer, intent(in) :: iNeigh

      !> Number of orbitals on iAt
      integer, intent(in) :: iNOrbAt

      !> Number of orbitals on neighbour atom
      integer, intent(in) :: iNOrbNeigh

      !> Local block
      real(dp), dimension(:), target, intent(inout) :: localBlock

      !> Pointer to local block
      real(dp), dimension(:,:), pointer, intent(out) :: pLocalBlock

      integer :: ind

      ind = iPair(iNeigh, iAt) + 1
      localBlock(1:iNOrbNeigh*iNOrbAt) = overSparse(ind:ind+iNOrbNeigh*iNOrbAt-1)
      pLocalBlock(1:iNOrbNeigh, 1:iNOrbAt) => localBlock(1:iNOrbNeigh*iNOrbAt)

    end subroutine copyOverlapBlock


    !> Copies density matrix block from sparse matrix.
    pure subroutine copyDensityBlock(desc1, desc2, localBlock, pLocalBlock)

      !> Start, end and range of first block
      integer, dimension(descLen), intent(in) :: desc1

      !> Start, end and range of second block
      integer, dimension(descLen), intent(in) :: desc2

      !> Local block in 1D format
      real(dp), dimension(:), target, intent(inout) :: localBlock

      !> Pointer to local block
      real(dp), dimension(:,:), pointer, intent(out) :: pLocalBlock

      pLocalBlock(1:desc1(iNOrb), 1:desc2(iNOrb)) => localBlock(1:desc1(iNOrb) * desc2(iNOrb))
      pLocalBlock(:,:) = tmpDRho(desc1(iStart):desc1(iEnd), desc2(iStart):desc2(iEnd))

    end subroutine copyDensityBlock


    !> Transposes a block.
    pure subroutine transposeBlock(orig, localBlock, pLocalBlock)

      !> Original matrix block
      real(dp), dimension(:,:), intent(in) :: orig

      !> Local copy in 1D
      real(dp), dimension(:), target, intent(out) :: localBlock

      !> Pointer to local copy
      real(dp), dimension(:,:), pointer, intent(out) :: pLocalBlock

      pLocalBlock(1:size(orig, dim=2), 1:size(orig, dim=1)) => localBlock(1:size(orig))
      pLocalBlock = transpose(orig)

    end subroutine transposeBlock


    !> Adds a contribution to a Hamiltonian block.
    subroutine updateHamiltonianBlock(descM, descN, pSma, pSbN, pPab)

      !> Start, end and range of row
      integer, dimension(descLen), intent(in) :: descM

      !> Start, end and range of column
      integer, dimension(descLen), intent(in) :: descN

      !> First overlap block
      real(dp), dimension(:,:), pointer, intent(in) :: pSma

      !> Second overlap block
      real(dp), dimension(:,:), pointer, intent(in) :: pSbN

      !> Density matrix block
      real(dp), dimension(:,:), pointer, intent(in) :: pPab

      real(dp), dimension(:,:), pointer :: pHmn

      pHmn => tmpHH(descM(iStart):descM(iEnd), descN(iStart):descN(iEnd))
      pHmn(:,:) = pHmn - gammaTot * matmul(matmul(pSma, pPab), pSbn)

    end subroutine updateHamiltonianBlock

  end subroutine addCamHamiltonianNeighbour_real


#:if WITH_SCALAPACK

  !> Update Hamiltonian with CAM range-separated contributions, using a matrix-matrix multiplication
  !! based algorithm.
  !! (real non-periodic and real Gamma-only version)
  !!
  !! Eq.(B3) of Phys. Rev. Materials 7, 063802 (DOI: 10.1103/PhysRevMaterials.7.063802)
  subroutine addCamHamiltonianMatrix_real_blacs(this, env, denseDesc, overlap, densSqr, HH)

    !> Class instance
    class(THybridXcFunc), intent(inout) :: this

    !> Environment settings
    type(TEnvironment), intent(inout) :: env

    !> Dense matrix descriptor
    type(TDenseDescr), intent(in) :: denseDesc

    !> Square (unpacked) overlap matrix.
    real(dp), intent(in) :: overlap(:,:)

    !> Square (unpacked) density matrix
    real(dp), intent(in) :: densSqr(:,:)

    !> Square (unpacked) Hamiltonian to be updated.
    real(dp), intent(inout) :: HH(:,:)

    !! Symmetrized, square (unpacked) Hamiltonian
    real(dp), allocatable :: Hcam(:,:)

    !! Square matrix filled with orbital-resolved gamma values
    !! Actually the diatomic gamma elements are just spread to all orbitals
    real(dp), allocatable :: camGammaAO(:,:)

    !! Number of local rows and columns
    integer :: nLocRow, nLocCol

    !! Local energy contribution
    real(dp) :: Etmp

    nLocRow = size(overlap, dim=1)
    nLocCol = size(overlap, dim=2)

    allocate(camGammaAO(nLocRow, nLocCol))
    allocate(Hcam(nLocRow, nLocCol))

    ! Compared to serial algorithm, no need to symmetrize Hamiltonian, overlap and density matrices
    ! as both triangles should be constructed if MPI parallel

    call initGamma(this, env, denseDesc, camGammaAO)
    call evaluateHamiltonian(this, denseDesc%blacsOrbSqr, overlap, densSqr, camGammaAO, Hcam)

    HH(:,:) = HH + Hcam

    ! Locally stored part of the energy, will collect over rank later:
    Etmp = evaluateEnergy_real(Hcam, densSqr)

    this%camEnergy = this%camEnergy + Etmp

  contains

    !> Get orbital-by-orbital gamma matrix
    subroutine initGamma(this, env, denseDesc, camGammaAO)

      !> Class instance
      class(THybridXcFunc), intent(inout) :: this

      !> Environment settings
      type(TEnvironment), intent(inout) :: env

      !> Dense matrix descriptor
      type(TDenseDescr), intent(in) :: denseDesc

      !> Symmetrized long-range gamma matrix
      real(dp), intent(out) :: camGammaAO(:,:)

      integer :: iAt1, iAt2, ii, jj, iOrb1, iOrb2

      ! Get CAM gamma variable
      camGammaAO(:,:) = 0.0_dp

      do jj = 1, size(camGammaAO, dim=2)
        iOrb2 = scalafx_indxl2g(jj, denseDesc%blacsOrbSqr(NB_), env%blacs%orbitalGrid%mycol,&
            & denseDesc%blacsOrbSqr(CSRC_), env%blacs%orbitalGrid%ncol)
        call bisection(iAt2, denseDesc%iAtomStart, iOrb2)
        do ii = 1, size(camGammaAO, dim=1)
          iOrb1 = scalafx_indxl2g(ii, denseDesc%blacsOrbSqr(MB_), env%blacs%orbitalGrid%myrow,&
              & denseDesc%blacsOrbSqr(RSRC_), env%blacs%orbitalGrid%nrow)
          call bisection(iAt1, denseDesc%iAtomStart, iOrb1)
          camGammaAO(ii, jj) = this%camGammaEval0(iAt1, iAt2)
        end do
      end do

    end subroutine initGamma


    !> Evaluates the Hamiltonian using PGEMM operations.
    subroutine evaluateHamiltonian(this, desc, Smat, Dmat, camGammaAO, Hcam)

      !> Class instance
      class(THybridXcFunc), intent(inout) :: this

      !> BLACS matrix descriptor
      integer, intent(in) :: desc(:)

      !> Symmetrized square overlap matrix
      real(dp), intent(in) :: Smat(:,:)

      !> Symmetrized square density matrix
      real(dp), intent(in) :: Dmat(:,:)

      !> Symmetrized long-range gamma matrix
      real(dp), intent(in) :: camGammaAO(:,:)

      !> Symmetrized long-range Hamiltonian matrix
      real(dp), intent(inout) :: Hcam(:,:)

      !!
      real(dp), allocatable :: Hmat(:,:)

      !!
      real(dp), allocatable :: tmpMat(:,:)

      !! Size of distributed matrices
      integer :: nRows, nCols

      nRows = size(Smat, dim=1)
      nCols = size(Smat, dim=2)

      allocate(Hmat(nRows, nCols), source=0.0_dp)
      allocate(tmpMat(nRows, nCols), source=0.0_dp)

      Hcam(:,:) = 0.0_dp

      call pblasfx_pgemm(Smat, desc, Dmat, desc, tmpMat, desc)
      call pblasfx_pgemm(tmpMat, desc, Smat, desc, Hcam, desc)
      Hcam(:,:) = Hcam * camGammaAO

      tmpMat(:,:) = tmpMat * camGammaAO
      call pblasfx_pgemm(tmpMat, desc, Smat, desc, Hcam, desc, alpha=1.0_dp, beta=1.0_dp)

      Hmat(:,:) = Dmat * camGammaAO
      call pblasfx_pgemm(Smat, desc, hMat, desc, tmpMat, desc)
      call pblasfx_pgemm(tmpMat, desc, Smat, desc, Hcam, desc, alpha=1.0_dp, beta=1.0_dp)

      call pblasfx_pgemm(Dmat, desc, Smat, desc, tmpMat, desc)
      tmpMat(:,:) = tmpMat * camGammaAO
      call pblasfx_pgemm(Smat, desc, tmpMat, desc, Hcam, desc, alpha=1.0_dp, beta=1.0_dp)

      if (this%tSpin .or. this%tREKS) then
        Hcam(:,:) = -0.25_dp * Hcam
      else
        Hcam(:,:) = -0.125_dp * Hcam
      end if

    end subroutine evaluateHamiltonian

  end subroutine addCamHamiltonianMatrix_real_blacs

#:else

  !> Update Hamiltonian with CAM range-separated contributions, using a matrix-matrix multiplication
  !! based algorithm.
  !! (real non-periodic and real Gamma-only version)
  !!
  !! Eq.(B3) of Phys. Rev. Materials 7, 063802 (DOI: 10.1103/PhysRevMaterials.7.063802)
  subroutine addCamHamiltonianMatrix_real(this, iSquare, overlap, densSqr, HH)

    !> Class instance
    class(THybridXcFunc), intent(inout) :: this

    !> Position of each atom in the rows/columns of the square matrices. Shape: (nAtom)
    integer, dimension(:), intent(in) :: iSquare

    !> Square (unpacked) overlap matrix
    real(dp), intent(in) :: overlap(:,:)

    !> Square (unpacked) density matrix
    real(dp), intent(in) :: densSqr(:,:)

    !> Square (unpacked) Hamiltonian to be updated
    real(dp), intent(inout) :: HH(:,:)

    !! Symmetrized, square (unpacked) overlap matrix
    real(dp), allocatable :: Smat(:,:)

    !! Symmetrized, square (unpacked) density matrix
    real(dp), allocatable :: Dmat(:,:)

    !! Symmetrized, square (unpacked) Hamiltonian
    real(dp), allocatable :: Hcam(:,:)

    !! Square matrix filled with orbital-resolved gamma values
    !! Actually the diatomic gamma elements are just spread to all orbitals
    real(dp), allocatable :: camGammaAO(:,:)

    !! Number of orbitals in square matrices
    integer :: nOrb

    nOrb = size(overlap, dim=1)

    allocate(Smat(nOrb, nOrb))
    allocate(Dmat(nOrb, nOrb))
    allocate(camGammaAO(nOrb, nOrb))
    allocate(Hcam(nOrb, nOrb))

    call allocateAndInit(this, iSquare, overlap, densSqr, HH, Smat, Dmat, camGammaAO)
    call evaluateHamiltonian(this, Smat, Dmat, camGammaAO, Hcam)

    HH(:,:) = HH + Hcam
    this%camEnergy = this%camEnergy + evaluateEnergy_real(Hcam, Dmat)

  contains

    !> Set up storage and get orbital-by-orbital gamma matrix
    subroutine allocateAndInit(this, iSquare, overlap, densSqr, HH, Smat, Dmat, camGammaAO)

      !> Class instance
      class(THybridXcFunc), intent(inout) :: this

      !> Position of each atom in the rows/columns of the square matrices. Shape: (nAtom)
      integer, intent(in) :: iSquare(:)

      !> Square (unpacked) overlap matrix
      real(dp), intent(in) :: overlap(:,:)

      !> Square (unpacked) density matrix
      real(dp), intent(in) :: densSqr(:,:)

      !> Square (unpacked) Hamiltonian to be updated
      real(dp), intent(inout) :: HH(:,:)

      !> Symmetrized square overlap matrix
      real(dp), intent(out) :: Smat(:,:)

      !> Symmetrized square density matrix
      real(dp), intent(out) :: Dmat(:,:)

      !> Symmetrized CAM gamma matrix
      real(dp), intent(out) :: camGammaAO(:,:)

      !! Number of atoms
      integer :: nAtom

      !! Indices iterating over atoms
      integer :: iAt, jAt

      nAtom = size(this%camGammaEval0, dim=1)

      ! Symmetrize Hamiltonian, overlap, density matrices
      call adjointLowerTriangle(HH)
      Smat(:,:) = overlap
      call adjointLowerTriangle(Smat)
      Dmat(:,:) = densSqr
      call adjointLowerTriangle(Dmat)

      ! Get CAM gamma variable
      camGammaAO(:,:) = 0.0_dp
      do iAt = 1, nAtom
        do jAt = 1, nAtom
          camGammaAO(iSquare(jAt):iSquare(jAt+1)-1, iSquare(iAt):iSquare(iAt+1)-1)&
              & = this%camGammaEval0(jAt, iAt)
        end do
      end do

    end subroutine allocateAndInit


    !> Evaluates the Hamiltonian, using GEMM operations.
    subroutine evaluateHamiltonian(this, Smat, Dmat, camGammaAO, Hcam)

      !> Class instance
      class(THybridXcFunc), intent(inout) :: this

      !> Symmetrized square overlap matrix
      real(dp), intent(in) :: Smat(:,:)

      !> Symmetrized square density matrix
      real(dp), intent(in) :: Dmat(:,:)

      !> Symmetrized CAM gamma matrix
      real(dp), intent(in) :: camGammaAO(:,:)

      !> Symmetrized CAM Hamiltonian matrix
      real(dp), intent(out) :: Hcam(:,:)

      real(dp), allocatable :: Hmat(:,:)
      real(dp), allocatable :: tmpMat(:,:)

      !! Number of orbitals in square matrices
      integer :: nOrb

      nOrb = size(Smat, dim=1)

      allocate(Hmat(nOrb, nOrb))
      allocate(tmpMat(nOrb, nOrb))

      call symm(tmpMat, 'l', Smat, Dmat)
      call symm(Hcam, 'r', Smat, tmpMat)
      Hcam(:,:) = Hcam * camGammaAO

      tmpMat(:,:) = tmpMat * camGammaAO
      call symm(Hcam, 'r', Smat, tmpMat, alpha=1.0_dp, beta=1.0_dp)

      Hmat(:,:) = Dmat * camGammaAO
      call symm(tmpMat, 'l', Smat, Hmat)
      call symm(Hcam, 'r', Smat, tmpMat, alpha=1.0_dp, beta=1.0_dp)

      call symm(tmpMat, 'l', Dmat, Smat)
      tmpMat(:,:) = tmpMat * camGammaAO
      call symm(Hcam, 'l', Smat, tmpMat, alpha=1.0_dp, beta=1.0_dp)

      if (this%tSpin .or. this%tREKS) then
        Hcam(:,:) = -0.25_dp * Hcam
      else
        Hcam(:,:) = -0.125_dp * Hcam
      end if

    end subroutine evaluateHamiltonian

  end subroutine addCamHamiltonianMatrix_real

#:endif


  !> Update Hamiltonian with CAM range-separated contributions, using a matrix-matrix multiplication
  !! based algorithm.
  !! (complex non-periodic and complex Gamma-only version)
  !!
  !! Eq.(B3) of Phys. Rev. Materials 7, 063802 (DOI: 10.1103/PhysRevMaterials.7.063802)
  subroutine addCamHamiltonianMatrix_cluster_cmplx(this, iSquare, overlap, densSqr, HH)

    !> Class instance
    class(THybridXcFunc), intent(inout) :: this

    !> Position of each atom in the rows/columns of the square matrices. Shape: (nAtom)
    integer, intent(in) :: iSquare(:)

    !> Square (unpacked) overlap matrix
    complex(dp), intent(in) :: overlap(:,:)

    !> Square (unpacked) density matrix
    complex(dp), intent(in) :: densSqr(:,:)

    !> Square (unpacked) Hamiltonian to be updated
    complex(dp), intent(inout) :: HH(:,:)

    !! Symmetrized, square (unpacked) overlap matrix
    complex(dp), allocatable :: Smat(:,:)

    !! Symmetrized, square (unpacked) density matrix
    complex(dp), allocatable :: Dmat(:,:)

    !! Symmetrized, square (unpacked) Hamiltonian
    complex(dp), allocatable :: Hcam(:,:)

    !! Square matrix filled with orbital-resolved gamma values
    !! Actually the diatomic gamma elements are just spread to all orbitals
    real(dp), allocatable :: camGammaAO(:,:)
    complex(dp), allocatable :: gammaCmplx(:,:)

    !! Number of orbitals in square matrices
    integer :: nOrb

    nOrb = size(overlap, dim=1)

    allocate(Smat(nOrb, nOrb))
    allocate(Dmat(nOrb, nOrb))
    allocate(camGammaAO(nOrb, nOrb))
    allocate(gammaCmplx(nOrb, nOrb))
    allocate(Hcam(nOrb, nOrb))

    call allocateAndInit(this, iSquare, overlap, densSqr, HH, Smat, Dmat, camGammaAO, gammaCmplx)

    call evaluateHamiltonian(this, Smat, Dmat, gammaCmplx, Hcam)

    HH(:,:) = HH + Hcam

    this%camEnergy = this%camEnergy + 0.5_dp * real(sum(Dmat * Hcam), dp)

  contains

    subroutine allocateAndInit(this, iSquare, overlap, densSqr, HH, Smat, Dmat, camGammaAO,&
        & gammaCmplx)

      !> Instance
      class(THybridXcFunc), intent(inout) :: this

      !> Position of each atom in the rows/columns of the square matrices. Shape: (nAtom)
      integer, intent(in) :: iSquare(:)

      !> Square (unpacked) overlap matrix
      complex(dp), intent(in) :: overlap(:,:)

      !> Square (unpacked) density matrix
      complex(dp), intent(in) :: densSqr(:,:)

      !> Square (unpacked) Hamiltonian to be updated
      complex(dp), intent(inout) :: HH(:,:)

      !> Symmetrized square overlap matrix
      complex(dp), intent(out) :: Smat(:,:)

      !> Symmetrized square density matrix
      complex(dp), intent(out) :: Dmat(:,:)

      !> Symmetrized CAM gamma matrix
      real(dp), intent(out) :: camGammaAO(:,:)

      !> Symmetrized CAM gamma matrix
      complex(dp), intent(out) :: gammaCmplx(:,:)

      !! Indices iterating over atoms
      integer :: iAt, jAt

      !! Number of atoms
      integer :: nAtom

      nAtom = size(this%camGammaEval0, dim=1)

      !! Symmetrize Hamiltonian, overlap, density matrices
      call adjointLowerTriangle(HH)
      Smat(:,:) = overlap
      call adjointLowerTriangle(Smat)
      Dmat(:,:) = densSqr
      call adjointLowerTriangle(Dmat)

      ! Get CAM gamma variable
      camGammaAO(:,:) = 0.0_dp
      do iAt = 1, nAtom
        do jAt = 1, nAtom
          camGammaAO(iSquare(jAt):iSquare(jAt+1)-1,iSquare(iAt):iSquare(iAt+1)-1) =&
              & this%camGammaEval0(jAt, iAt)
        end do
      end do
      gammaCmplx(:,:) = camGammaAO

    end subroutine allocateAndInit


    subroutine evaluateHamiltonian(this, Smat, Dmat, gammaCmplx, Hcam)

      !> Instance
      class(THybridXcFunc), intent(inout) :: this

      !> Symmetrized square overlap matrix
      complex(dp), intent(in) :: Smat(:,:)

      !> Symmetrized square density matrix
      complex(dp), intent(in) :: Dmat(:,:)

      !> Symmetrized CAM gamma matrix
      complex(dp), intent(in) :: gammaCmplx(:,:)

      !> Symmetrized CAM Hamiltonian matrix
      complex(dp), intent(out) :: Hcam(:,:)

      complex(dp), allocatable :: Hmat(:,:)
      complex(dp), allocatable :: tmpMat(:,:)

      !! Number of orbitals in square matrices
      integer :: nOrb

      nOrb = size(Smat, dim=1)

      allocate(Hmat(nOrb, nOrb))
      allocate(tmpMat(nOrb, nOrb))

      call hemm(tmpMat, 'l', Smat, Dmat)
      call hemm(Hcam, 'r', Smat, tmpMat)
      Hcam(:,:) = Hcam * gammaCmplx

      tmpMat(:,:) = tmpMat * gammaCmplx
      call hemm(Hcam, 'r', Smat, tmpMat, alpha=(1.0_dp,0.0_dp), beta=(1.0_dp,0.0_dp))

      Hmat(:,:) = Dmat * gammaCmplx
      call hemm(tmpMat, 'l', Smat, Hmat)
      call hemm(Hcam, 'r', Smat, tmpMat, alpha=(1.0_dp,0.0_dp), beta=(1.0_dp,0.0_dp))

      call hemm(tmpMat, 'l', Dmat, Smat)
      tmpMat(:,:) = tmpMat * gammaCmplx
      call hemm(Hcam, 'l', Smat, tmpMat, alpha=(1.0_dp,0.0_dp), beta=(1.0_dp,0.0_dp))

      if (this%tSpin) then
        Hcam(:,:) = -0.25_dp * Hcam
      else
        Hcam(:,:) = -0.125_dp * Hcam
      end if

    end subroutine evaluateHamiltonian

  end subroutine addCamHamiltonianMatrix_cluster_cmplx


  !> Adds range-separated contributions to Hamiltonian, using matrix based algorithm.
  !! (k-point version)
  subroutine addCamHamiltonianMatrix_kpts(this, env, denseDesc, ints, densityMatrix, neighbourList,&
      & nNeighbourSK, iCellVec, cellVec, iSparseStart, img2CentCell, kPoints, kWeights, HSqrCplxCam)

    !> Class instance
    class(THybridXcFunc), intent(inout) :: this

    !> Environment settings
    type(TEnvironment), intent(inout) :: env

    !> Dense matrix descriptor
    type(TDenseDescr), intent(in) :: denseDesc

    !> Integral container
    type(TIntegral), intent(in) :: ints

    !> Holds real and complex delta density matrices
    type(TDensityMatrix), intent(in), target :: densityMatrix

    !> List of neighbours for each atom
    type(TNeighbourList), intent(in) :: neighbourList

    !> Number of neighbours for each of the atoms
    integer, intent(in) :: nNeighbourSK(:)

    !> Index for which unit cell atoms are associated with
    integer, intent(in) :: iCellVec(:)

    !> Vectors (in units of the lattice constants) to cells of the lattice
    real(dp), intent(in) :: cellVec(:,:)

    !> Index array for the start of atomic blocks in sparse arrays
    integer, intent(in) :: iSparseStart(:,:)

    !> Map from image atoms to the original unique atom
    integer, intent(in) :: img2CentCell(:)

    !> The k-points in relative coordinates to calculate delta H(k) for
    real(dp), intent(in), target :: kPoints(:,:)

    !> The k-point weights
    real(dp), intent(in), target :: kWeights(:)

    !> Square (unpacked) Hamiltonian for all k-point/spin composite indices
    complex(dp), intent(out), allocatable :: HSqrCplxCam(:,:,:)

    !! Temporary y-storage
    complex(dp), allocatable :: gammaAO(:,:), gammaAOCc(:,:)

    !! Size of square matrices (e.g. Hamiltonian)
    integer :: squareSize

    !! Number of k-point-spin compound indices
    integer :: nKS

    !! K-point indices
    integer :: iK, iKPrime

    !! Spin index
    integer :: iS

    !! Number of k-points and spins
    integer :: nK, nKPrime, nS

    !! Global iKS for arrays present at every MPI rank
    integer :: iGlobalKS, iGlobalKPrimeS

    !! Temporary storage
    complex(dp), allocatable :: tmp(:,:), tmp2(:,:)
    complex(dp), dimension(:,:), allocatable :: Sp_c, dPp_c, Sp_dPp, dPp_Sp, Sp_dPp_cc, dPp_Sp_cc

    !! Composite index for two nested loops over central cell
    integer, allocatable :: iKSComposite(:,:)

    !! Dense, square, k-space overlap matrices S(k) on all MPI processes
    complex(dp), allocatable, target :: SSqrCplx(:,:,:)

    !! Pointer to dense, square, k-space overlap matrices S(k') on all MPI processes
    complex(dp), pointer :: SSqrCplxPrime(:,:,:)

    !! Buffer for dense, square, k-space overlap matrices S(k')
    complex(dp), allocatable, target :: SSqrCplxBuffer(:,:,:)

    !! Pointer to composite index for mapping iKPrime/iS --> iGlobalKPrimeS
    !! (for arrays present at every MPI rank)
    integer, pointer :: iKPrimeiSToiGlobalKPrimeS(:,:)

    !! Buffer for composite index for mapping iKPrime/iS --> iGlobalKPrimeS
    integer, allocatable, target :: iKPrimeiSToiGlobalKPrimeSBuffer(:,:)

    !! Pointer to k' k-points that are possibly different from the current ones in case of a
    !! bandstructure calculation
    real(dp), pointer :: kPointPrime(:,:)

    !> Pointer to weights of k' k-points
    real(dp), pointer :: kWeightPrime(:)

    !! Auxiliary variable
    integer :: ii

    squareSize = size(densityMatrix%deltaRhoInCplx, dim=1)
    nS = size(densityMatrix%iKiSToiGlobalKS, dim=2)
    nK = size(kPoints, dim=2)
    nKS = nK * nS

    ! check and initialize screening
    if (.not. this%tScreeningInited) then
      allocate(this%hprevCplxHS(squareSize, squareSize, nKS), source=(0.0_dp, 0.0_dp))
      this%tScreeningInited = .true.
    end if

    ! allocate exchange contribution to Hamiltonian
    allocate(HSqrCplxCam(squareSize, squareSize, nKS), source=(0.0_dp, 0.0_dp))

    ! skip whole procedure if delta density matrix is close to zero, e.g. in the first SCC iteration
    if (maxval(abs(densityMatrix%deltaRhoInCplx)) < 1e-16_dp) return

    call getDenseSqrDualSpaceOverlap(env, denseDesc, ints, neighbourList, nNeighbourSK, iCellVec,&
        & cellVec, iSparseStart, img2CentCell, kPoints, SSqrCplx)

    if (allocated(densityMatrix%kPointPrime)) then
      @:ASSERT(allocated(densityMatrix%kWeightPrime))
      nKPrime = size(densityMatrix%kPointPrime, dim=2)
      kPointPrime => densityMatrix%kPointPrime
      kWeightPrime => densityMatrix%kWeightPrime
      ! Build spin/k-point composite index for all spins and k-points (global)
      iGlobalKPrimeS = 1
      allocate(iKPrimeiSToiGlobalKPrimeSBuffer(nKPrime, nS))
      do iS = 1, nS
        do iKPrime = 1, nKPrime
          iKPrimeiSToiGlobalKPrimeSBuffer(iKPrime, iS) = iGlobalKPrimeS
          iGlobalKPrimeS = iGlobalKPrimeS + 1
        end do
      end do
      iKPrimeiSToiGlobalKPrimeS => iKPrimeiSToiGlobalKPrimeSBuffer
      call getDenseSqrDualSpaceOverlap(env, denseDesc, ints, neighbourList, nNeighbourSK, iCellVec,&
          & cellVec, iSparseStart, img2CentCell, densityMatrix%kPointPrime, SSqrCplxBuffer)
      SSqrCplxPrime => SSqrCplxBuffer
    else
      nKPrime = nK
      kPointPrime => kPoints
      kWeightPrime => kWeights
      iKPrimeiSToiGlobalKPrimeS => densityMatrix%iKiSToiGlobalKS
      SSqrCplxPrime => SSqrCplx
    end if

    allocate(tmp(squareSize, squareSize))
    allocate(tmp2(squareSize, squareSize))

    allocate(Sp_c(squareSize, squareSize))
    allocate(dPp_c(squareSize, squareSize))

    allocate(Sp_dPp(squareSize, squareSize))
    allocate(dPp_Sp(squareSize, squareSize))
    allocate(Sp_dPp_cc(squareSize, squareSize))
    allocate(dPp_Sp_cc(squareSize, squareSize))

    allocate(gammaAO(squareSize, squareSize))
    allocate(gammaAOCc(squareSize, squareSize))

    ! this is why the number of MPI groups may not exceed nS * nK * nKPrime processes
    ! (there wouldn't be any additional speedup)
    call getiKSiKSPrimeCompositeIndex(env, nS, nK, nKPrime, iKSComposite)

    do ii = 1, size(iKSComposite, dim=2)
      iS = iKSComposite(1, ii)
      iK = iKSComposite(2, ii)
      iGlobalKS = densityMatrix%iKiSToiGlobalKS(iK, iS)

      iKPrime = iKSComposite(3, ii)
      iGlobalKPrimeS = iKPrimeiSToiGlobalKPrimeS(iKPrime, iS)

      call hemm(Sp_dPp, 'l', SSqrCplxPrime(:,:, iKPrime),&
          & densityMatrix%deltaRhoInCplx(:,:, iGlobalKPrimeS))
      call hemm(dPp_Sp, 'l', densityMatrix%deltaRhoInCplx(:,:, iGlobalKPrimeS),&
          & SSqrCplxPrime(:,:, iKPrime))

      ! use conjg(S(k)) = transpose(S(k)) = S(-k)
      ! use conjg(dP(k)) = transpose(dP(k)) = dP(-k)
      Sp_c(:,:) = conjg(SSqrCplxPrime(:,:, iKPrime))
      dPp_c(:,:) = conjg(densityMatrix%deltaRhoInCplx(:,:, iGlobalKPrimeS))
      call hemm(Sp_dPp_cc, 'l', Sp_c, dPp_c)
      call hemm(dPp_Sp_cc, 'l', dPp_c, Sp_c)

      ! spin-polarized case: y-matrix constructed even if k and k' do not change
      call getCamGammaFourierAO(this, denseDesc%iAtomStart, this%cellVecsG, this%rCellVecsG,&
          & kPoints(:, iK), kPointPrime(:, iKPrime), gammaAO)
      call getCamGammaFourierAO(this, denseDesc%iAtomStart, this%cellVecsG, this%rCellVecsG,&
          & kPoints(:, iK), -kPointPrime(:, iKPrime), gammaAOCc)

      ! Term 1
      call hemm(tmp, 'r', SSqrCplxPrime(:,:, iKPrime), Sp_dPp)
      tmp(:,:) = tmp * gammaAO
      ! Term 1 (complex conjugated for inverse k-points)
      ! use conjg(S(k)) = transpose(S(k)) = S(-k)
      call hemm(tmp2, 'r', Sp_c, Sp_dPp_cc)
      tmp2(:,:) = tmp2 * gammaAOCc
      tmp(:,:) = tmp + tmp2

      ! Term 2
      call hemm(tmp, 'r', SSqrCplx(:,:, iK), Sp_dPp * gammaAO, beta=(1.0_dp, 0.0_dp))
      ! Term 2 (complex conjugated for inverse k-points)
      call hemm(tmp, 'r', SSqrCplx(:,:, iK), Sp_dPp_cc * gammaAOCc, beta=(1.0_dp, 0.0_dp))

      ! Term 3
      call hemm(tmp, 'l', SSqrCplx(:,:, iK), dPp_Sp * gammaAO, beta=(1.0_dp, 0.0_dp))
      ! Term 3 (complex conjugated for inverse k-points)
      call hemm(tmp, 'l', SSqrCplx(:,:, iK), dPp_Sp_cc * gammaAOCc, beta=(1.0_dp, 0.0_dp))

      ! Add terms 1-3
      ! (the factor 0.5 accounts for the additional -k' points)
      HSqrCplxCam(:,:, iGlobalKS) = HSqrCplxCam(:,:, iGlobalKS)&
          & + 0.5_dp * kWeightPrime(iKPrime) * tmp

      ! Term 4
      tmp(:,:) = densityMatrix%deltaRhoInCplx(:,:, iGlobalKPrimeS) * gammaAO
      tmp2(:,:) = tmp
      call hemm(tmp, 'r', SSqrCplx(:,:, iK), tmp2)
      tmp2(:,:) = tmp
      call hemm(tmp, 'l', SSqrCplx(:,:, iK), tmp2)

      ! Add term 4
      ! (the factor 0.5 accounts for the additional -k' points)
      HSqrCplxCam(:,:, iGlobalKS) = HSqrCplxCam(:,:, iGlobalKS)&
          & + 0.5_dp * kWeightPrime(iKPrime) * tmp

      ! Term 4 (complex conjugated for inverse k-points)
      ! use conjg(dP(k)) = transpose(dP(k)) = dP(-k)
      tmp(:,:) = dPp_c * gammaAOCc
      tmp2(:,:) = tmp
      call hemm(tmp, 'r', SSqrCplx(:,:, iK), tmp2)
      tmp2(:,:) = tmp
      call hemm(tmp, 'l', SSqrCplx(:,:, iK), tmp2)

      ! Add term 4
      ! (the factor 0.5 accounts for the additional -k' points)
      HSqrCplxCam(:,:, iGlobalKS) = HSqrCplxCam(:,:, iGlobalKS)&
          & + 0.5_dp * kWeightPrime(iKPrime) * tmp

    end do

    if (this%tSpin .or. this%tREKS) then
      HSqrCplxCam(:,:,:) = -0.25_dp * HSqrCplxCam
    else
      HSqrCplxCam(:,:,:) = -0.125_dp * HSqrCplxCam
    end if

  #:if WITH_MPI
    ! Sum up contributions of current MPI group
    call mpifx_allreduceip(env%mpi%interGroupComm, HSqrCplxCam, MPI_SUM)
  #:endif

    ! Safe for the serial case as well, since env%tGlobalLead = .true. for MPI-disabled binary
    if (env%tGlobalLead) then
      this%hprevCplxHS(:,:,:) = HSqrCplxCam
    end if

  #:if WITH_MPI
    call mpifx_bcast(env%mpi%globalComm, this%hprevCplxHS)
  #:endif

  end subroutine addCamHamiltonianMatrix_kpts


  !> Calculates the dense, square, dual-space overlap matrices S(k) on all MPI processes.
  subroutine getDenseSqrDualSpaceOverlap(env, denseDesc, ints, neighbourList, nNeighbourSK,&
      & iCellVec, cellVec, iSparseStart, img2CentCell, kPoint, SSqrCplx)

    !> Environment settings
    type(TEnvironment), intent(inout) :: env

    !> Dense matrix descriptor
    type(TDenseDescr), intent(in) :: denseDesc

    !> Integral container
    type(TIntegral), intent(in) :: ints

    !> List of neighbours for each atom
    type(TNeighbourList), intent(in) :: neighbourList

    !> Number of neighbours for each of the atoms
    integer, intent(in) :: nNeighbourSK(:)

    !> Index for which unit cell atoms are associated with
    integer, intent(in) :: iCellVec(:)

    !> Vectors (in units of the lattice constants) to cells of the lattice
    real(dp), intent(in) :: cellVec(:,:)

    !> Index array for the start of atomic blocks in sparse arrays
    integer, intent(in) :: iSparseStart(:,:)

    !> Map from image atoms to the original unique atom
    integer, intent(in) :: img2CentCell(:)

    !> The k-points
    real(dp), intent(in) :: kPoint(:,:)

    !> Dense, square, k-space overlap matrices S(k) on all MPI processes
    complex(dp), intent(out), allocatable :: SSqrCplx(:,:,:)

    !! Start and end index for MPI parallelization, if applicable
    integer :: iParallelStart, iParallelEnd

    !! Number of k-points
    integer :: nKpoints

    !! Iterates over tile of k-points
    integer :: iK

    nKpoints = size(kPoint, dim=2)
    allocate(SSqrCplx(denseDesc%nOrb, denseDesc%nOrb, nKpoints), source=(0.0_dp, 0.0_dp))

    call getStartAndEndIndex(env, nKpoints, iParallelStart, iParallelEnd)

    ! Pre-generate overlap matrix on all MPI processes
    do iK = iParallelStart, iParallelEnd
      call env%globalTimer%startTimer(globalTimers%sparseToDense)
      call unpackHS(SSqrCplx(:,:, iK), ints%overlap, kPoint(:, iK),&
          & neighbourList%iNeighbour, nNeighbourSK, iCellVec, cellVec, denseDesc%iAtomStart,&
          & iSparseStart, img2CentCell)
      call env%globalTimer%stopTimer(globalTimers%sparseToDense)
      call adjointLowerTriangle(SSqrCplx(:,:, iK))
    end do
  #:if WITH_SCALAPACK
    ! Distribute overlap matrices to all nodes via global summation
    call mpifx_allreduceip(env%mpi%globalComm, SSqrCplx, MPI_SUM)
  #:endif

  end subroutine getDenseSqrDualSpaceOverlap


  !> Adds CAM range-separated contributions to Hamiltonian, using neighbour-list based algorithm.
  !! (k-point version)
  !!
  !! Eq.(43) of Phys. Rev. Materials 7, 063802 (DOI: 10.1103/PhysRevMaterials.7.063802)
  subroutine addCamHamiltonianNeighbour_kpts_mic(this, env, deltaRhoSqr, symNeighbourList,&
      & nNeighbourCamSym, rCellVecs, cellVecs, iSquare, orb, kPoints, iKiSToiGlobalKS, HSqrCplxCam,&
      & errStatus)

    !> Class instance
    class(THybridXcFunc), intent(inout), target :: this

    !> Environment settings
    type(TEnvironment), intent(inout) :: env

    !> Square (unpacked) delta spin-density matrix at BvK real-space shifts
    real(dp), intent(in) :: deltaRhoSqr(:,:,:,:,:,:)

    !> List of neighbours for each atom (symmetric version)
    type(TSymNeighbourList), intent(in) :: symNeighbourList

    !> Nr. of neighbours for each atom
    integer, intent(in) :: nNeighbourCamSym(:)

    !> Vectors to neighboring unit cells in absolute units
    real(dp), intent(in) :: rCellVecs(:,:)

    !> Vectors to neighboring unit cells in relative units
    real(dp), intent(in) :: cellVecs(:,:)

    !> Position of each atom in the rows/columns of the square matrices. Shape: (nAtom)
    integer, intent(in) :: iSquare(:)

    !> Orbital information
    type(TOrbitals), intent(in) :: orb

    !> The k-points in relative coordinates to calculate delta H(k) for
    real(dp), intent(in) :: kPoints(:,:)

    !> Composite index for mapping iK/iS --> iGlobalKS for arrays present at every MPI rank
    integer, intent(in) :: iKiSToiGlobalKS(:,:)

    !> Square (unpacked) Hamiltonian for all k-point/spin composite indices to be updated
    complex(dp), intent(out), allocatable :: HSqrCplxCam(:,:,:)

    !> Error status
    type(TStatus), intent(inout) :: errStatus

    !! Dense matrix descriptor indices
    integer, parameter :: descLen = 3, iStart = 1, iEnd = 2, iNOrb = 3

    !! Atom blocks from sparse, real-space overlap matrices S_{\alpha\mu}, S_{\beta\nu}
    real(dp), pointer :: pSam(:,:), pSbn(:,:)
    real(dp), allocatable :: pSamT(:,:)

    !! \tilde{\gamma}_{\mu\nu}, \tilde{\gamma}_{\mu\beta},
    !! \tilde{\gamma}_{\alpha\nu}, \tilde{\gamma}_{\alpha\beta}
    real(dp), allocatable :: gammaMN(:), gammaMB(:), gammaAN(:), gammaAB(:)

    !! Density matrix block \alpha\beta
    real(dp), allocatable :: Pab(:,:,:,:,:)

    !! Stores start/end index and number of orbitals of square matrices
    integer :: descA(descLen), descB(descLen), descM(descLen), descN(descLen)

    !! Temporary storages
    real(dp), allocatable :: deltaDeltaRhoSqr(:,:,:,:,:,:)

    !! Number of atoms in central cell
    integer :: nAtom0

    !! Real-space \vec{h} and \vec{l} vectors in relative coordinates
    real(dp) :: vecH(3), vecL(3)

    !! Real-space \vec{h} and \vec{l} vectors in absolute coordinates
    real(dp) :: rVecH(3), rVecL(3)

    !! Temporary arrays for gemm operations
    real(dp), dimension(orb%mOrb, orb%mOrb) :: pSamT_Pab, pSamT_Pab_pSbn, Pab_Sbn
    real(dp), dimension(orb%mOrb, orb%mOrb) :: pSamT_Pab_gammaAB, pSamT_Pab_gammaMB_pSbn
    real(dp), dimension(orb%mOrb, orb%mOrb) :: pSamT_Pab_Sbn_gammaAN, pSamT_Pab_gammaAB_pSbn
    complex(dp) :: tot(orb%mOrb, orb%mOrb, size(kPoints, dim=2) * size(deltaRhoSqr, dim=6))

    !! K-point-spin compound index
    integer :: iGlobalKS

    !! K-point index
    integer :: iK

    !! Spin index
    integer :: iS

    !! Atom indices (central cell)
    integer :: iAtM, iAtN

    !! Neighbour indices (+corresponding atom indices)
    integer :: iAtA, iAtB

    !! Folded (to central cell) atom indices
    integer :: iAtAfold, iAtBfold

    !! Sorted (according to max overlap estimates) neighbour indices
    integer :: iNeighMsort, iNeighNsort

    !! Auxiliary variables for setting up 2D pointer to sparse overlap
    integer :: ind, nOrbAt, nOrbNeigh

    !! Size of square matrices (e.g. Hamiltonian)
    integer :: squareSize

    !! Max estimate for difference of square delta rho to previous SCC iteration
    real(dp) :: pMax

    !! Integer BvK index
    integer :: bvKIndex(3)

    !! Phase factor
    complex(dp) :: phase

    !! Species indices
    integer :: iSpM, iSpN, iSpA, iSpB

    !! Iterates over all BvK real-space vectors
    integer :: iG

    !! Composite index for four nested loops
    integer :: ii
    integer, allocatable :: compositeIndex(:,:)

    !! Dummy array with zeros
    real(dp) :: zeros(3)

    !! Number of k-points and spins
    integer :: nK, nS

    !! Number of k-point-spin compound indices
    integer :: nKS

    zeros(:) = 0.0_dp

    tot(:,:,:) = (0.0_dp, 0.0_dp)

    squareSize = size(deltaRhoSqr, dim=1)
    nAtom0 = size(this%species0)
    nS = size(deltaRhoSqr, dim=6)
    nK = size(kPoints, dim=2)
    nKS = nK * nS

    ! check and initialize screening
    if (.not. this%tScreeningInited) then
      allocate(this%hprevCplxHS(squareSize, squareSize, nKS))
      this%hprevCplxHS(:,:,:) = (0.0_dp, 0.0_dp)
      ! there is no previous delta density matrix, therefore just copy over
      deltaDeltaRhoSqr = deltaRhoSqr
      this%tScreeningInited = .true.
    else
      ! allocate and initialize difference of delta rho to previous SCC iteration
      deltaDeltaRhoSqr = deltaRhoSqr - this%dRhoPrevCplxHS
    end if

    pMax = maxval(abs(deltaDeltaRhoSqr))
    ! store delta density matrix
    this%dRhoPrevCplxHS = deltaRhoSqr

    ! allocate delta Hamiltonian
    allocate(HSqrCplxCam(squareSize, squareSize, nKS))
    HSqrCplxCam(:,:,:) = (0.0_dp, 0.0_dp)

    ! skip whole procedure if delta density matrix is close to zero, e.g. in the first SCC iteration
    if (pMax < 1e-16_dp) return

    call getFourLoopCompositeIndex(this, env, nNeighbourCamSym, pMax, compositeIndex, errStatus)
    @:PROPAGATE_ERROR(errStatus)

    allocate(gammaMN(size(this%cellVecsG, dim=2)))
    allocate(gammaMB(size(this%cellVecsG, dim=2)))
    allocate(gammaAN(size(this%cellVecsG, dim=2)))
    allocate(gammaAB(size(this%cellVecsG, dim=2)))

    loopMNBA: do ii = 1, size(compositeIndex, dim=2)
      ! Recover indices from composite
      iAtM = compositeIndex(1, ii)
      iAtN = compositeIndex(2, ii)
      iNeighNsort = compositeIndex(3, ii)
      iNeighMsort = compositeIndex(4, ii)

      ! Former loopMN
      iSpM = this%species0(iAtM)
      iSpN = this%species0(iAtN)
      descM = getDescriptor(iAtM, iSquare)
      descN = getDescriptor(iAtN, iSquare)

      ! Former loopB
      iAtB = symNeighbourList%neighbourList%iNeighbour(iNeighNsort, iAtN)
      iAtBfold = symNeighbourList%img2CentCell(iAtB)
      descB = getDescriptor(iAtBfold, iSquare)
      iSpB = this%species0(iAtBfold)
      ! get real-space \vec{l} for gamma arguments
      rVecL(:) = rCellVecs(:, symNeighbourList%iCellVec(iAtB))
      vecL(:) = cellVecs(:, symNeighbourList%iCellVec(iAtB))
      ! get 2D pointer to Sbn overlap block
      ind = symNeighbourList%iPair(iNeighNsort, iAtN) + 1
      nOrbAt = descN(iNOrb)
      nOrbNeigh = descB(iNOrb)
      pSbn(1:nOrbNeigh, 1:nOrbAt) => this%overSym(ind:ind + nOrbNeigh * nOrbAt - 1)

      ! Former loopA
      iAtA = symNeighbourList%neighbourList%iNeighbour(iNeighMsort, iAtM)
      iAtAfold = symNeighbourList%img2CentCell(iAtA)
      descA = getDescriptor(iAtAfold, iSquare)
      iSpA = this%species0(iAtAfold)
      ! get real-space \vec{h} for gamma arguments
      rVecH(:) = rCellVecs(:, symNeighbourList%iCellVec(iAtA))
      vecH(:) = cellVecs(:, symNeighbourList%iCellVec(iAtA))
      ! get 2D pointer to Sam overlap block
      ind = symNeighbourList%iPair(iNeighMsort, iAtM) + 1
      nOrbAt = descM(iNOrb)
      nOrbNeigh = descA(iNOrb)
      ! S_{\alpha\mu}
      pSam(1:nOrbNeigh, 1:nOrbAt) => this%overSym(ind:ind + nOrbNeigh * nOrbAt - 1)
      ! S_{\mu\alpha}
      pSamT = transpose(pSam(1:nOrbNeigh, 1:nOrbAt))

      gammaMN(:) = getCamGammaGResolved(this, iAtM, iAtN, iSpM, iSpN, this%rCellVecsG,&
          & rVecH - rVecL)
      gammaMB(:) = getCamGammaGResolved(this, iAtM, iAtBfold, iSpM, iSpB, this%rCellVecsG, rVecH)
      gammaAN(:) = getCamGammaGResolved(this, iAtAfold, iAtN, iSpA, iSpN, this%rCellVecsG, -rVecL)
      gammaAB(:) = getCamGammaGResolved(this, iAtAfold, iAtBfold, iSpA, iSpB, this%rCellVecsG,&
          & zeros)

      tot(:,:,:) = (0.0_dp, 0.0_dp)

      loopS: do iS = 1, nS

        ! get continuous 2D copy of Pab density matrix block
        Pab = deltaDeltaRhoSqr(descA(iStart):descA(iEnd), descB(iStart):descB(iEnd), :,:,:, iS)

        loopG: do iG = 1, size(this%cellVecsG, dim=2)
          bvKIndex(:) = this%foldToBvKIndex(-this%cellVecsG(:, iG))

          ! term #1/2
          pSamT_Pab(1:descM(iNOrb), 1:descB(iNOrb)) = matmul(pSamT, Pab(:,:, bvKIndex(1),&
              & bvKIndex(2), bvKIndex(3)))

          ! term #1
          pSamT_Pab_pSbn(1:descM(iNOrb), 1:descN(iNOrb)) = matmul(pSamT_Pab(1:descM(iNOrb),&
              & 1:descB(iNOrb)), pSbn)

          ! term #2
          pSamT_Pab_gammaMB_pSbn(1:descM(iNOrb), 1:descN(iNOrb))&
              & = matmul(pSamT_Pab(1:descM(iNOrb), 1:descB(iNOrb)), pSbn * gammaMB(iG))

          ! term #3
          Pab_Sbn(1:descA(iNOrb), 1:descN(iNOrb)) = matmul(Pab(:,:, bvKIndex(1), bvKIndex(2),&
              & bvKIndex(3)), pSbn)
          pSamT_Pab_Sbn_gammaAN(1:descM(iNOrb), 1:descN(iNOrb))&
              & = matmul(pSamT, Pab_Sbn(1:descA(iNOrb), 1:descN(iNOrb)) * gammaAN(iG))

          ! term #4
          pSamT_Pab_gammaAB(1:descM(iNorb), 1:descB(iNorb)) = matmul(pSamT, Pab(:,:,&
              & bvKIndex(1), bvKIndex(2), bvKIndex(3)) * gammaAB(iG))
          pSamT_Pab_gammaAB_pSbn(1:descM(iNOrb), 1:descN(iNOrb))&
              & = matmul(pSamT_Pab_gammaAB(1:descM(iNOrb), 1:descB(iNOrb)), pSbn)

          loopK: do iK = 1, nK

            iGlobalKS = iKiSToiGlobalKS(iK, iS)

            phase = exp(-imag * dot_product(2.0_dp * pi * kPoints(:, iK),&
                & this%cellVecsG(:, iG) - vecL + vecH))

            ! term #1
            tot(1:descM(iNOrb), 1:descN(iNOrb), iGlobalKS)&
                & = tot(1:descM(iNOrb), 1:descN(iNOrb), iGlobalKS)&
                & + cmplx(pSamT_Pab_pSbn(1:descM(iNOrb), 1:descN(iNOrb)) * gammaMN(iG), 0, dp)&
                & * phase

            ! term #2
            tot(1:descM(iNOrb), 1:descN(iNOrb), iGlobalKS)&
                & = tot(1:descM(iNOrb), 1:descN(iNOrb), iGlobalKS)&
                & + cmplx(pSamT_Pab_gammaMB_pSbn(1:descM(iNOrb), 1:descN(iNOrb)), 0, dp) * phase

            ! term #3
            tot(1:descM(iNOrb), 1:descN(iNOrb), iGlobalKS)&
                & = tot(1:descM(iNOrb), 1:descN(iNOrb), iGlobalKS)&
                & + cmplx(pSamT_Pab_Sbn_gammaAN(1:descM(iNOrb), 1:descN(iNOrb)), 0, dp) * phase

            ! term #4
            tot(1:descM(iNOrb), 1:descN(iNOrb), iGlobalKS)&
                & = tot(1:descM(iNOrb), 1:descN(iNOrb), iGlobalKS)&
                & + cmplx(pSamT_Pab_gammaAB_pSbn(1:descM(iNOrb), 1:descN(iNOrb)), 0, dp) * phase

          end do loopK

        end do loopG

      end do loopS

      HSqrCplxCam(descM(iStart):descM(iEnd), descN(iStart):descN(iEnd), :)&
          & = HSqrCplxCam(descM(iStart):descM(iEnd), descN(iStart):descN(iEnd), :)&
          & + tot(1:descM(iNOrb), 1:descN(iNOrb), :)

    end do loopMNBA

    if (this%tSpin .or. this%tREKS) then
      HSqrCplxCam(:,:,:) = -0.25_dp * HSqrCplxCam
    else
      HSqrCplxCam(:,:,:) = -0.125_dp * HSqrCplxCam
    end if

  #:if WITH_MPI
    ! Sum up contributions of current MPI group
    call mpifx_allreduceip(env%mpi%globalComm, HSqrCplxCam, MPI_SUM)
  #:endif

    ! Safe for the serial case as well, since env%tGlobalLead = .true. for MPI-disabled binary
    if (env%tGlobalLead) then
      this%hprevCplxHS(:,:,:) = this%hprevCplxHS + HSqrCplxCam
      HSqrCplxCam(:,:,:) = this%hprevCplxHS
    end if

  #:if WITH_MPI
    call mpifx_bcast(env%mpi%globalComm, this%hprevCplxHS)
    call mpifx_bcast(env%mpi%globalComm, HSqrCplxCam)
  #:endif

  end subroutine addCamHamiltonianNeighbour_kpts_mic


  !> Adds range-separated contributions to Hamiltonian, using neighbour-list based algorithm.
  !! (k-point version)
  !!
  !! Eq.(43) of Phys. Rev. Materials 7, 063802 (DOI: 10.1103/PhysRevMaterials.7.063802)
  subroutine addCamHamiltonianNeighbour_kpts_ct(this, env, deltaRhoSqr, symNeighbourList,&
      & nNeighbourCamSym, cellVecs, iSquare, orb, kPoints, iKiSToiGlobalKS, HSqrCplxCam, errStatus)

    !> Class instance
    class(THybridXcFunc), intent(inout), target :: this

    !> Environment settings
    type(TEnvironment), intent(inout) :: env

    !> Square (unpacked) delta spin-density matrix at BvK real-space shifts
    real(dp), intent(in) :: deltaRhoSqr(:,:,:,:,:,:)

    !> List of neighbours for each atom (symmetric version)
    type(TSymNeighbourList), intent(in) :: symNeighbourList

    !> Nr. of neighbours for each atom
    integer, intent(in) :: nNeighbourCamSym(:)

    !> Vectors to neighboring unit cells in relative units
    real(dp), intent(in) :: cellVecs(:,:)

    !> Position of each atom in the rows/columns of the square matrices. Shape: (nAtom)
    integer, intent(in) :: iSquare(:)

    !> Orbital information
    type(TOrbitals), intent(in) :: orb

    !> The k-points in relative coordinates to calculate delta H(k) for
    real(dp), intent(in) :: kPoints(:,:)

    !> Composite index for mapping iK/iS --> iGlobalKS for arrays present at every MPI rank
    integer, intent(in) :: iKiSToiGlobalKS(:,:)

    !> Square (unpacked) Hamiltonian for all k-point/spin composite indices
    complex(dp), intent(out), allocatable :: HSqrCplxCam(:,:,:)

    !> Error status
    type(TStatus), intent(inout) :: errStatus

    !! Dense matrix descriptor indices
    integer, parameter :: descLen = 3, iStart = 1, iEnd = 2, iNOrb = 3

    !! Atom blocks from sparse, real-space overlap matrices S_{\alpha\mu}, S_{\beta\nu}
    real(dp), pointer :: pSam(:,:), pSbn(:,:)
    real(dp), allocatable :: pSamT(:,:)

    !! Density matrix block \alpha\beta
    real(dp), allocatable :: Pab(:,:,:,:,:)

    !! Stores start/end index and number of orbitals of square matrices
    integer :: descA(descLen), descB(descLen), descM(descLen), descN(descLen)

    !! Temporary storages
    real(dp), allocatable :: deltaDeltaRhoSqr(:,:,:,:,:,:)

    !! Number of atoms in central cell
    integer :: nAtom0

    !! Real-space \vec{h} and \vec{l} vectors in relative coordinates
    real(dp) :: vecH(3), vecL(3)

    !! Temporary arrays for gemm operations
    real(dp), dimension(orb%mOrb, orb%mOrb) :: pSamT_Pab, pSamT_Pab_pSbn, Pab_Sbn
    real(dp), dimension(orb%mOrb, orb%mOrb) :: pSamT_Pab_gammaAB, pSamT_Pab_gammaMB_pSbn
    real(dp), dimension(orb%mOrb, orb%mOrb) :: pSamT_Pab_Sbn_gammaAN, pSamT_Pab_gammaAB_pSbn
    complex(dp) :: tot(orb%mOrb, orb%mOrb, size(kPoints, dim=2) * size(deltaRhoSqr, dim=6))

    !! K-point-spin compound index
    integer :: iGlobalKS

    !! K-point index
    integer :: iK

    !! Spin index
    integer :: iS

    !! Atom indices (central cell)
    integer :: iAtM, iAtN

    !! Neighbour indices (+corresponding atom indices)
    integer :: iAtA, iAtB

    !! Folded (to central cell) atom indices
    integer :: iAtAfold, iAtBfold

    !! Sorted (according to max overlap estimates) neighbour indices
    integer :: iNeighMsort, iNeighNsort

    !! Auxiliary variables for setting up 2D pointer to sparse overlap
    integer :: ind, nOrbAt, nOrbNeigh

    !! Size of square matrices (e.g. Hamiltonian)
    integer :: squareSize

    !! Max estimate for difference of square delta rho to previous SCC iteration
    real(dp) :: pMax

    !! Integer BvK index
    integer :: bvKIndex(3)

    !! Phase factor
    complex(dp) :: phase

    !! Iterates over all BvK real-space vectors
    integer :: iG, iGMN, iGMB, iGAN, iGAB

    !! Composite index for four nested loops
    integer :: ii
    integer, allocatable :: compositeIndex(:,:)

    !! Number of k-points and spins
    integer :: nK, nS

    !! Number of k-point-spin compound indices
    integer :: nKS

    tot(:,:,:) = (0.0_dp, 0.0_dp)

    squareSize = size(deltaRhoSqr, dim=1)
    nAtom0 = size(this%species0)
    nS = size(deltaRhoSqr, dim=6)
    nK = size(kPoints, dim=2)
    nKS = nK * nS

    ! check and initialize screening
    if (.not. this%tScreeningInited) then
      allocate(this%hprevCplxHS(squareSize, squareSize, nKS))
      this%hprevCplxHS(:,:,:) = (0.0_dp, 0.0_dp)
      ! there is no previous delta density matrix, therefore just copy over
      deltaDeltaRhoSqr = deltaRhoSqr
      this%tScreeningInited = .true.
    else
      ! allocate and initialize difference of delta rho to previous SCC iteration
      deltaDeltaRhoSqr = deltaRhoSqr - this%dRhoPrevCplxHS
    end if

    pMax = maxval(abs(deltaDeltaRhoSqr))
    ! store delta density matrix
    this%dRhoPrevCplxHS = deltaRhoSqr

    ! allocate delta Hamiltonian
    allocate(HSqrCplxCam(squareSize, squareSize, nKS))
    HSqrCplxCam(:,:,:) = (0.0_dp, 0.0_dp)

    ! skip whole procedure if delta density matrix is close to zero, e.g. in the first SCC iteration
    if (pMax < 1e-16_dp) return

    call getFourLoopCompositeIndex(this, env, nNeighbourCamSym, pMax, compositeIndex, errStatus)
    @:PROPAGATE_ERROR(errStatus)

    loopMNBA: do ii = 1, size(compositeIndex, dim=2)
      ! Recover indices from composite
      iAtM = compositeIndex(1, ii)
      iAtN = compositeIndex(2, ii)
      iNeighNsort = compositeIndex(3, ii)
      iNeighMsort = compositeIndex(4, ii)

      ! Former loopMN
      descM = getDescriptor(iAtM, iSquare)
      descN = getDescriptor(iAtN, iSquare)

      ! Former loopB
      iAtB = symNeighbourList%neighbourList%iNeighbour(iNeighNsort, iAtN)
      iAtBfold = symNeighbourList%img2CentCell(iAtB)
      descB = getDescriptor(iAtBfold, iSquare)
      ! get real-space \vec{l} for gamma arguments
      vecL(:) = cellVecs(:, symNeighbourList%iCellVec(iAtB))
      ! get 2D pointer to Sbn overlap block
      ind = symNeighbourList%iPair(iNeighNsort, iAtN) + 1
      nOrbAt = descN(iNOrb)
      nOrbNeigh = descB(iNOrb)
      pSbn(1:nOrbNeigh, 1:nOrbAt) => this%overSym(ind:ind + nOrbNeigh * nOrbAt - 1)

      ! Former loopA
      iAtA = symNeighbourList%neighbourList%iNeighbour(iNeighMsort, iAtM)
      iAtAfold = symNeighbourList%img2CentCell(iAtA)
      descA = getDescriptor(iAtAfold, iSquare)
      ! get real-space \vec{h} for gamma arguments
      vecH(:) = cellVecs(:, symNeighbourList%iCellVec(iAtA))
      ! get 2D pointer to Sam overlap block
      ind = symNeighbourList%iPair(iNeighMsort, iAtM) + 1
      nOrbAt = descM(iNOrb)
      nOrbNeigh = descA(iNOrb)
      ! S_{\alpha\mu}
      pSam(1:nOrbNeigh, 1:nOrbAt) => this%overSym(ind:ind + nOrbNeigh * nOrbAt - 1)
      pSamT = transpose(pSam(1:nOrbNeigh, 1:nOrbAt))

      tot(:,:,:) = (0.0_dp, 0.0_dp)

      loopS: do iS = 1, nS

        ! get continuous 2D copy of Pab density matrix block
        Pab = deltaDeltaRhoSqr(descA(iStart):descA(iEnd), descB(iStart):descB(iEnd), :,:,:, iS)

        loopGMN: do iG = 1, size(this%nNonZeroGammaG(iAtM, iAtN)%data)
          iGMN = this%nNonZeroGammaG(iAtM, iAtN)%data(iG)
          bvKIndex(:) = this%foldToBvKIndex(vecH - vecL - this%cellVecsG(:, iGMN))

          pSamT_Pab(1:descM(iNOrb), 1:descB(iNOrb)) = matmul(pSamT,&
              & Pab(:,:, bvKIndex(1), bvKIndex(2), bvKIndex(3)))
          pSamT_Pab_pSbn(1:descM(iNOrb), 1:descN(iNOrb)) = matmul(pSamT_Pab(1:descM(iNOrb),&
              & 1:descB(iNOrb)), pSbn)

          loopKMN: do iK = 1, nK
            iGlobalKS = iKiSToiGlobalKS(iK, iS)
            phase = exp(-imag * dot_product(2.0_dp * pi * kPoints(:, iK),&
                & this%cellVecsG(:, iGMN)))

            ! term #1
            tot(1:descM(iNOrb), 1:descN(iNOrb), iGlobalKS)&
                & = tot(1:descM(iNOrb), 1:descN(iNOrb), iGlobalKS)&
                & + cmplx(pSamT_Pab_pSbn(1:descM(iNOrb), 1:descN(iNOrb))&
                & * this%camGammaEvalG(iAtM, iAtN)%data(iG), 0, dp) * phase
          end do loopKMN
        end do loopGMN

        loopGMB: do iG = 1, size(this%nNonZeroGammaG(iAtM, iAtBfold)%data)
          iGMB = this%nNonZeroGammaG(iAtM, iAtBfold)%data(iG)
          bvKIndex(:) = this%foldToBvKIndex(vecH - this%cellVecsG(:, iGMB))

          pSamT_Pab(1:descM(iNOrb), 1:descB(iNOrb)) = matmul(pSamT,&
              & Pab(:,:, bvKIndex(1), bvKIndex(2), bvKIndex(3)))
          pSamT_Pab_gammaMB_pSbn(1:descM(iNOrb), 1:descN(iNOrb))&
              & = matmul(pSamT_Pab(1:descM(iNOrb), 1:descB(iNOrb)),&
              & pSbn * this%camGammaEvalG(iAtM, iAtBfold)%data(iG))

          loopKMB: do iK = 1, nK
            iGlobalKS = iKiSToiGlobalKS(iK, iS)
            phase = exp(-imag * dot_product(2.0_dp * pi * kPoints(:, iK),&
                & this%cellVecsG(:, iGMB) - vecL))

            ! term #2
            tot(1:descM(iNOrb), 1:descN(iNOrb), iGlobalKS)&
                & = tot(1:descM(iNOrb), 1:descN(iNOrb), iGlobalKS)&
                & + cmplx(pSamT_Pab_gammaMB_pSbn(1:descM(iNOrb), 1:descN(iNOrb)), 0, dp) * phase
          end do loopKMB
        end do loopGMB

        loopGAN: do iG = 1, size(this%nNonZeroGammaG(iAtAfold, iAtN)%data)
          iGAN = this%nNonZeroGammaG(iAtAfold, iAtN)%data(iG)
          bvKIndex(:) = this%foldToBvKIndex(-this%cellVecsG(:, iGAN) - vecL)

          Pab_Sbn(1:descA(iNOrb), 1:descN(iNOrb)) = matmul(Pab(:,:, bvKIndex(1), bvKIndex(2),&
              & bvKIndex(3)), pSbn)
          pSamT_Pab_Sbn_gammaAN(1:descM(iNOrb), 1:descN(iNOrb))&
              & = matmul(pSamT, Pab_Sbn(1:descA(iNOrb), 1:descN(iNOrb))&
              & * this%camGammaEvalG(iAtAfold, iAtN)%data(iG))

          loopKAN: do iK = 1, nK
            iGlobalKS = iKiSToiGlobalKS(iK, iS)
            phase = exp(-imag * dot_product(2.0_dp * pi * kPoints(:, iK),&
                & this%cellVecsG(:, iGAN) + vecH))

            ! term #3
            tot(1:descM(iNOrb), 1:descN(iNOrb), iGlobalKS)&
                & = tot(1:descM(iNOrb), 1:descN(iNOrb), iGlobalKS)&
                & + cmplx(pSamT_Pab_Sbn_gammaAN(1:descM(iNOrb), 1:descN(iNOrb)), 0, dp) * phase
          end do loopKAN
        end do loopGAN

        loopGAB: do iG = 1, size(this%nNonZeroGammaG(iAtAfold, iAtBfold)%data)
          iGAB = this%nNonZeroGammaG(iAtAfold, iAtBfold)%data(iG)
          bvKIndex(:) = this%foldToBvKIndex(-this%cellVecsG(:, iGAB))

          pSamT_Pab_gammaAB(1:descM(iNorb), 1:descB(iNorb)) = matmul(pSamT, Pab(:,:,&
              & bvKIndex(1), bvKIndex(2), bvKIndex(3))&
              & * this%camGammaEvalG(iAtAfold, iAtBfold)%data(iG))
          pSamT_Pab_gammaAB_pSbn(1:descM(iNOrb), 1:descN(iNOrb))&
              & = matmul(pSamT_Pab_gammaAB(1:descM(iNOrb), 1:descB(iNOrb)), pSbn)

          loopKAB: do iK = 1, nK
            iGlobalKS = iKiSToiGlobalKS(iK, iS)
            phase = exp(-imag * dot_product(2.0_dp * pi * kPoints(:, iK),&
                & this%cellVecsG(:, iGAB) - vecL + vecH))

            ! term #4
            tot(1:descM(iNOrb), 1:descN(iNOrb), iGlobalKS)&
                & = tot(1:descM(iNOrb), 1:descN(iNOrb), iGlobalKS)&
                & + cmplx(pSamT_Pab_gammaAB_pSbn(1:descM(iNOrb), 1:descN(iNOrb)), 0, dp) * phase
          end do loopKAB
        end do loopGAB

      end do loopS

      HSqrCplxCam(descM(iStart):descM(iEnd), descN(iStart):descN(iEnd), :)&
          & = HSqrCplxCam(descM(iStart):descM(iEnd), descN(iStart):descN(iEnd), :)&
          & + tot(1:descM(iNOrb), 1:descN(iNOrb), :)

    end do loopMNBA

    if (this%tSpin .or. this%tREKS) then
      HSqrCplxCam(:,:,:) = -0.25_dp * HSqrCplxCam
    else
      HSqrCplxCam(:,:,:) = -0.125_dp * HSqrCplxCam
    end if

  #:if WITH_MPI
    ! Sum up contributions of current MPI group
    call mpifx_allreduceip(env%mpi%globalComm, HSqrCplxCam, MPI_SUM)
  #:endif

    ! Safe for the serial case as well, since env%tGlobalLead = .true. for MPI-disabled binary
    if (env%tGlobalLead) then
      this%hprevCplxHS(:,:,:) = this%hprevCplxHS + HSqrCplxCam
      HSqrCplxCam(:,:,:) = this%hprevCplxHS
    end if

  #:if WITH_MPI
    call mpifx_bcast(env%mpi%globalComm, this%hprevCplxHS)
    call mpifx_bcast(env%mpi%globalComm, HSqrCplxCam)
  #:endif

  end subroutine addCamHamiltonianNeighbour_kpts_ct


  !> Adds the hybrid functional contribution to the total energy.
  subroutine addHybridEnergy_real(this, env, energy)

    !> Class instance
    class(THybridXcFunc), intent(inout) :: this

    !> Environment settings
    type(TEnvironment), intent(in) :: env

    !> Total energy
    real(dp), intent(inout) :: energy

  #:if WITH_MPI
    call mpifx_allreduceip(env%mpi%globalComm, this%camEnergy, MPI_SUM)
  #:endif

    energy = energy + this%camEnergy

    ! hack for spin unrestricted calculation
    ! probably broken for MPI groups over spin
    this%camEnergy = 0.0_dp

  end subroutine addHybridEnergy_real


  !> Adds the hybrid functional contribution to the total energy.
  subroutine addHybridEnergy_kpts(this, env, localKS, iKiSToiGlobalKS, kWeights, deltaRhoOutCplx,&
      & energy)

    !> Class instance
    class(THybridXcFunc), intent(inout) :: this

    !> Environment settings
    type(TEnvironment), intent(in) :: env

    !> The (K, S) tuples of the local processor group (localKS(1:2,iKS))
    !> Usage: iK = localKS(1, iKS); iS = localKS(2, iKS)
    integer, intent(in) :: localKS(:,:)

    !> Composite index for mapping iK/iS --> iGlobalKS for arrays present at every MPI rank
    integer, intent(in) :: iKiSToiGlobalKS(:,:)

    !> The k-point weights
    real(dp), intent(in) :: kWeights(:)

    !> Complex, dense, square k-space delta density matrix of all spins/k-points
    complex(dp), intent(in) :: deltaRhoOutCplx(:,:,:)

    !> Total energy
    real(dp), intent(inout) :: energy

    !! Spin and k-point indices
    integer :: iS, iK

    !! Local and global iKS for arrays present at every MPI rank
    integer :: iLocalKS, iGlobalKS, iDensMatKS

    do iLocalKS = 1, size(localKS, dim=2)
      iK = localKS(1, iLocalKS)
      iS = localKS(2, iLocalKS)
      iGlobalKS = iKiSToiGlobalKS(iK, iS)
      if (this%hybridXcAlg == hybridXcAlgo%matrixBased) then
        iDensMatKS = iGlobalKS
      else
        iDensMatKS = iLocalKS
      end if
      this%camEnergy = this%camEnergy&
          & + evaluateEnergy_cplx(this%hprevCplxHS(:,:, iGlobalKS), kWeights(iK),&
          & transpose(deltaRhoOutCplx(:,:, iDensMatKS)))
    end do

  #:if WITH_MPI
    call mpifx_allreduceip(env%mpi%globalComm, this%camEnergy, MPI_SUM)
  #:endif

    energy = energy + this%camEnergy

    ! hack for spin unrestricted calculation
    this%camEnergy = 0.0_dp

  end subroutine addHybridEnergy_kpts


  !> Finds location of relevant atomic block indices in a dense matrix.
  pure function getDescriptor(iAt, iSquare) result(desc)

    !> Relevant atom
    integer, intent(in) :: iAt

    !> Indexing array for start of atom orbitals
    integer, intent(in) :: iSquare(:)

    !> Resulting location ranges
    integer :: desc(3)

    desc(:) = [iSquare(iAt), iSquare(iAt + 1) - 1, iSquare(iAt + 1) - iSquare(iAt)]

  end function getDescriptor


  !> Evaluates energy from Hamiltonian and density matrix (real version).
  pure function evaluateEnergy_real(hamiltonian, densityMat) result(energy)

    !> Hamiltonian matrix
    real(dp), intent(in) :: hamiltonian(:,:)

    !> Density matrix
    real(dp), intent(in) :: densityMat(:,:)

    !> Resulting energy due to CAM contribution
    real(dp) :: energy

    integer :: mu

    energy = 0.5_dp * sum(hamiltonian * densityMat)

  end function evaluateEnergy_real


  !> Evaluates energy from the Hamiltonian and density matrix (complex version).
  pure function evaluateEnergy_cplx(hamiltonian, kWeight, densityMat) result(energy)

    !> Hamiltonian matrix
    complex(dp), intent(in) :: hamiltonian(:,:)

    !> The k-point weight (for energy contribution)
    real(dp), intent(in) :: kWeight

    !> Density matrix in k-space
    complex(dp), intent(in) :: densityMat(:,:)

    !> Resulting energy due to CAM contribution
    real(dp) :: energy

    energy = 0.5_dp * real(sum(hamiltonian * densityMat), dp) * kWeight

  end function evaluateEnergy_cplx


  !> Returns the value of a polynomial of 5th degree at x (or its derivative).
  !! The polynomial is created with the following boundary conditions:
  !! Its value, its 1st and 2nd derivatives are zero at x = rCut and agree with the provided values
  !! at x = rDamp, i.e. x = rCut - delta.
  !! WARNING: To avoid additional branches, there are no consistency checks, e.g. that rDamp < rCut
  pure function poly5zero(y0, y0p, y0pp, xx, rDamp, rCut, tDerivative) result(yy)

    !> Value of the polynom at x = rDamp
    real(dp), intent(in) :: y0

    !> Value of the 1st derivative at x = rDamp
    real(dp), intent(in) :: y0p

    !> Value of the 2nd derivative at x = rDamp
    real(dp), intent(in) :: y0pp

    !> Point where the polynomial should be calculated
    real(dp), intent(in) :: xx

    !> Point, where the polynomials value and first two derivatives should take the provided values
    real(dp), intent(in) :: rDamp

    !> Point, where the polynomial (and its 1st/2nd derivative) is supposed to be zero
    real(dp), intent(in) :: rCut

    !> True, if the derivative at xx is desired, otherwise the function value it returned
    logical, intent(in), optional :: tDerivative

    !! True, if the derivative at xx is desired, otherwise the function value it returned
    !! Default: .false.
    logical :: tPrime

    !! Value of the polynomial at xx (in general should satisfy rDamp <= x <= rCut)
    real(dp) :: yy

    !! Polynomial coefficients
    real(dp) :: aa, bb, cc, dd, ee, ff

    !! rCut - rDamp
    real(dp) :: delta

    if (present(tDerivative)) then
      tPrime = tDerivative
    else
      tPrime = .false.
    end if

    ! 5th order polynomial definition
    ! f(x) = ax^5 + bx^4 + cx^3 + dx^2 + ex + f
    ! f'(x) = 5ax^4 + 4bx^3 + 3cx^2 + 2dx + e
    ! f''(x) = 20ax^3 + 12bx^2 + 6cx + 2d

    ! boundary conditions:
    ! f(rCut) = 0, f'(rCut) = 0, f''(rCut) = 0
    ! f(rDamp) = y0, f'(rDamp) = y0p, f''(rDamp) = y0pp

    delta = rCut - rDamp

    aa = -(delta**2 * y0pp + 6.0_dp * delta * y0p + 12.0_dp * y0) / (2.0_dp * delta**5)

    bb = (rCut * (5.0_dp * delta**2 * y0pp + 30.0_dp * delta * y0p + 60.0_dp * y0)&
        & - 2.0_dp * delta**3 * y0pp - 14.0_dp * delta**2 * y0p - 30.0_dp * delta * y0)&
        & / (2.0_dp * delta**5)

    cc = -(rCut * (-8.0_dp * delta**3 * y0pp - 56.0_dp * delta**2 * y0p - 120.0_dp * delta * y0)&
        & + rCut**2 * (10.0_dp * delta**2 * y0pp + 60.0_dp * delta * y0p + 120.0_dp * y0)&
        & + delta**4 * y0pp + 8.0_dp * delta**3 * y0p + 20.0_dp * delta**2 * y0)&
        & / (2.0_dp * delta**5)

    dd = (rCut * (3.0_dp * delta**4 * y0pp + 24.0_dp * delta**3 * y0p + 60.0_dp * delta**2 * y0)&
        & + rCut**2 * (-12.0_dp * delta**3 * y0pp - 84.0_dp * delta**2 * y0p - 180.0_dp * delta&
        & * y0) + rCut**3 * (10.0_dp * delta**2 * y0pp + 60.0_dp * delta * y0p + 120.0_dp * y0))&
        & / (2.0_dp * delta**5)

    ee = -(rCut**2 * (3.0_dp * delta**4 * y0pp + 24.0_dp * delta**3 * y0p + 60.0_dp * delta**2&
        & * y0) + rCut**3*(-8.0_dp * delta**3 * y0pp - 56.0_dp * delta**2 * y0p - 120.0_dp * delta&
        & * y0) + rCut**4 * (5.0_dp * delta**2 * y0pp + 30.0_dp * delta * y0p + 60.0_dp * y0))&
        & / (2.0_dp * delta**5)

    ff = (rCut**3 * (delta**4 * y0pp + 8.0_dp * delta**3 * y0p + 20.0_dp * delta**2 * y0)&
        & + rCut**4 * (-2.0_dp * delta**3 * y0pp - 14.0_dp * delta**2 * y0p - 30.0_dp * delta * y0)&
        & + rCut**5 * (delta**2 * y0pp + 6.0_dp * delta * y0p + 12.0_dp * y0)) / (2.0_dp * delta**5)

    if (tPrime) then
      yy = 5.0_dp * aa * xx**4 + 4.0_dp * bb * xx**3 + 3.0_dp * cc * xx**2 + 2.0_dp * dd * xx + ee
    else
      yy = aa * xx**5 + bb * xx**4 + cc * xx**3 + dd * xx**2 + ee * xx + ff
    end if

  end function poly5zero


  !> Returns the numerical second derivative of long-range gamma, by means of a central finite
  !! difference.
  function getddLrNumericalGammaValue(this, iSp1, iSp2, dist, delta) result(ddGamma)

    !> Class instance
    class(THybridXcFunc), intent(in) :: this

    !> First species
    integer, intent(in) :: iSp1

    !> Second species
    integer, intent(in) :: iSp2

    !> Distance between atoms
    real(dp), intent(in) :: dist

    !> Delta for finite differences
    real(dp), intent(in) :: delta

    !> Numerical gamma derivative
    real(dp) :: ddGamma

    ddGamma = getddLrNumericalGammaValue_workhorse(this%hubbu(iSp1), this%hubbu(iSp2), this%omega,&
        & dist, delta)

  end function getddLrNumericalGammaValue


  !> Returns the numerical second derivative of full-range Hartree-Fock gamma, by means of a central
  !! finite difference.
  function getddHfNumericalGammaDeriv(this, iSp1, iSp2, dist, delta) result(ddGamma)

    !> Class instance
    class(THybridXcFunc), intent(in) :: this

    !> First species
    integer, intent(in) :: iSp1

    !> Second species
    integer, intent(in) :: iSp2

    !> Distance between atoms
    real(dp), intent(in) :: dist

    !> Delta for finite differences
    real(dp), intent(in) :: delta

    !> Numerical gamma derivative
    real(dp) :: ddGamma

    ddGamma = (&
        & getdHfAnalyticalGammaValue_workhorse(this%hubbu(iSp1), this%hubbu(iSp2), dist + delta)&
        & - getdHfAnalyticalGammaValue_workhorse(this%hubbu(iSp1), this%hubbu(iSp2), dist - delta))&
        & / (2.0_dp * delta)

  end function getddHfNumericalGammaDeriv


  !> Calculates analytical, truncated Coulomb, long-range gamma.
  function getLrTruncatedGammaValue(this, iSp1, iSp2, dist) result(gamma)

    !> Class instance
    class(THybridXcFunc_truncated), intent(in) :: this

    !> First species
    integer, intent(in) :: iSp1

    !> Second species
    integer, intent(in) :: iSp2

    !> Distance between atoms
    real(dp), intent(in) :: dist

    !> Resulting truncated gamma
    real(dp) :: gamma

    if (dist >= this%gammaCutoff) then
      gamma = 0.0_dp
    else
      gamma = getLrAnalyticalGammaValue_workhorse(this%hubbu(iSp1), this%hubbu(iSp2), this%omega,&
          & dist)
    end if

  end function getLrTruncatedGammaValue


  !> Calculates analytical, long-range gamma for MIC algorithm.
  function getLrMicGammaValue(this, iSp1, iSp2, dist) result(gamma)

    !> Class instance
    class(THybridXcFunc_mic), intent(in) :: this

    !> First species
    integer, intent(in) :: iSp1

    !> Second species
    integer, intent(in) :: iSp2

    !> Distance between atoms
    real(dp), intent(in) :: dist

    !> Resulting truncated gamma
    real(dp) :: gamma

    gamma = getLrAnalyticalGammaValue_workhorse(this%hubbu(iSp1), this%hubbu(iSp2), this%omega,&
        & dist)

  end function getLrMicGammaValue


  !> Calculates analytical, truncated and poly5zero damped Coulomb, long-range gamma.
  function getLrTruncatedAndDampedGammaValue(this, iSp1, iSp2, dist) result(gamma)

    !> Class instance
    class(THybridXcFunc_truncatedAndDamped), intent(in) :: this

    !> First species
    integer, intent(in) :: iSp1

    !> Second species
    integer, intent(in) :: iSp2

    !> Distance between atoms
    real(dp), intent(in) :: dist

    !> Resulting truncated gamma
    real(dp) :: gamma

    if (dist > this%gammaDamping .and. dist < this%gammaCutoff) then
      gamma = poly5zero(this%lrGammaAtDamping(iSp1, iSp2), this%lrdGammaAtDamping(iSp1, iSp2),&
          & this%lrddGammaAtDamping(iSp1, iSp2), dist, this%gammaDamping, this%gammaCutoff,&
          & tDerivative=.false.)
    elseif (dist >= this%gammaCutoff) then
      gamma = 0.0_dp
    else
      gamma = getLrAnalyticalGammaValue_workhorse(this%hubbu(iSp1), this%hubbu(iSp2), this%omega,&
          & dist)
    end if

  end function getLrTruncatedAndDampedGammaValue


  !> Calculates analytical, truncated Coulomb, full-range Hartree-Fock gamma.
  function getHfTruncatedGammaValue(this, iSp1, iSp2, dist) result(gamma)

    !> Class instance
    class(THybridXcFunc_truncated), intent(in) :: this

    !> First species
    integer, intent(in) :: iSp1

    !> Second species
    integer, intent(in) :: iSp2

    !> Distance between atoms
    real(dp), intent(in) :: dist

    !> Resulting truncated gamma
    real(dp) :: gamma

    if (dist >= this%gammaCutoff) then
      gamma = 0.0_dp
    else
      gamma = getHfAnalyticalGammaValue_workhorse(this%hubbu(iSp1), this%hubbu(iSp2), dist)
    end if

  end function getHfTruncatedGammaValue


  !> Calculates analytical, full-range Hartree-Fock gamma for MIC algorithm.
  function getHfMicGammaValue(this, iSp1, iSp2, dist) result(gamma)

    !> Class instance
    class(THybridXcFunc_mic), intent(in) :: this

    !> First species
    integer, intent(in) :: iSp1

    !> Second species
    integer, intent(in) :: iSp2

    !> Distance between atoms
    real(dp), intent(in) :: dist

    !> Resulting truncated gamma
    real(dp) :: gamma

    gamma = getHfAnalyticalGammaValue_workhorse(this%hubbu(iSp1), this%hubbu(iSp2), dist)

  end function getHfMicGammaValue


  !> Calculates analytical, truncated and poly5zero damped Coulomb, full-range Hartree-Fock gamma.
  function getHfTruncatedAndDampedGammaValue(this, iSp1, iSp2, dist) result(gamma)

    !> Class instance
    class(THybridXcFunc_truncatedAndDamped), intent(in) :: this

    !> First species
    integer, intent(in) :: iSp1

    !> Second species
    integer, intent(in) :: iSp2

    !> Distance between atoms
    real(dp), intent(in) :: dist

    !> Resulting truncated gamma
    real(dp) :: gamma

    if (dist > this%gammaDamping .and. dist < this%gammaCutoff) then
      gamma = poly5zero(this%hfGammaAtDamping(iSp1, iSp2), this%hfdGammaAtDamping(iSp1, iSp2),&
          & this%hfddGammaAtDamping(iSp1, iSp2), dist, this%gammaDamping, this%gammaCutoff,&
          & tDerivative=.false.)
    elseif (dist >= this%gammaCutoff) then
      gamma = 0.0_dp
    else
      gamma = getHfAnalyticalGammaValue_workhorse(this%hubbu(iSp1), this%hubbu(iSp2), dist)
    end if

  end function getHfTruncatedAndDampedGammaValue


  !> Returns g-sum of long-range and full-range Hartree-Fock gamma's.
  function getCamGammaGSum(this, iAt1, iAt2, iSp1, iSp2, rCellVecsG) result(gamma)

    !> Class instance
    class(THybridXcFunc), intent(in) :: this

    !> Index of first and second atom
    integer, intent(in) :: iAt1, iAt2

    !> Index of first and second species
    integer, intent(in) :: iSp1, iSp2

    !> Vectors to unit cells in absolute units
    real(dp), intent(in) :: rCellVecsG(:,:)

    !> Resulting Gamma, summed up for g-vectors
    real(dp) :: gamma

    !! Distance between the two atoms
    real(dp) :: dist

    !! Index of real-space \vec{g} summation
    integer :: iG

    gamma = 0.0_dp

    if (this%hybridXcType == hybridXcFunc%lc .or. this%hybridXcType == hybridXcFunc%cam) then
      loopGLr: do iG = 1, size(rCellVecsG, dim=2)
        dist = norm2(this%rCoords(:, iAt1) - (this%rCoords(:, iAt2) + rCellVecsG(:, iG)))
        gamma = gamma + this%camBeta * this%getLrGammaValue(iSp1, iSp2, dist)
      end do loopGLr
    end if

    if (this%hybridXcType == hybridXcFunc%hyb .or. this%hybridXcType == hybridXcFunc%cam) then
      loopGHf: do iG = 1, size(rCellVecsG, dim=2)
        dist = norm2(this%rCoords(:, iAt1) - (this%rCoords(:, iAt2) + rCellVecsG(:, iG)))
        gamma = gamma + this%camAlpha * this%getHfGammaValue(iSp1, iSp2, dist)
      end do loopGHf
    end if

  end function getCamGammaGSum


  !> Calculates the derivative of the square, dense, unpacked, dual-space overlap matrix w.r.t.
  !! given atom.
  subroutine getUnpackedOverlapPrime_kpts(iAtomPrime, skOverCont, orb, derivator, symNeighbourList,&
      & nNeighbourCamSym, iSquare, cellVec, rCoords, kPoint, overSqrPrime)

    !> Overlap derivative calculated w.r.t. this atom in the central cell
    integer, intent(in) :: iAtomPrime

    !> Sparse overlap container
    type(TSlakoCont), intent(in) :: skOverCont

    !> Orbital information.
    type(TOrbitals), intent(in) :: orb

    !> Differentiation object
    class(TNonSccDiff), intent(in) :: derivator

    !> List of neighbours for each atom (symmetric version)
    type(TSymNeighbourList), intent(in) :: symNeighbourList

    !> Nr. of neighbours for each atom.
    integer, intent(in) :: nNeighbourCamSym(:)

    !> Position of each atom in the rows/columns of the square matrices. Shape: (nAtom)
    integer, intent(in) :: iSquare(:)

    !> Vectors to neighboring unit cells in relative units
    real(dp), intent(in) :: cellVec(:,:)

    !> Atomic coordinates in absolute units, potentially including periodic images
    real(dp), intent(in) :: rCoords(:,:)

    !> Relative coordinates of the k-point where the overlap should be constructed
    real(dp), intent(in) :: kPoint(:)

    !> Overlap derivative
    complex(dp), intent(out) :: overSqrPrime(:,:,:)

    !! Temporary overlap derivative storage
    real(dp) :: overPrime(orb%mOrb, orb%mOrb, 3)

    !! Dense matrix descriptor indices
    integer, parameter :: descLen = 3, iStart = 1, iEnd = 2, iNOrb = 3

    !! Stores start/end index and number of orbitals of square matrices
    integer :: descAt1(descLen), descAt2(descLen)

    !! Atom to calculate energy gradient components for
    integer :: iAt2, iAt2fold, iNeigh

    !! Iterates over coordinates
    integer :: iCoord

    !! The k-point in relative coordinates, multiplied by 2pi
    real(dp) :: kPoint2p(3)

    !! Phase factor
    complex(dp) :: phase

    !! Auxiliary variables
    integer :: iVec

    overSqrPrime(:,:,:) = (0.0_dp, 0.0_dp)

    kPoint2p(:) = 2.0_dp * pi * kPoint

    descAt1 = getDescriptor(iAtomPrime, iSquare)
    do iNeigh = 0, nNeighbourCamSym(iAtomPrime)
      iAt2 = symNeighbourList%neighbourList%iNeighbour(iNeigh, iAtomPrime)
      iAt2fold = symNeighbourList%img2CentCell(iAt2)
      if (iAtomPrime == iAt2fold) cycle
      descAt2 = getDescriptor(iAt2fold, iSquare)
      iVec = symNeighbourList%iCellVec(iAt2)
      overPrime(:,:,:) = 0.0_dp
      call derivator%getFirstDeriv(overPrime, skOverCont, rCoords, symNeighbourList%species,&
          & iAtomPrime, iAt2, orb)
      phase = exp(imag * dot_product(kPoint2p, cellVec(:, iVec)))
      do iCoord = 1, 3
        overSqrPrime(iCoord, descAt2(iStart):descAt2(iStart) + descAt2(iNOrb) - 1,&
            & descAt1(iStart):descAt1(iStart) + descAt1(iNOrb) - 1)&
            & = overSqrPrime(iCoord, descAt2(iStart):descAt2(iStart) + descAt2(iNOrb) - 1,&
            & descAt1(iStart):descAt1(iStart) + descAt1(iNOrb) - 1)&
            & + overPrime(1:descAt2(iNOrb), 1:descAt1(iNOrb), iCoord) * phase
      end do
    end do

  end subroutine getUnpackedOverlapPrime_kpts


  !> Calculates pseudo Fourier transform of square CAM y-matrix with shape [nOrb, nOrb].
  subroutine getCamGammaFourierAO(this, iSquare, cellVecsG, rCellVecsG, kPoint, kPointPrime,&
      & camGammaSqr)

    !> Class instance
    class(THybridXcFunc), intent(in) :: this

    !> Position of each atom in the rows/columns of the square matrices. Shape: (nAtom)
    integer, intent(in) :: iSquare(:)

    !> Translation vectors to lattice cells in units of lattice constants
    real(dp), intent(in) :: cellVecsG(:,:)

    !> Vectors to unit cells in absolute units
    real(dp), intent(in) :: rCellVecsG(:,:)

    !> First and second k-point in relative coordinates
    real(dp), intent(in) :: kPoint(:), kPointPrime(:)

    !> Pseudo Fourier transform of square gamma matrix
    complex(dp), intent(out) :: camGammaSqr(:,:)

    !! Indices of shifted atom
    integer :: iAt1, iSp1

    !! Indices of atom in central cell
    integer :: iAt2, iSp2

    !! Number of atoms in central cell
    integer :: nAtom0

    !! Auxiliary variables
    integer :: iOrbStart1, iOrbEnd1, iOrbStart2, iOrbEnd2

    nAtom0 = size(this%species0)

    camGammaSqr(:,:) = (0.0_dp, 0.0_dp)

    do iAt2 = 1, nAtom0
      iSp2 = this%species0(iAt2)
      iOrbStart2 = iSquare(iAt2)
      iOrbEnd2 = iSquare(iAt2 + 1) - 1
      do iAt1 = 1, nAtom0
        iSp1 = this%species0(iAt1)
        iOrbStart1 = iSquare(iAt1)
        iOrbEnd1 = iSquare(iAt1 + 1) - 1
        camGammaSqr(iOrbStart1:iOrbEnd1, iOrbStart2:iOrbEnd2)&
            & = getCamGammaFourier(this, iAt1, iAt2, iSp1, iSp2, cellVecsG, rCellVecsG, kPoint,&
            & kPointPrime)
      end do
    end do

  end subroutine getCamGammaFourierAO


  !> Returns pseudo Fourier transform of long-range and full-range Hartree-Fock gammas.
  function getCamGammaFourier(this, iAt1, iAt2, iSp1, iSp2, cellVecsG, rCellVecsG, kPoint,&
      & kPointPrime) result(gamma)

    !> Class instance
    class(THybridXcFunc), intent(in) :: this

    !> Index of first and second atom
    integer, intent(in) :: iAt1, iAt2

    !> Index of first and second species
    integer, intent(in) :: iSp1, iSp2

    !> Translation vectors to lattice cells in units of lattice constants
    real(dp), intent(in) :: cellVecsG(:,:)

    !> Vectors to unit cells in absolute units
    real(dp), intent(in) :: rCellVecsG(:,:)

    !> First and second k-point in relative coordinates
    real(dp), intent(in) :: kPoint(:), kPointPrime(:)

    !> Resulting Gamma, summed up for g-vectors
    complex(dp) :: gamma

    !! Phase factor
    complex(dp) :: phase

    !! Distance between the two atoms
    real(dp) :: dist

    !! Difference (k - k') with both k-points in relative coordinates, multiplied by 2pi
    real(dp) :: kPoint2p(3)

    !! Index of real-space \vec{g} summation
    integer :: iG

    kPoint2p(:) = 2.0_dp * pi * (kPoint - kPointPrime)

    gamma = (0.0_dp, 0.0_dp)

    if (this%hybridXcType == hybridXcFunc%lc .or. this%hybridXcType == hybridXcFunc%cam) then
      loopGLr: do iG = 1, size(rCellVecsG, dim=2)
        dist = norm2(this%rCoords(:, iAt2) - (this%rCoords(:, iAt1) + rCellVecsG(:, iG)))
        phase = exp(imag * dot_product(kPoint2p, cellVecsG(:, iG)))
        gamma = gamma + this%camBeta * this%getLrGammaValue(iSp1, iSp2, dist) * phase
      end do loopGLr
    end if

    if (this%hybridXcType == hybridXcFunc%hyb .or. this%hybridXcType == hybridXcFunc%cam) then
      loopGHf: do iG = 1, size(rCellVecsG, dim=2)
        dist = norm2(this%rCoords(:, iAt2) - (this%rCoords(:, iAt1) + rCellVecsG(:, iG)))
        phase = exp(imag * dot_product(kPoint2p, cellVecsG(:, iG)))
        gamma = gamma + this%camAlpha * this%getHfGammaValue(iSp1, iSp2, dist) * phase
      end do loopGHf
    end if

  end function getCamGammaFourier


  !> Calculates derivatives w.r.t. given atom of the pseudo Fourier transformed square CAM y-matrix
  !! with shape [nOrb, nOrb].
  subroutine getCamGammaFourierAOPrime(this, iAtomPrime, iSquare, cellVecsG, rCellVecsG, kPoint,&
      & kPointPrime, camdGammaAO)

    !> Class instance
    class(THybridXcFunc), intent(in) :: this

    !> Derivative calculated w.r.t. this atom in the central cell
    integer, intent(in) :: iAtomPrime

    !> Position of each atom in the rows/columns of the square matrices. Shape: (nAtom)
    integer, intent(in) :: iSquare(:)

    !> Translation vectors to lattice cells in units of lattice constants
    real(dp), intent(in) :: cellVecsG(:,:)

    !> Vectors to unit cells in absolute units
    real(dp), intent(in) :: rCellVecsG(:,:)

    !> First and second k-point in relative coordinates
    real(dp), intent(in) :: kPoint(:), kPointPrime(:)

    !> Pseudo Fourier transform of square gamma matrix
    complex(dp), intent(out) :: camdGammaAO(:,:,:)

    !! Indices of atom to calculate derivative for
    integer :: iSpAtPrime, iOrbStartAtPrime, iOrbEndAtPrime

    !! Indices of atom in central cell
    integer :: iAt2, iSp2

    !! Number of atoms in central cell
    integer :: nAtom0

    !! Temporary storage for Gamma derivative
    complex(dp) :: dGamma(3)

    !! Iterates over coordinates
    integer :: iCoord

    !! Auxiliary variables
    integer :: iOrbStart2, iOrbEnd2

    nAtom0 = size(this%species0)

    camdGammaAO(:,:,:) = (0.0_dp, 0.0_dp)

    iSpAtPrime = this%species0(iAtomPrime)
    iOrbStartAtPrime = iSquare(iAtomPrime)
    iOrbEndAtPrime = iSquare(iAtomPrime + 1) - 1

    do iAt2 = 1, nAtom0
      if (iAtomPrime == iAt2) cycle
      iSp2 = this%species0(iAt2)
      iOrbStart2 = iSquare(iAt2)
      iOrbEnd2 = iSquare(iAt2 + 1) - 1
      dGamma(:) = getCamGammaFourierPrime(this, iAtomPrime, iAt2, iSpAtPrime, iSp2, cellVecsG,&
          & rCellVecsG, kPoint, kPointPrime)
      do iCoord = 1, 3
        camdGammaAO(iOrbStartAtPrime:iOrbEndAtPrime, iOrbStart2:iOrbEnd2, iCoord) = dGamma(iCoord)
      end do
    end do

  end subroutine getCamGammaFourierAOPrime


  !> Calculates derivative w.r.t. given atom of the pseudo Fourier transformed square CAM y-matrix.
  function getCamGammaFourierPrime(this, iAt1, iAt2, iSp1, iSp2, cellVecsG, rCellVecsG, kPoint,&
      & kPointPrime) result(dGamma)

    !> Class instance
    class(THybridXcFunc), intent(in) :: this

    !> Index of first and second atom
    integer, intent(in) :: iAt1, iAt2

    !> Index of first and second species
    integer, intent(in) :: iSp1, iSp2

    !> Translation vectors to lattice cells in units of lattice constants
    real(dp), intent(in) :: cellVecsG(:,:)

    !> Vectors to unit cells in absolute units
    real(dp), intent(in) :: rCellVecsG(:,:)

    !> First and second k-point in relative coordinates
    real(dp), intent(in) :: kPoint(:), kPointPrime(:)

    !> Resulting Gamma derivative, summed up for g-vectors
    complex(dp) :: dGamma(3)

    !! Phase factor
    complex(dp) :: phase

    !! Temporary distance vector
    real(dp) :: distVect(3)

    !! Distance between the two atoms
    real(dp) :: dist

    !! Difference (k - k') with both k-points in relative coordinates, multiplied by 2pi
    real(dp) :: kPoint2p(3)

    !! Index of real-space \vec{g} summation
    integer :: iG

    kPoint2p(:) = 2.0_dp * pi * (kPoint - kPointPrime)

    dGamma(:) = (0.0_dp, 0.0_dp)

    if (this%hybridXcType == hybridXcFunc%lc .or. this%hybridXcType == hybridXcFunc%cam) then
      loopGLr: do iG = 1, size(rCellVecsG, dim=2)
        distVect(:) = this%rCoords(:, iAt1) - (this%rCoords(:, iAt2) + rCellVecsG(:, iG))
        dist = norm2(distVect)
        distVect(:) = distVect / dist
        phase = exp(imag * dot_product(kPoint2p, cellVecsG(:, iG)))
        dGamma(:) = dGamma&
            & + distVect * this%camBeta * this%getLrGammaPrimeValue(iSp1, iSp2, dist) * phase
      end do loopGLr
    end if

    if (this%hybridXcType == hybridXcFunc%hyb .or. this%hybridXcType == hybridXcFunc%cam) then
      loopGHf: do iG = 1, size(rCellVecsG, dim=2)
        distVect(:) = this%rCoords(:, iAt1) - (this%rCoords(:, iAt2) + rCellVecsG(:, iG))
        dist = norm2(distVect)
        distVect(:) = distVect / dist
        phase = exp(imag * dot_product(kPoint2p, cellVecsG(:, iG)))
        dGamma(:) = dGamma&
            & + distVect * this%camAlpha * this%getHfGammaPrimeValue(iSp1, iSp2, dist) * phase
      end do loopGHf
    end if

  end function getCamGammaFourierPrime


  !> Returns g-resolved, long-range + full-range Hartree-Fock gamma.
  function getCamGammaGResolved(this, iAt1, iAt2, iSp1, iSp2, rCellVecsG, rShift)&
      & result(gammas)

    !> Class instance
    class(THybridXcFunc), intent(in) :: this

    !> Index of first and second atom
    integer, intent(in) :: iAt1, iAt2

    !> Index of first and second species
    integer, intent(in) :: iSp1, iSp2

    !> Vectors to unit cells in absolute units
    real(dp), intent(in) :: rCellVecsG(:,:)

    !> Absolute shift of atom position (Gamma arguments)
    real(dp), intent(in) :: rShift(:)

    !> Resulting Gammas, g-vector resolved
    real(dp) :: gammas(size(rCellVecsG, dim=2))

    !! Total shift of atom position (Gamma arguments) in absolute coordinates
    real(dp) :: rTotshift(3)

    !! Distance between the two atoms
    real(dp) :: dist

    !! Index of real-space \vec{g} summation
    integer :: iG

    gammas(:) = 0.0_dp

    if (this%hybridXcType == hybridXcFunc%lc .or. this%hybridXcType == hybridXcFunc%cam) then
      loopGLr: do iG = 1, size(rCellVecsG, dim=2)
        rTotshift(:) = rCellVecsG(:, iG) + rShift
        dist = norm2(this%rCoords(:, iAt1) - (this%rCoords(:, iAt2) + rTotshift))
        gammas(iG) = gammas(iG) + this%camBeta * this%getLrGammaValue(iSp1, iSp2, dist)
      end do loopGLr
    end if

    if (this%hybridXcType == hybridXcFunc%hyb .or. this%hybridXcType == hybridXcFunc%cam) then
      loopGHf: do iG = 1, size(rCellVecsG, dim=2)
        rTotshift(:) = rCellVecsG(:, iG) + rShift
        dist = norm2(this%rCoords(:, iAt1) - (this%rCoords(:, iAt2) + rTotshift))
        gammas(iG) = gammas(iG) + this%camAlpha * this%getHfGammaValue(iSp1, iSp2, dist)
      end do loopGHf
    end if

  end function getCamGammaGResolved


  !> Calculates analytical, unaltered long-range gamma derivative.
  function getdLrAnalyticalGammaValue(this, iSp1, iSp2, dist) result(dGamma)

    !> Class instance
    class(THybridXcFunc_full), intent(in) :: this

    !> First species
    integer, intent(in) :: iSp1

    !> Second species
    integer, intent(in) :: iSp2

    !> Distance between atoms
    real(dp), intent(in) :: dist

    !> Resulting gamma derivative
    real(dp) :: dGamma

    dGamma = getdLrAnalyticalGammaValue_workhorse(this%hubbu(iSp1), this%hubbu(iSp2), this%omega,&
        & dist)

  end function getdLrAnalyticalGammaValue


  !> Calculates analytical, unaltered full-range Hartree-Fock gamma derivative.
  function getdHfAnalyticalGammaValue(this, iSp1, iSp2, dist) result(dGamma)

    !> Class instance
    class(THybridXcFunc_full), intent(in) :: this

    !> First species
    integer, intent(in) :: iSp1

    !> Second species
    integer, intent(in) :: iSp2

    !> Distance between atoms
    real(dp), intent(in) :: dist

    !> Resulting gamma derivative
    real(dp) :: dGamma

    dGamma = getdHfAnalyticalGammaValue_workhorse(this%hubbu(iSp1), this%hubbu(iSp2), dist)

  end function getdHfAnalyticalGammaValue


  !> Calculates analytical, truncated Coulomb, long-range gamma derivative.
  function getdLrTruncatedGammaValue(this, iSp1, iSp2, dist) result(dGamma)

    !> Class instance
    class(THybridXcFunc_truncated), intent(in) :: this

    !> First species
    integer, intent(in) :: iSp1

    !> Second species
    integer, intent(in) :: iSp2

    !> Distance between atoms
    real(dp), intent(in) :: dist

    !> Resulting truncated gamma derivative
    real(dp) :: dGamma

    if (dist >= this%gammaCutoff) then
      dGamma = 0.0_dp
    else
      dGamma = getdLrAnalyticalGammaValue_workhorse(this%hubbu(iSp1), this%hubbu(iSp2), this%omega,&
          & dist)
    end if

  end function getdLrTruncatedGammaValue


  !> Calculates analytical, long-range gamma derivative for MIC algorithm.
  function getdLrMicGammaValue(this, iSp1, iSp2, dist) result(dGamma)

    !> Class instance
    class(THybridXcFunc_mic), intent(in) :: this

    !> First species
    integer, intent(in) :: iSp1

    !> Second species
    integer, intent(in) :: iSp2

    !> Distance between atoms
    real(dp), intent(in) :: dist

    !> Resulting truncated gamma derivative
    real(dp) :: dGamma

    dGamma = getdLrAnalyticalGammaValue_workhorse(this%hubbu(iSp1), this%hubbu(iSp2), this%omega,&
        & dist)

  end function getdLrMicGammaValue


  !> Calculates analytical, truncated Coulomb, full-range Hartree-Fock gamma derivative.
  function getdHfTruncatedGammaValue(this, iSp1, iSp2, dist) result(dGamma)

    !> Class instance
    class(THybridXcFunc_truncated), intent(in) :: this

    !> First species
    integer, intent(in) :: iSp1

    !> Second species
    integer, intent(in) :: iSp2

    !> Distance between atoms
    real(dp), intent(in) :: dist

    !> Resulting truncated gamma derivative
    real(dp) :: dGamma

    if (dist >= this%gammaCutoff) then
      dGamma = 0.0_dp
    else
      dGamma = getdHfAnalyticalGammaValue_workhorse(this%hubbu(iSp1), this%hubbu(iSp2), dist)
    end if

  end function getdHfTruncatedGammaValue


  !> Calculates analytical, full-range Hartree-Fock gamma derivative for MIC algorithm.
  function getdHfMicGammaValue(this, iSp1, iSp2, dist) result(dGamma)

    !> Class instance
    class(THybridXcFunc_mic), intent(in) :: this

    !> First species
    integer, intent(in) :: iSp1

    !> Second species
    integer, intent(in) :: iSp2

    !> Distance between atoms
    real(dp), intent(in) :: dist

    !> Resulting truncated gamma derivative
    real(dp) :: dGamma

    dGamma = getdHfAnalyticalGammaValue_workhorse(this%hubbu(iSp1), this%hubbu(iSp2), dist)

  end function getdHfMicGammaValue


  !> Calculates analytical, truncated and poly5zero damped Coulomb, long-range gamma derivative.
  function getdLrTruncatedAndDampedGammaValue(this, iSp1, iSp2, dist) result(dGamma)

    !> Class instance
    class(THybridXcFunc_truncatedAndDamped), intent(in) :: this

    !> First species
    integer, intent(in) :: iSp1

    !> Second species
    integer, intent(in) :: iSp2

    !> Distance between atoms
    real(dp), intent(in) :: dist

    !> Resulting truncated gamma derivative
    real(dp) :: dGamma

    if (dist > this%gammaDamping .and. dist < this%gammaCutoff) then
      dGamma = poly5zero(this%lrGammaAtDamping(iSp1, iSp2), this%lrdGammaAtDamping(iSp1, iSp2),&
          & this%lrddGammaAtDamping(iSp1, iSp2), dist, this%gammaDamping, this%gammaCutoff,&
          & tDerivative=.true.)
    elseif (dist >= this%gammaCutoff) then
      dGamma = 0.0_dp
    else
      dGamma = getdLrAnalyticalGammaValue_workhorse(this%hubbu(iSp1), this%hubbu(iSp2), this%omega,&
          & dist)
    end if

  end function getdLrTruncatedAndDampedGammaValue


  !> Calculates analytical, truncated and poly5zero damped Coulomb, full-range Hartree-Fock gamma
  !! derivative.
  function getdHfTruncatedAndDampedGammaValue(this, iSp1, iSp2, dist) result(dGamma)

    !> Class instance
    class(THybridXcFunc_truncatedAndDamped), intent(in) :: this

    !> First species
    integer, intent(in) :: iSp1

    !> Second species
    integer, intent(in) :: iSp2

    !> Distance between atoms
    real(dp), intent(in) :: dist

    !> Resulting truncated gamma derivative (1st)
    real(dp) :: dGamma

    if (dist > this%gammaDamping .and. dist < this%gammaCutoff) then
      dGamma = poly5zero(this%lrGammaAtDamping(iSp1, iSp2), this%lrdGammaAtDamping(iSp1, iSp2),&
          & this%lrddGammaAtDamping(iSp1, iSp2), dist, this%gammaDamping, this%gammaCutoff,&
          & tDerivative=.true.)
    elseif (dist >= this%gammaCutoff) then
      dGamma = 0.0_dp
    else
      dGamma = getdHfAnalyticalGammaValue_workhorse(this%hubbu(iSp1), this%hubbu(iSp2), dist)
    end if

  end function getdHfTruncatedAndDampedGammaValue


  !> Returns g-sum of long-range and full-range Hartree-Fock gamma derivatives.
  function getCamGammaPrimeGSum(this, iAt1, iAt2, iSp1, iSp2, rCellVecsG) result(dGamma)

    !> Class instance
    class(THybridXcFunc), intent(in) :: this

    !> Index of first and second atom
    integer, intent(in) :: iAt1, iAt2

    !> Index of first and second species
    integer, intent(in) :: iSp1, iSp2

    !> Vectors to unit cells in absolute units
    real(dp), intent(in) :: rCellVecsG(:,:)

    !> Resulting Gamma derivative (1st), summed up for g-vectors
    real(dp) :: dGamma(3)

    !! Temporary distance vector
    real(dp) :: distVect(3)

    !! Distance between the two atoms
    real(dp) :: dist

    !! Index of real-space \vec{g} summation
    integer :: iG

    dGamma(:) = 0.0_dp

    if (this%hybridXcType == hybridXcFunc%lc .or. this%hybridXcType == hybridXcFunc%cam) then
      loopGLr: do iG = 1, size(rCellVecsG, dim=2)
        distVect(:) = this%rCoords(:, iAt1) - (this%rCoords(:, iAt2) + rCellVecsG(:, iG))
        dist = norm2(distVect)
        distVect(:) = distVect / dist
        dGamma(:) = dGamma + distVect * this%camBeta * this%getLrGammaPrimeValue(iSp1, iSp2, dist)
      end do loopGLr
    end if

    if (this%hybridXcType == hybridXcFunc%hyb .or. this%hybridXcType == hybridXcFunc%cam) then
      loopGHf: do iG = 1, size(rCellVecsG, dim=2)
        distVect(:) = this%rCoords(:, iAt1) - (this%rCoords(:, iAt2) + rCellVecsG(:, iG))
        dist = norm2(distVect)
        distVect(:) = distVect / dist
        dGamma(:) = dGamma + distVect * this%camAlpha * this%getHfGammaPrimeValue(iSp1, iSp2, dist)
      end do loopGHf
    end if

  end function getCamGammaPrimeGSum


  !> Returns g-resolved long-range and full-range Hartree-Fock gamma derivatives.
  function getCamGammaPrimeGResolved(this, iAt1, iAt2, iSp1, iSp2, rCellVecsG) result(dGammas)

    !> Class instance
    class(THybridXcFunc), intent(in) :: this

    !> Index of first and second atom
    integer, intent(in) :: iAt1, iAt2

    !> Index of first and second species
    integer, intent(in) :: iSp1, iSp2

    !> Vectors to unit cells in absolute units
    real(dp), intent(in) :: rCellVecsG(:,:)

    !> Resulting Gamma derivative (1st), for each g-vector
    real(dp) :: dGammas(3, size(rCellVecsG, dim=2))

    !! Temporary distance vector
    real(dp) :: distVect(3)

    !! Distance between the two atoms
    real(dp) :: dist

    !! Index of real-space \vec{g} summation
    integer :: iG

    dGammas(:,:) = 0.0_dp

    if (this%hybridXcType == hybridXcFunc%lc .or. this%hybridXcType == hybridXcFunc%cam) then
      loopGLr: do iG = 1, size(rCellVecsG, dim=2)
        distVect(:) = this%rCoords(:, iAt1) - (this%rCoords(:, iAt2) + rCellVecsG(:, iG))
        dist = norm2(distVect)
        distVect(:) = distVect / dist
        dGammas(:, iG) = dGammas(:, iG)&
            & + distVect * this%camBeta * this%getLrGammaPrimeValue(iSp1, iSp2, dist)
      end do loopGLr
    end if

    if (this%hybridXcType == hybridXcFunc%hyb .or. this%hybridXcType == hybridXcFunc%cam) then
      loopGHf: do iG = 1, size(rCellVecsG, dim=2)
        distVect(:) = this%rCoords(:, iAt1) - (this%rCoords(:, iAt2) + rCellVecsG(:, iG))
        dist = norm2(distVect)
        distVect(:) = distVect / dist
        dGammas(:, iG) = dGammas(:, iG)&
            & + distVect * this%camAlpha * this%getHfGammaPrimeValue(iSp1, iSp2, dist)
      end do loopGHf
    end if

  end function getCamGammaPrimeGResolved


  !> Calculates analytical long-range gamma.
  function getLrAnalyticalGammaValue(this, iSp1, iSp2, dist) result(gamma)

    !> Class instance
    class(THybridXcFunc_full), intent(in) :: this

    !> First species
    integer, intent(in) :: iSp1

    !> Second species
    integer, intent(in) :: iSp2

    !> Distance between atoms
    real(dp), intent(in) :: dist

    !> Resulting gamma
    real(dp) :: gamma

    gamma = getLrAnalyticalGammaValue_workhorse(this%hubbu(iSp1), this%hubbu(iSp2), this%omega,&
        & dist)

  end function getLrAnalyticalGammaValue


  !> Returns the derivative of CAM range-separated gamma for iAtom1, iAtom2 (non-periodic version).
  subroutine getDirectionalCamGammaPrimeValue_cluster(this, grad, iAtom1, iAtom2)

    !> Class instance
    class(THybridXcFunc), intent(in) :: this

    !> Gradient of gamma between atoms
    real(dp), intent(out) :: grad(:)

    !> First atom
    integer, intent(in) :: iAtom1

    !> Second atom
    integer, intent(in) :: iAtom2

    !! Species index of first and second atom
    integer :: iSp1, iSp2

    !! Distance(-vector) of the two atoms
    real(dp) :: vect(3), dist

    grad(:) = 0.0_dp

    iSp1 = this%species0(iAtom1)
    iSp2 = this%species0(iAtom2)

    ! analytical derivatives
    vect(:) = this%rCoords(:, iAtom1) - this%rCoords(:, iAtom2)
    dist = sqrt(sum(vect**2))
    vect(:) = vect / dist

    if (this%hybridXcType == hybridXcFunc%lc .or. this%hybridXcType == hybridXcFunc%cam) then
      grad(:) = grad + vect * this%camBeta * this%getLrGammaPrimeValue(iSp1, iSp2, dist)
    end if

    if (this%hybridXcType == hybridXcFunc%hyb .or. this%hybridXcType == hybridXcFunc%cam) then
      grad(:) = grad + vect * this%camAlpha * this%getHfGammaPrimeValue(iSp1, iSp2, dist)
    end if

  end subroutine getDirectionalCamGammaPrimeValue_cluster


  !> Returns the derivative of CAM range-separated gamma for iAtom1, iAtom2 (periodic version).
  subroutine getDirectionalCamGammaPrimeValue_periodic(this, grad, iAtom1, iAtom2, img2CentCell)

    !> Class instance
    class(THybridXcFunc), intent(in) :: this

    !> Gradient of gamma between atoms
    real(dp), intent(out) :: grad(:)

    !> First atom
    integer, intent(in) :: iAtom1

    !> Second atom
    integer, intent(in) :: iAtom2

    !> Map images of atoms to the central cell
    integer, intent(in) :: img2CentCell(:)

    !! Species index of first and second atom
    integer :: iSp1, iSp2

    !! Distance(-vector) of the two atoms
    real(dp) :: vect(3), dist

    iSp1 = this%species0(img2CentCell(iAtom1))
    iSp2 = this%species0(img2CentCell(iAtom2))

    ! analytical derivatives
    vect(:) = this%rCoords(:, iAtom1) - this%rCoords(:, iAtom2)
    dist = sqrt(sum(vect**2))
    vect(:) = vect / dist

    if (this%hybridXcType == hybridXcFunc%lc .or. this%hybridXcType == hybridXcFunc%cam) then
      grad(:) = grad + vect * this%camBeta * this%getLrGammaPrimeValue(iSp1, iSp2, dist)
    end if

    if (this%hybridXcType == hybridXcFunc%hyb .or. this%hybridXcType == hybridXcFunc%cam) then
      grad(:) = grad + vect * this%camAlpha * this%getHfGammaPrimeValue(iSp1, iSp2, dist)
    end if

  end subroutine getDirectionalCamGammaPrimeValue_periodic


#:if WITH_SCALAPACK

  !> Interface routine to add gradients due to CAM range-separated contributions.
  !! (non-periodic and Gamma-only version)
  subroutine addCamGradients_real(this, env, parallelKS, deltaRhoSqr, SSqrReal, skOverCont,&
      & symNeighbourList, nNeighbourCamSym, orb, derivator, denseDesc, nSpin, tPeriodic, gradients,&
      & errStatus)

    !> Class instance
    class(THybridXcFunc), intent(inout) :: this

    !> Environment settings
    type(TEnvironment), intent(in) :: env

    !> The k-points and spins to process
    type(TParallelKS), intent(in) :: parallelKS

    !> Square (unpacked) delta density matrix
    real(dp), intent(in) :: deltaRhoSqr(:,:,:)

    !> Square (unpacked) overlap matrix
    real(dp), intent(in) :: SSqrReal(:,:)

    !> Sparse overlap part
    type(TSlakoCont), intent(in) :: skOverCont

    !> List of neighbours for each atom (symmetric version)
    type(TSymNeighbourList), intent(in) :: symNeighbourList

    !> Nr. of neighbours for each atom
    integer, intent(in) :: nNeighbourCamSym(:)

    !> Orbital information for system
    type(TOrbitals), intent(in) :: orb

    !> Differentiation object
    class(TNonSccDiff), intent(in) :: derivator

    !> Dense matrix descriptor
    type(TDenseDescr), intent(in) :: denseDesc

    !> Number of spin channels in calculation
    integer, intent(in) :: nSpin

    !> True, if system is periodic (i.e. Gamma-only)
    logical, intent(in) :: tPeriodic

    !> Energy gradients
    real(dp), intent(inout) :: gradients(:,:)

    !> Error status
    type(TStatus), intent(inout) :: errStatus

    if (tPeriodic) then
      call this%tabulateCamdGammaEval0_gamma()
    else
      call this%tabulateCamdGammaEval0_cluster()
    end if

    select case(this%hybridXcAlg)
    case (hybridXcAlgo%thresholdBased, hybridXcAlgo%neighbourBased)
      @:RAISE_ERROR(errStatus, -1, "HybridXc Module: MPI parallelized force evaluation not&
          & supported, choose matrix-based algorithm instead.")
    case (hybridXcAlgo%matrixBased)
      call addCamGradientsMatrix_real_blacs(this, env, parallelKS, deltaRhoSqr, SSqrReal,&
          & skOverCont, symNeighbourList, nNeighbourCamSym, orb, derivator, denseDesc, nSpin,&
          & gradients)
    end select

  end subroutine addCamGradients_real

#:else

  !> Interface routine to add gradients due to CAM range-separated contributions.
  !! (non-periodic and Gamma-only version)
  subroutine addCamGradients_real(this, deltaRhoSqr, SSqrReal, skOverCont, orb, iSquare,&
      & iNeighbour, nNeighbourSK, derivator, tPeriodic, gradients, symNeighbourList,&
      & nNeighbourCamSym)

    !> Class instance
    class(THybridXcFunc), intent(inout) :: this

    !> Square (unpacked) delta density matrix
    real(dp), intent(in) :: deltaRhoSqr(:,:,:)

    !> Square (unpacked) overlap matrix
    real(dp), intent(in) :: SSqrReal(:,:)

    !> Sparse overlap part
    type(TSlakoCont), intent(in) :: skOverCont

    !> Orbital information for system
    type(TOrbitals), intent(in) :: orb

    !> Index for dense arrays
    integer, intent(in) :: iSquare(:)

    !> Neighbours of atoms
    integer, intent(in) :: iNeighbour(0:,:)

    !> Number of atoms neighbouring each site where the overlap is non-zero
    integer, intent(in) :: nNeighbourSK(:)

    !> Differentiation object
    class(TNonSccDiff), intent(in) :: derivator

    !> True, if system is periodic (i.e. Gamma-only)
    logical, intent(in) :: tPeriodic

    !> Energy gradients
    real(dp), intent(inout) :: gradients(:,:)

    !> List of neighbours for each atom (symmetric version)
    type(TSymNeighbourList), intent(in), optional :: symNeighbourList

    !> Nr. of neighbours for each atom
    integer, intent(in), optional :: nNeighbourCamSym(:)

    if (tPeriodic) then
      call this%tabulateCamdGammaEval0_gamma()
    else
      call this%tabulateCamdGammaEval0_cluster()
    end if

    select case(this%hybridXcAlg)
    case (hybridXcAlgo%thresholdBased, hybridXcAlgo%neighbourBased)
      if (tPeriodic) then
        @:ASSERT(present(symNeighbourList) .and. present(nNeighbourCamSym))
        call addCamGradientsNeighbour_gamma(this, deltaRhoSqr, skOverCont, symNeighbourList,&
            & nNeighbourCamSym, iSquare, orb, derivator, gradients)
      else
        call addCamGradientsNeighbour_cluster(this, derivator, deltaRhoSqr, skOverCont, orb,&
            & iSquare, SSqrReal, iNeighbour, nNeighbourSK, gradients)
      end if
    case (hybridXcAlgo%matrixBased)
      @:ASSERT(present(symNeighbourList) .and. present(nNeighbourCamSym))
      call addCamGradientsMatrix_real(this, deltaRhoSqr, SSqrReal, skOverCont,&
          & symNeighbourList, nNeighbourCamSym, iSquare, orb, derivator, gradients)
    end select

  end subroutine addCamGradients_real

#:endif


  !> Adds gradients due to CAM range-separated contributions, using the neighbour-list based
  !! formulation (non-periodic and Gamma-only version).
  !!
  !! PhD thesis of Vitalij Lutsker (2015)
  !! "Range-Separated Hybrid Functionals in the Density Functional-Based Tight-Binding Method"
  !! Appendix, Section E, Eq.(E.1-E.9)
  subroutine addCamGradientsNeighbour_cluster(this, derivator, deltaRhoSqr, skOverCont, orb,&
      & iSquare, SSqrReal, iNeighbour, nNeighbourSK, gradients)

    !> Class instance
    class(THybridXcFunc), intent(in) :: this

    !> Density matrix difference from reference q0
    real(dp), intent(in) :: deltaRhoSqr(:,:,:)

    !> Sparse overlap part
    type(TSlakoCont), intent(in) :: skOverCont

    !> Orbital information for system
    type(TOrbitals), intent(in) :: orb

    !> Index for dense arrays
    integer, intent(in) :: iSquare(:)

    !> Overlap matrix
    real(dp), intent(in) :: SSqrReal(:,:)

    !> Neighbours of atoms
    integer, intent(in) :: iNeighbour(0:,:)

    !> Number of atoms neighbouring each site where the overlap is non-zero
    integer, intent(in) :: nNeighbourSK(:)

    !> Differentiation object
    class(TNonSccDiff), intent(in) :: derivator

    !> Energy gradients
    real(dp), intent(inout) :: gradients(:,:)

    integer :: nAtom0, iAtK, iNeighK, iAtB, iNeighB, iAtC, iAtA, kpa
    real(dp) :: tmpgamma1, tmpgamma2
    real(dp) :: tmpforce(3), tmpforce_r(3), tmpforce2, tmpmultvar1
    integer :: nSpin, iSpin, mu, alpha, beta, ccc, kkk
    real(dp) :: sPrimeTmp(orb%mOrb, orb%mOrb, 3)
    real(dp) :: sPrimeTmp2(orb%mOrb, orb%mOrb, 3)
    real(dp), allocatable :: tmpOvr(:,:), tmpRho(:,:,:), tmpderiv(:,:)

    nSpin = size(deltaRhoSqr, dim=3)
    call allocateAndInit(tmpOvr, tmpRho, tmpderiv)
    nAtom0 = size(this%species0)
    tmpderiv(:,:) = 0.0_dp
    ! sum K
    loopK: do iAtK = 1, nAtom0
      ! C >= K
      loopC: do iNeighK = 0, nNeighbourSK(iAtK)
        iAtC = iNeighbour(iNeighK, iAtK)
        ! evaluate the ovr_prime
        sPrimeTmp2(:,:,:) = 0.0_dp
        sPrimeTmp(:,:,:) = 0.0_dp
        if (iAtK /= iAtC) then
          call derivator%getFirstDeriv(sPrimeTmp, skOverCont, this%rCoords, this%species0, iAtK,&
              & iAtC, orb)
          call derivator%getFirstDeriv(sPrimeTmp2, skOverCont, this%rCoords, this%species0, iAtC,&
              & iAtK, orb)
        end if
        loopB: do iAtB = 1, nAtom0
          ! A > B
          loopA: do iNeighB = 0, nNeighbourSK(iAtB)
            iAtA = iNeighbour(iNeighB, iAtB)
            tmpgamma1 = this%camGammaEval0(iAtK, iAtB) + this%camGammaEval0(iAtC, iAtB)
            tmpgamma2 = tmpgamma1 + this%camGammaEval0(iAtK, iAtA) + this%camGammaEval0(iAtC, iAtA)
            tmpforce(:) = 0.0_dp
            tmpforce_r(:) = 0.0_dp
            tmpforce2 = 0.0_dp
            ccc = 0
            do mu = iSquare(iAtC), iSquare(iAtC + 1) - 1
              ccc = ccc + 1
              kkk = 0
              do kpa = iSquare(iAtK), iSquare(iAtK + 1) - 1
                kkk = kkk + 1
                tmpmultvar1 = 0.0_dp
                do iSpin = 1, nSpin
                  do alpha = iSquare(iAtA), iSquare(iAtA + 1) - 1
                    do beta = iSquare(iAtB), iSquare(iAtB + 1) - 1
                      tmpmultvar1 = tmpmultvar1 + tmpOvr(beta, alpha)&
                          & * (tmpRho(beta, kpa, iSpin) * tmpRho(alpha, mu, iSpin)&
                          & + tmpRho(alpha, kpa, iSpin) * tmpRho(beta, mu, iSpin))
                    end do
                  end do
                end do
                tmpforce(:) = tmpforce + tmpmultvar1 * sPrimeTmp(ccc, kkk, :)
                tmpforce_r(:) = tmpforce_r + tmpmultvar1 * sPrimeTmp2(kkk, ccc, :)
                tmpforce2 = tmpforce2 + tmpmultvar1 * tmpOvr(kpa, mu)
              end do
            end do

            ! C /= K
            if (iAtK /= iAtC) then
              if (iAtB /= iAtA) then
                tmpforce(:) = tmpforce * tmpgamma2
                tmpforce_r(:) = tmpforce_r * tmpgamma2
                tmpforce(:) = tmpforce + tmpforce2 * (this%camdGammaEval0(:, iAtK, iAtA)&
                    & + this%camdGammaEval0(:, iAtK, iAtB))
                tmpforce_r(:) = tmpforce_r + tmpforce2 * (this%camdGammaEval0(:, iAtC, iAtA)&
                    & + this%camdGammaEval0(:, iAtC, iAtB))
              else
                tmpforce(:) = tmpforce * tmpgamma1
                tmpforce_r(:) = tmpforce_r * tmpgamma1
                tmpforce(:) = tmpforce + tmpforce2 * this%camdGammaEval0(:, iAtK, iAtA)
                tmpforce_r(:) = tmpforce_r + tmpforce2 * this%camdGammaEval0(:, iAtC, iAtA)
              end if
            else
              if (iAtB /= iAtA) then
                tmpforce(:) = tmpforce + tmpforce2 * (this%camdGammaEval0(:, iAtK, iAtA)&
                    & + this%camdGammaEval0(:, iAtK, iAtB))
              else
                tmpforce(:) = tmpforce + tmpforce2 * (this%camdGammaEval0(:, iAtK, iAtA))
              end if
            end if
            tmpderiv(:, iAtK) = tmpderiv(:, iAtK) + tmpforce
            tmpderiv(:, iAtC) = tmpderiv(:, iAtC) + tmpforce_r
          end do loopA
        end do loopB
      end do loopC
    end do loopK

    if (this%tREKS) then
      gradients(:,:) = gradients - 0.5_dp * tmpderiv
    else
      gradients(:,:) = gradients - 0.25_dp * nSpin * tmpderiv
    end if


  contains

    !> Initialise.
    subroutine allocateAndInit(tmpOvr, tmpRho, tmpderiv)

      !> Storage for the overlap
      real(dp), allocatable, intent(inout) :: tmpOvr(:,:)

      !> Storage for density matrix
      real(dp), allocatable, intent(inout) :: tmpRho(:,:,:)

      !> Workspace for the derivatives
      real(dp), allocatable, intent(inout) :: tmpderiv(:,:)

      !! Holds long-range gamma derivatives of a single interaction
      real(dp) :: tmp(3)

      !! Atom indices
      integer :: iAt1, iAt2

      !! Spin channel index
      integer :: iSpin

      !! Number of atoms
      integer :: nAtom0

      nAtom0 = size(this%species0)

      allocate(tmpOvr(size(SSqrReal, dim=1), size(SSqrReal, dim=1)))
      allocate(tmpRho(size(deltaRhoSqr, dim=1), size(deltaRhoSqr, dim=1), size(deltaRhoSqr, dim=3)))
      allocate(tmpderiv(3, size(gradients, dim=2)))
      tmpOvr(:,:) = SSqrReal
      tmpRho(:,:,:) = deltaRhoSqr

      call adjointLowerTriangle(tmpOvr)

      do iSpin = 1, size(deltaRhoSqr, dim=3)
        call adjointLowerTriangle(tmpRho(:,:, iSpin))
      end do

    end subroutine allocateAndInit

  end subroutine addCamGradientsNeighbour_cluster


  !> Adds CAM gradients due to CAM range-separated contributions, using the neighbour-list based
  !! formulation (Gamma-only version).
  subroutine addCamGradientsNeighbour_gamma(this, deltaRhoSqr, skOverCont, symNeighbourList,&
      & nNeighbourCamSym, iSquare, orb, derivator, gradients)

    !> Class instance
    class(THybridXcFunc), intent(inout), target :: this

    !> Square (unpacked) delta density matrix
    real(dp), intent(in) :: deltaRhoSqr(:,:,:)

    !> Sparse overlap container
    type(TSlakoCont), intent(in) :: skOverCont

    !> List of neighbours for each atom (symmetric version)
    type(TSymNeighbourList), intent(in) :: symNeighbourList

    !> Nr. of neighbours for each atom.
    integer, intent(in) :: nNeighbourCamSym(:)

    !> Position of each atom in the rows/columns of the square matrices. Shape: (nAtom)
    integer, intent(in) :: iSquare(:)

    !> Orbital information.
    type(TOrbitals), intent(in) :: orb

    !> Differentiation object
    class(TNonSccDiff), intent(in) :: derivator

    !> Energy gradients
    real(dp), intent(inout) :: gradients(:,:)

    !! Dense matrix descriptor indices
    integer, parameter :: descLen = 3, iStart = 1, iEnd = 2, iNOrb = 3

    !! Atom blocks from sparse, real-space overlap matrices
    real(dp), pointer :: pSam(:,:), pSbk(:,:)

    !! Stores start/end index and number of orbitals of square matrices
    integer :: descA(descLen), descB(descLen), descM(descLen)

    !! Stores start/end index and number of orbitals of square matrices
    integer :: descK(descLen)

    !! Temporary storages
    real(dp), allocatable :: tmpDeltaRhoSqr(:,:,:), tmpGradients(:,:)

    !! Number of atoms in central cell
    integer :: nAtom0

    !! Overlap matrix elements
    real(dp) :: Sam, Sbk

    !! Density matrix elements
    real(dp) :: dPab, dPmk

    !! Overlap derivatives
    real(dp), dimension(orb%mOrb, orb%mOrb, 3) :: SbnPrimeKequalsN

    !! \tilde{\gamma}_{\mu\kappa}, \tilde{\gamma}_{\alpha\kappa},
    !! \tilde{\gamma}_{\alpha\beta}
    real(dp) :: gammaMK, gammaAK, gammaAB, gammaMKMB, gammaTot

    !! 1st derivatives of
    !! \tilde{\gamma}_{\mu\nu}, \tilde{\gamma}_{\mu\beta},
    !! \tilde{\gamma}_{\alpha\nu}, \tilde{\gamma}_{\alpha\beta}
    real(dp), dimension(3) :: dGammaMK, dGammaAK

    !! Atom to calculate energy gradient components for
    integer :: iAtK, iNeighK

    !! Index of atom in central cell
    integer :: iAtM

    !! Neighbour indices (+corresponding atom indices)
    integer :: iNeighM, iAtA, iAtB

    !! Folded (to central cell) atom indices
    integer :: iAtAfold, iAtBfold

    !! Auxiliary variables for setting up 2D pointer to sparse overlap
    integer :: ind, nOrbAt, nOrbNeigh

    !! Spin index and total number of spin channels
    integer :: iSpin, nSpin

    !! Orbital indices
    integer :: mu, kk, alpha, beta

    !! Product dPmk * Sam
    real(dp) :: dPmkSam

    !! Product dPmk * gammaTot
    real(dp) :: dPmkGammaTot

    !! Product dPmk * gammaTot * Sam
    real(dp) :: dPmkGammaTotSam

    !! Composite index iAtM/iAtN
    integer, allocatable :: iAtMN(:,:)

    nAtom0 = size(this%species0)
    nSpin = size(deltaRhoSqr, dim=3)

    call getTwoLoopCompositeIndex(nAtom0, nAtom0, iAtMN)

    ! allocate gradient contribution
    allocate(tmpGradients(3, size(gradients, dim=2)))
    tmpGradients(:,:) = 0.0_dp

    tmpDeltaRhoSqr = deltaRhoSqr
    do iSpin = 1, nSpin
      call adjointLowerTriangle(tmpDeltaRhoSqr(:,:, iSpin))
    end do

    loopK1: do iAtK = 1, nAtom0
      descK = getDescriptor(iAtK, iSquare)
      loopM1: do iAtM = 1, nAtom0
        descM = getDescriptor(iAtM, iSquare)
        gammaMK = this%camGammaEval0(iAtM, iAtK)
        loopB1: do iNeighK = 0, nNeighbourCamSym(iAtK)
          iAtB = symNeighbourList%neighbourList%iNeighbour(iNeighK, iAtK)
          iAtBfold = symNeighbourList%img2CentCell(iAtB)
          if (iAtK == iAtBfold) cycle
          descB = getDescriptor(iAtBfold, iSquare)
          ! \tilde{\gamma}_{\mu\beta}
          gammaMKMB = gammaMK + this%camGammaEval0(iAtM, iAtBfold)
          SbnPrimeKequalsN(:,:,:) = 0.0_dp
          call derivator%getFirstDeriv(SbnPrimeKequalsN, skOverCont, this%rCoords,&
              & symNeighbourList%species, iAtK, iAtB, orb)
          loopA1: do iNeighM = 0, nNeighbourCamSym(iAtM)
            iAtA = symNeighbourList%neighbourList%iNeighbour(iNeighM, iAtM)
            iAtAfold = symNeighbourList%img2CentCell(iAtA)
            descA = getDescriptor(iAtAfold, iSquare)
            ! \tilde{\gamma}_{\alpha\nu}
            gammaAK = this%camGammaEval0(iAtAfold, iAtK)
            ! \tilde{\gamma}_{\alpha\beta}
            gammaAB = this%camGammaEval0(iAtAfold, iAtBfold)
            ! get 2D pointer to Sam overlap block
            ind = symNeighbourList%iPair(iNeighM, iAtM) + 1
            nOrbAt = descM(iNOrb)
            nOrbNeigh = descA(iNOrb)
            pSam(1:nOrbNeigh, 1:nOrbAt) => this%overSym(ind:ind + nOrbNeigh * nOrbAt - 1)

            gammaTot = gammaMKMB + gammaAK + gammaAB

            do mu = 1, descM(iNOrb)
              do kk = 1, descK(iNOrb)
                do iSpin = 1, nSpin
                  dPmk = tmpDeltaRhoSqr(descM(iStart) + mu - 1, descK(iStart) + kk - 1, iSpin)
                  dPmkGammaTot = dPmk * gammaTot
                  do alpha = 1, descA(iNOrb)
                    Sam = pSam(alpha, mu)
                    dPmkGammaTotSam = dPmkGammaTot * Sam
                    do beta = 1, descB(iNOrb)
                      dPab = tmpDeltaRhoSqr(descA(iStart) + alpha - 1,&
                          & descB(iStart) + beta - 1, iSpin)

                      tmpGradients(:, iAtK) = tmpGradients(:, iAtK)&
                          & + dPmkGammaTotSam * dPab * SbnPrimeKequalsN(beta, kk, :)
                    end do
                  end do
                end do
              end do
            end do

          end do loopA1
        end do loopB1
      end do loopM1
    end do loopK1

    loopK2: do iAtK = 1, nAtom0
      descK = getDescriptor(iAtK, iSquare)
      loopM2: do iAtM = 1, nAtom0
        if (iAtK == iAtM) cycle
        descM = getDescriptor(iAtM, iSquare)
        dGammaMK(:) = -this%camdGammaEval0(:, iAtM, iAtK)
        loopB2: do iNeighK = 0, nNeighbourCamSym(iAtK)
          iAtB = symNeighbourList%neighbourList%iNeighbour(iNeighK, iAtK)
          iAtBfold = symNeighbourList%img2CentCell(iAtB)
          descB = getDescriptor(iAtBfold, iSquare)
          ! get 2D pointer to Sbk overlap block
          ind = symNeighbourList%iPair(iNeighK, iAtK) + 1
          nOrbAt = descK(iNOrb)
          nOrbNeigh = descB(iNOrb)
          pSbk(1:nOrbNeigh, 1:nOrbAt) => this%overSym(ind:ind + nOrbNeigh * nOrbAt - 1)
          loopA2: do iNeighM = 0, nNeighbourCamSym(iAtM)
            iAtA = symNeighbourList%neighbourList%iNeighbour(iNeighM, iAtM)
            iAtAfold = symNeighbourList%img2CentCell(iAtA)
            descA = getDescriptor(iAtAfold, iSquare)
            ! get 2D pointer to Sam overlap block
            ind = symNeighbourList%iPair(iNeighM, iAtM) + 1
            nOrbAt = descM(iNOrb)
            nOrbNeigh = descA(iNOrb)
            pSam(1:nOrbNeigh, 1:nOrbAt) => this%overSym(ind:ind + nOrbNeigh * nOrbAt - 1)

            do mu = 1, descM(iNOrb)
              do kk = 1, descK(iNOrb)
                do iSpin = 1, nSpin
                  dPmk = tmpDeltaRhoSqr(descM(iStart) + mu - 1, descK(iStart) + kk - 1, iSpin)
                  do alpha = 1, descA(iNOrb)
                    Sam = pSam(alpha, mu)
                    dPmkSam = dPmk * Sam
                    do beta = 1, descB(iNOrb)
                      Sbk = pSbk(beta, kk)
                      dPab = tmpDeltaRhoSqr(descA(iStart) + alpha - 1,&
                          & descB(iStart) + beta - 1, iSpin)

                      tmpGradients(:, iAtK) = tmpGradients(:, iAtK)&
                          & + dPmkSam * dPab * Sbk * dGammaMK
                    end do
                  end do
                end do
              end do
            end do

          end do loopA2
        end do loopB2
      end do loopM2
    end do loopK2

    loopK3: do iAtK = 1, nAtom0
      descK = getDescriptor(iAtK, iSquare)
      loopM3: do iAtM = 1, nAtom0
        descM = getDescriptor(iAtM, iSquare)
        loopA3: do iNeighM = 0, nNeighbourCamSym(iAtM)
          iAtA = symNeighbourList%neighbourList%iNeighbour(iNeighM, iAtM)
          iAtAfold = symNeighbourList%img2CentCell(iAtA)
          if (iAtK == iAtAfold) cycle
          dGammaAK(:) = -this%camdGammaEval0(:, iAtAfold, iAtK)
          descA = getDescriptor(iAtAfold, iSquare)
          ! get 2D pointer to Sam overlap block
          ind = symNeighbourList%iPair(iNeighM, iAtM) + 1
          nOrbAt = descM(iNOrb)
          nOrbNeigh = descA(iNOrb)
          pSam(1:nOrbNeigh, 1:nOrbAt) => this%overSym(ind:ind + nOrbNeigh * nOrbAt - 1)
          loopB3: do iNeighK = 0, nNeighbourCamSym(iAtK)
            iAtB = symNeighbourList%neighbourList%iNeighbour(iNeighK, iAtK)
            iAtBfold = symNeighbourList%img2CentCell(iAtB)
            descB = getDescriptor(iAtBfold, iSquare)
            ! get 2D pointer to Sbk overlap block
            ind = symNeighbourList%iPair(iNeighK, iAtK) + 1
            nOrbAt = descK(iNOrb)
            nOrbNeigh = descB(iNOrb)
            pSbk(1:nOrbNeigh, 1:nOrbAt) => this%overSym(ind:ind + nOrbNeigh * nOrbAt - 1)

            do mu = 1, descM(iNOrb)
              do kk = 1, descK(iNOrb)
                do iSpin = 1, nSpin
                  dPmk = tmpDeltaRhoSqr(descM(iStart) + mu - 1, descK(iStart) + kk - 1, iSpin)
                  do alpha = 1, descA(iNOrb)
                    Sam = pSam(alpha, mu)
                    dPmkSam = dPmk * Sam
                    do beta = 1, descB(iNOrb)
                      Sbk = pSbk(beta, kk)
                      dPab = tmpDeltaRhoSqr(descA(iStart) + alpha - 1,&
                          & descB(iStart) + beta - 1, iSpin)

                      tmpGradients(:, iAtK) = tmpGradients(:, iAtK)&
                          & + dPmkSam * dPab * Sbk * dGammaAK
                    end do
                  end do
                end do
              end do
            end do

          end do loopB3
        end do loopA3
      end do loopM3
    end do loopK3

    if (this%tREKS) then
      gradients(:,:) = gradients - 0.5_dp * tmpGradients
    else
      gradients(:,:) = gradients - 0.25_dp * nSpin * tmpGradients
    end if

  end subroutine addCamGradientsNeighbour_gamma


  !> Pre-tabulates summed gamma function derivatives (non-periodic version).
  subroutine tabulateCamdGammaEval0_cluster(this)

    !> Class instance
    class(THybridXcFunc), intent(inout) :: this

    !! Holds long-range gamma derivatives of a single interaction
    real(dp) :: tmp(3)

    !! Indices of interacting atoms in central cell
    integer :: iAt1, iAt2

    !! Number of atoms in central cell
    integer :: nAtom0

    nAtom0 = size(this%species0)

    if (allocated(this%camdGammaEval0)) deallocate(this%camdGammaEval0)
    allocate(this%camdGammaEval0(3, nAtom0, nAtom0))
    this%camdGammaEval0(:,:,:) = 0.0_dp

    do iAt1 = 1, nAtom0
      do iAt2 = 1, nAtom0
        if (iAt1 /= iAt2) then
          call getDirectionalCamGammaPrimeValue_cluster(this, tmp, iAt1, iAt2)
          this%camdGammaEval0(:, iAt1, iAt2) = tmp
        end if
      end do
    end do

  end subroutine tabulateCamdGammaEval0_cluster


  !> Pre-tabulates summed gamma function derivatives (Gamma-only version).
  subroutine tabulateCamdGammaEval0_gamma(this)

    !> Class instance
    class(THybridXcFunc), intent(inout) :: this

    !! Indices of interacting atoms in central cell, as well as their global species index
    integer :: iAtM, iSpM, iAtN, iSpN

    !! Number of atoms in central cell
    integer :: nAtom0

    nAtom0 = size(this%species0)

    if (allocated(this%camdGammaEval0)) deallocate(this%camdGammaEval0)
    allocate(this%camdGammaEval0(3, nAtom0, nAtom0))
    this%camdGammaEval0(:,:,:) = 0.0_dp

    do iAtM = 1, nAtom0
      iSpM = this%species0(iAtM)
      do iAtN = 1, nAtom0
        if (iAtM == iAtN) cycle
        iSpN = this%species0(iAtN)
        this%camdGammaEval0(:, iAtM, iAtN) = getCamGammaPrimeGSum(this, iAtM, iAtN, iSpM, iSpN,&
            & this%rCellVecsG)
        this%camdGammaEval0(:, iAtN, iAtM) = getCamGammaPrimeGSum(this, iAtN, iAtM, iSpN, iSpM,&
            & this%rCellVecsG)
      end do
    end do

  end subroutine tabulateCamdGammaEval0_gamma

#:if WITH_SCALAPACK

  !> Collects full, collected square matrix from distributed contributions.
  subroutine getFullFromDistributed(env, denseDesc, parallelKS, distrib, collected)

    !> Environment settings
    type(TEnvironment), intent(in) :: env

    !> Dense matrix descriptor
    type(TDenseDescr), intent(in) :: denseDesc

    !> The k-points and spins to process
    type(TParallelKS), intent(in) :: parallelKS

    !> Distributed matrix (source)
    real(dp), intent(in) :: distrib(:,:,:)

    !> Collected, full square matrix (target)
    real(dp), intent(out), allocatable :: collected(:,:,:)

    !! Spin and composite index
    integer :: nSpin, iSpin, iKS

    !! Auxiliary variables
    integer :: iGroup, iEig, nOrb

    !! Type for communicating a row or a column of a distributed matrix
    type(linecomm) :: collector

    !! Temporary storage for a single line of the collected, dense, square matrix
    real(dp), allocatable :: localLine(:)

    nOrb = denseDesc%fullSize
    nSpin = maxval(parallelKS%groupKS(2, :, 0:))
    allocate(localLine(nOrb))
    if (env%mpi%tGlobalLead) allocate(collected(nOrb, nOrb, nSpin), source=0.0_dp)

    call collector%init(env%blacs%orbitalGrid, denseDesc%blacsOrbSqr, "c")

    leadOrFollow: if (env%mpi%tGlobalLead) then
      do iKS = 1, parallelKS%maxGroupKS
        group: do iGroup = 0, env%mpi%nGroup - 1
          if (iKS > parallelKS%nGroupKS(iGroup)) then
            cycle group
          end if
          iSpin = parallelKS%groupKS(2, iKS, iGroup)
          do iEig = 1, nOrb
            if (iGroup == 0) then
              call collector%getline_lead(env%blacs%orbitalGrid, iEig, distrib(:,:,iKS),&
                  & localLine)
            else
              call mpifx_recv(env%mpi%interGroupComm, localLine, iGroup)
            end if
            collected(:, iEig, iSpin) = localLine
          end do
        end do group
      end do
    else
      do iKS = 1, parallelKS%nLocalKS
        do iEig = 1, nOrb
          if (env%mpi%tGroupLead) then
            call collector%getline_lead(env%blacs%orbitalGrid, iEig, distrib(:,:,iKS),&
                & localLine)
            call mpifx_send(env%mpi%interGroupComm, localLine, env%mpi%interGroupComm%leadrank)
          else
            call collector%getline_follow(env%blacs%orbitalGrid, iEig, distrib(:,:,iKS))
          end if
        end do
      end do
    end if leadOrFollow

  end subroutine getFullFromDistributed


  !> Scatters full, collected square matrix to distributed matrices.
  subroutine scatterFullToDistributed(env, denseDesc, parallelKS, collected, distrib)

    !> Environment settings
    type(TEnvironment), intent(in) :: env

    !> Dense matrix descriptor
    type(TDenseDescr), intent(in) :: denseDesc

    !> The k-points and spins to process
    type(TParallelKS), intent(in) :: parallelKS

    !> Collected, full square matrix (source)
    real(dp), intent(in), allocatable :: collected(:,:,:)

    !> Distributed matrix (target)
    real(dp), intent(out) :: distrib(:,:,:)

    !! Spin and composite index
    integer :: iSpin, iKS

    !! Auxiliary variables
    integer :: iGroup, iEig, nOrb

    !! Type for communicating a row or a column of a distributed matrix
    type(linecomm) :: collector

    !! Temporary storage for a single line of the collected, dense, square matrix
    real(dp), allocatable :: localLine(:)

    nOrb = denseDesc%fullSize
    allocate(localLine(nOrb))
    distrib(:,:,:) = 0.0_dp

    call collector%init(env%blacs%orbitalGrid, denseDesc%blacsOrbSqr, "c")

    leadOrFollow: if (env%mpi%tGlobalLead) then
      do iKS = 1, parallelKS%maxGroupKS
        group: do iGroup = 0, env%mpi%nGroup - 1
          if (iKS > parallelKS%nGroupKS(iGroup)) then
            cycle group
          end if
          iSpin = parallelKS%groupKS(2, iKS, iGroup)
          do iEig = 1, nOrb
            if (iGroup == 0) then
              call collector%setline_lead(env%blacs%orbitalGrid, iEig, collected(:, iEig, iSpin),&
                  & distrib(:,:,iKS))
            else
              call mpifx_send(env%mpi%interGroupComm, collected(:, iEig, iSpin), iGroup)
            end if
          end do
        end do group
      end do
    else
      do iKS = 1, parallelKS%nLocalKS
        do iEig = 1, nOrb
          if (env%mpi%tGroupLead) then
            call mpifx_recv(env%mpi%interGroupComm, localLine, env%mpi%interGroupComm%leadrank)
            call collector%setline_lead(env%blacs%orbitalGrid, iEig, localLine, distrib(:,:,iKS))
          else
            call collector%setline_follow(env%blacs%orbitalGrid, iEig, distrib(:,:,iKS))
          end if
        end do
      end do
    end if leadOrFollow

  end subroutine scatterFullToDistributed


  !> Adds CAM gradients due to CAM range-separated contributions, using a matrix-based formulation.
  !! (non-periodic and Gamma-only version)
  !!
  !! Eq.(B5) of Phys. Rev. Materials 7, 063802 (DOI: 10.1103/PhysRevMaterials.7.063802)
  subroutine addCamGradientsMatrix_real_blacs(this, env, parallelKS, deltaRhoSqr,&
      & overlap, skOverCont, symNeighbourList, nNeighbourCamSym, orb, derivator, denseDesc,&
      & nSpin, gradients)

    !> Class instance
    class(THybridXcFunc), intent(inout), target :: this

    !> Environment settings
    type(TEnvironment), intent(in) :: env

    !> The k-points and spins to process
    type(TParallelKS), intent(in) :: parallelKS

    !> Square (unpacked) delta density matrix
    real(dp), intent(in) :: deltaRhoSqr(:,:,:)

    !> Square (unpacked) overlap matrix
    real(dp), intent(in) :: overlap(:,:)

    !> Sparse overlap container
    type(TSlakoCont), intent(in) :: skOverCont

    !> List of neighbours for each atom (symmetric version)
    type(TSymNeighbourList), intent(in) :: symNeighbourList

    !> Nr. of neighbours for each atom.
    integer, intent(in) :: nNeighbourCamSym(:)

    !> Orbital information
    type(TOrbitals), intent(in) :: orb

    !> Differentiation object
    class(TNonSccDiff), intent(in) :: derivator

    !> Dense matrix descriptor
    type(TDenseDescr), intent(in) :: denseDesc

    !> Number of spin channels in calculation
    integer, intent(in) :: nSpin

    !> Energy gradients
    real(dp), intent(inout) :: gradients(:,:)

    !! Temporary energy gradients
    real(dp), allocatable :: tmpGradients(:,:)

    !! Number of atoms in central cell
    integer :: nAtom0

    !! Overlap derivative
    real(dp), allocatable :: overSqrPrime(:,:,:)

    !! Atom indices (central cell)
    integer :: iAtForce, iAt1, iAt2

    !! CAM gamma matrix, including periodic images
    real(dp), allocatable :: camGammaAO(:,:), camdGammaAO(:,:,:)

    !! Temporary storages
    real(dp), allocatable :: symSqrMat1(:,:,:), symSqrMat2(:,:,:), symSqrMatTmp(:,:)
    real(dp), allocatable :: deltaRhoOverlap(:,:,:)

    !! Composite index (iK, iS)
    integer :: iKS

    !! Iterates over coordinates
    integer :: iCoord

    !! Number of local rows and columns
    integer :: nLocRow, nLocCol

    !! Auxiliary variables
    integer :: ii, jj, iOrb1, iOrb2

    nLocCol = size(overlap, dim=1)
    nLocRow = size(overlap, dim=2)

    nAtom0 = size(this%species0)

    ! allocate CAM \tilde{gamma}
    allocate(camGammaAO(nLocCol, nLocRow))
    allocate(camdGammaAO(nLocCol, nLocRow, 3))

    ! get CAM \tilde{gamma} super-matrix
    do jj = 1, size(camGammaAO, dim=2)
      iOrb1 = scalafx_indxl2g(jj, denseDesc%blacsOrbSqr(NB_), env%blacs%orbitalGrid%mycol,&
          & denseDesc%blacsOrbSqr(CSRC_), env%blacs%orbitalGrid%ncol)
      call bisection(iAt1, denseDesc%iAtomStart, iOrb1)
      do ii = 1, size(camGammaAO, dim=1)
        iOrb2 = scalafx_indxl2g(ii, denseDesc%blacsOrbSqr(MB_), env%blacs%orbitalGrid%myrow,&
            & denseDesc%blacsOrbSqr(RSRC_), env%blacs%orbitalGrid%nrow)
        call bisection(iAt2, denseDesc%iAtomStart, iOrb2)
        camGammaAO(ii, jj) = this%camGammaEval0(iAt2, iAt1)
      end do
    end do

    ! pre-calculate deltaRho * overlap
    allocate(deltaRhoOverlap, mold=deltaRhoSqr)
    do iKS = 1, parallelKS%nLocalKS
      call pblasfx_psymm(deltaRhoSqr(:,:, iKS), denseDesc%blacsOrbSqr, overlap,&
          & denseDesc%blacsOrbSqr, deltaRhoOverlap(:,:, iKS), denseDesc%blacsOrbSqr, side="L")
    end do

    ! calculate first symmetrized square matrix of Eq.(B5)
    allocate(symSqrMat1, mold=deltaRhoSqr)
    allocate(symSqrMatTmp(size(symSqrMat1, dim=1), size(symSqrMat1, dim=2)))
    do iKS = 1, parallelKS%nLocalKS
      call pblasfx_pgemm(deltaRhoOverlap(:,:, iKS), denseDesc%blacsOrbSqr,&
          & deltaRhoSqr(:,:, iKS) * camGammaAO, denseDesc%blacsOrbSqr, symSqrMat1(:,:, iKS),&
          & denseDesc%blacsOrbSqr)
      call pblasfx_pgemm(deltaRhoOverlap(:,:, iKS) * camGammaAO, denseDesc%blacsOrbSqr,&
          & deltaRhoSqr(:,:, iKS), denseDesc%blacsOrbSqr, symSqrMat1(:,:, iKS),&
          & denseDesc%blacsOrbSqr, beta=1.0_dp)
      symSqrMatTmp(:,:) = symSqrMat1(:,:, iKS)
      ! symmetrize temporary storage
      call pblasfx_ptran(symSqrMatTmp, denseDesc%blacsOrbSqr, symSqrMat1(:,:, iKS),&
          & denseDesc%blacsOrbSqr, alpha=0.5_dp, beta=0.5_dp)
    end do

    ! free some memory
    deallocate(camGammaAO)

    ! calculate second symmetrized square matrix of Eq.(B5)
    allocate(symSqrMat2, mold=deltaRhoSqr)
        do iKS = 1, parallelKS%nLocalKS
      call pblasfx_ptran(deltaRhoOverlap(:,:, iKS), denseDesc%blacsOrbSqr,&
          & symSqrMat2(:,:, iKS), denseDesc%blacsOrbSqr, alpha=1.0_dp, beta=0.0_dp)
      symSqrMat2(:,:, iKS) = symSqrMat2(:,:, iKS) * deltaRhoOverlap(:,:, iKS)
      call pblasfx_pgemm(overlap, denseDesc%blacsOrbSqr, deltaRhoOverlap(:,:, iKS),&
          & denseDesc%blacsOrbSqr, symSqrMatTmp, denseDesc%blacsOrbSqr)
      symSqrMat2(:,:, iKS) = symSqrMat2(:,:, iKS) + symSqrMatTmp * deltaRhoSqr(:,:, iKS)
      ! symmetrize temporary storage
      symSqrMatTmp(:,:) = symSqrMat2(:,:, iKS)
      call pblasfx_ptran(symSqrMatTmp, denseDesc%blacsOrbSqr, symSqrMat2(:,:, iKS),&
          & denseDesc%blacsOrbSqr, alpha=0.5_dp, beta=0.5_dp)
    end do

    ! free some memory
    deallocate(symSqrMatTmp)
    deallocate(deltaRhoOverlap)

    allocate(overSqrPrime(3, nLocCol, nLocRow))
    allocate(tmpGradients, mold=gradients)
    tmpGradients(:,:) = 0.0_dp

    loopForceAtom: do iAtForce = 1, nAtom0
      call getUnpackedOverlapPrime_real(env, denseDesc, iAtForce, skOverCont, orb, derivator,&
          & symNeighbourList, nNeighbourCamSym, denseDesc%iAtomStart, this%rCoords, overSqrPrime)
      call getUnpackedCamGammaAOPrime(env, denseDesc, orb, iAtForce, this%camdGammaEval0,&
          & camdGammaAO)

      do iKS = 1, parallelKS%nLocalKS
        do iCoord = 1, 3
          ! first term of Eq.(B5)
          tmpGradients(iCoord, iAtForce) = tmpGradients(iCoord, iAtForce)&
              & - sum(overSqrPrime(iCoord, :,:) * symSqrMat1(:,:, iKS))
          ! second term of Eq.(B5)
          tmpGradients(iCoord, iAtForce) = tmpGradients(iCoord, iAtForce)&
              & - 0.5_dp * sum(camdGammaAO(:,:, iCoord) * symSqrMat2(:,:, iKS))
        end do
      end do
    end do loopForceAtom

    call mpifx_allreduceip(env%mpi%globalComm, tmpGradients, MPI_SUM)

    if (this%tREKS) then
      gradients(:,:) = gradients + tmpGradients
    else
      gradients(:,:) = gradients + 0.5_dp * nSpin * tmpGradients
    end if

  end subroutine addCamGradientsMatrix_real_blacs


  !> Calculates derivative of the (long-range + HF full-range) gamma interaction w.r.t. given atom.
  !!
  !! 2nd term of Eq.(B5) of Phys. Rev. Materials 7, 063802 (DOI: 10.1103/PhysRevMaterials.7.063802)
  subroutine getUnpackedCamGammaAOPrime(env, denseDesc, orb, iAtomPrime, camdGammaEval0,&
      & camdGammaAO)

    !> Environment settings
    type(TEnvironment), intent(in) :: env

    !> Dense matrix descriptor
    type(TDenseDescr), intent(in) :: denseDesc

    !> Orbital information
    type(TOrbitals), intent(in) :: orb

    !> Overlap derivative calculated w.r.t. this atom in the central cell
    integer, intent(in) :: iAtomPrime

    !> Pre-tabulated (long-range + HF full-range) gamma derivative of At1 and At2 (in central cell)
    real(dp), intent(in) :: camdGammaEval0(:,:,:)

    !> CAM gamma matrix, including periodic images
    real(dp), intent(out) :: camdGammaAO(:,:,:)

    !! Atom interacting with iAtomPrime
    integer :: iAt2

    !! Gamma derivatives of current diatomic block
    real(dp) :: tmp(orb%nOrb, orb%nOrb)

    !! Auxiliary variables
    integer :: iOrbStart1, iOrbEnd1, nOrbAtomPrime, iOrbStart2, iOrbEnd2, nOrb2, iCoord

    camdGammaAO(:,:,:) = 0.0_dp

    iOrbStart1 = denseDesc%iAtomStart(iAtomPrime)
    iOrbEnd1 = denseDesc%iAtomStart(iAtomPrime + 1) - 1
    nOrbAtomPrime = iOrbEnd1 - iOrbStart1 + 1

    do iAt2 = 1, size(camdGammaEval0, dim=3)
      iOrbStart2 = denseDesc%iAtomStart(iAt2)
      iOrbEnd2 = denseDesc%iAtomStart(iAt2 + 1) - 1
      nOrb2 = iOrbEnd2 - iOrbStart2 + 1
      do iCoord = 1, 3
        tmp(1:nOrbAtomPrime, 1:nOrb2) = camdGammaEval0(iCoord, iAtomPrime, iAt2)
        call scalafx_addl2g(env%blacs%orbitalGrid, tmp(1:nOrbAtomPrime, 1:nOrb2),&
            & denseDesc%blacsOrbSqr, iOrbStart1, iOrbStart2, camdGammaAO(:,:, iCoord))
      end do
    end do

  end subroutine getUnpackedCamGammaAOPrime


  !> Calculates the derivative of the square, dense, unpacked, Gamma-point overlap matrix w.r.t.
  !! given atom.
  !!
  !! 1st term of Eq.(B5) of Phys. Rev. Materials 7, 063802 (DOI: 10.1103/PhysRevMaterials.7.063802)
  subroutine getUnpackedOverlapPrime_real(env, denseDesc, iAtomPrime, skOverCont, orb, derivator,&
      & symNeighbourList, nNeighbourCamSym, iSquare, rCoords, overSqrPrime)

    !> Environment settings
    type(TEnvironment), intent(in) :: env

    !> Dense matrix descriptor
    type(TDenseDescr), intent(in) :: denseDesc

    !> Overlap derivative calculated w.r.t. this atom in the central cell
    integer, intent(in) :: iAtomPrime

    !> Sparse overlap container
    type(TSlakoCont), intent(in) :: skOverCont

    !> Orbital information.
    type(TOrbitals), intent(in) :: orb

    !> Differentiation object
    class(TNonSccDiff), intent(in) :: derivator

    !> List of neighbours for each atom (symmetric version)
    type(TSymNeighbourList), intent(in) :: symNeighbourList

    !> Nr. of neighbours for each atom.
    integer, intent(in) :: nNeighbourCamSym(:)

    !> Position of each atom in the rows/columns of the square matrices. Shape: (nAtom)
    integer, intent(in) :: iSquare(:)

    !> Atomic coordinates in absolute units, potentially including periodic images
    real(dp), intent(in) :: rCoords(:,:)

    !> Overlap derivative
    real(dp), intent(out) :: overSqrPrime(:,:,:)

    !! Temporary overlap derivative storage
    real(dp) :: overPrime(orb%mOrb, orb%mOrb, 3)

    !! Dense matrix descriptor indices
    integer, parameter :: descLen = 3, iStart = 1, iEnd = 2, iNOrb = 3

    !! Stores start/end index and number of orbitals of square matrices
    integer :: descAt1(descLen), descAt2(descLen)

    !! Atom to calculate energy gradient components for
    integer :: iAt2, iAt2fold, iNeigh

    !! Iterates over coordinates
    integer :: iCoord

    overSqrPrime(:,:,:) = 0.0_dp

    descAt1 = getDescriptor(iAtomPrime, iSquare)
    do iNeigh = 0, nNeighbourCamSym(iAtomPrime)
      iAt2 = symNeighbourList%neighbourList%iNeighbour(iNeigh, iAtomPrime)
      iAt2fold = symNeighbourList%img2CentCell(iAt2)
      if (iAtomPrime == iAt2fold) cycle
      descAt2 = getDescriptor(iAt2fold, iSquare)
      overPrime(:,:,:) = 0.0_dp
      call derivator%getFirstDeriv(overPrime, skOverCont, rCoords, symNeighbourList%species,&
          & iAtomPrime, iAt2, orb)
      do iCoord = 1, 3
        call scalafx_addl2g(env%blacs%orbitalGrid, overPrime(1:descAt2(iNOrb), 1:descAt1(iNOrb),&
            & iCoord), denseDesc%blacsOrbSqr, descAt2(iStart), descAt1(iStart),&
            & overSqrPrime(iCoord, :,:))
      end do
    end do

  end subroutine getUnpackedOverlapPrime_real

#:else

  !> Adds CAM gradients due to CAM range-separated contributions, using a matrix-based formulation.
  !! (non-periodic and Gamma-only version)
  !!
  !! Eq.(B5) of Phys. Rev. Materials 7, 063802 (DOI: 10.1103/PhysRevMaterials.7.063802)
  subroutine addCamGradientsMatrix_real(this, deltaRhoSqr, overlap, skOverCont,&
      & symNeighbourList, nNeighbourCamSym, iSquare, orb, derivator, gradients)

    !> Class instance
    class(THybridXcFunc), intent(inout), target :: this

    !> Square (unpacked) delta density matrix
    real(dp), intent(in) :: deltaRhoSqr(:,:,:)

    !> Square (unpacked) overlap matrix
    real(dp), intent(in) :: overlap(:,:)

    !> Sparse overlap container
    type(TSlakoCont), intent(in) :: skOverCont

    !> List of neighbours for each atom (symmetric version)
    type(TSymNeighbourList), intent(in) :: symNeighbourList

    !> Nr. of neighbours for each atom.
    integer, intent(in) :: nNeighbourCamSym(:)

    !> Position of each atom in the rows/columns of the square matrices. Shape: (nAtom)
    integer, intent(in) :: iSquare(:)

    !> Orbital information.
    type(TOrbitals), intent(in) :: orb

    !> Differentiation object
    class(TNonSccDiff), intent(in) :: derivator

    !> Energy gradients
    real(dp), intent(inout) :: gradients(:,:)

    !! Temporary energy gradients
    real(dp), allocatable :: tmpGradients(:,:)

    !! Dense matrix descriptor indices
    integer, parameter :: descLen = 3, iStart = 1, iEnd = 2, iNOrb = 3

    !! Number of atoms in central cell
    integer :: nAtom0

    !! Overlap derivative
    real(dp), allocatable :: overSqrPrime(:,:,:)

    !! Atom indices (central cell)
    integer :: iAtForce, iAt1, iAt2

    !! CAM gamma matrix, including periodic images
    real(dp), allocatable :: camGammaAO(:,:), camdGammaAO(:,:,:)

    !! Symmetrized, square (unpacked) overlap matrix
    real(dp), allocatable :: overlapSym(:,:)

    !! Symmetrized, square (unpacked) delta density matrix
    real(dp), allocatable :: deltaRhoSqrSym(:,:,:)

    !! Temporary storages
    real(dp), allocatable :: symSqrMat1(:,:,:), symSqrMat2(:,:,:), symSqrMat2Tmp(:,:)
    real(dp), allocatable :: deltaRhoOverlap(:,:,:)

    !! Spin index and total number of spin channels
    integer :: iSpin, nSpin

    !! Number of orbitals in square matrices
    integer :: nOrb

    !! Iterates over coordinates
    integer :: iCoord

    nAtom0 = size(this%species0)
    nSpin = size(deltaRhoSqr, dim=3)
    nOrb = size(overlap, dim=1)

    ! allocate CAM \tilde{gamma}
    allocate(camGammaAO(nOrb, nOrb))
    allocate(camdGammaAO(nOrb, nOrb, 3))

    ! symmetrize square overlap and density matrix
    overlapSym = overlap
    call adjointLowerTriangle(overlapSym)
    deltaRhoSqrSym = deltaRhoSqr
    do iSpin = 1, nSpin
      call adjointLowerTriangle(deltaRhoSqrSym(:,:, iSpin))
    end do

    ! get CAM \tilde{gamma} super-matrix
    do iAt2 = 1, nAtom0
      do iAt1 = 1, nAtom0
        camGammaAO(iSquare(iAt1):iSquare(iAt1 + 1) - 1, iSquare(iAt2):iSquare(iAt2 + 1) - 1)&
            & = this%camGammaEval0(iAt1, iAt2)
      end do
    end do

    ! pre-calculate deltaRho * overlap
    allocate(deltaRhoOverlap, mold=deltaRhoSqrSym)
    do iSpin = 1, nSpin
      call symm(deltaRhoOverlap(:,:, iSpin), "l", deltaRhoSqrSym(:,:, iSpin), overlapSym)
    end do

    ! calculate first symmetrized square matrix of Eq.(B5)
    allocate(symSqrMat1, mold=deltaRhoSqrSym)
    do iSpin = 1, nSpin
      call gemm(symSqrMat1(:,:, iSpin),&
          & deltaRhoOverlap(:,:, iSpin),&
          & deltaRhoSqrSym(:,:, iSpin) * camGammaAO)
      call gemm(symSqrMat1(:,:, iSpin),&
          & deltaRhoOverlap(:,:, iSpin) * camGammaAO,&
          & deltaRhoSqrSym(:,:, iSpin),&
          & beta=1.0_dp)
      ! symmetrize temporary storage
      symSqrMat1(:,:, iSpin) = 0.5_dp * (symSqrMat1(:,:, iSpin) + transpose(symSqrMat1(:,:, iSpin)))
    end do

    ! free some memory
    deallocate(camGammaAO)

    ! calculate second symmetrized square matrix of Eq.(B5)
    allocate(symSqrMat2, mold=deltaRhoSqrSym)
    allocate(symSqrMat2Tmp(size(symSqrMat2, dim=1), size(symSqrMat2, dim=1)))
    do iSpin = 1, nSpin
      symSqrMat2(:,:, iSpin) = transpose(deltaRhoOverlap(:,:, iSpin)) * deltaRhoOverlap(:,:, iSpin)
      call gemm(symSqrMat2Tmp, overlapSym, deltaRhoOverlap(:,:, iSpin))
      symSqrMat2(:,:, iSpin) = symSqrMat2(:,:, iSpin) + symSqrMat2Tmp * deltaRhoSqrSym(:,:, iSpin)
      ! symmetrize temporary storage
      symSqrMat2(:,:, iSpin) = 0.5_dp * (symSqrMat2(:,:, iSpin) + transpose(symSqrMat2(:,:, iSpin)))
    end do

    ! free some memory
    deallocate(overlapSym)
    deallocate(deltaRhoSqrSym)
    deallocate(symSqrMat2Tmp)
    deallocate(deltaRhoOverlap)

    allocate(overSqrPrime(3, nOrb, nOrb))
    allocate(tmpGradients, mold=gradients)
    tmpGradients(:,:) = 0.0_dp

    loopForceAtom: do iAtForce = 1, nAtom0
      call getUnpackedOverlapPrime_real(iAtForce, skOverCont, orb, derivator, symNeighbourList,&
          & nNeighbourCamSym, iSquare, this%rCoords, overSqrPrime)
      call getUnpackedCamGammaAOPrime(iAtForce, this%camdGammaEval0, iSquare, camdGammaAO)
      do iSpin = 1, nSpin
        do iCoord = 1, 3
          ! first term of Eq.(B5)
          tmpGradients(iCoord, iAtForce) = tmpGradients(iCoord, iAtForce)&
              & - sum(overSqrPrime(iCoord, :,:) * symSqrMat1(:,:, iSpin))
          ! second term of Eq.(B5)
          tmpGradients(iCoord, iAtForce) = tmpGradients(iCoord, iAtForce)&
              & - 0.5_dp * sum(camdGammaAO(:,:, iCoord) * symSqrMat2(:,:, iSpin))
        end do
      end do
    end do loopForceAtom

    if (this%tREKS) then
      gradients(:,:) = gradients + tmpGradients
    else
      gradients(:,:) = gradients + 0.5_dp * nSpin * tmpGradients
    end if

  end subroutine addCamGradientsMatrix_real


  !> Calculates derivative of the (long-range + HF full-range) gamma interaction w.r.t. given atom.
  !!
  !! 2nd term of Eq.(B5) of Phys. Rev. Materials 7, 063802 (DOI: 10.1103/PhysRevMaterials.7.063802)
  subroutine getUnpackedCamGammaAOPrime(iAtomPrime, camdGammaEval0, iSquare, camdGammaAO)

    !> Overlap derivative calculated w.r.t. this atom in the central cell
    integer, intent(in) :: iAtomPrime

    !> Pre-tabulated (long-range + HF full-range) gamma derivative of At1 and At2 (in central cell)
    real(dp), intent(in) :: camdGammaEval0(:,:,:)

    !> Position of each atom in the rows/columns of the square matrices. Shape: (nAtom)
    integer, intent(in) :: iSquare(:)

    !> CAM gamma matrix, including periodic images
    real(dp), intent(out) :: camdGammaAO(:,:,:)

    !! Atom interacting with iAtomPrime
    integer :: iAt2

    !! Iterates over coordinates
    integer :: iCoord

    camdGammaAO(:,:,:) = 0.0_dp

    do iAt2 = 1, size(camdGammaEval0, dim=3)
      do iCoord = 1, 3
        camdGammaAO(iSquare(iAtomPrime):iSquare(iAtomPrime + 1) - 1,&
            & iSquare(iAt2):iSquare(iAt2 + 1) - 1, iCoord)&
            & = camdGammaEval0(iCoord, iAtomPrime, iAt2)
      end do
    end do

  end subroutine getUnpackedCamGammaAOPrime


  !> Calculates the derivative of the square, dense, unpacked, Gamma-point overlap matrix w.r.t.
  !! given atom.
  !!
  !! 1st term of Eq.(B5) of Phys. Rev. Materials 7, 063802 (DOI: 10.1103/PhysRevMaterials.7.063802)
  subroutine getUnpackedOverlapPrime_real(iAtomPrime, skOverCont, orb, derivator, symNeighbourList,&
      & nNeighbourCamSym, iSquare, rCoords, overSqrPrime)

    !> Overlap derivative calculated w.r.t. this atom in the central cell
    integer, intent(in) :: iAtomPrime

    !> Sparse overlap container
    type(TSlakoCont), intent(in) :: skOverCont

    !> Orbital information.
    type(TOrbitals), intent(in) :: orb

    !> Differentiation object
    class(TNonSccDiff), intent(in) :: derivator

    !> List of neighbours for each atom (symmetric version)
    type(TSymNeighbourList), intent(in) :: symNeighbourList

    !> Nr. of neighbours for each atom.
    integer, intent(in) :: nNeighbourCamSym(:)

    !> Position of each atom in the rows/columns of the square matrices. Shape: (nAtom)
    integer, intent(in) :: iSquare(:)

    !> Atomic coordinates in absolute units, potentially including periodic images
    real(dp), intent(in) :: rCoords(:,:)

    !> Overlap derivative
    real(dp), intent(out) :: overSqrPrime(:,:,:)

    !! Temporary overlap derivative storage
    real(dp) :: overPrime(orb%mOrb, orb%mOrb, 3)

    !! Dense matrix descriptor indices
    integer, parameter :: descLen = 3, iStart = 1, iEnd = 2, iNOrb = 3

    !! Stores start/end index and number of orbitals of square matrices
    integer :: descAt1(descLen), descAt2(descLen)

    !! Atom to calculate energy gradient components for
    integer :: iAt2, iAt2fold, iNeigh

    !! Iterates over coordinates
    integer :: iCoord

    overSqrPrime(:,:,:) = 0.0_dp

    descAt1 = getDescriptor(iAtomPrime, iSquare)
    do iNeigh = 0, nNeighbourCamSym(iAtomPrime)
      iAt2 = symNeighbourList%neighbourList%iNeighbour(iNeigh, iAtomPrime)
      iAt2fold = symNeighbourList%img2CentCell(iAt2)
      if (iAtomPrime == iAt2fold) cycle
      descAt2 = getDescriptor(iAt2fold, iSquare)
      overPrime(:,:,:) = 0.0_dp
      call derivator%getFirstDeriv(overPrime, skOverCont, rCoords, symNeighbourList%species,&
            & iAtomPrime, iAt2, orb)
      do iCoord = 1, 3
        overSqrPrime(iCoord, descAt2(iStart):descAt2(iStart) + descAt2(iNOrb) - 1,&
            & descAt1(iStart):descAt1(iStart) + descAt1(iNOrb) - 1)&
            & = overSqrPrime(iCoord, descAt2(iStart):descAt2(iStart) + descAt2(iNOrb) - 1,&
            & descAt1(iStart):descAt1(iStart) + descAt1(iNOrb) - 1)&
            & + overPrime(1:descAt2(iNOrb), 1:descAt1(iNOrb), iCoord)
      end do
    end do

  end subroutine getUnpackedOverlapPrime_real

#:endif


  !> Interface routine to add gradients due to CAM range-separated contributions.
  !! (k-point version)
  subroutine addCamGradients_kpts(this, env, denseDesc, ints, orb, skOverCont, derivator,&
      & densityMatrix, neighbourList, nNeighbourSK, symNeighbourList, nNeighbourCamSym, cellVec,&
      & iCellVec, iSparseStart, img2CentCell, kPoints, kWeights, gradients, errStatus)

    !> Class instance
    class(THybridXcFunc), intent(inout), target :: this

    !> Environment settings
    type(TEnvironment), intent(inout) :: env

    !> Dense matrix descriptor
    type(TDenseDescr), intent(in) :: denseDesc

    !> Integral container
    type(TIntegral), intent(in) :: ints

    !> Orbital information.
    type(TOrbitals), intent(in) :: orb

    !> Sparse overlap container
    type(TSlakoCont), intent(in) :: skOverCont

    !> Differentiation object
    class(TNonSccDiff), intent(in) :: derivator

    !> Holds real and complex delta density matrices
    type(TDensityMatrix), intent(in) :: densityMatrix

    !> List of neighbours for each atom
    type(TNeighbourList), intent(in) :: neighbourList

    !> Number of neighbours for each of the atoms
    integer, intent(in) :: nNeighbourSK(:)

    !> List of neighbours for each atom (symmetric version)
    type(TSymNeighbourList), intent(in) :: symNeighbourList

    !> Nr. of neighbours for each atom.
    integer, intent(in) :: nNeighbourCamSym(:)

    !> Vectors to neighboring unit cells in relative units
    real(dp), intent(in) :: cellVec(:,:)

    !> Index for which unit cell atoms are associated with
    integer, intent(in) :: iCellVec(:)

    !> Index array for the start of atomic blocks in sparse arrays
    integer, intent(in) :: iSparseStart(:,:)

    !> Map from image atoms to the original unique atom
    integer, intent(in) :: img2CentCell(:)

    !> The k-points in relative coordinates to calculate delta H(k) for
    real(dp), intent(in) :: kPoints(:,:)

    !> The k-point weights (for energy contribution)
    real(dp), intent(in) :: kWeights(:)

    !> Energy gradients
    real(dp), intent(inout) :: gradients(:,:)

    !> Error status
    type(TStatus), intent(inout) :: errStatus

    select case(this%hybridXcAlg)
    case (hybridXcAlgo%thresholdBased)
      @:RAISE_ERROR(errStatus, -1, "Hybrid functionals don't yet support gradient calculations for&
          & periodic systems beyond the Gamma point using the thresholded HFX construction&
          & algorithm.")
    case (hybridXcAlgo%neighbourBased)
      if (this%gammaType == hybridXcGammaTypes%mic) then
        @:RAISE_ERROR(errStatus, -1, "Hybrid functionals don't yet support gradient calculations&
            & for periodic systems beyond the Gamma point for CoulombMatrix setting&
            & 'MinimumImage'.")
      else
        @:ASSERT(allocated(densityMatrix%deltaRhoInCplxHS)&
            & .and. allocated(densityMatrix%deltaRhoOutCplx))
        call addCamGradientsNeighbour_kpts_ct(this, densityMatrix%deltaRhoInCplxHS,&
            & densityMatrix%deltaRhoOutCplx, symNeighbourList, nNeighbourCamSym, cellVec,&
            & denseDesc%iAtomStart, orb, kPoints, kWeights, skOverCont, derivator, gradients)
      end if
    case (hybridXcAlgo%matrixBased)
      call addCamGradientsMatrix_kpts_ct(this, env, denseDesc, skOverCont, derivator, orb, ints,&
          & densityMatrix, neighbourList, nNeighbourSK, symNeighbourList, nNeighbourCamSym,&
          & iCellVec, cellVec, iSparseStart, img2CentCell, kPoints, kWeights, gradients)
    end select

  end subroutine addCamGradients_kpts


  !> Adds range-separated contributions to Hamiltonian, using matrix based algorithm.
  !! (k-point version)
  subroutine addCamGradientsMatrix_kpts_ct(this, env, denseDesc, skOverCont, derivator, orb, ints,&
      & densityMatrix, neighbourList, nNeighbourSK, symNeighbourList, nNeighbourCamSym, iCellVec,&
      & cellVec, iSparseStart, img2CentCell, kPoints, kWeights, gradients)

    !> Class instance
    class(THybridXcFunc), intent(inout) :: this

    !> Environment settings
    type(TEnvironment), intent(inout) :: env

    !> Dense matrix descriptor
    type(TDenseDescr), intent(in) :: denseDesc

    !> Sparse overlap container
    type(TSlakoCont), intent(in) :: skOverCont

    !> Differentiation object
    class(TNonSccDiff), intent(in) :: derivator

    !> Orbital information.
    type(TOrbitals), intent(in) :: orb

    !> Integral container
    type(TIntegral), intent(in) :: ints

    !> Holds real and complex delta density matrices
    type(TDensityMatrix), intent(in), target :: densityMatrix

    !> List of neighbours for each atom
    type(TNeighbourList), intent(in) :: neighbourList

    !> Number of neighbours for each of the atoms
    integer, intent(in) :: nNeighbourSK(:)

    !> List of neighbours for each atom (symmetric version)
    type(TSymNeighbourList), intent(in) :: symNeighbourList

    !> Nr. of neighbours for each atom.
    integer, intent(in) :: nNeighbourCamSym(:)

    !> Index for which unit cell atoms are associated with
    integer, intent(in) :: iCellVec(:)

    !> Vectors (in units of the lattice constants) to cells of the lattice
    real(dp), intent(in) :: cellVec(:,:)

    !> Index array for the start of atomic blocks in sparse arrays
    integer, intent(in) :: iSparseStart(:,:)

    !> Map from image atoms to the original unique atom
    integer, intent(in) :: img2CentCell(:)

    !> The k-points in relative coordinates to calculate delta H(k) for
    real(dp), intent(in), target :: kPoints(:,:)

    !> The k-point weights
    real(dp), intent(in), target :: kWeights(:)

    !> Energy gradients
    real(dp), intent(inout) :: gradients(:,:)

    !! Temporary energy gradients
    real(dp), allocatable :: tmpGradients(:,:)

    !! Temporary y-storage
    complex(dp), allocatable :: gammaAO(:,:), gammaAOCc(:,:)

    !! Temporary y-storage
    complex(dp), allocatable :: dGammaAO(:,:,:), dGammaAOCc(:,:,:)

    !! Temporary storage for derivative of S(k')
    complex(dp), allocatable :: overSqrPrime(:,:,:)

    !! Size of square matrices (e.g. Hamiltonian)
    integer :: squareSize

    !! Number of k-point-spin compound indices
    integer :: nKS

    !! K-point indices
    integer :: iK, iKPrime

    !! Spin index
    integer :: iS

    !! Number of k-points and spins
    integer :: nK, nS

    !! Global iKS for arrays present at every MPI rank
    integer :: iGlobalKS, iGlobalKPrimeS

    !! Index of atom to calculate derivative for
    integer :: iAtForce

    !! Iterates over coordinates
    integer :: iCoord

    !! Number of atoms in central cell
    integer :: nAtom0

    !! Temporary storage
    complex(dp), allocatable :: tmp(:,:)
    complex(dp), dimension(:,:), allocatable :: dPm, dP_S_cc, dP_S, dPp_Sp, dPp_Sp_T, Sp_dPp_Sp_T

    !! Composite index for two nested loops over central cell
    integer, allocatable :: iKSComposite(:,:)

    !! Dense, square, k-space overlap matrices S(k) on all MPI processes
    complex(dp), allocatable :: SSqrCplx(:,:,:)

    !! Auxiliary variables
    real(dp) :: wk_wkp
    integer :: ii

    nAtom0 = size(this%species0)
    squareSize = size(densityMatrix%deltaRhoOutCplx, dim=1)
    nS = size(densityMatrix%iKiSToiGlobalKS, dim=2)
    nK = size(kPoints, dim=2)
    nKS = nK * nS

    call getDenseSqrDualSpaceOverlap(env, denseDesc, ints, neighbourList, nNeighbourSK, iCellVec,&
        & cellVec, iSparseStart, img2CentCell, kPoints, SSqrCplx)

    allocate(tmpGradients, mold=gradients)
    tmpGradients(:,:) = 0.0_dp

    allocate(tmp(squareSize, squareSize))

    ! dP(k)* = dP(-k)
    allocate(dPm(squareSize, squareSize))
    ! [dP(k)@S(k)]* = dP(-k)@S(-k)
    allocate(dP_S_cc(squareSize, squareSize))

    ! dP(k)@S(k)
    allocate(dP_S(squareSize, squareSize))
    ! dP(k')@S(k')
    allocate(dPp_Sp(squareSize, squareSize))

    ! [dP(k')@S(k')]^T
    allocate(dPp_Sp_T(squareSize, squareSize))
    ! [S(k')@dP(k')@S(k')]^T = [S(k')@dP(k')@S(k')]*
    allocate(Sp_dPp_Sp_T(squareSize, squareSize))

    ! y(k'-k)
    allocate(gammaAO(squareSize, squareSize))
    ! y(k'+k)
    allocate(gammaAOCc(squareSize, squareSize))

    ! dy(k'-k)
    allocate(dGammaAO(squareSize, squareSize, 3))
    ! dy(k'+k)
    allocate(dGammaAOCc(squareSize, squareSize, 3))

    ! dS(k')
    allocate(overSqrPrime(3, squareSize, squareSize))

    ! this is why the number of MPI groups may not exceed nS * nK**2 processes
    ! (there wouldn't be any additional speedup)
    call getiKSiKSPrimeCompositeIndex(env, nS, nK, nK, iKSComposite)

    loopForceAtom: do iAtForce = 1, nAtom0
      do ii = 1, size(iKSComposite, dim=2)
        iS = iKSComposite(1, ii)
        iK = iKSComposite(2, ii)
        iGlobalKS = densityMatrix%iKiSToiGlobalKS(iK, iS)

        iKPrime = iKSComposite(3, ii)
        iGlobalKPrimeS = densityMatrix%iKiSToiGlobalKS(iKPrime, iS)

        wk_wkp = kWeights(iK) * kWeights(iKPrime)

        ! dP(k)@S(k)
        call hemm(dP_S, 'l', densityMatrix%deltaRhoOutCplx(:,:, iGlobalKS), SSqrCplx(:,:, iK))
        dP_S_cc(:,:) = conjg(dP_S)

        ! dP(k')@S(k')
        call hemm(dPp_Sp, 'l', densityMatrix%deltaRhoOutCplx(:,:, iGlobalKPrimeS),&
            & SSqrCplx(:,:, iKPrime))

        ! [dP(k')@S(k')]^T
        dPp_Sp_T(:,:) = transpose(dPp_Sp)

        ! [S(k')@dP(k')@S(k')]^T = [S(k')@dP(k')@S(k')]*
        call hemm(tmp, 'l', SSqrCplx(:,:, iKPrime), dPp_Sp)
        Sp_dPp_Sp_T(:,:) = conjg(tmp)

        ! dP(-k)
        dPm(:,:) = conjg(densityMatrix%deltaRhoOutCplx(:,:, iGlobalKS))

        ! spin-polarized case: y-matrix constructed even if k and k' do not change
        ! y(k'-k)
        call getCamGammaFourierAO(this, denseDesc%iAtomStart, this%cellVecsG, this%rCellVecsG,&
            & -kPoints(:, iK), -kPoints(:, iKPrime), gammaAO)
        ! y(k'+k)
        call getCamGammaFourierAO(this, denseDesc%iAtomStart, this%cellVecsG, this%rCellVecsG,&
            & kPoints(:, iK), -kPoints(:, iKPrime), gammaAOCc)
        ! dy(k'-k)
        call getCamGammaFourierAOPrime(this, iAtForce, denseDesc%iAtomStart, this%cellVecsG,&
            & this%rCellVecsG, -kPoints(:, iK), -kPoints(:, iKPrime), dGammaAO)
        ! dy(k'+k)
        call getCamGammaFourierAOPrime(this, iAtForce, denseDesc%iAtomStart, this%cellVecsG,&
            & this%rCellVecsG, kPoints(:, iK), -kPoints(:, iKPrime), dGammaAOCc)

        ! dS(k')
        call getUnpackedOverlapPrime_kpts(iAtForce, skOverCont, orb, derivator, symNeighbourList,&
            & nNeighbourCamSym, denseDesc%iAtomStart, cellVec, this%rCoords, kPoints(:, iKPrime),&
            & overSqrPrime)

        ! Term 1
        call hemm(tmp, 'r', densityMatrix%deltaRhoOutCplx(:,:, iGlobalKS) * gammaAO, dPp_Sp)
        call hemm(tmp, 'r', densityMatrix%deltaRhoOutCplx(:,:, iGlobalKPrimeS), dP_S * gammaAO,&
            & beta=(1.0_dp, 0.0_dp))
        ! Term 1 (complex conjugated for inverse k-points, k -> -k)
        call hemm(tmp, 'r', dPm * gammaAOCc, dPp_Sp, beta=(1.0_dp, 0.0_dp))
        call hemm(tmp, 'r', densityMatrix%deltaRhoOutCplx(:,:, iGlobalKPrimeS),&
            & dP_S_cc * gammaAOCc, beta=(1.0_dp, 0.0_dp))
        ! [...]^adj = 0.5 * [(...) + (...)^adj]
        ! we actually calculate ([...]^adj)^T here, so that the transpose in dS(k') * tmp^T can be
        ! omitted when adding the contribution to tmpGradients below
        ! (extra factor of 0.5 accounts for the additional -k points)
        tmp(:,:) = 0.25_dp * (transpose(tmp) + conjg(tmp))

        ! Add term 1
        ! (sum should run over dS(k') * tmp^T, the transpose is already considered in the
        ! hermitianization above)
        do iCoord = 1, 3
          tmpGradients(iCoord, iAtForce) = tmpGradients(iCoord, iAtForce)&
              & - wk_wkp * sum(overSqrPrime(iCoord, :,:) * tmp)
        end do

        ! Term 2
        tmp(:,:) = dP_S * dPp_Sp_T + Sp_dPp_Sp_T * densityMatrix%deltaRhoOutCplx(:,:, iGlobalKS)
        ! [...]^adj = 0.5 * [(...) + (...)^adj]
        ! we actually calculate ([...]^adj)^T here, so that the transpose in dS(k') * tmp^T can be
        ! omitted when adding the contribution to tmpGradients below
        ! (extra factor of 0.5 accounts for the additional -k points)
        tmp(:,:) = 0.25_dp * (transpose(tmp) + conjg(tmp))

        ! Add term 2
        ! (sum should run over dS(k') * tmp^T, the transpose is already considered in the
        ! hermitianization above)
        do iCoord = 1, 3
          tmpGradients(iCoord, iAtForce) = tmpGradients(iCoord, iAtForce)&
              & - 0.5_dp * wk_wkp * sum(dGammaAO(:,:, iCoord) * tmp)
        end do

        ! Term 2 (complex conjugated for inverse k-points, k -> -k)
        tmp(:,:) = dP_S_cc * dPp_Sp_T + Sp_dPp_Sp_T * dPm
        ! [...]^adj = 0.5 * [(...) + (...)^adj]
        ! we actually calculate ([...]^adj)^T here, so that the transpose in dS(k') * tmp^T can be
        ! omitted when adding the contribution to tmpGradients below
        ! (extra factor of 0.5 accounts for the additional -k points)
        tmp(:,:) = 0.25_dp * (transpose(tmp) + conjg(tmp))

        ! Add term 2
        ! (sum should run over dS(k') * tmp^T, the transpose is already considered in the
        ! hermitianization above)
        do iCoord = 1, 3
          tmpGradients(iCoord, iAtForce) = tmpGradients(iCoord, iAtForce)&
              & - 0.5_dp * wk_wkp * sum(dGammaAOCc(:,:, iCoord) * tmp)
        end do
      end do
    end do loopForceAtom

  #:if WITH_MPI
    call mpifx_allreduceip(env%mpi%globalComm, tmpGradients, MPI_SUM)
  #:endif

    if (this%tREKS) then
      gradients(:,:) = gradients + tmpGradients
    else
      gradients(:,:) = gradients + 0.5_dp * nS * tmpGradients
    end if

  end subroutine addCamGradientsMatrix_kpts_ct


  !> Adds CAM gradients due to CAM range-separated contributions (k-point version).
  subroutine addCamGradientsNeighbour_kpts_ct(this, deltaRhoSqr, deltaRhoOutSqrCplx,&
      & symNeighbourList, nNeighbourCamSym, cellVecs, iSquare, orb, kPoints, kWeights, skOverCont,&
      & derivator, gradients)

    !> Class instance
    class(THybridXcFunc), intent(inout), target :: this

    !> Square (unpacked) delta spin-density matrix at BvK real-space shifts
    real(dp), intent(in) :: deltaRhoSqr(:,:,:,:,:,:)

    !> Square (unpacked) delta spin-density matrix in k-space
    complex(dp), intent(in) :: deltaRhoOutSqrCplx(:,:,:)

    !> List of neighbours for each atom (symmetric version)
    type(TSymNeighbourList), intent(in) :: symNeighbourList

    !> Nr. of neighbours for each atom.
    integer, intent(in) :: nNeighbourCamSym(:)

    !> Vectors to neighboring unit cells in relative units
    real(dp), intent(in) :: cellVecs(:,:)

    !> Position of each atom in the rows/columns of the square matrices. Shape: (nAtom)
    integer, intent(in) :: iSquare(:)

    !> Orbital information.
    type(TOrbitals), intent(in) :: orb

    !> The k-points in relative coordinates to calculate delta H(k) for
    real(dp), intent(in) :: kPoints(:,:)

    !> The k-point weights (for energy contribution)
    real(dp), intent(in) :: kWeights(:)

    !> Sparse overlap container
    type(TSlakoCont), intent(in) :: skOverCont

    !> Differentiation object
    class(TNonSccDiff), intent(in) :: derivator

    !> Energy gradients
    real(dp), intent(inout) :: gradients(:,:)

    !! Dense matrix descriptor indices
    integer, parameter :: descLen = 3, iStart = 1, iEnd = 2, iNOrb = 3

    !! Atom blocks from sparse, real-space overlap matrices
    real(dp), pointer :: pSak(:,:), pSam(:,:), pSbk(:,:), pSbn(:,:)

    !! Stores start/end index and number of orbitals of square matrices
    integer :: descA(descLen), descB(descLen), descM(descLen), descN(descLen), descK(descLen)

    !! Temporary storages
    real(dp), allocatable :: tmpGradients(:,:)

    !! Overlap derivatives
    real(dp), dimension(orb%mOrb, orb%mOrb, 3) :: SbnPrimeKequalsN

    !! Number of atoms in central cell
    integer :: nAtom0

    !! Real-space \vec{h} and \vec{l} vectors in relative coordinates
    real(dp) :: vecH(3), vecL(3)

    !! K-point index
    integer :: iK

    !! Spin index
    integer :: iS

    !! Atom indices (central cell)
    integer :: iAtM, iAtN, iAtK

    !! Neighbour indices (+corresponding atom indices)
    integer :: iNeighK, iNeighM, iNeighN, iAtA, iAtB

    !! Folded (to central cell) atom indices
    integer :: iAtAfold, iAtBfold

    !! Auxiliary variables for setting up 2D pointer to sparse overlap
    integer :: ind, nOrbAt, nOrbNeigh

    !! Integer BvK index
    integer :: bvKIndex(3)

    !! Phase factor
    complex(dp) :: phase

    !! Iterates over all BvK real-space vectors
    integer :: iG, iGKB, iGMK, iGMB, iGAK, iGAB

    !! Orbital indices
    integer :: mu, nu, kk, alpha, beta

    !! Number of k-points and spins
    integer :: nK, nS

    !! Overlap matrix elements
    real(dp) :: Sam, Sak, Sbk, Sbn

    !! Density matrix elements
    real(dp) :: dPab
    complex(dp) :: dPkm, dPnk, dPnm

    !! Products phase * gammaMK, phase * gammaMB, phase * gammaAK, phase * gammaAB
    complex(dp) :: phaseGammaMK, phaseGammaMB, phaseGammaAK, phaseGammaAB

    !! Product dPkm * phase * gammaMK
    complex(dp) :: dPkmPhaseGammaMK

    !! Product dPkm * phase * gammaMB
    complex(dp) :: dPkmPhaseGammaMB

    !! Product dPkm * phase * gammaAK
    complex(dp) :: dPkmPhaseGammaAK

    !! Product dPkm * phase * gammaAB
    complex(dp) :: dPkmPhaseGammaAB

    !! Product dPkm * phase * gammaMK * Sam
    complex(dp) :: dPkmPhaseGammaMKSam

    !! Product dPkm * phase * gammaMB * Sam
    complex(dp) :: dPkmPhaseGammaMBSam

    !! Product dPkm * phase * gammaAK * Sam
    complex(dp) :: dPkmPhaseGammaAKSam

    !! Product dPkm * phase * gammaAB * Sam
    complex(dp) :: dPkmPhaseGammaABSam

    !! Directional derivatives
    real(dp), allocatable :: dGammaMK(:,:), dGammaAK(:,:), dGammaKB(:,:)

    !! Product phase * dGammaMK(:, iG)
    complex(dp) :: phasedGammaMK(3), phasedGammaKB(3), phasedGammaAK(3)

    complex(dp) :: dPkmPhasedGammaMK(3), dPkmPhasedGammaMKSam(3), dPnkPhasedGammaKB(3)
    complex(dp) :: dPnkPhasedGammaKBSak(3), dPkmPhasedGammaAK(3), dPkmPhasedGammaAKSam(3)
    complex(dp) :: dPnmPhasedGammaAK(3), dPnmPhasedGammaAKSam(3)

    !! Composite index for mapping iK/iS --> iGlobalKS
    integer, allocatable :: iKiSToiGlobalKS(:,:)
    integer :: iGlobalKS

    nAtom0 = size(this%species0)
    nS = size(deltaRhoSqr, dim=6)
    nK = size(kPoints, dim=2)

    ! Build spin/k-point composite index for all spins and k-points (global)
    iGlobalKS = 1
    allocate(iKiSToiGlobalKS(nK, nS))
    do iS = 1, nS
      do iK = 1, nK
        iKiSToiGlobalKS(iK, iS) = iGlobalKS
        iGlobalKS = iGlobalKS + 1
      end do
    end do

    ! allocate gradient contribution
    allocate(tmpGradients(3, size(gradients, dim=2)))
    tmpGradients(:,:) = 0.0_dp

    ! First terms of (48)
    loopK1: do iAtK = 1, nAtom0
      descK = getDescriptor(iAtK, iSquare)
      loopM1: do iAtM = 1, nAtom0
        descM = getDescriptor(iAtM, iSquare)
        loopB1: do iNeighK = 0, nNeighbourCamSym(iAtK)
          iAtB = symNeighbourList%neighbourList%iNeighbour(iNeighK, iAtK)
          iAtBfold = symNeighbourList%img2CentCell(iAtB)
          if (iAtBfold == iAtK) cycle
          descB = getDescriptor(iAtBfold, iSquare)
          vecL(:) = cellVecs(:, symNeighbourList%iCellVec(iAtB))
          SbnPrimeKequalsN(:,:,:) = 0.0_dp
          call derivator%getFirstDeriv(SbnPrimeKequalsN, skOverCont, this%rCoords,&
              & symNeighbourList%species, iAtK, iAtB, orb)
          loopA1: do iNeighM = 0, nNeighbourCamSym(iAtM)
            iAtA = symNeighbourList%neighbourList%iNeighbour(iNeighM, iAtM)
            iAtAfold = symNeighbourList%img2CentCell(iAtA)
            descA = getDescriptor(iAtAfold, iSquare)
            vecH(:) = cellVecs(:, symNeighbourList%iCellVec(iAtA))
            ! get 2D pointer to Sam overlap block
            ind = symNeighbourList%iPair(iNeighM, iAtM) + 1
            nOrbAt = descM(iNOrb)
            nOrbNeigh = descA(iNOrb)
            pSam(1:nOrbNeigh, 1:nOrbAt) => this%overSym(ind:ind + nOrbNeigh * nOrbAt - 1)

            loopGMK1: do iG = 1, size(this%nNonZeroGammaG(iAtM, iAtK)%data)
              iGMK = this%nNonZeroGammaG(iAtM, iAtK)%data(iG)
              bvKIndex(:) = this%foldToBvKIndex(vecH - vecL - this%cellVecsG(:, iGMK))

              loopKptsMK1: do iK = 1, nK
                phase = exp(imag * dot_product(2.0_dp * pi * kPoints(:, iK),&
                    & -this%cellVecsG(:, iGMK)))
                phaseGammaMK = phase * cmplx(this%camGammaEvalG(iAtM, iAtK)%data(iG)&
                    & * kWeights(iK), 0, dp)

                ! term #1
                do mu = 1, descM(iNOrb)
                  do kk = 1, descK(iNOrb)
                    do iS = 1, nS
                      dPkm = deltaRhoOutSqrCplx(descK(iStart) + kk - 1, descM(iStart) + mu - 1,&
                          & iKiSToiGlobalKS(iK, iS))
                      dPkmPhaseGammaMK = dPkm * phaseGammaMK
                      do alpha = 1, descA(iNOrb)
                        Sam = pSam(alpha, mu)
                        dPkmPhaseGammaMKSam = dPkmPhaseGammaMK * cmplx(Sam, 0, dp)
                        do beta = 1, descB(iNOrb)
                          dPab = deltaRhoSqr(descA(iStart) + alpha - 1, descB(iStart) + beta - 1,&
                              & bvKIndex(1), bvKIndex(2), bvKIndex(3), iS)

                          tmpGradients(:, iAtK) = tmpGradients(:, iAtK)&
                              & + real(2.0_dp * dPkmPhaseGammaMKSam&
                              & * cmplx(dPab * SbnPrimeKequalsN(beta, kk, :), 0, dp), dp)
                        end do
                      end do
                    end do
                  end do
                end do

              end do loopKptsMK1
            end do loopGMK1

            loopGMB: do iG = 1, size(this%nNonZeroGammaG(iAtM, iAtBfold)%data)
              iGMB = this%nNonZeroGammaG(iAtM, iAtBfold)%data(iG)
              bvKIndex(:) = this%foldToBvKIndex(vecH - this%cellVecsG(:, iGMB))

              loopKptsMB: do iK = 1, nK
                phase = exp(imag * dot_product(2.0_dp * pi * kPoints(:, iK),&
                    & vecL - this%cellVecsG(:, iGMB)))
                phaseGammaMB = phase * cmplx(this%camGammaEvalG(iAtM, iAtBfold)%data(iG)&
                    & * kWeights(iK), 0, dp)

                ! term #2
                do mu = 1, descM(iNOrb)
                  do kk = 1, descK(iNOrb)
                    do iS = 1, nS
                      dPkm = deltaRhoOutSqrCplx(descK(iStart) + kk - 1, descM(iStart) + mu - 1,&
                          & iKiSToiGlobalKS(iK, iS))
                      dPkmPhaseGammaMB = dPkm * phaseGammaMB
                      do alpha = 1, descA(iNOrb)
                        Sam = pSam(alpha, mu)
                        dPkmPhaseGammaMBSam = dPkmPhaseGammaMB * cmplx(Sam, 0, dp)
                        do beta = 1, descB(iNOrb)
                          dPab = deltaRhoSqr(descA(iStart) + alpha - 1, descB(iStart) + beta - 1,&
                              & bvKIndex(1), bvKIndex(2), bvKIndex(3), iS)

                          tmpGradients(:, iAtK) = tmpGradients(:, iAtK)&
                              & + real(2.0_dp * dPkmPhaseGammaMBSam&
                              & * cmplx(dPab * SbnPrimeKequalsN(beta, kk, :), 0, dp), dp)
                        end do
                      end do
                    end do
                  end do
                end do

              end do loopKptsMB
            end do loopGMB

            loopGAK1: do iG = 1, size(this%nNonZeroGammaG(iAtAfold, iAtK)%data)
              iGAK = this%nNonZeroGammaG(iAtAfold, iAtK)%data(iG)
              bvKIndex(:) = this%foldToBvKIndex(-this%cellVecsG(:, iGAK) - vecL)

              loopKptsAK1: do iK = 1, nK
                phase = exp(imag * dot_product(2.0_dp * pi * kPoints(:, iK),&
                    & -this%cellVecsG(:, iGAK) - vecH))
                phaseGammaAK = phase * cmplx(this%camGammaEvalG(iAtAfold, iAtK)%data(iG)&
                    & * kWeights(iK), 0, dp)

                ! term #3
                do mu = 1, descM(iNOrb)
                  do kk = 1, descK(iNOrb)
                    do iS = 1, nS
                      dPkm = deltaRhoOutSqrCplx(descK(iStart) + kk - 1, descM(iStart) + mu - 1,&
                          & iKiSToiGlobalKS(iK, iS))
                      dPkmPhaseGammaAK = dPkm * phaseGammaAK
                      do alpha = 1, descA(iNOrb)
                        Sam = pSam(alpha, mu)
                        dPkmPhaseGammaAKSam = dPkmPhaseGammaAK * cmplx(Sam, 0, dp)
                        do beta = 1, descB(iNOrb)
                          dPab = deltaRhoSqr(descA(iStart) + alpha - 1, descB(iStart) + beta - 1,&
                              & bvKIndex(1), bvKIndex(2), bvKIndex(3), iS)

                          tmpGradients(:, iAtK) = tmpGradients(:, iAtK)&
                              & + real(2.0_dp * dPkmPhaseGammaAKSam&
                              & * cmplx(dPab * SbnPrimeKequalsN(beta, kk, :), 0, dp), dp)
                        end do
                      end do
                    end do
                  end do
                end do

              end do loopKptsAK1
            end do loopGAK1

            loopGAB1: do iG = 1, size(this%nNonZeroGammaG(iAtAfold, iAtBfold)%data)
              iGAB = this%nNonZeroGammaG(iAtAfold, iAtBfold)%data(iG)
              bvKIndex(:) = this%foldToBvKIndex(-this%cellVecsG(:, iGAB))

              loopKptsAB1: do iK = 1, nK
                phase = exp(imag * dot_product(2.0_dp * pi * kPoints(:, iK),&
                    & -this%cellVecsG(:, iGAB) + vecL - vecH))
                phaseGammaAB = phase * cmplx(this%camGammaEvalG(iAtAfold, iAtBfold)%data(iG)&
                    & * kWeights(iK), 0, dp)

                ! term #4
                do mu = 1, descM(iNOrb)
                  do kk = 1, descK(iNOrb)
                    do iS = 1, nS
                      dPkm = deltaRhoOutSqrCplx(descK(iStart) + kk - 1, descM(iStart) + mu - 1,&
                          & iKiSToiGlobalKS(iK, iS))
                      dPkmPhaseGammaAB = dPkm * phaseGammaAB
                      do alpha = 1, descA(iNOrb)
                        Sam = pSam(alpha, mu)
                        dPkmPhaseGammaABSam = dPkmPhaseGammaAB * cmplx(Sam, 0, dp)
                        do beta = 1, descB(iNOrb)
                          dPab = deltaRhoSqr(descA(iStart) + alpha - 1, descB(iStart) + beta - 1,&
                              & bvKIndex(1), bvKIndex(2), bvKIndex(3), iS)

                          tmpGradients(:, iAtK) = tmpGradients(:, iAtK)&
                              & + real(2.0_dp * dPkmPhaseGammaABSam&
                              & * cmplx(dPab * SbnPrimeKequalsN(beta, kk, :), 0, dp), dp)
                        end do
                      end do
                    end do
                  end do
                end do

              end do loopKptsAB1
            end do loopGAB1

          end do loopA1
        end do loopB1
      end do loopM1
    end do loopK1

    ! term #1 with prefactor 1/8 of eq.(48)
    loopK2: do iAtK = 1, nAtom0
      descK = getDescriptor(iAtK, iSquare)
      loopM2: do iAtM = 1, nAtom0
        if (iAtK == iAtM) cycle
        descM = getDescriptor(iAtM, iSquare)
        dGammaMK = -this%camdGammaEvalG(iAtM, iAtK)%data
        loopB2: do iNeighK = 0, nNeighbourCamSym(iAtK)
          iAtB = symNeighbourList%neighbourList%iNeighbour(iNeighK, iAtK)
          iAtBfold = symNeighbourList%img2CentCell(iAtB)
          descB = getDescriptor(iAtBfold, iSquare)
          vecL(:) = cellVecs(:, symNeighbourList%iCellVec(iAtB))
          ! get 2D pointer to Sbk overlap block
          ind = symNeighbourList%iPair(iNeighK, iAtK) + 1
          nOrbAt = descK(iNOrb)
          nOrbNeigh = descB(iNOrb)
          pSbk(1:nOrbNeigh, 1:nOrbAt) => this%overSym(ind:ind + nOrbNeigh * nOrbAt - 1)
          loopA2: do iNeighM = 0, nNeighbourCamSym(iAtM)
            iAtA = symNeighbourList%neighbourList%iNeighbour(iNeighM, iAtM)
            iAtAfold = symNeighbourList%img2CentCell(iAtA)
            descA = getDescriptor(iAtAfold, iSquare)
            vecH(:) = cellVecs(:, symNeighbourList%iCellVec(iAtA))
            ! get 2D pointer to Sam overlap block
            ind = symNeighbourList%iPair(iNeighM, iAtM) + 1
            nOrbAt = descM(iNOrb)
            nOrbNeigh = descA(iNOrb)
            pSam(1:nOrbNeigh, 1:nOrbAt) => this%overSym(ind:ind + nOrbNeigh * nOrbAt - 1)

            loopGMK: do iG = 1, size(this%nNonZeroGammaG(iAtM, iAtK)%data)
              iGMK = this%nNonZeroGammaG(iAtM, iAtK)%data(iG)
              bvKIndex(:) = this%foldToBvKIndex(vecH - vecL - this%cellVecsG(:, iGMK))

              loopKptsMK: do iK = 1, nK
                phase = exp(imag * dot_product(2.0_dp * pi * kPoints(:, iK),&
                    & -this%cellVecsG(:, iGMK)))
                phasedGammaMK(:) = phase * cmplx(dGammaMK(:, iG) * kWeights(iK), 0, dp)

                do mu = 1, descM(iNOrb)
                  do kk = 1, descK(iNOrb)
                    do iS = 1, nS
                      dPkm = deltaRhoOutSqrCplx(descK(iStart) + kk - 1, descM(iStart) + mu - 1,&
                          & iKiSToiGlobalKS(iK, iS))
                      dPkmPhasedGammaMK(:) = dPkm * phasedGammaMK
                      do alpha = 1, descA(iNOrb)
                        Sam = pSam(alpha, mu)
                        dPkmPhasedGammaMKSam(:) = dPkmPhasedGammaMK * cmplx(Sam, 0, dp)
                        do beta = 1, descB(iNOrb)
                          Sbk = pSbk(beta, kk)
                          dPab = deltaRhoSqr(descA(iStart) + alpha - 1, descB(iStart) + beta - 1,&
                              & bvKIndex(1), bvKIndex(2), bvKIndex(3), iS)

                          tmpGradients(:, iAtK) = tmpGradients(:, iAtK)&
                              & + real(dPkmPhasedGammaMKSam * cmplx(Sbk * dPab, 0, dp), dp)
                        end do
                      end do
                    end do
                  end do
                end do

              end do loopKptsMK
            end do loopGMK

          end do loopA2
        end do loopB2
      end do loopM2
    end do loopK2

    ! term #2 with prefactor 1/8 of eq.(48)
    loopK3: do iAtK = 1, nAtom0
      descK = getDescriptor(iAtK, iSquare)
      loopN1: do iAtN = 1, nAtom0
        descN = getDescriptor(iAtN, iSquare)
        loopB3: do iNeighN = 0, nNeighbourCamSym(iAtN)
          iAtB = symNeighbourList%neighbourList%iNeighbour(iNeighN, iAtN)
          iAtBfold = symNeighbourList%img2CentCell(iAtB)
          if (iAtK == iAtBfold) cycle
          descB = getDescriptor(iAtBfold, iSquare)
          vecL(:) = cellVecs(:, symNeighbourList%iCellVec(iAtB))
          dGammaKB = this%camdGammaEvalG(iAtK, iAtBfold)%data
          ! get 2D pointer to Sbn overlap block
          ind = symNeighbourList%iPair(iNeighN, iAtN) + 1
          nOrbAt = descN(iNOrb)
          nOrbNeigh = descB(iNOrb)
          pSbn(1:nOrbNeigh, 1:nOrbAt) => this%overSym(ind:ind + nOrbNeigh * nOrbAt - 1)
          loopA3: do iNeighK = 0, nNeighbourCamSym(iAtK)
            iAtA = symNeighbourList%neighbourList%iNeighbour(iNeighK, iAtK)
            iAtAfold = symNeighbourList%img2CentCell(iAtA)
            descA = getDescriptor(iAtAfold, iSquare)
            vecH(:) = cellVecs(:, symNeighbourList%iCellVec(iAtA))
            ! get 2D pointer to Sak overlap block
            ind = symNeighbourList%iPair(iNeighK, iAtK) + 1
            nOrbAt = descK(iNOrb)
            nOrbNeigh = descA(iNOrb)
            pSak(1:nOrbNeigh, 1:nOrbAt) => this%overSym(ind:ind + nOrbNeigh * nOrbAt - 1)

            loopGKB: do iG = 1, size(this%nNonZeroGammaG(iAtK, iAtBfold)%data)
              iGKB = this%nNonZeroGammaG(iAtK, iAtBfold)%data(iG)
              bvKIndex(:) = this%foldToBvKIndex(vecH - this%cellVecsG(:, iGKB))

              loopKptsKB: do iK = 1, nK
                phase = exp(imag * dot_product(2.0_dp * pi * kPoints(:, iK),&
                    & vecL - this%cellVecsG(:, iGKB)))
                phasedGammaKB(:) = phase * cmplx(dGammaKB(:, iG) * kWeights(iK), 0, dp)

                do nu = 1, descN(iNOrb)
                  do kk = 1, descK(iNOrb)
                    do iS = 1, nS
                      dPnk = deltaRhoOutSqrCplx(descN(iStart) + nu - 1, descK(iStart) + kk - 1,&
                          & iKiSToiGlobalKS(iK, iS))
                      dPnkPhasedGammaKB(:) = dPnk * phasedGammaKB
                      do alpha = 1, descA(iNOrb)
                        Sak = pSak(alpha, kk)
                        dPnkPhasedGammaKBSak(:) = dPnkPhasedGammaKB * cmplx(Sak, 0, dp)
                        do beta = 1, descB(iNOrb)
                          Sbn = pSbn(beta, nu)
                          dPab = deltaRhoSqr(descA(iStart) + alpha - 1, descB(iStart) + beta - 1,&
                              & bvKIndex(1), bvKIndex(2), bvKIndex(3), iS)

                          tmpGradients(:, iAtK) = tmpGradients(:, iAtK)&
                              & + real(dPnkPhasedGammaKBSak * cmplx(Sbn * dPab, 0, dp), dp)
                        end do
                      end do
                    end do
                  end do
                end do

              end do loopKptsKB
            end do loopGKB

          end do loopA3
        end do loopB3
      end do loopN1
    end do loopK3

    ! term #3 with prefactor 1/8 of eq.(48)
    loopK4: do iAtK = 1, nAtom0
      descK = getDescriptor(iAtK, iSquare)
      loopM3: do iAtM = 1, nAtom0
        descM = getDescriptor(iAtM, iSquare)
        loopB4: do iNeighK = 0, nNeighbourCamSym(iAtK)
          iAtB = symNeighbourList%neighbourList%iNeighbour(iNeighK, iAtK)
          iAtBfold = symNeighbourList%img2CentCell(iAtB)
          descB = getDescriptor(iAtBfold, iSquare)
          vecL(:) = cellVecs(:, symNeighbourList%iCellVec(iAtB))
          ! get 2D pointer to Sbk overlap block
          ind = symNeighbourList%iPair(iNeighK, iAtK) + 1
          nOrbAt = descK(iNOrb)
          nOrbNeigh = descB(iNOrb)
          pSbk(1:nOrbNeigh, 1:nOrbAt) => this%overSym(ind:ind + nOrbNeigh * nOrbAt - 1)
          loopA4: do iNeighM = 0, nNeighbourCamSym(iAtM)
            iAtA = symNeighbourList%neighbourList%iNeighbour(iNeighM, iAtM)
            iAtAfold = symNeighbourList%img2CentCell(iAtA)
            if (iAtK == iAtAfold) cycle
            descA = getDescriptor(iAtAfold, iSquare)
            vecH(:) = cellVecs(:, symNeighbourList%iCellVec(iAtA))
            dGammaAK = -this%camdGammaEvalG(iAtAfold, iAtK)%data
            ! get 2D pointer to Sam overlap block
            ind = symNeighbourList%iPair(iNeighM, iAtM) + 1
            nOrbAt = descM(iNOrb)
            nOrbNeigh = descA(iNOrb)
            pSam(1:nOrbNeigh, 1:nOrbAt) => this%overSym(ind:ind + nOrbNeigh * nOrbAt - 1)

            loopGAK2: do iG = 1, size(this%nNonZeroGammaG(iAtAfold, iAtK)%data)
              iGAK = this%nNonZeroGammaG(iAtAfold, iAtK)%data(iG)
              bvKIndex(:) = this%foldToBvKIndex(-this%cellVecsG(:, iGAK) - vecL)

              loopKptsAK2: do iK = 1, nK
                phase = exp(imag * dot_product(2.0_dp * pi * kPoints(:, iK),&
                    & -this%cellVecsG(:, iGAK) - vecH))
                phasedGammaAK(:) = phase * cmplx(dGammaAK(:, iG) * kWeights(iK), 0, dp)

                do mu = 1, descM(iNOrb)
                  do kk = 1, descK(iNOrb)
                    do iS = 1, nS
                      dPkm = deltaRhoOutSqrCplx(descK(iStart) + kk - 1, descM(iStart) + mu - 1,&
                          & iKiSToiGlobalKS(iK, iS))
                      dPkmPhasedGammaAK(:) = dPkm * phasedGammaAK
                      do alpha = 1, descA(iNOrb)
                        Sam = pSam(alpha, mu)
                        dPkmPhasedGammaAKSam(:) = dPkmPhasedGammaAK * cmplx(Sam, 0, dp)
                        do beta = 1, descB(iNOrb)
                          Sbk = pSbk(beta, kk)
                          dPab = deltaRhoSqr(descA(iStart) + alpha - 1, descB(iStart) + beta - 1,&
                              & bvKIndex(1), bvKIndex(2), bvKIndex(3), iS)

                          tmpGradients(:, iAtK) = tmpGradients(:, iAtK)&
                              & + real(dPkmPhasedGammaAKSam * cmplx(Sbk * dPab, 0, dp), dp)
                        end do
                      end do
                    end do
                  end do
                end do

              end do loopKptsAK2
            end do loopGAK2

          end do loopA4
        end do loopB4
      end do loopM3
    end do loopK4

    ! term #4 with prefactor 1/8 of eq.(48)
    loopM4: do iAtM = 1, nAtom0
      descM = getDescriptor(iAtM, iSquare)
      loopN3: do iAtN = 1, nAtom0
        descN = getDescriptor(iAtN, iSquare)
        loopB5: do iNeighN = 0, nNeighbourCamSym(iAtN)
          iAtB = symNeighbourList%neighbourList%iNeighbour(iNeighN, iAtN)
          iAtBfold = symNeighbourList%img2CentCell(iAtB)
          iAtK = iAtBfold
          descK = getDescriptor(iAtK, iSquare)
          descB = getDescriptor(iAtBfold, iSquare)
          vecL(:) = cellVecs(:, symNeighbourList%iCellVec(iAtB))
          ! get 2D pointer to Sbn overlap block
          ind = symNeighbourList%iPair(iNeighN, iAtN) + 1
          nOrbAt = descN(iNOrb)
          nOrbNeigh = descB(iNOrb)
          pSbn(1:nOrbNeigh, 1:nOrbAt) => this%overSym(ind:ind + nOrbNeigh * nOrbAt - 1)
          loopA5: do iNeighM = 0, nNeighbourCamSym(iAtM)
            iAtA = symNeighbourList%neighbourList%iNeighbour(iNeighM, iAtM)
            iAtAfold = symNeighbourList%img2CentCell(iAtA)
            if (iAtK == iAtAfold) cycle
            descA = getDescriptor(iAtAfold, iSquare)
            vecH(:) = cellVecs(:, symNeighbourList%iCellVec(iAtA))
            dGammaAK = -this%camdGammaEvalG(iAtAfold, iAtK)%data
            ! get 2D pointer to Sam overlap block
            ind = symNeighbourList%iPair(iNeighM, iAtM) + 1
            nOrbAt = descM(iNOrb)
            nOrbNeigh = descA(iNOrb)
            pSam(1:nOrbNeigh, 1:nOrbAt) => this%overSym(ind:ind + nOrbNeigh * nOrbAt - 1)

            loopGAK3: do iG = 1, size(this%nNonZeroGammaG(iAtAfold, iAtK)%data)
              iGAK = this%nNonZeroGammaG(iAtAfold, iAtK)%data(iG)
              bvKIndex(:) = this%foldToBvKIndex(-this%cellVecsG(:, iGAK))

              loopKptsAK3: do iK = 1, nK
                phase = exp(imag * dot_product(2.0_dp * pi * kPoints(:, iK),&
                    & vecL - this%cellVecsG(:, iGAK) - vecH))
                phasedGammaAK(:) = phase * cmplx(dGammaAK(:, iG) * kWeights(iK), 0, dp)

                do mu = 1, descM(iNOrb)
                  do nu = 1, descN(iNOrb)
                    do iS = 1, nS
                      dPnm = deltaRhoOutSqrCplx(descN(iStart) + nu - 1, descM(iStart) + mu - 1,&
                          & iKiSToiGlobalKS(iK, iS))
                      dPnmPhasedGammaAK(:) = dPnm * phasedGammaAK
                      do alpha = 1, descA(iNOrb)
                        Sam = pSam(alpha, mu)
                        dPnmPhasedGammaAKSam(:) = dPnmPhasedGammaAK * cmplx(Sam, 0, dp)
                        do beta = 1, descB(iNOrb)
                          Sbn = pSbn(beta, nu)
                          dPab = deltaRhoSqr(descA(iStart) + alpha - 1, descB(iStart) + beta - 1,&
                              & bvKIndex(1), bvKIndex(2), bvKIndex(3), iS)

                          tmpGradients(:, iAtK) = tmpGradients(:, iAtK)&
                              & + real(dPnmPhasedGammaAKSam * cmplx(Sbn * dPab, 0, dp), dp)
                        end do
                      end do
                    end do
                  end do
                end do

              end do loopKptsAK3
            end do loopGAK3

          end do loopA5
        end do loopB5
      end do loopN3
    end do loopM4

    if (this%tREKS) then
      tmpGradients(:,:) = -0.25_dp * tmpGradients
    else
      tmpGradients(:,:) = -0.125_dp * nS * tmpGradients
    end if

    gradients(:,:) = gradients + tmpGradients

  end subroutine addCamGradientsNeighbour_kpts_ct


  !> Returns the array of atomic species in central cell.
  subroutine getCentralCellSpecies(this, species0)

    !> Class instance
    class(THybridXcFunc), intent(in) :: this

    !> 1D array for output, will be allocated
    integer, intent(out), allocatable :: species0(:)

    species0 = this%species0

  end subroutine getCentralCellSpecies


  !> Returns tabulated (long-range + full-range Hartree-Fock) gamma integrals.
  !! (non-periodic systems only)
  subroutine getCamGammaCluster(this, camGamma0)

    !> Class instance
    class(THybridXcFunc), intent(in) :: this

    !> Long-range + full-range Hartree-Fock gamma integrals in AO basis
    real(dp), intent(out) :: camGamma0(:,:)

    camGamma0(:,:) = this%camGammaEval0

  end subroutine getCamGammaCluster


  !> Calculates (long-range + full-range Hartree-Fock) gamma derivative integrals.
  !! (non-periodic systems only)
  subroutine getCamGammaDerivCluster(this, camGammaDeriv0)

    !> Class instance
    class(THybridXcFunc), intent(in) :: this

    !> Long-range + full-range Hartree-Fock gamma derivative integrals
    real(dp), intent(out) :: camGammaDeriv0(:,:,:)

    !! Holds long-range + full-range Hartree-Fock gamma derivatives of a single interaction
    real(dp) :: tmp(3)

    !! Number of atoms and indices of interacting atoms
    integer :: nAtom0, iAt1, iAt2

    nAtom0 = size(camGammaDeriv0, dim=1)

    do iAt1 = 1, nAtom0
      do iAt2 = 1, nAtom0
        if (iAt1 /= iAt2) then
          call getDirectionalCamGammaPrimeValue_cluster(this, tmp, iAt1, iAt2)
          camGammaDeriv0(iAt2, iAt1, :) = tmp
        end if
      end do
    end do

  end subroutine getCamGammaDerivCluster


  !> Calculates analytical full-range Hartree-Fock gamma.
  function getHfAnalyticalGammaValue(this, iSp1, iSp2, dist) result(gamma)

    !> Class instance
    class(THybridXcFunc_full), intent(in) :: this

    !> First species
    integer, intent(in) :: iSp1

    !> Second species
    integer, intent(in) :: iSp2

    !> Distance between atoms
    real(dp), intent(in) :: dist

    !> Resulting gamma
    real(dp) :: gamma

    gamma = getHfAnalyticalGammaValue_workhorse(this%hubbu(iSp1), this%hubbu(iSp2), dist)

  end function getHfAnalyticalGammaValue

end module dftbp_dftb_hybridxc
