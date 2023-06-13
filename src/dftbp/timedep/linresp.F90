!--------------------------------------------------------------------------------------------------!
!  DFTB+: general package for performing fast atomistic simulations                                !
!  Copyright (C) 2006 - 2023  DFTB+ developers group                                               !
!                                                                                                  !
!  See the LICENSE file for terms of usage and distribution.                                       !
!--------------------------------------------------------------------------------------------------!

#:include 'common.fypp'

!> Linear response formulation of TD-DFTB as developed by Niehaus et al.
!>
!> The functionality of the module has some limitation:
!> * Third order does not work.
!> * Periodic system do not work yet apart from Gamma point.
!> * Orbital potentials or spin-orbit does not work yet.
!> * Only for closed shell or colinear spin polarization (excitation energies only in that
!>   case).
!> * Onsite corrections are not included in this version
module dftbp_timedep_linresp
  use dftbp_common_accuracy, only : dp
  use dftbp_common_file, only : TFileDescr
  use dftbp_dftb_nonscc, only : TNonSccDiff
  use dftbp_dftb_scc, only : TScc
  use dftbp_dftb_slakocont, only : TSlakoCont
  use dftbp_extlibs_arpack, only : withArpack
  use dftbp_io_message, only : error
  use dftbp_io_taggedoutput, only : TTaggedWriter
  use dftbp_timedep_linrespgrad, only : LinRespGrad_old
  use dftbp_timedep_linresptypes, only : TLinResp, linrespSolverTypes
  use dftbp_type_commontypes, only : TOrbitals
  use dftbp_type_densedescr, only : TDenseDescr
  use dftbp_dftb_rangeseparated, only : TRangeSepFunc
  use dftbp_common_environment, only : TEnvironment
  implicit none

  private
  public :: TLinresp, TLinrespini
  public :: LinResp_init, linResp_calcExcitations, LinResp_addGradients

  !> Data type for initial values for linear response calculations
  type :: TLinrespini

    !> number of excitations to be found
    integer :: nExc

    !> Is an energy window being used?
    logical :: tEnergyWindow

    !> Is an oscillator window being used?
    logical :: tOscillatorWindow

    !> transitions to include that have energies above the single particle transition nExc
    real(dp) :: energyWindow

    !> single particle transitions to be included above the energy window if they are brighter than
    !> this state of interest
    real(dp) :: oscillatorWindow

    !> should transition charges be cached
    logical :: tCacheCharges

    !> number of excited states to find
    integer :: nStat

    !> symmetry of states being calculated
    character :: sym

    !> atom resolved Hubbard U
    real(dp), allocatable :: HubbardU(:)

    !> atom resolved spin constants
    real(dp), allocatable :: spinW(:)

    !> print excited state mulliken populations
    logical :: tMulliken

    !> Write expansion coefficients of excited states
    logical :: tCoeffs

    !> Include ground state in the excited state natural orbitals
    logical :: tGrndState

    !> Write natural orbitals for excited state
    logical :: tPrintEigVecs

    !> Should the density matrix be stored to disc?
    logical :: tWriteDensityMatrix

    !> write X+Y vector sqrt(wij) / sqrt(omega) * F^ia_I
    logical :: tXplusY

    !> Initial and final state for non-adiabatic coupling evaluation
    integer :: indNACouplings(2)

    !> write single particle transitions
    logical :: tSPTrans

    !> write more detail on excited state transitions
    logical :: tTrans

    !> write excited state transition charges
    logical :: tTransQ

    !> dipole strengths to excited states
    logical :: tTradip

    !> RPA solver choice
    integer :: iLinRespSolver

    !> subspace dimension factor Stratmann diagonaliser
    integer :: subSpaceFactorStratmann

    !> print state of Arnoldi solver
    logical :: tArnoldi

    !> diagnose output of Arnoldi solver
    logical :: tDiagnoseArnoldi

  end type TLinrespini


contains


  !> Initialize an internal data type for linear response excitations
  subroutine LinResp_init(this, ini, nAtom, nEl, nSpin, onSiteMatrixElements)

    !> data structure for linear response
    type(TLinResp), intent(out) :: this

    !> initial values for setting parameters
    type(TLinrespini), intent(inout) :: ini

    !> number of atoms in central cell
    integer, intent(in) :: nAtom

    !> number of electrons in total
    real(dp), intent(in) :: nEl

    !> Number of spin channels
    integer, intent(in) :: nSpin

    !> onsite corrections if in use
    real(dp), allocatable :: onSiteMatrixElements(:,:,:,:)

    integer :: dLev

    this%tinit = .false.

    this%iLinRespSolver = ini%iLinRespSolver

    if (any([linrespSolverTypes%Arpack, linrespSolverTypes%Stratmann] ==&
        & this%iLinRespSolver)) then
      this%tinit = .true.
    else
      call error('Internal error: Illegal routine call to LinResp_init.')
    end if

    if (this%iLinRespSolver == linrespSolverTypes%Arpack .and. .not. withArpack) then
      call error('This binary is buit without ARPACK support, but it is requested.')
    end if

    this%nExc = ini%nExc
    this%tEnergyWindow = ini%tEnergyWindow
    this%energyWindow = ini%energyWindow
    this%tOscillatorWindow = ini%tOscillatorWindow
    this%oscillatorWindow = ini%oscillatorWindow
    ! Final decision on value of tCacheChargesSame is made in linRespGrad
    this%tCacheChargesOccVir = ini%tCacheCharges
    this%tCacheChargesSame = ini%tCacheCharges
    this%nStat = ini%nStat
    this%symmetry = ini%sym

    this%tWriteDensityMatrix = ini%tWriteDensityMatrix

    if (nSpin == 1) then
      this%tSpin = .false.
    else if (nSpin == 2) then
      this%tSpin = .true.
    else
      call error("Unknown number of spin channels for excited state")
    end if

    if (this%tOscillatorWindow .and. this%OscillatorWindow <= 0.0_dp) then
      call error("Excited Oscillator window should be non-zero if used")
    end if
    if (this%tEnergyWindow .and. this%energyWindow <= 0.0_dp) then
      call error("Excited energy window should be non-zero if used")
    end if

    if(all(ini%indNACouplings == 0)) then
      this%tNaCoupling = .false.
    else
      if (any(ini%indNACouplings < 0)) then
        call error("StateCouplings: Indices must be positive.")
      end if
      if (ini%indNACouplings(1) >=  ini%indNACouplings(2)) then
        call error("StateCouplings: Second index must be larger than first one.")
      end if
      if (ini%nExc < ini%indNACouplings(2)) then
        call error('StateCouplings: Index must not exceed number of states to calculate.')
      end if
      if (this%tSpin) then
        call error('StateCouplings: Spin-polarized systems currently not available.')
      end if
      this%tNaCoupling = .true.
      this%indNACouplings = ini%indNACouplings
      dLev = ini%indNACouplings(2) - ini%indNACouplings(1)
    endif

    this%writeMulliken = ini%tMulliken
    this%writeCoeffs = ini%tCoeffs
    this%tGrndState = ini%tGrndState
    this%writeTrans = ini%tTrans
    this%writeTransQ = ini%tTransQ
    this%writeSPTrans = ini%tSPTrans
    this%writeXplusY = ini%tXplusY
    this%writeTransDip = ini%tTradip

    this%nAtom = nAtom
    this%nEl = nEl

    call move_alloc(ini%spinW, this%spinW)
    call move_alloc(ini%hubbardU, this%HubbardU)
    if (allocated(onSiteMatrixElements)) then
      allocate(this%onSiteMatrixElements(size(onSiteMatrixElements,dim=1),&
          & size(onSiteMatrixElements,dim=2), size(onSiteMatrixElements,dim=3),&
          & size(onSiteMatrixElements,dim=4)))
      this%onSiteMatrixElements(:,:,:,:) = onSiteMatrixElements
    end if

    select case(this%iLinRespSolver)
    case(linrespSolverTypes%Arpack)
      this%testArnoldi = ini%tDiagnoseArnoldi
      this%tArnoldi = ini%tArnoldi
    case (linrespSolverTypes%Stratmann)
      this%subSpaceFactorStratmann = ini%subSpaceFactorStratmann
    end select


  end subroutine LinResp_init


  !> Wrapper to call the actual linear response routine for excitation energies
  subroutine linResp_calcExcitations(env, this, tSpin, denseDesc, eigVec, eigVal, SSqrReal, filling,&
      & coords0, sccCalc, dqAt, species0, iNeighbour, img2CentCell, orb, tWriteTagged, fdTagged,&
      & taggedWriter, rangeSep, excEnergy, allExcEnergies)

    !> Environment settings
    type(TEnvironment), intent(in) :: env
    
    !> data structure with additional linear response values
    type(TLinresp), intent(inout) :: this

    !> is this a spin-polarized calculation
    logical, intent(in) :: tSpin

    !> Indexing array for dense H and S
    type(TDenseDescr), intent(in) :: denseDesc

    !> ground state eigenvectors
    real(dp), intent(in) :: eigVec(:,:,:)

    !> ground state eigenvalues
    real(dp), intent(in) :: eigVal(:,:)

    !> square overlap matrix
    real(dp), intent(in) :: SSqrReal(:,:)

    !> ground state occupations
    real(dp), intent(in) :: filling(:,:)

    !> central cell atomic coordinates
    real(dp), intent(in) :: coords0(:,:)

    !> Self-consistent charge module settings
    type(TScc), intent(in) :: sccCalc

    !> Gross Mulliken atomic charges for ground state
    real(dp), intent(in) :: dqAt(:,:)

    !> chemical type of atoms in central cell
    integer, intent(in) :: species0(:)

    !> index array for atomic neighbours
    integer, intent(in) :: img2CentCell(:)

    !> folding back to the central cell
    integer, intent(in) :: iNeighbour(0:,:)

    !> data type with atomic orbital information
    type(TOrbitals), intent(in) :: orb

    !> print tag information
    logical, intent(in) :: tWriteTagged

    !> file id for tagging information
    type(TFileDescr), intent(in) :: fdTagged

    !> tagged writer
    type(TTaggedWriter), intent(inout) :: taggedWriter

    !> Data for range separated calcualtion
    type(TRangeSepFunc), allocatable, intent(inout) :: rangeSep

    !> excitation energy (only when nStat /=0, othewise set numerically 0)
    real(dp), intent(out) :: excEnergy

    !> energies of all solved states
    real(dp), intent(inout), allocatable :: allExcEnergies(:)

    real(dp), pointer :: dummyPtr(:,:,:) => null()

    if (this%tInit) then
      @:ASSERT(size(orb%nOrbAtom) == this%nAtom)
      call LinRespGrad_old(env, this, denseDesc, eigVec, eigVal, sccCalc, dqAt, coords0,&
          & SSqrReal, filling, species0, iNeighbour, img2CentCell, orb, fdTagged, taggedWriter,&
          & rangeSep, excEnergy, allExcEnergies, dummyPtr)
    else
      call error('Internal error: Illegal routine call to LinResp_calcExcitations.')
    end if

  end subroutine linResp_calcExcitations


  !> Wrapper to call linear response calculations of excitations and forces in excited states
  subroutine LinResp_addGradients(env, tSpin, this, denseDesc, eigVec, eigVal, SSqrReal, filling,&
      & coords0, sccCalc, dqAt, species0, iNeighbour, img2CentCell, orb, skHamCont, skOverCont,&
      & fdTagged, taggedWriter, rangeSep, excEnergy, allExcEnergies, excgradient,&
      & derivator, rhoSqr, deltaRho, occNatural, naturalOrbs)

    !> Environment settings
    type(TEnvironment), intent(in) :: env
    
    !> is this a spin-polarized calculation
    logical, intent(in) :: tSpin

    !> data for the actual calculation
    type(TLinResp), intent(inout) :: this

    !> Indexing array for dense H and S
    type(TDenseDescr), intent(in) :: denseDesc

    !> ground state eigenvectors
    real(dp), intent(in) :: eigVec(:,:,:)

    !> ground state eigenvalues
    real(dp), intent(in) :: eigVal(:,:)

    !> square overlap matrix (must be symmetrised)
    real(dp), intent(in) :: SSqrReal(:,:)

    !> ground state occupations
    real(dp), intent(in) :: filling(:,:)

    !> central cell atomic coordinates
    real(dp), intent(in) :: coords0(:,:)

    !> Self-consistent charge module settings
    type(TScc), intent(in) :: sccCalc

    !> Gross atomic charges in ground state
    real(dp), intent(in) :: dqAt(:,:)

    !> chemical species of atoms in central cell
    integer, intent(in) :: species0(:)

    !> index array for atoms within cutoff distances
    integer, intent(in) :: iNeighbour(0:,:)

    !> folding to central cell (not really needed for non-periodic systems)
    integer, intent(in) :: img2CentCell(:)

    !> orbital data structure
    type(TOrbitals), intent(in) :: orb

    !> non-SCC H0 data
    type(TSlakoCont), intent(in) :: skHamCont

    !> overlap data
    type(TSlakoCont), intent(in) :: skOverCont

    !> method for calculating derivatives of S and H0 matrices
    class(TNonSccDiff), intent(in) :: derivator

    !> ground state density matrix (square matrix plus spin index)
    real(dp), intent(in) :: rhoSqr(:,:,:)

    !> difference density matrix (vs. uncharged atoms)
    real(dp), intent(inout), pointer :: deltaRho(:,:,:)

    !> file descriptor for tagged data
    type(TFileDescr), intent(in) :: fdTagged

    !> Tagged writer
    type(TTaggedWriter), intent(inout) :: taggedWriter

    !> Data for range-separated calculation
    type(TRangeSepFunc), allocatable, intent(inout) :: rangeSep

    !> energy of particular excited state
    real(dp), intent(out) :: excenergy

    !> energies of all solved states
    real(dp), intent(inout), allocatable :: allExcEnergies(:)

    !> contribution to forces from derivative of excited state energy
    real(dp), intent(inout), allocatable :: excgradient(:,:)

    !> occupations of the natural orbitals from the density matrix
    real(dp), intent(inout), allocatable :: occNatural(:)

    !> the natural orbitals of the excited state transition density matrix or the total density
    !> matrix in the excited state
    real(dp), intent(inout), allocatable :: naturalOrbs(:,:,:)

    real(dp), allocatable :: shiftPerAtom(:), shiftPerL(:,:)

    if (this%tInit) then
      @:ASSERT(this%nAtom == size(orb%nOrbAtom))
      @:ASSERT(allocated(occNatural) .eqv. allocated(naturalOrbs))

      allocate(shiftPerAtom(this%nAtom))
      allocate(shiftPerL(orb%mShell, this%nAtom))
      call sccCalc%getShiftPerAtom(shiftPerAtom)
      call sccCalc%getShiftPerL(shiftPerL)
      shiftPerAtom = shiftPerAtom + shiftPerL(1,:)

      if (allocated(occNatural)) then
        call LinRespGrad_old(env, this, denseDesc, eigVec, eigVal, sccCalc, dqAt, coords0, SSqrReal,&
            & filling, species0, iNeighbour, img2CentCell, orb, fdTagged, taggedWriter, rangeSep,&
            & excEnergy, allExcEnergies, deltaRho, shiftPerAtom, skHamCont, skOverCont,&
            & excgradient, derivator, rhoSqr, occNatural, naturalOrbs)
      else
        call LinRespGrad_old(env, this, denseDesc, eigVec, eigVal, sccCalc, dqAt, coords0, SSqrReal,&
            & filling, species0, iNeighbour, img2CentCell, orb, fdTagged, taggedWriter, rangeSep,&
            & excEnergy, allExcEnergies, deltaRho, shiftPerAtom, skHamCont, skOverCont,&
            & excgradient, derivator, rhoSqr)
      end if

    else
      call error('Internal error: Illegal routine call to LinResp_addGradients.')
    endif

  end subroutine LinResp_addGradients

end module dftbp_timedep_linresp
