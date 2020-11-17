!--------------------------------------------------------------------------------------------------!
!  DFTB+: general package for performing fast atomistic simulations                                !
!  Copyright (C) 2006 - 2020  DFTB+ developers group                                               !
!                                                                                                  !
!  See the LICENSE file for terms of usage and distribution.                                       !
!--------------------------------------------------------------------------------------------------!

#:include 'common.fypp'

!> Linear response formulation of TD-DFTB as developed by Niehaus et al.
!>
!> The functionality of the module has some limitation:
!> * Third order does not work.
!> * Periodic system do not work yet appart from Gamma point.
!> * Orbital potentials or spin-orbit does not work yet.
!> * Only for closed shell or colinear spin polarization (excitation energies only in that
!>   case).
!> * Onsite corrections are not included in this version
module dftbp_linresp
  use dftbp_assert
  use dftbp_accuracy
  use dftbp_message
  use dftbp_commontypes
  use dftbp_slakocont
  use dftbp_fileid
  use dftbp_scc, only : TScc
  use dftbp_nonscc, only : TNonSccDiff
  use dftbp_densedescr
  use dftbp_taggedoutput, only : TTaggedWriter
  use dftbp_linrespgrad
  use dftbp_arpack, only : withArpack
  use dftbp_linresptypes
  implicit none
  private

  public :: TLinresp, TLinrespini
  public :: LinResp_init, linResp_calcExcitations, addGradients

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

    !> write single particle transitions
    logical :: tSPTrans

    !> write more detail on excited state transitions
    logical :: tTrans

    !> write excited state transition charges
    logical :: tTransQ

    !> dipole strengths to excited states
    logical :: tTradip

    !> print state of Arnoldi solver
    logical :: tArnoldi

    !> diagnose output of Arnoldi solver
    logical :: tDiagnoseArnoldi

    !> Initialised data structure?
    logical :: tInit = .false.

  end type TLinrespini


  !> excitations plus forces and some other properties (excited
  interface addGradients
    module procedure LinResp_addGradients
  end interface addGradients

contains


  !> Initialize an internal data type for linear response excitations
  subroutine LinResp_init(this, ini, nAtom, nEl, onSiteMatrixElements)

    !> data structure for linear response
    type(TLinResp), intent(out) :: this

    !> initial values for setting parameters
    type(TLinrespini), intent(inout) :: ini

    !> number of atoms in central cell
    integer, intent(in) :: nAtom

    !> number of electrons in total
    real(dp), intent(in) :: nEl

    !> onsite corrections if in use
    real(dp), allocatable :: onSiteMatrixElements(:,:,:,:)

    this%tinit = .false.
    if (withArpack) then

      this%nExc = ini%nExc
      this%tEnergyWindow = ini%tEnergyWindow
      this%energyWindow = ini%energyWindow
      this%tOscillatorWindow = ini%tOscillatorWindow
      this%oscillatorWindow = ini%oscillatorWindow
      this%tCacheCharges = ini%tCacheCharges
      this%nStat = ini%nStat
      this%symmetry = ini%sym

      this%tWriteDensityMatrix = ini%tWriteDensityMatrix

      if (this%tOscillatorWindow .and. this%OscillatorWindow <= 0.0_dp) then
        call error("Excited Oscillator window should be non-zero if used")
      end if
      if (this%tEnergyWindow .and. this%energyWindow <= 0.0_dp) then
        call error("Excited energy window should be non-zero if used")
      end if

      if (ini%tMulliken) then
        this%fdMulliken = getFileId()
      else
        this%fdMulliken = -1
      end if
      if (ini%tCoeffs) then
        this%fdCoeffs = getFileId()
      else
        this%fdCoeffs = -1
      end if
      this%tGrndState = ini%tGrndState

      if (ini%tTrans) then
        this%fdTrans = getFileId()
      else
        this%fdTrans = -1
      end if

      if (ini%tTransQ) then
        this%fdTransQ = getFileId()
      else
        this%fdTransQ = -1
      end if

      if (ini%tSPTrans) then
        this%fdSPTrans = getFileId()
      else
        this%fdSPTrans = -1
      end if
      if (ini%tXplusY) then
        this%fdXplusY = getFileId()
      else
        this%fdXplusY = -1
      end if
      if (ini%tTradip) then
        this%fdTradip = getFileId()
      else
        this%fdTradip = -1
      end if

      this%nAtom = nAtom
      this%nEl = nEl
      this%fdExc = getFileId() ! file for excitations

      call move_alloc(ini%spinW, this%spinW)
      call move_alloc(ini%hubbardU, this%HubbardU)
      if (allocated(onSiteMatrixElements)) then
        allocate(this%onSiteMatrixElements(size(onSiteMatrixElements,dim=1),&
            & size(onSiteMatrixElements,dim=2), size(onSiteMatrixElements,dim=3),&
            & size(onSiteMatrixElements,dim=4)))
        this%onSiteMatrixElements(:,:,:,:) = onSiteMatrixElements
      end if

      this%tinit = .true.
    end if

    if (withArpack) then

      if (ini%tDiagnoseArnoldi) then
        this%fdArnoldiDiagnosis = getFileId()
      else
        this%fdArnoldiDiagnosis = -1
      end if
      this%tArnoldi = ini%tArnoldi
      this%fdArnoldi = getFileId()
      this%tinit = .true.

    else

      call error('Internal error: Illegal routine call to LinResp_init.')

    endif

  end subroutine LinResp_init


  !> Wrapper to call the actual linear response routine for excitation energies
  subroutine linResp_calcExcitations(this, tSpin, denseDesc, eigVec, eigVal, SSqrReal, filling,&
      & coords0, sccCalc, dqAt, species0, iNeighbour, img2CentCell, orb, tWriteTagged, fdTagged,&
      & taggedWriter, excEnergy, allExcEnergies)

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

    !> This-consistent charge module settings
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
    integer, intent(in) :: fdTagged

    !> tagged writer
    type(TTaggedWriter), intent(inout) :: taggedWriter

    !> excitation energy (only when nStat /=0, othewise set numerically 0)
    real(dp), intent(out) :: excEnergy

    !> energes of all solved states
    real(dp), intent(inout), allocatable :: allExcEnergies(:)

    if (withArpack) then
      @:ASSERT(this%tInit)
      @:ASSERT(size(orb%nOrbAtom) == this%nAtom)
      call LinRespGrad_old(tSpin, this, denseDesc%iAtomStart, eigVec, eigVal, sccCalc, dqAt,&
          & coords0, SSqrReal, filling, species0, iNeighbour, img2CentCell, orb, tWriteTagged,&
          & fdTagged, taggedWriter, excEnergy, allExcEnergies)
    else
      call error('Internal error: Illegal routine call to LinResp_calcExcitations')
    end if

  end subroutine linResp_calcExcitations


  !> Wrapper to call linear response calculations of excitations and forces in excited states
  subroutine LinResp_addGradients(tSpin, this, iAtomStart, eigVec, eigVal, SSqrReal, filling,&
      & coords0, sccCalc, dqAt, species0, iNeighbour, img2CentCell, orb, skHamCont, skOverCont,&
      & tWriteTagged, fdTagged, taggedWriter, excEnergy, allExcEnergies, excgradient, derivator,&
      & rhoSqr, occNatural, naturalOrbs)

    !> is this a spin-polarized calculation
    logical, intent(in) :: tSpin

    !> data for the actual calculation
    type(TLinResp), intent(inout) :: this

    !> indexing array for ground state square matrices
    integer, intent(in) :: iAtomStart(:)

    !> ground state eigenvectors
    real(dp), intent(in) :: eigVec(:,:,:)

    !> ground state eigenvalues
    real(dp), intent(in) :: eigVal(:,:)

    !> square overlap matrix (must be symmetriezed)
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

    !> print tag information
    logical, intent(in) :: tWriteTagged

    !> file descriptor for tagged data
    integer, intent(in) :: fdTagged

    !> Tagged writer
    type(TTaggedWriter), intent(inout) :: taggedWriter

    !> energy of particular excited state
    real(dp), intent(out) :: excenergy

    !> energes of all solved states
    real(dp), intent(inout), allocatable :: allExcEnergies(:)

    !> contribution to forces from derivative of excited state energy
    real(dp), intent(out) :: excgradient(:,:)

    !> occupations of the natural orbitals from the density matrix
    real(dp), intent(inout), allocatable :: occNatural(:)

    !> the natural orbitals of the excited state transition density matrix or the total density
    !> matrix in the excited state
    real(dp), intent(inout), allocatable :: naturalOrbs(:,:,:)

    real(dp), allocatable :: shiftPerAtom(:), shiftPerL(:,:)

    if (withArpack) then
      @:ASSERT(this%tInit)
      @:ASSERT(this%nAtom == size(orb%nOrbAtom))
      @:ASSERT(allocated(occNatural) .eqv. allocated(naturalOrbs))

      allocate(shiftPerAtom(this%nAtom))
      allocate(shiftPerL(orb%mShell, this%nAtom))
      call sccCalc%getShiftPerAtom(shiftPerAtom)
      call sccCalc%getShiftPerL(shiftPerL)
      shiftPerAtom = shiftPerAtom + shiftPerL(1,:)

      if (allocated(occNatural)) then
        call LinRespGrad_old(tSpin, this, iAtomStart, eigVec, eigVal, sccCalc, dqAt, coords0,&
            & SSqrReal, filling, species0, iNeighbour, img2CentCell, orb, tWriteTagged, fdTagged,&
            & taggedWriter, excEnergy, allExcEnergies, shiftPerAtom, skHamCont, skOverCont,&
            & excgradient, derivator, rhoSqr, occNatural, naturalOrbs)
      else
        call LinRespGrad_old(tSpin, this, iAtomStart, eigVec, eigVal, sccCalc, dqAt, coords0,&
            & SSqrReal, filling, species0, iNeighbour, img2CentCell, orb, tWriteTagged, fdTagged,&
            & taggedWriter, excEnergy, allExcEnergies, shiftPerAtom, skHamCont, skOverCont,&
            & excgradient, derivator, rhoSqr)
      end if

    else
      call error('Internal error: Illegal routine call to LinResp_addGradients.')
    endif

  end subroutine LinResp_addGradients

end module dftbp_linresp
