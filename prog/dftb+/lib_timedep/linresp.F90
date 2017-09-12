!--------------------------------------------------------------------------------------------------!
!  DFTB+: general package for performing fast atomistic simulations                                !
!  Copyright (C) 2017  DFTB+ developers group                                                      !
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
module linresp_module
  use assert
  use accuracy
  use message
  use commontypes
  use slakocont
  use fileid
  use scc, only : getShiftPerAtom, getShiftPerL, typSCC
  use nonscc, only : NonSccDiff
#:if WITH_ARPACK
  ! code is compiled with arpack available
  use linrespgrad
#:endif
  implicit none
  private

  public :: linresp, linrespini
  public :: init, calcExcitations, addGradients


  !> Data type for initial values for linear response calculations
  type :: linrespini

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

    !> write X+Y vector sqrt(wij) / sqrt(omega) * F^ia_I
    logical :: tXplusY

    !> write single particle transitions
    logical :: tSPTrans

    !> write more detail on excited state transitions
    logical :: tTrans

    !> dipole strengths to excited states
    logical :: tTradip

    !> print state and diagnose output of Arnoldi solver
    logical :: tArnoldi, tDiagnoseArnoldi

    !> Initialised data structure?
    logical :: tInit = .false.
  end type linrespini


  !> Data type for linear response internal settings
  type :: linresp
    integer :: nExc, nStat
    logical :: tEnergyWindow
    real(dp) :: energyWindow
    logical :: tOscillatorWindow
    real(dp) :: oscillatorWindow
    integer :: nOcc, nVir, nAtom
    real(dp) :: nEl
    character :: symmetry
    real(dp), allocatable :: spinW(:)
    real(dp), allocatable :: HubbardU(:)
    integer :: fdXplusY = -1
    integer :: fdCoeffs = -1
    logical :: tGrndState = .true.
    integer :: fdMulliken = -1
    integer :: fdTrans = -1
    integer :: fdSPTrans = -1
    integer :: fdExc = -1
    integer :: fdTradip = -1
    logical :: tArnoldi


    !> file unit for Arnoldi solver file unit for tests on output of Arnoldi solver
    integer :: fdArnoldi = -1

    integer :: fdArnoldiDiagnosis = -1
    logical :: tPrintEigVecs
    logical :: tInit = .false.
  end type linresp


  !> Initialise data structure
  interface init
    module procedure LinResp_init
  end interface init


  !> actually calculate excitations
  interface calcExcitations
    module procedure LinResp_calcExcitations
  end interface calcExcitations


  !> excitations plus forces and some other properties (excited
  interface addGradients
    module procedure LinResp_addGradients
  end interface addGradients

contains


  !> Initialize an internal data type for linear response excitations
  subroutine LinResp_init(self, ini, nAtom, nEl, orb)

    !> data structure for linear response
    type(linresp), intent(out) :: self

    !> initial values for setting parameters
    type(linrespini), intent(inout) :: ini

    !> number of atoms in central cell
    integer, intent(in) :: nAtom

    !> number of electrons in total
    real(dp), intent(in) :: nEl

    !> data on atomic orbitals
    type(TOrbitals), intent(in) :: orb

#:if WITH_ARPACK
    self%nExc = ini%nExc
    self%tEnergyWindow = ini%tEnergyWindow
    self%energyWindow = ini%energyWindow
    self%tOscillatorWindow = ini%tOscillatorWindow
    self%oscillatorWindow = ini%oscillatorWindow
    self%nStat = ini%nStat
    self%symmetry = ini%sym

    if (self%tOscillatorWindow .and. self%OscillatorWindow <= 0.0_dp) then
      call error("Excited Oscillator window should be non-zero if used")
    end if
    if (self%tEnergyWindow .and. self%energyWindow <= 0.0_dp) then
      call error("Excited energy window should be non-zero if used")
    end if

    if (ini%tMulliken) then
      self%fdMulliken = getFileId()
    else
      self%fdMulliken = -1
    end if
    if (ini%tCoeffs) then
      self%fdCoeffs = getFileId()
    else
      self%fdCoeffs = -1
    end if
    self%tGrndState = ini%tGrndState
    if (ini%tDiagnoseArnoldi) then
      self%fdArnoldiDiagnosis = getFileId()
    else
      self%fdArnoldiDiagnosis = -1
    end if
    if (ini%tTrans) then
      self%fdTrans = getFileId()
    else
      self%fdTrans = -1
    end if
    if (ini%tSPTrans) then
      self%fdSPTrans = getFileId()
    else
      self%fdSPTrans = -1
    end if
    if (ini%tXplusY) then
      self%fdXplusY = getFileId()
    else
      self%fdXplusY = -1
    end if
    if (ini%tTradip) then
      self%fdTradip = getFileId()
    else
      self%fdTradip = -1
    end if

    self%tArnoldi = ini%tArnoldi
    self%fdArnoldi = getFileId()
    self%nAtom = nAtom
    self%nEl = nEl
    self%nOcc = ceiling(nEl / 2.0_dp)
    self%nVir = orb%nOrb - self%nOcc
    self%fdExc = getFileId() ! file for excitations

    ! Write to disc
    self%tPrintEigVecs = ini%tPrintEigVecs

    call move_alloc(ini%spinW, self%spinW)
    call move_alloc(ini%hubbardU, self%HubbardU)
    self%tinit = .true.
#:else
    self%tinit = .false.
    call error('Internal error: Illegal routine call to LinResp_init.')
#:endif

  end subroutine LinResp_init


  !> Wrapper to call the actual linear response routine for excitation energies
  subroutine LinResp_calcExcitations(tSpin, self, iAtomStart, eigVec, eigVal, SSqrReal, filling, &
      & coords0, oSCC, dqAt, species0, iNeighbor, img2CentCell, orb, tWriteTagged, fdTagged, &
      & excEnergy)

    !> is this a spin-polarized calculation
    logical, intent(in) :: tSpin

    !> data structure with additional linear response values
    type(linresp), intent(inout) :: self

    !> indexing array for ground state square matrices
    integer, intent(in) :: iAtomStart(:)

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
    type(typSCC), intent(in) :: oSCC

    !> net Mulliken atomic charges for ground state
    real(dp), intent(in) :: dqAt(:)

    !> chemical type of atoms in central cell
    integer, intent(in) :: species0(:)

    !> index array for atomic neighbours
    integer, intent(in) :: img2CentCell(:)

    !> folding back to the central cell
    integer, intent(in) :: iNeighbor(0:,:)

    !> data type with atomic orbital information
    type(TOrbitals), intent(in) :: orb

    !> print tag information
    logical, intent(in) :: tWriteTagged

    !> file id for tagging information
    integer, intent(in) :: fdTagged

    !> excitation energy (only when nStat /=0, othewise set numerically 0)
    real(dp), intent(out) :: excEnergy

#:if WITH_ARPACK
    @:ASSERT(self%tInit)
    @:ASSERT(size(orb%nOrbAtom) == self%nAtom)
    call LinRespGrad_old(tSpin, self%nAtom, iAtomStart, eigVec, eigVal, oSCC, dqAt, coords0, &
        & self%nExc, self%nStat, self%symmetry, SSqrReal, filling, species0, self%HubbardU, &
        & self%spinW, self%nEl, iNeighbor, img2CentCell, orb, tWriteTagged, fdTagged, &
        & self%fdMulliken, self%fdCoeffs, self%tGrndState, self%fdXplusY, self%fdTrans, &
        & self%fdSPTrans, self%fdTradip, self%tArnoldi, self%fdArnoldi,  self%fdArnoldiDiagnosis, &
        & self%fdExc, self%tEnergyWindow, self%energyWindow, self%tOscillatorWindow, &
        & self%oscillatorWindow, excEnergy)

#:else
    call error('Internal error: Illegal routine call to &
        &LinResp_calcExcitations')
#:endif

  end subroutine LinResp_calcExcitations


  !> Wrapper to call linear response calculations of excitations and forces in excited states
  subroutine LinResp_addGradients(tSpin, self, iAtomStart, eigVec, eigVal, SSqrReal, filling, &
      & coords0, oSCC, dqAt, species0, iNeighbor, img2CentCell, orb, skHamCont, skOverCont, &
      & tWriteTagged, fdTagged, excEnergy, excgradient, derivator, rhoSqr, occNatural, &
      & naturalOrbs)

    !> is this a spin-polarized calculation
    logical, intent(in) :: tSpin

    !> data for the actual calculation
    type(linresp), intent(inout) :: self

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
    type(typSCC), intent(in) :: oSCC

    !> net atomic charges in ground state
    real(dp), intent(in) :: dqAt(:)

    !> chemical species of atoms in central cell
    integer, intent(in) :: species0(:)

    !> index array for atoms within cutoff distances
    integer, intent(in) :: iNeighbor(0:,:)

    !> folding to central cell (not really needed for non-periodic systems)
    integer, intent(in) :: img2CentCell(:)

    !> orbital data structure
    type(TOrbitals), intent(in) :: orb

    !> non-SCC H0 data
    type(OSlakoCont), intent(in) :: skHamCont

    !> overlap data
    type(OSlakoCont), intent(in) :: skOverCont

    !> print tag information
    logical, intent(in) :: tWriteTagged

    !> file descriptor for tagged data
    integer, intent(in) :: fdTagged

    !> energy of particular excited state
    real(dp), intent(out) :: excenergy

    !> contribution to forces from derivative of excited state energy
    real(dp), intent(out) :: excgradient(:,:)

    !> method for calculating derivatives of S and H0 matrices
    class(NonSccDiff), intent(in), optional :: derivator

    !> ground state density matrix (square matrix plus spin index)
    real(dp), intent(in), optional :: rhoSqr(:,:,:)

    !> occupations of the natural orbitals from the density matrix
    real(dp), intent(out), optional :: occNatural(:)

    !> the natural orbitals of the excited state transition density matrix or the total density
    !> matrix in the excited state
    real(dp), intent(out), optional :: naturalOrbs(:,:,:)

#:if WITH_ARPACK

    real(dp), allocatable :: shiftPerAtom(:), shiftPerL(:,:)
    @:ASSERT(self%tInit)
    @:ASSERT(self%nAtom == size(orb%nOrbAtom))
    ! BA: SCC is currently ugly, it gives back an array with an additional dimension (spin),
    ! however, fills always the first channel only!
    ALLOCATE(shiftPerAtom(self%nAtom))
    ALLOCATE(shiftPerL(orb%mShell, self%nAtom))
    call getShiftPerAtom(shiftPerAtom, oSCC)
    call getShiftPerL(shiftPerL, oSCC)
    shiftPerAtom = shiftPerAtom + shiftPerL(1,:)

    call LinRespGrad_old(tSpin, self%nAtom, iAtomStart, eigVec, eigVal, oSCC, dqAt, coords0, &
        & self%nExc, self%nStat, self%symmetry, SSqrReal, filling, species0, self%HubbardU, &
        & self%spinW, self%nEl, iNeighbor, img2CentCell, orb, tWriteTagged, fdTagged, &
        & self%fdMulliken, self%fdCoeffs, self%tGrndState, self%fdXplusY, self%fdTrans, &
        & self%fdSPTrans, self%fdTradip, self%tArnoldi, self%fdArnoldi, self%fdArnoldiDiagnosis, &
        & self%fdExc, self%tEnergyWindow, self%energyWindow, self%tOscillatorWindow, &
        & self%oscillatorWindow, excEnergy, shiftPerAtom, skHamCont, skOverCont, excgradient, &
        & derivator, rhoSqr, occNatural, naturalOrbs)

#:else
    call error('Internal error: Illegal routine call to LinResp_addGradients.')
#:endif

  end subroutine LinResp_addGradients

end module linresp_module
