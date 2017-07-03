!--------------------------------------------------------------------------------------------------!
!  DFTB+: general package for performing fast atomistic simulations                                !
!  Copyright (C) 2017  DFTB+ developers group                                                      !
!                                                                                                  !
!  See the LICENSE file for terms of usage and distribution.                                       !
!--------------------------------------------------------------------------------------------------!

#:include 'common.fypp'

!> Linear response formulation of TD-DFTB as developed by Niehaus et al.
!!
!! \note The functionality of the module has some limitation:
!!   o Third order does not work.
!!   o Periodic system do not work yet appart from Gamma point.
!!   o Orbital potentials or spin-orbit does not work yet.
!!   o Only for closed shell or colinear spin polarization (excitation
!!     energies only in that case).
!!   o Onsite corrections are not included in this version
!!
module linresp_module
  use assert
  use accuracy
  use message
  use commontypes
  use slakocont
  use fileid
  use scc, only : getShiftPerAtom, getShiftPerL
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
    integer :: nExc ! number of excitations to be found
    logical :: tEnergyWindow, tOscillatorWindow ! are windows being used?
    real(dp) :: energyWindow ! energies to include above nExc
    real(dp) :: oscillatorWindow ! single particle transitions to be
                                 ! included above the energy window if
                                 ! they are brighter than this
    ! state of interest (warning - this is also used as a signaling variable)
    integer :: nStat

    ! symmetry of states being calculated
    character :: sym
    ! homo- and hetero-spin Hubbard U derivatives as W
    real(dp), allocatable :: HubbardU(:)
    real(dp), allocatable :: spinW(:)
    ! print excited state mulliken populations
    logical :: tMulliken
    ! Expansion coefficients of excited states
    logical :: tCoeffs
    ! Include ground state in the excited state natural orbitals
    logical :: tGrndState
    ! Write natural orbitals for excited state
    logical :: tPrintEigVecs
    ! write X+Y vector sqrt(wij) / sqrt(omega) * F^ia_I
    logical :: tXplusY
    ! single particle transitions
    logical :: tSPTrans
    ! transitions to excited states
    logical :: tTrans
    ! dipole strengths to excited states
    logical :: tTradip
    ! print state and diagnose output of Arnoldi solver
    logical :: tArnoldi, tDiagnoseArnoldi
    ! Initialised data structure?
    logical :: tInit = .false.
  end type linrespini

  !> Data type for linear response internal settings
  type :: linresp
    integer   :: nExc, nStat
    logical   :: tEnergyWindow
    real(dp)  :: energyWindow
    logical   :: tOscillatorWindow
    real(dp)  :: oscillatorWindow
    integer   :: nOcc, nVir, nAtom
    real(dp)  :: nEl
    character :: symmetry
    real(dp), allocatable :: spinW(:)
    real(dp), allocatable :: HubbardU(:)
    integer   :: fdXplusY = -1
    integer   :: fdCoeffs = -1
    logical   :: tGrndState = .true.
    integer   :: fdMulliken = -1
    integer   :: fdTrans = -1
    integer   :: fdSPTrans = -1
    integer   :: fdExc = -1
    integer   :: fdTradip = -1
    logical   :: tArnoldi
    ! file unit for Arnoldi solver file unit for tests on output of Arnoldi
    ! solver
    integer   :: fdArnoldi = -1
    integer   :: fdArnoldiDiagnosis = -1
    logical   :: tPrintEigVecs
    logical   :: tInit = .false.
  end type linresp

  interface init
    module procedure LinResp_init
  end interface

  interface calcExcitations
    module procedure LinResp_calcExcitations
  end interface

  interface addGradients
    module procedure LinResp_addGradients
  end interface

contains

  !> Initialize an internal data type for linear response excitations
  !! \param self data structure for linear response
  !! \param ini initial values for setting parameters
  !! \param nAtom number of atoms in central cell
  !! \param nEl number of electrons in total
  !! \param orb data on atomic orbitals
  subroutine LinResp_init(self, ini, nAtom, nEl, orb)
    type(linresp), intent(out) :: self
    type(linrespini), intent(inout) :: ini
    integer, intent(in) :: nAtom
    real(dp), intent(in) :: nEl
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
  !! \param tSpin is this a spin-polarized calculation
  !! \param self data structure with additional linear response values
  !! \param iAtomStart indexing array for ground state square matrices
  !! \param eigVec ground state eigenvectors
  !! \param eigVal ground state eigenvalues
  !! \param SSqrReal square overlap matrix
  !! \param filling ground state occupations
  !! \param coords0 central cell atomic coordinates
  !! \param dqAt net Mulliken atomic charges for ground state
  !! \param species0 chemical type of atoms in central cell
  !! \param iNeighbor index array for atomic neighbours
  !! \param img2centcell folding back to the central cell
  !! \param orb data type with atomic orbital information
  !! \param tWriteTagged print tag information
  !! \param fdTagged file id for tagging information
  !! \param excEnergy excitation energy (only when nStat /=0, othewise
  !! set numerically 0)
  subroutine LinResp_calcExcitations(tSpin, self, iAtomStart, eigVec, eigVal, SSqrReal, filling, &
      & coords0, dqAt, species0, iNeighbor, img2CentCell, orb, tWriteTagged, fdTagged, excEnergy)
    logical, intent(in) :: tSpin
    type(linresp), intent(inout) :: self
    integer, intent(in) :: iAtomStart(:)
    real(dp), intent(in) :: eigVec(:,:,:)
    real(dp), intent(in) :: eigVal(:,:), SSqrReal(:,:)
    real(dp), intent(in) :: filling(:,:)
    real(dp), intent(in) :: coords0(:,:)
    real(dp), intent(in) :: dqAt(:)
    integer, intent(in) :: species0(:)
    integer, intent(in) :: iNeighbor(0:,:), img2CentCell(:)
    type(TOrbitals), intent(in) :: orb
    logical, intent(in) :: tWriteTagged
    integer, intent(in) :: fdTagged
    real(dp), intent(out) :: excEnergy

  #:if WITH_ARPACK
    @:ASSERT(self%tInit)
    @:ASSERT(size(orb%nOrbAtom) == self%nAtom)
    call LinRespGrad_old(tSpin, self%nAtom, iAtomStart, eigVec,&
        & eigVal, dqAt, coords0, self%nExc, self%nStat, self%symmetry,&
        & SSqrReal, filling, species0, self%HubbardU, self%spinW, self%nEl, iNeighbor, &
        & img2CentCell, orb, tWriteTagged, fdTagged, self%fdMulliken, self%fdCoeffs, &
        & self%tGrndState, self%fdXplusY, self%fdTrans, self%fdSPTrans, &
        & self%fdTradip, self%tArnoldi, self%fdArnoldi, &
        & self%fdArnoldiDiagnosis, self%fdExc, self%tEnergyWindow, self%energyWindow, &
        & self%tOscillatorWindow, self%oscillatorWindow, excEnergy, .false.)

  #:else
    call error('Internal error: Illegal routine call to &
        &LinResp_calcExcitations')
  #:endif

  end subroutine LinResp_calcExcitations


  !> Wrapper to call linear response calculations of excitations
  !! and forces in excited states
  !! \param tSpin is this a spin-polarized calculation
  !! \param self data for the actual calculation
  !! \param iAtomStart indexing array for ground state square matrices
  !! \param eigVec ground state eigenvectors
  !! \param eigVal ground state eigenvalues
  !! \param SSqrReal square overlap matrix (must be symmetriezed)
  !! \param filling ground state occupations
  !! \param coords0 central cell atomic coordinates
  !! \param qdAt net atomic charges in ground state
  !! \param species0 chemical species of atoms in central cell
  !! \param iNeighbor index array for atoms within cutoff distances
  !! \param img2CentCell folding to central cell (not really needed
  !! for non-periodic systems)
  !! \param orb orbital data structure
  !! \param skHamCont non-SCC H0 data
  !! \param skOverCont overlap data
  !! \param tWriteTagged print tag information
  !! \param fdTagged file descriptor for tagged data
  !! \param excEnergy energy of particular excited state
  !! \param derivator method for calculating derivatives of S and H0 matrices
  !! \param rhoSqr ground state density matrix (square matrix plus spin index)
  subroutine LinResp_addGradients(tSpin, self, iAtomStart, eigVec, eigVal, SSqrReal, filling, &
      & coords0, dqAt, species0, iNeighbor, img2CentCell, orb, skHamCont, skOverCont, &
      & tWriteTagged, fdTagged, excEnergy, tForces, excgradient, derivator, rhoSqr, occNatural, &
      & naturalOrbs)
    logical, intent(in) :: tSpin
    type(linresp), intent(inout) :: self
    integer, intent(in) :: iAtomStart(:)
    real(dp), intent(in) :: eigVec(:,:,:)
    real(dp), intent(in) :: eigVal(:,:), SSqrReal(:,:)
    real(dp), intent(in) :: filling(:,:)
    real(dp), intent(in) :: coords0(:,:)
    real(dp), intent(in) :: dqAt(:)
    integer, intent(in) :: species0(:)
    integer, intent(in) :: iNeighbor(0:,:), img2CentCell(:)
    type(TOrbitals), intent(in) :: orb
    type(OSlakoCont), intent(in) :: skHamCont, skOverCont
    logical, intent(in) :: tWriteTagged
    integer, intent(in) :: fdTagged
    real(dp), intent(out) :: excenergy
    logical, intent(in)   :: tForces
    real(dp), intent(out) :: excgradient(:,:)
    class(NonSccDiff), intent(in), optional :: derivator
    real(dp), intent(in), optional :: rhoSqr(:,:,:)
    real(dp), intent(out), optional :: occNatural(:)
    real(dp), intent(out), optional :: naturalOrbs(:,:)

  #:if WITH_ARPACK

    real(dp), allocatable :: shiftPerAtom(:,:), shiftPerL(:,:,:)
    @:ASSERT(self%tInit)
    @:ASSERT(self%nAtom == size(orb%nOrbAtom))
    ! BA: SCC is currently ugly, it gives back an array with an additional
    ! dimension (spin), however, fills always the first channel only!
    ALLOCATE(shiftPerAtom(self%nAtom, 1))
    ALLOCATE(shiftPerL(orb%mShell, self%nAtom, 1))
    call getShiftPerAtom(shiftPerAtom)
    call getShiftPerL(shiftPerL)
    shiftPerAtom = shiftPerAtom + shiftPerL(1,:,:)

    call LinRespGrad_old(tSpin, self%nAtom, iAtomStart, eigVec,&
        & eigVal, dqAt, coords0, self%nExc, self%nStat, self%symmetry,&
        & SSqrReal, filling, species0, self%HubbardU, self%spinW, self%nEl, iNeighbor, &
        & img2CentCell, orb, tWriteTagged, fdTagged, self%fdMulliken, self%fdCoeffs, &
        & self%tGrndState, self%fdXplusY, self%fdTrans, self%fdSPTrans, &
        & self%fdTradip, self%tArnoldi, self%fdArnoldi, &
        & self%fdArnoldiDiagnosis, self%fdExc, self%tEnergyWindow, self%energyWindow, &
        & self%tOscillatorWindow, self%oscillatorWindow, excEnergy, tForces, shiftPerAtom(:,1), &
        & skHamCont, skOverCont, excgradient, derivator, rhoSqr, &
        & occNatural, naturalOrbs)

  #:else
    call error('Internal error: Illegal routine call to LinResp_addGradients.')
  #:endif

  end subroutine LinResp_addGradients

end module linresp_module
