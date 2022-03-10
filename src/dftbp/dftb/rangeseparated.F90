!--------------------------------------------------------------------------------------------------!
!  DFTB+: general package for performing fast atomistic simulations                                !
!  Copyright (C) 2006 - 2022  DFTB+ developers group                                               !
!                                                                                                  !
!  See the LICENSE file for terms of usage and distribution.                                       !
!--------------------------------------------------------------------------------------------------!

#:include 'common.fypp'


!> Contains range separated related routines.
module dftbp_dftb_rangeseparated
  use dftbp_common_accuracy, only : dp, tolSameDist, MinHubDiff
  use dftbp_common_environment, only : TEnvironment, globalTimers
  use dftbp_common_globalenv, only : stdOut
  use dftbp_dftb_nonscc, only : TNonSccDiff
  use dftbp_dftb_slakocont, only : TSlakoCont
  use dftbp_dftb_sparse2dense, only : blockSymmetrizeHS, symmetrizeHS, hermitianSquareMatrix
  use dftbp_io_message, only : error
  use dftbp_math_blasroutines, only : gemm
  use dftbp_math_sorting, only : index_heap_sort
  use dftbp_type_commontypes, only : TOrbitals
  implicit none

  private
  public :: TRangeSepSKTag, TRangeSepFunc, RangeSepFunc_init, getGammaPrimeValue, rangeSepTypes


  type :: TRangeSepTypesEnum

    !> Neighbour based
    integer :: neighbour = 0

    !> Threshold based
    integer :: threshold = 1

    !> Matrix based
    integer :: matrixBased = 2

  end type TRangeSepTypesEnum


  !> Container for enumerated range separation types
  type(TRangeSepTypesEnum), parameter :: rangeSepTypes = TRangeSepTypesEnum()


  !> Slater-Koster file RangeSep tag structure
  type :: TRangeSepSKTag

    !> range-separation parameter
    real(dp) :: omega

    !> CAM alpha parameter
    real(dp) :: camAlpha

    !> CAM beta parameter
    real(dp) :: camBeta

  end type TRangeSepSKTag


  !> Range-Sep module dftbp_poisson_structure
  type :: TRangeSepFunc
    private

    !> coordinates of the atom
    real(dp), allocatable :: coords(:,:)

    !> evaluated long-range gamma, Atom1, Atom2 at each geometry step
    real(dp), allocatable :: lrGammaEval(:,:)

    !> evaluated Hartree-Fock gamma, Atom1, Atom2 at each geometry step
    real(dp), allocatable :: hfGammaEval(:,:)

    !> range-separation parameter
    real(dp) :: omega

    !> CAM alpha parameter
    real(dp) :: camAlpha

    !> CAM beta parameter
    real(dp) :: camBeta

    !> True, for CAM range-separation, otherwise LC (F)
    logical :: tCam

    !> Hubbard U values for atoms
    real(dp), allocatable :: hubbu(:)

    ! Hamiltonian Screening

    !> previous hamiltonian in screening by tolerance
    real(dp), allocatable :: hprev(:,:)

    !> previous delta density matrix in screening by tolerance
    real(dp), allocatable :: dRhoPrev(:,:)

    !> Is screening initialised
    logical :: tScreeningInited

    !> threshold for screening by value
    real(dp) :: pScreeningThreshold

    !> total full-range Hartree-Fock energy
    real(dp) :: hfEnergy

    !> total long-range energy
    real(dp) :: lrEnergy

    !> Is this spin restricted (F) or unrestricted (T)
    logical :: tSpin

    !> Is this DFTB/SSR formalism
    logical :: tREKS

    !> algorithm for range separation screening
    integer :: rsAlg

    !> species of atoms
    integer, allocatable :: species(:)

  contains

    procedure :: updateCoords
    procedure :: addCamHamiltonian
    procedure :: addLrHamiltonianMatrixCmplx
    procedure :: addCamEnergy
    procedure :: addCamGradients
    procedure :: evaluateLrEnergyDirect
    procedure :: getSpecies
    procedure :: getLrGamma
    procedure :: getLrGammaDeriv
    procedure :: getHartreeFockGamma
    procedure :: getHartreeFockGammaDeriv

  end type TRangeSepFunc


contains


  !> Intitializes the range-sep module.
  subroutine RangeSepFunc_init(this, nAtom, species, hubbu, screen, omega, camAlpha, camBeta,&
      & tSpin, tREKS, rsAlg)

    !> Instance
    type(TRangeSepFunc), intent(out) :: this

    !> Number of atoms
    integer, intent(in) :: nAtom

    !> List of all atomic species
    integer, intent(in) :: species(:)

    !> Atomic hubbards
    real(dp), intent(in) :: hubbu(:)

    !> Screening threshold value
    real(dp), intent(in) :: screen

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

    !> lr-hamiltonian construction algorithm
    integer, intent(in) :: rsAlg

    call initAndAllocate(this, nAtom, hubbu, species, screen, omega, camAlpha, camBeta, rsAlg,&
        & tSpin, tREKS)
    call checkRequirements(this)

  contains


    !> Intitializes data structures and allocate storage.
    subroutine initAndAllocate(this, nAtom, hubbu, species, screen, omega, camAlpha, camBeta,&
        & rsAlg, tSpin, tREKS)

      !> Instance
      type(TRangeSepFunc), intent(out) :: this

      !> Number of atoms
      integer, intent(in) :: nAtom

      !> Hubbard U values for atoms
      real(dp), intent(in) :: hubbu(:)

      !> Species of atoms
      integer, intent(in) :: species(:)

      !> screening cutoff if using appropriate method
      real(dp), intent(in) :: screen

      !> Range separation parameter
      real(dp), intent(in) :: omega

      !> CAM alpha parameter
      real(dp), intent(in) :: camAlpha

      !> CAM beta parameter
      real(dp), intent(in) :: camBeta

      !> Algorithm for range separation
      integer, intent(in) :: rsAlg

      !> Is this spin restricted (F) or unrestricted (T)
      logical, intent(in) :: tSpin

      !> Is this DFTB/SSR formalism
      logical, intent(in) :: tREKS

      this%tScreeningInited = .false.
      this%pScreeningThreshold = screen
      this%omega = omega
      this%lrEnergy = 0.0_dp
      this%hfEnergy = 0.0_dp
      this%rsAlg = rsAlg
      this%tSpin = tSpin
      this%tREKS = tREKS

      this%camAlpha = camAlpha
      this%camBeta = camBeta

      if ((abs(this%camAlpha) < 1.0e-16_dp) .and. (abs(this%camBeta - 1.0_dp) < 1.0e-16_dp)) then
        ! apparently this is a pure LC calculation
        this%tCam = .false.
      else
        this%tCam = .true.
      end if

      if (this%tREKS .and. this%tCam) then
        call error("General CAM functionals not currently implemented for REKS.")
      end if

      allocate(this%coords(3, nAtom))
      allocate(this%lrGammaEval(nAtom, nAtom))
      allocate(this%hfGammaEval(nAtom, nAtom))

      allocate(this%hubbu(size(hubbu)))
      this%hubbu(:) = hubbu
      allocate(this%species(nAtom))
      this%species(:) = species

    end subroutine initAndAllocate


    !> Test for option consistency
    subroutine checkRequirements(this)

      !> instance
      type(TRangeSepFunc), intent(inout) :: this

      ! Check for current restrictions
      if (this%tSpin .and. this%rsAlg == rangeSepTypes%threshold) then
        call error("Spin-unrestricted calculation for thresholded range separation not yet&
            & implemented!")
      end if

      if (this%tREKS .and. this%rsAlg == rangeSepTypes%threshold) then
        call error("REKS calculation with thresholded range separation not yet implemented!")
      end if

      if (.not. any([rangeSepTypes%neighbour, rangeSepTypes%threshold,&
            & rangeSepTypes%matrixBased] == this%rsAlg)) then
        call error("Unknown algorithm for screening the exchange in range separation!")
      end if

    end subroutine checkRequirements

  end subroutine RangeSepFunc_init


  !> Updates the rangeSep module on coordinate change.
  subroutine updateCoords(this, coords)

    !> Class instance
    class(TRangeSepFunc), intent(inout) :: this

    !> List of atomic coordinates
    real(dp), intent(in) :: coords(:,:)

    integer :: nAtom, iAtom1, iAtom2, iSp1, iSp2
    real(dp) :: dist

    this%coords(:,:) = coords
    nAtom = size(this%species)
    do iAtom1 = 1, nAtom
      do iAtom2 = 1, iAtom1
        iSp1 = this%species(iAtom1)
        iSp2 = this%species(iAtom2)
        dist = norm2(this%coords(:, iAtom1) - this%coords(:, iAtom2))
        this%lrGammaEval(iAtom1, iAtom2) = getAnalyticalGammaValue(this, iSp1, iSp2, dist)
        this%lrGammaEval(iAtom2, iAtom1) = this%lrGammaEval(iAtom1, iAtom2)
      end do
    end do

    if (this%tCam) then
      do iAtom1 = 1, nAtom
        do iAtom2 = 1, iAtom1
          iSp1 = this%species(iAtom1)
          iSp2 = this%species(iAtom2)
          dist = norm2(this%coords(:, iAtom1) - this%coords(:, iAtom2))
          this%hfGammaEval(iAtom1, iAtom2) = getAnalyticalHartreeFockGammaValue(this, iSp1, iSp2,&
              & dist)
          this%hfGammaEval(iAtom2, iAtom1) = this%hfGammaEval(iAtom1, iAtom2)
        end do
      end do
    end if

    if (this%tScreeningInited) then
      this%hprev(:,:) = 0.0_dp
      this%dRhoPrev(:,:) = 0.0_dp
      this%lrEnergy = 0.0_dp
      this%hfEnergy = 0.0_dp
    end if

  end subroutine updateCoords


  !> Interface routine for adding CAM range-separated contributions to the Hamiltonian.
  subroutine addCamHamiltonian(this, env, densSqr, over, iNeighbour, nNeighbourLC, iSquare, iPair,&
      & orb, HH, overlap)

    !> Class instance
    class(TRangeSepFunc), intent(inout) :: this

    !> Environment settings
    type(TEnvironment), intent(inout) :: env

    ! Neighbour based screening

    !> Square (unpacked) density matrix
    real(dp), intent(in), target :: densSqr(:,:)

    !> Sparse (packed) overlap matrix.
    real(dp), intent(in) :: over(:)

    !> Neighbour indices.
    integer, intent(in) :: iNeighbour(0:,:)

    !> Nr. of neighbours for each atom.
    integer, intent(in) :: nNeighbourLC(:)

    !> Position of each atom in the rows/columns of the square matrices. Shape: (nAtom)
    integer, intent(in) :: iSquare(:)

    !> iPair Position of each (neighbour, atom) pair in the sparse matrix.
    !> Shape: (0:maxNeighbour, nAtom)
    integer, intent(in) :: iPair(0:,:)

    !> Orbital information.
    type(TOrbitals), intent(in) :: orb

    !> Square (unpacked) Hamiltonian to be updated.
    real(dp), intent(inout), target :: HH(:,:)

    ! Threshold based screening

    !> square real overlap matrix
    real(dp), intent(in) :: overlap(:,:)

    call env%globalTimer%startTimer(globalTimers%rangeSeparatedH)

    ! always add the LR contribution
    call addLrHamiltonian(this, densSqr, over, iNeighbour, nNeighbourLC, iSquare, iPair, orb, HH,&
        & overlap)

    ! If xc-functional is more general, also add the full-range Hartree-Fock part.
    ! For pure LC, camAlpha would be zero anyway, but we want to save as much time as possible.
    if (this%tCam) then
      call addHartreeFockHamiltonian(this, densSqr, over, iNeighbour, nNeighbourLC, iSquare, iPair,&
          & orb, HH, overlap)
    end if

    call env%globalTimer%stopTimer(globalTimers%rangeSeparatedH)

  end subroutine addCamHamiltonian


  !> Interface routine for adding LC range-separated contributions to the Hamiltonian.
  subroutine addLrHamiltonian(this, densSqr, over, iNeighbour, nNeighbourLC, iSquare, iPair, orb,&
      & HH, overlap)

    !> Instance
    type(TRangeSepFunc), intent(inout) :: this

    ! Neighbour based screening

    !> Square (unpacked) density matrix
    real(dp), intent(in), target :: densSqr(:,:)

    !> Sparse (packed) overlap matrix.
    real(dp), intent(in) :: over(:)

    !> Neighbour indices.
    integer, intent(in) :: iNeighbour(0:,:)

    !> Nr. of neighbours for each atom.
    integer, intent(in) :: nNeighbourLC(:)

    !> Position of each atom in the rows/columns of the square matrices. Shape: (nAtom)
    integer, intent(in) :: iSquare(:)

    !> iPair Position of each (neighbour, atom) pair in the sparse matrix.
    !> Shape: (0:maxNeighbour, nAtom)
    integer, intent(in) :: iPair(0:,:)

    !> Orbital information.
    type(TOrbitals), intent(in) :: orb

    !> Square (unpacked) Hamiltonian to be updated.
    real(dp), intent(inout), target :: HH(:,:)

    ! Threshold based screening

    !> square real overlap matrix
    real(dp), intent(in) :: overlap(:,:)

    select case(this%rsAlg)
    case (rangeSepTypes%threshold)
      call addLrHamiltonianThreshold(this, overlap, densSqr, iNeighbour, nNeighbourLC, iSquare, HH,&
          & orb)
    case (rangeSepTypes%neighbour)
      call addLrHamiltonianNeighbour(this, densSqr, over, iNeighbour, nNeighbourLC, iSquare, iPair,&
          & orb, HH)
    case (rangeSepTypes%matrixBased)
      call addLrHamiltonianMatrix(this, iSquare, overlap, densSqr, HH)
    end select

  end subroutine addLrHamiltonian


  !> Adds the LR-exchange contribution to hamiltonian using the thresholding algorithm.
  subroutine addLrHamiltonianThreshold(this, overlap, deltaRho, iNeighbour, nNeighbourLC, iSquare,&
      & hamiltonian, orb)

    !> Class instance
    type(TRangeSepFunc), intent(inout) :: this

    !> Square real overlap matrix
    real(dp), intent(in) :: overlap(:,:)

    !> Square density matrix (deltaRho in DFTB terms)
    real(dp), intent(in) :: deltaRho(:,:)

    !> Neighbour indices.
    integer, dimension(0:,:), intent(in) :: iNeighbour

    !> Nr. of neighbours for each atom.
    integer, dimension(:), intent(in) :: nNeighbourLC

    !> Mapping atom_number -> number of the first basis function of the atomic block atom_number
    integer, intent(in) :: iSquare(:)

    !> Current hamiltonian
    real(dp), intent(inout) :: hamiltonian(:,:)

    !> Orbital information
    type(TOrbitals), intent(in) :: orb

    real(dp), allocatable :: tmpOvr(:,:), tmpDRho(:,:), testOvr(:,:), tmpDDRho(:,:), tmpDHam(:,:)
    integer, allocatable :: ovrInd(:,:)
    integer, parameter :: descLen = 3, iStart = 1, iEnd = 2, iNOrb = 3

    call allocateAndInit()
    call evaluateHamiltonian(tmpDHam)
    this%hprev(:,:) = this%hprev + tmpDHam
    hamiltonian(:,:) = hamiltonian + this%camBeta * this%hprev
    this%lrEnergy = this%lrEnergy + evaluateEnergy(this%hprev, tmpDRho)

  contains

    !> allocate and initialise some necessary arrays
    subroutine allocateAndInit()

      integer :: matrixSize, nAtom
      real(dp) :: tmp
      integer :: iAtMu, iAtNu, iNeigh

      matrixSize = size(hamiltonian, dim=1)
      nAtom = size(this%species)
      allocate(tmpOvr(matrixSize, matrixSize))
      allocate(tmpDHam(matrixSize, matrixSize))
      allocate(tmpDRho(matrixSize, matrixSize))
      allocate(tmpDDRho(matrixSize, matrixSize))
      allocate(testOvr(nAtom, nAtom))
      allocate(ovrInd(nAtom, nAtom))
      tmpOvr(:,:) = overlap
      call blockSymmetrizeHS(tmpOvr, iSquare)
      tmpDRho(:,:) = deltaRho
      call symmetrizeHS(tmpDRho)
      tmpDHam(:,:) = 0.0_dp
      call checkAndInitScreening(this, matrixSize, tmpDRho)
      tmpDDRho(:,:) = tmpDRho - this%dRhoPrev
      this%dRhoPrev(:,:) = tmpDRho

      testOvr(:,:) = 0.0_dp
      do iAtMu = 1, nAtom
        do iNeigh = 0, nNeighbourLC(iAtMu)
          iAtNu = iNeighbour(iNeigh, iAtMu)
          tmp = maxval(abs(tmpOvr(iSquare(iAtMu) : iSquare(iAtMu + 1) - 1,&
              & iSquare(iAtNu) : iSquare(iAtNu + 1) - 1)))
          testOvr(iAtMu, iAtNu) = tmp
          testOvr(iAtNu, iAtMu) = tmp
        end do
      end do
      do iAtMu = 1, nAtom
        call index_heap_sort(ovrInd(iAtMu,:), testOvr(iAtMu,:))
      end do

    end subroutine allocateAndInit


    !> Evaluate the update to Hamiltonian due to change in the DM.
    pure subroutine evaluateHamiltonian(tmpDHam)

      !> Update for the old hamiltonian on exit
      real(dp), intent(out) :: tmpDHam(:,:)

      integer :: nAtom
      real(dp) :: pbound, prb
      real(dp) :: tmpvec1(orb%mOrb), tmpvec2(orb%mOrb)
      real(dp) :: tmp, tstbound, gammabatch, gammabatchtmp
      integer :: iAtMu, iAtNu, iAt1, iAt2, iSp1, iSp2, nOrb1, nOrb2
      integer :: kk, ll, jj, ii, mu, nu
      integer, dimension(descLen) :: desc1, desc2, descM, descN

      nAtom = size(this%species)

      pbound = maxval(abs(tmpDDRho))
      tmpDHam = 0.0_dp
      loopMu: do iAtMu = 1, nAtom
        descM = getDescriptor(iAtMu, iSquare)
        loopKK: do kk = 1, nAtom
          iAt1 = ovrInd(iAtMu, nAtom + 1 - kk)
          desc1 = getDescriptor(iAt1, iSquare)
          iSp1 = this%species(iAt1)
          nOrb1 = orb%nOrbSpecies(iSp1)
          prb = pbound * testOvr(iAt1, iAtMu)
          if (abs(prb) < this%pScreeningThreshold) then
            exit loopKK
          end if
          loopNu: do iAtNu = 1, iAtMu
            descN = getDescriptor(iAtNu, iSquare)
            gammabatchtmp = this%lrGammaEval(iAtMu, iAtNu) + this%lrGammaEval(iAt1, iAtNu)
            loopLL: do ll = 1, nAtom
              iAt2 = ovrInd(iAtNu, nAtom + 1 - ll)
              iSp2 = this%species(iAt2)
              nOrb2 = orb%nOrbSpecies(iSp2)
              tstbound = prb * testOvr(iAt2, iAtNu)
              if (abs(tstbound) < this%pScreeningThreshold) then
                exit loopLL
              end if
              desc2 = getDescriptor(iAt2, iSquare)
              gammabatch = (this%lrGammaEval(iAtMu, iAt2) + this%lrGammaEval(iAt1, iAt2)&
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

      !> Instance
      type(TRangeSepFunc), intent(inout) :: this

      !> linear dimension of matrix
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

  end subroutine addLrHamiltonianThreshold


  !> Updates the Hamiltonian with the range separated contribution.
  subroutine addLrHamiltonianNeighbour(this, densSqr, over, iNeighbour, nNeighbourLC, iSquare,&
      & iPair, orb, HH)

    !> instance of object
    type(TRangeSepFunc), intent(inout) :: this

    !> Square (unpacked) density matrix
    real(dp), intent(in), target :: densSqr(:,:)

    !> Sparse (packed) overlap matrix.
    real(dp), intent(in) :: over(:)

    !> Neighbour indices.
    integer, intent(in) :: iNeighbour(0:,:)

    !> Nr. of neighbours for each atom.
    integer, intent(in) :: nNeighbourLC(:)

    !> Position of each atom in the rows/columns of the square matrices. Shape: (nAtom)
    integer, intent(in) :: iSquare(:)

    !> Position of each (neighbour, atom) pair in the sparse matrix. Shape: (0:maxNeighbour, nAtom)
    integer, intent(in) :: iPair(0:,:)

    !> Orbital information.
    type(TOrbitals), intent(in) :: orb

    !> Square (unpacked) Hamiltonian to be updated.
    real(dp), intent(inout), target :: HH(:,:)

    integer, parameter :: descLen = 3, iStart = 1, iEnd = 2, iNOrb = 3
    real(dp), dimension(orb%mOrb**2), target :: Sma, Sam, Snb, Sbn
    real(dp), dimension(orb%mOrb**2), target :: Pab, Pmb, Pan, Pmn
    real(dp), dimension(:,:), pointer :: pSma, pSam, pSnb, pSbn
    real(dp), dimension(:,:), pointer :: pPab, pPmb, pPan, pPmn
    real(dp) :: gamma1, gamma2, gammaTot
    integer :: nAtom
    integer :: iAtM, iAtN, iAtA, iAtB, iNeighN, iNeighA
    integer, dimension(descLen) :: descA, descB, descM, descN
    real(dp), dimension(:,:), allocatable, target :: tmpDRho
    real(dp), dimension(:,:), allocatable, target :: tmpHH

    call allocateAndInit(tmpHH, tmpDRho)
    call evaluateHamiltonian()
    call symmetrizeHS(tmpHH)

    HH(:,:) = HH + this%camBeta * tmpHH
    this%lrEnergy = this%lrEnergy + evaluateEnergy(tmpHH, tmpDRho)

  contains

    !> Allocate storage for mapping 1D<->2D array sections
    subroutine allocateAndInit(tmpHH, tmpDRho)

      !> density matrix case
      real(dp), dimension(:,:), allocatable, target, intent(inout) :: tmpDRho

      !> hamiltonian matrix case
      real(dp), dimension(:,:), allocatable, target, intent(inout) :: tmpHH

      allocate(tmpHH(size(HH, dim=1), size(HH, dim=2)))
      tmpHH(:,:) = 0.0_dp
      allocate(tmpDRho(size(densSqr, dim=1), size(densSqr, dim=1)))
      tmpDRho(:,:) = densSqr
      call symmetrizeHS(tmpDRho)

    end subroutine allocateAndInit


    !> actually evaluate the neighbour based cut-off hamiltonian
    subroutine evaluateHamiltonian()

      nAtom = size(this%species)

      loopN: do iAtN = 1, nAtom
        descN = getDescriptor(iAtN, iSquare)
        loopB: do iNeighN = 0, nNeighbourLC(iAtN)
          iAtB = iNeighbour(iNeighN, iAtN)
          descB = getDescriptor(iAtB, iSquare)
          call copyOverlapBlock(iAtN, iNeighN, descN(iNOrb), descB(iNOrb), Sbn, pSbn)
          call transposeBlock(pSbn, Snb, pSnb)
          loopA: do iAtA = 1, nAtom
            descA = getDescriptor(iAtA, iSquare)
            call copyDensityBlock(descA, descB, Pab, pPab)
            call copyDensityBlock(descA, descN, Pan, pPan)
            gamma1 = this%lrGammaEval(iAtA, iAtN) + this%lrGammaEval(iAtA, iAtB)
            loopM: do iNeighA = 0, nNeighbourLC(iAtA)
              iAtM = iNeighbour(iNeighA, iAtA)
              descM = getDescriptor(iAtM, iSquare)
              call copyOverlapBlock(iAtA, iNeighA, descA(iNOrb), descM(iNOrb), Sma, pSma)
              call transposeBlock(pSma, Sam, pSam)
              gamma2 = this%lrGammaEval(iAtM, iAtN) + this%lrGammaEval(iAtM, iAtB)
              gammaTot = gamma1 + gamma2
              !
              if (iAtM >= iAtN) then
                call updateHamiltonianBlock(descM, descN, pSma, pSbn, pPab)
              end if
              if (iAtA >= iAtN .and. iAtM /= iAtA) then
                call copyDensityBlock(descM, descB, Pmb, pPmb)
                call updateHamiltonianBlock(descA, descN, pSam, pSbn, pPmb)
              end if
              if (iAtM >= iAtB .and. iAtN /= iAtB) then
                call updateHamiltonianBlock(descM, descB, pSma, pSnb, pPan)
              end if
              if (iAtA >= iAtB .and. iAtM /= iAtA .and. iAtN /= iAtB) then
                call copyDensityBlock(descM, descN, Pmn, pPmn)
                call updateHamiltonianBlock(descA, descB, pSam, pSnb, pPmn)
              end if
            end do loopM
          end do loopA
        end do loopB
      end do loopN

    end subroutine evaluateHamiltonian


    !> copy atom block from sparse matrix
    pure subroutine copyOverlapBlock(iAt, iNeigh, iNOrbAt, iNOrbNeigh, localBlock, pLocalBlock)

      !> Atom for which this is a neighbour
      integer, intent(in) :: iAt

      !> Number of neighbour for this block
      integer, intent(in) :: iNeigh

      !> Number of orbitals on iAt
      integer, intent(in) :: iNOrbAt

      !> Number of orbitals on neighbour atom
      integer, intent(in) :: iNOrbNeigh

      !> local block
      real(dp), dimension(:), target, intent(inout) :: localBlock

      !> Pointer to local block
      real(dp), dimension(:,:), pointer, intent(out) :: pLocalBlock

      integer :: ind

      ind = iPair(iNeigh, iAt) + 1
      localBlock(1:iNOrbNeigh*iNOrbAt) = over(ind:ind+iNOrbNeigh*iNOrbAt-1)
      pLocalBlock(1:iNOrbNeigh, 1:iNOrbAt) => localBlock(1:iNOrbNeigh*iNOrbAt)

    end subroutine copyOverlapBlock


    !> copy a density matrix block from sparse matrix
    pure subroutine copyDensityBlock(desc1, desc2, localBlock, pLocalBlock)

      !> start, end and range of first block
      integer, dimension(descLen), intent(in) :: desc1

      !> start, end and range of second block
      integer, dimension(descLen), intent(in) :: desc2

      !> local block in 1D format
      real(dp), dimension(:), target, intent(inout) :: localBlock

      !> Pointer to local block
      real(dp), dimension(:,:), pointer, intent(out) :: pLocalBlock

      pLocalBlock(1:desc1(iNOrb), 1:desc2(iNOrb)) => localBlock(1:desc1(iNOrb) * desc2(iNOrb))
      pLocalBlock(:,:) = tmpDRho(desc1(iStart):desc1(iEnd), desc2(iStart):desc2(iEnd))

    end subroutine copyDensityBlock


    !> Transpose a block
    pure subroutine transposeBlock(orig, localBlock, pLocalBlock)

      !> Original matrix block
      real(dp), dimension(:,:), intent(in) :: orig

      !> local copy in 1D
      real(dp), dimension(:), target, intent(out) :: localBlock

      !> pointer to local copy
      real(dp), dimension(:,:), pointer, intent(out) :: pLocalBlock

      pLocalBlock(1:size(orig, dim=2), 1:size(orig, dim=1)) => localBlock(1:size(orig))
      pLocalBlock = transpose(orig)


    end subroutine transposeBlock


    !> Add a contribution to a Hamiltonian block
    subroutine updateHamiltonianBlock(descM, descN, pSma, pSbN, pPab)

      !> start, end and range of row
      integer, dimension(descLen), intent(in) :: descM

      !> start, end and range of column
      integer, dimension(descLen), intent(in) :: descN

      !> First overlap block
      real(dp), dimension(:,:), pointer, intent(in) :: pSma

      !> Second overlap block
      real(dp), dimension(:,:), pointer, intent(in) :: pSbN

      !> density matrix block
      real(dp), dimension(:,:), pointer, intent(in) :: pPab

      real(dp), dimension(:,:), pointer :: pHmn

      pHmn => tmpHH(descM(iStart):descM(iEnd), descN(iStart):descN(iEnd))
      if (this%tSpin .or. this%tREKS) then
        pHmn(:,:) = pHmn - 0.25_dp * gammaTot * matmul(matmul(pSma, pPab), pSbn)
      else
        pHmn(:,:) = pHmn - 0.125_dp * gammaTot * matmul(matmul(pSma, pPab), pSbn)
      end if

    end subroutine updateHamiltonianBlock

  end subroutine addLrHamiltonianNeighbour


  !> Update Hamiltonian with long-range contribution using matrix-matrix multiplications
  !>
  !> The routine provides a matrix-matrix multiplication based implementation of
  !> the 3rd term in Eq. 26 in https://doi.org/10.1063/1.4935095
  !>
  subroutine addLrHamiltonianMatrixCmplx(this, iSquare, overlap, densSqr, HH)

    !> Class instance
    class(TRangeSepFunc), intent(inout) :: this

    !> Position of each atom in the rows/columns of the square matrices. Shape: (nAtom)
    integer, dimension(:), intent(in) :: iSquare

    !> Square (unpacked) overlap matrix.
    complex(dp), intent(in) :: overlap(:,:)

    !> Square (unpacked) density matrix
    complex(dp), intent(in) :: densSqr(:,:)

    !> Square (unpacked) Hamiltonian to be updated.
    complex(dp), intent(inout) :: HH(:,:)

    complex(dp), allocatable :: Smat(:,:)
    complex(dp), allocatable :: Dmat(:,:)
    real(dp), allocatable :: lrGammaAO(:,:)
    complex(dp), allocatable :: gammaCmplx(:,:)
    complex(dp), allocatable :: Hlr(:,:)

    integer :: nOrb

    nOrb = size(overlap,dim=1)

    allocate(Smat(nOrb,nOrb))
    allocate(Dmat(nOrb,nOrb))
    allocate(lrGammaAO(nOrb,nOrb))
    allocate(gammaCmplx(nOrb,nOrb))
    allocate(Hlr(nOrb,nOrb))

    call allocateAndInit(this, iSquare, overlap, densSqr, HH, Smat, Dmat, lrGammaAO, gammaCmplx)

    call evaluateHamiltonian(this, Smat, Dmat, gammaCmplx, Hlr)

    HH(:,:) = HH + this%camBeta * Hlr

    this%lrenergy = this%lrenergy + 0.5_dp * real(sum(Dmat * Hlr), dp)

  contains

    subroutine allocateAndInit(this, iSquare, overlap, densSqr, HH, Smat, Dmat, lrGammaAO,&
        & gammaCmplx)

      !> instance
      type(TRangeSepFunc), intent(inout) :: this

      !> Position of each atom in the rows/columns of the square matrices. Shape: (nAtom)
      integer, dimension(:), intent(in) :: iSquare

      !> Square (unpacked) overlap matrix.
      complex(dp), intent(in) :: overlap(:,:)

      !> Square (unpacked) density matrix
      complex(dp), intent(in) :: densSqr(:,:)

      !> Square (unpacked) Hamiltonian to be updated.
      complex(dp), intent(inout) :: HH(:,:)

      !> Symmetrized square overlap matrix
      complex(dp), intent(out) :: Smat(:,:)

      !> Symmetrized square density matrix
      complex(dp), intent(out) :: Dmat(:,:)

      !> Symmetrized long-range gamma matrix
      real(dp), intent(out) :: lrGammaAO(:,:)

      !> Symmetrized long-range gamma matrix
      complex(dp), intent(out) :: gammaCmplx(:,:)

      integer :: nAtom, iAt, jAt

      nAtom = size(this%lrGammaEval, dim=1)

      !! Symmetrize Hamiltonian, overlap, density matrices
      call hermitianSquareMatrix(HH)
      Smat(:,:) = overlap
      call hermitianSquareMatrix(Smat)
      Dmat(:,:) = densSqr
      call hermitianSquareMatrix(Dmat)

      ! Get long-range gamma variable
      lrGammaAO(:,:) = 0.0_dp
      do iAt = 1, nAtom
        do jAt = 1, nAtom
          lrGammaAO(iSquare(jAt):iSquare(jAt+1)-1,iSquare(iAt):iSquare(iAt+1)-1) =&
              & this%lrGammaEval(jAt, iAt)
        end do
      end do
      gammaCmplx = lrGammaAO

    end subroutine allocateAndInit


    subroutine evaluateHamiltonian(this, Smat, Dmat, gammaCmplx, Hlr)

      !> instance
      type(TRangeSepFunc), intent(inout) :: this

      !> Symmetrized square overlap matrix
      complex(dp), intent(in) :: Smat(:,:)

      !> Symmetrized square density matrix
      complex(dp), intent(in) :: Dmat(:,:)

      !> Symmetrized long-range gamma matrix
      complex(dp), intent(in) :: gammaCmplx(:,:)

      !> Symmetrized long-range Hamiltonian matrix
      complex(dp), intent(out) :: Hlr(:,:)

      complex(dp), allocatable :: Hmat(:,:)
      complex(dp), allocatable :: tmpMat(:,:)

      integer :: nOrb

      nOrb = size(Smat,dim=1)

      allocate(Hmat(nOrb,nOrb))
      allocate(tmpMat(nOrb,nOrb))

      Hlr(:,:) = cmplx(0.0_dp,0.0_dp,dp)

      call gemm(tmpMat, Smat, Dmat)
      call gemm(Hlr, tmpMat, Smat)
      Hlr(:,:) = Hlr * gammaCmplx

      tmpMat(:,:) = tmpMat * gammaCmplx
      call gemm(Hlr, tmpMat, Smat, alpha=(1.0_dp,0.0_dp), beta=(1.0_dp,0.0_dp))

      Hmat(:,:) = Dmat * gammaCmplx
      call gemm(tmpMat, Smat, Hmat)
      call gemm(Hlr, tmpMat, Smat, alpha=(1.0_dp,0.0_dp), beta=(1.0_dp,0.0_dp))

      call gemm(tmpMat, Dmat, Smat)
      tmpMat(:,:) = tmpMat * gammaCmplx
      call gemm(Hlr, Smat, tmpMat, alpha=(1.0_dp,0.0_dp), beta=(1.0_dp,0.0_dp))

      if (this%tSpin) then
        Hlr(:,:) = -0.25_dp * Hlr
      else
        Hlr(:,:) = -0.125_dp * Hlr
      end if

    end subroutine evaluateHamiltonian

  end subroutine addLrHamiltonianMatrixCmplx


  !> Update Hamiltonian with long-range contribution using matrix-matrix multiplications
  !>
  !> The routine provides a matrix-matrix multiplication based implementation of
  !> the 3rd term in Eq. 26 in https://doi.org/10.1063/1.4935095
  !>
  subroutine addLrHamiltonianMatrix(this, iSquare, overlap, densSqr, HH)

    !> Class instance
    type(TRangeSepFunc), intent(inout) :: this

    !> Position of each atom in the rows/columns of the square matrices. Shape: (nAtom)
    integer, dimension(:), intent(in) :: iSquare

    !> Square (unpacked) overlap matrix.
    real(dp), intent(in) :: overlap(:,:)

    !> Square (unpacked) density matrix
    real(dp), intent(in) :: densSqr(:,:)

    !> Square (unpacked) Hamiltonian to be updated.
    real(dp), intent(inout) :: HH(:,:)

    real(dp), allocatable :: Smat(:,:)
    real(dp), allocatable :: Dmat(:,:)
    real(dp), allocatable :: lrGammaAO(:,:)
    real(dp), allocatable :: Hlr(:,:)

    integer :: nOrb

    nOrb = size(overlap,dim=1)

    allocate(Smat(nOrb, nOrb))
    allocate(Dmat(nOrb, nOrb))
    allocate(lrGammaAO(nOrb, nOrb))
    allocate(Hlr(nOrb, nOrb))

    call allocateAndInit(this, iSquare, overlap, densSqr, HH, Smat, Dmat, lrGammaAO)
    call evaluateHamiltonian(this, Smat, Dmat, lrGammaAO, Hlr)
    HH(:,:) = HH + this%camBeta * Hlr
    this%lrEnergy = this%lrEnergy + 0.5_dp * sum(Dmat * Hlr)

  contains

    !> Set up storage and get orbital-by-orbital gamma matrix
    subroutine allocateAndInit(this, iSquare, overlap, densSqr, HH, Smat, Dmat, lrGammaAO)

      !> Class instance
      type(TRangeSepFunc), intent(inout) :: this

      !> Position of each atom in the rows/columns of the square matrices. Shape: (nAtom)
      integer, dimension(:), intent(in) :: iSquare

      !> Square (unpacked) overlap matrix.
      real(dp), intent(in) :: overlap(:,:)

      !> Square (unpacked) density matrix
      real(dp), intent(in) :: densSqr(:,:)

      !> Square (unpacked) Hamiltonian to be updated.
      real(dp), intent(inout) :: HH(:,:)

      !> Symmetrized square overlap matrix
      real(dp), intent(out) :: Smat(:,:)

      !> Symmetrized square density matrix
      real(dp), intent(out) :: Dmat(:,:)

      !> Symmetrized long-range gamma matrix
      real(dp), intent(out) :: lrGammaAO(:,:)

      integer :: nAtom, iAt, jAt

      nAtom = size(this%lrGammaEval,dim=1)

      ! Symmetrize Hamiltonian, overlap, density matrices
      call symmetrizeHS(HH)
      Smat(:,:) = overlap
      call symmetrizeHS(Smat)
      Dmat(:,:) = densSqr
      call symmetrizeHS(Dmat)

      ! Get long-range gamma variable
      lrGammaAO(:,:) = 0.0_dp
      do iAt = 1, nAtom
        do jAt = 1, nAtom
          lrGammaAO(iSquare(jAt):iSquare(jAt+1)-1,iSquare(iAt):iSquare(iAt+1)-1) =&
              & this%lrGammaEval(jAt, iAt)
        end do
      end do

    end subroutine allocateAndInit


    !> Evaluate the hamiltonian using GEMM operations
    subroutine evaluateHamiltonian(this, Smat, Dmat, lrGammaAO, Hlr)

      !> Class instance
      type(TRangeSepFunc), intent(inout) :: this

      !> Symmetrized square overlap matrix
      real(dp), intent(in) :: Smat(:,:)

      !> Symmetrized square density matrix
      real(dp), intent(in) :: Dmat(:,:)

      !> Symmetrized long-range gamma matrix
      real(dp), intent(in) :: lrGammaAO(:,:)

      !> Symmetrized long-range Hamiltonian matrix
      real(dp), intent(out) :: Hlr(:,:)

      real(dp), allocatable :: Hmat(:,:)
      real(dp), allocatable :: tmpMat(:,:)

      integer :: nOrb

      nOrb = size(Smat, dim=1)

      allocate(Hmat(nOrb, nOrb))
      allocate(tmpMat(nOrb, nOrb))

      Hlr(:,:) = 0.0_dp

      call gemm(tmpMat, Smat, Dmat)
      call gemm(Hlr, tmpMat, Smat)
      Hlr(:,:) = Hlr * lrGammaAO

      tmpMat(:,:) = tmpMat * lrGammaAO
      call gemm(Hlr, tmpMat, Smat, alpha=1.0_dp, beta=1.0_dp)

      Hmat(:,:) = Dmat * lrGammaAO
      call gemm(tmpMat, Smat, Hmat)
      call gemm(Hlr, tmpMat, Smat, alpha=1.0_dp, beta=1.0_dp)

      call gemm(tmpMat, Dmat, Smat)
      tmpMat(:,:) = tmpMat * lrGammaAO
      call gemm(Hlr, Smat, tmpMat, alpha=1.0_dp, beta=1.0_dp)

      if (this%tSpin .or. this%tREKS) then
        Hlr(:,:) = -0.25_dp * Hlr
      else
        Hlr(:,:) = -0.125_dp * Hlr
      end if

    end subroutine evaluateHamiltonian

  end subroutine addLRHamiltonianMatrix


  !> Add the CAM-energy contribution to the total energy.
  subroutine addCamEnergy(this, energy)

    !> Class instance
    class(TRangeSepFunc), intent(inout) :: this

    !> Total energy
    real(dp), intent(inout) :: energy

    call addLrEnergy(this, energy)
    call addHartreeFockEnergy(this, energy)

  end subroutine addCamEnergy


  !> Add the LR-energy contribution to the total energy.
  subroutine addLrEnergy(this, energy)

    !> Instance
    type(TRangeSepFunc), intent(inout) :: this

    !> Total energy
    real(dp), intent(inout) :: energy

    energy = energy + this%camBeta * this%lrEnergy

    ! hack for spin unrestricted calculation
    this%lrEnergy = 0.0_dp

  end subroutine addLrenergy


  !> Finds location of relevant atomic block indices in a dense matrix.
  pure function getDescriptor(iAt, iSquare) result(desc)

    !> relevant atom
    integer, intent(in) :: iAt

    !> indexing array for start of atom orbitals
    integer, intent(in) :: iSquare(:)

    !> resulting location ranges
    integer :: desc(3)

    desc(:) = [iSquare(iAt), iSquare(iAt + 1) - 1, iSquare(iAt + 1) - iSquare(iAt)]

  end function getDescriptor


  !> Evaluates energy from triangles of the Hamiltonian and density matrix.
  pure function evaluateEnergy(hamiltonian, densityMat) result(energy)

    !> Hamiltonian matrix
    real(dp), intent(in) :: hamiltonian(:,:)

    !> Density matrix
    real(dp), intent(in) :: densityMat(:,:)

    !> Resulting energy
    real(dp) :: energy

    integer :: mu

    energy = 0.0_dp
    do mu = 1, size(hamiltonian, dim=2)
      energy = energy + hamiltonian(mu, mu) * densityMat(mu, mu)&
          & + 2.0_dp * sum(hamiltonian(mu + 1 :, mu) * densityMat(mu + 1 :, mu))
    end do
    energy = 0.5_dp * energy

  end function evaluateEnergy


  !> Calculates analytical long-range gamma.
  function getAnalyticalGammaValue(this, Sp1, Sp2, dist) result(gamma)

    !> Instance
    type(TRangeSepFunc), intent(in) :: this

    !> First species
    integer, intent(in) :: Sp1

    !> Second species
    integer, intent(in) :: Sp2

    !> Distance between atoms
    real(dp), intent(in) :: dist

    !> Resulting gamma
    real(dp) :: gamma

    real(dp) :: tauA, tauB, omega
    real(dp) :: prefac, tmp, tmp2, tau

    tauA = 3.2_dp * this%hubbu(Sp1)
    tauB = 3.2_dp * this%hubbu(Sp2)
    omega = this%omega

    if (dist < tolSameDist) then
      ! on-site case
      if (abs(tauA - tauB) < MinHubDiff) then
        tau = 0.5_dp * (tauA + tauB)
        tmp = 5.0_dp * tau**6 + 15.0_dp * tau**4 * omega**2 - 5.0_dp * tau**2 * omega**4 + omega**6
        tmp = tmp * 0.0625_dp / tau**5 - omega
        tmp = tmp * tau**8 / (tau**2 - omega**2)**4
        gamma = tau * 0.3125_dp - tmp
      else
        call error("Error(RangeSep): R = 0, Ua != Ub")
      end if
    else
      ! off-site case, Ua == Ub
      if (abs(tauA - tauB) < MinHubDiff ) then
        tauA = 0.5_dp * (tauA + tauB)
        tmp2 = ((dist * tauA)**3 / 48.0_dp + 0.1875_dp * (dist * tauA)**2 +&
            & 0.6875_dp * (dist * tauA) + 1.0_dp) * exp(-tauA * dist) / dist
        tmp = -tauA**8 / (tauA**2 - omega**2)**4 * (tmp2 + exp(-tauA*dist) * &
            & (dist**2 * (3.0_dp * tauA**4 * omega**4 - 3.0_dp * tauA**6 * omega**2 - &
            & tauA**2 * omega**6) + dist * (15.0_dp * tauA**3 * omega**4 - &
            & 21.0_dp * tauA**5 * omega**2 - 3.0_dp * tauA * omega**6) + &
            & (15.0_dp * tauA**2 * omega**4 - 45.0_dp * tauA**4 * omega**2 - &
            & 3.0_dp * omega**6)) / (48.0_dp * tauA**5))
        gamma = 1.0_dp/dist - tmp2 - (tauA**8 / (tauA**2 - omega**2)**4 *&
            & exp(-omega * dist) / dist + tmp)
      else
        ! off-site, Ua != Ub
        prefac = tauA**4 / (tauA * tauA - omega * omega)**2
        prefac = prefac * tauB**4 / (tauB * tauB - omega * omega)**2
        prefac = prefac * exp(-omega * dist) / dist
        tmp = prefac&
            & - getYGammaSubPart(tauA, tauB, dist, omega)&
            & - getYGammaSubPart(tauB, tauA, dist, omega)
        tmp = 1.0_dp / dist - tmp
        tmp = tmp&
            & - getYGammaSubPart(tauA, tauB, dist, 0.0_dp)&
            & - getYGammaSubPart(tauB, tauA, dist, 0.0_dp)
        gamma = tmp
      end if
    end if

  end function getAnalyticalGammaValue


  !> Returns the subexpression for the evaluation of the off-site Y-Gamma-integral.
  pure function getYGammaSubPart(tauA, tauB, R, omega) result(yGamma)

    !> decay constant site A
    real(dp), intent(in) :: tauA

    !> decay constant site B
    real(dp), intent(in) :: tauB

    !> separation of the sites A and B
    real(dp), intent(in) :: R

    !> range-separation parameter
    real(dp), intent(in) :: omega

    !> resulting off-site Y-Gamma-integral
    real(dp) :: yGamma

    !!
    real(dp) :: prefac, tmp

    tmp = (tauA - omega)
    tmp = tmp * (tauA + omega)
    prefac = tauA * tauA / tmp
    tmp = (tauB**6 - 3.0_dp * tauA * tauA * tauB**4 + 2.0_dp * omega * omega * tauB**4) / R
    tmp = tmp * prefac * prefac / (tauA * tauA - tauB * tauB)**3
    tmp = tauA * tauB**4 * 0.5_dp * prefac / (tauB * tauB - tauA * tauA)**2 - tmp
    yGamma = tmp * exp(-tauA * R)

  end function getYGammaSubPart


  !> Derivative of analytical long-range Gamma
  function getdAnalyticalGammaDeriv(this, iSp1, iSp2, dist) result(dAnalyticalGammaDeriv)

    !> Instance
    type(TRangeSepFunc), intent(in) :: this

    !> Species index of first and second atom
    integer, intent(in) :: iSp1, iSp2

    !> Distance between atoms
    real(dp), intent(in) :: dist

    !> Resulting d gamma / d dist
    real(dp) :: dAnalyticalGammaDeriv

    real(dp) :: tauA, tauB, omega
    real(dp) :: prefac, tmp, tmp2, dTmp, dTmp2

    tauA = 3.2_dp * this%hubbu(iSp1)
    tauB = 3.2_dp * this%hubbu(iSp2)
    omega = this%omega

    if (dist < tolSameDist) then
      ! on-site case
      if (abs(tauA - tauB) < MinHubDiff ) then
        dAnalyticalGammaDeriv = 0.0_dp
      else
        call error("Error(RangeSep): R = 0, Ua != Ub")
      end if
    else
      ! off-site case, Ua == Ub
      if (abs(tauA - tauB) < MinHubDiff ) then
        tauA = 0.5_dp * (tauA + tauB)

        tmp = dist**2 * (3.0_dp*tauA**4*omega**4 - 3.0_dp * tauA**6 * omega**2 - tauA**2*omega**6)&
            & + dist * (15.0_dp*tauA**3*omega**4 - 21.0_dp*tauA**5*omega**2 - 3.0_dp*tauA*omega**6)&
            & + (15.0_dp * tauA**2 * omega**4 - 45.0_dp * tauA**4 * omega**2 - 3.0_dp * omega**6)

        dTmp = 2.0_dp*dist*(3.0_dp*tauA**4*omega**4 - 3.0_dp*tauA**6*omega**2 - tauA**2*omega**6)&
            & + (15.0_dp*tauA**3*omega**4 - 21.0_dp*tauA**5*omega**2 - 3.0_dp*tauA*omega**6)

        dtmp = (dtmp*exp(-tauA*dist) -tmp*tauA*exp(-tauA*dist))/ (48.0_dp * tauA**5)

        tmp2 = ( dist**2 * tauA**3 / 48.0_dp + 0.1875_dp * dist * tauA**2 + 0.6875_dp * tauA&
            & + 1.0_dp / dist ) * exp(-tauA * dist)

        dTmp2 = &
            & (2.0_dp*dist*tauA**3/48.0_dp + 0.1875_dp*tauA**2 -1.0_dp/dist**2) * exp(-tauA*dist)&
            & -(dist**2*tauA**3/48.0_dp + 0.1875_dp*dist*tauA**2 + 0.6875_dp*tauA +1.0_dp/dist)&
            & * tauA * exp(-tauA * dist)

        dAnalyticalGammaDeriv = -1.0_dp/dist**2 -dtmp2&
            & + (tauA**8 / (tauA**2 - omega**2)**4) * (dtmp + dtmp2 + omega*exp(-omega * dist)/dist&
            & +exp(-omega * dist) / dist**2)

      else
        ! off-site, Ua != Ub
        prefac = tauA**4 / (tauA * tauA - omega * omega )**2
        prefac = prefac * tauB**4 / (tauB * tauB - omega * omega )**2
        prefac = prefac * (-omega * exp(-omega * dist) / dist - exp(-omega * dist) / dist**2)
        dAnalyticalGammaDeriv = -1.0_dp / (dist**2) - prefac&
            & + getdYGammaSubPart(tauA, tauB, dist, omega)&
            & + getdYGammaSubPart(tauB, tauA, dist, omega)&
            & - getdYGammaSubPart(tauA, tauB, dist, 0.0_dp)&
            & - getdYGammaSubPart(tauB, tauA, dist, 0.0_dp)
      end if
    end if

  end function getdAnalyticalGammaDeriv


  !> Returns the derivative of the subexpression for the evaluation of the off-site
  !! Y-Gamma-integral. Note that tauA /= tauB.
  pure function getdYGammaSubPart(tauA, tauB, R, omega) result(dYGammaSubPart)

    !> decay constant site A
    real(dp), intent(in) :: tauA

    !> decay constant site B
    real(dp), intent(in) :: tauB

    !> separation of the sites A and B
    real(dp), intent(in) :: R

    !> range-separation parameter
    real(dp), intent(in) :: omega

    !> resulting derivative of the subexpression
    real(dp) :: dYGammaSubPart

    !! auxiliary variables
    real(dp) :: prefac, tmp, tmp2, dtmp

    tmp = tauA**2 - omega**2
    prefac = tauA * tauA / tmp
    tmp = prefac * prefac / (tauA * tauA - tauB * tauB)**3
    dtmp = tmp * (tauB**6 - 3.0_dp * tauA * tauA * tauB**4 + 2.0_dp * omega * omega * tauB**4)/R**2
    tmp = tmp * (tauB**6 - 3.0_dp * tauA * tauA * tauB**4 + 2.0_dp * omega * omega * tauB**4) / R
    tmp2 = tauA * tauB**4 * 0.5_dp * prefac / (tauB * tauB - tauA * tauA )**2 - tmp

    dYGammaSubPart = (dtmp - tmp2 * tauA) * exp(-tauA * R)

  end function getdYGammaSubPart


  !> Returns the derivative of long-range gamma for iAtom1, iAtom2.
  subroutine getGammaPrimeValue(this, grad, iAtom1, iAtom2, coords, species)

    !> Instance
    type(TRangeSepFunc), intent(in) :: this

    !> Gradient of gamma between atoms
    real(dp), intent(out) :: grad(3)

    !> First atom
    integer, intent(in) :: iAtom1

    !> Second atom
    integer, intent(in) :: iAtom2

    !> Coordinates of atoms
    real(dp), intent(in) :: coords(:,:)

    !> List of all atomic species
    integer, intent(in) :: species(:)

    !! Species index of first and second atom
    integer :: iSp1, iSp2

    !! Distance(-vector) of the two atoms
    real(dp) :: vect(3), dist

    iSp1 = species(iAtom1)
    iSp2 = species(iAtom2)

    ! analytical derivatives
    vect(:) = coords(:, iAtom1) - coords(:, iAtom2)
    dist = sqrt(sum(vect**2))
    vect(:) = vect / dist
    grad(:) = vect * getdAnalyticalGammaDeriv(this, iSp1, iSp2, dist)

  end subroutine getGammaPrimeValue


  !> Adds CAM gradients due to full-/long-range HF-contributions.
  subroutine addCamGradients(this, gradients, derivator, deltaRho, skOverCont, coords, species,&
      & orb, iSquare, ovrlapMat, iNeighbour, nNeighbourSK)

    !> Class instance
    class(TRangeSepFunc), intent(in) :: this

    !> Energy gradients
    real(dp), intent(inout) :: gradients(:,:)

    !> Density matrix difference from reference q0
    real(dp), intent(in) :: deltaRho(:,:,:)

    !> Sparse overlap part
    type(TSlakoCont), intent(in) :: skOverCont

    !> Atomic coordinates
    real(dp), intent(in) :: coords(:,:)

    !> Chemical species of atoms
    integer, intent(in) :: species(:)

    !> Orbital information for system
    type(TOrbitals), intent(in) :: orb

    !> Index for dense arrays
    integer, intent(in) :: iSquare(:)

    !> Overlap matrix
    real(dp), intent(in) :: ovrlapMat(:,:)

    !> Neighbours of atoms
    integer, intent(in) :: iNeighbour(0:,:)

    !> Number of atoms neighbouring each site where the overlap is non-zero
    integer, intent(in) :: nNeighbourSK(:)

    !> Differentiation object
    class(TNonSccDiff), intent(in) :: derivator

    ! always add the LR contribution
    call addLrGradients(this, gradients, derivator, deltaRho, skOverCont, coords, species, orb,&
        & iSquare, ovrlapMat, iNeighbour, nNeighbourSK)

    ! If xc-functional is more general, also add the full-range Hartree-Fock part.
    ! For pure LC, camAlpha would be zero anyway, but we want to save as much time as possible.
    if (this%tCam) then
      call addHartreeFockGradients(this, gradients, derivator, deltaRho, skOverCont, coords,&
          & species, orb, iSquare, ovrlapMat, iNeighbour, nNeighbourSK)
    end if

  end subroutine addCamGradients


  !> Adds gradients due to long-range HF-contribution.
  subroutine addLrGradients(this, gradients, derivator, deltaRho, skOverCont, coords, species, orb,&
      & iSquare, ovrlapMat, iNeighbour, nNeighbourSK)

    !> Instance
    type(TRangeSepFunc), intent(in) :: this

    !> Energy gradients
    real(dp), intent(inout) :: gradients(:,:)

    !> Density matrix difference from reference q0
    real(dp), intent(in) :: deltaRho(:,:,:)

    !> Sparse overlap part
    type(TSlakoCont), intent(in) :: skOverCont

    !> Atomic coordinates
    real(dp), intent(in) :: coords(:,:)

    !> Chemical species of atoms
    integer, intent(in) :: species(:)

    !> Orbital information for system
    type(TOrbitals), intent(in) :: orb

    !> Index for dense arrays
    integer, intent(in) :: iSquare(:)

    !> Overlap matrix
    real(dp), intent(in) :: ovrlapMat(:,:)

    !> Neighbours of atoms
    integer, intent(in) :: iNeighbour(0:,:)

    !> Number of atoms neighbouring each site where the overlap is non-zero
    integer, intent(in) :: nNeighbourSK(:)

    !> Differentiation object
    class(TNonSccDiff), intent(in) :: derivator

    integer :: nAtom, iAtK, iNeighK, iAtB, iNeighB, iAtC, iAtA, kpa
    real(dp) :: tmpgamma1, tmpgamma2
    real(dp) :: tmpforce(3), tmpforce_r(3), tmpforce2, tmpmultvar1
    integer :: nSpin, iSpin, mu, alpha, beta, ccc, kkk
    real(dp) :: sPrimeTmp(orb%mOrb, orb%mOrb, 3)
    real(dp) :: sPrimeTmp2(orb%mOrb, orb%mOrb, 3)
    real(dp), allocatable :: gammaPrimeTmp(:,:,:), tmpOvr(:,:), tmpRho(:,:,:), tmpderiv(:,:)

    nSpin = size(deltaRho, dim=3)
    call allocateAndInit(tmpOvr, tmpRho, gammaPrimeTmp, tmpderiv)
    nAtom = size(this%species)
    tmpderiv(:,:) = 0.0_dp
    ! sum K
    loopK: do iAtK = 1, nAtom
      ! C >= K
      loopC: do iNeighK = 0, nNeighbourSK(iAtK)
        iAtC = iNeighbour(iNeighK, iAtK)
        ! evaluate the ovr_prime
        sPrimeTmp2(:,:,:) = 0.0_dp
        sPrimeTmp(:,:,:) = 0.0_dp
        if (iAtK /= iAtC) then
          call derivator%getFirstDeriv(sPrimeTmp, skOverCont, coords, species, iAtK, iAtC, orb)
          call derivator%getFirstDeriv(sPrimeTmp2, skOverCont, coords, species, iAtC, iAtK, orb)
        end if
        loopB: do iAtB = 1, nAtom
          ! A > B
          loopA: do iNeighB = 0, nNeighbourSK(iAtB)
            iAtA = iNeighbour(iNeighB, iAtB)
            tmpgamma1 = this%lrGammaEval(iAtK, iAtB) + this%lrGammaEval(iAtC, iAtB)
            tmpgamma2 = tmpgamma1 + this%lrGammaEval(iAtK, iAtA) + this%lrGammaEval(iAtC, iAtA)
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
                tmpforce(:) = tmpforce + tmpforce2 * (gammaPrimeTmp(:, iAtK, iAtA)&
                    & + gammaPrimeTmp(:, iAtK, iAtB))
                tmpforce_r(:) = tmpforce_r + tmpforce2 * (gammaPrimeTmp(:, iAtC, iAtA)&
                    & + gammaPrimeTmp(:, iAtC, iAtB))
              else
                tmpforce(:) = tmpforce * tmpgamma1
                tmpforce_r(:) = tmpforce_r * tmpgamma1
                tmpforce(:) = tmpforce + tmpforce2 * gammaPrimeTmp(:, iAtK, iAtA)
                tmpforce_r(:) = tmpforce_r + tmpforce2 * gammaPrimeTmp(:, iAtC, iAtA)
              end if
            else
              if (iAtB /= iAtA) then
                tmpforce(:) = tmpforce + tmpforce2 * (gammaPrimeTmp(:, iAtK, iAtA)&
                    & + gammaPrimeTmp(:, iAtK, iAtB))
              else
                tmpforce(:) = tmpforce + tmpforce2 * (gammaPrimeTmp(:, iAtK, iAtA))
              end if
            end if
            tmpderiv(:, iAtK) = tmpderiv(:, iAtK) + tmpforce
            tmpderiv(:, iAtC) = tmpderiv(:, iAtC) + tmpforce_r
          end do loopA
        end do loopB
      end do loopC
    end do loopK

    if (this%tREKS) then
      gradients(:,:) = gradients - 0.5_dp * this%camBeta * tmpderiv
    else
      gradients(:,:) = gradients - 0.25_dp * this%camBeta * nSpin * tmpderiv
    end if


  contains

    !> Initialise the
    subroutine allocateAndInit(tmpOvr, tmpRho, gammaPrimeTmp, tmpderiv)

      !> Storage for the overlap
      real(dp), allocatable, intent(inout) :: tmpOvr(:,:)

      !> storage for density matrix
      real(dp), allocatable, intent(inout) :: tmpRho(:,:,:)

      !> storage for derivative of gamma interaction
      real(dp), allocatable, intent(inout) :: gammaPrimeTmp(:,:,:)

      !> workspace for the derivatives
      real(dp), allocatable, intent(inout) :: tmpderiv(:,:)

      !! Holds long-range gamma derivatives of a single interaction
      real(dp) :: tmp(3)

      !!
      integer :: iSpin, iAt1, iAt2, nAtom

      nAtom = size(this%species)

      allocate(tmpOvr(size(ovrlapMat, dim=1), size(ovrlapMat, dim=1)))
      allocate(tmpRho(size(deltaRho, dim=1), size(deltaRho, dim=1), size(deltaRho, dim=3)))
      allocate(gammaPrimeTmp(3, nAtom, nAtom))
      allocate(tmpderiv(3, size(gradients, dim=2)))
      tmpOvr = ovrlapMat
      tmpRho = deltaRho

      call symmetrizeHS(tmpOvr)

      do iSpin = 1, size(deltaRho, dim=3)
        call symmetrizeHS(tmpRho(:,:, iSpin))
      end do

      ! precompute the gamma derivatives
      gammaPrimeTmp(:,:,:) = 0.0_dp
      do iAt1 = 1, nAtom
        do iAt2 = 1, nAtom
          if (iAt1 /= iAt2) then
            call getGammaPrimeValue(this, tmp, iAt1, iAt2, coords, species)
            gammaPrimeTmp(:, iAt1, iAt2) = tmp
          end if
        end do
      end do
    end subroutine allocateAndInit

  end subroutine addLrGradients


  !> evaluate the LR-Energy contribution directly. Very slow, use addLrEnergy instead.
  function evaluateLrEnergyDirect(this, env, deltaRho, overlap, iSquare) result(energy)

    !> instance of LR
    class(TRangeSepFunc), intent(in) :: this

    !> Environment settings
    type(TEnvironment), intent(inout) :: env

    !> square density matrix
    real(dp), intent(in) :: deltaRho(:,:)

    !> square overlap matrix
    real(dp), intent(in) :: overlap(:,:)

    !> Dense matrix atom indexing
    integer, intent(in) :: iSquare(:)

    !> resulting energy
    real(dp) :: energy

    integer :: iAt1, iAt2, nAtom, mu, nu, alpha, beta
    real(dp), allocatable :: tmpOvr(:,:), tmpDRho(:,:)
    real(dp) :: tmp

    call env%globalTimer%startTimer(globalTimers%energyEval)

    nAtom = size(this%species)

    allocate(tmpOvr(size(overlap, dim=1), size(overlap, dim=1)))
    allocate(tmpDRho(size(deltaRho, dim=1), size(deltaRho, dim=1)))
    tmpOvr(:,:) = overlap
    tmpDRho(:,:) = deltaRho

    call symmetrizeHS(tmpOvr)
    call symmetrizeHS(tmpDRho)

    energy = 0.0_dp
    do iAt1 = 1, nAtom
      do iAt2 = 1, nAtom
        tmp = 0.0_dp
        do mu = iSquare(iAt1), iSquare(iAt1 + 1) - 1
          do nu = iSquare(iAt2), iSquare(iAt2 + 1) - 1
            do alpha = 1, size(tmpOvr, dim=1)
              do beta = 1, size(tmpOvr, dim=1)
                tmp = tmp + (&
                    & tmpDRho(alpha, beta) * tmpDRho(mu, nu)&
                    & + tmpDRho(mu,beta) * tmpDRho(alpha,nu)) * tmpOvr(mu,alpha) * tmpOvr(nu,beta)
              end do
            end do
          end do
        end do
        energy = energy + tmp * this%lrGammaEval(iAt1, iAt2)
      end do
    end do
    energy = -energy / 8.0_dp

    call env%globalTimer%stopTimer(globalTimers%energyEval)

  end function evaluateLrEnergyDirect


  !> Obtains the array of atomic species.
  subroutine getSpecies(this, species)

    !> Class instance
    class(TRangeSepFunc), intent(in) :: this

    !> 1D array for output, will be allocated
    integer, intent(out), allocatable :: species(:)

    species = this%species

  end subroutine getSpecies


  !> Returns long-range gamma integrals.
  subroutine getLrGamma(this, lrGamma)

    !> Class instance
    class(TRangeSepFunc), intent(in) :: this

    !> Long-range gamma integrals in AO basis
    real(dp), intent(out) :: lrGamma(:,:)

    lrGamma(:,:) = this%lrGammaEval

  end subroutine getLrGamma


  !> Calculates long-range gamma derivative integrals.
  subroutine getLrGammaDeriv(this, coords, species, lrGammaDeriv)

    !> Class instance
    class(TRangeSepFunc), intent(in) :: this

    !> Atomic coordinates
    real(dp), intent(in) :: coords(:,:)

    !> Species of all atoms including images
    integer, intent(in) :: species(:)

    !> Long-range gamma derivative integrals
    real(dp), intent(out) :: lrGammaDeriv(:,:,:)

    !! Holds long-range gamma derivatives of a single interaction
    real(dp) :: tmp(3)

    !! Number of atoms and indices of interacting atoms
    integer :: nAtom, iAt1, iAt2

    nAtom = size(lrGammaDeriv, dim=1)

    do iAt1 = 1, nAtom
      do iAt2 = 1, nAtom
        if (iAt1 /= iAt2) then
          call getGammaPrimeValue(this, tmp, iAt1, iAt2, coords, species)
          lrGammaDeriv(iAt2, iAt1, :) = tmp
        end if
      end do
    end do

  end subroutine getLrGammaDeriv


  !!
  !! Full-range Hartree-Fock related Routines.
  !!

  !> Interface routine for the full-range Hartree-Fock contribution to the Hamiltonian.
  subroutine addHartreeFockHamiltonian(this, densSqr, over, iNeighbour, nNeighbourLC, iSquare,&
      & iPair, orb, HH, overlap)

    !> Instance
    type(TRangeSepFunc), intent(inout) :: this

    ! Neighbour based screening

    !> Square (unpacked) density matrix
    real(dp), intent(in), target :: densSqr(:,:)

    !> Sparse (packed) overlap matrix.
    real(dp), intent(in) :: over(:)

    !> Neighbour indices
    integer, intent(in) :: iNeighbour(0:,:)

    !> Nr. of neighbours for each atom.
    integer, intent(in) :: nNeighbourLC(:)

    !> Position of each atom in the rows/columns of the square matrices. Shape: (nAtom)
    integer, intent(in) :: iSquare(:)

    !> Position of each (neighbour, atom) pair in the sparse matrix.
    !> Shape: (0:maxNeighbour, nAtom)
    integer, intent(in) :: iPair(0:,:)

    !> Orbital information
    type(TOrbitals), intent(in) :: orb

    !> Square (unpacked) Hamiltonian to be updated
    real(dp), intent(inout), target :: HH(:,:)

    ! Threshold based screening

    !> Square real overlap matrix
    real(dp), intent(in) :: overlap(:,:)

    select case(this%rsAlg)
    case (rangeSepTypes%matrixBased)
      call addHartreeFockHamiltonianMatrix(this, iSquare, overlap, densSqr, HH)
    end select

  end subroutine addHartreeFockHamiltonian


  !> Updates Hamiltonian with Hartree-Fock contribution using matrix-matrix multiplications.
  subroutine addHartreeFockHamiltonianMatrix(this, iSquare, overlap, densSqr, HH)

    !> Instance
    type(TRangeSepFunc), intent(inout) :: this

    !> Position of each atom in the rows/columns of the square matrices. Shape: (nAtom)
    integer, intent(in) :: iSquare(:)

    !> Square (unpacked) overlap matrix.
    real(dp), intent(in) :: overlap(:,:)

    !> Square (unpacked) density matrix
    real(dp), intent(in) :: densSqr(:,:)

    !> Square (unpacked) Hamiltonian to be updated.
    real(dp), intent(inout) :: HH(:,:)

    real(dp), allocatable :: Smat(:,:)
    real(dp), allocatable :: Dmat(:,:)
    real(dp), allocatable :: hfGammaAO(:,:)
    real(dp), allocatable :: Hhf(:,:)

    integer :: nOrb

    nOrb = size(overlap, dim=1)

    allocate(Smat(nOrb, nOrb))
    allocate(Dmat(nOrb, nOrb))
    allocate(hfGammaAO(nOrb, nOrb))
    allocate(Hhf(nOrb, nOrb))

    call allocateAndInit(this, iSquare, overlap, densSqr, HH, Smat, Dmat, hfGammaAO)
    call evaluateHamiltonian(this, Smat, Dmat, hfGammaAO, Hhf)
    HH(:,:) = HH + this%camAlpha * Hhf
    this%hfEnergy = this%hfEnergy + 0.5_dp * sum(Dmat * Hhf)

  contains

    !> Set up storage and get orbital-by-orbital gamma matrix
    subroutine allocateAndInit(this, iSquare, overlap, densSqr, HH, Smat, Dmat, hfGammaAO)

      !> Instance
      type(TRangeSepFunc), intent(in) :: this

      !> Position of each atom in the rows/columns of the square matrices. Shape: (nAtom)
      integer, intent(in) :: iSquare(:)

      !> Square (unpacked) overlap matrix.
      real(dp), intent(in) :: overlap(:,:)

      !> Square (unpacked) density matrix
      real(dp), intent(in) :: densSqr(:,:)

      !> Square (unpacked) Hamiltonian to be updated.
      real(dp), intent(inout) :: HH(:,:)

      !> Symmetrized square overlap matrix
      real(dp), intent(out) :: Smat(:,:)

      !> Symmetrized square density matrix
      real(dp), intent(out) :: Dmat(:,:)

      !> Symmetrized long-range gamma matrix
      real(dp), intent(out) :: hfGammaAO(:,:)

      integer :: nAtom, iAt, jAt

      nAtom = size(this%hfGammaEval, dim=1)

      ! Symmetrize Hamiltonian, overlap, density matrices
      call symmetrizeHS(HH)
      Smat(:,:) = overlap
      call symmetrizeHS(Smat)
      Dmat(:,:) = densSqr
      call symmetrizeHS(Dmat)

      ! Get long-range gamma variable
      hfGammaAO(:,:) = 0.0_dp
      do iAt = 1, nAtom
        do jAt = 1, nAtom
          hfGammaAO(iSquare(jAt):iSquare(jAt + 1) - 1, iSquare(iAt):iSquare(iAt + 1) - 1) =&
              & this%hfGammaEval(jAt, iAt)
        end do
      end do

    end subroutine allocateAndInit


    !> Evaluate the Hamiltonian using GEMM operations.
    subroutine evaluateHamiltonian(this, Smat, Dmat, hfGammaAO, Hhf)

      !> Instance
      type(TRangeSepFunc), intent(in) :: this

      !> Symmetrized square overlap matrix
      real(dp), intent(in) :: Smat(:,:)

      !> Symmetrized square density matrix
      real(dp), intent(in) :: Dmat(:,:)

      !> Symmetrized long-range gamma matrix
      real(dp), intent(in) :: hfGammaAO(:,:)

      !> Symmetrized full-range Hartree-Fock Hamiltonian matrix
      real(dp), intent(out) :: Hhf(:,:)

      real(dp), allocatable :: Hmat(:,:)
      real(dp), allocatable :: tmpMat(:,:)

      integer :: nOrb

      nOrb = size(Smat,dim=1)

      allocate(Hmat(nOrb,nOrb))
      allocate(tmpMat(nOrb,nOrb))

      Hhf(:,:) = 0.0_dp

      call gemm(tmpMat, Smat, Dmat)
      call gemm(Hhf, tmpMat, Smat)
      Hhf(:,:) = Hhf * hfGammaAO

      tmpMat(:,:) = tmpMat * hfGammaAO
      call gemm(Hhf, tmpMat, Smat, alpha=1.0_dp, beta=1.0_dp)

      Hmat(:,:) = Dmat * hfGammaAO
      call gemm(tmpMat, Smat, Hmat)
      call gemm(Hhf, tmpMat, Smat, alpha=1.0_dp, beta=1.0_dp)

      call gemm(tmpMat, Dmat, Smat)
      tmpMat(:,:) = tmpMat * hfGammaAO
      call gemm(Hhf, Smat, tmpMat, alpha=1.0_dp, beta=1.0_dp)

      if (this%tSpin .or. this%tREKS) then
        Hhf(:,:) = -0.25_dp * Hhf
      else
        Hhf(:,:) = -0.125_dp * Hhf
      end if

    end subroutine evaluateHamiltonian

  end subroutine addHartreeFockHamiltonianMatrix


  !> Add the full-range Hartree-Fock Energy contribution to the total energy.
  subroutine addHartreeFockEnergy(this, energy)

    !> Instance
    type(TRangeSepFunc), intent(inout) :: this

    !> Total energy
    real(dp), intent(inout) :: energy

    energy = energy + this%camAlpha * this%hfEnergy

    ! hack for spin unrestricted calculation
    this%hfEnergy = 0.0_dp

  end subroutine addHartreeFockEnergy


  !> Calculates analytical full-range Hartree-Fock gamma.
  function getAnalyticalHartreeFockGammaValue(this, iSp1, iSp2, dist) result(gamma)

    !> Instance
    type(TRangeSepFunc), intent(in) :: this

    !> First species
    integer, intent(in) :: iSp1

    !> Second species
    integer, intent(in) :: iSp2

    !> Distance between atoms
    real(dp), intent(in) :: dist

    !> Resulting gamma
    real(dp) :: gamma

    real(dp) :: tauA, tauB
    real(dp) :: prefac, tmp, tau

    tauA = 3.2_dp * this%hubbu(iSp1)
    tauB = 3.2_dp * this%hubbu(iSp2)

    if (dist < tolSameDist) then
      ! on-site case
      if (abs(tauA - tauB) < MinHubDiff) then
        tau = 0.5_dp * (tauA + tauB)
        gamma = tau * 0.3125_dp
      else
        call error("Error(RangeSep): R = 0, Ua != Ub")
      end if
    else
      ! off-site case, Ua == Ub
      if (abs(tauA - tauB) < MinHubDiff ) then
        tauA = 0.5_dp * (tauA + tauB)
        tmp = ((dist * tauA)**3 / 48.0_dp + 0.1875_dp * (dist * tauA)**2 +&
            & 0.6875_dp * (dist * tauA) + 1.0_dp) * exp(-tauA * dist) / dist
        gamma = 1.0_dp / dist - tmp
      ! off-site, Ua != Ub
      else
        gamma = 1.0_dp / dist&
            & - getYGammaSubPart(tauA, tauB, dist, 0.0_dp)&
            & - getYGammaSubPart(tauB, tauA, dist, 0.0_dp)
      end if
    end if

  end function getAnalyticalHartreeFockGammaValue


  !> Derivative of analytical full-range Hartree-Fock gamma.
  function getdAnalyticalHartreeFockGammaDeriv(this, Sp1, Sp2, dist) result(dGamma)

    !> Instance
    type(TRangeSepFunc), intent(in) :: this

    !> First species
    integer, intent(in) :: Sp1

    !> Second species
    integer, intent(in) :: Sp2

    !> Distance between atoms
    real(dp), intent(in) :: dist

    !> Resulting d gamma / d dist
    real(dp) :: dGamma

    real(dp) :: tauA, tauB
    real(dp) :: dTmp

    tauA = 3.2_dp * this%hubbu(Sp1)
    tauB = 3.2_dp * this%hubbu(Sp2)

    if (dist < tolSameDist) then
      ! on-site case
      if (abs(tauA - tauB) < MinHubDiff) then
        dGamma = 0.0_dp
      else
        call error("Error(RangeSep): R = 0, Ua != Ub")
      end if
    else
      ! off-site case, Ua == Ub
      if (abs(tauA - tauB) < MinHubDiff) then
        tauA = 0.5_dp * (tauA + tauB)

        dTmp = &
            & (2.0_dp * dist * tauA**3 / 48.0_dp + 0.1875_dp * tauA**2 - 1.0_dp / dist**2)&
            & * exp(-tauA * dist) - (dist**2 * tauA**3 / 48.0_dp + 0.1875_dp * dist * tauA**2&
            & + 0.6875_dp * tauA + 1.0_dp / dist) * tauA * exp(-tauA * dist)

        dGamma = -1.0_dp / dist**2 - dtmp

      ! off-site, Ua != Ub
      else
        dGamma = -1.0_dp / (dist**2)&
            & - getdYGammaSubPart(tauA, tauB, dist, 0.0_dp)&
            & - getdYGammaSubPart(tauB, tauA, dist, 0.0_dp)
      end if
    end if

  end function getdAnalyticalHartreeFockGammaDeriv


  !> Returns the derivative of full-range Hartree-Fock gamma for iAtom1, iAtom2.
  subroutine getHartreeFockGammaPrimeValue(this, grad, iAtom1, iAtom2, coords, species)

    !> Instance
    type(TRangeSepFunc), intent(in) :: this

    !> Gradient of gamma between atoms
    real(dp), intent(out) :: grad(3)

    !> First atom
    integer, intent(in) :: iAtom1

    !> Second atom
    integer, intent(in) :: iAtom2

    !> Coordinates of atoms
    real(dp), intent(in) :: coords(:,:)

    !> List of all atomic species
    integer, intent(in) :: species(:)

    !! Species index of first and second atom
    integer :: iSp1, iSp2

    !! Distance(-vector) of the two atoms
    real(dp) :: vect(3), dist

    iSp1 = species(iAtom1)
    iSp2 = species(iAtom2)

    ! analytical derivatives
    vect(:) = coords(:, iAtom1) - coords(:, iAtom2)
    dist = sqrt(sum(vect**2))
    vect(:) = vect / dist
    grad(:) = vect * getdAnalyticalHartreeFockGammaDeriv(this, iSp1, iSp2, dist)

  end subroutine getHartreeFockGammaPrimeValue


  !> Returns full-range HartreeFock gamma integrals.
  pure subroutine getHartreeFockGamma(this, hfGamma)

    !> Class instance
    class(TRangeSepFunc), intent(in) :: this

    !> Full-range HartreeFock gamma integrals in AO basis
    real(dp), intent(out) :: hfGamma(:,:)

    hfGamma(:,:) = this%hfGammaEval

  end subroutine getHartreeFockGamma


  !> Calculate long-range gamma derivative integrals
  subroutine getHartreeFockGammaDeriv(this, coords, species, hfGammaDeriv)

    !> Class instance
    class(TRangeSepFunc), intent(in) :: this

    !> Atomic coordinates
    real(dp), intent(in) :: coords(:,:)

    !> Species of all atoms including images
    integer, intent(in) :: species(:)

    !> Full-range HartreeFock gamma derivative integrals
    real(dp), intent(out) :: hfGammaDeriv(:,:,:)

    !! Holds long-range gamma derivatives of a single interaction
    real(dp) :: tmp(3)

    !! Number of atoms and indices of interacting atoms
    integer :: nAtom, iAt1, iAt2

    nAtom = size(hfGammaDeriv, dim=1)

    do iAt1 = 1, nAtom
      do iAt2 = 1, nAtom
        if (iAt1 /= iAt2) then
          call getHartreeFockGammaPrimeValue(this, tmp, iAt1, iAt2, coords, species)
          hfGammaDeriv(iAt2, iAt1, :) = tmp
        end if
      end do
    end do

  end subroutine getHartreeFockGammaDeriv


  !> Adds gradients due to full-range Hartree-Fock contribution.
  subroutine addHartreeFockGradients(this, gradients, derivator, deltaRho, skOverCont, coords,&
      & species, orb, iSquare, overlapMat, iNeighbour, nNeighbourSK)

    !> Class instance
    type(TRangeSepFunc), intent(in) :: this

    !> Energy gradients
    real(dp), intent(inout) :: gradients(:,:)

    !> Density matrix difference from reference q0
    real(dp), intent(in) :: deltaRho(:,:,:)

    !> Sparse overlap part
    type(TSlakoCont), intent(in) :: skOverCont

    !> Atomic coordinates
    real(dp), intent(in) :: coords(:,:)

    !> Chemical species of atoms
    integer, intent(in) :: species(:)

    !> Orbital information for system
    type(TOrbitals), intent(in) :: orb

    !> Index for dense arrays
    integer, intent(in) :: iSquare(:)

    !> Overlap matrix
    real(dp), intent(in) :: overlapMat(:,:)

    !> Neighbours of atoms
    integer, intent(in) :: iNeighbour(0:,:)

    !> Number of atoms neighbouring each site where the overlap is non-zero
    integer, intent(in) :: nNeighbourSK(:)

    !> Differentiation object
    class(TNonSccDiff), intent(in) :: derivator

    integer :: nAtom, iAtK, iNeighK, iAtB, iNeighB, iAtC, iAtA, kpa
    real(dp) :: tmpgamma1, tmpgamma2
    real(dp) :: tmpforce(3), tmpforce_r(3), tmpforce2, tmpmultvar1
    integer :: nSpin, iSpin, mu, alpha, beta, ccc, kkk
    real(dp) :: sPrimeTmp(orb%mOrb, orb%mOrb, 3)
    real(dp) :: sPrimeTmp2(orb%mOrb, orb%mOrb, 3)
    real(dp), allocatable :: gammaPrimeTmp(:,:,:), tmpOvr(:,:), tmpRho(:,:,:), tmpderiv(:,:)

    nAtom = size(this%species)
    nSpin = size(deltaRho, dim=3)
    call allocateAndInit(tmpOvr, tmpRho, gammaPrimeTmp, tmpderiv)

    ! sum K
    loopK: do iAtK = 1, nAtom
      ! C >= K
      loopC: do iNeighK = 0, nNeighbourSK(iAtK)
        iAtC = iNeighbour(iNeighK, iAtK)
        ! evaluate the ovr_prime
        sPrimeTmp2(:,:,:) = 0.0_dp
        sPrimeTmp(:,:,:) = 0.0_dp
        if (iAtK /= iAtC) then
          call derivator%getFirstDeriv(sPrimeTmp, skOverCont, coords, species, iAtK, iAtC, orb)
          call derivator%getFirstDeriv(sPrimeTmp2, skOverCont, coords, species, iAtC, iAtK, orb)
        end if
        loopB: do iAtB = 1, nAtom
          ! A > B
          loopA: do iNeighB = 0, nNeighbourSK(iAtB)
            iAtA = iNeighbour(iNeighB, iAtB)
            tmpgamma1 = this%hfGammaEval(iAtK, iAtB) + this%hfGammaEval(iAtC, iAtB)
            tmpgamma2 = tmpgamma1 + this%hfGammaEval(iAtK, iAtA) + this%hfGammaEval(iAtC, iAtA)
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
                      tmpmultvar1 = tmpmultvar1&
                          & + tmpOvr(beta, alpha) * (tmpRho(beta, kpa, iSpin)&
                          & * tmpRho(alpha, mu, iSpin) + tmpRho(alpha, kpa, iSpin)&
                          & * tmpRho(beta, mu, iSpin))
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
                tmpforce(:) = tmpforce + tmpforce2 * (gammaPrimeTmp(:, iAtK, iAtA)&
                    & + gammaPrimeTmp(:, iAtK, iAtB))
                tmpforce_r(:) = tmpforce_r + tmpforce2 * (gammaPrimeTmp(:, iAtC, iAtA)&
                    & + gammaPrimeTmp(:, iAtC, iAtB))
              else
                tmpforce(:) = tmpforce * tmpgamma1
                tmpforce_r(:) = tmpforce_r * tmpgamma1
                tmpforce(:) = tmpforce + tmpforce2 * gammaPrimeTmp(:, iAtK, iAtA)
                tmpforce_r(:) = tmpforce_r + tmpforce2 * gammaPrimeTmp(:, iAtC, iAtA)
              end if
            else
              if (iAtB /= iAtA) then
                tmpforce(:) = tmpforce + tmpforce2 * (gammaPrimeTmp(:, iAtK, iAtA)&
                    & + gammaPrimeTmp(:, iAtK, iAtB))
              else
                tmpforce(:) = tmpforce + tmpforce2 * gammaPrimeTmp(:, iAtK, iAtA)
              end if
            end if
            tmpderiv(:, iAtK) = tmpderiv(:, iAtK) + tmpforce
            tmpderiv(:, iAtC) = tmpderiv(:, iAtC) + tmpforce_r
          end do loopA
        end do loopB
      end do loopC
    end do loopK

    if (this%tREKS) then
      gradients(:,:) = gradients - 0.5_dp * this%camAlpha * tmpderiv
    else
      gradients(:,:) = gradients - 0.25_dp * this%camAlpha * nSpin * tmpderiv
    end if

  contains

    !> Initialise.
    subroutine allocateAndInit(tmpOvr, tmpRho, gammaPrimeTmp, tmpderiv)

      !> Storage for the overlap
      real(dp), intent(inout), allocatable :: tmpOvr(:,:)

      !> Storage for density matrix
      real(dp), intent(inout), allocatable :: tmpRho(:,:,:)

      !> Storage for derivative of gamma interaction
      real(dp), intent(inout), allocatable :: gammaPrimeTmp(:,:,:)

      !> Workspace for the derivatives
      real(dp), intent(inout), allocatable :: tmpderiv(:,:)

      !! Holds gamma derivative for a single iAt1-iAt2 interaction
      real(dp) :: tmp(3)

      !! Spin index (up/down)
      integer :: iSpin

      !! Indices of interacting atoms
      integer :: iAt1, iAt2

      !! Number of atoms
      integer :: nAtom

      nAtom = size(this%species)
      allocate(tmpOvr(size(overlapMat, dim=1), size(overlapMat, dim=1)))
      allocate(tmpRho(size(deltaRho, dim=1), size(deltaRho, dim=1), size(deltaRho, dim=3)))
      allocate(gammaPrimeTmp(3, nAtom, nAtom))
      allocate(tmpderiv(3, size(gradients, dim=2)))

      tmpOvr(:,:) = overlapMat
      tmpRho(:,:,:) = deltaRho
      gammaPrimeTmp(:,:,:) = 0.0_dp
      tmpderiv(:,:) = 0.0_dp

      call symmetrizeHS(tmpOvr)
      do iSpin = 1, size(deltaRho, dim=3)
        call symmetrizeHS(tmpRho(:,:, iSpin))
      end do

      ! precompute the gamma derivatives
      do iAt1 = 1, nAtom
        do iAt2 = 1, nAtom
          if (iAt1 /= iAt2) then
            call getHartreeFockGammaPrimeValue(this, tmp, iAt1, iAt2, coords, species)
            gammaPrimeTmp(:, iAt1, iAt2) = tmp
          end if
        end do
      end do

    end subroutine allocateAndInit

  end subroutine addHartreeFockGradients

end module dftbp_dftb_rangeseparated
