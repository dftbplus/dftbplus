!--------------------------------------------------------------------------------------------------!
!  DFTB+: general package for performing fast atomistic simulations                                !
!  Copyright (C) 2006 - 2020  DFTB+ developers group                                               !
!                                                                                                  !
!  See the LICENSE file for terms of usage and distribution.                                       !
!--------------------------------------------------------------------------------------------------!

#:include 'common.fypp'


!> Contains range separated related routines.
module dftbp_rangeseparated
  use dftbp_accuracy
  use dftbp_environment
  use dftbp_assert
  use dftbp_message
  use dftbp_nonscc, only : TNonSccDiff
  use dftbp_slakocont, only : TSlakoCont
  use dftbp_commontypes
  use dftbp_sorting
  use dftbp_sparse2dense, only : blockSymmetrizeHS, symmetrizeHS, hermitianSquareMatrix
  use dftbp_globalenv, only : stdOut
  use dftbp_f08math
  use dftbp_blasroutines, only : gemm
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

    !> range separation parameter
    real(dp) :: omega

  end type TRangeSepSKTag


  !> Range-Sep module structure
  type :: TRangeSepFunc
    private

    !> coordinates of the atom
    real(dp), allocatable :: coords(:,:)

    !> evaluated gamma, Atom1, Atom2 at each geometry step
    real(dp), allocatable :: lrGammaEval(:,:)

    !> range-separation parameter
    real(dp) :: omega

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

    !> total long range energy
    real(dp) :: lrEnergy

    !> spin up part of energy
    real(dp) :: lrEnergyUp

    !> spin down part of energy
    real(dp) :: lrEnergyDn

    !> Is this spin restricted (F) or unrestricted (T)
    logical :: tSpin

    !> Is this DFTB/SSR formalism
    logical :: tREKS

    !> algorithm for range separation screening
    integer :: rsAlg

    !> species of atoms
    integer, allocatable :: species(:)

    !> Nr. of neighbours
    integer, allocatable :: nNeighbours(:)

  contains

    procedure :: updateCoords
    procedure :: addLrHamiltonian
    procedure :: addLrHamiltonianMatrixCmplx
    procedure :: addLrEnergy
    procedure :: addLrGradients
    procedure :: evaluateLrEnergyDirect
    procedure :: getSpecies
    procedure :: getLrGamma
    procedure :: getLrGammaDeriv

  end type TRangeSepFunc


contains


  !> Intitialize the range-sep module
  subroutine RangeSepFunc_init(this, nAtom, species, hubbu, screen, omega, tSpin, tREKS, rsAlg)

    !> class instance
    type(TRangeSepFunc), intent(out) :: this

    !> number of atoms
    integer, intent(in) :: nAtom

    !> list of all atomic species
    integer, intent(in) :: species(:)

    !> atomic hubbards
    real(dp), intent(in) :: hubbu(:)

    !> screening threshold value
    real(dp), intent(in) :: screen

    !> range separation parameter
    real(dp), intent(in) :: omega

    !> Is this spin restricted (F) or unrestricted (T)
    logical, intent(in) :: tSpin

    !> Is this DFTB/SSR formalism
    logical, intent(in) :: tREKS

    !> lr-hamiltonian construction algorithm
    integer, intent(in) :: rsAlg

    call initAndAllocate(this, nAtom, hubbu, species, screen, omega, rsAlg, tSpin, tREKS)
    call checkRequirements(this)

  contains


    !> initialise data structures and allocate storage
    subroutine initAndAllocate(this, nAtom, hubbu, species, screen, omega, rsAlg, tSpin, tREKS)

      !> Instance
      class(TRangeSepFunc), intent(out) :: this

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
      this%rsAlg = rsAlg
      this%tSpin = tSpin
      this%tREKS = tREKS

      allocate(this%coords(3, nAtom))
      allocate(this%species(nAtom))
      allocate(this%lrGammaEval(nAtom,nAtom))
      allocate(this%hubbu(size(hubbu)))
      this%hubbu(:) = hubbu
      this%species(:) = species

    end subroutine initAndAllocate


    !> Test for option consistency
    subroutine checkRequirements(this)

      !> instance
      class(TRangeSepFunc), intent(inout) :: this

      ! Check for current restrictions
      if (this%tSpin .and. this%rsAlg == rangeSepTypes%threshold) then
        call error("Spin-unrestricted calculation for thresholding algorithm not yet implemented!")
      end if

      if (this%tREKS .and. this%rsAlg == rangeSepTypes%threshold) then
        call error("REKS calculation with thresholding algorithm not yet implemented!")
      end if

      if (.not. any([rangeSepTypes%neighbour, rangeSepTypes%threshold,&
            & rangeSepTypes%matrixBased] == this%rsAlg)) then
        call error("Unknown algorithm for screening the exchange")
      end if

    end subroutine checkRequirements

  end subroutine RangeSepFunc_init


  !> update the rangeSep module on coordinate change
  subroutine updateCoords(this, coords)

    !> class instance
    class(TRangeSepFunc), intent(inout) :: this

    !> list of atomic coordinates
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

    if (this%tScreeningInited) then
      this%hprev(:,:) = 0.0_dp
      this%dRhoPrev(:,:) = 0.0_dp
      this%lrEnergy = 0.0_dp
    end if

  end subroutine updateCoords


  !> Interface routine.
  subroutine addLrHamiltonian(this, env, densSqr, over, iNeighbour, nNeighbourLC, iSquare, iPair,&
      & orb, HH, overlap)

    !> class instance
    class(TRangeSepFunc), intent(inout) :: this

    !> Environment settings
    type(TEnvironment), intent(inout) :: env

    ! Neighbour based screening

    !> Square (unpacked) density matrix
    real(dp), dimension(:,:), target, intent(in) :: densSqr

    !> Sparse (packed) overlap matrix.
    real(dp), dimension(:), intent(in) :: over

    !> Neighbour indices.
    integer, dimension(0:,:), intent(in) :: iNeighbour

    !> Nr. of neighbours for each atom.
    integer, dimension(:), intent(in) :: nNeighbourLC

    !> Position of each atom in the rows/columns of the square matrices. Shape: (nAtom)
    integer, dimension(:), intent(in) :: iSquare

    !> iPair Position of each (neighbour, atom) pair in the sparse matrix.
    !> Shape: (0:maxNeighbour,nAtom)
    integer, dimension(0:,:), intent(in) :: iPair

    !> Orbital information.
    type(TOrbitals), intent(in) :: orb

    !> Square (unpacked) Hamiltonian to be updated.
    real(dp), dimension(:,:), intent(inout), target :: HH

    ! Threshold based screening

    !> square real overlap matrix
    real(dp), intent(in) :: overlap(:,:)

    call env%globalTimer%startTimer(globalTimers%rangeSeparatedH)
    select case(this%rsAlg)
    case (rangeSepTypes%threshold)
      call addLrHamiltonianThreshold(this, env, overlap, densSqr, iNeighbour, nNeighbourLC,&
          & iSquare, HH, orb)
    case (rangeSepTypes%neighbour)
      call addLrHamiltonianNeighbour(this, env, densSqr, over, iNeighbour, nNeighbourLC, iSquare,&
          & iPair, orb, HH)
    case (rangeSepTypes%matrixBased)
      call addLrHamiltonianMatrix(this, iSquare, overlap, densSqr, HH)
    end select
    call env%globalTimer%stopTimer(globalTimers%rangeSeparatedH)

  end subroutine addLrHamiltonian


  !> Adds the LR-exchange contribution to hamiltonian using the thresholding algorithm
  subroutine addLrHamiltonianThreshold(this, env, overlap, deltaRho, iNeighbour, nNeighbourLC,&
      & iSquare, hamiltonian, orb)

    !> class instance
    type(TRangeSepFunc), intent(inout) :: this

    !> Environment settings
    type(TEnvironment), intent(inout) :: env

    !> square real overlap matrix
    real(dp), intent(in) :: overlap(:,:)

    !> square density matrix (deltaRho in DFTB terms)
    real(dp), intent(in) :: deltaRho(:,:)

    !> Neighbour indices.
    integer, dimension(0:,:), intent(in) :: iNeighbour

    !> Nr. of neighbours for each atom.
    integer, dimension(:), intent(in) :: nNeighbourLC

    !> mapping atom_number -> number of the first basis function of the atomic block atom_number
    integer, intent(in) :: iSquare(:)

    !> current hamiltonian
    real(dp), intent(inout) :: hamiltonian(:,:)

    !> orbital information
    type(TOrbitals), intent(in) :: orb

    real(dp), allocatable :: tmpOvr(:,:), tmpDRho(:,:), testOvr(:,:), tmpDDRho(:,:), tmpDHam(:,:)
    integer, allocatable :: ovrInd(:,:)
    integer, parameter :: descLen = 3, iStart = 1, iEnd = 2, iNOrb = 3

    call allocateAndInit()
    call evaluateHamiltonian(tmpDHam)
    this%hprev(:,:) = this%hprev + tmpDHam
    hamiltonian(:,:) = hamiltonian + this%hprev
    this%lrEnergy = this%lrEnergy + evaluateEnergy(this%hprev, tmpDRho)

  contains

    !> allocate and initialise some necessary arrays
    subroutine allocateAndInit()

      integer :: matrixSize, nAtom
      real(dp) :: tmp
      integer :: iAtMu, iAtNu, iNeigh

      matrixSize = size(hamiltonian, dim = 1)
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


    !> Evaluate the update to hamiltonian due to change the in the DM
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
          if(abs(prb) < this%pScreeningThreshold) then
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
              if(abs(tstbound) < this%pScreeningThreshold) then
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


    !> Initialise the screening matrices
    subroutine checkAndInitScreening(this, matrixSize, tmpDRho)

      !> Instance
      class(TRangeSepFunc), intent(inout) :: this

      !> linear dimension of matrix
      integer, intent(in) :: matrixSize

      !> Delta rho from iteration
      real(dp), allocatable, intent(in) :: tmpDRho(:,:)

      if(.not. this%tScreeningInited) then
        allocate(this%hprev(matrixSize, matrixSize))
        allocate(this%dRhoPrev(matrixSize, matrixSize))
        this%hprev(:,:) = 0.0_dp
        this%dRhoPrev(:,:) = tmpDRho
        this%tScreeningInited = .true.
      end if

    end subroutine checkAndInitScreening

  end subroutine addLrHamiltonianThreshold


  !> Updates the Hamiltonian with the range separated contribution.
  subroutine addLrHamiltonianNeighbour(this, env, densSqr, over, iNeighbour, nNeighbourLC, iSquare,&
      & iPair, orb, HH)

    !> instance of object
    type(TRangeSepFunc), intent(inout) :: this

    !> Environment settings
    type(TEnvironment), intent(inout) :: env

    !> Square (unpacked) density matrix
    real(dp), dimension(:,:), target, intent(in) :: densSqr

    !> Sparse (packed) overlap matrix.
    real(dp), dimension(:), intent(in) :: over

    !> Neighbour indices.
    integer, dimension(0:,:), intent(in) :: iNeighbour

    !> Nr. of neighbours for each atom.
    integer, dimension(:), intent(in) :: nNeighbourLC

    !> Position of each atom in the rows/columns of the square matrices. Shape: (nAtom)
    integer, dimension(:), intent(in) :: iSquare

    !> Position of each (neighbour, atom) pair in the sparse matrix. Shape: (0:maxNeighbour, nAtom)
    integer, dimension(0:,:), intent(in) :: iPair

    !> Orbital information.
    type(TOrbitals), intent(in) :: orb

    !> Square (unpacked) Hamiltonian to be updated.
    real(dp), dimension(:,:), intent(inout), target :: HH

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
    HH(:,:) = HH + tmpHH
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

    !> class instance
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
    real(dp), allocatable :: LRgammaAO(:,:)
    complex(dp), allocatable :: gammaCmplx(:,:)
    complex(dp), allocatable :: Hlr(:,:)

    integer :: nOrb

    nOrb = size(overlap,dim=1)

    allocate(Smat(nOrb,nOrb))
    allocate(Dmat(nOrb,nOrb))
    allocate(LRgammaAO(nOrb,nOrb))
    allocate(gammaCmplx(nOrb,nOrb))
    allocate(Hlr(nOrb,nOrb))

    call allocateAndInit(this, iSquare, overlap, densSqr, HH, Smat, Dmat, LRgammaAO, gammaCmplx)

    call evaluateHamiltonian(this, Smat, Dmat, gammaCmplx, Hlr)

    HH(:,:) = HH + Hlr

    this%lrenergy = this%lrenergy + 0.5_dp * real(sum(Dmat * Hlr), dp)

  contains

    subroutine allocateAndInit(this, iSquare, overlap, densSqr, HH, Smat, Dmat, LRgammaAO,&
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
      real(dp), intent(out) :: LRgammaAO(:,:)

      !> Symmetrized long-range gamma matrix
      complex(dp), intent(out) :: gammaCmplx(:,:)

      integer :: nAtom, iAt, jAt

      nAtom = size(this%lrGammaEval,dim=1)

      !! Symmetrize Hamiltonian, overlap, density matrices
      call hermitianSquareMatrix(HH)
      Smat(:,:) = overlap
      call hermitianSquareMatrix(Smat)
      Dmat(:,:) = densSqr
      call hermitianSquareMatrix(Dmat)

      ! Get long-range gamma variable
      LRgammaAO(:,:) = 0.0_dp
      do iAt = 1, nAtom
        do jAt = 1, nAtom
          LRgammaAO(iSquare(jAt):iSquare(jAt+1)-1,iSquare(iAt):iSquare(iAt+1)-1) =&
              & this%lrGammaEval(jAt,iAt)
        end do
      end do
      gammaCmplx = LRgammaAO

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

    !> class instance
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
    real(dp), allocatable :: LrGammaAO(:,:)
    real(dp), allocatable :: Hlr(:,:)

    integer :: nOrb

    nOrb = size(overlap,dim=1)

    allocate(Smat(nOrb,nOrb))
    allocate(Dmat(nOrb,nOrb))
    allocate(LrGammaAO(nOrb,nOrb))
    allocate(Hlr(nOrb,nOrb))

    call allocateAndInit(this, iSquare, overlap, densSqr, HH, Smat, Dmat, LrGammaAO)
    call evaluateHamiltonian(this, Smat, Dmat, LrGammaAO, Hlr)
    HH(:,:) = HH + Hlr
    this%lrenergy = this%lrenergy + 0.5_dp * sum(Dmat * Hlr)

  contains

    !> Set up storage and get orbital-by-orbital gamma matrix
    subroutine allocateAndInit(this, iSquare, overlap, densSqr, HH, Smat, Dmat, LrGammaAO)

      !> class instance
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
      real(dp), intent(out) :: LrGammaAO(:,:)

      integer :: nAtom, iAt, jAt

      nAtom = size(this%lrGammaEval,dim=1)

      ! Symmetrize Hamiltonian, overlap, density matrices
      call symmetrizeHS(HH)
      Smat(:,:) = overlap
      call symmetrizeHS(Smat)
      Dmat(:,:) = densSqr
      call symmetrizeHS(Dmat)

      ! Get long-range gamma variable
      LrGammaAO(:,:) = 0.0_dp
      do iAt = 1, nAtom
        do jAt = 1, nAtom
          LrGammaAO(iSquare(jAt):iSquare(jAt+1)-1,iSquare(iAt):iSquare(iAt+1)-1) =&
              & this%lrGammaEval(jAt,iAt)
        end do
      end do

    end subroutine allocateAndInit


    !> Evaluate the hamiltonian using GEMM operations
    subroutine evaluateHamiltonian(this, Smat, Dmat, LrGammaAO, Hlr)

      !> class instance
      type(TRangeSepFunc), intent(inout) :: this

      !> Symmetrized square overlap matrix
      real(dp), intent(in) :: Smat(:,:)

      !> Symmetrized square density matrix
      real(dp), intent(in) :: Dmat(:,:)

      !> Symmetrized long-range gamma matrix
      real(dp), intent(in) :: LrGammaAO(:,:)

      !> Symmetrized long-range Hamiltonian matrix
      real(dp), intent(out) :: Hlr(:,:)

      real(dp), allocatable :: Hmat(:,:)
      real(dp), allocatable :: tmpMat(:,:)

      integer :: nOrb

      nOrb = size(Smat,dim=1)

      allocate(Hmat(nOrb,nOrb))
      allocate(tmpMat(nOrb,nOrb))

      Hlr(:,:) = 0.0_dp

      call gemm(tmpMat, Smat, Dmat)
      call gemm(Hlr, tmpMat, Smat)
      Hlr(:,:) = Hlr * LrGammaAO

      tmpMat(:,:) = tmpMat * LrGammaAO
      call gemm(Hlr, tmpMat, Smat, alpha=1.0_dp, beta=1.0_dp)

      Hmat(:,:) = Dmat * LrGammaAO
      call gemm(tmpMat, Smat, Hmat)
      call gemm(Hlr, tmpMat, Smat, alpha=1.0_dp, beta=1.0_dp)

      call gemm(tmpMat, Dmat, Smat)
      tmpMat(:,:) = tmpMat * LrGammaAO
      call gemm(Hlr, Smat, tmpMat, alpha=1.0_dp, beta=1.0_dp)

      if (this%tSpin .or. this%tREKS) then
        Hlr(:,:) = -0.25_dp * Hlr
      else
        Hlr(:,:) = -0.125_dp * Hlr
      end if

    end subroutine evaluateHamiltonian

  end subroutine addLRHamiltonianMatrix


  !> Add the LR-Energy contribution to the total energy
  subroutine addLrEnergy(this, energy)

    !> RangeSep class instance
    class(TRangeSepFunc), intent(inout) :: this

    !> total energy
    real(dp), intent(inout) :: energy

    energy = energy + this%lrEnergy
    ! hack for spin unrestricted calculation
    this%lrEnergy = 0.0_dp

  end subroutine addLrenergy


  !> location of relevant atomic block indices in a dense matrix
  pure function getDescriptor(iAt, iSquare) result(desc)

    !> relevant atom
    integer, intent(in) :: iAt

    !> indexing array for start of atom orbitals
    integer, intent(in) :: iSquare(:)

    !> resulting location ranges
    integer :: desc(3)

    desc(:) = [iSquare(iAt), iSquare(iAt + 1) - 1, iSquare(iAt + 1) - iSquare(iAt)]

  end function getDescriptor


  !> evaluate energy from triangles of the hamiltonian and density matrix
  pure function evaluateEnergy(hamiltonian, densityMat) result(egy)

    !> hamiltonian matrix
    real(dp), intent(in) :: hamiltonian(:,:)

    !> density matrix
    real(dp), intent(in) :: densityMat(:,:)

    !> resulting energy
    real(dp) :: egy

    integer :: mu

    egy = 0.0_dp
    do mu = 1, size(hamiltonian, dim=2)
      egy = egy + hamiltonian(mu, mu) * densityMat(mu, mu)&
          & + 2.0_dp * sum(hamiltonian(mu + 1 :, mu) * densityMat(mu + 1 :, mu))
    end do
    egy = 0.5_dp * egy

  end function evaluateEnergy


  !> Analytical long-range Gamma
  function getAnalyticalGammaValue(this, Sp1, Sp2, dist)

    !> RangeSepFunc instance
    class(TRangeSepFunc), intent(inout) :: this

    !> first species
    integer, intent(in) :: Sp1

    !> second species
    integer, intent(in) :: Sp2

    !> distance between atoms
    real(dp), intent(in) :: dist

    !> resulting gamma
    real(dp) :: getAnalyticalGammaValue

    real(dp) :: tauA, tauB, omega
    real(dp) :: prefac, tmp, tmp2, tau

    tauA = 3.2_dp * this%hubbu(Sp1)
    tauB = 3.2_dp * this%hubbu(Sp2)
    omega = this%omega

    if (dist < tolSameDist) then
      ! on-site case
      if (abs(tauA - tauB) < MinHubDiff) then
        tau = 0.5_dp * (tauA + tauB)
        tmp = 5.0_dp * tau**6 + 15.0_dp * tau**4 * omega**2 - 5.0_dp * tau**2 * omega**4  + omega**6
        tmp = tmp * 0.0625_dp / tau**5 - omega
        tmp = tmp * tau**8 / (tau**2 - omega**2)**4
        getAnalyticalGammaValue = tau * 0.3125_dp - tmp
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
        getAnalyticalGammaValue = 1.0_dp/dist - tmp2 - (tauA**8 / (tauA**2 - omega**2)**4 *&
            & exp(-omega * dist) / dist + tmp)
      else
        ! off-site, Ua != Ub
        prefac = tauA**4 / (tauA * tauA - omega * omega )**2
        prefac = prefac * tauB**4 / (tauB * tauB - omega * omega )**2
        prefac = prefac * exp(-omega * dist) / dist
        tmp = prefac -getYGammaSubPart(tauA,tauB,dist,omega) -getYGammaSubPart(tauB,tauA,dist,omega)
        tmp = 1.0_dp / dist - tmp
        tmp = tmp -getYGammaSubPart(tauA,tauB,dist,0.0_dp) -getYGammaSubPart(tauB,tauA,dist,0.0_dp)
        getAnalyticalGammaValue = tmp
      end if
    end if

  end function getAnalyticalGammaValue


  !> returns the subexpression for the evaluation of the off-site Y-Gamma-integral
  pure function getYGammaSubPart(tauA, tauB, R, omega)

    !> decay constant site A
    real(dp), intent(in) :: tauA

    !> decay constant site B
    real(dp), intent(in) :: tauB

    !> separation of the sites A and B
    real(dp), intent(in) :: R

    !> range-separation parameter
    real(dp), intent(in) :: omega

    real(dp) :: getYGammaSubPart
    real(dp) :: prefac, tmp

    tmp = (tauA - omega)
    tmp = tmp * (tauA + omega)
    prefac = tauA * tauA / tmp
    tmp = (tauB**6 - 3.0_dp * tauA * tauA * tauB**4 + 2.0_dp * omega * omega * tauB**4) / R
    tmp = tmp * prefac * prefac / (tauA * tauA - tauB * tauB)**3
    tmp = tauA * tauB**4 * 0.5_dp * prefac / (tauB * tauB - tauA * tauA )**2 - tmp
    getYGammaSubPart = tmp * exp(-tauA * R)

  end function getYGammaSubPart


  !> Derivative of analytical long-range Gamma
  function getdAnalyticalGammaDeriv(this, Sp1, Sp2, dist)

    !> RangeSepFunc instance
    class(TRangeSepFunc), intent(inout) :: this

    !> first species
    integer, intent(in) :: Sp1

    !> second species
    integer, intent(in) :: Sp2

    !> distance between atoms
    real(dp), intent(in) :: dist

    !> resulting d gamma / d dist
    real(dp) :: getdAnalyticalGammaDeriv

    real(dp) :: tauA, tauB, omega
    real(dp) :: prefac, tmp, tmp2, dTmp, dTmp2

    tauA = 3.2_dp * this%hubbu(Sp1)
    tauB = 3.2_dp * this%hubbu(Sp2)
    omega = this%omega

    if (dist < tolSameDist) then
      ! on-site case
      if (abs(tauA - tauB) < MinHubDiff ) then
        getdAnalyticalGammaDeriv = 0.0_dp
      else
        call error("Error(RangeSep): R = 0, Ua != Ub")
      end if
    else
      ! off-site case, Ua == Ub
      if (abs(tauA - tauB) < MinHubDiff ) then
        tauA = 0.5_dp * (tauA + tauB)

        tmp = dist**2 * (3.0_dp*tauA**4*omega**4 - 3.0_dp*tauA**6*omega**2 - tauA**2*omega**6)&
            & + dist * (15.0_dp*tauA**3*omega**4 - 21.0_dp*tauA**5*omega**2 - 3.0_dp*tauA*omega**6)&
            & + (15.0_dp * tauA**2 * omega**4 - 45.0_dp * tauA**4 * omega**2 - 3.0_dp * omega**6)

        dTmp = 2.0_dp*dist*(3.0_dp*tauA**4*omega**4 - 3.0_dp*tauA**6*omega**2 - tauA**2*omega**6)&
            & + (15.0_dp*tauA**3*omega**4 - 21.0_dp*tauA**5*omega**2 - 3.0_dp*tauA*omega**6)

        dtmp = (dtmp*exp(-tauA*dist) -tmp*tauA*exp(-tauA*dist))/ (48.0_dp * tauA**5)

        tmp2 = ( dist**2 * tauA**3 / 48.0_dp + 0.1875_dp * dist * tauA**2 + 0.6875_dp * tauA&
            & + 1.0_dp / dist ) * exp(-tauA * dist)

        dTmp2 = &
            & (2.0_dp*dist*tauA**3/48.0_dp +0.1875_dp*tauA**2 -1.0_dp/dist**2) * exp(-tauA*dist)&
            & -(dist**2*tauA**3/48.0_dp + 0.1875_dp*dist*tauA**2 + 0.6875_dp*tauA +1.0_dp/dist)&
            & * tauA * exp(-tauA * dist)

        getdAnalyticalGammaDeriv = -1.0_dp/dist**2 -dtmp2&
            & + (tauA**8 / (tauA**2 - omega**2)**4) * (dtmp + dtmp2 + omega*exp(-omega * dist)/dist&
            & +exp(-omega * dist) / dist**2)

      else
        ! off-site, Ua != Ub
        prefac = tauA**4 / (tauA * tauA - omega * omega )**2
        prefac = prefac * tauB**4 / (tauB * tauB - omega * omega )**2
        prefac = prefac * ( -omega * exp(-omega * dist) / dist - exp(-omega * dist) / dist**2)
        getdAnalyticalGammaDeriv = -1.0_dp / (dist**2) - prefac&
            & + getdYGammaSubPart(tauA,tauB,dist,omega) +getdYGammaSubPart(tauB,tauA,dist,omega)&
            & -getdYGammaSubPart(tauA,tauB,dist,0.0_dp) -getdYGammaSubPart(tauB,tauA,dist,0.0_dp)
      end if
    end if

  end function getdAnalyticalGammaDeriv


  !> returns the derivative of the subexpression for the evaluation of the off-site
  !> Y-Gamma-integral. Note that tauA /= tauB
  pure function getdYGammaSubPart(tauA, tauB, R, omega)

    !> decay constant site A
    real(dp), intent(in) :: tauA

    !> decay constant site B
    real(dp), intent(in) :: tauB

    !> separation of the sites A and B
    real(dp), intent(in) :: R

    !> range-separation parameter
    real(dp), intent(in) :: omega

    real(dp) :: getdYGammaSubPart
    real(dp) :: prefac, tmp, tmp2, dtmp

    tmp = tauA**2 - omega**2
    prefac = tauA * tauA / tmp
    tmp = prefac * prefac / (tauA * tauA - tauB * tauB)**3
    dtmp = tmp * (tauB**6 - 3.0_dp * tauA * tauA * tauB**4 + 2.0_dp * omega * omega * tauB**4)/R**2
    tmp = tmp * (tauB**6 - 3.0_dp * tauA * tauA * tauB**4 + 2.0_dp * omega * omega * tauB**4) / R
    tmp2 = tauA * tauB**4 * 0.5_dp * prefac / (tauB * tauB - tauA * tauA )**2 - tmp

    getdYGammaSubPart = (dtmp -tmp2 * tauA) * exp(-tauA * R)

  end function getdYGammaSubPart


  !> Returns the derivative of lr-gamma for iAtom1, iAtom2
  subroutine getGammaPrimeValue(this, grad, iAtom1, iAtom2, coords, species)

    !> class instance
    class(TRangeSepFunc), intent(inout) :: this

    !> gradient of gamma between atoms
    real(dp), intent(out) :: grad(3)

    !> first atom
    integer, intent(in) :: iAtom1

    !> second atom
    integer, intent(in) :: iAtom2

    !> coordinates of atoms
    real(dp), intent(in) :: coords(:,:)

    !> list of all atomic species
    integer, intent(in) :: species(:)

    integer :: sp1, sp2
    real(dp) :: vect(3), dist

    sp1 = species(iAtom1)
    sp2 = species(iAtom2)

    ! analytical derivatives
    vect(:) = coords(:,iAtom1) - coords(:,iAtom2)
    dist = sqrt(sum(vect(:)**2))
    vect(:) = vect(:) / dist
    grad(:) = vect(:) * getdAnalyticalGammaDeriv(this, sp1, sp2, dist)

  end subroutine getGammaPrimeValue


  !> Adds gradients due to long-range HF-contribution
  subroutine addLrGradients(this, gradients, derivator, deltaRho, skOverCont, coords, species, orb,&
      & iSquare, ovrlapMat, iNeighbour, nNeighbourSK)

    !> class instance
    class(TRangeSepFunc), intent(inout) :: this

    !> energy gradients
    real(dp), intent(inout) :: gradients(:,:)

    !> density matrix difference from reference q0
    real(dp), intent(in) :: deltaRho(:,:,:)

    !> sparse overlap part
    type(TSlakoCont), intent(in) :: skOverCont

    !> atomic coordinates
    real(dp), intent(in) :: coords(:,:)

    !> chemical species of atoms
    integer, intent(in) :: species(:)

    !> orbital information for system
    type(TOrbitals), intent(in) :: orb

    !> index for dense arrays
    integer, intent(in) :: iSquare(:)

    !> overlap matrix
    real(dp), intent(in) :: ovrlapMat(:,:)

    !> neighbours of atoms
    integer, intent(in) :: iNeighbour(0:,:)

    !> number of atoms neighbouring each site where the overlap is non-zero
    integer, intent(in) :: nNeighbourSK(:)

    !> differentiation object
    class(TNonSccDiff), intent(in) :: derivator

    integer :: nAtom, iAtK, iNeighK, iAtB, iNeighB, iAtC, iAtA, kpa
    real(dp) :: tmpgamma1, tmpgamma2
    real(dp) :: tmpforce(3), tmpforce_r(3), tmpforce2, tmpmultvar1
    integer :: nSpin, iSpin, mu, alpha, beta, ccc, kkk
    real(dp) :: sPrimeTmp(orb%mOrb,orb%mOrb,3)
    real(dp) :: sPrimeTmp2(orb%mOrb,orb%mOrb,3)
    real(dp), allocatable :: gammaPrimeTmp(:,:,:), tmpOvr(:,:), tmpRho(:,:,:), tmpderiv(:,:)

    nSpin = size(deltaRho,dim=3)
    call allocateAndInit(tmpOvr, tmpRho, gammaPrimeTmp, tmpderiv)
    nAtom = size(this%species)
    tmpderiv = 0.0_dp
    ! sum K
    loopK: do iAtK = 1, nAtom
      ! C >= K
      loopC: do iNeighK = 0, nNeighbourSK(iAtK)
        iAtC = iNeighbour(iNeighK, iAtK)
        ! evaluate the ovr_prime
        sPrimeTmp2 = 0.0_dp
        sPrimeTmp = 0.0_dp
        if ( iAtK /= iAtC ) then
          call derivator%getFirstDeriv(sPrimeTmp, skOverCont, coords, species, iAtK, iAtC, orb)
          call derivator%getFirstDeriv(sPrimeTmp2, skOverCont, coords, species, iAtC, iAtK, orb)
        end if
        loopB: do iAtB = 1, nAtom
          ! A > B
          loopA: do iNeighB = 0, nNeighbourSK(iAtB)
            iAtA = iNeighbour(iNeighB, iAtB)
            tmpgamma1 = this%lrGammaEval(iAtK,iAtB) + this%lrGammaEval(iAtC,iAtB)
            tmpgamma2 = tmpgamma1 + this%lrGammaEval(iAtK,iAtA) + this%lrGammaEval(iAtC,iAtA)
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
                      tmpmultvar1 = tmpmultvar1 + tmpOvr(beta, alpha) * (tmpRho(beta,kpa,iSpin) &
                       & * tmpRho(alpha,mu,iSpin) + tmpRho(alpha,kpa,iSpin) * tmpRho(beta,mu,iSpin))
                    end do
                  end do
                end do
                tmpforce(:) = tmpforce(:) + tmpmultvar1 * (sPrimeTmp(ccc,kkk,:))
                tmpforce_r(:) = tmpforce_r(:) + tmpmultvar1 * (sPrimeTmp2(kkk,ccc,:))
                tmpforce2 = tmpforce2 + tmpmultvar1 * tmpOvr(kpa,mu)
              end do
            end do

            ! C /= K
            if( iAtK /= iAtC ) then
              if( iAtB /= iAtA) then
                tmpforce(:) = tmpforce(:) * tmpgamma2
                tmpforce_r(:) = tmpforce_r(:) * tmpgamma2
                tmpforce(:) = tmpforce(:) + tmpforce2 * (gammaPrimeTmp(:,iAtK,iAtA) &
                    & + gammaPrimeTmp(:,iAtK,iAtB))
                tmpforce_r(:) = tmpforce_r(:) + tmpforce2 * (gammaPrimeTmp(:,iAtC,iAtA) &
                    & + gammaPrimeTmp(:,iAtC,iAtB))
              else
                tmpforce(:) = tmpforce(:) * tmpgamma1
                tmpforce_r(:) = tmpforce_r(:) * tmpgamma1
                tmpforce(:) = tmpforce(:) + tmpforce2 * (gammaPrimeTmp(:,iAtK,iAtA))
                tmpforce_r(:) = tmpforce_r(:) + tmpforce2 * (gammaPrimeTmp(:,iAtC,iAtA))
              end if
            else
              if( iAtB /= iAtA) then
                tmpforce(:) = tmpforce(:) + tmpforce2 * (gammaPrimeTmp(:,iAtK,iAtA) &
                    & + gammaPrimeTmp(:,iAtK,iAtB))
              else
                tmpforce(:) = tmpforce(:) + tmpforce2 * (gammaPrimeTmp(:,iAtK,iAtA))
              end if
            end if
            tmpderiv(:,iAtK) = tmpderiv(:,iAtK) + tmpforce(:)
            tmpderiv(:,iAtC) = tmpderiv(:,iAtC) + tmpforce_r(:)
          end do loopA
        end do loopB
      end do loopC
    end do loopK

    if (this%tREKS) then
      gradients(:,:) = gradients - 0.5_dp * tmpderiv
    else
      gradients(:,:) = gradients - 0.25_dp * nSpin * tmpderiv
    end if

    deallocate(tmpOvr, tmpRho, gammaPrimeTmp, tmpderiv)

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

      real(dp) :: tmp(3)
      integer :: iSpin, iAt1, iAt2, nAtom

      nAtom = size(this%species)
      allocate(tmpOvr(size(ovrlapMat, dim = 1), size(ovrlapMat, dim = 1)))
      allocate(tmpRho(size(deltaRho, dim = 1), size(deltaRho, dim = 1), size(deltaRho, dim = 3)))
      allocate(gammaPrimeTmp(3, nAtom, nAtom))
      allocate(tmpderiv(3, size(gradients, dim = 2)))
      tmpOvr = ovrlapMat
      tmpRho = deltaRho
      call symmetrizeHS(tmpOvr)
      do iSpin = 1, size(deltaRho, dim = 3)
        call symmetrizeHS(tmpRho(:,:,iSpin))
      enddo
      ! precompute the gamma derivatives
      gammaPrimeTmp = 0.0_dp
      do iAt1 = 1, nAtom
        do iAt2 = 1, nAtom
          if(iAt1 /= iAt2) then
            call getGammaPrimeValue(this, tmp, iAt1, iAt2, coords, species)
            gammaPrimeTmp(:,iAt1, iAt2) = tmp(:)
          end if
        end do
      end do
    end subroutine allocateAndInit

  end subroutine addLrGradients


  !> evaluate the LR-Energy contribution directly. Very slow, use addLrEnergy instead.
  function evaluateLrEnergyDirect(this, env, deltaRho, ovrlap, iSquare) result (energy)

    !> instance of LR
    class(TRangeSepFunc), intent(inout) :: this

    !> Environment settings
    type(TEnvironment), intent(inout) :: env

    !> square density matrix
    real(dp), intent(in) :: deltaRho(:,:)

    !> square overlap matrix
    real(dp), intent(in) :: ovrlap(:,:)

    !> Dense matrix atom indexing
    integer, intent(in) :: iSquare(:)

    !> resulting energy
    real(dp) :: energy

    integer :: iAt1, iAt2, nAtom, mu, nu, alpha, beta
    real(dp), allocatable :: tmpOvr(:,:), tmpDRho(:,:)
    real(dp) :: tmp

    call env%globalTimer%startTimer(globalTimers%energyEval)
    nAtom = size(this%species)
    allocate(tmpOvr(size(ovrlap, dim = 1), size(ovrlap, dim = 1)))
    allocate(tmpDRho(size(deltaRho, dim = 1), size(deltaRho, dim = 1)))
    tmpOvr = ovrlap
    tmpDRho = deltaRho
    call symmetrizeHS(tmpOvr)
    call symmetrizeHS(tmpDRho)

    energy = 0.0_dp
    do iAt1 = 1, nAtom
      do iAt2 = 1, nAtom
        tmp = 0.0_dp
        do mu = iSquare(iAt1), iSquare(iAt1 + 1) - 1
          do nu = iSquare(iAt2), iSquare(iAt2 + 1) - 1
            do alpha = 1, size(tmpOvr, dim = 1)
              do beta = 1, size(tmpOvr, dim = 1)
                tmp = tmp + (&
                    & tmpDRho(alpha,beta) * tmpDRho(mu,nu)&
                    & +tmpDRho(mu,beta) * tmpDRho(alpha,nu)) * tmpOvr(mu,alpha) * tmpOvr(nu,beta)
              end do
            end do
          end do
        end do
        energy = energy + tmp * this%lrGammaEval(iAt1,iAt2)
      end do
    end do
    energy = -energy / 8.0_dp

    call env%globalTimer%stopTimer(globalTimers%energyEval)

  end function evaluateLrEnergyDirect


  !> obtain the array of atomic species
  subroutine getSpecies(this, targetArray)

    !> 1D array for output, will be allocated
    integer, allocatable, intent(out) :: targetArray(:)

    !> instance
    class(TRangeSepFunc), intent(in) :: this

    !> dimension of the species array
    integer :: dim1

    dim1 = size(this%species)
    allocate(targetArray(dim1))

    targetArray(:) = this%species

  end subroutine getSpecies


  !> Get long-range gamma integrals
  subroutine getLrGamma(this, LrGamma)

    !> class instance
    class(TRangeSepFunc), intent(inout) :: this

    !> long-range gamma integrals in AO basis
    real(dp), intent(out) :: LrGamma(:,:)

    LrGamma(:,:) = this%lrGammaEval

  end subroutine getLrGamma


  !> Calculate long-range gamma derivative integrals
  subroutine getLrGammaDeriv(this, coords, species, LrGammaDeriv)

    !> class instance
    class(TRangeSepFunc), intent(inout) :: this

    !> atomic coordinates
    real(dp), intent(in) :: coords(:,:)

    !> Species of all atoms including images
    integer, intent(in) :: species(:)

    !> long-range gamma derivative integrals
    real(dp), intent(out) :: LrGammaDeriv(:,:,:)

    real(dp) :: tmp(3)
    integer :: nAtom, iAt1, iAt2

    nAtom = size(LrGammaDeriv,dim=1)

    do iAt1 = 1, nAtom
      do iAt2 = 1, nAtom
        if (iAt1 /= iAt2) then
          call getGammaPrimeValue(this, tmp, iAt1, iAt2, coords, species)
          LrGammaDeriv(iAt2,iAt1,:) = tmp
        end if
      end do
    end do

  end subroutine getLrGammaDeriv

end module dftbp_rangeseparated
