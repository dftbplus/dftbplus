!--------------------------------------------------------------------------------------------------!
!  DFTB+: general package for performing fast atomistic simulations                                !
!  Copyright (C) 2018  DFTB+ developers group                                                      !
!                                                                                                  !
!  See the LICENSE file for terms of usage and distribution.                                       !
!--------------------------------------------------------------------------------------------------!

#:include 'common.fypp'


!> Contains range separated related routines.
module rangeseparated
  use accuracy
  use environment
  use assert
  use message
  use nonscc, only : NonSccDiff
  use SlakoCont, only : OSlakoCont
  use CommonTypes
  use sorting
  use sparse2dense, only : blockSymmetrizeHS
  use globalenv, only : stdOut
  implicit none
  private

  public :: RangeSepFunc, TRangeSepSKTag, RangeSep_init

  !> Slater-Koster file RangeSep tag structure
  type :: TRangeSepSKTag

    !> range separation parameter
    real(dp) :: omega

  end type TRangeSepSKTag

  !> Range-Sep module structure
  type :: RangeSepFunc

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
    real(dp), allocatable :: dRhoprev(:,:)

    !> Is screening initialised
    logical :: tScreeningInited

    !> threshold for screening by value
    real(dp) :: pScreeningThreshold

    ! lr-energy

    !> total long range energy
    real(dp) :: lrenergy

    !> spin up part of energy
    real(dp) :: lrenergyUp

    !> spin down part of energy
    real(dp) :: lrenergyDn

    !> is the system spin unrestricted or restricted
    logical :: tSpin

    !> algorithm for range separation screening
    character(lc) :: RSAlg

    !> species of atoms
    integer, allocatable :: species(:)

  contains

    procedure :: updateCoords
    procedure :: addLRHamiltonian
    procedure :: addLREnergy
    procedure :: addLRGradients
    procedure :: evaluateLREnergyDirect

  end type RangeSepFunc


contains


  !> Intitialize the range-sep module
  subroutine RangeSep_init(self, nAtom, species, speciesNames, hubbu, screen, omega, tSpin, RSAlg)

    !> class instance
    type(RangeSepFunc), intent(inout) :: self

    !> number of atoms
    integer, intent(in) :: nAtom

    !> list of all atomic species
    integer, intent(in) :: species(:)

    !> list of all atomic species names
    character(mc), intent(in) :: speciesNames(:)

    !> atomic hubbards
    real(dp), intent(in) :: hubbu(:)

    !> screening threshold value
    real(dp), intent(in) :: screen

    !> range separation parameter
    real(dp), intent(in) :: omega

    !> spin unrestricted?
    logical, intent(in) :: tSpin

    !> lr-hamiltonian construction algorithm
    character(lc), intent(in) :: RSAlg

    call initAndAllocate(self, nAtom, hubbu, species, screen, omega, RSAlg, tSpin)
    call printModuleInfoAndCheckReqs(self)

  contains


    !> initialise data structures and allocate storage
    subroutine initAndAllocate(self, nAtom, hubbu, species, screen, omega, RSAlg, tSpin)

      !> Instance
      class(RangeSepFunc), intent(inout) :: self

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
      character(lc), intent(in) :: RSAlg

      !> Is this spin restricted (F) or unrestricted (T)
      logical, intent(in) :: tSpin

      self%tScreeningInited = .false.
      self%pScreeningThreshold = screen
      self%omega = omega
      self%lrenergy = 0.0_dp
      self%RSAlg = RSAlg
      self%tSpin = tSpin

      allocate(self%coords(3, nAtom))
      allocate(self%species(nAtom))
      allocate(self%lrGammaEval(nAtom,nAtom))
      allocate(self%hubbu(size(hubbu(:))))
      self%hubbu = hubbu
      self%species(:) = species

    end subroutine initAndAllocate


    !> test for option consistency and print some information
    subroutine printModuleInfoAndCheckReqs(self)

      !> instance
      class(RangeSepFunc), intent(inout) :: self

      write(StdOut,'(a)') "================================"
      write(StdOut,'(a)') "Range-separated Hybrids in DFTB "
      write(StdOut,'(a)') ""
      write(StdOut,'(a)') "================================"
      write(StdOut,'(a)') "=> Initializing RangeSep module"

      ! Check for current restrictions
      if (self%tSpin .and. self%RSAlg == "tr") then
        call error("Spin-unrestricted calculation for thresholding algorithm not yet implemented!")
      end if

      ! summarise module settings
      write(StdOut,'(a,F10.3)') "  -> range-separation parameter [1/a0]:", self%omega

      if (self%tSpin) then
        write(StdOut,'(a)') "  -> spin-unrestricted calculation"
      else
        write(StdOut,'(a)') "  -> spin-restricted calculation"
      end if

      select case (self%RSAlg)
      case ("nb")
        write(StdOut,'(a)') "  -> using the neighbour list-based algorithm"
      case ("tr")
        write(StdOut,'(a)') "  -> using the thresholding algorithm"
        write(StdOut,'(a,E17.8)') "     -> Screening Threshold:", self%pScreeningThreshold
      case default
        call error("Invalid algorithm for screening exchange")
      end select

    end subroutine printModuleInfoAndCheckReqs

  end subroutine RangeSep_init


  !> update the rangeSep module on coordinate change
  subroutine updateCoords(self, coords)

    !> class instance
    class(RangeSepFunc), intent(inout) :: self

    !> list of atomic coordinates
    real(dp), intent(in) :: coords(:,:)

    integer :: nAtom, iAtom1,iAtom2,ii,iSp1,iSp2
    real(dp) :: dist

    @:ASSERT(all(shape(coords) == shape(self%coords)))
    self%coords(:,:) = coords
    nAtom = size(self%species)
    dist = 0.0_dp
    do iAtom1 = 1, nAtom
      do iAtom2 = 1, iAtom1
        iSp1 = self%species(iAtom1)
        iSp2 = self%species(iAtom2)
        dist = 0.0_dp
        do ii = 1, 3
          dist = dist + (self%coords(ii, iAtom1) - self%coords(ii, iAtom2))**2
        end do
        dist = sqrt(dist)
        self%lrGammaEval(iAtom1, iAtom2) = getAnalyticalGammaValue(self, iSp1, iSp2, dist)
        self%lrGammaEval(iAtom2, iAtom1) = self%lrGammaEval(iAtom1, iAtom2)
      end do
    end do
    ! reinitialise the screening
    if( self%tScreeningInited) then
      self%hprev = 0.0_dp
      self%dRhoprev = 0.0_dp
      !
      self%lrenergy = 0.0_dp
    end if

  end subroutine updateCoords


  !> Adds the LR-exchange contribution to hamiltonian using the thresholding algorithm
  subroutine addLRHamiltonian_tr(self, env, overlap, deltaRho, iSquare, hamiltonian, orb)

    !> class instance
    type(RangeSepFunc), intent(inout) :: self

    !> Environment settings
    type(TEnvironment), intent(inout) :: env

    !> square real overlap matrix
    real(dp), intent(in) :: overlap(:,:)

    !> square density matrix (deltaRho in DFTB terms)
    real(dp), intent(in) :: deltaRho(:,:)

    !> mapping atom_number -> number of the first basis function of the atomic block atom_number
    integer, intent(in) :: iSquare(:)

    !> current hamiltonian
    real(dp), intent(inout) :: hamiltonian(:,:)

    !> orbital information
    type(TOrbitals), intent(in) :: orb

    real(dp), allocatable :: tmpovr(:,:), tmpDRho(:,:), testovr(:,:), tmpDDRho(:,:), tmpDham(:,:)
    integer, allocatable :: ovrind(:,:)
    integer, parameter :: DESC_LEN = 3, ISTART = 1, IEND = 2, INORB = 3

    call allocateAndInit(tmpovr, tmpDham, tmpDRho, tmpDDRho, testovr, ovrind)
    call evaluateHamiltonian(tmpDHam)
    self%hprev = self%hprev + tmpDham
    hamiltonian = hamiltonian + self%hprev
    self%lrenergy = self%lrenergy + evaluateEnergy(self%hprev, tmpDRho)

  contains

    !> allocate and initialise some necessary arrays
    subroutine allocateAndInit(tmpovr, tmpDham, tmpDRho, tmpDDRho, testovr, ovrind)

      !> overlap matrix
      real(dp), allocatable, intent(inout) :: tmpovr(:,:)

      !> Update on hamiltonian from DM changes
      real(dp), allocatable, intent(inout) :: tmpDham(:,:)

      !> Density matrix minus reference density matrix
      real(dp), allocatable, intent(inout) :: tmpDRho(:,:)

      !> Change in tmpDRho between iteration, used for update
      real(dp), allocatable, intent(inout) :: tmpDDRho(:,:)

      !> matrix of test values for overlap (based on maximum overlap elements between atoms)
      real(dp), allocatable, intent(inout) :: testovr(:,:)

      !> sorted index array for maximal overlap elements between atom blocks
      integer, allocatable, intent(inout) :: ovrind(:,:)

      integer :: matrixSize, nAtom
      real(dp) :: tmp
      integer :: iAtMu, iAtNu

      matrixSize = size(hamiltonian, dim = 1)
      nAtom = size(self%species)
      allocate(tmpovr(matrixSize, matrixSize))
      allocate(tmpDham(matrixSize, matrixSize))
      allocate(tmpDRho(matrixSize, matrixSize))
      allocate(tmpDDRho(matrixSize, matrixSize))
      allocate(testovr(nAtom,nAtom))
      allocate(ovrind(nAtom,nAtom))
      tmpovr = overlap
      call blockSymmetrizeHS(tmpovr, iSquare)
      tmpDRho = deltaRho
      call symmetrizeSquareMatrix(tmpDRho)
      tmpDham = 0.0_dp
      call checkAndInitScreening(self, matrixSize, tmpDRho)
      tmpDDRho = tmpDRho - self%dRhoprev
      self%dRhoprev = tmpDRho
      do iAtMu = 1, nAtom
        do iAtNu = 1, nAtom
          tmp = maxval( abs( tmpovr(iSquare(iAtMu):iSquare(iAtMu + 1) - 1,&
              & iSquare(iAtNu):iSquare(iAtNu + 1) - 1) ) )
          testovr(iAtMu,iAtNu) = tmp
        end do
      end do
      do iAtMu = 1, nAtom
        call index_heap_sort(ovrind(iAtMu,:),testovr(iAtMu,:))
      end do

    end subroutine allocateAndInit


    !> Evaluate the update to hamiltonian due to change the in the DM
    pure subroutine evaluateHamiltonian(tmpDHam)

      !> Update for the old hamiltonian on exit
      real(dp), intent(out) :: tmpDHam(:,:)

      integer :: nAtom
      real(dp) :: pbound, prb
      real(dp) :: tmpvec1(orb%mOrb), tmpvec2(orb%mOrb), tmpvec3(orb%mOrb)
      real(dp) :: tmp, tstbound, gammabatch, gammabatchtmp
      integer :: iAtMu, iAtNu, iAt1, iAt2, iSp1, iSp2, nOrb1, nOrb2
      integer :: kk, ll, jj, ii, mu, nu
      integer, dimension(DESC_LEN) :: descA, descB, descM, descN

      nAtom = size(self%species)

      pbound = maxval(abs(tmpDDRho))
      tmpDham = 0.0_dp
      loopMu: do iAtMu = 1, nAtom
        descM = getDescriptor(iAtMu, iSquare)
        loopKK: do kk = 1, nAtom
          iAt1 = ovrind(iAtMu, nAtom + 1 - kk)
          descA = getDescriptor(iAt1, iSquare)
          iSp1 = self%species(iAt1)
          nOrb1 = orb%nOrbSpecies(iSp1)
          prb = pbound * testovr(iAt1, iAtMu)
          if(abs(prb) >= self%pScreeningThreshold) then
            loopNu: do iAtNu = 1, iAtMu
              descN = getDescriptor(iAtNu, iSquare)
              gammabatchtmp = self%lrGammaEval(iAtMu, iAtNu) + self%lrGammaEval(iAt1, iAtNu)
              loopLL: do ll = 1, nAtom
                iAt2 = ovrind(iAtNu, nAtom + 1 - ll)
                iSp2 = self%species(iAt2)
                nOrb2 = orb%nOrbSpecies(iSp2)
                ! screening condition
                tstbound = prb * testovr(iAt2, iAtNu)
                if(abs(tstbound) >= self%pScreeningThreshold) then
                  descB = getDescriptor(iAt2, iSquare)
                  gammabatch = (self%lrGammaEval(iAtMu, iAt2) + self%lrGammaEval(iAt1, iAt2)&
                      & + gammabatchtmp)
                  gammabatch = -0.125_dp * gammabatch
                  ! calculate the Q_AB
                  do nu = descN(ISTART), descN(IEND)
                    jj = 0
                    tmpvec2(1:nOrb2) = tmpovr(descB(ISTART):descB(IEND), nu)
                    do ii = descA(ISTART), descA(IEND)
                      jj = jj + 1
                      tmpvec1(jj) = sum(tmpvec2(1:nOrb2) * tmpDDRho(ii, descB(ISTART):descB(IEND)))
                    end do
                    tmp = 0.0_dp
                    do mu = descM(ISTART), descM(IEND)
                      tmp = sum(tmpovr(descA(ISTART):descA(IEND), mu) * tmpvec1(1:nOrb1))
                      tmpDham(mu, nu) = tmpDham(mu, nu) + gammabatch * tmp
                    end do
                  end do
                else
                  exit
                end if
              end do loopLL
            end do loopNu
          else
            exit
          end if
        end do loopKK
      end do loopMu

    end subroutine evaluateHamiltonian


    !> Initialise the screening matrices
    subroutine checkAndInitScreening(self, matrixSize, tmpDRho)

      !> Instance
      class(RangeSepFunc), intent(inout) :: self

      !> linear dimension of matrix
      integer, intent(in) :: matrixSize

      !> Delta rho from iteration
      real(dp), allocatable, intent(in) :: tmpDRho(:,:)

      if(.not. self%tScreeningInited) then
        allocate(self%hprev(matrixSize, matrixSize))
        allocate(self%dRhoprev(matrixSize, matrixSize))
        self%hprev = 0.0_dp
        self%dRhoprev = tmpDRho
        self%tScreeningInited = .true.
      end if

    end subroutine checkAndInitScreening

  end subroutine addLRHamiltonian_tr


  !> Updates the Hamiltonian with the range separated contribution.
  subroutine addLRHamiltonian_nb(self, env, densSqr, over, iNeighbour, nNeighbourLC, iSquare,&
      & iPair, orb, HH)

    !> instance of object
    type(RangeSepFunc), intent(inout) :: self

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

    integer, parameter :: DESC_LEN = 3, ISTART = 1, IEND = 2, INORB = 3
    real(dp), dimension(orb%mOrb**2), target :: Sma, Sam, Snb, Sbn
    real(dp), dimension(orb%mOrb**2), target :: Pab, Pmb, Pan, Pmn
    real(dp), dimension(:,:), pointer :: pSma, pSam, pSnb, pSbn, pHH
    real(dp), dimension(:,:), pointer :: pPab, pPmb, pPan, pPmn
    real(dp) :: gamma1, gamma2, gammaTot
    integer :: nAtom, ii, jj
    integer :: iAtM, iAtN, iAtA, iAtB, iNeighN, iNeighA
    integer, dimension(DESC_LEN) :: descA, descB, descM, descN
    real(dp), dimension(:,:), allocatable, target :: tmpDRho
    real(dp), dimension(:,:), allocatable, target :: tmpHH
    real(dp) :: tmp
    integer :: mu, nu

    call allocateAndInit(tmpHH, tmpDRho)
    call evaluateHamiltonian()
    HH = HH + tmpHH
    self%lrenergy = self%lrenergy + evaluateEnergy(tmpHH, tmpDRho)

  contains

    !> Allocate storage for mapping 1D<->2D array sections
    subroutine allocateAndInit(tmpHH, tmpDRho)

      !> density matrix case
      real(dp), dimension(:,:), allocatable, target, intent(inout) :: tmpDRho

      !> hamiltonian matrix case
      real(dp), dimension(:,:), allocatable, target, intent(inout) :: tmpHH

      allocate(tmpHH(size(HH, dim = 1), size(HH, dim = 2)))
      tmpHH = 0.0_dp
      allocate(tmpDRho(size(densSqr, dim = 1), size(densSqr, dim = 1)))
      tmpDRho = densSqr
      call symmetrizeSquareMatrix(tmpDRho)

    end subroutine allocateAndInit


    !> actually evaluate the neighbour based cut-off hamiltonian
    subroutine evaluateHamiltonian()

      nAtom = size(self%species)

      loopN: do iAtN = 1, nAtom
        descN = getDescriptor(iAtN, iSquare)
        loopB: do iNeighN = 0, nNeighbourLC(iAtN)
          iAtB = iNeighbour(iNeighN, iAtN)
          descB = getDescriptor(iAtB, iSquare)
          call copyOverlapBlock(iAtN, iNeighN, descN(INORB), descB(INORB), Sbn, pSbn)
          call transposeBlock(pSbn, Snb, pSnb)
          loopA: do iAtA = 1, nAtom
            descA = getDescriptor(iAtA, iSquare)
            call copyDensityBlock(descA, descB, Pab, pPab)
            call copyDensityBlock(descA, descN, Pan, pPan)
            gamma1 = self%lrGammaEval(iAtA, iAtN) + self%lrGammaEval(iAtA, iAtB)
            loopM: do iNeighA = 0, nNeighbourLC(iAtA)
              iAtM = iNeighbour(iNeighA, iAtA)
              descM = getDescriptor(iAtM, iSquare)
              call copyOverlapBlock(iAtA, iNeighA, descA(INORB), descM(INORB), Sma, pSma)
              call transposeBlock(pSma, Sam, pSam)
              gamma2 = self%lrGammaEval(iAtM, iAtN) + self%lrGammaEval(iAtM, iAtB)
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
    pure subroutine copyOverlapBlock(iAt, iNeigh, nOrbAt, nOrbNeigh, localBlock, pLocalBlock)

      !> Atom for which this is a neighbour
      integer, intent(in) :: iAt

      !> Number of neighbour for this block
      integer, intent(in) :: iNeigh

      !> Number of orbitals on iAt
      integer, intent(in) :: nOrbAt

      !> Number of orbitals on neighbour atom
      integer, intent(in) :: nOrbNeigh

      !> local block
      real(dp), dimension(:), target, intent(inout) :: localBlock

      !> Pointer to local block
      real(dp), dimension(:,:), pointer, intent(out) :: pLocalBlock

      integer :: ind

      ind = iPair(iNeigh, iAt) + 1
      localBlock(1:nOrbNeigh*nOrbAt) = over(ind:ind+nOrbNeigh*nOrbAt-1)
      pLocalBlock(1:nOrbNeigh, 1:nOrbAt) => localBlock(1:nOrbNeigh*nOrbAt)

    end subroutine copyOverlapBlock

    !> copy a density matrix block from sparse matrix
    pure subroutine copyDensityBlock(desc1, desc2, localBlock, pLocalBlock)

      !> start, end and range of first block
      integer, dimension(DESC_LEN), intent(in) :: desc1

      !> start, end and range of second block
      integer, dimension(DESC_LEN), intent(in) :: desc2

      !> local block in 1D format
      real(dp), dimension(:), target, intent(inout) :: localBlock

      !> Pointer to local block
      real(dp), dimension(:,:), pointer, intent(out) :: pLocalBlock

      pLocalBlock(1:desc1(INORB), 1:desc2(INORB)) => localBlock(1:desc1(INORB)*desc2(INORB))
      pLocalBlock(:,:) = tmpDRho(desc1(ISTART):desc1(IEND), desc2(ISTART):desc2(IEND))

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
      integer, dimension(DESC_LEN), intent(in) :: descM

      !> start, end and range of column
      integer, dimension(DESC_LEN), intent(in) :: descN

      !> First overlap block
      real(dp), dimension(:,:), pointer, intent(in) :: pSma

      !> Second overlap block
      real(dp), dimension(:,:), pointer, intent(in) :: pSbN

      !> density matrix block
      real(dp), dimension(:,:), pointer, intent(in) :: pPab

      real(dp), dimension(:,:), pointer :: pHmn

      pHmn => tmpHH(descM(ISTART):descM(IEND), descN(ISTART):descN(IEND))
      if(self%tSpin) then
        pHmn(:,:) = pHmn - 0.25_dp * gammaTot * matmul(matmul(pSma, pPab), pSbn)
      else
        pHmn(:,:) = pHmn - 0.125_dp * gammaTot * matmul(matmul(pSma, pPab), pSbn)
      end if

    end subroutine updateHamiltonianBlock

  end subroutine addLRHamiltonian_nb


  !> Add the LR-Energy contribution to the total energy
  subroutine addLREnergy(self, energy)

    !> RangeSep class instance
    class(RangeSepFunc), intent(inout) :: self

    !> total energy
    real(dp), intent(inout) :: energy

    energy = energy + self%lrenergy
    ! hack for spin unrestricted calculation
    self%lrenergy = 0.0_dp

  end subroutine addLRenergy


  !> copy lower triangle to upper for a square matrix
  subroutine symmetrizeSquareMatrix(matrix)

    !> matrix to symmetrize
    real(dp), allocatable, intent(inout) :: matrix(:,:)

    integer :: ii, matSize

    matSize = size(matrix, dim = 1)
  @:ASSERT(size(matrix, dim = 2) == matSize)
    do ii = 1, matSize
      matrix(ii, ii+1:matSize) = matrix(ii+1:matSize, ii)
    end do

  end subroutine symmetrizeSquareMatrix


  !> location of relevant atomic block indices in a dense matrix
  pure function getDescriptor(iAt, iSquare) result(desc)

    !> relevant atom
    integer, intent(in) :: iAt

    !> indexing array for start of atom orbitals
    integer, intent(in) :: iSquare(:)

    !> resulting location ranges
    integer :: desc(3)

    desc(:) = [ iSquare(iAt), iSquare(iAt + 1) - 1, iSquare(iAt + 1) - iSquare(iAt) ]

  end function getDescriptor


  !> evaluate energy from triangles of the hamitonian and density matrix
  pure function evaluateEnergy(hamiltonian, densityMat) result(egy)

    !> hamiltonian matrix
    real(dp), intent(in) :: hamiltonian(:,:)

    !> density matrix
    real(dp), intent(in) :: densityMat(:,:)

    !> resulting energy
    real(dp) :: egy

    integer :: mu

    egy = 0.0_dp
    do mu = 1, size(hamiltonian, dim = 2)
      egy = egy + hamiltonian(mu,mu)*densityMat(mu,mu)&
          & + 2.0_dp*sum(hamiltonian(mu+1:,mu)*densityMat(mu+1:,mu))
    end do
    egy = 0.5_dp * egy

  end function evaluateEnergy


  !> Interface routine.
  subroutine addLRHamiltonian(self, env, densSqr, over, iNeighbour, nNeighbourLC, iSquare, iPair,&
      & orb, HH, overlap)

    !> class instance
    class(RangeSepFunc), intent(inout) :: self

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
    select case(trim(self%RSAlg))
    case ("tr")
      call addLRHamiltonian_tr(self, env, overlap, densSqr, iSquare, HH, orb)
    case ("nb")
      call addLRHamiltonian_nb(self, env, densSqr, over, iNeighbour, nNeighbourLC, iSquare, iPair,&
          & orb, HH)
    case default
    end select
    call env%globalTimer%stopTimer(globalTimers%rangeSeparatedH)

  end subroutine addLRHamiltonian


  !> Analytical long-range Gamma
  function getAnalyticalGammaValue(self, Sp1, Sp2, dist)

    !> RangeSepFunc instance
    class(RangeSepFunc), intent(inout) :: self

    !> first species
    integer, intent(in) :: Sp1

    !> second species
    integer, intent(in) :: Sp2

    !> distance between atoms
    real(dp), intent(in) :: dist

    !> resulting gamma
    real(dp) :: getAnalyticalGammaValue

    integer :: ii
    real(dp) :: tauA, tauB, omega
    real(dp) :: prefac, tmp, tmp2, tau

    tauA = 3.2_dp * self%hubbu(Sp1)
    tauB = 3.2_dp * self%hubbu(Sp2)
    omega = self%omega

    if (dist < tolSameDist) then
      ! on-site case
      if (abs(tauA - tauB) < MinHubDiff ) then
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
  function getdAnalyticalGammaDeriv(self, Sp1, Sp2, dist)

    !> RangeSepFunc instance
    class(RangeSepFunc), intent(inout) :: self

    !> first species
    integer, intent(in) :: Sp1

    !> second species
    integer, intent(in) :: Sp2

    !> distance between atoms
    real(dp), intent(in) :: dist

    !> resulting d gamma / d dist
    real(dp) :: getdAnalyticalGammaDeriv

    integer :: ii
    real(dp) :: tauA, tauB, omega
    real(dp) :: prefac, tmp, tmp2, tau, dTmp, dTmp2

    tauA = 3.2_dp * self%hubbu(Sp1)
    tauB = 3.2_dp * self%hubbu(Sp2)
    omega = self%omega

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
  subroutine getGammaPrimeValue(self, grad, iAtom1, iAtom2, coords, species)

    !> class instance
    class(RangeSepFunc), intent(inout) :: self

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

    !!> finite difference choice
    !real(dp), parameter :: deltaXDiff = epsilon(1.0_dp)**0.25_dp

    integer :: sp1, sp2, jj, ii
    real(dp) :: vect(3), tmp(3),tmp2(3), dist
    real(dp) :: tauA, tauB, omega

    sp1 = species(iAtom1)
    sp2 = species(iAtom2)

    ! numerical finite difference
    !    do jj = 1, 3 ! x,y,z
    !      tmp(jj) = 0.0_dp
    !      do ii = 1, 2 ! +h, -h
    !        ! difference vector
    !        vect(:) = coords(:,iAtom2) - coords(:,iAtom1)
    !        vect(jj) = vect(jj) - real(2 * ii - 3, dp) * deltaXDiff
    !        dist = sqrt(sum(vect(:)**2))
    !        vect(:) = vect(:) / dist
    !        tmp(jj) = tmp(jj) + real(2 * ii - 3, dp)*getAnalyticalGammaValue(self, sp1, sp2, dist)
    !      end do
    !    end do
    !    tmp(:) = 0.5_dp * tmp(:) / deltaXDiff
    !    grad(:) = tmp(:)

    ! analytical derivatives
    vect(:) = coords(:,iAtom1) - coords(:,iAtom2)
    dist = sqrt(sum(vect(:)**2))
    vect(:) = vect(:) / dist
    grad(:) = vect(:) * getdAnalyticalGammaDeriv(self, sp1, sp2, dist)

  end subroutine getGammaPrimeValue


  !> Adds gradients due to long-range HF-contribution
  subroutine addLRGradients(self, gradients, derivator, deltaRho, skHamCont, skOverCont, coords,&
      & species, orb, iSquare, ovrlapMat, iNeighbour, nNeighbourSK)

    !> class instance
    class(RangeSepFunc), intent(inout) :: self

    !> energy gradients
    real(dp), intent(inout) :: gradients(:,:)

    !> density matrix difference from reference q0
    real(dp), intent(in) :: deltaRho(:,:)

    !> sparse hamiltonian (non-scc)
    type(OSlakoCont), intent(in) :: skHamCont

    !> sparse overlap part
    type(OSlakoCont), intent(in) :: skOverCont

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
    class(NonSccDiff), intent(in) :: derivator

    integer :: nAtom, iAtK, iNeighK, iAtB, iNeighB, iAtC, iAtA, kpa
    real(dp) :: tmpgamma1, tmpgamma2
    real(dp) :: tmpforce(3), tmpforce_r(3), tmpforce2, tmpmultvar1
    integer :: mu, nu, alpha, beta, ccc, kkk
    real(dp) :: dummy(orb%mOrb,orb%mOrb,3), sPrimeTmp(orb%mOrb,orb%mOrb,3)
    real(dp) :: sPrimeTmp2(orb%mOrb,orb%mOrb,3)
    real(dp), allocatable :: gammaPrimeTmp(:,:,:), tmpovr(:,:), tmpRho(:,:), tmpderiv(:,:)

    write(stdOut,'(a)') "rangeSep: addLRGradients"
    @:ASSERT(size(gradients,dim=1) == 3)
    call allocateAndInit(tmpovr, tmpRho, gammaPrimeTmp, tmpderiv)
    nAtom = size(self%species)
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
            tmpgamma1 = self%lrGammaEval(iAtK,iAtB) + self%lrGammaEval(iAtC,iAtB)
            tmpgamma2 = tmpgamma1 + self%lrGammaEval(iAtK,iAtA) + self%lrGammaEval(iAtC,iAtA)
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
                do alpha = iSquare(iAtA), iSquare(iAtA + 1) - 1
                  do beta = iSquare(iAtB), iSquare(iAtB + 1) - 1
                    tmpmultvar1 = tmpmultvar1 + tmpovr(beta, alpha) * (tmpRho(beta,kpa) &
                        & * tmpRho(alpha,mu) + tmpRho(alpha,kpa) * tmpRho(beta,mu))
                  end do
                end do
                tmpforce(:) = tmpforce(:) + tmpmultvar1 * (sPrimeTmp(ccc,kkk,:))
                tmpforce_r(:) = tmpforce_r(:) + tmpmultvar1 * (sPrimeTmp2(kkk,ccc,:))
                tmpforce2 = tmpforce2 + tmpmultvar1 * tmpovr(kpa,mu)
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

    gradients(:,:) = gradients -0.25_dp * tmpderiv

    deallocate(tmpovr, tmpRho, gammaPrimeTmp, tmpderiv)

  contains

    !> Initialise the
    subroutine allocateAndInit(tmpovr, tmpRho, gammaPrimeTmp, tmpderiv)

      !> Storage for the overlap
      real(dp), allocatable, intent(inout) :: tmpovr(:,:)

      !> storage for density matrix
      real(dp), allocatable, intent(inout) :: tmpRho(:,:)

      !> storage for derivative of gamma interaction
      real(dp), allocatable, intent(inout) :: gammaPrimeTmp(:,:,:)

      !> workspace for the derivatives
      real(dp), allocatable, intent(inout) :: tmpderiv(:,:)

      real(dp) :: tmp(3)
      integer :: iAt1, iAt2, nAtom

      nAtom = size(self%species)
      allocate(tmpovr(size(ovrlapMat, dim = 1), size(ovrlapMat, dim = 1)))
      allocate(tmpRho(size(deltaRho, dim = 1), size(deltaRho, dim = 1)))
      allocate(gammaPrimeTmp(3, nAtom, nAtom))
      allocate(tmpderiv(3, size(gradients, dim = 2)))
      tmpovr = ovrlapMat
      tmpRho = deltaRho
      call symmetrizeSquareMatrix(tmpovr)
      call symmetrizeSquareMatrix(tmpRho)
      ! precompute the gamma derivatives
      write(stdOut,'(a)') "precomputing the lr-gamma derivatives"
      gammaPrimeTmp = 0.0_dp
      do iAt1 = 1, nAtom
        do iAt2 = 1, nAtom
          if(iAt1 /= iAt2) then
            call getGammaPrimeValue(self, tmp, iAt1, iAt2, coords, species)
            gammaPrimeTmp(:,iAt1, iAt2) = tmp(:)
          end if
        end do
      end do
    end subroutine allocateAndInit

  end subroutine addLRGradients


  !> evaluate the LR-Energy contribution directly. Very slow, use addLREnergy instead.
  function evaluateLREnergyDirect(self, env, deltaRho, ovrlap, iSquare) result (energy)

    !> instance of LR
    class(RangeSepFunc), intent(inout) :: self

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
    real(dp), allocatable :: tmpovr(:,:), tmpDRho(:,:)
    real(dp) :: tmp

    call env%globalTimer%startTimer(globalTimers%energyEval)
    nAtom = size(self%species)
    allocate(tmpovr(size(ovrlap, dim = 1), size(ovrlap, dim = 1)))
    allocate(tmpDRho(size(deltaRho, dim = 1), size(deltaRho, dim = 1)))
    tmpovr = ovrlap
    tmpDRho = deltaRho
    call symmetrizeSquareMatrix(tmpovr)
    call symmetrizeSquareMatrix(tmpDRho)

    energy = 0.0_dp
    do iAt1 = 1, nAtom
      do iAt2 = 1, nAtom
        tmp = 0.0_dp
        do mu = iSquare(iAt1), iSquare(iAt1 + 1) - 1
          do nu = iSquare(iAt2), iSquare(iAt2 + 1) - 1
            do alpha = 1, size(tmpovr, dim = 1)
              do beta = 1, size(tmpovr, dim = 1)
                tmp = tmp + (&
                    & tmpDRho(alpha,beta) * tmpDRho(mu,nu)&
                    & +tmpDRho(mu,beta) * tmpDRho(alpha,nu)) * tmpovr(mu,alpha) * tmpovr(nu,beta)
              end do
            end do
          end do
        end do
        energy = energy + tmp * self%lrGammaEval(iAt1,iAt2)
      end do
    end do
    energy = -energy / 8.0_dp

    call env%globalTimer%stopTimer(globalTimers%energyEval)

  end function evaluateLREnergyDirect

end module rangeseparated
