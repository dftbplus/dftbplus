
#:include 'common.fypp'


!> Contains range separated related routines.
module rangeseparated 
  use accuracy
  use assert
  use nonscc, only : NonSccDiff
  use periodic, only : getLatticePoints
  use SlakoCont, only : OSlakoCont, getMIntegrals, getSKIntegrals
  use SK, only : rotateH0
  use CommonTypes
  use Interpolation
  use sorting
  use sparse2dense, only : blockSymmetrizeHS
  implicit none
  private

  public :: RangeSepFunc, TRangeSepSKTag
  real(dp), pointer :: gLatPoint_(:,:)     ! Lattice points for reciprocal Ewald

  !> Slater-Koster file RangeSep tag structure
  type :: TRangeSepSKTag
     character(lc) :: type 
     real(dp) :: omega
  end type TRangeSepSKTag

  !> Range-Sep module structure
  type :: RangeSepFunc 

     real(dp), allocatable :: coords(:,:)
     real(dp), allocatable :: lrGamma(:,:,:,:) ! lrGamma(spec1, spec2, dist, (1,2))
     real(dp), allocatable :: lrGammaEval(:,:) ! Atom1, Atom2 at each geometry step
     real(dp), allocatable :: grid(:)          ! 
     !
     real(dp) :: omega                         ! range-separation parameter
     real(dp) :: gamma_range
     integer :: num_mesh_points
     real(dp) :: delta
     real(dp), allocatable :: hubbu(:)
     !
     ! Hamiltonian Screening
     real(dp), allocatable :: hprev(:,:)
     real(dp), allocatable :: dRhoprev(:,:)
     logical :: tScreeningInited
     real(dp) :: pScreeningTreshold
     !
     ! lr-energy
     real(dp) :: lrenergy
     real(dp) :: lrenergyUp
     real(dp) :: lrenergyDn
     ! periodic boundary
     logical :: tPeriodic
     integer :: nKPoints
     !
     logical :: tSpin
     logical :: tTabGamma
     character(lc) :: RSAlg
     !
     integer, allocatable :: species(:)
   contains
     procedure :: initModule
     procedure :: updateCoords
     procedure :: addLRHamiltonian
     procedure :: addLRHamiltonian_tr
     procedure :: addLRHamiltonian_nb
     procedure :: addLREnergy
     procedure :: addLRGradients
     procedure :: getRSAlg
     procedure :: addLRHamiltonian_pbc
     procedure :: evaluateLREnergyDirect
     !
     ! procedure :: addSquareHamiltonian
     ! procedure :: addSquareHamiltonian_opt
     ! procedure :: addGradients
     ! procedure :: addGradients_tst
     ! procedure :: addGradients_tst2
     ! procedure :: getGammaValue
     ! procedure :: updateHamiltonian_thr
  end type RangeSepFunc

contains

  !> Return the RSAlg value
  !!  \param self, class instance
  function getRSALg(self) result(res)
    class(RangeSepFunc), intent(in) :: self
    character(lc) :: res
    res = self%RSAlg
  end function getRSALg

  !> Intitialize the range-sep module
  !!  \param self, class instance
  !!  \param nAtom, number of atoms
  !!  \param species, list of all atomic species
  !!  \param speciesNames, list of all atomic species names
  !!  \param hubbu, atomic hubbards
  !!  \param screen screening threshold value
  !!  \param omega range separation parameter
  !!  \param nKPts number of K-points for the calculation (1 for cluster)
  !!  \param tSpin spin unrestricted?
  !!  \param tTabGamma tabulated or analytical lr-gamma?
  !!  \param RSAlg lr-hamiltonian construction algorithm
  subroutine initModule(self, nAtom, species, speciesNames,hubbu,screen,omega,nKPts,tSpin,tTabGamma,RSAlg)
    class(RangeSepFunc), intent(inout) :: self
    integer, intent(in) :: nAtom
    integer, intent(in) :: species(:)
    character(mc), intent(in) :: speciesNames(:)
    ! 
    real(dp), intent(in) :: hubbu(:)
    real(dp), intent(in) :: screen
    real(dp), intent(in) :: omega
    integer, intent(in) :: nKPts
    logical, intent(in) :: tSpin
    logical, intent(in) :: tTabGamma
    character(lc), intent(in) :: RSAlg

    call initAndAllocate(self, nAtom, hubbu, species, screen, omega, RSAlg, tSpin, tTabGamma)
    call printModuleInfoAndCheckReqs(self, nKPts)
    if(tTabGamma) then
       call loadAndProcessTabulatedLRGammas(self, speciesNames, species)
    end if

  contains

    subroutine initAndAllocate(self, nAtom, hubbu, species, screen, omega, RSAlg, tSpin, tTabGamma)
      class(RangeSepFunc), intent(inout) :: self
      integer, intent(in) :: nAtom
      real(dp), intent(in) :: hubbu(:)
      integer, intent(in) :: species(:)
      real(dp), intent(in) :: screen
      real(dp), intent(in) :: omega
      character(lc), intent(in) :: RSAlg
      logical, intent(in) :: tSpin
      logical, intent(in) :: tTabGamma
      !
      self%tScreeningInited = .false.
      self%pScreeningTreshold = screen
      self%omega = omega
      self%lrenergy = 0.0_dp
      self%RSAlg = RSAlg
      self%tSpin = tSpin
      self%tTabGamma = tTabGamma

      allocate(self%coords(3, nAtom))
      allocate(self%species(nAtom))
      allocate(self%lrGammaEval(nAtom,nAtom))
      allocate(self%hubbu(size(hubbu(:))))
      self%hubbu = hubbu
      self%species(:) = species
    end subroutine initAndAllocate

    subroutine printModuleInfoAndCheckReqs(self, nKPts)
      class(RangeSepFunc), intent(inout) :: self
      integer, intent(in) :: nKPts
      !
      write(*,'(a)') "================================"
      write(*,'(a)') "Range-separated Hybrids in DFTB "
      write(*,'(a)') ""
      write(*,'(a)') "================================"
      write(*,'(a)') "=> Initializing RangeSep module"
      
      ! Current restrictions
      if(self%tSpin .and. self%RSAlg == "tr") then
         write(*,'(a)') "ERROR! Spin-unrestricted calculation for threshoding algortihm not yet implemented!"
         stop
      end if

      if(nKPts > 1) then
         write(*,'(a)') "ERROR! PBC not yet implemented!"
         stop
      end if

      ! summarize module settings
      write(*,'(a,F10.3)') "  -> range-separation parameter [1/a0]:", self%omega

      if(self%tSpin) then
         write(*,'(a)') "  -> spin-unrestricted calculation"
      else
         write(*,'(a)') "  -> spin-restricted calculation"
      end if

      if(self%tTabGamma) then
         write(*,'(a)') "  -> use the tabulated long-range Gamma"
      else
         write(*,'(a)') "  -> use analytical long-range Gamma"
      end if

      select case(self%RSAlg)
      case ("nb") 
         write(*,'(a)') "  -> using the neighbour list-based algorithm"
      case ("tr")
         write(*,'(a)') "  -> using the thresholding algorithm"
         write(*,'(a,E17.8)') "     -> Screening Threshold:", self%pScreeningTreshold
      case default
         write(*,'(a)') " ERROR! Invalid algorithm"
         stop
      end select

    end subroutine printModuleInfoAndCheckReqs

    subroutine loadAndProcessTabulatedLRGammas(self, speciesNames, species)
      class(RangeSepFunc), intent(inout) :: self
      character(mc), intent(in) :: speciesNames(:)
      integer, intent(in) :: species(:)

      character(lc) :: fname, dummy
      integer :: ii, jj, kk, st
      logical :: tFileExists
      integer :: num_mesh_points
      real(dp) :: delta,tmp, dist
      real(dp), allocatable :: test_array2(:)

      num_mesh_points=500
      delta=0.02_dp
      do ii = 1, size(speciesNames)
         do jj = 1, size(speciesNames)
            fname = "gamma-" // trim(speciesNames(ii)) // "-" &
                 & // trim(speciesNames(jj))
            print *, "process the heaader of '" // trim(fname) // ".dat'"
            fname = trim(fname) // ".dat"
            !
            inquire (file=fname, exist=tFileExists, iostat=st)
            if(tFileExists) then
               open(95,FILE=trim(fname),FORM='formatted',STATUS='unknown')
               read(95,*) dummy, self%omega, self%num_mesh_points, self%delta
               close(95)
               ! 
               write(*,'(a,E17.8)') "RS-Parameter k=", self%omega
               self%gamma_range=real(self%num_mesh_points,dp)*self%delta
               write(*,'(a,I6,2E17.8)') "gamma mesh parameters:", self&
                    &%num_mesh_points, self%delta, self%gamma_range
               num_mesh_points=self%num_mesh_points
               delta=self%delta
            else
               write(*,'(a)') "RangeSep: Tabulated long-range Gamma file does not exist!"
               stop
            end if
         end do
      end do
      !
      write(*,'(a,3I5)') "allocating the array for lr-gamma, dimensions:",&
           & size(speciesNames),size(speciesNames), num_mesh_points
      allocate(self%lrGamma(size(speciesNames),size(speciesNames)&
           &,num_mesh_points,2))
      allocate(self%grid(num_mesh_points))
      allocate(test_array2(num_mesh_points))
      !
      self%species(:) = species
      do ii = 1, size(speciesNames)
         do jj = 1, size(speciesNames)
            fname = "gamma-" // trim(speciesNames(ii)) // "-" &
                 & // trim(speciesNames(jj))
            print *, "Need to process '" // trim(fname) // ".dat'"
            fname = trim(fname) // ".dat"
            !
            open(95,FILE=trim(fname),FORM='formatted',STATUS='unknown')
            read(95,*)
            do kk=1,num_mesh_points
               read(95,*) self%grid(kk),self%lrGamma(ii,jj,kk,1)
            end do
            close(95)
            ! 
            write(*,'(a)') "prepare the interpolation"
            call set_cubic_spline(self%grid,self%lrGamma(ii,jj,:,1),test_array2)
            self%lrGamma(ii,jj,:,2)=test_array2
            write(*,'(a)') "Done!"
         end do
      end do

    end subroutine loadAndProcessTabulatedLRGammas

  end subroutine initModule

  !> update the rangeSep module on coordinate change
  !!  \param self, class instance
  !!  \param coords, list of atomic coordinates
  subroutine updateCoords(self, coords)
    class(RangeSepFunc), intent(inout) :: self
    real(dp), intent(in) :: coords(:,:)

    integer :: nAtom, iAtom1,iAtom2,ii,iSp1,iSp2
    real(dp) :: dist

    write(*,*) "rangesep update coords entry"
    @:ASSERT(all(shape(coords) == shape(self%coords)))
    self%coords(:,:) = coords
    nAtom = size(self%species)
    write(*,*) "Evaluating LR-GAMMA for atom pairs"
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
          self%lrGammaEval(iAtom1, iAtom2) = getGammaValue(self, iSp1, iSp2, dist)
          self%lrGammaEval(iAtom2, iAtom1) = self%lrGammaEval(iAtom1, iAtom2)
       end do
    end do
    print *, "Range separated coordinates updated"
    ! reinit the screening
    if( self%tScreeningInited) then
       self%hprev = 0.0_dp
       self%dRhoprev = 0.0_dp
       !
       self%lrenergy = 0.0_dp
    end if
  end subroutine updateCoords

  !> Adds the LR-exchange contribution to hamiltonian using the threshold algorithm
  !! \param self, class instance
  !! \param overlap, square real overlap matrix
  !! \param deltaRho, square density matrix (deltaRho in DFTB terms)
  !! \param iSquare, mapping atom_number -> number of the first basis function of the
  !!                 atomic block atom_number
  !! \param hamiltonian, current hamiltonian
  !! \param orb
  subroutine addLRHamiltonian_tr(self, overlap, deltaRho, iSquare, hamiltonian, orb)
    class(RangeSepFunc), intent(inout) :: self
    real(dp), intent(in) :: overlap(:,:), deltaRho(:,:)
    integer, intent(in) :: iSquare(:)
    real(dp), intent(inout) :: hamiltonian(:,:)
    type(TOrbitals), intent(in) :: orb
    !
    real(dp), allocatable :: tmpovr(:,:), tmpDRho(:,:), testovr(:,:), tmpDDRho(:,:), tmpDham(:,:)
    real :: start, finish
    integer, allocatable :: ovrind(:,:)
    integer, parameter :: DESC_LEN = 3, ISTART = 1, IEND = 2, INORB = 3

    call cpu_time(start)
    call allocateAndInit(tmpovr, tmpDham, tmpDRho, tmpDDRho, testovr, ovrind)
    call evaluateHamiltonian(tmpDHam)
    self%hprev = self%hprev + tmpDham
    hamiltonian = hamiltonian + self%hprev
    self%lrenergy = evaluateEnergy()
    call deallocateTempMatrices()
    call cpu_time(finish)
    !print '("--> Done. time = ", f10.4, " sec.")', finish - start  

  contains

    subroutine allocateAndInit(tmpovr, tmpDham, tmpDRho, tmpDDRho, testovr, ovrind)
      real(dp), allocatable, intent(inout) :: tmpovr(:,:), tmpDham(:,:),  tmpDRho(:,:), tmpDDRho(:,:)
      real(dp), allocatable, intent(inout) :: testovr(:,:)
      integer, allocatable :: ovrind(:,:)
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
            tmp = maxval(abs(tmpovr((iSquare(iAtMu)):(iSquare(iAtMu + 1) - 1)&
                 &,(iSquare(iAtNu)):(iSquare(iAtNu + 1) - 1))))
            testovr(iAtMu,iAtNu) = tmp
         end do
      end do
      do iAtMu = 1, nAtom
         call index_heap_sort(ovrind(iAtMu,:),testovr(iAtMu,:))
      end do
    end subroutine allocateAndInit

    subroutine evaluateHamiltonian(tmpDHam)
      real(dp), allocatable, intent(inout) :: tmpDHam(:,:)
      integer :: nAtom
      real(dp) :: pbound, prb
      real(dp) :: tmpvec1(orb%mOrb), tmpvec2(orb%mOrb), tmpvec3(orb%mOrb)
      real(dp) :: tmp, tstbound, gammabatch, gammabatchtmp
      integer :: iAtMu, iAtNu, iAt1, iAt2, iSp1, iSp2, nOrb1, nOrb2
      integer :: kk, ll, jj, ii, mu, nu
      integer, dimension(DESC_LEN) :: descA, descB, descM, descN

      nAtom = size(self%species)
      pbound = getMaxDensElement()
      tmpDham = 0.0_dp
      loopMu: do iAtMu = 1, nAtom
         descM = getDescriptor(iAtMu)
         loopKK: do kk = 1, nAtom
            iAt1 = ovrind(iAtMu, nAtom + 1 - kk)
            descA = getDescriptor(iAt1)
            iSp1 = self%species(iAt1)
            nOrb1 = orb%nOrbSpecies(iSp1)
            prb = pbound * testovr(iAt1, iAtMu)
            if(abs(prb) >= self%pScreeningTreshold) then
               loopNu: do iAtNu = 1, iAtMu
                  descN = getDescriptor(iAtNu)
                  gammabatchtmp = self%lrGammaEval(iAtMu, iAtNu) + self%lrGammaEval(iAt1, iAtNu)
                  loopLL: do ll = 1, nAtom
                     iAt2 = ovrind(iAtNu, nAtom + 1 - ll)
                     iSp2 = self%species(iAt2)
                     nOrb2 = orb%nOrbSpecies(iSp2)
                     ! screening condition
                     tstbound = prb * testovr(iAt2, iAtNu)
                     if(abs(tstbound) >= self%pScreeningTreshold) then
                        descB = getDescriptor(iAt2)
                        gammabatch = (self%lrGammaEval(iAtMu, iAt2)&
                             &+ self%lrGammaEval(iAt1, iAt2)&
                             &+ gammabatchtmp)
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

    function evaluateEnergy() result(energy)
      integer :: nAtom, iAtMu, iAtNu, mu, nu
      real(dp) :: tmp, energy

      nAtom = size(self%species)
      tmp = 0.0_dp
      do iAtMu = 1, nAtom
         do iAtNu = 1, iAtMu
            do mu = iSquare(iAtMu), iSquare(iAtMu + 1) - 1
               do nu = iSquare(iAtNu), iSquare(iAtNu + 1) - 1
                  if(iAtNu == iAtMu) then
                     tmp = tmp + self%hprev(mu, nu) * tmpDRho(mu, nu)
                  else
                     tmp = tmp + 2.0_dp * self%hprev(mu, nu) * tmpDRho(mu, nu)
                  end if
               end do
            end do
         end do
      end do
      energy = 0.5_dp * tmp
    end function evaluateEnergy

    subroutine deallocateTempMatrices()
      deallocate(tmpovr, tmpDham, tmpDRho,tmpDDRho)
      deallocate(testovr, ovrind)
    end subroutine deallocateTempMatrices

    subroutine checkAndInitScreening(self, matrixSize, tmpDRho)
      class(RangeSepFunc), intent(inout) :: self
      integer, intent(in) :: matrixSize
      real(dp), allocatable, intent(in) :: tmpDRho(:,:)

      if(.not. self%tScreeningInited) then
         write(*,'(a)') "RangeSep: Iinitialize Screening"
         self%tScreeningInited = .true.
         allocate(self%hprev(matrixSize, matrixSize))
         allocate(self%dRhoprev(matrixSize, matrixSize))
         self%hprev = 0.0_dp
         self%dRhoprev = tmpDRho
      end if
    end subroutine checkAndInitScreening

    function getMaxDensElement() result(pbound)
      real(dp) :: pbound
      pbound = maxval(abs(tmpDDRho))
    end function getMaxDensElement
    
    function getDescriptor(iAt) result(desc)
      integer, intent(in) :: iAt
      integer, dimension(DESC_LEN) :: desc

      desc(:) = [ iSquare(iAt), iSquare(iAt + 1) - 1, &
          & iSquare(iAt + 1) - iSquare(iAt) ]
    end function getDescriptor

  end subroutine addLRHamiltonian_tr

  !> Add the LR-Energy contribution to the total energy
  !! \param self, RangeSep class instance
  !! \param energy, total energy
  subroutine addLREnergy(self, energy)
    class(RangeSepFunc), intent(inout) :: self
    real(dp), intent(inout) :: energy    

    energy = energy + self%lrenergy
    ! hack for spin unrestricted calculation
    self%lrenergy = 0.0_dp
  end subroutine addLRenergy

  subroutine symmetrizeSquareMatrix(matrix)
    real(dp), allocatable, intent(inout) :: matrix(:,:)
    integer :: ii, jj, matSize

    matSize = size(matrix, dim = 1)
    if( matSize /= size(matrix, dim = 2)) then
       write(*,'(a)') "Error(RangeSep): not a square matrix"
       stop
    end if
    do ii = 1, matSize
       do jj = ii + 1, matSize
          matrix(ii, jj) = matrix(jj, ii)
       end do
    end do
  end subroutine symmetrizeSquareMatrix

  !> Updates the Hamiltonian with the range separated contribution.
  !!
  !! \param densSqr  Square (unpacked) density matrix
  !! \param over  Sparse (packed) overlap matrix.
  !! \param iNeighbor  Neighbor indices.
  !! \param nNeighbor  Nr. of neighbors for each atom.
  !! \param iSquare  Position of each atom in the rows/columns of the square
  !!     matrices. Shape: (nAtom)
  !! \param iPair  Position of each (neighbor, atom) pair in the sparse
  !!   matrix. Shape: (0:maxNeighbor, nAtom)
  !! \param orb  Orbital information.
  !! \param HH  Square (unpacked) Hamiltonian to be updated.
  subroutine addLRHamiltonian_nb(self, densSqr, over, iNeighbor, nNeighbor, iSquare, iPair, orb, HH)
    class(RangeSepFunc), intent(inout) :: self
    real(dp), dimension(:,:), target, intent(in) :: densSqr
    real(dp), dimension(:), intent(in) :: over
    integer, dimension(0:,:), intent(in) :: iNeighbor
    integer, dimension(:), intent(in) :: nNeighbor
    integer, dimension(:), intent(in) :: iSquare
    integer, dimension(0:,:), intent(in) :: iPair
    type(TOrbitals), intent(in) :: orb
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
    real :: start, finish

    call cpu_time(start)
    call allocateAndInit(tmpHH, tmpDRho)
    call evaluateHamiltonian()
    HH = HH + tmpHH
    call evaluateEnergy()
    call cpu_time(finish)
    print '("--> Done. time = ", f10.4, " sec.")', finish - start  
    
  contains

    subroutine allocateAndInit(tmpHH, tmpDRho)
      real(dp), dimension(:,:), allocatable, target, intent(inout) :: tmpDRho, tmpHH
      allocate(tmpHH(size(HH, dim = 1), size(HH, dim = 2)))
      tmpHH = 0.0_dp
      allocate(tmpDRho(size(densSqr, dim = 1), size(densSqr, dim = 1)))
      tmpDRho = densSqr
      call symmetrizeSquareMatrix(tmpDRho)
    end subroutine allocateAndInit

    subroutine evaluateHamiltonian()
      nAtom = size(self%species)
      loopN: do iAtN = 1, nAtom
         descN = getDescriptor(iAtN)
         loopB: do iNeighN = 0, nNeighbor(iAtN)
            iAtB = iNeighbor(iNeighN, iAtN)
            descB = getDescriptor(iAtB)
            call copyOverlapBlock(iAtN, iNeighN, descN(INORB), descB(INORB), &
                 & Sbn, pSbn)
            call transposeBlock(pSbn, Snb, pSnb)
            loopA: do iAtA = 1, nAtom
               descA = getDescriptor(iAtA)
               call copyDensityBlock(descA, descB, Pab, pPab)
               call copyDensityBlock(descA, descN, Pan, pPan)
               gamma1 = self%lrGammaEval(iAtA, iAtN) + self%lrGammaEval(iAtA, iAtB)
               loopM: do iNeighA = 0, nNeighbor(iAtA)
                  iAtM = iNeighbor(iNeighA, iAtA)
                  descM = getDescriptor(iAtM)
                  call copyOverlapBlock(iAtA, iNeighA, descA(INORB), descM(INORB), &
                       & Sma, pSma)
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

    subroutine evaluateEnergy()
      tmp = 0.0_dp
      do iAtM = 1, nAtom
         do iAtN = 1, iAtM
            do mu = iSquare(iAtM), iSquare(iAtM + 1) - 1
               do nu = iSquare(iAtN), iSquare(iAtN + 1) - 1
                  if(iAtN == iAtM) then
                     tmp = tmp + tmpHH(mu, nu) * tmpDRho(mu, nu)
                  else
                     tmp = tmp + 2.0_dp * tmpHH(mu, nu) * tmpDRho(mu, nu)
                  end if
               end do
            end do
         end do
      end do
      self%lrenergy = self%lrenergy + 0.5_dp * tmp
    end subroutine evaluateEnergy

    function getDescriptor(iAt) result(desc)
      integer, intent(in) :: iAt
      integer, dimension(DESC_LEN) :: desc

      desc(:) = [ iSquare(iAt), iSquare(iAt + 1) - 1, &
          & iSquare(iAt + 1) - iSquare(iAt) ]
    end function getDescriptor

    subroutine copyOverlapBlock(iAt, iNeigh, nOrbAt, nOrbNeigh, localBlock, &
        & pLocalBlock) 
      integer, intent(in) :: iAt, iNeigh, nOrbAt, nOrbNeigh
      real(dp), dimension(:), target, intent(inout) :: localBlock
      real(dp), dimension(:,:), pointer, intent(out) :: pLocalBlock
      
      integer :: ind

      ind = iPair(iNeigh, iAt) + 1
      localBlock(1:nOrbNeigh*nOrbAt) = over(ind:ind+nOrbNeigh*nOrbAt-1)
      pLocalBlock(1:nOrbNeigh, 1:nOrbAt) => localBlock(1:nOrbNeigh*nOrbAt)
    end subroutine copyOverlapBlock

    subroutine copyDensityBlock(desc1, desc2, localBlock, pLocalBlock)
      integer, dimension(DESC_LEN), intent(in) :: desc1, desc2
      real(dp), dimension(:), target, intent(inout) :: localBlock
      real(dp), dimension(:,:), pointer, intent(out) :: pLocalBlock
      
      pLocalBlock(1:desc1(INORB), 1:desc2(INORB)) => &
          & localBlock(1:desc1(INORB)*desc2(INORB))
!      pLocalBlock(:,:) = &
!          & densSqr(desc1(ISTART):desc1(IEND), desc2(ISTART):desc2(IEND))
      pLocalBlock(:,:) = &
          & tmpDRho(desc1(ISTART):desc1(IEND), desc2(ISTART):desc2(IEND))
    end subroutine copyDensityBlock

    subroutine transposeBlock(orig, localBlock, pLocalBlock)
      real(dp), dimension(:,:), intent(in) :: orig
      real(dp), dimension(:), target, intent(out) :: localBlock
      real(dp), dimension(:,:), pointer, intent(out) :: pLocalBlock

      pLocalBlock(1:size(orig, dim=2), 1:size(orig, dim=1)) => &
          & localBlock(1:size(orig))
      pLocalBlock = transpose(orig)
    end subroutine transposeBlock

    subroutine updateHamiltonianBlock(descM, descN, pSma, pSbN, pPab)
      integer, dimension(DESC_LEN), intent(in) :: descM, descN
      real(dp), dimension(:,:), pointer, intent(in) :: pSma, pSbN, pPab

      real(dp), dimension(:,:), pointer :: pHmn

!      pHmn => HH(descM(ISTART):descM(IEND), descN(ISTART):descN(IEND))
      pHmn => tmpHH(descM(ISTART):descM(IEND), descN(ISTART):descN(IEND))
      if(self%tSpin) then
         pHmn(:,:) = pHmn + gammaTot * matmul(matmul(pSma, pPab), pSbn) * (-0.25_dp) 
      else
         pHmn(:,:) = pHmn + gammaTot * matmul(matmul(pSma, pPab), pSbn) * (-0.125_dp)
      end if
    end subroutine updateHamiltonianBlock

  end subroutine addLRHamiltonian_nb

  !> Interface routine. 
  !!
  !! \param self, class instance
  !! \param densSqr  Square (unpacked) density matrix
  !! \param over  Sparse (packed) overlap matrix.
  !! \param iNeighbor  Neighbor indices.
  !! \param nNeighbor  Nr. of neighbors for each atom.
  !! \param iSquare  Position of each atom in the rows/columns of the square
  !!     matrices. Shape: (nAtom)
  !! \param iPair  Position of each (neighbor, atom) pair in the sparse
  !!   matrix. Shape: (0:maxNeighbor, nAtom)
  !! \param orb  Orbital information.
  !! \param HH  Square (unpacked) Hamiltonian to be updated.
  !! \param overlap, square real overlap matrix
  !! \param deltaRho, square density matrix (deltaRho in DFTB terms)
  subroutine addLRHamiltonian(self, densSqr, over, iNeighbor, nNeighbor, iSquare, iPair, orb, HH,&
    &overlap, deltaRho)
    class(RangeSepFunc), intent(inout) :: self
    ! NB
    real(dp), dimension(:,:), target, intent(in) :: densSqr
    real(dp), dimension(:), intent(in) :: over
    integer, dimension(0:,:), intent(in) :: iNeighbor
    integer, dimension(:), intent(in) :: nNeighbor
    integer, dimension(:), intent(in) :: iSquare
    integer, dimension(0:,:), intent(in) :: iPair
    type(TOrbitals), intent(in) :: orb
    real(dp), dimension(:,:), intent(inout), target :: HH
    ! TR
    real(dp), intent(in) :: overlap(:,:), deltaRho(:,:)

    select case(self%getRSAlg())
    case ("tr")
       call self%addLRHamiltonian_tr(overlap, deltaRho, iSquare, HH, orb)
    case ("nb")
       call self%addLRHamiltonian_nb(densSqr, over, iNeighbor, nNeighbor, iSquare, iPair, orb, HH)
    case default
    end select
  end subroutine addLRHamiltonian

  !!> Analytical long-range Gamma
  !! \param self, RangeSepFunc instatnce
  !! \param Sp1, Sp2, species
  !! \param dist, separation of sites 1 and 2
  function getAnalyticalGammaValue(self, Sp1, Sp2, dist)
    class(RangeSepFunc), intent(inout) :: self
    integer, intent(in) :: Sp1, Sp2
    real(dp), intent(in) :: dist
    real(dp) :: getAnalyticalGammaValue
    !
    integer :: ii
    real(dp) :: tauA, tauB, omega
    real(dp) :: prefac, tmp, tmp2, tau
    !
    tauA = 3.2_dp * self%hubbu(Sp1)
    tauB = 3.2_dp * self%hubbu(Sp2)
    omega = self%omega
    ! 
    if (dist < tolSameDist) then
       ! on-site case
       if (abs(tauA - tauB) < MinHubDiff ) then
          tau = 0.5_dp * (tauA + tauB)
          tmp = 5.0_dp * tau**6 + 15.0_dp * tau**4 * omega**2 - 5.0_dp * tau**2 * omega**4  + omega**6
          tmp = tmp * 0.0625_dp / tau**5 - omega
          tmp = tmp * tau**8 / (tau**2 - omega**2)**4 
          getAnalyticalGammaValue = tau * 0.3125_dp - tmp
       else
          write(*,'(a)') "Error(RangeSep): R = 0, Ua != Ub"
       end if
    else
       ! off-site case, Ua == Ub
       if (abs(tauA - tauB) < MinHubDiff ) then
          tauA = 0.5_dp * (tauA + tauB)
          tmp2 = ((dist * tauA)**3 / 48.0_dp + 0.1875_dp * (dist * tauA)**2 + &
               & 0.6875_dp * (dist * tauA) + 1.0_dp) * exp(-tauA * dist) / dist
          tmp = -tauA**8 / (tauA**2 - omega**2)**4 * (tmp2 + exp(-tauA*dist) * &
               & (dist**2 * (3.0_dp * tauA**4 * omega**4 - 3.0_dp * tauA**6 * omega**2 - &
               & tauA**2 * omega**6) + dist * (15.0_dp * tauA**3 * omega**4 - &
               & 21.0_dp * tauA**5 * omega**2 - 3.0_dp * tauA * omega**6) + &
               & (15.0_dp * tauA**2 * omega**4 - 45.0_dp * tauA**4 * omega**2 - &
               & 3.0_dp * omega**6)) / (48.0_dp * tauA**5))
          getAnalyticalGammaValue = 1/dist - tmp2 - (tauA**8 / (tauA**2 - omega**2)**4 *&
               & exp(-omega * dist) / dist + tmp)
       else
          ! off-site, Ua != Ub
          prefac = tauA**4 / (tauA * tauA - omega * omega )**2
          prefac = prefac * tauB**4 / (tauB * tauB - omega * omega )**2
          prefac = prefac * exp(-omega * dist) / dist
          tmp = prefac - getYGammaSubPart(tauA,tauB,dist,omega) - getYGammaSubPart(tauB,tauA,dist,omega)
          tmp = 1 / dist - tmp  
          tmp = tmp - getYGammaSubPart(tauA,tauB,dist,0.0_dp) - getYGammaSubPart(tauB,tauA,dist,0.0_dp)
          getAnalyticalGammaValue = tmp 
       end if
    end if
  end function getAnalyticalGammaValue

  !> returns the subexpression for the evaluation 
  !!  of the off-site Y-Gamma-integral
  !! \param tauA decay constant site A
  !! \param tauB decay constant site B
  !! \param R separation of the sites A and B
  !! \param omega, range-separation parameter
  function getYGammaSubPart(tauA, tauB, R, omega)
    real(dp), intent(in) :: tauA, tauB, R, omega

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

  !> returns the long-range gamma
  !! \param self, class instance
  !! \param Sp1, species 1
  !! \param Sp1, species 2
  !! \param dist, distance
  function getGammaValue(self, Sp1, Sp2, dist)
    class(RangeSepFunc), intent(inout) :: self
    integer, intent(in) :: Sp1, Sp2
    real(dp), intent(in) :: dist

    real(dp) :: getGammaValue
    real(dp) :: tmp, tmp2

    if(self%tTabGamma) then
       ! This option is used for test purposes only!
       ! In the normal case the analytical formula (getAnalyticalGammaValue) should be used
       if(abs(dist) .le. 1.0e-16_dp) then
          tmp2 = 0.02_dp
          call get_cubic_spline(self%grid, self%lrGamma(Sp1,Sp2,:,1), self%lrGamma(Sp1,Sp2,:,2), tmp2, tmp)
          getGammaValue = tmp * 1.0_dp
       else
          ! NOTE: gammas are generated to 10.0 Bohr,
          ! the interpolation for larger distances is undefined!!!
          if(abs(dist) .le. self%gamma_range) then
             call get_cubic_spline(self%grid, self%lrGamma(Sp1,Sp2,:,1), self%lrGamma(Sp1,Sp2,:,2), dist, tmp)
             getGammaValue = tmp * 1.0_dp
          else
             ! for large distances gamma is equal to (1-exp(-kR))/R
             getGammaValue = (1.0_dp - exp(-self%omega * dist)) / dist
          end if
       end if
    else
       ! Usual case is the analytical gamma formula:
       getGammaValue = getAnalyticalGammaValue(self, Sp1, Sp2, dist)
    end if
    

  end function getGammaValue

  !> Returns the numerical derivative of lr-gamma for iAtom1, iAtom2
  !! \param self, class instance
  !! \param iAtom1,
  !! \param iAtom2,
  !! \param coords,
  !! \param species,
  subroutine getGammaPrimeValue(self, grad, iAtom1, iAtom2, coords, species)
    class(RangeSepFunc), intent(inout) :: self
    real(dp), intent(out) :: grad(3)
    integer, intent(in) :: iAtom1, iAtom2
    real(dp), intent(in) :: coords(:,:)
    integer, intent(in) :: species(:)
    !
    real(dp), parameter :: deltaXDiff = epsilon(1.0_dp)**0.25_dp
    integer :: sp1, sp2, jj, ii
    real(dp) :: vect(3), tmp(3),tmp2(3), dist, deltah
    !
    sp1 = species(iAtom1)
    sp2 = species(iAtom2)
    deltah = deltaXDiff!*0.1_dp
    !    
    do jj = 1, 3 ! x,y,z
       tmp(jj) = 0.0_dp
       do ii = 1, 2 ! +h, -h
          ! difference vector
          vect(:) = coords(:,iAtom2) - coords(:,iAtom1)
          vect(jj) = vect(jj) - real(2 * ii - 3, dp) * deltah
          dist = sqrt(sum(vect(:)**2))
          vect(:) = vect(:) / dist
          !
          tmp(jj) = tmp(jj) + real(2 * ii - 3, dp) * getGammaValue(self, sp1, sp2, dist)
       end do
    end do
    tmp(:) = 0.5_dp * tmp(:) / deltah
    grad = tmp
  end subroutine getGammaPrimeValue

  !> Adds gradients due to long-range HF-contribution.
  !! \param self, class instance
  !! \param gradients 
  !! \param deltaRho, square difference DM (triangle form) 
  !! \param skHamCont
  !! \param skOverCont
  !! \param coords
  !! \param species
  !! \param orb
  !! \param iSquare
  !! \param ovrlapMat
  !! \param iNeighbor
  !! \param nNeighbor
  subroutine addLRGradients(self, gradients, derivator, deltaRho, skHamCont, skOverCont,&
      & coords, species, orb, iSquare, ovrlapMat, iNeighbor, nNeighbor)
    class(RangeSepFunc), intent(inout) :: self
    real(dp), intent(inout) :: gradients(:,:)
    class(NonSccDiff), intent(in) :: derivator
    real(dp), intent(in) :: deltaRho(:,:)
    type(OSlakoCont), intent(in) :: skHamCont, skOverCont
    real(dp), intent(in) :: coords(:,:)
    integer, intent(in) :: species(:)
    type(TOrbitals), intent(in) :: orb
    integer, intent(in) :: iSquare(:)
    real(dp), intent(in) :: ovrlapMat(:,:)
    integer, intent(in) :: iNeighbor(0:,:), nNeighbor(:)

    integer :: nAtom, iAtK, iNeighK, iAtB, iNeighB, iAtC, iAtA, kpa
    real(dp) :: tmpgamma1, tmpgamma2
    real(dp) :: tmpforce(3), tmpforce_r(3), tmpforce2, tmpmultvar1
    integer :: mu, nu, alpha, beta, ccc, kkk
    real(dp)  :: dummy(orb%mOrb,orb%mOrb,3), sPrimeTmp(orb%mOrb,orb%mOrb,3), sPrimeTmp2(orb%mOrb,orb%mOrb,3)
    real(dp), allocatable :: gammaPrimeTmp(:,:,:), tmpovr(:,:), tmpRho(:,:), tmpderiv(:,:)
    real :: start, finish

    @:ASSERT(size(gradients,dim=1) == 3)
    call cpu_time(start)
    call allocateAndInit(tmpovr, tmpRho, gammaPrimeTmp, tmpderiv)
    nAtom = size(self%species)
    tmpderiv = 0.0_dp
    ! sum K
    loopK: do iAtK = 1, nAtom
       ! C >= K
       loopC: do iNeighK = 0, nNeighbor(iAtK)
          iAtC = iNeighbor(iNeighK, iAtK)
          ! evaluate the ovr_prime
          sPrimeTmp2 = 0.0_dp
          sPrimeTmp = 0.0_dp
          if ( iAtK /= iAtC ) then
             call derivator%getFirstDeriv(dummy, skHamCont, coords, &
                 &species, iAtK, iAtC, orb)
             call derivator%getFirstDeriv(sPrimeTmp, skOverCont, coords, &
                 &species, iAtK, iAtC, orb)
             call derivator%getFirstDeriv(dummy, skHamCont, coords, &
                 &species, iAtC, iAtK, orb)
             call derivator%getFirstDeriv(sPrimeTmp2, skOverCont, coords, &
                 &species, iAtC, iAtK, orb)
             sPrimeTmp = 0.5_dp * sPrimeTmp
             sPrimeTmp2 = 0.5_dp * sPrimeTmp2
          end if
          loopB: do iAtB = 1, nAtom
             ! A > B
             loopA: do iNeighB = 0, nNeighbor(iAtB)
                iAtA = iNeighbor(iNeighB, iAtB)
                tmpgamma1 = self%lrGammaEval(iAtK,iAtB) + self%lrGammaEval(iAtC,iAtB)
                tmpgamma2 = tmpgamma1 + self%lrGammaEval(iAtK,iAtA) + self%lrGammaEval(iAtC,iAtA)
                tmpforce = 0.0_dp
                tmpforce_r = 0.0_dp
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

    gradients(:,:) = gradients + tmpderiv * (-0.25_dp)
    
    call cpu_time(finish)
    deallocate(tmpovr, tmpRho, gammaPrimeTmp, tmpderiv)

  contains 

    subroutine allocateAndInit(tmpovr, tmpRho, gammaPrimeTmp, tmpderiv)
      real(dp), allocatable, intent(inout) :: tmpovr(:,:), tmpRho(:,:), gammaPrimeTmp(:,:,:), tmpderiv(:,:)
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

  !> add the long-range HF exchange part to the hamiltonian
  !! \param self, class instance
  !! \param overlap, square overlap matrix (triangle form)
  !! \param deltaRho, square difference DM (triangle form)
  !! \param iSquare, atomic block index array
  !! \param hamiltonian, square hamiltonian (triangle form)
  !! \param species list of all atomic species
  !! \param iNeighbor neighbor list for atoms
  !! \param nNeighbor number of neighbors of each atom
  !! \param orb Information about the orbitals
  subroutine addLRHamiltonian_pbc(self, overlap, deltaRho, iSquare,&
      & hamiltonian, species, iNeighbor, nNeighbor, orb)
    class(RangeSepFunc), intent(inout) :: self
    real(dp), intent(in) :: overlap(:,:), deltaRho(:,:)
    integer, intent(in) :: iSquare(:)
    real(dp), intent(inout) :: hamiltonian(:,:)
    integer, intent(in) :: species(:)
    integer, intent(in) :: iNeighbor(0:,:)
    integer, intent(in) :: nNeighbor(:)
    type(TOrbitals), intent(in) :: orb

    write(*,'(a)') "Rangesep: update the Hamiltonian with PBC"

  end subroutine addLRHamiltonian_pbc
  
  !> evaluate the LR-Energy contribution directly. Very slow, use addLREnergy instead.
  !! \param deltaRho, square density matrix
  !! \param ovrlap, square overlap matrix
  !! \param iSquare
  function evaluateLREnergyDirect(self, deltaRho, ovrlap, iSquare) result (energy)
    class(RangeSepFunc), intent(inout) :: self
    real(dp), intent(in) :: deltaRho(:,:)
    real(dp), intent(in) :: ovrlap(:,:)
    integer, intent(in) :: iSquare(:)

    real(dp) :: energy, tmp
    integer :: iAt1, iAt2, nAtom, mu, nu, alpha, beta
    real(dp), allocatable :: tmpovr(:,:), tmpDRho(:,:)
    real :: start, finish
    !
    call cpu_time(start)
    nAtom = size(self%species)
    allocate(tmpovr(size(ovrlap, dim = 1), size(ovrlap, dim = 1)))
    allocate(tmpDRho(size(deltaRho, dim = 1), size(deltaRho, dim = 1)))
    tmpovr = ovrlap
    tmpDRho = deltaRho
    call symmetrizeSquareMatrix(tmpovr)
    call symmetrizeSquareMatrix(tmpDRho)
    !
    energy = 0.0_dp
    do iAt1 = 1, nAtom
       do iAt2 = 1, nAtom
          tmp = 0.0_dp
          do mu = iSquare(iAt1), iSquare(iAt1 + 1) - 1
             do nu = iSquare(iAt2), iSquare(iAt2 + 1) - 1
                do alpha = 1, size(tmpovr, dim = 1)
                   do beta = 1, size(tmpovr, dim = 1)
                   tmp = tmp + (tmpDRho(alpha,beta) * tmpDRho(mu,nu) + tmpDRho(mu,beta) * tmpDRho(alpha,nu)) &
                        & * tmpovr(mu,alpha) * tmpovr(nu,beta)
                   end do
                end do
             end do
          end do
          energy = energy + tmp * self%lrGammaEval(iAt1,iAt2)
       end do
    end do
    energy = -energy / 8.0_dp

    call cpu_time(finish)
    print '("-> rangeSep.evaluateLREnergyDirect done. time = ", f10.4, " sec.")', finish - start  

  end function evaluateLREnergyDirect


!!!======================================================
!!! Checked until here
!!!======================================================











  !> Returns analytical derivative of lr-gamma for iAtom1, iAtom2
  !! \param self, class instance
  !! \param iAtom1,
  !! \param iAtom2,
  !! \param coords,
  !! \param species,
  subroutine getAnalyticalGammaPrimeValue(self, grad, iAtom1, iAtom2, coords, species)
    class(RangeSepFunc), intent(inout) :: self
    real(dp), intent(out) :: grad(3)
    integer, intent(in) :: iAtom1, iAtom2
    real(dp), intent(in) :: coords(:,:)
    integer, intent(in) :: species(:)
    !
    !ToDO: implement the analytical derivatives, instead of finite differences

    !
  end subroutine getAnalyticalGammaPrimeValue

  !
  subroutine evalEnergy_tst(self,coords,species,orb,iSquare,skOverCont&
      &,iMovedAtom,chgvec,tmpDRho,nrgy)
    ! arguments
    class(RangeSepFunc), intent(inout) :: self
    real(dp), intent(in) :: coords(:,:)
    integer, intent(in) :: species(:)
    type(TOrbitals), intent(in) :: orb
    integer, intent(in) :: iSquare(:)
    type(OSlakoCont), intent(in) :: skOverCont
    integer, intent(in) :: iMovedAtom
    real(dp), intent(in) :: chgvec(3)
    real(dp), intent(in) :: tmpDRho(:,:)
    real(dp), intent(out) :: nrgy
    ! local variables
    real(dp) :: interSKOver(getMIntegrals(skOverCont))  ! interpolated S integs.
    integer :: sp1,sp2,jj,ii,mu,nu,alpha,beta
    integer :: iAtom1, iAtom2, nAtom
    real(dp) :: dist, vect(3), tmp(3), tmp22, lrenergy
    real(dp), allocatable :: coordsmod(:,:), tmpovr(:,:), gammatmp(:,:)
    real(dp) :: Stmp(orb%mOrb,orb%mOrb)
    real :: start, finish
    ! 
    !===============================
    ! get number of atoms
    nAtom = size(orb%nOrbAtom)
    !
    allocate(coordsmod(3,size(coords,dim=2)))
!    write(*,*) "size:",(iSquare(nAtom+1)-1)
    allocate(tmpovr((iSquare(nAtom+1)-1),(iSquare(nAtom+1)-1)))
    allocate(gammatmp(nAtom,nAtom))
    ! prepare the overlap

    !
    coordsmod = coords
    coordsmod(:,iMovedAtom) = coordsmod(:,iMovedAtom) + chgvec(:)
    !
    ! evaluate the overlap matrix

    ! vect(:) = coords(:,atomJ) - coords(:,atomI)
    ! vect(jj) = vect(jj) - real(2*ii-3,dp) * deltaXDiff
    ! dist = sqrt(sum(vect(:)**2))
    ! vect(:) = vect(:) / dist
    ! call getSKIntegrals(skOverCont, interSKOver, dist, sp1, sp2)
    ! call rotateH0(Stmp(:,:,ii,jj), interSKOver, vect(1), vect(2), &
    !     &vect(3), sp1, sp2, orb)
    ! write(*,*) "chgvec)"
    ! write(*,*) chgvec
    ! write(*,*) "Coords"
    ! write(*,*) coords
    ! write(*,*) "Coords_mod"
    ! write(*,*) coordsmod

    tmpovr = 0.0_dp
    do ii = 1, size(tmpovr,dim=1)
       tmpovr(ii,ii) = 1.0_dp
    end do
!    gammatmp = 0.0_dp
    do iAtom1 = 1, nAtom
       do iAtom2 = 1, nAtom
          if(iAtom1 .ne. iAtom2 ) then
 !            write(*,*) iAtom1, iAtom2
             ! normed difference vector
             vect(:) = coordsmod(:,iAtom2) - coordsmod(:,iAtom1)
             dist = sqrt(sum(vect(:)**2))
             vect(:) = vect(:) / dist
             !
             sp1 = species(iAtom1)
             sp2 = species(iAtom2)
             call getSKIntegrals(skOverCont, interSKOver, dist, sp1, sp2)
             call rotateH0(Stmp, interSKOver, vect(1), vect(2), &
              &vect(3), sp1, sp2, orb)
!             write(*,*) "orbAtom:", orb%nOrbAtom(iAtom1),orb%nOrbAtom(iAtom2)
!             write(*,*) Stmp(1,1),Stmp(2,1),Stmp(3,1),Stmp(4,1)
             tmpovr(iSquare(iAtom1):(iSquare(iAtom1+1)-1)&
                 &,iSquare(iAtom2):(iSquare(iAtom2+1)-1)) = transpose(Stmp(1:orb&
                 &%nOrbAtom(iAtom1),1:orb%nOrbAtom(iAtom2)))
             !
!             gammatmp(iAtom1,iAtom2) = getGammaValue(self, sp1, sp2, dist)
             !
          end if
       end do
    end do

    ! write(*,*) "OVERLAP:"
    ! do ii = 1, size(tmpovr,dim=1)
    !    do jj = 1, size(tmpovr,dim=1)
    !       write(*,'(E17.8)',advance='no') tmpovr(ii,jj)
    !    end do
    !    write(*,*)
    ! end do

    ! get gamma

    gammatmp = 0.0_dp
    do iAtom1 = 1, nAtom
       do iAtom2 = 1, nAtom
          !
          vect(:) = coordsmod(:,iAtom2) - coordsmod(:,iAtom1)
          dist = sqrt(sum(vect(:)**2))
          vect(:) = vect(:) / dist
          !
          sp1 = species(iAtom1)
          sp2 = species(iAtom2)
          gammatmp(iAtom1,iAtom2) = getGammaValue(self, sp1, sp2, dist)
       end do
    end do


    ! write(*,*) "SELF_GAMMA:"
    ! do ii = 1, nAtom
    !    do jj = 1, nAtom
    !       write(*,'(E17.8)',advance='no') self%lrGammaEval(ii,jj)
    !    end do
    !    write(*,*)
    ! end do

    ! write(*,*) "MOD_GAMMA:"
    ! do ii = 1, nAtom
    !    do jj = 1, nAtom
    !       write(*,'(E17.8)',advance='no') gammatmp(ii,jj)
    !    end do
    !    write(*,*)
    ! end do

    
!     sp1 = species(iAtom1)
!     sp2 = species(iAtom2)
!     !    
!     do jj = 1, 3 ! x,y,z
!        tmp(jj) = 0.0_dp
!        do ii = 1, 2 ! +h, -h
!           ! difference vector
!           vect(:) = coords(:,iAtom2) - coords(:,iAtom1)
!           vect(jj) = vect(jj) - real(2*ii-3,dp) * deltaXDiff
!           dist = sqrt(sum(vect(:)**2))
!           vect(:) = vect(:) / dist
!           !
!           tmp(jj) = tmp(jj) + real(2*ii-3,dp)*getGammaValue(self, sp1, sp2, dist)
!        end do
!     end do

!     tmp(:) = 0.5_dp * tmp(:) / deltaXDiff

    
    ! write(*,*) "MOD_DENS:"
    ! do ii = 1, nAtom
    !    do jj = 1, nAtom
    !       write(*,'(E17.8)',advance='no') tmpDRho(ii,jj)
    !    end do
    !    write(*,*)
    ! end do

    call cpu_time(start)
    
    ! evaluate the energy with modified overlap and gamma
    tmp22 = 0.0_dp
    lrenergy = 0.0_dp
    do iAtom1 = 1, nAtom
       do iAtom2 = 1, nAtom
          tmp22 = 0.0_dp
          do mu = iSquare(iAtom1),iSquare(iAtom1+1)-1
             do nu = iSquare(iAtom2),iSquare(iAtom2+1)-1
                do alpha = 1, size(tmpovr, dim = 1)
                   do beta = 1, size(tmpovr, dim = 1)
                      tmp22 = tmp22 + (tmpDRho(alpha,beta)*tmpDRho(mu,nu)&
                          &+tmpDRho(mu,beta)*tmpDRho(alpha,nu) )*tmpovr(mu&
                       &,alpha)*tmpovr(nu,beta)
                   end do
                end do
             end do
          end do
          lrenergy = lrenergy + tmp22 * gammatmp(iAtom1,iAtom2)!self%lrGammaEval(iAtom1,iAtom2)
       end do
    end do
    nrgy= -lrenergy/8.0_dp

    call cpu_time(finish)
    print '("-> Done. time = ", f10.4, " sec.")', finish&
        &-start  

    
  end subroutine evalEnergy_tst

  

  !!* Adds gradients due to long-range HF-contribution.
  !!* @param self, class instance
  !!* @param gradients 
  !!* @param deltaRho, square difference DM (triangle form) 
  subroutine addGradients_tst2(self, gradients, derivator, deltaRho, skHamCont, skOverCont,&
      & coords, species,orb,iSquare,ovrlapMat,iNeighbor,nNeighbor)
    class(RangeSepFunc), intent(inout) :: self
    real(dp), intent(inout) :: gradients(:,:)
    class(NonSccDiff), intent(in) :: derivator
    real(dp), intent(in) :: deltaRho(:,:)
    type(OSlakoCont), intent(in) :: skHamCont, skOverCont
    real(dp), intent(in) :: coords(:,:)
    integer, intent(in) :: species(:)
    type(TOrbitals), intent(in) :: orb
    integer, intent(in) :: iSquare(:)
    real(dp), intent(in) :: ovrlapMat(:,:)
    !
    integer, intent(in) :: iNeighbor(0:,:)
    integer, intent(in) :: nNeighbor(:)
    !
    !
!    real(dp) :: Stmp(size(Sprime,dim=1), size(Sprime,dim=2),2,3)
    !
    integer :: nAtom, iAtom1, iAtom2, sp1, sp2, ii, jj, kk
    integer :: iAtK, iNeighK, iAtB, iNeighB, iAtC, iAtA, kpa
    real(dp) :: tmpgamma1, tmpgamma2
    real(dp) :: tmpforce(3), tmpforce_r(3), tmpforce2, tmpmultvar1
    integer :: iAtom3, iAtom4, mu, nu, alpha, beta, iMovedAtom
    ! NOTE: dummy is not needed
    real(dp)  :: dummy(orb%mOrb,orb%mOrb,3), sPrimeTmp(orb%mOrb,orb%mOrb,3), sPrimeTmp2(orb%mOrb,orb%mOrb,3)
    real(dp) :: vect(3), tmp(3), lrenergy, tmp33, lrenergy2, tmp55
    real(dp) :: dist, rr, tmp1, tmp2(3), tmp3(3), lrgammatmp(3), tmpomat(3)
    real(dp), allocatable :: gammaPrimeTmp(:,:,:), ovrPrimeTmp(:,:,:),&
        & tmpovr(:,:), tmpRho(:,:), tmpderiv(:,:), tmpderiv2(:,:)
    ! arrays for screening
    real(dp), allocatable :: tstdensmat(:,:), tstovr(:,:), tstovrprime(:,:,:)
    real(dp) :: tsttmp_1, tsttmp_2, scrcond1, scrcond2, scrcondtst
    ! for cpu_start/finish routines
    real :: start, finish, start2, finish2
    !
    real(dp), allocatable :: sigmatst(:,:,:)
    integer :: ccc, kkk
    !
    write(*,*) "rangesep.addGradients_tst"
    !write(*,*) "deltaXDiff=", deltaXDiff
    !@:ASSERT(self%tInit)
    ! 
    @:ASSERT(size(gradients,dim=1) == 3)
    !

    !========================================== 
    ! orb%nOrbAtom(ii) => number of orbitals of atom ii
    nAtom = size(orb%nOrbAtom)
    ! deriv(:,:) = 0.0_dp

    !========================
    ! Brute force prototype
    !========================

    ! arrays for screening
!    allocate(tstdensmat(nAtom,nAtom))
!    allocate(tstovr(nAtom,nAtom))
!    allocate(tstovrprime(nAtom,nAtom,3))
    !
    allocate(tmpovr(size(ovrlapMat, dim=1),size(ovrlapMat, dim=1)))
    allocate(tmpRho(size(deltaRho, dim=1),size(deltaRho, dim=1)))

    allocate(gammaPrimeTmp(3,nAtom,nAtom))
!    allocate(ovrPrimeTmp(size(deltaRho,dim=1),size(deltaRho,dim=1),3))
    allocate(tmpderiv(3,size(gradients,dim=2)))
!    allocate(tmpderiv2(3,size(gradients,dim=2)))

    ! precompute values overlap times densmat
 !   allocate(sigmatst(size(ovrlapMat, dim=1),size(ovrlapMat, dim=1),nAtom))
    !

!    write(*,*) "OVERLAP:"
    do ii=1,size(ovrlapMat, dim=1)
       do jj=1,size(ovrlapMat, dim=1)
          if(jj .gt. ii) then
             tmpovr(ii,jj)=ovrlapMat(jj,ii)
          else
             tmpovr(ii,jj)=ovrlapMat(ii,jj)
          end if
!          write(*,'(E17.8)', advance='no') tmpovr(ii,jj)!overlap(ii,jj)
       end do
!              write(*,*)
    end do
!       write(*,*) "DELTAP "
    do ii=1,size(deltaRho, dim=1)
       do jj=1,size(deltaRho, dim=1)
          if(jj .gt. ii) then
             tmpRho(ii,jj)=deltaRho(jj,ii)
          else
             tmpRho(ii,jj)=deltaRho(ii,jj)
          end if
 !                  write(*,'(E17.8)', advance='no') tmpRho(ii,jj)!deltaRho(ii,jj)
       end do
 !            write(*,*)
    end do
    
    ! precompute the gamma derivatives
 !   ovrPrimeTmp = 0.0_dp
    gammaPrimeTmp = 0.0_dp
    do iAtom1 = 1, nAtom
       do iAtom2 = 1, nAtom
          if(iAtom1 .ne. iAtom2) then
             call getGammaPrimeValue(self, tmp, iAtom1, iAtom2, coords, species)
             gammaPrimeTmp(:,iAtom1,iAtom2) = tmp(:)
             !
             ! call H0Sprime(dummy, sPrimeTmp, skHamCont, skOverCont, coords, &
             !     &species, iAtom1, iAtom2, orb)
             ! do ii =1, 3
             !    ovrPrimeTmp(iSquare(iAtom1):(iSquare(iAtom1+1)-1)&
             !        &,iSquare(iAtom2):(iSquare(iAtom2+1)-1),ii) = transpose(sPrimeTmp(1:orb&
             !        &%nOrbAtom(iAtom1),1:orb%nOrbAtom(iAtom2),ii))
             ! end do
          end if
       end do
    end do
!    ovrPrimeTmp = ovrPrimeTmp * 0.5_dp
    !

    !
    ! do iAtom1=1, nAtom
    !    do iAtom2=1,nAtom
    !       ! overlap
    !       tsttmp_1=maxval(abs(tmpovr((iSquare(iAtom1)):(iSquare(iAtom1+1)-1)&
    !           &,(iSquare(iAtom2)):(iSquare(iAtom2+1)-1))))
    !       tstovr(iAtom1,iAtom2)=abs(tsttmp_1)
    !       ! overlap prime
    !       tsttmp_1=maxval(abs(ovrPrimeTmp((iSquare(iAtom1)):(iSquare(iAtom1+1)-1)&
    !           &,(iSquare(iAtom2)):(iSquare(iAtom2+1)-1),1)))
    !       tstovrprime(iAtom1,iAtom2,1)=abs(tsttmp_1)
    !       tsttmp_1=maxval(abs(ovrPrimeTmp((iSquare(iAtom1)):(iSquare(iAtom1+1)-1)&
    !           &,(iSquare(iAtom2)):(iSquare(iAtom2+1)-1),2)))
    !       tstovrprime(iAtom1,iAtom2,2)=abs(tsttmp_1)
    !       tsttmp_1=maxval(abs(ovrPrimeTmp((iSquare(iAtom1)):(iSquare(iAtom1+1)-1)&
    !           &,(iSquare(iAtom2)):(iSquare(iAtom2+1)-1),3)))
    !       tstovrprime(iAtom1,iAtom2,3)=abs(tsttmp_1)
    !       ! densmat 
    !       tsttmp_1=maxval(abs(tmpRho((iSquare(iAtom1)):(iSquare(iAtom1+1)-1)&
    !           &,(iSquare(iAtom2)):(iSquare(iAtom2+1)-1))))
    !       tstdensmat(iAtom1,iAtom2)=abs(tsttmp_1)
    !    end do
    ! end do
    !==============================
    ! NEIGHBOR LIST BASED
    !==============================
    call cpu_time(start2)
    tmpderiv = 0.0_dp
    ! sum K
    do iAtK = 1, nAtom
       ! C >= K
       do iNeighK = 0, nNeighbor(iAtK)
          iAtC = iNeighbor(iNeighK, iAtK)
          !
          ! evaluate the ovr_prime
          sPrimeTmp2 = 0.0_dp
          sPrimeTmp = 0.0_dp
          if ( iAtK /= iAtC ) then
             call derivator%getFirstDeriv(dummy, skHamCont, coords, &
                 &species, iAtK, iAtC, orb)
             call derivator%getFirstDeriv(sPrimeTmp, skOverCont, coords, &
                 &species, iAtK, iAtC, orb)
             call derivator%getFirstDeriv(dummy, skHamCont, coords, &
                 &species, iAtC, iAtK, orb)
             call derivator%getFirstDeriv(sPrimeTmp2, skOverCont, coords, &
                 &species, iAtC, iAtK, orb)
             sPrimeTmp = 0.5_dp * sPrimeTmp
             sPrimeTmp2 = 0.5_dp * sPrimeTmp2
          end if

          do iAtB = 1, nAtom
             ! A > B
             do iNeighB = 0, nNeighbor(iAtB)
                iAtA = iNeighbor(iNeighB, iAtB)
!             do iAtA = 1, nAtom
                ! mu, kpa, alpha, beta loop
                tmpgamma1 = self%lrGammaEval(iAtK,iAtB) + self&
                    &%lrGammaEval(iAtC,iAtB)
                tmpgamma2 = tmpgamma1 + self%lrGammaEval(iAtK,iAtA) + self&
                    &%lrGammaEval(iAtC,iAtA)
                tmpforce = 0.0_dp
                tmpforce_r = 0.0_dp
                tmpforce2 = 0.0_dp


                ccc = 0
                do mu = iSquare(iAtC), iSquare(iAtC+1)-1
                   ccc = ccc + 1
                   kkk = 0
                   do kpa = iSquare(iAtK), iSquare(iAtK+1)-1
                      kkk = kkk + 1 
                      tmpmultvar1 = 0.0_dp
                      do alpha = iSquare(iAtA), iSquare(iAtA+1)-1
                         do beta = iSquare(iAtB), iSquare(iAtB+1)-1
                            !
!                            tmpmultvar1 = tmpovr(alpha,beta)*(tmpRho(kpa,beta)&
!                                &*tmpRho(mu,alpha)+tmpRho(kpa,alpha)*tmpRho(mu&
!                                &,beta))
!                            tmpmultvar1 = tmpmultvar1 + tmpovr(alpha,beta)*(tmpRho(kpa,beta)&
!                                &*tmpRho(mu,alpha)+tmpRho(kpa,alpha)*tmpRho(mu&
!                                &,beta))
                            tmpmultvar1 = tmpmultvar1 + tmpovr(beta,alpha)&
                                &*(tmpRho(beta,kpa)&
                                &*tmpRho(alpha,mu)+tmpRho(alpha,kpa)&
                                &*tmpRho(beta,mu))
                            !                            
                            
                         end do
                      end do
                      
                      tmpforce(:) = tmpforce(:) + tmpmultvar1&
                          &*(sPrimeTmp(ccc,kkk,:))
                      tmpforce_r(:) = tmpforce_r(:) + tmpmultvar1&
                          &*(sPrimeTmp2(kkk,ccc,:))

                      
!                       tmpforce(:) = tmpforce(:) + tmpmultvar1*ovrPrimeTmp(kpa,mu,:)
!                       tmpforce_r(:) = tmpforce_r(:) + tmpmultvar1&
!                           &*ovrPrimeTmp(mu,kpa,:)
                      !
                      tmpforce2 = tmpforce2 + tmpmultvar1*tmpovr(kpa,mu)
                      
                   end do
                end do
                !
                ! C /= K
                if( iAtK /= iAtC ) then
                   if( iAtB /= iAtA) then
                      tmpforce(:) = tmpforce(:)*tmpgamma2
                      tmpforce_r(:) = tmpforce_r(:)*tmpgamma2
                      !

                      tmpforce(:) = tmpforce(:) + tmpforce2 * (gammaPrimeTmp(:&
                          &,iAtK,iAtA) + gammaPrimeTmp(:,iAtK,iAtB))
                      tmpforce_r(:) = tmpforce_r(:) + tmpforce2 * (gammaPrimeTmp(:&
                          &,iAtC,iAtA) + gammaPrimeTmp(:,iAtC,iAtB))
                   else
                       tmpforce(:) = tmpforce(:)*tmpgamma1
                       tmpforce_r(:) = tmpforce_r(:)*tmpgamma1

                      tmpforce(:) = tmpforce(:) + tmpforce2 * (gammaPrimeTmp(:&
                          &,iAtK,iAtA))
                      tmpforce_r(:) = tmpforce_r(:) + tmpforce2 * (gammaPrimeTmp(:&
                          &,iAtC,iAtA))
                      
                   end if
                else
                   !
                   if( iAtB /= iAtA) then
                      tmpforce(:) = tmpforce(:) + tmpforce2 * (gammaPrimeTmp(:&
                          &,iAtK,iAtA) + gammaPrimeTmp(:,iAtK,iAtB))
                   else
                       tmpforce(:) = tmpforce(:) + tmpforce2 * (gammaPrimeTmp(:&
                           &,iAtK,iAtA))
                   end if
                end if
                
                tmpderiv(:,iAtK) = tmpderiv(:,iAtK) + tmpforce(:)
                tmpderiv(:,iAtC) = tmpderiv(:,iAtC) + tmpforce_r(:)

             end do
          end do
       end do
    end do

    tmpderiv(:,:) = tmpderiv*(-0.25_dp)
    gradients(:,:)=gradients+tmpderiv
    
    call cpu_time(finish2)
 



!      write(*,*) "LR-FORCES 2:"
!      write(*,*) tmpderiv(:,:)
!     !==============================

    
!     call cpu_time(start)
!     !    
!     tmpderiv = 0.0_dp
!     tmp1 = 0.0_dp
!     do iMovedAtom = 1, nAtom
!        write(*,*) "Moving Atom", iMovedAtom
!              call cpu_time(start2)

!              do iAtom3 = 1, nAtom
!                 if(iAtom3 .ne. iMovedAtom) then
                   
!                    ! A=1, B=2, C=3, K=Moved
!                    !
!                    do iAtom1 = 1, nAtom
!                       do iAtom2 = 1, nAtom
!                          !
!                          scrcond1 = tstovr(iMovedAtom,iAtom2)*tstovr(iAtom1,iAtom3)&
!                              &*(tstdensmat(iMovedAtom,iAtom3)*tstdensmat(iAtom2&
!                              &,iAtom1)+tstdensmat(iMovedAtom,iAtom1)&
!                              &*tstdensmat(iAtom2,iAtom3))
!                          scrcond1 = 100.0_dp * maxval(abs(gammaPrimeTmp(:,iMovedAtom&
!                              &,iAtom3))) * scrcond1
!                          !
!                          scrcond2 = tstovr(iAtom1,iAtom2) *&
!                              & maxval(abs(tstovrprime(iMovedAtom,iAtom3,:)))&
!                              &*(tstdensmat(iMovedAtom,iAtom2)*tstdensmat(iAtom3&
!                              &,iAtom1)+tstdensmat(iMovedAtom,iAtom1)&
!                              &*tstdensmat(iAtom3,iAtom2))
!                          scrcond2 = scrcond2 * 100.0_dp *(self&
!                              &%lrGammaEval(iMovedAtom,iAtom2)+ self%lrGammaEval(iAtom3,iAtom2))

!                         scrcondtst = max(scrcond2,scrcond1)
!                         !
!  !                       if(iAtom1 /= iAtom2 .and. iAtom3 /= iMovedAtom) then
!                         if(scrcond1 .gt. 1.0e-16_dp) then
!                            !

!                            do beta=iSquare(iAtom2),iSquare(iAtom2+1)-1
!                               do nu=iSquare(iMovedAtom),iSquare(iMovedAtom+1)-1
!                                  do alpha=iSquare(iAtom1),iSquare(iAtom1+1)-1
!                                     do mu=iSquare(iAtom3),iSquare(iAtom3+1)-1
                                       
!                                        tmp1 = tmpRho(alpha,beta)*tmpRho(mu&
!                                            &,nu) +&
!                                            & tmpRho(mu,beta)*tmpRho(alpha,nu)
                                       
!                                        tmpderiv(:,iMovedAtom)=tmpderiv(:,iMovedAtom)+tmp1&
!                                            &*tmpovr(mu,alpha)*tmpovr(nu,beta)&
!                                            &*gammaPrimeTmp(:,iMovedAtom&
!                                            &,iAtom3)
                                       
!                                     end do
!                                  end do
!                               end do
!                            end do
!                         end if
!   !                   end if
!                             !
! !                            if(iAtom1 /= iAtom2) then
!                         if(scrcond2 .gt. 1.0e-16_dp) then
!                             !
!                             do beta=iSquare(iAtom2),iSquare(iAtom2+1)-1
!                                do nu=iSquare(iMovedAtom),iSquare(iMovedAtom+1)-1
!                                   do alpha=iSquare(iAtom1),iSquare(iAtom1+1)-1
!                                      do mu=iSquare(iAtom3),iSquare(iAtom3+1)-1

!                                         tmp1 = tmpRho(nu,beta)*tmpRho(mu,alpha)+tmpRho(mu&
!                                             &,beta)*tmpRho(nu,alpha)
!                                         !
!                                         tmpderiv(:,iMovedAtom) = tmpderiv(:,iMovedAtom) +&
!                                             & tmp1*(self&
!                                             &%lrGammaEval(iMovedAtom,iAtom2)&
!                                             &+self%lrGammaEval(iAtom3,iAtom2))&
!                                             &*tmpovr(alpha,beta)&
!                                             &*ovrPrimeTmp(nu,mu,:)
!                                      end do
!                                   end do
!                                end do
!                             end do
!  !                        end if
                            
!                          end if
                            
!                        end do
!                    end do
                   
!                 end if
!                 !
!              end do
!              call cpu_time(finish2)
! !             write(*,'(a,I5)') "Atom moved", iMovedAtom
! !             print '("-> done. time = ", f10.4, " sec.")', finish2&
! !                 &-start2  
             
!           end do
          
!           tmpderiv2 = -tmpderiv/4.0_dp
!     ! update the gradients
! !    gradients=gradients+tmpderiv2


!  !   write(*,*) "LR-FORCES 2:"
!  !   write(*,*) tmpderiv2(:,:)

    
!    call cpu_time(finish)
!    print '("-> addGradients_tst done. time = ", f10.4, " sec.")', finish&
!        &-start  
 
    deallocate(gammaPrimeTmp)
!    deallocate(tmpovr,tmpRho,ovrPrimeTmp,tmpderiv)
    deallocate(tmpovr,tmpRho,tmpderiv)
  end subroutine addGradients_tst2



  !!* Adds gradients due to long-range HF-contribution.
  !!* @param self, class instance
  !!* @param gradients 
  !!* @param deltaRho, square difference DM (triangle form) 
  subroutine addGradients_tst(self, gradients, derivator, deltaRho, skHamCont, skOverCont,&
      & coords, species,orb,iSquare,ovrlapMat,iNeighbor,nNeighbor)
    class(RangeSepFunc), intent(inout) :: self
    real(dp), intent(inout) :: gradients(:,:)
    class(NonSccDiff), intent(in) :: derivator
    real(dp), intent(in) :: deltaRho(:,:)
    type(OSlakoCont), intent(in) :: skHamCont, skOverCont
    real(dp), intent(in) :: coords(:,:)
    integer, intent(in) :: species(:)
    type(TOrbitals), intent(in) :: orb
    integer, intent(in) :: iSquare(:)
    real(dp), intent(in) :: ovrlapMat(:,:)
    !
    integer, intent(in) :: iNeighbor(0:,:)
    integer, intent(in) :: nNeighbor(:)
    !
    !
!    real(dp) :: Stmp(size(Sprime,dim=1), size(Sprime,dim=2),2,3)
    !
    integer :: nAtom, iAtom1, iAtom2, sp1, sp2, ii, jj, kk
    integer :: iAtK, iNeighK, iAtB, iNeighB, iAtC, iAtA, kpa
    real(dp) :: tmpgamma1, tmpgamma2
    real(dp) :: tmpforce(3), tmpforce_r(3), tmpforce2, tmpmultvar1
    integer :: iAtom3, iAtom4, mu, nu, alpha, beta, iMovedAtom
    ! NOTE: dummy is not needed
    real(dp)  :: dummy(orb%mOrb,orb%mOrb,3), sPrimeTmp(orb%mOrb,orb%mOrb,3)
    real(dp) :: vect(3), tmp(3), lrenergy, tmp33, lrenergy2, tmp55
    real(dp) :: dist, rr, tmp1, tmp2(3), tmp3(3), lrgammatmp(3), tmpomat(3)
    real(dp), allocatable :: gammaPrimeTmp(:,:,:), ovrPrimeTmp(:,:,:),&
        & tmpovr(:,:), tmpRho(:,:), tmpderiv(:,:), tmpderiv2(:,:)
    ! arrays for screening
    real(dp), allocatable :: tstdensmat(:,:), tstovr(:,:), tstovrprime(:,:,:)
    real(dp) :: tsttmp_1, tsttmp_2, scrcond1, scrcond2, scrcondtst
    ! for cpu_start/finish routines
    real :: start, finish, start2, finish2
    !
    real(dp), allocatable :: sigmatst(:,:,:)
    !
    write(*,*) "rangesep.addGradients_tst"
    !write(*,*) "deltaXDiff=", deltaXDiff
    !@:ASSERT(self%tInit)
    ! 
    @:ASSERT(size(gradients,dim=1) == 3)
    !

    !========================================== 
    ! orb%nOrbAtom(ii) => number of orbitals of atom ii
    nAtom = size(orb%nOrbAtom)
    ! deriv(:,:) = 0.0_dp

    !========================
    ! Brute force prototype
    !========================

    ! arrays for screening
    allocate(tstdensmat(nAtom,nAtom))
    allocate(tstovr(nAtom,nAtom))
    allocate(tstovrprime(nAtom,nAtom,3))
    !
    allocate(tmpovr(size(ovrlapMat, dim=1),size(ovrlapMat, dim=1)))
    allocate(tmpRho(size(deltaRho, dim=1),size(deltaRho, dim=1)))

    allocate(gammaPrimeTmp(3,nAtom,nAtom))
    allocate(ovrPrimeTmp(size(deltaRho,dim=1),size(deltaRho,dim=1),3))
    allocate(tmpderiv(3,size(gradients,dim=2)))
    allocate(tmpderiv2(3,size(gradients,dim=2)))

    ! precompute values overlap times densmat
    allocate(sigmatst(size(ovrlapMat, dim=1),size(ovrlapMat, dim=1),nAtom))
    !

!    write(*,*) "OVERLAP:"
    do ii=1,size(ovrlapMat, dim=1)
       do jj=1,size(ovrlapMat, dim=1)
          if(jj .gt. ii) then
             tmpovr(ii,jj)=ovrlapMat(jj,ii)
          else
             tmpovr(ii,jj)=ovrlapMat(ii,jj)
          end if
!          write(*,'(E17.8)', advance='no') tmpovr(ii,jj)!overlap(ii,jj)
       end do
!              write(*,*)
    end do
!       write(*,*) "DELTAP "
    do ii=1,size(deltaRho, dim=1)
       do jj=1,size(deltaRho, dim=1)
          if(jj .gt. ii) then
             tmpRho(ii,jj)=deltaRho(jj,ii)
          else
             tmpRho(ii,jj)=deltaRho(ii,jj)
          end if
 !                  write(*,'(E17.8)', advance='no') tmpRho(ii,jj)!deltaRho(ii,jj)
       end do
 !            write(*,*)
    end do
    
    ! precompute the gamma derivatives
    write(*,*) "precomputing the lr-gamma derivatives"
    ovrPrimeTmp = 0.0_dp
    gammaPrimeTmp = 0.0_dp
    do iAtom1 = 1, nAtom
       do iAtom2 = 1, nAtom
          if(iAtom1 .ne. iAtom2) then
             call getGammaPrimeValue(self, tmp, iAtom1, iAtom2, coords, species)
             gammaPrimeTmp(:,iAtom1,iAtom2) = tmp(:)
             !
             call derivator%getFirstDeriv(dummy, skHamCont, coords, &
                 &species, iAtom1, iAtom2, orb)
             call derivator%getFirstDeriv(sPrimeTmp, skOverCont, coords, &
                 &species, iAtom1, iAtom2, orb)
             do ii =1, 3
                ovrPrimeTmp(iSquare(iAtom1):(iSquare(iAtom1+1)-1)&
                    &,iSquare(iAtom2):(iSquare(iAtom2+1)-1),ii) = transpose(sPrimeTmp(1:orb&
                    &%nOrbAtom(iAtom1),1:orb%nOrbAtom(iAtom2),ii))
             end do
          end if
       end do
    end do
    ovrPrimeTmp = ovrPrimeTmp * 0.5_dp
    !

    !
    do iAtom1=1, nAtom
       do iAtom2=1,nAtom
          ! overlap
          tsttmp_1=maxval(abs(tmpovr((iSquare(iAtom1)):(iSquare(iAtom1+1)-1)&
              &,(iSquare(iAtom2)):(iSquare(iAtom2+1)-1))))
          tstovr(iAtom1,iAtom2)=abs(tsttmp_1)
          ! overlap prime
          tsttmp_1=maxval(abs(ovrPrimeTmp((iSquare(iAtom1)):(iSquare(iAtom1+1)-1)&
              &,(iSquare(iAtom2)):(iSquare(iAtom2+1)-1),1)))
          tstovrprime(iAtom1,iAtom2,1)=abs(tsttmp_1)
          tsttmp_1=maxval(abs(ovrPrimeTmp((iSquare(iAtom1)):(iSquare(iAtom1+1)-1)&
              &,(iSquare(iAtom2)):(iSquare(iAtom2+1)-1),2)))
          tstovrprime(iAtom1,iAtom2,2)=abs(tsttmp_1)
          tsttmp_1=maxval(abs(ovrPrimeTmp((iSquare(iAtom1)):(iSquare(iAtom1+1)-1)&
              &,(iSquare(iAtom2)):(iSquare(iAtom2+1)-1),3)))
          tstovrprime(iAtom1,iAtom2,3)=abs(tsttmp_1)
          ! densmat 
          tsttmp_1=maxval(abs(tmpRho((iSquare(iAtom1)):(iSquare(iAtom1+1)-1)&
              &,(iSquare(iAtom2)):(iSquare(iAtom2+1)-1))))
          tstdensmat(iAtom1,iAtom2)=abs(tsttmp_1)
       end do
    end do

    call cpu_time(start2)
    tmpderiv = 0.0_dp
    ! sum K
    do iAtK = 1, nAtom
       ! C > K
       do iNeighK = 1, nNeighbor(iAtK)
          iAtC = iNeighbor(iNeighK, iAtK)
          !
          do iAtB = 1, nAtom
             ! A > B
             do iNeighB = 0, nNeighbor(iAtB)
                iAtA = iNeighbor(iNeighB, iAtB)
!             do iAtA = 1, nAtom
                ! mu, kpa, alpha, beta loop
                tmpgamma1 = self%lrGammaEval(iAtK,iAtB) + self&
                    &%lrGammaEval(iAtC,iAtB)
                tmpgamma2 = tmpgamma1 + self%lrGammaEval(iAtK,iAtA) + self&
                    &%lrGammaEval(iAtC,iAtA)
                tmpforce = 0.0_dp
                tmpforce_r = 0.0_dp
                tmpforce2 = 0.0_dp
!                write(*,'(4I4)') iAtK, iAtC, iAtA, iAtB
                do mu = iSquare(iAtC), iSquare(iAtC+1)-1
                   do kpa = iSquare(iAtK), iSquare(iAtK+1)-1
                      do alpha = iSquare(iAtA), iSquare(iAtA+1)-1
                         do beta = iSquare(iAtB), iSquare(iAtB+1)-1
                            !
                            tmpmultvar1 = tmpovr(alpha,beta)*(tmpRho(kpa,beta)&
                                &*tmpRho(mu,alpha)+tmpRho(kpa,alpha)*tmpRho(mu&
                                &,beta))
                            
!                            tmpforce(:) = tmpforce(:) + tmpmultvar1*(ovrPrimeTmp(kpa,mu,:))*tmpgamma2
!                            tmpforce_r(:) = tmpforce_r(:) + tmpmultvar1*(ovrPrimeTmp(mu,kpa,:))*tmpgamma2
                            if(iAtB /= iAtA ) then
                               tmpforce(:) = tmpforce(:) + tmpmultvar1*(ovrPrimeTmp(kpa,mu,:))*tmpgamma2
                               tmpforce_r(:) = tmpforce_r(:) + tmpmultvar1&
                                   &*(ovrPrimeTmp(mu,kpa,:))*tmpgamma2
                            else
                               tmpforce(:) = tmpforce(:) + tmpmultvar1*(ovrPrimeTmp(kpa,mu,:))*tmpgamma1
                               tmpforce_r(:) = tmpforce_r(:) + tmpmultvar1&
                                   &*(ovrPrimeTmp(mu,kpa,:))*tmpgamma1
                            end if
                            tmpforce2 = tmpforce2 + tmpovr(kpa,alpha)&
                                &*tmpovr(beta,mu)*(tmpRho(kpa,mu)*tmpRho(alpha&
                                &,beta)+tmpRho(kpa,beta)*tmpRho(alpha,mu))
                            
                            tmpforce2 = tmpforce2 + tmpovr(kpa,beta)&
                                &*tmpovr(alpha,mu)*(tmpRho(kpa,mu)*tmpRho(beta&
                                &,alpha)+tmpRho(kpa,alpha)*tmpRho(beta,mu)) 
                            
                         end do
                      end do
                   end do
                end do
                write(*,'(a,4I3,a,3F17.8)') "###:", iAtK,iAtC,iAtA,iAtB,":",&
                    & tmpforce(:)
                !
                ! write(*,*) "tmpgamma1=", tmpgamma1
                ! write(*,*) "gammaPrimeTmp", gammaPrimeTmp(:,iAtK,iAtC)
                ! write(*,*) "gammaPrimeTmp", gammaPrimeTmp(:,iAtC,iAtK)
                ! write(*,*) "tmpforce(:)", tmpforce(:)
                ! write(*,*) "tmpforce2", tmpforce2
                !
                tmpderiv(:,iAtK) = tmpderiv(:,iAtK) + tmpforce(:) !+ tmpforce2&
!                    &*(gammaPrimeTmp(:,iAtK,iAtC) )
                !
                tmpderiv(:,iAtC) = tmpderiv(:,iAtC) + tmpforce_r(:) !+ tmpforce2&
!                    &*(gammaPrimeTmp(:,iAtC,iAtK))
             end do
!              ! A = B
!               tmpforce = 0.0_dp
!               tmpforce_r = 0.0_dp
!               tmpforce2 = 0.0_dp
!               do mu = iSquare(iAtC), iSquare(iAtC+1)-1
!                  do kpa = iSquare(iAtK), iSquare(iAtK+1)-1
!                     do alpha = iSquare(iAtB), iSquare(iAtB+1)-1
!                        do beta = iSquare(iAtB), iSquare(iAtB+1)-1
! ! !                      alpha = beta
! !                          !
! !                          ! tmpforce(:) = tmpforce(:) + tmpRho(kpa,beta)&
! !                          !     &*tmpRho(mu,beta)*ovrPrimeTmp(kpa,mu,:)
! !                          ! tmpforce_r(:) = tmpforce_r(:) + tmpRho(kpa,beta)&
! !                          !     &*tmpRho(mu,beta)*ovrPrimeTmp(mu,kpa,:)
                         
! !                          ! tmpforce2 = tmpforce2 + tmpovr(kpa,beta)*tmpovr(beta&
! !                          !     &,mu)*(tmpRho(kpa,mu)*tmpRho(beta,beta) +&
! !                          !     & tmpRho(kpa,beta)*tmpRho(beta,mu))
!                           tmpforce(:) = tmpforce(:) + tmpovr(alpha,beta)*(tmpRho(kpa,beta)&
!                               &*tmpRho(mu,alpha) + tmpRho(kpa,alpha)*tmpRho(mu&
!                               &,beta))*ovrPrimeTmp(kpa,mu,:)
                         
!                           tmpforce_r(:) = tmpforce_r(:) + tmpovr(alpha,beta)*(tmpRho(kpa,beta)&
!                               &*tmpRho(mu,alpha) + tmpRho(kpa,alpha)*tmpRho(mu&
!                               &,beta))*ovrPrimeTmp(mu,kpa,:)
                         
! !                          ! tmpforce2 = tmpforce2 + tmpovr(kpa,beta)*tmpovr(beta&
! !                          !     &,mu)*(tmpRho(kpa,mu)*tmpRho(beta,beta) +&
! !                          !     & tmpRho(kpa,beta)*tmpRho(beta,mu))
!                        end do
!                     end do
!                  end do
!               end do
!                tmpderiv(:,iAtK) = tmpderiv(:,iAtK) + tmpforce(:)*tmpgamma1
!                tmpderiv(:,iAtC) = tmpderiv(:,iAtC) + tmpforce_r(:)*tmpgamma1
             ! !
!              tmpderiv(:,iAtK) = tmpderiv(:,iAtK) + tmpforce2*(gammaPrimeTmp(:&
!                  &,iAtK,iAtC) )
             ! !
!              tmpderiv(:,iAtC) = tmpderiv(:,iAtC) + tmpforce2*(gammaPrimeTmp(:&
!                  &,iAtC,iAtK) )
             !
          end do
       end do
    end do

!     call cpu_time(start2)
!     tmpderiv = 0.0_dp
!     ! sum K
!     do iAtK = 1, nAtom
!        ! C > K
!        do iNeighK = 1, nNeighbor(iAtK)
!           iAtC = iNeighbor(iNeighK, iAtK)
!           !
!           do iAtB = 1, nAtom
!              ! A > B
!              do iNeighB = 1, nNeighbor(iAtB)
!                 iAtA = iNeighbor(iNeighB, iAtB)
!                 ! mu, kpa, alpha, beta loop
!                 tmpgamma1 = self%lrGammaEval(iAtK,iAtB) + self&
!                     &%lrGammaEval(iAtC,iAtB)
!                 tmpgamma2 = tmpgamma1 + self%lrGammaEval(iAtK,iAtA) + self&
!                     &%lrGammaEval(iAtC,iAtA)
!                 tmpforce = 0.0_dp
!                 tmpforce_r = 0.0_dp
!                 tmpforce2 = 0.0_dp
!                 write(*,'(4I4)') iAtK, iAtC, iAtA, iAtB
!                 do mu = iSquare(iAtC), iSquare(iAtC+1)-1
!                    do kpa = iSquare(iAtK), iSquare(iAtK+1)-1
!                       do alpha = iSquare(iAtA), iSquare(iAtA+1)-1
!                          do beta = iSquare(iAtB), iSquare(iAtB+1)-1
!                             !
!                             tmpforce(:) = tmpforce(:) + tmpovr(alpha,beta)&
!                                 &*(tmpRho(kpa,beta)*tmpRho(mu,alpha)&
!                                 &+tmpRho(kpa,alpha)*tmpRho(mu,beta))&
!                                 &*(ovrPrimeTmp(kpa,mu,:))*tmpgamma2
                            
!                             tmpforce_r(:) = tmpforce_r(:) + tmpovr(alpha,beta)&
!                                 &*(tmpRho(kpa,beta)*tmpRho(mu,alpha)&
!                                 &+tmpRho(kpa,alpha)*tmpRho(mu,beta))&
!                                 &*(ovrPrimeTmp(mu,kpa,:))*tmpgamma2
!                             !
                            
!                             ! tmpforce2 = tmpforce2 + tmpovr(kpa,alpha)&
!                             !     &*tmpovr(beta,mu)*(tmpRho(kpa,mu)*tmpRho(alpha&
!                             !     &,beta)+tmpRho(kpa,beta)*tmpRho(alpha,mu))
                            
!                             ! tmpforce2 = tmpforce2 + tmpovr(kpa,beta)&
!                             !     &*tmpovr(alpha,mu)*(tmpRho(kpa,mu)*tmpRho(beta&
!                             !     &,alpha)+tmpRho(kpa,alpha)*tmpRho(beta,mu)) 
                            
!                          end do
!                       end do
!                    end do
!                 end do
!                 !
!                 ! write(*,*) "tmpgamma1=", tmpgamma1
!                 ! write(*,*) "gammaPrimeTmp", gammaPrimeTmp(:,iAtK,iAtC)
!                 ! write(*,*) "gammaPrimeTmp", gammaPrimeTmp(:,iAtC,iAtK)
!                 ! write(*,*) "tmpforce(:)", tmpforce(:)
!                 ! write(*,*) "tmpforce2", tmpforce2
!                 !
!                 tmpderiv(:,iAtK) = tmpderiv(:,iAtK) + tmpforce(:) !+ tmpforce2&
! !                    &*(gammaPrimeTmp(:,iAtK,iAtC) )
!                 !
!                 tmpderiv(:,iAtC) = tmpderiv(:,iAtC) + tmpforce_r(:) !+ tmpforce2&
! !                    &*(gammaPrimeTmp(:,iAtC,iAtK))
!              end do
!              ! A = B
!              tmpforce = 0.0_dp
!              tmpforce_r = 0.0_dp
!              tmpforce2 = 0.0_dp
!              do mu = iSquare(iAtC), iSquare(iAtC+1)-1
!                 do kpa = iSquare(iAtK), iSquare(iAtK+1)-1
!                    do alpha = iSquare(iAtB), iSquare(iAtB+1)-1
!                       do beta = iSquare(iAtB), iSquare(iAtB+1)-1
! !                      alpha = beta
!                          !
!                          ! tmpforce(:) = tmpforce(:) + tmpRho(kpa,beta)&
!                          !     &*tmpRho(mu,beta)*ovrPrimeTmp(kpa,mu,:)
!                          ! tmpforce_r(:) = tmpforce_r(:) + tmpRho(kpa,beta)&
!                          !     &*tmpRho(mu,beta)*ovrPrimeTmp(mu,kpa,:)
                         
!                          ! tmpforce2 = tmpforce2 + tmpovr(kpa,beta)*tmpovr(beta&
!                          !     &,mu)*(tmpRho(kpa,mu)*tmpRho(beta,beta) +&
!                          !     & tmpRho(kpa,beta)*tmpRho(beta,mu))
!                          tmpforce(:) = tmpforce(:) + tmpovr(alpha,beta)*(tmpRho(kpa,beta)&
!                              &*tmpRho(mu,alpha) + tmpRho(kpa,alpha)*tmpRho(mu&
!                              &,beta))*ovrPrimeTmp(kpa,mu,:)
                         
!                          tmpforce_r(:) = tmpforce_r(:) + tmpovr(alpha,beta)*(tmpRho(kpa,beta)&
!                              &*tmpRho(mu,alpha) + tmpRho(kpa,alpha)*tmpRho(mu&
!                              &,beta))*ovrPrimeTmp(mu,kpa,:)
                         
!                          ! tmpforce2 = tmpforce2 + tmpovr(kpa,beta)*tmpovr(beta&
!                          !     &,mu)*(tmpRho(kpa,mu)*tmpRho(beta,beta) +&
!                          !     & tmpRho(kpa,beta)*tmpRho(beta,mu))
!                       end do
!                    end do
!                 end do
!              end do
!               tmpderiv(:,iAtK) = tmpderiv(:,iAtK) + tmpforce(:)*tmpgamma1
!               tmpderiv(:,iAtC) = tmpderiv(:,iAtC) + tmpforce_r(:)*tmpgamma1
!              ! !
! !              tmpderiv(:,iAtK) = tmpderiv(:,iAtK) + tmpforce2*(gammaPrimeTmp(:&
! !                  &,iAtK,iAtC) )
!              ! !
! !              tmpderiv(:,iAtC) = tmpderiv(:,iAtC) + tmpforce2*(gammaPrimeTmp(:&
! !                  &,iAtC,iAtK) )
!              !
!           end do
!        end do
!     end do

    
    call cpu_time(finish2)
    write(*,'(a)') "Neighbor list based eval:"
    print '("-> done. time = ", f10.4, " sec.")', finish2&
        &-start2  


    write(*,*) "LR-FORCES 2:"
    write(*,*) tmpderiv(:,:)*(-0.25_dp)

    
    call cpu_time(start)
    !    
    tmpderiv = 0.0_dp
    tmp1 = 0.0_dp
    do iMovedAtom = 1, nAtom
       write(*,*) "Moving Atom", iMovedAtom
             call cpu_time(start2)

             do iAtom3 = 1, nAtom
                if(iAtom3 .ne. iMovedAtom) then
                   
                   ! A=1, B=2, C=3, K=Moved
                   !
                   do iAtom1 = 1, nAtom
                      do iAtom2 = 1, nAtom
                         !
                         scrcond1 = tstovr(iMovedAtom,iAtom2)*tstovr(iAtom1,iAtom3)&
                             &*(tstdensmat(iMovedAtom,iAtom3)*tstdensmat(iAtom2&
                             &,iAtom1)+tstdensmat(iMovedAtom,iAtom1)&
                             &*tstdensmat(iAtom2,iAtom3))
                         scrcond1 = 100.0_dp * maxval(abs(gammaPrimeTmp(:,iMovedAtom&
                             &,iAtom3))) * scrcond1
                         !
                         scrcond2 = tstovr(iAtom1,iAtom2) *&
                             & maxval(abs(tstovrprime(iMovedAtom,iAtom3,:)))&
                             &*(tstdensmat(iMovedAtom,iAtom2)*tstdensmat(iAtom3&
                             &,iAtom1)+tstdensmat(iMovedAtom,iAtom1)&
                             &*tstdensmat(iAtom3,iAtom2))
                         scrcond2 = scrcond2 * 100.0_dp *(self&
                             &%lrGammaEval(iMovedAtom,iAtom2)+ self%lrGammaEval(iAtom3,iAtom2))

                        scrcondtst = max(scrcond2,scrcond1)
                        !
                        if(scrcond1 .gt. 1.0e-16_dp) then
                           !

                           do beta=iSquare(iAtom2),iSquare(iAtom2+1)-1
                              do nu=iSquare(iMovedAtom),iSquare(iMovedAtom+1)-1
                                 do alpha=iSquare(iAtom1),iSquare(iAtom1+1)-1
                                    do mu=iSquare(iAtom3),iSquare(iAtom3+1)-1
                                       
                                       tmp1 = tmpRho(alpha,beta)*tmpRho(mu&
                                           &,nu) +&
                                           & tmpRho(mu,beta)*tmpRho(alpha,nu)
                                       
                                       ! tmpderiv(:,iMovedAtom)=tmpderiv(:,iMovedAtom)+tmp1&
                                       !     &*tmpovr(mu,alpha)*tmpovr(nu,beta)&
                                       !     &*gammaPrimeTmp(:,iMovedAtom&
                                       !     &,iAtom3)
                                       
                                    end do
                                 end do
                              end do
                           end do
                            end if
                            !
!                            if(iAtom1 /= iAtom2) then
                        if(scrcond2 .gt. 1.0e-16_dp) then
                            !
                            do beta=iSquare(iAtom2),iSquare(iAtom2+1)-1
                               do nu=iSquare(iMovedAtom),iSquare(iMovedAtom+1)-1
                                  do alpha=iSquare(iAtom1),iSquare(iAtom1+1)-1
                                     do mu=iSquare(iAtom3),iSquare(iAtom3+1)-1

                                        tmp1 = tmpRho(nu,beta)*tmpRho(mu,alpha)+tmpRho(mu&
                                            &,beta)*tmpRho(nu,alpha)
                                        !
                                        tmpderiv(:,iMovedAtom) = tmpderiv(:,iMovedAtom) +&
                                            & tmp1*(self&
                                            &%lrGammaEval(iMovedAtom,iAtom2)&
                                            &+self%lrGammaEval(iAtom3,iAtom2))&
                                            &*tmpovr(alpha,beta)&
                                            &*ovrPrimeTmp(nu,mu,:)
                                     end do
                                  end do
                               end do
                            end do
 !                        end if
                            
                         end if
                            
                       end do
                   end do
                   
                end if
                !
             end do
             call cpu_time(finish2)
!             write(*,'(a,I5)') "Atom moved", iMovedAtom
!             print '("-> done. time = ", f10.4, " sec.")', finish2&
!                 &-start2  
             
          end do
          
          tmpderiv2 = -tmpderiv/4.0_dp
    ! update the gradients
    gradients=gradients+tmpderiv2


    write(*,*) "LR-FORCES 2:"
    write(*,*) tmpderiv2(:,:)

    
    call cpu_time(finish)
    print '("-> addGradients_tst done. time = ", f10.4, " sec.")', finish&
        &-start  
 
    deallocate(gammaPrimeTmp)
    deallocate(tmpovr,tmpRho,ovrPrimeTmp,tmpderiv)
  end subroutine addGradients_tst

  


    

  !!* Adds gradients due to long-range HF-contribution.
  !!* @param self, class instance
  !!* @param gradients 
  !!* @param deltaRho, square difference DM (triangle form) 
  subroutine addGradients(self, gradients, derivator, deltaRho, skHamCont, skOverCont,&
      & coords, species,orb,iSquare,ovrlapMat,iNeighbor,nNeighbor)
    class(RangeSepFunc), intent(inout) :: self
    real(dp), intent(inout) :: gradients(:,:)
    class(NonSccDiff), intent(in) :: derivator
    real(dp), intent(in) :: deltaRho(:,:)
    type(OSlakoCont), intent(in) :: skHamCont, skOverCont
    real(dp), intent(in) :: coords(:,:)
    integer, intent(in) :: species(:)
    type(TOrbitals), intent(in) :: orb
    integer, intent(in) :: iSquare(:)
    real(dp), intent(in) :: ovrlapMat(:,:)
    !
    integer, intent(in) :: iNeighbor(0:,:)
    integer, intent(in) :: nNeighbor(:)
    !
    !
!    real(dp) :: Stmp(size(Sprime,dim=1), size(Sprime,dim=2),2,3)
    !
    integer :: nAtom, iAtom1, iAtom2, sp1, sp2, ii, jj, kk
    integer :: iAtom3, iAtom4, mu, nu, alpha, beta, iMovedAtom
    ! NOTE: dummy is not needed
    real(dp)  :: dummy(orb%mOrb,orb%mOrb,3), sPrimeTmp(orb%mOrb,orb%mOrb,3)
    real(dp) :: vect(3), tmp(3), lrenergy, tmp33, lrenergy2, tmp55
    real(dp) :: dist, rr, tmp1, tmp2(3), tmp3(3), lrgammatmp(3), tmpomat(3)
    real(dp), allocatable :: gammaPrimeTmp(:,:,:), ovrPrimeTmp(:,:,:),&
        & tmpovr(:,:), tmpRho(:,:), tmpderiv(:,:), tmpderiv2(:,:)
    ! arrays for screening
    real(dp), allocatable :: tstdensmat(:,:), tstovr(:,:), tstovrprime(:,:,:)
    real(dp) :: tsttmp_1, tsttmp_2, scrcond1, scrcond2, scrcondtst
    ! for cpu_start/finish routines
    real :: start, finish, start2, finish2
    !
    real(dp), allocatable :: sigmatst(:,:,:)
    !
    write(*,*) "rangesep.addGradients"
    !write(*,*) "deltaXDiff=", deltaXDiff
    !@:ASSERT(self%tInit)
    ! 
    @:ASSERT(size(gradients,dim=1) == 3)
    !

    !========================================== 
    ! orb%nOrbAtom(ii) => number of orbitals of atom ii
    nAtom = size(orb%nOrbAtom)
    ! deriv(:,:) = 0.0_dp

    !========================
    ! Brute force prototype
    !========================
    call cpu_time(start)

    ! arrays for screening
    allocate(tstdensmat(nAtom,nAtom))
    allocate(tstovr(nAtom,nAtom))
    allocate(tstovrprime(nAtom,nAtom,3))
    !
    allocate(tmpovr(size(ovrlapMat, dim=1),size(ovrlapMat, dim=1)))
    allocate(tmpRho(size(deltaRho, dim=1),size(deltaRho, dim=1)))

    allocate(gammaPrimeTmp(3,nAtom,nAtom))
    allocate(ovrPrimeTmp(size(deltaRho,dim=1),size(deltaRho,dim=1),3))
    allocate(tmpderiv(3,size(gradients,dim=2)))
    allocate(tmpderiv2(3,size(gradients,dim=2)))

    ! precompute values overlap times densmat
    allocate(sigmatst(size(ovrlapMat, dim=1),size(ovrlapMat, dim=1),nAtom))
    !

!    write(*,*) "OVERLAP:"
    do ii=1,size(ovrlapMat, dim=1)
       do jj=1,size(ovrlapMat, dim=1)
          if(jj .gt. ii) then
             tmpovr(ii,jj)=ovrlapMat(jj,ii)
          else
             tmpovr(ii,jj)=ovrlapMat(ii,jj)
          end if
!          write(*,'(E17.8)', advance='no') tmpovr(ii,jj)!overlap(ii,jj)
       end do
!              write(*,*)
    end do
!       write(*,*) "DELTAP "
    do ii=1,size(deltaRho, dim=1)
       do jj=1,size(deltaRho, dim=1)
          if(jj .gt. ii) then
             tmpRho(ii,jj)=deltaRho(jj,ii)
          else
             tmpRho(ii,jj)=deltaRho(ii,jj)
          end if
 !                  write(*,'(E17.8)', advance='no') tmpRho(ii,jj)!deltaRho(ii,jj)
       end do
 !            write(*,*)
    end do


!!!!!!!!!!!!!


!!!!!!!!!!!!!
    sigmatst = 0.0_dp
    write(*,*) "precomputing the sigma"
    call cpu_time(start2)

    do iAtom1 = 1, nAtom
       do ii = 1, size(tmpovr,dim=1)
          do jj = 1, size(tmpovr,dim=1)
             !
             do mu = iSquare(iAtom1), iSquare(iAtom1+1)-1
sigmatst(jj,ii,iAtom1) = sigmatst(jj,ii,iAtom1) + tmpovr(jj,mu)*tmpRho(mu,ii)
             end do
          end do
       end do

    end do
    call cpu_time(finish2)
    print '("-> done. time = ", f10.4, " sec.")', finish2&
        &-start2  

    
    
    ! precompute the gamma derivatives
    write(*,*) "precomputing the lr-gamma derivatives"
    ovrPrimeTmp = 0.0_dp
    gammaPrimeTmp = 0.0_dp
    do iAtom1 = 1, nAtom
       do iAtom2 = 1, nAtom
          if(iAtom1 .ne. iAtom2) then
          call getGammaPrimeValue(self, tmp, iAtom1, iAtom2, coords, species)
          gammaPrimeTmp(:,iAtom1,iAtom2) = tmp(:)

          ! write(*,*) "##", iAtom1, iAtom2
           ! write(*,*) "GAMMAPRIMEX", gammaPrimeTmp(1,iAtom1,iAtom2)
           ! write(*,*) "GAMMAPRIMEY", gammaPrimeTmp(2,iAtom1,iAtom2)
           ! write(*,*) "GAMMAPRIMEZ", gammaPrimeTmp(3,iAtom1,iAtom2)
          ! write(*,*) "## gamma"
          ! write(*,*) self%lrGammaEval(iAtom1,iAtom2)
          !
           call derivator%getFirstDeriv(dummy, skHamCont, coords, &
               &species, iAtom1, iAtom2, orb)
           call derivator%getFirstDeriv(sPrimeTmp, skOverCont, coords, &
               &species, iAtom1, iAtom2, orb)
           ! write(*,*) "ORB1", orb%nOrbAtom(iAtom1)
           ! write(*,*) "ORB2", orb%nOrbAtom(iAtom2)
           ! write(*,*) "iSq, iSq+1", iSquare(iAtom1), iSquare(iAtom1+1)-1
           ! write(*,*) "iSq, iSq+1", iSquare(iAtom2), iSquare(iAtom2+1)-1

           ! do kk=1,3
           !    do ii=1, size(sPrimeTmp,dim=1)
           !       do jj=1, size(sPrimeTmp,dim=2)
           !          write(*,'(E17.8)', advance='no') sPrimeTmp(ii,jj,kk)
           !       end do
           !       write(*,*)
           !    end do
           !    write(*,*) "***"
           ! end do
           
           do ii =1, 3
           ovrPrimeTmp(iSquare(iAtom1):(iSquare(iAtom1+1)-1)&
               &,iSquare(iAtom2):(iSquare(iAtom2+1)-1),ii) = transpose(sPrimeTmp(1:orb&
               &%nOrbAtom(iAtom1),1:orb%nOrbAtom(iAtom2),ii))
        end do
        end if
       end do
    end do
    ovrPrimeTmp = ovrPrimeTmp * 0.5_dp


    !
    do iAtom1=1, nAtom
       do iAtom2=1,nAtom
          ! overlap
          tsttmp_1=maxval(abs(tmpovr((iSquare(iAtom1)):(iSquare(iAtom1+1)-1)&
              &,(iSquare(iAtom2)):(iSquare(iAtom2+1)-1))))
          tstovr(iAtom1,iAtom2)=abs(tsttmp_1)
          ! overlap prime
          tsttmp_1=maxval(abs(ovrPrimeTmp((iSquare(iAtom1)):(iSquare(iAtom1+1)-1)&
              &,(iSquare(iAtom2)):(iSquare(iAtom2+1)-1),1)))
          tstovrprime(iAtom1,iAtom2,1)=abs(tsttmp_1)
          tsttmp_1=maxval(abs(ovrPrimeTmp((iSquare(iAtom1)):(iSquare(iAtom1+1)-1)&
              &,(iSquare(iAtom2)):(iSquare(iAtom2+1)-1),2)))
          tstovrprime(iAtom1,iAtom2,2)=abs(tsttmp_1)
          tsttmp_1=maxval(abs(ovrPrimeTmp((iSquare(iAtom1)):(iSquare(iAtom1+1)-1)&
              &,(iSquare(iAtom2)):(iSquare(iAtom2+1)-1),3)))
          tstovrprime(iAtom1,iAtom2,3)=abs(tsttmp_1)
          ! densmat 
          tsttmp_1=maxval(abs(tmpRho((iSquare(iAtom1)):(iSquare(iAtom1+1)-1)&
              &,(iSquare(iAtom2)):(iSquare(iAtom2+1)-1))))
          tstdensmat(iAtom1,iAtom2)=abs(tsttmp_1)
       end do
    end do

    
    !===== TEST OUTPUT =====
    ! do iAtom1 = 1, nAtom
    !    do iAtom2 = 1, iAtom1-1
    !       !   
    !       write(*,*) "GAMMAPRIMEX", gammaPrimeTmp(1,iAtom1,iAtom2)
    !       write(*,*) "GAMMAPRIMEY", gammaPrimeTmp(2,iAtom1,iAtom2)
    !       write(*,*) "GAMMAPRIMEZ", gammaPrimeTmp(3,iAtom1,iAtom2)
    !    end do
    ! end do
    !===== TEST OUTPUT =====
    
    ! ovrPrimeTmp(1,4,:) = 0.0_dp
    ! ovrPrimeTmp(4,1,:) = 0.0_dp
    ! ovrPrimeTmp(5,8,:) = 0.0_dp
    ! ovrPrimeTmp(8,5,:) = 0.0_dp

    
    ! do kk=1,3
    !    write(*,*) "## ovrPrimeTmp", kk 
    !    do ii=1,size(ovrlapMat, dim=1)
    !       do jj=1,size(ovrlapMat, dim=1)
    !          write(*,'(E17.8)', advance='no') ovrPrimeTmp(ii,jj,kk)!overlap(ii,jj)
    !       end do
    !       write(*,*)
    !    end do
    ! end do


    
    
    ! orb%nOrbAtom(iAtom1) number orbitals atom1


!     tmpderiv = 0.0_dp
!     tmp1 = 0.0_dp
!     do iMovedAtom = 1, nAtom
!        write(*,*) "Moving Atom", iMovedAtom
!        do iAtom1 = 1, nAtom
!           do iAtom2 = 1, nAtom
! !             do ii = 1, nNeighbor(iMovedAtom)
! !                iAtom3 = iNeighbor(ii, iMovedAtom)
!              do iAtom3 = 1, nAtom
!                 if(iAtom3 .ne. iMovedAtom) then
!                 ! A=1, B=2, C=3, K=Moved
!                    !                   write(*,'(3I4)') iAtom1,iAtom2,iAtom3
!                    call cpu_time(start2)
!                    !
!                 do nu=iSquare(iMovedAtom),iSquare(iMovedAtom+1)-1
!                       do alpha=iSquare(iAtom1),iSquare(iAtom1+1)-1
!                          do beta=iSquare(iAtom2),iSquare(iAtom2+1)-1
!                             do mu=iSquare(iAtom3),iSquare(iAtom3+1)-1
!                             tmp1 = tmpRho(alpha,beta)*tmpRho(mu,nu) +&
!                                 & tmpRho(mu,beta)*tmpRho(alpha,nu)
!                             tmpderiv(:,iMovedAtom)=tmpderiv(:,iMovedAtom)+tmp1&
!                                 &*tmpovr(mu,alpha)*tmpovr(nu,beta)&
!                                 &*gammaPrimeTmp(:,iMovedAtom,iAtom3)

!                             tmp1 = tmpRho(nu,beta)*tmpRho(mu,alpha)+tmpRho(mu&
!                                 &,beta)*tmpRho(nu,alpha)
                            
!                             ! tmpderiv(:,iMovedAtom) = tmpderiv(:,iMovedAtom) +&
!                             !     & tmp1*(self%lrGammaEval(iMovedAtom,iAtom1)&
!                             !     &+self%lrGammaEval(iAtom3,iAtom1))&
!                             !     &*tmpovr(alpha,beta)*ovrPrimeTmp(nu,mu,:)
!                             tmpderiv(:,iMovedAtom) = tmpderiv(:,iMovedAtom) +&
!                                 & tmp1*(self%lrGammaEval(iMovedAtom,iAtom2)&
!                                 &+self%lrGammaEval(iAtom3,iAtom2))&
!                                 &*tmpovr(alpha,beta)*ovrPrimeTmp(nu,mu,:)
!                          end do
!                       end do
!                    end do
!                 end do
!                 call cpu_time(finish2)
!                 write(*,'(a,3I5)') "A,B,C", iAtom1, iAtom2, iAtom3
!                 print '("-> done. time = ", f10.4, " sec.")', finish2&
!                     &-start2  

!              end if
!                 !      
!              end do
!           end do
!        end do
!     end do
    
    tmpderiv = 0.0_dp
    tmp1 = 0.0_dp
    do iMovedAtom = 1, nAtom
       write(*,*) "Moving Atom", iMovedAtom
             call cpu_time(start2)

             do iAtom3 = 1, nAtom
                if(iAtom3 .ne. iMovedAtom) then
                   
                   ! A=1, B=2, C=3, K=Moved
                   !
                   do iAtom1 = 1, nAtom
                      do iAtom2 = 1, nAtom


                         scrcond1 = tstovr(iMovedAtom,iAtom2)*tstovr(iAtom1,iAtom3)&
                             &*(tstdensmat(iMovedAtom,iAtom3)*tstdensmat(iAtom2&
                             &,iAtom1)+tstdensmat(iMovedAtom,iAtom1)&
                             &*tstdensmat(iAtom2,iAtom3))
                         scrcond1 = 100.0_dp * maxval(abs(gammaPrimeTmp(:,iMovedAtom&
                             &,iAtom3))) * scrcond1
                         !                  write(*,'(a,I5,E17.8)') "scr1:", iMovedAtom, scrcond1

                         scrcond2 = tstovr(iAtom1,iAtom2) *&
                             & maxval(abs(tstovrprime(iMovedAtom,iAtom3,:)))&
                             &*(tstdensmat(iMovedAtom,iAtom2)*tstdensmat(iAtom3&
                             &,iAtom1)+tstdensmat(iMovedAtom,iAtom1)&
                             &*tstdensmat(iAtom3,iAtom2))
                         scrcond2 = scrcond2 * 100.0_dp *(self&
                             &%lrGammaEval(iMovedAtom,iAtom2)+ self%lrGammaEval(iAtom3,iAtom2))
                         !                   write(*,'(a,I5,E17.8)') "scr2:",
                         !                    iMovedAtom, scrcond2

                        scrcondtst = max(scrcond2,scrcond1)

                        if(scrcond1 .gt. 1.0e-6_dp) then
                           !

                           do beta=iSquare(iAtom2),iSquare(iAtom2+1)-1
                              do nu=iSquare(iMovedAtom),iSquare(iMovedAtom+1)-1
                                 do alpha=iSquare(iAtom1),iSquare(iAtom1+1)-1
                                    do mu=iSquare(iAtom3),iSquare(iAtom3+1)-1
                                       
                                       tmp1 = tmpRho(alpha,beta)*tmpRho(mu&
                                           &,nu) +&
                                           & tmpRho(mu,beta)*tmpRho(alpha,nu)
                                       
                                       tmpderiv(:,iMovedAtom)=tmpderiv(:,iMovedAtom)+tmp1&
                                           &*tmpovr(mu,alpha)*tmpovr(nu,beta)&
                                           &*gammaPrimeTmp(:,iMovedAtom&
                                           &,iAtom3)
                                       
                                    end do
                                 end do
                              end do
                           end do
                           ! do beta=iSquare(iAtom2),iSquare(iAtom2+1)-1
                           !    do nu=iSquare(iMovedAtom),iSquare(iMovedAtom+1)-1
                           !       do alpha=iSquare(iAtom1),iSquare(iAtom1+1)-1
                           !          do mu=iSquare(iAtom3),iSquare(iAtom3+1)-1
                           !             tmp1 = tmpRho(alpha,beta)*tmpRho(mu&
                           !                 &,nu) +&
                           !                 & tmpRho(mu,beta)*tmpRho(alpha,nu)
                           !             tmpderiv(:,iMovedAtom)=tmpderiv(:,iMovedAtom)+tmp1&
                           !                 &*tmpovr(mu,alpha)*tmpovr(nu,beta)&
                           !                 &*gammaPrimeTmp(:,iMovedAtom&
                           !                 &,iAtom3)
                                       
                           !          end do
                           !       end do
                           !    end do
                           ! end do
                        end if
                        !
                        if(scrcond2 .gt. 1.0e-6_dp) then
                            !
                            do beta=iSquare(iAtom2),iSquare(iAtom2+1)-1
                               do nu=iSquare(iMovedAtom),iSquare(iMovedAtom+1)-1
                                  do alpha=iSquare(iAtom1),iSquare(iAtom1+1)-1
                                     do mu=iSquare(iAtom3),iSquare(iAtom3+1)-1

                                        tmp1 = tmpRho(nu,beta)*tmpRho(mu,alpha)+tmpRho(mu&
                                            &,beta)*tmpRho(nu,alpha)
                                        !
                                        tmpderiv(:,iMovedAtom) = tmpderiv(:,iMovedAtom) +&
                                            & tmp1*(self&
                                            &%lrGammaEval(iMovedAtom,iAtom2)&
                                            &+self%lrGammaEval(iAtom3,iAtom2))&
                                            &*tmpovr(alpha,beta)&
                                            &*ovrPrimeTmp(nu,mu,:)
                                     end do
                                  end do
                               end do
                            end do
                         end if
                            
                           
                            
                            ! do beta=iSquare(iAtom2),iSquare(iAtom2+1)-1
                            !    do nu=iSquare(iMovedAtom),iSquare(iMovedAtom+1)-1
                            !       do alpha=iSquare(iAtom1),iSquare(iAtom1+1)-1
                            !          do mu=iSquare(iAtom3),iSquare(iAtom3+1)-1
                            !             tmp1 = tmpRho(alpha,beta)*tmpRho(mu&
                            !                 &,nu) +&
                            !                 & tmpRho(mu,beta)*tmpRho(alpha,nu)
                            !             tmpderiv(:,iMovedAtom)=tmpderiv(:,iMovedAtom)+tmp1&
                            !                 &*tmpovr(mu,alpha)*tmpovr(nu,beta)&
                            !                 &*gammaPrimeTmp(:,iMovedAtom&
                            !                 &,iAtom3)

                            !             tmp1 = tmpRho(nu,beta)*tmpRho(mu,alpha)+tmpRho(mu&
                            !                 &,beta)*tmpRho(nu,alpha)

                            !             ! tmpderiv(:,iMovedAtom) = tmpderiv(:,iMovedAtom) +&
                            !             !     & tmp1*(self%lrGammaEval(iMovedAtom,iAtom1)&
                            !             !     &+self%lrGammaEval(iAtom3,iAtom1))&
                            !             !     &*tmpovr(alpha,beta)*ovrPrimeTmp(nu,mu,:)
                            !             tmpderiv(:,iMovedAtom) = tmpderiv(:,iMovedAtom) +&
                            !                 & tmp1*(self&
                            !                 &%lrGammaEval(iMovedAtom,iAtom2)&
                            !                 &+self%lrGammaEval(iAtom3,iAtom2))&
                            !                 &*tmpovr(alpha,beta)&
                            !                 &*ovrPrimeTmp(nu,mu,:)
                            !          end do
                            !       end do
                            !    end do
                            ! end do
                   !
!                         end if
                      end do


                   end do
                   
                end if
                !
             end do
             call cpu_time(finish2)
             write(*,'(a,I5)') "Atom moved", iMovedAtom
             print '("-> done. time = ", f10.4, " sec.")', finish2&
                 &-start2  
             
          end do
          
          tmpderiv2 = -tmpderiv/4.0_dp

 !    write(*,*) "LR-FORCES 2:"
 !    write(*,*) tmpderiv2(:,:)

    !
    ! write(*,*) "NEIGHBORS:"
    ! write(*,*) "==================="
    ! do iAtom1 = 1, nAtom
    !    write(*,*) "iAtom1", iAtom1
    !    do ii = 1, nNeighbor(iAtom1)
    !       iAtom2 = iNeighbor(ii, iAtom1)
    !       write(*,*) ii, iAtom2
    !    end do
    ! end do

    ! write(*,*) "OVERLAP:"
    ! do ii=1,size(ovrlapMat, dim=1)
    !    do jj=1,size(ovrlapMat, dim=1)
    !       write(*,'(E17.8)', advance='no') tmpovr(ii,jj)!overlap(ii,jj)
    !    end do
    !    write(*,*)
    ! end do

    ! write(*,*) "DENS:"
    ! do ii = 1, nAtom
    !    do jj = 1, nAtom
    !       write(*,'(E17.8)',advance='no') tmpRho(ii,jj)
    !    end do
    !    write(*,*)
    ! end do


!==================
! FOR TEST ONLY
!==================
!     tmpderiv = 0.0_dp
!     do iMovedAtom = 1, nAtom
! !!!!!!!!!!!!!!!!!
!        ! +x
!        do ii=1,3
!           tmp55 = 0.0_dp
!           tmp = 0.0_dp
!           tmp(ii) = deltaXDiff!0.01_dp
!           lrenergy = 0.0_dp
! !       write(*,*) "============================================="
! !       write(*,*) "Energy_tst"
!           call evalEnergy_tst(self,coords,species,orb,iSquare,skOverCont,iMovedAtom,tmp&
!               &,tmpRho,lrenergy)
!        !    write(*,*) "ENERGY MOD+", lrenergy
!           tmp55 = lrenergy
!           ! -x
!           tmp = 0.0_dp
!           tmp(ii) = -deltaXDiff
!           lrenergy = 0.0_dp
!           call evalEnergy_tst(self,coords,species,orb,iSquare,skOverCont,iMovedAtom,tmp&
!               &,tmpRho,lrenergy)
!        !    write(*,*) "ENERGY MOD-", lrenergy
!           tmpderiv(ii,iMovedAtom) = 0.5_dp * (tmp55 - lrenergy) / deltaXDiff
!        end do
!     end do
!==================
    !

!    write(*,*) "LR-FORCES:"
!    write(*,*) tmpderiv(:,:)

    ! update the gradients
    gradients=gradients+tmpderiv2


    call cpu_time(finish)
    print '("-> addGradients done. time = ", f10.4, " sec.")', finish&
        &-start  
 
    deallocate(gammaPrimeTmp)
    deallocate(tmpovr,tmpRho,ovrPrimeTmp,tmpderiv)
  end subroutine addGradients







  !> Updates the Hamiltonian with the range separated contribution.
  !!
  !! \param densSqr  Square (unpacked) density matrix
  !! \param over  Sparse (packed) overlap matrix.
  !! \param gammas  Value of gamma for each atom pair. Shape: (nAtom, nAtom)
  !! \param iNeighbor  Neighbor indices.
  !! \param nNeighbor  Nr. of neighbors for each atom.
  !! \param iSquare  Position of each atom in the rows/columns of the square
  !!     matrices. Shape: (nAtom)
  !! \param iPair  Position of each (neighbor, atom) pair in the sparse
  !!   matrix. Shape: (0:maxNeighbor, nAtom)
  !! \param orb  Orbital information.
  !! \param HH  Square (unpacked) Hamiltonian to be updated.
  !!
!  subroutine updateHamiltonian(self,densSqr, over, gammas, iNeighbor, nNeighbor, &
!      & iSquare, iPair, orb, HH)
  subroutine updateHamiltonian_thr(self,densSqr, over, iNeighbor, nNeighbor, &
      & iSquare, iPair, orb, HH, overlap, species)
    class(RangeSepFunc), intent(inout) :: self
    !
    real(dp), dimension(:,:), target, intent(in) :: densSqr
    real(dp), dimension(:), intent(in) :: over
!    real(dp), dimension(:,:), intent(in) :: gammas
    integer, dimension(0:,:), intent(in) :: iNeighbor
    integer, dimension(:), intent(in) :: nNeighbor
    integer, dimension(:), intent(in) :: iSquare
    integer, dimension(0:,:), intent(in) :: iPair
    type(TOrbitals), intent(in) :: orb
    real(dp), dimension(:,:), intent(inout), target :: HH
    real(dp), intent(in) :: overlap(:,:)
    integer, intent(in) :: species(:)

    integer, parameter :: DESC_LEN = 3, ISTART = 1, IEND = 2, INORB = 3

    real(dp), dimension(orb%mOrb**2), target :: Sma, Sam, Snb, Sbn
    real(dp), dimension(orb%mOrb**2), target :: Pab, Pmb, Pan, Pmn
    real(dp), dimension(:,:), pointer :: pSma, pSam, pSnb, pSbn, pHH
    real(dp), dimension(:,:), pointer :: pPab, pPmb, pPan, pPmn
    real(dp) :: gamma1, gamma2, gammaTot

    integer :: nAtom, ii, jj
    integer :: iAtM, iAtN, iAtA, iAtB, iNeighN, iNeighA
    integer, dimension(DESC_LEN) :: descA, descB, descM, descN

    real(dp), dimension(:,:), allocatable :: gammas
    real(dp), dimension(:,:), allocatable, target :: tmpDRho, tmpDham
!!!
    real(dp), allocatable :: tmpovr(:,:), testovr(:,:)&
        &,tmpDDRho(:,:), tstham(:,:)
    integer :: iAtom_mu, iAtom_nu
    real(dp) :: tmp,tmp2,pbound,tstbound, gammabatch,prb!,hamtmp, hamtmp2,tmp3
    integer, allocatable :: ovrind(:,:), maxindex(:) !, tripindex(:,:,:)
    integer :: kk, iAtom1, iSp1, nOrb1
    real(dp), dimension(orb%mOrb,orb%mOrb), target :: Ovr_ma, Ovr_am, Ovr_nb, Ovr_bn
!!!
    real :: start, finish

    call cpu_time(start)
    nAtom = size(orb%nOrbAtom)
    allocate(gammas(size(self%lrGammaEval,dim=1),size(self%lrGammaEval,dim=1)))
    gammas(:,:) = self%lrGammaEval(:,:)
    allocate(tmpDRho(size(densSqr,dim=1),size(densSqr,dim=1)))
    !!!!!
    allocate(tmpovr(size(overlap, dim=1),size(overlap, dim=1)))
    allocate(tstham(size(overlap, dim=1),size(overlap, dim=1)))
    allocate(tmpDham(size(overlap, dim=1),size(overlap, dim=1)))
    allocate(tmpDDRho(size(densSqr, dim=1),size(densSqr, dim=1)))


    do ii=1,size(overlap, dim=1)
       do jj=1,size(overlap, dim=1)
          if(jj .gt. ii) then
             tmpovr(ii,jj)=overlap(jj,ii)
          else
             tmpovr(ii,jj)=overlap(ii,jj)
          end if
       end do
    end do

!    write(*,*) "DELTAP BEFORE"
    do ii=1,size(densSqr, dim=1)
       do jj=1,size(densSqr, dim=1)
           if(jj .gt. ii) then
              tmpDRho(ii,jj)=densSqr(jj,ii)
!              densSqr(ii,jj)=densSqr(jj,ii)
           else
              tmpDRho(ii,jj)=densSqr(ii,jj)
           end if
 !         write(*,'(E17.8)', advance='no') tmpDRho(ii,jj)!deltaRho(ii,jj)
       end do
  !     write(*,*)
    end do


    ! Hamiltonian Screening workaround
    if(.not. self%tScreeningInited) then
       write(*,*) "!!! Iinitialize Screening"
       self%tScreeningInited=.true.
       allocate(self%hprev(size(HH, dim=1),size(HH, dim=1)))
       allocate(self%dRhoprev(size(densSqr, dim=1),size(densSqr, dim=1)))
       self%hprev=0.0_dp
       self%dRhoprev=tmpDRho
    end if
    tmpDDRho=tmpDRho-self%dRhoprev
    tmpDham=0.0_dp
    self%dRhoprev=tmpDRho
    ! Hamiltonian Screening workaround
    nAtom=size(self%species)
    allocate(testovr(nAtom,nAtom))
    ! prepare the max overlaps
    do iAtom_mu=1, nAtom
       do iAtom_nu=1,nAtom
          tmp2=maxval(abs(tmpovr((iSquare(iAtom_mu)):(iSquare(iAtom_mu+1)-1)&
              &,(iSquare(iAtom_nu)):(iSquare(iAtom_nu+1)-1))))
          testovr(iAtom_mu,iAtom_nu)=tmp2
       end do
    end do


    ! sort 
    allocate(maxindex(nAtom))
    allocate(ovrind(nAtom,nAtom))
    do iAtom_mu=1, nAtom
       call index_heap_sort(ovrind(iAtom_mu,:),testovr(iAtom_mu,:))
    end do
    ! Hamiltonian Screening workaround
    pbound=maxval(abs(tmpDDRho))
    do iAtom_mu = 1, nAtom
       do kk=1, nAtom
          iAtom1=ovrind(iAtom_mu,nAtom+1-kk)
          iSp1 = species(iAtom1)
          nOrb1 = orb%nOrbSpecies(iSp1)
          prb=pbound*testovr(iAtom1,iAtom_mu)
          if(abs(prb) .le. self%pScreeningTreshold) then
             !
             maxindex(iAtom_mu) = kk
             exit
          end if
       end do
    end do

    call cpu_time(finish)
    print '("!!!preparation done. time = ", f10.4, " sec.")', finish&
        &-start  
    call cpu_time(start)

    tmpDham=0.0_dp
    loopN: do iAtN = 1, nAtom
       descN = getDescriptor(iAtN)
       loopB: do kk = 1, nAtom ! maxindex(iAtN)
          iAtB = ovrind(iAtN,nAtom+1-kk) !iNeighbor(iNeighN, iAtN)
          descB = getDescriptor(iAtB)
          prb=pbound*testovr(iAtB,iAtN)
          if(abs(prb) .ge. self%pScreeningTreshold) then
          loopA: do iAtA = 1, nAtom
             descA = getDescriptor(iAtA)
             gamma1 = gammas(iAtA, iAtN) + gammas(iAtA, iAtB)
             loopM: do jj = 1, nAtom
                iAtM = jj !iNeighbor(iNeighA, iAtA)
                descM = getDescriptor(iAtM)
                tstbound=prb*testovr(iAtA,iAtM)
                if(abs(tstbound) .ge. self%pScreeningTreshold) then
                   !
                   !
                   gamma2 = gammas(iAtM, iAtN) + gammas(iAtM, iAtB)
                   gammaTot = gamma1 + gamma2
                   !
                   call copyOverlapBlock_2(descB, descN, Sbn, pSbn)
                   call transposeBlock(pSbn, Snb, pSnb)
                   call copyDensityBlock(descA, descB, Pab, pPab)
                   call copyOverlapBlock_2(descM, descA, Sma, pSma)
                   call transposeBlock(pSma, Sam, pSam)
                   !
                   if (iAtM >= iAtN) then
                      call updateHamiltonianBlock(descM, descN, pSma, pSbn, pPab)
                   end if
                   if (iAtA >= iAtN .and. iAtM /= iAtA) then
                      call copyDensityBlock(descM, descB, Pmb, pPmb)
                      call updateHamiltonianBlock(descA, descN, pSam, pSbn, pPmb)
                   end if
                   if (iAtM >= iAtB .and. iAtN /= iAtB) then
                      call copyDensityBlock(descA, descN, Pan, pPan)
                      call updateHamiltonianBlock(descM, descB, pSma, pSnb, pPan)
                   end if
                   if (iAtA >= iAtB .and. iAtM /= iAtA .and. iAtN /= iAtB) then
                      call copyDensityBlock(descM, descN, Pmn, pPmn)
                      call updateHamiltonianBlock(descA, descB, pSam, pSnb, pPmn)
                   end if
                end if
             end do loopM
          end do loopA
       end if
       end do loopB
    end do loopN


!     tmpDham = 0.0_dp
!     loopM: do iAtM = 1, nAtom
!        descM = getDescriptor(iAtM)
!        loopA: do kk = 1, nAtom ! maxindex(iAtN)
!           iAtA = kk!ovrind(iAtN,nAtom+1-kk) !iNeighbor(iNeighN, iAtN)
!           descA = getDescriptor(iAtA)
! !          call transposeBlock(pSbn, Snb, pSnb)
!           loopN: do iAtN = 1, iAtM !nAtom
!              descN = getDescriptor(iAtN)
!              loopB: do jj = 1, nAtom
!                 iAtB = jj !iNeighbor(iNeighA, iAtA)
!                 descB = getDescriptor(iAtB)
!                 call copyDensityBlock(descA, descB, Pab, pPab)
!                 call copyOverlapBlock_2(descM, descA, Sma, pSma)
!                 call copyOverlapBlock_2(descB, descN, Sbn, pSbn)
! !                call transposeBlock(pSma, Sam, pSam)
! !                call transposeBlock(pSbn, Snb, pSnb)
!                 gamma1 = gammas(iAtA, iAtN) + gammas(iAtA, iAtB)
!                 gamma2 = gammas(iAtM, iAtN) + gammas(iAtM, iAtB)
!                 gammaTot = gamma1 + gamma2
!                 !
!                 call updateHamiltonianBlock(descM, descN, pSma, pSbn, pPab)
!              end do loopB
!           end do loopN
!        end do loopA
!     end do loopM

    !=================================
    ! add deltaH to H_LR 
    self%hprev = self%hprev + tmpDham
    ! add H_LR to current Hamiltonian
    HH = HH + self%hprev
    ! hamiltonian=hamiltonian+self%hprev*0.5_dp
    !=================================

    call cpu_time(finish)
    print '("!!!Balint Algorithm ham build done. time = ", f10.4, " sec.")', finish&
        &-start  
    
  contains

    function getDescriptor(iAt) result(desc)
      integer, intent(in) :: iAt
      integer, dimension(DESC_LEN) :: desc

      desc(:) = [ iSquare(iAt), iSquare(iAt + 1) - 1, &
          & iSquare(iAt + 1) - iSquare(iAt) ]

    end function getDescriptor


    subroutine copyOverlapBlock_2(desc1, desc2, localBlock, pLocalBlock) 
      integer, dimension(DESC_LEN), intent(in) :: desc1, desc2
      real(dp), dimension(:), target, intent(inout) :: localBlock
      real(dp), dimension(:,:), pointer, intent(out) :: pLocalBlock

      pLocalBlock(1:desc1(INORB), 1:desc2(INORB)) => &
          & localBlock(1:desc1(INORB)*desc2(INORB))
      pLocalBlock(:,:) = &
          & overlap(desc1(ISTART):desc1(IEND), desc2(ISTART):desc2(IEND))

    end subroutine copyOverlapBlock_2


    subroutine copyOverlapBlock(iAt, iNeigh, nOrbAt, nOrbNeigh, localBlock, &
        & pLocalBlock) 
      integer, intent(in) :: iAt, iNeigh, nOrbAt, nOrbNeigh
      real(dp), dimension(:), target, intent(inout) :: localBlock
      real(dp), dimension(:,:), pointer, intent(out) :: pLocalBlock
      
      integer :: ind

      ind = iPair(iNeigh, iAt) + 1
      localBlock(1:nOrbNeigh*nOrbAt) = over(ind:ind+nOrbNeigh*nOrbAt-1)
      pLocalBlock(1:nOrbNeigh, 1:nOrbAt) => localBlock(1:nOrbNeigh*nOrbAt)

    end subroutine copyOverlapBlock
      

    subroutine copyDensityBlock(desc1, desc2, localBlock, pLocalBlock)
      integer, dimension(DESC_LEN), intent(in) :: desc1, desc2
      real(dp), dimension(:), target, intent(inout) :: localBlock
      real(dp), dimension(:,:), pointer, intent(out) :: pLocalBlock
      
      pLocalBlock(1:desc1(INORB), 1:desc2(INORB)) => &
          & localBlock(1:desc1(INORB)*desc2(INORB))
!      pLocalBlock(:,:) = &
!          & densSqr(desc1(ISTART):desc1(IEND), desc2(ISTART):desc2(IEND))
      pLocalBlock(:,:) = &
          & tmpDDRho(desc1(ISTART):desc1(IEND), desc2(ISTART):desc2(IEND))

    end subroutine copyDensityBlock


    subroutine transposeBlock(orig, localBlock, pLocalBlock)
      real(dp), dimension(:,:), intent(in) :: orig
      real(dp), dimension(:), target, intent(out) :: localBlock
      real(dp), dimension(:,:), pointer, intent(out) :: pLocalBlock

      pLocalBlock(1:size(orig, dim=2), 1:size(orig, dim=1)) => &
          & localBlock(1:size(orig))
      pLocalBlock = transpose(orig)

    end subroutine transposeBlock


    subroutine updateHamiltonianBlock(descM, descN, pSma, pSbN, pPab)
      integer, dimension(DESC_LEN), intent(in) :: descM, descN
      real(dp), dimension(:,:), pointer, intent(in) :: pSma, pSbN, pPab

      real(dp), dimension(:,:), pointer :: pHmn

!      pHmn => HH(descM(ISTART):descM(IEND), descN(ISTART):descN(IEND))
      pHmn => tmpDHam(descM(ISTART):descM(IEND), descN(ISTART):descN(IEND))
      pHmn(:,:) = pHmn + gammaTot * matmul(matmul(pSma, pPab), pSbn) * (-0.125_dp)

    end subroutine updateHamiltonianBlock

  end subroutine updateHamiltonian_thr


 

  !!* add the long-range HF exchange part to the hamiltonian
  !!* @param self, class instance
  !!* @param overlap, square overlap matrix (triangle form)
  !!* @param deltaRho, square difference DM (triangle form)
  !!* @param iSquare, atomic block index array
  !!* @param hamiltonian, square hamiltonian (triangle form)
!  subroutine addSquareHamiltonian(self, overlap, deltaRho, iSquare,&
  !      & hamiltonian)
  !!* @param species list of all atomic species
  !!* @param iNeighbor neighbor list for atoms
  !!* @param nNeighbor number of neighbors of each atom
  !!* @param orb Information about the orbitals
  subroutine addSquareHamiltonian_opt(self, overlap, deltaRho, iSquare,&
      & hamiltonian,species,iNeighbor,nNeighbor,orb,iPair,over)
    class(RangeSepFunc), intent(inout) :: self
    real(dp), intent(in) :: overlap(:,:), deltaRho(:,:)
    integer, intent(in) :: iSquare(:)
    real(dp), intent(inout) :: hamiltonian(:,:)
    ! neighbor list variables
    integer, intent(in) :: species(:)
    integer, intent(in) :: iNeighbor(0:,:)
    integer, intent(in) :: nNeighbor(:)
    type(TOrbitals), intent(in) :: orb
    !
    integer, dimension(0:,:), intent(in) :: iPair
    real(dp), dimension(:), intent(in) :: over


    !
    integer :: iAtom1, iAtom2, nAtom!, iSp1, iSp2
    integer :: iAtom3, iAtom4, iNeigh1, iNeigh2, iSp1, iSp2, iSp3, iSp4
    integer :: nOrb1, nOrb2, nOrb3, nOrb4
    integer :: ind, nOrb12
    integer :: iAtom_mu, iAtom_nu!, iSp_mu, iSp_nu
    integer :: ii,jj,mu,nu,alpha,beta
    real(dp) :: tmp,tmp2,pbound,tstbound, gammabatch,prb!,hamtmp, hamtmp2,tmp3
!    integer :: tres2, tres3,tres4,tres5,tres6,tres7,tres8
!    real(dp) :: dist

    real(dp), allocatable :: tmpovr(:,:), tmpDRho(:,:), testovr(:,:)&
        &,tmpDDRho(:,:),tmpDham(:,:), tstham(:,:)
    !
    !    integer, allocatable :: countarray(:)
    real(dp), dimension(orb%mOrb**2), target :: TmpOvr1, TmpOvr2
    real(dp), dimension(:,:), pointer :: pTmpOvr1, pTmpOvr2
    real(dp) :: tmpvec1(orb%mOrb)
    !
    ! for cpu_start/finish routines
    real :: start, finish

    ! for energy
    real(dp) :: energy
integer :: iAtC, iAtD, iAtA, iAtB, iNeighC, iNeighD
    !
    !integer, parameter :: DESC_LEN = 3, ISTART = 1, IEND = 2, INORB = 3
    


    write(*,*) "Adding range separated contribution to Hamiltonian"

    ! start timing
    call cpu_time(start)
    !

    allocate(tmpovr(size(overlap, dim=1),size(overlap, dim=1)))
    allocate(tmpDham(size(overlap, dim=1),size(overlap, dim=1)))
    allocate(tmpDRho(size(deltaRho, dim=1),size(deltaRho, dim=1)))
    allocate(tmpDDRho(size(deltaRho, dim=1),size(deltaRho, dim=1)))
    !
    allocate(tstham(size(overlap, dim=1),size(overlap, dim=1)))

    !    write(*,*) "HAM BEFORE"
    do ii=1,size(hamiltonian, dim=1)
       do jj=1,size(hamiltonian, dim=1)
          !          write(*,'(E17.8)', advance='no') hamiltonian(ii,jj)
       end do
       !       write(*,*)
    end do
!       write(*,*) "OVR BEFORE"
    do ii=1,size(overlap, dim=1)
       do jj=1,size(overlap, dim=1)
          if(jj .gt. ii) then
             tmpovr(ii,jj)=overlap(jj,ii)
          else
             tmpovr(ii,jj)=overlap(ii,jj)
          end if
!                    write(*,'(E17.8)', advance='no') tmpovr(ii,jj)!overlap(ii,jj)
       end do
 !             write(*,*)
    end do
 !      write(*,*) "DELTAP BEFORE"
    do ii=1,size(deltaRho, dim=1)
       do jj=1,size(deltaRho, dim=1)
          if(jj .gt. ii) then
             tmpDRho(ii,jj)=deltaRho(jj,ii)
          else
             tmpDRho(ii,jj)=deltaRho(ii,jj)
          end if
  !                 write(*,'(E17.8)', advance='no') tmpDRho(ii,jj)!deltaRho(ii,jj)
       end do
  !           write(*,*)
    end do

    ! Hamiltonian Screening workaround
    if(.not. self%tScreeningInited) then
       write(*,*) "!!! Iinitialize Screening"
       self%tScreeningInited=.true.
       allocate(self%hprev(size(hamiltonian, dim=1),size(hamiltonian, dim=1)))
       allocate(self%dRhoprev(size(deltaRho, dim=1),size(deltaRho, dim=1)))
       self%hprev=0.0_dp
       self%dRhoprev=tmpDRho
    end if
    tmpDDRho=tmpDRho-self%dRhoprev
    tmpDham=0.0_dp
    self%dRhoprev=tmpDRho
    ! Hamiltonian Screening workaround

    
    nAtom=size(self%species)
    allocate(testovr(nAtom,nAtom))

    ! prepare the max overlaps
    do iAtom_mu=1, nAtom
       do iAtom_nu=1,nAtom
          tmp2=maxval(abs(tmpovr((iSquare(iAtom_mu)):(iSquare(iAtom_mu+1)-1)&
              &,(iSquare(iAtom_nu)):(iSquare(iAtom_nu+1)-1))))
          testovr(iAtom_mu,iAtom_nu)=tmp2
       end do
    end do

    ! Hamiltonian Screening workaround
    pbound=maxval(abs(tmpDDRho))
    ! Hamiltonian Screening workaround

    ! end timing
    call cpu_time(finish)
    print '("!!!preparation done. time = ", f10.4, " sec.")', finish&
        &-start  

    write(*,*) "size(tmpovr,dim=1):", size(tmpovr,dim=1)
    write(*,*) "nAtom:", nAtom

    ! start timing
    call cpu_time(start)

    ! get the position in packed overlap mat
    iAtom2 = 1
    iAtom1 = 1
    ind = iPair(iAtom2-1,iAtom1) + 1
    nOrb1 = ( iSquare(iAtom1+1) - iSquare(iAtom1) )
    nOrb2 = ( iSquare(iAtom2+1) - iSquare(iAtom2) )
    nOrb12 = nOrb1 * nOrb2
    TmpOvr1(1:nOrb12) = over(ind:ind+nOrb12-1)
    pTmpOvr1(1:nOrb2,1:nOrb1) => TmpOvr1(1:nOrb12)

    write(*,*) "OVERLAP BLOCK:", iAtom2, iAtom1, ind, nOrb2, nOrb1
    do ii = 1, nOrb2
       do jj = 1, nOrb1
          write(*,'(E17.8)',advance='no') pTmpOvr1(ii,jj)
       end do
       write(*,*)
    end do

    write(*,*) "======================================"
    write(*,*) " TEST:"

    !
    do iAtC = 1, nAtom !
       do iAtD = 1, nAtom
          !
          do iNeighC = 0, nNeighbor(iAtC)
             iAtA = iNeighbor(iNeighC, iAtC)
             !
             do iNeighD = 0, nNeighbor(iAtD)
                iAtB = iNeighbor(iNeighD, iAtD)
                !

                !
             end do
          end do
          !
       end do
    end do

    write(*,*) "======================================"
    
    
    
    ! Hamiltonian Screening workaround
    tmpDham=0.0_dp
    do iAtom_mu=1, nAtom
       do iAtom1=1, nAtom
          iSp1 = species(iAtom1)
          nOrb1 = orb%nOrbSpecies(iSp1)
          prb=pbound*testovr(iAtom1,iAtom_mu)
          if(abs(prb) .ge. self%pScreeningTreshold) then
             do iAtom_nu=1,iAtom_mu
                do iAtom2=1,nAtom
                   ! screening condition
                   tstbound=prb*testovr(iAtom2,iAtom_nu)
                   if(abs(tstbound) .ge. self%pScreeningTreshold) then
                      !
                      gammabatch=(self%lrGammaEval(iAtom_mu,iAtom2)&
                          &+self%lrGammaEval(iAtom1,iAtom2)&
                          &+self%lrGammaEval(iAtom_mu,iAtom_nu)&
                          &+self%lrGammaEval(iAtom1,iAtom_nu))
                      gammabatch=-0.125_dp*gammabatch
                      !
                      ! calculate the Q_AB
                      do nu=iSquare(iAtom_nu),iSquare(iAtom_nu+1)-1
                         !
                         jj = 0
                         do ii=iSquare(iAtom1),iSquare(iAtom1+1)-1
                            jj = jj + 1
                            tmpvec1(jj)=sum(tmpovr(iSquare(iAtom2):(iSquare(iAtom2+1)&
                                &-1),nu)*tmpDDRho(ii,iSquare(iAtom2):(iSquare(iAtom2+1)&
                                &-1)))
                         end do
                         tmp=0.0_dp
                         do mu=iSquare(iAtom_mu),iSquare(iAtom_mu+1)-1
                            tmp = sum(tmpovr(iSquare(iAtom1):iSquare(iAtom1+1)&
                                &-1,mu)*tmpvec1(1:nOrb1))
                            tmpDham(mu,nu)=tmpDham(mu,nu)+gammabatch*tmp
                         end do
                      end do
                      !
                   end if
                end do
             end do
          end if
       end do
    end do

    

    ! end timing
    call cpu_time(finish)
    print '("!!!!ham build done. time = ", f10.4, " sec.")', finish&
        &-start  
    
    self%hprev=self%hprev+tmpDham
    hamiltonian=hamiltonian+self%hprev
    
    ! evaluate the LR-energy from LR-Hamiltonian
    ! Note! it takes the current RhoIn, so it is not the final energy
    !
    tmp = 0.0_dp
    do iAtom_mu=1, nAtom
       do iAtom_nu=1, iAtom_mu
          do mu=iSquare(iAtom_mu), iSquare(iAtom_mu+1)-1
             do nu=iSquare(iAtom_nu), iSquare(iAtom_nu+1)-1
             !          write(*,'(E12.5)',advance='no') self%hprev(mu,nu)

                if(iAtom_nu .eq. iAtom_mu) then
                   tmp = tmp + self%hprev(mu,nu)*tmpDRho(mu,nu)
                else
                   tmp = tmp + 2.0_dp*self%hprev(mu,nu)*tmpDRho(mu,nu)
                end if
             end do
!             write(*,*)
          end do
       end do
    end do
    self%lrenergy = 0.5_dp*tmp
    write(*,*) "LR-ENERGY FROM HAMILTONIAN:", self%lrenergy
    
    write(*,'(a)') " -> Done."
    deallocate(tmpovr,tmpDRho)
  end subroutine addSquareHamiltonian_opt

 


  subroutine addSquareHamiltonian(self, overlap, deltaRho, iSquare,&
      & hamiltonian,species,iNeighbor,nNeighbor,orb)
    class(RangeSepFunc), intent(inout) :: self
    real(dp), intent(in) :: overlap(:,:), deltaRho(:,:)
    integer, intent(in) :: iSquare(:)
    real(dp), intent(inout) :: hamiltonian(:,:)
    ! neighbor list variables
    integer, intent(in) :: species(:)
    integer, intent(in) :: iNeighbor(0:,:)
    integer, intent(in) :: nNeighbor(:)
    type(TOrbitals), intent(in) :: orb

    !
    integer :: iAtom1, iAtom2, nAtom!, iSp1, iSp2
    integer :: iAtom3, iAtom4, iNeigh1, iNeigh2, iSp1, iSp2, iSp3, iSp4
    integer :: nOrb1, nOrb2, nOrb3, nOrb4
    integer :: iAtom_mu, iAtom_nu!, iSp_mu, iSp_nu
    integer :: ii,jj,mu,nu,alpha,beta, nnn
    real(dp) :: tmp,tmp2,pbound,tstbound, gammabatch,prb!,hamtmp, hamtmp2,tmp3
!    integer :: tres2, tres3,tres4,tres5,tres6,tres7,tres8
!    real(dp) :: dist

    real(dp), allocatable :: tmpovr(:,:), tmpDRho(:,:), testovr(:,:)&
        &,tmpDDRho(:,:),tmpDham(:,:), tstham(:,:)
    !
    !    integer, allocatable :: countarray(:)
    real(dp) :: tmpar1(orb%mOrb,orb%mOrb)
    real(dp) :: tmpar2(orb%mOrb,orb%mOrb)
    real(dp) :: tmpar3(orb%mOrb,orb%mOrb)
    real(dp) :: tmpar4(orb%mOrb,orb%mOrb)
    real(dp) :: tmpvec1(orb%mOrb)
    !
    ! for cpu_start/finish routines
    real :: start, finish

    ! for energy
    real(dp) :: energy
    

    write(*,*) "Adding range separated contribution to Hamiltonian"
    ! hamiltonian(:,:) = hamiltonian(:,:) + 0.0_dp

    ! start timing
    call cpu_time(start)
    !

    allocate(tmpovr(size(overlap, dim=1),size(overlap, dim=1)))
    allocate(tmpDham(size(overlap, dim=1),size(overlap, dim=1)))
    allocate(tmpDRho(size(deltaRho, dim=1),size(deltaRho, dim=1)))
    allocate(tmpDDRho(size(deltaRho, dim=1),size(deltaRho, dim=1)))
    !
    allocate(tstham(size(overlap, dim=1),size(overlap, dim=1)))

    !    write(*,*) "HAM BEFORE"
!    do ii=1,size(hamiltonian, dim=1)
!       do jj=1,size(hamiltonian, dim=1)
          !          write(*,'(E17.8)', advance='no') hamiltonian(ii,jj)
!       end do
       !       write(*,*)
!    end do
!       write(*,*) "OVR BEFORE"
    do ii=1,size(overlap, dim=1)
       do jj=1,size(overlap, dim=1)
          if(jj .gt. ii) then
             tmpovr(ii,jj)=overlap(jj,ii)
          else
             tmpovr(ii,jj)=overlap(ii,jj)
          end if
 !                   write(*,'(E17.8)', advance='no') tmpovr(ii,jj)!overlap(ii,jj)
       end do
 !             write(*,*)
    end do

    !TST only
    nnn=0
    do ii=1,size(overlap, dim=1)
       do jj=1,size(overlap, dim=1)
          if(tmpovr(ii,jj) .gt. 1.0e-16) then
             nnn=nnn+1
          end if
       end do
    end do

    write(*,*) 'NON-VANISH OVERLAP COUNT:', nnn
  !  write(*,'(a,I10)') 'NON-VANISH OVERLAP COUNT:', nnn
 !     write(*,*) "DELTAP BEFORE"
    do ii=1,size(deltaRho, dim=1)
       do jj=1,size(deltaRho, dim=1)
          if(jj .gt. ii) then
             tmpDRho(ii,jj)=deltaRho(jj,ii)
          else
             tmpDRho(ii,jj)=deltaRho(ii,jj)
          end if
!                   write(*,'(E17.8)', advance='no') tmpDRho(ii,jj)!deltaRho(ii,jj)
       end do
 !            write(*,*)
    end do

    ! Hamiltonian Screening workaround
    if(.not. self%tScreeningInited) then
       write(*,*) "!!! Iinitialize Screening"
       self%tScreeningInited=.true.
       allocate(self%hprev(size(hamiltonian, dim=1),size(hamiltonian, dim=1)))
       allocate(self%dRhoprev(size(deltaRho, dim=1),size(deltaRho, dim=1)))
       self%hprev=0.0_dp
       self%dRhoprev=tmpDRho
       !       self%dRhoprev=0.0_dp
    end if
    tmpDDRho=tmpDRho-self%dRhoprev
    tmpDham=0.0_dp
    self%dRhoprev=tmpDRho
    ! Hamiltonian Screening workaround

    
    nAtom=size(self%species)

    allocate(testovr(nAtom,nAtom))
!    allocate(countarray(16))

!    countarray=0

    ! prepare the max overlaps
    do iAtom_mu=1, nAtom
       do iAtom_nu=1,nAtom
          tmp2=maxval(abs(tmpovr((iSquare(iAtom_mu)):(iSquare(iAtom_mu+1)-1)&
              &,(iSquare(iAtom_nu)):(iSquare(iAtom_nu+1)-1))))
          testovr(iAtom_mu,iAtom_nu)=tmp2
       end do
    end do


    ! Hamiltonian Screening workaround
    !   pbound=maxval(abs(tmpDRho))
    pbound=maxval(abs(tmpDDRho))
    ! Hamiltonian Screening workaround

    ! end timing
    call cpu_time(finish)
    print '("!!!preparation done. time = ", f10.4, " sec.")', finish&
        &-start  

    write(*,*) "size(tmpovr,dim=1):", size(tmpovr,dim=1)
    write(*,*) "nAtom:", nAtom

    write(*,*) "!!! max(deltap):", maxval(abs(tmpDRho))
    write(*,*) "!!! max(ddeltap):", maxval(abs(tmpDDRho))

    ! start timing
!    call cpu_time(start)



!     ! Hamiltonian Screening workaround
!     !   pbound=maxval(abs(tmpDRho))
!     pbound=maxval(abs(tmpDRho))
!     ! Hamiltonian Screening workaround
    
!     tstham = 0.0_dp
!     !==================================
!     ! Hamiltonian from neighbour lists
!     do iAtom1 = 1, nAtom ! C
!        iSp1 = species(iAtom1)
!        nOrb1 = orb%nOrbSpecie(iSp1)
!        !
!        do iNeigh1 = 0, nNeighbor(iAtom1) ! C <= A
!           iAtom3 = iNeighbor(iNeigh1, iAtom1) ! A
!           iSp3 = species(iAtom3) 
!           nOrb3 = orb%nOrbSpecie(iSp3)
!           !
!           do iAtom2 = 1, nAtom ! D
!              iSp2 = species(iAtom2) 
!              nOrb2 = orb%nOrbSpecie(iSp2)

!              do iNeigh2 = 0, nNeighbor(iAtom2) !D <= B
!                 iAtom4 = iNeighbor(iNeigh2, iAtom2) ! B
!                 iSp4 = species(iAtom4) 
!                 nOrb4 = orb%nOrbSpecie(iSp4)
!                 ! 1 eq.
!                 ! get gamma
!                 tmp = (self%lrGammaEval(iAtom1,iAtom4)+self&
!                     &%lrGammaEval(iAtom1,iAtom2)+self&
!                     &%lrGammaEval(iAtom3,iAtom4)+self&
!                     &%lrGammaEval(iAtom3,iAtom2))
!                 !=================================                
!                 ! eq. 1
!                 do nu=iSquare(iAtom2),iSquare(iAtom2+1)-1
!                    ii=0
!                    do alpha=iSquare(iAtom3),iSquare(iAtom3+1)-1
!                       ii=ii+1
!                       tmpvec1(ii) = sum(tmpDRho(alpha,iSquare(iAtom4)&
!                           &:iSquare(iAtom4+1)-1)&
!                           &*tmpovr(iSquare(iAtom4):iSquare(iAtom4+1)&
!                           &-1,nu))
!                    end do
!                    do mu=iSquare(iAtom1),iSquare(iAtom1+1)-1
!                       if( mu .le. nu) then
!                          tstham(mu,nu) = tstham(mu,nu) +  sum(tmpovr(mu&
!                              &,iSquare(iAtom3):iSquare(iAtom3+1)-1)*tmpvec1(1:nOrb3))&
!                              &*tmp
!                       end if
!                    end do
!                 end do
!                 ! eq 2.
!                 if( iAtom4 /= iAtom2) then ! B .ne. D
!                    do beta=iSquare(iAtom4),iSquare(iAtom4+1)-1
!                       ii = 0
!                       do alpha=iSquare(iAtom3),iSquare(iAtom3+1)-1
!                          ii = ii + 1
!                          tmpvec1(ii) = sum(tmpDRho(alpha&
!                                       &,iSquare(iAtom2):iSquare(iAtom2+1)-1)&
!                                       &*tmpovr(iSquare(iAtom2):iSquare(iAtom2&
!                                       &+1)-1,beta))
                                  
!                       end do
!                       do mu=iSquare(iAtom1),iSquare(iAtom1+1)-1
!                          if( mu .le. beta) then
!                             tstham(mu,beta) = tstham(mu,beta) +&
!                                 & sum(tmpovr(mu,iSquare(iAtom3):iSquare(iAtom3&
!                                 &+1)-1)&
!                                 &*tmpvec1(1:nOrb3))*tmp
!                          end if
!                       end do
!                    end do
!                 end if


!                 ! 3 eq.
!                 if( iAtom3 /= iAtom1 ) then ! A .ne. C
!                    do nu=iSquare(iAtom2),iSquare(iAtom2+1)-1
!                       ii = 0
!                       do mu=iSquare(iAtom1),iSquare(iAtom1+1)-1
!                          ii = ii + 1
!                          !
!                          tmpvec1(ii) = sum(tmpDRho(mu,iSquare(iAtom4)&
!                              &:iSquare(iAtom4+1)-1)&
!                              &*tmpovr(iSquare(iAtom4)&
!                              &:iSquare(iAtom4+1)-1,nu))
                         
!                       end do
!                       do alpha=iSquare(iAtom3),iSquare(iAtom3+1)-1
!                          if( alpha .le. nu) then
!                             tstham(alpha,nu) = tstham(alpha,nu) +&
!                                 & sum(tmpovr(alpha&
!                                 &,iSquare(iAtom1):iSquare(iAtom1+1)-1)*tmpvec1(1:nOrb1))*tmp
!                          end if
!                       end do
!                    end do
!                    !
!                 end if
!                 ! 4 eq.
!                 if( (iAtom4 /= iAtom2) .and. (iAtom3 /= iAtom1 ) ) then
!                    !
!                    do beta=iSquare(iAtom4),iSquare(iAtom4+1)-1
!                       ii = 0
!                       do mu=iSquare(iAtom1),iSquare(iAtom1+1)-1
!                          ii = ii + 1
!                          !
!                          tmpvec1(ii) = sum(tmpDRho(mu,iSquare(iAtom2)&
!                              &:iSquare(iAtom2+1)-1)&
!                              &*tmpovr(iSquare(iAtom2)&
!                              &:iSquare(iAtom2+1)-1,beta))
                         
!                       end do
!                       do alpha=iSquare(iAtom3),iSquare(iAtom3+1)-1
!                          if( alpha .le. beta) then
                            
!                             tstham(alpha,beta) = tstham(alpha,beta) +&
!                                 & sum(tmpovr(alpha&
!                                 &,iSquare(iAtom1):iSquare(iAtom1+1)-1)*tmpvec1(1:nOrb1))*tmp
                            
                            
!                          end if
!                       end do
!                    end do
                   
!                 end if

                
                
! !                if(abs(tstbound) .ge. self%pScreeningTreshold ) then
!                 !                tmpar1(1:nOrb1,1:nOrb2) = 0.0_dp
! !====                
! !                 do nu=iSquare(iAtom2),iSquare(iAtom2+1)-1
! !                    do mu=iSquare(iAtom1),iSquare(iAtom1+1)-1
! !                       if( mu .le. nu) then
! !                          do alpha=iSquare(iAtom3),iSquare(iAtom3+1)-1
! !  !                           do beta=iSquare(iAtom4),iSquare(iAtom4+1)-1
! !                            ! tstham(mu,nu) = tstham(mu,nu) + tmpovr(mu,alpha)&
! !                             !     &*tmpDRho(alpha,beta)*tmpovr(beta,nu)*tmp
! !                             !===
! !                                tmp2 = 0.0_dp
! !                                tmp2 = sum(tmpDRho(alpha,iSquare(iAtom4)&
! !                                    &:iSquare(iAtom4+1)-1)&
! !                                    &*tmpovr(iSquare(iAtom4):iSquare(iAtom4+1)&
! !                                    &-1,nu))
! !                                tstham(mu,nu) = tstham(mu,nu) + tmpovr(mu,alpha)&
! !                                                                   &*tmp2*tmp
! !                                !===

! ! !                            end do
! !                          end do
! !                       end if
! !                    end do
! !                 end do
! !==== 
!                 !            end if
! !===                
! !                 ! 2 eq.
! !                 if( iAtom4 /= iAtom2) then ! B .ne. D
! ! !                   write(*,'(4I5)') iAtom1, iAtom3, iAtom2, iAtom4
! !                    !
! !                    do mu=iSquare(iAtom1),iSquare(iAtom1+1)-1
! !                       do beta=iSquare(iAtom4),iSquare(iAtom4+1)-1
! !                          if( mu .le. beta) then
! ! !                            do nu=iSquare(iAtom2),iSquare(iAtom2+1)-1
! !                                do alpha=iSquare(iAtom3),iSquare(iAtom3+1)-1
! !                                   ! tstham(mu,beta) = tstham(mu,beta) +&
! !                                   !     & tmpovr(mu,alpha)&
! !                                   !     &*tmpDRho(alpha,nu)*tmpovr(nu,beta)*tmp
! !                                   !
! !                                   tmp2 = 0.0_dp
! !                                   tmp2 = sum(tmpDRho(alpha&
! !                                       &,iSquare(iAtom2):iSquare(iAtom2+1)-1)&
! !                                       &*tmpovr(iSquare(iAtom2):iSquare(iAtom2&
! !                                       &+1)-1,beta))
! !                                   tstham(mu,beta) = tstham(mu,beta) +&
! !                                       & tmpovr(mu,alpha)&
! !                                       &*tmp2*tmp
                                  
! !                                end do
! ! !                            end do
! !                          end if
! !                       end do
! !                    end do
                   
! !                 end if

! !===
! !                 ! 3 eq.
! !                 if( iAtom3 /= iAtom1 ) then ! A .ne. C
! !  !                  write(*,'(4I5)') iAtom3, iAtom1, iAtom4, iAtom2
! !                    !
! !                       do nu=iSquare(iAtom2),iSquare(iAtom2+1)-1
! !                          do alpha=iSquare(iAtom3),iSquare(iAtom3+1)-1
! !                             if( alpha .le. nu) then
! !                                do mu=iSquare(iAtom1),iSquare(iAtom1+1)-1
! ! !                                  do beta=iSquare(iAtom4),iSquare(iAtom4+1)-1
! !                                       ! tstham(alpha,nu) = tstham(alpha,nu) +&
! !                                       !     & tmpovr(alpha,mu)&
! !                                       !     &*tmpDRho(mu,beta)*tmpovr(beta,nu)&
! !                                       !     &*tmp
! !                                       !
! !                                       tmp2=0.0_dp
! !                                       tmp2 = sum(tmpDRho(mu,iSquare(iAtom4)&
! !                                           &:iSquare(iAtom4+1)-1)&
! !                                           &*tmpovr(iSquare(iAtom4)&
! !                                           &:iSquare(iAtom4+1)-1,nu))
! !                                       tstham(alpha,nu) = tstham(alpha,nu) +&
! !                                           & tmpovr(alpha,mu)*tmp2*tmp

                                      
! !  !                                 end do
! !                                end do
! !                             end if
! !                       end do
! !                    end do
                   
! !                 end if
! !                 ! 4 eq.
! !                 if( (iAtom4 /= iAtom2) .and. (iAtom3 /= iAtom1 ) ) then
! !   !                 write(*,'(4I5)') iAtom3, iAtom1, iAtom2, iAtom4
! !                    !
! !                    do alpha=iSquare(iAtom3),iSquare(iAtom3+1)-1
! !                       do beta=iSquare(iAtom4),iSquare(iAtom4+1)-1
! !                          if( alpha .le. beta) then
! !                             do mu=iSquare(iAtom1),iSquare(iAtom1+1)-1
! ! !                               do nu=iSquare(iAtom2),iSquare(iAtom2+1)-1
! !                                    ! tstham(alpha,beta) = tstham(alpha,beta) +&
! !                                    !     & tmpovr(alpha,mu)&
! !                                    !     &*tmpDRho(mu,nu)*tmpovr(nu,beta)*tmp
! !                                    !
! !                                    tmp2=0.0_dp
! !                                    tmp2 = sum(tmpDRho(mu,iSquare(iAtom2)&
! !                                        &:iSquare(iAtom2+1)-1)&
! !                                        &*tmpovr(iSquare(iAtom2)&
! !                                        &:iSquare(iAtom2+1)-1,beta))
! !                                    tstham(alpha,beta) = tstham(alpha,beta) +&
! !                                        & tmpovr(alpha,mu)*tmp2*tmp

                                  
! !   !                             end do
! !                             end do
! !                          end if
! !                       end do
! !                    end do
                   
! !                end if
!              end do
!           end do
!        end do
!     end do
!     ! end timing
!     call cpu_time(finish)
!     print '("!!!!ham from neighbours build done. time = ", f10.4, " sec.")', finish&
!         &-start  

    
    ! write(*,*) "HAM, ALG1"
    ! do ii=1,size(tmpDRho,dim=1)
    !    do jj=1,size(tmpDRho,dim=1)
    !       write(*,'(E17.8)',advance='no') tstham(jj,ii)*(-0.125_dp)
    !    end do
    !    write(*,*)
    ! end do
    
    ! Hamiltonian Screening workaround
    !   pbound=maxval(abs(tmpDRho))
!    pbound=maxval(abs(tmpDDRho))
    ! Hamiltonian Screening workaround

    ! start timing
    call cpu_time(start)


    
    ! Hamiltonian Screening workaround
    tmpDham=0.0_dp
    do iAtom_mu=1, nAtom
       !      do iAtom_nu=1,iAtom_mu
       !
       do iAtom1=1, nAtom
          iSp1 = species(iAtom1)
          nOrb1 = orb%nOrbSpecies(iSp1)
          prb=pbound*testovr(iAtom1,iAtom_mu)
          if(abs(prb) .ge. self%pScreeningTreshold) then
             do iAtom_nu=1,iAtom_mu
                do iAtom2=1,nAtom
                   ! screening condition
                   tstbound=prb*testovr(iAtom2,iAtom_nu)
                   if(abs(tstbound) .ge. self%pScreeningTreshold) then
                      !
                      gammabatch=(self%lrGammaEval(iAtom_mu,iAtom2)&
                           &+self%lrGammaEval(iAtom1,iAtom2)&
                           &+self%lrGammaEval(iAtom_mu,iAtom_nu)&
                           &+self%lrGammaEval(iAtom1,iAtom_nu))
                      gammabatch=-0.125_dp*gammabatch
                      !
                      ! calculate the Q_AB
                      do nu=iSquare(iAtom_nu),iSquare(iAtom_nu+1)-1
                         !
                         jj = 0
                         do ii=iSquare(iAtom1),iSquare(iAtom1+1)-1
                            jj = jj + 1
                            tmpvec1(jj)=sum(tmpovr(iSquare(iAtom2):(iSquare(iAtom2+1)&
                                 &-1),nu)*tmpDDRho(ii,iSquare(iAtom2):(iSquare(iAtom2+1)&
                                 &-1)))
                         end do
                         tmp=0.0_dp
                         do mu=iSquare(iAtom_mu),iSquare(iAtom_mu+1)-1
                            tmp = sum(tmpovr(iSquare(iAtom1):iSquare(iAtom1+1)&
                                 &-1,mu)*tmpvec1(1:nOrb1))
                            tmpDham(mu,nu)=tmpDham(mu,nu)+gammabatch*tmp
                         end do
                      end do
                      !
                      !
                   end if

                end do
             end do
             end if
          end do
          !
!       end do
    end do

!     ! Hamiltonian Screening workaround
!     tmpDham=0.0_dp
!     do iAtom_mu=1, nAtom
!        do iAtom_nu=1,iAtom_mu
!           !
!           do iAtom1=1, nAtom
!              prb=pbound*testovr(iAtom1,iAtom_mu)
!              do iAtom2=1,nAtom
!                 ! screening condition
!                 tstbound=prb*testovr(iAtom2,iAtom_nu)
! !                if(abs(tstbound) .ge. 1.0e-3) then
!                 if(abs(tstbound) .ge. self%pScreeningTreshold) then
!                    !
!                    gammabatch=(self%lrGammaEval(iAtom_mu,iAtom2)&
!                        &+self%lrGammaEval(iAtom1,iAtom2)&
!                        &+self%lrGammaEval(iAtom_mu,iAtom_nu)&
!                        &+self%lrGammaEval(iAtom1,iAtom_nu))
!                    gammabatch=-0.125_dp*gammabatch
!                    do mu=iSquare(iAtom_mu),iSquare(iAtom_mu+1)-1
!                       do nu=iSquare(iAtom_nu),iSquare(iAtom_nu+1)-1
!                          !
!                          ! calculate the Q_AB
!                          tmp=0.0_dp
!                          do ii=iSquare(iAtom1),iSquare(iAtom1+1)-1
!                             tmp2=0.0_dp
!                             tmp2=sum(tmpovr(iSquare(iAtom2):(iSquare(iAtom2+1)&
!                                 &-1),nu)*tmpDDRho(ii,iSquare(iAtom2):(iSquare(iAtom2+1)&
!                                 &-1)))
!                             !                            do jj=iSquare(iAtom2),iSquare(iAtom2+1)-1
!                             !                               tmp2=tmp2+tmpovr(jj,nu)*tmpDDRho(ii,jj)
!                             !                            end do
!                             tmp=tmp+tmpovr(mu,ii)*tmp2
!                          end do
!                          tmpDham(mu,nu)=tmpDham(mu,nu)+gammabatch*tmp
!                          !                         tmpDham(mu,nu)=tmpDham(mu,nu)-0.125_dp&
!                          !                             &*gammabatch*tmp
!                          !                         hamiltonian(mu,nu)=hamiltonian(mu,nu)-0.125_dp&
!                          !                             &*gammabatch*tmp
!                       end do
!                    end do
!                    !
!                 end if

!              end do
!           end do
!           !
!        end do
!     end do

    

    !=================================
    ! add deltaH to H_LR 
    self%hprev=self%hprev+tmpDham
    ! add H_LR to current Hamiltonian
    hamiltonian=hamiltonian+self%hprev
!    hamiltonian=hamiltonian+self%hprev*0.5_dp
    !=================================

    ! end timing
    call cpu_time(finish)
!    write(*,'(2I6)') iAtom_mu,iAtom_nu

    print '("!!!!ham build done. time = ", f10.4, " sec.")', finish&
        &-start  


    
    !    hamiltonian = hamiltonian + 0.0_dp
    !====
!    hamiltonian=hamiltonian+self%hprev
 !   hamiltonian=hamiltonian+transpose(tstham(:,:))*(-0.125_dp)

    !
    ! write(*,*) "HAM, ALG2"
    ! do ii=1,size(tmpDRho,dim=1)
    !    do jj=1,size(tmpDRho,dim=1)
    !       if( ii < jj) then
    !          self%hprev(ii,jj) = 0.0_dp
    !       end if
    !          write(*,'(E17.8)',advance='no') self%hprev(ii,jj)
    !    end do
    !    write(*,*)
    ! end do

  !  write(*,*) "MAX ERROR:", maxval(self%hprev(:,:)-transpose(tstham(:,:))*(-0.125_dp))
    ! do ii=1,size(tmpDRho,dim=1)
    !    do jj=1,size(tmpDRho,dim=1)
    !       write(*,'(E17.8)',advance='no') self%hprev(ii,jj)-tstham(jj,ii)*(-0.125_dp)
    !    end do
    !    write(*,*)
    ! end do


    call cpu_time(start)

    
    ! evaluate the LR-energy from LR-Hamiltonian
    ! Note! it takes the current RhoIn, so it is not the final energy
    !
    tmp = 0.0_dp
    do iAtom_mu=1, nAtom
       do iAtom_nu=1, iAtom_mu
          do mu=iSquare(iAtom_mu), iSquare(iAtom_mu+1)-1
             do nu=iSquare(iAtom_nu), iSquare(iAtom_nu+1)-1
             !          write(*,'(E12.5)',advance='no') self%hprev(mu,nu)

                if(iAtom_nu .eq. iAtom_mu) then
                   tmp = tmp + self%hprev(mu,nu)*tmpDRho(mu,nu)
                else
                   tmp = tmp + 2.0_dp*self%hprev(mu,nu)*tmpDRho(mu,nu)
                end if
             end do
!             write(*,*)
          end do
       end do
    end do
    self%lrenergy = 0.5_dp*tmp
    write(*,*) "LR-ENERGY FROM HAMILTONIAN:", self%lrenergy


    call cpu_time(finish)
!    write(*,'(2I6)') iAtom_mu,iAtom_nu

    print '("!!!!energy eval done. time = ", f10.4, " sec.")', finish&
        &-start  


    
    !

    ! JUST FOR TESTING:
    ! evaluate the SCC-energy:
    ! energy = 0.0_dp
    ! do ii=1, size(hamiltonian,dim=1)
    !    do jj=1, ii!size(hamiltonian,dim=1)
    !       if(ii .eq. jj) then
    !          energy = energy + hamiltonian(ii,jj) * (tmpDRho(ii,jj) + 1.0_dp )
    !       else
    !          energy = energy + 2.0_dp * hamiltonian(ii,jj) * (tmpDRho(ii,jj) )
    !       end if
    !    end do
    ! end do

    !
    ! do ii=1,size(tmpDRho,dim=1)
    !    do jj=1,size(tmpDRho,dim=1)
    !       write(*,'(E17.8)',advance='no') tmpDRho(ii,jj)
    !    end do
    !    write(*,*)
    ! end do
    ! write(*,*)
    ! do ii=1,size(tmpDRho,dim=1)
    !    do jj=1,size(tmpDRho,dim=1)
    !       write(*,'(E17.8)',advance='no') hamiltonian(ii,jj)
    !    end do
    !    write(*,*)
    ! end do

    
    ! ! energy=energy * 0.5_dp
    ! write(*,*) "FULL-SCC-CONTRIBUTION TO ENERGY:", energy
    ! JUST FOR TESTING.

    write(*,'(a)') " -> Done."
    deallocate(tmpovr,tmpDRho)

  end subroutine addSquareHamiltonian



      subroutine set_cubic_spline(xx, fct, gama)
        real(dp), intent(in) :: xx(:)
        real(dp), intent(in) :: fct(:)
        real(dp), allocatable, intent(out) :: gama(:)
        !
        integer :: ii, kk, nn
        real(dp) :: p,qn,sig,un,yp1,ypn
        real(dp), allocatable :: u(:)
        !
        nn = size(xx)
        allocate(gama(nn))
        allocate(u(nn))
        !
        yp1=exp(-xx(1))
        ypn=exp(-xx(nn))
        ! natural spline
        gama(1)=0
        u(1)=0
        !
        !     gama(1)=-0.5_dp
        !     u(1)=(3.0_dp/(xx(2)-xx(1)))*((fct(2)-fct(1))/(xx(2)-xx(1))-yp1)
        !

        do ii=2,nn-1
           sig = (xx(ii)-xx(ii-1))/(xx(ii+1)-xx(ii-1))
           p=sig*gama(ii-1)+2.0_dp
           gama(ii)=(sig-1.0_dp)/p
           u(ii)=(6.0_dp*((fct(ii+1)-fct(ii))/(xx(ii+1)-xx(ii)) - (fct(ii)-fct(ii&
                &-1))/(xx(ii)-xx(ii-1)))/(xx(ii+1)-xx(ii-1))-sig*u(ii-1))/p
        end do
        ! natural spline
        qn=0.0_dp
        un=0.0_dp
        !
        !     qn=0.5_dp
        !     un=(3.0_dp/(xx(nn)-xx(nn-1)))*(ypn-(fct(nn)-fct(nn-1))/(xx(nn)-xx(nn-1)))
        !
        gama(nn)=(un-qn*u(nn-1))/(qn*gama(nn-1)+1.0_dp)
        do kk=nn-1,1,-1
           gama(kk)=gama(kk)*gama(kk+1)+u(kk)
        end do

        gama=0.0_dp

      end subroutine set_cubic_spline
  
  !> returns cubic spline value at specified point x
  !> Placed here, rather than in interpolation.F90, because only called here.
  !> (algorithm is based on NR)
  subroutine get_cubic_spline(xx, fct, dds, x, y)

    !> array with abscissae
    real(dp), intent(in) :: xx(:)

    !> array with ordinate
    real(dp), intent(in) :: fct(:)

    !> splines second derivatives, derived via set cubic splines
    real(dp), intent(in) :: dds(:)
   
    !> evaluation point
    real(dp), intent(in) :: x

    !> evaluated spline
    real(dp), intent(out) :: y
    !
    integer :: right, left, middle
    real(dp) :: step, A, B
    ! bisection
    left = 1
    right = size(xx)
    do
       if((right - left) <= 1) exit
       middle = (right + left)/2
       if( x >= xx(middle)) then
          left = middle
       else 
          right = middle
       end if
    end do
    step = xx(right) - xx(left)
    A = (xx(right) - x)/step
    B = (x - xx(left))/step
    ! calculate the spline value
    y = A*fct(left)+B*fct(right)+step*step/6.0_dp*(A*A - 1.0_dp)*A*dds(left)&
        &+step*step/6.0_dp*(B*B-1.0_dp)*B*dds(right)
  end subroutine get_cubic_spline



  ! !!* returns the long-range gamma
  ! !!* @param self, class instance
  ! !!* @param Sp1, species 1
  ! !!* @param Sp1, species 2
  ! !!* @param dist, distance
  ! function getGammaValue(self,Sp1,Sp2,dist)
  !   class(RangeSepFunc), intent(inout) :: self
  !   integer, intent(in) :: Sp1, Sp2
  !   real(dp), intent(in) :: dist
  !   real(dp) :: getGammaValue
  !   !
  !   real(dp) :: tmp, tmp2
  !   !
  !   if(abs(dist) .le. 1.0e-16_dp) then
  !      ! !!! THIS IS JUST AN APPROXIMATION !!!
  !      tmp2=0.02_dp
  !      call get_cubic_spline(self%grid,self%lrGamma(Sp1,Sp2,:,1),self%lrGamma(Sp1,Sp2,:,2),tmp2,tmp)
  !      getGammaValue=tmp*1.0_dp
  !   else
  !      ! !!! THIS IS JUST AN APPROXIMATION !!!
  !      ! note! gammas are generated to 10.0 Bohr,
  !      ! the interpolation for larger distances is undefined!!!
  !      if(abs(dist) .le. self%gamma_range) then
  !         call get_cubic_spline(self%grid,self%lrGamma(Sp1,Sp2,:,1),self%lrGamma(Sp1,Sp2,:,2),dist,tmp)
  !         getGammaValue=tmp*1.0_dp
  !      else
  !         ! for large distances gamma is equal to (1-exp(-kR))/R
  !         getGammaValue=(1.0_dp-exp(-self%omega*dist))/dist
  !      end if
       
  !   end if

  ! end function getGammaValue



  


!##############################
! dev. version:
!##############################
!   !
!   subroutine addSquareHamiltonian(self, overlap, deltaRho, iSquare,&
!       & hamiltonian)
!     class(RangeSepFunc), intent(inout) :: self
!     real(dp), intent(in) :: overlap(:,:), deltaRho(:,:)
!     integer, intent(in) :: iSquare(:)
!     real(dp), intent(inout) :: hamiltonian(:,:)
!     !
!     integer :: iAtom1, iAtom2, iSp1, iSp2, nAtom
!     integer :: iAtom_mu, iAtom_nu, iSp_mu, iSp_nu
!     integer :: ii,jj,mu,nu
!     real(dp) :: hamtmp, hamtmp2, tmp,tmp2,tmp3,pbound,tstbound, gammabatch,prb
!     integer :: tres2, tres3,tres4,tres5,tres6,tres7,tres8
!     real(dp) :: dist

!     real(dp), allocatable :: tmpovr(:,:), tmpDRho(:,:), testovr(:,:)&
!         &,tmpDDRho(:,:),tmpDham(:,:)
!     integer, allocatable :: countarray(:)

!     ! for cpu_start/finish routines
!     real :: start, finish

    
!     write(*,*) "Adding range separated contribution to Hamiltonian"
!     ! hamiltonian(:,:) = hamiltonian(:,:) + 0.0_dp

!     ! start timing
!     call cpu_time(start)
!     !


    
!     allocate(tmpovr(size(overlap, dim=1),size(overlap, dim=1)))
!     allocate(tmpDham(size(overlap, dim=1),size(overlap, dim=1)))
!     allocate(tmpDRho(size(deltaRho, dim=1),size(deltaRho, dim=1)))
!     allocate(tmpDDRho(size(deltaRho, dim=1),size(deltaRho, dim=1)))
   
! !    write(*,*) "HAM BEFORE"
!     do ii=1,size(hamiltonian, dim=1)
!        do jj=1,size(hamiltonian, dim=1)
! !          write(*,'(E17.8)', advance='no') hamiltonian(ii,jj)
!        end do
! !       write(*,*)
!     end do
!  !   write(*,*) "OVR BEFORE"
!     do ii=1,size(overlap, dim=1)
!        do jj=1,size(overlap, dim=1)
!           if(jj .gt. ii) then
!              tmpovr(ii,jj)=overlap(jj,ii)
!           else
!              tmpovr(ii,jj)=overlap(ii,jj)
!           end if
! !          write(*,'(E17.8)', advance='no') tmpovr(ii,jj)!overlap(ii,jj)
!        end do
! !       write(*,*)
!     end do
!  !   write(*,*) "DELTAP BEFORE"
!     do ii=1,size(deltaRho, dim=1)
!        do jj=1,size(deltaRho, dim=1)
!           if(jj .gt. ii) then
!              tmpDRho(ii,jj)=deltaRho(jj,ii)
!           else
!              tmpDRho(ii,jj)=deltaRho(ii,jj)
!           end if
!  !         write(*,'(E17.8)', advance='no') tmpDRho(ii,jj)!deltaRho(ii,jj)
!        end do
!  !      write(*,*)
!     end do

!     ! Hamiltonian Screening workaround
!     if(.not. self%tScreeningInited) then
!        write(*,*) "!!! Iinitialize Screening"
!        self%tScreeningInited=.true.
!        allocate(self%hprev(size(hamiltonian, dim=1),size(hamiltonian, dim=1)))
!        allocate(self%dRhoprev(size(deltaRho, dim=1),size(deltaRho, dim=1)))
!        self%hprev=0.0_dp
!        self%dRhoprev=tmpDRho
!     end if
!     tmpDDRho=tmpDRho-self%dRhoprev
!     tmpDham=0.0_dp
!     self%dRhoprev=tmpDRho
!     ! Hamiltonian Screening workaround

    
 
!     nAtom=size(self%species)

!     allocate(testovr(nAtom,nAtom))
!     allocate(countarray(16))

!     countarray=0
    
!     do iAtom_mu=1, nAtom
!        do iAtom_nu=1,nAtom
!           !
!           ! tmp3=0.0_dp
!           ! do mu=iSquare(iAtom_mu),iSquare(iAtom_mu+1)-1
!           !    do nu=iSquare(iAtom_nu),iSquare(iAtom_nu+1)-1
!           !       if(abs(tmp3) .le. abs(tmpovr(mu,nu))) then
!           !          tmp3=tmpovr(mu,nu)
!           !       end if
!           !    end do
!           ! end do

!            tmp2=maxval(abs(tmpovr((iSquare(iAtom_mu)):(iSquare(iAtom_mu+1)-1)&
!                &,(iSquare(iAtom_nu)):(iSquare(iAtom_nu+1)-1))))
!            testovr(iAtom_mu,iAtom_nu)=tmp2

!            ! if(tmp2 .lt. 1.0e-16_dp) then
!            !    countarray(1)=countarray(1)+1
!            ! end if
!            ! if(tmp2 .ge. 1.0e-16_dp .and. tmp2 .lt. 1.0e-14_dp) then
!            !    countarray(2)=countarray(2)+1
!            ! end if
!            ! if(tmp2 .ge. 1.0e-14_dp .and. tmp2 .lt. 1.0e-12_dp) then
!            !    countarray(3)=countarray(3)+1
!            ! end if
!            ! if(tmp2 .ge. 1.0e-12_dp .and. tmp2 .lt. 1.0e-10_dp) then
!            !    countarray(4)=countarray(4)+1
!            ! end if
!            ! if(tmp2 .ge. 1.0e-10_dp .and. tmp2 .lt. 1.0e-8_dp) then
!            !    countarray(5)=countarray(5)+1
!            ! end if
!            ! if(tmp2 .ge. 1.0e-8_dp .and. tmp2 .lt. 1.0e-6_dp) then
!            !    countarray(6)=countarray(6)+1
!            ! end if
!            ! if(tmp2 .ge. 1.0e-6_dp .and. tmp2 .lt. 1.0e-5_dp) then
!            !    countarray(7)=countarray(7)+1
!            ! end if
!            ! if(tmp2 .ge. 1.0e-5_dp .and. tmp2 .lt. 1.0e-4_dp) then
!            !    countarray(8)=countarray(8)+1
!            ! end if
!            ! if(tmp2 .ge. 1.0e-4_dp .and. tmp2 .lt. 1.0e-3_dp) then
!            !    countarray(9)=countarray(9)+1
!            ! end if
!            ! if(tmp2 .ge. 1.0e-3_dp .and. tmp2 .lt. 1.0e-2_dp) then
!            !    countarray(10)=countarray(10)+1
!            ! end if
!            ! if(tmp2 .ge. 1.0e-2_dp .and. tmp2 .lt. 1.0_dp) then
!            !    countarray(11)=countarray(11)+1
!            ! end if
!            ! if(tmp2 .ge. 1.0 ) then
!            !    countarray(12)=countarray(12)+1
!            ! end if

           
! !          write(*,*) iAtom_mu,iAtom_nu,tmp2,tmp3
! !             end do
! !          end do
!        end do
!     end do

    
!     ! Hamiltonian Screening workaround
!     !   pbound=maxval(abs(tmpDRho))
!     pbound=maxval(abs(tmpDDRho))
!      ! Hamiltonian Screening workaround

!     ! end timing
!     call cpu_time(finish)
!     print '("!!!preparation done. time = ", f10.4, " sec.")', finish&
!         &-start  
   
!     ! tmp=real(nAtom,dp)
!     ! tmp=tmp*tmp
!     ! write(*,*) "DIM:", tmp
!     ! write(*,*) "OVR COUNT:"
!     ! write(*,'(E17.8,E17.8, I10)') 1.0e-16_dp,real(countarray(1),dp)/tmp, countarray(1)
!     ! write(*,'(E17.8,E17.8, I10)') 1.0e-14_dp,real(countarray(2),dp)/tmp, countarray(2)
!     ! write(*,'(E17.8,E17.8, I10)') 1.0e-12_dp,real(countarray(3),dp)/tmp, countarray(3)
!     ! write(*,'(E17.8,E17.8, I10)') 1.0e-10_dp,real(countarray(4),dp)/tmp, countarray(4)
!     ! write(*,'(E17.8,E17.8, I10)') 1.0e-8_dp,real(countarray(5),dp)/tmp, countarray(5)
!     ! write(*,'(E17.8,E17.8, I10)') 1.0e-6_dp,real(countarray(6),dp)/tmp, countarray(6)
!     ! write(*,'(E17.8,E17.8, I10)') 1.0e-5_dp,real(countarray(7),dp)/tmp, countarray(7)
!     ! write(*,'(E17.8,E17.8, I10)') 1.0e-4_dp,real(countarray(8),dp)/tmp, countarray(8)
!     ! write(*,'(E17.8,E17.8, I10)') 1.0e-3_dp,real(countarray(9),dp)/tmp, countarray(9)
!     ! write(*,'(E17.8,E17.8, I10)') 1.0e-2_dp,real(countarray(10),dp)/tmp, countarray(10)
!     ! write(*,'(E17.8,E17.8, I10)') 1.0_dp,real(countarray(11),dp)/tmp, countarray(11)
!     ! write(*,'(E17.8,E17.8, I10)') 100.0_dp,real(countarray(12),dp)/tmp, countarray(12)

!     ! write(*,*) "OVR COUNT END!"

    
    
!   !   ! now evaluate the gammas for atom pairs
!   !   write(*,*) "Evaluating LR-GAMMA for atom pairs"
!   !   dist=0.0_dp
!   !   do iAtom1=1,nAtom
!   !      do iAtom2=1,iAtom1!nAtom
!   !         iSp1=self%species(iAtom1)
!   !         iSp2=self%species(iAtom2)
!   !         ! get distance
!   !         dist=0.0_dp
!   !         do ii=1,3
!   !            dist=dist+(self%coords(ii,iAtom1)-self%coords(ii&
!   !                &,iAtom2))**2
!   !         end do
!   !         dist=sqrt(dist)
!   !         ! evaluate gamma
!   !         self%lrGammaEval(iAtom1,iAtom2)=getGammaValue(self,iSp1,iSp2,dist)
!   !         self%lrGammaEval(iAtom2,iAtom1)=self%lrGammaEval(iAtom1,iAtom2)
!   ! !        write(*,'(2(a,E17.8))') "GAMMA:", self%lrGammaEval(iAtom1,iAtom2), ", dist=", dist
!   !      end do
!   !   end do

!     !
!     ! do iAtom_mu=1, nAtom
!     !    do iAtom_nu=1,nAtom
!     !       !
!     !       do mu=iSquare(iAtom_mu),iSquare(iAtom_mu+1)-1
!     !          do nu=iSquare(iAtom_nu),iSquare(iAtom_nu+1)-1
!     !             ! alpha, beta summation
!     !             hamtmp=0.0_dp
!     !             !
!     !             do iAtom1=1, nAtom
!     !                do iAtom2=1,nAtom
!     !                   !
!     !                   do ii=iSquare(iAtom1),iSquare(iAtom1+1)-1
!     !                      do jj=iSquare(iAtom2),iSquare(iAtom2+1)-1
!     !                         hamtmp=hamtmp+tmpovr(mu,ii)*tmpovr(jj,nu)*tmpDRho(ii&
!     !                             &,jj)*(self%lrGammaEval(iAtom_mu,iAtom2)&
!     !                             &+self%lrGammaEval(iAtom1,iAtom2)&
!     !                             &+self%lrGammaEval(iAtom_mu,iAtom_nu)&
!     !                             &+self%lrGammaEval(iAtom1,iAtom_nu))
!     !                      end do
!     !                   end do
!     !                end do
!     !             end do
!     !             hamiltonian(mu,nu)=hamiltonian(mu,nu)-1.0_dp/8.0_dp&
!     !                 &*hamtmp
!     !             !
!     !          end do
!     !       end do
!     !       !
!     !    end do
!     ! end do


!     ! do iAtom_mu=1, nAtom
!     !    do iAtom_nu=1,nAtom
!     !       !
!     !       do mu=iSquare(iAtom_mu),iSquare(iAtom_mu+1)-1
!     !          do nu=iSquare(iAtom_nu),iSquare(iAtom_nu+1)-1
!     !             ! alpha, beta summation
!     !             hamtmp=0.0_dp
!     !             !
!     !             do iAtom1=1, nAtom
!     !                do iAtom2=1,nAtom
!     !                   !
!     !                   tmp=0.0_dp
!     !                   do ii=iSquare(iAtom1),iSquare(iAtom1+1)-1
!     !                      do jj=iSquare(iAtom2),iSquare(iAtom2+1)-1
!     !                         tmp=tmp+tmpovr(mu,ii)*tmpovr(jj,nu)*tmpDRho(ii&
!     !                             &,jj)
!     !                      end do
!     !                   end do
!     !                   hamtmp=hamtmp+tmp*(self%lrGammaEval(iAtom_mu,iAtom2)&
!     !                             &+self%lrGammaEval(iAtom1,iAtom2)&
!     !                             &+self%lrGammaEval(iAtom_mu,iAtom_nu)&
!     !                             &+self%lrGammaEval(iAtom1,iAtom_nu))
!     !                end do
!     !             end do
!     !             hamiltonian(mu,nu)=hamiltonian(mu,nu)-1.0_dp/8.0_dp&
!     !                 &*hamtmp
!     !             !
!     !          end do
!     !       end do
!     !       !
!     !    end do
!     ! end do

!     write(*,*) "size(tmpovr,dim=1):", size(tmpovr,dim=1)
!     write(*,*) "nAtom:", nAtom

!     write(*,*) "!!! max(deltap):", maxval(abs(tmpDRho))
!     write(*,*) "!!! max(ddeltap):", maxval(abs(tmpDDRho))

    
!     ! start timing
!     call cpu_time(start)

!     tres2=0
!     tres3=0
!     tres4=0

    
! !     do iAtom_mu=1, nAtom


! !        do iAtom_nu=1,iAtom_mu!nAtom
! !           !

! !           do mu=iSquare(iAtom_mu),iSquare(iAtom_mu+1)-1
! !              do nu=iSquare(iAtom_nu),iSquare(iAtom_nu+1)-1
! !                 ! alpha, beta summation
! !                 hamtmp=0.0_dp
! !                 !
! !                 tmp3=0.0_dp

! !                 do iAtom1=1, nAtom
! !                    do iAtom2=1,nAtom
! !                       !
! !                       tmp=0.0_dp


! !                       ! do ii=iSquare(iAtom1),iSquare(iAtom1+1)-1
! !                       !    do jj=iSquare(iAtom2),iSquare(iAtom2+1)-1
! !                       !       tmp=tmp+tmpovr(mu,ii)*tmpovr(jj,nu)*tmpDRho(ii,jj)
! !                       !    end do
! !                       ! end do

! !                       tstbound=testovr(iAtom1,iAtom_mu)*testovr(iAtom2&
! !                           &,iAtom_nu)*pbound
! !                       if(abs(tstbound) .ge. 1.0e-6) then
                         
! !                          do ii=iSquare(iAtom1),iSquare(iAtom1+1)-1
! !                             tmp2=0.0_dp
! ! !                            if(abs(tmpovr(mu,ii)) .gt. 1.0e-16_dp ) then
! !                                !
! !                                do jj=iSquare(iAtom2),iSquare(iAtom2+1)-1
! !                                   tmp2=tmp2+tmpovr(jj,nu)*tmpDRho(ii,jj)
! !                                end do
! !                                tmp=tmp+tmpovr(mu,ii)*tmp2
! ! !                            end if
! !                          end do
                         
! !                       end if
! !                       hamtmp=hamtmp+tmp*(self%lrGammaEval(iAtom_mu,iAtom2)&
! !                           &+self%lrGammaEval(iAtom1,iAtom2)&
! !                           &+self%lrGammaEval(iAtom_mu,iAtom_nu)&
! !                           &+self%lrGammaEval(iAtom1,iAtom_nu))
! !                    end do
! !                 end do

! ! !                tmp3=-0.125_dp*hamtmp
! !                 ! if( abs(tmp3) .le. 1.0e-6_dp) then
! !                 !    tres2=tres2+1
! !                 ! end if
! !                 ! if( abs(tmp3) .le. 1.0e-8_dp) then
! !                 !    tres3=tres3+1
! !                 ! end if
! !                 ! if( abs(tmp3) .le. 1.0e-16_dp) then
! !                 !    tres4=tres4+1
! !                 ! end if
                
! !                 hamiltonian(mu,nu)=hamiltonian(mu,nu)-0.125_dp*hamtmp
! !                 !
                
! !              end do
! !           end do
! !           !

 
! !        end do
! !     end do


!     ! do iAtom_mu=1, nAtom
!     !    do iAtom_nu=1,iAtom_mu
!     !       !
!     !       do iAtom1=1, nAtom
!     !          prb=pbound*testovr(iAtom1,iAtom_mu)
!     !          do iAtom2=1,nAtom
!     !             ! screening condition
!     !             tstbound=prb*testovr(iAtom2,iAtom_nu)
!     !             if(abs(tstbound) .ge. 1.0e-8) then
!     !                !
!     !                gammabatch=(self%lrGammaEval(iAtom_mu,iAtom2)&
!     !                    &+self%lrGammaEval(iAtom1,iAtom2)&
!     !                    &+self%lrGammaEval(iAtom_mu,iAtom_nu)&
!     !                    &+self%lrGammaEval(iAtom1,iAtom_nu))
!     !                do mu=iSquare(iAtom_mu),iSquare(iAtom_mu+1)-1
!     !                   do nu=iSquare(iAtom_nu),iSquare(iAtom_nu+1)-1
!     !                      !
!     !                      ! calculate the Q_AB
!     !                      tmp=0.0_dp
!     !                      do ii=iSquare(iAtom1),iSquare(iAtom1+1)-1
!     !                         tmp2=0.0_dp
!     !                         do jj=iSquare(iAtom2),iSquare(iAtom2+1)-1
!     !                            tmp2=tmp2+tmpovr(jj,nu)*tmpDRho(ii,jj)
!     !                         end do
!     !                         tmp=tmp+tmpovr(mu,ii)*tmp2
!     !                      end do
!     !                      hamiltonian(mu,nu)=hamiltonian(mu,nu)-0.125_dp&
!     !                          &*gammabatch*tmp
!     !                   end do
!     !                end do
!     !                !
!     !             end if

!     !          end do
!     !       end do
!     !       !
!     !    end do
!     ! end do


!     ! Hamiltonian Screening workaround
!     tmpDham=0.0_dp
!    do iAtom_mu=1, nAtom
!        do iAtom_nu=1,iAtom_mu
!           !
!           do iAtom1=1, nAtom
!              prb=pbound*testovr(iAtom1,iAtom_mu)
!              do iAtom2=1,nAtom
!                 ! screening condition
!                 tstbound=prb*testovr(iAtom2,iAtom_nu)
!                 if(abs(tstbound) .ge. 1.0e-6) then
!                    !
!                    gammabatch=(self%lrGammaEval(iAtom_mu,iAtom2)&
!                        &+self%lrGammaEval(iAtom1,iAtom2)&
!                        &+self%lrGammaEval(iAtom_mu,iAtom_nu)&
!                        &+self%lrGammaEval(iAtom1,iAtom_nu))
!                    do mu=iSquare(iAtom_mu),iSquare(iAtom_mu+1)-1
!                       do nu=iSquare(iAtom_nu),iSquare(iAtom_nu+1)-1
!                          !
!                          ! calculate the Q_AB
!                          tmp=0.0_dp
!                          do ii=iSquare(iAtom1),iSquare(iAtom1+1)-1
!                             tmp2=0.0_dp
!                             do jj=iSquare(iAtom2),iSquare(iAtom2+1)-1
!                                tmp2=tmp2+tmpovr(jj,nu)*tmpDDRho(ii,jj)
!                             end do
!                             tmp=tmp+tmpovr(mu,ii)*tmp2
!                          end do
!                          tmpDham(mu,nu)=tmpDham(mu,nu)-0.125_dp&
!                              &*gammabatch*tmp
! !                         hamiltonian(mu,nu)=hamiltonian(mu,nu)-0.125_dp&
! !                             &*gammabatch*tmp
!                       end do
!                    end do
!                    !
!                 end if

!              end do
!           end do
!           !
!        end do
!     end do

!     self%hprev=self%hprev+tmpDham
    
!     hamiltonian=hamiltonian+self%hprev

    
    

! !    write(*,*) "ZERO COUNT:", tres2, tres3, tres4
    
!        ! end timing
!        call cpu_time(finish)
!        write(*,'(2I6)') iAtom_mu,iAtom_nu

!        print '("!!!!ham build done. time = ", f10.4, " sec.")', finish&
!            &-start  
!        !

    
!     write(*,'(a)') " -> Done."

!     ! write(*,*) "HAM AFTER"
!     ! do ii=1,size(hamiltonian, dim=1)
!     !    do jj=1,size(hamiltonian, dim=1)
!     !       write(*,'(E17.8)', advance='no') hamiltonian(ii,jj)
!     !    end do
!     !    write(*,*)
!     ! end do
    
!     deallocate(tmpovr,tmpDRho)
    
!   end subroutine addSquareHamiltonian


  
end module rangeseparated
