! Main Dynamics module

module timeprop_module

  use globalenv
  use commontypes ! for types
  use potentials ! for TPotentials type
  use scc
  use shift
  use accuracy
  use constants
  use sparse2dense
  use densitymatrix
  use blasroutines
  use lapackroutines
  use populations
  use blas
  use lapack
  use spin
  use repcont
  use energies, only: TEnergies, init
  use populations
  use repulsive
  use periodic
  use environment
  use timer
  use taggedoutput

  implicit none

  !> Data type to  initialize electronic dynamics variables from parser
  type TElecDynamicsInp
     real(dp) :: tdField, Dt, Omega, ReFieldPolVec(3), ImFieldPolVec(3)
     real(dp) :: Time0, Time1
     integer :: PolDir, Steps, writeFreq, PertType, EnvType, SpType, restartFreq
     logical :: tPopulations, tRestart, tWriteRestart
     real(dp) :: Phase
  end type TElecDynamicsInp

  !> Data type for electronic dynamics internal settings
  type TElecDynamics
     private
     real(dp) :: Field, Dt, Omega, Time0, Time1
     complex(cp) :: FieldDir(3)
     real(dp), allocatable :: TDFunction(:, :), phase
     integer :: Nsteps, writeFreq, PertType, EnvType, SpType
     integer :: nAtom, nOrbs, nSpin=1, currPolDir=1, restartFreq, fdAutotest
     integer, allocatable :: species(:), PolDirs(:)
     character(mc), allocatable :: speciesName(:)
     logical :: tPopulations, tSpinPol=.false.
     logical :: tRestart, tWriteRestart, tWriteAutotest
     logical :: tLaser = .false., tKick = .false., tEnvFromFile = .false.
     type(TScc), allocatable :: sccCalc
  end type TElecDynamics

  public :: runDynamics, TElecDynamics_init
  public :: TElecDynamicsInp, TElecDynamics
  public :: iKick, iLaser, iNoTDPert
  public :: iTDConstant, iTDGaussian, iTDSin2, iTDFromFile
  public :: iTDSinglet, iTDTriplet

  !> Enumerating available types of perturbation
  integer, parameter :: iKick = 1, iLaser = 2, iNoTDPert = 3

  !> Enumerating available types of envelope function
  integer, parameter :: iTDConstant = 1, iTDGaussian = 2, iTDSin2 = 3, iTDFromFile = 4

  !> Enumerating available types of spin polarized spectra
  integer, parameter :: iTDSinglet = 1, iTDTriplet = 2


contains

  !> Initialisation of input variables
  subroutine TElecDynamics_init(this, inp, species, speciesName, tWriteAutotest, fdAutotest, &
       &autotestTag)

    !> ElecDynamics instance
    type(TElecDynamics), intent(out) :: this

    !> ElecDynamicsInp instance
    type(TElecDynamicsInp), intent(in) :: inp

    !> species of all atoms in the system
    integer, allocatable, intent(in) :: species(:)

    !> label for each atomic chemical species
    character(mc), allocatable, intent(in) :: speciesName(:)

    !> produce tagged output?
    logical, intent(in) :: tWriteAutotest

    !> File unit for autotest data
    integer, intent(in) :: fdAutotest

    !> Tagged output files (machine readable)
    character(*), intent(in) :: autotestTag

    real(dp) :: norm

    this%Field = inp%tdField
    this%Dt = inp%Dt
    this%Nsteps = inp%Steps
    this%PertType = inp%PertType
    this%EnvType = inp%EnvType
    this%SpType = inp%SpType
    this%tPopulations = inp%tPopulations
    this%tRestart = inp%tRestart
    this%tWriteRestart = inp%tWriteRestart
    this%phase = inp%Phase
    this%writeFreq = inp%writeFreq
    this%restartFreq = inp%restartFreq
    allocate(this%species, source=species)
    allocate(this%speciesName, source=speciesName)
    allocate(this%sccCalc)

    if (inp%EnvType /= iTDConstant) then
       this%Time0 = inp%Time0
       this%Time1 = inp%Time1
    end if

    if (inp%PertType == iLaser) then
       this%tLaser = .true.
    else if (inp%PertType == iKick) then
       this%tKick = .true.
    end if

    if (this%tLaser) then
       this%Omega = inp%Omega
       this%FieldDir = inp%ReFieldPolVec + imag * inp%ImFieldPolVec
       norm = sqrt(sum(this%FieldDir(:)*conjg(this%FieldDir(:))))
       this%FieldDir = this%FieldDir/norm ! normalize polarization vector
       allocate(this%TDFunction(0:this%Nsteps, 3))
       if (this%EnvType == iTDFromFile) then
          this%tEnvFromFile = .true.
       end if
       call getTDFunction(this)
    end if

    if (this%tKick) then
       if (inp%PolDir == 4) then
          allocate(this%PolDirs(3))
          this%PolDirs(:) = (/ 1, 2, 3 /)
       else
          allocate(this%PolDirs(1))
          this%PolDirs(1) = inp%PolDir
       end if
    end if

    this%tWriteAutotest = tWriteAutotest
    if (tWriteAutotest) then
       this%fdAutotest = fdAutotest
       open(fdAutotest, file=autotestTag, position="append")
    end if

  end subroutine TElecDynamics_init


  !> Driver of time dependent propagation to calculate wither spectrum or laser
  subroutine runDynamics(this, Hsq, ham, H0, q0, over, filling, neighborList, nNeighbor,&
       &iSquare, iPair, img2CentCell, orb, coord, W, pRepCont, sccCalc, env)

    !> ElecDynamics instance
    type(TElecDynamics) :: this

    !> Eigenvectors
    real(dp), intent(inout) :: Hsq(:,:,:)

    !> Sparse storage for non-SCC hamitonian
    real(dp), intent(inout) :: H0(:)

    !> reference atomic occupations
    real(dp), intent(inout) :: q0(:,:,:)

    !> resulting hamitonian (sparse)
    real(dp), allocatable, intent(inout) :: ham(:,:)

    !> overlap (sparse)
    real(dp), allocatable, intent(inout) :: over(:)

    !> atomic coordinates
    real(dp), allocatable, intent(inout) :: coord(:,:)

    !> spin constants
    real(dp), intent(in) :: W(:,:,:)

    !> occupations
    real(dp), intent(in) :: filling(:,:,:)

    !> Number of neighbours for each of the atoms
    integer, intent(inout) :: nNeighbor(:)

    !> index array for location of atomic blocks in large sparse arrays
    integer, allocatable, intent(inout) :: iPair(:,:)

    !> image atoms to their equivalent in the central cell
    integer, allocatable, intent(inout) :: img2CentCell(:)

    !> Index array for start of atomic block in dense matrices
    integer, intent(in) :: iSquare(:)

    !> list of neighbours for each atom
    type(TNeighborList), intent(inout) :: neighborList

    !> repulsive information
    type(ORepCont), intent(in) :: pRepCont

    !> Atomic orbital information
    type(TOrbitals), intent(in) :: orb

    !> SCC module internal variables
    type(TScc), intent(in) :: sccCalc

    !> Environment settings
    type(TEnvironment), intent(inout) :: env

    integer :: iPol

    this%sccCalc = sccCalc

    this%nSpin = size(ham(:,:), dim=2)
    if (this%nSpin > 1) then
       this%tSpinPol = .true.
       call qm2ud(q0)
    end if

    this%nOrbs = size(Hsq, dim=1)
    this%nAtom = size(coord, dim=2)

    if (this%tKick) then
       do iPol = 1, size(this%PolDirs)
          this%currPolDir = this%PolDirs(iPol)
          call doDynamics(this,Hsq,ham,H0,q0,over,filling,neighborList,nNeighbor, &
               &iSquare,iPair,img2CentCell,orb,coord,W,pRepCont,env)
       end do
    else
       call doDynamics(this,Hsq,ham,H0,q0,over,filling,neighborList,nNeighbor, &
            &iSquare,iPair,img2CentCell,orb,coord,W,pRepCont,env)
    end if

  end subroutine runDynamics


  !> Runs the electronic dynamics of the system
  subroutine doDynamics(this,Hsq,ham,H0,q0,over,filling,neighborList,nNeighbor, &
       &iSquare,iPair,img2CentCell,orb,coord,W,pRepCont,env)

    !> ElecDynamics instance
    type(TElecDynamics) :: this

    !> Eigenvectors
    real(dp), intent(inout) :: Hsq(:,:,:)

    !> Sparse storage for non-SCC hamitonian
    real(dp), intent(inout) :: H0(:)

    !> reference atomic occupations
    real(dp), intent(inout) :: q0(:,:,:)

    !> resulting hamitonian (sparse)
    real(dp), allocatable, intent(inout) :: ham(:,:)

    !> overlap (sparse)
    real(dp), allocatable, intent(inout) :: over(:)

    !> atomic coordinates
    real(dp), allocatable, intent(inout) :: coord(:,:)

    !> spin constants
    real(dp), intent(in) :: W(:,:,:)

    !> occupations
    real(dp), intent(in) :: filling(:,:,:)

    !> Number of neighbours for each of the atoms
    integer, intent(inout) :: nNeighbor(:)

    !> index array for location of atomic blocks in large sparse arrays
    integer, allocatable, intent(inout) :: iPair(:,:)

    !> image atoms to their equivalent in the central cell
    integer, allocatable, intent(inout) :: img2CentCell(:)

    !> Index array for start of atomic block in dense matrices
    integer, intent(in) :: iSquare(:)

    !> list of neighbours for each atom
    type(TNeighborList), intent(inout) :: neighborList

    !> repulsive information
    type(ORepCont), intent(in) :: pRepCont

    !> Atomic orbital information
    type(TOrbitals), intent(in) :: orb

    !> Environment settings
    type(TEnvironment), intent(inout) :: env


    complex(cp) :: Ssqr(this%nOrbs,this%nOrbs), Sinv(this%nOrbs,this%nOrbs)
    complex(cp) :: Rho(this%nOrbs,this%nOrbs,this%nSpin), Rhoold(this%nOrbs,this%nOrbs,this%nSpin)
    complex(cp) :: H1(this%nOrbs,this%nOrbs,this%nSpin), T1(this%nOrbs,this%nOrbs)
    complex(cp) :: Rhonew(this%nOrbs,this%nOrbs,this%nSpin)
    complex(cp), allocatable :: Eiginv(:,:,:), EiginvAdj(:,:,:)
    real(dp) :: qq(orb%mOrb, this%nAtom, this%nSpin), deltaQ(this%nAtom,this%nSpin)
    real(dp) :: dipole(3,this%nSpin), chargePerShell(orb%mShell,this%nAtom,this%nSpin)
    real(dp), allocatable :: rhoPrim(:,:), ham0(:)
    real(dp) :: time, dTime, startTime = 0.0_dp, timeElec = 0.0_dp
    integer :: dipoleDat, qDat, energyDat, populDat(2)
    integer :: iStep = 0, iAtom, iSpin
    type(TPotentials) :: potential
    type(TEnergies) :: energy
    character(4) :: dumpIdx
    type(TTimer) :: loopTime

    call env%globalTimer%startTimer(globalTimers%elecDynInit)
    if (this%tRestart) then
       call readRestart(Rho, Ssqr, coord, startTime)
    end if

    call createMatrices(this, Rho, H1, Ssqr, Sinv, H0, ham0, &
         & over, ham, Hsq, filling, orb, rhoPrim, potential, &
         & neighborList%iNeighbor, nNeighbor, iSquare, iPair, img2CentCell, &
         & Eiginv, EiginvAdj, energy)

    call initTDOutput(this, dipoleDat, qDat, energyDat, populDat)

    call getChargeDipole(this, deltaQ, qq, dipole, q0, Rho, Ssqr, coord, iSquare) !q_0
    call updateH(this, H1, ham, over, ham0, qq, q0, coord, orb, potential, &  !H1_0
         &neighborList%iNeighbor, nNeighbor, iSquare, iPair, img2CentCell, &
         &iStep, chargePerShell, W, env)

    ! Apply kick to Rho if necessary
    if (this%tKick) then
       call kickDM(this, Rho, Ssqr, Sinv, iSquare, coord)
    end if

    ! Initialize electron dinamics
    Rhoold(:,:,:) = Rho
    call initializePropagator(this, this%Dt, Rho, Rhoold, H1, Sinv)

    call getTDEnergy(this, energy, rhoPrim, Rhoold, neighborList%iNeighbor, &
         & nNeighbor, orb, iSquare, iPair, img2CentCell, ham0, qq, q0, &
         & potential, chargePerShell, coord, pRepCont)

    call env%globalTimer%stopTimer(globalTimers%elecDynInit)
    ! End of initialization

    !! Main loop
    call env%globalTimer%startTimer(globalTimers%elecDynLoop)
    call loopTime%start()

    write(stdOut,"(A)")'Starting dynamics'

    do iStep = 0, this%Nsteps
       time = iStep * this%Dt + startTime

       if (.not. this%tRestart .or. iStep > 0) then
          call writeTDOutputs(this, dipoleDat, qDat, energyDat, time, energy, dipole, deltaQ, iStep)
       end if

       call getChargeDipole(this, deltaQ, qq, dipole, q0, Rho, Ssqr, coord, iSquare)
       call updateH(this, H1, ham, over, ham0, qq, q0, coord, orb, &
            & potential, neighborList%iNeighbor, nNeighbor, iSquare, iPair, &
            & img2CentCell, iStep, chargePerShell, W, env)

       if ((this%tWriteRestart) .and. (iStep > 0) .and. &
            & (mod(iStep, this%restartFreq) == 0)) then
          call writeRestart(Rho, Ssqr, coord, time)
       end if

       call getTDEnergy(this, energy, rhoPrim, Rho, neighborList%iNeighbor, &
            & nNeighbor, orb, iSquare, iPair, img2CentCell, ham0, qq, q0, &
            & potential, chargePerShell, coord, pRepCont)

       do iSpin = 1, this%nSpin
          call scal(H1(:,:,iSpin), imag)
          call propagateRho(this, Rhoold(:,:,iSpin), Rho(:,:,iSpin), &
               & H1(:,:,iSpin), Sinv, T1, 2.0_dp * this%Dt)
          call swap(Rhoold(:,:,iSpin), Rho(:,:,iSpin))

          if ((this%tPopulations) .and. (mod(iStep, this%writeFreq) == 0)) then
             call getTDPopulations(this, Rho, Eiginv, EiginvAdj, T1, populDat, time, iSpin)
          end if
       end do

       if (mod(iStep,this%Nsteps/10) == 0) then
          call loopTime%stop()
          timeElec  = loopTime%getWallClockTime()
          write(stdOut,"(A,2x,I6,2(2x,A,F10.6))") 'Step ', iStep, 'elapsed loop time: ', &
               & timeElec, 'average time per loop ', timeElec/(iStep+1)
       end if
    end do

    write(stdOut, "(A)") 'Dynamics finished OK!'
    call env%globalTimer%stopTimer(globalTimers%elecDynLoop)
    if (this%tWriteAutotest) then
       call writeTDAutotest(this, dipole, energy, deltaQ)
    end if
    call closeTDOutputs(this, dipoleDat, qDat, energyDat, populDat)
  end subroutine doDynamics



  !> Updates the hamiltonian with SCC and external TD field (if any) contributions
  subroutine updateH(this, H1, ham, over, ham0, qq, q0, coord, orb, potential, &
       &iNeighbor, nNeighbor, iSquare, iPair, img2CentCell, iStep, chargePerShell, W, env)

    !> ElecDynamics instance
    type(TElecDynamics) :: this

    !> Square hamiltonian
    complex(cp), intent(inout) :: H1(:,:,:)

    !> resulting hamitonian (sparse)
    real(dp), allocatable, intent(inout) :: ham(:,:)

    !> overlap (sparse)
    real(dp), allocatable, intent(inout) :: over(:)

    !> Sparse storage for non-SCC hamitonian
    real(dp), intent(in) :: ham0(:)

    !> atomic ocupations
    real(dp), intent(inout) :: qq(:,:,:)

    !> reference atomic occupations
    real(dp), intent(inout) :: q0(:,:,:)

    !> atomic coordinates
    real(dp), allocatable, intent(inout) :: coord(:,:)

    !> Atomic orbital information
    type(TOrbitals), intent(in) :: orb

    !> potential acting on the system
    type(TPotentials), intent(inout) :: potential

    !> Atomic neighbour data
    integer, intent(in) :: iNeighbor(0:,:)

    !> Number of neighbours for each of the atoms
    integer, intent(in) :: nNeighbor(:)

    !> Index array for start of atomic block in dense matrices
    integer, intent(in) :: iSquare(:)

    !> index array for location of atomic blocks in large sparse arrays
    integer, intent(in) :: iPair(0:,:)

    !> image atoms to their equivalent in the central cell
    integer, intent(in) :: img2CentCell(:)

    !> current step of the propagation
    integer, intent(in) :: iStep

    !> electrons in each atomi shell
    real(dp), intent(inout) :: chargePerShell(:,:,:)

    !> spin constants
    real(dp), intent(in) :: W(:,:,:)

    !> Environment settings
    type(TEnvironment), intent(inout) :: env

    real(dp) :: T2(this%nOrbs,this%nOrbs)
    real(dp) :: shellPot(orb%mShell, this%nAtom, this%nSpin), atomPot(this%nAtom, this%nSpin)
    integer :: iAtom, iSpin

    ham(:,1) = ham0(:)

    if (this%nSpin == 2) then
       ham(:,2) = 0.0_dp
       call ud2qm(qq)
       call ud2qm(q0)
       call getChargePerShell(qq, orb, this%species, chargePerShell)
    end if

    call clearPotential(potential)

    call this%sccCalc%updateCharges(env, qq, q0, orb, this%species, iNeighbor, img2CentCell)
    call this%sccCalc%getShiftPerAtom(atomPot(:,1))
    call this%sccCalc%getShiftPerL(shellPot(:,:,1))
    potential%intAtom(:,1) = potential%intAtom(:,1) + atomPot(:,1)
    potential%intShell(:,:,1) = potential%intShell(:,:,1) + shellPot(:,:,1)

    ! Build spin contribution (if necessary)
    if (this%tSpinPol) then
       call getSpinShift(shellPot, chargePerShell, this%species, orb, W)
       potential%intShell = potential%intShell + shellPot
    end if

    call total_shift(potential%intShell, potential%intAtom, orb, this%species)
    call total_shift(potential%intBlock, potential%intShell, orb, this%species)

    !! Add time dependent field if necessary
    if (this%tLaser) then
       do iAtom = 1, this%nAtom
          potential%extAtom(iAtom, 1) = dot_product(coord(:,iAtom), &
               & this%TDFunction(iStep, :))
       end do
       call total_shift(potential%extShell, potential%extAtom, orb, this%species)
       call total_shift(potential%extBlock, potential%extShell, orb, this%species)
       potential%intShell = potential%intShell + potential%extShell ! for SCC
       potential%intBlock = potential%intBlock + potential%extBlock ! for forces
    end if

    call add_shift(ham, over, nNeighbor, iNeighbor, &
         & this%species, orb, iPair, this%nAtom, img2CentCell, potential%intShell)

    if (this%nSpin == 2) then
       ham(:,:) = 2.0_dp * ham
       call qm2ud(ham)
       call qm2ud(q0)
    end if

    do iSpin=1,this%nSpin
       call unpackHS(T2,ham(:,iSpin),iNeighbor,nNeighbor,iSquare,iPair,img2CentCell)
       call blockSymmetrizeHS(T2,iSquare)
       H1(:,:,iSpin) = cmplx(T2, 0.0_dp, dp)
    end do

  end subroutine updateH


  !> Clear potential object
  subroutine clearPotential(potential)

    !> potential acting on the system
    type(TPotentials), intent(inout) :: potential

    potential%intAtom = 0.0_dp
    potential%intShell = 0.0_dp
    potential%intBlock = 0.0_dp
    potential%extAtom = 0.0_dp
    potential%extShell = 0.0_dp
    potential%extBlock = 0.0_dp
  end subroutine clearPotential


  !> Kick the density matrix for spectrum calculations
  subroutine kickDM(this, Rho, Ssqr, Sinv, iSquare, coord)
    !> ElecDynamics instance
    type(TElecDynamics), intent(in) :: this

    !> Square overlap
    complex(cp), intent(in) :: Ssqr(:,:)

    !> Square overlap inverse
    complex(cp), intent(in) :: Sinv(:,:)

    !> Density matrix
    complex(cp), intent(inout) :: Rho(:,:,:)

    !> atomic coordinates
    real(dp), allocatable, intent(in) :: coord(:,:)

    !> Index array for start of atomic block in dense matrices
    integer, intent(in) :: iSquare(:)

    complex(cp) :: T1(this%nOrbs, this%nOrbs, this%nSpin), T3(this%nOrbs, this%nOrbs, this%nSpin)
    complex(cp) :: T2(this%nOrbs, this%nOrbs), T4(this%nOrbs, this%nOrbs)
    integer :: iAt, iStart, iEnd, iSpin, iOrb
    real(dp) :: pkick(this%nSpin)

    pkick(1) = this%Field ! check units

    if (this%nSpin == 2 .and. this%SpType == iTDSinglet) pkick(2) = pkick(1)
    if (this%nSpin == 2 .and. this%SpType == iTDTriplet) pkick(2) = - pkick(1)

    T1 = 0.0_dp
    T2 = 0.0_dp
    T3 = 0.0_dp
    T4 = 0.0_dp

    do iSpin=1,this%nSpin
       do iAt = 1,this%nAtom
          iStart = iSquare(iAt)
          iEnd = iSquare(iAt+1)-1
          do iOrb = iStart, iEnd
             T1(iOrb, iOrb, iSpin) = Exp(cmplx(0, -pkick(iSpin)*coord(this%currPolDir, iAt), cp))
             T3(iOrb, iOrb, iSpin) = Exp(cmplx(0, pkick(iSpin)*coord(this%currPolDir, iAt), cp))
          end do
       end do
    end do

    do iSpin=1,this%nSpin
       call gemm(T2, T1(:,:,iSpin), Rho(:,:,iSpin), cmplx(1, 0, cp), cmplx(0, 0, cp))
       call gemm(T4, T2, Ssqr, cmplx(1, 0, cp), cmplx(0, 0, cp))
       call gemm(T2, T4, T3(:,:,iSpin), cmplx(1, 0, cp), cmplx(0, 0, cp))
       call gemm(Rho(:,:,iSpin), T2, Sinv, cmplx(0.5, 0, cp), cmplx(0, 0, cp), 'N', 'N')
       call gemm(Rho(:,:,iSpin), Sinv, T2, cmplx(0.5, 0, cp), cmplx(1, 0, cp), 'N', 'C')
    end do

    write(stdout,"(A)")'Density kicked!'

  end subroutine kickDM


  !> Creates array for external TD field
  subroutine getTDFunction(this)

    !> ElecDynamics instance
    type(TElecDynamics), intent(inout) :: this

    real(dp) :: midPulse, deltaT, angFreq, E0, time, envelope
    real(dp) :: tdfun(3)
    integer :: iStep, laserDat

    midPulse = (this%Time0 + this%Time1)/2.0_dp
    deltaT = this%Time1 - this%Time0
    angFreq = this%Omega
    E0 = this%Field
    this%TDFunction(:,:) = 0.0_dp

    open(newunit=laserDat, file='laser.dat')

    do iStep = 0,this%Nsteps
       time = iStep * this%Dt

       if (this%EnvType /= iTDConstant) then
          if (this%EnvType == iTDGaussian) then
             envelope = exp(-4.0_dp*pi*(time-midPulse)**2 / deltaT**2)
          else if (this%EnvType == iTDSin2 .and. (time >= this%Time0) .and. &
               & (time <= this%Time1)) then
             envelope = sin(pi*(time-this%Time0)/deltaT)**2
          end if
       else if (this%EnvType == iTDConstant) then
          envelope = 1.0_dp
       else
          envelope = 0.0_dp
       end if

       if (this%EnvType /= iTDFromFile) then
          this%TDFunction(iStep, :) = E0 * envelope * aimag(Exp(imag*(time*angFreq + this%phase)) &
               & * this%FieldDir)
       end if

       if (this%tEnvFromFile) then
          read(laserDat, *)time, tdfun(1), tdfun(2), tdfun(3)
          this%TDFunction(iStep, :) = this%TDFunction(iStep, :) + tdfun * (Bohr__AA / Hartree__eV)
       else
          write(laserDat, "(5F15.8)") time * au__fs, &
               & this%TDFunction(iStep, :) * (Hartree__eV / Bohr__AA)
       end if

    end do

    close(laserDat)
  end subroutine getTDFunction


  !> Calculate charges, dipole moments
  subroutine getChargeDipole(this, deltaQ, qq, dipole, q0, Rho, Ssqr, coord, iSquare)

    !> ElecDynamics instance
    type(TElecDynamics), intent(in) :: this

    !> Positive gross charge
    real(dp), intent(out) :: deltaQ(:,:)

    !> Dipole moment
    real(dp), intent(out) :: dipole(:,:)

    !> atomic ocupations
    real(dp), intent(out) :: qq(:,:,:)

    !> reference atomic occupations
    real(dp), intent(in) :: q0(:,:,:)

    !> atomic coordinates
    real(dp), intent(in) :: coord(:,:)

    !> Density matrix
    complex(cp), intent(in) :: Rho(:,:,:)

    !> Square overlap matrix
    complex(cp), intent(in) :: Ssqr(:,:)

    !> Index array for start of atomic block in dense matrices
    integer, intent(in) :: iSquare(:)

    integer :: iAt, iOrb, iSpin, iOrb2

    qq = 0.0_dp
    dipole = 0.0_dp

    do iSpin=1, this%nSpin
       do iAt = 1, this%nAtom
          iOrb = 0
          do iOrb2 = iSquare(iAt), iSquare(iAt+1)-1
             iOrb = iOrb + 1
             qq(iOrb,iAt,iSpin) = sum(Rho(iOrb2, :, iSpin) * Ssqr(iOrb2, :)) ! 1 <= iOrb <= orb%mOrb
          end do
          dipole(1,iSpin) = dipole(1,iSpin) + sum(q0(:, iAt, iSpin) - qq(:, iAt, iSpin))&
               & * coord(1, iAt)
          dipole(2,iSpin) = dipole(2,iSpin) + sum(q0(:, iAt, iSpin)-qq(:,iAt,iSpin))&
               & * coord(2,iAt)
          dipole(3,iSpin) = dipole(3,iSpin) + sum(q0(:,iAt,iSpin)-qq(:,iAt,iSpin))&
               & * coord(3,iAt)
       end do
    end do

    deltaQ(:,:) = sum((qq(:,:,:)-q0(:,:,:)), dim=1)
  end subroutine getChargeDipole


  !> Calculate energy
  subroutine getTDEnergy(this, energy, rhoPrim, Rho, iNeighbor, &
       & nNeighbor, orb, iSquare, iPair, img2CentCell, ham0, qq, q0, &
       & potential, chargePerShell, coord, pRepCont)

    !> ElecDynamics instance
    type(TElecDynamics), intent(inout) :: this

    !> data type for energy components and total
    type(TEnergies), intent(inout) :: energy

    !> sparse density matrix
    real(dp), allocatable, intent(inout) :: rhoPrim(:,:)

    !> Density matrix
    complex(cp), intent(in) :: Rho(:,:,:)

    !> Sparse storage for non-SCC hamitonian
    real(dp), intent(in) :: ham0(:)

    !> atomic ocupations
    real(dp), intent(inout) :: qq(:,:,:)

    !> reference atomic occupations
    real(dp), intent(in) :: q0(:,:,:)

    !> Atomic orbital information
    type(TOrbitals), intent(in) :: orb

    !> Atomic neighbour data
    integer, intent(in) :: iNeighbor(0:,:)

    !> Number of neighbours for each of the atoms
    integer, intent(in) :: nNeighbor(:)

    !> Index array for start of atomic block in dense matrices
    integer, intent(in) :: iSquare(:)

    !> index array for location of atomic blocks in large sparse arrays
    integer, intent(in) :: iPair(0:,:)

    !> image atoms to their equivalent in the central cell
    integer, intent(in) :: img2CentCell(:)

    !> electrons in each atomi shell
    real(dp), intent(in) :: chargePerShell(:,:,:)

    !> potential acting on the system
    type(TPotentials), intent(in) :: potential

    !> atomic coordinates
    real(dp), intent(in) :: coord(:,:)

    !> repulsive information
    type(ORepCont), intent(in) :: pRepCont

    integer :: iSpin

    iSpin = 1

    rhoPrim(:,iSpin) = 0.0_dp
    call packHS(rhoPrim(:,iSpin), real(Rho(:,:,iSpin), dp), iNeighbor, & ! should be complex
         &nNeighbor, orb%mOrb, iSquare, iPair, img2CentCell)

    energy%ETotal = 0.0_dp
    energy%EnonSCC = 0.0_dp
    energy%atomNonSCC(:) = 0.0_dp
    call mulliken(energy%atomNonSCC(:), rhoPrim(:,1), ham0, orb,&
         &iNeighbor, nNeighbor, img2CentCell, iPair)
    energy%EnonSCC =  sum(energy%atomNonSCC)

    if (this%tLaser) then ! energy in external field
       energy%atomExt = -sum(q0(:, :, 1) - qq(:, :, 1),dim=1) &
            & * potential%extAtom(:,1)
       energy%Eext =  sum(energy%atomExt)
    else
       energy%Eext = 0.0_dp
       energy%atomExt = 0.0_dp
    end if

    call this%sccCalc%getEnergyPerAtom(energy%atomSCC)
    energy%eSCC = sum(energy%atomSCC)

    if (this%nSpin > 1) then
       energy%atomSpin = 0.5_dp * sum(sum(potential%intShell(:,:,2:this%nSpin)&
            & * chargePerShell(:,:,2:this%nSpin), dim=1),dim=2)
       energy%Espin = sum(energy%atomSpin)
    else
       energy%atomSpin = 0.0_dp
       energy%eSpin = 0.0_dp
    end if

    !! Calculate repulsive energy
    call getERep(energy%atomRep, coord, nNeighbor, iNeighbor, &
         &this%species, pRepCont, img2CentCell)
    energy%Erep = sum(energy%atomRep)

    energy%Eelec = energy%EnonSCC + energy%eSCC + energy%Espin + energy%Eext
    energy%Etotal = energy%Eelec + energy%Erep + energy%eDisp

  end subroutine getTDEnergy


  !> Create all necessary matrices and instances for dynamics
  subroutine createMatrices(this, Rho, H1, Ssqr, Sinv, H0, ham0, &
       & over, ham, Hsq, filling, orb, rhoPrim, potential, &
       & iNeighbor, nNeighbor, iSquare, iPair, img2CentCell, Eiginv, &
       & EiginvAdj, energy)

    !> ElecDynamics instance
    type(TElecDynamics), intent(inout) :: this

    !> Eigenvectors
    real(dp), intent(inout) :: Hsq(:,:,:)

    !> overlap (sparse)
    real(dp), allocatable, intent(in) :: over(:)

    !> resulting hamitonian (sparse)
    real(dp), allocatable, intent(in) :: ham(:,:)

    !> occupations
    real(dp), intent(in) :: filling(:,:,:)

    !> Sparse storage for non-SCC hamitonian
    real(dp), intent(in) :: H0(:)

    !> Atomic neighbour data
    integer, intent(in) :: iNeighbor(0:,:)

    !> Number of neighbours for each of the atoms
    integer, intent(in) :: nNeighbor(:)

    !> Index array for start of atomic block in dense matrices
    integer, intent(in) :: iSquare(:)

    !> index array for location of atomic blocks in large sparse arrays
    integer, intent(in) :: iPair(0:,:)

    !> image atoms to their equivalent in the central cell
    integer, intent(in) :: img2CentCell(:)

    !> Atomic orbital information
    type(TOrbitals), intent(in) :: orb

    !> potential acting on the system
    type(TPotentials), intent(out) :: potential

    !> data type for energy components and total
    type(TEnergies), intent(out) :: energy

    !> sparse density matrix
    real(dp), allocatable, intent(out) :: rhoPrim(:,:)

    !> Sparse storage for non-SCC hamitonian
    real(dp), allocatable, intent(out) :: ham0(:)

    !> Square overlap matrix
    complex(cp), intent(out) :: Ssqr(:,:)

    !> Square overlap inverse
    complex(cp), intent(out) :: Sinv(:,:)

    !> Square hamiltonian
    complex(cp), intent(out) :: H1(:,:,:)

    !> Density matrix
    complex(cp), intent(out) :: Rho(:,:,:)

    !> Inverse of eigenvectors matrix (for populations)
    complex(cp), allocatable, intent(out), optional :: Eiginv(:,:,:)

    !> Adjoint of the inverse of eigenvectors matrix (for populations)
    complex(cp), allocatable, intent(out), optional :: EiginvAdj(:,:,:)

    real(dp) :: T2(this%nOrbs,this%nOrbs), T3(this%nOrbs, this%nOrbs)
    integer :: iSpin, iOrb, iOrb2

    allocate(rhoPrim(size(ham), this%nSpin))
    allocate(ham0(size(H0)))
    ham0(:) = H0

    T2 = 0.0_dp
    T3 = 0.0_dp
    call unpackHS(T2,over,iNeighbor,nNeighbor,iSquare,iPair,img2CentCell)
    call blockSymmetrizeHS(T2,iSquare)
    Ssqr(:,:) = cmplx(T2, 0, cp)

    do iSpin=1,this%nSpin
       call unpackHS(T3,ham(:,iSpin),iNeighbor,nNeighbor,iSquare,iPair,img2CentCell)
       call blockSymmetrizeHS(T3,iSquare)
       H1(:,:,iSpin) = cmplx(T3, 0, cp)
       T3 = 0.0_dp
    end do

    if (this%tPopulations) then
       allocate(Eiginv(this%nOrbs, this%nOrbs, this%nSpin))
       allocate(EiginvAdj(this%nOrbs, this%nOrbs, this%nSpin))
       do iSpin=1,this%nSpin
          call tdPopulInit(this, Eiginv, EiginvAdj, HSq, iSpin)
       end do
    end if

    do iOrb = 1, this%nOrbs
       T3(iOrb, iOrb) = 1.0_dp
    end do
    call gesv(T2, T3)
    Sinv(:,:) = cmplx(T3, 0, cp)
    write(stdOut,"(A)")'S inverted'

    do iSpin=1,this%nSpin
       T2 = 0.0_dp
       call makeDensityMatrix(T2,Hsq(:,:,iSpin),filling(:,1,iSpin))
       Rho(:,:,iSpin) = cmplx(T2, 0, cp)
       do iOrb = 1, this%nOrbs-1
          do iOrb2 = iOrb+1, this%nOrbs
             Rho(iOrb, iOrb2, iSpin) = Rho(iOrb2, iOrb, iSpin)
          end do
       end do
    end do

    call init(potential, orb, this%nAtom, this%nSpin)
    call init(energy, this%nAtom)
  end subroutine createMatrices


  !> Perfoms a step backwards to boot the dynamics using the Euler algorithm
  subroutine initializePropagator(this, step, Rhoold, Rho, H1, Sinv)
    !> ElecDynamics instance
    type(TElecDynamics), intent(in) :: this

    !> Density matrix
    complex(cp), intent(in) :: Rho(:,:,:)

    !> Square overlap inverse
    complex(cp), intent(in) :: Sinv(:,:)

    !> Square hamiltonian
    complex(cp), intent(inout) :: H1(:,:,:)

    !> Density matrix at previous step
    complex(cp), intent(out) :: Rhoold(:,:,:)

    !> Time step in atomic units (with sign, to perform step backwards or forwards)
    real(dp), intent(in) :: step

    complex(cp) :: T1(this%nOrbs,this%nOrbs)
    integer :: iSpin

    Rhoold(:,:,:) = Rho

    do iSpin=1,this%nSpin
       T1 = 0.0_cp
       H1(:,:,iSpin) = imag * H1(:,:,iSpin)
       call propagateRho(this, Rhoold(:,:,iSpin), Rho(:,:,iSpin), &
            & H1(:,:,iSpin), Sinv, T1, step)
    end do

  end subroutine initializePropagator


  !> Propagate Rho, notice that H = iH (coeficients are real)
  subroutine propagateRho(this, Rhoold, Rho, H1, Sinv, T1, step)

    !> ElecDynamics instance
    type(TElecDynamics), intent(in) :: this

    !> Density matrix at previous step
    complex(cp), intent(inout) :: Rhoold(:,:)

    !> Density matrix
    complex(cp), intent(in) :: Rho(:,:)

    !> Square hamiltonian
    complex(cp), intent(in) :: H1(:,:)

    !> Square overlap inverse
    complex(cp), intent(in) :: Sinv(:,:)

    !> Auxiliary matrix
    complex(cp), intent(out) :: T1(:,:)

    !> Time step in atomic units
    real(dp), intent(in) :: step

    call gemm(T1, Sinv, H1, cmplx(1, 0, cp), cmplx(0, 0, cp), &
         & 'N', 'N', this%nOrbs, this%nOrbs, this%nOrbs)
    call gemm(Rhoold, T1, Rho, cmplx(-step, 0, cp),&
         & cmplx(1, 0, cp), 'N', 'N', this%nOrbs, this%nOrbs, this%nOrbs)
    call gemm(Rhoold, Rho, T1, cmplx(-step, 0, cp),&
         & cmplx(1, 0, cp), 'N', 'C', this%nOrbs, this%nOrbs, this%nOrbs)
  end subroutine propagateRho


  !> Initialize output files
  subroutine initTDOutput(this, dipoleDat, qDat, energyDat, populDat)

    !> ElecDynamics instance
    type(TElecDynamics), intent(in) :: this

    !> Dipole output file ID
    integer, intent(out) :: dipoleDat

    !> Charge output file ID
    integer, intent(out) :: qDat

    !> Energy output file ID
    integer, intent(out) :: energyDat

    !> Populations  output file ID
    integer, intent(out) :: populDat(2)

    character(20) :: dipoleFileName
    character(1) :: strSpin
    integer :: iSpin
    logical :: exist

    if (this%tKick) then
       if (this%currPolDir == 1) then
          dipoleFileName = 'mux.dat'
       else if (this%currPolDir == 2) then
          dipoleFileName = 'muy.dat'
       else if (this%currPolDir == 3) then
          dipoleFileName = 'muz.dat'
       end if
    else
       dipoleFileName = 'mu.dat'
    end if

    call openFile(this, dipoleDat, dipoleFileName)
    call openFile(this, qDat, 'qsvst.dat')
    call openFile(this, energyDat, 'energyvst.dat')

    if (this%tPopulations) then
       do iSpin=1,this%nSpin
          write(strSpin,'(i1)')iSpin
          call openFile(this, populDat(iSpin), 'molPopul' // trim(strSpin) // '.dat')
       end do
    end if

  end subroutine initTDOutput


  !> Close output files
  subroutine closeTDOutputs(this, dipoleDat, qDat, energyDat, populDat)

    !> ElecDynamics instance
    type(TElecDynamics), intent(in) :: this

    !> Dipole output file ID
    integer, intent(out) :: dipoleDat

    !> Charge output file ID
    integer, intent(out) :: qDat

    !> Energy output file ID
    integer, intent(out) :: energyDat

    !> Populations output file ID
    integer, intent(in) :: populDat(2)

    integer :: iSpin

    close(dipoleDat)
    close(qDat)
    close(energyDat)

    if (this%tPopulations) then
       do iSpin = 1, this%nSpin
          close(populDat(iSpin))
       end do
    end if

    if (this%tWriteAutotest) then
       close(this%fdAutotest)
    end if
  end subroutine closeTDOutputs


  !> Open files in different ways deppending on their previous existance
  subroutine openFile(this, unitName, fileName)

    !> ElecDynamics instance
    type(TElecDynamics), intent(in) :: this

    !> File ID
    integer :: unitName

    !> Name of the file to open
    character(*) :: fileName

    logical :: exist=.false.

    if (this%tRestart) then
       inquire(file=fileName, exist=exist)
    end if

    if (exist) then
       open(newunit=unitName, file=fileName, status="old", position="append", action="write")
    else
       open(newunit=unitName, file=fileName, action="write")
    end if
  end subroutine openFile


  !> Write to and read from restart files
  subroutine writeRestart(Rho, Ssqr, coord, time, dumpName)

    complex(cp), intent(in) :: Rho(:,:,:), Ssqr(:,:)
    !> atomic coordinates
    real(dp), intent(in) :: coord(:,:)

    !> elapsed simulated time in atomic units
    real(dp), intent(in) :: time

    !> name of the dump file
    character(len=*), intent(in), optional :: dumpName

    integer :: dumpBin

    if (present(dumpName)) then
       open(newunit=dumpBin, file=dumpName, form='unformatted', access='stream', &
            & action='write')
    else
       open(newunit=dumpBin, file='tddump.bin', form='unformatted', access='stream', &
            & action='write')
    end if
    write(dumpBin) Rho, Ssqr, coord, time
    close(dumpBin)
  end subroutine writeRestart


  subroutine readRestart(Rho, Ssqr, coord, time)

    !> Density Matrix
    complex(cp), intent(out) :: Rho(:,:,:)

    !> Square overlap matrix
    complex(cp), intent(out) :: Ssqr(:,:)

    !> atomic coordinates
    real(dp), intent(out) :: coord(:,:)

    !> Previous simulation elapsed time until restart file writing
    real(dp), intent(out) :: time

    integer :: dumpBin

    open(newunit=dumpBin, file='tddump.bin', form='unformatted', access='stream', &
         & action='read')
    read(dumpBin) Rho, Ssqr, coord, time
    close(dumpBin)
  end subroutine readRestart


  !> Write results to file
  subroutine writeTDOutputs(this, dipoleDat, qDat, energyDat, &
       & time, energy, dipole, deltaQ, iStep)

    !> ElecDynamics instance
    type(TElecDynamics), intent(in) :: this

    !> data type for energy components and total
    type(TEnergies), intent(in) :: energy

    !> Dipole output file ID
    integer, intent(out) :: dipoleDat

    !> Charge output file ID
    integer, intent(out) :: qDat

    !> Energy output file ID
    integer, intent(out) :: energyDat

    !> Elapsed simulation time
    real(dp), intent(in) :: time

    !> Dipole moment
    real(dp), intent(in) :: dipole(:,:)

    !> Positive gross charge
    real(dp), intent(in) :: deltaQ(:,:)

    !> current step of the propagation
    integer, intent(in) :: iStep

    integer :: iAtom, iSpin, iDir

    write(energydat, '(9F25.15)') time * au__fs, energy%Etotal, energy%EnonSCC, &
         & energy%eSCC, energy%Espin, energy%Eext, energy%Erep

    write(dipoleDat, '(7F25.15)') time * au__fs, ((dipole(iDir, iSpin) * Bohr__AA, iDir=1, 3),&
         & iSpin=1, this%nSpin)

    if (mod(iStep, this%writeFreq) == 0) then
       write(qDat, '(*(2x,F25.15))') time * au__fs, sum(deltaQ),&
            & (sum(deltaQ(iAtom,:)), iAtom=1, this%nAtom)
    end if

  end subroutine writeTDOutputs


  !> Initialize matrices for populations
  subroutine tdPopulInit(this, Eiginv, EiginvAdj, HSq, iSpin)
    !> ElecDynamics instance
    type(TElecDynamics), intent(in) :: this

    !> Inverse of eigenvectors matrix (for populations)
    complex(cp), intent(out) :: Eiginv(:,:,:)

    !> Adjoint of the inverse of eigenvectors matrix (for populations)
    complex(cp), intent(out) :: EiginvAdj(:,:,:)

    !> Eigenvectors
    real(dp), intent(in) :: Hsq(:,:,:)

    !> Spin index
    integer, intent(in) :: iSpin

    real(dp), allocatable :: T2(:,:), T3(:,:)
    integer :: iOrb, iOrb2

    allocate(T2(this%nOrbs, this%nOrbs), T3(this%nOrbs, this%nOrbs))

    T2 = Hsq(:,:,iSpin)
    T3 = 0.0_dp
    do iOrb = 1, this%nOrbs
       T3(iOrb, iOrb) = 1.0_dp
    end do
    call gesv(T2,T3)
    Eiginv(:, :, iSpin) = cmplx(T3, 0, cp)

    T3 = 0.0_dp
    do iOrb = 1, this%nOrbs
       T3(iOrb, iOrb) = 1.0_dp
       do iOrb2 = 1, this%nOrbs
          T2(iOrb, iOrb2) = Hsq(iOrb2, iOrb, iSpin)
       end do
    end do
    call gesv(T2,T3)
    EiginvAdj(:, :, iSpin) = cmplx(T3, 0, cp)

    deallocate(T2, T3)

  end subroutine tdPopulInit


  !> Calculate populations at each time step
  subroutine getTDPopulations(this, Rho, Eiginv, EiginvAdj, T1, &
       & populDat, time, iSpin)

    !> ElecDynamics instance
    type(TElecDynamics), intent(in) :: this

    !> Density Matrix
    complex(cp), intent(in) :: Rho(:,:,:)
    !> Inverse of eigenvectors matrix (for populations)

    complex(cp), intent(inout) :: Eiginv(:,:,:)

    !> Adjoint of the inverse of eigenvectors matrix (for populations)
    complex(cp), intent(inout) :: EiginvAdj(:,:,:)

    !> Auxiliary matrix
    complex(cp), intent(inout) :: T1(:,:)

    !> Elapsed simulation time
    real(dp), intent(in) :: time

    !> Populations output file ID
    integer, intent(in) :: populDat(2)

    !> Spin index
    integer, intent(in) :: iSpin

    complex(cp) :: T2(this%nOrbs, this%nOrbs)
    real(dp) :: T3(this%nOrbs, this%nOrbs)
    integer :: iOrb

    T3 = 0.0_dp

    call gemm(T1, Eiginv(:,:,iSpin), Rho(:,:,iSpin))
    call gemm(T2, T1, EiginvAdj(:,:,iSpin))

    write(populDat(iSpin),'(*(2x,F25.15))') time * au__fs, &
         & (real(T2(iOrb, iOrb), dp), iOrb=1, this%nOrbs)

  end subroutine getTDPopulations

  subroutine writeTDAutotest(this, dipole, energy, deltaQ)
    !> ElecDynamics instance
    type(TElecDynamics), intent(in) :: this

    !> Dipole moment
    real(dp), intent(in) :: dipole(:,:)

    !> data type for energy components and total
    type(TEnergies), intent(in) :: energy

    !> Positive gross charge
    real(dp), intent(in) :: deltaQ(:,:)

    call writeTagged(this%fdAutotest, tag_tdenergy, energy%eSCC)
    call writeTagged(this%fdAutotest, tag_tddipole, dipole)
    call writeTagged(this%fdAutotest, tag_tdcharges, deltaQ)
  end subroutine writeTDAutotest

end module timeprop_module
