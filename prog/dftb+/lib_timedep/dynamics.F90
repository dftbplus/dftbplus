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

  implicit none

  type TElecDynamicsInp
     real(dp) :: tdField, Dt, Omega, ReFieldPolVec(3), ImFieldPolVec(3)
     real(dp) :: Time0, Time1
     integer :: PolDir, Steps, SaveEvery, PertType, EnvType, SpType, restartFreq
     logical :: Populations, Restart, WriteRestart
     real(dp) :: Phase
  end type TElecDynamicsInp

  type TElecDynamics
     private
     real(dp) :: Field, Dt, Omega, Time0, Time1
     complex(cp) :: FieldDir(3)
     real(dp), allocatable :: TDFunction(:, :), phase
     integer :: Nsteps, SaveEvery, PertType, EnvType, SpType
     integer :: nAtom, nOrbs, nSpin=1, currPolDir=1, restartEvery
     integer, allocatable :: species(:), PolDirs(:)
     character(mc), allocatable :: speciesName(:)
     logical :: Populations, SpinPol=.false.
     logical :: Restart, WriteRestart
     logical :: doLaser = .false., doKick = .false., EnvFromFile = .false.
     type(TScc), allocatable :: sccCalc
  end type TElecDynamics

  public :: runDynamics, TElecDynamics_init
  public :: TElecDynamicsInp, TElecDynamics
  public :: iKick, iLaser, iNoTDPert
  public :: iTDConstant, iTDGaussian, iTDSin2, iTDFromFile
  public :: iTDSinglet, iTDTriplet

  !! Enumerating available types of perturbation
  integer, parameter :: iKick = 1, iLaser = 2, iNoTDPert = 3

  !! Enumerating available types of envelope function
  integer, parameter :: iTDConstant = 1, iTDGaussian = 2, iTDSin2 = 3, iTDFromFile = 4

  !! Enumerating available types of spin polarized spectra
  integer, parameter :: iTDSinglet = 1, iTDTriplet = 2


contains

  !! This initializes the input variables
  subroutine TElecDynamics_init(this, inp, species, speciesName)
    type(TElecDynamics), intent(out) :: this
    type(TElecDynamicsInp), intent(in) :: inp
    integer, allocatable, intent(in) :: species(:)
    character(mc), allocatable, intent(in) :: speciesName(:)

    real(dp) :: norm

    this%Field = inp%tdField
    this%Dt = inp%Dt
    this%Nsteps = inp%Steps
    this%PertType = inp%PertType
    this%EnvType = inp%EnvType
    this%SpType = inp%SpType
    this%Populations = inp%Populations
    this%Restart = inp%Restart
    this%WriteRestart = inp%WriteRestart
    this%phase = inp%Phase
    this%SaveEvery = inp%SaveEvery
    this%restartEvery = int(inp%Steps/inp%restartFreq)
    allocate(this%species, source=species)
    allocate(this%speciesName, source=speciesName)
    allocate(this%sccCalc)

    if (inp%EnvType /= iTDConstant) then
       this%Time0 = inp%Time0
       this%Time1 = inp%Time1
    end if

    if (inp%PertType == iLaser) then
       this%doLaser = .true.
    else if (inp%PertType == iKick) then
       this%doKick = .true.
    end if

    if (this%doLaser) then
       this%Omega = inp%Omega
       this%FieldDir = inp%ReFieldPolVec + imag * inp%ImFieldPolVec
       norm = sqrt(sum(this%FieldDir(:)*conjg(this%FieldDir(:))))
       this%FieldDir = this%FieldDir/norm ! normalize polarization vector
       allocate(this%TDFunction(0:this%Nsteps, 3))
       if (this%EnvType == iTDFromFile) then
          this%EnvFromFile = .true.
       end if
       call getTDFunction(this)
    end if

    if (this%doKick) then
       if (inp%PolDir == 4) then
          allocate(this%PolDirs(3))
          this%PolDirs(:) = (/ 1, 2, 3 /)
       else
          allocate(this%PolDirs(1))
          this%PolDirs(1) = inp%PolDir
       end if
    end if
  end subroutine TElecDynamics_init


  !! Driver of calculation in order to perform a spectrum
  subroutine runDynamics(this, Hsq, ham, H0, q0, over, filling, neighborList, nNeighbor,&
       &iSquare, iPair, img2CentCell, orb, coord, W, pRepCont, sccCalc, env)
    type(TElecDynamics) :: this
    real(dp), intent(inout) :: Hsq(:,:,:), H0(:), q0(:,:,:)
    real(dp), allocatable, intent(inout) :: ham(:,:), over(:), coord(:,:)
    real(dp), intent(in) :: W(:,:,:), filling(:,:,:)
    integer, intent(inout) :: nNeighbor(:)
    integer, allocatable, intent(inout) :: iPair(:,:), img2CentCell(:)
    integer, intent(in) :: iSquare(:)
    type(TNeighborList), intent(inout) :: neighborList
    type(ORepCont), intent(in) :: pRepCont
    type(TOrbitals), intent(in) :: orb
    type(TScc), intent(in) :: sccCalc
    type(TEnvironment), intent(inout) :: env
    integer :: iPol

    this%sccCalc = sccCalc

    this%nSpin = size(ham(:,:), dim=2)
    if (this%nSpin > 1) then
       this%SpinPol = .true.
       call qm2ud(q0)
    end if

    this%nOrbs = size(Hsq, dim=1)
    this%nAtom = size(coord, dim=2)

    if (this%doKick) then
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


  !! Runs the electronic dynamics of the system
  subroutine doDynamics(this,Hsq,ham,H0,q0,over,filling,neighborList,nNeighbor, &
       &iSquare,iPair,img2CentCell,orb,coord,W,pRepCont,env)
    type(TElecDynamics) :: this
    real(dp), intent(inout) :: Hsq(:,:,:), H0(:), q0(:,:,:)
    real(dp), allocatable, intent(inout) :: ham(:,:), over(:), coord(:,:)
    real(dp), intent(in) :: W(:,:,:), filling(:,:,:)
    integer, intent(inout) :: nNeighbor(:)
    integer, allocatable, intent(inout) :: iPair(:,:), img2CentCell(:)
    integer, intent(in) :: iSquare(:)
    type(TNeighborList), intent(inout) :: neighborList
    type(ORepCont), intent(in) :: pRepCont
    type(TOrbitals), intent(in) :: orb
    type(TEnvironment), intent(inout) :: env

    complex(cp) :: Ssqr(this%nOrbs,this%nOrbs), Sinv(this%nOrbs,this%nOrbs), T1(this%nOrbs,this%nOrbs)
    complex(cp) :: Rho(this%nOrbs,this%nOrbs,this%nSpin), Rhoold(this%nOrbs,this%nOrbs,this%nSpin)
    complex(cp) :: H1(this%nOrbs,this%nOrbs,this%nSpin)
    complex(cp) :: Rhonew(this%nOrbs,this%nOrbs,this%nSpin)
    complex(cp), allocatable :: Eiginv(:,:,:), EiginvAdj(:,:,:)
    real(dp) :: qq(orb%mOrb, this%nAtom, this%nSpin), deltaQ(this%nAtom,this%nSpin), dipole(3,this%nSpin)
    real(dp) :: chargePerShell(orb%mShell,this%nAtom,this%nSpin)
    real(dp), allocatable :: rhoPrim(:,:), ham0(:)
    real(dp) :: time, dTime, startTime = 0.0_dp, timeInit = 0.0_dp, timeElec = 0.0_dp
    integer :: dipoleDat, qDat, energyDat, populDat(2)
    integer :: iStep = 0, iAtom, iSpin, iTimeInit, iTimeElec, iTimeLoop
    type(TPotentials) :: potential
    type(TEnergies) :: energy
    character(4) :: dumpIdx

    ! Initialize timer
    call env%globalTimer%startTimer(globalTimers%elecDynInit)
    if (this%Restart) then
       call readRestart(Rho, Ssqr, coord, startTime)
    end if

    ! Initialize stuff
    call createMatrices(this, Rho, H1, Ssqr, Sinv, H0, ham0, &
         & over, ham, Hsq, filling, orb, rhoPrim, potential, &
         & neighborList%iNeighbor, nNeighbor, iSquare, iPair, img2CentCell, &
         & Eiginv, EiginvAdj, energy)
    ! Initialize output files
    call initTDOutput(this, dipoleDat, qDat, energyDat, populDat)

    call getChargeDipole(this, deltaQ, qq, dipole, q0, Rho, Ssqr, coord, iSquare) !q_0
    call updateH(this, H1, ham, over, ham0, qq, q0, coord, orb, potential, &  !H1_0
         &neighborList%iNeighbor, nNeighbor, iSquare, iPair, img2CentCell, iStep, chargePerShell, W, env)

    ! Apply kick to Rho if necessary
    if (this%doKick) then
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
    call tic(iTimeLoop)

    write(stdOut,"(A)")'Starting dynamics'

    do iStep = 0, this%Nsteps
       time = iStep * this%Dt + startTime

       if (.not. this%Restart .or. iStep > 0) then
          call writeTDOutputs(this, dipoleDat, qDat, energyDat, time, energy, dipole, deltaQ, iStep)
       end if

       call getChargeDipole(this, deltaQ, qq, dipole, q0, Rho, Ssqr, coord, iSquare)
       call updateH(this, H1, ham, over, ham0, qq, q0, coord, orb, &
            & potential, neighborList%iNeighbor, nNeighbor, iSquare, iPair, &
            & img2CentCell, iStep, chargePerShell, W, env)

       if ((this%WriteRestart) .and. (iStep > 0) .and. &
            & (mod(iStep, this%restartEvery) == 0)) then
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

          if ((this%Populations) .and. (mod(iStep, this%SaveEvery) == 0)) then
             call getTDPopulations(this, Rho, Eiginv, EiginvAdj, T1, populDat, time, iSpin)
          end if
       end do

       if (mod(iStep,this%Nsteps/10) == 0) then
          call tac(dTime,iTimeLoop)
          write(stdOut,"(A,2x,I6,2(2x,A,F10.6))") 'Step ', iStep, 'elapsed loop time: ', dTime, &
               & 'average time per loop ', dTime/(iStep+1)
       end if
    end do

    write(stdOut, "(A)") 'Dynamics finished OK!'
    call env%globalTimer%stopTimer(globalTimers%elecDynLoop)
    call closeTDOutputs(this, dipoleDat, qDat, energyDat, populDat)
  end subroutine doDynamics



  !! Applies SCC to hamiltonian and adds TD perturbation (if any)
  subroutine updateH(this, H1, ham, over, ham0, qq, q0, coord, orb, potential, &
       &iNeighbor, nNeighbor, iSquare, iPair, img2CentCell, iStep, chargePerShell, W, env)
    type(TElecDynamics), intent(inout) :: this
    complex(cp), intent(inout) :: H1(:,:,:)
    real(dp), allocatable, intent(inout) :: ham(:,:), over(:), coord(:,:)
    real(dp), intent(in) :: ham0(:)
    real(dp), intent(inout) :: qq(:,:,:), q0(:,:,:)
    real(dp), intent(in) :: W(:,:,:)
    type(TOrbitals), intent(in) :: orb
    integer, intent(in) :: iNeighbor(0:,:),nNeighbor(:)
    integer, intent(in) :: iSquare(:), iPair(0:,:), img2CentCell(:)
    integer, intent(in) :: iStep
    real(dp), intent(inout) :: chargePerShell(:,:,:)
    real(dp) :: T2(this%nOrbs,this%nOrbs)
    real(dp) :: shellPot(orb%mShell, this%nAtom, this%nSpin), atomPot(this%nAtom, this%nSpin)
    type(TPotentials), intent(inout) :: potential
    integer :: iAtom, iSpin
    type(TEnvironment), intent(inout) :: env

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
    if (this%SpinPol) then
       call getSpinShift(shellPot, chargePerShell, this%species, orb, W)
       potential%intShell = potential%intShell + shellPot
    end if

    call total_shift(potential%intShell, potential%intAtom, orb, this%species)
    call total_shift(potential%intBlock, potential%intShell, orb, this%species)

    !! Add time dependent field if necessary
    if (this%doLaser) then
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


  !! Clear potential object
  subroutine clearPotential(potential)
    type(TPotentials), intent(inout) :: potential

    potential%intAtom = 0.0_dp
    potential%intShell = 0.0_dp
    potential%intBlock = 0.0_dp
    potential%extAtom = 0.0_dp
    potential%extShell = 0.0_dp
    potential%extBlock = 0.0_dp
  end subroutine clearPotential


  !! Kick the density matrix for spectrum calculations
  subroutine kickDM(this, Rho, Ssqr, Sinv, iSquare, coord)
    type(TElecDynamics), intent(in) :: this
    complex(cp), intent(in) :: Ssqr(:,:), Sinv(:,:)
    complex(cp), intent(inout) :: Rho(:,:,:)
    real(dp), allocatable, intent(in) :: coord(:,:)
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


  !! Calculate envelope function
  subroutine getTDFunction(this)
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

       if (this%EnvFromFile) then
          read(laserDat, *)time, tdfun(1), tdfun(2), tdfun(3)
          this%TDFunction(iStep, :) = this%TDFunction(iStep, :) + tdfun * (Bohr__AA / Hartree__eV)
       else
          write(laserDat, "(5F15.8)") time * au__fs, &
               & this%TDFunction(iStep, :) * (Hartree__eV / Bohr__AA)
       end if

    end do

    close(laserDat)
  end subroutine getTDFunction


  !! Calculate charges, dipole moments
  subroutine getChargeDipole(this, deltaQ, qq, dipole, q0, Rho, Ssqr, coord, iSquare)
    type(TElecDynamics), intent(in) :: this
    real(dp), intent(out) :: qq(:,:,:), dipole(:,:), deltaQ(:,:)
    real(dp), intent(in) :: coord(:,:), q0(:,:,:)
    complex(cp), intent(in) :: Rho(:,:,:), Ssqr(:,:)
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


  !! Calculate energy
  subroutine getTDEnergy(this, energy, rhoPrim, Rho, iNeighbor, &
       & nNeighbor, orb, iSquare, iPair, img2CentCell, ham0, qq, q0, &
       & potential, chargePerShell, coord, pRepCont)
    type(TElecDynamics), intent(inout) :: this
    type(TEnergies), intent(inout) :: energy
    real(dp), allocatable, intent(inout) :: rhoPrim(:,:)
    complex(cp), intent(in) :: Rho(:,:,:)
    real(dp), intent(in) :: ham0(:)
    real(dp), intent(inout) :: qq(:,:,:)
    real(dp), intent(in) :: q0(:,:,:)
    type(TOrbitals), intent(in) :: orb
    integer, intent(in) :: iNeighbor(0:,:),nNeighbor(:)
    integer, intent(in) :: iSquare(:), iPair(0:,:), img2CentCell(:)
    real(dp), intent(in) :: chargePerShell(:,:,:)
    type(TPotentials), intent(in) :: potential
    integer :: iSpin
    real(dp), intent(in) :: coord(:,:)
    type(ORepCont), intent(in) :: pRepCont

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

    if (this%doLaser) then ! energy in external field
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


  !! Create all necessary matrices for dynamics
  subroutine createMatrices(this, Rho, H1, Ssqr, Sinv, H0, ham0, &
       & over, ham, Hsq, filling, orb, rhoPrim, potential, &
       & iNeighbor, nNeighbor, iSquare, iPair, img2CentCell, Eiginv, &
       & EiginvAdj, energy)
    type(TElecDynamics), intent(inout) :: this
    real(dp), intent(inout) :: Hsq(:,:,:)
    real(dp), allocatable, intent(in) :: over(:), ham(:,:)
    real(dp), intent(in) :: filling(:,:,:), H0(:)
    integer, intent(in) :: iNeighbor(0:,:),nNeighbor(:)
    integer, intent(in) :: iSquare(:), iPair(0:,:), img2CentCell(:)
    type(TOrbitals), intent(in) :: orb
    type(TPotentials), intent(out) :: potential
    type(TEnergies), intent(out) :: energy

    real(dp), allocatable, intent(out) :: rhoPrim(:,:), ham0(:)
    complex(cp), intent(out) :: Ssqr(:,:), Sinv(:,:), H1(:,:,:)
    complex(cp), intent(out) :: Rho(:,:,:)
    complex(cp), allocatable, intent(out), optional :: Eiginv(:,:,:), EiginvAdj(:,:,:)
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

    if (this%Populations) then
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


  !! Perfoms a step backwards to boot the dynamics using the Euler algorithm
  subroutine initializePropagator(this, step, Rhoold, Rho, H1, Sinv)
    type(TElecDynamics), intent(in) :: this
    complex(cp), intent(in) :: Rho(:,:,:), Sinv(:,:)
    complex(cp), intent(inout) :: H1(:,:,:) ! Warning! H1 is modified
    complex(cp), intent(out) :: Rhoold(:,:,:)
    real(dp), intent(in) :: step ! Step with its sign, to perform Euler backwards (-) or forwards (+)
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


  !! Propagate Rho, notice that H = iH (coeficients are real)
  subroutine propagateRho(this, Rhoold, Rho, H1, Sinv, T1, step)
    type(TElecDynamics), intent(in) :: this
    complex(cp), intent(inout) :: Rhoold(:,:)
    complex(cp), intent(in) :: Rho(:,:), H1(:,:), Sinv(:,:)
    complex(cp), intent(out) :: T1(:,:)
    real(dp), intent(in) :: step

    call gemm(T1, Sinv, H1, cmplx(1, 0, cp), cmplx(0, 0, cp), &
         & 'N', 'N', this%nOrbs, this%nOrbs, this%nOrbs)
    call gemm(Rhoold, T1, Rho, cmplx(-step, 0, cp),&
         & cmplx(1, 0, cp), 'N', 'N', this%nOrbs, this%nOrbs, this%nOrbs)
    call gemm(Rhoold, Rho, T1, cmplx(-step, 0, cp),&
         & cmplx(1, 0, cp), 'N', 'C', this%nOrbs, this%nOrbs, this%nOrbs)
  end subroutine propagateRho


  !! Initialize output files
  subroutine initTDOutput(this, dipoleDat, qDat, energyDat, populDat)
    type(TElecDynamics), intent(in) :: this
    integer, intent(out) :: dipoleDat, qDat, energyDat, populDat(2)
    character(20) :: dipoleFileName
    character(1) :: strSpin
    integer :: iSpin
    logical :: exist

    if (this%doKick) then
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

    if (this%Populations) then
       do iSpin=1,this%nSpin
          write(strSpin,'(i1)')iSpin
          call openFile(this, populDat(iSpin), 'molPopul' // trim(strSpin) // '.dat')
       end do
    end if

  end subroutine initTDOutput

  subroutine closeTDOutputs(this, dipoleDat, qDat, energyDat, populDat)
    type(TElecDynamics), intent(in) :: this
    integer, intent(in) :: dipoleDat, qDat, energyDat, populDat(2)
    integer :: iSpin

    close(dipoleDat)
    close(qDat)
    close(energyDat)

    if (this%Populations) then
       do iSpin = 1, this%nSpin
          close(populDat(iSpin))
       end do
    end if
  end subroutine closeTDOutputs


  !! Open files in different ways deppending on their previous existance
  subroutine openFile(this, unitName, fileName)
    type(TElecDynamics), intent(in) :: this
    logical :: exist=.false.
    integer :: unitName
    character(*) :: fileName

    if (this%Restart) then
       inquire(file=fileName, exist=exist)
    end if

    if (exist) then
       open(newunit=unitName, file=fileName, status="old", position="append", action="write")
    else
       open(newunit=unitName, file=fileName, action="write")
    end if
  end subroutine openFile


  !! Write to and read from restart files
  subroutine writeRestart(Rho, Ssqr, coord, time, dumpName)
    complex(cp), intent(in) :: Rho(:,:,:), Ssqr(:,:)
    real(dp), intent(in) :: coord(:,:), time
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
    complex(cp), intent(out) :: Rho(:,:,:), Ssqr(:,:)
    real(dp), intent(out) :: coord(:,:), time
    integer :: dumpBin

    open(newunit=dumpBin, file='tddump.bin', form='unformatted', access='stream', &
         & action='read')
    read(dumpBin) Rho, Ssqr, coord, time
    close(dumpBin)
  end subroutine readRestart


  !! Write results to file
  subroutine writeTDOutputs(this, dipoleDat, qDat, energyDat, &
       & time, energy, dipole, deltaQ, iStep)
    type(TElecDynamics), intent(in) :: this
    type(TEnergies), intent(in) :: energy
    integer, intent(in) :: dipoleDat, qDat, energyDat
    real(dp), intent(in) :: time, dipole(:,:), deltaQ(:,:)
    integer, intent(in) :: iStep
    integer :: iAtom, iSpin, iDir

    write(energydat, '(9F25.15)') time * au__fs, energy%Etotal, energy%EnonSCC, &
         & energy%eSCC, energy%Espin, energy%Eext, energy%Erep
    write(dipoleDat, '(7F25.15)') time * au__fs, ((dipole(iDir, iSpin) * Bohr__AA, iDir=1, 3),&
         & iSpin=1, this%nSpin)
    if (mod(iStep, this%SaveEvery) == 0) then
       write(qDat, '(*(2x,F25.15))') time * au__fs, sum(deltaQ), (sum(deltaQ(iAtom,:)), iAtom=1, this%nAtom)
    end if

  end subroutine writeTDOutputs


  !! Initialize matrices for populations
  subroutine tdPopulInit(this, Eiginv, EiginvAdj, HSq, iSpin)
    type(TElecDynamics), intent(in) :: this
    complex(cp), intent(out) :: Eiginv(:,:,:), EiginvAdj(:,:,:)
    real(dp), intent(in) :: Hsq(:,:,:)
    real(dp), allocatable :: T2(:,:), T3(:,:)
    integer, intent(in) :: iSpin
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


  !! Calculate populations at each time step
  subroutine getTDPopulations(this, Rho, Eiginv, EiginvAdj, T1, &
       & populDat, time, iSpin)
    type(TElecDynamics), intent(in) :: this
    complex(cp), intent(in) :: Rho(:,:,:)
    complex(cp), intent(inout) :: Eiginv(:,:,:), EiginvAdj(:,:,:)
    complex(cp), intent(inout) :: T1(:,:)
    complex(cp) :: T2(this%nOrbs, this%nOrbs)  ! Does it allocate every time?
    real(dp) :: T3(this%nOrbs, this%nOrbs)
    real(dp), intent(in) :: time
    integer, intent(in) :: populDat(2), iSpin
    integer :: iOrb

    T3 = 0.0_dp

    call gemm(T1, Eiginv(:,:,iSpin), Rho(:,:,iSpin))
    call gemm(T2, T1, EiginvAdj(:,:,iSpin))

    write(populDat(iSpin),'(*(2x,F25.15))') time * au__fs, &
         & (real(T2(iOrb, iOrb), dp), iOrb=1, this%nOrbs)

  end subroutine getTDPopulations


  !! Timing functions
  subroutine tic(t)
    integer, intent(out) :: t
    call system_clock(t)
  end subroutine tic

  subroutine tac(deltat,ti)
    integer, intent(in) :: ti
    real(dp), intent(out) :: deltat
    integer :: tf,rate
    call system_clock(tf, rate)
    deltat = real(tf-ti, dp)/real(rate, dp)
  end subroutine tac

end module timeprop_module
