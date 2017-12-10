! Main Dynamics module

module timeprop_module

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

type ElecDynamicsInp
   real(dp) :: tdField, tdDt, tdOmega, tdReFieldPolVec(3), tdImFieldPolVec(3)
   real(dp) :: tdTime0, tdTime1
   integer :: tdPolDir, tdSteps, tdSaveEvery, tdType, tdEnvType, tdSpType
   logical :: tdPopulations, tdRestart, tdWriteRestart
   real(dp) :: tdPhase
end type ElecDynamicsInp

type ElecDynamics
   private
   real(dp) :: Field, Dt, Omega, Time0, Time1
   complex(cp) :: FieldDir(3)
   real(dp), allocatable :: TDFunction(:, :), phase
   integer :: PolDir=0, Nsteps, SaveEvery, Type, EnvType, SpType
   integer :: nAtom, nOrbs, nSpin=1
   integer, allocatable :: species(:)
   character(mc), allocatable :: speciesName(:)
   logical :: Populations, SpinPol=.false.
   logical :: Restart, WriteRestart, PumpProbe
   logical :: doLaser = .false., doKick = .false., EnvFromFile = .false.
   type(TScc), allocatable :: sccCalc
   type(TEnvironment) :: env
end type ElecDynamics

private
integer :: ii,jj,kk,ll

public :: runDynamics, initElecDynamics
public :: ElecDynamicsInp, ElecDynamics
public :: iKick, iLaser, iNoTDPert
public :: iTDConstant, iTDGaussian, iTDSin2, iTDFromFile
public :: iTDSinglet, iTDTriplet

!! Enumerating available types of perturbation
integer, parameter :: iKick = 1, iLaser = 2, iNoTDPert = 3

!! Enumerating available types of envelope function
integer, parameter :: iTDConstant = 1, iTDGaussian = 2, iTDSin2 = 3, iTDFromFile = 5

!! Enumerating available types of spin polarized spectra
integer, parameter :: iTDSinglet = 1, iTDTriplet = 2


contains

!! This initializes the input variables
  subroutine initElecDynamics(this, inp, species, speciesName, env)
    type(ElecDynamics), intent(out) :: this
    type(ElecDynamicsInp), intent(in) :: inp
    integer, allocatable, intent(in) :: species(:)
    character(mc), allocatable, intent(in) :: speciesName(:)
    type(TEnvironment), intent(in) :: env
   
    real(dp) :: norm
    complex(cp) :: im
    im = cmplx(0, 1, cp)

    this%Field = inp%tdField
    this%Dt = inp%tdDt
    this%Nsteps = inp%tdSteps
    this%Type = inp%tdType
    this%EnvType = inp%tdEnvType
    this%SpType = inp%tdSpType
    this%Populations = inp%tdPopulations
    this%Restart = inp%tdRestart
    this%WriteRestart = inp%tdWriteRestart
    this%phase = inp%tdPhase
    this%SaveEvery = inp%tdSaveEvery
    allocate(this%species, source=species)
    allocate(this%speciesName, source=speciesName)
    allocate(this%sccCalc)
    this%env = env

    if (inp%tdEnvType /= iTDConstant) then
       this%Time0 = inp%tdTime0
       this%Time1 = inp%tdTime1
    end if

    if (inp%tdType == iLaser) then
       this%doLaser = .true.
    else if (inp%tdType == iKick) then
       this%doKick = .true.
    end if

    if (this%doLaser) then
       this%Omega = inp%tdOmega
       this%FieldDir = inp%tdReFieldPolVec + im * inp%tdImFieldPolVec
       norm = sqrt(sum(this%FieldDir(:)*conjg(this%FieldDir(:))))
       this%FieldDir = this%FieldDir/norm ! normalize polarization vector
       allocate(this%TDFunction(0:this%Nsteps, 3))
       if (this%EnvType == iTDFromFile) then
          this%EnvFromFile = .true.
       end if
       call getTDFunction(this)
    end if

    if (this%doKick) then
       this%PolDir = inp%tdPolDir
    end if
  end subroutine initElecDynamics


  !! Driver of calculation in order to perform a spectrum
  subroutine runDynamics(sf,Hsq,ham,H0,q0,over,filling,neighborList,nNeighbor,&
       &iSquare,iPair,img2CentCell,orb,coord,W,pRepCont, sccCalc)
    type(ElecDynamics) :: sf
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
    integer :: iPol

    sf%sccCalc = sccCalc

    sf%nSpin = size(ham(:,:), dim=2)
    if (sf%nSpin > 1) then
       sf%SpinPol = .true.
       call qm2ud(q0)
    end if

    sf%nOrbs = size(Hsq, dim=1)
    sf%nAtom = size(coord, dim=2)

    if (sf%PolDir == 4) then
       do iPol = 1,3
          sf%PolDir = iPol
          call doDynamics(sf,Hsq,ham,H0,q0,over,filling,neighborList,nNeighbor, &
               &iSquare,iPair,img2CentCell,orb,coord,W,pRepCont)
       end do
    else
       call doDynamics(sf,Hsq,ham,H0,q0,over,filling,neighborList,nNeighbor, &
            &iSquare,iPair,img2CentCell,orb,coord,W,pRepCont) 
    end if

  end subroutine runDynamics


  !! Runs the electronic dynamics of the system
  subroutine doDynamics(sf,Hsq,ham,H0,q0,over,filling,neighborList,nNeighbor, & 
       &iSquare,iPair,img2CentCell,orb,coord,W,pRepCont)
    type(ElecDynamics) :: sf
    real(dp), intent(inout) :: Hsq(:,:,:), H0(:), q0(:,:,:)
    real(dp), allocatable, intent(inout) :: ham(:,:), over(:), coord(:,:)
    real(dp), intent(in) :: W(:,:,:), filling(:,:,:)
    integer, intent(inout) :: nNeighbor(:)
    integer, allocatable, intent(inout) :: iPair(:,:), img2CentCell(:)
    integer, intent(in) :: iSquare(:)
    type(TNeighborList), intent(inout) :: neighborList
    type(ORepCont), intent(in) :: pRepCont
    type(TOrbitals), intent(in) :: orb

    complex(cp) :: Ssqr(sf%nOrbs,sf%nOrbs), Sinv(sf%nOrbs,sf%nOrbs), T1(sf%nOrbs,sf%nOrbs)
    complex(cp) :: Rho(sf%nOrbs,sf%nOrbs,sf%nSpin), Rhoold(sf%nOrbs,sf%nOrbs,sf%nSpin)
    complex(cp) :: H1(sf%nOrbs,sf%nOrbs,sf%nSpin)
    complex(cp) :: Rhonew(sf%nOrbs,sf%nOrbs,sf%nSpin)
    complex(cp), allocatable :: Eiginv(:,:,:), EiginvAdj(:,:,:)
    real(dp) :: qInput(orb%mOrb, sf%nAtom, sf%nSpin), qq(sf%nAtom,sf%nSpin), mu(3,sf%nSpin)
    real(dp) :: chargePerShell(orb%mShell,sf%nAtom,sf%nSpin)
    real(dp), allocatable :: rhoPrim(:,:), ham0(:)
    real(dp) :: time, dTime, startTime, timeInit = 0.0_dp, timeElec = 0.0_dp
    integer :: muDat, qDat, energyDat, populDat(2)
    integer :: iStep = 0, iAtom, iSpin, iTimeInit, iTimeElec, iTimeLoop
    type(TPotentials) :: potential
    type(TEnergies) :: energy
    character(4) :: dumpIdx

    
    ! Initialize timer
    call tic(iTimeInit) 
    startTime = 0.0_dp   
    if (sf%Restart) then
       call readRestart(Rho, Ssqr, coord, startTime)
    end if

    ! Initialize stuff  
    call createMatrices(sf, Rho, H1, Ssqr, Sinv, H0, ham0, &
         & over, ham, Hsq, filling, orb, rhoPrim, potential, &
         & neighborList%iNeighbor, nNeighbor, iSquare, iPair, img2CentCell, & 
         & Eiginv, EiginvAdj)
    ! Initialize output files
    call initTDOutput(sf, muDat, qDat, energyDat, populDat)

    call getChargeMu(sf, qq, qInput, mu, q0, Rho, Ssqr, coord, iSquare) !q_0
    call addSCCandField(sf, H1, ham, over, ham0, qInput, q0, coord, orb, potential, &  !H1_0
         &neighborList%iNeighbor, nNeighbor, iSquare, iPair, img2CentCell, iStep, chargePerShell, W)

    ! Apply kick to Rho if necessary
    if (sf%doKick) then
       call kickDM(sf, Rho, Ssqr, Sinv, iSquare, coord)
    end if

    ! Initialize electron dinamics
    Rhoold(:,:,:) = Rho
    call euler(sf, sf%Dt, Rho, Rhoold, H1, Sinv)

    call getTDEnergy(sf, energy, rhoPrim, Rhoold, neighborList%iNeighbor, &
         & nNeighbor, orb, iSquare, iPair, img2CentCell, ham0, qInput, q0, &
         & potential, chargePerShell, coord, pRepCont)

    call tac(timeInit, iTimeInit)
    ! End of initialization

    !! Main loop
    call tic(iTimeLoop)
    write(*,*)'Starting dynamics'

    do iStep = 0, sf%Nsteps
       time = iStep * sf%Dt + startTime
       ! Write everything to file
       if ((sf%Restart) .and. (iStep > 0) .or. (.not. sf%Restart)) then
          call writeTDOutputs(sf, muDat, qDat, energyDat, &
               & time, energy, mu, qq, iStep)
       end if

       call getChargeMu(sf, qq, qInput, mu, q0, Rho, Ssqr, coord, iSquare) !q_1
       call addSCCandField(sf, H1, ham, over, ham0, qInput, q0, coord, orb, &
           & potential, neighborList%iNeighbor, nNeighbor, iSquare, iPair, &
           & img2CentCell, iStep, chargePerShell, W) ! H_1

       if ((sf%WriteRestart) .and. (iStep > 0) .and. &
            & (mod(iStep,sf%Nsteps/10) == 0)) then
          call writeRestart(Rho, Ssqr, coord, time)    
       end if

       call getTDEnergy(sf, energy, rhoPrim, Rho, neighborList%iNeighbor, &
            & nNeighbor, orb, iSquare, iPair, img2CentCell, ham0, qInput, q0, &
            & potential, chargePerShell, coord, pRepCont) ! E_1

       do iSpin = 1, sf%nSpin
          call tic(iTimeElec)          
          call zscal(sf%nOrbs*sf%nOrbs, cmplx(0, 1, cp), H1(:,:,iSpin), 1)                    
          call propagateRhoCPU(sf, Rhoold(:,:,iSpin), Rho(:,:,iSpin), &
               & H1(:,:,iSpin), Sinv, T1, 2.0_dp * sf%Dt)
          call zswap(sf%nOrbs*sf%nOrbs, Rhoold(:,:,iSpin), 1, Rho(:,:,iSpin), 1)
          call tac(dTime,iTimeElec)
          timeElec = timeElec + dTime

          if ((sf%Populations) .and. (mod(iStep, sf%SaveEvery) == 0)) then
             call getTDPopulations(sf, Rho, Eiginv, EiginvAdj, T1, &
                  & populDat, time, iSpin)
          end if
       end do

       if (mod(iStep,sf%Nsteps/10) == 0) then
          call tac(dTime,iTimeLoop)
          write(*,*)'Step ',iStep,'elapsed loop time:',dTime,'average time per loop',dTime/(iStep+1)
       end if
    end do
    
    write(*,*)'Dynamics OK!'
    call tac(dTime, iTimeInit)
    write(*,*)'Total time',dTime
    write(*,*)'Time spent in initialization',timeInit
    write(*,*)'Time spent in electronic propagation',timeElec

    close(energyDat)
    close(qDat)
    close(muDat)
  end subroutine doDynamics


  
  !! Applies SCC to hamiltonian and adds TD perturbation (if any)
  subroutine addSCCandField(sf, H1, ham, over, ham0, qInput, q0, coord, orb, potential, &
       &iNeighbor, nNeighbor, iSquare, iPair, img2CentCell, iStep, chargePerShell, W)
    type(ElecDynamics), intent(inout) :: sf
    complex(cp), intent(inout) :: H1(:,:,:)
    real(dp), allocatable, intent(inout) :: ham(:,:), over(:), coord(:,:)
    real(dp), intent(in) :: ham0(:)
    real(dp), intent(inout) :: qInput(:,:,:), q0(:,:,:)
    real(dp), intent(in) :: W(:,:,:)
    type(TOrbitals), intent(in) :: orb
    integer, intent(in) :: iNeighbor(0:,:),nNeighbor(:)
    integer, intent(in) :: iSquare(:), iPair(0:,:), img2CentCell(:)
    integer, intent(in) :: iStep
    real(dp), intent(inout) :: chargePerShell(:,:,:)
    real(dp) :: T2(sf%nOrbs,sf%nOrbs)
    real(dp) :: shellPot(orb%mShell, sf%nAtom, sf%nSpin), atomPot(sf%nAtom, sf%nSpin)
    type(TPotentials), intent(inout) :: potential
    integer :: iAtom, iSpin, iSp1, iSh1

    ham(:,1) = ham0(:)
    
    if (sf%nSpin == 2) then
       ham(:,2) = 0.0_dp
 
       call ud2qm(qInput)
       call ud2qm(q0)
 
       chargePerShell(:,:,:) = 0.0_dp ! hack for the moment to get charge and magnetization
       do iAtom = 1, sf%nAtom
          iSp1 = sf%species(iAtom)
          do iSh1 = 1, orb%nShell(iSp1)
             chargePerShell(iSh1,iAtom,1:sf%nSpin) = &
                  & chargePerShell(iSh1,iAtom,1:sf%nSpin) + &
                  & sum(qInput(orb%posShell(iSh1,iSp1): &
                  & orb%posShell(iSh1+1,iSp1)-1,iAtom,1:sf%nSpin),dim=1)
          end do
       end do
    end if

    potential%intAtom = 0.0_dp
    potential%intShell = 0.0_dp
    potential%intBlock = 0.0_dp
    potential%extAtom = 0.0_dp
    potential%extShell = 0.0_dp
    potential%extBlock = 0.0_dp

    call sf%sccCalc%updateCharges(sf%env, qInput, q0, orb, sf%species, iNeighbor, img2CentCell)
    call sf%sccCalc%getShiftPerAtom(atomPot(:,1))
    call sf%sccCalc%getShiftPerL(shellPot(:,:,1))
    potential%intAtom(:,1) = potential%intAtom(:,1) + atomPot(:,1)
    potential%intShell(:,:,1) = potential%intShell(:,:,1) + shellPot(:,:,1)

    !! Build spin contribution (if necessary)
    if (sf%SpinPol) then
       call getSpinShift(shellPot, chargePerShell, sf%species, orb, W)
       potential%intShell = potential%intShell + shellPot
    end if
    
    call total_shift(potential%intShell, potential%intAtom, orb, sf%species)
    call total_shift(potential%intBlock, potential%intShell, orb, sf%species)

    !! Add time dependent field if necessary
    if (sf%doLaser) then
       do iAtom = 1, sf%nAtom
          potential%extAtom(iAtom, 1) = dot_product(coord(:,iAtom), &
               & sf%TDFunction(iStep, :))
       end do
       call total_shift(potential%extShell, potential%extAtom, orb, sf%species)
       call total_shift(potential%extBlock, potential%extShell, orb, sf%species)
       potential%intShell = potential%intShell + potential%extShell ! for SCC
       potential%intBlock = potential%intBlock + potential%extBlock ! for forces
    end if

    call add_shift(ham, over, nNeighbor, iNeighbor, &
         & sf%species, orb, iPair, sf%nAtom, img2CentCell, potential%intShell) !or intBlock?

    if (sf%nSpin == 2) then
       ham(:,:) = 2.0_dp * ham
       call qm2ud(ham)
       call qm2ud(q0)
    end if

    do iSpin=1,sf%nSpin
       call unpackHS(T2,ham(:,iSpin),iNeighbor,nNeighbor,iSquare,iPair,img2CentCell)
       call blockSymmetrizeHS(T2,iSquare)
       H1(:,:,iSpin) = cmplx(T2, 0.0_dp, dp)
    end do

  end subroutine addSCCandField


  !! Kick the density matrix for spectrum calculations
  subroutine kickDM(sf, Rho, Ssqr, Sinv, iSquare, coord)
    type(ElecDynamics), intent(in) :: sf
    complex(cp), intent(in) :: Ssqr(:,:), Sinv(:,:)
    complex(cp), intent(inout) :: Rho(:,:,:)
    real(dp), allocatable, intent(in) :: coord(:,:)
    integer, intent(in) :: iSquare(:)
    complex(cp) :: T1(sf%nOrbs,sf%nOrbs,sf%nSpin),T3(sf%nOrbs,sf%nOrbs,sf%nSpin)
    complex(cp) :: T2(sf%nOrbs,sf%nOrbs),T4(sf%nOrbs,sf%nOrbs)
    integer :: iAt, iStart, iEnd, iSpin
    real(dp) :: pkick(sf%nSpin)

    pkick(1) = sf%Field ! check units

    if (sf%nSpin == 2 .and. sf%SpType == iTDSinglet) pkick(2) = pkick(1)
    if (sf%nSpin == 2 .and. sf%SpType == iTDTriplet) pkick(2) = - pkick(1)

    T1 = 0.0_dp
    T2 = 0.0_dp
    T3 = 0.0_dp
    T4 = 0.0_dp

    do iSpin=1,sf%nSpin
       do iAt = 1,sf%nAtom
          iStart = iSquare(iAt)
          iEnd = iSquare(iAt+1)-1
          do jj = iStart, iEnd
             T1(jj,jj,iSpin) = Exp(cmplx(0, -pkick(iSpin)*coord(sf%PolDir, iAt), cp))
             T3(jj,jj,iSpin) = Exp(cmplx(0, pkick(iSpin)*coord(sf%PolDir, iAt), cp))
          end do
       end do
    end do

    do iSpin=1,sf%nSpin
       call gemm(T2, T1(:,:,iSpin), Rho(:,:,iSpin), cmplx(1, 0, cp), cmplx(0, 0, cp))
       call gemm(T4, T2, Ssqr, cmplx(1, 0, cp), cmplx(0, 0, cp))
       call gemm(T2, T4, T3(:,:,iSpin), cmplx(1, 0, cp), cmplx(0, 0, cp))
       !! To account for non-commutating Rho and H, I replace the following line
       !call gemm(Rho(:,:,iSpin), T2, Sinv, cmplx(1, 0, cp), cmplx(0, 0, cp))
       !! By:
       call gemm(Rho(:,:,iSpin), T2, Sinv, cmplx(0.5, 0, cp), cmplx(0, 0, cp), 'N', 'N')
       call gemm(Rho(:,:,iSpin), Sinv, T2, cmplx(0.5, 0, cp), cmplx(1, 0, cp), 'N', 'C')       
    end do

    write(*,*)'Density kicked!'

  end subroutine kickDM


  !! Calculate envelope function
  subroutine getTDFunction(sf)
    type(ElecDynamics), intent(inout) :: sf
    real(dp) :: midPulse, deltaT, angFreq, E0, time, env, TDFuncMod
    real(dp) :: tdfun(3)
    integer :: iStep, laserDat
    complex(cp) :: im

    im = cmplx(0, 1, cp)
    midPulse = (sf%Time0 + sf%Time1)/2.0_dp
    deltaT = sf%Time1 - sf%Time0
    angFreq = sf%Omega
    E0 = sf%Field
    sf%TDFunction(:,:) = 0.0_dp

    open(newunit=laserDat, file='laser.dat')

    do iStep = 0,sf%Nsteps
       time = iStep * sf%Dt

       if (sf%EnvType /= iTDConstant .and. (time >= sf%Time0) .and. (time <= sf%Time1)) then
          if (sf%EnvType == iTDGaussian) then
             env = Exp(-(time-midPulse)**2/(2.*248.**2))
          else if (sf%EnvType == iTDSin2) then
             env = Sin(pi*(time+sf%Time0)/deltaT)**2
          end if
       else if (sf%EnvType == iTDConstant) then
          env = 1.0_dp
       else
          env = 0.0_dp
       end if

       if (sf%EnvType /= iTDFromFile) then
          sf%TDFunction(iStep, :) = E0 * env * aimag(Exp(im*(time*angFreq + sf%phase)) * sf%FieldDir)
          TDFuncMod = sqrt(sum(sf%TDFunction(iStep, :)**2))
       end if

       if (sf%EnvFromFile) then
          read(laserDat, *)time, TDFuncMod, tdfun(1), tdfun(2), tdfun(3)
          sf%TDFunction(iStep, :) = sf%TDFunction(iStep, :) + tdfun(:) * (Bohr__AA / Hartree__eV)
       else
          write(laserDat, "(5F15.8)") time * au__fs, &
               & TDFuncMod * (Hartree__eV / Bohr__AA), &
               & sf%TDFunction(iStep, :) * (Hartree__eV / Bohr__AA)
       end if
       
    end do

    close(laserDat)    
  end subroutine getTDFunction


  !! Calculate charges, dipole moments
    subroutine getChargeMu(sf, qq, qInput, mu, q0, Rho, Ssqr, coord, iSquare)
    type(ElecDynamics), intent(in) :: sf
    real(dp), intent(out) :: qInput(:,:,:), mu(:,:), qq(:,:)
    real(dp), intent(in) :: coord(:,:), q0(:,:,:)
    complex(cp), intent(in) :: Rho(:,:,:), Ssqr(:,:)
    integer, intent(in) :: iSquare(:)
    integer :: iAt, iOrb, iSpin

    qInput = 0.0_dp
    mu = 0.0_dp

    do iSpin=1,sf%nSpin
       do iAt = 1,sf%nAtom
          iOrb = 0
          do jj = iSquare(iAt),iSquare(iAt+1)-1
             iOrb = iOrb + 1
             qInput(iOrb,iAt,iSpin) = sum(Rho(jj,:,iSpin)*Ssqr(jj,:)) ! 1 <= iOrb <= orb%mOrb
          end do
          mu(1,iSpin) = mu(1,iSpin) + sum(q0(:,iAt,iSpin)-qInput(:,iAt,iSpin))&
               & * coord(1,iAt)
          mu(2,iSpin) = mu(2,iSpin) + sum(q0(:,iAt,iSpin)-qInput(:,iAt,iSpin))&
               & * coord(2,iAt)
          mu(3,iSpin) = mu(3,iSpin) + sum(q0(:,iAt,iSpin)-qInput(:,iAt,iSpin))&
               & * coord(3,iAt)
       end do
    end do

    qq(:,:) = sum((qInput(:,:,:)-q0(:,:,:)), dim=1)
  end subroutine getChargeMu


  !! Calculate energy
  subroutine getTDEnergy(sf, energy, rhoPrim, Rho, iNeighbor, &
       & nNeighbor, orb, iSquare, iPair, img2CentCell, ham0, qInput, q0, &
       & potential, chargePerShell, coord, pRepCont)
    type(ElecDynamics), intent(inout) :: sf
    type(TEnergies), intent(out) :: energy
    real(dp), allocatable, intent(inout) :: rhoPrim(:,:)
    complex(cp), intent(in) :: Rho(:,:,:)
    real(dp), intent(in) :: ham0(:)
    real(dp), intent(inout) :: qInput(:,:,:)
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

    call init(energy, sf%nAtom)
    energy%ETotal = 0.0_dp

    energy%EnonSCC = 0.0_dp
    energy%atomNonSCC(:) = 0.0_dp
    call mulliken(energy%atomNonSCC(:), rhoPrim(:,1), ham0, orb,&
         &iNeighbor, nNeighbor, img2CentCell, iPair)
    energy%EnonSCC =  sum(energy%atomNonSCC)

    if (sf%doLaser) then ! energy in external field
       energy%atomExt = -sum(q0(:, :, 1) - qInput(:, :, 1),dim=1) &
            & * potential%extAtom(:,1)
       energy%Eext =  sum(energy%atomExt)
    else
       energy%Eext = 0.0_dp
       energy%atomExt = 0.0_dp
    end if

    call sf%sccCalc%getEnergyPerAtom(energy%atomSCC)
    energy%eSCC = sum(energy%atomSCC)

    if (sf%nSpin > 1) then
       energy%atomSpin = 0.5_dp * sum(sum(potential%intShell(:,:,2:sf%nSpin)&
            & * chargePerShell(:,:,2:sf%nSpin), dim=1),dim=2)
       energy%Espin = sum(energy%atomSpin)
    else
       energy%atomSpin = 0.0_dp
       energy%eSpin = 0.0_dp
    end if

    !! Calculate repulsive energy
    call getERep(energy%atomRep, coord, nNeighbor, iNeighbor, &
         &sf%species, pRepCont, img2CentCell)
    energy%Erep = sum(energy%atomRep)
    
    energy%Eelec = energy%EnonSCC + energy%eSCC + energy%Espin + energy%Eext
    energy%Etotal = energy%Eelec + energy%Erep + energy%eDisp
    
  end subroutine getTDEnergy


  !! Create all necessary matrices for dynamics
  subroutine createMatrices(sf, Rho, H1, Ssqr, Sinv, H0, ham0, &
       & over, ham, Hsq, filling, orb, rhoPrim, potential, &
       & iNeighbor, nNeighbor, iSquare, iPair, img2CentCell, Eiginv, &
       & EiginvAdj)
    type(ElecDynamics), intent(inout) :: sf
    real(dp), intent(inout) :: Hsq(:,:,:)
    real(dp), allocatable, intent(in) :: over(:), ham(:,:)
    real(dp), intent(in) :: filling(:,:,:), H0(:)
    integer, intent(in) :: iNeighbor(0:,:),nNeighbor(:)
    integer, intent(in) :: iSquare(:), iPair(0:,:), img2CentCell(:)
    type(TOrbitals), intent(in) :: orb
    type(TPotentials), intent(out) :: potential

    real(dp), allocatable, intent(out) :: rhoPrim(:,:), ham0(:)
    complex(cp), intent(out) :: Ssqr(:,:), Sinv(:,:), H1(:,:,:)
    complex(cp), intent(out) :: Rho(:,:,:)
    complex(cp), allocatable, intent(out), optional :: Eiginv(:,:,:), EiginvAdj(:,:,:)
    real(dp) :: T2(sf%nOrbs,sf%nOrbs), T3(sf%nOrbs, sf%nOrbs)
    integer :: iSpin

    allocate(rhoPrim(size(ham), sf%nSpin))
    allocate(ham0(size(H0)))
    ham0(:) = H0
    
    T2 = 0.0_dp
    T3 = 0.0_dp
    call unpackHS(T2,over,iNeighbor,nNeighbor,iSquare,iPair,img2CentCell)
    call blockSymmetrizeHS(T2,iSquare)
    Ssqr(:,:) = cmplx(T2, 0, cp) ! Overlap

    do iSpin=1,sf%nSpin
       call unpackHS(T3,ham(:,iSpin),iNeighbor,nNeighbor,iSquare,iPair,img2CentCell)
       call blockSymmetrizeHS(T3,iSquare)
       H1(:,:,iSpin) = cmplx(T3, 0, cp) ! Hamiltonian
       T3 = 0.0_dp
    end do

    if (sf%Populations) then
       allocate(Eiginv(sf%nOrbs, sf%nOrbs, sf%nSpin))
       allocate(EiginvAdj(sf%nOrbs, sf%nOrbs, sf%nSpin))
       do iSpin=1,sf%nSpin
          call TDPopulInit(sf, Eiginv, EiginvAdj, HSq, iSpin)
       end do
    end if

    do ii = 1, sf%nOrbs
       T3(ii,ii) = 1.0_dp
    end do
    call gesv(T2,T3)
    Sinv(:,:) = cmplx(T3, 0, cp) ! Inverse overlap
    write(*,*)'S inverted'

    do iSpin=1,sf%nSpin
       T2 = 0.0_dp
       call makeDensityMatrix(T2,Hsq(:,:,iSpin),filling(:,1,iSpin))
       Rho(:,:,iSpin) = cmplx(T2, 0, cp) ! Density matrix
       do kk = 1, sf%nOrbs-1
          do ll = kk+1, sf%nOrbs
             Rho(kk,ll,iSpin) = Rho(ll,kk,iSpin) ! Symmetrization neccesary
          end do
       end do
    end do

    call init(potential, orb, sf%nAtom, sf%nSpin)
  end subroutine createMatrices

  
  !! Perfoms a step backwards to boot the dynamics using the Euler algorithm
  subroutine euler(sf, step, Rhoold, Rho, H1, Sinv)
    type(ElecDynamics), intent(in) :: sf
    complex(cp), intent(in) :: Rho(:,:,:), Sinv(:,:)
    complex(cp), intent(inout) :: H1(:,:,:) ! Warning! H1 is modified
    complex(cp), intent(out) :: Rhoold(:,:,:)
    real(dp), intent(in) :: step ! Step with its sign, to perform euler backwards (-) or forwards (+)
    complex(cp) :: T1(sf%nOrbs,sf%nOrbs)
    integer :: iSpin

    Rhoold(:,:,:) = Rho

    do iSpin=1,sf%nSpin
       T1 = 0.0_cp
       H1(:,:,iSpin) = cmplx(0, 1, cp) * H1(:,:,iSpin)
       call propagateRhoCPU(sf, Rhoold(:,:,iSpin), Rho(:,:,iSpin), &
            & H1(:,:,iSpin), Sinv, T1, step)
    end do

  end subroutine euler


  !! Propagate Rho, notice that H = iH (coeficients are real)
  subroutine propagateRhoCPU(sf, Rhoold, Rho, H1, Sinv, T1, step)
    type(ElecDynamics), intent(in) :: sf
    complex(cp), intent(inout) :: Rhoold(:,:)
    complex(cp), intent(in) :: Rho(:,:), H1(:,:), Sinv(:,:)
    complex(cp), intent(out) :: T1(:,:)
    real(dp), intent(in) :: step

    call gemm(T1, Sinv, H1, cmplx(1, 0, cp), cmplx(0, 0, cp), &
         & 'N', 'N', sf%nOrbs, sf%nOrbs, sf%nOrbs)
    call gemm(Rhoold, T1, Rho, cmplx(-step, 0, cp),&
         & cmplx(1, 0, cp), 'N', 'N', sf%nOrbs, sf%nOrbs, sf%nOrbs)
    call gemm(Rhoold, Rho, T1, cmplx(-step, 0, cp),&
         & cmplx(1, 0, cp), 'N', 'C', sf%nOrbs, sf%nOrbs, sf%nOrbs)
  end subroutine propagateRhoCPU


  !! Initialize output files
  subroutine initTDOutput(sf, muDat, qDat, energyDat, populDat)
    type(ElecDynamics), intent(in) :: sf
    integer, intent(out) :: muDat, qDat, energyDat, populDat(2)
    character(20) :: muFileName
    character(1) :: strSpin
    integer :: iSpin
    logical :: exist

    if (sf%doKick) then
       if (sf%PolDir == 1) then
          muFileName = 'mux.dat'
       else if (sf%PolDir == 2) then
          muFileName = 'muy.dat'
       else if (sf%PolDir == 3) then
          muFileName = 'muz.dat'
       end if
    else
       muFileName = 'mu.dat'
    end if

    call openFile(sf, muDat, muFileName)
    call openFile(sf, qDat, 'qsvst.dat')
    call openFile(sf, energyDat, 'energyvst.dat')

    if (sf%Populations) then
       do iSpin=1,sf%nSpin
          write(strSpin,'(i1)')iSpin
          call openFile(sf, populDat(iSpin), 'molPopul' // trim(strSpin) // '.dat')
       end do
    end if

  end subroutine initTDOutput

  
  !! Open files in different ways deppending on their previous existance
  subroutine openFile(sf, unitName, fileName)
    type(ElecDynamics), intent(in) :: sf
    logical :: exist=.false.
    integer :: unitName
    character(*) :: fileName
    if (sf%Restart) then
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
  subroutine writeTDOutputs(sf, muDat, qDat, energyDat, &
       & time, energy, mu, qq, iStep)
    type(ElecDynamics), intent(in) :: sf
    type(TEnergies), intent(in) :: energy
    integer, intent(in) :: muDat, qDat, energyDat
    real(dp), intent(in) :: time, mu(:,:), qq(:,:)
    integer, intent(in) :: iStep
    integer :: iAtom, iSpin

    write(energydat, '(9F25.15)') time * au__fs, energy%Etotal, energy%EnonSCC, &
         & energy%eSCC, energy%Espin, energy%Eext, energy%Erep
    write(muDat, '(7F25.15)') time * au__fs, ((mu(ii, iSpin) * Bohr__AA, ii=1, 3), iSpin=1, sf%nSpin)
    if (mod(iStep, sf%SaveEvery) == 0) then
       write(qDat, '(5000F25.15)') time * au__fs, sum(qq(:,:)), (sum(qq(iAtom,:)), iAtom=1,sf%nAtom)
    end if

  end subroutine writeTDOutputs


  !! Initialize matrices for populations
  subroutine TDPopulInit(sf, Eiginv, EiginvAdj, HSq, iSpin)
    type(ElecDynamics), intent(in) :: sf 
    complex(cp), intent(out) :: Eiginv(:,:,:), EiginvAdj(:,:,:)
    real(dp), intent(in) :: Hsq(:,:,:)
    real(dp), allocatable :: T2(:,:), T3(:,:)
    integer, intent(in) :: iSpin
    allocate(T2(sf%nOrbs, sf%nOrbs), T3(sf%nOrbs, sf%nOrbs))

    T2 = Hsq(:,:,iSpin)
    T3 = 0.0_dp
    do ii = 1, sf%nOrbs
       T3(ii,ii) = 1.0_dp
    end do
    call gesv(T2,T3)
    Eiginv(:, :, iSpin) = cmplx(T3, 0, cp)

    T3 = 0.0_dp
    do ii = 1, sf%nOrbs
       T3(ii, ii) = 1.0_dp
       do jj = 1, sf%nOrbs
          T2(ii, jj) = Hsq(jj, ii, iSpin)
       end do
    end do
    call gesv(T2,T3)
    EiginvAdj(:, :, iSpin) = cmplx(T3, 0, cp)

    deallocate(T2, T3)

  end subroutine TDPopulInit


  !! Calculate populations at each time step
  subroutine getTDPopulations(sf, Rho, Eiginv, EiginvAdj, T1, &
       & populDat, time, iSpin)
    type(ElecDynamics), intent(in) :: sf
    complex(cp), intent(in) :: Rho(:,:,:)
    complex(cp), intent(inout) :: Eiginv(:,:,:), EiginvAdj(:,:,:)
    complex(cp), intent(inout) :: T1(:,:)
    complex(cp) :: T2(sf%nOrbs, sf%nOrbs)  ! Does it allocate every time?
    real(dp) :: T3(sf%nOrbs, sf%nOrbs)
    real(dp), intent(in) :: time
    integer, intent(in) :: populDat(2), iSpin

    T3 = 0.0_dp

    call gemm(T1, Eiginv(:,:,iSpin), Rho(:,:,iSpin), cmplx(1, 0, cp), & 
         & cmplx(0, 0, cp), 'N', 'N', sf%nOrbs, sf%nOrbs, sf%nOrbs)
    call gemm(T2, T1, EiginvAdj(:,:,iSpin), cmplx(1, 0, cp), &
         & cmplx(0, 0, cp), 'N', 'N', sf%nOrbs, sf%nOrbs, sf%nOrbs)

    write(populDat(iSpin),'(5000F25.15)') time * au__fs, &
         & (real(T2(ii,ii), dp), ii=1, sf%nOrbs)

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
