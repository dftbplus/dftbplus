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
use forces
use repulsive
use slakocont
use repcont
use thermostat
use mdintegrator
use dummytherm
use mdcommon
use ranlux
use periodic 
use velocityverlet
use nonscc
use energies, only: TEnergies
use populations
use repulsive
use eigenvects
use sk
use dispiface

implicit none

type ElecDynamicsInp
   real(dp) :: tdField, tdDt, tdOmega, tdReFieldPolVec(3), tdImFieldPolVec(3)
   real(dp) :: tdTime0, tdTime1
   integer :: tdPolDir, tdSteps, tdSaveEvery, tdType, tdEnvType, tdSpType
   integer, allocatable :: indMovedAtom(:)
   integer :: nMovedAtom, tdEulerEvery, tdPPFrames
   logical :: tdPopulations, tdForces, tdIons, tReadMDVelocities, tddoEulers
   logical :: tdPairWise, tdRestart, tdWriteRestart, tdPumpProbe, tdOnsiteGradients
   real(dp) :: tempAtom, tdPhase
   real(dp), allocatable :: initialVelocities(:,:)
end type ElecDynamicsInp

type ElecDynamics
   private
   real(dp) :: Field, Dt, Omega, Time0, Time1
   complex(cp) :: FieldDir(3)
   real(dp), allocatable :: TDFunction(:, :)
   real(dp), allocatable :: initialVelocities(:,:)
   real(dp), allocatable :: movedVelo(:,:), movedMass(:,:), onsiteGrads(:,:,:,:)
   real(dp) :: mCutoff, skRepCutoff, phase
   real(dp) :: rCellVec(3,1)
   real(dp), allocatable :: atomEigVal(:,:)
   integer :: PolDir=0, Nsteps, SaveEvery, Type, EnvType, SpType
   integer :: nAtom, nOrbs, nSpin=1, nMovedAtom, nSparse, EulerEvery, PuProbeFrames
   integer, allocatable :: iCellVec(:)
   integer, allocatable :: species(:)
   integer, allocatable :: indMovedAtom(:)
   character(mc), allocatable :: speciesName(:)
   logical :: Populations, Ions, Forces, SpinPol=.false., tDispersion=.false.
   logical :: ReadMDVelocities, Restart, WriteRestart, PumpProbe
   logical :: FirstIonStep = .true.
   logical :: doEulers = .false.
   logical :: PairWiseEnergy = .false., CalcOnsiteGradients = .false.
   logical :: doLaser = .false., doKick = .false., EnvFromFile = .false.
   type(OThermostat), allocatable :: pThermostat
   type(OMDIntegrator), allocatable :: pMDIntegrator
   class(DispersionIface), allocatable :: dispersion
   type(NonSccDiff), allocatable :: derivator
   type(TScc), allocatable :: sccCalc
end type ElecDynamics

private
integer :: ii,jj,kk,ll

public :: runDynamics, initElecDynamics
public :: ElecDynamicsInp, ElecDynamics
public :: iKick, iLaser, iKickAndLaser, iNoTDPert, i2Lasers
public :: iTDConstant, iTDGaussian, iTDSin2, iTDCos, iTDFromFile
public :: iTDSinglet, iTDTriplet

!! Enumerating available types of perturbation
integer, parameter :: iKick = 1, iLaser = 2, iKickAndLaser = 3, iNoTDPert = 4, i2Lasers = 5

!! Enumerating available types of envelope function
integer, parameter :: iTDConstant = 1, iTDGaussian = 2, iTDSin2 = 3
integer, parameter :: iTDCos = 4, iTDFromFile = 5

!! Enumerating available types of spin polarized spectra
integer, parameter :: iTDSinglet = 1, iTDTriplet = 2


contains

!! This initializes the input variables
  subroutine initElecDynamics(this, inp, randomThermostat, & 
       &mass, nAtom, species, skRepCutoff, mCutoff, &
       &iCellVec, atomEigVal, speciesName, dispersion, nonSccDeriv, sccCalc)
    type(ElecDynamics), intent(out) :: this
    type(ElecDynamicsInp), intent(in) :: inp 
    real(dp), intent(in), allocatable :: atomEigVal(:,:)
    integer, intent(in) :: nAtom
    integer, allocatable, intent(in) :: species(:)
    character(mc), allocatable, intent(in) :: speciesName(:)
    integer, allocatable, intent(in) :: iCellVec(:)
    type(ORanlux), allocatable, intent(inout) :: randomThermostat  
    real(dp), intent(in) :: mCutoff, skRepCutoff, mass(:) 
    class(DispersionIface), allocatable, intent(inout) :: dispersion
    type(NonSccDiff), intent(in) :: nonSccDeriv
    type(TScc), intent(in) :: sccCalc
   
    type(ODummyThermostat), allocatable :: pDummyTherm
    type(OMDCommon), allocatable :: pMDFrame
    real(dp) :: norm, tempAtom
    logical :: tMDstill, tDispersion
    complex(cp) :: im
    im = cmplx(0, 1, cp)

    this%Field = inp%tdField
    this%Dt = inp%tdDt
    this%Nsteps = inp%tdSteps
    this%Type = inp%tdType
    this%EnvType = inp%tdEnvType
    this%SpType = inp%tdSpType
    this%Populations = inp%tdPopulations
    this%Ions = inp%tdIons
    this%Forces = inp%tdForces
    this%SaveEvery = inp%tdSaveEvery
    this%doEulers = inp%tddoEulers
    this%EulerEvery = inp%tdEulerEvery
    this%PairWiseEnergy = inp%tdPairWise
    this%Restart = inp%tdRestart
    this%WriteRestart = inp%tdWriteRestart
    this%PumpProbe = inp%tdPumpProbe
    this%PuProbeFrames = inp%tdPPFrames
    this%phase = inp%tdPhase
    this%CalcOnsiteGradients = inp%tdOnsiteGradients
    allocate(this%species, source=species)
    allocate(this%speciesName, source=speciesName)
    allocate(this%sccCalc)

    if (inp%tdEnvType /= iTDConstant) then
       this%Time0 = inp%tdTime0
       this%Time1 = inp%tdTime1
    end if

    if (inp%tdType == iLaser) then
       this%doLaser = .true.
    else if (inp%tdType == iKick) then
       this%doKick = .true.
    else if (inp%tdType == iKickAndLaser) then
       this%doLaser = .true.
       this%doKick = .true.
    else if (inp%tdType == i2Lasers) then
       this%doLaser = .true.
       this%EnvFromFile = .true.
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

    if (this%Forces .and. .not. this%Ions) then
       this%nMovedAtom = nAtom
       allocate(this%movedMass(3, this%nMovedAtom))
       this%movedMass(:,:) = spread(mass(this%species(:)),1,3)
    end if

    !! For forces
    if (this%Ions) then
       this%Forces = .true.
       allocate(this%indMovedAtom, source=inp%indMovedAtom)  
       this%nMovedAtom = inp%nMovedAtom
       tempAtom = inp%tempAtom
       tMDstill = .false.
       allocate(this%initialVelocities(3, this%nMovedAtom))
       
       if (allocated(inp%initialVelocities)) then
          this%initialVelocities(:,:) = inp%initialVelocities
       end if
         
       allocate(this%movedVelo(3, this%nMovedAtom))
       allocate(this%movedMass(3, this%nMovedAtom))
       this%ReadMDVelocities = inp%tReadMDVelocities
       this%movedMass(:,:) = spread(mass(this%indMovedAtom),1,3)
       allocate(this%pThermostat)
       allocate(pMDFrame)
       call init(pMDFrame, this%nMovedAtom, nAtom, tMDstill)
       allocate(pDummyTherm)
       call init(pDummyTherm, tempAtom, mass(this%indMovedAtom), randomThermostat, &
            &pMDFrame)
       call init(this%pThermostat, pDummyTherm)
       allocate(this%derivator, source = nonSccDeriv)
    else
       allocate(this%movedVelo(3, nAtom))
       this%movedVelo = 0.0_dp
     end if

     tDispersion = allocated(dispersion)
     if (tDispersion) then
        this%tDispersion = .true.
        allocate(this%dispersion, source=dispersion)
     end if

    this%skRepCutoff = skRepCutoff
    this%mCutoff = mCutoff
    allocate(this%iCellVec, source=iCellVec)
    this%rCellVec(:,1) = [0.0_dp, 0.0_dp, 0.0_dp]
    allocate(this%atomEigVal, source=atomEigVal)      
  end subroutine initElecDynamics


  !! Driver of calculation in order to perform a spectrum
  subroutine runDynamics(sf,Hsq,ham,H0,q0,over,filling,neighborList,nNeighbor,&
       &iSquare,iPair,img2CentCell,orb,coord,W,skHamCont,skOverCont,pRepCont, sccCalc)
    type(ElecDynamics) :: sf
    real(dp), intent(inout) :: Hsq(:,:,:), H0(:), q0(:,:,:)
    real(dp), allocatable, intent(inout) :: ham(:,:), over(:), coord(:,:)
    real(dp), intent(in) :: W(:,:,:), filling(:,:,:)
    integer, intent(inout) :: nNeighbor(:)
    integer, allocatable, intent(inout) :: iPair(:,:), img2CentCell(:)
    integer, intent(in) :: iSquare(:)
    type(TNeighborList), intent(inout) :: neighborList
    type(OSlakoCont), intent(in) :: skHamCont, skOverCont
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
               &iSquare,iPair,img2CentCell,orb,coord,W,skHamCont,skOverCont,pRepCont)
       end do
    else
       call doDynamics(sf,Hsq,ham,H0,q0,over,filling,neighborList,nNeighbor, &
            &iSquare,iPair,img2CentCell,orb,coord,W,skHamCont,skOverCont,pRepCont) 
    end if

  end subroutine runDynamics


  !! Runs the electronic dynamics of the system
  subroutine doDynamics(sf,Hsq,ham,H0,q0,over,filling,neighborList,nNeighbor, & 
       &iSquare,iPair,img2CentCell,orb,coord,W,skHamCont,skOverCont,pRepCont)
    type(ElecDynamics) :: sf
    real(dp), intent(inout) :: Hsq(:,:,:), H0(:), q0(:,:,:)
    real(dp), allocatable, intent(inout) :: ham(:,:), over(:), coord(:,:)
    real(dp), intent(in) :: W(:,:,:), filling(:,:,:)
    integer, intent(inout) :: nNeighbor(:)
    integer, allocatable, intent(inout) :: iPair(:,:), img2CentCell(:)
    integer, intent(in) :: iSquare(:)
    type(TNeighborList), intent(inout) :: neighborList
    type(OSlakoCont), intent(in) :: skHamCont, skOverCont
    type(ORepCont), intent(in) :: pRepCont
    type(TOrbitals), intent(in) :: orb

    complex(cp) :: Ssqr(sf%nOrbs,sf%nOrbs), Sinv(sf%nOrbs,sf%nOrbs), T1(sf%nOrbs,sf%nOrbs)
    complex(cp) :: Rho(sf%nOrbs,sf%nOrbs,sf%nSpin), Rhoold(sf%nOrbs,sf%nOrbs,sf%nSpin)
    complex(cp) :: H1(sf%nOrbs,sf%nOrbs,sf%nSpin)
    complex(cp) :: Rhonew(sf%nOrbs,sf%nOrbs,sf%nSpin), RdotSprime(sf%nOrbs,sf%nOrbs)
    complex(cp), allocatable :: Eiginv(:,:,:), EiginvAdj(:,:,:)
    real(dp) :: qInput(orb%mOrb, sf%nAtom, sf%nSpin), qq(sf%nAtom,sf%nSpin), mu(3,sf%nSpin)
    real(dp) :: chargePerShell(orb%mShell,sf%nAtom,sf%nSpin), coordNew(3, sf%nAtom)
    real(dp), allocatable :: rhoPrim(:,:), ErhoPrim(:), ham0(:)
    real(dp) :: time, totalForce(3, sf%nAtom), dTime, startTim, ePerBond(sf%nAtom, sf%nAtom) 
    real(dp) :: movedAccel(3, sf%nMovedAtom), energyKin, new3Coord(3, sf%nMovedAtom)
    real(dp) :: timeInit = 0.0_dp, timeElec = 0.0_dp, timeIon = 0.0_dp, timeInver = 0.0_dp, startTime
    integer :: muDat, qDat, forceDat, energyDat, populDat(2), coorDat, tmpDat, ePBondDat
    integer :: iStep = 0, iAtom, iSpin, iTimeInit, iTimeElec, iTimeIon, iTimeLoop
    type(TPotentials) :: potential
    type(TEnergies) :: energy
    character(4) :: dumpIdx

    
    ! Initialize timer
    call tic(iTimeInit)
    ! Initialize stuff  
    call createMatrices(sf, Rho, H1, Ssqr, Sinv, H0, ham0, &
         & over, ham, Hsq, filling, orb, rhoPrim, ErhoPrim, potential, &
         & neighborList%iNeighbor, nNeighbor, iSquare, iPair, img2CentCell, & 
         & Eiginv, EiginvAdj, skOverCont)
    ! Initialize output files
    call initTDOutput(sf, muDat, qDat, energyDat, forceDat, populDat, coorDat, ePBondDat)
    startTime = 0.0_dp

    if (sf%Restart) then
       call readRestart(Rho, coord, sf%movedVelo, startTime)
       call updateH0S(sf, Ssqr, Sinv, coord, &
            & orb, neighborList, nNeighbor, iSquare, iPair, img2CentCell, &
            & skHamCont, skOverCont, ham, ham0, over, timeInver)
       if (sf%Ions) then
          sf%initialVelocities(:,:) = sf%movedVelo
          sf%ReadMDVelocities = .true.
       end if
    end if

    ! Initialice forces and nuclear dynamics
    energyKin = 0.0_dp
    call getChargeMu(sf, qq, qInput, mu, q0, Rho, Ssqr, coord, iSquare) !q_0
    call addSCCandField(sf, H1, ham, over, ham0, qInput, q0, coord, orb, potential, &  !H1_0
         &neighborList%iNeighbor, nNeighbor, iSquare, iPair, img2CentCell, iStep, chargePerShell, W)
    if (sf%Forces) then
       totalForce = 0.0_dp
       call getForces(sf, movedAccel, totalForce, Rho, H1, Sinv, neighborList,&    !F_0
           & nNeighbor, img2CentCell, iPair, iSquare, potential, orb, skHamCont, &
           & skOverCont, qInput, q0, pRepCont, coord, rhoPrim, ErhoPrim)
       if (sf%Ions) then
          call initIonDynamics(sf, coordNew, coord, movedAccel)
       end if
    end if

    ! Apply kick to Rho if necessary
    if (sf%doKick) then
       call kickDM(sf, Rho, Ssqr, Sinv, iSquare, coord)
    end if

    ! Initialize electron dinamics
    Rhoold(:,:,:) = Rho
    call euler(sf, sf%Dt, Rho, Rhoold, H1, Sinv, coord, skHamCont, &
         & skOverCont, orb, neighborList, nNeighbor, img2CentCell, iSquare)

    call getTDEnergy(sf, energy, energyKin, rhoPrim, Rhoold, neighborList%iNeighbor, &
         & nNeighbor, orb, iSquare, iPair, img2CentCell, ham0, qInput, q0, &
         & potential, chargePerShell, coord, pRepCont)

    if (sf%PairWiseEnergy) then
       call pairWiseBondEO(sf, ePerBond, rhoPrim(:,1), ham0, iSquare, &
            & neighborList%iNeighbor, nNeighbor, img2CentCell, iPair)
    end if

    call tac(timeInit, iTimeInit)
    ! End of initialization

    !! Main loop
    call tic(iTimeLoop)
    write(*,*)'Starting dynamics'

    do iStep = 0, sf%Nsteps
       time = iStep * sf%Dt + startTime
       ! Write everything to file
       if ((sf%Restart) .and. (iStep > 0) .or. (.not. sf%Restart)) then
          call writeTDOutputs(sf, muDat, qDat, energyDat, forceDat, populDat, coorDat, &
               & time, energy, energyKin, mu, qq, coord, totalForce, iStep, & 
               & ePerBond, ePBondDat)
       end if

       call tic(iTimeIon)
       if (sf%Ions) then
          coord(:,:) = coordNew
          call updateH0S(sf, Ssqr, Sinv, coord, & !S_1, Sinv_1, ham0_1
               & orb, neighborList, nNeighbor, iSquare, iPair, img2CentCell, &
               & skHamCont, skOverCont, ham, ham0, over, timeInver)
       end if
       call tac(dTime,iTimeIon)
       timeIon = timeIon + dTime
       
       if ((sf%PumpProbe) .and. (mod(iStep, sf%Nsteps/sf%PuProbeFrames) == 0)) then
         write(dumpIdx,'(i4)')int(sf%PuProbeFrames*iStep/sf%Nsteps)
         call writeRestart(Rho, coord, sf%movedVelo, 0.0_dp, &
             & trim(dumpIdx) // 'ppdump.bin')
       end if

       call getChargeMu(sf, qq, qInput, mu, q0, Rho, Ssqr, coord, iSquare) !q_1
       call addSCCandField(sf, H1, ham, over, ham0, qInput, q0, coord, orb, &
           & potential, neighborList%iNeighbor, nNeighbor, iSquare, iPair, &
           & img2CentCell, iStep, chargePerShell, W) ! H_1

       if ((sf%WriteRestart) .and. (iStep > 0) .and. &
            & (mod(iStep,sf%Nsteps/10) == 0)) then
          call writeRestart(Rho, coord, sf%movedVelo, time)    
       end if

       call tic(iTimeIon)
       if (sf%Forces) then
          call getForces(sf, movedAccel, totalForce, Rho, H1, Sinv, neighborList,&  !F_1
               & nNeighbor, img2CentCell, iPair, iSquare, potential, orb, skHamCont, &
               & skOverCont, qInput, q0, pRepCont, coord, rhoPrim, ErhoPrim)
       end if
 
       if (sf%Ions) then
          new3Coord(:,:) = coordNew(:, sf%indMovedAtom)
          call next(sf%pMDIntegrator, movedAccel, new3Coord, sf%movedVelo) !v_1, x_2 saved for later
          coordNew(:, sf%indMovedAtom) = new3Coord
          call getRdotSprime(sf, RdotSprime, coord, &
               & skHamCont, skOverCont, orb, img2CentCell, neighborList, &
               & nNeighbor, iSquare)
       end if
       call tac(dTime,iTimeIon)
       timeIon = timeIon + dTime

       call getTDEnergy(sf, energy, energyKin, rhoPrim, Rho, neighborList%iNeighbor, &
            & nNeighbor, orb, iSquare, iPair, img2CentCell, ham0, qInput, q0, &
            & potential, chargePerShell, coord, pRepCont) ! E_1

       if (sf%PairWiseEnergy) then
          call pairWiseBondEO(sf, ePerBond, rhoPrim(:,1), ham0, iSquare, &
               & neighborList%iNeighbor, nNeighbor, img2CentCell, iPair)
       end if

       do iSpin = 1, sf%nSpin
          call tic(iTimeElec)          
          if (sf%Ions) then
             call zscal(sf%nOrbs*sf%nOrbs, cmplx(0, 1, cp), H1(:,:,iSpin), 1)
             call zaxpy(sf%nOrbs*sf%nOrbs, 1.0_dp, RdotSprime, 1, H1(:,:,iSpin), 1)
          else
             call zscal(sf%nOrbs*sf%nOrbs, cmplx(0, 1, cp), H1(:,:,iSpin), 1)          
          end if

          if (sf%doEulers.and.(iStep > 100).and.(mod(iStep,sf%EulerEvery) == 0)) then
             call zcopy(sf%nOrbs*sf%nOrbs, Rho(:,:,iSpin), 1, Rhoold(:,:,iSpin), 1)
             call propagateRhoCPU(sf, Rhoold(:,:,iSpin), Rho(:,:,iSpin), &
                  & H1(:,:,iSpin), Sinv, T1, sf%Dt)
          else
             call propagateRhoCPU(sf, Rhoold(:,:,iSpin), Rho(:,:,iSpin), &
                  & H1(:,:,iSpin), Sinv, T1, 2.0_dp * sf%Dt)
          end if

          call zswap(sf%nOrbs*sf%nOrbs, Rhoold(:,:,iSpin), 1, Rho(:,:,iSpin), 1)

          call tac(dTime,iTimeElec)
          timeElec = timeElec + dTime

          if ((sf%Populations) .and. (mod(iStep, sf%SaveEvery) == 0)) then
             call getTDPopulations(sf, Rho, Eiginv, EiginvAdj, T1, Hsq, &
                  & populDat, time, iSpin, ham, over, neighborList%iNeighbor, &
                  & nNeighbor, iSquare, iPair, img2CentCell)
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
    write(*,*)'Time spent in forces + ion dynamics',timeIon
    write(*,*)'Time spent in overlap inversion', timeInver

    close(energyDat)
    close(qDat)
    close(muDat)
  end subroutine doDynamics



  !! Initialize ion dynamics
  subroutine initIonDynamics(sf, coordNew, coord, movedAccel)
    type(ElecDynamics), intent(inout) :: sf
    type(OVelocityVerlet), allocatable :: pVelocityVerlet
    real(dp), intent(in) :: coord(:,:), movedAccel(:,:)
    real(dp), intent(out) :: coordNew(:,:)
    logical :: halfVelocities = .true.
    real(dp) :: velocities(3, sf%nMovedAtom)

    allocate(pVelocityVerlet)
    if (sf%ReadMDVelocities) then
       call init(pVelocityVerlet, sf%Dt, coord(:, sf%indMovedAtom),&
            & sf%pThermostat, sf%initialVelocities, halfVelocities)
       sf%movedVelo(:, sf%indMovedAtom) = sf%initialVelocities
    else
       call init(pVelocityVerlet, sf%Dt, coord(:, sf%indMovedAtom), &
            &sf%pThermostat, halfVelocities, velocities)
       sf%movedVelo(:,:) = velocities
    end if

    ! Euler step forward
    sf%movedVelo(:,:) = sf%movedVelo - 0.5_dp * movedAccel * sf%Dt ! Has to be done for good initialization 
    coordNew(:,:) = coord
    coordNew(:,sf%indMovedAtom) = coord(:,sf%indMovedAtom) &
         & + sf%movedVelo(:,:) * sf%Dt + 0.5_dp * movedAccel(:,:) * sf%Dt**2

    ! This re-initializes the VVerlet propagator with coordNew
    sf%movedVelo(:,:) = sf%movedVelo + 0.5_dp * movedAccel * sf%Dt
    call init(pVelocityVerlet, sf%Dt, coordNew(:, sf%indMovedAtom),&
         & sf%pThermostat, sf%movedVelo, halfVelocities)
    allocate(sf%pMDIntegrator)
    call init(sf%pMDIntegrator, pVelocityVerlet)
  end subroutine initIonDynamics


  !! Calculates properties perbond.
  !! If hamover = ham0 is energy
  !! If hamover = over is bond order
  subroutine pairWiseBondEO(sf, EObond, rhoPrim, hamover, iSquare, iNeighbor, nNeighbor, &
       & img2CentCell, iPair) 
    type(ElecDynamics), intent(in), target :: sf
    real(dp), intent(in) :: rhoPrim(:)
    real(dp), intent(in) :: hamover(:)
    real(dp), intent(out) :: EObond(sf%nAtom, sf%nAtom)
    integer, intent(in) :: iNeighbor(0:,:), nNeighbor(:)
    integer, allocatable, intent(in) :: iPair(:,:), img2CentCell(:)
    integer, intent(in) :: iSquare(:)
    integer :: iAt1, iAt2, iAt2f, nOrb1, nOrb2, iOrig, iStart, iEnd, iNeigh, mOrb

    EObond = 0.0_dp
    
    do iAt1 = 1, sf%nAtom
       ii = iSquare(iAt1)
       nOrb1 = iSquare(iAt1+1) - ii
       do iNeigh = 0, nNeighbor(iAt1)
          iOrig = iPair(iNeigh,iAt1) + 1
          iAt2 = iNeighbor(iNeigh, iAt1)
          iAt2f = img2CentCell(iAt2)
          jj = iSquare(iAt2f)
          nOrb2 = iSquare(iAt2f+1) - jj

          EObond(iAt2, iAt1) = &
               & sum(rhoPrim(iOrig:iOrig+nOrb1*nOrb2-1) * &
               & hamover(iOrig:iOrig+nOrb1*nOrb2-1))

          EObond(iAt1, iAt2) = EObond(iAt2, iAt1)

       end do
    end do
  end subroutine pairWiseBondEO


  !! Calculates forces on nuclei and updates positions
  subroutine updateH0S(sf, Ssqr, Sinv, coord, &
       & orb, neighborList, nNeighbor, iSquare, iPair, img2CentCell, &
       & skHamCont, skOverCont, ham, ham0, over, timeInver)
    type(ElecDynamics), intent(inout), target :: sf 
    complex(cp), intent(inout) :: Sinv(:,:), Ssqr(:,:)
    real(dp), intent(inout), allocatable :: ham0(:)
    real(dp), allocatable, intent(out) :: ham(:,:), over(:)
    real(dp), allocatable, intent(inout) :: coord(:,:)

    type(TNeighborList), intent(inout) :: neighborList    
    integer, intent(inout) :: nNeighbor(:)
    integer, allocatable, intent(inout) :: iPair(:,:), img2CentCell(:)
    integer, intent(in) :: iSquare(:)
    type(TOrbitals), intent(in) :: orb
    type(OSlakoCont), intent(in) :: skHamCont, skOverCont

    real(dp) :: Sreal(sf%nOrbs,sf%nOrbs), SinvReal(sf%nOrbs,sf%nOrbs)
    real(dp) :: coord0Fold(3,sf%nAtom)
    integer :: nAllAtom, nAllOrb
    integer, allocatable :: specie0(:)
    integer :: info, ipiv(sf%nOrbs, sf%nOrbs)
    real(dp), intent(inout) :: timeInver
    real(dp) :: dTime
    integer :: iSpin, iInverTime, sparseSize
    
    !! Calculate overlap, H0 for new geometry
    nAllAtom = sf%nAtom
    coord0Fold(:,:) = coord
    allocate(specie0, source=sf%species)

    call updateNeighborListAndSpecies(coord, sf%species, img2CentCell, sf%iCellVec, &
         &neighborList, nAllAtom, coord0Fold, specie0, sf%mCutoff, sf%rCellVec)
    nAllOrb = sum(orb%nOrbSpecies(sf%species(1:nAllAtom)))
    call getNrOfNeighborsForAll(nNeighbor, neighborList, sf%skRepCutoff)
    call getSparseDescriptor(neighborList%iNeighbor, nNeighbor, img2CentCell, orb, iPair,&
         & sparseSize)

    deallocate(ham)
    deallocate(over)
    deallocate(ham0)
    allocate(ham(sparseSize, sf%nSpin))
    allocate(over(sparseSize))
    allocate(ham0(sparseSize))
    sf%nSparse = sparseSize

    call sf%sccCalc%updateCoords(coord, sf%species, neighborList, img2CentCell)
!   call updateCoords_SCC(coord, sf%species, neighborList, img2CentCell)

   if (sf%tDispersion) then
      call sf%dispersion%updateCoords(neighborList, img2CentCell, coord, &
           & specie0)
!     call updateCoords(sf%dispersion, neighborList, img2CentCell, coord, &
!         &specie0)
   end if
   
   call buildH0(ham0, skHamCont, sf%atomEigVal, coord, nNeighbor, neighborList%iNeighbor, &
        &sf%species, iPair, orb)
   call buildS(over, skOverCont, coord, nNeighbor, neighborList%iNeighbor, sf%species,&
        & iPair, orb)
!   call buildH0S(ham0, over, skHamCont, skOverCont, sf%atomEigVal, coord, &
!        &nNeighbor, neighborList%iNeighbor, sf%species, iPair, orb)

   call tic(iInverTime)
   Sreal = 0.0_dp
   call unpackHS(Sreal,over,neighborList%iNeighbor,nNeighbor,iSquare,iPair,img2CentCell)
   call blockSymmetrizeHS(Sreal,iSquare)
   Ssqr(:,:) = cmplx(Sreal, 0, cp)

   SinvReal = 0.0_dp
   do ii = 1, sf%nOrbs
      SinvReal(ii,ii) = 1.0_dp
   end do
   call gesv(Sreal,SinvReal)
   Sinv(:,:) = cmplx(SinvReal, 0, cp)
   !    call zgesv(sf%nOrbs, sf%nOrbs, Sreal, sf%nOrbs, ipiv, Sinv, sf%nOrbs, info)
   call tac(dTime, iInverTime)
   timeInver = timeInver + dTime

  end subroutine updateH0S


  !! Calculates force
  subroutine getForces(sf, movedAccel, totalForce, Rho, H1, Sinv, neighborList,&
       & nNeighbor, img2CentCell, iPair, iSquare, potential, orb, skHamCont, &
       & skOverCont, qInput, q0, pRepCont, coord, rhoPrim, ErhoPrim)
    type(ElecDynamics), intent(inout), target :: sf
    complex(cp), intent(in) :: Rho(:,:,:), H1(:,:,:), Sinv(:,:)
    type(TNeighborList), intent(in) :: neighborList    
    integer, intent(in) :: nNeighbor(:)
    integer, allocatable, intent(in) :: iPair(:,:), img2CentCell(:)
    integer, intent(in) :: iSquare(:)
    real(dp), intent(out) :: totalForce(3, sf%nAtom)
    real(dp), intent(out) :: movedAccel(3, sf%nMovedAtom)
    type(TPotentials), intent(in) :: potential
    type(TOrbitals), intent(in) :: orb
    type(OSlakoCont), intent(in) :: skHamCont, skOverCont
    real(dp), intent(in) :: qInput(:,:,:), q0(:,:,:)
    type(ORepCont), intent(in) :: pRepCont
    real(dp), intent(in) :: coord(:,:)
    real(dp), allocatable, intent(inout) :: rhoPrim(:,:), ErhoPrim(:)

    real(dp) :: T1(sf%nOrbs,sf%nOrbs),T2(sf%nOrbs,sf%nOrbs)
    real(dp) :: derivs(3,sf%nAtom), repulsiveDerivs(3,sf%nAtom), totalDeriv(3, sf%nAtom)
    integer :: iSpin
    iSpin = 1

    ! Eigenvectors stored in HSqrReal are overwritten
    call gemm(T1, real(H1(:,:,iSpin), dp), real(Rho(:,:,iSpin), dp),&
         & 1.0_dp , 0.0_dp, 'N', 'N', sf%nOrbs, sf%nOrbs, sf%nOrbs)
    call gemm(T2,  real(Sinv, dp), T1, 0.5_dp, 0.0_dp, 'N', 'N', sf%nOrbs, sf%nOrbs, sf%nOrbs)
    call daxpy(sf%nOrbs*sf%nOrbs, 1.0_dp, transpose(T2), 1, T2, 1)
!    T1 = T2 + transpose(T2)

    if (size(rhoPrim, dim=1) /= sf%nSparse) then
       deallocate(rhoPrim, ErhoPrim)
       allocate(rhoPrim(sf%nSparse, sf%nSpin))
       allocate(ErhoPrim(sf%nSparse))
    end if

    rhoPrim(:,iSpin) = 0.0_dp
    ErhoPrim(:) = 0.0_dp

    call packHS(rhoPrim(:,iSpin), real(Rho(:,:,iSpin), dp), neighborList%iNeighbor, & ! should be complex
         &nNeighbor, orb%mOrb, iSquare, iPair, img2CentCell)
    call packHS(ErhoPrim, T2, neighborList%iNeighbor, & ! should be complex
         &nNeighbor, orb%mOrb, iSquare, iPair, img2CentCell)

    derivs(:,:) = 0.0_dp 
    repulsiveDerivs(:,:) = 0.0_dp

    call derivative_shift(derivs,sf%derivator,rhoPrim,ErhoPrim(:),skHamCont, &
         &skOverCont, coord, sf%species, neighborList%iNeighbor, nNeighbor, &
         &img2CentCell, iPair, orb, potential%intBlock)
    call sf%sccCalc%updateCharges(qInput, q0, orb, sf%species, neighborList%iNeighbor, &
         &img2CentCell)
!    call updateCharges_SCC(qInput, q0, orb, sf%species, neighborList%iNeighbor, &
!         &img2CentCell)
    call sf%sccCalc%addForceDc(derivs, sf%species, neighborList%iNeighbor, img2CentCell, coord)
!    call addForceDCSCC(derivs, sf%species, neighborList%iNeighbor, img2CentCell, coord)
    call getERepDeriv(repulsiveDerivs, coord, nNeighbor, &
         &neighborList%iNeighbor,sf%species,pRepCont, img2CentCell)

    totalDeriv(:,:) = repulsiveDerivs(:,:) + derivs(:,:)
    if (sf%tDispersion) then
       call sf%dispersion%addGradients(totalDeriv)
    end if

    totalForce(:,:) = - totalDeriv(:,:)

    if (sf%Ions) then
       movedAccel(:,:) = totalForce(:, sf%indMovedAtom) / sf%movedMass
    else
       movedAccel(:,:) = totalForce / sf%movedMass
    end if
    
  end subroutine getForces


  !! Calculates nonadiabatic matrix (Sprime) times velocities (Rdot)
  subroutine getRdotSprime(sf, RdotSprime, coord,&
       & skHamCont, skOverCont, orb, img2CentCell, neighborList, &
       &nNeighbor, iSquare)
    type(ElecDynamics), intent(in), target :: sf
    type(OSlakoCont), intent(in) :: skHamCont, skOverCont
    complex(cp), intent(out) :: RdotSprime(sf%nOrbs,sf%nOrbs)
    type(TOrbitals), intent(in) :: orb
    real(dp) :: sPrimeTmp(orb%mOrb,orb%mOrb,3)
    real(dp) :: sPrimeTmp2(orb%mOrb,orb%mOrb), dcoord(3,sf%nAtom)
    real(dp), intent(in) :: coord(:,:)
    integer :: iAtom1,iStart1,iEnd1,iSp1,nOrb1,iAtomAux
    integer :: iNeigh,iStart2,iEnd2,iAtom2,iAtom2f,iSp2,nOrb2
    type(TNeighborList), intent(in) :: neighborList
    integer, intent(in) :: nNeighbor(:)
    integer, intent(in) :: iSquare(:)
    integer, allocatable, intent(in) :: img2CentCell(:)

    dcoord = 0.0_dp
    do ii=1,sf%nMovedAtom
       dcoord(:,sf%indMovedAtom(ii)) = sf%movedVelo(:,ii)
    end do

    RdotSprime = 0.0_cp

    !$OMP PARALLEL DO PRIVATE(iAtom1,iStart1,iEnd1,iSp1,nOrb1,sPrimeTmp2,iNeigh,iAtom2, &
    !$OMP& iAtom2f,iStart2,iEnd2,iSp2,nOrb2,sPrimeTmp,ii) DEFAULT(SHARED) &
    !$OMP& SCHEDULE(RUNTIME)
    do iAtom1 = 1, sf%nAtom
       iStart1 = iSquare(iAtom1)
       iEnd1 = iSquare(iAtom1+1)-1
       iSp1 = sf%species(iAtom1)
       nOrb1 = orb%nOrbAtom(iAtom1)

       ! Onsite blocks
       if (sf%CalcOnsiteGradients) then
          sPrimeTmp2(:,:) = 0.0_dp
          do ii=1,3
             sPrimeTmp2(:,:) = sPrimeTmp2 + &
                  &sf%onsiteGrads(ii,iAtom1,:,:) * dcoord(ii,iAtom1)
          end do
          RdotSprime(iStart1:iEnd1,iStart1:iEnd1) = cmplx(sPrimeTmp2(1:nOrb1,1:nOrb1), 0, cp)
       end if
    
       ! Offsite blocks
       do iNeigh = 1, nNeighbor(iAtom1)
          iAtom2 = neighborList%iNeighbor(iNeigh, iAtom1)
          iAtom2f = img2CentCell(iAtom2)
          iStart2 = iSquare(iAtom2f)
          iEnd2 = iSquare(iAtom2f+1)-1
          iSp2 = sf%species(iAtom2f)
          nOrb2 = orb%nOrbAtom(iAtom2f)
          if (iAtom2f /= iAtom1) then
             call sf%derivator%getFirstDeriv(sPrimeTmp, skOverCont, coord, sf%species,&
                  & iAtom1, iAtom2, orb)

             sPrimeTmp2(:,:) = 0.0_dp
             do ii=1,3
                sPrimeTmp2(:,:) = sPrimeTmp2 + &
                     &sPrimeTmp(:,:,ii) * dcoord(ii,iAtom1) ! - dcoord(ii,iAtom2))
             end do
             RdotSprime(iStart2:iEnd2,iStart1:iEnd1) = &
                  & cmplx(sPrimeTmp2(1:nOrb2,1:nOrb1), 0, cp)

             sPrimeTmp2(:,:) = 0.0_dp
             do ii=1,3
                sPrimeTmp2(:,:) = sPrimeTmp2 - &
                     &sPrimeTmp(:,:,ii) * dcoord(ii,iAtom2) ! - dcoord(ii,iAtom1))
             end do
             RdotSprime(iStart1:iEnd1,iStart2:iEnd2) = &
                  & cmplx(transpose(sPrimeTmp2(1:nOrb2,1:nOrb1)), 0, cp)
          end if
       end do
    end do
    !$OMP END PARALLEL DO
    
  end subroutine getRdotSprime
  

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

    call sf%sccCalc%updateCharges(qInput, q0, orb, sf%species, iNeighbor, img2CentCell)
    !    call updateCharges_SCC(qInput, q0, orb, sf%species, iNeighbor, img2CentCell)
    call sf%sccCalc%getShiftPerAtom(atomPot(:,1))
    call sf%sccCalc%getShiftPerL(shellPot(:,:,1))
    potential%intAtom(:,1) = potential%intAtom(:,1) + atomPot(:,1)
    potential%intShell(:,:,1) = potential%intShell(:,:,1) + shellPot(:,:,1)
    !    call getShiftPerAtom(potential%intAtom)
    !    call getShiftPerL(potential%intShell)
    
    !! Build spin contribution (if necessary)
    if (sf%SpinPol) then
       call getSpinShift(shellPot, chargePerShell, sf%species, orb, W)
       potential%intShell = potential%intShell + shellPot
!       call addSpinShift(potential%intShell, chargePerShell, sf%species, orb, W)
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
         & sf%species, orb, iPair, sf%nAtom, img2CentCell, potential%intShell)

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
    integer :: iStep, laserDat, laserDat2
    complex(cp) :: im

    im = cmplx(0, 1, cp)
    midPulse = (sf%Time0 + sf%Time1)/2.0_dp
    deltaT = sf%Time1 - sf%Time0
    angFreq = sf%Omega
    E0 = sf%Field
    sf%TDFunction(:,:) = 0.0_dp

    open(newunit=laserDat, file='laser.dat')
    if (sf%EnvFromFile) then
       open(newunit=laserDat2, file='final-laser.dat')
    end if

    do iStep = 0,sf%Nsteps
       time = iStep * sf%Dt

       if (sf%EnvType /= iTDConstant .and. (time >= sf%Time0) .and. (time <= sf%Time1)) then
          if (sf%EnvType == iTDGaussian) then
             env = Exp(-(time-midPulse)**2/(2.*248.**2))
          else if (sf%EnvType == iTDSin2) then
             env = Sin(pi*(time+sf%Time0)/deltaT)**2
          else if (sf%EnvType == iTDCos) then
             stop "Cosine envelope not yet implemented"
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
          write(laserDat2, "(5F15.8)") time * au__fs, &
               & TDFuncMod * (Hartree__eV / Bohr__AA), &
               & sf%TDFunction(iStep, :) * (Hartree__eV / Bohr__AA)
       else
          write(laserDat, "(5F15.8)") time * au__fs, &
               & TDFuncMod * (Hartree__eV / Bohr__AA), &
               & sf%TDFunction(iStep, :) * (Hartree__eV / Bohr__AA)
       end if
       
    end do

    close(laserDat)
    if (sf%EnvType /= iTDFromFile .and. sf%EnvFromFile) then
       close(laserDat2)
    end if
    
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
  subroutine getTDEnergy(sf, energy, energyKin, rhoPrim, Rho, iNeighbor, &
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
    real(dp), intent(out) :: energyKin
    real(dp), intent(in) :: coord(:,:)
    logical :: tUseBuggyRepSum = .false.
    type(ORepCont), intent(in) :: pRepCont

    iSpin = 1

    if (size(rhoPrim, dim=1) /= sf%nSparse) then
       deallocate(rhoPrim)
       allocate(rhoPrim(sf%nSparse, sf%nSpin))
    end if
    rhoPrim(:,iSpin) = 0.0_dp
    call packHS(rhoPrim(:,iSpin), real(Rho(:,:,iSpin), dp), iNeighbor, & ! should be complex
         &nNeighbor, orb%mOrb, iSquare, iPair, img2CentCell)

    call init(energy, sf%nAtom)
!    call create(energy, sf%nAtom)
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
    
    if (sf%tDispersion) then
       call sf%dispersion%getEnergies(energy%atomDisp)
       energy%eDisp = sum(energy%atomDisp)
    else
       energy%atomDisp(:) = 0.0_dp
       energy%eDisp = 0.0_dp
    end if

    energy%Eelec = energy%EnonSCC + energy%eSCC + energy%Espin + energy%Eext
    energy%Etotal = energy%Eelec + energy%Erep + energy%eDisp
    
    if (sf%Ions) then
      energyKin = 0.5_dp * sum(sf%movedMass(:,:)*sf%movedVelo(:,:)**2)
      energy%Etotal = energy%Etotal + energyKin
    end if
  end subroutine getTDEnergy


  !! Create all necessary matrices for dynamics
  subroutine createMatrices(sf, Rho, H1, Ssqr, Sinv, H0, ham0, &
       & over, ham, Hsq, filling, orb, rhoPrim, ErhoPrim, potential, &
       & iNeighbor, nNeighbor, iSquare, iPair, img2CentCell, Eiginv, &
       & EiginvAdj, skOverCont)
    type(ElecDynamics), intent(inout) :: sf
    real(dp), intent(inout) :: Hsq(:,:,:)
    real(dp), allocatable, intent(in) :: over(:), ham(:,:)
    real(dp), intent(in) :: filling(:,:,:), H0(:)
    integer, intent(in) :: iNeighbor(0:,:),nNeighbor(:)
    integer, intent(in) :: iSquare(:), iPair(0:,:), img2CentCell(:)
    type(TOrbitals), intent(in) :: orb
    type(TPotentials), intent(out) :: potential
    type(OSlakoCont), intent(in) :: skOverCont

    real(dp), allocatable, intent(out) :: rhoPrim(:,:), ErhoPrim(:), ham0(:)
    complex(cp), intent(out) :: Ssqr(:,:), Sinv(:,:), H1(:,:,:)
    complex(cp), intent(out) :: Rho(:,:,:)
    complex(cp), allocatable, intent(out), optional :: Eiginv(:,:,:), EiginvAdj(:,:,:)
    real(dp) :: T2(sf%nOrbs,sf%nOrbs), T3(sf%nOrbs, sf%nOrbs)
    integer :: iSpin

    allocate(rhoPrim(size(ham), sf%nSpin), ErhoPrim(size(ham)))
    sf%nSparse = size(H0)
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

    if (sf%Ions .and. sf%CalcOnsiteGradients) then
       allocate(sf%onsiteGrads(3, sf%nAtom, orb%mOrb, orb%mOrb))
       call getOnsiteGrads(sf, skOverCont, orb)
    end if

    call init(potential, orb, sf%nAtom, sf%nSpin)
!    call create(potential, orb, sf%nAtom, sf%nSpin)
  end subroutine createMatrices

  
  !! Calculates onsite gradients for non-adiabatic coupling
  subroutine getOnsiteGrads(sf, skOverCont, orb)
    type(ElecDynamics), intent(inout) :: sf
    type(OSlakoCont), intent(in) :: skOverCont
    type(TOrbitals), intent(in) :: orb    
    real(dp) :: dist, uVects(3,3), vect(3), Stmp(2, orb%mOrb, orb%mOrb), Sder(orb%mOrb, orb%mOrb)
    real(dp) :: Stmp2(3, orb%mOrb, orb%mOrb), Shorsf(orb%mOrb, orb%mOrb)
    real(dp) :: interSKOver(getMIntegrals(skOverCont))
    integer :: iAt, iSp, nOrb, dir

    dist = 0.02_dp
    uVects(1,:) = (/ 1.0_dp, 0.0_dp, 0.0_dp /)
    uVects(2,:) = (/ 0.0_dp, 1.0_dp, 0.0_dp /)
    uVects(3,:) = (/ 0.0_dp, 0.0_dp, 1.0_dp /)
    sf%onsiteGrads = 0.0_dp

    do iAt = 1, sf%nAtom
       Sder(:,:) = 0.0_dp
       Shorsf(:,:) = 0.0_dp
       iSp = sf%species(iAt)
       nOrb = orb%nOrbSpecies(iSp)
       do dir = 1, 3
          do ii = 1, 2
             Stmp(ii, :, :) = 0.0_dp
             vect = (-1)**(ii+1) * uVects(dir, :)

             call getSKIntegrals(skOverCont, interSKOver, dist, iSp, iSp)
             call rotateH0(Stmp(ii, :, :), interSKOver, vect(1), vect(2), vect(3), &
                  &iSp, iSp, orb)   
          end do
          Sder(:,:) = (Stmp(1, :, :) - Stmp(2, :, :)) / (2.0_dp * dist)
          sf%onsiteGrads(dir, iAt, 1:nOrb, 1:nOrb) = Sder(1:nOrb, 1:nOrb)
       end do
    end do

  end subroutine getOnsiteGrads


  !! Perfoms a step backwards to boot the dynamics using the Euler algorithm
  subroutine euler(sf, step, Rhoold, Rho, H1, Sinv, coord, &
       &skHamCont, skOverCont, orb, neighborList, nNeighbor, img2CentCell, iSquare)
    type(ElecDynamics), intent(in) :: sf
    complex(cp), intent(in) :: Rho(:,:,:), Sinv(:,:)
    complex(cp), intent(inout) :: H1(:,:,:) ! Warning! H1 is modified
    complex(cp), intent(out) :: Rhoold(:,:,:)
    complex(cp) :: RdotSprime(sf%nOrbs,sf%nOrbs)
    complex(cp) :: T1(sf%nOrbs,sf%nOrbs)
    real(dp), intent(in) :: coord(:,:)
    type(OSlakoCont), intent(in) :: skHamCont, skOverCont
    type(TOrbitals), intent(in) :: orb
    type(TNeighborList), intent(in) :: neighborList    
    integer, intent(in) :: nNeighbor(:), iSquare(:)
    integer, allocatable, intent(in) :: img2CentCell(:)
    real(dp), intent(in) :: step ! Step with its sign, to perform euler backwards (-) or forwards (+)
    integer :: iSpin

    Rhoold(:,:,:) = Rho

    if (sf%Ions) then
       call getRdotSprime(sf, RdotSprime, coord, &
            & skHamCont, skOverCont, orb, img2CentCell, neighborList, &
            & nNeighbor, iSquare)
    end if

    do iSpin=1,sf%nSpin
       T1 = 0.0_cp
       if (sf%Ions) then
          H1(:,:,iSpin) = RdotSprime + cmplx(0, 1, cp) * H1(:,:,iSpin)
       else
          H1(:,:,iSpin) = cmplx(0, 1, cp) * H1(:,:,iSpin)
       end if
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
  subroutine initTDOutput(sf, muDat, qDat, energyDat, forceDat, populDat, coorDat, ePBondDat)
    type(ElecDynamics), intent(in) :: sf
    integer, intent(out) :: muDat, qDat, energyDat, forceDat, populDat(2), coorDat, ePBondDat
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

    if (sf%Forces) then
       call openFile(sf, forceDat, 'forcesvst.dat')
    end if

    if (sf%Ions) then
       call openFile(sf, coorDat, 'tdcoords.xyz')
    end if
    
    if (sf%PairWiseEnergy) then
       if (sf%Restart) then
          inquire(file='eperbond.dat', exist=exist)
       end if
       if (exist) then
          open(newunit=ePBondDat, file='eperbond.dat', status="old", &
               &position="append", form='unformatted', access='stream')
       else
          open(newunit=ePBondDat, file='eperbond.dat', form='unformatted', access='stream')
       end if
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
  subroutine writeRestart(Rho, coord, veloc, time, dumpName)
    complex(cp), intent(in) :: Rho(:,:,:) 
    real(dp), intent(in) :: coord(:,:), time
    real(dp), intent(in) :: veloc(:,:)
    character(len=*), intent(in), optional :: dumpName
    integer :: dumpBin

    if (present(dumpName)) then
       open(newunit=dumpBin, file=dumpName, form='unformatted', access='stream', &
            & action='write')
    else
       open(newunit=dumpBin, file='tddump.bin', form='unformatted', access='stream', &
            & action='write')
    end if
    write(dumpBin) Rho, coord, veloc, time
    close(dumpBin)
  end subroutine writeRestart

  
  subroutine readRestart(Rho, coord, veloc, time)
    complex(cp), intent(out) :: Rho(:,:,:)
    real(dp), intent(out) :: coord(:,:), time
    real(dp), intent(out) :: veloc(:,:)
    integer :: dumpBin

    open(newunit=dumpBin, file='tddump.bin', form='unformatted', access='stream', &
         & action='read')
    read(dumpBin) Rho, coord, veloc, time
    close(dumpBin)
  end subroutine readRestart

  
  !! Write results to file
  subroutine writeTDOutputs(sf, muDat, qDat, energyDat, forceDat, populDat, coorDat, &
       & time, energy, energyKin, mu, qq, coord, totalForce, iStep, &
       & ePerBond, ePBondDat)
    type(ElecDynamics), intent(in) :: sf
    type(TEnergies), intent(in) :: energy
    integer, intent(in) :: muDat, qDat, energyDat, forceDat, populDat(2), coorDat
    real(dp), intent(in) :: time, mu(:,:), qq(:,:), coord(:,:), ePerBond(:,:)
    real(dp), intent(in) :: energyKin, totalForce(:,:)
    integer, intent(in) :: iStep, ePBondDat
    integer :: iAtom, iSpin
    real(dp) :: auxVeloc(3, sf%nAtom)

    write(energydat, '(9F25.15)') time * au__fs, energy%Etotal, energy%EnonSCC, &
         & energy%eSCC, energy%Espin, energy%Eext, energy%Erep, energyKin, energy%eDisp
    write(muDat, '(7F25.15)') time * au__fs, ((mu(ii, iSpin) * Bohr__AA, ii=1, 3), iSpin=1, sf%nSpin)
    if (mod(iStep,sf%SaveEvery) == 0) then
       write(qDat, '(5000F25.15)') time * au__fs, sum(qq(:,:)), (sum(qq(iAtom,:)), iAtom=1,sf%nAtom)
    end if

    if (sf%Ions .and. (mod(iStep,sf%SaveEvery) == 0)) then
       auxVeloc = 0.0_dp
       auxVeloc(:, sf%indMovedAtom) = sf%movedVelo
       write(coorDat,'(I5)')sf%nAtom
       write(coorDat,*) 'MD step:', iStep, 'time', time * au__fs
       do iAtom=1,sf%nAtom
          write(coorDat, '(A2, 6F16.8)') trim(sf%speciesName(sf%species(iAtom))), &
               &coord(:, iAtom) * Bohr__AA, auxVeloc(:, iAtom) * Bohr__AA / au__fs
       end do
    endif

    if (sf%Forces .and. (mod(iStep,sf%SaveEvery) == 0)) then
       write(forceDat, '(10000F25.15)') time * au__fs, (totalForce(:,iAtom), iAtom=1,sf%nAtom) 
    end if

    if (sf%PairWiseEnergy .and. mod(iStep,sf%SaveEvery) == 0) then
       write(ePBondDat) time * au__fs, sum(ePerBond(:,:)), &
            & ((ePerBond(ii, jj), ii=1,sf%nAtom), jj=1,sf%nAtom)
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
  subroutine getTDPopulations(sf, Rho, Eiginv, EiginvAdj, T1, Hsq, &
       & populDat, time, iSpin, ham, over, iNeighbor, &
       & nNeighbor, iSquare, iPair, img2CentCell)
    type(ElecDynamics), intent(in) :: sf
    complex(cp), intent(in) :: Rho(:,:,:)
    complex(cp), intent(inout) :: Eiginv(:,:,:), EiginvAdj(:,:,:)
    complex(cp), intent(inout) :: T1(:,:)
    complex(cp) :: T2(sf%nOrbs, sf%nOrbs)  ! Does it allocate every time?
    real(dp) :: T3(sf%nOrbs, sf%nOrbs)
    real(dp), intent(in) :: time
    integer, intent(in) :: populDat(2), iSpin
    integer, intent(in) :: iNeighbor(0:,:),nNeighbor(:)
    integer, intent(in) :: iSquare(:), iPair(0:,:), img2CentCell(:)
    real(dp), allocatable, intent(in) :: over(:), ham(:,:)
    real(dp), intent(inout) :: Hsq(:,:,:)
    real(dp) :: eigen(sf%nOrbs)

    T3 = 0.0_dp

    if (sf%Ions) then
       call diagonalize(HSq(:,:,iSpin), T3, eigen, &
            &ham(:,iSpin), over, iNeighbor, nNeighbor, &
            &iSquare, iPair, img2CentCell, 2, 'V') ! solver=2 (devide and conquer) 
       
       call TDPopulInit(sf, Eiginv, EiginvAdj, HSq, iSpin)
    end if

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
