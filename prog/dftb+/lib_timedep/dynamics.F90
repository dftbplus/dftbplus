!--------------------------------------------------------------------------------------------------!
!  DFTB+: general package for performing fast atomistic simulations                                !
!  Copyright (C) 2017  DFTB+ developers group                                                      !
!                                                                                                  !
!  See the LICENSE file for terms of usage and distribution.                                       !
!--------------------------------------------------------------------------------------------------!

#:include 'common.fypp'

!> Implements real-time time-dependent DFTB by numerically propagating the electronic one-electron
!> density matrix of the system in the presence of an external perturbation (kick or laser)
!>
!> 1) Oviedo, M. B., Negre, C. F. A. & Sánchez, C. G. Physical Chemistry Chemical Physics,
!> 12(25), 6706 (2010) https://doi.org/10.1039/b926051j
!> 2) Negre, C. F. A, Fuertes, V. C., Oviedo, M. B., Oliva, F. Y. & Sanchez, C. G.
!> Journal of Physical Chemistry C, 116(28), 14748–14753 (2012) https://doi.org/10.1021/jp210248k

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
  use environment
  implicit none
  private
  integer :: ii,jj,kk,ll

  public :: runDynamics, TElecDynamics_init
  public :: TElecDynamicsInp, TElecDynamics
  public :: iKick, iLaser, iNoTDPert
  public :: iTDConstant, iTDGaussian, iTDSin2, iTDFromFile
  public :: iTDSinglet, iTDTriplet

  !> Data type to  initialize electronic dynamics variables from parser
  type TElecDynamicsInp
    !> External field peak intensity
    real(dp) :: tdField

    !> Timestep for propagation
    real(dp) :: dt

    !> Field energy in eV
    real(dp) :: omega

    !> Real part of the polarization direction for laser fields
    real(dp) :: reFieldPolVec(3)

    !> Imaginary part of the polarization direction for laser fields
    real(dp) :: imFieldPolVec(3)

    !> Initial time of pulse
    real(dp) :: time0

    !> Final time of pulse
    real(dp) :: time1

    !> Phase applied of the laser field
    real(dp) :: phase

    !> Polarization direction of kicks
    integer :: polDir

    !> Number of steps to propagate
    integer :: steps

    !> Frequency of output writing for charges and populations
    integer :: writeFreq

    !> Type of external perturbation applied
    integer :: pertType

    !> Envelope shape for laser perturbation
    integer :: envType

    !> Spin type of absorption spectrum for spin polarised calculations
    integer :: spType

    !> Frequency of variable dumping to restart file
    integer :: restartFreq

    !> If time-dependent populations should be calculated
    logical :: tPopulations

    !> If calculation should be restarted from dump file
    logical :: tRestart

    !> If dump file should be written during the dynamics
    logical :: tWriteRestart

   integer, allocatable :: indMovedAtom(:), indExcitedAtom(:)
   integer :: nMovedAtom, tdEulerEvery, tdPPFrames, nExcitedAtom
   logical :: tdForces, tdIons, tReadMDVelocities, tddoEulers, tdPairWise, tdPumpProbe
   logical :: tdOnsiteGradients
   real(dp) :: tempAtom
   real(dp), allocatable :: initialVelocities(:,:)
end type TElecDynamicsInp

type TElecDynamics
   private
   real(dp) :: field, dt, omega, time0, time1
   complex(dp) :: fieldDir(3)
   real(dp), allocatable :: tdFunction(:, :), phase
   integer :: nSteps, writeFreq, pertType, envType, spType
   integer :: nAtom, nOrbs, nSpin=1, currPolDir=1, restartFreq
   integer, allocatable :: species(:), polDirs(:)
   character(mc), allocatable :: speciesName(:)
   logical :: tPopulations, tSpinPol=.false.
   logical :: tRestart, tWriteRestart, tWriteAutotest
   logical :: tLaser = .false., tKick = .false., tEnvFromFile = .false.
   type(TScc), allocatable :: sccCalc
   character(mc) :: autotestTag

   real(dp), allocatable :: initialVelocities(:,:)
   real(dp), allocatable :: movedVelo(:,:), movedMass(:,:), onsiteGrads(:,:,:,:)
   real(dp) :: mCutoff, skRepCutoff, rCellVec(3,1)
   real(dp), allocatable :: atomEigVal(:,:)
   integer :: nExcitedAtom, nMovedAtom, nSparse, EulerEvery, PuProbeFrames
   integer, allocatable :: iCellVec(:)
   integer, allocatable :: indMovedAtom(:), indExcitedAtom(:)
   logical :: Ions, Forces, tDispersion=.false., ReadMDVelocities, PumpProbe
   logical :: FirstIonStep = .true., doEulers = .false.
   logical :: PairWiseEnergy = .false., CalcOnsiteGradients = .false.
   type(OThermostat), allocatable :: pThermostat
   type(OMDIntegrator), allocatable :: pMDIntegrator
   class(DispersionIface), allocatable :: dispersion
   type(NonSccDiff), allocatable :: derivator
end type TElecDynamics

  !> Enumerating available types of perturbation
  integer, parameter :: iKick = 1, iLaser = 2, iNoTDPert = 3

  !> Enumerating available types of envelope function
  integer, parameter :: iTDConstant = 1, iTDGaussian = 2, iTDSin2 = 3, iTDFromFile = 4

  !> Enumerating available types of spin polarized spectra
  integer, parameter :: iTDSinglet = 1, iTDTriplet = 2


contains

  !! This initializes the input variables
    !> Initialisation of input variables
  subroutine TElecDynamics_init(this, inp, species, speciesName, tWriteAutotest, autotestTag,
    randomThermostat, mass, nAtom, skRepCutoff, mCutoff, iCellVec, atomEigVal, dispersion, nonSccDeriv)

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

    !> Tagged output files (machine readable)
    character(*), intent(in) :: autotestTag

    real(dp) :: norm

    real(dp), intent(in), allocatable :: atomEigVal(:,:)
    integer, intent(in) :: nAtom
    integer, allocatable, intent(in) :: iCellVec(:)
    type(ORanlux), allocatable, intent(inout) :: randomThermostat
    real(dp), intent(in) :: mCutoff, skRepCutoff, mass(:)
    class(DispersionIface), allocatable, intent(inout) :: dispersion
    type(NonSccDiff), intent(in) :: nonSccDeriv

    type(ODummyThermostat), allocatable :: pDummyTherm
    type(OMDCommon), allocatable :: pMDFrame
    real(dp) :: norm, tempAtom
    logical :: tMDstill, tDispersion
    complex(cp) :: im

    this%field = inp%tdfield
    this%dt = inp%dt
    this%nSteps = inp%steps
    this%pertType = inp%pertType
    this%envType = inp%envType
    this%spType = inp%spType
    this%tPopulations = inp%tPopulations
    this%tRestart = inp%tRestart
    this%tWriteRestart = inp%tWriteRestart
    this%phase = inp%phase
    this%writeFreq = inp%writeFreq
    this%restartFreq = inp%restartFreq
    this%species = species
    this%speciesName = speciesName
    allocate(this%sccCalc)

    if (inp%envType /= iTDConstant) then
      this%time0 = inp%time0
      this%time1 = inp%time1
    end if

    if (inp%pertType == iLaser) then
      this%tLaser = .true.
    else if (inp%pertType == iKick) then
      this%tKick = .true.
    end if

    if (this%tLaser) then
      this%omega = inp%omega
      this%fieldDir = inp%reFieldPolVec + imag * inp%imFieldPolVec
      norm = sqrt(sum(this%fieldDir * conjg(this%fieldDir)))
      ! normalize polarization vector
      this%fieldDir = this%fieldDir / norm
      allocate(this%tdFunction(0:this%nSteps, 3))
      this%tEnvFromFile = (this%envType == iTDFromFile)
      call getTDFunction(this)
    end if

    if (this%tKick) then
      if (inp%polDir == 4) then
        this%polDirs = [1, 2, 3]
      else
        this%polDirs = [inp%polDir]
      end if
    end if

    this%tWriteAutotest = tWriteAutotest
    this%autotestTag = autotestTag

    this%Ions = inp%tdIons
    this%Forces = inp%tdForces
    this%doEulers = inp%tddoEulers
    this%EulerEvery = inp%tdEulerEvery
    this%PairWiseEnergy = inp%tdPairWise
    this%PumpProbe = inp%tdPumpProbe
    this%PuProbeFrames = inp%tdPPFrames
    this%CalcOnsiteGradients = inp%tdOnsiteGradients

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
       allocate(this%derivator, source=nonSccDeriv)
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
  end subroutine TElecDynamics_init

  !> Driver of time dependent propagation to calculate wither spectrum or laser
  subroutine runDynamics(this, Hsq, ham, H0, q0, over, filling, neighbourList, nNeighbourSK,&
      & iSquare, iPair, img2CentCell, orb, coord, W, pRepCont, sccCalc, env, skHamCont, skOverCont)

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
    real(dp), allocatable, intent(in) :: W(:,:,:)

    !> occupations
    real(dp), intent(in) :: filling(:,:,:)

    !> Number of neighbours for each of the atoms
    integer, intent(inout) :: nNeighbourSK(:)

    !> index array for location of atomic blocks in large sparse arrays
    integer, allocatable, intent(inout) :: iPair(:,:)

    !> image atoms to their equivalent in the central cell
    integer, allocatable, intent(inout) :: img2CentCell(:)

    !> Index array for start of atomic block in dense matrices
    integer, intent(in) :: iSquare(:)

    !> list of neighbours for each atom
    type(TNeighbourList), intent(inout) :: neighbourList

    !> repulsive information
    type(ORepCont), intent(in) :: pRepCont

    !> Atomic orbital information
    type(TOrbitals), intent(in) :: orb

    !> SCC module internal variables
    type(TScc), intent(in) :: sccCalc

    !> Environment settings
    type(TEnvironment), intent(inout) :: env

    integer :: iPol
    logical :: tWriteAutotest
    type(OSlakoCont), intent(in) :: skHamCont, skOverCont

    this%sccCalc = sccCalc

    this%nSpin = size(ham(:,:), dim=2)
    if (this%nSpin > 1) then
      this%tSpinPol = .true.
      call qm2ud(q0)
    end if

    this%nOrbs = size(Hsq, dim=1)
    this%nAtom = size(coord, dim=2)

    tWriteAutotest = this%tWriteAutotest
    if (this%tKick) then
      do iPol = 1, size(this%polDirs)
        this%currPolDir = this%polDirs(iPol)
        ! Make sure only last component enters autotest
        tWriteAutotest = tWriteAutotest .and. (iPol == size(this%polDirs))
        call doDynamics(this, Hsq, ham, H0, q0, over, filling, neighbourList, nNeighbourSK,&
            & iSquare, iPair, img2CentCell, orb, coord, W, pRepCont, env, skHamCont, skOverCont)
      end do
    else
      call doDynamics(this, Hsq, ham, H0, q0, over, filling, neighbourList, nNeighbourSK, iSquare,&
          & iPair, img2CentCell, orb, coord, W, pRepCont, env, skHamCont, skOverCont)
    end if

  end subroutine runDynamics


  !> Runs the electronic dynamics of the system
  subroutine doDynamics(this, Hsq, ham, H0, q0, over, filling, neighbourList, nNeighbourSK,&
      & iSquare, iPair, img2CentCell, orb, coord, W, pRepCont, env)

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
    real(dp), allocatable, intent(in) :: W(:,:,:)

    !> occupations
    real(dp), intent(in) :: filling(:,:,:)

    !> Number of neighbours for each of the atoms
    integer, intent(inout) :: nNeighbourSK(:)

    !> index array for location of atomic blocks in large sparse arrays
    integer, allocatable, intent(inout) :: iPair(:,:)

    !> image atoms to their equivalent in the central cell
    integer, allocatable, intent(inout) :: img2CentCell(:)

    !> Index array for start of atomic block in dense matrices
    integer, intent(in) :: iSquare(:)

    !> list of neighbours for each atom
    type(TNeighbourList), intent(inout) :: neighbourList

    !> repulsive information
    type(ORepCont), intent(in) :: pRepCont

    !> Atomic orbital information
    type(TOrbitals), intent(in) :: orb

    !> Environment settings
    type(TEnvironment), intent(inout) :: env

    type(OSlakoCont), intent(in) :: skHamCont, skOverCont

    complex(dp) :: Ssqr(this%nOrbs,this%nOrbs), Sinv(this%nOrbs,this%nOrbs)
    complex(dp) :: rho(this%nOrbs,this%nOrbs,this%nSpin), rhoOld(this%nOrbs,this%nOrbs,this%nSpin)
    complex(dp) :: H1(this%nOrbs,this%nOrbs,this%nSpin), T1(this%nOrbs,this%nOrbs)
    complex(dp), allocatable :: Eiginv(:,:,:), EiginvAdj(:,:,:)
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

    complex(cp) :: RdotSprime(this%nOrbs,this%nOrbs)
    real(dp) :: coordNew(3, this%nAtom), totalForce(3, this%nAtom), ePerBond(this%nAtom, this%nAtom)
    real(dp) :: movedAccel(3, this%nMovedAtom), energyKin, new3Coord(3, this%nMovedAtom)
    real(dp), allocatable :: ErhoPrim(:)
    integer :: coorDat, ePBondDat

    call env%globalTimer%startTimer(globalTimers%elecDynInit)

    call initializeTDVariables(this, rho, H1, Ssqr, Sinv, H0, ham0, over, ham, Hsq, filling, orb,&
        & rhoPrim, potential, neighbourList%iNeighbour, nNeighbourSK, iSquare, iPair, img2CentCell,&
        & Eiginv, EiginvAdj, energy, ErhoPrim, skOverCont)

    call initTDOutput(this, dipoleDat, qDat, energyDat, populDat, forceDat, coorDat, ePBondDat)

    if (this%tRestart) then
       call readRestart(rho, Ssqr, coord, startTime)
       call updateH0S(this, Ssqr, Sinv, coord, &
            & orb, neighborList, nNeighbor, iSquare, iPair, img2CentCell, &
            & skHamCont, skOverCont, ham, ham0, over, timeInver)
       if (this%Ions) then
          this%initialVelocities(:,:) = this%movedVelo
          this%ReadMDVelocities = .true.
       end if
    end if

    ! Initialice forces and nuclear dynamics
    energyKin = 0.0_dp
    call getChargeDipole(this, deltaQ, qq, dipole, q0, rho, Ssqr, coord, iSquare)
    call updateH(this, H1, ham, over, ham0, qq, q0, coord, orb, potential,&
         & neighbourList%iNeighbour, nNeighbourSK, iSquare, iPair, img2CentCell, iStep,&
         & chargePerShell, W, env)

    if (this%Forces) then
       totalForce = 0.0_dp
       call getForces(this, movedAccel, totalForce, Rho, H1, Sinv, neighborList,&    !F_0
           & nNeighbor, img2CentCell, iPair, iSquare, potential, orb, skHamCont, &
           & skOverCont, qInput, q0, pRepCont, coord, rhoPrim, ErhoPrim)
       if (this%Ions) then
          call initIonDynamics(this, coordNew, coord, movedAccel)
       end if
    end if

    ! Apply kick to rho if necessary
    if (this%tKick) then
      call kickDM(this, rho, Ssqr, Sinv, iSquare, coord)
    end if

    ! Initialize electron dinamics
    rhoOld(:,:,:) = rho
    call initializePropagator(this, this%dt, rho, rhoOld, H1, Sinv)

    call getTDEnergy(this, energy, rhoPrim, rhoOld, neighbourList%iNeighbour, nNeighbourSK, orb,&
         & iSquare, iPair, img2CentCell, ham0, qq, q0, potential, chargePerShell, coord, &
         & pRepCont, energyKin)

    if (this%PairWiseEnergy) then
       call pairWiseBondEO(this, ePerBond, rhoPrim(:,1), ham0, iSquare, &
            & neighborList%iNeighbor, nNeighbor, img2CentCell, iPair)
    end if

    call env%globalTimer%stopTimer(globalTimers%elecDynInit)

    ! Main loop
    call env%globalTimer%startTimer(globalTimers%elecDynLoop)
    call loopTime%start()

    write(stdOut, "(A)") 'Starting dynamics'

    do iStep = 0, this%nSteps
      time = iStep * this%dt + startTime

      if (.not. this%tRestart .or. iStep > 0) then
        call writeTDOutputs(this, dipoleDat, qDat, energyDat, time, energy, dipole, deltaQ, iStep)
     end if


     if (this%Ions) then
        coord(:,:) = coordNew
        call updateH0S(this, Ssqr, Sinv, coord, & !S_1, Sinv_1, ham0_1
             & orb, neighborList, nNeighbor, iSquare, iPair, img2CentCell, &
             & skHamCont, skOverCont, ham, ham0, over, timeInver)
     end if

     if ((this%PumpProbe) .and. (mod(iStep, this%Nsteps/this%PuProbeFrames) == 0)) then
        write(dumpIdx,'(i4)')int(this%PuProbeFrames*iStep/this%Nsteps)
        call writeRestart(Rho, coord, this%movedVelo, 0.0_dp, &
             & trim(dumpIdx) // 'ppdump.bin')
     end if

     call getChargeDipole(this, deltaQ, qq, dipole, q0, rho, Ssqr, coord, iSquare)
     call updateH(this, H1, ham, over, ham0, qq, q0, coord, orb, potential,&
          & neighbourList%iNeighbour, nNeighbourSK, iSquare, iPair, img2CentCell, iStep,&
          & chargePerShell, W, env)

     if ((this%tWriteRestart) .and. (iStep > 0) .and. (mod(iStep, this%restartFreq) == 0)) then
        call writeRestart(rho, Ssqr, coord, time)
     end if

     call tic(iTimeIon)
     if (this%Forces) then
        call getForces(this, movedAccel, totalForce, Rho, H1, Sinv, neighborList,&  !F_1
             & nNeighbor, img2CentCell, iPair, iSquare, potential, orb, skHamCont, &
             & skOverCont, qInput, q0, pRepCont, coord, rhoPrim, ErhoPrim)
     end if

     if (this%Ions) then
        new3Coord(:,:) = coordNew(:, this%indMovedAtom)
        call next(this%pMDIntegrator, movedAccel, new3Coord, this%movedVelo) !v_1, x_2 saved for later
        coordNew(:, this%indMovedAtom) = new3Coord
        call getRdotSprime(this, RdotSprime, coord, skOverCont, orb, img2CentCell, &
             &neighborList, nNeighbor, iSquare)
     end if

     call getTDEnergy(this, energy, rhoPrim, rho, neighbourList%iNeighbour, nNeighbourSK, orb,&
          & iSquare, iPair, img2CentCell, ham0, qq, q0, potential, chargePerShell, coord, &
          & pRepCont, energyKin)

     if (this%PairWiseEnergy) then
        call pairWiseBondEO(this, ePerBond, rhoPrim(:,1), ham0, iSquare, &
             & neighborList%iNeighbor, nNeighbor, img2CentCell, iPair)
     end if

     do iSpin = 1, this%nSpin
        call tic(iTimeElec)
        if (this%Ions) then
           call zscal(this%nOrbs*this%nOrbs, cmplx(0, 1, cp), H1(:,:,iSpin), 1)
           call zaxpy(this%nOrbs*this%nOrbs, 1.0_dp, RdotSprime, 1, H1(:,:,iSpin), 1)
        else
           call zscal(this%nOrbs*this%nOrbs, cmplx(0, 1, cp), H1(:,:,iSpin), 1)
        end if

        if (this%doEulers.and.(iStep > 100).and.(mod(iStep,this%EulerEvery) == 0)) then
           call zcopy(this%nOrbs*this%nOrbs, Rho(:,:,iSpin), 1, Rhoold(:,:,iSpin), 1)
           call propagateRhoCPU(this, Rhoold(:,:,iSpin), Rho(:,:,iSpin), &
                & H1(:,:,iSpin), Sinv, T1, this%Dt)
        else
           call propagateRhoCPU(this, Rhoold(:,:,iSpin), Rho(:,:,iSpin), &
                & H1(:,:,iSpin), Sinv, T1, 2.0_dp * this%Dt)
        end if

        call zswap(this%nOrbs*this%nOrbs, Rhoold(:,:,iSpin), 1, Rho(:,:,iSpin), 1)

        if ((this%tPopulations) .and. (mod(iStep, this%writeFreq) == 0)) then
           call getTDPopulations(this, rho, Eiginv, EiginvAdj, T1, populDat, time, iSpin, &
           & H1, Ssqr, ham, over, neighborList%iNeighbor, nNeighbor, iSquare, iPair, img2CentCell)
        end if

     end do

       if (mod(iStep,this%Nsteps/10) == 0) then
          call tac(dTime,iTimeLoop)
          write(*,*)'Step ',iStep,'elapsed loop time:',dTime,'average time per loop',dTime/(iStep+1)
       end if
    end do

    if (mod(iStep, this%nSteps / 10) == 0) then
       call loopTime%stop()
       timeElec  = loopTime%getWallClockTime()
       write(stdOut, "(A,2x,I6,2(2x,A,F10.6))") 'Step ', iStep, 'elapsed loop time: ',&
            & timeElec, 'average time per loop ', timeElec / (iStep + 1)
    end if

    write(stdOut, "(A)") 'Dynamics finished OK!'
    call env%globalTimer%stopTimer(globalTimers%elecDynLoop)

    if (this%tWriteAutotest) then
       call writeTDAutotest(this, dipole, energy, deltaQ)
    end if

    call closeTDOutputs(this, dipoleDat, qDat, energyDat, populDat)

    if (this%Ions) then
       deallocate(this%pMDIntegrator)
    end if

  end subroutine doDynamics


  !> Updates the hamiltonian with SCC and external TD field (if any) contributions
  subroutine updateH(this, H1, ham, over, ham0, qq, q0, coord, orb, potential, iNeighbour,&
       & nNeighbourSK, iSquare, iPair, img2CentCell, iStep, chargePerShell, W, env)

    !> ElecDynamics instance
    type(TElecDynamics) :: this

    !> Square hamiltonian
    complex(dp), intent(inout) :: H1(:,:,:)

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
    integer, intent(in) :: iNeighbour(0:,:)

    !> Number of neighbours for each of the atoms
    integer, intent(in) :: nNeighbourSK(:)

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
    real(dp), allocatable, intent(in) :: W(:,:,:)

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

    call this%sccCalc%updateCharges(env, qq, q0, orb, this%species, iNeighbour, img2CentCell)
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
          potential%extAtom(iAtom, 1) = dot_product(coord(:,iAtom), this%tdFunction(iStep, :))
       end do
       call total_shift(potential%extShell, potential%extAtom, orb, this%species)
       call total_shift(potential%extBlock, potential%extShell, orb, this%species)
       potential%intShell = potential%intShell + potential%extShell ! for SCC
       potential%intBlock = potential%intBlock + potential%extBlock ! for forces
    end if

    call add_shift(ham, over, nNeighbourSK, iNeighbour, this%species, orb, iPair, this%nAtom,&
         & img2CentCell, potential%intShell)

    if (this%nSpin == 2) then
       ham(:,:) = 2.0_dp * ham
       call qm2ud(ham)
       call qm2ud(q0)
    end if

    do iSpin=1,this%nSpin
       call unpackHS(T2,ham(:,iSpin),iNeighbour,nNeighbourSK,iSquare,iPair,img2CentCell)
       call blockSymmetrizeHS(T2,iSquare)
       H1(:,:,iSpin) = cmplx(T2, 0.0_dp, dp)
    end do

  end subroutine updateH


  !> Clear potential object
  subroutine clearPotential(potential)

    !> potential acting on the system
    type(TPotentials), intent(inout) :: potential

    potential%intAtom(:,:) = 0.0_dp
    potential%intShell(:,:,:) = 0.0_dp
    potential%intBlock(:,:,:,:) = 0.0_dp
    potential%extAtom(:,:) = 0.0_dp
    potential%extShell(:,:,:) = 0.0_dp
    potential%extBlock(:,:,:,:) = 0.0_dp
  end subroutine clearPotential


  !! Kick the density matrix for spectrum calculations
  subroutine kickDM(this, Rho, Ssqr, Sinv, iSquare, coord)
    real(dp) :: pkick(this%nSpin)
    !> ElecDynamics instance
    type(TElecDynamics), intent(in) :: this

    !> Square overlap
    complex(dp), intent(in) :: Ssqr(:,:)

    !> Square overlap inverse
    complex(dp), intent(in) :: Sinv(:,:)

    !> Density matrix
    complex(dp), intent(inout) :: rho(:,:,:)

    !> atomic coordinates
    real(dp), allocatable, intent(in) :: coord(:,:)

    !> Index array for start of atomic block in dense matrices
    integer, intent(in) :: iSquare(:)

    complex(dp) :: T1(this%nOrbs, this%nOrbs, this%nSpin), T3(this%nOrbs, this%nOrbs, this%nSpin)
    complex(dp) :: T2(this%nOrbs, this%nOrbs), T4(this%nOrbs, this%nOrbs)
    integer :: iAt, iStart, iEnd, iSpin, iOrb
    real(dp) :: pkick(this%nSpin)

    character(1), parameter :: localDir(3) = ['x', 'y', 'z']

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
          do jj = iStart, iEnd
             T1(jj,jj,iSpin) = Exp(cmplx(0, -pkick(iSpin)*coord(this%PolDir, iAt), cp))
             T3(jj,jj,iSpin) = Exp(cmplx(0, pkick(iSpin)*coord(this%PolDir, iAt), cp))
          end do
       end do
    end do

    do iSpin=1,this%nSpin
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
  subroutine getTDFunction(this)
    !> ElecDynamics instance
    type(TElecDynamics), intent(inout) :: this

    real(dp) :: midPulse, deltaT, angFreq, E0, time, envelope
    real(dp) :: tdfun(3)
    integer :: iStep, laserDat

    im = cmplx(0, 1, cp)
    midPulse = (this%Time0 + this%Time1)/2.0_dp
    deltaT = this%Time1 - this%Time0
    angFreq = this%Omega
    E0 = this%Field
    this%TDFunction(:,:) = 0.0_dp

    open(newunit=laserDat, file='laser.dat')
    if (this%EnvFromFile) then
       open(newunit=laserDat2, file='final-laser.dat')
    end if

    do iStep = 0,this%Nsteps
       time = iStep * this%Dt

       if (this%EnvType /= iTDConstant .and. (time >= this%Time0) .and. (time <= this%Time1)) then
          if (this%EnvType == iTDGaussian) then
             env = Exp(-(time-midPulse)**2/(2.*248.**2))
          else if (this%EnvType == iTDSin2) then
             env = Sin(pi*(time+this%Time0)/deltaT)**2
          end if
       else if (this%EnvType == iTDConstant) then
          env = 1.0_dp
       else
          env = 0.0_dp
       end if

       if (this%EnvType /= iTDFromFile) then
          this%TDFunction(iStep, :) = E0 * env * aimag(Exp(im*(time*angFreq + this%phase)) * this%FieldDir)
          TDFuncMod = sqrt(sum(this%TDFunction(iStep, :)**2))
       end if

       if (this%EnvFromFile) then
          read(laserDat, *)time, TDFuncMod, tdfun(1), tdfun(2), tdfun(3)
          this%TDFunction(iStep, :) = this%TDFunction(iStep, :) + tdfun(:) * (Bohr__AA / Hartree__eV)
          write(laserDat2, "(5F15.8)") time * au__fs, &
               & TDFuncMod * (Hartree__eV / Bohr__AA), &
               & this%TDFunction(iStep, :) * (Hartree__eV / Bohr__AA)
       else
          write(laserDat, "(5F15.8)") time * au__fs, &
               & TDFuncMod * (Hartree__eV / Bohr__AA), &
               & this%TDFunction(iStep, :) * (Hartree__eV / Bohr__AA)
       end if

    end do

    close(laserDat)
    if (this%EnvType /= iTDFromFile .and. this%EnvFromFile) then
       close(laserDat2)
    end if

  end subroutine getTDFunction



  !! Calculate charges, dipole moments
  subroutine getChargeDipole(this, qq, qInput, dipole, q0, Rho, Ssqr, coord, iSquare)
        !> ElecDynamics instance
    type(TElecDynamics), intent(in) :: this

    !> Negative gross charge
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
    complex(dp), intent(in) :: rho(:,:,:)

    !> Square overlap matrix
    complex(dp), intent(in) :: Ssqr(:,:)

    !> Index array for start of atomic block in dense matrices
    integer, intent(in) :: iSquare(:)

    integer :: iAt, iOrb, iSpin, iOrb2

    qInput = 0.0_dp
    dipole = 0.0_dp

    do iSpin=1,this%nSpin
       do iAt = 1,this%nAtom
          iOrb = 0
          do jj = iSquare(iAt),iSquare(iAt+1)-1
             iOrb = iOrb + 1
             qInput(iOrb,iAt,iSpin) = sum(Rho(jj,:,iSpin)*Ssqr(jj,:)) ! 1 <= iOrb <= orb%mOrb
          end do
          dipole(1,iSpin) = dipole(1,iSpin) + sum(q0(:,iAt,iSpin)-qInput(:,iAt,iSpin))&
               & * coord(1,iAt)
          dipole(2,iSpin) = dipole(2,iSpin) + sum(q0(:,iAt,iSpin)-qInput(:,iAt,iSpin))&
               & * coord(2,iAt)
          dipole(3,iSpin) = dipole(3,iSpin) + sum(q0(:,iAt,iSpin)-qInput(:,iAt,iSpin))&
               & * coord(3,iAt)
       end do
    end do

    qq(:,:) = sum((qInput(:,:,:)-q0(:,:,:)), dim=1)
  end subroutine getChargeDipole


  !! Calculate energy
  subroutine getTDEnergy(this, energy, energyKin, rhoPrim, Rho, iNeighbor, &
       & nNeighbor, orb, iSquare, iPair, img2CentCell, ham0, qInput, q0, &
       & potential, chargePerShell, coord, pRepCont)

    !> ElecDynamics instance
    type(TElecDynamics), intent(inout) :: this

    !> data type for energy components and total
    type(TEnergies), intent(inout) :: energy

    !> sparse density matrix
    real(dp), allocatable, intent(inout) :: rhoPrim(:,:)

    !> Density matrix
    complex(dp), intent(in) :: rho(:,:,:)

    !> Sparse storage for non-SCC hamitonian
    real(dp), intent(in) :: ham0(:)

    !> atomic ocupations
    real(dp), intent(inout) :: qq(:,:,:)

    !> reference atomic occupations
    real(dp), intent(in) :: q0(:,:,:)

    !> Atomic orbital information
    type(TOrbitals), intent(in) :: orb

    !> Atomic neighbour data
    integer, intent(in) :: iNeighbour(0:,:)

    !> Number of neighbours for each of the atoms
    integer, intent(in) :: nNeighbourSK(:)

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

    real(dp), intent(out) :: energyKin
    logical :: tUseBuggyRepSum = .false.

    if (size(rhoPrim, dim=1) /= this%nSparse) then
       deallocate(rhoPrim)
       allocate(rhoPrim(this%nSparse, this%nSpin))
    end if

    do iSpin = 1, this%nSpin

       rhoPrim(:,iSpin) = 0.0_dp
       call packHS(rhoPrim(:,iSpin), real(Rho(:,:,iSpin), dp), iNeighbor, & ! should be complex
            &nNeighbor, orb%mOrb, iSquare, iPair, img2CentCell)
    end do

    if (this%nSpin == 2) then
       call ud2qm(rhoPrim)
    end if

    call init(energy, this%nAtom)
    energy%ETotal = 0.0_dp

    energy%EnonSCC = 0.0_dp
    energy%atomNonSCC(:) = 0.0_dp
    call mulliken(energy%atomNonSCC(:), rhoPrim(:,1), ham0, orb,&
         &iNeighbor, nNeighbor, img2CentCell, iPair)
    energy%EnonSCC =  sum(energy%atomNonSCC)

    if (this%doLaser) then ! energy in external field
       energy%atomExt = -sum(q0(:, :, 1) - qInput(:, :, 1),dim=1) &
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

    if (this%tDispersion) then
       call this%dispersion%getEnergies(energy%atomDisp)
       energy%eDisp = sum(energy%atomDisp)
    else
       energy%atomDisp(:) = 0.0_dp
       energy%eDisp = 0.0_dp
    end if

    energy%Eelec = energy%EnonSCC + energy%eSCC + energy%Espin + energy%Eext
    energy%Etotal = energy%Eelec + energy%Erep + energy%eDisp

    if (this%Ions) then
      energyKin = 0.5_dp * sum(this%movedMass(:,:)*this%movedVelo(:,:)**2)
      energy%Etotal = energy%Etotal + energyKin
    end if
  end subroutine getTDEnergy


  !> Create all necessary matrices and instances for dynamics
  subroutine initializeTDVariables(this, rho, H1, Ssqr, Sinv, H0, ham0, over, ham, Hsq, filling,&
      & orb, rhoPrim, potential, iNeighbour, nNeighbourSK, iSquare, iPair, img2CentCell, Eiginv,&
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
    integer, intent(in) :: iNeighbour(0:,:)

    !> Number of neighbours for each of the atoms
    integer, intent(in) :: nNeighbourSK(:)

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
    complex(dp), intent(out) :: Ssqr(:,:)

    !> Square overlap inverse
    complex(dp), intent(out) :: Sinv(:,:)

    !> Square hamiltonian
    complex(dp), intent(out) :: H1(:,:,:)

    !> Density matrix
    complex(dp), intent(out) :: rho(:,:,:)

    !> Inverse of eigenvectors matrix (for populations)
    complex(dp), allocatable, intent(out), optional :: Eiginv(:,:,:)

    !> Adjoint of the inverse of eigenvectors matrix (for populations)
    complex(dp), allocatable, intent(out), optional :: EiginvAdj(:,:,:)

    real(dp) :: T2(this%nOrbs,this%nOrbs), T3(this%nOrbs, this%nOrbs)
    integer :: iSpin, iOrb, iOrb2

    type(OSlakoCont), intent(in) :: skOverCont
    real(dp), allocatable, intent(out) :: ErhoPrim(:)

    allocate(rhoPrim(size(ham), this%nSpin), ErhoPrim(size(ham)))
    this%nSparse = size(H0)
    allocate(ham0(size(H0)))
    ham0(:) = H0

    T2 = 0.0_dp
    T3 = 0.0_dp
    call unpackHS(T2,over,iNeighbor,nNeighbor,iSquare,iPair,img2CentCell)
    call blockSymmetrizeHS(T2,iSquare)
    Ssqr(:,:) = cmplx(T2, 0, cp) ! Overlap

    do iSpin=1,this%nSpin
       call unpackHS(T3,ham(:,iSpin),iNeighbor,nNeighbor,iSquare,iPair,img2CentCell)
       call blockSymmetrizeHS(T3,iSquare)
       H1(:,:,iSpin) = cmplx(T3, 0, cp) ! Hamiltonian
       T3 = 0.0_dp
    end do

    if (this%Populations) then
       allocate(Eiginv(this%nOrbs, this%nOrbs, this%nSpin))
       allocate(EiginvAdj(this%nOrbs, this%nOrbs, this%nSpin))
       do iSpin=1,this%nSpin
          call TDPopulInit(this, Eiginv, EiginvAdj, HSq, iSpin)
       end do
    end if

    do ii = 1, this%nOrbs
       T3(ii,ii) = 1.0_dp
    end do
    call gesv(T2,T3)
    Sinv(:,:) = cmplx(T3, 0, cp) ! Inverse overlap
    write(*,*)'S inverted'

    do iSpin=1,this%nSpin
       T2 = 0.0_dp
       call makeDensityMatrix(T2, Hsq(:,:,iSpin), filling(:,1,iSpin))
       Rho(:,:,iSpin) = cmplx(T2, 0, cp) ! Density matrix
       do kk = 1, this%nOrbs-1
          do ll = kk+1, this%nOrbs
             Rho(kk,ll,iSpin) = Rho(ll,kk,iSpin) ! Symmetrization neccesary
          end do
       end do
    end do

    if (this%Ions .and. this%CalcOnsiteGradients) then
       allocate(this%onsiteGrads(3, this%nAtom, orb%mOrb, orb%mOrb))
       call getOnsiteGrads(this, skOverCont, orb)
    end if

    call init(potential, orb, this%nAtom, this%nSpin)

  end subroutine initializeTDVariables


  !> Perfoms a step backwards to boot the dynamics using the Euler algorithm
  subroutine initializePropagator(this, step, rhoOld, rho, H1, Sinv, &
    & coord, skHamCont, skOverCont, orb, neighborList, nNeighbor, img2CentCell, iSquare)
    !> ElecDynamics instance
    type(TElecDynamics), intent(in) :: this

    !> Density matrix
    complex(dp), intent(in) :: rho(:,:,:)

    !> Square overlap inverse
    complex(dp), intent(in) :: Sinv(:,:)

    !> Square hamiltonian
    complex(dp), intent(inout) :: H1(:,:,:)

    !> Density matrix at previous step
    complex(dp), intent(out) :: rhoOld(:,:,:)

    !> Time step in atomic units (with sign, to perform step backwards or forwards)
    real(dp), intent(in) :: step

    complex(dp) :: T1(this%nOrbs,this%nOrbs)
    integer :: iSpin

    complex(cp) :: RdotSprime(this%nOrbs,this%nOrbs)
    real(dp), intent(in) :: coord(:,:)
    type(OSlakoCont), intent(in) :: skHamCont, skOverCont
    type(TOrbitals), intent(in) :: orb
    type(TNeighborList), intent(in) :: neighborList
    integer, intent(in) :: nNeighbor(:), iSquare(:)
    integer, allocatable, intent(in) :: img2CentCell(:)

    Rhoold(:,:,:) = Rho

    if (this%Ions) then
       call getRdotSprime(this, RdotSprime, coord, skOverCont, orb, img2CentCell, &
            &neighborList, nNeighbor, iSquare)
    end if

    do iSpin=1,this%nSpin
       T1 = 0.0_cp
       if (this%Ions) then
          H1(:,:,iSpin) = RdotSprime + cmplx(0, 1, cp) * H1(:,:,iSpin)
       else
          H1(:,:,iSpin) = cmplx(0, 1, cp) * H1(:,:,iSpin)
       end if
       call propagateRho(this, Rhoold(:,:,iSpin), Rho(:,:,iSpin), &
            & H1(:,:,iSpin), Sinv, T1, step)
    end do

  end subroutine initializePropagator


  !> Propagate rho, notice that H = iH (coeficients are real)
  subroutine propagateRho(this, rhoOld, rho, H1, Sinv, T1, step)

    !> ElecDynamics instance
    type(TElecDynamics), intent(in) :: this

    !> Density matrix at previous step
    complex(dp), intent(inout) :: rhoOld(:,:)

    !> Density matrix
    complex(dp), intent(in) :: rho(:,:)

    !> Square hamiltonian
    complex(dp), intent(in) :: H1(:,:)

    !> Square overlap inverse
    complex(dp), intent(in) :: Sinv(:,:)

    !> Auxiliary matrix
    complex(dp), intent(out) :: T1(:,:)

    !> Time step in atomic units
    real(dp), intent(in) :: step

    call gemm(T1, Sinv, H1, cmplx(1, 0, cp), cmplx(0, 0, cp), &
         & 'N', 'N', this%nOrbs, this%nOrbs, this%nOrbs)
    call gemm(Rhoold, T1, Rho, cmplx(-step, 0, cp),&
         & cmplx(1, 0, cp), 'N', 'N', this%nOrbs, this%nOrbs, this%nOrbs)
    call gemm(Rhoold, Rho, T1, cmplx(-step, 0, cp),&
         & cmplx(1, 0, cp), 'N', 'C', this%nOrbs, this%nOrbs, this%nOrbs)
  end subroutine propagateRho


  !! Initialize output files
  subroutine initTDOutput(this, dipoleDat, qDat, energyDat, forceDat, populDat, coorDat, ePBondDat)
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

    integer, intent(out) :: forceDat, coorDat, ePBondDat

    character(20) :: dipoleFileName
    character(1) :: strSpin
    integer :: iSpin
    logical :: exist

    if (this%doKick) then
       if (this%PolDir == 1) then
          dipoleFileName = 'mux.dat'
       else if (this%PolDir == 2) then
          dipoleFileName = 'muy.dat'
       else if (this%PolDir == 3) then
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

    if (this%Forces) then
       call openFile(this, forceDat, 'forcesvst.dat')
    end if

    if (this%Ions) then
       call openFile(this, coorDat, 'tdcoords.xyz')
    end if

    if (this%PairWiseEnergy) then
       if (this%Restart) then
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


  !> Close output files
  subroutine closeTDOutputs(this, dipoleDat, qDat, energyDat, populDat, forceDat, coorDat, ePBondDat)

    !> ElecDynamics instance
    type(TElecDynamics), intent(in) :: this

    !> Dipole output file ID
    integer, intent(in) :: dipoleDat

    !> Charge output file ID
    integer, intent(in) :: qDat

    !> Energy output file ID
    integer, intent(in) :: energyDat

    !> Populations output file ID
    integer, intent(in) :: populDat(2)

    integer, intent(in) :: forceDat, coorDat, ePBondDat

    integer :: iSpin

    close(dipoleDat)
    close(qDat)
    close(energyDat)

    if (this%tPopulations) then
       do iSpin = 1, this%nSpin
          close(populDat(iSpin))
       end do
    end if

  end subroutine closeTDOutputs


  !! Open files in different ways deppending on their previous existance
  subroutine openFile(this, unitName, fileName)
    !> ElecDynamics instance
    type(TElecDynamics), intent(in) :: this

    !> File ID
    integer :: unitName

    !> Name of the file to open
    character(*) :: fileName

    logical :: exist=.false.

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
  subroutine writeRestart(Rho, coord, veloc, time, dumpName)

    complex(dp), intent(in) :: rho(:,:,:), Ssqr(:,:)

    !> atomic coordinates
    real(dp), intent(in) :: coord(:,:)

    !> elapsed simulated time in atomic units
    real(dp), intent(in) :: time

    !> name of the dump file
    character(len=*), intent(in), optional :: dumpName

    integer :: dumpBin

    complex(cp), intent(in) :: Rho(:,:,:)
    real(dp), intent(in) :: coord(:,:), time
    real(dp), intent(in) :: veloc(:,:)
    character(len=*), intent(in), optional :: dumpName
    integer :: dumpBin

    if (present(dumpName)) then
       open(newunit=dumpBin, file=dumpName, form='unformatted', access='stream', action='write')
    else
       open(newunit=dumpBin, file='tddump.bin', form='unformatted', access='stream', action='write')
    end if
    write(dumpBin) rho, coord, veloc, time
    close(dumpBin)
  end subroutine writeRestart


  subroutine readRestart(Rho, coord, veloc, time)

    !> Density Matrix
    complex(dp), intent(out) :: rho(:,:,:)

    !> atomic coordinates
    real(dp), intent(out) :: coord(:,:)

    !> Previous simulation elapsed time until restart file writing
    real(dp), intent(out) :: time

    real(dp), intent(out) :: veloc(:,:)
    integer :: dumpBin

    open(newunit=dumpBin, file='tddump.bin', form='unformatted', access='stream', &
         & action='read')
    read(dumpBin) Rho, coord, veloc, time
    close(dumpBin)

  end subroutine readRestart


  !! Write results to file
  subroutine writeTDOutputs(this, dipoleDat, qDat, energyDat, forceDat, populDat, coorDat, &
       & time, energy, energyKin, dipole, qq, coord, totalForce, iStep, &
       & ePerBond, ePBondDat)

    !> ElecDynamics instance
    type(TElecDynamics), intent(in) :: this

    !> data type for energy components and total
    type(TEnergies), intent(in) :: energy

    !> Dipole output file ID
    integer, intent(in) :: dipoleDat

    !> Charge output file ID
    integer, intent(in) :: qDat

    !> Energy output file ID
    integer, intent(in) :: energyDat

    !> Elapsed simulation time
    real(dp), intent(in) :: time

    !> Dipole moment
    real(dp), intent(in) :: dipole(:,:)

    !> Negative gross charge
    real(dp), intent(in) :: deltaQ(:,:)

    !> current step of the propagation
    integer, intent(in) :: iStep

    integer :: iAtom, iSpin, iDir

    integer, intent(in) :: forceDat, coorDat, ePBondDat
    real(dp), intent(in) :: coord(:,:), ePerBond(:,:)
    real(dp), intent(in) :: energyKin, totalForce(:,:)
    real(dp) :: auxVeloc(3, this%nAtom)

    write(energydat, '(9F25.15)') time * au__fs, energy%Etotal, energy%EnonSCC, &
         & energy%eSCC, energy%Espin, energy%Eext, energy%Erep, energyKin, energy%eDisp
    write(dipoleDat, '(7F25.15)') time * au__fs, ((dipole(ii, iSpin) * Bohr__AA, ii=1, 3), iSpin=1, this%nSpin)
    if (mod(iStep,this%SaveEvery) == 0) then
       write(qDat, '(5000F25.15)') time * au__fs, sum(qq(:,:)), (sum(qq(iAtom,:)), iAtom=1,this%nAtom)
    end if

    if (this%Ions .and. (mod(iStep,this%SaveEvery) == 0)) then
       auxVeloc = 0.0_dp
       auxVeloc(:, this%indMovedAtom) = this%movedVelo
       write(coorDat,'(I5)')this%nAtom
       write(coorDat,*) 'MD step:', iStep, 'time', time * au__fs
       do iAtom=1,this%nAtom
          write(coorDat, '(A2, 6F16.8)') trim(this%speciesName(this%species(iAtom))), &
               &coord(:, iAtom) * Bohr__AA, auxVeloc(:, iAtom) * Bohr__AA / au__fs
       end do
    endif

    if (this%Forces .and. (mod(iStep,this%SaveEvery) == 0)) then
       write(forceDat, '(10000F25.15)') time * au__fs, (totalForce(:,iAtom), iAtom=1,this%nAtom)
    end if

    if (this%PairWiseEnergy .and. mod(iStep,this%SaveEvery) == 0) then
       write(ePBondDat) time * au__fs, sum(ePerBond(:,:)), &
            & ((ePerBond(ii, jj), ii=1,this%nAtom), jj=1,this%nAtom)
    end if

  end subroutine writeTDOutputs


  !> Initialize matrices for populations
  subroutine tdPopulInit(this, Eiginv, EiginvAdj, HSq, iSpin)
    !> ElecDynamics instance
    type(TElecDynamics), intent(in) :: this

    !> Inverse of eigenvectors matrix (for populations)
    complex(dp), intent(out) :: Eiginv(:,:,:)

    !> Adjoint of the inverse of eigenvectors matrix (for populations)
    complex(dp), intent(out) :: EiginvAdj(:,:,:)

    !> Eigenvectors
    real(dp), intent(in) :: Hsq(:,:,:)

    !> Spin index
    integer, intent(in) :: iSpin

    real(dp), allocatable :: T2(:,:), T3(:,:)
    integer :: iOrb, iOrb2

    allocate(T2(this%nOrbs, this%nOrbs), T3(this%nOrbs, this%nOrbs))

    T2 = Hsq(:,:,iSpin)
    T3 = 0.0_dp
    do ii = 1, this%nOrbs
       T3(ii,ii) = 1.0_dp
    end do
    call gesv(T2,T3)
    Eiginv(:, :, iSpin) = cmplx(T3, 0, cp)

    T3 = 0.0_dp
    do ii = 1, this%nOrbs
       T3(ii, ii) = 1.0_dp
       do jj = 1, this%nOrbs
          T2(ii, jj) = Hsq(jj, ii, iSpin)
       end do
    end do
    call gesv(T2,T3)
    EiginvAdj(:, :, iSpin) = cmplx(T3, 0, cp)

    deallocate(T2, T3)

  end subroutine TDPopulInit


  !! Calculate populations at each time step
  subroutine getTDPopulations(this, Rho, Eiginv, EiginvAdj, T1, H1, Ssqr, &
       & populDat, time, iSpin, ham, over, iNeighbor, &
       & nNeighbor, iSquare, iPair, img2CentCell)
        !> ElecDynamics instance
    type(TElecDynamics), intent(in) :: this

    !> Density Matrix
    complex(dp), intent(in) :: rho(:,:,:)
    !> Inverse of eigenvectors matrix (for populations)

    complex(dp), intent(inout) :: Eiginv(:,:,:)

    !> Adjoint of the inverse of eigenvectors matrix (for populations)
    complex(dp), intent(inout) :: EiginvAdj(:,:,:)

    !> Auxiliary matrix
    complex(dp), intent(inout) :: T1(:,:)

    !> Elapsed simulation time
    real(dp), intent(in) :: time

    !> Populations output file ID
    integer, intent(in) :: populDat(2)

    !> Spin index
    integer, intent(in) :: iSpin

    complex(dp) :: T2(this%nOrbs, this%nOrbs)
    integer :: iOrb

    real(dp) :: T3(this%nOrbs, this%nOrbs)
    integer, intent(in) :: iNeighbor(0:,:),nNeighbor(:)
    integer, intent(in) :: iSquare(:), iPair(0:,:), img2CentCell(:)
    real(dp), allocatable, intent(in) :: over(:), ham(:,:)
    complex(cp), intent(inout) :: H1(:,:,:), Ssqr(:,:)
    real(dp) :: eigen(this%nOrbs)

    T3 = 0.0_dp

    if (this%Ions) then
       call diagDenseMtx(2, 'V', H1(:,:,iSpin), SSqr, eigen)
       call TDPopulInit(this, Eiginv, EiginvAdj, real(H1), iSpin)
    end if

    call gemm(T1, Eiginv(:,:,iSpin), Rho(:,:,iSpin), cmplx(1, 0, cp), &
         & cmplx(0, 0, cp), 'N', 'N', this%nOrbs, this%nOrbs, this%nOrbs)
    call gemm(T2, T1, EiginvAdj(:,:,iSpin), cmplx(1, 0, cp), &
         & cmplx(0, 0, cp), 'N', 'N', this%nOrbs, this%nOrbs, this%nOrbs)

    write(populDat(iSpin),'(5000F25.15)') time * au__fs, &
         & (real(T2(ii,ii), dp), ii=1, this%nOrbs)

  end subroutine getTDPopulations


  !> Write time-dependent tagged information to autotestTag file
  subroutine writeTDAutotest(this, dipole, energy, deltaQ)

    !> ElecDynamics instance
    type(TElecDynamics), intent(in) :: this

    !> Dipole moment
    real(dp), intent(in) :: dipole(:,:)

    !> data type for energy components and total
    type(TEnergies), intent(in) :: energy

    !> Negative gross charge
    real(dp), intent(in) :: deltaQ(:,:)

    integer :: fdAutotest

    open(newunit=fdAutotest, file=trim(this%autotestTag), position="append")

    call writeTagged(fdAutotest, tag_tdenergy, energy%eSCC)
    call writeTagged(fdAutotest, tag_tddipole, dipole)
    call writeTagged(fdAutotest, tag_tdcharges, deltaQ)

    close(fdAutotest)

  end subroutine writeTDAutotest


  !! Initialize ion dynamics
  subroutine initIonDynamics(this, coordNew, coord, movedAccel)
    type(ElecDynamics), intent(inout) :: this
    type(OVelocityVerlet), allocatable :: pVelocityVerlet
    real(dp), intent(in) :: coord(:,:), movedAccel(:,:)
    real(dp), intent(out) :: coordNew(:,:)
    logical :: halfVelocities = .true.
    real(dp) :: velocities(3, this%nMovedAtom)

    allocate(pVelocityVerlet)
    if (this%ReadMDVelocities) then
       call init(pVelocityVerlet, this%Dt, coord(:, this%indMovedAtom),&
            & this%pThermostat, this%initialVelocities, halfVelocities)
       this%movedVelo(:, this%indMovedAtom) = this%initialVelocities
    else
       call init(pVelocityVerlet, this%Dt, coord(:, this%indMovedAtom), &
            &this%pThermostat, halfVelocities, velocities)
       this%movedVelo(:,:) = velocities
    end if

    ! Euler step forward
    this%movedVelo(:,:) = this%movedVelo - 0.5_dp * movedAccel * this%Dt ! Has to be done for good initialization
    coordNew(:,:) = coord
    coordNew(:,this%indMovedAtom) = coord(:,this%indMovedAtom) &
         & + this%movedVelo(:,:) * this%Dt + 0.5_dp * movedAccel(:,:) * this%Dt**2

    ! This re-initializes the VVerlet propagator with coordNew
    this%movedVelo(:,:) = this%movedVelo + 0.5_dp * movedAccel * this%Dt
    call init(pVelocityVerlet, this%Dt, coordNew(:, this%indMovedAtom),&
         & this%pThermostat, this%movedVelo, halfVelocities)
    allocate(this%pMDIntegrator)
    call init(this%pMDIntegrator, pVelocityVerlet)
  end subroutine initIonDynamics


  !! Calculates forces on nuclei and updates positions
  subroutine updateH0S(this, Ssqr, Sinv, coord, &
       & orb, neighborList, nNeighbor, iSquare, iPair, img2CentCell, &
       & skHamCont, skOverCont, ham, ham0, over, timeInver)
    type(ElecDynamics), intent(inout), target :: this
    complex(cp), intent(inout) :: Sinv(:,:), Ssqr(:,:)
    real(dp), allocatable, intent(inout) :: ham0(:)
    real(dp), allocatable, intent(inout) :: ham(:,:), over(:)
    real(dp), allocatable, intent(inout) :: coord(:,:)

    type(TNeighborList), intent(inout) :: neighborList
    integer, intent(inout) :: nNeighbor(:)
    integer, allocatable, intent(inout) :: iPair(:,:), img2CentCell(:)
    integer, intent(in) :: iSquare(:)
    type(TOrbitals), intent(in) :: orb
    type(OSlakoCont), intent(in) :: skHamCont, skOverCont

    real(dp) :: Sreal(this%nOrbs,this%nOrbs), SinvReal(this%nOrbs,this%nOrbs)
    real(dp) :: coord0Fold(3,this%nAtom)
    integer :: nAllAtom, nAllOrb
    integer, allocatable :: specie0(:)
    integer :: info, ipiv(this%nOrbs, this%nOrbs)
    real(dp), intent(inout) :: timeInver
    real(dp) :: dTime
    integer :: iSpin, iInverTime, sparseSize

    !! Calculate overlap, H0 for new geometry
    nAllAtom = this%nAtom
    coord0Fold(:,:) = coord
    allocate(specie0, source=this%species)

    call updateNeighborListAndSpecies(coord, this%species, img2CentCell, this%iCellVec, &
         &neighborList, nAllAtom, coord0Fold, specie0, this%mCutoff, this%rCellVec)
    nAllOrb = sum(orb%nOrbSpecies(this%species(1:nAllAtom)))
    call getNrOfNeighborthisorAll(nNeighbor, neighborList, this%skRepCutoff)
    call getSparseDescriptor(neighborList%iNeighbor, nNeighbor, img2CentCell, orb, iPair,&
         & sparseSize)

    deallocate(ham)
    deallocate(over)
    deallocate(ham0)
    allocate(ham(sparseSize, this%nSpin))
    allocate(over(sparseSize))
    allocate(ham0(sparseSize))
    this%nSparse = sparseSize

    call this%sccCalc%updateCoords(this%env, coord, this%species, neighborList, img2CentCell)

   if (this%tDispersion) then
      call this%dispersion%updateCoords(neighborList, img2CentCell, coord, &
           & specie0)
   end if

   call buildH0(ham0, skHamCont, this%atomEigVal, coord, nNeighbor, neighborList%iNeighbor, &
        &this%species, iPair, orb)
   call buildS(over, skOverCont, coord, nNeighbor, neighborList%iNeighbor, this%species,&
        & iPair, orb)

   call tic(iInverTime)
   Sreal = 0.0_dp
   call unpackHS(Sreal,over,neighborList%iNeighbor,nNeighbor,iSquare,iPair,img2CentCell)
   call blockSymmetrizeHS(Sreal,iSquare)
   Ssqr(:,:) = cmplx(Sreal, 0, cp)

   SinvReal = 0.0_dp
   do ii = 1, this%nOrbs
      SinvReal(ii,ii) = 1.0_dp
   end do
   call gesv(Sreal,SinvReal)
   Sinv(:,:) = cmplx(SinvReal, 0, cp)
   call tac(dTime, iInverTime)
   timeInver = timeInver + dTime

  end subroutine updateH0S


  !! Calculates force
  subroutine getForces(this, movedAccel, totalForce, Rho, H1, Sinv, neighborList,&
       & nNeighbor, img2CentCell, iPair, iSquare, potential, orb, skHamCont, &
       & skOverCont, qInput, q0, pRepCont, coord, rhoPrim, ErhoPrim)
    type(ElecDynamics), intent(inout), target :: this
    complex(cp), intent(in) :: Rho(:,:,:), H1(:,:,:), Sinv(:,:)
    type(TNeighborList), intent(in) :: neighborList
    integer, intent(in) :: nNeighbor(:)
    integer, allocatable, intent(in) :: iPair(:,:), img2CentCell(:)
    integer, intent(in) :: iSquare(:)
    real(dp), intent(out) :: totalForce(3, this%nAtom)
    real(dp), intent(out) :: movedAccel(3, this%nMovedAtom)
    type(TPotentials), intent(in) :: potential
    type(TOrbitals), intent(in) :: orb
    type(OSlakoCont), intent(in) :: skHamCont, skOverCont
    real(dp), intent(inout) :: qInput(:,:,:), q0(:,:,:)
    type(ORepCont), intent(in) :: pRepCont
    real(dp), intent(in) :: coord(:,:)
    real(dp), allocatable, intent(inout) :: rhoPrim(:,:), ErhoPrim(:)

    real(dp) :: T1(this%nOrbs,this%nOrbs),T2(this%nOrbs,this%nOrbs)
    real(dp) :: derivs(3,this%nAtom), repulsiveDerivs(3,this%nAtom), totalDeriv(3, this%nAtom)
    integer :: iSpin

    if (size(rhoPrim, dim=1) /= this%nSparse) then
       deallocate(rhoPrim, ErhoPrim)
       allocate(rhoPrim(this%nSparse, this%nSpin))
       allocate(ErhoPrim(this%nSparse))
    end if

    ErhoPrim(:) = 0.0_dp

    do iSpin = 1, this%nSpin

       call gemm(T1, real(H1(:,:,iSpin), dp), real(Rho(:,:,iSpin), dp),&
            & 1.0_dp , 0.0_dp, 'N', 'N', this%nOrbs, this%nOrbs, this%nOrbs)
       call gemm(T2,  real(Sinv, dp), T1, 0.5_dp, 0.0_dp, 'N', 'N', this%nOrbs, this%nOrbs, this%nOrbs)
       call daxpy(this%nOrbs*this%nOrbs, 1.0_dp, transpose(T2), 1, T2, 1)
       !    T1 = T2 + transpose(T2)

       rhoPrim(:,iSpin) = 0.0_dp

       call packHS(rhoPrim(:,iSpin), real(Rho(:,:,iSpin), dp), neighborList%iNeighbor, & ! should be complex
            &nNeighbor, orb%mOrb, iSquare, iPair, img2CentCell)
       call packHS(ErhoPrim, T2, neighborList%iNeighbor, & ! should be complex
            &nNeighbor, orb%mOrb, iSquare, iPair, img2CentCell)
    end do

    if (this%nSpin == 2) then
       call ud2qm(qInput)
       call ud2qm(q0)
       call ud2qm(rhoPrim)
    end if

    derivs(:,:) = 0.0_dp
    repulsiveDerivs(:,:) = 0.0_dp

    call derivative_shift(derivs,this%derivator,rhoPrim,ErhoPrim(:),skHamCont, &
         &skOverCont, coord, this%species, neighborList%iNeighbor, nNeighbor, &
         &img2CentCell, iPair, orb, potential%intBlock)
    call this%sccCalc%updateCharges(this%env, qInput, q0, orb, this%species, neighborList%iNeighbor, &
         &img2CentCell)
    call this%sccCalc%addForceDc(this%env, derivs, this%species, neighborList%iNeighbor, &
         & img2CentCell, coord)
    call getERepDeriv(repulsiveDerivs, coord, nNeighbor, &
         &neighborList%iNeighbor,this%species,pRepCont, img2CentCell)

    totalDeriv(:,:) = repulsiveDerivs(:,:) + derivs(:,:)
    if (this%tDispersion) then
       call this%dispersion%addGradients(totalDeriv)
    end if

    totalForce(:,:) = - totalDeriv(:,:)

    if (this%Ions) then
       movedAccel(:,:) = totalForce(:, this%indMovedAtom) / this%movedMass
    else
       movedAccel(:,:) = totalForce / this%movedMass
    end if

    if (this%nSpin == 2) then
       call qm2ud(qInput)
       call qm2ud(q0)
    end if

  end subroutine getForces


  !! Calculates nonadiabatic matrix (Sprime) times velocities (Rdot)
  subroutine getRdotSprime(this, RdotSprime, coord, skOverCont, orb, img2CentCell, &
       &neighborList, nNeighbor, iSquare)
    type(ElecDynamics), intent(in), target :: this
    type(OSlakoCont), intent(in) :: skOverCont
    complex(cp), intent(out) :: RdotSprime(this%nOrbs,this%nOrbs)
    type(TOrbitals), intent(in) :: orb
    real(dp) :: sPrimeTmp(orb%mOrb,orb%mOrb,3)
    real(dp) :: sPrimeTmp2(orb%mOrb,orb%mOrb), dcoord(3,this%nAtom)
    real(dp), intent(in) :: coord(:,:)
    integer :: iAtom1,iStart1,iEnd1,iSp1,nOrb1,iAtomAux
    integer :: iNeigh,iStart2,iEnd2,iAtom2,iAtom2f,iSp2,nOrb2
    type(TNeighborList), intent(in) :: neighborList
    integer, intent(in) :: nNeighbor(:)
    integer, intent(in) :: iSquare(:)
    integer, allocatable, intent(in) :: img2CentCell(:)

    dcoord = 0.0_dp
    do ii=1,this%nMovedAtom
       dcoord(:,this%indMovedAtom(ii)) = this%movedVelo(:,ii)
    end do

    sPrimeTmp = 0.0_dp
    RdotSprime = 0.0_cp

    !$OMP PARALLEL DO PRIVATE(iAtom1,iStart1,iEnd1,iSp1,nOrb1,sPrimeTmp2,iNeigh,iAtom2, &
    !$OMP& iAtom2f,iStart2,iEnd2,iSp2,nOrb2,sPrimeTmp,ii) DEFAULT(SHARED) &
    !$OMP& SCHEDULE(RUNTIME)
    do iAtom1 = 1, this%nAtom
       iStart1 = iSquare(iAtom1)
       iEnd1 = iSquare(iAtom1+1)-1
       iSp1 = this%species(iAtom1)
       nOrb1 = orb%nOrbAtom(iAtom1)

       ! Onsite blocks
       if (this%CalcOnsiteGradients) then
          sPrimeTmp2(:,:) = 0.0_dp
          do ii=1,3
             sPrimeTmp2(:,:) = sPrimeTmp2 + &
                  &this%onsiteGrads(ii,iAtom1,:,:) * dcoord(ii,iAtom1)
          end do
          RdotSprime(iStart1:iEnd1,iStart1:iEnd1) = cmplx(sPrimeTmp2(1:nOrb1,1:nOrb1), 0, cp)
       end if

       ! Offsite blocks
       do iNeigh = 1, nNeighbor(iAtom1)
          iAtom2 = neighborList%iNeighbor(iNeigh, iAtom1)
          iAtom2f = img2CentCell(iAtom2)
          iStart2 = iSquare(iAtom2f)
          iEnd2 = iSquare(iAtom2f+1)-1
          iSp2 = this%species(iAtom2f)
          nOrb2 = orb%nOrbAtom(iAtom2f)
          if (iAtom2f /= iAtom1) then
             call this%derivator%getFirstDeriv(sPrimeTmp, skOverCont, coord, this%species,&
                  & iAtom1, iAtom2, orb)

             sPrimeTmp2(:,:) = 0.0_dp
             do ii=1,3
                sPrimeTmp2(:,:) = sPrimeTmp2 + &
                     &sPrimeTmp(:,:,ii) * dcoord(ii,iAtom1)
             end do
             RdotSprime(iStart2:iEnd2,iStart1:iEnd1) = &
                  & cmplx(sPrimeTmp2(1:nOrb2,1:nOrb1), 0, cp)

             sPrimeTmp2(:,:) = 0.0_dp
             do ii=1,3
                sPrimeTmp2(:,:) = sPrimeTmp2 - &
                     &sPrimeTmp(:,:,ii) * dcoord(ii,iAtom2)
             end do
             RdotSprime(iStart1:iEnd1,iStart2:iEnd2) = &
                  & cmplx(transpose(sPrimeTmp2(1:nOrb2,1:nOrb1)), 0, cp)
          end if
       end do
    end do
    !$OMP END PARALLEL DO

  end subroutine getRdotSprime


  !! Calculates onsite gradients for non-adiabatic coupling
  subroutine getOnsiteGrads(this, skOverCont, orb)
    type(ElecDynamics), intent(inout) :: this
    type(OSlakoCont), intent(in) :: skOverCont
    type(TOrbitals), intent(in) :: orb
    real(dp) :: dist, uVects(3,3), vect(3), Stmp(2, orb%mOrb, orb%mOrb), Sder(orb%mOrb, orb%mOrb)
    real(dp) :: Stmp2(3, orb%mOrb, orb%mOrb)
    real(dp) :: interSKOver(getMIntegrals(skOverCont))
    integer :: iAt, iSp, nOrb, dir

    dist = 0.02_dp
    uVects(1,:) = (/ 1.0_dp, 0.0_dp, 0.0_dp /)
    uVects(2,:) = (/ 0.0_dp, 1.0_dp, 0.0_dp /)
    uVects(3,:) = (/ 0.0_dp, 0.0_dp, 1.0_dp /)
    this%onsiteGrads = 0.0_dp

    do iAt = 1, this%nAtom
       Sder(:,:) = 0.0_dp
       iSp = this%species(iAt)
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
          this%onsiteGrads(dir, iAt, 1:nOrb, 1:nOrb) = Sder(1:nOrb, 1:nOrb)
       end do
    end do

  end subroutine getOnsiteGrads


  !! Calculates properties perbond.
  !! If hamover = ham0 is energy
  !! If hamover = over is bond order
  subroutine pairWiseBondEO(this, EObond, rhoPrim, hamover, iSquare, iNeighbor, nNeighbor, &
       & img2CentCell, iPair)
    type(ElecDynamics), intent(in), target :: this
    real(dp), intent(in) :: rhoPrim(:)
    real(dp), intent(in) :: hamover(:)
    real(dp), intent(out) :: EObond(this%nAtom, this%nAtom)
    integer, intent(in) :: iNeighbor(0:,:), nNeighbor(:)
    integer, allocatable, intent(in) :: iPair(:,:), img2CentCell(:)
    integer, intent(in) :: iSquare(:)
    integer :: iAt1, iAt2, iAt2f, nOrb1, nOrb2, iOrig, iStart, iEnd, iNeigh, mOrb

    EObond = 0.0_dp

    do iAt1 = 1, this%nAtom
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

end module timeprop_module
