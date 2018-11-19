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
!>
module timeprop_module
  use globalenv
  use commontypes
  use potentials
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
  use evaluateenergies
  use thirdorder_module, only : ThirdOrder
  use populations
  use repulsive
  use periodic
  use environment
  use timer
  use taggedoutput
  use hamiltonian
  use solvertypes
  use onsitecorrection
  implicit none
  private

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

    !> Electric field angular frequency
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
  end type TElecDynamicsInp

  !> Data type for electronic dynamics internal settings
  type TElecDynamics
    private
    real(dp) :: field, dt, omega, time0, time1
    complex(dp) :: fieldDir(3)
    real(dp), allocatable :: tdFunction(:, :), phase
    integer :: nSteps, writeFreq, pertType, envType, spType
    integer :: nAtom, nOrbs, nSpin=1, currPolDir=1, restartFreq
    integer, allocatable :: polDirs(:)
    logical :: tPopulations
    logical :: tRestart, tWriteRestart, tWriteAutotest
    logical :: tLaser = .false., tKick = .false., tEnvFromFile = .false.
    type(TScc), allocatable :: sccCalc
    character(mc) :: autotestTag
  end type TElecDynamics

  !> Enumerating available types of perturbation
  integer, parameter :: iKick = 1, iLaser = 2, iNoTDPert = 3

  !> Enumerating available types of envelope function
  integer, parameter :: iTDConstant = 1, iTDGaussian = 2, iTDSin2 = 3, iTDFromFile = 4

  !> Enumerating available types of spin polarized spectra
  integer, parameter :: iTDSinglet = 1, iTDTriplet = 2


contains


  !> Initialisation of input variables
  subroutine TElecDynamics_init(this, inp, tWriteAutotest, autotestTag)

    !> ElecDynamics instance
    type(TElecDynamics), intent(out) :: this

    !> ElecDynamicsInp instance
    type(TElecDynamicsInp), intent(in) :: inp

    !> produce tagged output?
    logical, intent(in) :: tWriteAutotest

    !> Tagged output files (machine readable)
    character(*), intent(in) :: autotestTag

    real(dp) :: norm

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
      norm = sqrt(real(dot_product(this%fieldDir,this%fieldDir)))
      ! normalize polarization vector
      this%fieldDir = this%fieldDir / norm
      allocate(this%tdFunction(3, 0:this%nSteps))
      this%tEnvFromFile = (this%envType == iTDFromFile)
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

  end subroutine TElecDynamics_init


  !> Driver of time dependent propagation to calculate wither spectrum or laser
  subroutine runDynamics(this, Hsq, ham, H0, species, q0, over, filling, neighbourList,&
      & nNeighbourSK, iSquare, iSparseStart, img2CentCell, orb, coord, spinW, pRepCont, sccCalc,&
      & env, tDualSpinOrbit, xi, thirdOrd, nDftbUFunc, UJ, nUJ, iUJ, niUJ, iHam,&
      & iAtInCentralRegion, tFixEf, Ef, onSiteElements)

    !> ElecDynamics instance
    type(TElecDynamics) :: this

    !> Eigenvectors
    real(dp), intent(inout) :: Hsq(:,:,:)

    !> Sparse non-SCC hamitonian
    real(dp), intent(in) :: H0(:)

    !> species of all atoms in the system
    integer, intent(in) :: species(:)

    !> reference atomic occupations
    real(dp), intent(inout) :: q0(:,:,:)

    !> resulting hamitonian (sparse)
    real(dp), allocatable, intent(inout) :: ham(:,:)

    !> overlap (sparse)
    real(dp), allocatable, intent(inout) :: over(:)

    !> atomic coordinates
    real(dp), allocatable, intent(inout) :: coord(:,:)

    !> spin constants
    real(dp), allocatable, intent(in) :: spinW(:,:,:)

    !> occupations
    real(dp), intent(in) :: filling(:,:,:)

    !> Number of neighbours for each of the atoms
    integer, intent(inout) :: nNeighbourSK(:)

    !> index array for location of atomic blocks in large sparse arrays
    integer, allocatable, intent(inout) :: iSparseStart(:,:)

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

    !> Is dual spin orbit being used (block potentials)
    logical, intent(in) :: tDualSpinOrbit

    !> Spin orbit constants if required
    real(dp), allocatable, intent(in) :: xi(:,:)

    !> 3rd order settings
    type(ThirdOrder), intent(inout), allocatable :: thirdOrd

    !> which DFTB+U functional (if used)
    integer, intent(in), optional :: nDftbUFunc

    !> U-J prefactors in DFTB+U
    real(dp), intent(in), allocatable :: UJ(:,:)

    !> Number DFTB+U blocks of shells for each atom type
    integer, intent(in), allocatable :: nUJ(:)

    !> which shells are in each DFTB+U block
    integer, intent(in), allocatable :: iUJ(:,:,:)

    !> Number of shells in each DFTB+U block
    integer, intent(in), allocatable :: niUJ(:,:)

    !> Imaginary part of sparse hamiltonian storage
    real(dp), allocatable, intent(inout) :: iHam(:,:)

    !> Atoms over which to sum the total energies
    integer, intent(in) :: iAtInCentralRegion(:)

    !> Whether fixed Fermi level(s) should be used. (No charge conservation!)
    logical, intent(in) :: tFixEf

    !> If tFixEf is .true. contains reservoir chemical potential, otherwise the Fermi levels found
    !> from the given number of electrons
    real(dp), intent(inout) :: Ef(:)

    !> Corrections terms for on-site elements
    real(dp), intent(in), allocatable :: onSiteElements(:,:,:,:)

    integer :: iPol
    logical :: tWriteAutotest

    this%sccCalc = sccCalc

    this%nSpin = size(ham(:,:), dim=2)
    if (this%nSpin > 1) then
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
        call doDynamics(this, Hsq, ham, H0, species, q0, over, filling, neighbourList,&
            & nNeighbourSK, iSquare, iSparseStart, img2CentCell, orb, coord, spinW, pRepCont, env,&
            & tDualSpinOrbit, xi, thirdOrd, nDftbUFunc, UJ, nUJ, iUJ, niUJ, iHam,&
            & iAtInCentralRegion, tFixEf, Ef, onSiteElements)
      end do
    else
      call doDynamics(this, Hsq, ham, H0, species, q0, over, filling, neighbourList, nNeighbourSK,&
          & iSquare, iSparseStart, img2CentCell, orb, coord, spinW, pRepCont, env, tDualSpinOrbit,&
          & xi, thirdOrd, nDftbUFunc, UJ, nUJ, iUJ, niUJ, iHam, iAtInCentralRegion, tFixEf, Ef,&
          & onSiteElements)
    end if

  end subroutine runDynamics


  !> Runs the electronic dynamics of the system
  subroutine doDynamics(this, Hsq, ham, H0, species, q0, over, filling, neighbourList,&
      & nNeighbourSK, iSquare, iSparseStart, img2CentCell, orb, coord, spinW, pRepCont, env,&
      & tDualSpinOrbit, xi, thirdOrd, nDftbUFunc, UJ, nUJ, iUJ, niUJ, iHam, iAtInCentralRegion,&
      & tFixEf, Ef, onSiteElements)

    !> ElecDynamics instance
    type(TElecDynamics) :: this

    !> Eigenvectors
    real(dp), intent(inout) :: Hsq(:,:,:)

    !> Sparse storage for non-SCC hamitonian
    real(dp), intent(in) :: H0(:)

    !> species of all atoms in the system
    integer, intent(in) :: species(:)

    !> reference atomic occupations
    real(dp), intent(inout) :: q0(:,:,:)

    !> resulting hamitonian (sparse)
    real(dp), allocatable, intent(inout) :: ham(:,:)

    !> overlap (sparse)
    real(dp), allocatable, intent(inout) :: over(:)

    !> atomic coordinates
    real(dp), allocatable, intent(inout) :: coord(:,:)

    !> spin constants
    real(dp), allocatable, intent(in) :: spinW(:,:,:)

    !> occupations
    real(dp), intent(in) :: filling(:,:,:)

    !> Number of neighbours for each of the atoms
    integer, intent(inout) :: nNeighbourSK(:)

    !> index array for location of atomic blocks in large sparse arrays
    integer, allocatable, intent(inout) :: iSparseStart(:,:)

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

    !> Is dual spin orbit being used (block potentials)
    logical, intent(in) :: tDualSpinOrbit

    !> which DFTB+U functional (if used)
    integer, intent(in), optional :: nDftbUFunc

    !> U-J prefactors in DFTB+U
    real(dp), intent(in), allocatable :: UJ(:,:)

    !> Number DFTB+U blocks of shells for each atom type
    integer, intent(in), allocatable :: nUJ(:)

    !> which shells are in each DFTB+U block
    integer, intent(in), allocatable :: iUJ(:,:,:)

    !> Number of shells in each DFTB+U block
    integer, intent(in), allocatable :: niUJ(:,:)

    !> Spin orbit constants if required
    real(dp), allocatable, intent(in) :: xi(:,:)

    !> 3rd order settings
    type(ThirdOrder), intent(inout), allocatable :: thirdOrd

    !> Imaginary part of sparse hamiltonian storage
    real(dp), allocatable, intent(inout) :: iHam(:,:)

    !> Atoms over which to sum the total energies
    integer, intent(in) :: iAtInCentralRegion(:)

    !> Whether fixed Fermi level(s) should be used. (No charge conservation!)
    logical, intent(in) :: tFixEf

    !> If tFixEf is .true. contains reservoir chemical potential, otherwise the Fermi levels found
    !> from the given number of electrons
    real(dp), intent(inout) :: Ef(:)

    !> Corrections terms for on-site elements
    real(dp), intent(in), allocatable :: onSiteElements(:,:,:,:)

    complex(dp) :: Ssqr(this%nOrbs,this%nOrbs), Sinv(this%nOrbs,this%nOrbs)
    complex(dp) :: rho(this%nOrbs,this%nOrbs,this%nSpin), rhoOld(this%nOrbs,this%nOrbs,this%nSpin)
    complex(dp) :: H1(this%nOrbs,this%nOrbs,this%nSpin), T1(this%nOrbs,this%nOrbs)
    complex(dp), allocatable :: Eiginv(:,:,:), EiginvAdj(:,:,:)
    real(dp) :: qq(orb%mOrb, this%nAtom, this%nSpin), deltaQ(this%nAtom,this%nSpin)
    real(dp) :: dipole(3,this%nSpin), chargePerShell(orb%mShell,this%nAtom,this%nSpin)
    real(dp), allocatable :: rhoPrim(:,:)
    real(dp) :: time, startTime = 0.0_dp, timeElec = 0.0_dp
    integer :: dipoleDat, qDat, energyDat, populDat(2)
    integer :: iStep = 0, iSpin
    type(TPotentials) :: potential
    type(TEnergies) :: energy
    type(TTimer) :: loopTime
    real(dp) :: TS(this%nSpin)
    real(dp), allocatable :: qBlock(:,:,:,:), qiBlock(:,:,:,:)

    call env%globalTimer%startTimer(globalTimers%elecDynInit)

    if (allocated(UJ) .or. allocated(onSiteElements)) then
      allocate(qBlock(orb%mOrb, orb%mOrb, this%nAtom, this%nSpin))
    end if

    if (this%tRestart) then
      call readRestart(rho, rhoOld, Ssqr, coord, startTime)
    end if

    if (this%tLaser) then
      call getTDFunction(this, startTime)
    end if

    call initializeTDVariables(this, rho, H1, Ssqr, Sinv, over, ham, Hsq, filling, orb,&
        & rhoPrim, potential, neighbourList%iNeighbour, nNeighbourSK, iSquare, iSparseStart,&
        & img2CentCell, Eiginv, EiginvAdj, energy)

    ! Calculate repulsive energy
    call getERep(energy%atomRep, coord, nNeighbourSK, neighbourList%iNeighbour, species, pRepCont,&
        & img2CentCell)
    energy%Erep = sum(energy%atomRep)

    call initTDOutput(this, dipoleDat, qDat, energyDat, populDat)

    rhoPrim(:,:) = 0.0_dp
    do iSpin = 1, this%nSpin
      call packHS(rhoPrim(:,iSpin), real(rho(:,:,iSpin), dp), neighbourList%iNeighbour,&
          & nNeighbourSK, orb%mOrb, iSquare, iSparseStart, img2CentCell)
    end do

    !call getChargeDipole(this, deltaQ, qq, dipole, q0, rho, Ssqr, coord, iSquare)
    if (allocated(qBlock)) then
      qBlock(:,:,:,:) = 0.0_dp
      do iSpin = 1, this%nSpin
        call mulliken(qBlock(:,:,:,iSpin), over, rhoPrim(:,iSpin), orb, neighbourList%iNeighbour,&
            & nNeighbourSK, img2CentCell, iSparseStart)
      end do
    end if
    qq(:,:,:) = 0.0_dp
    do iSpin = 1, this%nSpin
      call mulliken(qq(:,:,iSpin), over, rhoPrim(:,iSpin), orb, neighbourList%iNeighbour,&
          & nNeighbourSK, img2CentCell, iSparseStart)
    end do
    deltaQ(:,:) = sum((qq - q0), dim=1)
    dipole(:,:) = -matmul(coord(:,:), deltaQ(:,:))

    call updateH(this, H1, ham, over, H0, species, qq, q0, coord, orb, potential,&
        & neighbourList, nNeighbourSK, iSquare, iSparseStart, img2CentCell, iStep,&
        & chargePerShell, spinW, env, tDualSpinOrbit, xi, thirdOrd, iHam, qBlock, qiBlock,&
        & nDftbUFunc, UJ, nUJ, iUJ, niUJ, onSiteElements)

    ! Apply kick to rho if necessary
    if (this%tKick) then
      call kickDM(this, rho, Ssqr, Sinv, iSquare, coord)
    end if

    ! had to add the "or tKick" option to override rhoOld if tRestart = yes, otherwise it will be
    ! badly initialised
    if (.not.this%tRestart .or. this%tKick) then
      ! Initialize electron dynamics
      rhoOld(:,:,:) = rho
      call initializePropagator(this, this%dt, rho, rhoOld, H1, Sinv)
    end if

    rhoPrim(:,:) = 0.0_dp
    do iSpin = 1, this%nSpin
      call packHS(rhoPrim(:,iSpin), real(rho(:,:,iSpin), dp), neighbourList%iNeighbour,&
          & nNeighbourSK, orb%mOrb, iSquare, iSparseStart, img2CentCell)
    end do

    TS = 0.0_dp
    call getEnergies(this%sccCalc, qq, q0, chargePerShell, species, this%tLaser, .false.,&
        & .false., tDualSpinOrbit, rhoPrim, H0, orb, neighbourList, nNeighbourSK, img2CentCell,&
        & iSparseStart, 0.0_dp, 0.0_dp, TS, potential, energy, thirdOrd, qBlock, qiBlock,&
        & nDftbUFunc, UJ, nUJ, iUJ, niUJ, xi, iAtInCentralRegion, tFixEf, Ef, onSiteElements)

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

      rhoPrim(:,:) = 0.0_dp
      do iSpin = 1, this%nSpin
        call packHS(rhoPrim(:,iSpin), real(rho(:,:,iSpin), dp), neighbourList%iNeighbour,&
            & nNeighbourSK, orb%mOrb, iSquare, iSparseStart, img2CentCell)
      end do

      !call getChargeDipole(this, deltaQ, qq, dipole, q0, rho, Ssqr, coord, iSquare)
      if (allocated(qBlock)) then
        qBlock(:,:,:,:) = 0.0_dp
        do iSpin = 1, this%nSpin
          call mulliken(qBlock(:,:,:,iSpin), over, rhoPrim(:,iSpin), orb, neighbourList%iNeighbour,&
            & nNeighbourSK, img2CentCell, iSparseStart)
        end do
      end if
      qq(:,:,:) = 0.0_dp
      do iSpin = 1, this%nSpin
        call mulliken(qq(:,:,iSpin), over, rhoPrim(:,iSpin), orb, neighbourList%iNeighbour,&
            & nNeighbourSK, img2CentCell, iSparseStart)
      end do
      deltaQ(:,:) = sum((qq - q0), dim=1)
      dipole(:,:) = -matmul(coord(:,:), deltaQ(:,:))

      call updateH(this, H1, ham, over, H0, species, qq, q0, coord, orb, potential,&
          & neighbourList, nNeighbourSK, iSquare, iSparseStart, img2CentCell, iStep,&
          & chargePerShell, spinW, env, tDualSpinOrbit, xi, thirdOrd, iHam, qBlock, qiBlock,&
          & nDftbUFunc, UJ, nUJ, iUJ, niUJ, onSiteElements)

      if ((this%tWriteRestart) .and. (iStep > 0) .and. (mod(iStep, this%restartFreq) == 0)) then
        call writeRestart(rho, rhoOld, Ssqr, coord, time)
      end if

      rhoPrim(:,:) = 0.0_dp
      do iSpin = 1, this%nSpin
        call packHS(rhoPrim(:,iSpin), real(rho(:,:,iSpin), dp), neighbourList%iNeighbour,&
            & nNeighbourSK, orb%mOrb, iSquare, iSparseStart, img2CentCell)
      end do

      TS = 0.0_dp
      call getEnergies(this%sccCalc, qq, q0, chargePerShell, species, this%tLaser, .false.,&
          & .false., tDualSpinOrbit, rhoPrim, H0, orb, neighbourList, nNeighbourSK, img2CentCell,&
          & iSparseStart, 0.0_dp, 0.0_dp, TS, potential, energy, thirdOrd, qBlock, qiBlock,&
          & nDftbUFunc, UJ, nUJ, iUJ, niUJ, xi, iAtInCentralRegion, tFixEf, Ef, onSiteElements)

      do iSpin = 1, this%nSpin
        call scal(H1(:,:,iSpin), imag)
        call propagateRho(rhoOld(:,:,iSpin), rho(:,:,iSpin), H1(:,:,iSpin), Sinv, T1,&
            & 2.0_dp * this%dt)
        call swap(rhoOld(:,:,iSpin), rho(:,:,iSpin))

        if ((this%tPopulations) .and. (mod(iStep, this%writeFreq) == 0)) then
          call getTDPopulations(rho, Eiginv, EiginvAdj, T1, populDat, time, iSpin)
        end if
      end do

      if (mod(iStep, this%nSteps / 10) == 0) then
        call loopTime%stop()
        timeElec  = loopTime%getWallClockTime()
        write(stdOut, "(A,2x,I6,2(2x,A,F10.6))") 'Step ', iStep, 'elapsed loop time: ',&
            & timeElec, 'average time per loop ', timeElec / (iStep + 1)
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
  subroutine updateH(this, H1, ham, over, H0, species, qq, q0, coord, orb, potential,&
      & neighbourList, nNeighbourSK, iSquare, iSparseStart, img2CentCell, iStep, chargePerShell,&
      & spinW, env, tDualSpinOrbit, xi, thirdOrd, iHam, qBlock, qiBlock, nDftbUFunc, UJ, nUJ, iUJ,&
      & niUJ, onSiteElements)

    !> ElecDynamics instance
    type(TElecDynamics) :: this

    !> Square hamiltonian
    complex(dp), intent(inout) :: H1(:,:,:)

    !> resulting hamitonian (sparse)
    real(dp), allocatable, intent(inout) :: ham(:,:)

    !> overlap (sparse)
    real(dp), allocatable, intent(inout) :: over(:)

    !> Sparse storage for non-SCC hamitonian
    real(dp), intent(in) :: H0(:)

    !> species of all atoms in the system
    integer, intent(in) :: species(:)

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

    !> list of neighbours for each atom
    type(TNeighbourList), intent(inout) :: neighbourList

    !> Number of neighbours for each of the atoms
    integer, intent(in) :: nNeighbourSK(:)

    !> Index array for start of atomic block in dense matrices
    integer, intent(in) :: iSquare(:)

    !> index array for location of atomic blocks in large sparse arrays
    integer, intent(in) :: iSparseStart(0:,:)

    !> image atoms to their equivalent in the central cell
    integer, intent(in) :: img2CentCell(:)

    !> current step of the propagation
    integer, intent(in) :: iStep

    !> electrons in each atomi shell
    real(dp), intent(inout) :: chargePerShell(:,:,:)

    !> spin constants
    real(dp), allocatable, intent(in) :: spinW(:,:,:)

    !> Environment settings
    type(TEnvironment), intent(inout) :: env

    !> Is dual spin orbit being used (block potentials)
    logical, intent(in) :: tDualSpinOrbit

    !> Spin orbit constants if required
    real(dp), allocatable, intent(in) :: xi(:,:)

    !> 3rd order settings
    type(ThirdOrder), intent(inout), allocatable :: thirdOrd

    !> Imaginary part of sparse hamiltonian storage
    real(dp), allocatable, intent(inout) :: iHam(:,:)

    !> block (dual) atomic populations
    real(dp), intent(in), allocatable :: qBlock(:,:,:,:)

    !> Imaginary part of block atomic populations
    real(dp), intent(in), allocatable :: qiBlock(:,:,:,:)

    !> which DFTB+U functional (if used)
    integer, intent(in), optional :: nDftbUFunc

    !> U-J prefactors in DFTB+U
    real(dp), intent(in), allocatable :: UJ(:,:)

    !> Number DFTB+U blocks of shells for each atom type
    integer, intent(in), allocatable :: nUJ(:)

    !> which shells are in each DFTB+U block
    integer, intent(in), allocatable :: iUJ(:,:,:)

    !> Number of shells in each DFTB+U block
    integer, intent(in), allocatable :: niUJ(:,:)

    !> Corrections terms for on-site elements
    real(dp), intent(in), allocatable :: onSiteElements(:,:,:,:)

    real(dp) :: T2(this%nOrbs,this%nOrbs)
    integer :: iAtom, iSpin
    logical :: tDFTBU, tImHam

    ! left over from Poisson shift upload from transport being messy
    real(dp), allocatable :: dummy(:,:)

    ham = 0.0_dp
    if (this%nSpin == 2) then
      call ud2qm(qq)
      call ud2qm(q0)
    end if

    tDFTBU = .false.
    if (allocated(UJ)) then
      if (size(UJ) > 0) then
        tDFTBU = .true.
      end if
    end if
    tImHam = .false. ! for the moment

    call resetExternalPotentials(potential)
    call resetInternalPotentials(tDualSpinOrbit, xi, orb, species, potential)

    call getChargePerShell(qq, orb, species, chargePerShell)
    call addChargePotentials(env, this%sccCalc, qq, q0, chargePerShell, orb, species,&
        & neighbourList, img2CentCell, spinW, thirdOrd, potential, gammaf, .false., .false., dummy)
    if (allocated(UJ)) then
      call addBlockChargePotentials(qBlock, qiBlock, tDftbU, tImHam, species, orb, nDftbUFunc, UJ,&
          & nUJ, iUJ, niUJ, potential)
    end if

    if (allocated(onSiteElements)) then
      call addOnsShift(potential%intBlock, potential%iOrbitalBlock, qBlock, qiBlock, q0,&
          & onSiteElements, species, orb)
    end if

    ! Add time dependent field if necessary
    if (this%tLaser) then
      do iAtom = 1, this%nAtom
        potential%extAtom(iAtom, 1) = dot_product(coord(:,iAtom), this%tdFunction(:, iStep))
      end do
      call total_shift(potential%extShell, potential%extAtom, orb, species)
      call total_shift(potential%extBlock, potential%extShell, orb, species)
    end if

    potential%intBlock = potential%intBlock + potential%extBlock

    call getSccHamiltonian(H0, over, nNeighbourSK, neighbourList, species, orb, iSparseStart,&
        & img2CentCell, potential, ham, iHam)

    ! Hack due to not using Pauli-type structure outside of this part of the routine
    if (this%nSpin == 2) then
      ham = 2.0_dp * ham
      call qm2ud(ham)
      call qm2ud(q0)
    end if

    do iSpin=1,this%nSpin
      call unpackHS(T2, ham(:,iSpin), neighbourList%iNeighbour, nNeighbourSK, iSquare,&
          & iSparseStart, img2CentCell)
      call blockSymmetrizeHS(T2,iSquare)
      H1(:,:,iSpin) = cmplx(T2, 0.0_dp, dp)
    end do

  end subroutine updateH


  !> Kick the density matrix for spectrum calculations
  subroutine kickDM(this, rho, Ssqr, Sinv, iSquare, coord)

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

    pkick(1) = this%field

    if (this%nSpin == 2) then
      select case(this%spType)
      case (iTDSinglet)
        pkick(2) = pkick(1)
      case(iTDTriplet)
        pkick(2) = - pkick(1)
      end select
    end if

    T1(:,:,:) = 0.0_dp
    T2(:,:) = 0.0_dp
    T3(:,:,:) = 0.0_dp
    T4(:,:) = 0.0_dp

    do iSpin= 1, this%nSpin
      do iAt = 1, this%nAtom
        iStart = iSquare(iAt)
        iEnd = iSquare(iAt + 1) - 1
        do iOrb = iStart, iEnd
          T1(iOrb, iOrb, iSpin) = exp(cmplx(0, -pkick(iSpin) * coord(this%currPolDir, iAt), dp))
          T3(iOrb, iOrb, iSpin) = exp(cmplx(0,  pkick(iSpin) * coord(this%currPolDir, iAt), dp))
        end do
      end do
    end do

    do iSpin=1,this%nSpin
      call gemm(T2, T1(:,:,iSpin), rho(:,:,iSpin))
      call gemm(T4, T2, Ssqr, cmplx(1, 0, dp))
      call gemm(T2, T4, T3(:,:,iSpin))
      call gemm(rho(:,:,iSpin), T2, Sinv, cmplx(0.5, 0, dp))
      call gemm(rho(:,:,iSpin), Sinv, T2, cmplx(0.5, 0, dp), cmplx(1, 0, dp), 'N', 'C')
    end do

    write(stdout,"(A)")'Density kicked along ' // localDir(this%currPolDir) //'!'

  end subroutine kickDM


  !> Creates array for external TD field
  subroutine getTDFunction(this, startTime)

    !> ElecDynamics instance
    type(TElecDynamics), intent(inout) :: this

    !> starting time of the simulation, if relevant
    real(dp), intent(in) :: startTime

    real(dp) :: midPulse, deltaT, angFreq, E0, time, envelope
    real(dp) :: tdfun(3)
    integer :: iStep, laserDat

    midPulse = (this%time0 + this%time1)/2.0_dp
    deltaT = this%time1 - this%time0
    angFreq = this%omega
    E0 = this%field
    this%tdFunction(:,:) = 0.0_dp

    open(newunit=laserDat, file='laser.dat')

    do iStep = 0,this%nSteps
      time = iStep * this%dt + startTime

      if (this%envType == iTDConstant) then
        envelope = 1.0_dp
      else if (this%envType == iTDGaussian) then
        envelope = exp(-4.0_dp*pi*(time-midPulse)**2 / deltaT**2)
      else if (this%envType == iTDSin2 .and. (time >= this%time0) .and. (time <= this%time1)) then
        envelope = sin(pi*(time-this%time0)/deltaT)**2
      else
        envelope = 0.0_dp
      end if

      if (this%tEnvFromFile) then
        read(laserDat, *)time, tdfun(1), tdfun(2), tdfun(3)
        this%tdFunction(:, iStep) = tdfun * (Bohr__AA / Hartree__eV)
      else
        this%tdFunction(:, iStep) = E0 * envelope * aimag(exp(imag*(time*angFreq + this%phase))&
             & * this%fieldDir)
        write(laserDat, "(5F15.8)") time * au__fs,&
            & this%tdFunction(:, iStep) * (Hartree__eV / Bohr__AA)
      end if

    end do

    close(laserDat)

  end subroutine getTDFunction


  !> Calculate charges, dipole moments
  subroutine getChargeDipole(this, deltaQ, qq, dipole, q0, rho, Ssqr, coord, iSquare)

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

    integer :: iAt, iSpin, iOrb1, iOrb2

    qq = 0.0_dp
    do iSpin=1, this%nSpin
      do iAt = 1, this%nAtom
        iOrb1 = iSquare(iAt)
        iOrb2 = iSquare(iAt+1)-1
        ! hermitian transpose used as real part only is needed
        qq(:iOrb2-iOrb1+1,iAt,iSpin) = real(sum(rho(:, iOrb1:iOrb2, iSpin) * Ssqr(:, iOrb1:iOrb2),&
            & dim=1), dp)
      end do
    end do

    deltaQ(:,:) = sum((qq - q0), dim=1)
    dipole(:,:) = -matmul(coord(:,:), deltaQ(:,:))

  end subroutine getChargeDipole


  !> Create all necessary matrices and instances for dynamics
  subroutine initializeTDVariables(this, rho, H1, Ssqr, Sinv, over, ham, Hsq, filling,&
      & orb, rhoPrim, potential, iNeighbour, nNeighbourSK, iSquare, iSparseStart, img2CentCell,&
      & Eiginv, EiginvAdj, energy)

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

    !> Atomic neighbour data
    integer, intent(in) :: iNeighbour(0:,:)

    !> Number of neighbours for each of the atoms
    integer, intent(in) :: nNeighbourSK(:)

    !> Index array for start of atomic block in dense matrices
    integer, intent(in) :: iSquare(:)

    !> index array for location of atomic blocks in large sparse arrays
    integer, intent(in) :: iSparseStart(0:,:)

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

    !> Square overlap matrix
    complex(dp), intent(out) :: Ssqr(:,:)

    !> Square overlap inverse
    complex(dp), intent(out) :: Sinv(:,:)

    !> Square hamiltonian
    complex(dp), intent(out) :: H1(:,:,:)

    !> Density matrix
    complex(dp), intent(inout) :: rho(:,:,:)

    !> Inverse of eigenvectors matrix (for populations)
    complex(dp), allocatable, intent(out), optional :: Eiginv(:,:,:)

    !> Adjoint of the inverse of eigenvectors matrix (for populations)
    complex(dp), allocatable, intent(out), optional :: EiginvAdj(:,:,:)

    real(dp) :: T2(this%nOrbs,this%nOrbs), T3(this%nOrbs, this%nOrbs)
    integer :: iSpin, iOrb, iOrb2

    allocate(rhoPrim(size(ham, dim=1), this%nSpin))

    T2 = 0.0_dp
    T3 = 0.0_dp
    call unpackHS(T2,over,iNeighbour,nNeighbourSK,iSquare,iSparseStart,img2CentCell)
    call blockSymmetrizeHS(T2,iSquare)
    Ssqr(:,:) = cmplx(T2, 0, dp)

    do iSpin=1,this%nSpin
      call unpackHS(T3,ham(:,iSpin),iNeighbour,nNeighbourSK,iSquare,iSparseStart,img2CentCell)
      call blockSymmetrizeHS(T3,iSquare)
      H1(:,:,iSpin) = cmplx(T3, 0, dp)
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
    Sinv(:,:) = cmplx(T3, 0, dp)
    write(stdOut,"(A)")'S inverted'

    if (.not.this%tRestart) then
      do iSpin=1,this%nSpin
        T2 = 0.0_dp
        call makeDensityMatrix(T2,Hsq(:,:,iSpin),filling(:,1,iSpin))
        rho(:,:,iSpin) = cmplx(T2, 0, dp)
        do iOrb = 1, this%nOrbs-1
          do iOrb2 = iOrb+1, this%nOrbs
            rho(iOrb, iOrb2, iSpin) = rho(iOrb2, iOrb, iSpin)
          end do
        end do
      end do
    end if

    call init(potential, orb, this%nAtom, this%nSpin)
    call init(energy, this%nAtom)

  end subroutine initializeTDVariables


  !> Perfoms a step backwards to boot the dynamics using the Euler algorithm
  subroutine initializePropagator(this, step, rhoOld, rho, H1, Sinv)

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

    rhoOld(:,:,:) = rho

    do iSpin=1,this%nSpin
      T1 = 0.0_dp
      H1(:,:,iSpin) = imag * H1(:,:,iSpin)
      call propagateRho(rhoOld(:,:,iSpin), rho(:,:,iSpin), H1(:,:,iSpin), Sinv, T1, step)
    end do

  end subroutine initializePropagator


  !> Propagate rho, notice that H = iH (coeficients are real)
  subroutine propagateRho(rhoOld, rho, H1, Sinv, T1, step)

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

    call gemm(T1, Sinv, H1)
    call gemm(rhoOld, T1, rho, cmplx(-step, 0, dp), cmplx(1, 0, dp))
    call gemm(rhoOld, rho, T1, cmplx(-step, 0, dp), cmplx(1, 0, dp), 'N', 'C')

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
    integer, intent(in) :: dipoleDat

    !> Charge output file ID
    integer, intent(in) :: qDat

    !> Energy output file ID
    integer, intent(in) :: energyDat

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

  end subroutine closeTDOutputs


  !> Open files in different ways deppending on their previous existance
  subroutine openFile(this, unitName, fileName)

    !> ElecDynamics instance
    type(TElecDynamics), intent(in) :: this

    !> File ID
    integer :: unitName

    !> Name of the file to open
    character(*) :: fileName

    character(30) :: newName

    character(1) :: strCount

    logical :: exist=.false.

    integer :: count

    ! changed the append by this block to rename the restarted output
    if (this%tRestart) then
      inquire(file=fileName, exist=exist)
      count = 1
      do while (exist)
        write(strCount,'(i1)') count
        newName = "rest" // trim(strCount) // "_" // fileName
        inquire(file=newName, exist=exist)
        count = count + 1
     end do
    else
      newName = fileName
    end if

    open(newunit=unitName, file=newName, action="write")

  end subroutine openFile


  !> Write to and read from restart files
  subroutine writeRestart(rho, rhoOld, Ssqr, coord, time, dumpName)

    complex(dp), intent(in) :: rho(:,:,:), rhoOld(:,:,:), Ssqr(:,:)
    !> atomic coordinates
    real(dp), intent(in) :: coord(:,:)

    !> elapsed simulated time in atomic units
    real(dp), intent(in) :: time

    !> name of the dump file
    character(len=*), intent(in), optional :: dumpName

    integer :: dumpBin

    if (present(dumpName)) then
      open(newunit=dumpBin, file=dumpName, form='unformatted', access='stream', action='write')
    else
      open(newunit=dumpBin, file='tddump.bin', form='unformatted', access='stream', action='write')
    end if
    write(dumpBin) rho, rhoOld, Ssqr, coord, time
    close(dumpBin)
  end subroutine writeRestart


  !> read a restart file containing density matrix, overlap, coordinates and time step
  subroutine readRestart(rho, rhoOld, Ssqr, coord, time)

    !> Density Matrix
    complex(dp), intent(out) :: rho(:,:,:)

    !> Previous density Matrix
    complex(dp), intent(out) :: rhoOld(:,:,:)

    !> Square overlap matrix
    complex(dp), intent(out) :: Ssqr(:,:)

    !> atomic coordinates
    real(dp), intent(out) :: coord(:,:)

    !> Previous simulation elapsed time until restart file writing
    real(dp), intent(out) :: time

    integer :: dumpBin

    open(newunit=dumpBin, file='tddump.bin', form='unformatted', access='stream', action='read')
    read(dumpBin) rho, rhoOld, Ssqr, coord, time
    close(dumpBin)

  end subroutine readRestart


  !> Write results to file
  subroutine writeTDOutputs(this, dipoleDat, qDat, energyDat, time, energy, dipole, deltaQ, iStep)

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

    write(energydat, '(9F25.15)') time * au__fs, energy%Etotal, energy%EnonSCC, energy%eSCC,&
        & energy%Espin, energy%Eext, energy%Erep

    write(dipoleDat, '(7F25.15)') time * au__fs, ((dipole(iDir, iSpin) * Bohr__AA, iDir=1, 3),&
        & iSpin=1, this%nSpin)

    if (mod(iStep, this%writeFreq) == 0) then
      write(qDat, "(2X,2F25.15)", advance="no") time * au__fs, -sum(deltaQ)
      do iAtom = 1, this%nAtom
        write(qDat, "(F25.15)", advance="no")-sum(deltaQ(iAtom,:))
      end do
      write(qDat,*)
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
    integer :: iOrb

    allocate(T2(this%nOrbs, this%nOrbs), T3(this%nOrbs, this%nOrbs))

    T2 = Hsq(:,:,iSpin)
    T3 = 0.0_dp
    do iOrb = 1, this%nOrbs
      T3(iOrb, iOrb) = 1.0_dp
    end do
    call gesv(T2,T3)
    Eiginv(:, :, iSpin) = cmplx(T3, 0, dp)

    T2(:,:) = transpose(Hsq(:,:,iSpin))
    T3 = 0.0_dp
    do iOrb = 1, this%nOrbs
      T3(iOrb, iOrb) = 1.0_dp
    end do
    call gesv(T2,T3)
    EiginvAdj(:, :, iSpin) = cmplx(T3, 0, dp)

    deallocate(T2, T3)

  end subroutine tdPopulInit


  !> Calculate populations at each time step
  subroutine getTDPopulations(rho, Eiginv, EiginvAdj, T1, populDat, time, iSpin)

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

    real(dp) :: occ(size(T1,dim=2))
    integer :: ii

    call gemm(T1, rho(:,:,iSpin), EiginvAdj(:,:,iSpin))
    T1 = transpose(Eiginv(:,:,iSpin)) * T1

    occ = real(sum(T1,dim=1))
    write(populDat(iSpin),'(*(2x,F25.15))', advance='no') time * au__fs
    do ii = 1, size(occ)
      write(populDat(iSpin),'(*(2x,F25.15))', advance='no')occ(ii)
    end do
    write(populDat(iSpin),*)

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

end module timeprop_module
