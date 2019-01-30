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
  use energies, only: TEnergies, init
  use evaluateenergies
  use thirdorder_module, only : ThirdOrder
  use populations
  use eigenvects
  use sk
  use dispiface
  use environment
  use repcont
  use timer
  use taggedoutput
  use hamiltonian
  use solvertypes
  use onsitecorrection
  implicit none
  private

  public :: runDynamics, TElecDynamics_init
  public :: TElecDynamicsInp, TElecDynamics
  public :: iKick, iLaser, iNoTDPert, iKickAndLaser
  public :: iTDConstant, iTDGaussian, iTDSin2, iTDFromFile
  public :: iTDSinglet, iTDTriplet

  !> Data type to  initialize electronic dynamics variables from parser
  type TElecDynamicsInp
    !> External field peak intensity
    real(dp) :: tdfield

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

    integer, allocatable :: indMovedAtom(:), indExcitedAtom(:)
    integer :: nMovedAtom, eulerFreq, tdPPFrames, nExcitedAtom
    logical :: tForces, tIons, tReadMDVelocities, tEulers, tPairWise, tPumpProbe
    logical :: tOnsiteGradients
    real(dp) :: tempAtom
    real(dp), allocatable :: initialVelocities(:,:)
end type TElecDynamicsInp

!> Data type for electronic dynamics internal settings
type TElecDynamics
   private
   real(dp) :: field, dt, omega, time0, time1
   complex(dp) :: fieldDir(3)
   real(dp), allocatable :: tdFunction(:, :), phase
   integer :: nSteps, writeFreq, pertType, envType, spType
   integer :: nAtom, nOrbs, nSpin=1, currPolDir=1, restartFreq
   integer, allocatable :: species(:), polDirs(:), species0(:)
   character(mc), allocatable :: speciesName(:)
   logical :: tPopulations, tSpinPol=.false.
   logical :: tRestart, tWriteRestart, tWriteAutotest
   logical :: tLaser = .false., tKick = .false., tEnvFromFile = .false.
   type(TScc), allocatable :: sccCalc
   character(mc) :: autotestTag

   real(dp), allocatable :: initialVelocities(:,:), movedVelo(:,:), movedMass(:,:)
   real(dp) :: mCutoff, skRepCutoff, rCellVec(3,1)
   real(dp), allocatable :: atomEigVal(:,:), onsiteGrads(:,:,:,:)
   integer :: nExcitedAtom, nMovedAtom, nSparse, eulerFreq, PuProbeFrames
   integer, allocatable :: iCellVec(:), indMovedAtom(:), indExcitedAtom(:)
   logical :: tIons, tForces, tDispersion=.false., ReadMDVelocities, tPumpProbe
   logical :: FirstIonStep = .true., tEulers = .false.
   logical :: tPairWiseEnergy = .false., tCalcOnsiteGradients = .false.
   type(OThermostat), allocatable :: pThermostat
   type(OMDIntegrator), allocatable :: pMDIntegrator
   class(DispersionIface), allocatable :: dispersion
   type(NonSccDiff), allocatable :: derivator
end type TElecDynamics

  !> Enumerating available types of perturbation
  integer, parameter :: iKick = 1, iLaser = 2, iNoTDPert = 3, iKickAndLaser = 4

  !> Enumerating available types of envelope function
  integer, parameter :: iTDConstant = 1, iTDGaussian = 2, iTDSin2 = 3, iTDFromFile = 4

  !> Enumerating available types of spin polarized spectra
  integer, parameter :: iTDSinglet = 1, iTDTriplet = 2

contains

  !> Initialisation of input variables
  subroutine TElecDynamics_init(this, inp, species, speciesName, tWriteAutotest, autotestTag, &
       &randomThermostat, mass, nAtom, skRepCutoff, mCutoff, iCellVec, atomEigVal, dispersion, &
       &nonSccDeriv)

    !> ElecDynamics instance
    type(TElecDynamics), intent(out) :: this

    !> ElecDynamicsInp instance
    type(TElecDynamicsInp), intent(in) :: inp

    !> label for each atomic chemical species
    character(mc), allocatable, intent(in) :: speciesName(:)

    !> produce tagged output?
    logical, intent(in) :: tWriteAutotest

    !> Tagged output files (machine readable)
    character(*), intent(in) :: autotestTag

    real(dp), intent(in), allocatable :: atomEigVal(:,:)
    integer, intent(in) :: nAtom
    integer, allocatable, intent(in) :: iCellVec(:)
    type(ORanlux), allocatable, intent(inout) :: randomThermostat
    real(dp), intent(in) :: mCutoff, skRepCutoff, mass(:)
    class(DispersionIface), allocatable, intent(inout) :: dispersion
    type(NonSccDiff), intent(in) :: nonSccDeriv
    integer, intent(in) :: species(:)

    type(ODummyThermostat), allocatable :: pDummyTherm
    type(OMDCommon), allocatable :: pMDFrame
    real(dp) :: norm, tempAtom
    logical :: tMDstill, tDispersion

    this%field = inp%tdField
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

    this%tIons = inp%tIons
    this%tForces = inp%tForces
    this%tEulers = inp%tEulers
    this%eulerFreq = inp%eulerFreq
    this%tPairWiseEnergy = inp%tPairWise
    this%tPumpProbe = inp%tPumpProbe
    this%PuProbeFrames = inp%tdPPFrames
    this%tCalcOnsiteGradients = inp%tOnsiteGradients

    if (this%tForces .and. .not. this%tIons) then
       this%nMovedAtom = nAtom
       allocate(this%movedMass(3, this%nMovedAtom))
       this%movedMass(:,:) = spread(mass(this%species(:)),1,3)
    end if

    !! For forces
    if (this%tIons) then
       this%tForces = .true.
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
       call init(pDummyTherm, tempAtom, mass(this%indMovedAtom), randomThermostat, pMDFrame)
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
       !call this%dispersion%updateCoords(neighbourList, img2CentCell, coord, species0)
    end if

    this%skRepCutoff = skRepCutoff
    this%mCutoff = mCutoff
    allocate(this%iCellVec, source=iCellVec)
    this%rCellVec(:,1) = [0.0_dp, 0.0_dp, 0.0_dp]
    allocate(this%atomEigVal, source=atomEigVal)
  end subroutine TElecDynamics_init


  !> Driver of time dependent propagation to calculate wither spectrum or laser
  subroutine runDynamics(this, Hsq, ham, H0, species, q0, over, filling, neighbourList,&
      & nNeighbourSK, iSquare, iSparseStart, img2CentCell, orb, coord, spinW, pRepCont, sccCalc,&
      & env, tDualSpinOrbit, xi, thirdOrd, nDftbUFunc, UJ, nUJ, iUJ, niUJ, iHam,&
      & iAtInCentralRegion, tFixEf, Ef, species0, coordAll, onSiteElements, skHamCont, skOverCont)

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

    !> central atomic coordinates
    real(dp), allocatable, intent(inout) :: coord(:,:)

    !> all atomic coordinates
    real(dp), allocatable, intent(inout) :: coordAll(:,:)

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

    type(OSlakoCont), intent(in) :: skHamCont, skOverCont

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

    !> species of central atoms
    integer, intent(in) :: species0(:)

    !> Corrections terms for on-site elements
    real(dp), intent(in), allocatable :: onSiteElements(:,:,:,:)

    integer :: iPol
    logical :: tWriteAutotest


    this%sccCalc = sccCalc
    allocate(this%species(size(species)))
    allocate(this%species0(size(species0)))
    this%species = species
    this%species0 = species0

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
            & iAtInCentralRegion, tFixEf, Ef, tWriteAutotest, coordAll, onSiteElements, skHamCont,&
            & skOverCont)
      end do
    else
      call doDynamics(this, Hsq, ham, H0, species, q0, over, filling, neighbourList, nNeighbourSK,&
          & iSquare, iSparseStart, img2CentCell, orb, coord, spinW, pRepCont, env, tDualSpinOrbit,&
          & xi, thirdOrd, nDftbUFunc, UJ, nUJ, iUJ, niUJ, iHam, iAtInCentralRegion, tFixEf, Ef,&
          & tWriteAutotest, coordAll, onSiteElements, skHamCont, skOverCont)
    end if

  end subroutine runDynamics


  !> Runs the electronic dynamics of the system
  subroutine doDynamics(this, Hsq, ham, H0, species, q0, over, filling, neighbourList,&
      & nNeighbourSK, iSquare, iSparseStart, img2CentCell, orb, coord, spinW, pRepCont, env,&
      & tDualSpinOrbit, xi, thirdOrd, nDftbUFunc, UJ, nUJ, iUJ, niUJ, iHam,&
      & iAtInCentralRegion, tFixEf, Ef, tWriteAutotest, coordAll, onSiteElements, skHamCont,&
      & skOverCont)

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

    !> all atomic coordinates
    real(dp), allocatable, intent(inout) :: coordAll(:,:)

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

    type(OSlakoCont), intent(in) :: skHamCont, skOverCont

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

    !> Should autotest data be written?
    logical, intent(in) :: tWriteAutotest

    !> Corrections terms for on-site elements
    real(dp), intent(in), allocatable :: onSiteElements(:,:,:,:)

    complex(dp) :: Ssqr(this%nOrbs,this%nOrbs), Sinv(this%nOrbs,this%nOrbs)
    complex(dp) :: rho(this%nOrbs,this%nOrbs,this%nSpin), rhoOld(this%nOrbs,this%nOrbs,this%nSpin)
    complex(dp) :: H1(this%nOrbs,this%nOrbs,this%nSpin), T1(this%nOrbs,this%nOrbs)
    complex(dp), allocatable :: Eiginv(:,:,:), EiginvAdj(:,:,:)
    real(dp) :: qq(orb%mOrb, this%nAtom, this%nSpin), deltaQ(this%nAtom,this%nSpin)
    real(dp) :: dipole(3,this%nSpin), chargePerShell(orb%mShell,this%nAtom,this%nSpin)
    real(dp), allocatable :: rhoPrim(:,:), iRhoPrim(:,:), ham0(:), ErhoPrim(:)
    real(dp) :: time, startTime = 0.0_dp, timeElec = 0.0_dp
    integer :: dipoleDat, qDat, energyDat, populDat(2), forceDat, coorDat, ePBondDat
    integer :: iStep = 0, iSpin
    type(TPotentials) :: potential
    type(TEnergies) :: energy
    type(TTimer) :: loopTime
    real(dp), allocatable :: qBlock(:,:,:,:), qiBlock(:,:,:,:)

    complex(dp) :: RdotSprime(this%nOrbs,this%nOrbs)
    real(dp) :: coordNew(3, this%nAtom), totalForce(3, this%nAtom), ePerBond(this%nAtom, this%nAtom)
    real(dp) :: movedAccel(3, this%nMovedAtom), energyKin, new3Coord(3, this%nMovedAtom)
    character(4) :: dumpIdx

    call env%globalTimer%startTimer(globalTimers%elecDynInit)

    if (this%tRestart) then
       call readRestart(rho, rhoOld, Ssqr, coord, this%movedVelo, startTime)
       call updateH0S(this, Ssqr, Sinv, coord, orb, neighbourList, nNeighbourSK, iSquare,&
            & iSparseStart, img2CentCell, skHamCont, skOverCont, ham, ham0, over, env)
       if (this%tIons) then
          this%initialVelocities(:,:) = this%movedVelo
          this%ReadMDVelocities = .true.
       end if
    end if

    if (this%tLaser) then
      call getTDFunction(this, startTime)
    end if

    call initializeTDVariables(this, rho, H1, Ssqr, Sinv, H0, ham0, over, ham, Hsq, filling, orb,&
        & rhoPrim, potential, neighbourList%iNeighbour, nNeighbourSK, iSquare, iSparseStart,&
        & img2CentCell, Eiginv, EiginvAdj, energy, ErhoPrim, skOverCont)

    if (allocated(UJ) .or. allocated(onSiteElements)) then
      allocate(qBlock(orb%mOrb, orb%mOrb, this%nAtom, this%nSpin))
      allocate(qiBlock(orb%mOrb, orb%mOrb, this%nAtom, this%nSpin))
      allocate(iRhoPrim(size(rhoPrim,dim=1), size(rhoPrim,dim=2)))
    end if

    call initTDOutput(this, dipoleDat, qDat, energyDat, populDat, forceDat, coorDat, ePBondDat)

    rhoPrim(:,:) = 0.0_dp
    if (allocated(iRhoPrim)) then
      iRhoPrim(:,:) = 0.0_dp
      do iSpin = 1, this%nSpin
        call packHS(rhoPrim(:,iSpin), iRhoPrim(:,iSpin), rho(:,:,iSpin), neighbourList%iNeighbour,&
            & nNeighbourSK, orb%mOrb, iSquare, iSparseStart, img2CentCell)
      end do
    else
      do iSpin = 1, this%nSpin
        call packHS(rhoPrim(:,iSpin), real(rho(:,:,iSpin), dp), neighbourList%iNeighbour,&
            & nNeighbourSK, orb%mOrb, iSquare, iSparseStart, img2CentCell)
      end do
    end if

    if (allocated(qBlock)) then
      qBlock(:,:,:,:) = 0.0_dp
      qiBlock(:,:,:,:) = 0.0_dp
      do iSpin = 1, this%nSpin
        call mulliken(qBlock(:,:,:,iSpin), over, rhoPrim(:,iSpin), orb, neighbourList%iNeighbour,&
            & nNeighbourSK, img2CentCell, iSparseStart)
        call skewMulliken(qiBlock(:,:,:,iSpin), over, iRhoPrim(:,iSpin), orb,&
            & neighbourList%iNeighbour, nNeighbourSK, img2CentCell, iSparseStart)
      end do
    end if
    call getChargeDipole(this, deltaQ, qq, dipole, q0, rho, Ssqr, coord, iSquare)
    !qq(:,:,:) = 0.0_dp
    !do iSpin = 1, this%nSpin
    !  call mulliken(qq(:,:,iSpin), over, rhoPrim(:,iSpin), orb, neighbourList%iNeighbour,&
    !      & nNeighbourSK, img2CentCell, iSparseStart)
    !end do
    !deltaQ(:,:) = sum((qq - q0), dim=1)
    !dipole(:,:) = -matmul(coord(:,:), deltaQ(:,:))

    call updateH(this, H1, ham, over, H0, species, qq, q0, coord, orb, potential,&
        & neighbourList, nNeighbourSK, iSquare, iSparseStart, img2CentCell, iStep,&
        & chargePerShell, spinW, env, tDualSpinOrbit, xi, thirdOrd, iHam, qBlock, qiBlock,&
        & nDftbUFunc, UJ, nUJ, iUJ, niUJ, onSiteElements)

    if (this%tForces) then
       totalForce = 0.0_dp
       call getForces(this, movedAccel, totalForce, rho, H1, Sinv, neighbourList, nNeighbourSK, &
            & img2CentCell, iSparseStart, iSquare, potential, orb, skHamCont, skOverCont, qq, q0, &
            & pRepCont, coord, rhoPrim, ErhoPrim, 0, env)
       if (this%tIons) then
          call initIonDynamics(this, coordNew, coord, movedAccel)
       end if
    end if

    ! Apply kick to rho if necessary
    if (this%tKick) then
      call kickDM(this, rho, Ssqr, Sinv, iSquare, coord)
    end if

    if (this%tPairWiseEnergy) then
       call pairWiseBondEO(this, ePerBond, rhoPrim(:,1), ham0, iSquare, &
            & neighbourList%iNeighbour, nNeighbourSK, img2CentCell, iSparseStart)
    end if

    ! had to add the "or tKick" option to override rhoOld if tRestart = yes, otherwise it will be
    ! badly initialised
    if (.not.this%tRestart .or. this%tKick) then
      ! Initialize electron dynamics
      rhoOld(:,:,:) = rho
      call initializePropagator(this, -this%dt, rhoOld, rho, H1, Sinv, coord, skHamCont,&
         & skOverCont, orb, neighbourList, nNeighbourSK, img2CentCell, iSquare)
    end if

    call getTDEnergy(this, energy, rhoPrim, rho, neighbourList, nNeighbourSK, orb,&
         & iSquare, iSparseStart, img2CentCell, ham0, qq, q0, potential, chargePerShell, coord, &
         & pRepCont, energyKin, tDualSpinOrbit, thirdOrd, qBlock, qiBlock, nDftbUFunc, UJ, nUJ,&
         & iUJ, niUJ, xi, iAtInCentralRegion, tFixEf, Ef, onSiteElements)

    call env%globalTimer%stopTimer(globalTimers%elecDynInit)

    ! Main loop
    call env%globalTimer%startTimer(globalTimers%elecDynLoop)
    call loopTime%start()

    write(stdOut, "(A)") 'Starting dynamics'

    do iStep = 0, this%nSteps
      time = iStep * this%dt + startTime

      if (iStep > 0) then
         call writeTDOutputs(this, dipoleDat, qDat, energyDat, forceDat, coorDat, time-this%dt,&
              & energy, energyKin, dipole, deltaQ, coord, totalForce, iStep-1, ePerBond, ePBondDat)
      end if

     if (this%tIons) then
        coord(:,:) = coordNew
        call updateH0S(this, Ssqr, Sinv, coord, orb, neighbourList, nNeighbourSK, iSquare,&
             &iSparseStart, img2CentCell, skHamCont, skOverCont, ham, ham0, over, env)
     end if

     if ((this%tPumpProbe) .and. (mod(iStep, this%nSteps/this%PuProbeFrames) == 0)) then
        write(dumpIdx,'(i4)')int(this%PuProbeFrames*iStep/this%nSteps)
        call writeRestart(rho, rhoOld, Ssqr, coord, this%movedVelo, 0.0_dp,&
             & trim(dumpIdx) // 'ppdump.bin')
     end if

      rhoPrim(:,:) = 0.0_dp
      if (allocated(iRhoPrim)) then
        iRhoPrim(:,:) = 0.0_dp
        do iSpin = 1, this%nSpin
          call packHS(rhoPrim(:,iSpin), iRhoPrim(:,iSpin), rho(:,:,iSpin),&
              & neighbourList%iNeighbour, nNeighbourSK, orb%mOrb, iSquare, iSparseStart,&
              & img2CentCell)
        end do
      else
        do iSpin = 1, this%nSpin
          call packHS(rhoPrim(:,iSpin), real(rho(:,:,iSpin), dp), neighbourList%iNeighbour,&
              & nNeighbourSK, orb%mOrb, iSquare, iSparseStart, img2CentCell)
        end do
      end if

      if (allocated(qBlock)) then
        qBlock(:,:,:,:) = 0.0_dp
        qiBlock(:,:,:,:) = 0.0_dp
        do iSpin = 1, this%nSpin
          call mulliken(qBlock(:,:,:,iSpin), over, rhoPrim(:,iSpin), orb, neighbourList%iNeighbour,&
              & nNeighbourSK, img2CentCell, iSparseStart)
          call skewMulliken(qiBlock(:,:,:,iSpin), over, iRhoPrim(:,iSpin), orb,&
              & neighbourList%iNeighbour, nNeighbourSK, img2CentCell, iSparseStart)
        end do
      end if
      call getChargeDipole(this, deltaQ, qq, dipole, q0, rho, Ssqr, coord, iSquare)
      !qq(:,:,:) = 0.0_dp
      !do iSpin = 1, this%nSpin
      !  call mulliken(qq(:,:,iSpin), over, rhoPrim(:,iSpin), orb, neighbourList%iNeighbour,&
      !      & nNeighbourSK, img2CentCell, iSparseStart)
      !end do
      !deltaQ(:,:) = sum((qq - q0), dim=1)
      !dipole(:,:) = -matmul(coord(:,:), deltaQ(:,:))

     call updateH(this, H1, ham, over, H0, species, qq, q0, coord, orb, potential,&
          & neighbourList, nNeighbourSK, iSquare, iSparseStart, img2CentCell, iStep,&
          & chargePerShell, spinW, env, tDualSpinOrbit, xi, thirdOrd, iHam, qBlock, qiBlock,&
          & nDftbUFunc, UJ, nUJ, iUJ, niUJ, onSiteElements)

     if ((this%tWriteRestart) .and. (iStep > 0) .and. (mod(iStep, this%restartFreq) == 0)) then
        call writeRestart(rho, rhoOld, Ssqr, coord, this%movedVelo, time)
     end if

     if (this%tForces) then
        call getForces(this, movedAccel, totalForce, rho, H1, Sinv, neighbourList,&  !F_1
             & nNeighbourSK, img2CentCell, iSparseStart, iSquare, potential, orb, skHamCont, &
             & skOverCont, qq, q0, pRepCont, coord, rhoPrim, ErhoPrim, iStep, env)
     end if

     if (this%tIons) then
        new3Coord(:,:) = coordNew(:, this%indMovedAtom)
        call next(this%pMDIntegrator, movedAccel, new3Coord, this%movedVelo)
        coordNew(:, this%indMovedAtom) = new3Coord
        call getRdotSprime(this, RdotSprime, coord, skOverCont, orb, img2CentCell, &
             &neighbourList, nNeighbourSK, iSquare)
     end if

     call getTDEnergy(this, energy, rhoPrim, rho, neighbourList, nNeighbourSK, orb,&
          & iSquare, iSparseStart, img2CentCell, ham0, qq, q0, potential, chargePerShell, coord, &
          & pRepCont, energyKin, tDualSpinOrbit, thirdOrd, qBlock, qiBlock, nDftbUFunc, UJ, nUJ,&
          & iUJ, niUJ, xi, iAtInCentralRegion, tFixEf, Ef, onSiteElements)

     if (this%tPairWiseEnergy) then
        call pairWiseBondEO(this, ePerBond, rhoPrim(:,1), ham0, iSquare, &
             & neighbourList%iNeighbour, nNeighbourSK, img2CentCell, iSparseStart)
     end if

     do iSpin = 1, this%nSpin
        if (this%tIons) then
           call scal(H1(:,:,iSpin), imag)
           call zaxpy(this%nOrbs*this%nOrbs, 1.0_dp, RdotSprime, 1, H1(:,:,iSpin), 1)

           if (this%tEulers .and. (iStep > 100) .and. (mod(iStep, this%eulerFreq) == 0)) then
              call zcopy(this%nOrbs*this%nOrbs, rho(:,:,iSpin), 1, rhoOld(:,:,iSpin), 1)
              call propagateRho(this, rhoOld(:,:,iSpin), rho(:,:,iSpin), H1(:,:,iSpin), Sinv,&
                   & this%dt)
           else
              call propagateRho(this, rhoOld(:,:,iSpin), rho(:,:,iSpin), H1(:,:,iSpin), Sinv,&
                   & 2.0_dp * this%dt)
           end if
        else
           !The following is commented for the fast popagate that considers a real H
           !call scal(H1(:,:,iSpin), imag)
           call propagateRhoRealH(this, rhoOld(:,:,iSpin), rho(:,:,iSpin), H1(:,:,iSpin), Sinv,&
                & 2.0_dp * this%dt)
        end if

        call swap(rhoOld(:,:,iSpin), rho(:,:,iSpin))

        if ((this%tPopulations) .and. (mod(iStep, this%writeFreq) == 0)) then
           call getTDPopulations(this, rho, Eiginv, EiginvAdj, T1, H1, Ssqr, populDat, time, iSpin)
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

    if (tWriteAutotest) then
      call writeTDAutotest(this, dipole, energy, deltaQ)
    end if

    call closeTDOutputs(this, dipoleDat, qDat, energyDat, populDat, forceDat, coorDat, ePBondDat)

    if (this%tIons) then
       deallocate(this%pMDIntegrator)
    end if

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
    real(dp), intent(inout), allocatable :: qBlock(:,:,:,:)

    !> Imaginary part of block atomic populations
    real(dp), intent(inout), allocatable :: qiBlock(:,:,:,:)

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

    ! Build spin contribution (if necessary)
!    if (this%tSpinPol) then
!       call getSpinShift(shellPot, chargePerShell, this%species, orb, W)
!       potential%intShell = potential%intShell + shellPot
    if (allocated(UJ) .or. allocated(onSiteElements)) then
      ! convert to qm representation
      call ud2qm(qBlock)
      call ud2qm(qiBlock)
    end if

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

    potential%intBlock = potential%intBlock + potential%extBlock ! is this necessary?
    potential%intShell = potential%intShell + potential%extShell ! for SCC

    call getSccHamiltonian(H0, over, nNeighbourSK, neighbourList, species, orb, iSparseStart,&
        & img2CentCell, potential, ham, iHam)

    ! Hack due to not using Pauli-type structure outside of this part of the routine
    if (this%nSpin == 2) then
      ham = 2.0_dp * ham
      call qm2ud(ham)
      call qm2ud(q0)
      call qm2ud(qq)
    end if

    do iSpin = 1, this%nSpin
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

    pkick(1) = this%field ! check units

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
    integer :: iStep, laserDat, laserDat2

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


  !> Calculate energy - modify to include new way to calculate energy
  subroutine getTDEnergy(this, energy, rhoPrim, rho, neighbourList, nNeighbourSK, orb,&
       & iSquare, iSparseStart, img2CentCell, ham0, qq, q0, potential, chargePerShell, coord, &
       & pRepCont, energyKin, tDualSpinOrbit, thirdOrd, qBlock, qiBlock, nDftbUFunc, UJ, nUJ, iUJ,&
       & niUJ, xi, iAtInCentralRegion, tFixEf, Ef, onSiteElements)

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

    !> neighbour list
    type(TNeighbourList), intent(in) :: neighbourList

    !> Number of neighbours for each of the atoms
    integer, intent(in) :: nNeighbourSK(:)

    !> Index array for start of atomic block in dense matrices
    integer, intent(in) :: iSquare(:)

    !> index array for location of atomic blocks in large sparse arrays
    integer, intent(in) :: iSparseStart(0:,:)

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

    !> kinetic energy
    real(dp), intent(out) :: energyKin

    !> Is dual spin orbit being used
    logical, intent(in) :: tDualSpinOrbit

        !> 3rd order settings
    type(ThirdOrder), intent(inout), allocatable :: thirdOrd

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

    !> Spin orbit constants
    real(dp), intent(in), allocatable :: xi(:,:)

    !> Atoms over which to sum the total energies
    integer, intent(in) :: iAtInCentralRegion(:)

    !> Whether fixed Fermi level(s) should be used. (No charge conservation!)
    logical, intent(in) :: tFixEf

    !> If tFixEf is .true. contains reservoir chemical potential, otherwise the Fermi levels found
    !> from the given number of electrons
    real(dp), intent(inout) :: Ef(:)

    !> Corrections terms for on-site elements
    real(dp), intent(in), allocatable :: onSiteElements(:,:,:,:)

    integer :: iSpin
    logical :: tUseBuggyRepSum = .false.
    real(dp) :: TS(this%nSpin)

    if (size(rhoPrim, dim=1) /= this%nSparse) then
       deallocate(rhoPrim)
       allocate(rhoPrim(this%nSparse, this%nSpin))
    end if

    rhoPrim(:,:) = 0.0_dp
    do iSpin = 1, this%nSpin
      call packHS(rhoPrim(:,iSpin), real(rho(:,:,iSpin), dp), neighbourList%iNeighbour,&
          & nNeighbourSK, orb%mOrb, iSquare, iSparseStart, img2CentCell)
    end do

    call ud2qm(rhoPrim)

    ! Calculate repulsive energy
    call getERep(energy%atomRep, coord, nNeighbourSK, neighbourList%iNeighbour, this%species,&
         & pRepCont, img2CentCell)
    energy%Erep = sum(energy%atomRep)

    ! Calculate dispersion component
    if (this%tDispersion) then
       call this%dispersion%getEnergies(energy%atomDisp)
       energy%eDisp = sum(energy%atomDisp)
    else
       energy%atomDisp(:) = 0.0_dp
       energy%eDisp = 0.0_dp
    end if

    TS = 0.0_dp
    call getEnergies(this%sccCalc, qq, q0, chargePerShell, this%species, this%tLaser, .false.,&
         & .false., tDualSpinOrbit, rhoPrim, ham0, orb, neighbourList, nNeighbourSK, img2CentCell,&
         & iSparseStart, 0.0_dp, 0.0_dp, TS, potential, energy, thirdOrd, qBlock, qiBlock,&
         & nDftbUFunc, UJ, nUJ, iUJ, niUJ, xi, iAtInCentralRegion, tFixEf, Ef, onSiteElements)
    ! getEnergies returns the total energy Etotal including repulsive and dispersions energies

    ! Calculate nuclear kinetic energy
    energyKin = 0.0_dp
    if (this%tIons) then
       energyKin = 0.5_dp * sum(this%movedMass(:,:) * this%movedVelo(:,:)**2)
       energy%Etotal = energy%Etotal + energyKin
    end if

  end subroutine getTDEnergy


  !> Create all necessary matrices and instances for dynamics
  subroutine initializeTDVariables(this, rho, H1, Ssqr, Sinv, H0, ham0, over, ham, Hsq, filling,&
      & orb, rhoPrim, potential, iNeighbour, nNeighbourSK, iSquare, iSparseStart, img2CentCell,&
      & Eiginv, EiginvAdj, energy, ErhoPrim, skOverCont)

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

    !> Non-SCC hamitonian
    real(dp), intent(in) :: H0(:)

    !> Local sparse storage for non-SCC hamitonian
    real(dp), allocatable, intent(out) :: ham0(:)

    type(OSlakoCont), intent(in) :: skOverCont
    real(dp), allocatable, intent(out) :: ErhoPrim(:)

    real(dp) :: T2(this%nOrbs,this%nOrbs), T3(this%nOrbs, this%nOrbs)
    integer :: iSpin, iOrb, iOrb2

    allocate(rhoPrim(size(ham, dim=1), this%nSpin), ErhoPrim(size(ham, dim=1)))
    this%nSparse = size(H0)
    allocate(ham0(size(H0)))
    ham0(:) = H0

    T2 = 0.0_dp
    T3 = 0.0_dp
    call unpackHS(T2,over,iNeighbour,nNeighbourSK,iSquare,iSparseStart,img2CentCell)
    call blockSymmetrizeHS(T2,iSquare)
    Ssqr(:,:) = cmplx(T2, 0, dp)

    do iSpin = 1, this%nSpin
      call unpackHS(T3, ham(:,iSpin), iNeighbour, nNeighbourSK, iSquare, iSparseStart, img2CentCell)
      call blockSymmetrizeHS(T3, iSquare)
      H1(:,:,iSpin) = cmplx(T3, 0, dp)
      T3 = 0.0_dp
    end do

    if (this%tPopulations) then
      allocate(Eiginv(this%nOrbs, this%nOrbs, this%nSpin))
      allocate(EiginvAdj(this%nOrbs, this%nOrbs, this%nSpin))
      do iSpin=1, this%nSpin
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

    if (this%tIons .and. this%tCalcOnsiteGradients) then
       allocate(this%onsiteGrads(3, this%nAtom, orb%mOrb, orb%mOrb))
       call getOnsiteGrads(this, skOverCont, orb)
    end if

  end subroutine initializeTDVariables


  !> Perfoms a step backwards to boot the dynamics using the Euler algorithm
  subroutine initializePropagator(this, step, rhoOld, rho, H1, Sinv, coord, skHamCont,&
         & skOverCont, orb, neighbourList, nNeighbourSK, img2CentCell, iSquare)

    !> ElecDynamics instance
    type(TElecDynamics), intent(inout) :: this

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

    real(dp), intent(in) :: coord(:,:)
    type(OSlakoCont), intent(in) :: skHamCont, skOverCont
    type(TOrbitals), intent(in) :: orb
    type(TNeighbourList), intent(in) :: neighbourList
    integer, intent(in) :: nNeighbourSK(:), iSquare(:)
    integer, allocatable, intent(in) :: img2CentCell(:)
    integer :: iSpin
    complex(dp) :: RdotSprime(this%nOrbs,this%nOrbs)

    rhoOld(:,:,:) = rho

    if (this%tIons) then
       call getRdotSprime(this, RdotSprime, coord, skOverCont, orb, img2CentCell, &
            &neighbourList, nNeighbourSK, iSquare)
    end if

    do iSpin=1,this%nSpin
       if (this%tIons) then
          H1(:,:,iSpin) = RdotSprime + imag * H1(:,:,iSpin)
          call propagateRho(this, rhoOld(:,:,iSpin), rho(:,:,iSpin), H1(:,:,iSpin), Sinv, step)
       else
          ! The following line is commented to make the fast propagate work since it needs a real H
          !H1(:,:,iSpin) = imag * H1(:,:,iSpin)
          call propagateRhoRealH(this, rhoOld(:,:,iSpin), rho(:,:,iSpin), H1(:,:,iSpin), Sinv, step)
       end if
    end do

  end subroutine initializePropagator


  !> Propagate rho, notice that H = iH (coeficients are real)
  subroutine propagateRho(this, rhoOld, rho, H1, Sinv, step)

    !> ElecDynamics instance
    type(TElecDynamics), intent(inout) :: this

    !> Density matrix at previous step
    complex(dp), intent(inout) :: rhoOld(:,:)

    !> Density matrix
    complex(dp), intent(in) :: rho(:,:)

    !> Square hamiltonian
    complex(dp), intent(in) :: H1(:,:)

    !> Square overlap inverse
    complex(dp), intent(in) :: Sinv(:,:)

    !> Time step in atomic units
    real(dp), intent(in) :: step

    complex(dp) :: T1(this%nOrbs,this%nOrbs)

    T1 = 0.0_dp
    call gemm(T1, Sinv, H1)
    call gemm(rhoOld, T1, rho, cmplx(-step, 0, dp), cmplx(1, 0, dp))
    call gemm(rhoOld, rho, T1, cmplx(-step, 0, dp), cmplx(1, 0, dp), 'N', 'C')

  end subroutine propagateRho


  !> Propagate rho for real Hamiltonian (used for frozen nuclei dynamics and gamma point periodic)
  subroutine propagateRhoRealH(this, rhoOld, rho, H1, Sinv, step)

    !> ElecDynamics instance
    type(TElecDynamics), intent(inout) :: this

    !> Density matrix at previous step
    complex(dp), intent(inout) :: rhoOld(:,:)

    !> Density matrix
    complex(dp), intent(in) :: rho(:,:)

    !> Square hamiltonian
    complex(dp), intent(in) :: H1(:,:)

    !> Square overlap inverse
    complex(dp), intent(in) :: Sinv(:,:)

    !> Time step in atomic units
    real(dp), intent(in) :: step

    real(dp) :: T1R(this%nOrbs,this%nOrbs),T2R(this%nOrbs,this%nOrbs)
    real(dp) :: T3R(this%nOrbs,this%nOrbs),T4R(this%nOrbs,this%nOrbs),T5R(this%nOrbs,this%nOrbs)
    integer :: i,j

    ! The code below takes into account that Sinv and H1 are real, this is twice as fast as the
    ! original above (propageteRho)

    ! get the real part of Sinv and H1
    T1R = real(H1)
    T2R = real(Sinv)
    call gemm(T3R,T2R,T1R)

    ! calculate the first term products for the real and imaginary parts independently
    T1R = real(rho)
    T2R = aimag(rho)
    call gemm(T4R,T3R,T1R)
    call gemm(T5R,T3R,T2R)

    ! build the commutator combining the real and imaginary parts of the previous result
    !$omp parallel do private(i,j)
    do i=1,this%nOrbs
      do j=1,this%nOrbs
        rhoOld(i,j) = rhoOld(i,j) + cmplx(0, -step, dp) * (T4R(i,j) + imag * T5R(i,j)) &
                      + cmplx(0, step, dp) * conjg(T4R(j,i) + imag * T5R(j,i))
      enddo
    enddo
    !$omp end parallel do

  end subroutine propagateRhoRealH


  !> Initialize output files
  subroutine initTDOutput(this, dipoleDat, qDat, energyDat, populDat, forceDat, coorDat, ePBondDat)
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
    logical :: exist=.false.

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

    if (this%tForces) then
       call openFile(this, forceDat, 'forcesvst.dat')
    end if

    if (this%tIons) then
       call openFile(this, coorDat, 'tdcoords.xyz')
    end if

    if (this%tPairWiseEnergy) then
       if (this%tRestart) then
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
  subroutine closeTDOutputs(this, dipoleDat, qDat, energyDat, populDat, forceDat, coorDat,&
       & ePBondDat)

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

    if (this%tIons) then
      close(coorDat)
    end if

    if (this%tForces) then
      close(forceDat)
   end if

   if (this%tPairWiseEnergy) then
      close(ePBondDat)
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
  subroutine writeRestart(rho, rhoOld, Ssqr, coord, veloc, time, dumpName)

    complex(dp), intent(in) :: rho(:,:,:), rhoOld(:,:,:), Ssqr(:,:)

    !> atomic coordinates
    real(dp), intent(in) :: coord(:,:)

    !> elapsed simulated time in atomic units
    real(dp), intent(in) :: time

    !> name of the dump file
    character(len=*), intent(in), optional :: dumpName

    !> atomic velocities
    real(dp), intent(in) :: veloc(:,:)

    integer :: dumpBin

    if (present(dumpName)) then
       open(newunit=dumpBin, file=dumpName, form='unformatted', access='stream', action='write')
    else
       open(newunit=dumpBin, file='tddump.bin', form='unformatted', access='stream', action='write')
    end if

    write(dumpBin) rho, rhoOld, Ssqr, coord, veloc, time
    close(dumpBin)
  end subroutine writeRestart


  !> read a restart file containing density matrix, overlap, coordinates and time step
  subroutine readRestart(rho, rhoOld, Ssqr, coord, veloc, time)

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

    !> atomic velocities
    real(dp), intent(out) :: veloc(:,:)
    integer :: dumpBin

    open(newunit=dumpBin, file='tddump.bin', form='unformatted', access='stream', action='read')
    read(dumpBin) rho, rhoOld, Ssqr, coord, veloc, time
    close(dumpBin)
  end subroutine readRestart

  !> Write results to file
  subroutine writeTDOutputs(this, dipoleDat, qDat, energyDat, forceDat, coorDat, &
       & time, energy, energyKin, dipole, deltaQ, coord, totalForce, iStep, ePerBond, ePBondDat)

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

    integer :: iAtom, iAtom2, iSpin, iDir

    integer, intent(in) :: forceDat, coorDat, ePBondDat
    real(dp), intent(in) :: coord(:,:), ePerBond(:,:)
    real(dp), intent(in) :: energyKin, totalForce(:,:)
    real(dp) :: auxVeloc(3, this%nAtom)

    write(energydat, '(9F25.15)') time * au__fs, energy%Etotal, energy%EnonSCC, energy%eSCC,&
         & energy%Espin, energy%Eext, energy%Erep, energyKin, energy%eDisp
    write(dipoleDat, '(7F25.15)') time * au__fs, ((dipole(iDir, iSpin) * Bohr__AA, iDir=1, 3),&
         & iSpin=1, this%nSpin)

    if (mod(iStep, this%writeFreq) == 0) then
      write(qDat, "(2X,2F25.15)", advance="no") time * au__fs, -sum(deltaQ)
      do iAtom = 1, this%nAtom
        write(qDat, "(F25.15)", advance="no")-sum(deltaQ(iAtom,:))
      end do
      write(qDat,*)
    end if

    if (this%tIons .and. (mod(iStep,this%writeFreq) == 0)) then
       auxVeloc = 0.0_dp
       auxVeloc(:, this%indMovedAtom) = this%movedVelo
       write(coorDat,'(I5)')this%nAtom
       write(coorDat,*) 'MD step:', iStep, 'time', time * au__fs
       do iAtom=1,this%nAtom
          write(coorDat, '(A2, 6F16.8)') trim(this%speciesName(this%species(iAtom))), &
               &coord(:, iAtom) * Bohr__AA, auxVeloc(:, iAtom) * Bohr__AA / au__fs
       end do
    endif

    if (this%tForces .and. (mod(iStep,this%writeFreq) == 0)) then
      write(forceDat, "(F25.15)", advance="no") time * au__fs
      do iAtom = 1, this%nAtom
        write(forceDat, "(3F25.15)", advance="no") totalForce(:,iAtom)
      end do
      write(forceDat,*)
    end if

!    if (this%tForces .and. (mod(iStep,this%writeFreq) == 0)) then
!       write(forceDat, '(10000F25.15)') time * au__fs, (totalForce(:,iAtom), iAtom=1,this%nAtom)
!    end if

    if (this%tPairWiseEnergy .and. mod(iStep,this%writeFreq) == 0) then
       write(ePBondDat) time * au__fs, sum(ePerBond(:,:)), &
            & ((ePerBond(iAtom, iAtom2), iAtom=1,this%nAtom), iAtom2=1,this%nAtom)
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
  subroutine getTDPopulations(this, rho, Eiginv, EiginvAdj, T1, H1, Ssqr, populDat, time, iSpin)

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

    real(dp) :: occ(size(T1,dim=2))
    integer :: ii

    complex(dp), intent(inout) :: H1(:,:,:), Ssqr(:,:)
    real(dp) :: eigen(this%nOrbs)

    if (this%tIons) then
       H1(:,:,iSpin) = -imag * H1(:,:,iSpin) ! change back to real H1
       ! only do previous step because of using other propagateRho for frozen nuclei
       call diagDenseMtx(2, 'V', H1(:,:,iSpin), SSqr, eigen)
       call tdPopulInit(this, Eiginv, EiginvAdj, real(H1), iSpin)
    end if

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


  !> Initialize ion dynamics
  subroutine initIonDynamics(this, coordNew, coord, movedAccel)
    type(TElecDynamics), intent(inout) :: this
    type(OVelocityVerlet), allocatable :: pVelocityVerlet
    real(dp), intent(in) :: coord(:,:), movedAccel(:,:)
    real(dp), intent(out) :: coordNew(:,:)
    logical :: halfVelocities = .true.
    real(dp) :: velocities(3, this%nMovedAtom)

    allocate(pVelocityVerlet)
    if (this%ReadMDVelocities) then
       call init(pVelocityVerlet, this%dt, coord(:, this%indMovedAtom),&
            & this%pThermostat, this%initialVelocities, halfVelocities)
       this%movedVelo(:, :) = this%initialVelocities
    else
       call init(pVelocityVerlet, this%dt, coord(:, this%indMovedAtom), &
            &this%pThermostat, halfVelocities, velocities)
       this%movedVelo(:,:) = velocities
    end if

    ! Euler step forward
    this%movedVelo(:,:) = this%movedVelo - 0.5_dp * movedAccel * this%dt ! Has ensures good initialization
    coordNew(:,:) = coord
    coordNew(:,this%indMovedAtom) = coord(:,this%indMovedAtom) &
         & + this%movedVelo(:,:) * this%dt + 0.5_dp * movedAccel(:,:) * this%dt**2

    ! This re-initializes the VVerlet propagator with coordNew
    this%movedVelo(:,:) = this%movedVelo + 0.5_dp * movedAccel * this%dt
    call init(pVelocityVerlet, this%dt, coordNew(:, this%indMovedAtom),&
         & this%pThermostat, this%movedVelo, halfVelocities)
    allocate(this%pMDIntegrator)
    call init(this%pMDIntegrator, pVelocityVerlet)
  end subroutine initIonDynamics


  !> Calculates forces on nuclei and updates positions
  subroutine updateH0S(this, Ssqr, Sinv, coord, orb, neighbourList, nNeighbourSK, iSquare, iSparseStart,&
       & img2CentCell, skHamCont, skOverCont, ham, ham0, over, env)
    type(TElecDynamics), intent(inout), target :: this
    complex(dp), intent(inout) :: Sinv(:,:), Ssqr(:,:)
    real(dp), allocatable, intent(inout) :: ham0(:)
    real(dp), allocatable, intent(inout) :: ham(:,:), over(:)
    real(dp), allocatable, intent(inout) :: coord(:,:)

    type(TNeighbourList), intent(inout) :: neighbourList
    integer, intent(inout) :: nNeighbourSK(:)
    integer, allocatable, intent(inout) :: iSparseStart(:,:), img2CentCell(:)
    integer, intent(in) :: iSquare(:)
    type(TOrbitals), intent(in) :: orb
    type(OSlakoCont), intent(in) :: skHamCont, skOverCont
    !> Environment settings
    type(TEnvironment), intent(inout) :: env

    real(dp) :: Sreal(this%nOrbs,this%nOrbs), SinvReal(this%nOrbs,this%nOrbs)
    real(dp) :: coord0Fold(3,this%nAtom)
    integer :: nAllAtom, nAllOrb
    integer, allocatable :: specie0(:)
    integer :: info, ipiv(this%nOrbs, this%nOrbs)
    integer :: iSpin, sparseSize, iOrb

    !! Calculate overlap, H0 for new geometry
    nAllAtom = this%nAtom
    coord0Fold(:,:) = coord
    allocate(specie0, source=this%species)

    call updateNeighbourListAndSpecies(coord, this%species, img2CentCell, this%iCellVec, &
         &neighbourList, nAllAtom, coord0Fold, specie0, this%mCutoff, this%rCellVec)
    nAllOrb = sum(orb%nOrbSpecies(this%species(1:nAllAtom)))
    call getNrOfNeighboursForAll(nNeighbourSK, neighbourList, this%skRepCutoff)
    call getSparseDescriptor(neighbourList%iNeighbour, nNeighbourSK, img2CentCell, orb, iSparseStart,&
         & sparseSize)

    deallocate(ham)
    deallocate(over)
    deallocate(ham0)
    allocate(ham(sparseSize, this%nSpin))
    allocate(over(sparseSize))
    allocate(ham0(sparseSize))
    this%nSparse = sparseSize

    call this%sccCalc%updateCoords(env, coord, this%species, neighbourList)

   if (this%tDispersion) then
      call this%dispersion%updateCoords(neighbourList, img2CentCell, coord, specie0)
   end if

   call buildH0(env, ham0, skHamCont, this%atomEigVal, coord, nNeighbourSK, &
        & neighbourList%iNeighbour, this%species, iSparseStart, orb)
   call buildS(env, over, skOverCont, coord, nNeighbourSK, neighbourList%iNeighbour, this%species,&
        & iSparseStart, orb)

   Sreal = 0.0_dp
   call unpackHS(Sreal,over,neighbourList%iNeighbour,nNeighbourSK,iSquare,iSparseStart,img2CentCell)
   call blockSymmetrizeHS(Sreal,iSquare)
   Ssqr(:,:) = cmplx(Sreal, 0, dp)

   SinvReal = 0.0_dp
   do iOrb = 1, this%nOrbs
      SinvReal(iOrb, iOrb) = 1.0_dp
   end do
   call gesv(Sreal,SinvReal)
   Sinv(:,:) = cmplx(SinvReal, 0, dp)
  end subroutine updateH0S


  !> Calculates force
  subroutine getForces(this, movedAccel, totalForce, rho, H1, Sinv, neighbourList, nNeighbourSK,&
       & img2CentCell, iSparseStart, iSquare, potential, orb, skHamCont, skOverCont, qq, q0,&
       & pRepCont, coord, rhoPrim, ErhoPrim, iStep, env)
    type(TElecDynamics), intent(inout), target :: this
    complex(dp), intent(in) :: rho(:,:,:), H1(:,:,:), Sinv(:,:)
    type(TNeighbourList), intent(in) :: neighbourList
    integer, intent(in) :: nNeighbourSK(:)
    integer, allocatable, intent(in) :: iSparseStart(:,:), img2CentCell(:)
    integer, intent(in) :: iSquare(:)
    real(dp), intent(out) :: totalForce(3, this%nAtom)
    real(dp), intent(out) :: movedAccel(3, this%nMovedAtom)
    type(TPotentials), intent(in) :: potential
    type(TOrbitals), intent(in) :: orb
    type(OSlakoCont), intent(in) :: skHamCont, skOverCont
    real(dp), intent(inout) :: qq(:,:,:), q0(:,:,:)
    type(ORepCont), intent(in) :: pRepCont
    real(dp), intent(in) :: coord(:,:)
    real(dp), allocatable, intent(inout) :: rhoPrim(:,:), ErhoPrim(:)
    integer, intent(in) :: iStep
    real(dp) :: T1(this%nOrbs,this%nOrbs),T2(this%nOrbs,this%nOrbs)
    real(dp) :: derivs(3,this%nAtom), repulsiveDerivs(3,this%nAtom), totalDeriv(3, this%nAtom)
    integer :: iSpin, iDir
    !> Environment settings
    type(TEnvironment), intent(inout) :: env

    if (size(rhoPrim, dim=1) /= this%nSparse) then
       deallocate(rhoPrim, ErhoPrim)
       allocate(rhoPrim(this%nSparse, this%nSpin))
       allocate(ErhoPrim(this%nSparse))
    end if

    ErhoPrim(:) = 0.0_dp
    rhoPrim(:,:) = 0.0_dp

    do iSpin = 1, this%nSpin
       call gemm(T1, real(H1(:,:,iSpin), dp), real(rho(:,:,iSpin), dp))
       call gemm(T2, real(Sinv, dp), T1, 0.5_dp, 0.0_dp)
       call daxpy(this%nOrbs*this%nOrbs, 1.0_dp, transpose(T2), 1, T2, 1)
       !    T1 = T2 + transpose(T2)

       call packHS(rhoPrim(:,iSpin), real(rho(:,:,iSpin), dp), neighbourList%iNeighbour, &
            &nNeighbourSK, orb%mOrb, iSquare, iSparseStart, img2CentCell)
       call packHS(ErhoPrim, T2, neighbourList%iNeighbour, nNeighbourSK, orb%mOrb, iSquare, iSparseStart,&
            & img2CentCell)
    end do

    if (this%nSpin == 2) then
       call ud2qm(qq)
       call ud2qm(q0)
       call ud2qm(rhoPrim)
    end if

    derivs(:,:) = 0.0_dp
    repulsiveDerivs(:,:) = 0.0_dp

    call derivative_shift(env, derivs, this%derivator, rhoPrim, ErhoPrim(:), skHamCont,&
         &skOverCont, coord, this%species, neighbourList%iNeighbour, nNeighbourSK, img2CentCell,&
         & iSparseStart, orb, potential%intBlock)
    call this%sccCalc%updateCharges(env, qq, q0, orb, this%species)
    call this%sccCalc%addForceDc(env, derivs, this%species, neighbourList%iNeighbour, &
         & img2CentCell)
    call getERepDeriv(repulsiveDerivs, coord, nNeighbourSK, neighbourList%iNeighbour, this%species,&
         & pRepCont, img2CentCell)

    if (this%tLaser) then
       do iDir = 1, 3
          derivs(iDir, :)=derivs(iDir,:)-sum(q0(:,:,1)-qq(:,:,1), dim=1)*this%TDFunction(iStep,iDir)
       end do
    end if

    totalDeriv(:,:) = repulsiveDerivs(:,:) + derivs(:,:)
    if (this%tDispersion) then
       call this%dispersion%addGradients(totalDeriv)
    end if

    totalForce(:,:) = - totalDeriv(:,:)

    if (this%tIons) then
       movedAccel(:,:) = totalForce(:, this%indMovedAtom) / this%movedMass
    else
       movedAccel(:,:) = totalForce / this%movedMass
    end if

    if (this%nSpin == 2) then
       call qm2ud(qq)
       call qm2ud(q0)
       call qm2ud(rhoPrim)
    end if

  end subroutine getForces


  !> Calculates nonadiabatic matrix (Sprime) times velocities (Rdot)
  subroutine getRdotSprime(this, RdotSprime, coord, skOverCont, orb, img2CentCell, &
       &neighbourList, nNeighbourSK, iSquare)
    type(TElecDynamics), intent(in), target :: this
    type(OSlakoCont), intent(in) :: skOverCont
    complex(dp), intent(out) :: RdotSprime(this%nOrbs,this%nOrbs)
    type(TOrbitals), intent(in) :: orb
    real(dp) :: sPrimeTmp(orb%mOrb,orb%mOrb,3)
    real(dp) :: sPrimeTmp2(orb%mOrb,orb%mOrb), dcoord(3,this%nAtom)
    real(dp), intent(in) :: coord(:,:)
    integer :: iAtom1,iStart1,iEnd1,iSp1,nOrb1,iAtomAux,iDir
    integer :: iNeigh,iStart2,iEnd2,iAtom2,iAtom2f,iSp2,nOrb2
    type(TNeighbourList), intent(in) :: neighbourList
    integer, intent(in) :: nNeighbourSK(:)
    integer, intent(in) :: iSquare(:)
    integer, allocatable, intent(in) :: img2CentCell(:)

    dcoord = 0.0_dp
    do iAtom1=1,this%nMovedAtom
       dcoord(:,this%indMovedAtom(iAtom1)) = this%movedVelo(:,iAtom1)
    end do

    sPrimeTmp = 0.0_dp
    RdotSprime = 0.0_dp

    !$OMP PARALLEL DO PRIVATE(iAtom1,iStart1,iEnd1,iSp1,nOrb1,sPrimeTmp2,iNeigh,iAtom2, &
    !$OMP& iAtom2f,iStart2,iEnd2,iSp2,nOrb2,sPrimeTmp,iDir) DEFAULT(SHARED) &
    !$OMP& SCHEDULE(RUNTIME)
    do iAtom1 = 1, this%nAtom
       iStart1 = iSquare(iAtom1)
       iEnd1 = iSquare(iAtom1+1)-1
       iSp1 = this%species(iAtom1)
       nOrb1 = orb%nOrbAtom(iAtom1)

       ! Onsite blocks
       if (this%tCalcOnsiteGradients) then
          sPrimeTmp2(:,:) = 0.0_dp
          do iDir = 1, 3
             sPrimeTmp2(:,:) = sPrimeTmp2 + &
                  &this%onsiteGrads(iDir, iAtom1, :, :) * dcoord(iDir, iAtom1)
          end do
          RdotSprime(iStart1:iEnd1,iStart1:iEnd1) = cmplx(sPrimeTmp2(1:nOrb1,1:nOrb1), 0, dp)
       end if

       ! Offsite blocks
       do iNeigh = 1, nNeighbourSK(iAtom1)
          iAtom2 = neighbourList%iNeighbour(iNeigh, iAtom1)
          iAtom2f = img2CentCell(iAtom2)
          iStart2 = iSquare(iAtom2f)
          iEnd2 = iSquare(iAtom2f+1)-1
          iSp2 = this%species(iAtom2f)
          nOrb2 = orb%nOrbAtom(iAtom2f)
          if (iAtom2f /= iAtom1) then
             call this%derivator%getFirstDeriv(sPrimeTmp, skOverCont, coord, this%species,&
                  & iAtom1, iAtom2, orb)

             sPrimeTmp2(:,:) = 0.0_dp
             do iDir=1,3
                sPrimeTmp2(:,:) = sPrimeTmp2 + &
                     &sPrimeTmp(:,:,iDir) * dcoord(iDir,iAtom1)
             end do
             RdotSprime(iStart2:iEnd2,iStart1:iEnd1) = &
                  & cmplx(sPrimeTmp2(1:nOrb2,1:nOrb1), 0, dp)

             sPrimeTmp2(:,:) = 0.0_dp
             do iDir=1,3
                sPrimeTmp2(:,:) = sPrimeTmp2 - &
                     &sPrimeTmp(:,:,iDir) * dcoord(iDir,iAtom2)
             end do
             RdotSprime(iStart1:iEnd1,iStart2:iEnd2) = &
                  & cmplx(transpose(sPrimeTmp2(1:nOrb2,1:nOrb1)), 0, dp)
          end if
       end do
    end do
    !$OMP END PARALLEL DO

  end subroutine getRdotSprime


  !> Calculates onsite gradients for non-adiabatic coupling
  subroutine getOnsiteGrads(this, skOverCont, orb)
    type(TElecDynamics), intent(inout) :: this
    type(OSlakoCont), intent(in) :: skOverCont
    type(TOrbitals), intent(in) :: orb
    real(dp) :: dist, uVects(3,3), vect(3), Stmp(2, orb%mOrb, orb%mOrb), Sder(orb%mOrb, orb%mOrb)
    real(dp) :: Stmp2(3, orb%mOrb, orb%mOrb)
    real(dp) :: interSKOver(getMIntegrals(skOverCont))
    integer :: iAt, iSp, nOrb, dir, iAux

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
          do iAux = 1, 2
             Stmp(iAux, :, :) = 0.0_dp
             vect = (-1)**(iAux+1) * uVects(dir, :)

             call getSKIntegrals(skOverCont, interSKOver, dist, iSp, iSp)
             call rotateH0(Stmp(iAux, :, :), interSKOver, vect(1), vect(2), vect(3), &
                  &iSp, iSp, orb)
          end do
          Sder(:,:) = (Stmp(1, :, :) - Stmp(2, :, :)) / (2.0_dp * dist)
          this%onsiteGrads(dir, iAt, 1:nOrb, 1:nOrb) = Sder(1:nOrb, 1:nOrb)
       end do
    end do

  end subroutine getOnsiteGrads

  !! Calculates properties per bond.
  !! If hamover = ham0 is energy
  !! If hamover = over is bond order
  subroutine pairWiseBondEO(this, EObond, rhoPrim, hamover, iSquare, iNeighbour, nNeighbourSK, &
       & img2CentCell, iSparseStart)
    type(TElecDynamics), intent(in), target :: this
    real(dp), intent(in) :: rhoPrim(:)
    real(dp), intent(in) :: hamover(:)
    real(dp), intent(out) :: EObond(this%nAtom, this%nAtom)
    integer, intent(in) :: iNeighbour(0:,:), nNeighbourSK(:)
    integer, allocatable, intent(in) :: iSparseStart(:,:), img2CentCell(:)
    integer, intent(in) :: iSquare(:)
    integer :: iAt1, iAt2, iAt2f, nOrb1, nOrb2, iOrig, iStart, iEnd, iNeigh, mOrb, iOrb, iOrb2

    EObond = 0.0_dp

    do iAt1 = 1, this%nAtom
       iOrb = iSquare(iAt1)
       nOrb1 = iSquare(iAt1+1) - iOrb
       do iNeigh = 0, nNeighbourSK(iAt1)
          iOrig = iSparseStart(iNeigh,iAt1) + 1
          iAt2 = iNeighbour(iNeigh, iAt1)
          iAt2f = img2CentCell(iAt2)
          iOrb2 = iSquare(iAt2f)
          nOrb2 = iSquare(iAt2f+1) - iOrb2

          EObond(iAt2, iAt1) = &
               & sum(rhoPrim(iOrig:iOrig+nOrb1*nOrb2-1) * &
               & hamover(iOrig:iOrig+nOrb1*nOrb2-1))
          EObond(iAt1, iAt2) = EObond(iAt2, iAt1)
       end do
    end do
  end subroutine pairWiseBondEO

end module timeprop_module
