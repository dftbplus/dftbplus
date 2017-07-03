!--------------------------------------------------------------------------------------------------!
!  DFTB+: general package for performing fast atomistic simulations                                !
!  Copyright (C) 2017  DFTB+ developers group                                                      !
!                                                                                                  !
!  See the LICENSE file for terms of usage and distribution.                                       !
!--------------------------------------------------------------------------------------------------!

#:include 'common.fypp'

!!* Global variables and initialization for the main program
!!* @todo Assignment (copy) operator for TNeighbors!!!
module initprogram
  use assert
  use inputdata_module
  use constants
  use periodic
  use accuracy
  use intrinsicpr
  use shortgamma
  use coulomb
  use message

  use mainio, only : SetEigVecsTxtOutput
  use mixer
  use simplemixer
  use andersonmixer
  use broydenmixer
  use diismixer

  use geoopt
  use conjgrad
  use steepdesc
  use gdiis

  use ranlux
  use mdcommon
  use mdintegrator
  use velocityverlet
  use thermostat
  use dummytherm
  use andersentherm
  use berendsentherm
  use nhctherm
  use tempprofile
  use numderivs2
  use lapackroutines
  use simplealgebra
  use nonscc
  use scc
  use sccinit
  use slakocont
  use repcont

  use fileid

  use spin, only: Spin_getOrbitalEquiv, ud2qm, qm2ud
  use dftbplusu

  use dispersions

  use thirdorder_module
  use linresp_module
  use stress
  use orbitalequiv
  use commontypes
  use sorting, only : heap_sort
  use fifo
  use linkedlist
  use xlbomd_module
  use etemp, only : Fermi
  use ipisocket
  implicit none


  character(*), parameter :: fChargeIn = "charges.bin"
  character(*), parameter :: fStopDriver = "stop_driver"
  character(*), parameter :: fStopSCC = "stop_scc"

  logical               :: tSCC            !* Is the calculation SCC?
  integer,  parameter   :: nCutoff = 1     !* Nr. of different cutoffs

  integer               :: nAtom           !* nr. of atoms
  integer               :: nAllAtom        !* nr. of all (image and orig) atoms
  integer, allocatable :: Img2CentCell(:) !* nr. of original atom in centre
  integer               :: nType           !* nr of different types (nAtom)

  type(TOrbitals), target :: orb

  integer               :: nOrb            !* nr. of orbitals in the system
  integer               :: nAllOrb         !* nr. of orbitals for all atoms
  integer, allocatable :: species(:)       !* types of the atoms (nAllAtom)
  integer,  allocatable :: species0(:)      !* type of the atoms (nAtom)
  real(dp), allocatable :: coord(:,:)      !* Coords of the atoms (3, nAllAtom)
  real(dp), allocatable, target :: coord0(:,:)   !* Coords (3, nAtom)
  real(dp), allocatable :: tmpCoords(:)    !* temporary array of coords
  real(dp), allocatable :: tmpWeight(:)    !* temporary weights
  real(dp), allocatable :: tmp3Coords(:,:) !* temporary array of coords (3,:)
  logical               :: tPeriodic       !* if calculation is periodic
  logical               :: tShowFoldedCoord!* Should central cell coordinates
  !* be output?
  integer :: forceType  ! How to calculate forces
  logical               :: tFracCoord      !* are atomic coordinates fractional?
  real(dp)              :: sccTol          !* Tollerance for SCC cycle

  real(dp), allocatable, target :: latVec(:,:)  !* lattice vectors as columns
  real(dp), allocatable, target :: recVec(:,:)  !* reciprocal vecs as columns

  !* original lattice vectors used for optimizing
  real(dp)              :: origLatVec(3,3)
  ! normalized vectors in those directions
  real(dp)              :: normOrigLatVec(3,3)

  real(dp), allocatable :: recVec2p(:,:)   !* reciprocal vectors in 2pi units
  real(dp)              :: CellVol          !* cell volume
  real(dp)              :: recCellVol       !* reciprocal cell volume
  real(dp), allocatable :: cellVec(:,:)    !* translation vecs for interacting
                                           !* image cells (3, nImgCell + 1)
  real(dp), allocatable :: rCellVec(:,:)   !* cell vectors in absolute units
  integer, allocatable :: iCellVec(:)     !* index in cellVec for each atom

  !!* ADT for neighbor parameters
  type(TNeighborList), allocatable, save :: neighborList
  integer,  allocatable :: nNeighbor(:)    !* nr. of neighbors for SK + rep
  integer, allocatable :: iPair(:,:)      !* H/S indexing array

  integer,  allocatable :: iAtomStart(:)   !* atom start pos for squared H/S

  real(dp), allocatable, target :: hubbU(:,:)      !* Hubbard Us (orbital, atom)
  real(dp), allocatable :: atomEigVal(:,:) !* self energy (orbital, atom)
  real(dp), allocatable :: referenceN0(:,:)!* reference n_0 charges for each
                                           !* atom
  real(dp), allocatable :: mass(:)         !* list of atomic masses
  real(dp), allocatable :: speciesMass(:)  !* list of atomic masses for each species

  type(OSlakoCont)  :: skHamCont
  type(OSlakoCont) :: skOverCont
  type(ORepCont) :: pRepCont
  real(dp) :: skCutoff
  real(dp) :: skRepCutoff

  real(dp)              :: mCutoff        !* longest pair interaction

  real(dp), allocatable :: ham(:,:)       !* Hamiltonian
  real(dp), allocatable :: iHam(:,:)      !* imaginary part of the Hamiltonian
  real(dp), allocatable :: chargePerShell(:,:,:)
  real(dp), allocatable :: over(:)        !* Overlap


  integer               :: nKPoint        !* nr. of K-points
  real(dp), allocatable :: kPoint(:,:)    !* K-points
  real(dp), allocatable :: KWeight(:)     !* weight of the K-Points

  real(dp)              :: pressure        !* external pressure if periodic
  logical               :: tBarostat
  real(dp)              :: BarostatStrength

  logical               :: tRealHS        !* H and S are real

  real(dp), allocatable :: nEl(:)         !* nr. of electrons
  real(dp)              :: nEl0           !* Nr. of all electrons if neutral
  real(dp), allocatable :: spinW(:,:,:)       !* Spin W's    !'
  real(dp), allocatable :: xi(:,:)        !* Spin orbit constants

  logical               :: tDFTBU         !* is this a DFTB+U calculation?
  integer               :: nDFTBUfunc     !* Choice of orbital functional
  real(dp), allocatable :: UJ(:,:)        !* list of U-J for species
  integer, allocatable  :: nUJ(:)         !* How many U-J for each species
  integer, allocatable  :: niUJ(:,:)      !* number of l-values of U-J for each
  !* block
  integer, allocatable  :: iUJ(:,:,:)     !* l-values of U-J for each block

  real(dp) :: tempElec                    !* electron temperature
  logical :: tFillKSep                    !* If K points should filled
                                          !* separately
  logical :: tFixEf                       !* Fix Fermi energy at specified value
  real(dp) :: Ef(2)                       !* Fermi energy

  logical :: tSetFillingTemp              !* Filling temp updated by MD.
  integer  :: iDistribFn = 0              !* Choice of electron distribution
                                          !* function, defaults to Fermi
  real(dp) :: tempAtom                    !* atomic temperature
  real(dp) :: deltaT                      !* MD stepsize
  integer :: solver                       !* eigensolver
  integer :: nSCCIter                     !* number of SCC iterations
  integer :: nSpin                        !* Number of spin components, 1
                                          !* is unpolarised, 2 is polarised, 4
  !* is noncolinear / spin-orbit
  logical :: tSpin                        !* is this a spin polarized
  !* calculation?
  logical :: tSpinOrbit                   !* is there spin-orbit coupling
  logical :: tDualSpinOrbit               !* Use block like dual representation
  logical :: tImHam                       !* complex hamiltonian in real space
  logical :: t2Component                  !* is this a two component
  !* calculation (spin orbit or non-collinear spin)

  logical :: tSpinSharedEf ! Common Fermi level accross spin channels

  real(dp) :: almix

  logical               :: tGeoOpt          !* Geometry optimization needed?
  logical               :: tCoordOpt        !* optimize internal coordinates?
  logical               :: tLatOpt          !* optimize lattice constants?
  logical               :: tLatOptFixAng    !* Fix angles between lattice
  !!* vectors when optimizing?
  logical               :: tLatOptFixLen(3)
  logical               :: tLatOptIsotropic
  logical               :: tMD              !* Is this a MD calculation?
  logical               :: tDerivs          !* Is this a derivatives calc?
  logical               :: tMulliken        !* Do we need Mulliken charges?
  logical               :: tLocalise        !* Calculate localised orbitals?
  logical               :: tPipekMezey      !* Use PipekMezey localisation?
  logical               :: tPipekDense      !* use a dense algorithm
                                            !* for Pipek-Mezey
                                            !* localisation?
  real(dp), allocatable :: sparsePipekTols(:) !* tollerances if
                                              !* instead using a
                                              !* sparse version
  real(dp)              :: PipekTol         ! halting tollerance for
                                            ! localisation
  integer               :: PipekMaxIter     ! number of localisation iterations
  logical               :: tPrintMulliken   !* Do we need to show
                                            !* Mulliken charges?
  logical               :: tDipole          !* calculate an electric dipole?
  logical               :: tAtomicEnergy    !* Do we need atom resolved E?
  logical               :: tPrintEigVecs    !* Print out eigenvectors?
  logical               :: tPrintEigVecsTxt !* Store as a text file
  logical               :: tProjEigenvecs   !* Print eigenvector projections?
  logical               :: tForces          !* Do we need forces?
  logical               :: tPrintForces     !* are forces being returned
  integer               :: nMovedAtom       !* Number of moved atoms
  integer, allocatable  :: indMovedAtom(:)  !* Index of the moved atoms
  integer               :: nMovedCoord      !* Nr. of moved coordinates
  integer               :: nGeoSteps        !* Nr. of geo movements to do
  integer               :: nGeoConstr       !* Nr. of geometry constraints
  integer,  allocatable :: conAtom(:)       !* Index of constrained atoms
  real(dp), allocatable :: conVec(:,:)      !* Constraint vectors

  logical :: tSocket         !* use commands from a socket
  type(IpiSocketComm), allocatable :: socket


  character(lc) :: geoOutFile     !* File containing end geometry

  logical               :: tAppendGeo       !* Append geometries in the output?
  logical :: tConvrgForces
  character(mc), allocatable :: speciesName(:)

  real(dp)              :: random_pool(10)   !* pool of initial random numbers
  ! for future use. See comment in code at create(pRanlux, in this routine.

  logical            :: tInitialized = .false.
  private :: tInitialized

  type(OGeoOpt), allocatable :: pGeoCoordOpt  !* General geometry optimizer
  type(OGeoOpt), allocatable :: pGeoLatOpt    !* Geometry optimizer for lattice
                                          !* consts

  type(OMixer), allocatable :: pChrgMixer    !* Charge mixer

  type(ORanlux), allocatable, target :: randomGenerator !* Random number generator
  type(ORanlux), pointer :: pRandomGenerator

  type(OMDCommon), allocatable :: pMDFrame  !* MD Framework
  type(OMDIntegrator), allocatable :: pMDIntegrator !* MD integrator
  type(OTempProfile), allocatable, target :: temperatureProfile

  type(OnumDerivs), allocatable, target :: derivDriver

  !! Charge related variables
  real(dp), allocatable    :: q0(:, :, :)
  real(dp), allocatable    :: qShell0(:,:)
  real(dp), allocatable    :: qInput(:, :, :)
  real(dp), allocatable    :: qOutput(:, :, :)
  real(dp), allocatable    :: qBlockIn(:, :, :, :)
  real(dp), allocatable    :: qBlockOut(:, :, :, :)
  real(dp), allocatable    :: qiBlockIn(:, :, :, :)
  real(dp), allocatable    :: qiBlockOut(:, :, :, :)
  real(dp), allocatable    :: qInpRed(:), qOutRed(:), qDiffRed(:)
  integer, allocatable     :: iEqOrbitals(:,:,:)  !* Orbital equiv. relations
  integer :: nIneqOrb !* nr. of inequivalent orbitals
  integer :: nMixElements !* nr. of elements to go through the mixer - may
  ! contain reduced orbitals and also orbital blocks (if tDFTBU)
  integer, allocatable :: iEqBlockDFTBU(:,:,:,:) !* Orbital equivalency for
  !* orbital blocks
  integer, allocatable :: iEqBlockDFTBULS(:,:,:,:) !* Orbital equivalency for
  !* orbital blocks with spin-orbit


  !! External charges
  integer :: nExtChrg   !* Nr. of external charges
  logical :: tExtChrg   !* If external charges must be considered

  logical  :: tEField = .false. ! external electric field
  real(dp) :: EFieldStrength = 0.0_dp ! field strength
  real(dp) :: EfieldVector(3) = 0.0_dp ! field direction
  logical  :: tTDEfield = .false. ! time dependent
  real(dp) :: EfieldOmega = 0.0_dp ! angular frequency
  integer  :: EfieldPhase = 0 ! phase of field at step 0

  ! PDOS projection
  type(listIntR1), save :: iOrbRegion
  type(listCharLc), save :: regionLabels
  integer, allocatable, save :: fdProjEig(:) ! file units for results

  !! Third order
  logical :: t3rd, t3rdFull
  type(ThirdOrder) :: thirdOrd

  !! Linear response
  logical :: tLinResp, tLinRespZVect
  logical :: tPrintExcitedEigVecs = .false.
  type(linresp), save :: lresp

  !! Other stuff
  logical :: tReadChrg    !* If initial charges/dens mtx. from external file.
  logical :: tWriteTagged !* produce tagged output?
  logical :: tWriteDetailedXML !* Produce detailed.xml
  logical :: tWriteResultsTag !* Produce detailed.tag
  logical :: tWriteDetailedOut !* Produce detailed.out
  logical :: tWriteBandDat !* Produce band.dat
  logical :: tWriteHS, tWriteRealHS  !* Should HS (square and real) be printed?
  logical :: tMinMemory, tStoreEigvecs


  integer :: runId !* Program run id

  integer :: restartFreq    !* Frequency for saving restart info

  logical :: tDispersion  !* If dispersion should be calculated
  logical :: tStress = .true. !* Can stress be calculated? - start by
  !* assuming it can
  class(DispersionIface), allocatable :: dispersion

  type(OFifoRealR2), allocatable :: storeEigvecsReal(:)
  type(OFifoCplxR2), allocatable :: storeEigvecsCplx(:)

  ! XLBOMD related parameters
  type(Xlbomd), allocatable :: xlbomdIntegrator
  logical :: tXlbomd

  type(NonSccDiff), save :: nonSccDeriv

  integer, parameter :: nInitNeighbor = 40  !* First guess for nr. of neighbors.
  private :: nInitNeighbor



contains


  !!* Initializes the variables in the module based on the parsed input
  !!* @param input Holds the parsed input data.
  subroutine initProgramVariables(input)
    type(inputData), intent(inout), target :: input

    !! Mixer related local variables
    integer  :: nGeneration
    real(dp) :: mixParam
    integer :: iMixer        !* mixer number
    type(OSimpleMixer), allocatable :: pSimpleMixer
    type(OAndersonMixer), allocatable :: pAndersonMixer
    type(OBroydenMixer), allocatable :: pBroydenMixer
    type(ODIISMixer), allocatable :: pDIISMixer

    !! Geometry optimizer related local variables
    type(OConjGrad), allocatable :: pConjGrad    !* Conjugate gradient driver
    type(OSteepDesc), allocatable :: pSteepDesc  !* Steepest descent driver
    type(OConjGrad), allocatable :: pConjGradLat   !* Conjugate gradient driver
    type(OSteepDesc), allocatable :: pSteepDescLat !* Steepest descent driver
    type(ODIIS), allocatable :: pDIIS    !* gradient DIIS driver

    !! MD related local variables
    type(OThermostat), allocatable :: pThermostat
    type(ODummyThermostat), allocatable :: pDummyTherm
    type(OAndersenThermostat), allocatable :: pAndersenTherm
    type(OBerendsenThermostat), allocatable :: pBerendsenTherm
    type(ONHCThermostat), allocatable :: pNHCTherm

    type(OVelocityVerlet), allocatable :: pVelocityVerlet
    type(OTempProfile), pointer :: pTempProfile

    integer :: ind, ii, jj, kk, iS, iAt, iSp, iSh, iOrb
    integer :: iStart, iEnd

    !! Dispersion
    type(DispSlaKirk), allocatable :: slaKirk
    type(DispUFF), allocatable :: uff
  #:if WITH_DFTD3
    type(DispDftD3), allocatable :: dftd3
  #:endif

    character(lc) :: strTmp, strTmp2
    logical :: tFirst ! flag to check for first cycle through a loop
    integer  :: iSeed, timeValues(8)
    real(dp) :: rTmp

    logical :: tExist ! Flag if some files do exist or not

    !! Orbital equivalency for SCC and Spin
    integer, allocatable :: iEqOrbSCC(:,:,:), iEqOrbSpin(:,:,:)
    !! Orbital equivalency for orbital potentials
    integer, allocatable :: iEqOrbDFTBU(:,:,:)

    !! Damped interactions
    logical, allocatable, target :: tDampedShort(:)
    type(ThirdOrderInp) :: thirdInp

    !! PDOS stuff
    integer :: iReg, nAtomRegion, nOrbRegion, iTmp
    integer, allocatable :: iAtomRegion(:)
    integer :: valshape(1)
    character(lc) :: tmpStr
    integer, allocatable :: tmpir1(:)

    type(TSCCInit), allocatable :: sccInit

    ! Used for indexing linear response
    integer :: homoLoc(1)

    @:ASSERT(input%tInitialized)

    !! Basic variables
    tSCC = input%ctrl%tScc
    tDFTBU = input%ctrl%tDFTBU
    tSpin = input%ctrl%tSpin
    if (tSpin) then
      nSpin = 2
    else
      nSpin = 1
    end if
    tSpinSharedEf = input%ctrl%tSpinSharedEf
    tSpinOrbit = input%ctrl%tSpinOrbit
    tDualSpinOrbit = input%ctrl%tDualSpinOrbit
    t2Component = input%ctrl%t2Component


    if (t2Component) then
      nSpin = 4
    end if

    if (nSpin /= 2 .and. tSpinSharedEf) then
      call error("Colinear spin polarization required for shared Ef over spin&
          & channels")
    end if


    sccTol = input%ctrl%sccTol
    nAtom = input%geom%nAtom
    nType = input%geom%nSpecies
    tPeriodic = input%geom%tPeriodic
    tShowFoldedCoord = input%ctrl%tShowFoldedCoord
    if (tShowFoldedCoord .and. .not. tPeriodic) then
      call error("Folding coordinates back into the central cell is&
          & meaningless for molecular boundary conditions!")
    end if
    tFracCoord = input%geom%tFracCoord
    solver = input%ctrl%iSolver
    if (tSCC) then
      nSCCIter = input%ctrl%maxIter
    else
      nSCCIter = 1
    end if

    if (tPeriodic) then
      allocate(latVec(3, 3))
      @:ASSERT(all(shape(input%geom%latVecs) == shape(latVec)))
      latVec(:,:) = input%geom%latVecs(:,:)
      allocate(recVec(3, 3))
      allocate(recVec2p(3, 3))
      recVec2p = latVec(:,:)
      call matinv(recVec2p)
      recVec2p = reshape(recVec2p, (/3, 3/), order=(/2, 1/))
      recVec = 2.0_dp * pi * recVec2p
      CellVol = abs(determinant33(latVec))
      recCellVol = abs(determinant33(recVec))
    else
      allocate(latVec(0, 0))
      allocate(recVec(0, 0))
      allocate(recVec2p(0, 0))
      CellVol = 0.0_dp
      recCellVol = 0.0_dp
    end if

    orb = input%slako%orb

    !! Slater-Koster tables
    skHamCont = input%slako%skHamCont
    skOverCont = input%slako%skOverCont
    pRepCont = input%slako%repCont

    allocate(atomEigVal(orb%mShell, nType))
    @:ASSERT(size(input%slako%skSelf, dim=1) == orb%mShell)
    @:ASSERT(size(input%slako%skSelf, dim=2) == size(atomEigVal, dim=2))
    atomEigVal(:,:) = input%slako%skSelf(1:orb%mShell, :)

    @:ASSERT(size(input%slako%skOcc, dim=1) >= orb%mShell)
    allocate(referenceN0(orb%mShell, nType))
    referenceN0(:,:) = input%slako%skOcc(1:orb%mShell, :)
    @:ASSERT(size(input%slako%mass) == nType)
    allocate(speciesMass(nType))
    speciesMass(:) = input%slako%mass(:)

    ! Spin W's !'
    if (allocated(input%ctrl%spinW)) then
      allocate(spinW(orb%mShell, orb%mShell, nType))
      spinW(:,:,:) = 0.0_dp
      do iSp = 1, nType
        do jj = 1, orb%nShell(iSp)
          do kk = 1, orb%nShell(iSp)
            spinW(jj, kk, iSp) = input%ctrl%spinW(jj, kk, iSp)
          end do
        end do
      end do
    else
      allocate(spinW(0,0,0))
    end if

    if (tSpinOrbit) then
      allocate(xi(orb%mShell,nType))
      xi(:,:) = 0.0_dp
      do iSp=1,nType
        do jj=1, orb%nShell(iSp)
          xi(jj,iSp)=input%ctrl%xi(jj,iSp)
        end do
      end do
    else
      allocate(xi(0,0))
    end if

    ! DFTB+U parameters
    if (tDFTBU) then
      nDFTBUfunc = input%ctrl%DFTBUfunc
      allocate(UJ(size(input%ctrl%UJ,dim=1),size(input%ctrl%UJ,dim=2)))
      allocate(nUJ(size(input%ctrl%nUJ)))
      allocate(niUJ(size(input%ctrl%niUJ,dim=1),size(input%ctrl%niUJ,dim=2)))
      allocate(iUJ(size(input%ctrl%iUJ,dim=1), size(input%ctrl%iUJ,dim=2),&
          & size(input%ctrl%iUJ,dim=3)))

      UJ(:,:) = input%ctrl%UJ(:,:)
      nUJ(:) = input%ctrl%nUJ(:)
      niUJ(:,:) = input%ctrl%niUJ(:,:)
      iUJ(:,:,:) = input%ctrl%iUJ(:,:,:)
      do iSp = 1, nType
        do jj = 1, nUJ(iSp)
          if (niUJ(jj,iSp)>1) then
            call heap_sort(iUJ(1:niUJ(jj,iSp),jj,iSp))
          end if
        end do
      end do
    else
      allocate(UJ(0,0))
      allocate(nUJ(0))
      allocate(niUJ(0,0))
      allocate(iUJ(0,0,0))
    end if

    !! Cutoffs from SlaKo and repulsive
    skCutoff = max(getCutoff(skHamCont), getCutoff(skOverCont))
    skRepCutoff = max(skCutoff, getCutoff(pRepCont))
    mCutoff = skRepCutoff

    !! Get species names and output file
    geoOutFile = input%ctrl%outFile
    allocate(speciesName(size(input%geom%speciesNames)))
    speciesName(:) = input%geom%speciesNames(:)

    do iSp = 1, nType
      do jj = iSp+1, nType
        if (speciesName(iSp) == speciesName(jj)) then
          write(tmpStr,"('Duplicate identical species labels in the geometry:&
              & ',A)")speciesName(iSp)
              call error(tmpStr)
        end if
      end do
    end do


    !! Initialise the SCC module (the two copies of the Hubbard Us are rather
    !! artifical, since the copy for the main program is only used for dumping
    !! into the tagged format for autotest)
    allocate(hubbU(orb%mShell, nType))
    @:ASSERT(size(input%slako%skHubbU, dim=1) >= orb%mShell)
    @:ASSERT(size(input%slako%skHubbU, dim=2) == nType)
    hubbU(:,:) = input%slako%skHubbU(1:orb%mShell, :)
    if (allocated(input%ctrl%hubbU)) then
      where (input%ctrl%hubbU > 0.0_dp)
        hubbU = input%ctrl%hubbU
      end where
    end if
    if (tSCC) then
      allocate(sccInit)
      sccInit%orb => orb
      if (tPeriodic) then
        sccInit%latVecs = latVec
        sccInit%recVecs = recVec
        sccInit%volume = CellVol
      end if
      sccInit%hubbU = hubbU
      allocate(tDampedShort(nType))
      if (input%ctrl%tDampH) then
        tDampedShort = (speciesMass < 3.5_dp * amu__au)
        !tDampedShort(:) = (speciesName == "H" .or. speciesName == "h")
      else
        tDampedShort(:) = .false.
      end if
      sccInit%tDampedShort = tDampedShort
      sccInit%dampExp = input%ctrl%dampExp
      nExtChrg = input%ctrl%nExtChrg
      tExtChrg = (nExtChrg > 0)
      if (tExtChrg) then
        if (.not.tSCC) then
          call error("External charges can only be used in an SCC calculation")
        end if
        tStress = .false. ! Stress calculations not allowed
        @:ASSERT(size(input%ctrl%extChrg, dim=1) == 4)
        @:ASSERT(size(input%ctrl%extChrg, dim=2) == nExtChrg)
        sccInit%extCharges = input%ctrl%extChrg
        if (allocated(input%ctrl%extChrgBlurWidth)) then
          sccInit%blurWidths = input%ctrl%extChrgblurWidth
        end if
      end if
      if (allocated(input%ctrl%chrgConstr)) then
        @:ASSERT(all(shape(input%ctrl%chrgConstr) == (/ nAtom, 2 /)))
        if (any(abs(input%ctrl%chrgConstr(:,2)) > epsilon(1.0_dp))) then
          sccInit%chrgConstraints = input%ctrl%chrgConstr
        end if
      end if

      if (allocated(input%ctrl%thirdOrderOn)) then
        @:ASSERT(tSCC)
        @:ASSERT(all(shape(input%ctrl%thirdOrderOn) == (/ nAtom, 2 /)))
        sccInit%thirdOrderOn = input%ctrl%thirdOrderOn
      end if

      sccInit%ewaldAlpha = input%ctrl%ewaldAlpha
      call init_SCC(sccInit)
      deallocate(sccInit)
      mCutoff = max(mCutoff, getSCCCutoff())

      if (input%ctrl%t3rd .and. input%ctrl%tOrbResolved) then
        call error("Onsite third order DFTB only compatible with orbital non&
            & resolved SCC")
      end if

      ! Initialize full 3rd order module
      t3rd = input%ctrl%t3rd
      t3rdFull = input%ctrl%t3rdFull
      if (t3rdFull) then
        @:ASSERT(tSCC)
        thirdInp%orb => orb
        thirdInp%hubbUs = hubbU
        thirdInp%hubbUDerivs = input%ctrl%hubDerivs
        allocate(thirdInp%damped(nType))
        thirdInp%damped(:) = tDampedShort
        thirdInp%dampExp = input%ctrl%dampExp
        thirdInp%shellResolved = input%ctrl%tOrbResolved
        call ThirdOrder_init(thirdOrd, thirdInp)
        mCutoff = max(mCutoff, thirdOrd%getCutoff())
      end if
    end if

    !! Initial coordinates
    allocate(coord0(3, nAtom))
    @:ASSERT(all(shape(coord0) == shape(input%geom%coords)))
    coord0(:,:) = input%geom%coords(:,:)
    allocate(species0(nAtom))
    @:ASSERT(all(shape(species0) == shape(input%geom%species)))
    species0(:) = input%geom%species(:)

    allocate(mass(nAtom))
    mass = speciesMass(species0)
    if (allocated(input%ctrl%masses)) then
      @:ASSERT(size(input%ctrl%masses) == nAtom)
      where (input%ctrl%masses >= 0.0_dp)
        mass = input%ctrl%masses
      end where
    end if

    if (tPeriodic) then
      !! Make some guess for the nr. of all interacting atoms
      nAllAtom = int((real(nAtom, dp)**(1.0_dp/3.0_dp) + 3.0_dp)**3)
    else
      nAllAtom = nAtom
    end if
    allocate(coord(3, nAllAtom))
    allocate(species(nAllAtom))
    allocate(img2CentCell(nAllAtom))
    allocate(iCellVec(nAllAtom))
    allocate(iAtomStart(nAtom + 1))
    call buildSquaredAtomIndex(iAtomStart, orb)

    !! Intialize Hamilton and overlap
    if (tSCC) then
      allocate(chargePerShell(orb%mShell,nAtom,nSpin))
    else
       allocate(chargePerShell(0,0,0))
    end if
    allocate(ham(0, nSpin))
    allocate(iHam(0, nSpin))
    allocate(over(0))
    allocate(iPair(0, nAtom))

    !! Brillouin zone sampling
    if (tPeriodic) then
      nKPoint = input%ctrl%nKPoint
      allocate(kPoint(3, nKPoint))
      allocate(kWeight(nKPoint))
      @:ASSERT(all(shape(kPoint) == shape(input%ctrl%KPoint)))
      @:ASSERT(all(shape(kWeight) == shape(input%ctrl%kWeight)))
      kPoint(:,:) = input%ctrl%KPoint(:,:)
      if (sum(input%ctrl%kWeight(:)) < epsilon(1.0_dp)) then
        call error("Sum of k-point weights should be greater than zero!")
      end if
      kWeight(:) = input%ctrl%kWeight(:) / sum(input%ctrl%kWeight(:))
    else
      nKPoint = 1
      allocate(kPoint(3, nKPoint))
      allocate(kWeight(nKpoint))
      kPoint(:,1) = 0.0_dp
      kWeight(1) = 1.0_dp
    end if

    if ((.not. tPeriodic) .or. (nKPoint == 1 .and. &
        &all(kPoint(:, 1) == (/ 0.0_dp, 0.0_dp, 0.0_dp /)))) then
      tRealHS = .true.
    else
      tRealHS = .false.
    end if

    !! Other usefull quantities
    nOrb = orb%nOrb

    if (nSpin == 4) then
      allocate(nEl(1))
    else
      allocate(nEl(nSpin))
    end if

    nEl0 = 0.0_dp
    do ii = 1, nAtom
      nEl0 = nEl0 + sum(input%slako%skOcc(1:orb%nShell(species0(ii)), &
          & species0(ii)))
    end do
    nEl(:) = 0.0_dp
    if (nSpin == 1 .or. nSpin == 4) then
      nEl(1) = nEl0 - input%ctrl%nrChrg
      if(ceiling(nEl(1)) > 2.0_dp*nOrb) then
        call error("More electrons than basis functions!")
      end if
    else
      nEl(1) = 0.5_dp * (nEl0 - input%ctrl%nrChrg + input%ctrl%nrSpinPol)
      nEl(2) = 0.5_dp * (nEl0 - input%ctrl%nrChrg - input%ctrl%nrSpinPol)
      if (any(ceiling(nEl(:)) > nOrb)) then
        call error("More electrons than basis functions!")
      end if
    end if

    if (.not.all(nEl(:) >= 0.0_dp)) then
      call error("Less than 0 electrons!")
    end if

    iDistribFn = input%ctrl%iDistribFn
    tempElec = input%ctrl%tempElec
    tFixEf = input%ctrl%tFixEf
    if (tFixEf) then
      Ef = input%ctrl%Ef
    else
      Ef = 0.0_dp
    end if
    tSetFillingTemp = input%ctrl%tSetFillingTemp
    tFillKSep = input%ctrl%tFillKSep
    tempAtom = input%ctrl%tempAtom
    deltaT = input%ctrl%deltaT

    tImHam = tDualSpinOrbit .or. (tSpinOrbit .and. tDFTBU) ! .or. tBField

    !! Create equivalency relations
    if (tSCC) then
      allocate(iEqOrbitals(orb%mOrb, nAtom, nSpin))
      allocate(iEqOrbSCC(orb%mOrb, nAtom, nSpin))
      call SCC_getOrbitalEquiv(orb, species0, iEqOrbSCC)
      if (nSpin == 1) then
        iEqOrbitals(:,:,:) = iEqOrbSCC(:,:,:)
      else
        allocate(iEqOrbSpin(orb%mOrb, nAtom, nSpin))
        call Spin_getOrbitalEquiv(orb, species0, iEqOrbSpin)
        call OrbitalEquiv_merge(iEqOrbSCC, iEqOrbSpin, orb, iEqOrbitals)
        deallocate(iEqOrbSpin)
      end if
      deallocate(iEqOrbSCC)
      nIneqOrb = maxval(iEqOrbitals)
      nMixElements = nIneqOrb
      if (tDFTBU) then
        allocate(iEqOrbSpin(orb%mOrb, nAtom, nSpin))
        allocate(iEqOrbDFTBU(orb%mOrb, nAtom, nSpin))
        call DFTBplsU_getOrbitalEquiv(iEqOrbDFTBU,orb, species0, nUJ, niUJ, iUJ)
        call OrbitalEquiv_merge(iEqOrbitals, iEqOrbDFTBU, orb, iEqOrbSpin)
        iEqOrbitals(:,:,:) = iEqOrbSpin(:,:,:)
        nIneqOrb = maxval(iEqOrbitals)
        deallocate(iEqOrbSpin)
        deallocate(iEqOrbDFTBU)
        allocate(iEqBlockDFTBU(orb%mOrb, orb%mOrb, nAtom, nSpin))
        call DFTBU_blockIndx(iEqBlockDFTBU, nIneqOrb, orb, species0, &
            & nUJ, niUJ, iUJ)
        nMixElements = max(nMixElements,maxval(iEqBlockDFTBU)) ! as
        !  iEqBlockDFTBU does not include diagonal elements, so in the case of
        !  a purely s-block DFTB+U calculation, maxval(iEqBlockDFTBU) would
        !  return 0
        if (tImHam) then
          allocate(iEqBlockDFTBULS(orb%mOrb, orb%mOrb, nAtom, nSpin))
          call DFTBU_blockIndx(iEqBlockDFTBULS,nMixElements , orb, species0, &
            & nUJ, niUJ, iUJ)
          nMixElements = max(nMixElements,maxval(iEqBlockDFTBULS))
        end if
      end if
    else
      nIneqOrb = nOrb
      nMixElements = 0
    end if

    if (.not.tDFTBU) then
      allocate(iEqBlockDFTBU(0, 0, 0, 0))
    end if
    if (.not.(tDFTBU.and.tImHam)) then
      allocate(iEqBlockDFTBULS(0, 0, 0, 0))
    end if


    !! Initialize mixer
    !! (at the moment, the mixer does not need to know about the size of the
    !! vector to mix.)
    if (tSCC) then
      allocate(pChrgMixer)
      iMixer = input%ctrl%iMixSwitch
      nGeneration = input%ctrl%iGenerations
      mixParam = input%ctrl%almix
      select case (iMixer)
      case (1)
        allocate(pSimplemixer)
        call init(pSimpleMixer, mixParam)
        call init(pChrgMixer, pSimpleMixer)
      case (2)
        allocate(pAndersonMixer)
        if (input%ctrl%andersonNrDynMix > 0) then
          call init(pAndersonMixer, nGeneration, mixParam, &
              &input%ctrl%andersonInitMixing, input%ctrl%andersonDynMixParams, &
              &input%ctrl%andersonOmega0)
        else
          call init(pAndersonMixer, nGeneration, mixParam, &
              &input%ctrl%andersonInitMixing, omega0=input%ctrl%andersonOmega0)
        end if
        call init(pChrgMixer, pAndersonMixer)
      case (3)
        allocate(pBroydenMixer)
        call init(pBroydenMixer, nSCCIter, mixParam, input%ctrl%broydenOmega0,&
            & input%ctrl%broydenMinWeight, input%ctrl%broydenMaxWeight, input%ctrl%broydenWeightFac)
        call init(pChrgMixer, pBroydenMixer)
      case(4)
        allocate(pDIISMixer)
        call init(pDIISMixer,nGeneration, mixParam, input%ctrl%tFromStart)
        call init(pChrgMixer, pDIISMixer)
      case default
        call error("Unknown charge mixer type.")
      end select
    end if

    !! initialise in cases where atoms move
    tGeoOpt = input%ctrl%tGeoOpt
    tCoordOpt = input%ctrl%tCoordOpt
    tLatOpt = (input%ctrl%tLatOpt .and. tPeriodic)
    if (tLatOpt) then
      if (tExtChrg) then
        ! Stop as not sure, what to do with the coordinates of the
        ! external charges, when the lattice changes.
        call error("External charges and lattice optimisation can not be used &
            &together.")
      end if
    end if
    if (tLatOpt) then
      tLatOptFixAng = input%ctrl%tLatOptFixAng
      tLatOptFixLen = input%ctrl%tLatOptFixLen
      tLatOptIsotropic = input%ctrl%tLatOptIsotropic
      if (tLatOptFixAng .or. any(tLatOptFixLen) .or. tLatOptIsotropic) then
        origLatVec(:,:) = latVec(:,:)
        do ii = 1, 3
           normOrigLatVec(:,ii) = origLatVec(:,ii) &
                & / sqrt(sum(origLatVec(:,ii)**2))
        end do
      end if
    end if
    pressure = input%ctrl%pressure
    tBarostat = input%ctrl%tBarostat
    BarostatStrength = input%ctrl%BarostatStrength

    tSocket = allocated(input%ctrl%socketInput)
    if (tSocket) then
      write(*,*) "Initialising for socket communication to host ", &
          & trim(input%ctrl%socketInput%host)
      input%ctrl%socketInput%nAtom = nAtom
      socket = IpiSocketComm(input%ctrl%socketInput)
      tForces = .true.
      tGeoOpt = .false.
      tMD = .false.
    end if

    tAppendGeo = input%ctrl%tAppendGeo
    tConvrgForces = (input%ctrl%tConvrgForces .and. tSCC) ! no point if not SCC
    tMD = input%ctrl%tMD
    tDerivs = input%ctrl%tDerivs
    tPrintMulliken = input%ctrl%tPrintMulliken
    tEField = input%ctrl%tEfield ! external electric field
    tMulliken = input%ctrl%tMulliken .or. tPrintMulliken .or. tEField
    tAtomicEnergy = input%ctrl%tAtomicEnergy
    tPrintEigVecs = input%ctrl%tPrintEigVecs

    ! false if not set anywhere else
    call SetEigVecsTxtOutput(input%ctrl%tPrintEigVecsTxt &
        & .or. input%ctrl%tPipekMezey)

    ! Projection of eigenstates onto specific regions of the system
    tProjEigenvecs = input%ctrl%tProjEigenvecs
    if (tProjEigenvecs) then
      call init(iOrbRegion)
      call init(regionLabels)
      do iReg = 1, size(input%ctrl%tShellResInRegion)
        call elemShape(input%ctrl%iAtInRegion, valshape, iReg)
        nAtomRegion = valshape(1)
        allocate(iAtomRegion(nAtomRegion))
        call intoArray(input%ctrl%iAtInRegion, iAtomRegion, iTmp, iReg)
        if (input%ctrl%tOrbResInRegion(iReg) .or. input%ctrl%tShellResInRegion(iReg)) then

          if (input%ctrl%tOrbResInRegion(iReg)) then
            iSp = species0(iAtomRegion(1)) ! all atoms the same in the region
            @:ASSERT(all(species0(iAtomRegion) == iSp))
            nOrbRegion = nAtomRegion
            ! Create orbital index.
            allocate(tmpir1(nOrbRegion))
            do iOrb = 1, orb%nOrbSpecies(iSp)
              tmpir1 = 0
              ind = 1
              do iAt = 1, nAtomRegion
                tmpir1(ind) = iAtomStart(iAt) + iOrb - 1
                ind = ind + 1
              end do
              call append(iOrbRegion, tmpir1)
              write(tmpStr, "(A,'.',I0,'.',I0,'.out')") &
                  & trim(input%ctrl%RegionLabel(iReg)), orb%iShellOrb(iOrb,iSp), &
                  & iOrb-orb%posShell(orb%iShellOrb(iOrb,iSp),iSp) &
                  & -orb%angShell(orb%iShellOrb(iOrb,iSp),iSp)
              call append(regionLabels, tmpStr)
            end do
            deallocate(tmpir1)
          end if

          if (input%ctrl%tShellResInRegion(iReg)) then
            iSp = species0(iAtomRegion(1)) ! all atoms the same in the region
            @:ASSERT(all(species0(iAtomRegion) == iSp))
            ! Create a separate region for each shell. It will contain
            ! the orbitals of that given shell for each atom in the region.
            do iSh = 1, orb%nShell(iSp)
              nOrbRegion = nAtomRegion &
                  &* (orb%posShell(iSh + 1, iSp) - orb%posShell(iSh, iSp))
              ind = 1
              ! Create orbital index.
              allocate(tmpir1(nOrbRegion))
              do ii = 1, nAtomRegion
                iAt = iAtomRegion(ii)
                do jj = orb%posShell(iSh, iSp), orb%posShell(iSh + 1, iSp) - 1
                  tmpir1(ind) = iAtomStart(iAt) + jj - 1
                  ind = ind + 1
                end do
              end do
              call append(iOrbRegion, tmpir1)
              deallocate(tmpir1)
              write(tmpStr, "(A,'.',I0,'.out')") &
                  &trim(input%ctrl%RegionLabel(iReg)), iSh
              call append(regionLabels, tmpStr)
            end do
          end if

        else
          ! We take all orbitals from all atoms.
          nOrbRegion = 0
          do ii = 1, nAtomRegion
            nOrbRegion = nOrbRegion + orb%nOrbAtom(iAtomRegion(ii))
          end do
          ind = 1
          allocate(tmpir1(nOrbRegion))
          ! Create an index of the orbitals
          do ii = 1, nAtomRegion
            iAt = iAtomRegion(ii)
            do jj = 1, orb%nOrbAtom(iAt)
              tmpir1(ind) = iAtomStart(iAt) + jj - 1
              ind = ind + 1
            end do
          end do
          call append(iOrbRegion, tmpir1)
          deallocate(tmpir1)
          write(tmpStr, "(A,'.out')") trim(input%ctrl%RegionLabel(iReg))
          call append(regionLabels, tmpStr)
        end if
        deallocate(iAtomRegion)
      end do

      allocate(fdProjEig(len(iOrbRegion)))
      do ii = 1, len(iOrbRegion)
        fdProjEig(ii) = getFileId()
      end do
    else
      allocate(fdProjEig(0))
    end if

    tPrintForces = input%ctrl%tPrintForces
    tForces = input%ctrl%tForces .or. tPrintForces
    if (tSCC) then
      forceType = input%ctrl%forceType
    else
      if (input%ctrl%forceType /= 0) then
        call error("Invalid force evaluation method for non-SCC calculations.")
      end if
    end if
    if (forceType == 2 .and. tempElec > minTemp) then
       call error("This ForceEvaluation method requires the electron&
           & temperature to be zero")
    end if
    if (tForces) then
      select case(input%ctrl%iDerivMethod)
      case (1)
        ! set step size from input
        if (input%ctrl%deriv1stDelta < epsilon(1.0_dp)) then
          write(tmpStr, "(A,E12.4)") 'Too small value for finite difference &
              &step :', input%ctrl%deriv1stDelta
          call error(tmpStr)
        end if
        call NonSccDiff_init(nonSccDeriv, diffTypes%finiteDiff, &
            & input%ctrl%deriv1stDelta)
      case (2)
        call NonSccDiff_init(nonSccDeriv, diffTypes%richardson)
      end select
    end if

    ! requires stress to already be possible and it being a periodic calculation
    ! with forces
    tStress = ((tPeriodic .and. tForces).and.tStress)

    nMovedAtom = input%ctrl%nrMoved
    nMovedCoord = 3 * nMovedAtom

    if (input%ctrl%maxRun == -1) then
      nGeoSteps = huge(1)
    else
      nGeoSteps = input%ctrl%maxRun
    end if

    if (nMovedAtom > 0) then
      allocate(indMovedAtom(size(input%ctrl%indMovedAtom)))
      indMovedAtom(:) = input%ctrl%indMovedAtom(:)
    else
      allocate(indMovedAtom(0))
    end if

    allocate(pGeoCoordOpt)
    if (tCoordOpt) then
      allocate(tmpCoords(nMovedCoord))
      tmpCoords(1:nMovedCoord) = reshape(coord0(:, indMovedAtom), &
          & (/ nMovedCoord /))
      select case (input%ctrl%iGeoOpt)
      case(1)
        allocate(tmpWeight(nMovedCoord))
        tmpWeight(1:nMovedCoord) = 0.5_dp * deltaT**2 &
            & / reshape(spread(mass(indMovedAtom), 1, 3), &
            & (/nMovedCoord/))
        allocate(pSteepDesc)
        call init(pSteepDesc, size(tmpCoords), input%ctrl%maxForce, &
             & input%ctrl%maxAtomDisp,tmpWeight )
        deallocate(tmpWeight)
        call init(pGeoCoordOpt, pSteepDesc)
      case (2)
        allocate(pConjGrad)
        call init(pConjGrad, size(tmpCoords), input%ctrl%maxForce, &
            & input%ctrl%maxAtomDisp)
        call init(pGeoCoordOpt, pConjGrad)
      case (3)
        allocate(pDIIS)
        call init(pDIIS, size(tmpCoords), input%ctrl%maxForce, &
            & input%ctrl%deltaGeoOpt, input%ctrl%iGenGeoOpt)
        call init(pGeoCoordOpt, pDIIS)
      end select
      call reset(pGeoCoordOpt, tmpCoords)
    end if

    allocate(pGeoLatOpt)
    if (tLatOpt) then
      select case (input%ctrl%iGeoOpt)
      case(1)
        allocate(tmpWeight(9))
        tmpWeight = 1.0_dp
        allocate(pSteepDescLat)
        call init(pSteepDescLat, 9, input%ctrl%maxForce, &
            & input%ctrl%maxLatDisp, tmpWeight )
        deallocate(tmpWeight)
        call init(pGeoLatOpt, pSteepDescLat)
      case(2,3) ! use CG lattice for both DIIS and CG
        allocate(pConjGradLat)
        call init(pConjGradLat, 9, input%ctrl%maxForce, &
            & input%ctrl%maxLatDisp)
         call init(pGeoLatOpt, pConjGradLat)
      end select
      if (tLatOptIsotropic ) then ! optimization uses scaling factor
                                  ! of unit cell
        call reset( pGeoLatOpt, &
            &(/1.0_dp,0.0_dp,0.0_dp,0.0_dp,0.0_dp,0.0_dp,0.0_dp,0.0_dp,0.0_dp/))
      else if (tLatOptFixAng) then
        call reset( pGeoLatOpt, & ! optimization uses scaling factor
                                  ! of lattice vectors
            &(/1.0_dp,1.0_dp,1.0_dp,0.0_dp,0.0_dp,0.0_dp,0.0_dp,0.0_dp,0.0_dp/))
      else
        call reset( pGeoLatOpt, reshape(latVec, (/ 9 /)) )
      end if
    end if

    if (.not.(tGeoOpt.or.tMD.or.tSocket)) then
      nGeoSteps = 0
    end if

    !! Initialize constraints
    nGeoConstr = input%ctrl%nrConstr
    if (nGeoConstr > 0) then
      allocate(conAtom(input%ctrl%nrConstr))
      allocate(conVec(3, input%ctrl%nrConstr))
      @:ASSERT(all(shape(conAtom) == shape(input%ctrl%conAtom)))
      @:ASSERT(all(shape(conVec) == shape(input%ctrl%conVec)))
      conAtom(:) = input%ctrl%conAtom(:)
      conVec(:,:) = input%ctrl%conVec(:,:)
      do ii = 1, nGeoConstr
        conVec(:,ii) = conVec(:,ii) / sqrt(sum(conVec(:,ii)**2))
      end do
    else
      allocate(conAtom(0))
      allocate(conVec(3, 0))
    end if

    !! Dispersion
    tDispersion = allocated(input%ctrl%dispInp)
    if (tDispersion) then
      if (allocated(input%ctrl%dispInp%slakirk)) then
        tStress = .false.
        if (tLatOpt) then
          call error("Sorry, lattice optimisation and Slater-Kirkwood type&
              & dispersion can not be used together")
        end if
        if (tBarostat) then
          call error("Sorry, barostatic MD and Slater-Kirkwood type&
              & dispersion can not be used together")
        end if
        allocate(slaKirk)
        if (tPeriodic) then
          call DispSlaKirk_init(slaKirk, input%ctrl%dispInp%slakirk, &
              & latVec)
        else
          call DispSlaKirk_init(slaKirk, input%ctrl%dispInp%slakirk)
        end if
        call move_alloc(slaKirk, dispersion)

      elseif (allocated(input%ctrl%dispInp%uff)) then
        allocate(uff)
        if (tPeriodic) then
          call DispUff_init(uff, input%ctrl%dispInp%uff, nAtom, species0, &
              & latVec)
        else
          call DispUff_init(uff, input%ctrl%dispInp%uff, nAtom)
        end if
        call move_alloc(uff, dispersion)

    #:if WITH_DFTD3
      elseif (allocated(input%ctrl%dispInp%dftd3)) then
        allocate(dftd3)
        if (tPeriodic) then
          call DispDftD3_init(dftd3, input%ctrl%dispInp%dftd3, nAtom, &
              & species0, speciesName, latVec)
        else
          call DispDftD3_init(dftd3, input%ctrl%dispInp%dftd3, nAtom, &
              & species0, speciesName)
        end if
        call move_alloc(dftd3, dispersion)
    #:endif
      end if
      mCutoff = max(mCutoff, dispersion%getRCutoff())

    end if

    if (input%ctrl%nrChrg == 0.0_dp .and. (.not.tPeriodic) .and. tMulliken) then
      tDipole = .true.
    else
      tDipole = .false.
    end if

    tLocalise = input%ctrl%tLocalise
    if (tLocalise .and. tStoreEigvecs) then
      call error("Localisation of electronic states currently unsupported when&
          & storing data on disk due to fifo limitations")
    end if
    if (tLocalise .and. (nSpin > 2 .or. t2Component)) then
      call error("Localisation of electronic states currently unsupported for&
          & non-collinear and spin orbit calculations")
    end if
    tPipekMezey = input%ctrl%tPipekMezey
    PipekTol =  input%ctrl%PipekTol
    PipekMaxIter =  input%ctrl%PipekMaxIter
    tPipekDense = input%ctrl%tPipekDense
    if (.not.tPipekDense.and.tPipekMezey) then
      allocate(sparsePipekTols(size(input%ctrl%sparsePipekTols)))
      sparsePipekTols(:) = input%ctrl%sparsePipekTols(:)
      if (any(sparsePipekTols < epsilon(0.0_dp))) then
        call error('Tollerances for sparse Pipek-Mezey localisation too small.')
      end if
    else
      allocate(sparsePipekTols(0))
    end if

    ! Linear response
    tLinResp = input%ctrl%lrespini%tInit

    if (tLinResp) then

      ! input sanity checking
    #:if not WITH_ARPACK
      call error("This binary has been compiled without support for linear response calculations.")
    #:endif
      if (.not. tSCC) then
        call error("Linear response excitation requires SCC=Yes")
      end if
      if (nspin > 2) then
        call error("Linear reponse does not work with non-colinear spin polarization yet")
      elseif (tSpin .and. tForces) then
        call error("excited state relaxation is not implemented yet for spin-polarized systems")
      elseif (tPeriodic .and. .not.tRealHS) then
        call error("Linear response only works with non-periodic or gamma-point molecular&
            & crystals")
      elseif (tSpinOrbit) then
        call error("Linear response does not support spin orbit coupling at the moment.")
      elseif (tDFTBU) then
        call error("Linear response does not support LDA+U yet")
      elseif (input%ctrl%tOrbResolved) then
        call error("Linear response does not support orbital resolved scc yet")
      end if
      if (tempElec > 0.0_dp .and. tForces) then
        write(tmpStr, "(A,E12.4,A)")"Excited state forces are not implemented yet for fractional&
            & occupations, kT=", tempElec/Boltzmann,"K"
        call warning(tmpStr)
      end if

      if (input%ctrl%lrespini%nstat == 0) then
        if (input%ctrl%lrespini%tMulliken) then
          call error("Muliken analysis only available for StateOfInterest non zero.")
        end if
        if (tForces) then
          call error("Excited forces only available for StateOfInterest non zero.")
        end if
        if (input%ctrl%lrespini%tPrintEigVecs .or. input%ctrl%lrespini%tCoeffs) then
          call error("Excited eigenvectors only available for StateOfInterest non zero.")
        end if
      end if
      if (input%ctrl%lrespini%energyWindow < 0.0_dp) then
        call error("Negative energy window for excitations")
      end if

      ! Hubbard U and spin constants for excitations (W only needed for triplet/spin polarised)
      allocate(input%ctrl%lrespini%HubbardU(nType))
      allocate(input%ctrl%lrespini%spinW(nType))
      input%ctrl%lrespini%HubbardU = 0.0_dp
      input%ctrl%lrespini%spinW = 0.0_dp

      ! calculate linear response Gamma values from HOAO shell Hubbard U (non
      ! shell resolved)
      do iSp = 1, nType
        homoLoc = maxloc(atomEigVal(:orb%nShell(iSp), iSp), &
            & mask=input%slako%skOcc(:orb%nShell(iSp), iSp) > 0.0_dp)
        input%ctrl%lrespini%HubbardU(iSp) = hubbU(homoLoc(1), iSp)
      end do

      ! and atomic HOAO spin W value if needed
      input%ctrl%lrespini%spinW(:) = 0.0_dp
      select case(input%ctrl%lrespini%sym)
      case("S")
        ! Singlet case, no need for spin constants
      case("T","B"," ")
        ! triplet or spin-polarised
        do iSp = 1, nType
          homoLoc = maxloc(atomEigVal(:orb%nShell(iSp), iSp), &
              & mask=input%slako%skOcc(:orb%nShell(iSp), iSp) > 0.0_dp)
          input%ctrl%lrespini%spinW(iSp) = spinW(homoLoc(1), homoLoc(1), iSp)
        end do
      case default
        call error("Unknown excitation type requested")
      end select

      tPrintExcitedEigVecs = input%ctrl%lrespini%tPrintEigVecs
      tLinRespZVect = (input%ctrl%lrespini%tMulliken .or. tForces &
          & .or. input%ctrl%lrespini%tCoeffs .or. tPrintExcitedEigVecs)

      call init(lresp, input%ctrl%lrespini, nAtom, nEl(1), orb)

    end if

    !! Generate a random id for the run. Seed with system time if possible.
    call system_clock(iSeed)
    !! Try date_and_time if system_clock does not work properly.
    if (iSeed < 1) then
      call date_and_time(values=timeValues)
      if (timeValues(5) >= 0) then
        iSeed = 1000 * (60 * (60 * timeValues(5) + timeValues(6)) &
            &+ timeValues(7)) + timeValues(8)
      end if
    end if
    !! Try  default seeder if attempts to use the clock failed.
    if (iSeed < 1) then
      call random_seed
      call random_number(rTmp)
      ! Make sure seed > 0.
      iSeed = int(real(huge(iSeed) - 1, dp) * rTmp) + 1
    end if
    allocate(randomGenerator)
    call init(randomGenerator, initSeed=iSeed)
    call getRandom(randomGenerator, rTmp)
    runId = int(real(huge(runId) - 1, dp) * rTmp) + 1
    deallocate(randomGenerator)

    !! Create random generator and pull off first 10
    !! random numbers to avoid disturbing the subsequent sequence.
    !! random_pool can then be used to seed other generators if needed,
    !! with a supply of random numbers controlled from the initial seed.
    !! If the size of random_pool is changed then reproducibility of
    !! the random numbers if initialised from a seed is lost.
    iSeed = input%ctrl%iSeed
    if (iSeed < 1) then
      iSeed = runId     ! No seed specified, use random runId
    end if
    allocate(randomGenerator)
    call init(randomGenerator, 3, iSeed)
    call getRandom(randomGenerator,random_pool(:))
    pRandomGenerator => randomGenerator


    !! MD stuff
    if (tMD) then
      !! Create MD framework.
      allocate(pMDFrame)
      call init(pMDFrame, nMovedAtom, nAtom, input%ctrl%tMDstill)

      !! Create temperature profile, if thermostat is not the dummy one
      if (input%ctrl%iThermostat /= 0) then
        allocate(temperatureProfile)
        call init(temperatureProfile, input%ctrl%tempMethods, input%ctrl%tempSteps,&
            & input%ctrl%tempValues)
        pTempProfile => temperatureProfile
      else
        nullify(pTempProfile)
      end if

      !! Create thermostat
      allocate(pThermostat)
      select case (input%ctrl%iThermostat)
      case (0) ! No thermostat
        allocate(pDummyTherm)
        call init(pDummyTherm, tempAtom, mass(indMovedAtom), pRandomGenerator, pMDFrame)
        call init(pThermostat, pDummyTherm)
      case (1) ! Andersen thermostat
        allocate(pAndersenTherm)
        call init(pAndersenTherm, pRandomGenerator, mass(indMovedAtom), pTempProfile,&
            & input%ctrl%tRescale, input%ctrl%wvScale, pMDFrame)
        call init(pThermostat, pAndersenTherm)
      case (2) ! Berendsen thermostat
        allocate(pBerendsenTherm)
        call init(pBerendsenTherm, pRandomGenerator, mass(indMovedAtom), pTempProfile,&
            & input%ctrl%wvScale, pMDFrame)
        call init(pThermostat, pBerendsenTherm)
      case (3) ! Nose-Hoover-Chain thermostat
        allocate(pNHCTherm)
        if (input%ctrl%tInitNHC) then
          call init(pNHCTherm, pRandomGenerator, mass(indMovedAtom), &
              & pTempProfile, input%ctrl%wvScale, pMDFrame, input%ctrl%deltaT, &
              & input%ctrl%nh_npart, input%ctrl%nh_nys, input%ctrl%nh_nc, &
              & input%ctrl%xnose, input%ctrl%vnose, input%ctrl%gnose)
        else
          call init(pNHCTherm, pRandomGenerator, mass(indMovedAtom), pTempProfile,&
              & input%ctrl%wvScale, pMDFrame, input%ctrl%deltaT, &
              & input%ctrl%nh_npart, input%ctrl%nh_nys, input%ctrl%nh_nc)
        end if
        call init(pThermostat, pNHCTherm)
      end select

      !! Create MD integrator
      allocate(pVelocityVerlet)
      if (input%ctrl%tReadMDVelocities) then
        if (tBarostat) then
          call init(pVelocityVerlet, deltaT, coord0(:,indMovedAtom),&
              & pThermostat,input%ctrl%initialVelocities, &
              & BarostatStrength,pressure,input%ctrl%tIsotropic)
        else
          call init(pVelocityVerlet, deltaT, coord0(:,indMovedAtom),&
              & pThermostat,input%ctrl%initialVelocities)
        end if
      else
        if (tBarostat) then
          call init(pVelocityVerlet, deltaT, coord0(:,indMovedAtom),&
              & pThermostat, BarostatStrength,pressure,input%ctrl%tIsotropic)
        else
          call init(pVelocityVerlet, deltaT, coord0(:,indMovedAtom), pThermostat)
        end if
      end if
      allocate(pMDIntegrator)
      call init(pMDIntegrator, pVelocityVerlet)
    end if

    ! Check for extended Born-Oppenheimer MD
    tXlbomd = allocated(input%ctrl%xlbomd)
    if (tXlbomd) then
      if (input%ctrl%iThermostat /= 0) then
        call error("XLBOMD does not work with thermostats yet")
      elseif (tBarostat) then
        call error("XLBOMD does not work with barostats yet")
      elseif (nSpin /= 1 .or. tDFTBU) then
        call error("XLBOMD does not work for spin or DFTB+U yet")
      elseif (forceType /= 2 .and. forceType /= 3) then
        call error("Force evaluation method incompatible with XLBOMD")
      elseif (iDistribFn /= Fermi) then
        call error("Filling function incompatible with XLBOMD")
      end if
      allocate(xlbomdIntegrator)
      call Xlbomd_init(xlbomdIntegrator, input%ctrl%xlbomd, nIneqOrb)
    end if

    if (tDerivs) then
      allocate(tmp3Coords(3,nMovedAtom))
      tmp3Coords = coord0(:,indMovedAtom)
      call create(derivDriver, tmp3Coords, input%ctrl%deriv2ndDelta)
      coord0(:,indMovedAtom) = tmp3Coords
      deallocate(tmp3Coords)
      nGeoSteps = 2 * 3 * nMovedAtom - 1
    end if

    if (tEField) then
      EFieldStrength = input%ctrl%EFieldStrength
      EfieldVector(:) = input%ctrl%EfieldVector(:)
      tTDEfield = input%ctrl%tTDEfield
      EfieldOmega = input%ctrl%EfieldOmega
      EfieldPhase = input%ctrl%EfieldPhase
      if (tTDEfield .and. .not. tMD) then
        call error ("Time dependent electric fields only possible for MD!")
      end if
      ! parser should catch all of these:
      @:ASSERT(.not.tTDEfield .or. tMD)
    else
      tEField = .false.
      EFieldStrength = 0.0_dp
      EfieldVector(:) = 0.0_dp
      tTDEfield = .false.
      EfieldOmega = 0.0_dp
      EfieldPhase = 0
    end if

    !! Allocate charge arrays
    if (tMulliken) then ! automatically true if tSCC
      allocate(q0(orb%mOrb, nAtom, nSpin))
      q0(:,:,:) = 0.0_dp

      allocate(qShell0(orb%mShell, nAtom))
      qShell0(:,:) = 0.0_dp
    else
      allocate(q0(0,0,0))
      allocate(qShell0(0,0))
    end if

    allocate(qInput(orb%mOrb, nAtom, nSpin))
    allocate(qOutput(orb%mOrb, nAtom, nSpin))
    qInput(:,:,:) = 0.0_dp
    qOutput(:,:,:) = 0.0_dp

    if (tDFTBU) then
      allocate(qBlockIn(orb%mOrb, orb%mOrb, nAtom, nSpin))
      allocate(qBlockOut(orb%mOrb, orb%mOrb, nAtom, nSpin))
      qBlockIn = 0.0_dp
      qBlockOut = 0.0_dp
      if (tImHam) then
        allocate(qiBlockIn(orb%mOrb, orb%mOrb, nAtom, nSpin))
        qiBlockIn = 0.0_dp
      else
        allocate(qiBlockIn(0, 0, 0, 0))
        qiBlockIn = 0.0_dp
      end if
    else
      allocate(qBlockIn(0, 0, 0, 0))
      allocate(qBlockOut(0, 0, 0, 0))
      allocate(qiBlockIn(0, 0, 0, 0))
      qiBlockIn = 0.0_dp
      qBlockIn = 0.0_dp
      qBlockOut = 0.0_dp
    end if

    if (tImHam) then
      allocate(qiBlockOut(orb%mOrb, orb%mOrb, nAtom, nSpin))
      qiBlockOut = 0.0_dp
    end if

    if (tSCC) then
      allocate(qDiffRed(nMixElements))
      allocate(qInpRed(nMixElements))
      allocate(qOutRed(nMixElements))
      qDiffRed = 0.0_dp
      qInpRed = 0.0_dp
      qOutRed = 0.0_dp
    end if

    !! Initialize Mulliken charges
    if (tMulliken .or. tLinResp) then
      call initQFromShellChrg(q0, referenceN0, species0, orb)
    end if


    tReadChrg = input%ctrl%tReadChrg
    if (tSCC) then
      do iAt = 1, nAtom
        iSp = species0(iAt)
        do iSh = 1, orb%nShell(iSp)
          qShell0 (iSh,iAt) = sum(q0(orb%posShell(iSh,iSp): &
              & orb%posShell(iSh+1,iSp)-1,iAt,1))
        end do
      end do
      if (tReadChrg) then
        if (tDFTBU) then
          if (nSpin == 2) then
            if (tFixEf) then ! do not check charge or magnetisation from file
              call initQFromFile(qInput, fChargeIn, orb, qBlock=qBlockIn)
            else
              call initQFromFile(qInput, fChargeIn, orb, nEl = sum(nEl), &
                  & magnetisation=nEl(1)-nEl(2), qBlock=qBlockIn)
            end if
          else
            if (tImHam) then
              if (tFixEf) then
                call initQFromFile(qInput, fChargeIn, orb, &
                    & qBlock=qBlockIn,qiBlock=qiBlockIn)
              else
                call initQFromFile(qInput, fChargeIn, orb, nEl = nEl(1), &
                    & qBlock=qBlockIn,qiBlock=qiBlockIn)
              end if
            else
              if (tFixEf) then
                call initQFromFile(qInput, fChargeIn, orb, qBlock=qBlockIn)
              else
                call initQFromFile(qInput, fChargeIn, orb, nEl = nEl(1), &
                    & qBlock=qBlockIn)
              end if
            end if
          end if
        else
          ! hack again caused by going from up/down to q and M
          if (nSpin == 2) then
            if (tFixEf) then
              call initQFromFile(qInput, fChargeIn, orb)
            else
              call initQFromFile(qInput, fChargeIn, orb, nEl = sum(nEl),&
                  & magnetisation=nEl(1)-nEl(2))
            end if
          else
            if (tFixEf) then
              call initQFromFile(qInput, fChargeIn, orb)
            else
              call initQFromFile(qInput, fChargeIn, orb, nEl = nEl(1))
            end if
          end if
        end if
      else
        if (allocated(input%ctrl%initialCharges)) then
          if (abs(sum(input%ctrl%initialCharges) - input%ctrl%nrChrg) &
              &> 1e-4_dp) then
            write(strTmp, "(A,G13.6,A,G13.6,A,A)") "Sum of initial charges&
                & does not match specified total charge. (", &
                & sum(input%ctrl%initialCharges), " vs. ", &
                &input%ctrl%nrChrg, ") ", &
                & "Your initial charge distribution will be rescaled."
            call warning(strTmp)
          end if
          call initQFromAtomChrg(qInput, input%ctrl%initialCharges, &
              &referenceN0, species0, speciesName, orb)
        else
          qInput(:,:,:) = q0
        end if
        ! Rescaling to ensure correct number of electrons in the system
        qInput(:,:,1) = qInput(:,:,1) *  sum(nEl) / sum(qInput(:,:,1))


        select case (nSpin)
        case (1)
          ! nothing to do
        case (2)
          if (allocated(input%ctrl%initialSpins)) then
            do ii = 1, nAtom
              !! does not actually matter if additional spin polarization pushes
              !! charges to <0 as the initial charges are not mixed in to later
              !! iterations
              qInput(1:orb%nOrbAtom(ii),ii,2) = &
                  & qInput(1:orb%nOrbAtom(ii),ii,1) * &
                  & input%ctrl%initialSpins(1,ii) / &
                  & sum(qInput(1:orb%nOrbAtom(ii),ii,1))
            end do
          else
            do ii = 1, nAtom
              qInput(1:orb%nOrbAtom(ii),ii,2) = &
                  & qInput(1:orb%nOrbAtom(ii),ii,1) * &
                  & (nEl(1)-nEl(2))/sum(qInput(:,:,1))
            end do
          end if
        case (4)
          if (tSpin) then
            if (.not. allocated(input%ctrl%initialSpins)) then
              call error("Missing initial spins!")
            end if
            if (any(shape(input%ctrl%initialSpins)/=(/3,nAtom/))) then
              call error("Incorrect shape initialSpins array!")
            end if
            do ii = 1, nAtom
              do jj = 1, 3
                qInput(1:orb%nOrbAtom(ii),ii,jj+1) = &
                    & qInput(1:orb%nOrbAtom(ii),ii,1) * &
                    & input%ctrl%initialSpins(jj,ii) &
                    & / sum(qInput(1:orb%nOrbAtom(ii),ii,1))
              end do
            end do
          end if
        end select
        if (tDFTBU) then
          qBlockIn = 0.0_dp
          do iS = 1, nSpin
            do iAt = 1, nAtom
              iSp = species0(iAt)
              do iSh = 1, orb%nShell(iSp)
                iStart = orb%posShell(iSh,iSp)
                iEnd = orb%posShell(iSh+1,iSp)-1
                rTmp = sum(qInput(iStart:iEnd,iAt,iS))
                rTmp = rTmp / real(iEnd+1-iStart,dp)
                do ii = iStart, iEnd
                  qBlockIn(ii,ii,iAt,iS) = rTmp
                end do
              end do
            end do
          end do
          if (tImHam) then
            qiBlockIn = 0.0_dp
          end if
        end if
      end if

      qInpRed = 0.0_dp
      if (nSpin == 2) then
        call qm2ud(qInput)
        if (tDFTBU) then
          call qm2ud(qBlockIn)
        end if
      end if

      call OrbitalEquiv_reduce(qInput, iEqOrbitals, orb, qInpRed(1:nIneqOrb))
      if (tDFTBU) then
        call AppendBlock_reduce( qBlockIn,iEqBlockDFTBU, orb, &
            & qInpRed )
        if (tImHam) then
          call AppendBlock_reduce( qiBlockIn,iEqBlockDFTBULS, orb, &
              & qInpRed, skew=.true. )
        end if
      end if

      if (nSpin == 2) then
        call ud2qm(qInput)
        if (tDFTBU) then
          call ud2qm(qBlockIn)
        end if
      end if
    end if

    !! Initalize images (translations)
    if (tPeriodic) then
      call getCellTranslations(cellVec, rCellVec, latVec, recVec2p, mCutoff)
    else
      allocate(cellVec(3, 1))
      allocate(rCellVec(3, 1))
      cellVec(:,1) = [0.0_dp, 0.0_dp, 0.0_dp]
      rCellVec(:,1) = [0.0_dp, 0.0_dp, 0.0_dp]
    end if

    !! Initialize neighborlist.
    allocate(neighborList)
    call init(neighborList, nAtom, nInitNeighbor)
    allocate(nNeighbor(nAtom))


    !! Set various options
    tWriteTagged = input%ctrl%tWriteTagged
    tWriteDetailedXML = input%ctrl%tWriteDetailedXML
    tWriteResultsTag = input%ctrl%tWriteResultsTag
    tWriteDetailedOut = input%ctrl%tWriteDetailedOut
    tWriteBandDat = input%ctrl%tWriteBandDat
    tWriteHS = input%ctrl%tWriteHS
    tWriteRealHS = input%ctrl%tWriteRealHS

    !! Minimize memory usage?
    tMinMemory = input%ctrl%tMinMemory
    tStoreEigvecs = tMinMemory .and. (nKPoint > 1 .or. nSpin == 2 )
    if (tStoreEigvecs) then
      if (tRealHS.and.(.not.t2Component)) then
        allocate(storeEigvecsReal(nSpin))
        allocate(storeEigvecsCplx(0))
        do iS = 1, nSpin
          call init(storeEigvecsReal(iS), 0, "tmp_eigvr_")
        end do
      else
        if (t2Component) then
          allocate(storeEigvecsCplx(1))
          allocate(storeEigvecsReal(0))
          call init(storeEigvecsCplx(1), 0, "tmp_eigvc_")
        else
          allocate(storeEigvecsCplx(nSpin))
          allocate(storeEigvecsReal(0))
          do iS = 1, nSpin
            call init(storeEigvecsCplx(iS), 0, "tmp_eigvc_")
          end do
        end if
      end if
    else
      allocate(storeEigvecsReal(0))
      allocate(storeEigvecsCplx(0))
    end if

    !! Check if stopfiles already exist and quit if yes
    inquire(file=fStopSCC, exist=tExist)
    if (tExist) then
      call error("Stop file '" // fStopSCC // "' already present at startup")
    end if
    inquire(file=fStopDriver, exist=tExist)
    if (tExist) then
      call error("Stop file '" // fStopDriver // "' already present at startup")
    end if

    restartFreq = input%ctrl%restartFreq

    tInitialized = .true.

    if (input%ctrl%tMD) then
      select case(input%ctrl%iThermostat)
      case (0)
        if (tBarostat) then
          write(*,"('Mode:',T30,A,/,T30,A)") 'MD without scaling of &
              &velocities', '(a.k.a. "NPE" ensemble)'
        else
          write(*,"('Mode:',T30,A,/,T30,A)") 'MD without scaling of &
              &velocities', '(a.k.a. NVE ensemble)'
        end if
      case (1)
        if (tBarostat) then
          write(*,"('Mode:',T30,A,/,T30,A)") &
              & "MD with re-selection of velocities according to temperature", &
              & "(a.k.a. NPT ensemble using Andersen thermostating + Berensen&
              & barostat)"
        else
          write(*,"('Mode:',T30,A,/,T30,A)") &
              & "MD with re-selection of velocities according to temperature", &
              & "(a.k.a. NVT ensemble using Andersen thermostating)"
        end if
      case(2)
        if (tBarostat) then
          write(*,"('Mode:',T30,A,/,T30,A)") &
              & "MD with scaling of velocities according to temperature", &
              & "(a.k.a. 'not' NVP ensemble using Berendsen thermostating and&
              & barostat)"
        else
          write(*,"('Mode:',T30,A,/,T30,A)") &
              & "MD with scaling of velocities according to temperature", &
              & "(a.k.a. 'not' NVT ensemble using Berendsen thermostating)"
        end if
      case(3)
        if (tBarostat) then
          write(*,"('Mode:',T30,A,/,T30,A)") &
              & "MD with scaling of velocities according to", &
              & "Nose-Hoover-Chain thermostat + Berensen barostat"
        else
          write(*,"('Mode:',T30,A,/,T30,A)") &
              & "MD with scaling of velocities according to", &
              & "Nose-Hoover-Chain thermostat"
        end if

      case default
        call error("Unknown thermostat mode")
      end select
    elseif (tGeoOpt) then
      if (nGeoConstr > 0) then
        strTmp = "with constraints"
      else
        strTmp = ""
      end if
      select case (input%ctrl%iGeoOpt)
      case (1)
        write(*,"('Mode:',T30,A)")'Steepest descent' // trim(strTmp)
      case (2)
        write(*,"('Mode:',T30,A)") 'Conjugate gradient relaxation' &
            &// trim(strTmp)
      case (3)
        write(*,"('Mode:',T30,A)") 'Modified gDIIS relaxation' &
            &// trim(strTmp)
      case default
        call error("Unknown optimisation mode")
      end select
    elseif (tDerivs) then
      write (*,"('Mode:',T30,A)") "2nd derivatives calculation"
      write (*,"('Mode:',T30,A)") "Calculated for atoms:"
      write (*,*) indMovedAtom
    elseif (tSocket) then
      write (*,"('Mode:',T30,A)") "Socket controlled calculation"
    else
      write (*,"('Mode:',T30,A)") "Static calculation"
    end if

    if (tSCC) then
      write (*, "(A,':',T30,A)") "Self consistent charges", "Yes"
      write (*, "(A,':',T30,E14.6)") "SCC-tolerance", sccTol
      write (*, "(A,':',T30,I14)") "Max. scc iterations", nSCCIter
      write (*, "(A,':',T30,E14.6)") "Ewald alpha parameter", getSCCEwaldPar()
      if (tDFTBU) then
        write (*, "(A,':',T35,A)") "Orbitally dependant functional", "Yes"
        write (*, "(A,':',T30,I14)") "Orbital functional number",nDFTBUfunc !
        !  use module to reverse look up name
      end if
    else
      write (*, "(A,':',T30,A)") "Self consistent charges", "No"
    end if

    select case (nSpin)
    case(1)
      write (*,"(A,':',T30,A)") "Spin polarisation", "No"
      write (*, "(A,':',T30,F12.6,/,A,':',T30,F12.6)") "Nr. of up electrons", &
          &0.5_dp*nEl(1), "Nr. of down electrons", 0.5_dp*nEl(1)
    case(2)
      write (*,"(A,':',T30,A)") "Spin polarisation", "Yes"
      write (*, "(A,':',T30,F12.6,/,A,':',T30,F12.6)") "Nr. of up electrons", &
          &nEl(1), "Nr. of down electrons", nEl(2)
    case(4)
      write (*,"(A,':',T30,A)") "Non-collinear calculation", "Yes"
      write (*, "(A,':',T30,F12.6)") "Nr. of electrons", nEl(1)
    end select

    if (tPeriodic) then
      write (*, "(A,':',T30,A)") "Periodic boundaries", "Yes"
      if (tLatOpt) then
        write (*, "(A,':',T30,A)") "Lattice optimisation", "Yes"
        write (*, "(A,':',T30,f12.6)") "Pressure", pressure
      end if
    else
      write (*, "(A,':',T30,A)") "Periodic boundaries", "No"
    end if

    select case (solver)
    case(1)
      write (strTmp, "(A)") "Standard"
    case(2)
      write (strTmp, "(A)") "Divide and Conquer"
    case(3)
      write (strTmp, "(A)") "Relatively robust (version 1)"
    case(4)
      write (strTmp, "(A)") "Relativel robust (version 2)"
    case default
      call error("Unknown eigensolver!")
    end select
    write (*, "(A,':',T30,A)") "Diagonalizer", trim(strTmp)

    if (tSCC) then
      select case (iMixer)
      case(1)
        write (strTmp, "(A)") "Simple"
      case(2)
        write (strTmp, "(A)") "Anderson"
      case(3)
        write (strTmp, "(A)") "Broyden"
      case(4)
        write (strTmp, "(A)") "DIIS"
      end select
      write (*, "(A,':',T30,A,' ',A)") "Mixer", trim(strTmp), "mixer"
      write (*, "(A,':',T30,F14.6)") "Mixing parameter", mixParam
      write (*, "(A,':',T30,I14)") "Maximal SCC-cycles", nSCCIter
      select case (iMixer)
      case(2)
        write (*, "(A,':',T30,I14)") "Nr. of chrg. vectors to mix", nGeneration
      case(3)
        write (*, "(A,':',T30,I14)") "Nr. of chrg. vec. in memory", nGeneration
      case(4)
        write (*, "(A,':',T30,I14)") "Nr. of chrg. vectors to mix", nGeneration
      end select
    end if

    if (tCoordOpt) then
      write (*, "(A,':',T30,I14)") "Nr. of moved atoms", nMovedAtom
    end if
    if (tGeoOpt) then
      write (*, "(A,':',T30,I14)") "Max. nr. of geometry steps", nGeoSteps
      write (*, "(A,':',T30,E14.6)") "Force tolerance", input%ctrl%maxForce
      if (input%ctrl%iGeoOpt == 1) then
        write (*, "(A,':',T30,E14.6)") "Step size", deltaT
      end if
    end if

    if (tForces) then
      select case (forceType)
      case(0)
        strTmp = "Traditional"
      case(2)
        strTmp = "Dynamics, zero electronic temp."
      case(3)
        strTmp = "Dynamics, finite electronic temp."
      end select
      write(*, "(A,':',T30,A)") "Force evaluation method", strTmp
    end if

    tFirst = .true.
    if (nGeoConstr > 0) then
      do ii = 1, nAtom
        do jj = 1, nGeoConstr
          if (conAtom(jj) == ii) then
            if (tFirst) then
              write(strTmp, "(A,':')") "Geometry constraints"
              tFirst = .false.
            else
              write(strTmp, "(A)") ""
            end if
            write(*, "(A,T30,'At',I4,': ',3F10.6)") trim(strTmp), &
                &ii, (conVec(kk,jj), kk=1,3)
          end if
        end do
      end do
    end if

    if (.not.input%ctrl%tSetFillingTemp) then
      write (*, "(A,':',T30,E14.6)") "Electronic temperature", tempElec
    end if
    if (tMD) then
      write (*, "(A,':',T30,E14.6)") "Time step", deltaT
      if (input%ctrl%iThermostat == 0 .and. .not.input%ctrl%tReadMDVelocities)&
          & then
        write (*, "(A,':',T30,E14.6)") "Temperature", tempAtom
      end if
      write (*, "(A,':',T30,I14)") "Random seed", iSeed
      if (input%ctrl%tRescale) then
        write (*, "(A,':',T30,E14.6)") "Rescaling probability", &
            &input%ctrl%wvScale
      end if
    end if

    if (tSCC) then
      if (input%ctrl%tReadChrg) then
        write (strTmp, "(A,A,A)") "Read in from '", trim(fChargeIn), "'"
      else
        write (strTmp, "(A,E11.3,A)") "Set automatically (system chrg: ", &
            &input%ctrl%nrChrg, ")"
      end if
      write (*, "(A,':',T30,A)") "Initial charges", trim(strTmp)
    end if

    do iSp = 1, nType
      if (iSp == 1) then
        write (strTmp, "(A,':')") "Included shells"
      else
        write (strTmp, "(A)") ""
      end if
      do jj = 1, orb%nShell(iSp)
        if (jj == 1) then
          strTmp2 = trim(orbitalNames(orb%angShell(jj, iSp) + 1))
        else
          strTmp2 = trim(strTmp2) // ", " &
              &// trim(orbitalNames(orb%angShell(jj, iSp) + 1))
        end if
      end do
      write (*, "(A,T29,A2,':  ',A)") trim(strTmp), trim(speciesName(iSp)), &
          &trim(strTmp2)
    end do

    if (tPeriodic) then
      do ii = 1, nKPoint
        if (ii == 1) then
          write(strTmp, "(A,':')") "K-points and weights"
        else
          write(strTmp, "(A)") ""
        end if
        write(*,"(A,T28,I6,':',3F10.6,3X,F10.6)") trim(strTmp), ii, &
            & (kPoint(jj, ii), jj=1, 3), kWeight(ii)
      end do
    end if

    if (tDispersion) then
      select type (dispersion)
      type is (DispSlaKirk)
        write(*,"(A)") "Using Slater-Kirkwood dispersion corrections"
      type is (DispUff)
        write(*,"(A)") "Using Lennard-Jones dispersion corrections"
    #:if WITH_DFTD3
      type is (DispDftD3)
        write(*,"(A)") "Using DFT-D3 dispersion corrections"
    #:endif
      class default
        call error("Unknown dispersion model - this should not happen!")
      end select
    end if

    if (tSCC) then
      ! Have the SK values of U been replaced?
      if (allocated(input%ctrl%hubbU)) then
        strTmp = ""
        tFirst = .true.
        do iSp = 1, nType
          if (all(input%ctrl%hubbU(1:orb%nShell(iSp), iSp) == 0.0_dp)) then
            cycle
          end if
          do jj = 1, orb%nShell(iSp)
            if (tFirst) then
              write(strTmp, "(A,':')") "Non-default Hubbard U"
              tFirst = .false.
            else
              write(strTmp, "(A)") ""
            end if
            write(*, "(A,T30,A2,2X,I1,'(',A1,'): ',E14.6)") trim(strTmp), speciesName(iSp), jj, &
                & orbitalNames(orb%angShell(jj, iSp)+1), hubbU(jj, iSp)
          end do
        end do
      end if
    end if

    tFirst = .true.
    if (tSpin) then
      do iSp = 1, nType
        do jj = 1, orb%nShell(iSp)
          do kk = 1, orb%nShell(iSp)
            if (tFirst) then
              write(strTmp, "(A,':')") "Spin coupling constants"
              tFirst = .false.
            else
              write(strTmp, "(A)") ""
            end if
            write(*, "(A,T30,A2,2X,I1,'(',A1,')-',I1,'(',A1,'): ',E14.6)") &
                &trim(strTmp), speciesName(iSp), &
                &jj, orbitalNames(orb%angShell(jj, iSp)+1), &
                &kk, orbitalNames(orb%angShell(kk, iSp)+1), &
                &spinW(kk, jj, iSp)
          end do
        end do
      end do
    end if

    tFirst = .true.
    if (tSpinOrbit) then
      if (tDualSpinOrbit) then
        write(*,"(A)")"Dual representation spin orbit"
      end if
      do iSp = 1, nType
        do jj = 1, orb%nShell(iSp)
          if (tFirst) then
            write(strTmp, "(A,':')") "Spin orbit constants"
            tFirst = .false.
          else
            write(strTmp, "(A)") ""
          end if
          write(*, "(A,T30,A2,2X,I1,'(',A1,'): ',E14.6)") &
                &trim(strTmp), speciesName(iSp), &
                &jj, orbitalNames(orb%angShell(jj, iSp)+1), &
                &xi(jj, iSp)
          if (xi(jj, iSp) /= 0.0_dp .and. orb%angShell(jj, iSp) == 0) then
            call error("Program halt due to non-zero s-orbital spin-orbit &
                &coupling constant!")
          end if
        end do
      end do
    end if

    if (tSCC) then
      if (t3rdFull) then
        write (*, "(A,T30,A)") "Full 3rd order correction", "Yes"
        if (input%ctrl%tOrbResolved) then
          write (*, "(A,T30,A)") "Orbital-resolved 3rd order", "Yes"
          write(*, "(A30)") "Shell-resolved Hubbard derivs:"
          write(*, "(A)") "        s-shell   p-shell   d-shell   f-shell"
          do iSp = 1, nType
            write(*,"(A3,A3,4F10.4)") "  ", trim(speciesName(iSp)),&
                & input%ctrl%hubDerivs(:orb%nShell(iSp),iSp)
          end do
        end if
      end if

      if (any(tDampedShort)) then
        write(*, "(A,T30,A)") "Damped SCC", "Yes"
        ii = count(tDampedShort)
        write(strTmp, "(A,I0,A)") "(A,T30,", ii, "(A,1X))"
        write(*, strTmp) "Damped species(s):", pack(speciesName, tDampedShort)
        deallocate(tDampedShort)
      end if
    end if

    write (*, "(A,':')") "Extra options"
    if (tPrintMulliken) then
      write (*, "(T30,A)") "Mulliken analysis"
    end if
    if (tPrintForces .and. .not. (tMD .or. tGeoOpt .or. tDerivs)) then
      write (*, "(T30,A)") "Force calculation"
    end if
    if (tPrintEigVecs) then
      write (*, "(T30,A)") "Eigenvector printing"
    end if
    if (tExtChrg) then
      write (*, "(T30,A)") "External charges specified"
    end if

    if (tEField) then
      if (tTDEfield) then
        write (*, "(T30,A)") "External electric field specified"
        write (*, "(A,':',T30,E14.6)") "Angular frequency", EfieldOmega
      else
        write (*, "(T30,A)") "External static electric field specified"
      end if
      write (*, "(A,':',T30,E14.6)") "Field strength", EFieldStrength
      write (*, "(A,':',T30,3F9.6)") "Direction", EfieldVector
      if (tPeriodic) then
        call warning("Saw tooth potential used for periodic geometry &
            &- make sure there is a vacuum region!")
      end if
    end if

    if (tDFTBU) then
      do iSp = 1, nType
        if (nUJ(iSp)>0) then
          write(strTmp, "(A,':')") "U-J coupling constants"
          write(*, "(A,T25,A2)")trim(strTmp), speciesName(iSp)
          do jj = 1, nUJ(iSp)
            write(strTmp, "(A,I1,A)")'(A,',niUJ(jj,iSp),'I2,T25,A,F6.4)'
            write(*,trim(strTmp))'Shells:',iUJ(1:niUJ(jj,iSp),jj,iSp),'UJ:', &
                & UJ(jj,iSp)
          end do
        end if
      end do

    end if

    if (tSpinOrbit) then
      if (tDualSpinOrbit) then
        if ( (tEField .or. tExtChrg) .and. tForces) then
          call error("Currently there is a force bug for dual spin orbit &
              &coupling")
        end if
      end if
    end if

    select case (forceType)
    case(0)
      write(*, "(A,T30,A)") "Force type", "original"
    case(1)
      write(*, "(A,T30,A)") "Force type", "erho with re-diagonalized &
          & eigenvalues"
    case(2)
      write(*, "(A,T30,A)") "Force type", "erho with DHD-product (T_elec = 0K)"
    case(3)
      write(*, "(A,T30,A)") "Force type", "erho with S^-1 H D (Te <> 0K)"
    end select

    if ((tSpinOrbit .and. tDFTBU) .and. tForces)  then
      call error("Currently there is a force bug for dual DFTB+U with spin &
          &orbit coupling")
    end if

    if (.not.tStress) then
      if (tBarostat) then
        call error("Sorry, MD with a barostat requires stress evaluation")
      end if
      if (tLatOpt) then
        call error("Sorry, lattice optimization requires stress tensor&
            & evaluation")
      end if
    end if

    if (tSpinOrbit .and. (tWriteHS .or.(tWriteRealHS.and..not.tDualSpinOrbit)))&
        & then
      call error("Writing of Hamiltonian currently not possible with spin orbit&
          & coupling enabled.")
    end if

    if (tLinResp) then
      if (tMinMemory) then
        call error("Linear response is not compatible with MinimiseMemoryUsage&
            & yet")
      end if

      if (tDFTBU) then
        call error("Linear response is not compatible with Orbitally dependant&
            & functionals yet")
      end if

      if (tForces .and. nSpin > 1) then
        call error("Linear response is not available for spin polarised forces&
            & yet")
      end if

      if (t2Component) then
        call error("Linear response is not compatibile with this spinor&
            & Hamiltonian")
      end if

      if (tStress) then
        call error("Excited state stresses not implemented")
      end if

      if (.not.tRealHS) then
        call error("Linear response does not support k-points")
      end if

    end if

    tInitialized = .true.

  end subroutine initProgramVariables



  subroutine destructProgramVariables()

    integer :: ii

    if (tProjEigenvecs) then
      call destruct(iOrbRegion)
      call destruct(RegionLabels)
    end if
    do ii = 1, size(storeEigvecsReal)
      call destruct(storeEigvecsReal(ii))
    end do
    do ii = 1, size(storeEigvecsCplx)
      call destruct(storeEigvecsCplx(ii))
    end do

  end subroutine destructProgramVariables


end module initprogram
