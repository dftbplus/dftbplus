!--------------------------------------------------------------------------------------------------!
!  DFTB+: general package for performing fast atomistic simulations                                !
!  Copyright (C) 2017  DFTB+ developers group                                                      !
!                                                                                                  !
!  See the LICENSE file for terms of usage and distribution.                                       !
!--------------------------------------------------------------------------------------------------!

!!* The main dftb+ program
program dftbplus
#include "allocate.h"
#include "assert.h"
  use constants
  use initprogram
  use inputdata_module
  use nonscc
  use eigenvects
  use repulsive
  use etemp
  use populations
  use densitymatrix
  use forces
  use stress
  use lapackroutines, only : matinv ! reguired to calculate lattice derivs
  ! from stress tensors
  use simplealgebra, only : determinant33 ! required for cell volumes
  use taggedoutput
  use scc
  use externalcharges
  use mixer
  use geoopt
  use numderivs2
  use spin
  use dftbplusu
  use fileid
  use formatout
  use mdcommon
  use energies
  use potentials
  use orbitalequiv
  use parser
  use sparse2dense, only : unpackHS, packHS, packERho, iPackHS,&
      & blockSymmetrizeHS, blockHermitianHS
  use blasroutines, only : hemv, gemm, symm, hemm
  use hsdutils
  use charmanip
  use shift
  use spinorbit
  use angmomentum
  use elecconstraints
  use pmlocalisation
  use linresp_module
!  use emfields
  use mainio
  use xmlf90
  implicit none

  !! Revision control strings
  character(len=*), parameter :: RELEASE_VERSION = '17.1'
  integer, parameter :: RELEASE_YEAR = 2017

  type(inputData)          :: input             ! Contains the parsed input

  integer                  :: nk, iEgy, nSpin2, nK2, iSpin2, iK2
  complex(dp), allocatable :: HSqrCplx(:,:,:,:), SSqrCplx(:,:), HSqrCplx2(:,:)
  real(dp),    allocatable :: HSqrReal(:,:,:), SSqrReal(:,:), HSqrReal2(:,:)
  real(dp),    allocatable :: eigen(:,:,:), eigen2(:,:,:)
  real(dp), allocatable    :: rhoPrim(:,:)
  real(dp), allocatable    :: iRhoPrim(:,:)
  real(dp), allocatable    :: ERhoPrim(:), ERhoPrim2(:)
  real(dp), allocatable    :: h0(:)

  ! variables for derivatives using the Hellmann-Feynman theorem:
  real(dp), allocatable    :: hprime(:,:) ! for derivatives of H wrt external
  real(dp), allocatable    :: potentialDerivative(:,:) ! for derivatives of V
  real(dp), allocatable    :: dipoleTmp(:,:) ! temporary dipole data

  real(dp), allocatable    :: filling(:,:,:)
  real(dp), allocatable :: Eband(:), TS(:), E0(:), Eold

  type(TEnergies)          :: energy
  type(TPotentials)        :: potential


  real(dp), allocatable    :: derivs(:,:),repulsiveDerivs(:,:),totalDeriv(:,:)
  real(dp), allocatable    :: chrgForces(:,:)
  real(dp), allocatable    :: excitedDerivs(:,:) ! excited state force addition

  ! Stress tensors for various contribution in periodic calculations
  real(dp) :: elecStress(3,3), repulsiveStress(3,3), kineticStress(3,3)
  real(dp) :: dispStress(3,3), totalStress(3,3)

  ! Derivatives of lattice vectors in periodic calculations
  real(dp) :: elecLatDeriv(3,3), repulsiveLatDeriv(3,3)
  real(dp) :: dispLatDeriv(3,3), totalLatDeriv(3,3)
  ! derivative of cell volume wrt to lattice vectors, needed for pV term
  real(dp) :: derivCellVol(3,3)
  
  real(dp) :: dipoleMoment(3)
  real(dp) :: angularMomentum(3) ! hold total angular momentum vector

  integer                  :: ii, jj, kk

  logical                  :: tConverged

  logical, parameter       :: tDensON2 = .false.  ! O(N^2) density mtx creation
  logical, parameter       :: tAppendDetailedOut = .false.

  character(len=*), parameter :: formatEnergy = '(8f12.5)'
  character(len=*), parameter :: formatEigen = '(8f14.8)'
  character(len=*), parameter :: formatHessian = '(4f16.10)'
  character(len=*), parameter :: formatGeoOut = "(I5,F16.8,F16.8,F16.8)"
  ! formats for data with 1 or two units, and exponential notation form:
  character(len=*), parameter :: format1U = "(' ',A,':',T32,F18.10,T51,A)"
  character(len=*), parameter :: format2U = &
      &"(' ',A,':',T32,F18.10,T51,A,T54,F16.4,T71,A)"
  character(len=*), parameter :: format1Ue = "(' ',A,':',T37,E13.6,T51,A)"
  character(len=*), parameter :: format2Ue = &
      &"(' ',A,':',T37,E13.6,T51,A,T57,E13.6,T71,A)"
  character(len=*), parameter :: format1U1e = &
      &"(' ',A,':',T32,F18.10,T51,A,T57,E13.6,T71,A)"
  real(dp) :: cellPressure

  !! Variables for the geometry optimization
  integer :: iGeoStep                      !* Geometry steps so far
  integer :: iLatGeoStep                   !* Lattice geometry steps so far
  logical :: tGeomEnd                      !* Do we have the final geometry?
  logical :: tCoordEnd                     !* Has this completed?
  logical :: tCoordStep                    !* do we take an optimization step
  !* on the lattice or the internal coordinates if optimizing both in a
  !* periodic geometry
  real(dp) :: invLatVec(3,3)
  real(dp), allocatable, target :: coord0Fold(:,:) !* Folded coords (3, nAtom)
  real(dp), pointer :: pCoord0Out(:,:)  ! Coordinates to print out
  real(dp), allocatable :: new3Coord(:,:)     !* New coordinates returned by
  !* the MD routines
  real(dp) :: tmpLatVecs(9), newLatVecs(9) !* lattice vectors returned by
  ! the optimizer
  real(dp) :: tmpLat3Vecs(3,3)
  real(dp), allocatable :: velocities(:,:) !* MD velocities
  real(dp), allocatable :: movedVelo(:,:)  !* MD velocities for moved atoms
  real(dp), allocatable :: movedAccel(:,:) !* MD acceleration for moved atoms
  real(dp), allocatable :: movedMass(:,:)  !* Mass of the moved atoms
  real(dp) :: KE                           !* MD Kinetic energy
  real(dp) :: kT                           !* MD instantaneous thermal energy

  real(dp) :: Efield(3), absEfield !* external electric field

  real(dp) :: diffGeo                      !* Difference between last calculated
  !* and new geometry.

  !!* Loop variables
  integer :: iSCCIter, iSpin, iAtom, iNeigh

  integer :: fdTagged  !!* File descriptor for the tagged writer
  integer :: fdUser    !!* File descriptor for the human readable output
  integer :: fdBand    !!* File descriptor for the band structure output
  integer :: fdEigvec  !!* File descriptor for the eigenvector output
  integer :: fdResultsTag !!* File descriptor for detailed.tag
  integer :: fdMD      !!* File descriptor for extra MD output
  integer :: fdHessian !!* File descriptor for numerical Hessian

  !!* Name of the human readable file
  character(*), parameter :: taggedOut = "autotest.tag"
  character(*), parameter :: userOut = "detailed.out"
  character(*), parameter :: bandOut = "band.out"
  character(*), parameter :: mdOut = "md.out"
  character(*), parameter :: resultsTag = "results.tag"
  character(*), parameter :: hessianOut = "hessian.out"


  real(dp) :: sccErrorQ           !* Charge error in the last iterations
  real(dp) :: rTmp
  real(dp), allocatable :: tmpDerivs(:)
  real(dp), allocatable :: tmpMatrix(:,:)
  real(dp), allocatable :: orbitalL(:,:,:), orbitalLPart(:,:,:)
  real(dp), allocatable    :: rVecTemp(:)
  character(lc) :: lcTmp

  character(lc) :: tmpStr  !!* temporary character variable

  real(dp), pointer :: pDynMatrix(:,:)

  logical :: tWriteRestart = .false. !* flag to write out geometries (and
  !* charge data if scc) when moving atoms about - in the case of conjugate
  !* gradient/steepest descent the geometries are written anyway
  integer :: minSCCIter                   !* Minimal number of SCC iterations

  type(xmlf_t) :: xf
  real(dp), allocatable :: bufferRealR2(:,:)
  logical :: tStopSCC, tStopDriver   ! if scf/geometry driver should be stopped
  integer :: ang, iSh1, iSp1

  real(dp), allocatable :: shift3rd(:)
  real(dp), allocatable :: orbresshift3rd(:,:)
  
  real(dp), allocatable :: dqAtom(:) ! net charge on each atom  
  
  real(dp), allocatable :: rhoSqrReal(:,:,:) ! density matrix

  ! Natural orbitals for excited state density matrix, if requested
  real(dp), allocatable :: naturalOrbs(:,:,:), occNatural(:,:)
  
  real(dp), allocatable :: invJacobian(:,:)

  real(dp) :: localisation ! locality measure for the wavefunction

  integer :: nFilledLev ! temporary variable for number of occupied levels

  integer :: nSpinHams  ! Nr. of different spin Hamiltonians
  integer :: sqrHamSize  ! Size of the sqr Hamiltonian

  call printDFTBHeader(RELEASE_VERSION, RELEASE_YEAR)
  write (*,'(/A/)') "***  Parsing and initializing"

  !! Parse input and set the variables in the local scope according the input.
  !! These variables are defined in the initprogram module.

  call parseHSDInput(input)
  write (*,"(/A)") "Starting initialization..."
  write (*,"(A80)") repeat("-", 80)
  call initProgramVariables(input)
  call destroy(input)
  write (*,*)

  elecStress = 0.0_dp
  repulsiveStress = 0.0_dp
  kineticStress = 0.0_dp
  dispStress = 0.0_dp
  totalStress = 0.0_dp
  elecLatDeriv = 0.0_dp
  repulsiveLatDeriv = 0.0_dp
  dispLatDeriv = 0.0_dp
  totalLatDeriv = 0.0_dp
  derivCellVol = 0.0_dp
  if (tWriteTagged.or.tWriteResultsTag) then
    !! Initialize tagged writer and human readable output
    call initTaggedWriter()
  end if
  if (tWriteTagged) then
    fdTagged = getFileId()
    open(fdTagged, file=taggedOut, position="rewind", status="replace")
    close(fdTagged)
  end if
  if (tWriteResultsTag) then
    fdResultsTag = getFileId()
  end if

  if (tWriteBandDat) then
    fdBand = getFileId()
  end if
  fdEigvec = getFileId()

  if (tDerivs) then
    fdHessian = getFileId()
    open(fdHessian, file=hessianOut, position="rewind", status="replace")
  end if

  ! initially empty file
  if (tWriteDetailedOut) then
    fdUser = getFileId()
    open(fdUser, file=userOut, position="rewind", status="replace")
  end if

  ! initially open to file to be empty
  if (tMD) then
    fdMD = getFileId()
    open(fdMD, file=mdOut, position="rewind", status="replace")
  end if

  ALLOCATE_(rhoPrim, (0, nSpin))
  ALLOCATE_(h0, (0))
  ALLOCATE_(iRhoPrim, (0, nSpin))

  ALLOCATE_(excitedDerivs, (0, 0))

  if (tForces) then
    ALLOCATE_(ERhoPrim, (0))
    ALLOCATE_(ERhoPrim2, (0))
  end if

  if (tForces) then
    ALLOCATE_(derivs,(3,nAtom))
    ALLOCATE_(repulsiveDerivs,(3,nAtom))
    ALLOCATE_(totalDeriv, (3,nAtom))
    if (tExtChrg) then
      ALLOCATE_(chrgForces, (3, nExtChrg))
    end if
    if (tLinResp) then
      DEALLOCATE_(excitedDerivs)
      ALLOCATE_(excitedDerivs, (3, nAtom))
    end if
  end if
  
  call create(energy,nAtom)

  call create(potential,orb,nAtom,nSpin)
  ALLOCATE_(shift3rd, (nAtom))
  ALLOCATE_(orbresshift3rd, (orb%mShell,nAtom))
  
  ! Nr. of independent spin Hamiltonians
  select case (nSpin)
  case (1)
    nSpinHams = 1
    sqrHamSize = nOrb
  case (2)
    nSpinHams = 2
    sqrHamSize = nOrb
  case (4)
    nSpinHams = 1
    sqrHamSize = 2 * nOrb
  end select

  ALLOCATE_(TS, (nSpinHams))
  ALLOCATE_(E0, (nSpinHams))
  ALLOCATE_(Eband, (nSpinHams))
  ALLOCATE_(eigen, (sqrHamSize, nKPoint, nSpinHams))
  ALLOCATE_(eigen2, (sqrHamSize, nKPoint, nSpinHams))
  ALLOCATE_(filling,(sqrHamSize, nKpoint, nSpinHams))

  ALLOCATE_(coord0Fold, (3, nAtom))
  if (tShowFoldedCoord) then
    pCoord0Out => coord0Fold
  else
    pCoord0Out => coord0
  end if


  if (tMD.or.tDerivs) then
    ALLOCATE_(new3Coord, (3, nMovedAtom))
  end if

  if (tCoordOpt) then
    ALLOCATE_(tmpDerivs,(size(tmpCoords)))
  else
    ALLOCATE_(tmpDerivs,(0))
  end if

  if ((tMulliken .and. tSpinOrbit) .or. tImHam) then
    ALLOCATE_(orbitalL,(3,orb%mShell,nAtom))
    orbitalL = 0.0_dp
  else
    ALLOCATE_(orbitalL,(0,0,0))
  end if

  if ((tMulliken .and. tSpinOrbit) .and. .not.  tDualSpinOrbit) then
    ALLOCATE_(orbitalLPart,(3,orb%mShell,nAtom))
    orbitalLPart = 0.0_dp
  else
    ALLOCATE_(orbitalLPart,(0,0,0))
  end if
  eigen(:,:,:) = 0.0_dp
  eigen2(:,:,:) = 0.0_dp

  if (tStoreEigvecs) then
    nSpin2 = 1
    nK2 = 1
  else
    nSpin2 = nSpin
    nK2 = nKPoint
  end if

  ! If only H/S should be printed, no allocation for square HS is needed
  if (.not. (tWriteRealHS .or. tWriteHS)) then
    if (t2Component) then
      ALLOCATE_(HSqrCplx, (sqrHamSize, sqrHamSize, nK2, 1))
      ALLOCATE_(SSqrCplx, (sqrHamSize, sqrHamSize))
    elseif (tRealHS) then
      ALLOCATE_(HSqrReal, (sqrHamSize, sqrHamSize, nSpin2))
      if (any(forceType == [ 1, 2, 3 ])) then
        ALLOCATE_(HSqrReal2, (sqrHamSize, sqrHamSize))
      end if
      ALLOCATE_(SSqrReal, (sqrHamSize, sqrHamSize))
    else
      ALLOCATE_(HSqrCplx, (sqrHamSize, sqrHamSize, nK2, nSpin2))
      if (any(forceType == [ 1, 2, 3 ])) then
        ALLOCATE_(HSqrCplx2, (sqrHamSize, sqrHamSize))
      end if
      ALLOCATE_(SSqrCplx, (sqrHamSize, sqrHamSize))
    end if
  end if

  ALLOCATE_(rhoSqrReal, (0,0,0))
  ALLOCATE_(dqAtom, (0))
  if (tLinResp) then
    DEALLOCATE_(dqAtom)
    ALLOCATE_(dqAtom, (nAtom))
    if (tLinRespZVect) then
      DEALLOCATE_(rhoSqrReal)
      ALLOCATE_(rhoSqrReal, (sqrHamSize, sqrHamSize, nSpin))
    end if    
  end if
  
  if (tLinResp .and. tPrintExcitedEigVecs) then
    ALLOCATE(naturalOrbs(nOrb,nOrb,1))
    ALLOCATE(occNatural(nOrb,1))
  else
    ALLOCATE(naturalOrbs(0,0,0))
    ALLOCATE(occNatural(0,0))
  end if
  naturalOrbs = 0.0_dp
  occNatural = 0.0_dp
  
  if (tMD) then
    ALLOCATE_(velocities,(3,nAtom))
    ALLOCATE_(movedVelo, (3, nMovedAtom))
    ALLOCATE_(movedAccel, (3, nMovedAtom))
    ALLOCATE_(movedMass, (3, nMovedAtom))
    movedMass(:,:) = spread(mass(indMovedAtom),1,3)
    velocities(:,:) = 0.0_dp
  end if


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !! Geometry loop
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  tGeomEnd = nGeoSteps == 0

  tCoordStep = .false.
  if (tCoordOpt) then
    tCoordStep = .true.
    tCoordEnd = .false.
  end if

  iGeoStep = 0
  iLatGeoStep = 0
  tStopDriver = .false.

  if (tSCC) then
    if (tDFTBU) then
      minSCCIter = 2
    else
      if (nSpin == 1) then
        minSCCIter = 1
      else
        minSCCIter = 2
      end if
    end if
  else
    minSCCIter = 1
  end if

  if (tXlbomd) then
    call xlbomdIntegrator%setDefaultSCCParameters(minSCCiter, nSCCiter, sccTol)
  end if

  lpGeomOpt: do while (iGeoStep <= nGeoSteps)
        
    if (tSocket) then
      call socket%receive(coord0, latvec)
      cellVol = determinant33(latVec)
      recVec2p = latVec(:,:)
      call matinv(recVec2p)
      recVec2p = reshape(recVec2p, (/3, 3/), order=(/2, 1/))
      recVec = 2.0_dp * pi * recVec2p
      recCellVol = determinant33(recVec)

      if (tSCC) then
        call updateLatVecs_SCC(latVec, recVec, cellVol)
        mCutoff = max(mCutoff, getSCCCutoff())
      end if
      if (tDispersion) then
        call dispersion%updateLatVecs(latVec)
        mCutoff = max(mCutoff, dispersion%getRCutoff())
      end if
      call getCellTranslations(cellVec, rCellVec, latVec, recVec2p, &
          & mCutoff)
    end if

    if (restartFreq > 0 .and. (tGeoOpt .or. tMD)) then
      tWriteRestart = (iGeoStep == nGeoSteps .or. &
          & (mod(iGeoStep, restartFreq) == 0 ))
    else
      tWriteRestart = .false.
    end if

    if (tMD.and.tWriteRestart) then
      write(fdMD,*)"MD step:",iGeoStep
      call state(pMDIntegrator,fdMD)
    end if

    !! Write out geometry information
    write(*, '(/,A)') repeat('-', 80)
    if (tCoordOpt .and. tLatOpt) then
      write (*, "(/,'***  Geometry step: ',I0,', Lattice step: ',I0,/)") iGeoStep,iLatGeoStep
    else
      write (*, "(/,'***  Geometry step: ',I0,/)") iGeoStep
    end if
    
    
    if (tPeriodic) then
      invLatVec = transpose(latVec)
      call matinv(invLatVec)
      CellVol = abs(determinant33(latVec))
      
      ! derivative of pV term in Gibbs energy
      if (tStress.and.pressure/=0.0_dp) then
        call derivDeterminant33(derivCellVol,latVec)
        derivCellVol(:,:) = pressure * derivCellVol(:,:)
      end if
      
    end if

    !! Save old coordinates and fold coords to unit cell
    coord0Fold(:,:) = coord0
    if (tPeriodic) then
      call foldCoordToUnitCell(coord0Fold, latVec, recVec2p)
    end if

    !! Initialize neighborlists
    call updateNeighborListAndSpecies(coord, species, img2CentCell, iCellVec, &
        &neighborList, nAllAtom, coord0Fold, species0, mCutoff, rCellVec)
    nAllOrb = sum(orb%nOrbSpecies(species(1:nAllAtom)))

    !! Calculate neighborlist for SK and repulsive calculation
    call getNrOfNeighborsForAll(nNeighbor, neighborList, skRepCutoff)

    !! Reallocate Hamiltonian and overlap based on the new neighbor list
    call reallocateHS(ham, over, iPair, neighborList%iNeighbor, nNeighbor, &
        &orb, img2CentCell)

    !! Reallocate density matrixes if necessary
    if (size(ham, dim=1) > size(rhoPrim, dim=1)) then
      DEALLOCATE_(H0)
      ALLOCATE_(H0,(size(ham,dim=1)))
      DEALLOCATE_(rhoPrim)
      ALLOCATE_(rhoPrim,(size(ham,dim=1),nSpin))
      if (tImHam) then
        DEALLOCATE_(iRhoPrim)
        ALLOCATE_(iRhoPrim,(size(ham,dim=1),nSpin))
        DEALLOCATE_PARR(iHam)
        INITALLOCATE_PARR(iHam,(size(ham,dim=1),nSpin))
      end if
      if (tForces) then
        DEALLOCATE_(ERhoPrim)
        ALLOCATE_(ERhoPrim,(size(ham,dim=1)))
        DEALLOCATE_(ERhoPrim2)
        ALLOCATE_(ERhoPrim2, (size(ham, dim=1)))
      end if
    end if

    !! (Re)Initialize mixer
    if (tSCC) then
      call reset(pChrgMixer, nMixElements)
    end if

    !! Notify various modules about coordinate changes
    if (tSCC) then
      call updateCoords_SCC(coord, species, neighborList, img2CentCell)
    end if
    if (tDispersion) then
      call dispersion%updateCoords(neighborList, img2CentCell, coord, &
          & species0)
    end if
    if (t3rdFull) then
      call thirdOrd%updateCoords(neighborList, species)
    end if

    !! Build non-scc Hamiltonian and overlap
    call buildH0(H0, skHamCont, atomEigVal, coord, nNeighbor,&
        &  neighborList%iNeighbor, species, iPair, orb)
    call buildS(over, skOverCont, coord, nNeighbor, neighborList%iNeighbor,&
        & species, iPair, orb)

    !! Adapt electron temperature to MD, if necessary
    if (tSetFillingTemp) then
      call getTemperature(pTempProfile, tempElec)
    end if

    if (tXlbomd) then
      call xlbomdIntegrator%getSCCParameters(minSCCIter, nSCCiter, sccTol)
    end if

    if (tSCC .and. (.not. tAppendDetailedOut)) then
      write(*,"(' ',A5,A18,A18,A18)") "iSCC", " Total electronic ", &
          & "  Diff electronic ", "     SCC error    "
    end if

    tConverged = .false.

    energy%ETotal = 0.0_dp
    energy%atomTotal(:) = 0.0_dp

    !! Calculate repulsive energy
    call getERep(energy%atomRep, coord, nNeighbor, neighborList%iNeighbor, &
        &species, pRepCont, img2CentCell)
    energy%Erep = sum(energy%atomRep)
    if (tDispersion) then
      call dispersion%getEnergies(energy%atomDisp)
      energy%eDisp = sum(energy%atomDisp)
    else
      energy%atomDisp(:) = 0.0_dp
    end if

    potential%extAtom = 0.0_dp
    potential%extShell = 0.0_dp
    potential%extBlock = 0.0_dp


    if (tEField) then
      Efield(:) = EFieldStrength * EfieldVector(:)
      if (tTDEfield) then
        Efield(:) = Efield(:) &
            & * sin(EfieldOmega*deltaT*real(iGeoStep+EfieldPhase,dp))
      end if
      absEfield = sqrt(sum(Efield**2))
      if (tPeriodic) then
        do iAtom = 1, nAtom
          do iNeigh = 1, nNeighbor(iAtom)
            ii = neighborList%iNeighbor(iNeigh,iAtom)
            if (iCellVec(ii) /= 0) then ! overlap between atom in central
              !  cell and non-central cell
              if (abs(dot_product(cellVec(:,iCellVec(ii)),EfieldVector))&
                  & /= 0.0_dp) then ! component of electric field projects
                ! onto vector between cells
                write(tmpStr,"('Interaction between atoms ',I0,' and ', I0,&
                    &' crosses the saw-tooth discontinuity in the &
                    &electric field.')")iAtom,img2centcell(ii)
                call error(tmpStr)
              end if
            end if
          end do
        end do
        do iAtom = 1, nAtom
          potential%extAtom(iAtom,1)=dot_product(coord0Fold(:,iAtom),Efield)
        end do
      else
        do iAtom = 1, nAtom
          potential%extAtom(iAtom,1)=dot_product(coord(:,iAtom),Efield)
        end do
      end if
    else
      Efield = 0.0_dp
    end if

    call total_shift(potential%extShell, potential%extAtom, orb, species)
    call total_shift(potential%extBlock, potential%extShell, orb, species)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !! SCC-loop
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    iSCCIter = 1
    tStopSCC = .false.

    lpSCC: do while (iSCCiter <= nSCCIter)
      rhoPrim(:,:) = 0.0_dp
      if (tImHam) then
        iRhoPrim(:,:) = 0.0_dp
      end if

      ham(:,:) = 0.0_dp
      do ii = 1, size(H0)
        ham(ii,1) = h0(ii)
      end do

      !! Build various contribution to the Hamiltonian

      ! Reset this in case DFTB+U terms are added to
      ! potential%iorbitalBlock later in loop:
      potential%iorbitalBlock = 0.0_dp
      if (tDualSpinOrbit) then
        call shiftLS(potential%iorbitalBlock,xi,orb,species)
      end if

      if (.not. tSCC) then
        tConverged = .true.
        potential%intBlock = 0.0_dp
      else
        chargePerShell(:,:,:) = 0.0_dp ! hack for the moment to get charge
        ! and magnetization
        do iAtom = 1, nAtom
          iSp1 = species(iAtom)
          do iSh1 = 1, orb%nShell(iSp1)
            chargePerShell(iSh1,iAtom,1:nSpin) = &
                & chargePerShell(iSh1,iAtom,1:nSpin) + &
                & sum(qInput(orb%posShell(iSh1,iSp1): &
                & orb%posShell(iSh1+1,iSp1)-1,iAtom,1:nSpin),dim=1)
          end do
        end do

        potential%intAtom = 0.0_dp
        potential%intShell = 0.0_dp
        potential%intBlock = 0.0_dp
        call updateCharges_SCC(qInput, q0, orb, species, &
            &neighborList%iNeighbor, img2CentCell)
        call getShiftPerAtom(potential%intAtom)
        call getShiftPerL(potential%intShell)

        if (t3rdFull) then
          call thirdOrd%updateCharges(species0, neighborList, qInput, q0,&
              & img2CentCell, orb)
          call thirdOrd%getShiftPerAtom(shift3rd)
          potential%intAtom(:,1) = potential%intAtom(:,1) + shift3rd
          call thirdOrd%getShiftPerShell(orbresshift3rd)
          potential%intShell(:,:,1) = potential%intShell(:,:,1)&
              & + orbresshift3rd(:,:)
        end if

        call total_shift(potential%intShell, potential%intAtom, orb, species)
        !! Build spin contribution (if necessary)
        if (tSpin) then
          call addSpinShift(potential%intShell,chargePerShell,species,orb,spinW)
        end if

        call total_shift(potential%intBlock, potential%intShell, orb, species)

        if (tDFTBU) then !! Apply LDA+U correction (if necessary)
          potential%orbitalBlock = 0.0_dp
          if (tImHam) then
            call shift_DFTBU(potential%orbitalBlock,potential%iorbitalBlock, &
                & qBlockIn, qiBlockIn, species,orb, nDFTBUfunc, &
                & UJ, nUJ, niUJ, iUJ)
          else
            call shift_DFTBU(potential%orbitalBlock,qBlockIn,species,orb, &
                & nDFTBUfunc, UJ, nUJ, niUJ, iUJ)
          end if
          potential%intBlock = potential%intBlock + potential%orbitalBlock
        end if

      end if

      potential%intBlock = potential%intBlock + potential%extBlock

      call add_shift(ham,over,nNeighbor, neighborList%iNeighbor, &
          & species,orb,iPair,nAtom,img2CentCell,potential%intBlock)

      if (tImHam) then
        iHam = 0.0_dp
        call add_shift(iHam,over,nNeighbor, neighborList%iNeighbor, &
            & species,orb,iPair,nAtom,img2CentCell,potential%iorbitalBlock)
        iHam(:,:) = 2.0_dp*iHam(:,:)
      end if

      ! hack due to not using Pauli-type structure for diagonalisation
      ! etc.
      if (nSpin > 1) then
        ham(:,:) = 2.0_dp * ham
      end if

      if (nSpin /= 4) then

        if (nSpin == 2) then
          call qm2ud(ham)
        end if

        if (tWriteRealHS .or. tWriteHS) then
          ! Write out matrices if necessary and quit.
          if (tImHam) then
            call writeHS(tWriteHS, tWriteRealHS, ham, over, &
                & neighborList%iNeighbor, nNeighbor, iAtomStart, iPair, &
                & img2CentCell, kPoint, iCellVec, cellVec, iHam=iHam)
          else
            call writeHS(tWriteHS, tWriteRealHS, ham, over, &
                & neighborList%iNeighbor, nNeighbor, iAtomStart, iPair, &
                & img2CentCell, kPoint, iCellVec, cellVec)
          end if
          write (*, "(A)") "Hamilton/Overlap written, exiting program."
          stop
        end if

        spinDiag: do iSpin = 1, nSpin
          if (tStoreEigvecs) then
            iSpin2 = 1
          else
            iSpin2 = iSpin
          end if

          ! Solve eigenproblem for real or complex Hamiltionian
          if (tRealHS) then
            call diagonalize(HSqrReal(:,:,iSpin2), SSqrReal, eigen(:,1,iSpin), ham(:,iSpin), over,&
                & neighborList%iNeighbor, nNeighbor, iAtomStart, iPair, img2CentCell, solver, 'V')
            if (tStoreEigvecs) then
              call reset(storeEigvecsReal(iSpin), [sqrHamSize, sqrHamSize])
              call push(storeEigvecsReal(iSpin), HSqrReal(:,:,iSpin2))
            end if
          else
            do nk = 1, nKPoint
              if (tStoreEigvecs) then
                iK2 = 1
              else
                iK2 = nK
              end if
              call diagonalize(HSqrCplx(:,:,iK2,iSpin2), SSqrCplx, eigen(:,nk,iSpin), ham(:,iSpin),&
                  & over,kPoint(:,nk), neighborList%iNeighbor, nNeighbor, iCellVec, cellVec,&
                  & iAtomStart, iPair, img2CentCell, solver, 'V')
              if (tStoreEigvecs) then
                call reset(storeEigvecsCplx(iSpin), [sqrHamSize, sqrHamSize])
                call push(storeEigvecsCplx(iSpin), HSqrCplx(:,:,iK2,iSpin2))
              end if
            end do
          end if
        end do spinDiag

        call getFillingsAndBandEnergies(eigen, nEl, nSpin, tempElec, kWeight, tSpinSharedEf,&
            & tFillKSep, tFixEf, iDistribFn, Ef, filling, Eband, TS, E0)

        ! Create density matrices
        spinDensMat: do iSpin = 1, nSpin
          if (tStoreEigvecs) then
            iSpin2 = 1
          else
            iSpin2 = iSpin
          end if

          if (tRealHS) then
            if (tStoreEigvecs) then
              call get(storeEigvecsReal(iSpin), HSqrReal(:,:,iSpin2))
            end if

            if (tDensON2) then
              call makeDensityMatrix(SSqrReal, HSqrReal(:,:,iSpin2), filling(:,1,iSpin),&
                  & neighborlist%iNeighbor, nNeighbor, orb, iAtomStart, img2CentCell)
            else
              call makeDensityMatrix(SSqrReal, HSqrReal(:,:,iSpin2), filling(:,1,iSpin))
            end if
            
            if (tLinResp .and. tLinRespZVect) then
              rhoSqrReal(:,:,iSpin) = SSqrReal
            end if
            
            call packHS(rhoPrim(:,iSpin), SSqrReal, neighborlist%iNeighbor, nNeighbor, orb%mOrb,&
                  & iAtomStart, iPair, img2CentCell)

            ! Store density matrix for later use if using linear response
            if (tLinResp .and. tForces) then
              rhoSqrReal(:,:,iSpin) = SSqrReal
            end if

          else

            do nk = 1, nKPoint
              if (tStoreEigvecs) then
                iK2 = 1
                call get(storeEigvecsCplx(iSpin), HSqrCplx(:,:,iK2, iSpin2))
              else
                iK2 = nK
              end if

              if (tDensON2) then
                call makeDensityMatrix(SSqrCplx, HSqrCplx(:,:,iK2,iSpin2), filling(:,nK,iSpin),&
                    & neighborlist%iNeighbor, nNeighbor, orb, iAtomStart, img2CentCell)
              else
                call makeDensityMatrix(SSqrCplx, HSqrCplx(:,:,iK2,iSpin2), filling(:,nK,iSpin))
              end if
              call packHS(rhoPrim(:,iSpin), SSqrCplx, kPoint(:,nK), kWeight(nk),&
                  & neighborList%iNeighbor, nNeighbor, orb%mOrb, iCellVec, cellVec, iAtomStart,&
                  & iPair, img2CentCell)
            end do
          end if
        end do spinDensMat

        do iSpin = 1, nSpin
          if (tWriteBandDat) then
            if (iSpin == 1) then
              open(unit=fdBand, file=bandOut, action="write")
            end if
            do nk=1,nKPoint
              write(fdBand,*)'KPT ',nk,' SPIN ', iSpin, &
                  &' KWEIGHT ', kweight(nK)
              do iEgy=1,sqrHamSize
                write(fdBand, formatEnergy) Hartree__eV*eigen(iEgy,nk,iSpin),&
                    & filling(iEgy,nk,iSpin)
              end do
              write(fdBand,*)
            end do
            if (iSpin == nSpin) then
              close(fdBand)
            end if
          end if
        end do

        call ud2qm(rhoPrim)

      else !  (nSpin == 4) then

        if (tRealHS) then
          if (tImHam) then
            call diagonalize(HSqrCplx(:,:,1,1), SSqrCplx, eigen(:,1,1), ham, over,&
                & neighborList%iNeighbor, nNeighbor, iAtomStart, iPair, img2CentCell, solver, 'V', &
                & iHam=iHam)
          else if (tSpinOrbit) then
            call diagonalize(HSqrCplx(:,:,1,1), SSqrCplx, eigen(:,1,1), ham, over,&
                & neighborList%iNeighbor, nNeighbor, iAtomStart, iPair, img2CentCell, solver, 'V', &
                & xi,orb,species)
          else
            call diagonalize(HSqrCplx(:,:,1,1), SSqrCplx, eigen(:,1,1), ham, over,&
                & neighborList%iNeighbor, nNeighbor, iAtomStart, iPair, img2CentCell, solver, 'V')
          end if
          if (tStoreEigvecs) then
            call reset(storeEigvecsCplx(1), [sqrHamSize, sqrHamSize])
            call push(storeEigvecsCplx(1), HSqrCplx(:,:,1,1))
          end if
        else
          do nk = 1, nKPoint
            if (tStoreEigvecs) then
              iK2 = 1
            else
              iK2 = nK
            end if
            if (tImHam) then
              call diagonalize(HSqrCplx(:,:,iK2,1), SSqrCplx, eigen(:,nk,1), ham, over,&
                  & kPoint(:,nk), neighborList%iNeighbor, nNeighbor, iCellVec, cellVec, iAtomStart,&
                  & iPair, img2CentCell, solver, 'V', iHam=iHam)
            else if (tSpinOrbit) then
              call diagonalize(HSqrCplx(:,:,iK2,1), SSqrCplx, eigen(:,nK,1), ham, over,&
                  & kPoint(:,nK), neighborList%iNeighbor, nNeighbor, iCellVec, cellVec, iAtomStart,&
                  & iPair, img2CentCell, solver, 'V', xi,orb,species)
            else
              call diagonalize(HSqrCplx(:,:,iK2,1), SSqrCplx, eigen(:,nK,1), ham, over,&
                  & kPoint(:,nK), neighborList%iNeighbor, nNeighbor, iCellVec, cellVec, iAtomStart,&
                  & iPair, img2CentCell, solver, 'V')
            end if
            if (tStoreEigvecs) then
              call reset(storeEigvecsCplx(1), [sqrHamSize, sqrHamSize])
              call push(storeEigvecsCplx(1), HSqrCplx(:,:,iK2,1))
            end if
          end do
        end if

        call getFillingsAndBandEnergies(eigen, nEl, nSpin, tempElec, kWeight, tSpinSharedEf,&
            & tFillKSep, tFixEf, iDistribFn, Ef, filling, Eband, TS, E0)

        ! Pauli structure of eigenvectors
        filling(:,1:nKPoint,1) = 2.0_dp * filling(:,1:nKPoint,1)

        SSqrCplx = 0.0_dp
        if (tSpinOrbit) then
          energy%atomLS = 0.0_dp
          if (.not.tDualSpinOrbit) then
            ALLOCATE_(rVecTemp,(nAtom))
          end if
        end if
        if (tImHam .and. tMulliken) then
          orbitalL = 0.0_dp
        end if
        nkLoop4: do nk = 1, nKPoint
          if (tStoreEigvecs) then
            iK2 = 1
            call get(storeEigvecsCplx(1), HSqrCplx(:,:,iK2,1))
          else
            iK2 = nK
          end if
          if (tDensON2) then
            call error("Currently missing.")
          else
            call makeDensityMatrix(SSqrCplx, HSqrCplx(:,:,iK2,1), &
                &filling(:,nK,1))
          end if

          if (tSpinOrbit .and. .not. tDualSpinOrbit) then
            rVecTemp = 0.0_dp
            call getEnergySpinOrbit(rVecTemp,SSqrCplx,iAtomStart, &
                & xi, orb, species)
            energy%atomLS = energy%atomLS + kWeight(nk)*rVecTemp
            if (tMulliken) then
              orbitalLPart = 0.0_dp
              call getL(orbitalLPart,SSqrCplx,iAtomStart, orb, species)
              orbitalL = orbitalL + kWeight(nk) * orbitalLPart
            end if
          end if

          if (tRealHS) then
            call packHS(rhoPrim(:,:), SSqrCplx, neighborlist%iNeighbor, &
                &nNeighbor, orb%mOrb, iAtomStart, iPair, img2CentCell)
            if (tImHam) then
              call iPackHS(iRhoPrim(:,:), SSqrCplx, neighborlist%iNeighbor, &
                  &nNeighbor, orb%mOrb, iAtomStart, iPair, img2CentCell)
            end if
          else
            call packHS(rhoPrim, SSqrCplx, kPoint(:,nK), &
                &kWeight(nk), neighborList%iNeighbor, nNeighbor, orb%mOrb, &
                &iCellVec, cellVec, iAtomStart, iPair, img2CentCell)
            if (tImHam) then
              call iPackHS(iRhoPrim, SSqrCplx, kPoint(:,nK), &
                  & kWeight(nk), neighborlist%iNeighbor, &
                  & nNeighbor, orb%mOrb, iCellVec, cellVec, iAtomStart, &
                  & iPair, img2CentCell)
            end if
          end if

        end do nkLoop4
        if (tSpinOrbit .and. .not. tDualSpinOrbit) then
          DEALLOCATE_(rVecTemp)
          energy%ELS = sum(energy%atomLS(:))
        end if
        filling(:,1:nKPoint,1) = 0.5_dp * filling(:,1:nKPoint,1)

        if (tWriteBandDat) then
          open(unit=fdBand, file=bandOut, action="write")
          do nk=1,nKPoint
            write(fdBand, *) 'KPT ',nk,' SPIN ', 1, ' KWEIGHT ', kweight(nK)
            do iEgy=1,sqrHamSize
              write(fdBand, formatEnergy) Hartree__eV*eigen(iEgy,nk,1),&
                  & filling(iEgy,nk,1)
            end do
            write(fdBand,*)
          end do
          close(fdBand)
        end if

      end if ! end of nSpin == 4 case


      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      !! Mulliken analysis
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      if (tMulliken) then
        qOutput(:,:,:) = 0.0_dp
        do iSpin = 1, nSpin
          call mulliken(qOutput(:,:,iSpin), over, rhoPrim(:,iSpin), &
              &orb, neighborList%iNeighbor, nNeighbor, img2CentCell, iPair)
        end do
      end if

      if (tImHam) then
        qiBlockOut(:,:,:,:) = 0.0_dp
        energy%atomLS = 0.0_dp
        do iSpin = 1, nSpin
          call skewMulliken(qiBlockOut(:,:,:,iSpin), over, iRhoPrim(:,iSpin), &
              &orb, neighborList%iNeighbor, nNeighbor, img2CentCell, iPair)
        end do
        call getL(orbitalL,qiBlockOut,orb,species)
        if (tDualSpinOrbit) then
          call getEnergySpinOrbit(energy%atomLS,qiBlockOut,xi,orb,species)
          energy%ELS = sum(energy%atomLS(:))
        end if
        qBlockOut(:,:,:,:) = 0.0_dp
      end if

      if (tDFTBU) then
        qBlockOut(:,:,:,:) = 0.0_dp
        do iSpin = 1, nSpin
          call mulliken(qBlockOut(:,:,:,iSpin), over, rhoPrim(:,iSpin), &
              &orb, neighborList%iNeighbor, nNeighbor, img2CentCell, iPair)
        end do
      end if

      ! Note: if XLBOMD is active, potential created with ingoing charges
      ! is needed later, therefore it should not be zeroed out.
      if (tSCC .and. .not. tXlbomd) then

        potential%intAtom = 0.0_dp
        potential%intShell = 0.0_dp
        potential%intBlock = 0.0_dp

        chargePerShell(:,:,:) = 0.0_dp ! hack for the moment to get charge
        ! and magnetization
        do iAtom = 1, nAtom
          iSp1 = species(iAtom)
          do iSh1 = 1, orb%nShell(iSp1)
            chargePerShell(iSh1,iAtom,1:nSpin) = &
                & chargePerShell(iSh1,iAtom,1:nSpin) + &
                & sum(qOutput(orb%posShell(iSh1,iSp1): &
                & orb%posShell(iSh1+1,iSp1)-1,iAtom,1:nSpin),dim=1)
          end do
        end do

        ! recalculate the SCC shifts for the output charge.

        !! SCC contribution is calculated with the output charges.
        call updateCharges_SCC(qOutput, q0, orb, species, &
            &neighborList%iNeighbor, img2CentCell)

        call getShiftPerAtom(potential%intAtom)
        call getShiftPerL(potential%intShell)
        if (t3rdFull) then
          call thirdOrd%updateCharges(species0, neighborList, qOutput, q0,&
              & img2CentCell, orb)
          call thirdOrd%getShiftPerAtom(shift3rd)
          potential%intAtom(:,1) = potential%intAtom(:,1) + shift3rd
          call thirdOrd%getShiftPerShell(orbresshift3rd)
          potential%intShell(:,:,1) = potential%intShell(:,:,1)&
              & + orbresshift3rd(:,:)
        end if

        call total_shift(potential%intShell, potential%intAtom, orb, species)

        !! Build spin contribution (if necessary)
        if (tSpin) then
          call addSpinShift(potential%intShell,chargePerShell,species,orb,spinW)
        end if

        call total_shift(potential%intBlock, potential%intShell, orb, species)
      end if


      !! Calculate energies

      ! non-SCC part
      energy%EnonSCC = 0.0_dp
      energy%atomNonSCC(:) = 0.0_dp

      call mulliken(energy%atomNonSCC(:), rhoPrim(:,1), H0,orb,&
          &neighborList%iNeighbor, nNeighbor, img2CentCell, iPair)
      energy%EnonSCC =  sum(energy%atomNonSCC)

      if (tEfield) then ! energy in external field
        energy%atomExt = -sum( q0(:, :, 1) - qOutput(:, :, 1),dim=1) &
            & * potential%extAtom(:,1)
        energy%Eext =  sum(energy%atomExt)
      else
        energy%Eext = 0.0_dp
        energy%atomExt = 0.0_dp
      end if

      if (tSCC) then
        if (tXlbomd) then
          call getEnergyPerAtom_SCC_Xlbomd(species, orb, qOutput, q0, &
              & energy%atomSCC)
        else
          call getEnergyPerAtom_SCC(energy%atomSCC)
        end if
        energy%eSCC = sum(energy%atomSCC)
        if (t3rdFull) then
          if (tXlbomd) then
            call thirdOrd%getEnergyPerAtomXlbomd(qOutput, q0, species, orb,&
                & energy%atom3rd)
          else
            call thirdOrd%getEnergyPerAtom(energy%atom3rd)
          end if
          energy%e3rd = sum(energy%atom3rd)
        end if

        if (nSpin > 1) then
          energy%atomSpin(:) = 0.5_dp * sum(sum(potential%intShell(:,:,2:nSpin)&
              & * chargePerShell(:,:,2:nSpin), dim=1),dim=2)
          energy%Espin = sum(energy%atomSpin)
        else
          energy%atomSpin(:) = 0.0_dp
          energy%eSpin = 0.0_dp
        end if
      end if


      potential%iorbitalBlock = 0.0_dp
      if (tDualSpinOrbit) then
        call shiftLS(potential%iorbitalBlock,xi,orb,species)
      end if

      if (tDFTBU) then
        energy%atomDftbu(:) = 0.0_dp
        if (.not. tImHam) then
          call E_DFTBU(energy%atomDftbu,qBlockOut,species,orb, &
              & nDFTBUfunc, UJ, nUJ, niUJ, iUJ)
        else
          call E_DFTBU(energy%atomDftbu,qBlockOut,species,orb, &
              & nDFTBUfunc, UJ, nUJ, niUJ, iUJ, qiBlockOut)
        end if
        energy%Edftbu = sum(energy%atomDftbu(:))
        potential%orbitalBlock = 0.0_dp

        if (tImHam) then
          call shift_DFTBU(potential%orbitalBlock,potential%iorbitalBlock, &
              & qBlockOut,qiBlockOut, species,orb, nDFTBUfunc, UJ, nUJ, niUJ, &
              & iUJ)
        else
          call shift_DFTBU(potential%orbitalBlock,qBlockOut,species,orb, &
              & nDFTBUfunc, UJ, nUJ, niUJ, iUJ)
        end if
      else
        energy%Edftbu = 0.0_dp
      end if


      energy%Eelec = energy%EnonSCC + energy%ESCC + energy%Espin &
          & + energy%ELS + energy%Edftbu + energy%Eext + energy%e3rd

      energy%atomElec(:) = energy%atomNonSCC(:) &
          & + energy%atomSCC(:) + energy%atomSpin(:) + energy%atomDftbu(:) &
          & + energy%atomLS(:) + energy%atomExt(:) + energy%atom3rd(:)

      energy%atomTotal(:) = energy%atomElec(:) + energy%atomRep(:) + &
          & energy%atomDisp(:)

      energy%Etotal = energy%Eelec + energy%Erep + energy%eDisp

      !! Stop SCC if appropriate stop file is present (We need this query here
      !! since the following block contains a check iSCCIter /= nSCCIter)
      inquire(file=fStopSCC, exist=tStopSCC)
      if (tStopSCC) then
        write (*,*) "Stop file '" // fStopSCC // "' found."
        nSCCIter = iSCCIter
        write (*,*) "Setting max number of scc cycles to current cycle."
      end if


      !! Mix charges
      if (tSCC) then
        qOutRed = 0.0_dp
        if (nSpin == 2) then
          call qm2ud(qOutput)
          if (tDFTBU) then
            call qm2ud(qBlockOut)
          end if
        end if
        call OrbitalEquiv_reduce(qOutput, iEqOrbitals, orb, &
            & qOutRed(1:nIneqOrb))
        if (tDFTBU) then
          call AppendBlock_reduce( qBlockOut, iEqBlockDFTBU, orb, &
              & qOutRed )
          if (tImHam) then
            call AppendBlock_reduce( qiBlockOut, iEqBlockDFTBULS, orb, &
                & qOutRed, skew=.true. )
          end if
        end if
        if (nSpin == 2) then
          call ud2qm(qOutput)
          if (tDFTBU) then
            call ud2qm(qBlockOut)
          end if
        end if

        qDiffRed(:) = qOutRed(:) - qInpRed(:)
        sccErrorQ = maxval(abs(qDiffRed))

        tConverged = (sccErrorQ < sccTol) .and. &
            & (iSCCiter >= minSCCIter .or. tReadChrg .or. iGeoStep > 0)
        if ((.not. tConverged) .and. iSCCiter /= nSCCiter) then
          !! Avoid mixing of spin unpolarised density for spin polarised
          !! cases, this is only a problem in iteration 1, as there is
          !! only the (spin unpolarised!) atomic input density at that
          !! point. (Unless charges had been initialized externally)
          if ((iSCCIter + iGeoStep) == 1 .and. (nSpin > 1.or.tDFTBU) &
              & .and. .not.tReadChrg) then
            qInput(:,:,:) = qOutput(:,:,:)
            qInpRed(:) = qOutRed(:)
            if (tDFTBU) then
              qBlockIn(:,:,:,:) = qBlockOut(:,:,:,:)
              if (tSpinOrbit) then
                qiBlockIn(:,:,:,:) = qiBlockOut(:,:,:,:)
              end if
            end if
          else
            call mix(pChrgMixer, qInpRed, qDiffRed)
            call OrbitalEquiv_expand(qInpRed(1:nIneqOrb), iEqOrbitals, &
                & orb, qInput)
            if (tDFTBU) then
              qBlockIn = 0.0_dp
              call Block_expand( qInpRed ,iEqBlockDFTBU, orb, &
                  & qBlockIn, species0, nUJ, niUJ, iUJ, orbEquiv=iEqOrbitals )
              if (tSpinOrbit) then
                call Block_expand( qInpRed ,iEqBlockDFTBULS, orb, &
                    & qiBlockIn, species0, nUJ, niUJ, iUJ, skew=.true. )
              end if
            end if
            if (nSpin == 2) then
              call ud2qm(qInput)
              if (tDFTBU) then
                call ud2qm(qBlockIn)
              end if
            end if
          end if
        end if

      end if


      !! Clear detailed.out if necessary
      if (tWriteDetailedOut .and. .not. tAppendDetailedOut) then
        close(fdUser)
        open(fdUser, file=userOut, position="rewind", status="replace")
        select case(iDistribFn)
        case(0)
          write(fdUser,*)'Fermi distribution function'
        case(1)
          write(fdUser,*)'Gaussian distribution function'
        case default
          write(fdUser,*)'Methfessel-Paxton distribution function order'&
              &,iDistribFn
        end select
        write(fdUser,*) ""
      end if

      if (tSCC) then
        if (iSCCiter > 1) then
          rTmp = (energy%Eelec-Eold)
        else
          rTmp = 0.0_dp
        end if
        Eold = energy%Eelec

        if (tWriteDetailedOut) then
          if (nGeoSteps > 0) then
            if (tMD) then
              write(fdUser,"('MD step: ',I0)") iGeoStep
            elseif (tDerivs) then
              write(fdUser,"('Difference derivative step: ',I0)") iGeoStep
            else
              if (tCoordOpt .and. tLatOpt) then
                write (fdUser, "('Geometry optimization step: ',I0, &
                    & ', Lattice step: ',I0)") &
                    & iGeoStep,iLatGeoStep
              else
                write(fdUser,"('Geometry optimization step: ',I0)") iGeoStep
              end if
            end if
          else
            write(fdUser,*) "Calculation with static geometry"
          end if
          write(fdUser,*)''
          write (fdUser,'(/A)') repeat("*", 80)
          write(fdUser,"(' ',A5,A18,A18,A18)") "iSCC", " Total electronic ", &
              & "  Diff electronic ", "     SCC error    "
          write(fdUser,"(I5,E18.8,E18.8,E18.8,E18.8)") iSCCIter, &
              & energy%Eelec, rTmp, sccErrorQ
          write (fdUser,'(A)') repeat("*", 80)
          write(fdUser,*)""
        end if

        if (tDFTBU) then
          write(*,"(I5,E18.8,E18.8,E18.8)") iSCCIter, energy%Eelec, rTmp, &
              &sccErrorQ
        else
          write(*,"(I5,E18.8,E18.8,E18.8)") iSCCIter, energy%Eelec, rTmp, &
              & sccErrorQ
        end if
      end if


      !! Not writing any restarting info if not converged and minimal number of
      !! SCC iterations not done.
      if (restartFreq > 0 .and. .not.(tMD .or. tGeoOpt .or. tDerivs) .and. nSCCIter > 1) then
        if (tConverged .or. ((iSCCIter >= minSCCIter &
            & .or. tReadChrg .or. iGeoStep > 0) &
            &.and. (iSCCIter == nSCCIter .or. mod(iSCCIter, restartFreq) ==&
            & 0))) then
          if (tMulliken.and.tSCC) then
            if (tDFTBU) then
              if (tSpinOrbit) then
                call writeQToFile(qInput, fChargeIn, orb, qBlockIn, qiBlockIn)
              else
                call writeQToFile(qInput, fChargeIn, orb, qBlockIn)
              end if
            else
              call writeQToFile(qInput, fChargeIn, orb)
            end if
            print "('>> Charges saved for restart in ',A)", fChargeIn
          end if
        end if
      end if

      if (tWriteDetailedOut) then
        if (nMovedAtom > 0 .and. .not. tDerivs) then
          write (fdUser,*) "Coordinates of moved atoms (au):"
          do ii = 1, nMovedAtom
            write(fdUser,formatGeoOut) indMovedAtom(ii), &
                &pCoord0Out(:, indMovedAtom(ii))
          end do
          write (fdUser,*) ""
        end if

        !! Write out atomic charges
        if (tPrintMulliken) then
          write (fdUser, "(/,A)") " Net atomic charges (e)"
          write (fdUser, "(1X,A5,1X,A16)")" Atom", " Net charge"
          do ii = 1, nAtom
            write (fdUser, "(1X,I5,1X,F16.8)") ii, &
                &sum(q0(:, ii, 1) - qOutput(:, ii, 1))
          end do
          write(fdUser,*) ""
        end if
        lpSpinPrint: do iSpin = 1, nSpinHams
          if (nSpin == 2) then
            write(fdUser,*) 'COMPONENT = ',trim(spinName(iSpin))
          else
            write(fdUser,*) 'COMPONENT = ',trim(quaternionName(iSpin))
          end if
          write(fdUser,*) ' '
          write(fdUser,*)'Eigenvalues /H'
          do iEgy = 1, sqrHamSize
            write(fdUser, formatEigen) (eigen(iEgy,ii,iSpin),ii=1,nKPoint)
          end do
          write(fdUser,*) ""
          write(fdUser,*)'Eigenvalues /eV'
          do iEgy = 1, sqrHamSize
            write(fdUser, formatEigen) &
                &(Hartree__eV*eigen(iEgy,ii,iSpin),ii=1,nKPoint)
          end do
          write (fdUser,*) ''
          write(fdUser,*)'Fillings'
          do iEgy = 1, sqrHamSize
            write(fdUser, formatEnergy) (filling(iEgy,ii,iSpin),ii=1,nKPoint)
          end do
          write (fdUser,*) ""
        end do lpSpinPrint
        if (nSpin == 4) then
          if (tPrintMulliken) then
            do jj = 1, 4
              write (fdUser,"(' Nr. of electrons (',A,'):',F16.8)") &
                  & quaternionName(jj),sum(qOutput(:, :,jj))
              write (fdUser,*) ""
              write (fdUser, "(' Atom populations (',A,')')") quaternionName(jj)
              write (fdUser, "(1X,A5,1X,A16)")" Atom", " Population"

              do ii = 1, nAtom
                write (fdUser, "(1X,I5,1X,F16.8)") ii, sum(qOutput(:, ii, jj))
              end do
              write (fdUser,*) ''
              write (fdUser, "(' l-shell populations (',A,')')") quaternionName(jj)
              write (fdUser, "(1X,A5,1X,A3,1X,A3,1X,A16)")" Atom", "Sh.", &
                  &"  l", " Population"
              do ii = 1, nAtom
                iSp1 = species(ii)
                do iSh1 = 1, orb%nShell(iSp1)
                  write (fdUser, "(1X,I5,1X,I3,1X,I3,1X,F16.8)") ii, iSh1, &
                      &orb%angShell(iSh1,iSp1), &
                      &sum(qOutput(orb%posShell(iSh1,iSp1)&
                      &:orb%posShell(iSh1+1,iSp1)-1, ii, jj))
                end do
              end do
              write (fdUser,*) ''
              write (fdUser, "(' Orbital populations (',A,')')") &
                  & quaternionName(jj)
              write (fdUser, "(1X,A5,1X,A3,1X,A3,1X,A3,1X,A16)") " Atom", &
                  & "Sh.","  l","  m", " Population"
              do ii = 1, nAtom
                iSp1 = species(ii)
                do iSh1 = 1, orb%nShell(iSp1)
                  ang = orb%angShell(iSh1, iSp1)
                  do kk = 0, 2 * ang
                    write (fdUser, "(' ',I5,1X,I3,1X,I3,1X,I3,1X,F16.8)") &
                        &ii, iSh1, ang, kk - ang, &
                        &qOutput(orb%posShell(iSh1,iSp1)+kk, ii, jj)
                  end do
                end do
              end do
              write (fdUser,*) ''
            end do
          end if
          if (tDFTBU) then
            do jj = 1, 4
              write (fdUser, "(' Block populations (',A,')')") &
                  & quaternionName(jj)
              do ii = 1, nAtom
                iSp1 = species(ii)
                write(fdUser,*)'Atom',ii
                do kk = 1, orb%nOrbSpecies(iSp1)
                  write(fdUser,"(16F8.4)") &
                      & qBlockOut(1:orb%nOrbSpecies(iSp1),kk,ii,jj)
                end do
                write (fdUser,*) ''
              end do
            end do
          end if
          if (tImHam .and. tPrintMulliken) then
            write (fdUser,*) ''
            write (fdUser,*) ' Electron angular momentum (mu_B/hbar)'
            write (fdUser, "(2X,A5,T10,A3,T14,A1,T20,A1,T35,A9)")"Atom", "Sh.",&
                &"l", "S", "Momentum"
            do ii = 1, nAtom
              iSp1 = species(ii)
              do iSh1 = 1, orb%nShell(iSp1)
                write (fdUser, "(1X,I5,1X,I3,1X,I3,1X,F14.8,' :',3F14.8)") &
                    & ii, iSh1, orb%angShell(iSh1,iSp1), &
                    & 0.5_dp*sqrt(sum(sum(qOutput(orb%posShell(iSh1,iSp1)&
                    & :orb%posShell(iSh1+1,iSp1)-1, ii, 2:4),dim=1)**2)), &
                    & -gfac*0.25_dp*sum(qOutput(orb%posShell(iSh1,iSp1)&
                    & :orb%posShell(iSh1+1,iSp1)-1, ii, 2:4),dim=1)
              end do
            end do
            write (fdUser,*) ''
            write (fdUser,*) ' Orbital angular momentum (mu_B/hbar)'
            write (fdUser, "(2X,A5,T10,A3,T14,A1,T20,A1,T35,A9)")"Atom", "Sh.",&
                &"l", "L", "Momentum"
            do ii = 1, nAtom
              iSp1 = species(ii)
              do iSh1 = 1, orb%nShell(iSp1)
                write (fdUser, "(1X,I5,1X,I3,1X,I3,1X,F14.8,' :',3F14.8)") &
                    & ii, iSh1, orb%angShell(iSh1,iSp1), &
                    & sqrt(sum(orbitalL(1:3,iSh1,ii)**2)),-orbitalL(1:3,iSh1,ii)
              end do
            end do
            write (fdUser,*) ''
            write (fdUser,*) ' Total angular momentum (mu_B/hbar)'
            write (fdUser, "(2X,A5,T10,A3,T14,A1,T20,A1,T35,A9)")"Atom", "Sh.",&
                &"l", "J", "Momentum"
            angularMomentum = 0.0_dp
            do ii = 1, nAtom
              iSp1 = species(ii)
              do iSh1 = 1, orb%nShell(iSp1)
                write (fdUser, "(1X,I5,1X,I3,1X,I3,1X,F14.8,' :',3F14.8)") &
                    & ii, iSh1, orb%angShell(iSh1,iSp1), sqrt(sum((&
                    & orbitalL(1:3,iSh1,ii) &
                    & +sum(0.5_dp*qOutput(orb%posShell(iSh1,iSp1)&
                    & :orb%posShell(iSh1+1,iSp1)-1, ii, 2:4),dim=1))**2)), &
                    & -orbitalL(1:3,iSh1,ii) &
                    & -gfac*0.25_dp*sum(qOutput(orb%posShell(iSh1,iSp1)&
                    & :orb%posShell(iSh1+1,iSp1)-1, ii, 2:4),dim=1)
                angularMomentum(1:3) = angularMomentum(1:3) &
                    &  -orbitalL(1:3,iSh1,ii) &
                    & -gfac*0.25_dp*sum(qOutput(orb%posShell(iSh1,iSp1)&
                    & :orb%posShell(iSh1+1,iSp1)-1, ii, 2:4),dim=1)
              end do
            end do
            write (fdUser,*) ''
          end if
        else
          if (nSpin == 2) then
            call qm2ud(qOutput)
            if (tDFTBU) then
              call qm2ud(qBlockOut)
            end if
          end if
          lpSpinPrint2: do iSpin = 1, nSpin
            if (tPrintMulliken) then
              write (fdUser,"(' Nr. of electrons (',A,'):',F16.8)") &
                  & trim(spinName(iSpin)),sum(qOutput(:, :,iSpin))
              write (fdUser, "(' Atom populations (',A,')')")&
                  & trim(spinName(iSpin))
              write (fdUser, "(1X,A5,1X,A16)")" Atom", " Population"

              do ii = 1, nAtom
                write (fdUser, "(1X,I5,1X,F16.8)")ii, sum(qOutput(:, ii, iSpin))
              end do
              write (fdUser,*) ''
              write (fdUser, "(' l-shell populations (',A,')')") &
                  & trim(spinName(iSpin))
              write (fdUser, "(1X,A5,1X,A3,1X,A3,1X,A16)")" Atom", "Sh.", &
                  &"  l", " Population"
              do ii = 1, nAtom
                iSp1 = species(ii)
                do iSh1 = 1, orb%nShell(iSp1)
                  write (fdUser, "(1X,I5,1X,I3,1X,I3,1X,F16.8)") ii, iSh1, &
                      &orb%angShell(iSh1,iSp1), &
                      &sum(qOutput(orb%posShell(iSh1,iSp1)&
                      &:orb%posShell(iSh1+1,iSp1)-1, ii, iSpin))
                end do
              end do
              write (fdUser,*) ''
              write (fdUser, "(' Orbital populations (',A,')')") &
                  & trim(spinName(iSpin))
              write (fdUser, "(1X,A5,1X,A3,1X,A3,1X,A3,1X,A16)")" Atom", "Sh.",&
                  &"  l","  m", " Population"
              do ii = 1, nAtom
                iSp1 = species(ii)
                do iSh1 = 1, orb%nShell(iSp1)
                  ang = orb%angShell(iSh1, iSp1)
                  do kk = 0, 2 * ang
                    write (fdUser, "(' ',I5,1X,I3,1X,I3,1X,I3,1X,F16.8)") &
                        &ii, iSh1, ang, kk - ang, &
                        &qOutput(orb%posShell(iSh1,iSp1)+kk, ii, iSpin)
                  end do
                end do
              end do
              write (fdUser,*) ''
            end if
            if (tDFTBU) then
              write (fdUser, "(' Block populations (',A,')')") &
                  & trim(spinName(iSpin))
              do ii = 1, nAtom
                iSp1 = species(ii)
                write(fdUser,*)'Atom',ii
                do kk = 1, orb%nOrbSpecies(iSp1)
                  write(fdUser,"(16F8.4)") &
                      & qBlockOut(1:orb%nOrbSpecies(iSp1),kk,ii,iSpin)
                end do
              end do
              write (fdUser,*) ''
            end if
          end do lpSpinPrint2
          if (nSpin == 2) then
            call ud2qm(qOutput)
            if (tDFTBU) then
              call ud2qm(qBlockOut)
            end if
          end if
        end if

        if (nSpin == 2) then
          call qm2ud(qOutput)
          call qm2ud(qInput)
        end if
        lpSpinPrint3: do iSpin = 1, nSpinHams
          if (nSpin == 2) then
            write(fdUser,*)'Spin ',trim(spinName(iSpin))
          end if
          write(fdUser,format2U) 'Fermi level', Ef(iSpin),"H", &
              & Hartree__eV*Ef(iSpin),'eV'
          write(fdUser,format2U) 'Band energy', Eband(iSpin),"H", &
              &Hartree__eV*Eband(iSpin),'eV'
          write(fdUser,format2U)'TS', TS(iSpin),"H", Hartree__eV*TS(iSpin),'eV'
          write(fdUser,format2U) 'Band free energy (E-TS)',&
              &Eband(iSpin)-TS(iSpin),"H",&
              &Hartree__eV*(Eband(iSpin)-TS(iSpin)),'eV'
          write(fdUser,format2U)'Extrapolated E(0K)',E0(iSpin),"H",&
              &Hartree__eV*(E0(iSpin)),'eV'
          if (tPrintMulliken) then
            if (nSpin == 2) then
              write(fdUser, &
                  & "(' Input/Output electrons (',A,'):',F16.8,F16.8)") &
                  &trim(spinName(iSpin)), sum(qInput(:, :, iSpin)), &
                  & sum(qOutput(:, :,iSpin))
            else
              write(fdUser, &
                  & "(' Input/Output electrons (',A,'):',F16.8,F16.8)") &
                  & quaternionName(iSpin), sum(qInput(:, :, iSpin)), &
                  & sum(qOutput(:, :,iSpin))
            end if
          end if
          write (fdUser,*) ''

        end do lpSpinPrint3
        if (nSpin == 2) then
          call ud2qm(qOutput)
          call ud2qm(qInput)
        end if

        write(fdUser,format2U) 'Energy H0', energy%EnonSCC,'H',  &
            & energy%EnonSCC*Hartree__eV,'eV'

        if (tSCC) then
          write (fdUser,format2U) 'Energy SCC', energy%ESCC,'H', &
              & energy%ESCC*Hartree__eV,'eV'
          if (tSpin) then
            write (fdUser,format2U) 'Energy SPIN', energy%Espin,'H', &
                & energy%Espin*Hartree__eV,'eV'
          end if
          if (tDFTBU) then
            write (fdUser,format2U) 'Energy DFTB+U', energy%Edftbu,'H', &
                &energy%Edftbu*Hartree__eV,'eV'
          end if
        end if
        if (tSpinOrbit) then
          write(fdUser,format2U) 'Energy L.S', energy%ELS,'H', &
              & energy%ELS*Hartree__eV,'eV'
        end if
        if (tEfield) then
          write(fdUser,format2U) 'Energy ext. field', energy%Eext,'H', &
              & energy%Eext*Hartree__eV,'eV'
        end if

        write(fdUser,format2U) 'Total Electronic energy', energy%Eelec,'H', &
            & energy%Eelec*Hartree__eV,'eV'

        !! Write out repulsive related data

        write(fdUser,format2U) 'Repulsive energy', energy%Erep,'H', &
            & energy%Erep*Hartree__eV,'eV'
        if (tDispersion) then
          write (fdUser,format2U) 'Dispersion energy', energy%eDisp,'H', &
              &energy%eDisp * Hartree__eV,'eV'
        end if
        write(fdUser,format2U) 'Total energy', energy%Etotal,'H', &
            &(energy%Etotal) * Hartree__eV,'eV'
        write(fdUser,format2U) 'Total Mermin free energy', energy%Etotal -&
            & sum(TS),'H', (energy%Etotal - sum(TS)) * Hartree__eV,'eV'
        if (tPeriodic.and.pressure/=0.0_dp) then
          write(fdUser,format2U) 'Gibbs free energy', energy%Etotal -&
              & sum(TS)+CellVol*pressure,'H', Hartree__eV * &
              & (energy%Etotal - sum(TS)+CellVol*pressure),'eV'
        end if
        write (fdUser,*) ''

        if (tAtomicEnergy) then
          write (fdUser, "(' Atom resolved electronic energies ')")
          do ii = 1, nAtom
            write (fdUser, "(1X,I5,F16.8,' H',F16.6,' eV')") ii, &
                & energy%atomElec(ii), Hartree__eV*energy%atomElec(ii)
          end do
          write (fdUser,*) ''
        end if

        if (tAtomicEnergy) then
          write (fdUser, "(' Atom resolved repulsive energies ')")
          do ii = 1, nAtom
            write (fdUser, "(1X,I5,F16.8,' H',F16.6,' eV')") ii, &
                & energy%atomRep(ii), Hartree__eV*energy%atomRep(ii)
          end do
          write (fdUser,*) ''
          write (fdUser, "(' Atom resolved total energies ')")
          do ii = 1, nAtom
            write (fdUser, "(1X,I5,F16.8,' H',F16.6,' eV')") ii, &
                & energy%atomTotal(ii), Hartree__eV*energy%atomTotal(ii)
          end do
          write (fdUser,*) ''
        end if

      end if

      if (tConverged) then
        exit lpSCC
      end if

      iSCCIter = iSCCIter + 1

    end do lpSCC

    !! Linear response
    energy%Eexcited = 0.0_dp
    excitedDerivs = 0.0_dp
    if (tLinResp) then
      ASSERT(.not. t3rd .and. tRealHS)
      dqAtom = sum( qOutput(:,:,1) - q0(:,:,1) , dim=1)
      call unpackHS(SSqrReal, over, neighborList%iNeighbor, nNeighbor,&
          & iAtomStart, iPair, img2CentCell)
      call blockSymmetrizeHS(SSqrReal, iAtomStart)
      if (tForces) then
        ASSERT(.not. tPeriodic)
        do iSpin = 1, nSpin
          call blockSymmetrizeHS(rhoSqrReal(:,:,iSpin), iAtomStart)
        end do
      end if
      if (tWriteTagged) then
        open(fdTagged, file=taggedOut, position="append")
      end if
      
      if (tLinRespZVect) then
        if (tPrintExcitedEigVecs) then
          
          call addGradients(tSpin, lresp, iAtomStart, &
              & HSqrReal, eigen(:,1,:), SSqrReal, filling(:,1,:), coord0, &
              & dqAtom, species0, neighborList%iNeighbor, &
              & img2CentCell, orb, skHamCont, skOverCont, tWriteTagged, &
              & fdTagged, energy%Eexcited, tForces, excitedDerivs, &
              & nonSccDeriv, rhoSqrReal, occNatural=occNatural(:,1), &
              & naturalOrbs=naturalOrbs(:,:,1))
          
          call writeEigvecs(fdEigvec, runId, nAtom, nSpin, neighborList, &
              & nNeighbor, iAtomStart, iPair, img2CentCell, orb, species, &
              & speciesName, over, naturalOrbs(:,:,1:1), SSqrReal, &
              & fileName="excitedOrbs")
          
        else
          call addGradients(tSpin, lresp, iAtomStart, &
              & HSqrReal, eigen(:,1,:), SSqrReal, filling(:,1,:), coord0, &
              & dqAtom, species0, neighborList%iNeighbor, &
              & img2CentCell, orb, skHamCont, skOverCont, tWriteTagged, &
              & fdTagged, energy%Eexcited, tForces, excitedDerivs, &
              & nonSccDeriv, rhoSqrReal)
        end if
      else
        call calcExcitations(tSpin, lresp, iAtomStart,&
            & HSqrReal, eigen(:,1,:), SSqrReal, filling(:,1,:), coord0,&
            & dqAtom, species0, neighborList%iNeighbor,&
            & img2CentCell, orb, tWriteTagged, fdTagged, energy%Eexcited)
      end if
      energy%Etotal = energy%Etotal + energy%Eexcited
      if (tWriteTagged) then
        close(fdTagged)
      end if
    end if

    if (tXlbomd) then
      if (xlbomdIntegrator%needsInverseJacobian()) then
        write(*, "(A)") ">> Updating XLBOMD Inverse Jacobian"
        ALLOCATE_(invJacobian, (nIneqOrb, nIneqOrb))
        call getInverseJacobian(pChrgMixer, invJacobian)
        call xlbomdIntegrator%setInverseJacobian(invJacobian)
        DEALLOCATE_(invJacobian)
      end if

      call xlbomdIntegrator%getNextCharges(qOutRed(1:nIneqOrb), &
          & qInpRed(1:nIneqOrb))
      call OrbitalEquiv_expand(qInpRed(1:nIneqOrb), iEqOrbitals, orb, qInput)
      if (tDFTBU) then
        qBlockIn = 0.0_dp
        call Block_expand(qInpRed ,iEqBlockDFTBU, orb, qBlockIn, species0, nUJ, &
            & niUJ, iUJ, orbEquiv=iEqOrbitals)
        if (tSpinOrbit) then
          call Block_expand(qInpRed ,iEqBlockDFTBULS, orb, qiBlockIn, species0, &
              & nUJ, niUJ, iUJ, skew=.true.)
        end if
      end if
      if (nSpin == 2) then
        call ud2qm(qInput)
        if (tDFTBU) then
          call ud2qm(qBlockIn)
        end if
      end if
    end if

    if (tPrintEigVecs) then
      if (tRealHS) then
        if (nSpin == 4) then
          if (tStoreEigvecs) then
            call writeEigvecs(fdEigvec, runId, nAtom, nSpin, neighborList,&
                & nNeighbor, cellVec, iCellVec, iAtomStart, iPair, &
                & img2CentCell, orb, species, speciesName, over, &
                & reshape([0.0_dp,0.0_dp,0.0_dp],(/3,1/)), HSqrCplx, SSqrCplx,&
                & storeEigvecsCplx)
          else
            call writeEigvecs(fdEigvec, runId, nAtom, nSpin, neighborList,&
                & nNeighbor, cellVec, iCellVec, iAtomStart, iPair, &
                & img2CentCell, orb, species, speciesName, over, &
                & reshape([0.0_dp,0.0_dp,0.0_dp],(/3,1/)), HSqrCplx, SSqrCplx)
          end if
        else
          if (tStoreEigvecs) then
            call writeEigvecs(fdEigvec, runId, nAtom, nSpin, neighborList, &
                &nNeighbor, iAtomStart, iPair, img2CentCell, orb, species, &
                &speciesName, over, HSqrReal, SSqrReal, storeEigvecsReal)
          else
            call writeEigvecs(fdEigvec, runId, nAtom, nSpin, neighborList, &
                &nNeighbor, iAtomStart, iPair, img2CentCell, orb, species, &
                &speciesName, over, HSqrReal, SSqrReal)
          end if
        end if
      else
        if (tStoreEigvecs) then
          call writeEigvecs(fdEigvec, runId, nAtom, nSpin, neighborList, &
              &nNeighbor, cellVec, iCellVec, iAtomStart, iPair, img2CentCell,&
              &orb, species, speciesName, over, kpoint, HSqrCplx, SSqrCplx, &
              &storeEigvecsCplx)
        else
          call writeEigvecs(fdEigvec, runId, nAtom, nSpin, neighborList,&
              &nNeighbor, cellVec, iCellVec, iAtomStart, iPair, img2CentCell,&
              &orb, species, speciesName, over, kpoint, HSqrCplx, SSqrCplx)
        end if
      end if
    end if




    if (tProjEigenvecs) then
      if (.not.tRealHS .or. (nSpin == 4)) then
        if (tStoreEigvecs) then
          call writeProjEigvecs(regionLabels, fdProjEig, eigen, nSpin, neighborList, &
              & nNeighbor, cellVec, iCellVec, iAtomStart, iPair, &
              & img2CentCell, orb, over, kpoint, kWeight, HSqrCplx, &
              & SSqrCplx, iOrbRegion, storeEigvecsCplx)
        else
          call writeProjEigvecs(regionLabels, fdProjEig, eigen, nSpin, neighborList, &
              & nNeighbor, cellVec, iCellVec, iAtomStart, iPair, &
              & img2CentCell, orb, over, kpoint, kWeight, HSqrCplx, &
              & SSqrCplx, iOrbRegion)
        end if
      else
        if (tStoreEigvecs) then
          call writeProjEigvecs(regionLabels, fdProjEig, eigen, nSpin, neighborList, &
              & nNeighbor, iAtomStart, iPair, img2CentCell, orb, over, &
              & HSqrReal, SSqrReal, iOrbRegion, storeEigvecsReal)
        else
          call writeProjEigvecs(regionLabels, fdProjEig, eigen, nSpin, neighborList, &
              & nNeighbor, iAtomStart, iPair, img2CentCell, orb, over, &
              & HSqrReal, SSqrReal, iOrbRegion)
        end if
      end if
    end if
    
    if (tGeoOpt .or. tMD) then
      if (iGeoStep == 0) then
        write (lcTmp, "(A,A)") trim(geoOutFile), ".gen"
        call clearFile(trim(lcTmp))
        write (lcTmp, "(A,A)") trim(geoOutFile), ".xyz"
        call clearFile(trim(lcTmp))
      end if
      write (lcTmp, "(A,A)") trim(geoOutFile), ".xyz"
    end if
    
    if (tGeoOpt) then      
      if (.not. tAppendGeo) then
        call clearFile(trim(lcTmp))
      end if
      if (tLatOpt) then
        write (tmpStr, "('** Geometry step: ',I0,', Lattice step: ',I0)") &
          & iGeoStep,iLatGeoStep
      else
        write(tmpStr,"('Geometry Step: ',i0)")iGeoStep
      end if
      ! save geometry in gen format
      call writeGenGeometry()
      
      if (tPrintMulliken) then
        if (nSpin == 4) then
          ALLOCATE_(tmpMatrix,(3,nAtom))
          do jj = 1, nAtom
            do ii = 1, 3
              tmpMatrix(ii,jj) = sum(qOutput(:,jj,ii+1))
            end do
          end do
          ! convert by the inverse of the scaling used in writeXYZFormat :
          tmpMatrix(:,:) = tmpMatrix(:,:) * au__fs / (1000_dp * Bohr__AA)
          call writeXYZFormat(trim(lcTmp), pCoord0Out, species0, speciesName, &
              &charges=sum(qOutput(:,:,1),dim=1), velocities = tmpMatrix, &
              & comment=trim(tmpStr))
          DEALLOCATE_(tmpMatrix)
        else
          call writeXYZFormat(trim(lcTmp), pCoord0Out, species0, speciesName, &
              &charges=sum(qOutput(:,:,1),dim=1),comment=trim(tmpStr))
        end if
      else
        call writeXYZFormat(trim(lcTmp), pCoord0Out, species0, speciesName, &
            &comment=trim(tmpStr))
      end if
    end if

    write (*,*)
    write (*, format2U) "Total Energy", energy%Etotal,"H", &
        & Hartree__eV*energy%Etotal,"eV"
    write (*, format2U) "Total Mermin free energy", &
        & energy%Etotal - sum(TS),"H", &
        & Hartree__eV*(energy%Etotal - sum(TS)),"eV"

    if (tDipole) then
      dipoleMoment(:) = 0.0_dp
      do iAtom = 1, nAtom
        dipoleMoment(:) = dipoleMoment(:) &
            & + sum(q0(:, iAtom, 1) - qOutput(:, iAtom, 1)) * coord(:,iAtom)
      end do
#if DEBUG >= 1
      ! extra test for the potential in the code, does the dipole from
      ! charge positions match the derivative of energy wrt an external E field?
      ALLOCATE_(hprime,(size(h0),1))
      ALLOCATE_(dipoleTmp,(size(qOutput,dim=1),nAtom))
      ALLOCATE_(potentialDerivative,(nAtom,1))
      write(*,"(A)",advance='no')'Hellmann Feynman dipole:'
      do ii = 1, 3 ! loop over directions
        potentialDerivative = 0.0_dp
        potentialDerivative(:,1) = -coord(ii,:) ! Potential from dH/dE
        hprime = 0.0_dp
        dipoleTmp = 0.0_dp
        call add_shift(hprime,over,nNeighbor, neighborList%iNeighbor, &
            & species,orb,iPair,nAtom,img2CentCell,potentialDerivative)
        ! evaluate <psi| dH/dE | psi>
        call mulliken(dipoleTmp, hprime(:,1), rhoPrim(:,1), &
            &orb, neighborList%iNeighbor, nNeighbor, img2CentCell, iPair)
        ! add nuclei term for derivative wrt E
        do iAtom = 1, nAtom
          dipoleTmp(1,iAtom) = dipoleTmp(1,iAtom) &
              & + sum(q0(:,iAtom,1))*coord0(ii,iAtom)
        end do
        write(*,"(f12.8)",advance='no')sum(dipoleTmp)
      end do
      write(*,*)" au"
      DEALLOCATE_(potentialDerivative)
      DEALLOCATE_(hprime)
      DEALLOCATE_(dipoleTmp)
#endif
    else
      dipoleMoment(:) = 0.0_dp
    end if

    !! Calculate energy weighted density matrix
    if (tForces) then

      !if (tXLBOMD) then
      !  if (nSpin == 4 .or. tDFTBU) then
      !    call warning("XLBOMD does not work correctly with noncollinear spin&
      !        & or DFTB+U so far. Your forces will be incorrect !!!")
      !
      !  end if
      !  ! Rebuild Hamiltonian
      !  ! (ONLY WORKS FOR SCC AND SPIN BUT NOT FOR DFTB+U as the only
      !  ! those shift vectors are updated from qOutput)
      !  ham(:,:) = 0.0_dp
      !  ham(:,1) = h0
      !  call add_shift(ham, over, nNeighbor, neighborList%iNeighbor, &
      !      & species, orb, iPair, nAtom, img2CentCell, &
      !      & potential%intBlock + potential%extBlock)
      !end if

      !! Calculate the identity part of the energy weighted density matrix
      ERhoPrim(:) = 0.0_dp

      if (nSpin == 4) then
        do nK = 1, nKPoint
          ! get eigenvectors out o storage if neccessary
          if (tStoreEigvecs) then
            iK2 = 1
            call get(storeEigvecsCplx(1), HSqrCplx(:,:,iK2, 1))
          else
            iK2 = nK
          end if

          if (tDensON2) then
            call error("Currently missing.")
          else
            call makeDensityMatrix(SSqrCplx, HSqrCplx(:,:,iK2,1), &
                &filling(:,nk,1), eigen(:,nk,1))
          end if
          if (tRealHS) then
            call packERho(ERhoPrim(:), SSqrCplx, neighborList%iNeighbor, &
                &nNeighbor, orb%mOrb, iAtomStart, iPair, img2CentCell)
          else
            call packERho(ERhoPrim(:), SSqrCplx, kPoint(:,nk), &
                &kWeight(nk), neighborList%iNeighbor, nNeighbor, orb%mOrb, &
                &iCellVec, cellVec, iAtomStart, iPair, img2CentCell)
          end if
        end do
      else
        if (tRealHS) then
          do iSpin = 1, nSpin

            if (tStoreEigvecs) then
              iSpin2 = 1
              call get(storeEigvecsReal(mod(iSpin,3)), &
                  & HSqrReal(:,:,mod(iSpin2,3)))
            else
              iSpin2 = iSpin
            end if

            ! Build energy weighted density matrix
            if (tDensON2) then
              call makeDensityMatrix(SSqrReal, HSqrReal(:,:,iSpin2), &
                  &filling(:,1,iSpin), eigen(:,1,iSpin), &
                  & neighborlist%iNeighbor, nNeighbor, orb, iAtomStart, &
                  & img2CentCell)
            else
              select case (forceType)
              case(0)
                ! Original (nonconsistent) scheme
                call makeDensityMatrix(SSqrReal, HSqrReal(:,:,iSpin2), &
                    & filling(:,1,iSpin), eigen(:,1,iSpin))
              case(1)
                ! Recreate eigenvalues for a consistent energy weighted
                ! density matrix (yields, however, incorrect forces for XLBOMD)
                call diagonalize(HSqrReal2, SSqrReal, &
                    & eigen2(:,1,iSpin), ham(:,iSpin), over, &
                    & neighborList%iNeighbor, nNeighbor, &
                    & iAtomStart, iPair, img2CentCell, solver, 'N')
                call makeDensityMatrix(SSqrReal, HSqrReal(:,:,iSpin2), &
                    & filling(:,1,iSpin), eigen2(:,1,iSpin))
              case(2)
                !! Correct force for XLBOMD for T=0K (DHD)
                !! Eigenvectors stored in HSqrReal are overwritten
                call unpackHS(SSqrReal, ham(:,iSpin), neighborlist%iNeighbor, &
                    & nNeighbor, iAtomStart, iPair, img2CentCell)
                call blockSymmetrizeHS(SSqrReal, iAtomStart)
                call makeDensityMatrix(HSqrReal2, HSqrReal(:,:,iSpin2), &
                    &filling(:,1,iSpin))
                ! D H
                call symm(HSqrReal(:,:,iSpin2), "L", HSqrReal2, SSqrReal)
                ! (D H) D
                call symm(SSqrReal, "R", HSqrReal2, HSqrReal(:,:,iSpin2), &
                    &alpha=0.5_dp)
              case(3)
                ! Correct force for XLBOMD for T <> 0K (DHS^-1 + S^-1HD)
                ! Eigenvectors stored in HSqrReal are overwritten
                call makeDensityMatrix(SSqrReal, HSqrReal(:,:,iSpin2), &
                    &filling(:,1,iSpin))
                call unpackHS(HSqrReal2, ham(:,iSpin), &
                    & neighborlist%iNeighbor, nNeighbor, iAtomStart, iPair, &
                    & img2CentCell)
                call blocksymmetrizeHS(HSqrReal2, iAtomStart)
                call symm(HSqrReal(:,:,iSpin2), "L", SSqrReal, HSqrReal2)
                call unpackHS(SSqrReal, over, neighborlist%iNeighbor, &
                    & nNeighbor, iAtomStart, iPair, img2CentCell)
                call symmatinv(SSqrReal)
                call symm(HSqrReal2, "R", SSqrReal, HSqrReal(:,:,iSpin2), &
                    & alpha=0.5_dp)
                SSqrReal = HSqrReal2 + transpose(HSqrReal2)
              end select
            end if
            call packHS(ERhoPrim, SSqrReal, neighborList%iNeighbor, &
                &nNeighbor, orb%mOrb, iAtomStart, iPair, img2CentCell)
          end do
        else
          do iSpin = 1, nSpin
            if (tStoreEigvecs) then
              iSpin2 = 1
            else
              iSpin2 = iSpin
            end if

            do nK = 1, nKPoint
              !! Calculate eigenvectors, if necessary. Eigenvectors for the last
              !! spin in the last k-points are still there, so use those
              !! directly
              if (tStoreEigvecs) then
                iK2 = 1
                call get(storeEigvecsCplx(iSpin), HSqrCplx(:,:,iK2, iSpin2))
              else
                iK2 = nK
              end if

              if (tDensON2) then
                call makeDensityMatrix(SSqrCplx, HSqrCplx(:,:,iK2,iSpin2), &
                    &filling(:,nK,iSpin), eigen(:,nK, iSpin), &
                    &neighborlist%iNeighbor, nNeighbor, orb, iAtomStart,&
                    &img2CentCell)
              else
                select case (forceType)
                case(0)
                  ! Original (nonconsistent) scheme
                  call makeDensityMatrix(SSqrCplx, HSqrCplx(:,:,iK2,iSpin2), &
                      &filling(:,nK,iSpin), eigen(:,nK, iSpin))
                case (1)
                  call error("Force type 1 not implemented for complex H")
                case(2)
                  !! Correct force for XLBOMD for T=0K (DHD)
                  !! Eigenvectors stored in HSqrCplx are overwritten
                  call makeDensityMatrix(HSqrCplx2, HSqrCplx(:,:,iK2,iSpin2), &
                      &filling(:,nK,iSpin))
                  call unpackHS(SSqrCplx, ham(:,iSpin), kPoint(:,nK), &
                      & neighborlist%iNeighbor, nNeighbor, iCellVec, cellVec, &
                      & iAtomStart, iPair, img2CentCell)
                  call blockHermitianHS(SSqrCplx, iAtomStart)
                  call hemm(HSqrCplx(:,:,iK2,iSpin2), "L", HSqrCplx2, SSqrCplx)
                  call hemm(SSqrCplx, "R", HSqrCplx2, &
                      & HSqrCplx(:,:,iK2,iSpin2), alpha=(0.5_dp, 0.0_dp))
                case(3)
                  ! Correct force for XLBOMD for T <> 0K (DHS^-1 + S^-1HD)
                  ! Eigenvectors stored in HSqrReal are overwritten
                  call makeDensityMatrix(SSqrCplx, HSqrCplx(:,:,iK2,iSpin2), &
                      &filling(:,nK,iSpin))
                  call unpackHS(HSqrCplx2, ham(:,iSpin), kPoint(:,nK), &
                      & neighborlist%iNeighbor, nNeighbor, iCellVec, cellVec, &
                      & iAtomStart, iPair, img2CentCell)
                  call blockHermitianHS(HSqrCplx2, iAtomStart)
                  call hemm(HSqrCplx(:,:,iK2,iSpin2), "L", SSqrCplx, HSqrCplx2)
                  call unpackHS(SSqrCplx, over, kPoint(:,nK), &
                      & neighborlist%iNeighbor, nNeighbor, iCellVec, cellVec, &
                      & iAtomStart, iPair, img2CentCell)
                  call hermatinv(SSqrCplx)
                  call hemm(HSqrCplx2, "R", SSqrCplx, &
                      & HSqrCplx(:,:,iK2,iSpin2), alpha=(0.5_dp, 0.0_dp))
                  SSqrCplx = HSqrCplx2 + transpose(conjg(HSqrCplx2))
                end select
              end if
              call packHS(ERhoPrim(:), SSqrCplx, kPoint(:,nk), &
                  &kWeight(nk), neighborList%iNeighbor, nNeighbor, orb%mOrb, &
                  &iCellVec, cellVec, iAtomStart, iPair, img2CentCell)
            end do
          end do
        end if
      end if

      derivs(:,:) = 0.0_dp
      if (tExtChrg) then
        chrgForces(:,:) = 0.0_dp
      end if

      if (.not. (tSCC.or.tEField)) then ! no external or internal potentials
        if (tImHam) then
          call derivative_shift(derivs, nonSccDeriv, rhoPrim, iRhoPrim,&
              & erhoPrim, skHamCont, skOverCont, coord, species, &
              & neighborList%iNeighbor, nNeighbor, img2CentCell, iPair, orb,&
              & potential%intBlock, potential%iorbitalBlock)
        else
          call derivative_nonscc(derivs, nonSccDeriv, rhoPrim(:,1), ERhoPrim,&
              & skHamCont, skOverCont, coord, species, neighborList%iNeighbor,&
              & nNeighbor, img2CentCell, iPair, orb)
        end if
      else
        if (tSCC) then
          potential%intBlock = potential%intBlock + potential%extBlock
        else
          potential%intBlock = potential%extBlock
        end if

        if (tDFTBU) then
          potential%intBlock = potential%orbitalBlock + potential%intBlock
        end if

        if (tImHam) then
          call derivative_shift(derivs, nonSccDeriv, rhoPrim, iRhoPrim,&
              & erhoPrim, skHamCont, skOverCont, coord, species,&
              & neighborList%iNeighbor, nNeighbor, img2CentCell, iPair, orb,&
              & potential%intBlock, potential%iorbitalBlock)
        else
          call derivative_shift(derivs, nonSccDeriv, rhoPrim, erhoPrim,&
              & skHamCont, skOverCont, coord, species, neighborList%iNeighbor,&
              & nNeighbor, img2CentCell, iPair, orb, potential%intBlock)
        end if

        ! add double counting terms in force :
        if (tExtChrg) then
          if (tXlbomd) then
            call error("XLBOMD does not work with external charges yet!")
            !call addForceDCSCC_XLBOMD(species, orb, neighborList%iNeighbor, &
            !    & img2CentCell, coord, qOutput, q0, derivs)
          else
            call addForceDCSCC(derivs, species, neighborList%iNeighbor, &
                & img2CentCell, coord, chrgForces)
          end if
        elseif (tSCC) then
          if (tXlbomd) then
            call addForceDCSCC_Xlbomd(species, orb, neighborList%iNeighbor, &
                & img2CentCell, coord, qOutput, q0, derivs)
          else
            call addForceDCSCC(derivs, species, neighborList%iNeighbor, &
                & img2CentCell, coord)
          end if
        end if
        if (tEField) then
          do ii = 1, 3
            derivs(ii,:) = derivs(ii,:) - &
                & sum(q0(:, :, 1)-qOutput(:,:,1),dim=1)*EField(ii)
          end do
        end if
        if (t3rdFull) then
          if (tXlbomd) then
            call thirdOrd%addGradientDcXlbomd(neighborList, species, coord, &
                & img2CentCell, qOutput, q0, orb, derivs)
          else
            call thirdOrd%addGradientDc(neighborList, species, coord, &
                & img2CentCell, derivs)
          end if
        end if
      end if

      call getERepDeriv(repulsiveDerivs, coord, nNeighbor, &
          &neighborList%iNeighbor, species,pRepCont, img2CentCell)

      totalDeriv(:,:) = repulsiveDerivs(:,:) + derivs(:,:)

      if (tLinResp) then
        totalDeriv(:,:) = totalDeriv(:,:) + excitedDerivs(:,:)
      end if

      if (tDispersion) then
        call dispersion%addGradients(totalDeriv)
      end if

      if (tStress) then
        repulsiveStress = 0.0_dp
        elecStress = 0.0_dp
        dispStress = 0.0_dp
        call getRepulsiveStress(repulsiveStress, coord, nNeighbor, &
            & neighborList%iNeighbor, species, img2CentCell, pRepCont, CellVol)
        if (tSCC) then
          if (tImHam) then
            call getBlockiStress(elecStress, nonSccDeriv, rhoPrim, iRhoPrim,&
                & ERhoPrim, skHamCont, skOverCont, coord, species,&
                & neighborList%iNeighbor, nNeighbor, img2CentCell, iPair, orb,&
                & potential%intBlock, potential%iorbitalBlock, CellVol)
          else
            call getBlockStress(elecStress, nonSccDeriv, rhoPrim, ERhoPrim,&
                & skHamCont, skOverCont, coord, species, neighborList%iNeighbor,&
                & nNeighbor, img2CentCell, iPair, orb, potential%intBlock,&
                & cellVol)
          end if

          call addStressDCSCC(elecStress,species,neighborList%iNeighbor, &
              & img2CentCell,coord)
          
        else
          if (tImHam) then
            call getBlockiStress(elecStress, nonSccDeriv, rhoPrim, iRhoPrim,&
                & ERhoPrim, skHamCont, skOverCont, coord, species, &
                & neighborList%iNeighbor, nNeighbor, img2CentCell, iPair, orb,&
                & potential%intBlock, potential%iorbitalBlock, cellVol)
          else
            call getNonSCCStress(elecStress, nonSccDeriv, rhoPrim(:,1),&
                & ERhoPrim, skHamCont, skOverCont, coord, species,&
                & neighborList%iNeighbor, nNeighbor, img2CentCell, iPair, orb,&
                & cellVol)
          end if
        end if

        if (tDispersion) then
          call dispersion%getStress(dispStress)
          dispLatDeriv = -CellVol * matmul(dispStress,invLatVec)
        end if
        
        if (tEField) then
          elecLatDeriv = 0.0_dp
          call cart2frac(coord0,latVec)
          do iAtom = 1, nAtom
            do ii = 1, 3
              do jj = 1, 3
                elecLatDeriv(jj,ii) =  elecLatDeriv(jj,ii) - &
                    & sum(q0(:, iAtom, 1)-qOutput(:,iAtom,1),dim=1) &
                    & * EField(ii) * coord0(jj,iAtom)
              end do
            end do
          end do
          call frac2cart(coord0,latVec)
          elecStress = elecStress &
              & -matmul(elecLatDeriv,transpose(latVec))/CellVol
        end if
        
        totalStress = repulsiveStress + elecStress + dispStress 
        
        cellPressure = ( totalStress(1,1) + totalStress(2,2) &
            & + totalStress(3,3) )/3.0_dp
        
        repulsiveLatDeriv = -CellVol * matmul(repulsiveStress,invLatVec)
        elecLatDeriv = -CellVol * matmul(elecStress,invLatVec)        
        totalLatDeriv = repulsiveLatDeriv + elecLatDeriv + dispLatDeriv

        write(*,format2Ue)'Volume',CellVol,'au^3',(Bohr__AA**3)*CellVol,'A^3'
        
      end if
      
    end if
    
    ! MD case includes atomic kinetic energy contribution, so print that later
    if (tStress .and. .not. tMD) then
      write(*,format2Ue)'Pressure',cellPressure,'au',&
          & cellPressure* au__pascal, 'Pa'
      if (pressure/=0.0_dp) then
        write (*, format2U) "Gibbs free energy", &
            & energy%Etotal - sum(TS(:)) + CellVol*pressure,'H', &
            & Hartree__eV*(energy%Etotal - sum(TS(:)) + CellVol*pressure), 'eV'
      end if
    end if

    !! Write out information after the end of the SCC loop, but within the
    !! geometry loop

    if (tWriteDetailedOut) then
      if (tSCC) then
        if (tConverged) then
          write (fdUser,*) "SCC converged"
          write (fdUser,*) ""
        else
          if (.not. tXlbomd) then
            write (fdUser,*)"SCC is NOT converged, maximal SCC&
                & iterations exceeded"
            write (fdUser,*) ""
            if (tConvrgForces) then
              call error("SCC is NOT converged, maximal SCC iterations&
                  & exceeded")
            else
              call warning("SCC is NOT converged, maximal SCC iterations&
                  & exceeded")
            end if
          end if
        end if
      else
        write (fdUser,*) "Non-SCC calculation"
        write (fdUser,*) ""
      end if

      if (tLinResp.and.energy%Eexcited /= 0.0_dp) then
        write (fdUser, format2U) "Excitation Energy", energy%Eexcited,"H", &
            & Hartree__eV*energy%Eexcited,"eV"
        write (fdUser,*)
      end if

      if (tGeoOpt .or. tMD) then
        write(fdUser,*) "Full geometry written in ",trim(geoOutFile), &
            & ".{xyz|gen}"
        write(fdUser,*)''
      end if

      !! Write out forces
      if (tPrintForces) then

        write(fdUser,*)'Total Forces'
        do ii = 1, nAtom
          write(fdUser,*) -totalDeriv(:,ii)
        end do
        write(fdUser,*)''

        if (tStress.and. .not.tMD) then
          write(fdUser,*)'Total stress tensor'
          do ii = 1, 3
            write(fdUser,"(3F20.12)")totalStress(:,ii)
          end do
          write(fdUser,*)
          write(fdUser,*)'Total lattice derivs'
          do ii = 1, 3
            write(fdUser,"(3F20.12)")totalLatDeriv(:,ii)
          end do
          write(fdUser,*)''
        end if

        write(fdUser,format1Ue) "Maximal derivative component", &
            & maxval(abs(totalderiv)),'au'
        if (nMovedAtom > 0) then
          write (fdUser,format1Ue) "Max force for moved atoms:", &
              &maxval(abs(totalDeriv(:,indMovedAtom))),'au'
        end if
        write(fdUser,*)''

        if (tExtChrg) then
          write (fdUser,*) "Forces on external charges"
          do ii = 1, nExtChrg
            write (fdUser, *) -chrgForces(:,ii)
          end do
          write(fdUser,*)''
        end if

        if (tPeriodic .and. .not. tMD) then
          write(fdUser,format1Ue)'Volume',CellVol,'au^3'
          if (tStress) then
            write(fdUser,format2Ue)'Pressure',cellPressure, 'au', &
                & cellPressure * au__pascal, 'Pa'
          end if
          write(fdUser,*)''
        end if

      end if
    end if

    if (tForces) then
      !! Set force components along constraint vectors zero
      do ii = 1, nGeoConstr
        iAtom = conAtom(ii)
        totalDeriv(:,iAtom) = totalDeriv(:,iAtom) &
            &- conVec(:,ii) * dot_product(conVec(:,ii), totalDeriv(:,iAtom))
      end do

      if (tCoordOpt) then
        tmpDerivs(1:nMovedCoord) = &
            & reshape(totalDeriv(:,indMovedAtom),(/ nMovedCoord /))
        write (*,"(' ',A,':',T30,E20.6)") "Maximal force component", &
            &maxval(abs(tmpDerivs))
      end if

      if (tLatOpt) then
        ! Only include the derivCellVol contribution if not MD, as the barostat
        ! would otherwise take care of this, hence add it here rather than to
        ! totalLatDeriv itself
        tmpLat3Vecs(:,:) = totalLatDeriv + derivCellVol
        tmpLatVecs(1:9) = reshape(tmpLat3Vecs, [9])

        if (tLatOptFixAng) then ! project forces to be along original lattice
          tmpLat3Vecs = tmpLat3Vecs * normOrigLatVec
          tmpLatVecs(:) = 0.0_dp
          if (any(tLatOptFixLen)) then
            tmpLatVecs(:) = 0.0_dp
            do ii = 1, 3
              if (.not.tLatOptFixLen(ii)) then
                tmpLatVecs(ii) = sum(tmpLat3Vecs(:,ii))
              end if
            end do
          else
            tmpLatVecs(1:3) = sum(tmpLat3Vecs,dim=1)
          end if
        elseif (tLatOptIsotropic) then
          tmpLat3Vecs = tmpLat3Vecs * normOrigLatVec
          tmpLatVecs(:) = 0.0_dp
          tmpLatVecs(1) = sum(tmpLat3Vecs)
        end if
        write (*,format1Ue) "Maximal Lattice force component", &
            & maxval(abs(tmpLatVecs)),'au'
      end if

      if (tSocket) then
        ! stress was computed above in the force evaluation block
        call socket%send(energy%ETotal - sum(TS), -totalDeriv, &
            & totalStress * cellVol)
      end if

      !! If geometry minimizer finished and the last calculated geometry is the
      !! minimal one (not necessary the case, depends on the optimizer!)
      !! -> we are finished.
      !! Otherwise we have to recalc everything in the converged geometry.

      if (tGeomEnd) then
        exit lpGeomOpt
      else
        if (tWriteRestart .and. tMulliken .and. tSCC .and. .not. tDerivs .and. nSCCIter > 1) then
          if (tDFTBU) then
            if (tSpinOrbit) then
              call writeQToFile(qInput, fChargeIn, orb, qBlockIn, qiBlockIn)
            else
              call writeQToFile(qInput, fChargeIn, orb, qBlockIn)
            end if
          else
            call writeQToFile(qInput, fChargeIn, orb)
          end if
          print "('>> Charges saved for restart in ',A)", fChargeIn
        end if
        if (tDerivs) then
          call next(pDerivDriver,new3Coord,totalDeriv(:,indMovedAtom), tGeomEnd)
          coord0(:,indMovedAtom) = new3Coord(:,:)
          if (tGeomEnd) exit lpGeomOpt
        elseif (tGeoOpt) then
          if (tCoordStep) then
            call next(pGeoCoordOpt, energy%Etotal - sum(TS), tmpDerivs, &
                & tmpCoords,tCoordEnd)
            if (.not.tLatOpt) tGeomEnd = tCoordEnd
          else
            call next(pGeoLatOpt, energy%Etotal - sum(TS) + CellVol*pressure, &
                & tmpLatVecs, newLatVecs,tGeomEnd)
            if (tLatOptFixAng) then ! optimization uses scaling factor of
              !  lattice vectors
              if (any(tLatOptFixLen)) then
                do ii = 3, 1, -1
                  if (.not.tLatOptFixLen(ii)) then
                    newLatVecs(3*ii-2:3*ii) =  newLatVecs(ii)*origLatVec(:,ii)
                  else
                    newLatVecs(3*ii-2:3*ii) =  origLatVec(:,ii)
                  end if
                end do
              else
                newLatVecs(7:9) =  newLatVecs(3)*origLatVec(:,3)
                newLatVecs(4:6) =  newLatVecs(2)*origLatVec(:,2)
                newLatVecs(1:3) =  newLatVecs(1)*origLatVec(:,1)
              end if
            else if (tLatOptIsotropic) then ! optimization uses scaling factor
              !  unit cell
              do ii = 3,1, -1 ! loop downwards as reusing newLatVecs
                newLatVecs(3*ii-2:3*ii) =  newLatVecs(1)*origLatVec(:,ii)
              end do
            end if
            iLatGeoStep = iLatGeoStep + 1
          end if
        elseif(tMD) then
          movedAccel(:,:) = -totalDeriv(:,indMovedAtom) / movedMass
          call next(pMDIntegrator, movedAccel ,new3Coord, movedVelo)
          if (associated(pTempProfile)) then
            call next(pTempProfile)
          end if
          call evalKE(KE, movedVelo, movedMass(1,:))
          call evalkT(pMDFrame, kT, movedVelo, movedMass(1,:))
          velocities(:, indMovedAtom) = movedVelo(:,:)
          if (tWriteRestart) then
            write(tmpStr,"('MD iter: ',i0)")iGeoStep
            ! save geometry in gen format
            call writeGenGeometry()
            if (tMulliken) then
              call writeXYZFormat(trim(lcTmp), pCoord0Out, species0, &
                  &speciesName, charges=sum(qOutput(:,:,1), dim=1),&
                  &velocities=velocities, comment=trim(tmpStr))
            else
              call writeXYZFormat(trim(lcTmp), pCoord0Out, species0, &
                  &speciesName, velocities=velocities, &
                  &comment=trim(tmpStr))
            end if
          end if
          
          if (tStress) then

            ! contribution from kinetic energy in MD, now that velocities for
            ! this geometry step are available
            call getKineticStress(kineticStress, mass, species0, velocities, &
                & CellVol)
            
            totalStress = totalStress + kineticStress
            cellPressure = ( totalStress(1,1) + totalStress(2,2) &
                & + totalStress(3,3) )/3.0_dp
            
            totalLatDeriv = -CellVol * matmul(totalStress,invLatVec)
            
          end if
          
          if (tStress .and. tWriteDetailedOut .and. tPrintForces) then
            write(fdUser,*)'Total stress tensor'
            do ii = 1, 3
              write(fdUser,"(3F20.12)")totalStress(:,ii)
            end do
            write(fdUser,*)''
            write(fdUser,*)'Total lattice derivs'
            do ii = 1, 3
              write(fdUser,"(3F20.12)")totalLatDeriv(:,ii)
            end do
            write(fdUser,*)''
          end if

          if (tSetFillingTemp) then
            write(*,format2U)'Electronic Temperature:',tempElec,'H',&
                &tempElec/Boltzmann,'K'
          end if
          if (tEfield) then
            write(*,format1U1e)'External E field', absEField, 'au', &
                & absEField * au__V_m, 'V/m'
          end if
          write(*,format2U)"MD Temperature:",kT,"H",kT/Boltzmann,"K"
          !write(*,format1U)"MD Kinetic Energy", KE,"H"
          write(*,format2U)"MD Kinetic Energy", KE,"H",Hartree__eV*KE,"eV"
          write(*,format2U)"Total MD Energy",KE + energy%Etotal - sum(TS),"H", &
              & Hartree__eV*(KE + energy%Etotal - sum(TS)),"eV"
          if (tPeriodic) then
            write(*,format2Ue)'Pressure',cellPressure,'au',&
                & cellPressure* au__pascal, 'Pa'
            if (pressure/=0.0_dp) then
              write(*,format2U) 'Gibbs free energy including KE', KE + &
                  & energy%Etotal -sum(TS(:))+CellVol*pressure, &
                  & 'H', Hartree__eV * (KE+energy%Etotal - sum(TS(:)) &
                  & + CellVol*pressure),'eV'
            end if
          end if
          if (tWriteDetailedOut) then
            if (tSetFillingTemp) then
              write (fdUser, format2U)"Electronic Temperature", tempElec,'au',&
                  &tempElec * Hartree__eV,'eV'
            end if
            write(fdUser,format1U)"MD Kinetic Energy", KE,"H"
            write(fdUser,format1U)"Total MD Energy", &
                & KE + energy%Etotal - sum(TS),"H"
            if (tPeriodic) then
              write(fdUser,format2Ue)'Pressure',cellPressure,'au',&
                  & cellPressure* au__pascal, 'Pa'
              if (pressure/=0.0_dp) then
                write(fdUser,format2U) 'Gibbs free energy including KE', KE + &
                    & energy%Etotal -sum(TS(:))+CellVol*pressure, &
                    & 'H', Hartree__eV * (KE+energy%Etotal - sum(TS(:)) &
                    & + CellVol*pressure),'eV'
              end if
            end if
            write(fdUser,format2U)"MD Temperature",kT,"H",kT/Boltzmann,"K"
          end if
        end if

        if (tGeomEnd.and.tGeoOpt) then
          diffGeo = 0.0_dp
          if (tLatOpt) then
            diffGeo = max( maxval(abs(reshape(latVec,(/9/))- newLatVecs)), &
                & diffGeo)
          end if

          if (tCoordOpt) then
            diffGeo = max(maxval(abs(reshape(coord0(:,indMovedAtom), &
                & (/nMovedCoord/))- tmpCoords)),diffGeo)
          end if
          if (diffGeo < tolSameDist) then
            tGeomEnd = .true.
            exit lpGeomOpt
          end if
        end if

        if (.not. tGeomEnd .and. .not. tSocket) then
          if (tGeoOpt) then
            if (tCoordStep) then
              if (tCoordEnd) then
                diffGeo = maxval(abs(reshape(coord0(:,indMovedAtom), &
                    & (/nMovedCoord/))- tmpCoords))
                if (diffGeo < tolSameDist) then
                  tCoordStep = .false.
                  if (tLatOpt) then
                    tCoordEnd = .false.
                  end if
                end if
              end if
              if (nMovedCoord > 0) then
                ! Workaround for NAG
                !coord0(:,indMovedAtom) = reshape(tmpCoords(:), (/3, nMovedAtom&
                !    &/))
                do ii = 1, nMovedAtom
                  coord0(:,indMovedAtom(ii)) = tmpCoords((ii-1)*3+1:ii*3)
                end do
              end if
            else
              call cart2frac(coord0,latVec)
              latVec = reshape(newLatVecs, (/3,3/))
              call frac2cart(coord0,latVec)
              recVec2p = latVec(:,:)
              call matinv(recVec2p)
              recVec2p = reshape(recVec2p, (/3, 3/), order=(/2, 1/))
              recVec = 2.0_dp * pi * recVec2p
              CellVol = abs(determinant33(latVec))
              recCellVol = abs(determinant33(recVec))
              if (tSCC) then
                call updateLatVecs_SCC(latVec, recVec, CellVol)
                mCutoff = max(mCutoff, getSCCCutoff())
              end if
              if (tDispersion) then
                call dispersion%updateLatVecs(latVec)
                mCutoff = max(mCutoff, dispersion%getRCutoff())
              end if
              call getCellTranslations(cellVec, rCellVec, latVec, recVec2p, &
                  & mCutoff)
              if (tCoordOpt) then
                tCoordStep = .true.
                tCoordEnd = .false.
                tmpCoords(1:nMovedCoord) = reshape(coord0(:, indMovedAtom), &
                    & (/ nMovedCoord /))
                call reset(pGeoCoordOpt, tmpCoords)
              end if
            end if
          elseif (tMD) then

            coord0(:,indMovedAtom) = new3Coord(:,:)

            if (tBarostat) then ! apply a Barostat
              call rescale(pMDIntegrator,coord0,latVec,totalStress)
              recVec2p = latVec(:,:)
              call matinv(recVec2p)
              recVec2p = reshape(recVec2p, (/3, 3/), order=(/2, 1/))
              recVec = 2.0_dp * pi * recVec2p
              recCellVol = abs(determinant33(recVec))
              if (tSCC) then
                call updateLatVecs_SCC(latVec, recVec, CellVol)
                mCutoff = max(mCutoff, getSCCCutoff())
              end if
              if (tDispersion) then
                call dispersion%updateLatVecs(latVec)
                mCutoff = max(mCutoff, dispersion%getRCutoff())
              end if
              call getCellTranslations(cellVec, rCellVec, latVec, recVec2p, &
                  & mCutoff)
            end if

            if (tWriteRestart) then
              if (tStress) then
                if (tBarostat) then
                  write(fdMD,*)'Lattice vectors (A)'
                  do ii = 1, 3
                    write(fdMD,*)latVec(:,ii)*Bohr__AA
                  end do
                  write(fdMD,format2Ue)'Volume',CellVol,'au^3', &
                      & (Bohr__AA**3)*CellVol,'A^3'
                end if
                write(fdMD,format2Ue)'Pressure',cellPressure,'au',&
                    & cellPressure * au__pascal, 'Pa'
                if (pressure/=0.0_dp) then
                  write(fdMD,format2U) 'Gibbs free energy', &
                      & energy%Etotal -sum(TS)+CellVol*pressure, &
                      & 'H', Hartree__eV * (energy%Etotal - sum(TS) &
                      & + CellVol*pressure),'eV'
                  write(fdMD,format2U) 'Gibbs free energy including KE', KE + &
                      & energy%Etotal -sum(TS)+CellVol*pressure, &
                      & 'H', Hartree__eV * (KE+energy%Etotal - sum(TS) &
                      & + CellVol*pressure),'eV'
                end if
              end if
              if (tLinResp .and. energy%Eexcited /= 0.0_dp) then ! need to
                ! sort out properly where this is actually evaluated
                write (fdMD, format2U) "Excitation Energy", &
                    & energy%Eexcited,"H", Hartree__eV*energy%Eexcited,"eV"
              end if
              write(fdMD,format2U)'Potential Energy', &
                  & energy%Etotal - sum(TS),'H', &
                  & (energy%Etotal - sum(TS))*Hartree__eV,'eV'
              write(fdMD,format2U)'MD Kinetic Energy',KE,'H',KE*Hartree__eV,'eV'
              write(fdMD,format2U)'Total MD Energy', &
                  & KE+energy%Etotal - sum(TS),'H', &
                  & (KE+energy%Etotal - sum(TS))*Hartree__eV,'eV'
              write(fdMD,format2U)'MD Temperature', kT,'au',kT/Boltzmann,'K'
              if (tEfield) then
                write(fdMD,format1U1e)'External E field', absEField, 'au', &
                    &absEField * au__V_m, 'V/m'
              end if
              if (tFixEf .and. tPrintMulliken) then
                write(fdMD,"(' Net charge     : ',F14.8)") &
                    & sum(q0(:, ii, 1) - qOutput(:, ii, 1))
              end if
              if (tDipole) then
                write(fdMD,"(' Dipole moment  :',3f14.8,' au')")dipoleMoment
                write(fdMD,"(' Dipole moment  :',3f14.8,' Debye')") &
                    & dipoleMoment*au__Debye
              end if
            end if
          end if
        end if
      end if

      if (tWriteDetailedOut.and.tMD) then
        write(fdUser,format1U)"MD Kinetic Energy", KE,"H"
        write(fdUser,format2U)"Total MD Energy", &
            & KE + energy%Etotal - sum(TS),"H", &
            & Hartree__eV*(KE + energy%Etotal - sum(TS)), "eV"
        write(fdUser,format2U)"MD Temperature",kT,"H",kT/Boltzmann,"K"
        write(fdUser,*) ""
      end if
    end if

    if (tWriteDetailedOut .and. tEfield) then
      write(fdUser,format1U1e)'External E field', absEField,'au', &
          & absEField * au__V_m, 'V/m'
    end if

    !! Stop reading of initial charges/block populations again
    tReadChrg = .false.

    !! Stop SCC if appropriate stop file is present
    if (.not. tStopSCC) then
      inquire(file=fStopDriver, exist=tStopDriver)
      if (tStopDriver) then
        write (*,*) "Stop file '" // fStopDriver // "' found."
      end if
    end if
    if (tStopSCC .or. tStopDriver) then
      nGeoSteps = iGeoStep
      write (*,*) "Setting max number of geometry steps to current step number."
    end if

    iGeoStep = iGeoStep + 1
  end do lpGeomOpt

  if (tSocket) then
    call socket%shutdown()
  end if

  if (tWriteDetailedOut.and.tDipole) then
    write(fdUser,"(' Dipole moment  :',3f14.8,' au')")dipoleMoment
    write(fdUser,"(' Dipole moment  :',3f14.8,' Debye')") &
        & dipoleMoment*au__Debye
    write(fdUser,*)''
  end if

  tGeomEnd = tMD .or. tGeomEnd .or. tDerivs

  if (tWriteDetailedOut) then
    if (tGeoOpt) then
      if (tGeomEnd) then
        write (fdUser,*) "Geometry converged"
      else
        write (fdUser, *) "!!! Geometry did NOT converge!"
      end if
    elseif (tMD) then
      if (tGeomEnd) then
        write (fdUser,*) "Molecular dynamics completed"
      else
        write (fdUser, *) "!!! Molecular dynamics terminated abnormally!"
      end if
    elseif (tDerivs) then
      if (tGeomEnd) then
        write (fdUser,*) "Second derivatives completed"
      else
        write (fdUser, *) "!!! Second derivatives terminated abnormally!"
      end if
    end if
    write(fdUser,*)''
    close(fdUser)
  end if

  if (tGeoOpt) then
    if (tGeomEnd) then
      write (*,*)
      write (*,*) "Geometry converged"
    else
      call warning("!!! Geometry did NOT converge!")
    end if
  elseif (tMD) then
    if (tGeomEnd) then
      write (*,*)
      write (*,*) "Molecular dynamics completed"
    else
      call warning("!!! Molecular dynamics terminated abnormally!")
    end if
  elseif (tDerivs) then
    if (tGeomEnd) then
      write (*,*)
      write (*,*) "Second derivatives completed"
    else
      call warning("!!! Second derivatives terminated abnormally!")
    end if
  end if

  if (tMD) then
    write(*,*)'MD information accumulated in ',mdOut
    close(fdMD)
  end if

  if (tDerivs) then
    call getHessianMatrix(pDerivDriver,pDynMatrix)
    write(*,*)'Hessian matrix written to ',hessianOut
    do ii = 1, size(pDynMatrix,dim=2)
      write(fdHessian,formatHessian)pDynMatrix(:,ii)
    end do
    close(fdHessian)
  end if

  if (tLocalise) then ! warning the canonical DFTB ground state
                      ! orbitals are over-written after this point :
    if (tPipekMezey) then

      if (tStoreEigvecs) then
        call error("Pipek-Mezey localisation not implemented for stored &
            &eigenvectors")
      end if

      if (nSpin > 2) then
        call error("Pipek-Mezey localisation not implemented for &
            &non-colinear DFTB")
      end if

      if (any( abs(mod(filling,real(3-nSpin,dp))) > elecTolMax)) then
        call warning("Fractional occupations present for electron &
            &localisation")
      end if

      if (tRealHS) then

        call unpackHS(SSqrReal,over,neighborList%iNeighbor, nNeighbor, &
            &iAtomStart, iPair, img2CentCell)
        do iSpin = 1, nSpin
          nFilledLev = floor(nEl(iSpin)/real(3-nSpin,dp))
          localisation = PipekMezeyLocalisation(HSqrReal(:,1:nFilledLev,iSpin),&
              & SSqrReal,iAtomStart)

          write(*,*)'Original localisation',localisation

          if (tPipekDense) then
            call PipekMezey(HSqrReal(:,1:nFilledLev,iSpin), SSqrReal, &
                & iAtomStart,PipekTol,PipekMaxIter)
          else
            do ii = 1, size(sparsePipekTols)
              call PipekMezey(HSqrReal(:,1:nFilledLev,iSpin), SSqrReal, &
                  & iAtomStart,PipekTol,PipekMaxIter,sparsePipekTols(ii))
            end do
          end if

          localisation = PipekMezeyLocalisation(HSqrReal(:,1:nFilledLev,iSpin),&
              & SSqrReal,iAtomStart)

          write(*,*)'Final localisation ',localisation

        end do

        call writeEigvecs(fdEigvec, runId, nAtom, nSpin, neighborList, &
            & nNeighbor, iAtomStart, iPair, img2CentCell, orb, species, &
            & speciesName, over, HSqrReal, SSqrReal, fileName="localOrbs")

      else

        do iSpin = 1, nSpin

          nFilledLev = floor(nEl(iSpin)/real(3-nSpin,dp))

          localisation = sum(PipekMezeyLocalisation( &
              & HSqrCplx(:,:nFilledLev,:,iSpin), SSqrCplx, over, kpoint, &
              & kweight, neighborList%iNeighbor, nNeighbor, iCellVec, cellVec, &
              & iAtomStart, iPair, img2CentCell))

          write(*,*)'Original localisation',localisation

          call PipekMezey(HSqrCplx(:,:nFilledLev,:,iSpin), SSqrCplx, &
              & over, kpoint, kweight, neighborList%iNeighbor, nNeighbor, &
              & iCellVec, cellVec, iAtomStart, iPair, img2CentCell, PipekTol, &
              & PipekMaxIter)

          localisation = sum(PipekMezeyLocalisation( &
              & HSqrCplx(:,:nFilledLev,:,iSpin), SSqrCplx, over, kpoint, &
              & kweight, neighborList%iNeighbor, nNeighbor, iCellVec, cellVec, &
              & iAtomStart, iPair, img2CentCell))

          write(*,*)'Final localisation',localisation

        end do

        call writeEigvecs(fdEigvec, runId, nAtom, nSpin, neighborList, &
            & nNeighbor, cellVec, iCellVec, iAtomStart, iPair, img2CentCell, &
            & orb, species, speciesName, over, kpoint, HSqrCplx, SSqrCplx, &
            & fileName="localOrbs")

      end if

    end if
  end if

  if (tWriteTagged) then
    !! Write out final results in tagged format
    open(fdTagged, file=taggedOut, position="append")
    if (tPeriodic) then
      CellVol = abs(determinant33(latVec))
      call writeTagged(fdTagged, tag_volume, CellVol)
    end if
    if (tPrintMulliken) then
      call qm2ud(qOutput)
      call writeTagged(fdTagged, tag_qOutput, qOutput(:,:,1))
      call ud2qm(qOutput)
    end if
    if (tForces) then
      call writeTagged(fdTagged, tag_forceTot, -totalDeriv)
      if (tExtChrg) then
        call writeTagged(fdTagged, tag_chrgForces, -chrgForces)
      end if
      if (tLinResp) then
        call writeTagged(fdTagged, tag_excForce, -excitedDerivs)
      end if
      if (tStress) then
        call writeTagged(fdTagged, tag_stressTot, totalStress)
      end if
    end if
    if (tDerivs) then
      call writeTagged(fdTagged, tag_HessianNum, pDynMatrix)
    end if
    call writeTagged(fdTagged, tag_freeEgy, energy%Etotal - sum(TS))
    if (pressure/=0.0_dp) then
      call writeTagged(fdTagged, tag_Gibbsfree, &
          & energy%Etotal - sum(TS) + CellVol*pressure)
    end if
    call writeTagged(fdTagged, tag_endCoord, coord0)
    if (tLocalise) then
      call writeTagged(fdTagged, tag_pmlocalise, localisation)
    end if
    close(fdTagged)
  end if

  if (tWriteResultsTag) then
    open(fdResultsTag, file=resultsTag, position="rewind", status="replace")
    call writeTagged(fdResultsTag, tag_egyTotal, energy%Etotal)
    if (tAtomicEnergy) then
      call writeTagged(fdResultsTag, tag_egyTotalAt, energy%atomTotal)
    end if
    call writeTagged(fdResultsTag, tag_forces, tForces)
    if (tForces) then
      call writeTagged(fdResultsTag, tag_forceTot, -totalDeriv)
      if (tExtChrg) then
        call writeTagged(fdResultsTag, tag_chrgForces, -chrgForces)
      end if
    end if
    if (tStress) then
      call writeTagged(fdResultsTag, tag_stressTot, totalStress)
    end if
    if (tDerivs) then
      call writeTagged(fdResultsTag, tag_HessianNum, pDynMatrix)
    end if
    call writeTagged(fdResultsTag, tag_scc, tSCC)
    if (tSCC) then
      call writeTagged(fdResultsTag, tag_nSCC, iSCCIter)
      call writeTagged(fdResultsTag, tag_sccConv, tConverged)
    end if
    if (tPrintMulliken) then
      call writeTagged(fdResultsTag, tag_qOutputAt, sum(qOutput, dim=1))
      call writeTagged(fdResultsTag, tag_qOutAtNet, &
          &sum(q0(:,:,1) - qOutput(:,:,1), dim=1))
    end if
    call writeTagged(fdResultsTag, tag_eigenVal, eigen)
    call writeTagged(fdResultsTag, tag_filling, filling)
    call writeTagged(fdResultsTag, tag_efermi, Ef)
    if (size(nEl) == 1) then
      call writeTagged(fdResultsTag, tag_nElUp, 0.5_dp*nEl(1))
      call writeTagged(fdResultsTag, tag_nElDown, 0.5_dp*nEl(1))
    else
      call writeTagged(fdResultsTag, tag_nElUp, nEl(1))
      call writeTagged(fdResultsTag, tag_nElDown, nEl(2))
    end if
     if (tPeriodic) then
      CellVol = abs(determinant33(latVec))
      call writeTagged(fdResultsTag, tag_volume, CellVol)
    end if

    close(fdResultsTag)
  end if


  if (tWriteDetailedXML) then
    !! Ugly hack for printing out xml info, will be removed later
    call xml_OpenFile("detailed.xml", xf, indent=.true.)
    call xml_ADDXMLDeclaration(xf)
    call xml_NewElement(xf, "detailedout")
    call writeChildValue(xf, "identity", runId)
    call xml_NewElement(xf, "geometry")
    call writeChildValue(xf, "typenames", speciesName)
    call writeChildValue(xf, "typesandcoordinates", &
        &reshape(species0, (/ 1, size(species0) /)), pCoord0Out)
    call writeChildValue(xf, "periodic", tPeriodic)
    if (tPeriodic) then
      call writeChildValue(xf, "latticevectors", latVec)
    end if
    call xml_EndElement(xf, "geometry")
    call writeChildValue(xf, "real", tRealHS)
    call writeChildValue(xf, "nrofkpoints", nKPoint)
    call writeChildValue(xf, "nrofspins", nSpin)
    call writeChildValue(xf, "nrofstates", size(eigen, dim=1))
    call writeChildValue(xf, "nroforbitals", nOrb)
    ALLOCATE_(bufferRealR2, (4, nKPoint))
    bufferRealR2(1:3, :) = kPoint(:,:)
    bufferRealR2(4, :) = kWeight(:)
    call writeChildValue(xf, "kpointsandweights", bufferRealR2)
    DEALLOCATE_(bufferRealR2)
    call xml_NewElement(xf, "occupations")
    do ii = 1, nSpin
      call xml_NewElement(xf, "spin" // i2c(ii))
      do jj = 1, nKpoint
        call writeChildValue(xf, "k" // i2c(jj), filling(:, jj, mod(ii,3)))
      end do
      call xml_EndElement(xf, "spin" // i2c(ii))
    end do
    call xml_EndElement(xf, "occupations")
    if (tLinResp .and. tPrintExcitedEigVecs) then
      call xml_NewElement(xf, "excitedoccupations")
      call xml_NewElement(xf, "spin" // i2c(1))
      call writeChildValue(xf, "k" // i2c(1), occNatural)
      call xml_EndElement(xf, "spin" // i2c(1))
      call xml_EndElement(xf, "excitedoccupations")
    end if

    call xml_EndElement(xf, "detailedout")
    call xml_Close(xf)
  end if

  DEALLOCATE_(orbitalL)
  DEALLOCATE_(orbitalLPart)
  DEALLOCATE_(new3Coord)
  DEALLOCATE_(tmpDerivs)

  !! Deallocate arrays
  DEALLOCATE_(filling)

  if (tForces) then
    DEALLOCATE_(derivs)
    DEALLOCATE_(repulsiveDerivs)
    DEALLOCATE_(chrgForces)
  end if

  call destroy(energy)
  call destroy(potential)

  if (tMulliken) then
    DEALLOCATE_(qOutput)
    DEALLOCATE_(qInput)
    DEALLOCATE_(q0)
  end if
  DEALLOCATE_(rhoPrim)
  DEALLOCATE_(iRhoPrim)
  if (tForces) then
    DEALLOCATE_(ERhoPrim)
  end if

  if (tMD) then
    DEALLOCATE_(velocities)
    DEALLOCATE_(movedVelo)
    DEALLOCATE_(movedAccel)
    DEALLOCATE_(movedMass)
  end if

  DEALLOCATE_(HSqrCplx)
  DEALLOCATE_(SSqrCplx)
  DEALLOCATE_(HSqrReal)
  DEALLOCATE_(SSqrReal)
  DEALLOCATE_(eigen)

  if (tSCC) then
    call destruct_SCC()
    call destroy(pChrgMixer)
  end if

  call destroyProgramVariables()

contains

  ! Invokes the writing routines for the Hamiltonian and overlap matrices.
  subroutine writeHS(tWriteHS, tWriteRealHS, ham, over, iNeighbor, &
      &nNeighbor, iAtomStart, iPair, img2CentCell, kPoint, iCellVec, &
      &cellVec, iHam)
    logical, intent(in) :: tWriteHS, tWriteRealHS
    real(dp), intent(in) :: ham(:,:), over(:)
    integer, intent(in) :: iNeighbor(0:,:), nNeighbor(:)
    integer, intent(in) :: iAtomStart(:), iPair(0:,:), img2CentCell(:)
    real(dp), intent(in) :: kPoint(:,:)
    integer, intent(in) :: iCellVec(:)
    real(dp), intent(in) :: cellVec(:,:)
    real(dp), intent(in), optional :: iHam(:,:)

    integer :: iS, nSpin

    nSpin = size(ham, dim=2)

    if (tWriteRealHS) then
      do iS = 1, nSpin
        call writeSparse("hamreal" // i2c(iS) // ".dat", ham(:,iS), iNeighbor, &
            &nNeighbor, iAtomStart, iPair, img2CentCell, iCellVec, cellVec)
        if (present(iHam)) then
          call writeSparse("hamimag" // i2c(iS) // ".dat", iHam(:,iS),&
              & iNeighbor, nNeighbor, iAtomStart, iPair, img2CentCell,iCellVec,&
              & cellVec)
        end if
      end do
      call writeSparse("overreal.dat", over, iNeighbor, &
          &nNeighbor, iAtomStart, iPair, img2CentCell, iCellVec, cellVec)
    end if
    if (tWriteHS) then
      if (tRealHS) then
        do iS = 1, nSpin
          call writeSparseAsSquare("hamsqr" // i2c(iS) // ".dat", ham(:,iS), &
              &iNeighbor, nNeighbor, iAtomStart, iPair, img2CentCell)
        end do
        call writeSparseAsSquare("oversqr.dat", over, iNeighbor, nNeighbor, &
            &iAtomStart, iPair, img2CentCell)
      else
        do iS = 1, nSpin
          call writeSparseAsSquare("hamsqr" // i2c(iS) // ".dat", ham(:,iS), &
              &kPoint, iNeighbor, nNeighbor, iAtomStart, iPair, img2CentCell, &
              &iCellVec, cellVec)
        end do
        call writeSparseAsSquare("oversqr.dat", over, kPoint, iNeighbor, &
            &nNeighbor, iAtomStart, iPair, img2CentCell, iCellVec, cellVec)
      end if
    end if

  end subroutine writeHS


  !> Calculates electron fillings and resulting band energy terms.
  !!
  subroutine getFillingsAndBandEnergies(eigvals, nElectrons, nSpinBlocks, tempElec, kWeights,&
      & tSpinSharedEf, tFillKSep, tFixEf, iDistribFn, Ef, fillings, Eband, TS, E0)

    !> Eigenvalue of each level, kpoint and spin channel
    real(dp), intent(in) :: eigvals(:,:,:)

    !> Nr. of electrons for each spin channel
    real(dp), intent(in) :: nElectrons(:)

    !> Nr. of spin blocks in the Hamiltonian (1 - spin avg, 2 - colinear, 4 - non-colinear)
    integer, intent(in) :: nSpinBlocks

    !> Electronic temperature
    real(dp), intent(in) :: tempElec

    !> Weight of the k-points.
    real(dp), intent(in) :: kWeights(:)

    !> Whether for colinear spin a common Fermi level for both spin channels should be used
    logical, intent(in) :: tSpinSharedEf

    !> Whether each K-point should be filled separately (individual Fermi-level for each k-point)
    logical, intent(in) ::  tFillKSep

    !> Whether fixed Fermi level(s) should be used. (No charge conservation!)
    logical, intent(in) :: tFixEf

    !> Selector for the distribution function
    integer, intent(in) :: iDistribFn

    !> Fixed Fermi levels on entry, if tFixEf is .true., otherwise the Fermi levels found for the
    !! given number of electrons on exit
    real(dp), intent(inout) :: Ef(:)

    !> Fillings
    real(dp), intent(out) :: fillings(:,:,:)

    !> Band energies
    real(dp), intent(out) :: Eband(:)

    !> Band entropies
    real(dp), intent(out) :: TS(:)

    !> Band energies extrapolated to zero Kelvin
    real(dp), intent(out) :: E0(:)

    real(dp) :: EbandTmp(1), TSTmp(1), E0Tmp(1)
    real(dp) :: EfTmp
    real(dp) :: nElecFill(2)
    integer :: nSpinHams, nKPoints
    integer :: iS, iK

    nKPoints = size(fillings, dim=2)
    nSpinHams = size(fillings, dim=3)

    if (nSpinBlocks == 1) then
      ! Filling functions assume one electron per level, but for spin unpolarised we have two
      nElecFill(1) = nElectrons(1) / 2.0_dp
    else
      nElecFill(1:nSpinHams) = nElectrons(1:nSpinHams)
    end if

    if (tFixEf) then
      ! Fixed Fermi level for each spin channel
      do iS = 1, nSpinHams
        call electronFill(Eband(iS:iS), fillings(:,:,iS:iS), TS(iS:iS), E0(iS:iS), Ef(iS),&
            & eigvals(:,:,iS:iS), tempElec, iDistribFn, kWeights)
      end do
    else if (nSpinHams == 2 .and. tSpinSharedEf) then
      ! Common Fermi level across two colinear spin channels
      call Efilling(Eband, Ef(1), TS, E0, fillings, eigvals, sum(nElecFill), tempElec, kWeights,&
          & iDistribFn)
      Ef(2) = Ef(1)
    else if (tFillKSep) then
      ! Every spin channel and every k-point filled up individually.
      Eband(:) = 0.0_dp
      Ef(:) = 0.0_dp
      TS(:) = 0.0_dp
      E0(:) = 0.0_dp
      do iS = 1, nSpinHams
        do iK = 1, nKPoints
          call Efilling(EbandTmp, EfTmp, TSTmp, E0Tmp, fillings(:, iK:iK, iS:iS),&
              & eigvals(:, iK:iK, iS:iS), nElecFill(iS), tempElec, [1.0_dp], iDistribFn)
          Eband(iS) = Eband(iS) + EbandTmp(1) * kWeights(iK)
          Ef(iS) = Ef(iS) + EfTmp * kWeights(iK)
          TS(iS) = TS(iS) + TSTmp(1) * kWeights(iK)
          E0(iS) = E0(iS) + E0Tmp(1) * kWeights(iK)
        end do
      end do
    else
      ! Every spin channel (but no the k-points) filled up individually
      do iS = 1, nSpinHams
        call Efilling(Eband(iS:iS), Ef(iS), TS(iS:iS), E0(iS:iS), fillings(:,:,iS:iS),&
            & eigvals(:,:,iS:iS), nElecFill(iS), tempElec, kWeights, iDistribFn)
      end do
    end if

    if (nSpinBlocks == 1) then
      ! Prefactor 2 for spin unpolarised calculations
      Eband(:) = 2.0_dp * Eband
      E0(:) = 2.0_dp * E0
      TS(:) = 2.0_dp * TS
      fillings(:,:,:) = 2.0_dp * fillings
    end if

  end subroutine getFillingsAndBandEnergies

  !> Write out geometry in gen format if needed
  subroutine writeGenGeometry()
    character(lc) :: lcTmpLocal
    if (tGeoOpt .or. tMD) then
      if (tWriteRestart) then
        write (lcTmpLocal, "(A,A)") trim(geoOutFile), ".gen"
        call clearFile(trim(lcTmpLocal))
        if (tPeriodic) then
          call writeGenFormat(trim(lcTmpLocal), pCoord0Out, species0, speciesName, &
              &latVec, tFracCoord)
        else
          call writeGenFormat(trim(lcTmpLocal), coord0, species0, speciesName)
        end if
      end if
    end if
  end subroutine writeGenGeometry
  
end program dftbplus
