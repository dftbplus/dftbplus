!--------------------------------------------------------------------------------------------------!
!  DFTB+: general package for performing fast atomistic simulations                                !
!  Copyright (C) 2017  DFTB+ developers group                                                      !
!                                                                                                  !
!  See the LICENSE file for terms of usage and distribution.                                       !
!--------------------------------------------------------------------------------------------------!

#:include 'common.fypp'

!> The main routines for DFTB+
module main
  use assert
  use constants
  use io
  use inputdata_module
  use nonscc
  use eigenvects
  use repulsive
  use etemp
  use populations
  use densitymatrix
  use forces
  use stress
  use taggedoutput
  use scc
  use sccinit
  use externalcharges
  use periodic
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
  use sparse2dense
  use blasroutines, only : symm, hemm
  use hsdutils
  use charmanip
  use shift
  use spinorbit
  use angmomentum
  use elecconstraints
  use pmlocalisation, only : TPipekMezey
  use linresp_module
  use mainio
  use commontypes
  use dispersions, only : DispersionIface
  use xmlf90
  use ipisocket, only : IpiSocketComm
  use thirdorder_module, only : ThirdOrder
  use simplealgebra
  use message
  use repcont
  use xlbomd_module
  use fifo
  use slakocont
  use linkedlist
  use lapackroutines
  use mdcommon
  use mdintegrator
  use tempprofile
  implicit none
  private

  public :: runDftbPlus

  ! Name of the human readable files

  !> Tagged output files (machine readable)
  character(*), parameter :: autotestTag = "autotest.tag"

  !> Detailed user output
  character(*), parameter :: userOut = "detailed.out"

  !> band structure and filling information
  character(*), parameter :: bandOut = "band.out"

  !> File accumulating data during an MD run
  character(*), parameter :: mdOut = "md.out"

  !> Machine readable tagged output
  character(*), parameter :: resultsTag = "results.tag"

  !> Second derivative of the energy with respect to atomic positions
  character(*), parameter :: hessianOut = "hessian.out"

  !> O(N^2) density matrix creation
  logical, parameter :: tDensON2 = .false.

  !> Should further output be appended to detailed.out?
  logical, parameter :: tAppendDetailedOut = .false.


contains

  !> The main DFTB program itself
  subroutine runDftbPlus()
    use initprogram

    !> Square dense hamiltonian storage for cases with k-points
    complex(dp), allocatable :: HSqrCplx(:,:,:,:)

    !> Square dense overlap storage for cases with k-points
    complex(dp), allocatable :: SSqrCplx(:,:)

    !> Square dense hamiltonian storage
    real(dp), allocatable :: HSqrReal(:,:,:)

    !> Square dense overlap storage
    real(dp), allocatable :: SSqrReal(:,:)

    !> Eigenvalues
    real(dp), allocatable :: eigen(:,:,:)

    !> Sparse storage of density matrix
    real(dp), allocatable :: rhoPrim(:,:)

    !> Imaginary part of density matrix in sparse storage
    real(dp), allocatable :: iRhoPrim(:,:)

    !> Energy weighted density matrix
    real(dp), allocatable :: ERhoPrim(:)

    !> Non-SCC part of the hamiltonian in sparse storage
    real(dp), allocatable :: h0(:)

    !> electronic filling
    real(dp), allocatable :: filling(:,:,:)

    !> band structure energy
    real(dp), allocatable :: Eband(:)

    !> entropy of electrons at temperature T
    real(dp), allocatable :: TS(:)

    !> zero temperature electronic energy
    real(dp), allocatable :: E0(:)

    !> energy in previous scc cycles
    real(dp) :: Eold

    !> Total energy components
    type(TEnergies) :: energy

    !> Potentials for orbitals
    type(TPotentials) :: potential

    !> Energy derivative with respect to atomic positions
    real(dp), allocatable :: totalDeriv(:,:)

    !> Forces on any external charges
    real(dp), allocatable :: chrgForces(:,:)

    !> excited state force addition
    real(dp), allocatable :: excitedDerivs(:,:)

    !> Stress tensors for various contribution in periodic calculations
    real(dp) :: totalStress(3,3), totalLatDeriv(3,3)

    !> derivative of cell volume wrt to lattice vectors, needed for pV term
    real(dp) :: derivCellVol(3,3)

    !> dipole moments when available
    real(dp), allocatable :: dipoleMoment(:)

    !> whether scc converged
    logical :: tConverged

    !> internal pressure within the cell
    real(dp) :: cellPressure

    !> Geometry steps so far
    integer :: iGeoStep

    !> Lattice geometry steps so far
    integer :: iLatGeoStep

    !> Do we have the final geometry?
    logical :: tGeomEnd

    !> Has this completed?
    logical :: tCoordEnd

    !> do we take an optimization step on the lattice or the internal coordinates if optimizing both
    !> in a periodic geometry
    logical :: tCoordStep

    !> Folded coords (3, nAtom)
    real(dp), allocatable, target :: coord0Fold(:,:)

    !> Coordinates to print out
    real(dp), pointer :: pCoord0Out(:,:)

    !> New coordinates returned by the MD routines
    real(dp), allocatable :: newCoords(:,:)

    !> lattice vectors returned by the optimizer
    real(dp) :: constrLatDerivs(9)

    !> MD velocities
    real(dp), allocatable :: velocities(:,:)

    !> MD velocities for moved atoms
    real(dp), allocatable :: movedVelo(:,:)

    !> MD acceleration for moved atoms
    real(dp), allocatable :: movedAccel(:,:)

    !> Mass of the moved atoms
    real(dp), allocatable :: movedMass(:,:)

    !> MD instantaneous thermal energy
    real(dp) :: kT

    !> external electric field
    real(dp) :: Efield(3), absEfield

    !> Difference between last calculated and new geometry.
    real(dp) :: diffGeo

    !> Loop variables
    integer :: iSCCIter

    !> File descriptor for the tagged writer
    integer :: fdAutotest

    !> File descriptor for the human readable output
    integer :: fdUser

    !> File descriptor for the band structure output
    integer :: fdBand

    !> File descriptor for the eigenvector output
    integer :: fdEigvec

    !> File descriptor for detailed.tag
    integer :: fdResultsTag

    !> File descriptor for extra MD output
    integer :: fdMD

    !> File descriptor for numerical Hessian
    integer :: fdHessian

    !> Charge error in the last iterations
    real(dp) :: sccErrorQ, diffElec
    real(dp), allocatable :: orbitalL(:,:,:)

    !> Dynamical (Hessian) matrix
    real(dp), pointer :: pDynMatrix(:,:)

    !> flag to write out geometries (and charge data if scc) when moving atoms about - in the case
    !> of conjugate gradient/steepest descent the geometries are written anyway
    logical :: tWriteRestart = .false.

    !> Minimal number of SCC iterations
    integer :: minSCCIter

    !> if scc/geometry driver should be stopped
    logical :: tStopSCC, tStopDriver

    !> Whether scc restart info should be written in current iteration
    logical :: tWriteSccRestart

    !> Whether charges should be written
    logical :: tWriteCharges

    !> net charge on each atom
    real(dp), allocatable :: dqAtom(:)

    !> density matrix
    real(dp), allocatable :: rhoSqrReal(:,:,:)

    !> Natural orbitals for excited state density matrix, if requested
    real(dp), allocatable, target :: occNatural(:)

    !> locality measure for the wavefunction
    real(dp) :: localisation

    ! set up output files
    call initOutputFiles(tWriteAutotest, tWriteResultsTag, tWriteBandDat, tDerivs,&
        & tWriteDetailedOut, tMd, tGeoOpt, fdAutotest, fdResultsTag, fdBand, fdEigvec, fdHessian,&
        & fdUser, fdMd, geoOutFile)

    ! set up larger arrays
    call initArrays(tForces, tExtChrg, tLinResp, tLinRespZVect, tMd, tMulliken, tSpinOrbit, tImHam,&
        & tStoreEigvecs, tWriteRealHS, tWriteHS, t2Component, tRealHS, tPrintExcitedEigvecs,&
        & tDipole, orb, nAtom, nMovedAtom, nKPoint, nSpin, nExtChrg, indMovedAtom, mass, rhoPrim,&
        & h0, iRhoPrim, excitedDerivs, ERhoPrim, totalDeriv, chrgForces, energy, potential, TS, E0,&
        & Eband, eigen, filling, coord0Fold, newCoords, orbitalL, HSqrCplx, SSqrCplx, HSqrReal,&
        & SSqrReal, rhoSqrReal, dqAtom, chargePerShell, occNatural, velocities, movedVelo,&
        & movedAccel, movedMass, dipoleMoment)

    if (tShowFoldedCoord) then
      pCoord0Out => coord0Fold
    else
      pCoord0Out => coord0
    end if

    call initGeoOptParameters(tCoordOpt, nGeoSteps, tGeomEnd, tCoordStep, tStopDriver, iGeoStep,&
        & iLatGeoStep)

    minSccIter = getMinSccIters(tScc, tDftbU, nSpin)

    if (tXlbomd) then
      call xlbomdIntegrator%setDefaultSCCParameters(minSCCiter, maxSccIter, sccTol)
    end if

    ! If the geometry is periodic, need to update lattice information in geometry loop
    tLatticeChanged = tPeriodic

    ! As first geometry iteration, require updates for coordinates in dependent routines
    tCoordsChanged = .true.

    ! Belongs in initprogram where initial geometry is set up
    if (tSocket) then
      call receiveGeometryFromSocket(socket, tPeriodic, coord0, latVec, tCoordsChanged,&
          & tLatticeChanged, tStopDriver)
    end if

    ! Main geometry loop
    lpGeomOpt: do iGeoStep = 0, nGeoSteps

      call printGeoStepInfo(tCoordOpt, tLatOpt, iLatGeoStep, iGeoStep)

      tWriteRestart = needsRestartWriting(tGeoOpt, tMd, iGeoStep, nGeoSteps, restartFreq)
      if (tMD .and. tWriteRestart) then
        call writeMdOut1(fdMd, mdOut, iGeoStep, pMDIntegrator)
      end if

      if (tLatticeChanged) then
        call handleLatticeChange(latVec, tScc, tStress, pressure, mCutoff, dispersion,&
            & recVec, invLatVec, cellVol, recCellVol, derivCellVol, cellVec, rCellVec)
      end if

      if (tCoordsChanged) then
        call handleCoordinateChange(coord0, latVec, invLatVec, species0, mCutoff, skRepCutoff, orb,&
            & tPeriodic, tScc, tDispersion, dispersion, thirdOrd, img2CentCell, iCellVec,&
            & neighborList, nAllAtom, coord0Fold, coord, species, rCellVec, nAllOrb, nNeighbor,&
            & ham, over, H0, rhoPrim, iRhoPrim, iHam, ERhoPrim, iPair)
      end if

      if (tSCC) then
        call reset(pChrgMixer, nMixElements)
      end if

      call buildH0(H0, skHamCont, atomEigVal, coord, nNeighbor, neighborList%iNeighbor, species,&
          & iPair, orb)
      call buildS(over, skOverCont, coord, nNeighbor, neighborList%iNeighbor, species, iPair, orb)

      if (tSetFillingTemp) then
        call getTemperature(temperatureProfile, tempElec)
      end if

      call calcRepulsiveEnergy(coord, species, img2CentCell, nNeighbor, neighborList, pRepCont,&
          & energy%atomRep, energy%ERep)

      if (tDispersion) then
        call calcDispersionEnergy(dispersion, energy%atomDisp, energy%Edisp)
      end if

      call resetExternalPotentials(potential)
      if (tEField) then
        call setUpExternalElectricField(tTDEField, tPeriodic, EFieldStrength, EFieldVector,&
            & EFieldOmega, EFieldPhase, neighborList, nNeighbor, iCellVec, img2CentCell, cellVec,&
            & deltaT, iGeoStep, coord0Fold, coord, EField, potential%extAtom(:,1), absEField)
      end if

      call mergeExternalPotentials(orb, species, potential)

      call initSccLoop(tScc, xlbomdIntegrator, minSccIter, maxSccIter, sccTol, tConverged)

      lpSCC: do iSccIter = 1, maxSccIter

        call resetInternalPotentials(tDualSpinOrbit, xi, orb, species, potential)

        if (tScc) then
          call getChargePerShell(qInput, orb, species, chargePerShell)
          call addChargePotentials(qInput, q0, chargePerShell, orb, species, species0,&
              & neighborList, img2CentCell, spinW, thirdOrd, potential)
          call addBlockChargePotentials(qBlockIn, qiBlockIn, tDftbU, tImHam, species, orb,&
              & nDftbUFunc, UJ, nUJ, iUJ, niUJ, potential)
        end if
        potential%intBlock = potential%intBlock + potential%extBlock

        call getSccHamiltonian(H0, over, nNeighbor, neighborList, species, orb, iPair,&
            & img2CentCell, potential, ham, iHam)

        if (tWriteRealHS .or. tWriteHS) then
          call writeHSAndStop(tWriteHS, tWriteRealHS, tRealHS, over, neighborList, nNeighbor,&
              & iAtomStart, iPair, img2CentCell, kPoint, iCellVec, cellVec, ham, iHam)
        end if

        call getDensity(ham, over, neighborList, nNeighbor, iAtomStart, iPair, img2CentCell,&
            & iCellVec, cellVec, kPoint, kWeight, orb, species, solver, tRealHS, tSpinSharedEf,&
            & tSpinOrbit, tDualSpinOrbit, tFillKSep, tFixEf, tMulliken, iDistribFn, tempElec, nEl,&
            & Ef, energy, eigen, filling, rhoPrim, Eband, TS, E0, iHam, xi, orbitalL, HSqrReal,&
            & SSqrReal, iRhoPrim, HSqrCplx, SSqrCplx, rhoSqrReal, storeEigvecsReal,&
            & storeEigvecsCplx)

        if (tWriteBandDat) then
          call writeBandOut(fdBand, bandOut, eigen, filling, kWeight)
        end if

        if (tMulliken) then
          call getMullikenPopulation(rhoPrim, over, orb, neighborList, nNeighbor, img2CentCell,&
              & iPair, qOutput, iRhoPrim=iRhoPrim, qBlock=qBlockOut, qiBlock=qiBlockOut)
        end if

        ! For non-dual spin-orbit orbitalL is determined during getDensity() call above
        if (tDualSpinOrbit) then
          call getL(orbitalL, qiBlockOut, orb, species)
        end if

        ! Note: if XLBOMD is active, potential created with input charges is needed later,
        ! therefore it should not be overwritten here.
        if (tSCC .and. .not. tXlbomd) then
          call resetInternalPotentials(tDualSpinOrbit, xi, orb, species, potential)
          call getChargePerShell(qOutput, orb, species, chargePerShell)
          call addChargePotentials(qOutput, q0, chargePerShell, orb, species, species0,&
              & neighborList, img2CentCell, spinW, thirdOrd, potential)
          call addBlockChargePotentials(qBlockOut, qiBlockOut, tDftbU, tImHam, species, orb,&
              & nDftbUFunc, UJ, nUJ, iUJ, niUJ, potential)
          potential%intBlock = potential%intBlock + potential%extBlock
        end if

        call getEnergies(qOutput, q0, chargePerShell, species, tEField, tScc, tXlbomd, tDftbU,&
            & tDualSpinOrbit, rhoPrim, H0, orb, neighborList, nNeighbor, img2CentCell, iPair,&
            & cellVol, pressure, TS, potential, energy, thirdOrd, qBlockOut, qiBlockOut,&
            & nDftbUFunc, UJ, nUJ, iUJ, niUJ, xi)

        tStopScc = hasStopFile(fStopScc)

        if (tScc) then
          call getNextInputCharges(pChrgMixer, qOutput, qOutRed, orb, nIneqOrb, iEqOrbitals,&
              & iGeoStep, iSccIter, minSccIter, maxSccIter, sccTol, tStopScc, tDftbU, tReadChrg,&
              & qInput, qInpRed, sccErrorQ, tConverged, qBlockOut, iEqBlockDftbU, qBlockIn,&
              & qiBlockOut, iEqBlockDftbULS, species0, nUJ, iUJ, niUJ, qiBlockIn)
          call getSccInfo(iSccIter, energy%Eelec, Eold, diffElec)
          call printSccInfo(tDftbU, iSccIter, energy%Eelec, diffElec, sccErrorQ)
          tWriteSccRestart = needsSccRestartWriting(restartFreq, iGeoStep, iSccIter, minSccIter,&
              & maxSccIter, tMd, tGeoOpt, tDerivs, tConverged, tReadChrg, tStopScc)
          if (tWriteSccRestart) then
            call writeRestart(fChargeIn, orb, qInput, qBlockIn, qiBlockIn)
          end if
        end if

        if (tWriteDetailedOut) then
          call writeDetailedOut1(fdUser, userOut, tAppendDetailedOut, iDistribFn, nGeoSteps,&
              & iGeoStep, tMD, tDerivs, tCoordOpt, tLatOpt, iLatGeoStep, iSccIter, energy,&
              & diffElec, sccErrorQ, indMovedAtom, pCoord0Out, q0, qInput, qOutput, eigen, filling,&
              & orb, species, tDFTBU, tImHam, tPrintMulliken, orbitalL, qBlockOut, Ef, Eband, TS,&
              & E0, pressure, cellVol, tAtomicEnergy, tDispersion, tEField, tPeriodic, nSpin,&
              & tSpinOrbit, tScc)
        end if

        if (tConverged .or. tStopScc) then
          exit lpSCC
        end if

      end do lpSCC

      if (tLinResp) then
        call ensureLinRespConditions(t3rd, tRealHS, tPeriodic, tForces)
        call calculateLinRespExcitations(lresp, qOutput, q0, over, HSqrReal, eigen(:,1,:),&
            & filling(:,1,:), coord0, species, species0, speciesName, orb, skHamCont, skOverCont,&
            & fdAutotest, fdEigvec, runId, neighborList, nNeighbor, iAtomStart, iPair,&
            & img2CentCell, tWriteAutotest, tForces, tLinRespZVect, tPrintExcitedEigvecs,&
            & nonSccDeriv, energy, SSqrReal, rhoSqrReal, excitedDerivs, occNatural)
      end if

      if (tXlbomd) then
        call getXlbomdCharges(xlbomdIntegrator, qOutRed, pChrgMixer, orb, nIneqOrb, iEqOrbitals,&
            & qInput, qInpRed, iEqBlockDftbU, qBlockIn, species0, nUJ, iUJ, niUJ, iEqBlockDftbuLs,&
            & qiBlockIn)
      end if

      if (tDipole) then
        call getDipoleMoment(qOutput, q0, coord, dipoleMoment)
      #:call DEBUG_CODE
        call checkDipoleViaHellmannFeynman(size(h0), rhoPrim, q0, coord0, over, orb, neighborList,&
            & nNeighbor, species, iPair, img2CentCell)
      #:endcall DEBUG_CODE
      end if

      if (tPrintEigVecs) then
        call writeEigenvectors(nSpin, fdEigvec, runId, nAtom, neighborList, nNeighbor, cellVec,&
            & iCellVEc, iAtomStart, iPair, img2CentCell, species, speciesName, orb, kPoint, over,&
            & HSqrReal, SSqrReal, HSqrCplx, SSqrCplx, storeEigvecsReal, storeEigvecsCplx)
      end if

      if (tProjEigenvecs) then
        call writeProjectedEigenvectors(regionLabels, fdProjEig, eigen, nSpin, neighborList,&
            & nNeighbor, cellVec, iCellVec, iAtomStart, iPair, img2CentCell, orb, over, kPoint,&
            & kWeight, iOrbRegion, HSqrReal, SSqrReal, HSqrCplx, SSqrCplx, storeEigvecsReal,&
            & storeEigvecsCplx)
      end if

      ! MD geometry files are written only later, once velocities for the current geometry are known
      if (tGeoOpt .and. tWriteRestart) then
        call writeCurrentGeometry(geoOutFile, pCoord0Out, tLatOpt, tMd, tAppendGeo, tFracCoord,&
            & tPeriodic, tPrintMulliken, species0, speciesName, latVec, iGeoStep, iLatGeoStep,&
            & nSpin, qOutput, velocities)
      end if

      call printEnergies(energy)

      if (tForces) then
        call getEnergyWeightedDensityMtx(forceType, filling, eigen, kPoint, kWeight, neighborList,&
            & nNeighbor, orb, iAtomStart, iPair, img2CentCell, iCellVEc, cellVec, tRealHS, ham,&
            & over, solver, ERhoPrim, HSqrReal, SSqrReal, HSqrCplx, SSqrCplx, storeEigvecsReal,&
            & storeEigvecsCplx)
        call getGradients(tScc, tEField, tXlbomd, nonSccDeriv, Efield, rhoPrim, ERhoPrim, qOutput,&
            & q0, skHamCont, skOverCont, pRepCont, neighborList, nNeighbor, species, img2CentCell,&
            & iPair, orb, potential, coord, dispersion, totalDeriv, iRhoPrim, thirdOrd, chrgForces)
        if (tLinResp) then
          totalDeriv(:,:) = totalDeriv(:,:) + excitedDerivs(:,:)
        end if

        if (tStress) then
          call getStress(tScc, tEField, nonSccDeriv, EField, rhoPrim, ERhoPrim, qOutput, q0,&
              & skHamCont, skOverCont, pRepCont, neighborList, nNeighbor, species, img2CentCell,&
              & iPair, orb, potential, coord, latVec, invLatVec, cellVol, coord0, dispersion,&
              & totalStress, totalLatDeriv, cellPressure, iRhoPrim)
          call printVolume(cellVol)
          ! MD case includes the atomic kinetic energy contribution, so print that later
          if (.not. tMD) then
            call printPressureAndFreeEnergy(pressure, cellPressure, energy%EGibbs)
          end if
        end if
      end if

      if (tWriteDetailedOut) then
        call writeDetailedOut2(fdUser, tScc, tConverged, tXlbomd, tLinResp, tGeoOpt, tMD,&
            & tPrintForces, tStress, tPeriodic, energy, totalStress, totalLatDeriv, totalDeriv, &
            & chrgForces, indMovedAtom, cellVol, cellPressure, geoOutFile)
      end if

      if (tScc .and. .not. tXlbomd .and. .not. tConverged) then
        if (tConvrgForces) then
          call error("SCC is NOT converged, maximal SCC iterations exceeded")
        else
          call warning("SCC is NOT converged, maximal SCC iterations exceeded")
        end if
      end if

      if (tForces) then
        if (allocated(conAtom)) then
          call constrainForces(conAtom, conVec, totalDeriv)
        end if

        if (tCoordOpt) then
          call printMaxForce(maxval(abs(totalDeriv(:, indMovedAtom))))
        end if

        if (tLatOpt) then
          ! Only include the derivCellVol contribution if not MD, as the barostat would otherwise
          ! take care of this, hence add it here rather than to totalLatDeriv itself
          call constrainLatticeDerivs(totalLatDeriv + derivCellVol, normOrigLatVec,&
              & tLatOptFixAng, tLatOptFixLen, tLatOptIsotropic, constrLatDerivs)
          call printMaxLatticeForce(maxval(abs(constrLatDerivs)))
        end if

        if (tSocket) then
          ! stress was computed above in the force evaluation block or is 0 if aperiodic
          call socket%send(energy%ETotal - sum(TS), -totalDeriv, totalStress * cellVol)
        end if

        ! If geometry minimizer finished and the last calculated geometry is the minimal one (not
        ! necessarily the case, depends on the optimizer!) we are finished.  Otherwise we have to
        ! recalculate everything at the converged geometry.
        
        if (tGeomEnd) then
          exit lpGeomOpt
        end if

        tWriteCharges = tWriteRestart .and. tMulliken .and. tSCC .and. .not. tDerivs&
            & .and. maxSccIter > 1
        if (tWriteCharges) then
          call writeCharges(qInput, fChargeIn, orb, qBlockIn, qiBlockIn)
        end if

        ! initially assume coordinates are not being updated
        tCoordsChanged = .false.
        tLatticeChanged = .false.

        if (tDerivs) then
          call getNextDerivStep(derivDriver, totalDeriv, indMovedAtom, coord0, tGeomEnd)
          if (tGeomEnd) then
            exit lpGeomOpt
          end if
          tCoordsChanged = .true.
        else if (tGeoOpt) then
          tCoordsChanged = .true.
          if (tCoordStep) then
            call getNextCoordinateOptStep(pGeoCoordOpt, energy%EMermin, totalDeriv, indMovedAtom,&
                & coord0, diffGeo, tCoordEnd)
            if (.not. tLatOpt) then
              tGeomEnd = tCoordEnd
            end if
            if (.not. tGeomEnd .and. tCoordEnd .and. diffGeo < tolSameDist) then
              tCoordStep = .false.
            end if
          else
            call getNextLatticeOptStep(pGeoLatOpt, energy%EGibbs, constrLatDerivs, origLatVec,&
                & tLatOptFixAng, tLatOptFixLen, tLatOptIsotropic, indMovedAtom, latVec, coord0,&
                & diffGeo, tGeomEnd)
            iLatGeoStep = iLatGeoStep + 1
            tLatticeChanged = .true.
            if (.not. tGeomEnd .and. tCoordOpt) then
              tCoordStep = .true.
              call reset(pGeoCoordOpt, reshape(coord0(:, indMovedAtom), [nMovedCoord]))
            end if
          end if
          if (tGeomEnd .and. diffGeo < tolSameDist) then
            exit lpGeomOpt
          end if
        else if (tMD) then
          ! New MD coordinates saved in a temporary variable, as writeCurrentGeometry() below
          ! needs the old ones to write out consistent geometries and velocities.
          newCoords(:,:) = coord0
          call getNextMdStep(pMdIntegrator, pMdFrame, temperatureProfile, totalDeriv, movedMass,&
              & mass, cellVol, invLatVec, species0, indMovedAtom, tStress, tBarostat, energy,&
              & newCoords, latVec, cellPressure, totalStress, totalLatDeriv, velocities, kT)
          tCoordsChanged = .true.
          tLatticeChanged = tBarostat
          call printMdInfo(tSetFillingTemp, tEField, tPeriodic, tempElec, absEField, kT,&
              & cellPressure, pressure, energy)
          if (tWriteRestart) then
            if (tPeriodic) then
              cellVol = abs(determinant33(latVec))
              energy%EGibbs = energy%EMermin + pressure * cellVol
            end if
            call writeMdOut2(fdMd, tStress, tBarostat, tLinResp, tEField, tFixEf, tPrintMulliken,&
                & energy, latVec, cellVol, cellPressure, pressure, kT, absEField, qOutput, q0,&
                & dipoleMoment)
            call writeCurrentGeometry(geoOutFile, pCoord0Out, .false., .true., .true., tFracCoord,&
                & tPeriodic, tPrintMulliken, species0, speciesName, latVec, iGeoStep, iLatGeoStep,&
                & nSpin, qOutput, velocities)
          end if
          coord0(:,:) = newCoords
          if (tWriteDetailedOut) then
            call writeDetailedOut3(fdUser, tPrintForces, tSetFillingTemp, tPeriodic, tStress,&
                & totalStress, totalLatDeriv, energy, tempElec, pressure, cellPressure, kT)
          end if
        else if (tSocket .and. iGeoStep < nGeoSteps) then
          ! Only receive geometry from socket, if there are still geometry iterations left
          call receiveGeometryFromSocket(socket, tPeriodic, coord0, latVec, tCoordsChanged,&
              & tLatticeChanged, tStopDriver)
        end if
      end if

      if (tWriteDetailedOut) then
        call writeDetailedOut4(fdUser, tMD, energy, kT)
      end if

      tStopDriver = tStopScc .or. tStopDriver .or. hasStopFile(fStopDriver)
      if (tStopDriver) then
        exit lpGeomOpt
      end if

    end do lpGeomOpt

    if (tSocket) then
      call socket%shutdown()
    end if

    tGeomEnd = tMD .or. tGeomEnd .or. tDerivs

    if (tWriteDetailedOut) then
      call writeDetailedOut5(fdUser, tGeoOpt, tGeomEnd, tMd, tDerivs, tEField, absEField,&
          & dipoleMoment)
    end if

    call writeFinalDriverStatus(tGeoOpt, tGeomEnd, tMd, tDerivs)

    if (tMD) then
      call writeMdOut3(fdMd, mdOut)
    end if

    if (tDerivs) then
      call getHessianMatrix(derivDriver, pDynMatrix)
      call writeHessianOut(fdHessian, hessianOut, pDynMatrix)
    else
      nullify(pDynMatrix)
    end if


    if (allocated(pipekMezey)) then
      ! NOTE: the canonical DFTB ground state orbitals are over-written after this point
      if (tStoreEigvecs) then
        call error("Pipek-Mezey localisation not implemented for stored eigenvectors")
      end if
      if (nSpin > 2) then
        call error("Pipek-Mezey localisation not implemented for non-colinear DFTB")
      end if
      call calcPipekMezeyLocalisation(pipekMezey, nEl, filling, over, kPoint, kWeight,&
          & neighborList, nNeighbor, iAtomStart, iPair, img2CentCell, iCellVec, cellVec, fdEigvec,&
          & runId, orb, species, speciesName, localisation, HSqrReal, SSqrReal, HsqrCplx, SSqrCplx)
    end if

    if (tWriteAutotest) then
      if (tPeriodic) then
        cellVol = abs(determinant33(latVec))
        energy%EGibbs = energy%EMermin + pressure * cellVol
      end if
      call writeAutotestTag(fdAutotest, autotestTag, tPeriodic, cellVol, tMulliken, qOutput,&
          & totalDeriv, chrgForces, excitedDerivs, tStress, totalStress, pDynMatrix,&
          & energy%EMermin, pressure, energy%EGibbs, coord0, tLocalise, localisation)
    end if
    if (tWriteResultsTag) then
      call writeResultsTag(fdResultsTag, resultsTag, totalDeriv, chrgForces, tStress, totalStress,&
          & pDynMatrix, tPeriodic, cellVol)
    end if
    if (tWriteDetailedXML) then
      call writeDetailedXml(runId, speciesName, species0, pCoord0Out, tPeriodic, latVec, tRealHS,&
          & nKPoint, nSpin, size(eigen, dim=1), nOrb, kPoint, kWeight, filling, occNatural)
    end if

    call destructProgramVariables()

  end subroutine runDftbPlus


  !> Initialises (clears) output files.
  subroutine initOutputFiles(tWriteAutotest, tWriteResultsTag, tWriteBandDat, tDerivs,&
      & tWriteDetailedOut, tMd, tGeoOpt, fdAutotest, fdResultsTag, fdBand, fdEigvec, fdHessian,&
      & fdUser, fdMd, geoOutFile)

    !> Should tagged regression test data be printed
    logical, intent(in) :: tWriteAutotest

    !> Write tagged output for machine readable results
    logical, intent(in) :: tWriteResultsTag

    !> Band structure and fillings
    logical, intent(in) :: tWriteBandDat

    !> Finite different second derivatives
    logical, intent(in) :: tDerivs

    !> Write detail on the calculation to file
    logical, intent(in) :: tWriteDetailedOut

    !> Is this a molecular dynamics calculation
    logical, intent(in) :: tMd

    !> Are atomic coodinates being optimised
    logical, intent(in) :: tGeoOpt

    !> File unit for autotest data
    integer, intent(out) :: fdAutotest

    !> File unit for tagged results data
    integer, intent(out) :: fdResultsTag

    !> File unit for band structure
    integer, intent(out) :: fdBand

    !> File unit for eigenvectors
    integer, intent(out) :: fdEigvec

    !> File unit for second derivatives information
    integer, intent(out) :: fdHessian

    !> File unit for detailed.out
    integer, intent(out) :: fdUser

    !> File unit for information during molecular dynamics
    integer, intent(out) :: fdMd

    !> Filename for geometry output
    character(*), intent(in) :: geoOutFile

    call initTaggedWriter()
    if (tWriteAutotest) then
      call initOutputFile(autotestTag, fdAutotest)
    end if
    if (tWriteResultsTag) then
      call initOutputFile(resultsTag, fdResultsTag)
    end if
    if (tWriteBandDat) then
      call initOutputFile(bandOut, fdBand)
    end if
    fdEigvec = getFileId()
    if (tDerivs) then
      call initOutputFile(hessianOut, fdHessian)
    end if
    if (tWriteDetailedOut) then
      call initOutputFile(userOut, fdUser)
    end if
    if (tMD) then
      call initOutputFile(mdOut, fdMD)
    end if
    if (tGeoOpt .or. tMD) then
      call clearFile(trim(geoOutFile) // ".gen")
      call clearFile(trim(geoOutFile) // ".xyz")
    end if

  end subroutine initOutputFiles


  !> Allocates most of the large arrays needed during the DFTB run.
  subroutine initArrays(tForces, tExtChrg, tLinResp, tLinRespZVect, tMd, tMulliken, tSpinOrbit,&
      & tImHam, tStoreEigvecs, tWriteRealHS, tWriteHS, t2Component, tRealHS, tPrintExcitedEigvecs,&
      & tDipole, orb, nAtom, nMovedAtom, nKPoint, nSpin, nExtChrg, indMovedAtom, mass, rhoPrim, h0,&
      & iRhoPrim, excitedDerivs, ERhoPrim, totalDeriv, chrgForces, energy, potential, TS, E0,&
      & Eband, eigen, filling, coord0Fold, newCoords, orbitalL, HSqrCplx, SSqrCplx, HSqrReal,&
      & SSqrReal, rhoSqrReal, dqAtom, chargePerShell, occNatural, velocities, movedVelo,&
      & movedAccel, movedMass, dipoleMoment)

    !> Are forces required
    logical, intent(in) :: tForces

    !> Are the external charges
    logical, intent(in) :: tExtChrg

    !> Are excitation energies being calculated
    logical, intent(in) :: tLinResp

    !> Are excited state properties being calculated
    logical, intent(in) :: tLinRespZVect

    !> Is this a molecular dynamics calculation
    logical, intent(in) :: tMd

    !> Is Mulliken analysis being performed
    logical, intent(in) :: tMulliken

    !> Are there spin orbit interactions
    logical, intent(in) :: tSpinOrbit

    !> Are there imaginary parts to the hamiltonian
    logical, intent(in) :: tImHam

    !> Should eigenvectors be stored on disc
    logical, intent(in) :: tStoreEigvecs

    !> Should the sparse hamiltonian and overlap be writen to disc?
    logical, intent(in) :: tWriteRealHS

    !> Should the hamiltonian and overlap be written out as dense matrices
    logical, intent(in) :: tWriteHS

    !> Is the hamiltonian for a two component (Pauli) system
    logical, intent(in) :: t2Component

    !> Is the hamiltonian real?
    logical, intent(in) :: tRealHS

    !> Print the excited state eigenvectors
    logical, intent(in) :: tPrintExcitedEigvecs

    !> Print the dipole moment
    logical, intent(in) :: tDipole

    !> data structure with atomic orbital information
    type(TOrbitals), intent(in) :: orb

    !> Number of atoms
    integer, intent(in) :: nAtom

    !> Number of atoms moved about during the calculation
    integer, intent(in) :: nMovedAtom

    !> Number of k-points
    integer, intent(in) :: nKPoint

    !> Number of spin channels
    integer, intent(in) :: nSpin

    !> Number of external charges
    integer, intent(in) :: nExtChrg

    !> Indices for any moving atoms
    integer, intent(in) :: indMovedAtom(:)

    !> Masses of atoms
    real(dp), intent(in) :: mass(:)

    !> Sparse storage density matrix
    real(dp), intent(out), allocatable :: rhoPrim(:,:)

    !> Non-scc part of the hamiltonian
    real(dp), intent(out), allocatable :: h0(:)

    !> Imaginary part of the sparse density matrix
    real(dp), intent(out), allocatable :: iRhoPrim(:,:)

    !> Excitation energy derivatives with respect to atomic coordinates
    real(dp), intent(out), allocatable :: excitedDerivs(:,:)

    !> Energy weighted density matrix
    real(dp), intent(out), allocatable :: ERhoPrim(:)

    !> Derivatives of total energy with respect to atomic coordinates
    real(dp), intent(out), allocatable :: totalDeriv(:,:)

    !> Forces on (any) external charges
    real(dp), intent(out), allocatable :: chrgForces(:,:)

    !> Energy terms
    type(TEnergies), intent(out) :: energy

    !> Potentials acting on the system
    type(TPotentials), intent(out) :: potential

    !> Electron entropy contribution at T
    real(dp), intent(out), allocatable :: TS(:)

    !> zero temperature extrapolated electronic energy
    real(dp), intent(out), allocatable :: E0(:)

    !> band  energy
    real(dp), intent(out), allocatable :: Eband(:)

    !> single particle energies (band structure)
    real(dp), intent(out), allocatable :: eigen(:,:,:)

    !> Occupations of single particle states
    real(dp), intent(out), allocatable :: filling(:,:,:)

    !> Coordinates in central cell
    real(dp), intent(out), allocatable :: coord0Fold(:,:)

    !> Updated coordinates
    real(dp), intent(out), allocatable :: newCoords(:,:)

    !> Orbital angular momentum
    real(dp), intent(out), allocatable :: orbitalL(:,:,:)

    !> Complex dense hamiltonian
    complex(dp), intent(out), allocatable :: HSqrCplx(:,:,:,:)

    !> overlap matrix dense storage
    complex(dp), intent(out), allocatable :: SSqrCplx(:,:)

    !> real dense hamiltonian
    real(dp), intent(out), allocatable :: HSqrReal(:,:,:)

    !> overlap matrix dense storage
    real(dp), intent(out), allocatable :: SSqrReal(:,:)

    !> density matrix dense storage
    real(dp), intent(out), allocatable :: rhoSqrReal(:,:,:)

    !> Atomic net charges
    real(dp), intent(out), allocatable :: dqAtom(:)

    !> Number of electron in each atomic shell
    real(dp), intent(out), allocatable :: chargePerShell(:,:,:)

    !> Occupations for natural orbitals
    real(dp), intent(out), allocatable :: occNatural(:)

    !> Atomic velocities
    real(dp), intent(out), allocatable :: velocities(:,:)

    !> Array for moving atom velocities
    real(dp), intent(out), allocatable :: movedVelo(:,:)

    !> moving atom accelerations
    real(dp), intent(out), allocatable :: movedAccel(:,:)

    !> moving atoms masses
    real(dp), intent(out), allocatable :: movedMass(:,:)

    !> system dipole moment
    real(dp), intent(out), allocatable :: dipoleMoment(:)

    integer :: nSpinHams, sqrHamSize
    integer :: nSpinStored, nKPointStored

    allocate(rhoPrim(0, nSpin))
    allocate(h0(0))
    if (tImHam) then
      allocate(iRhoPrim(0, nSpin))
    end if

    if (tForces) then
      allocate(ERhoPrim(0))
      allocate(totalDeriv(3, nAtom))
      if (tExtChrg) then
        allocate(chrgForces(3, nExtChrg))
      end if
    end if
    if (tLinRespZVect) then
      allocate(excitedDerivs(3, nAtom))
    end if

    call init(energy, nAtom)
    call init(potential, orb, nAtom, nSpin)

    ! Nr. of independent spin Hamiltonians
    select case (nSpin)
    case (1)
      nSpinHams = 1
      sqrHamSize = orb%nOrb
    case (2)
      nSpinHams = 2
      sqrHamSize = orb%nOrb
    case (4)
      nSpinHams = 1
      sqrHamSize = 2 * orb%nOrb
    end select

    allocate(TS(nSpinHams))
    allocate(E0(nSpinHams))
    allocate(Eband(nSpinHams))
    allocate(eigen(sqrHamSize, nKPoint, nSpinHams))
    allocate(filling(sqrHamSize, nKpoint, nSpinHams))

    allocate(coord0Fold(3, nAtom))

    if (tMD) then
      allocate(newCoords(3, nAtom))
    end if

    if ((tMulliken .and. tSpinOrbit) .or. tImHam) then
      allocate(orbitalL(3, orb%mShell, nAtom))
    end if

    if (tStoreEigvecs) then
      nSpinStored = 1
      nKPointStored = 1
    else
      nSpinStored = nSpin
      nKPointStored = nKPoint
    end if

    ! If only H/S should be printed, no allocation for square HS is needed
    if (.not. (tWriteRealHS .or. tWriteHS)) then
      if (t2Component) then
        allocate(HSqrCplx(sqrHamSize, sqrHamSize, nKPointStored, 1))
        allocate(SSqrCplx(sqrHamSize, sqrHamSize))
      elseif (tRealHS) then
        allocate(HSqrReal(sqrHamSize, sqrHamSize, nSpinStored))
        allocate(SSqrReal(sqrHamSize, sqrHamSize))
      else
        allocate(HSqrCplx(sqrHamSize, sqrHamSize, nKPointStored, nSpinStored))
        allocate(SSqrCplx(sqrHamSize, sqrHamSize))
      end if
    end if

    if (tLinResp) then
      allocate(dqAtom(nAtom))
      if (tLinRespZVect) then
        allocate(rhoSqrReal(sqrHamSize, sqrHamSize, nSpin))
      end if
    end if
    allocate(chargePerShell(orb%mShell, nAtom, nSpin))

    if (tLinResp .and. tPrintExcitedEigVecs) then
      allocate(occNatural(orb%nOrb))
    end if

    if (tMD) then
      allocate(velocities(3, nAtom))
      allocate(movedVelo(3, nMovedAtom))
      allocate(movedAccel(3, nMovedAtom))
      allocate(movedMass(3, nMovedAtom))
      movedMass(:,:) = spread(mass(indMovedAtom),1,3)
    end if

    if (tDipole) then
      allocate(dipoleMoment(3))
    end if

  end subroutine initArrays


  !> Initialises some parameters before geometry loop starts.
  subroutine initGeoOptParameters(tCoordOpt, nGeoSteps, tGeomEnd, tCoordStep, tStopDriver,&
      & iGeoStep, iLatGeoStep)

    !> Are atomic coordinates changing
    logical, intent(in) :: tCoordOpt

    !> Number of geometry steps
    integer, intent(in) :: nGeoSteps

    !> Have the geometry changes terminated
    logical, intent(out) :: tGeomEnd

    !> Are the atomic coordinates changing
    logical, intent(out) :: tCoordStep

    !> Should the geometry driver stop
    logical, intent(out) :: tStopDriver

    !> Step of the geometry driver
    integer, intent(out) :: iGeoStep

    !> Number of steps changing the lattice vectors
    integer, intent(out) :: iLatGeoStep

    tGeomEnd = (nGeoSteps == 0)

    tCoordStep = .false.
    if (tCoordOpt) then
      tCoordStep = .true.
    end if
    tStopDriver = .false.

    iGeoStep = 0
    iLatGeoStep = 0

  end subroutine initGeoOptParameters


  !> Receives the geometry from socket communication.
  subroutine receiveGeometryFromSocket(socket, tPeriodic, coord0, latVecs, tCoordsChanged,&
      & tLatticeChanged, tStopDriver)

    !> Socket communication object
    type(IpiSocketComm), allocatable, intent(inout) :: socket

    !> Is the system periodic
    logical, intent(in) :: tPeriodic

    !> Coordinates for atoms
    real(dp), intent(inout) :: coord0(:,:)

    !> Lattice vectors for the unit cell (not referenced if not periodic)
    real(dp), intent(inout) :: latVecs(:,:)

    !> Have the atomic coordinates changed
    logical, intent(out) :: tCoordsChanged

    !> Have the lattice vectors changed
    logical, intent(out) :: tLatticeChanged

    !> Stop the geometry driver if true
    logical, intent(out) :: tStopDriver

    real(dp) :: tmpLatVecs(3, 3)

    call socket%receive(coord0, tmpLatVecs, tStopDriver)
    tCoordsChanged = .true.
    if (tPeriodic .and. .not. tStopDriver) then
      latVecs(:,:) = tmpLatVecs
    end if
    tLatticeChanged = tPeriodic

  end subroutine receiveGeometryFromSocket


  !> Initialises SCC related parameters before geometry loop starts
  function getMinSccIters(tScc, tDftbU, nSpin) result(minSccIter)

    !> Is this a self consistent calculation
    logical, intent(in) :: tScc

    !> Are there orbital potentials present
    logical, intent(in) :: tDftbU

    !> Number of spin channels
    integer, intent(in) :: nSpin

    !> Minimum possible number of self consistent iterations
    integer :: minSccIter

    if (tScc) then
      if (tDftbU) then
        minSccIter = 2
      else
        if (nSpin == 1) then
          minSccIter = 1
        else
          minSccIter = 2
        end if
      end if
    else
      minSccIter = 1
    end if

  end function getMinSccIters


  !> Does the operations that are necessary after a lattice vector update
  subroutine handleLatticeChange(latVecs, tScc, tStress, pressure, mCutoff, dispersion,&
      & recVecs, recVecs2p, cellVol, recCellVol, derivCellVol, cellVecs, rCellVecs)

    !> lattice vectors
    real(dp), intent(in) :: latVecs(:,:)

    !> self consistent calculation?
    logical, intent(in) :: tScc

    !> evaluate stress
    logical, intent(in) :: tStress

    !> External presure
    real(dp), intent(in) :: pressure

    !> Maximum distance for interactions
    real(dp), intent(inout) :: mCutoff

    !> Dispersion interactions object
    class(DispersionIface), allocatable, intent(inout) :: dispersion

    !> Reciprocal lattice vectors
    real(dp), intent(out) :: recVecs(:,:)

    !> Reciprocal lattice vectors in units of 2 pi
    real(dp), intent(out) :: recVecs2p(:,:)

    !> Unit cell volume
    real(dp), intent(out) :: cellVol

    !> reciprocal lattice unit cell volume
    real(dp), intent(out) :: recCellVol

    !> derivative of pV term
    real(dp), intent(out) :: derivCellVol(:,:)

    !> translation vectors to lattice cells in units of lattice constants
    real(dp), allocatable, intent(out) :: cellVecs(:,:)

    !> Vectors to unit cells in absolute units
    real(dp), allocatable, intent(out) :: rCellVecs(:,:)

    cellVol = abs(determinant33(latVecs))
    recVecs2p(:,:) = latVecs
    call matinv(recVecs2p)
    recVecs2p = transpose(recVecs2p)
    recVecs = 2.0_dp * pi * recVecs2p
    recCellVol = abs(determinant33(recVecs))
    if (tStress) then
      call derivDeterminant33(derivCellVol, latVecs)
      derivCellVol(:,:) = pressure * derivCellVol
    end if
    if (tSCC) then
      call updateLatVecs_SCC(latVecs, recVecs, cellVol)
      mCutoff = max(mCutoff, getSCCCutoff())
    end if
    if (allocated(dispersion)) then
      call dispersion%updateLatVecs(latVecs)
      mCutoff = max(mCutoff, dispersion%getRCutoff())
    end if
    call getCellTranslations(cellVecs, rCellVecs, latVecs, recVecs2p, mCutoff)

  end subroutine handleLatticeChange


  !> Does the operations that are necessary after atomic coordinates change
  subroutine handleCoordinateChange(coord0, latVec, invLatVec, species0, mCutoff, skRepCutoff, &
      & orb, tPeriodic, tScc, tDispersion, dispersion, thirdOrd, img2CentCell, iCellVec,&
      & neighborList, nAllAtom, coord0Fold, coord, species, rCellVec, nAllOrb, nNeighbor, ham,&
      & over, H0, rhoPrim, iRhoPrim, iHam, ERhoPrim, iPair)

    !> Central cell coordinates
    real(dp), intent(in) :: coord0(:,:)

    !> Lattice vectors if periodic
    real(dp), intent(in) :: latVec(:,:)

    !> Inverse of the lattice vectors
    real(dp), intent(in) :: invLatVec(:,:)

    !> chemical species of central cell atoms
    integer, intent(in) :: species0(:)

    !> Longest cutoff distance that neighbour maps are generated
    real(dp), intent(in) :: mCutoff

    !> Cutoff for repulsive interaction from SK data
    real(dp), intent(in) :: skRepCutoff

    !> Atomic orbital information
    type(TOrbitals), intent(in) :: orb

    !> Is the geometry periodic
    logical, intent(in) :: tPeriodic

    !> Charge self consistent mode?
    logical, intent(in) :: tScc

    !> Is dispersion included
    logical, intent(in) :: tDispersion

    !> Dispersion interactions
    class(DispersionIface), allocatable, intent(inout) :: dispersion

    !> Third order SCC interactions
    type(ThirdOrder), allocatable, intent(inout) :: thirdOrd

    !> image atoms to their equivalent in the central cell
    integer, allocatable, intent(inout) :: img2CentCell(:)

    !> Index for which unit cell an atom is in
    integer, allocatable, intent(inout) :: iCellVec(:)

    !> List of neighbouring atoms
    type(TNeighborList), intent(inout) :: neighborList

    !> Total number of atoms including images
    integer, intent(out) :: nAllAtom

    !> Total number of atomic orbitals including image atoms
    integer, intent(out) :: nAllOrb

    !> Central cell atomic coordinates, folded inside the central cell
    real(dp), intent(out) :: coord0Fold(:,:)

    !> Coordinates of all atoms including images
    real(dp), allocatable, intent(inout) :: coord(:,:)

    !> Species of all atoms including images
    integer, allocatable, intent(inout) :: species(:)

    !> Vectors to units cells in absolute units
    real(dp), allocatable, intent(inout) :: rCellVec(:,:)

    !> Number of neighbours of each real atom
    integer, intent(out) :: nNeighbor(:)

    !> Sparse hamiltonian storage
    real(dp), allocatable, intent(inout) :: ham(:,:)

    !> sparse overlap storage
    real(dp), allocatable, intent(inout) :: over(:)

    !> Non-SCC hamitonian storage
    real(dp), allocatable, intent(inout) :: h0(:)

    !> Sparse density matrix storage
    real(dp), allocatable, intent(inout) :: rhoPrim(:,:)

    !> Imaginary part of sparse density matrix storage
    real(dp), allocatable, intent(inout) :: iRhoPrim(:,:)

    !> Imaginary part of sparse hamiltonian storage
    real(dp), allocatable, intent(inout) :: iHam(:,:)

    !> energy weighted density matrix storage
    real(dp), allocatable, intent(inout) :: ERhoPrim(:)

    !> index array for location of atomic blocks in large sparse arrays
    integer, allocatable, intent(inout) :: iPair(:,:)

    !> Total size of orbitals in the sparse data structures, where the decay of the overlap sets the
    !> sparsity pattern
    integer :: sparseSize

    coord0Fold(:,:) = coord0
    if (tPeriodic) then
      call foldCoordToUnitCell(coord0Fold, latVec, invLatVec)
    end if

    call updateNeighborListAndSpecies(coord, species, img2CentCell, iCellVec, &
        &neighborList, nAllAtom, coord0Fold, species0, mCutoff, rCellVec)
    nAllOrb = sum(orb%nOrbSpecies(species(1:nAllAtom)))
    call getNrOfNeighborsForAll(nNeighbor, neighborList, skRepCutoff)
    call getSparseDescriptor(neighborList%iNeighbor, nNeighbor, img2CentCell, orb, iPair,&
        & sparseSize)
    call reallocateSparseArrays(sparseSize, ham, over, H0, rhoPrim, iHam, iRhoPrim, ERhoPrim)

    ! Notify various modules about coordinate changes
    if (tSCC) then
      call updateCoords_SCC(coord, species, neighborList, img2CentCell)
    end if
    if (tDispersion) then
      call dispersion%updateCoords(neighborList, img2CentCell, coord, &
          & species0)
    end if
    if (allocated(thirdOrd)) then
      call thirdOrd%updateCoords(neighborList, species)
    end if

  end subroutine handleCoordinateChange


  !> Decides, whether restart file should be written during the run.
  function needsRestartWriting(tGeoOpt, tMd, iGeoStep, nGeoSteps, restartFreq) result(tWriteRestart)
    
    !> Are geometries being optimised
    logical, intent(in) :: tGeoOpt

    !> Is this a molecular dynamics run
    logical, intent(in) :: tMd

    !> Current geometry step
    integer, intent(in) :: iGeoStep

    !> Number of geometry steps in total
    integer, intent(in) :: nGeoSteps

    !> Frequency of restart in geometry steps
    integer, intent(in) :: restartFreq

    !> Should a restart file be written?
    logical :: tWriteRestart

    if (restartFreq > 0 .and. (tGeoOpt .or. tMD)) then
      tWriteRestart = (iGeoStep == nGeoSteps .or. (mod(iGeoStep, restartFreq) == 0))
    else
      tWriteRestart = .false.
    end if

  end function needsRestartWriting


  !> Ensures that sparse array have enough storage to hold all necessary elements.
  subroutine reallocateSparseArrays(sparseSize, ham, over, H0, rhoPrim, iHam, iRhoPrim, ERhoPrim)

    !> Size of the sparse overlap
    integer, intent(in) :: sparseSize

    !> Sparse storage for hamitonian (sparseSize,nSpin)
    real(dp), allocatable, intent(inout) :: ham(:,:)

    !> Sparse storage for overlap
    real(dp), allocatable, intent(inout) :: over(:)

    !> Sparse storage for non-SCC hamitonian
    real(dp), allocatable, intent(inout) :: H0(:)

    !> Sparse storage for density matrix
    real(dp), allocatable, intent(inout) :: rhoPrim(:,:)

    !> Sparse storage for imaginary hamitonian (not reallocated if not initially allocated)
    real(dp), allocatable, intent(inout) :: iHam(:,:)

    !> Sparse storage for imaginary part of density matrix (not reallocated if not initially
    !> allocated)
    real(dp), allocatable, intent(inout) :: iRhoPrim(:,:)

    !> Sparse storage for energy weighted density matrix (not reallocated if not initially
    !> allocated)
    real(dp), allocatable, intent(inout) :: ERhoPrim(:)

    integer :: nSpin

    #:call ASSERT_CODE
      @:ASSERT(size(over) == size(ham, dim=1))
      @:ASSERT(size(H0) == size(ham, dim=1))
      @:ASSERT(all(shape(rhoPrim) == shape(ham)))
      if (allocated(iRhoPrim)) then
        @:ASSERT(all(shape(iRhoPrim) == shape(ham)))
        @:ASSERT(all(shape(iHam) == shape(ham)))
      end if
      if (allocated(ERhoPrim)) then
        @:ASSERT(all(shape(ERhoPrim) == shape(ham)))
      end if
    #:endcall ASSERT_CODE

    if (size(ham, dim=1) >= sparseSize) then
      ! Sparse matrices are big enough
      return
    end if

    nSpin = size(ham, dim=2)
    deallocate(ham)
    deallocate(over)
    deallocate(H0)
    deallocate(rhoPrim)
    allocate(ham(sparseSize, nSpin))
    allocate(over(sparseSize))
    allocate(H0(sparseSize))
    allocate(rhoPrim(sparseSize, nSpin))
    if (allocated(iRhoPrim)) then
      deallocate(iRhoPrim)
      deallocate(iHam)
      allocate(iRhoPrim(sparseSize, nSpin))
      allocate(iHam(sparseSize, nSpin))
    end if
    if (allocated(ERhoPrim)) then
      deallocate(ERhoPrim)
      allocate(ERhoPrim(sparseSize))
    end if

  end subroutine reallocateSparseArrays


  !> Calculates repulsive energy for current geometry
  subroutine calcRepulsiveEnergy(coord, species, img2CentCell, nNeighbor, neighborList,&
      & pRepCont, Eatom, Etotal)

    !> All atomic coordinates
    real(dp), intent(in) :: coord(:,:)

    !> All atoms chemical species 
    integer, intent(in) :: species(:)

    !> Image atom indices to central cell atoms
    integer, intent(in) :: img2CentCell(:)

    !> Number of neighbours for each actual atom
    integer, intent(in) :: nNeighbor(:)

    !> List of neighbours for each atom
    type(TNeighborList), intent(in) :: neighborList

    !> Repulsive interaction data
    type(ORepCont), intent(in) :: pRepCont

    !> Energy for each atom
    real(dp), intent(out) :: Eatom(:)

    !> Total energy
    real(dp), intent(out) :: Etotal

    call getERep(Eatom, coord, nNeighbor, neighborList%iNeighbor, species, pRepCont, img2CentCell)
    Etotal = sum(Eatom)

  end subroutine calcRepulsiveEnergy


  !> Calculates dispersion energy for current geometry.
  subroutine calcDispersionEnergy(dispersion, Eatom, Etotal)

    !> dispersion interactions
    class(DispersionIface), intent(inout) :: dispersion

    !> energy per atom
    real(dp), intent(out) :: Eatom(:)

    !> total energy
    real(dp), intent(out) :: Etotal

    call dispersion%getEnergies(Eatom)
    Etotal = sum(Eatom)

  end subroutine calcDispersionEnergy


  !> Sets the external potential components to zero
  subroutine resetExternalPotentials(potential)

    !> Potential contributions
    type(TPotentials), intent(inout) :: potential

    potential%extAtom(:,:) = 0.0_dp
    potential%extShell(:,:,:) = 0.0_dp
    potential%extBlock(:,:,:,:) = 0.0_dp

  end subroutine resetExternalPotentials


  !> Merges atomic and shell resolved external potentials into blocked one
  subroutine mergeExternalPotentials(orb, species, potential)

    !> Atomic orbital information
    type(TOrbitals), intent(in) :: orb

    !> species for atoms
    integer, intent(in) :: species(:)

    !> Potential energy contributions
    type(TPotentials), intent(inout) :: potential

    call total_shift(potential%extShell, potential%extAtom, orb, species)
    call total_shift(potential%extBlock, potential%extShell, orb, species)

  end subroutine mergeExternalPotentials


  !> Sets up electric external field
  subroutine setUpExternalElectricField(tTimeDepEField, tPeriodic, EFieldStrength, EFieldVector,&
      & EFieldOmega, EFieldPhase, neighborList, nNeighbor, iCellVec, img2CentCell, cellVec, deltaT,&
      & iGeoStep, coord0Fold, coord, EField, extAtomPot, absEField)

    !> Is there an electric field that varies with geometry step during MD?
    logical, intent(in) :: tTimeDepEField

    !> Is this a periodic geometry
    logical, intent(in) :: tPeriodic

    !> What is the field strength
    real(dp), intent(in) :: EFieldStrength

    !> What is the field direction
    real(dp), intent(in) :: EFieldVector(:)

    !> Is there an angular frequency for the applied field
    real(dp), intent(in) :: EFieldOmega

    !> What is the phase of the field
    integer, intent(in) :: EFieldPhase

    !> Atomi neighbours
    type(TNeighborList), intent(in) :: neighborList

    !> Number of neighbours for each atom
    integer, intent(in) :: nNeighbor(:)

    !> Index for unit cells
    integer, intent(in) :: iCellVec(:)

    !> Image atom to central cell atom number
    integer, intent(in) :: img2CentCell(:)

    !> Vectors to image unit cells
    real(dp), intent(in) :: cellVec(:,:)

    !> Time step in MD
    real(dp), intent(in) :: deltaT

    !> Number of the geometry step
    integer, intent(in) :: iGeoStep

    !> Atomic coordinates in central cell
    real(dp), allocatable, intent(in) :: coord0Fold(:,:)

    !> all coordinates
    real(dp), intent(in) :: coord(:,:)

    !> Resulting electric field
    real(dp), intent(out) :: EField(:)

    !> Potentials on atomic sites
    real(dp), intent(out) :: extAtomPot(:)

    !> Magnitude of the field
    real(dp), intent(out) :: absEField

    integer :: nAtom
    integer :: iAt1, iAt2, iNeigh
    character(lc) :: tmpStr

    nAtom = size(nNeighbor)

    Efield(:) = EFieldStrength * EfieldVector
    if (tTimeDepEField) then
      Efield(:) = Efield * sin(EfieldOmega * deltaT * real(iGeoStep + EfieldPhase, dp))
    end if
    absEfield = sqrt(sum(Efield**2))
    if (tPeriodic) then
      do iAt1 = 1, nAtom
        do iNeigh = 1, nNeighbor(iAt1)
          iAt2 = neighborList%iNeighbor(iNeigh, iAt1)
          ! overlap between atom in central cell and non-central cell
          if (iCellVec(iAt2) /= 0) then
            ! component of electric field projects onto vector between cells
            if (abs(dot_product(cellVec(:, iCellVec(iAt2)), EfieldVector)) > epsilon(1.0_dp)) then
              write(tmpStr, "(A, I0, A, I0, A)") 'Interaction between atoms ', iAt1, ' and ',&
                  & img2CentCell(iAt2),&
                  & ' crosses the saw-tooth discontinuity in the electric field.'
              call error(tmpStr)
            end if
          end if
        end do
      end do
      do iAt1 = 1, nAtom
        extAtomPot(iAt1) = dot_product(coord0Fold(:, iAt1), Efield)
      end do
    else
      do iAt1 = 1, nAtom
        extAtomPot(iAt1) = dot_product(coord(:, iAt1), Efield)
      end do
    end if

  end subroutine setUpExternalElectricField


  !> Initialise basic variables before the scc loop.
  subroutine initSccLoop(tScc, xlbomdIntegrator, minSccIter, maxSccIter, sccTol, tConverged)

    !> Is this an SCC calculation?
    logical, intent(in) :: tScc

    !> Details for extended Lagrange integrator (of used)
    type(Xlbomd), allocatable, intent(inout) :: xlbomdIntegrator

    !> Minimum number of SCC cycles that can be used
    integer, intent(inout) :: minSccIter

    !> Maximum number of SCC cycles
    integer, intent(inout) :: maxSccIter

    !> Tollerance for SCC convergence
    real(dp), intent(inout) :: sccTol

    !> Has SCC convergence been achieved?
    logical, intent(out) :: tConverged

    if (allocated(xlbomdIntegrator)) then
      call xlbomdIntegrator%getSCCParameters(minSCCIter, maxSccIter, sccTol)
    end if

    tConverged = (.not. tScc)

    if (tScc) then
      call printSccHeader()
    end if

  end subroutine initSccLoop


  !> Calculate the number of charges per shell.
  subroutine getChargePerShell(qq, orb, species, chargePerShell)

    !> charges in each orbital, for each atom and spin channel
    real(dp), intent(in) :: qq(:,:,:)

    !> orbital information
    type(TOrbitals), intent(in) :: orb

    !> species of each atom
    integer, intent(in) :: species(:)

    !> Resulting charges in atomic shells
    real(dp), intent(out) :: chargePerShell(:,:,:)

    integer :: iAt, iSp, iSh
    integer :: nAtom, nSpin

    nAtom = size(chargePerShell, dim=2)
    nSpin = size(chargePerShell, dim=3)
    chargePerShell(:,:,:) = 0.0_dp
    do iAt = 1, nAtom
      iSp = species(iAt)
      do iSh = 1, orb%nShell(iSp)
        chargePerShell(iSh, iAt, 1:nSpin) = chargePerShell(iSh, iAt, 1:nSpin)&
            & + sum(qq(orb%posShell(iSh, iSp) : orb%posShell(iSh + 1, iSp) - 1, iAt, 1:nSpin),&
            & dim=1)
      end do
    end do

  end subroutine getChargePerShell


  !> Reset internal potential related quantities
  subroutine resetInternalPotentials(tDualSpinOrbit, xi, orb, species, potential)

    !> Is dual spin orbit being used (block potentials)
    logical, intent(in) :: tDualSpinOrbit

    !> Spin orbit constants if required
    real(dp), allocatable, intent(in) :: xi(:,:)

    !> atomic orbital information
    type(TOrbitals), intent(in) :: orb

    !> chemical species
    integer, intent(in) :: species(:)

    !> potentials in the system
    type(TPotentials), intent(inout) :: potential

    @:ASSERT(tDualSpinOrbit .eqv. allocate(xi))

    potential%intAtom(:,:) = 0.0_dp
    potential%intShell(:,:,:) = 0.0_dp
    potential%intBlock(:,:,:,:) = 0.0_dp
    potential%orbitalBlock(:,:,:,:) = 0.0_dp
    potential%iOrbitalBlock(:,:,:,:) = 0.0_dp
    if (tDualSpinOrbit) then
      call shiftLS(potential%iOrbitalBlock, xi, orb, species)
    end if

  end subroutine resetInternalPotentials


  !> Add potentials comming from point charges.
  subroutine addChargePotentials(qInput, q0, chargePerShell, orb, species, species0, neighborList,&
      & img2CentCell, spinW, thirdOrd, potential)

    !> Input atomic populations
    real(dp), intent(in) :: qInput(:,:,:)

    !> reference atomic occupations
    real(dp), intent(in) :: q0(:,:,:)

    !> charges per atomic shell
    real(dp), intent(in) :: chargePerShell(:,:,:)

    !> atomic orbital information
    type(TOrbitals), intent(in) :: orb

    !> species of all atoms
    integer, intent(in) :: species(:)
    
    !> species of central cell atoms
    integer, intent(in) :: species0(:)

    !> neighbours to atoms
    type(TNeighborList), intent(in) :: neighborList

    !> map from image atom to real atoms
    integer, intent(in) :: img2CentCell(:)

    !> spin constants
    real(dp), intent(in), allocatable :: spinW(:,:,:)

    !> third order SCC interactions
    type(ThirdOrder), allocatable, intent(inout) :: thirdOrd

    !> Potentials acting
    type(TPotentials), intent(inout) :: potential

    real(dp), allocatable :: atomPot(:,:)
    real(dp), allocatable :: shellPot(:,:,:)
    integer :: nAtom, nSpin

    nAtom = size(qInput, dim=2)
    nSpin = size(qInput, dim=3)

    allocate(atomPot(nAtom, nSpin))
    allocate(shellPot(orb%mShell, nAtom, nSpin))

    call updateCharges_SCC(qInput, q0, orb, species, neighborList%iNeighbor, img2CentCell)
    call getShiftPerAtom(atomPot(:,1))
    call getShiftPerL(shellPot(:,:,1))
    potential%intAtom(:,1) = potential%intAtom(:,1) + atomPot(:,1)
    potential%intShell(:,:,1) = potential%intShell(:,:,1) + shellPot(:,:,1)

    if (allocated(thirdOrd)) then
      call thirdOrd%updateCharges(species0, neighborList, qInput, q0, img2CentCell, orb)
      call thirdOrd%getShifts(atomPot(:,1), shellPot(:,:,1))
      potential%intAtom(:,1) = potential%intAtom(:,1) + atomPot(:,1)
      potential%intShell(:,:,1) = potential%intShell(:,:,1) + shellPot(:,:,1)
    end if

    if (allocated(spinW)) then
      call getSpinShift(shellPot, chargePerShell, species, orb, spinW)
      potential%intShell = potential%intShell + shellPot
    end if

    call total_shift(potential%intShell, potential%intAtom, orb, species)
    call total_shift(potential%intBlock, potential%intShell, orb, species)

  end subroutine addChargePotentials


  !> Add potentials comming from on-site block of the dual density matrix.
  subroutine addBlockChargePotentials(qBlockIn, qiBlockIn, tDftbU, tImHam, species, orb,&
      & nDftbUFunc, UJ, nUJ, iUJ, niUJ, potential)

    !> block input charges
    real(dp), allocatable, intent(in) :: qBlockIn(:,:,:,:)

    !> imaginary part
    real(dp), allocatable, intent(in) :: qiBlockIn(:,:,:,:)

    !> is this a +U calculation
    logical, intent(in) :: tDftbU

    !> does the hamitonian have an imaginary part in real space?
    logical, intent(in) :: tImHam

    !> chemical species of all atoms
    integer, intent(in) :: species(:)

    !> Orbital information
    type(TOrbitals), intent(in) :: orb

    !> choice of +U functional 
    integer, intent(in) :: nDftbUFunc

    !> prefactor for +U potential
    real(dp), allocatable, intent(in) :: UJ(:,:)

    !> Number DFTB+U blocks of shells for each atom type
    integer, intent(in), optional :: nUJ(:)

    !> which shells are in each DFTB+U block
    integer, intent(in), optional :: iUJ(:,:,:)

    !> Number of shells in each DFTB+U block
    integer, intent(in), optional :: niUJ(:,:)

    !> potentials acting in system
    type(TPotentials), intent(inout) :: potential


    if (tDFTBU) then
      if (tImHam) then
        call getDftbUShift(potential%orbitalBlock, potential%iorbitalBlock, qBlockIn, qiBlockIn,&
            & species,orb, nDFTBUfunc, UJ, nUJ, niUJ, iUJ)
      else
        call getDftbUShift(potential%orbitalBlock, qBlockIn, species, orb, nDFTBUfunc, UJ, nUJ,&
            & niUJ, iUJ)
      end if
      potential%intBlock = potential%intBlock + potential%orbitalBlock
    end if

  end subroutine addBlockChargePotentials



  !> Returns the Hamiltonian for the given scc iteration
  subroutine getSccHamiltonian(H0, over, nNeighbor, neighborList, species, orb, iPair,&
      & img2CentCell, potential, ham, iHam)

    !> non-SCC hamitonian (sparse)
    real(dp), intent(in) :: H0(:)

    !> overlap (sparse)
    real(dp), intent(in) :: over(:)

    !> Number of atomic neighbours
    integer, intent(in) :: nNeighbor(:)

    !> list of atomic neighbours
    type(TNeighborList), intent(in) :: neighborList

    !> species of atoms
    integer, intent(in) :: species(:)

    !> atomic orbital information
    type(TOrbitals), intent(in) :: orb

    !> Index for atomic blocks in sparse data
    integer, intent(in) :: iPair(:,:)

    !> image atoms to central cell atoms
    integer, intent(in) :: img2CentCell(:)

    !> potential acting on sustem
    type(TPotentials), intent(in) :: potential

    !> resulting hamitonian (sparse)
    real(dp), intent(out) :: ham(:,:)

    !> imaginary part of hamitonian (if required, signalled by being allocated)
    real(dp), allocatable, intent(inout) :: iHam(:,:)

    integer :: nAtom

    nAtom = size(orb%nOrbAtom)

    ham(:,:) = 0.0_dp
    ham(:,1) = h0
    call add_shift(ham, over, nNeighbor, neighborList%iNeighbor, species, orb, iPair, nAtom,&
        & img2CentCell, potential%intBlock)

    if (allocated(iHam)) then
      iHam(:,:) = 0.0_dp
      call add_shift(iHam, over, nNeighbor, neighborList%iNeighbor, species, orb, iPair, nAtom,&
          & img2CentCell, potential%iorbitalBlock)
    end if

  end subroutine getSccHamiltonian


  !> Writes Hamiltonian and overlap matrices and stops program execution.
  subroutine writeHSAndStop(tWriteHS, tWriteRealHS, tRealHS, over, neighborList, nNeighbor,&
      & iAtomStart, iPair, img2CentCell, kPoint, iCellVec, cellVec, ham, iHam)

    !> Write dense hamiltonian and overlap matrices
    logical, intent(in) :: tWriteHS

    !> write sparse hamitonian and overlap matrices
    logical, intent(in) :: tWriteRealHS

    !> Is the hamitonian real?
    logical, intent(in) :: tRealHS

    !> overlap in sparse storage
    real(dp), intent(in) :: over(:)

    !> atomic neighbours
    type(TNeighborList), intent(in) :: neighborList

    !> number of neighbours for each central cell atom
    integer, intent(in) :: nNeighbor(:)

    !> dense matrix indexing for atomic blocks
    integer, intent(in) :: iAtomStart(:)

    !> sparse matrix indexing for atomic blocks
    integer, intent(in) :: iPair(:,:)

    !> Image atoms to central cell
    integer, intent(in) :: img2CentCell(:)

    !> k-points
    real(dp), intent(in) :: kPoint(:,:)

    !> index  for which unit cell an atom is in
    integer, intent(in) :: iCellVec(:)

    !> vectors to unit cells, in lattice constant units
    real(dp), intent(in) :: cellVec(:,:)

    !> sparse hamitonian
    real(dp), intent(inout) :: ham(:,:)

    !> imaginary part of hamitonian (used if allocated)
    real(dp), allocatable, intent(inout) :: iHam(:,:)

    integer :: nSpin

    nSpin = size(ham, dim=2)

    ! Sanity check, although this should have been caught in initprogram already.
    if (nSpin == 4) then
      call error('Internal error: Hamiltonian writing for Pauli-Hamiltoninan not implemented')
    end if

    if (nSpin == 2) then
      call qm2ud(ham)
    end if

    ! Write out matrices if necessary and quit.
    if (allocated(iHam)) then
      call writeHS(tWriteHS, tWriteRealHS, tRealHS, ham, over, neighborList%iNeighbor, nNeighbor,&
          & iAtomStart, iPair, img2CentCell, kPoint, iCellVec, cellVec, iHam=iHam)
    else
      call writeHS(tWriteHS, tWriteRealHS, tRealHS, ham, over, neighborList%iNeighbor, nNeighbor,&
          & iAtomStart, iPair, img2CentCell, kPoint, iCellVec, cellVec)
    end if
    write(stdOut, "(A)") "Hamilton/Overlap written, exiting program."
    stop

  end subroutine writeHSAndStop


  !> Returns the sparse density matrix.
  !>
  !> All operations (e.g. non-dual spin orbit coupling), which need access to full (unpacked)
  !> Hamiltonian or the full (unpacked) density matrix, must also invoked from within this routine,
  !> as those unpacked quantities do not exist elsewhere.
  !>
  subroutine getDensity(ham, over, neighborList, nNeighbor, iAtomStart, iPair, img2CentCell, &
      & iCellVec, cellVec, kPoint, kWeight, orb, species, solver, tRealHS, tSpinSharedEf,&
      & tSpinOrbit, tDualSpinOrbit, tFillKSep, tFixEf, tMulliken, iDistribFn, tempElec, nEl, Ef, &
      & energy, eigen, filling, rhoPrim, Eband, TS, E0, iHam, xi, orbitalL, HSqrReal, SSqrReal, &
      & iRhoPrim, HSqrCplx, SSqrCplx, rhoSqrReal, storeEigvecsReal, storeEigvecsCplx)
    real(dp), intent(inout) :: ham(:,:)
    real(dp), intent(in) :: over(:)
    type(TNeighborList), intent(in) :: neighborList
    integer, intent(in) :: nNeighbor(:)
    integer, intent(in) :: iAtomStart(:)
    integer, intent(in) :: iPair(:,:)
    integer, intent(in) :: img2CentCell(:)
    integer, intent(in) :: iCellVec(:)
    real(dp), intent(in) :: cellVec(:,:)
    real(dp), intent(in) :: kPoint(:,:)
    real(dp), intent(in) :: kWeight(:)
    type(TOrbitals), intent(in) :: orb
    integer, intent(in) :: species(:)
    integer, intent(in) :: solver
    logical, intent(in) :: tRealHS
    logical, intent(in) :: tSpinSharedEf
    logical, intent(in) :: tSpinOrbit
    logical, intent(in) :: tDualSpinOrbit
    logical, intent(in) :: tFillKSep
    logical, intent(in) :: tFixEf    
    logical, intent(in) :: tMulliken
    integer, intent(in) :: iDistribFn
    real(dp), intent(in) :: tempElec, nEl(:)
    real(dp), intent(inout) :: Ef(:)
    type(TEnergies), intent(inout) :: energy
    real(dp), intent(out) :: eigen(:,:,:)
    real(dp), intent(out) :: filling(:,:,:)
    real(dp), intent(out) :: rhoPrim(:,:)
    real(dp), intent(out) :: Eband(:)
    real(dp), intent(out) :: TS(:)
    real(dp), intent(out) :: E0(:)
    real(dp), intent(inout), optional :: iHam(:,:)
    real(dp), intent(in), optional :: xi(:,:)
    real(dp), intent(out), optional :: orbitalL(:,:,:)
    real(dp), intent(out), optional :: iRhoPrim(:,:)
    real(dp), intent(out), optional :: HSqrReal(:,:,:)
    real(dp), intent(out), optional :: SSqrReal(:,:)
    complex(dp), intent(out), optional :: HSqrCplx(:,:,:,:)
    complex(dp), intent(out), optional :: SSqrCplx(:,:)
    real(dp), intent(out), optional :: rhoSqrReal(:,:,:)
    type(OFifoRealR2), intent(inout), optional :: storeEigvecsReal(:)
    type(OFifoCplxR2), intent(inout), optional :: storeEigvecsCplx(:)

    integer :: nSpin

    nSpin = size(ham, dim=2)

    ! Hack due to not using Pauli-type structure for diagonalisation
    if (nSpin > 1) then
      ham(:,:) = 2.0_dp * ham
      if (present(iHam)) then
        iHam(:,:) = 2.0_dp * iHam
      end if
    end if

    if (nSpin /= 4) then
      call qm2ud(ham)
      if (tRealHS) then
        call buildAndDiagDenseHam(ham, over, neighborList, nNeighbor, iAtomStart, iPair,&
            & img2CentCell, solver, HSqrReal, SSqrReal, eigen(:,1,:))
      else
        call buildAndDiagDenseKDepHam(ham, over, kPoint, neighborList, nNeighbor, iAtomStart,&
            & iPair, img2CentCell, iCellVec, cellVec, solver, HSqrCplx, SSqrCplx, eigen)
      end if
    else
      if (tRealHS) then
        call buildAndDiagDensePauliHam(ham, over, neighborList, nNeighbor, iAtomStart, iPair,&
            & img2CentCell, solver, eigen(:,1,1), HSqrCplx(:,:,1,1), SSqrCplx, iHam, xi, species,&
            & orb, storeEigvecsCplx)
      else
        call buildAndDiagDenseKDepPauliHam(ham, over, kPoint, neighborList, nNeighbor, iAtomStart,&
            & iPair, img2CentCell, iCellVec, cellVec, solver, eigen(:,:,1), HSqrCplx(:,:,:,1),&
            & SSqrCplx, iHam, xi, species, orb, storeEigvecsCplx)
      end if
    end if

    call getFillingsAndBandEnergies(eigen, nEl, nSpin, tempElec, kWeight, tSpinSharedEf,&
        & tFillKSep, tFixEf, iDistribFn, Ef, filling, Eband, TS, E0)

    if (nSpin /= 4) then
      if (tRealHS) then
        call getDensityFromEigvecs(filling(:,1,:), neighborList, nNeighbor, iPair, iAtomStart,&
            & img2CentCell, orb, HSqrReal, rhoPrim, SSqrReal, rhoSqrReal, storeEigvecsReal)
      else
        call getDensityFromKDepEigvecs(filling, kPoint, kWeight, neighborList, nNeighbor, iPair,&
            & iAtomStart, img2CentCell, iCellVec, cellVec, orb, HSqrCplx, rhoPrim, SSqrCplx,&
            & storeEigvecsCplx)
      end if
      call ud2qm(rhoPrim)
    else
      ! Pauli structure of eigenvectors
      filling(:,:,1) = 2.0_dp * filling(:,:,1)
      call getDensityFromPauliEigvecs(tRealHS, tSpinOrbit, tDualSpinOrbit, tMulliken, kPoint,&
          & kWeight, filling(:,:,1), neighborList, nNeighbor, orb, iAtomStart, iPair,&
          & img2CentCell, iCellVec, cellVec, species, HSqrCplx(:,:,:,1), SSqrCplx, energy, rhoPrim,&
          & xi, orbitalL, iRhoPrim, storeEigvecsCplx)
      filling(:,:,1) = 0.5_dp * filling(:,:,1)
    end if

  end subroutine getDensity


  !> Builds and diagonalises dense Hamiltonians.
  subroutine buildAndDiagDenseHam(ham, over, neighborList, nNeighbor, iAtomStart, iPair,&
      & img2CentCell, solver, HSqrReal, SSqrReal, eigen, storeEigvecsReal)
    real(dp), intent(in) :: ham(:,:)
    real(dp), intent(in) :: over(:)
    type(TNeighborList), intent(in) :: neighborList
    integer, intent(in) :: nNeighbor(:)
    integer, intent(in) :: iAtomStart(:)
    integer, intent(in) :: iPair(:,:)
    integer, intent(in) :: img2CentCell(:)
    integer, intent(in) :: solver
    real(dp), intent(out) :: HSqrReal(:,:,:)
    real(dp), intent(out) :: SSqrReal(:,:)
    real(dp), intent(out) :: eigen(:,:)
    type(OFifoRealR2), intent(inout), optional :: storeEigvecsReal(:)

    logical :: tStoreEigvecs
    integer :: nSpin
    integer :: iSpin, iSpin2

    nSpin = size(ham, dim=2)
    tStoreEigvecs = present(storeEigvecsReal)

    #:call ASSERT_CODE
      @:ASSERT(size(HSqrReal, dim=3) == nSpin .or. (tStoreEigvecs .and. size(HSqrReal, dim=3) == 1))
      if (tStoreEigvecs) then
        @:ASSERT(size(storeEigvecsReal) == nSpin)
      end if
    #:endcall ASSERT_CODE

    do iSpin = 1, nSpin
      if (tStoreEigvecs) then
        iSpin2 = 1
      else
        iSpin2 = iSpin
      end if
      call unpackHS(HSqrReal(:,:,iSpin2), ham(:,iSpin), neighborList%iNeighbor, nNeighbor,&
          & iAtomStart, iPair, img2CentCell)
      call unpackHS(SSqrReal, over, neighborList%iNeighbor, nNeighbor, iAtomStart, iPair,&
          & img2CentCell)
      call diagonalizeDenseMtx(solver, 'V', HSqrReal(:,:,iSpin2), SSqrReal, eigen(:,iSpin))
      if (tStoreEigvecs) then
        call reset(storeEigvecsReal(iSpin), [size(HSqrReal, dim=1), size(HSqrReal, dim=2)])
        call push(storeEigvecsReal(iSpin), HSqrReal(:,:,iSpin2))
      end if
    end do

  end subroutine buildAndDiagDenseHam


  !> Builds and diagonalises dense k-point dependent Hamiltonians.
  subroutine buildAndDiagDenseKDepHam(ham, over, kPoint, neighborList, nNeighbor, iAtomStart,&
      & iPair, img2CentCell, iCellVec, cellVec, solver, HSqrCplx, SSqrCplx, eigen, storeEigvecsCplx)
    real(dp), intent(in) :: ham(:,:)
    real(dp), intent(in) :: over(:)
    real(dp), intent(in) :: kPoint(:,:)
    type(TNeighborList), intent(in) :: neighborList
    integer, intent(in) :: nNeighbor(:)
    integer, intent(in) :: iAtomStart(:)
    integer, intent(in) :: iPair(:,:)
    integer, intent(in) :: img2CentCell(:)
    integer, intent(in) :: iCellVec(:)
    real(dp), intent(in) :: cellVec(:,:)
    integer, intent(in) :: solver
    complex(dp), intent(out) :: HSqrCplx(:,:,:,:)
    complex(dp), intent(out) :: SSqrCplx(:,:)
    real(dp), intent(out) :: eigen(:,:,:)
    type(OFifoCplxR2), intent(inout), optional :: storeEigvecsCplx(:)

    integer :: nSpin, nKPoint
    integer :: iK, iK2, iSpin, iSpin2
    logical :: tStoreEigvecs

    nSpin = size(ham, dim=2)
    nKPoint = size(kPoint, dim=2)
    tStoreEigvecs = present(storeEigvecsCplx)

    #:call ASSERT_CODE
      @:ASSERT((size(HSqrCplx, dim=2) == nKPoint .and. size(HSqrReal, dim=3) == nSpin)&
          & .or. (tStoreEigvecs .and. size(HSqrReal, dim=2) == 1 .and. size(HSqrReal, dim=3) == 1))
      if (tStoreEigvecs) then
        @:ASSERT(size(storeEigvecsCplx) == nSpin)
      end if
    #:endcall ASSERT_CODE

    do iSpin = 1, nSpin
      do iK = 1, nKPoint
        if (tStoreEigvecs) then
          iSpin2 = 1
          iK2 = 1
        else
          iSpin2 = iSpin
          iK2 = iK
        end if
        call unpackHS(HSqrCplx(:,:,iK2,iSpin2), ham(:,iSpin), kPoint(:,iK), neighborList%iNeighbor,&
            & nNeighbor, iCellVec, cellVec, iAtomStart, iPair, img2CentCell)
        call unpackHS(SSqrCplx, over, kPoint(:,iK), neighborList%iNeighbor, nNeighbor, iCellVec,&
            & cellVec, iAtomStart, iPair, img2CentCell)
        call diagonalizeDenseMtx(solver, 'V', HSqrCplx(:,:,iK2,iSpin2), SSqrCplx, eigen(:,iK,iSpin))
        if (tStoreEigvecs) then
          call reset(storeEigvecsCplx(iSpin), [size(HSqrCplx, dim=1), size(HSqrCplx, dim=2)])
          call push(storeEigvecsCplx(iSpin), HSqrCplx(:,:,iK2,iSpin2))
        end if
      end do
    end do

  end subroutine buildAndDiagDenseKDepHam


  !> Builds and diagonalizes Pauli two-component Hamiltonians.
  subroutine buildAndDiagDensePauliHam(ham, over, neighborList, nNeighbor, iAtomStart, iPair,&
      & img2CentCell, solver, eigen, HSqrCplx, SSqrCplx, iHam, xi, species, orb, storeEigvecsCplx)
    real(dp), intent(in) :: ham(:,:)
    real(dp), intent(in) :: over(:)
    type(TNeighborList), intent(in) :: neighborList
    integer, intent(in) :: nNeighbor(:)
    integer, intent(in) :: iAtomStart(:)
    integer, intent(in) :: iPair(:,:)
    integer, intent(in) :: img2CentCell(:)
    integer, intent(in) :: solver
    real(dp), intent(out) :: eigen(:)
    complex(dp), intent(out) :: HSqrCplx(:,:)
    complex(dp), intent(out) :: SSqrCplx(:,:)
    real(dp), intent(in), optional :: iHam(:,:)
    real(dp), intent(in), optional :: xi(:,:)
    integer, intent(in), optional :: species(:)
    type(TOrbitals), intent(in), optional :: orb
    type(OFifoCplxR2), intent(inout), optional :: storeEigvecsCplx(:)

    logical :: tStoreEigvecs

    tStoreEigvecs = present(storeEigvecsCplx)

    #:call ASSERT_CODE
      @:ASSERT(.not. present(xi) .or. (present(species) .and. present(orb)))
      if (tStoreEigvecs) then
        @:ASSERT(size(storeEigvecsCplx) == 1)
      end if
    #:endcall ASSERT_CODE

    call unpackHSPauli(ham, over, neighborList%iNeighbor, nNeighbor, iPair, iAtomStart,&
        & img2CentCell, HSqrCplx, SSqrCplx, iHam=iHam)
    if (present(xi) .and. .not. present(iHam)) then
      call addOnsiteSpinOrbitContrib(xi, species, orb, iAtomStart, HSqrCplx)
    end if
    call diagonalizeDenseMtx(solver, 'V', HSqrCplx, SSqrCplx, eigen)
    if (tStoreEigvecs) then
      call reset(storeEigvecsCplx(1), shape(HSqrCplx))
      call push(storeEigvecsCplx(1), HSqrCplx)
    end if

  end subroutine buildAndDiagDensePauliHam


  !> Builds and diagonalizes k-dependent Pauli two-component Hamiltonians.
  subroutine buildAndDiagDenseKDepPauliHam(ham, over, kPoint, neighborList, nNeighbor, iAtomStart,&
      & iPair, img2CentCell, iCellVec, cellVec, solver, eigen, HSqrCplx, SSqrCplx, iHam, xi,&
      & species, orb, storeEigvecsCplx)
    real(dp), intent(in) :: ham(:,:)
    real(dp), intent(in) :: over(:)
    real(dp), intent(in) :: kPoint(:,:)
    type(TNeighborList), intent(in) :: neighborList
    integer, intent(in) :: nNeighbor(:)
    integer, intent(in) :: iAtomStart(:)
    integer, intent(in) :: iPair(:,:)
    integer, intent(in) :: img2CentCell(:)
    integer, intent(in) :: iCellVec(:)
    real(dp), intent(in) :: cellVec(:,:)
    integer, intent(in) :: solver
    real(dp), intent(out) :: eigen(:,:)
    complex(dp), intent(out) :: HSqrCplx(:,:,:)
    complex(dp), intent(out) :: SSqrCplx(:,:)
    real(dp), intent(in), optional :: iHam(:,:)
    real(dp), intent(in), optional :: xi(:,:)
    integer, intent(in), optional :: species(:)
    type(TOrbitals), intent(in), optional :: orb
    type(OFifoCplxR2), intent(inout), optional :: storeEigvecsCplx(:)

    integer :: nKPoint
    integer :: iK, iK2
    logical :: tStoreEigvecs

    nKPoint = size(kPoint, dim=2)
    tStoreEigvecs = present(storeEigvecsCplx)

    #:call ASSERT_CODE
      @:ASSERT(.not. present(xi) .or. (present(species) .and. present(orb)))
      if (tStoreEigvecs) then
        @:ASSERT(size(storeEigvecsCplx) == 1)
      end if
    #:endcall ASSERT_CODE

    do iK = 1, nKPoint
      if (tStoreEigvecs) then
        iK2 = 1
      else
        iK2 = iK
      end if
      call unpackHSPauliK(ham, over, kPoint(:,iK), neighborList%iNeighbor, nNeighbor,&
          & iAtomStart, iPair, img2CentCell, iCellVec, cellVec, HSqrCplx(:,:,iK2),&
          & SSqrCplx, iHam=iHam)
      if (present(xi) .and. .not. present(iHam)) then
        call addOnsiteSpinOrbitContrib(xi, species, orb, iAtomStart, HSqrCplx(:,:,iK))
      end if
      call diagonalizeDenseMtx(solver, 'V', HSqrCplx(:,:,iK2), SSqrCplx, eigen(:,iK))
      if (tStoreEigvecs) then
        call reset(storeEigvecsCplx(1), [size(HSqrCplx, dim=1), size(HSqrCplx, dim=2)])
        call push(storeEigvecsCplx(1), HSqrCplx(:,:,iK2))
      end if
    end do

  end subroutine buildAndDiagDenseKDepPauliHam


  !> Creates sparse density matrix from real eigenvectors.
  subroutine getDensityFromEigvecs(filling, neighborList, nNeighbor, iPair, iAtomStart,&
      & img2CentCell, orb, eigvecs, rhoPrim, work, rhoSqrReal, storeEigvecsReal)
    real(dp), intent(in) :: filling(:,:)
    type(TNeighborList), intent(in) :: neighborList
    integer, intent(in) :: nNeighbor(:)
    integer, intent(in) :: iPair(:,:)
    integer, intent(in) :: iAtomStart(:)
    integer, intent(in) :: img2CentCell(:)
    type(TOrbitals), intent(in) :: orb
    real(dp), intent(inout) :: eigvecs(:,:,:)
    real(dp), intent(out) :: rhoPrim(:,:)
    real(dp), intent(out) :: work(:,:)
    real(dp), intent(out), optional  :: rhoSqrReal(:,:,:)
    type(OFifoRealR2), intent(inout), optional :: storeEigvecsReal(:)

    integer :: nSpin
    integer :: iSpin, iSpin2
    logical :: tStoreEigvecs

    nSpin = size(filling, dim=2)
    tStoreEigvecs = present(storeEigvecsReal)

    rhoPrim(:,:) = 0.0_dp
    do iSpin = 1, nSpin
      if (tStoreEigvecs) then
        iSpin2 = 1
        call get(storeEigvecsReal(iSpin), eigvecs(:,:,iSpin2))
      else
        iSpin2 = iSpin
      end if

      if (tDensON2) then
        call makeDensityMatrix(work, eigvecs(:,:,iSpin2), filling(:,iSpin),&
            & neighborlist%iNeighbor, nNeighbor, orb, iAtomStart, img2CentCell)
      else
        call makeDensityMatrix(work, eigvecs(:,:,iSpin2), filling(:,iSpin))
      end if
      call packHS(rhoPrim(:,iSpin), work, neighborlist%iNeighbor, nNeighbor, orb%mOrb, iAtomStart,&
          & iPair, img2CentCell)

      if (present(rhoSqrReal)) then
        rhoSqrReal(:,:,iSpin) = work
      end if
    end do

  end subroutine getDensityFromEigvecs


  !> Creates sparse density matrix from complex eigenvectors.
  subroutine getDensityFromKDepEigvecs(filling, kPoint, kWeight, neighborList, nNeighbor, iPair,&
      & iAtomStart, img2CentCell, iCellVec, cellVec, orb, eigvecs, rhoPrim, work, storeEigvecsCplx)
    real(dp), intent(in) :: filling(:,:,:)
    real(dp), intent(in) :: kPoint(:,:)
    real(dp), intent(in) :: kWeight(:)
    type(TNeighborList), intent(in) :: neighborList
    integer, intent(in) :: nNeighbor(:)
    integer, intent(in) :: iPair(:,:)
    integer, intent(in) :: iAtomStart(:)
    integer, intent(in) :: img2CentCell(:)
    integer, intent(in) :: iCellVec(:)
    real(dp), intent(in) :: cellVec(:,:)
    type(TOrbitals), intent(in) :: orb
    complex(dp), intent(inout) :: eigvecs(:,:,:,:)
    real(dp), intent(out) :: rhoPrim(:,:)
    complex(dp), intent(out) :: work(:,:)
    type(OFifoCplxR2), intent(inout), optional :: storeEigvecsCplx(:)

    integer :: nSpin
    integer :: iK, iK2, iSpin, iSpin2
    logical :: tStoreEigvecs

    nSpin = size(filling, dim=3)
    tStoreEigvecs = present(storeEigvecsCplx)

    rhoPrim(:,:) = 0.0_dp

    do iSpin = 1, nSpin
      do iK = 1, size(kPoint, dim=2)
        if (tStoreEigvecs) then
          iSpin2 = 1
          iK2 = 1
          call get(storeEigvecsCplx(iSpin), eigvecs(:,:,iK2,iSpin2))
        else
          iSpin2 = iSpin
          iK2 = iK
        end if

        if (tDensON2) then
          call makeDensityMatrix(work, eigvecs(:,:,iK2,iSpin2), filling(:,iK,iSpin),&
              & neighborlist%iNeighbor, nNeighbor, orb, iAtomStart, img2CentCell)
        else
          call makeDensityMatrix(work, eigvecs(:,:,iK2,iSpin2), filling(:,iK,iSpin))
        end if
        call packHS(rhoPrim(:,iSpin), work, kPoint(:,iK), kWeight(iK), neighborList%iNeighbor,&
            & nNeighbor, orb%mOrb, iCellVec, cellVec, iAtomStart, iPair, img2CentCell)
      end do
    end do

  end subroutine getDensityFromKDepEigvecs


  !> Creates sparse density matrix from two component complex eigenvectors.
  subroutine getDensityFromPauliEigvecs(tRealHS, tSpinOrbit, tDualSpinOrbit, tMulliken,&
      & kPoint, kWeight, filling, neighborList, nNeighbor, orb, iAtomStart, iPair,&
      & img2CentCell, iCellVec, cellVec, species, eigvecs, work, &
      & energy, rhoPrim, xi, orbitalL, iRhoPrim, eigvecsFifo)
    logical, intent(in) :: tRealHS
    logical, intent(in) :: tSpinOrbit
    logical, intent(in) :: tDualSpinOrbit
    logical, intent(in) :: tMulliken
    real(dp), intent(in) :: kPoint(:,:)
    real(dp), intent(in) :: kWeight(:)
    real(dp), intent(in) :: filling(:,:)
    type(TNeighborList), intent(in) :: neighborList
    integer, intent(in) :: nNeighbor(:)
    type(TOrbitals), intent(in) :: orb
    integer, intent(in) :: iAtomStart(:)
    integer, intent(in) :: iPair(:,:)
    integer, intent(in) :: img2CentCell(:)
    integer, intent(in) :: iCellVec(:)
    real(dp), intent(in) :: cellVec(:,:)
    integer, intent(in) :: species(:)
    complex(dp), intent(inout) :: eigvecs(:,:,:)
    complex(dp), intent(inout) :: work(:,:)
    type(TEnergies), intent(inout) :: energy
    real(dp), intent(out) :: rhoPrim(:,:)
    real(dp), intent(in), optional :: xi(:,:)
    real(dp), intent(out), optional :: orbitalL(:,:,:)
    real(dp), intent(out), optional :: iRhoPrim(:,:)
    type(OFifoCplxR2), intent(inout), optional :: eigvecsFifo(:)


    real(dp), allocatable :: rVecTemp(:), orbitalLPart(:,:,:)
    integer :: nKPoint, nAtom
    integer :: iK, iK2
    logical :: tStoreEigvecs, tImHam

    nAtom = size(orb%nOrbAtom)
    tImHam = present(iRhoPrim)
    tStoreEigvecs = present(eigvecsFifo)
    nKPoint = size(kWeight)

    rhoPrim(:,:) = 0.0_dp
    if (present(iRhoPrim)) then
      iRhoPrim(:,:) = 0.0_dp
    end if

    work(:,:) = 0.0_dp
    if (tSpinOrbit .and. .not. tDualSpinOrbit) then
      energy%atomLS(:) = 0.0_dp
      allocate(rVecTemp(nAtom))
    end if

    if (tMulliken .and. tSpinOrbit .and. .not. tDualSpinOrbit) then
      allocate(orbitalLPart(3, orb%mShell, nAtom))
      orbitalL(:,:,:) = 0.0_dp
    end if

    do iK = 1, nKPoint
      if (tStoreEigvecs) then
        iK2 = 1
        call get(eigvecsFifo(1), eigvecs(:,:,iK2))
      else
        iK2 = iK
      end if
      call makeDensityMatrix(work, eigvecs(:,:,iK2), filling(:,iK))
      if (tSpinOrbit .and. .not. tDualSpinOrbit) then
        rVecTemp(:) = 0.0_dp
        call getEnergySpinOrbit(rVecTemp, work, iAtomStart, xi, orb, species)
        energy%atomLS = energy%atomLS + kWeight(iK) * rVecTemp
        if (tMulliken) then
          orbitalLPart(:,:,:) = 0.0_dp
          call getL(orbitalLPart, work, iAtomStart, orb, species)
          orbitalL(:,:,:) = orbitalL + kWeight(iK) * orbitalLPart
        end if
      end if

      if (tRealHS) then
        call packHSPauli(rhoPrim, work, neighborlist%iNeighbor, nNeighbor, orb%mOrb, iAtomStart,&
            & iPair, img2CentCell)
        if (tImHam) then
          call packHSPauliImag(iRhoPrim, work, neighborlist%iNeighbor, nNeighbor, orb%mOrb,&
              & iAtomStart, iPair, img2CentCell)
        end if
      else
        call packHS(rhoPrim, work, kPoint(:,iK), kWeight(iK), neighborList%iNeighbor, nNeighbor,&
            & orb%mOrb, iCellVec, cellVec, iAtomStart, iPair, img2CentCell)
        if (tImHam) then
          call iPackHS(iRhoPrim, work, kPoint(:,iK), kWeight(iK), neighborlist%iNeighbor, &
              & nNeighbor, orb%mOrb, iCellVec, cellVec, iAtomStart, iPair, img2CentCell)
        end if
      end if
    end do

    if (tSpinOrbit .and. .not. tDualSpinOrbit) then
      energy%ELS = sum(energy%atomLS)
    end if

  end subroutine getDensityFromPauliEigvecs


  !> Calculates electron fillings and resulting band energy terms.
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
    logical, intent(in) :: tFillKSep

    !> Whether fixed Fermi level(s) should be used. (No charge conservation!)
    logical, intent(in) :: tFixEf

    !> Selector for the distribution function
    integer, intent(in) :: iDistribFn

    !> Fixed Fermi levels on entry, if tFixEf is .true., otherwise the Fermi levels found for the
    !> given number of electrons on exit
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
      ! Fixed value of the Fermi level for each spin channel
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


  !> Adds spin-orbit contribution to dense Hamiltonian (for non-dual spin-orbit model).
  subroutine addOnsiteSpinOrbitContrib(xi, species, orb, iAtomStart, HSqrCplx)

    !> Spin orbit constants for each species
    real(dp), intent(in) :: xi(:,:)

    !> chemical species
    integer, intent(in) :: species(:)

    !> atomic orbital information
    type(TOrbitals), intent(in) :: orb

    !> index array for atomic blocks in dense matrices
    integer, intent(in) :: iAtomStart(:)

    !> Dense hamitonian matrix (2 component)
    complex(dp), intent(inout) :: HSqrCplx(:,:)

    complex(dp), allocatable :: atomZ(:,:,:), atomPlus(:,:,:), Lz(:,:), Lplus(:,:)
    integer :: nAtom, nSpecies, nOrb
    integer :: iSp, iOrb, iAt, ll

    nAtom = size(orb%nOrbAtom)
    nSpecies = maxval(species(1:nAtom))
    nOrb = size(HSqrCplx, dim=1) / 2
    allocate(atomZ(orb%mOrb, orb%mOrb, nSpecies))
    allocate(atomPlus(orb%mOrb, orb%mOrb, nSpecies))
    atomZ(:,:,:) = 0.0_dp
    atomPlus(:,:,:) = 0.0_dp
    allocate(Lz(orb%mOrb, orb%mOrb))
    allocate(Lplus(orb%mOrb, orb%mOrb))

    do iSp = 1, nSpecies
      do iOrb = 1, orb%nShell(iSp)
        Lz(:,:) = 0.0_dp
        Lplus(:,:) = 0.0_dp
        ll = orb%angShell(iOrb,iSp)
        call loperators(Lplus(1:2*ll+1,1:2*ll+1), Lz(1:2*ll+1,1:2*ll+1), ll)
        atomZ(orb%posShell(iOrb, iSp) : orb%posShell(iOrb + 1, iSp) - 1, &
            & orb%posShell(iOrb, iSp) : orb%posShell(iOrb + 1, iSp) - 1, iSp) =&
            & 0.5_dp * xi(iOrb, iSp) * Lz(1 : 2 * ll + 1, 1 : 2 * ll + 1)
        atomPlus(orb%posShell(iOrb, iSp) : orb%posShell(iOrb + 1, iSp) - 1, &
            & orb%posShell(iOrb, iSp) : orb%posShell(iOrb + 1, iSp) - 1, iSp) =&
            & 0.5_dp * xi(iOrb, iSp) * Lplus(1 : 2 * ll + 1, 1 : 2 * ll + 1)
      end do
    end do

    do iAt = 1, nAtom
      iSp = species(iAt)
      HSqrCplx(iAtomStart(iAt) : iAtomStart(iAt + 1) - 1, &
          & iAtomStart(iAt) : iAtomStart(iAt + 1) - 1) = &
          & HSqrCplx(iAtomStart(iAt) : iAtomStart(iAt + 1) - 1, &
          & iAtomStart(iAt) : iAtomStart(iAt + 1) - 1) &
          & + atomZ(1 : orb%nOrbSpecies(iSp), 1 : orb%nOrbSpecies(iSp), iSp)
      HSqrCplx(nOrb + iAtomStart(iAt) : nOrb + iAtomStart(iAt + 1) - 1, &
          & nOrb + iAtomStart(iAt) : nOrb + iAtomStart(iAt + 1) - 1) = &
          & HSqrCplx(nOrb + iAtomStart(iAt) : nOrb + iAtomStart(iAt + 1) - 1, &
          & nOrb + iAtomStart(iAt) : nOrb + iAtomStart(iAt + 1) - 1) &
          & - atomZ(1 : orb%nOrbSpecies(iSp), 1 : orb%nOrbSpecies(iSp), iSp)
      HSqrCplx(nOrb+ iAtomStart(iAt) : nOrb + iAtomStart(iAt + 1) - 1, &
          & iAtomStart(iAt) : iAtomStart(iAt + 1) - 1) = &
          & HSqrCplx(nOrb + iAtomStart(iAt) : nOrb + iAtomStart(iAt + 1) - 1, &
          & iAtomStart(iAt) : iAtomStart(iAt + 1) - 1) &
          & + atomPlus(1 : orb%nOrbSpecies(iSp), 1 : orb%nOrbSpecies(iSp), iSp)
    end do

  end subroutine addOnsiteSpinOrbitContrib


  !> Calculate Mulliken population from sparse density matrix.
  subroutine getMullikenPopulation(rhoPrim, over, orb, neighborList, nNeighbor, img2CentCell,&
      & iPair, qOrb, iRhoPrim, qBlock, qiBlock)

    !> sparse density matrix
    real(dp), intent(in) :: rhoPrim(:,:)

    !> sparse overlap matrix
    real(dp), intent(in) :: over(:)

    !> Atomic orbital information
    type(TOrbitals), intent(in) :: orb

    !> Atomic neighbours
    type(TNeighborList), intent(in) :: neighborList

    !> Number of neighbours for each atom within overlap distance
    integer, intent(in) :: nNeighbor(:)

    !> image to actual atom indexing
    integer, intent(in) :: img2CentCell(:)

    !> sparse matrix indexing array
    integer, intent(in) :: iPair(:,:)

    !> orbital charges
    real(dp), intent(out) :: qOrb(:,:,:)

    !> imaginary part of density matrix
    real(dp), intent(in), optional :: iRhoPrim(:,:)

    !> Dual atomic charges
    real(dp), intent(out), optional :: qBlock(:,:,:,:)

    !> Imaginary part of dual atomic charges
    real(dp), intent(out), optional :: qiBlock(:,:,:,:)

    integer :: iSpin

    qOrb(:,:,:) = 0.0_dp
    do iSpin = 1, size(qOrb, dim=3)
      call mulliken(qOrb(:,:,iSpin), over, rhoPrim(:,iSpin), orb, neighborList%iNeighbor,&
          & nNeighbor, img2CentCell, iPair)
    end do

    if (present(qBlock)) then
      qBlock(:,:,:,:) = 0.0_dp
      do iSpin = 1, size(qBlock, dim=4)
        call mulliken(qBlock(:,:,:,iSpin), over, rhoPrim(:,iSpin), orb, neighborList%iNeighbor,&
            & nNeighbor, img2CentCell, iPair)
      end do
    end if

    if (present(qiBlock)) then
      qiBlock(:,:,:,:) = 0.0_dp
      do iSpin = 1, size(qiBlock, dim=4)
        call skewMulliken(qiBlock(:,:,:,iSpin), over, iRhoPrim(:,iSpin), orb,&
            & neighborList%iNeighbor, nNeighbor, img2CentCell, iPair)
      end do
    end if

  end subroutine getMullikenPopulation


  !> Calculates various energy contributions
  subroutine getEnergies(qOrb, q0, chargePerShell, species, tEField, tScc, tXlbomd, tDftbU,&
      & tDualSpinOrbit, rhoPrim, H0, orb, neighborList, nNeighbor, img2CentCell, iPair, cellVol,&
      & pressure, TS, potential, energy, thirdOrd, qBlock, qiBlock, nDftbUFunc, UJ, nUJ, iUJ, niUJ,&
      & xi)

    !> Electrons in each atomic orbital
    real(dp), intent(in) :: qOrb(:,:,:)

    !> reference charges
    real(dp), intent(in) :: q0(:,:,:)

    !> electrons in each atomi shell
    real(dp), intent(in) :: chargePerShell(:,:,:)

    !> chemical species
    integer, intent(in) :: species(:)

    !> is an external electric field present
    logical, intent(in) :: tEField

    !> is this a self-consistent calculation
    logical, intent(in) :: tScc

    !> Is the extended Lagrangian being used for MD
    logical, intent(in) :: tXlbomd

    !> Are there orbital potentials present
    logical, intent(in) :: tDftbU

    !> Is dual spin orbit being used
    logical, intent(in) :: tDualSpinOrbit

    !> density matrix in sparse storage
    real(dp), intent(in) :: rhoPRim(:,:)

    !> non-self-consistent hamiltonian
    real(dp), intent(in) :: H0(:)

    !> atomic orbital information
    type(TOrbitals), intent(in) :: orb

    !> neighbour list
    type(TNeighborList), intent(in) :: neighborList

    !> Number of neighbours within cutoff for each atom
    integer, intent(in) :: nNeighbor(:)

    !> image to real atom mapping
    integer, intent(in) :: img2CentCell(:)

    !> index for sparse large matrices
    integer, intent(in) :: iPair(:,:)

    !> unit cell volume
    real(dp), intent(in) :: cellVol

    !> external pressure
    real(dp), intent(in) :: pressure

    !> electron entropy contribution
    real(dp), intent(in) :: TS(:)

    !> potentials acting
    type(TPotentials), intent(in) :: potential

    !> energy contributions
    type(TEnergies), intent(inout) :: energy

    !> 3rd order settings
    type(ThirdOrder), intent(inout), optional :: thirdOrd

    !> block (dual) atomic populations
    real(dp), intent(in), optional :: qBlock(:,:,:,:)

    !> Imaginary part of block atomic populations
    real(dp), intent(in), optional :: qiBlock(:,:,:,:)

    !> which DFTB+U functional (if used)
    integer, intent(in), optional :: nDftbUFunc

    !> U-J prefactors in DFTB+U
    real(dp), intent(in), optional :: UJ(:,:)

    !> Number DFTB+U blocks of shells for each atom type
    integer, intent(in), optional :: nUJ(:)

    !> which shells are in each DFTB+U block
    integer, intent(in), optional :: iUJ(:,:,:)

    !> Number of shells in each DFTB+U block
    integer, intent(in), optional :: niUJ(:,:)

    !> Spin orbit constants
    real(dp), intent(in), optional :: xi(:,:)

    integer :: nSpin

    nSpin = size(rhoPrim, dim=2)

    ! Tr[H0 * Rho] can be done with the same algorithm as Mulliken-analysis
    energy%atomNonSCC(:) = 0.0_dp
    call mulliken(energy%atomNonSCC, rhoPrim(:,1), H0, orb, neighborList%iNeighbor, nNeighbor,&
        & img2CentCell, iPair)
    energy%EnonSCC =  sum(energy%atomNonSCC)

    if (tEfield) then
      energy%atomExt = sum(qOrb(:,:,1) - q0(:,:,1), dim=1) * potential%extAtom(:,1)
      energy%Eext = sum(energy%atomExt)
    end if

    if (tSCC) then
      if (tXlbomd) then
        call getEnergyPerAtom_SCC_Xlbomd(species, orb, qOrb, q0, energy%atomSCC)
      else
        call getEnergyPerAtom_SCC(energy%atomSCC)
      end if
      energy%Escc = sum(energy%atomSCC)
      if (present(thirdOrd)) then
        if (tXlbomd) then
          call thirdOrd%getEnergyPerAtomXlbomd(qOrb, q0, species, orb, energy%atom3rd)
        else
          call thirdOrd%getEnergyPerAtom(energy%atom3rd)
        end if
        energy%e3rd = sum(energy%atom3rd)
      end if

      if (nSpin > 1) then
        energy%atomSpin(:) = 0.5_dp * sum(sum(potential%intShell(:,:,2:nSpin)&
            & * chargePerShell(:,:,2:nSpin), dim=1), dim=2)
        energy%Espin = sum(energy%atomSpin)
      end if

    end if

    if (tDftbU) then
      call E_DFTBU(energy%atomDftbu, qBlock, species, orb, nDFTBUfunc, UJ, nUJ, niUJ, iUJ, qiBlock)
      energy%Edftbu = sum(energy%atomDftbu)
    end if

    if (tDualSpinOrbit) then
      energy%atomLS(:) = 0.0_dp
      call getEnergySpinOrbit(energy%atomLS, qiBlock, xi, orb, species)
      energy%ELS = sum(energy%atomLS)
    end if

    energy%Eelec = energy%EnonSCC + energy%ESCC + energy%Espin + energy%ELS + energy%Edftbu&
        & + energy%Eext + energy%e3rd
    energy%atomElec(:) = energy%atomNonSCC + energy%atomSCC + energy%atomSpin + energy%atomDftbu&
        & + energy%atomLS + energy%atomExt + energy%atom3rd
    energy%atomTotal(:) = energy%atomElec + energy%atomRep + energy%atomDisp
    energy%Etotal = energy%Eelec + energy%Erep + energy%eDisp
    energy%EMermin = energy%Etotal - sum(TS)
    energy%EGibbs = energy%EMermin + cellVol * pressure

  end subroutine getEnergies


  !> Checks for the presence of a stop file on disc.
  function hasStopFile(fileName) result(tStop)

    !> name of file to check for
    character(*), intent(in) :: fileName

    !> Is the file present
    logical :: tStop

    inquire(file=fileName, exist=tStop)
    if (tStop) then
      write(stdOut, "(3A)") "Stop file '" // fileName // "' found."
    end if

  end function hasStopFile


  !> Returns input charges for next SCC iteration.
  subroutine getNextInputCharges(pChrgMixer, qOutput, qOutRed, orb, nIneqOrb, iEqOrbitals,&
      & iGeoStep, iSccIter, minSccIter, maxSccIter, sccTol, tStopScc, tDftbU, tReadChrg, qInput,&
      & qInpRed, sccErrorQ, tConverged, qBlockOut, iEqBlockDftbU, qBlockIn, qiBlockOut,&
      & iEqBlockDftbuLS, species0, nUJ, iUJ, niUJ, qiBlockIn)

    !> Charge mixing object
    type(OMixer), intent(inout) :: pChrgMixer

    !> Output electrons
    real(dp), intent(inout) :: qOutput(:,:,:)

    !> Output electrons reduced by unique orbital types
    real(dp), intent(inout) :: qOutRed(:)

    !> Atomic orbital data
    type(TOrbitals), intent(in) :: orb

    !> Total number of inequivalent atomic orbitals
    integer, intent(in) :: nIneqOrb

    !> Equivalence relations between orbitals
    integer, intent(in) :: iEqOrbitals(:,:,:)

    !> Number of current geometry step
    integer, intent(in) :: iGeoStep

    !> Number of current SCC step
    integer, intent(in) :: iSccIter

    !> minumum number of SCC iterations to perform
    integer, intent(in) :: minSccIter

    !> maximum number of SCC iterations before terminating loop
    integer, intent(in) :: maxSccIter

    !> Tolerance on SCC charges between input and output
    real(dp), intent(in) :: sccTol

    !> Should the SCC loop stop
    logical, intent(in) :: tStopScc

    !> are orbital potentials being used
    logical, intent(in) :: tDftbU

    !> Were intial charges read from disc?
    logical, intent(in) :: tReadChrg

    !> Resulting input charges for next SCC iteration
    real(dp), intent(inout) :: qInput(:,:,:)

    !> Equivalence reduced input charges
    real(dp), intent(inout) :: qInpRed(:)

    !> SCC error
    real(dp), intent(out) :: sccErrorQ

    !> Has the calculation converged>
    logical, intent(out) :: tConverged

    !> Dual output charges
    real(dp), intent(inout), optional :: qBlockOut(:,:,:,:)

    !> equivalence mapping for dual charge blocks
    integer, intent(in), optional :: iEqBlockDftbu(:,:,:,:)

    !> block charge input (if needed for orbital potentials)
    real(dp), intent(out), optional ::qBlockIn(:,:,:,:)

    !> Imaginary part of block charges
    real(dp), intent(in), optional :: qiBlockOut(:,:,:,:)

    !> Equivalence mappings in the case of spin orbit and DFTB+U
    integer, intent(in), optional :: iEqBlockDftbuLS(:,:,:,:)

    !> atomic species for atoms
    integer, intent(in), optional :: species0(:)

    !> Number DFTB+U blocks of shells for each atom type
    integer, intent(in), optional :: nUJ(:)

    !> which shells are in each DFTB+U block
    integer, intent(in), optional :: iUJ(:,:,:)

    !> Number of shells in each DFTB+U block
    integer, intent(in), optional :: niUJ(:,:)

    !> Imaginary part of block atomic input populations
    real(dp), intent(out), optional :: qiBlockIn(:,:,:,:)

    real(dp), allocatable :: qDiffRed(:)
    integer :: nSpin

    nSpin = size(qOutput, dim=3)
    call reduceCharges(orb, nIneqOrb, iEqOrbitals, qOutput, qOutRed, qBlockOut, iEqBlockDftbu,&
        & qiBlockOut, iEqBlockDftbuLS)
    qDiffRed = qOutRed - qInpRed
    sccErrorQ = maxval(abs(qDiffRed))
    tConverged = (sccErrorQ < sccTol)&
        & .and. (iSCCiter >= minSCCIter .or. tReadChrg .or. iGeoStep > 0)
    if ((.not. tConverged) .and. (iSCCiter /= maxSccIter .and. .not. tStopScc)) then
      ! Avoid mixing of spin unpolarised density for spin polarised cases, this is only a problem in
      ! iteration 1, as there is only the (spin unpolarised!) atomic input density at that
      ! point. (Unless charges had been initialized externally)
      if ((iSCCIter + iGeoStep) == 1 .and. (nSpin > 1 .or. tDFTBU) .and. .not. tReadChrg) then
        qInpRed(:) = qOutRed
        qInput(:,:,:) = qOutput
        if (present(qBlockIn)) then
          qBlockIn(:,:,:,:) = qBlockOut
          if (present(qiBlockIn)) then
            qiBlockIn(:,:,:,:) = qiBlockOut
          end if
        end if
      else
        call mix(pChrgMixer, qInpRed, qDiffRed)
        call expandCharges(qInpRed, orb, nIneqOrb, iEqOrbitals, qInput, qBlockIn, iEqBlockDftbu,&
            & species0, nUJ, iUJ, niUJ, qiBlockIn, iEqBlockDftbuLS)
      end if
    end if

  end subroutine getNextInputCharges


  !> Reduce charges according to orbital equivalency rules.
  subroutine reduceCharges(orb, nIneqOrb, iEqOrbitals, qOrb, qRed, qBlock, iEqBlockDftbu,&
      & qiBlock, iEqBlockDftbuLS)

    !> Atomic orbital information
    type(TOrbitals), intent(in) :: orb

    !> Number of unique types of atomic orbitals
    integer, intent(in) :: nIneqOrb

    !> equivalence index
    integer, intent(in) :: iEqOrbitals(:,:,:)

    !> Electrons in atomic orbitals
    real(dp), intent(inout) :: qOrb(:,:,:)

    !> Reduction of atomic populations
    real(dp), intent(out) :: qRed(:)

    !> Block (dual) populations, if also being reduced
    real(dp), intent(inout), optional :: qBlock(:,:,:,:)

    !> equivalences for block charges
    integer, intent(in), optional :: iEqBlockDftbu(:,:,:,:)

    !> Imaginary part of block charges if present
    real(dp), intent(in), optional :: qiBlock(:,:,:,:)

    !> Equivalences for spin orbit if needed
    integer, intent(in), optional :: iEqBlockDftbuLS(:,:,:,:)

    qRed(:) = 0.0_dp
    call qm2ud(qOrb)
    if (present(qBlock)) then
      call qm2ud(qBlock)
    end if
    call OrbitalEquiv_reduce(qOrb, iEqOrbitals, orb, qRed(1:nIneqOrb))
    if (present(qBlock)) then
      call AppendBlock_reduce(qBlock, iEqBlockDFTBU, orb, qRed)
      if (present(qiBlock)) then
        call AppendBlock_reduce(qiBlock, iEqBlockDFTBULS, orb, qRed, skew=.true.)
      end if
    end if
    call ud2qm(qOrb)
    if (present(qBlock)) then
      call ud2qm(qBlock)
    end if

  end subroutine reduceCharges


  !> Expand reduced charges according orbital equivalency rules.
  subroutine expandCharges(qRed, orb, nIneqOrb, iEqOrbitals, qOrb, qBlock, iEqBlockDftbu, species0,&
      & nUJ, iUJ, niUJ, qiBlock, iEqBlockDftbuLS)

    !> Reduction of atomic populations
    real(dp), intent(in) :: qRed(:)

    !> Atomic orbital information
    type(TOrbitals), intent(in) :: orb

    !> Number of unique types of atomic orbitals
    integer, intent(in) :: nIneqOrb

    !> equivalence index
    integer, intent(in) :: iEqOrbitals(:,:,:)

    !> Electrons in atomic orbitals
    real(dp), intent(out) :: qOrb(:,:,:)

    !> Block (dual) populations, if also stored in reduced form
    real(dp), intent(inout), optional :: qBlock(:,:,:,:)

    !> equivalences for block charges
    integer, intent(in), optional :: iEqBlockDftbU(:,:,:,:)

    !> species of central cell atoms
    integer, intent(in), optional :: species0(:)

    !> Number DFTB+U blocks of shells for each atom type
    integer, intent(in), optional :: nUJ(:)

    !> which shells are in each DFTB+U block
    integer, intent(in), optional :: iUJ(:,:,:)

    !> Number of shells in each DFTB+U block
    integer, intent(in), optional :: niUJ(:,:)

    !> Imaginary part of block atomic populations
    real(dp), intent(inout), optional :: qiBlock(:,:,:,:)

    !> Equivalences for spin orbit if needed
    integer, intent(in), optional :: iEqBlockDftbULS(:,:,:,:)

    integer :: nSpin

    @:ASSERT(present(qBlock) .eqv. present(iEqBlockDftbU))
    @:ASSERT(present(qBlock) .eqv. present(species0))
    @:ASSERT(present(qBlock) .eqv. present(nUJ))
    @:ASSERT(present(qBlock) .eqv. present(iUJ))
    @:ASSERT(present(qBlock) .eqv. present(niUJ))
    @:ASSERT(.not. present(qiBlock) .or. present(qBlock))
    @:ASSERT(present(qiBlock) .eqv. present(iEqBlockDftbuLS))

    nSpin = size(qOrb, dim=3)
    call OrbitalEquiv_expand(qRed(1:nIneqOrb), iEqOrbitals, orb, qOrb)
    if (present(qBlock)) then
      qBlock(:,:,:,:) = 0.0_dp
      call Block_expand(qRed, iEqBlockDftbu, orb, qBlock, species0, nUJ, niUJ, iUJ,&
          & orbEquiv=iEqOrbitals)
      if (present(qiBlock)) then
        call Block_expand(qRed, iEqBlockDftbuLS, orb, qiBlock, species0, nUJ, niUJ, iUJ,&
            & skew=.true.)
      end if
    end if
    if (nSpin == 2) then
      call ud2qm(qOrb)
      if (present(qBlock)) then
        call ud2qm(qBlock)
      end if
    end if

  end subroutine expandCharges


  !> Get some info about scc convergence.
  subroutine getSccInfo(iSccIter, Eelec, EelecOld, diffElec)

    !> Iteration number
    integer, intent(in) :: iSccIter

    !> Electronic energy
    real(dp), intent(in) :: Eelec

    !> old electronic energy, overwritten on exit with current value
    real(dp), intent(inout) :: EelecOld

    !> difference in electronic energies between iterations
    real(dp), intent(out) :: diffElec

    if (iScciter > 1) then
      diffElec = Eelec - EelecOld
    else
      diffElec = 0.0_dp
    end if
    EelecOld = Eelec

  end subroutine getSccInfo


  !> Prints info about scc convergence.
  subroutine printSccInfo(tDftbU, iSccIter, Eelec, diffElec, sccErrorQ)

    !> Are orbital potentials being used
    logical, intent(in) :: tDftbU

    !> Iteration count
    integer, intent(in) :: iSccIter

    !> electronic energy
    real(dp), intent(in) :: Eelec

    !> Difference in electronic energy between this iteration and the last
    real(dp), intent(in) :: diffElec

    !> Maximum charge difference between input and output
    real(dp), intent(in) :: sccErrorQ

    if (tDFTBU) then
      write(stdOut, "(I5,E18.8,E18.8,E18.8)") iSCCIter, Eelec, diffElec, sccErrorQ
    else
      write(stdOut, "(I5,E18.8,E18.8,E18.8)") iSCCIter, Eelec, diffElec, sccErrorQ
    end if

  end subroutine printSccInfo


  !> Whether restart information needs to be written in the current scc loop.
  function needsSccRestartWriting(restartFreq, iGeoStep, iSccIter, minSccIter, maxSccIter, tMd,&
      & tGeoOpt, tDerivs, tConverged, tReadChrg, tStopScc) result(tRestart)

    !> frequency of charge  write out
    integer, intent(in) :: restartFreq

    !> current geometry step
    integer, intent(in) :: iGeoStep

    !> current SCC step
    integer, intent(in) :: iSccIter

    !> minimum number of SCC cycles to perform
    integer, intent(in) :: minSccIter

    !> maximum number of SCC cycles to perform
    integer, intent(in) :: maxSccIter

    !> is this molecular dynamics
    logical, intent(in) :: tMd

    !> Is there geometry optimisation
    logical, intent(in) :: tGeoOpt

    !> are finite difference changes happening
    logical, intent(in) :: tDerivs

    !> Is this converged SCC
    logical, intent(in) :: tConverged

    !> have the charges been read from disc
    logical, intent(in) :: tReadChrg

    !> Has the SCC cycle been stopped?
    logical, intent(in) :: tStopScc

    !> resulting decision as to whether to write charges to disc
    logical :: tRestart

    tRestart = ((restartFreq > 0 .and. .not. (tMD .or. tGeoOpt .or. tDerivs) .and. maxSccIter > 1)&
        & .and. (tConverged&
        & .or. ((iSCCIter >= minSCCIter .or. tReadChrg .or. iGeoStep > 0)&
        & .and. ((iSCCIter == maxSccIter .or. tStopScc)&
        & .or. mod(iSCCIter, restartFreq) == 0))))

  end function needsSccRestartWriting


  !> Write restart information.
  subroutine writeRestart(fChargeIn, orb, qInput, qBlockIn, qiBlockIn)

    !> File name for charge output
    character(*), intent(in) :: fChargeIn

    !> Atomic orbital information
    type(TOrbitals), intent(in) :: orb

    !> Input charges
    real(dp), intent(in) :: qInput(:,:,:)

    !> input block charges
    real(dp), intent(in), optional :: qBlockIn(:,:,:,:)

    !> imaginary part of input block charges
    real(dp), intent(in), optional :: qiBlockIn(:,:,:,:)

    if (present(qBlockIn)) then
      if (present(qiBlockIn)) then
        call writeQToFile(qInput, fChargeIn, orb, qBlockIn, qiBlockIn)
      else
        call writeQToFile(qInput, fChargeIn, orb, qBlockIn)
      end if
    else
      call writeQToFile(qInput, fChargeIn, orb)
    end if
    write(stdout, "(2A)") ">> Charges saved for restart in ", fChargeIn

  end subroutine writeRestart


  !> Stop if linear response module can not be invoked due to unimplemented combinations of
  !> features.
  subroutine ensureLinRespConditions(t3rd, tRealHS, tPeriodic, tForces)

    !> 3rd order hamiltonian contributions included
    logical, intent(in) :: t3rd

    !> a real hamiltonian
    logical, intent(in) :: tRealHs

    !> periodic boundary conditions
    logical, intent(in) :: tPeriodic

    !> forces being evaluated
    logical, intent(in) :: tForces

    if (t3rd) then
      call error("Third order currently incompatible with excited state")
    end if
    if (.not. tRealHS) then
      call error("Only real systems are supported for excited state calculations")
    end if
    if (tPeriodic .and. tForces) then
      call error("Forces in the excited state for periodic geometries are currently&
          & unavailable")
    end if

  end subroutine ensureLinRespConditions


  !> Do the linear response excitation calculation.
  subroutine calculateLinRespExcitations(lresp, qOutput, q0, over, HSqrReal, eigen, filling,&
      & coord0, species, species0, speciesName, orb, skHamCont, skOverCont, fdAutotest, fdEigvec,&
      & runId, neighborList, nNeighbor, iAtomStart, iPair, img2CentCell, tWriteAutotest, tForces,&
      & tLinRespZVect, tPrintExcitedEigvecs, nonSccDeriv, energy, SSqrReal, rhoSqrReal,&
      & excitedDerivs, occNatural)
    type(LinResp), intent(inout) :: lresp
    real(dp), intent(in) :: qOutput(:,:,:)
    real(dp), intent(in) :: q0(:,:,:)
    real(dp), intent(in) :: over(:)
    real(dp), intent(in) :: HSqrReal(:,:,:)
    real(dp), intent(in) :: eigen(:,:)
    real(dp), intent(in) :: filling(:,:)
    real(dp), intent(in) :: coord0(:,:)
    integer, intent(in) :: species(:)
    integer, intent(in) :: species0(:)
    character(*), intent(in) :: speciesName(:)
    type(TOrbitals), intent(in) :: orb
    type(OSlakoCont), intent(in) :: skHamCont
    type(OSlakoCont), intent(in) :: skOverCont
    integer, intent(in) :: fdAutotest
    integer, intent(in) :: fdEigvec
    integer, intent(in) :: runId
    type(TNeighborList), intent(in) :: neighborList
    integer, intent(in) :: nNeighbor(:)
    integer, intent(in) :: iAtomStart(:)
    integer, intent(in) :: iPair(:,:)
    integer, intent(in) :: img2CentCell(:)
    logical, intent(in) :: tWriteAutotest
    logical, intent(in) :: tForces
    logical, intent(in) :: tLinRespZVect
    logical, intent(in) :: tPrintExcitedEigvecs
    type(NonSccDiff), intent(in) :: nonSccDeriv
    type(TEnergies), intent(inout) :: energy
    real(dp), intent(out) :: SSqrReal(:,:)
    real(dp), intent(inout), optional :: rhoSqrReal(:,:,:)
    real(dp), intent(out), optional :: excitedDerivs(:,:)
    real(dp), intent(out), optional :: occNatural(:)

    real(dp), allocatable :: dQAtom(:)
    real(dp), allocatable, target :: naturalOrbs(:)
    real(dp), pointer :: pNaturalOrbs2(:,:), pNaturalOrbs3(:,:,:)
    integer :: iSpin, nSpin, nAtom
    logical :: tSpin

    nAtom = size(qOutput, dim=2)
    nSpin = size(eigen, dim=2)
    tSpin = (nSpin == 2)

    energy%Eexcited = 0.0_dp
    allocate(dQAtom(nAtom))
    dQAtom(:) = sum(qOutput(:,:,1) - q0(:,:,1), dim=1)
    call unpackHS(SSqrReal, over, neighborList%iNeighbor, nNeighbor, iAtomStart, iPair,&
        & img2CentCell)
    call blockSymmetrizeHS(SSqrReal, iAtomStart)
    if (tForces) then
      do iSpin = 1, nSpin
        call blockSymmetrizeHS(rhoSqrReal(:,:,iSpin), iAtomStart)
      end do
    end if
    if (tWriteAutotest) then
      open(fdAutotest, file=autotestTag, position="append")
    end if

    if (tLinRespZVect) then
      if (tPrintExcitedEigVecs) then
        allocate(naturalOrbs(orb%nOrb * orb%nOrb))
        pNaturalOrbs2(1:orb%nOrb, 1:orb%nOrb) => naturalOrbs
        pNaturalOrbs3(1:orb%nOrb, 1:orb%nOrb, 1:1) => naturalOrbs
      else
        pNaturalOrbs2 => null()
        pNaturalOrbs3 => null()
      end if
      call addGradients(tSpin, lresp, iAtomStart, HSqrReal, eigen, SSqrReal,&
          & filling, coord0, dQAtom, species0, neighborList%iNeighbor, img2CentCell, orb,&
          & skHamCont, skOverCont, tWriteAutotest, fdAutotest, energy%Eexcited, excitedDerivs, &
          & nonSccDeriv, rhoSqrReal, occNatural=occNatural, naturalOrbs=pNaturalOrbs2)
      if (tPrintExcitedEigvecs) then
        call writeEigvecs(fdEigvec, runId, nAtom, nSpin, neighborList, nNeighbor, iAtomStart,&
            & iPair, img2CentCell, orb, species, speciesName, over, pNaturalOrbs3, SSqrReal,&
            & fileName="excitedOrbs")
      end if
    else
      call calcExcitations(tSpin, lresp, iAtomStart, HSqrReal, eigen, SSqrReal, filling, coord0,&
          & dQAtom, species0, neighborList%iNeighbor, img2CentCell, orb, tWriteAutotest,&
          & fdAutotest, energy%Eexcited)
    end if
    energy%Etotal = energy%Etotal + energy%Eexcited
    energy%EMermin = energy%EMermin + energy%Eexcited
    energy%EGibbs = energy%EGibbs + energy%Eexcited
    if (tWriteAutotest) then
      close(fdAutotest)
    end if

  end subroutine calculateLinRespExcitations


  !> Get the XLBOMD charges for the current geometry.
  subroutine getXlbomdCharges(xlbomdIntegrator, qOutRed, pChrgMixer, orb, nIneqOrb, iEqOrbitals,&
      & qInput, qInpRed, iEqBlockDftbu, qBlockIn, species0, nUJ, iUJ, niUJ, iEqBlockDftbuLS,&
      & qiBlockIn)
    type(Xlbomd), intent(inout) :: xlbomdIntegrator
    real(dp), intent(in) :: qOutRed(:)
    type(OMixer), intent(inout) :: pChrgMixer
    type(TOrbitals), intent(in) :: orb
    integer, intent(in) :: nIneqOrb
    integer, intent(in) :: iEqOrbitals(:,:,:)
    real(dp), intent(out) :: qInput(:,:,:)
    real(dp), intent(out) :: qInpRed(:)
    integer, intent(in), optional :: iEqBlockDftbU(:,:,:,:)
    integer, intent(in), optional :: species0(:)
    real(dp), intent(out), optional :: qBlockIn(:,:,:,:)
    integer, intent(in), optional :: nUJ(:)
    integer, intent(in), optional :: iUJ(:,:,:)
    integer, intent(in), optional :: niUJ(:,:)
    integer, intent(in), optional :: iEqBlockDftbuLS(:,:,:,:)
    real(dp), intent(out), optional :: qiBlockIn(:,:,:,:)

    real(dp), allocatable :: invJacobian(:,:)

    if (xlbomdIntegrator%needsInverseJacobian()) then
      write(stdOut, "(A)") ">> Updating XLBOMD Inverse Jacobian"
      allocate(invJacobian(nIneqOrb, nIneqOrb))
      call getInverseJacobian(pChrgMixer, invJacobian)
      call xlbomdIntegrator%setInverseJacobian(invJacobian)
      deallocate(invJacobian)
    end if
    call xlbomdIntegrator%getNextCharges(qOutRed(1:nIneqOrb), qInpRed(1:nIneqOrb))
    call expandCharges(qInpRed, orb, nIneqOrb, iEqOrbitals, qInput, qBlockIn, iEqBlockDftbu,&
        & species0, nUJ, iUJ, niUJ, qiBlockIn, iEqBlockDftbuLS)

  end subroutine getXlbomdCharges


  !> Writes the eigenvectors to disc.
  subroutine writeEigenvectors(nSpin, fdEigvec, runId, nAtom, neighborList, nNeighbor, cellVec,&
      & iCellVec, iAtomStart, iPair, img2CentCell, species, speciesName, orb, kPoint, over,&
      & HSqrReal, SSqrReal, HSqrCplx, SSqrCplx, storeEigvecsReal, storeEigvecsCplx)
    integer, intent(in) :: nSpin
    integer, intent(in) :: fdEigvec
    integer, intent(in) :: runId
    integer, intent(in) :: nAtom
    type(TNeighborList), intent(in) :: neighborList
    integer, intent(in) :: nNeighbor(:)
    integer, intent(in) :: iCellVec(:)
    integer, intent(in) :: img2CentCell(:)
    integer, intent(in) :: iAtomStart(:)
    integer, intent(in) :: iPair(:,:)
    real(dp), intent(in) :: cellVec(:,:)
    integer, intent(in) :: species(:)
    character(*), intent(in) :: speciesName(:)
    type(TOrbitals), intent(in) :: orb
    real(dp), intent(in) :: kPoint(:,:)
    real(dp), intent(in) :: over(:)
    real(dp), intent(inout), optional :: HSqrReal(:,:,:)
    real(dp), intent(inout), optional :: SSqrReal(:,:)
    complex(dp), intent(inout), optional :: HSqrCplx(:,:,:,:)
    complex(dp), intent(inout), optional :: SSqrCplx(:,:)
    type(OFifoRealR2), intent(inout), optional :: storeEigvecsReal(:)
    type(OFifoCplxR2), intent(inout), optional :: storeEigvecsCplx(:)

    @:ASSERT(present(HSqrReal) .neqv. present(HSqrReal))
    @:ASSERT(present(SSqrReal) .neqv. present(SSqrReal))
    @:ASSERT(.not. present(storeEigvecsReal) .or. present(HSqrReal))
    @:ASSERT(.not. present(storeEigvecsCplx) .or. present(HSqrCplx))

    if (present(HSqrCplx)) then
      !> Complex Pauli-Hamiltonian without k-points
      call writeEigvecs(fdEigvec, runId, nAtom, nSpin, neighborList, nNeighbor, cellVec, iCellVec,&
          & iAtomStart, iPair, img2CentCell, orb, species, speciesName, over, kPoint, HSqrCplx,&
          & SSqrCplx, storeEigvecsCplx)
    else
      !> Real Hamiltonian
      call writeEigvecs(fdEigvec, runId, nAtom, nSpin, neighborList, nNeighbor, iAtomStart, iPair,&
          & img2CentCell, orb, species, speciesName, over, HSqrReal, SSqrReal, storeEigvecsReal)
    end if

  end subroutine writeEigenvectors


  !> Write projected eigenvectors.
  subroutine writeProjectedEigenvectors(regionLabels, fdProjEig, eigen, nSpin, neighborList,&
      & nNeighbor, cellVec, iCellVec, iAtomStart, iPair, img2CentCell, orb, over, kPoint, kWeight,&
      & iOrbRegion, HSqrReal, SSqrReal, HSqrCplx, SSqrCplx, storeEigvecsReal, storeEigvecsCplx)
    type(ListCharLc), intent(inout) :: regionLabels
    integer, intent(in) :: fdProjEig(:)
    real(dp), intent(in) :: eigen(:,:,:)
    integer, intent(in) :: nSpin
    type(TNeighborList), intent(in) :: neighborList
    integer, intent(in) :: nNeighbor(:)
    real(dp), intent(in) :: cellVec(:,:)
    integer, intent(in) :: iCellVec(:)
    integer, intent(in) :: iAtomStart(:)
    integer, intent(in) :: iPair(:,:)
    integer, intent(in) :: img2CentCell(:)
    type(TOrbitals), intent(in) :: orb
    real(dp), intent(in) :: over(:)
    real(dp), intent(in) :: kPoint(:,:)
    real(dp), intent(in) :: kWeight(:)
    type(ListIntR1), intent(inout) :: iOrbRegion
    real(dp), intent(inout), optional :: HSqrReal(:,:,:)
    real(dp), intent(inout), optional :: SSqrReal(:,:)
    complex(dp), intent(inout), optional :: HSqrCplx(:,:,:,:)
    complex(dp), intent(inout), optional :: SSqrCplx(:,:)
    type(OFifoRealR2), intent(inout), optional :: storeEigvecsReal(:)
    type(OFifoCplxR2), intent(inout), optional :: storeEigvecsCplx(:)

    @:ASSERT(present(HSqrReal) .neqv. present(HSqrReal))
    @:ASSERT(present(SSqrReal) .neqv. present(SSqrReal))
    @:ASSERT(.not. present(storeEigvecsReal) .or. present(HSqrReal))
    @:ASSERT(.not. present(storeEigvecsCplx) .or. present(HSqrCplx))


    if (present(SSqrCplx)) then
      call writeProjEigvecs(regionLabels, fdProjEig, eigen, nSpin, neighborList, nNeighbor,&
          & cellVec, iCellVec, iAtomStart, iPair, img2CentCell, orb, over, kpoint, kWeight,&
          & HSqrCplx, SSqrCplx, iOrbRegion, storeEigvecsCplx)
    else
      call writeProjEigvecs(regionLabels, fdProjEig, eigen, nSpin, neighborList, nNeighbor,&
          & iAtomStart, iPair, img2CentCell, orb, over, HSqrReal, SSqrReal, iOrbRegion,&
          & storeEigvecsReal)
    end if

  end subroutine writeProjectedEigenvectors


  !> Write current geometry to disc
  subroutine writeCurrentGeometry(geoOutFile, pCoord0Out, tLatOpt, tMd, tAppendGeo, tFracCoord,&
      & tPeriodic, tPrintMulliken, species0, speciesName, latVec, iGeoStep, iLatGeoStep, nSpin,&
      & qOutput, velocities)
    character(*), intent(in) :: geoOutFile
    real(dp), intent(in) :: pCoord0Out(:,:)
    logical, intent(in) :: tLatOpt
    logical, intent(in) :: tMd
    logical, intent(in) :: tAppendGeo
    logical, intent(in) :: tFracCoord
    logical, intent(in) :: tPeriodic
    logical, intent(in) :: tPrintMulliken
    integer, intent(in) :: species0(:)
    character(*), intent(in) :: speciesName(:)
    real(dp), intent(in) :: latVec(:,:)
    integer, intent(in) :: iGeoStep
    integer, intent(in) :: iLatGeoStep
    integer, intent(in) :: nSpin
    real(dp), intent(in), optional :: qOutput(:,:,:)
    real(dp), intent(in), optional :: velocities(:,:)

    real(dp), allocatable :: tmpMatrix(:,:)
    integer :: nAtom
    integer :: ii, jj
    character(lc) :: comment, fname

    nAtom = size(pCoord0Out, dim=2)

    fname = trim(geoOutFile) // ".gen"
    if (tPeriodic) then
      call writeGenFormat(fname, pCoord0Out, species0, speciesName, latVec, tFracCoord)
    else
      call writeGenFormat(fname, pCoord0Out, species0, speciesName)
    end if

    fname = trim(geoOutFile) // ".xyz"
    if (tLatOpt) then
      write (comment, "(A, I0, A, I0)") '** Geometry step: ', iGeoStep, ', Lattice step: ',&
          & iLatGeoStep
    elseif (tMD) then
      write(comment, "(A, I0)") 'MD iter: ', iGeoStep
    else
      write(comment,"(A, I0)") 'Geometry Step: ', iGeoStep
    end if

    if (tPrintMulliken) then
      ! For non-colinear spin without velocities write magnetisation into the velocity field
      if (nSpin == 4 .and. .not. present(velocities)) then
        allocate(tmpMatrix(3, nAtom))
        do jj = 1, nAtom
          do ii = 1, 3
            tmpMatrix(ii,jj) = sum(qOutput(:, jj, ii + 1))
          end do
        end do
        ! convert by the inverse of the scaling used in writeXYZFormat
        tmpMatrix(:,:) = tmpMatrix * au__fs / (1000_dp * Bohr__AA)
        call writeXYZFormat(fname, pCoord0Out, species0, speciesName,&
            & charges=sum(qOutput(:,:,1), dim=1), velocities=tmpMatrix, comment=comment,&
            & append=tAppendGeo)
      else
        call writeXYZFormat(fname, pCoord0Out, species0, speciesName,&
            & charges=sum(qOutput(:,:,1),dim=1), velocities=velocities, comment=comment,&
            & append=tAppendGeo)
      end if
    else
      call writeXYZFormat(fname, pCoord0Out, species0, speciesName, velocities=velocities,&
          & comment=comment, append=tAppendGeo)
    end if

  end subroutine writeCurrentGeometry


  !> Prints current total energies
  subroutine printEnergies(energy)

    !> energy components
    type(TEnergies), intent(in) :: energy

    write(stdOut, *)
    write(stdOut, format2U) "Total Energy", energy%Etotal,"H", Hartree__eV * energy%Etotal,"eV"
    write(stdOut, format2U) "Total Mermin free energy", energy%EMermin, "H",&
        & Hartree__eV * energy%EMermin," eV"

  end subroutine printEnergies


  !> Calculates dipole moment.
  subroutine getDipoleMoment(qOutput, q0, coord, dipoleMoment)
    real(dp), intent(in) :: qOutput(:,:,:)
    real(dp), intent(in) :: q0(:,:,:)
    real(dp), intent(in) :: coord(:,:)
    real(dp), intent(out) :: dipoleMoment(:)

    integer :: nAtom, iAtom

    nAtom = size(qOutput, dim=2)
    dipoleMoment(:) = 0.0_dp
    do iAtom = 1, nAtom
      dipoleMoment(:) = dipoleMoment(:) &
          & + sum(q0(:, iAtom, 1) - qOutput(:, iAtom, 1)) * coord(:,iAtom)
    end do

  end subroutine getDipoleMoment


  !> Prints dipole moment calcululated by the derivative of H with respect of the external field.
  subroutine checkDipoleViaHellmannFeynman(sparseSize, rhoPrim, q0, coord0, over, orb,&
      & neighborList, nNeighbor, species, iPair, img2CentCell)
    integer, intent(in) :: sparseSize
    real(dp), intent(in) :: rhoPrim(:,:), q0(:,:,:), coord0(:,:), over(:)
    type(TOrbitals), intent(in) :: orb
    type(TNeighborList), intent(in) :: neighborList
    integer, intent(in) :: nNeighbor(:), species(:), iPair(:,:), img2CentCell(:)

    real(dp), allocatable :: hprime(:,:), dipole(:,:), potentialDerivative(:,:)
    integer :: nAtom
    integer :: iAt, ii

    nAtom = size(q0, dim=2)
    allocate(hprime(sparseSize, 1))
    allocate(dipole(size(q0, dim=1), nAtom))
    allocate(potentialDerivative(nAtom, 1))
    write(stdOut, "(A)", advance='no') 'Hellmann Feynman dipole:'

    ! loop over directions
    do ii = 1, 3
      potentialDerivative(:,:) = 0.0_dp
      ! Potential from dH/dE
      potentialDerivative(:,1) = -coord0(ii,:)
      hprime(:,:) = 0.0_dp
      dipole(:,:) = 0.0_dp
      call add_shift(hprime, over, nNeighbor, neighborList%iNeighbor, species, orb, iPair, nAtom,&
          & img2CentCell, potentialDerivative)

      ! evaluate <psi| dH/dE | psi>
      call mulliken(dipole, hprime(:,1), rhoPrim(:,1), orb, neighborList%iNeighbor, nNeighbor,&
          & img2CentCell, iPair)

      ! add nuclei term for derivative wrt E
      do iAt = 1, nAtom
        dipole(1, iAt) = dipole(1, iAt) + sum(q0(:, iAt, 1)) * coord0(ii, iAt)
      end do
      write(stdOut, "(F12.8)", advance='no') sum(dipole)
    end do
    write(stdOut, *) " au"

  end subroutine checkDipoleViaHellmannFeynman


  !> Calculate the energy weighted density matrix
  !>
  !> NOTE: Dense eigenvector and overlap matrices are overwritten.
  !>
  subroutine getEnergyWeightedDensityMtx(forceType, filling, eigen, kPoint, kWeight, neighborList,&
      & nNeighbor, orb, iAtomStart, iPair, img2CentCell, iCellVEc, cellVec, tRealHS, ham, over,&
      & solver, ERhoPrim, HSqrReal, SSqrReal, HSqrCplx, SSqrCplx, storeEigvecsReal,&
      & storeEigvecsCplx)
    integer, intent(in) :: forceType
    real(dp), intent(in) :: filling(:,:,:), eigen(:,:,:), kPoint(:,:), kWeight(:)
    type(TNeighborList), intent(in) :: neighborList
    integer, intent(in) :: nNeighbor(:)
    type(TOrbitals), intent(in) :: orb
    integer, intent(in) :: iAtomStart(:), iPair(:,:), img2CentCell(:), iCellVec(:)
    real(dp), intent(in) :: cellVec(:,:)
    logical, intent(in) :: tRealHS
    real(dp), intent(in) :: ham(:,:), over(:)
    integer, intent(in) :: solver
    real(dp), intent(out) :: ERhoPrim(:)
    real(dp), intent(inout), optional :: HSqrReal(:,:,:), SSqrReal(:,:)
    complex(dp), intent(inout), optional :: HSqrCplx(:,:,:,:), SSqrCplx(:,:)
    type(OFifoRealR2), intent(inout), optional :: storeEigvecsReal(:)
    type(OFifoCplxR2), intent(inout), optional :: storeEigvecsCplx(:)

    integer :: nSpin

    nSpin = size(ham, dim=2)

    if (nSpin == 4) then
      call getEDensityMtxFromPauliEigvecs(filling, eigen, kPoint, kWeight, neighborList, nNeighbor,&
          & orb, iAtomStart, iPair, img2CentCell, iCellVec, cellVec, tRealHS, HSqrCplx, SSqrCplx,&
          & ERhoPrim, storeEigvecsCplx)
    else if (tRealHS) then
      call getEDensityMtxFromRealEigvecs(forceType, filling, eigen, neighborList, nNeighbor, orb,&
          & iAtomStart, iPair, img2CentCell, ham, over, solver, HSqrReal, SSqrReal, ERhoPrim,&
          & storeEigvecsReal)
    else
      call getEDensityMtxFromComplexEigvecs(forceType, filling, eigen, kPoint, kWeight,&
          & neighborList, nNeighbor, orb, iAtomStart, iPair, img2CentCell, iCellVEc, cellVec, ham,&
          & over, HSqrCplx, SSqrCplx, ERhoPrim, storeEigvecsCplx)
    end if

  end subroutine getEnergyWeightedDensityMtx


  !> Calculates density matrix from real eigenvectors.
  subroutine getEDensityMtxFromRealEigvecs(forceType, filling, eigen, neighborList, nNeighbor, orb,&
      & iAtomStart, iPair, img2CentCell, ham, over, solver, HSqrReal, SSqrReal, ERhoPrim,&
      & storeEigvecsReal)
    integer, intent(in) :: forceType
    real(dp), intent(in) :: filling(:,:,:), eigen(:,:,:)
    type(TNeighborList), intent(in) :: neighborList
    integer, intent(in) :: nNeighbor(:)
    type(TOrbitals), intent(in) :: orb
    integer, intent(in) :: iAtomStart(:), iPair(:,:), img2CentCell(:)
    real(dp), intent(in) :: ham(:,:), over(:)
    integer, intent(in) :: solver
    real(dp), intent(inout) :: HSqrReal(:,:,:), SSqrReal(:,:)
    real(dp), intent(out) :: ERhoPrim(:)
    type(OFifoRealR2), intent(inout), optional :: storeEigvecsReal(:)

    real(dp), allocatable :: HSqrReal2(:,:), eigen2(:,:,:)
    integer :: nSpin, sqrHamSize
    integer :: iS, iS2

    nSpin = size(eigen, dim=3)
    sqrHamSize = size(HSqrReal, dim=1)
    if (forceType == 1 .or. forceType == 2 .or. forceType == 3) then
      allocate(HSqrReal2(sqrHamSize, sqrHamSize))
      if (forceType == 2) then
        allocate(eigen2(sqrHamSize, 1, nSpin))
      end if
    end if

    ERhoPrim(:) = 0.0_dp
    do iS = 1, nSpin
      if (present(storeEigvecsReal)) then
        call get(storeEigvecsReal(iS), HSqrReal(:,:,1))
        iS2 = 1
      else
        iS2 = iS
      end if

      select case (forceType)

      case(0)
        ! Original (nonconsistent) scheme
        if (tDensON2) then
          call makeDensityMatrix(SSqrReal, HSqrReal(:,:,iS2), filling(:,1,iS), eigen(:,1,iS),&
              & neighborlist%iNeighbor, nNeighbor, orb, iAtomStart, img2CentCell)
        else
          call makeDensityMatrix(SSqrReal, HSqrReal(:,:,iS2), filling(:,1,iS), eigen(:,1,iS))
        end if

      case(1)
        ! Recreate eigenvalues for a consistent energy weighted density matrix
        ! (yields, however, incorrect forces for XLBOMD)
        call diagonalize(HSqrReal2, SSqrReal, eigen2(:,1,iS), ham(:,iS), over, &
            & neighborList%iNeighbor, nNeighbor, iAtomStart, iPair, img2CentCell, solver, 'N')
        call makeDensityMatrix(SSqrReal, HSqrReal(:,:,iS2), filling(:,1,iS), eigen2(:,1,iS))

      case(2)
        ! Correct force for XLBOMD for T=0K (DHD)
        call unpackHS(SSqrReal, ham(:,iS), neighborlist%iNeighbor, nNeighbor, iAtomStart,&
            & iPair, img2CentCell)
        call blockSymmetrizeHS(SSqrReal, iAtomStart)
        call makeDensityMatrix(HSqrReal2, HSqrReal(:,:,iS2), filling(:,1,iS))
        ! D H
        call symm(HSqrReal(:,:,iS2), "L", HSqrReal2, SSqrReal)
        ! (D H) D
        call symm(SSqrReal, "R", HSqrReal2, HSqrReal(:,:,iS2), alpha=0.5_dp)

      case(3)
        ! Correct force for XLBOMD for T <> 0K (DHS^-1 + S^-1HD)
        call makeDensityMatrix(SSqrReal, HSqrReal(:,:,iS2), filling(:,1,iS))
        call unpackHS(HSqrReal2, ham(:,iS), neighborlist%iNeighbor, nNeighbor, iAtomStart, iPair,&
            & img2CentCell)
        call blocksymmetrizeHS(HSqrReal2, iAtomStart)
        call symm(HSqrReal(:,:,iS2), "L", SSqrReal, HSqrReal2)
        call unpackHS(SSqrReal, over, neighborlist%iNeighbor, nNeighbor, iAtomStart, iPair,&
            & img2CentCell)
        call symmatinv(SSqrReal)
        call symm(HSqrReal2, "R", SSqrReal, HSqrReal(:,:,iS2), alpha=0.5_dp)
        SSqrReal(:,:) = HSqrReal2 + transpose(HSqrReal2)
      end select
      call packHS(ERhoPrim, SSqrReal, neighborList%iNeighbor, nNeighbor, orb%mOrb, iAtomStart,&
          & iPair, img2CentCell)
    end do

  end subroutine getEDensityMtxFromRealEigvecs


  !> Calculates density matrix from complex eigenvectors.
  subroutine getEDensityMtxFromComplexEigvecs(forceType, filling, eigen, kPoint, kWeight,&
      & neighborList, nNeighbor, orb, iAtomStart, iPair, img2CentCell, iCellVec, cellVec, ham,&
      & over, HSqrCplx, SSqrCplx, ERhoPrim, storeEigvecsCplx)
    integer, intent(in) :: forceType
    real(dp), intent(in) :: filling(:,:,:), eigen(:,:,:), kPoint(:,:), kWeight(:)
    type(TNeighborList), intent(in) :: neighborList
    integer, intent(in) :: nNeighbor(:)
    type(TOrbitals), intent(in) :: orb
    integer, intent(in) :: iAtomStart(:), iPair(:,:), img2CentCell(:), iCellVec(:)
    real(dp), intent(in) :: cellVec(:,:)
    real(dp), intent(in) :: ham(:,:), over(:)
    complex(dp), intent(inout) :: HSqrCplx(:,:,:,:), SSqrCplx(:,:)
    real(dp), intent(out) :: ERhoPrim(:)
    type(OFifoCplxR2), intent(inout), optional :: storeEigvecsCplx(:)

    complex(dp), allocatable :: HSqrCplx2(:,:)
    integer :: nSpin, nKPoint, sqrHamSize
    integer :: iS, iK, iS2, iK2

    nKPoint = size(eigen, dim=2)
    nSpin = size(eigen, dim=3)
    sqrHamSize = size(HSqrCplx, dim=1)

    if (forceType == 1 .or. forceType == 2 .or. forceType == 3) then
      allocate(HSqrCplx2(sqrHamSize, sqrHamSize))
    end if

    ERhoPrim(:) = 0.0_dp
    do iS = 1, nSpin
      do iK = 1, nKPoint
        if (present(storeEigvecsCplx)) then
          iK2 = 1
          iS2 = 1
          call get(storeEigvecsCplx(iS), HSqrCplx(:,:,iK2, iS2))
        else
          iS2 = iS
          iK2 = iK
        end if

        select case (forceType)

        case(0)
          ! Original (nonconsistent) scheme
          if (tDensON2) then
            call makeDensityMatrix(SSqrCplx, HSqrCplx(:,:,iK2,iS2), filling(:,iK,iS),&
                & eigen(:,iK, iS), neighborlist%iNeighbor, nNeighbor, orb, iAtomStart,&
                & img2CentCell)
          else
            call makeDensityMatrix(SSqrCplx, HSqrCplx(:,:,iK2,iS2), filling(:,iK,iS),&
                & eigen(:,iK, iS))
          end if

        case (1)
          call error("Force type 1 not implemented for complex H")

        case(2)
            ! Correct force for XLBOMD for T=0K (DHD)
          call makeDensityMatrix(HSqrCplx2, HSqrCplx(:,:,iK2,iS2), filling(:,iK,iS))
          call unpackHS(SSqrCplx, ham(:,iS), kPoint(:,iK), neighborlist%iNeighbor, nNeighbor,&
              & iCellVec, cellVec, iAtomStart, iPair, img2CentCell)
          call blockHermitianHS(SSqrCplx, iAtomStart)
          call hemm(HSqrCplx(:,:,iK2,iS2), "L", HSqrCplx2, SSqrCplx)
          call hemm(SSqrCplx, "R", HSqrCplx2, HSqrCplx(:,:,iK2,iS2), alpha=(0.5_dp, 0.0_dp))

        case(3)
            ! Correct force for XLBOMD for T <> 0K (DHS^-1 + S^-1HD)
          call makeDensityMatrix(SSqrCplx, HSqrCplx(:,:,iK2,iS2), filling(:,iK,iS))
          call unpackHS(HSqrCplx2, ham(:,iS), kPoint(:,iK), neighborlist%iNeighbor, nNeighbor,&
              & iCellVec, cellVec, iAtomStart, iPair, img2CentCell)
          call blockHermitianHS(HSqrCplx2, iAtomStart)
          call hemm(HSqrCplx(:,:,iK2,iS2), "L", SSqrCplx, HSqrCplx2)
          call unpackHS(SSqrCplx, over, kPoint(:,iK), neighborlist%iNeighbor, nNeighbor, iCellVec,&
              & cellVec, iAtomStart, iPair, img2CentCell)
          call hermatinv(SSqrCplx)
          call hemm(HSqrCplx2, "R", SSqrCplx, HSqrCplx(:,:,iK2,iS2), alpha=(0.5_dp, 0.0_dp))
          SSqrCplx(:,:) = HSqrCplx2 + transpose(conjg(HSqrCplx2))
        end select

        call packHS(ERhoPrim, SSqrCplx, kPoint(:,iK), kWeight(iK), neighborList%iNeighbor,&
            & nNeighbor, orb%mOrb, iCellVec, cellVec, iAtomStart, iPair, img2CentCell)
      end do
    end do

  end subroutine getEDensityMtxFromComplexEigvecs


  !> Calculates density matrix from Pauli-type two component eigenvectors.
  subroutine getEDensityMtxFromPauliEigvecs(filling, eigen, kPoint, kWeight, neighborList,&
      & nNeighbor, orb, iAtomStart, iPair, img2CentCell, iCellVec, cellVec, tRealHS, HSqrCplx,&
      & SSqrCplx, ERhoPrim, storeEigvecsCplx)
    real(dp), intent(in) :: filling(:,:,:), eigen(:,:,:), kPoint(:,:), kWeight(:)
    type(TNeighborList), intent(in) :: neighborList
    integer, intent(in) :: nNeighbor(:)
    type(TOrbitals), intent(in) :: orb
    integer, intent(in) :: iAtomStart(:), iPair(:,:), img2CentCell(:), iCellVec(:)
    real(dp), intent(in) :: cellVec(:,:)
    logical, intent(in) :: tRealHS
    complex(dp), intent(inout) :: HSqrCplx(:,:,:,:), SSqrCplx(:,:)
    real(dp), intent(out) :: ERhoPrim(:)
    type(OFifoCplxR2), intent(inout), optional :: storeEigvecsCplx(:)

    integer :: nKPoint
    integer :: iK, iK2

    nKPoint = size(eigen, dim=2)

    ERhoPrim(:) = 0.0_dp
    do iK = 1, nKPoint
      if (present(storeEigvecsCplx)) then
        iK2 = 1
        call get(storeEigvecsCplx(1), HSqrCplx(:,:,iK2,1))
      else
        iK2 = iK
      end if
      call makeDensityMatrix(SSqrCplx, HSqrCplx(:,:,iK2,1), filling(:,iK,1), eigen(:,iK,1))
      if (tRealHS) then
        call packERho(ERhoPrim, SSqrCplx, neighborList%iNeighbor, nNeighbor, orb%mOrb, iAtomStart,&
            & iPair, img2CentCell)
      else
        call packERho(ERhoPrim, SSqrCplx, kPoint(:,iK), kWeight(iK), neighborList%iNeighbor,&
            & nNeighbor, orb%mOrb, iCellVec, cellVec, iAtomStart, iPair, img2CentCell)
      end if
    end do

  end subroutine getEDensityMtxFromPauliEigvecs


  !> Calculates the gradients
  subroutine getGradients(tScc, tEField, tXlbomd, nonSccDeriv, Efield, rhoPrim, ERhoPrim, qOutput,&
      & q0, skHamCont, skOverCont, pRepCont, neighborList,&
      & nNeighbor, species, img2CentCell, iPair, orb, potential, coord, dispersion, &
      & derivs, iRhoPrim, thirdOrd, chrgForces)
    logical, intent(in) :: tScc, tEField, tXlbomd
    type(NonSccDiff), intent(in) :: nonSccDeriv
    real(dp), intent(in) :: Efield(:), rhoPrim(:,:), ERhoPrim(:), qOutput(:,:,:), q0(:,:,:)
    type(OSlakoCont), intent(in) :: skHamCont, skOverCont
    type(ORepCont), intent(in) :: pRepCont
    type(TNeighborList), intent(in) :: neighborList
    integer, intent(in) :: nNeighbor(:), species(:), img2CentCell(:), iPair(:,:)
    type(TOrbitals), intent(in) :: orb
    type(TPotentials), intent(in) :: potential
    real(dp), intent(in) :: coord(:,:)
    ! Workaround:ifort 17.0: Pass as allocatable instead of optional to prevent segfault
    class(DispersionIface), intent(inout), allocatable :: dispersion
    real(dp), intent(out) :: derivs(:,:)
    real(dp), intent(in), optional :: iRhoPrim(:,:)
    type(ThirdOrder), intent(inout), optional :: thirdOrd
    real(dp), intent(out), optional :: chrgForces(:,:)

    real(dp), allocatable :: tmpDerivs(:,:)
    logical :: tImHam, tExtChrg
    integer :: ii

    tImHam = present(iRhoPrim)
    tExtChrg = present(chrgForces)

    derivs(:,:) = 0.0_dp

    if (.not. (tSCC .or. tEField)) then
      ! No external or internal potentials
      if (tImHam) then
        call derivative_shift(derivs, nonSccDeriv, rhoPrim, iRhoPrim, ERhoPrim, skHamCont,&
            & skOverCont, coord, species, neighborList%iNeighbor, nNeighbor, img2CentCell,&
            & iPair, orb, potential%intBlock, potential%iorbitalBlock)
      else
        call derivative_shift(derivs, nonSccDeriv, rhoPrim(:,1), ERhoPrim, skHamCont, skOverCont,&
            & coord, species, neighborList%iNeighbor, nNeighbor, img2CentCell, iPair, orb)
      end if
    else
      if (tImHam) then
        call derivative_shift(derivs, nonSccDeriv, rhoPrim, iRhoPrim, ERhoPrim, skHamCont,&
            & skOverCont, coord, species, neighborList%iNeighbor, nNeighbor, img2CentCell, iPair,&
            & orb, potential%intBlock, potential%iorbitalBlock)
      else
        call derivative_shift(derivs, nonSccDeriv, rhoPrim, ERhoPrim, skHamCont, skOverCont, coord,&
            & species, neighborList%iNeighbor, nNeighbor, img2CentCell, iPair, orb,&
            & potential%intBlock)
      end if

      if (tExtChrg) then
        chrgForces(:,:) = 0.0_dp
        if (tXlbomd) then
          call error("XLBOMD does not work with external charges yet!")
        else
          call addForceDCSCC(derivs, species, neighborList%iNeighbor, img2CentCell, coord,&
              & chrgForces)
        end if
      elseif (tSCC) then
        if (tXlbomd) then
          call addForceDCSCC_Xlbomd(species, orb, neighborList%iNeighbor, img2CentCell, coord,&
              & qOutput, q0, derivs)
        else
          call addForceDCSCC(derivs, species, neighborList%iNeighbor, img2CentCell, coord)
        end if
      end if

      if (present(thirdOrd)) then
        if (tXlbomd) then
          call thirdOrd%addGradientDcXlbomd(neighborList, species, coord, img2CentCell, qOutput,&
              & q0, orb, derivs)
        else
          call thirdOrd%addGradientDc(neighborList, species, coord, img2CentCell, derivs)
        end if
      end if

      if (tEField) then
        do ii = 1, 3
          derivs(ii,:) = derivs(ii,:) - sum(q0(:,:,1) - qOutput(:,:,1), dim=1) * EField(ii)
        end do
      end if
    end if

    if (allocated(dispersion)) then
      call dispersion%addGradients(derivs)
    end if

    allocate(tmpDerivs(3, size(q0, dim=2)))
    call getERepDeriv(tmpDerivs, coord, nNeighbor, neighborList%iNeighbor, species, pRepCont,&
        & img2CentCell)
    derivs(:,:) = derivs + tmpDerivs

  end subroutine getGradients


  !> Calculates stress tensor and lattice derivatives.
  subroutine getStress(tScc, tEField, nonSccDeriv, EField, rhoPrim, ERhoPrim, qOutput, q0,&
      & skHamCont, skOverCont, pRepCont, neighborList, nNeighbor, species, img2CentCell, iPair,&
      & orb, potential, coord, latVec, invLatVec, cellVol, coord0, dispersion, totalStress,&
      & totalLatDeriv, cellPressure, iRhoPrim)
    logical, intent(in) :: tScc, tEField
    type(NonSccDiff), intent(in) :: nonSccDeriv
    real(dp), intent(in) :: Efield(:), rhoPrim(:,:), ERhoPrim(:), qOutput(:,:,:), q0(:,:,:)
    type(OSlakoCont), intent(in) :: skHamCont, skOverCont
    type(ORepCont), intent(in) :: pRepCont
    type(TNeighborList), intent(in) :: neighborList
    integer, intent(in) :: nNeighbor(:), species(:), img2CentCell(:), iPair(:,:)
    type(TOrbitals), intent(in) :: orb
    type(TPotentials), intent(in) :: potential
    real(dp), intent(in) :: coord(:,:), latVec(:,:), invLatVec(:,:), cellVol
    real(dp), intent(inout) :: coord0(:,:)
    ! Workaround:ifort 17.0: Pass as allocatable instead of optional to prevent segfault
    class(DispersionIface), allocatable, intent(inout) :: dispersion
    real(dp), intent(out) :: totalStress(:,:), totalLatDeriv(:,:), cellPressure
    real(dp), intent(in), optional :: iRhoPrim(:,:)

    real(dp) :: tmpStress(3, 3)
    logical :: tImHam

    tImHam = present(iRhoPrim)

    if (tSCC) then
      if (tImHam) then
        call getBlockiStress(totalStress, nonSccDeriv, rhoPrim, iRhoPrim, ERhoPrim, skHamCont,&
            & skOverCont, coord, species, neighborList%iNeighbor, nNeighbor, img2CentCell, iPair,&
            & orb, potential%intBlock, potential%iorbitalBlock, cellVol)
      else
        call getBlockStress(totalStress, nonSccDeriv, rhoPrim, ERhoPrim, skHamCont, skOverCont,&
            & coord, species, neighborList%iNeighbor, nNeighbor, img2CentCell, iPair, orb,&
            & potential%intBlock, cellVol)
      end if
      call addStressDCSCC(totalStress,species,neighborList%iNeighbor, img2CentCell,coord)
    else
      if (tImHam) then
        call getBlockiStress(totalStress, nonSccDeriv, rhoPrim, iRhoPrim, ERhoPrim, skHamCont,&
            & skOverCont, coord, species, neighborList%iNeighbor, nNeighbor, img2CentCell, iPair,&
            & orb, potential%intBlock, potential%iorbitalBlock, cellVol)
      else
        call getNonSCCStress(totalStress, nonSccDeriv, rhoPrim(:,1), ERhoPrim, skHamCont,&
            & skOverCont, coord, species, neighborList%iNeighbor, nNeighbor, img2CentCell, iPair,&
            & orb, cellVol)
      end if
    end if

    if (allocated(dispersion)) then
      call dispersion%getStress(tmpStress)
      totalStress(:,:) = totalStress + tmpStress
    end if

    if (tEField) then
      call getEFieldStress(latVec, cellVol, q0, qOutput, Efield, coord0, tmpStress)
      totalStress(:,:) = totalStress + tmpStress
    end if

    call getRepulsiveStress(tmpStress, coord, nNeighbor, neighborList%iNeighbor, species,&
        & img2CentCell, pRepCont, cellVol)
    totalStress(:,:) = totalStress + tmpStress

    cellPressure = (totalStress(1,1) + totalStress(2,2) + totalStress(3,3)) / 3.0_dp
    totalLatDeriv(:,:) = -cellVol * matmul(totalStress, invLatVec)

  end subroutine getStress


  !> Calculates stress from external electric field.
  subroutine getEFieldStress(latVec, cellVol, q0, qOutput, Efield, coord0, stress)
    real(dp), intent(in) :: latVec(:,:), cellVol, q0(:,:,:), qOutput(:,:,:), Efield(:)
    real(dp), intent(inout) :: coord0(:,:)
    real(dp), intent(out) :: stress(:,:)

    real(dp) :: latDerivs(3,3)
    integer :: nAtom
    integer :: iAtom, ii, jj

    nAtom = size(coord0, dim=2)

    latDerivs(:,:) = 0.0_dp
    call cart2frac(coord0, latVec)
    do iAtom = 1, nAtom
      do ii = 1, 3
        do jj = 1, 3
          latDerivs(jj,ii) =  latDerivs(jj,ii)&
              & - sum(q0(:,iAtom,1) - qOutput(:,iAtom,1), dim=1) * EField(ii) * coord0(jj,iAtom)
        end do
      end do
    end do
    call frac2cart(coord0, latVec)
    stress(:,:) = -matmul(latDerivs, transpose(latVec)) / cellVol

  end subroutine getEFieldStress


  !> Prints cell volume.
  subroutine printVolume(cellVol)
    real(dp), intent(in) :: cellVol

    write(stdOut, format2Ue) 'Volume', cellVol, 'au^3', (Bohr__AA**3) * cellVol, 'A^3'
  end subroutine printVolume


  !> Prints pressure and free energy.
  subroutine printPressureAndFreeEnergy(pressure, cellPressure, EGibbs)
    real(dp), intent(in) :: pressure, cellPressure, EGibbs

    write(stdOut, format2Ue) 'Pressure', cellPressure, 'au', cellPressure * au__pascal, 'Pa'
    if (abs(pressure) > epsilon(1.0_dp)) then
      write(stdOut, format2U) "Gibbs free energy", EGibbs, 'H', Hartree__eV * EGibbs, 'eV'
    end if

  end subroutine printPressureAndFreeEnergy


  !> Removes forces components along constraint directions
  subroutine constrainForces(conAtom, conVec, totalDerivs)
    integer, intent(in) :: conAtom(:)
    real(dp), intent(in) :: conVec(:,:)
    real(dp), intent(inout) :: totalDerivs(:,:)

    integer :: ii, iAtom

    ! Set force components along constraint vectors zero
    do ii = 1, size(conAtom)
      iAtom = conAtom(ii)
      totalDerivs(:,iAtom) = totalDerivs(:,iAtom)&
          & - conVec(:,ii) * dot_product(conVec(:,ii), totalDerivs(:,iAtom))
    end do

  end subroutine constrainForces


  !> Writes maximal force component.
  subroutine printMaxForce(maxForce)
    real(dp), intent(in) :: maxForce

    write(stdOut, "(A, ':', T30, E20.6)") "Maximal force component", maxForce

  end subroutine printMaxForce


  !> Flattens lattice components and applies lattice optimisation constraints.
  subroutine constrainLatticeDerivs(totalLatDerivs, normLatVecs, tLatOptFixAng,&
      & tLatOptFixLen, tLatOptIsotropic, constrLatDerivs)
    real(dp), intent(in) :: totalLatDerivs(:,:), normLatVecs(:,:)
    logical, intent(in) :: tLatOptFixAng, tLatOptFixLen(:), tLatOptIsotropic
    real(dp), intent(out) :: constrLatDerivs(:)

    real(dp) :: tmpLatDerivs(3, 3)
    integer :: ii

    tmpLatDerivs(:,:) = totalLatDerivs
    constrLatDerivs = reshape(tmpLatDerivs, [9])
    if (tLatOptFixAng) then
      ! project forces to be along original lattice
      tmpLatDerivs(:,:) = tmpLatDerivs * normLatVecs
      constrLatDerivs(:) = 0.0_dp
      if (any(tLatOptFixLen)) then
        do ii = 1, 3
          if (.not. tLatOptFixLen(ii)) then
            constrLatDerivs(ii) = sum(tmpLatDerivs(:,ii))
          end if
        end do
      else
        constrLatDerivs(1:3) = sum(tmpLatDerivs, dim=1)
      end if
    elseif (tLatOptIsotropic) then
      tmpLatDerivs(:,:) = tmpLatDerivs * normLatVecs
      constrLatDerivs(:) = 0.0_dp
      constrLatDerivs(1) = sum(tmpLatDerivs)
    end if

  end subroutine constrainLatticeDerivs


  !> Print maximal lattice force component
  subroutine printMaxLatticeForce(maxLattForce)
    real(dp), intent(in) :: maxLattForce

    write(stdOut, format1Ue) "Maximal Lattice force component", maxLattForce, 'au'

  end subroutine printMaxLatticeForce


  !> Write out charges.
  subroutine writeCharges(qInput, fChargeIn, orb, qBlockIn, qiBlockIn)
    real(dp), intent(in) :: qInput(:,:,:)
    character(*), intent(in) :: fChargeIn
    type(TOrbitals), intent(in) :: orb
    real(dp), intent(in), optional :: qBlockIn(:,:,:,:), qiBlockIn(:,:,:,:)

    call writeQToFile(qInput, fChargeIn, orb, qBlockIn, qiBlockIn)
    write(stdOut, "(A,A)") '>> Charges saved for restart in ', trim(fChargeIn)

  end subroutine writeCharges


  !> Unfold contrained lattice vectors to full one.
  subroutine unconstrainLatticeVectors(constrLatVecs, origLatVecs, tLatOptFixAng, tLatOptFixLen,&
      & tLatOptIsotropic, newLatVecs)
    real(dp), intent(in) :: constrLatVecs(:), origLatVecs(:,:)
    logical, intent(in) :: tLatOptFixAng, tLatOptFixLen(:), tLatOptIsotropic
    real(dp), intent(out) :: newLatVecs(:,:)

    real(dp) :: tmpLatVecs(9)
    integer :: ii

    tmpLatVecs(:) = constrLatVecs
    if (tLatOptFixAng) then
      ! Optimization uses scaling factor of lattice vectors
      if (any(tLatOptFixLen)) then
        do ii = 3, 1, -1
          if (.not. tLatOptFixLen(ii)) then
            tmpLatVecs(3 * ii - 2 : 3 * ii) =  tmpLatVecs(ii) * origLatVecs(:,ii)
          else
            tmpLatVecs(3 * ii - 2 : 3 * ii) =  origLatVecs(:,ii)
          end if
        end do
      else
        tmpLatVecs(7:9) =  tmpLatVecs(3) * origLatVecs(:,3)
        tmpLatVecs(4:6) =  tmpLatVecs(2) * origLatVecs(:,2)
        tmpLatVecs(1:3) =  tmpLatVecs(1) * origLatVecs(:,1)
      end if
    else if (tLatOptIsotropic) then
      ! Optimization uses scaling factor unit cell
      do ii = 3, 1, -1
        tmpLatVecs(3 * ii - 2 : 3 * ii) =  tmpLatVecs(1) * origLatVecs(:,ii)
      end do
    end if
    newLatVecs(:,:) = reshape(tmpLatVecs, [3, 3])

  end subroutine unconstrainLatticeVectors


  !> Returns the coordinates for the next Hessian calculation step.
  subroutine getNextDerivStep(derivDriver, derivs, indMovedAtoms, coords, tGeomEnd)
    type(OnumDerivs), intent(inout) :: derivDriver
    real(dp), intent(in) :: derivs(:,:)
    integer, intent(in) :: indMovedAtoms(:)
    real(dp), intent(out) :: coords(:,:)
    logical, intent(out) :: tGeomEnd

    real(dp) :: newCoords(3, size(indMovedAtoms))

    call next(derivDriver, newCoords, derivs(:, indMovedAtoms), tGeomEnd)
    coords(:, indMovedAtoms) = newCoords

  end subroutine getNextDerivStep


  !> Returns the coordinates for the next coordinate optimisation step.
  subroutine getNextCoordinateOptStep(pGeoCoordOpt, EMermin, totalDerivs, indMovedAtom, coords0,&
      & diffGeo, tCoordEnd)
    type(OGeoOpt), intent(inout) :: pGeoCoordOpt
    real(dp), intent(in) :: EMermin, totalDerivs(:,:)
    integer, intent(in) :: indMovedAtom(:)
    real(dp), intent(inout) :: coords0(:,:)
    real(dp), intent(out) :: diffGeo
    logical, intent(out) :: tCoordEnd

    real(dp) :: totalDerivsMoved(3 * size(indMovedAtom))
    real(dp), target :: newCoordsMoved(3 * size(indMovedAtom))
    real(dp), pointer :: pNewCoordsMoved(:,:)

    totalDerivsMoved(:) = reshape(totalDerivs(:, indMovedAtom), [3 * size(indMovedAtom)])
    call next(pGeoCoordOpt, EMermin, totalDerivsMoved, newCoordsMoved, tCoordEnd)
    pNewCoordsMoved(1:3, 1:size(indMovedAtom)) => newCoordsMoved(1 : 3 * size(indMovedAtom))
    diffGeo = maxval(abs(pNewCoordsMoved - coords0(:, indMovedAtom)))
    coords0(:, indMovedAtom) = pNewCoordsMoved

  end subroutine getNextCoordinateOptStep


  !> Returns the coordinates and lattice vectors for the next lattice optimisation step.
  subroutine getNextLatticeOptStep(pGeoLatOpt, EGibbs, constrLatDerivs, origLatVec, tLatOptFixAng,&
      & tLatOptFixLen, tLatOptIsotropic, indMovedAtom, latVec, coord0, diffGeo, tGeomEnd)
    type(OGeoOpt), intent(inout) :: pGeoLatOpt
    real(dp), intent(in) :: EGibbs, constrLatDerivs(:), origLatVec(:,:)
    logical, intent(in) :: tLatOptFixAng, tLatOptFixLen(:), tLatOptIsotropic
    integer, intent(in) :: indMovedAtom(:)
    real(dp), intent(inout) :: latVec(:,:), coord0(:,:)
    real(dp), intent(out) :: diffGeo
    logical, intent(out) :: tGeomEnd

    real(dp) :: newLatVecsFlat(9), newLatVecs(3, 3), oldMovedCoords(3, size(indMovedAtom))

    call next(pGeoLatOpt, EGibbs, constrLatDerivs, newLatVecsFlat,tGeomEnd)
    call unconstrainLatticeVectors(newLatVecsFlat, origLatVec, tLatOptFixAng, tLatOptFixLen,&
        & tLatOptIsotropic, newLatVecs)
    oldMovedCoords(:,:) = coord0(:, indMovedAtom)
    call cart2frac(coord0, latVec)
    latVec(:,:) = newLatVecs
    call frac2cart(coord0, latVec)
    diffGeo = max(maxval(abs(newLatVecs - latVec)),&
        & maxval(abs(oldMovedCoords - coord0(:, indMovedAtom))))

  end subroutine getNextLatticeOptStep


  !> Delivers data for next MD step (and updates data depending on velocities of current step)
  subroutine getNextMdStep(pMdIntegrator, pMdFrame, temperatureProfile, totalDeriv, movedMass,&
      & mass, cellVol, invLatVec, species0, indMovedAtom, tStress, tBarostat, energy, coord0,&
      & latVec, cellPressure, totalStress, totalLatDeriv, velocities, kT)
    type(OMdIntegrator), intent(inout) :: pMdIntegrator
    type(OMdCommon), intent(in) :: pMdFrame
    type(OTempProfile), allocatable, intent(inout) :: temperatureProfile
    real(dp), intent(in) :: totalDeriv(:,:), movedMass(:,:), mass(:), cellVol, invLatVec(:,:)
    integer, intent(in) :: species0(:), indMovedAtom(:)
    logical, intent(in) :: tStress, tBarostat
    type(TEnergies), intent(inout) :: energy
    real(dp), intent(inout) :: coord0(:,:), latVec(:,:)
    real(dp), intent(inout) :: cellPressure, totalStress(:,:), totalLatDeriv(:,:)
    real(dp), intent(out) :: velocities(:,:), kT

    real(dp) :: movedAccel(3, size(indMovedAtom)), movedVelo(3, size(indMovedAtom))
    real(dp) :: movedCoords(3, size(indMovedAtom))
    real(dp) :: kineticStress(3, 3)

    movedAccel(:,:) = -totalDeriv(:, indMovedAtom) / movedMass
    call next(pMdIntegrator, movedAccel, movedCoords, movedVelo)
    coord0(:, indMovedAtom) = movedCoords
    velocities(:,:) = 0.0_dp
    velocities(:, indMovedAtom) = movedVelo(:,:)

    if (allocated(temperatureProfile)) then
      call next(temperatureProfile)
    end if
    call evalKE(energy%Ekin, movedVelo, movedMass(1,:))
    call evalkT(pMdFrame, kT, movedVelo, movedMass(1,:))
    energy%EMerminKin = energy%EMermin + energy%Ekin
    energy%EGibbsKin = energy%EGibbs + energy%Ekin

    if (tStress) then
      ! contribution from kinetic energy in MD, now that velocities for this geometry step are
      ! available
      call getKineticStress(kineticStress, mass, species0, velocities, cellVol)
      totalStress = totalStress + kineticStress
      cellPressure = (totalStress(1,1) + totalStress(2,2) + totalStress(3,3)) / 3.0_dp
      totalLatDeriv = -cellVol * matmul(totalStress, invLatVec)
    end if

    if (tBarostat) then
      call rescale(pMDIntegrator, coord0, latVec, totalStress)
    end if

  end subroutine getNextMdStep


  !> Prints out info about current MD step.
  subroutine printMdInfo(tSetFillingTemp, tEField, tPeriodic, tempElec, absEField, kT,&
      & cellPressure, pressure, energy)
    logical, intent(in) :: tSetFillingTemp, tEFIeld, tPeriodic
    real(dp), intent(in) :: tempElec, absEField, kT, cellPressure, pressure
    type(TEnergies), intent(in) :: energy

    if (tSetFillingTemp) then
      write(stdOut, format2U) 'Electronic Temperature:', tempElec, 'H', tempElec / Boltzmann, 'K'
    end if
    if (tEfield) then
      write(stdOut, format1U1e) 'External E field', absEField, 'au', absEField * au__V_m, 'V/m'
    end if
    write(stdOut, format2U) "MD Temperature:", kT, "H", kT / Boltzmann, "K"
    write(stdOut, format2U) "MD Kinetic Energy", energy%Ekin, "H", Hartree__eV * energy%Ekin, "eV"
    write(stdOut, format2U) "Total MD Energy", energy%EMerminKin, "H",&
        & Hartree__eV * energy%EMerminKin, "eV"
    if (tPeriodic) then
      write(stdOut, format2Ue) 'Pressure', cellPressure, 'au', cellPressure * au__pascal, 'Pa'
      if (abs(pressure) < epsilon(1.0_dp)) then
        write(stdOut, format2U) 'Gibbs free energy including KE', energy%EGibbsKin, 'H',&
            & Hartree__eV * energy%EGibbsKin, 'eV'
      end if
    end if

  end subroutine printMdInfo


  !> Write out final status of the geometry driver.
  subroutine writeFinalDriverStatus(tGeoOpt, tGeomEnd, tMd, tDerivs)
    logical, intent(in) :: tGeoOpt, tGeomEnd, tMd, tDerivs

    if (tGeoOpt) then
      if (tGeomEnd) then
        write(stdOut, "(/, A)") "Geometry converged"
      else
        call warning("!!! Geometry did NOT converge!")
      end if
    elseif (tMD) then
      if (tGeomEnd) then
        write(stdOut, "(/, A)") "Molecular dynamics completed"
      else
        call warning("!!! Molecular dynamics terminated abnormally!")
      end if
    elseif (tDerivs) then
      if (tGeomEnd) then
        write(stdOut, "(/, A)") "Second derivatives completed"
      else
        call warning("!!! Second derivatives terminated abnormally!")
      end if
    end if

  end subroutine writeFinalDriverStatus


  !> Calculates and prints Pipek-Mezey localisation
  subroutine calcPipekMezeyLocalisation(pipekMezey, nEl, filling, over, kPoint, kWeight,&
      & neighborList, nNeighbor, iAtomStart, iPair, img2CentCell, iCellVec, cellVec, fdEigvec,&
      & runId, orb, species, speciesName, localisation, HSqrReal, SSqrReal, HsqrCplx, SSqrCplx)
    type(TPipekMezey), intent(in) :: pipekMezey
    real(dp), intent(in) :: nEl(:), filling(:,:,:), over(:), kPoint(:,:), kWeight(:)
    type(TNeighborList), intent(in) :: neighborList
    integer, intent(in) :: nNeighbor(:), iAtomStart(:), iPair(:,:), img2CentCell(:), iCellVec(:)
    real(dp), intent(in) :: cellVec(:,:)
    integer, intent(in) :: fdEigvec, runId
    type(TOrbitals), intent(in) :: orb
    integer, intent(in) :: species(:)
    character(*), intent(in) :: speciesName(:)
    real(dp), intent(out) :: localisation
    real(dp), intent(inout), optional :: HSqrReal(:,:,:), SSqrReal(:,:)
    complex(dp), intent(inout), optional :: HSqrCplx(:,:,:,:), SSqrCplx(:,:)

    integer :: nFilledLev, nAtom, nSpin
    integer :: iSpin

    nAtom = size(orb%nOrbAtom)
    nSpin = size(nEl)

    if (any(abs(mod(filling, real(3 - nSpin, dp))) > elecTolMax)) then
      call warning("Fractional occupations present for electron localisation")
    end if

    if (present(HSqrReal)) then
      call unpackHS(SSqrReal,over,neighborList%iNeighbor, nNeighbor, iAtomStart, iPair,&
          & img2CentCell)
      do iSpin = 1, nSpin
        nFilledLev = floor(nEl(iSpin) / real(3 - nSpin, dp))
        localisation = pipekMezey%getLocalisation(HSqrReal(:, 1:nFilledLev, iSpin), SSqrReal,&
            & iAtomStart)
        write(stdOut, "(A, E15.8)") 'Original localisation', localisation
        call pipekMezey%calcCoeffs(HSqrReal(:, 1:nFilledLev, iSpin), SSqrReal, iAtomStart)
        localisation = pipekMezey%getLocalisation(HSqrReal(:,1:nFilledLev,iSpin), SSqrReal,&
            & iAtomStart)
        write(stdOut, "(A, E20.12)") 'Final localisation ', localisation
      end do

      call writeEigvecs(fdEigvec, runId, nAtom, nSpin, neighborList, nNeighbor, iAtomStart, iPair,&
          & img2CentCell, orb, species, speciesName, over, HSqrReal, SSqrReal, fileName="localOrbs")
    else
      do iSpin = 1, nSpin
        nFilledLev = floor(nEl(iSpin) / real( 3 - nSpin, dp))
        localisation = sum(pipekMezey%getLocalisation(&
            & HSqrCplx(:,:nFilledLev,:,iSpin), SSqrCplx, over, kpoint, kweight, neighborList,&
            & nNeighbor, iCellVec, cellVec, iAtomStart, iPair, img2CentCell))
        write(stdOut, "(A, E20.12)") 'Original localisation', localisation
        call pipekMezey%calcCoeffs(HSqrCplx(:,:nFilledLev,:,iSpin), SSqrCplx, over, kpoint,&
            & kweight, neighborList, nNeighbor, iCellVec, cellVec, iAtomStart, iPair, img2CentCell)
        localisation = sum(pipekMezey%getLocalisation(&
            & HSqrCplx(:,:nFilledLev,:,iSpin), SSqrCplx, over, kpoint, kweight, neighborList,&
            & nNeighbor, iCellVec, cellVec, iAtomStart, iPair, img2CentCell))
        write(stdOut, "(A, E20.12)") 'Final localisation', localisation
      end do

      call writeEigvecs(fdEigvec, runId, nAtom, nSpin, neighborList, nNeighbor, cellVec, iCellVec,&
          & iAtomStart, iPair, img2CentCell, orb, species, speciesName, over, kpoint, HSqrCplx,&
          & SSqrCplx, fileName="localOrbs")
    end if

  end subroutine calcPipekMezeyLocalisation


end module main
