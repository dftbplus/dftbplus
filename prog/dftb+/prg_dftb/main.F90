!--------------------------------------------------------------------------------------------------!
!  DFTB+: general package for performing fast atomistic simulations                                !
!  Copyright (C) 2018  DFTB+ developers group                                                      !
!                                                                                                  !
!  See the LICENSE file for terms of usage and distribution.                                       !
!--------------------------------------------------------------------------------------------------!

#:include 'common.fypp'

!> The main routines for DFTB+
module main
#:if WITH_MPI
  use mpifx
#:endif
#:if WITH_SCALAPACK
  use scalapackfx
  use scalafxext
#:endif
#:if WITH_PROGRESS
  use sp2progress
#:endif
  use assert
  use constants
  use globalenv
  use environment
  use densedescr
  use inputdata_module
  use nonscc
  use eigenvects
  use repulsive
  use etemp
  use populations
  use densitymatrix
  use forces
  use stress
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
  use elstatpot, only : TElStatPotentials
  use dftbp_solvers
  implicit none
  private

  public :: runDftbPlus

  !> O(N^2) density matrix creation
  logical, parameter :: tDensON2 = .false.

  !> Should further output be appended to detailed.out?
  logical, parameter :: tAppendDetailedOut = .false.


contains

  !> The main DFTB program itself
  subroutine runDftbPlus(env)
    use initprogram

    !> Environment settings
    type(TEnvironment), intent(inout) :: env

    !> energy in previous scc cycles
    real(dp) :: Eold

    !> Stress tensors for various contribution in periodic calculations
    !> Sign convention: Positive diagonal elements expand the supercell
    real(dp) :: totalStress(3,3)

    !> Derivative of total energy with respect to lattice vectors
    !> Sign convention: This is in the uphill energy direction for the lattice vectors (each row
    !> pertaining to a separate lattice vector), i.e. opposite to the force.
    !>
    !> The component of a derivative vector that is orthogonal to the plane containing the other two
    !> lattice vectors will expand (contract) the supercell if it is on the opposite (same) same
    !> side of the plane as its associated lattice vector.
    !>
    !> In the special case of cartesian axis aligned orthorhombic lattice vectors, negative diagonal
    !> elements expand the supercell.
    real(dp) :: totalLatDeriv(3,3)

    !> derivative of cell volume wrt to lattice vectors, needed for pV term
    real(dp) :: extLatDerivs(3,3)

    !> whether scc converged
    logical :: tConverged

    !> internal pressure within the cell
    real(dp) :: intPressure

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

    !> lattice vectors returned by the optimizer
    real(dp) :: constrLatDerivs(9)

    !> MD instantaneous thermal energy
    real(dp) :: tempIon

    !> external electric field
    real(dp) :: Efield(3), absEfield

    !> Difference between last calculated and new geometry.
    real(dp) :: diffGeo

    !> Loop variables
    integer :: iSCCIter

    !> Charge error in the last iterations
    real(dp) :: sccErrorQ, diffElec

    !> flag to write out geometries (and charge data if scc) when moving atoms about - in the case
    !> of conjugate gradient/steepest descent the geometries are written anyway
    logical :: tWriteRestart

    !> Minimal number of SCC iterations
    integer :: minSCCIter

    !> if scc/geometry driver should be stopped
    logical :: tStopSCC, tStopDriver

    !> Whether scc restart info should be written in current iteration
    logical :: tWriteSccRestart

    !> Whether charges should be written
    logical :: tWriteCharges

    !> locality measure for the wavefunction
    real(dp) :: localisation

    call initGeoOptParameters(tCoordOpt, nGeoSteps, tGeomEnd, tCoordStep, tStopDriver, iGeoStep,&
        & iLatGeoStep)

    minSccIter = getMinSccIters(tSccCalc, tDftbU, nSpin)

    if (tXlbomd) then
      call xlbomdIntegrator%setDefaultSCCParameters(minSCCiter, maxSccIter, sccTol)
    end if

    ! If the geometry is periodic, need to update lattice information in geometry loop
    tLatticeChanged = tPeriodic

    ! As first geometry iteration, require updates for coordinates in dependent routines
    tCoordsChanged = .true.

    ! Main geometry loop
    lpGeomOpt: do iGeoStep = 0, nGeoSteps
      call env%globalTimer%startTimer(globalTimers%preSccInit)

      call printGeoStepInfo(tCoordOpt, tLatOpt, iLatGeoStep, iGeoStep)

      tWriteRestart = env%tGlobalMaster&
          & .and. needsRestartWriting(tGeoOpt, tMd, iGeoStep, nGeoSteps, restartFreq)
      if (tMD .and. tWriteRestart) then
        call writeMdOut1(fdMd, mdOut, iGeoStep, pMDIntegrator)
      end if

      if (tLatticeChanged) then
        call handleLatticeChange(latVec, sccCalc, tStress, extPressure, mCutoff,&
            & dispersion, recVec, invLatVec, cellVol, recCellVol, extLatDerivs, cellVec, rCellVec)
      end if

      if (tCoordsChanged) then
        call handleCoordinateChange(env, coord0, latVec, invLatVec, species0, mCutoff, skRepCutoff,&
            & orb, tPeriodic, sccCalc, dispersion, thirdOrd, img2CentCell, iCellVec, neighborList,&
            & nAllAtom, coord0Fold, coord, species, rCellVec, nAllOrb, nNeighbor, ham, over, H0,&
            & rhoPrim, iRhoPrim, iHam, ERhoPrim, iSparseStart)
      end if

      if (tSccCalc) then
        call reset(pChrgMixer, nMixElements)
      end if

      call env%globalTimer%startTimer(globalTimers%sparseH0S)
      call buildH0(env, H0, skHamCont, atomEigVal, coord, nNeighbor, neighborList%iNeighbor,&
          & species, iSparseStart, orb)
      call buildS(env, over, skOverCont, coord, nNeighbor, neighborList%iNeighbor, species,&
          & iSparseStart, orb)
      call env%globalTimer%stopTimer(globalTimers%sparseH0S)

    #:if WITH_PROGRESS
      ! The following will build the inverse overlap congruence transformation for SP2
      if (tCoordsChanged .and. solver%solverType == solverTypes%progressSp2) then
        call solver%sp2Solver%buildZMatrix(over, neighborList, nNeighbor, iSparseStart,&
            & img2CentCell, denseDesc, orb)
      end if
    #:endif

      if (tSetFillingTemp) then
        call getTemperature(temperatureProfile, tempElec)
      end if

      call calcRepulsiveEnergy(coord, species, img2CentCell, nNeighbor, neighborList, pRepCont,&
          & energy%atomRep, energy%ERep)

      if (tDispersion) then
        call calcDispersionEnergy(dispersion, energy%atomDisp, energy%Edisp)
      end if

      call resetExternalPotentials(potential)
      call setUpExternalElectricField(tEField, tTDEField, tPeriodic, EFieldStrength,&
          & EFieldVector, EFieldOmega, EFieldPhase, neighborList, nNeighbor, iCellVec,&
          & img2CentCell, cellVec, deltaT, iGeoStep, coord0Fold, coord, EField,&
          & potential%extAtom(:,1), absEField)
    
      call mergeExternalPotentials(orb, species, potential)

      call initSccLoop(tSccCalc, xlbomdIntegrator, minSccIter, maxSccIter, sccTol, tConverged)

      call env%globalTimer%stopTimer(globalTimers%preSccInit)

      call env%globalTimer%startTimer(globalTimers%scc)
      lpSCC: do iSccIter = 1, maxSccIter
        call resetInternalPotentials(tDualSpinOrbit, xi, orb, species, potential)
        if (tSccCalc) then
          call getChargePerShell(qInput, orb, species, chargePerShell)
          call addChargePotentials(env, sccCalc, qInput, q0, chargePerShell, orb, species,&
              & neighborList, img2CentCell, spinW, thirdOrd, potential)
          call addBlockChargePotentials(qBlockIn, qiBlockIn, tDftbU, tImHam, species, orb,&
              & nDftbUFunc, UJ, nUJ, iUJ, niUJ, potential)
        end if
        potential%intBlock = potential%intBlock + potential%extBlock

        call getSccHamiltonian(H0, over, nNeighbor, neighborList, species, orb, iSparseStart,&
            & img2CentCell, potential, ham, iHam)

        if (tWriteRealHS .or. tWriteHS) then
          if (withMpi) then
            call error("Writing of HS not working with MPI yet")
          else
            call writeHSAndStop(tWriteHS, tWriteRealHS, tRealHS, over, neighborList, nNeighbor,&
                & denseDesc%iAtomStart, iSparseStart, img2CentCell, kPoint, iCellVec, cellVec,&
                & ham, iHam)
          end if
        end if

        call getDensity(env, denseDesc, ham, over, neighborList, nNeighbor, iSparseStart,&
            & img2CentCell, iCellVec, cellVec, kPoint, kWeight, orb, species, solver, tRealHS,&
            & tSpinSharedEf, tSpinOrbit, tDualSpinOrbit, tFillKSep, tFixEf, tMulliken, iDistribFn,&
            & tempElec, nEl, parallelKS, Ef, energy, eigen, filling, rhoPrim, Eband, TS, E0, iHam,&
            & xi, orbitalL, HSqrReal, SSqrReal, eigvecsReal, iRhoPrim, HSqrCplx, SSqrCplx,&
            & eigvecsCplx, rhoSqrReal)

        if (tWriteBandDat) then
          call writeBandOut(fdBand, bandOut, eigen, filling, kWeight)
        end if

        if (tMulliken) then
          call getMullikenPopulation(rhoPrim, over, orb, neighborList, nNeighbor, img2CentCell,&
              & iSparseStart, qOutput, iRhoPrim=iRhoPrim, qBlock=qBlockOut, qiBlock=qiBlockOut)
        end if

        ! For non-dual spin-orbit orbitalL is determined during getDensity() call above
        if (tDualSpinOrbit) then
          call getLDual(orbitalL, qiBlockOut, orb, species)
        end if

        ! Note: if XLBOMD is active, potential created with input charges is needed later,
        ! therefore it should not be overwritten here.
        if (tSccCalc .and. .not. tXlbomd) then
          call resetInternalPotentials(tDualSpinOrbit, xi, orb, species, potential)
          call getChargePerShell(qOutput, orb, species, chargePerShell)
          call addChargePotentials(env, sccCalc, qOutput, q0, chargePerShell, orb, species,&
              & neighborList, img2CentCell, spinW, thirdOrd, potential)
          call addBlockChargePotentials(qBlockOut, qiBlockOut, tDftbU, tImHam, species, orb,&
              & nDftbUFunc, UJ, nUJ, iUJ, niUJ, potential)
          potential%intBlock = potential%intBlock + potential%extBlock

        end if

        call getEnergies(sccCalc, qOutput, q0, chargePerShell, species, tEField, tXlbomd,&
            & tDftbU, tDualSpinOrbit, rhoPrim, H0, orb, neighborList, nNeighbor, img2CentCell,&
            & iSparseStart, cellVol, extPressure, TS, potential, energy, thirdOrd, qBlockOut,&
            & qiBlockOut, nDftbUFunc, UJ, nUJ, iUJ, niUJ, xi)

        tStopScc = hasStopFile(fStopScc)

        if (tSccCalc) then
          call getNextInputCharges(env, pChrgMixer, qOutput, qOutRed, orb, nIneqOrb, iEqOrbitals,&
              & iGeoStep, iSccIter, minSccIter, maxSccIter, sccTol, tStopScc, tDftbU, tReadChrg,&
              & qInput, qInpRed, sccErrorQ, tConverged, qBlockOut, iEqBlockDftbU, qBlockIn,&
              & qiBlockOut, iEqBlockDftbULS, species0, nUJ, iUJ, niUJ, qiBlockIn)
          call getSccInfo(iSccIter, energy%Eelec, Eold, diffElec)
          call printSccInfo(tDftbU, iSccIter, energy%Eelec, diffElec, sccErrorQ)
          tWriteSccRestart = env%tGlobalMaster .and. &
              & needsSccRestartWriting(restartFreq, iGeoStep, iSccIter, minSccIter, maxSccIter,&
              & tMd, tGeoOpt, tDerivs, tConverged, tReadChrg, tStopScc)
          if (tWriteSccRestart) then
            call writeCharges(fCharges, fdCharges, tWriteChrgAscii, orb, qInput, qBlockIn,&
                & qiBlockIn)
          end if
        end if

        if (tWriteDetailedOut) then
          call writeDetailedOut1(fdDetailedOut, userOut, tAppendDetailedOut, iDistribFn, nGeoSteps,&
              & iGeoStep, tMD, tDerivs, tCoordOpt, tLatOpt, iLatGeoStep, iSccIter, energy,&
              & diffElec, sccErrorQ, indMovedAtom, pCoord0Out, q0, qInput, qOutput, eigen, filling,&
              & orb, species, tDFTBU, tImHam, tPrintMulliken, orbitalL, qBlockOut, Ef, Eband, TS,&
              & E0, extPressure, cellVol, tAtomicEnergy, tDispersion, tEField, tPeriodic, nSpin,&
              & tSpinOrbit, tSccCalc, invLatVec, kPoint)
        end if

        if (tConverged .or. tStopScc) then
          exit lpSCC
        end if

      end do lpSCC
      call env%globalTimer%stopTimer(globalTimers%scc)

      call env%globalTimer%startTimer(globalTimers%postSCC)
      if (tLinResp) then
        if (withMpi) then
          call error("Linear response calc. does not work with MPI yet")
        end if
        call ensureLinRespConditions(t3rd, tRealHS, tPeriodic, tForces)
        call calculateLinRespExcitations(env, lresp, parallelKS, sccCalc, qOutput, q0, over,&
            & eigvecsReal, eigen(:,1,:), filling(:,1,:), coord0, species, speciesName, orb,&
            & skHamCont, skOverCont, fdAutotest, autotestTag, fdEigvec, runId, neighborList,&
            & nNeighbor, denseDesc, iSparseStart, img2CentCell, tWriteAutotest, tForces,&
            & tLinRespZVect, tPrintExcitedEigvecs, tPrintEigvecsTxt, nonSccDeriv, energy, SSqrReal,&
            & rhoSqrReal, excitedDerivs, occNatural)
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
            & nNeighbor, species, iSparseStart, img2CentCell)
      #:endcall DEBUG_CODE
      end if

      call env%globalTimer%startTimer(globalTimers%eigvecWriting)
      if (tPrintEigVecs) then
        call writeEigenvectors(env, fdEigvec, runId, neighborList, nNeighbor, cellVec, iCellVec,&
            & denseDesc, iSparseStart, img2CentCell, species, speciesName, orb, kPoint, over,&
            & parallelKS, tPrintEigvecsTxt, eigvecsReal, SSqrReal, eigvecsCplx, SSqrCplx)
      end if

      if (tProjEigenvecs) then
        call writeProjectedEigenvectors(env, regionLabels, fdProjEig, eigen, neighborList,&
              & nNeighbor, cellVec, iCellVec, denseDesc, iSparseStart, img2CentCell, orb, over,&
              & kPoint, kWeight, iOrbRegion, parallelKS, eigvecsReal, SSqrReal, eigvecsCplx,&
              & SSqrCplx)
      end if
      call env%globalTimer%stopTimer(globalTimers%eigvecWriting)

      ! MD geometry files are written only later, once velocities for the current geometry are known
      if (tGeoOpt .and. tWriteRestart) then
        call writeCurrentGeometry(geoOutFile, pCoord0Out, tLatOpt, tMd, tAppendGeo, tFracCoord,&
            & tPeriodic, tPrintMulliken, species0, speciesName, latVec, iGeoStep, iLatGeoStep,&
            & nSpin, qOutput, velocities)
      end if

      call printEnergies(energy)

      if (tForces) then
        call env%globalTimer%startTimer(globalTimers%forceCalc)
        call env%globalTimer%startTimer(globalTimers%energyDensityMatrix)
        call getEnergyWeightedDensity(env, denseDesc, forceType, filling, eigen, kPoint, kWeight,&
            & neighborList, nNeighbor, orb, iSparseStart, img2CentCell, iCellVec, cellVec,&
            & tRealHS, ham, over, parallelKS, ERhoPrim, eigvecsReal, SSqrReal, eigvecsCplx,&
            & SSqrCplx, rhoSqrReal)
        call env%globalTimer%stopTimer(globalTimers%energyDensityMatrix)
        call getGradients(env, sccCalc, tEField, tXlbomd, nonSccDeriv, Efield, rhoPrim, ERhoPrim,&
            & qOutput, q0, skHamCont, skOverCont, pRepCont, neighborList, nNeighbor, species,&
            & img2CentCell, iSparseStart, orb, potential, coord, derivs, iRhoPrim, thirdOrd,&
            & chrgForces, dispersion)
        if (tLinResp) then
          derivs(:,:) = derivs + excitedDerivs
        end if
        call env%globalTimer%stopTimer(globalTimers%forceCalc)

        if (tStress) then
          call env%globalTimer%startTimer(globalTimers%stressCalc)
          call getStress(env, sccCalc, tEField, nonSccDeriv, EField, rhoPrim, ERhoPrim, qOutput,&
              & q0, skHamCont, skOverCont, pRepCont, neighborList, nNeighbor, species,&
              & img2CentCell, iSparseStart, orb, potential, coord, latVec, invLatVec, cellVol,&
              & coord0, totalStress, totalLatDeriv, intPressure, iRhoPrim, dispersion)
          call env%globalTimer%stopTimer(globalTimers%stressCalc)
          call printVolume(cellVol)
          ! MD case includes the atomic kinetic energy contribution, so print that later
          if (.not. tMD) then
            call printPressureAndFreeEnergy(extPressure, intPressure, energy%EGibbs)
          end if
        end if

      end if

      if (tWriteDetailedOut) then
        call writeDetailedOut2(fdDetailedOut, tSccCalc, tConverged, tXlbomd, tLinResp, tGeoOpt,&
            & tMD, tPrintForces, tStress, tPeriodic, energy, totalStress, totalLatDeriv, derivs, &
            & chrgForces, indMovedAtom, cellVol, intPressure, geoOutFile)
      end if

      if (tSccCalc .and. .not. tXlbomd .and. .not. tConverged) then
        call warning("SCC is NOT converged, maximal SCC iterations exceeded")
        if (tUseConvergedForces) then
          call env%shutdown()
        end if
      end if

      if (tSccCalc .and. allocated(esp) .and. (.not. (tGeoOpt .or. tMD) .or. &
          & needsRestartWriting(tGeoOpt, tMd, iGeoStep, nGeoSteps, restartFreq))) then
        call esp%evaluate(env, sccCalc, EField)
        call writeEsp(esp, env, iGeoStep, nGeoSteps)
      end if

      if (tForces) then
        if (allocated(conAtom)) then
          call constrainForces(conAtom, conVec, derivs)
        end if

        if (tCoordOpt) then
          call printMaxForce(maxval(abs(derivs(:, indMovedAtom))))
        end if

        if (tLatOpt) then
          ! Only include the extLatDerivs contribution if not MD, as the barostat would otherwise
          ! take care of this, hence add it here rather than to totalLatDeriv itself
          call constrainLatticeDerivs(totalLatDeriv + extLatDerivs, normOrigLatVec,&
              & tLatOptFixAng, tLatOptFixLen, tLatOptIsotropic, constrLatDerivs)
          call printMaxLatticeForce(maxval(abs(constrLatDerivs)))
        end if

        if (tSocket .and. env%tGlobalMaster) then
          ! stress was computed above in the force evaluation block or is 0 if aperiodic
        #:if WITH_SOCKETS
          call socket%send(energy%ETotal - sum(TS), -derivs, totalStress * cellVol)
        #:else
          call error("Should not be here - compiled without socket support")
        #:endif
        end if

        ! If geometry minimizer finished and the last calculated geometry is the minimal one (not
        ! necessarily the case, depends on the optimizer!) we are finished.  Otherwise we have to
        ! recalculate everything at the converged geometry.

        if (tGeomEnd) then
          call env%globalTimer%stopTimer(globalTimers%postSCC)
          exit lpGeomOpt
        end if

        tWriteCharges = tWriteRestart .and. tMulliken .and. tSccCalc .and. .not. tDerivs&
            & .and. maxSccIter > 1
        if (tWriteCharges) then
          call writeCharges(fCharges, fdCharges, tWriteChrgAscii, orb, qInput, qBlockIn, qiBlockIn)
        end if

        ! initially assume coordinates are not being updated
        tCoordsChanged = .false.
        tLatticeChanged = .false.

        if (tDerivs) then
          call getNextDerivStep(derivDriver, derivs, indMovedAtom, coord0, tGeomEnd)
          if (tGeomEnd) then
            call env%globalTimer%stopTimer(globalTimers%postSCC)
            exit lpGeomOpt
          end if
          tCoordsChanged = .true.
        else if (tGeoOpt) then
          tCoordsChanged = .true.
          if (tCoordStep) then
            call getNextCoordinateOptStep(pGeoCoordOpt, energy%EMermin, derivs, indMovedAtom,&
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
            call env%globalTimer%stopTimer(globalTimers%postSCC)
            exit lpGeomOpt
          end if
        else if (tMD) then
          ! New MD coordinates saved in a temporary variable, as writeCurrentGeometry() below
          ! needs the old ones to write out consistent geometries and velocities.
          newCoords(:,:) = coord0
          call getNextMdStep(pMdIntegrator, pMdFrame, temperatureProfile, derivs, movedMass,&
              & mass, cellVol, invLatVec, species0, indMovedAtom, tStress, tBarostat, energy,&
              & newCoords, latVec, intPressure, totalStress, totalLatDeriv, velocities, tempIon)
          tCoordsChanged = .true.
          tLatticeChanged = tBarostat
          call printMdInfo(tSetFillingTemp, tEField, tPeriodic, tempElec, absEField, tempIon,&
              & intPressure, extPressure, energy)
          if (tWriteRestart) then
            if (tPeriodic) then
              cellVol = abs(determinant33(latVec))
              energy%EGibbs = energy%EMermin + extPressure * cellVol
            end if
            call writeMdOut2(fdMd, tStress, tBarostat, tLinResp, tEField, tFixEf, tPrintMulliken,&
                & energy, latVec, cellVol, intPressure, extPressure, tempIon, absEField, qOutput,&
                & q0, dipoleMoment)
            call writeCurrentGeometry(geoOutFile, pCoord0Out, .false., .true., .true., tFracCoord,&
                & tPeriodic, tPrintMulliken, species0, speciesName, latVec, iGeoStep, iLatGeoStep,&
                & nSpin, qOutput, velocities)
          end if
          coord0(:,:) = newCoords
          if (tWriteDetailedOut) then
            call writeDetailedOut3(fdDetailedOut, tPrintForces, tSetFillingTemp, tPeriodic,&
                & tStress, totalStress, totalLatDeriv, energy, tempElec, extPressure, intPressure,&
                & tempIon)
          end if
        else if (tSocket .and. iGeoStep < nGeoSteps) then
          ! Only receive geometry from socket, if there are still geometry iterations left
        #:if WITH_SOCKETS
          call receiveGeometryFromSocket(env, socket, tPeriodic, coord0, latVec, tCoordsChanged,&
              & tLatticeChanged, tStopDriver)
        #:else
          call error("Internal error: code compiled without socket support")
        #:endif
        end if
      end if

      if (tWriteDetailedOut) then
        call writeDetailedOut4(fdDetailedOut, tMD, energy, tempIon)
      end if

      tStopDriver = tStopScc .or. tStopDriver .or. hasStopFile(fStopDriver)
      if (tStopDriver) then
        call env%globalTimer%stopTimer(globalTimers%postSCC)
        exit lpGeomOpt
      end if

      call env%globalTimer%stopTimer(globalTimers%postSCC)

    end do lpGeomOpt

    call env%globalTimer%startTimer(globalTimers%postGeoOpt)

  #:if WITH_SOCKETS
    if (tSocket .and. env%tGlobalMaster) then
      call socket%shutdown()
    end if
  #:endif

    tGeomEnd = tMD .or. tGeomEnd .or. tDerivs

    if (env%tGlobalMaster) then
      if (tWriteDetailedOut) then
        call writeDetailedOut5(fdDetailedOut, tGeoOpt, tGeomEnd, tMd, tDerivs, tEField, absEField,&
            & dipoleMoment)
      end if

      call writeFinalDriverStatus(tGeoOpt, tGeomEnd, tMd, tDerivs)

      if (tMD) then
        call writeMdOut3(fdMd, mdOut)
      end if
    end if

    if (env%tGlobalMaster .and. tDerivs) then
      call getHessianMatrix(derivDriver, pDynMatrix)
      call writeHessianOut(fdHessian, hessianOut, pDynMatrix)
    else
      nullify(pDynMatrix)
    end if


    if (allocated(pipekMezey)) then
      ! NOTE: the canonical DFTB ground state orbitals are over-written after this point
      if (withMpi) then
        call error("Pipek-Mezey localisation does not yet work with MPI")
      end if
      if (nSpin > 2) then
        call error("Pipek-Mezey localisation not implemented for non-colinear DFTB")
      end if
      call calcPipekMezeyLocalisation(env, pipekMezey, tPrintEigvecsTxt, nEl, filling, over,&
          & kPoint, neighborList, nNeighbor, denseDesc, iSparseStart, img2CentCell, iCellVec,&
          & cellVec, fdEigvec, runId, orb, species, speciesName, parallelKS, localisation,&
          & eigvecsReal, SSqrReal, eigvecsCplx, SSqrCplx)
    end if

    if (tWriteAutotest) then
      if (tPeriodic) then
        cellVol = abs(determinant33(latVec))
        energy%EGibbs = energy%EMermin + extPressure * cellVol
      end if
      call writeAutotestTag(fdAutotest, autotestTag, tPeriodic, cellVol, tMulliken, qOutput,&
          & derivs, chrgForces, excitedDerivs, tStress, totalStress, pDynMatrix,&
          & energy%EMermin, extPressure, energy%EGibbs, coord0, tLocalise, localisation, esp)
    end if
    if (tWriteResultsTag) then
      call writeResultsTag(fdResultsTag, resultsTag, derivs, chrgForces, tStress, totalStress,&
          & pDynMatrix, tPeriodic, cellVol)
    end if
    if (tWriteDetailedXML) then
      call writeDetailedXml(runId, speciesName, species0, pCoord0Out, tPeriodic, latVec, tRealHS,&
          & nKPoint, nSpin, size(eigen, dim=1), nOrb, kPoint, kWeight, filling, occNatural)
    end if

    call env%globalTimer%startTimer(globalTimers%postGeoOpt)

    call destructProgramVariables()

  end subroutine runDftbPlus


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


  !> Initialises SCC related parameters before geometry loop starts
  function getMinSccIters(tSccCalc, tDftbU, nSpin) result(minSccIter)

    !> Is this a self consistent calculation
    logical, intent(in) :: tSccCalc

    !> Are there orbital potentials present
    logical, intent(in) :: tDftbU

    !> Number of spin channels
    integer, intent(in) :: nSpin

    !> Minimum possible number of self consistent iterations
    integer :: minSccIter

    if (tSccCalc) then
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
  subroutine handleLatticeChange(latVecs, sccCalc, tStress, extPressure, mCutoff,&
      & dispersion, recVecs, recVecs2p, cellVol, recCellVol, extLatDerivs, cellVecs, rCellVecs)

    !> lattice vectors
    real(dp), intent(in) :: latVecs(:,:)

    !> Module variables
    type(TScc), allocatable, intent(inout) :: sccCalc

    !> evaluate stress
    logical, intent(in) :: tStress

    !> External presure
    real(dp), intent(in) :: extPressure

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
    real(dp), intent(out) :: extLatDerivs(:,:)

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
      call derivDeterminant33(extLatDerivs, latVecs)
      extLatDerivs(:,:) = extPressure * extLatDerivs
    end if
    if (allocated(sccCalc)) then
      call sccCalc%updateLatVecs(latVecs, recVecs, cellVol)
      mCutoff = max(mCutoff, sccCalc%getCutoff())
    end if
    if (allocated(dispersion)) then
      call dispersion%updateLatVecs(latVecs)
      mCutoff = max(mCutoff, dispersion%getRCutoff())
    end if
    call getCellTranslations(cellVecs, rCellVecs, latVecs, recVecs2p, mCutoff)

  end subroutine handleLatticeChange


  !> Does the operations that are necessary after atomic coordinates change
  subroutine handleCoordinateChange(env, coord0, latVec, invLatVec, species0, mCutoff,&
      & skRepCutoff, orb, tPeriodic, sccCalc, dispersion, thirdOrd, img2CentCell, iCellVec,&
      & neighborList, nAllAtom, coord0Fold, coord, species, rCellVec, nAllOrb, nNeighbor, ham,&
      & over, H0, rhoPrim, iRhoPrim, iHam, ERhoPrim, iSparseStart)

    !> Environment settings
    type(TEnvironment), intent(in) :: env

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

    !> SCC module internal variables
    type(TScc), allocatable, intent(inout) :: sccCalc

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
    integer, allocatable, intent(inout) :: iSparseStart(:,:)

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
    call getSparseDescriptor(neighborList%iNeighbor, nNeighbor, img2CentCell, orb, iSparseStart,&
        & sparseSize)
    call reallocateSparseArrays(sparseSize, ham, over, H0, rhoPrim, iHam, iRhoPrim, ERhoPrim)

    ! Notify various modules about coordinate changes
    if (allocated(sccCalc)) then
      call sccCalc%updateCoords(env, coord, species, neighborList)
    end if
    if (allocated(dispersion)) then
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
        @:ASSERT(size(ERhoPrim) == size(ham, dim=1))
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
  subroutine setUpExternalElectricField(tEfield, tTimeDepEField, tPeriodic, EFieldStrength,&
      & EFieldVector, EFieldOmega, EFieldPhase, neighborList, nNeighbor, iCellVec, img2CentCell,&
      & cellVec, deltaT, iGeoStep, coord0Fold, coord, EField, extAtomPot, absEField)

    !> Whether electric field should be considered at all
    logical, intent(in) :: tEfield

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

    !> Vectors (in units of the lattice constants) to cells of the lattice
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

    if (.not. tEField) then
      EField(:) = 0.0_dp
      absEField = 0.0_dp
      extAtomPot(:) = 0.0_dp
      return
    end if

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
  subroutine initSccLoop(tSccCalc, xlbomdIntegrator, minSccIter, maxSccIter, sccTol, tConverged)

    !> Is this an SCC calculation?
    logical, intent(in) :: tSccCalc

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

    tConverged = (.not. tSccCalc)

    if (tSccCalc) then
      call printSccHeader()
    end if

  end subroutine initSccLoop


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

    @:ASSERT(.not. tDualSpinOrbit .or. allocated(xi))

    potential%intAtom(:,:) = 0.0_dp
    potential%intShell(:,:,:) = 0.0_dp
    potential%intBlock(:,:,:,:) = 0.0_dp
    potential%orbitalBlock(:,:,:,:) = 0.0_dp
    potential%iOrbitalBlock(:,:,:,:) = 0.0_dp
    if (tDualSpinOrbit) then
      call getDualSpinOrbitShift(potential%iOrbitalBlock, xi, orb, species)
    end if

  end subroutine resetInternalPotentials


  !> Add potentials comming from point charges.
  subroutine addChargePotentials(env, sccCalc, qInput, q0, chargePerShell, orb, species,&
      & neighborList, img2CentCell, spinW, thirdOrd, potential)

    !> Environment settings
    type(TEnvironment), intent(in) :: env

    !> SCC module internal variables
    type(TScc), intent(inout) :: sccCalc

    !> Input atomic populations
    real(dp), intent(in) :: qInput(:,:,:)

    !> reference atomic occupations
    real(dp), intent(in) :: q0(:,:,:)

    !> charges per atomic shell
    real(dp), intent(in) :: chargePerShell(:,:,:)

    !> atomic orbital information
    type(TOrbitals), intent(in) :: orb

    !> species of all atoms
    integer, target, intent(in) :: species(:)

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
    integer, pointer :: pSpecies0(:)
    integer :: nAtom, nSpin

    nAtom = size(qInput, dim=2)
    nSpin = size(qInput, dim=3)
    pSpecies0 => species(1:nAtom)

    allocate(atomPot(nAtom, nSpin))
    allocate(shellPot(orb%mShell, nAtom, nSpin))

    call sccCalc%updateCharges(env, qInput, q0, orb, species, neighborList%iNeighbor, img2CentCell)
    call sccCalc%getShiftPerAtom(atomPot(:,1))
    call sccCalc%getShiftPerL(shellPot(:,:,1))
    potential%intAtom(:,1) = potential%intAtom(:,1) + atomPot(:,1)
    potential%intShell(:,:,1) = potential%intShell(:,:,1) + shellPot(:,:,1)

    if (allocated(thirdOrd)) then
      call thirdOrd%updateCharges(pSpecies0, neighborList, qInput, q0, img2CentCell, orb)
      call thirdOrd%getShifts(atomPot(:,1), shellPot(:,:,1))
      potential%intAtom(:,1) = potential%intAtom(:,1) + atomPot(:,1)
      potential%intShell(:,:,1) = potential%intShell(:,:,1) + shellPot(:,:,1)
    end if

    if (nSpin /= 1 .and. allocated(spinW)) then
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
    integer, intent(in), allocatable :: nUJ(:)

    !> which shells are in each DFTB+U block
    integer, intent(in), allocatable :: iUJ(:,:,:)

    !> Number of shells in each DFTB+U block
    integer, intent(in), allocatable :: niUJ(:,:)

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
  subroutine getSccHamiltonian(H0, over, nNeighbor, neighborList, species, orb, iSparseStart,&
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
    integer, intent(in) :: iSparseStart(:,:)

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
    call add_shift(ham, over, nNeighbor, neighborList%iNeighbor, species, orb, iSparseStart, nAtom,&
        & img2CentCell, potential%intBlock)

    if (allocated(iHam)) then
      iHam(:,:) = 0.0_dp
      call add_shift(iHam, over, nNeighbor, neighborList%iNeighbor, species, orb, iSparseStart,&
          & nAtom, img2CentCell, potential%iorbitalBlock)
    end if

  end subroutine getSccHamiltonian


  !> Returns the sparse density matrix.
  !>
  !> All operations (e.g. non-dual spin orbit coupling), which need access to full (unpacked)
  !> Hamiltonian or the full (unpacked) density matrix, must also invoked from within this routine,
  !> as those unpacked quantities do not exist elsewhere.
  !>
  subroutine getDensity(env, denseDesc, ham, over, neighborList, nNeighbor, iSparseStart,&
      & img2CentCell, iCellVec, cellVec, kPoint, kWeight, orb, species, solver, tRealHS,&
      & tSpinSharedEf, tSpinOrbit, tDualSpinOrbit, tFillKSep, tFixEf, tMulliken, iDistribFn,&
      & tempElec, nEl, parallelKS, Ef, energy, eigen, filling, rhoPrim, Eband, TS, E0, iHam, xi,&
      & orbitalL, HSqrReal, SSqrReal, eigvecsReal, iRhoPrim, HSqrCplx, SSqrCplx, eigvecsCplx,&
      & rhoSqrReal)

    !> Environment settings
    type(TEnvironment), intent(inout) :: env

    !> Dense matrix descriptor
    type(TDenseDescr), intent(in) :: denseDesc

    !> hamiltonian in sparse storage
    real(dp), intent(inout) :: ham(:,:)

    !> sparse overlap matrix
    real(dp), intent(in) :: over(:)

    !> list of neighbours for each atom
    type(TNeighborList), intent(in) :: neighborList

    !> Number of neighbours for each of the atoms
    integer, intent(in) :: nNeighbor(:)

    !> Index array for the start of atomic blocks in sparse arrays
    integer, intent(in) :: iSparseStart(:,:)

    !> map from image atoms to the original unique atom
    integer, intent(in) :: img2CentCell(:)

    !> Index for which unit cell atoms are associated with
    integer, intent(in) :: iCellVec(:)

    !> Vectors (in units of the lattice constants) to cells of the lattice
    real(dp), intent(in) :: cellVec(:,:)

    !> k-points
    real(dp), intent(in) :: kPoint(:,:)

    !> Weights for k-points
    real(dp), intent(in) :: kWeight(:)

    !> Atomic orbital information
    type(TOrbitals), intent(in) :: orb

    !> species of all atoms in the system
    integer, intent(in) :: species(:)

    !> Eigensolver choice
    type(TSolver), intent(inout) :: solver

    !> Is the hamitonian real (no k-points/molecule/gamma point)?
    logical, intent(in) :: tRealHS

    !> Is the Fermi level common accross spin channels?
    logical, intent(in) :: tSpinSharedEf

    !> Are spin orbit interactions present
    logical, intent(in) :: tSpinOrbit

    !> Are block population spin orbit interactions present
    logical, intent(in) :: tDualSpinOrbit

    !> Fill k-points separately if true (no charge transfer accross the BZ)
    logical, intent(in) :: tFillKSep

    !> Whether fixed Fermi level(s) should be used. (No charge conservation!)
    logical, intent(in) :: tFixEf

    !> Should Mulliken populations be generated/output
    logical, intent(in) :: tMulliken

    !> occupation function for electronic states
    integer, intent(in) :: iDistribFn

    !> Electronic temperature
    real(dp), intent(in) :: tempElec

    !> Number of electrons
    real(dp), intent(in) :: nEl(:)

    !> K-points and spins to process
    type(TParallelKS), intent(in) :: parallelKS

    !> Fermi level(s)
    real(dp), intent(inout) :: Ef(:)

    !> Energy contributions and total
    type(TEnergies), intent(inout) :: energy

    !> eigenvalues
    real(dp), intent(out) :: eigen(:,:,:)

    !> occupations
    real(dp), intent(out) :: filling(:,:,:)

    !> sparse density matrix
    real(dp), intent(out) :: rhoPrim(:,:)

    !> band structure energy
    real(dp), intent(out) :: Eband(:)

    !> electronic entropy times temperature
    real(dp), intent(out) :: TS(:)

    !> extrapolated 0 temperature band energy
    real(dp), intent(out) :: E0(:)

    !> imaginary part of hamitonian
    real(dp), intent(inout), allocatable :: iHam(:,:)

    !> spin orbit constants
    real(dp), intent(in), allocatable :: xi(:,:)

    !> orbital moments of atomic shells
    real(dp), intent(inout), allocatable :: orbitalL(:,:,:)

    !> imaginary part of density matrix
    real(dp), intent(inout), allocatable :: iRhoPrim(:,:)

    !> dense real hamiltonian storage
    real(dp), intent(inout), allocatable :: HSqrReal(:,:)

    !> dense real overlap storage
    real(dp), intent(inout), allocatable :: SSqrReal(:,:)

    !> real eigenvectors on exit
    real(dp), intent(inout), allocatable :: eigvecsReal(:,:,:)

    !> dense complex (k-points) hamiltonian storage
    complex(dp), intent(inout), allocatable :: HSqrCplx(:,:)

    !> dense complex (k-points) overlap storage
    complex(dp), intent(inout), allocatable :: SSqrCplx(:,:)

    !> complex eigenvectors on exit
    complex(dp), intent(inout), allocatable :: eigvecsCplx(:,:,:)

    !> Dense density matrix
    real(dp), intent(inout), allocatable :: rhoSqrReal(:,:,:)

    integer :: nSpin

    nSpin = size(ham, dim=2)

    ! Hack due to not using Pauli-type structure for diagonalisation
    if (nSpin > 1) then
      ham(:,:) = 2.0_dp * ham
      if (allocated(iHam)) then
        iHam(:,:) = 2.0_dp * iHam
      end if
    end if

    call qm2ud(ham)

    select case (solver%solverType)

    case (solverTypes%lapackQr, solverTypes%lapackDivAndConq, solverTypes%lapackRelRobust)
      call getDensityByDenseDiag(env, denseDesc, ham, over, neighborList, nNeighbor,&
          & iSparseStart, img2CentCell, iCellVec, cellVec, kPoint, kWeight, orb, species, solver,&
          & tRealHS, tSpinSharedEf, tSpinOrbit, tDualSpinOrbit, tFillKSep, tFixEf, tMulliken,&
          & iDistribFn, tempElec, nEl, parallelKS, Ef, energy, eigen, filling, rhoPrim, Eband, TS,&
          & E0, iHam, xi, orbitalL, HSqrReal, SSqrReal, eigvecsReal, iRhoPrim, HSqrCplx, SSqrCplx,&
          & eigvecsCplx, rhoSqrReal)

    case (solverTypes%progressSp2)
      if (.not. tRealHS .or. nSpin == 4) then
        call error("Internal error: SP2 solver does not work with complex Hamiltonian or&
            & non-colinear spin")
      end if
      call getDensityBySp2(env, denseDesc, ham, neighborList, nNeighbor, iSparseStart,&
          & img2CentCell, orb, solver, nEl, parallelKS, rhoPrim, rhoSqrReal)

    case default
      call error("Internal error: invalid solver type in getDensity")

    end select

    call ud2qm(rhoPrim)

  end subroutine getDensity


  subroutine getDensityByDenseDiag(env, denseDesc, ham, over, neighborList, nNeighbor,&
      & iSparseStart, img2CentCell, iCellVec, cellVec, kPoint, kWeight, orb, species, solver,&
      & tRealHS, tSpinSharedEf, tSpinOrbit, tDualSpinOrbit, tFillKSep, tFixEf, tMulliken,&
      & iDistribFn, tempElec, nEl, parallelKS, Ef, energy, eigen, filling, rhoPrim, Eband, TS, E0,&
      & iHam, xi, orbitalL, HSqrReal, SSqrReal, eigvecsReal, iRhoPrim, HSqrCplx, SSqrCplx,&
      & eigvecsCplx, rhoSqrReal)

    !> Environment settings
    type(TEnvironment), intent(inout) :: env

    !> Dense matrix descriptor
    type(TDenseDescr), intent(in) :: denseDesc

    !> hamiltonian in sparse storage
    real(dp), intent(in) :: ham(:,:)

    !> sparse overlap matrix
    real(dp), intent(in) :: over(:)

    !> list of neighbours for each atom
    type(TNeighborList), intent(in) :: neighborList

    !> Number of neighbours for each of the atoms
    integer, intent(in) :: nNeighbor(:)

    !> Index array for the start of atomic blocks in sparse arrays
    integer, intent(in) :: iSparseStart(:,:)

    !> map from image atoms to the original unique atom
    integer, intent(in) :: img2CentCell(:)

    !> Index for which unit cell atoms are associated with
    integer, intent(in) :: iCellVec(:)

    !> Vectors (in units of the lattice constants) to cells of the lattice
    real(dp), intent(in) :: cellVec(:,:)

    !> k-points
    real(dp), intent(in) :: kPoint(:,:)

    !> Weights for k-points
    real(dp), intent(in) :: kWeight(:)

    !> Atomic orbital information
    type(TOrbitals), intent(in) :: orb

    !> species of all atoms in the system
    integer, intent(in) :: species(:)

    !> Eigensolver choice
    type(TSolver), intent(in) :: solver

    !> Is the hamitonian real (no k-points/molecule/gamma point)?
    logical, intent(in) :: tRealHS

    !> Is the Fermi level common accross spin channels?
    logical, intent(in) :: tSpinSharedEf

    !> Are spin orbit interactions present
    logical, intent(in) :: tSpinOrbit

    !> Are block population spin orbit interactions present
    logical, intent(in) :: tDualSpinOrbit

    !> Fill k-points separately if true (no charge transfer accross the BZ)
    logical, intent(in) :: tFillKSep

    !> Whether fixed Fermi level(s) should be used. (No charge conservation!)
    logical, intent(in) :: tFixEf

    !> Should Mulliken populations be generated/output
    logical, intent(in) :: tMulliken

    !> occupation function for electronic states
    integer, intent(in) :: iDistribFn

    !> Electronic temperature
    real(dp), intent(in) :: tempElec

    !> Number of electrons
    real(dp), intent(in) :: nEl(:)

    !> K-points and spins to process
    type(TParallelKS), intent(in) :: parallelKS

    !> Fermi level(s)
    real(dp), intent(inout) :: Ef(:)

    !> Energy contributions and total
    type(TEnergies), intent(inout) :: energy

    !> eigenvalues
    real(dp), intent(out) :: eigen(:,:,:)

    !> occupations
    real(dp), intent(out) :: filling(:,:,:)

    !> sparse density matrix
    real(dp), intent(out) :: rhoPrim(:,:)

    !> band structure energy
    real(dp), intent(out) :: Eband(:)

    !> electronic entropy times temperature
    real(dp), intent(out) :: TS(:)

    !> extrapolated 0 temperature band energy
    real(dp), intent(out) :: E0(:)

    !> imaginary part of hamitonian
    real(dp), intent(in), allocatable :: iHam(:,:)

    !> spin orbit constants
    real(dp), intent(in), allocatable :: xi(:,:)

    !> orbital moments of atomic shells
    real(dp), intent(inout), allocatable :: orbitalL(:,:,:)

    !> imaginary part of density matrix
    real(dp), intent(inout), allocatable :: iRhoPrim(:,:)

    !> dense real hamiltonian storage
    real(dp), intent(inout), allocatable :: HSqrReal(:,:)

    !> dense real overlap storage
    real(dp), intent(inout), allocatable :: SSqrReal(:,:)

    !> real eigenvectors on exit
    real(dp), intent(inout), allocatable :: eigvecsReal(:,:,:)

    !> dense complex (k-points) hamiltonian storage
    complex(dp), intent(inout), allocatable :: HSqrCplx(:,:)

    !> dense complex (k-points) overlap storage
    complex(dp), intent(inout), allocatable :: SSqrCplx(:,:)

    !> complex eigenvectors on exit
    complex(dp), intent(inout), allocatable :: eigvecsCplx(:,:,:)

    !> Dense density matrix
    real(dp), intent(inout), allocatable :: rhoSqrReal(:,:,:)

    integer :: nSpin

    nSpin = size(ham, dim=2)

    call env%globalTimer%startTimer(globalTimers%diagonalization)

    if (nSpin /= 4) then
      if (tRealHS) then
        call buildAndDiagDenseRealHam(env, denseDesc, ham, over, neighborList, nNeighbor,&
            & iSparseStart, img2CentCell, solver, parallelKS, HSqrReal, SSqrReal, eigVecsReal,&
            & eigen(:,1,:))
      else
        call buildAndDiagDenseCplxHam(env, denseDesc, ham, over, kPoint, neighborList, nNeighbor,&
            & iSparseStart, img2CentCell, iCellVec, cellVec, solver, parallelKS, HSqrCplx,&
            & SSqrCplx, eigVecsCplx, eigen)
      end if
    else
      call buildAndDiagDensePauliHam(env, denseDesc, ham, over, kPoint, neighborList, nNeighbor,&
          & iSparseStart, img2CentCell, iCellVec, cellVec, orb, solver, parallelKS, eigen(:,:,1),&
          & HSqrCplx, SSqrCplx, eigVecsCplx, iHam, xi, species)
    end if

    call env%globalTimer%stopTimer(globalTimers%diagonalization)
    
    call getFillingsAndBandEnergies(eigen, nEl, nSpin, tempElec, kWeight, tSpinSharedEf,&
        & tFillKSep, tFixEf, iDistribFn, Ef, filling, Eband, TS, E0)

    call env%globalTimer%startTimer(globalTimers%densityMatrix)
    if (nSpin /= 4) then
      if (tRealHS) then
        call getDensityFromRealEigvecs(env, denseDesc, filling(:,1,:), neighborList, nNeighbor,&
            & iSparseStart, img2CentCell, orb, eigVecsReal, parallelKS, rhoPrim, SSqrReal,&
            & rhoSqrReal)
      else
        call getDensityFromCplxEigvecs(env, denseDesc, filling, kPoint, kWeight, neighborList,&
            & nNeighbor, iSparseStart, img2CentCell, iCellVec, cellVec, orb, parallelKS,&
            & eigvecsCplx, rhoPrim, SSqrCplx)
      end if
    else
      ! Pauli structure of eigenvectors
      filling(:,:,1) = 2.0_dp * filling(:,:,1)
      call getDensityFromPauliEigvecs(env, denseDesc, tRealHS, tSpinOrbit, tDualSpinOrbit,&
          & tMulliken, kPoint, kWeight, filling(:,:,1), neighborList, nNeighbor, orb, iSparseStart,&
          & img2CentCell, iCellVec, cellVec, species, parallelKS, eigVecsCplx, SSqrCplx, energy,&
          & rhoPrim, xi, orbitalL, iRhoPrim)
      filling(:,:,1) = 0.5_dp * filling(:,:,1)
    end if

    call env%globalTimer%stopTimer(globalTimers%densityMatrix)

  end subroutine getDensityByDenseDiag


  !> Obtains the density by using the SP2 method
  subroutine getDensityBySp2(env, denseDesc, ham, neighborList, nNeighbor, iSparseStart,&
      & img2CentCell, orb, solver, nEl, parallelKS, rhoPrim, rhoSqrReal)

    !> Environment settings
    type(TEnvironment), intent(inout) :: env

    !> Dense matrix descriptor
    type(TDenseDescr), intent(in) :: denseDesc

    !> hamiltonian in sparse storage
    real(dp), intent(in) :: ham(:,:)

    !> list of neighbours for each atom
    type(TNeighborList), intent(in) :: neighborList

    !> Number of neighbours for each of the atoms
    integer, intent(in) :: nNeighbor(:)

    !> Index array for the start of atomic blocks in sparse arrays
    integer, intent(in) :: iSparseStart(:,:)

    !> map from image atoms to the original unique atom
    integer, intent(in) :: img2CentCell(:)

    !> Atomic orbital information
    type(TOrbitals), intent(in) :: orb

    !> Eigensolver choice
    type(TSolver), intent(inout) :: solver

    !> Number of electrons
    real(dp), intent(in) :: nEl(:)

    !> K-points and spins to process
    type(TParallelKS), intent(in) :: parallelKS

    !> sparse density matrix
    real(dp), intent(out) :: rhoPrim(:,:)

    !> Dense density matrix to be used for later purposes
    real(dp), allocatable, intent(inout) :: rhoSqrReal(:,:,:)

    call env%globalTimer%startTimer(globalTimers%densityMatrix)
    
  #:if WITH_PROGRESS
    call solver%sp2Solver%getDensity(ham, neighborList, nNeighbor, iSparseStart,&
        & img2CentCell, denseDesc, orb, nEl, parallelKS, rhoPrim, rhoSqrReal)
  #:else
    call error("Internal error: getDensityBySp2 called despite missing Progress support")
  #:endif
      
    call env%globalTimer%stopTimer(globalTimers%densityMatrix)

  end subroutine getDensityBySp2


  !> Builds and diagonalises dense Hamiltonians.
  subroutine buildAndDiagDenseRealHam(env, denseDesc, ham, over, neighborList, nNeighbor,&
      & iSparseStart, img2CentCell, solver, parallelKS, HSqrReal, SSqrReal, eigvecsReal, eigen)

    !> Environment settings
    type(TEnvironment), intent(inout) :: env

    !> Dense matrix descriptor
    type(TDenseDescr), intent(in) :: denseDesc

    !> hamiltonian in sparse storage
    real(dp), intent(in) :: ham(:,:)

    !> sparse overlap matrix
    real(dp), intent(in) :: over(:)

    !> list of neighbours for each atom
    type(TNeighborList), intent(in) :: neighborList

    !> Number of neighbours for each of the atoms
    integer, intent(in) :: nNeighbor(:)

    !> Index array for the start of atomic blocks in sparse arrays
    integer, intent(in) :: iSparseStart(:,:)

    !> map from image atoms to the original unique atom
    integer, intent(in) :: img2CentCell(:)

    !> Eigensolver choice
    type(TSolver), intent(in) :: solver

    !> K-points and spins to be handled
    type(TParallelKS), intent(in) :: parallelKS

    !> dense hamitonian matrix
    real(dp), intent(out) :: HSqrReal(:,:)

    !> dense overlap matrix
    real(dp), intent(out) :: SSqrReal(:,:)

    !> Eigenvectors on eixt
    real(dp), intent(out) :: eigvecsReal(:,:,:)

    !> eigenvalues
    real(dp), intent(out) :: eigen(:,:)

    integer :: iKS, iSpin

    eigen(:,:) = 0.0_dp
    do iKS = 1, parallelKS%nLocalKS
      iSpin = parallelKS%localKS(2, iKS)
    #:if WITH_SCALAPACK
      call env%globalTimer%startTimer(globalTimers%sparseToDense)
      call unpackHSRealBlacs(env%blacs, ham(:,iSpin), neighborList%iNeighbor, nNeighbor,&
          & iSparseStart, img2CentCell, denseDesc, HSqrReal)
      call unpackHSRealBlacs(env%blacs, over, neighborList%iNeighbor, nNeighbor, iSparseStart,&
          & img2CentCell, denseDesc, SSqrReal)
      call env%globalTimer%stopTimer(globalTimers%sparseToDense)
      call diagDenseMtxBlacs(solver, 'V', denseDesc%blacsOrbSqr, HSqrReal, SSqrReal,&
          & eigen(:,iSpin), eigvecsReal(:,:,iKS))
    #:else
      call env%globalTimer%startTimer(globalTimers%sparseToDense)
      call unpackHS(HSqrReal, ham(:,iSpin), neighborList%iNeighbor, nNeighbor,&
          & denseDesc%iAtomStart, iSparseStart, img2CentCell)
      call unpackHS(SSqrReal, over, neighborList%iNeighbor, nNeighbor, denseDesc%iAtomStart,&
          & iSparseStart, img2CentCell)
      call env%globalTimer%stopTimer(globalTimers%sparseToDense)
      call diagDenseMtx(solver, 'V', HSqrReal, SSqrReal, eigen(:,iSpin))
      eigvecsReal(:,:,iKS) = HSqrReal
    #:endif
    end do

  #:if WITH_SCALAPACK
    ! Distribute all eigenvalues to all nodes via global summation
    call mpifx_allreduceip(env%mpi%interGroupComm, eigen, MPI_SUM)
  #:endif

  end subroutine buildAndDiagDenseRealHam


  !> Builds and diagonalises dense k-point dependent Hamiltonians.
  subroutine buildAndDiagDenseCplxHam(env, denseDesc, ham, over, kPoint, neighborList, nNeighbor,&
      & iSparseStart, img2CentCell, iCellVec, cellVec, solver, parallelKS, HSqrCplx, SSqrCplx,&
      & eigvecsCplx, eigen)

    !> Environment settings
    type(TEnvironment), intent(inout) :: env

    !> Dense matrix descriptor
    type(TDenseDescr), intent(in) :: denseDesc

    !> hamiltonian in sparse storage
    real(dp), intent(in) :: ham(:,:)

    !> sparse overlap matrix
    real(dp), intent(in) :: over(:)

    !> k-points
    real(dp), intent(in) :: kPoint(:,:)

    !> list of neighbours for each atom
    type(TNeighborList), intent(in) :: neighborList

    !> Number of neighbours for each of the atoms
    integer, intent(in) :: nNeighbor(:)

    !> Index array for the start of atomic blocks in sparse arrays
    integer, intent(in) :: iSparseStart(:,:)

    !> map from image atoms to the original unique atom
    integer, intent(in) :: img2CentCell(:)

    !> Index for which unit cell atoms are associated with
    integer, intent(in) :: iCellVec(:)

    !> Vectors (in units of the lattice constants) to cells of the lattice
    real(dp), intent(in) :: cellVec(:,:)

    !> Eigensolver choice
    type(TSolver), intent(in) :: solver

    !> K-points and spins to be handled
    type(TParallelKS), intent(in) :: parallelKS

    !> dense hamitonian matrix
    complex(dp), intent(out) :: HSqrCplx(:,:)

    !> dense overlap matrix
    complex(dp), intent(out) :: SSqrCplx(:,:)

    !> Complex eigenvectors
    complex(dp), intent(out) :: eigvecsCplx(:,:,:)

    !> eigenvalues
    real(dp), intent(out) :: eigen(:,:,:)

    integer :: iKS, iK, iSpin

    eigen(:,:,:) = 0.0_dp
    do iKS = 1, parallelKS%nLocalKS
      iK = parallelKS%localKS(1, iKS)
      iSpin = parallelKS%localKS(2, iKS)
    #:if WITH_SCALAPACK
      call env%globalTimer%startTimer(globalTimers%sparseToDense)
      call unpackHSCplxBlacs(env%blacs, ham(:,iSpin), kPoint(:,iK), neighborList%iNeighbor,&
          & nNeighbor, iCellVec, cellVec, iSparseStart, img2CentCell, denseDesc, HSqrCplx)
      call unpackHSCplxBlacs(env%blacs, over, kPoint(:,iK), neighborList%iNeighbor, nNeighbor,&
          & iCellVec, cellVec, iSparseStart, img2CentCell, denseDesc, SSqrCplx)
      call env%globalTimer%stopTimer(globalTimers%sparseToDense)
      call diagDenseMtxBlacs(solver, 'V', denseDesc%blacsOrbSqr, HSqrCplx, SSqrCplx,&
          & eigen(:,iK,iSpin), eigvecsCplx(:,:,iKS))
    #:else
      call env%globalTimer%startTimer(globalTimers%sparseToDense)
      call unpackHS(HSqrCplx, ham(:,iSpin), kPoint(:,iK), neighborList%iNeighbor, nNeighbor,&
          & iCellVec, cellVec, denseDesc%iAtomStart, iSparseStart, img2CentCell)
      call unpackHS(SSqrCplx, over, kPoint(:,iK), neighborList%iNeighbor, nNeighbor, iCellVec,&
          & cellVec, denseDesc%iAtomStart, iSparseStart, img2CentCell)
      call env%globalTimer%stopTimer(globalTimers%sparseToDense)
      call diagDenseMtx(solver, 'V', HSqrCplx, SSqrCplx, eigen(:,iK,iSpin))
      eigvecsCplx(:,:,iKS) = HSqrCplx
    #:endif
    end do

  #:if WITH_SCALAPACK
    call mpifx_allreduceip(env%mpi%interGroupComm, eigen, MPI_SUM)
  #:endif

  end subroutine buildAndDiagDenseCplxHam


  !> Builds and diagonalizes Pauli two-component Hamiltonians.
  subroutine buildAndDiagDensePauliHam(env, denseDesc, ham, over, kPoint, neighborList,&
      & nNeighbor, iSparseStart, img2CentCell, iCellVec, cellVec, orb, solver, parallelKS, eigen,&
      & HSqrCplx, SSqrCplx, eigvecsCplx, iHam, xi, species)

    !> Environment settings
    type(TEnvironment), intent(inout) :: env

    !> Dense matrix descriptor
    type(TDenseDescr), intent(in) :: denseDesc

    !> hamiltonian in sparse storage
    real(dp), intent(in) :: ham(:,:)

    !> sparse overlap matrix
    real(dp), intent(in) :: over(:)

    !> k-points
    real(dp), intent(in) :: kPoint(:,:)

    !> list of neighbours for each atom
    type(TNeighborList), intent(in) :: neighborList

    !> Number of neighbours for each of the atoms
    integer, intent(in) :: nNeighbor(:)

    !> Index array for the start of atomic blocks in sparse arrays
    integer, intent(in) :: iSparseStart(:,:)

    !> map from image atoms to the original unique atom
    integer, intent(in) :: img2CentCell(:)

    !> Index for which unit cell atoms are associated with
    integer, intent(in) :: iCellVec(:)

    !> Vectors (in units of the lattice constants) to cells of the lattice
    real(dp), intent(in) :: cellVec(:,:)

    !> atomic orbital information
    type(TOrbitals), intent(in) :: orb

    !> Eigensolver choice
    type(TSolver), intent(in) :: solver

    !> K-points and spins to be handled
    type(TParallelKS), intent(in) :: parallelKS

    !> eigenvalues
    real(dp), intent(out) :: eigen(:,:)

    !> dense hamitonian matrix
    complex(dp), intent(out) :: HSqrCplx(:,:)

    !> dense overlap matrix
    complex(dp), intent(out) :: SSqrCplx(:,:)

    !> eigenvectors
    complex(dp), intent(out) :: eigvecsCplx(:,:,:)

    !> imaginary part of the hamiltonian
    real(dp), intent(in), allocatable :: iHam(:,:)

    !> spin orbit constants
    real(dp), intent(in), allocatable :: xi(:,:)

    !> species of atoms
    integer, intent(in), optional :: species(:)

    integer :: iKS, iK

    eigen(:,:) = 0.0_dp
    do iKS = 1, parallelKS%nLocalKS
      iK = parallelKS%localKS(1, iKS)
      call env%globalTimer%startTimer(globalTimers%sparseToDense)
    #:if WITH_SCALAPACK
      if (allocated(iHam)) then
        call unpackHPauliBlacs(env%blacs, ham, kPoint(:,iK), neighborList%iNeighbor, nNeighbor,&
            & iCellVec, cellVec, iSparseStart, img2CentCell, orb%mOrb, denseDesc, HSqrCplx,&
            & iorig=iHam)
      else
        call unpackHPauliBlacs(env%blacs, ham, kPoint(:,iK), neighborList%iNeighbor, nNeighbor,&
            & iCellVec, cellVec, iSparseStart, img2CentCell, orb%mOrb, denseDesc, HSqrCplx)
      end if
      call unpackSPauliBlacs(env%blacs, over, kPoint(:,iK), neighborList%iNeighbor, nNeighbor,&
          & iCellVec, cellVec, iSparseStart, img2CentCell, orb%mOrb, denseDesc, SSqrCplx)
    #:else
      if (allocated(iHam)) then
        call unpackHPauli(ham, kPoint(:,iK), neighborList%iNeighbor, nNeighbor, iSparseStart, &
            & denseDesc%iAtomStart, img2CentCell, iCellVec, cellVec, HSqrCplx, iHam=iHam)
      else
        call unpackHPauli(ham, kPoint(:,iK), neighborList%iNeighbor, nNeighbor, iSparseStart, &
            & denseDesc%iAtomStart, img2CentCell, iCellVec, cellVec, HSqrCplx)
      end if
      call unpackSPauli(over, kPoint(:,iK), neighborList%iNeighbor, nNeighbor,&
          & denseDesc%iAtomStart, iSparseStart, img2CentCell, iCellVec, cellVec, SSqrCplx)
    #:endif
      if (allocated(xi) .and. .not. allocated(iHam)) then
        call addOnsiteSpinOrbitHam(env, xi, species, orb, denseDesc, HSqrCplx)
      end if
      call env%globalTimer%stopTimer(globalTimers%sparseToDense)
    #:if WITH_SCALAPACK
      call diagDenseMtxBlacs(solver, 'V', denseDesc%blacsOrbSqr, HSqrCplx, SSqrCplx, eigen(:,iK),&
          & eigvecsCplx(:,:,iKS))
    #:else
      call diagDenseMtx(solver, 'V', HSqrCplx, SSqrCplx, eigen(:,iK))
      eigvecsCplx(:,:,iKS) = HSqrCplx
    #:endif
    end do

  #:if WITH_SCALAPACK
    call mpifx_allreduceip(env%mpi%interGroupComm, eigen, MPI_SUM)
  #:endif

  end subroutine buildAndDiagDensePauliHam


  !> Creates sparse density matrix from real eigenvectors.
  subroutine getDensityFromRealEigvecs(env, denseDesc, filling, neighborList, nNeighbor,&
      & iSparseStart, img2CentCell, orb, eigvecs, parallelKS, rhoPrim, work, rhoSqrReal)

    !> Environment settings
    type(TEnvironment), intent(inout) :: env

    !> Dense matrix descriptor
    type(TDenseDescr), intent(in) :: denseDesc

    !> Filling
    real(dp), intent(in) :: filling(:,:)

    !> list of neighbours for each atom
    type(TNeighborList), intent(in) :: neighborList

    !> Number of neighbours for each of the atoms
    integer, intent(in) :: nNeighbor(:)

    !> Index array for the start of atomic blocks in sparse arrays
    integer, intent(in) :: iSparseStart(:,:)

    !> map from image atoms to the original unique atom
    integer, intent(in) :: img2CentCell(:)

    !> Atomic orbital information
    type(TOrbitals), intent(in) :: orb

    !> K-points and spins to process
    type(TParallelKS), intent(in) :: parallelKS

    !> eigenvectors
    real(dp), intent(inout) :: eigvecs(:,:,:)

    !> sparse density matrix
    real(dp), intent(out) :: rhoPrim(:,:)

    !> work space array
    real(dp), intent(out) :: work(:,:)

    !> Dense density matrix if needed
    real(dp), intent(inout), allocatable  :: rhoSqrReal(:,:,:)

    integer :: iKS, iSpin

    rhoPrim(:,:) = 0.0_dp
    do iKS = 1, parallelKS%nLocalKS
      iSpin = parallelKS%localKS(2, iKS)

    #:if WITH_SCALAPACK
      call makeDensityMtxRealBlacs(env%blacs%orbitalGrid, denseDesc%blacsOrbSqr, filling(:,iSpin),&
          & eigvecs(:,:,iKS), work)
      call env%globalTimer%startTimer(globalTimers%denseToSparse)
      call packRhoRealBlacs(env%blacs, denseDesc, work, neighborList%iNeighbor, nNeighbor,&
          & orb%mOrb, iSparseStart, img2CentCell, rhoPrim(:,iSpin))
      call env%globalTimer%stopTimer(globalTimers%denseToSparse)
    #:else
      if (tDensON2) then
        call makeDensityMatrix(work, eigvecs(:,:,iKS), filling(:,iSpin),&
            & neighborlist%iNeighbor, nNeighbor, orb, denseDesc%iAtomStart, img2CentCell)
      else
        call makeDensityMatrix(work, eigvecs(:,:,iKS), filling(:,iSpin))
      end if
      call env%globalTimer%startTimer(globalTimers%denseToSparse)
      call packHS(rhoPrim(:,iSpin), work, neighborlist%iNeighbor, nNeighbor, orb%mOrb,&
          & denseDesc%iAtomStart, iSparseStart, img2CentCell)
      call env%globalTimer%stopTimer(globalTimers%denseToSparse)
    #:endif

      if (allocated(rhoSqrReal)) then
        rhoSqrReal(:,:,iSpin) = work
      end if
    end do

  #:if WITH_SCALAPACK
    ! Add up and distribute density matrix contribution from each group
    call mpifx_allreduceip(env%mpi%globalComm, rhoPrim, MPI_SUM)
  #:endif

  end subroutine getDensityFromRealEigvecs


  !> Creates sparse density matrix from complex eigenvectors.
  subroutine getDensityFromCplxEigvecs(env, denseDesc, filling, kPoint, kWeight, neighborList,&
      & nNeighbor, iSparseStart, img2CentCell, iCellVec, cellVec, orb, parallelKS, eigvecs,&
      & rhoPrim, work)

    !> Environment settings
    type(TEnvironment), intent(inout) :: env

    !> Dense matrix descriptor
    type(TDenseDescr), intent(in) :: denseDesc

    !> Occupations of single particle states in the ground state
    real(dp), intent(in) :: filling(:,:,:)

    !> k-points
    real(dp), intent(in) :: kPoint(:,:)

    !> Weights for k-points
    real(dp), intent(in) :: kWeight(:)

    !> list of neighbours for each atom
    type(TNeighborList), intent(in) :: neighborList

    !> Number of neighbours for each of the atoms
    integer, intent(in) :: nNeighbor(:)

    !> Index array for the start of atomic blocks in sparse arrays
    integer, intent(in) :: iSparseStart(:,:)

    !> map from image atoms to the original unique atom
    integer, intent(in) :: img2CentCell(:)

    !> Index for which unit cell atoms are associated with
    integer, intent(in) :: iCellVec(:)

    !> Vectors (in units of the lattice constants) to cells of the lattice
    real(dp), intent(in) :: cellVec(:,:)

    !> Atomic orbital information
    type(TOrbitals), intent(in) :: orb

    !> K-points and spins to process
    type(TParallelKS), intent(in) :: parallelKS

    !> eigenvectors of the system
    complex(dp), intent(inout) :: eigvecs(:,:,:)

    !> density matrix in sparse storage
    real(dp), intent(out) :: rhoPrim(:,:)

    !> workspace array
    complex(dp), intent(out) :: work(:,:)

    integer :: iKS, iK, iSpin

    rhoPrim(:,:) = 0.0_dp

    do iKS = 1, parallelKS%nLocalKS
      iK = parallelKS%localKS(1, iKS)
      iSpin = parallelKS%localKS(2, iKS)
    #:if WITH_SCALAPACK
      call makeDensityMtxCplxBlacs(env%blacs%orbitalGrid, denseDesc%blacsOrbSqr,&
          & filling(:,iK,iSpin), eigvecs(:,:,iKS), work)
      call env%globalTimer%startTimer(globalTimers%denseToSparse)
      call packRhoCplxBlacs(env%blacs, denseDesc, work, kPoint(:,iK), kWeight(iK),&
          & neighborList%iNeighbor, nNeighbor, orb%mOrb, iCellVec, cellVec, iSparseStart,&
          & img2CentCell, rhoPrim(:,iSpin))
      call env%globalTimer%stopTimer(globalTimers%denseToSparse)
    #:else
      if (tDensON2) then
        call makeDensityMatrix(work, eigvecs(:,:,iKS), filling(:,iK,iSpin), neighborlist%iNeighbor,&
            & nNeighbor, orb, denseDesc%iAtomStart, img2CentCell)
      else
        call makeDensityMatrix(work, eigvecs(:,:,iKS), filling(:,iK,iSpin))
      end if
      call env%globalTimer%startTimer(globalTimers%denseToSparse)
      call packHS(rhoPrim(:,iSpin), work, kPoint(:,iK), kWeight(iK), neighborList%iNeighbor,&
          & nNeighbor, orb%mOrb, iCellVec, cellVec, denseDesc%iAtomStart, iSparseStart,&
          & img2CentCell)
      call env%globalTimer%stopTimer(globalTimers%denseToSparse)
    #:endif
    end do

  #:if WITH_SCALAPACK
    ! Add up and distribute density matrix contribution from each group
    call mpifx_allreduceip(env%mpi%globalComm, rhoPrim, MPI_SUM)
  #:endif

  end subroutine getDensityFromCplxEigvecs


  !> Creates sparse density matrix from two component complex eigenvectors.
  subroutine getDensityFromPauliEigvecs(env, denseDesc, tRealHS, tSpinOrbit, tDualSpinOrbit,&
      & tMulliken, kPoint, kWeight, filling, neighborList, nNeighbor, orb, iSparseStart,&
      & img2CentCell, iCellVec, cellVec, species, parallelKS, eigvecs, work, energy, rhoPrim, xi,&
      & orbitalL, iRhoPrim)

    !> Environment settings
    type(TEnvironment), intent(inout) :: env

    !> Dense matrix descriptor
    type(TDenseDescr), intent(in) :: denseDesc

    !> Is the hamitonian real (no k-points/molecule/gamma point)?
    logical, intent(in) :: tRealHS

    !> Are spin orbit interactions present
    logical, intent(in) :: tSpinOrbit

    !> Are block population spin orbit interactions present
    logical, intent(in) :: tDualSpinOrbit

    !> Should Mulliken populations be generated/output
    logical, intent(in) :: tMulliken

    !> k-points
    real(dp), intent(in) :: kPoint(:,:)

    !> Weights for k-points
    real(dp), intent(in) :: kWeight(:)

    !> occupations of molecular orbitals/Bloch states
    real(dp), intent(in) :: filling(:,:)

    !> list of neighbours for each atom
    type(TNeighborList), intent(in) :: neighborList

    !> Number of neighbours for each of the atoms
    integer, intent(in) :: nNeighbor(:)

    !> Atomic orbital information
    type(TOrbitals), intent(in) :: orb

    !> Index array for the start of atomic blocks in sparse arrays
    integer, intent(in) :: iSparseStart(:,:)

    !> map from image atoms to the original unique atom
    integer, intent(in) :: img2CentCell(:)

    !> Index for which unit cell atoms are associated with
    integer, intent(in) :: iCellVec(:)

    !> Vectors (in units of the lattice constants) to cells of the lattice
    real(dp), intent(in) :: cellVec(:,:)

    !> species of all atoms in the system
    integer, intent(in) :: species(:)

    !> K-points and spins to process
    type(TParallelKS), intent(in) :: parallelKS

    !> eigenvectors
    complex(dp), intent(inout) :: eigvecs(:,:,:)

    !> work space array
    complex(dp), intent(inout) :: work(:,:)

    !> Energy contributions and total
    type(TEnergies), intent(inout) :: energy

    !> sparse stored density matrix
    real(dp), intent(out) :: rhoPrim(:,:)

    !> spin orbit constants
    real(dp), intent(in), allocatable :: xi(:,:)

    !> Angular momentum of atomic shells
    real(dp), intent(inout), allocatable :: orbitalL(:,:,:)

    !> imaginary part of density matrix  if required
    real(dp), intent(inout), allocatable :: iRhoPrim(:,:)


    real(dp), allocatable :: rVecTemp(:), orbitalLPart(:,:,:)
    integer :: nKPoint, nAtom
    integer :: iKS, iK
    logical :: tImHam

    nAtom = size(orb%nOrbAtom)
    tImHam = allocated(iRhoPrim)
    nKPoint = size(kWeight)

    rhoPrim(:,:) = 0.0_dp
    if (allocated(iRhoPrim)) then
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

    do iKS = 1, parallelKS%nLocalKS
      iK = parallelKS%localKS(1, iKS)

    #:if WITH_SCALAPACK
      call makeDensityMtxCplxBlacs(env%blacs%orbitalGrid, denseDesc%blacsOrbSqr, filling(:,iK),&
          & eigvecs(:,:,iKS), work)
    #:else
      call makeDensityMatrix(work, eigvecs(:,:,iKS), filling(:,iK))
    #:endif
      if (tSpinOrbit .and. .not. tDualSpinOrbit) then
        call getOnsiteSpinOrbitEnergy(env, rVecTemp, work, denseDesc, xi, orb, species)
        energy%atomLS = energy%atomLS + kWeight(iK) * rVecTemp
        if (tMulliken) then
          orbitalLPart(:,:,:) = 0.0_dp
          call getLOnsite(env, orbitalLPart, work, denseDesc, orb, species)
          orbitalL(:,:,:) = orbitalL + kWeight(iK) * orbitalLPart
        end if
      end if

      call env%globalTimer%startTimer(globalTimers%denseToSparse)
    #:if WITH_SCALAPACK
      if (tImHam) then
        call packRhoPauliBlacs(env%blacs, denseDesc, work, kPoint(:,iK), kWeight(iK),&
            & neighborList%iNeighbor, nNeighbor, orb%mOrb, iCellVec, cellVec, iSparseStart,&
            & img2CentCell, rhoPrim, iRhoPrim)
      else
        call packRhoPauliBlacs(env%blacs, denseDesc, work, kPoint(:,iK), kWeight(iK),&
            & neighborList%iNeighbor, nNeighbor, orb%mOrb, iCellVec, cellVec, iSparseStart,&
            & img2CentCell, rhoPrim)
      end if
    #:else
      if (tRealHS) then
        call packHSPauli(rhoPrim, work, neighborlist%iNeighbor, nNeighbor, orb%mOrb,&
            & denseDesc%iAtomStart, iSparseStart, img2CentCell)
        if (tImHam) then
          call packHSPauliImag(iRhoPrim, work, neighborlist%iNeighbor, nNeighbor, orb%mOrb,&
              & denseDesc%iAtomStart, iSparseStart, img2CentCell)
        end if
      else
        call packHS(rhoPrim, work, kPoint(:,iK), kWeight(iK), neighborList%iNeighbor, nNeighbor,&
            & orb%mOrb, iCellVec, cellVec, denseDesc%iAtomStart, iSparseStart, img2CentCell)
        if (tImHam) then
          call iPackHS(iRhoPrim, work, kPoint(:,iK), kWeight(iK), neighborlist%iNeighbor, &
              & nNeighbor, orb%mOrb, iCellVec, cellVec, denseDesc%iAtomStart, iSparseStart,&
              & img2CentCell)
        end if
      end if
    #:endif
      call env%globalTimer%stopTimer(globalTimers%denseToSparse)
    end do

  #:if WITH_SCALAPACK
    call env%globalTimer%startTimer(globalTimers%denseToSparse)
    ! Add up and distribute contributions from each group
    call mpifx_allreduceip(env%mpi%globalComm, rhoPrim, MPI_SUM)
    if (allocated(iRhoPrim)) then
      call mpifx_allreduceip(env%mpi%globalComm, iRhoPrim, MPI_SUM)
    end if
    call mpifx_allreduceip(env%mpi%globalComm, energy%atomLS, MPI_SUM)
    if (tMulliken .and. tSpinOrbit .and. .not. tDualSpinOrbit) then
      call mpifx_allreduceip(env%mpi%globalComm, orbitalL, MPI_SUM)
    end if
    call env%globalTimer%stopTimer(globalTimers%denseToSparse)
  #:endif
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


  !> Calculate Mulliken population from sparse density matrix.
  subroutine getMullikenPopulation(rhoPrim, over, orb, neighborList, nNeighbor, img2CentCell,&
      & iSparseStart, qOrb, iRhoPrim, qBlock, qiBlock)

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
    integer, intent(in) :: iSparseStart(:,:)

    !> orbital charges
    real(dp), intent(out) :: qOrb(:,:,:)

    !> imaginary part of density matrix
    real(dp), intent(in), allocatable :: iRhoPrim(:,:)

    !> Dual atomic charges
    real(dp), intent(inout), allocatable :: qBlock(:,:,:,:)

    !> Imaginary part of dual atomic charges
    real(dp), intent(inout), allocatable :: qiBlock(:,:,:,:)

    integer :: iSpin

    qOrb(:,:,:) = 0.0_dp
    do iSpin = 1, size(qOrb, dim=3)
      call mulliken(qOrb(:,:,iSpin), over, rhoPrim(:,iSpin), orb, neighborList%iNeighbor,&
          & nNeighbor, img2CentCell, iSparseStart)
    end do

    if (allocated(qBlock)) then
      qBlock(:,:,:,:) = 0.0_dp
      do iSpin = 1, size(qBlock, dim=4)
        call mulliken(qBlock(:,:,:,iSpin), over, rhoPrim(:,iSpin), orb, neighborList%iNeighbor,&
            & nNeighbor, img2CentCell, iSparseStart)
      end do
    end if

    if (allocated(qiBlock)) then
      qiBlock(:,:,:,:) = 0.0_dp
      do iSpin = 1, size(qiBlock, dim=4)
        call skewMulliken(qiBlock(:,:,:,iSpin), over, iRhoPrim(:,iSpin), orb,&
            & neighborList%iNeighbor, nNeighbor, img2CentCell, iSparseStart)
      end do
    end if

  end subroutine getMullikenPopulation


  !> Calculates various energy contributions
  subroutine getEnergies(sccCalc, qOrb, q0, chargePerShell, species, tEField, tXlbomd,&
      & tDftbU, tDualSpinOrbit, rhoPrim, H0, orb, neighborList, nNeighbor, img2CentCell,&
      & iSparseStart, cellVol, extPressure, TS, potential, energy, thirdOrd, qBlock, qiBlock,&
      & nDftbUFunc, UJ, nUJ, iUJ, niUJ, xi)

    !> SCC module internal variables
    type(TScc), allocatable, intent(in) :: sccCalc

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
    integer, intent(in) :: iSparseStart(:,:)

    !> unit cell volume
    real(dp), intent(in) :: cellVol

    !> external pressure
    real(dp), intent(in) :: extPressure

    !> electron entropy contribution
    real(dp), intent(in) :: TS(:)

    !> potentials acting
    type(TPotentials), intent(in) :: potential

    !> energy contributions
    type(TEnergies), intent(inout) :: energy

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

    integer :: nSpin

    nSpin = size(rhoPrim, dim=2)

    ! Tr[H0 * Rho] can be done with the same algorithm as Mulliken-analysis
    energy%atomNonSCC(:) = 0.0_dp
    call mulliken(energy%atomNonSCC, rhoPrim(:,1), H0, orb, neighborList%iNeighbor, nNeighbor,&
        & img2CentCell, iSparseStart)
    energy%EnonSCC =  sum(energy%atomNonSCC)

    if (tEfield) then
      energy%atomExt = sum(qOrb(:,:,1) - q0(:,:,1), dim=1) * potential%extAtom(:,1)
      energy%Eext = sum(energy%atomExt)
    end if

    if (allocated(sccCalc)) then
      if (tXlbomd) then
        call sccCalc%getEnergyPerAtomXlbomd(species, orb, qOrb, q0, energy%atomSCC)
      else
        call sccCalc%getEnergyPerAtom(energy%atomSCC)
      end if
      energy%Escc = sum(energy%atomSCC)
      if (nSpin > 1) then
        energy%atomSpin(:) = 0.5_dp * sum(sum(potential%intShell(:,:,2:nSpin)&
            & * chargePerShell(:,:,2:nSpin), dim=1), dim=2)
        energy%Espin = sum(energy%atomSpin)
      end if
    end if
    if (allocated(thirdOrd)) then
      if (tXlbomd) then
        call thirdOrd%getEnergyPerAtomXlbomd(qOrb, q0, species, orb, energy%atom3rd)
      else
        call thirdOrd%getEnergyPerAtom(energy%atom3rd)
      end if
      energy%e3rd = sum(energy%atom3rd)
    end if

    if (tDftbU) then
      if (allocated(qiBlock)) then
        call E_DFTBU(energy%atomDftbu, qBlock, species, orb, nDFTBUfunc, UJ, nUJ, niUJ, iUJ,&
            & qiBlock)
      else
        call E_DFTBU(energy%atomDftbu, qBlock, species, orb, nDFTBUfunc, UJ, nUJ, niUJ, iUJ)
      end if
      energy%Edftbu = sum(energy%atomDftbu)
    end if

    if (tDualSpinOrbit) then
      energy%atomLS(:) = 0.0_dp
      call getDualSpinOrbitEnergy(energy%atomLS, qiBlock, xi, orb, species)
      energy%ELS = sum(energy%atomLS)
    end if

    energy%Eelec = energy%EnonSCC + energy%ESCC + energy%Espin + energy%ELS + energy%Edftbu&
        & + energy%Eext + energy%e3rd
    energy%atomElec(:) = energy%atomNonSCC + energy%atomSCC + energy%atomSpin + energy%atomDftbu&
        & + energy%atomLS + energy%atomExt + energy%atom3rd
    energy%atomTotal(:) = energy%atomElec + energy%atomRep + energy%atomDisp
    energy%Etotal = energy%Eelec + energy%Erep + energy%eDisp
    energy%EMermin = energy%Etotal - sum(TS)
    energy%EGibbs = energy%EMermin + cellVol * extPressure

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
  subroutine getNextInputCharges(env, pChrgMixer, qOutput, qOutRed, orb, nIneqOrb, iEqOrbitals,&
      & iGeoStep, iSccIter, minSccIter, maxSccIter, sccTol, tStopScc, tDftbU, tReadChrg, qInput,&
      & qInpRed, sccErrorQ, tConverged, qBlockOut, iEqBlockDftbU, qBlockIn, qiBlockOut,&
      & iEqBlockDftbuLS, species0, nUJ, iUJ, niUJ, qiBlockIn)

    !> Environment settings
    type(TEnvironment), intent(in) :: env

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
    real(dp), intent(inout), allocatable :: qBlockOut(:,:,:,:)

    !> equivalence mapping for dual charge blocks
    integer, intent(in), allocatable :: iEqBlockDftbu(:,:,:,:)

    !> block charge input (if needed for orbital potentials)
    real(dp), intent(inout), allocatable ::qBlockIn(:,:,:,:)

    !> Imaginary part of block charges
    real(dp), intent(in), allocatable :: qiBlockOut(:,:,:,:)

    !> Equivalence mappings in the case of spin orbit and DFTB+U
    integer, intent(in), allocatable :: iEqBlockDftbuLS(:,:,:,:)

    !> atomic species for atoms
    integer, intent(in), allocatable :: species0(:)

    !> Number DFTB+U blocks of shells for each atom type
    integer, intent(in), allocatable :: nUJ(:)

    !> which shells are in each DFTB+U block
    integer, intent(in), allocatable :: iUJ(:,:,:)

    !> Number of shells in each DFTB+U block
    integer, intent(in), allocatable :: niUJ(:,:)

    !> Imaginary part of block atomic input populations
    real(dp), intent(inout), allocatable :: qiBlockIn(:,:,:,:)

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
        if (allocated(qBlockIn)) then
          qBlockIn(:,:,:,:) = qBlockOut
          if (allocated(qiBlockIn)) then
            qiBlockIn(:,:,:,:) = qiBlockOut
          end if
        end if
      else
        call mix(pChrgMixer, qInpRed, qDiffRed)
      #:if WITH_MPI
        ! Synchronise charges in order to avoid mixers that store a history drifting apart
        call mpifx_allreduceip(env%mpi%globalComm, qInpRed, MPI_SUM)
        qInpRed(:) = qInpRed / env%mpi%globalComm%size
      #:endif
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
    real(dp), intent(in) :: qOrb(:,:,:)

    !> Reduction of atomic populations
    real(dp), intent(out) :: qRed(:)

    !> Block (dual) populations, if also being reduced
    real(dp), intent(in), allocatable :: qBlock(:,:,:,:)

    !> equivalences for block charges
    integer, intent(in), allocatable :: iEqBlockDftbu(:,:,:,:)

    !> Imaginary part of block charges if present
    real(dp), intent(in), allocatable :: qiBlock(:,:,:,:)

    !> Equivalences for spin orbit if needed
    integer, intent(in), allocatable :: iEqBlockDftbuLS(:,:,:,:)

    real(dp), allocatable :: qOrbUpDown(:,:,:), qBlockUpDown(:,:,:,:)

    qRed(:) = 0.0_dp
    qOrbUpDown = qOrb
    call qm2ud(qOrbUpDown)
    call orbitalEquiv_reduce(qOrbUpDown, iEqOrbitals, orb, qRed(1:nIneqOrb))
    if (allocated(qBlock)) then
      qBlockUpDown = qBlock
      call qm2ud(qBlockUpDown)
      call appendBlock_reduce(qBlockUpDown, iEqBlockDFTBU, orb, qRed)
      if (allocated(qiBlock)) then
        call appendBlock_reduce(qiBlock, iEqBlockDFTBULS, orb, qRed, skew=.true.)
      end if
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
    real(dp), intent(inout), allocatable :: qBlock(:,:,:,:)

    !> equivalences for block charges
    integer, intent(in), allocatable :: iEqBlockDftbU(:,:,:,:)

    !> species of central cell atoms
    integer, intent(in), allocatable :: species0(:)

    !> Number DFTB+U blocks of shells for each atom type
    integer, intent(in), allocatable :: nUJ(:)

    !> which shells are in each DFTB+U block
    integer, intent(in), allocatable :: iUJ(:,:,:)

    !> Number of shells in each DFTB+U block
    integer, intent(in), allocatable :: niUJ(:,:)

    !> Imaginary part of block atomic populations
    real(dp), intent(inout), allocatable :: qiBlock(:,:,:,:)

    !> Equivalences for spin orbit if needed
    integer, intent(in), allocatable :: iEqBlockDftbULS(:,:,:,:)

    integer :: nSpin

    @:ASSERT(allocated(qBlock) .eqv. allocated(iEqBlockDftbU))
    @:ASSERT(.not. allocated(qBlock) .or. allocated(species0))
    @:ASSERT(.not. allocated(qBlock) .or. allocated(nUJ))
    @:ASSERT(.not. allocated(qBlock) .or. allocated(iUJ))
    @:ASSERT(.not. allocated(qBlock) .or. allocated(niUJ))
    @:ASSERT(.not. allocated(qiBlock) .or. allocated(qBlock))
    @:ASSERT(allocated(qiBlock) .eqv. allocated(iEqBlockDftbuLS))

    nSpin = size(qOrb, dim=3)
    call OrbitalEquiv_expand(qRed(1:nIneqOrb), iEqOrbitals, orb, qOrb)
    if (allocated(qBlock)) then
      qBlock(:,:,:,:) = 0.0_dp
      call Block_expand(qRed, iEqBlockDftbu, orb, qBlock, species0, nUJ, niUJ, iUJ,&
          & orbEquiv=iEqOrbitals)
      if (allocated(qiBlock)) then
        call Block_expand(qRed, iEqBlockDftbuLS, orb, qiBlock, species0, nUJ, niUJ, iUJ,&
            & skew=.true.)
      end if
    end if
    if (nSpin == 2) then
      call ud2qm(qOrb)
      if (allocated(qBlock)) then
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

    logical :: tEnoughIters, tRestartIter

    ! Do we need restart at all?
    tRestart = (restartFreq > 0 .and. .not. (tMD .or. tGeoOpt .or. tDerivs) .and. maxSccIter > 1)
    if (tRestart) then

      ! Do we have enough iterations already?
      tEnoughIters = (iSccIter >= minSccIter .or. tReadChrg .or. iGeoStep > 0)

      ! Is current iteration the right one for writing a restart file?
      tRestartIter = (iSccIter == maxSccIter .or. tStopScc .or. mod(iSccIter, restartFreq) == 0)

      tRestart = (tConverged .or. (tEnoughIters .and. tRestartIter))
    end if

  end function needsSccRestartWriting


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
  subroutine calculateLinRespExcitations(env, lresp, parallelKS, sccCalc, qOutput, q0, over,&
      & eigvecsReal, eigen, filling, coord0, species, speciesName, orb, skHamCont, skOverCont,&
      & fdAutotest, autotestTag, fdEigvec, runId, neighborList, nNeighbor, denseDesc, iSparseStart,&
      & img2CentCell, tWriteAutotest, tForces, tLinRespZVect, tPrintExcEigvecs,&
      & tPrintExcEigvecsTxt, nonSccDeriv, energy, work, rhoSqrReal, excitedDerivs, occNatural)

    !> Environment settings
    type(TEnvironment), intent(in) :: env

    !> excited state settings
    type(LinResp), intent(inout) :: lresp

    !> K-points and spins to process
    type(TParallelKS), intent(in) :: parallelKS

    !> SCC module internal variables
    type(TScc), intent(in) :: sccCalc

    !> electrons in atomic orbitals
    real(dp), intent(in) :: qOutput(:,:,:)

    !> reference atomic orbital occupations
    real(dp), intent(in) :: q0(:,:,:)

    !> sparse overlap matrix
    real(dp), intent(in) :: over(:)

    !> ground state eigenvectors
    real(dp), intent(in) :: eigvecsReal(:,:,:)

    !> ground state eigenvalues
    real(dp), intent(in) :: eigen(:,:)

    !> ground state fillings
    real(dp), intent(in) :: filling(:,:)

    !> central cell coordinates
    real(dp), intent(in) :: coord0(:,:)

    !> species of all atoms in the system
    integer, target, intent(in) :: species(:)

    !> label for each atomic chemical species
    character(*), intent(in) :: speciesName(:)

    !> Atomic orbital information
    type(TOrbitals), intent(in) :: orb

    !> non-SCC hamiltonian information
    type(OSlakoCont), intent(in) :: skHamCont

    !> overlap information
    type(OSlakoCont), intent(in) :: skOverCont

    !> file ID for regression data
    integer, intent(in) :: fdAutotest

    !> File name for regression data
    character(*), intent(in) :: autotestTag

    !> File ID for ground state eigenvectors
    integer, intent(in) :: fdEigvec

    !> Job ID for future identification
    integer, intent(in) :: runId

    !> list of neighbours for each atom
    type(TNeighborList), intent(in) :: neighborList

    !> Number of neighbours for each of the atoms
    integer, intent(in) :: nNeighbor(:)

    !> Dense matrix descriptor
    type(TDenseDescr), intent(in) :: denseDesc

    !> Index array for the start of atomic blocks in sparse arrays
    integer, intent(in) :: iSparseStart(:,:)

    !> map from image atoms to the original unique atom
    integer, intent(in) :: img2CentCell(:)

    !> should regression test data be written
    logical, intent(in) :: tWriteAutotest

    !> forces to be calculated
    logical, intent(in) :: tForces

    !> require the Z vector for excited state properties
    logical, intent(in) :: tLinRespZVect

    !> print natural orbitals of the excited state
    logical, intent(in) :: tPrintExcEigvecs

    !> print natural orbitals also in text form?
    logical, intent(in) :: tPrintExcEigvecsTxt

    !> method for calculating derivatives of S and H0
    type(NonSccDiff), intent(in) :: nonSccDeriv

    !> Energy contributions and total
    type(TEnergies), intent(inout) :: energy

    !> Working array of the size of the dense matrices.
    real(dp), intent(out) :: work(:,:)

    !> density matrix in dense form
    real(dp), intent(inout), allocatable :: rhoSqrReal(:,:,:)

    !> excited state energy derivative with respect to atomic coordinates
    real(dp), intent(inout), allocatable :: excitedDerivs(:,:)

    !> natural orbital occupation numbers
    real(dp), intent(inout), allocatable :: occNatural(:)

    real(dp), allocatable :: dQAtom(:)
    real(dp), allocatable :: naturalOrbs(:,:,:)
    integer, pointer :: pSpecies0(:)
    integer :: iSpin, nSpin, nAtom
    logical :: tSpin

    nAtom = size(qOutput, dim=2)
    nSpin = size(eigen, dim=2)
    tSpin = (nSpin == 2)
    pSpecies0 => species(1:nAtom)

    energy%Eexcited = 0.0_dp
    allocate(dQAtom(nAtom))
    dQAtom(:) = sum(qOutput(:,:,1) - q0(:,:,1), dim=1)
    call unpackHS(work, over, neighborList%iNeighbor, nNeighbor, denseDesc%iAtomStart,&
        & iSparseStart, img2CentCell)
    call blockSymmetrizeHS(work, denseDesc%iAtomStart)
    if (tForces) then
      do iSpin = 1, nSpin
        call blockSymmetrizeHS(rhoSqrReal(:,:,iSpin), denseDesc%iAtomStart)
      end do
    end if
    if (tWriteAutotest) then
      open(fdAutotest, file=autotestTag, position="append")
    end if

    if (tLinRespZVect) then
      if (tPrintExcEigVecs) then
        allocate(naturalOrbs(orb%nOrb, orb%nOrb, 1))
      end if
      call addGradients(tSpin, lresp, denseDesc%iAtomStart, eigvecsReal, eigen, work, filling,&
          & coord0, sccCalc, dQAtom, pSpecies0, neighborList%iNeighbor, img2CentCell, orb,&
          & skHamCont, skOverCont, tWriteAutotest, fdAutotest, energy%Eexcited, excitedDerivs,&
          & nonSccDeriv, rhoSqrReal, occNatural, naturalOrbs)
      if (tPrintExcEigvecs) then
        call writeRealEigvecs(env, fdEigvec, runId, neighborList, nNeighbor, denseDesc,&
            & iSparseStart, img2CentCell, pSpecies0, speciesName, orb, over, parallelKS,&
            & tPrintExcEigvecsTxt, naturalOrbs, work, fileName="excitedOrbs")
      end if
    else
      call calcExcitations(lresp, tSpin, denseDesc, eigvecsReal, eigen, work, filling, coord0,&
          & sccCalc, dQAtom, pSpecies0, neighborList%iNeighbor, img2CentCell, orb, tWriteAutotest,&
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

    !> integrator for the extended Lagrangian
    type(Xlbomd), intent(inout) :: xlbomdIntegrator

    !> output charges, reduced by equivalences
    real(dp), intent(in) :: qOutRed(:)

    !> SCC mixer
    type(OMixer), intent(inout) :: pChrgMixer

    !> Atomic orbital information
    type(TOrbitals), intent(in) :: orb

    !> number of inequivalent orbitals
    integer, intent(in) :: nIneqOrb

    !> equivalence map
    integer, intent(in) :: iEqOrbitals(:,:,:)

    !> input charges
    real(dp), intent(out) :: qInput(:,:,:)

    !> input charges reduced by equivalences
    real(dp), intent(out) :: qInpRed(:)

    !> +U equivalences
    integer, intent(in), allocatable :: iEqBlockDftbU(:,:,:,:)

    !> central cell species
    integer, intent(in), allocatable :: species0(:)

    !> block input charges
    real(dp), intent(inout), allocatable :: qBlockIn(:,:,:,:)

    !> Number DFTB+U blocks of shells for each atom type
    integer, intent(in), allocatable :: nUJ(:)

    !> which shells are in each DFTB+U block
    integer, intent(in), allocatable :: iUJ(:,:,:)

    !> Number of shells in each DFTB+U block
    integer, intent(in), allocatable :: niUJ(:,:)

    !> equivalences for spin orbit
    integer, intent(in), allocatable :: iEqBlockDftbuLS(:,:,:,:)

    !> imaginary part of dual charges
    real(dp), intent(inout), allocatable :: qiBlockIn(:,:,:,:)

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


  !> Calculates dipole moment.
  subroutine getDipoleMoment(qOutput, q0, coord, dipoleMoment)

    !> electrons in orbitals
    real(dp), intent(in) :: qOutput(:,:,:)

    !> reference atomic charges
    real(dp), intent(in) :: q0(:,:,:)

    !> atomic coordinates
    real(dp), intent(in) :: coord(:,:)

    !> resulting dipole moment
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
      & neighborList, nNeighbor, species, iSparseStart, img2CentCell)
    integer, intent(in) :: sparseSize
    real(dp), intent(in) :: rhoPrim(:,:)
    real(dp), intent(in) :: q0(:,:,:)
    real(dp), intent(in) :: coord0(:,:)
    real(dp), intent(in) :: over(:)

    !> Atomic orbital information
    type(TOrbitals), intent(in) :: orb

    !> list of neighbours for each atom
    type(TNeighborList), intent(in) :: neighborList

    !> Number of neighbours for each of the atoms
    integer, intent(in) :: nNeighbor(:)

    !> species of all atoms in the system
    integer, intent(in) :: species(:)

    !> Index array for the start of atomic blocks in sparse arrays
    integer, intent(in) :: iSparseStart(:,:)

    !> map from image atoms to the original unique atom
    integer, intent(in) :: img2CentCell(:)

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
      call add_shift(hprime, over, nNeighbor, neighborList%iNeighbor, species, orb, iSparseStart,&
          & nAtom, img2CentCell, potentialDerivative)

      ! evaluate <psi| dH/dE | psi>
      call mulliken(dipole, hprime(:,1), rhoPrim(:,1), orb, neighborList%iNeighbor, nNeighbor,&
          & img2CentCell, iSparseStart)

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
  subroutine getEnergyWeightedDensity(env, denseDesc, forceType, filling, eigen, kPoint, kWeight,&
      & neighborList, nNeighbor, orb, iSparseStart, img2CentCell, iCellVEc, cellVec, tRealHS, ham,&
      & over, parallelKS, ERhoPrim, HSqrReal, SSqrReal, HSqrCplx, SSqrCplx, rhoSqrReal)

    !> Environment settings
    type(TEnvironment), intent(in) :: env

    !> Dense matrix descriptor
    type(TDenseDescr), intent(in) :: denseDesc

    !> Force type
    integer, intent(in) :: forceType

    !> Occupations of single particle states in the ground state
    real(dp), intent(in) :: filling(:,:,:)

    !> Eigenvalues
    real(dp), intent(in) :: eigen(:,:,:)

    !> K-points
    real(dp), intent(in) :: kPoint(:,:)

    !> Weights for k-points
    real(dp), intent(in) :: kWeight(:)

    !> list of neighbours for each atom
    type(TNeighborList), intent(in) :: neighborList

    !> Number of neighbours for each of the atoms
    integer, intent(in) :: nNeighbor(:)

    !> Atomic orbital information
    type(TOrbitals), intent(in) :: orb

    !> Index array for the start of atomic blocks in sparse arrays
    integer, intent(in) :: iSparseStart(:,:)

    !> map from image atoms to the original unique atom
    integer, intent(in) :: img2CentCell(:)

    !> Index for which unit cell atoms are associated with
    integer, intent(in) :: iCellVec(:)

    !> Vectors (in units of the lattice constants) to cells of the lattice
    real(dp), intent(in) :: cellVec(:,:)

    !> Is the hamitonian real (no k-points/molecule/gamma point)?
    logical, intent(in) :: tRealHS

    !> Sparse Hamiltonian
    real(dp), intent(in) :: ham(:,:)

    !> Sparse overlap
    real(dp), intent(in) :: over(:)

    !> K-points and spins to process
    type(TParallelKS), intent(in) :: parallelKS

    !> Energy weighted sparse matrix
    real(dp), intent(out) :: ERhoPrim(:)

    !> Storage for dense hamiltonian matrix
    real(dp), intent(inout), allocatable :: HSqrReal(:,:,:)

    !> Storage for dense overlap matrix
    real(dp), intent(inout), allocatable :: SSqrReal(:,:)

    !> Storage for dense hamitonian matrix (complex case)
    complex(dp), intent(inout), allocatable :: HSqrCplx(:,:,:)

    !> Storage for dense overlap matrix (complex case)
    complex(dp), intent(inout), allocatable :: SSqrCplx(:,:)

    !> Density matrix in case it has been stored
    real(dp), allocatable, intent(inout) :: rhoSqrReal(:,:,:)

    integer :: nSpin

    nSpin = size(ham, dim=2)

    if (nSpin == 4) then
      call getEDensityMtxFromPauliEigvecs(env, denseDesc, filling, eigen, kPoint, kWeight,&
          & neighborList, nNeighbor, orb, iSparseStart, img2CentCell, iCellVec, cellVec, tRealHS,&
          & parallelKS, HSqrCplx, SSqrCplx, ERhoPrim)
    else if (tRealHS) then
      call getEDensityMtxFromRealEigvecs(env, denseDesc, forceType, filling, eigen, neighborList,&
          & nNeighbor, orb, iSparseStart, img2CentCell, ham, over, parallelKS, HSqrReal, SSqrReal,&
          & ERhoPrim, rhoSqrReal)
    else
      call getEDensityMtxFromComplexEigvecs(env, denseDesc, forceType, filling, eigen, kPoint,&
          & kWeight, neighborList, nNeighbor, orb, iSparseStart, img2CentCell, iCellVec, cellVec,&
          & ham, over, parallelKS, HSqrCplx, SSqrCplx, ERhoPrim)
    end if

  end subroutine getEnergyWeightedDensity


  !> Calculates density matrix from real eigenvectors.
  subroutine getEDensityMtxFromRealEigvecs(env, denseDesc, forceType, filling, eigen, neighborList,&
      & nNeighbor, orb, iSparseStart, img2CentCell, ham, over, parallelKS, eigvecsReal, work,&
      & ERhoPrim, rhoSqrReal)

    !> Environment settings
    type(TEnvironment), intent(in) :: env

    !> Dense matrix descriptor
    type(TDenseDescr), intent(in) :: denseDesc

    !> How to calculate the force
    integer, intent(in) :: forceType

    !> Occupations of single particle states in the ground state
    real(dp), intent(in) :: filling(:,:,:)

    !> Eigenvalues
    real(dp), intent(in) :: eigen(:,:,:)

    !> list of neighbours for each atom
    type(TNeighborList), intent(in) :: neighborList

    !> Number of neighbours for each of the atoms
    integer, intent(in) :: nNeighbor(:)

    !> Atomic orbital information
    type(TOrbitals), intent(in) :: orb

    !> Index array for the start of atomic blocks in sparse arrays
    integer, intent(in) :: iSparseStart(:,:)

    !> map from image atoms to the original unique atom
    integer, intent(in) :: img2CentCell(:)

    !> Sparse Hamiltonian
    real(dp), intent(in) :: ham(:,:)

    !> Sparse overlap
    real(dp), intent(in) :: over(:)

    !> K-points and spins to process
    type(TParallelKS), intent(in) :: parallelKS

    !> Eigenvectors (NOTE: they will be rewritten with work data on exit!)
    real(dp), intent(inout) :: eigvecsReal(:,:,:)

    !> Work array for storing temporary data
    real(dp), intent(out) :: work(:,:)

    !> Energy weighted density matrix
    real(dp), intent(out) :: ERhoPrim(:)

    !> Density matrix in case it has been stored
    real(dp), allocatable, intent(inout) :: rhoSqrReal(:,:,:)

    real(dp), allocatable :: work2(:,:)
    integer :: nSpin, nLocalRows, nLocalCols, nOrb
    integer :: iKS, iS

    nSpin = size(eigen, dim=3)
    nLocalRows = size(eigvecsReal, dim=1)
    nLocalCols = size(eigvecsReal, dim=2)
    nOrb = denseDesc%iAtomStart(size(denseDesc%iAtomStart)) - 1
    if (forceType == 2 .or. forceType == 3) then
      allocate(work2(nLocalRows, nLocalCols))
    end if

    ERhoPrim(:) = 0.0_dp
    do iKS = 1, parallelKS%nLocalKS
      iS = parallelKS%localKS(2, iKS)

      select case (forceType)

      case(0)
        ! Original (non-consistent) scheme
      #:if WITH_SCALAPACK
        call makeDensityMtxRealBlacs(env%blacs%orbitalGrid, denseDesc%blacsOrbSqr, filling(:,1,iS),&
            & eigvecsReal(:,:,iKS), work, eigen(:,1,iS))
      #:else
        if (tDensON2) then
          call makeDensityMatrix(work, eigvecsReal(:,:,iKS), filling(:,1,iS), eigen(:,1,iS),&
              & neighborlist%iNeighbor, nNeighbor, orb, denseDesc%iAtomStart, img2CentCell)
        else
          call makeDensityMatrix(work, eigvecsReal(:,:,iKS), filling(:,1,iS), eigen(:,1,iS))
        end if
      #:endif

      case(2)
        ! Correct force for XLBOMD for T=0K (DHD)
      #:if WITH_SCALAPACK
        call unpackHSRealBlacs(env%blacs, ham(:,iS), neighborList%iNeighbor, nNeighbor,&
            & iSparseStart, img2CentCell, denseDesc, work)
        call makeDensityMtxRealBlacs(env%blacs%orbitalGrid, denseDesc%blacsOrbSqr, filling(:,1,iS),&
            & eigVecsReal(:,:,iKS), work2)
        call pblasfx_psymm(work2, denseDesc%blacsOrbSqr, work, denseDesc%blacsOrbSqr,&
            & eigvecsReal(:,:,iKS), denseDesc%blacsOrbSqr, side="L")
        call pblasfx_psymm(work2, denseDesc%blacsOrbSqr, eigvecsReal(:,:,iKS),&
            & denseDesc%blacsOrbSqr, work, denseDesc%blacsOrbSqr, side="R", alpha=0.5_dp)
      #:else
        call unpackHS(work, ham(:,iS), neighborlist%iNeighbor, nNeighbor, denseDesc%iAtomStart,&
            & iSparseStart, img2CentCell)
        call blockSymmetrizeHS(work, denseDesc%iAtomStart)
        if (allocated(rhoSqrReal)) then
          work2(:,:) = rhoSqrReal(:,:,iS)
        else
          call makeDensityMatrix(work2, eigvecsReal(:,:,iKS), filling(:,1,iS))
        end if
        ! D H
        call symm(eigvecsReal(:,:,iKS), "L", work2, work)
        ! (D H) D
        call symm(work, "R", work2, eigvecsReal(:,:,iKS), alpha=0.5_dp)
      #:endif

      case(3)
        ! Correct force for XLBOMD for T <> 0K (DHS^-1 + S^-1HD)
      #:if WITH_SCALAPACK
        call makeDensityMtxRealBlacs(env%blacs%orbitalGrid, denseDesc%blacsOrbSqr, filling(:,1,iS),&
            & eigVecsReal(:,:,iKS), work)
        call unpackHSRealBlacs(env%blacs, ham(:,iS), neighborlist%iNeighbor, nNeighbor,&
            & iSparseStart, img2CentCell, denseDesc, work2)
        call pblasfx_psymm(work, denseDesc%blacsOrbSqr, work2, denseDesc%blacsOrbSqr,&
            & eigvecsReal(:,:,iKS), denseDesc%blacsOrbSqr, side="L")
        call unpackHSRealBlacs(env%blacs, over, neighborlist%iNeighbor, nNeighbor, iSparseStart,&
            & img2CentCell, denseDesc, work)
        call psymmatinv(denseDesc%blacsOrbSqr, work)
        call pblasfx_psymm(work, denseDesc%blacsOrbSqr, eigvecsReal(:,:,iKS),&
            & denseDesc%blacsOrbSqr, work2, denseDesc%blacsOrbSqr, side="R", alpha=0.5_dp)
        work(:,:) = work2
        call pblasfx_ptran(work2, denseDesc%blacsOrbSqr, work, denseDesc%blacsOrbSqr, alpha=1.0_dp,&
            & beta=1.0_dp)
      #:else
        call makeDensityMatrix(work, eigvecsReal(:,:,iKS), filling(:,1,iS))
        call unpackHS(work2, ham(:,iS), neighborlist%iNeighbor, nNeighbor, denseDesc%iAtomStart,&
            & iSparseStart, img2CentCell)
        call blocksymmetrizeHS(work2, denseDesc%iAtomStart)
        call symm(eigvecsReal(:,:,iKS), "L", work, work2)
        call unpackHS(work, over, neighborlist%iNeighbor, nNeighbor, denseDesc%iAtomStart,&
            & iSparseStart, img2CentCell)
        call symmatinv(work)
        call symm(work2, "R", work, eigvecsReal(:,:,iKS), alpha=0.5_dp)
        work(:,:) = work2 + transpose(work2)
      #:endif
      end select

    #:if WITH_SCALAPACK
      call packRhoRealBlacs(env%blacs, denseDesc, work, neighborList%iNeighbor, nNeighbor,&
          & orb%mOrb, iSparseStart, img2CentCell, ERhoPrim)
    #:else
      call packHS(ERhoPrim, work, neighborList%iNeighbor, nNeighbor, orb%mOrb,&
          & denseDesc%iAtomStart, iSparseStart, img2CentCell)
    #:endif
    end do

  #:if WITH_SCALAPACK
    ! Add up and distribute energy weighted density matrix contribution from each group
    call mpifx_allreduceip(env%mpi%globalComm, ERhoPrim, MPI_SUM)
  #:endif

  end subroutine getEDensityMtxFromRealEigvecs


  !> Calculates density matrix from complex eigenvectors.
  subroutine getEDensityMtxFromComplexEigvecs(env, denseDesc, forceType, filling, eigen, kPoint,&
      & kWeight, neighborList, nNeighbor, orb, iSparseStart, img2CentCell, iCellVec, cellVec,&
      & ham, over, parallelKS, eigvecsCplx, work, ERhoPrim)

    !> Environment settings
    type(TEnvironment), intent(in) :: env

    !> Dense matrix descriptor
    type(TDenseDescr), intent(in) :: denseDesc

    integer, intent(in) :: forceType

    !> Occupations of single particle states in the ground state
    real(dp), intent(in) :: filling(:,:,:)

    !> eigen-values of the system
    real(dp), intent(in) :: eigen(:,:,:)

    !> k-points of the system
    real(dp), intent(in) :: kPoint(:,:)

    !> Weights for k-points
    real(dp), intent(in) :: kWeight(:)

    !> list of neighbours for each atom
    type(TNeighborList), intent(in) :: neighborList

    !> Number of neighbours for each of the atoms
    integer, intent(in) :: nNeighbor(:)

    !> Atomic orbital information
    type(TOrbitals), intent(in) :: orb

    !> Index array for the start of atomic blocks in sparse arrays
    integer, intent(in) :: iSparseStart(:,:)

    !> map from image atoms to the original unique atom
    integer, intent(in) :: img2CentCell(:)

    !> Index for which unit cell atoms are associated with
    integer, intent(in) :: iCellVec(:)

    !> Vectors (in units of the lattice constants) to cells of the lattice
    real(dp), intent(in) :: cellVec(:,:)

    !> Sparse Hamiltonian
    real(dp), intent(in) :: ham(:,:)

    !> Sparse overlap
    real(dp), intent(in) :: over(:)

    !> K-points and spins to process
    type(TParallelKS), intent(in) :: parallelKS

    !> Eigenvectors of the system
    complex(dp), intent(inout) :: eigvecsCplx(:,:,:)

    !> work array (sized like overlap matrix)
    complex(dp), intent(inout) :: work(:,:)

    !> Energy weighted sparse density matrix (charge only part)
    real(dp), intent(out) :: ERhoPrim(:)

    complex(dp), allocatable :: work2(:,:)
    integer :: nLocalRows, nLocalCols, nOrb
    integer :: iKS, iS, iK

    nLocalRows = size(eigvecsCplx, dim=1)
    nLocalCols = size(eigvecsCplx, dim=2)
    nOrb = denseDesc%iAtomStart(size(denseDesc%iAtomStart)) - 1

    if (forceType == 2 .or. forceType == 3) then
      allocate(work2(nLocalRows, nLocalCols))
    end if

    ERhoPrim(:) = 0.0_dp
    do iKS = 1, parallelKS%nLocalKS
      iK = parallelKS%localKS(1, iKS)
      iS = parallelKS%localKS(2, iKS)

      select case (forceType)

      case(0)
        ! Original (non-consistent) scheme
      #:if WITH_SCALAPACK
        call makeDensityMtxCplxBlacs(env%blacs%orbitalGrid, denseDesc%blacsOrbSqr,&
            & filling(:,iK,iS), eigvecsCplx(:,:,iKS), work, eigen(:,iK,iS))
      #:else
        if (tDensON2) then
          call makeDensityMatrix(work, eigvecsCplx(:,:,iKS), filling(:,iK,iS),&
              & eigen(:,iK, iS), neighborlist%iNeighbor, nNeighbor, orb, denseDesc%iAtomStart,&
              & img2CentCell)
        else
          call makeDensityMatrix(work, eigvecsCplx(:,:,iKS), filling(:,iK,iS), eigen(:,iK, iS))
        end if
      #:endif

      case(2)
        ! Correct force for XLBOMD for T=0K (DHD)
      #:if WITH_SCALAPACK
        call unpackHSCplxBlacs(env%blacs, ham(:,iS), kPoint(:,iK), neighborList%iNeighbor,&
            & nNeighbor, iCellVec, cellVec, iSparseStart, img2CentCell, denseDesc, work)
        call makeDensityMtxCplxBlacs(env%blacs%orbitalGrid, denseDesc%blacsOrbSqr, filling(:,1,iS),&
            & eigvecsCplx(:,:,iKS), work2)
        call pblasfx_phemm(work2, denseDesc%blacsOrbSqr, work, denseDesc%blacsOrbSqr,&
            & eigvecsCplx(:,:,iKS), denseDesc%blacsOrbSqr, side="L")
        call pblasfx_phemm(work2, denseDesc%blacsOrbSqr, eigvecsCplx(:,:,iKS),&
            & denseDesc%blacsOrbSqr, work, denseDesc%blacsOrbSqr, side="R", alpha=(0.5_dp, 0.0_dp))
      #:else
        call makeDensityMatrix(work2, eigvecsCplx(:,:,iKS), filling(:,iK,iS))
        call unpackHS(work, ham(:,iS), kPoint(:,iK), neighborlist%iNeighbor, nNeighbor,&
            & iCellVec, cellVec, denseDesc%iAtomStart, iSparseStart, img2CentCell)
        call blockHermitianHS(work, denseDesc%iAtomStart)
        call hemm(eigvecsCplx(:,:,iKS), "L", work2, work)
        call hemm(work, "R", work2, eigvecsCplx(:,:,iKS), alpha=(0.5_dp, 0.0_dp))
      #:endif

      case(3)
        ! Correct force for XLBOMD for T <> 0K (DHS^-1 + S^-1HD)
      #:if WITH_SCALAPACK
        call makeDensityMtxCplxBlacs(env%blacs%orbitalGrid, denseDesc%blacsOrbSqr,&
            & filling(:,iK,iS), eigVecsCplx(:,:,iKS), work)
        call unpackHSCplxBlacs(env%blacs, ham(:,iS), kPoint(:,iK), neighborlist%iNeighbor,&
            & nNeighbor, iCellVec, cellVec, iSparseStart, img2CentCell, denseDesc, work2)
        call pblasfx_phemm(work, denseDesc%blacsOrbSqr, work2, denseDesc%blacsOrbSqr,&
            & eigvecsCplx(:,:,iKS), denseDesc%blacsOrbSqr, side="L")
        call unpackHSCplxBlacs(env%blacs, over, kPoint(:,iK), neighborlist%iNeighbor,&
            & nNeighbor, iCellVec, cellVec, iSparseStart, img2CentCell, denseDesc, work)
        call phermatinv(denseDesc%blacsOrbSqr, work)
        call pblasfx_phemm(work, denseDesc%blacsOrbSqr, eigvecsCplx(:,:,iKS),&
            & denseDesc%blacsOrbSqr, work2, denseDesc%blacsOrbSqr, side="R", alpha=(0.5_dp, 0.0_dp))
        work(:,:) = work2
        call pblasfx_ptranc(work2, denseDesc%blacsOrbSqr, work, denseDesc%blacsOrbSqr,&
            & alpha=(1.0_dp, 0.0_dp), beta=(1.0_dp, 0.0_dp))
      #:else
        call makeDensityMatrix(work, eigvecsCplx(:,:,iKS), filling(:,iK,iS))
        call unpackHS(work2, ham(:,iS), kPoint(:,iK), neighborlist%iNeighbor, nNeighbor,&
            & iCellVec, cellVec, denseDesc%iAtomStart, iSparseStart, img2CentCell)
        call blockHermitianHS(work2, denseDesc%iAtomStart)
        call hemm(eigvecsCplx(:,:,iKS), "L", work, work2)
        call unpackHS(work, over, kPoint(:,iK), neighborlist%iNeighbor, nNeighbor, iCellVec,&
            & cellVec, denseDesc%iAtomStart, iSparseStart, img2CentCell)
        call hermatinv(work)
        call hemm(work2, "R", work, eigvecsCplx(:,:,iKS), alpha=(0.5_dp, 0.0_dp))
        work(:,:) = work2 + transpose(conjg(work2))
      #:endif
      end select

    #:if WITH_SCALAPACK
      call packRhoCplxBlacs(env%blacs, denseDesc, work, kPoint(:,iK), kWeight(iK),&
          &neighborList%iNeighbor, nNeighbor, orb%mOrb, iCellVec, cellVec, iSparseStart,&
          & img2CentCell, ERhoPrim)
    #:else
      call packHS(ERhoPrim, work, kPoint(:,iK), kWeight(iK), neighborList%iNeighbor,&
          & nNeighbor, orb%mOrb, iCellVec, cellVec, denseDesc%iAtomStart, iSparseStart,&
          & img2CentCell)
    #:endif
    end do

  #:if WITH_SCALAPACK
    ! Add up and distribute energy weighted density matrix contribution from each group
    call mpifx_allreduceip(env%mpi%globalComm, ERhoPrim, MPI_SUM)
  #:endif

  end subroutine getEDensityMtxFromComplexEigvecs


  !> Calculates density matrix from Pauli-type two component eigenvectors.
  subroutine getEDensityMtxFromPauliEigvecs(env, denseDesc, filling, eigen, kPoint, kWeight,&
      & neighborList, nNeighbor, orb, iSparseStart, img2CentCell, iCellVec, cellVec, tRealHS,&
      & parallelKS, eigvecsCplx, work, ERhoPrim)

    !> Environment settings
    type(TEnvironment), intent(in) :: env

    !> Dense matrix descriptor
    type(TDenseDescr), intent(in) :: denseDesc

    !> Occupations of single particle states in the ground state
    real(dp), intent(in) :: filling(:,:,:)

    !> Eigenvalues
    real(dp), intent(in) :: eigen(:,:,:)

    !> K-points
    real(dp), intent(in) :: kPoint(:,:)

    !> Weights for k-points
    real(dp), intent(in) :: kWeight(:)

    !> list of neighbours for each atom
    type(TNeighborList), intent(in) :: neighborList

    !> Number of neighbours for each of the atoms
    integer, intent(in) :: nNeighbor(:)

    !> Atomic orbital information
    type(TOrbitals), intent(in) :: orb

    !> Index array for the start of atomic blocks in sparse arrays
    integer, intent(in) :: iSparseStart(:,:)

    !> map from image atoms to the original unique atom
    integer, intent(in) :: img2CentCell(:)

    !> Index for which unit cell atoms are associated with
    integer, intent(in) :: iCellVec(:)

    !> Vectors (in units of the lattice constants) to cells of the lattice
    real(dp), intent(in) :: cellVec(:,:)

    !> Is the hamitonian real (no k-points/molecule/gamma point)?
    logical, intent(in) :: tRealHS

    !> K-points and spins to process
    type(TParallelKS), intent(in) :: parallelKS

    !> Eigenvectors
    complex(dp), intent(inout) :: eigvecsCplx(:,:,:)

    !> Work array
    complex(dp), intent(out) :: work(:,:)

    !> Sparse energy weighted density matrix
    real(dp), intent(out) :: ERhoPrim(:)

    integer :: iKS, iK

    ERhoPrim(:) = 0.0_dp
    do iKS = 1, parallelKS%nLocalKS
      iK = parallelKS%localKS(1, iKS)
    #:if WITH_SCALAPACK
      call makeDensityMtxCplxBlacs(env%blacs%orbitalGrid, denseDesc%blacsOrbSqr, filling(:,iK,1),&
          & eigvecsCplx(:,:,iKS), work, eigen(:,iK,1))
      call packERhoPauliBlacs(env%blacs, denseDesc, work, kPoint(:,iK), kWeight(iK),&
          & neighborList%iNeighbor, nNeighbor, orb%mOrb, iCellVec, cellVec, iSparseStart,&
          & img2CentCell, ERhoPrim)
    #:else
      call makeDensityMatrix(work, eigvecsCplx(:,:,iKS), filling(:,iK,1), eigen(:,iK,1))
      if (tRealHS) then
        call packERho(ERhoPrim, work, neighborList%iNeighbor, nNeighbor, orb%mOrb,&
            & denseDesc%iAtomStart, iSparseStart, img2CentCell)
      else
        call packERho(ERhoPrim, work, kPoint(:,iK), kWeight(iK), neighborList%iNeighbor,&
            & nNeighbor, orb%mOrb, iCellVec, cellVec, denseDesc%iAtomStart, iSparseStart,&
            & img2CentCell)
      end if
    #:endif
    end do

  #:if WITH_SCALAPACK
    ! Add up and distribute energy weighted density matrix contribution from each group
    call mpifx_allreduceip(env%mpi%globalComm, ERhoPrim, MPI_SUM)
  #:endif

  end subroutine getEDensityMtxFromPauliEigvecs


  !> Calculates the gradients
  subroutine getGradients(env, sccCalc, tEField, tXlbomd, nonSccDeriv, Efield, rhoPrim, ERhoPrim,&
      & qOutput, q0, skHamCont, skOverCont, pRepCont, neighborList, nNeighbor, species,&
      & img2CentCell, iSparseStart, orb, potential, coord, derivs, iRhoPrim, thirdOrd, chrgForces,&
      & dispersion)

    !> Environment settings
    type(TEnvironment), intent(in) :: env

    !> SCC module internal variables
    type(TScc), allocatable, intent(in) :: sccCalc

    !> external electric field
    logical, intent(in) :: tEField

    !> extended Lagrangian active?
    logical, intent(in) :: tXlbomd

    !> method for calculating derivatives of S and H0
    type(NonSccDiff), intent(in) :: nonSccDeriv

    !> Any applied electric field
    real(dp), intent(in) :: Efield(:)

    !> sparse density matrix
    real(dp), intent(in) :: rhoPrim(:,:)

    !> energy  weighted density matrix
    real(dp), intent(in) :: ERhoPrim(:)

    !> electron populations (may be unallocated for non-scc case)
    real(dp), allocatable, intent(in) :: qOutput(:,:,:)

    !> reference atomic charges (may be unallocated for non-scc case)
    real(dp), allocatable, intent(in) :: q0(:,:,:)

    !> non-SCC hamiltonian information
    type(OSlakoCont), intent(in) :: skHamCont

    !> overlap information
    type(OSlakoCont), intent(in) :: skOverCont

    !> repulsive information
    type(ORepCont), intent(in) :: pRepCont

    !> list of neighbours for each atom
    type(TNeighborList), intent(in) :: neighborList

    !> Number of neighbours for each of the atoms
    integer, intent(in) :: nNeighbor(:)

    !> species of all atoms in the system
    integer, intent(in) :: species(:)

    !> map from image atoms to the original unique atom
    integer, intent(in) :: img2CentCell(:)

    !> Index array for the start of atomic blocks in sparse arrays
    integer, intent(in) :: iSparseStart(:,:)

    !> Atomic orbital information
    type(TOrbitals), intent(in) :: orb

    !>  potential acting on the system
    type(TPotentials), intent(in) :: potential

    !> atomic coordinates
    real(dp), intent(in) :: coord(:,:)

    !> derivatives of energy wrt to atomic positions
    real(dp), intent(out) :: derivs(:,:)

    !> imaginary part of density matrix
    real(dp), intent(in), allocatable :: iRhoPrim(:,:)

    !> Is 3rd order SCC being used
    type(ThirdOrder), intent(inout), allocatable :: thirdOrd

    !> forces on external charges
    real(dp), intent(inout), allocatable :: chrgForces(:,:)

    !> dispersion interactions
    class(DispersionIface), intent(inout), allocatable :: dispersion

    real(dp), allocatable :: tmpDerivs(:,:)
    logical :: tImHam, tExtChrg, tSccCalc
    integer :: nAtom
    integer :: ii

    tSccCalc = allocated(sccCalc)
    tImHam = allocated(iRhoPrim)
    tExtChrg = allocated(chrgForces)
    nAtom = size(derivs, dim=2)

    derivs(:,:) = 0.0_dp

    if (.not. (tSccCalc .or. tEField)) then
      ! No external or internal potentials
      if (tImHam) then
        call derivative_shift(env, derivs, nonSccDeriv, rhoPrim, iRhoPrim, ERhoPrim, skHamCont,&
            & skOverCont, coord, species, neighborList%iNeighbor, nNeighbor, img2CentCell,&
            & iSparseStart, orb, potential%intBlock, potential%iorbitalBlock)
      else
        call derivative_shift(env, derivs, nonSccDeriv, rhoPrim(:,1), ERhoPrim, skHamCont,&
            & skOverCont, coord, species, neighborList%iNeighbor, nNeighbor, img2CentCell,&
            & iSparseStart, orb)
      end if
    else
      if (tImHam) then
        call derivative_shift(env, derivs, nonSccDeriv, rhoPrim, iRhoPrim, ERhoPrim, skHamCont,&
            & skOverCont, coord, species, neighborList%iNeighbor, nNeighbor, img2CentCell,&
            & iSparseStart, orb, potential%intBlock, potential%iorbitalBlock)
      else
        call derivative_shift(env, derivs, nonSccDeriv, rhoPrim, ERhoPrim, skHamCont, skOverCont,&
            & coord, species, neighborList%iNeighbor, nNeighbor, img2CentCell, iSparseStart, orb,&
            & potential%intBlock)
      end if

      if (tExtChrg) then
        chrgForces(:,:) = 0.0_dp
        if (tXlbomd) then
          call error("XLBOMD does not work with external charges yet!")
        else
          call sccCalc%addForceDc(env, derivs, species, neighborList%iNeighbor, img2CentCell, &
              & chrgForces)
        end if
      else if (tSccCalc) then
        if (tXlbomd) then
          call sccCalc%addForceDcXlbomd(env, species, orb, neighborList%iNeighbor, img2CentCell,&
              & qOutput, q0, derivs)
        else
          call sccCalc%addForceDc(env, derivs, species, neighborList%iNeighbor, img2CentCell)
        end if
      end if

      if (allocated(thirdOrd)) then
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

    allocate(tmpDerivs(3, nAtom))
    call getERepDeriv(tmpDerivs, coord, nNeighbor, neighborList%iNeighbor, species, pRepCont,&
        & img2CentCell)
    derivs(:,:) = derivs + tmpDerivs

  end subroutine getGradients


  !> Calculates stress tensor and lattice derivatives.
  subroutine getStress(env, sccCalc, tEField, nonSccDeriv, EField, rhoPrim, ERhoPrim, qOutput,&
      & q0, skHamCont, skOverCont, pRepCont, neighborList, nNeighbor, species, img2CentCell,&
      & iSparseStart, orb, potential, coord, latVec, invLatVec, cellVol, coord0, totalStress,&
      & totalLatDeriv, intPressure, iRhoPrim, dispersion)

    !> Environment settings
    type(TEnvironment), intent(in) :: env

    !> SCC module internal variables
    type(TScc), allocatable, intent(in) :: sccCalc

    !> External electric field
    logical, intent(in) :: tEField

    !> method for calculating derivatives of S and H0
    type(NonSccDiff), intent(in) :: nonSccDeriv

    !> external electric field
    real(dp), intent(in) :: Efield(:)

    !> density matrix
    real(dp), intent(in) :: rhoPrim(:,:)

    !> energy weighted density matrix
    real(dp), intent(in) :: ERhoPrim(:)

    !> electrons in orbitals
    real(dp), intent(in) :: qOutput(:,:,:)

    !> refernce charges
    real(dp), intent(in) :: q0(:,:,:)

    !> non-SCC hamitonian information
    type(OSlakoCont), intent(in) :: skHamCont

    !> overlap information
    type(OSlakoCont), intent(in) :: skOverCont

    !> repulsive information
    type(ORepCont), intent(in) :: pRepCont

    !> list of neighbours for each atom
    type(TNeighborList), intent(in) :: neighborList

    !> Number of neighbours for each of the atoms
    integer, intent(in) :: nNeighbor(:)

    !> species of all atoms in the system
    integer, intent(in) :: species(:)

    !> map from image atoms to the original unique atom
    integer, intent(in) :: img2CentCell(:)

    !> Index array for the start of atomic blocks in sparse arrays
    integer, intent(in) :: iSparseStart(:,:)

    !> Atomic orbital information
    type(TOrbitals), intent(in) :: orb

    !> potentials acting
    type(TPotentials), intent(in) :: potential

    !> coordinates of all atoms
    real(dp), intent(in) :: coord(:,:)

    !> lattice vectors
    real(dp), intent(in) :: latVec(:,:)

    !> inverse of the lattice vectors
    real(dp), intent(in) :: invLatVec(:,:)

    !> unit cell volume
    real(dp), intent(in) :: cellVol

    !> central cell coordinates of atoms
    real(dp), intent(inout) :: coord0(:,:)

    !> stress tensor
    real(dp), intent(out) :: totalStress(:,:)

    !> energy derivatives with respect to lattice vectors
    real(dp), intent(out) :: totalLatDeriv(:,:)

    !> internal pressure in cell
    real(dp), intent(out) :: intPressure

    !> imaginary part of the density matrix (if present)
    real(dp), intent(in), allocatable :: iRhoPrim(:,:)

    !> dispersion interactions
    class(DispersionIface), allocatable, intent(inout) :: dispersion

    real(dp) :: tmpStress(3, 3)
    logical :: tImHam

    tImHam = allocated(iRhoPrim)

    if (allocated(sccCalc)) then
      if (tImHam) then
        call getBlockiStress(env, totalStress, nonSccDeriv, rhoPrim, iRhoPrim, ERhoPrim, skHamCont,&
            & skOverCont, coord, species, neighborList%iNeighbor, nNeighbor, img2CentCell,&
            & iSparseStart, orb, potential%intBlock, potential%iorbitalBlock, cellVol)
      else
        call getBlockStress(env, totalStress, nonSccDeriv, rhoPrim, ERhoPrim, skHamCont,&
            & skOverCont, coord, species, neighborList%iNeighbor, nNeighbor, img2CentCell,&
            & iSparseStart, orb, potential%intBlock, cellVol)
      end if
      call sccCalc%addStressDc(totalStress, env, species, neighborList%iNeighbor, img2CentCell)
    else
      if (tImHam) then
        call getBlockiStress(env, totalStress, nonSccDeriv, rhoPrim, iRhoPrim, ERhoPrim, skHamCont,&
            & skOverCont, coord, species, neighborList%iNeighbor, nNeighbor, img2CentCell,&
            & iSparseStart, orb, potential%intBlock, potential%iorbitalBlock, cellVol)
      else
        call getNonSCCStress(env, totalStress, nonSccDeriv, rhoPrim(:,1), ERhoPrim, skHamCont,&
            & skOverCont, coord, species, neighborList%iNeighbor, nNeighbor, img2CentCell,&
            & iSparseStart, orb, cellVol)
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

    intPressure = (totalStress(1,1) + totalStress(2,2) + totalStress(3,3)) / 3.0_dp
    totalLatDeriv(:,:) = -cellVol * matmul(totalStress, invLatVec)

  end subroutine getStress


  !> Calculates stress from external electric field.
  subroutine getEFieldStress(latVec, cellVol, q0, qOutput, Efield, coord0, stress)

    !> lattice vectors
    real(dp), intent(in) :: latVec(:,:)

    !> unit cell volume
    real(dp), intent(in) :: cellVol

    !> reference atomic charges
    real(dp), intent(in) :: q0(:,:,:)

    !> number of electrons in each orbital
    real(dp), intent(in) :: qOutput(:,:,:)

    !> external electric field
    real(dp), intent(in) :: Efield(:)

    !> central cell coordinates of atoms
    real(dp), intent(inout) :: coord0(:,:)

    !> Stress tensor
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


  !> Removes forces components along constraint directions
  subroutine constrainForces(conAtom, conVec, derivss)

    !> atoms being constrained
    integer, intent(in) :: conAtom(:)

    !> vector to project out forces
    real(dp), intent(in) :: conVec(:,:)

    !> on input energy derivatives, on exit resulting projected derivatives
    real(dp), intent(inout) :: derivss(:,:)

    integer :: ii, iAtom

    ! Set force components along constraint vectors zero
    do ii = 1, size(conAtom)
      iAtom = conAtom(ii)
      derivss(:,iAtom) = derivss(:,iAtom)&
          & - conVec(:,ii) * dot_product(conVec(:,ii), derivss(:,iAtom))
    end do

  end subroutine constrainForces


  !> Flattens lattice components and applies lattice optimisation constraints.
  subroutine constrainLatticeDerivs(totalLatDerivs, normLatVecs, tLatOptFixAng,&
      & tLatOptFixLen, tLatOptIsotropic, constrLatDerivs)

    !> energy derivative with respect to lattice vectors
    real(dp), intent(in) :: totalLatDerivs(:,:)

    !> unit normals parallel to lattice vectors
    real(dp), intent(in) :: normLatVecs(:,:)

    !> Are the angles of the lattice being fixed during optimisation?
    logical, intent(in) :: tLatOptFixAng

    !> Are the magnitude of the lattice vectors fixed
    logical, intent(in) :: tLatOptFixLen(:)

    !> is the optimisation isotropic
    logical, intent(in) :: tLatOptIsotropic

    !> lattice vectors returned by the optimizer
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


  !> Unfold contrained lattice vectors to full one.
  subroutine unconstrainLatticeVectors(constrLatVecs, origLatVecs, tLatOptFixAng, tLatOptFixLen,&
      & tLatOptIsotropic, newLatVecs)

    !> packaged up lattice vectors (depending on optimisation mode)
    real(dp), intent(in) :: constrLatVecs(:)

    !> original vectors at start
    real(dp), intent(in) :: origLatVecs(:,:)

    !> Are the angles of the lattice vectors fixed
    logical, intent(in) :: tLatOptFixAng

    !> are the magnitudes of the lattice vectors fixed
    logical, intent(in) :: tLatOptFixLen(:)

    !> is the optimisation isotropic
    logical, intent(in) :: tLatOptIsotropic

    !> resulting lattice vectors
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

    !> Driver for the finite difference second derivatives
    type(OnumDerivs), intent(inout) :: derivDriver

    !> first derivatives of energy at the current coordinates
    real(dp), intent(in) :: derivs(:,:)

    !> moving atoms
    integer, intent(in) :: indMovedAtoms(:)

    !> atomic coordinates
    real(dp), intent(out) :: coords(:,:)

    !> has the process terminated
    logical, intent(out) :: tGeomEnd

    real(dp) :: newCoords(3, size(indMovedAtoms))

    call next(derivDriver, newCoords, derivs(:, indMovedAtoms), tGeomEnd)
    coords(:, indMovedAtoms) = newCoords

  end subroutine getNextDerivStep


  !> Returns the coordinates for the next coordinate optimisation step.
  subroutine getNextCoordinateOptStep(pGeoCoordOpt, EMermin, derivss, indMovedAtom, coords0,&
      & diffGeo, tCoordEnd)

    !> optimiser for atomic coordinates
    type(OGeoOpt), intent(inout) :: pGeoCoordOpt

    !> electronic free energy U -TS
    real(dp), intent(in) :: EMermin

    !> Derivative of energy with respect to atomic coordinates
    real(dp), intent(in) :: derivss(:,:)

    !> numbers of the moving atoms
    integer, intent(in) :: indMovedAtom(:)

    !> central cell atomic coordinates
    real(dp), intent(inout) :: coords0(:,:)

    !> largest change in atomic coordinates
    real(dp), intent(out) :: diffGeo

    !> has the geometry optimisation finished
    logical, intent(out) :: tCoordEnd

    real(dp) :: derivssMoved(3 * size(indMovedAtom))
    real(dp), target :: newCoordsMoved(3 * size(indMovedAtom))
    real(dp), pointer :: pNewCoordsMoved(:,:)

    derivssMoved(:) = reshape(derivss(:, indMovedAtom), [3 * size(indMovedAtom)])
    call next(pGeoCoordOpt, EMermin, derivssMoved, newCoordsMoved, tCoordEnd)
    pNewCoordsMoved(1:3, 1:size(indMovedAtom)) => newCoordsMoved(1 : 3 * size(indMovedAtom))
    diffGeo = maxval(abs(pNewCoordsMoved - coords0(:, indMovedAtom)))
    coords0(:, indMovedAtom) = pNewCoordsMoved

  end subroutine getNextCoordinateOptStep


  !> Returns the coordinates and lattice vectors for the next lattice optimisation step.
  subroutine getNextLatticeOptStep(pGeoLatOpt, EGibbs, constrLatDerivs, origLatVec, tLatOptFixAng,&
      & tLatOptFixLen, tLatOptIsotropic, indMovedAtom, latVec, coord0, diffGeo, tGeomEnd)

    !> lattice vector optimising object
    type(OGeoOpt), intent(inout) :: pGeoLatOpt

    !> Gibbs free energy (U -TS_elec + pV)
    real(dp), intent(in) :: EGibbs

    !> lattice vectors returned by the optimizer
    real(dp), intent(in) :: constrLatDerivs(:)

    !> Starting lattice vectors
    real(dp), intent(in) :: origLatVec(:,:)

    !> Fix angles between lattice vectors
    logical, intent(in) :: tLatOptFixAng

    !> Fix the magnitudes of lattice vectors
    logical, intent(in) :: tLatOptFixLen(:)

    !> Optimise isotropically
    logical, intent(in) :: tLatOptIsotropic

    !> numbers of the moving atoms
    integer, intent(in) :: indMovedAtom(:)

    !> lattice vectors
    real(dp), intent(inout) :: latVec(:,:)

    !> central cell coordinates of atoms
    real(dp), intent(inout) :: coord0(:,:)

    !> Maximum change in geometry at this step
    real(dp), intent(out) :: diffGeo

    !> has the geometry optimisation finished
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
  subroutine getNextMdStep(pMdIntegrator, pMdFrame, temperatureProfile, derivs, movedMass,&
      & mass, cellVol, invLatVec, species0, indMovedAtom, tStress, tBarostat, energy, coord0,&
      & latVec, intPressure, totalStress, totalLatDeriv, velocities, tempIon)

    !> Molecular dynamics integrator
    type(OMdIntegrator), intent(inout) :: pMdIntegrator

    !> Molecular dynamics reference frame information
    type(OMdCommon), intent(in) :: pMdFrame

    !> Temperature profile in MD
    type(OTempProfile), allocatable, intent(inout) :: temperatureProfile

    !> Energy derivative wrt to atom positions
    real(dp), intent(in) :: derivs(:,:)

    !> Masses of moving atoms
    real(dp), intent(in) :: movedMass(:,:)

    !> Masses of each chemical species
    real(dp), intent(in) :: mass(:)

    !> unit cell volume
    real(dp), intent(in) :: cellVol

    !> inverse of the lattice vectors
    real(dp), intent(in) :: invLatVec(:,:)

    !> species of atoms in the central cell
    integer, intent(in) :: species0(:)

    !> numbers of the moving atoms
    integer, intent(in) :: indMovedAtom(:)

    !> Is stress being evaluated?
    logical, intent(in) :: tStress

    !> Is there a barostat
    logical, intent(in) :: tBarostat

    !> Energy contributions and total
    type(TEnergies), intent(inout) :: energy

    !> central cell coordinates of atoms
    real(dp), intent(inout) :: coord0(:,:)

    !> lattice vectors
    real(dp), intent(inout) :: latVec(:,:)

    !> Internal pressure in the unit cell
    real(dp), intent(inout) :: intPressure

    !> Stress tensor
    real(dp), intent(inout) :: totalStress(:,:)

    !> Derivative of energy with respect to lattice vectors
    real(dp), intent(inout) :: totalLatDeriv(:,:)

    !> Atomic velocities
    real(dp), intent(out) :: velocities(:,:)

    !> Atomic kinetic energy
    real(dp), intent(out) :: tempIon

    real(dp) :: movedAccel(3, size(indMovedAtom)), movedVelo(3, size(indMovedAtom))
    real(dp) :: movedCoords(3, size(indMovedAtom))
    real(dp) :: kineticStress(3, 3)

    movedAccel(:,:) = -derivs(:, indMovedAtom) / movedMass
    call next(pMdIntegrator, movedAccel, movedCoords, movedVelo)
    coord0(:, indMovedAtom) = movedCoords
    velocities(:,:) = 0.0_dp
    velocities(:, indMovedAtom) = movedVelo(:,:)

    if (allocated(temperatureProfile)) then
      call next(temperatureProfile)
    end if
    call evalKE(energy%Ekin, movedVelo, movedMass(1,:))
    call evalkT(pMdFrame, tempIon, movedVelo, movedMass(1,:))
    energy%EMerminKin = energy%EMermin + energy%Ekin
    energy%EGibbsKin = energy%EGibbs + energy%Ekin

    if (tStress) then
      ! contribution from kinetic energy in MD, now that velocities for this geometry step are
      ! available
      call getKineticStress(kineticStress, mass, species0, velocities, cellVol)
      totalStress = totalStress + kineticStress
      intPressure = (totalStress(1,1) + totalStress(2,2) + totalStress(3,3)) / 3.0_dp
      totalLatDeriv = -cellVol * matmul(totalStress, invLatVec)
    end if

    if (tBarostat) then
      call rescale(pMDIntegrator, coord0, latVec, totalStress)
    end if

  end subroutine getNextMdStep


  !> Calculates and prints Pipek-Mezey localisation
  subroutine calcPipekMezeyLocalisation(env, pipekMezey, tPrintEigvecsTxt, nEl, filling, over,&
      & kPoint, neighborList, nNeighbor, denseDesc,  iSparseStart, img2CentCell, iCellVec,&
      & cellVec, fdEigvec, runId, orb, species, speciesName, parallelKS, localisation, eigvecsReal,&
      & SSqrReal, eigvecsCplx, SSqrCplx)

    !> Environment settings
    type(TEnvironment), intent(in) :: env

    !> Localisation methods for single electron states (if used)
    type(TPipekMezey), intent(in) :: pipekMezey

    !> Store eigenvectors as a text file
    logical, intent(in) :: tPrintEigVecsTxt

    !> Number of electrons
    real(dp), intent(in) :: nEl(:)

    !> Occupations of single particle states in the ground state
    real(dp), intent(in) :: filling(:,:,:)

    !> sparse overlap matrix
    real(dp), intent(in) :: over(:)

    !> k-points in the system (0,0,0) if molecular
    real(dp), intent(in) :: kPoint(:,:)

    !> list of neighbours for each atom
    type(TNeighborList), intent(in) :: neighborList

    !> Number of neighbours for each of the atoms
    integer, intent(in) :: nNeighbor(:)

    !> Dense matrix descriptor
    type(TDenseDescr), intent(in) :: denseDesc

    !> Index array for the start of atomic blocks in sparse arrays
    integer, intent(in) :: iSparseStart(:,:)

    !> map from image atoms to the original unique atom
    integer, intent(in) :: img2CentCell(:)

    !> Index for which unit cell atoms are associated with
    integer, intent(in) :: iCellVec(:)

    !> Vectors (in units of the lattice constants) to cells of the lattice
    real(dp), intent(in) :: cellVec(:,:)

    !> File ID for ground state eigenvectors
    integer, intent(in) :: fdEigvec

    !> Job ID for future identification
    integer, intent(in) :: runId

    !> Atomic orbital information
    type(TOrbitals), intent(in) :: orb

    !> species of all atoms in the system
    integer, intent(in) :: species(:)

    !> label for each atomic chemical species
    character(*), intent(in) :: speciesName(:)

    !> K-points and spins to process
    type(TParallelKS), intent(in) :: parallelKS

    !> Localisation measure of single particle states
    real(dp), intent(out) :: localisation

    !> Storage for dense hamiltonian matrix
    real(dp), intent(inout), allocatable :: eigvecsReal(:,:,:)

    !> Storage for dense overlap matrix
    real(dp), intent(inout), allocatable :: SSqrReal(:,:)

    !> Storage for dense hamitonian matrix (complex case)
    complex(dp), intent(inout), allocatable :: eigvecsCplx(:,:,:)

    !> Storage for dense overlap matrix (complex case)
    complex(dp), intent(inout), allocatable :: SSqrCplx(:,:)

    integer :: nFilledLev, nAtom, nSpin
    integer :: iSpin, iKS, iK

    nAtom = size(orb%nOrbAtom)
    nSpin = size(nEl)

    if (any(abs(mod(filling, real(3 - nSpin, dp))) > elecTolMax)) then
      call warning("Fractional occupations allocated for electron localisation")
    end if

    if (allocated(eigvecsReal)) then
      call unpackHS(SSqrReal,over,neighborList%iNeighbor, nNeighbor, denseDesc%iAtomStart,&
          & iSparseStart, img2CentCell)
      do iKS = 1, parallelKS%nLocalKS
        iSpin = parallelKS%localKS(2, iKS)
        nFilledLev = floor(nEl(iSpin) / real(3 - nSpin, dp))
        localisation = pipekMezey%getLocalisation(eigvecsReal(:, 1:nFilledLev, iKS), SSqrReal,&
            & denseDesc%iAtomStart)
        write(stdOut, "(A, E15.8)") 'Original localisation', localisation
        call pipekMezey%calcCoeffs(eigvecsReal(:, 1:nFilledLev, iKS), SSqrReal,&
            & denseDesc%iAtomStart)
        localisation = pipekMezey%getLocalisation(eigvecsReal(:,1:nFilledLev,iKS), SSqrReal,&
            & denseDesc%iAtomStart)
        write(stdOut, "(A, E20.12)") 'Final localisation ', localisation
      end do

      call writeRealEigvecs(env, fdEigvec, runId, neighborList, nNeighbor, denseDesc, iSparseStart,&
          & img2CentCell, species(:nAtom), speciesName, orb, over, parallelKS, tPrintEigvecsTxt,&
          & eigvecsReal, SSqrReal, fileName="localOrbs")
    else

      localisation = 0.0_dp
      do iKS = 1, parallelKS%nLocalKS
        iK = parallelKS%localKS(1, iKS)
        iSpin = parallelKS%localKS(2, iKS)
        nFilledLev = floor(nEl(iSpin) / real( 3 - nSpin, dp))
        localisation = localisation + pipekMezey%getLocalisation( &
            & eigvecsCplx(:,:nFilledLev,iKS), SSqrCplx, over, kpoint(:,iK), neighborList,&
            & nNeighbor, iCellVec, cellVec, denseDesc%iAtomStart, iSparseStart, img2CentCell)
      end do
      write(stdOut, "(A, E20.12)") 'Original localisation', localisation

      ! actual localisation calls
      do iKS = 1, parallelKS%nLocalKS
        iK = parallelKS%localKS(1, iKS)
        iSpin = parallelKS%localKS(2, iKS)
        nFilledLev = floor(nEl(iSpin) / real( 3 - nSpin, dp))
        call pipekMezey%calcCoeffs(eigvecsCplx(:,:nFilledLev,iKS), SSqrCplx, over, kpoint(:,iK),&
            & neighborList, nNeighbor, iCellVec, cellVec, denseDesc%iAtomStart, iSparseStart,&
            & img2CentCell)
      end do

      localisation = 0.0_dp
      do iKS = 1, parallelKS%nLocalKS
        iK = parallelKS%localKS(1, iKS)
        iSpin = parallelKS%localKS(2, iKS)
        nFilledLev = floor(nEl(iSpin) / real( 3 - nSpin, dp))
        localisation = localisation + pipekMezey%getLocalisation( &
            & eigvecsCplx(:,:nFilledLev,iKS), SSqrCplx, over, kpoint(:,iK), neighborList,&
            & nNeighbor, iCellVec, cellVec, denseDesc%iAtomStart, iSparseStart, img2CentCell)
      end do
      write(stdOut, "(A, E20.12)") 'Final localisation', localisation

      call writeCplxEigvecs(env, fdEigvec, runId, neighborList, nNeighbor, cellVec, iCellVec,&
          & denseDesc, iSparseStart, img2CentCell, species, speciesName, orb, kPoint, over,&
          & parallelKS, tPrintEigvecsTxt, eigvecsCplx, SSqrCplx, fileName="localOrbs")

    end if

  end subroutine calcPipekMezeyLocalisation


end module main
