!--------------------------------------------------------------------------------------------------!
!  DFTB+: general package for performing fast atomistic simulations                                !
!  Copyright (C) 2006 - 2020  DFTB+ developers group                                               !
!                                                                                                  !
!  See the LICENSE file for terms of usage and distribution.                                       !
!--------------------------------------------------------------------------------------------------!

#:include 'common.fypp'

!> The main routines for DFTB+
module dftbp_main
#:if WITH_MPI
  use dftbp_mpifx
#:endif
#:if WITH_SCALAPACK
  use dftbp_scalapackfx
  use dftbp_scalafxext
#:endif
#:if WITH_SOCKETS
  use dftbp_ipisocket, only : IpiSocketComm
#:endif
  use dftbp_elecsolvers, only : TElectronicSolver, electronicSolverTypes
  use dftbp_assert
  use dftbp_constants
  use dftbp_globalenv
  use dftbp_environment
  use dftbp_densedescr
  use dftbp_inputdata
  use dftbp_nonscc
  use dftbp_eigenvects
  use dftbp_repulsive
  use dftbp_etemp
  use dftbp_populations
  use dftbp_densitymatrix
  use dftbp_forces
  use dftbp_stress
  use dftbp_scc
  use dftbp_sccinit
  use dftbp_onsitecorrection
  use dftbp_externalcharges
  use dftbp_periodic
  use dftbp_mixer
  use dftbp_geoopt
  use dftbp_numderivs2
  use dftbp_spin
  use dftbp_dftbplusu
  use dftbp_mdcommon
  use dftbp_energies
  use dftbp_potentials
  use dftbp_orbitalequiv
  use dftbp_parser
  use dftbp_sparse2dense
#:if not WITH_SCALAPACK
  use dftbp_blasroutines, only : symm, hemm
#:endif
  use dftbp_hsdutils
  use dftbp_charmanip
  use dftbp_shift
  use dftbp_spinorbit
  use dftbp_angmomentum
  use dftbp_elecconstraints
  use dftbp_pmlocalisation, only : TPipekMezey
  use dftbp_linresp
  use dftbp_mainio
  use dftbp_commontypes
  use dftbp_dispersions, only : DispersionIface
  use dftbp_xmlf90
  use dftbp_thirdorder, only : ThirdOrder
  use dftbp_rangeseparated, only : RangeSepFunc
  use dftbp_simplealgebra
  use dftbp_message
  use dftbp_repcont
  use dftbp_halogenx
  use dftbp_xlbomd
  use dftbp_slakocont
  use dftbp_linkedlist
  use dftbp_lapackroutines
  use dftbp_mdcommon
  use dftbp_mdintegrator
  use dftbp_tempprofile
  use dftbp_elstatpot, only : TElStatPotentials
  use dftbp_elstattypes, only : elstatTypes
  use dftbp_forcetypes, only : forceTypes
  use dftbp_initprogram, only : TRefExtPot
  use dftbp_qdepextpotproxy, only : TQDepExtPotProxy
  use dftbp_taggedoutput, only : TTaggedWriter
  use dftbp_plumed, only : TPlumedCalc, TPlumedCalc_final
#:if WITH_TRANSPORT
  use libnegf_vars, only : TTransPar
  use negf_int
  use poisson_init
#:endif
  use dftbp_transportio

  implicit none
  private

  public :: runDftbPlus
  public :: processGeometry

  !> O(N^2) density matrix creation
  logical, parameter :: tDensON2 = .false.

  !> Should further output be appended to detailed.out?
  logical, parameter :: tAppendDetailedOut = .false.


contains

  !> The main DFTB program itself
  subroutine runDftbPlus(env)
    use dftbp_initprogram

    !> Environment settings
    type(TEnvironment), intent(inout) :: env

    !> Geometry steps so far
    integer :: iGeoStep

    !> Lattice geometry steps so far
    integer :: iLatGeoStep

    !> Do we have the final geometry?
    logical :: tGeomEnd

    !> do we take an optimization step on the lattice or the internal coordinates if optimizing both
    !> in a periodic geometry
    logical :: tCoordStep

    !> if scc/geometry driver should be stopped
    logical :: tStopScc, tStopDriver

    !> locality measure for the wavefunction
    real(dp) :: localisation

    !> flag to write out geometries (and charge data if scc) when moving atoms about - in the case
    !> of drivers like conjugate gradient/steepest descent the geometries are written anyway
    logical :: tWriteRestart

    !> lattice vectors returned by the optimizer
    real(dp) :: constrLatDerivs(9)

    !> MD instantaneous thermal energy
    real(dp) :: tempIon

    !> Whether charges should be written
    logical :: tWriteCharges

    logical :: tExitGeoOpt


    call initGeoOptParameters(tCoordOpt, nGeoSteps, tGeomEnd, tCoordStep, tStopDriver, iGeoStep,&
        & iLatGeoStep)

    ! If the geometry is periodic, need to update lattice information in geometry loop
    tLatticeChanged = tPeriodic

    ! As first geometry iteration, require updates for coordinates in dependent routines
    tCoordsChanged = .true.

    ! Main geometry loop
    geoOpt: do iGeoStep = 0, nGeoSteps
      tWriteRestart = env%tGlobalMaster&
          & .and. needsRestartWriting(tGeoOpt, tMd, iGeoStep, nGeoSteps, restartFreq)
      call printGeoStepInfo(tCoordOpt, tLatOpt, iLatGeoStep, iGeoStep)
      call processGeometry(env, iGeoStep, iLatGeoStep, tWriteRestart, tStopDriver, tStopScc,&
          & tExitGeoOpt)
      if (tExitGeoOpt) then
        exit geoOpt
      end if
      call postProcessDerivs(derivs, conAtom, conVec, tLatOpt, totalLatDeriv, extLatDerivs,&
          & normOrigLatVec, tLatOptFixAng, tLatOptFixLen, tLatOptIsotropic, constrLatDerivs)
      call printMaxForces(derivs, constrLatDerivs, tCoordOpt, tLatOpt, indMovedAtom)
    #:if WITH_SOCKETS
      if (tSocket) then
        call sendEnergyAndForces(env, socket, energy, TS, derivs, totalStress, cellVol)
      end if
    #:endif
      tWriteCharges = tWriteRestart .and. tMulliken .and. tSccCalc .and. .not. tDerivs&
          & .and. maxSccIter > 1
      if (tWriteCharges) then
        call writeCharges(fCharges, tWriteChrgAscii, orb, qInput, qBlockIn, qiBlockIn, deltaRhoIn)
      end if

      if (tForces) then
        call getNextGeometry(env, iGeoStep, tWriteRestart, constrLatDerivs, tCoordStep, tGeomEnd,&
            & tStopDriver, iLatGeoStep, tempIon, tExitGeoOpt)
        if (tExitGeoOpt) then
          exit geoOpt
        end if
      end if

      if (tWriteDetailedOut .and. tMd) then
        call writeDetailedOut4(fdDetailedOut, energy, tempIon)
      end if

      if (tGeomEnd) then
        call env%globalTimer%stopTimer(globalTimers%postSCC)
        exit geoOpt
      end if

      tStopDriver = tStopScc .or. tStopDriver .or. hasStopFile(fStopDriver)
      if (tStopDriver) then
        call env%globalTimer%stopTimer(globalTimers%postSCC)
        exit geoOpt
      end if
      call env%globalTimer%stopTimer(globalTimers%postSCC)
    end do geoOpt

    call env%globalTimer%startTimer(globalTimers%postGeoOpt)

  #:if WITH_SOCKETS
    if (tSocket .and. env%tGlobalMaster) then
      call socket%shutdown()
    end if
  #:endif

    if (allocated(plumedCalc)) then
      call TPlumedCalc_final(plumedCalc)
    end if

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
      call writeHessianOut(hessianOut, pDynMatrix)
    else
      nullify(pDynMatrix)
    end if

    if (tWriteShifts) then
      call writeShifts(fShifts, orb, potential%intShell)
    endif

  #:if WITH_TRANSPORT
    if (tContCalc) then
      ! Note: shift and charges are saved in QM representation (not UD)
      call writeContactShifts(transpar%contacts(transpar%taskContInd)%output, orb, &
          & potential%intShell, qOutput, Ef)
    end if

    if (tLocalCurrents) then
      call local_currents(env, parallelKS%localKS, ham, over, neighbourList, nNeighbourSK,&
          & cutOff%skCutoff, denseDesc%iAtomStart, iSparseStart, img2CentCell, iCellVec, cellVec,&
          & rCellVec, orb, kPoint, kWeight, coord0Fold, species0, speciesName, mu, lCurrArray)
    end if

    if (tTunn) then
      call calc_current(env, parallelKS%localKS, ham, over, neighbourList%iNeighbour, nNeighbourSK,&
          & densedesc%iAtomStart, iSparseStart, img2CentCell, iCellVec, cellVec, orb, kPoint,&
          & kWeight, tunneling, current, ldos, leadCurrents, writeTunn, tWriteLDOS,&
          & regionLabelLDOS, mu)
    end if

  #:endif

    if (allocated(pipekMezey)) then
      ! NOTE: the canonical DFTB ground state orbitals are over-written after this point
      if (withMpi) then
        call error("Pipek-Mezey localisation does not yet work with MPI")
      end if
      if (nSpin > 2) then
        call error("Pipek-Mezey localisation not implemented for non-colinear DFTB")
      end if
      call calcPipekMezeyLocalisation(env, pipekMezey, tPrintEigvecsTxt, nEl, filling, over,&
          & kPoint, neighbourList, nNeighbourSk, denseDesc, iSparseStart, img2CentCell, iCellVec,&
          & cellVec, runId, orb, species, speciesName, parallelKS, localisation, eigvecsReal,&
          & SSqrReal, eigvecsCplx, SSqrCplx, tHelical, coord)
    end if

    if (tWriteAutotest) then
      if (tPeriodic) then
        cellVol = abs(determinant33(latVec))
        energy%EGibbs = energy%EMermin + extPressure * cellVol
      end if
      call writeAutotestTag(autotestTag, tPeriodic, cellVol, tMulliken, qOutput, derivs,&
          & chrgForces, excitedDerivs, tStress, totalStress, pDynMatrix, energy, extPressure,&
          & coord0, tLocalise, localisation, esp, taggedWriter, tunneling, ldos, tDefinedFreeE,&
          & lCurrArray)
    end if
    if (tWriteResultsTag) then
      call writeResultsTag(resultsTag, energy, derivs, chrgForces, electronicSolver, tStress,&
          & totalStress, pDynMatrix, tPeriodic, cellVol, tMulliken, qOutput, q0, taggedWriter,&
          & tDefinedFreeE)
    end if
    if (tWriteDetailedXML) then
      call writeDetailedXml(runId, speciesName, species0, pCoord0Out, tPeriodic, latVec, tRealHS,&
          & nKPoint, nSpin, size(eigen, dim=1), nOrb, kPoint, kWeight, filling, occNatural)
    end if

    call env%globalTimer%stopTimer(globalTimers%postGeoOpt)

  #:if WITH_TRANSPORT
    if (tPoisson) then
      call poiss_destroy()
    end if
    if (electronicSolver%iSolver == electronicSolverTypes%GF) then
      call negf_destroy()
    end if
  #:endif

  end subroutine runDftbPlus


  !> Process current geometry
  subroutine processGeometry(env, iGeoStep, iLatGeoStep, tWriteRestart, tStopDriver, tStopScc,&
      & tExitGeoOpt)
    use dftbp_initprogram

    !> Environment settings
    type(TEnvironment), intent(inout) :: env

    !> Current geometry step
    integer, intent(in) :: iGeoStep

    !> Current lattice step
    integer, intent(in) :: iLatGeoStep

    !> flag to write out geometries (and charge data if scc)
    logical, intent(in) :: tWriteRestart

    !> if scc/geometry driver should be stopped
    logical, intent(inout) :: tStopDriver

    !> if scc driver should be stopped
    logical, intent(out) :: tStopScc

    !> Whether main code should exit the geometry optimisation loop
    logical, intent(out) :: tExitGeoOpt

    ! Charge error in the last iterations
    real(dp) :: sccErrorQ, diffElec

    ! Loop variables
    integer :: iSccIter

    ! energy in previous scc cycles
    real(dp) :: Eold

    ! whether scc converged
    logical :: tConverged

    ! Whether scc restart info should be written in current iteration
    logical :: tWriteSccRestart

    ! Charge difference
    real(dp), allocatable :: dQ(:,:,:)

    ! loop index
    integer :: iSpin

    call env%globalTimer%startTimer(globalTimers%preSccInit)

    if (allocated(qDepExtPot)) then
      allocate(dQ(orb%mShell, nAtom, nSpin))
    end if

    call electronicSolver%reset()
    tExitGeoOpt = .false.

    if (tMD .and. tWriteRestart) then
      call writeMdOut1(fdMd, mdOut, iGeoStep, pMDIntegrator)
    end if

    if (tLatticeChanged) then
      call handleLatticeChange(latVec, sccCalc, tStress, extPressure, cutOff%mCutOff, dispersion,&
          & recVec, invLatVec, cellVol, recCellVol, extLatDerivs, cellVec, rCellVec)
    end if

    if (tCoordsChanged) then
      call handleCoordinateChange(env, coord0, latVec, invLatVec, species0, cutOff, orb,&
          & tPeriodic, tHelical, sccCalc, dispersion, thirdOrd, rangeSep, img2CentCell, iCellVec,&
          & neighbourList, nAllAtom, coord0Fold, coord, species, rCellVec, nNeighbourSk,&
          & nNeighbourRep, nNeighbourLC, ham, over, H0, rhoPrim, iRhoPrim, iHam, ERhoPrim,&
          & iSparseStart, tPoisson)
    end if

    #:if WITH_TRANSPORT
      if (tNegf) then
        call initNegfStuff(denseDesc, transpar, ginfo, neighbourList, nNeighbourSK, img2CentCell,&
            & orb)
      end if
    #:endif

    if (tSccCalc) then
      call reset(pChrgMixer, nMixElements)
    end if

    if (electronicSolver%isElsiSolver .and. .not. tLargeDenseMatrices) then
      call electronicSolver%elsi%updateGeometry(env, neighbourList, nNeighbourSK,&
          & denseDesc%iAtomStart, iSparseStart, img2CentCell)
    end if

    call env%globalTimer%startTimer(globalTimers%sparseH0S)
    call buildH0(env, H0, skHamCont, atomEigVal, coord, nNeighbourSk, neighbourList%iNeighbour,&
        & species, iSparseStart, orb)
    call buildS(env, over, skOverCont, coord, nNeighbourSk, neighbourList%iNeighbour, species,&
        & iSparseStart, orb)
    call env%globalTimer%stopTimer(globalTimers%sparseH0S)

    if (tSetFillingTemp) then
      call getTemperature(temperatureProfile, tempElec)
    end if

    call electronicSolver%updateElectronicTemp(tempElec)

    call calcRepulsiveEnergy(coord, species, img2CentCell, nNeighbourRep, neighbourList,&
        & pRepCont, energy%atomRep, energy%ERep, iAtInCentralRegion)

    if (tDispersion) then
      call calcDispersionEnergy(dispersion, energy%atomDisp, energy%Edisp, iAtInCentralRegion)
    end if

    if (allocated(halogenXCorrection)) then
      call halogenXCorrection%getEnergies(energy%atomHalogenX, coord, species, neighbourList,&
          & img2CentCell)
      energy%EHalogenX = sum(energy%atomHalogenX(iAtInCentralRegion))
    end if

    call resetExternalPotentials(refExtPot, potential)

    if (tReadShifts) then
      call readShifts(fShifts, orb, nAtom, nSpin, potential%extShell)
    end if

    call setUpExternalElectricField(tEField, tTDEField, tPeriodic, EFieldStrength,&
        & EFieldVector, EFieldOmega, EFieldPhase, neighbourList, nNeighbourSk, iCellVec,&
        & img2CentCell, cellVec, deltaT, iGeoStep, coord0Fold, coord, potential%extAtom(:,1),&
        & potential%extGrad, EField, absEField)

    call mergeExternalPotentials(orb, species, potential)

    ! For non-scc calculations with transport only, jump out of geometry loop
    if (electronicSolver%iSolver == electronicSolverTypes%OnlyTransport) then
      if (tWriteDetailedOut) then
        call openDetailedOut(fdDetailedOut, userOut, tAppendDetailedOut, iGeoStep, 1)
      end if
      ! We need to define hamltonian by adding the potential
      call getSccHamiltonian(H0, over, nNeighbourSK, neighbourList, species, orb, iSparseStart,&
          & img2CentCell, potential, ham, iHam)
      tExitGeoOpt = .true.
      return
    end if

    if (electronicSolver%iSolver == electronicSolverTypes%pexsi) then
      call electronicSolver%elsi%initPexsiDeltaVRanges(tSccCalc, potential)
    end if

    call initSccLoop(tSccCalc, xlbomdIntegrator, minSccIter, maxSccIter, sccTol, tConverged, tNegf)

    call env%globalTimer%stopTimer(globalTimers%preSccInit)

    call env%globalTimer%startTimer(globalTimers%scc)

    lpSCC: do iSccIter = 1, maxSccIter

      call resetInternalPotentials(tDualSpinOrbit, xi, orb, species, potential)

      if (tSccCalc) then

        call getChargePerShell(qInput, orb, species, chargePerShell)

      #:if WITH_TRANSPORT
        ! Overrides input charges with uploaded contact charges
        if (tUpload) then
          call overrideContactCharges(qInput, chargeUp, transpar)
        end if
      #:endif

        call addChargePotentials(env, sccCalc, qInput, q0, chargePerShell, orb, species,&
            & neighbourList, img2CentCell, spinW, thirdOrd, potential, electrostatics, tPoisson,&
            & tUpload, shiftPerLUp)

        call addBlockChargePotentials(qBlockIn, qiBlockIn, tDftbU, tImHam, species, orb,&
            & nDftbUFunc, UJ, nUJ, iUJ, niUJ, potential)

        if (allocated(onSiteElements) .and. (iSCCIter > 1 .or. tReadChrg)) then
          call addOnsShift(potential%intBlock, potential%iOrbitalBlock, qBlockIn, qiBlockIn, q0,&
              & onSiteElements, species, orb)
        end if

      end if

      ! All potentials are added up into intBlock
      potential%intBlock = potential%intBlock + potential%extBlock

      if (allocated(qDepExtPot)) then
        call getChargePerShell(qInput, orb, species, dQ, qRef=q0)
        call qDepExtPot%addPotential(sum(dQ(:,:,1), dim=1), dQ(:,:,1), orb, species,&
            & potential%intBlock)
      end if

      if (electronicSolver%iSolver == electronicSolverTypes%pexsi .and. tSccCalc) then
        call electronicSolver%elsi%updatePexsiDeltaVRanges(potential)
      end if

      call getSccHamiltonian(H0, over, nNeighbourSK, neighbourList, species, orb, iSparseStart,&
          & img2CentCell, potential, ham, iHam)

      if (tWriteRealHS .or. tWriteHS .and. any(electronicSolver%iSolver ==&
          & [electronicSolverTypes%qr, electronicSolverTypes%divideandconquer,&
          & electronicSolverTypes%relativelyrobust, electronicSolverTypes%magma_gvd])) then
        call writeHSAndStop(env, tWriteHS, tWriteRealHS, tRealHS, over, neighbourList,&
            & nNeighbourSK, denseDesc%iAtomStart, iSparseStart, img2CentCell, kPoint, iCellVec,&
            & cellVec, ham, iHam)
      end if

      call convertToUpDownRepr(ham, iHam)

      call getDensity(env, iSccIter, denseDesc, ham, over, neighbourList, nNeighbourSk,&
          & iSparseStart, img2CentCell, iCellVec, cellVec, kPoint, kWeight, orb, tHelical, coord,&
          & species, electronicSolver, tRealHS, tSpinSharedEf, tSpinOrbit, tDualSpinOrbit,&
          & tFillKSep, tFixEf, tMulliken, iDistribFn, tempElec, nEl, parallelKS, Ef, mu, energy,&
          & rangeSep, eigen, filling, rhoPrim, Eband, TS, E0, iHam, xi, orbitalL, HSqrReal,&
          & SSqrReal, eigvecsReal, iRhoPrim, HSqrCplx, SSqrCplx, eigvecsCplx, rhoSqrReal,&
          & deltaRhoInSqr, deltaRhoOutSqr, qOutput, nNeighbourLC, tLargeDenseMatrices)

      !> For rangeseparated calculations deduct atomic charges from deltaRho
      if (tRangeSep) then
        select case(nSpin)
        case(2)
          do iSpin = 1, 2
            call denseSubtractDensityOfAtoms(q0, denseDesc%iAtomStart, deltaRhoOutSqr, iSpin)
          end do
        case(1)
          call denseSubtractDensityOfAtoms(q0, denseDesc%iAtomStart, deltaRhoOutSqr)
        case default
          call error("Range separation not implemented for non-colinear spin")
        end select
      end if

      if (tWriteBandDat) then
        call writeBandOut(bandOut, eigen, filling, kWeight)
      end if

      if (tMulliken) then
        call getMullikenPopulation(rhoPrim, over, orb, neighbourList, nNeighbourSk, img2CentCell,&
            & iSparseStart, qOutput, iRhoPrim=iRhoPrim, qBlock=qBlockOut, qiBlock=qiBlockOut)
      end if

      #:if WITH_TRANSPORT
        ! Override charges with uploaded contact charges
        if (tUpload) then
          call overrideContactCharges(qOutput, chargeUp, transpar)
        end if
      #:endif

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
            & neighbourList, img2CentCell, spinW, thirdOrd, potential, electrostatics,&
            & tPoissonTwice, tUpload, shiftPerLUp)

        call addBlockChargePotentials(qBlockOut, qiBlockOut, tDftbU, tImHam, species, orb,&
            & nDftbUFunc, UJ, nUJ, iUJ, niUJ, potential)

        if (allocated(onSiteElements)) then
          call addOnsShift(potential%intBlock, potential%iOrbitalBlock, qBlockOut, qiBlockOut,&
              & q0, onSiteElements, species, orb)
        end if

        potential%intBlock = potential%intBlock + potential%extBlock
      end if

      if (allocated(qDepExtPot)) then
        call getChargePerShell(qOutput, orb, species, dQ, qRef=q0)
        call qDepExtPot%addPotential(sum(dQ(:,:,1), dim=1), dQ(:,:,1), orb, species,&
            & potential%intBlock)
      end if

      call getEnergies(sccCalc, qOutput, q0, chargePerShell, species, tExtField, tXlbomd,&
          & tDftbU, tDualSpinOrbit, rhoPrim, H0, orb, neighbourList, nNeighbourSk, img2CentCell,&
          & iSparseStart, cellVol, extPressure, TS, potential, energy, thirdOrd, rangeSep,&
          & qDepExtPot, qBlockOut, qiBlockOut, nDftbUFunc, UJ, nUJ, iUJ, niUJ, xi,&
          & iAtInCentralRegion, tFixEf, Ef, onSiteElements)

      tStopScc = hasStopFile(fStopScc)

      ! Mix charges Input/Output
      if (tSccCalc) then
        if(.not. tRangeSep) then
          call getNextInputCharges(env, pChrgMixer, qOutput, qOutRed, orb, nIneqOrb, iEqOrbitals,&
              & iGeoStep, iSccIter, minSccIter, maxSccIter, sccTol, tStopScc, tMixBlockCharges,&
              & tReadChrg, qInput, qInpRed, sccErrorQ, tConverged, qBlockOut, iEqBlockDftbU,&
              & qBlockIn, qiBlockOut, iEqBlockDftbULS, species0, nUJ, iUJ, niUJ, qiBlockIn,&
              & iEqBlockOnSite, iEqBlockOnSiteLS)
        else
          call getNextInputDensity(SSqrReal, over, neighbourList, nNeighbourSK,&
              & denseDesc%iAtomStart, iSparseStart, img2CentCell, pChrgMixer, qOutput, orb,&
              & tHelical, species, coord, iGeoStep, iSccIter, minSccIter, maxSccIter, sccTol,&
              & tStopScc, tReadChrg, q0, qInput, sccErrorQ, tConverged, deltaRhoOut, deltaRhoIn,&
              & deltaRhoDiff, qBlockIn, qBlockOut)
        end if

        call getSccInfo(iSccIter, energy%Eelec, Eold, diffElec)
        if (tNegf) then
          call printSccHeader()
        end if
        call printSccInfo(tDftbU, iSccIter, energy%Eelec, diffElec, sccErrorQ)

        if (tNegf) then
          call printBlankLine()
        end if

        tWriteSccRestart = env%tGlobalMaster .and. &
            & needsSccRestartWriting(restartFreq, iGeoStep, iSccIter, minSccIter, maxSccIter,&
            & tMd, tGeoOpt, tDerivs, tConverged, tReadChrg, tStopScc)
        if (tWriteSccRestart) then
          call writeCharges(fCharges, tWriteChrgAscii, orb, qInput, qBlockIn, qiBlockIn, deltaRhoIn)
        end if
      end if

      if (tWriteDetailedOut) then
        call openDetailedOut(fdDetailedOut, userOut, tAppendDetailedOut, iGeoStep, iSccIter)
        call writeDetailedOut1(fdDetailedOut, iDistribFn, nGeoSteps, iGeoStep, tMD, tDerivs,&
            & tCoordOpt, tLatOpt, iLatGeoStep, iSccIter, energy, diffElec, sccErrorQ, indMovedAtom,&
            & pCoord0Out, q0, qInput, qOutput, eigen, filling, orb, species, tDFTBU,&
            & tImHam.or.tSpinOrbit, tPrintMulliken, orbitalL, qBlockOut, Ef, Eband, TS, E0,&
            & extPressure, cellVol, tAtomicEnergy, tDispersion, tEField, tPeriodic, nSpin, tSpin,&
            & tSpinOrbit, tSccCalc, allocated(onSiteElements), tNegf, invLatVec, kPoint,&
            & iAtInCentralRegion, electronicSolver, tDefinedFreeE, allocated(halogenXCorrection),&
            & tRangeSep, allocated(thirdOrd))
      end if

      if (tConverged .or. tStopScc) then
        exit lpSCC
      end if

    end do lpSCC

    call env%globalTimer%stopTimer(globalTimers%scc)

    #:if WITH_TRANSPORT
      if (tPoisson) then
        call poiss_savepotential()
      end if
    #:endif

    call env%globalTimer%startTimer(globalTimers%postSCC)

    if (tLinResp) then
      if (withMpi) then
        call error("Linear response calc. does not work with MPI yet")
      end if
      call ensureLinRespConditions(t3rd, tRealHS, tPeriodic, tCasidaForces)
      call calculateLinRespExcitations(env, lresp, parallelKS, sccCalc, qOutput, q0, over,&
          & eigvecsReal, eigen(:,1,:), filling(:,1,:), coord, species, speciesName, orb, tHelical,&
          & skHamCont, skOverCont, autotestTag, taggedWriter, runId, neighbourList, nNeighbourSK,&
          & denseDesc, iSparseStart, img2CentCell, tWriteAutotest, tCasidaForces, tLinRespZVect,&
          & tPrintExcitedEigvecs, tPrintEigvecsTxt, nonSccDeriv, energy, energiesCasida, SSqrReal,&
          & rhoSqrReal, excitedDerivs, occNatural)
    end if

    if (tXlbomd) then
      call getXlbomdCharges(xlbomdIntegrator, qOutRed, pChrgMixer, orb, nIneqOrb, iEqOrbitals,&
          & qInput, qInpRed, iEqBlockDftbU, qBlockIn, species0, nUJ, iUJ, niUJ, iEqBlockDftbuLs,&
          & qiBlockIn, iEqBlockOnSite, iEqBlockOnSiteLS)
    end if

    if (tDipole) then
      call getDipoleMoment(qOutput, q0, coord, dipoleMoment, iAtInCentralRegion)
    #:call DEBUG_CODE
      call checkDipoleViaHellmannFeynman(rhoPrim, q0, coord0, over, orb, neighbourList,&
          & nNeighbourSk, species, iSparseStart, img2CentCell)
    #:endcall DEBUG_CODE
    end if

    call env%globalTimer%startTimer(globalTimers%eigvecWriting)

    if (tPrintEigVecs) then
      call writeEigenvectors(env, runId, neighbourList, nNeighbourSk, cellVec, iCellVec, denseDesc,&
          & iSparseStart, img2CentCell, species, speciesName, orb, kPoint, over, parallelKS,&
          & tPrintEigvecsTxt, eigvecsReal, SSqrReal, eigvecsCplx, SSqrCplx)
    end if

    if (tProjEigenvecs) then
      call writeProjectedEigenvectors(env, regionLabels, eigen, neighbourList, nNeighbourSk,&
          & cellVec, iCellVec, denseDesc, iSparseStart, img2CentCell, orb, over, kPoint, kWeight,&
          & iOrbRegion, parallelKS, eigvecsReal, SSqrReal, eigvecsCplx, SSqrCplx)
    end if
    call env%globalTimer%stopTimer(globalTimers%eigvecWriting)

    ! MD geometry files are written only later, once velocities for the current geometry are known
    if (tGeoOpt .and. tWriteRestart) then
      call writeCurrentGeometry(geoOutFile, pCoord0Out, tLatOpt, tMd, tAppendGeo, tFracCoord,&
          & tPeriodic, tHelical, tPrintMulliken, species0, speciesName, latVec,&
          & iGeoStep, iLatGeoStep, nSpin, qOutput, velocities)
    end if

    call printEnergies(energy, TS, electronicSolver, tDefinedFreeE)

    if (tForces) then
      call env%globalTimer%startTimer(globalTimers%forceCalc)
      call env%globalTimer%startTimer(globalTimers%energyDensityMatrix)
      call getEnergyWeightedDensity(env, electronicSolver, denseDesc, forceType, filling, eigen,&
          & kPoint, kWeight, neighbourList, nNeighbourSk, orb, iSparseStart, img2CentCell,&
          & iCellVec, cellVec, tRealHS, ham, over, parallelKS, tHelical, species, coord, iSccIter,&
          & mu, ERhoPrim, eigvecsReal, SSqrReal, eigvecsCplx, SSqrCplx)
      call env%globalTimer%stopTimer(globalTimers%energyDensityMatrix)
      call getGradients(env, sccCalc, tExtField, tXlbomd, nonSccDeriv, EField, rhoPrim, ERhoPrim,&
          & qOutput, q0, skHamCont, skOverCont, pRepCont, neighbourList, nNeighbourSk,&
          & nNeighbourRep, species, img2CentCell, iSparseStart, orb, potential, coord, derivs,&
          & iRhoPrim, thirdOrd, qDepExtPot, chrgForces, dispersion, rangeSep, SSqrReal, over,&
          & denseDesc, deltaRhoOutSqr, tPoisson, halogenXCorrection, tHelical, coord0)

      if (tCasidaForces) then
        derivs(:,:) = derivs + excitedDerivs
      end if
      call env%globalTimer%stopTimer(globalTimers%forceCalc)

      call updateDerivsByPlumed(env, plumedCalc, nAtom, iGeoStep, derivs, energy%EMermin, coord0,&
          & mass, tPeriodic, latVec)

      if (tStress) then
        if (tHelical) then
          call env%globalTimer%startTimer(globalTimers%stressCalc)
          call getHelicalBCDerivs(env, latvec(:2,1), sccCalc, thirdOrd, tExtField, nonSccDeriv,&
              & rhoPrim, ERhoPrim, qOutput, q0, skHamCont, skOverCont, pRepCont, neighbourList,&
              & nNeighbourSk, nNeighbourRep, species, img2CentCell, iSparseStart, orb, potential,&
              & coord, iRhoPrim, dispersion, halogenXCorrection)
          call env%globalTimer%stopTimer(globalTimers%stressCalc)
        else
          call env%globalTimer%startTimer(globalTimers%stressCalc)
          call getStress(env, sccCalc, thirdOrd, tExtField, nonSccDeriv, rhoPrim, ERhoPrim,&
              & qOutput, q0, skHamCont, skOverCont, pRepCont, neighbourList, nNeighbourSk,&
              & nNeighbourRep, species, img2CentCell, iSparseStart, orb, potential, coord, latVec,&
              & invLatVec, cellVol, coord0, totalStress, totalLatDeriv, intPressure, iRhoPrim,&
              & dispersion, halogenXCorrection)
          call env%globalTimer%stopTimer(globalTimers%stressCalc)
          call printVolume(cellVol)
        end if

        ! MD case includes the atomic kinetic energy contribution, so print that later
        if (.not. (tMD .or. tHelical)) then
          call printPressureAndFreeEnergy(extPressure, intPressure, energy%EGibbs)
        end if
      end if

    end if

    if (tWriteDetailedOut) then
      call writeDetailedOut2(fdDetailedOut, tSccCalc, tConverged, tXlbomd, tLinResp, tGeoOpt,&
          & tMD, tPrintForces, tStress, tPeriodic, energy, totalStress, totalLatDeriv, derivs, &
          & chrgForces, indMovedAtom, cellVol, intPressure, geoOutFile, iAtInCentralRegion)
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

  end subroutine processGeometry


  subroutine postprocessDerivs(derivs, conAtom, conVec, tLatOpt, totalLatDerivs,&
      & extLatDerivs, normLatVecs, tLatOptFixAng, tLatOptFixLen, tLatOptIsotropic,&
      & constrLatDerivs)

    !> On input energy derivatives, on exit resulting projected derivatives
    real(dp), intent(inout), allocatable :: derivs(:,:)

    !> Atoms being constrained
    integer, allocatable, intent(in) :: conAtom(:)

    !> Vector to project out forces
    real(dp), allocatable, intent(in) :: conVec(:,:)

    !> Whether lattice optimisation is on
    logical, intent(in) :: tLatOpt

    !> Derivative of total energy with respect to lattice vectors
    real(dp) :: totalLatDerivs(:,:)

    !> derivative of cell volume wrt to lattice vectors, needed for pV term
    real(dp), intent(in) :: extLatDerivs(:,:)

    !> Unit normals parallel to lattice vectors
    real(dp), intent(in) :: normLatVecs(:,:)

    !> Are the angles of the lattice being fixed during optimisation?
    logical, intent(in) :: tLatOptFixAng

    !> Are the magnitude of the lattice vectors fixed
    logical, intent(in) :: tLatOptFixLen(:)

    !> Is the optimisation isotropic
    logical, intent(in) :: tLatOptIsotropic

    !> Lattice vectors returned by the optimizer
    real(dp), intent(out) :: constrLatDerivs(:)

    if (allocated(conAtom)) then
      call constrainForces(conAtom, conVec, derivs)
    end if

    if (tLatOpt) then
      ! Only include the extLatDerivs contribution if not MD, as the barostat would otherwise
      ! take care of this, hence add it here rather than to totalLatDeriv itself
      call constrainLatticeDerivs(totalLatDerivs + extLatDerivs, normLatVecs, tLatOptFixAng,&
          & tLatOptFixLen, tLatOptIsotropic, constrLatDerivs)
    end if

  end subroutine postprocessDerivs


  subroutine getNextGeometry(env, iGeoStep, tWriteRestart, constrLatDerivs, tCoordStep, tGeomEnd,&
      & tStopDriver, iLatGeoStep, tempIon, tExitGeoOpt)
    use dftbp_initprogram

    !> Environment settings
    type(TEnvironment), intent(inout) :: env

    !> Current geometry step
    integer, intent(in) :: iGeoStep

    !> flag to write out geometries (and charge data if scc)
    logical, intent(in) :: tWriteRestart

    !> lattice vectors returned by the optimizer
    real(dp), intent(in) :: constrLatDerivs(:)

    !> do we take an optimization step on the lattice or the internal coordinates if optimizing both
    !> in a periodic geometry
    logical, intent(inout) :: tCoordStep

    !> Do we have the final geometry?
    logical, intent(inout) :: tGeomEnd

    !> If geometry driver should be stopped
    logical, intent(inout) :: tStopDriver

    !> Current lattice step
    integer, intent(inout) :: iLatGeoStep

    !> MD instantaneous thermal energy
    real(dp), intent(out) :: tempIon

    !> Whether geometry optimisation should be stop
    logical, intent(out) :: tExitGeoOpt


    !> Difference between last calculated and new geometry.
    real(dp) :: diffGeo

    !> Has this completed?
    logical :: tCoordEnd

    ! initially assume that coordinates and lattice vectors won't be updated
    tCoordsChanged = .false.
    tLatticeChanged = .false.

    tExitGeoOpt = .false.

    if (tDerivs) then
      call getNextDerivStep(derivDriver, derivs, indMovedAtom, coord0, tGeomEnd)
      if (tGeomEnd) then
        call env%globalTimer%stopTimer(globalTimers%postSCC)
        tExitGeoOpt = .true.
        return
      end if
      tCoordsChanged = .true.
    else if (tGeoOpt) then
      tCoordsChanged = .true.
      if (tCoordStep) then
        call getNextCoordinateOptStep(pGeoCoordOpt, energy, derivs, indMovedAtom, coord0, diffGeo,&
            & tCoordEnd, .not. tCasidaForces)
        if (.not. tLatOpt) then
          tGeomEnd = tCoordEnd
        end if
        if (.not. tGeomEnd .and. tCoordEnd .and. diffGeo < tolSameDist) then
          tCoordStep = .false.
        end if
      else
        call getNextLatticeOptStep(pGeoLatOpt, energy, constrLatDerivs, origLatVec, tLatOptFixAng,&
            & tLatOptFixLen, tLatOptIsotropic, indMovedAtom, latVec, coord0, diffGeo, tGeomEnd)
        iLatGeoStep = iLatGeoStep + 1
        tLatticeChanged = .true.
        if (.not. tGeomEnd .and. tCoordOpt) then
          tCoordStep = .true.
          call reset(pGeoCoordOpt, reshape(coord0(:, indMovedAtom), [nMovedCoord]))
        end if
      end if
      if (tGeomEnd .and. diffGeo < tolSameDist) then
        call env%globalTimer%stopTimer(globalTimers%postSCC)
        tExitGeoOpt = .true.
        return
      end if
    else if (tMD) then
      ! New MD coordinates saved in a temporary variable, as writeCurrentGeometry() below
      ! needs the old ones to write out consistent geometries and velocities.
      newCoords(:,:) = coord0
      call getNextMdStep(pMdIntegrator, pMdFrame, temperatureProfile, derivs, movedMass, mass,&
          & cellVol, invLatVec, species0, indMovedAtom, tStress, tBarostat, energy, newCoords,&
          & latVec, intPressure, totalStress, totalLatDeriv, velocities, tempIon)
      tCoordsChanged = .true.
      tLatticeChanged = tBarostat
      call printMdInfo(tSetFillingTemp, tEField, tPeriodic, tempElec, absEField, tempIon,&
          & intPressure, extPressure, energy)
      if (tWriteRestart) then
        if (tPeriodic) then
          cellVol = abs(determinant33(latVec))
          energy%EGibbs = energy%EMermin + extPressure * cellVol
        end if
        call writeMdOut2(fdMd, tStress, tPeriodic, tBarostat, tLinResp, tEField, tFixEf,&
            & tPrintMulliken, energy, energiesCasida, latVec, cellVol, intPressure, extPressure,&
            & tempIon, absEField, qOutput, q0, dipoleMoment)
        call writeCurrentGeometry(geoOutFile, pCoord0Out, .false., .true., .true., tFracCoord,&
            & tPeriodic, tHelical, tPrintMulliken, species0, speciesName, latVec, iGeoStep,&
            & iLatGeoStep, nSpin, qOutput, velocities)
      end if
      coord0(:,:) = newCoords
      if (tWriteDetailedOut) then
        call writeDetailedOut3(fdDetailedOut, tPrintForces, tSetFillingTemp, tPeriodic, tStress,&
            & totalStress, totalLatDeriv, energy, tempElec, extPressure, intPressure, tempIon)
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

  end subroutine getNextGeometry



  !> Initialises some parameters before geometry loop starts.
  subroutine initGeoOptParameters(tCoordOpt, nGeoSteps, tGeomEnd, tCoordStep, tStopDriver, iGeoStep&
      &, iLatGeoStep)

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


  !> Does the operations that are necessary after a lattice vector update
  subroutine handleLatticeChange(latVecs, sccCalc, tStress, extPressure, mCutOff, dispersion,&
      & recVecs, recVecs2p, cellVol, recCellVol, extLatDerivs, cellVecs, rCellVecs)

    !> lattice vectors
    real(dp), intent(in) :: latVecs(:,:)

    !> Module variables
    type(TScc), allocatable, intent(inout) :: sccCalc

    !> evaluate stress
    logical, intent(in) :: tStress

    !> External presure
    real(dp), intent(in) :: extPressure

    !> Maximum distance for interactions
    real(dp), intent(inout) :: mCutOff

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
      mCutOff = max(mCutOff, sccCalc%getCutOff())
    end if
    if (allocated(dispersion)) then
      call dispersion%updateLatVecs(latVecs)
      mCutOff = max(mCutOff, dispersion%getRCutOff())
    end if
    call getCellTranslations(cellVecs, rCellVecs, latVecs, recVecs2p, mCutOff)

  end subroutine handleLatticeChange


  !> Does the operations that are necessary after atomic coordinates change
  subroutine handleCoordinateChange(env, coord0, latVec, invLatVec, species0, cutOff, orb,&
      & tPeriodic, tHelical, sccCalc, dispersion, thirdOrd, rangeSep, img2CentCell, iCellVec,&
      & neighbourList, nAllAtom, coord0Fold, coord, species, rCellVec, nNeighbourSK, nNeighbourRep,&
      & nNeighbourLC, ham, over, H0, rhoPrim, iRhoPrim, iHam, ERhoPrim, iSparseStart, tPoisson)

    use dftbp_initprogram, only : OCutoffs

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

    !> Longest cut-off distances that neighbour maps are generated for
    type(OCutoffs), intent(in) :: cutOff

    !> Atomic orbital information
    type(TOrbitals), intent(in) :: orb

    !> Is the geometry periodic
    logical, intent(in) :: tPeriodic

    !> Is the geometry helical
    logical, intent(in) :: tHelical

    !> SCC module internal variables
    type(TScc), allocatable, intent(inout) :: sccCalc

    !> Dispersion interactions
    class(DispersionIface), allocatable, intent(inout) :: dispersion

    !> Third order SCC interactions
    type(ThirdOrder), allocatable, intent(inout) :: thirdOrd

    !> Range separation contributions
    type(RangeSepFunc), allocatable, intent(inout) :: rangeSep

    !> Image atoms to their equivalent in the central cell
    integer, allocatable, intent(inout) :: img2CentCell(:)

    !> Index for which unit cell an atom is in
    integer, allocatable, intent(inout) :: iCellVec(:)

    !> List of neighbouring atoms
    type(TNeighbourList), intent(inout) :: neighbourList

    !> Total number of atoms including images
    integer, intent(out) :: nAllAtom

    !> Central cell atomic coordinates, folded inside the central cell
    real(dp), intent(out) :: coord0Fold(:,:)

    !> Coordinates of all atoms including images
    real(dp), allocatable, intent(inout) :: coord(:,:)

    !> Species of all atoms including images
    integer, allocatable, intent(inout) :: species(:)

    !> Vectors to units cells in absolute units
    real(dp), allocatable, intent(in) :: rCellVec(:,:)

    !> Number of neighbours of each real atom
    integer, intent(out) :: nNeighbourSK(:)

    !> Number of neighbours of each real atom close enough for repulsive interactions
    integer, intent(out) :: nNeighbourRep(:)

    !> Number of neighbours for each of the atoms for the exchange contributions in the long range
    !> functional
    integer, intent(inout), allocatable :: nNeighbourLC(:)

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

    !> Transport variables
    logical, intent(in) :: tPoisson

    !> Total size of orbitals in the sparse data structures, where the decay of the overlap sets the
    !> sparsity pattern
    integer :: sparseSize

    coord0Fold(:,:) = coord0
    if (tPeriodic .or. tHelical) then
      call foldCoordToUnitCell(coord0Fold, latVec, invLatVec)
    end if

    if (tHelical) then
      call updateNeighbourListAndSpecies(coord, species, img2CentCell, iCellVec, neighbourList,&
          & nAllAtom, coord0Fold, species0, cutoff%mCutoff, rCellVec, latVec=latVec)
    else
      call updateNeighbourListAndSpecies(coord, species, img2CentCell, iCellVec, neighbourList,&
          & nAllAtom, coord0Fold, species0, cutoff%mCutOff, rCellVec)
    end if

    call getNrOfNeighboursForAll(nNeighbourSK, neighbourList, cutoff%skCutOff)

    call getSparseDescriptor(neighbourList%iNeighbour, nNeighbourSK, img2CentCell, orb,&
        & iSparseStart, sparseSize)
    call reallocateSparseArrays(sparseSize, ham, over, H0, rhoPrim, iHam, iRhoPrim, ERhoPrim)

    ! count neighbours for repulsive interactions between atoms
    call getNrOfNeighboursForAll(nNeighbourRep, neighbourList, cutoff%repCutOff)

    if (allocated(nNeighbourLC)) then
      ! count neighbours for repulsive interactions between atoms
      call getNrOfNeighboursForAll(nNeighbourLC, neighbourList, cutoff%lcCutOff)
    end if

    ! Notify various modules about coordinate changes
  #:if WITH_TRANSPORT
    if (tPoisson) then
      !! TODO: poiss_updcoords pass coord0 and not coord0Fold because the
      !! folding can mess up the contact position. Could we have the supercell
      !! centered on the input atomic structure?
      call poiss_updcoords(coord0)
    end if
  #:endif

    if (allocated(sccCalc)) then
      call sccCalc%updateCoords(env, coord, species, neighbourList)
    end if

    if (allocated(dispersion)) then
      call dispersion%updateCoords(neighbourList, img2CentCell, coord, species0)
    end if
    if (allocated(thirdOrd)) then
      call thirdOrd%updateCoords(neighbourList, species)
    end if
    if (allocated(rangeSep)) then
       call rangeSep%updateCoords(coord0)
    end if


  end subroutine handleCoordinateChange


#:if WITH_TRANSPORT

  !> Initialise transport
  subroutine initNegfStuff(denseDescr, transpar, ginfo, neighbourList, nNeighbourSK, img2CentCell,&
      & orb)

    !> Dense matrix descriptor
    type(TDenseDescr), intent(in) :: denseDescr

    !> Transport settings
    type(TTransPar), intent(in) :: transpar

    !> libNEGF data
    type(TNEGFInfo), intent(in) :: ginfo

    !> Atomic orbital information
    type(TOrbitals), intent(in) :: orb

    !> Image atoms to their equivalent in the central cell
    integer, intent(in) :: img2CentCell(:)

    !> List of neighbouring atoms
    type(TNeighbourList), intent(in) :: neighbourList

    !> Number of neighbours of each real atom
    integer, intent(in) :: nNeighbourSK(:)

    ! known issue about the PLs: We need an automatic partitioning
    call negf_init_csr(denseDescr%iAtomStart, neighbourList%iNeighbour, nNeighbourSK, img2CentCell,&
        & orb)

    call negf_init_str(denseDescr, transpar, ginfo%greendens, neighbourList%iNeighbour,&
        & nNeighbourSK, img2CentCell)

    call negf_init_dephasing(ginfo%tundos)  !? why tundos

  end subroutine initNegfStuff

#:endif


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
  subroutine calcRepulsiveEnergy(coord, species, img2CentCell, nNeighbourRep, neighbourList,&
      & pRepCont, Eatom, Etotal, iAtInCentralRegion)

    !> All atomic coordinates
    real(dp), intent(in) :: coord(:,:)

    !> All atoms chemical species
    integer, intent(in) :: species(:)

    !> Image atom indices to central cell atoms
    integer, intent(in) :: img2CentCell(:)

    !> Number of neighbours for each atom within the repulsive distance
    integer, intent(in) :: nNeighbourRep(:)

    !> List of neighbours for each atom
    type(TNeighbourList), intent(in) :: neighbourList

    !> Repulsive interaction data
    type(ORepCont), intent(in) :: pRepCont

    !> Energy for each atom
    real(dp), intent(out) :: Eatom(:)

    !> Total energy
    real(dp), intent(out) :: Etotal

    !> atoms in the central cell (or device region if transport)
    integer, intent(in) :: iAtInCentralRegion(:)

    call getERep(Eatom, coord, nNeighbourRep, neighbourList%iNeighbour, species, pRepCont,&
        & img2CentCell)
    Etotal = sum(Eatom(iAtInCentralRegion))

  end subroutine calcRepulsiveEnergy


  !> Calculates dispersion energy for current geometry.
  subroutine calcDispersionEnergy(dispersion, Eatom, Etotal, iAtInCentralRegion)

    !> dispersion interactions
    class(DispersionIface), intent(inout) :: dispersion

    !> energy per atom
    real(dp), intent(out) :: Eatom(:)

    !> total energy
    real(dp), intent(out) :: Etotal

    !> atoms in the central cell (or device region if transport)
    integer, intent(in) :: iAtInCentralRegion(:)

    call dispersion%getEnergies(Eatom)
    Etotal = sum(Eatom(iAtInCentralRegion))

  end subroutine calcDispersionEnergy


  !> Sets the external potential components to zero
  subroutine resetExternalPotentials(refExtPot, potential)

    !> Reference external potential (usually set via API)
    type(TRefExtPot), intent(in) :: refExtPot

    !> Potential contributions
    type(TPotentials), intent(inout) :: potential

    if (allocated(refExtPot%atomPot)) then
      potential%extAtom(:,:) = refExtPot%atomPot
    else
      potential%extAtom(:,:) = 0.0_dp
    end if
    if (allocated(refExtPot%shellPot)) then
      potential%extShell(:,:,:) = refExtPot%shellPot
    else
      potential%extShell(:,:,:) = 0.0_dp
    end if
    potential%extBlock(:,:,:,:) = 0.0_dp
    if (allocated(refExtPot%potGrad)) then
      potential%extGrad(:,:) = refExtPot%potGrad
    else
      potential%extGrad(:,:) = 0.0_dp
    end if

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
      & EFieldVector, EFieldOmega, EFieldPhase, neighbourList, nNeighbourSK, iCellVec,&
      & img2CentCell, cellVec, deltaT, iGeoStep, coord0Fold, coord, extAtomPot, extPotGrad, EField,&
      & absEField)

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

    !> Atomic neighbours
    type(TNeighbourList), intent(in) :: neighbourList

    !> Number of neighbours for each atom
    integer, intent(in) :: nNeighbourSK(:)

    !> Index for unit cells
    integer, intent(in) :: iCellVec(:)

    !> Image atom to central cell atom number
    integer, intent(in) :: img2CentCell(:)

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

    !> Potentials on atomic sites
    real(dp), intent(inout) :: extAtomPot(:)

    !> Gradient of potential on atomic sites with respect of nucleus positions. Shape: (3, nAtom)
    real(dp), intent(inout) :: extPotGrad(:,:)

    !> Resulting electric field
    real(dp), intent(out) :: EField(:)

    !> Magnitude of the field
    real(dp), intent(out) :: absEField

    integer :: nAtom
    integer :: iAt1, iAt2, iNeigh
    character(lc) :: tmpStr

    if (.not. tEField) then
      EField(:) = 0.0_dp
      absEField = 0.0_dp
      return
    end if

    nAtom = size(nNeighbourSK)

    Efield(:) = EFieldStrength * EfieldVector
    if (tTimeDepEField) then
      Efield(:) = Efield * sin(EfieldOmega * deltaT * real(iGeoStep + EfieldPhase, dp))
    end if
    absEfield = sqrt(sum(Efield**2))
    if (tPeriodic) then
      do iAt1 = 1, nAtom
        do iNeigh = 1, nNeighbourSK(iAt1)
          iAt2 = neighbourList%iNeighbour(iNeigh, iAt1)
          ! overlap between atom in central cell and non-central cell
          if (iCellVec(iAt2) /= 0) then
            ! component of electric field projects onto vector between cells
            if (abs(dot_product(cellVec(:, iCellVec(iAt2)), EfieldVector)) > epsilon(1.0_dp)) then
              write(tmpStr, "(A, I0, A, I0, A)") 'Interaction between atoms ', iAt1, ' and ',&
                  & img2CentCell(iAt2), ' crosses the saw-tooth discontinuity in the electric&
                  & field.'
              call error(tmpStr)
            end if
          end if
        end do
      end do
      do iAt1 = 1, nAtom
        extAtomPot(iAt1) = extAtomPot(iAt1) + dot_product(coord0Fold(:, iAt1), Efield)
      end do
    else
      do iAt1 = 1, nAtom
        extAtomPot(iAt1) = extAtomPot(iAt1) + dot_product(coord(:, iAt1), Efield)
      end do
    end if
    extPotGrad(:,:) = extPotGrad + spread(EField, 2, nAtom)

  end subroutine setUpExternalElectricField


  !> Initialise basic variables before the scc loop.
  subroutine initSccLoop(tSccCalc, xlbomdIntegrator, minSccIter, maxSccIter, sccTol, tConverged,&
      & tNegf)

    !> Is this an SCC calculation?
    logical, intent(in) :: tSccCalc

    !> Details for extended Lagrange integrator (of used)
    type(Xlbomd), allocatable, intent(inout) :: xlbomdIntegrator

    !> Minimum number of SCC cycles that can be used
    integer, intent(inout) :: minSccIter

    !> Maximum number of SCC cycles
    integer, intent(inout) :: maxSccIter

    !> Tolerance for SCC convergence
    real(dp), intent(inout) :: sccTol

    !> Has SCC convergence been achieved?
    logical, intent(out) :: tConverged

    !> Is this a transport calculation?
    logical, intent(in) :: tNegf

    if (allocated(xlbomdIntegrator)) then
      call xlbomdIntegrator%getSCCParameters(minSccIter, maxSccIter, sccTol)
    end if

    tConverged = (.not. tSccCalc)

    if (tSccCalc .and. .not. tNegf) then
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


#:if WITH_TRANSPORT

  !> Replace charges with those from the stored contact values
  subroutine overrideContactCharges(qInput, chargeUp, transpar)
    !> input charges
    real(dp), intent(inout) :: qInput(:,:,:)

    !> uploaded charges
    real(dp), intent(in) :: chargeUp(:,:,:)

    !> Transport parameters
    type(TTransPar), intent(in) :: transpar

    integer :: ii, iStart, iEnd

    do ii = 1, transpar%ncont
      iStart = transpar%contacts(ii)%idxrange(1)
      iEnd = transpar%contacts(ii)%idxrange(2)
      qInput(:,iStart:iEnd,:) = chargeUp(:,iStart:iEnd,:)
    end do

  end subroutine overrideContactCharges

#:endif


  !> Add potentials comming from point charges.
  subroutine addChargePotentials(env, sccCalc, qInput, q0, chargePerShell, orb, species,&
      & neighbourList, img2CentCell, spinW, thirdOrd, potential, electrostatics, tPoisson,&
      & tUpload, shiftPerLUp)

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
    type(TNeighbourList), intent(in) :: neighbourList

    !> map from image atom to real atoms
    integer, intent(in) :: img2CentCell(:)

    !> spin constants
    real(dp), intent(in), allocatable :: spinW(:,:,:)

    !> third order SCC interactions
    type(ThirdOrder), allocatable, intent(inout) :: thirdOrd

    !> Potentials acting
    type(TPotentials), intent(inout) :: potential

    !> electrostatic solver (poisson or gamma-functional)
    integer, intent(in) :: electrostatics

    !> whether Poisson is solved (used with tPoissonTwice)
    logical, intent(in) :: tPoisson

    !> whether contacts are uploaded
    logical, intent(in) :: tUpload

    !> uploded potential per shell per atom
    real(dp), allocatable, intent(in) :: shiftPerLUp(:,:)

    ! local variables
    real(dp), allocatable :: atomPot(:,:)
    real(dp), allocatable :: shellPot(:,:,:)
    real(dp), allocatable, save :: shellPotBk(:,:)
    integer, pointer :: pSpecies0(:)
    integer :: nAtom, nSpin

    nAtom = size(qInput, dim=2)
    nSpin = size(qInput, dim=3)
    pSpecies0 => species(1:nAtom)

    allocate(atomPot(nAtom, nSpin))
    allocate(shellPot(orb%mShell, nAtom, nSpin))

    call sccCalc%updateCharges(env, qInput, q0, orb, species)

    select case(electrostatics)

    case(elstatTypes%gammaFunc)

      call sccCalc%updateShifts(env, orb, species, neighbourList%iNeighbour, img2CentCell)
      call sccCalc%getShiftPerAtom(atomPot(:,1))
      call sccCalc%getShiftPerL(shellPot(:,:,1))

    case(elstatTypes%poisson)

    #:if WITH_TRANSPORT
      ! NOTE: charge-magnetization representation is used
      !       iSpin=1 stores total charge
      ! Logic of calls order:
      ! shiftPerLUp      is 0.0 on the device region,
      ! poiss_getshift() updates only the device region
      if (tPoisson) then
        if (tUpload) then
          shellPot(:,:,1) = shiftPerLUp
        else
          ! Potentials for non-existing angular momenta must be 0 for later summations
          shellPot(:,:,1) = 0.0_dp
        end if
        call poiss_updcharges(qInput(:,:,1), q0(:,:,1))
        call poiss_getshift(shellPot(:,:,1))
        if (.not.allocated(shellPotBk)) then
          allocate(shellPotBk(orb%mShell, nAtom))
        end if
        shellPotBk = shellPot(:,:,1)
      else
        shellPot(:,:,1) = shellPotBk
      end if
      atomPot(:,:) = 0.0_dp
      call sccCalc%setShiftPerAtom(atomPot(:,1))
      call sccCalc%setShiftPerL(shellPot(:,:,1))
    #:else
      call error("poisson solver used without transport modules")
    #:endif

    end select

    potential%intAtom(:,1) = potential%intAtom(:,1) + atomPot(:,1)
    potential%intShell(:,:,1) = potential%intShell(:,:,1) + shellPot(:,:,1)

    if (allocated(thirdOrd)) then
      call thirdOrd%updateCharges(pSpecies0, neighbourList, qInput, q0, img2CentCell, orb)
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
  subroutine addBlockChargePotentials(qBlockIn, qiBlockIn, tDftbU, tImHam, species, orb, nDftbUFunc&
      &, UJ, nUJ, iUJ, niUJ, potential)

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
  subroutine getSccHamiltonian(H0, over, nNeighbourSK, neighbourList, species, orb, iSparseStart,&
      & img2CentCell, potential, ham, iHam)

    !> non-SCC hamitonian (sparse)
    real(dp), intent(in) :: H0(:)

    !> overlap (sparse)
    real(dp), intent(in) :: over(:)

    !> Number of atomic neighbours
    integer, intent(in) :: nNeighbourSK(:)

    !> list of atomic neighbours
    type(TNeighbourList), intent(in) :: neighbourList

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
    call add_shift(ham, over, nNeighbourSK, neighbourList%iNeighbour, species, orb, iSparseStart,&
        & nAtom, img2CentCell, potential%intBlock)

    if (allocated(iHam)) then
      iHam(:,:) = 0.0_dp
      call add_shift(iHam, over, nNeighbourSK, neighbourList%iNeighbour, species, orb,&
          & iSparseStart, nAtom, img2CentCell, potential%iorbitalBlock)
    end if

  end subroutine getSccHamiltonian


  !> Transform the hamiltonian from QM to UD representation
  !> Hack due to not using Pauli-type structure for diagonalisation
  !> For collinear spin, qm2ud will produce the right potential:
  !> (Vq, uB*Bz*\sigma_z) -> (Vq + uB*Bz*\sigma_z, Vq - uB*Bz*\sigma_z)
  !> For non-collinear spin-orbit, all blocks are multiplied by 1/2:
  !> (Vq/2, uL* Lx*\sigma_x/2, uL* Ly*\sigma_y/2, uL* Lz*\sigma_z/2)
  subroutine convertToUpDownRepr(Ham, iHam)
    real(dp), intent(inout) :: Ham(:,:)
    real(dp), intent(inout), allocatable :: iHam(:,:)

    integer :: nSpinBlocks

    nSpinBlocks = size(ham, dim=2)

    if (nSpinBlocks > 1) then
      ham = 2.0_dp * ham
      if (allocated(iHam)) then
        iHam = 2.0_dp * iHam
      end if
    end if

    if (nSpinBlocks /= 4) then
      call qm2ud(ham)
      if (allocated(iHam)) then
        call qm2ud(iHam)
      end if
    end if

  end subroutine convertToUpDownRepr


  !> Returns the sparse density matrix.
  !>
  !> All operations (e.g. non-dual spin orbit coupling), which need access to full (unpacked)
  !> Hamiltonian or the full (unpacked) density matrix, must also invoked from within this routine,
  !> as those unpacked quantities do not exist elsewhere.
  !>
  subroutine getDensity(env, iScc, denseDesc, ham, over, neighbourList, nNeighbourSK, iSparseStart,&
      & img2CentCell, iCellVec, cellVec, kPoint, kWeight, orb, tHelical, coord, species,&
      & electronicSolver, tRealHS, tSpinSharedEf, tSpinOrbit, tDualSpinOrbit, tFillKSep, tFixEf,&
      & tMulliken, iDistribFn, tempElec, nEl, parallelKS, Ef, mu, energy, rangeSep, eigen, filling,&
      & rhoPrim, Eband, TS, E0, iHam, xi, orbitalL, HSqrReal, SSqrReal, eigvecsReal, iRhoPrim,&
      & HSqrCplx, SSqrCplx, eigvecsCplx, rhoSqrReal, deltaRhoInSqr, deltaRhoOutSqr, qOutput,&
      & nNeighbourLC, tLargeDenseMatrices)

    !> Environment settings
    type(TEnvironment), intent(inout) :: env

    !> SCC iteration counter (needed by GF)
    integer, intent(in) :: iSCC

    !> Dense matrix descriptor
    type(TDenseDescr), intent(in) :: denseDesc

    !> hamiltonian in sparse storage
    real(dp), intent(in) :: ham(:,:)

    !> sparse overlap matrix
    real(dp), intent(in) :: over(:)

    !> list of neighbours for each atom
    type(TNeighbourList), intent(in) :: neighbourList

    !> Number of neighbours for each of the atoms
    integer, intent(in) :: nNeighbourSK(:)

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

    !> Is the geometry helical
    logical, intent(in) :: tHelical

    !> Coordinates of all atoms including images
    real(dp), allocatable, intent(inout) :: coord(:,:)

    !> species of all atoms in the system
    integer, intent(in) :: species(:)

    !> Electronic solver information
    type(TElectronicSolver), intent(inout) :: electronicSolver

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

    !> Electrochemical potentials (contact, spin)
    real(dp), allocatable, intent(in) :: mu(:,:)

    !> Energy contributions and total
    type(TEnergies), intent(inout) :: energy

    !> Data for rangeseparated calculation
    type(RangeSepFunc), allocatable, intent(inout) :: rangeSep

    !> eigenvalues (level, kpoint, spin)
    real(dp), intent(out) :: eigen(:,:,:)

    !> occupations (level, kpoint, spin)
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

    !> Change in density matrix during last SCC iteration
    real(dp), pointer, intent(inout) :: deltaRhoInSqr(:,:,:)

    !> Change in density matrix after SCC step
    real(dp), pointer, intent(inout) :: deltaRhoOutSqr(:,:,:)

    !> Output electrons
    real(dp), intent(inout) :: qOutput(:,:,:)

    !> Number of neighbours for each of the atoms for the exchange contributions in the long range
    !> functional
    integer, intent(in), allocatable :: nNeighbourLC(:)

    !> Are dense matrices for H, S, etc. being used
    logical, intent(in) :: tLargeDenseMatrices

    integer :: nSpin, iKS, iSp, iK, nAtom
    complex(dp), allocatable :: rhoSqrCplx(:,:)
    logical :: tImHam
    real(dp), allocatable :: rVecTemp(:)

    nSpin = size(ham, dim=2)
    tImHam = allocated(iRhoPrim)

    select case (electronicSolver%iSolver)

    case (electronicSolverTypes%GF)

      call env%globalTimer%startTimer(globalTimers%densityMatrix)
    #:if WITH_TRANSPORT
      call calcdensity_green(iSCC, env, parallelKS%localKS, ham, over, neighbourlist%iNeighbour,&
          & nNeighbourSK, denseDesc%iAtomStart, iSparseStart, img2CentCell, iCellVec, cellVec, orb,&
          & kPoint, kWeight, mu, rhoPrim, Eband, Ef, E0, TS)
    #:else
      call error("Internal error: getDensity : GF-solver although code compiled without transport")
    #:endif
      call ud2qm(rhoPrim)
      call env%globalTimer%stopTimer(globalTimers%densityMatrix)

    case (electronicSolverTypes%onlyTransport)

      call error("OnlyTransport solver cannot calculate the density matrix")

    case(electronicSolverTypes%qr, electronicSolverTypes%divideandconquer,&
        & electronicSolverTypes%relativelyrobust, electronicSolverTypes%elpa,&
        & electronicSolverTypes%magma_gvd)

      call getDensityFromDenseDiag(env, denseDesc, ham, over, neighbourList, nNeighbourSK,&
          & iSparseStart, img2CentCell, iCellVec, cellVec, kPoint, kWeight, orb,&
          & denseDesc%iAtomStart, tHelical, coord, species, electronicSolver, tRealHS,&
          & tSpinSharedEf, tSpinOrbit, tDualSpinOrbit, tFillKSep, tFixEf, tMulliken, iDistribFn,&
          & tempElec, nEl, parallelKS, Ef, energy, rangeSep, eigen, filling, rhoPrim, Eband, TS,&
          & E0, iHam, xi, orbitalL, HSqrReal, SSqrReal, eigvecsReal, iRhoPrim, HSqrCplx, SSqrCplx,&
          & eigvecsCplx, rhoSqrReal, deltaRhoInSqr, deltaRhoOutSqr, qOutput, nNeighbourLC)

    case(electronicSolverTypes%omm, electronicSolverTypes%pexsi, electronicSolverTypes%ntpoly)

      call env%globalTimer%startTimer(globalTimers%densityMatrix)
      call electronicSolver%elsi%getDensity(env, denseDesc, ham, over, neighbourList, nNeighbourSK,&
          & iSparseStart, img2CentCell, iCellVec, cellVec, kPoint, kWeight, tHelical, orb, species,&
          & coord, tRealHS, tSpinSharedEf, tSpinOrbit, tDualSpinOrbit, tMulliken, parallelKS, Ef,&
          & energy, rhoPrim, Eband, TS, iHam, xi, orbitalL, HSqrReal, SSqrReal, iRhoPrim, HSqrCplx,&
          & SSqrCplx)
      call env%globalTimer%stopTimer(globalTimers%densityMatrix)

    end select

  end subroutine getDensity


  !> Returns the density matrix using dense diagonalisation.
  subroutine getDensityFromDenseDiag(env, denseDesc, ham, over, neighbourList, nNeighbourSK,&
      & iSparseStart, img2CentCell, iCellVec, cellVec, kPoint, kWeight, orb, iAtomStart, tHelical, coord,&
      & species, electronicSolver, tRealHS, tSpinSharedEf, tSpinOrbit, tDualSpinOrbit, tFillKSep,&
      & tFixEf, tMulliken, iDistribFn, tempElec, nEl, parallelKS, Ef, energy, rangeSep, eigen,&
      & filling, rhoPrim, Eband, TS, E0, iHam, xi, orbitalL, HSqrReal, SSqrReal, eigvecsReal,&
      & iRhoPrim, HSqrCplx, SSqrCplx, eigvecsCplx, rhoSqrReal, deltaRhoInSqr, deltaRhoOutSqr,&
      & qOutput, nNeighbourLC)

    !> Environment settings
    type(TEnvironment), intent(inout) :: env

    !> Dense matrix descriptor
    type(TDenseDescr), intent(in) :: denseDesc

    !> hamiltonian in sparse storage
    real(dp), intent(in) :: ham(:,:)

    !> sparse overlap matrix
    real(dp), intent(in) :: over(:)

    !> list of neighbours for each atom
    type(TNeighbourList), intent(in) :: neighbourList

    !> Number of neighbours for each of the atoms
    integer, intent(in) :: nNeighbourSK(:)

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

    !>Start of atomic blocks in dense arrays
    integer, allocatable, intent(in) :: iAtomStart(:)

    !> Is the geometry helical
    logical, intent(in) :: tHelical

    !> Coordinates of all atoms including images
    real(dp), allocatable, intent(inout) :: coord(:,:)

    !> species of all atoms in the system
    integer, intent(in) :: species(:)

    !> Electronic solver information
    type(TElectronicSolver), intent(inout) :: electronicSolver

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

    !> Data for rangeseparated calculation
    type(RangeSepFunc), allocatable, intent(inout) :: rangeSep

    !> eigenvalues (level, kpoint, spin)
    real(dp), intent(out) :: eigen(:,:,:)

    !> occupations (level, kpoint, spin)
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

    !> Change in density matrix during last rangesep SCC cycle
    real(dp), pointer, intent(in) :: deltaRhoInSqr(:,:,:)

    !> Change in density matrix during this SCC step for rangesep
    real(dp), pointer, intent(inout) :: deltaRhoOutSqr(:,:,:)

    !> Output electrons
    real(dp), intent(inout) :: qOutput(:,:,:)

    !> Number of neighbours for each of the atoms for the exchange contributions in the long range
    !> functional
    integer, intent(in), allocatable :: nNeighbourLC(:)

    integer :: nSpin

    nSpin = size(ham, dim=2)
    call env%globalTimer%startTimer(globalTimers%diagonalization)
    if (nSpin /= 4) then
      if (tRealHS) then
        call buildAndDiagDenseRealHam(env, denseDesc, ham, over, species, neighbourList,&
            & nNeighbourSK, iSparseStart, img2CentCell, orb, iAtomStart, tHelical, coord,&
            & electronicSolver, parallelKS, rangeSep, deltaRhoInSqr, qOutput, nNeighbourLC,&
            & HSqrReal, SSqrReal, eigVecsReal, eigen(:,1,:))
      else
        call buildAndDiagDenseCplxHam(env, denseDesc, ham, over, kPoint, neighbourList,&
            & nNeighbourSK, iSparseStart, img2CentCell, iCellVec, cellVec, electronicSolver,&
            & parallelKS, tHelical, orb, species, coord, HSqrCplx, SSqrCplx, eigVecsCplx, eigen)
      end if
    else
      call buildAndDiagDensePauliHam(env, denseDesc, ham, over, kPoint, neighbourList,&
          & nNeighbourSK, iSparseStart, img2CentCell, iCellVec, cellVec, orb, electronicSolver,&
          & parallelKS, eigen(:,:,1), HSqrCplx, SSqrCplx, eigVecsCplx, iHam, xi, species)
    end if
    call env%globalTimer%stopTimer(globalTimers%diagonalization)

    call getFillingsAndBandEnergies(eigen, nEl, nSpin, tempElec, kWeight, tSpinSharedEf,&
        & tFillKSep, tFixEf, iDistribFn, Ef, filling, Eband, TS, E0)

    call env%globalTimer%startTimer(globalTimers%densityMatrix)
    if (nSpin /= 4) then
      if (tRealHS) then
        call getDensityFromRealEigvecs(env, denseDesc, filling(:,1,:), neighbourList, nNeighbourSK,&
            & iSparseStart, img2CentCell, orb, species, iAtomStart, coord, tHelical, eigVecsReal,&
            & parallelKS, rhoPrim, SSqrReal, rhoSqrReal, deltaRhoOutSqr)
      else
        call getDensityFromCplxEigvecs(env, denseDesc, filling, kPoint, kWeight, neighbourList,&
            & nNeighbourSK, iSparseStart, img2CentCell, iCellVec, cellVec, orb,&
            & parallelKS, tHelical, species, coord, eigvecsCplx, rhoPrim, SSqrCplx)
      end if
      call ud2qm(rhoPrim)
    else
      ! Pauli structure of eigenvectors
      filling(:,:,1) = 2.0_dp * filling(:,:,1)
      call getDensityFromPauliEigvecs(env, denseDesc, tRealHS, tSpinOrbit, tDualSpinOrbit,&
          & tMulliken, kPoint, kWeight, filling(:,:,1), neighbourList, nNeighbourSK, orb,&
          & iSparseStart, img2CentCell, iCellVec, cellVec, species, parallelKS, eigVecsCplx,&
          & SSqrCplx, energy, rhoPrim, xi, orbitalL, iRhoPrim)
      filling(:,:,1) = 0.5_dp * filling(:,:,1)
    end if
    call env%globalTimer%stopTimer(globalTimers%densityMatrix)

  end subroutine getDensityFromDenseDiag


  !> Builds and diagonalises dense Hamiltonians.
  subroutine buildAndDiagDenseRealHam(env, denseDesc, ham, over, species, neighbourList,&
      & nNeighbourSK, iSparseStart, img2CentCell, orb, iAtomStart, tHelical, coord,&
      & electronicSolver, parallelKS, rangeSep, deltaRhoInSqr, qOutput, nNeighbourLC, HSqrReal,&
      & SSqrReal, eigvecsReal, eigen)

    !> Environment settings
    type(TEnvironment), intent(inout) :: env

    !> Dense matrix descriptor
    type(TDenseDescr), intent(in) :: denseDesc

    !> hamiltonian in sparse storage
    real(dp), intent(in) :: ham(:,:)

    !> sparse overlap matrix
    real(dp), intent(in) :: over(:)

    !> list of neighbours for each atom
    type(TNeighbourList), intent(in) :: neighbourList

    !> Number of neighbours for each of the atoms
    integer, intent(in) :: nNeighbourSK(:)

    !> Index array for the start of atomic blocks in sparse arrays
    integer, intent(in) :: iSparseStart(:,:)

    !> map from image atoms to the original unique atom
    integer, intent(in) :: img2CentCell(:)

    !> Atomic orbital information
    type(TOrbitals), intent(in) :: orb

    !> species of all atoms in the system
    integer, intent(in) :: species(:)

    !>Start of atomic blocks in dense arrays
    integer, allocatable, intent(in) :: iAtomStart(:)

    !> Is the geometry helical
    logical, intent(in) :: tHelical

    !> Coordinates of all atoms including images
    real(dp), allocatable, intent(inout) :: coord(:,:)

    !> Electronic solver information
    type(TElectronicSolver), intent(inout) :: electronicSolver

    !> K-points and spins to be handled
    type(TParallelKS), intent(in) :: parallelKS

    !>Data for rangeseparated calcualtion
    type(RangeSepFunc), allocatable, intent(inout) :: rangeSep

    !> Change in density matrix during last rangesep SCC cycle
    real(dp), pointer, intent(in) :: deltaRhoInSqr(:,:,:)

    !> Output electrons
    real(dp), intent(inout) :: qOutput(:,:,:)

    !> Number of neighbours for each of the atoms for the exchange contributions in the long range
    !> functional
    integer, intent(in), allocatable :: nNeighbourLC(:)

    !> dense hamitonian matrix
    real(dp), intent(out) :: HSqrReal(:,:)

    !> dense overlap matrix
    real(dp), intent(out) :: SSqrReal(:,:)

    !> Eigenvectors on eixt
    real(dp), intent(out) :: eigvecsReal(:,:,:)

    !> eigenvalues
    real(dp), intent(out) :: eigen(:,:)

    integer :: iKS, iSpin, ii

    eigen(:,:) = 0.0_dp
    do iKS = 1, parallelKS%nLocalKS
      iSpin = parallelKS%localKS(2, iKS)
    #:if WITH_SCALAPACK
      call env%globalTimer%startTimer(globalTimers%sparseToDense)
      if (tHelical) then
        call unpackHSHelicalRealBlacs(env%blacs, ham(:,iSpin), neighbourList%iNeighbour,&
            & nNeighbourSK, iSparseStart, img2CentCell, orb, species, coord, denseDesc, HSqrReal)
        if (.not. electronicSolver%hasCholesky(1)) then
          call unpackHSHelicalRealBlacs(env%blacs, over, neighbourList%iNeighbour,&
              & nNeighbourSK, iSparseStart, img2CentCell, orb, species, coord, denseDesc, SSqrReal)
        end if
      else
        call unpackHSRealBlacs(env%blacs, ham(:,iSpin), neighbourList%iNeighbour, nNeighbourSK,&
            & iSparseStart, img2CentCell, denseDesc, HSqrReal)
        if (.not. electronicSolver%hasCholesky(1)) then
          call unpackHSRealBlacs(env%blacs, over, neighbourList%iNeighbour, nNeighbourSK,&
              & iSparseStart, img2CentCell, denseDesc, SSqrReal)
        end if
      end if
      call env%globalTimer%stopTimer(globalTimers%sparseToDense)
      call diagDenseMtxBlacs(electronicSolver, 1, 'V', denseDesc%blacsOrbSqr, HSqrReal, SSqrReal,&
          & eigen(:,iSpin), eigvecsReal(:,:,iKS))
    #:else
      call env%globalTimer%startTimer(globalTimers%sparseToDense)
      if (tHelical) then
        call unpackHS(HSqrReal, ham(:,iSpin), neighbourList%iNeighbour, nNeighbourSK,&
            & denseDesc%iAtomStart, iSparseStart, img2CentCell, orb, species, coord)
        call unpackHS(SSqrReal, over, neighbourList%iNeighbour, nNeighbourSK, denseDesc%iAtomStart,&
            & iSparseStart, img2CentCell, orb, species, coord)
      else
        call unpackHS(HSqrReal, ham(:,iSpin), neighbourList%iNeighbour, nNeighbourSK,&
            & denseDesc%iAtomStart, iSparseStart, img2CentCell)
        call unpackHS(SSqrReal, over, neighbourList%iNeighbour, nNeighbourSK, denseDesc%iAtomStart,&
            & iSparseStart, img2CentCell)
      end if
      call env%globalTimer%stopTimer(globalTimers%sparseToDense)

      ! Add rangeseparated contribution
      ! Assumes deltaRhoInSqr only used by rangeseparation
      ! Should this be used elsewhere, need to pass tRangeSep
      if (allocated(rangeSep)) then
        call denseMulliken(deltaRhoInSqr, SSqrReal, denseDesc%iAtomStart, qOutput)
        call rangeSep%addLRHamiltonian(env, deltaRhoInSqr(:,:,iSpin), over,&
            & neighbourList%iNeighbour,  nNeighbourLC, denseDesc%iAtomStart, iSparseStart,&
            & orb, HSqrReal, SSqrReal)
      end if

      call diagDenseMtx(electronicSolver, 'V', HSqrReal, SSqrReal, eigen(:,iSpin))
      eigvecsReal(:,:,iKS) = HSqrReal
    #:endif
    end do

  #:if WITH_SCALAPACK
    ! Distribute all eigenvalues to all nodes via global summation
    call mpifx_allreduceip(env%mpi%interGroupComm, eigen, MPI_SUM)
  #:endif

  end subroutine buildAndDiagDenseRealHam


  !> Builds and diagonalises dense k-point dependent Hamiltonians.
  subroutine buildAndDiagDenseCplxHam(env, denseDesc, ham, over, kPoint, neighbourList,&
      & nNeighbourSK, iSparseStart, img2CentCell, iCellVec, cellVec, electronicSolver, parallelKS,&
      & tHelical, orb, species, coord, HSqrCplx, SSqrCplx, eigvecsCplx, eigen)

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
    type(TNeighbourList), intent(in) :: neighbourList

    !> Number of neighbours for each of the atoms
    integer, intent(in) :: nNeighbourSK(:)

    !> Index array for the start of atomic blocks in sparse arrays
    integer, intent(in) :: iSparseStart(:,:)

    !> map from image atoms to the original unique atom
    integer, intent(in) :: img2CentCell(:)

    !> Index for which unit cell atoms are associated with
    integer, intent(in) :: iCellVec(:)

    !> Vectors (in units of the lattice constants) to cells of the lattice
    real(dp), intent(in) :: cellVec(:,:)

    !> Electronic solver information
    type(TElectronicSolver), intent(inout) :: electronicSolver

    !> K-points and spins to be handled
    type(TParallelKS), intent(in) :: parallelKS

    !> Is the geometry helical
    logical, intent(in) :: tHelical

    !> Atomic orbital information
    type(TOrbitals), intent(in) :: orb

    !> species of all atoms in the system
    integer, intent(in) :: species(:)

    !> atomic coordinates
    real(dp), intent(in) :: coord(:,:)

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
      if (tHelical) then
        call unpackHSHelicalCplxBlacs(env%blacs, ham(:,iSpin), kPoint(:,iK),&
            & neighbourList%iNeighbour, nNeighbourSK, iCellVec, cellVec, iSparseStart,&
            & img2CentCell, orb, species, coord, denseDesc, HSqrCplx)
        if (.not. electronicSolver%hasCholesky(iKS)) then
          call unpackHSHelicalCplxBlacs(env%blacs, over, kPoint(:,iK), neighbourList%iNeighbour,&
              & nNeighbourSK, iCellVec, cellVec, iSparseStart, img2CentCell, orb, species, coord,&
              & denseDesc, SSqrCplx)
        end if
      else
        call unpackHSCplxBlacs(env%blacs, ham(:,iSpin), kPoint(:,iK), neighbourList%iNeighbour,&
            & nNeighbourSK, iCellVec, cellVec, iSparseStart, img2CentCell, denseDesc, HSqrCplx)
        if (.not. electronicSolver%hasCholesky(iKS)) then
          call unpackHSCplxBlacs(env%blacs, over, kPoint(:,iK), neighbourList%iNeighbour,&
              & nNeighbourSK, iCellVec, cellVec, iSparseStart, img2CentCell, denseDesc, SSqrCplx)
        end if
      end if
      call env%globalTimer%stopTimer(globalTimers%sparseToDense)
      call diagDenseMtxBlacs(electronicSolver, iKS, 'V', denseDesc%blacsOrbSqr, HSqrCplx, SSqrCplx,&
          & eigen(:,iK,iSpin), eigvecsCplx(:,:,iKS))
    #:else
      call env%globalTimer%startTimer(globalTimers%sparseToDense)
      if (tHelical) then
        call unpackHS(HSqrCplx, ham(:,iSpin), kPoint(:,iK), neighbourList%iNeighbour, nNeighbourSK,&
            & iCellVec, cellVec, denseDesc%iAtomStart, iSparseStart, img2CentCell, orb, species,&
            & coord)
        call unpackHS(SSqrCplx, over, kPoint(:,iK), neighbourList%iNeighbour, nNeighbourSK,&
            & iCellVec, cellVec, denseDesc%iAtomStart, iSparseStart, img2CentCell, orb, species,&
            & coord)
      else
        call unpackHS(HSqrCplx, ham(:,iSpin), kPoint(:,iK), neighbourList%iNeighbour, nNeighbourSK,&
            & iCellVec, cellVec, denseDesc%iAtomStart, iSparseStart, img2CentCell)
        call unpackHS(SSqrCplx, over, kPoint(:,iK), neighbourList%iNeighbour, nNeighbourSK,&
            & iCellVec, cellVec, denseDesc%iAtomStart, iSparseStart, img2CentCell)
      end if
      call env%globalTimer%stopTimer(globalTimers%sparseToDense)
      call diagDenseMtx(electronicSolver, 'V', HSqrCplx, SSqrCplx, eigen(:,iK,iSpin))
      eigvecsCplx(:,:,iKS) = HSqrCplx
      write(345,*)
      write(345,*)eigen(:,iK,iSpin)
      write(345,*)
    #:endif
    end do
  
  #:if WITH_SCALAPACK
    call mpifx_allreduceip(env%mpi%interGroupComm, eigen, MPI_SUM)
  #:endif

  end subroutine buildAndDiagDenseCplxHam


  !> Builds and diagonalizes Pauli two-component Hamiltonians.
  subroutine buildAndDiagDensePauliHam(env, denseDesc, ham, over, kPoint, neighbourList,&
      & nNeighbourSK, iSparseStart, img2CentCell, iCellVec, cellVec, orb, electronicSolver,&
      & parallelKS, eigen, HSqrCplx, SSqrCplx, eigvecsCplx, iHam, xi, species)

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
    type(TNeighbourList), intent(in) :: neighbourList

    !> Number of neighbours for each of the atoms
    integer, intent(in) :: nNeighbourSK(:)

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

    !> Electronic solver information
    type(TElectronicSolver), intent(inout) :: electronicSolver

    !> K-points and spins to be handled
    type(TParallelKS), intent(in) :: parallelKS

    !> eigenvalues (orbital, kpoint)
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
        call unpackHPauliBlacs(env%blacs, ham, kPoint(:,iK), neighbourList%iNeighbour,&
            & nNeighbourSK, iCellVec, cellVec, iSparseStart, img2CentCell, orb%mOrb, denseDesc,&
            & HSqrCplx, iorig=iHam)
      else
        call unpackHPauliBlacs(env%blacs, ham, kPoint(:,iK), neighbourList%iNeighbour,&
            & nNeighbourSK, iCellVec, cellVec, iSparseStart, img2CentCell, orb%mOrb, denseDesc,&
            & HSqrCplx)
      end if
      if (.not. electronicSolver%hasCholesky(iKS)) then
        call unpackSPauliBlacs(env%blacs, over, kPoint(:,iK), neighbourList%iNeighbour,&
            & nNeighbourSK, iCellVec, cellVec, iSparseStart, img2CentCell, orb%mOrb, denseDesc,&
            & SSqrCplx)
      end if
    #:else
      if (allocated(iHam)) then
        call unpackHPauli(ham, kPoint(:,iK), neighbourList%iNeighbour, nNeighbourSK, iSparseStart,&
            & denseDesc%iAtomStart, img2CentCell, iCellVec, cellVec, HSqrCplx, iHam=iHam)
      else
        call unpackHPauli(ham, kPoint(:,iK), neighbourList%iNeighbour, nNeighbourSK, iSparseStart,&
            & denseDesc%iAtomStart, img2CentCell, iCellVec, cellVec, HSqrCplx)
      end if
      call unpackSPauli(over, kPoint(:,iK), neighbourList%iNeighbour, nNeighbourSK,&
          & denseDesc%iAtomStart, iSparseStart, img2CentCell, iCellVec, cellVec, SSqrCplx)
    #:endif
      if (allocated(xi) .and. .not. allocated(iHam)) then
        call addOnsiteSpinOrbitHam(env, xi, species, orb, denseDesc, HSqrCplx)
      end if
      call env%globalTimer%stopTimer(globalTimers%sparseToDense)
    #:if WITH_SCALAPACK
      call diagDenseMtxBlacs(electronicSolver, iKS, 'V', denseDesc%blacsOrbSqr, HSqrCplx, SSqrCplx,&
          & eigen(:,iK), eigvecsCplx(:,:,iKS))
    #:else
      call diagDenseMtx(electronicSolver, 'V', HSqrCplx, SSqrCplx, eigen(:,iK))
      eigvecsCplx(:,:,iKS) = HSqrCplx
    #:endif
    end do

  #:if WITH_SCALAPACK
    call mpifx_allreduceip(env%mpi%interGroupComm, eigen, MPI_SUM)
  #:endif

  end subroutine buildAndDiagDensePauliHam


  !> Creates sparse density matrix from real eigenvectors.
  subroutine getDensityFromRealEigvecs(env, denseDesc, filling, neighbourList, nNeighbourSK,&
      & iSparseStart, img2CentCell, orb, species, iAtomStart, coord, tHelical, eigvecs, parallelKS,&
      & rhoPrim, work, rhoSqrReal, deltaRhoOutSqr)

    !> Environment settings
    type(TEnvironment), intent(inout) :: env

    !> Dense matrix descriptor
    type(TDenseDescr), intent(in) :: denseDesc

    !> Filling
    real(dp), intent(in) :: filling(:,:)

    !> list of neighbours for each atom
    type(TNeighbourList), intent(in) :: neighbourList

    !> Number of neighbours for each of the atoms
    integer, intent(in) :: nNeighbourSK(:)

    !> Index array for the start of atomic blocks in sparse arrays
    integer, intent(in) :: iSparseStart(:,:)

    !> map from image atoms to the original unique atom
    integer, intent(in) :: img2CentCell(:)

    !> Atomic orbital information
    type(TOrbitals), intent(in) :: orb

    !> species of atoms
    integer, intent(in), optional :: species(:)

    !>Start of atomic blocks in dense arrays
    integer, allocatable, intent(in) :: iAtomStart(:)

    !> K-points and spins to process
    type(TParallelKS), intent(in) :: parallelKS

    !> all coordinates
    real(dp), intent(in) :: coord(:,:)

    !> Is the geometry helical
    logical, intent(in) :: tHelical

    !> eigenvectors
    real(dp), intent(inout) :: eigvecs(:,:,:)

    !> sparse density matrix
    real(dp), intent(out) :: rhoPrim(:,:)

    !> work space array
    real(dp), intent(out) :: work(:,:)

    !> Dense density matrix if needed
    real(dp), intent(inout), allocatable  :: rhoSqrReal(:,:,:)

    !> Change in density matrix during this SCC step for rangesep
    real(dp), pointer, intent(inout) :: deltaRhoOutSqr(:,:,:)

    integer :: iKS, iSpin

    rhoPrim(:,:) = 0.0_dp
    do iKS = 1, parallelKS%nLocalKS
      iSpin = parallelKS%localKS(2, iKS)

    #:if WITH_SCALAPACK
      call makeDensityMtxRealBlacs(env%blacs%orbitalGrid, denseDesc%blacsOrbSqr, filling(:,iSpin),&
          & eigvecs(:,:,iKS), work)
      call env%globalTimer%startTimer(globalTimers%denseToSparse)
      if (tHelical) then
        call packRhoHelicalRealBlacs(env%blacs, denseDesc, work, neighbourList%iNeighbour,&
            & nNeighbourSK, iSparseStart, img2CentCell, orb, species, coord, rhoPrim(:,iSpin))
      else
        call packRhoRealBlacs(env%blacs, denseDesc, work, neighbourList%iNeighbour, nNeighbourSK,&
            & orb%mOrb, iSparseStart, img2CentCell, rhoPrim(:,iSpin))
      end if
      call env%globalTimer%stopTimer(globalTimers%denseToSparse)
    #:else
      !> Either pack density matrix or delta density matrix
      if(.not. associated(deltaRhoOutSqr)) then
        if (tDensON2) then
          call makeDensityMatrix(work, eigvecs(:,:,iKS), filling(:,iSpin),&
              & neighbourlist%iNeighbour, nNeighbourSK, orb, denseDesc%iAtomStart, img2CentCell)
        else
          call makeDensityMatrix(work, eigvecs(:,:,iKS), filling(:,iSpin))
        end if

        call env%globalTimer%startTimer(globalTimers%denseToSparse)
        if (tHelical) then
          call packHS(rhoPrim(:,iSpin), work, neighbourlist%iNeighbour, nNeighbourSK,&
              & denseDesc%iAtomStart, iSparseStart, img2CentCell, orb, species, coord)
        else
          call packHS(rhoPrim(:,iSpin), work, neighbourlist%iNeighbour, nNeighbourSK, orb%mOrb,&
              & denseDesc%iAtomStart, iSparseStart, img2CentCell)
        end if
        call env%globalTimer%stopTimer(globalTimers%denseToSparse)
      else
        ! Rangeseparated case: pack delta density matrix
        call makeDensityMatrix(deltaRhoOutSqr(:,:,iSpin),&
            & eigvecs(:,:,iKS), filling(:,iSpin))
        call env%globalTimer%startTimer(globalTimers%denseToSparse)
        if (tHelical) then
          call packHS(rhoPrim(:,iSpin), deltaRhoOutSqr(:,:,iSpin), neighbourlist%iNeighbour,&
              & nNeighbourSK, denseDesc%iAtomStart, iSparseStart, img2CentCell, orb, species, coord)
        else
          call packHS(rhoPrim(:,iSpin), deltaRhoOutSqr(:,:,iSpin), neighbourlist%iNeighbour,&
              & nNeighbourSK, orb%mOrb, denseDesc%iAtomStart, iSparseStart, img2CentCell)
        end if
        call env%globalTimer%stopTimer(globalTimers%denseToSparse)
      end if
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
  subroutine getDensityFromCplxEigvecs(env, denseDesc, filling, kPoint, kWeight, neighbourList,&
      & nNeighbourSK, iSparseStart, img2CentCell, iCellVec, cellVec, orb, parallelKS, tHelical,&
      & species, coord, eigvecs, rhoPrim, work)

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
    type(TNeighbourList), intent(in) :: neighbourList

    !> Number of neighbours for each of the atoms
    integer, intent(in) :: nNeighbourSK(:)

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

    !> Is the geometry helical
    logical, intent(in) :: tHelical

    !> species for atoms
    integer, intent(in) :: species(:)

    !> all coordinates
    real(dp), intent(in) :: coord(:,:)

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
      call makeDensityMtxCplxBlacs(env%blacs%orbitalGrid, denseDesc%blacsOrbSqr, filling(:,iK&
          &,iSpin), eigvecs(:,:,iKS), work)
      call env%globalTimer%startTimer(globalTimers%denseToSparse)
      if (tHelical) then
        call packRhoHelicalCplxBlacs(env%blacs, denseDesc, work, kPoint(:,iK), kWeight(iK),&
            & neighbourList%iNeighbour, nNeighbourSK, iCellVec, cellVec, iSparseStart,&
            & img2CentCell, orb, species, coord, rhoPrim(:,iSpin))
      else
        call packRhoCplxBlacs(env%blacs, denseDesc, work, kPoint(:,iK), kWeight(iK),&
            & neighbourList%iNeighbour, nNeighbourSK, orb%mOrb, iCellVec, cellVec, iSparseStart,&
            & img2CentCell, rhoPrim(:,iSpin))
      end if
      call env%globalTimer%stopTimer(globalTimers%denseToSparse)
    #:else
      if (tDensON2) then
        call makeDensityMatrix(work, eigvecs(:,:,iKS), filling(:,iK,iSpin),&
            & neighbourlist%iNeighbour, nNeighbourSK, orb, denseDesc%iAtomStart, img2CentCell)
      else
        call makeDensityMatrix(work, eigvecs(:,:,iKS), filling(:,iK,iSpin))
      end if
      call env%globalTimer%startTimer(globalTimers%denseToSparse)
      if (tHelical) then
        call packHS(rhoPrim(:,iSpin), work, kPoint(:,iK), kWeight(iK), neighbourList%iNeighbour,&
            & nNeighbourSK, orb%mOrb, iCellVec, cellVec, denseDesc%iAtomStart, iSparseStart,&
            & img2CentCell, orb, species, coord)
      else
        call packHS(rhoPrim(:,iSpin), work, kPoint(:,iK), kWeight(iK), neighbourList%iNeighbour,&
            & nNeighbourSK, orb%mOrb, iCellVec, cellVec, denseDesc%iAtomStart, iSparseStart,&
            & img2CentCell)
      end if
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
      & tMulliken, kPoint, kWeight, filling, neighbourList, nNeighbourSK, orb, iSparseStart,&
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
    type(TNeighbourList), intent(in) :: neighbourList

    !> Number of neighbours for each of the atoms
    integer, intent(in) :: nNeighbourSK(:)

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
    integer :: nAtom
    integer :: iKS, iK
    logical :: tImHam

    nAtom = size(orb%nOrbAtom)
    tImHam = allocated(iRhoPrim)

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
            & neighbourList%iNeighbour, nNeighbourSK, orb%mOrb, iCellVec, cellVec, iSparseStart,&
            & img2CentCell, rhoPrim, iRhoPrim)
      else
        call packRhoPauliBlacs(env%blacs, denseDesc, work, kPoint(:,iK), kWeight(iK),&
            & neighbourList%iNeighbour, nNeighbourSK, orb%mOrb, iCellVec, cellVec, iSparseStart,&
            & img2CentCell, rhoPrim)
      end if
    #:else
      if (tRealHS) then
        call packHSPauli(rhoPrim, work, neighbourlist%iNeighbour, nNeighbourSK, orb%mOrb,&
            & denseDesc%iAtomStart, iSparseStart, img2CentCell)
        if (tImHam) then
          call packHSPauliImag(iRhoPrim, work, neighbourlist%iNeighbour, nNeighbourSK, orb%mOrb,&
              & denseDesc%iAtomStart, iSparseStart, img2CentCell)
        end if
      else
        call packHS(rhoPrim, work, kPoint(:,iK), kWeight(iK), neighbourList%iNeighbour,&
            & nNeighbourSK, orb%mOrb, iCellVec, cellVec, denseDesc%iAtomStart, iSparseStart,&
            & img2CentCell)
        if (tImHam) then
          call iPackHS(iRhoPrim, work, kPoint(:,iK), kWeight(iK), neighbourlist%iNeighbour,&
              & nNeighbourSK, orb%mOrb, iCellVec, cellVec, denseDesc%iAtomStart, iSparseStart,&
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

    !> Fillings (orbital, kpoint, spin)
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
      ! Every spin channel (but not the k-points) filled up individually
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
  subroutine getMullikenPopulation(rhoPrim, over, orb, neighbourList, nNeighbourSK, img2CentCell,&
      & iSparseStart, qOrb, iRhoPrim, qBlock, qiBlock)

    !> sparse density matrix
    real(dp), intent(in) :: rhoPrim(:,:)

    !> sparse overlap matrix
    real(dp), intent(in) :: over(:)

    !> Atomic orbital information
    type(TOrbitals), intent(in) :: orb

    !> Atomic neighbours
    type(TNeighbourList), intent(in) :: neighbourList

    !> Number of neighbours for each atom within overlap distance
    integer, intent(in) :: nNeighbourSK(:)

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
      call mulliken(qOrb(:,:,iSpin), over, rhoPrim(:,iSpin), orb, neighbourList%iNeighbour,&
          & nNeighbourSK, img2CentCell, iSparseStart)
    end do

    if (allocated(qBlock)) then
      qBlock(:,:,:,:) = 0.0_dp
      do iSpin = 1, size(qBlock, dim=4)
        call mulliken(qBlock(:,:,:,iSpin), over, rhoPrim(:,iSpin), orb, neighbourList%iNeighbour,&
            & nNeighbourSK, img2CentCell, iSparseStart)
      end do
    end if

    if (allocated(qiBlock)) then
      qiBlock(:,:,:,:) = 0.0_dp
      do iSpin = 1, size(qiBlock, dim=4)
        call skewMulliken(qiBlock(:,:,:,iSpin), over, iRhoPrim(:,iSpin), orb,&
            & neighbourList%iNeighbour, nNeighbourSK, img2CentCell, iSparseStart)
      end do
    end if

  end subroutine getMullikenPopulation


  !> Calculates various energy contribution that can potentially update for the same geometry
  subroutine getEnergies(sccCalc, qOrb, q0, chargePerShell, species, tExtField, tXlbomd, tDftbU,&
      & tDualSpinOrbit, rhoPrim, H0, orb, neighbourList, nNeighbourSK, img2CentCell, iSparseStart,&
      & cellVol, extPressure, TS, potential, energy, thirdOrd, rangeSep, qDepExtPot, qBlock,&
      & qiBlock, nDftbUFunc, UJ, nUJ, iUJ, niUJ, xi, iAtInCentralRegion, tFixEf, Ef, onSiteElements)

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
    logical, intent(in) :: tExtField

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
    type(TNeighbourList), intent(in) :: neighbourList

    !> Number of neighbours within cut-off for each atom
    integer, intent(in) :: nNeighbourSK(:)

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

    !> Data from rangeseparated calculations
    type(RangeSepFunc), intent(inout), allocatable ::rangeSep

    !> Proxy for querying Q-dependant external potentials
    type(TQDepExtPotProxy), intent(inout), allocatable :: qDepExtPot

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

    integer :: nSpin
    real(dp) :: nEl(2)

    nSpin = size(rhoPrim, dim=2)

    ! Tr[H0 * Rho] can be done with the same algorithm as Mulliken-analysis
    energy%atomNonSCC(:) = 0.0_dp
    call mulliken(energy%atomNonSCC, rhoPrim(:,1), H0, orb, neighbourList%iNeighbour, nNeighbourSK,&
        & img2CentCell, iSparseStart)
    energy%EnonSCC = sum(energy%atomNonSCC(iAtInCentralRegion(:)))

    energy%atomExt(:) = 0.0_dp
    if (tExtField) then
      energy%atomExt(:) = energy%atomExt&
          & + sum(qOrb(:,:,1) - q0(:,:,1), dim=1) * potential%extAtom(:,1)
    end if
    if (allocated(qDepExtPot)) then
      call qDepExtPot%addEnergy(energy%atomExt)
    end if
    energy%Eext = sum(energy%atomExt)

    if (allocated(sccCalc)) then
      if (tXlbomd) then
        call sccCalc%getEnergyPerAtomXlbomd(species, orb, qOrb, q0, energy%atomSCC)
      else
        call sccCalc%getEnergyPerAtom(energy%atomSCC)
      end if
      energy%Escc = sum(energy%atomSCC(iAtInCentralRegion(:)))
      if (nSpin > 1) then
        energy%atomSpin(:) = 0.5_dp * sum(sum(potential%intShell(:,:,2:nSpin)&
            & * chargePerShell(:,:,2:nSpin), dim=1), dim=2)
        energy%Espin = sum(energy%atomSpin(iAtInCentralRegion(:)))
      end if
    end if

    if (allocated(thirdOrd)) then
      if (tXlbomd) then
        call thirdOrd%getEnergyPerAtomXlbomd(qOrb, q0, species, orb, energy%atom3rd)
      else
        call thirdOrd%getEnergyPerAtom(energy%atom3rd)
      end if
      energy%e3rd = sum(energy%atom3rd(iAtInCentralRegion(:)))
    end if

    if (allocated(onSiteElements)) then
      call getEons(energy%atomOnSite, qBlock, qiBlock, q0, onSiteElements, species, orb)
      energy%eOnSite = sum(energy%atomOnSite)
    end if

    if (tDftbU) then
      if (allocated(qiBlock)) then
        call E_DFTBU(energy%atomDftbu, qBlock, species, orb, nDFTBUfunc, UJ, nUJ, niUJ, iUJ,&
            & qiBlock)
      else
        call E_DFTBU(energy%atomDftbu, qBlock, species, orb, nDFTBUfunc, UJ, nUJ, niUJ, iUJ)
      end if
      energy%Edftbu = sum(energy%atomDftbu(iAtInCentralRegion(:)))
    end if

    if (tDualSpinOrbit) then
      energy%atomLS(:) = 0.0_dp
      call getDualSpinOrbitEnergy(energy%atomLS, qiBlock, xi, orb, species)
      energy%ELS = sum(energy%atomLS(iAtInCentralRegion(:)))
    end if

    energy%Eelec = energy%EnonSCC + energy%ESCC + energy%Espin + energy%ELS + energy%Edftbu&
        & + energy%Eext + energy%e3rd + energy%eOnSite

    !> Add exchange conribution for range separated calculations
    if (allocated(rangeSep)) then
      energy%Efock = 0.0_dp
      call rangeSep%addLREnergy(energy%Efock)
      energy%Eelec = energy%Eelec + energy%Efock
    end if

    energy%atomElec(:) = energy%atomNonSCC + energy%atomSCC + energy%atomSpin + energy%atomDftbu&
        & + energy%atomLS + energy%atomExt + energy%atom3rd + energy%atomOnSite
    energy%atomTotal(:) = energy%atomElec + energy%atomRep + energy%atomDisp + energy%atomHalogenX
    energy%Etotal = energy%Eelec + energy%Erep + energy%eDisp + energy%eHalogenX
    energy%EMermin = energy%Etotal - sum(TS)
    ! extrapolated to 0 K
    energy%Ezero = energy%Etotal - 0.5_dp * sum(TS)
    energy%EGibbs = energy%EMermin + cellVol * extPressure

    energy%EForceRelated = energy%EGibbs
    if (tFixEf) then
      if (nSpin == 2) then
        nEl(:) = sum(sum(qOrb(:,iAtInCentralRegion(:),:),dim=1),dim=1)
        nEl(1) = 0.5_dp * ( nEl(1) + nEl(2) )
        nEl(2) = nEl(1) - nEl(2)
        ! negative sign due to electron charge
        energy%EForceRelated = energy%EForceRelated  - sum(nEl(:2) * Ef(:2))
      else
        nEl = 0.0_dp
        nEl(1) = sum(qOrb(:,iAtInCentralRegion(:),1))
        ! negative sign due to electron charge
        energy%EForceRelated = energy%EForceRelated  - nEl(1) * Ef(1)
      end if
    end if

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
      & iGeoStep, iSccIter, minSccIter, maxSccIter, sccTol, tStopScc, tMixBlockCharges, tReadChrg,&
      & qInput, qInpRed, sccErrorQ, tConverged, qBlockOut, iEqBlockDftbU, qBlockIn, qiBlockOut,&
      & iEqBlockDftbuLS, species0, nUJ, iUJ, niUJ, qiBlockIn, iEqBlockOnSite, iEqBlockOnSiteLS)

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
    logical, intent(in) :: tMixBlockCharges

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

    !> Equivalences for onsite block corrections if needed
    integer, intent(in), allocatable :: iEqBlockOnSite(:,:,:,:)

    !> Equivalences for onsite block corrections if needed for imaginary elements
    integer, intent(in), allocatable :: iEqBlockOnSiteLS(:,:,:,:)

    real(dp), allocatable :: qDiffRed(:)
    integer :: nSpin

    nSpin = size(qOutput, dim=3)

    call reduceCharges(orb, nIneqOrb, iEqOrbitals, qOutput, qOutRed, qBlockOut, iEqBlockDftbu,&
        & qiBlockOut, iEqBlockDftbuLS, iEqBlockOnSite, iEqBlockOnSiteLS)
    qDiffRed = qOutRed - qInpRed
    sccErrorQ = maxval(abs(qDiffRed))
    tConverged = (sccErrorQ < sccTol)&
        & .and. (iSccIter >= minSccIter .or. tReadChrg .or. iGeoStep > 0)
    if ((.not. tConverged) .and. (iSccIter /= maxSccIter .and. .not. tStopScc)) then
      ! Avoid mixing of spin unpolarised density for spin polarised cases, this is only a problem in
      ! iteration 1, as there is only the (spin unpolarised!) atomic input density at that
      ! point. (Unless charges had been initialized externally)
      if ((iSCCIter + iGeoStep) == 1 .and. (nSpin > 1 .or. tMixBlockCharges) .and. .not. tReadChrg)&
          & then
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
            & species0, nUJ, iUJ, niUJ, qiBlockIn, iEqBlockDftbuLS, iEqBlockOnSite,&
            & iEqBlockOnSiteLS)
      end if
    end if

  end subroutine getNextInputCharges


  !> Update delta density matrix rather than merely q for rangeseparation
  subroutine getNextInputDensity(SSqrReal, over, neighbourList, nNeighbourSK, iAtomStart,&
      & iSparseStart, img2CentCell, pChrgMixer, qOutput, orb, tHelical, species, coord, iGeoStep,&
      & iSccIter, minSccIter, maxSccIter, sccTol, tStopScc, tReadChrg, q0, qInput, sccErrorQ,&
      & tConverged, deltaRhoOut, deltaRhoIn, deltaRhoDiff, qBlockIn, qBlockOut)

    !> Square dense overlap storage
    real(dp), allocatable, intent(inout) :: SSqrReal(:,:)

    !> sparse overlap matrix
    real(dp), intent(in) :: over(:)

    !> list of neighbours for each atom
    type(TNeighbourList), intent(in) :: neighbourList

    !> Number of neighbours for each of the atoms
    integer, intent(in) :: nNeighbourSK(:)

    !>Start of atomic blocks in dense arrays
    integer, allocatable, intent(in) :: iAtomStart(:)

    !> Index array for the start of atomic blocks in sparse arrays
    integer, intent(in) :: iSparseStart(:,:)

    !> map from image atoms to the original unique atom
    integer, intent(in) :: img2CentCell(:)

    !> Charge mixing object
    type(OMixer), intent(inout) :: pChrgMixer

    !> Output electrons
    real(dp), intent(inout) :: qOutput(:,:,:)

    !> Atomic orbital data
    type(TOrbitals), intent(in) :: orb

    !> Is the geometry helical
    logical, intent(in) :: tHelical

    !> species for atoms
    integer, intent(in) :: species(:)

    !> all coordinates
    real(dp), intent(in) :: coord(:,:)

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

    !> Were intial charges read from disc?
    logical, intent(in) :: tReadChrg

    !> reference charges
    real(dp), intent(in) :: q0(:,:,:)

    !> Resulting input charges for next SCC iteration
    real(dp), intent(inout) :: qInput(:,:,:)

    !> SCC error
    real(dp), intent(out) :: sccErrorQ

    !> Has the calculation converged>
    logical, intent(out) :: tConverged

    !> delta density matrix for rangeseparated calculations
    real(dp), intent(inout) :: deltaRhoOut(:)

    !> delta density matrix as inpurt for next SCC cycle
    real(dp), target, intent(inout) :: deltaRhoIn(:)

    !> difference of delta density matrix in and out
    real(dp), intent(inout) :: deltaRhoDiff(:)

    !> block charge input (if needed for orbital potentials)
    real(dp), intent(inout), allocatable :: qBlockIn(:,:,:,:)

    !> Dual output charges
    real(dp), intent(inout), allocatable :: qBlockOut(:,:,:,:)


    integer :: nSpin, iSpin, iAt, iOrb
    real(dp), pointer :: deltaRhoInSqr(:,:,:)

    nSpin = size(qOutput, dim=3)

    deltaRhoDiff(:) = deltaRhoOut - deltaRhoIn
    sccErrorQ = maxval(abs(deltaRhoDiff))
    tConverged = (sccErrorQ < sccTol)&
        & .and. (iSCCiter >= minSCCIter .or. tReadChrg .or. iGeoStep > 0)

    if ((.not. tConverged) .and. (iSCCiter /= maxSccIter .and. .not. tStopScc)) then
      if ((iSCCIter + iGeoStep) == 1 .and. (nSpin > 1 .and. .not. tReadChrg)) then
        deltaRhoIn(:) = deltaRhoOut
        qInput(:,:,:) = qOutput
        if (allocated(qBlockIn)) then
          qBlockIn(:,:,:,:) = qBlockOut
        end if
      else
        call mix(pChrgMixer, deltaRhoIn, deltaRhoDiff)
        if (tHelical) then
          call unpackHS(SSqrReal, over, neighbourList%iNeighbour, nNeighbourSK, iAtomStart,&
              & iSparseStart, img2CentCell, orb, species, coord)
        else
          call unpackHS(SSqrReal, over, neighbourList%iNeighbour, nNeighbourSK, iAtomStart,&
              & iSparseStart, img2CentCell)
        end if
        deltaRhoInSqr(1:orb%nOrb, 1:orb%nOrb, 1:nSpin) => deltaRhoIn
        call denseMulliken(deltaRhoInSqr, SSqrReal, iAtomStart, qInput)

        ! RangeSep: for spin-unrestricted calculation the initial guess should be equally
        ! distributed to alpha and beta density matrices
        if(nSpin == 2) then
          qInput(:,:,1) = qInput(:,:,1) + q0(:,:,1) * 0.5_dp
          qInput(:,:,2) = qInput(:,:,2) + q0(:,:,1) * 0.5_dp
        else
          qInput(:,:,:) = qInput + q0
        end if

        if (allocated(qBlockIn)) then
          call denseBlockMulliken(deltaRhoInSqr, SSqrReal, iAtomStart, qBlockIn)
          do iSpin = 1, nSpin
            do iAt = 1, size(qInput, dim=2)
              do iOrb = 1, size(qInput, dim=1)
                qBlockIn(iOrb, iOrb, iAt, iSpin) = qInput(iOrb, iAt, iSpin)
              end do
            end do
          end do
        end if

        call ud2qm(qInput)
        if (allocated(qBlockIn)) then
          call ud2qm(qBlockIn)
        end if
      end if
    end if

  end subroutine getNextInputDensity


  !> Reduce charges according to orbital equivalency rules.
  subroutine reduceCharges(orb, nIneqOrb, iEqOrbitals, qOrb, qRed, qBlock, iEqBlockDftbu, qiBlock,&
      & iEqBlockDftbuLS, iEqBlockOnSite, iEqBlockOnSiteLS)

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

    !> Equivalences for onsite block corrections if needed
    integer, intent(in), allocatable :: iEqBlockOnSite(:,:,:,:)

    !> Equivalences for onsite block corrections if needed for imaginary part
    integer, intent(in), allocatable :: iEqBlockOnSiteLS(:,:,:,:)

    real(dp), allocatable :: qOrbUpDown(:,:,:), qBlockUpDown(:,:,:,:)

    qRed(:) = 0.0_dp
    qOrbUpDown = qOrb
    call qm2ud(qOrbUpDown)
    call orbitalEquiv_reduce(qOrbUpDown, iEqOrbitals, orb, qRed(1:nIneqOrb))
    if (allocated(qBlock)) then
      qBlockUpDown = qBlock
      call qm2ud(qBlockUpDown)
      if (allocated(iEqBlockOnSite)) then
        ! all blocks are full of unique elements
        call onsBlock_reduce(qBlockUpDown, iEqBlockOnSite, orb, qRed)
        if (allocated(qiBlock)) then
          call onsBlock_reduce(qiBlock, iEqBlockOnSiteLS, orb, qRed, skew=.true.)
        end if
      else
        ! only a subset of blocks are covered in +U type operations
        call appendBlock_reduce(qBlockUpDown, iEqBlockDFTBU, orb, qRed)
        if (allocated(qiBlock)) then
          call appendBlock_reduce(qiBlock, iEqBlockDFTBULS, orb, qRed, skew=.true.)
        end if
      end if
    end if

  end subroutine reduceCharges


  !> Expand reduced charges according orbital equivalency rules.
  subroutine expandCharges(qRed, orb, nIneqOrb, iEqOrbitals, qOrb, qBlock, iEqBlockDftbu, species0,&
      & nUJ, iUJ, niUJ, qiBlock, iEqBlockDftbuLS, iEqBlockOnSite, iEqBlockOnSiteLS)

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

    !> Equivalences for onsite block corrections if needed
    integer, intent(in), allocatable :: iEqBlockOnSite(:,:,:,:)

    !> Equivalences for onsite block corrections if needed for imaginary part
    integer, intent(in), allocatable :: iEqBlockOnSiteLS(:,:,:,:)

    integer :: nSpin

    @:ASSERT(allocated(qBlock) .eqv. (allocated(iEqBlockDftbU) .or. allocated(iEqBlockOnSite)))
    @:ASSERT(.not. allocated(qBlock) .or. allocated(species0))
    @:ASSERT(.not. allocated(qBlock) .or. allocated(nUJ))
    @:ASSERT(.not. allocated(qBlock) .or. allocated(iUJ))
    @:ASSERT(.not. allocated(qBlock) .or. allocated(niUJ))
    @:ASSERT(.not. allocated(qiBlock) .or. allocated(qBlock))
    @:ASSERT(allocated(qiBlock) .eqv. (allocated(iEqBlockDftbuLS) .or. allocated(iEqBlockOnSiteLS)))

    nSpin = size(qOrb, dim=3)
    call OrbitalEquiv_expand(qRed(1:nIneqOrb), iEqOrbitals, orb, qOrb)
    if (allocated(qBlock)) then
      qBlock(:,:,:,:) = 0.0_dp
      if (allocated(iEqBlockOnSite)) then
        ! all blocks are full of unique elements
        call Onsblock_expand(qRed, iEqBlockOnSite, orb, qBlock, orbEquiv=iEqOrbitals)
        if (allocated(qiBlock)) then
          call Onsblock_expand(qRed, iEqBlockOnSiteLS, orb, qiBlock, skew=.true.)
        end if
      else
        ! only a subset of blocks are covered in +U type operations
        call Block_expand(qRed, iEqBlockDftbu, orb, qBlock, species0, nUJ, niUJ, iUJ,&
            & orbEquiv=iEqOrbitals)
        if (allocated(qiBlock)) then
          call Block_expand(qRed, iEqBlockDftbuLS, orb, qiBlock, species0, nUJ, niUJ, iUJ,&
              & skew=.true.)
        end if
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

    if (iSccIter > 1) then
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

    !> forces being evaluated in the excited state
    logical, intent(in) :: tForces

    if (t3rd) then
      call error("Third order currently incompatible with excited state")
    end if
    if (.not. tRealHS) then
      call error("Only real systems are supported for excited state calculations")
    end if
    if (tPeriodic .and. tForces) then
      call error("Forces in the excited state for periodic geometries are currently unavailable")
    end if

  end subroutine ensureLinRespConditions


  !> Do the linear response excitation calculation.
  subroutine calculateLinRespExcitations(env, lresp, parallelKS, sccCalc, qOutput, q0, over,&
      & eigvecsReal, eigen, filling, coord, species, speciesName, orb, tHelical, skHamCont,&
      & skOverCont, autotestTag, taggedWriter, runId, neighbourList, nNeighbourSk, denseDesc,&
      & iSparseStart, img2CentCell, tWriteAutotest, tForces, tLinRespZVect, tPrintExcEigvecs,&
      & tPrintExcEigvecsTxt, nonSccDeriv, energy, energies, work, rhoSqrReal, excitedDerivs,&
      & occNatural)

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

    !> ground state eigenvalues (orbital, kpoint)
    real(dp), intent(in) :: eigen(:,:)

    !> ground state fillings (orbital, kpoint)
    real(dp), intent(in) :: filling(:,:)

    !> all atomic coordinates
    real(dp), intent(in) :: coord(:,:)

    !> species of all atoms in the system
    integer, target, intent(in) :: species(:)

    !> label for each atomic chemical species
    character(*), intent(in) :: speciesName(:)

    !> Atomic orbital information
    type(TOrbitals), intent(in) :: orb

    !> Is the geometry helical
    logical, intent(in) :: tHelical

    !> non-SCC hamiltonian information
    type(OSlakoCont), intent(in) :: skHamCont

    !> overlap information
    type(OSlakoCont), intent(in) :: skOverCont

    !> File name for regression data
    character(*), intent(in) :: autotestTag

    !> Tagged writer
    type(TTaggedWriter), intent(inout) :: taggedWriter

    !> Job ID for future identification
    integer, intent(in) :: runId

    !> list of neighbours for each atom
    type(TNeighbourList), intent(in) :: neighbourList

    !> Number of neighbours for each of the atoms
    integer, intent(in) :: nNeighbourSK(:)

    !> Dense matrix descriptor
    type(TDenseDescr), intent(in) :: denseDesc

    !> Index array for the start of atomic blocks in sparse arrays
    integer, intent(in) :: iSparseStart(:,:)

    !> map from image atoms to the original unique atom
    integer, intent(in) :: img2CentCell(:)

    !> should regression test data be written
    logical, intent(in) :: tWriteAutotest

    !> forces to be calculated in the excited state
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

    !> energes of all solved states
    real(dp), intent(inout), allocatable :: energies(:)

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
    integer :: iSpin, nSpin, nAtom, fdAutotest
    logical :: tSpin

    nAtom = size(qOutput, dim=2)
    nSpin = size(eigen, dim=2)
    tSpin = (nSpin == 2)
    pSpecies0 => species(1:nAtom)

    energy%Eexcited = 0.0_dp
    allocate(dQAtom(nAtom))
    dQAtom(:) = sum(qOutput(:,:,1) - q0(:,:,1), dim=1)
    if (tHelical) then
      call unpackHS(work, over, neighbourList%iNeighbour, nNeighbourSK, denseDesc%iAtomStart,&
          & iSparseStart, img2CentCell, orb, species, coord)
    else
      call unpackHS(work, over, neighbourList%iNeighbour, nNeighbourSK, denseDesc%iAtomStart,&
          & iSparseStart, img2CentCell)
    end if
    call blockSymmetrizeHS(work, denseDesc%iAtomStart)
    if (tForces) then
      do iSpin = 1, nSpin
        call blockSymmetrizeHS(rhoSqrReal(:,:,iSpin), denseDesc%iAtomStart)
      end do
    end if
    if (tWriteAutotest) then
      open(newUnit=fdAutotest, file=autotestTag, position="append")
    end if

    if (tLinRespZVect) then
      if (tPrintExcEigVecs) then
        allocate(naturalOrbs(orb%nOrb, orb%nOrb, 1))
      end if
      call addGradients(tSpin, lresp, denseDesc%iAtomStart, eigvecsReal, eigen, work, filling,&
          & coord(:,:nAtom), sccCalc, dQAtom, pSpecies0, neighbourList%iNeighbour, img2CentCell,&
          & orb, skHamCont, skOverCont, tWriteAutotest, fdAutotest, taggedWriter, energy%Eexcited,&
         & energies, excitedDerivs, nonSccDeriv, rhoSqrReal, occNatural, naturalOrbs)
      if (tPrintExcEigvecs) then
        call writeRealEigvecs(env, runId, neighbourList, nNeighbourSK, denseDesc, iSparseStart,&
            & img2CentCell, pSpecies0, speciesName, orb, over, parallelKS, tPrintExcEigvecsTxt,&
            & naturalOrbs, work, fileName="excitedOrbs")
      end if
    else
      call calcExcitations(lresp, tSpin, denseDesc, eigvecsReal, eigen, work, filling,&
          & coord(:,:nAtom), sccCalc, dQAtom, pSpecies0, neighbourList%iNeighbour, img2CentCell,&
          & orb, tWriteAutotest, fdAutotest, taggedWriter, energy%Eexcited, energies)
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
      & qiBlockIn, iEqBlockOnSite, iEqBlockOnSiteLS)

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

    !> Equivalences for onsite block corrections if needed
    integer, intent(inout), allocatable :: iEqBlockOnSite(:,:,:,:)

    !> Equivalences for onsite block corrections if needed for imaginary part
    integer, intent(inout), allocatable :: iEqBlockOnSiteLS(:,:,:,:)

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
        & species0, nUJ, iUJ, niUJ, qiBlockIn, iEqBlockDftbuLS, iEqBlockOnSite, iEqBlockOnSiteLS)

  end subroutine getXlbomdCharges


  !> Calculates dipole moment.
  subroutine getDipoleMoment(qOutput, q0, coord, dipoleMoment, iAtInCentralRegion)

    !> electrons in orbitals
    real(dp), intent(in) :: qOutput(:,:,:)

    !> reference atomic charges
    real(dp), intent(in) :: q0(:,:,:)

    !> atomic coordinates
    real(dp), intent(in) :: coord(:,:)

    !> resulting dipole moment
    real(dp), intent(out) :: dipoleMoment(:)

    !> atoms in the central cell (or device region if transport)
    integer, intent(in) :: iAtInCentralRegion(:)

    integer :: nAtom, ii, iAtom

    nAtom = size(qOutput, dim=2)
    dipoleMoment(:) = 0.0_dp
    do ii = 1, size(iAtInCentralRegion)
      iAtom = iAtInCentralRegion(ii)
      dipoleMoment(:) = dipoleMoment(:)&
          & + sum(q0(:, iAtom, 1) - qOutput(:, iAtom, 1)) * coord(:,iAtom)
    end do

  end subroutine getDipoleMoment


  !> Prints dipole moment calcululated by the derivative of H with respect of the external field.
  subroutine checkDipoleViaHellmannFeynman(rhoPrim, q0, coord0, over, orb, neighbourList,&
      & nNeighbourSK, species, iSparseStart, img2CentCell)

    !> Density matrix in sparse storage
    real(dp), intent(in) :: rhoPrim(:,:)

    !> Reference orbital charges
    real(dp), intent(in) :: q0(:,:,:)

    !> Central cell atomic coordinates
    real(dp), intent(in) :: coord0(:,:)

    !> Overlap matrix in sparse storage
    real(dp), intent(in) :: over(:)

    !> Atomic orbital information
    type(TOrbitals), intent(in) :: orb

    !> list of neighbours for each atom
    type(TNeighbourList), intent(in) :: neighbourList

    !> Number of neighbours for each of the atoms
    integer, intent(in) :: nNeighbourSK(:)

    !> species of all atoms in the system
    integer, intent(in) :: species(:)

    !> Index array for the start of atomic blocks in sparse arrays
    integer, intent(in) :: iSparseStart(:,:)

    !> map from image atoms to the original unique atom
    integer, intent(in) :: img2CentCell(:)

    real(dp), allocatable :: hprime(:,:), dipole(:,:), potentialDerivative(:,:)
    integer :: nAtom, sparseSize, iAt, ii

    sparseSize = size(over)
    nAtom = size(q0, dim=2)
    allocate(hprime(sparseSize, 1))
    allocate(dipole(size(q0, dim=1), nAtom))
    allocate(potentialDerivative(nAtom, 1))
    write(stdOut,*)
    write(stdOut, "(A)", advance='no') 'Hellmann Feynman dipole:'

    ! loop over directions
    do ii = 1, 3
      potentialDerivative(:,:) = 0.0_dp
      ! Potential from dH/dE
      potentialDerivative(:,1) = -coord0(ii,:)
      hprime(:,:) = 0.0_dp
      dipole(:,:) = 0.0_dp
      call add_shift(hprime, over, nNeighbourSK, neighbourList%iNeighbour, species, orb,&
          & iSparseStart, nAtom, img2CentCell, potentialDerivative)

      ! evaluate <psi| dH/dE | psi>
      call mulliken(dipole, hprime(:,1), rhoPrim(:,1), orb, neighbourList%iNeighbour, nNeighbourSK,&
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
  subroutine getEnergyWeightedDensity(env, electronicSolver, denseDesc, forceType, filling,&
      & eigen, kPoint, kWeight, neighbourList, nNeighbourSK, orb, iSparseStart, img2CentCell,&
      & iCellVEc, cellVec, tRealHS, ham, over, parallelKS, tHelical, species, coord, iSCC, mu,&
      & ERhoPrim, HSqrReal, SSqrReal, HSqrCplx, SSqrCplx)

    !> Environment settings
    type(TEnvironment), intent(inout) :: env

    !> Electronic solver information
    type(TElectronicSolver), intent(inout) :: electronicSolver

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
    type(TNeighbourList), intent(in) :: neighbourList

    !> Number of neighbours for each of the atoms
    integer, intent(in) :: nNeighbourSK(:)

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

    !> Is the geometry helical
    logical, intent(in) :: tHelical

    !> species of atoms
    integer, intent(in), optional :: species(:)

    !> all coordinates
    real(dp), intent(in) :: coord(:,:)

    !> iteration counter
    integer, intent(in) :: iSCC

    !> Electrochemical potentials per contact and spin
    real(dp), allocatable, intent(in) :: mu(:,:)

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

    integer :: nSpin

    nSpin = size(ham, dim=2)

    call env%globalTimer%startTimer(globalTimers%energyDensityMatrix)

    select case (electronicSolver%iSolver)

    case (electronicSolverTypes%GF)

      if (forceType /= forceTypes%orig) then
        call error("Alternative force evaluation methods are not supported by this electronic&
            & solver.")
      end if

    #:if WITH_TRANSPORT
      if (electronicSolver%iSolver == electronicSolverTypes%GF) then
        call calcEdensity_green(iSCC, env, parallelKS%localKS, ham, over, neighbourlist%iNeighbour,&
            & nNeighbourSK, denseDesc%iAtomStart, iSparseStart, img2CentCell, iCellVec, cellVec,&
            & orb, kPoint, kWeight, mu, ERhoPrim)
      end if
    #:endif

    case (electronicSolverTypes%onlyTransport)

      call error("The OnlyTransport solver cannot calculate the energy weighted density matrix")

    case (electronicSolverTypes%qr, electronicSolverTypes%divideandconquer,&
        & electronicSolverTypes%relativelyrobust, electronicSolverTypes%elpa, &
        & electronicSolverTypes%magma_gvd)

      call getEDensityMtxFromEigvecs(env, denseDesc, forceType, filling, eigen, kPoint, kWeight,&
          & neighbourList, nNeighbourSK, orb, iSparseStart, img2CentCell, iCellVec, cellVec,&
          & tRealHS, ham, over, parallelKS, tHelical, species, coord, ERhoPrim, HSqrReal, SSqrReal,&
          & HSqrCplx, SSqrCplx)

    case (electronicSolverTypes%omm, electronicSolverTypes%pexsi, electronicSolverTypes%ntpoly)

      if (forceType /= forceTypes%orig) then
        call error("Alternative force evaluation methods are not supported by this electronic&
            & solver.")
      end if

      call electronicSolver%elsi%getEDensity(env, denseDesc, nSpin, kPoint, kWeight, neighbourList,&
          & nNeighbourSK, tHelical, orb, species, coord, iSparseStart, img2CentCell, iCellVec,&
          & cellVec, tRealHS, parallelKS, ERhoPrim, SSqrReal, SSqrCplx)

    end select

    call env%globalTimer%stopTimer(globalTimers%energyDensityMatrix)

  end subroutine getEnergyWeightedDensity


  !> Calculates the energy weighted density matrix using eigenvectors
  subroutine getEDensityMtxFromEigvecs(env, denseDesc, forceType, filling, eigen, kPoint, kWeight,&
      & neighbourList, nNeighbourSK, orb, iSparseStart, img2CentCell, iCellVec, cellVec, tRealHS,&
      & ham, over, parallelKS, tHelical, species, coord, ERhoPrim, HSqrReal, SSqrReal, HSqrCplx,&
      & SSqrCplx)

    !> Environment settings
    type(TEnvironment), intent(inout) :: env

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
    type(TNeighbourList), intent(in) :: neighbourList

    !> Number of neighbours for each of the atoms
    integer, intent(in) :: nNeighbourSK(:)

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

    !> Is the geometry helical
    logical, intent(in) :: tHelical

    !> species of atoms
    integer, intent(in), optional :: species(:)

    !> all coordinates
    real(dp), intent(in) :: coord(:,:)

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

    integer :: nSpin

    nSpin = size(ham, dim=2)

    if (nSpin == 4) then
      call getEDensityMtxFromPauliEigvecs(env, denseDesc, filling, eigen, kPoint, kWeight,&
          & neighbourList, nNeighbourSK, orb, iSparseStart, img2CentCell, iCellVec, cellVec,&
          & tRealHS, parallelKS, HSqrCplx, SSqrCplx, ERhoPrim)
    else
      if (tRealHS) then
        call getEDensityMtxFromRealEigvecs(env, denseDesc, forceType, filling, eigen,&
            & neighbourList, nNeighbourSK, orb, iSparseStart, img2CentCell, ham, over,&
            & parallelKS, tHelical, species, coord, HSqrReal, SSqrReal, ERhoPrim)
      else
        call getEDensityMtxFromComplexEigvecs(env, denseDesc, forceType, filling, eigen, kPoint,&
            & kWeight, neighbourList, nNeighbourSK, orb, iSparseStart, img2CentCell, iCellVec,&
            & cellVec, ham, over, parallelKS, tHelical, species, coord, HSqrCplx, SSqrCplx,&
            & ERhoPrim)
      end if
    end if

  end subroutine getEDensityMtxFromEigvecs


  !> Calculates density matrix from real eigenvectors.
  subroutine getEDensityMtxFromRealEigvecs(env, denseDesc, forceType, filling, eigen,&
      & neighbourList, nNeighbourSK, orb, iSparseStart, img2CentCell, ham, over, parallelKS,&
      & tHelical, species, coord, eigvecsReal, work, ERhoPrim)

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
    type(TNeighbourList), intent(in) :: neighbourList

    !> Number of neighbours for each of the atoms
    integer, intent(in) :: nNeighbourSK(:)

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

    !> Is the geometry haelical
    logical, intent(in) :: tHelical

    !> species of atoms
    integer, intent(in), optional :: species(:)

    !> all coordinates
    real(dp), intent(in) :: coord(:,:)

    !> Eigenvectors (NOTE: they will be rewritten with work data on exit!)
    real(dp), intent(inout) :: eigvecsReal(:,:,:)

    !> Work array for storing temporary data
    real(dp), intent(out) :: work(:,:)

    !> Energy weighted density matrix
    real(dp), intent(out) :: ERhoPrim(:)

    real(dp), allocatable :: work2(:,:)
    integer :: nLocalRows, nLocalCols
    integer :: iKS, iS

    nLocalRows = size(eigvecsReal, dim=1)
    nLocalCols = size(eigvecsReal, dim=2)
    if (forceType == forceTypes%dynamicT0 .or. forceType == forceTypes%dynamicTFinite) then
      allocate(work2(nLocalRows, nLocalCols))
    end if

    ERhoPrim(:) = 0.0_dp
    do iKS = 1, parallelKS%nLocalKS
      iS = parallelKS%localKS(2, iKS)

      select case (forceType)

      case(forceTypes%orig)
        ! Original (non-consistent) scheme
      #:if WITH_SCALAPACK
        call makeDensityMtxRealBlacs(env%blacs%orbitalGrid, denseDesc%blacsOrbSqr, filling(:,1,iS),&
            & eigvecsReal(:,:,iKS), work, eigen(:,1,iS))
      #:else
        if (tDensON2) then
          call makeDensityMatrix(work, eigvecsReal(:,:,iKS), filling(:,1,iS), eigen(:,1,iS),&
              & neighbourlist%iNeighbour, nNeighbourSK, orb, denseDesc%iAtomStart, img2CentCell)
        else
          call makeDensityMatrix(work, eigvecsReal(:,:,iKS), filling(:,1,iS), eigen(:,1,iS))
        end if
      #:endif

      case(forceTypes%dynamicT0)
        ! Correct force for XLBOMD for T=0K (DHD)
      #:if WITH_SCALAPACK
        if (tHelical) then
          call unpackHSHelicalRealBlacs(env%blacs, ham(:,iS), neighbourList%iNeighbour,&
              & nNeighbourSK, iSparseStart, img2CentCell, orb, species, coord, denseDesc, work)
        else
          call unpackHSRealBlacs(env%blacs, ham(:,iS), neighbourList%iNeighbour, nNeighbourSK,&
              & iSparseStart, img2CentCell, denseDesc, work)
        end if
        call makeDensityMtxRealBlacs(env%blacs%orbitalGrid, denseDesc%blacsOrbSqr, filling(:,1,iS),&
            & eigVecsReal(:,:,iKS), work2)
        call pblasfx_psymm(work2, denseDesc%blacsOrbSqr, work, denseDesc%blacsOrbSqr,&
            & eigvecsReal(:,:,iKS), denseDesc%blacsOrbSqr, side="L")
        call pblasfx_psymm(work2, denseDesc%blacsOrbSqr, eigvecsReal(:,:,iKS),&
            & denseDesc%blacsOrbSqr, work, denseDesc%blacsOrbSqr, side="R", alpha=0.5_dp)
      #:else
        if (tHelical) then
          call unpackHS(work, ham(:,iS), neighbourlist%iNeighbour, nNeighbourSK,&
              & denseDesc%iAtomStart, iSparseStart, img2CentCell, orb, species, coord)
        else
          call unpackHS(work, ham(:,iS), neighbourlist%iNeighbour, nNeighbourSK,&
              & denseDesc%iAtomStart, iSparseStart, img2CentCell)
        end if
        call blockSymmetrizeHS(work, denseDesc%iAtomStart)
        call makeDensityMatrix(work2, eigvecsReal(:,:,iKS), filling(:,1,iS))
        ! D H
        call symm(eigvecsReal(:,:,iKS), "L", work2, work)
        ! (D H) D
        call symm(work, "R", work2, eigvecsReal(:,:,iKS), alpha=0.5_dp)
      #:endif

      case(forceTypes%dynamicTFinite)
        ! Correct force for XLBOMD for T <> 0K (DHS^-1 + S^-1HD)
      #:if WITH_SCALAPACK
        call makeDensityMtxRealBlacs(env%blacs%orbitalGrid, denseDesc%blacsOrbSqr, filling(:,1,iS),&
            & eigVecsReal(:,:,iKS), work)
        if (tHelical) then
          call unpackHSHelicalRealBlacs(env%blacs, ham(:,iS), neighbourlist%iNeighbour,&
              & nNeighbourSK, iSparseStart, img2CentCell, orb, species, coord, denseDesc, work2)
        else
          call unpackHSRealBlacs(env%blacs, ham(:,iS), neighbourlist%iNeighbour, nNeighbourSK,&
              & iSparseStart, img2CentCell, denseDesc, work2)
        end if
        call pblasfx_psymm(work, denseDesc%blacsOrbSqr, work2, denseDesc%blacsOrbSqr,&
            & eigvecsReal(:,:,iKS), denseDesc%blacsOrbSqr, side="L")
        if (tHelical) then
          call unpackHSHelicalRealBlacs(env%blacs, over, neighbourlist%iNeighbour, nNeighbourSK,&
              & iSparseStart, img2CentCell, orb, species, coord, denseDesc, work)
        else
          call unpackHSRealBlacs(env%blacs, over, neighbourlist%iNeighbour, nNeighbourSK,&
              & iSparseStart, img2CentCell, denseDesc, work)
        end if
        call psymmatinv(denseDesc%blacsOrbSqr, work)
        call pblasfx_psymm(work, denseDesc%blacsOrbSqr, eigvecsReal(:,:,iKS),&
            & denseDesc%blacsOrbSqr, work2, denseDesc%blacsOrbSqr, side="R", alpha=0.5_dp)
        work(:,:) = work2
        call pblasfx_ptran(work2, denseDesc%blacsOrbSqr, work, denseDesc%blacsOrbSqr, alpha=1.0_dp,&
            & beta=1.0_dp)
      #:else
        call makeDensityMatrix(work, eigvecsReal(:,:,iKS), filling(:,1,iS))
        if (tHelical) then
          call unpackHS(work2, ham(:,iS), neighbourlist%iNeighbour, nNeighbourSK,&
              & denseDesc%iAtomStart, iSparseStart, img2CentCell, orb, species, coord)
        else
          call unpackHS(work2, ham(:,iS), neighbourlist%iNeighbour, nNeighbourSK,&
              & denseDesc%iAtomStart, iSparseStart, img2CentCell)
        end if
        call blocksymmetrizeHS(work2, denseDesc%iAtomStart)
        call symm(eigvecsReal(:,:,iKS), "L", work, work2)
        if (tHelical) then
          call unpackHS(work, over, neighbourlist%iNeighbour, nNeighbourSK, denseDesc%iAtomStart,&
              & iSparseStart, img2CentCell, orb, species, coord)
        else
          call unpackHS(work, over, neighbourlist%iNeighbour, nNeighbourSK, denseDesc%iAtomStart,&
              & iSparseStart, img2CentCell)
        end if
        call symmatinv(work)
        call symm(work2, "R", work, eigvecsReal(:,:,iKS), alpha=0.5_dp)
        work(:,:) = work2 + transpose(work2)
      #:endif
      end select

    #:if WITH_SCALAPACK
      if (tHelical) then
        call packRhoHelicalRealBlacs(env%blacs, denseDesc, work, neighbourList%iNeighbour,&
            & nNeighbourSK, iSparseStart, img2CentCell, orb, species, coord, ERhoPrim)
      else
        call packRhoRealBlacs(env%blacs, denseDesc, work, neighbourList%iNeighbour, nNeighbourSK,&
            & orb%mOrb, iSparseStart, img2CentCell, ERhoPrim)
      end if
    #:else
      if (tHelical) then
        call packHS(ERhoPrim, work, neighbourList%iNeighbour, nNeighbourSK,&
            & denseDesc%iAtomStart, iSparseStart, img2CentCell, orb, species, coord)
      else
        call packHS(ERhoPrim, work, neighbourList%iNeighbour, nNeighbourSK, orb%mOrb,&
            & denseDesc%iAtomStart, iSparseStart, img2CentCell)
      end if
    #:endif
    end do

  #:if WITH_SCALAPACK
    ! Add up and distribute energy weighted density matrix contribution from each group
    call mpifx_allreduceip(env%mpi%globalComm, ERhoPrim, MPI_SUM)
  #:endif

  end subroutine getEDensityMtxFromRealEigvecs


  !> Calculates density matrix from complex eigenvectors.
  subroutine getEDensityMtxFromComplexEigvecs(env, denseDesc, forceType, filling, eigen, kPoint,&
      & kWeight, neighbourList, nNeighbourSK, orb, iSparseStart, img2CentCell, iCellVec, cellVec,&
      & ham, over, parallelKS, tHelical, species, coord, eigvecsCplx, work, ERhoPrim)

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
    type(TNeighbourList), intent(in) :: neighbourList

    !> Number of neighbours for each of the atoms
    integer, intent(in) :: nNeighbourSK(:)

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

        !> Is the geometry haelical
    logical, intent(in) :: tHelical

    !> species of atoms
    integer, intent(in), optional :: species(:)

    !> all coordinates
    real(dp), intent(in) :: coord(:,:)

    !> Eigenvectors of the system
    complex(dp), intent(inout) :: eigvecsCplx(:,:,:)

    !> work array (sized like overlap matrix)
    complex(dp), intent(inout) :: work(:,:)

    !> Energy weighted sparse density matrix (charge only part)
    real(dp), intent(out) :: ERhoPrim(:)

    complex(dp), allocatable :: work2(:,:)
    integer :: nLocalRows, nLocalCols
    integer :: iKS, iS, iK

    nLocalRows = size(eigvecsCplx, dim=1)
    nLocalCols = size(eigvecsCplx, dim=2)

    if (forceType == forceTypes%dynamicT0 .or. forceType == forceTypes%dynamicTFinite) then
      allocate(work2(nLocalRows, nLocalCols))
    end if

    ERhoPrim(:) = 0.0_dp
    do iKS = 1, parallelKS%nLocalKS
      iK = parallelKS%localKS(1, iKS)
      iS = parallelKS%localKS(2, iKS)

      select case (forceType)

      case(forceTypes%orig)
        ! Original (non-consistent) scheme
      #:if WITH_SCALAPACK
        call makeDensityMtxCplxBlacs(env%blacs%orbitalGrid, denseDesc%blacsOrbSqr,&
            & filling(:,iK,iS), eigvecsCplx(:,:,iKS), work, eigen(:,iK,iS))
      #:else
        if (tDensON2) then
          call makeDensityMatrix(work, eigvecsCplx(:,:,iKS), filling(:,iK,iS),&
              & eigen(:,iK, iS), neighbourlist%iNeighbour, nNeighbourSK, orb, denseDesc%iAtomStart,&
              & img2CentCell)
        else
          call makeDensityMatrix(work, eigvecsCplx(:,:,iKS), filling(:,iK,iS), eigen(:,iK, iS))
        end if
      #:endif

      case(forceTypes%dynamicT0)
        ! Correct force for XLBOMD for T=0K (DHD)
      #:if WITH_SCALAPACK
        if (tHelical) then
          call unpackHSHelicalCplxBlacs(env%blacs, ham(:,iS), kPoint(:,iK),&
              & neighbourList%iNeighbour, nNeighbourSK, iCellVec, cellVec, iSparseStart,&
              & img2CentCell, orb, species, coord, denseDesc, work)
        else
          call unpackHSCplxBlacs(env%blacs, ham(:,iS), kPoint(:,iK), neighbourList%iNeighbour,&
              & nNeighbourSK, iCellVec, cellVec, iSparseStart, img2CentCell, denseDesc, work)
        end if
        call makeDensityMtxCplxBlacs(env%blacs%orbitalGrid, denseDesc%blacsOrbSqr, filling(:,1,iS),&
            & eigvecsCplx(:,:,iKS), work2)
        call pblasfx_phemm(work2, denseDesc%blacsOrbSqr, work, denseDesc%blacsOrbSqr,&
            & eigvecsCplx(:,:,iKS), denseDesc%blacsOrbSqr, side="L")
        call pblasfx_phemm(work2, denseDesc%blacsOrbSqr, eigvecsCplx(:,:,iKS),&
            & denseDesc%blacsOrbSqr, work, denseDesc%blacsOrbSqr, side="R", alpha=(0.5_dp, 0.0_dp))
      #:else
        call makeDensityMatrix(work2, eigvecsCplx(:,:,iKS), filling(:,iK,iS))
        if (tHelical) then
          call unpackHS(work, ham(:,iS), kPoint(:,iK), neighbourlist%iNeighbour, nNeighbourSK,&
              & iCellVec, cellVec, denseDesc%iAtomStart, iSparseStart, img2CentCell, orb, species,&
              & coord)
        else
          call unpackHS(work, ham(:,iS), kPoint(:,iK), neighbourlist%iNeighbour, nNeighbourSK,&
              & iCellVec, cellVec, denseDesc%iAtomStart, iSparseStart, img2CentCell)
        end if
        call blockHermitianHS(work, denseDesc%iAtomStart)
        call hemm(eigvecsCplx(:,:,iKS), "L", work2, work)
        call hemm(work, "R", work2, eigvecsCplx(:,:,iKS), alpha=(0.5_dp, 0.0_dp))
      #:endif

      case(forceTypes%dynamicTFinite)
        ! Correct force for XLBOMD for T <> 0K (DHS^-1 + S^-1HD)
      #:if WITH_SCALAPACK
        call makeDensityMtxCplxBlacs(env%blacs%orbitalGrid, denseDesc%blacsOrbSqr,&
            & filling(:,iK,iS), eigVecsCplx(:,:,iKS), work)
        if (tHelical) then
          call unpackHSHelicalCplxBlacs(env%blacs, ham(:,iS), kPoint(:,iK),&
              & neighbourlist%iNeighbour, nNeighbourSK, iCellVec, cellVec, iSparseStart,&
              & img2CentCell, orb, species, coord, denseDesc, work2)
        else
          call unpackHSCplxBlacs(env%blacs, ham(:,iS), kPoint(:,iK), neighbourlist%iNeighbour,&
              & nNeighbourSK, iCellVec, cellVec, iSparseStart, img2CentCell, denseDesc, work2)
        end if
        call pblasfx_phemm(work, denseDesc%blacsOrbSqr, work2, denseDesc%blacsOrbSqr,&
            & eigvecsCplx(:,:,iKS), denseDesc%blacsOrbSqr, side="L")
        if (tHelical) then
          call unpackHSHelicalCplxBlacs(env%blacs, over, kPoint(:,iK), neighbourlist%iNeighbour,&
              & nNeighbourSK, iCellVec, cellVec, iSparseStart, img2CentCell, orb, species, coord,&
              & denseDesc, work)
        else
          call unpackHSCplxBlacs(env%blacs, over, kPoint(:,iK), neighbourlist%iNeighbour,&
              & nNeighbourSK, iCellVec, cellVec, iSparseStart, img2CentCell, denseDesc, work)
        end if
        call phermatinv(denseDesc%blacsOrbSqr, work)
        call pblasfx_phemm(work, denseDesc%blacsOrbSqr, eigvecsCplx(:,:,iKS),&
            & denseDesc%blacsOrbSqr, work2, denseDesc%blacsOrbSqr, side="R", alpha=(0.5_dp, 0.0_dp))
        work(:,:) = work2
        call pblasfx_ptranc(work2, denseDesc%blacsOrbSqr, work, denseDesc%blacsOrbSqr,&
            & alpha=(1.0_dp, 0.0_dp), beta=(1.0_dp, 0.0_dp))
      #:else
        call makeDensityMatrix(work, eigvecsCplx(:,:,iKS), filling(:,iK,iS))
        if (tHelical) then
          call unpackHS(work2, ham(:,iS), kPoint(:,iK), neighbourlist%iNeighbour, nNeighbourSK,&
              & iCellVec, cellVec, denseDesc%iAtomStart, iSparseStart, img2CentCell, orb, species,&
              & coord)
        else
          call unpackHS(work2, ham(:,iS), kPoint(:,iK), neighbourlist%iNeighbour, nNeighbourSK,&
              & iCellVec, cellVec, denseDesc%iAtomStart, iSparseStart, img2CentCell)
        end if
        call blockHermitianHS(work2, denseDesc%iAtomStart)
        call hemm(eigvecsCplx(:,:,iKS), "L", work, work2)
        if  (tHelical) then
          call unpackHS(work, over, kPoint(:,iK), neighbourlist%iNeighbour, nNeighbourSK, iCellVec,&
              & cellVec, denseDesc%iAtomStart, iSparseStart, img2CentCell, orb, species, coord)
        else
          call unpackHS(work, over, kPoint(:,iK), neighbourlist%iNeighbour, nNeighbourSK, iCellVec,&
              & cellVec, denseDesc%iAtomStart, iSparseStart, img2CentCell)
        end if
        call hermatinv(work)
        call hemm(work2, "R", work, eigvecsCplx(:,:,iKS), alpha=(0.5_dp, 0.0_dp))
        work(:,:) = work2 + transpose(conjg(work2))
      #:endif
      end select

    #:if WITH_SCALAPACK
      if (tHelical) then
        call packRhoHelicalCplxBlacs(env%blacs, denseDesc, work, kPoint(:,iK), kWeight(iK),&
            & neighbourList%iNeighbour, nNeighbourSK, iCellVec, cellVec, iSparseStart,&
            & img2CentCell, orb, species, coord, ERhoPrim)
      else
        call packRhoCplxBlacs(env%blacs, denseDesc, work, kPoint(:,iK), kWeight(iK),&
            &neighbourList%iNeighbour, nNeighbourSK, orb%mOrb, iCellVec, cellVec, iSparseStart,&
            & img2CentCell, ERhoPrim)
      end if
    #:else
      if (tHelical) then
        call packHS(ERhoPrim, work, kPoint(:,iK), kWeight(iK), neighbourList%iNeighbour,&
            & nNeighbourSK, orb%mOrb, iCellVec, cellVec, denseDesc%iAtomStart, iSparseStart,&
            & img2CentCell, orb, species, coord)
      else
        call packHS(ERhoPrim, work, kPoint(:,iK), kWeight(iK), neighbourList%iNeighbour,&
            & nNeighbourSK, orb%mOrb, iCellVec, cellVec, denseDesc%iAtomStart, iSparseStart,&
            & img2CentCell)
      end if
    #:endif
    end do

  #:if WITH_SCALAPACK
    ! Add up and distribute energy weighted density matrix contribution from each group
    call mpifx_allreduceip(env%mpi%globalComm, ERhoPrim, MPI_SUM)
  #:endif

  end subroutine getEDensityMtxFromComplexEigvecs


  !> Calculates density matrix from Pauli-type two component eigenvectors.
  subroutine getEDensityMtxFromPauliEigvecs(env, denseDesc, filling, eigen, kPoint, kWeight,&
      & neighbourList, nNeighbourSK, orb, iSparseStart, img2CentCell, iCellVec, cellVec, tRealHS,&
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
    type(TNeighbourList), intent(in) :: neighbourList

    !> Number of neighbours for each of the atoms
    integer, intent(in) :: nNeighbourSK(:)

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
          & neighbourList%iNeighbour, nNeighbourSK, orb%mOrb, iCellVec, cellVec, iSparseStart,&
          & img2CentCell, ERhoPrim)
    #:else
      call makeDensityMatrix(work, eigvecsCplx(:,:,iKS), filling(:,iK,1), eigen(:,iK,1))
      if (tRealHS) then
        call packERho(ERhoPrim, work, neighbourList%iNeighbour, nNeighbourSK, orb%mOrb,&
            & denseDesc%iAtomStart, iSparseStart, img2CentCell)
      else
        call packERho(ERhoPrim, work, kPoint(:,iK), kWeight(iK), neighbourList%iNeighbour,&
            & nNeighbourSK, orb%mOrb, iCellVec, cellVec, denseDesc%iAtomStart, iSparseStart,&
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
  subroutine getGradients(env, sccCalc, tExtField, tXlbomd, nonSccDeriv, EField, rhoPrim, ERhoPrim,&
      & qOutput, q0, skHamCont, skOverCont, pRepCont, neighbourList, nNeighbourSK, nNeighbourRep,&
      & species, img2CentCell, iSparseStart, orb, potential, coord, derivs, iRhoPrim, thirdOrd,&
      & qDepExtPot, chrgForces, dispersion, rangeSep, SSqrReal, over, denseDesc, deltaRhoOutSqr,&
      & tPoisson, halogenXCorrection, tHelical, coord0)

    use dftbp_quaternions, only : rotate3

    !> Environment settings
    type(TEnvironment), intent(in) :: env

    !> SCC module internal variables
    type(TScc), allocatable, intent(in) :: sccCalc

    !> external electric field
    logical, intent(in) :: tExtField

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
    type(TNeighbourList), intent(in) :: neighbourList

    !> Number of neighbours for each of the atoms
    integer, intent(in) :: nNeighbourSK(:)

    !> Number of neighbours for each of the atoms closer than the repulsive cut-off
    integer, intent(in) :: nNeighbourRep(:)

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

    !> Population dependant external potential
    type(TQDepExtPotProxy), intent(inout), allocatable :: qDepExtPot

    !> forces on external charges
    real(dp), intent(inout), allocatable :: chrgForces(:,:)

    !> dispersion interactions
    class(DispersionIface), intent(inout), allocatable :: dispersion

    !> Data from rangeseparated calculations
    type(RangeSepFunc), intent(inout), allocatable ::rangeSep

    !> dense overlap matrix, required for rangeSep
    real(dp), intent(inout), allocatable :: SSqrReal(:,:)

    !> sparse overlap matrix, required for rangeSep
    real(dp), intent(in) :: over(:)

    !> Dense matrix descriptor,required for rangeSep
    type(TDenseDescr), intent(in) :: denseDesc

    !> Change in density matrix during this SCC step for rangesep
    real(dp), pointer, intent(in) :: deltaRhoOutSqr(:,:,:)

    !> whether Poisson solver is used
    logical, intent(in) :: tPoisson

    !> Correction for halogen bonds
    type(THalogenX), allocatable, intent(inout) :: halogenXCorrection

    !> Is the geometry helical
    logical, intent(in) :: tHelical

    !> Central cell atomic coordinates
    real(dp), intent(in) :: coord0(:,:)

    ! Locals
    real(dp), allocatable :: tmpDerivs(:,:)
    real(dp), allocatable :: dummyArray(:,:)
    real(dp), allocatable :: dQ(:,:,:)
    logical :: tImHam, tExtChrg, tSccCalc
    integer :: nAtom, iAt
    integer :: ii

    real(dp), parameter :: zAxis(3) = (/0.0_dp,0.0_dp,1.0_dp/)


    tSccCalc = allocated(sccCalc)
    tImHam = allocated(iRhoPrim)
    tExtChrg = allocated(chrgForces)
    nAtom = size(derivs, dim=2)

    allocate(tmpDerivs(3, nAtom))
    if (tPoisson) allocate(dummyArray(orb%mshell, nAtom))
    derivs(:,:) = 0.0_dp

    if (.not. (tSccCalc .or. tExtField)) then
      ! No external or internal potentials
      if (tImHam) then
        call derivative_shift(env, derivs, nonSccDeriv, rhoPrim, iRhoPrim, ERhoPrim, skHamCont,&
            & skOverCont, coord, species, neighbourList%iNeighbour, nNeighbourSK, img2CentCell,&
            & iSparseStart, orb, potential%intBlock, potential%iorbitalBlock)
      else
        call derivative_shift(env, derivs, nonSccDeriv, rhoPrim(:,1), ERhoPrim, skHamCont,&
            & skOverCont, coord, species, neighbourList%iNeighbour, nNeighbourSK, img2CentCell,&
            & iSparseStart, orb, tHelical)
      end if
    else
      if (tImHam) then
        call derivative_shift(env, derivs, nonSccDeriv, rhoPrim, iRhoPrim, ERhoPrim, skHamCont,&
            & skOverCont, coord, species, neighbourList%iNeighbour, nNeighbourSK, img2CentCell,&
            & iSparseStart, orb, potential%intBlock, potential%iorbitalBlock)
      else
        call derivative_shift(env, derivs, nonSccDeriv, rhoPrim, ERhoPrim, skHamCont, skOverCont,&
            & coord, species, neighbourList%iNeighbour, nNeighbourSK, img2CentCell, iSparseStart,&
            & orb, potential%intBlock)
      end if

      if (tPoisson) then
        tmpDerivs = 0.0_dp
      #:if WITH_TRANSPORT
        call poiss_getshift(dummyArray, tmpDerivs)
      #:endif
        derivs(:,:) = derivs + tmpDerivs
      else

        if (tExtChrg) then
          chrgForces(:,:) = 0.0_dp
          if (tXlbomd) then
            call error("XLBOMD does not work with external charges yet!")
          else
            call sccCalc%addForceDc(env, derivs, species, neighbourList%iNeighbour, img2CentCell,&
                & chrgForces)
          end if
        else if (tSccCalc) then
          if (tXlbomd) then
            call sccCalc%addForceDcXlbomd(env, species, orb, neighbourList%iNeighbour,&
                & img2CentCell, qOutput, q0, derivs)
          else
            call sccCalc%addForceDc(env, derivs, species, neighbourList%iNeighbour, img2CentCell)

          end if
        endif

        if (allocated(thirdOrd)) then
          if (tXlbomd) then
            call thirdOrd%addGradientDcXlbomd(neighbourList, species, coord, img2CentCell, qOutput,&
                & q0, orb, derivs)
          else
            call thirdOrd%addGradientDc(neighbourList, species, coord, img2CentCell, derivs)
          end if
        end if

      if (allocated(qDepExtPot)) then
        allocate(dQ(orb%mShell, nAtom, size(qOutput, dim=3)))
        call getChargePerShell(qOutput, orb, species, dQ, qRef=q0)
        call qDepExtPot%addGradientDc(sum(dQ(:,:,1), dim=1), dQ(:,:,1), derivs)
      end if

      if (tExtField) then
        do iAt = 1, nAtom
          derivs(:, iAt) = derivs(:, iAt)&
              & + sum(qOutput(:, iAt, 1) - q0(:, iAt, 1)) * potential%extGrad(:, iAt)
          end do
        end if

      end if
    end if

    if (allocated(dispersion)) then
      call dispersion%addGradients(derivs)
    end if

    if (allocated(halogenXCorrection)) then
      call halogenXCorrection%addGradients(derivs, coord, species, neighbourList, img2CentCell)
    end if

    if (allocated(rangeSep)) then
      if (tHelical) then
        call unpackHS(SSqrReal, over, neighbourList%iNeighbour, nNeighbourSK, denseDesc%iAtomStart,&
            & iSparseStart, img2CentCell, orb, species, coord)
      else
        call unpackHS(SSqrReal, over, neighbourList%iNeighbour, nNeighbourSK, denseDesc%iAtomStart,&
            & iSparseStart, img2CentCell)
      end if
      if (size(deltaRhoOutSqr, dim=3) > 2) then
        call error("Range separated forces do not support non-colinear spin")
      else
        call rangeSep%addLRGradients(derivs, nonSccDeriv, deltaRhoOutSqr, skHamCont, skOverCont,&
            & coord, species, orb, denseDesc%iAtomStart, SSqrReal, neighbourList%iNeighbour,&
            & nNeighbourSK)
      end if
    end if

    call getERepDeriv(tmpDerivs, coord, nNeighbourRep, neighbourList%iNeighbour, species, pRepCont,&
        & img2CentCell, tHelical)

    derivs(:,:) = derivs + tmpDerivs

    if (tHelical) then
      ! correct for folding to unit cell
      do ii = 1, nAtom
        call rotate3( derivs(:,ii), -atan2(coord(2,ii),coord(1,ii))&
            & +atan2(coord0(2,ii),coord0(1,ii)), zAxis)
      end do
    end if

  end subroutine getGradients


  !> use plumed to update derivatives
  subroutine updateDerivsByPlumed(env, plumedCalc, nAtom, iGeoStep, derivs, energy, coord0, mass,&
      & tPeriodic, latVecs)

    !> Environment settings
    type(TEnvironment), intent(in) :: env

    !> PLUMED calculator
    type(TPlumedCalc), allocatable, intent(inout) :: plumedCalc

    !> number of atoms
    integer, intent(in) :: nAtom

    !> steps taken during simulation
    integer, intent(in) :: iGeoStep

    !> the derivatives array
    real(dp), intent(inout), target, contiguous :: derivs(:,:)

    !> current energy
    real(dp), intent(in) :: energy

    !> current atomic positions
    real(dp), intent(in), target, contiguous :: coord0(:,:)

    !> atomic masses array
    real(dp), intent(in), target, contiguous :: mass(:)

    !> periodic?
    logical, intent(in) :: tPeriodic

    !> lattice vectors
    real(dp), intent(in), target, contiguous :: latVecs(:,:)

    if (.not. allocated(plumedCalc)) then
      return
    end if
    derivs(:,:) = -derivs
    call plumedCalc%sendCmdVal("setStep", iGeoStep)
    call plumedCalc%sendCmdPtr("setForces", derivs)
    call plumedCalc%sendCmdVal("setEnergy", energy)
    call plumedCalc%sendCmdPtr("setPositions", coord0)
    call plumedCalc%sendCmdPtr("setMasses", mass)
    if (tPeriodic) then
      call plumedCalc%sendCmdPtr("setBox", latVecs)
    end if
    call plumedCalc%sendCmdVal("calc", 0)
    derivs(:,:) = -derivs

  end subroutine updateDerivsByPlumed


  !> Calculates stress tensor and lattice derivatives.
  subroutine getStress(env, sccCalc, thirdOrd, tExtField, nonSccDeriv, rhoPrim, ERhoPrim, qOutput,&
      & q0, skHamCont, skOverCont, pRepCont, neighbourList, nNeighbourSk, nNeighbourRep, species,&
      & img2CentCell, iSparseStart, orb, potential, coord, latVec, invLatVec, cellVol, coord0,&
      & totalStress, totalLatDeriv, intPressure, iRhoPrim, dispersion, halogenXCorrection)

    !> Environment settings
    type(TEnvironment), intent(in) :: env

    !> SCC module internal variables
    type(TScc), allocatable, intent(in) :: sccCalc

    !> Is 3rd order SCC being used
    type(ThirdOrder), intent(inout), allocatable :: thirdOrd

    !> External field
    logical, intent(in) :: tExtField

    !> method for calculating derivatives of S and H0
    type(NonSccDiff), intent(in) :: nonSccDeriv

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
    type(TNeighbourList), intent(in) :: neighbourList

    !> Number of neighbours for each of the atoms
    integer, intent(in) :: nNeighbourSK(:)

    !> Number of neighbours for each of the atoms closer than the repulsive cut-off
    integer, intent(in) :: nNeighbourRep(:)

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

    !> Correction for halogen bonds
    type(THalogenX), allocatable, intent(inout) :: halogenXCorrection

    real(dp) :: tmpStress(3, 3)
    logical :: tImHam

    tImHam = allocated(iRhoPrim)

    if (allocated(sccCalc)) then
      if (tImHam) then
        call getBlockiStress(env, totalStress, nonSccDeriv, rhoPrim, iRhoPrim, ERhoPrim, skHamCont,&
            & skOverCont, coord, species, neighbourList%iNeighbour, nNeighbourSK, img2CentCell,&
            & iSparseStart, orb, potential%intBlock, potential%iorbitalBlock, cellVol)
      else
        call getBlockStress(env, totalStress, nonSccDeriv, rhoPrim, ERhoPrim, skHamCont,&
            & skOverCont, coord, species, neighbourList%iNeighbour, nNeighbourSK, img2CentCell,&
            & iSparseStart, orb, potential%intBlock, cellVol)
      end if
      call sccCalc%addStressDc(totalStress, env, species, neighbourList%iNeighbour, img2CentCell)
      if (allocated(thirdOrd)) then
        call thirdOrd%addStressDc(neighbourList, species, coord, img2CentCell, cellVol, totalStress)
      end if
    else
      if (tImHam) then
        call getBlockiStress(env, totalStress, nonSccDeriv, rhoPrim, iRhoPrim, ERhoPrim, skHamCont,&
            & skOverCont, coord, species, neighbourList%iNeighbour, nNeighbourSK, img2CentCell,&
            & iSparseStart, orb, potential%intBlock, potential%iorbitalBlock, cellVol)
      else
        call getNonSCCStress(env, totalStress, nonSccDeriv, rhoPrim(:,1), ERhoPrim, skHamCont,&
            & skOverCont, coord, species, neighbourList%iNeighbour, nNeighbourSK, img2CentCell,&
            & iSparseStart, orb, cellVol)
      end if
    end if

    if (allocated(dispersion)) then
      call dispersion%getStress(tmpStress)
      totalStress(:,:) = totalStress + tmpStress
    end if

    if (allocated(halogenXCorrection)) then
      call halogenXCorrection%getStress(tmpStress, coord, neighbourList, species, img2CentCell,&
          & cellVol)
      totalStress(:,:) = totalStress + tmpStress
    end if

    if (tExtField) then
      call getExtFieldStress(latVec, cellVol, q0, qOutput, potential%extGrad, coord0, tmpStress)
      totalStress(:,:) = totalStress + tmpStress
    end if

    call getRepulsiveStress(tmpStress, coord, nNeighbourRep, neighbourList%iNeighbour, species,&
        & img2CentCell, pRepCont, cellVol)
    totalStress(:,:) = totalStress + tmpStress

    intPressure = (totalStress(1,1) + totalStress(2,2) + totalStress(3,3)) / 3.0_dp
    totalLatDeriv(:,:) = -cellVol * matmul(totalStress, invLatVec)

  end subroutine getStress


  !> Calculates stress from external electric field.
  subroutine getExtFieldStress(latVec, cellVol, q0, qOutput, extPotGrad, coord0, stress)

    !> lattice vectors
    real(dp), intent(in) :: latVec(:,:)

    !> unit cell volume
    real(dp), intent(in) :: cellVol

    !> reference atomic charges
    real(dp), intent(in) :: q0(:,:,:)

    !> number of electrons in each orbital
    real(dp), intent(in) :: qOutput(:,:,:)

    !> Gradient of the external field
    real(dp), intent(in) :: extPotGrad(:,:)

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
              & - sum(q0(:,iAtom,1) - qOutput(:,iAtom,1), dim=1) * extPotGrad(ii, iAtom)&
              & * coord0(jj, iAtom)
        end do
      end do
    end do
    call frac2cart(coord0, latVec)
    stress(:,:) = -matmul(latDerivs, transpose(latVec)) / cellVol

  end subroutine getExtFieldStress


  !> Calculates derivatives with respect to helical bounday conditions
  subroutine getHelicalBCDerivs(env, boundaryConds, sccCalc, thirdOrd, tExtField, nonSccDeriv,&
      & rhoPrim, ERhoPrim, qOutput, q0, skHamCont, skOverCont, pRepCont, neighbourList,&
      & nNeighbourSk, nNeighbourRep, species, img2CentCell, iSparseStart, orb, potential, coord,&
      & iRhoPrim, dispersion, halogenXCorrection)

    !> Environment settings
    type(TEnvironment), intent(in) :: env

    !> Helical boundary conditions (translation and angle)
    real(dp),intent(in) :: boundaryConds(2)

    !> SCC module internal variables
    type(TScc), allocatable, intent(in) :: sccCalc

    !> Is 3rd order SCC being used
    type(ThirdOrder), intent(inout), allocatable :: thirdOrd

    !> External field
    logical, intent(in) :: tExtField

    !> method for calculating derivatives of S and H0
    type(NonSccDiff), intent(in) :: nonSccDeriv

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
    type(TNeighbourList), intent(in) :: neighbourList

    !> Number of neighbours for each of the atoms
    integer, intent(in) :: nNeighbourSK(:)

    !> Number of neighbours for each of the atoms closer than the repulsive cut-off
    integer, intent(in) :: nNeighbourRep(:)

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

    !> imaginary part of the density matrix (if present)
    real(dp), intent(in), allocatable :: iRhoPrim(:,:)

    !> dispersion interactions
    class(DispersionIface), allocatable, intent(inout) :: dispersion

    !> Correction for halogen bonds
    type(THalogenX), allocatable, intent(inout) :: halogenXCorrection

    real(dp) :: derivTotal(2), derivRep(2), derivNonSCC(2)
    logical :: tImHam

    tImHam = allocated(iRhoPrim)

    derivTotal(:) = 0.0_dp
    derivRep(:) = 0.0_dp
    derivNonSCC(:) = 0.0_dp

    if (allocated(sccCalc)) then
      call error("Currently missing")
    else
      if (tImHam) then
        call error("Currently missing")
      else
        call getHelicalNonSCCLatForce(derivNonSCC, nonSccDeriv, rhoPrim(:,1), ERhoPrim, skHamCont,&
            & skOverCont, coord, species, neighbourList%iNeighbour, nNeighbourSK, img2CentCell,&
            & iSparseStart, orb, boundaryConds)
      end if
    end if

    if (allocated(dispersion)) then
      call error("Currently missing")
    end if

    if (allocated(halogenXCorrection)) then
      call error("Currently missing")
    end if

    if (tExtField) then
      call error("Currently missing")
    end if

    call getHelicalRepLatForce(derivRep, coord, nNeighbourRep, neighbourList%iNeighbour, species,&
        & img2CentCell, pRepCont, boundaryConds)

    derivTotal(:) = derivRep + derivNonSCC

    write(stdOut,*)derivRep
    write(stdOut,*)derivNonSCC
    write(stdOut,*)'Helical BC forces', -derivTotal

  end subroutine getHelicalBCDerivs


  !> Removes forces components along constraint directions
  subroutine constrainForces(conAtom, conVec, derivs)

    !> atoms being constrained
    integer, intent(in) :: conAtom(:)

    !> vector to project out forces
    real(dp), intent(in) :: conVec(:,:)

    !> on input energy derivatives, on exit resulting projected derivatives
    real(dp), intent(inout) :: derivs(:,:)

    integer :: ii, iAtom

    ! Set force components along constraint vectors zero
    do ii = 1, size(conAtom)
      iAtom = conAtom(ii)
      derivs(:,iAtom) = derivs(:,iAtom)&
          & - conVec(:,ii) * dot_product(conVec(:,ii), derivs(:,iAtom))
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
  subroutine getNextDerivStep(derivDriver, derivs, indMovedAtoms, coord, tGeomEnd)

    !> Driver for the finite difference second derivatives
    type(OnumDerivs), intent(inout) :: derivDriver

    !> first derivatives of energy at the current coordinates
    real(dp), intent(in) :: derivs(:,:)

    !> moving atoms
    integer, intent(in) :: indMovedAtoms(:)

    !> atomic coordinates
    real(dp), intent(out) :: coord(:,:)

    !> has the process terminated
    logical, intent(out) :: tGeomEnd

    real(dp) :: newCoords(3, size(indMovedAtoms))

    call next(derivDriver, newCoords, derivs(:, indMovedAtoms), tGeomEnd)
    coord(:, indMovedAtoms) = newCoords

  end subroutine getNextDerivStep


  !> Returns the coordinates for the next coordinate optimisation step.
  subroutine getNextCoordinateOptStep(pGeoCoordOpt, energy, derivss, indMovedAtom, coord0,&
      & diffGeo, tCoordEnd, tRemoveExcitation)

    !> optimiser for atomic coordinates
    type(OGeoOpt), intent(inout) :: pGeoCoordOpt

    !> energies
    type(TEnergies), intent(in) :: energy

    !> Derivative of energy with respect to atomic coordinates
    real(dp), intent(in) :: derivss(:,:)

    !> numbers of the moving atoms
    integer, intent(in) :: indMovedAtom(:)

    !> central cell atomic coordinates
    real(dp), intent(inout) :: coord0(:,:)

    !> largest change in atomic coordinates
    real(dp), intent(out) :: diffGeo

    !> has the geometry optimisation finished
    logical, intent(out) :: tCoordEnd

    !> remove excited state energy from the step, to be consistent with the forces
    logical, intent(in) :: tRemoveExcitation

    real(dp) :: derivssMoved(3 * size(indMovedAtom))
    real(dp), target :: newCoordsMoved(3 * size(indMovedAtom))
    real(dp), pointer :: pNewCoordsMoved(:,:)

    derivssMoved(:) = reshape(derivss(:, indMovedAtom), [3 * size(indMovedAtom)])
    if (tRemoveExcitation) then
      call next(pGeoCoordOpt, energy%EForceRelated, derivssMoved, newCoordsMoved, tCoordEnd)
    else
      call next(pGeoCoordOpt, energy%EForceRelated + energy%Eexcited, derivssMoved, newCoordsMoved,&
          & tCoordEnd)
    end if
    pNewCoordsMoved(1:3, 1:size(indMovedAtom)) => newCoordsMoved(1 : 3 * size(indMovedAtom))
    diffGeo = maxval(abs(pNewCoordsMoved - coord0(:, indMovedAtom)))
    coord0(:, indMovedAtom) = pNewCoordsMoved

  end subroutine getNextCoordinateOptStep


  !> Returns the coordinates and lattice vectors for the next lattice optimisation step.
  subroutine getNextLatticeOptStep(pGeoLatOpt, energy, constrLatDerivs, origLatVec, tLatOptFixAng,&
      & tLatOptFixLen, tLatOptIsotropic, indMovedAtom, latVec, coord0, diffGeo, tGeomEnd)

    !> lattice vector optimising object
    type(OGeoOpt), intent(inout) :: pGeoLatOpt

    !> Energy contributions and total
    type(TEnergies), intent(inout) :: energy

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

    call next(pGeoLatOpt, energy%EForceRelated, constrLatDerivs, newLatVecsFlat,tGeomEnd)
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
    energy%EForceRelated = energy%EForceRelated + energy%Ekin

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
      & kPoint, neighbourList, nNeighbourSK, denseDesc,  iSparseStart, img2CentCell, iCellVec,&
      & cellVec, runId, orb, species, speciesName, parallelKS, localisation, eigvecsReal, SSqrReal,&
      & eigvecsCplx, SSqrCplx, tHelical, coord)

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
    type(TNeighbourList), intent(in) :: neighbourList

    !> Number of neighbours for each of the atoms
    integer, intent(in) :: nNeighbourSK(:)

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

    !> Is the geometry helical
    logical, intent(in) :: tHelical

    !> atomic coordinates
    real(dp), intent(in) :: coord(:,:)

    integer :: nFilledLev, nAtom, nSpin
    integer :: iSpin, iKS, iK

    nAtom = size(orb%nOrbAtom)
    nSpin = size(nEl)

    if (any(abs(mod(filling, real(3 - nSpin, dp))) > elecTolMax)) then
      call warning("Fractional occupations allocated for electron localisation")
    end if

    if (allocated(eigvecsReal)) then
      if (tHelical) then
        call unpackHS(SSqrReal,over,neighbourList%iNeighbour, nNeighbourSK, denseDesc%iAtomStart,&
            & iSparseStart, img2CentCell, orb, species, coord)
      else
        call unpackHS(SSqrReal,over,neighbourList%iNeighbour, nNeighbourSK, denseDesc%iAtomStart,&
            & iSparseStart, img2CentCell)
      end if
      do iKS = 1, parallelKS%nLocalKS
        iSpin = parallelKS%localKS(2, iKS)
        nFilledLev = nint(nEl(iSpin) / real(3 - nSpin, dp))
        localisation = pipekMezey%getLocalisation(eigvecsReal(:, 1:nFilledLev, iKS), SSqrReal,&
            & denseDesc%iAtomStart)
        write(stdOut, "(A, E15.8)") 'Original localisation', localisation
        call pipekMezey%calcCoeffs(eigvecsReal(:, 1:nFilledLev, iKS), SSqrReal,&
            & denseDesc%iAtomStart)
        localisation = pipekMezey%getLocalisation(eigvecsReal(:,1:nFilledLev,iKS), SSqrReal,&
            & denseDesc%iAtomStart)
        write(stdOut, "(A, E20.12)") 'Final localisation ', localisation
      end do

      call writeRealEigvecs(env, runId, neighbourList, nNeighbourSK, denseDesc, iSparseStart,&
          & img2CentCell, species(:nAtom), speciesName, orb, over, parallelKS, tPrintEigvecsTxt,&
          & eigvecsReal, SSqrReal, fileName="localOrbs")
    else

      localisation = 0.0_dp
      do iKS = 1, parallelKS%nLocalKS
        iK = parallelKS%localKS(1, iKS)
        iSpin = parallelKS%localKS(2, iKS)
        nFilledLev = nint(nEl(iSpin) / real( 3 - nSpin, dp))
        localisation = localisation + pipekMezey%getLocalisation(&
            & eigvecsCplx(:,:nFilledLev,iKS), SSqrCplx, over, kpoint(:,iK), neighbourList,&
            & nNeighbourSK, iCellVec, cellVec, denseDesc%iAtomStart, iSparseStart, img2CentCell)
      end do
      write(stdOut, "(A, E20.12)") 'Original localisation', localisation

      ! actual localisation calls
      do iKS = 1, parallelKS%nLocalKS
        iK = parallelKS%localKS(1, iKS)
        iSpin = parallelKS%localKS(2, iKS)
        nFilledLev = nint(nEl(iSpin) / real( 3 - nSpin, dp))
        call pipekMezey%calcCoeffs(eigvecsCplx(:,:nFilledLev,iKS), SSqrCplx, over, kpoint(:,iK),&
            & neighbourList, nNeighbourSK, iCellVec, cellVec, denseDesc%iAtomStart, iSparseStart,&
            & img2CentCell)
      end do

      localisation = 0.0_dp
      do iKS = 1, parallelKS%nLocalKS
        iK = parallelKS%localKS(1, iKS)
        iSpin = parallelKS%localKS(2, iKS)
        nFilledLev = nint(nEl(iSpin) / real( 3 - nSpin, dp))
        localisation = localisation + pipekMezey%getLocalisation(&
            & eigvecsCplx(:,:nFilledLev,iKS), SSqrCplx, over, kpoint(:,iK), neighbourList,&
            & nNeighbourSK, iCellVec, cellVec, denseDesc%iAtomStart, iSparseStart, img2CentCell)
      end do
      write(stdOut, "(A, E20.12)") 'Final localisation', localisation

      call writeCplxEigvecs(env, runId, neighbourList, nNeighbourSK, cellVec, iCellVec, denseDesc,&
          & iSparseStart, img2CentCell, species, speciesName, orb, kPoint, over, parallelKS,&
          & tPrintEigvecsTxt, eigvecsCplx, SSqrCplx, fileName="localOrbs")

    end if

  end subroutine calcPipekMezeyLocalisation

  subroutine printMaxForces(derivs, constrLatDerivs, tCoordOpt, tLatOpt, indMovedAtoms)
    real(dp), intent(in), allocatable :: derivs(:,:)
    real(dp), intent(in) :: constrLatDerivs(:)
    logical, intent(in) :: tCoordOpt
    logical, intent(in) :: tLatOpt
    integer, intent(in) :: indMovedAtoms(:)

    if (tCoordOpt) then
      call printMaxForce(maxval(abs(derivs(:, indMovedAtoms))))
    end if
    if (tLatOpt) then
      call printMaxLatticeForce(maxval(abs(constrLatDerivs)))
    end if

  end subroutine printMaxForces


#:if WITH_SOCKETS

  subroutine sendEnergyAndForces(env, socket, energy, TS, derivs, totalStress, cellVol)
    type(TEnvironment), intent(in) :: env
    ! Socket may be unallocated (as on slave processes)
    type(IpiSocketComm), allocatable, intent(inout) :: socket
    type(TEnergies), intent(in) :: energy
    real(dp), intent(in) :: TS(:)
    real(dp), intent(in) :: derivs(:,:)
    real(dp), intent(in) :: totalStress(:,:)
    real(dp), intent(in) :: cellVol

    if (env%tGlobalMaster) then
      ! stress was computed above in the force evaluation block or is 0 if aperiodic
      call socket%send(energy%ETotal - sum(TS), -derivs, totalStress * cellVol)
    end if
  end subroutine sendEnergyAndForces

#:endif

end module dftbp_main
