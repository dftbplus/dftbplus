!--------------------------------------------------------------------------------------------------!
!  DFTB+: general package for performing fast atomistic simulations                                !
!  Copyright (C) 2006 - 2022  DFTB+ developers group                                               !
!                                                                                                  !
!  See the LICENSE file for terms of usage and distribution.                                       !
!--------------------------------------------------------------------------------------------------!

#:include 'common.fypp'
#:include 'error.fypp'

!> The main routines for DFTB+
module dftbp_dftbplus_main
  use dftbp_common_accuracy, only : dp, elecTolMax, tolSameDist
  use dftbp_common_constants, only : pi
  use dftbp_common_environment, only : TEnvironment, globalTimers
  use dftbp_common_globalenv, only : stdOut, withMpi
  use dftbp_common_hamiltoniantypes, only : hamiltonianTypes
  use dftbp_common_status, only : TStatus
  use dftbp_derivs_numderivs2, only : TNumderivs, next, getHessianMatrix
  use dftbp_derivs_perturb, only : TResponse
  use dftbp_dftb_blockpothelper, only : appendBlockReduced
  use dftbp_dftb_boundarycond, only : TBoundaryConditions
  use dftbp_dftb_densitymatrix, only : makeDensityMatrix
  use dftbp_dftb_determinants, only : TDftbDeterminants, TDftbDeterminants_init, determinants
  use dftbp_dftb_dftbplusu, only : TDftbU
  use dftbp_dftb_dispersions, only : TDispersionIface
  use dftbp_dftb_elstatpot, only : TElStatPotentials
  use dftbp_dftb_energytypes, only : TEnergies
  use dftbp_dftb_etemp, only : electronFill, Efilling
  use dftbp_dftb_extfields, only : addUpExternalField
  use dftbp_dftb_forces, only : derivative_shift
  use dftbp_dftb_getenergies, only : calcEnergies, calcDispersionEnergy, sumEnergies
  use dftbp_dftb_halogenx, only : THalogenX
  use dftbp_dftb_hamiltonian, only : resetInternalPotentials, addChargePotentials,&
      & getSccHamiltonian, mergeExternalPotentials, resetExternalPotentials,&
      & addBlockChargePotentials
  use dftbp_dftb_nonscc, only : TNonSccDiff, buildS, buildH0
  use dftbp_dftb_onsitecorrection, only : Onsblock_expand, onsBlock_reduce, addOnsShift
  use dftbp_dftb_orbitalequiv, only : OrbitalEquiv_expand, orbitalEquiv_reduce
  use dftbp_dftb_periodic, only : TNeighbourList, updateNeighbourListAndSpecies, cart2frac,&
      & frac2cart, getNrOfNeighboursForAll, getCellTranslations
  use dftbp_dftb_pmlocalisation, only : TPipekMezey
  use dftbp_dftb_populations, only : getChargePerShell, denseSubtractDensityOfAtoms, mulliken,&
      & denseMulliken, denseBlockMulliken, skewMulliken, getOnsitePopulation, &
      & getAtomicMultipolePopulation
  use dftbp_dftb_potentials, only : TPotentials
  use dftbp_dftb_rangeseparated, only : TRangeSepFunc
  use dftbp_dftb_repulsive_repulsive, only : TRepulsive
  use dftbp_dftb_scc, only : TScc
  use dftbp_dftb_shift, only : addShift
  use dftbp_dftb_slakocont, only : TSlakoCont
  use dftbp_dftb_sparse2dense, only : unpackHPauli, unpackHS, blockSymmetrizeHS, packHS,&
      & blockSymmetrizeHS, packHS, SymmetrizeHS, unpackHelicalHS, packerho, blockHermitianHS,&
      & packHSPauli, packHelicalHS, packHSPauliImag, iPackHS, unpackSPauli, getSparseDescriptor
  use dftbp_dftb_spin, only : ud2qm, qm2ud
  use dftbp_dftb_spinorbit, only : addOnsiteSpinOrbitHam, getOnsiteSpinOrbitEnergy
  use dftbp_dftb_stress, only : getkineticstress, getBlockStress, getBlockiStress, getNonSCCStress
  use dftbp_dftb_thirdorder, only : TThirdOrder
  use dftbp_dftbplus_eigenvects, only : diagDenseMtx
  use dftbp_dftbplus_elstattypes, only : elstatTypes
  use dftbp_dftbplus_forcetypes, only : forceTypes
  use dftbp_dftbplus_initprogram, only : TDftbPlusMain, TCutoffs, TNegfInt
  use dftbp_dftbplus_inputdata, only : TNEGFInfo
  use dftbp_dftbplus_mainio, only : writeRealEigvecs, writeCplxEigVecs, readEigenVecs,&
      & printMaxForce, printMaxLatticeForce, printReksSccHeader, printSccHeader, printMdInfo,&
      & writeMdOut2, writeDetailedOut5, writeMdOut1, openOutputFile, printReksSccInfo,&
      & writeReksDetailedOut1, writebandout, writehsandstop, printSccInfo, printBlankLine,&
      & writeCharges, writeDetailedOut1, writeDetailedOut2, writeDetailedOut3,&
      & writeEigenVectors, writeProjectedEigenvectors, writeCurrentGeometry, writeDetailedOut4,&
      & writeEsp, printGeostepInfo, writeDetailedOut2dets, printEnergies, printVolume,&
      & printPressureAndFreeEnergy, writeDetailedOut6, writeDetailedOut7,&
      & writeFinalDriverstatus, writeHessianout, writeAutotestTag, writeResultsTag,&
      & writeDetailedXml, writeCosmoFile, printForceNorm, printLatticeForceNorm, writeDerivBandOut,&
      & writeDetailedOut8, writeDetailedOut9
  use dftbp_dftbplus_outputfiles, only : autotestTag, bandOut, fCharges, fShifts, fStopScc, mdOut,&
      & userOut, fStopDriver, hessianOut, resultsTag, derivEBandOut
  use dftbp_dftbplus_qdepextpotproxy, only : TQDepExtPotProxy
  use dftbp_dftbplus_transportio, only : readShifts, writeShifts, writeContactShifts
  use dftbp_elecsolvers_elecsolvers, only : TElectronicSolver
  use dftbp_elecsolvers_elecsolvertypes, only : electronicSolverTypes
  use dftbp_extlibs_plumed, only : TPlumedCalc, TPlumedCalc_final
  use dftbp_extlibs_tblite, only : TTBLite
  use dftbp_geoopt_geoopt, only : TGeoOpt, next, reset
  use dftbp_io_message, only : error, warning
  use dftbp_io_taggedoutput, only : TTaggedWriter
  use dftbp_math_angmomentum, only : getLOnsite, getLDual
  use dftbp_math_blasroutines, only : hemm, symm
  use dftbp_math_lapackroutines, only : hermatinv, matinv, symmatinv
  use dftbp_math_simplealgebra, only : determinant33, derivDeterminant33
  use dftbp_md_mdcommon, only : TMdCommon, evalKE, evalKT
  use dftbp_md_mdintegrator, only : TMdIntegrator, next, rescale
  use dftbp_md_tempprofile, only : TTempProfile
  use dftbp_md_xlbomd, only : TXLBOMD
  use dftbp_mixer_mixer, only : TMixer, reset, mix, getInverseJacobian
  use dftbp_reks_reks, only : TReksCalc, guessneweigvecs, optimizeFONs, calcweights, activeorbswap,&
      & getfilling, calcsareksenergy, printsareksenergy, qm2udl, printreksmicrostates, qmexpandl,&
      & ud2qml, constructmicrostates, checkgammapoint, getfockanddiag, printrekssainfo,&
      & getstateinteraction, getreksenproperties, getreksgradients, getreksgradproperties,&
      & getReksStress
  use dftbp_solvation_cm5, only : TChargeModel5
  use dftbp_solvation_fieldscaling, only : TScaleExtEField
  use dftbp_solvation_solvation, only : TSolvation
  use dftbp_timedep_linresp, only : TLinResp, linResp_calcExcitations, LinResp_addGradients
  use dftbp_timedep_pprpa, only : ppRpaEnergies
  use dftbp_timedep_timeprop, only : runDynamics
  use dftbp_type_commontypes, only : TOrbitals, TParallelKS
  use dftbp_type_densedescr, only : TDenseDescr
  use dftbp_type_integral, only : TIntegral
  use dftbp_type_multipole, only : TMultipole
#:if WITH_SCALAPACK
  use dftbp_dftb_densitymatrix, only : makeDensityMtxRealBlacs, makeDensityMtxCplxBlacs
  use dftbp_dftb_sparse2dense, only : packRhoRealBlacs, packRhoCplxBlacs, packRhoPauliBlacs,&
      & packRhoHelicalRealBlacs, packRhoHelicalCplxBlacs, packERhoPauliBlacs, unpackHSRealBlacs,&
      & unpackHSCplxBlacs, unpackHPauliBlacs, unpackSPauliBlacs, unpackHSHelicalRealBlacs,&
      & unpackHSHelicalCplxBlacs
  use dftbp_dftbplus_eigenvects, only : diagDenseMtxBlacs
  use dftbp_extlibs_mpifx, only : MPI_SUM, mpifx_allreduceip
  use dftbp_extlibs_scalapackfx, only : pblasfx_phemm, pblasfx_psymm, pblasfx_ptran,&
      & pblasfx_ptranc, blacsfx_gemr2d
  use dftbp_math_scalafxext, only : phermatinv, psymmatinv
#:endif
#:if WITH_SOCKETS
  use dftbp_io_ipisocket, only : IpiSocketComm
  use dftbp_dftbplus_mainio, only : receiveGeometryFromSocket
#:endif
#:if WITH_TRANSPORT
  use dftbp_dftbplus_initprogram, only : overrideContactCharges
  use dftbp_transport_negfint, only : TNegfInt_final
  use dftbp_transport_negfvars, only : TTransPar
#:endif
  implicit none

  private
  public :: runDftbPlus
  public :: processGeometry

  !> Should further output be appended to detailed.out?
  logical, parameter :: tAppendDetailedOut = .false.


contains

  !> The main DFTB program itself
  subroutine runDftbPlus(this, env)

    !> Global variables
    type(TDftbPlusMain), intent(inout) :: this

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

    !> Should the geometry loop be stopped early
    logical :: tExitGeoOpt

    !> Which state is being calculated in the determinant loop?
    integer :: iDet
    logical :: isUnReduced

    type(TStatus) :: errStatus

    call initGeoOptParameters(this%tCoordOpt, this%nGeoSteps, tGeomEnd, tCoordStep, tStopDriver,&
        & iGeoStep, iLatGeoStep)

    ! If the geometry is periodic, need to update lattice information in geometry loop
    this%tLatticeChanged = this%tPeriodic

    ! As first geometry iteration, require updates for coordinates in dependent routines
    this%tCoordsChanged = .true.

    ! Main geometry loop
    geoOpt: do iGeoStep = 0, this%nGeoSteps
      tWriteRestart = env%tGlobalLead .and. needsRestartWriting(this%isGeoOpt .or.&
          & allocated(this%geoOpt), this%tMd, iGeoStep, this%nGeoSteps, this%restartFreq)

      if (.not. this%tRestartNoSC) then
        call printGeoStepInfo(this%tCoordOpt, this%tLatOpt, iLatGeoStep, iGeoStep)
      end if

      ! DFTB Determinant Loop
      ! Will pass though loop once, unless specified in input to perform multiple determinants
      lpDets : do iDet = 1, this%nDets

        this%deltaDftb%iDeterminant = iDet

        call preDetCharges(isUnReduced, iDet, this%nDets, iGeoStep, this%deltaDftb, this%qInput,&
            & this%qDets, this%qBlockIn, this%qBlockDets, this%deltaRhoIn, this%deltaRhoDets)
        if (isUnReduced) then
          call reduceCharges(this%orb, this%nIneqOrb, this%iEqOrbitals, this%qInput, this%qInpRed,&
              & this%qBlockIn, this%dftbU, this%iEqBlockDftbu, this%qiBlockIn,&
              & this%iEqBlockDftbuLS, this%iEqBlockOnSite, this%iEqBlockOnSiteLS)
        end if

        call processGeometry(this, env, iGeoStep, iLatGeoStep, tWriteRestart, tStopScc,&
            & tExitGeoOpt, errStatus)
        if (errStatus%hasError()) then
          call error(errStatus%message)
        end if

        call postDetCharges(iDet, this%nDets, this%qOutput, this%qDets, this%qBlockDets,&
            & this%qBlockOut, this%deltaRhoDets, this%deltaRhoOut)

      end do lpDets

      call this%deltaDftb%postProcessDets(this%dftbEnergy, this%qOutput, this%qDets,&
          & this%qBlockOut, this%qBlockDets, this%dipoleMoment, this%totalStress,&
          & this%tripletStress, this%mixedStress, this%derivs, this%tripletderivs, this%mixedderivs)

      if (this%tWriteDetailedOut .and. this%deltaDftb%nDeterminant() > 1) then
        call writeDetailedOut2Dets(this%fdDetailedOut, userOut, tAppendDetailedOut,&
            & this%dftbEnergy, this%electronicSolver, this%deltaDftb, this%q0, this%orb,&
            & this%qOutput, this%qDets, this%qBlockDets, this%species, this%iAtInCentralRegion,&
            & this%tPrintMulliken, this%cm5Cont)
      end if

      if (.not.this%tRestartNoSC .and.&
          & this%electronicSolver%iSolver /= electronicSolverTypes%OnlyTransport) then
        call printEnergies(this%dftbEnergy, this%electronicSolver, this%deltaDftb)
      end if

      if (this%tStress) then

        call printVolume(this%cellVol)

        ! MD case includes the atomic kinetic energy contribution, so print that later
        if (.not. (this%tMD .or. this%tHelical)) then
          call printPressureAndFreeEnergy(this%extPressure, this%intPressure,&
              & this%dftbEnergy(this%deltaDftb%iDeterminant)%EGibbs)
        end if
      end if

      call postprocessDerivs(this%derivs, this%conAtom, this%conVec, this%tLatOpt,&
          & this%totalLatDeriv, this%extLatDerivs, this%normOrigLatVec, this%tLatOptFixAng,&
          & this%tLatOptFixLen, this%tLatOptIsotropic, constrLatDerivs)

      if (tExitGeoOpt) then
        exit geoOpt
      end if

      call printMaxForces(this%derivs, constrLatDerivs, this%tCoordOpt, this%tLatOpt,&
          & this%indMovedAtom)
    #:if WITH_SOCKETS
      if (this%tSocket) then
        call sendEnergyAndForces(env, this%socket, this%dftbEnergy(this%deltaDftb%iFinal),&
            & this%derivs, this%totalStress, this%cellVol)
      end if
    #:endif
      tWriteCharges =  allocated(this%qInput) .and. tWriteRestart .and. this%tMulliken&
          & .and. this%tSccCalc .and. .not. this%tDerivs&
          & .and. this%maxSccIter > 1 .and. this%deltaDftb%nDeterminant() == 1
      if (tWriteCharges) then
        call writeCharges(fCharges, this%tWriteChrgAscii, this%orb, this%qInput, this%qBlockIn,&
            & this%qiBlockIn, this%deltaRhoIn, size(this%iAtInCentralRegion), this%multipoleInp)
      end if

      if (this%tForces) then
        call getNextGeometry(this, env, iGeoStep, tWriteRestart, constrLatDerivs, tCoordStep,&
            & tGeomEnd, tStopDriver, iLatGeoStep, tempIon, tExitGeoOpt)
        if (tExitGeoOpt) then
          exit geoOpt
        end if
      end if

      if (this%tWriteDetailedOut .and. this%tMd) then
        call writeDetailedOut6(this%fdDetailedOut%unit, this%dftbEnergy(this%deltaDftb%iFinal),&
            & tempIon)
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
    if (this%tSocket .and. env%tGlobalLead) then
      call this%socket%shutdown()
    end if
  #:endif

    if (allocated(this%plumedCalc)) then
      call TPlumedCalc_final(this%plumedCalc)
    end if

    tGeomEnd = this%tMD .or. tGeomEnd .or. this%tDerivs

    if (env%tGlobalLead) then
      if (this%tWriteDetailedOut) then
        call writeDetailedOut7(this%fdDetailedOut%unit,&
            & this%isGeoOpt .or. allocated(this%geoOpt), tGeomEnd, this%tMd, this%tDerivs,&
            & this%eField, this%dipoleMoment, this%deltaDftb, this%eFieldScaling,&
            & this%dipoleMessage)
      end if

      call writeFinalDriverStatus(this%isGeoOpt .or. allocated(this%geoOpt), tGeomEnd, this%tMd,&
          & this%tDerivs)

      if (this%tMD) then
        deallocate(this%fdMd)
        write(stdOut, "(2A)") 'MD information accumulated in ', mdOut
      end if
    end if

    if (env%tGlobalLead .and. this%tDerivs) then
      call getHessianMatrix(this%derivDriver, this%pDynMatrix)
      call writeHessianOut(hessianOut, this%pDynMatrix, this%indMovedAtom)
    else
      nullify(this%pDynMatrix)
    end if

    if (this%tWriteShifts) then
      call writeShifts(fShifts, this%orb, this%potential%intShell)
    endif

    ! Here time propagation is called
    if (allocated(this%electronDynamics)) then
      call runDynamics(this%electronDynamics, this%boundaryCond, this%eigvecsReal, this%H0,&
          & this%species, this%q0, this%referenceN0, this%ints, this%filling, this%neighbourList,&
          & this%nNeighbourSK, this%nNeighbourLC, this%denseDesc%iAtomStart, this%iSparseStart,&
          & this%img2CentCell, this%orb, this%coord0, this%spinW, this%repulsive, env,&
          & this%tDualSpinOrbit, this%xi, this%thirdOrd, this%solvation, this%eFieldScaling,&
          & this%rangeSep, this%qDepExtPot, this%dftbU, this%iAtInCentralRegion, this%tFixEf,&
          & this%Ef, this%coord, this%onsiteElements, this%skHamCont, this%skOverCont, this%latVec,&
          & this%invLatVec, this%iCellVec, this%rCellVec, this%cellVec, this%electronicSolver,&
          & this%eigvecsCplx, this%taggedWriter, this%refExtPot, errStatus)
      if (errStatus%hasError()) then
        call error(errStatus%message)
      end if
    end if

  #:if WITH_TRANSPORT

    if (this%isAContactCalc) then
      ! Note: shift and charges are saved in QM representation (not UD)
      call writeContactShifts(this%transpar%contacts(this%transpar%taskContInd)%name, this%orb,&
          & this%potential%coulombShell, this%qOutput, this%Ef, this%qBlockOut,&
          & .not.this%transpar%tWriteBinShift)
    end if

    if (this%tLocalCurrents) then
      call this%negfInt%local_currents(env, this%parallelKS%localKS, this%ints%hamiltonian,&
          & this%ints%overlap, this%neighbourList, this%nNeighbourSK, this%cutOff%skCutoff,&
          & this%denseDesc%iAtomStart, this%iSparseStart, this%img2CentCell, this%iCellVec,&
          & this%cellVec, this%rCellVec, this%orb, this%kPoint, this%kWeight, this%coord0Fold,&
          & this%species0, this%speciesName, this%mu, this%lCurrArray)
    end if

    if (this%tTunn) then
      call this%negfInt%calc_current(env, this%parallelKS%localKS, this%ints%hamiltonian,&
          & this%ints%overlap, this%neighbourList%iNeighbour, this%nNeighbourSK,&
          & this%densedesc%iAtomStart, this%iSparseStart, this%img2CentCell, this%iCellVec,&
          & this%cellVec, this%orb, this%kPoint, this%kWeight, this%tunneling, this%current,&
          & this%ldos, this%leadCurrents, this%writeTunn, this%tWriteLDOS, this%regionLabelLDOS,&
          & this%mu)
    end if

  #:endif

    if (this%isDFTBPT) then

      if (this%isEResp) then
        call this%response%wrtEField(env, this%parallelKS, this%filling, this%eigen,&
            & this%eigVecsReal, this%eigvecsCplx, this%ints%hamiltonian, this%ints%overlap,&
            & this%orb, this%nAtom, this%species, this%neighbourList, this%nNeighbourSK,&
            & this%denseDesc, this%iSparseStart, this%img2CentCell, this%coord, this%scc,&
            & this%maxSccIter, this%sccTol, this%isSccConvRequired, this%nMixElements,&
            & this%nIneqOrb, this%iEqOrbitals, this%tempElec, this%Ef, this%spinW,&
            & this%thirdOrd, this%dftbU, this%iEqBlockDftbu, this%onSiteElements,&
            & this%iEqBlockOnSite, this%rangeSep, this%nNeighbourLC, this%pChrgMixer, this%kPoint,&
            & this%kWeight, this%iCellVec, this%cellVec, this%polarisability, this%dEidE,&
            & this%dqOut, this%neFermi, this%dEfdE, errStatus, this%dynRespEFreq)
        if (errStatus%hasError()) then
          call error(errStatus%message)
        end if
        if (this%tWriteBandDat) then
          call writeDerivBandOut(derivEBandOut, this%dEidE, this%kWeight)
        end if
        if (env%tGlobalLead .and. this%tWriteDetailedOut) then
          call writeDetailedOut8(this%fdDetailedOut%unit, this%neFermi)
          call writeDetailedOut9(this%fdDetailedOut%unit, this%orb, this%polarisability,&
              & this%dqOut, this%dEfdE)
        end if
      end if

      if (this%isKernelResp) then
        call this%response%wrtVAtom(env, this%parallelKS, this%tWriteAutotest, autotestTag,&
            & this%tWriteResultsTag, resultsTag, this%taggedWriter, this%tWriteBandDat,&
            & this%fdDetailedOut, this%filling, this%eigen, this%eigVecsReal, this%eigvecsCplx,&
            & this%ints%hamiltonian, this%ints%overlap, this%orb, this%nAtom, this%species,&
            & this%neighbourList, this%nNeighbourSK, this%denseDesc, this%iSparseStart,&
            & this%img2CentCell, this%isRespKernelRPA, this%scc, this%maxSccIter, this%sccTol,&
            & this%isSccConvRequired, this%nMixElements, this%nIneqOrb, this%iEqOrbitals,&
            & this%tempElec, this%Ef, this%spinW, this%thirdOrd, this%dftbU, this%iEqBlockDftbu,&
            & this%onSiteElements, this%iEqBlockOnSite, this%rangeSep, this%nNeighbourLC,&
            & this%pChrgMixer, this%kPoint, this%kWeight, this%iCellVec, this%cellVec,&
            & this%neFermi, errStatus, this%dynKernelFreq, this%tHelical, this%coord)
        if (errStatus%hasError()) then
          call error(errStatus%message)
        end if
        if (env%tGlobalLead .and. this%tWriteDetailedOut) then
          call writeDetailedOut8(this%fdDetailedOut%unit, this%neFermi)
        end if
      end if

    end if

    if (env%tGlobalLead .and. this%tWriteDetailedOut) then
      deallocate(this%fdDetailedOut)
    end if

    if (allocated(this%pipekMezey)) then
      ! NOTE: the canonical DFTB ground state orbitals are over-written after this point
      if (withMpi) then
        call error("Pipek-Mezey localisation does not yet work with MPI")
      end if
      if (this%nSpin > 2) then
        call error("Pipek-Mezey localisation not implemented for non-colinear DFTB")
      end if
      call calcPipekMezeyLocalisation(env, this%pipekMezey, this%tPrintEigvecsTxt, this%nEl,&
          & this%filling, this%ints, this%kPoint, this%neighbourList, this%nNeighbourSk,&
          & this%denseDesc, this%iSparseStart, this%img2CentCell, this%iCellVec, this%cellVec,&
          & this%runId, this%orb, this%species, this%speciesName, this%parallelKS, localisation,&
          & this%eigvecsReal, this%SSqrReal, this%eigvecsCplx, this%SSqrCplx, this%tHelical,&
          & this%coord)
    end if

    if (this%tWriteAutotest.and..not.this%tRestartNoSC) then
      if (this%tPeriodic) then
        this%cellVol = abs(determinant33(this%latVec))
        this%dftbEnergy(this%deltaDftb%iFinal)%EGibbs =&
            & this%dftbEnergy(this%deltaDftb%iFinal)%EMermin + this%extPressure * this%cellVol
      end if
      call writeAutotestTag(autotestTag, this%electronicSolver, this%tPeriodic, this%cellVol,&
          & this%tMulliken, this%qOutput, this%derivs, this%chrgForces, this%excitedDerivs,&
          & this%tStress, this%totalStress, this%pDynMatrix,&
          & this%dftbEnergy(this%deltaDftb%iFinal), this%extPressure, this%coord0, this%tLocalise,&
          & localisation, this%electrostatPot, this%taggedWriter, this%tunneling, this%ldos,&
          & this%lCurrArray, this%polarisability, this%dEidE, this%dipoleMoment, this%eFieldScaling)
    end if
    if (this%tWriteResultsTag) then
      call writeResultsTag(resultsTag, this%dftbEnergy(this%deltaDftb%iFinal), this%derivs,&
          & this%chrgForces, this%nEl, this%Ef, this%eigen, this%filling, this%electronicSolver,&
          & this%tStress, this%totalStress, this%pDynMatrix, this%tPeriodic, this%cellVol,&
          & this%tMulliken, this%qOutput, this%q0, this%taggedWriter, this%cm5Cont,&
          & this%polarisability, this%dEidE, this%dqOut, this%neFermi, this%dEfdE,&
          & this%coord0, this%dipoleMoment, this%multipoleOut, this%eFieldScaling)
    end if
    if (this%tWriteCosmoFile .and. allocated(this%solvation)) then
      call writeCosmoFile(this%solvation, this%species0, this%speciesName, this%coord0, &
          & this%dftbEnergy(this%deltaDftb%iFinal)%EMermin)
    end if
    if (this%tWriteDetailedXML) then
      call writeDetailedXml(this%runId, this%speciesName, this%species0, this%pCoord0Out,&
          & this%tPeriodic, this%tHelical, this%latVec, this%origin, this%tRealHS, this%nKPoint,&
          & this%nSpin, size(this%eigen, dim=1), this%nOrb, this%kPoint, this%kWeight,&
          & this%filling, this%occNatural)
    end if

    call env%globalTimer%stopTimer(globalTimers%postGeoOpt)

  #:if WITH_TRANSPORT
    if (this%electronicSolver%iSolver == electronicSolverTypes%GF .or. &
      & this%electronicSolver%iSolver == electronicSolverTypes%OnlyTransport) then
      call TNegfInt_final(this%negfInt)
    end if
  #:endif

  call this%neighbourList%finalize()

  end subroutine runDftbPlus


  !> Set up charges before determinant calculations
  subroutine preDetCharges(isUnReduced, iDet, nDets, iGeoStep, deltaDftb, qInput, qDets, qBlockIn,&
      & qBlockDets, deltaRhoIn, deltaRhoDets)

    !> Are charge reductions required after this routine
    logical, intent(out) :: isUnReduced

    !> The current determinant being processed
    integer, intent(in) :: iDet

    !> Total number of determinants
    integer, intent(in) :: nDets

    !> Geometry step
    integer, intent(in) :: iGeoStep

    !> Determinant derived type
    type(TDftbDeterminants), intent(in) :: deltaDftb

    !> input charges
    real(dp), intent(inout) :: qInput(:,:,:)

    !> input charges from multiple determinants
    real(dp), intent(inout), allocatable :: qDets(:,:,:,:)

    !> block charge input (if needed for orbital potentials)
    real(dp), intent(inout), allocatable :: qBlockIn(:,:,:,:)

    !> block charge input (if needed for orbital potentials), from multiple determinants
    real(dp), intent(inout), allocatable :: qBlockDets(:,:,:,:,:)

    !> delta density matrix as input for next SCC cycle (if needed for range sep. potentials)
    real(dp), intent(inout), allocatable :: deltaRhoIn(:)

    !> delta density matrix (if needed for range sep. potentials), from multiple determinants
    real(dp), intent(inout), allocatable :: deltaRhoDets(:,:)

    isUnReduced = .false.
    if (nDets > 1) then
      if (iGeoStep == 0) then
        if (deltaDftb%iGround > 0 .and. iDet /= deltaDftb%iGround) then
          qInput(:,:,:) = qDets(:,:,:,deltaDftb%iGround)
          if (allocated(qBlockIn)) then
            qBlockIn(:,:,:,:) = qBlockDets(:,:,:,:,deltaDftb%iGround)
          end if
          if (allocated(deltaRhoIn)) then
            deltaRhoIn(:) = deltaRhoDets(:,deltaDftb%iGround)
          end if
          isUnReduced = .true.
        end if
      else
        qInput(:,:,:) = qDets(:,:,:,iDet)
        if (allocated(qBlockIn)) then
          qBlockIn(:,:,:,:) = qBlockDets(:,:,:,:,iDet)
        end if
        if (allocated(deltaRhoIn)) then
          deltaRhoIn(:) = deltaRhoDets(:,iDet)
        end if
        isUnReduced = .true.
      end if
    end if

  end subroutine preDetCharges


  !> Store (if necessary) charges after determinant calculations
  subroutine postDetCharges(iDet, nDets, qOutput, qDets, qBlockDets, qBlockOut, deltaRhoDets,&
      & deltaRhoOut)

    !> The current determinant being processed
    integer, intent(in) :: iDet

    !> Total number of determinants
    integer, intent(in) :: nDets

    !> output charges
    real(dp), intent(inout) :: qOutput(:,:,:)

    !> output charges from multiple determinants
    real(dp), intent(inout), allocatable :: qDets(:,:,:,:)

    !> block charge output (if needed for orbital potentials), from multiple determinants
    real(dp), intent(inout), allocatable :: qBlockDets(:,:,:,:,:)

    !> block charge output (if needed for orbital potentials)
    real(dp), intent(inout), allocatable :: qBlockOut(:,:,:,:)

    !> delta density matrix (if needed for range sep. potentials), from multiple determinants
    real(dp), intent(inout), allocatable :: deltaRhoDets(:,:)

    !> delta density matrix as input for next SCC cycle (if needed for range sep. potentials)
    real(dp), intent(inout), allocatable :: deltaRhoOut(:)

    if (nDets > 1) then
      qDets(:,:,:,iDet) = qOutput(:,:,:)
      if (allocated(qBlockOut)) then
        qBlockDets(:,:,:,:,iDet) = qBlockOut
      end if
      if (allocated(deltaRhoDets)) then
        deltaRhoDets(:,iDet) = deltaRhoOut
      end if
    end if

  end subroutine postDetCharges


  !> Process current geometry
  subroutine processGeometry(this, env, iGeoStep, iLatGeoStep, tWriteRestart, tStopScc,&
      & tExitGeoOpt, errStatus)

    !> Global variables
    type(TDftbPlusMain), intent(inout) :: this

    !> Environment settings
    type(TEnvironment), intent(inout) :: env

    !> Current geometry step
    integer, intent(in) :: iGeoStep

    !> Current lattice step
    integer, intent(in) :: iLatGeoStep

    !> flag to write out geometries (and charge data if scc)
    logical, intent(in) :: tWriteRestart

    !> if scc driver should be stopped
    logical, intent(out) :: tStopScc

    !> Whether main code should exit the geometry optimisation loop
    logical, intent(out) :: tExitGeoOpt

    !> Status of operation
    type(TStatus), intent(out) :: errStatus

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
    integer :: iSpin, iKS

    real(dp), allocatable :: dipoleTmp(:)

    if (this%tDipole) then
      allocate(dipoleTmp(3))
    end if

    call env%globalTimer%startTimer(globalTimers%preSccInit)

    if (allocated(this%qDepExtPot)) then
      allocate(dQ(this%orb%mShell, this%nAtom, this%nSpin))
    end if

    call this%electronicSolver%reset()
    tExitGeoOpt = .false.

    if (this%tMD .and. tWriteRestart) then
      if (iGeoStep == 0) then
        call openOutputFile(mdOut, .false., this%fdMd)
      end if
      call writeMdOut1(this%fdMd%unit, iGeoStep, this%pMDIntegrator)
    end if

    if (this%tLatticeChanged) then
      call handleLatticeChange(this%latVec, this%scc, this%tblite, this%tStress, this%extPressure,&
          & this%cutOff%mCutOff, this%repulsive, this%dispersion, this%solvation, this%cm5Cont,&
          & this%recVec, this%invLatVec, this%cellVol, this%recCellVol, this%extLatDerivs,&
          & this%cellVec, this%rCellVec, this%boundaryCond)
    end if

    if (this%tCoordsChanged) then
      call handleCoordinateChange(env, this%boundaryCond, this%coord0, this%latVec, this%invLatVec,&
          & this%species0, this%cutOff, this%orb, this%tPeriodic, this%tHelical, this%scc,&
          & this%tblite, this%repulsive, this%dispersion,this%solvation, this%thirdOrd,&
          & this%rangeSep, this%reks, this%img2CentCell, this%iCellVec, this%neighbourList,&
          & this%nAllAtom, this%coord0Fold, this%coord,this%species, this%rCellVec,&
          & this%nNeighbourSk, this%nNeighbourLC, this%ints, this%H0, this%rhoPrim, this%iRhoPrim,&
          & this%ERhoPrim, this%iSparseStart, this%cm5Cont, errStatus)
        @:PROPAGATE_ERROR(errStatus)
    end if

    #:if WITH_TRANSPORT
      if (this%tNegf) then
        call setupNegfStuff(this%negfInt, this%denseDesc, this%transpar, this%ginfo,&
            & this%neighbourList, this%nNeighbourSK, this%img2CentCell, this%orb)
      end if
    #:endif

    if (this%tSccCalc .and. .not. allocated(this%reks) .and. .not.&
        & this%tRestartNoSC) then
      call reset(this%pChrgMixer, this%nMixElements)
    end if

    if (this%electronicSolver%isElsiSolver .and. .not. this%tLargeDenseMatrices) then
      call this%electronicSolver%elsi%updateGeometry(env, this%neighbourList, this%nNeighbourSK,&
          & this%denseDesc%iAtomStart, this%iSparseStart, this%img2CentCell)
    end if

    call env%globalTimer%startTimer(globalTimers%sparseH0S)
    select case(this%hamiltonianType)
    case default
      call error("Invalid Hamiltonian")
    case(hamiltonianTypes%dftb)
      call buildH0(env, this%H0, this%skHamCont, this%atomEigVal, this%coord, this%nNeighbourSk,&
          & this%neighbourList%iNeighbour, this%species, this%iSparseStart, this%orb)
      call buildS(env, this%ints%overlap, this%skOverCont, this%coord, this%nNeighbourSk,&
          & this%neighbourList%iNeighbour, this%species, this%iSparseStart, this%orb)
    case(hamiltonianTypes%xtb)
      @:ASSERT(allocated(this%tblite))
      call this%tblite%buildSH0(env, this%species, this%coord, this%nNeighbourSk, &
          & this%neighbourList%iNeighbour, this%img2CentCell, this%iSparseStart, &
          & this%orb, this%H0, this%ints%overlap, this%ints%dipoleBra, this%ints%dipoleKet, &
          & this%ints%quadrupoleBra, this%ints%quadrupoleKet)
    end select
    call env%globalTimer%stopTimer(globalTimers%sparseH0S)

    if (this%tSetFillingTemp) then
      call this%temperatureProfile%getTemperature(this%tempElec)
    end if

    call this%electronicSolver%updateElectronicTemp(this%tempElec)

    if (allocated(this%repulsive)) then
      call this%repulsive%getEnergy(this%coord, this%species, this%img2CentCell,&
          & this%neighbourList,this%dftbEnergy(this%deltaDftb%iDeterminant)%atomRep,&
          & this%dftbEnergy(this%deltaDftb%iDeterminant)%ERep,&
          & iAtInCentralRegion=this%iAtInCentralRegion)
    end if

    if (allocated(this%halogenXCorrection)) then
      call this%halogenXCorrection%getEnergies(this%dftbEnergy(&
          & this%deltaDftb%iDeterminant)%atomHalogenX, this%coord, this%species,&
          & this%neighbourList, this%img2CentCell)
      this%dftbEnergy(this%deltaDftb%iDeterminant)%EHalogenX =&
          & sum(this%dftbEnergy(this%deltaDftb%iDeterminant)%atomHalogenX(this%iAtInCentralRegion))
    end if

    call resetExternalPotentials(this%refExtPot, this%potential)

    if (this%tReadShifts) then
      call readShifts(fShifts, this%orb, this%nAtom, this%nSpin, this%potential%extShell)
    end if

    call addUpExternalField(this%eField, this%tPeriodic, this%neighbourList, this%nNeighbourSk,&
        & this%iCellVec, this%img2CentCell, this%cellVec, this%deltaT, iGeoStep, this%coord0Fold,&
        & this%coord, this%potential)

    call mergeExternalPotentials(this%orb, this%species, this%potential)

    ! For non-scc calculations with transport only, jump out of geometry loop
    if (this%electronicSolver%iSolver == electronicSolverTypes%OnlyTransport) then
      if (this%tWriteDetailedOut) then
        call openOutputFile(userOut, tAppendDetailedOut, this%fdDetailedOut)
      end if
      ! We need to define hamiltonian by adding the potential
      call getSccHamiltonian(this%H0, this%ints, this%nNeighbourSK, this%neighbourList,&
          & this%species, this%orb, this%iSparseStart, this%img2CentCell, this%potential,&
          & allocated(this%reks), this%ints%hamiltonian, this%ints%iHamiltonian)
      tExitGeoOpt = .true.
      return
    end if

    if (this%electronicSolver%iSolver == electronicSolverTypes%pexsi) then
      call this%electronicSolver%elsi%initPexsiDeltaVRanges(this%tSccCalc, this%potential)
    end if

    if (.not.this%tRestartNoSC) then
      call initSccLoop(this%tSccCalc, this%xlbomdIntegrator, this%minSccIter, this%maxSccIter,&
          & this%sccTol, tConverged, this%tNegf, this%reks)
    else
      tConverged = .true.
    end if

    call env%globalTimer%stopTimer(globalTimers%preSccInit)

    call env%globalTimer%startTimer(globalTimers%scc)

    REKS_SCC: if (allocated(this%reks)) then

      lpSCC_REKS: do iSccIter = 1, this%maxSccIter

        if (iSccIter == 1) then
          call getReksInitialSettings(env, this%denseDesc, this%h0, this%ints, this%neighbourList,&
              & this%nNeighbourSK, this%iSparseStart, this%img2CentCell, this%electronicSolver,&
              & iGeoStep, this%HSqrReal, this%SSqrReal, this%eigvecsReal, this%eigen, this%reks,&
              & errStatus)
          @:PROPAGATE_ERROR(errStatus)
        end if

        call getDensityMatrixL(env, this%denseDesc, this%neighbourList, this%nNeighbourSK,&
            & this%iSparseStart, this%img2CentCell, this%orb, this%species, this%coord,&
            & this%tHelical, this%eigvecsReal, this%parallelKS, this%rhoPrim, this%SSqrReal,&
            & this%rhoSqrReal, this%q0, this%deltaRhoOutSqr, this%reks)
        call getMullikenPopulationL(env, this%denseDesc, this%neighbourList, this%nNeighbourSK,&
            & this%img2CentCell, this%iSparseStart, this%orb, this%rhoPrim, this%ints,&
            & this%iRhoPrim, this%qBlockOut, this%qiBlockOut, this%qNetAtom, this%reks)

        call getHamiltonianLandEnergyL(env, this%denseDesc, this%scc, this%tblite, this%orb,&
            & this%species, this%neighbourList, this%nNeighbourSK, this%iSparseStart,&
            & this%img2CentCell, this%H0, this%ints, this%spinW, this%cellVol, this%extPressure,&
            & this%dftbEnergy(1), this%q0, this%iAtInCentralRegion, this%solvation, this%thirdOrd,&
            & this%potential, this%rangeSep, this%nNeighbourLC, this%tDualSpinOrbit, this%xi,&
            & this%isExtField, this%isXlbomd, this%dftbU, this%dftbEnergy(1)%TS, this%qDepExtPot,&
            & this%qBlockOut, this%qiBlockOut, this%tFixEf, this%Ef, this%rhoPrim,&
            & this%onSiteElements, this%dispersion, tConverged, this%species0, this%referenceN0,&
            & this%qNetAtom, this%multipoleOut, this%reks)
        call optimizeFONsAndWeights(this%eigvecsReal, this%filling, this%dftbEnergy(1), this%reks)

        call getFockandDiag(env, this%denseDesc, this%neighbourList, this%nNeighbourSK,&
            & this%iSparseStart, this%img2CentCell, this%eigvecsReal, this%electronicSolver,&
            & this%eigen, this%reks)

        ! Creates (delta) density matrix for averaged state from real eigenvectors.
        call getDensityFromRealEigvecs(env, this%denseDesc, this%filling(:,1,:),&
            & this%neighbourList, this%nNeighbourSK, this%iSparseStart, this%img2CentCell,&
            & this%orb, this%species, this%denseDesc%iAtomStart, this%coord, this%tHelical,&
            & this%eigVecsReal, this%parallelKS, this%rhoPrim, this%SSqrReal, this%rhoSqrReal,&
            & this%deltaRhoOutSqr)
        ! For rangeseparated calculations deduct atomic charges from deltaRho
        if (this%isRangeSep) then
          call denseSubtractDensityOfAtoms(this%q0, this%denseDesc%iAtomStart, this%deltaRhoOutSqr)
        end if
        call getMullikenPopulation(this%rhoPrim, this%ints, this%orb, this%neighbourList,&
            & this%nNeighbourSK, this%img2CentCell, this%iSparseStart, this%qOutput,&
            & iRhoPrim=this%iRhoPrim, qBlock=this%qBlockOut, qiBlock=this%qiBlockOut,&
            & qNetAtom=this%qNetAtom, multipoles=this%multipoleOut)

        ! Check charge convergence and guess new eigenvectors
        tStopScc = hasStopFile(fStopScc)
        if (this%isRangeSep) then
          call getReksNextInputDensity(sccErrorQ, this%sccTol, tConverged, iSccIter,&
              & this%minSccIter, this%maxSccIter, iGeoStep, tStopScc, this%eigvecsReal,&
              & this%deltaRhoOut, this%deltaRhoIn, this%deltaRhoDiff, this%reks)
        else
          call getReksNextInputCharges(this%qInput, this%qOutput, this%qDiff, sccErrorQ,&
              & this%sccTol, tConverged, iSccIter, this%minSccIter, this%maxSccIter, iGeoStep,&
              & tStopScc, this%eigvecsReal, this%reks)
        end if

        call getSccInfo(iSccIter, this%dftbEnergy(1)%Eavg, Eold, diffElec)
        call printReksSccInfo(iSccIter, this%dftbEnergy(1)%Eavg, diffElec, sccErrorQ,&
            & this%reks)

        if (tConverged .or. tStopScc) then

          call printReksSAInfo(this%reks, this%dftbEnergy(1)%Etotal)

          call getStateInteraction(env, this%denseDesc, this%neighbourList, this%nNeighbourSK,&
              & this%iSparseStart, this%img2CentCell, this%coord, this%iAtInCentralRegion,&
              & this%eigvecsReal, this%electronicSolver, this%eigen, this%qOutput, this%q0,&
              & this%tDipole, dipoleTmp, this%reks)
          call assignDipoleMoment(dipoleTmp, this%dipoleMoment, this%deltaDftb%iDeterminant,&
              & this%tDipole, this%reks, isSingleState=.true.)

          call getReksEnProperties(env, this%denseDesc, this%neighbourList, this%nNeighbourSK,&
              & this%img2CentCell, this%iSparseStart, this%eigvecsReal, this%coord0, this%reks)

          if (this%tWriteDetailedOut .and. this%deltaDftb%nDeterminant() == 1) then
            ! In this routine the correct Etotal is evaluated.
            ! If TargetStateL > 0, certain microstate
            ! is optimized. If not, SSR state is optimized.
            call openOutputFile(userOut, tAppendDetailedOut, this%fdDetailedOut)
            call writeReksDetailedOut1(this%fdDetailedOut%unit, this%nGeoSteps, iGeoStep,&
                & this%tMD, this%tDerivs, this%tCoordOpt, this%tLatOpt, iLatGeoStep, iSccIter,&
                & this%dftbEnergy(1), diffElec, sccErrorQ, this%indMovedAtom, this%pCoord0Out,&
                & this%q0, this%qOutput, this%orb, this%species, this%tPrintMulliken,&
                & this%extPressure, this%cellVol, this%dftbEnergy(1)%TS, this%tAtomicEnergy,&
                & this%dispersion, this%tPeriodic, this%tSccCalc, this%invLatVec, this%kPoint,&
                & this%iAtInCentralRegion, this%electronicSolver, this%reks,&
                & allocated(this%thirdOrd), this%isRangeSep, qNetAtom=this%qNetAtom)
          end if
          if (this%tWriteBandDat) then
            call writeBandOut(bandOut, this%eigen, this%filling, this%kWeight)
          end if

          exit lpSCC_REKS
        end if
      end do lpSCC_REKS

    else ! not REKS_SCC

      ! Standard spin free or unrestricted DFTB

      lpSCC: do iSccIter = 1, this%maxSccIter

        call resetInternalPotentials(this%tDualSpinOrbit, this%xi, this%orb, this%species,&
            & this%potential)

        if (this%tSccCalc) then

          call getChargePerShell(this%qInput, this%orb, this%species, this%chargePerShell)

        #:if WITH_TRANSPORT
          ! Overrides input charges with uploaded contact charges
          if (this%tUpload) then
            call overrideContactCharges(this%qInput, this%chargeUp, this%transpar, this%qBlockIn,&
                & this%blockUp)
          end if
        #:endif

          call addChargePotentials(env, this%scc, this%tblite, .true., this%qInput, this%q0,&
              & this%chargePerShell, this%orb, this%multipoleInp, this%species, this%neighbourList,&
              & this%img2CentCell, this%spinW, this%solvation, this%thirdOrd, this%dispersion,&
              & this%potential)

          call addBlockChargePotentials(this%qBlockIn, this%qiBlockIn, this%dftbU, this%tImHam,&
              & this%species, this%orb, this%potential)

          if (allocated(this%onSiteElements) .and. (iSCCIter > 1 .or. this%tReadChrg)) then
            call addOnsShift(this%potential%intBlock, this%potential%iOrbitalBlock, this%qBlockIn,&
                & this%qiBlockIn, this%onSiteElements, this%species, this%orb, this%q0)
          end if

        end if

        ! All potentials are added up into intBlock
        this%potential%intBlock = this%potential%intBlock + this%potential%extBlock

        if (allocated(this%qDepExtPot)) then
          call getChargePerShell(this%qInput, this%orb, this%species, dQ, qRef=this%q0)
          call this%qDepExtPot%addPotential(sum(dQ(:,:,1), dim=1), dQ(:,:,1), this%orb,&
              & this%species, this%potential%intBlock)
        end if

        if (this%electronicSolver%iSolver == electronicSolverTypes%pexsi .and. this%tSccCalc) then
          call this%electronicSolver%elsi%updatePexsiDeltaVRanges(this%potential)
        end if

        call getSccHamiltonian(this%H0, this%ints, this%nNeighbourSK, this%neighbourList,&
            & this%species, this%orb, this%iSparseStart, this%img2CentCell, this%potential,&
            & allocated(this%reks), this%ints%hamiltonian, this%ints%iHamiltonian)

        if (this%tWriteRealHS .or. this%tWriteHS&
            & .and. any(this%electronicSolver%iSolver&
            & == [electronicSolverTypes%qr, electronicSolverTypes%divideandconquer,&
            & electronicSolverTypes%relativelyrobust, electronicSolverTypes%magma_gvd])) then
          call writeHSAndStop(env, this%tWriteHS, this%tWriteRealHS, this%tRealHS,&
              & this%ints%overlap, this%neighbourList, this%nNeighbourSK,&
              & this%denseDesc%iAtomStart, this%iSparseStart, this%img2CentCell, this%kPoint,&
              & this%iCellVec, this%cellVec, this%ints%hamiltonian, this%ints%iHamiltonian)
        end if

        call convertToUpDownRepr(this%ints%hamiltonian, this%ints%iHamiltonian)

        call getDensity(env, this%negfInt, iSccIter, this%denseDesc, this%ints,&
            & this%neighbourList, this%nNeighbourSk, this%iSparseStart, this%img2CentCell,&
            & this%iCellVec, this%cellVec, this%kPoint, this%kWeight, this%orb, this%tHelical,&
            & this%coord, this%species, this%electronicSolver, this%tRealHS, this%tSpinSharedEf,&
            & this%tSpinOrbit, this%tDualSpinOrbit, this%tFillKSep, this%tFixEf, this%tMulliken,&
            & this%iDistribFn, this%tempElec, this%nEl, this%parallelKS, this%Ef, this%mu,&
            & this%dftbEnergy(this%deltaDftb%iDeterminant), this%rangeSep, this%eigen,&
            & this%filling, this%rhoPrim, this%xi, this%orbitalL, this%HSqrReal,&
            & this%SSqrReal, this%eigvecsReal, this%iRhoPrim, this%HSqrCplx, this%SSqrCplx,&
            & this%eigvecsCplx, this%rhoSqrReal, this%deltaRhoInSqr, this%deltaRhoOutSqr,&
            & this%qOutput, this%nNeighbourLC, this%tLargeDenseMatrices, this%deltaDftb, errStatus)
        if (errStatus%hasError()) then
          call error(errStatus%message)
        end if

        !> For rangeseparated calculations deduct atomic charges from deltaRho
        if (this%isRangeSep) then
          select case(this%nSpin)
          case(2)
            do iSpin = 1, 2
              call denseSubtractDensityOfAtoms(this%q0, this%denseDesc%iAtomStart,&
                  & this%deltaRhoOutSqr, iSpin)
            end do
          case(1)
            call denseSubtractDensityOfAtoms(this%q0, this%denseDesc%iAtomStart,&
                & this%deltaRhoOutSqr)
          case default
            call error("Range separation not implemented for non-colinear spin")
          end select
        end if

        if (this%tWriteBandDat .and. this%deltaDftb%nDeterminant() == 1) then
          call writeBandOut(bandOut, this%eigen, this%filling, this%kWeight)
        end if

        if (this%tMulliken) then
          call getMullikenPopulation(this%rhoPrim, this%ints, this%orb, this%neighbourList,&
              & this%nNeighbourSk, this%img2CentCell, this%iSparseStart, this%qOutput,&
              & iRhoPrim=this%iRhoPrim, qBlock=this%qBlockOut, qiBlock=this%qiBlockOut,&
              & qNetAtom=this%qNetAtom, multipoles=this%multipoleOut)

          if (this%tSpinSharedEf .or. this%tFixEf .or.&
              & this%electronicSolver%iSolver == electronicSolverTypes%GF) then
            this%nEl(:) = sum(sum(this%qOutput(:, this%iAtInCentralRegion, :size(this%nEl)),dim=1),&
                & dim=1)
            call qm2ud(this%nEl)
          end if

        end if

      #:if WITH_TRANSPORT
        ! Override charges with uploaded contact charges
        if (this%tUpload) then
          call overrideContactCharges(this%qOutput, this%chargeUp, this%transpar, this%qBlockIn,&
              & this%blockUp)
        end if
      #:endif

        ! For non-dual spin-orbit orbitalL is determined during getDensity() call above
        if (this%tDualSpinOrbit) then
          call getLDual(this%orbitalL, this%qiBlockOut, this%orb, this%species)
        end if

        ! Note: if XLBOMD is active, potential created with input charges is needed later,
        ! therefore it should not be overwritten here.
        if (this%tSccCalc .and. .not. this%isXlbomd) then
          call resetInternalPotentials(this%tDualSpinOrbit, this%xi, this%orb, this%species,&
              & this%potential)
          call getChargePerShell(this%qOutput, this%orb, this%species, this%chargePerShell)

          call addChargePotentials(env, this%scc, this%tblite, this%updateSccAfterDiag,&
              & this%qOutput, this%q0, this%chargePerShell, this%orb, this%multipoleOut,&
              & this%species, this%neighbourList, this%img2CentCell, this%spinW, this%solvation,&
              & this%thirdOrd, this%dispersion, this%potential)

          call addBlockChargePotentials(this%qBlockOut, this%qiBlockOut, this%dftbU, this%tImHam,&
              & this%species, this%orb, this%potential)

          if (allocated(this%onSiteElements)) then
            call addOnsShift(this%potential%intBlock, this%potential%iOrbitalBlock, this%qBlockOut,&
                & this%qiBlockOut, this%onSiteElements, this%species, this%orb, this%q0)
          end if

          this%potential%intBlock = this%potential%intBlock + this%potential%extBlock
        end if

        if (allocated(this%qDepExtPot)) then
          call getChargePerShell(this%qOutput, this%orb, this%species, dQ, qRef=this%q0)
          call this%qDepExtPot%addPotential(sum(dQ(:,:,1), dim=1), dQ(:,:,1), this%orb,&
              & this%species, this%potential%intBlock)
        end if

        call calcEnergies(this%scc, this%tblite, this%qOutput, this%q0, this%chargePerShell,&
            & this%multipoleOut, this%species, this%isExtField, this%isXlbomd, this%dftbU,&
            & this%tDualSpinOrbit, this%rhoPrim, this%H0, this%orb, this%neighbourList,&
            & this%nNeighbourSk, this%img2CentCell, this%iSparseStart, this%cellVol,&
            & this%extPressure, this%dftbEnergy(this%deltaDftb%iDeterminant)%TS, this%potential,&
            & this%dftbEnergy(this%deltaDftb%iDeterminant), this%thirdOrd, this%solvation,&
            & this%rangeSep, this%reks, this%qDepExtPot, this%qBlockOut, this%qiBlockOut,&
            & this%xi, this%iAtInCentralRegion, this%tFixEf, this%Ef, this%onSiteElements,&
            & this%qNetAtom, this%potential%intOnSiteAtom, this%potential%extOnSiteAtom)

        tStopScc = hasStopFile(fStopScc)

        ! Mix charges Input/Output
        if (this%tSccCalc) then
          if(.not. this%isRangeSep) then
            call getNextInputCharges(env, this%pChrgMixer, this%qOutput, this%qOutRed, this%orb,&
                & this%nIneqOrb, this%iEqOrbitals, iGeoStep, iSccIter, this%minSccIter,&
                & this%maxSccIter, this%sccTol, tStopScc, this%tMixBlockCharges, this%tReadChrg,&
                & this%qInput, this%qInpRed, sccErrorQ, tConverged, this%dftbU, this%qBlockOut,&
                & this%iEqBlockDftbU, this%qBlockIn, this%qiBlockOut, this%iEqBlockDftbULS,&
                & this%species0, this%qiBlockIn, this%iEqBlockOnSite, this%iEqBlockOnSiteLS,&
                & this%nIneqDip, this%nIneqQuad, this%iEqDipole, this%iEqQuadrupole, &
                & this%multipoleOut, this%multipoleInp)
          else
            call getNextInputDensity(this%SSqrReal, this%ints, this%neighbourList,&
                & this%nNeighbourSK, this%denseDesc%iAtomStart, this%iSparseStart,&
                & this%img2CentCell, this%pChrgMixer, this%qOutput, this%orb, this%tHelical,&
                & this%species, this%coord, iGeoStep, iSccIter, this%minSccIter, this%maxSccIter,&
                & this%sccTol, tStopScc, this%tReadChrg, this%q0, this%qInput, sccErrorQ,&
                & tConverged, this%deltaRhoOut, this%deltaRhoIn, this%deltaRhoDiff, this%qBlockIn,&
                & this%qBlockOut)
          end if

          call getSccInfo(iSccIter, this%dftbEnergy(this%deltaDftb%iDeterminant)%Eelec, Eold,&
              & diffElec)
          if (this%tNegf) then
            call printSccHeader()
          end if
          call printSccInfo(allocated(this%dftbU), iSccIter,&
              & this%dftbEnergy(this%deltaDftb%iDeterminant)%Eelec, diffElec, sccErrorQ)

          if (this%tNegf) then
            call printBlankLine()
          end if

          tWriteSccRestart = env%tGlobalLead .and. needsSccRestartWriting(this%restartFreq,&
              & iGeoStep, iSccIter, this%minSccIter, this%maxSccIter, this%tMd, &
              & this%isGeoOpt .or. allocated(this%geoOpt),&
              & this%tDerivs, tConverged, this%tReadChrg, tStopScc)
          if (tWriteSccRestart) then
            call writeCharges(fCharges, this%tWriteChrgAscii, this%orb, this%qInput, this%qBlockIn,&
                & this%qiBlockIn, this%deltaRhoIn, size(this%iAtInCentralRegion), this%multipoleInp)
          end if
        end if

        if (allocated(this%dispersion) .and. .not. tConverged) then
          call this%dispersion%updateOnsiteCharges(this%qNetAtom, this%orb, this%referenceN0,&
              & this%species0, tConverged)
          call calcDispersionEnergy(this%dispersion,&
              & this%dftbEnergy(this%deltaDftb%iDeterminant)%atomDisp,&
              & this%dftbEnergy(this%deltaDftb%iDeterminant)%Edisp, this%iAtInCentralRegion)
        end if
        call sumEnergies(this%dftbEnergy(this%deltaDftb%iDeterminant))

        if (this%tWriteDetailedOut .and. this%deltaDftb%nDeterminant() == 1) then
          call openOutputFile(userOut, tAppendDetailedOut, this%fdDetailedOut)
          call writeDetailedOut1(this%fdDetailedOut%unit, this%iDistribFn, this%nGeoSteps,&
              & iGeoStep, this%tMD, this%tDerivs, this%tCoordOpt, this%tLatOpt, iLatGeoStep,&
              & iSccIter, this%dftbEnergy(this%deltaDftb%iDeterminant), diffElec, sccErrorQ,&
              & this%indMovedAtom, this%pCoord0Out, this%tPeriodic, this%tSccCalc, this%tNegf,&
              & this%invLatVec, this%kPoint)
          call writeDetailedOut2(this%fdDetailedOut%unit, this%q0, this%qInput, this%qOutput,&
              & this%orb, this%species, allocated(this%dftbU), this%tImHam .or. this%tSpinOrbit,&
              & this%tPrintMulliken, this%orbitalL, this%qBlockOut, this%nSpin,&
              & allocated(this%onSiteElements), this%iAtInCentralRegion, this%cm5Cont,&
              & this%qNetAtom)
          call writeDetailedOut3(this%fdDetailedOut%unit, this%qInput, this%qOutput,&
              & this%dftbEnergy(this%deltaDftb%iDeterminant), this%species, allocated(this%dftbU),&
              & this%tPrintMulliken, this%Ef, this%extPressure, this%cellVol, this%tAtomicEnergy,&
              & this%dispersion, allocated(this%eField), this%tPeriodic, this%nSpin, this%tSpin,&
              & this%tSpinOrbit, this%tSccCalc, allocated(this%onSiteElements), this%tNegf,&
              & this%iAtInCentralRegion, this%electronicSolver, allocated(this%halogenXCorrection),&
              & this%isRangeSep, allocated(this%thirdOrd), allocated(this%solvation))
        end if

        if (tConverged .or. tStopScc) then
          exit lpSCC
        end if

      end do lpSCC

    end if REKS_SCC

  #:if WITH_SCALAPACK
    if (this%isSparseReorderRequired) then
      if (allocated(this%eigvecsReal)) then
        do iKS = 1, this%parallelKS%nLocalKS
          call blacsfx_gemr2d(this%denseDesc%nOrb, this%denseDesc%nOrb,&
              & this%eigvecsReal(:,:,iKS), 1, 1, this%denseDesc%blacsOrbSqr,&
              & this%eigVecsRealReordered(:,:,iKS), 1, 1, this%denseDesc%blacsColumnSqr,&
              & env%blacs%orbitalGrid%ctxt)
        end do
      else if (allocated(this%eigvecsCplx)) then
        do iKS = 1, this%parallelKS%nLocalKS
          call blacsfx_gemr2d(this%denseDesc%nOrb, this%denseDesc%nOrb,&
              & this%eigvecsCplx(:,:,iKS), 1, 1, this%denseDesc%blacsOrbSqr,&
              & this%eigVecsCplxReordered(:,:,iKS), 1, 1, this%denseDesc%blacsColumnSqr,&
              & env%blacs%orbitalGrid%ctxt)
        end do
      end if
    end if
  #:endif

    if (allocated(this%dispersion) .and. .not.allocated(this%reks)) then
      ! If we get to this point for a dispersion model, if it is charge dependent it may require
      ! evaluation post-hoc if SCC was not achieved but the input settings are to proceed with
      ! non-converged SCC.
      call this%dispersion%updateOnsiteCharges(this%qNetAtom, this%orb, this%referenceN0,&
          & this%species0, tConverged .or. .not. this%isSccConvRequired)
      call calcDispersionEnergy(this%dispersion,&
          & this%dftbEnergy(this%deltaDftb%iDeterminant)%atomDisp,&
          & this%dftbEnergy(this%deltaDftb%iDeterminant)%Edisp,&
          & this%iAtInCentralRegion)
      call sumEnergies(this%dftbEnergy(this%deltaDftb%iDeterminant))
    end if

    if (this%tWriteDetailedOut .and. this%deltaDftb%nDeterminant() == 1) then
      deallocate(this%fdDetailedOut)
      call openOutputFile(userOut, tAppendDetailedOut, this%fdDetailedOut)
      if (allocated(this%reks)) then
        call writeReksDetailedOut1(this%fdDetailedOut%unit, this%nGeoSteps, iGeoStep, this%tMD,&
            & this%tDerivs, this%tCoordOpt, this%tLatOpt, iLatGeoStep, iSccIter,&
            & this%dftbEnergy(1), diffElec, sccErrorQ, this%indMovedAtom, this%pCoord0Out,&
            & this%q0, this%qOutput, this%orb, this%species, this%tPrintMulliken,&
            & this%extPressure, this%cellVol, this%dftbEnergy(1)%TS, this%tAtomicEnergy,&
            & this%dispersion, this%tPeriodic, this%tSccCalc, this%invLatVec, this%kPoint,&
            & this%iAtInCentralRegion, this%electronicSolver, this%reks,&
            & allocated(this%thirdOrd), this%isRangeSep, qNetAtom=this%qNetAtom)
      else
        call writeDetailedOut1(this%fdDetailedOut%unit, this%iDistribFn, this%nGeoSteps,&
            & iGeoStep, this%tMD, this%tDerivs, this%tCoordOpt, this%tLatOpt, iLatGeoStep,&
            & iSccIter, this%dftbEnergy(this%deltaDftb%iDeterminant), diffElec, sccErrorQ,&
            & this%indMovedAtom, this%pCoord0Out, this%tPeriodic, this%tSccCalc, this%tNegf,&
            & this%invLatVec, this%kPoint)
        call writeDetailedOut2(this%fdDetailedOut%unit, this%q0, this%qInput, this%qOutput,&
            & this%orb, this%species, allocated(this%dftbU), this%tImHam.or.this%tSpinOrbit,&
            & this%tPrintMulliken, this%orbitalL, this%qBlockOut, this%nSpin,&
            & allocated(this%onSiteElements), this%iAtInCentralRegion, this%cm5Cont, this%qNetAtom)
        call writeDetailedOut3(this%fdDetailedOut%unit, this%qInput, this%qOutput,&
            & this%dftbEnergy(this%deltaDftb%iDeterminant), this%species, allocated(this%dftbU),&
            & this%tPrintMulliken, this%Ef, this%extPressure, this%cellVol, this%tAtomicEnergy,&
            & this%dispersion, allocated(this%eField), this%tPeriodic, this%nSpin, this%tSpin,&
            & this%tSpinOrbit, this%tSccCalc, allocated(this%onSiteElements), this%tNegf,&
            & this%iAtInCentralRegion, this%electronicSolver, allocated(this%halogenXCorrection),&
            & this%isRangeSep, allocated(this%thirdOrd), allocated(this%solvation))
      end if
    end if

    call env%globalTimer%stopTimer(globalTimers%scc)

    if (allocated(this%scc)) then
      call this%scc%finishSccLoop(env)
    end if

    if (allocated(this%dispersion)) then
      if (.not.this%dispersion%energyAvailable()) then
        call warning("Dispersion contributions are not included in the energy")
      end if
    end if

    if (this%tSccCalc .and. .not. this%isXlbomd .and. .not. tConverged&
        & .and. .not. this%tRestartNoSC) then
      call warning("SCC is NOT converged, maximal SCC iterations exceeded")
      if (this%isSccConvRequired) then
        call env%shutdown()
      end if
    end if

    call env%globalTimer%startTimer(globalTimers%postSCC)

    if (this%isLinResp) then
      call calculateLinRespExcitations(env, this%linearResponse, this%parallelKS, this%scc,&
          & this%qOutput, this%q0, this%ints, this%eigvecsReal, this%eigen(:,1,:),&
          & this%filling(:,1,:), this%coord, this%species, this%speciesName, this%orb,&
          & this%skHamCont, this%skOverCont, autotestTag, this%taggedWriter, this%runId,&
          & this%neighbourList, this%nNeighbourSK, this%denseDesc, this%iSparseStart,&
          & this%img2CentCell, this%tWriteAutotest, this%tCasidaForces, this%tLinRespZVect,&
          & this%tPrintExcitedEigvecs, this%tPrintEigvecsTxt, this%nonSccDeriv,&
          & this%dftbEnergy(1), this%energiesCasida, this%SSqrReal, this%rhoSqrReal,&
          & this%deltaRhoOutSqr, this%excitedDerivs, this%dQAtomEx, this%occNatural, this%rangeSep)
    end if

    if (allocated(this%ppRPA)) then
      call unpackHS(this%SSqrReal, this%ints%overlap, this%neighbourList%iNeighbour,&
          & this%nNeighbourSK, this%denseDesc%iAtomStart, this%iSparseStart, this%img2CentCell)
      call blockSymmetrizeHS(this%SSqrReal, this%denseDesc%iAtomStart)
      if (withMpi) then
        call error("pp-RPA calc. does not work with MPI yet")
      end if
      call ppRPAenergies(this%ppRPA, this%denseDesc, this%eigvecsReal, this%eigen(:,1,:),&
          & this%scc, this%SSqrReal, this%species0, this%nEl(1), this%neighbourList%iNeighbour,&
          & this%img2CentCell, this%orb, this%tWriteAutotest, autotestTag, this%taggedWriter)
    end if

    if (this%isXlbomd) then
      call getXlbomdCharges(this%xlbomdIntegrator, this%qOutRed, this%pChrgMixer, this%orb,&
          & this%nIneqOrb, this%iEqOrbitals, this%qInput, this%qInpRed, this%dftbU,&
          & this%iEqBlockDftbU, this%qBlockIn, this%species0, this%iEqBlockDftbuLs, this%qiBlockIn,&
          & this%iEqBlockOnSite, this%iEqBlockOnSiteLS)
    end if

    if (this%tDipole .and. .not. allocated(this%reks) .and. .not. this%tRestartNoSC) then
      call getDipoleMoment(this%qOutput, this%q0, this%multipoleOut%dipoleAtom, this%coord,&
          & this%dipoleMoment(:,this%deltaDftb%iDeterminant), this%iAtInCentralRegion)
    #:block DEBUG_CODE
      if (this%hamiltonianType == hamiltonianTypes%dftb) then
        call checkDipoleViaHellmannFeynman(this%rhoPrim, this%q0, this%coord0, this%ints, this%orb,&
            & this%neighbourList, this%nNeighbourSk, this%species, this%iSparseStart,&
            & this%img2CentCell, this%eFieldScaling)
      end if
    #:endblock DEBUG_CODE
    end if

    call env%globalTimer%startTimer(globalTimers%eigvecWriting)

    if (this%tPrintEigVecs) then
      call writeEigenvectors(env, this%runId, this%neighbourList, this%nNeighbourSk, this%cellVec,&
          & this%iCellVec, this%denseDesc, this%iSparseStart, this%img2CentCell, this%species,&
          & this%speciesName, this%orb, this%kPoint, this%ints%overlap, this%parallelKS,&
          & this%tPrintEigvecsTxt, this%eigvecsReal, this%SSqrReal, this%eigvecsCplx, this%SSqrCplx)
    end if

    if (this%tProjEigenvecs) then
      call writeProjectedEigenvectors(env, this%regionLabels, this%eigen, this%neighbourList,&
          & this%nNeighbourSk, this%cellVec, this%iCellVec, this%denseDesc, this%iSparseStart,&
          & this%img2CentCell, this%orb, this%ints%overlap, this%kPoint, this%kWeight,&
          & this%iOrbRegion, this%parallelKS, this%eigvecsReal, this%SSqrReal, this%eigvecsCplx,&
          & this%SSqrCplx)
    end if
    call env%globalTimer%stopTimer(globalTimers%eigvecWriting)

    ! MD geometry files are written only later, once velocities for the current geometry are known
    if ((this%isGeoOpt .or. allocated(this%geoOpt)) .and. tWriteRestart) then
      if (.not. (this%deltaDftb%isSpinPurify .and.&
          & this%deltaDftb%iDeterminant == determinants%triplet)) then
        call writeCurrentGeometry(this%geoOutFile, this%pCoord0Out, this%tLatOpt, this%tMd,&
            & this%tAppendGeo.and.iGeoStep>0, this%tFracCoord, this%tPeriodic, this%tHelical,&
            & this%tPrintMulliken, this%species0, this%speciesName, this%latVec, this%origin,&
            & iGeoStep, iLatGeoStep, this%nSpin, this%qOutput, this%velocities)
      endif
    end if

    if (this%tForces) then
      call env%globalTimer%startTimer(globalTimers%forceCalc)
      if (allocated(this%reks)) then
        call getReksGradients(env, this%denseDesc, this%scc, this%rangeSep, this%dispersion,&
            & this%neighbourList, this%nNeighbourSK, this%iSparseStart, this%img2CentCell,&
            & this%orb, this%nonSccDeriv, this%skHamCont, this%skOverCont, this%repulsive,&
            & this%coord, this%coord0, this%species, this%q0, this%eigvecsReal,&
            & this%chrgForces, this%ints%overlap, this%spinW, this%derivs, this%tWriteAutotest,&
            & autotestTag, this%taggedWriter, this%reks)
        call getReksGradProperties(env, this%denseDesc, this%neighbourList, this%nNeighbourSK,&
            & this%iSparseStart, this%img2CentCell, this%eigvecsReal, this%orb,&
            & this%iAtInCentralRegion, this%coord, this%coord0, this%ints%overlap, this%rhoPrim,&
            & this%qOutput, this%q0, this%tDipole, dipoleTmp, this%chrgForces, this%reks)
        call assignDipoleMoment(dipoleTmp, this%dipoleMoment, this%deltaDftb%iDeterminant,&
            & this%tDipole, this%reks, isSingleState=.false.)
      else
        call env%globalTimer%startTimer(globalTimers%energyDensityMatrix)
        call getEnergyWeightedDensity(env, this%negfInt, this%electronicSolver, this%denseDesc,&
            & this%forceType, this%filling, this%eigen, this%kPoint, this%kWeight,&
            & this%neighbourList, this%nNeighbourSK, this%orb, this%iSparseStart,&
            & this%img2CentCell, this%iCellVec, this%cellVec, this%tRealHS, this%ints,&
            & this%parallelKS, this%tHelical, this%species, this%coord, iSccIter, this%mu,&
            & this%ERhoPrim, this%eigvecsReal, this%SSqrReal, this%eigvecsCplx, this%SSqrCplx,&
            & errStatus)
        @:PROPAGATE_ERROR(errStatus)
        call env%globalTimer%stopTimer(globalTimers%energyDensityMatrix)
        call getGradients(env, this%boundaryCond, this%scc, this%tblite, this%isExtField,&
            & this%isXlbomd, this%nonSccDeriv, this%rhoPrim, this%ERhoPrim, this%qOutput, this%q0,&
            & this%skHamCont, this%skOverCont, this%repulsive, this%neighbourList,&
            & this%nNeighbourSk, this%species, this%img2CentCell, this%iSparseStart,&
            & this%orb, this%potential, this%coord, this%derivs, this%groundDerivs,&
            & this%tripletderivs, this%mixedderivs, this%iRhoPrim, this%thirdOrd,&
            & this%solvation, this%qDepExtPot, this%chrgForces, this%dispersion,&
            & this%rangeSep, this%SSqrReal, this%ints, this%denseDesc, this%deltaRhoOutSqr,&
            & this%halogenXCorrection, this%tHelical, this%coord0, this%deltaDftb)

        if (this%tCasidaForces) then
          this%derivs(:,:) = this%derivs + this%excitedDerivs
        end if
      end if

      call env%globalTimer%stopTimer(globalTimers%forceCalc)

      call updateDerivsByPlumed(env, this%plumedCalc, this%nAtom, iGeoStep, this%derivs,&
          & this%dftbEnergy(this%deltaDftb%iDeterminant)%EMermin, this%coord0, this%mass,&
          & this%tPeriodic, this%latVec)

      if (this%tStress) then
        call env%globalTimer%startTimer(globalTimers%stressCalc)
        if (allocated(this%reks)) then
          call getReksStress(env, this%denseDesc, this%scc, this%nonSccDeriv, this%skHamCont,&
              & this%skOverCont, this%repulsive, this%neighbourList, this%nNeighbourSk,&
              & this%species, this%img2CentCell, this%iSparseStart, this%orb,&
              & this%dispersion, this%coord, this%q0, this%invLatVec, this%cellVol,&
              & this%totalStress, this%totalLatDeriv, this%intPressure, this%reks)
        else
          call getStress(env, this%scc, this%tblite, this%thirdOrd, this%isExtField,&
              & this%nonSccDeriv, this%rhoPrim, this%ERhoPrim, this%qOutput, this%q0,&
              & this%skHamCont, this%skOverCont, this%repulsive, this%neighbourList,&
              & this%nNeighbourSk, this%species, this%img2CentCell, this%iSparseStart,&
              & this%orb, this%potential, this%coord, this%latVec, this%invLatVec,&
              & this%cellVol, this%coord0, this%totalStress, this%totalLatDeriv,&
              & this%intPressure, this%iRhoPrim, this%solvation, this%dispersion,&
              & this%halogenXCorrection, this%deltaDftb, this%tripletStress, this%mixedStress)
        end if
        call env%globalTimer%stopTimer(globalTimers%stressCalc)

      end if

    end if

    if (this%tWriteDetailedOut  .and. this%deltaDftb%nDeterminant() == 1) then
      call writeDetailedOut4(this%fdDetailedOut%unit, this%tSccCalc, tConverged, this%isXlbomd,&
          & this%isLinResp, this%isGeoOpt .or. allocated(this%geoOpt), this%tMD,&
          & this%tPrintForces, this%tStress,&
          & this%tPeriodic, this%dftbEnergy(this%deltaDftb%iDeterminant), this%totalStress,&
          & this%totalLatDeriv, this%derivs, this%chrgForces, this%indMovedAtom, this%cellVol,&
          & this%intPressure, this%geoOutFile, this%iAtInCentralRegion)
    end if

    if (this%tSccCalc .and. allocated(this%electrostatPot)&
        & .and. (.not. (this%isGeoOpt .or. allocated(this%geoOpt) .or. this%tMD)&
        & .or. needsRestartWriting(this%isGeoOpt .or. allocated(this%geoOpt), this%tMd,&
        & iGeoStep, this%nGeoSteps, this%restartFreq))) then
      call this%electrostatPot%evaluate(env, this%scc, this%eField)
      call writeEsp(this%electrostatPot, env, iGeoStep, this%nGeoSteps)
    end if

  end subroutine processGeometry


  !> Process geometry for constrains
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
    real(dp), intent(in) :: totalLatDerivs(:,:)

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


  !> Next geometry step from driver
  subroutine getNextGeometry(this, env, iGeoStep, tWriteRestart, constrLatDerivs,&
      & tCoordStep, tGeomEnd, tStopDriver, iLatGeoStep, tempIon, tExitGeoOpt)

    !> Global variables
    type(TDftbPlusMain), intent(inout) :: this

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
    logical :: tCoordEnd, converged

    ! initially assume that coordinates and lattice vectors won't be updated
    this%tCoordsChanged = .false.
    this%tLatticeChanged = .false.

    tExitGeoOpt = .false.

    if (this%tDerivs) then
      call getNextDerivStep(this%derivDriver, this%derivs, this%indMovedAtom, &
           & this%indDerivAtom, this%coord0, tGeomEnd)
      if (tGeomEnd) then
        call env%globalTimer%stopTimer(globalTimers%postSCC)
        tExitGeoOpt = .true.
        return
      end if
      this%tCoordsChanged = .true.
    else if (allocated(this%geoOpt)) then
      block
        real(dp) :: energy, ediff, dnorm, damax, gnorm, gamax
        logical :: econv, dconv, gconv

        dnorm = norm2(this%displ)
        damax = maxval(abs(this%displ))

        call this%filter%transformDerivative(this%coord0, this%latVec, &
            & this%derivs, this%totalStress, this%gcurr)
        call this%geoOpt%step(this%dftbEnergy(this%deltaDftb%iFinal)%Emermin, &
            & this%gcurr, this%displ)
        call this%filter%transformStructure(this%coord0, this%latVec, this%displ)

        energy = this%dftbEnergy(this%deltaDftb%iFinal)%Emermin
        ediff = energy - this%elast
        gnorm = norm2(this%gcurr)
        gamax = maxval(abs(this%gcurr))
        write(stdOut, '(a)') stepSummary(energy, ediff, gnorm, gamax, dnorm, damax, 2)

        econv = ediff <= epsilon(0.0_dp) .and. abs(ediff) < this%optTol%energy
        dconv = dnorm / this%filter%nvar < this%optTol%dispNorm .and. damax < this%optTol%dispElem
        gconv = gnorm / this%filter%nvar < this%optTol%gradNorm .and. gamax < this%optTol%gradElem
        converged = econv .and. gconv .and. dconv
      end block

      this%elast = this%dftbEnergy(this%deltaDftb%iFinal)%Emermin

      if (tGeomEnd) then
        call env%globalTimer%stopTimer(globalTimers%postSCC)
        tExitGeoOpt = .true.
        return
      end if
      tGeomEnd = converged
      this%tCoordsChanged = .true.
      this%tLatticeChanged = this%tPeriodic .and. this%filter%lattice
    else if (this%isGeoOpt) then
      this%tCoordsChanged = .true.
      if (tCoordStep) then
        call getNextCoordinateOptStep(this%pGeoCoordOpt, this%dftbEnergy(this%deltaDftb%iFinal),&
            & this%derivs, this%indMovedAtom, this%coord0, diffGeo, tCoordEnd,&
            & .not. this%tCasidaForces)
        if (.not. this%tLatOpt) then
          tGeomEnd = tCoordEnd
        end if
        if (.not. tGeomEnd .and. tCoordEnd .and. diffGeo < tolSameDist) then
          tCoordStep = .false.
        end if
      else
        call getNextLatticeOptStep(this%pGeoLatOpt, this%dftbEnergy(this%deltaDftb%iDeterminant),&
            & constrLatDerivs, this%origLatVec, this%tLatOptFixAng, this%tLatOptFixLen,&
            & this%tLatOptIsotropic, this%indMovedAtom, this%latVec, this%coord0, diffGeo, tGeomEnd)
        iLatGeoStep = iLatGeoStep + 1
        this%tLatticeChanged = .true.
        if (.not. tGeomEnd .and. this%tCoordOpt) then
          tCoordStep = .true.
          call reset(this%pGeoCoordOpt,&
              & reshape(this%coord0(:, this%indMovedAtom), [this%nMovedCoord]))
        end if
      end if
      if (tGeomEnd .and. diffGeo < tolSameDist) then
        call env%globalTimer%stopTimer(globalTimers%postSCC)
        tExitGeoOpt = .true.
        return
      end if
    else if (this%tMD) then
      ! New MD coordinates saved in a temporary variable, as writeCurrentGeometry() below
      ! needs the old ones to write out consistent geometries and velocities.
      this%newCoords(:,:) = this%coord0
      call getNextMdStep(this%pMdIntegrator, this%pMdFrame, this%temperatureProfile, this%derivs,&
          & this%movedMass, this%mass, this%cellVol, this%invLatVec, this%species0,&
          & this%indMovedAtom, this%tStress, this%tBarostat,&
          & this%dftbEnergy(this%deltaDftb%iDeterminant), this%newCoords, this%latVec,&
          & this%intPressure, this%totalStress, this%totalLatDeriv, this%velocities, tempIon)
      this%tCoordsChanged = .true.
      this%tLatticeChanged = this%tBarostat
      call printMdInfo(this%tSetFillingTemp, this%eField, this%tPeriodic, this%tempElec,&
          & tempIon, this%intPressure, this%extPressure,&
          & this%dftbEnergy(this%deltaDftb%iDeterminant))
      if (tWriteRestart) then
        if (this%tPeriodic) then
          this%cellVol = abs(determinant33(this%latVec))
          this%dftbEnergy(this%deltaDftb%iDeterminant)%EGibbs =&
              & this%dftbEnergy(this%deltaDftb%iDeterminant)%EMermin&
              & + this%extPressure * this%cellVol
        end if
        call writeMdOut2(this%fdMd%unit, this%tStress, this%tBarostat, this%tPeriodic,&
            & this%isLinResp, this%eField, this%tFixEf, this%tPrintMulliken,&
            & this%dftbEnergy(this%deltaDftb%iDeterminant), this%energiesCasida, this%latVec,&
            & this%cellVol, this%intPressure, this%extPressure, tempIon, this%qOutput, this%q0,&
            & this%dipoleMoment, this%eFieldScaling, this%dipoleMessage)
        call writeCurrentGeometry(this%geoOutFile, this%pCoord0Out, .false., .true., .true.,&
            & this%tFracCoord, this%tPeriodic, this%tHelical, this%tPrintMulliken, this%species0,&
            & this%speciesName, this%latVec, this%origin, iGeoStep, iLatGeoStep, this%nSpin,&
            & this%qOutput, this%velocities)
      end if
      this%coord0(:,:) = this%newCoords
      if (this%tWriteDetailedOut  .and. this%deltaDftb%nDeterminant() == 1) then
        call writeDetailedOut5(this%fdDetailedOut%unit, this%tPrintForces, this%tSetFillingTemp,&
            & this%tPeriodic, this%tStress, this%totalStress, this%totalLatDeriv,&
            & this%dftbEnergy(this%deltaDftb%iDeterminant), this%tempElec, this%extPressure,&
            & this%intPressure, tempIon)
      end if
    else if (this%tSocket .and. iGeoStep < this%nGeoSteps) then
      ! Only receive geometry from socket, if there are still geometry iterations left
    #:if WITH_SOCKETS
      call receiveGeometryFromSocket(env, this%socket, this%tPeriodic, this%coord0, this%latVec,&
          & this%tCoordsChanged, this%tLatticeChanged, tStopDriver)
    #:else
      call error("Internal error: code compiled without socket support")
    #:endif
    end if

  end subroutine getNextGeometry


  function stepSummary(energy, ediff, gnorm, gamax, dnorm, damax, prlevel) result(str)
    real(dp), intent(in) :: energy, ediff, gnorm, gamax, dnorm, damax
    integer, intent(in) :: prlevel
    character(len=:), allocatable :: str
    character(len=*), parameter :: nl = new_line('a')

    str = ""
    if (prlevel <= 0) return

    if (prlevel > 1) then
      str = str//nl//"total energy  "//format_string(energy, '(es14.7)')//" H    "//&
        &            "   energy change "//format_string(ediff, '(es14.7)')//" H"
      str = str//nl//"gradient norm "//format_string(gnorm, '(es14.7)')//" H/a0 "//&
        &            "   max. gradient "//format_string(gamax, '(es14.7)')//" H/a0"
      str = str//nl//"step length   "//format_string(dnorm, '(es14.7)')//" a0   "//&
        &            "   max. step     "//format_string(damax, '(es14.7)')//" a0"
    else
      str = str//nl//"total energy  "//format_string(energy, '(es14.7)')//" H"
    end if

  contains
    function format_string(val, format) result(str)
      real(dp), intent(in) :: val
      character(len=*), intent(in) :: format
      character(len=:), allocatable :: str

      character(len=128) :: buffer
      integer :: stat

      write(buffer, format, iostat=stat) val
      if (stat == 0) then
        str = trim(buffer)
      else
        str = "*"
      end if
    end function format_string
  end function stepSummary

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


  !> Does the operations that are necessary after a lattice vector update
  subroutine handleLatticeChange(latVecs, sccCalc, tblite, tStress, extPressure, mCutOff,&
      & repulsive, dispersion, solvation, cm5Cont, recVecs, recVecs2p, cellVol, recCellVol,&
      & extLatDerivs, cellVecs, rCellVecs, boundaryCond)

    !> lattice vectors
    real(dp), intent(in) :: latVecs(:,:)

    !> Module variables
    type(TScc), allocatable, intent(inout) :: sccCalc

    !> Library interface handler
    type(TTBLite), allocatable, intent(inout) :: tblite

    !> evaluate stress
    logical, intent(in) :: tStress

    !> External pressure
    real(dp), intent(in) :: extPressure

    !> Maximum distance for interactions
    real(dp), intent(inout) :: mCutOff

    !> Repulsive interaction
    class(TRepulsive), allocatable, intent(inout) :: repulsive

    !> Dispersion interactions object
    class(TDispersionIface), allocatable, intent(inout) :: dispersion

    !> Solvation model
    class(TSolvation), allocatable, intent(inout) :: solvation

    !> Charge model 5
    type(TChargeModel5), allocatable, intent(inout) :: cm5Cont

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

    !> Boundary conditions on the calculation
    type(TBoundaryConditions), intent(in) :: boundaryCond

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
      call sccCalc%updateLatVecs(latVecs, recVecs, boundaryCond, cellVol)
      mCutOff = max(mCutOff, sccCalc%getCutOff())
    end if
    if (allocated(tblite)) then
      call tblite%updateLatVecs(latVecs)
      mCutOff = max(mCutOff, tblite%getRCutOff())
    end if
    if (allocated(repulsive)) then
      call repulsive%updateLatVecs(latVecs)
    end if
    if (allocated(dispersion)) then
      call dispersion%updateLatVecs(latVecs)
      mCutOff = max(mCutOff, dispersion%getRCutOff())
    end if
    if (allocated(solvation)) then
      call solvation%updateLatVecs(latVecs)
      mCutOff = max(mCutOff, solvation%getRCutOff())
    end if
    if (allocated(cm5Cont)) then
       call cm5Cont%updateLatVecs(latVecs)
       mCutoff = max(mCutOff, cm5Cont%getRCutOff())
    end if
    call getCellTranslations(cellVecs, rCellVecs, latVecs, recVecs2p, mCutOff)

  end subroutine handleLatticeChange


  !> Does the operations that are necessary after atomic coordinates change
  subroutine handleCoordinateChange(env, boundaryCond, coord0, latVec, invLatVec, species0, cutOff,&
      & orb, tPeriodic, tHelical, sccCalc, tblite, repulsive, dispersion, solvation, thirdOrd,&
      & rangeSep, reks, img2CentCell, iCellVec, neighbourList, nAllAtom, coord0Fold, coord,&
      & species, rCellVec, nNeighbourSK, nNeighbourLC, ints, H0, rhoPrim, iRhoPrim,&
      & ERhoPrim, iSparseStart, cm5Cont, errStatus)

    !> Environment settings
    type(TEnvironment), intent(in) :: env

    !> Boundary conditions on the calculation
    type(TBoundaryConditions), intent(in) :: boundaryCond

    !> Central cell coordinates
    real(dp), intent(in) :: coord0(:,:)

    !> Lattice vectors if periodic
    real(dp), intent(in) :: latVec(:,:)

    !> Inverse of the lattice vectors
    real(dp), intent(in) :: invLatVec(:,:)

    !> chemical species of central cell atoms
    integer, intent(in) :: species0(:)

    !> Longest cut-off distances that neighbour maps are generated for
    type(TCutoffs), intent(in) :: cutOff

    !> Atomic orbital information
    type(TOrbitals), intent(in) :: orb

    !> Is the geometry periodic
    logical, intent(in) :: tPeriodic

    !> Is the geometry helical
    logical, intent(in) :: tHelical

    !> SCC module internal variables
    type(TScc), allocatable, intent(inout) :: sccCalc

    !> Library interface handler
    type(TTBLite), allocatable, intent(inout) :: tblite

    !> Repulsive
    class(TRepulsive), allocatable, intent(inout) :: repulsive

    !> Dispersion interactions
    class(TDispersionIface), allocatable, intent(inout) :: dispersion

    !> Solvation model
    class(TSolvation), allocatable, intent(inout) :: solvation

    !> Third order SCC interactions
    type(TThirdOrder), allocatable, intent(inout) :: thirdOrd

    !> Range separation contributions
    type(TRangeSepFunc), allocatable, intent(inout) :: rangeSep

    !> data type for REKS
    type(TReksCalc), allocatable, intent(inout) :: reks

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

    !> Number of neighbours for each of the atoms for the exchange contributions in the long range
    !> functional
    integer, intent(inout), allocatable :: nNeighbourLC(:)

    !> Integral container
    type(TIntegral), intent(inout) :: ints

    !> Non-SCC hamiltonian storage
    real(dp), allocatable, intent(inout) :: h0(:)

    !> Sparse density matrix storage
    real(dp), allocatable, intent(inout) :: rhoPrim(:,:)

    !> Imaginary part of sparse density matrix storage
    real(dp), allocatable, intent(inout) :: iRhoPrim(:,:)

    !> energy weighted density matrix storage
    real(dp), allocatable, intent(inout) :: ERhoPrim(:)

    !> index array for location of atomic blocks in large sparse arrays
    integer, allocatable, intent(inout) :: iSparseStart(:,:)

    !> Charge model 5
    type(TChargeModel5), allocatable, intent(inout) :: cm5Cont

    !> Status of operation
    type(TStatus), intent(out) :: errStatus

    !> Total size of orbitals in the sparse data structures, where the decay of the overlap sets the
    !> sparsity pattern
    integer :: sparseSize

    coord0Fold(:,:) = coord0
    call boundaryCond%foldCoordsToCell(coord0Fold, latVec)

    if (tHelical) then
      call updateNeighbourListAndSpecies(env, coord, species, img2CentCell, iCellVec,&
          & neighbourList, nAllAtom, coord0Fold, species0, cutoff%mCutoff, rCellVec,&
          & helicalBoundConds=latVec)
    else
      call updateNeighbourListAndSpecies(env, coord, species, img2CentCell, iCellVec,&
          & neighbourList, nAllAtom, coord0Fold, species0, cutoff%mCutOff, rCellVec)
    end if

    call getNrOfNeighboursForAll(nNeighbourSK, neighbourList, cutoff%skCutOff)

    call getSparseDescriptor(neighbourList%iNeighbour, nNeighbourSK, img2CentCell, orb,&
        & iSparseStart, sparseSize)
    call reallocateSparseArrays(sparseSize, reks, ints, H0,&
        & rhoPrim, iRhoPrim, ERhoPrim)

    if (allocated(nNeighbourLC)) then
      call getNrOfNeighboursForAll(nNeighbourLC, neighbourList, cutoff%lcCutOff)
    end if

    if (allocated(sccCalc)) then
      call sccCalc%updateCoords(env, coord0, coord, species, neighbourList)
    end if

    if (allocated(tblite)) then
      call tblite%updateCoords(env, neighbourList, img2CentCell, coord, species)
    end if

    if (allocated(repulsive)) then
      call repulsive%updateCoords(coord, species, img2CentCell, neighbourList)
    end if

    if (allocated(dispersion)) then
      call dispersion%updateCoords(env, neighbourList, img2CentCell, coord, species0, errStatus)
      @:PROPAGATE_ERROR(errStatus)
    end if
    if (allocated(solvation)) then
      call solvation%updateCoords(env, neighbourList, img2CentCell, coord, species0)
    end if
    if (allocated(thirdOrd)) then
      call thirdOrd%updateCoords(neighbourList, species)
    end if
    if (allocated(rangeSep)) then
      call rangeSep%updateCoords(coord0)
    end if
    if (allocated(cm5Cont)) then
       call cm5Cont%updateCoords(neighbourList, img2CentCell, coord, species)
    end if


  end subroutine handleCoordinateChange


#:if WITH_TRANSPORT

  !> Initialise transport
  subroutine setupNegfStuff(negfInt, denseDescr, transpar, ginfo, neighbourList, nNeighbourSK,&
      & img2CentCell, orb)

    !> NEGF interface
    type(TNegfInt), intent(inout) :: negfInt

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
    call negfInt%setup_csr(denseDescr%iAtomStart, neighbourList%iNeighbour, nNeighbourSK,&
        & img2CentCell, orb)

    call negfInt%setup_str(denseDescr, transpar, ginfo%greendens, neighbourList%iNeighbour,&
        & nNeighbourSK, img2CentCell)

    call negfInt%setup_dephasing(ginfo%tundos)  !? why tundos

  end subroutine setupNegfStuff

#:endif


  !> Decides, whether restart file should be written during the run.
  function needsRestartWriting(isGeoOpt, tMd, iGeoStep, nGeoSteps, restartFreq)&
      & result(tWriteRestart)

    !> Are geometries being optimised
    logical, intent(in) :: isGeoOpt

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

    if (restartFreq > 0 .and. (isGeoOpt .or. tMD)) then
      tWriteRestart = (iGeoStep == nGeoSteps .or. (mod(iGeoStep, restartFreq) == 0))
    else
      tWriteRestart = .false.
    end if

  end function needsRestartWriting


  !> Ensures that sparse array have enough storage to hold all necessary elements.
  subroutine reallocateSparseArrays(sparseSize, reks, ints,&
      & H0, rhoPrim, iRhoPrim, ERhoPrim)

    !> Size of the sparse overlap
    integer, intent(in) :: sparseSize

    !> data type for REKS
    type(TReksCalc), allocatable, intent(inout) :: reks

    !> Integral container
    type(TIntegral), intent(inout) :: ints

    !> Sparse storage for non-SCC hamiltonian
    real(dp), allocatable, intent(inout) :: H0(:)

    !> Sparse storage for density matrix
    real(dp), allocatable, intent(inout) :: rhoPrim(:,:)

    !> Sparse storage for imaginary part of density matrix (not reallocated if not initially
    !> allocated)
    real(dp), allocatable, intent(inout) :: iRhoPrim(:,:)

    !> Sparse storage for energy weighted density matrix (not reallocated if not initially
    !> allocated)
    real(dp), allocatable, intent(inout) :: ERhoPrim(:)

    integer :: nSpin, nDipole, nQuadrupole

    #:block DEBUG_CODE
      @:ASSERT(size(H0) == size(ints%overlap))
      if (.not. allocated(reks)) then
        @:ASSERT(size(ints%hamiltonian, dim=1) == size(ints%overlap))
        @:ASSERT(all(shape(rhoPrim) == shape(ints%hamiltonian)))
        if (allocated(iRhoPrim)) then
          @:ASSERT(all(shape(iRhoPrim) == shape(rhoPrim)))
          @:ASSERT(all(shape(ints%iHamiltonian) == shape(ints%hamiltonian)))
        end if
        if (allocated(ERhoPrim)) then
          @:ASSERT(size(ERhoPrim) == size(rhoPrim, dim=1))
        end if
      end if
      @:ASSERT(allocated(ints%dipoleKet) .eqv. allocated(ints%dipoleBra))
      if (allocated(ints%dipoleKet)) then
        @:ASSERT(size(ints%overlap) == size(ints%dipoleKet, 2))
        @:ASSERT(all(shape(ints%dipoleKet) == shape(ints%dipoleBra)))
      end if
      @:ASSERT(allocated(ints%quadrupoleKet) .eqv. allocated(ints%quadrupoleBra))
      if (allocated(ints%quadrupoleKet)) then
        @:ASSERT(size(ints%overlap) == size(ints%quadrupoleKet, 2))
        @:ASSERT(all(shape(ints%quadrupoleKet) == shape(ints%quadrupoleBra)))
      end if
    #:endblock DEBUG_CODE

    if (allocated(reks)) then
      if (size(ints%overlap, dim=1) == sparseSize) then
        ! When the size of sparse matrices are different,
        ! phase of MOs can affect gradient of REKS
        return
      end if
    else
      if (size(ints%hamiltonian, dim=1) >= sparseSize) then
        ! Sparse matrices are big enough
        return
      end if
    end if

    nSpin = size(rhoPrim, dim=2)
    if (.not. allocated(reks)) then
      deallocate(ints%hamiltonian)
    end if
    deallocate(ints%overlap)
    deallocate(H0)
    deallocate(rhoPrim)
    if (.not. allocated(reks)) then
      allocate(ints%hamiltonian(sparseSize, nSpin))
    end if
    allocate(ints%overlap(sparseSize))
    allocate(H0(sparseSize))
    allocate(rhoPrim(sparseSize, nSpin))
    if (allocated(iRhoPrim)) then
      deallocate(iRhoPrim)
      deallocate(ints%iHamiltonian)
      allocate(iRhoPrim(sparseSize, nSpin))
      allocate(ints%iHamiltonian(sparseSize, nSpin))
    end if
    if (allocated(ERhoPrim)) then
      deallocate(ERhoPrim)
      allocate(ERhoPrim(sparseSize))
    end if
    if (allocated(reks)) then
      call reks%reallocate(sparseSize)
    end if
    if (allocated(ints%dipoleKet)) then
      nDipole = size(ints%dipoleKet, 1)
      deallocate(ints%dipoleBra, ints%dipoleKet)
      allocate(ints%dipoleKet(nDipole, sparseSize))
      allocate(ints%dipoleBra(nDipole, sparseSize))
    end if
    if (allocated(ints%quadrupoleKet)) then
      nQuadrupole = size(ints%quadrupoleKet, 1)
      deallocate(ints%quadrupoleBra, ints%quadrupoleKet)
      allocate(ints%quadrupoleKet(nQuadrupole, sparseSize))
      allocate(ints%quadrupoleBra(nQuadrupole, sparseSize))
    end if

  end subroutine reallocateSparseArrays


  !> Initialise basic variables before the scc loop.
  subroutine initSccLoop(tSccCalc, xlbomdIntegrator, minSccIter, maxSccIter, sccTol, tConverged,&
      & tNegf, reks)

    !> Is this an SCC calculation?
    logical, intent(in) :: tSccCalc

    !> Details for extended Lagrange integrator (of used)
    type(TXLBOMD), allocatable, intent(inout) :: xlbomdIntegrator

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

    !> data type for REKS
    type(TReksCalc), allocatable, intent(inout) :: reks

    if (allocated(xlbomdIntegrator)) then
      call xlbomdIntegrator%getSCCParameters(minSccIter, maxSccIter, sccTol)
    end if

    tConverged = (.not. tSccCalc)

    if (allocated(reks)) then
      if (tSccCalc) then
        call printReksSccHeader(reks)
      end if
    else
      if (tSccCalc .and. .not. tNegf) then
        call printSccHeader()
      end if
    end if

  end subroutine initSccLoop


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
  subroutine getDensity(env, negfInt, iScc, denseDesc, ints, neighbourList, nNeighbourSK,&
      & iSparseStart, img2CentCell, iCellVec, cellVec, kPoint, kWeight, orb, tHelical, coord,&
      & species, electronicSolver, tRealHS, tSpinSharedEf, tSpinOrbit, tDualSpinOrbit, tFillKSep,&
      & tFixEf, tMulliken, iDistribFn, tempElec, nEl, parallelKS, Ef, mu, energy, rangeSep, eigen,&
      & filling, rhoPrim, xi, orbitalL, HSqrReal, SSqrReal, eigvecsReal, iRhoPrim, HSqrCplx,&
      & SSqrCplx, eigvecsCplx, rhoSqrReal, deltaRhoInSqr, deltaRhoOutSqr, qOutput, nNeighbourLC,&
      & tLargeDenseMatrices, deltaDftb, errStatus)

    !> Environment settings
    type(TEnvironment), intent(inout) :: env

    !> NEGF interface
    type(TNegfInt), intent(inout) :: negfInt

    !> SCC iteration counter (needed by GF)
    integer, intent(in) :: iSCC

    !> Dense matrix descriptor
    type(TDenseDescr), intent(in) :: denseDesc

    !> Integral container
    type(TIntegral), intent(in) :: ints

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

    !> Is the hamiltonian real (no k-points/molecule/gamma point)?
    logical, intent(in) :: tRealHS

    !> Is the Fermi level common across spin channels?
    logical, intent(in) :: tSpinSharedEf

    !> Are spin orbit interactions present
    logical, intent(in) :: tSpinOrbit

    !> Are block population spin orbit interactions present
    logical, intent(in) :: tDualSpinOrbit

    !> Fill k-points separately if true (no charge transfer across the BZ)
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
    type(TRangeSepFunc), allocatable, intent(inout) :: rangeSep

    !> eigenvalues (level, kpoint, spin)
    real(dp), intent(out) :: eigen(:,:,:)

    !> occupations (level, kpoint, spin)
    real(dp), intent(out) :: filling(:,:,:)

    !> sparse density matrix
    real(dp), intent(out) :: rhoPrim(:,:)

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

    !> Determinant derived type
    type(TDftbDeterminants), intent(inout) :: deltaDftb

    !> Status of operation
    type(TStatus), intent(out) :: errStatus

    integer :: nSpin, iKS, iSp, iK, nAtom
    complex(dp), allocatable :: rhoSqrCplx(:,:)
    logical :: tImHam

    nSpin = size(ints%hamiltonian, dim=2)
    tImHam = allocated(iRhoPrim)

    select case (electronicSolver%iSolver)

    case (electronicSolverTypes%GF)

      call env%globalTimer%startTimer(globalTimers%densityMatrix)
    #:if WITH_TRANSPORT
      call negfInt%calcdensity_green(iSCC, env, parallelKS%localKS, ints%hamiltonian, ints%overlap,&
          & neighbourlist%iNeighbour, nNeighbourSK, denseDesc%iAtomStart, iSparseStart,&
          & img2CentCell, iCellVec, cellVec, orb, kPoint, kWeight, mu, rhoPrim, energy%Eband, Ef,&
          & energy%E0, energy%TS)
    #:else
      call error("Internal error: getDensity : GF solver called although code compiled without&
          & transport")
    #:endif
      call ud2qm(rhoPrim)
      call env%globalTimer%stopTimer(globalTimers%densityMatrix)

    case (electronicSolverTypes%onlyTransport)

      call error("OnlyTransport solver cannot calculate the density matrix")

    case(electronicSolverTypes%qr, electronicSolverTypes%divideandconquer,&
        & electronicSolverTypes%relativelyrobust, electronicSolverTypes%elpa,&
        & electronicSolverTypes%magma_gvd)

      call getDensityFromDenseDiag(env, denseDesc, ints, neighbourList, nNeighbourSK,&
          & iSparseStart, img2CentCell, iCellVec, cellVec, kPoint, kWeight, orb,&
          & denseDesc%iAtomStart, tHelical, coord, species, electronicSolver, tRealHS,&
          & tSpinSharedEf, tSpinOrbit, tDualSpinOrbit, tFillKSep, tFixEf, tMulliken, iDistribFn,&
          & tempElec, nEl, parallelKS, Ef, energy, rangeSep, eigen, filling, rhoPrim, xi,&
          & orbitalL, HSqrReal, SSqrReal, eigvecsReal, iRhoPrim, HSqrCplx, SSqrCplx, eigvecsCplx,&
          & rhoSqrReal, deltaRhoInSqr, deltaRhoOutSqr, qOutput, nNeighbourLC, deltaDftb, errStatus)
      @:PROPAGATE_ERROR(errStatus)

    case(electronicSolverTypes%omm, electronicSolverTypes%pexsi, electronicSolverTypes%ntpoly,&
        &electronicSolverTypes%elpadm)

      call env%globalTimer%startTimer(globalTimers%densityMatrix)

      call electronicSolver%elsi%getDensity(env, denseDesc, ints%hamiltonian, ints%overlap,&
          & neighbourList, nNeighbourSK, iSparseStart, img2CentCell, iCellVec, cellVec, kPoint,&
          & kWeight, tHelical, orb, species, coord, tRealHS, tSpinSharedEf, tSpinOrbit,&
          & tDualSpinOrbit, tMulliken, parallelKS, Ef, energy, rhoPrim, energy%Eband, energy%TS,&
          & ints%iHamiltonian, xi, orbitalL, HSqrReal, SSqrReal, iRhoPrim, HSqrCplx, SSqrCplx)
      call env%globalTimer%stopTimer(globalTimers%densityMatrix)

    end select

  end subroutine getDensity


  !> Returns the density matrix using dense diagonalisation.
  subroutine getDensityFromDenseDiag(env, denseDesc, ints, neighbourList, nNeighbourSK,&
      & iSparseStart, img2CentCell, iCellVec, cellVec, kPoint, kWeight, orb, iAtomStart, tHelical,&
      & coord, species, electronicSolver, tRealHS, tSpinSharedEf, tSpinOrbit, tDualSpinOrbit,&
      & tFillKSep, tFixEf, tMulliken, iDistribFn, tempElec, nEl, parallelKS, Ef, energy, rangeSep,&
      & eigen, filling, rhoPrim, xi, orbitalL, HSqrReal, SSqrReal, eigvecsReal, iRhoPrim,&
      & HSqrCplx, SSqrCplx, eigvecsCplx, rhoSqrReal, deltaRhoInSqr, deltaRhoOutSqr, qOutput,&
      & nNeighbourLC, deltaDftb, errStatus)

    !> Environment settings
    type(TEnvironment), intent(inout) :: env

    !> Dense matrix descriptor
    type(TDenseDescr), intent(in) :: denseDesc

    !> Integral container
    type(TIntegral), intent(in) :: ints

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

    !> Start of atomic blocks in dense arrays
    integer, allocatable, intent(in) :: iAtomStart(:)

    !> Is the geometry helical
    logical, intent(in) :: tHelical

    !> Coordinates of all atoms including images
    real(dp), allocatable, intent(inout) :: coord(:,:)

    !> species of all atoms in the system
    integer, intent(in) :: species(:)

    !> Electronic solver information
    type(TElectronicSolver), intent(inout) :: electronicSolver

    !> Is the hamiltonian real (no k-points/molecule/gamma point)?
    logical, intent(in) :: tRealHS

    !> Is the Fermi level common across spin channels?
    logical, intent(in) :: tSpinSharedEf

    !> Are spin orbit interactions present
    logical, intent(in) :: tSpinOrbit

    !> Are block population spin orbit interactions present
    logical, intent(in) :: tDualSpinOrbit

    !> Fill k-points separately if true (no charge transfer across the BZ)
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
    type(TRangeSepFunc), allocatable, intent(inout) :: rangeSep

    !> eigenvalues (level, kpoint, spin)
    real(dp), intent(out) :: eigen(:,:,:)

    !> occupations (level, kpoint, spin)
    real(dp), intent(out) :: filling(:,:,:)

    !> sparse density matrix
    real(dp), intent(out) :: rhoPrim(:,:)

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

    !> Determinant derived type
    type(TDftbDeterminants), intent(inout) :: deltaDftb

    !> Status of operation
    type(TStatus), intent(out) :: errStatus

    integer :: nSpin

    nSpin = size(ints%hamiltonian, dim=2)
    call env%globalTimer%startTimer(globalTimers%diagonalization)
    if (nSpin /= 4) then
      if (tRealHS) then
        call buildAndDiagDenseRealHam(env, denseDesc, ints, species, neighbourList,&
            & nNeighbourSK, iSparseStart, img2CentCell, orb, iAtomStart, tHelical, coord,&
            & electronicSolver, parallelKS, rangeSep, deltaRhoInSqr, qOutput, nNeighbourLC,&
            & HSqrReal, SSqrReal, eigVecsReal, eigen(:,1,:), errStatus)
        @:PROPAGATE_ERROR(errStatus)
      else
        call buildAndDiagDenseCplxHam(env, denseDesc, ints, kPoint, neighbourList,&
            & nNeighbourSK, iSparseStart, img2CentCell, iCellVec, cellVec, electronicSolver,&
            & parallelKS, tHelical, orb, species, coord, HSqrCplx, SSqrCplx, eigVecsCplx, eigen,&
            & errStatus)
        @:PROPAGATE_ERROR(errStatus)
      end if
    else
      call buildAndDiagDensePauliHam(env, denseDesc, ints, kPoint, neighbourList,&
          & nNeighbourSK, iSparseStart, img2CentCell, iCellVec, cellVec, orb, electronicSolver,&
          & parallelKS, eigen(:,:,1), HSqrCplx, SSqrCplx, eigVecsCplx, errStatus, xi, species)
      @:PROPAGATE_ERROR(errStatus)
    end if
    call env%globalTimer%stopTimer(globalTimers%diagonalization)

    call getFillingsAndBandEnergies(eigen, nEl, nSpin, tempElec, kWeight, tSpinSharedEf,&
        & tFillKSep, tFixEf, iDistribFn, Ef, filling, energy%Eband, energy%TS, energy%E0, deltaDftb)

    call env%globalTimer%startTimer(globalTimers%densityMatrix)
    if (nSpin /= 4) then
      if (tRealHS) then
        call getDensityFromRealEigvecs(env, denseDesc, filling(:,1,:), neighbourList, nNeighbourSK,&
            & iSparseStart, img2CentCell, orb, species, denseDesc%iAtomStart, coord, tHelical,&
            & eigVecsReal, parallelKS, rhoPrim, SSqrReal, rhoSqrReal, deltaRhoOutSqr)
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
          & iSparseStart, img2CentCell, iCellVec, cellVec, species, parallelKS, deltaDftb,&
          & eigVecsCplx, SSqrCplx, energy, rhoPrim, xi, orbitalL, iRhoPrim)
      filling(:,:,1) = 0.5_dp * filling(:,:,1)
    end if
    call env%globalTimer%stopTimer(globalTimers%densityMatrix)

  end subroutine getDensityFromDenseDiag


  !> Builds and diagonalises dense Hamiltonians.
  subroutine buildAndDiagDenseRealHam(env, denseDesc, ints, species, neighbourList,&
      & nNeighbourSK, iSparseStart, img2CentCell, orb, iAtomStart, tHelical, coord,&
      & electronicSolver, parallelKS, rangeSep, deltaRhoInSqr, qOutput, nNeighbourLC, HSqrReal,&
      & SSqrReal, eigvecsReal, eigen, errStatus)

    !> Environment settings
    type(TEnvironment), intent(inout) :: env

    !> Dense matrix descriptor
    type(TDenseDescr), intent(in) :: denseDesc

    !> Integral container
    type(TIntegral), intent(in) :: ints

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

    !> Start of atomic blocks in dense arrays
    integer, allocatable, intent(in) :: iAtomStart(:)

    !> Is the geometry helical
    logical, intent(in) :: tHelical

    !> Coordinates of all atoms including images
    real(dp), allocatable, intent(inout) :: coord(:,:)

    !> Electronic solver information
    type(TElectronicSolver), intent(inout) :: electronicSolver

    !> K-points and spins to be handled
    type(TParallelKS), intent(in) :: parallelKS

    !>Data for rangeseparated calculation
    type(TRangeSepFunc), allocatable, intent(inout) :: rangeSep

    !> Change in density matrix during last rangesep SCC cycle
    real(dp), pointer, intent(in) :: deltaRhoInSqr(:,:,:)

    !> Output electrons
    real(dp), intent(inout) :: qOutput(:,:,:)

    !> Number of neighbours for each of the atoms for the exchange contributions in the long range
    !> functional
    integer, intent(in), allocatable :: nNeighbourLC(:)

    !> dense hamiltonian matrix
    real(dp), intent(out) :: HSqrReal(:,:)

    !> dense overlap matrix
    real(dp), intent(out) :: SSqrReal(:,:)

    !> Eigenvectors on eixt
    real(dp), intent(out) :: eigvecsReal(:,:,:)

    !> eigenvalues
    real(dp), intent(out) :: eigen(:,:)

    !> Status of operation
    type(TStatus), intent(out) :: errStatus

    integer :: iKS, iSpin

    eigen(:,:) = 0.0_dp
    do iKS = 1, parallelKS%nLocalKS
      iSpin = parallelKS%localKS(2, iKS)
    #:if WITH_SCALAPACK
      call env%globalTimer%startTimer(globalTimers%sparseToDense)
      if (tHelical) then
        call unpackHSHelicalRealBlacs(env%blacs, ints%hamiltonian(:,iSpin),&
            & neighbourList%iNeighbour, nNeighbourSK, iSparseStart, img2CentCell, orb, species,&
            & coord, denseDesc, HSqrReal)
        if (.not. electronicSolver%hasCholesky(1)) then
          call unpackHSHelicalRealBlacs(env%blacs, ints%overlap, neighbourList%iNeighbour,&
              & nNeighbourSK, iSparseStart, img2CentCell, orb, species, coord, denseDesc, SSqrReal)
        end if
      else
        call unpackHSRealBlacs(env%blacs, ints%hamiltonian(:,iSpin), neighbourList%iNeighbour,&
            & nNeighbourSK, iSparseStart, img2CentCell, denseDesc, HSqrReal)
        if (.not. electronicSolver%hasCholesky(1)) then
          call unpackHSRealBlacs(env%blacs, ints%overlap, neighbourList%iNeighbour, nNeighbourSK,&
              & iSparseStart, img2CentCell, denseDesc, SSqrReal)
        end if
      end if
      call env%globalTimer%stopTimer(globalTimers%sparseToDense)
      call diagDenseMtxBlacs(electronicSolver, 1, 'V', denseDesc%blacsOrbSqr, HSqrReal, SSqrReal,&
          & eigen(:,iSpin), eigvecsReal(:,:,iKS), errStatus)
      @:PROPAGATE_ERROR(errStatus)
    #:else
      call env%globalTimer%startTimer(globalTimers%sparseToDense)
      if (tHelical) then
        call unpackHelicalHS(HSqrReal, ints%hamiltonian(:,iSpin), neighbourList%iNeighbour,&
            & nNeighbourSK, denseDesc%iAtomStart, iSparseStart, img2CentCell, orb, species, coord)
        call unpackHelicalHS(SSqrReal, ints%overlap, neighbourList%iNeighbour, nNeighbourSK,&
            & denseDesc%iAtomStart, iSparseStart, img2CentCell, orb, species, coord)
      else
        call unpackHS(HSqrReal, ints%hamiltonian(:,iSpin), neighbourList%iNeighbour, nNeighbourSK,&
            & denseDesc%iAtomStart, iSparseStart, img2CentCell)
        call unpackHS(SSqrReal, ints%overlap, neighbourList%iNeighbour, nNeighbourSK,&
            & denseDesc%iAtomStart, iSparseStart, img2CentCell)
      end if
      call env%globalTimer%stopTimer(globalTimers%sparseToDense)

      ! Add rangeseparated contribution
      ! Assumes deltaRhoInSqr only used by rangeseparation
      ! Should this be used elsewhere, need to pass isRangeSep
      if (allocated(rangeSep)) then
        call denseMulliken(deltaRhoInSqr, SSqrReal, denseDesc%iAtomStart, qOutput)
        call rangeSep%addLRHamiltonian(env, deltaRhoInSqr(:,:,iSpin), ints%overlap,&
            & neighbourList%iNeighbour,  nNeighbourLC, denseDesc%iAtomStart, iSparseStart,&
            & orb, HSqrReal, SSqrReal)
      end if

      call diagDenseMtx(env, electronicSolver, 'V', HSqrReal, SSqrReal, eigen(:,iSpin), errStatus)
      @:PROPAGATE_ERROR(errStatus)
      eigvecsReal(:,:,iKS) = HSqrReal
    #:endif
    end do

  #:if WITH_SCALAPACK
    ! Distribute all eigenvalues to all nodes via global summation
    call mpifx_allreduceip(env%mpi%interGroupComm, eigen, MPI_SUM)
  #:endif

  end subroutine buildAndDiagDenseRealHam


  !> Builds and diagonalises dense k-point dependent Hamiltonians.
  subroutine buildAndDiagDenseCplxHam(env, denseDesc, ints, kPoint, neighbourList,&
      & nNeighbourSK, iSparseStart, img2CentCell, iCellVec, cellVec, electronicSolver, parallelKS,&
      & tHelical, orb, species, coord, HSqrCplx, SSqrCplx, eigvecsCplx, eigen, errStatus)

    !> Environment settings
    type(TEnvironment), intent(inout) :: env

    !> Dense matrix descriptor
    type(TDenseDescr), intent(in) :: denseDesc

    !> Integral container
    type(TIntegral), intent(in) :: ints

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

    !> dense hamiltonian matrix
    complex(dp), intent(out) :: HSqrCplx(:,:)

    !> dense overlap matrix
    complex(dp), intent(out) :: SSqrCplx(:,:)

    !> Complex eigenvectors
    complex(dp), intent(out) :: eigvecsCplx(:,:,:)

    !> eigenvalues
    real(dp), intent(out) :: eigen(:,:,:)

    !> Status of operation
    type(TStatus), intent(out) :: errStatus

    integer :: iKS, iK, iSpin

    eigen(:,:,:) = 0.0_dp
    do iKS = 1, parallelKS%nLocalKS
      iK = parallelKS%localKS(1, iKS)
      iSpin = parallelKS%localKS(2, iKS)
    #:if WITH_SCALAPACK
      call env%globalTimer%startTimer(globalTimers%sparseToDense)
      if (tHelical) then
        call unpackHSHelicalCplxBlacs(env%blacs, ints%hamiltonian(:,iSpin), kPoint(:,iK),&
            & neighbourList%iNeighbour, nNeighbourSK, iCellVec, cellVec, iSparseStart,&
            & img2CentCell, orb, species, coord, denseDesc, HSqrCplx)
        if (.not. electronicSolver%hasCholesky(iKS)) then
          call unpackHSHelicalCplxBlacs(env%blacs, ints%overlap, kPoint(:,iK),&
              & neighbourList%iNeighbour, nNeighbourSK, iCellVec, cellVec, iSparseStart,&
              & img2CentCell, orb, species, coord, denseDesc, SSqrCplx)
        end if
      else
        call unpackHSCplxBlacs(env%blacs, ints%hamiltonian(:,iSpin), kPoint(:,iK),&
            & neighbourList%iNeighbour, nNeighbourSK, iCellVec, cellVec, iSparseStart,&
            & img2CentCell, denseDesc, HSqrCplx)
        if (.not. electronicSolver%hasCholesky(iKS)) then
          call unpackHSCplxBlacs(env%blacs, ints%overlap, kPoint(:,iK), neighbourList%iNeighbour,&
              & nNeighbourSK, iCellVec, cellVec, iSparseStart, img2CentCell, denseDesc, SSqrCplx)
        end if
      end if
      call env%globalTimer%stopTimer(globalTimers%sparseToDense)
      call diagDenseMtxBlacs(env, electronicSolver, iKS, 'V', denseDesc%blacsOrbSqr, HSqrCplx,&
          & SSqrCplx, eigen(:,iK,iSpin), eigvecsCplx(:,:,iKS), errStatus)
      @:PROPAGATE_ERROR(errStatus)
    #:else
      call env%globalTimer%startTimer(globalTimers%sparseToDense)
      if (tHelical) then
        call unpackHelicalHS(HSqrCplx, ints%hamiltonian(:,iSpin), kPoint(:,iK),&
            & neighbourList%iNeighbour, nNeighbourSK, iCellVec, cellVec, denseDesc%iAtomStart,&
            & iSparseStart, img2CentCell, orb, species, coord)
        call unpackHelicalHS(SSqrCplx, ints%overlap, kPoint(:,iK), neighbourList%iNeighbour,&
            & nNeighbourSK, iCellVec, cellVec, denseDesc%iAtomStart, iSparseStart, img2CentCell,&
            & orb, species, coord)
      else
        call unpackHS(HSqrCplx, ints%hamiltonian(:,iSpin), kPoint(:,iK), neighbourList%iNeighbour,&
            & nNeighbourSK, iCellVec, cellVec, denseDesc%iAtomStart, iSparseStart, img2CentCell)
        call unpackHS(SSqrCplx, ints%overlap, kPoint(:,iK), neighbourList%iNeighbour, nNeighbourSK,&
            & iCellVec, cellVec, denseDesc%iAtomStart, iSparseStart, img2CentCell)
      end if
      call env%globalTimer%stopTimer(globalTimers%sparseToDense)
      call diagDenseMtx(env, electronicSolver, 'V', HSqrCplx, SSqrCplx, eigen(:,iK,iSpin),&
          & errStatus)
      @:PROPAGATE_ERROR(errStatus)
      eigvecsCplx(:,:,iKS) = HSqrCplx
    #:endif
    end do

  #:if WITH_SCALAPACK
    call mpifx_allreduceip(env%mpi%interGroupComm, eigen, MPI_SUM)
  #:endif

  end subroutine buildAndDiagDenseCplxHam


  !> Builds and diagonalizes Pauli two-component Hamiltonians.
  subroutine buildAndDiagDensePauliHam(env, denseDesc, ints, kPoint, neighbourList,&
      & nNeighbourSK, iSparseStart, img2CentCell, iCellVec, cellVec, orb, electronicSolver,&
      & parallelKS, eigen, HSqrCplx, SSqrCplx, eigvecsCplx, errStatus, xi, species)

    !> Environment settings
    type(TEnvironment), intent(inout) :: env

    !> Dense matrix descriptor
    type(TDenseDescr), intent(in) :: denseDesc

    !> Integral container
    type(TIntegral), intent(in) :: ints

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

    !> dense hamiltonian matrix
    complex(dp), intent(out) :: HSqrCplx(:,:)

    !> dense overlap matrix
    complex(dp), intent(out) :: SSqrCplx(:,:)

    !> eigenvectors
    complex(dp), intent(out) :: eigvecsCplx(:,:,:)

    !> Status of operation
    type(TStatus), intent(out) :: errStatus

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
      if (allocated(ints%iHamiltonian)) then
        call unpackHPauliBlacs(env%blacs, ints%hamiltonian, kPoint(:,iK), neighbourList%iNeighbour,&
            & nNeighbourSK, iCellVec, cellVec, iSparseStart, img2CentCell, orb%mOrb, denseDesc,&
            & HSqrCplx, iorig=ints%iHamiltonian)
      else
        call unpackHPauliBlacs(env%blacs, ints%hamiltonian, kPoint(:,iK), neighbourList%iNeighbour,&
            & nNeighbourSK, iCellVec, cellVec, iSparseStart, img2CentCell, orb%mOrb, denseDesc,&
            & HSqrCplx)
      end if
      if (.not. electronicSolver%hasCholesky(iKS)) then
        call unpackSPauliBlacs(env%blacs, ints%overlap, kPoint(:,iK), neighbourList%iNeighbour,&
            & nNeighbourSK, iCellVec, cellVec, iSparseStart, img2CentCell, orb%mOrb, denseDesc,&
            & SSqrCplx)
      end if
    #:else
      if (allocated(ints%iHamiltonian)) then
        call unpackHPauli(ints%hamiltonian, kPoint(:,iK), neighbourList%iNeighbour, nNeighbourSK,&
            & iSparseStart, denseDesc%iAtomStart, img2CentCell, iCellVec, cellVec, HSqrCplx,&
            & iHam=ints%iHamiltonian)
      else
        call unpackHPauli(ints%hamiltonian, kPoint(:,iK), neighbourList%iNeighbour, nNeighbourSK,&
            & iSparseStart, denseDesc%iAtomStart, img2CentCell, iCellVec, cellVec, HSqrCplx)
      end if
      call unpackSPauli(ints%overlap, kPoint(:,iK), neighbourList%iNeighbour, nNeighbourSK,&
          & denseDesc%iAtomStart, iSparseStart, img2CentCell, iCellVec, cellVec, SSqrCplx)
    #:endif
      if (allocated(xi) .and. .not. allocated(ints%iHamiltonian)) then
        call addOnsiteSpinOrbitHam(env, xi, species, orb, denseDesc, HSqrCplx)
      end if
      call env%globalTimer%stopTimer(globalTimers%sparseToDense)
    #:if WITH_SCALAPACK
      call diagDenseMtxBlacs(env, electronicSolver, iKS, 'V', denseDesc%blacsOrbSqr, HSqrCplx,&
          & SSqrCplx, eigen(:,iK), eigvecsCplx(:,:,iKS), errStatus)
      @:PROPAGATE_ERROR(errStatus)
    #:else
      call diagDenseMtx(env, electronicSolver, 'V', HSqrCplx, SSqrCplx, eigen(:,iK), errStatus)
      @:PROPAGATE_ERROR(errStatus)
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

    !> Start of atomic blocks in dense arrays
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
        call makeDensityMatrix(work, eigvecs(:,:,iKS), filling(:,iSpin))
        call env%globalTimer%startTimer(globalTimers%denseToSparse)
        if (tHelical) then
          call packHelicalHS(rhoPrim(:,iSpin), work, neighbourlist%iNeighbour, nNeighbourSK,&
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
          call packHelicalHS(rhoPrim(:,iSpin), deltaRhoOutSqr(:,:,iSpin), neighbourlist%iNeighbour,&
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
      call makeDensityMatrix(work, eigvecs(:,:,iKS), filling(:,iK,iSpin))
      call env%globalTimer%startTimer(globalTimers%denseToSparse)
      if (tHelical) then
        call packHelicalHS(rhoPrim(:,iSpin), work, kPoint(:,iK), kWeight(iK),&
            & neighbourList%iNeighbour, nNeighbourSK, orb%mOrb, iCellVec, cellVec,&
            & denseDesc%iAtomStart, iSparseStart, img2CentCell, orb, species, coord)
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
      & img2CentCell, iCellVec, cellVec, species, parallelKS, deltaDftb, eigvecs, work, dftbEnergy,&
      & rhoPrim, xi, orbitalL, iRhoPrim)

    !> Environment settings
    type(TEnvironment), intent(inout) :: env

    !> Dense matrix descriptor
    type(TDenseDescr), intent(in) :: denseDesc

    !> Is the hamiltonian real (no k-points/molecule/gamma point)?
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

    !> Determinant derived type
    type(TDftbDeterminants), intent(inout) :: deltaDftb

    !> eigenvectors
    complex(dp), intent(inout) :: eigvecs(:,:,:)

    !> work space array
    complex(dp), intent(inout) :: work(:,:)

    !> Energy contributions and total
    type(TEnergies), intent(inout) :: dftbEnergy

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
      dftbEnergy%atomLS(:) = 0.0_dp
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
        dftbEnergy%atomLS = dftbEnergy%atomLS + kWeight(iK) * rVecTemp
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
    call mpifx_allreduceip(env%mpi%globalComm, dftbEnergy%atomLS, MPI_SUM)
    if (tMulliken .and. tSpinOrbit .and. .not. tDualSpinOrbit) then
      call mpifx_allreduceip(env%mpi%globalComm, orbitalL, MPI_SUM)
    end if
    call env%globalTimer%stopTimer(globalTimers%denseToSparse)
  #:endif
    if (tSpinOrbit .and. .not. tDualSpinOrbit) then
      dftbEnergy%ELS = sum(dftbEnergy%atomLS)
    end if

  end subroutine getDensityFromPauliEigvecs


  !> Calculates electron fillings and resulting band energy terms.
  subroutine getFillingsAndBandEnergies(eigvals, nElectrons, nSpinBlocks, tempElec, kWeights,&
      & tSpinSharedEf, tFillKSep, tFixEf, iDistribFn, Ef, fillings, Eband, TS, E0, deltaDftb)

    !> Eigenvalue of each level, kpoint and spin channel
    real(dp), intent(inout) :: eigvals(:,:,:)

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

    !> Determinant derived type
    type(TDftbDeterminants), intent(inout) :: deltaDftb

    real(dp) :: EbandTmp(2), TSTmp(2), E0Tmp(2), EfTmp(2), nElecFill(2), kWeightTmp(2)
    integer :: nSpinHams, nKPoints, nLevels, iS, iK, iConfig

    nLevels = size(fillings, dim=1)
    nKPoints = size(fillings, dim=2)
    nSpinHams = size(fillings, dim=3)

    if (nSpinBlocks == 1) then
      ! Filling functions assume one electron per level, but for spin unpolarised we have two
      nElecFill(1) = 0.5_dp * nElectrons(1)
    else
      nElecFill(:nSpinHams) = nElectrons(:nSpinHams)
    end if

    if (tFixEf) then
      ! Fixed value of the Fermi level for each spin channel
      do iS = 1, nSpinHams
        call electronFill(Eband(iS:iS), fillings(:,:,iS:iS), TS(iS:iS), E0(iS:iS), Ef(iS),&
            & eigvals(:,:,iS:iS), tempElec, iDistribFn, kWeights)
      end do
    else if (nSpinHams == 2 .and. tSpinSharedEf) then
      ! Common Fermi level across two colinear spin channels
      TS(:) = 0.0_dp
      E0(:) = 0.0_dp
      Eband(:) = 0.0_dp
      call Efilling(Eband, Ef(1), TS, E0, fillings, eigvals, sum(nElecFill), tempElec, kWeights,&
          & iDistribFn)
      Ef(2) = Ef(1)
    else if (tFillKSep) then
      ! Every spin channel and every k-point filled up individually.
      Eband(:) = 0.0_dp
      Ef(:) = 0.0_dp
      TS(:) = 0.0_dp
      E0(:) = 0.0_dp
      kWeightTmp(:) = 1.0_dp
      do iK = 1, nKPoints
        call deltaDftb%detFilling(fillings(:,iK:iK,:), EBandTmp(:nSpinHams), EfTmp, TSTmp, E0Tmp, nElecFill,&
            & eigVals(:,iK:iK,:), tempElec, kWeightTmp(:nSpinHams), iDistribFn)
        Eband(:) = Eband + EbandTmp(:nSpinHams) * kWeights(iK)
        Ef(:) = Ef + EfTmp(:nSpinHams) * kWeights(iK)
        TS(:) = TS + TSTmp(:nSpinHams) * kWeights(iK)
        E0(:) = E0 + E0Tmp(:nSpinHams) * kWeights(iK)
      end do
    else
      call deltaDftb%detFilling(fillings, EBand, Ef, TS, E0, nElecFill, eigVals, tempElec,&
          & kWeights, iDistribFn)
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
  subroutine getMullikenPopulation(rhoPrim, ints, orb, neighbourList, nNeighbourSK, img2CentCell,&
      & iSparseStart, qOrb, iRhoPrim, qBlock, qiBlock, qNetAtom, multipoles)

    !> sparse density matrix
    real(dp), intent(in) :: rhoPrim(:,:)

    !> Integral container
    type(TIntegral), intent(in) :: ints

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

    !> Onsite Mulliken charges per atom
    real(dp), intent(inout), allocatable :: qNetAtom(:)

    !> Multipole moments
    type(TMultipole), intent(inout), optional :: multipoles

    integer :: iSpin

    qOrb(:,:,:) = 0.0_dp
    do iSpin = 1, size(rhoPrim, dim=2)
      call mulliken(qOrb(:,:,iSpin), ints%overlap, rhoPrim(:,iSpin), orb, neighbourList%iNeighbour,&
          & nNeighbourSK, img2CentCell, iSparseStart)
    end do

    if (allocated(qBlock)) then
      qBlock(:,:,:,:) = 0.0_dp
      do iSpin = 1, size(rhoPrim, dim=2)
        call mulliken(qBlock(:,:,:,iSpin), ints%overlap, rhoPrim(:,iSpin), orb,&
            & neighbourList%iNeighbour, nNeighbourSK, img2CentCell, iSparseStart)
      end do
    end if

    if (allocated(qiBlock)) then
      qiBlock(:,:,:,:) = 0.0_dp
      do iSpin = 1, size(iRhoPrim, dim=2)
        call skewMulliken(qiBlock(:,:,:,iSpin), ints%overlap, iRhoPrim(:,iSpin), orb,&
            & neighbourList%iNeighbour, nNeighbourSK, img2CentCell, iSparseStart)
      end do
    end if

    if (allocated(qNetAtom)) then
      call getOnsitePopulation(rhoPrim(:,1), orb, iSparseStart, qNetAtom)
    end if

    if (present(multipoles)) then

      if (allocated(multipoles%dipoleAtom)) then
        call getAtomicMultipolePopulation(multipoles%dipoleAtom, ints%dipoleBra, ints%dipoleKet, &
            & rhoPrim, orb, neighbourList%iNeighbour, nNeighbourSK, img2CentCell, &
            & iSparseStart)
      end if

      if (allocated(multipoles%quadrupoleAtom)) then
        call getAtomicMultipolePopulation(multipoles%quadrupoleAtom, ints%quadrupoleBra,&
            & ints%quadrupoleKet, rhoPrim, orb, neighbourList%iNeighbour, nNeighbourSK,&
            & img2CentCell, iSparseStart)
      end if

    end if

  end subroutine getMullikenPopulation


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
      & qInput, qInpRed, sccErrorQ, tConverged, dftbU, qBlockOut, iEqBlockDftbU, qBlockIn,&
      & qiBlockOut, iEqBlockDftbuLS, species0, qiBlockIn, iEqBlockOnSite, iEqBlockOnSiteLS,&
      & nIneqDip, nIneqQuad, iEqDipole, iEqQuadrupole, multipoleOut, multipoleInp)

    !> Environment settings
    type(TEnvironment), intent(in) :: env

    !> Charge mixing object
    type(TMixer), intent(inout) :: pChrgMixer

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

    !> minimum number of SCC iterations to perform
    integer, intent(in) :: minSccIter

    !> maximum number of SCC iterations before terminating loop
    integer, intent(in) :: maxSccIter

    !> Tolerance on SCC charges between input and output
    real(dp), intent(in) :: sccTol

    !> Should the SCC loop stop
    logical, intent(in) :: tStopScc

    !> are orbital potentials being used
    logical, intent(in) :: tMixBlockCharges

    !> Were initial charges read from disc?
    logical, intent(in) :: tReadChrg

    !> Resulting input charges for next SCC iteration
    real(dp), intent(inout) :: qInput(:,:,:)

    !> Equivalence reduced input charges
    real(dp), intent(inout) :: qInpRed(:)

    !> SCC error
    real(dp), intent(out) :: sccErrorQ

    !> Has the calculation converged
    logical, intent(out) :: tConverged

    !> Are there orbital potentials present
    type(TDftbU), intent(in), allocatable :: dftbU

    !> Dual output charges
    real(dp), intent(inout), allocatable :: qBlockOut(:,:,:,:)

    !> equivalence mapping for dual charge blocks
    integer, intent(in), allocatable :: iEqBlockDftbu(:,:,:,:)

    !> block charge input (if needed for orbital potentials)
    real(dp), intent(inout), allocatable :: qBlockIn(:,:,:,:)

    !> Imaginary part of block charges
    real(dp), intent(in), allocatable :: qiBlockOut(:,:,:,:)

    !> Equivalence mappings in the case of spin orbit and DFTB+U
    integer, intent(in), allocatable :: iEqBlockDftbuLS(:,:,:,:)

    !> atomic species for atoms
    integer, intent(in), allocatable :: species0(:)

    !> Imaginary part of block atomic input populations
    real(dp), intent(inout), allocatable :: qiBlockIn(:,:,:,:)

    !> Equivalences for onsite block corrections if needed
    integer, intent(in), allocatable :: iEqBlockOnSite(:,:,:,:)

    !> Equivalences for onsite block corrections if needed for imaginary elements
    integer, intent(in), allocatable :: iEqBlockOnSiteLS(:,:,:,:)

    !> Total number of inequivalent cumulative atomic dipole moments
    integer, intent(in) :: nIneqDip

    !> Total number of inequivalent cumulative atomic quadrupole moments
    integer, intent(in) :: nIneqQuad

    !> Equivalence relations between cumulative atomic dipole moments
    integer, intent(in) :: iEqDipole(:,:)

    !> Equivalence relations between cumulative atomic quadrupole moments
    integer, intent(in) :: iEqQuadrupole(:,:)

    !> Multipole moments
    type(TMultipole), intent(inout) :: multipoleInp

    !> Multipole moments
    type(TMultipole), intent(inout) :: multipoleOut

    real(dp), allocatable :: qDiffRed(:)
    integer :: nSpin, nMix

    nSpin = size(qOutput, dim=3)

    call reduceCharges(orb, nIneqOrb, iEqOrbitals, qOutput, qOutRed, qBlockOut, dftbU,&
        & iEqBlockDftbu, qiBlockOut, iEqBlockDftbuLS, iEqBlockOnSite, iEqBlockOnSiteLS)
    if (nIneqDip > 0) then
      ! FIXME: Assumes we always mix all dipole moments
      nMix = nIneqOrb
      @:ASSERT(allocated(multipoleOut%dipoleAtom))
      qOutRed(nMix+1:nMix+nIneqDip) = reshape(multipoleOut%dipoleAtom, [nIneqDip])
    end if
    if (nIneqQuad > 0) then
      ! FIXME: Assumes we always mix all quadrupole moments
      nMix = nIneqOrb + nIneqDip
      @:ASSERT(allocated(multipoleOut%quadrupoleAtom))
      qOutRed(nMix+1:nMix+nIneqQuad) = reshape(multipoleOut%quadrupoleAtom, [nIneqQuad])
    end if
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
        multipoleInp = multipoleOut
      else
        call mix(pChrgMixer, qInpRed, qDiffRed)
      #:if WITH_MPI
        ! Synchronise charges in order to avoid mixers that store a history drifting apart
        call mpifx_allreduceip(env%mpi%globalComm, qInpRed, MPI_SUM)
        qInpRed(:) = qInpRed / env%mpi%globalComm%size
      #:endif
        call expandCharges(qInpRed, orb, nIneqOrb, iEqOrbitals, qInput, dftbU, qBlockIn,&
            & iEqBlockDftbu, species0, qiBlockIn, iEqBlockDftbuLS, iEqBlockOnSite, iEqBlockOnSiteLS)
        if (nIneqDip > 0) then
          ! FIXME: Assumes we always mix all dipole moments
          nMix = nIneqOrb
          multipoleInp%dipoleAtom(:, :, :) = reshape(qInpRed(nMix+1:nMix+nIneqDip), &
              & shape(multipoleInp%dipoleAtom))
        end if
        if (nIneqQuad > 0) then
          ! FIXME: Assumes we always mix all quadrupole moments
          nMix = nIneqOrb + nIneqDip
          multipoleInp%quadrupoleAtom(:, :, :) = reshape(qInpRed(nMix+1:nMix+nIneqQuad),&
              & shape(multipoleInp%quadrupoleAtom))
        end if
      end if
    end if

  end subroutine getNextInputCharges


  !> Update delta density matrix rather than merely q for rangeseparation
  subroutine getNextInputDensity(SSqrReal, ints, neighbourList, nNeighbourSK, iAtomStart,&
      & iSparseStart, img2CentCell, pChrgMixer, qOutput, orb, tHelical, species, coord, iGeoStep,&
      & iSccIter, minSccIter, maxSccIter, sccTol, tStopScc, tReadChrg, q0, qInput, sccErrorQ,&
      & tConverged, deltaRhoOut, deltaRhoIn, deltaRhoDiff, qBlockIn, qBlockOut)

    !> Square dense overlap storage
    real(dp), allocatable, intent(inout) :: SSqrReal(:,:)

    !> Integral container
    type(TIntegral), intent(in) :: ints

    !> list of neighbours for each atom
    type(TNeighbourList), intent(in) :: neighbourList

    !> Number of neighbours for each of the atoms
    integer, intent(in) :: nNeighbourSK(:)

    !> Start of atomic blocks in dense arrays
    integer, allocatable, intent(in) :: iAtomStart(:)

    !> Index array for the start of atomic blocks in sparse arrays
    integer, intent(in) :: iSparseStart(:,:)

    !> map from image atoms to the original unique atom
    integer, intent(in) :: img2CentCell(:)

    !> Charge mixing object
    type(TMixer), intent(inout) :: pChrgMixer

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

    !> minimum number of SCC iterations to perform
    integer, intent(in) :: minSccIter

    !> maximum number of SCC iterations before terminating loop
    integer, intent(in) :: maxSccIter

    !> Tolerance on SCC charges between input and output
    real(dp), intent(in) :: sccTol

    !> Should the SCC loop stop
    logical, intent(in) :: tStopScc

    !> Were initial charges read from disc?
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

    !> delta density matrix as input for next SCC cycle
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
          call unpackHelicalHS(SSqrReal, ints%overlap, neighbourList%iNeighbour, nNeighbourSK,&
              & iAtomStart, iSparseStart, img2CentCell, orb, species, coord)
        else
          call unpackHS(SSqrReal, ints%overlap, neighbourList%iNeighbour, nNeighbourSK, iAtomStart,&
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
  subroutine reduceCharges(orb, nIneqOrb, iEqOrbitals, qOrb, qRed, qBlock, dftbU, iEqBlockDftbu,&
      & qiBlock, iEqBlockDftbuLS, iEqBlockOnSite, iEqBlockOnSiteLS)

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

    !> Are there orbital potentials present
    type(TDftbU), intent(in), allocatable :: dftbU

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
          call onsBlock_reduce(qiBlock, iEqBlockOnSiteLS, orb, qRed, isSkew=.true.)
        end if
      else
        ! only a subset of blocks are covered in +U type operations
        call appendBlockReduced(qBlockUpDown, iEqBlockDFTBU, orb, qRed)
        if (allocated(qiBlock)) then
          call appendBlockReduced(qiBlock, iEqBlockDFTBULS, orb, qRed, isSkew=.true.)
        end if
      end if
    end if

  end subroutine reduceCharges


  !> Expand reduced charges according orbital equivalency rules.
  subroutine expandCharges(qRed, orb, nIneqOrb, iEqOrbitals, qOrb, dftbU, qBlock, iEqBlockDftbu,&
      & species0, qiBlock, iEqBlockDftbuLS, iEqBlockOnSite, iEqBlockOnSiteLS)

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

    !> Are there orbital potentials present
    type(TDftbU), intent(in), allocatable :: dftbU

    !> Block (dual) populations, if also stored in reduced form
    real(dp), intent(inout), allocatable :: qBlock(:,:,:,:)

    !> equivalences for block charges
    integer, intent(in), allocatable :: iEqBlockDftbU(:,:,:,:)

    !> species of central cell atoms
    integer, intent(in), allocatable :: species0(:)

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
          call Onsblock_expand(qRed, iEqBlockOnSiteLS, orb, qiBlock, isSkew=.true.)
        end if
      else
        ! only a subset of blocks are covered in +U type operations
        call dftbU%expandBlock(qRed, iEqBlockDftbu, orb, qBlock, species0, orbEquiv=iEqOrbitals)
        if (allocated(qiBlock)) then
          call dftbU%expandBlock(qRed, iEqBlockDftbuLS, orb, qiBlock, species0, isSkew=.true.)
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
      & isGeoOpt, tDerivs, tConverged, tReadChrg, tStopScc) result(tRestart)

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
    logical, intent(in) :: isGeoOpt

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
    tRestart = (restartFreq > 0 .and. .not. (tMD .or. isGeoOpt .or. tDerivs) .and. maxSccIter > 1)
    if (tRestart) then

      ! Do we have enough iterations already?
      tEnoughIters = (iSccIter >= minSccIter .or. tReadChrg .or. iGeoStep > 0)

      ! Is current iteration the right one for writing a restart file?
      tRestartIter = (iSccIter == maxSccIter .or. tStopScc .or. mod(iSccIter, restartFreq) == 0)

      tRestart = (tConverged .or. (tEnoughIters .and. tRestartIter))
    end if

  end function needsSccRestartWriting

  !> Do the linear response excitation calculation.
  subroutine calculateLinRespExcitations(env, linearResponse, parallelKS, sccCalc, qOutput, q0,&
      & ints, eigvecsReal, eigen, filling, coord, species, speciesName, orb, skHamCont,&
      & skOverCont, autotestTag, taggedWriter, runId, neighbourList, nNeighbourSk, denseDesc,&
      & iSparseStart, img2CentCell, tWriteAutotest, tForces, tLinRespZVect, tPrintExcEigvecs,&
      & tPrintExcEigvecsTxt, nonSccDeriv, dftbEnergy, energies, work, rhoSqrReal, deltaRhoOutSqr,&
      & excitedDerivs, dQAtomEx, occNatural, rangeSep)

    !> Environment settings
    type(TEnvironment), intent(in) :: env

    !> excited state settings
    type(TLinResp), intent(inout), allocatable :: linearResponse

    !> K-points and spins to process
    type(TParallelKS), intent(in) :: parallelKS

    !> SCC module internal variables
    type(TScc), intent(in) :: sccCalc

    !> electrons in atomic orbitals
    real(dp), intent(in) :: qOutput(:,:,:)

    !> reference atomic orbital occupations
    real(dp), intent(in) :: q0(:,:,:)

    !> Integral container
    type(TIntegral), intent(in) :: ints

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

    !> non-SCC hamiltonian information
    type(TSlakoCont), intent(in) :: skHamCont

    !> overlap information
    type(TSlakoCont), intent(in) :: skOverCont

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
    type(TNonSccDiff), intent(in) :: nonSccDeriv

    !> Energy contributions and total
    type(TEnergies), intent(inout) :: dftbEnergy

    !> energies of all solved states
    real(dp), intent(inout), allocatable :: energies(:)

    !> Working array of the size of the dense matrices.
    real(dp), intent(out) :: work(:,:)

    !> density matrix in dense form
    real(dp), intent(inout), allocatable :: rhoSqrReal(:,:,:)

    !> difference density matrix (vs. uncharged atoms) in dense form
    real(dp), pointer, intent(inout) :: deltaRhoOutSqr(:,:,:)

    !> excited state energy derivative with respect to atomic coordinates
    real(dp), intent(inout), allocatable :: excitedDerivs(:,:)

    !> Mulliken charges in excited state
    real(dp), intent(out) :: dQAtomEx(:)

    !> natural orbital occupation numbers
    real(dp), intent(inout), allocatable :: occNatural(:)

    !> Data for range-separated calculation
    type(TRangeSepFunc), allocatable, intent(inout) :: rangeSep

    real(dp), allocatable :: dQAtom(:,:)
    real(dp), allocatable :: naturalOrbs(:,:,:)
    integer, pointer :: pSpecies0(:)
    integer :: iSpin, nSpin, nAtom, fdAutotest
    logical :: tSpin

    nAtom = size(qOutput, dim=2)
    nSpin = size(eigen, dim=2)
    tSpin = (nSpin == 2)
    pSpecies0 => species(1:nAtom)

    dftbEnergy%Eexcited = 0.0_dp
    allocate(dQAtom(nAtom, nSpin))
    dQAtom(:,:) = sum(qOutput(:,:,:) - q0(:,:,:), dim=1)
    call unpackHS(work, ints%overlap, neighbourList%iNeighbour, nNeighbourSK, denseDesc%iAtomStart,&
        & iSparseStart, img2CentCell)
    call blockSymmetrizeHS(work, denseDesc%iAtomStart)
    if (allocated(rhoSqrReal)) then
      do iSpin = 1, nSpin
        call blockSymmetrizeHS(rhoSqrReal(:,:,iSpin), denseDesc%iAtomStart)
      end do
    end if
    if (tForces .and. allocated(rangeSep)) then
      do iSpin = 1, nSpin
        call blockSymmetrizeHS(deltaRhoOutSqr(:,:,iSpin), denseDesc%iAtomStart)
      end do
    end if
    if (tWriteAutotest) then
      open(newUnit=fdAutotest, file=autotestTag, position="append")
    end if

    if (tLinRespZVect) then
      if (tPrintExcEigVecs) then
        allocate(naturalOrbs(orb%nOrb, orb%nOrb, 1))
      end if
      call LinResp_addGradients(tSpin, linearResponse, denseDesc%iAtomStart, eigvecsReal, eigen,&
          & work, filling, coord(:,:nAtom), sccCalc, dQAtom, pSpecies0, neighbourList%iNeighbour,&
          & img2CentCell, orb, skHamCont, skOverCont, tWriteAutotest, fdAutotest, taggedWriter,&
          & rangeSep, dftbEnergy%Eexcited, energies, excitedDerivs, nonSccDeriv,&
          & rhoSqrReal, deltaRhoOutSqr, occNatural, naturalOrbs)
      if (tPrintExcEigvecs) then
        call writeRealEigvecs(env, runId, neighbourList, nNeighbourSK, denseDesc, iSparseStart,&
            & img2CentCell, pSpecies0, speciesName, orb, ints%overlap, parallelKS, &
            & tPrintExcEigvecsTxt, naturalOrbs, work, fileName="excitedOrbs")
      end if
    else
      call linResp_calcExcitations(linearResponse, tSpin, denseDesc, eigvecsReal, eigen, work,&
          & filling, coord(:,:nAtom), sccCalc, dQAtom, pSpecies0, neighbourList%iNeighbour,&
          & img2CentCell, orb, tWriteAutotest, fdAutotest, taggedWriter,&
          & rangeSep, dftbEnergy%Eexcited, energies)
    end if
    dftbEnergy%Etotal = dftbEnergy%Etotal + dftbEnergy%Eexcited
    dftbEnergy%EMermin = dftbEnergy%EMermin + dftbEnergy%Eexcited
    dftbEnergy%EGibbs = dftbEnergy%EGibbs + dftbEnergy%Eexcited
    if (tWriteAutotest) then
      close(fdAutotest)
    end if

  end subroutine calculateLinRespExcitations

  !> Get the XLBOMD charges for the current geometry.
  subroutine getXlbomdCharges(xlbomdIntegrator, qOutRed, pChrgMixer, orb, nIneqOrb, iEqOrbitals,&
      & qInput, qInpRed, dftbU, iEqBlockDftbu, qBlockIn, species0, iEqBlockDftbuLS, qiBlockIn,&
      & iEqBlockOnSite, iEqBlockOnSiteLS)

    !> integrator for the extended Lagrangian
    type(TXLBOMD), intent(inout) :: xlbomdIntegrator

    !> output charges, reduced by equivalences
    real(dp), intent(in) :: qOutRed(:)

    !> SCC mixer
    type(TMixer), intent(inout) :: pChrgMixer

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

    !> Are there orbital potentials present
    type(TDftbU), intent(in), allocatable :: dftbU

    !> +U equivalences
    integer, intent(in), allocatable :: iEqBlockDftbU(:,:,:,:)

    !> central cell species
    integer, intent(in), allocatable :: species0(:)

    !> block input charges
    real(dp), intent(inout), allocatable :: qBlockIn(:,:,:,:)

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
    call expandCharges(qInpRed, orb, nIneqOrb, iEqOrbitals, qInput, dftbU, qBlockIn, iEqBlockDftbu,&
        & species0, qiBlockIn, iEqBlockDftbuLS, iEqBlockOnSite, iEqBlockOnSiteLS)

  end subroutine getXlbomdCharges


  !> Calculates dipole moment.
  subroutine getDipoleMoment(qOutput, q0, dipAtom, coord, dipoleMoment, iAtInCentralRegion)

    !> electrons in orbitals
    real(dp), intent(in) :: qOutput(:,:,:)

    !> reference atomic charges
    real(dp), intent(in) :: q0(:,:,:)

    !> Dipole populations for each atom
    real(dp), intent(in), optional :: dipAtom(:,:,:)

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

    if (present(dipAtom)) then
      do ii = 1, size(iAtInCentralRegion)
        iAtom = iAtInCentralRegion(ii)
        dipoleMoment(:) = dipoleMoment(:) - dipAtom(:, iAtom, 1)
      end do
    end if

  end subroutine getDipoleMoment


  !> Prints dipole moment calculated by the derivative of H with respect to the external field.
  subroutine checkDipoleViaHellmannFeynman(rhoPrim, q0, coord0, ints, orb, neighbourList,&
      & nNeighbourSK, species, iSparseStart, img2CentCell, eFieldScaling)

    !> Density matrix in sparse storage
    real(dp), intent(in) :: rhoPrim(:,:)

    !> Reference orbital charges
    real(dp), intent(in) :: q0(:,:,:)

    !> Central cell atomic coordinates
    real(dp), intent(in) :: coord0(:,:)

    !> Integral container
    type(TIntegral), intent(in) :: ints

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

    !> Instance of electric/dipole scaling due to any dielectric media effects
    class(TScaleExtEField), intent(in) :: eFieldScaling

    real(dp), allocatable :: hprime(:,:), dipole(:,:), potentialDerivative(:,:)
    integer :: nAtom, sparseSize, iAt, iCart

    sparseSize = size(ints%overlap)
    nAtom = size(q0, dim=2)
    allocate(hprime(sparseSize, 1))
    allocate(dipole(size(q0, dim=1), nAtom))
    allocate(potentialDerivative(nAtom, 1))
    write(stdOut,*)
    write(stdOut, "(A)", advance='no') 'Hellmann Feynman dipole:'

    ! loop over directions
    do iCart = 1, 3
      potentialDerivative(:,:) = 0.0_dp
      ! Potential from dH/dE
      potentialDerivative(:,1) = -eFieldScaling%scaledExtEField(coord0(iCart,:))
      hprime(:,:) = 0.0_dp
      dipole(:,:) = 0.0_dp
      call addShift(hprime, ints%overlap, nNeighbourSK, neighbourList%iNeighbour, species, orb,&
          & iSparseStart, nAtom, img2CentCell, potentialDerivative)

      ! evaluate <psi| dH/dE | psi> = Tr_part rho dH/dE
      call mulliken(dipole, hprime(:,1), rhoPrim(:,1), orb, neighbourList%iNeighbour, nNeighbourSK,&
          & img2CentCell, iSparseStart)

      ! add nuclei term for derivative wrt E
      do iAt = 1, nAtom
        dipole(1, iAt) = dipole(1, iAt) + sum(q0(:, iAt, 1)) *&
            & eFieldScaling%scaledExtEField(coord0(iCart, iAt))
      end do
      write(stdOut, "(F16.8)", advance='no') sum(dipole)
    end do
    write(stdOut, *) " au"

  end subroutine checkDipoleViaHellmannFeynman


  !> Calculate the energy weighted density matrix
  !>
  !> NOTE: Dense eigenvector and overlap matrices are overwritten.
  !>
  subroutine getEnergyWeightedDensity(env, negfInt, electronicSolver, denseDesc, forceType,&
      & filling, eigen, kPoint, kWeight, neighbourList, nNeighbourSK, orb, iSparseStart,&
      & img2CentCell, iCellVEc, cellVec, tRealHS, ints, parallelKS, tHelical, species, coord,&
      & iSCC, mu, ERhoPrim, HSqrReal, SSqrReal, HSqrCplx, SSqrCplx, errStatus)

    !> Environment settings
    type(TEnvironment), intent(inout) :: env

    !> NEGF interface
    type(TNegfInt), intent(inout) :: negfInt

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

    !> Is the hamiltonian real (no k-points/molecule/gamma point)?
    logical, intent(in) :: tRealHS

    !> Integral container
    type(TIntegral), intent(in) :: ints

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

    !> Storage for dense hamiltonian matrix (complex case)
    complex(dp), intent(inout), allocatable :: HSqrCplx(:,:,:)

    !> Storage for dense overlap matrix (complex case)
    complex(dp), intent(inout), allocatable :: SSqrCplx(:,:)

    !> Status of operation
    type(TStatus), intent(out) :: errStatus

    integer :: nSpin

    nSpin = size(ints%hamiltonian, dim=2)

    call env%globalTimer%startTimer(globalTimers%energyDensityMatrix)

    select case (electronicSolver%iSolver)

    case (electronicSolverTypes%GF)

      if (forceType /= forceTypes%orig) then
        @:RAISE_ERROR(errStatus, -1, "Alternative force evaluation methods are not supported by&
            & this electronic solver.")
      end if

    #:if WITH_TRANSPORT
      if (electronicSolver%iSolver == electronicSolverTypes%GF) then
        call negfInt%calcEdensity_green(iSCC, env, parallelKS%localKS, ints%hamiltonian,&
            & ints%overlap, neighbourlist%iNeighbour, nNeighbourSK, denseDesc%iAtomStart,&
            & iSparseStart, img2CentCell, iCellVec, cellVec, orb, kPoint, kWeight, mu, ERhoPrim)
      end if
    #:endif

    case (electronicSolverTypes%onlyTransport)

      @:RAISE_ERROR(errStatus, -1, "The OnlyTransport solver cannot calculate the energy weighted&
          & density matrix")

    case (electronicSolverTypes%qr, electronicSolverTypes%divideandconquer,&
        & electronicSolverTypes%relativelyrobust, electronicSolverTypes%elpa, &
        & electronicSolverTypes%magma_gvd)

      call getEDensityMtxFromEigvecs(env, denseDesc, forceType, filling, eigen, kPoint, kWeight,&
          & neighbourList, nNeighbourSK, orb, iSparseStart, img2CentCell, iCellVec, cellVec,&
          & tRealHS, ints, parallelKS, tHelical, species, coord, ERhoPrim, HSqrReal, SSqrReal,&
          & HSqrCplx, SSqrCplx, errStatus)
      @:PROPAGATE_ERROR(errStatus)

    case (electronicSolverTypes%omm, electronicSolverTypes%pexsi, electronicSolverTypes%ntpoly,&
        &electronicSolverTypes%elpadm)

      if (forceType /= forceTypes%orig) then
        @:RAISE_ERROR(errStatus, -1, "Alternative force evaluation methods are not supported by&
            & this electronic solver.")
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
      & ints, parallelKS, tHelical, species, coord, ERhoPrim, HSqrReal, SSqrReal, HSqrCplx,&
      & SSqrCplx, errStatus)

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

    !> Is the hamiltonian real (no k-points/molecule/gamma point)?
    logical, intent(in) :: tRealHS

    !> Integral container
    type(TIntegral), intent(in) :: ints

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

    !> Storage for dense hamiltonian matrix (complex case)
    complex(dp), intent(inout), allocatable :: HSqrCplx(:,:,:)

    !> Storage for dense overlap matrix (complex case)
    complex(dp), intent(inout), allocatable :: SSqrCplx(:,:)

    !> Status of operation
    type(TStatus), intent(out) :: errStatus

    integer :: nSpin

    nSpin = size(ints%hamiltonian, dim=2)

    if (nSpin == 4) then
      call getEDensityMtxFromPauliEigvecs(env, denseDesc, forceType, filling, eigen, kPoint,&
          & kWeight, neighbourList, nNeighbourSK, orb, iSparseStart, img2CentCell, iCellVec,&
          & cellVec, tRealHS, parallelKS, HSqrCplx, SSqrCplx, ERhoPrim, errStatus)
    else
      if (tRealHS) then
        call getEDensityMtxFromRealEigvecs(env, denseDesc, forceType, filling, eigen,&
            & neighbourList, nNeighbourSK, orb, iSparseStart, img2CentCell, ints,&
            & parallelKS, tHelical, species, coord, HSqrReal, SSqrReal, ERhoPrim, errStatus)

      else
        call getEDensityMtxFromComplexEigvecs(env, denseDesc, forceType, filling, eigen, kPoint,&
            & kWeight, neighbourList, nNeighbourSK, orb, iSparseStart, img2CentCell, iCellVec,&
            & cellVec, ints, parallelKS, tHelical, species, coord, HSqrCplx, SSqrCplx, ERhoPrim,&
            & errStatus)
      end if
    end if
    @:PROPAGATE_ERROR(errStatus)

  end subroutine getEDensityMtxFromEigvecs


  !> Calculates density matrix from real eigenvectors.
  subroutine getEDensityMtxFromRealEigvecs(env, denseDesc, forceType, filling, eigen,&
      & neighbourList, nNeighbourSK, orb, iSparseStart, img2CentCell, ints, parallelKS,&
      & tHelical, species, coord, eigvecsReal, work, ERhoPrim, errStatus)

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

    !> Integral container
    type(TIntegral), intent(in) :: ints

    !> K-points and spins to process
    type(TParallelKS), intent(in) :: parallelKS

    !> Is the geometry helical
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

    !> Status of operation
    type(TStatus), intent(out) :: errStatus

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
        call makeDensityMatrix(work, eigvecsReal(:,:,iKS), filling(:,1,iS), eigen(:,1,iS))
      #:endif

      case(forceTypes%dynamicT0)
        ! Correct force for XLBOMD for T=0K (DHD)
      #:if WITH_SCALAPACK
        if (tHelical) then
          call unpackHSHelicalRealBlacs(env%blacs, ints%hamiltonian(:,iS),&
              & neighbourList%iNeighbour, nNeighbourSK, iSparseStart, img2CentCell, orb, species,&
              & coord, denseDesc, work)
        else
          call unpackHSRealBlacs(env%blacs, ints%hamiltonian(:,iS), neighbourList%iNeighbour,&
              & nNeighbourSK, iSparseStart, img2CentCell, denseDesc, work)
        end if
        call makeDensityMtxRealBlacs(env%blacs%orbitalGrid, denseDesc%blacsOrbSqr, filling(:,1,iS),&
            & eigVecsReal(:,:,iKS), work2)
        call pblasfx_psymm(work2, denseDesc%blacsOrbSqr, work, denseDesc%blacsOrbSqr,&
            & eigvecsReal(:,:,iKS), denseDesc%blacsOrbSqr, side="L")
        call pblasfx_psymm(work2, denseDesc%blacsOrbSqr, eigvecsReal(:,:,iKS),&
            & denseDesc%blacsOrbSqr, work, denseDesc%blacsOrbSqr, side="R", alpha=0.5_dp)
      #:else
        if (tHelical) then
          call unpackHelicalHS(work, ints%hamiltonian(:,iS), neighbourlist%iNeighbour,&
              & nNeighbourSK, denseDesc%iAtomStart, iSparseStart, img2CentCell, orb, species, coord)
        else
          call unpackHS(work, ints%hamiltonian(:,iS), neighbourlist%iNeighbour, nNeighbourSK,&
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
          call unpackHSHelicalRealBlacs(env%blacs, ints%hamiltonian(:,iS),&
              & neighbourlist%iNeighbour, nNeighbourSK, iSparseStart, img2CentCell, orb, species,&
              & coord, denseDesc, work2)
        else
          call unpackHSRealBlacs(env%blacs, ints%hamiltonian(:,iS), neighbourlist%iNeighbour,&
              & nNeighbourSK, iSparseStart, img2CentCell, denseDesc, work2)
        end if
        call pblasfx_psymm(work, denseDesc%blacsOrbSqr, work2, denseDesc%blacsOrbSqr,&
            & eigvecsReal(:,:,iKS), denseDesc%blacsOrbSqr, side="L")
        if (tHelical) then
          call unpackHSHelicalRealBlacs(env%blacs, ints%overlap, neighbourlist%iNeighbour,&
              & nNeighbourSK, iSparseStart, img2CentCell, orb, species, coord, denseDesc, work)
        else
          call unpackHSRealBlacs(env%blacs, ints%overlap, neighbourlist%iNeighbour, nNeighbourSK,&
              & iSparseStart, img2CentCell, denseDesc, work)
        end if
        call psymmatinv(denseDesc%blacsOrbSqr, work, errStatus)
        @:PROPAGATE_ERROR(errStatus)
        call pblasfx_psymm(work, denseDesc%blacsOrbSqr, eigvecsReal(:,:,iKS),&
            & denseDesc%blacsOrbSqr, work2, denseDesc%blacsOrbSqr, side="R", alpha=0.5_dp)
        work(:,:) = work2
        call pblasfx_ptran(work2, denseDesc%blacsOrbSqr, work, denseDesc%blacsOrbSqr, alpha=1.0_dp,&
            & beta=1.0_dp)
      #:else
        call makeDensityMatrix(work, eigvecsReal(:,:,iKS), filling(:,1,iS))
        if (tHelical) then
          call unpackHelicalHS(work2, ints%hamiltonian(:,iS), neighbourlist%iNeighbour,&
              & nNeighbourSK, denseDesc%iAtomStart, iSparseStart, img2CentCell, orb, species, coord)
        else
          call unpackHS(work2, ints%hamiltonian(:,iS), neighbourlist%iNeighbour, nNeighbourSK,&
              & denseDesc%iAtomStart, iSparseStart, img2CentCell)
        end if
        call blocksymmetrizeHS(work2, denseDesc%iAtomStart)
        call symm(eigvecsReal(:,:,iKS), "L", work, work2)
        if (tHelical) then
          call unpackHelicalHS(work, ints%overlap, neighbourlist%iNeighbour, nNeighbourSK,&
              & denseDesc%iAtomStart, iSparseStart, img2CentCell, orb, species, coord)
        else
          call unpackHS(work, ints%overlap, neighbourlist%iNeighbour, nNeighbourSK,&
              & denseDesc%iAtomStart, iSparseStart, img2CentCell)
        end if
        call symmatinv(work, errStatus)
        @:PROPAGATE_ERROR(errStatus)
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
        call packHelicalHS(ERhoPrim, work, neighbourList%iNeighbour, nNeighbourSK,&
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
      & ints, parallelKS, tHelical, species, coord, eigvecsCplx, work, ERhoPrim, errStatus)

    !> Environment settings
    type(TEnvironment), intent(in) :: env

    !> Dense matrix descriptor
    type(TDenseDescr), intent(in) :: denseDesc

    !> Force evaluation method
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

    !> Integral container
    type(TIntegral), intent(in) :: ints

    !> K-points and spins to process
    type(TParallelKS), intent(in) :: parallelKS

    !> Is the geometry helical
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

    !> Status of operation
    type(TStatus), intent(out) :: errStatus

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
        call makeDensityMatrix(work, eigvecsCplx(:,:,iKS), filling(:,iK,iS), eigen(:,iK, iS))
      #:endif

      case(forceTypes%dynamicT0)
        ! Correct force for XLBOMD for T=0K (DHD)
      #:if WITH_SCALAPACK
        if (tHelical) then
          call unpackHSHelicalCplxBlacs(env%blacs, ints%hamiltonian(:,iS), kPoint(:,iK),&
              & neighbourList%iNeighbour, nNeighbourSK, iCellVec, cellVec, iSparseStart,&
              & img2CentCell, orb, species, coord, denseDesc, work)
        else
          call unpackHSCplxBlacs(env%blacs, ints%hamiltonian(:,iS), kPoint(:,iK),&
              & neighbourList%iNeighbour, nNeighbourSK, iCellVec, cellVec, iSparseStart,&
              & img2CentCell, denseDesc, work)
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
          call unpackHelicalHS(work, ints%hamiltonian(:,iS), kPoint(:,iK),&
              & neighbourlist%iNeighbour, nNeighbourSK, iCellVec, cellVec, denseDesc%iAtomStart,&
              & iSparseStart, img2CentCell, orb, species, coord)
        else
          call unpackHS(work, ints%hamiltonian(:,iS), kPoint(:,iK), neighbourlist%iNeighbour,&
              & nNeighbourSK, iCellVec, cellVec, denseDesc%iAtomStart, iSparseStart, img2CentCell)
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
          call unpackHSHelicalCplxBlacs(env%blacs, ints%hamiltonian(:,iS), kPoint(:,iK),&
              & neighbourlist%iNeighbour, nNeighbourSK, iCellVec, cellVec, iSparseStart,&
              & img2CentCell, orb, species, coord, denseDesc, work2)
        else
          call unpackHSCplxBlacs(env%blacs, ints%hamiltonian(:,iS), kPoint(:,iK),&
              & neighbourlist%iNeighbour, nNeighbourSK, iCellVec, cellVec, iSparseStart,&
              & img2CentCell, denseDesc, work2)
        end if
        call pblasfx_phemm(work, denseDesc%blacsOrbSqr, work2, denseDesc%blacsOrbSqr,&
            & eigvecsCplx(:,:,iKS), denseDesc%blacsOrbSqr, side="L")
        if (tHelical) then
          call unpackHSHelicalCplxBlacs(env%blacs, ints%overlap, kPoint(:,iK),&
              & neighbourlist%iNeighbour, nNeighbourSK, iCellVec, cellVec, iSparseStart,&
              & img2CentCell, orb, species, coord, denseDesc, work)
        else
          call unpackHSCplxBlacs(env%blacs, ints%overlap, kPoint(:,iK), neighbourlist%iNeighbour,&
              & nNeighbourSK, iCellVec, cellVec, iSparseStart, img2CentCell, denseDesc, work)
        end if
        call phermatinv(denseDesc%blacsOrbSqr, work, errStatus)
        @:PROPAGATE_ERROR(errStatus)
        call pblasfx_phemm(work, denseDesc%blacsOrbSqr, eigvecsCplx(:,:,iKS),&
            & denseDesc%blacsOrbSqr, work2, denseDesc%blacsOrbSqr, side="R", alpha=(0.5_dp, 0.0_dp))
        work(:,:) = work2
        call pblasfx_ptranc(work2, denseDesc%blacsOrbSqr, work, denseDesc%blacsOrbSqr,&
            & alpha=(1.0_dp, 0.0_dp), beta=(1.0_dp, 0.0_dp))
      #:else
        call makeDensityMatrix(work, eigvecsCplx(:,:,iKS), filling(:,iK,iS))
        if (tHelical) then
          call unpackHelicalHS(work2, ints%hamiltonian(:,iS), kPoint(:,iK),&
              & neighbourlist%iNeighbour, nNeighbourSK, iCellVec, cellVec, denseDesc%iAtomStart,&
              & iSparseStart, img2CentCell, orb, species, coord)
        else
          call unpackHS(work2, ints%hamiltonian(:,iS), kPoint(:,iK), neighbourlist%iNeighbour,&
              & nNeighbourSK, iCellVec, cellVec, denseDesc%iAtomStart, iSparseStart, img2CentCell)
        end if
        call blockHermitianHS(work2, denseDesc%iAtomStart)
        call hemm(eigvecsCplx(:,:,iKS), "L", work, work2)
        if  (tHelical) then
          call unpackHelicalHS(work, ints%overlap, kPoint(:,iK), neighbourlist%iNeighbour,&
              & nNeighbourSK, iCellVec, cellVec, denseDesc%iAtomStart, iSparseStart, img2CentCell,&
              & orb, species, coord)
        else
          call unpackHS(work, ints%overlap, kPoint(:,iK), neighbourlist%iNeighbour, nNeighbourSK,&
              & iCellVec, cellVec, denseDesc%iAtomStart, iSparseStart, img2CentCell)
        end if
        call hermatinv(work, errStatus)
        @:PROPAGATE_ERROR(errStatus)
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
        call packHelicalHS(ERhoPrim, work, kPoint(:,iK), kWeight(iK), neighbourList%iNeighbour,&
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
  subroutine getEDensityMtxFromPauliEigvecs(env, denseDesc, forceType, filling, eigen, kPoint,&
      & kWeight, neighbourList, nNeighbourSK, orb, iSparseStart, img2CentCell, iCellVec, cellVec,&
      & tRealHS, parallelKS, eigvecsCplx, work, ERhoPrim, errStatus)

    !> Environment settings
    type(TEnvironment), intent(in) :: env

    !> Dense matrix descriptor
    type(TDenseDescr), intent(in) :: denseDesc

    !> Force evaluation method
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

    !> Is the hamiltonian real (no k-points/molecule/gamma point)?
    logical, intent(in) :: tRealHS

    !> K-points and spins to process
    type(TParallelKS), intent(in) :: parallelKS

    !> Eigenvectors
    complex(dp), intent(inout) :: eigvecsCplx(:,:,:)

    !> Work array
    complex(dp), intent(out) :: work(:,:)

    !> Sparse energy weighted density matrix
    real(dp), intent(out) :: ERhoPrim(:)

    !> Status of operation
    type(TStatus), intent(out) :: errStatus

    integer :: iKS, iK

    if (forceType /= forceTypes%orig) then
      @:RAISE_ERROR(errStatus, -1, "Alternative force evaluation methods are not supported for&
          & Pauli hamiltonians at the moment.")
    end if

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
  subroutine getGradients(env, boundaryConds, sccCalc, tblite, isExtField, isXlbomd, nonSccDeriv,&
      & rhoPrim, ERhoPrim, qOutput, q0, skHamCont, skOverCont, repulsive, neighbourList,&
      & nNeighbourSK, species, img2CentCell, iSparseStart, orb, potential, coord, derivs,&
      & groundDerivs, tripletderivs, mixedderivs, iRhoPrim, thirdOrd, solvation, qDepExtPot,&
      & chrgForces, dispersion, rangeSep, SSqrReal, ints, denseDesc, deltaRhoOutSqr,&
      & halogenXCorrection, tHelical, coord0, deltaDftb)

    !> Environment settings
    type(TEnvironment), intent(inout) :: env

    !> Boundary conditions on the geometry
    type(TBoundaryConditions), intent(in) :: boundaryConds

    !> SCC module internal variables
    type(TScc), allocatable, intent(inout) :: sccCalc

    !> Library communication interface
    type(TTBLite), allocatable, intent(inout) :: tblite

    !> external electric field
    logical, intent(in) :: isExtField

    !> extended Lagrangian active?
    logical, intent(in) :: isXlbomd

    !> method for calculating derivatives of S and H0
    type(TNonSccDiff), intent(in) :: nonSccDeriv

    !> sparse density matrix
    real(dp), intent(in) :: rhoPrim(:,:)

    !> energy  weighted density matrix
    real(dp), intent(in) :: ERhoPrim(:)

    !> electron populations (may be unallocated for non-scc case)
    real(dp), allocatable, intent(in) :: qOutput(:,:,:)

    !> reference atomic charges (may be unallocated for non-scc case)
    real(dp), allocatable, intent(in) :: q0(:,:,:)

    !> non-SCC hamiltonian information
    type(TSlakoCont), intent(in) :: skHamCont

    !> overlap information
    type(TSlakoCont), intent(in) :: skOverCont

    !> repulsive information
    class(TRepulsive), allocatable, intent(in) :: repulsive

    !> list of neighbours for each atom
    type(TNeighbourList), intent(in) :: neighbourList

    !> Number of neighbours for each of the atoms
    integer, intent(in) :: nNeighbourSK(:)

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

    !> derivatives of ground state energy wrt to atomic positions
    real(dp), intent(inout), allocatable :: groundDerivs(:,:)

    !> derivatives of triplet energy wrt to atomic positions (TI-DFTB excited states)
    real(dp), intent(inout), allocatable :: tripletDerivs(:,:)

    !> derivatives of mixed energy wrt to atomic positions (TI-DFTB excited states)
    real(dp), intent(inout), allocatable :: mixedDerivs(:,:)

    !> imaginary part of density matrix
    real(dp), intent(in), allocatable :: iRhoPrim(:,:)

    !> Is 3rd order SCC being used
    type(TThirdOrder), intent(inout), allocatable :: thirdOrd

    !> Solvation model
    class(TSolvation), allocatable, intent(inout) :: solvation

    !> Population dependant external potential
    type(TQDepExtPotProxy), intent(inout), allocatable :: qDepExtPot

    !> forces on external charges
    real(dp), intent(inout), allocatable :: chrgForces(:,:)

    !> dispersion interactions
    class(TDispersionIface), intent(inout), allocatable :: dispersion

    !> Data from rangeseparated calculations
    type(TRangeSepFunc), intent(inout), allocatable :: rangeSep

    !> dense overlap matrix, required for rangeSep
    real(dp), intent(inout), allocatable :: SSqrReal(:,:)

    !> Integral container
    type(TIntegral), intent(in) :: ints

    !> Dense matrix descriptor,required for rangeSep
    type(TDenseDescr), intent(in) :: denseDesc

    !> Change in density matrix during this SCC step for rangesep
    real(dp), pointer, intent(in) :: deltaRhoOutSqr(:,:,:)

    !> Correction for halogen bonds
    type(THalogenX), allocatable, intent(inout) :: halogenXCorrection

    !> Is the geometry helical
    logical, intent(in) :: tHelical

    !> Central cell atomic coordinates
    real(dp), intent(in) :: coord0(:,:)

    !> Determinant derived type
    type(TDftbDeterminants), intent(in) :: deltaDftb

    real(dp), allocatable :: tmpDerivs(:,:)
    real(dp), allocatable :: dQ(:,:,:)
    real(dp), allocatable :: dipoleAtom(:,:)
    logical :: tImHam, tExtChrg, tSccCalc
    integer :: nAtom, iAt

    tSccCalc = allocated(sccCalc)
    tImHam = allocated(iRhoPrim)
    tExtChrg = allocated(chrgForces)
    nAtom = size(derivs, dim=2)

    if (allocated(tblite)) then
      dipoleAtom = potential%dipoleAtom
      if (allocated(potential%extDipoleAtom)) then
        dipoleAtom(:, :) = dipoleAtom + potential%extDipoleAtom
      end if
    end if

    allocate(tmpDerivs(3, nAtom))
    derivs(:,:) = 0.0_dp

    if (.not. (tSccCalc .or. isExtField)) then
      ! No external or internal potentials
      if (allocated(tblite)) then
        call tblite%buildDerivativeShift(env, rhoPrim, ERhoPrim, coord, species, &
          & nNeighbourSK, neighbourList%iNeighbour, img2CentCell, iSparseStart, orb, &
          & potential%intBlock, dipoleAtom, potential%quadrupoleAtom)
        call tblite%addGradients(env, neighbourList, species, coord, img2centCell, derivs)
      else
        if (tImHam) then
          call derivative_shift(env, derivs, nonSccDeriv, rhoPrim, iRhoPrim, ERhoPrim, skHamCont,&
              & skOverCont, coord, species, neighbourList%iNeighbour, nNeighbourSK, img2CentCell,&
              & iSparseStart, orb, potential%intBlock, potential%iorbitalBlock)
        else
          call derivative_shift(env, derivs, nonSccDeriv, rhoPrim(:,1), ERhoPrim, skHamCont,&
              & skOverCont, coord, species, neighbourList%iNeighbour, nNeighbourSK, img2CentCell,&
              & iSparseStart, orb, tHelical)
        end if
      end if
    else
      if (allocated(tblite)) then
        call tblite%buildDerivativeShift(env, rhoPrim, ERhoPrim, coord, species, &
          & nNeighbourSK, neighbourList%iNeighbour, img2CentCell, iSparseStart, orb, &
          & potential%intBlock, dipoleAtom, potential%quadrupoleAtom)
        call tblite%addGradients(env, neighbourList, species, coord, img2centCell, derivs)
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
      end if

      if (tExtChrg) then
        chrgForces(:,:) = 0.0_dp
        if (isXlbomd) then
          call error("XLBOMD does not work with external charges yet!")
        else
          call sccCalc%addForceDc(env, derivs, species, neighbourList%iNeighbour, img2CentCell,&
              & chrgForces)
        end if
      else if (tSccCalc) then
        if (isXlbomd) then
          call sccCalc%addForceDcXlbomd(env, species, orb, neighbourList%iNeighbour,&
              & img2CentCell, qOutput, q0, derivs)
        else
          call sccCalc%addForceDc(env, derivs, species, neighbourList%iNeighbour, img2CentCell)
        end if
      end if

      if (allocated(thirdOrd)) then
        if (isXlbomd) then
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

      if (isExtField) then
        do iAt = 1, nAtom
          derivs(:, iAt) = derivs(:, iAt)&
              & + sum(qOutput(:, iAt, 1) - q0(:, iAt, 1)) * potential%extGrad(:, iAt)
        end do
      end if

    end if

    if (allocated(solvation)) then
      if (isXlbomd) then
        call error("XLBOMD does not work with solvation yet!")
      else
        call solvation%addGradients(env, neighbourList, species, coord, img2CentCell, derivs)
      end if
    end if

    if (allocated(dispersion)) then
      call dispersion%addGradients(env, neighbourList, img2CentCell, coord, species, derivs)
    end if

    if (allocated(halogenXCorrection)) then
      call halogenXCorrection%addGradients(derivs, coord, species, neighbourList, img2CentCell)
    end if

    if (allocated(rangeSep)) then
      if (tHelical) then
        call unpackHelicalHS(SSqrReal, ints%overlap, neighbourList%iNeighbour, nNeighbourSK,&
            & denseDesc%iAtomStart, iSparseStart, img2CentCell, orb, species, coord)
      else
        call unpackHS(SSqrReal, ints%overlap, neighbourList%iNeighbour, nNeighbourSK,&
            & denseDesc%iAtomStart, iSparseStart, img2CentCell)
      end if
      if (size(deltaRhoOutSqr, dim=3) > 2) then
        call error("Range separated forces do not support non-colinear spin")
      else
        call rangeSep%addLRGradients(derivs, nonSccDeriv, deltaRhoOutSqr, skOverCont, coord,&
            & species, orb, denseDesc%iAtomStart, SSqrReal, neighbourList%iNeighbour, nNeighbourSK)
      end if
    end if

    if (allocated(repulsive)) then
      call repulsive%getGradients(coord, species, img2CentCell, neighbourList, tmpDerivs)
    else
      tmpDerivs(:,:) = 0.0_dp
    end if

    derivs(:,:) = derivs + tmpDerivs

    call boundaryConds%alignVectorCentralCell(derivs, coord, coord0, nAtom)

    if(deltaDftb%isNonAufbau) then
      select case (deltaDftb%whichDeterminant(deltaDftb%iDeterminant))
      case (determinants%ground)
        groundDerivs(:,:) = derivs
      case (determinants%triplet)
        tripletDerivs(:,:) = derivs
      case (determinants%mixed)
        mixedDerivs(:,:) = derivs
      end select
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
  subroutine getStress(env, sccCalc, tblite, thirdOrd, isExtField, nonSccDeriv, rhoPrim,&
      & ERhoPrim, qOutput, q0, skHamCont, skOverCont, repulsive, neighbourList, nNeighbourSk,&
      & species, img2CentCell, iSparseStart, orb, potential, coord, latVec, invLatVec,&
      & cellVol, coord0, totalStress, totalLatDeriv, intPressure, iRhoPrim, solvation,&
      & dispersion, halogenXCorrection, deltaDftb, tripletStress, mixedStress)

    !> Environment settings
    type(TEnvironment), intent(in) :: env

    !> SCC module internal variables
    type(TScc), allocatable, intent(in) :: sccCalc

    !> Library interface handler
    type(TTBLite), allocatable, intent(inout) :: tblite

    !> Is 3rd order SCC being used
    type(TThirdOrder), intent(inout), allocatable :: thirdOrd

    !> External field
    logical, intent(in) :: isExtField

    !> method for calculating derivatives of S and H0
    type(TNonSccDiff), intent(in) :: nonSccDeriv

    !> density matrix
    real(dp), intent(in) :: rhoPrim(:,:)

    !> energy weighted density matrix
    real(dp), intent(in) :: ERhoPrim(:)

    !> electrons in orbitals
    real(dp), intent(in) :: qOutput(:,:,:)

    !> reference charges
    real(dp), intent(in) :: q0(:,:,:)

    !> non-SCC hamiltonian information
    type(TSlakoCont), intent(in) :: skHamCont

    !> overlap information
    type(TSlakoCont), intent(in) :: skOverCont

    !> repulsive information
    class(TRepulsive), allocatable, intent(in) :: repulsive

    !> list of neighbours for each atom
    type(TNeighbourList), intent(in) :: neighbourList

    !> Number of neighbours for each of the atoms
    integer, intent(in) :: nNeighbourSK(:)

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

    !> Solvation model
    class(TSolvation), allocatable, intent(inout) :: solvation

    !> dispersion interactions
    class(TDispersionIface), allocatable, intent(inout) :: dispersion

    !> Correction for halogen bonds
    type(THalogenX), allocatable, intent(inout) :: halogenXCorrection

    !> Determinant derived type
    type(TDftbDeterminants), intent(in) :: deltaDftb

    !> Stress tensor in triplet state (TI-DFTB excited states)
    real(dp), intent(inout), optional :: tripletStress(:,:)

    !> Stress tensor in mixed state (TI-DFTB excited states)
    real(dp), intent(inout), optional :: mixedStress(:,:)

    real(dp) :: tmpStress(3, 3)
    logical :: tImHam

    tImHam = allocated(iRhoPrim)
    totalStress(:, :) = 0.0_dp

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
      if (allocated(tblite)) then
        call tblite%getStress(tmpStress)
        totalStress(:, :) = totalStress + tmpStress
      else
        if (tImHam) then
          call getBlockiStress(env, totalStress, nonSccDeriv, rhoPrim, iRhoPrim, ERhoPrim,&
              & skHamCont, skOverCont, coord, species, neighbourList%iNeighbour, nNeighbourSK,&
              & img2CentCell, iSparseStart, orb, potential%intBlock, potential%iorbitalBlock,&
              & cellVol)
        else
          call getNonSCCStress(env, totalStress, nonSccDeriv, rhoPrim(:,1), ERhoPrim, skHamCont,&
              & skOverCont, coord, species, neighbourList%iNeighbour, nNeighbourSK, img2CentCell,&
              & iSparseStart, orb, cellVol)
        end if
      end if
    end if

    if (allocated(solvation)) then
      call solvation%getStress(tmpStress)
      totalStress(:,:) = totalStress + tmpStress
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

    if (isExtField) then
      call getExtFieldStress(latVec, cellVol, q0, qOutput, potential%extGrad, coord0, tmpStress)
      totalStress(:,:) = totalStress + tmpStress
    end if

    if (allocated(repulsive)) then
      call repulsive%getStress(coord, species, img2CentCell, neighbourList, cellVol, tmpStress)
    else
      tmpStress(:,:) = 0.0_dp
    end if
    totalStress(:,:) = totalStress + tmpStress

    if(deltaDftb%isNonAufbau) then
      select case (deltaDftb%whichDeterminant(deltaDftb%iDeterminant))
      case (determinants%triplet)
        tripletStress(:,:) = totalStress
      case (determinants%mixed)
        mixedStress(:,:) = totalStress
      end select
    end if

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


  !> Unfold constrained lattice vectors to full one.
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
  subroutine getNextDerivStep(derivDriver, derivs, indMovedAtoms, indDerivAtoms, coord, tGeomEnd)

    !> Driver for the finite difference second derivatives
    type(TNumDerivs), intent(inout) :: derivDriver

    !> first derivatives of energy at the current coordinates
    real(dp), intent(in) :: derivs(:,:)

    !> indices of moving atoms
    integer, intent(in) :: indMovedAtoms(:)

    !> indices of atoms for which 2nd derivatives should be calculated
    integer, intent(in) :: indDerivAtoms(:)

    !> atomic coordinates
    real(dp), intent(inout) :: coord(:,:)

    !> has the process terminated
    logical, intent(out) :: tGeomEnd

    real(dp) :: newCoords(3, size(indMovedAtoms))

    call next(derivDriver, newCoords, derivs(:, indDerivAtoms), tGeomEnd)
    coord(:, indMovedAtoms) = newCoords

  end subroutine getNextDerivStep


  !> Returns the coordinates for the next coordinate optimisation step.
  subroutine getNextCoordinateOptStep(pGeoCoordOpt, energy, derivss, indMovedAtom, coord0,&
      & diffGeo, tCoordEnd, tRemoveExcitation)

    !> optimiser for atomic coordinates
    type(TGeoOpt), intent(inout) :: pGeoCoordOpt

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
    type(TGeoOpt), intent(inout) :: pGeoLatOpt

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
    type(TMdIntegrator), intent(inout) :: pMdIntegrator

    !> Molecular dynamics reference frame information
    type(TMdCommon), intent(in) :: pMdFrame

    !> Temperature profile in MD
    type(TTempProfile), allocatable, intent(inout) :: temperatureProfile

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
      call temperatureProfile%next()
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
  subroutine calcPipekMezeyLocalisation(env, pipekMezey, tPrintEigvecsTxt, nEl, filling, ints,&
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

    !> Integral container
    type(TIntegral), intent(in) :: ints

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

    !> Storage for dense hamiltonian matrix (complex case)
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
        call unpackHelicalHS(SSqrReal, ints%overlap, neighbourList%iNeighbour, nNeighbourSK,&
            & denseDesc%iAtomStart, iSparseStart, img2CentCell, orb, species, coord)
      else
        call unpackHS(SSqrReal, ints%overlap, neighbourList%iNeighbour, nNeighbourSK,&
            & denseDesc%iAtomStart, iSparseStart, img2CentCell)
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
          & img2CentCell, species(:nAtom), speciesName, orb, ints%overlap, parallelKS,&
          & tPrintEigvecsTxt, eigvecsReal, SSqrReal, fileName="localOrbs")
    else

      localisation = 0.0_dp
      do iKS = 1, parallelKS%nLocalKS
        iK = parallelKS%localKS(1, iKS)
        iSpin = parallelKS%localKS(2, iKS)
        nFilledLev = nint(nEl(iSpin) / real( 3 - nSpin, dp))
        localisation = localisation + pipekMezey%getLocalisation(&
            & eigvecsCplx(:,:nFilledLev,iKS), SSqrCplx, ints%overlap, kpoint(:,iK), neighbourList,&
            & nNeighbourSK, iCellVec, cellVec, denseDesc%iAtomStart, iSparseStart, img2CentCell)
      end do
      write(stdOut, "(A, E20.12)") 'Original localisation', localisation

      ! actual localisation calls
      do iKS = 1, parallelKS%nLocalKS
        iK = parallelKS%localKS(1, iKS)
        iSpin = parallelKS%localKS(2, iKS)
        nFilledLev = nint(nEl(iSpin) / real( 3 - nSpin, dp))
        call pipekMezey%calcCoeffs(eigvecsCplx(:,:nFilledLev,iKS), SSqrCplx, ints%overlap,&
            & kpoint(:,iK), neighbourList, nNeighbourSK, iCellVec, cellVec, denseDesc%iAtomStart,&
            & iSparseStart, img2CentCell)
      end do

      localisation = 0.0_dp
      do iKS = 1, parallelKS%nLocalKS
        iK = parallelKS%localKS(1, iKS)
        iSpin = parallelKS%localKS(2, iKS)
        nFilledLev = nint(nEl(iSpin) / real( 3 - nSpin, dp))
        localisation = localisation + pipekMezey%getLocalisation(&
            & eigvecsCplx(:,:nFilledLev,iKS), SSqrCplx, ints%overlap, kpoint(:,iK), neighbourList,&
            & nNeighbourSK, iCellVec, cellVec, denseDesc%iAtomStart, iSparseStart, img2CentCell)
      end do
      write(stdOut, "(A, E20.12)") 'Final localisation', localisation

      call writeCplxEigvecs(env, runId, neighbourList, nNeighbourSK, cellVec, iCellVec, denseDesc,&
          & iSparseStart, img2CentCell, species, speciesName, orb, kPoint, ints%overlap,&
          & parallelKS, tPrintEigvecsTxt, eigvecsCplx, SSqrCplx, fileName="localOrbs")

    end if

  end subroutine calcPipekMezeyLocalisation


  !> Prints information about maximal forces in the system.
  subroutine printMaxForces(derivs, constrLatDerivs, tCoordOpt, tLatOpt, indMovedAtoms)

    !> Gradients on atoms ]3, nAtom]
    real(dp), intent(in), allocatable :: derivs(:,:)

    !> Lattice derivatives. Shape: [9]
    real(dp), intent(in) :: constrLatDerivs(:)

    !> Whether coordinate optimization is on.
    logical, intent(in) :: tCoordOpt

    !> Whether lattice optimization is on.
    logical, intent(in) :: tLatOpt

    !> Indices of moved atoms. Shape [nMovedAtoms].
    integer, intent(in) :: indMovedAtoms(:)

    real(dp) :: normedDeriv

    if (tCoordOpt) then
      call printMaxForce(maxval(abs(derivs(:, indMovedAtoms))))
      normedDeriv = norm2(derivs(:, indMovedAtoms)) / sqrt(real(size(indMovedAtoms), dp))
      call printForceNorm(normedDeriv)
    end if
    if (tLatOpt) then
      call printMaxLatticeForce(maxval(abs(constrLatDerivs)))
      normedDeriv = norm2(constrLatDerivs) / sqrt(3.0_dp)
      call printLatticeForceNorm(normedDeriv)
    end if

  end subroutine printMaxForces


#:if WITH_SOCKETS

  subroutine sendEnergyAndForces(env, socket, energy, derivs, totalStress, cellVol)

    !> Environment
    type(TEnvironment), intent(in) :: env

    !> Socket may be unallocated (as on follower processes)
    type(ipiSocketComm), allocatable, intent(inout) :: socket

    !> energy structure
    type(TEnergies), intent(in) :: energy

    !> energy derivatives
    real(dp), intent(in) :: derivs(:,:)

    !> stress tensor
    real(dp), intent(in) :: totalStress(:,:)

    !> cell volume
    real(dp), intent(in) :: cellVol

    if (env%tGlobalLead) then
      ! stress was computed above in the force evaluation block or is 0 if aperiodic
      call socket%send(energy%ETotal - sum(energy%TS), -derivs, totalStress * cellVol)
    end if
  end subroutine sendEnergyAndForces

#:endif


  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !!!! REKS subroutines
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Diagonalize H0 to obtain initial guess of eigenvectors
  !> or read eigenvectors in REKS
  !> Save dense overlap matrix elements
  !> Check Gamma point condition and set filling information
  subroutine getReksInitialSettings(env, denseDesc, h0, ints, neighbourList, nNeighbourSK,&
      & iSparseStart, img2CentCell, electronicSolver, iGeoStep, HSqrReal, SSqrReal, eigvecsReal,&
      & eigen, reks, errStatus)

    !> Environment settings
    type(TEnvironment), intent(inout) :: env

    !> Dense matrix descriptor
    type(TDenseDescr), intent(in) :: denseDesc

    !> hamiltonian in sparse storage
    real(dp), intent(in) :: h0(:)

    !> Integral container
    type(TIntegral), intent(in) :: ints

    !> list of neighbours for each atom
    type(TNeighbourList), intent(in) :: neighbourList

    !> Number of neighbours for each of the atoms
    integer, intent(in) :: nNeighbourSK(:)

    !> Index array for the start of atomic blocks in sparse arrays
    integer, intent(in) :: iSparseStart(:,:)

    !> map from image atoms to the original unique atom
    integer, intent(in) :: img2CentCell(:)

    !> Electronic solver information
    type(TElectronicSolver), intent(inout) :: electronicSolver

    !> Number of current geometry step
    integer, intent(in) :: iGeoStep

    !> dense hamiltonian matrix
    real(dp), intent(out) :: HSqrReal(:,:)

    !> dense overlap matrix
    real(dp), intent(out) :: SSqrReal(:,:)

    !> Eigenvectors on eixt
    real(dp), intent(inout) :: eigvecsReal(:,:,:)

    !> eigenvalues
    real(dp), intent(out) :: eigen(:,:,:)

    !> data type for REKS
    type(TReksCalc), intent(inout) :: reks

    !> Status of operation
    type(TStatus), intent(out) :: errStatus

    call env%globalTimer%startTimer(globalTimers%sparseToDense)
    call unpackHS(SSqrReal, ints%overlap, neighbourList%iNeighbour, nNeighbourSK, &
        & denseDesc%iAtomStart, iSparseStart, img2CentCell)
    call env%globalTimer%stopTimer(globalTimers%sparseToDense)

    reks%overSqr(:,:) = SSqrReal
    call blockSymmetrizeHS(reks%overSqr, denseDesc%iAtomStart)

    if (iGeoStep == 0) then

      if (.not. reks%tReadMO) then

        call env%globalTimer%startTimer(globalTimers%sparseToDense)
        call unpackHS(HSqrReal, h0, neighbourList%iNeighbour, nNeighbourSK, &
            & denseDesc%iAtomStart, iSparseStart, img2CentCell)
        call env%globalTimer%stopTimer(globalTimers%sparseToDense)

        eigen(:,:,:) = 0.0_dp
        call env%globalTimer%startTimer(globalTimers%diagonalization)
        call diagDenseMtx(env, electronicSolver, 'V', HSqrReal, SSqrReal, eigen(:,1,1), errStatus)
        @:PROPAGATE_ERROR(errStatus)
        call env%globalTimer%stopTimer(globalTimers%diagonalization)
        eigvecsReal(:,:,1) = HSqrReal

      else

        call readEigenvecs(eigvecsReal(:,:,1))
        call renormalizeEigenvecs(env, electronicSolver, eigvecsReal, reks)

      end if

      call constructMicrostates(reks)

    else

      call renormalizeEigenvecs(env, electronicSolver, eigvecsReal, reks)

    end if

    call checkGammaPoint(denseDesc, neighbourList%iNeighbour, &
        & nNeighbourSK, iSparseStart, img2CentCell, ints%overlap, reks)

  end subroutine getReksInitialSettings


  !> Normalize eigenvectors with unitary transformation
  subroutine renormalizeEigenvecs(env, electronicSolver, eigvecsReal, reks)
    use dftbp_math_blasroutines, only : gemm
    use dftbp_math_eigensolver, only : heev

    !> Environment settings
    type(TEnvironment), intent(inout) :: env

    !> Electronic solver information
    type(TElectronicSolver), intent(inout) :: electronicSolver

    !> Eigenvectors on eixt
    real(dp), intent(inout) :: eigvecsReal(:,:,:)

    !> data type for REKS
    type(TReksCalc), intent(inout) :: reks

    real(dp), allocatable :: tmpMat(:,:)
    real(dp), allocatable :: tmpS(:,:)
    real(dp), allocatable :: tmpC(:,:)
    real(dp), allocatable :: tmpEigen(:)
    real(dp), allocatable :: unitaryMat(:,:)

    integer :: nOrb, ii

    nOrb = size(reks%overSqr,dim=1)

    allocate(tmpMat(nOrb,nOrb))
    allocate(tmpS(nOrb,nOrb))
    allocate(tmpC(nOrb,nOrb))
    allocate(tmpEigen(nOrb))
    allocate(unitaryMat(nOrb,nOrb))

    ! Calculate CSC = C^T_old * S_new * C_old
    tmpMat(:,:) = 0.0_dp
    call gemm(tmpMat, eigvecsReal(:,:,1), reks%overSqr, transA='T')
    tmpC(:,:) = 0.0_dp
    call gemm(tmpC, tmpMat, eigvecsReal(:,:,1))

    ! Diagonalize CSC to obtain a unitary matrix, U = CSC^(-1/2)

    tmpEigen(:) = 0.0_dp
    call env%globalTimer%startTimer(globalTimers%diagonalization)
    ! tmpC becomes eigenvectors (X) of CSC
    call heev(tmpC, tmpEigen, 'U', 'V')
    call env%globalTimer%stopTimer(globalTimers%diagonalization)

    ! Make inverse square root matrix consisting of eigenvalues (s) of CSC
    tmpS(:,:) = 0.0_dp
    do ii = 1, nOrb
      tmpS(ii,ii) = 1.0_dp / max(sqrt(tmpEigen(ii)), epsilon(0.0_dp))
    end do

    ! Calculate a unitary matrix U = CSC^(-1/2) = X * s^(-1/2) * X^T
    tmpMat(:,:) = 0.0_dp
    call gemm(tmpMat, tmpS, tmpC, transB='T')
    unitaryMat(:,:) = 0.0_dp
    call gemm(unitaryMat, tmpC, tmpMat)

    ! C_new = C_old * U
    tmpC(:,:) = 0.0_dp
    call gemm(tmpC, eigvecsReal(:,:,1), unitaryMat)

    eigvecsReal(:,:,1) = tmpC

  end subroutine renormalizeEigenvecs


  !> Creates (delta) density matrix for each microstate from real eigenvectors.
  subroutine getDensityMatrixL(env, denseDesc, neighbourList, nNeighbourSK, iSparseStart,&
      & img2CentCell, orb, species, coord, tHelical, eigvecs, parallelKS, rhoPrim, work,&
      & rhoSqrReal, q0, deltaRhoOutSqr, reks)

    !> Environment settings
    type(TEnvironment), intent(inout) :: env

    !> Dense matrix descriptor
    type(TDenseDescr), intent(in) :: denseDesc

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

    !> species of all atoms
    integer, target, intent(in) :: species(:)

    !> Coordinates of all atoms including images
    real(dp), allocatable, intent(inout) :: coord(:,:)

    !> Is the geometry helical
    logical, intent(in) :: tHelical

    !> eigenvectors
    real(dp), intent(inout) :: eigvecs(:,:,:)

    !> K-points and spins to process
    type(TParallelKS), intent(in) :: parallelKS

    !> sparse density matrix
    real(dp), intent(inout) :: rhoPrim(:,:)

    !> work space array
    real(dp), intent(inout) :: work(:,:)

    !> Dense density matrix if needed
    real(dp), intent(inout), allocatable :: rhoSqrReal(:,:,:)

    !> reference atomic occupations
    real(dp), intent(in) :: q0(:,:,:)

    !> Change in density matrix during this SCC step for rangesep
    real(dp), pointer, intent(inout) :: deltaRhoOutSqr(:,:,:)

    !> data type for REKS
    type(TReksCalc), intent(inout) :: reks

    integer :: iL

    call env%globalTimer%startTimer(globalTimers%densityMatrix)

    if (reks%tForces) then
      reks%rhoSqrL(:,:,:,:) = 0.0_dp
    else
      reks%rhoSpL(:,:,:) = 0.0_dp
    end if

    do iL = 1, reks%Lmax

      call getDensityFromRealEigvecs(env, denseDesc, reks%fillingL(:,:,iL), neighbourList,&
          & nNeighbourSK, iSparseStart, img2CentCell, orb, species, denseDesc%iAtomStart, coord,&
          & tHelical, eigvecs, parallelKS, rhoPrim, work, rhoSqrReal, deltaRhoOutSqr)

      if (reks%tForces) then
        ! reks%rhoSqrL has (my_ud) component
        if (reks%isRangeSep) then
          reks%rhoSqrL(:,:,1,iL) = deltaRhoOutSqr(:,:,1)
        else
          reks%rhoSqrL(:,:,1,iL) = work
        end if
      else
        ! reks%rhoSpL has (my_ud) component
        reks%rhoSpL(:,1,iL) = rhoPrim(:,1)
      end if

      if (reks%isRangeSep) then
        ! reks%deltaRhoSqrL has (my_ud) component
        reks%deltaRhoSqrL(:,:,1,iL) = deltaRhoOutSqr(:,:,1)
      end if

      if (reks%tForces) then
        call symmetrizeHS(reks%rhoSqrL(:,:,1,iL))
      end if
      if (reks%isRangeSep) then
        call symmetrizeHS(reks%deltaRhoSqrL(:,:,1,iL))
        call denseSubtractDensityOfAtoms(q0, denseDesc%iAtomStart, reks%deltaRhoSqrL(:,:,:,iL), 1)
      end if

    end do

    if (reks%tForces) then
      ! reks%rhoSqrL has (my_qm) component
      call ud2qmL(reks%rhoSqrL, reks%Lpaired)
    else
      ! reks%rhoSpL has (my_qm) component
      call ud2qmL(reks%rhoSpL, reks%Lpaired)
    end if

    call env%globalTimer%stopTimer(globalTimers%densityMatrix)

  end subroutine getDensityMatrixL


  !> Calculate Mulliken population for each microstate from sparse density matrix.
  subroutine getMullikenPopulationL(env, denseDesc, neighbourList, nNeighbourSK, &
      & img2CentCell, iSparseStart, orb, rhoPrim, ints, iRhoPrim, qBlock, &
      & qiBlock, qNetAtom, reks)

    !> Environment settings
    type(TEnvironment), intent(inout) :: env

    !> Dense matrix descriptor
    type(TDenseDescr), intent(in) :: denseDesc

    !> Atomic neighbours
    type(TNeighbourList), intent(in) :: neighbourList

    !> Number of neighbours for each atom within overlap distance
    integer, intent(in) :: nNeighbourSK(:)

    !> image to actual atom indexing
    integer, intent(in) :: img2CentCell(:)

    !> sparse matrix indexing array
    integer, intent(in) :: iSparseStart(:,:)

    !> Atomic orbital information
    type(TOrbitals), intent(in) :: orb

    !> sparse density matrix
    real(dp), intent(inout) :: rhoPrim(:,:)

    !> Integral container
    type(TIntegral), intent(in) :: ints

    !> imaginary part of density matrix
    real(dp), intent(in), allocatable :: iRhoPrim(:,:)

    !> Dual atomic charges
    real(dp), intent(inout), allocatable :: qBlock(:,:,:,:)

    !> Imaginary part of dual atomic charges
    real(dp), intent(inout), allocatable :: qiBlock(:,:,:,:)

    !> Onsite Mulliken charges per atom
    real(dp), intent(inout), allocatable :: qNetAtom(:)

    !> data type for REKS
    type(TReksCalc), intent(inout) :: reks

    integer :: iL

    do iL = 1, reks%Lmax

      if (reks%tForces) then
        rhoPrim(:,1) = 0.0_dp
        call env%globalTimer%startTimer(globalTimers%denseToSparse)
        call packHS(rhoPrim(:,1), reks%rhoSqrL(:,:,1,iL), neighbourlist%iNeighbour, &
            & nNeighbourSK, orb%mOrb, denseDesc%iAtomStart, iSparseStart, img2CentCell)
        call env%globalTimer%stopTimer(globalTimers%denseToSparse)
      else
        rhoPrim(:,1) = reks%rhoSpL(:,1,iL)
      end if

      ! reks%qOutputL & reks%qNetAtomL has (my_qm) component
      reks%qOutputL(:,:,:,iL) = 0.0_dp
      call getMullikenPopulation(rhoPrim, ints, orb, neighbourList, nNeighbourSK, &
          & img2CentCell, iSparseStart, reks%qOutputL(:,:,:,iL), iRhoPrim=iRhoPrim, &
          & qBlock=qBlock, qiBlock=qiBlock, qNetAtom=qNetAtom)

      ! Get correct net charge per atom
      ! Note that qNetAtomL does not have spin dependency so it does not
      ! correspond to (my_qm) or (my_ud) representation
      if (reks%isQNetAllocated) then
        if (iL > reks%Lpaired) then
          if (mod(iL,2) == 0) then
            reks%qNetAtomL(:,iL) = reks%qNetAtomL(:,iL-1)
          else
            reks%qNetAtomL(:,iL) = qNetAtom
          end if
        else
          reks%qNetAtomL(:,iL) = qNetAtom
        end if
      end if

    end do

    ! reks%qOutputL has (qm) component
    call qmExpandL(reks%qOutputL, reks%Lpaired)

  end subroutine getMullikenPopulationL


  !> Build L, spin dependent Hamiltonian with various contributions
  !> and compute the energy of microstates
  subroutine getHamiltonianLandEnergyL(env, denseDesc, sccCalc, tblite, orb, species,&
      & neighbourList, nNeighbourSK, iSparseStart, img2CentCell, H0, ints, spinW, cellVol,&
      & extPressure, energy, q0, iAtInCentralRegion, solvation, thirdOrd, potential, rangeSep,&
      & nNeighbourLC, tDualSpinOrbit, xi, isExtField, isXlbomd, dftbU, TS, qDepExtPot, qBlock,&
      & qiBlock, tFixEf, Ef, rhoPrim, onSiteElements, dispersion, tConverged, species0,&
      & referenceN0, qNetAtom, multipole, reks)

    !> Environment settings
    type(TEnvironment), intent(inout) :: env

    !> Dense matrix descriptor
    type(TDenseDescr), intent(in) :: denseDesc

    !> SCC module internal variables
    type(TScc), allocatable, intent(inout) :: sccCalc

    !> Library interface handler
    type(TTBLite), allocatable, intent(inout) :: tblite

    !> atomic orbital information
    type(TOrbitals), intent(in) :: orb

    !> species of all atoms
    integer, target, intent(in) :: species(:)

    !> neighbours to atoms
    type(TNeighbourList), intent(in) :: neighbourList

    !> Number of atomic neighbours
    integer, intent(in) :: nNeighbourSK(:)

    !> Index for atomic blocks in sparse data
    integer, intent(in) :: iSparseStart(:,:)

    !> map from image atom to real atoms
    integer, intent(in) :: img2CentCell(:)

    !> non-SCC hamiltonian (sparse)
    real(dp), intent(in) :: H0(:)

    !> Integral container
    type(TIntegral), intent(inout) :: ints

    !> spin constants
    real(dp), allocatable, intent(in) :: spinW(:,:,:)

    !> unit cell volume
    real(dp), intent(in) :: cellVol

    !> external pressure
    real(dp), intent(in) :: extPressure

    !> energy contributions
    type(TEnergies), intent(inout) :: energy

    !> reference atomic occupations
    real(dp), intent(in) :: q0(:,:,:)

    !> Atoms over which to sum the total energies
    integer, intent(in) :: iAtInCentralRegion(:)

    !> Solvation mode
    class(TSolvation), allocatable, intent(inout) :: solvation

    !> third order SCC interactions
    type(TThirdOrder), allocatable, intent(inout) :: thirdOrd

    !> potentials acting
    type(TPotentials), intent(inout) :: potential

    !> Data for rangeseparated calculation
    type(TRangeSepFunc), allocatable, intent(inout) :: rangeSep

    !> Nr. of neighbours for each atom in the long-range functional.
    integer, allocatable, intent(in) :: nNeighbourLC(:)

    !> Is dual spin orbit being used (block potentials)
    logical, intent(in) :: tDualSpinOrbit

    !> Spin orbit constants if required
    real(dp), allocatable, intent(in) :: xi(:,:)

    !> is an external electric field present
    logical, intent(in) :: isExtField

    !> Is the extended Lagrangian being used for MD
    logical, intent(in) :: isXlbomd

    !> Are there orbital potentials present
    type(TDftbU), intent(in), allocatable :: dftbU

    !> electron entropy contribution
    real(dp), intent(in) :: TS(:)

    !> Proxy for querying Q-dependant external potentials
    type(TQDepExtPotProxy), intent(inout), allocatable :: qDepExtPot

    !> block (dual) atomic populations
    real(dp), intent(in), allocatable :: qBlock(:,:,:,:)

    !> Imaginary part of block atomic populations
    real(dp), intent(in), allocatable :: qiBlock(:,:,:,:)

    !> Whether fixed Fermi level(s) should be used. (No charge conservation!)
    logical, intent(in) :: tFixEf

    !> If tFixEf is .true. contains reservoir chemical potential, otherwise the Fermi levels found
    !> from the given number of electrons
    real(dp), intent(inout) :: Ef(:)

    !> sparse density matrix
    real(dp), intent(inout) :: rhoPrim(:,:)

    !> Corrections terms for on-site elements
    real(dp), intent(in), allocatable :: onSiteElements(:,:,:,:)

    !> dispersion interactions
    class(TDispersionIface), allocatable, intent(inout) :: dispersion

    !> Has the calculation converged>
    logical, intent(in) :: tConverged

    !> species of atoms in the central cell
    integer, intent(in) :: species0(:)

    !> reference n_0 charges for each atom
    real(dp), intent(in) :: referenceN0(:,:)

    !> Onsite Mulliken charges per atom
    real(dp), intent(inout), allocatable :: qNetAtom(:)

    !> Multipole container
    type(TMultipole), intent(inout) :: multipole

    !> data type for REKS
    type(TReksCalc), allocatable, intent(inout) :: reks

    real(dp), allocatable :: tmpHamSp(:,:)
    real(dp), allocatable :: tmpEn(:)

    integer, pointer :: pSpecies0(:)
    integer :: sparseSize, nAtom, nSpin, iL, tmpL, rsL

    sparseSize = size(ints%overlap,dim=1)
    nAtom = size(reks%qOutputL,dim=2)
    nSpin = size(reks%qOutputL,dim=3)
    pSpecies0 => species(1:nAtom)

    allocate(tmpHamSp(sparseSize,1))
    if (reks%isRangeSep) then
      allocate(tmpEn(reks%Lmax))
    end if

    ! Calculate contribution to Hamiltonian except rangeseparated part
    reks%intShellL(:,:,:,:) = 0.0_dp
    reks%intBlockL(:,:,:,:,:) = 0.0_dp
    do iL = 1, reks%Lmax

      ! reks%chargePerShellL has (qm) component
      call getChargePerShell(reks%qOutputL(:,:,:,iL), orb, species,&
          & reks%chargePerShellL(:,:,:,iL))
      call resetInternalPotentials(tDualSpinOrbit, xi, orb, species, potential)
      call addChargePotentials(env, sccCalc, tblite, .true., reks%qOutputL(:,:,:,iL), q0,&
          & reks%chargePerShellL(:,:,:,iL), orb, multipole, species, neighbourList,&
          & img2CentCell, spinW, solvation, thirdOrd, dispersion, potential)

      ! reks%intShellL, reks%intBlockL has (qm) component
      reks%intShellL(:,:,:,iL) = potential%intShell
      reks%intBlockL(:,:,:,:,iL) = potential%intBlock

      ! qm representation is converted to my_qm representation
      if (iL <= reks%Lpaired) then
        ! If iL = 1, then qm = 1u + 1d, 1u - 1d and my_qm = 1u + 1d
        ! calculate charge part
        potential%intBlock(:,:,:,1) = reks%intBlockL(:,:,:,1,iL)
        tmpHamSp(:,1) = h0
      else
        if (mod(iL,2) == 1) then
          ! If iL = 3, then qm = 3u + 3d, 3u - 3d and my_qm = 3u + 3d
          ! calculate charge part
          potential%intBlock(:,:,:,1) = reks%intBlockL(:,:,:,1,iL)
          tmpHamSp(:,1) = h0
        else
          ! If iL = 4, then qm = 4u + 4d, 4u - 4d and my_qm = 3u - 3d (= -(4u - 4d))
          ! calculate magnetization part
          potential%intBlock(:,:,:,1) = -reks%intBlockL(:,:,:,2,iL)
          tmpHamSp(:,1) = 0.0_dp
        end if
      end if

      ! tmpHamSp has (my_qm) component
      call getSccHamiltonian(H0, ints, nNeighbourSK, neighbourList, species, orb,&
          & iSparseStart, img2CentCell, potential, allocated(reks), tmpHamSp, ints%iHamiltonian)
      tmpHamSp(:,1) = 2.0_dp * tmpHamSp(:,1)

      if (reks%isRangeSep) then
        ! reks%hamSqrL has (my_qm) component
        reks%hamSqrL(:,:,1,iL) = 0.0_dp
        call env%globalTimer%startTimer(globalTimers%sparseToDense)
        call unpackHS(reks%hamSqrL(:,:,1,iL), tmpHamSp(:,1), neighbourList%iNeighbour, &
            & nNeighbourSK, denseDesc%iAtomStart, iSparseStart, img2CentCell)
        call env%globalTimer%stopTimer(globalTimers%sparseToDense)
        call blockSymmetrizeHS(reks%hamSqrL(:,:,1,iL), denseDesc%iAtomStart)
      else
        ! reks%hamSpL has (my_qm) component
        reks%hamSpL(:,1,iL) = tmpHamSp(:,1)
      end if

    end do

    ! Calculate contribution to Hamiltonian including rangeseparated part
    if (.not. reks%isRangeSep) then
      ! reks%hamSpL has (my_ud) component
      call qm2udL(reks%hamSpL, reks%Lpaired)
    else
      ! reks%hamSqrL has (my_ud) component
      call qm2udL(reks%hamSqrL, reks%Lpaired)
      tmpEn(:) = 0.0_dp
      do iL = 1, reks%Lmax
        ! Add rangeseparated contribution
        call rangeSep%addLRHamiltonian(env, reks%deltaRhoSqrL(:,:,1,iL), ints%overlap, &
            & neighbourList%iNeighbour, nNeighbourLC, denseDesc%iAtomStart, &
            & iSparseStart, orb, reks%hamSqrL(:,:,1,iL), reks%overSqr)
        ! Calculate the long-range exchange energy for up spin
        call rangeSep%addLREnergy(tmpEn(iL))
      end do
    end if


    ! Calculate energy contribution corresponding to upper Hamiltonian
    do iL = 1, reks%Lmax

      ! Get microstate index for non-SCC and rangeseparation energy contribution
      if (iL <= reks%Lpaired) then
        tmpL = iL
        rsL = iL
      else
        if (mod(iL,2) == 1) then
          tmpL = iL
          rsL = iL + 1
        else
          tmpL = iL - 1
          rsL = iL - 1
        end if
      end if
      ! Set the long-range corrected energy contribution
      if (reks%isRangeSep) then
        energy%Efock = tmpEn(iL) + tmpEn(rsL)
      end if

      if (reks%tForces) then
        rhoPrim(:,1) = 0.0_dp
        call env%globalTimer%startTimer(globalTimers%denseToSparse)
        ! reks%rhoSqrL has (my_qm) component
        call packHS(rhoPrim(:,1), reks%rhoSqrL(:,:,1,tmpL), &
            & neighbourList%iNeighbour, nNeighbourSK, orb%mOrb, &
            & denseDesc%iAtomStart, iSparseStart, img2CentCell)
        call env%globalTimer%stopTimer(globalTimers%denseToSparse)
      else
        ! reks%rhoSpL has (my_qm) component
        rhoPrim(:,1) = reks%rhoSpL(:,1,tmpL)
      end if

      ! Calculate correct charge contribution for each microstate
      call sccCalc%updateCharges(env, reks%qOutputL(:,:,:,iL), orb, species, q0)
      call sccCalc%updateShifts(env, orb, species, neighbourList%iNeighbour, img2CentCell)
      potential%intShell(:,:,:) = reks%intShellL(:,:,:,iL)
      if (allocated(thirdOrd)) then
        call thirdOrd%updateCharges(pSpecies0, neighbourList, reks%qOutputL(:,:,:,iL), q0,&
            & img2CentCell, orb)
      end if

      call calcEnergies(sccCalc, tblite, reks%qOutputL(:,:,:,iL), q0,&
          & reks%chargePerShellL(:,:,:,iL), multipole, species, isExtField, isXlbomd, dftbU,&
          & tDualSpinOrbit, rhoPrim, H0, orb, neighbourList, nNeighbourSk, img2CentCell,&
          & iSparseStart, cellVol, extPressure, TS, potential, energy, thirdOrd, solvation,&
          & rangeSep, reks, qDepExtPot, qBlock, qiBlock, xi, iAtInCentralRegion, tFixEf, Ef,&
          & onSiteElements)

      if (allocated(dispersion)) then
        ! For dftd4 dispersion, update charges
        call dispersion%updateCharges(env, pSpecies0, neighbourList, reks%qOutputL(:,:,:,iL),&
            & q0, img2CentCell, orb)
        ! For MBD/TS dispersion, update onsite charges
        ! TODO : Currently, reks%qNetAtomL does not affect Hamiltonian
        if (reks%isQNetAllocated) then
          qNetAtom(:) = reks%qNetAtomL(:,iL)
        end if
        call dispersion%updateOnsiteCharges(qNetAtom, orb, referenceN0,&
            & species0, tConverged)
        call calcDispersionEnergy(dispersion, energy%atomDisp, energy%Edisp,&
            & iAtInCentralRegion)
      end if
      call sumEnergies(energy)

      ! Assign energy contribution of each microstate
      reks%enLnonSCC(iL) = energy%EnonSCC
      reks%enLSCC(iL) = energy%Escc
      reks%enLspin(iL) = energy%Espin
      if (allocated(thirdOrd)) then
        reks%enL3rd(iL) = energy%e3rd
      end if
      if (reks%isRangeSep) then
        reks%enLfock(iL) = energy%Efock
      end if
      if (reks%isDispersion) then
        reks%enLdisp(iL) = energy%Edisp
      end if
      reks%enLtot(iL) = energy%Etotal

      ! REKS is not affected by filling, so TS becomes 0
      energy%EMermin = energy%Etotal
      ! extrapolated to 0 K
      energy%Ezero = energy%Etotal
      energy%EGibbs = energy%EMermin + cellVol * extPressure
      energy%EForceRelated = energy%EGibbs

    end do

    if (reks%Plevel >= 2) then
      call printReksMicrostates(reks, energy%Erep)
    end if

  end subroutine getHamiltonianLandEnergyL


  !> Optimize the fractional occupation numbers (FONs) and weights
  !> Swap the active orbitals when fa < fb
  !> Compute the several energy contributions
  subroutine optimizeFONsAndWeights(eigvecs, filling, energy, reks)

    !> eigenvectors
    real(dp), intent(inout) :: eigvecs(:,:,:)

    !> occupations (level, kpoint, spin)
    real(dp), intent(out) :: filling(:,:,:)

    !> energy contributions
    type(TEnergies), intent(inout) :: energy

    !> data type for REKS
    type(TReksCalc), intent(inout) :: reks

    call optimizeFons(reks)
    call calcWeights(reks)

    call activeOrbSwap(reks, eigvecs(:,:,1))
    call getFilling(reks, filling(:,1,1))

    call calcSaReksEnergy(reks, energy)

    if (reks%Plevel >= 2) then
      call printSaReksEnergy(reks)
    end if

  end subroutine optimizeFONsAndWeights


  !> Returns input charges for next SCC iteration.
  subroutine getReksNextInputCharges(qInput, qOutput, qDiff, sccErrorQ, sccTol, tConverged,&
      & iSccIter, minSccIter, maxSccIter, iGeoStep, tStopScc, eigvecs, reks)

    !> input charges (for potentials)
    real(dp), intent(inout) :: qInput(:, :, :)

    !> Output electrons
    real(dp), intent(inout) :: qOutput(:,:,:)

    !> charge differences between input and output charges
    real(dp), intent(inout) :: qDiff(:,:,:)

    !> SCC error
    real(dp), intent(out) :: sccErrorQ

    !> Tolerance on SCC charges between input and output
    real(dp), intent(in) :: sccTol

    !> Has the calculation converged>
    logical, intent(out) :: tConverged

    !> Number of current SCC step
    integer, intent(in) :: iSccIter

    !> minimum number of SCC iterations to perform
    integer, intent(in) :: minSccIter

    !> maximum number of SCC iterations before terminating loop
    integer, intent(in) :: maxSccIter

    !> Number of current geometry step
    integer, intent(in) :: iGeoStep

    !> Should the SCC loop stop
    logical, intent(in) :: tStopScc

    !> Eigenvectors on eixt
    real(dp), intent(inout) :: eigvecs(:,:,:)

    !> data type for REKS
    type(TReksCalc), intent(inout) :: reks

    qDiff(:,:,:) = qOutput - qInput
    sccErrorQ = maxval(abs(qDiff))

    tConverged = (sccErrorQ < sccTol) &
        & .and. (iSccIter >= minSccIter .or. reks%tReadMO .or. iGeoStep > 0)
    if ((.not. tConverged) .and. (iSccIter /= maxSccIter .and. .not. tStopScc)) then
      qInput(:,:,:) = qOutput
      call guessNewEigvecs(eigvecs(:,:,1), reks%eigvecsFock)
    end if

  end subroutine getReksNextInputCharges


  !> Update delta density matrix rather than merely q for rangeseparation
  subroutine getReksNextInputDensity(sccErrorQ, sccTol, tConverged, &
      & iSccIter, minSccIter, maxSccIter, iGeoStep, tStopScc, &
      & eigvecs, deltaRhoOut, deltaRhoIn, deltaRhoDiff, reks)

    !> SCC error
    real(dp), intent(out) :: sccErrorQ

    !> Tolerance on SCC charges between input and output
    real(dp), intent(in) :: sccTol

    !> Has the calculation converged>
    logical, intent(out) :: tConverged

    !> Number of current SCC step
    integer, intent(in) :: iSccIter

    !> minimum number of SCC iterations to perform
    integer, intent(in) :: minSccIter

    !> maximum number of SCC iterations before terminating loop
    integer, intent(in) :: maxSccIter

    !> Number of current geometry step
    integer, intent(in) :: iGeoStep

    !> Should the SCC loop stop
    logical, intent(in) :: tStopScc

    !> Eigenvectors on eixt
    real(dp), intent(inout) :: eigvecs(:,:,:)

    !> delta density matrix for rangeseparated calculations
    real(dp), intent(inout) :: deltaRhoOut(:)

    !> delta density matrix as inpurt for next SCC cycle
    real(dp), target, intent(inout) :: deltaRhoIn(:)

    !> difference of delta density matrix in and out
    real(dp), intent(inout) :: deltaRhoDiff(:)

    !> data type for REKS
    type(TReksCalc), intent(inout) :: reks

    deltaRhoDiff(:) = deltaRhoOut - deltaRhoIn
    sccErrorQ = maxval(abs(deltaRhoDiff))

    tConverged = (sccErrorQ < sccTol) &
        & .and. (iSccIter >= minSccIter .or. reks%tReadMO .or. iGeoStep > 0)
    if ((.not. tConverged) .and. (iSccIter /= maxSccIter .and. .not. tStopScc)) then
      deltaRhoIn(:) = deltaRhoOut
      call guessNewEigvecs(eigvecs(:,:,1), reks%eigvecsFock)
    end if

  end subroutine getReksNextInputDensity


  !> Set correct dipole moment according to type of REKS calculation
  subroutine assignDipoleMoment(dipoleTmp, dipoleMoment, iDet, tDipole, reks, isSingleState)

    !> resulting temporary dipole moment
    real(dp), allocatable, intent(in) :: dipoleTmp(:)

    !> resulting dipole moment
    real(dp), allocatable, intent(inout) :: dipoleMoment(:,:)

    !> Which state is being calculated in the determinant loop?
    integer, intent(in) :: iDet

    !> calculate an electric dipole?
    logical, intent(in) :: tDipole

    !> data type for REKS
    type(TReksCalc), intent(inout) :: reks

    !> calculate a single-state REKS?
    logical, intent(in) :: isSingleState

    ! Set correct dipole moment to this%dipoleMoment
    if (tDipole) then
      if (isSingleState) then
        ! For single-state REKS case, see getStateInteraction routine.
        if (reks%Efunction == 1) then
          dipoleMoment(:,iDet) = dipoleTmp
        end if
      else
        ! (SI)-SA-REKS case, see getReksGradProperties routine.
        if (reks%Efunction > 1) then
          dipoleMoment(:,iDet) = dipoleTmp
        end if
      end if
    end if

  end subroutine assignDipoleMoment

end module dftbp_dftbplus_main
