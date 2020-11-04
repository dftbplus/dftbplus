!--------------------------------------------------------------------------------------------------!
!  DFTB+: general package for performing fast atomistic simulations                                !
!  Copyright (C) 2006 - 2020  DFTB+ developers group                                               !
!                                                                                                  !
!  See the LICENSE file for terms of usage and distribution.                                       !
!--------------------------------------------------------------------------------------------------!

#:include 'common.fypp'
#:include 'error.fypp'

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
  use dftbp_hamiltoniantypes
  use dftbp_nonscc
  use dftbp_eigenvects
  use dftbp_repulsive
  use dftbp_etemp
  use dftbp_populations
  use dftbp_densitymatrix
  use dftbp_forces
  use dftbp_stress
  use dftbp_scc
  use dftbp_hamiltonian
  use dftbp_getenergies, only : calcEnergies, calcRepulsiveEnergy, calcDispersionEnergy, sumEnergies
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
  use dftbp_energytypes, only : TEnergies
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
  use dftbp_linresptypes
  use dftbp_pprpa, only : ppRPAenergies
#:if WITH_ARPACK
  use dftbp_RS_LinearResponse
#:endif
  use dftbp_mainio
  use dftbp_commontypes
  use dftbp_dispersions, only : TDispersionIface
  use dftbp_solvation, only : TSolvation
  use dftbp_cm5, only : TChargeModel5
  use dftbp_xmlf90
  use dftbp_thirdorder, only : TThirdOrder
  use dftbp_rangeseparated, only : TRangeSepFunc
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
  use dftbp_timeprop
  use dftbp_qdepextpotproxy, only : TQDepExtPotProxy
  use dftbp_taggedoutput, only : TTaggedWriter
  use dftbp_reks
  use dftbp_plumed, only : TPlumedCalc, TPlumedCalc_final
  use dftbp_determinants
#:if WITH_TRANSPORT
  use libnegf_vars, only : TTransPar
  use negf_int
#:endif
  use poisson_init
  use dftbp_transportio
  use dftbp_initprogram

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
  subroutine runDftbPlus(env, globalData)

    !> Environment settings
    type(TEnvironment), intent(inout) :: env

    !> Global variables
    type(TGlobalData), intent(inout) :: globalData

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

    call initGeoOptParameters(globalData%tCoordOpt, globalData%nGeoSteps, tGeomEnd, tCoordStep,&
        & tStopDriver, iGeoStep, iLatGeoStep)

    ! If the geometry is periodic, need to update lattice information in geometry loop
    globalData%tLatticeChanged = globalData%tPeriodic

    ! As first geometry iteration, require updates for coordinates in dependent routines
    globalData%tCoordsChanged = .true.

    ! Main geometry loop
    geoOpt: do iGeoStep = 0, globalData%nGeoSteps
      tWriteRestart = env%tGlobalLead .and. needsRestartWriting(globalData%isGeoOpt,&
          & globalData%tMd, iGeoStep, globalData%nGeoSteps, globalData%restartFreq)

      if (.not. globalData%tRestartNoSC) then
        call printGeoStepInfo(globalData%tCoordOpt, globalData%tLatOpt, iLatGeoStep, iGeoStep)
      end if

      ! DFTB Determinant Loop
      ! Will pass though loop once, unless specified in input to perform multiple determiants
      lpDets : do iDet = 1, globalData%nDets

        globalData%deltaDftb%iDeterminant = iDet

        call preDetCharges(isUnReduced, iDet, globalData%nDets, iGeoStep, globalData%deltaDftb,&
            & globalData%qInput, globalData%qDets, globalData%qBlockIn, globalData%qBlockDets,&
            & globalData%deltaRhoIn, globalData%deltaRhoDets)
        if (isUnReduced) then
          call reduceCharges(globalData%orb, globalData%nIneqOrb, globalData%iEqOrbitals,&
              & globalData%qInput, globalData%qInpRed, globalData%qBlockIn,&
              & globalData%iEqBlockDftbu, globalData%qiBlockIn, globalData%iEqBlockDftbuLS,&
              & globalData%iEqBlockOnSite, globalData%iEqBlockOnSiteLS)
        end if

        call processGeometry(env, globalData, iGeoStep, iLatGeoStep, tWriteRestart, tStopScc,&
            & tExitGeoOpt)

        call postDetCharges(iDet, globalData%nDets, globalData%qOutput, globalData%qDets,&
            & globalData%qBlockDets, globalData%qBlockOut, globalData%deltaRhoDets,&
            & globalData%deltaRhoOut)

      end do lpDets

      call globalData%deltaDftb%postProcessDets(globalData%dftbEnergy, globalData%qOutput,&
          & globalData%qDets, globalData%qBlockOut, globalData%qBlockDets,&
          & globalData%dipoleMoment, globalData%totalStress, globalData%tripletStress,&
          & globalData%mixedStress, globalData%derivs, globalData%tripletderivs,&
          & globalData%mixedderivs)

      if (globalData%tWriteDetailedOut .and. globalData%deltaDftb%nDeterminant() > 1) then
        call writeDetailedOut2Dets(globalData%fdDetailedOut, userOut,&
            & tAppendDetailedOut, globalData%dftbEnergy, globalData%electronicSolver,&
            & globalData%deltaDftb, globalData%q0, globalData%orb, globalData%qOutput,&
            & globalData%qDets, globalData%qBlockDets, globalData%species,&
            & globalData%iAtInCentralRegion, globalData%tPrintMulliken, globalData%cm5Cont)
      end if

      if (.not.globalData%tRestartNoSC) then
        call printEnergies(globalData%dftbEnergy, globalData%electronicSolver,&
            & globalData%deltaDftb)
      end if

      if (globalData%tStress) then
        globalData%intPressure = (globalData%totalStress(1,1) + globalData%totalStress(2,2) +&
            & globalData%totalStress(3,3)) / 3.0_dp
        globalData%totalLatDeriv(:,:) = -globalData%cellVol * matmul(globalData%totalStress,&
            & globalData%invLatVec)

        call printVolume(globalData%cellVol)

        ! MD case includes the atomic kinetic energy contribution, so print that later
        if (.not. (globalData%tMD .or. globalData%tHelical)) then
          call printPressureAndFreeEnergy(globalData%extPressure, globalData%intPressure,&
              & globalData%dftbEnergy(globalData%deltaDftb%iDeterminant)%EGibbs)
        end if
      end if

      call postprocessDerivs(globalData%derivs, globalData%conAtom, globalData%conVec,&
          & globalData%tLatOpt, globalData%totalLatDeriv, globalData%extLatDerivs,&
          & globalData%normOrigLatVec, globalData%tLatOptFixAng, globalData%tLatOptFixLen,&
          & globalData%tLatOptIsotropic, constrLatDerivs)

      if (tExitGeoOpt) then
        exit geoOpt
      end if

      call printMaxForces(globalData%derivs, constrLatDerivs, globalData%tCoordOpt,&
          & globalData%tLatOpt, globalData%indMovedAtom)
    #:if WITH_SOCKETS
      if (globalData%tSocket) then
        call sendEnergyAndForces(env, globalData%socket,&
            & globalData%dftbEnergy(globalData%deltaDftb%iFinal), globalData%derivs,&
            & globalData%totalStress, globalData%cellVol)
      end if
    #:endif
      tWriteCharges =  allocated(globalData%qInput) .and. tWriteRestart .and. globalData%tMulliken&
          & .and. globalData%tSccCalc .and. .not. globalData%tDerivs .and.&
          & globalData%maxSccIter > 1 .and. globalData%deltaDftb%nDeterminant() == 1
      if (tWriteCharges) then
        call writeCharges(fCharges, globalData%tWriteChrgAscii, globalData%orb,&
            & globalData%qInput, globalData%qBlockIn, globalData%qiBlockIn, globalData%deltaRhoIn)
      end if

      if (globalData%tForces) then
        call getNextGeometry(env, globalData, iGeoStep, tWriteRestart, constrLatDerivs,&
            & tCoordStep, tGeomEnd, tStopDriver, iLatGeoStep, tempIon, tExitGeoOpt)
        if (tExitGeoOpt) then
          exit geoOpt
        end if
      end if

      if (globalData%tWriteDetailedOut .and. globalData%tMd) then
        call writeDetailedOut6(globalData%fdDetailedOut,&
            & globalData%dftbEnergy(globalData%deltaDftb%iFinal), tempIon)
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
    if (globalData%tSocket .and. env%tGlobalLead) then
      call globalData%socket%shutdown()
    end if
  #:endif

    if (allocated(globalData%plumedCalc)) then
      call TPlumedCalc_final(globalData%plumedCalc)
    end if

    tGeomEnd = globalData%tMD .or. tGeomEnd .or. globalData%tDerivs

    if (env%tGlobalLead) then
      if (globalData%tWriteDetailedOut) then
        call writeDetailedOut7(globalData%fdDetailedOut, globalData%isGeoOpt, tGeomEnd,&
            & globalData%tMd, globalData%tDerivs, globalData%tEField, globalData%absEField,&
            & globalData%dipoleMoment, globalData%deltaDftb)
      end if

      call writeFinalDriverStatus(globalData%isGeoOpt, tGeomEnd, globalData%tMd,&
          & globalData%tDerivs)

      if (globalData%tMD) then
        call writeMdOut3(globalData%fdMd, mdOut)
      end if
    end if

    if (env%tGlobalLead .and. globalData%tDerivs) then
      call getHessianMatrix(globalData%derivDriver, globalData%pDynMatrix)
      call writeHessianOut(hessianOut, globalData%pDynMatrix)
    else
      nullify(globalData%pDynMatrix)
    end if

    if (globalData%tWriteShifts) then
      call writeShifts(fShifts, globalData%orb, globalData%potential%intShell)
    endif

    ! Here time propagation is called
    if (allocated(globalData%electronDynamics)) then
      call runDynamics(globalData%electronDynamics, globalData%eigvecsReal, globalData%ham,&
          & globalData%H0, globalData%species, globalData%q0, globalData%referenceN0,&
          & globalData%over, globalData%filling, globalData%neighbourList,&
          & globalData%nNeighbourSK, globalData%nNeighbourLC, globalData%denseDesc%iAtomStart,&
          & globalData%iSparseStart, globalData%img2CentCell, globalData%orb, globalData%coord0,&
          & globalData%spinW, globalData%pRepCont, globalData%sccCalc, env,&
          & globalData%tDualSpinOrbit, globalData%xi, globalData%thirdOrd, globalData%solvation,&
          & globalData%rangeSep, globalData%qDepExtPot, globalData%nDftbUFunc, globalData%UJ,&
          & globalData%nUJ, globalData%iUJ, globalData%niUJ, globalData%iAtInCentralRegion,&
          & globalData%tFixEf, globalData%Ef, globalData%coord, globalData%onsiteElements,&
          & globalData%skHamCont, globalData%skOverCont, globalData%latVec, globalData%invLatVec,&
          & globalData%iCellVec, globalData%rCellVec, globalData%cellVec,&
          & globalData%electronicSolver, globalData%eigvecsCplx, globalData%taggedWriter,&
          & globalData%refExtPot)
    end if

  #:if WITH_TRANSPORT
    if (globalData%tContCalc) then
      ! Note: shift and charges are saved in QM representation (not UD)
      call writeContactShifts(globalData%transpar%contacts(globalData%transpar%taskContInd)%name,&
          & globalData%orb, globalData%potential%intShell, globalData%qOutput, globalData%Ef,&
          & globalData%potential%intBlock, globalData%qBlockOut,&
          & .not.globalData%transpar%tWriteBinShift)
    end if

    if (globalData%tLocalCurrents) then
      call local_currents(env, globalData%parallelKS%localKS, globalData%ham, globalData%over,&
          & globalData%neighbourList, globalData%nNeighbourSK, globalData%cutOff%skCutoff,&
          & globalData%denseDesc%iAtomStart, globalData%iSparseStart, globalData%img2CentCell,&
          & globalData%iCellVec, globalData%cellVec, globalData%rCellVec, globalData%orb,&
          & globalData%kPoint, globalData%kWeight, globalData%coord0Fold, globalData%species0,&
          & globalData%speciesName, globalData%mu, globalData%lCurrArray)
    end if

    if (globalData%tTunn) then
      call calc_current(env, globalData%parallelKS%localKS, globalData%ham, globalData%over,&
          & globalData%neighbourList%iNeighbour, globalData%nNeighbourSK,&
          & globalData%densedesc%iAtomStart, globalData%iSparseStart, globalData%img2CentCell,&
          & globalData%iCellVec, globalData%cellVec, globalData%orb, globalData%kPoint,&
          & globalData%kWeight, globalData%tunneling, globalData%current, globalData%ldos,&
          & globalData%leadCurrents, globalData%writeTunn, globalData%tWriteLDOS,&
          & globalData%regionLabelLDOS, globalData%mu)
    end if

  #:endif

    if (allocated(globalData%pipekMezey)) then
      ! NOTE: the canonical DFTB ground state orbitals are over-written after this point
      if (withMpi) then
        call error("Pipek-Mezey localisation does not yet work with MPI")
      end if
      if (globalData%nSpin > 2) then
        call error("Pipek-Mezey localisation not implemented for non-colinear DFTB")
      end if
      call calcPipekMezeyLocalisation(env, globalData%pipekMezey, globalData%tPrintEigvecsTxt,&
          & globalData%nEl, globalData%filling, globalData%over, globalData%kPoint,&
          & globalData%neighbourList, globalData%nNeighbourSk, globalData%denseDesc,&
          & globalData%iSparseStart, globalData%img2CentCell, globalData%iCellVec,&
          & globalData%cellVec, globalData%runId, globalData%orb, globalData%species,&
          & globalData%speciesName, globalData%parallelKS, localisation, globalData%eigvecsReal,&
          & globalData%SSqrReal, globalData%eigvecsCplx, globalData%SSqrCplx, globalData%tHelical,&
          & globalData%coord)
    end if

    if (globalData%tWriteAutotest.and..not.globalData%tRestartNoSC) then
      if (globalData%tPeriodic) then
        globalData%cellVol = abs(determinant33(globalData%latVec))
        globalData%dftbEnergy(globalData%deltaDftb%iFinal)%EGibbs =&
            & globalData%dftbEnergy(globalData%deltaDftb%iFinal)%EMermin +&
            & globalData%extPressure * globalData%cellVol
      end if
      call writeAutotestTag(autotestTag, globalData%electronicSolver,&
          & globalData%tPeriodic, globalData%cellVol, globalData%tMulliken, globalData%qOutput,&
          & globalData%derivs, globalData%chrgForces, globalData%excitedDerivs,&
          & globalData%tStress, globalData%totalStress, globalData%pDynMatrix,&
          & globalData%dftbEnergy(globalData%deltaDftb%iFinal), globalData%extPressure,&
          & globalData%coord0, globalData%tLocalise, localisation, globalData%esp,&
          & globalData%taggedWriter, globalData%tunneling, globalData%ldos, globalData%lCurrArray)
    end if
    if (globalData%tWriteResultsTag) then
      call writeResultsTag(resultsTag,&
          & globalData%dftbEnergy(globalData%deltaDftb%iFinal), globalData%derivs,&
          & globalData%chrgForces, globalData%nEl, globalData%Ef, globalData%eigen,&
          & globalData%filling, globalData%electronicSolver, globalData%tStress,&
          & globalData%totalStress, globalData%pDynMatrix, globalData%tPeriodic,&
          & globalData%cellVol, globalData%tMulliken, globalData%qOutput, globalData%q0,&
          & globalData%taggedWriter, globalData%cm5Cont)
    end if
    if (globalData%tWriteDetailedXML) then
      call writeDetailedXml(globalData%runId, globalData%speciesName, globalData%species0,&
          & globalData%pCoord0Out, globalData%tPeriodic, globalData%tHelical, globalData%latVec,&
          & globalData%origin, globalData%tRealHS, globalData%nKPoint, globalData%nSpin,&
          & size(globalData%eigen, dim=1), globalData%nOrb, globalData%kPoint, globalData%kWeight,&
          & globalData%filling, globalData%occNatural)
    end if

    call env%globalTimer%stopTimer(globalTimers%postGeoOpt)

    if (globalData%tPoisson) then
      call poiss_destroy(env)
    end if
  #:if WITH_TRANSPORT
    if (globalData%electronicSolver%iSolver == electronicSolverTypes%GF .or. &
      & globalData%electronicSolver%iSolver == electronicSolverTypes%OnlyTransport) then
      call negf_destroy()
    end if
  #:endif

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
  subroutine processGeometry(env, globalData, iGeoStep, iLatGeoStep, tWriteRestart, tStopScc,&
      & tExitGeoOpt, stat)

    !> Environment settings
    type(TEnvironment), intent(inout) :: env

    !> Global variables
    type(TGlobalData), intent(inout) :: globalData

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
    integer, intent(out), optional :: stat

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

    real(dp), allocatable :: dipoleTmp(:)

    if (globalData%tDipole) then
      allocate(dipoleTmp(3))
    end if

    call env%globalTimer%startTimer(globalTimers%preSccInit)

    if (allocated(globalData%qDepExtPot)) then
      allocate(dQ(globalData%orb%mShell, globalData%nAtom, globalData%nSpin))
    end if

    call globalData%electronicSolver%reset()
    tExitGeoOpt = .false.

    if (globalData%tMD .and. tWriteRestart) then
      call writeMdOut1(globalData%fdMd, mdOut, iGeoStep, globalData%pMDIntegrator)
    end if

    if (globalData%tLatticeChanged) then
      call handleLatticeChange(globalData%latVec, globalData%sccCalc, globalData%tStress,&
          & globalData%extPressure, globalData%cutOff%mCutOff, globalData%dispersion,&
          & globalData%solvation, globalData%cm5Cont, globalData%recVec, globalData%invLatVec,&
          & globalData%cellVol, globalData%recCellVol, globalData%extLatDerivs, globalData%cellVec,&
          & globalData%rCellVec)
    end if

    if (globalData%tCoordsChanged) then
      call handleCoordinateChange(env, globalData%coord0, globalData%latVec, globalData%invLatVec,&
          & globalData%species0, globalData%cutOff, globalData%orb, globalData%tPeriodic,&
          & globalData%tHelical, globalData%sccCalc, globalData%dispersion, globalData%solvation,&
          & globalData%thirdOrd, globalData%rangeSep, globalData%reks, globalData%img2CentCell,&
          & globalData%iCellVec, globalData%neighbourList, globalData%nAllAtom,&
          & globalData%coord0Fold, globalData%coord, globalData%species, globalData%rCellVec,&
          & globalData%nNeighbourSk, globalData%nNeighbourRep, globalData%nNeighbourLC,&
          & globalData%ham, globalData%over, globalData%H0, globalData%rhoPrim,&
          & globalData%iRhoPrim, globalData%iHam, globalData%ERhoPrim, globalData%iSparseStart,&
          & globalData%tPoisson, globalData%cm5Cont, stat)
        @:HANDLE_ERROR(stat)
    end if

    #:if WITH_TRANSPORT
      if (globalData%tNegf) then
        call initNegfStuff(globalData%denseDesc, globalData%transpar, globalData%ginfo,&
            & globalData%neighbourList, globalData%nNeighbourSK, globalData%img2CentCell,&
            & globalData%orb)
      end if
    #:endif

    if (globalData%tSccCalc .and. .not. allocated(globalData%reks) .and. .not.&
        & globalData%tRestartNoSC) then
      call reset(globalData%pChrgMixer, globalData%nMixElements)
    end if

    if (globalData%electronicSolver%isElsiSolver .and. .not. globalData%tLargeDenseMatrices) then
      call globalData%electronicSolver%elsi%updateGeometry(env, globalData%neighbourList,&
          & globalData%nNeighbourSK, globalData%denseDesc%iAtomStart, globalData%iSparseStart,&
          & globalData%img2CentCell)
    end if

    call env%globalTimer%startTimer(globalTimers%sparseH0S)
    select case(globalData%hamiltonianType)
    case default
      call error("Invalid Hamiltonian")
    case(hamiltonianTypes%dftb)
      call buildH0(env, globalData%H0, globalData%skHamCont, globalData%atomEigVal,&
          & globalData%coord, globalData%nNeighbourSk, globalData%neighbourList%iNeighbour,&
          & globalData%species, globalData%iSparseStart, globalData%orb)
      call buildS(env, globalData%over, globalData%skOverCont, globalData%coord,&
          & globalData%nNeighbourSk, globalData%neighbourList%iNeighbour, globalData%species,&
          & globalData%iSparseStart, globalData%orb)
    case(hamiltonianTypes%xtb)
      ! TODO
      call error("xTB calculation currently not supported")
    end select
    call env%globalTimer%stopTimer(globalTimers%sparseH0S)

    if (globalData%tSetFillingTemp) then
      call globalData%temperatureProfile%getTemperature(globalData%tempElec)
    end if

    call globalData%electronicSolver%updateElectronicTemp(globalData%tempElec)

    call calcRepulsiveEnergy(globalData%coord, globalData%species, globalData%img2CentCell,&
        & globalData%nNeighbourRep, globalData%neighbourList, globalData%pRepCont,&
        & globalData%dftbEnergy(globalData%deltaDftb%iDeterminant)%atomRep,&
        & globalData%dftbEnergy(globalData%deltaDftb%iDeterminant)%ERep,&
        & globalData%iAtInCentralRegion)

    if (allocated(globalData%halogenXCorrection)) then
      call globalData%halogenXCorrection%getEnergies(globalData%dftbEnergy(&
          & globalData%deltaDftb%iDeterminant)%atomHalogenX, globalData%coord, globalData%species,&
          & globalData%neighbourList, globalData%img2CentCell)
      globalData%dftbEnergy(globalData%deltaDftb%iDeterminant)%EHalogenX =&
          & sum(globalData%dftbEnergy(globalData%deltaDftb%iDeterminant)%atomHalogenX(&
          & globalData%iAtInCentralRegion))
    end if

    call resetExternalPotentials(globalData%refExtPot, globalData%potential)

    if (globalData%tReadShifts) then
      call readShifts(fShifts, globalData%orb, globalData%nAtom, globalData%nSpin,&
          & globalData%potential%extShell)
    end if

    call setUpExternalElectricField(globalData%tEField, globalData%tTDEField, globalData%tPeriodic,&
        & globalData%EFieldStrength, globalData%EFieldVector, globalData%EFieldOmega,&
        & globalData%EFieldPhase, globalData%neighbourList, globalData%nNeighbourSk,&
        & globalData%iCellVec, globalData%img2CentCell, globalData%cellVec, globalData%deltaT,&
        & iGeoStep, globalData%coord0Fold, globalData%coord, globalData%potential%extAtom(:,1),&
        & globalData%potential%extGrad, globalData%EField, globalData%absEField)

    call mergeExternalPotentials(globalData%orb, globalData%species, globalData%potential)

    ! For non-scc calculations with transport only, jump out of geometry loop
    if (globalData%electronicSolver%iSolver == electronicSolverTypes%OnlyTransport) then
      if (globalData%tWriteDetailedOut) then
        call openDetailedOut(globalData%fdDetailedOut, userOut,&
            & tAppendDetailedOut)
      end if
      ! We need to define hamltonian by adding the potential
      call getSccHamiltonian(globalData%H0, globalData%over, globalData%nNeighbourSK,&
          & globalData%neighbourList, globalData%species, globalData%orb, globalData%iSparseStart,&
          & globalData%img2CentCell, globalData%potential, allocated(globalData%reks),&
          & globalData%ham, globalData%iHam)
      tExitGeoOpt = .true.
      return
    end if

    if (globalData%electronicSolver%iSolver == electronicSolverTypes%pexsi) then
      call globalData%electronicSolver%elsi%initPexsiDeltaVRanges(globalData%tSccCalc,&
          & globalData%potential)
    end if

    if (.not.globalData%tRestartNoSC) then
      call initSccLoop(globalData%tSccCalc, globalData%xlbomdIntegrator, globalData%minSccIter,&
          & globalData%maxSccIter, globalData%sccTol, tConverged, globalData%tNegf, globalData%reks)
    else
      tConverged = .true.
    end if

    call env%globalTimer%stopTimer(globalTimers%preSccInit)

    call env%globalTimer%startTimer(globalTimers%scc)

    REKS_SCC: if (allocated(globalData%reks)) then

      lpSCC_REKS: do iSccIter = 1, globalData%maxSccIter

        if (iSccIter == 1) then
          call getReksInitialSettings(env, globalData%denseDesc, globalData%h0, globalData%over,&
              & globalData%neighbourList, globalData%nNeighbourSK, globalData%iSparseStart,&
              & globalData%img2CentCell, globalData%electronicSolver, iGeoStep,&
              & globalData%HSqrReal, globalData%SSqrReal, globalData%eigvecsReal, globalData%eigen,&
              & globalData%reks)
        end if

        call getDensityMatrixL(env, globalData%denseDesc, globalData%neighbourList,&
            & globalData%nNeighbourSK, globalData%iSparseStart, globalData%img2CentCell,&
            & globalData%orb, globalData%species, globalData%coord, globalData%tHelical,&
            & globalData%eigvecsReal, globalData%parallelKS, globalData%rhoPrim,&
            & globalData%SSqrReal, globalData%rhoSqrReal, globalData%q0, globalData%deltaRhoOutSqr,&
            & globalData%reks)
        call getMullikenPopulationL(env, globalData%denseDesc, globalData%neighbourList,&
            & globalData%nNeighbourSK, globalData%img2CentCell, globalData%iSparseStart,&
            & globalData%orb, globalData%rhoPrim, globalData%over, globalData%iRhoPrim,&
            & globalData%qBlockOut, globalData%qiBlockOut, globalData%reks)

        call getHamiltonianLandEnergyL(env, globalData%denseDesc, globalData%sccCalc,&
            & globalData%orb, globalData%species, globalData%neighbourList,&
            & globalData%nNeighbourSK, globalData%iSparseStart, globalData%img2CentCell,&
            & globalData%H0, globalData%over, globalData%spinW, globalData%cellVol,&
            & globalData%extPressure, globalData%dftbEnergy(1), globalData%q0,&
            & globalData%iAtInCentralRegion, globalData%solvation, globalData%thirdOrd,&
            & globalData%potential, globalData%electrostatics,  globalData%tPoisson,&
            & globalData%tUpload, globalData%shiftPerLUp, globalData%rangeSep,&
            & globalData%nNeighbourLC, globalData%tDualSpinOrbit, globalData%xi,&
            & globalData%tExtField, globalData%isXlbomd, globalData%tDftbU,&
            & globalData%dftbEnergy(1)%TS, globalData%qDepExtPot, globalData%qBlockOut,&
            & globalData%qiBlockOut, globalData%nDftbUFunc, globalData%UJ, globalData%nUJ,&
            & globalData%iUJ, globalData%niUJ, globalData%tFixEf, globalData%Ef,&
            & globalData%rhoPrim, globalData%onSiteElements, globalData%iHam,&
            & globalData%dispersion, globalData%reks)
        call optimizeFONsAndWeights(globalData%eigvecsReal, globalData%filling,&
            & globalData%dftbEnergy(1), globalData%reks)

        call getFockandDiag(env, globalData%denseDesc, globalData%neighbourList,&
            & globalData%nNeighbourSK, globalData%iSparseStart, globalData%img2CentCell,&
            & globalData%eigvecsReal, globalData%electronicSolver, globalData%eigen,&
            & globalData%reks)

        ! Creates (delta) density matrix for averaged state from real eigenvectors.
        call getDensityFromRealEigvecs(env, globalData%denseDesc, globalData%filling(:,1,:),&
            & globalData%neighbourList, globalData%nNeighbourSK, globalData%iSparseStart,&
            & globalData%img2CentCell, globalData%orb, globalData%species,&
            & globalData%denseDesc%iAtomStart, globalData%coord, globalData%tHelical,&
            & globalData%eigVecsReal, globalData%parallelKS, globalData%rhoPrim,&
            & globalData%SSqrReal, globalData%rhoSqrReal, globalData%deltaRhoOutSqr)
        ! For rangeseparated calculations deduct atomic charges from deltaRho
        if (globalData%isRangeSep) then
          call denseSubtractDensityOfAtoms(globalData%q0, globalData%denseDesc%iAtomStart,&
              & globalData%deltaRhoOutSqr)
        end if
        call getMullikenPopulation(globalData%rhoPrim, globalData%over, globalData%orb,&
            & globalData%neighbourList, globalData%nNeighbourSK, globalData%img2CentCell,&
            & globalData%iSparseStart, globalData%qOutput, iRhoPrim=globalData%iRhoPrim,&
            & qBlock=globalData%qBlockOut, qiBlock=globalData%qiBlockOut,&
            & qNetAtom=globalData%qNetAtom)

        ! Check charge convergece and guess new eigenvectors
        tStopScc = hasStopFile(fStopScc)
        if (globalData%isRangeSep) then
          call getReksNextInputDensity(sccErrorQ, globalData%sccTol, tConverged, iSccIter,&
              & globalData%minSccIter, globalData%maxSccIter, iGeoStep, tStopScc,&
              & globalData%eigvecsReal, globalData%deltaRhoOut, globalData%deltaRhoIn,&
              & globalData%deltaRhoDiff, globalData%reks)
        else
          call getReksNextInputCharges(globalData%qInput, globalData%qOutput, globalData%qDiff,&
              & sccErrorQ, globalData%sccTol, tConverged, iSccIter, globalData%minSccIter,&
              & globalData%maxSccIter, iGeoStep, tStopScc, globalData%eigvecsReal, globalData%reks)
        end if

        if (allocated(globalData%dispersion)) then
          call globalData%dispersion%updateOnsiteCharges(globalData%qNetAtom, globalData%orb,&
              & globalData%referenceN0, globalData%species0, tConverged)
          call calcDispersionEnergy(globalData%dispersion, globalData%dftbEnergy(1)%atomDisp,&
              & globalData%dftbEnergy(1)%Edisp, globalData%iAtInCentralRegion)
          call sumEnergies(globalData%dftbEnergy(1))
        end if

        call getSccInfo(iSccIter, globalData%dftbEnergy(1)%Etotal, Eold, diffElec)
        call printReksSccInfo(iSccIter, globalData%dftbEnergy(1)%Etotal, diffElec, sccErrorQ,&
            & globalData%reks)

        if (tConverged .or. tStopScc) then

          call printReksSAInfo(globalData%reks, globalData%dftbEnergy(1)%Etotal)

          call getStateInteraction(env, globalData%denseDesc, globalData%neighbourList,&
              & globalData%nNeighbourSK, globalData%iSparseStart, globalData%img2CentCell,&
              & globalData%coord, globalData%iAtInCentralRegion, globalData%eigvecsReal,&
              & globalData%electronicSolver, globalData%eigen, globalData%qOutput, globalData%q0,&
              & globalData%tDipole, dipoleTmp, globalData%reks)

          call getReksEnProperties(globalData%eigvecsReal, globalData%coord0, globalData%reks)

          if (globalData%tWriteDetailedOut) then
            ! In this routine the correct Etotal is evaluated.
            ! If TargetStateL > 0, certain microstate
            ! is optimized. If not, SSR state is optimized.
            call openDetailedOut(globalData%fdDetailedOut, userOut,&
                & tAppendDetailedOut)
            call writeReksDetailedOut1(globalData%fdDetailedOut, globalData%nGeoSteps, iGeoStep,&
                & globalData%tMD, globalData%tDerivs, globalData%tCoordOpt, globalData%tLatOpt,&
                & iLatGeoStep, iSccIter, globalData%dftbEnergy(1), diffElec, sccErrorQ,&
                & globalData%indMovedAtom, globalData%pCoord0Out, globalData%q0,&
                & globalData%qOutput, globalData%orb, globalData%species,&
                & globalData%tPrintMulliken, globalData%extPressure, globalData%cellVol,&
                & globalData%dftbEnergy(1)%TS, globalData%tAtomicEnergy, globalData%dispersion,&
                & globalData%tPeriodic, globalData%tSccCalc, globalData%invLatVec,&
                & globalData%kPoint, globalData%iAtInCentralRegion, globalData%electronicSolver,&
                & globalData%reks, allocated(globalData%thirdOrd), globalData%isRangeSep)
          end if
          if (globalData%tWriteBandDat) then
            call writeBandOut(bandOut, globalData%eigen, globalData%filling,&
                & globalData%kWeight)
          end if

          exit lpSCC_REKS
        end if
      end do lpSCC_REKS

    else ! not REKS_SCC

      ! Standard spin free or unrestricted DFTB

      lpSCC: do iSccIter = 1, globalData%maxSccIter

        call resetInternalPotentials(globalData%tDualSpinOrbit, globalData%xi, globalData%orb,&
            & globalData%species, globalData%potential)

        if (globalData%tSccCalc) then

          call getChargePerShell(globalData%qInput, globalData%orb, globalData%species,&
              & globalData%chargePerShell)

        #:if WITH_TRANSPORT
          ! Overrides input charges with uploaded contact charges
          if (globalData%tUpload) then
            call overrideContactCharges(globalData%qInput, globalData%chargeUp,&
                & globalData%transpar, globalData%qBlockIn, globalData%blockUp)
          end if
        #:endif

          call addChargePotentials(env, globalData%sccCalc, globalData%qInput, globalData%q0,&
              & globalData%chargePerShell, globalData%orb, globalData%species,&
              & globalData%neighbourList, globalData%img2CentCell, globalData%spinW,&
              & globalData%solvation, globalData%thirdOrd, globalData%potential,&
              & globalData%electrostatics, globalData%tPoisson, globalData%tUpload,&
              & globalData%shiftPerLUp, globalData%dispersion)

          call addBlockChargePotentials(globalData%qBlockIn, globalData%qiBlockIn,&
              & globalData%tDftbU, globalData%tImHam, globalData%species, globalData%orb,&
              & globalData%nDftbUFunc, globalData%UJ, globalData%nUJ, globalData%iUJ,&
              & globalData%niUJ, globalData%potential)

          if (allocated(globalData%onSiteElements) .and. (iSCCIter > 1 .or. globalData%tReadChrg))&
              & then
            call addOnsShift(globalData%potential%intBlock, globalData%potential%iOrbitalBlock,&
                & globalData%qBlockIn, globalData%qiBlockIn, globalData%q0,&
                & globalData%onSiteElements, globalData%species, globalData%orb)
          end if

        end if

        ! All potentials are added up into intBlock
        globalData%potential%intBlock = globalData%potential%intBlock +&
            & globalData%potential%extBlock

        if (allocated(globalData%qDepExtPot)) then
          call getChargePerShell(globalData%qInput, globalData%orb, globalData%species,&
              & dQ, qRef=globalData%q0)
          call globalData%qDepExtPot%addPotential(sum(dQ(:,:,1), dim=1), dQ(:,:,1), globalData%orb,&
              & globalData%species, globalData%potential%intBlock)
        end if

        if (globalData%electronicSolver%iSolver == electronicSolverTypes%pexsi .and.&
            & globalData%tSccCalc) then
          call globalData%electronicSolver%elsi%updatePexsiDeltaVRanges(globalData%potential)
        end if

        call getSccHamiltonian(globalData%H0, globalData%over, globalData%nNeighbourSK,&
            & globalData%neighbourList, globalData%species, globalData%orb,&
            & globalData%iSparseStart, globalData%img2CentCell, globalData%potential,&
            & allocated(globalData%reks), globalData%ham, globalData%iHam)

        if (globalData%tWriteRealHS .or. globalData%tWriteHS .and.&
            & any(globalData%electronicSolver%iSolver ==&
            & [electronicSolverTypes%qr, electronicSolverTypes%divideandconquer,&
            & electronicSolverTypes%relativelyrobust, electronicSolverTypes%magma_gvd])) then
          call writeHSAndStop(env, globalData%tWriteHS, globalData%tWriteRealHS,&
              & globalData%tRealHS, globalData%over, globalData%neighbourList,&
              & globalData%nNeighbourSK, globalData%denseDesc%iAtomStart, globalData%iSparseStart,&
              & globalData%img2CentCell, globalData%kPoint, globalData%iCellVec,&
              & globalData%cellVec, globalData%ham, globalData%iHam)
        end if

        call convertToUpDownRepr(globalData%ham, globalData%iHam)

        call getDensity(env, iSccIter, globalData%denseDesc, globalData%ham, globalData%over,&
            & globalData%neighbourList, globalData%nNeighbourSk, globalData%iSparseStart,&
            & globalData%img2CentCell, globalData%iCellVec, globalData%cellVec, globalData%kPoint,&
            & globalData%kWeight, globalData%orb, globalData%tHelical, globalData%coord,&
            & globalData%species, globalData%electronicSolver, globalData%tRealHS,&
            & globalData%tSpinSharedEf, globalData%tSpinOrbit, globalData%tDualSpinOrbit,&
            & globalData%tFillKSep, globalData%tFixEf, globalData%tMulliken, globalData%iDistribFn,&
            & globalData%tempElec, globalData%nEl, globalData%parallelKS, globalData%Ef,&
            & globalData%mu, globalData%dftbEnergy(globalData%deltaDftb%iDeterminant),&
            & globalData%rangeSep, globalData%eigen, globalData%filling, globalData%rhoPrim,&
            & globalData%iHam, globalData%xi, globalData%orbitalL, globalData%HSqrReal,&
            & globalData%SSqrReal, globalData%eigvecsReal, globalData%iRhoPrim,&
            & globalData%HSqrCplx, globalData%SSqrCplx, globalData%eigvecsCplx,&
            & globalData%rhoSqrReal, globalData%deltaRhoInSqr, globalData%deltaRhoOutSqr,&
            & globalData%qOutput, globalData%nNeighbourLC, globalData%tLargeDenseMatrices,&
            & globalData%deltaDftb)

        !> For rangeseparated calculations deduct atomic charges from deltaRho
        if (globalData%isRangeSep) then
          select case(globalData%nSpin)
          case(2)
            do iSpin = 1, 2
              call denseSubtractDensityOfAtoms(globalData%q0, globalData%denseDesc%iAtomStart,&
                  & globalData%deltaRhoOutSqr, iSpin)
            end do
          case(1)
            call denseSubtractDensityOfAtoms(globalData%q0, globalData%denseDesc%iAtomStart,&
                & globalData%deltaRhoOutSqr)
          case default
            call error("Range separation not implemented for non-colinear spin")
          end select
        end if

        if (globalData%tWriteBandDat .and. globalData%deltaDftb%nDeterminant() == 1) then
          call writeBandOut(bandOut, globalData%eigen, globalData%filling,&
              & globalData%kWeight)
        end if

        if (globalData%tMulliken) then
          call getMullikenPopulation(globalData%rhoPrim, globalData%over, globalData%orb,&
              & globalData%neighbourList, globalData%nNeighbourSk, globalData%img2CentCell,&
              & globalData%iSparseStart, globalData%qOutput, iRhoPrim=globalData%iRhoPrim,&
              & qBlock=globalData%qBlockOut, qiBlock=globalData%qiBlockOut,&
              & qNetAtom=globalData%qNetAtom)
        end if

      #:if WITH_TRANSPORT
        ! Override charges with uploaded contact charges
        if (globalData%tUpload) then
          call overrideContactCharges(globalData%qOutput, globalData%chargeUp, globalData%transpar,&
              & globalData%qBlockIn, globalData%blockUp)
        end if
      #:endif

        ! For non-dual spin-orbit orbitalL is determined during getDensity() call above
        if (globalData%tDualSpinOrbit) then
          call getLDual(globalData%orbitalL, globalData%qiBlockOut, globalData%orb,&
              & globalData%species)
        end if

        ! Note: if XLBOMD is active, potential created with input charges is needed later,
        ! therefore it should not be overwritten here.
        if (globalData%tSccCalc .and. .not. globalData%isXlbomd) then
          call resetInternalPotentials(globalData%tDualSpinOrbit, globalData%xi, globalData%orb,&
              & globalData%species, globalData%potential)
          call getChargePerShell(globalData%qOutput, globalData%orb, globalData%species,&
              & globalData%chargePerShell)

          call addChargePotentials(env, globalData%sccCalc, globalData%qOutput, globalData%q0,&
              & globalData%chargePerShell, globalData%orb, globalData%species,&
              & globalData%neighbourList, globalData%img2CentCell, globalData%spinW,&
              & globalData%solvation, globalData%thirdOrd, globalData%potential,&
              & globalData%electrostatics, globalData%tPoissonTwice, globalData%tUpload,&
              & globalData%shiftPerLUp, globalData%dispersion)

          call addBlockChargePotentials(globalData%qBlockOut, globalData%qiBlockOut,&
              & globalData%tDftbU, globalData%tImHam, globalData%species, globalData%orb,&
              & globalData%nDftbUFunc, globalData%UJ, globalData%nUJ, globalData%iUJ,&
              & globalData%niUJ, globalData%potential)

          if (allocated(globalData%onSiteElements)) then
            call addOnsShift(globalData%potential%intBlock, globalData%potential%iOrbitalBlock,&
                & globalData%qBlockOut, globalData%qiBlockOut, globalData%q0,&
                & globalData%onSiteElements, globalData%species, globalData%orb)
          end if

          globalData%potential%intBlock = globalData%potential%intBlock +&
              & globalData%potential%extBlock
        end if

        if (allocated(globalData%qDepExtPot)) then
          call getChargePerShell(globalData%qOutput, globalData%orb, globalData%species,&
              & dQ, qRef=globalData%q0)
          call globalData%qDepExtPot%addPotential(sum(dQ(:,:,1), dim=1), dQ(:,:,1), globalData%orb,&
              & globalData%species, globalData%potential%intBlock)
        end if

        call calcEnergies(globalData%sccCalc, globalData%qOutput, globalData%q0,&
            & globalData%chargePerShell, globalData%species, globalData%tExtField,&
            & globalData%isXlbomd, globalData%tDftbU, globalData%tDualSpinOrbit,&
            & globalData%rhoPrim, globalData%H0, globalData%orb, globalData%neighbourList,&
            & globalData%nNeighbourSk, globalData%img2CentCell, globalData%iSparseStart,&
            & globalData%cellVol, globalData%extPressure,&
            & globalData%dftbEnergy(globalData%deltaDftb%iDeterminant)%TS, globalData%potential,&
            & globalData%dftbEnergy(globalData%deltaDftb%iDeterminant), globalData%thirdOrd,&
            & globalData%solvation, globalData%rangeSep, globalData%reks, globalData%qDepExtPot,&
            & globalData%qBlockOut, globalData%qiBlockOut, globalData%nDftbUFunc, globalData%UJ,&
            & globalData%nUJ, globalData%iUJ, globalData%niUJ, globalData%xi,&
            & globalData%iAtInCentralRegion, globalData%tFixEf, globalData%Ef,&
            & globalData%onSiteElements)

        tStopScc = hasStopFile(fStopScc)

        ! Mix charges Input/Output
        if (globalData%tSccCalc) then
          if(.not. globalData%isRangeSep) then
            call getNextInputCharges(env, globalData%pChrgMixer, globalData%qOutput,&
                & globalData%qOutRed, globalData%orb, globalData%nIneqOrb, globalData%iEqOrbitals,&
                & iGeoStep, iSccIter, globalData%minSccIter, globalData%maxSccIter,&
                & globalData%sccTol, tStopScc, globalData%tMixBlockCharges, globalData%tReadChrg,&
                & globalData%qInput, globalData%qInpRed, sccErrorQ, tConverged,&
                & globalData%qBlockOut, globalData%iEqBlockDftbU, globalData%qBlockIn,&
                & globalData%qiBlockOut, globalData%iEqBlockDftbULS, globalData%species0,&
                & globalData%nUJ, globalData%iUJ, globalData%niUJ, globalData%qiBlockIn,&
                & globalData%iEqBlockOnSite, globalData%iEqBlockOnSiteLS)
          else
            call getNextInputDensity(globalData%SSqrReal, globalData%over,&
                & globalData%neighbourList, globalData%nNeighbourSK,&
                & globalData%denseDesc%iAtomStart, globalData%iSparseStart,&
                & globalData%img2CentCell, globalData%pChrgMixer, globalData%qOutput,&
                & globalData%orb, globalData%tHelical, globalData%species, globalData%coord,&
                & iGeoStep, iSccIter, globalData%minSccIter, globalData%maxSccIter,&
                & globalData%sccTol, tStopScc, globalData%tReadChrg, globalData%q0,&
                & globalData%qInput, sccErrorQ, tConverged, globalData%deltaRhoOut,&
                & globalData%deltaRhoIn, globalData%deltaRhoDiff, globalData%qBlockIn,&
                & globalData%qBlockOut)
          end if

          call getSccInfo(iSccIter, globalData%dftbEnergy(globalData%deltaDftb%iDeterminant)%Eelec,&
              & Eold, diffElec)
          if (globalData%tNegf) then
            call printSccHeader()
          end if
          call printSccInfo(globalData%tDftbU, iSccIter,&
              & globalData%dftbEnergy(globalData%deltaDftb%iDeterminant)%Eelec, diffElec, sccErrorQ)

          if (globalData%tNegf) then
            call printBlankLine()
          end if

          tWriteSccRestart = env%tGlobalLead .and. needsSccRestartWriting(globalData%restartFreq,&
              & iGeoStep, iSccIter, globalData%minSccIter, globalData%maxSccIter, globalData%tMd,&
              & globalData%isGeoOpt, globalData%tDerivs, tConverged, globalData%tReadChrg, tStopScc)
          if (tWriteSccRestart) then
            call writeCharges(fCharges, globalData%tWriteChrgAscii, globalData%orb,&
                & globalData%qInput, globalData%qBlockIn, globalData%qiBlockIn,&
                & globalData%deltaRhoIn)
          end if
        end if

        if (allocated(globalData%dispersion)) then
          call globalData%dispersion%updateOnsiteCharges(globalData%qNetAtom, globalData%orb,&
              & globalData%referenceN0, globalData%species0, tConverged)
          call calcDispersionEnergy(globalData%dispersion,&
              & globalData%dftbEnergy(globalData%deltaDftb%iDeterminant)%atomDisp,&
              & globalData%dftbEnergy(globalData%deltaDftb%iDeterminant)%Edisp,&
              & globalData%iAtInCentralRegion)
        end if

        call sumEnergies(globalData%dftbEnergy(globalData%deltaDftb%iDeterminant))

        if (globalData%tWriteDetailedOut .and. globalData%deltaDftb%nDeterminant() == 1) then
          call openDetailedOut(globalData%fdDetailedOut, userOut,&
              & tAppendDetailedOut)
          call writeDetailedOut1(globalData%fdDetailedOut, globalData%iDistribFn,&
              & globalData%nGeoSteps, iGeoStep, globalData%tMD, globalData%tDerivs,&
              & globalData%tCoordOpt, globalData%tLatOpt, iLatGeoStep, iSccIter,&
              & globalData%dftbEnergy(globalData%deltaDftb%iDeterminant), diffElec, sccErrorQ,&
              & globalData%indMovedAtom, globalData%pCoord0Out, globalData%tPeriodic,&
              & globalData%tSccCalc, globalData%tNegf, globalData%invLatVec, globalData%kPoint)
          call writeDetailedOut2(globalData%fdDetailedOut, globalData%q0, globalData%qInput,&
              & globalData%qOutput, globalData%orb, globalData%species, globalData%tDFTBU,&
              & globalData%tImHam .or. globalData%tSpinOrbit, globalData%tPrintMulliken,&
              & globalData%orbitalL, globalData%qBlockOut, globalData%nSpin,&
              & allocated(globalData%onSiteElements), globalData%iAtInCentralRegion,&
              & globalData%cm5Cont, globalData%qNetAtom)
          call writeDetailedOut3(globalData%fdDetailedOut, globalData%qInput, globalData%qOutput,&
              & globalData%dftbEnergy(globalData%deltaDftb%iDeterminant), globalData%species,&
              & globalData%tDFTBU, globalData%tPrintMulliken, globalData%Ef,&
              & globalData%extPressure, globalData%cellVol, globalData%tAtomicEnergy,&
              & globalData%dispersion, globalData%tEField, globalData%tPeriodic, globalData%nSpin,&
              & globalData%tSpin, globalData%tSpinOrbit, globalData%tSccCalc,&
              & allocated(globalData%onSiteElements), globalData%tNegf,&
              & globalData%iAtInCentralRegion, globalData%electronicSolver,&
              & allocated(globalData%halogenXCorrection), globalData%isRangeSep,&
              & allocated(globalData%thirdOrd), allocated(globalData%solvation))
        end if

        if (tConverged .or. tStopScc) then
          exit lpSCC
        end if

      end do lpSCC

    end if REKS_SCC

    if (allocated(globalData%dispersion)) then
      ! If we get to this point for a dispersion model, if it is charge dependent it may require
      ! evaluation post-hoc if SCC was not achieved but the input settings are to proceed with
      ! non-converged SCC.
      call globalData%dispersion%updateOnsiteCharges(globalData%qNetAtom, globalData%orb,&
          & globalData%referenceN0, globalData%species0, tConverged .or. .not.&
          & globalData%isSccConvRequired)
      call calcDispersionEnergy(globalData%dispersion,&
          & globalData%dftbEnergy(globalData%deltaDftb%iDeterminant)%atomDisp,&
          & globalData%dftbEnergy(globalData%deltaDftb%iDeterminant)%Edisp,&
          & globalData%iAtInCentralRegion)
      call sumEnergies(globalData%dftbEnergy(globalData%deltaDftb%iDeterminant))

      if (globalData%tWriteDetailedOut .and. globalData%deltaDftb%nDeterminant() == 1) then
        close(globalData%fdDetailedOut)
        call openDetailedOut(globalData%fdDetailedOut, userOut,&
            & tAppendDetailedOut)
        if (allocated(globalData%reks)) then
          call writeReksDetailedOut1(globalData%fdDetailedOut, globalData%nGeoSteps, iGeoStep,&
              & globalData%tMD, globalData%tDerivs, globalData%tCoordOpt, globalData%tLatOpt ,&
              & iLatGeoStep, iSccIter, globalData%dftbEnergy(1), diffElec, sccErrorQ,&
              & globalData%indMovedAtom, globalData%pCoord0Out, globalData%q0, globalData%qOutput,&
              & globalData%orb, globalData%species, globalData%tPrintMulliken,&
              & globalData%extPressure, globalData%cellVol, globalData%dftbEnergy(1)%TS,&
              & globalData%tAtomicEnergy, globalData%dispersion, globalData%tPeriodic,&
              & globalData%tSccCalc, globalData%invLatVec, globalData%kPoint,&
              & globalData%iAtInCentralRegion, globalData%electronicSolver, globalData%reks,&
              & allocated(globalData%thirdOrd), globalData%isRangeSep)
        else
          call writeDetailedOut1(globalData%fdDetailedOut, globalData%iDistribFn,&
              & globalData%nGeoSteps, iGeoStep, globalData%tMD, globalData%tDerivs,&
              & globalData%tCoordOpt, globalData%tLatOpt, iLatGeoStep, iSccIter,&
              & globalData%dftbEnergy(globalData%deltaDftb%iDeterminant), diffElec, sccErrorQ,&
              & globalData%indMovedAtom, globalData%pCoord0Out, globalData%tPeriodic,&
              & globalData%tSccCalc, globalData%tNegf, globalData%invLatVec, globalData%kPoint)
          call writeDetailedOut2(globalData%fdDetailedOut, globalData%q0, globalData%qInput,&
              & globalData%qOutput, globalData%orb, globalData%species, globalData%tDFTBU,&
              & globalData%tImHam.or.globalData%tSpinOrbit, globalData%tPrintMulliken,&
              & globalData%orbitalL, globalData%qBlockOut, globalData%nSpin,&
              & allocated(globalData%onSiteElements), globalData%iAtInCentralRegion,&
              & globalData%cm5Cont, globalData%qNetAtom)
          call writeDetailedOut3(globalData%fdDetailedOut, globalData%qInput, globalData%qOutput,&
              & globalData%dftbEnergy(globalData%deltaDftb%iDeterminant), globalData%species,&
              & globalData%tDFTBU, globalData%tPrintMulliken, globalData%Ef,&
              & globalData%extPressure, globalData%cellVol, globalData%tAtomicEnergy,&
              & globalData%dispersion, globalData%tEField, globalData%tPeriodic, globalData%nSpin,&
              & globalData%tSpin, globalData%tSpinOrbit, globalData%tSccCalc,&
              & allocated(globalData%onSiteElements), globalData%tNegf,&
              & globalData%iAtInCentralRegion, globalData%electronicSolver,&
              & allocated(globalData%halogenXCorrection), globalData%isRangeSep,&
              & allocated(globalData%thirdOrd), allocated(globalData%solvation))
        end if
      end if

    end if

    call env%globalTimer%stopTimer(globalTimers%scc)

    if (globalData%tPoisson) then
      call poiss_savepotential(env)
    end if

    call env%globalTimer%startTimer(globalTimers%postSCC)

    if (globalData%isLinResp) then
      if (.not. globalData%isRS_LinResp) then
        call calculateLinRespExcitations(env, globalData%linearResponse, globalData%parallelKS,&
            & globalData%sccCalc, globalData%qOutput, globalData%q0, globalData%over,&
            & globalData%eigvecsReal, globalData%eigen(:,1,:), globalData%filling(:,1,:),&
            & globalData%coord, globalData%species, globalData%speciesName, globalData%orb,&
            & globalData%skHamCont, globalData%skOverCont, autotestTag,&
            & globalData%taggedWriter, globalData%runId, globalData%neighbourList,&
            & globalData%nNeighbourSK, globalData%denseDesc, globalData%iSparseStart,&
            & globalData%img2CentCell, globalData%tWriteAutotest, globalData%tCasidaForces,&
            & globalData%tLinRespZVect, globalData%tPrintExcitedEigvecs,&
            & globalData%tPrintEigvecsTxt, globalData%nonSccDeriv, globalData%dftbEnergy(1),&
            & globalData%energiesCasida, globalData%SSqrReal, globalData%rhoSqrReal,&
            & globalData%excitedDerivs, globalData%dQAtomEx, globalData%occNatural)
      else
        call calculateLinRespExcitations_RS(env, globalData%linearResponse, globalData%parallelKS,&
            & globalData%sccCalc, globalData%qOutput, globalData%q0, globalData%over,&
            & globalData%eigvecsReal, globalData%eigen(:,1,:), globalData%filling(:,1,:),&
            & globalData%coord0, globalData%species, globalData%speciesName, globalData%orb,&
            & globalData%skHamCont, globalData%skOverCont, autotestTag,&
            & globalData%taggedWriter, globalData%runId, globalData%neighbourList,&
            & globalData%nNeighbourSK, globalData%denseDesc, globalData%iSparseStart,&
            & globalData%img2CentCell, globalData%tWriteAutotest, globalData%tCasidaForces,&
            & globalData%tLinRespZVect, globalData%tPrintExcitedEigvecs,&
            & globalData%tPrintEigvecsTxt, globalData%nonSccDeriv, globalData%dftbEnergy(1),&
            & globalData%energiesCasida, globalData%SSqrReal, globalData%deltaRhoOutSqr,&
            & globalData%excitedDerivs, globalData%dQAtomEx, globalData%occNatural,&
            & globalData%rangeSep)
      end if
    end if

    if (allocated(globalData%ppRPA)) then
      call unpackHS(globalData%SSqrReal, globalData%over, globalData%neighbourList%iNeighbour,&
          & globalData%nNeighbourSK, globalData%denseDesc%iAtomStart, globalData%iSparseStart,&
          & globalData%img2CentCell)
      call blockSymmetrizeHS(globalData%SSqrReal, globalData%denseDesc%iAtomStart)
      if (withMpi) then
        call error("pp-RPA calc. does not work with MPI yet")
      end if
      call ppRPAenergies(globalData%ppRPA, globalData%denseDesc, globalData%eigvecsReal,&
          & globalData%eigen(:,1,:), globalData%sccCalc, globalData%SSqrReal, globalData%species0,&
          & globalData%nEl(1), globalData%neighbourList%iNeighbour, globalData%img2CentCell,&
          & globalData%orb, globalData%tWriteAutotest, autotestTag,&
          & globalData%taggedWriter)
    end if

    if (globalData%isXlbomd) then
      call getXlbomdCharges(globalData%xlbomdIntegrator, globalData%qOutRed, globalData%pChrgMixer,&
          & globalData%orb, globalData%nIneqOrb, globalData%iEqOrbitals, globalData%qInput,&
          & globalData%qInpRed, globalData%iEqBlockDftbU, globalData%qBlockIn, globalData%species0,&
          & globalData%nUJ, globalData%iUJ, globalData%niUJ, globalData%iEqBlockDftbuLs,&
          & globalData%qiBlockIn, globalData%iEqBlockOnSite, globalData%iEqBlockOnSiteLS)
    end if

    if (globalData%tDipole .and. .not. allocated(globalData%reks) .and. .not.&
        & globalData%tRestartNoSC) then
      call getDipoleMoment(globalData%qOutput, globalData%q0, globalData%coord,&
          & globalData%dipoleMoment(:,globalData%deltaDftb%iDeterminant),&
          & globalData%iAtInCentralRegion)
    #:block DEBUG_CODE
      call checkDipoleViaHellmannFeynman(globalData%rhoPrim, globalData%q0, globalData%coord0,&
          & globalData%over, globalData%orb, globalData%neighbourList, globalData%nNeighbourSk,&
          & globalData%species, globalData%iSparseStart, globalData%img2CentCell)
    #:endblock DEBUG_CODE
    end if

    call env%globalTimer%startTimer(globalTimers%eigvecWriting)

    if (globalData%tPrintEigVecs) then
      call writeEigenvectors(env, globalData%runId, globalData%neighbourList,&
          & globalData%nNeighbourSk, globalData%cellVec, globalData%iCellVec, globalData%denseDesc,&
          & globalData%iSparseStart, globalData%img2CentCell, globalData%species,&
          & globalData%speciesName, globalData%orb, globalData%kPoint, globalData%over,&
          & globalData%parallelKS, globalData%tPrintEigvecsTxt, globalData%eigvecsReal,&
          & globalData%SSqrReal, globalData%eigvecsCplx, globalData%SSqrCplx)
    end if

    if (globalData%tProjEigenvecs) then
      call writeProjectedEigenvectors(env, globalData%regionLabels, globalData%eigen,&
          & globalData%neighbourList, globalData%nNeighbourSk, globalData%cellVec,&
          & globalData%iCellVec, globalData%denseDesc, globalData%iSparseStart,&
          & globalData%img2CentCell, globalData%orb, globalData%over, globalData%kPoint,&
          & globalData%kWeight, globalData%iOrbRegion, globalData%parallelKS,&
          & globalData%eigvecsReal, globalData%SSqrReal, globalData%eigvecsCplx,&
          & globalData%SSqrCplx)
    end if
    call env%globalTimer%stopTimer(globalTimers%eigvecWriting)

    ! MD geometry files are written only later, once velocities for the current geometry are known
    if (globalData%isGeoOpt .and. tWriteRestart) then
      if (.not. (globalData%deltaDftb%isSpinPurify .and.&
          & globalData%deltaDftb%iDeterminant == determinants%triplet)) then
        call writeCurrentGeometry(globalData%geoOutFile, globalData%pCoord0Out, globalData%tLatOpt,&
            & globalData%tMd, globalData%tAppendGeo, globalData%tFracCoord, globalData%tPeriodic,&
            & globalData%tHelical, globalData%tPrintMulliken, globalData%species0,&
            & globalData%speciesName, globalData%latVec, globalData%origin, iGeoStep, iLatGeoStep,&
            & globalData%nSpin, globalData%qOutput, globalData%velocities)
      endif
    end if

    if (globalData%tForces) then
      call env%globalTimer%startTimer(globalTimers%forceCalc)
      if (allocated(globalData%reks)) then
        call getReksGradients(env, globalData%denseDesc, globalData%sccCalc, globalData%rangeSep,&
            & globalData%dispersion, globalData%neighbourList, globalData%nNeighbourSK,&
            & globalData%nNeighbourRep, globalData%iSparseStart, globalData%img2CentCell,&
            & globalData%orb, globalData%nonSccDeriv, globalData%skHamCont, globalData%skOverCont,&
            & globalData%pRepCont, globalData%coord, globalData%coord0, globalData%species,&
            & globalData%q0, globalData%eigvecsReal, globalData%chrgForces, globalData%over,&
            & globalData%spinW, globalData%derivs, globalData%tWriteAutotest,&
            & autotestTag, globalData%taggedWriter, globalData%reks)
        call getReksGradProperties(env, globalData%denseDesc, globalData%neighbourList,&
            & globalData%nNeighbourSK, globalData%iSparseStart, globalData%img2CentCell,&
            & globalData%eigvecsReal, globalData%orb, globalData%iAtInCentralRegion,&
            & globalData%coord, globalData%coord0, globalData%over, globalData%rhoPrim,&
            & globalData%qOutput, globalData%q0, globalData%tDipole, dipoleTmp,&
            & globalData%chrgForces, globalData%reks)
      else
        call env%globalTimer%startTimer(globalTimers%energyDensityMatrix)
        call getEnergyWeightedDensity(env, globalData%electronicSolver, globalData%denseDesc,&
            & globalData%forceType, globalData%filling, globalData%eigen, globalData%kPoint,&
            & globalData%kWeight, globalData%neighbourList, globalData%nNeighbourSK, globalData%orb,&
            & globalData%iSparseStart, globalData%img2CentCell, globalData%iCellVec,&
            & globalData%cellVec, globalData%tRealHS, globalData%ham, globalData%over,&
            & globalData%parallelKS, globalData%tHelical, globalData%species, globalData%coord,&
            & iSccIter, globalData%mu, globalData%ERhoPrim, globalData%eigvecsReal,&
            & globalData%SSqrReal, globalData%eigvecsCplx, globalData%SSqrCplx)
        call env%globalTimer%stopTimer(globalTimers%energyDensityMatrix)
        call getGradients(env, globalData%sccCalc, globalData%tExtField, globalData%isXlbomd,&
            & globalData%nonSccDeriv, globalData%EField, globalData%rhoPrim, globalData%ERhoPrim,&
            & globalData%qOutput, globalData%q0, globalData%skHamCont, globalData%skOverCont,&
            & globalData%pRepCont, globalData%neighbourList, globalData%nNeighbourSk,&
            & globalData%nNeighbourRep, globalData%species, globalData%img2CentCell,&
            & globalData%iSparseStart, globalData%orb, globalData%potential, globalData%coord,&
            & globalData%derivs, globalData%groundDerivs, globalData%tripletderivs,&
            & globalData%mixedderivs, globalData%iRhoPrim, globalData%thirdOrd,&
            & globalData%solvation, globalData%qDepExtPot, globalData%chrgForces,&
            & globalData%dispersion, globalData%rangeSep, globalData%SSqrReal, globalData%over,&
            & globalData%denseDesc, globalData%deltaRhoOutSqr, globalData%tPoisson,&
            & globalData%halogenXCorrection, globalData%tHelical, globalData%coord0,&
            & globalData%deltaDftb)

        if (globalData%tCasidaForces) then
          globalData%derivs(:,:) = globalData%derivs + globalData%excitedDerivs
        end if
      end if

      call env%globalTimer%stopTimer(globalTimers%forceCalc)

      call updateDerivsByPlumed(env, globalData%plumedCalc, globalData%nAtom, iGeoStep,&
          & globalData%derivs, globalData%dftbEnergy(globalData%deltaDftb%iDeterminant)%EMermin,&
          & globalData%coord0, globalData%mass, globalData%tPeriodic, globalData%latVec)

      if (globalData%tStress) then
        call env%globalTimer%startTimer(globalTimers%stressCalc)
        if (allocated(globalData%reks)) then
          call getReksStress(env, globalData%denseDesc, globalData%sccCalc, globalData%nonSccDeriv,&
              & globalData%skHamCont, globalData%skOverCont, globalData%pRepCont,&
              & globalData%neighbourList, globalData%nNeighbourSk, globalData%nNeighbourRep,&
              & globalData%species, globalData%img2CentCell, globalData%iSparseStart,&
              & globalData%orb, globalData%dispersion, globalData%coord, globalData%q0,&
              & globalData%invLatVec, globalData%cellVol, globalData%totalStress,&
              & globalData%totalLatDeriv, globalData%intPressure, globalData%reks)
        else
          call getStress(env, globalData%sccCalc, globalData%thirdOrd, globalData%tExtField,&
              & globalData%nonSccDeriv, globalData%rhoPrim, globalData%ERhoPrim,&
              & globalData%qOutput, globalData%q0, globalData%skHamCont, globalData%skOverCont,&
              & globalData%pRepCont, globalData%neighbourList, globalData%nNeighbourSk,&
              & globalData%nNeighbourRep, globalData%species, globalData%img2CentCell,&
              & globalData%iSparseStart, globalData%orb, globalData%potential, globalData%coord,&
              & globalData%latVec, globalData%invLatVec, globalData%cellVol, globalData%coord0,&
              & globalData%totalStress, globalData%totalLatDeriv, globalData%intPressure,&
              & globalData%iRhoPrim, globalData%solvation, globalData%dispersion,&
              & globalData%halogenXCorrection, globalData%deltaDftb, globalData%tripletStress,&
              & globalData%mixedStress)
        end if
        call env%globalTimer%stopTimer(globalTimers%stressCalc)

      end if

    end if

    if (globalData%tWriteDetailedOut  .and. globalData%deltaDftb%nDeterminant() == 1) then
      call writeDetailedOut4(globalData%fdDetailedOut, globalData%tSccCalc, tConverged,&
          & globalData%isXlbomd, globalData%isLinResp, globalData%isGeoOpt, globalData%tMD,&
          & globalData%tPrintForces, globalData%tStress, globalData%tPeriodic,&
          & globalData%dftbEnergy(globalData%deltaDftb%iDeterminant), globalData%totalStress,&
          & globalData%totalLatDeriv, globalData%derivs, globalData%chrgForces,&
          & globalData%indMovedAtom, globalData%cellVol, globalData%intPressure,&
          & globalData%geoOutFile, globalData%iAtInCentralRegion)
    end if

    if (allocated(globalData%dispersion)) then
      if (.not.globalData%dispersion%energyAvailable()) then
        call warning("Dispersion contributions are not included in the energy")
      end if
    end if

    if (globalData%tSccCalc .and. .not. globalData%isXlbomd .and. .not. tConverged .and. .not.&
        & globalData%tRestartNoSC) then
      call warning("SCC is NOT converged, maximal SCC iterations exceeded")
      if (globalData%isSccConvRequired) then
        call env%shutdown()
      end if
    end if

    if (globalData%tSccCalc .and. allocated(globalData%esp) .and. (.not. (globalData%isGeoOpt .or.&
        & globalData%tMD) .or. needsRestartWriting(globalData%isGeoOpt, globalData%tMd, iGeoStep,&
        & globalData%nGeoSteps, globalData%restartFreq))) then
      call globalData%esp%evaluate(env, globalData%sccCalc, globalData%EField)
      call writeEsp(globalData%esp, env, iGeoStep, globalData%nGeoSteps)
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


  !> Next geometry step from driver
  subroutine getNextGeometry(env, globalData, iGeoStep, tWriteRestart, constrLatDerivs,&
      & tCoordStep, tGeomEnd, tStopDriver, iLatGeoStep, tempIon, tExitGeoOpt)

    !> Environment settings
    type(TEnvironment), intent(inout) :: env

    !> Global variables
    type(TGlobalData), intent(inout) :: globalData

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
    globalData%tCoordsChanged = .false.
    globalData%tLatticeChanged = .false.

    tExitGeoOpt = .false.

    if (globalData%tDerivs) then
      call getNextDerivStep(globalData%derivDriver, globalData%derivs, globalData%indMovedAtom,&
          & globalData%coord0, tGeomEnd)
      if (tGeomEnd) then
        call env%globalTimer%stopTimer(globalTimers%postSCC)
        tExitGeoOpt = .true.
        return
      end if
      globalData%tCoordsChanged = .true.
    else if (globalData%isGeoOpt) then
      globalData%tCoordsChanged = .true.
      if (tCoordStep) then
        call getNextCoordinateOptStep(globalData%pGeoCoordOpt, globalData%dftbEnergy(&
            & globalData%deltaDftb%iFinal), globalData%derivs, globalData%indMovedAtom,&
            & globalData%coord0, diffGeo, tCoordEnd, .not. globalData%tCasidaForces)
        if (.not. globalData%tLatOpt) then
          tGeomEnd = tCoordEnd
        end if
        if (.not. tGeomEnd .and. tCoordEnd .and. diffGeo < tolSameDist) then
          tCoordStep = .false.
        end if
      else
        call getNextLatticeOptStep(globalData%pGeoLatOpt, globalData%dftbEnergy(&
            & globalData%deltaDftb%iDeterminant), constrLatDerivs, globalData%origLatVec,&
            & globalData%tLatOptFixAng, globalData%tLatOptFixLen, globalData%tLatOptIsotropic,&
            & globalData%indMovedAtom, globalData%latVec, globalData%coord0, diffGeo, tGeomEnd)
        iLatGeoStep = iLatGeoStep + 1
        globalData%tLatticeChanged = .true.
        if (.not. tGeomEnd .and. globalData%tCoordOpt) then
          tCoordStep = .true.
          call reset(globalData%pGeoCoordOpt, reshape(globalData%coord0(&
              & :, globalData%indMovedAtom), [globalData%nMovedCoord]))
        end if
      end if
      if (tGeomEnd .and. diffGeo < tolSameDist) then
        call env%globalTimer%stopTimer(globalTimers%postSCC)
        tExitGeoOpt = .true.
        return
      end if
    else if (globalData%tMD) then
      ! New MD coordinates saved in a temporary variable, as writeCurrentGeometry() below
      ! needs the old ones to write out consistent geometries and velocities.
      globalData%newCoords(:,:) = globalData%coord0
      call getNextMdStep(globalData%pMdIntegrator, globalData%pMdFrame,&
          & globalData%temperatureProfile, globalData%derivs, globalData%movedMass,&
          & globalData%mass, globalData%cellVol, globalData%invLatVec, globalData%species0,&
          & globalData%indMovedAtom, globalData%tStress, globalData%tBarostat,&
          & globalData%dftbEnergy(globalData%deltaDftb%iDeterminant), globalData%newCoords,&
          & globalData%latVec, globalData%intPressure, globalData%totalStress,&
          & globalData%totalLatDeriv, globalData%velocities, tempIon)
      globalData%tCoordsChanged = .true.
      globalData%tLatticeChanged = globalData%tBarostat
      call printMdInfo(globalData%tSetFillingTemp, globalData%tEField, globalData%tPeriodic,&
          & globalData%tempElec, globalData%absEField, tempIon, globalData%intPressure,&
          & globalData%extPressure, globalData%dftbEnergy(globalData%deltaDftb%iDeterminant))
      if (tWriteRestart) then
        if (globalData%tPeriodic) then
          globalData%cellVol = abs(determinant33(globalData%latVec))
          globalData%dftbEnergy(globalData%deltaDftb%iDeterminant)%EGibbs =&
              & globalData%dftbEnergy(globalData%deltaDftb%iDeterminant)%EMermin +&
              & globalData%extPressure * globalData%cellVol
        end if
        call writeMdOut2(globalData%fdMd, globalData%tStress, globalData%tBarostat,&
            & globalData%tPeriodic, globalData%isLinResp, globalData%tEField, globalData%tFixEf,&
            & globalData%tPrintMulliken, globalData%dftbEnergy(globalData%deltaDftb%iDeterminant),&
            & globalData%energiesCasida, globalData%latVec, globalData%cellVol,&
            & globalData%intPressure, globalData%extPressure, tempIon, globalData%absEField,&
            & globalData%qOutput, globalData%q0, globalData%dipoleMoment)
        call writeCurrentGeometry(globalData%geoOutFile, globalData%pCoord0Out, .false., .true.,&
            & .true., globalData%tFracCoord, globalData%tPeriodic, globalData%tHelical,&
            & globalData%tPrintMulliken, globalData%species0, globalData%speciesName,&
            & globalData%latVec, globalData%origin, iGeoStep, iLatGeoStep, globalData%nSpin,&
            & globalData%qOutput, globalData%velocities)
      end if
      globalData%coord0(:,:) = globalData%newCoords
      if (globalData%tWriteDetailedOut  .and. globalData%deltaDftb%nDeterminant() == 1) then
        call writeDetailedOut5(globalData%fdDetailedOut, globalData%tPrintForces,&
            & globalData%tSetFillingTemp, globalData%tPeriodic, globalData%tStress,&
            & globalData%totalStress, globalData%totalLatDeriv,&
            & globalData%dftbEnergy(globalData%deltaDftb%iDeterminant), globalData%tempElec,&
            & globalData%extPressure, globalData%intPressure, tempIon)
      end if
    else if (globalData%tSocket .and. iGeoStep < globalData%nGeoSteps) then
      ! Only receive geometry from socket, if there are still geometry iterations left
    #:if WITH_SOCKETS
      call receiveGeometryFromSocket(env, globalData%socket, globalData%tPeriodic,&
          & globalData%coord0, globalData%latVec, globalData%tCoordsChanged, tLatticeChanged,&
          & tStopDriver)
    #:else
      call error("Internal error: code compiled without socket support")
    #:endif
    end if

  end subroutine getNextGeometry



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
  subroutine handleLatticeChange(latVecs, sccCalc, tStress, extPressure, mCutOff, dispersion,&
      & solvation, cm5Cont, recVecs, recVecs2p, cellVol, recCellVol, extLatDerivs, cellVecs,&
      & rCellVecs)

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
  subroutine handleCoordinateChange(env, coord0, latVec, invLatVec, species0, cutOff,&
      & orb, tPeriodic, tHelical, sccCalc, dispersion, solvation, thirdOrd, rangeSep, reks,&
      & img2CentCell, iCellVec, neighbourList, nAllAtom, coord0Fold, coord, species, rCellVec,&
      & nNeighbourSK, nNeighbourRep, nNeighbourLC, ham, over, H0, rhoPrim, iRhoPrim, iHam,&
      & ERhoPrim, iSparseStart, tPoisson, cm5Cont, stat)

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
    type(TCutoffs), intent(in) :: cutOff

    !> Atomic orbital information
    type(TOrbitals), intent(in) :: orb

    !> Is the geometry periodic
    logical, intent(in) :: tPeriodic

    !> Is the geometry helical
    logical, intent(in) :: tHelical

    !> SCC module internal variables
    type(TScc), allocatable, intent(inout) :: sccCalc

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

    !> Number of neighbours of each real atom close enough for repulsive interactions
    integer, intent(out) :: nNeighbourRep(:)

    !> Number of neighbours for each of the atoms for the exchange contributions in the long range
    !> functional
    integer, intent(inout), allocatable :: nNeighbourLC(:)

    !> Sparse hamiltonian storage
    real(dp), allocatable, intent(inout) :: ham(:,:)

    !> sparse overlap storage
    real(dp), allocatable, intent(inout) :: over(:)

    !> Non-SCC hamiltonian storage
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

    !> Charge model 5
    type(TChargeModel5), allocatable, intent(inout) :: cm5Cont

    !> Status of operation
    integer, intent(out), optional :: stat

    !> Total size of orbitals in the sparse data structures, where the decay of the overlap sets the
    !> sparsity pattern
    integer :: sparseSize

    coord0Fold(:,:) = coord0
    if (tPeriodic .or. tHelical) then
      call foldCoordToUnitCell(coord0Fold, latVec, invLatVec)
    end if

    if (tHelical) then
      call updateNeighbourListAndSpecies(coord, species, img2CentCell, iCellVec, neighbourList,&
          & nAllAtom, coord0Fold, species0, cutoff%mCutoff, rCellVec, helicalBoundConds=latVec)
    else
      call updateNeighbourListAndSpecies(coord, species, img2CentCell, iCellVec, neighbourList,&
          & nAllAtom, coord0Fold, species0, cutoff%mCutOff, rCellVec)
    end if

    call getNrOfNeighboursForAll(nNeighbourSK, neighbourList, cutoff%skCutOff)

    call getSparseDescriptor(neighbourList%iNeighbour, nNeighbourSK, img2CentCell, orb,&
        & iSparseStart, sparseSize)
    call reallocateSparseArrays(sparseSize, reks, ham, over, H0,&
        & rhoPrim, iHam, iRhoPrim, ERhoPrim)

    ! count neighbours for repulsive interactions between atoms
    call getNrOfNeighboursForAll(nNeighbourRep, neighbourList, cutoff%repCutOff)

    if (allocated(nNeighbourLC)) then
      ! count neighbours for repulsive interactions between atoms
      call getNrOfNeighboursForAll(nNeighbourLC, neighbourList, cutoff%lcCutOff)
    end if

    ! Notify various modules about coordinate changes
    if (tPoisson) then
      !! TODO: poiss_updcoords pass coord0 and not coord0Fold because the
      !! folding can mess up the contact position. Could we have the supercell
      !! centered on the input atomic structure?
      call poiss_updcoords(coord0)
    end if

    if (allocated(sccCalc)) then
      call sccCalc%updateCoords(env, coord, species, neighbourList)
    end if

    if (allocated(dispersion)) then
      call dispersion%updateCoords(env, neighbourList, img2CentCell, coord, species0, stat)
      @:HANDLE_ERROR(stat)
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
  subroutine reallocateSparseArrays(sparseSize, reks, ham, over,&
      & H0, rhoPrim, iHam, iRhoPrim, ERhoPrim)

    !> Size of the sparse overlap
    integer, intent(in) :: sparseSize

    !> data type for REKS
    type(TReksCalc), allocatable, intent(inout) :: reks

    !> Sparse storage for hamiltonian (sparseSize,nSpin)
    real(dp), allocatable, intent(inout) :: ham(:,:)

    !> Sparse storage for overlap
    real(dp), allocatable, intent(inout) :: over(:)

    !> Sparse storage for non-SCC hamiltonian
    real(dp), allocatable, intent(inout) :: H0(:)

    !> Sparse storage for density matrix
    real(dp), allocatable, intent(inout) :: rhoPrim(:,:)

    !> Sparse storage for imaginary hamiltonian (not reallocated if not initially allocated)
    real(dp), allocatable, intent(inout) :: iHam(:,:)

    !> Sparse storage for imaginary part of density matrix (not reallocated if not initially
    !> allocated)
    real(dp), allocatable, intent(inout) :: iRhoPrim(:,:)

    !> Sparse storage for energy weighted density matrix (not reallocated if not initially
    !> allocated)
    real(dp), allocatable, intent(inout) :: ERhoPrim(:)

    integer :: nSpin

    #:block DEBUG_CODE
      @:ASSERT(size(H0) == size(over))
      if (.not. allocated(reks)) then
        @:ASSERT(size(ham, dim=1) == size(over))
        @:ASSERT(all(shape(rhoPrim) == shape(ham)))
        if (allocated(iRhoPrim)) then
          @:ASSERT(all(shape(iRhoPrim) == shape(rhoPrim)))
          @:ASSERT(all(shape(iHam) == shape(ham)))
        end if
        if (allocated(ERhoPrim)) then
          @:ASSERT(size(ERhoPrim) == size(rhoPrim, dim=1))
        end if
      end if
    #:endblock DEBUG_CODE

    if (allocated(reks)) then
      if (size(over, dim=1) == sparseSize) then
        ! When the size of sparse matrices are different,
        ! phase of MOs can affect gradient of REKS
        return
      end if
    else
      if (size(ham, dim=1) >= sparseSize) then
        ! Sparse matrices are big enough
        return
      end if
    end if

    nSpin = size(rhoPrim, dim=2)
    if (.not. allocated(reks)) then
      deallocate(ham)
    end if
    deallocate(over)
    deallocate(H0)
    deallocate(rhoPrim)
    if (.not. allocated(reks)) then
      allocate(ham(sparseSize, nSpin))
    end if
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
    if (allocated(reks)) then
      call reks%reallocate(sparseSize)
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


#:if WITH_TRANSPORT

  !> Replace charges with those from the stored contact values
  subroutine overrideContactCharges(qInput, chargeUp, transpar, qBlockInput, blockUp)

    !> input charges
    real(dp), intent(inout) :: qInput(:,:,:)

    !> uploaded charges
    real(dp), intent(in) :: chargeUp(:,:,:)

    !> Transport parameters
    type(TTransPar), intent(in) :: transpar

    !> block charges, for example from DFTB+U
    real(dp), allocatable, intent(inout) :: qBlockInput(:,:,:,:)

    !> uploaded block charges
    real(dp), allocatable, intent(in) :: blockUp(:,:,:,:)

    integer :: ii, iStart, iEnd

    do ii = 1, transpar%ncont
      iStart = transpar%contacts(ii)%idxrange(1)
      iEnd = transpar%contacts(ii)%idxrange(2)
      qInput(:,iStart:iEnd,:) = chargeUp(:,iStart:iEnd,:)
    end do

  @:ASSERT(allocated(qBlockInput) .eqv. allocated(blockUp))
    if (allocated(qBlockInput)) then
      do ii = 1, transpar%ncont
        iStart = transpar%contacts(ii)%idxrange(1)
        iEnd = transpar%contacts(ii)%idxrange(2)
        qBlockInput(:,:,iStart:iEnd,:) = blockUp(:,:,iStart:iEnd,:)
      end do
    end if

  end subroutine overrideContactCharges

#:endif


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
      & rhoPrim, iHam, xi, orbitalL, HSqrReal, SSqrReal, eigvecsReal, iRhoPrim, HSqrCplx, SSqrCplx,&
      & eigvecsCplx, rhoSqrReal, deltaRhoInSqr, deltaRhoOutSqr, qOutput, nNeighbourLC,&
      & tLargeDenseMatrices, deltaDftb)

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

    !> Is the hamiltonian real (no k-points/molecule/gamma point)?
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
    type(TRangeSepFunc), allocatable, intent(inout) :: rangeSep

    !> eigenvalues (level, kpoint, spin)
    real(dp), intent(out) :: eigen(:,:,:)

    !> occupations (level, kpoint, spin)
    real(dp), intent(out) :: filling(:,:,:)

    !> sparse density matrix
    real(dp), intent(out) :: rhoPrim(:,:)

    !> imaginary part of hamiltonian
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

    !> Determinant derived type
    type(TDftbDeterminants), intent(inout) :: deltaDftb

    integer :: nSpin, iKS, iSp, iK, nAtom
    complex(dp), allocatable :: rhoSqrCplx(:,:)
    logical :: tImHam

    nSpin = size(ham, dim=2)
    tImHam = allocated(iRhoPrim)

    select case (electronicSolver%iSolver)

    case (electronicSolverTypes%GF)

      call env%globalTimer%startTimer(globalTimers%densityMatrix)
    #:if WITH_TRANSPORT
      call calcdensity_green(iSCC, env, parallelKS%localKS, ham, over, neighbourlist%iNeighbour,&
          & nNeighbourSK, denseDesc%iAtomStart, iSparseStart, img2CentCell, iCellVec, cellVec, orb,&
          & kPoint, kWeight, mu, rhoPrim, energy%Eband, Ef, energy%E0, energy%TS)
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
          & tempElec, nEl, parallelKS, Ef, energy, rangeSep, eigen, filling, rhoPrim, iHam, xi,&
          & orbitalL, HSqrReal, SSqrReal, eigvecsReal, iRhoPrim, HSqrCplx, SSqrCplx, eigvecsCplx,&
          & rhoSqrReal, deltaRhoInSqr, deltaRhoOutSqr, qOutput, nNeighbourLC, deltaDftb)

    case(electronicSolverTypes%omm, electronicSolverTypes%pexsi, electronicSolverTypes%ntpoly,&
        &electronicSolverTypes%elpadm)

      call env%globalTimer%startTimer(globalTimers%densityMatrix)

      call electronicSolver%elsi%getDensity(env, denseDesc, ham, over, neighbourList, nNeighbourSK,&
          & iSparseStart, img2CentCell, iCellVec, cellVec, kPoint, kWeight, tHelical, orb, species,&
          & coord, tRealHS, tSpinSharedEf, tSpinOrbit, tDualSpinOrbit, tMulliken, parallelKS, Ef,&
          & energy, rhoPrim, energy%Eband, energy%TS, iHam, xi, orbitalL, HSqrReal, SSqrReal,&
          & iRhoPrim, HSqrCplx, SSqrCplx)
      call env%globalTimer%stopTimer(globalTimers%densityMatrix)

    end select

  end subroutine getDensity


  !> Returns the density matrix using dense diagonalisation.
  subroutine getDensityFromDenseDiag(env, denseDesc, ham, over, neighbourList, nNeighbourSK,&
      & iSparseStart, img2CentCell, iCellVec, cellVec, kPoint, kWeight, orb, iAtomStart, tHelical,&
      & coord, species, electronicSolver, tRealHS, tSpinSharedEf, tSpinOrbit, tDualSpinOrbit,&
      & tFillKSep, tFixEf, tMulliken, iDistribFn, tempElec, nEl, parallelKS, Ef, energy, rangeSep,&
      & eigen, filling, rhoPrim, iHam, xi, orbitalL, HSqrReal, SSqrReal, eigvecsReal, iRhoPrim,&
      & HSqrCplx, SSqrCplx, eigvecsCplx, rhoSqrReal, deltaRhoInSqr, deltaRhoOutSqr, qOutput,&
      & nNeighbourLC, deltaDftb)

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
    type(TRangeSepFunc), allocatable, intent(inout) :: rangeSep

    !> eigenvalues (level, kpoint, spin)
    real(dp), intent(out) :: eigen(:,:,:)

    !> occupations (level, kpoint, spin)
    real(dp), intent(out) :: filling(:,:,:)

    !> sparse density matrix
    real(dp), intent(out) :: rhoPrim(:,:)

    !> imaginary part of hamiltonian
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

    !> Determinant derived type
    type(TDftbDeterminants), intent(inout) :: deltaDftb

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

    !>Data for rangeseparated calcualtion
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

    integer :: iKS, iSpin

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
        call unpackHelicalHS(HSqrReal, ham(:,iSpin), neighbourList%iNeighbour, nNeighbourSK,&
            & denseDesc%iAtomStart, iSparseStart, img2CentCell, orb, species, coord)
        call unpackHelicalHS(SSqrReal, over, neighbourList%iNeighbour, nNeighbourSK,&
            & denseDesc%iAtomStart, iSparseStart, img2CentCell, orb, species, coord)
      else
        call unpackHS(HSqrReal, ham(:,iSpin), neighbourList%iNeighbour, nNeighbourSK,&
            & denseDesc%iAtomStart, iSparseStart, img2CentCell)
        call unpackHS(SSqrReal, over, neighbourList%iNeighbour, nNeighbourSK, denseDesc%iAtomStart,&
            & iSparseStart, img2CentCell)
      end if
      call env%globalTimer%stopTimer(globalTimers%sparseToDense)

      ! Add rangeseparated contribution
      ! Assumes deltaRhoInSqr only used by rangeseparation
      ! Should this be used elsewhere, need to pass isRangeSep
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

    !> dense hamiltonian matrix
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
        call unpackHelicalHS(HSqrCplx, ham(:,iSpin), kPoint(:,iK), neighbourList%iNeighbour,&
            & nNeighbourSK, iCellVec, cellVec, denseDesc%iAtomStart, iSparseStart, img2CentCell,&
            & orb, species, coord)
        call unpackHelicalHS(SSqrCplx, over, kPoint(:,iK), neighbourList%iNeighbour, nNeighbourSK,&
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

    !> dense hamiltonian matrix
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
        if (tDensON2) then
          call makeDensityMatrix(work, eigvecs(:,:,iKS), filling(:,iSpin),&
              & neighbourlist%iNeighbour, nNeighbourSK, orb, denseDesc%iAtomStart, img2CentCell)
        else
          call makeDensityMatrix(work, eigvecs(:,:,iKS), filling(:,iSpin))
        end if
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
      if (tDensON2) then
        call makeDensityMatrix(work, eigvecs(:,:,iKS), filling(:,iK,iSpin),&
            & neighbourlist%iNeighbour, nNeighbourSK, orb, denseDesc%iAtomStart, img2CentCell)
      else
        call makeDensityMatrix(work, eigvecs(:,:,iKS), filling(:,iK,iSpin))
      end if
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
        call deltaDftb%detFilling(fillings(:,iK:iK,:), EBand, EfTmp, TSTmp, E0Tmp, nElecFill,&
            & eigVals(:,iK:iK,:), tempElec, kWeightTmp(:nSpinHams), iDistribFn)
        Eband(:) = Eband + EbandTmp(1) * kWeights(iK)
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
  subroutine getMullikenPopulation(rhoPrim, over, orb, neighbourList, nNeighbourSK, img2CentCell,&
      & iSparseStart, qOrb, iRhoPrim, qBlock, qiBlock, qNetAtom)

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

    !> Onsite Mulliken charges per atom
    real(dp), intent(inout), optional :: qNetAtom(:)

    integer :: iSpin

    qOrb(:,:,:) = 0.0_dp
    do iSpin = 1, size(rhoPrim, dim=2)
      call mulliken(qOrb(:,:,iSpin), over, rhoPrim(:,iSpin), orb, neighbourList%iNeighbour,&
          & nNeighbourSK, img2CentCell, iSparseStart)
    end do

    if (allocated(qBlock)) then
      qBlock(:,:,:,:) = 0.0_dp
      do iSpin = 1, size(rhoPrim, dim=2)
        call mulliken(qBlock(:,:,:,iSpin), over, rhoPrim(:,iSpin), orb, neighbourList%iNeighbour,&
            & nNeighbourSK, img2CentCell, iSparseStart)
      end do
    end if

    if (allocated(qiBlock)) then
      qiBlock(:,:,:,:) = 0.0_dp
      do iSpin = 1, size(iRhoPrim, dim=2)
        call skewMulliken(qiBlock(:,:,:,iSpin), over, iRhoPrim(:,iSpin), orb,&
            & neighbourList%iNeighbour, nNeighbourSK, img2CentCell, iSparseStart)
      end do
    end if

    if (present(qNetAtom)) then
      call getOnsitePopulation(rhoPrim(:,1), orb, iSparseStart, qNetAtom)
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
      & qInput, qInpRed, sccErrorQ, tConverged, qBlockOut, iEqBlockDftbU, qBlockIn, qiBlockOut,&
      & iEqBlockDftbuLS, species0, nUJ, iUJ, niUJ, qiBlockIn, iEqBlockOnSite, iEqBlockOnSiteLS)

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
    real(dp), intent(inout), allocatable :: qBlockIn(:,:,:,:)

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
          call unpackHelicalHS(SSqrReal, over, neighbourList%iNeighbour, nNeighbourSK, iAtomStart,&
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
      & over, eigvecsReal, eigen, filling, coord, species, speciesName, orb, skHamCont,&
      & skOverCont, autotestTag, taggedWriter, runId, neighbourList, nNeighbourSk, denseDesc,&
      & iSparseStart, img2CentCell, tWriteAutotest, tForces, tLinRespZVect, tPrintExcEigvecs,&
      & tPrintExcEigvecsTxt, nonSccDeriv, dftbEnergy, energies, work, rhoSqrReal, excitedDerivs,&
      & dQAtomEx, occNatural)

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

    !> energes of all solved states
    real(dp), intent(inout), allocatable :: energies(:)

    !> Working array of the size of the dense matrices.
    real(dp), intent(out) :: work(:,:)

    !> density matrix in dense form
    real(dp), intent(inout), allocatable :: rhoSqrReal(:,:,:)

    !> excited state energy derivative with respect to atomic coordinates
    real(dp), intent(inout), allocatable :: excitedDerivs(:,:)

    !> Mulliken charges in excited state
    real(dp), intent(out) :: dQAtomEx(:)

    !> natural orbital occupation numbers
    real(dp), intent(inout), allocatable :: occNatural(:)

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
    call unpackHS(work, over, neighbourList%iNeighbour, nNeighbourSK, denseDesc%iAtomStart,&
        & iSparseStart, img2CentCell)
    call blockSymmetrizeHS(work, denseDesc%iAtomStart)
    if (allocated(rhoSqrReal)) then
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
      call addGradients(tSpin, linearResponse, denseDesc%iAtomStart, eigvecsReal, eigen, work,&
          & filling, coord(:,:nAtom), sccCalc, dQAtom, pSpecies0, neighbourList%iNeighbour,&
          & img2CentCell, orb, skHamCont, skOverCont, tWriteAutotest, fdAutotest, taggedWriter,&
          & dftbEnergy%Eexcited, energies, excitedDerivs, nonSccDeriv,&
          & rhoSqrReal, occNatural, naturalOrbs)
      if (tPrintExcEigvecs) then
        call writeRealEigvecs(env, runId, neighbourList, nNeighbourSK, denseDesc, iSparseStart,&
            & img2CentCell, pSpecies0, speciesName, orb, over, parallelKS, tPrintExcEigvecsTxt,&
            & naturalOrbs, work, fileName="excitedOrbs")
      end if
    else
      call linResp_calcExcitations(linearResponse, tSpin, denseDesc, eigvecsReal, eigen, work,&
          & filling, coord(:,:nAtom), sccCalc, dQAtom, pSpecies0, neighbourList%iNeighbour,&
          & img2CentCell, orb, tWriteAutotest, fdAutotest, taggedWriter,&
          & dftbEnergy%Eexcited, energies)
    end if
    dftbEnergy%Etotal = dftbEnergy%Etotal + dftbEnergy%Eexcited
    dftbEnergy%EMermin = dftbEnergy%EMermin + dftbEnergy%Eexcited
    dftbEnergy%EGibbs = dftbEnergy%EGibbs + dftbEnergy%Eexcited
    if (tWriteAutotest) then
      close(fdAutotest)
    end if

  end subroutine calculateLinRespExcitations


  !> Do the linear response excitation calculation with range-separated Hamiltonian.
  subroutine calculateLinRespExcitations_RS(env, linearResponse, parallelKS, sccCalc, qOutput, q0,&
      & over, eigvecsReal, eigen, filling, coord0, species, speciesName, orb, skHamCont,&
      & skOverCont, autotestTag, taggedWriter, runId, neighbourList, nNeighbourSk, denseDesc,&
      & iSparseStart, img2CentCell, tWriteAutotest, tForces, tLinRespZVect, tPrintExcEigvecs,&
      & tPrintExcEigvecsTxt, nonSccDeriv, dftbEnergy, energies, work, deltaRhoOutSqr,&
      & excitedDerivs, dQAtomEx, occNatural, rangeSep)

    !> Environment settings
    type(TEnvironment), intent(in) :: env

    !> excited state settings
    type(TLinResp), intent(inout) :: linearResponse

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
    real(dp), intent(inout) :: eigvecsReal(:,:,:)

    !> ground state eigenvalues (orbital, kpoint)
    real(dp), intent(in) :: eigen(:,:)

    !> ground state fillings (orbital, kpoint)
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

    !> energes of all solved states
    real(dp), intent(inout), allocatable :: energies(:)

    !> Working array of the size of the dense matrices.
    real(dp), intent(out) :: work(:,:)

    !> density matrix in dense form
   !real(dp), intent(inout), allocatable :: rhoSqrReal(:,:,:)
    real(dp), intent(inout) :: deltaRhoOutSqr(:,:,:)

    !> excited state energy derivative with respect to atomic coordinates
    real(dp), intent(inout), allocatable :: excitedDerivs(:,:)

    !> Mulliken charges in excited state
    real(dp), intent(out) :: dQAtomEx(:)

    !> natural orbital occupation numbers
    real(dp), intent(inout), allocatable :: occNatural(:)

    !> Data from rangeseparated calculations
    type(TRangeSepFunc), intent(inout), allocatable :: rangeSep


    real(dp), allocatable :: dQAtom(:)
    real(dp), allocatable :: naturalOrbs(:,:,:)
    integer, pointer :: pSpecies0(:)
    integer :: iSpin, nSpin, nAtom, fdAutotest
    logical :: tSpin
    ! Onsite corrections -- remain dummy as they are not calculated
    logical, parameter :: tOnsite = .false.
    integer :: norb

    nAtom = size(qOutput, dim=2)
    nSpin = size(eigen, dim=2)
    tSpin = (nSpin == 2)
    pSpecies0 => species(1:nAtom)

    norb = size(eigvecsReal, dim=1)

    dftbEnergy%Eexcited = 0.0_dp
    allocate(dQAtom(nAtom))
    dQAtom(:) = sum(qOutput(:,:,1) - q0(:,:,1), dim=1)
    call unpackHS(work, over, neighbourList%iNeighbour, nNeighbourSK, denseDesc%iAtomStart,&
        & iSparseStart, img2CentCell)
    call blockSymmetrizeHS(work, denseDesc%iAtomStart)
    !call unpackHS(SSqrReal, over, neighbourList%iNeighbour, nNeighbour,&
    !     & denseDesc%iAtomStart, iPair, img2CentCell)
    !call blockSymmetrizeHS(SSqrReal, iAtomStart)
    if (tForces) then
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
      ! WITH FORCES
      excitedDerivs = 0.0_dp
    #:if WITH_ARPACK
      call linRespCalcExcitationsRS(tSpin, tOnsite, linearResponse, denseDesc%iAtomStart,&
          & eigvecsReal, eigen, sccCalc, work, filling, coord0, dQAtom, pSpecies0,&
          & linearResponse%HubbardU, neighbourList%iNeighbour, img2CentCell, orb, rangeSep,&
          & tWriteAutotest, fdAutotest, taggedWriter, dftbEnergy%Eexcited, skHamCont,&
          & skOverCont, nonSccDeriv, deltaRhoOutSqr(:,:,1), excitedDerivs, dQAtomEx)
    #:else
      call error("Should not be here - compiled without ARPACK")
    #:endif
      if (tPrintExcEigvecs) then
        call writeRealEigvecs(env, runId, neighbourList, nNeighbourSK, denseDesc, iSparseStart,&
            & img2CentCell, pSpecies0, speciesName, orb, over, parallelKS, tPrintExcEigvecsTxt,&
            & naturalOrbs, work, fileName="excitedOrbs")
      end if
    else
      ! NO FORCES
    #:if WITH_ARPACK
      call linRespCalcExcitationsRS(tSpin, tOnsite, linearResponse, denseDesc%iAtomStart,&
          & eigvecsReal, eigen, sccCalc, work, filling, coord0, dQAtom, pSpecies0,&
          & linearResponse%HubbardU, neighbourList%iNeighbour, img2CentCell, orb, rangeSep,&
          & tWriteAutotest, fdAutotest, taggedWriter, dftbEnergy%Eexcited)
    #:else
      call error("Should not be here - compiled without ARPACK")
    #:endif
    end if

    dftbEnergy%Etotal = dftbEnergy%Etotal + dftbEnergy%Eexcited
    dftbEnergy%EMermin = dftbEnergy%EMermin + dftbEnergy%Eexcited
    dftbEnergy%EGibbs = dftbEnergy%EGibbs + dftbEnergy%Eexcited
    if (tWriteAutotest) then
      close(fdAutotest)
    end if

  end subroutine calculateLinRespExcitations_RS

  !> Get the XLBOMD charges for the current geometry.
  subroutine getXlbomdCharges(xlbomdIntegrator, qOutRed, pChrgMixer, orb, nIneqOrb, iEqOrbitals,&
      & qInput, qInpRed, iEqBlockDftbu, qBlockIn, species0, nUJ, iUJ, niUJ, iEqBlockDftbuLS,&
      & qiBlockIn, iEqBlockOnSite, iEqBlockOnSiteLS)

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


  !> Prints dipole moment calculated by the derivative of H with respect to the external field.
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

    !> Is the hamiltonian real (no k-points/molecule/gamma point)?
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

    !> Storage for dense hamiltonian matrix (complex case)
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

    case (electronicSolverTypes%omm, electronicSolverTypes%pexsi, electronicSolverTypes%ntpoly,&
        &electronicSolverTypes%elpadm)

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

    !> Is the hamiltonian real (no k-points/molecule/gamma point)?
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

    !> Storage for dense hamiltonian matrix (complex case)
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
          call unpackHelicalHS(work, ham(:,iS), neighbourlist%iNeighbour, nNeighbourSK,&
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
          call unpackHelicalHS(work2, ham(:,iS), neighbourlist%iNeighbour, nNeighbourSK,&
              & denseDesc%iAtomStart, iSparseStart, img2CentCell, orb, species, coord)
        else
          call unpackHS(work2, ham(:,iS), neighbourlist%iNeighbour, nNeighbourSK,&
              & denseDesc%iAtomStart, iSparseStart, img2CentCell)
        end if
        call blocksymmetrizeHS(work2, denseDesc%iAtomStart)
        call symm(eigvecsReal(:,:,iKS), "L", work, work2)
        if (tHelical) then
          call unpackHelicalHS(work, over, neighbourlist%iNeighbour, nNeighbourSK,&
              & denseDesc%iAtomStart, iSparseStart, img2CentCell, orb, species, coord)
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
          call unpackHelicalHS(work, ham(:,iS), kPoint(:,iK), neighbourlist%iNeighbour,&
              & nNeighbourSK, iCellVec, cellVec, denseDesc%iAtomStart, iSparseStart, img2CentCell,&
              & orb, species, coord)
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
          call unpackHelicalHS(work2, ham(:,iS), kPoint(:,iK), neighbourlist%iNeighbour,&
              & nNeighbourSK, iCellVec, cellVec, denseDesc%iAtomStart, iSparseStart, img2CentCell,&
              & orb, species, coord)
        else
          call unpackHS(work2, ham(:,iS), kPoint(:,iK), neighbourlist%iNeighbour, nNeighbourSK,&
              & iCellVec, cellVec, denseDesc%iAtomStart, iSparseStart, img2CentCell)
        end if
        call blockHermitianHS(work2, denseDesc%iAtomStart)
        call hemm(eigvecsCplx(:,:,iKS), "L", work, work2)
        if  (tHelical) then
          call unpackHelicalHS(work, over, kPoint(:,iK), neighbourlist%iNeighbour, nNeighbourSK,&
              & iCellVec, cellVec, denseDesc%iAtomStart, iSparseStart, img2CentCell, orb, species,&
              & coord)
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
  subroutine getGradients(env, sccCalc, tExtField, isXlbomd, nonSccDeriv, EField, rhoPrim,&
      & ERhoPrim, qOutput, q0, skHamCont, skOverCont, pRepCont, neighbourList, nNeighbourSK,&
      & nNeighbourRep, species, img2CentCell, iSparseStart, orb, potential, coord, derivs,&
      & groundDerivs, tripletderivs, mixedderivs, iRhoPrim, thirdOrd, solvation, qDepExtPot,&
      & chrgForces, dispersion, rangeSep, SSqrReal, over, denseDesc, deltaRhoOutSqr, tPoisson,&
      & halogenXCorrection, tHelical, coord0, deltaDftb)

    !> Environment settings
    type(TEnvironment), intent(inout) :: env

    !> SCC module internal variables
    type(TScc), allocatable, intent(in) :: sccCalc

    !> external electric field
    logical, intent(in) :: tExtField

    !> extended Lagrangian active?
    logical, intent(in) :: isXlbomd

    !> method for calculating derivatives of S and H0
    type(TNonSccDiff), intent(in) :: nonSccDeriv

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
    type(TSlakoCont), intent(in) :: skHamCont

    !> overlap information
    type(TSlakoCont), intent(in) :: skOverCont

    !> repulsive information
    type(TRepCont), intent(in) :: pRepCont

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

    !> Determinant derived type
    type(TDftbDeterminants), intent(in) :: deltaDftb

    ! Locals
    real(dp), allocatable :: tmpDerivs(:,:)
    real(dp), allocatable :: dummyArray(:,:)
    real(dp), allocatable :: dQ(:,:,:)
    logical :: tImHam, tExtChrg, tSccCalc
    integer :: nAtom, iAt


    tSccCalc = allocated(sccCalc)
    tImHam = allocated(iRhoPrim)
    tExtChrg = allocated(chrgForces)
    nAtom = size(derivs, dim=2)

    allocate(tmpDerivs(3, nAtom))
    if (tPoisson) then
      allocate(dummyArray(orb%mshell, nAtom))
    end if
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
        call poiss_getshift(env, dummyArray, tmpDerivs)
        derivs(:,:) = derivs + tmpDerivs

      else

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
        endif

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

        if (tExtField) then
          do iAt = 1, nAtom
            derivs(:, iAt) = derivs(:, iAt)&
                & + sum(qOutput(:, iAt, 1) - q0(:, iAt, 1)) * potential%extGrad(:, iAt)
          end do
        end if

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
      call dispersion%addGradients(derivs)
    end if

    if (allocated(halogenXCorrection)) then
      call halogenXCorrection%addGradients(derivs, coord, species, neighbourList, img2CentCell)
    end if

    if (allocated(rangeSep)) then
      if (tHelical) then
        call unpackHelicalHS(SSqrReal, over, neighbourList%iNeighbour, nNeighbourSK,&
            & denseDesc%iAtomStart, iSparseStart, img2CentCell, orb, species, coord)
      else
        call unpackHS(SSqrReal, over, neighbourList%iNeighbour, nNeighbourSK, denseDesc%iAtomStart,&
            & iSparseStart, img2CentCell)
      end if
      if (size(deltaRhoOutSqr, dim=3) > 2) then
        call error("Range separated forces do not support non-colinear spin")
      else
        call rangeSep%addLRGradients(derivs, nonSccDeriv, deltaRhoOutSqr, skOverCont, coord,&
            & species, orb, denseDesc%iAtomStart, SSqrReal, neighbourList%iNeighbour, nNeighbourSK)
      end if
    end if


    call getERepDeriv(tmpDerivs, coord, nNeighbourRep, neighbourList%iNeighbour, species, pRepCont,&
        & img2CentCell, tHelical)

    derivs(:,:) = derivs + tmpDerivs

    call helicalTwistFolded(derivs, coord, coord0, nAtom, tHelical)

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


  !> Correct for z folding into central unit cell requiring a twist in helical cases
  pure subroutine helicalTwistFolded(derivs, coord, coord0, nAtom, tHelical)

    use dftbp_quaternions, only : rotate3
    use dftbp_boundarycond, only : zAxis

    !> Derivatives
    real(dp), intent(inout) :: derivs(:,:)

    !> Unfolded atoms
    real(dp), intent(in) :: coord(:,:)

    !> Central cell atoms
    real(dp), intent(in) :: coord0(:,:)

    !> number of atoms
    integer, intent(in) :: nAtom

    !> Is this a helical geometry
    logical, intent(in) :: tHelical

    integer :: ii
    real(dp) :: deltaTheta

    if (tHelical) then
      do ii = 1, nAtom
        deltaTheta = atan2(coord0(2,ii),coord0(1,ii)) - atan2(coord(2,ii),coord(1,ii))
        call rotate3(derivs(:,ii), deltaTheta, zAxis)
      end do
    end if

  end subroutine helicalTwistFolded


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
      & totalStress, totalLatDeriv, intPressure, iRhoPrim, solvation, dispersion,&
      & halogenXCorrection, deltaDftb, tripletStress, mixedStress)

    !> Environment settings
    type(TEnvironment), intent(in) :: env

    !> SCC module internal variables
    type(TScc), allocatable, intent(in) :: sccCalc

    !> Is 3rd order SCC being used
    type(TThirdOrder), intent(inout), allocatable :: thirdOrd

    !> External field
    logical, intent(in) :: tExtField

    !> method for calculating derivatives of S and H0
    type(TNonSccDiff), intent(in) :: nonSccDeriv

    !> density matrix
    real(dp), intent(in) :: rhoPrim(:,:)

    !> energy weighted density matrix
    real(dp), intent(in) :: ERhoPrim(:)

    !> electrons in orbitals
    real(dp), intent(in) :: qOutput(:,:,:)

    !> refernce charges
    real(dp), intent(in) :: q0(:,:,:)

    !> non-SCC hamiltonian information
    type(TSlakoCont), intent(in) :: skHamCont

    !> overlap information
    type(TSlakoCont), intent(in) :: skOverCont

    !> repulsive information
    type(TRepCont), intent(in) :: pRepCont

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

    if (tExtField) then
      call getExtFieldStress(latVec, cellVol, q0, qOutput, potential%extGrad, coord0, tmpStress)
      totalStress(:,:) = totalStress + tmpStress
    end if

    call getRepulsiveStress(tmpStress, coord, nNeighbourRep, neighbourList%iNeighbour, species,&
        & img2CentCell, pRepCont, cellVol)
    totalStress(:,:) = totalStress + tmpStress

    if(deltaDftb%isNonAufbau) then
      select case (deltaDftb%whichDeterminant(deltaDftb%iDeterminant))
      case (determinants%triplet)
        tripletStress(:,:) = totalStress
      case (determinants%mixed)
        mixedStress(:,:) = totalStress
      end select
    end if

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
    type(TNumDerivs), intent(inout) :: derivDriver

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
        call unpackHelicalHS(SSqrReal,over,neighbourList%iNeighbour, nNeighbourSK,&
            & denseDesc%iAtomStart, iSparseStart, img2CentCell, orb, species, coord)
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
  subroutine getReksInitialSettings(env, denseDesc, h0, over, neighbourList, &
      & nNeighbourSK, iSparseStart, img2CentCell, electronicSolver, iGeoStep, &
      & HSqrReal, SSqrReal, eigvecsReal, eigen, reks)

    !> Environment settings
    type(TEnvironment), intent(inout) :: env

    !> Dense matrix descriptor
    type(TDenseDescr), intent(in) :: denseDesc

    !> hamiltonian in sparse storage
    real(dp), intent(in) :: h0(:)

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

    call env%globalTimer%startTimer(globalTimers%sparseToDense)
    call unpackHS(SSqrReal, over, neighbourList%iNeighbour, nNeighbourSK, &
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
        call diagDenseMtx(electronicSolver, 'V', HSqrReal, SSqrReal, eigen(:,1,1))
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
        & nNeighbourSK, iSparseStart, img2CentCell, over, reks)

  end subroutine getReksInitialSettings


  !> Normalize eigenvectors with unitary transformation
  subroutine renormalizeEigenvecs(env, electronicSolver, eigvecsReal, reks)

    use dftbp_blasroutines, only : gemm
    use dftbp_eigensolver, only : heev

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
    real(dp), intent(inout), allocatable  :: rhoSqrReal(:,:,:)

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
      & img2CentCell, iSparseStart, orb, rhoPrim, over, iRhoPrim, qBlock, &
      & qiBlock, reks)

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

    !> sparse overlap matrix
    real(dp), intent(in) :: over(:)

    !> imaginary part of density matrix
    real(dp), intent(in), allocatable :: iRhoPrim(:,:)

    !> Dual atomic charges
    real(dp), intent(inout), allocatable :: qBlock(:,:,:,:)

    !> Imaginary part of dual atomic charges
    real(dp), intent(inout), allocatable :: qiBlock(:,:,:,:)

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

      ! reks%qOutputL has (my_qm) component
      reks%qOutputL(:,:,:,iL) = 0.0_dp
      call getMullikenPopulation(rhoPrim, over, orb, neighbourList, nNeighbourSK, &
          & img2CentCell, iSparseStart, reks%qOutputL(:,:,:,iL), iRhoPrim=iRhoPrim, &
          & qBlock=qBlock, qiBlock=qiBlock)

    end do

    ! reks%qOutputL has (qm) component
    call qmExpandL(reks%qOutputL, reks%Lpaired)

  end subroutine getMullikenPopulationL


  !> Build L, spin dependent Hamiltonian with various contributions
  !> and compute the energy of microstates
  subroutine getHamiltonianLandEnergyL(env, denseDesc, sccCalc, orb, species, neighbourList, &
      & nNeighbourSK, iSparseStart, img2CentCell, H0, over, spinW, cellVol, extPressure, &
      & energy, q0, iAtInCentralRegion, solvation, thirdOrd, potential, electrostatics, &
      & tPoisson, tUpload, shiftPerLUp, rangeSep, nNeighbourLC, tDualSpinOrbit, xi, tExtField, &
      & isXlbomd, tDftbU, TS, qDepExtPot, qBlock, qiBlock, nDftbUFunc, UJ, nUJ, iUJ, niUJ,&
      & tFixEf, Ef, rhoPrim, onSiteElements, iHam, dispersion, reks)

    !> Environment settings
    type(TEnvironment), intent(inout) :: env

    !> Dense matrix descriptor
    type(TDenseDescr), intent(in) :: denseDesc

    !> SCC module internal variables
    type(TScc), allocatable, intent(inout) :: sccCalc

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

    !> sparse overlap matrix
    real(dp), intent(in) :: over(:)

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

    !> potentials in the system
    type(TPotentials), intent(inout) :: potential

    !> electrostatic solver (poisson or gamma-functional)
    integer, intent(in) :: electrostatics

    !> whether Poisson is solved (used with tPoissonTwice)
    logical, intent(in) :: tPoisson

    !> whether contacts are uploaded
    logical, intent(in) :: tUpload

    !> uploded potential per shell per atom
    real(dp), allocatable, intent(in) :: shiftPerLUp(:,:)

    !> Data for rangeseparated calculation
    type(TRangeSepFunc), allocatable, intent(inout) :: rangeSep

    !> Nr. of neighbours for each atom in the long-range functional.
    integer, allocatable, intent(in) :: nNeighbourLC(:)

    !> Is dual spin orbit being used (block potentials)
    logical, intent(in) :: tDualSpinOrbit

    !> Spin orbit constants if required
    real(dp), allocatable, intent(in) :: xi(:,:)

    !> is an external electric field present
    logical, intent(in) :: tExtField

    !> Is the extended Lagrangian being used for MD
    logical, intent(in) :: isXlbomd

    !> Are there orbital potentials present
    logical, intent(in) :: tDftbU

    !> electron entropy contribution
    real(dp), intent(in) :: TS(:)

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

    !> Whether fixed Fermi level(s) should be used. (No charge conservation!)
    logical, intent(in) :: tFixEf

    !> If tFixEf is .true. contains reservoir chemical potential, otherwise the Fermi levels found
    !> from the given number of electrons
    real(dp), intent(inout) :: Ef(:)

    !> sparse density matrix
    real(dp), intent(inout) :: rhoPrim(:,:)

    !> Corrections terms for on-site elements
    real(dp), intent(in), allocatable :: onSiteElements(:,:,:,:)

    !> imaginary part of hamiltonian (if required, signalled by being allocated)
    real(dp), allocatable, intent(inout) :: iHam(:,:)

    !> dispersion interactions
    class(TDispersionIface), allocatable, intent(in) :: dispersion

    !> data type for REKS
    type(TReksCalc), allocatable, intent(inout) :: reks

    real(dp), allocatable :: tmpHamSp(:,:)
    real(dp), allocatable :: tmpEn(:)

    integer, pointer :: pSpecies0(:)
    integer :: sparseSize, nAtom, nSpin, iL, tmpL, rsL

    sparseSize = size(over,dim=1)
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
      call addChargePotentials(env, sccCalc, reks%qOutputL(:,:,:,iL), q0, &
          & reks%chargePerShellL(:,:,:,iL), orb, species, neighbourList, &
          & img2CentCell, spinW, solvation, thirdOrd, potential, electrostatics, &
          & tPoisson, tUpload, shiftPerLUp, dispersion)

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
      call getSccHamiltonian(H0, over, nNeighbourSK, neighbourList, species, orb,&
          & iSparseStart, img2CentCell, potential, allocated(reks), tmpHamSp, iHam)
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
        call rangeSep%addLRHamiltonian(env, reks%deltaRhoSqrL(:,:,1,iL), over, &
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
      call sccCalc%updateCharges(env, reks%qOutputL(:,:,:,iL), q0, orb, species)
      call sccCalc%updateShifts(env, orb, species, neighbourList%iNeighbour, img2CentCell)
      potential%intShell(:,:,:) = reks%intShellL(:,:,:,iL)
      if (allocated(thirdOrd)) then
        call thirdOrd%updateCharges(pSpecies0, neighbourList, &
            & reks%qOutputL(:,:,:,iL), q0, img2CentCell, orb)
      end if

      call calcEnergies(sccCalc, reks%qOutputL(:,:,:,iL), q0, reks%chargePerShellL(:,:,:,iL),&
          & species, tExtField, isXlbomd, tDftbU, tDualSpinOrbit, rhoPrim, H0, orb,&
          & neighbourList, nNeighbourSk, img2CentCell, iSparseStart, cellVol, extPressure, TS,&
          & potential, energy, thirdOrd, solvation, rangeSep, reks, qDepExtPot, qBlock, qiBlock,&
          & nDftbUFunc, UJ, nUJ, iUJ, niUJ, xi, iAtInCentralRegion, tFixEf, Ef, onSiteElements)
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

    !> minumum number of SCC iterations to perform
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

    !> minumum number of SCC iterations to perform
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


end module dftbp_main
