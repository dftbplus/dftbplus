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
  use lapackroutines, only : matinv
  use taggedoutput
  use scc
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
  use sparse2dense, only : unpackHS, packHS, packERho, iPackHS,&
      & blockSymmetrizeHS, blockHermitianHS
  use blasroutines, only : symm, hemm
  use hsdutils
  use charmanip
  use shift
  use spinorbit
  use angmomentum
  use elecconstraints
  use pmlocalisation
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
  implicit none
  private

  public :: runDftbPlus

  ! Name of the human readable files
  character(*), parameter :: autotestTag = "autotest.tag"
  character(*), parameter :: userOut = "detailed.out"
  character(*), parameter :: bandOut = "band.out"
  character(*), parameter :: mdOut = "md.out"
  character(*), parameter :: resultsTag = "results.tag"
  character(*), parameter :: hessianOut = "hessian.out"

  
contains

  subroutine runDftbPlus()
    use initprogram

    real(dp), allocatable :: orbResShift3rd(:,:)
    real(dp), allocatable :: shift3rd(:)
    
    
    integer                  :: nk, iSpin2, iK2

    complex(dp), allocatable :: HSqrCplx(:,:,:,:), SSqrCplx(:,:), HSqrCplx2(:,:)
    real(dp),    allocatable :: HSqrReal(:,:,:), SSqrReal(:,:), HSqrReal2(:,:)
    real(dp),    allocatable :: eigen(:,:,:), eigen2(:,:,:)
    real(dp), allocatable :: rhoPrim(:,:)
    real(dp), allocatable :: iRhoPrim(:,:)
    real(dp), allocatable :: ERhoPrim(:)
    real(dp), allocatable :: h0(:)

    ! variables for derivatives using the Hellmann-Feynman theorem:

    !> for derivatives of H wrt external perturbation
    real(dp), allocatable :: hprime(:,:)

    !> for derivatives of V
    real(dp), allocatable :: potentialDerivative(:,:)

    !> temporary dipole data
    real(dp), allocatable :: dipoleTmp(:,:)


    !> electronic filling
    real(dp), allocatable :: filling(:,:,:)

    !> band structure energy
    real(dp), allocatable :: Eband(:)

    !> entropy of electrons at temperature T
    real(dp), allocatable :: TS(:)

    !> zero temperature electronic energy
    real(dp), allocatable :: E0(:)

    !> energy in previous scc cycles
    real(dp), allocatable :: Eold


    !> Total energy components
    type(TEnergies) :: energy

    !> Potentials for orbitals
    type(TPotentials) :: potential

    real(dp), allocatable :: derivs(:,:),repulsiveDerivs(:,:),totalDeriv(:,:)
    real(dp), allocatable :: chrgForces(:,:)

    !> excited state force addition
    real(dp), allocatable :: excitedDerivs(:,:)


    !> Stress tensors for various contribution in periodic calculations
    real(dp) :: elecStress(3,3), repulsiveStress(3,3), kineticStress(3,3)
    real(dp) :: dispStress(3,3), totalStress(3,3)


    !> Derivatives of lattice vectors in periodic calculations
    real(dp) :: elecLatDeriv(3,3), repulsiveLatDeriv(3,3)
    real(dp) :: dispLatDeriv(3,3), totalLatDeriv(3,3)

    !> derivative of cell volume wrt to lattice vectors, needed for pV term
    real(dp) :: derivCellVol(3,3)


    !> dipole moments when available
    real(dp) :: dipoleMoment(3)

    integer :: ii, jj

    logical :: tConverged

    logical, parameter :: tDensON2 = .false.  ! O(N^2) density mtx creation
    logical, parameter :: tAppendDetailedOut = .false.

    real(dp) :: cellPressure

    ! Variables for the geometry optimization

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
    real(dp), allocatable :: new3Coord(:,:)

    !> lattice vectors returned by the optimizer
    real(dp) :: tmpLatVecs(9), newLatVecs(9)
    real(dp) :: tmpLat3Vecs(3,3)

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
    integer :: iSCCIter, iSpin, iAtom

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
    real(dp), allocatable :: tmpDerivs(:)
    real(dp), allocatable :: tmpMatrix(:,:)
    real(dp), allocatable :: orbitalL(:,:,:), orbitalLPart(:,:,:)
    real(dp), allocatable :: rVecTemp(:)
    character(lc) :: lcTmp


    !> temporary character variable
    character(lc) :: tmpStr

    real(dp), pointer :: pDynMatrix(:,:)


    !> flag to write out geometries (and charge data if scc) when moving atoms about - in the case
    !> of conjugate gradient/steepest descent the geometries are written anyway
    logical :: tWriteRestart = .false.

    !> Minimal number of SCC iterations
    integer :: minSCCIter

    !> if scf/geometry driver should be stopped
    logical :: tStopSCC, tStopDriver

    !> net charge on each atom
    real(dp), allocatable :: dqAtom(:)


    !> density matrix
    real(dp), allocatable :: rhoSqrReal(:,:,:)


    !> Natural orbitals for excited state density matrix, if requested
    real(dp), allocatable :: naturalOrbs(:,:,:), occNatural(:,:)

    real(dp), allocatable :: invJacobian(:,:)


    !> locality measure for the wavefunction
    real(dp) :: localisation

    !> temporary variable for number of occupied levels
    integer :: nFilledLev

    call initOutputFiles(tWriteAutotest, tWriteResultsTag, tWriteBandDat, tDerivs,&
        & tWriteDetailedOut, tMd, fdAutotest, fdResultsTag, fdBand, fdEigvec, fdHessian, fdUser,&
        & fdMd)

    call initArrays(tForces, tExtChrg, tLinResp, tLinRespZVect, t3rdFull, tMd, tDerivs,&
      & tCoordOpt, tMulliken, tSpinOrbit, tDualSpinOrbit, tImHam, tStoreEigvecs, tWriteRealHS,&
      & tWriteHS, t2Component, tRealHS, tPrintExcitedEigvecs, orb, nAtom, nMovedAtom, nKPoint,&
      & nSpin, nExtChrg, forceType, indMovedAtom, mass, rhoPrim, h0, iRhoPrim, excitedDerivs,&
      & ERhoPrim, derivs, repulsiveDerivs, totalDeriv, chrgForces, energy, potential,&
      & shift3rd, orbResShift3rd, TS, E0, Eband, eigen, eigen2, filling, coord0Fold, new3Coord,&
      & tmpDerivs, orbitalL, orbitalLPart, HSqrCplx, HSqrCplx2, SSqrCplx, HSqrReal, HsqrReal2,&
      & SSqrReal, rhoSqrReal, dqAtom, chargePerShell, naturalOrbs, occNatural, velocities,&
      & movedVelo, movedAccel, movedMass)

    if (tShowFoldedCoord) then
      pCoord0Out => coord0Fold
    else
      pCoord0Out => coord0
    end if

    call initGeoOptParameters(tCoordOpt, nGeoSteps, tGeomEnd, tCoordStep, tCoordEnd, tStopDriver,&
        & iGeoStep, iLatGeoStep)

    minSccIter = getMinSccIters(tScc, tDftbU, nSpin)
    if (tXlbomd) then
      call xlbomdIntegrator%setDefaultSCCParameters(minSCCiter, maxSccIter, sccTol)
    end if

    tLatticeChanged = tPeriodic
    tCoordsChanged = .true.

    ! Belongs into initprogram where initial geometry is set up
    if (tSocket) then
      call receiveGeometryFromSocket(socket, tPeriodic, coord0, latVec, tCoordsChanged,&
          & tLatticeChanged)
    end if

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

      call initSccLoop(tScc, xlbomdIntegrator, minSccIter, maxSccIter, sccTol, tConverged, tStopScc)
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
        
        if (nSpin /= 4) then

          if (nSpin == 2) then
            call qm2ud(ham)
          end if

          if (tWriteRealHS .or. tWriteHS) then
            ! Write out matrices if necessary and quit.
            if (tImHam) then
              call writeHS(tWriteHS, tWriteRealHS, tRealHS, ham, over, &
                  & neighborList%iNeighbor, nNeighbor, iAtomStart, iPair, &
                  & img2CentCell, kPoint, iCellVec, cellVec, iHam=iHam)
            else
              call writeHS(tWriteHS, tWriteRealHS, tRealHS, ham, over, &
                  & neighborList%iNeighbor, nNeighbor, iAtomStart, iPair, &
                  & img2CentCell, kPoint, iCellVec, cellVec)
            end if
            write(stdOut, "(A)") "Hamilton/Overlap written, exiting program."
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
              call diagonalize(HSqrReal(:,:,iSpin2), SSqrReal, eigen(:,1,iSpin), ham(:,iSpin),&
                  & over, neighborList%iNeighbor, nNeighbor, iAtomStart, iPair, img2CentCell,&
                  & solver, 'V')
              if (tStoreEigvecs) then
                call reset(storeEigvecsReal(iSpin), [size(HSqrReal, dim=1), size(HSqrReal, dim=2)])
                call push(storeEigvecsReal(iSpin), HSqrReal(:,:,iSpin2))
              end if
            else
              do nk = 1, nKPoint
                if (tStoreEigvecs) then
                  iK2 = 1
                else
                  iK2 = nK
                end if
                call diagonalize(HSqrCplx(:,:,iK2,iSpin2), SSqrCplx, eigen(:,nk,iSpin),&
                    & ham(:,iSpin), over,kPoint(:,nk), neighborList%iNeighbor, nNeighbor,&
                    & iCellVec, cellVec, iAtomStart, iPair, img2CentCell, solver, 'V')
                if (tStoreEigvecs) then
                  call reset(storeEigvecsCplx(iSpin),&
                      & [size(HSqrCplx, dim=1), size(HSqrCplx, dim=2)])
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
            rhoPrim(:, iSpin) = 0.0_dp
            if (tImHam) then
              iRhoPrim(:, iSpin) = 0.0_dp
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

          if (tWriteBandDat) then
            call writeBandOut(fdBand, bandOut, eigen, filling, kWeight)
          end if


          call ud2qm(rhoPrim)

        else !  (nSpin == 4) then

          rhoPrim(:,:) = 0.0_dp
          if (tImHam) then
            iRhoPrim(:,:) = 0.0_dp
          end if


          if (tRealHS) then
            if (tImHam) then
              call diagonalize(HSqrCplx(:,:,1,1), SSqrCplx, eigen(:,1,1), ham, over,&
                  & neighborList%iNeighbor, nNeighbor, iAtomStart, iPair, img2CentCell, solver,&
                  & 'V', iHam=iHam)
            else if (tSpinOrbit) then
              call diagonalize(HSqrCplx(:,:,1,1), SSqrCplx, eigen(:,1,1), ham, over,&
                  & neighborList%iNeighbor, nNeighbor, iAtomStart, iPair, img2CentCell, solver,&
                  & 'V', xi,orb,species)
            else
              call diagonalize(HSqrCplx(:,:,1,1), SSqrCplx, eigen(:,1,1), ham, over,&
                  & neighborList%iNeighbor, nNeighbor, iAtomStart, iPair, img2CentCell, solver, 'V')
            end if
            if (tStoreEigvecs) then
              call reset(storeEigvecsCplx(1), [size(HSqrCplx, dim=1), size(HSqrCplx, dim=2)])
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
                    & kPoint(:,nk), neighborList%iNeighbor, nNeighbor, iCellVec, cellVec,&
                    & iAtomStart, iPair, img2CentCell, solver, 'V', iHam=iHam)
              else if (tSpinOrbit) then
                call diagonalize(HSqrCplx(:,:,iK2,1), SSqrCplx, eigen(:,nK,1), ham, over,&
                    & kPoint(:,nK), neighborList%iNeighbor, nNeighbor, iCellVec, cellVec,&
                    & iAtomStart, iPair, img2CentCell, solver, 'V', xi,orb,species)
              else
                call diagonalize(HSqrCplx(:,:,iK2,1), SSqrCplx, eigen(:,nK,1), ham, over,&
                    & kPoint(:,nK), neighborList%iNeighbor, nNeighbor, iCellVec, cellVec,&
                    & iAtomStart, iPair, img2CentCell, solver, 'V')
              end if
              if (tStoreEigvecs) then
                call reset(storeEigvecsCplx(1), [size(HSqrCplx, dim=1), size(HSqrCplx, dim=2)])
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
              allocate(rVecTemp(nAtom))
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
            deallocate(rVecTemp)
            energy%ELS = sum(energy%atomLS(:))
          end if
          filling(:,1:nKPoint,1) = 0.5_dp * filling(:,1:nKPoint,1)

          if (tWriteBandDat) then
            call writeBandOut(fdBand, bandOut, eigen, filling, kWeight)
          end if

        end if ! end of nSpin == 4 case


        ! Mulliken analysis


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
          call resetInternalPotentials(tDualSpinOrbit, xi, orb, species, potential)
          call getChargePerShell(qOutput, orb, species, chargePerShell)
          call addChargePotentials(qOutput, q0, chargePerShell, orb, species, species0,&
              & neighborList, img2CentCell, spinW, thirdOrd, potential)
          call addBlockChargePotentials(qBlockOut, qiBlockOut, tDftbU, tImHam, species, orb,&
              & nDftbUFunc, UJ, nUJ, iUJ, niUJ, potential)
          potential%intBlock = potential%intBlock + potential%extBlock
        end if

        ! Calculate energies

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
        energy%EMermin = energy%Etotal - sum(TS)
        energy%EGibbs = energy%EMermin + cellVol * pressure

        ! Stop SCC if appropriate stop file is present (We need this query here
        ! since the following block contains a check iSCCIter /= maxSccIter)
        inquire(file=fStopSCC, exist=tStopSCC)
        if (tStopSCC) then
          write(stdOut, "(3A)") "Stop file '" // fStopSCC // "' found."
          maxSccIter = iSCCIter
          write(stdOut, "(A)") "Setting max number of scc cycles to current cycle."
        end if

        ! Mix charges
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
          if ((.not. tConverged) .and. iSCCiter /= maxSccIter) then
            ! Avoid mixing of spin unpolarised density for spin polarised
            ! cases, this is only a problem in iteration 1, as there is
            ! only the (spin unpolarised!) atomic input density at that
            ! point. (Unless charges had been initialized externally)
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

        if (tSCC) then
          if (iSCCiter > 1) then
            diffElec = energy%Eelec - Eold
          else
            diffElec = 0.0_dp
          end if
          Eold = energy%Eelec
        end if

        if (tWriteDetailedOut) then
          call writeDetailedOut1(fdUser, userOut, tAppendDetailedOut, iDistribFn, nGeoSteps,&
              & iGeoStep, tMD, tDerivs, tCoordOpt, tLatOpt, iLatGeoStep, iSccIter, energy,&
              & diffElec, sccErrorQ, indMovedAtom, pCoord0Out, q0, qInput, qOutput, eigen, filling,&
              & orb, species, tDFTBU, tImHam, tPrintMulliken, orbitalL, qBlockOut, Ef, Eband, TS,&
              & E0, pressure, cellVol, tAtomicEnergy, tDispersion, tEField, tPeriodic, nSpin,&
              & tSpinOrbit, tScc)
        end if

        if (tSCC) then
          if (tDFTBU) then
            write(stdOut, "(I5,E18.8,E18.8,E18.8)") iSCCIter, energy%Eelec, diffElec, sccErrorQ
          else
            write(stdOut, "(I5,E18.8,E18.8,E18.8)") iSCCIter, energy%Eelec, diffElec, sccErrorQ
          end if
        end if

        ! Not writing any restarting info if not converged and minimal number of
        ! SCC iterations not done.
        if (restartFreq > 0 .and. .not.(tMD .or. tGeoOpt .or. tDerivs) .and. maxSccIter > 1) then
          if (tConverged .or. ((iSCCIter >= minSCCIter &
              & .or. tReadChrg .or. iGeoStep > 0) &
              &.and. (iSCCIter == maxSccIter .or. mod(iSCCIter, restartFreq) ==&
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

        if (tConverged) then
          exit lpSCC
        end if

      end do lpSCC

      ! Linear response
      energy%Eexcited = 0.0_dp
      if (tLinResp) then
        if (t3rd) then
          call error("Third order currently incompatible with excited state")
        end if
        if (.not.tRealHS) then
          call error("Only real systems are supported for excited state calculations")
        end if
        if (tPeriodic .and. tForces) then
          call error("Forces in the excited state for periodic geometries are currently&
              & unavailable")
        end if

        dqAtom = sum( qOutput(:,:,1) - q0(:,:,1) , dim=1)
        call unpackHS(SSqrReal, over, neighborList%iNeighbor, nNeighbor,&
            & iAtomStart, iPair, img2CentCell)
        call blockSymmetrizeHS(SSqrReal, iAtomStart)
        if (tForces) then
          @:ASSERT(.not. tPeriodic)
          do iSpin = 1, nSpin
            call blockSymmetrizeHS(rhoSqrReal(:,:,iSpin), iAtomStart)
          end do
        end if
        if (tWriteAutotest) then
          open(fdAutotest, file=autotestTag, position="append")
        end if

        if (tLinRespZVect) then
          if (tPrintExcitedEigVecs) then

            call addGradients(tSpin, lresp, iAtomStart, &
                & HSqrReal, eigen(:,1,:), SSqrReal, filling(:,1,:), coord0, &
                & dqAtom, species0, neighborList%iNeighbor, &
                & img2CentCell, orb, skHamCont, skOverCont, tWriteAutotest, &
                & fdAutotest, energy%Eexcited, excitedDerivs, &
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
                & img2CentCell, orb, skHamCont, skOverCont, tWriteAutotest, &
                & fdAutotest, energy%Eexcited, excitedDerivs, &
                & nonSccDeriv, rhoSqrReal)
          end if
        else
          call calcExcitations(tSpin, lresp, iAtomStart,&
              & HSqrReal, eigen(:,1,:), SSqrReal, filling(:,1,:), coord0,&
              & dqAtom, species0, neighborList%iNeighbor,&
              & img2CentCell, orb, tWriteAutotest, fdAutotest, energy%Eexcited)
        end if
        energy%Etotal = energy%Etotal + energy%Eexcited
        energy%EMermin = energy%EMermin + energy%Eexcited
        energy%EGibbs = energy%EGibbs + energy%Eexcited
        if (tWriteAutotest) then
          close(fdAutotest)
        end if
      end if

      if (tXlbomd) then
        if (xlbomdIntegrator%needsInverseJacobian()) then
          write(stdOut, "(A)") ">> Updating XLBOMD Inverse Jacobian"
          allocate(invJacobian(nIneqOrb, nIneqOrb))
          call getInverseJacobian(pChrgMixer, invJacobian)
          call xlbomdIntegrator%setInverseJacobian(invJacobian)
          deallocate(invJacobian)
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
          write (tmpStr, "(A, I0, A, I0)") '** Geometry step: ', iGeoStep, ', Lattice step: ',&
              & iLatGeoStep
        else
          write(tmpStr,"(A, I0)") 'Geometry Step: ', iGeoStep
        end if
        ! save geometry in gen format
        call writeGenGeometry(tGeoOpt, tMd, tWriteRestart, tFracCoord, tPeriodic, geoOutFile,&
            & pCoord0Out, species0, speciesName, latVec)

        if (tPrintMulliken) then
          if (nSpin == 4) then
            allocate(tmpMatrix(3,nAtom))
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
            deallocate(tmpMatrix)
          else
            call writeXYZFormat(trim(lcTmp), pCoord0Out, species0, speciesName, &
                &charges=sum(qOutput(:,:,1),dim=1),comment=trim(tmpStr))
          end if
        else
          call writeXYZFormat(trim(lcTmp), pCoord0Out, species0, speciesName, &
              &comment=trim(tmpStr))
        end if
      end if

      write(stdOut, *)
      write(stdOut, format2U) "Total Energy", energy%Etotal,"H", Hartree__eV * energy%Etotal,"eV"
      write(stdOut, format2U) "Total Mermin free energy", energy%EMermin, "H",&
          & Hartree__eV * energy%EMermin,"eV"

      if (tDipole) then
        dipoleMoment(:) = 0.0_dp
        do iAtom = 1, nAtom
          dipoleMoment(:) = dipoleMoment(:) &
              & + sum(q0(:, iAtom, 1) - qOutput(:, iAtom, 1)) * coord(:,iAtom)
        end do

#:call DEBUG_CODE
        ! extra test for the potential in the code, does the dipole from
        ! charge positions match the derivative of energy wrt an external E field?
        allocate(hprime(size(h0),1))
        allocate(dipoleTmp(size(qOutput,dim=1),nAtom))
        allocate(potentialDerivative(nAtom,1))
        write(stdOut, "(A)", advance='no') 'Hellmann Feynman dipole:'
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
          write(stdOut, "(F12.8)", advance='no') sum(dipoleTmp)
        end do
        write(stdOut, *) " au"
        deallocate(potentialDerivative)
        deallocate(hprime)
        deallocate(dipoleTmp)
#:endcall DEBUG_CODE
      else
        dipoleMoment(:) = 0.0_dp
      end if

      ! Calculate energy weighted density matrix
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

        ! Calculate the identity part of the energy weighted density matrix
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
                  ! Correct force for XLBOMD for T=0K (DHD)
                  ! Eigenvectors stored in HSqrReal are overwritten
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
                ! Calculate eigenvectors, if necessary. Eigenvectors for the last
                ! spin in the last k-points are still there, so use those
                ! directly
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
                    ! Correct force for XLBOMD for T=0K (DHD)
                    ! Eigenvectors stored in HSqrCplx are overwritten
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
          !if (tSCC) then
          !  potential%intBlock = potential%intBlock + potential%extBlock
          !else
          !  potential%intBlock = potential%extBlock
          !end if

          !if (tDFTBU) then
          !  potential%intBlock = potential%orbitalBlock + potential%intBlock
          !end if

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
            dispLatDeriv = -CellVol * matmul(dispStress, invLatVec)
          else
            dispStress(:,:) = 0.0_dp
            dispLatDeriv(:,:) = 0.0_dp
          end if

          if (tEField) then
            elecLatDeriv(:,:) = 0.0_dp
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

          write(stdOut, format2Ue) 'Volume', cellVol, 'au^3', (Bohr__AA**3) * cellVol, 'A^3'

        end if

      end if

      ! MD case includes atomic kinetic energy contribution, so print that later
      if (tStress .and. .not. tMD) then
        write(stdOut, format2Ue) 'Pressure', cellPressure, 'au', cellPressure * au__pascal, 'Pa'
        if (pressure /= 0.0_dp) then
          write(stdOut, format2U) "Gibbs free energy", energy%EGibbs, 'H',&
              & Hartree__eV * energy%EGibbs, 'eV'
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
        ! Set force components along constraint vectors zero
        do ii = 1, nGeoConstr
          iAtom = conAtom(ii)
          totalDeriv(:,iAtom) = totalDeriv(:,iAtom) &
              &- conVec(:,ii) * dot_product(conVec(:,ii), totalDeriv(:,iAtom))
        end do

        if (tCoordOpt) then
          tmpDerivs(1:nMovedCoord) = reshape(totalDeriv(:,indMovedAtom),(/ nMovedCoord /))
          write(stdOut, "(A, ':', T30, E20.6)") "Maximal force component", maxval(abs(tmpDerivs))
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
          write(stdOut, format1Ue) "Maximal Lattice force component", maxval(abs(tmpLatVecs)), 'au'
        end if

        if (tSocket) then
          ! stress was computed above in the force evaluation block or is 0 if aperiodic
          call socket%send(energy%ETotal - sum(TS), -totalDeriv, &
              & totalStress * cellVol)
        end if

        ! If geometry minimizer finished and the last calculated geometry is the
        ! minimal one (not necessary the case, depends on the optimizer!)
        ! -> we are finished.
        ! Otherwise we have to recalc everything in the converged geometry.

        if (tGeomEnd) then
          exit lpGeomOpt
        else
          tCoordsChanged = .false.
          tLatticeChanged = .false.
          if (tWriteRestart .and. tMulliken .and. tSCC .and. .not. tDerivs .and. maxSccIter > 1)&
              & then
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
            call next(derivDriver, new3Coord, totalDeriv(:,indMovedAtom), tGeomEnd)
            coord0(:,indMovedAtom) = new3Coord(:,:)
            if (tGeomEnd) exit lpGeomOpt
            tCoordsChanged = .true.
          else if (tGeoOpt) then
            if (tCoordStep) then
              call next(pGeoCoordOpt, energy%EMermin, tmpDerivs, tmpCoords,tCoordEnd)
              tCoordsChanged = .true.
              if (.not.tLatOpt) tGeomEnd = tCoordEnd
            else
              call next(pGeoLatOpt, energy%EGibbs, tmpLatVecs, newLatVecs,tGeomEnd)
              tLatticeChanged = .true.
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
          else if (tMD) then
            movedAccel(:,:) = -totalDeriv(:,indMovedAtom) / movedMass
            call next(pMDIntegrator, movedAccel ,new3Coord, movedVelo)
            tCoordsChanged = .true.
            if (allocated(temperatureProfile)) then
              call next(temperatureProfile)
            end if
            call evalKE(energy%Ekin, movedVelo, movedMass(1,:))
            call evalkT(pMDFrame, kT, movedVelo, movedMass(1,:))
            velocities(:,:) = 0.0_dp
            velocities(:, indMovedAtom) = movedVelo(:,:)
            energy%EMerminKin = energy%EMermin + energy%Ekin
            energy%EGibbsKin = energy%EGibbs + energy%Ekin
            if (tWriteRestart) then
              write(tmpStr, "(A, I0)") 'MD iter: ', iGeoStep
              ! save geometry in gen format
              call writeGenGeometry(tGeoOpt, tMd, tWriteRestart, tFracCoord, tPeriodic, geoOutFile,&
                  & pCoord0Out, species0, speciesName, latVec)
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
              call getKineticStress(kineticStress, mass, species0, velocities, CellVol)
              totalStress = totalStress + kineticStress
              cellPressure = ( totalStress(1,1) + totalStress(2,2) &
                  & + totalStress(3,3) )/3.0_dp

              totalLatDeriv = -CellVol * matmul(totalStress,invLatVec)
            end if

            if (tSetFillingTemp) then
              write(stdOut, format2U) 'Electronic Temperature:', tempElec, 'H',&
                  & tempElec/Boltzmann, 'K'
            end if
            if (tEfield) then
              write(stdOut, format1U1e) 'External E field', absEField, 'au',&
                  & absEField * au__V_m, 'V/m'
            end if
            write(stdOut, format2U) "MD Temperature:", kT, "H", kT / Boltzmann, "K"
            write(stdOut, format2U) "MD Kinetic Energy", energy%Ekin, "H",&
                & Hartree__eV * energy%Ekin, "eV"
            write(stdOut, format2U) "Total MD Energy", energy%EMerminKin, "H", &
                & Hartree__eV * energy%EMerminKin, "eV"
            if (tPeriodic) then
              write(stdOut, format2Ue) 'Pressure', cellPressure, 'au',&
                  & cellPressure * au__pascal, 'Pa'
              if (pressure /= 0.0_dp) then
                write(stdOut, format2U) 'Gibbs free energy including KE', energy%EGibbsKin, 'H',&
                    & Hartree__eV * energy%EGibbsKin,'eV'
              end if
            end if

            if (tWriteDetailedOut) then
              call writeDetailedOut3(fdUser, tPrintForces, tSetFillingTemp, tPeriodic, tStress,&
                  & totalStress, totalLatDeriv, energy, tempElec, pressure, cellPressure, kT)
            end if
          else if (tSocket .and. iGeoStep < nGeoSteps) then
            ! Only receive geometry from socket, if there are still geometry iterations left
            call receiveGeometryFromSocket(socket, tPeriodic, coord0, latVec, tCoordsChanged,&
                & tLatticeChanged)
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
                  ! Workaround NAG (6.1)
                  !coord0(:,indMovedAtom) = reshape(tmpCoords, [3, nMovedAtom])
                  do ii = 1, nMovedAtom
                    coord0(:,indMovedAtom(ii)) = tmpCoords((ii-1)*3+1:ii*3)
                  end do                  
                end if
              else
                call cart2frac(coord0, latVec)
                latVec(:,:) = reshape(newLatVecs, [3, 3])
                call frac2cart(coord0, latVec)
                tCoordsChanged = .true.
                if (tCoordOpt) then
                  tCoordStep = .true.
                  tCoordEnd = .false.
                  tmpCoords(1:nMovedCoord) = reshape(coord0(:, indMovedAtom), [nMovedCoord])
                  call reset(pGeoCoordOpt, tmpCoords)
                end if
              end if
            elseif (tMD) then

              coord0(:,indMovedAtom) = new3Coord(:,:)

              if (tBarostat) then ! apply a Barostat
                call rescale(pMDIntegrator,coord0,latVec,totalStress)
                tCoordsChanged = .true.
                tLatticeChanged = .true.
              end if

              if (tWriteRestart) then
                if (tPeriodic) then
                  cellVol = abs(determinant33(latVec))
                  energy%EGibbs = energy%EMermin + pressure * cellVol
                end if
                call writeMdOut2(fdMd, tStress, tBarostat, tLinResp, tEField, tFixEf,&
                    & tPrintMulliken, tDipole, energy, latVec, cellVol, cellPressure, pressure, kT,&
                    & absEField, dipoleMoment, qOutput, q0)
              end if
            end if
          end if
        end if
      end if
      
      if (tWriteDetailedOut) then
        call writeDetailedOut4(fdUser, tMD, energy, kT)
      end if

      ! Stop reading of initial charges/block populations again
      tReadChrg = .false.

      ! Stop SCC if appropriate stop file is present
      if (.not. tStopSCC) then
        inquire(file=fStopDriver, exist=tStopDriver)
        if (tStopDriver) then
          write(stdOut, "(3A)") "Stop file '" // fStopDriver // "' found."
        end if
      end if
      if (tStopSCC .or. tStopDriver) then
        nGeoSteps = iGeoStep
        write(stdOut, "(A)") "Setting max number of geometry steps to current step number."
      end if

    end do lpGeomOpt

    if (tSocket) then
      call socket%shutdown()
    end if

    tGeomEnd = tMD .or. tGeomEnd .or. tDerivs

    if (tWriteDetailedOut) then
      call writeDetailedOut5(fdUser, tGeoOpt, tGeomEnd, tMd, tDerivs, tEField, tDipole, absEField,&
          & dipoleMoment)
    end if

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

    if (tMD) then
      call writeMdOut3(fdMd)
      write(stdOut, "(2A)") 'MD information accumulated in ', mdOut
    end if

    if (tDerivs) then
      call getHessianMatrix(derivDriver, pDynMatrix)
      write(stdOut, "(2A)") 'Hessian matrix written to ', hessianOut
      call writeHessianOut(fdHessian, hessianOut, pDynMatrix)
    else
      nullify(pDynMatrix)
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

            write(stdOut, *) 'Original localisation', localisation

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

            write(stdOut, "(A, E20.12)") 'Final localisation ', localisation

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

            write(stdOut, "(A, E20.12)") 'Original localisation', localisation

            call PipekMezey(HSqrCplx(:,:nFilledLev,:,iSpin), SSqrCplx, &
                & over, kpoint, kweight, neighborList%iNeighbor, nNeighbor, &
                & iCellVec, cellVec, iAtomStart, iPair, img2CentCell, PipekTol, &
                & PipekMaxIter)

            localisation = sum(PipekMezeyLocalisation( &
                & HSqrCplx(:,:nFilledLev,:,iSpin), SSqrCplx, over, kpoint, &
                & kweight, neighborList%iNeighbor, nNeighbor, iCellVec, cellVec, &
                & iAtomStart, iPair, img2CentCell))

            write(stdOut, "(A, E20.12)") 'Final localisation', localisation

          end do

          call writeEigvecs(fdEigvec, runId, nAtom, nSpin, neighborList, &
              & nNeighbor, cellVec, iCellVec, iAtomStart, iPair, img2CentCell, &
              & orb, species, speciesName, over, kpoint, HSqrCplx, SSqrCplx, &
              & fileName="localOrbs")

        end if

      end if
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
      call writeResultsTag(fdResultsTag, resultsTag, energy, tAtomicEnergy, totalDeriv, chrgForces,&
          & tStress, totalStress, pDynMatrix, tScc, iSccIter, tConverged, tPrintMulliken, qOutput,&
          & q0, eigen, filling, Ef, nEl, tPeriodic, cellVol)
    end if
    if (tWriteDetailedXML) then
      call writeDetailedXml(runId, speciesName, species0, pCoord0Out, tPeriodic, latVec, tRealHS,&
          & nKPoint, nSpin, size(eigen, dim=1), nOrb, kPoint, kWeight, filling, occNatural)
    end if

    call destructProgramVariables()

  end subroutine runDftbPlus


  subroutine initOutputFiles(tWriteAutotest, tWriteResultsTag, tWriteBandDat, tDerivs,&
      & tWriteDetailedOut, tMd, fdAutotest, fdResultsTag, fdBand, fdEigvec, fdHessian, fdUser, fdMd)
    logical, intent(in) :: tWriteAutotest, tWriteResultsTag, tWriteBandDat, tDerivs
    logical, intent(in) :: tWriteDetailedOut, tMd
    integer, intent(out) :: fdAutotest, fdResultsTag, fdBand, fdEigvec, fdHessian, fdUser, fdMd
    
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

  end subroutine initOutputFiles


  subroutine initArrays(tForces, tExtChrg, tLinResp, tLinRespZVect, t3rdFull, tMd, tDerivs,&
      & tCoordOpt, tMulliken, tSpinOrbit, tDualSpinOrbit, tImHam, tStoreEigvecs, tWriteRealHS,&
      & tWriteHS, t2Component, tRealHS, tPrintExcitedEigvecs, orb, nAtom, nMovedAtom, nKPoint,&
      & nSpin, nExtChrg, forceType, indMovedAtom, mass, rhoPrim, h0, iRhoPrim, excitedDerivs,&
      & ERhoPrim, derivs, repulsiveDerivs, totalDeriv, chrgForces, energy, potential,&
      & shift3rd, orbResShift3rd, TS, E0, Eband, eigen, eigen2, filling, coord0Fold, new3Coord,&
      & tmpDerivs, orbitalL, orbitalLPart, HSqrCplx, HSqrCplx2, SSqrCplx, HSqrReal, HsqrReal2,&
      & SSqrReal, rhoSqrReal, dqAtom, chargePerShell, naturalOrbs, occNatural, velocities,&
      & movedVelo, movedAccel, movedMass)
    logical, intent(in) :: tForces, tExtChrg, tLinResp, tLinRespZVect, t3rdFull, tMd, tDerivs
    logical, intent(in) :: tCoordOpt, tMulliken, tSpinOrbit, tDualSpinOrbit, tImHam, tStoreEigvecs
    logical, intent(in) :: tWriteRealHS, tWriteHS, t2Component, tRealHS, tPrintExcitedEigvecs
    type(TOrbitals), intent(in) :: orb
    integer, intent(in) :: nAtom, nMovedAtom, nKPoint, nSpin, nExtChrg, forceType
    integer, intent(in) :: indMovedAtom(:)
    real(dp), intent(in) :: mass(:)
    real(dp), intent(out), allocatable :: rhoPrim(:,:), h0(:), iRhoPrim(:,:), excitedDerivs(:,:)
    real(dp), intent(out), allocatable :: ERhoPrim(:)
    real(dp), intent(out), allocatable :: derivs(:,:), repulsiveDerivs(:,:), totalDeriv(:,:)
    real(dp), intent(out), allocatable :: chrgForces(:,:)
    type(TEnergies), intent(out) :: energy
    type(TPotentials), intent(out) :: potential
    real(dp), intent(out), allocatable :: shift3rd(:), orbResShift3rd(:,:), TS(:), E0(:), Eband(:)
    real(dp), intent(out), allocatable :: eigen(:,:,:), eigen2(:,:,:), filling(:,:,:)
    real(dp), intent(out), allocatable :: coord0Fold(:,:), new3Coord(:,:), tmpDerivs(:)
    real(dp), intent(out), allocatable :: orbitalL(:,:,:), orbitalLPart(:,:,:)
    complex(dp), intent(out), allocatable :: HSqrCplx(:,:,:,:), HSqrCplx2(:,:), SSqrCplx(:,:)
    real(dp), intent(out), allocatable :: HSqrReal(:,:,:), HSqrReal2(:,:), SSqrReal(:,:)
    real(dp), intent(out), allocatable :: rhoSqrReal(:,:,:)
    real(dp), intent(out), allocatable :: dqAtom(:), chargePerShell(:,:,:)
    real(dp), intent(out), allocatable :: naturalOrbs(:,:,:), occNatural(:,:)
    real(dp), intent(out), allocatable :: velocities(:,:), movedVelo(:,:), movedAccel(:,:)
    real(dp), intent(out), allocatable :: movedMass(:,:)

    integer :: nSpinHams, sqrHamSize
    integer :: nSpinStored, nKPointStored

    allocate(rhoPrim(0, nSpin))
    allocate(h0(0))
    if (tImHam) then
      allocate(iRhoPrim(0, nSpin))
    end if

    if (tForces) then
      allocate(ERhoPrim(0))
      allocate(derivs(3, nAtom))
      allocate(repulsiveDerivs(3, nAtom))
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

    if (t3rdFull) then
      allocate(shift3rd(nAtom))
      allocate(orbresshift3rd(orb%mShell,nAtom))
    end if

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
    allocate(eigen2(sqrHamSize, nKPoint, nSpinHams))
    allocate(filling(sqrHamSize, nKpoint, nSpinHams))

    allocate(coord0Fold(3, nAtom))

    if (tMD.or.tDerivs) then
      allocate(new3Coord(3, nMovedAtom))
    end if

    if (tCoordOpt) then
      allocate(tmpDerivs(3 * nMovedAtom))
    end if

    if ((tMulliken .and. tSpinOrbit) .or. tImHam) then
      allocate(orbitalL(3, orb%mShell, nAtom))
    end if
    if ((tMulliken .and. tSpinOrbit) .and. .not.  tDualSpinOrbit) then
      allocate(orbitalLPart(3, orb%mShell, nAtom))
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
        if (any(forceType == [ 1, 2, 3 ])) then
          allocate(HSqrReal2(sqrHamSize, sqrHamSize))
        end if
        allocate(SSqrReal(sqrHamSize, sqrHamSize))
      else
        allocate(HSqrCplx(sqrHamSize, sqrHamSize, nKPointStored, nSpinStored))
        if (any(forceType == [ 1, 2, 3 ])) then
          allocate(HSqrCplx2(sqrHamSize, sqrHamSize))
        end if
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
      allocate(naturalOrbs(orb%nOrb, orb%nOrb, 1))
      allocate(occNatural(orb%nOrb, 1))
    end if

    if (tMD) then
      allocate(velocities(3, nAtom))
      allocate(movedVelo(3, nMovedAtom))
      allocate(movedAccel(3, nMovedAtom))
      allocate(movedMass(3, nMovedAtom))
      movedMass(:,:) = spread(mass(indMovedAtom),1,3)
    end if

  end subroutine initArrays


  !> Initialises some parameters before geometry loop starts.
  subroutine initGeoOptParameters(tCoordOpt, nGeoSteps, tGEomEnd, tCoordStep, tCoordEnd,&
      & tStopDriver, iGeoStep, iLatGeoStep)
    logical, intent(in) :: tCoordOpt
    integer, intent(in) :: nGeoSteps
    logical, intent(out) :: tGeomEnd, tCoordStep, tCoordEnd, tStopDriver
    integer, intent(out) :: iGeoStep, iLatGeoStep
    
    tGeomEnd = (nGeoSteps == 0)

    tCoordStep = .false.
    if (tCoordOpt) then
      tCoordStep = .true.
      tCoordEnd = .false.
    end if

    iGeoStep = 0
    iLatGeoStep = 0
    tStopDriver = .false.

  end subroutine initGeoOptParameters


  !> Receives geometry from socket.
  subroutine receiveGeometryFromSocket(socket, tPeriodic, coord0, latVecs, tCoordsChanged,&
      & tLatticeChanged)
    type(IpiSocketComm), allocatable, intent(inout) :: socket
    logical, intent(in) :: tPeriodic
    real(dp), intent(out) :: coord0(:,:)
    real(dp), intent(out) :: latVecs(:,:)
    logical, intent(out) :: tCoordsChanged, tLatticeChanged

    real(dp) :: tmpLatVecs(3, 3)
    
    call socket%receive(coord0, tmpLatVecs)
    tCoordsChanged = .true.
    if (tPeriodic) then
      latVecs(:,:) = tmpLatVecs
    end if
    tLatticeChanged = tPeriodic

  end subroutine receiveGeometryFromSocket


  !> Initialises SCC related parameters before geometry loop starts
  function getMinSccIters(tScc, tDftbU, nSpin) result(minSccIter)
    logical, intent(in) :: tScc, tDftbU
    integer, intent(in) :: nSpin
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


  !> Does the steps necessary after a lattice vector update
  subroutine handleLatticeChange(latVecs, tScc, tStress, pressure, mCutoff, dispersion,&
      & recVecs, recVecs2p, cellVol, recCellVol, derivCellVol, cellVecs, rCellVecs)
    real(dp), intent(in) :: latVecs(:,:)
    logical, intent(in) :: tScc, tStress
    real(dp), intent(in) :: pressure
    real(dp), intent(inout) :: mCutoff
    class(DispersionIface), allocatable, intent(inout) :: dispersion
    real(dp), intent(out) :: recVecs(:,:), recVecs2p(:,:)
    real(dp), intent(out) :: cellVol, recCellVol
    real(dp), intent(out) :: derivCellVol(:,:)
    real(dp), allocatable, intent(out) :: cellVecs(:,:), rCellVecs(:,:)

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


  !> Does necessary steps when ionic coordinates change
  subroutine handleCoordinateChange(coord0, latVec, invLatVec, species0, mCutoff, skRepCutoff, &
      & orb, tPeriodic, tScc, tDispersion, dispersion, thirdOrd, img2CentCell, iCellVec,&
      & neighborList, nAllAtom, coord0Fold, coord, species, rCellVec, nAllOrb, nNeighbor, ham,&
      & over, H0, rhoPrim, iRhoPrim, iHam, ERhoPrim, iPair)
    real(dp), intent(in) :: coord0(:,:), latVec(:,:), invLatVec(:,:)
    integer, intent(in) :: species0(:)
    real(dp), intent(in) :: mCutoff, skRepCutoff
    type(TOrbitals), intent(in) :: orb
    logical, intent(in) :: tPeriodic, tScc, tDispersion
    class(DispersionIface), allocatable, intent(inout) :: dispersion
    type(ThirdOrder), allocatable, intent(inout) :: thirdOrd
    integer, allocatable, intent(inout) :: img2CentCell(:), iCellVec(:)
    type(TNeighborList), intent(inout) :: neighborList
    integer, intent(out) :: nAllAtom, nAllOrb
    real(dp), intent(out) :: coord0Fold(:,:)
    real(dp), allocatable, intent(inout) :: coord(:,:)
    integer, allocatable, intent(inout) :: species(:)
    real(dp), allocatable, intent(inout) :: rCellVec(:,:)
    integer, intent(out) :: nNeighbor(:)
    real(dp), allocatable, intent(inout) :: ham(:,:), over(:), h0(:)
    real(dp), allocatable, intent(inout) :: rhoPrim(:,:), iRhoPrim(:,:), iHam(:,:), ERhoPrim(:)
    integer, allocatable, intent(inout) :: iPair(:,:)

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
    logical, intent(in) :: tGeoOpt, tMd
    integer, intent(in) :: iGeoStep, nGeoSteps, restartFreq
    logical :: tWriteRestart
    
    if (restartFreq > 0 .and. (tGeoOpt .or. tMD)) then
      tWriteRestart = (iGeoStep == nGeoSteps .or. (mod(iGeoStep, restartFreq) == 0))
    else
      tWriteRestart = .false.
    end if

  end function needsRestartWriting


  !> Ensures that sparse array have enough storage to hold all necessary elements.
  subroutine reallocateSparseArrays(sparseSize, ham, over, H0, rhoPrim, iHam, iRhoPrim, ERhoPrim)
    integer, intent(in) :: sparseSize
    real(dp), allocatable, intent(inout) :: ham(:,:), over(:), H0(:), rhoPrim(:,:)
    real(dp), allocatable, intent(inout) :: iHam(:,:), iRhoPrim(:,:), ERhoPrim(:)

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
    real(dp), intent(in) :: coord(:,:)
    integer, intent(in) :: species(:), img2CentCell(:), nNeighbor(:)
    type(TNeighborList), intent(in) :: neighborList
    type(ORepCont), intent(in) :: pRepCont
    real(dp), intent(out) :: Eatom(:), Etotal
    
    call getERep(Eatom, coord, nNeighbor, neighborList%iNeighbor, species, pRepCont, img2CentCell)
    Etotal = sum(Eatom)

  end subroutine calcRepulsiveEnergy


  !> Calculates dispersion energy for current geometry.
  subroutine calcDispersionEnergy(dispersion, Eatom, Etotal)
    class(DispersionIface), intent(inout) :: dispersion
    real(dp), intent(out) :: Eatom(:), Etotal

    call dispersion%getEnergies(Eatom)
    Etotal = sum(Eatom)

  end subroutine calcDispersionEnergy


  !> Sets the external potential components to zero
  subroutine resetExternalPotentials(potential)
    type(TPotentials), intent(inout) :: potential
    
    potential%extAtom(:,:) = 0.0_dp
    potential%extShell(:,:,:) = 0.0_dp
    potential%extBlock(:,:,:,:) = 0.0_dp

  end subroutine resetExternalPotentials


  !> Merges atomic and shell resolved external potentials into blocked one
  subroutine mergeExternalPotentials(orb, species, potential)
    type(TOrbitals), intent(in) :: orb
    integer, intent(in) :: species(:)
    type(TPotentials), intent(inout) :: potential
    
    call total_shift(potential%extShell, potential%extAtom, orb, species)
    call total_shift(potential%extBlock, potential%extShell, orb, species)

  end subroutine mergeExternalPotentials


  !> Sets up electric external field
  subroutine setUpExternalElectricField(tTimeDepEField, tPeriodic, EFieldStrength, EFieldVector,&
      & EFieldOmega, EFieldPhase, neighborList, nNeighbor, iCellVec, img2CentCell, cellVec, deltaT,&
      & iGeoStep, coord0Fold, coord, EField, extAtomPot, absEField)
    logical, intent(in) :: tTimeDepEField, tPeriodic
    real(dp), intent(in) :: EFieldStrength, EFieldVector(:), EFieldOmega
    integer, intent(in) :: EFieldPhase
    type(TNeighborList), intent(in) :: neighborList
    integer, intent(in) :: nNeighbor(:), iCellVec(:), img2CentCell(:)
    real(dp), intent(in) :: cellVec(:,:)
    real(dp), intent(in) :: deltaT
    integer, intent(in) :: iGeoStep
    real(dp), allocatable, intent(in) :: coord0Fold(:,:)
    real(dp), intent(in) :: coord(:,:)
    real(dp), intent(out) :: EField(:), extAtomPot(:)
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
  subroutine initSccLoop(tScc, xlbomdIntegrator, minSccIter, maxSccIter, sccTol, tConverged,&
      & tStopScc)
    logical, intent(in) :: tScc
    type(Xlbomd), allocatable, intent(inout) :: xlbomdIntegrator
    integer, intent(inout) :: minSccIter, maxSccIter
    real(dp), intent(inout) :: sccTol
    logical, intent(out) :: tConverged, tStopScc
    
    if (allocated(xlbomdIntegrator)) then
      call xlbomdIntegrator%getSCCParameters(minSCCIter, maxSccIter, sccTol)
    end if

    tConverged = (.not. tScc)
    tStopSCC = .false.

    if (tScc) then
      call printSccHeader()
    end if
    
  end subroutine initSccLoop


  !> Calculate the number of charges per shell.
  subroutine getChargePerShell(qq, orb, species, chargePerShell)
    real(dp), intent(in) :: qq(:,:,:)
    type(TOrbitals), intent(in) :: orb
    integer, intent(in) :: species(:)
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
    logical, intent(in) :: tDualSpinOrbit
    real(dp), intent(in) :: xi(:,:)
    type(TOrbitals), intent(in) :: orb
    integer, intent(in) :: species(:)
    type(TPotentials), intent(inout) :: potential

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
    real(dp), intent(in) :: qInput(:,:,:), q0(:,:,:), chargePerShell(:,:,:)
    type(TOrbitals), intent(in) :: orb
    integer, intent(in) :: species(:), species0(:)
    type(TNeighborList), intent(in) :: neighborList
    integer, intent(in) :: img2CentCell(:)
    real(dp), intent(in), allocatable :: spinW(:,:,:)
    type(ThirdOrder), allocatable, intent(inout) :: thirdOrd
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
    real(dp), allocatable, intent(in) :: qBlockIn(:,:,:,:), qiBlockIn(:,:,:,:)
    logical, intent(in) :: tDftbU, tImHam
    integer, intent(in) :: species(:)
    type(TOrbitals), intent(in) :: orb
    integer, intent(in) :: nDftbUFunc
    real(dp), allocatable, intent(in) :: UJ(:,:)
    integer, allocatable, intent(in) :: nUJ(:), iUJ(:,:,:), niUJ(:,:)
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
    real(dp), intent(in) :: H0(:), over(:)
    integer, intent(in) :: nNeighbor(:)
    type(TNeighborList), intent(in) :: neighborList
    integer, intent(in) :: species(:)
    type(TOrbitals), intent(in) :: orb
    integer, intent(in) :: iPair(:,:), img2CentCell(:)
    type(TPotentials), intent(in) :: potential
    real(dp), intent(out) :: ham(:,:)
    real(dp), allocatable, intent(inout) :: iHam(:,:)

    integer :: nAtom, nSpin

    nAtom = size(orb%nOrbAtom)
    nSpin = size(ham, dim=2)

    ham(:,:) = 0.0_dp
    ham(:, 1) = h0
    call add_shift(ham, over, nNeighbor, neighborList%iNeighbor, species, orb, iPair, nAtom,&
        & img2CentCell, potential%intBlock)
    ! hack due to not using Pauli-type structure for diagonalisation
    if (nSpin > 1) then
      ham(:,:) = 2.0_dp * ham
    end if

    if (allocated(iHam)) then
      iHam(:,:) = 0.0_dp
      call add_shift(iHam, over, nNeighbor, neighborList%iNeighbor, species, orb, iPair, nAtom,&
          & img2CentCell, potential%iorbitalBlock)
      iHam(:,:) = 2.0_dp * iHam
    end if
    
  end subroutine getSccHamiltonian



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


end module main
