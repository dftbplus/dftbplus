!--------------------------------------------------------------------------------------------------!
!  DFTB+: general package for performing fast atomistic simulations                                !
!  Copyright (C) 2018  DFTB+ developers group                                                      !
!                                                                                                  !
!  See the LICENSE file for terms of usage and distribution.                                       !
!--------------------------------------------------------------------------------------------------!

#:include 'common.fypp'

!> Module for linear response derivative calculations using perturbation methods
module dftbp_perturbderivs
  use dftbp_accuracy
  use dftbp_constants
  use dftbp_globalenv
  use dftbp_message
  use dftbp_commontypes
  use dftbp_potentials
  use dftbp_scc
  use dftbp_orbitalequiv
  use dftbp_populations
  use dftbp_spin
  use dftbp_thirdorder_module, only : ThirdOrder
  use dftbp_dftbplusu
  use dftbp_rangeseparated, only : RangeSepFunc
  use dftbp_onsitecorrection
  use dftbp_mainio
  use dftbp_shift
  use dftbp_mixer
  use dftbp_finitethelper
  use dftbp_scalapackfx
  use dftbp_environment
  use dftbp_periodic
  use dftbp_densedescr
  use dftbp_sparse2dense
  use dftbp_taggedoutput
  use dftbp_degeneratePerturb
#:if WITH_MPI
  use dftbp_mpifx
#:endif
#:if WITH_SCALAPACK
  use dftbp_scalafxext
#:else
  use dftbp_blasroutines
#:endif

  implicit none

  private
  public :: staticPerturWrtE

  !> small complex value for frequency dependent cases
  complex(dp), parameter :: eta = (0.0_dp,1.0E-8_dp)

  !> Direction labels
  character(len=1), parameter :: direction(3) = ['x','y','z']

contains

  !> Static (frequency independent) perturbation at q=0
  subroutine staticPerturWrtE(env, parallelKS, filling, eigvals, eigVecsReal, eigVecsCplx, ham,&
      & over, orb, nAtom, species, speciesnames, neighbourList, nNeighbourSK, denseDesc,&
      & iSparseStart, img2CentCell, coord, sccCalc, maxSccIter, sccTol, nMixElements,&
      & nIneqMixElements, iEqOrbitals, tempElec, Ef, tFixEf, spinW, thirdOrd, tDFTBU, UJ, nUJ,&
      & iUJ, niUJ, iEqBlockDftbu, onsMEs, iEqBlockOnSite, rangeSep, nNeighbourLC, pChrgMixer,&
      & taggedWriter, tWriteAutoTest, autoTestTagFile, tWriteTaggedOut, taggedResultsFile,&
      & tWriteDetailedOut, fdDetailedOut, kPoint, kWeight, iCellVec, cellVec, tPeriodic, omegas)

    !> Environment settings
    type(TEnvironment), intent(inout) :: env

    !> K-points and spins to process
    type(TParallelKS), intent(in) :: parallelKS

    !> Filling
    real(dp), intent(in) :: filling(:,:,:)

    !> Eigenvalue of each level, kpoint and spin channel
    real(dp), intent(in) :: eigvals(:,:,:)

    !> ground state eigenvectors
    real(dp), intent(in), allocatable :: eigVecsReal(:,:,:)

    !> ground state complex eigenvectors
    complex(dp), intent(in), allocatable :: eigvecsCplx(:,:,:)

    !> Sparse Hamiltonian
    real(dp), intent(in) :: ham(:,:)

    !> sparse overlap matrix
    real(dp), intent(in) :: over(:)

    !> Atomic orbital information
    type(TOrbitals), intent(in) :: orb

    !> Number of central cell atoms
    integer, intent(in) :: nAtom

    !> chemical species
    integer, intent(in) :: species(:)

    !> label for each atomic chemical species
    character(mc), intent(in) :: speciesnames(:)

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

    !> atomic coordinates
    real(dp), intent(in) :: coord(:,:)

    !> SCC module internal variables
    type(TScc), intent(inout) :: sccCalc

    !> maximal number of SCC iterations
    integer, intent(in) :: maxSccIter

    !> Tolerance for SCC convergence
    real(dp), intent(in) :: sccTol

    !> nr. of elements to go through the mixer - may contain reduced orbitals and also orbital
    !> blocks (if tDFTBU or onsite corrections)
    integer, intent(in) :: nMixElements

    !> nr. of inequivalent charges
    integer, intent(in) :: nIneqMixElements

    !> Equivalence relations between orbitals
    integer, intent(in) :: iEqOrbitals(:,:,:)

    !> onsite matrix elements for shells (elements between s orbitals on the same shell are ignored)
    real(dp), intent(in), allocatable :: onsMEs(:,:,:,:)

    !> Equivalences for onsite block corrections if needed
    integer, intent(in), allocatable :: iEqBlockOnSite(:,:,:,:)

    !> Electron temperature
    real(dp), intent(in) :: tempElec

    !> Fermi level(s)
    real(dp), intent(in) :: Ef(:)

    !> Whether fixed Fermi level(s) should be used. (No charge conservation!)
    logical, intent(in) :: tFixEf

    !> spin constants
    real(dp), intent(in), allocatable :: spinW(:,:,:)

    !> Third order SCC interactions
    type(ThirdOrder), allocatable, intent(inout) :: thirdOrd

    !> is this a +U calculation
    logical, intent(in) :: tDftbU

    !> prefactor for +U potential
    real(dp), allocatable, intent(in) :: UJ(:,:)

    !> Number DFTB+U blocks of shells for each atom type
    integer, intent(in), allocatable :: nUJ(:)

    !> which shells are in each DFTB+U block
    integer, intent(in), allocatable :: iUJ(:,:,:)

    !> Number of shells in each DFTB+U block
    integer, intent(in), allocatable :: niUJ(:,:)

    !> equivalence mapping for dual charge blocks
    integer, intent(in), allocatable :: iEqBlockDftbu(:,:,:,:)

    !> Data for range-separated calculation
    type(RangeSepFunc), allocatable, intent(inout) :: rangeSep

    !> Number of neighbours for each of the atoms for the exchange contributions in the long range
    !> functional
    integer, intent(inout), allocatable :: nNeighbourLC(:)

    !> Charge mixing object
    type(OMixer), intent(inout) :: pChrgMixer

    !> Tagged writer object
    type(TTaggedWriter), intent(inout) :: taggedWriter

    !> should regression test data be written
    logical, intent(in) :: tWriteAutoTest

    !> File name for regression data
    character(*), intent(in) :: autoTestTagFile

    !> should machine readable output data be written
    logical, intent(in) :: tWriteTaggedOut

    !> File name for machine readable results data
    character(*), intent(in) :: taggedResultsFile

    !> should detailed.out be written to
    logical, intent(in) :: tWriteDetailedOut

    !> File id for detailed.out
    integer, intent(in) :: fdDetailedOut

    !> k-points
    real(dp), intent(in) :: kPoint(:,:)

    !> Weights for k-points
    real(dp), intent(in) :: kWeight(:)

    !> Vectors (in units of the lattice constants) to cells of the lattice
    real(dp), intent(in) :: cellVec(:,:)

    !> Index for which unit cell atoms are associated with
    integer, intent(in) :: iCellVec(:)

    !> Is this a periodic geometry
    logical, intent(in) :: tPeriodic

    !> Driving frequencies for dynamic polarization calculations (values of 0 or not allocated being
    !> the static case)
    real(dp), intent(in), allocatable :: omegas(:)

    integer :: iS, iK, iKS, iAt, iNeigh, iCart, iSCC, iLev, iSh, iSp, jAt, jAtf

    integer :: nSpin, nKpts, nOrbs, nIndepHam

    ! maximum allowed number of electrons in a single particle state
    real(dp) :: maxFill

    integer, allocatable :: nFilled(:,:), nEmpty(:,:)

    integer :: ii, jj, iGlob, jGlob, iEmpty, iFilled
    integer :: iSCCIter
    logical :: tStopSCC

    ! matrices for derivatives of terms in hamiltonian and outputs
    real(dp), allocatable :: dHam(:,:), idHam(:,:)
    real(dp) :: drho(size(over),size(ham, dim=2))
    real(dp) :: drhoExtra(size(over),size(ham, dim=2))
    real(dp), allocatable :: idRho(:,:), idRhoExtra(:,:)
    real(dp) :: dqIn(orb%mOrb,nAtom,size(ham, dim=2))
    real(dp) :: dqOut(orb%mOrb, nAtom, size(ham, dim=2)), dqUpDown(orb%mOrb, size(ham, dim=2))
    real(dp) :: dqInpRed(nMixElements), dqOutRed(nMixElements)
    real(dp) :: dqDiffRed(nMixElements), sccErrorQ
    real(dp) :: dqPerShell(orb%mShell,nAtom,size(ham, dim=2))

    real(dp), allocatable :: dqBlockIn(:,:,:,:), SSqrReal(:,:)
    real(dp), allocatable :: dqBlockOut(:,:,:,:)
    real(dp), allocatable :: dummy(:,:,:,:)

    ! derivative of potentials
    type(TPotentials) :: dPotential

    real(dp), allocatable :: shellPot(:,:,:), atomPot(:,:)

    real(dp) :: nF(size(Ef)), dEf(size(Ef))

    logical :: tSccCalc, tMetallic, tConverged

    real(dp), allocatable :: polarisability(:,:,:)

    real(dp), allocatable :: dEi(:,:,:,:)
    real(dp), allocatable :: dPsiReal(:,:,:,:)
    complex(dp), allocatable :: dPsiCmplx(:,:,:,:,:)

    integer :: fdResults, fdResponses

    ! dynamical perturbation
    integer :: iOmega, nOmega
    real(dp) :: omega

    ! used for range separated contributions, note this stays in the up/down representation
    ! throughout if spin polarised
    real(dp), pointer :: dRhoOutSqr(:,:,:), dRhoInSqr(:,:,:)
    real(dp), allocatable, target :: dRhoOut(:), dRhoIn(:)

  #:if WITH_SCALAPACK
    ! need distributed matrix descriptors
    integer :: desc(DLEN_), nn

    nn = denseDesc%fullSize
    call scalafx_getdescriptor(env%blacs%orbitalGrid, nn, nn, env%blacs%rowBlockSize,&
        & env%blacs%columnBlockSize, desc)
  #:endif

    if (allocated(omegas)) then
      nOmega = size(omegas)
    else
      nOmega = 1
    end if

    allocate(polarisability(3,3,nOmega))

    if (tPeriodic) then
      call error("Electric field polarizability not currently implemented for periodic systems")
    end if

    if (tFixEf) then
      call error("Perturbation expressions not currently implemented for fixed Fermi energy")
    end if

    write(stdOut,*)
    write(stdOut,*)'Perturbation calculation of electric polarisability'
    write(stdOut,*)

    nSpin = size(ham, dim=2)
    select case(nSpin)
    case(1,4)
      nIndepHam = 1
    case(2)
      nIndepHam = 2
    end select
    select case(nSpin)
    case(1)
      maxFill = 2.0_dp
    case(2,4)
      maxFill = 1.0_dp
    end select

    nOrbs = size(filling,dim=1)

    nKpts = size(filling,dim=2)

    allocate(dHam(size(ham,dim=1),nSpin))

    if (allocated(rangeSep)) then
      allocate(SSqrReal(nOrbs, nOrbs))
      SSqrReal(:,:) = 0.0_dp
      call unpackHS(SSqrReal, over, neighbourList%iNeighbour, nNeighbourSK, denseDesc%iAtomStart,&
          & iSparseStart, img2CentCell)
      allocate(dRhoOut(nOrbs * nOrbs * nSpin))
      dRhoOutSqr(1:nOrbs, 1:nOrbs, 1:nSpin) => dRhoOut(:nOrbs*nOrbs*nSpin)
      allocate(dRhoIn(nOrbs * nOrbs * nSpin))
      dRhoInSqr(1:nOrbs, 1:nOrbs, 1:nSpin) => dRhoIn(:nOrbs*nOrbs*nSpin)
    else
      dRhoInSqr => null()
      dRhoOutSqr => null()
    end if

    tSccCalc = (maxSccIter > 1)

    if (tDFTBU .or. allocated(onsMEs)) then
      allocate(dqBlockIn(orb%mOrb,orb%mOrb,nAtom,nSpin))
      allocate(dqBlockOut(orb%mOrb,orb%mOrb,nAtom,nSpin))
    end if

    call init(dPotential,orb,nAtom,nSpin)

    allocate(nFilled(nIndepHam, nKpts))
    allocate(nEmpty(nIndepHam, nKpts))

    if (allocated(spinW) .or. allocated(thirdOrd)) then
      allocate(shellPot(orb%mShell, nAtom, nSpin))
    end if
    if (allocated(thirdOrd)) then
      allocate(atomPot(nAtom, nSpin))
    end if


    ! If derivatives of eigenvalues are needed
    !allocate(dEi(nOrbs, nKpts, nIndepHam, 3))
    !dEi(:,:,:,:) = 0.0_dp

    nFilled(:,:) = -1
    do iS = 1, nIndepHam
      do iK = 1, nKPts
        do iLev = 1, nOrbs
          if ( filling(iLev,iK,iS) < epsilon(1.0) ) then
            nFilled(iS,iK) = iLev - 1
            exit
          end if
        end do
        if (nFilled(iS, iK) < 0) then
          nFilled(iS, iK) = nOrbs
        end if
      end do
    end do
    nEmpty(:,:) = -1
    do iS = 1, nIndepHam
      do iK = 1, nKpts
        do iLev = 1, nOrbs
          if ( abs( filling(iLev,iK,iS) - maxFill ) > epsilon(1.0)) then
            nEmpty(iS, iK) = iLev
            exit
          end if
        end do
        if (nEmpty(iS, iK) < 0) then
          nEmpty(iS, iK) = 1
        end if
      end do
    end do

    tMetallic = (.not.all(nFilled == nEmpty -1))

    ! if derivatives of valence wavefunctions needed
    !if (.not.tMetallic) then
    !  if (allocated(eigVecsReal)) then
    !    allocate(dPsiReal(size(eigVecsReal,dim=1), size(eigVecsReal,dim=2), nIndepHam, 3))
    !  else
    !    allocate(dPsiCmplx(size(eigvecsCplx,dim=1), size(eigvecsCplx,dim=2), nKpts, nIndepHam, 3))
    !  end if
    !end if

    if (tMetallic) then
      write(stdOut,*)'Metallic system'
    else
      write(stdOut,*)'Non-metallic system'
    end if

    if (tMetallic) then
      ! Density of electrons at the Fermi energy, required to correct later for shift in Fermi level
      ! at q=0 in metals
      do iS = 1, nIndepHam
        nf(iS) = 0.0_dp
        do iK = 1, nKpts
          do ii = nEmpty(iS, iK), nFilled(iS, iK)
            nf(iS) = nf(iS) + kWeight(iK) * deltamn(Ef(iS), eigvals(ii,iK,iS), tempElec)
          end do
        end do
        nf(iS) = maxFill * nf(iS)
      end do
      write(stdOut,*)'Density of states at the Fermi energy Nf:',nF
    end if

    dqOut(:,:,:) = 0.0_dp

    ! polarisation direction
    ! note, could MPI parallelise over this
    lpCart: do iCart = 1, 3

      dqIn(:,:,:) = 0.0_dp
      if (tDFTBU .or. allocated(onsMEs)) then
        dqBlockIn(:,:,:,:) = 0.0_dp
        dqBlockOut(:,:,:,:) = 0.0_dp
      end if

      dPotential%extAtom(:,:) = 0.0_dp
      dPotential%extShell(:,:,:) = 0.0_dp
      dPotential%extBlock(:,:,:,:) = 0.0_dp

      ! derivative wrt to electric field as a perturbation
      do iAt = 1, nAtom
        dPotential%extAtom(iAt,1) = coord(iCart,iAt)
      end do
      call total_shift(dPotential%extShell, dPotential%extAtom, orb, species)
      call total_shift(dPotential%extBlock, dPotential%extShell, orb, species)

      lpOmega: do iOmega = 1, nOmega

        if (allocated(omegas)) then
          omega = omegas(iOmega)
        else
          omega = 0.0_dp
        end if

        if (tSccCalc) then
          call reset(pChrgMixer, nMixElements)
          if (iOmega == 1) then
            dqInpRed(:) = 0.0_dp
            dqPerShell(:,:,:) = 0.0_dp
            if (allocated(rangeSep)) then
              dRhoIn(:) = 0.0_dp
              dRhoOut(:) = 0.0_dp
            end if
          end if
        end if

        if (tSccCalc) then
          write(stdOut,"(1X,A,T12,A)")'SCC Iter','Error'
        end if

        iSCCIter = 1
        tStopSCC = .false.
        lpSCC: do while (iSCCiter <= maxSccIter)

          dPotential%intAtom(:,:) = 0.0_dp
          dPotential%intShell(:,:,:) = 0.0_dp
          dPotential%intBlock(:,:,:,:) = 0.0_dp

          if (tDFTBU .or. allocated(onsMEs)) then
            dPotential%orbitalBlock(:,:,:,:) = 0.0_dp
          end if

          if (tSccCalc .and. iSCCiter>1) then
            call sccCalc%updateCharges(env, dqIn, orb, species)
            call sccCalc%updateShifts(env, orb, species, neighbourList%iNeighbour, img2CentCell)
            call sccCalc%getShiftPerAtom(dPotential%intAtom(:,1))
            call sccCalc%getShiftPerL(dPotential%intShell(:,:,1))

            if (allocated(spinW)) then
              call getChargePerShell(dqIn, orb, species, dqPerShell)
              call getSpinShift(shellPot, dqPerShell, species, orb, spinW)
              dPotential%intShell(:,:,:) = dPotential%intShell + shellPot
            end if

            if (allocated(thirdOrd)) then
              atomPot(:,:) = 0.0_dp
              shellPot(:,:,:) = 0.0_dp
              call thirdOrd%getdShiftdQ(atomPot(:,1), shellPot(:,:,1), species, neighbourList,&
                  & dqIn, img2CentCell, orb)
              dPotential%intAtom(:,1) = dPotential%intAtom(:,1) + atomPot(:,1)
              dPotential%intShell(:,:,1) = dPotential%intShell(:,:,1)&
                  & + shellPot(:,:,1)
            end if

            if (tDFTBU) then
              ! note the derivatives of both FLL and pSIC are the same (pSIC, i.e. case 2 in module)
              call getDftbUShift(dPotential%orbitalBlock, dqBlockIn, species, orb, 2,&
                  & UJ, nUJ, niUJ, iUJ)
            end if
            if (allocated(onsMEs)) then
              ! onsite corrections
              call addOnsShift(dPotential%orbitalBlock, dPotential%iOrbitalBlock,&
                  & dqBlockIn, dummy, onsMEs, species, orb)
            end if

          end if

          call total_shift(dPotential%intShell,dPotential%intAtom, orb,species)
          call total_shift(dPotential%intBlock,dPotential%intShell, orb,species)
          if (tDFTBU .or. allocated(onsMEs)) then
            dPotential%intBlock(:,:,:,:) = dPotential%intBlock&
                & + dPotential%orbitalBlock
          end if
          dPotential%intBlock(:,:,:,:) = dPotential%intBlock&
              & + dPotential%extBlock

          dHam(:,:) = 0.0_dp
          call add_shift(dHam, over, nNeighbourSK, neighbourList%iNeighbour, species, orb,&
              & iSparseStart, nAtom, img2CentCell, dPotential%intBlock)

          if (nSpin > 1) then
            dHam(:,:) = 2.0_dp * dHam(:,:)
            if (allocated(idHam)) then
              idHam(:,:) = 2.0_dp * idHam(:,:)
            end if
          end if
          call qm2ud(dHam)
          if (allocated(idHam)) then
            call qm2ud(idHam)
          end if

          dRho(:,:) = 0.0_dp
          if (allocated(idRho)) then
            idRho(:,:) = 0.0_dp
          end if

          ! evaluate derivative of density matrix
          if (allocated(eigVecsReal)) then

            do iKS = 1, parallelKS%nLocalKS

              iS = parallelKS%localKS(2, iKS)

              if (allocated(dRhoOut)) then
                ! replace with case that will get updated in dRhoStaticReal
                dRhoOutSqr(:,:,iS) = dRhoInSqr(:,:,iS)
              end if

              call dRhoStaticReal(env, dHam, omega, neighbourList, nNeighbourSK, iSparseStart,&
                  & img2CentCell, denseDesc, iKS, parallelKS, nFilled(:,1), nEmpty(:,1),&
                  & eigVecsReal, eigVals, Ef, tempElec, orb, drho(:,iS), iCart, dRhoOutSqr,&
                  & rangeSep, over, nNeighbourLC,&
                #:if WITH_SCALAPACK
                  & desc,&
                #:endif
                  & dEi, dPsiReal)
            end do

          elseif (nSpin > 2) then

            do iKS = 1, parallelKS%nLocalKS

              iK = parallelKS%localKS(1, iKS)

              call dRhoStaticPauli(env, dHam, idHam, omega, neighbourList, nNeighbourSK,&
                  & iSparseStart, img2CentCell, denseDesc, parallelKS, nFilled(:, iK),&
                  & nEmpty(:, iK), eigvecsCplx, eigVals, Ef, tempElec, orb, dRho, idRho, kPoint,&
                  & kWeight, iCellVec, cellVec, iKS, iCart,&
                #:if WITH_SCALAPACK
                  & desc,&
                #:endif
                  & dEi, dPsiCmplx)

            end do

          else

            call error("Shouldn't be here")

          end if

        #:if WITH_SCALAPACK
          ! Add up and distribute density matrix contributions from each group
          call mpifx_allreduceip(env%mpi%globalComm, dRho, MPI_SUM)
        #:endif

          dRhoExtra = 0.0_dp
          if (tMetallic) then
            ! correct for Fermi level shift for q=0 fields

            do iKS = 1, parallelKS%nLocalKS
              iK = parallelKS%localKS(1, iKS)
              iS = parallelKS%localKS(2, iKS)

              dqOut(:,:,iS) = 0.0_dp
              call mulliken(dqOut(:,:,iS), over, drho(:,iS), orb, &
                  & neighbourList%iNeighbour, nNeighbourSK, img2CentCell, iSparseStart)

              dEf(iS) = -sum(dqOut(:, :, iS)) / nF(iS)

              if (abs(dEf(iS)) > 10.0_dp*epsilon(1.0_dp)) then
                ! Fermi level changes, so need to correct for the change in the number of charges

                if (allocated(eigVecsReal)) then

                  ! real case, no k-points
                  call dRhoFermiChangeStaticReal(dRhoExtra(:, iS), env, parallelKS, iKS,&
                      & neighbourList, nNeighbourSK, img2CentCell, iSparseStart, dEf, Ef,&
                      & nFilled(:,iK), nEmpty(:,iK), eigVecsReal, orb, denseDesc, tempElec,&
                      & eigVals, dRhoOutSqr&
                    #:if WITH_SCALAPACK
                      &, desc&
                    #:endif
                      &)

                elseif (nSpin > 2) then

                  ! two component wavefunction cases
                  call dRhoFermiChangeStaticPauli(dRhoExtra, idRhoExtra, env, parallelKS, iKS,&
                      & kPoint, kWeight, iCellVec, cellVec, neighbourList, nNEighbourSK,&
                      & img2CentCell, iSparseStart, dEf, Ef, nFilled(:,iK), nEmpty(:,iK),&
                      & eigVecsCplx, orb, denseDesc, tempElec, eigVals&
                    #:if WITH_SCALAPACK
                      &, desc&
                    #:endif
                      &)

                else

                  ! k-points
                  call error("Not added yet")

                end if

              end if

            end do

          #:if WITH_SCALAPACK
            ! Add up and distribute density matrix contribution from each group
            call mpifx_allreduceip(env%mpi%globalComm, dRhoExtra, MPI_SUM)
          #:endif
            dRho(:,:) = dRho + dRhoExtra

          end if

          dRho(:,:) = maxFill * drho
          if (allocated(dRhoOut)) then
            dRhoOut(:) = maxFill * dRhoOut
          end if
          call ud2qm(dRho)

          if (allocated(idRho)) then
            idRho(:,:) = maxFill * drho
            call ud2qm(idRho)
          end if

          dqOut(:,:,:) = 0.0_dp
          do iS = 1, nSpin
            call mulliken(dqOut(:,:,iS), over, drho(:,iS), orb, &
                & neighbourList%iNeighbour, nNeighbourSK, img2CentCell, iSparseStart)
            if (tDFTBU .or. allocated(onsMEs)) then
              dqBlockOut(:,:,:,iS) = 0.0_dp
              call mulliken(dqBlockOut(:,:,:,iS), over, drho(:,iS), orb, neighbourList%iNeighbour,&
                  & nNeighbourSK, img2CentCell, iSparseStart)
            end if
          end do

          if (tSccCalc) then

            if (allocated(rangeSep)) then
              dqDiffRed(:) = dRhoOut - dRhoIn
            else
              dqOutRed = 0.0_dp
              call OrbitalEquiv_reduce(dqOut, iEqOrbitals, orb,&
                  & dqOutRed(:nIneqMixElements))
              if (tDFTBU) then
                call AppendBlock_reduce(dqBlockOut, iEqBlockDFTBU, orb, dqOutRed )
              end if
              if (allocated(onsMEs)) then
                call onsBlock_reduce(dqBlockOut, iEqBlockOnSite, orb, dqOutRed)
              end if
              dqDiffRed(:) = dqOutRed - dqInpRed
            end if
            sccErrorQ = maxval(abs(dqDiffRed))

            write(stdOut,"(1X,I0,T10,E20.12)")iSCCIter, sccErrorQ
            tConverged = (sccErrorQ < sccTol)

            if ((.not. tConverged) .and. iSCCiter /= maxSccIter) then
              if (iSCCIter == 1) then
                if (allocated(rangeSep)) then
                  dRhoIn(:) = dRhoOut
                  call denseMulliken(dRhoInSqr, SSqrReal, denseDesc%iAtomStart, dqIn)
                else
                  dqIn(:,:,:) = dqOut
                  dqInpRed(:) = dqOutRed(:)
                  if (tDFTBU .or. allocated(onsMEs)) then
                    dqBlockIn(:,:,:,:) = dqBlockOut(:,:,:,:)
                  end if
                end if

              else

                if (allocated(rangeSep)) then
                  call mix(pChrgMixer, dRhoIn, dqDiffRed)
                  call denseMulliken(dRhoInSqr, SSqrReal, denseDesc%iAtomStart, dqIn)
                else
                  call mix(pChrgMixer, dqInpRed, dqDiffRed)
                #:if WITH_MPI
                  ! Synchronise charges in order to avoid mixers that store a history drifting apart
                  call mpifx_allreduceip(env%mpi%globalComm, dqInpRed, MPI_SUM)
                  dqInpRed(:) = dqInpRed / env%mpi%globalComm%size
                #:endif

                  call OrbitalEquiv_expand(dqInpRed(:nIneqMixElements), iEqOrbitals, orb, dqIn)

                  if (tDFTBU .or. allocated(onsMEs)) then
                    dqBlockIn(:,:,:,:) = 0.0_dp
                    if (tDFTBU) then
                      call Block_expand( dqInpRed ,iEqBlockDFTBU, orb, dqBlockIn, species(:nAtom),&
                          & nUJ, niUJ, iUJ, orbEquiv=iEqOrbitals )
                    else
                      call Onsblock_expand(dqInpRed, iEqBlockOnSite, orb, dqBlockIn,&
                          & orbEquiv=iEqOrbitals)
                    end if
                  end if
                end if

              end if

              if (allocated(rangeSep)) then
                call ud2qm(dqIn)
              end if

            end if
          else
            tConverged = .true.
          end if

          if (tConverged) then
            exit lpSCC
          end if

          if (allocated(spinW)) then
            dqPerShell = 0.0_dp
            do iAt = 1, nAtom
              iSp = species(iAt)
              do iSh = 1, orb%nShell(iSp)
                dqPerShell(iSh,iAt,:nSpin) = dqPerShell(iSh,iAt,:nSpin) +&
                    & sum(dqIn(orb%posShell(iSh,iSp): orb%posShell(iSh+1,iSp)-1,iAt,:nSpin),dim=1)
              end do
            end do

          end if

          iSCCIter = iSCCIter +1

        end do lpSCC

        if (tMetallic) then
          write(stdOut,*)
          write(stdOut,"(A,2E20.12)")'d E_f / d E_'//trim(direction(iCart))//':', dEf
          write(stdOut,*)
        end if

        do ii = 1, 3
          polarisability(ii, iCart, iOmega) =&
              & -sum(sum(dqOut(:,:nAtom,1),dim=1)*coord(ii,:nAtom))
        end do

      end do lpOmega

    end do lpCart

  #:if WITH_SCALAPACK
    if (allocated(dEi)) then
      call mpifx_allreduceip(env%mpi%globalComm, dEi, MPI_SUM)
    end if
  #:endif

    if (allocated(omegas)) then
      if (omegas(1) == 0.0_dp) then
        write(stdOut,*)'Static polarisability'
        do iCart = 1, 3
          write(stdOut,"(3E20.12)")polarisability(:, iCart, 1)
        end do
        write(stdOut,*)
      end if
      if (tIoProc) then
        open(newunit=fdResponses, file="polarisabilities.dat", action="write", status="replace")
        do iOmega = 1, nOmega
          write(fdResponses,"(10E20.12)")omegas(iOmega)*Hartree__eV, polarisability(:, :, iOmega)
        end do
        close(fdResponses)
      end if
    else
      write(stdOut,*)'Static polarisability'
      do iCart = 1, 3
        write(stdOut,"(3E20.12)")polarisability(:, iCart, 1)
      end do
      write(stdOut,*)
    end if

    if (tWriteAutoTest) then
      open(newunit=fdResults, file=trim(autoTestTagFile), position="append")
      call taggedWriter%write(fdResults, tagLabels%dmudEPerturb, polarisability(:,:,nOmega))
      close(fdResults)
    end if
    if (tWriteTaggedOut) then
      open(newunit=fdResults, file=trim(taggedResultsFile), action="write", status="old",&
          & position="append")
      call taggedWriter%write(fdResults, tagLabels%dmudEPerturb, polarisability)
      call taggedWriter%write(fdResults, tagLabels%dqdEPerturb, sum(dqout(:,:natom,:), dim = 1))
      close(fdResults)
    end if

    if (tWriteDetailedOut) then
      write(fdDetailedOut,*)'Linear response to electric field'
      if (allocated(dEi)) then

        write(fdDetailedOut,*)"Warning! Gauge dependent property:"
        do iCart = 1, 3
          write(fdDetailedOut,"(A)")'d e_i / d E_'//trim(direction(iCart))//':'
          do iS = 1, nIndepHam
            do iK = 1, nKpts
              do iLev = 1, nOrbs
                write(fdDetailedOut,*)iS,iLev,dEi(iLev,iK,iS,iCart)
              end do
            end do
          end do
          write(fdDetailedOut,*)
        end do
      end if

      write(fdDetailedOut,*)'Atomic charge derivatives (a.u.)'
      do iCart = 1, 3
        write(fdDetailedOut,"(A)")'d q / d E_'//trim(direction(iCart))//':'
        do iAt = 1, nAtom
          dqUpDown(:orb%nOrbAtom(iAt), :) = dqOut(:orb%nOrbAtom(iAt),iAt,:)
          call qm2ud(dqUpDown)
          do iS = 1, nSpin
            write(fdDetailedOut,"(I4,4E20.12)")iAt, sum(dqUpDown(:orb%nOrbAtom(iAt), :))
          end do
        end do
      end do

      write(fdDetailedOut,*)
      write(fdDetailedOut,*)'Polarisability (a.u.)'
      do iOmega = 1, nOmega
        do iCart = 1, 3
          write(fdDetailedOut,"(3E20.12)")polarisability(:, iCart, iOmega)
        end do
        if (nOmega > 1) then
          write(fdDetailedOut,*)
        end if
      end do
    end if
    write(fdDetailedOut,*)


  end subroutine staticPerturWrtE


  !> Calculate the derivative of density matrix from derivative of hamiltonian in static case at
  !> q=0, k=0
  subroutine dRhoStaticReal(env, dHam, omega, neighbourList, nNeighbourSK, iSparseStart,&
      & img2CentCell, denseDesc, iKS, parallelKS, nFilled, nEmpty, eigVecsReal, eigVals, Ef,&
      & tempElec, orb, dRhoSparse, iCart, dRhoSqr, rangeSep, over, nNeighbourLC,&
    #:if WITH_SCALAPACK
      & desc,&
    #:endif
      & dEi, dPsi)

    !> Environment settings
    type(TEnvironment), intent(inout) :: env

    !> Derivative of the hamiltonian
    real(dp), intent(in) :: dHam(:,:)

    !> Driving frequency (if non-zero, is required for time dependent DFTB)
    real(dp), intent(out) :: omega

    !> list of neighbours for each atom
    type(TNeighbourList), intent(in) :: neighbourList

    !> Number of neighbours for each of the atoms
    integer, intent(in) :: nNeighbourSK(:)

    !> Index array for the start of atomic blocks in sparse arrays
    integer, intent(in) :: iSparseStart(:,:)

    !> map from image atoms to the original unique atom
    integer, intent(in) :: img2CentCell(:)

    !> Dense matrix descriptor
    type(TDenseDescr), intent(in) :: denseDesc

    !> Particular spin/k-point
    integer, intent(in) :: iKS

    !> K-points and spins to process
    type(TParallelKS), intent(in) :: parallelKS

    !> ground state eigenvectors
    real(dp), intent(in) :: eigVecsReal(:,:,:)

    !> Eigenvalue of each level, kpoint and spin channel
    real(dp), intent(in) :: eigvals(:,:,:)

    !> Fermi level(s)
    real(dp), intent(in) :: Ef(:)

    !> Last (partly) filled level in each spin channel
    integer, intent(in) :: nFilled(:)

    !> First (partly) empty level in each spin channel
    integer, intent(in) :: nEmpty(:)

    !> Electron temperature
    real(dp), intent(in) :: tempElec

    !> Atomic orbital information
    type(TOrbitals), intent(in) :: orb

    !> returning dRhoSparse on exit
    real(dp), intent(out) :: dRhoSparse(:)

    !> Cartesian direction of perturbation
    integer, intent(in) :: iCart

    !> Derivative of rho as a square matrix, if needed
    real(dp), pointer :: dRhoSqr(:,:,:)

    !> Data for range-separated calculation
    type(RangeSepFunc), allocatable, intent(inout) :: rangeSep

    !> sparse overlap matrix
    real(dp), intent(in) :: over(:)

    !> Number of neighbours for each of the atoms for the exchange contributions in the long range
    !> functional
    integer, intent(inout), allocatable :: nNeighbourLC(:)

  #:if WITH_SCALAPACK
    !> BLACS matrix descriptor
    integer, intent(in) :: desc(DLEN_)
  #:endif

    !> Derivative of single particle eigenvalues
    real(dp), allocatable, intent(inout) :: dEi(:,:,:,:)

    !> Optional derivatives of single particle wavefunctions
    real(dp), allocatable, intent(inout) :: dPsi(:,:,:,:)

    integer :: ii, jj, iGlob, jGlob, iFilled, iEmpty, iS, iK, nOrb
    real(dp) :: workLocal(size(eigVecsReal,dim=1), size(eigVecsReal,dim=2))
    complex(dp), allocatable :: cWorkLocal(:, :)
    real(dp), allocatable :: dRho(:,:)

    integer :: iSignOmega

    iK = parallelKS%localKS(1, iKS)
    iS = parallelKS%localKS(2, iKS)

    if (omega /= 0.0_dp) then
      allocate(cWorkLocal(size(eigVecsReal,dim=1), size(eigVecsReal,dim=2)))
    end if

    if (allocated(dEi)) then
      dEi(:, iK, iS, iCart) = 0.0_dp
    end if
    if (allocated(dPsi)) then
      dPsi(:, :, iS, iCart) = 0.0_dp
    end if

    workLocal(:,:) = 0.0_dp
    allocate(dRho(size(eigVecsReal,dim=1), size(eigVecsReal,dim=2)))
    dRho(:,:) = 0.0_dp

  #:if WITH_SCALAPACK

    ! dH in square form
    call unpackHSRealBlacs(env%blacs, dHam(:,iS), neighbourList%iNeighbour, nNeighbourSK,&
        & iSparseStart, img2CentCell, denseDesc, workLocal)

    ! dH times c_i
    call pblasfx_psymm(workLocal, denseDesc%blacsOrbSqr, eigVecsReal(:,:,iKS),&
        & denseDesc%blacsOrbSqr, dRho, denseDesc%blacsOrbSqr) !, mm=nFilled(iS))

    ! c_i times dH times c_i
    call pblasfx_pgemm(eigVecsReal(:,:,iKS), denseDesc%blacsOrbSqr, dRho,&
        & denseDesc%blacsOrbSqr, workLocal, denseDesc%blacsOrbSqr, transa="T")

    ! derivative of eigenvalues stored diagonal of matrix workLocal, from <c|h'|c>
    if (allocated(dEi)) then
      if (iGlob == jGlob) then
        !if (iGlob == jGlob) then workLocal(ii,jj) contains a derivative of an eigenvalue
        dEi(iGlob, iK, iS, iCart) = workLocal(ii,jj)
      end if
    end if


    ! weight with inverse of energy differences
    do jj = 1, size(workLocal,dim=2)
      jGlob = scalafx_indxl2g(jj, desc(NB_), env%blacs%orbitalGrid%mycol, desc(CSRC_),&
          & env%blacs%orbitalGrid%ncol)
      do ii = 1, size(workLocal,dim=1)
        iGlob = scalafx_indxl2g(ii, desc(MB_), env%blacs%orbitalGrid%myrow, desc(RSRC_),&
            & env%blacs%orbitalGrid%nrow)

        workLocal(ii,jj) = workLocal(ii,jj) * &
            & invDiff(eigvals(jGlob,1,iS),eigvals(iGlob,1,iS),Ef(iS),tempElec)&
            & * theta(eigvals(jGlob,1,iS),eigvals(iGlob,1,iS),tempElec)

      end do
    end do

    ! Derivatives of states
    call pblasfx_pgemm(eigVecsReal(:,:,iKS), denseDesc%blacsOrbSqr, workLocal,&
        & denseDesc%blacsOrbSqr, dRho, denseDesc%blacsOrbSqr)

    if (allocated(dPsi)) then
      dPsi(:, :, iS, iCart) = workLocal
    end if

    ! Form derivative of occupied density matrix
    call pblasfx_pgemm(dRho, denseDesc%blacsOrbSqr,eigVecsReal(:,:,iKS),&
        & denseDesc%blacsOrbSqr, workLocal, denseDesc%blacsOrbSqr, transb="T",&
        & kk=nFilled(iS))

    dRho(:,:) = workLocal
    ! and symmetrize
    call pblasfx_ptran(workLocal, denseDesc%blacsOrbSqr, dRho, denseDesc%blacsOrbSqr,&
        & beta=1.0_dp)

  #:else

    ! serial case
    nOrb = size(dRho, dim = 1)

    ! dH matrix in square form
    call unpackHS(dRho, dHam(:,iS), neighbourList%iNeighbour, nNeighbourSK,&
        & denseDesc%iAtomStart, iSparseStart, img2CentCell)

    if (allocated(rangeSep)) then
      call unpackHS(workLocal, over, neighbourList%iNeighbour, nNeighbourSK,&
          & denseDesc%iAtomStart, iSparseStart, img2CentCell)
      call rangeSep%addLRHamiltonian(env, dRhoSqr(:,:,iS), over, neighbourList%iNeighbour,&
          & nNeighbourLC, denseDesc%iAtomStart, iSparseStart, orb, dRho, workLocal)
    end if

    ! form |c> H' <c|
    call symm(workLocal, 'l', dRho, eigVecsReal(:,:,iS))
    workLocal(:,:) = matmul(transpose(eigVecsReal(:,:,iS)), workLocal)

    ! diagonal elements of workLocal are now derivatives of eigenvalues if needed
    if (allocated(dEi)) then
      do ii = 1, nOrb
        dEi(ii, iK, iS, iCart) = workLocal(ii,ii)
      end do
    end if

    if (omega == 0.0_dp) then
      ! static case

      ! Form actual perturbation U matrix for eigenvectors (potentially at finite T) by weighting
      ! the elements
      do iFilled = 1, nFilled(iS)
        do iEmpty = nEmpty(iS), nOrb
          workLocal(iEmpty, iFilled) = workLocal(iEmpty, iFilled) * &
              & invDiff(eigvals(iFilled, iK, iS), eigvals(iEmpty, iK, iS), Ef(iS), tempElec)&
              & *theta(eigvals(iFilled, iK, iS), eigvals(iEmpty, iK, iS), tempElec)
        end do
      end do

      ! calculate the derivatives of the eigenvectors
      workLocal(:, :nFilled(iS)) =&
          & matmul(eigVecsReal(:, nEmpty(iS):, iS), workLocal(nEmpty(iS):, :nFilled(iS)))

      if (allocated(dPsi)) then
        dPsi(:, :, iS, iCart) = workLocal
      end if

      ! zero the uncalculated virtual states
      workLocal(:, nFilled(iS)+1:) = 0.0_dp

      ! form the derivative of the density matrix
      dRho(:,:) = matmul(workLocal(:, :nFilled(iS)), transpose(eigVecsReal(:, :nFilled(iS), iS)))&
          & + matmul(eigVecsReal(:, :nFilled(iS), iS), transpose(workLocal(:, :nFilled(iS))))

    else
      ! dynamic case

      dRho(:,:) = 0.0_dp
      ! loop over real and imaginary frequencies of driving
      do iSignOmega = -1, 1, 2

        ! Form actual perturbation U matrix for eigenvectors (potentially at finite T) by weighting
        ! the elements
        do iFilled = 1, nFilled(iS)
          do iEmpty = nEmpty(iS), nOrb
            cWorkLocal(iEmpty, iFilled) = workLocal(iEmpty, iFilled) * &
                & ( theta(eigvals(iEmpty, iK, iS), Ef(iS), tempElec)&
                & - theta(eigvals(iFilled, iK, iS), Ef(iS), tempElec) ) &
                & * theta(eigvals(iFilled, iK, iS), eigvals(iEmpty, iK, iS), tempElec)&
                & / ( eigvals(iEmpty, iK, iS) - eigvals(iFilled, iK, iS) + iSignOmega * omega + eta)
          end do
        end do

        ! calculate the derivatives of the eigenvectors
        cWorkLocal(:, :nFilled(iS)) =&
            & matmul(eigVecsReal(:, nEmpty(iS):, iS), cWorkLocal(nEmpty(iS):, :nFilled(iS)))

        ! form the derivative of the density matrix
        dRho(:,:) = dRho + 0.5_dp*real(&
            & matmul(cWorkLocal(:, :nFilled(iS)), transpose(eigVecsReal(:, :nFilled(iS), iS)))&
            & + matmul(eigVecsReal(:, :nFilled(iS), iS), transpose(cWorkLocal(:, :nFilled(iS)))) )

      end do

    end if

  #:endif

    dRhoSparse(:) = 0.0_dp
  #:if WITH_SCALAPACK
    call packRhoRealBlacs(env%blacs, denseDesc, dRho, neighbourList%iNeighbour, nNeighbourSK,&
        & orb%mOrb, iSparseStart, img2CentCell, dRhoSparse)
  #:else
    call packHS(dRhoSparse, dRho, neighbourList%iNeighbour, nNeighbourSK, orb%mOrb,&
        & denseDesc%iAtomStart, iSparseStart, img2CentCell)
  #:endif

    if (associated(dRhoSqr)) then
      dRhoSqr(:,:,iS) = dRho
    end if

  end subroutine dRhoStaticReal


  !> Calculate the change in the density matrix due to shift in the Fermi energy
  subroutine dRhoFermiChangeStaticReal(dRhoExtra, env, parallelKS, iKS, neighbourList,&
      & nNEighbourSK, img2CentCell, iSparseStart, dEf, Ef, nFilled, nEmpty, eigVecsReal, orb,&
      & denseDesc, tempElec, eigVals, dRhoSqr&
    #:if WITH_SCALAPACK
      &, desc&
    #:endif
      &)

    !> Additional contribution to the density matrix to cancel effect of Fermi energy change
    real(dp), intent(out) :: dRhoExtra(:)

    !> Environment settings
    type(TEnvironment), intent(inout) :: env

    !> K-points and spins to process
    type(TParallelKS), intent(in) :: parallelKS

    !> spin/kpoint channel
    integer, intent(in) :: iKS

    !> list of neighbours for each atom
    type(TNeighbourList), intent(in) :: neighbourList

    !> Number of neighbours for each of the atoms
    integer, intent(in) :: nNeighbourSK(:)

    !> map from image atoms to the original unique atom
    integer, intent(in) :: img2CentCell(:)

    !> Index array for the start of atomic blocks in sparse arrays
    integer, intent(in) :: iSparseStart(:,:)

    !> Fermi level derivative
    real(dp), intent(in) :: dEf(:)

    !> Fermi level
    real(dp), intent(in) :: Ef(:)

    !> Last (partly) filled level in each spin channel
    integer, intent(in) :: nFilled(:)

    !> First (partly) empty level in each spin channel
    integer, intent(in) :: nEmpty(:)

    !> ground state eigenvectors
    real(dp), intent(in) :: eigVecsReal(:,:,:)

    !> Atomic orbital information
    type(TOrbitals), intent(in) :: orb

    !> Dense matrix descriptor
    type(TDenseDescr), intent(in) :: denseDesc

    !> Electron temperature
    real(dp), intent(in) :: tempElec

    !> Eigenvalue of each level, kpoint and spin channel
    real(dp), intent(in) :: eigvals(:,:,:)

    !> Derivative of rho as a square matrix, if needed
    real(dp), pointer :: dRhoSqr(:,:,:)

  #:if WITH_SCALAPACK
    !> BLACS matrix descriptor
    integer, intent(in) :: desc(DLEN_)
  #:endif

    integer :: iFilled, jj, jGlob, nSpin, iS
    real(dp) :: workReal(size(eigVecsReal, dim=1), size(eigVecsReal, dim=2))

    iS = parallelKS%localKS(2, iKS)
    nSpin = size(eigVecsReal, dim=3)

    workReal(:,:) = 0.0_dp

  #:if WITH_SCALAPACK

    do jj = 1, size(workReal,dim=2)
      jGlob = scalafx_indxl2g(jj,desc(NB_), env%blacs%orbitalGrid%mycol, desc(CSRC_),&
          & env%blacs%orbitalGrid%ncol)
      if (jGlob >= nEmpty(iS) .and. jGlob <= nFilled(iS)) then
        workReal(:,jj) = eigVecsReal(:,jj,iKS) * &
            & deltamn(eigVals(jGlob, 1, iKS), Ef(iS), tempElec) * dEf(iS)
      end if
    end do
    call pblasfx_pgemm(workReal, denseDesc%blacsOrbSqr,eigVecsReal(:,:,iKS),&
        & denseDesc%blacsOrbSqr, workReal, denseDesc%blacsOrbSqr, transb="T",&
        & alpha=real(3-nSpin,dp))

  #:else

    do iFilled = nEmpty(iS), nFilled(iS)
      workReal(:, iFilled) = eigVecsReal(:, iFilled, iS) * &
          & deltamn(eigvals(iFilled, 1, iS), Ef(iS), tempElec) * dEf(iS)
    end do
    workReal(:, :) = real(3-nSpin,dp)&
        & * matmul(workReal(:, nEmpty(iS):nFilled(iS)),&
        & transpose(eigVecsReal(:, nEmpty(iS):nFilled(iS), iS)))

  #:endif

    ! pack extra term into density matrix
  #:if WITH_SCALAPACK
    call packRhoRealBlacs(env%blacs, denseDesc, workReal, neighbourList%iNeighbour, nNeighbourSK,&
        & orb%mOrb, iSparseStart, img2CentCell, drhoExtra)
  #:else
    call packHS(drhoExtra, workReal, neighbourList%iNeighbour, nNeighbourSK, orb%mOrb,&
        & denseDesc%iAtomStart, iSparseStart, img2CentCell)
  #:endif

    if  (associated(dRhoSqr)) then
      dRhoSqr(:,:,iS) = dRhoSqr(:,:,iS) + workReal
    end if

  end subroutine dRhoFermiChangeStaticReal


  !> Calculate the derivative of density matrix from derivative of hamiltonian in static case at q=0
  subroutine dRhoStaticPauli(env, dHam, idHam, omega, neighbourList, nNeighbourSK, iSparseStart,&
      & img2CentCell, denseDesc, parallelKS, nFilled, nEmpty, eigVecsCplx, eigVals, Ef, tempElec,&
      & orb, dRhoSparse, idRhoSparse, kPoint, kWeight, iCellVec, cellVec, iKS, iCart,&
    #:if WITH_SCALAPACK
      & desc,&
    #:endif
      & dEi, dPsi)

    !> Environment settings
    type(TEnvironment), intent(in) :: env

    !> Derivative of the hamiltonian
    real(dp), intent(in) :: dHam(:,:)

    !> Derivative of the imaginary part of the hamiltonian
    real(dp), intent(in), allocatable :: idHam(:,:)

    !> Driving frequency (if non-zero, is required for time dependent DFTB)
    real(dp), intent(out) :: omega

    !> list of neighbours for each atom
    type(TNeighbourList), intent(in) :: neighbourList

    !> Number of neighbours for each of the atoms
    integer, intent(in) :: nNeighbourSK(:)

    !> Index array for the start of atomic blocks in sparse arrays
    integer, intent(in) :: iSparseStart(:,:)

    !> map from image atoms to the original unique atom
    integer, intent(in) :: img2CentCell(:)

    !> Dense matrix descriptor
    type(TDenseDescr), intent(in) :: denseDesc

    !> K-points and spins to process
    type(TParallelKS), intent(in) :: parallelKS

    !> ground state eigenvectors
    complex(dp), intent(in) :: eigVecsCplx(:,:,:)

    !> Eigenvalue of each level, kpoint and spin channel
    real(dp), intent(in) :: eigvals(:,:,:)

    !> Fermi level(s)
    real(dp), intent(in) :: Ef(:)

    !> Last (partly) filled level in each spin channel
    integer, intent(in) :: nFilled(:)

    !> First (partly) empty level in each spin channel
    integer, intent(in) :: nEmpty(:)

    !> Electron temperature
    real(dp), intent(in) :: tempElec

    !> Atomic orbital information
    type(TOrbitals), intent(in) :: orb

    !> returning dRhoSparse on exit
    real(dp), intent(out) :: dRhoSparse(:,:)

    !> returning imaginary part of dRhoSparse on exit
    real(dp), intent(out), allocatable :: idRhoSparse(:,:)

    !> k-points
    real(dp), intent(in) :: kPoint(:,:)

    !> Weights for k-points
    real(dp), intent(in) :: kWeight(:)

    !> Vectors (in units of the lattice constants) to cells of the lattice
    real(dp), intent(in) :: cellVec(:,:)

    !> Index for which unit cell atoms are associated with
    integer, intent(in) :: iCellVec(:)

    !> spin/kpoint channel
    integer, intent(in) :: iKS

    !> Cartesian direction of perturbation
    integer, intent(in) :: iCart

  #:if WITH_SCALAPACK
    !> BLACS matrix descriptor
    integer, intent(in) :: desc(DLEN_)
  #:endif

    !> Derivative of single particle eigenvalues
    real(dp), allocatable, intent(inout) :: dEi(:,:,:,:)

    !> Optional derivatives of single particle wavefunctions
    complex(dp), allocatable, intent(inout) :: dPsi(:,:,:,:,:)

    integer :: ii, jj, iGlob, jGlob, iFilled, iEmpty, iK, iS, nOrb
    complex(dp) :: workLocal(size(eigVecsCplx,dim=1), size(eigVecsCplx,dim=2))
    complex(dp), allocatable :: cWorkLocal(:, :)
    complex(dp) :: dRho(size(eigVecsCplx,dim=1), size(eigVecsCplx,dim=2))

    integer :: iSignOmega

    iK = parallelKS%localKS(1, iKS)
    iS = parallelKS%localKS(2, iKS)

    if (omega /= 0.0_dp) then
      allocate(cWorkLocal(size(eigVecsCplx,dim=1), size(eigVecsCplx,dim=2)))
    end if

    if (allocated(dEi)) then
      dEi(:, iK, iS, iCart) = 0.0_dp
    end if
    if (allocated(dPsi)) then
      dPsi(:, :, iK, iS, iCart) = cmplx(0,0,dp)
    end if

    workLocal(:,:) = cmplx(0,0,dp)
    dRho(:,:) = cmplx(0,0,dp)

  #:if WITH_SCALAPACK

    ! dH in square form
    if (allocated(idHam)) then
      call unpackHPauliBlacs(env%blacs, dHam, kPoint(:,iK), neighbourList%iNeighbour,&
          & nNeighbourSK, iCellVec, cellVec, iSparseStart, img2CentCell, orb%mOrb, denseDesc,&
          & workLocal, iorig=idHam)
    else
      call unpackHPauliBlacs(env%blacs, dHam, kPoint(:,iK), neighbourList%iNeighbour,&
          & nNeighbourSK, iCellVec, cellVec, iSparseStart, img2CentCell, orb%mOrb, denseDesc,&
          & workLocal)
    end if

    ! dH times c_i
    call pblasfx_phemm(workLocal, denseDesc%blacsOrbSqr, eigVecsCplx(:,:,iKS),&
        & denseDesc%blacsOrbSqr, dRho, denseDesc%blacsOrbSqr) !, mm=nFilled(iS))

    ! c_i times dH times c_i
    call pblasfx_pgemm(eigVecsCplx(:,:,iKS), denseDesc%blacsOrbSqr, dRho,&
        & denseDesc%blacsOrbSqr, workLocal, denseDesc%blacsOrbSqr, transa="C")

    ! derivative of eigenvalues stored diagonal of matrix workLocal, from <c|h'|c>
    if (allocated(dEi)) then
      do jj = 1, size(workLocal,dim=2)
        jGlob = scalafx_indxl2g(jj, desc(NB_), env%blacs%orbitalGrid%mycol, desc(CSRC_),&
            & env%blacs%orbitalGrid%ncol)
        do ii = 1, size(workLocal,dim=1)
          iGlob = scalafx_indxl2g(ii, desc(MB_), env%blacs%orbitalGrid%myrow, desc(RSRC_),&
              & env%blacs%orbitalGrid%nrow)
          if (iGlob == jGlob) then
            dEi(iGlob, iK, iS, iCart) = real(workLocal(ii,jj),dp)
          end if
        end do
      end do
    end if

    ! weight matrix with inverse of energy differences
    do jj = 1, size(workLocal,dim=2)
      jGlob = scalafx_indxl2g(jj, desc(NB_), env%blacs%orbitalGrid%mycol, desc(CSRC_),&
          & env%blacs%orbitalGrid%ncol)
      if (jGlob > nFilled(1)) then
        workLocal(:, jj) = 0.0_dp
        cycle
      end if
      do ii = 1, size(workLocal,dim=1)
        iGlob = scalafx_indxl2g(ii, desc(MB_), env%blacs%orbitalGrid%myrow, desc(RSRC_),&
            & env%blacs%orbitalGrid%nrow)
        if (iGlob < nEmpty(1)) then
          workLocal(ii, :) = 0.0_dp
          cycle
        end if
        workLocal(ii, jj) = workLocal(ii, jj) * &
            & invDiff(eigvals(jGlob,iK,iS),eigvals(iGlob,iK,iS),Ef(iS),tempElec)&
            & * theta(eigvals(jGlob,iK,iS),eigvals(iGlob,iK,iS),tempElec)
      end do
    end do

    ! Derivatives of states
    call pblasfx_pgemm(eigVecsCplx(:,:,iKS), denseDesc%blacsOrbSqr, workLocal,&
        & denseDesc%blacsOrbSqr, dRho, denseDesc%blacsOrbSqr)

    if (allocated(dPsi)) then
      dPsi(:, :, iK, iS, iCart) = workLocal
    end if

    ! Form derivative of occupied density matrix
    call pblasfx_pgemm(dRho, denseDesc%blacsOrbSqr,eigVecsCplx(:,:,iKS),&
        & denseDesc%blacsOrbSqr, workLocal, denseDesc%blacsOrbSqr, transb="C",&
        & kk=nFilled(iS))
    dRho(:,:) = workLocal
    ! and hermitize
    call pblasfx_ptranc(workLocal, denseDesc%blacsOrbSqr, dRho, denseDesc%blacsOrbSqr,&
        & beta=(1.0_dp,0.0_dp))

  #:else

    ! serial case
    nOrb = size(dRho, dim = 1)

    if (allocated(idHam)) then
      call unpackHPauli(dHam, kPoint(:,iK), neighbourList%iNeighbour, nNeighbourSK, iSparseStart,&
          & denseDesc%iAtomStart, img2CentCell, iCellVec, cellVec, dRho, iHam=idHam)
    else
      call unpackHPauli(dHam, kPoint(:,iK), neighbourList%iNeighbour, nNeighbourSK, iSparseStart,&
          & denseDesc%iAtomStart, img2CentCell, iCellVec, cellVec, dRho)
    end if

    ! form |c> H' <c|
    call hemm(workLocal, 'l', dRho, eigVecsCplx(:,:,iKS))
    workLocal(:,:) = matmul(transpose(conjg(eigVecsCplx(:,:,iKS))), workLocal)

    ! diagonal elements of workLocal are now derivatives of eigenvalues if needed
    if (allocated(dEi)) then
      do ii = 1, nOrb
        dEi(ii, iK, iS, iCart) = real(workLocal(ii,ii),dp)
      end do
    end if

    if (omega == 0.0_dp) then
      ! static case

      ! Form actual perturbation U matrix for eigenvectors (potentially at finite T) by
      ! weighting the elements
      do iFilled = 1, nFilled(1)
        do iEmpty = nEmpty(1), nOrb
          workLocal(iEmpty, iFilled) = workLocal(iEmpty, iFilled) * &
              & invDiff(eigvals(iFilled, iK, 1), eigvals(iEmpty, iK, 1), Ef(1), tempElec)&
              & *theta(eigvals(iFilled, iK, 1), eigvals(iEmpty, iK, 1), tempElec)
        end do
      end do

      ! calculate the derivatives of ci
      workLocal(:, :nFilled(1)) =&
          & matmul(eigVecsCplx(:, nEmpty(1):, iKS), workLocal(nEmpty(1):, :nFilled(1)))

      if (allocated(dPsi)) then
        dPsi(:, :, iK, iS, iCart) = workLocal
      end if

      ! zero the uncalculated virtual states
      workLocal(:, nFilled(1)+1:) = 0.0_dp

      ! form the derivative of the density matrix
      dRho(:,:) = matmul(workLocal(:, :nFilled(1)),&
          & transpose(conjg(eigVecsCplx(:, :nFilled(1), iKS))) )&
          & + matmul(eigVecsCplx(:, :nFilled(1), iKS),&
          & transpose(conjg(workLocal(:, :nFilled(iKS)))) )

    else

      ! dynamic case

      dRho(:,:) = 0.0_dp
      ! loop over real and imaginary frequencies of driving
      do iSignOmega = -1, 1, 2

        ! Form actual perturbation U matrix for eigenvectors (potentially at finite T) by weighting
        ! the elements
        do iFilled = 1, nFilled(iS)
          do iEmpty = nEmpty(iS), nOrb
            cWorkLocal(iEmpty, iFilled) = workLocal(iEmpty, iFilled) * &
                & ( theta(eigvals(iEmpty, iK, iS), Ef(iS), tempElec)&
                & - theta(eigvals(iFilled, iK, iS), Ef(iS), tempElec) ) &
                & * theta(eigvals(iFilled, iK, iS), eigvals(iEmpty, iK, iS), tempElec)&
                & / ( eigvals(iEmpty, iK, iS) - eigvals(iFilled, iK, iS) + iSignOmega * omega + eta)
          end do
        end do

        ! calculate the derivatives of ci
        cWorkLocal(:, :nFilled(1)) =&
            & matmul(eigVecsCplx(:, nEmpty(1):, iKS), cWorkLocal(nEmpty(1):, :nFilled(1)))

        ! form the derivative of the density matrix
        dRho(:,:) = dRho + 0.5_dp * ( matmul(cWorkLocal(:, :nFilled(1)),&
            & transpose(conjg(eigVecsCplx(:, :nFilled(1), iKS))) )&
            & + matmul(eigVecsCplx(:, :nFilled(1), iKS),&
            & transpose(conjg(cWorkLocal(:, :nFilled(iKS)))) ) )

      end do

    end if

  #:endif

    dRhoSparse(:,:) = 0.0_dp
    if (allocated(idRhoSparse)) then
      idRhoSparse(:,:) = 0.0_dp
    end if

  #:if WITH_SCALAPACK
    if (allocated(idRhoSparse)) then
      call packRhoPauliBlacs(env%blacs, denseDesc, dRho, kPoint(:,iK), kWeight(iK),&
          & neighbourList%iNeighbour, nNeighbourSK, orb%mOrb, iCellVec, cellVec, iSparseStart,&
          & img2CentCell, dRhoSparse, idRhoSparse)
      else
        call packRhoPauliBlacs(env%blacs, denseDesc, dRho, kPoint(:,iK), kWeight(iK),&
            & neighbourList%iNeighbour, nNeighbourSK, orb%mOrb, iCellVec, cellVec, iSparseStart,&
            & img2CentCell, dRhoSparse)
      end if
  #:else
      call packHSPauli(dRhoSparse, dRho, neighbourlist%iNeighbour, nNeighbourSK, orb%mOrb,&
          & denseDesc%iAtomStart, iSparseStart, img2CentCell)
      if (allocated(idRhoSparse)) then
        call packHSPauliImag(idRhoSparse, dRho, neighbourlist%iNeighbour, nNeighbourSK,&
            & orb%mOrb, denseDesc%iAtomStart, iSparseStart, img2CentCell)
      end if
  #:endif

      ! adjustment from Pauli to charge/spin
      dRhoSparse(:,:) = 2.0_dp * dRhoSparse

  end subroutine dRhoStaticPauli


  !> Calculate the change in the density matrix due to shift in the Fermi energy
  subroutine dRhoFermiChangeStaticPauli(dRhoExtra, idRhoExtra, env, parallelKS, iKS, kPoint,&
      & kWeight, iCellVec, cellVec, neighbourList, nNEighbourSK, img2CentCell, iSparseStart, dEf,&
      & Ef, nFilled, nEmpty, eigVecsCplx, orb, denseDesc, tempElec, eigVals&
    #:if WITH_SCALAPACK
      &, desc&
    #:endif
      &)

    !> Additional contribution to the density matrix to cancel effect of Fermi energy change
    real(dp), intent(out) :: dRhoExtra(:,:)

    !> Imaginary part of additional contribution to the density matrix to cancel effect of Fermi
    !> energy change
    real(dp), intent(inout), allocatable :: idRhoExtra(:,:)

    !> Environment settings
    type(TEnvironment), intent(inout) :: env

    !> K-points and spins to process
    type(TParallelKS), intent(in) :: parallelKS

    !> spin/kpoint channel
    integer, intent(in) :: iKS

    !> k-points
    real(dp), intent(in) :: kPoint(:,:)

    !> Weights for k-points
    real(dp), intent(in) :: kWeight(:)

    !> Vectors (in units of the lattice constants) to cells of the lattice
    real(dp), intent(in) :: cellVec(:,:)

    !> Index for which unit cell atoms are associated with
    integer, intent(in) :: iCellVec(:)

    !> list of neighbours for each atom
    type(TNeighbourList), intent(in) :: neighbourList

    !> Number of neighbours for each of the atoms
    integer, intent(in) :: nNeighbourSK(:)

    !> map from image atoms to the original unique atom
    integer, intent(in) :: img2CentCell(:)

    !> Index array for the start of atomic blocks in sparse arrays
    integer, intent(in) :: iSparseStart(:,:)

    !> Fermi level derivative
    real(dp), intent(in) :: dEf(:)

    !> Fermi level
    real(dp), intent(in) :: Ef(:)

    !> Last (partly) filled level in each spin channel
    integer, intent(in) :: nFilled(:)

    !> First (partly) empty level in each spin channel
    integer, intent(in) :: nEmpty(:)

    !> ground state eigenvectors
    complex(dp), intent(in) :: eigVecsCplx(:,:,:)

    !> Atomic orbital information
    type(TOrbitals), intent(in) :: orb

    !> Dense matrix descriptor
    type(TDenseDescr), intent(in) :: denseDesc

    !> Electron temperature
    real(dp), intent(in) :: tempElec

    !> Eigenvalue of each level, kpoint and spin channel
    real(dp), intent(in) :: eigvals(:,:,:)

  #:if WITH_SCALAPACK
    !> BLACS matrix descriptor
    integer, intent(in) :: desc(DLEN_)
  #:endif

    integer :: iFilled, jj, jGlob, iK, iS
    complex(dp) :: workLocal(size(eigVecsCplx, dim=1), size(eigVecsCplx, dim=2))
  #:if WITH_SCALAPACK
    complex(dp) :: workLocal2(size(eigVecsCplx, dim=1), size(eigVecsCplx, dim=2))
    complex(dp) :: workLocal3(size(eigVecsCplx, dim=1), size(eigVecsCplx, dim=2))
  #:endif

    workLocal(:,:) = cmplx(0,0,dp)

    iK = parallelKS%localKS(1, iKS)
    iS = parallelKS%localKS(2, iKS)

  #:if WITH_SCALAPACK

    do jj = 1, size(workLocal,dim=2)
      jGlob = scalafx_indxl2g(jj, desc(NB_), env%blacs%orbitalGrid%mycol, desc(CSRC_),&
          & env%blacs%orbitalGrid%ncol)
      if (jGlob >= nEmpty(iS) .and. jGlob <= nFilled(iS)) then
        workLocal(:, jj) = eigVecsCplx(:, jj, iKS) * &
            & deltamn(eigVals(jGlob, iK, iS), Ef(iS), tempElec) * dEf(iS)
      end if
    end do

    workLocal3(:,:) = eigVecsCplx(:,:,iKS)
    call pblasfx_pgemm(workLocal, denseDesc%blacsOrbSqr, workLocal3,&
        & denseDesc%blacsOrbSqr, workLocal2, denseDesc%blacsOrbSqr, transb="C")
    workLocal(:,:) = workLocal2
    ! Hermitian symmetry
    call pblasfx_ptranc(workLocal2, denseDesc%blacsOrbSqr, workLocal, denseDesc%blacsOrbSqr,&
        & alpha=(0.5_dp,0.0_dp), beta=(0.5_dp,0.0_dp))
  #:else

    do iFilled = nEmpty(1), nFilled(1)
      workLocal(:, iFilled) = eigVecsCplx(:, iFilled, 1) * &
          & deltamn(eigvals(iFilled, 1, 1), Ef(1), tempElec) * dEf(1)
    end do

    workLocal(:, :) = matmul(workLocal(:, nEmpty(1):nFilled(1)),&
        & transpose(conjg(eigVecsCplx(:, nEmpty(1):nFilled(1), 1))))

  #:endif

    if (allocated(idRhoExtra)) then
      idRhoExtra(:,:) = 0.0_dp
    end if

    ! pack extra term into density matrix
  #:if WITH_SCALAPACK
    if (allocated(idRhoExtra)) then
      call packRhoPauliBlacs(env%blacs, denseDesc, workLocal, kPoint(:,iK), kWeight(iK),&
          & neighbourList%iNeighbour, nNeighbourSK, orb%mOrb, iCellVec, cellVec, iSparseStart,&
          & img2CentCell, dRhoExtra, idRhoExtra)
    else
      call packRhoPauliBlacs(env%blacs, denseDesc, workLocal, kPoint(:,iK), kWeight(iK),&
          & neighbourList%iNeighbour, nNeighbourSK, orb%mOrb, iCellVec, cellVec, iSparseStart,&
          & img2CentCell, dRhoExtra)
    end if
  #:else
    call packHSPauli(dRhoExtra, workLocal, neighbourlist%iNeighbour, nNeighbourSK, orb%mOrb,&
        & denseDesc%iAtomStart, iSparseStart, img2CentCell)
    if (allocated(idRhoExtra)) then
      call packHSPauliImag(idRhoExtra, workLocal, neighbourlist%iNeighbour, nNeighbourSK,&
          & orb%mOrb, denseDesc%iAtomStart, iSparseStart, img2CentCell)
    end if

  #:endif

    ! adjustment from Pauli to charge/spin
    dRhoExtra(:,:) = 2.0_dp * dRhoExtra
    if (allocated(idRhoExtra)) then
      idRhoExtra(:,:) = 2.0_dp * idRhoExtra
    end if

  end subroutine dRhoFermiChangeStaticPauli

end module dftbp_Perturbderivs
