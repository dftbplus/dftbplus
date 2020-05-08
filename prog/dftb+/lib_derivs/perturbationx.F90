!--------------------------------------------------------------------------------------------------!
!  DFTB+: general package for performing fast atomistic simulations                                !
!  Copyright (C) 2018  DFTB+ developers group                                                      !
!                                                                                                  !
!  See the LICENSE file for terms of usage and distribution.                                       !
!--------------------------------------------------------------------------------------------------!

#:include 'common.fypp'

!> Module for linear response derivative calculations using perturbation methods
module dftbp_perturbxderivs
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
  use dftbp_thirdorder, only : TThirdOrder
  use dftbp_dftbplusu
  use dftbp_rangeseparated, only : TRangeSepFunc
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
  use dftbp_rotateDegenerateOrbs
  use dftbp_slakocont
  use dftbp_nonscc, only : TNonSccDiff
  use dftbp_taggedoutput

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
  public :: dPsidx

  !> Direction labels
  character(len=1), parameter :: direction(3) = ['x','y','z']

contains

  !> Static (frequency independent) perturbation at q=0
  subroutine dPsidx(env, parallelKS, filling, eigvals, eigVecsReal, eigVecsCplx, rhoPrim,&
      & potential, qOrb, q0, ham, over, skHamCont, skOverCont, nonSccDeriv, orb, nAtom, species,&
      & speciesnames, neighbourList, nNeighbourSK, denseDesc, iSparseStart, img2CentCell, coord,&
      & sccCalc, maxSccIter, sccTol, nMixElements, nIneqMixElements, iEqOrbitals, tempElec, Ef,&
      & tFixEf, spinW, thirdOrd, tDFTBU,UJ, nUJ, iUJ, niUJ, iEqBlockDftbu, onsMEs, iEqBlockOnSite,&
      & rangeSep, nNeighbourLC, pChrgMixer, taggedWriter, tWriteAutoTest, autoTestTagFile,&
      & tWriteTaggedOut, taggedResultsFile, tWriteDetailedOut, fdDetailedOut, kPoint, kWeight,&
      & iCellVec, cellVec, tPeriodic, tMulliken)

    !> Environment settings
    type(TEnvironment), intent(inout) :: env

    !> K-points and spins to process
    type(TParallelKS), intent(in) :: parallelKS

    !> Fillings of unperturbed system
    real(dp), intent(in) :: filling(:,:,:)

    !> Eigenvalue of each level, kpoint and spin channel
    real(dp), intent(in) :: eigvals(:,:,:)

    !> ground state eigenvectors
    real(dp), intent(in), allocatable :: eigVecsReal(:,:,:)

    !> ground state complex eigenvectors
    complex(dp), intent(in), allocatable :: eigvecsCplx(:,:,:)

    !> Unperturbed density matrix in sparse format
    real(dp), intent(in) :: rhoPrim(:,:)

    !> Unperturbed potentials
    type(TPotentials), intent(in) :: potential

    !> Electrons in each atomic orbital
    real(dp), intent(in) :: qOrb(:,:,:)

    !> reference charges
    real(dp), intent(in) :: q0(:,:,:)

    !> Sparse Hamiltonian
    real(dp), intent(in) :: ham(:,:)

    !> sparse overlap matrix
    real(dp), intent(in) :: over(:)

    !> Container for SK Hamiltonian integrals
    type(TSlakoCont), intent(in) :: skHamCont

    !> Container for SK overlap integrals
    type(TSlakoCont), intent(in) :: skOverCont

    !> method for calculating derivatives of S and H0
    type(TNonSccDiff), intent(in) :: nonSccDeriv

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
    type(TScc), intent(inout), allocatable :: sccCalc

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
    type(TThirdOrder), allocatable, intent(inout) :: thirdOrd

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
    type(TRangeSepFunc), allocatable, intent(inout) :: rangeSep

    !> Number of neighbours for each of the atoms for the exchange contributions in the long range
    !> functional
    integer, intent(inout), allocatable :: nNeighbourLC(:)

    !> Charge mixing object
    type(TMixer), intent(inout) :: pChrgMixer

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

    !> Should Mulliken populations be generated/output
    logical, intent(in) :: tMulliken

    integer :: iS, iK, iKS, iAt, iNeigh, iCart, iSCC, iLev, iSh, iSp, jAt, jAtf, iOrb, jCart

    integer :: nSpin, nKpts, nOrbs, nIndepHam

    ! maximum allowed number of electrons in a single particle state
    real(dp) :: maxFill

    integer, allocatable :: nFilled(:,:), nEmpty(:,:)

    integer :: ii, jj, iGlob, jGlob
    integer :: iSCCIter
    logical :: tStopSCC

    ! matrices for derivatives of terms in hamiltonian and outputs
    real(dp), allocatable :: dHam(:,:), idHam(:,:), dOver(:,:), dH0(:,:)

    ! overlap derivative terms in potential, omega dS + d(delta q gammma) S
    real(dp), allocatable :: sOmega(:,:,:)

    real(dp) :: drho(size(over),size(ham, dim=2))
    real(dp) :: drhoExtra(size(over),size(ham, dim=2))
    real(dp), allocatable :: idRho(:,:), idRhoExtra(:,:)
    real(dp) :: dqIn(orb%mOrb,nAtom,size(ham, dim=2))
    real(dp) :: dqOut(orb%mOrb, nAtom, size(ham, dim=2), 3, nAtom)
    real(dp) :: dqUpDown(orb%mOrb, size(ham, dim=2))
    real(dp) :: dqInpRed(nMixElements), dqOutRed(nMixElements)
    real(dp) :: dqDiffRed(nMixElements), sccErrorQ
    real(dp) :: dqPerShell(orb%mShell,nAtom,size(ham, dim=2))

    ! eigenvalue weighted vectors
    real(dp), allocatable :: eCiReal(:, :, :)
    complex(dp), allocatable :: eCiCplx(:, :, :)

    real(dp), allocatable :: Vat(:,:), vdgamma(:,:,:)

    real(dp), allocatable :: dqBlockIn(:,:,:,:), SSqrReal(:,:)
    real(dp), allocatable :: dqBlockOut(:,:,:,:)
    real(dp), allocatable :: dummy(:,:,:,:)

    ! derivative of potentials
    type(TPotentials) :: dPotential

    real(dp), allocatable :: shellPot(:,:,:), atomPot(:,:)

    logical :: tSccCalc, tConverged
    logical, allocatable :: tMetallic(:)

    real(dp), allocatable :: dEi(:,:,:,:)
    real(dp), allocatable :: dPsiReal(:,:,:,:)
    complex(dp), allocatable :: dPsiCmplx(:,:,:,:,:)

    integer :: fdResults, fdResponses

    ! used for range separated contributions, note this stays in the up/down representation
    ! throughout if spin polarised
    real(dp), pointer :: dRhoOutSqr(:,:,:), dRhoInSqr(:,:,:)
    real(dp), allocatable, target :: dRhoOut(:), dRhoIn(:)

    ! non-variational part of charge derivative
    real(dp), allocatable :: dqNonVariational(:,:,:), dqNonVariationalBlock(:,:,:,:)

    real(dp) :: dDipole(3)

  #:if WITH_SCALAPACK
    ! need distributed matrix descriptors
    integer :: desc(DLEN_), nn

    type(blocklist) :: blocks
    integer :: blockSize, iLoc

    nn = denseDesc%fullSize
    call scalafx_getdescriptor(env%blacs%orbitalGrid, nn, nn, env%blacs%rowBlockSize,&
        & env%blacs%columnBlockSize, desc)
    call blocks%init(env%blacs%orbitalGrid, desc, "c")

  #:endif

    if (tFixEf) then
      call error("Perturbation expressions not currently implemented for fixed Fermi energy")
    end if

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

    allocate(tMetallic(nIndepHam))

    nOrbs = size(filling,dim=1)

    nKpts = size(filling,dim=2)

    allocate(dEi(nOrbs, nAtom, nSpin, 3))

    allocate(dqNonVariational(orb%mOrb,nAtom,nSpin))

  #:if WITH_SCALAPACK

    if (allocated(eigVecsReal)) then

      allocate(eCiReal(size(eigVecsReal,dim=1), size(eigVecsReal,dim=2), size(eigVecsReal,dim=3)))

      ! e_i |c_i>
      do iKS = 1, parallelKS%nLocalKS
        iK = parallelKS%localKS(1, iKS)
        iS = parallelKS%localKS(2, iKS)
        do ii = 1, size(blocks)
          call blocks%getblock(ii, iGlob, iLoc, blockSize)
          iK = parallelKS%localKS(1, iKS)
          iS = parallelKS%localKS(2, iKS)
          do jj = 0, min(blockSize - 1, nOrbs - iGlob)
            eCiReal(:, iLoc + jj, iKS) = eigVals(iGlob + jj, iK, iS)&
                & * eigVecsReal(:, iLoc + jj, iKS)
          end do
        end do
      end do

    else

      allocate(eCiCplx(size(eigVecsCplx,dim=1), size(eigVecsCplx,dim=2), size(eigVecsCplx,dim=3)))

      ! e_i |c_i>
      do iKS = 1, parallelKS%nLocalKS
        iK = parallelKS%localKS(1, iKS)
        iS = parallelKS%localKS(2, iKS)
        do ii = 1, size(blocks)
          call blocks%getblock(ii, iGlob, iLoc, blockSize)
          do jj = 0, min(blockSize - 1, nOrbs - iGlob)
            eCiCplx(:, iLoc + jj, iKS) = eigVals(iGlob + jj, iK, iS)&
                & * eigVecsCplx(:, iLoc + jj, iKS)
          end do
        end do
      end do

    end if

  #:else

    if (allocated(eigVecsReal)) then
      allocate(eCiReal(size(eigVecsReal,dim=1), size(eigVecsReal,dim=2), size(eigVecsReal,dim=3)))
      do iKS = 1, size(eigVecsReal,dim=3)
        do iOrb = 1, size(eigVecsReal,dim=2)
          eCiReal(:,iOrb, iKS) = eigVecsReal(:,iOrb, iKS) * eigVals(iOrb, 1, iKS)
        end do
      end do
    else
      allocate(eCiCplx(size(eigVecsCplx,dim=1), size(eigVecsCplx,dim=2), size(eigVecsCplx,dim=3)))
      do iKS = 1, parallelKS%nLocalKS
        iK = parallelKS%localKS(1, iKS)
        iS = parallelKS%localKS(2, iKS)
        do iOrb = 1, size(eigVecsCplx,dim=2)
          eCiCplx(:,iOrb, iKS) = eigVecsCplx(:,iOrb, iKS) * eigVals(iOrb, iK, iS)
        end do
      end do
    end if

  #:endif

    allocate(dHam(size(ham,dim=1),nSpin))
    allocate(dOver(size(ham,dim=1),3))
    allocate(dH0(size(ham,dim=1),3))

    tSccCalc = allocated(sccCalc)

    ! terms v S' and v' S
    if (tSccCalc) then
      allocate(sOmega(size(ham,dim=1),nSpin,2))
      allocate(Vat(nAtom,nSpin))
      allocate(vdgamma(orb%mShell,nAtom,nSpin))
    end if

    if (allocated(rangeSep)) then
    #:if WITH_SCALAPACK
      call error("Range separation not supported for MPI at the moment")
    #:endif
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

    if (tDFTBU .or. allocated(onsMEs)) then
      allocate(dqBlockIn(orb%mOrb,orb%mOrb,nAtom,nSpin))
      allocate(dqBlockOut(orb%mOrb,orb%mOrb,nAtom,nSpin))
      allocate(dqNonVariationalBlock(orb%mOrb,orb%mOrb,nAtom,nSpin))
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

    do iS = 1, nIndepHam
      tMetallic(iS) = (.not.all(nFilled(iS,:) == nEmpty(iS,:) -1))
      !write(stdOut,*)'Fractionally filled range'
      !do iK = 1, nKpts
      !  write(stdOut,*) nEmpty(:,iK), ':', nFilled(:,iK)
      !end do
    end do

    dqOut(:,:,:,:,:) = 0.0_dp
    dEi(:,:,:,:) = 0.0_dp

    ! Displaced atom to differentiate wrt
    lpAtom: do iAt = 1, nAtom

      call nonSccDeriv%getFirstDeriv(dOver, env, skOverCont, coord, species, iAt, orb,&
          & nNeighbourSK, neighbourList%iNeighbour, iSparseStart, img2centcell)

      call nonSccDeriv%getFirstDeriv(dH0, env, skHamCont, coord, species, iAt, orb, nNeighbourSK,&
          & neighbourList%iNeighbour, iSparseStart, img2centcell)

      ! perturbation direction
      lpCart: do iCart = 1, 3

        write(stdOut,*)'Calculating derivative for displacement along ', &
            & trim(direction(iCart)),' for atom', iAt

        if (tSccCalc) then
          sOmega(:,:,:) = 0.0_dp
          ! First part, omega dS
          call add_shift(sOmega(:,:,1), dOver(:,iCart), nNeighbourSK, neighbourList%iNeighbour,&
              & species, orb, iSparseStart, nAtom, img2CentCell, potential%intBlock)
        end if

        if (tSccCalc) then

          vdgamma(:,:,:) = 0.0_dp
          vat(:,:) = 0.0_dp

          call sccCalc%updateCoords(env, coord, species, neighbourList)
          call sccCalc%updateCharges(env, qOrb, orb, species, q0)
          call sccCalc%addPotentialDeriv(env, vAt(:,1), vdgamma, species, neighbourList%iNeighbour,&
              & img2CentCell, coord, orb, iCart, iAt)
          call total_shift(vdgamma, vAt, orb, species)

        end if

        ! non-variational part of the charge change due to basis derivatives
        if (tMulliken) then
          dqNonVariational(:,:,:) = 0.0_dp
          do iS = 1, nSpin
            call mulliken(dqNonVariational(:,:,iS), dOver(:,iCart), rhoPrim(:,iS), orb,&
                & neighbourList%iNeighbour, nNeighbourSK, img2CentCell, iSparseStart)
          end do
          if (tDFTBU .or. allocated(onsMEs)) then
            dqNonVariationalBlock(:,:,:,:) = 0.0_dp
            do iS = 1, nSpin
              call mulliken(dqNonVariationalBlock(:,:,:,iS), dOver(:,iCart), rhoPrim(:,iS), orb,&
                  & neighbourList%iNeighbour, nNeighbourSK, img2CentCell, iSparseStart)
            end do
          end if
        end if

        dqIn(:,:,:) = 0.0_dp
        if (tDFTBU .or. allocated(onsMEs)) then
          dqBlockIn(:,:,:,:) = 0.0_dp
          dqBlockOut(:,:,:,:) = 0.0_dp
        end if

        dPotential%extAtom(:,:) = 0.0_dp
        dPotential%extShell(:,:,:) = 0.0_dp
        dPotential%extBlock(:,:,:,:) = 0.0_dp

        if (tSccCalc) then
          call reset(pChrgMixer, nMixElements)
          dqInpRed(:) = 0.0_dp
          dqPerShell(:,:,:) = 0.0_dp
          if (allocated(rangeSep)) then
            dRhoIn(:) = 0.0_dp
            dRhoOut(:) = 0.0_dp
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
            call sccCalc%updateCharges(env, dqIn+dqNonVariational, orb, species)
            call sccCalc%updateShifts(env, orb, species, neighbourList%iNeighbour, img2CentCell)
            call sccCalc%getShiftPerAtom(dPotential%intAtom(:,1))
            call sccCalc%getShiftPerL(dPotential%intShell(:,:,1))

            if (allocated(spinW)) then
              call getChargePerShell(dqIn+dqNonVariational, orb, species, dqPerShell)
              call getSpinShift(shellPot, dqPerShell, species, orb, spinW)
              dPotential%intShell(:,:,:) = dPotential%intShell + shellPot
            end if

            if (allocated(thirdOrd)) then
              atomPot(:,:) = 0.0_dp
              shellPot(:,:,:) = 0.0_dp
              call thirdOrd%getdShiftdQ(atomPot(:,1), shellPot(:,:,1), species, neighbourList,&
                  & dqIn+dqNonVariational, img2CentCell, orb)
              dPotential%intAtom(:,1) = dPotential%intAtom(:,1) + atomPot(:,1)
              dPotential%intShell(:,:,1) = dPotential%intShell(:,:,1)&
                  & + shellPot(:,:,1)
            end if

            if (tDFTBU) then
              ! note the derivatives of both FLL and pSIC are the same (pSIC, i.e. case 2 in module)
              call getDftbUShift(dPotential%orbitalBlock, dqBlockIn+dqNonVariationalBlock, species,&
                  & orb, 2, UJ, nUJ, niUJ, iUJ)
            end if
            if (allocated(onsMEs)) then
              ! onsite corrections
              call addOnsShift(dPotential%orbitalBlock, dPotential%iOrbitalBlock,&
                  & dqBlockIn + dqNonVariationalBlock, dummy, onsMEs, species, orb)
            end if

          end if

          call total_shift(dPotential%intShell,dPotential%intAtom, orb,species)
          call total_shift(dPotential%intBlock,dPotential%intShell, orb,species)
          if (tDFTBU .or. allocated(onsMEs)) then
            dPotential%intBlock(:,:,:,:) = dPotential%intBlock + dPotential%orbitalBlock
          end if
          dPotential%intBlock(:,:,:,:) = dPotential%intBlock + dPotential%extBlock

          if (tSccCalc) then
            sOmega(:,:,2) = 0.0_dp
            ! add the (Delta q) * d gamma / dx term
            call add_shift(sOmega(:,:,2), over, nNeighbourSK, neighbourList%iNeighbour,&
                & species, orb, iSparseStart, nAtom, img2CentCell, vdgamma)
            ! and add gamma * d (Delta q) / dx
            call add_shift(sOmega(:,:,2), over, nNeighbourSK, neighbourList%iNeighbour, species,&
                & orb, iSparseStart, nAtom, img2CentCell, dpotential%intBlock)
          end if

          dHam = 0.0_dp

          dHam(:,1) = dH0(:,iCart)

          if (tSccCalc) then
            dHam(:,:) = dHam + sOmega(:,:,1) + sOmega(:,:,2)
          end if

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
            drho(:,:) = 0.0_dp
            do iKS = 1, parallelKS%nLocalKS

              iS = parallelKS%localKS(2, iKS)

              if (allocated(dRhoOut)) then
                ! replace with case that will get updated in dRhoReal
                dRhoOutSqr(:,:,iS) = dRhoInSqr(:,:,iS)
              end if

              call dRhoReal(env, dHam, dOver(:,iCart), neighbourList, nNeighbourSK,&
                  & iSparseStart, img2CentCell, denseDesc, iKS, parallelKS, nFilled(:,1),&
                  & nEmpty(:,1), eigVecsReal, eigVals, Ef, tempElec, orb, drho(:,iS), iCart,&
                  & dRhoOutSqr, rangeSep, over, nNeighbourLC, eCiReal, tMetallic,&
                  & filling / maxFill,&
                #:if WITH_SCALAPACK
                  & desc,&
                #:endif
                  & dEi, dPsiReal, iAt)
            end do

          elseif (nSpin > 2) then

            do iKS = 1, parallelKS%nLocalKS

              iK = parallelKS%localKS(1, iKS)

              call dRhoPauli(env, dHam, idHam, dOver(:,iCart), neighbourList, nNeighbourSK,&
                  & iSparseStart, img2CentCell, denseDesc, parallelKS, nFilled(:, iK),&
                  & nEmpty(:, iK), eigvecsCplx, eigVals, Ef, tempElec, orb, dRho, idRho, kPoint,&
                  & kWeight, iCellVec, cellVec, iKS, iCart,&
                #:if WITH_SCALAPACK
                  & desc,&
                #:endif
                  & dEi, dPsiCmplx, iAt)

            end do

          else

            call error("Shouldn't be here")

          end if

        #:if WITH_SCALAPACK
          ! Add up and distribute density matrix contributions from each group
          call mpifx_allreduceip(env%mpi%globalComm, dRho, MPI_SUM)
        #:endif


          dRho(:,:) = maxFill * drho
          if (allocated(dRhoOut)) then
            dRhoOut(:) = maxFill * dRhoOut
          end if
          call ud2qm(dRho)

          if (allocated(idRho)) then
            idRho(:,:) = maxFill * drho
            call ud2qm(idRho)
          end if

          dqOut(:, :, :, iCart, iAt) = 0.0_dp
          do iS = 1, nSpin
            call mulliken(dqOut(:, :, iS, iCart, iAt), over, drho(:,iS), orb,&
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
              call OrbitalEquiv_reduce(dqOut(:, :, :, iCart, iAt), iEqOrbitals, orb,&
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
                  dqIn(:,:,:) = dqOut(:, :, :, iCart, iAt)
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
            do jAt = 1, nAtom
              iSp = species(jAt)
              do iSh = 1, orb%nShell(iSp)
                dqPerShell(iSh,jAt,:nSpin) = dqPerShell(iSh,jAt,:nSpin) +&
                    & sum(dqIn(orb%posShell(iSh,iSp): orb%posShell(iSh+1,iSp)-1,jAt,:nSpin),dim=1)
              end do
            end do

          end if

          iSCCIter = iSCCIter +1

        end do lpSCC

        dqOut(:, :, :, iCart, iAt) = dqOut(:, :, :, iCart, iAt) + dqNonVariational

      end do lpCart

    end do lpAtom

  #:if WITH_SCALAPACK
    call mpifx_allreduceip(env%mpi%globalComm, dEi, MPI_SUM)
  #:endif

    write(stdOut, *)'dEi'
    do iCart = 1, 3
      write(stdOut, *)iCart
      do iS = 1, nSpin
        do iAt = 1, nAtom
          write(stdOut, *) dEi(:, iAt, iS, iCart) ! * Hartree__eV
        end do
      end do
    end do

    if (tMulliken .or. tSccCalc) then
      write(stdOut, *)
      write(stdOut, *)'Charge derivatives'
      do iAt = 1, nAtom
        write(stdOut,"(A,I0)")'/d Atom_',iAt
        do iS = 1, nSpin
          do jAt = 1, nAtom
            write(stdOut, *)jAt, -sum(dqOut(:, jAt, iS, :, iAt), dim=1)
          end do
          write(stdOut, *)
        end do
      end do
      write(stdOut, *)

      write(stdOut, *)'Born effective charges'
      ! i.e. derivative of dipole moment wrt to atom positions, or equivalently derivative of forces
      ! wrt to a homogeneous electric field
      do iAt = 1, nAtom
        do iCart = 1, 3
          do jCart = 1, 3
            dDipole(jCart) = -sum(sum(dqOut(:, : , 1, iCart, iAt), dim=1) * coord(jCart, :))
          end do
          dDipole(iCart) = dDipole(iCart) -sum(qOrb(:,iAt,1) - q0(:,iAt,1))
          write(stdOut,*)dDipole
        end do
        write(stdOut, *)
      end do
      write(stdOut, *)

    end if

    if (tWriteAutoTest) then
      open(newunit=fdResults, file=autoTestTagFile, position="append")
      call taggedWriter%write(fdResults, tagLabels%dqdx, dqOut)
      close(fdResults)
    end if
    if (tWriteTaggedOut) then
      open(newunit=fdResults, file=taggedResultsFile, position="append")
      call taggedWriter%write(fdResults, tagLabels%dqdx, dqOut)
      close(fdResults)
    end if

  end subroutine dPsidx


  !> Calculate the derivative of density matrix from derivative of hamiltonian in static case at
  !> q=0, k=0
  subroutine dRhoReal(env, dHam, dOver, neighbourList, nNeighbourSK, iSparseStart,&
      & img2CentCell, denseDesc, iKS, parallelKS, nFilled, nEmpty, eigVecsReal, eigVals, Ef,&
      & tempElec, orb, dRhoSparse, iCart, dRhoSqr, rangeSep, over, nNeighbourLC, eCiReal,&
      & tMetallic, filling,&
    #:if WITH_SCALAPACK
      & desc,&
    #:endif
      & dEi, dPsi, iAtom)

    !> Environment settings
    type(TEnvironment), intent(inout) :: env

    !> Derivative of the hamiltonian
    real(dp), intent(in) :: dHam(:,:)

    !> Derivative of the overlap
    real(dp), intent(in) :: dOver(:)

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
    type(TRangeSepFunc), allocatable, intent(inout) :: rangeSep

    !> sparse overlap matrix
    real(dp), intent(in) :: over(:)

    !> Number of neighbours for each of the atoms for the exchange contributions in the long range
    !> functional
    integer, intent(inout), allocatable :: nNeighbourLC(:)

    !> Eigenvalue weighted eigenvectors
    real(dp), intent(in) :: eCiReal(:,:,:)

    !> Is this a metallic system
    logical, intent(in) :: tMetallic(:)

    !> Fillings of unperturbed system
    real(dp), intent(in) :: filling(:,:,:)

  #:if WITH_SCALAPACK
    !> BLACS matrix descriptor
    integer, intent(in) :: desc(DLEN_)
  #:endif

    !> Derivative of single particle eigenvalues
    real(dp), intent(out) :: dEi(:,:,:,:)

    !> Optional derivatives of single particle wavefunctions
    real(dp), allocatable, intent(inout) :: dPsi(:,:,:,:)

    !> Atom with which the the derivative is being calculated
    integer, intent(in) :: iAtom

    integer :: ii, jj, iGlob, jGlob, iFilled, iEmpty, iS, iK, nOrb
    real(dp) :: workLocal(size(eigVecsReal,dim=1), size(eigVecsReal,dim=2))
    real(dp) :: work2Local(size(eigVecsReal,dim=1), size(eigVecsReal,dim=2))
    real(dp) :: work3Local(size(eigVecsReal,dim=1), size(eigVecsReal,dim=2))
  #:if WITH_SCALAPACK
    real(dp) :: work4Local(size(eigVecsReal,dim=1), size(eigVecsReal,dim=2))
  #:endif
    real(dp), allocatable :: dRho(:,:)
    type(TDegeneracyTransform) :: transform

    real(dp), allocatable :: dFilling(:)

    iK = parallelKS%localKS(1, iKS)
    iS = parallelKS%localKS(2, iKS)

    if (tMetallic(iS)) then
      allocate(dFilling(size(dEi, dim=1)))
    end if

    call transform%init()

    dEi(:, iAtom, iS, iCart) = 0.0_dp
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
        & denseDesc%blacsOrbSqr, dRho, denseDesc%blacsOrbSqr)

    ! dS in square form
    call unpackHSRealBlacs(env%blacs, dOver, neighbourList%iNeighbour, nNeighbourSK,&
        & iSparseStart, img2CentCell, denseDesc, work2Local)

    ! H' - e S' <c|
    call pblasfx_psymm(work2Local, denseDesc%blacsOrbSqr, eCiReal(:,:,iKS), denseDesc%blacsOrbSqr,&
        & dRho, denseDesc%blacsOrbSqr, alpha=-1.0_dp, beta=1.0_dp)

    ! c_i times dH times c_i
    call pblasfx_pgemm(eigVecsReal(:,:,iKS), denseDesc%blacsOrbSqr, dRho,&
        & denseDesc%blacsOrbSqr, workLocal, denseDesc%blacsOrbSqr, transa="T")

    ! |c> S' <c|, note not fully efficient, as could replace second operation with pointwise product
    ! and sum along first index (distributed)
    call pblasfx_psymm(work2Local, denseDesc%blacsOrbSqr, eigVecsReal(:,:,iKS),&
        & denseDesc%blacsOrbSqr, work3local, denseDesc%blacsOrbSqr)
    call pblasfx_pgemm(eigVecsReal(:,:,iKS), denseDesc%blacsOrbSqr, work3Local,&
        & denseDesc%blacsOrbSqr, work4Local, denseDesc%blacsOrbSqr, transa="T")

    ! weight with inverse of energy differences
    do jj = 1, size(workLocal,dim=2)
      jGlob = scalafx_indxl2g(jj, desc(NB_), env%blacs%orbitalGrid%mycol, desc(CSRC_),&
          & env%blacs%orbitalGrid%ncol)
      do ii = 1, size(workLocal,dim=1)
        iGlob = scalafx_indxl2g(ii, desc(MB_), env%blacs%orbitalGrid%myrow, desc(RSRC_),&
            & env%blacs%orbitalGrid%nrow)
        ! derivative of eigenvalues stored in diagonal of matrix workLocal, from <c|h'|c>
        if (iGlob == jGlob) then
          !if (iGlob == jGlob) then workLocal(ii,jj) contains a derivative of an eigenvalue
          dEi(iGlob, iAtom, iS, iCart) = workLocal(ii,jj)
        end if
        if (iGlob == jGlob) then
          workLocal(ii,jj) = -0.5_dp * work4Local(ii, jj)
        else
          workLocal(ii,jj) = workLocal(ii,jj) / (eigvals(jGlob,1,iS) - eigvals(iGlob,1,iS))
        end if
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

    dRho(:,:) = 0.0_dp
    call unpackHS(dRho, dHam(:,iS), neighbourList%iNeighbour, nNeighbourSK,&
        & denseDesc%iAtomStart, iSparseStart, img2CentCell)

    if (allocated(rangeSep)) then
      call unpackHS(workLocal, over, neighbourList%iNeighbour, nNeighbourSK,&
          & denseDesc%iAtomStart, iSparseStart, img2CentCell)
      call rangeSep%addLRHamiltonian(env, dRhoSqr(:,:,iS), over, neighbourList%iNeighbour,&
          & nNeighbourLC, denseDesc%iAtomStart, iSparseStart, orb, dRho, workLocal)
    end if

    ! form H' |c>
    call symm(workLocal, 'l', dRho, eigVecsReal(:,:,iS))

    ! form H' - e S' |c>
    dRho(:,:) = 0.0_dp
    call unpackHS(dRho, dOver, neighbourList%iNeighbour, nNeighbourSK, denseDesc%iAtomStart,&
        & iSparseStart, img2CentCell)
    call symm(workLocal, 'l', dRho, eCiReal(:,:,iS), alpha=-1.0_dp, beta=1.0_dp)

    ! form |c> H' - e S' <c|
    workLocal(:,:) = matmul(transpose(eigVecsReal(:,:,iS)), workLocal)

    call transform%generateUnitary(workLocal, eigvals(:,iK,iS))
    call transform%degenerateTransform(workLocal)

    ! diagonal elements of workLocal are now derivatives of eigenvalues
    do ii = 1, nOrb
      dEi(ii, iAtom, iS, iCart) = workLocal(ii,ii)
    end do

    if  (tMetallic(iS)) then
      call dEida(dFilling, filling(:,iK,iS), dEi(:,iAtom, iS, iCart), tempElec)
      !write(stdOut,*)'dEf', dEfda(filling(:,iK,iS), dEi(:,iAtom, iS, iCart))
      !write(stdOut,*)dFilling
    end if

    work3Local = eigVecsReal(:,:,iS)
    call transform%applyUnitary(work3Local)

    dRho(:,:) = 0.0_dp
    call unpackHS(dRho, dOver, neighbourList%iNeighbour, nNeighbourSK, denseDesc%iAtomStart,&
        & iSparseStart, img2CentCell)
    call symm(work2Local, 'l', dRho, work3Local)
    work2Local(:,:) = work2Local * work3Local

    ! Form actual perturbation U matrix for eigenvectors by weighting the elements
    do iFilled = 1, nFilled(iS)
      do iEmpty = 1, nOrb
        if (iFilled == iEmpty) then
          workLocal(iFilled, iFilled) = -0.5_dp * sum(work2Local(:, iFilled))
        else
          if (.not.transform%degenerate(iFilled,iEmpty)) then
            workLocal(iEmpty, iFilled) = workLocal(iEmpty, iFilled)&
                & / (eigvals(iFilled, iK, iS) - eigvals(iEmpty, iK, iS))
          else
            workLocal(iEmpty, iFilled) = 0.0_dp
            workLocal(iFilled, iEmpty) = 0.0_dp
          end if
        end if
      end do
    end do

    ! calculate the derivatives of the eigenvectors
    workLocal(:, :nFilled(iS)) = matmul(work3Local, workLocal(:, :nFilled(iS)))

    if (allocated(dPsi)) then
      dPsi(:, :, iS, iCart) = work3Local
    end if

    do iFilled = 1, nOrb
      workLocal(:, iFilled) = workLocal(:, iFilled) * filling(iFilled, iK, iS)
    end do

    ! form the derivative of the density matrix
    dRho(:,:) = matmul(workLocal, transpose(work3Local)) + matmul(work3Local, transpose(workLocal))

    if (tMetallic(iS)) then
      ! extra contribution from change in Fermi level leading to change in occupations
      do iFilled = nEmpty(iS), nFilled(iS)
        workLocal(:, iFilled) = work3Local(:, iFilled) * dFilling(iFilled)
      end do
      dRho(:,:) = dRho + 0.5_dp * matmul(workLocal(:, nEmpty(iS):nFilled(iS)),&
          & transpose(work3Local(:, nEmpty(iS):nFilled(iS))))&
          & + 0.5 * matmul(work3Local(:, nEmpty(iS):nFilled(iS)),&
          & transpose(workLocal(:, nEmpty(iS):nFilled(iS))))
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

    call transform%destroy()

  end subroutine dRhoReal


  !> Calculate the derivative of density matrix from derivative of hamiltonian in static case at q=0
  subroutine dRhoPauli(env, dHam, idHam, dOver, neighbourList, nNeighbourSK, iSparseStart,&
      & img2CentCell, denseDesc, parallelKS, nFilled, nEmpty, eigVecsCplx, eigVals, Ef, tempElec,&
      & orb, dRhoSparse, idRhoSparse, kPoint, kWeight, iCellVec, cellVec, iKS, iCart,&
    #:if WITH_SCALAPACK
      & desc,&
    #:endif
      & dEi, dPsi, iAtom)

    !> Environment settings
    type(TEnvironment), intent(in) :: env

    !> Derivative of the hamiltonian
    real(dp), intent(in) :: dHam(:,:)

    !> Derivative of the imaginary part of the hamiltonian
    real(dp), intent(in), allocatable :: idHam(:,:)

    !> Derivative of the overlap matrix
    real(dp), intent(in) :: dOver(:)

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
    real(dp), intent(out) :: dEi(:,:,:,:)

    !> Optional derivatives of single particle wavefunctions
    complex(dp), allocatable, intent(inout) :: dPsi(:,:,:,:,:)


    integer, intent(in) :: iAtom

    integer :: ii, jj, iGlob, jGlob, iFilled, iEmpty, iK, iS, nOrb
    complex(dp) :: workLocal(size(eigVecsCplx,dim=1), size(eigVecsCplx,dim=2))
    complex(dp) :: dRho(size(eigVecsCplx,dim=1), size(eigVecsCplx,dim=2))

    iK = parallelKS%localKS(1, iKS)
    iS = parallelKS%localKS(2, iKS)

    dEi(:, iAtom, iS, iCart) = 0.0_dp
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

    ! derivative of eigenvalues stored in diagonal of matrix workLocal, from <c|h'|c>
    do jj = 1, size(workLocal,dim=2)
      jGlob = scalafx_indxl2g(jj, desc(NB_), env%blacs%orbitalGrid%mycol, desc(CSRC_),&
          & env%blacs%orbitalGrid%ncol)
      do ii = 1, size(workLocal,dim=1)
        iGlob = scalafx_indxl2g(ii, desc(MB_), env%blacs%orbitalGrid%myrow, desc(RSRC_),&
            & env%blacs%orbitalGrid%nrow)
        if (iGlob == jGlob) then
          dEi(iGlob, iAtom, iS, iCart) = real(workLocal(ii,jj),dp)
        end if
      end do
    end do

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
    do ii = 1, nOrb
      dEi(ii, iAtom, iS, iCart) = real(workLocal(ii,ii),dp)
    end do

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

  end subroutine dRhoPauli

end module dftbp_perturbxderivs
