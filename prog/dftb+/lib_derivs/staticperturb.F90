!--------------------------------------------------------------------------------------------------!
!  DFTB+: general package for performing fast atomistic simulations                                !
!  Copyright (C) 2006 - 2021  DFTB+ developers group                                               !
!                                                                                                  !
!  See the LICENSE file for terms of usage and distribution.                                       !
!--------------------------------------------------------------------------------------------------!

#:include 'common.fypp'

!> Module for static linear response derivative calculations using perturbation methods
module dftbp_staticperturb
  use dftbp_accuracy, only : dp, mc
  use dftbp_constants, only : Hartree__eV, quaternionName
  use dftbp_globalenv, only : stdOut
  use dftbp_message, only : error, warning
  use dftbp_commontypes, only : TOrbitals
  use dftbp_potentials, only : TPotentials, TPotentials_init
  use dftbp_scc, only : TScc
  use dftbp_orbitalequiv, only : OrbitalEquiv_reduce, OrbitalEquiv_expand
  use dftbp_populations, only : mulliken, densemulliken, getchargepershell
  use dftbp_spin, only : getSpinShift, ud2qm, qm2ud
  use dftbp_thirdorder, only : TThirdOrder
  use dftbp_dftbplusu, only : TDftbU, TDftbU_init, plusUFunctionals
  use dftbp_rangeseparated, only : TRangeSepFunc
  use dftbp_onsitecorrection, only : addOnsShift, onsblock_expand
  use dftbp_shift, only : add_shift, total_shift
  use dftbp_mixer, only : TMixer, mix, reset
  use dftbp_fermihelper, only : theta, deltamn, invDiff
  use dftbp_environment, only : TEnvironment
  use dftbp_periodic, only : TNeighbourList
  use dftbp_densedescr, only : TDenseDescr
  use dftbp_rotatedegen, only : TRotateDegen, TRotateDegen_init
  use dftbp_parallelks, only : TParallelKS
  use dftbp_blockpothelper, only : appendBlockReduced
#:if WITH_MPI
  use dftbp_mpifx
#:endif
#:if WITH_SCALAPACK
  use dftbp_sparse2dense
  use dftbp_scalapackfx
  use dftbp_scalafxext
#:else
  use dftbp_sparse2dense, only : unpackHS, packHS, packHSPauli, unpackHPauli, packHSPauliImag
  use dftbp_blasroutines, only : hemm, symm
#:endif

  implicit none

  private
  public :: staticPerturWrtE

  !> small complex value for frequency dependent cases
  complex(dp), parameter :: eta = (0.0_dp,1.0E-8_dp)

contains

  !> Static (frequency independent) perturbation at dq=0
  subroutine staticPerturWrtE(env, parallelKS, filling, eigvals, eigVecsReal, eigVecsCplx, ham,&
      & over, orb, nAtom, species, speciesnames, neighbourList, nNeighbourSK, denseDesc,&
      & iSparseStart, img2CentCell, coord, sccCalc, maxSccIter, sccTol, isSccConvRequired,&
      & nMixElements, nIneqMixElements, iEqOrbitals, tempElec, Ef, tFixEf, spinW, thirdOrd, dftbU,&
      & iEqBlockDftbu, onsMEs, iEqBlockOnSite, rangeSep, nNeighbourLC, pChrgMixer, kPoint, kWeight,&
      & iCellVec, cellVec, tPeriodic, polarisability, dEi, dqOut, neFermi, dEfdE)

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
    type(TScc), intent(inout), allocatable :: sccCalc

    !> maximal number of SCC iterations
    integer, intent(in) :: maxSccIter

    !> Tolerance for SCC convergence
    real(dp), intent(in) :: sccTol

    !> Use converged derivatives of charges
    logical, intent(in) :: isSccConvRequired

    !> nr. of elements to go through the mixer - may contain reduced orbitals and also orbital
    !> blocks (if a DFTB+U or onsite correction calculation)
    integer, intent(in) :: nMixElements

    !> nr. of inequivalent charges
    integer, intent(in) :: nIneqMixElements

    !> Equivalence relations between orbitals
    integer, intent(in), allocatable :: iEqOrbitals(:,:,:)

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

    !> Are there orbital potentials present
    type(TDftbU), intent(in), allocatable :: dftbU

    !> equivalence mapping for dual charge blocks
    integer, intent(in), allocatable :: iEqBlockDftbu(:,:,:,:)

    !> Data for range-separated calculation
    type(TRangeSepFunc), allocatable, intent(inout) :: rangeSep

    !> Number of neighbours for each of the atoms for the exchange contributions in the long range
    !> functional
    integer, intent(inout), allocatable :: nNeighbourLC(:)

    !> Charge mixing object
    type(TMixer), intent(inout), allocatable :: pChrgMixer

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

    !> Static electric polarisability
    real(dp), intent(out) :: polarisability(:,:)

    !> Derivatives of eigenvalues, if required
    real(dp), allocatable, intent(inout) :: dEi(:,:,:,:)

    !> Derivatives of Mulliken charges, if required
    real(dp), allocatable, intent(inout) :: dqOut(:,:,:,:)

    !> Number of electrons at the Fermi energy (if metallic)
    real(dp), allocatable, intent(inout) :: neFermi(:)

    !> Derivative of the Fermi energy (if metallic)
    real(dp), allocatable, intent(inout) :: dEfdE(:,:)

    integer :: iS, iK, iKS, iAt, iNeigh, iCart, iSCC, iLev, iSh, iSp, jAt, jAtf

    integer :: nSpin, nKpts, nOrbs, nIndepHam

    ! maximum allowed number of electrons in a single particle state
    real(dp) :: maxFill

    integer, allocatable :: nFilled(:,:), nEmpty(:,:)

    integer :: ii, jj, iGlob, jGlob, iEmpty, iFilled
    integer :: iSCCIter

    ! matrices for derivatives of terms in hamiltonian and outputs
    real(dp), allocatable :: dHam(:,:), idHam(:,:)
    real(dp) :: drho(size(over),size(ham, dim=2))
    real(dp) :: drhoExtra(size(over),size(ham, dim=2))
    real(dp), allocatable :: idRho(:,:), idRhoExtra(:,:)
    real(dp) :: dqIn(orb%mOrb,nAtom,size(ham, dim=2))
    real(dp) :: dqInpRed(nMixElements), dqOutRed(nMixElements)
    real(dp) :: dqDiffRed(nMixElements), sccErrorQ
    real(dp) :: dqPerShell(orb%mShell,nAtom,size(ham, dim=2))

    real(dp), allocatable :: dqBlockIn(:,:,:,:), SSqrReal(:,:)
    real(dp), allocatable :: dqBlockOut(:,:,:,:)
    real(dp), allocatable :: dummy(:,:,:,:)

    ! derivative of potentials
    type(TPotentials) :: dPotential

    real(dp), allocatable :: shellPot(:,:,:), atomPot(:,:)

    logical :: tSccCalc, tMetallic, tConverged

    real(dp), allocatable :: dPsiReal(:,:,:,:)
    complex(dp), allocatable :: dPsiCmplx(:,:,:,:,:)

    integer :: fdResults

    ! used for range separated contributions, note this stays in the up/down representation
    ! throughout if spin polarised
    real(dp), pointer :: dRhoOutSqr(:,:,:), dRhoInSqr(:,:,:)
    real(dp), allocatable, target :: dRhoOut(:), dRhoIn(:)

    !> For transformation in the  case of degeneracies
    type(TRotateDegen), allocatable :: transform(:)

  #:if WITH_SCALAPACK
    ! need distributed matrix descriptors
    integer :: desc(DLEN_), nn

    nn = denseDesc%fullSize
    call scalafx_getdescriptor(env%blacs%orbitalGrid, nn, nn, env%blacs%rowBlockSize,&
        & env%blacs%columnBlockSize, desc)
  #:endif

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

    allocate(transform(nIndepHam))
    do ii = 1, nIndepHam
      call TRotateDegen_init(transform(ii))
    end do

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

    tSccCalc = allocated(sccCalc)

    if (allocated(dftbU) .or. allocated(onsMEs)) then
      allocate(dqBlockIn(orb%mOrb,orb%mOrb,nAtom,nSpin))
      allocate(dqBlockOut(orb%mOrb,orb%mOrb,nAtom,nSpin))
    end if

    call TPotentials_init(dPotential,orb,nAtom,nSpin)

    allocate(nFilled(nIndepHam, nKpts))
    allocate(nEmpty(nIndepHam, nKpts))

    if (allocated(spinW) .or. allocated(thirdOrd)) then
      allocate(shellPot(orb%mShell, nAtom, nSpin))
    end if
    if (allocated(thirdOrd)) then
      allocate(atomPot(nAtom, nSpin))
    end if

    ! If derivatives of eigenvalues are needed
    if (allocated(dEi)) then
      dEi(:,:,:,:) = 0.0_dp
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

    ! should really check for each spin channel separately:
    tMetallic = (.not.all(nFilled == nEmpty -1))

    ! if derivatives of valence wavefunctions needed. Note these will have an arbitrary set of
    ! global phases
    ! if (allocated(eigVecsReal)) then
    !   allocate(dPsiReal(size(eigVecsReal,dim=1), size(eigVecsReal,dim=2), nIndepHam, 3))
    ! else
    !   allocate(dPsiCmplx(size(eigvecsCplx,dim=1), size(eigvecsCplx,dim=2), nKpts, nIndepHam, 3))
    ! end if

    if (tMetallic) then
      write(stdOut,*)'Metallic system'
    else
      write(stdOut,*)'Non-metallic system'
    end if

    if (tMetallic) then
      ! Density of electrons at the Fermi energy, required to correct later for shift in Fermi level
      ! at q=0 in metals
      if (allocated(neFermi)) then
        deallocate(neFermi)
        deallocate(dEfdE)
      end if
      allocate(neFermi(size(Ef)))
      allocate(dEfdE(size(Ef),3))

      do iS = 1, nIndepHam
        neFermi(iS) = 0.0_dp
        do iK = 1, nKpts
          do ii = nEmpty(iS, iK), nFilled(iS, iK)
            neFermi(iS) = neFermi(iS) + kWeight(iK) * deltamn(Ef(iS), eigvals(ii,iK,iS), tempElec)
          end do
        end do
        neFermi(iS) = maxFill * neFermi(iS)
      end do
      write(stdOut,*)'Density of states at the Fermi energy Nf (a.u.):', neFermi
    end if

    dqOut(:,:,:,:) = 0.0_dp

    ! polarisation direction
    ! note, could MPI parallelise over this
    lpCart: do iCart = 1, 3

      dqIn(:,:,:) = 0.0_dp
      if (allocated(dftbU) .or. allocated(onsMEs)) then
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

      lpSCC: do iSCCIter = 1, maxSccIter

        dPotential%intAtom(:,:) = 0.0_dp
        dPotential%intShell(:,:,:) = 0.0_dp
        dPotential%intBlock(:,:,:,:) = 0.0_dp

        if (allocated(dftbU) .or. allocated(onsMEs)) then
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

          if (allocated(dftbU)) then
            ! note the derivatives of both FLL and pSIC are the same (i.e. pSIC)
            call dftbU%getDftbUShift(dPotential%orbitalBlock, dqBlockIn, species, orb,&
                & plusUFunctionals%pSIC)
          end if
          if (allocated(onsMEs)) then
            ! onsite corrections
            call addOnsShift(dPotential%orbitalBlock, dPotential%iOrbitalBlock, dqBlockIn, dummy,&
                & onsMEs, species, orb)
          end if

        end if

        call total_shift(dPotential%intShell,dPotential%intAtom, orb,species)
        call total_shift(dPotential%intBlock,dPotential%intShell, orb,species)
        if (allocated(dftbU) .or. allocated(onsMEs)) then
          dPotential%intBlock(:,:,:,:) = dPotential%intBlock + dPotential%orbitalBlock
        end if
        dPotential%intBlock(:,:,:,:) = dPotential%intBlock + dPotential%extBlock

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

            call dRhoStaticReal(env, dHam, neighbourList, nNeighbourSK, iSparseStart,&
                & img2CentCell, denseDesc, iKS, parallelKS, nFilled(:,1), nEmpty(:,1),&
                & eigVecsReal, eigVals, Ef, tempElec, orb, drho(:,iS), iCart, dRhoOutSqr,&
                & rangeSep, over, nNeighbourLC, transform(iKS), &
                #:if WITH_SCALAPACK
                & desc,&
                #:endif
                & dEi, dPsiReal)
          end do

        elseif (nSpin > 2) then

          do iKS = 1, parallelKS%nLocalKS

            iK = parallelKS%localKS(1, iKS)

            call dRhoStaticPauli(env, dHam, idHam, neighbourList, nNeighbourSK,&
                & iSparseStart, img2CentCell, denseDesc, parallelKS, nFilled(:, iK),&
                & nEmpty(:, iK), eigvecsCplx, eigVals, Ef, tempElec, orb, dRho, idRho, kPoint,&
                & kWeight, iCellVec, cellVec, iKS, iCart, transform(iKS), &
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

            dqOut(:,:,iS,iCart) = 0.0_dp
            call mulliken(dqOut(:,:,iS,iCart), over, drho(:,iS), orb, &
                & neighbourList%iNeighbour, nNeighbourSK, img2CentCell, iSparseStart)

            dEfdE(iS, iCart) = -sum(dqOut(:, :, iS, iCart))
            if (abs(dEfdE(iS, iCart)) > epsilon(0.0_dp) .and. abs(neFermi(iS)) >  epsilon(0.0_dp))&
                & then
              dEfdE(iS, iCart) = -sum(dqOut(:, :, iS, iCart)) / neFermi(iS)
            else
              dEfdE(iS, iCart) = 0.0_dp
            end if

            if (abs(dEfdE(iS, iCart)) > 10.0_dp*epsilon(1.0_dp)) then
              ! Fermi level changes, so need to correct for the change in the number of charges

              if (allocated(eigVecsReal)) then

                ! real case, no k-points
                call dRhoFermiChangeStaticReal(dRhoExtra(:, iS), env, parallelKS, iKS,&
                    & neighbourList, nNeighbourSK, img2CentCell, iSparseStart, dEfdE(:, iCart),&
                    & Ef, nFilled(:,iK), nEmpty(:,iK), eigVecsReal, orb, denseDesc, tempElec,&
                    & eigVals, dRhoOutSqr&
                    #:if WITH_SCALAPACK
                    &, desc&
                    #:endif
                    &)

              elseif (nSpin > 2) then

                ! two component wavefunction cases
                call dRhoFermiChangeStaticPauli(dRhoExtra, idRhoExtra, env, parallelKS, iKS,&
                    & kPoint, kWeight, iCellVec, cellVec, neighbourList, nNEighbourSK,&
                    & img2CentCell, iSparseStart, dEfdE(:, iCart), Ef, nFilled(:,iK),&
                    & nEmpty(:,iK), eigVecsCplx, orb, denseDesc, tempElec, eigVals&
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

        dqOut(:,:,:,iCart) = 0.0_dp
        do iS = 1, nSpin
          call mulliken(dqOut(:,:,iS,iCart), over, drho(:,iS), orb, &
              & neighbourList%iNeighbour, nNeighbourSK, img2CentCell, iSparseStart)
          if (allocated(dftbU) .or. allocated(onsMEs)) then
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
            call OrbitalEquiv_reduce(dqOut(:,:,:,iCart), iEqOrbitals, orb,&
                & dqOutRed(:nIneqMixElements))
            if (allocated(dftbU)) then
              call AppendBlockReduced(dqBlockOut, iEqBlockDFTBU, orb, dqOutRed)
            end if
            if (allocated(onsMEs)) then
              call AppendBlockReduced(dqBlockOut, iEqBlockOnsite, orb, dqOutRed)
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
                dqIn(:,:,:) = dqOut(:,:,:,iCart)
                dqInpRed(:) = dqOutRed(:)
                if (allocated(dftbU) .or. allocated(onsMEs)) then
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

                if (allocated(dftbU) .or. allocated(onsMEs)) then
                  dqBlockIn(:,:,:,:) = 0.0_dp
                  if (allocated(dftbU)) then
                    call dftbU%expandBlock(dqInpRed, iEqBlockDFTBU, orb, dqBlockIn, species,&
                        & orbEquiv=iEqOrbitals)
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

      end do lpSCC

      if (tSccCalc .and. .not.tConverged) then
        call warning("SCC in perturbation is NOT converged, maximal SCC iterations exceeded")
        if (isSccConvRequired) then
          call env%shutdown()
        end if
      end if

      if (tMetallic) then
        write(stdOut,*)
        write(stdOut,"(A,2E20.12)")'d E_f / d E_'//trim(quaternionName(iCart+1))//':',&
            & dEfdE(:,iCart)
        write(stdOut,*)
      end if

      do ii = 1, 3
        polarisability(ii, iCart) = -sum(sum(dqOut(:,:nAtom,1,iCart),dim=1)*coord(ii,:nAtom))
      end do

    end do lpCart

  #:if WITH_SCALAPACK
    if (allocated(dEi)) then
      call mpifx_allreduceip(env%mpi%globalComm, dEi, MPI_SUM)
    end if
  #:endif

    write(stdOut,*)
    write(stdOut,*)'Static polarisability (a.u.)'
    do iCart = 1, 3
      write(stdOut,"(3E20.12)")polarisability(:, iCart)
    end do
    write(stdOut,*)

  end subroutine staticPerturWrtE


  !> Calculate the derivative of density matrix from derivative of hamiltonian in static case at
  !> q=0, k=0
  subroutine dRhoStaticReal(env, dHam, neighbourList, nNeighbourSK, iSparseStart, img2CentCell,&
      & denseDesc, iKS, parallelKS, nFilled, nEmpty, eigVecsReal, eigVals, Ef, tempElec, orb,&
      & dRhoSparse, iCart, dRhoSqr, rangeSep, over, nNeighbourLC, transform,&
    #:if WITH_SCALAPACK
      & desc,&
    #:endif
      & dEi, dPsi)

    !> Environment settings
    type(TEnvironment), intent(inout) :: env

    !> Derivative of the hamiltonian
    real(dp), intent(in) :: dHam(:,:)

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

    !> Transformation structure for degenerate orbitals
    type(TRotateDegen), intent(inout) :: transform

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
    real(dp), allocatable :: dRho(:,:)
    real(dp), allocatable :: eigVecsTransformed(:,:)
    logical :: isTransformed

    iK = parallelKS%localKS(1, iKS)
    iS = parallelKS%localKS(2, iKS)

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
        & denseDesc%blacsOrbSqr, dRho, denseDesc%blacsOrbSqr)

    ! c_i times dH times c_i
    call pblasfx_pgemm(eigVecsReal(:,:,iKS), denseDesc%blacsOrbSqr, dRho, denseDesc%blacsOrbSqr,&
        & workLocal, denseDesc%blacsOrbSqr, transa="T")

    eigvecsTransformed = eigVecsReal(:,:,iKS)
    call transform%generateUnitary(env, worklocal, eigvals(:,1,iS), eigVecsTransformed, denseDesc,&
        & isTransformed)
    ! now have states orthogonalised agains the operator in degenerate cases, |c~>

    if (isTransformed) then
      ! re-form |c~> H' <c~| with the transformed vectors

      dRho(:,:) = 0.0_dp
      workLocal(:,:) = 0.0_dp

      ! dH in square form
      call unpackHSRealBlacs(env%blacs, dHam(:,iS), neighbourList%iNeighbour, nNeighbourSK,&
          & iSparseStart, img2CentCell, denseDesc, workLocal)

      ! dH times c_i
      call pblasfx_psymm(workLocal, denseDesc%blacsOrbSqr, eigVecsTransformed,&
          & denseDesc%blacsOrbSqr, dRho, denseDesc%blacsOrbSqr)

      ! c_i times dH times c_i
      call pblasfx_pgemm(eigVecsTransformed, denseDesc%blacsOrbSqr, dRho, denseDesc%blacsOrbSqr,&
          & workLocal, denseDesc%blacsOrbSqr, transa="T")

    end if

    ! derivative of eigenvalues stored in diagonal of matrix workLocal, from <c|h'|c>
    if (allocated(dEi)) then
      do jj = 1, size(workLocal,dim=2)
        jGlob = scalafx_indxl2g(jj, desc(NB_), env%blacs%orbitalGrid%mycol, desc(CSRC_),&
            & env%blacs%orbitalGrid%ncol)
        do ii = 1, size(workLocal,dim=1)
          iGlob = scalafx_indxl2g(ii, desc(MB_), env%blacs%orbitalGrid%myrow, desc(RSRC_),&
              & env%blacs%orbitalGrid%nrow)
          if (iGlob == jGlob) then
            !if (iGlob == jGlob) then workLocal(ii,jj) contains a derivative of an eigenvalue
            dEi(iGlob, iK, iS, iCart) = workLocal(ii,jj)
          end if
        end do
      end do
    end if

    ! weight matrix with inverse of energy differences
    do jj = 1, size(workLocal,dim=2)
      jGlob = scalafx_indxl2g(jj, desc(NB_), env%blacs%orbitalGrid%mycol, desc(CSRC_),&
          & env%blacs%orbitalGrid%ncol)
      do ii = 1, size(workLocal,dim=1)
        iGlob = scalafx_indxl2g(ii, desc(MB_), env%blacs%orbitalGrid%myrow, desc(RSRC_),&
            & env%blacs%orbitalGrid%nrow)
        if (abs(eigvals(jGlob,iK,iS) - eigvals(iGlob,iK,iS)) < epsilon(0.0_dp)&
            & .and. iGlob /= jGlob) then
          ! degenerate, so no contribution
          workLocal(ii,jj) = 0.0_dp
        else
          workLocal(ii,jj) = workLocal(ii,jj) * &
              & invDiff(eigvals(jGlob,1,iS),eigvals(iGlob,1,iS),Ef(iS),tempElec)&
              & * theta(eigvals(jGlob,1,iS),eigvals(iGlob,1,iS),tempElec)
        end if
      end do
    end do

    ! Derivatives of states
    call pblasfx_pgemm(eigvecsTransformed, denseDesc%blacsOrbSqr, workLocal, denseDesc%blacsOrbSqr,&
        & dRho, denseDesc%blacsOrbSqr)

    if (allocated(dPsi)) then
      dPsi(:, :, iS, iCart) = workLocal
    end if

    ! Form derivative of occupied density matrix
    call pblasfx_pgemm(dRho, denseDesc%blacsOrbSqr, eigvecsTransformed, denseDesc%blacsOrbSqr,&
        & workLocal, denseDesc%blacsOrbSqr, transb="T", kk=nFilled(iS))

    dRho(:,:) = workLocal
    ! and symmetrize
    call pblasfx_ptran(workLocal, denseDesc%blacsOrbSqr, dRho, denseDesc%blacsOrbSqr, beta=1.0_dp)

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

    ! orthogonalise degenerate states against perturbation, producing |c~> H' <c~|
    call transform%generateUnitary(workLocal, eigvals(:,iK,iS))
    call transform%degenerateTransform(workLocal)

    ! diagonal elements of workLocal are now derivatives of eigenvalues if needed
    if (allocated(dEi)) then
      do ii = 1, nOrb
        dEi(ii, iK, iS, iCart) = workLocal(ii,ii)
      end do
    end if

    ! Form actual perturbation U matrix for eigenvectors (potentially at finite T) by weighting
    ! the elements
    do iFilled = 1, nFilled(iS)
      do iEmpty = nEmpty(iS), nOrb
        if (.not.transform%degenerate(iFilled,iEmpty) .or. iEmpty == iFilled) then
          workLocal(iEmpty, iFilled) = workLocal(iEmpty, iFilled) * &
              & invDiff(eigvals(iFilled, iK, iS), eigvals(iEmpty, iK, iS), Ef(iS), tempElec)&
              & *theta(eigvals(iFilled, iK, iS), eigvals(iEmpty, iK, iS), tempElec)
        else
          ! rotation should already have set these elements to zero
          workLocal(iEmpty, iFilled) = 0.0_dp
        end if
      end do
    end do

    eigvecsTransformed = eigVecsReal(:,:,iS)
    call transform%applyUnitary(eigvecsTransformed)

    ! calculate the derivatives of the eigenvectors
    workLocal(:, :nFilled(iS)) =&
        & matmul(eigvecsTransformed(:, nEmpty(iS):), workLocal(nEmpty(iS):, :nFilled(iS)))

    if (allocated(dPsi)) then
      dPsi(:, :, iS, iCart) = workLocal
    end if

    ! zero the uncalculated virtual states
    workLocal(:, nFilled(iS)+1:) = 0.0_dp

    ! form the derivative of the density matrix
    dRho(:,:) = matmul(workLocal(:, :nFilled(iS)), transpose(eigvecsTransformed(:, :nFilled(iS))))&
        & + matmul(eigvecsTransformed(:, :nFilled(iS)), transpose(workLocal(:, :nFilled(iS))))

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
      & nNEighbourSK, img2CentCell, iSparseStart, dEfdE, Ef, nFilled, nEmpty, eigVecsReal, orb,&
      & denseDesc, tempElec, eigVals, dRhoSqr&
    #:if WITH_SCALAPACK
      &, desc&
    #:endif
      &)

    !> Additional contribution to the density matrix to cancel effect of Fermi energy change
    real(dp), intent(inout) :: dRhoExtra(:)

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
    real(dp), intent(in) :: dEfdE(:)

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
            & deltamn(eigVals(jGlob, 1, iKS), Ef(iS), tempElec) * dEfdE(iS)
      end if
    end do
    call pblasfx_pgemm(workReal, denseDesc%blacsOrbSqr,eigVecsReal(:,:,iKS),&
        & denseDesc%blacsOrbSqr, workReal, denseDesc%blacsOrbSqr, transb="T",&
        & alpha=real(3-nSpin,dp))

  #:else

    do iFilled = nEmpty(iS), nFilled(iS)
      workReal(:, iFilled) = eigVecsReal(:, iFilled, iS) * &
          & deltamn(eigvals(iFilled, 1, iS), Ef(iS), tempElec) * dEfdE(iS)
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
  subroutine dRhoStaticPauli(env, dHam, idHam, neighbourList, nNeighbourSK, iSparseStart,&
      & img2CentCell, denseDesc, parallelKS, nFilled, nEmpty, eigVecsCplx, eigVals, Ef, tempElec,&
      & orb, dRhoSparse, idRhoSparse, kPoint, kWeight, iCellVec, cellVec, iKS, iCart, transform,&
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

    !> Transformation structure for degenerate orbitals
    type(TRotateDegen), intent(inout) :: transform

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
    complex(dp) :: dRho(size(eigVecsCplx,dim=1), size(eigVecsCplx,dim=2))
    complex(dp), allocatable :: eigVecsTransformed(:,:)
    logical :: isTransformed

    iK = parallelKS%localKS(1, iKS)
    iS = parallelKS%localKS(2, iKS)

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
        & denseDesc%blacsOrbSqr, dRho, denseDesc%blacsOrbSqr)

    ! c_i times dH times c_i
    call pblasfx_pgemm(eigVecsCplx(:,:,iKS), denseDesc%blacsOrbSqr, dRho,&
        & denseDesc%blacsOrbSqr, workLocal, denseDesc%blacsOrbSqr, transa="C")

    eigvecsTransformed = eigVecsCplx(:,:,iKS)
    call transform%generateUnitary(env, worklocal, eigvals(:,1,iS), eigVecsTransformed, denseDesc,&
        & isTransformed)
    ! now have states orthogonalised agains the operator in degenerate cases, |c~>

    if (isTransformed) then
      ! re-form |c~> H' <c~| with the transformed vectors

      dRho(:,:) = 0.0_dp
      workLocal(:,:) = 0.0_dp

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
      call pblasfx_phemm(workLocal, denseDesc%blacsOrbSqr, eigVecsTransformed,&
          & denseDesc%blacsOrbSqr, dRho, denseDesc%blacsOrbSqr)

      ! c_i times dH times c_i
      call pblasfx_pgemm(eigVecsTransformed, denseDesc%blacsOrbSqr, dRho,&
          & denseDesc%blacsOrbSqr, workLocal, denseDesc%blacsOrbSqr, transa="C")

    end if

    ! derivative of eigenvalues stored in diagonal of matrix workLocal, from <c|h'|c>
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
        if (abs(eigvals(jGlob,iK,iS) - eigvals(iGlob,iK,iS)) < epsilon(0.0_dp)&
            & .and. iGlob /= jGlob) then
          ! degenerate, so no contribution
          workLocal(ii,jj) = 0.0_dp
        else
          workLocal(ii, jj) = workLocal(ii, jj) * &
              & invDiff(eigvals(jGlob,iK,iS),eigvals(iGlob,iK,iS),Ef(iS),tempElec)&
              & * theta(eigvals(jGlob,iK,iS),eigvals(iGlob,iK,iS),tempElec)
        end if
      end do
    end do

    ! Derivatives of states
    call pblasfx_pgemm(eigvecsTransformed, denseDesc%blacsOrbSqr, workLocal, denseDesc%blacsOrbSqr,&
        & dRho, denseDesc%blacsOrbSqr)

    if (allocated(dPsi)) then
      dPsi(:, :, iK, iS, iCart) = workLocal
    end if

    ! Form derivative of occupied density matrix
    call pblasfx_pgemm(dRho, denseDesc%blacsOrbSqr,eigvecsTransformed, denseDesc%blacsOrbSqr,&
        & workLocal, denseDesc%blacsOrbSqr, transb="C", kk=nFilled(iS))
    dRho(:,:) = workLocal
    ! and hermitize
    call pblasfx_ptranc(workLocal, denseDesc%blacsOrbSqr, dRho, denseDesc%blacsOrbSqr,&
        & beta=(1.0_dp,0.0_dp))

  #:else

    ! serial case
    nOrb = size(dRho, dim = 1)

    call unpackHPauli(dHam, kPoint(:,iK), neighbourList%iNeighbour, nNeighbourSK, iSparseStart,&
        & denseDesc%iAtomStart, img2CentCell, iCellVec, cellVec, dRho, idHam)

    ! form |c> H' <c|
    call hemm(workLocal, 'l', dRho, eigVecsCplx(:,:,iKS))
    workLocal(:,:) = matmul(transpose(conjg(eigVecsCplx(:,:,iKS))), workLocal)

    ! orthogonalise degenerate states against perturbation
    call transform%generateUnitary(workLocal, eigvals(:,iK,iS))
    call transform%degenerateTransform(workLocal)

    ! diagonal elements of workLocal are now derivatives of eigenvalues if needed
    if (allocated(dEi)) then
      do ii = 1, nOrb
        dEi(ii, iK, iS, iCart) = real(workLocal(ii,ii),dp)
      end do
    end if

    ! static case

    ! Form actual perturbation U matrix for eigenvectors (potentially at finite T) by
    ! weighting the elements
    do iFilled = 1, nFilled(1)
      do iEmpty = nEmpty(1), nOrb
        if (.not.transform%degenerate(iFilled,iEmpty) .or. iEmpty == iFilled) then
          workLocal(iEmpty, iFilled) = workLocal(iEmpty, iFilled) * &
              & invDiff(eigvals(iFilled, iK, 1), eigvals(iEmpty, iK, 1), Ef(1), tempElec)&
              & *theta(eigvals(iFilled, iK, 1), eigvals(iEmpty, iK, 1), tempElec)
        else
          ! rotation should already have set these elements to zero
          workLocal(iEmpty, iFilled) = 0.0_dp
        end if
      end do
    end do

    eigvecsTransformed = eigVecsCplx(:,:,iS)
    call transform%applyUnitary(eigvecsTransformed)

    ! calculate the derivatives of the eigenvectors
    workLocal(:, :nFilled(1)) =&
        & matmul(eigvecsTransformed(:, nEmpty(1):), workLocal(nEmpty(1):, :nFilled(1)))

    if (allocated(dPsi)) then
      dPsi(:, :, iK, iS, iCart) = workLocal
    end if

    ! zero the uncalculated virtual states
    workLocal(:, nFilled(1)+1:) = 0.0_dp

    ! form the derivative of the density matrix
    dRho(:,:) = matmul(workLocal(:, :nFilled(1)),&
        & transpose(conjg(eigvecsTransformed(:, :nFilled(1)))) )&
        & + matmul(eigvecsTransformed(:, :nFilled(1)),&
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

    end subroutine dRhoStaticPauli


  !> Calculate the change in the density matrix due to shift in the Fermi energy
  subroutine dRhoFermiChangeStaticPauli(dRhoExtra, idRhoExtra, env, parallelKS, iKS, kPoint,&
      & kWeight, iCellVec, cellVec, neighbourList, nNEighbourSK, img2CentCell, iSparseStart, dEfdE,&
      & Ef, nFilled, nEmpty, eigVecsCplx, orb, denseDesc, tempElec, eigVals&
    #:if WITH_SCALAPACK
      &, desc&
    #:endif
      &)

    !> Additional contribution to the density matrix to cancel effect of Fermi energy change
    real(dp), intent(inout) :: dRhoExtra(:,:)

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
    real(dp), intent(in) :: dEfdE(:)

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
            & deltamn(eigVals(jGlob, iK, iS), Ef(iS), tempElec) * dEfdE(iS)
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
          & deltamn(eigvals(iFilled, 1, 1), Ef(1), tempElec) * dEfdE(1)
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

end module dftbp_staticperturb
