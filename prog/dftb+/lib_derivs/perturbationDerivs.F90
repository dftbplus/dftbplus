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

  !> Static (frequency independent) perturbation
  subroutine staticPerturWrtE(env, parallelKS, filling, SSqrReal, eigvals, eigvecs, ham, over, orb,&
      & nAtom, species, speciesnames, neighbourList, nNeighbourSK, denseDesc, iSparseStart,&
      & img2CentCell, coord, sccCalc, maxSccIter, sccTol, nMixElements, nIneqMixElements,&
      & iEqOrbitals, tempElec, Ef, tFixEf, tSpinSharedEf, spinW, thirdOrd, tDFTBU, UJ, nUJ, iUJ,&
      & niUJ, iEqBlockDftbu, onsMEs, iEqBlockOnSite, pChrgMixer, taggedWriter, tWriteAutoTest,&
      & autoTestTagFile, tWriteTaggedOut, taggedResultsFile, tWriteDetailedOut, fdDetailedOut)

    !> Environment settings
    type(TEnvironment), intent(inout) :: env

    !> K-points and spins to process
    type(TParallelKS), intent(in) :: parallelKS

    !> Filling
    real(dp), intent(in) :: filling(:,:,:)

    !> dense real overlap storage
    real(dp), intent(inout) :: SSqrReal(:,:)

    !> Eigenvalue of each level, kpoint and spin channel
    real(dp), intent(in) :: eigvals(:,:,:)

    !> ground state eigenvectors
    real(dp), intent(in) :: eigvecs(:,:,:)

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

    !> Is the Fermi level common accross spin channels?
    logical, intent(in) :: tSpinSharedEf

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


    integer :: iS, iK, iKS, iAt, iCart, iSCC, iLev, iSh, iSp
    integer :: nSpin, nOrbs
    integer, allocatable :: nFilled(:), nEmpty(:)

    integer :: ii, jj, iGlob, jGlob, iEmpty, iFilled
    integer :: iSCCIter
    logical :: tStopSCC

    ! matrices for derivatives of hamiltonian, eigenvalues and eigenvectors
    real(dp), allocatable :: dham(:,:)
    real(dp) :: drho(size(over),size(ham, dim=2))
    real(dp) :: drhoExtra(size(over),size(ham, dim=2))
    real(dp) :: dqIn(orb%mOrb,nAtom,size(ham, dim=2))
    real(dp) :: dqOut(orb%mOrb,nAtom,size(ham, dim=2))
    real(dp) :: dqInpRed(nMixElements), dqOutRed(nMixElements)
    real(dp) :: dqDiffRed(nMixElements), sccErrorQ
    real(dp) :: dqPerShell(orb%mShell,nAtom,size(ham, dim=2))

    real(dp), allocatable :: dqBlockIn(:,:,:,:)
    real(dp), allocatable :: dqBlockOut(:,:,:,:)
    real(dp), allocatable :: dummy(:,:,:,:)

    ! derivative of potentials
    type(TPotentials) :: dpotential

    real(dp), allocatable :: shellPot(:,:,:)

    real(dp) :: dci(size(eigvecs,dim=1),size(eigvecs,dim=2), &
        & size(eigvecs,dim=3))
    real(dp) :: deigvals(size(eigvals,dim=1),size(eigvals,dim=2), &
        & size(eigvals,dim=3))
    real(dp) :: work(size(eigvecs,dim=1),size(eigvecs,dim=2), &
        & size(eigvecs,dim=3))
    real(dp) :: work2(size(eigvecs,dim=1),size(eigvecs,dim=2), &
        & size(eigvecs,dim=3))

    real(dp) :: nF(size(ham, dim=2)), dEf(size(ham, dim=2))

    logical :: tSccCalc, tMetallic, tConverged

    real(dp) :: polarisability(3,3)

    integer :: fdResults

  #:if WITH_SCALAPACK
    ! need distributed matrix descriptors
    integer :: desc(DLEN_), nn

    nn = denseDesc%fullSize
    call scalafx_getdescriptor(env%blacs%orbitalGrid, nn, nn, env%blacs%rowBlockSize,&
        & env%blacs%columnBlockSize, desc)
  #:endif

    if (tFixEf) then
      call error("Perturbation expressions not currently implemented for fixed Fermi energy")
    end if
    if (tSpinSharedEf) then
      call error("Perturbation expressions not currently implemented for shared Fermi energy&
          & between spins")
    end if
    if (allocated(thirdOrd)) then
      call error("Perturbation expressions not currently implemented for 3rd order model")
    end if

    write(stdOut,*)
    write(stdOut,*)'Perturbation calculation of polarisation'
    write(stdOut,*)

    nSpin = size(ham, dim=2)
    nOrbs = size(filling,dim=1)
    allocate(dham(size(ham,dim=1),nSpin))

    tSccCalc = (maxSccIter > 1)

    if (tDFTBU .or. allocated(onsMEs)) then
      allocate(dqBlockIn(orb%mOrb,orb%mOrb,nAtom,nSpin))
      allocate(dqBlockOut(orb%mOrb,orb%mOrb,nAtom,nSpin))
    end if

    call init(dpotential,orb,nAtom,nSpin)

    allocate(nFilled(nSpin))
    allocate(nEmpty(nSpin))

    if (allocated(spinW)) then
      allocate(shellPot(orb%mShell, nAtom, nSpin))
    end if

    nFilled = -1
    do iS = 1, nSpin
      do iLev = 1, nOrbs
        if ( filling(iLev,1,iS) < epsilon(1.0) ) then
          nFilled(iS) = iLev - 1
          exit
        end if
      end do
      if (nFilled(iS) < 0) then
        nFilled(iS) = nOrbs
      end if
    end do
    nEmpty = -1
    do iS = 1, nSpin
      do iLev = 1, nOrbs
        if ( abs( filling(iLev,1,iS) - real(3-nSpin,dp) ) &
            & > epsilon(1.0)) then
          nEmpty(iS) = iLev
          exit
        end if
      end do
      if (nEmpty(iS) < 0) then
        nEmpty(iS) = 1
      end if
    end do

    tMetallic = (.not.all(nFilled == nEmpty -1))

    if (nSpin == 1) then
      write(stdOut,"(1X,A,T40,I0)")'Fully or partly filled states end at :', nFilled
      write(stdOut,"(1X,A,T40,I0)")'Fully or partly empty states start at:', nEmpty
    else
      write(stdOut,"(1X,A,T40,I0,1X,I0)")'Fully or partly filled states end at :', nFilled
      write(stdOut,"(1X,A,T40,I0,1X,I0)")'Fully or partly empty states start at:', nEmpty
    end if
    if (tMetallic) then
      write(stdOut,*)'Metallic system'
    else
      write(stdOut,*)'Non-metallic system'
    end if

    if (tMetallic) then

      ! Density of electrons at the Fermi energy
      do iS = 1, nSpin
        nf(iS) = 0.0_dp
        do ii = nEmpty(iS), nFilled(iS)
          nf(iS) = nf(iS) &
              & + deltamn(Ef(iS),eigvals(ii,1,iS),tempElec)
        end do
        nf(iS) = real(3-nSpin,dp)*nf(iS)
      end do
      write(stdOut,*)'Density of states at the Fermi energy Nf:',nF
    end if

    ! polarisation direction
    ! note, could MPI parallelise over this
    do iCart = 1, 3

      if (tSccCalc) then
        write(stdOut,*)
        write(stdOut,"(1X,A,1X,A,1X,A)")'Calculating',direction(iCart),'direction field'
      end if

      dqIn(:,:,:) = 0.0_dp
      dqOut(:,:,:) = 0.0_dp
      if (tDFTBU .or. allocated(onsMEs)) then
        dqBlockIn(:,:,:,:) = 0.0_dp
        dqBlockOut(:,:,:,:) = 0.0_dp
      end if

      ! derivative of E.x
      dpotential%extAtom(:,:) = 0.0_dp
      do iAt = 1, nAtom
        dpotential%extAtom(iAt,1) = coord(iCart,iAt)
      end do
      dpotential%extShell(:,:,:) = 0.0_dp
      dpotential%extBlock(:,:,:,:) = 0.0_dp
      call total_shift(dpotential%extShell, dpotential%extAtom, orb, species)
      call total_shift(dpotential%extBlock, dpotential%extShell, orb, species)

      if (tSccCalc) then
        call reset(pChrgMixer, nMixElements)
        dqInpRed(:) = 0.0_dp
        dqPerShell(:,:,:) = 0.0_dp
      end if

      if (tSccCalc) then
        write(stdOut,"(1X,A,T12,A)")'SCC Iter','Error'
      end if

      iSCCIter = 1
      tStopSCC = .false.
      lpSCC: do while (iSCCiter <= maxSccIter)

        dpotential%intAtom(:,:) = 0.0_dp
        dpotential%intShell(:,:,:) = 0.0_dp
        dpotential%intBlock(:,:,:,:) = 0.0_dp

        if (tDFTBU .or. allocated(onsMEs)) then
          dpotential%orbitalBlock(:,:,:,:) = 0.0_dp
        end if

        if (tSccCalc .and. iSCCiter>1) then
          call sccCalc%updateCharges(env, dqIn, orb, species)
          call sccCalc%updateShifts(env, orb, species, neighbourList%iNeighbour, img2CentCell)
          call sccCalc%getShiftPerAtom(dpotential%intAtom(:,1))
          call sccCalc%getShiftPerL(dpotential%intShell(:,:,1))

          if (allocated(spinW)) then
            call getChargePerShell(dqIn, orb, species, dqPerShell)
            call getSpinShift(shellPot, dqPerShell, species, orb, spinW)
            dpotential%intShell(:,:,:) = dpotential%intShell + shellPot
          end if

          if (tDFTBU) then
            ! note the derivatives of both FLL and pSIC are pSIC, which is case 2
            call getDftbUShift(dpotential%orbitalBlock, dqBlockIn, species, orb, 2,&
                & UJ, nUJ, niUJ, iUJ)
          end if
          if (allocated(onsMEs)) then
            call addOnsShift(dpotential%orbitalBlock, dpotential%iOrbitalBlock, dqBlockIn,&
                & dummy, onsMEs, species, orb)
          end if


        end if

        call total_shift(dpotential%intShell,dpotential%intAtom, orb,species)
        call total_shift(dpotential%intBlock,dpotential%intShell, orb,species)
        dpotential%intBlock(:,:,:,:) = dpotential%intBlock + dpotential%extBlock
        if (tDFTBU .or. allocated(onsMEs)) then
          dpotential%intBlock(:,:,:,:) = dpotential%intBlock + dpotential%orbitalBlock
        end if

        dham(:,:) = 0.0_dp
        call add_shift(dham, over, nNeighbourSK, neighbourList%iNeighbour, species, orb,&
            & iSparseStart, nAtom, img2CentCell, dpotential%intBlock)

        if (nSpin == 2) then
          dham(:,:) = 2.0_dp * dham(:,:)
        end if
        call qm2ud(dham)

        drho(:,:) = 0.0_dp

        do iKS = 1, parallelKS%nLocalKS
          iK = parallelKS%localKS(1, iKS)
          iS = parallelKS%localKS(2, iKS)

        #:if WITH_SCALAPACK

          call unpackHSRealBlacs(env%blacs, dHam(:,iS), neighbourList%iNeighbour, nNeighbourSK,&
              & iSparseStart, img2CentCell, denseDesc, work(:,:,iKS))

          ! dH times c_i
          call pblasfx_psymm(work(:,:,iKS), denseDesc%blacsOrbSqr, eigvecs(:,:,iKS),&
              & denseDesc%blacsOrbSqr, work2(:,:,iKS), denseDesc%blacsOrbSqr) !, mm=nFilled(iS))

          ! c_i times dH times c_i
          call pblasfx_pgemm(eigvecs(:,:,iKS), denseDesc%blacsOrbSqr, work2(:,:,iKS),&
              & denseDesc%blacsOrbSqr,  work(:,:,iKS), denseDesc%blacsOrbSqr, transa="T")

          ! derivative of eigenvalues stored diagonal of matrix work, from <c|h'|c>
          do jj = 1, size(work,dim=2)
            jGlob = scalafx_indxl2g(jj, desc(NB_), env%blacs%orbitalGrid%mycol, desc(CSRC_),&
                &  env%blacs%orbitalGrid%ncol)
            do ii = 1, size(work,dim=1)
              iGlob = scalafx_indxl2g(ii, desc(MB_), env%blacs%orbitalGrid%myrow, desc(RSRC_),&
                  & env%blacs%orbitalGrid%nrow)

              !if (iGlob == jGlob) then work(ii,jj,iKS) contains a derivative of an eigenvalue

              ! weight with inverse of energy differences
              work(ii,jj,iKS) = work(ii,jj,iKS) * &
                  & invDiff(eigvals(jGlob,1,iS),eigvals(iGlob,1,iS),Ef(iS),tempElec)&
                  & * theta(eigvals(jGlob,1,iS),eigvals(iGlob,1,iS),tempElec)
            end do
          end do

          ! Derivatives of states
          call pblasfx_pgemm(eigvecs(:,:,iKS), denseDesc%blacsOrbSqr, work(:,:,iKS),&
              & denseDesc%blacsOrbSqr, work2(:,:,iKS), denseDesc%blacsOrbSqr)

          ! Form derivative of density matrix
          call pblasfx_pgemm(work2(:,:,iKS), denseDesc%blacsOrbSqr,eigvecs(:,:,iKS),&
              & denseDesc%blacsOrbSqr, work(:,:,iKS), denseDesc%blacsOrbSqr, transb="T",&
              & kk=nFilled(iS))
          work2(:,:,iKS) = work(:,:,iKS)
          ! and symmetrize
          call pblasfx_ptran(work(:,:,iKS), denseDesc%blacsOrbSqr, work2(:,:,iKS),&
              & denseDesc%blacsOrbSqr, beta=1.0_dp)

        #:else

          work2(:,:,iS) = 0.0_dp
          call unpackHS(work2(:,:,iS), dHam(:,iS), neighbourList%iNeighbour, nNeighbourSK,&
              & denseDesc%iAtomStart, iSparseStart, img2CentCell)

          ! form |c> H' <c|
          call symm(work(:,:,iS), 'l', work2(:,:,iS), eigvecs(:,:,iS))
          work(:,:,iS) = matmul(transpose(eigvecs(:,:,iS)), work(:,:,iS))

          ! diagonal elements of work(:,:,iS) are now derivatives of eigenvalues if needed

          ! Form actual perturbation U matrix for eigenvectors (potentially at finite T) by
          ! weighting the elements
          do iFilled = 1, nFilled(iS)
            do iEmpty = nEmpty(iS), nOrbs
              work(iEmpty, iFilled, iS) = work(iEmpty, iFilled, iS) * &
                  & invDiff(eigvals(iFilled, iK, iS), eigvals(iEmpty, iK, iS), Ef(iS), tempElec)&
                  & *theta(eigvals(iFilled, iK, iS), eigvals(iEmpty, iK, iS), tempElec)
            end do
          end do

          ! calculate the derivatives of ci
          work(:, :nFilled(iS), iKS) =&
              & matmul(eigvecs(:, nEmpty(iS):, iS), work(nEmpty(iS):, :nFilled(iS), iS))
          ! zero the uncalculated virtual states
          work(:, nFilled(iS)+1:, iS) = 0.0_dp

          ! form the derivative of the density matrix
          work2(:,:,iS) = matmul(work(:, :nFilled(iS), iS),&
              & transpose(eigvecs(:, :nFilled(iS), iS)))
          work2(:,:,iS) = work2(:,:,iS) +&
              & matmul(eigvecs(:, :nFilled(iS), iS), transpose(work(:, :nFilled(iS), iS)))

        #:endif


        #:if WITH_SCALAPACK
          call packRhoRealBlacs(env%blacs, denseDesc, work2(:,:,iKS), neighbourList%iNeighbour,&
              & nNeighbourSK, orb%mOrb, iSparseStart, img2CentCell, drho(:,iS))
        #:else
          call packHS(drho(:,iS), work2(:,:,iS), neighbourList%iNeighbour, nNeighbourSK,&
              & orb%mOrb, denseDesc%iAtomStart, iSparseStart, img2CentCell)
        #:endif

        end do

      #:if WITH_SCALAPACK
        ! Add up and distribute density matrix contribution from each group
        call mpifx_allreduceip(env%mpi%globalComm, dRho, MPI_SUM)
      #:endif

        drhoextra = 0.0_dp
        if (tMetallic) then
          ! correct for Fermi level shift for q=0 fields

          do iKS = 1, parallelKS%nLocalKS
            iK = parallelKS%localKS(1, iKS)
            iS = parallelKS%localKS(2, iKS)

            dqOut(:,:,iS) = 0.0_dp
            call mulliken(dqOut(:,:,iS), over, drho(:,iS), orb, &
                & neighbourList%iNeighbour, nNeighbourSK, img2CentCell, iSparseStart)

            dEf(iS) = -sum(dqOut(:,:,iS))/nF(iS)

            if (abs(dEf(iS)) > 10.0_dp*epsilon(1.0_dp)) then
              ! Fermi level changes, so need to correct for the change in the number of charges
              work(:,:,iKS) = 0.0_dp

            #:if WITH_SCALAPACK

              do jj = 1, size(work,dim=2)
                jGlob = scalafx_indxl2g(jj,desc(NB_), env%blacs%orbitalGrid%mycol, desc(CSRC_),&
                    & env%blacs%orbitalGrid%ncol)
                if (jGlob >= nEmpty(iS) .and. jGlob <= nFilled(iS)) then
                  work(:,jj,iKS) = eigvecs(:,jj,iKS) * &
                      & deltamn(eigVals(jGlob, 1, iKS), Ef(iS), tempElec) * dEf(iS)
                end if
              end do
              call pblasfx_pgemm(work(:,:,iKS), denseDesc%blacsOrbSqr,eigvecs(:,:,iKS),&
                  & denseDesc%blacsOrbSqr, work2(:,:,iKS), denseDesc%blacsOrbSqr, transb="T",&
                  & alpha=real(3-nSpin,dp))

            #:else

              work(:,:,iS) = 0.0_dp
              do iFilled = nEmpty(iS), nFilled(iS)
                work(:, iFilled, iS) = eigVecs(:, iFilled, iS) * &
                    & deltamn(eigvals(iFilled, iK, iS), Ef(iS), tempElec) * dEf(iS)
              end do
              work2(:, :, iS) = real(3-nSpin,dp) * matmul(work(:, nEmpty(iS):nFilled(iS), iS),&
                  & transpose(eigvecs(:, nEmpty(iS):nFilled(iS), iS)))

            #:endif

              ! pack extra term into density matrix
           #:if WITH_SCALAPACK
              call packRhoRealBlacs(env%blacs, denseDesc, work2(:,:,iKS), neighbourList%iNeighbour,&
                  & nNeighbourSK, orb%mOrb, iSparseStart, img2CentCell, drhoExtra(:,iS))
           #:else
              call packHS(drhoExtra(:,iS), work2(:,:,iS), neighbourList%iNeighbour, nNeighbourSK,&
                  & orb%mOrb, denseDesc%iAtomStart, iSparseStart, img2CentCell)
           #:endif

            end if

          end do

        #:if WITH_SCALAPACK
          ! Add up and distribute density matrix contribution from each group
          call mpifx_allreduceip(env%mpi%globalComm, dRhoExtra, MPI_SUM)
        #:endif
          drho = drho + drhoextra

        end if

        if (nSpin == 1) then
          drho = 2.0_dp * drho
        end if

        call ud2qm(drho)

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
          dqOutRed = 0.0_dp
          call OrbitalEquiv_reduce(dqOut, iEqOrbitals, orb, &
              & dqOutRed(:nIneqMixElements))
          if (tDFTBU) then
            call AppendBlock_reduce(dqBlockOut, iEqBlockDFTBU, orb, dqOutRed )
          end if
          if (allocated(onsMEs)) then
            call onsBlock_reduce(dqBlockOut, iEqBlockOnSite, orb, dqOutRed)
          end if

          dqDiffRed(:) = dqOutRed(:) - dqInpRed(:)
          sccErrorQ = maxval(abs(dqDiffRed))

          write(StdOut,"(1X,I0,T10,E20.12)")iSCCIter, sccErrorQ
          tConverged = (sccErrorQ < sccTol)

          if ((.not. tConverged) .and. iSCCiter /= maxSccIter) then
            if (iSCCIter == 1) then
              dqIn(:,:,:) = dqOut(:,:,:)
              dqInpRed(:) = dqOutRed(:)
              if (tDFTBU .or. allocated(onsMEs)) then
                dqBlockIn(:,:,:,:) = dqBlockOut(:,:,:,:)
              end if
            else

              call mix(pChrgMixer, dqInpRed, dqDiffRed)
            #:if WITH_MPI
              ! Synchronise charges in order to avoid mixers that store a history drifting apart
              call mpifx_allreduceip(env%mpi%globalComm, dqInpRed, MPI_SUM)
              dqInpRed(:) = dqInpRed / env%mpi%globalComm%size
            #:endif

              call OrbitalEquiv_expand(dqInpRed(:nIneqMixElements), iEqOrbitals, &
                  & orb, dqIn)
              if (tDFTBU .or. allocated(onsMEs)) then
                dqBlockIn(:,:,:,:) = 0.0_dp
                if (tDFTBU) then
                  call Block_expand( dqInpRed ,iEqBlockDFTBU, orb, dqBlockIn, species(:nAtom), nUJ,&
                      & niUJ, iUJ, orbEquiv=iEqOrbitals )
                else
                  call Onsblock_expand(dqInpRed, iEqBlockOnSite, orb, dqBlockIn,&
                      & orbEquiv=iEqOrbitals)
                end if
              end if

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
                  & sum(dqIn(orb%posShell(iSh,iSp): &
                  & orb%posShell(iSh+1,iSp)-1,iAt,:nSpin),dim=1)
            end do
          end do

        end if

        iSCCIter = iSCCIter +1

      end do lpSCC

      do ii = 1, 3
        polarisability(ii, iCart) = -sum(sum(dqOut(:,:nAtom,1),dim=1)*coord(ii,:nAtom))
      end do

    end do

    write(stdOut,*)
    write(stdOut,*)'Polarisability'
    do iCart = 1, 3
      write(stdOut,"(3E20.12)")polarisability(:, iCart)
    end do

    if (tWriteAutoTest) then
      open(newunit=fdResults, file=trim(autoTestTagFile), position="append")
      call taggedWriter%write(fdResults, tagLabels%dmudEPerturb, polarisability)
      close(fdResults)
    end if
    if (tWriteTaggedOut) then
      open(newunit=fdResults, file=trim(taggedResultsFile), position="append")
      call taggedWriter%write(fdResults, tagLabels%dmudEPerturb, polarisability)
      close(fdResults)
    end if
    if (tWriteDetailedOut) then
      write(fdDetailedOut,*)'Polarisability (a.u.)'
      do iCart = 1, 3
        write(fdDetailedOut,"(3E20.12)")polarisability(:, iCart)
      end do
    end if
    write(fdDetailedOut,*)


  end subroutine staticPerturWrtE


!  subroutine perturbDyn(mympi, allproc, grpproc, gridAtom, &
!      & groupKS, nGroups, groupsize, desc, filling, nElectrons, SSqrReal, eigvals, &
!      & eigvecs, ham, over, orb, nAtom, species, speciesnames, iNeighbour, &
!      & nNeighbourSK, iAtomStart, iSparseStart, img2CentCell, coord, maxSccIter, sccTol, &
!      & nMixElements, nIneqMixElements, iEqOrbitals, tempElec, Ef, W, pChrgMixer, omega)
!    type(mpifx_comm), intent(in) :: mympi
!    type(blacsgrid), intent(in) :: allproc
!    type(blacsgrid), intent(in) :: grpproc
!    type(blacsgrid), intent(in) :: gridAtom
!    integer, intent(in) :: groupKS(:,:)
!    integer, intent(in) :: nGroups, groupsize
!    integer, intent(in) :: desc(DLEN_)
!    real(dp), intent(in) :: filling(:,:,:)
!    real(dp), intent(in) :: nElectrons(:)
!    real(dp), intent(inout) :: SSqrReal(:,:)
!    real(dp), intent(in) :: eigvals(:,:,:)
!    real(dp), intent(in) :: eigvecs(:,:,:)
!    real(dp), intent(in) :: ham(:,:), over(:)
!    type(TOrbitals), pointer :: orb
!    integer, intent(in) :: nAtom, species(:)
!    character(mc), intent(in) :: speciesnames(:)
!    integer, intent(in) :: iNeighbour(0:,:), nNeighbourSK(:)
!    integer, intent(in) :: denseDesc%iAtomStart(:), iSparseStart(0:,:)
!    integer, intent(in) :: img2CentCell(:)
!    real(dp), intent(in) :: coord(:,:)
!    integer, intent(in) :: maxSccIter
!    real(dp), intent(in) :: sccTol
!    integer, intent(in) :: nMixElements, nIneqMixElements
!    integer, intent(in) :: iEqOrbitals(:,:,:)
!    real(dp), intent(in) :: tempElec, Ef(:)
!    real(dp), intent(in) :: W(:,:,:)
!    type(omixer), pointer, intent(inout) :: pChrgMixer
!    real(dp), intent(in) :: omega(:)
!
!    integer :: iS, iK, iKS, iAt, iCart, iSCC, iLev, iSh, iSp
!    integer :: nSpin, nOrbs, nSparse
!    integer, allocatable :: nFilled(:), nEmpty(:)
!
!    integer :: ii, jj, iGlob, jGlob, iSignOmega, iOmega, nOmega
!    integer :: iSCCIter
!    logical :: tStopSCC
!
!    ! derivatives of hamiltonian, eigenvalues and eigenvectors
!    real(dp), allocatable :: dham(:,:)
!    real(dp) :: drho(size(over),size(ham, dim=2))
!    real(dp) :: drhoExtra(size(over),size(ham, dim=2))
!    real(dp) :: dqIn(orb%mOrb,nAtom,size(ham, dim=2))
!    real(dp) :: dqOut(orb%mOrb,nAtom,size(ham, dim=2))
!    real(dp) :: dqInpRed(nMixElements), dqOutRed(nMixElements)
!    real(dp) :: dqDiffRed(nMixElements), sccErrorQ, polarisability(3,3)
!    real(dp) :: dqPerShell(orb%mOrb,nAtom,size(ham, dim=2))
!    real(dp) :: tmpFill(size(filling,dim=1))
!    type(TPotentials) :: dpotential
!
!    real(dp) :: dci(size(eigvecs,dim=1),size(eigvecs,dim=2), &
!        & size(eigvecs,dim=3))
!    real(dp) :: deigvals(size(eigvals,dim=1),size(eigvals,dim=2), &
!        & size(eigvals,dim=3))
!    real(dp) :: work(size(eigvecs,dim=1),size(eigvecs,dim=2), &
!        & size(eigvecs,dim=3))
!    real(dp) :: work2(size(eigvecs,dim=1),size(eigvecs,dim=2), &
!        & size(eigvecs,dim=3))
!    complex(dp) :: cwork1(size(eigvecs,dim=1),size(eigvecs,dim=2))
!    complex(dp) :: tmpCeigvecs(size(eigvecs,dim=1),size(eigvecs,dim=2))
!    complex(dp) :: cwork2(size(eigvecs,dim=1),size(eigvecs,dim=2))
!
!    real(dp) :: nF(size(ham, dim=2)), dEf(size(ham, dim=2))
!
!    logical :: tSccCalc, tSpin, tMetallic, tConverged
!
!    nSpin = size(ham, dim=2)
!    nSparse = size(ham,dim=1)
!    nOrbs = size(filling,dim=1)
!    nOmega = size(omega)
!
!    tSccCalc = (maxSccIter > 1)
!    tSpin = (nSpin > 1)
!
!    ASSERT(size(eigvecs,dim=3) == nSpin)
!    ASSERT(all(shape(coord) == [3,nAtom]))
!
!      do iCart = 1, 3
!        open(1001+iCart,file="pol_"//direction(iCart),position="rewind", &
!            & status="replace")
!      end do
!
!    call create(dpotential,orb,nAtom,nSpin)
!
!    allocate(nFilled(nSpin))
!    allocate(nEmpty(nSpin))
!
!    nFilled = -1
!    do iS = 1, nSpin
!      do iLev = 1, nOrbs
!        if ( filling(iLev,1,iS) <= epsilon(1.0_dp) ) then
!          nFilled(iS) = iLev - 1
!          exit
!        end if
!      end do
!      if (nFilled(iS) < 0) then
!        nFilled(iS) = nOrbs
!      end if
!    end do
!    nEmpty = -1
!    do iS = 1, nSpin
!      do iLev = 1, nOrbs
!        if ( abs( filling(iLev,1,iS) - real(3-nSpin,dp) ) &
!            & >= epsilon(1.0_dp)) then
!          nEmpty(iS) = iLev
!          exit
!        end if
!      end do
!      if (nEmpty(iS) < 0) then
!        nEmpty(iS) = 1
!      end if
!    end do
!
!    tMetallic = (.not.all(nFilled == nEmpty -1))
!
!      write(stdOut,*)'Fully or partly filled states',nFilled
!      write(stdOut,*)'Fully or partly empty states',nEmpty
!      if (tMetallic) then
!        write(stdOut,*)'Metallic system'
!      else
!        write(stdOut,*)'Non-metallic system'
!      end if
!
!    if (tMetallic) then
!
!      ! Number of electrons at the Fermi energy
!      do iS = 1, nSpin
!        nf(iS) = 0.0_dp
!        do ii = nEmpty(iS), nFilled(iS)
!          nf(iS) = nf(iS) &
!              & + deltamn(Ef(iS),eigvals(ii,1,iS),tempElec)
!        end do
!        nf(iS) = real(3-nSpin,dp)*nf(iS)
!      end do
!        write(stdOut,*)'Electrons at the Fermi energy Nf:',nF
!      end if
!
!    allocate(dham(nSparse,nSpin))
!
!    do iCart = 1, 3 ! polarisation direction
!
!      dqIn = 0.0_dp
!      dqOut = 0.0_dp
!      ! derivative of E.x
!      dpotential%extAtom = 0.0_dp
!      do iAt = 1, nAtom
!        dpotential%extAtom(iAt,1) = coord(iCart,iAt)
!      end do
!
!      dpotential%extBlock = 0.0_dp
!      call total_shift(dpotential%extBlock, dpotential%extAtom, orb, species)
!
!      do iOmega = 1, nOmega
!
!        if (tSccCalc) then
!          call reset(pChrgMixer, nMixElements)
!        end if
!        if (iOmega == 1) then
!          dqInpRed = 0.0_dp
!          dqPerShell = 0.0_dp
!        end if
!
!          write(stdOut,*)'At',omega(iOmega),'au',omega(iOmega)*Hartree__eV,'eV'
!
!        iSCCIter = 1
!        tStopSCC = .false.
!        lpSCC: do while (iSCCiter <= maxSccIter)
!
!          dpotential%intAtom = 0.0_dp
!          dpotential%intBlock = 0.0_dp
!          dpotential%intShell = 0.0_dp
!
!          if (tSccCalc .and. ( iSCCiter>1 .or. iOmega > 1)) then
!            call updateCharges_SCC(allproc, gridAtom, dqIn, orb, species, &
!                & neighbourList%iNeighbour, img2CentCell)
!            call getShiftPerAtom(dpotential%intAtom)
!            call getShiftPerL(dpotential%intShell)
!
!            if (tSpin) then
!              call addSpinShift(dpotential%intShell,dqPerShell,species,orb,W)
!            end if
!
!            !if (tDFTBU) then
!            !  call shift_DFTBU(dpotential%orbitalBlock,qBlockIn,specie,orb, &
!            !      & 2, UJ, nUJ, niUJ, iUJ)
!            !end if
!
!          end if
!
!          call total_shift(dpotential%intShell,dpotential%intAtom, orb,species)
!          call total_shift(dpotential%intBlock,dpotential%intShell, orb,species)
!
!          dpotential%intBlock = dpotential%intBlock + dpotential%extBlock
!
!          dham = 0.0_dp
!          call add_shift(dham, over, nNeighbourSK, neighbourList%iNeighbour, species, orb, iSparseStart, &
!              & nAtom, img2CentCell, dpotential%intBlock)
!
!          if (nSpin == 2) then
!            dham(:,:) = 2.0_dp *  dham(:,:)
!          end if
!          call qm2ud(dham)
!
!          drho = 0.0_dp
!
!          if (grpproc%iproc /= -1) then
!
  !            do iKS = 1, parallelKS%nLocalKS
  !  iK = parallelKS%localKS(1, iKS)
  !  iS = parallelKS%localKS(2, iKS)
!
!              call unpackhs_parallel_real(grpproc, dham(:,iS), neighbourList%iNeighbour, &
!                  & nNeighbourSK, denseDesc%iAtomStart, iSparseStart, img2CentCell, desc, &
!                  & work(:,:,iKS))
!
!              ! dH times c_i
!              call pblasfx_psymm(work(:,:,iKS), desc, eigvecs(:,:,iKS), desc, &
!                  & work2(:,:,iKS), desc) !, mm=nFilled(iS))
!
!              ! c_i times dH times c_i
!              call pblasfx_pgemm(eigvecs(:,:,iKS), desc, work2(:,:,iKS), desc, &
!                  & work(:,:,iKS), desc, transa="T")
!              ! derivative of eigenvalues stored diagonal of matrix work,
!              ! from <c|h'|c>
!
!              tmpCeigvecs = eigvecs(:,:,iKS) ! complex copy of eigenvectors
!
!              do iSignOmega = -1, 1, 2 ! loop over positive and negative frequencies
!
!                do jj = 1, size(work,dim=2)
!                  jGlob = scalafx_indxl2g(jj,desc(NB_),grpproc%mycol,desc(CSRC_),&
!                      & grpproc%ncol)
!                  do ii = 1, size(work,dim=1)
!                    iGlob = scalafx_indxl2g(ii,desc(MB_),grpproc%myrow,desc(RSRC_),&
!                        & grpproc%nrow)
!
!                    ! weight with inverse of energy differences
!                    cwork1(ii,jj) = work(ii,jj,iKS) * &
!                        & ( theta(eigvals(jGlob,1,iS),Ef(iS),tempElec) &
!                        & - theta(eigvals(iGlob,1,iS),Ef(iS),tempElec) ) &
!                        & * theta(eigvals(jGlob,1,iS),eigvals(iGlob,1,iS),tempElec) &
!                        & / (eigvals(jGlob,1,iS)-eigvals(iGlob,1,iS) &
!                        &    + iSignOmega * omega(iOmega) + eta)
!                  end do
!                end do
!
!                ! Derivatives of eigen-state vectors
!                call pblasfx_pgemm(tmpCeigvecs, desc, cwork1, desc, &
!                    & cwork2, desc)
!
!                ! Form derivative of density matrix and then hermitian symmetrise
!                call pblasfx_pgemm(cwork2, desc,tmpCeigvecs, desc, &
!                    & cwork1, desc, transb="T", kk=nFilled(iS))
!
!                cwork2 = cwork1
!
!                call pblasfx_ptranc(cwork2,desc,cwork1,desc, &
!                    & alpha=(0.5_dp,0.0_dp), beta=(0.5_dp,0.0_dp))
!
!                call packrho_parallel_real(grpproc, desc, real(cwork1), &
!                    & neighbourList%iNeighbour, nNeighbourSK, orb%mOrb, denseDesc%iAtomStart, iSparseStart, &
!                    & img2CentCell, drho(:,iS) )
!
!              end do
!
!            end do
!
!          end if
!
!          call blacsfx_gsum(allproc, drho, rdest=-1, cdest=-1)
!
!          if (tMetallic.and.grpproc%iproc /= -1) then
!            ! correct for Fermi level shift for q=0 fields
!            drhoextra = 0.0_dp
!            do iKS = 1, parallelKS%nLocalKS
!  iK = parallelKS%localKS(1, iKS)
!  iS = parallelKS%localKS(2, iKS)
!
!              dqOut = 0.0_dp
!              call mulliken(dqOut(:,:,iS), over, drho(:,iS), orb, &
!                  & neighbourList%iNeighbour, nNeighbourSK, img2CentCell, iSparseStart)
!
!              dEf(iS) = -sum(dqOut(:,:,iS))/nF(iS)
!
!              if (abs(dEf(iS)) > epsilon(1.0_dp)) then ! Fermi level changes,
!                ! so need to correct for change in the number of charges
!                work(:,:,iKS) = 0.0_dp
!                do jj = 1, size(work,dim=2)
!                  jGlob = scalafx_indxl2g(jj,desc(NB_), &
!                      & grpproc%mycol,desc(CSRC_), grpproc%ncol)
!                  if (jGlob >= nEmpty(iS) .and. jGlob <= nFilled(iS)) then
!                    work(:,jj,iKS) = eigvecs(:,jj,iKS) * &
!                        & deltamn(eigVals(jGlob,1,iKS),Ef(iS),tempElec)*dEf(iS)
!                  end if
!                end do
!                call pblasfx_pgemm(work(:,:,iKS), desc,eigvecs(:,:,iKS), desc, &
!                    & work2(:,:,iKS), desc, transb="T",alpha=2.0_dp)
!                call packrho_parallel_real(grpproc, desc, work2(:,:,iKS), &
!                    & neighbourList%iNeighbour, nNeighbourSK, orb%mOrb, denseDesc%iAtomStart, iSparseStart, &
!                    & img2CentCell, drhoextra(:,iS) )
!              end if
!
!            end do
!
!            call blacsfx_gsum(allproc, drhoextra, rdest=-1, cdest=-1)
!            drho = drho + drhoextra
!          end if
!
!          if (nSpin == 1) then
!            drho = 2.0_dp * drho
!          end if
!
!          call ud2qm(drho)
!
!          dqOut = 0.0_dp
!          do iS = 1, nSpin
!            call mulliken(dqOut(:,:,iS), over, drho(:,iS), orb, &
!                & neighbourList%iNeighbour, nNeighbourSK, img2CentCell, iSparseStart)
!            !if (tDFTBU) then
!            !  qBlockOut(:,:,:,iS) = 0.0_dp
!            !  call mulliken(qBlockOut(:,:,:,iS), over, drho(:,iS), &
!            !      &orb, neighbourList%iNeighbour, nNeighbourSK, img2CentCell, iSparseStart)
!            !end if
!          end do
!
!          if (tSccCalc) then
!            dqOutRed = 0.0_dp
!            call OrbitalEquiv_reduce(dqOut, iEqOrbitals, orb, &
!                & dqOutRed(:nIneqMixElements))
!            !if (tDFTBU) then
!            !  call AppendBlock_reduce( qBlockOut, iEqBlockDFTBU, orb, &
!            !      & qOutRed )
!            !end if
!
!            dqDiffRed(:) = dqOutRed(:) - dqInpRed(:)
!            sccErrorQ = maxval(abs(dqDiffRed))
!
!              write(stdOut,*)'Iter',iSCCIter,'Error',sccErrorQ
!            tConverged = (sccErrorQ < sccTol)
!
!            if ((.not. tConverged) .and. iSCCiter /= maxSccIter) then
!              if (iSCCIter == 1 .and. iOmega == 1) then
!                dqIn(:,:,:) = dqOut(:,:,:)
!                dqInpRed(:) = dqOutRed(:)
!                !if (tDFTBU) then
!                !  qBlockIn(:,:,:,:) = qBlockOut(:,:,:,:)
!                !end if
!              else
!                ! Mixing only done on master node as mixer may require IO.
!                if (ioproc) then
!                  call mix(pChrgMixer, dqInpRed, dqDiffRed)
!                end if
!                call mpifx_bcast(mympi, dqInpRed, root=ioprocid)
!
!                call OrbitalEquiv_expand(dqInpRed(:nIneqMixElements), iEqOrbitals, &
!                    & orb, dqIn)
!                !if (tDFTBU) then
!                !  qBlockIn = 0.0_dp
!                !  call Block_expand( qInpRed ,iEqBlockDFTBU, orb, &
!                !      & qBlockIn, specie0, nUJ, niUJ, iUJ, orbEquiv=iEqOrbitals )
!                !end if
!              end if
!            end if
!          else
!            tConverged = .true.
!          end if
!
!          if (tConverged) then
!            exit lpSCC
!          end if
!
!          if (tSpin) then
!            dqPerShell = 0.0_dp
!            do iAt = 1, nAtom
!              iSp = species(iAt)
!              do iSh = 1, orb%nShell(iSp)
!                dqPerShell(iSh,iAt,:nSpin) = dqPerShell(iSh,iAt,:nSpin) +&
!                    & sum(dqIn(orb%posShell(iSh,iSp): &
!                    & orb%posShell(iSh+1,iSp)-1,iAt,:nSpin),dim=1)
!              end do
!            end do
!
!          end if
!
!          iSCCIter = iSCCIter +1
!
!        end do lpSCC
!
!        if (ioproc) then
!          if (iCart==1) then
!            write(stdOut,*)'Polarisability'
!          end if
!          do ii = 1, 3
!            write(stdOut,"(E20.12)",advance ='no') &
!                & -sum(sum(dqOut(:,:nAtom,1),dim=1)*coord(ii,:nAtom))
!            polarisability(ii,iCart) = &
!                & -sum(sum(dqOut(:,:nAtom,1),dim=1)*coord(ii,:nAtom))
!          end do
!          write(stdOut,*)
!          write(1001+icart,"(5E20.12)")omega(iOmega)*Hartree__eV, &
!              & polarisability(:,iCart), sqrt(sum((dqOut(:,:,1)*dqOut(:,:,1))))
!          call flush(1001+icart)
!        end if
!
!      end do
!
!    end do
!
!    if (ioproc) then
!      do ii = 1, 3
!        close(1001+icart)
!      end do
!    end if
!    call destroy(dpotential)
!    deallocate(dham)
!    deallocate(nFilled)
!    deallocate(nEmpty)
!
!  end subroutine perturbDyn


end module dftbp_Perturbderivs
