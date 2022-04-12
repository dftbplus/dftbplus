!--------------------------------------------------------------------------------------------------!
!  DFTB+: general package for performing fast atomistic simulations                                !
!  Copyright (C) 2006 - 2022  DFTB+ developers group                                               !
!                                                                                                  !
!  See the LICENSE file for terms of usage and distribution.                                       !
!--------------------------------------------------------------------------------------------------!

#:include 'common.fypp'

!> Linear response excitations and gradients with respect to atomic coordinates
!>
!> Note: This module is NOT instance safe it uses a common block to communicate with ARPACK
!>
module dftbp_timedep_linrespgrad
  use dftbp_common_accuracy, only : dp, elecTolMax, lc, rsp
  use dftbp_common_constants, only : Hartree__eV, au__Debye, cExchange
  use dftbp_dftb_nonscc, only : TNonSccDiff
  use dftbp_dftb_rangeseparated, only : TRangeSepFunc, getGammaPrimeValue
  use dftbp_dftb_scc, only : TScc
  use dftbp_dftb_shortgammafuncs, only : expGammaPrime
  use dftbp_dftb_sk, only : rotateH0
  use dftbp_dftb_slakocont, only : TSlakoCont, getMIntegrals, getSKIntegrals
  use dftbp_extlibs_arpack, only : withArpack, saupd, seupd
  use dftbp_io_message, only : error
  use dftbp_io_taggedoutput, only : TTaggedWriter, tagLabels
  use dftbp_math_blasroutines, only : gemm, hemv, symm, herk
  use dftbp_math_degeneracy, only : TDegeneracyFind
  use dftbp_math_eigensolver, only : heev
  use dftbp_math_qm, only : makeSimilarityTrans
  use dftbp_math_sorting, only : index_heap_sort, merge_sort
  use dftbp_timedep_linrespcommon, only : excitedDipoleOut, excitedQOut, twothird,&
      & oscillatorStrength, indxoo, indxov, indxvv, rindxov_array, &
      & getSPExcitations, calcTransitionDipoles, dipselect, transitionDipole, writeSPExcitations,&
      & getExcSpin, writeExcMulliken, actionAplusB, actionAminusB, intialSubSpaceMatrixApmB,&
      & calcMatrixSqrt, incMemStratmann, orthonormalizeVectors, getSqrOcc
  use dftbp_timedep_linresptypes, only : TLinResp
  use dftbp_timedep_transcharges, only : TTransCharges, transq, TTransCharges_init
  use dftbp_type_commontypes, only : TOrbitals

  implicit none

  private
  public :: LinRespGrad_old

  !> Output files for results
  character(*), parameter :: transitionsOut = "TRA.DAT"
  character(*), parameter :: XplusYOut = "XplusY.DAT"
  character(*), parameter :: excitedCoefsOut = "COEF.DAT"
  character(*), parameter :: excitationsOut = "EXC.DAT"
  character(*), parameter :: transDipOut = "TDP.DAT"
  character(*), parameter :: naCouplingOut = "NACV.DAT"


  ! Solver related variables

  !> Tolerance for ARPACK solver.
  real(dp), parameter :: ARTOL = epsilon(1.0_rsp)

  !> Threshold for Stratmann solver
  real(dp), parameter :: CONV_THRESH_STRAT = epsilon(1.0_rsp)

  !> Maximal allowed iteration in the ARPACK solver.
  integer, parameter :: MAX_AR_ITER = 300

  !> Names of output files
  character(*), parameter :: arpackOut = "ARPACK.DAT"
  character(*), parameter :: testArpackOut = "TEST_ARPACK.DAT"


contains

  !> This subroutine analytically calculates excitations and gradients of excited state energies
  !> based on Time Dependent DFRT
  subroutine LinRespGrad_old(tSpin, this, iAtomStart, grndEigVecs, grndEigVal, sccCalc, dq, coord0,&
      & SSqr, filling, species0, iNeighbour, img2CentCell, orb, tWriteTagged, fdTagged,&
      & taggedWriter, rangeSep, omega, allOmega, deltaRho, shift, skHamCont, skOverCont, excgrad,&
      & derivator, rhoSqr, occNatural, naturalOrbs)

    !> spin polarized calculation
    logical, intent(in) :: tSpin

    type(TLinResp), intent(inout) :: this

    !> index vector for S and H matrices
    integer, intent(in) :: iAtomStart(:)

    !> ground state MO-coefficients
    real(dp), intent(in) :: grndEigVecs(:,:,:)

    !> ground state MO-energies
    real(dp), intent(in) :: grndEigVal(:,:)

    !> Self-consistent charge module settings
    type(TScc), intent(in) :: sccCalc

    !> converged ground state Mulliken gross charges - atomic charges
    real(dp), intent(in) :: dq(:,:)

    !> atomic positions
    real(dp), intent(in) :: coord0(:,:)

    !> square overlap matrix between basis functions, both triangles required
    real(dp), intent(in) :: SSqr(:,:)

    !> occupations for the states
    real(dp), intent(in) :: filling(:,:)

    !> chemical species of each atom
    integer, intent(in) :: species0(:)

    !> Atomic neighbour lists
    integer, intent(in) :: iNeighbour(0:,:)

    !> Mapping of atom number to central cell atom number
    integer, intent(in) :: img2CentCell(:)

    !> data type for atomic orbital information
    type(TOrbitals), intent(in) :: orb

    !> print tag information
    logical, intent(in) :: tWriteTagged

    !> file descriptor for the tagged data output
    integer, intent(in) :: fdTagged

    !> tagged writer
    type(TTaggedWriter), intent(inout) :: taggedWriter

    !> Data for range-separated calculation
    type(TRangeSepFunc), allocatable, intent(inout) :: rangeSep

    !> excitation energy of state nStat
    real(dp), intent(out) :: omega

    !> excitation energy of all states that have been solved
    real(dp), allocatable, intent(inout) :: allOmega(:)

    !> difference density matrix (vs. uncharged atoms)
    real(dp), intent(inout), pointer :: deltaRho(:,:,:)

    !> shift vector for potentials in the ground state
    real(dp), intent(in), optional :: shift(:)

    !> non-SCC hamiltonian data
    type(TSlakoCont), intent(in), optional :: skHamCont

    !> overlap data
    type(TSlakoCont), intent(in), optional :: skOverCont

    !> excitation energy gradient with respect to atomic positions
    real(dp), intent(out), optional :: excgrad(:,:)

    !> Differentiator for H0 and S matrices.
    class(TNonSccDiff), intent(in), optional :: derivator

    !> ground state density matrix
    real(dp), intent(in), optional :: rhoSqr(:,:,:)

    !> Occupation numbers for natural orbitals from the excited state density matrix
    real(dp), intent(out), optional :: occNatural(:)

    !> the single particle eigenvectors themselves for the excited state density matrix.
    real(dp), intent(out), optional :: naturalOrbs(:,:,:)


    real(dp) :: Ssq(this%nExc), omegaDif, omegaAvg
    real(dp), allocatable :: gammaMat(:,:), lrGamma(:,:), snglPartTransDip(:,:)
    real(dp), allocatable :: ovrXev(:,:,:), wij(:)
    real(dp), allocatable :: dqex(:,:), sposz(:), osz(:), pc(:,:,:)
    real(dp), allocatable :: xpy(:,:), xmy(:,:), sqrOccIA(:), xpym(:),xpyn(:),xmyn(:),xmym(:)
    real(dp), allocatable :: t(:,:,:), rhs(:), woo(:,:), wvv(:,:), wov(:)
    real(dp), allocatable :: rhsm(:), woom(:,:), wvvm(:,:), wovm(:)
    real(dp), allocatable :: eval(:),transitionDipoles(:,:), nacv(:,:)
    integer, allocatable :: win(:), getIA(:,:), getIJ(:,:), getAB(:,:)

    !> array from pairs of single particles states to compound index - should replace with a more
    !> compact data structure in the cases where there are oscilator windows
    integer, allocatable :: iatrans(:,:,:)

    character, allocatable :: symmetries(:)

    integer :: nxoo_max, nxvv_max
    integer, allocatable :: nocc_ud(:), nvir_ud(:)
    integer :: mHOMO, mLUMO
    integer :: nxov, nxov_ud(2), nxov_r, nxov_d, nxov_rd, nxoo_ud(2), nxvv_ud(2)
    integer :: norb, nxoo, nxvv
    integer :: i, j, iSpin, isym, iLev, nStartLev, nEndLev
    integer :: nCoupLev, mCoupLev
    integer :: aa, bb, ss, ab, qq
    integer :: nSpin
    character :: sym
    character(lc) :: tmpStr

    real(dp) :: energyThreshold

    integer :: nStat

    !> control variables
    logical :: tZVector, tFracOcc
    logical :: tRangeSep = .false.
    logical :: tNaCoupling = .true.

    !> should gradients be calculated
    logical :: tForces

    !> transition charges, either cached or evaluated on demand
    type(TTransCharges) :: transChrg

    integer :: fdMulliken, fdTrans, fdTransDip, fdArnoldi, fdXPlusY, fdExc

    !> Communication with ARPACK for progress information
    integer :: logfil, ndigit, mgetv0
    integer :: msaupd, msaup2, msaitr, mseigt, msapps, msgets, mseupd
    integer :: mnaupd, mnaup2, mnaitr, mneigh, mnapps, mngets, mneupd
    integer :: mcaupd, mcaup2, mcaitr, mceigh, mcapps, mcgets, mceupd

    !> Common block of ARPACK variables
    common /debug/ logfil, ndigit, mgetv0,&
        &    msaupd, msaup2, msaitr, mseigt, msapps, msgets, mseupd,&
        &    mnaupd, mnaup2, mnaitr, mneigh, mnapps, mngets, mneupd,&
        &    mcaupd, mcaup2, mcaitr, mceigh, mcapps, mcgets, mceupd

    if (withArpack) then

      ! ARPACK library variables
      ndigit = -3
      ! Output unit:
      logfil = 1
      msgets = 0
      msaitr = 0
      msapps = 0
      mseigt = 0
      mseupd = 0
      if(this%tArnoldi) then
        msaupd = 1
        msaup2 = 1
      else
        msaupd = 0
        msaup2 = 0
      endif
      ! End of ARPACK communication variables

    end if

    if (this%writeMulliken) then
      open(newunit=fdMulliken, file=excitedQOut,position="rewind", status="replace")
      close(fdMulliken)
      open(newunit=fdMulliken, file=excitedDipoleOut, position="rewind", status="replace")
      close(fdMulliken)
    end if

    if (this%tArnoldi) then
      open(newunit=fdArnoldi, file=arpackOut, position="rewind", status="replace")
    else
      fdArnoldi = -1
    end if

    nSpin = size(grndEigVal, dim=2)
    @:ASSERT(nSpin > 0 .and. nSpin <=2)

    norb = orb%nOrb

  #:block DEBUG_CODE
    if (present(excgrad)) then
      @:ASSERT(present(rhoSqr))
      @:ASSERT(present(shift))
      @:ASSERT(present(skHamCont))
      @:ASSERT(present(skOverCont))
      @:ASSERT(present(derivator))
    end if
  #:endblock DEBUG_CODE
    @:ASSERT(present(occNatural) .eqv. present(naturalOrbs))

    ! Should possibly not use allocation status but have a dedicated derived type variable?
    if(allocated(rangeSep)) then
       tRangeSep = .true.
       allocate(lrGamma(this%nAtom, this%nAtom))
       call rangeSep%getLrGamma(lrGamma)
    endif

    ! Try to detect fractional occupations
    tFracOcc = .false.
    do iSpin = 1, nSpin
      do i = 1, norb
        if (filling(i,iSpin) > elecTolMax .and. 2.0_dp/nSpin - filling(i,iSpin) >  elecTolMax) then
          tFracOcc = .true.
          exit
        end if
      end do
    end do
    if (tFracOcc .and. tRangeSep) then
      call error('Fractional occupations not implemented for TD-LC-DFTB.')
    end if

    ! count initial number of transitions from occupied to empty states
    nxov_ud = 0
    do iSpin = 1, nSpin
      !$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(i,j) SCHEDULE(RUNTIME) REDUCTION(+:nxov_ud)
      do i = 1, norb - 1
        do j = i, norb
          if (filling(i,iSpin) > filling(j,iSpin) + elecTolMax) then
            nxov_ud(iSpin) = nxov_ud(iSpin) + 1
          end if
        end do
      end do
      !$OMP  END PARALLEL DO
    end do
    nxov = sum(nxov_ud)

    ! # occupied/virtual states per spin channel
    ALLOCATE(nocc_ud(nSpin))
    ALLOCATE(nvir_ud(nSpin))
    nocc_ud = 0
    nvir_ud = 0
    do iSpin = 1, nSpin
      do i = 1, norb
        if (filling(i,iSpin) > elecTolMax) then
          nocc_ud(iSpin) = nocc_ud(iSpin) + 1
        else
          nvir_ud(iSpin) = nvir_ud(iSpin) + 1
        end if
      end do
    end do

    mHOMO = maxval(nocc_ud)
    mLUMO = minval(nocc_ud) + 1

    ! Dimension of getIJ and getAB
    nxoo_ud(:) = 0
    nxvv_ud(:) = 0
    do iSpin = 1, nSpin
       nxoo_ud(iSpin) = (nocc_ud(iSpin) * (nocc_ud(iSpin) + 1))/2
       nxvv_ud(iSpin) = (nvir_ud(iSpin) * (nvir_ud(iSpin) + 1))/2
    enddo
    nxoo = sum(nxoo_ud)
    nxvv = sum(nxvv_ud)
    nxoo_max = maxval(nxoo_ud)
    nxvv_max = maxval(nxvv_ud)

    if (this%nExc + 1 >= nxov) then
      write(tmpStr,"(' Insufficent single particle excitations, ',I0,&
          & ', for required number of excited states ',I0)")nxov, this%nExc
      call error(tmpStr)
    end if

    tForces = .false.
    ! are gradients required?
    if (present(excgrad)) then
      if (size(excgrad) > 0) then
        tForces = .true.
      end if
    end if

    !> is a z vector required?
    tZVector = tForces .or. this%writeMulliken .or. this%writeCoeffs .or. present(naturalOrbs) .or.&
        & this%tWriteDensityMatrix

    !> occ-occ/vir-vir charges only required for Z-vector/forces or TD-LC-DFTB
    if((.not. tZVector) .and. this%tCacheChargesSame) then
       this%tCacheChargesSame = .false.
    endif

    ! Sanity checks
    nstat = this%nStat
    if (nstat < 0 .and. this%symmetry /= "S") then
      call error("Linresp: Brightest mode only available for singlets.")
    end if
    if (nstat /= 0 .and. this%symmetry == "B") then
      call error("Linresp: Both symmetries not allowed if a specific state is excited")
    end if
    if (tZVector .and. this%nExc > nxov - 1) then
      call error("Linresp: With gradients/properties, nexc can be greater than the number of&
          & occupied-virtual excitations")
    end if

    ! Select symmetries to process
    if (.not. tSpin) then
      select case (this%symmetry)
      case ("B")
        ALLOCATE(symmetries(2))
        symmetries(:) = [ "T", "S" ]
      case ("S")
        ALLOCATE(symmetries(1))
        symmetries(:) = [ "S" ]
      case ("T")
        ALLOCATE(symmetries(1))
        symmetries(:) = [ "T" ]
      end select
    else
      ! ADG: temporary solution for spin polarized case.
      ALLOCATE(symmetries(1))
      symmetries(:) = [ " " ]
    end if
    ! Allocation for general arrays
    ALLOCATE(gammaMat(this%nAtom, this%nAtom))
    ALLOCATE(snglPartTransDip(nxov, 3))
    ALLOCATE(ovrXev(norb, norb, nSpin))
    ALLOCATE(wij(nxov))
    ALLOCATE(win(nxov))
    ALLOCATE(sqrOccIA(nxov))
    ALLOCATE(eval(this%nExc))
    ALLOCATE(getIA(nxov, 3))
    ALLOCATE(getIJ(nxoo, 3))
    ALLOCATE(getAB(nxvv, 3))
    ALLOCATE(transitionDipoles(this%nExc, 3))
    ALLOCATE(sposz(nxov))

    ! Overlap times wave function coefficients - most routines in DFTB+ use lower triangle (would
    ! remove the need to symmetrize the overlap and ground state density matrix in the main code if
    ! this could be used everywhere in these routines)
    do iSpin = 1, nSpin
      call symm(ovrXev(:,:,iSpin), "L", SSqr, grndEigVecs(:,:,iSpin))
    end do

    ! ground state Hubbard U softened coulombic interactions
    call sccCalc%getAtomicGammaMatrix(gammaMat, iNeighbour, img2CentCell)

    ! Oscillator strengths for exited states, when needed.
    ALLOCATE(osz(this%nExc))

    ! Find all single particle transitions and KS energy differences for cases that go from filled
    ! to empty states, create index arrays for ov,oo,vv
    call getSPExcitations(nocc_ud, nvir_ud, grndEigVal, filling, wij, getIA, getIJ, getAB)

    ! put them in ascending energy order
    if (this%tOscillatorWindow) then
      ! use a stable sort so that degenerate transitions from the same single particle state are
      ! grouped together in the results, allowing these to be selected together (since how intensity
      ! is shared out over degenerate transitions is arbitrary between eigensolvers/platforms).
      call merge_sort(win, wij, 1.0E-3_dp*epsilon(0.0_rsp))
    else
      ! do not require stability, use the usual routine to sort, saving an O(N) workspace
      call index_heap_sort(win, wij)
    end if
    wij = wij(win)

    ! Build square root of occupation difference between virtual and occupied states
    call getSqrOcc(filling, win, nxov_ud(1), nxov, getIA, tSpin, sqrOccIA)

    ! dipole strength of transitions between K-S states
    call calcTransitionDipoles(coord0, win, nxov_ud(1), getIA, iAtomStart, ovrXev, grndEigVecs,&
        & snglPartTransDip)

    ! single particle excitation oscillator strengths
    sposz(:) = twothird * wij(:) * sum(snglPartTransDip**2, dim=2)

    if (this%tOscillatorWindow .and. tZVector ) then
      call error("Incompabilitity between excited state property evaluation and an oscillator&
          & strength window at the moment.")
    end if

    if (this%tOscillatorWindow .or. this%tEnergyWindow) then

      if (.not. this%tEnergyWindow) then

        ! find transitions that are strongly dipole allowed (> oscillatorWindow)
        call dipselect(wij, sposz, win, snglPartTransDip, nxov_rd, this%oscillatorWindow,&
            & grndEigVal, getIA)

      else

        ! energy window above the lowest nexc single particle transitions
        energyThreshold = wij(this%nExc) + this%energyWindow
        nxov_r = count(wij <= energyThreshold)

        nxov_d = 0
        if (this%tOscillatorWindow) then

          ! find transitions that are strongly dipole allowed (> oscillatorWindow)
          if (nxov_r < nxov) then
            ! find transitions that are strongly dipole allowed (> oscillatorWindow)
            call dipselect(wij(nxov_r+1:), sposz(nxov_r+1:), win(nxov_r+1:),&
                & snglPartTransDip(nxov_r+1:,:),nxov_d, this%oscillatorWindow,&
                & grndEigVal, getIA)
          end if

        end if

        nxov_rd = nxov_r + nxov_d

      end if
    else

      nxov_rd = nxov

    end if

    if (withArpack) then
      ! just in case energy/dipole windows add no extra states, and is due to an arpack solver
      ! requirement combined with the need to get at least nexc states
      nxov_rd = max(nxov_rd,min(this%nExc+1,nxov))
    else
      nxov_rd = max(nxov_rd,min(this%nExc,nxov))
    end if

    call TTransCharges_init(transChrg, iAtomStart, ovrXev, grndEigVecs, nxov_rd, nxov_ud(1), &
        & nxoo_ud, nxvv_ud, getIA, getIJ, getAB, win, this%tCacheChargesOccVir,             &
        & this%tCacheChargesSame)

    if (this%writeXplusY) then
      open(newunit=fdXPlusY, file=XplusYOut, position="rewind", status="replace")
    else
      fdXPlusY = -1
    end if

    if(this%writeTrans) then
      open(newunit=fdTrans, file=transitionsOut, position="rewind", status="replace")
      write(fdTrans,*)
    else
      fdTrans = -1
    end if

    ! Many-body transition dipole file to excited states
    if (this%writeTransDip) then
      open(newunit=fdTransDip, file=transDipOut, position="rewind", status="replace")
      write(fdTransDip, *)
      write(fdTransDip, '(5x,a,5x,a,2x,a)') "#", 'w [eV]', 'Transition dipole (x,y,z) [Debye]'
      write(fdTransDip, *)
      write(fdTransDip, '(1x,60("="))')
      write(fdTransDip, *)
    else
      fdTransDip = -1
    end if

    ! excitation energies
    open(newunit=fdExc, file=excitationsOut, position="rewind", status="replace")
    write(fdExc,*)
    if (tSpin) then
      write(fdExc,'(5x,a,7x,a,9x,a,9x,a,6x,a,4x,a)') 'w [eV]', 'Osc.Str.', 'Transition',&
          & 'Weight', 'KS [eV]','D<S*S>'
    else
      write(fdExc,'(5x,a,7x,a,9x,a,9x,a,6x,a,4x,a)') 'w [eV]','Osc.Str.', 'Transition',&
          & 'Weight', 'KS [eV]','Sym.'
    end if

    write(fdExc,*)
    write(fdExc,'(1x,80("="))')
    write(fdExc,*)

    ! single particle excitations (output file and tagged file if needed).  Was used for nxov_rd =
    ! size(wij), but now for just states that are actually included in the excitation calculation.
    call writeSPExcitations(wij, win, nxov_ud(1), getIA, this%writeSPTrans, sposz, nxov_rd, tSpin)

    ALLOCATE(xpy(nxov_rd, this%nExc))
    if (tZVector .or. tRangeSep) then
      ALLOCATE(xmy(nxov_rd, this%nExc))
    end if

    ! set up transition indexing
    ALLOCATE(iatrans(norb, norb, nSpin))
    call rindxov_array(win, nxov, nxoo, nxvv, getIA, getIJ, getAB, iatrans)

    if (this%tUseArpack .and. tRangeSep) then
      call error("Range separation requires the Stratmann solver for excitations")
    end if

    do isym = 1, size(symmetries)

      sym = symmetries(isym)
      if (withArpack .and. this%tUseArpack .and. (.not. tRangeSep)) then
        call buildAndDiagExcMatrixArpack(tSpin, wij(:nxov_rd), sym, win, nocc_ud, nvir_ud,&
            & nxoo_ud, nxvv_ud, nxov_ud, nxov_rd, iaTrans, getIA, getIJ, getAB, iAtomStart,&
            & ovrXev, grndEigVecs, filling, sqrOccIA(:nxov_rd), gammaMat, species0, this%spinW,&
            & transChrg, this%testArnoldi, eval, xpy, xmy, this%onSiteMatrixElements, orb,&
            & tRangeSep, tZVector)
      else
        call buildAndDiagExcMatrixStratmann(tSpin, this%subSpaceFactorStratmann, wij(:nxov_rd),&
            & sym, win, nocc_ud, nvir_ud, nxoo_ud, nxvv_ud, nxov_ud, nxov_rd, iaTrans, getIA,&
            & getIJ, getAB, iAtomStart, ovrXev, grndEigVecs, filling, sqrOccIA(:nxov_rd), gammaMat,&
            & species0, this%spinW, transChrg, eval, xpy, xmy, this%onSiteMatrixElements, orb,&
            & tRangeSep, lrGamma, tZVector)
      end if

      ! Excitation oscillator strengths for resulting states
      call getOscillatorStrengths(sym, tSpin, snglPartTransDip(1:nxov_rd,:), eval, xpy, &
            & sqrOccIA(:nxov_rd), nstat, osz, this%writeTransDip, transitionDipoles)

      if (tSpin) then
        call getExcSpin(Ssq, nxov_ud(1), getIA, win, eval, xpy, filling, ovrXev, grndEigVecs)
        call writeExcitations(sym, osz, this%nExc, nxov_ud(1), getIA, win, eval, xpy,&
            & wij(:nxov_rd), fdXPlusY, fdTrans, fdTransDip,&
            & transitionDipoles,  tWriteTagged, fdTagged, taggedWriter, fdExc, Ssq)
      else
        call writeExcitations(sym, osz, this%nExc, nxov_ud(1), getIA, win, eval, xpy,&
            & wij(:nxov_rd), fdXPlusY, fdTrans, fdTransDip,&
            & transitionDipoles, tWriteTagged, fdTagged, taggedWriter, fdExc)
      end if

      if (allocated(allOmega)) then
        if (size(allOmega) /= size(symmetries) * this%nExc) then
          deallocate(allOmega)
        end if
      end if
      if (.not. allocated(allOmega)) then
        allocate(allOmega(size(symmetries) * this%nExc))
      end if
      allOmega(1+(iSym-1)*this%nExc:iSym*this%nExc) = sqrt(eval)

    end do

    if (fdArnoldi /= -1) close(fdArnoldi)
    if (fdTrans /= -1) close(fdTrans)
    if (fdXPlusY /= -1) close(fdXPlusY)
    if (fdTransDip /= -1) close(fdTransDip)
    close(fdExc)

    ! Remove some un-used memory
    deallocate(snglPartTransDip)
    deallocate(transitionDipoles)
    deallocate(sposz)

    if (.not. tZVector) then
      if (nstat == 0) then
        omega = 0.0_dp
      else
        omega = sqrt(eval(nstat))
      end if
    else
      ! calculate Furche vectors and transition density matrix for various properties

      if (nstat == 0) then
        nStartLev = 1
        nEndLev = this%nExc

        if (tForces) then
          call error("Forces currently not available unless a single excited state is specified")
        end if

      else
        nStartLev = nstat
        nEndLev = nstat
      end if

      if (tSpin) then
        if (any( abs(filling) > elecTolMax .and. abs(filling-1.0_dp) > elecTolMax ) ) then
          call error("Fractional fillings not currently possible for excited state property&
              & calculations")
        end if
      else
        if (any( abs(filling) > elecTolMax .and. abs(filling-2.0_dp) > elecTolMax ) ) then
          call error("Fractional fillings not currently possible for excited state property&
              & calculations")
        end if
      end if

      ! Arrays needed for Z vector
      ALLOCATE(t(norb, norb, nSpin))
      ALLOCATE(rhs(nxov_rd))
      ALLOCATE(woo(nxoo_max, nSpin))
      ALLOCATE(wvv(nxvv_max, nSpin))
      ALLOCATE(wov(nxov_rd))

      ! Arrays for gradients and Mulliken analysis
      if (tZVector) then
        ALLOCATE(dqex(this%nAtom, nSpin))
        ALLOCATE(pc(norb, norb, nSpin))
      end if

      do iLev = nStartLev, nEndLev

        omega = sqrt(eval(iLev))

        ! solve for Z and W to get excited state density matrix
        call getZVectorEqRHS(tRangeSep, xpy(:,iLev), xmy(:,iLev), win, iAtomStart, nocc_ud,&
            & nxov_ud(1), transChrg, getIA, getIJ, getAB, iatrans, this%nAtom, species0, &
            & grndEigVal, ovrXev, grndEigVecs, gammaMat, lrGamma, this%spinW, omega, sym, rhs, t,&
            & wov, woo, wvv)

        call solveZVectorPrecond(rhs, tSpin, wij(:nxov_rd), sym, win, nocc_ud, nvir_ud, nxoo_ud, &
            & nxvv_ud, nxov_ud, nxov_rd, iaTrans, getIA, getIJ, getAB, this%nAtom, iAtomStart, &
            & ovrXev, grndEigVecs, filling, sqrOccIA(:nxov_rd), gammaMat, species0, this%spinW, &
            & this%onSiteMatrixElements, orb, transChrg, tRangeSep, lrGamma)

        call calcWVectorZ(rhs, win, nocc_ud, nxov_ud(1), getIA, getIJ, getAB, iaTrans, iAtomStart,&
            & ovrXev, grndEigVecs, gammaMat, grndEigVal, wov, woo, wvv, transChrg, species0, &
            & this%spinW, tRangeSep, lrGamma)

        call calcPMatrix(t, rhs, win, getIA, pc)

        call writeCoeffs(pc, grndEigVecs, filling, this%writeCoeffs, this%tGrndState, occNatural,&
            & naturalOrbs)

        do iSpin = 1, nSpin
          ! Make MO to AO transformation of the excited density matrix
          call makeSimilarityTrans(pc(:,:,iSpin), grndEigVecs(:,:,iSpin))
          call getExcMulliken(iAtomStart, pc(:,:,iSpin), SSqr, dqex(:,iSpin))
        end do

        if (this%tWriteDensityMatrix) then
          call writeDM(iLev, pc, rhoSqr)
        end if

        if (this%writeMulliken) then
          !> for now, only total Mulliken charges
          call writeExcMulliken(sym, iLev, dq(:,1), sum(dqex,dim=2), coord0)
        end if

        if (tForces) then
          call addGradients(sym, nxov_rd, this%nAtom, species0, iAtomStart, norb, nocc_ud,&
              & getIA, getIJ, getAB, win, grndEigVecs, pc, ovrXev, dq, dqex, gammaMat, &
              & lrGamma, this%HubbardU, this%spinW, shift, woo, wov, wvv, transChrg, xpy(:,iLev), &
              & xmy(:,iLev), coord0, orb, skHamCont, skOverCont, derivator, rhoSqr, deltaRho,  &
              & tRangeSep, rangeSep, excgrad)
        end if

        if (tNaCoupling) then
          ! NACV computations require +/- variants of the T,W,Z matrices used for ex. gradients
          ! This overwrites T, RHS and W
          ALLOCATE(rhsm(nxov_rd))
          ALLOCATE(woom(nxoo_max, nSpin))
          ALLOCATE(wvvm(nxvv_max, nSpin))
          ALLOCATE(wovm(nxov_rd))
          ALLOCATE(nacv, mold = excgrad)
          ALLOCATE(xpyn, mold = xpy(:,1))
          ALLOCATE(xpym, mold = xpy(:,1))
          ALLOCATE(xmyn, mold = xpy(:,1))
          ALLOCATE(xmym, mold = xpy(:,1))
          woo = 0
          woom = 0
          wvvm = 0
          wvv = 0

          !!xmy = 0
          pc = 0
          t = 0 
          rhsm = 0
          rhs = 0
          open(67, file='nacv.out')
          do qq = 1,2
             if(qq==1) then
                nCoupLev = 1
                mCoupLev = 2
             else
                nCoupLev = 2
                mCoupLev = 1
             endif

          xpyn = xpy(:,nCoupLev)
          xmyn = xmy(:,nCoupLev)
          xpym = xpy(:,mCoupLev)
          xmym = xmy(:,mCoupLev)

          nacv = 0
          omegaDif = sqrt(eval(nCoupLev)) - sqrt(eval(mCoupLev))
          omegaAvg = 0.5_dp * (sqrt(eval(nCoupLev)) + sqrt(eval(mCoupLev)))
          print *,'TN: omegas', 27.2114* omegaDif,27.2114*omegaAvg,27.2114*sqrt(eval(nCoupLev)),27.2114*sqrt(eval(mCoupLev))

          ! compute + component of RHS for Z-vector eq. in the NaCoupling case
          ! also computes the + components of W and T
          call getNadiaZvectorEqRHS(tRangeSep, xpy(:,nCoupLev), xmy(:,nCoupLev), xpy(:,mCoupLev),  & 
            & xmy(:,mCoupLev), win, iAtomStart, nocc_ud, nxov_ud(1), transChrg, getIA, getIJ,      &
            & getAB, iatrans, this%nAtom, species0, grndEigVal, ovrXev, grndEigVecs, gammaMat,     &
            & lrGamma, this%spinW, omegaDif, omegaAvg, sym, rhs, rhsm, t, wov, wovm, woo,          &
            & woom, wvv, wvvm)

          call solveZVectorPrecond(rhs, tSpin, wij(:nxov_rd), sym, win, nocc_ud, nvir_ud, nxoo_ud, &
            & nxvv_ud, nxov_ud, nxov_rd, iaTrans, getIA, getIJ, getAB, this%nAtom, iAtomStart,     &
            & ovrXev, grndEigVecs, filling, sqrOccIA(:nxov_rd), gammaMat, species0, this%spinW,    &
            & this%onSiteMatrixElements, orb, transChrg, tRangeSep, lrGamma)

          call solveZVectorPrecondMinus(rhsm, tSpin, wij(:nxov_rd), sym, win, nocc_ud, nvir_ud,    &
            & nxoo_ud, nxvv_ud, nxov_ud, nxov_rd, iaTrans, getIA, getIJ, getAB, this%nAtom,        &
            & iAtomStart, ovrXev, grndEigVecs, filling, sqrOccIA(:nxov_rd), gammaMat, species0,    &
            & this%spinW, this%onSiteMatrixElements, orb, transChrg, tRangeSep, lrGamma)

          call calcNadiaWVectorZ(rhs, rhsm, win, nocc_ud, nxov_ud(1), getIA, getIJ, getAB, iaTrans,&
            & iAtomStart, ovrXev, grndEigVecs, gammaMat, grndEigVal, wov, wovm, woo, woom, wvvm,   & 
            & transChrg, species0, this%spinW, tRangeSep, lrGamma, omegaDif)      

          !!rhs = 0 
          !!rhsm = 0

          call calcNadiaPMatrix(t, rhs, rhsm, win, getIA, pc)
!!$          do i = 1,norb
!!$             do j = 1,norb
!!$                write(67,'(2x,i3,2x,i3,f20.16)') i,j,t(i,j,1)
!!$             enddo
!!$          enddo
!!$          write(67,*)
          !!pc = 0
          do iSpin = 1, nSpin
             ! Make MO to AO transformation of the excited density matrix
             call makeSimilarityTrans(pc(:,:,iSpin), grndEigVecs(:,:,iSpin))
             call getExcMulliken(iAtomStart, pc(:,:,iSpin), SSqr, dqex(:,iSpin))
          end do
          print *,dqex,'dqex'
!!$          woo = 0
          woom = 0
          wvvm = 0
!!$          wov = 0
          wovm = 0
!!$          wvv = 0
!!$          xpym = 0
!!$          xmym = 0
!!$          xpyn = 0
!!$          xmyn = 0
!!$          !!pc=0
!!$          nacv = 0

          call addNadiaGradients(sym, nxov_rd, this%nAtom, species0, iAtomStart, norb, nocc_ud,    &
              & getIA, getIJ, getAB, win, grndEigVecs, pc, ovrXev, dq, dqex, gammaMat, lrGamma,    &
              & this%HubbardU, this%spinW, shift, woo, woom, wov, wovm, wvv, wvvm, transChrg,      &
              & xpyn, xmyn, xpym, xmym, coord0, orb,   &   
              & skHamCont, skOverCont, derivator, rhoSqr, deltaRho, tRangeSep, rangeSep, nacv)
!!$         call addNadiaGradients(sym, nxov_rd, this%nAtom, species0, iAtomStart, norb, nocc_ud,    &
!!$              & getIA, getIJ, getAB, win, grndEigVecs, pc, ovrXev, dq, dqex, gammaMat, lrGamma,    &
!!$              & this%HubbardU, this%spinW, shift, woo, woom, wov, wovm, wvv, wvvm, transChrg,      &
!!$              & xpy(:,nCoupLev), xmy(:,nCoupLev), xpy(:,mCoupLev), xmy(:,mCoupLev), coord0, orb,   &   
!!$              & skHamCont, skOverCont, derivator, rhoSqr, deltaRho, tRangeSep, rangeSep, nacv)
          
          do i= 1, size(nacv(1,:))
             write(67,'(3(f16.8,2x))') nacv(1,i), nacv(2,i), nacv(3,i)
          enddo
          write(67,*)
          enddo
          close(67)

        end if
      end do

      if (nstat == 0) then
        omega = 0.0_dp
      end if

    end if

  end subroutine LinRespGrad_old

  !> Solves the RPA equations in their hermitian form (valid for local functionals) at finite T
  !>
  !>  [A  B] X   =    [C  0] X
  !>                w
  !>  [B  A] Y   =    [0 -C] Y
  !>
  !> (see definitions in Marc Casida, in Recent Advances in Density Functional Methods,
  !>  World Scientific, 1995, Part I, p. 155.)
  !>
  !> The hermitian EV problem is given by \Omega F = w^2 F, with
  !>  S = -C (A-B)^{-1} C, \Omega = - S^{-1/2} (A+B) S^{-1/2} and F = (X+Y) * sqrt(w/wia)
  !>
  !> In this routine \Omega is diagonalised by the iterative ARPACK diagonaliser.
  !> The code deals with closed shell systems by diagonalising dedicated singlet/triplet
  !> submatrices.
  !> See Dominguez JCTC 9 4901 (2013)
  !>
  subroutine buildAndDiagExcMatrixArpack(tSpin, wij, sym, win, nocc_ud, nvir_ud,&
      & nxoo_ud, nxvv_ud, nxov_ud, nxov_rd, iaTrans, getIA, getIJ, getAB, iAtomStart, ovrXev,&
      & grndEigVecs, filling, sqrOccIA, gammaMat, species0, spinW, transChrg, testArnoldi,&
      & eval, xpy, xmy, onsMEs, orb, tRangeSep, tZVector)

    !> spin polarisation?
    logical, intent(in) :: tSpin

    !> single particle excitation energies
    real(dp), intent(in) :: wij(:)

    !> symmetry to calculate transitions
    character, intent(in) :: sym

    !> index array for single particle excitations
    integer, intent(in) :: win(:)

    !> occupied orbitals per spin channel
    integer, intent(in) :: nocc_ud(:)

    !> virtual orbitals per spin channel
    integer, intent(in) :: nvir_ud(:)

    !> number of occ-occ transitions per spin channel
    integer, intent(in) :: nxoo_ud(:)

    !> number of vir-vir transitions per spin channel
    integer, intent(in) :: nxvv_ud(:)

    !> number of occ-vir transitions per spin channel
    integer, intent(in) :: nxov_ud(:)

    !> number of occupied-virtual transitions (possibly reduced by windowing)
    integer, intent(in) :: nxov_rd

    !> array from pairs of single particles states to compound index
    integer, intent(in) :: iaTrans(:,:,:)

    !> index array for occ-vir single particle excitations
    integer, intent(in) :: getIA(:,:)

    !> index array for occ-occ single particle excitations
    integer, intent(in) :: getIJ(:,:)

    !> index array for vir-vir single particle excitations
    integer, intent(in) :: getAB(:,:)

    !> indexing array for square matrices
    integer, intent(in) :: iAtomStart(:)

    !> overlap times ground state eigenvectors
    real(dp), intent(in) :: ovrXev(:,:,:)

    !> ground state eigenvectors
    real(dp), intent(in) :: grndEigVecs(:,:,:)

    !> occupation numbers
    real(dp), intent(in) :: filling(:,:)

    ! Square root of occupation difference between vir and occ states
    real(dp), intent(in) :: sqrOccIA(:)

    !> electrostatic matrix
    real(dp), intent(in) :: gammaMat(:,:)

    !> central cell chemical species
    integer, intent(in) :: species0(:)

    !> file handle for ARPACK eigenstate tests
    logical, intent(in) :: testArnoldi

    !> atomic resolved spin constants
    real(dp), intent(in) :: spinW(:)

    !> machinery for transition charges between single particle levels
    type(TTransCharges), intent(in) :: transChrg

    !> resulting eigenvalues for transitions (w^2)
    real(dp), intent(out) :: eval(:)

    !> eigenvectors (X+Y)
    real(dp), intent(out) :: xpy(:,:)

    !> eigenvectors (X-Y), only evaluated if Z-vector is needed
    real(dp), intent(inout), allocatable :: xmy(:,:)

    !> onsite corrections if in use
    real(dp), allocatable :: onsMEs(:,:,:,:)

    !> data type for atomic orbital information
    type(TOrbitals), intent(in) :: orb

    !> is calculation range-separated?
    logical, intent(in) :: tRangeSep

    !> is the Z-vector equation to be solved later?
    logical, intent(in) :: tZVector

    real(dp), allocatable :: workl(:), workd(:), resid(:), vv(:,:), qij(:)
    real(dp) :: sigma, omega
    integer :: iparam(11), ipntr(11), ii
    integer :: ido, ncv, lworkl, info
    logical, allocatable :: selection(:)
    logical :: rvec
    integer :: nexc, natom

    integer :: iState
    real(dp), allocatable :: Hv(:), orthnorm(:,:)
    character(lc) :: tmpStr
    integer :: fdArnoldiTest

    nexc = size(eval)
    natom = size(gammaMat, dim=1)

    @:ASSERT(all(shape(xpy) == [ nxov_rd, nexc ]))
    @:ASSERT(tRangeSep .eqv. .false.)

    ! Three times more Lanczos vectors than desired eigenstates
    ncv = min(3 * nexc, nxov_rd)

    lworkl = ncv * (ncv + 8)

    ALLOCATE(workl(lworkl))
    ALLOCATE(workd(3 * nxov_rd))
    ALLOCATE(resid(nxov_rd))
    ALLOCATE(selection(ncv))
    ALLOCATE(vv(nxov_rd, ncv))
    ALLOCATE(qij(natom))

    resid(:) = 0.0_dp
    workd(:) = 0.0_dp

    ! random initial vector used for dsaupd ARPACK call
    info = 0
    ! IDO must be zero on the first  call
    ido = 0
    ! restarting the iteration with a starting vector that is a linear combination of Ritz vectors
    ! associated with the "wanted" Ritz values.
    iparam(1) = 1
    ! maximum iterations of solver
    iparam(3) = MAX_AR_ITER
    ! solve A*x = lambda*x, with A symmetric
    iparam(7) = 1

    ! loop until exit
    do

      ! call the reverse communication interface from arpack
      call saupd (ido, "I", nxov_rd, "SM", nexc, ARTOL, resid, ncv, vv, nxov_rd, iparam, ipntr,&
          & workd, workl, lworkl, info)

      if (ido == 99) then
        ! has terminated normally, exit loop
        exit
      end if

      ! still running, test for an error return
      if (abs(ido) /= 1) then
        write(tmpStr,"(' Unexpected return from arpack routine saupd, IDO ',I0, ' INFO ',I0)")&
            & ido, info
        call error(tmpStr)
      end if

      ! Action of excitation supermatrix on supervector
      call actionAplusB(tSpin, wij, sym, win, nocc_ud, nvir_ud, nxoo_ud, nxvv_ud, nxov_ud,&
          & nxov_rd, iaTrans, getIA, getIJ, getAB, iAtomStart, ovrXev, grndEigVecs, filling,&
          & sqrOccIA, gammaMat, species0, spinW, onsMEs, orb, .false., transChrg, &
          & workd(ipntr(1):ipntr(1)+nxov_rd-1), workd(ipntr(2):ipntr(2)+nxov_rd-1), tRangeSep)

    end do

    ! check returned info flag for errors
    if (info < 0) then
      write(tmpStr,"(' Error with ARPACK routine saupd, info = ',I0)")info
      call error(tmpStr)
    else if (info  ==  1) then
      call error("Maximum number of iterations reached. Increase the number of excited states to&
          & solve for (NrOfExcitations).")
    else

      ! now want Ritz vectors
      rvec = .true.

      ! everything after the first 6 variables are passed directly to DSEUPD following the last call
      ! to DSAUPD.  These arguments MUST NOT BE MODIFIED between the the last call to DSAUPD and the
      ! call to DSEUPD.
      ! Note: At this point xpy holds the hermitian eigenvectors F
      call seupd (rvec, "All", selection, eval, xpy, nxov_rd, sigma, "I", nxov_rd, "SM", nexc, ARTOL,&
          & resid, ncv, vv, nxov_rd, iparam, ipntr, workd, workl, lworkl, info)

      ! check for error on return
      if (info  /=  0) then
        write(tmpStr,"(' Error with ARPACK routine seupd, info = ',I0)")info
        call error(tmpStr)
      end if

    end if

    if (testArnoldi) then
      ! tests for quality of returned eigenpairs
      open(newunit=fdArnoldiTest, file=testArpackOut, position="rewind", status="replace")
      allocate(Hv(nxov_rd))
      allocate(orthnorm(nxov_rd,nxov_rd))
      orthnorm = matmul(transpose(xpy(:,:nExc)),xpy(:,:nExc))

      write(fdArnoldiTest,"(A)")'State Ei deviation    Evec deviation  Norm deviation  Max&
          & non-orthog'
      do iState = 1, nExc
        call actionAplusB(tSpin, wij, sym, win, nocc_ud, nvir_ud, nxoo_ud, nxvv_ud, nxov_ud,&
            & nxov_rd, iaTrans, getIA, getIJ, getAB, iAtomStart, ovrXev, grndEigVecs, filling,&
            & sqrOccIA, gammaMat, species0, spinW, onsMEs, orb, .false., transChrg, xpy(:,iState),&
            & Hv, .false.)
        write(fdArnoldiTest,"(I4,4E16.8)")iState,&
            & dot_product(Hv,xpy(:,iState))-eval(iState),&
            & sqrt(sum( (Hv-xpy(:,iState)*eval(iState) )**2 )), orthnorm(iState,iState) - 1.0_dp,&
            & max(maxval(orthnorm(:iState-1,iState)), maxval(orthnorm(iState+1:,iState)))
      end do
      close(fdArnoldiTest)
    end if

    if (tZVector) then
      xmy(:,:) = 0.0_dp
    end if

    ! Conversion from eigenvectors of the hermitian problem (F) to (X+Y)
    do iState = 1, nExc
      omega = sqrt(eval(iState))
      xpy(:nxov_rd,iState) = xpy(:nxov_rd,iState) * sqrt(wij(:nxov_rd) / omega)
      if (tZVector) then
        xmy(:nxov_rd,iState) = xpy(:nxov_rd,iState) * omega / wij(:nxov_rd)
      end if
    end do

  end subroutine buildAndDiagExcMatrixArpack

  !> Solves the RPA equations in their standard form at finite T
  !>
  !>  [A  B] X   =    [C  0] X
  !>                w
  !>  [B  A] Y   =    [0 -C] Y
  !>
  !> (see definitions in Marc Casida, in Recent Advances in Density Functional Methods,
  !>  World Scientific, 1995, Part I, p. 155.)
  !>
  !> The RPA eqs are diagonalised by the Stratmann algorithm (JCP 109 8218 (1998).
  !> See also Dominguez JCTC 9 4901 (2013), Kranz JCTC 13 1737 (2017)
  !>
  !> Returns w^2 and (X+Y) (to be consistent with ARPACK diagonaliser)
  !>
  subroutine buildAndDiagExcMatrixStratmann(tSpin, subSpaceFactor, wij, sym, win, nocc_ud, nvir_ud,&
      & nxoo_ud, nxvv_ud, nxov_ud, nxov_rd, iaTrans, getIA, getIJ, getAB, iAtomStart, ovrXev,&
      & grndEigVecs, filling, sqrOccIA, gammaMat, species0, spinW, transChrg, eval, xpy, xmy,&
      & onsMEs, orb, tRangeSep, lrGamma, tZVector)

    !> spin polarisation?
    logical, intent(in) :: tSpin

    !> initial subspace is this factor times number of excited states
    integer :: subSpaceFactor

    !> single particle excitation energies
    real(dp), intent(in) :: wij(:)

    !> symmetry to calculate transitions
    character, intent(in) :: sym

    !> index array for single particle excitations
    integer, intent(in) :: win(:)

    !> occupied orbitals per spin channel
    integer, intent(in) :: nocc_ud(:)

    !> virtual orbitals per spin channel
    integer, intent(in) :: nvir_ud(:)

    !> number of occ-occ transitions per spin channel
    integer, intent(in) :: nxoo_ud(:)

    !> number of vir-vir transitions per spin channel
    integer, intent(in) :: nxvv_ud(:)

    !> number of occ-vir transitions per spin channel
    integer, intent(in) :: nxov_ud(:)

    !> number of occupied-virtual transitions (possibly reduced by windowing)
    integer, intent(in) :: nxov_rd

    !> array from pairs of single particles states to compound index
    integer, intent(in) :: iaTrans(:,:,:)

    !> index array for occ-vir single particle excitations
    integer, intent(in) :: getIA(:,:)

    !> index array for occ-occ single particle excitations
    integer, intent(in) :: getIJ(:,:)

    !> index array for vir-vir single particle excitations
    integer, intent(in) :: getAB(:,:)

    !> indexing array for square matrices
    integer, intent(in) :: iAtomStart(:)

    !> overlap times ground state eigenvectors
    real(dp), intent(in) :: ovrXev(:,:,:)

    !> ground state eigenvectors
    real(dp), intent(in) :: grndEigVecs(:,:,:)

    !> occupation numbers
    real(dp), intent(in) :: filling(:,:)

    ! Square root of occupation difference between vir and occ states
    real(dp), intent(in) :: sqrOccIA(:)

    !> electrostatic matrix
    real(dp), intent(in) :: gammaMat(:,:)

    !> central cell chemical species
    integer, intent(in) :: species0(:)

    !> atomic resolved spin constants
    real(dp), intent(in) :: spinW(:)

    !> machinery for transition charges between single particle levels
    type(TTransCharges), intent(in) :: transChrg

    !> resulting eigenvalues for transitions (actually w^2)
    real(dp), intent(out) :: eval(:)

    !> eigenvectors (X+Y)
    real(dp), intent(out) :: xpy(:,:)

    !> eigenvectors (X-Y), only evaluated if Z-vector is needed
    real(dp), intent(inout), allocatable :: xmy(:,:)

    !> onsite corrections if in use
    real(dp), allocatable :: onsMEs(:,:,:,:)

    !> data type for atomic orbital information
    type(TOrbitals), intent(in) :: orb

    !> is calculation range-separated?
    logical, intent(in) :: tRangeSep

    !> electrostatic matrix, long-range corrected
    real(dp), allocatable, intent(in) :: lrGamma(:,:)

    !> is the Z-vector equation to be solved later?
    logical, intent(in) :: tZVector

    real(dp), allocatable :: vecB(:,:) ! basis of subspace
    real(dp), allocatable :: evecL(:,:), evecR(:,:) ! left and right eigenvectors of Mnh
    real(dp), allocatable :: vP(:,:), vM(:,:) ! vec. for (A+B)b_i, (A-B)b_i
    ! matrices M_plus, M_minus, M_minus^(1/2), M_minus^(-1/2) and M_herm~=resp. mat on subsapce
    real(dp), allocatable :: mP(:,:), mM(:,:), mMsqrt(:,:), mMsqrtInv(:,:), mH(:,:)
    real(dp), allocatable :: evalInt(:) ! store eigenvectors within routine
    real(dp), allocatable :: dummyM(:,:), workArray(:)
    real(dp), allocatable :: vecNorm(:) ! will hold norms of residual vectors
    real(dp) :: dummyReal

    integer :: nExc, nAtom, info, dummyInt, newVec, iterStrat
    integer :: subSpaceDim, memDim, workDim, prevSubSpaceDim
    integer :: ii, jj, ia, ij, ab
    character(lc) :: tmpStr

    logical :: didConverge

    if (allocated(onsMEs)) then
      write(tmpStr,'(A)') 'Onsite corrections not available in Stratmann diagonaliser.'
      call error(tmpStr)
    endif

    ! Number of excited states to solve for
    nExc = size(eval)
    nAtom = size(gammaMat, dim=1)
    @:ASSERT(all(shape(xpy) == [ nxov_rd, nexc ]))

    ! Small subSpaceDim is faster but leads to convergence problems
    ! if large number of excited states is needed
    if (subSpaceFactor < 2) then
      write(tmpStr,'(A)') 'SubSpaceStratmann must be larger than one.'
      call error(tmpStr)
    endif
    subSpaceDim = min(subSpaceFactor * nExc, nxov_rd)
    iterStrat = 1
    write(*,'(A)')
    write(*,'(A)') '>> Stratmann diagonalisation of response matrix'
    write(*,'(3x,A,i6,A,i6)') 'Total dimension of A+B: ', nxov_rd, ' inital subspace: ',&
        & subSpaceDim
    ! Memory available for subspace calcs
    memDim = min(subSpaceDim + 6 * nExc, nxov_rd)
    workDim = 3 * memDim + 1
    allocate(vecB(nxov_rd, memDim))
    allocate(vP(nxov_rd, memDim))
    allocate(vM(nxov_rd, memDim))
    allocate(mP(memDim, memDim))
    allocate(mM(memDim, memDim))
    allocate(mMsqrt(memDim, memDim))
    allocate(mMsqrtInv(memDim, memDim))
    allocate(mH(memDim, memDim))
    allocate(dummyM(memDim, memDim))
    allocate(evalInt(memDim))
    allocate(evecL(memDim, nExc))
    allocate(evecR(memDim, nExc))
    allocate(workArray(3 * memDim + 1))
    allocate(vecNorm(2 * memDim))

    ! set initial bs
    vecB(:,:) = 0.0_dp
    do ii = 1, subSpaceDim
      vecB(ii, ii) = 1.0_dp
    end do

    if (tZVector) then
      xmy(:,:) = 0.0_dp
    end if

    prevSubSpaceDim = 0
    didConverge = .false.

    ! Solve the linear response problem. Iterative expansion of subspace:
    solveLinResp: do

      if (prevSubSpaceDim > 0) then

        ! Extend subspace matrices:
        do ii = prevSubSpaceDim + 1, subSpaceDim
          call actionAplusB(tSpin, wij, sym, win, nocc_ud, nvir_ud, nxoo_ud, nxvv_ud, nxov_ud,&
            & nxov_rd, iaTrans, getIA, getIJ, getAB, iAtomStart, ovrXev, grndEigVecs, filling,&
            & sqrOccIA, gammaMat, species0, spinW, onsMEs, orb, .true., transChrg, vecB(:,ii),&
            & vP(:,ii), tRangeSep, lrGamma)
          call actionAminusB(tSpin, wij, win, nocc_ud, nvir_ud, nxoo_ud, nxvv_ud, nxov_ud, nxov_rd,&
            & iaTrans, getIA, getIJ, getAB, iAtomStart, ovrXev, grndEigVecs, filling, sqrOccIA,&
            & transChrg, vecB(:,ii), vM(:,ii), tRangeSep, lrGamma)
        end do

        do ii = prevSubSpaceDim + 1, subSpaceDim
          do jj = 1, ii
            mP(ii,jj) = dot_product(vecB(:,jj), vP(:,ii))
            mP(jj,ii) = mP(ii,jj)
            mM(ii,jj) = dot_product(vecB(:,jj), vM(:,ii))
            mM(jj,ii) = mM(ii,jj)
          end do
        end do

      else
        ! We need (A+B)_iajb. Could be realized by calls to actionAplusB.
        ! Specific routine for this task is more effective
        call intialSubSpaceMatrixApmB(transChrg, subSpaceDim, wij, sym, win, &
            & nxov_ud(1), iAtomStart, ovrXev, grndEigVecs, filling, sqrOccIA, getIA, getIJ, getAB,&
            & iaTrans, gammaMat, lrGamma, species0, spinW, tSpin, tRangeSep, vP, vM, mP, mM)
      end if

      call calcMatrixSqrt(mM, subSpaceDim, memDim, workArray, workDim, mMsqrt, mMsqrtInv)
      call dsymm('L', 'U', subSpaceDim, subSpaceDim, 1.0_dp, mP, memDim, mMsqrt, memDim,&
          & 0.0_dp, dummyM, memDim)
      call dsymm('L', 'U', subSpaceDim, subSpaceDim, 1.0_dp, mMsqrt, memDim, dummyM, memDim,&
          & 0.0_dp, mH, memDim)

      ! Diagonalise in subspace
      call dsyev('V', 'U', subSpaceDim, mH, memDim, evalInt, workArray, workDim, info)
      if (info /= 0) then
        write(tmpStr,'(A)') 'TDDFT diagonalisation. Increase SubSpaceStratmann.'
        call error(tmpStr)
      endif

      ! This yields T=(A-B)^(-1/2)|X+Y>.
      ! Calc. |R_n>=|X+Y>=(A-B)^(1/2)T and |L_n>=|X-Y>=(A-B)^(-1/2)T.
      ! Transformation preserves orthonormality.
      ! Only compute up to nExc index, because only that much needed.
      call dsymm('L', 'U', subSpaceDim, nExc, 1.0_dp, Mmsqrt, memDim, Mh, memDim, 0.0_dp,&
          & evecR, memDim)
      call dsymm('L', 'U', subSpaceDim, nExc, 1.0_dp, Mmsqrtinv, memDim, Mh, memDim, 0.0_dp,&
          & evecL, memDim)

      ! Need |X-Y>=sqrt(w)(A-B)^(-1/2)T, |X+Y>=(A-B)^(1/2)T/sqrt(w) for proper solution to original
      ! EV problem, only use first nExc vectors
      do ii = 1, nExc
        dummyReal = sqrt(sqrt(evalInt(ii)))
        evecR(:,ii) = evecR(:,ii) / dummyReal
        evecL(:,ii) = evecL(:,ii) * dummyReal
      end do

      !see if more memory is needed to save extended basis. If so increase amount of memory.
      if (subSpaceDim + 2 * nExc > memDim) then
        call incMemStratmann(memDim, workDim, vecB, vP, vM, mP, mM, mH, mMsqrt, mMsqrtInv, &
             &  dummyM, evalInt, workArray, evecL, evecR, vecNorm)
      end if

      ! Calculate the residual vectors
      !   calcs. all |R_n>
      call dgemm('N', 'N', nxov_rd, nExc, subSpaceDim, 1.0_dp, vecB, nxov_rd, evecR, memDim,&
          & 0.0_dp, vecB(1,subSpaceDim+1), nxov_rd)
      !   calcs. all |L_n>
      call dgemm('N', 'N', nxov_rd, nExc, subSpaceDim, 1.0_dp, vecB, nxov_rd, evecL, memDim,&
          & 0.0_dp, vecB(1,subSpaceDim+1+nExc), nxov_rd)

      do ii = 1, nExc
        dummyReal = -sqrt(evalInt(ii))
        vecB(:,subSpaceDim + ii) = dummyReal * vecB(:, subSpaceDim + ii)
        vecB(:,subSpaceDim + nExc + ii) = dummyReal * vecB(:, subSpaceDim + nExc + ii)
      end do

      ! (A-B)|L_n> for all n=1,..,nExc
      call dgemm('N', 'N', nxov_rd, nExc, subSpaceDim, 1.0_dp, vM, nxov_rd, evecL, memDim, 1.0_dp,&
          & vecB(1, subSpaceDim + 1), nxov_rd)
      ! (A+B)|R_n> for all n=1,..,nExc
      call dgemm('N', 'N', nxov_rd, nExc, subSpaceDim, 1.0_dp, vP, nxov_rd, evecR, memDim, 1.0_dp,&
          & vecB(1, subSpaceDim + 1 + nExc), nxov_rd)

      ! calc. norms of residual vectors to check for convergence
      didConverge = .true.
      do ii = subSpaceDim + 1, subSpaceDim + nExc
        vecNorm(ii-subSpaceDim) = dot_product(vecB(:,ii), vecB(:,ii))
        if (vecNorm(ii-subSpaceDim) .gt. CONV_THRESH_STRAT) then
          didConverge = .false.
        end if
      end do

      if (didConverge) then
        do ii = subSpaceDim + nExc + 1, subSpaceDim + 2 * nExc
          vecNorm(ii-subSpaceDim) = dot_product(vecB(:,ii), vecB(:,ii))
          if (vecNorm(ii-subSpaceDim) .gt. CONV_THRESH_STRAT) then
            didConverge = .false.
          end if
        end do
      end if

      if ((.not. didConverge) .and. (subSpaceDim + 2 * nExc > nxov_rd)) then
        write(tmpStr,'(A)') 'Linear Response calculation in subspace did not converge!&
             & Increase SubspaceFactor.'
        call error(tmpStr)
      end if

      ! if converged then exit loop:
      if (didConverge) then
        eval(:) = evalInt(1:nExc)
        ! Calc. X+Y
        xpy(:,:) = matmul(vecB(:,1:subSpaceDim), evecR(1:subSpaceDim,:))
        ! Calc. X-Y, only when needed
        if (tZVector) then
          xmy(:,:) = matmul(vecB(:,1:subSpaceDim), evecL(1:subSpaceDim,:))
        end if
        write(*,'(A)') '>> Stratmann converged'
        exit solveLinResp ! terminate diag. routine
      end if

      ! Otherwise calculate new basis vectors and extend subspace with them
      ! only include new vectors if they add meaningful residue component
      newVec = 0
      do ii = 1, nExc
        if (vecNorm(ii) .gt. CONV_THRESH_STRAT) then
          newVec = newVec + 1
          dummyReal = sqrt(evalInt(ii))
          info = subSpaceDim + ii
          dummyInt = subSpaceDim + newVec
          do jj = 1, nxov_rd
            vecB(jj,dummyInt) = vecB(jj,info) / (dummyReal - wij(jj))
          end do
        end if
      end do

      do ii = 1, nExc
        if (vecNorm(nExc+ii) .gt. CONV_THRESH_STRAT) then
          newVec = newVec + 1
          info = subSpaceDim + nExc + ii
          dummyInt = subSpaceDim + newVec
          do jj = 1, nxov_rd
            vecB(jj,dummyInt) = vecB(jj,info) / (dummyReal - wij(jj))
          end do

        end if
      end do

      prevSubSpaceDim = subSpaceDim
      subSpaceDim = subSpaceDim + newVec
      if(iterStrat == 1) then
         write(*,'(3x,A)') 'Iteration  Subspace dimension'
      end if

      write(*,'(3x,i6,10x,i6)') iterStrat, subSpaceDim
      iterStrat = iterStrat + 1

      ! create orthogonal basis
      call orthonormalizeVectors(prevSubSpaceDim + 1, subSpaceDim, vecB)

    end do solveLinResp

  end subroutine buildAndDiagExcMatrixStratmann


  !> Calculate oscillator strength for a given excitation between KS states
  subroutine getOscillatorStrengths(sym, tSpin, snglPartTransDip, eval, xpy, sqrOccIA, istat, osz, &
     & tTradip, transitionDipoles)

    !> symmetry of transition
    character, intent(in) :: sym

    !> spin polarisation?
    logical, intent(in) :: tSpin

    !> dipole moments for single particle transtions
    real(dp), intent(in) :: snglPartTransDip(:,:)

    !> Low lying eigenvalues of Casida eqn (Omega^2)
    real(dp), intent(in) :: eval(:)

    !> eigenvectors of Casida eqn (X+Y)
    real(dp), intent(in) :: xpy(:,:)

    !> square root of KS occupation differences
    real(dp), intent(in) :: sqrOccIA(:)

    !> write transition dipole
    logical :: tTradip

    !> flag wich if <-1 on entry is returned as the brightest state
    integer, intent(inout) :: istat

    !> Oscilator strengths of transitions
    real(dp), intent(out) :: osz(:)

    !> resulting transition dipoles
    real(dp), intent(out) :: transitionDipoles(:,:)

    integer :: ii, nmat, oszLoc(1)

    nmat = size(xpy, dim=1)

    transitionDipoles(:,:) = 0.0_dp
    osz = 0.0_dp

    ! Triplet oscillator strength and transition dipole is zero for
    ! closed shell ground state
    if ((.not. tSpin) .and. (sym == "T")) then
      return
    end if

    !$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(ii) SCHEDULE(RUNTIME)
    do ii = 1, size(xpy, dim=2)
      osz(ii) = oscillatorStrength(tSpin, snglPartTransDip, sqrt(eval(ii)), xpy(:,ii), sqrOccIA)
    end do
    !$OMP  END PARALLEL DO

    if (istat < 0) then
      ! find largest transition dipole transition
      oszLoc = maxloc(osz)
      istat = oszLoc(1)
    end if

    if (tTradip) then
      call transitionDipole(tSpin, snglPartTransDip, xpy, sqrOccIA, transitionDipoles)
    end if

  end subroutine getOscillatorStrengths


  !> Build right hand side of the equation for the Z-vector and those parts of the W-vectors which
  !> do not depend on Z.
  subroutine getZVectorEqRHS(tRangeSep, xpy, xmy, win, iAtomStart, homo, nmatup, transChrg, getIA, &
      & getIJ, getAB, iatrans, natom, species0, grndEigVal, ovrXev, grndEigVecs, gammaMat, lrGamma,&
      & spinW, omega, sym, rhs, t, wov, woo, wvv)

    !> is calculation range-separated?
    logical, intent(in) :: tRangeSep

    !> X+Y Furche term
    real(dp), intent(in) :: xpy(:)

    !> X-Y Furche term
    real(dp), intent(in) :: xmy(:)

    !> index array for single particle transitions
    integer, intent(in) :: win(:)

    !> index vector for S and H matrices
    integer, intent(in) :: iAtomStart(:)

    !> highest occupied level
    integer, intent(in) :: homo(:)

    !> number of same spin excitations
    integer, intent(in) :: nmatup

    !> machinery for transition charges between single particle levels
    type(TTransCharges), intent(in) :: transChrg

    !> index array between transitions in square and 1D representations
    integer, intent(in) :: getIA(:,:)

    !> index array for vir-vir transitions
    integer, intent(in) :: getIJ(:,:)

    !> index array for occ-occ transitions
    integer, intent(in) :: getAB(:,:)

    !> index array from orbital pairs to compound index
    integer, intent(in) :: iatrans(:,:,:)

    !> number of central cell atoms
    integer, intent(in) :: natom

    !> central cell chemical species
    integer, intent(in) :: species0(:)

    !> ground state wavefunctions
    real(dp), intent(in) :: grndEigVal(:,:)

    !> overlap times ground state wavefunctions
    real(dp), intent(in) :: ovrXev(:,:,:)

    !> ground state wavefunctions
    real(dp), intent(in) :: grndEigVecs(:,:,:)

    !> softened coulomb matrix
    real(dp), intent(in) :: gammaMat(:,:)

    !> softened coulomb matrix, long-range corrected
    real(dp), allocatable, intent(in) :: lrGamma(:,:)

    !> ground state spin derivatives for each species
    real(dp), intent(in) :: spinW(:)

    !> Excitation energies
    real(dp), intent(in) :: omega

    !> Symmetry of the transitions
    character, intent(in) :: sym

    !> Right hand side for the Furche solution
    real(dp), intent(out) :: rhs(:)

    !> T matrix
    real(dp), intent(out) :: t(:,:,:)

    !> W vector occupied-virtual part
    real(dp), intent(out) :: wov(:)

    !> W vector occupied part
    real(dp), intent(out) :: woo(:,:)

    !> W vector virtual part
    real(dp), intent(out) :: wvv(:,:)

    real(dp), allocatable :: xpyq(:), qTr(:), gamxpyq(:), qgamxpyq(:,:), gamqt(:)
    real(dp), allocatable :: xpyqds(:), gamxpyqds(:)
    real(dp), allocatable :: vecHvvXpY(:), vecHvvXmY(:), vecHooXpY(:), vecHooXmY(:)
    real(dp), allocatable :: vecHovT(:), vecHooT(:), vecHvv(:)
    integer :: nxov
    integer, allocatable :: nxoo(:), nxvv(:), nvir(:)
    integer :: i, j, a, b, ias, ibs, abs, ij, ab, jas, ijs, s, nSpin, soo(2), svv(2), nOrb
    real(dp) :: tmp1, tmp2, fact
    logical :: tSpin

    nxov = size(rhs)
    nOrb = size(ovrXev, dim=1)

    ALLOCATE(xpyq(natom))
    ALLOCATE(qTr(natom))
    ALLOCATE(gamxpyq(natom))
    ALLOCATE(gamqt(natom))

    t(:,:,:) = 0.0_dp
    rhs(:) = 0.0_dp
    wov(:) = 0.0_dp
    woo(:,:) = 0.0_dp
    wvv(:,:) = 0.0_dp

    nSpin = size(t, dim=3)

    ALLOCATE(nxoo(nSpin))
    ALLOCATE(nxvv(nSpin))
    ALLOCATE(nvir(nSpin))

    nxoo(:) = (homo(:)*(homo(:)+1))/2
    nvir(:) = size(t, dim=1) - homo(:)
    nxvv(:) = (nvir(:)*(nvir(:)+1))/2

    !! transition charges use compound index ijs = ij + soo(s)
    soo(:) = (/ 0, nxoo(1) /)
    svv(:) = (/ 0, nxvv(1) /)

    ALLOCATE(qgamxpyq(max(maxval(nxoo), maxval(nxvv)), size(homo)))

    if (nSpin == 2) then
      tSpin = .true.
      ALLOCATE(xpyqds(natom))
      ALLOCATE(gamxpyqds(natom))
    else
      tSpin = .false.
    end if

    ! Build t_ab = 0.5 * sum_i (X+Y)_ia (X+Y)_ib + (X-Y)_ia (X-Y)_ib
    ! and w_ab = Q_ab with Q_ab as in (B16) but with corrected sign.
    ! factor 1 / (1 + delta_ab) follows later
    do ias = 1, nxov
      call indxov(win, ias, getIA, i, a, s)

      ! BA: is T_aa = 0?
      do b = homo(s) + 1, a
        ibs = iatrans(i, b, s)
        ab = iaTrans(a, b, s) - svv(s)
        tmp1 = xpy(ias) * xpy(ibs) + xmy(ias) * xmy(ibs)
        tmp2 = omega * (xpy(ias) * xmy(ibs)+ xmy(ias) * xpy(ibs))
        t(a,b,s) = t(a,b,s) + 0.5_dp * tmp1
        ! to prevent double counting
        if (a /= b) then
          t(b,a,s) = t(b,a,s) + 0.5_dp * tmp1
        end if
        ! Note: diagonal elements will be multiplied by 0.5 later.
        wvv(ab,s) = wvv(ab,s) + grndEigVal(i,s) * tmp1 + tmp2
      end do

      ! Build t_ij = 0.5 * sum_a (X+Y)_ia (X+Y)_ja + (X-Y)_ia (X-Y)_ja and 1 / (1 + delta_ij) Q_ij
      ! with Q_ij as in eq. (B9) (1st part of w_ij)
      do j = i, homo(s)
        jas = iatrans(j,a,s)

        ! ADG: assume no constraint on occ space atm (nocc_r = nocc)
        ! otherwise first argument should be nocc - nocc_r
        ij = iatrans(i, j, s) - soo(s)
        tmp1 = xpy(ias) * xpy(jas) + xmy(ias) * xmy(jas)
        tmp2 = omega * (xpy(ias) * xmy(jas) + xmy(ias) * xpy(jas))
        ! Note, there is a typo in Heringer et al. J. Comp Chem 28, 2589.
        ! The sign must be negative see Furche, J. Chem. Phys, 117 7433 (2002).
        t(i,j,s) = t(i,j,s) - 0.5_dp * tmp1
        ! to prevent double counting
        if (i /= j) then
          t(j,i,s) = t(j,i,s) - 0.5_dp * tmp1
        end if
        woo(ij,s) = woo(ij,s) - grndEigVal(a,s) * tmp1 + tmp2
      end do

    end do

    ! xpyq = Q * xpy
    xpyq(:) = 0.0_dp
    call transChrg%qMatVec(iAtomStart, ovrXev, grndEigVecs, getIA, win, xpy, xpyq)

    if (.not. tSpin) then  ! ---- spin-unpolarized case ----
      ! qgamxpyq(ab) = sum_jc K_ab,jc (X+Y)_jc
      if (sym == "S") then
        call hemv(gamxpyq, gammaMat,  xpyq)
        do ab = 1, nxvv(1)
          qTr(:) = transChrg%qTransAB(ab, iAtomStart, ovrXev, grndEigVecs, getAB)
          qgamxpyq(ab, 1) = 2.0_dp * sum(qTr * gamxpyq)
        end do
      else ! triplet case
        do ab = 1, nxvv(1)
          qTr(:) = transChrg%qTransAB(ab, iAtomStart, ovrXev, grndEigVecs, getAB)
          qgamxpyq(ab, 1) = 2.0_dp * sum(qTr * xpyq * spinW(species0))
        end do
      end if

    else  ! ---- spin-polarized case -----

      xpyqds(:) = 0.0_dp
      call transChrg%qMatVecDs(iAtomStart, ovrXev, grndEigVecs, getIA, win, xpy, xpyqds)

      call hemv(gamxpyq, gammaMat,  xpyq)
      do s = 1, 2
        if (s == 1) then
          fact = 1.0_dp
        else
          fact = -1.0_dp
        end if
        do ab = 1, nxvv(s)
          qTr(:) = transChrg%qTransAB(ab + svv(s), iAtomStart, ovrXev, grndEigVecs, getAB)
          qgamxpyq(ab, s) = sum(qTr * gamxpyq)
          !magnetization part
          qgamxpyq(ab, s) = qgamxpyq(ab, s) + fact * sum(qTr * xpyqds * spinW(species0))
        end do
      end do

    end if

    !! Debug TN   
!!$    ALLOCATE(vecHvv(sum(nxvv))) 
!!$    call getHvvXYfr(sym, nXvv, nAtom, iaTrans, getIA, getAB, win,&
!!$      & iAtomStart, species0, ovrXev, grndEigVecs, gammaMat, spinW, transChrg, xpy, vecHvv)
    ! rhs(ia) -= Qia = sum_b (X+Y)_ib * qgamxpyq(ab))
    do ias = 1, nxov
      call indxov(win, ias, getIA, i, a, s)

      do b = homo(s) + 1, a
        ab = iatrans(a, b, s) - svv(s)
        ibs = iatrans(i, b, s)
        !!> TN print *,iatrans(a, b, s),qgamxpyq(ab, s),vecHvv(iatrans(a, b, s)),'Grunfeld'
        rhs(ias) = rhs(ias) - 2.0_dp * xpy(ibs) * qgamxpyq(ab, s)
        ! Since qgamxpyq has only upper triangle
        if (a /= b) then
          rhs(ibs) = rhs(ibs) - 2.0_dp * xpy(ias) * qgamxpyq(ab, s)
        end if
      end do
    end do

    ! -rhs = -rhs - sum_j (X + Y)_ja H + _ij[X + Y]
    if (.not. tSpin) then  ! ---- spin-unpolarized case ----

      if (sym == "S") then
        do ij = 1, nxoo(1)
          qgamxpyq(ij, 1) = 0.0_dp
          qTr(:) = transChrg%qTransIJ(ij, iAtomStart, ovrXev, grndEigVecs, getIJ)
          ! qgamxpyq(ij) = sum_kb K_ij,kb (X+Y)_kb
          qgamxpyq(ij, 1) = 2.0_dp * sum(qTr * gamxpyq)
        end do
      else
        do ij = 1, nxoo(1)
          qgamxpyq(ij, 1) = 0.0_dp
          qTr(:) = transChrg%qTransIJ(ij, iAtomStart, ovrXev, grndEigVecs, getIJ)
          qgamxpyq(ij, 1) = 2.0_dp * sum(qTr * xpyq * spinW(species0))
        end do
      end if

    else  ! ---- spin-polarized case -----

      do s = 1, 2
        if (s == 1) then
          fact = 1.0_dp
        else
          fact = -1.0_dp
        end if
        do ij = 1, nxoo(s)
          qgamxpyq(ij, s) = 0.0_dp
          qTr(:) = transChrg%qTransIJ(ij + soo(s), iAtomStart, ovrXev, grndEigVecs, getIJ)
          qgamxpyq(ij, s) = sum(qTr * gamxpyq)
          !magnetization part
          qgamxpyq(ij, s) = qgamxpyq(ij, s) + fact * sum(qTr * xpyqds * spinW(species0))
        end do
      end do

    end if

    ! rhs(ia) += Qai = sum_j (X+Y)_ja qgamxpyq(ij)
    ! add Qai to Wia as well.
    do ias = 1, nxov
      call indxov(win, ias, getIA, i, a, s)
      do j = i, homo(s)
        jas = iatrans(j, a, s)
        ij = iatrans(i, j, s) - soo(s)
        tmp1 = 2.0_dp * xpy(jas) * qgamxpyq(ij, s)
        rhs(ias) = rhs(ias) + tmp1
        wov(ias) = wov(ias) + tmp1
        if (i /= j) then
          tmp2 = 2.0_dp * xpy(ias) * qgamxpyq(ij, s)
          rhs(jas) = rhs(jas) + tmp2
          wov(jas) = wov(jas) + tmp2
        end if
      end do
    end do

    ! gamxpyq(iAt2) = sum_ij q_ij(iAt2) T_ij
    gamxpyq(:) = 0.0_dp
    if (tSpin) then
      gamxpyqds(:) = 0.0_dp
    end if

    do s = 1, nSpin
      if (s == 1) then
        fact = 1.0_dp
      else
        fact = -1.0_dp
      end if
      do ij = 1, nxoo(s)
        i = getIJ(ij + soo(s), 1)
        j = getIJ(ij + soo(s), 2)
        qTr(:) = transChrg%qTransIJ(ij + soo(s), iAtomStart, ovrXev, grndEigVecs, getIJ)
        if (i == j) then
          gamxpyq(:) = gamxpyq(:) + t(i,j,s) * qTr(:)
          if (tSpin) then
            gamxpyqds(:) = gamxpyqds(:) + t(i,j,s) * qTr(:) * fact
          end if
        else
          ! factor 2 because of symmetry of the matrix
          gamxpyq(:) = gamxpyq(:) + 2.0_dp  * t(i,j,s) * qTr(:)
          if (tSpin) then
            gamxpyqds(:) = gamxpyqds(:) + 2.0_dp * t(i,j,s) * qTr(:) * fact
          end if
        end if
      end do

      ! gamxpyq(iAt2) += sum_ab q_ab(iAt2) T_ab
      do ab = 1, nxvv(s)
        a = getAB(ab + svv(s), 1)
        b = getAB(ab + svv(s), 2)
        qTr(:) = transChrg%qTransAB(ab + svv(s), iAtomStart, ovrXev, grndEigVecs, getAB)
        if (a == b) then
          gamxpyq(:) = gamxpyq(:) + t(a,b,s) * qTr(:)
          if (tSpin) then
            gamxpyqds(:) = gamxpyqds(:) + t(a,b,s) * qTr(:) * fact
          end if
        else
          ! factor 2 because of symmetry of the matrix
          gamxpyq(:) = gamxpyq(:) + 2.0_dp * t(a,b,s) * qTr(:)
          if (tSpin) then
            gamxpyqds(:) = gamxpyqds(:) + 2.0_dp * t(a,b,s) * qTr(:) * fact
          end if
        end if
      end do

    end do

    ! gamqt(iAt1) = sum_iAt2 gamma_iAt1,iAt2 gamxpyq(iAt2)
    call hemv(gamqt, gammaMat, gamxpyq)

    ! rhs -= sum_q^ia(iAt1) gamxpyq(iAt1)
    if (.not. tSpin) then
      call transChrg%qVecMat(iAtomStart, ovrXev, grndEigVecs, getIA, win, -4.0_dp*gamqt, rhs)
    else
      call transChrg%qVecMat(iAtomStart, ovrXev, grndEigVecs, getIA, win, -2.0_dp*gamqt, rhs)
      call transChrg%qVecMatDs(iAtomStart, ovrXev, grndEigVecs, getIA, win, &
           & -2.0_dp*gamxpyqds*spinW(species0), rhs)
    end if

    ! Furche vectors
    do s = 1, nSpin
      if (s == 1) then
        fact = 1.0_dp
      else
        fact = -1.0_dp
      end if
      do ij = 1, nxoo(s)
        qTr(:) = transChrg%qTransIJ(ij + soo(s), iAtomStart, ovrXev, grndEigVecs, getIJ)
        if (.not. tSpin) then
          woo(ij,s) = woo(ij,s) + 4.0_dp * sum(qTr * gamqt)
        else
          woo(ij,s) = woo(ij,s) + 2.0_dp * sum(qTr * gamqt)
          woo(ij,s) = woo(ij,s) + 2.0_dp * fact * sum(qTr * gamxpyqds * spinW(species0))
        end if
      end do
    end do

    ! Contributions due to range-separation
    if (tRangeSep) then

      allocate(vecHvvXpY(sum(nxvv)))
      allocate(vecHvvXmY(sum(nxvv)))
      allocate(vecHooXpY(sum(nxoo)))
      allocate(vecHooXmY(sum(nxoo)))
      allocate(vecHovT(nxov))
      allocate(vecHooT(sum(nxoo)))

      call getHvvXY( 1, nxvv, homo, natom, iatrans, getIA, getAB, win, iAtomStart,&
          & ovrXev, grndEigVecs, lrGamma, transChrg, xpy, vecHvvXpY)

      call getHvvXY(-1, nxvv, homo, natom, iatrans, getIA, getAB, win, iAtomStart,&
          & ovrXev, grndEigVecs, lrGamma, transChrg, xmy, vecHvvXmY)

      call getHooXY( 1, nxoo, homo, natom, iatrans, getIA, getIJ, win, iAtomStart,&
          & ovrXev, grndEigVecs, lrGamma, transChrg, xpy, vecHooXpY)

      call getHooXY(-1, nxoo, homo, natom, iatrans, getIA, getIJ, win, iAtomStart,&
          & ovrXev, grndEigVecs, lrGamma, transChrg, xmy, vecHooXmY)

      call getHovT(nxoo, nxvv, homo, natom, iatrans, getIA, getIJ, getAB, win,&
        & iAtomStart, ovrXev, grndEigVecs, lrGamma, transChrg, t, vecHovT)

      do ias = 1, nxov

        call indXov(win, ias, getIA, i, a, s)
        do b = homo(s) + 1, nOrb
          ibs = iaTrans(i, b, s)
          abs = iaTrans(a, b, s)
          rhs(ias) = rhs(ias) - cExchange * xpy(ibs) * vecHvvXpY(abs)
          if (a >= b) then
            rhs(ias) = rhs(ias) - cExchange * xmy(ibs) * vecHvvXmY(abs)
          else
            rhs(ias) = rhs(ias) + cExchange * xmy(ibs) * vecHvvXmY(abs)
          end if
        end do

        do j = 1, homo(s)
          jas = iaTrans(j, a, s)
          ijs = iaTrans(i, j, s)
          rhs(ias) = rhs(ias) + cExchange * xpy(jas) * vecHooXpY(ijs)
          wov(ias) = wov(ias) + cExchange * xpy(jas) * vecHooXpY(ijs)
          if (i >= j) then
            rhs(ias) = rhs(ias) + cExchange * xmy(jas) * vecHooXmY(ijs)
            wov(ias) = wov(ias) + cExchange * xmy(jas) * vecHooXmY(ijs)
          else
            rhs(ias) = rhs(ias) - cExchange * xmy(jas) * vecHooXmY(ijs)
            wov(ias) = wov(ias) - cExchange * xmy(jas) * vecHooXmY(ijs)
          end if
        end do
       rhs(ias) = rhs(ias) - cExchange * vecHovT(ias)

      end do

      call getHooT(nxov, nxoo, homo, natom, iatrans, getIA, getIJ, win, iAtomStart, &
       & ovrXev, grndEigVecs, lrGamma, transChrg, t, vecHooT)

      !! woo should be made 1D
      do s = 1, nSpin
        do ij = 1, nxoo(s)
          i = getIJ(ij + soo(s), 1)
          j = getIJ(ij + soo(s), 2)
          ijs = iaTrans(i, j, s)
          woo(ij,s) = woo(ij,s) + cExchange * vecHooT(ijs)
        end do
      end do

    endif

  end subroutine getZVectorEqRHS

  !> Build right hand side of the equation for the Z-vector and those parts of the W-vectors which
  !> do not depend on Z. Modified version of getZVectorEqRHS for NA couplings. 
  !> Here the +/- part of the RHS and W is computed. The non-symmetric T is constructed.  
  subroutine getNadiaZvectorEqRHS(tRangeSep, xpyn, xmyn, xpym, xmym, win, iAtomStart, homo,     & 
      & nmatup, transChrg, getIA, getIJ, getAB, iatrans, natom, species0, grndEigVal, ovrXev,&
      & grndEigVecs, gammaMat, lrGamma, spinW, omegaDif, omegaAvg, sym, rhsp, rhsm, t, wovp, &
      & wovm, woop, woom, wvvp, wvvm)

    !> is calculation range-separated?
    logical, intent(in) :: tRangeSep

    !> X+Y Furche term for excited state n 
    real(dp), intent(in) :: xpyn(:)

    !> X-Y Furche term for excited state n 
    real(dp), intent(in) :: xmyn(:)

    !> X+Y Furche term for excited state m 
    real(dp), intent(in) :: xpym(:)

    !> X-Y Furche term for excited state m 
    real(dp), intent(in) :: xmym(:)

    !> index array for single particle transitions
    integer, intent(in) :: win(:)

    !> index vector for S and H matrices
    integer, intent(in) :: iAtomStart(:)

    !> highest occupied level
    integer, intent(in) :: homo(:)

    !> number of same spin excitations
    integer, intent(in) :: nmatup

    !> machinery for transition charges between single particle levels
    type(TTransCharges), intent(in) :: transChrg

    !> index array between transitions in square and 1D representations
    integer, intent(in) :: getIA(:,:)

    !> index array for vir-vir transitions
    integer, intent(in) :: getIJ(:,:)

    !> index array for occ-occ transitions
    integer, intent(in) :: getAB(:,:)

    !> index array from orbital pairs to compound index
    integer, intent(in) :: iatrans(:,:,:)

    !> number of central cell atoms
    integer, intent(in) :: natom

    !> central cell chemical species
    integer, intent(in) :: species0(:)

    !> ground state wavefunctions
    real(dp), intent(in) :: grndEigVal(:,:)

    !> overlap times ground state wavefunctions
    real(dp), intent(in) :: ovrXev(:,:,:)

    !> ground state wavefunctions
    real(dp), intent(in) :: grndEigVecs(:,:,:)

    !> softened coulomb matrix
    real(dp), intent(in) :: gammaMat(:,:)

    !> softened coulomb matrix, long-range corrected
    real(dp), allocatable, intent(in) :: lrGamma(:,:)

    !> ground state spin derivatives for each species
    real(dp), intent(in) :: spinW(:)

    !> Excitation energy difference between states n and m 
    real(dp), intent(in) :: omegaDif

    !> Average excitation energy of states n and m 
    real(dp), intent(in) :: omegaAvg

    !> Symmetry of the transitions
    character, intent(in) :: sym

    !> Right hand side (P+Q)
    real(dp), intent(out) :: rhsp(:)

    !> Right hand side (P-Q)
    real(dp), intent(out) :: rhsm(:)

    !> T matrix (non-symmetric)
    real(dp), intent(out) :: t(:,:,:)

    !> W^+ vector occupied-virtual part
    real(dp), intent(out) :: wovp(:)

    !> W^- vector occupied-virtual part
    real(dp), intent(out) :: wovm(:)

    !> W^+ vector occupied part
    real(dp), intent(out) :: woop(:,:)

    !> W^- vector occupied part
    real(dp), intent(out) :: woom(:,:)

    !> W^+ vector virtual part
    real(dp), intent(out) :: wvvp(:,:)

    !> W^- vector virtual part
    real(dp), intent(out) :: wvvm(:,:)


    real(dp), allocatable :: xpyq(:), qTr(:), gamxpyq(:), qgamxpyq(:,:), gamqt(:)
    real(dp), allocatable :: xpyqds(:), gamxpyqds(:)
    real(dp), allocatable :: vecHvvXorY(:), vecHooXorY(:), vecHovT(:), vecHvvT(:)
    real(dp), allocatable :: vecHvvXpY(:), vecHvvXmY(:), vecHooXpY(:), vecHooXmY(:)
    real(dp), allocatable :: vecHovTp(:), vecHovTm(:), vecHooT(:), tp(:,:,:), tm(:,:,:)
    integer :: nxov
    integer, allocatable :: nxoo(:), nxvv(:), nvir(:)
    integer :: i, j, a, b, ias, ibs, abs, ij, ab, jas, ijs, s, nSpin, soo(2), svv(2), nOrb
    integer :: iOrb, jOrb
    real(dp) :: fact, ptmp1, mtmp1, ptmp2, mtmp2, tmp3, tmp4, tmp1
    logical :: tSpin

    nxov = size(rhsp)
    nOrb = size(ovrXev, dim=1)

    ALLOCATE(xpyq(natom))
    ALLOCATE(qTr(natom))
    ALLOCATE(gamxpyq(natom))
    ALLOCATE(gamqt(natom))


    t(:,:,:) = 0.0_dp
    rhsp(:) = 0.0_dp
    rhsm(:) = 0.0_dp
    wovp(:) = 0.0_dp
    wovm(:) = 0.0_dp
    woop(:,:) = 0.0_dp
    woom(:,:) = 0.0_dp
    wvvp(:,:) = 0.0_dp
    wvvm(:,:) = 0.0_dp

    nSpin = size(t, dim=3)

    ALLOCATE(nxoo(nSpin))
    ALLOCATE(nxvv(nSpin))
    ALLOCATE(nvir(nSpin))

    nxoo(:) = (homo(:)*(homo(:)+1))/2
    nvir(:) = size(t, dim=1) - homo(:)
    nxvv(:) = (nvir(:)*(nvir(:)+1))/2

    !! transition charges use compound index ijs = ij + soo(s)
    soo(:) = (/ 0, nxoo(1) /)
    svv(:) = (/ 0, nxvv(1) /)

    ALLOCATE(qgamxpyq(max(maxval(nxoo), maxval(nxvv)), size(homo)))
    ALLOCATE(vecHooXorY(sum(nxoo)))
    ALLOCATE(vecHvvXorY(sum(nxvv)))   

    if (nSpin == 2) then
      tSpin = .true.
      ALLOCATE(xpyqds(natom))
      ALLOCATE(gamxpyqds(natom))
    else
      tSpin = .false.
    end if

    ! Build t_ab = 0.5 * sum_i (X+Y)_ia (X+Y)_ib + (X-Y)_ia (X-Y)_ib
    ! and w_ab = Q_ab with Q_ab as in (B16) but with corrected sign.
    ! factor 1 / (1 + delta_ab) follows later
    do ias = 1, nxov
      call indxov(win, ias, getIA, i, a, s)

      ! BA: is T_aa = 0?
      do b = homo(s) + 1, a
        ibs = iatrans(i, b, s)
        ab = iaTrans(a, b, s) - svv(s)
        ptmp1 = xpyn(ias) * xpym(ibs) + xmyn(ias) * xmym(ibs) + &
              & xpym(ias) * xpyn(ibs) + xmym(ias) * xmyn(ibs)
   
        mtmp1 = xpym(ias) * xmyn(ibs) + xmym(ias) * xpyn(ibs) - &
              & xpyn(ias) * xmym(ibs) + xmyn(ias) * xpym(ibs)

        ptmp2 = (xpyn(ias) * xmym(ibs) + xmyn(ias) * xpym(ibs) + &
              &  xpym(ias) * xmyn(ibs) + xmym(ias) * xpyn(ibs)) 

        mtmp2 = (xpyn(ias) * xmym(ibs) + xmyn(ias) * xpym(ibs) - &
              &  xpym(ias) * xmyn(ibs) + xmym(ias) * xpyn(ibs)) 

        tmp3 = xpyn(ias) * xpym(ibs) + xmyn(ias) * xmym(ibs)
        tmp4 = xpyn(ibs) * xpym(ias) + xmyn(ibs) * xmym(ias)

        t(a,b,s) = t(a,b,s) + 0.5_dp * tmp3
        wvvp(ab,s) = wvvp(ab,s) + 0.25_dp * grndEigVal(i,s) * ptmp1 / omegaDif  &
                   & + 0.25_dp * omegaAvg * ptmp2 / omegaDif 
                            
        ! to prevent double counting
        if (a /= b) then
          t(b,a,s) = t(b,a,s) + 0.5_dp * tmp4
          wvvm(ab,s) = wvvm(ab,s) + 0.25_dp * grndEigVal(i,s) * mtmp1 / omegaDif  &
                     & + 0.25_dp * omegaAvg * mtmp2 / omegaDif 
        end if
        
      end do

      ! Build t_ij = 0.5 * sum_a (X+Y)_ia (X+Y)_ja + (X-Y)_ia (X-Y)_ja and 1 / (1 + delta_ij) Q_ij
      ! with Q_ij as in eq. (B9) (1st part of w_ij)
      do j = i, homo(s)
        jas = iatrans(j,a,s)

        ! ADG: assume no constraint on occ space atm (nocc_r = nocc)
        ! otherwise first argument should be nocc - nocc_r
        ij = iatrans(i, j, s) - soo(s)

        ptmp1 = (xpyn(ias) * xpym(jas) + xmyn(ias) * xmym(jas) + & 
              &  xpym(ias) * xpyn(jas) + xmym(ias) * xmyn(jas))

        mtmp1 = (xpyn(ias) * xmym(jas) + xmyn(ias) * xpym(jas) - & 
              &  xpym(ias) * xmyn(jas) + xmym(ias) * xpyn(jas))

        ptmp2 = (xpyn(ias) * xmym(jas) + xmyn(ias) * xpym(jas) + &
              &  xpym(ias) * xmyn(jas) + xmym(ias) * xpyn(jas)) 

        mtmp2 = (xpyn(ias) * xpym(jas) + xmyn(ias) * xmym(jas) - &
              &  xpym(ias) * xpyn(jas) + xmym(ias) * xmyn(jas))

        tmp3 = xpyn(ias) * xpym(jas) + xmyn(ias) * xmym(jas)
        tmp4 = xpyn(jas) * xpym(ias) + xmyn(jas) * xmym(ias)

        t(i,j,s) = t(i,j,s) - 0.5_dp * tmp3
        woop(ij,s) = woop(ij,s) - 0.25_dp * grndEigVal(a,s) * ptmp1 / omegaDif  &
                   & + 0.25_dp * omegaAvg * ptmp2 / omegaDif 

        ! to prevent double counting
        if (i /= j) then
          t(j,i,s) = t(j,i,s) - 0.5_dp * tmp4
          woom(ij,s) = woom(ij,s) - 0.25_dp * grndEigVal(a,s) * mtmp1 / omegaDif  &
                   & + 0.25_dp * omegaAvg * mtmp2 / omegaDif 
        end if

      end do

    end do

    ! Terms for (P+-Q) of form (X+Y)^m_ib H^+_ab[(X+Y)^n]  
    call getHplusXYfr(sym, nxoo, nxvv, nAtom, iaTrans, getIA, getIJ, getAB, win, iAtomStart,  &
      & species0, ovrXev, grndEigVecs, gammaMat, spinW, transChrg, xpyn, vecHooXorY, vecHvvXorY)

    do ias = 1, nxov
      call indxov(win, ias, getIA, i, a, s)
      do b = homo(s) + 1, a
        abs = iatrans(a, b, s) 
        ibs = iatrans(i, b, s)
        ! For the forces, we have a factor of 2 here 
        rhsp(ias) = rhsp(ias) - xpym(ibs) * vecHvvXorY(abs)
        rhsm(ias) = rhsm(ias) + xmym(ibs) * vecHvvXorY(abs)
        ! Since vecHvvXpY has only upper triangle
        if (a /= b) then
          rhsp(ibs) = rhsp(ibs) - xpym(ias) * vecHvvXorY(abs)
          rhsm(ibs) = rhsm(ibs) + xmym(ias) * vecHvvXorY(abs)
        end if
      end do

      do j = i, homo(s)
        jas = iatrans(j, a, s)
        ijs = iatrans(i, j, s) 
        ! For the forces, we have a factor of 2 here
        tmp1 = xpym(jas) * vecHooXorY(ijs)
        rhsp(ias) = rhsp(ias) + tmp1
        wovp(ias) = wovp(ias) + 0.5_dp * tmp1 / omegaDif
        tmp1 = xmym(jas) * vecHooXorY(ijs)
        rhsm(ias) = rhsm(ias) - tmp1
        wovm(ias) = wovm(ias) - 0.5_dp * tmp1 / omegaDif       
        if (i /= j) then
           tmp1 = xpym(ias) * vecHooXorY(ijs)
           rhsp(jas) = rhsp(jas) + tmp1
           wovp(jas) = wovp(jas) + 0.5_dp * tmp1 / omegaDif
           tmp1 = xmym(ias) * vecHooXorY(ijs)
           rhsm(jas) = rhsm(jas) - tmp1
           wovm(jas) = wovm(jas) - 0.5_dp * tmp1 / omegaDif
        end if
      end do

    end do

    ! Now m <-> n
    call getHplusXYfr(sym, nxoo, nxvv, nAtom, iaTrans, getIA, getIJ, getAB, win, iAtomStart,  &
      & species0, ovrXev, grndEigVecs, gammaMat, spinW, transChrg, xpym, vecHooXorY, vecHvvXorY)

    do ias = 1, nxov
      call indxov(win, ias, getIA, i, a, s)
      do b = homo(s) + 1, a
        abs = iatrans(a, b, s) 
        ibs = iatrans(i, b, s)
        ! For the forces, we have a factor of 2 here 
        rhsp(ias) = rhsp(ias) - xpyn(ibs) * vecHvvXorY(abs)
        rhsm(ias) = rhsm(ias) - xmyn(ibs) * vecHvvXorY(abs)
        ! Since vecHvvXpY has only upper triangle
        if (a /= b) then
          rhsp(ibs) = rhsp(ibs) - xpyn(ias) * vecHvvXorY(abs)
          rhsm(ibs) = rhsm(ibs) - xmyn(ias) * vecHvvXorY(abs)
        end if
      end do

      do j = i, homo(s)
        jas = iatrans(j, a, s)
        ijs = iatrans(i, j, s) 
        ! For the forces, we have a factor of 2 here
        tmp1 = xpyn(jas) * vecHooXorY(ijs)
        rhsp(ias) = rhsp(ias) + tmp1
        wovp(ias) = wovp(ias) + 0.5_dp * tmp1 / omegaDif
        tmp1 = xmyn(jas) * vecHooXorY(ijs)
        rhsm(ias) = rhsm(ias) + tmp1
        wovm(ias) = wovm(ias) + 0.5_dp * tmp1 / omegaDif       
        if (i /= j) then
           tmp1 = xpyn(ias) * vecHooXorY(ijs)
           rhsp(jas) = rhsp(jas) + tmp1
           wovp(jas) = wovp(jas) + 0.5_dp * tmp1 / omegaDif
           tmp1 = xmyn(ias) * vecHooXorY(ijs)
           rhsm(jas) = rhsm(jas) + tmp1
           wovm(jas) = wovm(jas) + 0.5_dp * tmp1 / omegaDif
        end if
      end do

    end do

    allocate(tp, mold=t)
    allocate(tm, mold=t)
    allocate(vecHovT(nxov))
    allocate(vecHooT(sum(nxoo)))
    allocate(vecHvvT(sum(nxvv)))

    do s = 1, nSpin
      tp(:,:,s) = 0.5_dp * (t(:,:,s) + transpose(t(:,:,s)))
      tm(:,:,s) = 0.5_dp * (t(:,:,s) - transpose(t(:,:,s)))
    end do

    !!> -RHS^+ += - H^+_ia[T^+]
    call getHplusMfr(3, sym, nxoo, nxvv, nxov, nAtom, iaTrans, getIA, getIJ, getAB, win,   &
      & iAtomStart, species0, ovrXev, grndEigVecs, gammaMat, spinW, transChrg,             &
      & tp, vecHovT)

    rhsp = rhsp - vecHovT

    !!> Woo^+ += 0.5 * H^+_ij[T+Z] / Omega_mn, Z part computed later 
    call getHplusMfr(1, sym, nxoo, nxvv, nxov, nAtom, iaTrans, getIA, getIJ, getAB, win,   &
      & iAtomStart, species0, ovrXev, grndEigVecs, gammaMat, spinW, transChrg,             &
      & t, vecHooT)

    do s = 1, nSpin
      do ij = 1, nxoo(s)
        ijs = ij + soo(s)
        woop(ij,s) = woop(ij,s) + 0.5_dp * vecHooT(ijs) / omegaDif
      end do
    end do

    ! Contributions due to range-separation
    if (tRangeSep) then
      allocate(vecHvvXpY(sum(nxvv)))
      allocate(vecHvvXmY(sum(nxvv)))
      allocate(vecHooXpY(sum(nxoo)))
      allocate(vecHooXmY(sum(nxoo)))

      !> Long-range part of H^+[(X+Y)^n] or H^-[(X-Y)^n] for occ-occ and vir-vir comp. of H
      call getHvvXY( 1, nxvv, homo, natom, iatrans, getIA, getAB, win, iAtomStart,&
          & ovrXev, grndEigVecs, lrGamma, transChrg, xpyn, vecHvvXpY)

      call getHvvXY(-1, nxvv, homo, natom, iatrans, getIA, getAB, win, iAtomStart,&
          & ovrXev, grndEigVecs, lrGamma, transChrg, xmyn, vecHvvXmY)

      call getHooXY( 1, nxoo, homo, natom, iatrans, getIA, getIJ, win, iAtomStart,&
          & ovrXev, grndEigVecs, lrGamma, transChrg, xpyn, vecHooXpY)

      call getHooXY(-1, nxoo, homo, natom, iatrans, getIA, getIJ, win, iAtomStart,&
          & ovrXev, grndEigVecs, lrGamma, transChrg, xmyn, vecHooXmY)

      do ias = 1, nxov

        call indXov(win, ias, getIA, i, a, s)
        do b = homo(s) + 1, nOrb
          ibs = iaTrans(i, b, s)
          abs = iaTrans(a, b, s)
          rhsp(ias) = rhsp(ias) - cExchange * 0.5_dp * xpym(ibs) * vecHvvXpY(abs)
          rhsm(ias) = rhsm(ias) + cExchange * 0.5_dp * xmym(ibs) * vecHvvXpY(abs)
          if (a >= b) then
            rhsp(ias) = rhsp(ias) - cExchange * 0.5_dp * xmym(ibs) * vecHvvXmY(abs)
            rhsm(ias) = rhsm(ias) + cExchange * 0.5_dp * xpym(ibs) * vecHvvXmY(abs)
          !> Only a>b is stored in vecHvvXmY, which is anti-symmetric
          else
            rhsp(ias) = rhsp(ias) + cExchange * 0.5_dp * xmym(ibs) * vecHvvXmY(abs)
            rhsm(ias) = rhsm(ias) - cExchange * 0.5_dp * xpym(ibs) * vecHvvXmY(abs)
          end if
        end do

        do j = 1, homo(s)
          jas = iaTrans(j, a, s)
          ijs = iaTrans(i, j, s)
          rhsp(ias) = rhsp(ias) + cExchange * 0.50_dp * xpym(jas) * vecHooXpY(ijs)
          wovp(ias) = wovp(ias) + cExchange * 0.25_dp * xpym(jas) * vecHooXpY(ijs) / omegaDif
          wovm(ias) = wovm(ias) - cExchange * 0.25_dp * xmym(jas) * vecHooXpY(ijs) / omegaDif
          if (i >= j) then
            rhsp(ias) = rhsp(ias) + cExchange * 0.50_dp * xmym(jas) * vecHooXmY(ijs)
            wovp(ias) = wovp(ias) + cExchange * 0.25_dp * xmym(jas) * vecHooXmY(ijs) / omegaDif
            wovm(ias) = wovm(ias) - cExchange * 0.25_dp * xpym(jas) * vecHooXmY(ijs) / omegaDif
          else
            rhsp(ias) = rhsp(ias) - cExchange * 0.50_dp * xmym(jas) * vecHooXmY(ijs)
            wovp(ias) = wovp(ias) - cExchange * 0.25_dp * xmym(jas) * vecHooXmY(ijs) / omegaDif
            wovm(ias) = wovm(ias) + cExchange * 0.25_dp * xpym(jas) * vecHooXmY(ijs) / omegaDif
          end if
        end do

      end do

      !!> Now n <-> m 

      !> Long-range part of H^+[(X+Y)^n] or H^-[(X-Y)^n] for occ-occ and vir-vir comp. of H
      call getHvvXY( 1, nxvv, homo, natom, iatrans, getIA, getAB, win, iAtomStart,&
          & ovrXev, grndEigVecs, lrGamma, transChrg, xpym, vecHvvXpY)

      call getHvvXY(-1, nxvv, homo, natom, iatrans, getIA, getAB, win, iAtomStart,&
          & ovrXev, grndEigVecs, lrGamma, transChrg, xmym, vecHvvXmY)

      call getHooXY( 1, nxoo, homo, natom, iatrans, getIA, getIJ, win, iAtomStart,&
          & ovrXev, grndEigVecs, lrGamma, transChrg, xpym, vecHooXpY)

      call getHooXY(-1, nxoo, homo, natom, iatrans, getIA, getIJ, win, iAtomStart,&
          & ovrXev, grndEigVecs, lrGamma, transChrg, xmym, vecHooXmY)

      do ias = 1, nxov

        call indXov(win, ias, getIA, i, a, s)
        do b = homo(s) + 1, nOrb
          ibs = iaTrans(i, b, s)
          abs = iaTrans(a, b, s)
          rhsp(ias) = rhsp(ias) - cExchange * 0.5_dp * xpyn(ibs) * vecHvvXpY(abs)
          rhsm(ias) = rhsm(ias) - cExchange * 0.5_dp * xmyn(ibs) * vecHvvXpY(abs)
          if (a >= b) then
            rhsp(ias) = rhsp(ias) - cExchange * 0.5_dp * xmyn(ibs) * vecHvvXmY(abs)
            rhsm(ias) = rhsm(ias) - cExchange * 0.5_dp * xpyn(ibs) * vecHvvXmY(abs)
          !> Only a>b is stored in vecHvvXmY, which is anti-symmetric
          else
            rhsp(ias) = rhsp(ias) + cExchange * 0.5_dp * xmyn(ibs) * vecHvvXmY(abs)
            rhsm(ias) = rhsm(ias) + cExchange * 0.5_dp * xpyn(ibs) * vecHvvXmY(abs)
          end if
        end do

        do j = 1, homo(s)
          jas = iaTrans(j, a, s)
          ijs = iaTrans(i, j, s)
          rhsp(ias) = rhsp(ias) + cExchange * 0.50_dp * xpyn(jas) * vecHooXpY(ijs)
          wovp(ias) = wovp(ias) + cExchange * 0.25_dp * xpyn(jas) * vecHooXpY(ijs) / omegaDif
          wovm(ias) = wovm(ias) + cExchange * 0.25_dp * xmyn(jas) * vecHooXpY(ijs) / omegaDif
          if (i >= j) then
            rhsp(ias) = rhsp(ias) + cExchange * 0.50_dp * xmyn(jas) * vecHooXmY(ijs)
            wovp(ias) = wovp(ias) + cExchange * 0.25_dp * xmyn(jas) * vecHooXmY(ijs) / omegaDif
            wovm(ias) = wovm(ias) + cExchange * 0.25_dp * xpyn(jas) * vecHooXmY(ijs) / omegaDif
          else
            rhsp(ias) = rhsp(ias) - cExchange * 0.50_dp * xmyn(jas) * vecHooXmY(ijs)
            wovp(ias) = wovp(ias) - cExchange * 0.25_dp * xmyn(jas) * vecHooXmY(ijs) / omegaDif
            wovm(ias) = wovm(ias) - cExchange * 0.25_dp * xpyn(jas) * vecHooXmY(ijs) / omegaDif
          end if
        end do

      end do

      !!> -RHS^- += H^-_ia[T^-]
      call getHovT(nxoo, nxvv, homo, natom, iatrans, getIA, getIJ, getAB, win,&
        & iAtomStart, ovrXev, grndEigVecs, lrGamma, transChrg, tm, vecHovT)

      rhsm = rhsm - cExchange * vecHovT

      !!> -RHS^+ += - H^+_ia[T^+]
      call getHovT(nxoo, nxvv, homo, natom, iatrans, getIA, getIJ, getAB, win,&
        & iAtomStart, ovrXev, grndEigVecs, lrGamma, transChrg, tp, vecHovT)

      rhsp = rhsp - cExchange * vecHovT

      !!> Woo^- += 0.5 * H^-_ij[T+Z] / Omega_mn, Z part computed later 
      call getHooT(nxov, nxoo, homo, natom, iatrans, getIA, getIJ, win, iAtomStart, &
       & ovrXev, grndEigVecs, lrGamma, transChrg, tm, vecHooT)

      do s = 1, nSpin
        do ij = 1, nxoo(s)
          ijs = ij + soo(s)
          woom(ij,s) = woom(ij,s) - cExchange * 0.5_dp * vecHooT(ijs) / omegaDif
        end do
      end do

      !!> Woo^+ += 0.5 * H^+_ij[T+Z] / Omega_mn, Z part computed later 
      call getHooT(nxov, nxoo, homo, natom, iatrans, getIA, getIJ, win, iAtomStart, &
       & ovrXev, grndEigVecs, lrGamma, transChrg, tp, vecHooT)

      do s = 1, nSpin
        do ij = 1, nxoo(s)
          ijs = ij + soo(s)
          woop(ij,s) = woop(ij,s) + cExchange * 0.5_dp * vecHooT(ijs) / omegaDif
        end do
      end do

      !!> Wvv^- += 0.5 * H^+_ab[T+Z] / Omega_mn, Z part computed later 
      call getHvvT(nxov, nxvv, homo, natom, iatrans, getIA, getAB, win, iAtomStart, &
       & ovrXev, grndEigVecs, lrGamma, transChrg, tm, vecHvvT)

      do s = 1, nSpin
        do ab = 1, nxvv(s)
          abs = ab + svv(s)
          wvvm(ab,s) = wvvm(ab,s) - cExchange * 0.5_dp * vecHvvT(abs) / omegaDif
        end do
      end do

    endif
    print *,'Finally got all W'

  end subroutine getNadiaZvectorEqRHS


  !> Solving the (A+B) Z = -R equation via diagonally preconditioned conjugate gradient
  subroutine solveZVectorPrecond(rhs, tSpin, wij, sym, win, nocc_ud, nvir_ud, nxoo_ud, nxvv_ud, &
      & nxov_ud, nxov_rd, iaTrans, getIA, getIJ, getAB, natom, iAtomStart, ovrXev, grndEigVecs, &
      & occNr, sqrOccIA, gammaMat, species0, spinW, onsMEs, orb, transChrg, tRangeSep, lrGamma)

    !> on entry -R, on exit Z
    real(dp), intent(inout) :: rhs(:)

    !> logical spin polarization
    logical, intent(in) :: tSpin

    !> excitation energies (wij = epsion_j - epsilon_i)
    real(dp), intent(in) :: wij(:)

    !> symmetry flag (singlet or triplet)
    character, intent(in) :: sym

    !> sorting index of the excitation energies.
    integer, intent(in) :: win(:)

    !> occupied orbitals per spin channel
    integer, intent(in) :: nocc_ud(:)

    !> virtual orbitals per spin channel
    integer, intent(in) :: nvir_ud(:)

    !> number of occ-occ transitions per spin channel
    integer, intent(in) :: nxoo_ud(:)

    !> number of vir-vir transitions per spin channel
    integer, intent(in) :: nxvv_ud(:)

    !> number of occ-vir transitions per spin channel
    integer, intent(in) :: nxov_ud(:)

    !> number of occupied-virtual transitions (possibly reduced by windowing)
    integer, intent(in) :: nxov_rd

    !> array from pairs of single particles states to compound index
    integer, intent(in) :: iaTrans(:,:,:)

    !> index array for occ-vir single particle excitations
    integer, intent(in) :: getIA(:,:)

    !> index array for occ-occ single particle excitations
    integer, intent(in) :: getIJ(:,:)

    !> index array for vir-vir single particle excitations
    integer, intent(in) :: getAB(:,:)

    !> number of atoms
    integer, intent(in) :: natom

    !> starting position of each atom in the list of orbitals.
    integer, intent(in) :: iAtomStart(:)

    !> overlap times eigenvector. (nOrb, nOrb)
    real(dp), intent(in) :: ovrXev(:,:,:)

    !> eigenvectors (nOrb, nOrb)
    real(dp), intent(in) :: grndEigVecs(:,:,:)

    !> occupation numbers
    real(dp), intent(in) :: occNr(:,:)

    ! Square root of occupation difference between vir and occ states
    real(dp), intent(in) :: sqrOccIA(:)

    !> DFTB gamma matrix (nAtm, nAtom)
    real(dp), intent(in) :: gammaMat(:,:)

    !> chemical species of the atoms
    integer, intent(in) :: species0(:)

    !> ground state spin constants for each species
    real(dp), intent(in) :: spinW(:)

    !> onsite matrix elements for shells (elements between s orbitals on the same shell are ignored)
    real(dp), intent(in), allocatable :: onsMEs(:,:,:,:)

    !> data type for atomic orbital information
    type(TOrbitals), intent(in) :: orb

    !> machinery for transition charges between single particle levels
    type(TTransCharges), intent(in) :: transChrg

    !> is calculation range-separated?
    logical, intent(in) :: tRangeSep

    !> long-range Gamma
    real(dp), allocatable, intent(in) :: lrGamma(:,:)

    integer :: nxov
    integer :: ia, kk, i, a, s, ias, iis, aas
    real(dp) :: rhs2(size(rhs)), rkm1(size(rhs)), zkm1(size(rhs)), pkm1(size(rhs)), apk(size(rhs))
    real(dp) :: qTmp(nAtom), rs, alphakm1, tmp1, tmp2, bkm1
    real(dp), allocatable :: qTr(:), P(:)

    nxov = nxov_rd
    allocate(qTr(nAtom))

    ! diagonal preconditioner
    ! P^-1 = 1 / (A+B)_ia,ia (diagonal of the supermatrix sum A+B)
    allocate(P(nxov))
    do ia = 1, nxov
      qTr = transChrg%qTransIA(ia, iAtomStart, ovrXev, grndEigVecs, getIA, win)
      call hemv(qTmp, gammaMat, qTr)
      if (.not. tSpin) then
        rs = 4.0_dp * dot_product(qTr, qTmp) + wij(ia)
      else
        rs = 2.0_dp * dot_product(qTr, qTmp) + wij(ia)
        rs = rs + 2.0_dp * sum(qTr * qTr * spinW(species0))
      end if

      !! Possibly reorder spin case
      if (tRangeSep) then
        call hemv(qTmp, lrGamma, qTr)
        rs = rs - cExchange * dot_product(qTr, qTmp)
        call indXov(win, ia, getIA, i, a, s)
        iis = iaTrans(i, i, s)
        qTr = transChrg%qTransIJ(iis, iAtomStart, ovrXev, grndEigVecs, getIJ)
        call hemv(qTmp, lrGamma, qTr)
        aas = iaTrans(a, a, s)
        qTr = transChrg%qTransAB(aas, iAtomStart, ovrXev, grndEigVecs, getAB)
        rs = rs - cExchange * dot_product(qTr, qTmp)
      end if

      P(ia) = 1.0_dp / rs
    end do

    ! Free some space, before entering the actionAplusB routine
    deallocate(qTr)

    ! unit vector as initial guess solution
    rhs2(:) = 1.0_dp / sqrt(real(nxov,dp))

    ! action of matrix on vector
    ! we need the singlet action even for triplet excitations!
    call actionAplusB(tSpin, wij, 'S', win, nocc_ud, nvir_ud, nxoo_ud, nxvv_ud, nxov_ud,&
      & nxov_rd, iaTrans, getIA, getIJ, getAB, iAtomStart, ovrXev, grndEigVecs, occNr, sqrOccIA,&
      & gammaMat, species0, spinW, onsMEs, orb, .true., transChrg, rhs2, rkm1, tRangeSep, lrGamma)

    rkm1(:) = rhs - rkm1
    zkm1(:) = P * rkm1
    pkm1(:) = zkm1

    ! Iteration: should be convergent in at most nxov steps for a quadradic surface, so set higher
    do kk = 1, nxov**2

      ! action of matrix on vector
      call actionAplusB(tSpin, wij, 'S', win, nocc_ud, nvir_ud, nxoo_ud, nxvv_ud, nxov_ud,&
         & nxov_rd, iaTrans, getIA, getIJ, getAB, iAtomStart, ovrXev, grndEigVecs, occNr, sqrOccIA,&
         & gammaMat, species0, spinW, onsMEs, orb, .true., transChrg, pkm1, apk, tRangeSep, lrGamma)

      tmp1 = dot_product(rkm1, zkm1)
      tmp2 = dot_product(pkm1, apk)
      alphakm1 = tmp1 / tmp2

      rhs2 = rhs2 + alphakm1 * pkm1

      rkm1 = rkm1 -alphakm1 * apk

      tmp2 = dot_product(rkm1, rkm1)

      ! residual
      if (tmp2 <= epsilon(1.0_dp)**2) then
        exit
      end if

      if (kk == nxov**2) then
        call error("solveZVectorEq : Z vector not converged!")
      end if

      zkm1(:) = P * rkm1

      tmp2 = dot_product(zkm1, rkm1)

      ! Fletcher-Reeves update
      bkm1 = tmp2 / tmp1

      pkm1 = zkm1 + bkm1 * pkm1

    end do

    rhs(:) = rhs2(:)

  end subroutine solveZVectorPrecond

  !> Solving the (A-B) Z = -R equation via diagonally preconditioned conjugate gradient
  subroutine solveZVectorPrecondMinus(rhs, tSpin, wij, sym, win, nocc_ud, nvir_ud, nxoo_ud, nxvv_ud, &
      & nxov_ud, nxov_rd, iaTrans, getIA, getIJ, getAB, natom, iAtomStart, ovrXev, grndEigVecs, &
      & occNr, sqrOccIA, gammaMat, species0, spinW, onsMEs, orb, transChrg, tRangeSep, lrGamma)

    !> on entry -R, on exit Z
    real(dp), intent(inout) :: rhs(:)

    !> logical spin polarization
    logical, intent(in) :: tSpin

    !> excitation energies (wij = epsion_j - epsilon_i)
    real(dp), intent(in) :: wij(:)

    !> symmetry flag (singlet or triplet)
    character, intent(in) :: sym

    !> sorting index of the excitation energies.
    integer, intent(in) :: win(:)

    !> occupied orbitals per spin channel
    integer, intent(in) :: nocc_ud(:)

    !> virtual orbitals per spin channel
    integer, intent(in) :: nvir_ud(:)

    !> number of occ-occ transitions per spin channel
    integer, intent(in) :: nxoo_ud(:)

    !> number of vir-vir transitions per spin channel
    integer, intent(in) :: nxvv_ud(:)

    !> number of occ-vir transitions per spin channel
    integer, intent(in) :: nxov_ud(:)

    !> number of occupied-virtual transitions (possibly reduced by windowing)
    integer, intent(in) :: nxov_rd

    !> array from pairs of single particles states to compound index
    integer, intent(in) :: iaTrans(:,:,:)

    !> index array for occ-vir single particle excitations
    integer, intent(in) :: getIA(:,:)

    !> index array for occ-occ single particle excitations
    integer, intent(in) :: getIJ(:,:)

    !> index array for vir-vir single particle excitations
    integer, intent(in) :: getAB(:,:)

    !> number of atoms
    integer, intent(in) :: natom

    !> starting position of each atom in the list of orbitals.
    integer, intent(in) :: iAtomStart(:)

    !> overlap times eigenvector. (nOrb, nOrb)
    real(dp), intent(in) :: ovrXev(:,:,:)

    !> eigenvectors (nOrb, nOrb)
    real(dp), intent(in) :: grndEigVecs(:,:,:)

    !> occupation numbers
    real(dp), intent(in) :: occNr(:,:)

    ! Square root of occupation difference between vir and occ states
    real(dp), intent(in) :: sqrOccIA(:)

    !> DFTB gamma matrix (nAtm, nAtom)
    real(dp), intent(in) :: gammaMat(:,:)

    !> chemical species of the atoms
    integer, intent(in) :: species0(:)

    !> ground state spin constants for each species
    real(dp), intent(in) :: spinW(:)

    !> onsite matrix elements for shells (elements between s orbitals on the same shell are ignored)
    real(dp), intent(in), allocatable :: onsMEs(:,:,:,:)

    !> data type for atomic orbital information
    type(TOrbitals), intent(in) :: orb

    !> machinery for transition charges between single particle levels
    type(TTransCharges), intent(in) :: transChrg

    !> is calculation range-separated?
    logical, intent(in) :: tRangeSep

    !> long-range Gamma
    real(dp), allocatable, intent(in) :: lrGamma(:,:)

    integer :: nxov
    integer :: ia, kk, i, a, s, ias, iis, aas
    real(dp) :: rhs2(size(rhs)), rkm1(size(rhs)), zkm1(size(rhs)), pkm1(size(rhs)), apk(size(rhs))
    real(dp) :: qTmp(nAtom), rs, alphakm1, tmp1, tmp2, bkm1
    real(dp), allocatable :: qTr(:), P(:)

    nxov = nxov_rd
    allocate(qTr(nAtom))

    ! diagonal preconditioner
    ! P^-1 = 1 / (A-B)_ia,ia (diagonal of the supermatrix sum A+B)
    allocate(P(nxov))

    do ia = 1, nxov
      rs = wij(ia)
      !! Possibly reorder spin case
      if (tRangeSep) then
        qTr = transChrg%qTransIA(ia, iAtomStart, ovrXev, grndEigVecs, getIA, win)
        call hemv(qTmp, gammaMat, qTr)
        call hemv(qTmp, lrGamma, qTr)
        rs = rs + cExchange * dot_product(qTr, qTmp)
        call indXov(win, ia, getIA, i, a, s)
        iis = iaTrans(i, i, s)
        qTr = transChrg%qTransIJ(iis, iAtomStart, ovrXev, grndEigVecs, getIJ)
        call hemv(qTmp, lrGamma, qTr)
        aas = iaTrans(a, a, s)
        qTr = transChrg%qTransAB(aas, iAtomStart, ovrXev, grndEigVecs, getAB)
        rs = rs - cExchange * dot_product(qTr, qTmp)
      end if
      P(ia) = 1.0_dp / rs
    end do

    ! Free some space, before entering the actionAplusB routine
    deallocate(qTr)

    ! unit vector as initial guess solution
    rhs2(:) = 1.0_dp / sqrt(real(nxov,dp))

    ! action of matrix on vector
    ! we need the singlet action even for triplet excitations!
    call actionAminusB(tSpin, wij, win, nocc_ud, nvir_ud, nxoo_ud, nxvv_ud, nxov_ud, nxov_rd,&
         & iaTrans, getIA, getIJ, getAB, iAtomStart, ovrXev, grndEigVecs, occNr, sqrOccIA,&
         & transChrg, rhs2, rkm1, tRangeSep, lrGamma)

    rkm1(:) = rhs - rkm1
    zkm1(:) = P * rkm1
    pkm1(:) = zkm1

    ! Iteration: should be convergent in at most nxov steps for a quadradic surface, so set higher
    do kk = 1, nxov**2

      ! action of matrix on vector
      call actionAminusB(tSpin, wij, win, nocc_ud, nvir_ud, nxoo_ud, nxvv_ud, nxov_ud, nxov_rd,&
           & iaTrans, getIA, getIJ, getAB, iAtomStart, ovrXev, grndEigVecs, occNr, sqrOccIA,&
           & transChrg, pkm1, apk, tRangeSep, lrGamma)

      tmp1 = dot_product(rkm1, zkm1)
      tmp2 = dot_product(pkm1, apk)
      alphakm1 = tmp1 / tmp2

      rhs2 = rhs2 + alphakm1 * pkm1

      rkm1 = rkm1 -alphakm1 * apk

      tmp2 = dot_product(rkm1, rkm1)

      ! residual
      if (tmp2 <= epsilon(1.0_dp)**2) then
        exit
      end if

      if (kk == nxov**2) then
        call error("solveZVectorEq : Z vector not converged!")
      end if

      zkm1(:) = P * rkm1

      tmp2 = dot_product(zkm1, rkm1)

      ! Fletcher-Reeves update
      bkm1 = tmp2 / tmp1

      pkm1 = zkm1 + bkm1 * pkm1

    end do

    rhs(:) = rhs2(:)

  end subroutine solveZVectorPrecondMinus


  !> Calculate Z-dependent parts of the W-vectors and divide diagonal elements of W_ij and W_ab by
  !> 2.
  subroutine calcWvectorZ(zz, win, homo, nmatup, getIA, getIJ, getAB, iaTrans, iAtomStart,&
      & ovrXev, grndEigVecs, gammaMat, grndEigVal, wov, woo, wvv, transChrg, species0, spinW, &
      & tRangeSep, lrGamma)

    !> Z vector
    real(dp), intent(in) :: zz(:)

    !> index array for single particle transitions
    integer, intent(in) :: win(:)

    !> highest occupied level
    integer, intent(in) :: homo(:)

    !> number of same spin excitations
    integer, intent(in) :: nmatup

    !> index array between occ-vir transitions in square and 1D representations
    integer, intent(in) :: getIA(:,:)

    !> index array between occ-occ transitions
    integer, intent(in) :: getIJ(:,:)

    !> index array between vir-vir transitions
    integer, intent(in) :: getAB(:,:)

    !> array from pairs of single particles states to compound index
    integer, intent(in) :: iaTrans(:,:,:)

    !> index array for S and H0 ground state square matrices
    integer, intent(in) :: iAtomStart(:)

    !> overlap times ground state wavefunctions
    real(dp), intent(in) :: ovrXev(:,:,:)

    !> ground state wavefunctions
    real(dp), intent(in) :: grndEigVecs(:,:,:)

    !> softened coulomb matrix
    real(dp), intent(in) :: gammaMat(:,:)

    !> ground state MO-energies
    real(dp), intent(in) :: grndEigVal(:,:)

    !> W vector occupied-virtual part
    real(dp), intent(inout) :: wov(:)

    !> W vector occupied part
    real(dp), intent(inout) :: woo(:,:)

    !> W vector virtual part
    real(dp), intent(inout) :: wvv(:,:)

    !> machinery for transition charges between single particle levels
    type(TTransCharges), intent(in) :: transChrg

    !> central cell chemical species
    integer, intent(in) :: species0(:)

    !> ground state spin derivatives for each species
    real(dp), intent(in) :: spinW(:)

    !> is calculation range-separated?
    logical, intent(in) :: tRangeSep

    !> long-range Gamma
    real(dp), allocatable, intent(in) :: lrGamma(:,:)

    integer :: nxov, natom, nSpin, soo(2), svv(2)
    integer, allocatable :: nxoo(:), nxvv(:), nvir(:)
    integer :: ij, ias, ijs, ab, i, j, a, b, s, iAt1
    real(dp) :: fact
    real(dp), allocatable :: qTr(:), gamxpyq(:), zq(:), zqds(:), vecHooZ(:)
    logical :: tSpin

    nxov = size(zz)
    natom = size(gammaMat, dim=1)
    nSpin = size(grndEigVal, dim=2)

    ALLOCATE(qTr(natom))
    ALLOCATE(gamxpyq(natom))
    ALLOCATE(zq(natom))
    ALLOCATE(nxoo(nSpin))
    ALLOCATE(nxvv(nSpin))
    ALLOCATE(nvir(nSpin))

    nxoo(:) = (homo(:)*(homo(:)+1))/2
    nvir(:) = size(grndEigVecs, dim=1) - homo(:)
    nxvv(:) = (nvir(:)*(nvir(:)+1))/2

    soo(:) = (/ 0, nxoo(1) /)
    svv(:) = (/ 0, nxvv(1) /)

    if ( nSpin == 2 ) then
      tSpin = .true.
      ALLOCATE(zqds(natom))
    else
      tSpin = .false.
    end if

    ! Adding missing epsilon_i * Z_ia term to W_ia
    do ias = 1, nxov
      call indxov(win, ias, getIA, i, a, s)
      wov(ias) = wov(ias) + zz(ias) * grndEigVal(i, s)
    end do

    ! Missing sum_kb 4 K_ijkb Z_kb term in W_ij: zq(iAt1) = sum_kb q^kb(iAt1) Z_kb
    zq(:) = 0.0_dp
    call transChrg%qMatVec(iAtomStart, ovrXev, grndEigVecs, getIA, win, zz, zq)
    call hemv(gamxpyq, gammaMat, zq)

    if (tSpin) then
      zqds(:) = 0.0_dp
      call transChrg%qMatVecDs(iAtomStart, ovrXev, grndEigVecs, getIA, win, zz, zqds)
    end if

    ! sum_iAt1 qTr(iAt1) gamxpyq(iAt1)
    do s = 1, nSpin
      if (s == 1) then
        fact = 1.0_dp
      else
        fact = -1.0_dp
      end if
      do ij = 1, nxoo(s)
        qTr(:) = transChrg%qTransIJ(ij + soo(s), iAtomStart, ovrXev, grndEigVecs, getIJ)
        ! W contains 1/2 for i == j.
        if (.not. tSpin) then
          woo(ij,s) = woo(ij,s) + 4.0_dp * sum(qTr * gamxpyq)
        else
          woo(ij,s) = woo(ij,s) + 2.0_dp * sum(qTr * gamxpyq)
          woo(ij,s) = woo(ij,s) + 2.0_dp * fact * sum(qTr * zqds * spinW(species0))
        end if
      end do
    end do

    if (tRangeSep) then

      allocate(vecHooZ(sum(nxoo)))
      call getHooXY(1, nxoo, homo, natom, iaTrans, getIA, getIJ, win,&
      & iAtomStart, ovrXev, grndEigVecs, lrGamma, transChrg, zz, vecHooZ)

      !! woo should be made 1D
      do s = 1, nSpin
        do ij = 1, nxoo(s)
          i = getIJ(ij + soo(s), 1)
          j = getIJ(ij + soo(s), 2)
          ijs = iaTrans(i, j, s)
          woo(ij,s) = woo(ij,s) + cExchange * vecHooZ(ijs)
        end do
      end do

    end if

    ! Divide diagonal elements of W_ij by 2.
    do s = 1, nSpin
      do ij = 1, nxoo(s)
        i = getIJ(ij + soo(s), 1)
        j = getIJ(ij + soo(s), 2)
        if (i == j) then
          woo(ij,s) = 0.5_dp * woo(ij,s)
        end if
      end do
    end do

    ! Divide diagonal elements of W_ab by 2.
    do s = 1, nSpin
      do ab = 1, nxvv(s)
        a = getAB(ab + svv(s), 1)
        b = getAB(ab + svv(s), 2)
        if (a == b) then
          wvv(ab,s) = 0.5_dp * wvv(ab,s)
        end if
      end do
    end do

  end subroutine calcWvectorZ

  !> Calculate Z-dependent parts of the W-vectors for the NAC vectors
  subroutine calcNadiaWvectorZ(zzp, zzm, win, homo, nmatup, getIA, getIJ, getAB, iaTrans, &
      & iAtomStart, ovrXev, grndEigVecs, gammaMat, grndEigVal, wovp, wovm, woop, woom,    & 
      & wvvm, transChrg, species0, spinW, tRangeSep, lrGamma, omegaDif)

    !> Z^+ vector
    real(dp), intent(in) :: zzp(:)

    !> Z^- vector
    real(dp), intent(in) :: zzm(:)

    !> index array for single particle transitions
    integer, intent(in) :: win(:)

    !> highest occupied level
    integer, intent(in) :: homo(:)

    !> number of same spin excitations
    integer, intent(in) :: nmatup

    !> index array between occ-vir transitions in square and 1D representations
    integer, intent(in) :: getIA(:,:)

    !> index array between occ-occ transitions
    integer, intent(in) :: getIJ(:,:)

    !> index array between vir-vir transitions
    integer, intent(in) :: getAB(:,:)

    !> array from pairs of single particles states to compound index
    integer, intent(in) :: iaTrans(:,:,:)

    !> index array for S and H0 ground state square matrices
    integer, intent(in) :: iAtomStart(:)

    !> overlap times ground state wavefunctions
    real(dp), intent(in) :: ovrXev(:,:,:)

    !> ground state wavefunctions
    real(dp), intent(in) :: grndEigVecs(:,:,:)

    !> softened coulomb matrix
    real(dp), intent(in) :: gammaMat(:,:)

    !> ground state MO-energies
    real(dp), intent(in) :: grndEigVal(:,:)

    !> W^+ vector occupied-virtual part
    real(dp), intent(inout) :: wovp(:)

    !> W^- vector occupied-virtual part
    real(dp), intent(inout) :: wovm(:)

    !> W^+ vector occupied part
    real(dp), intent(inout) :: woop(:,:)

    !> W^- vector occupied part
    real(dp), intent(inout) :: woom(:,:)

    !> W^- vector virtual part
    real(dp), intent(inout) :: wvvm(:,:)

    !> machinery for transition charges between single particle levels
    type(TTransCharges), intent(in) :: transChrg

    !> central cell chemical species
    integer, intent(in) :: species0(:)

    !> ground state spin derivatives for each species
    real(dp), intent(in) :: spinW(:)

    !> is calculation range-separated?
    logical, intent(in) :: tRangeSep

    !> long-range Gamma
    real(dp), allocatable, intent(in) :: lrGamma(:,:)

    real(dp), intent(in) ::omegaDif

    integer :: nxov, natom, nSpin, soo(2), svv(2)
    integer, allocatable :: nxoo(:), nxvv(:), nvir(:)
    integer :: ij, ias, ijs, ab, i, j, a, b, s, abs, iAt1
    real(dp) :: fact
    real(dp), allocatable :: qTr(:), gamxpyq(:), zq(:), zqds(:), vecHooZ(:), vecHvvZ(:)
    logical :: tSpin

    nxov = size(zzp)
    natom = size(gammaMat, dim=1)
    nSpin = size(grndEigVal, dim=2)

    ALLOCATE(qTr(natom))
    ALLOCATE(gamxpyq(natom))
    ALLOCATE(zq(natom))
    ALLOCATE(nxoo(nSpin))
    ALLOCATE(nxvv(nSpin))
    ALLOCATE(nvir(nSpin))

    nxoo(:) = (homo(:)*(homo(:)+1))/2
    nvir(:) = size(grndEigVecs, dim=1) - homo(:)
    nxvv(:) = (nvir(:)*(nvir(:)+1))/2

    soo(:) = (/ 0, nxoo(1) /)
    svv(:) = (/ 0, nxvv(1) /)

    if ( nSpin == 2 ) then
      tSpin = .true.
      ALLOCATE(zqds(natom))
    else
      tSpin = .false.
    end if

    ! Adding missing epsilon_i * Z_ia term to W_ia
    do ias = 1, nxov
      call indxov(win, ias, getIA, i, a, s)
      wovp(ias) = wovp(ias) + zzp(ias) * grndEigVal(i, s) / omegaDif
      wovm(ias) = wovm(ias) + zzm(ias) * grndEigVal(i, s) / omegaDif
    end do

    ! Missing sum_kb 4 K_ijkb Z_kb term in W_ij: zq(iAt1) = sum_kb q^kb(iAt1) Z_kb
    zq(:) = 0.0_dp
    call transChrg%qMatVec(iAtomStart, ovrXev, grndEigVecs, getIA, win, zzp, zq)
    call hemv(gamxpyq, gammaMat, zq)

    if (tSpin) then
      zqds(:) = 0.0_dp
      call transChrg%qMatVecDs(iAtomStart, ovrXev, grndEigVecs, getIA, win, zzp, zqds)
    end if

    ! sum_iAt1 qTr(iAt1) gamxpyq(iAt1)
    do s = 1, nSpin
      if (s == 1) then
        fact = 1.0_dp
      else
        fact = -1.0_dp
      end if
      do ij = 1, nxoo(s)
        qTr(:) = transChrg%qTransIJ(ij + soo(s), iAtomStart, ovrXev, grndEigVecs, getIJ)
        ! W contains 1/2 for i == j.
        ! Force contains factors of 4/2 here
        if (.not. tSpin) then
          woop(ij,s) = woop(ij,s) + 2.0_dp * sum(qTr * gamxpyq) / omegaDif
        else
          woop(ij,s) = woop(ij,s) + 1.0_dp * sum(qTr * gamxpyq) / omegaDif
          woop(ij,s) = woop(ij,s) + 1.0_dp * fact * sum(qTr * zqds * spinW(species0)) / omegaDif
        end if
      end do
    end do

    if (tRangeSep) then

      allocate(vecHooZ(sum(nxoo)))
      allocate(vecHvvZ(sum(nxvv)))

      call getHooXY(1, nxoo, homo, natom, iaTrans, getIA, getIJ, win,&
      & iAtomStart, ovrXev, grndEigVecs, lrGamma, transChrg, zzp, vecHooZ)

      !! woo should be made 1D
      do s = 1, nSpin
        do ij = 1, nxoo(s)
          i = getIJ(ij + soo(s), 1)
          j = getIJ(ij + soo(s), 2)
          ijs = iaTrans(i, j, s)
          woop(ij,s) = woop(ij,s) + cExchange * 0.5_dp * vecHooZ(ijs) / omegaDif
        end do
      end do

      call getHooXY(-1, nxoo, homo, natom, iaTrans, getIA, getIJ, win,&
      & iAtomStart, ovrXev, grndEigVecs, lrGamma, transChrg, zzm, vecHooZ)

      !! woo should be made 1D
      do s = 1, nSpin
        do ij = 1, nxoo(s)
          i = getIJ(ij + soo(s), 1)
          j = getIJ(ij + soo(s), 2)
          ijs = iaTrans(i, j, s)
          woom(ij,s) = woom(ij,s) + cExchange * 0.5_dp * vecHooZ(ijs) / omegaDif
        end do
      end do

      call getHvvXY(-1, nxvv, homo, natom, iatrans, getIA, getAB, win, iAtomStart,&
          & ovrXev, grndEigVecs, lrGamma, transChrg, zzm, vecHvvZ)

      do s = 1, nSpin
        do ab = 1, nxvv(s)
          a = getAB(ab + svv(s), 1)
          b = getAB(ab + svv(s), 2)
          abs = iaTrans(a, b, s)
          wvvm(ab,s) = wvvm(ab,s) + cExchange * 0.5_dp * vecHvvZ(abs) / omegaDif
        end do
      end do

    end if

    ! Division of diagonal elements done in addNadiaGradients

  end subroutine calcNadiaWvectorZ


  !> Write out density matrix, full if rhoSqr is present
  subroutine writeDM(iLev, pc, rhoSqr)

    !> Lable for excited state level
    integer, intent(in) :: iLev

    !> transition density matrix
    real(dp), intent(in) :: pc(:,:,:)

    !> ground state density matrix
    real(dp), intent(in), optional :: rhoSqr(:,:,:)

    integer :: fdUnit, iErr
    integer :: iSpin, nSpin
    character(lc) :: tmpStr, error_string

    nSpin = size(pc, dim=3)

    write(tmpStr, "(A,I0,A)")"DM", iLev, ".dat"

    open(newunit=fdUnit, file=trim(tmpStr), position="rewind", status="replace",&
        & form='unformatted',iostat=iErr)
    if (iErr /= 0) then
      write(error_string, *) "Failure to open density matrix"
      call error(error_string)
    end if

    ! size and spin channels
    do iSpin = 1, nSpin
      write(fdUnit)size(pc, dim=1), iSpin

      if (present(rhoSqr)) then
        write(fdUnit)cmplx(pc(:,:,iSpin)+rhoSqr(:,:,iSpin), 0.0_dp, dp)
      else
        write(fdUnit)cmplx(pc(:,:,iSpin), 0.0_dp, dp)
      end if
    end do

    close(fdUnit)

  end subroutine writeDM


  !> Mulliken population for a square density matrix and overlap
  !> Note: assumes both triangles of both square matrices are filled
  subroutine getExcMulliken(iAtomStart, pc, s, dqex)

    !> indexing array for atoms
    integer, intent(in) :: iAtomStart(:)

    !> density matrix
    real(dp), intent(in) :: pc(:,:)

    !> overlap matrix
    real(dp), intent(in) :: s(:,:)

    !> output atomic charges
    real(dp), intent(out) :: dqex(:)

    real(dp) :: tmp(size(pc,dim=1))
    integer :: iAt1

    @:ASSERT(all(shape(pc)==shape(s)))

    tmp = sum(pc * s,dim=2)
    dqex(:) = 0.0_dp
    do iAt1 = 1, size(dqex)
      dqex(iAt1) = sum(tmp(iAtomStart(iAt1):iAtomStart(iAt1 + 1) -1))
    end do

  end subroutine getExcMulliken


  !> Calculation of force from derivatives of excitation energy
  !> 1. we need the ground and excited Mulliken charges
  !> 2. we need P,(T,Z),W, X + Y from linear response
  !> 3. calculate dsmndr, dhmndr (dS/dR, dh/dR), dgabda (dGamma_{IAt1,IAt2}/dR_{IAt1}),
  !> dgext (dGamma-EXT_{IAt1,k}/dR_{IAt1})
  subroutine addGradients(sym, nxov, natom, species0, iAtomStart, norb, homo, getIA,&
      & getIJ, getAB, win, grndEigVecs, pc, ovrXev, dq_ud, dqex, gammaMat, lrGamma, HubbardU, &
      & spinW, shift, woo, wov, wvv, transChrg, xpy, xmy, coord0, orb, skHamCont, skOverCont, &
      & derivator, rhoSqr, deltaRho, tRangeSep, rangeSep, excgrad)

    !> symmetry of the transition
    character, intent(in) :: sym

    !> number of single particle transitions to include
    integer, intent(in) :: nxov

    !> number of central cell atoms
    integer, intent(in) :: natom

    !> central cell chemical species
    integer, intent(in) :: species0(:)

    !> index array for S and H0 ground state square matrices
    integer, intent(in) :: iAtomStart(:)

    !> number of orbitals for ground state system
    integer, intent(in) :: norb

    !> number of highest occupied state in ground state
    integer, intent(in) :: homo(:)

    !> index array from composite occ-vir transition index to specific single particle states
    integer, intent(in) :: getIA(:,:)

    !> index array from composite occ-occ transition index to specific single particle states
    integer, intent(in) :: getIJ(:,:)

    !> index array from composite vir-vir transition index to specific single particle states
    integer, intent(in) :: getAB(:,:)

    !> single particle transition index
    integer, intent(in) :: win(:)

    !> ground state eigenvectors
    real(dp), intent(in) :: grndEigVecs(:,:,:)

    !> transition density matrix
    real(dp), intent(in) :: pc(:,:,:)

    !> overlap times ground state eigenvectors
    real(dp), intent(in) :: ovrXev(:,:,:)

    !> ground state gross charges
    real(dp), intent(in) :: dq_ud(:,:)

    !> charge differences from ground to excited state
    real(dp), intent(in) :: dqex(:,:)

    !> softened coulomb matrix
    real(dp), intent(in) :: gammaMat(:,:)

    !> electrostatic matrix, long-range corrected
    real(dp), allocatable, intent(in) :: lrGamma(:,:)

    !> ground state Hubbard U values
    real(dp), intent(in) :: HubbardU(:)

    !> ground state spin derivatives for each species
    real(dp), intent(in) :: spinW(:)

    !> ground state potentials (shift vector)
    real(dp), intent(in) :: shift(:)

    !> W vector occupied part
    real(dp), intent(in) :: woo(:,:)

    !> W vector occupied-virtual part
    real(dp), intent(in) :: wov(:)

    !> W vector virtual part
    real(dp), intent(in) :: wvv(:,:)

    !> machinery for transition charges between single particle levels
    type(TTransCharges), intent(in) :: transChrg

    !> X+Y Furche term
    real(dp), intent(in) :: xpy(:)

    !> X-Y Furche term
    real(dp), intent(in) :: xmy(:)

    !> central cell atomic coordinates
    real(dp), intent(in) :: coord0(:,:)

    !> data type for atomic orbital information
    type(TOrbitals), intent(in) :: orb

    !> H0 data
    type(TSlakoCont), intent(in) :: skHamCont

    !> overlap data
    type(TSlakoCont), intent(in) :: skOverCont

    !> Differentiator for the non-scc matrices
    class(TNonSccDiff), intent(in) :: derivator

    !> ground state density matrix
    real(dp), intent(in) :: rhoSqr(:,:,:)

    !> difference density matrix (vs. uncharged atoms)
    real(dp), intent(inout), pointer :: deltaRho(:,:,:)

    !> is calculation range-separated?
    logical, intent(in) :: tRangeSep

    !> Data for range-separated calculation
    type(TRangeSepFunc), allocatable, intent(inout) :: rangeSep

    !> resulting excited state gradient
    real(dp), intent(out) :: excgrad(:,:)


    real(dp), allocatable :: shift_excited(:,:), xpyq(:), xpyqds(:)
    real(dp), allocatable :: shxpyq(:,:), xpycc(:,:,:), wcc(:,:,:), tmp5(:), tmp7(:), tmp11(:)
    real(dp), allocatable :: qTr(:), temp(:), dq(:), dm(:), dsigma(:)
    real(dp), allocatable :: dH0(:,:,:), dSo(:,:,:)
    real(dp), allocatable :: Dens(:,:), SpinDens(:,:)
    real(dp), allocatable :: xmycc(:,:,:), xpyas(:,:,:), xmyas(:,:,:)
    real(dp), allocatable :: overlap(:,:), lrGammaOrb(:,:), gammaLongRangePrime(:,:,:)
    real(dp), allocatable :: PS(:,:,:), DS(:,:,:), SPS(:,:,:), SDS(:,:,:), SX(:,:,:)
    real(dp), allocatable :: XS(:,:,:), SXS(:,:,:), SY(:,:,:), YS(:,:,:), SYS(:,:,:)
    integer :: ia, i, j, a, b, ab, ij, m, n, mu, nu, xyz, iAt1, iAt2, s, ka
    integer :: indalpha, indalpha1, indbeta, indbeta1, soo(2), svv(2)
    integer :: iSp1, iSp2, iSpin, nSpin
    real(dp) :: tmp1, tmp2, tmp3, tmp4, tmp6, tmp8, tmp9, tmp10, rab
    real(dp) :: diffvec(3), dgab(3), tmpVec(3), tmp3a, tmp3b, tmprs, tmprs2, tmps(2)
    real(dp) :: spinFactor
    integer, allocatable :: nxoo(:), nxvv(:), nvir(:), species(:)
    logical :: tSpin

    nSpin = size(grndEigVecs, dim=3)
    tSpin = (nSpin == 2)

    ALLOCATE(shift_excited(natom, nSpin))
    ALLOCATE(xpyq(natom))
    ALLOCATE(shxpyq(natom, nSpin))
    ALLOCATE(xpycc(norb, norb, nSpin))
    ALLOCATE(wcc(norb, norb, nSpin))
    ALLOCATE(qTr(natom))
    ALLOCATE(temp(norb))
    ALLOCATE(tmp5(nSpin))
    ALLOCATE(tmp7(nSpin))

    ALLOCATE(Dens(norb,norb))
    !! TO CHANGE: For tRangeSep density from call seems to be incorrect, have
    !! to recreate it from eigenvectors.
    Dens = 0._dp
    if (tRangeSep) then
      Dens = 0._dp
      call herk(Dens, grndEigVecs(:,1:homo(1),1), alpha=2.0_dp)
    else
      Dens(:,:) = sum(rhoSqr, dim=3)
    endif

    ALLOCATE(dH0(orb%mOrb, orb%mOrb, 3))
    ALLOCATE(dSo(orb%mOrb, orb%mOrb, 3))

    ALLOCATE(nxoo(nSpin))
    ALLOCATE(nxvv(nSpin))
    ALLOCATE(nvir(nSpin))

    nxoo(:) = (homo(:)*(homo(:)+1))/2
    nvir(:) = norb - homo(:)
    nxvv(:) = (nvir(:)*(nvir(:)+1))/2

    soo(:) = (/ 0, nxoo(1) /)
    svv(:) = (/ 0, nxvv(1) /)

    ALLOCATE(dq(natom))
    dq(:) = dq_ud(:,1)

    if (tSpin) then
      ALLOCATE(dm(natom))
      ALLOCATE(xpyqds(natom))
      ALLOCATE(tmp11(nSpin))

      ALLOCATE(SpinDens(norb,norb))
      SpinDens(:,:) = rhoSqr(:,:,1) - rhoSqr(:,:,2)

      ALLOCATE(dsigma(2))
      dsigma(1) = 1.0_dp
      dsigma(2) = -1.0_dp
      dm(:) = dq_ud(:,2)
    end if

    if (tRangeSep) then
      ALLOCATE(xmycc(norb, norb, nSpin))
      ALLOCATE(xpyas(norb, norb, nSpin))
      ALLOCATE(xmyas(norb, norb, nSpin))
      ALLOCATE(PS(norb, norb, nSpin))
      ALLOCATE(DS(norb, norb, nSpin))
      ALLOCATE(SPS(norb, norb, nSpin))
      ALLOCATE(SDS(norb, norb, nSpin))
      ALLOCATE(SX(norb, norb, nSpin))
      ALLOCATE(XS(norb, norb, nSpin))
      ALLOCATE(SXS(norb, norb, nSpin))
      ALLOCATE(SY(norb, norb, nSpin))
      ALLOCATE(YS(norb, norb, nSpin))
      ALLOCATE(SYS(norb, norb, nSpin))
      ALLOCATE(overlap(norb, norb))
      ALLOCATE(lrGammaOrb(norb, norb))
      ALLOCATE(gammaLongRangePrime(3, nAtom, nAtom))

      ! Symmetrize deltaRho
      do mu = 1, norb
        do nu = mu + 1, norb
          deltaRho(mu,nu,:) = deltaRho(nu,mu,:)
        end do
      end do

      ! Compute long-range gamma derivative
      gammaLongRangePrime(:,:,:) = 0._dp
      call rangeSep%getSpecies(species)
      do iAt1 = 1, nAtom
        do iAt2 = 1, nAtom
          if(iAt1 /= iAt2) then
            call getGammaPrimeValue(rangeSep, tmpVec, iAt1, iAt2, coord0, species)
            gammaLongRangePrime(:, iAt1, iAt2) = tmpVec
          end if
        end do
      end do

      ! Symmetrize S (can't we get S from caller?)
      call getSqrS(coord0, nAtom, skOverCont, orb, iAtomStart, species0, overlap)
      call getSqrGamma(nAtom, lrGamma, iAtomStart, lrGammaOrb)

    end if

    excgrad = 0.0_dp

    ! excited state potentials at atomic sites
    do iSpin = 1, nSpin
      call hemv(shift_excited(:,iSpin), gammaMat, dqex(:,iSpin))
    end do

    ! xypq(alpha) = sum_ia (X+Y)_ia q^ia(alpha)
    ! complexity norb * norb * norb
    xpyq(:) = 0.0_dp
    call transChrg%qMatVec(iAtomStart, ovrXev, grndEigVecs, getIA, win, xpy, xpyq)

    ! complexity norb * norb
    shxpyq(:,:) = 0.0_dp
    if (.not. tSpin) then
      if (sym == "S") then
        call hemv(shxpyq(:,1), gammaMat, xpyq)
      else
        shxpyq(:,1) = xpyq(:) * spinW(species0)
      end if
    else
      xpyqds(:) = 0.0_dp
      call transChrg%qMatVecDs(iAtomStart, ovrXev, grndEigVecs, getIA, win, xpy, xpyqds)
      do iSpin = 1, nSpin
        call hemv(shxpyq(:,iSpin), gammaMat, xpyq)
        shxpyq(:,iSpin) = shxpyq(:,iSpin) + dsigma(iSpin) * spinW(species0) * xpyqds
        shxpyq(:,iSpin) = 0.5_dp * shxpyq(:,iSpin)
      end do
    end if

    ! calculate xpycc
    ! (xpycc)_{mu nu} = sum_{ia} (X + Y)_{ia} (grndEigVecs(mu,i)grndEigVecs(nu,a)
    ! + grndEigVecs(nu,i)grndEigVecs(mu,a))
    ! complexity norb * norb * norb
    !
    ! xpycc(mu,nu) = sum_ia (X+Y)_ia grndEigVecs(mu,i) grndEigVecs(nu,a)
    ! xpycc(mu, nu) += sum_ia (X+Y)_ia grndEigVecs(mu,a) grndEigVecs(nu,i)
    xpycc(:,:,:) = 0.0_dp
    do ia = 1, nxov
      call indxov(win, ia, getIA, i, a, iSpin)
      ! should replace with DSYR2 call :
      do nu = 1, norb
        do mu = 1, norb
          xpycc(mu,nu,iSpin) = xpycc(mu,nu,iSpin) + xpy(ia) *&
              & ( grndEigVecs(mu,i,iSpin)*grndEigVecs(nu,a,iSpin)&
              & + grndEigVecs(mu,a,iSpin)*grndEigVecs(nu,i,iSpin) )
        end do
      end do
    end do

    if (tRangeSep) then

      ! Asymmetric contribution: xmycc_as = sum_ias (X-Y)_ias c_mas c_nis
      xmycc(:,:,:) = 0.0_dp
      xpyas(:,:,:) = 0.0_dp
      xmyas(:,:,:) = 0.0_dp
      do ia = 1, nxov
        call indxov(win, ia, getIA, i, a, iSpin)
        ! should replace with DSYR2 call :
        do nu = 1, norb
          do mu = 1, norb
            xmycc(mu,nu,iSpin) = xmycc(mu,nu,iSpin) + xmy(ia) * &
                & ( grndEigVecs(mu,i,iSpin) * grndEigVecs(nu,a,iSpin) &
                & + grndEigVecs(mu,a,iSpin) * grndEigVecs(nu,i,iSpin) )
            xpyas(mu,nu,iSpin) = xpyas(mu,nu,iSpin) + xpy(ia) * &
                &  grndEigVecs(mu,i,iSpin) * grndEigVecs(nu,a,iSpin)
            xmyas(mu,nu,iSpin) = xmyas(mu,nu,iSpin) + xmy(ia) * &
                &  grndEigVecs(mu,i,iSpin) * grndEigVecs(nu,a,iSpin)
          end do
        end do
      end do

      ! Account for normalization of S/T versus spin-polarized X+/-Y
      ! We have (X+Y)^S = 1/sqrt(2) [(X+Y)_up + (X+Y)_dn]
      if (tSpin) then
        xmycc = xmycc / sqrt(2._dp)
        xpyas = xpyas / sqrt(2._dp)
        xmyas = xmyas / sqrt(2._dp)
      end if

      do iSpin = 1, nSpin
        call symm(PS(:,:,iSpin), 'R', overlap, pc(:,:,iSpin), 'U', 1.0_dp, 0.0_dp, nOrb, nOrb)
        call symm(SPS(:,:,iSpin), 'L', overlap, PS(:,:,iSpin), 'U', 1.0_dp, 0.0_dp, nOrb, nOrb)
        call symm(DS(:,:,iSpin), 'R', overlap, deltaRho(:,:,iSpin), 'U', 1.0_dp, 0.0_dp, nOrb, nOrb)
        call symm(SDS(:,:,iSpin), 'L', overlap, DS(:,:,iSpin), 'U', 1.0_dp, 0.0_dp, nOrb, nOrb)
        call symm(XS(:,:,iSpin), 'R', overlap, xpyas(:,:,iSpin), 'U', 1.0_dp, 0.0_dp, nOrb, nOrb)
        call symm(SX(:,:,iSpin), 'L', overlap, xpyas(:,:,iSpin), 'U', 1.0_dp, 0.0_dp, nOrb, nOrb)
        call symm(SXS(:,:,iSpin), 'L', overlap, XS(:,:,iSpin), 'U', 1.0_dp, 0.0_dp, nOrb, nOrb)
        call symm(YS(:,:,iSpin), 'R', overlap, xmyas(:,:,iSpin), 'U', 1.0_dp, 0.0_dp, nOrb, nOrb)
        call symm(SY(:,:,iSpin), 'L', overlap, xmyas(:,:,iSpin), 'U', 1.0_dp, 0.0_dp, nOrb, nOrb)
        call symm(SYS(:,:,iSpin), 'L', overlap, YS(:,:,iSpin), 'U', 1.0_dp, 0.0_dp, nOrb, nOrb)
      end do
    end if

    ! calculate wcc = c_mu,i * W_ij * c_j,nu. We have only W_ab b > a and W_ij j > i:
    ! wcc(m,n) = sum_{pq, p <= q} w_pq (grndEigVecs(mu,p)grndEigVecs(nu,q)
    ! + grndEigVecs(nu,p)grndEigVecs(mu,q))
    ! complexity norb * norb * norb

    ! calculate the occ-occ part
    wcc(:,:,:) = 0.0_dp

    do iSpin = 1, nSpin
      do ij = 1, nxoo(iSpin)
        i = getIJ(ij + soo(iSpin), 1)
        j = getIJ(ij + soo(iSpin), 2)
        ! replace with DSYR2 call :
        do mu = 1, norb
          do nu = 1, norb
            wcc(mu,nu,iSpin) = wcc(mu,nu,iSpin) + woo(ij,iSpin) *&
                & ( grndEigVecs(mu,i,iSpin)*grndEigVecs(nu,j,iSpin)&
                & + grndEigVecs(mu,j,iSpin)*grndEigVecs(nu,i,iSpin) )
          end do
        end do

      end do
    end do

    ! calculate the occ-virt part : the same way as for xpycc
    do ia = 1, nxov
      call indxov(win, ia, getIA, i, a, iSpin)
      ! again replace with DSYR2 call :
      do nu = 1, norb
        do mu = 1, norb
          wcc(mu,nu,iSpin) = wcc(mu,nu,iSpin) + wov(ia) *&
              & ( grndEigVecs(mu,i,iSpin)*grndEigVecs(nu,a,iSpin)&
              & + grndEigVecs(mu,a,iSpin)*grndEigVecs(nu,i,iSpin) )
        end do
      end do
    end do

    ! calculate the virt - virt part
    do iSpin = 1, nSpin
      do ab = 1, nxvv(iSpin)
        a = getAB(ab + svv(iSpin), 1)
        b = getAB(ab + svv(iSpin), 2)
        ! replace with DSYR2 call :
        do mu = 1, norb
          do nu = 1, norb
            wcc(mu,nu,iSpin) = wcc(mu,nu,iSpin) + wvv(ab,iSpin) *&
                & ( grndEigVecs(mu,a,iSpin)*grndEigVecs(nu,b,iSpin)&
                & + grndEigVecs(mu,b,iSpin)*grndEigVecs(nu,a,iSpin) )
          end do
        end do

      end do
    end do

    ! now calculating the force complexity : norb * norb * 3

    ! as have already performed norb**3 operation to get here,
    ! calculate for all atoms

    ! BA: only for non-periodic systems!
    do iAt1 = 1, nAtom
      indalpha = iAtomStart(iAt1)
      indalpha1 = iAtomStart(iAt1 + 1) -1
      iSp1 = species0(iAt1)

      do iAt2 = 1, iAt1 - 1
        indbeta = iAtomStart(iAt2)
        indbeta1 = iAtomStart(iAt2 + 1) -1
        iSp2 = species0(iAt2)

        diffvec = coord0(:,iAt1) - coord0(:,iAt2)
        rab = sqrt(sum(diffvec**2))

        ! now holds unit vector in direction
        diffvec = diffvec / rab

        ! calculate the derivative of gamma
        dgab(:) = diffvec(:) * (-1.0_dp/rab**2 - expGammaPrime(rab, HubbardU(iSp1), HubbardU(iSp2)))

        tmp3a = 0.0_dp
        do iSpin = 1, nSpin
          tmp3a = tmp3a + dq(iAt1) * dqex(iAt2,iSpin) + dqex(iAt1,iSpin) * dq(iAt2)
        end do

        if (.not. tSpin) then
          if (sym == "S") then
            tmp3b = 4.0_dp * xpyq(iAt1) * xpyq(iAt2)
          else
            tmp3b = 0.0_dp
          end if
        else
          tmp3b = 2.0_dp * xpyq(iAt1) * xpyq(iAt2)
        end if

        excgrad(:,iAt1) = excgrad(:,iAt1) + dgab(:) * ( tmp3a + tmp3b )
        excgrad(:,iAt2) = excgrad(:,iAt2) - dgab(:) * ( tmp3a + tmp3b )

        tmp5(:) = shift_excited(iAt1,:) + shift_excited(iAt2,:)
        tmp7(:) = 2.0_dp * ( shxpyq(iAt1,:) + shxpyq(iAt2,:) )

        if (tSpin) then
          tmp9 = spinW(iSp1) * dm(iAt1) + spinW(iSp2) * dm(iAt2)
          tmp11(:) = spinW(iSp1) * dqex(iAt1,:) + spinW(iSp2) * dqex(iAt2,:)
        end if

        if (tRangeSep) then
          tmprs = 0.0_dp
          tmps(:) = 0.0_dp
          do iSpin = 1, nSpin
            do mu = indAlpha, indAlpha1
              do nu = indBeta, indBeta1
                tmprs = tmprs +&
          & ( 2.0_dp * (PS(mu,nu,iSpin) * DS(nu,mu,iSpin) + PS(nu,mu,iSpin) * DS(mu,nu,iSpin)) +&
          &   SPS(mu,nu,iSpin) * deltaRho(mu,nu,iSpin) + SPS(nu,mu,iSpin) * deltaRho(nu,mu,iSpin) +&
          &   pc(mu,nu,iSpin) * SDS(mu,nu,iSpin) + pc(nu,mu,iSpin) * SDS(nu,mu,iSpin) )
                tmprs = tmprs + 2.0_dp *&
          & ( xpyas(mu,nu,iSpin) * SXS(mu,nu,iSpin) + xpyas(nu,mu,iSpin) * SXS(nu,mu,iSpin) +&
          &   SX(mu,nu,iSpin) * XS(mu,nu,iSpin) + SX(nu,mu,iSpin) * XS(nu,mu,iSpin) )
                tmprs = tmprs +&
          & ( XS(mu,nu,iSpin) * XS(nu,mu,iSpin) + XS(nu,mu,iSpin) * XS(mu,nu,iSpin) +&
          &   SXS(mu,nu,iSpin) * xpyas(nu,mu,iSpin) + SXS(nu,mu,iSpin) * xpyas(mu,nu,iSpin) +&
          &   xpyas(mu,nu,iSpin) * SXS(nu,mu,iSpin) + xpyas(nu,mu,iSpin) * SXS(mu,nu,iSpin) +&
          &   SX(mu,nu,iSpin) * SX(nu,mu,iSpin) + SX(nu,mu,iSpin) * SX(mu,nu,iSpin) )
                tmprs = tmprs + 2.0_dp *&
          & ( xmyas(mu,nu,iSpin) * SYS(mu,nu,iSpin) + xmyas(nu,mu,iSpin) * SYS(nu,mu,iSpin) +&
          &   SY(mu,nu,iSpin) * YS(mu,nu,iSpin) + SY(nu,mu,iSpin) * YS(nu,mu,iSpin) )
                tmprs = tmprs -&
          & ( YS(mu,nu,iSpin) * YS(nu,mu,iSpin) + YS(nu,mu,iSpin) * YS(mu,nu,iSpin) +&
          &   SYS(mu,nu,iSpin) * xmyas(nu,mu,iSpin) + SYS(nu,mu,iSpin) * xmyas(mu,nu,iSpin) +&
          &   xmyas(mu,nu,iSpin) * SYS(nu,mu,iSpin) + xmyas(nu,mu,iSpin) * SYS(mu,nu,iSpin) +&
          &   SY(mu,nu,iSpin) * SY(nu,mu,iSpin) + SY(nu,mu,iSpin) * SY(mu,nu,iSpin) )
              end do
            end do
          end do
          ! Factor of two for spin-polarized calculation
          tmprs = cExchange * nSpin * tmprs

          excGrad(:,iAt1) = excGrad(:,iAt1) - 0.125_dp * tmprs * gammaLongRangePrime(:,iAt1,iAt2)
          excGrad(:,iAt2) = excGrad(:,iAt2) + 0.125_dp * tmprs * gammaLongRangePrime(:,iAt1,iAt2)
        end if

        call derivator%getFirstDeriv(dH0, skHamCont, coord0, species0,&
            & iAt1, iAt2, orb)
        call derivator%getFirstDeriv(dSo, skOverCont, coord0, species0,&
            & iAt1, iAt2, orb)

        do xyz = 1, 3

          tmp1 = 0.0_dp
          tmp2 = 0.0_dp
          tmp3 = 0.0_dp
          tmp4 = 0.0_dp
          tmp6 = 0.0_dp
          tmp8 = 0.0_dp
          tmp10 = 0.0_dp
          tmprs2 = 0.0_dp

          do iSpin = 1, nSpin

            do mu = indalpha, indalpha1
              do nu = indbeta, indbeta1
                m = mu - indalpha + 1
                n = nu - indbeta + 1

                tmp1 = tmp1 + 2.0_dp * dH0(n,m,xyz) * pc(mu,nu,iSpin)
                tmp2 = tmp2 + dSo(n,m,xyz) * pc(mu,nu,iSpin) * (shift(iAt1)+shift(iAt2))
                tmp3 = tmp3 - dSo(n,m,xyz) * wcc(mu,nu,iSpin)
                tmp4 = tmp4 + tmp5(iSpin) * dSo(n,m,xyz) * Dens(mu,nu)
                tmp6 = tmp6 + tmp7(iSpin) * dSo(n,m,xyz) * xpycc(mu,nu,iSpin)

                if (tSpin) then
                  tmp8 = tmp8 + tmp9 * dSo(n,m,xyz) * dsigma(iSpin) * pc(mu,nu,iSpin)
                  tmp10 = tmp10 + tmp11(iSpin) * dSo(n,m,xyz) * dsigma(iSpin) * SpinDens(mu,nu)
                end if

                if (tRangeSep) then
                  tmprs = 0.0_dp
                  do ka = 1, nOrb
                    tmprs = tmprs +&
            & ( PS(mu,ka,iSpin) * deltaRho(nu,ka,iSpin) + PS(nu,ka,iSpin) * deltaRho(mu,ka,iSpin) +&
            &   pc(mu,ka,iSpin) * DS(nu,ka,iSpin) + pc(nu,ka,iSpin) * DS(mu,ka,iSpin) ) *&
            &  (lrGammaOrb(mu,ka) + lrGammaOrb(nu,ka))
                    tmprs = tmprs +&
            & ( xpyas(mu,ka,iSpin) * XS(nu,ka,iSpin) + xpyas(ka,mu,iSpin) * SX(ka,nu,iSpin) +&
            &   xpyas(nu,ka,iSpin) * XS(mu,ka,iSpin) + xpyas(ka,nu,iSpin) * SX(ka,mu,iSpin) )*&
            &  (lrGammaOrb(mu,ka) + lrGammaOrb(nu,ka))
                    tmprs = tmprs +&
            & ( xmyas(mu,ka,iSpin) * YS(nu,ka,iSpin) + xmyas(ka,mu,iSpin) * SY(ka,nu,iSpin) +&
            &   xmyas(nu,ka,iSpin) * YS(mu,ka,iSpin) + xmyas(ka,nu,iSpin) * SY(ka,mu,iSpin) ) *&
            &  (lrGammaOrb(mu,ka) + lrGammaOrb(nu,ka))
                    tmprs = tmprs +&
            & ( XS(mu,ka,iSpin) * xpyas(ka,nu,iSpin) + XS(nu,ka,iSpin) * xpyas(ka,mu,iSpin) +&
            &   xpyas(mu,ka,iSpin) * SX(ka,nu,iSpin) + xpyas(nu,ka,iSpin) * SX(ka,mu,iSpin)) *&
            &  (lrGammaOrb(mu,ka) + lrGammaOrb(nu,ka))
                    tmprs = tmprs -&
            & ( YS(mu,ka,iSpin) * xmyas(ka,nu,iSpin) + YS(nu,ka,iSpin) * xmyas(ka,mu,iSpin) +&
            &   xmyas(mu,ka,iSpin) * SY(ka,nu,iSpin) + xmyas(nu,ka,iSpin) * SY(ka,mu,iSpin)) *&
            &  (lrGammaOrb(mu,ka) + lrGammaOrb(nu,ka))
                  end do
                  ! Factor of 2 for spin-polarized calculations
                  tmprs2 = tmprs2 + cExchange * nSpin * dSo(n,m,xyz) * tmprs
                end if

              end do
            end do

          end do

          excgrad(xyz,iAt1) = excgrad(xyz,iAt1)&
              & + tmp1 + tmp2 + tmp4 + tmp6 + tmp3 + tmp8 + tmp10 - 0.25_dp * tmprs2
          excgrad(xyz,iAt2) = excgrad(xyz,iAt2)&
              & - tmp1 - tmp2 - tmp4 - tmp6 - tmp3 - tmp8 - tmp10 + 0.25_dp * tmprs2
        end do
      end do
    end do

  end subroutine addGradients

  !> Calculation of nacv using gradient routine
  subroutine addNadiaGradients(sym, nxov, natom, species0, iAtomStart, norb, homo, getIA,     &
      & getIJ, getAB, win, grndEigVecs, pc, ovrXev, dq_ud, dqex, gammaMat, lrGamma, HubbardU, &
      & spinW, shift, woop, woom, wovp, wovm, wvvp, wvvm, transChrg, xpyn, xmyn, xpym, xmym,  &
      & coord0, orb, skHamCont, skOverCont, derivator, rhoSqr, deltaRho, tRangeSep, rangeSep, &
      & nacv)

    !> symmetry of the transition
    character, intent(in) :: sym

    !> number of single particle transitions to include
    integer, intent(in) :: nxov

    !> number of central cell atoms
    integer, intent(in) :: natom

    !> central cell chemical species
    integer, intent(in) :: species0(:)

    !> index array for S and H0 ground state square matrices
    integer, intent(in) :: iAtomStart(:)

    !> number of orbitals for ground state system
    integer, intent(in) :: norb

    !> number of highest occupied state in ground state
    integer, intent(in) :: homo(:)

    !> index array from composite occ-vir transition index to specific single particle states
    integer, intent(in) :: getIA(:,:)

    !> index array from composite occ-occ transition index to specific single particle states
    integer, intent(in) :: getIJ(:,:)

    !> index array from composite vir-vir transition index to specific single particle states
    integer, intent(in) :: getAB(:,:)

    !> single particle transition index
    integer, intent(in) :: win(:)

    !> ground state eigenvectors
    real(dp), intent(in) :: grndEigVecs(:,:,:)

    !> transition density matrix
    real(dp), intent(in) :: pc(:,:,:)

    !> overlap times ground state eigenvectors
    real(dp), intent(in) :: ovrXev(:,:,:)

    !> ground state gross charges
    real(dp), intent(in) :: dq_ud(:,:)

    !> charge differences from ground to excited state
    real(dp), intent(in) :: dqex(:,:)

    !> softened coulomb matrix
    real(dp), intent(in) :: gammaMat(:,:)

    !> electrostatic matrix, long-range corrected
    real(dp), allocatable, intent(in) :: lrGamma(:,:)

    !> ground state Hubbard U values
    real(dp), intent(in) :: HubbardU(:)

    !> ground state spin derivatives for each species
    real(dp), intent(in) :: spinW(:)

    !> ground state potentials (shift vector)
    real(dp), intent(in) :: shift(:)

    !> W vector occupied part
    real(dp), intent(inout) :: woop(:,:)

    !> W vector occupied part
    real(dp), intent(in) :: woom(:,:)

    !> W vector occupied-virtual part
    real(dp), intent(inout) :: wovp(:)

    !> W vector occupied-virtual part
    real(dp), intent(in) :: wovm(:)

    !> W vector virtual part
    real(dp), intent(inout) :: wvvp(:,:)

    !> W vector virtual part
    real(dp), intent(in) :: wvvm(:,:)

    !> machinery for transition charges between single particle levels
    type(TTransCharges), intent(in) :: transChrg

    !> X+Y Furche term (state n)
    real(dp), intent(in) :: xpyn(:)

    !> X-Y Furche term (state n)
    real(dp), intent(in) :: xmyn(:)

    !> X+Y Furche term (state m)
    real(dp), intent(in) :: xpym(:)

    !> X-Y Furche term (state m)
    real(dp), intent(in) :: xmym(:)

    !> central cell atomic coordinates
    real(dp), intent(in) :: coord0(:,:)

    !> data type for atomic orbital information
    type(TOrbitals), intent(in) :: orb

    !> H0 data
    type(TSlakoCont), intent(in) :: skHamCont

    !> overlap data
    type(TSlakoCont), intent(in) :: skOverCont

    !> Differentiator for the non-scc matrices
    class(TNonSccDiff), intent(in) :: derivator

    !> ground state density matrix
    real(dp), intent(in) :: rhoSqr(:,:,:)

    !> difference density matrix (vs. uncharged atoms)
    real(dp), intent(inout), pointer :: deltaRho(:,:,:)

    !> is calculation range-separated?
    logical, intent(in) :: tRangeSep

    !> Data for range-separated calculation
    type(TRangeSepFunc), allocatable, intent(inout) :: rangeSep

    !> resulting non-adiabatic coupling
    real(dp), intent(out) :: nacv(:,:)


    real(dp), allocatable :: shift_excited(:,:), xpyq(:,:), xpyqds(:,:)
    real(dp), allocatable :: shxpyq(:,:,:), xpycc(:,:,:,:), wcc(:,:,:), tmp5(:), tmp7(:), tmp11(:)
    real(dp), allocatable :: qTr(:), temp(:), dq(:), dm(:), dsigma(:)
    real(dp), allocatable :: dH0(:,:,:), dSo(:,:,:)
    real(dp), allocatable :: Dens(:,:), SpinDens(:,:)
    real(dp), allocatable :: xmycc(:,:,:,:), xpyas(:,:,:,:), xmyas(:,:,:,:)
    real(dp), allocatable :: overlap(:,:), lrGammaOrb(:,:), gammaLongRangePrime(:,:,:)
    real(dp), allocatable :: PS(:,:,:), DS(:,:,:), SPS(:,:,:), SDS(:,:,:), SX(:,:,:,:)
    real(dp), allocatable :: XS(:,:,:,:), SXS(:,:,:,:), SY(:,:,:,:), YS(:,:,:,:), SYS(:,:,:,:)
    real(dp), allocatable :: xpy(:,:), xmy(:,:)
    integer :: ia, i, j, a, b, ab, ij, m, n, mu, nu, xyz, iAt1, iAt2, s, ka
    integer :: indalpha, indalpha1, indbeta, indbeta1, soo(2), svv(2)
    integer :: iSp1, iSp2, iSpin, nSpin, iState
    real(dp) :: tmp1, tmp2, tmp3, tmp4, tmp6, tmp8, tmp9, tmp10, rab
    real(dp) :: diffvec(3), dgab(3), tmpVec(3), tmp3a, tmp3b, tmprs, tmprs2, tmps(2)
    real(dp) :: spinFactor
    integer, allocatable :: nxoo(:), nxvv(:), nvir(:), species(:)
    logical :: tSpin

    nSpin = size(grndEigVecs, dim=3)
    tSpin = (nSpin == 2)

    ALLOCATE(shift_excited(natom, nSpin))
    ALLOCATE(xpyq(natom,2))
    ALLOCATE(shxpyq(natom, nSpin,2))
    ALLOCATE(xpycc(norb, norb, nSpin,2))
    ALLOCATE(wcc(norb, norb, nSpin))
    ALLOCATE(qTr(natom))
    ALLOCATE(temp(norb))
    ALLOCATE(tmp5(nSpin))
    ALLOCATE(tmp7(nSpin))

    !! This should be changed to save memory
    ALLOCATE(xpy(nxov,2))
    ALLOCATE(xmy(nxov,2))
    xpy(:,1) = xpyn
    xpy(:,2) = xpym
    xmy(:,1) = xmyn
    xmy(:,2) = xmym

    ALLOCATE(Dens(norb,norb))
    !! TO CHANGE: For tRangeSep density from call seems to be incorrect, have
    !! to recreate it from eigenvectors.
    Dens = 0._dp
    if (tRangeSep) then
      Dens = 0._dp
      call herk(Dens, grndEigVecs(:,1:homo(1),1), alpha=2.0_dp)
    else
      Dens(:,:) = sum(rhoSqr, dim=3)
    endif

    ALLOCATE(dH0(orb%mOrb, orb%mOrb, 3))
    ALLOCATE(dSo(orb%mOrb, orb%mOrb, 3))

    ALLOCATE(nxoo(nSpin))
    ALLOCATE(nxvv(nSpin))
    ALLOCATE(nvir(nSpin))

    nxoo(:) = (homo(:)*(homo(:)+1))/2
    nvir(:) = norb - homo(:)
    nxvv(:) = (nvir(:)*(nvir(:)+1))/2

    soo(:) = (/ 0, nxoo(1) /)
    svv(:) = (/ 0, nxvv(1) /)

    ALLOCATE(dq(natom))
    dq(:) = dq_ud(:,1)

    if (tSpin) then
      ALLOCATE(dm(natom))
      ALLOCATE(xpyqds(natom,2))
      ALLOCATE(tmp11(nSpin))

      ALLOCATE(SpinDens(norb,norb))
      SpinDens(:,:) = rhoSqr(:,:,1) - rhoSqr(:,:,2)

      ALLOCATE(dsigma(2))
      dsigma(1) = 1.0_dp
      dsigma(2) = -1.0_dp
      dm(:) = dq_ud(:,2)
    end if

    if (tRangeSep) then
      ALLOCATE(xmycc(norb, norb, nSpin, 2))
      ALLOCATE(xpyas(norb, norb, nSpin, 2))
      ALLOCATE(xmyas(norb, norb, nSpin, 2))
      ALLOCATE(PS(norb, norb, nSpin))
      ALLOCATE(DS(norb, norb, nSpin))
      ALLOCATE(SPS(norb, norb, nSpin))
      ALLOCATE(SDS(norb, norb, nSpin))
      ALLOCATE(SX(norb, norb, nSpin, 2))
      ALLOCATE(XS(norb, norb, nSpin, 2))
      ALLOCATE(SXS(norb, norb, nSpin, 2))
      ALLOCATE(SY(norb, norb, nSpin, 2))
      ALLOCATE(YS(norb, norb, nSpin, 2))
      ALLOCATE(SYS(norb, norb, nSpin, 2))
      ALLOCATE(overlap(norb, norb))
      ALLOCATE(lrGammaOrb(norb, norb))
      ALLOCATE(gammaLongRangePrime(3, nAtom, nAtom))

      ! Symmetrize deltaRho
      do mu = 1, norb
        do nu = mu + 1, norb
          deltaRho(mu,nu,:) = deltaRho(nu,mu,:)
        end do
      end do

      ! Compute long-range gamma derivative
      gammaLongRangePrime(:,:,:) = 0._dp
      call rangeSep%getSpecies(species)
      do iAt1 = 1, nAtom
        do iAt2 = 1, nAtom
          if(iAt1 /= iAt2) then
            call getGammaPrimeValue(rangeSep, tmpVec, iAt1, iAt2, coord0, species)
            gammaLongRangePrime(:, iAt1, iAt2) = tmpVec
          end if
        end do
      end do

      ! Symmetrize S (can't we get S from caller?)
      call getSqrS(coord0, nAtom, skOverCont, orb, iAtomStart, species0, overlap)
      call getSqrGamma(nAtom, lrGamma, iAtomStart, lrGammaOrb)

    end if

    nacv = 0.0_dp

    ! excited state potentials at atomic sites
    do iSpin = 1, nSpin
      call hemv(shift_excited(:,iSpin), gammaMat, dqex(:,iSpin))
    end do

    ! xypq(alpha) = sum_ia (X+Y)_ia q^ia(alpha)
    ! complexity norb * norb * norb
    xpyq = 0.0_dp
    do iState = 1, 2
      call transChrg%qMatVec(iAtomStart, ovrXev, grndEigVecs, getIA, win, &
           & xpy(:,iState), xpyq(:,iState))
      
      ! complexity norb * norb
      shxpyq = 0.0_dp
      if (.not. tSpin) then
        if (sym == "S") then
          call hemv(shxpyq(:,1,iState), gammaMat, xpyq(:,iState))
        else
          shxpyq(:,1,iState) = xpyq(:,iState) * spinW(species0)
        end if
      else
        xpyqds(:,iState) = 0.0_dp
        call transChrg%qMatVecDs(iAtomStart, ovrXev, grndEigVecs, getIA, win, &
             & xpy(:,iState), xpyqds(:,iState))
        do iSpin = 1, nSpin
          call hemv(shxpyq(:,iSpin,iState), gammaMat, xpyq(:,iState))
          shxpyq(:,iSpin,iState) = shxpyq(:,iSpin,iState) + dsigma(iSpin) &
               & * spinW(species0) * xpyqds(:,iState)
          shxpyq(:,iSpin,iState) = 0.5_dp * shxpyq(:,iSpin,iState)
        end do
      end if

      ! calculate xpycc
      ! (xpycc)_{mu nu} = sum_{ia} (X + Y)_{ia} (grndEigVecs(mu,i)grndEigVecs(nu,a)
      ! + grndEigVecs(nu,i)grndEigVecs(mu,a))
      ! complexity norb * norb * norb
      !
      ! xpycc(mu,nu) = sum_ia (X+Y)_ia grndEigVecs(mu,i) grndEigVecs(nu,a)
      ! xpycc(mu, nu) += sum_ia (X+Y)_ia grndEigVecs(mu,a) grndEigVecs(nu,i)
      xpycc = 0.0_dp
      do ia = 1, nxov
        call indxov(win, ia, getIA, i, a, iSpin)
        ! should replace with DSYR2 call :
        do nu = 1, norb
          do mu = 1, norb
            xpycc(mu,nu,iSpin,iState) = xpycc(mu,nu,iSpin,iState) + xpy(ia,iState) *&
              & ( grndEigVecs(mu,i,iSpin)*grndEigVecs(nu,a,iSpin)&
              & + grndEigVecs(mu,a,iSpin)*grndEigVecs(nu,i,iSpin) )
          end do
        end do
      end do
    end do

    if (tRangeSep) then

      do iState = 1,2
        ! Asymmetric contribution: xmycc_as = sum_ias (X-Y)_ias c_mas c_nis
        xmycc = 0.0_dp
        xpyas = 0.0_dp
        xmyas = 0.0_dp
        do ia = 1, nxov
          call indxov(win, ia, getIA, i, a, iSpin)
          ! should replace with DSYR2 call :
          do nu = 1, norb
            do mu = 1, norb
               xmycc(mu,nu,iSpin,iState) = xmycc(mu,nu,iSpin,iState) + xmy(ia,iState) * &
                & ( grndEigVecs(mu,i,iSpin) * grndEigVecs(nu,a,iSpin) &
                & + grndEigVecs(mu,a,iSpin) * grndEigVecs(nu,i,iSpin) )
               xpyas(mu,nu,iSpin,iState) = xpyas(mu,nu,iSpin,iState) + xpy(ia,iState) * &
                &  grndEigVecs(mu,i,iSpin) * grndEigVecs(nu,a,iSpin)
               xmyas(mu,nu,iSpin,iState) = xmyas(mu,nu,iSpin,iState) + xmy(ia,iState) * &
                &  grndEigVecs(mu,i,iSpin) * grndEigVecs(nu,a,iSpin)
            end do
          end do
        end do
   
        ! Account for normalization of S/T versus spin-polarized X+/-Y
        ! We have (X+Y)^S = 1/sqrt(2) [(X+Y)_up + (X+Y)_dn]
        if (tSpin) then
          xmycc = xmycc / sqrt(2._dp)
          xpyas = xpyas / sqrt(2._dp)
          xmyas = xmyas / sqrt(2._dp)
        end if

        do iSpin = 1, nSpin
          call symm(PS(:,:,iSpin), 'R', overlap, pc(:,:,iSpin), 'U', 1.0_dp, 0.0_dp, nOrb, nOrb)
          call symm(SPS(:,:,iSpin), 'L', overlap, PS(:,:,iSpin), 'U', 1.0_dp, 0.0_dp, nOrb, nOrb)
          call symm(DS(:,:,iSpin), 'R', overlap, deltaRho(:,:,iSpin), 'U', 1.0_dp, 0.0_dp, nOrb, nOrb)
          call symm(SDS(:,:,iSpin), 'L', overlap, DS(:,:,iSpin), 'U', 1.0_dp, 0.0_dp, nOrb, nOrb)
          call symm(XS(:,:,iSpin,iState), 'R', overlap, xpyas(:,:,iSpin,iState), 'U', 1.0_dp, 0.0_dp, nOrb, nOrb)
          call symm(SX(:,:,iSpin,iState), 'L', overlap, xpyas(:,:,iSpin,iState), 'U', 1.0_dp, 0.0_dp, nOrb, nOrb)
          call symm(SXS(:,:,iSpin,iState), 'L', overlap, XS(:,:,iSpin,iState), 'U', 1.0_dp, 0.0_dp, nOrb, nOrb)
          call symm(YS(:,:,iSpin,iState), 'R', overlap, xmyas(:,:,iSpin,iState), 'U', 1.0_dp, 0.0_dp, nOrb, nOrb)
          call symm(SY(:,:,iSpin,iState), 'L', overlap, xmyas(:,:,iSpin,iState), 'U', 1.0_dp, 0.0_dp, nOrb, nOrb)
          call symm(SYS(:,:,iSpin,iState), 'L', overlap, YS(:,:,iSpin,iState), 'U', 1.0_dp, 0.0_dp, nOrb, nOrb)
        end do
      end do

    end if

    ! calculate wcc = c_mu,i * W_ij * c_j,nu. We have only W_ab b > a and W_ij j > i:
    ! wcc(m,n) = sum_{pq, p <= q} w_pq (grndEigVecs(mu,p)grndEigVecs(nu,q)
    ! + grndEigVecs(nu,p)grndEigVecs(mu,q))
    ! complexity norb * norb * norb

    ! calculate the occ-occ part
    wcc(:,:,:) = 0.0_dp

    do iSpin = 1, nSpin
      do ij = 1, nxoo(iSpin)
        i = getIJ(ij + soo(iSpin), 1)
        j = getIJ(ij + soo(iSpin), 2)
        ! For NACV diagonal elements have not yet been divided by 2 
        if(i==j) then
          woop(ij,iSpin) = 0.5_dp *  woop(ij,iSpin)
        end if
        ! replace with DSYR2 call :
        do mu = 1, norb
          do nu = 1, norb
            wcc(mu,nu,iSpin) = wcc(mu,nu,iSpin) + woop(ij,iSpin) *                      &
                & ( grndEigVecs(mu,i,iSpin)*grndEigVecs(nu,j,iSpin)                     &
                & + grndEigVecs(mu,j,iSpin)*grndEigVecs(nu,i,iSpin) )                   &
                & + woom(ij,iSpin) * ( grndEigVecs(mu,i,iSpin)*grndEigVecs(nu,j,iSpin)  &
                & - grndEigVecs(mu,j,iSpin)*grndEigVecs(nu,i,iSpin) )                              
          end do
        end do

      end do
    end do

    ! calculate the occ-virt part : the same way as for xpycc
    do ia = 1, nxov
      call indxov(win, ia, getIA, i, a, iSpin)
      ! again replace with DSYR2 call :
      do nu = 1, norb
        do mu = 1, norb
          wcc(mu,nu,iSpin) = wcc(mu,nu,iSpin) + wovp(ia) *                     &
              & ( grndEigVecs(mu,i,iSpin)*grndEigVecs(nu,a,iSpin)              &
              & + grndEigVecs(mu,a,iSpin)*grndEigVecs(nu,i,iSpin) )            &
              & + wovm(ia) * ( grndEigVecs(mu,i,iSpin)*grndEigVecs(nu,a,iSpin) &
              & - grndEigVecs(mu,a,iSpin)*grndEigVecs(nu,i,iSpin) ) 
        end do
      end do
    end do

    ! calculate the virt - virt part
    do iSpin = 1, nSpin
      do ab = 1, nxvv(iSpin)
        a = getAB(ab + svv(iSpin), 1)
        b = getAB(ab + svv(iSpin), 2)
        ! For NACV diagonal elements have not yet been divided by 2 
        if(a==b) then
          wvvp(ab,iSpin) = 0.5_dp *  wvvp(ab,iSpin)
        end if       
        ! replace with DSYR2 call :
        do mu = 1, norb
          do nu = 1, norb
            wcc(mu,nu,iSpin) = wcc(mu,nu,iSpin) + wvvp(ab,iSpin) *                     &
                & ( grndEigVecs(mu,a,iSpin)*grndEigVecs(nu,b,iSpin)                    &
                & + grndEigVecs(mu,b,iSpin)*grndEigVecs(nu,a,iSpin) )                  &
                & + wvvm(ab,iSpin) * ( grndEigVecs(mu,a,iSpin)*grndEigVecs(nu,b,iSpin) &
                & - grndEigVecs(mu,b,iSpin)*grndEigVecs(nu,a,iSpin) )
          end do
        end do

      end do
    end do

    ! now calculating the force complexity : norb * norb * 3

    ! as have already performed norb**3 operation to get here,
    ! calculate for all atoms

    ! BA: only for non-periodic systems!
    do iAt1 = 1, nAtom
      indalpha = iAtomStart(iAt1)
      indalpha1 = iAtomStart(iAt1 + 1) -1
      iSp1 = species0(iAt1)

      do iAt2 = 1, iAt1 - 1
        indbeta = iAtomStart(iAt2)
        indbeta1 = iAtomStart(iAt2 + 1) -1
        iSp2 = species0(iAt2)

        diffvec = coord0(:,iAt1) - coord0(:,iAt2)
        rab = sqrt(sum(diffvec**2))

        ! now holds unit vector in direction
        diffvec = diffvec / rab

        ! calculate the derivative of gamma
        dgab(:) = diffvec(:) * (-1.0_dp/rab**2 - expGammaPrime(rab, HubbardU(iSp1), HubbardU(iSp2)))

        tmp3a = 0.0_dp
        do iSpin = 1, nSpin
          tmp3a = tmp3a + dq(iAt1) * dqex(iAt2,iSpin) + dqex(iAt1,iSpin) * dq(iAt2)
        end do

        if (.not. tSpin) then
          if (sym == "S") then
            tmp3b = 4.0_dp * xpyq(iAt1,1) * xpyq(iAt2,2)
          else
            tmp3b = 0.0_dp
          end if
        else
          tmp3b = 2.0_dp * xpyq(iAt1,1) * xpyq(iAt2,2)
        end if

        nacv(:,iAt1) = nacv(:,iAt1) + dgab(:) * ( tmp3a + tmp3b )
        nacv(:,iAt2) = nacv(:,iAt2) - dgab(:) * ( tmp3a + tmp3b )

        tmp5(:) = shift_excited(iAt1,:) + shift_excited(iAt2,:)
        tmp7(:) = 2.0_dp * ( shxpyq(iAt1,:,2) + shxpyq(iAt2,:,2) )

        if (tSpin) then
          tmp9 = spinW(iSp1) * dm(iAt1) + spinW(iSp2) * dm(iAt2)
          tmp11(:) = spinW(iSp1) * dqex(iAt1,:) + spinW(iSp2) * dqex(iAt2,:)
        end if

        if (tRangeSep) then
          tmprs = 0.0_dp
          tmps(:) = 0.0_dp
          do iSpin = 1, nSpin
            do mu = indAlpha, indAlpha1
              do nu = indBeta, indBeta1
                tmprs = tmprs +&
          & ( 2.0_dp * (PS(mu,nu,iSpin) * DS(nu,mu,iSpin) + PS(nu,mu,iSpin) * DS(mu,nu,iSpin)) +&
          &   SPS(mu,nu,iSpin) * deltaRho(mu,nu,iSpin) + SPS(nu,mu,iSpin) * deltaRho(nu,mu,iSpin) +&
          &   pc(mu,nu,iSpin) * SDS(mu,nu,iSpin) + pc(nu,mu,iSpin) * SDS(nu,mu,iSpin) )
                tmprs = tmprs + 2.0_dp *&
          & ( xpyas(mu,nu,iSpin,1) * SXS(mu,nu,iSpin,2) + xpyas(nu,mu,iSpin,2) * SXS(nu,mu,iSpin,1) +&
          &   SX(mu,nu,iSpin,1) * XS(mu,nu,iSpin,2) + SX(nu,mu,iSpin,2) * XS(nu,mu,iSpin,1) )
                tmprs = tmprs +&
          & ( XS(mu,nu,iSpin,1) * XS(nu,mu,iSpin,2) + XS(nu,mu,iSpin,2) * XS(mu,nu,iSpin,1) +&
          &   SXS(mu,nu,iSpin,1) * xpyas(nu,mu,iSpin,2) + SXS(nu,mu,iSpin,2) * xpyas(mu,nu,iSpin,1) +&
          &   xpyas(mu,nu,iSpin,1) * SXS(nu,mu,iSpin,2) + xpyas(nu,mu,iSpin,2) * SXS(mu,nu,iSpin,1) +&
          &   SX(mu,nu,iSpin,1) * SX(nu,mu,iSpin,2) + SX(nu,mu,iSpin,2) * SX(mu,nu,iSpin,1) )
                tmprs = tmprs + 2.0_dp *&
          & ( xmyas(mu,nu,iSpin,1) * SYS(mu,nu,iSpin,2) + xmyas(nu,mu,iSpin,2) * SYS(nu,mu,iSpin,1) +&
          &   SY(mu,nu,iSpin,1) * YS(mu,nu,iSpin,2) + SY(nu,mu,iSpin,2) * YS(nu,mu,iSpin,1) )
                tmprs = tmprs -&
          & ( YS(mu,nu,iSpin,1) * YS(nu,mu,iSpin,2) + YS(nu,mu,iSpin,2) * YS(mu,nu,iSpin,1) +&
          &   SYS(mu,nu,iSpin,1) * xmyas(nu,mu,iSpin,2) + SYS(nu,mu,iSpin,2) * xmyas(mu,nu,iSpin,1) +&
          &   xmyas(mu,nu,iSpin,1) * SYS(nu,mu,iSpin,2) + xmyas(nu,mu,iSpin,2) * SYS(mu,nu,iSpin,1) +&
          &   SY(mu,nu,iSpin,1) * SY(nu,mu,iSpin,2) + SY(nu,mu,iSpin,2) * SY(mu,nu,iSpin,1) )
              end do
            end do
          end do
          ! Factor of two for spin-polarized calculation
          tmprs = cExchange * nSpin * tmprs

          nacv(:,iAt1) = nacv(:,iAt1) - 0.125_dp * tmprs * gammaLongRangePrime(:,iAt1,iAt2)
          nacv(:,iAt2) = nacv(:,iAt2) + 0.125_dp * tmprs * gammaLongRangePrime(:,iAt1,iAt2)
        end if

        call derivator%getFirstDeriv(dH0, skHamCont, coord0, species0,&
            & iAt1, iAt2, orb)
        call derivator%getFirstDeriv(dSo, skOverCont, coord0, species0,&
            & iAt1, iAt2, orb)

        do xyz = 1, 3

          tmp1 = 0.0_dp
          tmp2 = 0.0_dp
          tmp3 = 0.0_dp
          tmp4 = 0.0_dp
          tmp6 = 0.0_dp
          tmp8 = 0.0_dp
          tmp10 = 0.0_dp
          tmprs2 = 0.0_dp

          do iSpin = 1, nSpin

            do mu = indalpha, indalpha1
              do nu = indbeta, indbeta1
                m = mu - indalpha + 1
                n = nu - indbeta + 1

                tmp1 = tmp1 + 2.0_dp * dH0(n,m,xyz) * pc(mu,nu,iSpin)
                tmp2 = tmp2 + dSo(n,m,xyz) * pc(mu,nu,iSpin) * (shift(iAt1)+shift(iAt2))
                tmp3 = tmp3 - dSo(n,m,xyz) * wcc(mu,nu,iSpin)
                tmp4 = tmp4 + tmp5(iSpin) * dSo(n,m,xyz) * Dens(mu,nu)
                tmp6 = tmp6 + tmp7(iSpin) * dSo(n,m,xyz) * xpycc(mu,nu,iSpin,1)

                if (tSpin) then
                  tmp8 = tmp8 + tmp9 * dSo(n,m,xyz) * dsigma(iSpin) * pc(mu,nu,iSpin)
                  tmp10 = tmp10 + tmp11(iSpin) * dSo(n,m,xyz) * dsigma(iSpin) * SpinDens(mu,nu)
                end if

                if (tRangeSep) then
                  tmprs = 0.0_dp
                  do ka = 1, nOrb
                    tmprs = tmprs +&
            & ( PS(mu,ka,iSpin) * deltaRho(nu,ka,iSpin) + PS(nu,ka,iSpin) * deltaRho(mu,ka,iSpin) +&
            &   pc(mu,ka,iSpin) * DS(nu,ka,iSpin) + pc(nu,ka,iSpin) * DS(mu,ka,iSpin) ) *&
            &  (lrGammaOrb(mu,ka) + lrGammaOrb(nu,ka))
                    tmprs = tmprs +&
            & ( xpyas(mu,ka,iSpin,1) * XS(nu,ka,iSpin,2) + xpyas(ka,mu,iSpin,2) * SX(ka,nu,iSpin,1) +&
            &   xpyas(nu,ka,iSpin,1) * XS(mu,ka,iSpin,2) + xpyas(ka,nu,iSpin,2) * SX(ka,mu,iSpin,1) )*&
            &  (lrGammaOrb(mu,ka) + lrGammaOrb(nu,ka))
                    tmprs = tmprs +&
            & ( xmyas(mu,ka,iSpin,1) * YS(nu,ka,iSpin,2) + xmyas(ka,mu,iSpin,2) * SY(ka,nu,iSpin,1) +&
            &   xmyas(nu,ka,iSpin,1) * YS(mu,ka,iSpin,2) + xmyas(ka,nu,iSpin,2) * SY(ka,mu,iSpin,1) ) *&
            &  (lrGammaOrb(mu,ka) + lrGammaOrb(nu,ka))
                    tmprs = tmprs +&
            & ( XS(mu,ka,iSpin,1) * xpyas(ka,nu,iSpin,2) + XS(nu,ka,iSpin,1) * xpyas(ka,mu,iSpin,2) +&
            &   xpyas(mu,ka,iSpin,1) * SX(ka,nu,iSpin,2) + xpyas(nu,ka,iSpin,1) * SX(ka,mu,iSpin,2)) *&
            &  (lrGammaOrb(mu,ka) + lrGammaOrb(nu,ka))
                    tmprs = tmprs -&
            & ( YS(mu,ka,iSpin,1) * xmyas(ka,nu,iSpin,2) + YS(nu,ka,iSpin,1) * xmyas(ka,mu,iSpin,2) +&
            &   xmyas(mu,ka,iSpin,1) * SY(ka,nu,iSpin,2) + xmyas(nu,ka,iSpin,1) * SY(ka,mu,iSpin,2)) *&
            &  (lrGammaOrb(mu,ka) + lrGammaOrb(nu,ka))
                  end do
                  ! Factor of 2 for spin-polarized calculations
                  tmprs2 = tmprs2 + cExchange * nSpin * dSo(n,m,xyz) * tmprs
                end if

              end do
            end do

          end do
          !!print *,tmp1,tmp2,tmp4,tmp6,tmp3,tmp8,tmp10,tmprs2,'tmp1 + tmp2 + tmp4 + tmp6 + tmp3 + tmp8 + tmp10 - 0.25_dp * tmprs2'
          nacv(xyz,iAt1) = nacv(xyz,iAt1)&
              & + tmp1 + tmp2 + tmp4 + tmp6 + tmp3 + tmp8 + tmp10 - 0.25_dp * tmprs2
          nacv(xyz,iAt2) = nacv(xyz,iAt2)&
              & - tmp1 - tmp2 - tmp4 - tmp6 - tmp3 - tmp8 - tmp10 + 0.25_dp * tmprs2
        end do
      end do
    end do

  end subroutine addNadiaGradients


  !> Write out excitations projected onto ground state
  subroutine writeCoeffs(tt, grndEigVecs, occ, tCoeffs, tIncGroundState,&
      & occNatural, naturalOrbs)

    !> T part of the matrix
    real(dp), intent(in) :: tt(:,:,:)

    !> ground state eigenvectors
    real(dp), intent(in) :: grndEigVecs(:,:,:)

    !> ground state occupations
    real(dp), intent(in) :: occ(:,:)

    !> save the coefficients of the natural orbitals
    logical, intent(in) :: tCoeffs

    !> include the ground state as well as the transition part
    logical, intent(in) :: tIncGroundState

    !> Natural orbital occupation numbers
    real(dp), intent(out), optional :: occNatural(:)

    !> Natural orbitals
    real(dp), intent(out), optional :: naturalOrbs(:,:,:)

    real(dp), allocatable :: t2(:,:,:), occtmp(:,:)
    integer :: norb, nSpin, ii, jj, mm, iSpin
    logical :: tSpin

    integer :: fdCoeffs

    norb = size(tt, dim=1)
    nSpin = size(tt, dim=3)
    tSpin = (nSpin == 2)

    if (present(occNatural).or.tCoeffs) then

      ALLOCATE(t2(norb, norb, nSpin))
      t2 = tt
      if (tIncGroundState) then
        do ii = 1, norb
          t2(ii,ii,:) = t2(ii,ii,:) + occ(ii,:)
        end do
      end if

      if (present(occNatural)) then
        naturalOrbs = t2
        call evalCoeffs(naturalOrbs(:,:,1), occNatural, grndEigVecs(:,:,1))
        if (tCoeffs) then
          ALLOCATE(occtmp(size(occ), nSpin))
          occTmp(:,1) = occNatural
        end if
      else
        ALLOCATE(occtmp(size(occ), nSpin))
        occtmp = 0.0_dp
        do iSpin = 1, nSpin
          call evalCoeffs(t2(:,:,iSpin), occtmp(:,iSpin), grndEigVecs(:,:,iSpin))
        end do
      end if

      ! Better to get this by post-processing DFTB+ output, but here for
      ! compatibility at the moment
      if (tCoeffs) then
        open(newunit=fdCoeffs, file=excitedCoefsOut, position="append")
        write(fdCoeffs,*) 'T F'
        if (.not. tSpin) then
          do ii = 1, norb
            jj = norb - ii + 1
            write(fdCoeffs, '(1x,i3,1x,f13.10,1x,f13.10)') ii, occtmp(jj,1), 2.0_dp
            write(fdCoeffs, '(6(f13.10,1x))') (cmplx(t2(mm,jj,1), kind=dp),&
                & mm = 1, norb)
          end do
        else
          do iSpin = 1, nSpin
            write(fdCoeffs,*)
            write(fdCoeffs, '(1x,a,1x,i1)') 'SPIN', iSpin
            do ii = 1, norb
              jj = norb - ii + 1
              write(fdCoeffs, '(1x,i3,1x,f13.10,1x,f13.10)') ii, occtmp(jj,iSpin), 1.0_dp
              write(fdCoeffs, '(6(f13.10,1x))') (cmplx(t2(mm,jj,iSpin), kind=dp),&
                  & mm = 1, norb)
            end do
          end do
        end if

        close(fdCoeffs)
      end if

    end if

  end subroutine writeCoeffs


  !> Project MO density matrix onto ground state orbitals
  subroutine evalCoeffs(t2, occ, eig)

    !> density matrix
    real(dp), intent(inout) :: t2(:,:)

    !> resulting natural orbital occupations
    real(dp), intent(out) :: occ(:)

    !> 'natural' eigenvectors
    real(dp), intent(in) :: eig(:,:)

    real(dp), allocatable :: coeffs(:,:)

    ALLOCATE(coeffs(size(occ),size(occ)))

    call heev(t2, occ, 'U', 'V')
    call gemm(coeffs, eig, t2)
    t2 = coeffs

  end subroutine evalCoeffs


  !> Write out transitions from ground to excited state along with single particle transitions and
  !> dipole strengths
  subroutine writeExcitations(sym, osz, nexc, nmatup, getIA, win, eval, xpy, wij, &
      & fdXPlusY, fdTrans, fdTransDip, transitionDipoles, tWriteTagged,&
      & fdTagged, taggedWriter, fdExc, Ssq)

    !> Symmetry label for the type of transition
    character, intent(in) :: sym

    !> oscillator strengths for transitions from ground to excited states
    real(dp), intent(in) :: osz(:)

    !> number of excited states to solve for
    integer, intent(in) :: nexc

    !> number of same spin excitations
    integer, intent(in) :: nmatup

    !> index array between transitions in square and 1D representations
    integer, intent(in) :: getIA(:,:)

    !> index array for single particle excitations
    integer, intent(in) :: win(:)

    !> excitation energies
    real(dp), intent(in) :: eval(:)

    !> eigenvectors of excited states (X+Y)
    real(dp), intent(in) :: xpy(:,:)

    !> single particle excitation energies
    real(dp), intent(in) :: wij(:)

    !> single particle transition dipole moments
    real(dp), intent(in) :: transitionDipoles(:,:)

    !> should tagged information be written out
    logical, intent(in) :: tWriteTagged

    !> file unit for transition dipoles
    integer, intent(in) :: fdTransDip

    !> file unit for X+Y data
    integer, intent(in) :: fdXPlusY

    !> file unit for transitions
    integer, intent(in) :: fdTrans

    !> file unit for tagged output (> -1 for write out)
    integer, intent(in) :: fdTagged

    !> tagged writer
    type(TTaggedWriter), intent(inout) :: taggedWriter

    !> file unit for excitation energies
    integer, intent(in) :: fdExc

    !> For spin polarized systems, measure of spin
    real(dp), intent(in), optional :: Ssq(:)

    integer :: nmat
    integer :: ii, jj, iweight, indo, m, n, s
    real(dp), allocatable :: wvec(:)
    integer, allocatable :: wvin(:)
    real(dp) :: weight, wvnorm
    logical :: updwn, tSpin
    character :: sign
    type(TDegeneracyFind) :: DegeneracyFind
    logical :: tDegenerate
    integer, allocatable :: degenerate(:,:)
    real(dp), allocatable :: oDeg(:)

    tSpin = present(Ssq)
    nmat = size(wij)

    allocate(wvec(nmat))
    allocate(wvin(nmat))
    wvec(:) = 0.0_dp
    wvin(:) = 0

    if(fdXplusY /= -1) then
      write(fdXPlusY,*) nmat, nexc
    end if

    do ii = 1, nexc
      if (eval(ii) > 0.0_dp) then

        ! calculate weight of single particle transitions
        wvec(:) = xpy(:,ii)**2
        wvnorm = 1.0_dp / sqrt(sum(wvec**2))
        wvec(:) = wvec * wvnorm

        ! find largest coefficient in CI - should use maxloc
        call index_heap_sort(wvin,wvec)
        wvin = wvin(size(wvin):1:-1)
        wvec = wvec(wvin)

        weight = wvec(1)
        iweight = wvin(1)

        call indxov(win, iweight, getIA, m, n, s)
        sign = sym
        if (tSpin) then
          sign = " "
          write(fdExc,&
              & '(1x,f10.3,4x,f14.8,2x,i5,3x,a,1x,i5,7x,f6.3,2x,f10.3,4x,&
              & f6.3)')&
              & Hartree__eV * sqrt(eval(ii)), osz(ii), m, '->', n, weight,&
              & Hartree__eV * wij(iWeight), Ssq(ii)
        else
          write(fdExc,&
              & '(1x,f10.3,4x,f14.8,5x,i5,3x,a,1x,i5,7x,f6.3,2x,f10.3,6x,a)')&
              & Hartree__eV * sqrt(eval(ii)), osz(ii), m, '->', n, weight,&
              & Hartree__eV * wij(iWeight), sign
        end if

        if(fdXplusY /= -1) then
          if (tSpin) then
            updwn = (win(iweight) <= nmatup)
            sign = "D"
            if (updwn) sign = "U"
          end if
          write(fdXPlusY,'(1x,i5,3x,a,3x,ES17.10)') ii, sign, sqrt(eval(ii))
          write(fdXPlusY,'(6(1x,ES17.10))') xpy(:,ii)
        endif

        if (fdTrans /= -1) then
          write(fdTrans, '(2x,a,T12,i5,T21,ES17.10,1x,a,2x,a)')&
              & 'Energy ', ii,  Hartree__eV * sqrt(eval(ii)), 'eV', sign
          write(fdTrans,*)
          write(fdTrans,'(2x,a,9x,a,8x,a)')'Transition', 'Weight', 'KS [eV]'
          write(fdTrans,'(1x,45("="))')

          sign = " "
          do jj = 1, nmat
            !if (wvec(jj) < 1e-4_dp) exit ! ??????
            indo = wvin(jj)
            call indxov(win, indo, getIA, m, n, s)
            if (tSpin) then
              updwn = (win(indo) <= nmatup)
              sign = "D"
              if (updwn) sign = "U"
            end if
            write(fdTrans, '(i5,3x,a,1x,i5,1x,1a,T22,f10.8,T33,f14.8)')&
                & m, '->', n, sign, wvec(jj), Hartree__eV * wij(wvin(jj))
          end do
          write(fdTrans,*)
        end if

        if (fdTransDip /= -1) then
          write(fdTransDip, '(1x,i5,1x,f10.3,2x,3(ES14.6))')&
              & ii, Hartree__eV * sqrt(eval(ii)), (transitionDipoles(ii,jj)&
              & * au__Debye, jj=1,3)
        end if

      else

        ! find largest coefficient in CI - should use maxloc
        call index_heap_sort(wvin,wvec)
        wvin = wvin(size(wvin):1:-1)
        wvec = wvec(wvin)

        weight = wvec(1)
        iweight = wvin(1)
        call indxov(win, iWeight, getIA, m, n, s)
        sign = sym

        if (tSpin) then
          sign = " "
          write(fdExc,&
              & '(6x,A,T12,4x,f14.8,2x,i5,3x,a,1x,i5,7x,A,2x,f10.3,4x,f6.3)')&
              & '< 0', osz(ii), m, '->', n, '-', Hartree__eV * wij(iWeight),&
              & Ssq(ii)
        else
          write(fdExc,&
              & '(6x,A,T12,4x,f14.8,2x,i5,3x,a,1x,i5,7x,f6.3,2x,f10.3,6x,a)')&
              & '< 0', osz(ii), m, '->', n, weight, Hartree__eV * wij(iWeight), sign
        end if

        if(fdXplusY /= -1) then
          if (tSpin) then
            updwn = (win(iweight) <= nmatup)
            sign = "D"
            if (updwn) sign = "U"
          end if
          write(fdXPlusY,'(1x,i5,3x,a,3x,A)') ii,sign, '-'
        endif

        if (fdTrans /= -1) then
          write(fdTrans, '(2x,a,1x,i5,5x,a,1x,a,3x,a)')&
              & 'Energy ', ii,  '-', 'eV', sign
          write(fdTrans,*)
        end if

        if (fdTransDip /= -1) then
          write(fdTransDip, '(1x,i5,1x,A)') ii, '-'
        endif

      end if

    end do

    deallocate(wvec)
    deallocate(wvin)

    if (tWriteTagged) then

      call degeneracyFind%init(elecTolMax)
      call degeneracyFind%degeneracyTest(eval, tDegenerate)
      if (.not.tDegenerate) then
        call taggedWriter%write(fdTagged, tagLabels%excEgy, eval)
        call taggedWriter%write(fdTagged, tagLabels%excOsc, osz)
        ! Since the transition dipole file exists, transition dipoles had been calculated
        if (fdTransDip /= -1) then
          call taggedWriter%write(fdTagged, tagLabels%excDipole,&
              & sqrt(sum(transitionDipoles**2,dim=2)))
        end if
      else
        degenerate = DegeneracyFind%degenerateRanges()
        call taggedWriter%write(fdTagged, tagLabels%excEgy, eval(degenerate(1,:)))
        ! sum oscillator strength over any degenerate levels
        allocate(oDeg(DegeneracyFind%degenerateGroups()))
        do ii = 1, size(oDeg)
          oDeg(ii) = sum(osz(degenerate(1,ii):degenerate(2,ii)))
        end do
        call taggedWriter%write(fdTagged, tagLabels%excOsc, oDeg)
        ! Since the transition dipole file exists, transition dipoles had been calculated
        if (fdTransDip /= -1) then
          oDeg(:) = 0.0_dp
          do ii = 1, size(oDeg)
            oDeg(ii) = sqrt(sum(transitionDipoles(degenerate(1,ii):degenerate(2,ii),:)**2))
          end do
          call taggedWriter%write(fdTagged, tagLabels%excDipole, oDeg)
        end if
      end if

    end if

  end subroutine writeExcitations


  !> Create transition density matrix in MO basis P = T + 1/2 Z symmetric (paper has T + Z
  !> asymmetric) (Zab = Zij = 0, Tia = 0)
  subroutine calcPMatrix(t, rhs, win, getIA, pc)

    !> T matrix
    real(dp), intent(in) :: t(:,:,:)

    !> Z matrix
    real(dp), intent(in) :: rhs(:)

    !> index array for single particle transitions
    integer, intent(in) :: win(:)

    !> array of the occupied->virtual pairs (nTransitions,occ 1 or virtual 2)
    integer, intent(in) :: getIA(:,:)

    !> resulting excited state density matrix
    real(dp), intent(out) :: pc(:,:,:)

    integer :: ias, i, a, s, nSpin

    nSpin = size(pc, dim=3)

    pc = 0.0_dp
    do ias = 1, size(rhs)
      call indxov(win, ias, getIA, i, a, s)
      pc(i,a,s) = rhs(ias)
    end do

    do s = 1, nSpin
      pc(:,:,s) = 0.5_dp * ( pc(:,:,s) + transpose(pc(:,:,s)) )
    end do

    pc = pc + t

  end subroutine calcPMatrix

  !> Create transition density matrix in MO basis P = T + Z 
  !> non-symmetric for NA
  subroutine calcNadiaPMatrix(t, rhsp, rhsm, win, getIA, pc)

    !> T matrix
    real(dp), intent(in) :: t(:,:,:)

    !> Z^+ matrix
    real(dp), intent(in) :: rhsp(:)

    !> Z^- matrix
    real(dp), intent(in) :: rhsm(:)

    !> index array for single particle transitions
    integer, intent(in) :: win(:)

    !> array of the occupied->virtual pairs (nTransitions,occ 1 or virtual 2)
    integer, intent(in) :: getIA(:,:)

    !> resulting excited state density matrix
    real(dp), intent(out) :: pc(:,:,:)

    integer :: ias, i, a, s, nSpin

    nSpin = size(pc, dim=3)

    pc = 0.0_dp
    do ias = 1, size(rhsp)
      call indxov(win, ias, getIA, i, a, s)
      pc(i,a,s) = t(i,a,s) + rhsp(ias) + rhsm(ias)
      pc(a,i,s) = t(a,i,s) + rhsp(ias) - rhsm(ias)
    end do

    !! pc = pc + t

    do s = 1, nSpin
      pc(:,:,s) = 0.5_dp * ( pc(:,:,s) + transpose(pc(:,:,s)) )
    end do

  end subroutine calcNadiaPMatrix

  !> Computes H^+/-_pq [V] as defined in Furche JCP 117 7433 (2002) eq. 20
  !> Here p/q are virtual orbitals and V is either X+Y or X-Y
  subroutine getHvvXY(ipm, nXvv, homo, nAtom, iaTrans, getIA, getAB, win,&
      & iAtomStart, ovrXev, grndEigVecs, lrGamma, transChrg, XorY, vecHvv)

    !> sign s of H in H^(s)[V]
    integer, intent(in) :: ipm

    !> number of vir-vir transitions per spin channel
    integer, intent(in) :: nXvv(:)

    !> occupied orbitals per spin channel
    integer, intent(in) :: homo(:)

    !> number of atoms
    integer, intent(in) :: nAtom

    !> array from pairs of single particles states to compound index
    integer, intent(in) :: iaTrans(:,:,:)

    !> index array for occ-vir single particle excitations
    integer, intent(in) :: getIA(:,:)

    !> index array for vir-vir single particle excitations
    integer, intent(in) :: getAB(:,:)

    !> index array for single particle excitations
    integer, intent(in) :: win(:)

    !> indexing array for square matrices
    integer, intent(in) :: iAtomStart(:)

    !> overlap times ground state eigenvectors
    real(dp), intent(in) :: ovrXev(:,:,:)

    !> ground state eigenvectors
    real(dp), intent(in) :: grndEigVecs(:,:,:)

    !> electrostatic matrix, long-range corrected
    real(dp), allocatable, intent(in) :: lrGamma(:,:)

    !> machinery for transition charges between single particle levels
    type(TTransCharges), intent(in) :: transChrg

    !> RPA eigenvectors, either (X+Y) or (X-Y)
    real(dp), intent(in) :: XorY(:)

    !> Output vector H[V] virtual-virtual
    real(dp), intent(out) :: vecHvv(:)

    real(dp), allocatable :: qIJ(:), gqIJ(:), qX(:,:), Gq(:,:)
    integer :: i, a, b, s, ias, ibs, abs, nOrb, nXov

    nOrb = size(ovrXev, dim=1)
    nXov = size(XorY)

    allocate(qIJ(nAtom))
    allocate(gqIJ(nAtom))
    allocate(qX(nAtom, nXov))
    allocate(Gq(nAtom, nXov))

    qX(:,:) = 0.0_dp
    do ias = 1, nXov
      call indXov(win, ias, getIA, i, a, s)
      do b = homo(s) + 1, nOrb
        ibs = iaTrans(i, b, s)
        abs = iaTrans(a, b, s)
        qIJ = transChrg%qTransAB(abs, iAtomStart, ovrXev, grndEigVecs, getAB)
        qX(:,ias) = qX(:,ias) + qIJ * XorY(ibs)
      end do
    end do

    Gq(:,:)  = 0.0_dp
    do ias = 1, nXov
      qIJ = transChrg%qTransIA(ias, iAtomStart, ovrXev, grndEigVecs, getIA, win)
      call dsymv('U', nAtom, 1.0_dp, lrGamma, nAtom, qIJ, 1, 0.0_dp, gqIJ, 1)
      Gq(:,ias) = gqIJ(:)
    end do

    vecHvv(:) = 0.0_dp
    do abs = 1, sum(nXvv)
      a = getAB(abs, 1)
      b = getAB(abs, 2)
      s = getAB(abs, 3)
      do i = 1, homo(s)
        ias = iaTrans(i, a, s)
        ibs = iaTrans(i, b, s)
        vecHvv(abs) = vecHvv(abs) - ipm * (dot_product(qX(:,ias), Gq(:,ibs))&
            & + ipm * dot_product(Gq(:,ias), qX(:,ibs)))
      end do
    end do

  end subroutine getHvvXY

  !> Computes full range part of H^+_pq [V] as defined in Furche JCP 117 7433 (2002) eq. 20
  !> Here p/q are both virtual or both occupied orbitals and V is either X+Y or X-Y
  !> Note: The full range part of H^- is zero! 
  subroutine getHplusXYfr(sym, nXoo, nXvv, nAtom, iaTrans, getIA, getIJ, getAB, win,&
      & iAtomStart, species0, ovrXev, grndEigVecs, frGamma, spinW, transChrg, XorY, vecHoo, vecHvv)

    !> symmetry of the transition
    character, intent(in) :: sym

    !> number of occ-occ transitions per spin channel
    integer, intent(in) :: nXoo(:)

    !> number of vir-vir transitions per spin channel
    integer, intent(in) :: nXvv(:)

    !> number of atoms
    integer, intent(in) :: nAtom

    !> array from pairs of single particles states to compound index
    integer, intent(in) :: iaTrans(:,:,:)

    !> index array for occ-vir single particle excitations
    integer, intent(in) :: getIA(:,:)

    !> index array for occ-occ single particle excitations
    integer, intent(in) :: getIJ(:,:)

    !> index array for vir-vir single particle excitations
    integer, intent(in) :: getAB(:,:)

    !> index array for single particle excitations
    integer, intent(in) :: win(:)

    !> indexing array for square matrices
    integer, intent(in) :: iAtomStart(:)

    !> central cell chemical species
    integer, intent(in) :: species0(:)

    !> overlap times ground state eigenvectors
    real(dp), intent(in) :: ovrXev(:,:,:)

    !> ground state eigenvectors
    real(dp), intent(in) :: grndEigVecs(:,:,:)

    !> electrostatic matrix, full-range
    real(dp), intent(in) :: frGamma(:,:)

    !> ground state spin derivatives for each species
    real(dp), intent(in) :: spinW(:)

    !> machinery for transition charges between single particle levels
    type(TTransCharges), intent(in) :: transChrg

    !> RPA eigenvectors, either (X+Y) or (X-Y)
    real(dp), intent(in) :: XorY(:)

    !> Output vector H[V] occupied-occupied
    real(dp), intent(out) :: vecHoo(:)

    !> Output vector H[V] virtual-virtual
    real(dp), intent(out) :: vecHvv(:)

    integer :: nSpin, ab, s, abs, svv(2), ij, ijs, soo(2)
    real(dp) :: fact
    real(dp), allocatable  :: xpyq(:), gamxpyq(:), qTr(:), xpyqds(:), gamxpyqds(:)
    logical :: tSpin

    nSpin = size(grndEigVecs, dim=3)
    tSpin = (nSpin == 2)
    ALLOCATE(xpyq(nAtom))
    ALLOCATE(qTr(nAtom))
    ALLOCATE(gamxpyq(nAtom))
    soo(:) = (/ 0, nXoo(1) /)
    svv(:) = (/ 0, nXvv(1) /)

    vecHoo(:) = 0.0_dp
    vecHvv(:) = 0.0_dp
    xpyq(:) = 0.0_dp
    call transChrg%qMatVec(iAtomStart, ovrXev, grndEigVecs, getIA, win, XorY, xpyq)

    if (.not. tSpin) then  ! ---- spin-unpolarized case ----
      ! vecHvv(ab) = sum_jc K_ab,jc (X+Y)_jc
      if (sym == "S") then
        call hemv(gamxpyq, frGamma,  xpyq)
        do ab = 1, nXvv(1)
          qTr(:) = transChrg%qTransAB(ab, iAtomStart, ovrXev, grndEigVecs, getAB)
          vecHvv(ab) = 2.0_dp * sum(qTr * gamxpyq)
        end do
        do ij = 1, nXoo(1)
          qTr(:) = transChrg%qTransIJ(ij, iAtomStart, ovrXev, grndEigVecs, getIJ)
          ! vecHoo(ij) = sum_kb K_ij,kb (X+Y)_kb
          vecHoo(ij) = 2.0_dp * sum(qTr * gamxpyq)
        end do
      else ! triplet case
        do ab = 1, nXvv(1)
          qTr(:) = transChrg%qTransAB(ab, iAtomStart, ovrXev, grndEigVecs, getAB)
          vecHvv(ab) = 2.0_dp * sum(qTr * xpyq * spinW(species0))
        end do
        do ij = 1, nXoo(1)
          qTr(:) = transChrg%qTransIJ(ij, iAtomStart, ovrXev, grndEigVecs, getIJ)
          vecHoo(ij) = 2.0_dp * sum(qTr * xpyq * spinW(species0))
        end do
      end if

    else  ! ---- spin-polarized case -----

      ALLOCATE(xpyqds(nAtom))
      ALLOCATE(gamxpyqds(nAtom))
      xpyqds(:) = 0.0_dp
      call transChrg%qMatVecDs(iAtomStart, ovrXev, grndEigVecs, getIA, win, XorY, xpyqds)

      call hemv(gamxpyq, frGamma,  xpyq)
      do s = 1, 2
        if (s == 1) then
          fact = 1.0_dp
        else
          fact = -1.0_dp
        end if
        do ab = 1, nXvv(s)
          abs = ab + svv(s)
          qTr(:) = transChrg%qTransAB(abs, iAtomStart, ovrXev, grndEigVecs, getAB)
          vecHvv(abs) = sum(qTr * gamxpyq)
          !magnetization part
          vecHvv(abs) = vecHvv(abs) + fact * sum(qTr * xpyqds * spinW(species0))
        end do
        do ij = 1, nXoo(s)
          ijs = ij + soo(s)
          qTr(:) = transChrg%qTransIJ(ijs, iAtomStart, ovrXev, grndEigVecs, getIJ)
          vecHoo(ijs) = sum(qTr * gamxpyq)
          !magnetization part
          vecHoo(ijs) = vecHoo(ijs) + fact * sum(qTr * xpyqds * spinW(species0))
        end do
      end do

    end if

  end subroutine getHplusXYfr

  !> Computes full range part of H^+_pq [M] as defined in Furche JCP 117 7433 (2002) eq. 20
  !> Here pq are arbitrary orbitals and M is a general matrix with ov,oo,vv components
  !> iMode = 1: returns oo components of H
  !> iMode = 2: returns vv components of H
  !> iMode = 3: returns ov components of H
  !> Note: The full range part of H^- is zero! 
  subroutine getHplusMfr(iMode, sym, nXoo, nXvv, nXov, nAtom, iaTrans, getIA, getIJ, getAB, win,&
      & iAtomStart, species0, ovrXev, grndEigVecs, frGamma, spinW, transChrg, matM, vecH)

    !> type of return vector (oo, vv, ov)
    integer , intent(in)  :: iMode 

    !> symmetry of the transition
    character, intent(in) :: sym

    !> number of occ-occ transitions per spin channel
    integer, intent(in) :: nXoo(:)

    !> number of vir-vir transitions per spin channel
    integer, intent(in) :: nXvv(:)

    !> number of occ-vir transitions (both spins)
    integer, intent(in) :: nXov

    !> number of atoms
    integer, intent(in) :: nAtom

    !> array from pairs of single particles states to compound index
    integer, intent(in) :: iaTrans(:,:,:)

    !> index array for occ-vir single particle excitations
    integer, intent(in) :: getIA(:,:)

    !> index array for occ-occ single particle excitations
    integer, intent(in) :: getIJ(:,:)

    !> index array for vir-vir single particle excitations
    integer, intent(in) :: getAB(:,:)

    !> index array for single particle excitations
    integer, intent(in) :: win(:)

    !> indexing array for square matrices
    integer, intent(in) :: iAtomStart(:)

    !> central cell chemical species
    integer, intent(in) :: species0(:)

    !> overlap times ground state eigenvectors
    real(dp), intent(in) :: ovrXev(:,:,:)

    !> ground state eigenvectors
    real(dp), intent(in) :: grndEigVecs(:,:,:)

    !> electrostatic matrix, full-range
    real(dp), intent(in) :: frGamma(:,:)

    !> ground state spin derivatives for each species
    real(dp), intent(in) :: spinW(:)

    !> machinery for transition charges between single particle levels
    type(TTransCharges), intent(in) :: transChrg

    !> input Matrix spin-resolved
    real(dp), intent(in) :: matM(:,:,:)

    !> output vector H[M] 
    real(dp), intent(out) :: vecH(:)

    integer :: nSpin, ab, i, j, a, b, s, abs, svv(2), ij, ijs, soo(2), ias
    real(dp), dimension(2) :: spinFactor = (/ 1.0_dp, -1.0_dp /)
    real(dp), allocatable  :: xpyq(:), gamxpyq(:), qTr(:), xpyqds(:), gamxpyqds(:)
    real(dp), allocatable  :: gamqt(:)
    logical :: tSpin

    if(iMode == 1) then
      @:ASSERT(size(vecH) == sum(nXoo))
    else if(iMode == 2) then
      @:ASSERT(size(vecH) == sum(nXvv))
    else
      @:ASSERT(size(vecH) == nXov)
    end if

    nSpin = size(grndEigVecs, dim=3)
    tSpin = (nSpin == 2)
    if (tSpin) then
      ALLOCATE(xpyqds(natom))
      ALLOCATE(gamxpyqds(natom))
    endif       
    ALLOCATE(xpyq(nAtom))
    ALLOCATE(qTr(nAtom))
    ALLOCATE(gamxpyq(nAtom))
    ALLOCATE(gamqt(nAtom))
    soo(:) = (/ 0, nXoo(1) /)
    svv(:) = (/ 0, nXvv(1) /)

    vecH(:) = 0.0_dp
   ! gamxpyq(iAt2) = sum_ij q_ij(iAt2) M_ij
    gamxpyq(:) = 0.0_dp
    if (tSpin) then
      gamxpyqds(:) = 0.0_dp
    end if

    do s = 1, nSpin
      do ij = 1, nxoo(s)
        i = getIJ(ij + soo(s), 1)
        j = getIJ(ij + soo(s), 2)
        qTr(:) = transChrg%qTransIJ(ij + soo(s), iAtomStart, ovrXev, grndEigVecs, getIJ)
        if (i == j) then
          gamxpyq(:) = gamxpyq(:) + matM(i,j,s) * qTr(:)
          if (tSpin) then
            gamxpyqds(:) = gamxpyqds(:) + matM(i,j,s) * qTr(:) * spinFactor(s)
          end if
        else
          gamxpyq(:) = gamxpyq(:) + (matM(i,j,s) + matM(j,i,s)) * qTr(:)
          if (tSpin) then
            gamxpyqds(:) = gamxpyqds(:) + (matM(i,j,s) + matM(j,i,s)) * qTr(:) * spinFactor(s)
          end if
        end if
      end do
     

      ! gamxpyq(iAt2) += sum_ab q_ab(iAt2) M_ab
      do ab = 1, nxvv(s)
        a = getAB(ab + svv(s), 1)
        b = getAB(ab + svv(s), 2)
        qTr(:) = transChrg%qTransAB(ab + svv(s), iAtomStart, ovrXev, grndEigVecs, getAB)
        if (a == b) then
          gamxpyq(:) = gamxpyq(:) + matM(a,b,s) * qTr(:)
          if (tSpin) then
            gamxpyqds(:) = gamxpyqds(:) + matM(a,b,s) * qTr(:) * spinFactor(s)
          end if
        else
          ! factor 2 because of symmetry of the matrix
          gamxpyq(:) = gamxpyq(:) + (matM(a,b,s) + matM(b,a,s)) * qTr(:)
          if (tSpin) then
            gamxpyqds(:) = gamxpyqds(:) + (matM(a,b,s) + matM(b,a,s)) * qTr(:) * spinFactor(s)
          end if
        end if
      end do
    end do

    ! gamxpyq(iAt2) += sum_ab q_ab(iAt2) M_ia
    do ias = 1, nxov
      i = getIA(ias, 1)
      a = getIA(ias, 2)
      s = getIA(ias, 3)
      qTr(:) = transChrg%qTransIA(ias, iAtomStart, ovrXev, grndEigVecs, getIA, win)
      gamxpyq(:) = gamxpyq(:) + (matM(i,a,s) + matM(a,i,s)) * qTr(:)
      if (tSpin) then
         gamxpyqds(:) = gamxpyqds(:) + (matM(i,a,s) + matM(a,i,s)) * qTr(:) * spinFactor(s)
      end if
    end do

    ! gamqt(iAt1) = sum_iAt2 gamma_iAt1,iAt2 gamxpyq(iAt2)
    call hemv(gamqt, frGamma, gamxpyq)
    
    if(iMode == 1) then
      do s = 1, nSpin
        do ij = 1, nxoo(s)
          i = getIJ(ij + soo(s), 1)
          j = getIJ(ij + soo(s), 2)
          qTr(:) = transChrg%qTransIJ(ij + soo(s), iAtomStart, ovrXev, grndEigVecs, getIJ)
          if (.not. tSpin) then
            vecH(ij + soo(s)) = 2.0_dp * dot_product(gamqt,qTr)
          else 
            vecH(ij + soo(s)) =  dot_product(gamqt,qTr)
            vecH(ij + soo(s)) = vecH(ij + soo(s)) + &
                      & spinFactor(s) * dot_product(gamxpyqds * spinW(species0), qTr)
          end if
        end do
      end do
    else if (iMode == 2) then
      do s = 1, nSpin
        do ab = 1, nxvv(s)
          a = getAB(ab + svv(s), 1)
          b = getAB(ab + svv(s), 2)
          qTr(:) = transChrg%qTransAB(ab + svv(s), iAtomStart, ovrXev, grndEigVecs, getAB)
          if (.not. tSpin) then
            vecH(ab + svv(s)) = 2.0_dp * dot_product(gamqt,qTr)
          else 
            vecH(ab + svv(s)) = dot_product(gamqt,qTr)
            vecH(ab + svv(s)) = vecH(ab + svv(s))  + &
                      & spinFactor(s) * dot_product(gamxpyqds * spinW(species0), qTr)
          end if
        end do
      end do
    else
      do ias = 1, nxov
        i = getIA(ias, 1)
        a = getIA(ias, 2)
        s = getIA(ias, 3)
        qTr(:) = transChrg%qTransIA(ias, iAtomStart, ovrXev, grndEigVecs, getIA, win)
        if (.not. tSpin) then
          vecH(ias) = 2.0_dp * dot_product(gamqt,qTr)
        else 
          vecH(ias) =  dot_product(gamqt,qTr)
          vecH(ias) = vecH(ias)  + &
                      & spinFactor(s) * dot_product(gamxpyqds * spinW(species0), qTr)
        end if
      end do
    end if
  
  end subroutine getHplusMfr


  !> Computes long-range H^+/-_pq [V] as defined in Furche JCP 117 7433 (2002) eq. 20
  !> or for DFTB in Sokolov J. Chem. Theory Comput. 2021, 17, 2266−2282
  !> Here p/q are occupied orbitals and V is an occ-vir vector
  subroutine getHooXY(ipm, nXoo, homo, nAtom, iaTrans, getIA, getIJ, win,&
      & iAtomStart, ovrXev, grndEigVecs, lrGamma, transChrg, XorY, vecHoo)

    !> sign s of H in H^(s)[V]
    integer, intent(in) :: ipm

    !> number of occ-occ transitions per spin channel
    integer, intent(in) :: nXoo(:)

    !> occupied orbitals per spin channel
    integer, intent(in) :: homo(:)

    !> number of atoms
    integer, intent(in) :: nAtom

    !> array from pairs of single particles states to compound index
    integer, intent(in) :: iaTrans(:,:,:)

    !> index array for occ-vir single particle excitations
    integer, intent(in) :: getIA(:,:)

    !> index array for occ-occ single particle excitations
    integer, intent(in) :: getIJ(:,:)

    !> index array for single particle excitations
    integer, intent(in) :: win(:)

    !> indexing array for square matrices
    integer, intent(in) :: iAtomStart(:)

    !> overlap times ground state eigenvectors
    real(dp), intent(in) :: ovrXev(:,:,:)

    !> ground state eigenvectors
    real(dp), intent(in) :: grndEigVecs(:,:,:)

    !> electrostatic matrix, long-range corrected
    real(dp), allocatable, intent(in) :: lrGamma(:,:)

    !> machinery for transition charges between single particle levels
    type(TTransCharges), intent(in) :: transChrg

    !> RPA eigenvectors, either (X+Y) or (X-Y)
    real(dp), intent(in) :: XorY(:)

    !> Output vector H[V] occ-occ
    real(dp), intent(out) :: vecHoo(:)

    real(dp), allocatable :: qIJ(:), gqIJ(:), qX(:,:), Gq(:,:)
    integer :: i, j, a, s, ias, jas, ijs, nOrb, nXov

    nOrb = size(ovrXev, dim=1)
    nXov = size(XorY)

    allocate(qIJ(nAtom))
    allocate(gqIJ(nAtom))
    allocate(qX(nAtom, nXov))
    allocate(Gq(nAtom, nXov))

    qX(:,:) = 0.0_dp
    do ias = 1, nXov
      call indXov(win, ias, getIA, i, a, s)
      do j = 1, homo(s)
        jas = iaTrans(j, a, s)
        ijs = iaTrans(i, j, s)
        qIJ = transChrg%qTransIJ(ijs, iAtomStart, ovrXev, grndEigVecs, getIJ)
        qX(:,ias) = qX(:,ias) + qIJ * XorY(jas)
      end do
    end do

    Gq(:,:)  = 0.0_dp
    do ias = 1, nXov
      call indXov(win, ias, getIA, i, a, s)
      qIJ = transChrg%qTransIA(ias, iAtomStart, ovrXev, grndEigVecs, getIA, win)
      call dsymv('U', nAtom, 1.0_dp, lrGamma, nAtom, qIJ, 1, 0.0_dp, gqIJ, 1)
      Gq(:,ias) = gqIJ
    end do

    vecHoo(:) = 0.0_dp
    do ijs = 1, sum(nXoo)
      i = getIJ(ijs, 1)
      j = getIJ(ijs, 2)
      s = getIJ(ijs, 3)
      do a = homo(s) + 1, nOrb
        ias = iaTrans(i, a, s)
        jas = iaTrans(j, a, s)
        vecHoo(ijs) = vecHoo(ijs) - ipm * (dot_product(qX(:,ias), Gq(:,jas))&
            & + ipm * dot_product(Gq(:,ias), qX(:,jas)))
      end do
    end do

  end subroutine getHooXY

  !> Computes long-range H^+_pq [T] as defined in Furche JCP 117 7433 (2002) eq. 20
  !> or for DFTB in Sokolov J. Chem. Theory Comput. 2021, 17, 2266−2282
  !> Here p is an occupied MO and q is a virtual one
  !> Note that T must be symmetrical and must have no occ-vir components.
  !> Returns also  - H^-_pq [T] if T is anti-symmetric
  subroutine getHovT(nXoo, nXvv, homo, nAtom, iaTrans, getIA, getIJ, getAB, win,&
      & iAtomStart, ovrXev, grndEigVecs, lrGamma, transChrg, t, vecHovT)

    !> number of occ-occ transitions per spin channel
    integer, intent(in) :: nXoo(:)

    !> number of vir-vir transitions per spin channel
    integer, intent(in) :: nXvv(:)

    !> occupied orbitals per spin channel
    integer, intent(in) :: homo(:)

    !> number of atoms
    integer, intent(in) :: nAtom

    !> array from pairs of single particles states to compound index
    integer, intent(in) :: iaTrans(:,:,:)

    !> index array for occ-vir single particle excitations
    integer, intent(in) :: getIA(:,:)

    !> index array for occ-occ single particle excitations
    integer, intent(in) :: getIJ(:,:)

    !> index array for vir-vir single particle excitations
    integer, intent(in) :: getAB(:,:)

    !> index array for single particle excitations
    integer, intent(in) :: win(:)

    !> indexing array for square matrices
    integer, intent(in) :: iAtomStart(:)

    !> overlap times ground state eigenvectors
    real(dp), intent(in) :: ovrXev(:,:,:)

    !> ground state eigenvectors
    real(dp), intent(in) :: grndEigVecs(:,:,:)

    !> electrostatic matrix, long-range corrected
    real(dp), allocatable, intent(in) :: lrGamma(:,:)

    !> machinery for transition charges between single particle levels
    type(TTransCharges), intent(in) :: transChrg

    !> excited state density matrix
    real(dp), intent(in) :: t(:,:,:)

    !> Output vector H[T] occ-vir
    real(dp), intent(out) :: vecHovT(:)

    real(dp), allocatable :: qIJ(:), gqIJ(:), qX(:,:), Gq(:,:)
    integer :: i, j, a, b, s, ias, ibs, abs, ijs, jas, nOrb, nXov, iMx

    nOrb = size(ovrXev, dim=1)
    nXov = size(vecHovT)

    allocate(qIJ(nAtom))
    allocate(gqIJ(nAtom))
    allocate(qX(nAtom, nXov))
    iMx = max(sum(nXoo), sum(nXvv))
    allocate(Gq(nAtom, iMx))

    qX(:,:) = 0.0_dp
    do ias = 1, nXov
      call indXov(win, ias, getIA, i, a, s)
      do b = homo(s) + 1, nOrb
        ibs = iaTrans(i, b, s)
        qIJ = transChrg%qTransIA(ibs, iAtomStart, ovrXev, grndEigVecs, getIA, win)
        qX(:,ias) = qX(:,ias) + qIJ * t(a,b,s)
      end do
    end do

    Gq(:,:)  = 0.0_dp
    do abs = 1, sum(nXvv)
      qIJ = transChrg%qTransAB(abs, iAtomStart, ovrXev, grndEigVecs, getAB)
      call dsymv('U', nAtom, 1.0_dp, lrGamma, nAtom, qIJ, 1, 0.0_dp, gqIJ, 1)
      Gq(:,abs) = gqIJ
    end do

    vecHovT(:) = 0.0_dp
    do ias = 1, nXov
      call indXov(win, ias, getIA, i, a, s)
      do b = homo(s) + 1, nOrb
        ibs = iaTrans(i, b, s)
        abs = iaTrans(a, b, s)
        vecHovT(ias) = vecHovT(ias) - 2.0_dp * dot_product(qX(:,ibs), Gq(:,abs))
      end do
    end do

    qX(:,:) = 0.0_dp
    do ias = 1, nXov
      call indXov(win, ias, getIA, i, a, s)
      do j = 1, homo(s)
        jas = iaTrans(j, a, s)
        qIJ = transChrg%qTransIA(jas, iAtomStart, ovrXev, grndEigVecs, getIA, win)
        qX(:,ias) = qX(:,ias) + qIJ * t(i,j,s)
      end do
    end do

    Gq(:,:)  = 0.0_dp
    do ijs = 1, sum(nXoo)
      i = getIJ(ijs, 1)
      j = getIJ(ijs, 2)
      s = getIJ(ijs, 3)
      qIJ = transChrg%qTransIJ(ijs, iAtomStart, ovrXev, grndEigVecs, getIJ)
      call dsymv('U', nAtom, 1.0_dp, lrGamma, nAtom, qIJ, 1, 0.0_dp, gqIJ, 1)
      Gq(:,ijs) = gqIJ
    end do

    do ias = 1, nXov
      call indXov(win, ias, getIA, i, a, s)
      do j = 1, homo(s)
        jas = iaTrans(j, a, s)
        ijs = iaTrans(i, j, s)
        vecHovT(ias) = vecHovT(ias) - 2.0_dp * dot_product(qX(:,jas), Gq(:,ijs))
      end do
    end do

  end subroutine getHovT

  !> Computes long-range H^+_pq [T] as defined in Furche JCP 117 7433 (2002) eq. 20
  !> or for DFTB in Sokolov J. Chem. Theory Comput. 2021, 17, 2266−2282
  !> Here p and q are virtual 
  !> Note that T must be symmetrical and must have no occ-vir components.
  !> Returns also  - H^-_pq [T] if T is anti-symmetric
  subroutine getHvvT(nXov, nXvv, homo, nAtom, iaTrans, getIA, getAB, win,&
      & iAtomStart, ovrXev, grndEigVecs, lrGamma, transChrg, t, vecHvvT)

    !> number of occ-vir transitions per spin channel
    integer, intent(in) :: nXov

    !> number of vir-vir transitions per spin channel
    integer, intent(in) :: nXvv(:)

    !> occupied orbitals per spin channel
    integer, intent(in) :: homo(:)

    !> number of atoms
    integer, intent(in) :: nAtom

    !> array from pairs of single particles states to compound index
    integer, intent(in) :: iaTrans(:,:,:)

    !> index array for occ-vir single particle excitations
    integer, intent(in) :: getIA(:,:)

    !> index array for occ-occ single particle excitations
    integer, intent(in) :: getAB(:,:)

    !> index array for single particle excitations
    integer, intent(in) :: win(:)

    !> indexing array for square matrices
    integer, intent(in) :: iAtomStart(:)

    !> overlap times ground state eigenvectors
    real(dp), intent(in) :: ovrXev(:,:,:)

    !> ground state eigenvectors
    real(dp), intent(in) :: grndEigVecs(:,:,:)

    !> electrostatic matrix, long-range corrected
    real(dp), allocatable, intent(in) :: lrGamma(:,:)

    !> machinery for transition charges between single particle levels
    type(TTransCharges), intent(in) :: transChrg

    !> excited state density matrix
    real(dp), intent(in) :: t(:,:,:)

    !> Output vector H[T] vir-vir
    real(dp), intent(out) :: vecHvvT(:)

    real(dp), allocatable :: qIJ(:), gqIJ(:), qX(:,:), Gq(:,:), qXa(:,:,:)
    integer :: nOrb, iSpin, nSpin, iMx, svv(2)
    integer :: j, b, s, jbs, i, ibs, jas, abs, a, c, cbs, d, ads, bcs, ab

    nOrb = size(ovrXev, dim=1)
    nSpin = size(t, dim=3)
    svv(:) = (/ 0, nXvv(1) /)

    allocate(qIJ(nAtom))
    allocate(gqIJ(nAtom))
    iMx = max(sum(nXvv), nXov)
    allocate(qX(nAtom, iMx))
    allocate(Gq(nAtom, iMx))

    qX(:,:) = 0.0_dp
    do jbs = 1, nXov
      call indXov(win, jbs, getIA, j, b, s)
      do i = 1, homo(s) 
        ibs = iaTrans(i, b, s)
        qIJ = transChrg%qTransIA(ibs, iAtomStart, ovrXev, grndEigVecs, getIA, win)
        qX(:,jbs) = qX(:,jbs) + qIJ * t(i,j,s)
      end do
    end do

    Gq(:,:)  = 0.0_dp
    do jas = 1, nXov
      qIJ = transChrg%qTransIA(jas, iAtomStart, ovrXev, grndEigVecs, getIA, win)
      call dsymv('U', nAtom, 1.0_dp, lrGamma, nAtom, qIJ, 1, 0.0_dp, gqIJ, 1)
      Gq(:,jas) = gqIJ
    end do

    vecHvvT(:) = 0.0_dp
    do abs = 1, sum(nXvv)
      a = getAB(abs, 1)
      b = getAB(abs, 2)
      s = getAB(abs, 3)
      do j = 1, homo(s)
        jas = iaTrans(j, a, s)
        jbs = iaTrans(j, b, s)
        vecHvvT(abs) = vecHvvT(abs) - 2.0_dp * dot_product(qX(:,jbs), Gq(:,jas))
      end do
    end do

    deallocate(qX)

    Gq(:,:)  = 0.0_dp
    do ads = 1, sum(nXvv)
      qIJ = transChrg%qTransAB(ads, iAtomStart, ovrXev, grndEigVecs, getAB)
      call dsymv('U', nAtom, 1.0_dp, lrGamma, nAtom, qIJ, 1, 0.0_dp, gqIJ, 1)
      Gq(:,ads) = gqIJ(:)
    end do

    !! For qXa_abs = sum_k q_acs t(b,c,s), we need both qXa_abs and qXa_bas
    !! Need for a spin loop, don't think this can be simplified
    do iSpin = 1, nSpin

      allocate(qXa(nAtom, nOrb-homo(iSpin), nOrb-homo(iSpin)))
      qXa(:,:,:) = 0.0_dp
      do b = homo(iSpin) + 1, nOrb
        do c = homo(iSpin) + 1, nOrb
          bcs = iaTrans(b, c, iSpin)
          qIJ = transChrg%qTransAB(bcs, iAtomStart, ovrXev, grndEigVecs, getAB)
          do d = homo(iSpin) + 1, nOrb
            qXa(:,b-homo(iSpin),d-homo(iSpin)) = qXa(:,b-homo(iSpin),d-homo(iSpin)) &
                 & + qIJ * t(c,d,iSpin)
          end do
        end do
      end do

      do ab = 1, nXvv(iSpin)
        a = getAB(ab + svv(iSpin), 1)
        b = getAB(ab + svv(iSpin), 2)
        abs = iaTrans(a, b, iSpin)
        do d = homo(iSpin) + 1, nOrb
          ads = iaTrans(a, d, iSpin)
          vecHvvT(abs) = vecHvvT(abs) &
               & - 2.0_dp * dot_product(qXa(:,b-homo(iSpin),d-homo(iSpin)), Gq(:,ads))
        end do
      end do
      deallocate(qXa)

    end do

  end subroutine getHvvT

  !> Computes long-range H^+_pq [T] as defined in Furche JCP 117 7433 (2002) eq. 20
  !> or for DFTB in Sokolov J. Chem. Theory Comput. 2021, 17, 2266−2282
  !> Here p and q are occupied 
  !> Note that T must be symmetrical and must have no occ-vir components.
  !> Returns also  - H^-_pq [T] if T is anti-symmetric
  subroutine getHooT(nXov, nXoo, homo, nAtom, iaTrans, getIA, getIJ, win,&
      & iAtomStart, ovrXev, grndEigVecs, lrGamma, transChrg, t, vecHooT)

    !> number of occ-vir transitions per spin channel
    integer, intent(in) :: nXov

    !> number of occ-occ transitions per spin channel
    integer, intent(in) :: nXoo(:)

    !> occupied orbitals per spin channel
    integer, intent(in) :: homo(:)

    !> number of atoms
    integer, intent(in) :: nAtom

    !> array from pairs of single particles states to compound index
    integer, intent(in) :: iaTrans(:,:,:)

    !> index array for occ-vir single particle excitations
    integer, intent(in) :: getIA(:,:)

    !> index array for occ-occ single particle excitations
    integer, intent(in) :: getIJ(:,:)

    !> index array for single particle excitations
    integer, intent(in) :: win(:)

    !> indexing array for square matrices
    integer, intent(in) :: iAtomStart(:)

    !> overlap times ground state eigenvectors
    real(dp), intent(in) :: ovrXev(:,:,:)

    !> ground state eigenvectors
    real(dp), intent(in) :: grndEigVecs(:,:,:)

    !> electrostatic matrix, long-range corrected
    real(dp), allocatable, intent(in) :: lrGamma(:,:)

    !> machinery for transition charges between single particle levels
    type(TTransCharges), intent(in) :: transChrg

    !> excited state density matrix
    real(dp), intent(in) :: t(:,:,:)

    !> Output vector H[T] occ-occ
    real(dp), intent(out) :: vecHooT(:)

    real(dp), allocatable :: qIJ(:), gqIJ(:), qX(:,:), Gq(:,:), qXa(:,:,:)
    integer :: nOrb, iSpin, nSpin, iMx, soo(2)
    integer :: i, j, k, a, b, s, ij, ias, ibs, ijs, jas, iks, jks

    nOrb = size(ovrXev, dim=1)
    nSpin = size(t, dim=3)
    soo(:) = (/ 0, nXoo(1) /)

    allocate(qIJ(nAtom))
    allocate(gqIJ(nAtom))
    iMx = max(sum(nXoo), nXov)
    allocate(qX(nAtom, iMx))
    allocate(Gq(nAtom, iMx))

    qX(:,:) = 0.0_dp
    do ias = 1, nXov
      call indXov(win, ias, getIA, i, a, s)
      do b = homo(s) + 1, nOrb
        ibs = iaTrans(i, b, s)
        qIJ = transChrg%qTransIA(ibs, iAtomStart, ovrXev, grndEigVecs, getIA, win)
        qX(:,ias) = qX(:,ias) + qIJ * t(a,b,s)
      end do
    end do

    Gq(:,:)  = 0.0_dp
    do ias = 1, nXov
      qIJ = transChrg%qTransIA(ias, iAtomStart, ovrXev, grndEigVecs, getIA, win)
      call dsymv('U', nAtom, 1.0_dp, lrGamma, nAtom, qIJ, 1, 0.0_dp, gqIJ, 1)
      Gq(:,ias) = gqIJ
    end do

    vecHooT(:) = 0.0_dp
    do ijs = 1, sum(nXoo)
      i = getIJ(ijs, 1)
      j = getIJ(ijs, 2)
      s = getIJ(ijs, 3)
      do a = homo(s) + 1, nOrb
        ias = iaTrans(i, a, s)
        jas = iaTrans(j, a, s)
        vecHooT(ijs) = vecHooT(ijs) - 2.0_dp * dot_product(qX(:,ias), Gq(:,jas))
      end do
    end do

    deallocate(qX)

    Gq(:,:)  = 0.0_dp
    do ijs = 1, sum(nXoo)
      qIJ = transChrg%qTransIJ(ijs, iAtomStart, ovrXev, grndEigVecs, getIJ)
      call dsymv('U', nAtom, 1.0_dp, lrGamma, nAtom, qIJ, 1, 0.0_dp, gqIJ, 1)
      Gq(:,ijs) = gqIJ(:)
    end do

    !! For qXa_ijs = sum_k q_iks t(j,k,s), we need both qXa_ijs and qXa_jis
    !! Need for a spin loop, don't think this can be simplified
    do iSpin = 1, nSpin

      allocate(qXa(nAtom, homo(iSpin), homo(iSpin)))
      qXa(:,:,:) = 0.0_dp
      do i = 1, homo(iSpin)
        do k = 1, homo(iSpin)
          iks = iaTrans(i, k, iSpin)
          qIJ = transChrg%qTransIJ(iks, iAtomStart, ovrXev, grndEigVecs, getIJ)
          do j = 1, homo(iSpin)
            qXa(:,i,j) = qXa(:,i,j) + qIJ * t(j,k,iSpin)
          end do
        end do
      end do

      do ij = 1, nXoo(iSpin)
        i = getIJ(ij + soo(iSpin), 1)
        j = getIJ(ij + soo(iSpin), 2)
        ijs = iaTrans(i, j, iSpin)
        do k = 1, homo(iSpin)
          jks = iaTrans(j, k, iSpin)
          vecHooT(ijs) = vecHooT(ijs) - 2.0_dp * dot_product(qXa(:,i,k), Gq(:,jks))
        end do
      end do
      deallocate(qXa)

    end do

  end subroutine getHooT

  !> Constructs the full overlap matrix S
  subroutine getSqrS(coord, nAtom, skOverCont, orb, iAtomStart, species0, S)
    real(dp), intent(in) :: coord(:,:)
    integer,intent(in) :: nAtom, iAtomStart(:), species0(:)
    type(TSlakoCont), intent(in) :: skOverCont
    type(TOrbitals), intent(in) :: orb
    real(dp), intent(out) :: S(:,:)

    real(dp) :: SBlock(9,9)
    integer :: iAt1, iAt2, mu, nu, m, n

    S(:,:) = 0.0_dp

    do iAt1 = 1, nAtom
      do iAt2 = 1, iAt1-1

        call getSOffsite(coord(:,iAt1), coord(:,iAt2), species0(iAt1), species0(iAt2), orb,&
            & skOverCont, SBlock)

        do mu = iAtomStart(iAt1), iAtomStart(iAt1+1) - 1
          m = mu - iAtomStart(iAt1) + 1
          do nu = iAtomStart(iAt2), iAtomStart(iAt2+1) - 1
            n = nu - iAtomStart(iAt2) + 1
            S(mu,nu) = SBlock(n,m)
            S(nu,mu) = S(mu,nu)
          end do
        end do

      end do
    end do

    do mu = 1, size(S, dim=1)
      S(mu,mu) = 1.0_dp !Diagonal entries
    end do

  end subroutine getSqrS

  !> Constructs a Gamma-Matrix of dimension nOrb instead of nAtoms
  subroutine getSqrGamma(nAtom, lrGamma, iAtomStart, lrGammaOrb)
    real(dp), intent(in) :: lrGamma(:,:)
    integer,intent(in) :: nAtom, iAtomStart(:)
    real(dp), intent(out) :: lrGammaOrb(:,:)
    integer :: at1, at2, mu, nu, indAt1, indAt1p1, indAt2, indAt2p1

    lrGammaOrb(:,:) = 0.0_dp

    do at1 = 1, nAtom
      indAt1 = iAtomStart(at1)
      indAt1p1 = iAtomStart(at1+1) - 1
      do at2 = 1, at1
        indAt2 = iAtomStart(at2)
        indAt2p1 = iAtomStart(at2+1) - 1
        do mu = indAt1, indAt1p1
          do nu = indAt2, indAt2p1
            lrGammaOrb(mu, nu) = lrGamma(at1, at2)
            lrGammaOrb(nu, mu) = lrGammaOrb(mu, nu)
          end do
        end do
      end do
    end do

  end subroutine getSqrGamma

  !> Helper routine to construct overlap
  subroutine getSOffsite(coords1, coords2, iSp1, iSp2, orb, skOverCont, Sblock)
    real(dp), intent(in) :: coords1(:), coords2(:)
    integer, intent(in) :: iSp1, iSp2
    type(TOrbitals), intent(in) :: orb
    type(TSlakoCont), intent(in) :: skOverCont
    real(dp), intent(out) :: Sblock(:,:)

    real(dp) :: interSKOver(getMIntegrals(skOverCont))
    real(dp) :: vect(3), dist

    @:ASSERT(size(coords1) == 3)
    @:ASSERT(size(coords2) == 3)
    @:ASSERT(all(shape(Sblock) >= [orb%mOrb, orb%mOrb]))

    vect(:) = coords2 - coords1
    dist = sqrt(sum(vect**2))
    vect(:) = vect(:) / dist
    call getSKIntegrals(skOverCont, interSKOver, dist, iSp1, iSp2)
    call rotateH0(Sblock, interSKOver, vect(1), vect(2), vect(3), iSp1, iSp2, orb)

  end subroutine getSOffsite

end module dftbp_timedep_linrespgrad
