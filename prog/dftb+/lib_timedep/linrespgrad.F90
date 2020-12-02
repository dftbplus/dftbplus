!--------------------------------------------------------------------------------------------------!
!  DFTB+: general package for performing fast atomistic simulations                                !
!  Copyright (C) 2006 - 2020  DFTB+ developers group                                               !
!                                                                                                  !
!  See the LICENSE file for terms of usage and distribution.                                       !
!--------------------------------------------------------------------------------------------------!

#:include 'common.fypp'

!> Linear response excitations and gradients with respect to atomic coordinates
module dftbp_linrespgrad
  use dftbp_assert
  use dftbp_arpack
  use dftbp_linrespcommon
  use dftbp_commontypes
  use dftbp_slakocont
  use dftbp_shortgamma
  use dftbp_accuracy
  use dftbp_constants, only : Hartree__eV, au__Debye
  use dftbp_nonscc, only : TNonSccDiff
  use dftbp_scc, only : TScc
  use dftbp_blasroutines
  use dftbp_eigensolver
  use dftbp_lapackroutines
  use dftbp_message
  use dftbp_taggedoutput, only : TTaggedWriter, tagLabels
  use dftbp_sorting
  use dftbp_qm
  use dftbp_transcharges
  use dftbp_linresptypes
  use dftbp_degeneracyfind
  implicit none
  private

  public :: LinRespGrad_old

  !> Output files for results
  character(*), parameter :: transitionsOut = "TRA.DAT"
  character(*), parameter :: XplusYOut = "XplusY.DAT"
  character(*), parameter :: excitedCoefsOut = "COEF.DAT"
  character(*), parameter :: excitationsOut = "EXC.DAT"
  character(*), parameter :: transDipOut = "TDP.DAT"


  ! ARPACK related variables

  !> Tolerance for ARPACK solver.
  real(dp), parameter :: ARTOL = epsilon(1.0_rsp)

  !> Maximal allowed iteration in the ARPACK solver.
  integer, parameter :: MAX_AR_ITER = 300

  !> Names of output files
  character(*), parameter :: arpackOut = "ARPACK.DAT"
  character(*), parameter :: testArpackOut = "TEST_ARPACK.DAT"

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

contains

  !> This subroutine analytically calculates excitations and gradients of excited state energies
  !> based on Time Dependent DFRT
  subroutine LinRespGrad_old(tSpin, this, iAtomStart, grndEigVecs, grndEigVal, sccCalc, dq, coord0,&
      & SSqr, filling, species0, iNeighbour, img2CentCell, orb, tWriteTagged, fdTagged,&
      & taggedWriter, omega, allOmega, shift, skHamCont, skOverCont, excgrad, derivator, rhoSqr,&
      & occNatural, naturalOrbs)

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

    !> excitation energy of state nStat
    real(dp), intent(out) :: omega

    !> excitation energy of all states that have been solved
    real(dp), allocatable, intent(inout) :: allOmega(:)

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


    real(dp) :: Ssq(this%nExc)
    real(dp), allocatable :: gammaMat(:,:), snglPartTransDip(:,:)
    real(dp), allocatable :: stimc(:,:,:), wij(:)
    real(dp), allocatable :: dqex(:,:), sposz(:), osz(:), xpy(:), xmy(:), pc(:,:,:)
    real(dp), allocatable :: t(:,:,:), rhs(:), woo(:,:), wvv(:,:), wov(:)
    real(dp), allocatable :: evec(:,:), eval(:), transitionDipoles(:,:)
    integer, allocatable :: win(:), getij(:,:)

    !> array from pairs of single particles states to compound index - should replace with a more
    !> compact data structure in the cases where there are oscilator windows
    integer, allocatable :: iatrans(:,:,:)

    character, allocatable :: symmetries(:)

    integer :: mnvir, nxoo_max, nxvv_max
    integer, allocatable :: nocc_ud(:), nvir_ud(:)
    integer :: mHOMO, mLUMO
    integer :: nxov, nxov_ud(2), nxov_r, nxov_d, nxov_rd
    integer :: norb
    integer :: i, j, iSpin, isym, iLev, nStartLev, nEndLev
    integer :: nSpin
    character :: sym
    character(lc) :: tmpStr

    real(dp) :: energyThreshold

    integer :: nStat

    !> control variables
    logical :: tZVector, tCoeffs, tTradip

    !> printing data
    logical :: tMulliken

    !> should gradients be calculated
    logical :: tForces

    !> transition charges, either cached or evaluated on demand
    type(TTransCharges) :: transChrg


    if (withArpack) then

      ! ARPACK library variables
      ndigit = -3
      ! Output unit:
      logfil = this%fdArnoldi
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

    @:ASSERT(this%fdExc > 0)

    ! work out which data files are required, based on whether they have valid file IDs (>0)
    tMulliken = (this%fdMulliken > 0)
    tCoeffs = (this%fdCoeffs > 0)
    tTradip = (this%fdTradip > 0)

    if (tMulliken) then
      open(this%fdMulliken, file=excitedQOut,position="rewind", status="replace")
      close(this%fdMulliken)
      open(this%fdMulliken, file=excitedDipoleOut, position="rewind", status="replace")
      close(this%fdMulliken)
    end if

    @:ASSERT(this%fdArnoldi > 0)
    if (this%tArnoldi) then
      open(this%fdArnoldi, file=arpackOut, position="rewind", status="replace")
    end if

    nSpin = size(grndEigVal, dim=2)
    @:ASSERT(nSpin > 0 .and. nSpin <=2)

    norb = orb%nOrb

    @:ASSERT(present(excgrad) .eqv. present(shift))
    @:ASSERT(present(excgrad) .eqv. present(skHamCont))
    @:ASSERT(present(excgrad) .eqv. present(skOverCont))
    @:ASSERT(present(excgrad) .eqv. present(derivator))
  #:block DEBUG_CODE
    if (present(excgrad)) then
    @:ASSERT(present(rhoSqr))
    end if
  #:endblock DEBUG_CODE
    @:ASSERT(present(occNatural) .eqv. present(naturalOrbs))

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
    tZVector = tForces .or. tMulliken .or. tCoeffs .or. present(naturalOrbs) .or.&
        & this%tWriteDensityMatrix

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
    ALLOCATE(stimc(norb, norb, nSpin))
    ALLOCATE(wij(nxov))
    ALLOCATE(win(nxov))
    ALLOCATE(eval(this%nExc))
    ALLOCATE(getij(nxov, 3))
    ALLOCATE(transitionDipoles(this%nExc, 3))
    ALLOCATE(sposz(nxov))
    ALLOCATE(nocc_ud(nSpin))
    ALLOCATE(nvir_ud(nSpin))

    ! Overlap times wave function coefficients - most routines in DFTB+ use lower triangle (would
    ! remove the need to symmetrize the overlap and ground state density matrix in the main code if
    ! this could be used everywhere in these routines)
    do iSpin = 1, nSpin
      call symm(stimc(:,:,iSpin), "L", SSqr, grndEigVecs(:,:,iSpin))
    end do

    ! ground state Hubbard U softened coulombic interactions
    call sccCalc%getAtomicGammaMatrix(gammaMat, iNeighbour, img2CentCell)

    ! Oscillator strengths for exited states, when needed.
    ALLOCATE(osz(this%nExc))

    ! Find all single particle transitions and KS energy differences for cases that go from filled
    ! to empty states
    call getSPExcitations(grndEigVal, filling, wij, getij)

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

    ! dipole strength of transitions between K-S states
    call calcTransitionDipoles(coord0, win, nxov_ud(1), getij, iAtomStart, stimc, grndEigVecs,&
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
            & grndEigVal, getij)

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
                & grndEigVal, getij)
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

    call TTransCharges_init(transChrg, iAtomStart, stimc, grndEigVecs, nxov_rd, nxov_ud(1), getij,&
        & win, this%tCacheCharges)


    if (this%fdXplusY >  0) then
      open(this%fdXplusY, file=XplusYOut, position="rewind", status="replace")
    end if

    if(this%fdTrans>0) then
      open(this%fdTrans, file=transitionsOut, position="rewind", status="replace")
      write(this%fdTrans,*)
    endif

    ! Many-body transition dipole file to excited states
    if (this%fdTradip > 0) then
      open(this%fdTradip, file=transDipOut, position="rewind", status="replace")
      write(this%fdTradip,*)
      write(this%fdTradip,'(5x,a,5x,a,2x,a)') "#", 'w [eV]', 'Transition dipole (x,y,z) [Debye]'
      write(this%fdTradip,*)
      write(this%fdTradip,'(1x,60("="))')
      write(this%fdTradip,*)
    endif

    ! excitation energies
    open(this%fdExc, file=excitationsOut, position="rewind", status="replace")
    write(this%fdExc,*)
    if (tSpin) then
      write(this%fdExc,'(5x,a,7x,a,9x,a,9x,a,6x,a,4x,a)') 'w [eV]', 'Osc.Str.', 'Transition',&
          & 'Weight', 'KS [eV]','D<S*S>'
    else
      write(this%fdExc,'(5x,a,7x,a,9x,a,9x,a,6x,a,4x,a)') 'w [eV]','Osc.Str.', 'Transition',&
          & 'Weight', 'KS [eV]','Sym.'
    end if

    write(this%fdExc,*)
    write(this%fdExc,'(1x,80("="))')
    write(this%fdExc,*)

    ! single particle excitations (output file and tagged file if needed).  Was used for nxov_rd =
    ! size(wij), but now for just states that are actually included in the excitation calculation.
    call writeSPExcitations(wij, win, nxov_ud(1), getij, this%fdSPTrans, sposz, nxov_rd, tSpin)
    ALLOCATE(evec(nxov_rd, this%nExc))

    do isym = 1, size(symmetries)

      sym = symmetries(isym)
      if (withArpack) then
        call buildAndDiagExcMatrixArpack(tSpin, wij(:nxov_rd), sym, win, nxov_ud(1), nxov_rd,&
            & iAtomStart, stimc, grndEigVecs, filling, getij, gammaMat, species0, this%spinW,&
            & transChrg, this%fdArnoldiDiagnosis, eval, evec, this%onSiteMatrixElements, orb)
      else
        call error("No suitable eigensolver was compiled with this binary")
      end if

      ! Excitation oscillator strengths for resulting states
      call getOscillatorStrengths(sym, snglPartTransDip(1:nxov_rd,:), wij(:nxov_rd), eval, evec,&
          & filling, win, nxov_ud(1), getij, nstat, osz, tTradip, transitionDipoles)

      if (tSpin) then
        call getExcSpin(Ssq, nxov_ud(1), getij, win, eval, evec, wij(:nxov_rd), filling, stimc,&
            & grndEigVecs)
        call writeExcitations(sym, osz, this%nExc, nxov_ud(1), getij, win, eval, evec,&
            & wij(:nxov_rd), this%fdXplusY, this%fdTrans, this%fdTradip, transitionDipoles,&
            & tWriteTagged, fdTagged, taggedWriter, this%fdExc, Ssq)
      else
        call writeExcitations(sym, osz, this%nExc, nxov_ud(1), getij, win, eval, evec,&
            & wij(:nxov_rd), this%fdXplusY, this%fdTrans, this%fdTradip, transitionDipoles,&
            & tWriteTagged, fdTagged, taggedWriter, this%fdExc)
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

    if (this%tArnoldi) then
      close(this%fdArnoldi)
    end if

    if (this%fdTrans > 0) close(this%fdTrans)
    if (this%fdXplusY > 0) close(this%fdXplusY)
    if (this%fdExc > 0) close(this%fdExc)
    if (this%fdTradip > 0) close(this%fdTradip)

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

      mnvir = maxval(nvir_ud)
      nxoo_max = (mHOMO * (mHOMO + 1)) / 2
      nxvv_max = (mnvir * (mnvir + 1)) / 2

      ! Arrays needed for Z vector
      ALLOCATE(xpy(nxov_rd))
      ALLOCATE(xmy(nxov_rd))
      ALLOCATE(t(norb, norb, nSpin))
      ALLOCATE(rhs(nxov_rd))
      ALLOCATE(woo(nxoo_max, nSpin))
      ALLOCATE(wvv(nxvv_max, nSpin))
      ALLOCATE(wov(nxov_rd))
      ALLOCATE(iatrans(1:mHOMO, mLUMO:norb, nSpin))

      ! Arrays for gradients and Mulliken analysis
      if (tZVector) then
        ALLOCATE(dqex(this%nAtom, nSpin))
        ALLOCATE(pc(norb, norb, nSpin))
      end if

      ! set up transition indexing
      call rindxov_array(win, mLUMO, nxov, getij, iatrans)

      do iLev = nStartLev, nEndLev
        omega = sqrt(eval(iLev))
        ! Furche terms: X+Y, X-Y
        xpy(:nxov_rd) = evec(:nxov_rd,iLev) * sqrt(wij(:nxov_rd) / omega)
        xmy(:nxov_rd) = evec(:nxov_rd,iLev) * sqrt(omega / wij(:nxov_rd))

        ! solve for Z and W to get excited state density matrix
        call getZVectorEqRHS(xpy, xmy, win, iAtomStart, nocc_ud,&
            & nxov_ud(1), getij, iatrans, this%nAtom, species0,grndEigVal,&
            & stimc, grndEigVecs, gammaMat, this%spinW, omega, sym, rhs, t,&
            & wov, woo, wvv, transChrg)
        call solveZVectorPrecond(rhs, win, nxov_ud(1), getij, this%nAtom, iAtomStart,&
            & stimc, gammaMat, wij(:nxov_rd), grndEigVecs, transChrg, species0, this%spinW)
        call calcWVectorZ(rhs, win, nocc_ud, nxov_ud(1), getij, iAtomStart,&
            & stimc, grndEigVecs, gammaMat, grndEigVal, wov, woo, wvv, transChrg, species0, this%spinW)
        call calcPMatrix(t, rhs, win, getij, pc)

        call writeCoeffs(pc, grndEigVecs, filling, this%fdCoeffs,&
            & tCoeffs, this%tGrndState, occNatural, naturalOrbs)

        do iSpin = 1, nSpin
          ! Make MO to AO transformation of the excited density matrix
          call makeSimilarityTrans(pc(:,:,iSpin), grndEigVecs(:,:,iSpin))
          call getExcMulliken(iAtomStart, pc(:,:,iSpin), SSqr, dqex(:,iSpin))
        end do

        if (this%tWriteDensityMatrix) then
          call writeDM(iLev, pc, rhoSqr)
        end if

        if (tMulliken) then
          !> for now, only total Mulliken charges
          call writeExcMulliken(sym, iLev, dq(:,1), sum(dqex,dim=2), coord0, this%fdMulliken)
        end if

        if (tForces) then
          call addGradients(sym, nxov_rd, this%nAtom, species0, iAtomStart, norb, nocc_ud,&
              & getij, win, grndEigVecs, pc, stimc, dq, dqex, gammaMat, this%HubbardU,&
              & this%spinW, shift, woo, wov, wvv, transChrg, xpy, coord0, orb, skHamCont,&
              & skOverCont, derivator, rhoSqr, excgrad)
        end if

      end do

      if (nstat == 0) then
        omega = 0.0_dp
      end if

    end if

  end subroutine LinRespGrad_old


  !> Builds and diagonalizes the excitation matrix via iterative technique.
  subroutine buildAndDiagExcMatrixArpack(tSpin, wij, sym, win, nmatup, nxov, iAtomStart, stimc,&
      & grndEigVecs, filling, getij, gammaMat, species0, spinW, transChrg, fdArnoldiDiagnosis,&
      & eval, evec, onsMEs, orb)

    !> spin polarisation?
    logical, intent(in) :: tSpin

    !> single particle excitation energies
    real(dp), intent(in) :: wij(:)

    !> symmetry to calculate transitions
    character, intent(in) :: sym

    !> index array for single particle excitions
    integer, intent(in) :: win(:)

    !> number of same spin excitations
    integer, intent(in) :: nmatup

    !> number of occupied-virtual transitions
    integer, intent(in) :: nxov

    !> indexing array for square matrices
    integer, intent(in) :: iAtomStart(:)

    !> overlap times ground state eigenvectors
    real(dp), intent(in) :: stimc(:,:,:)

    !> ground state eigenvectors
    real(dp), intent(in) :: grndEigVecs(:,:,:)

    !> occupation numbers
    real(dp), intent(in) :: filling(:,:)

    !> electrostatic matrix
    real(dp), intent(in) :: gammaMat(:,:)

    !> index array for for single particle excitations
    integer, intent(in) :: getij(:,:)

    !> central cell chemical species
    integer, intent(in) :: species0(:)

    !> file handle for ARPACK eigenstate tests
    integer, intent(in) :: fdArnoldiDiagnosis

    !> atomic resolved spin constants
    real(dp), intent(in) :: spinW(:)

    !> machinery for transition charges between single particle levels
    type(TTransCharges), intent(in) :: transChrg

    !> resulting eigenvalues for transitions
    real(dp), intent(out) :: eval(:)

    !> eigenvectors for transitions
    real(dp), intent(out) :: evec(:,:)

    !> onsite corrections if in use
    real(dp), allocatable :: onsMEs(:,:,:,:)

    !> data type for atomic orbital information
    type(TOrbitals), intent(in) :: orb

    real(dp), allocatable :: workl(:), workd(:), resid(:), vv(:,:), qij(:)
    real(dp) :: sigma
    integer :: iparam(11), ipntr(11)
    integer :: ido, ncv, lworkl, info
    logical, allocatable :: selection(:)
    logical :: rvec
    integer :: nexc, natom

    integer :: iState
    real(dp), allocatable :: Hv(:), orthnorm(:,:)
    character(lc) :: tmpStr

    nexc = size(eval)
    natom = size(gammaMat, dim=1)
    @:ASSERT(all(shape(evec) == [ nxov, nexc ]))

    ! Three times more Lanczos vectors than desired eigenstates
    ncv = min(3 * nexc, nxov)

    lworkl = ncv * (ncv + 8)

    ALLOCATE(workl(lworkl))
    ALLOCATE(workd(3 * nxov))
    ALLOCATE(resid(nxov))
    ALLOCATE(selection(ncv))
    ALLOCATE(vv(nxov, ncv))
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
      call saupd (ido, "I", nxov, "SM", nexc, ARTOL, resid, ncv, vv, nxov, iparam, ipntr, workd,&
          & workl, lworkl, info)

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
      call omegatvec(tSpin, workd(ipntr(1):ipntr(1)+nxov-1), workd(ipntr(2):ipntr(2)+nxov-1),&
          & wij, sym, win, nmatup, iAtomStart, stimc, grndEigVecs, filling, getij, gammaMat,&
          & species0, spinW, onsMEs, orb, transChrg)

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
      call seupd (rvec, "All", selection, eval, evec, nxov, sigma, "I", nxov, "SM", nexc, ARTOL,&
          & resid, ncv, vv, nxov, iparam, ipntr, workd, workl, lworkl, info)

      ! check for error on return
      if (info  /=  0) then
        write(tmpStr,"(' Error with ARPACK routine seupd, info = ',I0)")info
        call error(tmpStr)
      end if

    end if

    if (fdArnoldiDiagnosis > 0) then
      ! tests for quality of returned eigenpairs
      open(fdArnoldiDiagnosis, file=testArpackOut, position="rewind", status="replace")
      ALLOCATE(Hv(nxov))
      ALLOCATE(orthnorm(nxov,nxov))
      orthnorm = matmul(transpose(evec(:,:nExc)),evec(:,:nExc))

      write(fdArnoldiDiagnosis,"(A)")'State Ei deviation    Evec deviation  Norm deviation  Max&
          & non-orthog'
      do iState = 1, nExc
        call omegatvec(tSpin, evec(:,iState), Hv, wij, sym, win, nmatup, iAtomStart, stimc,&
            & grndEigVecs, filling, getij, gammaMat, species0, spinW, onsMEs, orb, transChrg)
        write(fdArnoldiDiagnosis,"(I4,4E16.8)")iState, dot_product(Hv,evec(:,iState))-eval(iState),&
            & sqrt(sum( (Hv-evec(:,iState)*eval(iState) )**2 )), orthnorm(iState,iState) - 1.0_dp,&
            & max(maxval(orthnorm(:iState-1,iState)), maxval(orthnorm(iState+1:,iState)))
      end do
      close(fdArnoldiDiagnosis)
    end if

  end subroutine buildAndDiagExcMatrixArpack


  !> Calculate oscillator strength for a given excitation between KS states
  subroutine getOscillatorStrengths(sym, snglPartTransDip, wij, eval, evec, filling, win, nmatup,&
      & getij, istat, osz, tTradip, transitionDipoles)

    !> symmetry of transition
    character, intent(in) :: sym

    !> dipole moments for single particle transtions
    real(dp), intent(in) :: snglPartTransDip(:,:)

    !> energies for single particle transitions
    real(dp), intent(in) :: wij(:)

    !> Low lying eigenvalues of Casida eqn
    real(dp), intent(in) :: eval(:)

    !> eigenvectors of Casida eqn
    real(dp), intent(in) :: evec(:,:)

    !> Single particle occupations in the ground state
    real(dp), intent(in) :: filling(:,:)

    !> index for transitions
    integer, intent(in) :: win(:)

    !> number of up spin transitions before the down spin start
    integer, intent(in) :: nmatup

    !> index from single particle excitation to specific pair of single particle states involved
    integer, intent(in) :: getij(:,:)

    !> write transition dipole
    logical :: tTradip

    !> flag wich if <-1 on entry is returned as the brightest state
    integer, intent(inout) :: istat

    !> Oscilator strengths of transitions
    real(dp), intent(out) :: osz(:)

    !> resulting transition dipoles
    real(dp), intent(out) :: transitionDipoles(:,:)

    integer :: ii, nmat, oszLoc(1)
    real(dp) :: wnij(size(evec, dim=1))
    logical :: tSpin

    nmat = size(evec, dim=1)

    if (size(filling, dim=2) == 2) then
      tSpin = .true.
    else
      tSpin = .false.
    end if

    transitionDipoles(:,:) = 0.0_dp
    osz = 0.0_dp

    ! Triplet oscillator strength and transition dipole is zero for
    ! closed shell ground state
    if ((.not. tSpin) .and. (sym == "T")) then
      return
    end if

    call wtdn(wij, filling, win, nmatup, nmat, getij, wnij)

    !$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(ii) SCHEDULE(RUNTIME)
    do ii = 1, size(evec, dim=2)
      osz(ii) = oscillatorStrength(snglPartTransDip, wnij, evec(:,ii))
    end do
    !$OMP  END PARALLEL DO

    if (istat < 0) then
      ! find largest transition dipole transition
      oszLoc = maxloc(osz)
      istat = oszLoc(1)
    end if

    if (tTradip) then
      call transitionDipole(snglPartTransDip, wnij, eval, evec, transitionDipoles)
    end if

  end subroutine getOscillatorStrengths


  !> Build right hand side of the equation for the Z-vector and those parts of the W-vectors which
  !> do not depend on Z.
  subroutine getZVectorEqRHS(xpy, xmy, win, iAtomStart, homo, nmatup, getij, iatrans, natom,&
      & species0, grndEigVal, stimc, c, gammaMat, spinW, omega, sym, rhs, t, wov, woo, wvv,&
      & transChrg)

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

    !> index array between transitions in square and 1D representations
    integer, intent(in) :: getij(:,:)

    !> index array from orbital pairs to compound index
    integer, intent(in) :: iatrans(:,minval(homo)+1:,:)

    !> number of central cell atoms
    integer, intent(in) :: natom

    !> central cell chemical species
    integer, intent(in) :: species0(:)

    !> ground state wavefunctions
    real(dp), intent(in) :: grndEigVal(:,:)

    !> overlap times ground state wavefunctions
    real(dp), intent(in) :: stimc(:,:,:)

    !> ground state wavefunctions
    real(dp), intent(in) :: c(:,:,:)

    !> softened coulomb matrix
    real(dp), intent(in) :: gammaMat(:,:)

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

    !> machinery for transition charges between single particle levels
    type(TTransCharges), intent(in) :: transChrg

    real(dp), allocatable :: xpyq(:), qij(:), gamxpyq(:), qgamxpyq(:,:), gamqt(:)
    real(dp), allocatable :: xpyqds(:), gamxpyqds(:)
    integer :: nxov
    integer, allocatable :: nxoo(:), nxvv(:), nvir(:)
    integer :: i, j, a, b, ias, ibs, ij, ab, jas, s, nSpin
    real(dp) :: tmp1, tmp2, fact
    logical :: updwn, tSpin

    nxov = size(rhs)

    ALLOCATE(xpyq(natom))
    ALLOCATE(qij(natom))
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
      call indxov(win, ias, getij, i, a, s)

      ! BA: is T_aa = 0?
      do b = homo(s) + 1, a
        ibs = iatrans(i, b, s)
        call rindxvv(homo(s), a, b, ab)
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
        call rindxvv(0, j, i, ij)
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
    call transChrg%qMatVec(iAtomStart, stimc, c, getij, win, xpy, xpyq)

    if (.not. tSpin) then  ! ---- spin-unpolarized case ----

      ! qgamxpyq(ab) = sum_jc K_ab,jc (X+Y)_jc
      if (sym == "S") then
        call hemv(gamxpyq, gammaMat,  xpyq)
        do ab = 1, nxvv(1)
          call indxvv(homo(1), ab, a, b)
          qij(:) = transq(a, b, iAtomStart, .true., stimc, c)
          qgamxpyq(ab, 1) = 2.0_dp * sum(qij * gamxpyq)
        end do
      else ! triplet case
        do ab = 1, nxvv(1)
          call indxvv(homo(1), ab, a, b)
          qij(:) = transq(a, b, iAtomStart, .true., stimc, c)
          qgamxpyq(ab, 1) = 2.0_dp * sum(qij * xpyq * spinW(species0))
        end do
      end if

    else  ! ---- spin-polarized case -----

      xpyqds(:) = 0.0_dp
      call transChrg%qMatVecDs(iAtomStart, stimc, c, getij, win, xpy, xpyqds)

      call hemv(gamxpyq, gammaMat,  xpyq)
      do s = 1, 2
        if (s == 1) then
          updwn = .true.
          fact = 1.0_dp
        else
          updwn = .false.
          fact = -1.0_dp
        end if
        do ab = 1, nxvv(s)
          call indxvv(homo(s), ab, a, b)
          qij(:) = transq(a, b, iAtomStart, updwn, stimc, c)
          qgamxpyq(ab, s) = sum(qij * gamxpyq)
          !magnetization part
          qgamxpyq(ab, s) = qgamxpyq(ab, s) + fact * sum(qij * xpyqds * spinW(species0))
        end do
      end do

    end if

    ! rhs(ia) -= Qia = sum_b (X+Y)_ib * qgamxpyq(ab))
    do ias = 1, nxov
      call indxov(win, ias, getij, i, a, s)

      do b = homo(s) + 1, a
        call rindxvv(homo(s), a, b, ab)
        ibs = iatrans(i,b,s)
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
          call indxoo(ij, i, j)
          qij(:) = transq(i, j, iAtomStart, .true., stimc, c)
          ! qgamxpyq(ij) = sum_kb K_ij,kb (X+Y)_kb
          qgamxpyq(ij, 1) = 2.0_dp * sum(qij * gamxpyq)
        end do
      else
        do ij = 1, nxoo(1)
          qgamxpyq(ij, 1) = 0.0_dp
          call indxoo(ij, i, j)
          qij(:) = transq(i, j, iAtomStart, .true., stimc, c)
          qgamxpyq(ij, 1) = 2.0_dp * sum(qij * xpyq * spinW(species0))
        end do
      end if

    else  ! ---- spin-polarized case -----

      do s = 1, 2
        if (s == 1) then
          updwn = .true.
          fact = 1.0_dp
        else
          updwn = .false.
          fact = -1.0_dp
        end if
        do ij = 1, nxoo(s)
          qgamxpyq(ij, s) = 0.0_dp
          call indxoo(ij, i, j)
          qij(:) = transq(i, j, iAtomStart, updwn, stimc, c)
          qgamxpyq(ij, s) = sum(qij * gamxpyq)
          !magnetization part
          qgamxpyq(ij, s) = qgamxpyq(ij, s) + fact * sum(qij * xpyqds * spinW(species0))
        end do
      end do

    end if

    ! rhs(ia) += Qai = sum_j (X+Y)_ja qgamxpyq(ij)
    ! add Qai to Wia as well.
    do ias = 1, nxov
      call indxov(win, ias, getij, i, a, s)
      do j = i, homo(s)
        jas = iatrans(j, a, s)
        !ij = i-homo+nocc + ((j-homo+nocc - 1) * (j-homo+nocc)) / 2
        call rindxvv(0, j, i, ij)
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
        updwn = .true.
        fact = 1.0_dp
      else
        updwn = .false.
        fact = -1.0_dp
      end if
      do ij = 1, nxoo(s)
        call indxoo(ij, i, j)
        qij = transq(i, j, iAtomStart, updwn, stimc, c)
        if (i == j) then
          gamxpyq(:) = gamxpyq(:) + t(i,j,s) * qij(:)
          if (tSpin) then
            gamxpyqds(:) = gamxpyqds(:) + t(i,j,s) * qij(:) * fact
          end if
        else
          ! factor 2 because of symmetry of the matrix
          gamxpyq(:) = gamxpyq(:) + 2.0_dp  * t(i,j,s) * qij(:)
          if (tSpin) then
            gamxpyqds(:) = gamxpyqds(:) + 2.0_dp * t(i,j,s) * qij(:) * fact
          end if
        end if
      end do

      ! gamxpyq(iAt2) += sum_ab q_ab(iAt2) T_ab
      do ab = 1, nxvv(s)
        call indxvv(homo(s), ab, a, b)
        qij(:) = transq(a, b, iAtomStart, updwn, stimc, c)
        if (a == b) then
          gamxpyq(:) = gamxpyq(:) + t(a,b,s) * qij(:)
          if (tSpin) then
            gamxpyqds(:) = gamxpyqds(:) + t(a,b,s) * qij(:) * fact
          end if
        else
          ! factor 2 because of symmetry of the matrix
          gamxpyq(:) = gamxpyq(:) + 2.0_dp * t(a,b,s) * qij(:)
          if (tSpin) then
            gamxpyqds(:) = gamxpyqds(:) + 2.0_dp * t(a,b,s) * qij(:) * fact
          end if
        end if
      end do

    end do

    ! gamqt(iAt1) = sum_iAt2 gamma_iAt1,iAt2 gamxpyq(iAt2)
    call hemv(gamqt, gammaMat, gamxpyq)

    ! rhs -= sum_q^ia(iAt1) gamxpyq(iAt1)
    if (.not. tSpin) then
      call transChrg%qVecMat(iAtomStart, stimc, c, getij, win, -4.0_dp*gamqt, rhs)
    else
      call transChrg%qVecMat(iAtomStart, stimc, c, getij, win, -2.0_dp*gamqt, rhs)
      call transChrg%qVecMatDs(iAtomStart, stimc, c, getij, win, -2.0_dp*gamxpyqds*spinW(species0), rhs)
    end if

    ! Furche vectors
    do s = 1, nSpin
      if (s == 1) then
        updwn = .true.
        fact = 1.0_dp
      else
        updwn = .false.
        fact = -1.0_dp
      end if
      do ij = 1, nxoo(s)
        call indxoo(ij, i, j)
        qij(:) = transq(i, j, iAtomStart, updwn, stimc, c)
        if (.not. tSpin) then
          woo(ij,s) = woo(ij,s) + 4.0_dp * sum(qij * gamqt)
        else
          woo(ij,s) = woo(ij,s) + 2.0_dp * sum(qij * gamqt)
          woo(ij,s) = woo(ij,s) + 2.0_dp * fact * sum(qij * gamxpyqds * spinW(species0))
        end if
      end do
    end do

  end subroutine getZVectorEqRHS


  !> Solving the (A+B) Z = -R equation via diagonally preconditioned conjugate gradient
  subroutine solveZVectorPrecond(rhs, win, nmatup, getij, natom, iAtomStart, stimc, gammaMat, wij,&
      & c, transChrg, species0, spinW)

    !> on entry -R, on exit Z
    real(dp), intent(inout) :: rhs(:)

    !> index for single particle excitations
    integer, intent(in) :: win(:)

    !> number of transitions between only up states
    integer, intent(in) :: nmatup

    !> index array from composite index to specific filled-empty transition
    integer, intent(in) :: getij(:,:)

    !> number of atoms
    integer, intent(in) :: natom

    !> index vector for S and H0 matrices
    integer, intent(in) :: iAtomStart(:)

    !> overlap times ground state mo-coefficients
    real(dp), intent(in) :: stimc(:,:,:)

    !> Softened coulomb matrix
    real(dp), intent(in) :: gammaMat(:,:)

    !> single particle excitation energies
    real(dp), intent(in) :: wij(:)

    !> ground state mo-coefficients
    real(dp), intent(in) :: c(:,:,:)

    !> machinery for transition charges between single particle levels
    type(TTransCharges), intent(in) :: transChrg

    !> central cell chemical species
    integer, intent(in) :: species0(:)

    !> ground state spin derivatives for each species
    real(dp), intent(in) :: spinW(:)

    integer :: nxov
    integer :: ia, kk
    real(dp) :: rhs2(size(rhs)), rkm1(size(rhs)), zkm1(size(rhs)), pkm1(size(rhs)), apk(size(rhs))
    real(dp) :: qTmp(nAtom), rs, alphakm1, tmp1, tmp2, bkm1
    real(dp), allocatable :: qij(:), P(:)
    logical :: tSpin

    nxov = size(rhs)
    allocate(qij(nAtom))

    tSpin = (nxov > nmatup)

    ! diagonal preconditioner
    ! P^-1 = 1 / (A+B)_ia,ia (diagonal of the supermatrix sum A+B)
    allocate(P(nxov))
    do ia = 1, nxov
      qij = transChrg%qTransIJ(ia, iAtomStart, stimc, c, getij, win)
      call hemv(qTmp, gammaMat, qij)
      if (.not. tSpin) then
        rs = 4.0_dp * dot_product(qij, qTmp) + wij(ia)
      else
        rs = 2.0_dp * dot_product(qij, qTmp) + wij(ia)
        rs = rs + 2.0_dp * sum(qij * qij * spinW(species0))
      end if
      P(ia) = 1.0_dp / rs
    end do

    ! Free some space, before entering the apbw routine
    deallocate(qij)

    ! unit vector as initial guess solution
    rhs2(:) = 1.0_dp / sqrt(real(nxov,dp))

    ! action of matrix on vector
    call apbw(rkm1, rhs2, wij, nxov, natom, win, nmatup, getij, iAtomStart, stimc, c, gammaMat,&
        & transChrg, species0, spinW)

    rkm1(:) = rhs - rkm1
    zkm1(:) = P * rkm1
    pkm1(:) = zkm1

    ! Iteration: should be convergent in at most nxov steps for a quadradic surface, so set higher
    do kk = 1, nxov**2

      ! action of matrix on vector
      call apbw(apk, pkm1, wij, nxov, natom, win, nmatup, getij, iAtomStart, stimc, c, gammaMat,&
          & transChrg, species0, spinW)

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


  !> Calculate Z-dependent parts of the W-vectors and divide diagonal elements of W_ij and W_ab by
  !> 2.
  subroutine calcWvectorZ(zz, win, homo, nmatup, getij, iAtomStart, stimc, c, gammaMat,&
      & grndEigVal, wov, woo, wvv, transChrg, species0, spinW)

    !> Z vector
    real(dp), intent(in) :: zz(:)

    !> index array for single particle transitions
    integer, intent(in) :: win(:)

    !> highest occupied level
    integer, intent(in) :: homo(:)

    !> number of same spin excitations
    integer, intent(in) :: nmatup

    !> index array between transitions in square and 1D representations
    integer, intent(in) :: getij(:,:)

    !> index array for S and H0 ground state square matrices
    integer, intent(in) :: iAtomStart(:)

    !> overlap times ground state wavefunctions
    real(dp), intent(in) :: stimc(:,:,:)

    !> ground state wavefunctions
    real(dp), intent(in) :: c(:,:,:)

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

    integer :: nxov, natom, nSpin
    integer, allocatable :: nxoo(:), nxvv(:), nvir(:)
    integer :: ij, ias, ab, i, j, a, b, s, iAt1
    real(dp) :: fact
    real(dp), allocatable :: qij(:), gamxpyq(:), zq(:), zqds(:)
    logical :: updwn, tSpin

    nxov = size(zz)
    natom = size(gammaMat, dim=1)
    nSpin = size(grndEigVal, dim=2)

    ALLOCATE(qij(natom))
    ALLOCATE(gamxpyq(natom))
    ALLOCATE(zq(natom))
    ALLOCATE(nxoo(nSpin))
    ALLOCATE(nxvv(nSpin))
    ALLOCATE(nvir(nSpin))

    nxoo(:) = (homo(:)*(homo(:)+1))/2
    nvir(:) = size(c, dim=1) - homo(:)
    nxvv(:) = (nvir(:)*(nvir(:)+1))/2

    if ( nSpin == 2 ) then
      tSpin = .true.
      ALLOCATE(zqds(natom))
    else
      tSpin = .false.
    end if

    ! Adding missing epsilon_i * Z_ia term to W_ia
    do ias = 1, nxov
      call indxov(win, ias, getij, i, a, s)
      wov(ias) = wov(ias) + zz(ias) * grndEigVal(i, s)
    end do

    ! Missing sum_kb 4 K_ijkb Z_kb term in W_ij: zq(iAt1) = sum_kb q^kb(iAt1) Z_kb
    zq(:) = 0.0_dp
    call transChrg%qMatVec(iAtomStart, stimc, c, getij, win, zz, zq)
    call hemv(gamxpyq, gammaMat, zq)

    if (tSpin) then
      zqds(:) = 0.0_dp
      call transChrg%qMatVecDs(iAtomStart, stimc, c, getij, win, zz, zqds)
    end if

    ! sum_iAt1 qij(iAt1) gamxpyq(iAt1)
    do s = 1, nSpin
      if (s == 1) then
        updwn = .true.
        fact = 1.0_dp
      else
        updwn = .false.
        fact = -1.0_dp
      end if
      do ij = 1, nxoo(s)
        call indxoo(ij, i, j)
        qij(:) = transq(i, j, iAtomStart, updwn, stimc, c)
        ! W contains 1/2 for i == j.
        if (.not. tSpin) then
          woo(ij,s) = woo(ij,s) + 4.0_dp * sum(qij * gamxpyq)
        else
          woo(ij,s) = woo(ij,s) + 2.0_dp * sum(qij * gamxpyq)
          woo(ij,s) = woo(ij,s) + 2.0_dp * fact * sum(qij * zqds * spinW(species0))
        end if
      end do
    end do

    ! Divide diagonal elements of W_ij by 2.
    do s = 1, nSpin
      do ij = 1, nxoo(s)
        call indxoo(ij, i, j)
        if (i == j) then
          woo(ij,s) = 0.5_dp * woo(ij,s)
        end if
      end do
    end do

    ! Divide diagonal elements of W_ab by 2.
    do s = 1, nSpin
      do ab = 1, nxvv(s)
        call indxvv(homo(s), ab, a, b)
        if (a == b) then
          wvv(ab,s) = 0.5_dp * wvv(ab,s)
        end if
      end do
    end do

  end subroutine calcWvectorZ


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
  subroutine addGradients(sym, nxov, natom, species0, iAtomStart, norb, homo, getij,&
      & win, grndEigVecs, pc, stimc, dq_ud, dqex, gammaMat, HubbardU, spinW, shift, woo, wov, wvv,&
      & transChrg, xpy, coord0, orb, skHamCont, skOverCont, derivator, rhoSqr, excgrad)

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

    !> index array from composite transition index to specific single particle states
    integer, intent(in) :: getij(:,:)

    !> single particle transition index
    integer, intent(in) :: win(:)

    !> ground state eigenvectors
    real(dp), intent(in) :: grndEigVecs(:,:,:)

    !> transition density matrix
    real(dp), intent(in) :: pc(:,:,:)

    !> overlap times ground state eigenvectors
    real(dp), intent(in) :: stimc(:,:,:)

    !> ground state gross charges
    real(dp), intent(in) :: dq_ud(:,:)

    !> charge differences from ground to excited state
    real(dp), intent(in) :: dqex(:,:)

    !> softened coulomb matrix
    real(dp), intent(in) :: gammaMat(:,:)

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

    !> central cell atomic coordinates
    real(dp), intent(in) :: coord0(:,:)

    !> data type for atomic orbital information
    type(TOrbitals), intent(in) :: orb

    !> H0 data
    type(TSlakoCont), intent(in) :: skHamCont

    !> overlap data
    type(TSlakoCont), intent(in) :: skOverCont

    !> Differentiatior for the non-scc matrices
    class(TNonSccDiff), intent(in) :: derivator

    !> ground state density matrix
    real(dp), intent(in) :: rhoSqr(:,:,:)

    !> resulting excited state gradient
    real(dp), intent(out) :: excgrad(:,:)

    real(dp), allocatable :: shift_excited(:,:), xpyq(:), xpyqds(:)
    real(dp), allocatable :: shxpyq(:,:), xpycc(:,:,:), wcc(:,:,:), tmp5(:), tmp7(:), tmp11(:)
    real(dp), allocatable :: qij(:), temp(:), dq(:), dm(:), dsigma(:)
    real(dp), allocatable :: dH0(:,:,:), dS(:,:,:)
    real(dp), allocatable :: Dens(:,:), SpinDens(:,:)
    integer :: ia, i, j, a, b, ab, ij, m, n, mu, nu, xyz, iAt1, iAt2, s
    integer :: indalpha, indalpha1, indbeta, indbeta1
    integer :: iSp1, iSp2, iSpin, nSpin
    real(dp) :: tmp1, tmp2, tmp3, tmp4, tmp6, tmp8, tmp9, tmp10, rab
    real(dp) :: diffvec(3), dgab(3), tmp3a, tmp3b
    integer, allocatable :: nxoo(:), nxvv(:), nvir(:)
    logical :: tSpin

    nSpin = size(grndEigVecs, dim=3)
    tSpin = (nSpin == 2)

    ALLOCATE(shift_excited(natom, nSpin))
    ALLOCATE(xpyq(natom))
    ALLOCATE(shxpyq(natom, nSpin))
    ALLOCATE(xpycc(norb, norb, nSpin))
    ALLOCATE(wcc(norb, norb, nSpin))
    ALLOCATE(qij(natom))
    ALLOCATE(temp(norb))
    ALLOCATE(tmp5(nSpin))
    ALLOCATE(tmp7(nSpin))

    ALLOCATE(Dens(norb,norb))
    Dens(:,:) = sum(rhoSqr, dim=3)

    ALLOCATE(dH0(orb%mOrb, orb%mOrb, 3))
    ALLOCATE(dS(orb%mOrb, orb%mOrb, 3))

    ALLOCATE(nxoo(nSpin))
    ALLOCATE(nxvv(nSpin))
    ALLOCATE(nvir(nSpin))

    nxoo(:) = (homo(:)*(homo(:)+1))/2
    nvir(:) = norb - homo(:)
    nxvv(:) = (nvir(:)*(nvir(:)+1))/2

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

    excgrad = 0.0_dp

    ! excited state potentials at atomic sites
    do iSpin = 1, nSpin
      call hemv(shift_excited(:,iSpin), gammaMat, dqex(:,iSpin))
    end do

    ! xypq(alpha) = sum_ia (X+Y)_ia q^ia(alpha)
    ! complexity norb * norb * norb
    xpyq(:) = 0.0_dp
    call transChrg%qMatVec(iAtomStart, stimc, grndEigVecs, getij, win, xpy,xpyq)

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
      call transChrg%qMatVecDs(iAtomStart, stimc, grndEigVecs, getij, win, xpy, xpyqds)
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
      call indxov(win, ia, getij, i, a, iSpin)
      ! should replace with DSYR2 call :
      do nu = 1, norb
        do mu = 1, norb
          xpycc(mu,nu,iSpin) = xpycc(mu,nu,iSpin) + xpy(ia) *&
              & ( grndEigVecs(mu,i,iSpin)*grndEigVecs(nu,a,iSpin)&
              & + grndEigVecs(mu,a,iSpin)*grndEigVecs(nu,i,iSpin) )
        end do
      end do
    end do

    ! calculate wcc = c_mu,i * W_ij * c_j,nu. We have only W_ab b > a and W_ij j > i:
    ! wcc(m,n) = sum_{pq, p <= q} w_pq (grndEigVecs(mu,p)grndEigVecs(nu,q)
    ! + grndEigVecs(nu,p)grndEigVecs(mu,q))
    ! complexity norb * norb * norb

    ! calculate the occ-occ part
    wcc(:,:,:) = 0.0_dp

    do iSpin = 1, nSpin
      do ij = 1, nxoo(iSpin)
        call indxoo(ij, i, j)
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
      call indxov(win, ia, getij, i, a, iSpin)
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
      do ab =1, nxvv(iSpin)
        call indxvv(homo(iSpin), ab, a, b)
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

        call derivator%getFirstDeriv(dH0, skHamCont, coord0, species0,&
            & iAt1, iAt2, orb)
        call derivator%getFirstDeriv(dS, skOverCont, coord0, species0,&
            & iAt1, iAt2, orb)

        do xyz = 1, 3

          tmp1 = 0.0_dp
          tmp2 = 0.0_dp
          tmp3 = 0.0_dp
          tmp4 = 0.0_dp
          tmp6 = 0.0_dp
          tmp8 = 0.0_dp
          tmp10 = 0.0_dp

          do iSpin = 1, nSpin

            do mu = indalpha, indalpha1
              do nu = indbeta, indbeta1
                m = mu - indalpha + 1
                n = nu - indbeta + 1

                tmp1 = tmp1 + 2.0_dp * dH0(n,m,xyz) * pc(mu,nu,iSpin)
                tmp2 = tmp2 + dS(n,m,xyz) * pc(mu,nu,iSpin) * (shift(iAt1)+shift(iAt2))
                tmp3 = tmp3 - dS(n,m,xyz) * wcc(mu,nu,iSpin)
                tmp4 = tmp4 + tmp5(iSpin) * dS(n,m,xyz) * Dens(mu,nu)
                tmp6 = tmp6 + tmp7(iSpin) * dS(n,m,xyz) * xpycc(mu,nu,iSpin)

                if (tSpin) then
                  tmp8 = tmp8 + tmp9 * dS(n,m,xyz) * dsigma(iSpin) * pc(mu,nu,iSpin)
                  tmp10 = tmp10 + tmp11(iSpin) * dS(n,m,xyz) * dsigma(iSpin) * SpinDens(mu,nu)
                end if

              end do
            end do

          end do

          excgrad(xyz,iAt1) = excgrad(xyz,iAt1)&
              & + tmp1 + tmp2 + tmp4 + tmp6 + tmp3 + tmp8 + tmp10
          excgrad(xyz,iAt2) = excgrad(xyz,iAt2)&
              & - tmp1 - tmp2 - tmp4 - tmp6 - tmp3 - tmp8 - tmp10
        end do
      end do
    end do

  end subroutine addGradients


  !> Write out excitations projected onto ground state
  subroutine writeCoeffs(tt, grndEigVecs, occ, fdCoeffs, tCoeffs, tIncGroundState,&
      & occNatural, naturalOrbs)

    !> T part of the matrix
    real(dp), intent(in) :: tt(:,:,:)

    !> ground state eigenvectors
    real(dp), intent(in) :: grndEigVecs(:,:,:)

    !> ground state occupations
    real(dp), intent(in) :: occ(:,:)

    !> file descriptor to write data into
    integer, intent(in) :: fdCoeffs

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
        open(fdCoeffs, file=excitedCoefsOut, position="append")
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
  subroutine writeExcitations(sym, osz, nexc, nmatup, getij, win, eval, evec, wij, fdXplusY,&
      & fdTrans, fdTradip, transitionDipoles, tWriteTagged, fdTagged, taggedWriter, fdExc, Ssq)

    !> Symmetry label for the type of transition
    character, intent(in) :: sym

    !> oscillator strengths for transitions from ground to excited states
    real(dp), intent(in) :: osz(:)

    !> number of excited states to solve for
    integer, intent(in) :: nexc

    !> number of same spin excitations
    integer, intent(in) :: nmatup

    !> index array between transitions in square and 1D representations
    integer, intent(in) :: getij(:,:)

    !> index array for single particle excitions
    integer, intent(in) :: win(:)

    !> excitation energies
    real(dp), intent(in) :: eval(:)

    !> eigenvectors of excited states
    real(dp), intent(in) :: evec(:,:)

    !> single particle excitation energies
    real(dp), intent(in) :: wij(:)

    !> single particle transition dipole moments
    real(dp), intent(in) :: transitionDipoles(:,:)

    !> should tagged information be written out
    logical, intent(in) :: tWriteTagged

    !> file unit for transition dipoles
    integer, intent(in) :: fdTradip

    !> file unit for X+Y data
    integer, intent(in) :: fdXplusY

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
    real(dp), allocatable :: xply(:)
    integer, allocatable :: wvin(:)
    real(dp) :: rsqw, weight, wvnorm
    logical :: updwn, tSpin
    character :: sign
    type(TDegeneracyFind) :: DegeneracyFind
    logical :: tDegenerate
    integer, allocatable :: degenerate(:,:)
    real(dp), allocatable :: oDeg(:)

    @:ASSERT(fdExc > 0)

    tSpin = present(Ssq)
    nmat = size(wij)

    allocate(wvec(nmat))
    allocate(wvin(nmat))
    allocate(xply(nmat))
    wvec(:) = 0.0_dp
    wvin(:) = 0
    xply(:) = 0.0_dp

    if(fdXplusY > 0) then
      write(fdXplusY,*) nmat, nexc
    end if

    do ii = 1, nexc
      if (eval(ii) > 0.0_dp) then

        ! calculate weight of single particle transitions
        rsqw = 1.0_dp / sqrt(eval(ii))
        ! (X+Y)^ia_I = sqrt(wij) / sqrt(omega) * F^ia_I
        xply(:) = sqrt(rsqw) * sqrt(wij(:)) * evec(:,ii)
        wvec(:) = xply**2
        wvnorm = 1.0_dp / sqrt(sum(wvec**2))
        wvec(:) = wvec * wvnorm

        ! find largest coefficient in CI - should use maxloc
        call index_heap_sort(wvin,wvec)
        wvin = wvin(size(wvin):1:-1)
        wvec = wvec(wvin)

        weight = wvec(1)
        iweight = wvin(1)

        call indxov(win, iweight, getij, m, n, s)
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

        if(fdXplusY > 0) then
          if (tSpin) then
            updwn = (win(iweight) <= nmatup)
            sign = "D"
            if (updwn) sign = "U"
          end if
          write(fdXplusY,'(1x,i5,3x,a,3x,ES17.10)') ii, sign, sqrt(eval(ii))
          write(fdXplusY,'(6(1x,ES17.10))') xply
        endif

        if (fdTrans > 0) then
          write(fdTrans, '(2x,a,T12,i5,T21,ES17.10,1x,a,2x,a)')&
              & 'Energy ', ii,  Hartree__eV * sqrt(eval(ii)), 'eV', sign
          write(fdTrans,*)
          write(fdTrans,'(2x,a,9x,a,8x,a)')'Transition', 'Weight', 'KS [eV]'
          write(fdTrans,'(1x,45("="))')

          sign = " "
          do jj = 1, nmat
            !if (wvec(jj) < 1e-4_dp) exit ! ??????
            indo = wvin(jj)
            call indxov(win, indo, getij, m, n, s)
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

        if (fdTradip > 0) then
          write(fdTradip, '(1x,i5,1x,f10.3,2x,3(ES14.6))')&
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
        call indxov(win, iWeight, getij, m, n, s)
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

        if(fdXplusY > 0) then
          if (tSpin) then
            updwn = (win(iweight) <= nmatup)
            sign = "D"
            if (updwn) sign = "U"
          end if
          write(fdXplusY,'(1x,i5,3x,a,3x,A)') ii,sign, '-'
        endif

        if (fdTrans > 0) then
          write(fdTrans, '(2x,a,1x,i5,5x,a,1x,a,3x,a)')&
              & 'Energy ', ii,  '-', 'eV', sign
          write(fdTrans,*)
        end if

        if(fdTradip > 0) then
          write(fdTradip, '(1x,i5,1x,A)') ii, '-'
        endif

      end if

    end do

    deallocate(wvec)
    deallocate(wvin)
    deallocate(xply)

    if (tWriteTagged) then

      call DegeneracyFind%init(elecTolMax)
      call DegeneracyFind%degeneracyTest(eval, tDegenerate)
      if (.not.tDegenerate) then
        call taggedWriter%write(fdTagged, tagLabels%excEgy, eval)
        call taggedWriter%write(fdTagged, tagLabels%excOsc, osz)
        if (fdTradip > 0) then
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
        if (fdTradip > 0) then
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
  subroutine calcPMatrix(t, rhs, win, getij, pc)

    !> T matrix
    real(dp), intent(in) :: t(:,:,:)

    !> Z matrix
    real(dp), intent(in) :: rhs(:)

    !> index array for single particle transitions
    integer, intent(in) :: win(:)

    !> array of the occupied->virtual pairs (nTransitions,occ 1 or virtual 2)
    integer, intent(in) :: getij(:,:)

    !> resulting excited state density matrix
    real(dp), intent(out) :: pc(:,:,:)

    integer :: ias, i, a, s, nSpin

    nSpin = size(pc, dim=3)

    pc = 0.0_dp
    do ias = 1, size(rhs)
      call indxov(win, ias, getij, i, a, s)
      pc(i,a,s) = rhs(ias)
    end do

    do s = 1, nSpin
      pc(:,:,s) = 0.5_dp * ( pc(:,:,s) + transpose(pc(:,:,s)) )
    end do

    pc = pc + t

  end subroutine calcPMatrix


end module dftbp_linrespgrad
