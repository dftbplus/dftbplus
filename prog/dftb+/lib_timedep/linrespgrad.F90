!--------------------------------------------------------------------------------------------------!
!  DFTB+: general package for performing fast atomistic simulations                                !
!  Copyright (C) 2018  DFTB+ developers group                                                      !
!                                                                                                  !
!  See the LICENSE file for terms of usage and distribution.                                       !
!--------------------------------------------------------------------------------------------------!

#:include 'common.fypp'

!> Linear response excitations and gradients with respect to atomic coordinates
module linrespgrad
  use assert
  use arpack
  use linrespcommon
  use commontypes
  use slakocont
  use shortgamma
  use accuracy
  use constants, only : Hartree__eV, au__Debye
  use nonscc, only : NonSccDiff
  use scc, only : TScc
  use blasroutines
  use eigensolver
  use message
  use taggedoutput
  use sorting
  use qm
  implicit none
  private

  public :: LinRespGrad_old


  !> Tolerance for ARPACK solver.
  real(dp), parameter :: ARTOL = epsilon(1.0_dp)


  !> Maximal allowed iteration in the ARPACK solver.
  integer, parameter :: MAX_AR_ITER = 300

  character(lc) :: tmpStr


  !> Names of output files
  character(*), parameter :: arpackOut = "ARPACK.DAT"
  character(*), parameter :: testArpackOut = "TEST_ARPACK.DAT"
  character(*), parameter :: transitionsOut = "TRA.DAT"
  character(*), parameter :: XplusYOut = "XplusY.DAT"
  character(*), parameter :: excitedQOut = "XCH.DAT"
  character(*), parameter :: excitedDipoleOut = "XREST.DAT"
  character(*), parameter :: excitedCoefsOut = "COEF.DAT"
  character(*), parameter :: excitationsOut = "EXC.DAT"
  character(*), parameter :: transDipOut = "TDP.DAT"
  character(*), parameter :: singlePartOut = "SPX.DAT"


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
  subroutine LinRespGrad_old(tSpin, natom, iAtomStart, grndEigVecs, grndEigVal, sccCalc, dq,&
      & coord0, nexc, nstat0, symc, SSqr, filling, species0, HubbardU, spinW, rnel, iNeighbour,&
      & img2CentCell, orb, tWriteTagged, fdTagged, fdMulliken, fdCoeffs, tGrndState, fdXplusY,&
      & fdTrans, fdSPTrans, fdTradip, tArnoldi, fdArnoldi, fdArnoldiDiagnosis, fdExc,&
      & tEnergyWindow, energyWindow,tOscillatorWindow, oscillatorWindow, omega, shift, skHamCont,&
      & skOverCont, excgrad, derivator, rhoSqr, occNatural, naturalOrbs)

    !> spin polarized calculation
    logical, intent(in) :: tSpin

    !> number of atoms
    integer, intent(in) :: natom

    !> index vector for S and H matrices
    integer, intent(in) :: iAtomStart(:)

    !> ground state MO-coefficients
    real(dp), intent(in) :: grndEigVecs(:,:,:)

    !> ground state MO-energies
    real(dp), intent(in) :: grndEigVal(:,:)

    !> Self-consistent charge module settings
    type(TScc), intent(in) :: sccCalc

    !> converged ground state Mulliken gross charges - atomic charges
    real(dp), intent(in) :: dq(:)

    !> atomic positions
    real(dp), intent(in) :: coord0(:,:)

    !> number of excited states to solve for
    integer, intent(in) :: nexc

    !> state of interest (< 0 find brightest, 0 calculate all nexc states, > 0 that specific state)
    integer, intent(in) :: nstat0

    !> symmetry required singlet ('S'), triplet ("T") or both ("B")
    character, intent(in) :: symc

    !> square overlap matrix between basis functions, both triangles required
    real(dp), intent(in) :: SSqr(:,:)

    !> occupations for the states
    real(dp), intent(in) :: filling(:,:)

    !> chemical species of each atom
    integer, intent(in) :: species0(:)

    !> ground state Hubbard U values for each species
    real(dp), intent(in) :: HubbardU(:)

    !> ground state spin derivatives for each species
    real(dp), intent(in) :: spinW(:)

    !> real number of electrons in system
    real(dp), intent(in) :: rnel

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

    !> file unit for excited Mulliken populations?
    integer, intent(in) :: fdMulliken

    !> file unit if the coefficients for the excited states should be written to disc
    integer, intent(in) :: fdCoeffs

    !> Add the ground state to the excited state transition density matrix when determining the
    !> natural orbitals
    logical, intent(in) :: tGrndState

    !> file for X+Y data
    integer, intent(in) :: fdXplusY

    !> File unit for single particle (KS) transitions if required
    integer, intent(in) :: fdTrans

    !> File unit for single particle transition dipole strengths
    integer, intent(in) :: fdSPTrans

    !> File unit for transition dipole data
    integer, intent(in) :: fdTradip

    !> write state of Arnoldi solver to disc
    logical, intent(in) :: tArnoldi

    !> file unit for Arnoldi write out
    integer, intent(in) :: fdArnoldi

    !> file unit for Arnoldi solver tests, if this is < 1 no tests are performed
    integer, intent(in) :: fdArnoldiDiagnosis

    !> file handle for excitation energies
    integer, intent(in) :: fdExc

    !> is an energy window specified
    logical, intent(in) :: tEnergyWindow

    !> energy window for transitions above energy of nexc-th single particle transtion
    real(dp), intent(in) :: energyWindow

    !> is an oscillator window specified
    logical, intent(in) :: tOscillatorWindow

    !> the window for transitions not included in nexc and energy window (if used)
    real(dp), intent(in) :: oscillatorWindow

    !> excitation energy of state nstat0
    real(dp), intent(out) :: omega

    !> shift vector for potentials in the ground state
    real(dp), intent(in), optional :: shift(:)

    !> non-SCC hamitonian data
    type(OSlakoCont), intent(in), optional :: skHamCont

    !> overlap data
    type(OSlakoCont), intent(in), optional :: skOverCont

    !> excitation energy gradient with respect to atomic positions
    real(dp), intent(out), optional :: excgrad(:,:)

    !> Differentiator for H0 and S matrices.
    class(NonSccDiff), intent(in), optional :: derivator

    !> ground state square density matrix
    real(dp), intent(in), optional :: rhoSqr(:,:,:)

    !> Occupation numbers for natural orbitals from the excited state density matrix
    real(dp), intent(out), optional :: occNatural(:)

    !> the single particle eigenvectors themselves for the excited state density matrix.
    real(dp), intent(out), optional :: naturalOrbs(:,:,:)

    real(dp) :: Ssq(nexc)
    real(dp), allocatable :: gammaMat(:,:), snglPartTransDip(:,:)
    real(dp), allocatable :: stimc(:,:,:), wij(:)
    real(dp), allocatable :: dqex(:), sposz(:), osz(:), xpy(:), xmy(:), pc(:,:)
    real(dp), allocatable :: t(:,:), rhs(:), woo(:), wvv(:), wov(:)
    real(dp), allocatable :: evec(:,:), eval(:), transitionDipoles(:,:)
    integer, allocatable :: win(:), getij(:,:)


    !> array from pairs of single particles states to compound index - should replace with a more
    !> compact data structure in the cases where there are oscilator windows
    integer, allocatable :: iatrans(:,:)

    character, allocatable :: symmetries(:)

    integer :: nocc, nocc_r, nvir_r, nxoo_r, nxvv_r
    integer :: nxov, nxov_ud(2), nxov_r, nxov_d, nxov_rd
    integer :: norb
    integer :: i, j, iSpin, isym, iLev, nStartLev, nEndLev
    integer :: nSpin
    character :: sym

    real(dp) :: energyThreshold

    integer :: nStat


    !> control variables
    logical :: tZVector, tCoeffs, tTradip


    !> printing data
    logical :: tMulliken


    !> should gradients be calculated
    logical :: tForces

    ! ARPACK library variables
    ndigit = -3
    ! Output unit:
    logfil = fdArnoldi
    msgets = 0
    msaitr = 0
    msapps = 0
    mseigt = 0
    mseupd = 0
    if(tArnoldi) then
      msaupd = 1
      msaup2 = 1
    else
      msaupd = 0
      msaup2 = 0
    endif
    ! End of ARPACK communication variables

    @:ASSERT(fdExc > 0)

    ! work out which data files are required, based on whether they have valid file IDs (>0)
    tMulliken = (fdMulliken > 0)
    tCoeffs = (fdCoeffs > 0)
    tTradip = (fdTradip > 0)

    if (tMulliken) then
      open(fdMulliken, file=excitedQOut,position="rewind", status="replace")
      close(fdMulliken)
      open(fdMulliken, file=excitedDipoleOut, position="rewind", status="replace")
      close(fdMulliken)
    end if

    @:ASSERT(fdArnoldi > 0)
    if (tArnoldi) then
      open(fdArnoldi, file=arpackOut, position="rewind", status="replace")
    end if

    nSpin = size(grndEigVal, dim=2)
    @:ASSERT(nSpin > 0 .and. nSpin <=2)

    norb = orb%nOrb

    @:ASSERT(present(excgrad) .eqv. present(shift))
    @:ASSERT(present(excgrad) .eqv. present(skHamCont))
    @:ASSERT(present(excgrad) .eqv. present(skOverCont))
    @:ASSERT(present(excgrad) .eqv. present(rhoSqr))
    @:ASSERT(present(excgrad) .eqv. present(derivator))

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

    if (nexc + 1 >= nxov) then
      write(tmpStr,"(' Insufficent single particle excitations, ',I0,&
          & ', for required number of excited states ',I0)")nxov, nexc
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
    tZVector = tForces .or. tMulliken .or. tCoeffs .or. present(naturalOrbs)

    ! Sanity checks
    nstat = nstat0
    if (nstat < 0 .and. symc /= "S") then
      call error("Linresp: Brightest mode only available for singlets.")
    end if
    if (nstat /= 0 .and. symc == "B") then
      call error("Linresp: Both symmetries not allowed if a specific state is excited")
    end if
    if (tZVector .and. nexc > nxov - 1) then
      call error("Linresp: With gradients/properties, nexc can be greater than the number of&
          & occupied-virtual excitations")
    end if

    ! Select symmetries to process
    if (.not. tSpin) then
      select case (symc)
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
    ALLOCATE(gammaMat(natom, natom))
    ALLOCATE(snglPartTransDip(nxov, 3))
    ALLOCATE(stimc(norb, norb, nSpin))
    ALLOCATE(wij(nxov))
    ALLOCATE(win(nxov))
    ALLOCATE(eval(nexc))
    ALLOCATE(getij(nxov, 2))
    ALLOCATE(transitionDipoles(nxov, 3))
    ALLOCATE(sposz(nxov))

    ! Overlap times wave function coefficients - most routines in DFTB+ use lower triangle (would
    ! remove the need to symmetrize the overlap and ground state density matrix in the main code if
    ! this could be used everywhere in these routines)
    do iSpin = 1, nSpin
      call symm(stimc(:,:,iSpin), "L", SSqr, grndEigVecs(:,:,iSpin))
    end do

    ! ground state Hubbard U softened coulombic interactions
    call sccCalc%getAtomicGammaMatrix(gammaMat, iNeighbour, img2CentCell)

    ! Oscillator strengths for exited states, when needed.
    ALLOCATE(osz(nexc))

    ! Find all single particle transitions and KS energy differences for cases that go from filled
    ! to empty states
    call getSPExcitations(grndEigVal, filling, wij, getij)

    ! put them in ascending energy order
    if (tOscillatorWindow) then
      ! use a stable sort so that degenerate transitions from the same single particle state are
      ! grouped together in the results, allowing these to be selected together (since how intensity
      ! is shared out over degenerate transitions is arbitrary between eigensolvers/platforms).
      call merge_sort(win,wij, 1.0_dp*epsilon(1.0))
    else
      ! do not require stability, use the usual routine to sort, saving an O(N) workspace
      call index_heap_sort(win,wij)
    end if
    wij = wij(win)

    ! dipole strength of transitions between K-S states
    call calcTransitionDipoles(coord0, win, nxov_ud(1), getij, iAtomStart, stimc, grndEigVecs,&
        & snglPartTransDip)

    ! single particle excitation oscillator strengths
    sposz(:) = twothird * wij(:) * sum(snglPartTransDip**2, dim=2)

    if (tOscillatorWindow .and. tZVector ) then
      call error("Incompabilitity between excited state property evaluation and an oscillator&
          & strength window at the moment.")
    end if

    if (tOscillatorWindow .or. tEnergyWindow) then

      if (.not. tEnergyWindow) then

        ! find transitions that are strongly dipole allowed (> oscillatorWindow)
        call dipselect(wij, sposz, win, snglPartTransDip,nxov_rd, oscillatorWindow, grndEigVal,&
            & getij)

      else

        ! energy window above the lowest nexc single particle transitions
        energyThreshold = wij(nexc) + energyWindow
        nxov_r = count(wij <= energyThreshold)

        nxov_d = 0
        if (tOscillatorWindow) then

          ! find transitions that are strongly dipole allowed (> oscillatorWindow)
          if (nxov_r < nxov) then
            ! find transitions that are strongly dipole allowed (> oscillatorWindow)
            call dipselect(wij(nxov_r+1:), sposz(nxov_r+1:), win(nxov_r+1:),&
                & snglPartTransDip(nxov_r+1:,:),nxov_d, oscillatorWindow,&
                & grndEigVal, getij)
          end if

        end if

        nxov_rd = nxov_r + nxov_d

      end if
    else

      nxov_rd = nxov

    end if

    ! just in case energy/dipole windows add no extra states, and is due to an arpack solver
    ! requirement combined with the need to get at least nexc states
    nxov_rd = max(nxov_rd,min(nexc+1,nxov))

    if (fdXplusY >  0) then
      open(fdXplusY, file=XplusYOut, position="rewind", status="replace")
    end if

    if(fdTrans>0) then
      open(fdTrans, file=transitionsOut, position="rewind", status="replace")
      write(fdTrans,*)
    endif

    ! single particle transition dipole file
    if (fdTradip > 0) then
      open(fdTradip, file=transDipOut, position="rewind", status="replace")
      write(fdTradip,*)
      write(fdTradip,'(5x,a,5x,a,2x,a)') "#", 'w [eV]', 'Transition dipole (x,y,z) [Debye]'
      write(fdTradip,*)
      write(fdTradip,'(1x,57("="))')
      write(fdTradip,*)
    endif

    ! excitation energies
    open(fdExc, file=excitationsOut, position="rewind", status="replace")
    write(fdExc,*)
    if (tSpin) then
      write(fdExc,'(5x,a,7x,a,9x,a,9x,a,6x,a,4x,a)') 'w [eV]', 'Osc.Str.', 'Transition','Weight',&
          & 'KS [eV]','D<S*S>'
    else
      write(fdExc,'(5x,a,7x,a,9x,a,9x,a,6x,a,4x,a)') 'w [eV]','Osc.Str.', 'Transition','Weight',&
          & 'KS [eV]','Sym.'
    end if

    write(fdExc,*)
    write(fdExc,'(1x,80("="))')
    write(fdExc,*)

    ! single particle excitations (output file and tagged file if needed).  Was used for nxov_rd =
    ! size(wij), but now for just states that are actually included in the excitation calculation.
    call writeSPExcitations(wij, win, nxov_ud(1), getij, fdSPTrans, sposz, nxov_rd, tSpin)
    ALLOCATE(evec(nxov_rd, nexc))

    do isym = 1, size(symmetries)

      sym = symmetries(isym)
      call buildAndDiagExcMatrix(tSpin, wij(:nxov_rd), sym, win, nxov_ud(1), nxov_rd, iAtomStart,&
          & stimc, grndEigVecs, filling, getij, gammaMat, species0, spinW, fdArnoldiDiagnosis,&
          & eval, evec )

      ! Excitation oscillator strengths for resulting states
      call getOscillatorStrengths(sym, snglPartTransDip(1:nxov_rd,:), wij(:nxov_rd), eval, evec,&
          & filling, win, nxov_ud(1), getij, nstat, osz, tTradip, transitionDipoles)

      if (tSpin) then
        call getExcSpin(Ssq, nxov_ud(1), getij, win, eval, evec, wij(:nxov_rd), filling, stimc,&
            & grndEigVecs)
        call writeExcitations(sym, osz, nexc, nxov_ud(1), getij, win, eval, evec, wij(:nxov_rd),&
            & fdXplusY, fdTrans, fdTradip, transitionDipoles, tWriteTagged, fdTagged, fdExc, Ssq)
      else
        call writeExcitations(sym, osz, nexc, nxov_ud(1), getij, win, eval, evec, wij(:nxov_rd),&
            & fdXplusY, fdTrans, fdTradip, transitionDipoles, tWriteTagged, fdTagged, fdExc)
      end if

    end do

    if (tArnoldi) then
      close(fdArnoldi)
    end if

    if (fdTrans > 0) close(fdTrans)
    if (fdXplusY > 0) close(fdXplusY)
    if (fdExc > 0) close(fdExc)
    if (fdTradip > 0) close(fdTradip)

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
        nEndLev = nexc

        if (tForces) then
          call error("Forces currently not available unless a single excited state is specified")
        end if

      else
        nStartLev = nstat
        nEndLev = nstat
      end if

      if (tSpin) then
        call error("Z vector evaluation does not currently support spin polarization.")
      end if

      if (any( abs(filling) > elecTolMax .and. abs(filling-2.0_dp) > elecTolMax ) ) then
        call error("Fractional fillings not currently possible for excited state property&
            & calculations")
      end if

      ! redefine if needed (generalize it for spin-polarized and fractional occupancy)
      nocc = int(rnel) / 2

      ! count virtual and occupied states
      call getNorb_r(nxov_rd, win, getij, nocc, nocc_r, nvir_r)

      ! size of occ-occ and virt-virt blocks
      nxoo_r = (nocc_r * (nocc_r + 1)) / 2
      nxvv_r = (nvir_r * (nvir_r + 1)) / 2

      ! Arrays needed for Z vector
      ALLOCATE(xpy(nxov_rd))
      ALLOCATE(xmy(nxov_rd))
      ALLOCATE(t(norb, norb))
      ALLOCATE(rhs(nxov_rd))
      ALLOCATE(woo(nxoo_r))
      ALLOCATE(wvv(nxvv_r))
      ALLOCATE(wov(nxov_rd))
      ALLOCATE(iatrans(1:nocc, nocc+1:norb))

      ! Arrays for gradients and Mulliken analysis
      if (tZVector) then
        ALLOCATE(dqex(natom))
        ALLOCATE(pc(norb, norb))
      end if

      ! Furche terms: X+Y, X-Y
      xpy(:nxov_rd) = sqrt(wij(:nxov_rd)) / sqrt(omega) * evec(:nxov_rd,nstat)
      xmy(:nxov_rd) = sqrt(omega) / sqrt(wij(:nxov_rd)) * evec(:nxov_rd,nstat)

      ! set up transition indexing
      call rindxov_array(win, nocc, nxov, getij, iatrans)

      do iLev = nStartLev, nEndLev
        omega = sqrt(eval(iLev))
        ! Furche terms: X+Y, X-Y
        xpy(:nxov_rd) = evec(:nxov_rd,iLev) * sqrt(wij(:nxov_rd) / omega)
        xmy(:nxov_rd) = evec(:nxov_rd,iLev) * sqrt(omega / wij(:nxov_rd))

        ! solve for Z and W to get excited state density matrix
        call getZVectorEqRHS(xpy, xmy, win, iAtomStart, nocc, nocc_r,&
            & nxov_ud(1), getij, iatrans, natom, species0,grndEigVal(:,1),&
            & stimc, grndEigVecs, gammaMat, spinW, omega, sym, rhs, t,&
            & wov, woo, wvv)
        call solveZVectorEq(rhs, win, nxov_ud(1), getij, natom, iAtomStart,&
            & stimc, gammaMat, wij(:nxov_rd), grndEigVecs)
        call calcWVectorZ(rhs, win, nocc, nocc_r, nxov_ud(1), getij, iAtomStart,&
            & stimc, grndEigVecs, gammaMat, grndEigVal(:,1), wov, woo, wvv)
        call calcPMatrix(t, rhs, win, getij, pc)

        call writeCoeffs(pc, grndEigVecs, filling, nocc, fdCoeffs,&
            & tCoeffs, tGrndState, occNatural, naturalOrbs)

        ! Make MO to AO transformation of the excited density matrix
        call makeSimiliarityTrans(pc, grndEigVecs(:,:,1))

        call getExcMulliken(iAtomStart, pc, SSqr, dqex)
        if (tMulliken) then
          call writeExcMulliken(sym, iLev, dq, dqex, coord0, fdMulliken)
        end if

        if (tForces) then
          call addGradients(sym, nxov_rd, natom, species0, iAtomStart, norb,&
              & nocc, nocc_r, nxov_ud(1), getij, win, grndEigVecs, pc, stimc,&
              & dq, dqex, gammaMat, HubbardU, spinW, shift, woo, wov, wvv,&
              & xpy, coord0, orb, skHamCont, skOverCont, derivator,&
              & rhoSqr(:,:,1), excgrad)
        end if

      end do

      if (nstat == 0) then
        omega = 0.0_dp
      end if

    end if

  end subroutine LinRespGrad_old


  !> Builds and diagonalizes the excitation matrix via iterative technique.
  subroutine buildAndDiagExcMatrix(tSpin, wij, sym, win, nmatup, nxov, iAtomStart, stimc,&
      & grndEigVecs, filling, getij, gammaMat, species0, spinW, fdArnoldiDiagnosis, eval, evec)

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

    !> resulting eigenvalues for transitions
    real(dp), intent(out) :: eval(:)

    !> eigenvectors for transitions
    real(dp), intent(out) :: evec(:,:)

    real(dp), allocatable :: workl(:), workd(:), resid(:), vv(:,:), qij(:)
    real(dp) :: sigma
    integer :: iparam(11), ipntr(11)
    integer :: ido, ncv, lworkl, info
    logical, allocatable :: selection(:)
    logical :: rvec
    integer :: nexc, natom

    integer :: iState
    real(dp), allocatable :: Hv(:), orthnorm(:,:)

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

    resid = 0.0_dp
    workd = 0.0_dp

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
            & ido,info
        call error(tmpStr)
      end if

      ! Action of excitation supermatrix on supervector
      call omegatvec(tSpin, workd(ipntr(1):ipntr(1)+nxov-1), workd(ipntr(2):ipntr(2)+nxov-1),&
          & wij, sym, win, nmatup, iAtomStart, stimc, grndEigVecs, filling, getij, gammaMat,&
          & species0, spinW)

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
            & grndEigVecs, filling, getij, gammaMat, species0, spinW)
        write(fdArnoldiDiagnosis,"(I4,4E16.8)")iState, dot_product(Hv,evec(:,iState))-eval(iState),&
            & sqrt(sum( (Hv-evec(:,iState)*eval(iState) )**2 )), orthnorm(iState,iState) - 1.0_dp,&
            & max(maxval(orthnorm(:iState-1,iState)), maxval(orthnorm(iState+1:,iState)))
      end do
      close(fdArnoldiDiagnosis)
    end if

  end subroutine buildAndDiagExcMatrix


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

    transitionDipoles = 0.0_dp
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
      call transitionDipole(snglPartTransDip, wnij, eval, evec,&
          & transitionDipoles)
    end if

  end subroutine getOscillatorStrengths


  !> Calculate <S^2> as a measure of spin contamination (smaller magnitudes are better, 0.5 is
  !> considered an upper threshold for reliability according to Garcia thesis)
  subroutine getExcSpin(Ssq, nmatup, getij, win, eval, evec, wij, filling, stimc, grndEigVecs)

    !> spin contamination
    real(dp), intent(out) :: Ssq(:)

    !> number of spin up excitations
    integer, intent(in) :: nmatup

    !> index for composite excitations to specific occupied and empty states
    integer, intent(in) :: getij(:,:)

    !> single particle excitations
    integer, intent(in) :: win(:)

    !> Casida exitation energies
    real(dp), intent(in) :: eval(:)

    !> Casida excited eigenvectors
    real(dp), intent(in) :: evec(:,:)

    !> single particle excitation energies
    real(dp), intent(in) :: wij(:)

    !> occupations in ground state
    real(dp), intent(in) :: filling(:,:)

    !> Overlap times ground state eigenvectors
    real(dp), intent(in) :: stimc(:,:,:)

    !> Ground state eigenvectors
    real(dp), intent(in) :: grndEigVecs(:,:,:)

    integer:: i, k, l, m, ia, jb, ii, aa, jj, bb
    integer:: nmat, nexc, nup, ndwn
    real(dp) :: rsqw, TDvnorm
    real(dp), allocatable :: TDvec(:), TDvec_sq(:)
    integer, allocatable :: TDvin(:)
    logical :: ud_ia, ud_jb
    real(dp) :: s_iaja, s_iaib, s_iajb, tmp
    real(dp) :: wnij(size(wij))

    nmat = size(evec, dim=1)
    nexc = size(Ssq)
    nup = ceiling(sum(filling(:,1)))
    ndwn = ceiling(sum(filling(:,2)))
    ALLOCATE(TDvec(nmat))
    ALLOCATE(TDvec_sq(nmat))
    ALLOCATE(TDvin(nmat))

    call wtdn(wij, filling, win, nmatup, nmat, getij, wnij)

    do i = 1, nexc
      rsqw = 1.0_dp / sqrt(sqrt(eval(i)))
      TDvec(:) = sqrt(wnij(:)) * rsqw * evec(:,i)
      TDvnorm = 1.0_dp / sqrt(sum(TDvec**2))
      TDvec(:) = TDvec(:) * TDvnorm
      TDvec_sq = TDvec**2

      ! put these transition dipoles in order of descending magnitude
      call index_heap_sort(TDvin, TDvec_sq)
      TDvin = TDvin(nmat:1:-1)
      TDvec_sq = TDvec_sq(TDvin)

      ! S_{ia,ja}
      s_iaja = 0.0_dp
      do k = 1, nmat
        ia = TDvin(k)
        call indxov(win, ia, getij, ii, aa)
        ud_ia = (win(ia) <= nmatup)
        do l = 1, nmat
          jb = TDvin(l)
          call indxov(win, jb, getij, jj, bb)
          ud_jb = (win(jb) <= nmatup)

          if ( (bb /= aa) .or. (ud_jb .neqv. ud_ia) ) then
            cycle
          end if

          tmp = 0.0_dp
          if (ud_ia) then
            do m = 1,ndwn
              tmp = tmp + MOoverlap(ii,m,stimc,grndEigVecs) * MOoverlap(jj,m,stimc,grndEigVecs)
            end do
          else
            do m = 1,nup
              tmp = tmp + MOoverlap(m,ii,stimc,grndEigVecs) * MOoverlap(m,jj,stimc,grndEigVecs)
            end do
          end if

          s_iaja = s_iaja + TDvec(ia)*TDvec(jb)*tmp

        end do
      end do

      ! S_{ia,ib}
      s_iaib = 0.0_dp
      do k = 1, nmat
        ia = TDvin(k)
        call indxov(win, ia, getij, ii, aa)
        ud_ia = (win(ia) <= nmatup)
        do l = 1, nmat
          jb = TDvin(l)
          call indxov(win, jb, getij, jj, bb)
          ud_jb = (win(jb) <= nmatup)

          if ( (ii /= jj) .or. (ud_jb .neqv. ud_ia) ) then
            cycle
          end if

          tmp = 0.0_dp
          if (ud_ia) then
            do m = 1,ndwn
              tmp = tmp + MOoverlap(aa,m,stimc,grndEigVecs) * MOoverlap(bb,m,stimc,grndEigVecs)
            end do
          else
            do m = 1,nup
              tmp = tmp + MOoverlap(m,aa,stimc,grndEigVecs) * MOoverlap(m,bb,stimc,grndEigVecs)
            end do
          end if

          s_iaib = s_iaib + TDvec(ia)*TDvec(jb)*tmp
        end do
      end do

      ! S_{ia,jb}
      s_iajb = 0.0_dp
      do k = 1, nmat
        ia = TDvin(k)
        call indxov(win, ia, getij, ii, aa)
        ud_ia = (win(ia) <= nmatup)
        if (.not. ud_ia ) then
          cycle
        end if
        do l = 1, nmat
          jb = TDvin(l)
          call indxov(win, jb, getij, jj, bb)
          ud_jb = (win(jb) <= nmatup)

          if ( ud_jb ) cycle

          s_iajb = s_iajb + TDvec(ia)*TDvec(jb) * MOoverlap(aa,bb,stimc,grndEigVecs)&
              & * MOoverlap(ii,jj,stimc,grndEigVecs)

        end do
      end do

      Ssq(i) =  s_iaja - s_iaib - 2.0_dp*s_iajb

    end do

  end subroutine getExcSpin


  !> Build right hand side of the equation for the Z-vector and those parts of the W-vectors which
  !> do not depend on Z.
  subroutine getZVectorEqRHS(xpy, xmy, win, iAtomStart, homo, nocc, nmatup, getij, iatrans, natom,&
      & species0, grndEigVal, stimc, c, gammaMat, spinW, omega, sym, rhs, t, wov, woo, wvv)

    !> X+Y Furche term
    real(dp), intent(in) :: xpy(:)

    !> X-Y Furche term
    real(dp), intent(in) :: xmy(:)

    !> index array for single particle transitions
    integer, intent(in) :: win(:)

    !> index vector for S and H matrices
    integer, intent(in) :: iAtomStart(:)

    !> highest occupied level
    integer, intent(in) :: homo

    !> number of filled states
    integer, intent(in) :: nocc

    !> number of same spin excitations
    integer, intent(in) :: nmatup

    !> index array between transitions in square and 1D representations
    integer, intent(in) :: getij(:,:)

    !> index array from orbital pairs to compound index
    integer, intent(in) :: iatrans(:,homo+1:)

    !> number of central cell atoms
    integer, intent(in) :: natom

    !> central cell chemical species
    integer, intent(in) :: species0(:)

    !> ground state wavefunctions
    real(dp), intent(in) :: grndEigVal(:)

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
    real(dp), intent(out) :: t(:,:)

    !> W vector occupied-virtual part
    real(dp), intent(out) :: wov(:)

    !> W vector occupied part
    real(dp), intent(out) :: woo(:)

    !> W vector virtual part
    real(dp), intent(out) :: wvv(:)

    real(dp), allocatable :: xpyq(:), qij(:), gamxpyq(:), qgamxpyq(:), gamqt(:)
    integer :: nxov, nxoo, nxvv
    integer :: i, j, a, b, ia, ib, ij, ab, ja
    real(dp) :: tmp1, tmp2
    logical :: updwn

    nxov = size(rhs)
    nxoo = size(woo)
    nxvv = size(wvv)

    ALLOCATE(xpyq(natom))
    ALLOCATE(qij(natom))
    ALLOCATE(gamxpyq(natom))
    ALLOCATE(gamqt(natom))
    ALLOCATE(qgamxpyq(max(nxoo, nxvv)))

    t(:,:) = 0.0_dp
    rhs(:) = 0.0_dp
    wov(:) = 0.0_dp
    woo(:) = 0.0_dp
    wvv(:) = 0.0_dp
    xpyq(:) = 0.0_dp

    ! Build t_ab = 0.5 * sum_i (X+Y)_ia (X+Y)_ib + (X-Y)_ia (X-Y)_ib
    ! and w_ab = Q_ab with Q_ab as in (B16) but with corrected sign.
    ! factor 1 / (1 + delta_ab) follows later
    do ia = 1, nxov
      call indxov(win, ia, getij, i, a)

      ! BA: is T_aa = 0?
      do b = homo + 1, a
        ib = iatrans(i, b)
        call rindxvv(homo, a, b, ab)
        tmp1 = xpy(ia) * xpy(ib) + xmy(ia) * xmy(ib)
        tmp2 = omega * (xpy(ia) * xmy(ib)+ xmy(ia) * xpy(ib))
        t(a,b) = t(a,b) + 0.5_dp * tmp1
        ! to prevent double counting
        if (a /= b) then
          t(b,a) = t(b,a) + 0.5_dp * tmp1
        end if
        ! Note: diagonal elements will be multiplied by 0.5 later.
        wvv(ab) = wvv(ab) + grndEigVal(i) * tmp1 + tmp2
      end do

      ! Build t_ij = 0.5 * sum_a (X+Y)_ia (X+Y)_ja + (X-Y)_ia (X-Y)_ja and 1 / (1 + delta_ij) Q_ij
      ! with Q_ij as in eq. (B9) (1st part of w_ij)
      do j = i, homo
        ja = iatrans(j,a)
        call rindxvv(homo-nocc, j, i, ij)
        tmp1 = xpy(ia) * xpy(ja) + xmy(ia) * xmy(ja)
        tmp2 = omega * (xpy(ia) * xmy(ja) + xmy(ia) * xpy(ja))
        ! Note, there is a typo in Heringer et al. J. Comp Chem 28, 2589.
        ! The sign must be negative see Furche, J. Chem. Phys, 117 7433 (2002).
        t(i,j) = t(i,j) - 0.5_dp * tmp1
        ! to prevent double counting
        if (i /= j) then
          t(j,i) = t(j,i) - 0.5_dp * tmp1
        end if
        woo(ij) = woo(ij) - grndEigVal(a) * tmp1 + tmp2
      end do
    end do

    ! Build xpyq = sum_ia (X+Y)_ia
    do ia = 1, nxov
      call indxov(win, ia, getij, i, a)
      updwn = (win(ia) <= nmatup)
      call transq(i, a, iAtomStart, updwn, stimc, c, qij)
      xpyq(:) = xpyq + xpy(ia) * qij
    end do

    ! qgamxpyq(ab) = sum_jc K_ab,jc (X+Y)_jc
    if (sym == "S") then
      call hemv(gamxpyq, gammaMat,  xpyq)
      do ab = 1, nxvv
        call indxvv(homo, ab, a, b)
        call transq(a, b, iAtomStart, updwn, stimc, c, qij)
        qgamxpyq(ab) = 2.0_dp * sum(qij * gamxpyq)
      end do
    else
      do ab = 1, nxvv
        call indxvv(homo, ab, a, b)
        call transq(a, b, iAtomStart, updwn, stimc, c, qij)
        qgamxpyq(ab) = 2.0_dp * sum(qij * xpyq * spinW(species0))
      end do
    end if

    ! rhs(ia) -= Qia = sum_b (X+Y)_ib * qgamxpyq(ab))
    do ia = 1, nxov
      call indxov(win, ia, getij, i, a)
      do b = homo + 1, a
        call rindxvv(homo, a, b, ab)
        ib = iatrans(i,b)
        rhs(ia) = rhs(ia) - 2.0_dp * xpy(ib) * qgamxpyq(ab)
        ! Since qgamxpyq has only upper triangle
        if (a /= b) then
          rhs(ib) = rhs(ib) - 2.0_dp * xpy(ia) * qgamxpyq(ab)
        end if
      end do
    end do

    ! -rhs = -rhs - sum_j (X + Y)_ja H + _ij[X + Y]
    if (sym == "S") then
      do ij = 1, nxoo
        qgamxpyq(ij) = 0.0_dp
        call indxoo(homo, nocc, ij, i, j)
        call transq(i, j, iAtomStart, updwn, stimc, c, qij)
        ! qgamxpyq(ij) = sum_kb K_ij,kb (X+Y)_kb
        qgamxpyq(ij) = 2.0_dp * sum(qij * gamxpyq)
      end do
    else
      do ij = 1, nxoo
        qgamxpyq(ij) = 0.0_dp
        call indxoo(homo, nocc, ij, i, j)
        call transq(i, j, iAtomStart, updwn, stimc, c, qij)
        qgamxpyq(ij) = 2.0_dp * sum(qij * xpyq * spinW(species0))
      end do
    end if

    ! rhs(ia) += Qai = sum_j (X+Y)_ja qgamxpyq(ij)
    ! add Qai to Wia as well.
    do ia = 1, nxov
      call indxov(win, ia, getij, i, a)
      do j = i, homo
        ja = iatrans(j, a)
        ij = i-homo+nocc + ((j-homo+nocc - 1) * (j-homo+nocc)) / 2
        tmp1 = 2.0_dp * xpy(ja) * qgamxpyq(ij)
        rhs(ia) = rhs(ia) + tmp1
        wov(ia) = wov(ia) + tmp1
        if (i /= j) then
          tmp2 = 2.0_dp * xpy(ia) * qgamxpyq(ij)
          rhs(ja) = rhs(ja) + tmp2
          wov(ja) = wov(ja) + tmp2
        end if
      end do
    end do

    ! gamxpyq(iAt2) = sum_ij q_ij(iAt2) T_ij
    gamxpyq(:) = 0.0_dp
    do ij = 1, nxoo
      call indxoo(homo, nocc, ij, i, j)
      call transq(i, j, iAtomStart, updwn, stimc, c, qij)
      if (i == j) then
        gamxpyq(:) = gamxpyq(:) + t(i,j) * qij(:)
      else
        ! factor 2 because of symmetry of the matrix
        gamxpyq(:) = gamxpyq(:) + 2.0_dp  * t(i,j) * qij(:)
      end if
    end do

    ! gamxpyq(iAt2) += sum_ab q_ab(iAt2) T_ab
    do ab = 1, nxvv
      call indxvv(homo, ab, a, b)
      call transq(a, b, iAtomStart, updwn, stimc, c, qij)
      if (a == b) then
        gamxpyq(:) = gamxpyq(:) + t(a,b) * qij(:)
      else
        ! factor 2 because of symmetry of the matrix
        gamxpyq(:) = gamxpyq(:) + 2.0_dp * t(a,b) * qij(:)
      end if
    end do

    ! gamqt(iAt1) = sum_iAt2 gamma_iAt1,iAt2 gamxpyq(iAt2)
    call hemv(gamqt, gammaMat, gamxpyq)

    ! rhs -= sum_q^ia(iAt1) gamxpyq(iAt1)
    do ia = 1, nxov
      call indxov(win, ia, getij, i, a)
      updwn = (win(ia) <= nmatup)
      call transq(i, a, iAtomStart, updwn, stimc, c, qij)
      rhs(ia) = rhs(ia) - 4.0_dp * sum(qij * gamqt)
    end do

    ! Furche vectors
    do ij = 1, nxoo
      call indxoo(homo, nocc, ij, i, j)
      call transq(i, j, iAtomStart, updwn, stimc, c, qij)
      woo(ij) = woo(ij) + 4.0_dp * sum(qij * gamqt)
    end do

  end subroutine getZVectorEqRHS


  !> Solving the (A+B) Z = -R equation via conjugate gradient
  subroutine solveZVectorEq(rhs, win, nmatup, getij, natom, iAtomStart, stimc, gammaMat, wij, c)

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

    integer :: nxov
    integer :: ia, i, a, k
    real(dp) :: rhs2(size(rhs)),rkm1(size(rhs)),pkm1(size(rhs)),apk(size(rhs))
    real(dp) :: qTmp(nAtom), rs, alphakm1, tmp1, tmp2, bkm1
    real(dp), allocatable :: qij(:)
    logical :: updwn

    nxov = size(rhs)
    allocate(qij(nAtom))

    ! Choosing a start value
    ! rhs2 = rhs / (A+B)_ia,ia (diagonal of the supermatrix sum A+B)
    do ia = 1, nxov
      call indxov(win, ia, getij, i, a)
      updwn = (win(ia) <= nmatup)
      call transq(i, a, iAtomStart, updwn, stimc, c, qij)
      call hemv(qTmp, gammaMat, qij)
      rs = 4.0_dp * dot_product(qij, qTmp) + wij(ia)
      rhs2(ia) = rhs(ia) / rs
    end do

    ! unit vector
    rhs2 = 1.0_dp / sqrt(real(nxov,dp))

    ! Free some space, before entering the apbw routine
    deallocate(qij)

    ! action of matrix on vector
    call apbw(rkm1, rhs2, wij, nxov, natom, win, nmatup, getij, iAtomStart,&
        & stimc, c, gammaMat)

    rkm1 = rhs - rkm1
    pkm1 = rkm1

    ! Iteration: should be convergent in at most nxov steps for a quadradic surface, so set higher
    do k = 1, nxov**2

      ! action of matrix on vector
      call apbw(apk, pkm1, wij, nxov, natom,&
          & win, nmatup, getij, iAtomStart, stimc, c, gammaMat)

      tmp1 = dot_product(rkm1, rkm1)
      tmp2 = dot_product(pkm1, apk)
      alphakm1 = tmp1 / tmp2

      rhs2 = rhs2 + alphakm1 * pkm1

      rkm1 = rkm1 -alphakm1 * apk

      tmp2 = dot_product(rkm1, rkm1)

      ! residual
      if (tmp2 <= epsilon(1.0_dp)**2) then
        exit
      end if

      if (k == nxov**2) then
        call error("solveZVectorEq : Z vector not converged!")
      end if

      bkm1 = tmp2 / tmp1

      pkm1 = rkm1 + bkm1 * pkm1

    end do

    rhs(:) = rhs2(:)

  end subroutine solveZVectorEq


  !> Calculate Z-dependent parts of the W-vectors and divide diagonal elements of W_ij and W_ab by
  !> 2.
  subroutine calcWvectorZ(zz, win, homo, nocc, nmatup, getij, iAtomStart, stimc, c, gammaMat,&
      & grndEigVal, wov, woo, wvv)

    !> Z vector
    real(dp), intent(in) :: zz(:)

    !> index array for transitions
    integer, intent(in) :: win(:)

    !> highest occupied level
    integer, intent(in) :: homo

    !> number of filled levels
    integer, intent(in) :: nocc

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
    real(dp), intent(in) :: grndEigVal(:)

    !> W vector occupied-virtual part
    real(dp), intent(inout) :: wov(:)

    !> W vector occupied part
    real(dp), intent(inout) :: woo(:)

    !> W vector virtual part
    real(dp), intent(inout) :: wvv(:)

    integer :: nxov, nxoo, nxvv, natom
    integer :: ij, ia, ab, i, j, a, b, iAt1
    real(dp), allocatable :: qij(:), gamxpyq(:), zq(:)
    logical :: updwn

    nxov = size(zz)
    natom = size(gammaMat, dim=1)
    nxoo = size(woo)
    nxvv = size(wvv)
    ALLOCATE(qij(natom))
    ALLOCATE(gamxpyq(natom))
    ALLOCATE(zq(natom))

    ! Adding missing epsilon_i * Z_ia term to W_ia
    do ia = 1, nxov
      call indxov(win, ia, getij, i, a)
      wov(ia) = wov(ia) + zz(ia) * grndEigVal(i)
    end do

    ! Missing sum_kb 4 K_ijkb Z_kb term in W_ij: zq(iAt1) = sum_kb q^kb(iAt1) Z_kb
    do iAt1 = 1, natom
      zq(iAt1) = 0.0_dp
      do ia = 1, nxov
        call indxov(win, ia, getij, i, a)
        updwn = (win(ia) <= nmatup)
        call transq(i, a, iAtomStart, updwn, stimc, c, qij)
        zq(iAt1) = zq(iAt1) + zz(ia) * qij(iAt1)
      end do
    end do

    call hemv(gamxpyq, gammaMat, zq)

    ! sum_iAt1 qij(iAt1) gamxpyq(iAt1)
    do ij = 1, nxoo
      call indxoo(homo, nocc, ij, i, j)
      call transq(i, j, iAtomStart, updwn, stimc, c, qij)
      do iAt1 = 1, natom
        ! W contains 1/2 for i == j.
        woo(ij) = woo(ij) + 4.0_dp * qij(iAt1) * gamxpyq(iAt1)
      end do
    end do

    ! Divide diagonal elements of W_ij by 2.
    do ij = 1, nxoo
      call indxoo(homo, nocc, ij, i, j)
      if (i == j) then
        woo(ij) = 0.5_dp * woo(ij)
      end if
    end do

    ! Divide diagonal elements of W_ab by 2.
    do ab = 1, nxvv
      call indxvv(homo, ab, a, b)
      if (a == b) then
        wvv(ab) = 0.5_dp * wvv(ab)
      end if
    end do

  end subroutine calcWvectorZ


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


  !> Excited state Mulliken charges and dipole moments written to disc
  subroutine writeExcMulliken(sym, nstat, dq, dqex, coord0, fdMulliken)

    !> symmetry label
    character, intent(in) :: sym

    !> state index
    integer, intent(in) :: nstat

    !> ground state gross charge
    real(dp), intent(in) :: dq(:)

    !> change in atomic charges from ground to excited state
    real(dp), intent(in) :: dqex(:)

    !> central cell coordinates
    real(dp), intent(in) :: coord0(:,:)

    !> file unit for Mulliken data
    integer, intent(in) :: fdMulliken

    integer :: natom, m
    real(dp) :: dipol(3), dipabs

    natom = size(dq)

    @:ASSERT(size(dq) == size(dqex))
    @:ASSERT(all(shape(coord0) == [3,nAtom]))

    ! Output of excited state Mulliken charges
    open(fdMulliken, file=excitedQOut,position="append")
    write(fdMulliken, "(a,a,i2)") "# MULLIKEN CHARGES of excited state ",&
        & sym, nstat
    write(fdMulliken, "(a,2x,A,i4)") "#", 'Natoms =',natom
    write(fdMulliken, "('#',1X,A4,T15,A)")'Atom','netCharge'
    write(fdMulliken,'("#",41("="))')
    do m = 1,  natom
      write(fdMulliken,"(i5,1x,f16.8)") m, -dq(m) - dqex(m)
    end do
    close(fdMulliken)

    ! Calculation of excited state dipole moment
    dipol(:) = -1.0_dp * matmul(coord0, dq + dqex)
    dipabs = sqrt(sum(dipol**2))

    open(fdMulliken, file=excitedDipoleOut, position="append")
    write(fdMulliken, "(a,a,i2)") "Mulliken analysis of excited state ",&
        & sym, nstat
    write(fdMulliken, '(42("="))')
    write(fdMulliken, "(a)") " "
    write(fdMulliken, "(a)") "Mulliken exc. state dipole moment [Debye]"
    write(fdMulliken, '(42("="))')
    write(fdMulliken, "(3f14.8)") (dipol(m) * au__Debye, m = 1, 3)
    write(fdMulliken, "(a)") " "
    write(fdMulliken, "(a)") "Norm of exc. state dipole moment [Debye]"
    write(fdMulliken, '(42("="))')
    write(fdMulliken, "(e20.12)") dipabs * au__Debye
    write(fdMulliken, *)
    close(fdMulliken)

  end subroutine writeExcMulliken


  !> Calculate transition moments for transitions between Kohn-Sham states, including spin-flipping
  !> transitions
  subroutine calcTransitionDipoles(coord0, win, nmatup, getij, iAtomStart, stimc, grndEigVecs,&
      & snglPartTransDip)

    !> Atomic positions
    real(dp), intent(in) :: coord0(:,:)

    !> transition energies
    integer, intent(in) :: win(:)

    !> number of same-spin transitions
    integer, intent(in) :: nmatup

    !> index array for ground state square matrices
    integer, intent(in) :: iAtomStart(:)

    !> index array for excitation pairs
    integer, intent(in) :: getij(:,:)

    !> overlap times ground state wavefunctions
    real(dp), intent(in) :: stimc(:,:,:)

    !> ground state wavefunctions
    real(dp), intent(in) :: grndEigVecs(:,:,:)

    !> resulting transition dipoles
    real(dp), intent(out) :: snglPartTransDip(:,:)

    integer :: nxov, natom
    integer :: indm, ii, jj
    real(dp), allocatable :: qij(:)
    logical :: updwn

    nxov = size(win)
    natom = size(coord0, dim=2)

    ALLOCATE(qij(natom))

    ! Calculate transition dipole elements
    do indm = 1, nxov
      call indxov(win, indm, getij, ii, jj)
      updwn = (win(indm) <= nmatup)
      call transq(ii, jj, iAtomStart, updwn, stimc, grndEigVecs, qij)
      snglPartTransDip(indm, :) = matmul(coord0, qij)
    end do

  end subroutine calcTransitionDipoles


  !> Calculation of force from derivatives of excitation energy
  !> 1. we need the ground and excited Mulliken charges
  !> 2. we need P,(T,Z),W, X + Y from linear response
  !> 3. calculate dsmndr, dhmndr (dS/dR, dh/dR), dgabda (dGamma_{IAt1,IAt2}/dR_{IAt1}),
  !> dgext (dGamma-EXT_{IAt1,k}/dR_{IAt1})
  subroutine addGradients(sym, nxov, natom, species0, iAtomStart, norb, homo,&
      & nocc, nmatup, getij, win, grndEigVecs, pc, stimc, dq, dqex, gammaMat,&
      & HubbardU, spinW, shift, woo, wov, wvv, xpy, coord0, orb,&
      & skHamCont, skOverCont, derivator, rhoSqr, excgrad)

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
    integer, intent(in) :: homo

    !> number of occupied states in calculation (not neccessarily same as HOMO in the case of
    !> windowing)
    integer, intent(in) :: nocc

    !> single particle excitation energies
    integer, intent(in) :: win(:)

    !> number of up->up transitions
    integer, intent(in) :: nmatup

    !> index array from composite transition index to specific single particle states
    integer, intent(in) :: getij(:,:)

    !> ground state eigenvectors
    real(dp), intent(in) :: grndEigVecs(:,:,:)

    !> transition density matrix
    real(dp), intent(in) :: pc(:,:)

    !> overlap times ground state eigenvectors
    real(dp), intent(in) :: stimc(:,:,:)

    !> ground state gross charges
    real(dp), intent(in) :: dq(:)

    !> charge differences from ground to excited state
    real(dp), intent(in) :: dqex(:)

    !> softened coulomb matrix
    real(dp), intent(in) :: gammaMat(:,:)

    !> ground state Hubbard U values
    real(dp), intent(in) :: HubbardU(:)

    !> ground state spin derivatives for each species
    real(dp), intent(in) :: spinW(:)

    !> ground state potentials (shift vector)
    real(dp), intent(in) :: shift(:)

    !> W vector occupied part
    real(dp), intent(in) :: woo(:)

    !> W vector occupied-virtual part
    real(dp), intent(in) :: wov(:)

    !> W vector virtual part
    real(dp), intent(in) :: wvv(:)

    !> X+Y Furche term
    real(dp), intent(in) :: xpy(:)

    !> central cell atomic coordinates
    real(dp), intent(in) :: coord0(:,:)

    !> data type for atomic orbital information
    type(TOrbitals), intent(in) :: orb

    !> H0 data
    type(OSlakoCont), intent(in) :: skHamCont

    !> overlap data
    type(OSlakoCont), intent(in) :: skOverCont

    !> Differentiatior for the non-scc matrices
    class(NonSccDiff), intent(in) :: derivator

    !> ground state density matrix for spin-free case
    real(dp), intent(in) :: rhoSqr(:,:)

    !> resulting excited state gradient
    real(dp), intent(out) :: excgrad(:,:)

    real(dp), allocatable :: shift_excited(:), xpyq(:)
    real(dp), allocatable :: shxpyq(:), xpycc(:,:), wcc(:,:)
    real(dp), allocatable :: qij(:), temp(:)
    real(dp), allocatable :: dH0(:,:,:), dS(:,:,:)
    integer :: ia, i, j, a, b, ab, ij, m, n, mu, nu, xyz, iAt1, iAt2
    integer :: indalpha, indalpha1, indbeta, indbeta1
    integer :: iSp1, iSp2
    real(dp) :: tmp1, tmp2, tmp3, tmp4, tmp5, tmp6, tmp7, rab
    real(dp) :: diffvec(3), dgab(3), tmp3a, tmp3b

    logical :: updwn

    integer :: nxoo, nxvv

    ALLOCATE(shift_excited(natom))
    ALLOCATE(xpyq(natom))
    ALLOCATE(shxpyq(natom))
    ALLOCATE(xpycc(norb, norb))
    ALLOCATE(wcc(norb, norb))
    ALLOCATE(qij(natom))
    ALLOCATE(temp(norb))

    ALLOCATE(dH0(orb%mOrb, orb%mOrb, 3))
    ALLOCATE(dS(orb%mOrb, orb%mOrb, 3))

    nxoo = size(woo)
    nxvv = size(wvv)

    excgrad = 0.0_dp

    ! excited state potentials at atomic sites
    call hemv(shift_excited, gammaMat, dqex)

    ! xypq(alpha) = sum_ia (X+Y)_ia q^ia(alpha)
    ! complexity norb * norb * norb
    xpyq(:) = 0.0_dp
    do ia = 1, nxov
      call indxov(win, ia, getij, i, a)
      updwn = (win(ia) <= nmatup)
      call transq(i, a, iAtomStart, updwn, stimc, grndEigVecs, qij)
      xpyq(:) = xpyq(:) + xpy(ia) * qij(:)
    end do

    ! complexity norb * norb
    shxpyq(:) = 0.0_dp
    if (sym == "S") then
      call hemv(shxpyq, gammaMat, xpyq)
    else
      shxpyq(:) = xpyq(:) * spinW(species0)
    end if

    ! calculate xpycc
    ! (xpycc)_{mu nu} = sum_{ia} (X + Y)_{ia} (grndEigVecs(mu,i)grndEigVecs(nu,a)
    ! + grndEigVecs(nu,i)grndEigVecs(mu,a))
    ! complexity norb * norb * norb
    xpycc(:,:) = 0.0_dp

    ! xpycc(mu,nu) = sum_ia (X+Y)_ia grndEigVecs(mu,i) grndEigVecs(nu,a)
    ! xpycc(mu, nu) += sum_ia (X+Y)_ia grndEigVecs(mu,a) grndEigVecs(nu,i)
    do ia = 1, nxov
      call indxov(win, ia, getij, i, a)
      do nu = 1, norb
        do mu = 1, norb
          xpycc(mu,nu) = xpycc(mu,nu) + xpy(ia) *&
              & ( grndEigVecs(mu,i,1)*grndEigVecs(nu,a,1)&
              & + grndEigVecs(mu,a,1)*grndEigVecs(nu,i,1) )
        end do
      end do
    end do

    ! calculate wcc = c_mu,i * W_ij * c_j,nu. We have only W_ab b > a and W_ij j > i:
    ! wcc(m,n) = sum_{pq, p <= q} w_pq (grndEigVecs(mu,p)grndEigVecs(nu,q)
    ! + grndEigVecs(nu,p)grndEigVecs(mu,q))
    ! complexity norb * norb * norb

    ! calculate the occ-occ part
    wcc(:,:) = 0.0_dp

    do ij = 1, nxoo
      call indxoo(homo, nocc, ij, i, j)
      do mu = 1, norb
        do nu = 1, norb
          wcc(mu,nu) = wcc(mu,nu) + woo(ij) *&
              & ( grndEigVecs(mu,i,1)*grndEigVecs(nu,j,1)&
              & + grndEigVecs(mu,j,1)*grndEigVecs(nu,i,1) )
        end do
      end do

    end do

    ! calculate the occ-virt part : the same way as for xpycc
    do ia = 1, nxov
      call indxov(win, ia, getij, i, a)
      do nu = 1, norb
        do mu = 1, norb
          wcc(mu,nu) = wcc(mu,nu) + wov(ia) *&
              & ( grndEigVecs(mu,i,1)*grndEigVecs(nu,a,1)&
              & + grndEigVecs(mu,a,1)*grndEigVecs(nu,i,1) )
        end do
      end do
    end do

    ! calculate the virt - virt part
    do ab =1, nxvv
      call indxvv(homo, ab, a, b)
      do mu = 1, norb
        do nu = 1, norb
          wcc(mu,nu) = wcc(mu,nu) + wvv(ab) *&
              & ( grndEigVecs(mu,a,1)*grndEigVecs(nu,b,1)&
              & + grndEigVecs(mu,b,1)*grndEigVecs(nu,a,1) )
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

        tmp3a = dq(iAt1) * dqex(iAt2) + dqex(iAt1) * dq(iAt2)

        if (sym == "S") then
          tmp3b = 4.0_dp * xpyq(iAt1) * xpyq(iAt2)
        else
          tmp3b = 0.0_dp
        end if

        excgrad(:,iAt1) = excgrad(:,iAt1) + dgab(:) * ( tmp3a + tmp3b )
        excgrad(:,iAt2) = excgrad(:,iAt2) - dgab(:) * ( tmp3a + tmp3b )

        tmp5 = shift_excited(iAt1) + shift_excited(iAt2)
        tmp7 = 2.0_dp * ( shxpyq(iAt1) + shxpyq(iAt2) )

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

          do mu = indalpha, indalpha1
            do nu = indbeta, indbeta1
              m = mu - indalpha + 1
              n = nu - indbeta + 1

              tmp1 = tmp1 + 2.0_dp * dH0(n,m,xyz) * pc(mu,nu)
              tmp2 = tmp2 + dS(n,m,xyz) * pc(mu,nu) * (shift(iAt1)+shift(iAt2))
              tmp3 = tmp3 - dS(n,m,xyz) * wcc(mu,nu)
              tmp4 = tmp4 + tmp5 * dS(n,m,xyz) * rhoSqr(mu,nu)
              tmp6 = tmp6 + tmp7 * dS(n,m,xyz) * xpycc(mu,nu)
            end do
          end do
          excgrad(xyz,iAt1) = excgrad(xyz,iAt1)&
              & + tmp1 + tmp2 + tmp4 + tmp6 + tmp3
          excgrad(xyz,iAt2) = excgrad(xyz,iAt2)&
              & - tmp1 - tmp2 - tmp4 - tmp6 - tmp3
        end do
      end do
    end do

  end subroutine addGradients


  !> Write out excitations projected onto ground state
  subroutine writeCoeffs(tt, grndEigVecs, occ, nocc, fdCoeffs, tCoeffs, tIncGroundState,&
      & occNatural, naturalOrbs)

    !> T part of the matrix
    real(dp), intent(in) :: tt(:,:)

    !> ground state eigenvectors
    real(dp), intent(in) :: grndEigVecs(:,:,:)

    !> ground state occupations
    real(dp), intent(in) :: occ(:,:)

    !> number of filled states
    integer, intent(in) :: nocc

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

    real(dp), allocatable :: t2(:,:), occtmp(:)
    integer :: norb, ii, jj, mm

    norb = size(tt, dim=1)

    if (present(occNatural).or.tCoeffs) then

      ALLOCATE(t2(norb, norb))
      t2 = tt
      if (tIncGroundState) then
        do ii = 1, nocc
          t2(ii,ii) = t2(ii,ii) + occ(ii,1)
        end do
      end if

      if (present(occNatural)) then
        naturalOrbs(:,:,1) = t2
        call evalCoeffs(naturalOrbs(:,:,1) ,occNatural,grndEigVecs(:,:,1))
        if (tCoeffs) then
          ALLOCATE(occtmp(size(occ)))
          occTmp = occNatural
        end if
      else
        ALLOCATE(occtmp(size(occ)))
        occtmp = 0.0_dp
        call evalCoeffs(t2,occNatural,grndEigVecs(:,:,1))
      end if

      ! Better to get this by post-processing DFTB+ output, but here for
      ! compatibility at the moment
      if (tCoeffs) then
        open(fdCoeffs, file=excitedCoefsOut, position="append")
        write(fdCoeffs,*) 'T F'
        do ii = 1, norb
          jj = norb - ii + 1
          write(fdCoeffs, '(1x,i3,1x,f13.10,1x,f13.10)') ii, occtmp(jj), 2.0_dp
          write(fdCoeffs, '(6(f13.10,1x))') (cmplx(t2(mm,jj), kind=dp),&
              & mm = 1, norb)
        end do
        close(fdCoeffs)
      end if

    end if

  end subroutine writeCoeffs


  !> Project MO density matrix onto ground state orbitals
  subroutine evalCoeffs(t2,occ,eig)

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
      & fdTrans, fdTradip, transitionDipoles, tWriteTagged, fdTagged, fdExc, Ssq)

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

    !> file unit for excitation energies
    integer, intent(in) :: fdExc

    !> For spin polarized systems, measure of spin
    real(dp), intent(in), optional :: Ssq(:)

    integer :: nmat
    integer :: i, j, iweight, indo, m, n
    integer :: iDeg
    real(dp), allocatable :: wvec(:)
    real(dp), allocatable :: xply(:)
    real(dp), allocatable :: eDeg(:)
    real(dp), allocatable :: oDeg(:)
    integer, allocatable :: wvin(:)
    real(dp) :: rsqw, weight, wvnorm
    logical :: updwn, tSpin
    character :: sign

    @:ASSERT(fdExc > 0)

    tSpin = present(Ssq)
    nmat = size(wij)
    ALLOCATE(wvec(nmat))
    ALLOCATE(wvin(nmat))
    ALLOCATE(xply(nmat))
    ALLOCATE(eDeg(nexc))
    ALLOCATE(oDeg(nexc))
    wvec = 0.0_dp
    wvin = 0
    xply = 0.0_dp
    eDeg = 0.0_dp
    oDeg = 0.0_dp

    if(fdXplusY > 0) then
      write(fdXplusY,*) nmat, nexc
    end if

    do i = 1, nexc
      if (eval(i) > 0.0_dp) then

        ! calculate weight of single particle transitions
        rsqw = 1.0_dp / sqrt(eval(i))
        ! (X+Y)^ia_I = sqrt(wij) / sqrt(omega) * F^ia_I
        xply(:) = sqrt(rsqw) * sqrt(wij(:)) * evec(:,i)
        wvec(:) = xply(:)**2
        wvnorm = 1.0_dp / sqrt(sum(wvec**2))
        wvec(:) = wvec(:) * wvnorm

        ! find largest coefficient in CI - should use maxloc
        call index_heap_sort(wvin,wvec)
        wvin = wvin(size(wvin):1:-1)
        wvec = wvec(wvin)

        weight = wvec(1)
        iweight = wvin(1)

        call indxov(win, iweight, getij, m, n)
        sign = sym
        if (tSpin) then
          sign = " "
          write(fdExc,&
              & '(1x,f10.3,4x,f14.8,2x,i5,3x,a,1x,i5,7x,f6.3,2x,f10.3,4x,&
              & f6.3)')&
              & Hartree__eV * sqrt(eval(i)), osz(i), m, '->', n, weight,&
              & Hartree__eV * wij(iweight), Ssq(i)
        else
          write(fdExc,&
              & '(1x,f10.3,4x,f14.8,5x,i5,3x,a,1x,i5,7x,f6.3,2x,f10.3,6x,a)')&
              & Hartree__eV * sqrt(eval(i)), osz(i), m, '->', n, weight,&
              & Hartree__eV * wij(iweight), sign
        end if

        if(fdXplusY > 0) then
          if (tSpin) then
            updwn = (win(iweight) <= nmatup)
            sign = "D"
            if (updwn) sign = "U"
          end if
          write(fdXplusY,'(1x,i5,3x,a,3x,ES17.10)') i,sign, sqrt(eval(i))
          write(fdXplusY,'(6(1x,ES17.10))') xply(:)
        endif

        if (fdTrans > 0) then
          write(fdTrans, '(2x,a,T12,i5,T21,ES17.10,1x,a,2x,a)')&
              & 'Energy ', i,  Hartree__eV * sqrt(eval(i)), 'eV', sign
          write(fdTrans,*)
          write(fdTrans,'(2x,a,9x,a,8x,a)')&
              & 'Transition', 'Weight', 'KS [eV]'
          write(fdTrans,'(1x,45("="))')

          sign = " "
          do j = 1, nmat
            !if (wvec(j) < 1e-4_dp) exit ! ??????
            indo = wvin(j)
            call indxov(win, indo, getij, m, n)
            if (tSpin) then
              updwn = (win(indo) <= nmatup)
              sign = "D"
              if (updwn) sign = "U"
            end if
            write(fdTrans,&
                & '(i5,3x,a,1x,i5,1x,1a,T22,f10.8,T33,f14.8)')&
                & m, '->', n, sign, wvec(j), Hartree__eV * wij(wvin(j))
          end do
          write(fdTrans,*)
        end if

        if(fdTradip > 0) then
          write(fdTradip, '(1x,i5,1x,f10.3,2x,3(ES13.6))')&
              & i, Hartree__eV * sqrt(eval(i)), (transitionDipoles(i,j)&
              & * au__Debye, j=1,3)
        endif
      else

        ! find largest coefficient in CI - should use maxloc
        call index_heap_sort(wvin,wvec)
        wvin = wvin(size(wvin):1:-1)
        wvec = wvec(wvin)

        weight = wvec(1)
        iweight = wvin(1)
        call indxov(win, iweight, getij, m, n)
        sign = sym

        if (tSpin) then
          sign = " "
          write(fdExc,&
              & '(6x,A,T12,4x,f14.8,2x,i5,3x,a,1x,i5,7x,A,2x,f10.3,4x,f6.3)')&
              & '< 0', osz(i), m, '->', n, '-', Hartree__eV * wij(iweight),&
              & Ssq(i)
        else
          write(fdExc,&
              & '(6x,A,T12,4x,f14.8,2x,i5,3x,a,1x,i5,7x,f6.3,2x,f10.3,6x,a)')&
              & '< 0', osz(i), m, '->', n, weight,&
              & Hartree__eV * wij(iweight), sign
        end if

        if(fdXplusY > 0) then
          if (tSpin) then
            updwn = (win(iweight) <= nmatup)
            sign = "D"
            if (updwn) sign = "U"
          end if
          write(fdXplusY,'(1x,i5,3x,a,3x,A)') i,sign, '-'
        endif

        if (fdTrans > 0) then
          write(fdTrans, '(2x,a,1x,i5,5x,a,1x,a,3x,a)')&
              & 'Energy ', i,  '-', 'eV', sign
          write(fdTrans,*)
        end if

        if(fdTradip > 0) then
          write(fdTradip, '(1x,i5,1x,A)') i, '-'
        endif

      end if

    end do

    ! Determine degenerate levels and sum oscillator strength over any degenerate levels
    iDeg = 1
    eDeg(1) = eval(1)
    oDeg(1) = osz(1)
    do i = 2, nexc
      if(abs(eval(i)-eval(i-1)) < elecTolMax) then
        oDeg(iDeg) = oDeg(iDeg) + osz(i)
      else
        iDeg = iDeg + 1
        eDeg(iDeg) = eval(i)
        oDeg(iDeg) = osz(i)
      endif
    end do
    if (tWriteTagged) then
      call writeTagged(fdTagged, tag_excEgy, eDeg(:iDeg))
      call writeTagged(fdTagged, tag_excOsc, oDeg(:iDeg))
    end if

  end subroutine writeExcitations


  !> Create transition density matrix in MO basis P = T + 1/2 Z symmetric (paper has T + Z
  !> asymmetric) (Zab = Zij = 0, Tia = 0)
  subroutine calcPMatrix(t, rhs, win, getij, pc)

    !> T matrix
    real(dp), intent(in) :: t(:,:)

    !> Z matrix
    real(dp), intent(in) :: rhs(:)

    !> index array for single particle transitions
    integer, intent(in) :: win(:)

    !> array of the occupied->virtual pairs (nTransitions,occ 1 or virtual 2)
    integer, intent(in) :: getij(:,:)

    !> resulting excited state density matrix
    real(dp), intent(out) :: pc(:,:)

    integer :: ia, i, a

    pc = 0.0_dp
    do ia = 1, size(rhs)
      call indxov(win, ia, getij, i, a)
      pc(i,a) = rhs(ia)
    end do
    pc = 0.5_dp * ( pc + transpose(pc) )

    pc = pc + t

  end subroutine calcPMatrix


  !> Write single particle excitations to a file as well as potentially to tagged output file (in
  !> that case, summing over degeneracies)
  subroutine writeSPExcitations(wij, win, nmatup, getij, fdSPTrans, sposz, nxov, tSpin)

    !> single particle excitation energies
    real(dp), intent(in) :: wij(:)

    !> index array for single particle transitions
    integer, intent(in) :: win(:)

    !> number of transitions within same spin channel
    integer, intent(in) :: nmatup

    !> index from composite index to occupied and virtual single particle states
    integer, intent(in) :: getij(:,:)

    !> file descriptor for the single particle excitation data
    integer, intent(in) :: fdSPTrans

    !> single particle oscilation strengths
    real(dp), intent(in) :: sposz(:)

    !> Number of included single particle excitations to print out (assumes that win and wij are
    !> sorted so that the wanted transitions are first in the array)
    integer, intent(in) :: nxov

    !> is this a spin-polarized calculation?
    logical, intent(in) :: tSpin

    integer :: indm, m, n
    logical :: updwn
    character :: sign

    @:ASSERT(size(sposz)>=nxov)

    if (fdSPTrans > 0) then
      ! single particle excitations
      open(fdSPTrans, file=singlePartOut, position="rewind", status="replace")
      write(fdSPTrans,*)
      write(fdSPTrans,'(7x,a,7x,a,8x,a)') '#      w [eV]',&
          & 'Osc.Str.', 'Transition'
      write(fdSPTrans,*)
      write(fdSPTrans,'(1x,58("="))')
      write(fdSPTrans,*)
      do indm = 1, nxov
        call indxov(win, indm, getij, m, n)
        sign = " "
        if (tSpin) then
          updwn = (win(indm) <= nmatup)
          if (updwn) then
            sign = "U"
          else
            sign = "D"
          end if
        end if
        write(fdSPTrans,&
            & '(1x,i7,3x,f8.3,3x,f13.7,4x,i5,3x,a,1x,i5,1x,1a)')&
            & indm, Hartree__eV * wij(indm), sposz(indm), m, '->', n, sign
      end do
      write(fdSPTrans,*)
      close(fdSPTrans)
    end if

  end subroutine writeSPExcitations

end module linrespgrad
