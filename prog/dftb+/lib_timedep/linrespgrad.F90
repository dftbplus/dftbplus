!--------------------------------------------------------------------------------------------------!
!  DFTB+: general package for performing fast atomistic simulations                                !
!  Copyright (C) 2017  DFTB+ developers group                                                      !
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
  use scc, only : getAtomicGammaMatrix
  use blasroutines
  use eigensolver
  use message
  use taggedoutput
  use sorting
  use qm
  !use bisect
  implicit none
  private

  public :: LinRespGrad_old

  ! Tolerance for ARPACK solver.
  real(dp), parameter :: ARTOL = epsilon(1.0_dp)

  ! Maximal allowed iteration in the ARPACK solver.
  integer, parameter :: MAX_AR_ITER = 300

  character(lc) :: tmpStr

  ! Names of output files
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

  !! Communication with ARPACK (debug.h) for progress information
  integer :: logfil, ndigit, mgetv0
  integer :: msaupd, msaup2, msaitr, mseigt, msapps, msgets, mseupd
  integer :: mnaupd, mnaup2, mnaitr, mneigh, mnapps, mngets, mneupd
  integer :: mcaupd, mcaup2, mcaitr, mceigh, mcapps, mcgets, mceupd
  common /debug/ logfil, ndigit, mgetv0,&
      &    msaupd, msaup2, msaitr, mseigt, msapps, msgets, mseupd,&
      &    mnaupd, mnaup2, mnaitr, mneigh, mnapps, mngets, mneupd,&
      &    mcaupd, mcaup2, mcaitr, mceigh, mcapps, mcgets, mceupd

contains

  !> This subroutine analytically calculates excitations and gradients
  !! of excited state energies based on Time Dependent DFRT
  !! \param tSpin spin polarized calculation
  !! \param natom number of atoms
  !! \param iAtomStart index vector for S and H matrices
  !! \param ndim dimension of S,H
  !! \param grndEigVecs ground state MO-coefficients
  !! \param grndEigVal ground state MO-energies
  !! \param dq converged ground state Mulliken net charges - atomic charges
  !! \param coord0 atomic positions
  !! \param nexc number of excited states to solve for
  !! \param nstat0 state of interest (< 0 find brightest, 0 calculate
  !! all nexc states, > 0 that specific state)
  !! \param symc symmetry required singlet ('S'), triplet ("T") or both ("B")
  !! \param SSqr square overlap matrix between basis functions, both
  !! triangles required
  !! \param species0 chemical species of each atom
  !! \param HubbardU ground state Hubbard U values for each species
  !! \param spinW ground state spin derivatives for each species
  !! \param rnel real number of electrons in system
  !! \param iNeighbor Atomic neighbour lists
  !! \param img2CentCell Mapping of atom number to central cell atom number
  !! \param orb data type for atomic orbital information
  !! \param tWriteTagged print tag information
  !! \param fdTagged file descriptor for the tagged data output
  !! \param fdMulliken file unit for excited Mulliken populations?
  !! \param fdCoeffs file unit if the coefficients for the excited
  !! states should be written to disc
  !! \param fdXplusY file for X+Y data
  !! \param fdTrans File unit for transitions if required
  !! \param fdSPTrans File unit for single particle (KS) transitions
  !! if required
  !! \param fdTradip File unit for transition dipole strengths
  !! \param tArnoldi write state of Arnoldi solver to disc
  !! \param fdArnoldi file unit for Arnoldi write out
  !! \param fdArnoldiDiagnosis file unit for Arnoldi solver tests, if
  !! this is < 1 no tests are performed
  !! \param tEnergyWindow is an energy window specified
  !! \param energyWindow window for transitions above nstat0 (if nstat0 > 0)
  !! \param tOscillatorWindow is an oscillator window window specified
  !! \param oscillatorWindow window for transitions above nstat0 (if nstat0 > 0)
  !! \param omega energy of state nstat0 if nstat0 > 0
  !! \param shift shift vector for potentials
  !! \param skHamCont non-SCC hamitonian data
  !! \param skOverCont overlap data
  !! \param excgrad excitation energy gradient
  !! \param derivator Differentiatior for the non-scc components.
  !! \param rhoSqr ground state square density matrix
  !! \param occNatural Occupation numbers for natural orbitals from excited
  !! state density matrix
  !! \param naturalOrbs the single particle eigenvectors themselves.
  subroutine LinRespGrad_old(tSpin, natom, iAtomStart, &
      & grndEigVecs, grndEigVal, dq, coord0, nexc, nstat0, symc, SSqr, filling,&
      & species0, HubbardU, spinW, rnel, iNeighbor, img2CentCell, &
      & orb, tWriteTagged, fdTagged, fdMulliken, fdCoeffs, tGrndState, &
      & fdXplusY, fdTrans, fdSPTrans, fdTradip, tArnoldi, fdArnoldi, &
      & fdArnoldiDiagnosis, fdExc,tEnergyWindow, energyWindow,tOscillatorWindow, oscillatorWindow, &
      & omega, tGrads, shift, skHamCont, skOverCont, excgrad, derivator, rhoSqr, &
      & occNatural, naturalOrbs)
    logical, intent(in) :: tSpin
    integer, intent(in) :: natom, iAtomStart(:)
    real(dp), intent(in) :: grndEigVecs(:,:,:), grndEigVal(:,:)
    real(dp), intent(in) :: dq(:), coord0(:,:)
    integer, intent(in) :: nexc, nstat0
    character, intent(in) :: symc
    real(dp), intent(in) :: SSqr(:,:)
    real(dp), intent(in) :: filling(:,:)
    integer, intent(in) :: species0(:)
    real(dp), intent(in) :: HubbardU(:), spinW(:), rnel
    integer, intent(in) :: iNeighbor(0:,:), img2CentCell(:)
    type(TOrbitals), intent(in) :: orb
    logical, intent(in) :: tWriteTagged
    integer, intent(in) :: fdTagged, fdMulliken, fdCoeffs
    logical, intent(in) :: tGrndState
    integer, intent(in) :: fdTradip
    integer, intent(in) :: fdXplusY
    integer, intent(in) :: fdTrans
    integer, intent(in) :: fdSPTrans
    logical, intent(in) :: tArnoldi
    integer, intent(in) :: fdArnoldi
    integer, intent(in) :: fdArnoldiDiagnosis
    integer, intent(in) :: fdExc
    logical, intent(in)  :: tEnergyWindow
    real(dp), intent(in) :: energyWindow
    logical, intent(in)  :: tOscillatorWindow
    real(dp), intent(in) :: oscillatorWindow
    real(dp), intent(out) :: omega
    logical, intent(in)  :: tGrads
    real(dp), intent(in), optional :: shift(:)
    type(OSlakoCont), intent(in), optional :: skHamCont, skOverCont
    real(dp), intent(out), optional :: excgrad(:,:)
    class(NonSccDiff), intent(in), optional :: derivator
    real(dp), intent(in), optional :: rhoSqr(:,:,:)
    real(dp), intent(out), optional :: occNatural(:)
    real(dp), intent(out), optional :: naturalOrbs(:,:)

    real(dp) :: Ssq(nexc)
    real(dp), allocatable :: gammaMat(:,:), snglPartTransDip(:,:)
    real(dp), allocatable :: stimc(:,:,:), wij(:)
    real(dp), allocatable :: dqex(:), sposz(:), osz(:), xpy(:), xmy(:), pc(:,:)
    real(dp), allocatable :: t(:,:), rhs(:), woo(:), wvv(:), wov(:)
    real(dp), allocatable :: evec(:,:), eval(:), transitionDipoles(:,:)
    integer, allocatable :: win(:), iatrans(:,:), getij(:,:)
    character, allocatable :: symmetries(:)

    integer :: nocc, nocc_r, nvir_r, nxoo_r, nxvv_r
    integer :: nxov, nxov_ud(2), nxov_r, nxov_d, nxov_rd
    integer :: norb
    integer :: i, j, iSpin, isym
    integer :: nSpin
    character :: sym

    real(dp) :: energyThreshold

    integer :: nStat

    ! control variables
    logical :: tZVector, tCoeffs, tTradip

    ! printing data
    logical :: tMulliken

    ndigit = -3
    !! Output unit:
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
    !! End of ARPACK communication

    @:ASSERT(fdExc > 0)

    ! work out which data files are required, based on whether they
    ! have valid file IDs (>0)
    tMulliken = (fdMulliken > 0)
    tCoeffs = (fdCoeffs > 0)
    tTradip = (fdTradip > 0)

    @:ASSERT(fdArnoldi > 0)
    if (tArnoldi) then
      open(fdArnoldi, file=arpackOut, position="rewind", status="replace")
    end if

    nSpin = size(grndEigVal, dim=2)
    @:ASSERT(nSpin > 0 .and. nSpin <=2)

    norb = orb%nOrb

    @:ASSERT(present(excgrad) .eqv. present(shift))
    @:ASSERT(present(shift) .eqv. present(skHamCont))
    @:ASSERT(present(excgrad) .eqv. present(rhoSqr))
    @:ASSERT(present(skHamCont) .eqv. present(skOverCont))
    @:ASSERT(present(excgrad) .eqv. present(derivator))

    @:ASSERT(present(occNatural) .eqv. present(naturalOrbs))

    ! count initial number of transitions from occupied to empty states
    nxov_ud = 0
    do iSpin = 1, nSpin
      do i = 1, norb - 1
        do j = i, norb
          if (filling(i,iSpin) > filling(j,iSpin) + elecTolMax) then
            nxov_ud(iSpin) = nxov_ud(iSpin) + 1
          end if
        end do
      end do
    end do
    nxov = sum(nxov_ud)

    if (nexc + 1 >= nxov) then
      write(tmpStr,"(' Insufficent single particle excitations, ',I0, &
          &', for required number of excited states ',I0)")nxov, nexc
      call error(tmpStr)
    end if

    tZVector = tGrads .or. tMulliken .or. tCoeffs .or. &
        & present(naturalOrbs)

    ! Sanity checks
    nstat = nstat0
    if (nstat < 0 .and. symc /= "S") then
      call error("Linresp: Brightest mode only available for singlets.")
    end if
    if (nstat /= 0 .and. symc == "B") then
      call error("Linresp: Both symmetries not allowed if a specific state is excited")
    end if
    if (nstat == 0 .and. tZVector) then
      call error("Linresp: Properties calculation only available  with one selected excited state.")
    end if
    if (nstat == 0 .and. tMulliken) then
      call error("Linresp: Excited charges available only with one selected excited state.")
    end if
    if (nstat == 0 .and. tCoeffs) then
      call error("Linresp: Coefficients only available with one selected excited state.")
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
    else ! ADG: temporary solution for spin polarized case.
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

    ! Overlap times wave function coefficients - most routines in
    ! DFTB+ use lower triangle (would remove need to symmetrize
    ! overlap in main code)
    do iSpin = 1, nSpin
      call symm(stimc(:,:,iSpin), "L", SSqr, grndEigVecs(:,:,iSpin))
    end do

    ! ground state Hubbard U softened coulombic interactions
    call getAtomicGammaMatrix(gammaMat, iNeighbor, img2CentCell)

    ! Oscillator strength, when needed.
    ALLOCATE(osz(nexc))

    ! Find all transitions and KS energy differences for cases that go
    ! from filled to empty states
    call getSPExcitations(grndEigVal, filling, wij, getij)

    ! put them in ascending energy order
    if (tOscillatorWindow) then
      ! use a stable sort so that degenerate transitions from the same
      ! single particle state are grouped together in the results
      call merge_sort(win,wij, 1.0_dp*epsilon(1.0))
    else
      ! do not require stability, use usual routine to sort
      call index_heap_sort(win,wij)
    end if
    wij = wij(win)

    ! dipole strength of transitions between K-S states
    call calcTransitionDipoles(coord0, win, nxov_ud(1), getij, &
        & iAtomStart, stimc, grndEigVecs, snglPartTransDip)

    sposz(:) = twothird * wij(:) * sum(snglPartTransDip**2, dim=2)

    if (tOscillatorWindow .and. tZVector ) then
      call error("Incompabilitity between excited state property&
          & evaluation and an oscillator strength window at the moment.")
    end if

    if (tOscillatorWindow .or. tEnergyWindow) then

      if (.not. tEnergyWindow) then

        ! find transitions that are strongly dipole allowed (> oscillatorWindow)
        call dipselect(wij, sposz, win, snglPartTransDip,nxov_rd, &
            & oscillatorWindow, grndEigVal, getij)

      else

        energyThreshold = wij(nexc) + energyWindow
        nxov_r = count(wij <= energyThreshold)

        nxov_d = 0
        if (tOscillatorWindow) then

          ! find transitions that are strongly dipole allowed (> oscillatorWindow)
          call dipselect(wij(nxov_r+1:), sposz(nxov_r+1:), win(nxov_r+1:), &
              & snglPartTransDip(nxov_r+1:,:),nxov_d, oscillatorWindow, &
              & grndEigVal, getij)

        end if

        nxov_rd = nxov_r + nxov_d

      end if
    else

      nxov_rd = nxov

    end if

    ! just in case energy/dipole windows add no extra states, and is due to an
    ! arpack solver requirement combined with the need to get at least nexc
    ! states
    nxov_rd = max(nxov_rd,min(nexc+1,nxov))

    if (fdXplusY >  0) then
      open(fdXplusY, file=XplusYOut, position="rewind", status="replace")
    end if

    if(fdTrans>0) then
      open(fdTrans, file=transitionsOut, position="rewind", status="replace")
      write(fdTrans,*)
    endif

    if (fdTradip > 0) then ! single particle transition dipole file
      open(fdTradip, file=transDipOut, position="rewind", status="replace")
      write(fdTradip,*)
      write(fdTradip,'(5x,a,5x,a,7x,a,6x,a,6x,a,6x,a)') "#", 'w [eV]',&
          & 'Transition dipole (x,y,z) [Debye]'
      write(fdTradip,*)
      write(fdTradip,'(1x,65("="))')
      write(fdTradip,*)
    endif

    ! excitation energies
    open(fdExc, file=excitationsOut, position="rewind", status="replace")
    write(fdExc,*)
    if (tSpin) then
      write(fdExc,'(5x,a,7x,a,9x,a,9x,a,6x,a,4x,a)') 'w [eV]',&
          & 'Osc.Str.', 'Transition','Weight','KS [eV]','D<S*S>'
    else
      write(fdExc,'(5x,a,7x,a,9x,a,9x,a,6x,a,4x,a)') 'w [eV]',&
          & 'Osc.Str.', 'Transition','Weight','KS [eV]','Sym.'
    end if

    write(fdExc,*)
    write(fdExc,'(1x,80("="))')
    write(fdExc,*)

    ! single particle excitations (output file and tagged file if
    ! needed).  Was used for nxov_rd = size(wij), but now for just
    ! states that are actually included in the excitation calculation.
    call writeSPExcitations(wij, win, nxov_ud(1), getij, fdSPTrans, sposz, nxov_rd, tSpin)
    ALLOCATE(evec(nxov_rd, nexc))

    do isym = 1, size(symmetries)

      sym = symmetries(isym)
      call buildAndDiagExcMatrix(tSpin, wij(:nxov_rd), sym, win, &
          & nxov_ud(1), nxov_rd, iAtomStart, stimc, grndEigVecs, filling, &
          & getij, gammaMat, species0, spinW, fdArnoldiDiagnosis, eval, evec )
      ! Excitation oscillator strengths for resulting states
      call getOscillatorStrengths(sym, snglPartTransDip(1:nxov_rd,:), &
          & wij(:nxov_rd), eval, evec, filling, win, nxov_ud(1), getij, &
          & nstat, osz, tTradip, transitionDipoles)
      if (tSpin) then
        call getExcSpin(Ssq, nxov_ud(1), getij, win, eval, evec, &
            & wij(:nxov_rd), filling, stimc, grndEigVecs)
        call writeExcitations(sym, osz, nexc, nxov_ud(1), getij, win, &
            & eval, evec, wij(:nxov_rd), fdXplusY, fdTrans, &
            & fdTradip, transitionDipoles, tWriteTagged, fdTagged, fdExc, Ssq)
      else
        call writeExcitations(sym, osz, nexc, nxov_ud(1), getij, win, &
            & eval, evec, wij(:nxov_rd), fdXplusY, fdTrans, &
            & fdTradip, transitionDipoles, tWriteTagged, fdTagged, fdExc)
      end if
    end do


    if (tArnoldi) then
      close(fdArnoldi)
    end if

    if (fdTrans > 0) close(fdTrans)
    if (fdXplusY > 0) close(fdXplusY)
    if (fdExc > 0) close(fdExc)
    if (fdTradip > 0) close(fdTradip)

    if (nstat == 0) then
      omega = 0.0_dp
      return
    end if

    ! Remove some un-used memory
    deallocate(snglPartTransDip)
    deallocate(transitionDipoles)
    deallocate(sposz)

    omega = sqrt(eval(nstat))


    if (tZVector) then ! calculate Furche vectors and transition
                       ! density matrix for various properties

      if (tSpin) then
        call error("Z vector evaluation does not currently support spin&
            & polarization.")
      end if

      if (any( abs(filling) > elecTolMax .and. &
          & abs(filling-2.0_dp) > elecTolMax ) ) then
        call error("Fractional fillings not currently possible for excited &
            &state property calculations")
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


      ! solve for Z and W to get excited state density matrix
      call getZVectorEqRHS(xpy, xmy, win, iAtomStart, nocc, nocc_r, &
          & nxov_ud(1), getij, iatrans, natom, species0,grndEigVal(:,1), &
          & stimc, grndEigVecs, gammaMat, spinW, omega, sym, rhs, t, &
          & wov, woo, wvv)
      call solveZVectorEq(rhs, win, nxov_ud(1), getij, natom, iAtomStart, &
          & stimc, gammaMat, wij(:nxov_rd), grndEigVecs)
      call calcWVectorZ(rhs, win, nocc, nocc_r, nxov_ud(1), getij, iAtomStart, &
          & stimc, grndEigVecs, gammaMat, grndEigVal(:,1), wov, woo, wvv)
      call calcPMatrix(t, rhs, win, getij, pc)

      call writeCoeffs(pc, grndEigVecs, filling, nocc, fdCoeffs, &
          & tCoeffs, tGrndState, occNatural, naturalOrbs)

      ! Make MO to AO transformation of the excited density matrix
      call unitary(pc, grndEigVecs(:,:,1))

      call getExcMulliken(iAtomStart, pc, SSqr, dqex)
      if (tMulliken) then
        call writeExcMulliken(sym, nstat, dq, dqex, coord0, fdMulliken)
      end if


      if (tGrads) then
        call addGradients(sym, nxov_rd, natom, species0, iAtomStart, norb, &
            & nocc, nocc_r, nxov_ud(1), getij, win, grndEigVecs, pc, stimc, &
            & dq, dqex, gammaMat, HubbardU, spinW, shift, woo, wov, wvv, &
            & xpy, coord0, orb, skHamCont, skOverCont, derivator, &
            & rhoSqr(:,:,1), excgrad)
      end if

    end if

  end subroutine LinRespGrad_old

  !> Builds and diagonalizes the excitation matrix via iterative technique.
  !! \param tSpin spin polarisation?
  !! \param wij single particle excitation energies
  !! \param sym symmetry to calculate transitions
  !! \param win index array for single particle excitions
  !! \param nmatup number of same spin excitations
  !! \param nxov number of occupied-virtual transitions
  !! \param iAtomStart square atom array indexing
  !! \param stimc overlap times ground state eigenvectors
  !! \param grndEigVecs ground state eigenvectors
  !! \param filling occupation numbers
  !! \param getij index array between transitions in square and 1D
  !! representations
  !! \param gammaMat electrostatic matrix
  !! \param species0 central cell chemical species
  !! \param spinW ground state spin derivatives for each species
  !! \param eval resulting eigenvalues for transitions
  !! \param evec eigenvectors for transitions
  subroutine buildAndDiagExcMatrix(tSpin, wij, sym, win, nmatup, &
      & nxov, iAtomStart, stimc, grndEigVecs, filling, getij, gammaMat, &
      & species0, spinW, fdArnoldiDiagnosis, &
      & eval, evec)
    logical, intent(in) :: tSpin
    real(dp), intent(in) :: wij(:)
    character, intent(in) :: sym
    integer, intent(in) :: win(:), nmatup, nxov, iAtomStart(:)
    real(dp), intent(in) :: stimc(:,:,:), grndEigVecs(:,:,:)
    real(dp), intent(in) :: filling(:,:)
    real(dp), intent(in) :: gammaMat(:,:)
    integer, intent(in) :: getij(:,:), species0(:)
    integer, intent(in) :: fdArnoldiDiagnosis
    real(dp), intent(in) :: spinW(:)
    real(dp), intent(out) :: eval(:), evec(:,:)

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

    info = 0 ! random initial vector used for dsaupd
    ido = 0 ! IDO must be zero on the first  call
    iparam(1) = 1 ! restarting the iteration with a starting vector
                  ! that is a linear combination of Ritz vectors
                  ! associated with the "wanted" Ritz values.
    iparam(3) = MAX_AR_ITER ! maximum iterations of solver
    iparam(7) = 1 ! solve A*x = lambda*x, with A symmetric

    do ! loop until exit

      ! call the reverse communication interface from arpack
      call saupd (ido, "I", nxov, "SM", nexc, ARTOL, resid, ncv, vv, nxov,&
          & iparam, ipntr, workd, workl, lworkl, info)

      if (ido == 99) then ! has terminated normally
        exit
      end if
      if (abs(ido) /= 1) then
        write(tmpStr,"(' Unexpected return from arpack routine saupd, IDO ',I0,&
            & ' INFO ',I0)") ido,info
        call error(tmpStr)
      end if

      ! Action of excitation supermatrix on supervector
      call omegatvec(tSpin, workd(ipntr(1):ipntr(1)+nxov-1),&
          & workd(ipntr(2):ipntr(2)+nxov-1),&
          & wij, sym, win, nmatup, iAtomStart, stimc,&
          & grndEigVecs, filling, getij, gammaMat, species0,&
          & spinW)

    end do

    if (info < 0) then
      write(tmpStr,"(' Error with ARPACK routine saupd, info = ',I0)")info
      call error(tmpStr)
    else if (info  ==  1) then
      call error("Maximum number of iterations reached.&
          & Increase the number of excited states to solve for&
          & (NrOfExcitations).")
    else
      rvec = .true. ! want Ritz vectors

      ! everything after the first 6 variables are passed directly to
      ! DSEUPD following the last call to DSAUPD.  These arguments
      ! MUST NOT BE MODIFIED between the the last call to DSAUPD and
      ! the call to DSEUPD.
      call seupd (rvec, "All", selection, eval, evec, nxov, sigma, "I", nxov,&
          & "SM", nexc, ARTOL, resid, ncv, vv, nxov, iparam, ipntr, workd, &
          & workl, lworkl, info)

      if (info  /=  0) then
        write(tmpStr,"(' Error with ARPACK routine seupd, info = ',I0)")info
        call error(tmpStr)
      end if

    end if

    ! tests for quality of returned eigenpairs
    if (fdArnoldiDiagnosis > 0) then
      open(fdArnoldiDiagnosis, file=testArpackOut, position="rewind", &
          & status="replace")
      ALLOCATE(Hv(nxov))
      ALLOCATE(orthnorm(nxov,nxov))
      orthnorm = matmul(transpose(evec(:,:nExc)),evec(:,:nExc))

      write(fdArnoldiDiagnosis,"(A)")'State Ei deviation    Evec deviation &
          & Norm deviation  Max non-orthog'
      do iState = 1, nExc
        call omegatvec(tSpin, evec(:,iState),&
            & Hv, wij, sym, win, nmatup, iAtomStart, stimc,&
            & grndEigVecs, filling, getij, gammaMat, species0,&
            & spinW)
        write(fdArnoldiDiagnosis,"(I4,4E16.8)")iState, &
            & dot_product(Hv,evec(:,iState))-eval(iState), &
            & sqrt(sum( (Hv-evec(:,iState)*eval(iState) )**2 )), &
            & orthnorm(iState,iState) - 1.0_dp, &
            & max(maxval(orthnorm(:iState-1,iState)), &
            & maxval(orthnorm(iState+1:,iState)))
      end do
      close(fdArnoldiDiagnosis)
    end if

  end subroutine buildAndDiagExcMatrix

  ! Calculate oscillator strength for a given excitation between KS states
  subroutine getOscillatorStrengths(sym, snglPartTransDip, wij, eval, evec, &
      & filling, win, nmatup, getij, istat, osz, tTradip, transitionDipoles)
    character, intent(in) :: sym
    real(dp), intent(in) :: snglPartTransDip(:,:), wij(:), eval(:), evec(:,:)
    real(dp), intent(in) :: filling(:,:)
    integer, intent(in) :: win(:), nmatup, getij(:,:)
    logical :: tTradip
    integer, intent(inout) :: istat
    real(dp), intent(out) :: osz(:), transitionDipoles(:,:)

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

    do ii = 1, size(evec, dim=2)
      osz(ii) = oscillatorStrength(snglPartTransDip, wnij, evec(:,ii))
    end do

    if (istat < 0) then ! find largest transition dipole transition
      oszLoc = maxloc(osz)
      istat = oszLoc(1)
    end if

    if (tTradip) then
      call transitionDipole(snglPartTransDip, wnij, eval, evec, &
          & transitionDipoles)
    end if

  end subroutine getOscillatorStrengths

  !> Calculate excitation energies for closed and open shell systems
  !! \param Ssq <S^2> as a measure of spin contamination (smaller
  !! magnitudes are better, 0.5 is considered an upper threshold for
  !! reliability according to Garcia thesis)
  !! \param nmatup number of same spin excitations
  !! \param getij index for composite excitations to specific occupied
  !! and empty states
  !! \param win index for single particle excitations
  !! \param eval exitation energies
  !! \param evec excited eigenvectors
  !! \param wij single particle excitation energies
  !! \param filling occupations in ground state
  !! \param stimc Overlap times ground state eigenvectors
  !! \param grndEigVecs Ground state eigenvectors
  subroutine getExcSpin(Ssq, nmatup, getij, win, eval, evec, wij, filling, &
      & stimc, grndEigVecs)
    real(dp), intent(out) :: Ssq(:)
    integer, intent(in) :: nmatup, getij(:,:), win(:)
    real(dp), intent(in) :: eval(:), evec(:,:), wij(:)
    real(dp), intent(in) :: filling(:,:)
    real(dp), intent(in) :: stimc(:,:,:), grndEigVecs(:,:,:)

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
        !if (TDvec_sq(k) < 1e-4_dp) exit ! ??????
        ia = TDvin(k)
        call indxov(win, ia, getij, ii, aa)
        ud_ia = (win(ia) <= nmatup)
        do l = 1, nmat
          !if (TDvec_sq(l) < 1e-4_dp) exit ! ??????
          jb = TDvin(l)
          call indxov(win, jb, getij, jj, bb)
          ud_jb = (win(jb) <= nmatup)

          if ( (bb /= aa) .or. (ud_jb .neqv. ud_ia) ) then
            cycle
          end if

          tmp = 0.0_dp
          if (ud_ia) then
            do m = 1,ndwn
              tmp = tmp + MOoverlap(ii,m,stimc,grndEigVecs) &
                  &     * MOoverlap(jj,m,stimc,grndEigVecs)
            end do
          else
            do m = 1,nup
              tmp = tmp + MOoverlap(m,ii,stimc,grndEigVecs) &
                  &     * MOoverlap(m,jj,stimc,grndEigVecs)
            end do
          end if

          s_iaja = s_iaja + TDvec(ia)*TDvec(jb)*tmp

        end do
      end do

      ! S_{ia,ib}
      s_iaib = 0.0_dp
      do k = 1, nmat
        !if (TDvec_sq(k) < 1e-4_dp) exit ! ??????
        ia = TDvin(k)
        call indxov(win, ia, getij, ii, aa)
        ud_ia = (win(ia) <= nmatup)
        do l = 1, nmat
          !if (TDvec_sq(l) < 1e-4_dp) exit ! ??????
          jb = TDvin(l)
          call indxov(win, jb, getij, jj, bb)
          ud_jb = (win(jb) <= nmatup)

          if ( (ii /= jj) .or. (ud_jb .neqv. ud_ia) ) then
            cycle
          end if

          tmp = 0.0_dp
          if (ud_ia) then
            do m = 1,ndwn
              tmp = tmp + MOoverlap(aa,m,stimc,grndEigVecs) &
                  &     * MOoverlap(bb,m,stimc,grndEigVecs)
            end do
          else
            do m = 1,nup
              tmp = tmp + MOoverlap(m,aa,stimc,grndEigVecs) &
                  &     * MOoverlap(m,bb,stimc,grndEigVecs)
            end do
          end if

          s_iaib = s_iaib + TDvec(ia)*TDvec(jb)*tmp
        end do
      end do

      ! S_{ia,jb}
      s_iajb = 0.0_dp
      do k = 1, nmat
        !if (TDvec_sq(k) < 1e-4_dp) exit ! ??????
        ia = TDvin(k)
        call indxov(win, ia, getij, ii, aa)
        ud_ia = (win(ia) <= nmatup)
        if (.not. ud_ia ) then
          cycle
        end if
        do l = 1, nmat
          !if (TDvec_sq(l) < 1e-4_dp) exit ! ??????
          jb = TDvin(l)
          call indxov(win, jb, getij, jj, bb)
          ud_jb = (win(jb) <= nmatup)

          if ( ud_jb ) cycle

          s_iajb = s_iajb + TDvec(ia)*TDvec(jb)&
              & *MOoverlap(aa,bb,stimc,grndEigVecs) &
              & *MOoverlap(ii,jj,stimc,grndEigVecs)

        end do
      end do

      Ssq(i) =  s_iaja - s_iaib - 2.0_dp*s_iajb

    end do

  end subroutine getExcSpin

  ! Build right hand side of the equation for the Z-vector and those parts
  ! of the W-vectors which do not depend on Z.
  subroutine getZVectorEqRHS(xpy, xmy, win, iAtomStart, homo, nocc, nmatup, &
      & getij, iatrans, natom, species0, ev, stimc, c, gammaMat, spinW, &
      & omega, sym, rhs, t, wov, woo, wvv)
    real(dp), intent(in) :: xpy(:), xmy(:)
    integer, intent(in) :: win(:), iAtomStart(:)
    integer, intent(in) :: homo, nocc, nmatup
    integer, intent(in) :: getij(:,:), iatrans(1:,homo+1:), natom, species0(:)
    real(dp), intent(in) :: ev(:), stimc(:,:,:), c(:,:,:), gammaMat(:,:)
    real(dp), intent(in) :: spinW(:), omega
    character, intent(in) :: sym
    real(dp), intent(out) :: rhs(:), t(:,:)
    real(dp), intent(out) :: wov(:), woo(:), wvv(:)

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
        !ab = (b - homo) + (((a - homo) - 1) * (a - homo))/2
        call rindxvv(homo, a, b, ab)
        tmp1 = xpy(ia) * xpy(ib) + xmy(ia) * xmy(ib)
        tmp2 = omega * (xpy(ia) * xmy(ib)+ xmy(ia) * xpy(ib))
        t(a,b) = t(a,b) + 0.5_dp * tmp1
        ! to prevent double counting
        if (a /= b) then
          t(b,a) = t(b,a) + 0.5_dp * tmp1
        end if
        ! Note: diagonal elements will be multiplied by 0.5 later.
        wvv(ab) = wvv(ab) + ev(i) * tmp1 + tmp2
      end do

      ! Build t_ij = 0.5 * sum_a (X+Y)_ia (X+Y)_ja + (X-Y)_ia (X-Y)_ja
      ! and 1 / (1 + delta_ij) Q_ij with Q_ij as in eq. (B9) (1st part of w_ij)
      do j = i, homo
        ja = iatrans(j,a)
        !ij = i-homo+nocc + ((j-homo+nocc - 1) * (j-homo+nocc)) / 2
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
        woo(ij) = woo(ij) - ev(a) * tmp1 + tmp2
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
    ! qgamxpyq(ij) = sum_kb K_ij,kb (X+Y)_kb
    if (sym == "S") then
      do ij = 1, nxoo
        qgamxpyq(ij) = 0.0_dp
        call indxoo(homo, nocc, ij, i, j)
        call transq(i, j, iAtomStart, updwn, stimc, c, qij)
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

  ! Solving the (A+B) Z = -R equation via conjugate gradient optimization.
  !! \param rhs on entry -R, on exit Z
  !! \param win index for single particle excitations
  !! \param nmatup number of transitions between only up states
  !! \param getij index array from composite index to specific
  !! filled-empty transition
  !! \param natom number of atoms
  !! \param iAtomStart index vector for S and H0 matrices
  !! \param stimc overlap times ground state mo-coefficients
  !! \param gammaMat Softened coulomb matrix
  !! \param wij single particle excitation energies
  !! \param c ground state mo-coefficients
  subroutine solveZVectorEq(rhs, win, nmatup, getij, natom, iAtomStart, &
      & stimc, gammaMat, wij, c)
    real(dp), intent(inout) :: rhs(:)
    integer, intent(in) :: win(:), nmatup, getij(:,:), natom, iAtomStart(:)
    real(dp), intent(in) :: stimc(:,:,:), gammaMat(:,:), wij(:), c(:,:,:)

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
    call apbw(rkm1, rhs2, wij, nxov, natom, win, nmatup, getij, iAtomStart, &
        & stimc, c, gammaMat)

    rkm1 = rhs - rkm1
    pkm1 = rkm1

    ! Iteration: should be convergent in at most nxov steps for
    ! quadradic surface, so set higher
    do k = 1, nxov**2
      ! print *, "iteration = ", k, "of", nxov

      ! action of matrix on vector
      call apbw(apk, pkm1, wij, nxov, natom,&
          & win, nmatup, getij, iAtomStart, stimc, c, gammaMat)

      tmp1 = dot_product(rkm1, rkm1)
      tmp2 = dot_product(pkm1, apk)
      alphakm1 = tmp1 / tmp2

      !call axpy(alphakm1, pkm1, rhs2)
      rhs2 = rhs2 + alphakm1 * pkm1
      !call axpy(-alphakm1, apk, rkm1)
      rkm1 = rkm1 -alphakm1 * apk

      tmp2 = dot_product(rkm1, rkm1)
      ! print *, "residuo", tmp2
      if (tmp2 <= epsilon(1.0_dp)**2) then ! 1.0e-14_dp) then
        exit
      end if
      if (k == nxov**2) then
        call error("LrespoGrad : Z vector not converged!")
      end if

      bkm1 = tmp2 / tmp1

      !call scal(bkm1, pkm1)
      !call axpy(1.0_dp, rkm1, pkm1)
      pkm1 = rkm1 + bkm1 * pkm1

    end do

    rhs(:) = rhs2(:)

  end subroutine solveZVectorEq

  ! Calculate Z-dependent parts of the W-vectors and divide diagonal
  ! elements of W_ij and W_ab by 2.
  subroutine calcWvectorZ(zz, win, homo, nocc, nmatup, getij, iAtomStart, &
      & stimc, c, gammaMat, ev, wov, woo, wvv)
    real(dp), intent(in) :: zz(:)
    integer, intent(in) :: win(:), homo, nocc, nmatup, getij(:,:), iAtomStart(:)
    real(dp), intent(in) :: stimc(:,:,:), c(:,:,:), gammaMat(:,:), ev(:)
    real(dp), intent(inout) :: wov(:), woo(:), wvv(:)

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
      wov(ia) = wov(ia) + zz(ia) * ev(i)
    end do

    ! Missing sum_kb 4 K_ijkb Z_kb term in W_ij:
    ! zq(iAt1) = sum_kb q^kb(iAt1) Z_kb
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
        woo(ij) = woo(ij)&
            & + 4.0_dp * qij(iAt1) * gamxpyq(iAt1)
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
  !! \param iAtomStart indexing array for atoms
  !! \param pc density matrix
  !! \param s overlap matrix
  !! \param dqex output atomic charges
  !! \note assumes both triangles of both square matrices are filled
  subroutine getExcMulliken(iAtomStart, pc, s, dqex)
    integer, intent(in) :: iAtomStart(:)
    real(dp), intent(in) :: pc(:,:), s(:,:)
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
  !! \param sym symmetry label
  !! \param nstat state index
  !! \param dq ground state net charge
  !! \param dqex change in atomic charges from ground to excited state
  !! \param coord0 central cell coordinates
  !! \param fdMulliken file unit for Mulliken data
  !! \param fdTagged unit for autotest tagging data
  subroutine writeExcMulliken(sym, nstat, dq, dqex, coord0, fdMulliken)
    character, intent(in) :: sym
    integer, intent(in)   :: nstat
    real(dp), intent(in)  :: dq(:), dqex(:), coord0(:,:)
    integer, intent(in) :: fdMulliken

    integer :: natom, m
    real(dp) :: dipol(3), dipabs

    natom = size(dq)

    @:ASSERT(size(dq) == size(dqex))
    @:ASSERT(all(shape(coord0) == [3,nAtom]))

    ! Output of excited state Mulliken charges
    open(fdMulliken, file=excitedQOut,position="rewind", status="replace")
    write(fdMulliken, "(a,a,i2)") "# MULLIKEN CHARGES of excited state ",&
        & sym, nstat
    write(fdMulliken, "(a,2x,A,i4)") "#", 'Natoms =',natom
    write(fdMulliken, "(a)") "# Atom     netCharge  "
    write(fdMulliken, "(a)") "#============================================"
    do m = 1,  natom
      write(fdMulliken,"(i5,5x,f10.6)") m, -dq(m) - dqex(m)
    end do
    close(fdMulliken)

    ! Calculation of excited state dipole moment
    dipol(:) = -1.0_dp * matmul(coord0, dq + dqex)
    dipabs = sqrt(sum(dipol**2))

    open(fdMulliken, file=excitedDipoleOut, position="rewind", status="replace")
    write(fdMulliken, "(a,a,i2)") "Mulliken analysis of excited state ",&
        & sym, nstat
    write(fdMulliken, "(a)") "=============================================="
    write(fdMulliken, "(a)") " "
    write(fdMulliken, "(a)") "Mulliken exc. state dipole moment [Debye]"
    write(fdMulliken, "(a)") "=============================================="
    write(fdMulliken, "(3(1x,f20.12))") (dipol(m) * au__Debye, m = 1, 3)
    write(fdMulliken, "(a)") " "
    write(fdMulliken, "(a)") "Norm of exc. state dipole moment [Debye]"
    write(fdMulliken, "(a)") "=============================================="
    write(fdMulliken, "(1x,f20.12)") dipabs * au__Debye
    write(fdMulliken, "(a)")
    close(fdMulliken)

  end subroutine writeExcMulliken

  !> Calculate transition moments for transitions between Kohn-Sham
  !> states, including spin-flipping transitions
  !! \param coord0 Atomic positions
  !! \param win transition energies
  !! \param nmatup number of same-spin transitions
  !! \param getij index array
  !! \param iAtomStart index array for ground state square matrices
  !! \param stimc overlap times ground state wavefunctions
  !! \param grndEigVecs ground state wavefunctions
  !! \param snglPartTransDip resulting transition dipoles
  subroutine calcTransitionDipoles(coord0, win, nmatup, getij,&
      & iAtomStart, stimc, grndEigVecs, snglPartTransDip)
    real(dp), intent(in) :: coord0(:,:)
    integer, intent(in) :: win(:), nmatup, iAtomStart(:), getij(:,:)
    real(dp), intent(in) :: stimc(:,:,:), grndEigVecs(:,:,:)
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

  !> Calculation of force from derivative of excitation energy
  !! 1. we need the ground and excited Mulliken charges
  !! 2. we need P,(T,Z),W, X + Y from linear response
  !! 3. calculate dsmndr, dhmndr (dS/dR, dh/dR),
  !! dgabda (dGamma_{IAt1,IAt2}/dR_{IAt1}),
  !! dgext (dGamma-EXT_{IAt1,k}/dR_{IAt1})
  !! \param sym symmetry of the transition
  !! \param nxov number of single particle transitions to include
  !! \param natom number of central cell atoms
  !! \param species0 central cell chemical species
  !! \param iAtomStart index array for S and H0 ground state square matrices
  !! \param norb number of orbitals for ground state system
  !! \param homo number of highest occupied state in ground state
  !! \param nocc number of occupied states in calculation (not
  !! neccessarily same as HOMO in the case of windowing)
  !! \param nmatup number of up->up transitions
  !! \param getij index array from composite transition index to specific
  !! single particle states
  !! \param win single particle excitation energies
  !! \param grndEigVecs ground state eigenvectors
  !! \param pc transition density matrix
  !! \param stimc overlap times ground state eigenvectors
  !! \param dq ground state net charges
  !! \param dqex charge differences from ground to excited state
  !! \param gammaMat softened coulomb matrix
  !! \param HubbardU ground state Hubbard U values
  !! \param spinW ground state spin derivatives for each species
  !! \param shift ground state potentials (shift vector)
  !! \param woo W vector occupied part
  !! \param wov W vector occupied-virtual part
  !! \param wvv W vector virtual part
  !! \param xpy X+Y Furche term
  !! \param coord0 central cell atomic coordinates
  !! \param orb data type for atomic orbital information
  !! \param skHamCont H0 data
  !! \param skOverCont overlap data
  !! \param derivator Differentiatior for the non-scc matrices
  !! \param rhoSqr ground state density matrix for spin-free case
  !! \param excgrad resulting excited state gradient
  subroutine addGradients(sym, nxov, natom, species0, iAtomStart, norb, homo, &
      & nocc, nmatup, getij, win, grndEigVecs, pc, stimc, dq, dqex, gammaMat, &
      & HubbardU, spinW, shift, woo, wov, wvv, xpy, coord0, orb, &
      & skHamCont, skOverCont, derivator, rhoSqr, excgrad)
    character, intent(in) :: sym
    integer, intent(in)   :: nxov, natom, species0(:), iAtomStart(:)
    integer, intent(in)   :: norb, homo, nocc, win(:)
    integer, intent(in)   :: nmatup, getij(:,:)
    real(dp), intent(in)  :: grndEigVecs(:,:,:), pc(:,:), stimc(:,:,:)
    real(dp), intent(in)  :: dq(:), dqex(:)
    real(dp), intent(in)  :: gammaMat(:,:)
    real(dp), intent(in)  :: HubbardU(:), spinW(:), shift(:)
    real(dp), intent(in)  :: woo(:), wov(:), wvv(:), xpy(:)
    real(dp), intent(in)  :: coord0(:,:)
    type(TOrbitals), intent(in)   :: orb
    type(OSlakoCont), intent(in)  :: skHamCont, skOverCont
    class(NonSccDiff), intent(in) :: derivator
    real(dp), intent(in)  :: rhoSqr(:,:)
    real(dp), intent(out) :: excgrad(:,:)

    real(dp), allocatable :: shift_excited(:), xpyq(:)
    real(dp), allocatable :: shxpyq(:), xpycc(:,:), wcc(:,:)
    real(dp), allocatable :: qij(:), temp(:)
    real(dp), allocatable :: dH0(:,:,:), dS(:,:,:)
    integer :: ia, i, j, a, b, ab, ij, m, n, mu, nu, xyz, iAt1, iAt2
    integer :: indalpha, indalpha1, indbeta, indbeta1
    integer :: nOrb1, nOrb2, iSp1, iSp2
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
    ! (xpycc)_{mu nu} =
    ! =  sum_{ia} (X + Y)_{ia} (grndEigVecs(mu,i)grndEigVecs(nu,a)
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
              &( grndEigVecs(mu,i,1)*grndEigVecs(nu,a,1) &
              & + grndEigVecs(mu,a,1)*grndEigVecs(nu,i,1) )
        end do
      end do
    end do

    ! calculate wcc = c_mu,i * W_ij * c_j,nu
    ! We have only W_ab b > a and W_ij j > i:
    ! wcc(m,n) = sum_{pq, p <= q} w_pq (grndEigVecs(mu,p)grndEigVecs(nu,q)
    ! + grndEigVecs(nu,p)grndEigVecs(mu,q))
    ! complexity norb * norb * norb

    ! calculate the occ-occ part
    ! BA: Does not the diagonal contain twice as much as needed?
    wcc(:,:) = 0.0_dp

    do ij = 1, nxoo
      call indxoo(homo, nocc, ij, i, j)
      do mu = 1, norb
        do nu = 1, norb
          wcc(mu,nu) = wcc(mu,nu) + woo(ij) * &
              & ( grndEigVecs(mu,i,1)*grndEigVecs(nu,j,1) &
              & + grndEigVecs(mu,j,1)*grndEigVecs(nu,i,1) )
        end do
      end do

    end do

    ! calculate the occ-virt part : the same way as for xpycc

    do ia = 1, nxov
      call indxov(win, ia, getij, i, a)
      do nu = 1, norb
        do mu = 1, norb
          wcc(mu,nu) = wcc(mu,nu) + wov(ia) * &
              & ( grndEigVecs(mu,i,1)*grndEigVecs(nu,a,1) &
              & + grndEigVecs(mu,a,1)*grndEigVecs(nu,i,1) )
        end do
      end do
    end do

    ! calculate the virt - virt part
    do ab =1, nxvv
      call indxvv(homo, ab, a, b)
      do mu = 1, norb
        do nu = 1, norb
          wcc(mu,nu) = wcc(mu,nu) + wvv(ab) * &
              &( grndEigVecs(mu,a,1)*grndEigVecs(nu,b,1) &
              &+ grndEigVecs(mu,b,1)*grndEigVecs(nu,a,1) )
        end do
      end do
    end do

    !! now calculating the force !
    !! complexity : norb * norb * 3

    ! as have already performed norb**3 operation to get here,
    ! calculate for all atoms

    ! BA: only for non-periodic systems!
    do iAt1 = 1, nAtom
      indalpha = iAtomStart(iAt1)
      nOrb1 = iAtomStart(iAt1+1) - indalpha
      indalpha1 = iAtomStart(iAt1 + 1) -1
      iSp1 = species0(iAt1)

      do iAt2 = 1, iAt1 - 1
        indbeta = iAtomStart(iAt2)
        nOrb2 = iAtomStart(iAt2+1) - indbeta
        indbeta1 = iAtomStart(iAt2 + 1) -1
        iSp2 = species0(iAt2)

        diffvec = coord0(:,iAt1) - coord0(:,iAt2)
        rab = sqrt(sum(diffvec**2))

        diffvec = diffvec / rab ! now holds unit vector in direction

        ! calculate the derivative of gamma
        dgab(:) = diffvec(:) * &
            & (-1.0_dp / rab**2 - expGammaPrime(rab, HubbardU(iSp1), HubbardU(iSp2)))

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
          excgrad(xyz,iAt1) = excgrad(xyz,iAt1) &
              & + tmp1 + tmp2 + tmp4 + tmp6 + tmp3
          excgrad(xyz,iAt2) = excgrad(xyz,iAt2) &
              & - tmp1 - tmp2 - tmp4 - tmp6 - tmp3
        end do
      end do
    end do

  end subroutine addGradients

  !> Write out excitations projected onto ground state
  !! \param tt density matrix in the MO basis
  !! \param grndEigVecs ground state eigenvectors
  !! \param occ ground state occupations
  !! \param nocc number of filled states
  !! \param fdCoeffs file descriptor to write data into
  subroutine writeCoeffs(tt, grndEigVecs, occ, nocc, fdCoeffs, tCoeffs, &
      & tIncGroundState, occNatural, naturalOrbs)
    real(dp), intent(in) :: tt(:,:), grndEigVecs(:,:,:)
    real(dp), intent(in) :: occ(:,:)
    integer, intent(in)  :: nocc, fdCoeffs
    logical, intent(in)  :: tCoeffs
    logical, intent(in)  :: tIncGroundState
    real(dp), intent(out), optional :: occNatural(:), naturalOrbs(:,:)

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
        naturalOrbs = t2
        call evalCoeffs(naturalOrbs,occNatural,grndEigVecs(:,:,1))
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
        open(fdCoeffs, file=excitedCoefsOut, position="rewind", &
            & status="replace")
        write(fdCoeffs,*) 'T F'
        do ii = 1, norb
          jj = norb - ii + 1
          write(fdCoeffs, '(1x,i3,1x,f13.10,1x,f13.10)') ii, occtmp(jj), 2.0_dp
          write(fdCoeffs, '(6(f13.10,1x))') (cmplx(t2(mm,jj), kind=dp), &
              & mm = 1, norb)
        end do
        close(fdCoeffs)
      end if

    end if

  end subroutine writeCoeffs

  !> Project MO density matrix onto ground state orbitals
  !! \param t2 density matrix
  !! \param occ resulting natural orbital occupations
  !! \param eig 'natural' eigenvectors
  subroutine evalCoeffs(t2,occ,eig)
    real(dp), intent(inout) :: t2(:,:)
    real(dp), intent(out)   :: occ(:)
    real(dp), intent(in)    :: eig(:,:)

    real(dp), allocatable :: coeffs(:,:)

    ALLOCATE(coeffs(size(occ),size(occ)))

    call heev(t2, occ, 'U', 'V')
    call gemm(coeffs, eig, t2)
    t2 = coeffs

  end subroutine evalCoeffs

  !> Write out transitions from ground to excited state along with
  !> single particle transitions and dipole strengths
  !! \param sym Symmetry label for the type of transition
  !! \param osz oscillator strengths for transitions from ground to
  !! excited states
  !! \param nexc number of excited states to solve for
  !! \param nmatup number of same spin excitations
  !! \param getij index array between transitions in square and 1D
  !! representations
  !! \param win index array for single particle excitions
  !! \param eval excitation energies
  !! \param evec eigenvectors of excited states
  !! \param wij single particle excitation energies
  !! \param fdXplusY file unit for X+Y data
  !! \param fdTrans file unit for transitions
  !! \param fdTradip file unit for transition dipoles
  !! \param transitionDipoles
  !! \param tWriteTagged print tag information
  !! \param fdTagged file unit for tagged output (> -1 for write out)
  !! \param fdExc file unit for excitation energies
  !! \param Ssq For spin polarized systems, measure of spin
  subroutine writeExcitations(sym, osz, nexc, nmatup, getij, win, eval, evec, &
      & wij, fdXplusY, fdTrans, fdTradip, transitionDipoles, tWriteTagged, &
      & fdTagged, fdExc, Ssq)
    character, intent(in) :: sym
    real(dp), intent(in) :: osz(:)
    integer, intent(in) :: nexc, nmatup, getij(:,:), win(:)
    real(dp), intent(in) :: eval(:), evec(:,:), wij(:), transitionDipoles(:,:)
    logical, intent(in) :: tWriteTagged
    integer, intent(in) :: fdTradip, fdXplusY, fdTrans, fdTagged, fdExc
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
          write(fdExc, &
              & '(1x,f10.3,4x,f12.7,4x,i5,3x,a,1x,i5,7x,f6.3,2x,f10.3,4x, &
              & f6.3)') &
              & Hartree__eV * sqrt(eval(i)), osz(i), m, '->', n, weight, &
              & Hartree__eV * wij(iweight), Ssq(i)
        else
          write(fdExc, &
              & '(1x,f10.3,4x,f12.7,4x,i5,3x,a,1x,i5,7x,f6.3,2x,f10.3,6x,a)') &
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
          write(fdTrans,*)
          write(fdTrans, '(2x,a,1x,i5,5x,f10.3,1x,a,3x,a)') &
              & 'Energy ', i,  Hartree__eV * sqrt(eval(i)), 'eV', sign
          write(fdTrans,*)
          write(fdTrans,'(2x,a,9x,a,7x,a)') &
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
                & '(i5,3x,a,1x,i5,1x,1a,1x,f8.4,2x,f10.3)') &
                & m, '->', n, sign, wvec(j), Hartree__eV * wij(wvin(j))
          end do
        end if

        if(fdTradip > 0) then
          write(fdTradip, '(1x,i5,1x,f10.3,2x,3(3x,f10.6))') &
              & i, Hartree__eV * sqrt(eval(i)), (transitionDipoles(i,j) &
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
          write(fdExc, &
              & '(6x,A,T12,4x,f12.7,4x,i5,3x,a,1x,i5,7x,A,2x,f10.3,4x,f6.3)') &
              & '< 0', osz(i), m, '->', n, '-', Hartree__eV * wij(iweight), &
              & Ssq(i)
        else
          write(fdExc, &
              & '(6x,A,T12,4x,f12.7,4x,i5,3x,a,1x,i5,7x,f6.3,2x,f10.3,6x,a)') &
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
          write(fdTrans,*)
          write(fdTrans, '(2x,a,1x,i5,5x,a,1x,a,3x,a)') &
              & 'Energy ', i,  '-', 'eV', sign
          write(fdTrans,*)
        end if

        if(fdTradip > 0) then
          write(fdTradip, '(1x,i5,1x,A)') i, '-'
        endif

      end if

    end do

    ! Determine degenerate levels and sum oscillator strength over any
    ! degenerate levels
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

  !> Create transition density matrix in MO basis
  !! P = T + 1/2 Z symmetric (paper has T + Z asymmetric)
  !! (Zab = Zij = 0, Tia = 0)
  !! \param t
  !! \param rhs
  !! \param win index array for single particle transitions
  !! \param getij
  !! \param pc resulting excited state density matrix
  subroutine calcPMatrix(t, rhs, win, getij, pc)
    real(dp), intent(in) :: t(:,:), rhs(:)
    integer, intent(in) :: win(:), getij(:,:)
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

  !> Write single particle excitations to a file as well as
  !! potentially to tagged output file (in that case, summing over
  !! degeneracies)
  !! \param wij single particle excitation energies
  !! \param win index array
  !! \param nmatup number of transitions within same spin channel
  !! \param getij index from composite index to occupied and virtual
  !! single particle states
  !! \param fdTagged file descriptor for the tagged data output
  !! \param fgTagged tagged output file
  !! \param fdSPTrans file for transitions
  !! \param sposz single particle oscilation strengths
  !! \param nxov Number of included single particle excitations to
  !! print out (assumes that win and wij are sorted so that these are
  !! first
  !! \param tSpin is this a spin-polarized calculation?
  subroutine writeSPExcitations(wij, win, nmatup, getij, fdSPTrans, sposz, nxov, tSpin)
    real(dp), intent(in)  :: wij(:)
    integer, intent(in)   :: win(:), nmatup, getij(:,:)
    integer, intent(in)   :: fdSPTrans
    real(dp), intent(in)  :: sposz(:)
    integer, intent(in)   :: nxov
    logical, intent(in)   :: tSpin

    integer :: indm, m, n
    logical :: updwn
    character :: sign

    @:ASSERT(size(sposz)>=nxov)


    if (fdSPTrans > 0) then
      ! single particle excitations
      open(fdSPTrans, file=singlePartOut, position="rewind", status="replace")
      write(fdSPTrans,*)
      write(fdSPTrans,'(7x,a,7x,a,8x,a)') '#      w [eV]', &
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
        write(fdSPTrans, &
            & '(1x,i7,3x,f8.3,3x,f13.7,4x,i5,3x,a,1x,i5,1x,1a)') &
            & indm, Hartree__eV * wij(indm), sposz(indm), m, '->', n, sign
      end do
      write(fdSPTrans,*)
      close(fdSPTrans)
    end if

  end subroutine writeSPExcitations

end module linrespgrad
