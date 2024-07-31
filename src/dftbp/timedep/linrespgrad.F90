!--------------------------------------------------------------------------------------------------!
!  DFTB+: general package for performing fast atomistic simulations                                !
!  Copyright (C) 2006 - 2023  DFTB+ developers group                                               !
!                                                                                                  !
!  See the LICENSE file for terms of usage and distribution.                                       !
!--------------------------------------------------------------------------------------------------!

#:include 'common.fypp'

!> Linear response excitations and gradients with respect to atomic coordinates.
!!
!! Note: This module is NOT instance safe it uses a common block to communicate with ARPACK
module dftbp_timedep_linrespgrad
  use dftbp_common_accuracy, only : dp, elecTolMax, lc, rsp
  use dftbp_common_constants, only : Hartree__eV, au__Debye, cExchange
  use dftbp_common_schedule, only : distributeRangeInChunks, assembleChunks
  use dftbp_io_commonformats, only : format2U
  use dftbp_common_globalenv, only : stdOut
  use dftbp_common_file, only : TFileDescr, openFile, closeFile, clearFile
  use dftbp_dftb_nonscc, only : TNonSccDiff
  use dftbp_dftb_hybridxc, only : THybridXcFunc, getDirectionalCamGammaPrimeValue
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
      & getExcSpin, writeExcMulliken, actionAplusB, actionAminusB, initialSubSpaceMatrixApmB,&
      & calcMatrixSqrt, incMemStratmann, orthonormalizeVectors, getSqrOcc
  use dftbp_timedep_linresptypes, only : TLinResp, linrespSolverTypes
  use dftbp_timedep_transcharges, only : TTransCharges, transq, TTransCharges_init
  use dftbp_type_commontypes, only : TOrbitals
  use dftbp_type_densedescr, only : TDenseDescr
  use dftbp_common_environment, only : TEnvironment, globalTimers
  
#:if WITH_SCALAPACK

  use dftbp_timedep_linrespcommon, only : getExcSpin_MPI
  use dftbp_extlibs_scalapackfx, only : pblasfx_psymm
  use dftbp_extlibs_mpifx, only : MPI_SUM, mpifx_allreduceip
  use dftbp_math_scalafxext, only : distrib2replicated

#:endif
  
  implicit none

  private
  public :: LinRespGrad_old, conicalIntersectionOptimizer

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
  real(dp), parameter :: CONV_THRESH_STRAT = 0.1*epsilon(1.0_rsp)

  !> Maximal allowed iteration in the ARPACK solver.
  integer, parameter :: MAX_AR_ITER = 300

  !> Names of output files
  character(*), parameter :: arpackOut = "ARPACK.DAT"
  character(*), parameter :: testArpackOut = "TEST_ARPACK.DAT"

  !> Treshold for near-identical NACV values to fix phase
  real(dp), parameter :: nacTol = 1.0e-6_dp

contains

  !> This subroutine analytically calculates excitations and gradients of excited state energies
  !> based on Time Dependent DFRT
  subroutine LinRespGrad_old(env, this, denseDesc, grndEigVecs, grndEigVal, sccCalc, dq, coord0,&
      & SSqr, filling, species0, iNeighbour, img2CentCell, orb, fdTagged, taggedWriter, hybridXc,&
      & omega, allOmega, deltaRho, shift, skHamCont, skOverCont, excgrad, nacv, derivator, rhoSqr,&
      & occNatural, naturalOrbs)

    !> Environment settings
    type(TEnvironment), intent(inout) :: env
    
    type(TLinResp), intent(inout) :: this

    !> Dense matrix descriptor
    type(TDenseDescr), intent(in) :: denseDesc 

    !> Ground state MO-coefficients
    real(dp), intent(in) :: grndEigVecs(:,:,:)

    !> Ground state MO-energies
    real(dp), intent(in) :: grndEigVal(:,:)

    !> Self-consistent charge module settings
    type(TScc), intent(in) :: sccCalc

    !> Converged ground state Mulliken gross charges - atomic charges
    real(dp), intent(in) :: dq(:,:)

    !> Atomic positions
    real(dp), intent(in) :: coord0(:,:)

    !> Square overlap matrix between basis functions, both triangles required
    real(dp), intent(in) :: SSqr(:,:)

    !> Occupations for the states
    real(dp), intent(in) :: filling(:,:)

    !> Chemical species of each atom
    integer, intent(in) :: species0(:)

    !> Atomic neighbour lists
    integer, intent(in) :: iNeighbour(0:,:)

    !> Mapping of atom number to central cell atom number
    integer, intent(in) :: img2CentCell(:)

    !> Data type for atomic orbital information
    type(TOrbitals), intent(in) :: orb

    !> File descriptor for the tagged data output
    type(TFileDescr), intent(in) :: fdTagged

    !> Tagged writer
    type(TTaggedWriter), intent(inout) :: taggedWriter

    !> Data for hybrid xc-functional calculation
    class(THybridXcFunc), allocatable, intent(inout) :: hybridXc

    !> Excitation energy of state nStat
    real(dp), intent(out) :: omega

    !> Excitation energy of all states that have been solved
    real(dp), allocatable, intent(inout) :: allOmega(:)

    !> Difference density matrix (vs. uncharged atoms)
    real(dp), intent(inout), optional :: deltaRho(:,:,:)

    !> Shift vector for potentials in the ground state
    real(dp), intent(in), optional :: shift(:)

    !> Non-SCC hamiltonian data
    type(TSlakoCont), intent(in), optional :: skHamCont

    !> Overlap data
    type(TSlakoCont), intent(in), optional :: skOverCont

    !> Excitation energy gradients with respect to atomic positions
    real(dp), intent(out), optional :: excgrad(:,:,:)

    !> Non-adiabatic coupling vectors
    real(dp), intent(out), optional :: nacv(:,:,:)

    !> Differentiator for H0 and S matrices
    class(TNonSccDiff), intent(in), optional :: derivator

    !> Ground state density matrix
    real(dp), intent(in), optional :: rhoSqr(:,:,:)

    !> Occupation numbers for natural orbitals from the excited state density matrix
    real(dp), intent(out), optional :: occNatural(:)

    !> The single particle eigenvectors themselves for the excited state density matrix
    real(dp), intent(out), optional :: naturalOrbs(:,:,:)


    real(dp) :: Ssq(this%nExc), omegaDif, omegaAvg
    real(dp), allocatable :: gammaMat(:,:), lrGamma(:,:), snglPartTransDip(:,:)
    real(dp), allocatable :: ovrXev(:,:,:), wij(:)
    real(dp), allocatable :: dqex(:,:), sposz(:), osz(:), pc(:,:,:)
    real(dp), allocatable :: xpy(:,:), xmy(:,:), sqrOccIA(:)
    real(dp), allocatable :: xpym(:), xpyn(:), xmyn(:), xmym(:)
    real(dp), allocatable :: t(:,:,:), rhs(:), woo(:,:), wvv(:,:), wov(:)
    real(dp), allocatable :: eval(:), transitionDipoles(:,:)
    real(dp), allocatable :: eigVecGlb(:,:,:), ovrXevGlb(:,:,:)
    integer, allocatable :: win(:), getIA(:,:), getIJ(:,:), getAB(:,:)

    !> Array from pairs of single particles states to compound index - should replace with a more
    !> compact data structure in the cases where there are oscilator windows
    integer, allocatable :: iatrans(:,:,:)

    character, allocatable :: symmetries(:)

    integer :: nxoo_max, nxvv_max
    integer, allocatable :: nocc_ud(:), nvir_ud(:)
    integer :: mHOMO, mLUMO
    integer :: nxov, nxov_ud(2), nxov_r, nxov_d, nxov_rd, nxoo_ud(2), nxvv_ud(2)
    integer :: norb, nxoo, nxvv
    integer :: i, j, ss, iSpin, isym, iLev, iSav, nStartLev, nEndLev
    integer :: nCoupLev, mCoupLev, iNac, iGlobal, fGlobal
    integer :: nSpin
    character :: sym
    character(lc) :: tmpStr

    real(dp) :: energyThreshold

    integer :: nStat

    !> Control variables
    logical :: tZVector, doAllZVectors, tFracOcc, doVanillaZvector
    logical :: tHybridXc = .false.

    !> Should gradients be calculated
    logical :: tForces

    !> Transition charges, either cached or evaluated on demand
    type(TTransCharges) :: transChrg

    type(TFileDescr) :: fdTrans, fdTransDip, fdArnoldi, fdXPlusY, fdExc

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

    call env%globalTimer%startTimer(globalTimers%lrSetup)

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
      if (this%tArnoldi) then
        msaupd = 1
        msaup2 = 1
      else
        msaupd = 0
        msaup2 = 0
      endif
      ! End of ARPACK communication variables

    end if

    if (this%writeMulliken) then
      call clearFile(excitedQOut)
      call clearFile(excitedDipoleOut)
    end if

    if (this%tArnoldi) then
      call openFile(fdArnoldi, arpackOut, mode="w")
    end if

    nstat = this%nstat
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
    if(allocated(hybridXc)) then
       tHybridXc = .true.
       call env%globalTimer%startTimer(globalTimers%lrCoulomb)
       allocate(lrGamma(this%nAtom, this%nAtom))
       call hybridXc%getCamGammaCluster(lrGamma)
       call env%globalTimer%stopTimer(globalTimers%lrCoulomb)
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
    if (tFracOcc .and. tHybridXc) then
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
    allocate(nocc_ud(nSpin))
    allocate(nvir_ud(nSpin))
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
    end do
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

    !> Is a z vector required?
    tZVector = tForces .or. this%writeMulliken .or. this%writeCoeffs .or. present(naturalOrbs) .or.&
        & this%tWriteDensityMatrix .or. this%tNaCoupling
    doAllZVectors = tZVector .and. (nstat == 0) .and. (.not. this%isCIopt) .and. &
        & (.not. this%tNaCoupling)

    !> Occ-occ/vir-vir charges only required for Z-vector/forces or TD-LC-DFTB
    if (this%tCacheChargesOccVir) then
      if (tZVector .or. tHybridXc) this%tCacheChargesSame = .true.
    end if

    ! Sanity checks
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
    if (.not. this%tSpin) then
      select case (this%symmetry)
      case ("B")
        allocate(symmetries(2))
        symmetries(:) = [ "T", "S" ]
      case ("S")
        allocate(symmetries(1))
        symmetries(:) = [ "S" ]
      case ("T")
        allocate(symmetries(1))
        symmetries(:) = [ "T" ]
      end select
    else
      ! ADG: temporary solution for spin polarized case.
      allocate(symmetries(1))
      symmetries(:) = [ " " ]
    end if
    ! Allocation for general arrays
    allocate(gammaMat(this%nAtom, this%nAtom))
    allocate(snglPartTransDip(nxov, 3))
    allocate(ovrXev(size(grndEigVecs,dim=1), size(grndEigVecs,dim=2), nSpin))
    allocate(wij(nxov))
    allocate(win(nxov))
    allocate(sqrOccIA(nxov))
    allocate(eval(this%nExc))
    allocate(getIA(nxov, 3))
    allocate(getIJ(nxoo, 3))
    allocate(getAB(nxvv, 3))
    allocate(transitionDipoles(this%nExc, 3))
    allocate(sposz(nxov))

    ! Overlap times wave function coefficients - most routines in DFTB+ use lower triangle (would
    ! remove the need to symmetrize the overlap and ground state density matrix in the main code if
    ! this could be used everywhere in these routines)
  #:if WITH_SCALAPACK

    do iSpin = 1, nSpin
      call pblasfx_psymm(SSqr, denseDesc%blacsOrbSqr, grndEigVecs(:,:,iSpin), denseDesc%blacsOrbSqr,&
           & ovrXev(:,:,iSpin), denseDesc%blacsOrbSqr, side="L")
    end do

    ! The routines getExcSpin and the onsite corrections can not easily be parallelized. Global
    ! versions of the BLACS distributed arrays are required.
    if (this%tSpin .or. allocated(this%onSiteMatrixElements)) then
      allocate(eigVecGlb(norb,norb,nSpin))
      allocate(ovrXevGlb(norb,norb,nSpin))
      do ss = 1, nSpin
        call distrib2replicated(env%blacs%orbitalGrid, denseDesc%blacsOrbSqr, grndEigVecs(:,:,ss),&
            & eigVecGlb(:,:,ss))
        call distrib2replicated(env%blacs%orbitalGrid, denseDesc%blacsOrbSqr, ovrXev(:,:,ss),&
            & ovrXevGlb(:,:,ss))
      end do
    end if

    call env%globalTimer%startTimer(globalTimers%lrCoulomb)
    call sccCalc%getAtomicGammaMatrixBlacs(gammaMat, iNeighbour, img2CentCell, env)
    call env%globalTimer%stopTimer(globalTimers%lrCoulomb)
    
  #:else
    
    do iSpin = 1, nSpin
      call symm(ovrXev(:,:,iSpin), "L", SSqr, grndEigVecs(:,:,iSpin))
    end do

    if (this%tSpin .or. allocated(this%onSiteMatrixElements)) then
      allocate(eigVecGlb(norb,norb,nSpin))
      allocate(ovrXevGlb(norb,norb,nSpin))
      eigVecGlb = grndEigVecs
      ovrXevGlb = ovrXev
    end if
    
    call env%globalTimer%startTimer(globalTimers%lrCoulomb)
    call sccCalc%getAtomicGammaMatrix(gammaMat, iNeighbour, img2CentCell)
    call env%globalTimer%stopTimer(globalTimers%lrCoulomb)
   
  #:endif

    ! Oscillator strengths for exited states, when needed.
    allocate(osz(this%nExc))

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
    call getSqrOcc(filling, win, nxov_ud(1), nxov, getIA, this%tSpin, sqrOccIA)

    call env%globalTimer%startTimer(globalTimers%lrTransCharges)

    ! First call to initialize charges for all occ-vir transitions
    call TTransCharges_init(transChrg, env, denseDesc, ovrXev, grndEigVecs, norb, nxov,&
        & nxov_ud(1), nxoo_ud, nxvv_ud, getIA, getIJ, getAB, win, this%tCacheChargesOccVir,&
        & this%tCacheChargesSame, .true.)

    ! dipole strength of transitions between K-S states
    call calcTransitionDipoles(coord0, win, norb, nSpin, nxov_ud(1), getIA, transChrg, env, &
        & denseDesc, ovrXev, grndEigVecs, snglPartTransDip)

    call env%globalTimer%stopTimer(globalTimers%lrTransCharges)

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
    
    ! Recompute occ-vir transition charges, since win/wij and number has changed
    if (nxov_rd /= nxov .or. this%tOscillatorWindow .or. this%tEnergyWindow) then
      call env%globalTimer%startTimer(globalTimers%lrTransCharges)
      call TTransCharges_init(transChrg, env, denseDesc, ovrXev, grndEigVecs, norb, nxov_rd,&
        & nxov_ud(1), nxoo_ud, nxvv_ud, getIA, getIJ, getAB, win, this%tCacheChargesOccVir,&
        & this%tCacheChargesSame, .false.)
      call env%globalTimer%stopTimer(globalTimers%lrTransCharges)
    end if

    ! MPI distribution of RPA vectors according to these indices 
    call distributeRangeInChunks(env, 1, nxov_rd, iGlobal, fGlobal)

    if (this%writeXplusY) then
      call openfile(fdXPlusY, XplusYOut, mode="w")
    end if

    if (this%writeTrans) then
      call openFile(fdTrans, transitionsOut, mode="w")
      write(fdTrans%unit,*)
    end if

    ! Many-body transition dipole file to excited states
    if (this%writeTransDip) then
      call openFile(fdTransDip, transDipOut, mode="w")
      write(fdTransDip%unit, *)
      write(fdTransDip%unit, '(5x,a,5x,a,2x,a)') "#", 'w [eV]', 'Transition dipole (x,y,z) [Debye]'
      write(fdTransDip%unit, *)
      write(fdTransDip%unit, '(1x,60("="))')
      write(fdTransDip%unit, *)
    end if

    ! excitation energies
    call openFile(fdExc, excitationsOut, mode="w")
    write(fdExc%unit, *)
    if (this%tSpin) then
      write(fdExc%unit, '(5x,a,7x,a,9x,a,9x,a,6x,a,4x,a)') 'w [eV]', 'Osc.Str.', 'Transition',&
          & 'Weight', 'KS [eV]','D<S*S>'
    else
      write(fdExc%unit, '(5x,a,7x,a,9x,a,9x,a,6x,a,4x,a)') 'w [eV]','Osc.Str.', 'Transition',&
          & 'Weight', 'KS [eV]','Sym.'
    end if

    write(fdExc%unit, *)
    write(fdExc%unit, '(1x,80("="))')
    write(fdExc%unit, *)

    ! single particle excitations (output file and tagged file if needed).  Was used for nxov_rd =
    ! size(wij), but now for just states that are actually included in the excitation calculation.
    call writeSPExcitations(wij, win, nxov_ud(1), getIA, this%writeSPTrans, sposz, nxov_rd,&
        & this%tSpin)

    allocate(xpy(nxov_rd, this%nExc))
    if (tZVector .or. tHybridXc) then
      allocate(xmy(nxov_rd, this%nExc))
    end if

    ! set up transition indexing
    allocate(iatrans(norb, norb, nSpin))
    call rindxov_array(win, nxov, nxoo, nxvv, getIA, getIJ, getAB, iatrans)

    if (this%iLinRespSolver /= linrespSolverTypes%stratmann .and. tHybridXc) then
      call error("Range separation requires the Stratmann solver for excitations")
    end if

    call env%globalTimer%stopTimer(globalTimers%lrSetup)

    do isym = 1, size(symmetries)

      sym = symmetries(isym)
      call env%globalTimer%startTimer(globalTimers%lrSolver)
      select case (this%iLinRespSolver)
      case (linrespSolverTypes%arpack)
        call buildAndDiagExcMatrixArpack(iGlobal, fGlobal, this%tSpin, wij(:nxov_rd), sym, win,&
            & nocc_ud, nvir_ud, nxoo_ud, nxvv_ud, nxov_ud, nxov_rd, iaTrans, getIA, getIJ, getAB,&
            & env, denseDesc, ovrXev, ovrXevGlb, grndEigVecs, eigVecGlb, filling,&
            & sqrOccIA(:nxov_rd), gammaMat, species0, this%spinW, transChrg, this%testArnoldi,&
            & eval, xpy, xmy, this%onSiteMatrixElements, orb, tHybridXc, tZVector)
      case (linrespSolverTypes%stratmann)
        call buildAndDiagExcMatrixStratmann(iGlobal, fGlobal, this%tSpin, &
            & this%subSpaceFactorStratmann, wij(:nxov_rd), sym, win, nocc_ud, nvir_ud, nxoo_ud,&
            & nxvv_ud, nxov_ud, nxov_rd, iaTrans, getIA, getIJ, getAB, env, denseDesc, ovrXev,&
            & ovrXevGlb, grndEigVecs, eigVecGlb, filling, sqrOccIA(:nxov_rd), gammaMat, species0,&
            & this%spinW, transChrg, eval, xpy, xmy, this%onSiteMatrixElements, orb, tHybridXc,&
            & lrGamma, tZVector)
      end select
      call env%globalTimer%stopTimer(globalTimers%lrSolver)

      ! Excitation oscillator strengths for resulting states
      call getOscillatorStrengths(sym, this%tSpin, snglPartTransDip(1:nxov_rd,:), eval, xpy,&
            & sqrOccIA(:nxov_rd), nstat, osz, this%writeTransDip, transitionDipoles)

      if (this%tSpin) then

      #:if WITH_SCALAPACK

        call getExcSpin(Ssq, nxov_ud(1), getIA, win, eval, xpy, filling, ovrXevGlb, eigVecGlb)

      #:else
        
        call getExcSpin(Ssq, nxov_ud(1), getIA, win, eval, xpy, filling, ovrXev, grndEigVecs)

      #:endif
        
        call writeExcitations(sym, osz, this%nExc, nxov_ud(1), getIA, win, eval, xpy,&
            & wij(:nxov_rd), fdXPlusY, fdTrans, fdTransDip,&
            & transitionDipoles,  fdTagged, taggedWriter, fdExc, Ssq)
      else
        call writeExcitations(sym, osz, this%nExc, nxov_ud(1), getIA, win, eval, xpy,&
            & wij(:nxov_rd), fdXPlusY, fdTrans, fdTransDip,&
            & transitionDipoles, fdTagged, taggedWriter, fdExc)
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

    call closeFile(fdArnoldi)
    call closeFile(fdTrans)
    call closeFile(fdXPlusY)
    call closeFile(fdTransDip)
    call closeFile(fdExc)

    ! Remove some un-used memory
    deallocate(snglPartTransDip)
    deallocate(transitionDipoles)
    deallocate(sposz)

    ! Calculate Furche vectors and transition density matrix for various properties
    if (tZVector) then

      ! Differentiates between a standard Z-vector equation for transition densities and forces and
      ! a specific one for non-adiabatic couplings
      doVanillaZvector = .true.
      if (doAllZVectors) then

        call env%globalTimer%startTimer(globalTimers%lrZVector)

        nStartLev = 1
        nEndLev = this%nExc

        if (tForces) then
          call error("Forces currently not available unless a single excited state is specified")
        end if

      else if (this%isCIopt) then

        if (this%indNACouplings(1) == 0) then
          nStartLev = this%indNACouplings(1) + 1
        else
          nStartLev = this%indNACouplings(1)
        end if
        nEndLev = this%indNACouplings(2)

      else

        nStartLev = nstat
        nEndLev = nstat
        if (nstat == 0) then
          doVanillaZvector = .false.
        end if

      end if

      if (this%tSpin) then
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
      allocate(t(norb, norb, nSpin))
      allocate(rhs(nxov_rd))
      allocate(woo(nxoo_max, nSpin))
      allocate(wvv(nxvv_max, nSpin))
      allocate(wov(nxov_rd))

      ! Arrays for gradients and Mulliken analysis
      if (tZVector) then
        allocate(dqex(this%nAtom, nSpin))
        allocate(pc(norb, norb, nSpin))
      end if

      if (doVanillaZvector) then
        
        call env%globalTimer%startTimer(globalTimers%lrGradients)
        do iLev = nStartLev, nEndLev

          omega = sqrt(eval(iLev))

         ! solve for Z and W to get excited state density matrix
         call getZVectorEqRHS(tHybridXc, xpy(:,iLev), xmy(:,iLev), win, env, denseDesc,&
            & nocc_ud, transChrg, getIA, getIJ, getAB, iatrans, this%nAtom, species0, grndEigVal,&
            & ovrXev, grndEigVecs, gammaMat, lrGamma, this%spinW, omega, sym, rhs, t,&
            & wov, woo, wvv)

          call solveZVectorPrecond(rhs, this%tSpin, wij(:nxov_rd), win, nocc_ud, nvir_ud, nxoo_ud,&
            & nxvv_ud, nxov_ud, nxov_rd, iaTrans, getIA, getIJ, getAB, this%nAtom, env, denseDesc,&
            & ovrXev, ovrXevGlb, grndEigVecs, eigVecGlb, filling, sqrOccIA(:nxov_rd), gammaMat,&
            & species0, this%spinW, this%onSiteMatrixElements, orb, transChrg, tHybridXc, lrGamma)

          call calcWVectorZ(rhs, win, nocc_ud, getIA, getIJ, getAB, iaTrans, env, denseDesc,&
            & ovrXev, grndEigVecs, gammaMat, grndEigVal, wov, woo, wvv, transChrg, species0, &
            & this%spinW, tHybridXc, lrGamma)

          call calcPMatrix(t, rhs, win, getIA, pc)

          call writeCoeffs(pc, grndEigVecs, filling, this%writeCoeffs, this%tGrndState, occNatural,&
            & naturalOrbs)

          do iSpin = 1, nSpin
            ! Make MO to AO transformation of the excited density matrix
            call makeSimilarityTrans(pc(:,:,iSpin), grndEigVecs(:,:,iSpin))
            call getExcMulliken(denseDesc, pc(:,:,iSpin), SSqr, dqex(:,iSpin))
          end do

          if (this%tWriteDensityMatrix) then
            call writeDM(iLev, pc, rhoSqr)
          end if

          if (this%writeMulliken) then
            !> For now, only total Mulliken charges
            call writeExcMulliken(sym, iLev, dq(:,1), sum(dqex,dim=2), coord0)
          end if

          if (tForces) then
            iSav = iLev - nStartLev + 1
            call addGradients(sym, nxov_rd, this%nAtom, species0, env, denseDesc, norb, nocc_ud,&
                & getIA, getIJ, getAB, win, grndEigVecs, pc, ovrXev, dq, dqex, gammaMat,&
                & lrGamma, this%HubbardU, this%spinW, shift, woo, wov, wvv, transChrg, xpy(:,iLev),&
                & xmy(:,iLev), coord0, orb, skHamCont, skOverCont, derivator, rhoSqr, tHybridXc,&
                & hybridXc, excgrad(:,:,iSav), deltaRho=deltaRho)
          end if
        end do
        call env%globalTimer%stopTimer(globalTimers%lrGradients)      
      end if

      if (this%tNaCoupling) then

        call env%globalTimer%startTimer(globalTimers%lrNAC)

        ! This overwrites T, RHS and W
        allocate(xpyn, mold=xpy(:,1))
        allocate(xpym, mold=xpy(:,1))
        allocate(xmyn, mold=xpy(:,1))
        allocate(xmym, mold=xpy(:,1))

        iNac = 0
        nacv(:,:,:) = 0.0_dp
        do nCoupLev = this%indNACouplings(1), this%indNACouplings(2)-1
          do mCoupLev = nCoupLev+1, this%indNACouplings(2)

            iNac = iNac + 1
            woo(:,:) = 0.0_dp
            wvv(:,:) = 0.0_dp
            pc(:,:,:) = 0.0_dp
            t(:,:,:) = 0.0_dp
            rhs(:) = 0.0_dp

            !> Ground-to-excited NACV
            if (nCoupLev == 0) then
              xpym(:) = xpy(:,mCoupLev)
              xmym(:) = xmy(:,mCoupLev)
              omegaDif = sqrt(eval(mCoupLev))
              call grndToExcDensityMatrices(tHybridXc, xpym, xmym, win, env, denseDesc, nocc_ud,&
                 & transChrg, getIA, getIJ, getAB, iatrans, this%nAtom, species0, grndEigVal,&
                 & ovrXev, grndEigVecs, gammaMat, lrGamma, this%spinW, omegaDif, sym, pc, wov,&
                 & woo)

              do iSpin = 1, nSpin
                ! Make MO to AO transformation of the excited density matrix
                call makeSimilarityTrans(pc(:,:,iSpin), grndEigVecs(:,:,iSpin))
                call getExcMulliken(denseDesc, pc(:,:,iSpin), SSqr, dqex(:,iSpin))
              end do

              !> For 0-n couplings, the standard force routine can be used, where
              !> X+Y and W_ab are zeroed out
              wvv(:,:)  = 0.0_dp
              xpym(:) = 0.0_dp
              xmym(:) = 0.0_dp
              xpyn(:) = 0.0_dp
              xmyn(:) = 0.0_dp
              call addGradients(sym, nxov_rd, this%nAtom, species0, env, denseDesc, norb,&
                & nocc_ud, getIA, getIJ, getAB, win, grndEigVecs, pc, ovrXev, dq, dqex, gammaMat,&
                & lrGamma, this%HubbardU, this%spinW, shift, woo, wov, wvv, transChrg, xpym,&
                & xmym, coord0, orb, skHamCont, skOverCont, derivator, rhoSqr, tHybridXc, hybridXc,&
                & nacv(:,:,iNac), deltaRho=deltaRho)

            else

              xpyn(:) = xpy(:,nCoupLev)
              xmyn(:) = xmy(:,nCoupLev)
              xpym(:) = xpy(:,mCoupLev)
              xmym(:) = xmy(:,mCoupLev)
              omegaDif = sqrt(eval(nCoupLev)) - sqrt(eval(mCoupLev))
              omegaAvg = 0.5_dp * (sqrt(eval(nCoupLev)) + sqrt(eval(mCoupLev)))

              ! compute + component of RHS for Z-vector eq. in the NaCoupling case
              ! also computes the + components of W and T
              call getNadiaZvectorEqRHS(tHybridXc, xpy(:,nCoupLev), xmy(:,nCoupLev),&
                 & xpy(:,mCoupLev), xmy(:,mCoupLev), win, env, denseDesc, nocc_ud, transChrg,&
                 & getIA, getIJ, getAB, iatrans, this%nAtom, species0, grndEigVal, ovrXev,&
                 & grndEigVecs, gammaMat, lrGamma, this%spinW, omegaAvg, sym, rhs, t,&
                 & wov, woo, wvv)

              call solveZVectorPrecond(rhs, this%tSpin, wij(:nxov_rd), win, nocc_ud, nvir_ud,&
                  & nxoo_ud, nxvv_ud, nxov_ud, nxov_rd, iaTrans, getIA, getIJ, getAB, this%nAtom,&
                  & env, denseDesc, ovrXev, ovrXevGlb, grndEigVecs, eigVecGlb, filling,&
                  & sqrOccIA(:nxov_rd), gammaMat, species0, this%spinW, this%onSiteMatrixElements,&
                  & orb, transChrg, tHybridXc, lrGamma)

              call calcWVectorZ(rhs, win, nocc_ud, getIA, getIJ, getAB, iaTrans,&
                 & env, denseDesc, ovrXev, grndEigVecs, gammaMat, grndEigVal, wov, woo, wvv,&
                 & transChrg, species0, this%spinW, tHybridXc, lrGamma)

              call calcPMatrix(t, rhs, win, getIA, pc)

              do iSpin = 1, nSpin
                ! Make MO to AO transformation of the excited density matrix
                call makeSimilarityTrans(pc(:,:,iSpin), grndEigVecs(:,:,iSpin))
                call getExcMulliken(denseDesc, pc(:,:,iSpin), SSqr, dqex(:,iSpin))
              end do

              call addNadiaGradients(sym, nxov_rd, this%nAtom, species0, env, denseDesc, norb,&
                & nocc_ud, getIA, getIJ, getAB, win, grndEigVecs, pc, ovrXev, dq, dqex, gammaMat, lrGamma,&
                & this%HubbardU, this%spinW, shift, woo, wov, wvv, transChrg, xpyn, xmyn, xpym,&
                & xmym, coord0, orb, skHamCont, skOverCont, derivator, rhoSqr, deltaRho, tHybridXc,&
                & hybridXc, nacv(:,:,iNac))

            end if

            !> P and W have not yet been divided by excitation energy difference
            nacv(:,:,iNac) = nacv(:,:,iNac) / omegaDif

          end do
        end do

        !> Convention to determine arbitrary phase
        call fixNACVPhase(nacv)

        call writeNACV(this%indNACouplings(1), this%indNACouplings(2), fdTagged, taggedWriter, nacv)

        call env%globalTimer%stopTimer(globalTimers%lrNAC)

      end if

    end if

    call env%globalTimer%stopTimer(globalTimers%lrZVector)
    
    !> Omega has possibly been overwritten for CI optimization or NA couplings, but should always
    !! refer to nstat
    if (nstat == 0) then
      omega = 0.0_dp
    else
      omega = sqrt(eval(nstat))
    end if

  end subroutine LinRespGrad_old


  !> Solves the RPA equations in their hermitian form (valid for local functionals) at finite T
  !!
  !!  [A  B] X   =    [C  0] X
  !!                W
  !!  [B  A] Y   =    [0 -C] Y
  !!
  !! (see definitions by Mark Casida, in Recent Advances in Density Functional Methods,
  !!  World Scientific, 1995, Part I, p. 155.)
  !!
  !! The hermitian EV problem is given by \Omega F = w^2 F, with
  !!  S = -C (A-B)^{-1} C, \Omega = - S^{-1/2} (A+B) S^{-1/2} and F = (X+Y) * sqrt(w/wia)
  !!
  !! In this routine \Omega is diagonalised by the iterative ARPACK diagonaliser.
  !! The code deals with closed shell systems by diagonalising dedicated singlet/triplet
  !! submatrices.
  !! See Dominguez JCTC 9 4901 (2013)
  subroutine buildAndDiagExcMatrixArpack(iGlobal, fGlobal, tSpin, wij, sym, win, nocc_ud, nvir_ud,&
      & nxoo_ud, nxvv_ud, nxov_ud, nxov_rd, iaTrans, getIA, getIJ, getAB, env, denseDesc, ovrXev,&
      & ovrXevGlb, grndEigVecs, eigVecGlb, filling, sqrOccIA, gammaMat, species0, spinW, transChrg,&
      & testArnoldi, eval, xpy, xmy, onsMEs, orb, tHybridXc, tZVector)

    !> Starting index of current rank in global RPA vectors
    integer, intent(in) :: iGlobal

    !> End index of current rank in global RPA vectors
    integer, intent(in) :: fGlobal
    
    !> Spin polarisation?
    logical, intent(in) :: tSpin

    !> Single particle excitation energies
    real(dp), intent(in) :: wij(:)

    !> Symmetry to calculate transitions
    character, intent(in) :: sym

    !> Index array for single particle excitations
    integer, intent(in) :: win(:)

    !> Occupied orbitals per spin channel
    integer, intent(in) :: nocc_ud(:)

    !> Virtual orbitals per spin channel
    integer, intent(in) :: nvir_ud(:)

    !> Number of occ-occ transitions per spin channel
    integer, intent(in) :: nxoo_ud(:)

    !> Number of vir-vir transitions per spin channel
    integer, intent(in) :: nxvv_ud(:)

    !> Number of occ-vir transitions per spin channel
    integer, intent(in) :: nxov_ud(:)

    !> Number of occupied-virtual transitions (possibly reduced by windowing)
    integer, intent(in) :: nxov_rd

    !> array from pairs of single particles states to compound index
    integer, intent(in) :: iaTrans(:,:,:)

    !> Index array for occ-vir single particle excitations
    integer, intent(in) :: getIA(:,:)

    !> Index array for occ-occ single particle excitations
    integer, intent(in) :: getIJ(:,:)

    !> Index array for vir-vir single particle excitations
    integer, intent(in) :: getAB(:,:)

    !> Environment settings
    type(TEnvironment), intent(inout) :: env

    !> Dense matrix descriptor
    type(TDenseDescr), intent(in) :: denseDesc

    !> overlap times ground state eigenvectors (local in case of MPI)
    real(dp), intent(in) :: ovrXev(:,:,:)

    !> overlap times ground state eigenvectors (global in case of MPI)  
    real(dp), allocatable, intent(in) :: ovrXevGlb(:,:,:)

    !> ground state eigenvectors (local in case of MPI)
    real(dp), intent(in) :: grndEigVecs(:,:,:)

    !> ground state eigenvectors (global in case of MPI)
    real(dp), allocatable, intent(in) :: eigVecGlb(:,:,:)
    
    !> occupation numbers
    real(dp), intent(in) :: filling(:,:)

    ! Square root of occupation difference between vir and occ states
    real(dp), intent(in) :: sqrOccIA(:)

    !> Electrostatic matrix
    real(dp), intent(in) :: gammaMat(:,:)

    !> Central cell chemical species
    integer, intent(in) :: species0(:)

    !> File handle for ARPACK eigenstate tests
    logical, intent(in) :: testArnoldi

    !> Atomic resolved spin constants
    real(dp), intent(in) :: spinW(:)

    !> Machinery for transition charges between single particle levels
    type(TTransCharges), intent(in) :: transChrg

    !> Resulting eigenvalues for transitions (w^2)
    real(dp), intent(out) :: eval(:)

    !> Eigenvectors (X+Y)
    real(dp), intent(out) :: xpy(:,:)

    !> Eigenvectors (X-Y), only evaluated if Z-vector is needed
    real(dp), intent(inout), allocatable :: xmy(:,:)

    !> Onsite corrections if in use
    real(dp), allocatable :: onsMEs(:,:,:,:)

    !> Data type for atomic orbital information
    type(TOrbitals), intent(in) :: orb

    !> Is calculation range-separated?
    logical, intent(in) :: tHybridXc

    !> Is the Z-vector equation to be solved later?
    logical, intent(in) :: tZVector

    real(dp), allocatable :: workl(:), workd(:), resid(:), vv(:,:), qij(:)
    real(dp) :: sigma, omega
    integer :: iparam(11), ipntr(11)
    integer :: ido, ncv, lworkl, info
    logical, allocatable :: selection(:)
    logical :: rvec
    integer :: nexc, natom, nLoc

    integer :: iState, comm
    real(dp), allocatable :: Hv(:), orthnorm(:,:)
    character(lc) :: tmpStr
    type(TFileDescr) :: fdArnoldiTest

  #:if WITH_SCALAPACK and WITH_ARPACK 
    external pdsaupd, pdseupd
    comm = env%mpi%globalComm%id 
  #:endif

    !! Local chunk of RPA vectors have this size under MPI
    nLoc = fGlobal - iGlobal + 1
    
    nexc = size(eval)
    natom = size(gammaMat, dim=1)

    @:ASSERT(all(shape(xpy) == [ nxov_rd, nexc ]))
    @:ASSERT(tHybridXc .eqv. .false.)

    ! Three times more Lanczos vectors than desired eigenstates
    ncv = min(3 * nexc, nxov_rd)

    lworkl = ncv * (ncv + 8)

    allocate(workl(lworkl))
    allocate(qij(natom))
    allocate(selection(ncv))        
    allocate(workd(3 * nLoc))
    allocate(resid(nLoc))
    allocate(vv(nLoc, ncv))

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
    
    do

      ! call the reverse communication interface from arpack
    #:if WITH_SCALAPACK and WITH_ARPACK

      call pdsaupd(comm, ido, "I", nLoc, "SM", nexc, ARTOL, resid, ncv, vv, nLoc, iparam,&
          & ipntr, workd, workl, lworkl, info)
 
    #:else
      
      call saupd(ido, "I", nxov_rd, "SM", nexc, ARTOL, resid, ncv, vv, nxov_rd, iparam,&
          & ipntr, workd, workl, lworkl, info)
      
    #:endif

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
      call actionAplusB(iGlobal, fGlobal, tSpin, wij, sym, win, nocc_ud, nvir_ud, nxoo_ud, nxvv_ud,&
          & nxov_ud, nxov_rd, iaTrans, getIA, getIJ, getAB, env, denseDesc, ovrXev, ovrXevGlb,&
          & grndEigVecs, eigVecGlb, filling, sqrOccIA, gammaMat, species0, spinW, onsMEs, orb,&
          & .false., transChrg, workd(ipntr(1):ipntr(1)+nLoc-1), workd(ipntr(2):ipntr(2)+nLoc-1),&
          & tHybridXc)
      
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
      ! call to DSEUPD.s
      ! Note: At this point xpy holds the hermitian eigenvectors F
    #:if WITH_SCALAPACK and WITH_ARPACK

      call pdseupd (comm, rvec, "All", selection, eval, vv, nLoc, sigma, "I", nLoc,& 
          & "SM", nexc, ARTOL, resid, ncv, vv, nLoc, iparam, ipntr, workd, workl, lworkl, info)

      xpy(:,:) = 0.0_dp
      xpy(iGlobal:fGlobal,:nexc) = vv(:,:nexc)
      call assembleChunks(env, xpy)

     #:else
       
      call seupd(rvec, "All", selection, eval, xpy, nxov_rd, sigma, "I", nxov_rd, "SM",&
          & nexc, ARTOL, resid, ncv, vv, nxov_rd, iparam, ipntr, workd, workl, lworkl, info)

    #:endif
            
      ! check for error on return
      if (info  /=  0) then
        write(tmpStr,"(' Error with ARPACK routine seupd, info = ',I0)")info
        call error(tmpStr)
      end if

    end if

    if (testArnoldi) then
      ! tests for quality of returned eigenpairs
      call openFile(fdArnoldiTest, testArpackOut, mode="w")
      allocate(Hv(nxov_rd))
      allocate(orthnorm(nxov_rd,nxov_rd))
      orthnorm = matmul(transpose(xpy(:,:nExc)),xpy(:,:nExc))

      write(fdArnoldiTest%unit,"(A)")'State Ei deviation    Evec deviation  Norm deviation  Max&
          & non-orthog'
      do iState = 1, nExc

        call actionAplusB(iGlobal, fGlobal, tSpin, wij, sym, win, nocc_ud, nvir_ud, nxoo_ud,&
            & nxvv_ud, nxov_ud, nxov_rd, iaTrans, getIA, getIJ, getAB, env, denseDesc, ovrXev,&
            & ovrXevGlb, grndEigVecs, eigVecGlb, filling, sqrOccIA, gammaMat, species0, spinW,&
            & onsMEs, orb, .false., transChrg, xpy(iGlobal:fGlobal,iState), Hv(iGlobal:fGlobal),& 
            & .false.)
        call assembleChunks(env, Hv)
        
        write(fdArnoldiTest%unit,"(I4,4E16.8)")iState,&
            & dot_product(Hv,xpy(:,iState))-eval(iState),&
            & sqrt(sum( (Hv-xpy(:,iState)*eval(iState) )**2 )), orthnorm(iState,iState) - 1.0_dp,&
            & max(maxval(orthnorm(:iState-1,iState)), maxval(orthnorm(iState+1:,iState)))
      end do
      call closeFile(fdArnoldiTest)
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
  !!
  !!  [A  B] X   =    [C  0] X
  !!                W
  !!  [B  A] Y   =    [0 -C] Y
  !!
  !! (see definitions by Mark Casida, in Recent Advances in Density Functional Methods,
  !!  World Scientific, 1995, Part I, p. 155.)
  !!
  !! The RPA eqs are diagonalised by the Stratmann algorithm (JCP 109 8218 (1998).
  !! See also Dominguez JCTC 9 4901 (2013), Kranz JCTC 13 1737 (2017)
  !!
  !! Returns w^2 and (X+Y) (to be consistent with ARPACK diagonaliser)
  !!
  subroutine buildAndDiagExcMatrixStratmann(iGlobal, fGlobal, tSpin, subSpaceFactor, wij, sym, win,&
      & nocc_ud, nvir_ud, nxoo_ud, nxvv_ud, nxov_ud, nxov_rd, iaTrans, getIA, getIJ, getAB, env,&
      & denseDesc, ovrXev, ovrXevGlb, grndEigVecs, eigVecGlb, filling, sqrOccIA, gammaMat, species0,&
      & spinW, transChrg, eval, xpy, xmy, onsMEs, orb, tHybridXc, lrGamma, tZVector)

    !> Starting index of current rank in global RPA vectors
    integer, intent(in) :: iGlobal

    !> End index of current rank in global RPA vectors
    integer, intent(in) :: fGlobal
    
    !> Spin polarisation?
    logical, intent(in) :: tSpin

    !> Initial subspace is this factor times number of excited states
    integer :: subSpaceFactor

    !> Single particle excitation energies
    real(dp), intent(in) :: wij(:)

    !> Symmetry to calculate transitions
    character, intent(in) :: sym

    !> Index array for single particle excitations
    integer, intent(in) :: win(:)

    !> Occupied orbitals per spin channel
    integer, intent(in) :: nocc_ud(:)

    !> Virtual orbitals per spin channel
    integer, intent(in) :: nvir_ud(:)

    !> Number of occ-occ transitions per spin channel
    integer, intent(in) :: nxoo_ud(:)

    !> Number of vir-vir transitions per spin channel
    integer, intent(in) :: nxvv_ud(:)

    !> Number of occ-vir transitions per spin channel
    integer, intent(in) :: nxov_ud(:)

    !> Number of occupied-virtual transitions (possibly reduced by windowing)
    integer, intent(in) :: nxov_rd
    
    !> Array from pairs of single particles states to compound index
    integer, intent(in) :: iaTrans(:,:,:)

    !> Index array for occ-vir single particle excitations
    integer, intent(in) :: getIA(:,:)

    !> Index array for occ-occ single particle excitations
    integer, intent(in) :: getIJ(:,:)

    !> Index array for vir-vir single particle excitations
    integer, intent(in) :: getAB(:,:)

    !> Environment settings
    type(TEnvironment), intent(inout) :: env

    !> Dense matrix descriptor
    type(TDenseDescr), intent(in) :: denseDesc 

    !> overlap times ground state eigenvectors (local in case of MPI)
    real(dp), intent(in) :: ovrXev(:,:,:)

    !> overlap times ground state eigenvectors (global in case of MPI)  
    real(dp), allocatable, intent(in) :: ovrXevGlb(:,:,:)

    !> ground state eigenvectors (local in case of MPI)
    real(dp), intent(in) :: grndEigVecs(:,:,:)

    !> ground state eigenvectors (global in case of MPI)
    real(dp), allocatable, intent(in) :: eigVecGlb(:,:,:)

    !> occupation numbers
    real(dp), intent(in) :: filling(:,:)

    ! Square root of occupation difference between vir and occ states
    real(dp), intent(in) :: sqrOccIA(:)

    !> Electrostatic matrix
    real(dp), intent(in) :: gammaMat(:,:)

    !> Central cell chemical species
    integer, intent(in) :: species0(:)

    !> Atomic resolved spin constants
    real(dp), intent(in) :: spinW(:)

    !> Machinery for transition charges between single particle levels
    type(TTransCharges), intent(in) :: transChrg

    !> Resulting eigenvalues for transitions (actually w^2)
    real(dp), intent(out) :: eval(:)

    !> Eigenvectors (X+Y)
    real(dp), intent(out) :: xpy(:,:)

    !> Eigenvectors (X-Y), only evaluated if Z-vector is needed
    real(dp), intent(inout), allocatable :: xmy(:,:)

    !> Onsite corrections if in use
    real(dp), allocatable :: onsMEs(:,:,:,:)

    !> Data type for atomic orbital information
    type(TOrbitals), intent(in) :: orb

    !> Is calculation range-separated?
    logical, intent(in) :: tHybridXc

    !> Electrostatic matrix, long-range corrected
    real(dp), allocatable, intent(in) :: lrGamma(:,:)

    !> Is the Z-vector equation to be solved later?
    logical, intent(in) :: tZVector

    real(dp), allocatable :: vecB(:,:) ! basis of subspace
    real(dp), allocatable :: evecL(:,:), evecR(:,:) ! left and right eigenvectors of Mnh
    real(dp), allocatable :: vP(:,:), vM(:,:) ! vec. for (A+B)b_i, (A-B)b_i
    ! matrices M_plus, M_minus, M_minus^(1/2), M_minus^(-1/2) and M_herm~=resp. mat on subspace
    real(dp), allocatable :: mP(:,:), mM(:,:), mMsqrt(:,:), mMsqrtInv(:,:), mH(:,:) 
    ! Residual vectors
    real(dp), allocatable :: resR(:,:), resL(:,:)
    real(dp), allocatable :: evalInt(:) ! store eigenvectors within routine
    real(dp), allocatable :: dummyM(:,:)
    real(dp), allocatable :: vecNorm(:) ! will hold norms of residual vectors
    real(dp) :: dummyReal

    integer :: nExc, nAtom, dummyInt, newVec, iVec, info, iterStrat, nLoc
    integer :: subSpaceDim, prevSubSpaceDim
    integer :: ii, jj, myjj, myii
    character(lc) :: tmpStr

    logical :: didConverge
    
  #:if WITH_SCALAPACK
    external pdsaupd, pdseupd
  #:endif
    
    !! Local chunk of RPA vectors have this size under MPI
    nLoc = fGlobal - iGlobal + 1
    
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
      write(tmpStr,'(A)') 'SubSpaceFactor for Stratmann solver must be larger than one.'
      call error(tmpStr)
    endif
    subSpaceDim = min(subSpaceFactor * nExc, nxov_rd)
    iterStrat = 1

    write(stdOut,'(A)')
    write(stdOut,'(A)') '>> Stratmann diagonalisation of response matrix'
    write(stdOut,'(3x,A,i6,A,i6)') 'Total dimension of A+B: ', nxov_rd, ' inital subspace: ',&
      & subSpaceDim

    allocate(mP(subSpaceDim, subSpaceDim))
    allocate(mM(subSpaceDim, subSpaceDim))
    allocate(mMsqrt(subSpaceDim, subSpaceDim))
    allocate(mMsqrtInv(subSpaceDim, subSpaceDim))
    allocate(mH(subSpaceDim, subSpaceDim))
    allocate(dummyM(subSpaceDim, subSpaceDim))
    allocate(evalInt(subSpaceDim))
    allocate(evecL(subSpaceDim, nExc))
    allocate(evecR(subSpaceDim, nExc))
    allocate(vecNorm(2*nExc))

    allocate(vecB(nLoc, subSpaceDim))
    allocate(vP(nLoc, subSpaceDim))
    allocate(vM(nLoc, subSpaceDim))
    allocate(resL(nLoc, nExc))
    allocate(resR(nLoc, nExc))   

    ! set initial bs
    vecB(:,:) = 0.0_dp
    do ii = iGlobal, fGlobal
      myii = ii - iGlobal + 1
      if(ii > subSpaceDim) exit
      vecB(myii, ii) = 1.0_dp
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

          call actionAplusB(iGlobal, fGlobal, tSpin, wij, sym, win, nocc_ud, nvir_ud, &
            & nxoo_ud, nxvv_ud, nxov_ud, nxov_rd, iaTrans, getIA, getIJ, getAB, env, denseDesc,&
            & ovrXev, ovrXevGlb, grndEigVecs, eigVecGlb, filling, sqrOccIA, gammaMat, species0,&
            & spinW, onsMEs, orb, .true., transChrg, vecB(:,ii), vP(:,ii), tHybridXc, lrGamma)
          call actionAminusB(iGlobal, fGlobal, tSpin, wij, win, nocc_ud, nvir_ud, nxoo_ud,&
            & nxvv_ud, nxov_ud, nxov_rd, iaTrans, getIA, getIJ, getAB, env, denseDesc, ovrXev,&
            & grndEigVecs, filling, sqrOccIA, transChrg, vecB(:,ii), vM(:,ii), tHybridXc,& 
            & lrGamma)

        end do
 
       do ii = prevSubSpaceDim + 1, subSpaceDim
          do jj = 1, ii
            dummyReal = dot_product(vecB(:,jj), vP(:,ii))
            call assembleChunks(env, dummyReal)
            mP(ii,jj) = dummyReal
            mP(jj,ii) = mP(ii,jj)
            dummyReal = dot_product(vecB(:,jj), vM(:,ii))
            call assembleChunks(env, dummyReal)
            mM(ii,jj) = dummyReal
            mM(jj,ii) = mM(ii,jj)
          end do
        end do

      else
        ! We need (A+B)_iajb. Could be realized by calls to actionAplusB.
        ! Specific routine for this task is more effective
        call initialSubSpaceMatrixApmB(iGlobal, fGlobal, transChrg, subSpaceDim, wij, sym, win,&
            & nxov_ud(1), env, denseDesc, ovrXev, grndEigVecs, filling, sqrOccIA, getIA, getIJ,&
            & getAB, iaTrans, gammaMat, lrGamma, species0, spinW, tSpin, tHybridXc, vP, vM, mP, mM)

      end if

      call calcMatrixSqrt(mM, subSpaceDim, mMsqrt, mMsqrtInv)
      
      call symm(dummyM, 'L', mP, mMsqrt, uplo='U')
      call symm(mH, 'L', mMsqrt, dummyM, uplo='U')

      ! Diagonalise in subspace
      call heev(mH, evalInt, 'U', 'V', info)
      if (info /= 0) then
        if (subSpaceFactor * nExc < nxov_rd) then
          write(tmpStr,'(A)') 'TDDFT diagonalisation failure. Increase SubSpaceFactor.'
        else
          write(tmpStr,'(A)') 'TDDFT diagonalisation failure. Insufficient transitions available to&
              & converge.'
        end if
        call error(tmpStr)
      endif

      ! This yields T=(A-B)^(-1/2)|X+Y>.
      ! Calc. |R_n>=|X+Y>=(A-B)^(1/2)T and |L_n>=|X-Y>=(A-B)^(-1/2)T.
      ! Transformation preserves orthonormality.
      ! Only compute up to nExc index, because only that much needed.
      call symm(evecR, 'L', Mmsqrt, Mh, uplo='U')
      call symm(evecL, 'L', Mmsqrtinv, Mh, uplo='U')

      ! Need |X-Y>=sqrt(w)(A-B)^(-1/2)T, |X+Y>=(A-B)^(1/2)T/sqrt(w) for proper solution to original
      ! EV problem, only use first nExc vectors
      do ii = 1, nExc
        dummyReal = sqrt(sqrt(evalInt(ii)))
        evecR(:,ii) = evecR(:,ii) / dummyReal
        evecL(:,ii) = evecL(:,ii) * dummyReal
      end do

      ! Calculate the residual vectors
      !   calcs. all |R_n>
      call gemm(resR, vecB, evecR)
      !   calcs. all |L_n>
      call gemm(resL, vecB, evecL)

      do ii = 1, nExc
        dummyReal = -sqrt(evalInt(ii))
        resR(:,ii) = dummyReal * resR(:,ii)
        resL(:,ii) = dummyReal * resL(:,ii)
      end do

      ! (A-B)|L_n> for all n=1,..,nExc
      call gemm(resR, vM, evecL, beta=1.0_dp)
      ! (A+B)|R_n> for all n=1,..,nExc
      call gemm(resL, vP, evecR, beta=1.0_dp)

      ! calc. norms of residual vectors to check for convergence
      do ii = 1, nExc
        dummyReal = dot_product(resR(:,ii), resR(:,ii))
        call assembleChunks(env, dummyReal)
        vecNorm(ii) = dummyReal
        dummyReal = dot_product(resL(:,ii), resL(:,ii))
        call assembleChunks(env, dummyReal)
        vecNorm(nExc+ii) = dummyReal
      end do
      didConverge = all(vecNorm .lt. CONV_THRESH_STRAT)
      
      if ((.not. didConverge) .and. (subSpaceDim > nxov_rd)) then
        write(tmpStr,'(A)') 'Linear Response calculation in subspace did not converge!&
             & Increase SubspaceFactor.'
        call error(tmpStr)
      end if

      ! if converged then exit loop:
      if (didConverge) then
        eval(:) = evalInt(1:nExc)

        ! Calc. X+Y
        xpy(:,:) = 0.0_dp
        xpy(iGlobal:fGlobal,:) = matmul(vecB, evecR)
        call assembleChunks(env, xpy)
        
        ! Calc. X-Y, only when needed
        if (tZVector) then
          xmy(iGlobal:fGlobal,:) = matmul(vecB, evecL)
          call assembleChunks(env, xmy)
        end if

        write(stdOut,'(A)') '>> Stratmann converged'
        exit solveLinResp ! terminate diag. routine
      end if

      ! Otherwise calculate new basis vectors and extend subspace with them
      ! only include new vectors if they add meaningful residue component
      newVec = 0
      do ii = 1, 2*nExc
        if (vecNorm(ii) .gt. CONV_THRESH_STRAT) then
          newVec = newVec + 1
        endif
      enddo
      
      call incMemStratmann(subSpaceDim, subSpaceDim + newVec, vecB, vP, vM, mP, mM, mH, mMsqrt,&
            &  mMsqrtInv, dummyM, evalInt, evecL, evecR)

      iVec = 0
      do ii = 1, nExc
        if (vecNorm(ii) .gt. CONV_THRESH_STRAT) then
          iVec = iVec + 1
          dummyReal = sqrt(evalInt(ii))
          dummyInt = subSpaceDim + iVec

          do jj = iGlobal, fGlobal
            myjj = jj - iGlobal + 1
            vecB(myjj,dummyInt) = resR(myjj,ii) / (dummyReal - wij(jj))
          end do
          
        end if
      end do

      do ii = 1, nExc
        if (vecNorm(nExc+ii) .gt. CONV_THRESH_STRAT) then
          iVec = iVec + 1
          dummyInt = subSpaceDim + iVec

          do jj = iGlobal, fGlobal
            myjj = jj - iGlobal + 1
            vecB(myjj,dummyInt) = resL(myjj,ii) / (dummyReal - wij(jj))
          end do

        end if
      end do
      
      prevSubSpaceDim = subSpaceDim
      subSpaceDim = subSpaceDim + newVec

      if(iterStrat == 1) then
        write(stdOut,'(3x,A)') 'Iteration  Subspace dimension'
      end if
      write(stdOut,'(3x,i6,10x,i6)') iterStrat, subSpaceDim

      iterStrat = iterStrat + 1

      ! create orthogonal basis
      call orthonormalizeVectors(env, prevSubSpaceDim + 1, subSpaceDim, vecB)

    end do solveLinResp

  end subroutine buildAndDiagExcMatrixStratmann


  !> Calculate oscillator strength for a given excitation between KS states.
  subroutine getOscillatorStrengths(sym, tSpin, snglPartTransDip, eval, xpy, sqrOccIA, istat, osz,&
      & tTradip, transitionDipoles)

    !> Symmetry of transition
    character, intent(in) :: sym

    !> Spin polarisation?
    logical, intent(in) :: tSpin

    !> Dipole moments for single particle transtions
    real(dp), intent(in) :: snglPartTransDip(:,:)

    !> Low lying eigenvalues of Casida eqn (Omega^2)
    real(dp), intent(in) :: eval(:)

    !> Eigenvectors of Casida eqn (X+Y)
    real(dp), intent(in) :: xpy(:,:)

    !> Square root of KS occupation differences
    real(dp), intent(in) :: sqrOccIA(:)

    !> Write transition dipole
    logical :: tTradip

    !> Flag wich if <-1 on entry is returned as the brightest state
    integer, intent(inout) :: istat

    !> Oscilator strengths of transitions
    real(dp), intent(out) :: osz(:)

    !> Resulting transition dipoles
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
   subroutine getZVectorEqRHS(tHybridXc, xpy, xmy, win, env, denseDesc, homo, transChrg, getIA, getIJ,&
      & getAB, iatrans, natom, species0, grndEigVal, ovrXev, grndEigVecs, gammaMat, lrGamma,&
      & spinW, omega, sym, rhs, t, wov, woo, wvv)

    !> Is calculation range-separated?
    logical, intent(in) :: tHybridXc

    !> X+Y Furche term
    real(dp), intent(in) :: xpy(:)

    !> X-Y Furche term
    real(dp), intent(in) :: xmy(:)

    !> Index array for single particle transitions
    integer, intent(in) :: win(:)

    !> Environment settings
    type(TEnvironment), intent(inout) :: env

    !> Dense matrix descriptor
    type(TDenseDescr), intent(in) :: denseDesc 

    !> Highest occupied level
    integer, intent(in) :: homo(:)

    !> Machinery for transition charges between single particle levels
    type(TTransCharges), intent(in) :: transChrg

    !> Index array between transitions in square and 1D representations
    integer, intent(in) :: getIA(:,:)

    !> Index array for vir-vir transitions
    integer, intent(in) :: getIJ(:,:)

    !> Index array for occ-occ transitions
    integer, intent(in) :: getAB(:,:)

    !> Index array from orbital pairs to compound index
    integer, intent(in) :: iatrans(:,:,:)

    !> Number of central cell atoms
    integer, intent(in) :: natom

    !> Central cell chemical species
    integer, intent(in) :: species0(:)

    !> Ground state wavefunctions
    real(dp), intent(in) :: grndEigVal(:,:)

    !> Overlap times ground state wavefunctions
    real(dp), intent(in) :: ovrXev(:,:,:)

    !> Ground state wavefunctions
    real(dp), intent(in) :: grndEigVecs(:,:,:)

    !> Softened coulomb matrix
    real(dp), intent(in) :: gammaMat(:,:)

    !> Softened coulomb matrix, long-range corrected
    real(dp), allocatable, intent(in) :: lrGamma(:,:)

    !> Ground state spin derivatives for each species
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
    real(dp), allocatable :: vecHovT(:), vecHooT(:)
    integer :: nxov
    integer, allocatable :: nxoo(:), nxvv(:), nvir(:)
    integer :: i, j, a, b, ias, ibs, abs, ij, ab, jas, ijs, s, nSpin, soo(2), svv(2), nOrb
    real(dp) :: tmp1, tmp2, fact
    logical :: tSpin

    nxov = size(rhs)
    nOrb = size(ovrXev, dim=1)

    allocate(xpyq(natom))
    allocate(qTr(natom))
    allocate(gamxpyq(natom))
    allocate(gamqt(natom))

    t(:,:,:) = 0.0_dp
    rhs(:) = 0.0_dp
    wov(:) = 0.0_dp
    woo(:,:) = 0.0_dp
    wvv(:,:) = 0.0_dp

    nSpin = size(t, dim=3)

    allocate(nxoo(nSpin))
    allocate(nxvv(nSpin))
    allocate(nvir(nSpin))

    nxoo(:) = (homo(:)*(homo(:)+1))/2
    nvir(:) = size(t, dim=1) - homo(:)
    nxvv(:) = (nvir(:)*(nvir(:)+1))/2

    !! transition charges use compound index ijs = ij + soo(s)
    soo(:) = [0, nxoo(1)]
    svv(:) = [0, nxvv(1)]

    allocate(qgamxpyq(max(maxval(nxoo), maxval(nxvv)), size(homo)))

    if (nSpin == 2) then
      tSpin = .true.
      allocate(xpyqds(natom))
      allocate(gamxpyqds(natom))
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
    call transChrg%qMatVec(env, denseDesc, ovrXev, grndEigVecs, getIA, win, xpy, xpyq)

    if (.not. tSpin) then  ! ---- spin-unpolarized case ----
      ! qgamxpyq(ab) = sum_jc K_ab,jc (X+Y)_jc
      if (sym == "S") then
        call hemv(gamxpyq, gammaMat,  xpyq)
        do ab = 1, nxvv(1)
          qTr(:) = transChrg%qTransAB(ab, env, denseDesc, ovrXev, grndEigVecs, getAB)
          qgamxpyq(ab, 1) = 2.0_dp * sum(qTr * gamxpyq)
        end do
      else ! triplet case
        do ab = 1, nxvv(1)
          qTr(:) = transChrg%qTransAB(ab, env, denseDesc, ovrXev, grndEigVecs, getAB)
          qgamxpyq(ab, 1) = 2.0_dp * sum(qTr * xpyq * spinW(species0))
        end do
      end if

    else  ! ---- spin-polarized case -----

      xpyqds(:) = 0.0_dp
      call transChrg%qMatVecDs(env, denseDesc, ovrXev, grndEigVecs, getIA, win, xpy, xpyqds)

      call hemv(gamxpyq, gammaMat,  xpyq)
      do s = 1, 2
        if (s == 1) then
          fact = 1.0_dp
        else
          fact = -1.0_dp
        end if
        do ab = 1, nxvv(s)
          qTr(:) = transChrg%qTransAB(ab + svv(s), env, denseDesc, ovrXev, grndEigVecs, getAB)
          qgamxpyq(ab, s) = sum(qTr * gamxpyq)
          !magnetization part
          qgamxpyq(ab, s) = qgamxpyq(ab, s) + fact * sum(qTr * xpyqds * spinW(species0))
        end do
      end do

    end if

    ! rhs(ia) -= Qia = sum_b (X+Y)_ib * qgamxpyq(ab))
    do ias = 1, nxov
      call indxov(win, ias, getIA, i, a, s)

      do b = homo(s) + 1, a
        ab = iatrans(a, b, s) - svv(s)
        ibs = iatrans(i, b, s)
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
          qTr(:) = transChrg%qTransIJ(ij, env, denseDesc, ovrXev, grndEigVecs, getIJ)
          ! qgamxpyq(ij) = sum_kb K_ij,kb (X+Y)_kb
          qgamxpyq(ij, 1) = 2.0_dp * sum(qTr * gamxpyq)
        end do
      else
        do ij = 1, nxoo(1)
          qgamxpyq(ij, 1) = 0.0_dp
          qTr(:) = transChrg%qTransIJ(ij, env, denseDesc, ovrXev, grndEigVecs, getIJ)
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
          qTr(:) = transChrg%qTransIJ(ij + soo(s), env, denseDesc, ovrXev, grndEigVecs, getIJ)
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
        qTr(:) = transChrg%qTransIJ(ij + soo(s), env, denseDesc, ovrXev, grndEigVecs, getIJ)
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
        qTr(:) = transChrg%qTransAB(ab + svv(s), env, denseDesc, ovrXev, grndEigVecs, getAB)
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
      call transChrg%qVecMat(env, denseDesc, ovrXev, grndEigVecs, getIA, win, -4.0_dp*gamqt, rhs)
    else
      call transChrg%qVecMat(env, denseDesc, ovrXev, grndEigVecs, getIA, win, -2.0_dp*gamqt, rhs)
      call transChrg%qVecMatDs(env, denseDesc, ovrXev, grndEigVecs, getIA, win, &
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
        qTr(:) = transChrg%qTransIJ(ij + soo(s), env, denseDesc, ovrXev, grndEigVecs, getIJ)
        if (.not. tSpin) then
          woo(ij,s) = woo(ij,s) + 4.0_dp * sum(qTr * gamqt)
        else
          woo(ij,s) = woo(ij,s) + 2.0_dp * sum(qTr * gamqt)
          woo(ij,s) = woo(ij,s) + 2.0_dp * fact * sum(qTr * gamxpyqds * spinW(species0))
        end if
      end do
    end do

    ! Contributions due to range-separation
    if (tHybridXc) then

      allocate(vecHvvXpY(sum(nxvv)))
      allocate(vecHvvXmY(sum(nxvv)))
      allocate(vecHooXpY(sum(nxoo)))
      allocate(vecHooXmY(sum(nxoo)))
      allocate(vecHovT(nxov))
      allocate(vecHooT(sum(nxoo)))

      call getHvvXY( 1, nxvv, homo, natom, iatrans, getIA, getAB, win, env,&
          & denseDesc, ovrXev, grndEigVecs, lrGamma, transChrg, xpy, vecHvvXpY)

      call getHvvXY(-1, nxvv, homo, natom, iatrans, getIA, getAB, win, env,&
          & denseDesc, ovrXev, grndEigVecs, lrGamma, transChrg, xmy, vecHvvXmY)

      call getHooXY( 1, nxoo, homo, natom, iatrans, getIA, getIJ, win, env,&
          & denseDesc, ovrXev, grndEigVecs, lrGamma, transChrg, xpy, vecHooXpY)

      call getHooXY(-1, nxoo, homo, natom, iatrans, getIA, getIJ, win, env,&
          & denseDesc, ovrXev, grndEigVecs, lrGamma, transChrg, xmy, vecHooXmY)

      call getHovT(nxoo, nxvv, homo, natom, iatrans, getIA, getIJ, getAB, win, env,&
        & denseDesc, ovrXev, grndEigVecs, lrGamma, transChrg, t, vecHovT)

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

      call getHooT(nxov, nxoo, homo, natom, iatrans, getIA, getIJ, win, env, &
       & denseDesc, ovrXev, grndEigVecs, lrGamma, transChrg, t, vecHooT)

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


  !> Solving the (A+B) Z = -R equation via diagonally preconditioned conjugate gradient.
  subroutine solveZVectorPrecond(rhs, tSpin, wij, win, nocc_ud, nvir_ud, nxoo_ud, nxvv_ud,&
      & nxov_ud, nxov_rd, iaTrans, getIA, getIJ, getAB, natom, env, denseDesc, ovrXev, ovrXevGlb,&
      & grndEigVecs, eigVecGlb, occNr, sqrOccIA, gammaMat, species0, spinW, onsMEs, orb, transChrg,& 
      & tHybridXc, lrGamma)

    !> On entry -R, on exit Z
    real(dp), intent(inout) :: rhs(:)

    !> Logical spin polarization
    logical, intent(in) :: tSpin

    !> Excitation energies (wij = epsion_j - epsilon_i)
    real(dp), intent(in) :: wij(:)

    !> Sorting index of the excitation energies
    integer, intent(in) :: win(:)

    !> Occupied orbitals per spin channel
    integer, intent(in) :: nocc_ud(:)

    !> Virtual orbitals per spin channel
    integer, intent(in) :: nvir_ud(:)

    !> Number of occ-occ transitions per spin channel
    integer, intent(in) :: nxoo_ud(:)

    !> Number of vir-vir transitions per spin channel
    integer, intent(in) :: nxvv_ud(:)

    !> Number of occ-vir transitions per spin channel
    integer, intent(in) :: nxov_ud(:)

    !> Number of occupied-virtual transitions (possibly reduced by windowing)
    integer, intent(in) :: nxov_rd

    !> Array from pairs of single particles states to compound index
    integer, intent(in) :: iaTrans(:,:,:)

    !> Index array for occ-vir single particle excitations
    integer, intent(in) :: getIA(:,:)

    !> Index array for occ-occ single particle excitations
    integer, intent(in) :: getIJ(:,:)

    !> Index array for vir-vir single particle excitations
    integer, intent(in) :: getAB(:,:)

    !> Number of atoms
    integer, intent(in) :: natom

    !> Environment settings
    type(TEnvironment), intent(inout) :: env

    !> Dense matrix descriptor
    type(TDenseDescr), intent(in) :: denseDesc 

    !> Overlap times eigenvector (nOrb, nOrb)
    real(dp), intent(in) :: ovrXev(:,:,:)
    
    !> overlap times ground state eigenvectors (global in case of MPI)  
    real(dp), allocatable, intent(in) :: ovrXevGlb(:,:,:)   

    !> Eigenvectors (nOrb, nOrb)
    real(dp), intent(in) :: grndEigVecs(:,:,:)

    !> ground state eigenvectors (global in case of MPI)
    real(dp), allocatable, intent(in) :: eigVecGlb(:,:,:)

    !> Occupation numbers
    real(dp), intent(in) :: occNr(:,:)

    ! Square root of occupation difference between vir and occ states
    real(dp), intent(in) :: sqrOccIA(:)

    !> DFTB gamma matrix (nAtm, nAtom)
    real(dp), intent(in) :: gammaMat(:,:)

    !> Chemical species of the atoms
    integer, intent(in) :: species0(:)

    !> Ground state spin constants for each species
    real(dp), intent(in) :: spinW(:)

    !> Onsite matrix elements for shells (elements between s orbitals on the same shell are ignored)
    real(dp), intent(in), allocatable :: onsMEs(:,:,:,:)

    !> Data type for atomic orbital information
    type(TOrbitals), intent(in) :: orb

    !> Machinery for transition charges between single particle levels
    type(TTransCharges), intent(in) :: transChrg

    !> Is calculation range-separated?
    logical, intent(in) :: tHybridXc

    !> Long-range Gamma
    real(dp), allocatable, intent(in) :: lrGamma(:,:)

    integer :: nxov, iGlobal, fGlobal
    integer :: ia, kk, i, a, s, iis, aas
    real(dp) :: rhs2(size(rhs)), rkm1(size(rhs)), zkm1(size(rhs)), pkm1(size(rhs)), apk(size(rhs))
    real(dp) :: qTmp(nAtom), rs, alphakm1, tmp1, tmp2, bkm1
    real(dp), allocatable :: qTr(:), P(:)

    ! Need to change this for MPI, currently routine works on full RPA vectors
    iGlobal = 1
    fGlobal = nxov_rd
    
    nxov = nxov_rd
    allocate(qTr(nAtom))

    ! diagonal preconditioner
    ! P^-1 = 1 / (A+B)_ia,ia (diagonal of the supermatrix sum A+B)
    allocate(P(nxov))
    do ia = 1, nxov
      qTr = transChrg%qTransIA(ia, env, denseDesc, ovrXev, grndEigVecs, getIA, win)
      call hemv(qTmp, gammaMat, qTr)
      if (.not. tSpin) then
        rs = 4.0_dp * dot_product(qTr, qTmp) + wij(ia)
      else
        rs = 2.0_dp * dot_product(qTr, qTmp) + wij(ia)
        rs = rs + 2.0_dp * sum(qTr * qTr * spinW(species0))
      end if

      !! Possibly reorder spin case
      if (tHybridXc) then
        call hemv(qTmp, lrGamma, qTr)
        rs = rs - cExchange * dot_product(qTr, qTmp)
        call indXov(win, ia, getIA, i, a, s)
        iis = iaTrans(i, i, s)
        qTr = transChrg%qTransIJ(iis, env, denseDesc, ovrXev, grndEigVecs, getIJ)
        call hemv(qTmp, lrGamma, qTr)
        aas = iaTrans(a, a, s)
        qTr = transChrg%qTransAB(aas, env, denseDesc, ovrXev, grndEigVecs, getAB)
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
    call actionAplusB(iGlobal, fGlobal, tSpin, wij, 'S', win, nocc_ud, nvir_ud, nxoo_ud, nxvv_ud,&
      & nxov_ud, nxov_rd, iaTrans, getIA, getIJ, getAB, env, denseDesc, ovrXev, ovrXevGlb, &
      & grndEigVecs, eigVecGlb, occNr, sqrOccIA, gammaMat, species0, spinW, onsMEs, orb, .true.,& 
      & transChrg, rhs2, rkm1, tHybridXc, lrGamma)
    
    rkm1(:) = rhs - rkm1
    zkm1(:) = P * rkm1
    pkm1(:) = zkm1

    ! Iteration: should be convergent in at most nxov steps for a quadradic surface, so set higher
    do kk = 1, nxov**2

      ! action of matrix on vector
      call actionAplusB(iGlobal, fGlobal, tSpin, wij, 'S', win, nocc_ud, nvir_ud, nxoo_ud, nxvv_ud,&
         & nxov_ud, nxov_rd, iaTrans, getIA, getIJ, getAB, env, denseDesc, ovrXev, ovrXevGlb,&
         & grndEigVecs, eigVecGlb, occNr, sqrOccIA, gammaMat, species0, spinW, onsMEs, orb, .true.,& 
         & transChrg, pkm1, apk, tHybridXc, lrGamma)
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
  subroutine calcWvectorZ(zz, win, homo, getIA, getIJ, getAB, iaTrans, env, denseDesc, ovrXev,&
      & grndEigVecs, gammaMat, grndEigVal, wov, woo, wvv, transChrg, species0, spinW,&
      & tHybridXc, lrGamma)

    !> Z vector
    real(dp), intent(in) :: zz(:)

    !> Index array for single particle transitions
    integer, intent(in) :: win(:)

    !> Highest occupied level
    integer, intent(in) :: homo(:)

    !> Index array between occ-vir transitions in square and 1D representations
    integer, intent(in) :: getIA(:,:)

    !> Index array between occ-occ transitions
    integer, intent(in) :: getIJ(:,:)

    !> Index array between vir-vir transitions
    integer, intent(in) :: getAB(:,:)

    !> Array from pairs of single particles states to compound index
    integer, intent(in) :: iaTrans(:,:,:)

    !> Environment settings
    type(TEnvironment), intent(inout) :: env

    !> Dense matrix descriptor
    type(TDenseDescr), intent(in) :: denseDesc 

    !> Overlap times ground state wavefunctions
    real(dp), intent(in) :: ovrXev(:,:,:)

    !> Ground state wavefunctions
    real(dp), intent(in) :: grndEigVecs(:,:,:)

    !> Softened coulomb matrix
    real(dp), intent(in) :: gammaMat(:,:)

    !> Ground state MO-energies
    real(dp), intent(in) :: grndEigVal(:,:)

    !> W vector occupied-virtual part
    real(dp), intent(inout) :: wov(:)

    !> W vector occupied part
    real(dp), intent(inout) :: woo(:,:)

    !> W vector virtual part
    real(dp), intent(inout) :: wvv(:,:)

    !> Machinery for transition charges between single particle levels
    type(TTransCharges), intent(in) :: transChrg

    !> Central cell chemical species
    integer, intent(in) :: species0(:)

    !> Ground state spin derivatives for each species
    real(dp), intent(in) :: spinW(:)

    !> Is calculation range-separated?
    logical, intent(in) :: tHybridXc

    !> Long-range Gamma
    real(dp), allocatable, intent(in) :: lrGamma(:,:)

    integer :: nxov, natom, nSpin, soo(2), svv(2)
    integer, allocatable :: nxoo(:), nxvv(:), nvir(:)
    integer :: ij, ias, ijs, ab, i, j, a, b, s
    real(dp) :: fact
    real(dp), allocatable :: qTr(:), gamxpyq(:), zq(:), zqds(:), vecHooZ(:)
    logical :: tSpin

    nxov = size(zz)
    natom = size(gammaMat, dim=1)
    nSpin = size(grndEigVal, dim=2)

    allocate(qTr(natom))
    allocate(gamxpyq(natom))
    allocate(zq(natom))
    allocate(nxoo(nSpin))
    allocate(nxvv(nSpin))
    allocate(nvir(nSpin))

    nxoo(:) = (homo(:)*(homo(:)+1))/2
    nvir(:) = size(grndEigVecs, dim=1) - homo(:)
    nxvv(:) = (nvir(:)*(nvir(:)+1))/2

    soo(:) = [0, nxoo(1)]
    svv(:) = [0, nxvv(1)]

    if ( nSpin == 2 ) then
      tSpin = .true.
      allocate(zqds(natom))
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
    call transChrg%qMatVec(env, denseDesc, ovrXev, grndEigVecs, getIA, win, zz, zq)
    call hemv(gamxpyq, gammaMat, zq)

    if (tSpin) then
      zqds(:) = 0.0_dp
      call transChrg%qMatVecDs(env, denseDesc, ovrXev, grndEigVecs, getIA, win, zz, zqds)
    end if

    ! sum_iAt1 qTr(iAt1) gamxpyq(iAt1)
    do s = 1, nSpin
      if (s == 1) then
        fact = 1.0_dp
      else
        fact = -1.0_dp
      end if
      do ij = 1, nxoo(s)
        qTr(:) = transChrg%qTransIJ(ij + soo(s), env, denseDesc, ovrXev, grndEigVecs, getIJ)
        ! W contains 1/2 for i == j.
        if (.not. tSpin) then
          woo(ij,s) = woo(ij,s) + 4.0_dp * sum(qTr * gamxpyq)
        else
          woo(ij,s) = woo(ij,s) + 2.0_dp * sum(qTr * gamxpyq)
          woo(ij,s) = woo(ij,s) + 2.0_dp * fact * sum(qTr * zqds * spinW(species0))
        end if
      end do
    end do

    if (tHybridXc) then

      allocate(vecHooZ(sum(nxoo)))
      call getHooXY(1, nxoo, homo, natom, iaTrans, getIA, getIJ, win,&
      & env, denseDesc, ovrXev, grndEigVecs, lrGamma, transChrg, zz, vecHooZ)

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


  !> Write out density matrix, full matrix if rhoSqr is present, but transition part only if not.
  subroutine writeDM(iLev, pc, rhoSqr)

    !> Lable for excited state level
    integer, intent(in) :: iLev

    !> Transition density matrix
    real(dp), intent(in) :: pc(:,:,:)

    !> Ground state density matrix
    real(dp), intent(in), optional :: rhoSqr(:,:,:)

    type(TFileDescr) :: fd
    integer :: iErr
    integer :: iSpin, nSpin
    character(lc) :: tmpStr, error_string

    nSpin = size(pc, dim=3)

    write(tmpStr, "(A,I0,A)")"DM", iLev, ".dat"

    call openFile(fd, trim(tmpStr), mode="wb", ioStat=iErr)
    if (iErr /= 0) then
      write(error_string, *) "Failure to open density matrix"
      call error(error_string)
    end if

    ! size and spin channels
    do iSpin = 1, nSpin
      write(fd%unit) size(pc, dim=1), iSpin

      if (present(rhoSqr)) then
        write(fd%unit) cmplx(pc(:,:,iSpin) + rhoSqr(:,:,iSpin), 0.0_dp, dp)
      else
        write(fd%unit) cmplx(pc(:,:,iSpin), 0.0_dp, dp)
      end if
    end do

    call closeFile(fd)

  end subroutine writeDM


  !> Mulliken population for a square density matrix and overlap
  !> Note: assumes both triangles of both square matrices are filled
  subroutine getExcMulliken(denseDesc, pc, s, dqex)

    !> Dense matrix descriptor
    type(TDenseDescr), intent(in) :: denseDesc
    
    !> density matrix
    real(dp), intent(in) :: pc(:,:)

    !> Overlap matrix
    real(dp), intent(in) :: s(:,:)

    !> Output atomic charges
    real(dp), intent(out) :: dqex(:)

    real(dp) :: tmp(size(pc,dim=1))
    integer :: iAt1

    @:ASSERT(all(shape(pc)==shape(s)))

    tmp = sum(pc * s,dim=2)
    dqex(:) = 0.0_dp
    do iAt1 = 1, size(dqex)
      dqex(iAt1) = sum(tmp(denseDesc%iAtomStart(iAt1):denseDesc%iAtomStart(iAt1 + 1) -1))
    end do

  end subroutine getExcMulliken


  !> Calculation of force from derivatives of excitation energy
  !! 1. we need the ground and excited Mulliken charges
  !! 2. we need P,(T,Z),W, X + Y from linear response
  !! 3. calculate dsmndr, dhmndr (dS/dR, dh/dR), dgabda (dGamma_{IAt1,IAt2}/dR_{IAt1}),
  !! dgext (dGamma-EXT_{IAt1,k}/dR_{IAt1})
  subroutine addGradients(sym, nxov, natom, species0, env, denseDesc, norb, homo, getIA,&
      & getIJ, getAB, win, grndEigVecs, pc, ovrXev, dq_ud, dqex, gammaMat, lrGamma, HubbardU, &
      & spinW, shift, woo, wov, wvv, transChrg, xpy, xmy, coord0, orb, skHamCont, skOverCont, &
      & derivator, rhoSqr, tHybridXc, hybridXc, excgrad, deltaRho)

    !> Symmetry of the transition
    character, intent(in) :: sym

    !> Number of single particle transitions to include
    integer, intent(in) :: nxov

    !> Number of central cell atoms
    integer, intent(in) :: natom

    !> Central cell chemical species
    integer, intent(in) :: species0(:)

    !> Environment settings
    type(TEnvironment), intent(inout) :: env

    !> Dense matrix descriptor
    type(TDenseDescr), intent(in) :: denseDesc
    
    !> number of orbitals for ground state system
    integer, intent(in) :: norb

    !> Number of highest occupied state in ground state
    integer, intent(in) :: homo(:)

    !> Index array from composite occ-vir transition index to specific single particle states
    integer, intent(in) :: getIA(:,:)

    !> Index array from composite occ-occ transition index to specific single particle states
    integer, intent(in) :: getIJ(:,:)

    !> Index array from composite vir-vir transition index to specific single particle states
    integer, intent(in) :: getAB(:,:)

    !> Single particle transition index
    integer, intent(in) :: win(:)

    !> Ground state eigenvectors
    real(dp), intent(in) :: grndEigVecs(:,:,:)

    !> Transition density matrix
    real(dp), intent(in) :: pc(:,:,:)

    !> Overlap times ground state eigenvectors
    real(dp), intent(in) :: ovrXev(:,:,:)

    !> Ground state gross charges
    real(dp), intent(in) :: dq_ud(:,:)

    !> Charge differences from ground to excited state
    real(dp), intent(in) :: dqex(:,:)

    !> Softened coulomb matrix
    real(dp), intent(in) :: gammaMat(:,:)

    !> Electrostatic matrix, long-range corrected
    real(dp), allocatable, intent(in) :: lrGamma(:,:)

    !> Ground state Hubbard U values
    real(dp), intent(in) :: HubbardU(:)

    !> Ground state spin derivatives for each species
    real(dp), intent(in) :: spinW(:)

    !> Ground state potentials (shift vector)
    real(dp), intent(in) :: shift(:)

    !> W vector occupied part
    real(dp), intent(in) :: woo(:,:)

    !> W vector occupied-virtual part
    real(dp), intent(in) :: wov(:)

    !> W vector virtual part
    real(dp), intent(in) :: wvv(:,:)

    !> Machinery for transition charges between single particle levels
    type(TTransCharges), intent(in) :: transChrg

    !> X+Y Furche term
    real(dp), intent(in) :: xpy(:)

    !> X-Y Furche term
    real(dp), intent(in) :: xmy(:)

    !> Central cell atomic coordinates
    real(dp), intent(in) :: coord0(:,:)

    !> Data type for atomic orbital information
    type(TOrbitals), intent(in) :: orb

    !> H0 data
    type(TSlakoCont), intent(in) :: skHamCont

    !> Overlap data
    type(TSlakoCont), intent(in) :: skOverCont

    !> Differentiator for the non-scc matrices
    class(TNonSccDiff), intent(in) :: derivator

    !> Ground state density matrix
    real(dp), intent(in) :: rhoSqr(:,:,:)

    !> Is calculation range-separated?
    logical, intent(in) :: tHybridXc

    !> Data for range-separated calculation
    class(THybridXcFunc), allocatable, intent(inout) :: hybridXc

    !> Resulting excited state gradient
    real(dp), intent(out) :: excgrad(:,:)

    !> Difference density matrix (vs. uncharged atoms)
    real(dp), intent(inout), optional :: deltaRho(:,:,:)


    real(dp), allocatable :: shift_excited(:,:), xpyq(:), xpyqds(:)
    real(dp), allocatable :: shxpyq(:,:), xpycc(:,:,:), wcc(:,:,:), tmp5(:), tmp7(:), tmp11(:)
    real(dp), allocatable :: qTr(:), temp(:), dq(:), dm(:), dsigma(:)
    real(dp), allocatable :: dH0(:,:,:), dSo(:,:,:)
    real(dp), allocatable :: Dens(:,:), SpinDens(:,:)
    real(dp), allocatable :: xmycc(:,:,:), xpyas(:,:,:), xmyas(:,:,:)
    real(dp), allocatable :: overlap(:,:), lrGammaOrb(:,:), gammaLongRangePrime(:,:,:)
    real(dp), allocatable :: PS(:,:,:), DS(:,:,:), SPS(:,:,:), SDS(:,:,:), SX(:,:,:)
    real(dp), allocatable :: XS(:,:,:), SXS(:,:,:), SY(:,:,:), YS(:,:,:), SYS(:,:,:)
    integer :: ia, i, j, a, b, ab, ij, m, n, mu, nu, xyz, iAt1, iAt2, ka
    integer :: indalpha, indalpha1, indbeta, indbeta1, soo(2), svv(2)
    integer :: iSp1, iSp2, iSpin, nSpin
    real(dp) :: tmp1, tmp2, tmp3, tmp4, tmp6, tmp8, tmp9, tmp10, rab
    real(dp) :: diffvec(3), dgab(3), tmpVec(3), tmp3a, tmp3b, tmprs, tmprs2, tmps(2)
    integer, allocatable :: nxoo(:), nxvv(:), nvir(:), species(:)
    logical :: tSpin

    nSpin = size(grndEigVecs, dim=3)
    tSpin = (nSpin == 2)

    allocate(shift_excited(natom, nSpin))
    allocate(xpyq(natom))
    allocate(shxpyq(natom, nSpin))
    allocate(xpycc(norb, norb, nSpin))
    allocate(wcc(norb, norb, nSpin))
    allocate(qTr(natom))
    allocate(temp(norb))
    allocate(tmp5(nSpin))
    allocate(tmp7(nSpin))

    allocate(Dens(norb, norb))
    Dens(:,:) = sum(rhoSqr, dim=3)

    allocate(dH0(orb%mOrb, orb%mOrb, 3))
    allocate(dSo(orb%mOrb, orb%mOrb, 3))

    allocate(nxoo(nSpin))
    allocate(nxvv(nSpin))
    allocate(nvir(nSpin))

    nxoo(:) = (homo(:)*(homo(:)+1))/2
    nvir(:) = norb - homo(:)
    nxvv(:) = (nvir(:)*(nvir(:)+1))/2

    soo(:) = [0, nxoo(1)]
    svv(:) = [0, nxvv(1)]

    allocate(dq(natom))
    dq(:) = dq_ud(:,1)

    if (tSpin) then
      allocate(dm(natom))
      allocate(xpyqds(natom))
      allocate(tmp11(nSpin))

      allocate(SpinDens(norb,norb))
      SpinDens(:,:) = rhoSqr(:,:,1) - rhoSqr(:,:,2)

      allocate(dsigma(2))
      dsigma(1) = 1.0_dp
      dsigma(2) = -1.0_dp
      dm(:) = dq_ud(:,2)
    end if

    if (tHybridXc) then
      allocate(xmycc(norb, norb, nSpin))
      allocate(xpyas(norb, norb, nSpin))
      allocate(xmyas(norb, norb, nSpin))
      allocate(PS(norb, norb, nSpin))
      allocate(DS(norb, norb, nSpin))
      allocate(SPS(norb, norb, nSpin))
      allocate(SDS(norb, norb, nSpin))
      allocate(SX(norb, norb, nSpin))
      allocate(XS(norb, norb, nSpin))
      allocate(SXS(norb, norb, nSpin))
      allocate(SY(norb, norb, nSpin))
      allocate(YS(norb, norb, nSpin))
      allocate(SYS(norb, norb, nSpin))
      allocate(overlap(norb, norb))
      allocate(lrGammaOrb(norb, norb))
      allocate(gammaLongRangePrime(3, nAtom, nAtom))

      ! Symmetrize deltaRho
      do mu = 1, norb
        do nu = mu + 1, norb
          deltaRho(mu,nu,:) = deltaRho(nu,mu,:)
        end do
      end do

      ! Compute long-range gamma derivative
      gammaLongRangePrime(:,:,:) = 0._dp
      call hybridXc%getCentralCellSpecies(species)
      do iAt1 = 1, nAtom
        do iAt2 = 1, nAtom
          if (iAt1 /= iAt2) then
            call getDirectionalCamGammaPrimeValue(hybridXc, tmpVec, iAt1, iAt2)
            gammaLongRangePrime(:, iAt1, iAt2) = tmpVec
          end if
        end do
      end do

      ! Symmetrize S (can't we get S from caller?)
      call getSqrS(coord0, nAtom, skOverCont, orb, denseDesc%iAtomStart, species0, overlap)
      call getSqrGamma(nAtom, lrGamma, denseDesc%iAtomStart, lrGammaOrb)

    end if

    excgrad = 0.0_dp

    ! excited state potentials at atomic sites
    do iSpin = 1, nSpin
      call hemv(shift_excited(:,iSpin), gammaMat, dqex(:,iSpin))
    end do

    ! xypq(alpha) = sum_ia (X+Y)_ia q^ia(alpha)
    ! complexity norb * norb * norb
    xpyq(:) = 0.0_dp
    call transChrg%qMatVec(env, denseDesc, ovrXev, grndEigVecs, getIA, win, xpy, xpyq)

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
      call transChrg%qMatVecDs(env, denseDesc, ovrXev, grndEigVecs, getIA, win, xpy, xpyqds)
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

    if (tHybridXc) then

      ! Asymmetric contribution: xmycc_as = sum_ias (X-Y)_ias c_mas c_nis
      xmycc(:,:,:) = 0.0_dp
      xpyas(:,:,:) = 0.0_dp
      xmyas(:,:,:) = 0.0_dp
      do ia = 1, nxov
        call indxov(win, ia, getIA, i, a, iSpin)
        ! should replace with DSYR2 call:
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
      indalpha = denseDesc%iAtomStart(iAt1)
      indalpha1 = denseDesc%iAtomStart(iAt1 + 1) -1
      iSp1 = species0(iAt1)

      do iAt2 = 1, iAt1 - 1
        indbeta = denseDesc%iAtomStart(iAt2)
        indbeta1 = denseDesc%iAtomStart(iAt2 + 1) -1
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

        if (tHybridXc) then
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

                if (tHybridXc) then
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


  !> Write out excitations projected onto ground state.
  subroutine writeCoeffs(tt, grndEigVecs, occ, tCoeffs, tIncGroundState, occNatural, naturalOrbs)

    !> T part of the matrix
    real(dp), intent(in) :: tt(:,:,:)

    !> Ground state eigenvectors
    real(dp), intent(in) :: grndEigVecs(:,:,:)

    !> Ground state occupations
    real(dp), intent(in) :: occ(:,:)

    !> Save the coefficients of the natural orbitals
    logical, intent(in) :: tCoeffs

    !> Include the ground state as well as the transition part
    logical, intent(in) :: tIncGroundState

    !> Natural orbital occupation numbers
    real(dp), intent(out), optional :: occNatural(:)

    !> Natural orbitals
    real(dp), intent(out), optional :: naturalOrbs(:,:,:)

    real(dp), allocatable :: t2(:,:,:), occtmp(:,:)
    integer :: norb, nSpin, ii, jj, mm, iSpin
    logical :: tSpin

    type(TFileDescr) :: fdCoeffs

    norb = size(tt, dim=1)
    nSpin = size(tt, dim=3)
    tSpin = (nSpin == 2)

    if (present(occNatural).or.tCoeffs) then

      allocate(t2(norb, norb, nSpin))
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
          allocate(occtmp(size(occ), nSpin))
          occTmp(:,1) = occNatural
        end if
      else
        allocate(occtmp(size(occ), nSpin))
        occtmp = 0.0_dp
        do iSpin = 1, nSpin
          call evalCoeffs(t2(:,:,iSpin), occtmp(:,iSpin), grndEigVecs(:,:,iSpin))
        end do
      end if

      ! Better to get this by post-processing DFTB+ output, but here for
      ! compatibility at the moment
      if (tCoeffs) then
        call openFile(fdCoeffs, excitedCoefsOut, mode="a")
        write(fdCoeffs%unit,*) 'T F'
        if (.not. tSpin) then
          do ii = 1, norb
            jj = norb - ii + 1
            write(fdCoeffs%unit, '(1x,i3,1x,f13.10,1x,f13.10)') ii, occtmp(jj,1), 2.0_dp
            write(fdCoeffs%unit, '(6(f13.10,1x))') (cmplx(t2(mm,jj,1), kind=dp),&
                & mm = 1, norb)
          end do
        else
          do iSpin = 1, nSpin
            write(fdCoeffs%unit,*)
            write(fdCoeffs%unit, '(1x,a,1x,i1)') 'SPIN', iSpin
            do ii = 1, norb
              jj = norb - ii + 1
              write(fdCoeffs%unit, '(1x,i3,1x,f13.10,1x,f13.10)') ii, occtmp(jj,iSpin), 1.0_dp
              write(fdCoeffs%unit, '(6(f13.10,1x))') (cmplx(t2(mm,jj,iSpin), kind=dp),&
                  & mm = 1, norb)
            end do
          end do
        end if

        call closeFile(fdCoeffs)
      end if

    end if

  end subroutine writeCoeffs


  !> Project MO density matrix onto ground state orbitals.
  subroutine evalCoeffs(t2, occ, eig)

    !> Density matrix
    real(dp), intent(inout) :: t2(:,:)

    !> Resulting natural orbital occupations
    real(dp), intent(out) :: occ(:)

    !> 'natural' eigenvectors
    real(dp), intent(in) :: eig(:,:)

    real(dp), allocatable :: coeffs(:,:)

    allocate(coeffs(size(occ),size(occ)))

    call heev(t2, occ, 'U', 'V')
    call gemm(coeffs, eig, t2)
    t2 = coeffs

  end subroutine evalCoeffs


  !> Write out transitions from ground to excited state along with single particle transitions and
  !! dipole strengths.
  subroutine writeExcitations(sym, osz, nexc, nmatup, getIA, win, eval, xpy, wij, fdXPlusY,&
      & fdTrans, fdTransDip, transitionDipoles, fdTagged, taggedWriter, fdExc, Ssq)

    !> Symmetry label for the type of transition
    character, intent(in) :: sym

    !> Oscillator strengths for transitions from ground to excited states
    real(dp), intent(in) :: osz(:)

    !> Number of excited states to solve for
    integer, intent(in) :: nexc

    !> Number of same spin excitations
    integer, intent(in) :: nmatup

    !> Index array between transitions in square and 1D representations
    integer, intent(in) :: getIA(:,:)

    !> Index array for single particle excitations
    integer, intent(in) :: win(:)

    !> Excitation energies
    real(dp), intent(in) :: eval(:)

    !> Eigenvectors of excited states (X+Y)
    real(dp), intent(in) :: xpy(:,:)

    !> Single particle excitation energies
    real(dp), intent(in) :: wij(:)

    !> Single particle transition dipole moments
    real(dp), intent(in) :: transitionDipoles(:,:)

    !> File unit for transition dipoles
    type(TFileDescr), intent(in) :: fdTransDip

    !> File unit for X+Y data
    type(TFileDescr), intent(in) :: fdXPlusY

    !> File unit for transitions
    type(TFileDescr), intent(in) :: fdTrans

    !> File unit for tagged output (> -1 for write out)
    type(TFileDescr), intent(in) :: fdTagged

    !> Tagged writer
    type(TTaggedWriter), intent(inout) :: taggedWriter

    !> File unit for excitation energies
    type(TFileDescr), intent(in) :: fdExc

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

    if (fdXplusY%isConnected()) then
      write(fdXPlusY%unit, *) nmat, nexc
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
          write(fdExc%unit,&
              & '(1x,f10.3,4x,f14.8,2x,i5,3x,a,1x,i5,7x,f6.3,2x,f10.3,4x,&
              & f6.3)')&
              & Hartree__eV * sqrt(eval(ii)), osz(ii), m, '->', n, weight,&
              & Hartree__eV * wij(iWeight), Ssq(ii)
        else
          write(fdExc%unit,&
              & '(1x,f10.3,4x,f14.8,5x,i5,3x,a,1x,i5,7x,f6.3,2x,f10.3,6x,a)')&
              & Hartree__eV * sqrt(eval(ii)), osz(ii), m, '->', n, weight,&
              & Hartree__eV * wij(iWeight), sign
        end if

        if (fdXplusY%isConnected()) then
          if (tSpin) then
            updwn = (win(iweight) <= nmatup)
            sign = "D"
            if (updwn) sign = "U"
          end if
          write(fdXPlusY%unit, '(1x,i5,3x,a,3x,ES17.10)') ii, sign, sqrt(eval(ii))
          write(fdXPlusY%unit, '(6(1x,ES17.10))') xpy(:,ii)
        endif

        if (fdTrans%isConnected()) then
          write(fdTrans%unit, '(2x,a,T12,i5,T21,ES17.10,1x,a,2x,a)')&
              & 'Energy ', ii,  Hartree__eV * sqrt(eval(ii)), 'eV', sign
          write(fdTrans%unit,*)
          write(fdTrans%unit,'(2x,a,9x,a,8x,a)')'Transition', 'Weight', 'KS [eV]'
          write(fdTrans%unit,'(1x,45("="))')

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
            write(fdTrans%unit, '(i5,3x,a,1x,i5,1x,1a,T22,f10.8,T33,f14.8)')&
                & m, '->', n, sign, wvec(jj), Hartree__eV * wij(wvin(jj))
          end do
          write(fdTrans%unit,*)
        end if

        if (fdTransDip%isConnected()) then
          write(fdTransDip%unit, '(1x,i5,1x,f10.3,2x,3(ES14.6))')&
              & ii, Hartree__eV * sqrt(eval(ii)), (transitionDipoles(ii,jj)&
              & * au__Debye, jj=1,3)
        end if

      else

        ! find largest coefficient in CI - should use maxloc
        call index_heap_sort(wvin, wvec)
        wvin = wvin(size(wvin):1:-1)
        wvec = wvec(wvin)

        weight = wvec(1)
        iweight = wvin(1)
        call indxov(win, iWeight, getIA, m, n, s)
        sign = sym

        if (tSpin) then
          sign = " "
          write(fdExc%unit,&
              & '(6x,A,T12,4x,f14.8,2x,i5,3x,a,1x,i5,7x,A,2x,f10.3,4x,f6.3)')&
              & '< 0', osz(ii), m, '->', n, '-', Hartree__eV * wij(iWeight),&
              & Ssq(ii)
        else
          write(fdExc%unit,&
              & '(6x,A,T12,4x,f14.8,2x,i5,3x,a,1x,i5,7x,f6.3,2x,f10.3,6x,a)')&
              & '< 0', osz(ii), m, '->', n, weight, Hartree__eV * wij(iWeight), sign
        end if

        if (fdXplusY%isConnected()) then
          if (tSpin) then
            updwn = (win(iweight) <= nmatup)
            sign = "D"
            if (updwn) sign = "U"
          end if
          write(fdXPlusY%unit, '(1x,i5,3x,a,3x,A)') ii,sign, '-'
        endif

        if (fdTrans%isConnected()) then
          write(fdTrans%unit, '(2x,a,1x,i5,5x,a,1x,a,3x,a)') 'Energy ', ii,  '-', 'eV', sign
          write(fdTrans%unit,*)
        end if

        if (fdTransDip%isConnected()) then
          write(fdTransDip%unit, '(1x,i5,1x,A)') ii, '-'
        endif

      end if

    end do

    deallocate(wvec)
    deallocate(wvin)

    if (fdTagged%isConnected()) then

      call degeneracyFind%init(elecTolMax)
      call degeneracyFind%degeneracyTest(eval, tDegenerate)
      if (.not.tDegenerate) then
        call taggedWriter%write(fdTagged%unit, tagLabels%excEgy, eval)
        call taggedWriter%write(fdTagged%unit, tagLabels%excOsc, osz)
        ! Since the transition dipole file exists, transition dipoles had been calculated
        if (fdTransDip%isConnected()) then
          call taggedWriter%write(fdTagged%unit, tagLabels%excDipole,&
              & sqrt(sum(transitionDipoles**2,dim=2)))
        end if
      else
        degenerate = DegeneracyFind%degenerateRanges()
        call taggedWriter%write(fdTagged%unit, tagLabels%excEgy, eval(degenerate(1,:)))
        ! sum oscillator strength over any degenerate levels
        allocate(oDeg(DegeneracyFind%degenerateGroups()))
        do ii = 1, size(oDeg)
          oDeg(ii) = sum(osz(degenerate(1,ii):degenerate(2,ii)))
        end do
        call taggedWriter%write(fdTagged%unit, tagLabels%excOsc, oDeg)
        ! Since the transition dipole file exists, transition dipoles had been calculated
        if (fdTransDip%isConnected()) then
          oDeg(:) = 0.0_dp
          do ii = 1, size(oDeg)
            oDeg(ii) = sqrt(sum(transitionDipoles(degenerate(1,ii):degenerate(2,ii),:)**2))
          end do
          call taggedWriter%write(fdTagged%unit, tagLabels%excDipole, oDeg)
        end if
      end if
    end if

  end subroutine writeExcitations


  !> Create transition density matrix in MO basis P = T + 1/2 Z symmetric (paper has T + Z
  !! asymmetric) (Zab = Zij = 0, Tia = 0).
  subroutine calcPMatrix(t, rhs, win, getIA, pc)

    !> T matrix
    real(dp), intent(in) :: t(:,:,:)

    !> Z matrix
    real(dp), intent(in) :: rhs(:)

    !> Index array for single particle transitions
    integer, intent(in) :: win(:)

    !> Array of the occupied->virtual pairs (nTransitions,occ 1 or virtual 2)
    integer, intent(in) :: getIA(:,:)

    !> Resulting excited state density matrix
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


  !> Computes H^+/-_pq [V] as defined in Furche JCP 117 7433 (2002) eq. 20
  !! Here p/q are virtual orbitals and V is either X+Y or X-Y
  subroutine getHvvXY(ipm, nXvv, homo, nAtom, iaTrans, getIA, getAB, win,&
      & env, denseDesc, ovrXev, grndEigVecs, lrGamma, transChrg, XorY, vecHvv)

    !> Sign s of H in H^(s)[V]
    integer, intent(in) :: ipm

    !> Number of vir-vir transitions per spin channel
    integer, intent(in) :: nXvv(:)

    !> Occupied orbitals per spin channel
    integer, intent(in) :: homo(:)

    !> Number of atoms
    integer, intent(in) :: nAtom

    !> Array from pairs of single particles states to compound index
    integer, intent(in) :: iaTrans(:,:,:)

    !> Index array for occ-vir single particle excitations
    integer, intent(in) :: getIA(:,:)

    !> Index array for vir-vir single particle excitations
    integer, intent(in) :: getAB(:,:)

    !> Index array for single particle excitations
    integer, intent(in) :: win(:)

    !> Environment settings
    type(TEnvironment), intent(inout) :: env

    !> Dense matrix descriptor
    type(TDenseDescr), intent(in) :: denseDesc
    
    !> overlap times ground state eigenvectors
    real(dp), intent(in) :: ovrXev(:,:,:)

    !> Ground state eigenvectors
    real(dp), intent(in) :: grndEigVecs(:,:,:)

    !> Electrostatic matrix, long-range corrected
    real(dp), allocatable, intent(in) :: lrGamma(:,:)

    !> Machinery for transition charges between single particle levels
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
        qIJ = transChrg%qTransAB(abs, env, denseDesc, ovrXev, grndEigVecs, getAB)
        qX(:,ias) = qX(:,ias) + qIJ * XorY(ibs)
      end do
    end do

    Gq(:,:) = 0.0_dp
    do ias = 1, nXov
      qIJ = transChrg%qTransIA(ias, env, denseDesc, ovrXev, grndEigVecs, getIA, win)
      call hemv(gqIJ, lrGamma, qIJ, uplo='U')
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


  !> Computes H^+/-_pq [V] as defined in Furche JCP 117 7433 (2002) eq. 20
  !! Here p/q are occupied orbitals and V is either X+Y or X-Y
  subroutine getHooXY(ipm, nXoo, homo, nAtom, iaTrans, getIA, getIJ, win,&
      & env, denseDesc, ovrXev, grndEigVecs, lrGamma, transChrg, XorY, vecHoo)

    !> Sign s of H in H^(s)[V]
    integer, intent(in) :: ipm

    !> Number of occ-occ transitions per spin channel
    integer, intent(in) :: nXoo(:)

    !> Occupied orbitals per spin channel
    integer, intent(in) :: homo(:)

    !> Number of atoms
    integer, intent(in) :: nAtom

    !> Array from pairs of single particles states to compound index
    integer, intent(in) :: iaTrans(:,:,:)

    !> Index array for occ-vir single particle excitations
    integer, intent(in) :: getIA(:,:)

    !> Index array for occ-occ single particle excitations
    integer, intent(in) :: getIJ(:,:)

    !> Index array for single particle excitations
    integer, intent(in) :: win(:)

    !> Environment settings
    type(TEnvironment), intent(inout) :: env

    !> Dense matrix descriptor
    type(TDenseDescr), intent(in) :: denseDesc
    
    !> overlap times ground state eigenvectors
    real(dp), intent(in) :: ovrXev(:,:,:)

    !> Ground state eigenvectors
    real(dp), intent(in) :: grndEigVecs(:,:,:)

    !> Electrostatic matrix, long-range corrected
    real(dp), allocatable, intent(in) :: lrGamma(:,:)

    !> Machinery for transition charges between single particle levels
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
        qIJ = transChrg%qTransIJ(ijs, env, denseDesc, ovrXev, grndEigVecs, getIJ)
        qX(:,ias) = qX(:,ias) + qIJ * XorY(jas)
      end do
    end do

    Gq(:,:) = 0.0_dp
    do ias = 1, nXov
      call indXov(win, ias, getIA, i, a, s)
      qIJ = transChrg%qTransIA(ias, env, denseDesc, ovrXev, grndEigVecs, getIA, win)
      call hemv(gqIJ, lrGamma, qIJ, uplo='U')
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


  !> Computes H^+/-_pq [T] as defined in Furche JCP 117 7433 (2002) eq. 20
  !! Here p is an occupied MO and q is a virtual one, T is the relaxed difference density
  subroutine getHovT(nXoo, nXvv, homo, nAtom, iaTrans, getIA, getIJ, getAB, win,&
      & env, denseDesc, ovrXev, grndEigVecs, lrGamma, transChrg, t, vecHovT)

    !> Number of occ-occ transitions per spin channel
    integer, intent(in) :: nXoo(:)

    !> Number of vir-vir transitions per spin channel
    integer, intent(in) :: nXvv(:)

    !> Occupied orbitals per spin channel
    integer, intent(in) :: homo(:)

    !> Number of atoms
    integer, intent(in) :: nAtom

    !> Array from pairs of single particles states to compound index
    integer, intent(in) :: iaTrans(:,:,:)

    !> Index array for occ-vir single particle excitations
    integer, intent(in) :: getIA(:,:)

    !> Index array for occ-occ single particle excitations
    integer, intent(in) :: getIJ(:,:)

    !> Index array for vir-vir single particle excitations
    integer, intent(in) :: getAB(:,:)

    !> Index array for single particle excitations
    integer, intent(in) :: win(:)

    !> Environment settings
    type(TEnvironment), intent(inout) :: env

    !> Dense matrix descriptor
    type(TDenseDescr), intent(in) :: denseDesc 

    !> Overlap times ground state eigenvectors
    real(dp), intent(in) :: ovrXev(:,:,:)

    !> Ground state eigenvectors
    real(dp), intent(in) :: grndEigVecs(:,:,:)

    !> Electrostatic matrix, long-range corrected
    real(dp), allocatable, intent(in) :: lrGamma(:,:)

    !> Machinery for transition charges between single particle levels
    type(TTransCharges), intent(in) :: transChrg

    !> Excited state density matrix
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
        qIJ = transChrg%qTransIA(ibs, env, denseDesc, ovrXev, grndEigVecs, getIA, win)
        qX(:,ias) = qX(:,ias) + qIJ * t(a,b,s)
      end do
    end do

    Gq(:,:) = 0.0_dp
    do abs = 1, sum(nXvv)
      qIJ = transChrg%qTransAB(abs, env, denseDesc, ovrXev, grndEigVecs, getAB)
      call hemv(gqIJ, lrGamma, qIJ, uplo='U')
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
        qIJ = transChrg%qTransIA(jas, env, denseDesc, ovrXev, grndEigVecs, getIA, win)
        qX(:,ias) = qX(:,ias) + qIJ * t(i,j,s)
      end do
    end do

    Gq(:,:) = 0.0_dp
    do ijs = 1, sum(nXoo)
      i = getIJ(ijs, 1)
      j = getIJ(ijs, 2)
      s = getIJ(ijs, 3)
      qIJ = transChrg%qTransIJ(ijs, env, denseDesc, ovrXev, grndEigVecs, getIJ)
      call hemv(gqIJ, lrGamma, qIJ, uplo='U')
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


  !> Computes H^+/-_pq [T] as defined in Furche JCP 117 7433 (2002) eq. 20
  !! Here p/q are occupied MO, T is the relaxed difference density
  subroutine getHooT(nXov, nXoo, homo, nAtom, iaTrans, getIA, getIJ, win,&
      & env, denseDesc, ovrXev, grndEigVecs, lrGamma, transChrg, t, vecHooT)

    !> Number of occ-vir transitions per spin channel
    integer, intent(in) :: nXov

    !> Number of occ-occ transitions per spin channel
    integer, intent(in) :: nXoo(:)

    !> Occupied orbitals per spin channel
    integer, intent(in) :: homo(:)

    !> Number of atoms
    integer, intent(in) :: nAtom

    !> Array from pairs of single particles states to compound index
    integer, intent(in) :: iaTrans(:,:,:)

    !> Index array for occ-vir single particle excitations
    integer, intent(in) :: getIA(:,:)

    !> Index array for occ-occ single particle excitations
    integer, intent(in) :: getIJ(:,:)

    !> Index array for single particle excitations
    integer, intent(in) :: win(:)

    !> Environment settings
    type(TEnvironment), intent(inout) :: env

    !> Dense matrix descriptor
    type(TDenseDescr), intent(in) :: denseDesc 

    !> Overlap times ground state eigenvectors
    real(dp), intent(in) :: ovrXev(:,:,:)

    !> Ground state eigenvectors
    real(dp), intent(in) :: grndEigVecs(:,:,:)

    !> Electrostatic matrix, long-range corrected
    real(dp), allocatable, intent(in) :: lrGamma(:,:)

    !> Machinery for transition charges between single particle levels
    type(TTransCharges), intent(in) :: transChrg

    !> Excited state density matrix
    real(dp), intent(in) :: t(:,:,:)

    !> Output vector H[T] occ-occ
    real(dp), intent(out) :: vecHooT(:)

    real(dp), allocatable :: qIJ(:), gqIJ(:), qX(:,:), Gq(:,:), qXa(:,:,:)
    integer :: nOrb, iSpin, nSpin, iMx, soo(2)
    integer :: i, j, k, a, b, s, ij, ias, ibs, ijs, jas, iks, jks

    nOrb = size(ovrXev, dim=1)
    nSpin = size(t, dim=3)
    soo(:) = [0, nXoo(1)]

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
        qIJ = transChrg%qTransIA(ibs, env, denseDesc, ovrXev, grndEigVecs, getIA, win)
        qX(:,ias) = qX(:,ias) + qIJ * t(a,b,s)
      end do
    end do

    Gq(:,:) = 0.0_dp
    do ias = 1, nXov
      qIJ = transChrg%qTransIA(ias, env, denseDesc, ovrXev, grndEigVecs, getIA, win)
      call hemv(gqIJ, lrGamma, qIJ, uplo='U')
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

    Gq(:,:) = 0.0_dp
    do ijs = 1, sum(nXoo)
      qIJ = transChrg%qTransIJ(ijs, env, denseDesc, ovrXev, grndEigVecs, getIJ)
      call hemv(gqIJ, lrGamma, qIJ, uplo='U')
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
          qIJ = transChrg%qTransIJ(iks, env, denseDesc, ovrXev, grndEigVecs, getIJ)
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


  !> Constructs the full overlap matrix S.
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
      ! Diagonal entries
      S(mu, mu) = 1.0_dp
    end do

  end subroutine getSqrS


  !> Constructs a Gamma-Matrix of dimension nOrb instead of nAtoms.
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


  !> Helper routine to construct overlap.
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
    vect(:) = vect / dist
    call getSKIntegrals(skOverCont, interSKOver, dist, iSp1, iSp2)
    call rotateH0(Sblock, interSKOver, vect(1), vect(2), vect(3), iSp1, iSp2, orb)

  end subroutine getSOffsite

  !> Compute (fake) transition density matrix and W for ground-to-excited state couplings
  !! See TCA 140 34 (2020) and JCP 132 044107 (2010)
  !! Actually omega * W is computed
  !! TODO: Spin-polarized systems 
  subroutine grndToExcDensityMatrices(tHybridXc, xpy, xmy, win, env, denseDesc, homo, transChrg,& 
      & getIA, getIJ, getAB, iaTrans, nAtom, species0, grndEigVal, ovrXev, grndEigVecs, frGamma,&
      & lrGamma, spinW, omega, sym, pc, wov, woo)

    !> Is calculation range-separated?
    logical, intent(in) :: tHybridXc

    !> X+Y Furche term for excited state n
    real(dp), intent(in) :: xpy(:)

    !> X-Y Furche term for excited state n
    real(dp), intent(in) :: xmy(:)

    !> Index array for single particle transitions
    integer, intent(in) :: win(:)
    
    !> Environment settings
    type(TEnvironment), intent(inout) :: env

    !! Dense matrix descriptor
    type(TDenseDescr), intent(in) :: denseDesc 

    !> Highest occupied level
    integer, intent(in) :: homo(:)

    !> Machinery for transition charges between single particle levels
    type(TTransCharges), intent(in) :: transChrg

    !> Index array between transitions in square and 1D representations
    integer, intent(in) :: getIA(:,:)

    !> Index array for vir-vir transitions
    integer, intent(in) :: getIJ(:,:)

    !> Index array for occ-occ transitions
    integer, intent(in) :: getAB(:,:)

    !> Index array from orbital pairs to compound index
    integer, intent(in) :: iaTrans(:,:,:)

    !> Number of central cell atoms
    integer, intent(in) :: nAtom

    !> Central cell chemical species
    integer, intent(in) :: species0(:)

    !> Ground state wavefunctions
    real(dp), intent(in) :: grndEigVal(:,:)

    !> Overlap times ground state wavefunctions
    real(dp), intent(in) :: ovrXev(:,:,:)

    !> Ground state wavefunctions
    real(dp), intent(in) :: grndEigVecs(:,:,:)

    !> Softened coulomb matrix
    real(dp), intent(in) :: frGamma(:,:)

    !> Softened coulomb matrix, long-range corrected
    real(dp), allocatable, intent(in) :: lrGamma(:,:)

    !> Ground state spin derivatives for each species
    real(dp), intent(in) :: spinW(:)

    !> Excitation energy of states n
    real(dp), intent(in) :: omega

    !> Symmetry of the transitions
    character, intent(in) :: sym

    !> P matrix (symmetric)
    real(dp), intent(out) :: pc(:,:,:)

    !> W^+ vector occupied-virtual part
    real(dp), intent(out) :: wov(:)

    !> W^+ vector occupied part
    real(dp), intent(out) :: woo(:,:)

    real(dp), allocatable :: p(:), vecHoo(:)
    integer, allocatable :: nxoo(:), nxvv(:), nvir(:)
    integer :: soo(2), i, a, s, ias, j, ij, ijs, nxov, nSpin, nOrb

    nxov = size(wov)
    nOrb = size(ovrXev, dim=1)
    nSpin = size(homo)

    allocate(p(nxov))
    allocate(nxoo(nSpin))
    allocate(nxvv(nSpin))
    allocate(nvir(nSpin))

    nxoo(:) = (homo * (homo + 1)) / 2
    nvir(:) = nOrb - homo
    nxvv(:) = (nvir * (nvir + 1)) / 2

    allocate(vecHoo(sum(nxoo)))

    ! Transition charges use compound index ijs = ij + soo(s)
    soo(:) = [0, nxoo(1)]

    p(:) = 0.0_dp
    wov(:) = 0.0_dp
    woo(:,:) = 0.0_dp

    ! "Fake" density matrix for non-adiabatic coupling [Furche JCP 132 044107 (2010)]
    ! Restricted KS: P = 2 P^up ; (X+Y) = sqrt(2) (X+Y)^up
    do ias = 1, nxov
      call indxov(win, ias, getIA, i, a, s)
      p(ias) = sqrt(2.0_dp) * xpy(ias)
      wov(ias) = grndEigVal(i, s) * p(ias) + omega * xmy(ias) / sqrt(2.0_dp)
    end do

    ! Define P symmetrically (similar to treatment of excited state gradients)
    pc(:,:,:) = 0.0_dp
    do ias = 1, nxov
      call indxov(win, ias, getIA, i, a, s)
      pc(i,a,s) = 0.5_dp * p(ias)
      pc(a,i,s) = 0.5_dp * p(ias)
    end do
    call getHplusXYfr(sym, nXoo, nXvv, nAtom, getIA, getIJ, getAB, win, env, denseDesc,&
      & species0, ovrXev, grndEigVecs, frGamma, spinW, transChrg, p, vecHoo=vecHoo)

    do s = 1, nSpin
      do ij = 1, nXoo(s)
        i = getIJ(ij + soo(s), 1)
        j = getIJ(ij + soo(s), 2)
        ijs = iaTrans(i, j, s)
        !> GetHplusXYfr used with P instead of X+Y yields half the desired result
        woo(ij,s) = 2.0_dp * vecHoo(ijs)
      end do
    end do

    if (tHybridXc) then
      call getHooXY(1, nXoo, homo, nAtom, iaTrans, getIA, getIJ, win, env, denseDesc,&
       & ovrXev, grndEigVecs, lrGamma, transChrg, p, vecHoo)
      do s = 1, nSpin
        do ij = 1, nXoo(s)
          i = getIJ(ij + soo(s), 1)
          j = getIJ(ij + soo(s), 2)
          ijs = iaTrans(i, j, s)
          woo(ij,s) = woo(ij,s) + vecHoo(ijs)
        end do
      end do
    end if

    ! Divide diagonal elements of W_ij by 2.
    do s = 1, nSpin
      do ij = 1, nXoo(s)
        i = getIJ(ij + soo(s), 1)
        j = getIJ(ij + soo(s), 2)
        if (i == j) then
          woo(ij,s) = 0.5_dp * woo(ij,s)
        end if
      end do
    end do

  end subroutine grndToExcDensityMatrices


  !> Build right hand side of the equation for the Z-vector and those parts of the W-vectors which
  !! do not depend on Z. Modified version of getZVectorEqRHS for state-to-state NA couplings.
  !! Furche PCCP 21 18999 (2019)
  !! Here the + (symmetric) part of RHS, T and (omega_m-omega_n) * W (stored as W) is computed. 
  subroutine getNadiaZvectorEqRHS(tHybridXc, xpyn, xmyn, xpym, xmym, win, env, denseDesc, homo,& 
      & transChrg, getIA, getIJ, getAB, iatrans, natom, species0, grndEigVal, ovrXev, grndEigVecs,&
      &  gammaMat, lrGamma, spinW, omegaAvg, sym, rhs, t, wov, woo, wvv)

    !> Is calculation range-separated?
    logical, intent(in) :: tHybridXc

    !> X+Y Furche term for excited state n
    real(dp), intent(in) :: xpyn(:)

    !> X-Y Furche term for excited state n
    real(dp), intent(in) :: xmyn(:)

    !> X+Y Furche term for excited state m
    real(dp), intent(in) :: xpym(:)

    !> X-Y Furche term for excited state m
    real(dp), intent(in) :: xmym(:)

    !> Index array for single particle transitions
    integer, intent(in) :: win(:)

    !> Environment settings
    type(TEnvironment), intent(inout) :: env

    !> Dense matrix descriptor
    type(TDenseDescr), intent(in) :: denseDesc 

    !> Highest occupied level
    integer, intent(in) :: homo(:)

    !> Machinery for transition charges between single particle levels
    type(TTransCharges), intent(in) :: transChrg

    !> Index array between transitions in square and 1D representations
    integer, intent(in) :: getIA(:,:)

    !> Index array for vir-vir transitions
    integer, intent(in) :: getIJ(:,:)

    !> Index array for occ-occ transitions
    integer, intent(in) :: getAB(:,:)

    !> Index array from orbital pairs to compound index
    integer, intent(in) :: iatrans(:,:,:)

    !> Number of central cell atoms
    integer, intent(in) :: natom

    !> Central cell chemical species
    integer, intent(in) :: species0(:)

    !> Ground state wavefunctions
    real(dp), intent(in) :: grndEigVal(:,:)

    !> Overlap times ground state wavefunctions
    real(dp), intent(in) :: ovrXev(:,:,:)

    !> Ground state wavefunctions
    real(dp), intent(in) :: grndEigVecs(:,:,:)

    !> Softened coulomb matrix
    real(dp), intent(in) :: gammaMat(:,:)

    !> Softened coulomb matrix, long-range corrected
    real(dp), allocatable, intent(in) :: lrGamma(:,:)

    !> Ground state spin derivatives for each species
    real(dp), intent(in) :: spinW(:)

    !> Average excitation energy of states n and m
    real(dp), intent(in) :: omegaAvg

    !> Symmetry of the transitions
    character, intent(in) :: sym

    !> Right hand side (P+Q)
    real(dp), intent(out) :: rhs(:)

    !> T matrix (non-symmetric)
    real(dp), intent(out) :: t(:,:,:)

    !> W^+ vector occupied-virtual part
    real(dp), intent(out) :: wov(:)

    !> W^+ vector occupied part
    real(dp), intent(out) :: woo(:,:)

    !> W^+ vector virtual part
    real(dp), intent(out) :: wvv(:,:)

    real(dp), allocatable :: xpyq(:), qTr(:), gamxpyq(:), qgamxpyq(:,:), gamqt(:)
    real(dp), allocatable :: xpyqds(:), gamxpyqds(:)
    real(dp), allocatable :: vecHvvXorY(:), vecHooXorY(:), vecHovT(:), vecHvvT(:)
    real(dp), allocatable :: vecHvvXpY(:), vecHvvXmY(:), vecHooXpY(:), vecHooXmY(:)
    real(dp), allocatable :: vecHooT(:)
    integer :: nxov
    integer, allocatable :: nxoo(:), nxvv(:), nvir(:)
    integer :: i, j, a, b, ias, ibs, abs, ij, ab, jas, ijs, s, nSpin, soo(2), svv(2), nOrb
    real(dp) :: ptmp1, ptmp2, tmp3, tmp4, tmp1
    logical :: tSpin

    nxov = size(rhs)
    nOrb = size(ovrXev, dim=1)

    allocate(xpyq(natom))
    allocate(qTr(natom))
    allocate(gamxpyq(natom))
    allocate(gamqt(natom))

    t(:,:,:) = 0.0_dp
    rhs(:) = 0.0_dp
    wov(:) = 0.0_dp
    woo(:,:) = 0.0_dp
    wvv(:,:) = 0.0_dp

    nSpin = size(t, dim=3)

    allocate(nxoo(nSpin))
    allocate(nxvv(nSpin))
    allocate(nvir(nSpin))

    nxoo(:) = (homo * (homo + 1)) / 2
    nvir(:) = size(t, dim=1) - homo
    nxvv(:) = (nvir * (nvir + 1)) / 2

    !! transition charges use compound index ijs = ij + soo(s)
    soo(:) = [0, nxoo(1)]
    svv(:) = [0, nxvv(1)]

    allocate(qgamxpyq(max(maxval(nxoo), maxval(nxvv)), size(homo)))
    allocate(vecHooXorY(sum(nxoo)))
    allocate(vecHvvXorY(sum(nxvv)))

    if (nSpin == 2) then
      tSpin = .true.
      allocate(xpyqds(natom))
      allocate(gamxpyqds(natom))
    else
      tSpin = .false.
    end if

    ! Build state-to-state 1TDM and W (eq. 42 in Furche PCCP)
    ! We are symmetrizing the non-symmetric T of Furche
    ! Factor 1 / (1 + delta_ab) for W follows later
    do ias = 1, nxov
      call indxov(win, ias, getIA, i, a, s)

      ! BA: is T_aa = 0?
      do b = homo(s) + 1, a
        ibs = iatrans(i, b, s)
        ab = iaTrans(a, b, s) - svv(s)
        ptmp1 = xpyn(ias) * xpym(ibs) + xmyn(ias) * xmym(ibs)&
              & + xpym(ias) * xpyn(ibs) + xmym(ias) * xmyn(ibs)

        ptmp2 = (xpyn(ias) * xmym(ibs) + xmyn(ias) * xpym(ibs)&
              & + xpym(ias) * xmyn(ibs) + xmym(ias) * xpyn(ibs))

        tmp3 = xpyn(ias) * xpym(ibs) + xmyn(ias) * xmym(ibs)
        tmp4 = xpyn(ibs) * xpym(ias) + xmyn(ibs) * xmym(ias)

        ! Set t(a,b,s) = t(a,b,s) + 0.5_dp * tmp3 for asymmetric T
        t(a,b,s) = t(a,b,s) + 0.25_dp * (tmp3 + tmp4)
        wvv(ab,s) = wvv(ab,s) + 0.5_dp * grndEigVal(i,s) * ptmp1&
                   & + 0.5_dp * omegaAvg * ptmp2

        ! To prevent double counting
        if (a /= b) then
          ! Set t(b,a,s) = t(b,a,s) + 0.5_dp * tmp4 for asymmetric T
          t(b,a,s) = t(b,a,s) + 0.25_dp * (tmp3 + tmp4)
        end if

      end do

      do j = i, homo(s)
        jas = iatrans(j,a,s)

        ij = iatrans(i, j, s) - soo(s)

        ptmp1 = (xpyn(ias) * xpym(jas) + xmyn(ias) * xmym(jas)&
              & + xpym(ias) * xpyn(jas) + xmym(ias) * xmyn(jas))

        ptmp2 = (xpyn(ias) * xmym(jas) + xmyn(ias) * xpym(jas)&
              & + xpym(ias) * xmyn(jas) + xmym(ias) * xpyn(jas))

        tmp3 = xpyn(ias) * xpym(jas) + xmyn(ias) * xmym(jas)
        tmp4 = xpyn(jas) * xpym(ias) + xmyn(jas) * xmym(ias)

        t(i,j,s) = t(i,j,s) - 0.25_dp * (tmp3 + tmp4)
        woo(ij,s) = woo(ij,s) - 0.5_dp * grndEigVal(a,s) * ptmp1 + 0.5_dp * omegaAvg * ptmp2

        ! To prevent double counting
        if (i /= j) then
          t(j,i,s) = t(j,i,s) - 0.25_dp * (tmp3 + tmp4)
        end if

      end do

    end do

    ! Terms for (P+-Q) of form (X+Y)^m_ib H^+_ab[(X+Y)^n]  
    call getHplusXYfr(sym, nxoo, nxvv, nAtom, getIA, getIJ, getAB, win, env, denseDesc, &
      & species0, ovrXev, grndEigVecs, gammaMat, spinW, transChrg, xpyn, vecHoo=vecHooXorY,&
      & vecHvv=vecHvvXorY)

    do ias = 1, nxov
      call indxov(win, ias, getIA, i, a, s)
      do b = homo(s) + 1, a
        abs = iatrans(a, b, s)
        ibs = iatrans(i, b, s)
        ! For the forces, we have a factor of 2 here
        rhs(ias) = rhs(ias) - xpym(ibs) * vecHvvXorY(abs)
        ! Since vecHvvXpY has only upper triangle
        if (a /= b) then
          rhs(ibs) = rhs(ibs) - xpym(ias) * vecHvvXorY(abs)
        end if
      end do

      do j = i, homo(s)
        jas = iatrans(j, a, s)
        ijs = iatrans(i, j, s)
        ! For the forces, we have a factor of 2 here
        tmp1 = xpym(jas) * vecHooXorY(ijs)
        rhs(ias) = rhs(ias) + tmp1
        wov(ias) = wov(ias) + tmp1
        if (i /= j) then
           tmp1 = xpym(ias) * vecHooXorY(ijs)
           rhs(jas) = rhs(jas) + tmp1
           wov(jas) = wov(jas) + tmp1
        end if
      end do

    end do

    ! Now m <-> n
    call getHplusXYfr(sym, nxoo, nxvv, nAtom, getIA, getIJ, getAB, win, env, denseDesc,&
      & species0, ovrXev, grndEigVecs, gammaMat, spinW, transChrg, xpym, vecHoo=vecHooXorY,&
      & vecHvv=vecHvvXorY)

    do ias = 1, nxov
      call indxov(win, ias, getIA, i, a, s)
      do b = homo(s) + 1, a
        abs = iatrans(a, b, s)
        ibs = iatrans(i, b, s)
        ! For the forces, we have a factor of 2 here
        rhs(ias) = rhs(ias) - xpyn(ibs) * vecHvvXorY(abs)
        ! Since vecHvvXpY has only upper triangle
        if (a /= b) then
          rhs(ibs) = rhs(ibs) - xpyn(ias) * vecHvvXorY(abs)
        end if
      end do

      do j = i, homo(s)
        jas = iatrans(j, a, s)
        ijs = iatrans(i, j, s)
        ! For the forces, we have a factor of 2 here
        tmp1 = xpyn(jas) * vecHooXorY(ijs)
        rhs(ias) = rhs(ias) + tmp1
        wov(ias) = wov(ias) + tmp1
        if (i /= j) then
           tmp1 = xpyn(ias) * vecHooXorY(ijs)
           rhs(jas) = rhs(jas) + tmp1
           wov(jas) = wov(jas) + tmp1
        end if
      end do

    end do

    allocate(vecHovT(nxov))
    allocate(vecHooT(sum(nxoo)))
    allocate(vecHvvT(sum(nxvv)))

    ! -RHS^+ += - H^+_ia[T^+]
    call getHplusMfr(3, nxoo, nxvv, nxov, nAtom, getIA, getIJ, getAB, win, env, denseDesc,&
      & species0, ovrXev, grndEigVecs, gammaMat, spinW, transChrg, t, vecHovT)

    rhs = rhs - vecHovT

    ! Woo^+ += 0.5 * H^+_ij[T+Z] / Omega_mn, Z part computed later 
    call getHplusMfr(1, nxoo, nxvv, nxov, nAtom, getIA, getIJ, getAB, win, env, denseDesc,&
      &  species0, ovrXev, grndEigVecs, gammaMat, spinW, transChrg, t, vecHooT)

    do s = 1, nSpin
      do ij = 1, nxoo(s)
        ijs = ij + soo(s)
        woo(ij,s) = woo(ij,s) + vecHooT(ijs)
      end do
    end do

    ! Contributions due to range-separation
    if (tHybridXc) then
      allocate(vecHvvXpY(sum(nxvv)))
      allocate(vecHvvXmY(sum(nxvv)))
      allocate(vecHooXpY(sum(nxoo)))
      allocate(vecHooXmY(sum(nxoo)))

      !> Long-range part of H^+[(X+Y)^n] or H^-[(X-Y)^n] for occ-occ and vir-vir comp. of H
      call getHvvXY(1, nxvv, homo, natom, iatrans, getIA, getAB, win, env, denseDesc,&
          & ovrXev, grndEigVecs, lrGamma, transChrg, xpyn, vecHvvXpY)

      call getHvvXY(-1, nxvv, homo, natom, iatrans, getIA, getAB, win, env, denseDesc,&
          & ovrXev, grndEigVecs, lrGamma, transChrg, xmyn, vecHvvXmY)

      call getHooXY( 1, nxoo, homo, natom, iatrans, getIA, getIJ, win, env, denseDesc,&
          & ovrXev, grndEigVecs, lrGamma, transChrg, xpyn, vecHooXpY)

      call getHooXY(-1, nxoo, homo, natom, iatrans, getIA, getIJ, win, env, denseDesc,&
          & ovrXev, grndEigVecs, lrGamma, transChrg, xmyn, vecHooXmY)

      do ias = 1, nxov

        call indXov(win, ias, getIA, i, a, s)
        do b = homo(s) + 1, nOrb
          ibs = iaTrans(i, b, s)
          abs = iaTrans(a, b, s)
          rhs(ias) = rhs(ias) - cExchange * 0.5_dp * xpym(ibs) * vecHvvXpY(abs)
          if (a >= b) then
            rhs(ias) = rhs(ias) - cExchange * 0.5_dp * xmym(ibs) * vecHvvXmY(abs)
          ! Only a>b is stored in vecHvvXmY, which is anti-symmetric
          else
            rhs(ias) = rhs(ias) + cExchange * 0.5_dp * xmym(ibs) * vecHvvXmY(abs)
          end if
        end do

        do j = 1, homo(s)
          jas = iaTrans(j, a, s)
          ijs = iaTrans(i, j, s)
          rhs(ias) = rhs(ias) + cExchange * 0.5_dp * xpym(jas) * vecHooXpY(ijs)
          wov(ias) = wov(ias) + cExchange * 0.5_dp * xpym(jas) * vecHooXpY(ijs)
          if (i >= j) then
            rhs(ias) = rhs(ias) + cExchange * 0.5_dp * xmym(jas) * vecHooXmY(ijs)
            wov(ias) = wov(ias) + cExchange * 0.5_dp * xmym(jas) * vecHooXmY(ijs)
          else
            rhs(ias) = rhs(ias) - cExchange * 0.5_dp * xmym(jas) * vecHooXmY(ijs)
            wov(ias) = wov(ias) - cExchange * 0.5_dp * xmym(jas) * vecHooXmY(ijs)
          end if
        end do

      end do

      ! Now n <-> m

      !> Long-range part of H^+[(X+Y)^n] or H^-[(X-Y)^n] for occ-occ and vir-vir comp. of H
      call getHvvXY( 1, nxvv, homo, natom, iatrans, getIA, getAB, win, env, denseDesc,&
          & ovrXev, grndEigVecs, lrGamma, transChrg, xpym, vecHvvXpY)

      call getHvvXY(-1, nxvv, homo, natom, iatrans, getIA, getAB, win, env, denseDesc,&
          & ovrXev, grndEigVecs, lrGamma, transChrg, xmym, vecHvvXmY)

      call getHooXY( 1, nxoo, homo, natom, iatrans, getIA, getIJ, win, env, denseDesc,&
          & ovrXev, grndEigVecs, lrGamma, transChrg, xpym, vecHooXpY)

      call getHooXY(-1, nxoo, homo, natom, iatrans, getIA, getIJ, win, env, denseDesc,&
          & ovrXev, grndEigVecs, lrGamma, transChrg, xmym, vecHooXmY)

      do ias = 1, nxov

        call indXov(win, ias, getIA, i, a, s)
        do b = homo(s) + 1, nOrb
          ibs = iaTrans(i, b, s)
          abs = iaTrans(a, b, s)
          rhs(ias) = rhs(ias) - cExchange * 0.5_dp * xpyn(ibs) * vecHvvXpY(abs)
          if (a >= b) then
            rhs(ias) = rhs(ias) - cExchange * 0.5_dp * xmyn(ibs) * vecHvvXmY(abs)
          ! Only a>b is stored in vecHvvXmY, which is anti-symmetric
          else
            rhs(ias) = rhs(ias) + cExchange * 0.5_dp * xmyn(ibs) * vecHvvXmY(abs)
          end if
        end do

        do j = 1, homo(s)
          jas = iaTrans(j, a, s)
          ijs = iaTrans(i, j, s)
          rhs(ias) = rhs(ias) + cExchange * 0.5_dp * xpyn(jas) * vecHooXpY(ijs)
          wov(ias) = wov(ias) + cExchange * 0.5_dp * xpyn(jas) * vecHooXpY(ijs)
          if (i >= j) then
            rhs(ias) = rhs(ias) + cExchange * 0.5_dp * xmyn(jas) * vecHooXmY(ijs)
            wov(ias) = wov(ias) + cExchange * 0.5_dp * xmyn(jas) * vecHooXmY(ijs)
          else
            rhs(ias) = rhs(ias) - cExchange * 0.5_dp * xmyn(jas) * vecHooXmY(ijs)
            wov(ias) = wov(ias) - cExchange * 0.5_dp * xmyn(jas) * vecHooXmY(ijs)
          end if
        end do

      end do

      ! -RHS^+ += - H^+_ia[T^+]
      call getHovT(nxoo, nxvv, homo, natom, iatrans, getIA, getIJ, getAB, win,&
        & env, denseDesc, ovrXev, grndEigVecs, lrGamma, transChrg, t, vecHovT)

      rhs = rhs - cExchange * vecHovT

      ! Woo^+ += 0.5 * H^+_ij[T+Z] / Omega_mn, Z part computed later 
      call getHooT(nxov, nxoo, homo, natom, iatrans, getIA, getIJ, win, env, denseDesc,&
       & ovrXev, grndEigVecs, lrGamma, transChrg, t, vecHooT)

      do s = 1, nSpin
        do ij = 1, nxoo(s)
          ijs = ij + soo(s)
          woo(ij, s) = woo(ij, s) + cExchange * vecHooT(ijs)
        end do
      end do

    endif

  end subroutine getNadiaZvectorEqRHS


  !> Write out non-adiabatic coupling vectors.
  subroutine writeNACV(iLev, jLev, fdTagged, taggedWriter, nacv)

    !> Start level for coupling
    integer, intent(in) :: iLev

    !> End level for coupling
    integer, intent(in) :: jLev

    !> File descriptor for the tagged data output
    type(TFileDescr), intent(in) :: fdTagged

    !> Tagged writer
    type(TTaggedWriter), intent(inout) :: taggedWriter

    !> Non-adiabatic coupling vector
    real(dp), intent(in) :: nacv(:,:,:)

    type(TFileDescr) :: fdNaCoupl
    integer :: ii, iNac, nCoupLev, mCoupLev

    call openFile(fdNaCoupl, naCouplingOut, mode="w")

    iNac = 0
    do nCoupLev = iLev, jLev-1
      do mCoupLev = nCoupLev+1, jLev
         iNac = iNac + 1
         write(fdNaCoupl%unit, '(2(I4,2x))') nCoupLev, mCoupLev
         do ii = 1, size(nacv, dim=2)
           write(fdNaCoupl%unit, '(3(E20.12,2x))') nacv(1,ii,iNac), nacv(2,ii,iNac), nacv(3,ii,iNac)
         end do
      end do
    end do

    call closeFile(fdNaCoupl)

    if (fdTagged%isConnected()) call taggedWriter%write(fdTagged%unit, tagLabels%nacv, nacv)

  end subroutine writeNACV


  !> Calculation of nacv using gradient routine
  subroutine addNadiaGradients(sym, nxov, natom, species0, env, denseDesc, norb, homo, getIA,       &
      & getIJ, getAB, win, grndEigVecs, pc, ovrXev, dq_ud, dqex, gammaMat, lrGamma, HubbardU,   &
      & spinW, shift, woo, wov, wvv, transChrg, xpyn, xmyn, xpym, xmym, coord0, orb, skHamCont, &
      & skOverCont, derivator, rhoSqr, deltaRho, tHybridXc, hybridXc, nacv)

    !> Symmetry of the transition
    character, intent(in) :: sym

    !> Number of single particle transitions to include
    integer, intent(in) :: nxov

    !> Number of central cell atoms
    integer, intent(in) :: natom

    !> Central cell chemical species
    integer, intent(in) :: species0(:)

    !> Environment settings
    type(TEnvironment), intent(inout) :: env

    !> Dense matrix descriptor
    type(TDenseDescr), intent(in) :: denseDesc
    
    !> number of orbitals for ground state system
    integer, intent(in) :: norb

    !> Number of highest occupied state in ground state
    integer, intent(in) :: homo(:)

    !> Index array from composite occ-vir transition index to specific single particle states
    integer, intent(in) :: getIA(:,:)

    !> Index array from composite occ-occ transition index to specific single particle states
    integer, intent(in) :: getIJ(:,:)

    !> Index array from composite vir-vir transition index to specific single particle states
    integer, intent(in) :: getAB(:,:)

    !> Single particle transition index
    integer, intent(in) :: win(:)

    !> Ground state eigenvectors
    real(dp), intent(in) :: grndEigVecs(:,:,:)

    !> Transition density matrix
    real(dp), intent(in) :: pc(:,:,:)

    !> Overlap times ground state eigenvectors
    real(dp), intent(in) :: ovrXev(:,:,:)

    !> Ground state gross charges
    real(dp), intent(in) :: dq_ud(:,:)

    !> Charge differences from ground to excited state
    real(dp), intent(in) :: dqex(:,:)

    !> Softened coulomb matrix
    real(dp), intent(in) :: gammaMat(:,:)

    !> Electrostatic matrix, long-range corrected
    real(dp), allocatable, intent(in) :: lrGamma(:,:)

    !> Ground state Hubbard U values
    real(dp), intent(in) :: HubbardU(:)

    !> Ground state spin derivatives for each species
    real(dp), intent(in) :: spinW(:)

    !> Ground state potentials (shift vector)
    real(dp), intent(in) :: shift(:)

    !> W vector occupied part
    real(dp), intent(inout) :: woo(:,:)

    !> W vector occupied-virtual part
    real(dp), intent(inout) :: wov(:)

    !> W vector virtual part
    real(dp), intent(inout) :: wvv(:,:)

    !> Machinery for transition charges between single particle levels
    type(TTransCharges), intent(in) :: transChrg

    !> X+Y Furche term (state n)
    real(dp), intent(in) :: xpyn(:)

    !> X-Y Furche term (state n)
    real(dp), intent(in) :: xmyn(:)

    !> X+Y Furche term (state m)
    real(dp), intent(in) :: xpym(:)

    !> X-Y Furche term (state m)
    real(dp), intent(in) :: xmym(:)

    !> Central cell atomic coordinates
    real(dp), intent(in) :: coord0(:,:)

    !> Data type for atomic orbital information
    type(TOrbitals), intent(in) :: orb

    !> H0 data
    type(TSlakoCont), intent(in) :: skHamCont

    !> Overlap data
    type(TSlakoCont), intent(in) :: skOverCont

    !> Differentiator for the non-scc matrices
    class(TNonSccDiff), intent(in) :: derivator

    !> Ground state density matrix
    real(dp), intent(in) :: rhoSqr(:,:,:)

    !> Difference density matrix (vs. uncharged atoms)
    real(dp), intent(inout), optional :: deltaRho(:,:,:)

    !> Is calculation range-separated?
    logical, intent(in) :: tHybridXc

    !> Data for hybrid xc-functional calculation
    class(THybridXcFunc), allocatable, intent(inout) :: hybridXc

    !> Resulting non-adiabatic coupling
    real(dp), intent(out) :: nacv(:,:)

    real(dp), allocatable :: shift_excited(:,:), xpyq(:,:), xpyqds(:,:)
    real(dp), allocatable :: shxpyq(:,:,:), xpycc(:,:,:,:), wcc(:,:,:), tmp5(:), tmp7(:,:), tmp11(:)
    real(dp), allocatable :: qTr(:), temp(:), dq(:), dm(:), dsigma(:)
    real(dp), allocatable :: dH0(:,:,:), dSo(:,:,:)
    real(dp), allocatable :: Dens(:,:), SpinDens(:,:)
    real(dp), allocatable :: xmycc(:,:,:,:), xpyas(:,:,:,:), xmyas(:,:,:,:)
    real(dp), allocatable :: overlap(:,:), lrGammaOrb(:,:), gammaLongRangePrime(:,:,:)
    real(dp), allocatable :: PS(:,:,:), DS(:,:,:), SPS(:,:,:), SDS(:,:,:), SX(:,:,:,:)
    real(dp), allocatable :: XS(:,:,:,:), SXS(:,:,:,:), SY(:,:,:,:), YS(:,:,:,:), SYS(:,:,:,:)
    real(dp), allocatable :: xpy(:,:), xmy(:,:)
    integer :: ia, i, j, a, b, ab, ij, m, n, mu, nu, xyz, iAt1, iAt2, ka
    integer :: indalpha, indalpha1, indbeta, indbeta1, soo(2), svv(2)
    integer :: iSp1, iSp2, iSpin, nSpin, iState
    real(dp) :: tmp1, tmp2, tmp3, tmp4, tmp6, tmp8, tmp9, tmp10, rab
    real(dp) :: diffvec(3), dgab(3), tmpVec(3), tmp3a, tmp3b, tmprs, tmprs2, tmps(2)
    integer, allocatable :: nxoo(:), nxvv(:), nvir(:), species(:)
    logical :: tSpin

    nSpin = size(grndEigVecs, dim=3)
    tSpin = (nSpin == 2)

    allocate(shift_excited(natom, nSpin))
    allocate(xpyq(natom, 2))
    allocate(shxpyq(natom, nSpin, 2))
    allocate(xpycc(norb, norb, nSpin, 2))
    allocate(wcc(norb, norb, nSpin))
    allocate(qTr(natom))
    allocate(temp(norb))
    allocate(tmp5(nSpin))
    allocate(tmp7(nSpin, 2))

    !! This should be changed to save memory
    allocate(xpy(nxov, 2))
    allocate(xmy(nxov, 2))
    xpy(:,1) = xpyn
    xpy(:,2) = xpym
    xmy(:,1) = xmyn
    xmy(:,2) = xmym

    allocate(Dens(norb, norb))
    Dens(:,:) = sum(rhoSqr, dim=3)

    allocate(dH0(orb%mOrb, orb%mOrb, 3))
    allocate(dSo(orb%mOrb, orb%mOrb, 3))

    allocate(nxoo(nSpin))
    allocate(nxvv(nSpin))
    allocate(nvir(nSpin))

    nxoo(:) = (homo * (homo + 1)) / 2
    nvir(:) = norb - homo
    nxvv(:) = (nvir * (nvir + 1)) / 2

    soo(:) = [0, nxoo(1)]
    svv(:) = [0, nxvv(1)]

    allocate(dq(natom))
    dq(:) = dq_ud(:,1)

    if (tSpin) then
      allocate(dm(natom))
      allocate(xpyqds(natom, 2))
      xpyqds = 0.0_dp
      allocate(tmp11(nSpin))

      allocate(SpinDens(norb, norb))
      SpinDens(:,:) = rhoSqr(:,:,1) - rhoSqr(:,:,2)

      allocate(dsigma(2))
      dsigma(1) = 1.0_dp
      dsigma(2) = -1.0_dp
      dm(:) = dq_ud(:,2)
    end if

    if (tHybridXc) then
      allocate(xmycc(norb, norb, nSpin, 2))
      allocate(xpyas(norb, norb, nSpin, 2))
      allocate(xmyas(norb, norb, nSpin, 2))
      allocate(PS(norb, norb, nSpin))
      allocate(DS(norb, norb, nSpin))
      allocate(SPS(norb, norb, nSpin))
      allocate(SDS(norb, norb, nSpin))
      allocate(SX(norb, norb, nSpin, 2))
      allocate(XS(norb, norb, nSpin, 2))
      allocate(SXS(norb, norb, nSpin, 2))
      allocate(SY(norb, norb, nSpin, 2))
      allocate(YS(norb, norb, nSpin, 2))
      allocate(SYS(norb, norb, nSpin, 2))
      allocate(overlap(norb, norb))
      allocate(lrGammaOrb(norb, norb))
      allocate(gammaLongRangePrime(3, nAtom, nAtom))

      ! Symmetrize deltaRho
      do mu = 1, norb
        do nu = mu + 1, norb
          deltaRho(mu,nu,:) = deltaRho(nu,mu,:)
        end do
      end do

      ! Compute long-range gamma derivative
      gammaLongRangePrime(:,:,:) = 0.0_dp
      call hybridXc%getCentralCellSpecies(species)
      do iAt1 = 1, nAtom
        do iAt2 = 1, nAtom
          if (iAt1 /= iAt2) then
            call getDirectionalCamGammaPrimeValue(hybridXc, tmpVec, iAt1, iAt2)
            gammaLongRangePrime(:, iAt1, iAt2) = tmpVec
          end if
        end do
      end do

      ! Symmetrize S (can't we get S from caller?)
      call getSqrS(coord0, nAtom, skOverCont, orb, denseDesc%iAtomStart, species0, overlap)
      call getSqrGamma(nAtom, lrGamma, denseDesc%iAtomStart, lrGammaOrb)

    end if

    nacv = 0.0_dp

    ! Excited state potentials at atomic sites
    do iSpin = 1, nSpin
      call hemv(shift_excited(:,iSpin), gammaMat, dqex(:,iSpin))
    end do

    ! xypq(alpha) = sum_ia (X+Y)_ia q^ia(alpha)
    ! Complexity norb * norb * norb
    xpyq = 0.0_dp
    do iState = 1, 2
      call transChrg%qMatVec(env, denseDesc, ovrXev, grndEigVecs, getIA, win,&
           & xpy(:,iState), xpyq(:,iState))
      ! complexity norb * norb
      shxpyq(:,:,iState) = 0.0_dp
      if (.not. tSpin) then
        if (sym == "S") then
          call hemv(shxpyq(:,1,iState), gammaMat, xpyq(:,iState))
        else
          shxpyq(:,1,iState) = xpyq(:,iState) * spinW(species0)
        end if
      else
        call transChrg%qMatVecDs(env, denseDesc, ovrXev, grndEigVecs, getIA, win,&
             & xpy(:,iState), xpyqds(:,iState))
        do iSpin = 1, nSpin
          call hemv(shxpyq(:,iSpin,iState), gammaMat, xpyq(:,iState))
          shxpyq(:,iSpin,iState) = shxpyq(:,iSpin,iState) + dsigma(iSpin)&
               & * spinW(species0) * xpyqds(:,iState)
          shxpyq(:,iSpin,iState) = 0.5_dp * shxpyq(:,iSpin,iState)
        end do
      end if

      ! Calculate xpycc
      ! (xpycc)_{mu nu} = sum_{ia} (X + Y)_{ia} (grndEigVecs(mu,i)grndEigVecs(nu,a)
      ! + grndEigVecs(nu,i)grndEigVecs(mu,a))
      ! Complexity norb * norb * norb
      !
      ! xpycc(mu,nu) = sum_ia (X+Y)_ia grndEigVecs(mu,i) grndEigVecs(nu,a)
      ! xpycc(mu, nu) += sum_ia (X+Y)_ia grndEigVecs(mu,a) grndEigVecs(nu,i)
      xpycc(:,:,:,iState) = 0.0_dp
      do ia = 1, nxov
        call indxov(win, ia, getIA, i, a, iSpin)
        ! Should replace with DSYR2 call:
        do nu = 1, norb
          do mu = 1, norb
            xpycc(mu,nu,iSpin,iState) = xpycc(mu,nu,iSpin,iState) + xpy(ia,iState)&
                & * (grndEigVecs(mu,i,iSpin)*grndEigVecs(nu,a,iSpin)&
                & + grndEigVecs(mu,a,iSpin)*grndEigVecs(nu,i,iSpin))
          end do
        end do
      end do
    end do

    if (tHybridXc) then

      xmycc = 0.0_dp
      xpyas = 0.0_dp
      xmyas = 0.0_dp

      do iState = 1,2
        ! Asymmetric contribution: xmycc_as = sum_ias (X-Y)_ias c_mas c_nis
        do ia = 1, nxov
          call indxov(win, ia, getIA, i, a, iSpin)
          ! Should replace with DSYR2 call:
          do nu = 1, norb
            do mu = 1, norb
               xmycc(mu,nu,iSpin,iState) = xmycc(mu,nu,iSpin,iState) + xmy(ia,iState) *&
                & ( grndEigVecs(mu,i,iSpin) * grndEigVecs(nu,a,iSpin)&
                & + grndEigVecs(mu,a,iSpin) * grndEigVecs(nu,i,iSpin) )
               xpyas(mu,nu,iSpin,iState) = xpyas(mu,nu,iSpin,iState) + xpy(ia,iState) *&
                &  grndEigVecs(mu,i,iSpin) * grndEigVecs(nu,a,iSpin)
               xmyas(mu,nu,iSpin,iState) = xmyas(mu,nu,iSpin,iState) + xmy(ia,iState) *&
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
          call symm(DS(:,:,iSpin), 'R', overlap, deltaRho(:,:,iSpin), 'U', 1.0_dp, 0.0_dp, nOrb,&
              & nOrb)
          call symm(SDS(:,:,iSpin), 'L', overlap, DS(:,:,iSpin), 'U', 1.0_dp, 0.0_dp, nOrb, nOrb)
          call symm(XS(:,:,iSpin,iState), 'R', overlap, xpyas(:,:,iSpin,iState), 'U', 1.0_dp,&
              & 0.0_dp, nOrb, nOrb)
          call symm(SX(:,:,iSpin,iState), 'L', overlap, xpyas(:,:,iSpin,iState), 'U', 1.0_dp,&
              & 0.0_dp, nOrb, nOrb)
          call symm(SXS(:,:,iSpin,iState), 'L', overlap, XS(:,:,iSpin,iState), 'U', 1.0_dp, 0.0_dp,&
              & nOrb, nOrb)
          call symm(YS(:,:,iSpin,iState), 'R', overlap, xmyas(:,:,iSpin,iState), 'U', 1.0_dp,&
              & 0.0_dp, nOrb, nOrb)
          call symm(SY(:,:,iSpin,iState), 'L', overlap, xmyas(:,:,iSpin,iState), 'U', 1.0_dp,&
              & 0.0_dp, nOrb, nOrb)
          call symm(SYS(:,:,iSpin,iState), 'L', overlap, YS(:,:,iSpin,iState), 'U', 1.0_dp, 0.0_dp,&
              & nOrb, nOrb)
        end do
      end do

    end if

    ! Calculate wcc = c_mu,i * W_ij * c_j,nu. We have only W_ab b > a and W_ij j > i:
    ! wcc(m,n) = sum_{pq, p <= q} w_pq (grndEigVecs(mu,p)grndEigVecs(nu,q)
    ! + grndEigVecs(nu,p)grndEigVecs(mu,q))
    ! Complexity norb * norb * norb

    ! Calculate the occ-occ part
    wcc(:,:,:) = 0.0_dp

    do iSpin = 1, nSpin
      do ij = 1, nxoo(iSpin)
        i = getIJ(ij + soo(iSpin), 1)
        j = getIJ(ij + soo(iSpin), 2)
        ! Replace with DSYR2 call:
        do mu = 1, norb
          do nu = 1, norb
            wcc(mu,nu,iSpin) = wcc(mu,nu,iSpin) + woo(ij,iSpin) *&
                & ( grndEigVecs(mu,i,iSpin)*grndEigVecs(nu,j,iSpin)&
                & + grndEigVecs(mu,j,iSpin)*grndEigVecs(nu,i,iSpin) )
          end do
        end do

      end do
    end do

    ! Calculate the occ-virt part: the same way as for xpycc
    do ia = 1, nxov
      call indxov(win, ia, getIA, i, a, iSpin)
      ! Again replace with DSYR2 call:
      do nu = 1, norb
        do mu = 1, norb
          wcc(mu,nu,iSpin) = wcc(mu,nu,iSpin) + wov(ia) *&
              & ( grndEigVecs(mu,i,iSpin)*grndEigVecs(nu,a,iSpin)&
              & + grndEigVecs(mu,a,iSpin)*grndEigVecs(nu,i,iSpin) )
        end do
      end do
    end do

    ! Calculate the virt - virt part
    do iSpin = 1, nSpin
      do ab = 1, nxvv(iSpin)
        a = getAB(ab + svv(iSpin), 1)
        b = getAB(ab + svv(iSpin), 2)
        ! Replace with DSYR2 call:
        do mu = 1, norb
          do nu = 1, norb
            wcc(mu,nu,iSpin) = wcc(mu,nu,iSpin) + wvv(ab,iSpin) *&
                & ( grndEigVecs(mu,a,iSpin)*grndEigVecs(nu,b,iSpin)&
                & + grndEigVecs(mu,b,iSpin)*grndEigVecs(nu,a,iSpin) )
          end do
        end do

      end do
    end do

    ! Now calculating the force complexity : norb * norb * 3

    ! As have already performed norb**3 operation to get here, calculate for all atoms

    ! BA: only for non-periodic systems!
    do iAt1 = 1, nAtom
      indalpha = denseDesc%iAtomStart(iAt1)
      indalpha1 = denseDesc%iAtomStart(iAt1 + 1) -1
      iSp1 = species0(iAt1)

      do iAt2 = 1, iAt1 - 1
        indbeta = denseDesc%iAtomStart(iAt2)
        indbeta1 = denseDesc%iAtomStart(iAt2 + 1) -1
        iSp2 = species0(iAt2)

        diffvec = coord0(:,iAt1) - coord0(:,iAt2)
        rab = sqrt(sum(diffvec**2))

        ! Now holds unit vector in direction
        diffvec = diffvec / rab

        ! Calculate the derivative of gamma
        dgab(:) = diffvec * (-1.0_dp / rab**2 - expGammaPrime(rab, HubbardU(iSp1), HubbardU(iSp2)))

        tmp3a = 0.0_dp
        do iSpin = 1, nSpin
          tmp3a = tmp3a + dq(iAt1) * dqex(iAt2,iSpin) + dqex(iAt1,iSpin) * dq(iAt2)
        end do

        ! This term is symmetrized for NACV
        if (.not. tSpin) then
          if (sym == "S") then
            tmp3b = 2.0_dp * (xpyq(iAt1,1) * xpyq(iAt2,2) + xpyq(iAt2,1) * xpyq(iAt1,2))
          else
            tmp3b = 0.0_dp
          end if
        else
          tmp3b = xpyq(iAt1,1) * xpyq(iAt2,2) + xpyq(iAt2,1) * xpyq(iAt1,2)
        end if

        nacv(:,iAt1) = nacv(:,iAt1) + dgab * (tmp3a + tmp3b)
        nacv(:,iAt2) = nacv(:,iAt2) - dgab * (tmp3a + tmp3b)

        tmp5(:) = shift_excited(iAt1,:) + shift_excited(iAt2,:)
        tmp7(:,:) = 2.0_dp * (shxpyq(iAt1,:,:) + shxpyq(iAt2,:,:))

        if (tSpin) then
          tmp9 = spinW(iSp1) * dm(iAt1) + spinW(iSp2) * dm(iAt2)
          tmp11(:) = spinW(iSp1) * dqex(iAt1,:) + spinW(iSp2) * dqex(iAt2,:)
        end if

        if (tHybridXc) then
          tmprs = 0.0_dp
          tmps(:) = 0.0_dp
          do iSpin = 1, nSpin
            do mu = indAlpha, indAlpha1
              do nu = indBeta, indBeta1
                tmprs = tmprs + &
                    & ( 2.0_dp * (PS(mu,nu,iSpin) * DS(nu,mu,iSpin)&
                    & + PS(nu,mu,iSpin) * DS(mu,nu,iSpin))&
                    & + SPS(mu,nu,iSpin) * deltaRho(mu,nu,iSpin)&
                    & + SPS(nu,mu,iSpin) * deltaRho(nu,mu,iSpin)&
                    & + pc(mu,nu,iSpin) * SDS(mu,nu,iSpin)&
                    & + pc(nu,mu,iSpin) * SDS(nu,mu,iSpin) )

                tmprs = tmprs + &
                    & ( xpyas(mu,nu,iSpin,1) * SXS(mu,nu,iSpin,2)&
                    & + xpyas(nu,mu,iSpin,2) * SXS(nu,mu,iSpin,1)&
                    & + SX(mu,nu,iSpin,1) * XS(mu,nu,iSpin,2)&
                    & + SX(nu,mu,iSpin,2) * XS(nu,mu,iSpin,1) )
                tmprs = tmprs + &
                    & ( xpyas(mu,nu,iSpin,2) * SXS(mu,nu,iSpin,1)&
                    & + xpyas(nu,mu,iSpin,1) * SXS(nu,mu,iSpin,2)&
                    & + SX(mu,nu,iSpin,2) * XS(mu,nu,iSpin,1)&
                    & + SX(nu,mu,iSpin,1) * XS(nu,mu,iSpin,2) )

                tmprs = tmprs + 0.5_dp * &
                    & ( XS(mu,nu,iSpin,1) * XS(nu,mu,iSpin,2)&
                    & + XS(nu,mu,iSpin,2) * XS(mu,nu,iSpin,1)&
                    & + SXS(mu,nu,iSpin,1) * xpyas(nu,mu,iSpin,2)&
                    & + SXS(nu,mu,iSpin,2) * xpyas(mu,nu,iSpin,1)&
                    & + xpyas(mu,nu,iSpin,1) * SXS(nu,mu,iSpin,2)&
                    & + xpyas(nu,mu,iSpin,2) * SXS(mu,nu,iSpin,1)&
                    & + SX(mu,nu,iSpin,1) * SX(nu,mu,iSpin,2)&
                    & + SX(nu,mu,iSpin,2) * SX(mu,nu,iSpin,1) )
                tmprs = tmprs + 0.5_dp * &
                    & ( XS(mu,nu,iSpin,2) * XS(nu,mu,iSpin,1)&
                    & + XS(nu,mu,iSpin,1) * XS(mu,nu,iSpin,2)&
                    & + SXS(mu,nu,iSpin,2) * xpyas(nu,mu,iSpin,1)&
                    & + SXS(nu,mu,iSpin,1) * xpyas(mu,nu,iSpin,2)&
                    & + xpyas(mu,nu,iSpin,2) * SXS(nu,mu,iSpin,1)&
                    & + xpyas(nu,mu,iSpin,1) * SXS(mu,nu,iSpin,2)&
                    & + SX(mu,nu,iSpin,2) * SX(nu,mu,iSpin,1)&
                    & + SX(nu,mu,iSpin,1) * SX(mu,nu,iSpin,2) )

                tmprs = tmprs + &
                    & ( xmyas(mu,nu,iSpin,1) * SYS(mu,nu,iSpin,2)&
                    & + xmyas(nu,mu,iSpin,2) * SYS(nu,mu,iSpin,1)&
                    & + SY(mu,nu,iSpin,1) * YS(mu,nu,iSpin,2)&
                    & + SY(nu,mu,iSpin,2) * YS(nu,mu,iSpin,1) )
                tmprs = tmprs + &
                    & ( xmyas(mu,nu,iSpin,2) * SYS(mu,nu,iSpin,1)&
                    & + xmyas(nu,mu,iSpin,1) * SYS(nu,mu,iSpin,2)&
                    & + SY(mu,nu,iSpin,2) * YS(mu,nu,iSpin,1)&
                    & + SY(nu,mu,iSpin,1) * YS(nu,mu,iSpin,2) )

                tmprs = tmprs - 0.5_dp * &
                    & ( YS(mu,nu,iSpin,1) * YS(nu,mu,iSpin,2)&
                    & + YS(nu,mu,iSpin,2) * YS(mu,nu,iSpin,1)&
                    & + SYS(mu,nu,iSpin,1) * xmyas(nu,mu,iSpin,2)&
                    & + SYS(nu,mu,iSpin,2) * xmyas(mu,nu,iSpin,1)&
                    & + xmyas(mu,nu,iSpin,1) * SYS(nu,mu,iSpin,2)&
                    & + xmyas(nu,mu,iSpin,2) * SYS(mu,nu,iSpin,1)&
                    & + SY(mu,nu,iSpin,1) * SY(nu,mu,iSpin,2)&
                    & + SY(nu,mu,iSpin,2) * SY(mu,nu,iSpin,1) )
                tmprs = tmprs - 0.5_dp * &
                    & ( YS(mu,nu,iSpin,2) * YS(nu,mu,iSpin,1)&
                    & + YS(nu,mu,iSpin,1) * YS(mu,nu,iSpin,2)&
                    & + SYS(mu,nu,iSpin,2) * xmyas(nu,mu,iSpin,1)&
                    & + SYS(nu,mu,iSpin,1) * xmyas(mu,nu,iSpin,2)&
                    & + xmyas(mu,nu,iSpin,2) * SYS(nu,mu,iSpin,1)&
                    & + xmyas(nu,mu,iSpin,1) * SYS(mu,nu,iSpin,2)&
                    & + SY(mu,nu,iSpin,2) * SY(nu,mu,iSpin,1)&
                    & + SY(nu,mu,iSpin,1) * SY(mu,nu,iSpin,2) )
              end do
            end do
          end do
          ! Factor of two for spin-polarized calculation
          tmprs = cExchange * nSpin * tmprs

          nacv(:,iAt1) = nacv(:,iAt1) - 0.125_dp * tmprs * gammaLongRangePrime(:,iAt1,iAt2)
          nacv(:,iAt2) = nacv(:,iAt2) + 0.125_dp * tmprs * gammaLongRangePrime(:,iAt1,iAt2)
        end if

        call derivator%getFirstDeriv(dH0, skHamCont, coord0, species0, iAt1, iAt2, orb)
        call derivator%getFirstDeriv(dSo, skOverCont, coord0, species0, iAt1, iAt2, orb)

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
                ! tmp6 generalization could be wrong
                tmp6 = tmp6 + 0.5_dp * dSo(n,m,xyz) * (tmp7(iSpin,2) * xpycc(mu,nu,iSpin,1) +&
                         & tmp7(iSpin,1) * xpycc(mu,nu,iSpin,2))

                if (tSpin) then
                  tmp8 = tmp8 + tmp9 * dSo(n,m,xyz) * dsigma(iSpin) * pc(mu,nu,iSpin)
                  tmp10 = tmp10 + tmp11(iSpin) * dSo(n,m,xyz) * dsigma(iSpin) * SpinDens(mu,nu)
                end if

                if (tHybridXc) then
                  tmprs = 0.0_dp
                  do ka = 1, nOrb
                    tmprs = tmprs + (lrGammaOrb(mu,ka) + lrGammaOrb(nu,ka)) *&
                        & ( PS(mu,ka,iSpin) * deltaRho(nu,ka,iSpin)&
                        & + PS(nu,ka,iSpin) * deltaRho(mu,ka,iSpin)&
                        & + pc(mu,ka,iSpin) * DS(nu,ka,iSpin)&
                        & + pc(nu,ka,iSpin) * DS(mu,ka,iSpin) )
                    tmprs = tmprs + 0.5_dp * (lrGammaOrb(mu,ka) + lrGammaOrb(nu,ka)) *&
                        & ( xpyas(mu,ka,iSpin,1) * XS(nu,ka,iSpin,2)&
                        & + xpyas(ka,mu,iSpin,2) * SX(ka,nu,iSpin,1)&
                        & + xpyas(nu,ka,iSpin,1) * XS(mu,ka,iSpin,2)&
                        & + xpyas(ka,nu,iSpin,2) * SX(ka,mu,iSpin,1) )
                    tmprs = tmprs + 0.5_dp * (lrGammaOrb(mu,ka) + lrGammaOrb(nu,ka)) *&
                        & ( xpyas(mu,ka,iSpin,2) * XS(nu,ka,iSpin,1)&
                        & + xpyas(ka,mu,iSpin,1) * SX(ka,nu,iSpin,2)&
                        & + xpyas(nu,ka,iSpin,2) * XS(mu,ka,iSpin,1)&
                        & + xpyas(ka,nu,iSpin,1) * SX(ka,mu,iSpin,2) )
                    tmprs = tmprs + 0.5_dp * (lrGammaOrb(mu,ka) + lrGammaOrb(nu,ka)) *&
                        & ( xmyas(mu,ka,iSpin,1) * YS(nu,ka,iSpin,2)&
                        & + xmyas(ka,mu,iSpin,2) * SY(ka,nu,iSpin,1)&
                        & + xmyas(nu,ka,iSpin,1) * YS(mu,ka,iSpin,2)&
                        & + xmyas(ka,nu,iSpin,2) * SY(ka,mu,iSpin,1) )
                    tmprs = tmprs + 0.5_dp * (lrGammaOrb(mu,ka) + lrGammaOrb(nu,ka)) *&
                        & ( xmyas(mu,ka,iSpin,2) * YS(nu,ka,iSpin,1)&
                        & + xmyas(ka,mu,iSpin,1) * SY(ka,nu,iSpin,2)&
                        & + xmyas(nu,ka,iSpin,2) * YS(mu,ka,iSpin,1)&
                        & + xmyas(ka,nu,iSpin,1) * SY(ka,mu,iSpin,2) )
                    tmprs = tmprs + 0.5_dp * (lrGammaOrb(mu,ka) + lrGammaOrb(nu,ka)) *&
                        & ( XS(mu,ka,iSpin,1) * xpyas(ka,nu,iSpin,2)&
                        & + XS(nu,ka,iSpin,1) * xpyas(ka,mu,iSpin,2)&
                        & + xpyas(mu,ka,iSpin,1) * SX(ka,nu,iSpin,2)&
                        & + xpyas(nu,ka,iSpin,1) * SX(ka,mu,iSpin,2))
                    tmprs = tmprs + 0.5_dp * (lrGammaOrb(mu,ka) + lrGammaOrb(nu,ka)) *&
                        & ( XS(mu,ka,iSpin,2) * xpyas(ka,nu,iSpin,1)&
                        & + XS(nu,ka,iSpin,2) * xpyas(ka,mu,iSpin,1)&
                        & + xpyas(mu,ka,iSpin,2) * SX(ka,nu,iSpin,1)&
                        & + xpyas(nu,ka,iSpin,2) * SX(ka,mu,iSpin,1))
                    tmprs = tmprs - 0.5_dp * (lrGammaOrb(mu,ka) + lrGammaOrb(nu,ka)) *&
                        & ( YS(mu,ka,iSpin,1) * xmyas(ka,nu,iSpin,2)&
                        & + YS(nu,ka,iSpin,1) * xmyas(ka,mu,iSpin,2)&
                        & + xmyas(mu,ka,iSpin,1) * SY(ka,nu,iSpin,2)&
                        & + xmyas(nu,ka,iSpin,1) * SY(ka,mu,iSpin,2))
                    tmprs = tmprs - 0.5_dp * (lrGammaOrb(mu,ka) + lrGammaOrb(nu,ka)) *&
                        & ( YS(mu,ka,iSpin,2) * xmyas(ka,nu,iSpin,1)&
                        & + YS(nu,ka,iSpin,2) * xmyas(ka,mu,iSpin,1)&
                        & + xmyas(mu,ka,iSpin,2) * SY(ka,nu,iSpin,1)&
                        & + xmyas(nu,ka,iSpin,2) * SY(ka,mu,iSpin,1))
                  end do
                  ! Factor of 2 for spin-polarized calculations
                  tmprs2 = tmprs2 + cExchange * nSpin * dSo(n,m,xyz) * tmprs
                end if

              end do
            end do

          end do
          nacv(xyz,iAt1) = nacv(xyz,iAt1)&
              & + tmp1 + tmp2 + tmp4 + tmp6 + tmp3 + tmp8 + tmp10 - 0.25_dp * tmprs2
          nacv(xyz,iAt2) = nacv(xyz,iAt2)&
              & - tmp1 - tmp2 - tmp4 - tmp6 - tmp3 - tmp8 - tmp10 + 0.25_dp * tmprs2
        end do
      end do
    end do

  end subroutine addNadiaGradients


  !> Computes full range part of H^+_pq [X+-Y] as defined in Furche JCP 117 7433 (2002) eq. 20
  !! Here p/q are both virtual or both occupied orbitals and V is either X+Y or X-Y
  !! Note: The full range part of H^- is zero! 
  !! Note: This routine is specific for X+Y, for other quantities factors of 2 arise
  subroutine getHplusXYfr(sym, nXoo, nXvv, nAtom, getIA, getIJ, getAB, win, env, denseDesc, &
      & species0, ovrXev, grndEigVecs, frGamma, spinW, transChrg, XorY, vecHoo, vecHvv)

    !> Symmetry of the transition
    character, intent(in) :: sym

    !> Number of occ-occ transitions per spin channel
    integer, intent(in) :: nXoo(:)

    !> Number of vir-vir transitions per spin channel
    integer, intent(in) :: nXvv(:)

    !> Number of atoms
    integer, intent(in) :: nAtom

    !> Index array for occ-vir single particle excitations
    integer, intent(in) :: getIA(:,:)

    !> Index array for occ-occ single particle excitations
    integer, intent(in) :: getIJ(:,:)

    !> Index array for vir-vir single particle excitations
    integer, intent(in) :: getAB(:,:)

    !> Index array for single particle excitations
    integer, intent(in) :: win(:)

    !> Environment settings
    type(TEnvironment), intent(inout) :: env

    !> Dense matrix descriptor
    type(TDenseDescr), intent(in) :: denseDesc 

    !> Central cell chemical species
    integer, intent(in) :: species0(:)

    !> Overlap times ground state eigenvectors
    real(dp), intent(in) :: ovrXev(:,:,:)

    !> Ground state eigenvectors
    real(dp), intent(in) :: grndEigVecs(:,:,:)

    !> Electrostatic matrix, full-range
    real(dp), intent(in) :: frGamma(:,:)

    !> Ground state spin derivatives for each species
    real(dp), intent(in) :: spinW(:)

    !> Machinery for transition charges between single particle levels
    type(TTransCharges), intent(in) :: transChrg

    !> RPA eigenvectors, either (X+Y) or (X-Y)
    real(dp), intent(in) :: XorY(:)

    !> Output vector H[V] occupied-occupied
    real(dp), optional, intent(out) :: vecHoo(:)

    !> Output vector H[V] virtual-virtual
    real(dp), optional, intent(out) :: vecHvv(:)

    integer :: nSpin, ab, s, abs, svv(2), ij, ijs, soo(2)
    real(dp) :: fact
    real(dp), allocatable  :: xpyq(:), gamxpyq(:), qTr(:), xpyqds(:), gamxpyqds(:)
    logical :: tSpin

    nSpin = size(grndEigVecs, dim=3)
    tSpin = (nSpin == 2)
    allocate(xpyq(nAtom))
    allocate(qTr(nAtom))
    allocate(gamxpyq(nAtom))
    soo(:) = [0, nXoo(1)]
    svv(:) = [0, nXvv(1)]

    if (present(vecHoo)) then
      vecHoo(:) = 0.0_dp
    end if
    if (present(vecHvv)) then
      vecHvv(:) = 0.0_dp
    end if

    xpyq(:) = 0.0_dp
    call transChrg%qMatVec(env, denseDesc, ovrXev, grndEigVecs, getIA, win, XorY, xpyq)

    if (.not. tSpin) then  ! ---- spin-unpolarized case ----
      ! vecHvv(ab) = sum_jc K_ab,jc (X+Y)_jc
      if (sym == "S") then
        call hemv(gamxpyq, frGamma,  xpyq)
        if (present(vecHvv)) then
          do ab = 1, nXvv(1)
            qTr(:) = transChrg%qTransAB(ab, env, denseDesc, ovrXev, grndEigVecs, getAB)
            vecHvv(ab) = 2.0_dp * sum(qTr * gamxpyq)
          end do
        end if
        if (present(vecHoo)) then
          do ij = 1, nXoo(1)
            qTr(:) = transChrg%qTransIJ(ij, env, denseDesc, ovrXev, grndEigVecs, getIJ)
            ! vecHoo(ij) = sum_kb K_ij,kb (X+Y)_kb
            vecHoo(ij) = 2.0_dp * sum(qTr * gamxpyq)
          end do
        end if
      else ! Triplet case
        if (present(vecHvv)) then
          do ab = 1, nXvv(1)
            qTr(:) = transChrg%qTransAB(ab, env, denseDesc, ovrXev, grndEigVecs, getAB)
            vecHvv(ab) = 2.0_dp * sum(qTr * xpyq * spinW(species0))
          end do
        end if
        if (present(vecHoo)) then
          do ij = 1, nXoo(1)
            qTr(:) = transChrg%qTransIJ(ij, env, denseDesc, ovrXev, grndEigVecs, getIJ)
            vecHoo(ij) = 2.0_dp * sum(qTr * xpyq * spinW(species0))
          end do
        end if
      end if

    else  ! ---- spin-polarized case -----

      allocate(xpyqds(nAtom))
      allocate(gamxpyqds(nAtom))
      xpyqds(:) = 0.0_dp
      call transChrg%qMatVecDs(env, denseDesc, ovrXev, grndEigVecs, getIA, win, XorY, xpyqds)

      call hemv(gamxpyq, frGamma,  xpyq)
      do s = 1, 2
        if (s == 1) then
          fact = 1.0_dp
        else
          fact = -1.0_dp
        end if
        if (present(vecHvv)) then
          do ab = 1, nXvv(s)
            abs = ab + svv(s)
            qTr(:) = transChrg%qTransAB(abs, env, denseDesc, ovrXev, grndEigVecs, getAB)
            vecHvv(abs) = sum(qTr * gamxpyq)
            ! Magnetization part
            vecHvv(abs) = vecHvv(abs) + fact * sum(qTr * xpyqds * spinW(species0))
          end do
        end if
        if (present(vecHoo)) then
          do ij = 1, nXoo(s)
            ijs = ij + soo(s)
            qTr(:) = transChrg%qTransIJ(ijs, env, denseDesc, ovrXev, grndEigVecs, getIJ)
            vecHoo(ijs) = sum(qTr * gamxpyq)
            !magnetization part
            vecHoo(ijs) = vecHoo(ijs) + fact * sum(qTr * xpyqds * spinW(species0))
          end do
        end if
      end do

    end if

  end subroutine getHplusXYfr


  !> Computes full range part of H^+_pq [M] as defined in Furche JCP 117 7433 (2002) eq. 20
  !! Here pq are arbitrary orbitals and M is a general matrix with ov,oo,vv components
  !! iMode = 1: returns oo components of H
  !! iMode = 2: returns vv components of H
  !! iMode = 3: returns ov components of H
  !! Note: The full range part of H^- is zero! 
  !! Routine currently does not work for M_ia /= 0 on entry!
  subroutine getHplusMfr(iMode, nXoo, nXvv, nXov, nAtom, getIA, getIJ, getAB, win, env, denseDesc, &
      & species0, ovrXev, grndEigVecs, frGamma, spinW, transChrg, matM, vecH)

    !> Type of return vector (oo, vv, ov)
    integer , intent(in)  :: iMode

    !> Number of occ-occ transitions per spin channel
    integer, intent(in) :: nXoo(:)

    !> Number of vir-vir transitions per spin channel
    integer, intent(in) :: nXvv(:)

    !> Number of occ-vir transitions (both spins)
    integer, intent(in) :: nXov

    !> Number of atoms
    integer, intent(in) :: nAtom

    !> Index array for occ-vir single particle excitations
    integer, intent(in) :: getIA(:,:)

    !> Index array for occ-occ single particle excitations
    integer, intent(in) :: getIJ(:,:)

    !> Index array for vir-vir single particle excitations
    integer, intent(in) :: getAB(:,:)

    !> Index array for single particle excitations
    integer, intent(in) :: win(:)

    !> Environment settings
    type(TEnvironment), intent(inout) :: env

    !> Dense matrix descriptor
    type(TDenseDescr), intent(in) :: denseDesc
    
    !> central cell chemical species
    integer, intent(in) :: species0(:)

    !> Overlap times ground state eigenvectors
    real(dp), intent(in) :: ovrXev(:,:,:)

    !> Ground state eigenvectors
    real(dp), intent(in) :: grndEigVecs(:,:,:)

    !> Electrostatic matrix, full-range
    real(dp), intent(in) :: frGamma(:,:)

    !> Ground state spin derivatives for each species
    real(dp), intent(in) :: spinW(:)

    !> Machinery for transition charges between single particle levels
    type(TTransCharges), intent(in) :: transChrg

    !> Input Matrix spin-resolved
    real(dp), intent(in) :: matM(:,:,:)

    !> Output vector H[M]
    real(dp), intent(out) :: vecH(:)

    integer :: nSpin, ab, i, j, a, b, s, svv(2), ij, soo(2), ias
    real(dp), dimension(2) :: spinFactor = [1.0_dp, -1.0_dp]
    real(dp), allocatable  :: xpyq(:), gamxpyq(:), qTr(:), xpyqds(:), gamxpyqds(:)
    real(dp), allocatable  :: gamqt(:)
    logical :: tSpin

    if (iMode == 1) then
      @:ASSERT(size(vecH) == sum(nXoo))
    else if (iMode == 2) then
      @:ASSERT(size(vecH) == sum(nXvv))
    else
      @:ASSERT(size(vecH) == nXov)
    end if

    nSpin = size(grndEigVecs, dim=3)
    tSpin = (nSpin == 2)
    if (tSpin) then
      allocate(xpyqds(natom))
      allocate(gamxpyqds(natom))
    endif
    allocate(xpyq(nAtom))
    allocate(qTr(nAtom))
    allocate(gamxpyq(nAtom))
    allocate(gamqt(nAtom))
    soo(:) = [0, nXoo(1)]
    svv(:) = [0, nXvv(1)]

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
        qTr(:) = transChrg%qTransIJ(ij + soo(s), env, denseDesc, ovrXev, grndEigVecs, getIJ)
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
        qTr(:) = transChrg%qTransAB(ab + svv(s), env, denseDesc, ovrXev, grndEigVecs, getAB)
        if (a == b) then
          gamxpyq(:) = gamxpyq(:) + matM(a,b,s) * qTr(:)
          if (tSpin) then
            gamxpyqds(:) = gamxpyqds(:) + matM(a,b,s) * qTr(:) * spinFactor(s)
          end if
        else
          ! Factor 2 because of symmetry of the matrix
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
      qTr(:) = transChrg%qTransIA(ias, env, denseDesc, ovrXev, grndEigVecs, getIA, win)
      gamxpyq(:) = gamxpyq(:) + (matM(i,a,s) + matM(a,i,s)) * qTr(:)
      if (tSpin) then
         gamxpyqds(:) = gamxpyqds(:) + (matM(i,a,s) + matM(a,i,s)) * qTr(:) * spinFactor(s)
      end if
    end do

    ! gamqt(iAt1) = sum_iAt2 gamma_iAt1,iAt2 gamxpyq(iAt2)
    call hemv(gamqt, frGamma, gamxpyq)

    if (iMode == 1) then
      do s = 1, nSpin
        do ij = 1, nxoo(s)
          i = getIJ(ij + soo(s), 1)
          j = getIJ(ij + soo(s), 2)
          qTr(:) = transChrg%qTransIJ(ij + soo(s), env, denseDesc, ovrXev, grndEigVecs, getIJ)
          if (.not. tSpin) then
            vecH(ij + soo(s)) = 4.0_dp * dot_product(gamqt,qTr)
          else
            vecH(ij + soo(s)) = 2.0_dp * dot_product(gamqt,qTr)
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
          qTr(:) = transChrg%qTransAB(ab + svv(s), env, denseDesc, ovrXev, grndEigVecs, getAB)
          if (.not. tSpin) then
            vecH(ab + svv(s)) = 4.0_dp * dot_product(gamqt,qTr)
          else
            vecH(ab + svv(s)) = 2.0_dp * dot_product(gamqt,qTr)
            vecH(ab + svv(s)) = vecH(ab + svv(s)) + &
                & spinFactor(s) * dot_product(gamxpyqds * spinW(species0), qTr)
          end if
        end do
      end do
    else
      do ias = 1, nxov
        i = getIA(ias, 1)
        a = getIA(ias, 2)
        s = getIA(ias, 3)
        qTr(:) = transChrg%qTransIA(ias, env, denseDesc, ovrXev, grndEigVecs, getIA, win)
        if (.not. tSpin) then
           vecH(ias) = 4.0_dp * dot_product(gamqt,qTr)
        else
          vecH(ias) =  2.0_dp * dot_product(gamqt,qTr)
          vecH(ias) = vecH(ias) + &
              & spinFactor(s) * dot_product(gamxpyqds * spinW(species0), qTr)
        end if
      end do
    end if

  end subroutine getHplusMfr


  !> Convention to fix phase of NAC vector: Largest value is positive. If several entries are closer
  !! than nacTol to the maximum, the lowest index is chosen.
  subroutine fixNACVPhase(nacv)

    !> Non-adiabatic coupling matrix for (3, atoms, transitions)
    real(dp), intent(inout) :: nacv(:,:,:)

    real(dp) :: max
    integer :: ii, jj, iNac, nAtoms, phase

    nAtoms = size(nacv, dim=2)

    do iNac = 1, size(nacv, dim=3)
      max = maxval(abs(nacv(:,:,iNac)))
      outer: do ii = 1, nAtoms
        do jj = 1, 3
          if (abs(abs(nacv(jj,ii,iNac)) - max) < nacTol) exit outer
        end do
      end do outer

      phase = int(nacv(jj, ii, iNac) / abs(nacv(jj, ii, iNac)))
      nacv(:,:,iNac) = phase * nacv(:,:,iNac)
    end do

  end subroutine fixNACVPhase


  !> Implements the CI optimizer of Bearpark et al. Chem. Phys. Lett. 223 269 (1994) with
  !! modifications introduced by Harabuchi/Hatanaka/Maeda CPL X 2019.
  !! Previous published results [Niehaus JCP 158 054103 (2023), TCA 140 34 (2021)] were obtained
  !! with a differing version that assumed orthogonal X1 and X2 vectors, which leads to poor
  !! convergence.
  subroutine conicalIntersectionOptimizer(derivs, excDerivs, indNACouplings, energyShift,&
      & naCouplings, excEnergies)

    !> Ground state gradient (overwritten)
    real(dp), intent(inout) :: derivs(:,:)

    !> Gradients of the excitation energy
    real(dp), intent(in) :: excDerivs(:,:,:)

    !> States between which CI is optimized
    integer, intent(in) :: indNACouplings(2)

    !> Shift of excited state PES (Harabuchi/Hatanaka/Maeda CPL X 2019)
    real(dp), intent(in) :: energyShift

    !> Nonadiabatic coupling vectors
    real(dp), intent(in) :: naCouplings(:,:,:)

    !> Sn-S0 excitation energy
    real(dp), intent(in) :: excEnergies(:)

    integer :: nAtoms, nexcGrad, nCoupl
    real(dp), allocatable :: X1(:), X2(:), dE2(:), gpf(:)
    real(dp) :: dpX1, dpX2, dp12, prj1, prj2, deltaE, normGP
    real(dp) :: alpa, beta
    character(len=*), parameter :: format2U = "(A, ':', T32, F18.10, T51, A, T54, F16.4, T71, A)"

    nAtoms = size(derivs, dim=2)
    nexcGrad = indNACouplings(2) - indNACouplings(1) + 1
    nCoupl = nexcGrad * (nexcGrad - 1) / 2
    if (indNACouplings(1) == 0) then
      nexcGrad = nexcGrad - 1
    end if
    @:ASSERT(nexcGrad == size(excDerivs, dim=3))
    @:ASSERT(nCoupl == size(naCouplings, dim=3))

    allocate(X1(3 * nAtoms))
    allocate(X2(3 * nAtoms))
    allocate(dE2(3 * nAtoms))
    allocate(gpf(3 * nAtoms))

    if (indNACouplings(1) == 0) then
      X1 = reshape(excDerivs(:,:,nexcGrad), [3 * nAtoms])
    else
      X1 = reshape(excDerivs(:,:,nexcGrad) - excDerivs(:,:,1), [3 * nAtoms])
    end if
    ! Last entry of array holds coupling between states indNACouplings(1/2)
    X2 = reshape(naCouplings(:,:,nCoupl), [3 * nAtoms])

    ! Original Bearbark suggestion would be dS1/dq, Harabuchi chooses 0.5 (dS0/dq + dS1/dq)
    dE2 = reshape(derivs + 0.5_dp * (excDerivs(:,:,nexcGrad) + excDerivs(:,:,1)), [3 * nAtoms])

    dpX1 = dot_product(X1, X1)
    dpX2 = dot_product(X2, X2)
    dp12 = dot_product(X1, X2)
    prj1 = dot_product(dE2, X1)
    prj2 = dot_product(dE2, X2)

    alpa = (prj2 * dp12 - prj1 * dpX2) / (dpX1 * dpX2 - dp12 * dp12)
    beta = (prj2 * dpX1 - prj1 * dp12) / (dp12 * dp12 - dpX1 * dpX2)

    ! Eq.(5) in Bearpark et al.
    gpf(:) = dE2 + alpa * X1 + beta * X2
    normGP = norm2(gpf)

    ! Eq.(4) in Bearpark et al. with modifications by Harabuchi
    ! Yields approximate CI without running into SCF problems too early
    ! Shift should be brought to zero
    if (indNACouplings(1) == 0) then
      deltaE = excEnergies(indNACouplings(2))
    else
      deltaE = excEnergies(indNACouplings(2)) - excEnergies(indNACouplings(1))
    end if
    gpf(:) = gpf + 2.0_dp * (deltaE - energyShift) * X1 / sqrt(dpX1)

    derivs(:,:) = reshape(gpf, [3, nAtoms])

    write(stdOut, format2U) "Energy gap CI", deltaE, 'H', Hartree__eV * deltaE, 'eV'

  end subroutine conicalIntersectionOptimizer

end module dftbp_timedep_linrespgrad
