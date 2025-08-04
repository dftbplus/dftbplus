!--------------------------------------------------------------------------------------------------!
!  DFTB+: general package for performing fast atomistic simulations                                !
!  Copyright (C) 2006 - 2025  DFTB+ developers group                                               !
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
  use dftbp_extlibs_arpack, only : psaupd, pseupd, saupd, seupd, withArpack
  use dftbp_io_message, only : error
  use dftbp_io_taggedoutput, only : TTaggedWriter, tagLabels
  use dftbp_math_blasroutines, only : gemm, hemv, symm, herk
  use dftbp_math_degeneracy, only : TDegeneracyFind
  use dftbp_math_eigensolver, only : heev
  use dftbp_math_matrixops, only : orthonormalizeVectors
  use dftbp_math_qm, only : makeSimilarityTrans
  use dftbp_math_sorting, only : index_heap_sort, merge_sort
  use dftbp_timedep_linrespcommon, only : excitedDipoleOut, excitedQOut, twothird,&
      & oscillatorStrength, indxoo, indxov, indxvv, rindxov_array,&
      & getSPExcitations, calcTransitionDipoles, dipselect, transitionDipole, writeSPExcitations,&
      & getExcSpin, writeExcMulliken, actionAplusB, actionAminusB, initialSubSpaceMatrixApmB,&
      & calcMatrixSqrt, incMemStratmann, getSqrOcc
  use dftbp_timedep_linresptypes, only : TLinResp, linrespSolverTypes, TCasidaParameter,&
      & TCasidaParameter_init
  use dftbp_timedep_transcharges, only : TTransCharges, transq, TTransCharges_init
  use dftbp_type_commontypes, only : TOrbitals
  use dftbp_type_densedescr, only : TDenseDescr
  use dftbp_common_environment, only : TEnvironment, globalTimers

#:if WITH_SCALAPACK

  use dftbp_extlibs_scalapackfx, only : pblasfx_psymm
  use dftbp_extlibs_mpifx, only : MPI_SUM, mpifx_allreduceip, mpifx_bcast
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
  character(*), parameter :: transChrgOut = "ATQ.DAT"


  ! Solver related variables

  !> Tolerance for ARPACK solver.
  real(dp), parameter :: arTol = epsilon(1.0_rsp)

  !> Threshold for Stratmann solver
  real(dp), parameter :: convThreshStrat = 0.1*epsilon(1.0_rsp)

  !> Maximal allowed iteration in the ARPACK solver.
  integer, parameter :: maxArIter = 300

  !> Names of output files
  character(*), parameter :: arpackOut = "ARPACK.DAT"
  character(*), parameter :: testArpackOut = "TEST_ARPACK.DAT"

  !> Threshold for near-identical NACV values to fix phase
  real(dp), parameter :: nacTol = 1.0e-6_dp

contains

  !> This subroutine analytically calculates excitations and gradients of excited state energies
  !! based on Time Dependent DFRT
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
    integer, allocatable :: win(:), getIA(:,:), getIJ(:,:), getAB(:,:)

    !> MPI Global array
    real(dp), allocatable :: VecGlb(:,:)

    !> Array from pairs of single particles states to compound index - should replace with a more
    !> compact data structure in the cases where there are oscilator windows
    integer, allocatable :: iatrans(:,:,:)

    character, allocatable :: symmetries(:)

    integer :: nxoo_max, nxvv_max
    integer, allocatable :: nocc_ud(:), nvir_ud(:), nxoo_ud(:), nxvv_ud(:), nxov_ud(:)
    integer :: mHOMO, mLUMO
    integer :: nxov, nxov_r, nxov_d, nxov_rd
    integer :: norb, nxoo, nxvv
    integer :: i, j, iSpin, isym, iLev, iSav, nStartLev, nEndLev
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

    !> Casida parameters (number of transitions, index arrays and alike)
    type(TCasidaParameter) :: rpa

    type(TFileDescr) :: fdTrans, fdTransDip, fdArnoldi, fdXPlusY, fdExc, fdTransQ

    !> Communication with ARPACK for progress information
    integer :: logfil, ndigit, mgetv0
    integer :: msaupd, msaup2, msaitr, mseigt, msapps, msgets, mseupd
    integer :: mnaupd, mnaup2, mnaitr, mneigh, mnapps, mngets, mneupd
    integer :: mcaupd, mcaup2, mcaitr, mceigh, mcapps, mcgets, mceupd

    !> Common block of ARPACK variables
    common /debug/ logfil, ndigit, mgetv0,&
        & msaupd, msaup2, msaitr, mseigt, msapps, msgets, mseupd,&
        & mnaupd, mnaup2, mnaitr, mneigh, mnapps, mngets, mneupd,&
        & mcaupd, mcaup2, mcaitr, mceigh, mcapps, mcgets, mceupd

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
    allocate(nxov_ud(nSpin))
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
    nocc_ud(:) = 0
    nvir_ud(:) = 0
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
    allocate(nxoo_ud(nSpin))
    allocate(nxvv_ud(nSpin))
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

    ! Is a z vector required?
    tZVector = tForces .or. this%writeMulliken .or. this%writeCoeffs .or. present(naturalOrbs) .or.&
        & this%tWriteDensityMatrix .or. this%tNaCoupling
    doAllZVectors = tZVector .and. (nstat == 0) .and. (.not. this%isCIopt) .and.&
        & (.not. this%tNaCoupling)

    ! Occ-occ/vir-vir charges only required for Z-vector/forces or TD-LC-DFTB
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
      call pblasfx_psymm(SSqr, denseDesc%blacsOrbSqr, grndEigVecs(:,:,iSpin),&
          & denseDesc%blacsOrbSqr, ovrXev(:,:,iSpin), denseDesc%blacsOrbSqr, side="L")
    end do

    call env%globalTimer%startTimer(globalTimers%lrCoulomb)
    call sccCalc%getAtomicGammaMatrixBlacs(gammaMat, iNeighbour, img2CentCell, env)
    call env%globalTimer%stopTimer(globalTimers%lrCoulomb)

  #:else

    do iSpin = 1, nSpin
      call symm(ovrXev(:,:,iSpin), "L", SSqr, grndEigVecs(:,:,iSpin))
    end do

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
    wij(:) = wij(win)

    ! Build square root of occupation difference between virtual and occupied states
    call getSqrOcc(filling, win, nxov_ud(1), nxov, getIA, this%tSpin, sqrOccIA)

    call env%globalTimer%startTimer(globalTimers%lrTransCharges)

    ! First call to initialize charges for all occ-vir transitions
    call TTransCharges_init(transChrg, env, denseDesc, ovrXev, grndEigVecs, norb, nxov,&
        & nxov_ud(1), nxoo_ud, nxvv_ud, getIA, getIJ, getAB, win, this%tCacheChargesOccVir,&
        & this%tCacheChargesSame, .true.)

    ! dipole strength of transitions between K-S states
    call calcTransitionDipoles(coord0, win, getIA, transChrg, env, denseDesc, ovrXev, grndEigVecs,&
        & snglPartTransDip)

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

    ! set up transition indexing
    allocate(iatrans(norb, norb, nSpin))
    call rindxov_array(win, nxov, nxoo, nxvv, getIA, getIJ, getAB, iatrans)

    ! MPI distribution of RPA vectors according to these indices
    call distributeRangeInChunks(env, 1, nxov_rd, iGlobal, fGlobal)

    ! All relevant run time parameters of Casida are stored in a derived type
    ! Input arrays are deallocated on return
    call TCasidaParameter_init(rpa, nocc_ud, nvir_ud, nxoo_ud, nxvv_ud, nxov_ud, nxov_rd,&
        & iaTrans, getIA, getIJ, getAB, win, wij, sqrOccIA, tHybridXc, tZVector)

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
    if (this%writeExcitations) then
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
    end if

    ! single particle excitations output file
    call writeSPExcitations(this, rpa, sposz)

    allocate(xpy(rpa%nxov_rd, this%nExc))
    if (tZVector .or. tHybridXc) then
      allocate(xmy(rpa%nxov_rd, this%nExc))
    end if

    if (this%iLinRespSolver /= linrespSolverTypes%stratmann .and. tHybridXc) then
      call error("Range separation requires the Stratmann solver for excitations")
    end if

    call env%globalTimer%stopTimer(globalTimers%lrSetup)

    do isym = 1, size(symmetries)

      sym = symmetries(isym)
      call env%globalTimer%startTimer(globalTimers%lrSolver)
      select case (this%iLinRespSolver)
      case (linrespSolverTypes%arpack)
        call buildAndDiagExcMatrixArpack(iGlobal, fGlobal, env, orb, this, rpa, transChrg,&
            & denseDesc, ovrXev, grndEigVecs, gammaMat, species0, eval, sym, xpy, xmy)

      case (linrespSolverTypes%stratmann)
        call buildAndDiagExcMatrixStratmann(iGlobal, fGlobal, env, orb, this, rpa, transChrg,&
            & denseDesc, ovrXev, grndEigVecs, gammaMat, lrGamma, species0, eval, sym, xpy, xmy)
      end select
      call env%globalTimer%stopTimer(globalTimers%lrSolver)

      ! Excitation oscillator strengths for resulting states
      call getOscillatorStrengths(this, rpa, sym, eval, xpy, snglPartTransDip(1:rpa%nxov_rd,:),&
          & nstat, osz, transitionDipoles)

      ! Transition charges for state nstat
      if (this%writeTransQ) then
        call getAndWriteTransitionCharges(env, this, rpa, transChrg, sym, denseDesc, ovrXev,&
            & grndEigVecs, xpy, fdTransQ, fdTagged, taggedWriter)
      end if

      if (this%tSpin) then

        call getExcSpin(env, orb, rpa, denseDesc, Ssq, xpy, filling, ovrXev, grndEigVecs)

        call writeExcitations(this, rpa, sym, osz, eval, xpy, fdXPlusY, fdTrans, fdTransDip,&
            & transitionDipoles, fdTagged, taggedWriter, fdExc, Ssq)

      else

        call writeExcitations(this, rpa, sym, osz, eval, xpy, fdXPlusY, fdTrans, fdTransDip,&
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

      call env%globalTimer%startTimer(globalTimers%lrZVector)

      ! Differentiates between a standard Z-vector equation for transition densities and forces and
      ! a specific one for non-adiabatic couplings
      doVanillaZvector = .true.
      if (doAllZVectors) then

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
          call getZVectorEqRHS(env, orb, this, rpa, transChrg, sym, denseDesc, species0, grndEigVal,&
            & ovrXev, grndEigVecs, gammaMat, lrGamma, omega, xpy(:,iLev), xmy(:,iLev), rhs, t, wov,&
            & woo, wvv)

          call solveZVectorPrecond(env, orb, this, rpa, transChrg, denseDesc, species0, ovrXev,&
            & grndEigVecs, gammaMat, lrGamma, rhs)

          call calcWVectorZ(env, orb, this, rpa, transChrg, denseDesc, species0, ovrXev, grndEigVecs,&
            & grndEigVal, gammaMat, lrGamma, rhs, wov, woo, wvv)

          call calcPMatrix(env, rpa, t, rhs, pc)

          call writeCoeffs(pc, grndEigVecs, filling, this%writeCoeffs, this%tGrndState, occNatural,&
            & naturalOrbs)

          do iSpin = 1, nSpin
            ! Make MO to AO transformation of the excited density matrix
            allocate(VecGlb(norb,norb))
            call distrib2replicated(env%blacs%orbitalGrid, denseDesc%blacsOrbSqr, grndEigVecs(:,:,iSpin),&
                                    & VecGlb(:,:))

            call makeSimilarityTrans(pc(:,:,iSpin), VecGlb(:,:))

            call distrib2replicated(env%blacs%orbitalGrid, denseDesc%blacsOrbSqr, SSqr,&
                                    & VecGlb)

            call getExcMulliken(denseDesc, pc(:,:,iSpin), VecGlb, dqex(:,iSpin))
            deallocate(VecGlb)
          end do

          if (this%tWriteDensityMatrix) then
            call writeDM(iLev, pc, rhoSqr)
          end if

          if (this%writeMulliken) then
            ! This prints charges for an excited state from the relaxed transition density
            call writeExcMulliken(sym, iLev, dq(:,1), sum(dqex,dim=2), coord0)
          end if

          if (tForces) then
            iSav = iLev - nStartLev + 1
            call addGradients(env, orb, this, rpa, transChrg, hybridXc, denseDesc, sym, species0,&
                & ovrXev, grndEigVecs, gammaMat, lrGamma, coord0, dq, dqex, shift, xpy(:,iLev),&
                & xmy(:,iLev), woo, wov, wvv, skHamCont, skOverCont, derivator, rhoSqr, pc,&
                & excgrad(:,:,iSav), deltaRho=deltaRho)
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

            ! Ground-to-excited NACV
            if (nCoupLev == 0) then
              xpym(:) = xpy(:,mCoupLev)
              xmym(:) = xmy(:,mCoupLev)
              omegaDif = sqrt(eval(mCoupLev))
              call grndToExcDensityMatrices(env, orb, this, rpa, transChrg, denseDesc, sym, species0,&
                  & ovrXev, grndEigVecs, grndEigVal, gammaMat, lrGamma, omegaDif, pc, xpym, xmym,&
                  & wov, woo)

              do iSpin = 1, nSpin
                ! Make MO to AO transformation of the excited density matrix
                allocate(VecGlb(norb,norb))
                call distrib2replicated(env%blacs%orbitalGrid, denseDesc%blacsOrbSqr, grndEigVecs(:,:,iSpin),&
                                        & VecGlb(:,:))

                call makeSimilarityTrans(pc(:,:,iSpin), VecGlb(:,:))

                call distrib2replicated(env%blacs%orbitalGrid, denseDesc%blacsOrbSqr, SSqr,&
                                        & VecGlb)

                call getExcMulliken(denseDesc, pc(:,:,iSpin), VecGlb, dqex(:,iSpin))
                deallocate(VecGlb)
              end do

              ! For 0-n couplings, the standard force routine can be used, where
              ! X+Y and W_ab are zeroed out
              wvv(:,:)  = 0.0_dp
              xpym(:) = 0.0_dp
              xmym(:) = 0.0_dp
              xpyn(:) = 0.0_dp
              xmyn(:) = 0.0_dp

              call addGradients(env, orb, this, rpa, transChrg, hybridXc, denseDesc, sym, species0,&
                  & ovrXev, grndEigVecs, gammaMat, lrGamma, coord0, dq, dqex, shift, xpym, xmym,&
                  & woo, wov, wvv, skHamCont, skOverCont, derivator, rhoSqr, pc, nacv(:,:,iNac),&
                  & deltaRho=deltaRho)

            else

              xpyn(:) = xpy(:,nCoupLev)
              xmyn(:) = xmy(:,nCoupLev)
              xpym(:) = xpy(:,mCoupLev)
              xmym(:) = xmy(:,mCoupLev)
              omegaDif = sqrt(eval(nCoupLev)) - sqrt(eval(mCoupLev))
              omegaAvg = 0.5_dp * (sqrt(eval(nCoupLev)) + sqrt(eval(mCoupLev)))

              ! compute + component of RHS for Z-vector eq. in the NaCoupling case
              ! also computes the + components of W and T
              call getNadiaZvectorEqRHS(env, orb, this, rpa, transChrg, sym, denseDesc, species0,&
                  & grndEigVal, ovrXev, grndEigVecs, gammaMat, lrGamma, omegaAvg, xpy(:,nCoupLev),&
                  & xmy(:,nCoupLev), xpy(:,mCoupLev), xmy(:,mCoupLev), rhs, t, wov, woo, wvv)

              call solveZVectorPrecond(env, orb, this, rpa, transChrg, denseDesc, species0, ovrXev,&
                  & grndEigVecs, gammaMat, lrGamma, rhs)

              call calcWVectorZ(env, orb, this, rpa, transChrg, denseDesc, species0, ovrXev,&
                  & grndEigVecs, grndEigVal, gammaMat, lrGamma, rhs, wov, woo, wvv)

              call calcPMatrix(env, rpa, t, rhs, pc)

              do iSpin = 1, nSpin
                ! Make MO to AO transformation of the excited density matrix
                allocate(VecGlb(norb,norb))
                call distrib2replicated(env%blacs%orbitalGrid, denseDesc%blacsOrbSqr, grndEigVecs(:,:,iSpin),&
                                        & VecGlb(:,:))

                call makeSimilarityTrans(pc(:,:,iSpin), VecGlb(:,:))

                call distrib2replicated(env%blacs%orbitalGrid, denseDesc%blacsOrbSqr, SSqr,&
                                        & VecGlb)

                call getExcMulliken(denseDesc, pc(:,:,iSpin), VecGlb, dqex(:,iSpin))
                deallocate(VecGlb)
              end do

              call addNadiaGradients(env, orb, this, rpa, transChrg, hybridXc, denseDesc, sym,&
                  & species0, ovrXev, grndEigVecs, gammaMat, lrGamma, coord0, dq, dqex, shift,&
                  & xpyn, xmyn, xpym, xmym, woo, wov, wvv, skHamCont, skOverCont, derivator,&
                  & rhoSqr, pc, nacv(:,:,iNac), deltaRho)

            end if

            ! P and W have not yet been divided by excitation energy difference
            nacv(:,:,iNac) = nacv(:,:,iNac) / omegaDif

          end do
        end do

        ! Convention to determine arbitrary phase
        call fixNACVPhase(nacv)

        if (this%writeNacv) then
          call writeNACV(this%indNACouplings(1), this%indNACouplings(2), fdTagged, taggedWriter,&
              & nacv)
        end if

        call env%globalTimer%stopTimer(globalTimers%lrNAC)

      end if

      call env%globalTimer%stopTimer(globalTimers%lrZVector)

    end if

    call env%globalTimer%stopTimer(globalTimers%lrZVector)

    ! Omega has possibly been overwritten for CI optimization or NA couplings, but should always
    ! refer to nstat
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
  subroutine buildAndDiagExcMatrixArpack(iGlobal, fGlobal, env, orb, lr, rpa, transChrg,&
      & denseDesc, ovrXev, grndEigVecs, gammaMat, species0, eval, sym, xpy, xmy)

    !> Starting index of current rank in global RPA vectors
    integer, intent(in) :: iGlobal

    !> End index of current rank in global RPA vectors
    integer, intent(in) :: fGlobal

    !> Environment settings
    type(TEnvironment), intent(inout) :: env

    !> Data type for atomic orbital information
    type(TOrbitals), intent(in) :: orb

    !> Data structure for linear response
    type(TLinResp), intent(in) :: lr

    !> Run time parameters of the Casida routine
    type(TCasidaParameter), intent(in) :: rpa

    !> Machinery for transition charges between single particle levels
    type(TTransCharges), intent(in) :: transChrg

    !> Dense matrix descriptor
    type(TDenseDescr), intent(in) :: denseDesc

    !> overlap times ground state eigenvectors (local in case of MPI)
    real(dp), intent(in) :: ovrXev(:,:,:)

    !> ground state eigenvectors (local in case of MPI)
    real(dp), intent(in) :: grndEigVecs(:,:,:)

    !> Electrostatic matrix
    real(dp), intent(in) :: gammaMat(:,:)

    !> Central cell chemical species
    integer, intent(in) :: species0(:)

    !> Resulting eigenvalues for transitions (w^2)
    real(dp), intent(out) :: eval(:)

    !> Symmetry to calculate transitions
    character, intent(in) :: sym

    !> Eigenvectors (X+Y)
    real(dp), intent(out) :: xpy(:,:)

    !> Eigenvectors (X-Y), only evaluated if Z-vector is needed
    real(dp), intent(inout), allocatable :: xmy(:,:)

    integer :: iparam(11), ipntr(11)
    integer :: ido, ncv, lworkl, info
    integer :: nexc, natom, nLoc
    integer :: iState, comm
    real(dp), allocatable :: workl(:), workd(:), resid(:), vv(:,:), qij(:)
    real(dp), allocatable :: Hv(:), orthnorm(:,:)
    real(dp) :: sigma, omega
    logical, allocatable :: selection(:)
    logical :: rvec
    character(lc) :: tmpStr
    type(TFileDescr) :: fdArnoldiTest

  #:if WITH_PARPACK
    comm = env%mpi%globalComm%id
  #:endif

    ! Local chunk of RPA vectors have this size under MPI
    nLoc = fGlobal - iGlobal + 1

    nexc = size(eval)
    natom = size(gammaMat, dim=1)

    @:ASSERT(all(shape(xpy) == [ rpa%nxov_rd, nexc ]))
    @:ASSERT(rpa%tHybridXc .eqv. .false.)

    ! Three times more Lanczos vectors than desired eigenstates
    ncv = min(3 * nexc, rpa%nxov_rd)

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
    iparam(3) = maxArIter
    ! solve A*x = lambda*x, with A symmetric
    iparam(7) = 1

    do

      ! call the reverse communication interface from arpack
    #:if WITH_PARPACK
      call psaupd(comm, ido, "I", nLoc, "SM", nexc, arTol, resid, ncv, vv, nLoc, iparam,&
          & ipntr, workd, workl, lworkl, info)
    #:else
      call saupd(ido, "I", rpa%nxov_rd, "SM", nexc, arTol, resid, ncv, vv, rpa%nxov_rd, iparam,&
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
      call actionAplusB(iGlobal, fGlobal, env, orb, lr, rpa, transChrg, sym, denseDesc, species0,&
          & ovrXev, grndEigVecs, gammaMat, .false., workd(ipntr(1):ipntr(1)+nLoc-1),&
          & workd(ipntr(2):ipntr(2)+nLoc-1))

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
    #:if WITH_PARPACK

      call pseupd (comm, rvec, "All", selection, eval, vv, nLoc, sigma, "I", nLoc,&
          & "SM", nexc, arTol, resid, ncv, vv, nLoc, iparam, ipntr, workd, workl, lworkl, info)

      xpy(:,:) = 0.0_dp
      xpy(iGlobal:fGlobal,:nexc) = vv(:,:nexc)
      call assembleChunks(env, xpy)

     #:else

      call seupd(rvec, "All", selection, eval, xpy, rpa%nxov_rd, sigma, "I", rpa%nxov_rd, "SM",&
          & nexc, arTol, resid, ncv, vv, rpa%nxov_rd, iparam, ipntr, workd, workl, lworkl, info)

    #:endif

      ! check for error on return
      if (info  /=  0) then
        write(tmpStr,"(' Error with ARPACK routine seupd, info = ',I0)")info
        call error(tmpStr)
      end if

    end if

    if (lr%testArnoldi) then
      ! tests for quality of returned eigenpairs
      call openFile(fdArnoldiTest, testArpackOut, mode="w")
      allocate(Hv(rpa%nxov_rd))
      allocate(orthnorm(rpa%nxov_rd,rpa%nxov_rd))
      orthnorm = matmul(transpose(xpy(:,:nExc)),xpy(:,:nExc))

      write(fdArnoldiTest%unit,"(A)")'State Ei deviation    Evec deviation  Norm deviation  Max&
          & non-orthog'
      do iState = 1, nExc

        call actionAplusB(iGlobal, fGlobal, env, orb, lr, rpa, transChrg, sym, denseDesc, species0,&
          & ovrXev, grndEigVecs, gammaMat, .false., xpy(iGlobal:fGlobal,iState),&
          & Hv(iGlobal:fGlobal))

        call assembleChunks(env, Hv)

        write(fdArnoldiTest%unit,"(I4,4E16.8)")iState,&
            & dot_product(Hv,xpy(:,iState))-eval(iState),&
            & sqrt(sum( (Hv-xpy(:,iState)*eval(iState) )**2 )), orthnorm(iState,iState) - 1.0_dp,&
            & max(maxval(orthnorm(:iState-1,iState)), maxval(orthnorm(iState+1:,iState)))
      end do
      call closeFile(fdArnoldiTest)
    end if

    if (rpa%tZVector) then
      xmy(:,:) = 0.0_dp
    end if

    ! Conversion from eigenvectors of the hermitian problem (F) to (X+Y)
    do iState = 1, nExc
      omega = sqrt(eval(iState))
      xpy(:rpa%nxov_rd,iState) = xpy(:rpa%nxov_rd,iState) * sqrt(rpa%wij(:rpa%nxov_rd) / omega)
      if (rpa%tZVector) then
        xmy(:rpa%nxov_rd,iState) = xpy(:rpa%nxov_rd,iState) * omega / rpa%wij(:rpa%nxov_rd)
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
  subroutine buildAndDiagExcMatrixStratmann(iGlobal, fGlobal, env, orb, lr, rpa, transChrg,&
    & denseDesc, ovrXev, grndEigVecs, gammaMat, lrGamma, species0, eval, sym, xpy, xmy)

    !> Starting index of current rank in global RPA vectors
    integer, intent(in) :: iGlobal

    !> End index of current rank in global RPA vectors
    integer, intent(in) :: fGlobal

    !> Environment settings
    type(TEnvironment), intent(inout) :: env

    !> Data type for atomic orbital information
    type(TOrbitals), intent(in) :: orb

    !> Data structure for linear response
    type(TLinResp), intent(in) :: lr

    !> Run time parameters of the Casida routine
    type(TCasidaParameter), intent(in) :: rpa

    !> Machinery for transition charges between single particle levels
    type(TTransCharges), intent(in) :: transChrg

    !> Dense matrix descriptor
    type(TDenseDescr), intent(in) :: denseDesc

    !> overlap times ground state eigenvectors (local in case of MPI)
    real(dp), intent(in) :: ovrXev(:,:,:)

    !> ground state eigenvectors (local in case of MPI)
    real(dp), intent(in) :: grndEigVecs(:,:,:)

    !> Electrostatic matrix
    real(dp), intent(in) :: gammaMat(:,:)

    !> Electrostatic matrix, long-range corrected
    real(dp), allocatable, intent(in) :: lrGamma(:,:)

    !> Central cell chemical species
    integer, intent(in) :: species0(:)

    !> Resulting eigenvalues for transitions (w^2)
    real(dp), intent(out) :: eval(:)

    !> Symmetry to calculate transitions
    character, intent(in) :: sym

    !> Eigenvectors (X+Y)
    real(dp), intent(out) :: xpy(:,:)

    !> Eigenvectors (X-Y), only evaluated if Z-vector is needed
    real(dp), intent(inout), allocatable :: xmy(:,:)

    real(dp), allocatable :: vecB(:,:) ! basis of subspace
    real(dp), allocatable :: evecL(:,:), evecR(:,:) ! left and right eigenvectors of Mnh
    real(dp), allocatable :: vP(:,:), vM(:,:) ! vec. for (A+B)b_i, (A-B)b_i
    ! matrices M_plus, M_minus, M_minus^(1/2), M_minus^(-1/2) and M_herm~=resp. mat on subspace
    real(dp), allocatable :: mP(:,:), mM(:,:), mMsqrt(:,:), mMsqrtInv(:,:), mH(:,:)
    ! Residual vectors
    real(dp), allocatable :: resR(:,:), resL(:,:), dummyM(:,:)
    real(dp), allocatable :: evalInt(:) ! store eigenvectors within routine
    real(dp), allocatable :: vecNorm(:) ! will hold norms of residual vectors
    real(dp) :: dummyReal

    integer :: nExc, nAtom, dummyInt, newVec, iVec, info, iterStrat, nLoc
    integer :: subSpaceDim, prevSubSpaceDim
    integer :: ii, jj, myjj, myii
    character(lc) :: tmpStr

    logical :: didConverge

    ! Local chunk of RPA vectors have this size under MPI
    nLoc = fGlobal - iGlobal + 1

    if (allocated(lr%onSiteMatrixElements)) then
      write(tmpStr,'(A)') 'Onsite corrections not available in Stratmann diagonaliser.'
      call error(tmpStr)
    endif

    ! Number of excited states to solve for
    nExc = size(eval)
    nAtom = size(gammaMat, dim=1)
    @:ASSERT(all(shape(xpy) == [ rpa%nxov_rd, nexc ]))

    ! Small subSpaceDim is faster but leads to convergence problems
    ! if large number of excited states is needed
    if (lr%subSpaceFactorStratmann < 2) then
      write(tmpStr,'(A)') 'SubSpaceFactor for Stratmann solver must be larger than one.'
      call error(tmpStr)
    endif
    subSpaceDim = min(lr%subSpaceFactorStratmann * nExc, rpa%nxov_rd)
    iterStrat = 1

    write(stdOut,'(A)')
    write(stdOut,'(A)') '>> Stratmann diagonalisation of response matrix'
    write(stdOut,'(3x,A,i6,A,i6)') 'Total dimension of A+B: ', rpa%nxov_rd, ' inital subspace: ',&
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

    if (rpa%tZVector) then
      xmy(:,:) = 0.0_dp
    end if

    prevSubSpaceDim = 0
    didConverge = .false.

    ! Solve the linear response problem. Iterative expansion of subspace:
    solveLinResp: do

      if (prevSubSpaceDim > 0) then

        ! Extend subspace matrices:
        do ii = prevSubSpaceDim + 1, subSpaceDim

          call actionAplusB(iGlobal, fGlobal, env, orb, lr, rpa, transChrg, sym, denseDesc,&
              & species0, ovrXev, grndEigVecs, gammaMat, .true., vecB(:,ii), vP(:,ii), lrGamma)
          call actionAminusB(iGlobal, fGlobal, env, orb, lr, rpa, transChrg, denseDesc,&
              & ovrXev, grndEigVecs, vecB(:,ii), vM(:,ii), lrGamma)

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
        call initialSubSpaceMatrixApmB(iGlobal, fGlobal, env, lr, rpa, transChrg, sym, denseDesc,&
            & species0, ovrXev, grndEigVecs, gammaMat, lrGamma, subSpaceDim, vP, vM, mP, mM)

      end if

      call calcMatrixSqrt(mM, subSpaceDim, mMsqrt, mMsqrtInv)

      call symm(dummyM, 'L', mP, mMsqrt, uplo='U')
      call symm(mH, 'L', mMsqrt, dummyM, uplo='U')

      ! Diagonalise in subspace
      call heev(mH, evalInt, 'U', 'V', info)

    #:if WITH_SCALAPACK
      ! required for coherence between processors, as slices of eigenvectors are processed on
      ! separate ranks, so need the same global phase convention:
      call mpifx_bcast(env%mpi%globalComm, mH)
    #:endif

      if (info /= 0) then
        if (lr%subSpaceFactorStratmann * nExc < rpa%nxov_rd) then
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
      didConverge = all(vecNorm < convThreshStrat)

      if ((.not. didConverge) .and. (subSpaceDim > rpa%nxov_rd)) then
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
        if (rpa%tZVector) then
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
        if (vecNorm(ii) > convThreshStrat) then
          newVec = newVec + 1
        endif
      enddo

      call incMemStratmann(subSpaceDim, subSpaceDim + newVec, vecB, vP, vM, mP, mM, mH, mMsqrt,&
            & mMsqrtInv, dummyM, evalInt, evecL, evecR)

      iVec = 0
      do ii = 1, nExc
        if (vecNorm(ii) > convThreshStrat) then
          iVec = iVec + 1
          dummyReal = sqrt(evalInt(ii))
          dummyInt = subSpaceDim + iVec

          do jj = iGlobal, fGlobal
            myjj = jj - iGlobal + 1
            vecB(myjj,dummyInt) = resR(myjj,ii) / (dummyReal - rpa%wij(jj))
          end do

        end if
      end do

      do ii = 1, nExc
        if (vecNorm(nExc+ii) > convThreshStrat) then
          iVec = iVec + 1
          dummyInt = subSpaceDim + iVec

          do jj = iGlobal, fGlobal
            myjj = jj - iGlobal + 1
            vecB(myjj,dummyInt) = resL(myjj,ii) / (dummyReal - rpa%wij(jj))
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
  subroutine getOscillatorStrengths(lr, rpa, sym, eval, xpy, snglPartTransDip, istat, osz,&
      & transitionDipoles)

    !> Data structure for linear response
    type(TLinResp), intent(in) :: lr

    !> Run time parameters of the Casida routine
    type(TCasidaParameter), intent(in) :: rpa

    !> Symmetry of transition
    character, intent(in) :: sym

    !> Low lying eigenvalues of Casida eqn (Omega^2)
    real(dp), intent(in) :: eval(:)

    !> Eigenvectors of Casida eqn (X+Y)
    real(dp), intent(in) :: xpy(:,:)

    !> Dipole moments for single particle transtions
    real(dp), intent(in) :: snglPartTransDip(:,:)

    !> Flag wich if <-1 on entry is returned as the brightest state
    integer, intent(inout) :: istat

    !> Oscilator strengths of transitions
    real(dp), intent(out) :: osz(:)

    !> Resulting transition dipoles
    real(dp), intent(out) :: transitionDipoles(:,:)

    integer :: ii, nmat, oszLoc(1)

    nmat = size(xpy, dim=1)

    transitionDipoles(:,:) = 0.0_dp
    osz(:) = 0.0_dp

    ! Triplet oscillator strength and transition dipole is zero for
    ! closed shell ground state
    if ((.not. lr%tSpin) .and. (sym == "T")) then
      return
    end if

    !$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(ii) SCHEDULE(RUNTIME)
    do ii = 1, size(xpy, dim=2)
      osz(ii) = oscillatorStrength(lr%tSpin, snglPartTransDip, sqrt(eval(ii)), xpy(:,ii),&
          & rpa%sqrOccIA)
    end do
    !$OMP  END PARALLEL DO

    if (istat < 0) then
      ! find largest transition dipole transition
      oszLoc = maxloc(osz)
      istat = oszLoc(1)
    end if

    if (lr%writeTransDip) then
      call transitionDipole(lr%tSpin, snglPartTransDip, xpy, rpa%sqrOccIA, transitionDipoles)
    end if

  end subroutine getOscillatorStrengths


  !> Build right hand side of the equation for the Z-vector and those parts of the W-vectors which
  !! do not depend on Z.
  subroutine getZVectorEqRHS(env, orb, lr, rpa, transChrg, sym, denseDesc, species0, grndEigVal,&
      & ovrXev, grndEigVecs, gammaMat, lrGamma, omega, xpy, xmy, rhs, t, wov, woo, wvv)

    !> Environment settings
    type(TEnvironment), intent(inout) :: env

    !> Data type for atomic orbital information
    type(TOrbitals), intent(in) :: orb

    !> Data structure for linear response
    type(TLinResp), intent(in) :: lr

    !> Run time parameters of the Casida routine
    type(TCasidaParameter), intent(in) :: rpa

    !> Machinery for transition charges between single particle levels
    type(TTransCharges), intent(in) :: transChrg

    !> Symmetry of the transitions
    character, intent(in) :: sym

    !> Dense matrix descriptor
    type(TDenseDescr), intent(in) :: denseDesc

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

    !> Excitation energies
    real(dp), intent(in) :: omega

    !> X+Y Furche term
    real(dp), intent(in) :: xpy(:)

    !> X-Y Furche term
    real(dp), intent(in) :: xmy(:)

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

    integer :: i, j, a, b, ias, ibs, abs, ij, ab, jas, ijs, s, nSpin, soo(2), svv(2), nOrb, nxov
    integer :: iGlobal, fGlobal
    real(dp), allocatable :: xpyq(:), qTr(:), gamxpyq(:), qgamxpyq(:,:), gamqt(:)
    real(dp), allocatable :: xpyqds(:), gamxpyqds(:)
    real(dp), allocatable :: vecHvvXpY(:), vecHvvXmY(:), vecHooXpY(:), vecHooXmY(:)
    real(dp), allocatable :: vecHovT(:), vecHooT(:)
    real(dp) :: tmp1, tmp2, fact
    logical :: tSpin

    nxov = size(rhs)
    nOrb = orb%nOrb

    call distributeRangeInChunks(env, 1, nxov, iGlobal, fGlobal)

    allocate(xpyq(lr%nAtom))
    allocate(qTr(lr%nAtom))
    allocate(gamxpyq(lr%nAtom))
    allocate(gamqt(lr%nAtom))

    t(:,:,:) = 0.0_dp
    rhs(:) = 0.0_dp
    wov(:) = 0.0_dp
    woo(:,:) = 0.0_dp
    wvv(:,:) = 0.0_dp

    nSpin = size(t, dim=3)

    ! Transition charges use compound index ijs = ij + soo(s)
    soo(:) = [0, rpa%nxoo_ud(1)]
    svv(:) = [0, rpa%nxvv_ud(1)]

    allocate(qgamxpyq(max(maxval(rpa%nxoo_ud), maxval(rpa%nxvv_ud)), size(rpa%nocc_ud)))
    qgamxpyq(:,:) = 0.0_dp

    if (nSpin == 2) then
      tSpin = .true.
      allocate(xpyqds(lr%nAtom))
      allocate(gamxpyqds(lr%nAtom))
    else
      tSpin = .false.
    end if

    ! Build t_ab = 0.5 * sum_i (X+Y)_ia (X+Y)_ib + (X-Y)_ia (X-Y)_ib
    ! and w_ab = Q_ab with Q_ab as in (B16) but with corrected sign.
    ! factor 1 / (1 + delta_ab) follows later
    do ias = iGlobal, fGlobal
      call indxov(rpa%win, ias, rpa%getIA, i, a, s)

      ! BA: is T_aa = 0?
      do b = rpa%nocc_ud(s) + 1, a
        ibs = rpa%iaTrans(i, b, s)
        ab = rpa%iaTrans(a, b, s) - svv(s)
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
      do j = i, rpa%nocc_ud(s)
        jas = rpa%iaTrans(j,a,s)

        ! ADG: assume no constraint on occ space atm (nocc_r = nocc)
        ! otherwise first argument should be nocc - nocc_r
        ij = rpa%iaTrans(i, j, s) - soo(s)
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

    call assembleChunks(env, t)
    call assembleChunks(env, wvv)

    ! xpyq = Q * xpy
    xpyq(:) = 0.0_dp
    call transChrg%qMatVec(env, denseDesc, ovrXev, grndEigVecs, rpa%getIA, rpa%win, xpy, xpyq)

    if (.not. tSpin) then  ! ---- spin-unpolarized case ----
      call distributeRangeInChunks(env, 1, rpa%nxvv_ud(1), iGlobal, fGlobal)

      ! qgamxpyq(ab) = sum_jc K_ab,jc (X+Y)_jc
      if (sym == "S") then
        call hemv(gamxpyq, gammaMat,  xpyq)
        do ab = iGlobal, fGlobal
          qTr(:) = transChrg%qTransAB(ab, env, denseDesc, ovrXev, grndEigVecs, rpa%getAB)
          qgamxpyq(ab, 1) = 2.0_dp * sum(qTr * gamxpyq)
        end do
      else ! triplet case
        do ab = iGlobal, fGlobal
          qTr(:) = transChrg%qTransAB(ab, env, denseDesc, ovrXev, grndEigVecs, rpa%getAB)
          qgamxpyq(ab, 1) = 2.0_dp * sum(qTr * xpyq * lr%spinW(species0))
        end do
      end if

      call assembleChunks(env, qgamxpyq)

    else  ! ---- spin-polarized case -----

      xpyqds(:) = 0.0_dp
      call transChrg%qMatVecDs(env, denseDesc, ovrXev, grndEigVecs, rpa%getIA, rpa%win, xpy, xpyqds)

      call hemv(gamxpyq, gammaMat,  xpyq)
      do s = 1, 2
        if (s == 1) then
          fact = 1.0_dp
        else
          fact = -1.0_dp
        end if

        call distributeRangeInChunks(env, 1, rpa%nxvv_ud(s), iGlobal, fGlobal)

        do ab = iGlobal, fGlobal
          qTr(:) = transChrg%qTransAB(ab + svv(s), env, denseDesc, ovrXev, grndEigVecs, rpa%getAB)
          qgamxpyq(ab, s) = sum(qTr * gamxpyq)
          !magnetization part
          qgamxpyq(ab, s) = qgamxpyq(ab, s) + fact * sum(qTr * xpyqds * lr%spinW(species0))
        end do

        call assembleChunks(env, qgamxpyq(:,s))
      end do

    end if

    call distributeRangeInChunks(env, 1, nxov, iGlobal, fGlobal)

    ! rhs(ia) -= Qia = sum_b (X+Y)_ib * qgamxpyq(ab))
    do ias = iGlobal, fGlobal
      call indxov(rpa%win, ias, rpa%getIA, i, a, s)

      do b = rpa%nocc_ud(s) + 1, a
        ab = rpa%iaTrans(a, b, s) - svv(s)
        ibs = rpa%iaTrans(i, b, s)
        rhs(ias) = rhs(ias) - 2.0_dp * xpy(ibs) * qgamxpyq(ab, s)
        ! Since qgamxpyq has only upper triangle
        if (a /= b) then
          rhs(ibs) = rhs(ibs) - 2.0_dp * xpy(ias) * qgamxpyq(ab, s)
        end if
      end do
    end do

    ! -rhs = -rhs - sum_j (X + Y)_ja H + _ij[X + Y]
    qgamxpyq(:,:) = 0.0_dp
    if (.not. tSpin) then  ! ---- spin-unpolarized case ----
      call distributeRangeInChunks(env, 1, rpa%nxoo_ud(1), iGlobal, fGlobal)

      if (sym == "S") then
        do ij = iGlobal, fGlobal
          qgamxpyq(ij, 1) = 0.0_dp
          qTr(:) = transChrg%qTransIJ(ij, env, denseDesc, ovrXev, grndEigVecs, rpa%getIJ)
          ! qgamxpyq(ij) = sum_kb K_ij,kb (X+Y)_kb
          qgamxpyq(ij, 1) = 2.0_dp * sum(qTr * gamxpyq)
        end do
      else
        do ij = iGlobal, fGlobal
          qgamxpyq(ij, 1) = 0.0_dp
          qTr(:) = transChrg%qTransIJ(ij, env, denseDesc, ovrXev, grndEigVecs, rpa%getIJ)
          qgamxpyq(ij, 1) = 2.0_dp * sum(qTr * xpyq * lr%spinW(species0))
        end do
      end if

      call assembleChunks(env, qgamxpyq)

    else  ! ---- spin-polarized case -----

      do s = 1, 2
        if (s == 1) then
          fact = 1.0_dp
        else
          fact = -1.0_dp
        end if
        call distributeRangeInChunks(env, 1, rpa%nxoo_ud(s), iGlobal, fGlobal)
        do ij = iGlobal, fGlobal
          qgamxpyq(ij, s) = 0.0_dp
          qTr(:) = transChrg%qTransIJ(ij + soo(s), env, denseDesc, ovrXev, grndEigVecs, rpa%getIJ)
          qgamxpyq(ij, s) = sum(qTr * gamxpyq)
          !magnetization part
          qgamxpyq(ij, s) = qgamxpyq(ij, s) + fact * sum(qTr * xpyqds * lr%spinW(species0))
        end do
        call assembleChunks(env, qgamxpyq(:,s))
      end do

    end if

    call distributeRangeInChunks(env, 1, nxov, iGlobal, fGlobal)

    ! rhs(ia) += Qai = sum_j (X+Y)_ja qgamxpyq(ij)
    ! add Qai to Wia as well.
    do ias = iGlobal, fGlobal
      call indxov(rpa%win, ias, rpa%getIA, i, a, s)
      do j = i, rpa%nocc_ud(s)
        jas = rpa%iaTrans(j, a, s)
        ij = rpa%iaTrans(i, j, s) - soo(s)
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

    call assembleChunks(env, rhs)
    call assembleChunks(env, wov)

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

      call distributeRangeInChunks(env, 1, rpa%nxoo_ud(s), iGlobal, fGlobal)

      do ij = iGlobal, fGlobal
        i = rpa%getIJ(ij + soo(s), 1)
        j = rpa%getIJ(ij + soo(s), 2)
        qTr(:) = transChrg%qTransIJ(ij + soo(s), env, denseDesc, ovrXev, grndEigVecs, rpa%getIJ)
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

      call distributeRangeInChunks(env, 1, rpa%nxvv_ud(s), iGlobal, fGlobal)

      ! gamxpyq(iAt2) += sum_ab q_ab(iAt2) T_ab
      do ab = iGlobal, fGlobal
        a = rpa%getAB(ab + svv(s), 1)
        b = rpa%getAB(ab + svv(s), 2)
        qTr(:) = transChrg%qTransAB(ab + svv(s), env, denseDesc, ovrXev, grndEigVecs, rpa%getAB)
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

    call assembleChunks(env, gamxpyq)
    if (tSpin) then
      call assembleChunks(env, gamxpyqds)
    end if


    ! gamqt(iAt1) = sum_iAt2 gamma_iAt1,iAt2 gamxpyq(iAt2)
    call hemv(gamqt, gammaMat, gamxpyq)

    ! rhs -= sum_q^ia(iAt1) gamxpyq(iAt1)
    if (.not. tSpin) then
      call transChrg%qVecMat(env, denseDesc, ovrXev, grndEigVecs, rpa%getIA, rpa%win,&
          & -4.0_dp * gamqt, rhs)
    else
      call transChrg%qVecMat(env, denseDesc, ovrXev, grndEigVecs, rpa%getIA, rpa%win,&
          & -2.0_dp * gamqt, rhs)
      call transChrg%qVecMatDs(env, denseDesc, ovrXev, grndEigVecs, rpa%getIA, rpa%win,&
           & -2.0_dp * gamxpyqds * lr%spinW(species0), rhs)
    end if

    ! Furche vectors
    do s = 1, nSpin
      if (s == 1) then
        fact = 1.0_dp
      else
        fact = -1.0_dp
      end if

      call distributeRangeInChunks(env, 1, rpa%nxoo_ud(s), iGlobal, fGlobal)

      do ij = iGlobal, fGlobal
        qTr(:) = transChrg%qTransIJ(ij + soo(s), env, denseDesc, ovrXev, grndEigVecs, rpa%getIJ)
        if (.not. tSpin) then
          woo(ij,s) = woo(ij,s) + 4.0_dp * sum(qTr * gamqt)
        else
          woo(ij,s) = woo(ij,s) + 2.0_dp * sum(qTr * gamqt)
          woo(ij,s) = woo(ij,s) + 2.0_dp * fact * sum(qTr * gamxpyqds * lr%spinW(species0))
        end if
      end do
    end do

    ! Contributions due to range-separation
    if (rpa%tHybridXc) then

      allocate(vecHvvXpY(sum(rpa%nxvv_ud)))
      allocate(vecHvvXmY(sum(rpa%nxvv_ud)))
      allocate(vecHooXpY(sum(rpa%nxoo_ud)))
      allocate(vecHooXmY(sum(rpa%nxoo_ud)))
      allocate(vecHovT(nxov))
      allocate(vecHooT(sum(rpa%nxoo_ud)))

      call getHvvXY(env, orb, lr, rpa, transChrg, denseDesc, ovrXev, grndEigVecs, lrGamma,  1, xpy,&
          & vecHvvXpY)

      call getHvvXY(env, orb, lr, rpa, transChrg, denseDesc, ovrXev, grndEigVecs, lrGamma, -1, xmy,&
          & vecHvvXmY)

      call getHooXY(env, orb, lr, rpa, transChrg, denseDesc, ovrXev, grndEigVecs, lrGamma,  1, xpy,&
          & vecHooXpY)

      call getHooXY(env, orb, lr, rpa, transChrg, denseDesc, ovrXev, grndEigVecs, lrGamma, -1, xmy,&
          & vecHooXmY)

      call getHovT(env, orb, lr, rpa, transChrg, denseDesc, ovrXev, grndEigVecs, lrGamma, t, vecHovT)

      !TODO: i can not parallelized this because i had to do the assemble of rhs before
      ! in order to get gamqt (see before), if i dont do the assemble then gamqt is wrong.
      ! A possible solution is to move all the block of gamqt and the rhs before that to the end
      ! after the hybridxc section, then do the assemble of rsh and then calculate gamqt and the other eq of rhs
      ! call distributeRangeInChunks(env, 1, nxov, iGlobal, fGlobal)

      do ias = 1, nxov

        call indXov(rpa%win, ias, rpa%getIA, i, a, s)
        do b = rpa%nocc_ud(s) + 1, nOrb
          ibs = rpa%iaTrans(i, b, s)
          abs = rpa%iaTrans(a, b, s)
          rhs(ias) = rhs(ias) - cExchange * xpy(ibs) * vecHvvXpY(abs)
          if (a >= b) then
            rhs(ias) = rhs(ias) - cExchange * xmy(ibs) * vecHvvXmY(abs)
          else
            rhs(ias) = rhs(ias) + cExchange * xmy(ibs) * vecHvvXmY(abs)
          end if
        end do

        do j = 1, rpa%nocc_ud(s)
          jas = rpa%iaTrans(j, a, s)
          ijs = rpa%iaTrans(i, j, s)
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

      call getHooT(env, orb, lr, rpa, transChrg, denseDesc, ovrXev, grndEigVecs, lrGamma, t, vecHooT)

      call distributeRangeInChunks(env, 1, rpa%nxoo_ud(s), iGlobal, fGlobal)

      ! Array woo should be made 1D
      do s = 1, nSpin
        do ij = iGlobal, fGlobal
          i = rpa%getIJ(ij + soo(s), 1)
          j = rpa%getIJ(ij + soo(s), 2)
          ijs = rpa%iaTrans(i, j, s)
          woo(ij,s) = woo(ij,s) + cExchange * vecHooT(ijs)
        end do
      end do

      call assembleChunks(env, woo)

    endif

  end subroutine getZVectorEqRHS


  !> Solving the (A+B) Z = -R equation via diagonally preconditioned conjugate gradient.
  subroutine solveZVectorPrecond(env, orb, lr, rpa, transChrg, denseDesc, species0, ovrXev,&
       & grndEigVecs, gammaMat, lrGamma, rhs)

    !> Environment settings
    type(TEnvironment), intent(inout) :: env

    !> Data type for atomic orbital information
    type(TOrbitals), intent(in) :: orb

    !> Data structure for linear response
    type(TLinResp), intent(in) :: lr

    !> Run time parameters of the Casida routine
    type(TCasidaParameter), intent(in) :: rpa

    !> machinery for transition charges between single particle levels
    type(TTransCharges), intent(in) :: transChrg

    !> Dense matrix descriptor
    type(TDenseDescr), intent(in) :: denseDesc

    !> Chemical species of the atoms
    integer, intent(in) :: species0(:)

    !> Overlap times eigenvector (nOrb, nOrb)
    real(dp), intent(in) :: ovrXev(:,:,:)

    !> Eigenvectors (nOrb, nOrb)
    real(dp), intent(in) :: grndEigVecs(:,:,:)

    !> DFTB gamma matrix (nAtm, nAtom)
    real(dp), intent(in) :: gammaMat(:,:)

    !> Long-range Gamma
    real(dp), allocatable, intent(in) :: lrGamma(:,:)

    !> On entry -R, on exit Z
    real(dp), intent(inout) :: rhs(:)

    integer :: nxov, iGlobal, fGlobal
    integer :: ia, kk, i, a, s, iis, aas
    real(dp), allocatable :: qTr(:), qTmp(:), P(:)
    real(dp) :: rhs2(size(rhs)), rkm1(size(rhs)), zkm1(size(rhs)), pkm1(size(rhs)), apk(size(rhs))
    real(dp) :: rs, alphakm1, tmp1, tmp2, bkm1

    call distributeRangeInChunks(env, 1, rpa%nxov_rd, iGlobal, fGlobal)

    nxov = rpa%nxov_rd
    allocate(qTr(lr%nAtom))
    allocate(qTmp(lr%nAtom))

    ! diagonal preconditioner
    ! P^-1 = 1 / (A+B)_ia,ia (diagonal of the supermatrix sum A+B)
    allocate(P(nxov)); P(:) = 0.0_dp
    do ia = iGlobal, fGlobal
      qTr(:) = transChrg%qTransIA(ia, env, denseDesc, ovrXev, grndEigVecs, rpa%getIA, rpa%win)
      call hemv(qTmp, gammaMat, qTr)
      if (.not. lr%tSpin) then
        rs = 4.0_dp * dot_product(qTr, qTmp) + rpa%wij(ia)
      else
        rs = 2.0_dp * dot_product(qTr, qTmp) + rpa%wij(ia)
        rs = rs + 2.0_dp * sum(qTr * qTr * lr%spinW(species0))
      end if

      ! Possibly reorder spin case
      if (rpa%tHybridXc) then
        call hemv(qTmp, lrGamma, qTr)
        rs = rs - cExchange * dot_product(qTr, qTmp)
        call indXov(rpa%win, ia, rpa%getIA, i, a, s)
        iis = rpa%iaTrans(i, i, s)
        qTr(:) = transChrg%qTransIJ(iis, env, denseDesc, ovrXev, grndEigVecs, rpa%getIJ)
        call hemv(qTmp, lrGamma, qTr)
        aas = rpa%iaTrans(a, a, s)
        qTr(:) = transChrg%qTransAB(aas, env, denseDesc, ovrXev, grndEigVecs, rpa%getAB)
        rs = rs - cExchange * dot_product(qTr, qTmp)
      end if

      P(ia) = 1.0_dp / rs
    end do

    call assembleChunks(env, P)

    ! Free some space, before entering the actionAplusB routine
    deallocate(qTr)

    ! unit vector as initial guess solution
    rhs2(:) = 1.0_dp / sqrt(real(nxov,dp))
    rkm1(:) = 0.0_dp

    ! action of matrix on vector
    ! we need the singlet action even for triplet excitations!
    call actionAplusB(iGlobal, fGlobal, env, orb, lr, rpa, transChrg, 'S', denseDesc, species0,&
        & ovrXev, grndEigVecs, gammaMat, .true., rhs2(iGlobal:fGlobal), rkm1(iGlobal:fGlobal),&
        & lrGamma)

    call assembleChunks(env, rkm1)

    rkm1(:) = rhs - rkm1
    zkm1(:) = P * rkm1
    pkm1(:) = zkm1

    ! Iteration: should be convergent in at most nxov steps for a quadradic surface, so set higher
    do kk = 1, nxov**2

      ! action of matrix on vector
      apk = 0.0_dp
      call actionAplusB(iGlobal, fGlobal, env, orb, lr, rpa, transChrg, 'S', denseDesc, species0,&
          & ovrXev, grndEigVecs, gammaMat, .true., pkm1(iGlobal:fGlobal), apk(iGlobal:fGlobal), lrGamma)

      call assembleChunks(env, apk)

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
  !! 2.
  subroutine calcWvectorZ(env, orb, lr, rpa, transChrg, denseDesc, species0, ovrXev, grndEigVecs,&
      & grndEigVal, gammaMat, lrGamma, zz, wov, woo, wvv)

    !> Environment settings
    type(TEnvironment), intent(inout) :: env

    !> Data type for atomic orbital information
    type(TOrbitals), intent(in) :: orb

    !> Data structure for linear response
    type(TLinResp), intent(in) :: lr

    !> Run time parameters of the Casida routine
    type(TCasidaParameter), intent(in) :: rpa

    !> machinery for transition charges between single particle levels
    type(TTransCharges), intent(in) :: transChrg

    !> Dense matrix descriptor
    type(TDenseDescr), intent(in) :: denseDesc

    !> Central cell chemical species
    integer, intent(in) :: species0(:)

    !> Overlap times ground state wavefunctions
    real(dp), intent(in) :: ovrXev(:,:,:)

    !> Ground state wavefunctions
    real(dp), intent(in) :: grndEigVecs(:,:,:)

    !> Ground state MO-energies
    real(dp), intent(in) :: grndEigVal(:,:)

    !> Softened coulomb matrix
    real(dp), intent(in) :: gammaMat(:,:)

    !> Long-range Gamma
    real(dp), allocatable, intent(in) :: lrGamma(:,:)

    !> Z vector
    real(dp), intent(in) :: zz(:)

    !> W vector occupied-virtual part
    real(dp), intent(inout) :: wov(:)

    !> W vector occupied part
    real(dp), intent(inout) :: woo(:,:)

    !> W vector virtual part
    real(dp), intent(inout) :: wvv(:,:)

    integer :: nSpin, soo(2), svv(2)
    integer :: ij, ias, ijs, ab, i, j, a, b, s, iGlobal, fGlobal
    real(dp), allocatable :: qTr(:), gamxpyq(:), zq(:), zqds(:), vecHooZ(:)
    real(dp), allocatable :: wovLoc(:), wooLoc(:,:), wvvLoc(:,:)
    real(dp) :: fact

    nSpin = size(grndEigVal, dim=2)

    allocate(qTr(lr%nAtom))
    allocate(gamxpyq(lr%nAtom))
    allocate(zq(lr%nAtom))

    soo(:) = [0, rpa%nxoo_ud(1)]
    svv(:) = [0, rpa%nxvv_ud(1)]

    if (nSpin == 2) then
      allocate(zqds(lr%nAtom))
    end if

    ! MPI local arrays
    call distributeRangeInChunks(env, 1, rpa%nxov_rd, iGlobal, fGlobal)
    allocate(wovLoc(rpa%nxov_rd)); wovLoc = 0.0_dp

    ! Adding missing epsilon_i * Z_ia term to W_ia
    do ias = iGlobal, fGlobal
      call indxov(rpa%win, ias, rpa%getIA, i, a, s)
      wovLoc(ias) = zz(ias) * grndEigVal(i, s)
    end do
    call assembleChunks(env, wovLoc)
    wov = wov + wovLoc
    deallocate(wovLoc)

    ! Missing sum_kb 4 K_ijkb Z_kb term in W_ij: zq(iAt1) = sum_kb q^kb(iAt1) Z_kb
    zq(:) = 0.0_dp
    call transChrg%qMatVec(env, denseDesc, ovrXev, grndEigVecs, rpa%getIA, rpa%win, zz, zq)
    call hemv(gamxpyq, gammaMat, zq)

    if (lr%tSpin) then
      zqds(:) = 0.0_dp
      call transChrg%qMatVecDs(env, denseDesc, ovrXev, grndEigVecs, rpa%getIA, rpa%win, zz, zqds)
    end if

    ! MPI local arrays
    allocate(wooLoc(size(woo, dim=1),nSpin)); wooLoc = 0.0_dp

    ! sum_iAt1 qTr(iAt1) gamxpyq(iAt1)
    do s = 1, nSpin
      if (s == 1) then
        fact = 1.0_dp
      else
        fact = -1.0_dp
      end if
      call distributeRangeInChunks(env, 1, rpa%nxoo_ud(s), iGlobal, fGlobal)

      do ij = iGlobal, fGlobal
        qTr(:) = transChrg%qTransIJ(ij + soo(s), env, denseDesc, ovrXev, grndEigVecs, rpa%getIJ)
        ! W contains 1/2 for i == j.
        if (.not. lr%tSpin) then
          wooLoc(ij,s) = woo(ij,s) + 4.0_dp * sum(qTr * gamxpyq)
        else
          wooLoc(ij,s) = woo(ij,s) + 2.0_dp * sum(qTr * gamxpyq)
          wooLoc(ij,s) = wooLoc(ij,s) + 2.0_dp * fact * sum(qTr * zqds * lr%spinW(species0))
        end if
      end do
    end do

    if (rpa%tHybridXc) then

      allocate(vecHooZ(sum(rpa%nxoo_ud)))
      call getHooXY(env, orb, lr, rpa, transChrg, denseDesc, ovrXev, grndEigVecs, lrGamma, 1, zz,&
          & vecHooZ)

      ! Array woo should be made 1D
      do s = 1, nSpin
        call distributeRangeInChunks(env, 1, rpa%nxoo_ud(s), iGlobal, fGlobal)
        do ij = iGlobal, fGlobal
          i = rpa%getIJ(ij + soo(s), 1)
          j = rpa%getIJ(ij + soo(s), 2)
          ijs = rpa%iaTrans(i, j, s)
          wooLoc(ij,s) = wooLoc(ij,s) + cExchange * vecHooZ(ijs)
        end do
      end do

    end if

    ! Divide diagonal elements of W_ij by 2.
    do s = 1, nSpin
      call distributeRangeInChunks(env, 1, rpa%nxoo_ud(s), iGlobal, fGlobal)
      do ij = iGlobal, fGlobal
        i = rpa%getIJ(ij + soo(s), 1)
        j = rpa%getIJ(ij + soo(s), 2)
        if (i == j) then
          wooLoc(ij,s) = 0.5_dp * wooLoc(ij,s)
        end if
      end do
    end do

    call assembleChunks(env, wooLoc)
    woo = wooLoc
    deallocate(wooLoc)

    ! MPI local arrays
    allocate(wvvLoc(size(wvv, dim=1),nSpin)); wvvLoc = 0.0_dp

    ! Divide diagonal elements of W_ab by 2.
    do s = 1, nSpin
      call distributeRangeInChunks(env, 1, rpa%nxvv_ud(s), iGlobal, fGlobal)

      do ab = iGlobal, fGlobal
        a = rpa%getAB(ab + svv(s), 1)
        b = rpa%getAB(ab + svv(s), 2)
        wvvLoc(ab,s) = wvv(ab,s)
        if (a == b) then
          wvvLoc(ab,s) = 0.5_dp * wvv(ab,s)
        end if
      end do
    end do

    call assembleChunks(env, wvvLoc)
    wvv = wvvLoc
    deallocate(wvvLoc)

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
  !! Note: assumes both triangles of both square matrices are filled
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
  subroutine addGradients(env, orb, lr, rpa, transChrg, hybridXc, denseDesc, sym, species0, ovrXev,&
      & grndEigVecs, gammaMat, lrGamma, coord0, dq_ud, dqex, shift, xpy, xmy, woo, wov, wvv,&
      & skHamCont, skOverCont, derivator, rhoSqr, pc, excgrad, deltaRho)

    !> Environment settings
    type(TEnvironment), intent(inout) :: env

    !> Data type for atomic orbital information
    type(TOrbitals), intent(in) :: orb

    !> Data structure for linear response
    type(TLinResp), intent(in) :: lr

    !> Run time parameters of the Casida routine
    type(TCasidaParameter), intent(in) :: rpa

    !> machinery for transition charges between single particle levels
    type(TTransCharges), intent(in) :: transChrg

    !> Data for range-separated calculation
    class(THybridXcFunc), allocatable, intent(inout) :: hybridXc

    !> Dense matrix descriptor
    type(TDenseDescr), intent(in) :: denseDesc

    !> Symmetry of the transition
    character, intent(in) :: sym

    !> Central cell chemical species
    integer, intent(in) :: species0(:)

    !> Overlap times ground state eigenvectors
    real(dp), intent(in) :: ovrXev(:,:,:)

    !> Ground state eigenvectors
    real(dp), intent(in) :: grndEigVecs(:,:,:)

    !> Softened coulomb matrix
    real(dp), intent(in) :: gammaMat(:,:)

    !> Electrostatic matrix, long-range corrected
    real(dp), allocatable, intent(in) :: lrGamma(:,:)

    !> Central cell atomic coordinates
    real(dp), intent(in) :: coord0(:,:)

    !> Ground state gross charges
    real(dp), intent(in) :: dq_ud(:,:)

    !> Charge differences from ground to excited state
    real(dp), intent(in) :: dqex(:,:)

    !> Ground state potentials (shift vector)
    real(dp), intent(in) :: shift(:)

    !> X+Y Furche term
    real(dp), intent(in) :: xpy(:)

    !> X-Y Furche term
    real(dp), intent(in) :: xmy(:)

    !> W vector occupied part
    real(dp), intent(in) :: woo(:,:)

    !> W vector occupied-virtual part
    real(dp), intent(in) :: wov(:)

    !> W vector virtual part
    real(dp), intent(in) :: wvv(:,:)

    !> H0 data
    type(TSlakoCont), intent(in) :: skHamCont

    !> Overlap data
    type(TSlakoCont), intent(in) :: skOverCont

    !> Differentiator for the non-scc matrices
    class(TNonSccDiff), intent(in) :: derivator

    !> Ground state density matrix
    real(dp), intent(in) :: rhoSqr(:,:,:)

    !> Transition density matrix
    real(dp), intent(in) :: pc(:,:,:)

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
    real(dp), allocatable :: deltaRhoGlobal(:,:,:), grndEigVecsGlobal(:,:,:), DensGlobal(:,:,:)
    real(dp) :: tmp1, tmp2, tmp3, tmp4, tmp6, tmp8, tmp9, tmp10, rab
    real(dp) :: diffvec(3), dgab(3), tmpVec(3), tmp3a, tmp3b, tmprs, tmprs2, tmps(2)
    integer, allocatable :: species(:)
    integer :: ia, i, j, a, b, ab, ij, m, n, mu, nu, xyz, iAt1, iAt2, ka
    integer :: indalpha, indalpha1, indbeta, indbeta1, soo(2), svv(2)
    integer :: iSp1, iSp2, iSpin, nSpin, nOrb, iGlobal, fGlobal

    nSpin = size(grndEigVecs, dim=3)
    nOrb = orb%nOrb

    allocate(shift_excited(lr%nAtom, nSpin))
    allocate(xpyq(lr%nAtom))
    allocate(shxpyq(lr%nAtom, nSpin))
    allocate(xpycc(nOrb, nOrb, nSpin))
    allocate(wcc(nOrb, nOrb, nSpin))
    allocate(qTr(lr%nAtom))
    allocate(temp(nOrb))
    allocate(tmp5(nSpin))
    allocate(tmp7(nSpin))

    allocate(Dens(nOrb, nOrb), DensGlobal(nOrb,nOrb,size(rhoSqr,dim=3)))
    do iSpin = 1, size(rhoSqr, dim=3)
      call distrib2replicated(env%blacs%orbitalGrid, denseDesc%blacsOrbSqr, &
                           &  rhoSqr(:,:,iSpin), DensGlobal(:,:,iSpin))
    enddo
    Dens(:,:) = sum(DensGlobal, dim=3)
    deallocate(DensGlobal)

    ! NOTE: probably this is not necessary, we need to check
    ! Symmetrize RhoSqr
    do mu = 1, size(Dens, dim=1)
      do nu = mu + 1, size(Dens, dim=2)
        Dens(mu,nu) = Dens(nu,mu)
      end do
    end do

    allocate(dH0(orb%mOrb, orb%mOrb, 3))
    allocate(dSo(orb%mOrb, orb%mOrb, 3))

    soo(:) = [0, rpa%nxoo_ud(1)]
    svv(:) = [0, rpa%nxvv_ud(1)]

    allocate(dq(lr%nAtom))
    dq(:) = dq_ud(:,1)

    if (lr%tSpin) then
      allocate(dm(lr%nAtom))
      allocate(xpyqds(lr%nAtom))
      allocate(tmp11(nSpin))

      !FIXME: here nOrb is the global value but rhoSqr has dimension of nOrb local
      !TODO: The test NH forces does not run even for single process.
      allocate(SpinDens(nOrb,nOrb))
      SpinDens(:,:) = rhoSqr(:,:,1) - rhoSqr(:,:,2)

      allocate(dsigma(2))
      dsigma(1) = 1.0_dp
      dsigma(2) = -1.0_dp
      dm(:) = dq_ud(:,2)
    end if

    if (rpa%tHybridXc) then
      allocate(xmycc(nOrb, nOrb, nSpin))
      allocate(xpyas(nOrb, nOrb, nSpin))
      allocate(xmyas(nOrb, nOrb, nSpin))
      allocate(PS(nOrb, nOrb, nSpin))
      allocate(DS(nOrb, nOrb, nSpin))
      allocate(SPS(nOrb, nOrb, nSpin))
      allocate(SDS(nOrb, nOrb, nSpin))
      allocate(SX(nOrb, nOrb, nSpin))
      allocate(XS(nOrb, nOrb, nSpin))
      allocate(SXS(nOrb, nOrb, nSpin))
      allocate(SY(nOrb, nOrb, nSpin))
      allocate(YS(nOrb, nOrb, nSpin))
      allocate(SYS(nOrb, nOrb, nSpin))
      allocate(overlap(nOrb, nOrb))
      allocate(lrGammaOrb(nOrb, nOrb))
      allocate(gammaLongRangePrime(3, lr%nAtom, lr%nAtom))

      ! Convert local arrays to global
      allocate(deltaRhoGlobal(norb,norb,size(deltaRho,dim=3)))
      allocate(grndEigVecsGlobal(norb,norb,size(grndEigVecs,dim=3)))
      do iSpin = 1, nSpin
        call distrib2replicated(env%blacs%orbitalGrid, denseDesc%blacsOrbSqr, &
                             &  deltaRho(:,:,iSpin), deltaRhoGlobal(:,:,iSpin))
        call distrib2replicated(env%blacs%orbitalGrid, denseDesc%blacsOrbSqr, &
                             &  grndEigVecs(:,:,iSpin), grndEigVecsGlobal(:,:,iSpin))
      enddo

      ! Symmetrize deltaRhoGlobal
      do mu = 1, size(deltaRhoGlobal, dim=1)
        do nu = mu + 1, size(deltaRhoGlobal, dim=2)
          deltaRhoGlobal(mu,nu,:) = deltaRhoGlobal(nu,mu,:)
        end do
      end do

      ! Compute long-range gamma derivative
      call distributeRangeInChunks(env, 1, lr%nAtom, iGlobal, fGlobal)
      gammaLongRangePrime(:,:,:) = 0.0_dp
      do iAt1 = iGlobal, fGlobal
        do iAt2 = 1, lr%nAtom
          if (iAt1 /= iAt2) then
            call getDirectionalCamGammaPrimeValue(hybridXc, tmpVec, iAt1, iAt2)
            gammaLongRangePrime(:, iAt1, iAt2) = tmpVec
          end if
        end do
      end do
      call assembleChunks(env,gammaLongRangePrime)

      ! Symmetrize S (can't we get S from caller?)
      call getSqrS(coord0, lr%nAtom, skOverCont, orb, denseDesc%iAtomStart, species0, overlap)
      call getSqrGamma(lr%nAtom, lrGamma, denseDesc%iAtomStart, lrGammaOrb)

    end if

    excgrad(:,:) = 0.0_dp

    ! excited state potentials at atomic sites
    do iSpin = 1, nSpin
      call hemv(shift_excited(:,iSpin), gammaMat, dqex(:,iSpin))
    end do

    ! xypq(alpha) = sum_ia (X+Y)_ia q^ia(alpha)
    ! complexity nOrb * nOrb * nOrb
    xpyq(:) = 0.0_dp
    call transChrg%qMatVec(env, denseDesc, ovrXev, grndEigVecs, rpa%getIA, rpa%win, xpy, xpyq)

    ! complexity nOrb * nOrb
    shxpyq(:,:) = 0.0_dp
    if (.not. lr%tSpin) then
      if (sym == "S") then
        call hemv(shxpyq(:,1), gammaMat, xpyq)
      else
        shxpyq(:,1) = xpyq(:) * lr%spinW(species0)
      end if
    else
      xpyqds(:) = 0.0_dp
      call transChrg%qMatVecDs(env, denseDesc, ovrXev, grndEigVecs, rpa%getIA, rpa%win, xpy, xpyqds)
      do iSpin = 1, nSpin
        call hemv(shxpyq(:,iSpin), gammaMat, xpyq)
        shxpyq(:,iSpin) = shxpyq(:,iSpin) + dsigma(iSpin) * lr%spinW(species0) * xpyqds
        shxpyq(:,iSpin) = 0.5_dp * shxpyq(:,iSpin)
      end do
    end if

    ! calculate xpycc
    ! (xpycc)_{mu nu} = sum_{ia} (X + Y)_{ia} (grndEigVecs(mu,i)grndEigVecs(nu,a)
    ! + grndEigVecs(nu,i)grndEigVecs(mu,a))
    ! complexity nOrb * nOrb * nOrb
    !
    ! xpycc(mu,nu) = sum_ia (X+Y)_ia grndEigVecs(mu,i) grndEigVecs(nu,a)
    ! xpycc(mu, nu) += sum_ia (X+Y)_ia grndEigVecs(mu,a) grndEigVecs(nu,i)

    call distributeRangeInChunks(env, 1, rpa%nxov_rd, iGlobal, fGlobal)
    xpycc(:,:,:) = 0.0_dp
    do ia = iGlobal, fGlobal
      call indxov(rpa%win, ia, rpa%getIA, i, a, iSpin)
      ! should replace with DSYR2 call :
      do nu = 1, nOrb
        do mu = 1, nOrb
          xpycc(mu,nu,iSpin) = xpycc(mu,nu,iSpin) + xpy(ia) *&
              & ( grndEigVecsGlobal(mu,i,iSpin)*grndEigVecsGlobal(nu,a,iSpin)&
              & + grndEigVecsGlobal(mu,a,iSpin)*grndEigVecsGlobal(nu,i,iSpin) )
        end do
      end do
    end do
    call assembleChunks(env, xpycc)


    if (rpa%tHybridXc) then

      ! Asymmetric contribution: xmycc_as = sum_ias (X-Y)_ias c_mas c_nis
      xmycc(:,:,:) = 0.0_dp
      xpyas(:,:,:) = 0.0_dp
      xmyas(:,:,:) = 0.0_dp
      call distributeRangeInChunks(env, 1, rpa%nxov_rd, iGlobal, fGlobal)
      do ia = iGlobal, fGlobal
        call indxov(rpa%win, ia, rpa%getIA, i, a, iSpin)
        ! should replace with DSYR2 call:
        do nu = 1, nOrb
          do mu = 1, nOrb
            xmycc(mu,nu,iSpin) = xmycc(mu,nu,iSpin) + xmy(ia) *&
                & ( grndEigVecsGlobal(mu,i,iSpin) * grndEigVecsGlobal(nu,a,iSpin)&
                & + grndEigVecsGlobal(mu,a,iSpin) * grndEigVecsGlobal(nu,i,iSpin) )
            xpyas(mu,nu,iSpin) = xpyas(mu,nu,iSpin) + xpy(ia) *&
                & grndEigVecsGlobal(mu,i,iSpin) * grndEigVecsGlobal(nu,a,iSpin)
            xmyas(mu,nu,iSpin) = xmyas(mu,nu,iSpin) + xmy(ia) *&
                & grndEigVecsGlobal(mu,i,iSpin) * grndEigVecsGlobal(nu,a,iSpin)
          end do
        end do
      end do
      call assembleChunks(env, xmycc)
      call assembleChunks(env, xpyas)
      call assembleChunks(env, xmyas)

      ! Account for normalization of S/T versus spin-polarized X+/-Y
      ! We have (X+Y)^S = 1/sqrt(2) [(X+Y)_up + (X+Y)_dn]
      if (lr%tSpin) then
        xmycc(:,:,:) = xmycc / sqrt(2._dp)
        xpyas(:,:,:) = xpyas / sqrt(2._dp)
        xmyas(:,:,:) = xmyas / sqrt(2._dp)
      end if

      do iSpin = 1, nSpin
        call symm(PS(:,:,iSpin), 'R', overlap, pc(:,:,iSpin), 'U', 1.0_dp, 0.0_dp, nOrb, nOrb)
        call symm(SPS(:,:,iSpin), 'L', overlap, PS(:,:,iSpin), 'U', 1.0_dp, 0.0_dp, nOrb, nOrb)
        call symm(DS(:,:,iSpin), 'R', overlap, deltaRhoGlobal(:,:,iSpin), 'U', 1.0_dp, 0.0_dp, nOrb, nOrb)
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
    ! complexity nOrb * nOrb * nOrb

    ! calculate the occ-occ part
    wcc(:,:,:) = 0.0_dp
    do iSpin = 1, nSpin
      call distributeRangeInChunks(env, 1, rpa%nxoo_ud(iSpin), iGlobal, fGlobal)
      do ij = iGlobal, fGlobal
        i = rpa%getIJ(ij + soo(iSpin), 1)
        j = rpa%getIJ(ij + soo(iSpin), 2)
        ! replace with DSYR2 call :
        do mu = 1, nOrb
          do nu = 1, nOrb
            wcc(mu,nu,iSpin) = wcc(mu,nu,iSpin) + woo(ij,iSpin) *&
                & ( grndEigVecsGlobal(mu,i,iSpin)*grndEigVecsGlobal(nu,j,iSpin)&
                & + grndEigVecsGlobal(mu,j,iSpin)*grndEigVecsGlobal(nu,i,iSpin) )
          end do
        end do

      end do
    end do

    ! calculate the occ-virt part : the same way as for xpycc
    call distributeRangeInChunks(env, 1, rpa%nxov_rd, iGlobal, fGlobal)
    do ia = iGlobal, fGlobal
      call indxov(rpa%win, ia, rpa%getIA, i, a, iSpin)
      ! again replace with DSYR2 call :
      do nu = 1, nOrb
        do mu = 1, nOrb
          wcc(mu,nu,iSpin) = wcc(mu,nu,iSpin) + wov(ia) *&
              & ( grndEigVecsGlobal(mu,i,iSpin)*grndEigVecsGlobal(nu,a,iSpin)&
              & + grndEigVecsGlobal(mu,a,iSpin)*grndEigVecsGlobal(nu,i,iSpin) )
        end do
      end do
    end do

    ! calculate the virt - virt part
    do iSpin = 1, nSpin
      call distributeRangeInChunks(env, 1, rpa%nxvv_ud(iSpin), iGlobal, fGlobal)
      do ab = iGlobal, fGlobal
        a = rpa%getAB(ab + svv(iSpin), 1)
        b = rpa%getAB(ab + svv(iSpin), 2)
        ! replace with DSYR2 call :
        do mu = 1, nOrb
          do nu = 1, nOrb
            wcc(mu,nu,iSpin) = wcc(mu,nu,iSpin) + wvv(ab,iSpin) *&
                & ( grndEigVecsGlobal(mu,a,iSpin)*grndEigVecsGlobal(nu,b,iSpin)&
                & + grndEigVecsGlobal(mu,b,iSpin)*grndEigVecsGlobal(nu,a,iSpin) )
          end do
        end do
      end do
    end do
    call assembleChunks(env, wcc)

    ! now calculating the force complexity : nOrb * nOrb * 3

    ! as have already performed nOrb**3 operation to get here,
    ! calculate for all atoms

    ! BA: only for non-periodic systems!
    do iAt1 = 1, lr%nAtom
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
        dgab(:) = diffvec(:) * (-1.0_dp/rab**2 - expGammaPrime(rab, lr%HubbardU(iSp1),&
            & lr%HubbardU(iSp2)))

        tmp3a = 0.0_dp
        do iSpin = 1, nSpin
          tmp3a = tmp3a + dq(iAt1) * dqex(iAt2,iSpin) + dqex(iAt1,iSpin) * dq(iAt2)
        end do

        if (.not. lr%tSpin) then
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

        if (lr%tSpin) then
          tmp9 = lr%spinW(iSp1) * dm(iAt1) + lr%spinW(iSp2) * dm(iAt2)
          tmp11(:) = lr%spinW(iSp1) * dqex(iAt1,:) + lr%spinW(iSp2) * dqex(iAt2,:)
        end if

        if (rpa%tHybridXc) then
          tmprs = 0.0_dp
          tmps(:) = 0.0_dp
          do iSpin = 1, nSpin
            do mu = indAlpha, indAlpha1
              do nu = indBeta, indBeta1
                tmprs = tmprs +&
          & ( 2.0_dp * (PS(mu,nu,iSpin) * DS(nu,mu,iSpin) + PS(nu,mu,iSpin) * DS(mu,nu,iSpin)) +&
          &   SPS(mu,nu,iSpin) * deltaRhoGlobal(mu,nu,iSpin) + SPS(nu,mu,iSpin) * deltaRhoGlobal(nu,mu,iSpin) +&
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

                if (lr%tSpin) then
                  tmp8 = tmp8 + tmp9 * dSo(n,m,xyz) * dsigma(iSpin) * pc(mu,nu,iSpin)
                  tmp10 = tmp10 + tmp11(iSpin) * dSo(n,m,xyz) * dsigma(iSpin) * SpinDens(mu,nu)
                end if

                if (rpa%tHybridXc) then
                  tmprs = 0.0_dp
                  do ka = 1, nOrb
                    tmprs = tmprs +&
            & ( PS(mu,ka,iSpin) * deltaRhoGlobal(nu,ka,iSpin) + PS(nu,ka,iSpin) * deltaRhoGlobal(mu,ka,iSpin) +&
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
  subroutine writeExcitations(lr, rpa, sym, osz, eval, xpy, fdXPlusY, fdTrans, fdTransDip,&
      & transitionDipoles, fdTagged, taggedWriter, fdExc, Ssq)

    !> Data structure for linear response
    type(TLinResp), intent(in) :: lr

    !> Run time parameters of the Casida routine
    type(TCasidaParameter), intent(in) :: rpa

    !> Symmetry label for the type of transition
    character, intent(in) :: sym

    !> Oscillator strengths for transitions from ground to excited states
    real(dp), intent(in) :: osz(:)

    !> Excitation energies
    real(dp), intent(in) :: eval(:)

    !> Eigenvectors of excited states (X+Y)
    real(dp), intent(in) :: xpy(:,:)

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

    integer, allocatable :: wvin(:), degenerate(:,:)
    integer :: nmat, ii, jj, iweight, indo, m, n, s
    real(dp), allocatable :: wvec(:), oDeg(:)
    real(dp) :: weight, wvnorm
    logical :: updwn, tSpin,tDegenerate
    character :: sign
    type(TDegeneracyFind) :: DegeneracyFind

    tSpin = present(Ssq)
    nmat = size(rpa%wij)

    allocate(wvec(nmat))
    allocate(wvin(nmat))
    wvec(:) = 0.0_dp
    wvin(:) = 0

    if (fdXplusY%isConnected()) then
      write(fdXPlusY%unit, *) nmat, lr%nExc
    end if

    do ii = 1, lr%nExc
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

        call indxov(rpa%win, iweight, rpa%getIA, m, n, s)
        sign = sym
        if (fdExc%isConnected()) then
          if (lr%tSpin) then
            sign = " "
            write(fdExc%unit,&
                & '(1x,f10.3,4x,f14.8,2x,i5,3x,a,1x,i5,7x,f6.3,2x,f10.3,4x,&
                & f6.3)')&
                & Hartree__eV * sqrt(eval(ii)), osz(ii), m, '->', n, weight,&
                & Hartree__eV * rpa%wij(iWeight), Ssq(ii)
          else
            write(fdExc%unit,&
                & '(1x,f10.3,4x,f14.8,5x,i5,3x,a,1x,i5,7x,f6.3,2x,f10.3,6x,a)')&
                & Hartree__eV * sqrt(eval(ii)), osz(ii), m, '->', n, weight,&
                & Hartree__eV * rpa%wij(iWeight), sign
          end if
        end if

        if (fdXplusY%isConnected()) then
          if (tSpin) then
            updwn = (rpa%win(iweight) <= rpa%nxov_ud(1))
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
            call indxov(rpa%win, indo, rpa%getIA, m, n, s)
            if (tSpin) then
              updwn = (rpa%win(indo) <= rpa%nxov_ud(1))
              sign = "D"
              if (updwn) sign = "U"
            end if
            write(fdTrans%unit, '(i5,3x,a,1x,i5,1x,1a,T22,f10.8,T33,f14.8)')&
                & m, '->', n, sign, wvec(jj), Hartree__eV * rpa%wij(wvin(jj))
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
        call indxov(rpa%win, iWeight, rpa%getIA, m, n, s)
        sign = sym

        if (fdExc%isConnected()) then
          if (lr%tSpin) then
            sign = " "
            write(fdExc%unit,&
                & '(6x,A,T12,4x,f14.8,2x,i5,3x,a,1x,i5,7x,A,2x,f10.3,4x,f6.3)')&
                & '< 0', osz(ii), m, '->', n, '-', Hartree__eV * rpa%wij(iWeight),&
                & Ssq(ii)
          else
            write(fdExc%unit,&
                & '(6x,A,T12,4x,f14.8,2x,i5,3x,a,1x,i5,7x,f6.3,2x,f10.3,6x,a)')&
                & '< 0', osz(ii), m, '->', n, weight, Hartree__eV * rpa%wij(iWeight), sign
          end if
        end if

        if (fdXplusY%isConnected()) then
          if (lr%tSpin) then
            updwn = (rpa%win(iweight) <= rpa%nxov_ud(1))
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
  subroutine calcPMatrix(env, rpa, t, rhs, pc)

    !> Environment settings
    type(TEnvironment), intent(inout) :: env

    !> Run time parameters of the Casida routine
    type(TCasidaParameter), intent(in) :: rpa

    !> T matrix
    real(dp), intent(in) :: t(:,:,:)

    !> Z matrix
    real(dp), intent(in) :: rhs(:)

    !> Resulting excited state density matrix
    real(dp), intent(out) :: pc(:,:,:)

    integer :: ias, i, a, s, nSpin, iGlobal, fGlobal

    nSpin = size(pc, dim=3)

    call distributeRangeInChunks(env, 1, size(rhs), iGlobal, fGlobal)

    pc(:,:,:) = 0.0_dp
    do ias = iGlobal, fGlobal
      call indxov(rpa%win, ias, rpa%getIA, i, a, s)
      pc(i,a,s) = rhs(ias)
    end do
    call assembleChunks(env, pc)

    do s = 1, nSpin
      pc(:,:,s) = 0.5_dp * ( pc(:,:,s) + transpose(pc(:,:,s)) )
    end do

    pc(:,:,:) = pc + t

  end subroutine calcPMatrix


  !> Computes H^+/-_pq [V] as defined in Furche JCP 117 7433 (2002) eq. 20
  !! Here p/q are virtual orbitals and V is either X+Y or X-Y
  subroutine getHvvXY(env, orb, lr, rpa, transChrg, denseDesc, ovrXev, grndEigVecs, lrGamma, ipm, XorY,&
      & vecHvv)

    !> Environment settings
    type(TEnvironment), intent(inout) :: env

    !> Data type for atomic orbital information
    type(TOrbitals), intent(in) :: orb

    !> Data structure for linear response
    type(TLinResp), intent(in) :: lr

    !> Run time parameters of the Casida routine
    type(TCasidaParameter), intent(in) :: rpa

    !> machinery for transition charges between single particle levels
    type(TTransCharges), intent(in) :: transChrg

    !> Dense matrix descriptor
    type(TDenseDescr), intent(in) :: denseDesc

    !> overlap times eigenvector. (nOrb, nOrb) [distributed]
    real(dp), intent(in) :: ovrXev(:,:,:)

    !> eigenvectors (nOrb, nOrb) [distributed]
    real(dp), intent(in) :: grndEigVecs(:,:,:)

    !> long-range Gamma if in use
    real(dp), allocatable, intent(in) :: lrGamma(:,:)

    !> Sign s of H in H^(s)[V]
    integer, intent(in) :: ipm

    !> RPA eigenvectors, either (X+Y) or (X-Y)
    real(dp), intent(in) :: XorY(:)

    !> Output vector H[V] virtual-virtual
    real(dp), intent(out) :: vecHvv(:)

    real(dp), allocatable :: qIJ(:), gqIJ(:), qX(:,:), Gq(:,:)
    integer :: i, a, b, s, ias, ibs, abs, nOrb, iGlobal, fGlobal

    nOrb = orb%nOrb

    allocate(qIJ(lr%nAtom))
    allocate(gqIJ(lr%nAtom))
    allocate(qX(lr%nAtom, rpa%nxov_rd))
    allocate(Gq(lr%nAtom, rpa%nxov_rd))

    call distributeRangeInChunks(env, 1, rpa%nxov_rd, iGlobal, fGlobal)

    qX(:,:) = 0.0_dp
    do ias = iGlobal, fGlobal
      call indXov(rpa%win, ias, rpa%getIA, i, a, s)
      do b = rpa%nocc_ud(s) + 1, nOrb
        ibs = rpa%iaTrans(i, b, s)
        abs = rpa%iaTrans(a, b, s)
        qIJ(:) = transChrg%qTransAB(abs, env, denseDesc, ovrXev, grndEigVecs, rpa%getAB)
        qX(:,ias) = qX(:,ias) + qIJ * XorY(ibs)
      end do
    end do

    Gq(:,:) = 0.0_dp
    do ias = iGlobal, fGlobal
      qIJ(:) = transChrg%qTransIA(ias, env, denseDesc, ovrXev, grndEigVecs, rpa%getIA, rpa%win)
      call hemv(gqIJ, lrGamma, qIJ, uplo='U')
      Gq(:,ias) = gqIJ(:)
    end do

    call assembleChunks(env,qX)
    call assembleChunks(env,Gq)

    call distributeRangeInChunks(env, 1, sum(rpa%nxvv_ud), iGlobal, fGlobal)

    vecHvv(:) = 0.0_dp
    do abs = iGlobal, fGlobal
      a = rpa%getAB(abs, 1)
      b = rpa%getAB(abs, 2)
      s = rpa%getAB(abs, 3)
      do i = 1, rpa%nocc_ud(s)
        ias = rpa%iaTrans(i, a, s)
        ibs = rpa%iaTrans(i, b, s)
        vecHvv(abs) = vecHvv(abs) - ipm * (dot_product(qX(:,ias), Gq(:,ibs))&
            & + ipm * dot_product(Gq(:,ias), qX(:,ibs)))
      end do
    end do

    call assembleChunks(env,vecHvv)

  end subroutine getHvvXY


  !> Computes H^+/-_pq [V] as defined in Furche JCP 117 7433 (2002) eq. 20
  !! Here p/q are occupied orbitals and V is either X+Y or X-Y
  subroutine getHooXY(env, orb, lr, rpa, transChrg, denseDesc, ovrXev, grndEigVecs, lrGamma, ipm, XorY,&
      & vecHoo)

    !> Environment settings
    type(TEnvironment), intent(inout) :: env

    !> Data type for atomic orbital information
    type(TOrbitals), intent(in) :: orb

    !> Data structure for linear response
    type(TLinResp), intent(in) :: lr

    !> Run time parameters of the Casida routine
    type(TCasidaParameter), intent(in) :: rpa

    !> Machinery for transition charges between single particle levels
    type(TTransCharges), intent(in) :: transChrg

    !> Dense matrix descriptor
    type(TDenseDescr), intent(in) :: denseDesc

    !> Overlap times eigenvector. (nOrb, nOrb) [distributed]
    real(dp), intent(in) :: ovrXev(:,:,:)

    !> Eigenvectors (nOrb, nOrb) [distributed]
    real(dp), intent(in) :: grndEigVecs(:,:,:)

    !> Long-range Gamma if in use
    real(dp), allocatable, intent(in) :: lrGamma(:,:)

    !> Sign s of H in H^(s)[V]
    integer, intent(in) :: ipm

    !> RPA eigenvectors, either (X+Y) or (X-Y)
    real(dp), intent(in) :: XorY(:)

    !> Output vector H[V] occ-occ
    real(dp), intent(out) :: vecHoo(:)

    real(dp), allocatable :: qIJ(:), gqIJ(:), qX(:,:), Gq(:,:)
    integer :: i, j, a, s, ias, jas, ijs, nOrb, iGlobal, fGlobal

    nOrb = orb%nOrb

    allocate(qIJ(lr%nAtom))
    allocate(gqIJ(lr%nAtom))
    allocate(qX(lr%nAtom, rpa%nxov_rd))
    allocate(Gq(lr%nAtom, rpa%nxov_rd))

    call distributeRangeInChunks(env, 1, rpa%nxov_rd, iGlobal, fGlobal)

    qX(:,:) = 0.0_dp
    do ias = iGlobal, fGlobal
      call indXov(rpa%win, ias, rpa%getIA, i, a, s)
      do j = 1, rpa%nocc_ud(s)
        jas = rpa%iaTrans(j, a, s)
        ijs = rpa%iaTrans(i, j, s)
        qIJ(:) = transChrg%qTransIJ(ijs, env, denseDesc, ovrXev, grndEigVecs, rpa%getIJ)
        qX(:,ias) = qX(:,ias) + qIJ * XorY(jas)
      end do
    end do

    Gq(:,:) = 0.0_dp
    do ias = iGlobal, fGlobal
      call indXov(rpa%win, ias, rpa%getIA, i, a, s)
      qIJ(:) = transChrg%qTransIA(ias, env, denseDesc, ovrXev, grndEigVecs, rpa%getIA, rpa%win)
      call hemv(gqIJ, lrGamma, qIJ, uplo='U')
      Gq(:,ias) = gqIJ
    end do

    call assembleChunks(env,qX)
    call assembleChunks(env,Gq)

    call distributeRangeInChunks(env, 1, sum(rpa%nxoo_ud), iGlobal, fGlobal)

    vecHoo(:) = 0.0_dp
    do ijs = iGlobal, fGlobal
      i = rpa%getIJ(ijs, 1)
      j = rpa%getIJ(ijs, 2)
      s = rpa%getIJ(ijs, 3)
      do a = rpa%nocc_ud(s) + 1, nOrb
        ias = rpa%iaTrans(i, a, s)
        jas = rpa%iaTrans(j, a, s)
        vecHoo(ijs) = vecHoo(ijs) - ipm * (dot_product(qX(:,ias), Gq(:,jas))&
            & + ipm * dot_product(Gq(:,ias), qX(:,jas)))
      end do
    end do

    call assembleChunks(env,vecHoo)

  end subroutine getHooXY


  !> Computes H^+/-_pq [T] as defined in Furche JCP 117 7433 (2002) eq. 20
  !! Here p is an occupied MO and q is a virtual one, T is the relaxed difference density
  subroutine getHovT(env, orb, lr, rpa, transChrg, denseDesc, ovrXev, grndEigVecs, lrGamma, t, vecHovT)

    !> Environment settings
    type(TEnvironment), intent(inout) :: env

    !> Data type for atomic orbital information
    type(TOrbitals), intent(in) :: orb

    !> Data structure for linear response
    type(TLinResp), intent(in) :: lr

    !> Run time parameters of the Casida routine
    type(TCasidaParameter), intent(in) :: rpa

    !> Machinery for transition charges between single particle levels
    type(TTransCharges), intent(in) :: transChrg

    !> Dense matrix descriptor
    type(TDenseDescr), intent(in) :: denseDesc

    !> Overlap times eigenvector. (nOrb, nOrb) [distributed]
    real(dp), intent(in) :: ovrXev(:,:,:)

    !> Eigenvectors (nOrb, nOrb) [distributed]
    real(dp), intent(in) :: grndEigVecs(:,:,:)

    !> Long-range Gamma if in use
    real(dp), allocatable, intent(in) :: lrGamma(:,:)

    !> Excited state density matrix
    real(dp), intent(in) :: t(:,:,:)

    !> Output vector H[T] occ-vir
    real(dp), intent(out) :: vecHovT(:)

    real(dp), allocatable :: qIJ(:), gqIJ(:), qX(:,:), Gq(:,:)
    integer :: i, j, a, b, s, ias, ibs, abs, ijs, jas, nOrb, iMx, iGlobal, fGlobal

    nOrb = orb%nOrb

    allocate(qIJ(lr%nAtom))
    allocate(gqIJ(lr%nAtom))
    allocate(qX(lr%nAtom, rpa%nxov_rd))
    iMx = max(sum(rpa%nxoo_ud), sum(rpa%nxvv_ud))
    allocate(Gq(lr%nAtom, iMx))

    call distributeRangeInChunks(env, 1, rpa%nxov_rd, iGlobal, fGlobal)

    qX(:,:) = 0.0_dp
    do ias = iGlobal, fGlobal
      call indXov(rpa%win, ias, rpa%getIA, i, a, s)
      do b = rpa%nocc_ud(s) + 1, nOrb
        ibs = rpa%iaTrans(i, b, s)
        qIJ(:) = transChrg%qTransIA(ibs, env, denseDesc, ovrXev, grndEigVecs, rpa%getIA, rpa%win)
        qX(:,ias) = qX(:,ias) + qIJ * t(a,b,s)
      end do
    end do

    call assembleChunks(env,qX)
    call distributeRangeInChunks(env, 1, sum(rpa%nxvv_ud), iGlobal, fGlobal)

    Gq(:,:) = 0.0_dp
    do abs = iGlobal, fGlobal
      qIJ(:) = transChrg%qTransAB(abs, env, denseDesc, ovrXev, grndEigVecs, rpa%getAB)
      call hemv(gqIJ, lrGamma, qIJ, uplo='U')
      Gq(:,abs) = gqIJ
    end do

    call assembleChunks(env,Gq)
    call distributeRangeInChunks(env, 1, rpa%nxov_rd, iGlobal, fGlobal)

    vecHovT(:) = 0.0_dp
    do ias = iGlobal, fGlobal
      call indXov(rpa%win, ias, rpa%getIA, i, a, s)
      do b = rpa%nocc_ud(s) + 1, nOrb
        ibs = rpa%iaTrans(i, b, s)
        abs = rpa%iaTrans(a, b, s)
        vecHovT(ias) = vecHovT(ias) - 2.0_dp * dot_product(qX(:,ibs), Gq(:,abs))
      end do
    end do

    qX(:,:) = 0.0_dp
    do ias = iGlobal, fGlobal
      call indXov(rpa%win, ias, rpa%getIA, i, a, s)
      do j = 1, rpa%nocc_ud(s)
        jas = rpa%iaTrans(j, a, s)
        qIJ(:) = transChrg%qTransIA(jas, env, denseDesc, ovrXev, grndEigVecs, rpa%getIA, rpa%win)
        qX(:,ias) = qX(:,ias) + qIJ * t(i,j,s)
      end do
    end do

    call assembleChunks(env,qX)
    call distributeRangeInChunks(env, 1, sum(rpa%nxoo_ud), iGlobal, fGlobal)

    Gq(:,:) = 0.0_dp
    do ijs = iGlobal, fGlobal
      i = rpa%getIJ(ijs, 1)
      j = rpa%getIJ(ijs, 2)
      s = rpa%getIJ(ijs, 3)
      qIJ(:) = transChrg%qTransIJ(ijs, env, denseDesc, ovrXev, grndEigVecs, rpa%getIJ)
      call hemv(gqIJ, lrGamma, qIJ, uplo='U')
      Gq(:,ijs) = gqIJ
    end do

    call assembleChunks(env,Gq)
    call distributeRangeInChunks(env, 1, rpa%nxov_rd, iGlobal, fGlobal)

    do ias = iGlobal, fGlobal
      call indXov(rpa%win, ias, rpa%getIA, i, a, s)
      do j = 1, rpa%nocc_ud(s)
        jas = rpa%iaTrans(j, a, s)
        ijs = rpa%iaTrans(i, j, s)
        vecHovT(ias) = vecHovT(ias) - 2.0_dp * dot_product(qX(:,jas), Gq(:,ijs))
      end do
    end do

    call assembleChunks(env,vecHovT)

  end subroutine getHovT


  !> Computes H^+/-_pq [T] as defined in Furche JCP 117 7433 (2002) eq. 20
  !! Here p/q are occupied MO, T is the relaxed difference density
  subroutine getHooT(env, orb, lr, rpa, transChrg, denseDesc, ovrXev, grndEigVecs, lrGamma, t, vecHooT)

    !> Environment settings
    type(TEnvironment), intent(inout) :: env

    !> Data type for atomic orbital information
    type(TOrbitals), intent(in) :: orb

    !> Data structure for linear response
    type(TLinResp), intent(in) :: lr

    !> Run time parameters of the Casida routine
    type(TCasidaParameter), intent(in) :: rpa

    !> Machinery for transition charges between single particle levels
    type(TTransCharges), intent(in) :: transChrg

    !> Dense matrix descriptor
    type(TDenseDescr), intent(in) :: denseDesc

    !> Overlap times eigenvector. (nOrb, nOrb) [distributed]
    real(dp), intent(in) :: ovrXev(:,:,:)

    !> Eigenvectors (nOrb, nOrb) [distributed]
    real(dp), intent(in) :: grndEigVecs(:,:,:)

    !> Long-range Gamma if in use
    real(dp), allocatable, intent(in) :: lrGamma(:,:)

    !> Excited state density matrix
    real(dp), intent(in) :: t(:,:,:)

    !> Output vector H[T] occ-occ
    real(dp), intent(out) :: vecHooT(:)

    real(dp), allocatable :: qIJ(:), gqIJ(:), qX(:,:), Gq(:,:), qXa(:,:,:)
    integer :: nOrb, iSpin, nSpin, iMx, soo(2)
    integer :: i, j, k, a, b, s, ij, ias, ibs, ijs, jas, iks, jks, iGlobal, fGlobal

    nOrb = orb%nOrb
    nSpin = size(t, dim=3)
    soo(:) = [0, rpa%nxoo_ud(1)]

    allocate(qIJ(lr%nAtom))
    allocate(gqIJ(lr%nAtom))
    iMx = max(sum(rpa%nxoo_ud), rpa%nxov_rd)
    allocate(qX(lr%nAtom, iMx))
    allocate(Gq(lr%nAtom, iMx))

    call distributeRangeInChunks(env, 1, rpa%nxov_rd, iGlobal, fGlobal)

    qX(:,:) = 0.0_dp
    do ias = iGlobal, fGlobal
      call indXov(rpa%win, ias, rpa%getIA, i, a, s)
      do b = rpa%nocc_ud(s) + 1, nOrb
        ibs = rpa%iaTrans(i, b, s)
        qIJ(:) = transChrg%qTransIA(ibs, env, denseDesc, ovrXev, grndEigVecs, rpa%getIA, rpa%win)
        qX(:,ias) = qX(:,ias) + qIJ * t(a,b,s)
      end do
    end do

    Gq(:,:) = 0.0_dp
    do ias = iGlobal, fGlobal
      qIJ(:) = transChrg%qTransIA(ias, env, denseDesc, ovrXev, grndEigVecs, rpa%getIA, rpa%win)
      call hemv(gqIJ, lrGamma, qIJ, uplo='U')
      Gq(:,ias) = gqIJ
    end do

    call assembleChunks(env,qX)
    call assembleChunks(env,Gq)

    call distributeRangeInChunks(env, 1, sum(rpa%nxoo_ud), iGlobal, fGlobal)

    vecHooT(:) = 0.0_dp
    do ijs = iGlobal, fGlobal
      i = rpa%getIJ(ijs, 1)
      j = rpa%getIJ(ijs, 2)
      s = rpa%getIJ(ijs, 3)
      do a = rpa%nocc_ud(s) + 1, nOrb
        ias = rpa%iaTrans(i, a, s)
        jas = rpa%iaTrans(j, a, s)
        vecHooT(ijs) = vecHooT(ijs) - 2.0_dp * dot_product(qX(:,ias), Gq(:,jas))
      end do
    end do

    deallocate(qX)

    Gq(:,:) = 0.0_dp
    do ijs = iGlobal, fGlobal
      qIJ = transChrg%qTransIJ(ijs, env, denseDesc, ovrXev, grndEigVecs, rpa%getIJ)
      call hemv(gqIJ, lrGamma, qIJ, uplo='U')
      Gq(:,ijs) = gqIJ(:)
    end do

    call assembleChunks(env,Gq)

    ! For qXa_ijs = sum_k q_iks t(j,k,s), we need both qXa_ijs and qXa_jis
    ! Need for a spin loop, don't think this can be simplified
    do iSpin = 1, nSpin

      call distributeRangeInChunks(env, 1, rpa%nocc_ud(iSpin), iGlobal, fGlobal)

      allocate(qXa(lr%nAtom, rpa%nocc_ud(iSpin), rpa%nocc_ud(iSpin)))
      qXa(:,:,:) = 0.0_dp
      do i = iGlobal, fGlobal
        do k = 1, rpa%nocc_ud(iSpin)
          iks = rpa%iaTrans(i, k, iSpin)
          qIJ(:) = transChrg%qTransIJ(iks, env, denseDesc, ovrXev, grndEigVecs, rpa%getIJ)
          do j = 1, rpa%nocc_ud(iSpin)
            qXa(:,i,j) = qXa(:,i,j) + qIJ * t(j,k,iSpin)
          end do
        end do
      end do

      call assembleChunks(env,qXa)

      call distributeRangeInChunks(env, 1, rpa%nxoo_ud(iSpin), iGlobal, fGlobal)

      do ij = iGlobal, fGlobal
        i = rpa%getIJ(ij + soo(iSpin), 1)
        j = rpa%getIJ(ij + soo(iSpin), 2)
        ijs = rpa%iaTrans(i, j, iSpin)
        do k = 1, rpa%nocc_ud(iSpin)
          jks = rpa%iaTrans(j, k, iSpin)
          vecHooT(ijs) = vecHooT(ijs) - 2.0_dp * dot_product(qXa(:,i,k), Gq(:,jks))
        end do
      end do
      deallocate(qXa)

    end do

    call assembleChunks(env,vecHooT)

  end subroutine getHooT


  !> Constructs the full overlap matrix S.
  subroutine getSqrS(coord, nAtom, skOverCont, orb, iAtomStart, species0, S)

    !> Atom coordinates
    real(dp), intent(in) :: coord(:,:)

    !> Number of atoms
    integer,intent(in) :: nAtom

    !> Starting position of each atom in the list of orbitals
    integer,intent(in) :: iAtomStart(:)

    !> chemical species of the atoms
    integer,intent(in) :: species0(:)

    !> Overlap data
    type(TSlakoCont), intent(in) :: skOverCont

    !> Data type for atomic orbital information
    type(TOrbitals), intent(in) :: orb

    !> Overlap matrix
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


  !> Constructs a Gamma-Matrix of dimension nOrb instead of nAtoms (i.e., atomic orbitals).
  subroutine getSqrGamma(nAtom, lrGamma, iAtomStart, lrGammaOrb)

    !> Long-range Gamma
    real(dp), intent(in) :: lrGamma(:,:)

    !> Number of atoms
    integer,intent(in) :: nAtom

    !> Starting position of each atom in the list of orbitals
    integer,intent(in) :: iAtomStart(:)

    !> Resulting gamma matrix
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


  !> Helper routine to construct diatomic overlap.
  subroutine getSOffsite(coords1, coords2, iSp1, iSp2, orb, skOverCont, Sblock)

    !> First atom coordinates
    real(dp), intent(in) :: coords1(:)

    !> Second atom coordinates
    real(dp), intent(in) :: coords2(:)

    !> First atom species
    integer, intent(in) :: iSp1

    !> Second atom species
    integer, intent(in) :: iSp2

    !> Data type for atomic orbital information
    type(TOrbitals), intent(in) :: orb

    !> Overlap data
    type(TSlakoCont), intent(in) :: skOverCont

    !> Diatomic block
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
  subroutine grndToExcDensityMatrices(env, orb, lr, rpa, transChrg, denseDesc, sym, species0, ovrXev,&
      & grndEigVecs, grndEigVal, frGamma, lrGamma, omega, pc, xpy, xmy, wov, woo)

    !> Environment settings
    type(TEnvironment), intent(inout) :: env

    !> Data type for atomic orbital information
    type(TOrbitals), intent(in) :: orb

    !> Data structure for linear response
    type(TLinResp), intent(in) :: lr

    !> Run time parameters of the Casida routine
    type(TCasidaParameter), intent(in) :: rpa

    !> Machinery for transition charges between single particle levels
    type(TTransCharges), intent(in) :: transChrg

    !> Dense matrix descriptor
    type(TDenseDescr), intent(in) :: denseDesc

    !> Symmetry of the transition
    character, intent(in) :: sym

    !> Central cell chemical species
    integer, intent(in) :: species0(:)

    !> Overlap times ground state eigenvectors
    real(dp), intent(in) :: ovrXev(:,:,:)

    !> Ground state eigenvectors
    real(dp), intent(in) :: grndEigVecs(:,:,:)

    !> Ground state wavefunctions
    real(dp), intent(in) :: grndEigVal(:,:)

    !> Softened coulomb matrix
    real(dp), intent(in) :: frGamma(:,:)

    !> Electrostatic matrix, long-range corrected
    real(dp), allocatable, intent(in) :: lrGamma(:,:)

    !> Excitation energy of states n
    real(dp), intent(in) :: omega

    !> P matrix (symmetric)
    real(dp), intent(out) :: pc(:,:,:)

    !> X+Y Furche term for excited state n
    real(dp), intent(in) :: xpy(:)

    !> X-Y Furche term for excited state n
    real(dp), intent(in) :: xmy(:)

    !> W^+ vector occupied-virtual part
    real(dp), intent(out) :: wov(:)

    !> W^+ vector occupied part
    real(dp), intent(out) :: woo(:,:)

    real(dp), allocatable :: p(:), vecHoo(:)
    integer :: soo(2), i, a, s, ias, j, ij, ijs, nSpin, iGlobal, fGlobal


    nSpin = size(rpa%nocc_ud)

    allocate(p(rpa%nxov_rd))
    allocate(vecHoo(sum(rpa%nxoo_ud)))

    ! Transition charges use compound index ijs = ij + soo(s)
    soo(:) = [0, rpa%nxoo_ud(1)]

    p(:) = 0.0_dp
    wov(:) = 0.0_dp
    woo(:,:) = 0.0_dp

    ! MPI distribution indexes
    call distributeRangeInChunks(env, 1, rpa%nxov_rd, iGlobal, fGlobal)

    ! "Fake" density matrix for non-adiabatic coupling [Furche JCP 132 044107 (2010)]
    ! Restricted KS: P = 2 P^up ; (X+Y) = sqrt(2) (X+Y)^up
    do ias = iGlobal, fGlobal
      call indxov(rpa%win, ias, rpa%getIA, i, a, s)
      p(ias) = sqrt(2.0_dp) * xpy(ias)
      wov(ias) = grndEigVal(i, s) * p(ias) + omega * xmy(ias) / sqrt(2.0_dp)
    end do

    ! Define P symmetrically (similar to treatment of excited state gradients)
    pc(:,:,:) = 0.0_dp
    do ias = iGlobal, fGlobal
      call indxov(rpa%win, ias, rpa%getIA, i, a, s)
      pc(i,a,s) = 0.5_dp * p(ias)
      pc(a,i,s) = 0.5_dp * p(ias)
    end do

    call assembleChunks(env,p)
    call assembleChunks(env,pc)
    call assembleChunks(env,wov)

    call getHplusXYfr(env, lr, rpa, transChrg, denseDesc, sym, species0, ovrXev, grndEigVecs,&
        & frGamma, p, vecHoo=vecHoo)

    do s = 1, nSpin
      call distributeRangeInChunks(env, 1, rpa%nxoo_ud(s), iGlobal, fGlobal)
      do ij = iGlobal, fGlobal
        i = rpa%getIJ(ij + soo(s), 1)
        j = rpa%getIJ(ij + soo(s), 2)
        ijs = rpa%iaTrans(i, j, s)
        ! GetHplusXYfr used with P instead of X+Y yields half the desired result
        woo(ij,s) = 2.0_dp * vecHoo(ijs)
      end do
    end do

    if (rpa%tHybridXc) then
      call getHooXY(env, orb, lr, rpa, transChrg, denseDesc, ovrXev, grndEigVecs, lrGamma, 1, p, vecHoo)

      do s = 1, nSpin
        call distributeRangeInChunks(env, 1, rpa%nxoo_ud(s), iGlobal, fGlobal)
        do ij = iGlobal, fGlobal
          i = rpa%getIJ(ij + soo(s), 1)
          j = rpa%getIJ(ij + soo(s), 2)
          ijs = rpa%iaTrans(i, j, s)
          woo(ij,s) = woo(ij,s) + vecHoo(ijs)
        end do
      end do
    end if

    ! Divide diagonal elements of W_ij by 2.
    do s = 1, nSpin
      call distributeRangeInChunks(env, 1, rpa%nxoo_ud(s), iGlobal, fGlobal)
      do ij = iGlobal, fGlobal
        i = rpa%getIJ(ij + soo(s), 1)
        j = rpa%getIJ(ij + soo(s), 2)
        if (i == j) then
          woo(ij,s) = 0.5_dp * woo(ij,s)
        end if
      end do
    end do
    call assembleChunks(env,woo)

  end subroutine grndToExcDensityMatrices


  !> Build right hand side of the equation for the Z-vector and those parts of the W-vectors which
  !! do not depend on Z. Modified version of getZVectorEqRHS for state-to-state NA couplings.
  !! Furche PCCP 21 18999 (2019)
  !! Here the + (symmetric) part of RHS, T and (omega_m-omega_n) * W (stored as W) is computed.
  subroutine getNadiaZVectorEqRHS(env, orb, lr, rpa, transChrg, sym, denseDesc, species0, grndEigVal,&
      & ovrXev, grndEigVecs, gammaMat, lrGamma, omegaAvg, xpyn, xmyn, xpym, xmym, rhs, t, wov, woo,&
      & wvv)

    !> Environment settings
    type(TEnvironment), intent(inout) :: env

    !> Data type for atomic orbital information
    type(TOrbitals), intent(in) :: orb

    !> Data structure for linear response
    type(TLinResp), intent(in) :: lr

    !> Run time parameters of the Casida routine
    type(TCasidaParameter), intent(in) :: rpa

    !> Machinery for transition charges between single particle levels
    type(TTransCharges), intent(in) :: transChrg

    !> Symmetry of the transitions
    character, intent(in) :: sym

    !> Dense matrix descriptor
    type(TDenseDescr), intent(in) :: denseDesc

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

    !> Average excitation energy of states n and m
    real(dp), intent(in) :: omegaAvg

    !> X+Y Furche term for excited state n
    real(dp), intent(in) :: xpyn(:)

    !> X-Y Furche term for excited state n
    real(dp), intent(in) :: xmyn(:)

    !> X+Y Furche term for excited state m
    real(dp), intent(in) :: xpym(:)

    !> X-Y Furche term for excited state m
    real(dp), intent(in) :: xmym(:)

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
    real(dp), allocatable :: vecHvvXorY(:), vecHooXorY(:), vecHovT(:), vecHvvT(:)
    real(dp), allocatable :: vecHvvXpY(:), vecHvvXmY(:), vecHooXpY(:), vecHooXmY(:)
    real(dp), allocatable :: vecHooT(:)
    integer :: nxov, iGlobal, fGlobal
    integer :: i, j, a, b, ias, ibs, abs, ij, ab, jas, ijs, s, nSpin, soo(2), svv(2), nOrb
    real(dp) :: ptmp1, ptmp2, tmp3, tmp4, tmp1

    nxov = size(rhs)
    nOrb = orb%nOrb

    allocate(xpyq(lr%nAtom))
    allocate(qTr(lr%nAtom))
    allocate(gamxpyq(lr%nAtom))
    allocate(gamqt(lr%nAtom))

    t(:,:,:) = 0.0_dp
    rhs(:) = 0.0_dp
    wov(:) = 0.0_dp
    woo(:,:) = 0.0_dp
    wvv(:,:) = 0.0_dp

    nSpin = size(t, dim=3)

    ! Transition charges use compound index ijs = ij + soo(s)
    soo(:) = [0, rpa%nxoo_ud(1)]
    svv(:) = [0, rpa%nxvv_ud(1)]

    allocate(qgamxpyq(max(maxval(rpa%nxoo_ud), maxval(rpa%nxvv_ud)), size(rpa%nocc_ud)))
    allocate(vecHooXorY(sum(rpa%nxoo_ud)))
    allocate(vecHvvXorY(sum(rpa%nxvv_ud)))

    if (nSpin == 2) then
      allocate(xpyqds(lr%nAtom))
      allocate(gamxpyqds(lr%nAtom))
    end if

    call distributeRangeInChunks(env, 1, nxov, iGlobal, fGlobal)

    ! Build state-to-state 1TDM and W (eq. 42 in Furche PCCP)
    ! We are symmetrizing the non-symmetric T of Furche
    ! Factor 1 / (1 + delta_ab) for W follows later
    do ias = iGlobal, fGlobal
      call indxov(rpa%win, ias, rpa%getIA, i, a, s)

      ! BA: is T_aa = 0?
      do b = rpa%nocc_ud(s) + 1, a
        ibs = rpa%iaTrans(i, b, s)
        ab = rpa%iaTrans(a, b, s) - svv(s)
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

      do j = i, rpa%nocc_ud(s)
        jas = rpa%iaTrans(j,a,s)

        ij = rpa%iaTrans(i, j, s) - soo(s)

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

    call assembleChunks(env, t)

    ! Terms for (P+-Q) of form (X+Y)^m_ib H^+_ab[(X+Y)^n]
    call getHplusXYfr(env, lr, rpa, transChrg, denseDesc, sym, species0, ovrXev, grndEigVecs,&
        & gammaMat, xpyn, vecHoo=vecHooXorY, vecHvv=vecHvvXorY)

    do ias = iGlobal, fGlobal
      call indxov(rpa%win, ias, rpa%getIA, i, a, s)
      do b = rpa%nocc_ud(s) + 1, a
        abs = rpa%iaTrans(a, b, s)
        ibs = rpa%iaTrans(i, b, s)
        ! For the forces, we have a factor of 2 here
        rhs(ias) = rhs(ias) - xpym(ibs) * vecHvvXorY(abs)
        ! Since vecHvvXpY has only upper triangle
        if (a /= b) then
          rhs(ibs) = rhs(ibs) - xpym(ias) * vecHvvXorY(abs)
        end if
      end do

      do j = i, rpa%nocc_ud(s)
        jas = rpa%iaTrans(j, a, s)
        ijs = rpa%iaTrans(i, j, s)
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
    call getHplusXYfr(env, lr, rpa, transChrg, denseDesc, sym, species0, ovrXev, grndEigVecs,&
        & gammaMat, xpym, vecHoo=vecHooXorY, vecHvv=vecHvvXorY)

    do ias = iGlobal, fGlobal
      call indxov(rpa%win, ias, rpa%getIA, i, a, s)
      do b = rpa%nocc_ud(s) + 1, a
        abs = rpa%iaTrans(a, b, s)
        ibs = rpa%iaTrans(i, b, s)
        ! For the forces, we have a factor of 2 here
        rhs(ias) = rhs(ias) - xpyn(ibs) * vecHvvXorY(abs)
        ! Since vecHvvXpY has only upper triangle
        if (a /= b) then
          rhs(ibs) = rhs(ibs) - xpyn(ias) * vecHvvXorY(abs)
        end if
      end do

      do j = i, rpa%nocc_ud(s)
        jas = rpa%iaTrans(j, a, s)
        ijs = rpa%iaTrans(i, j, s)
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
    allocate(vecHooT(sum(rpa%nxoo_ud)))
    allocate(vecHvvT(sum(rpa%nxvv_ud)))

    ! -RHS^+ += - H^+_ia[T^+]
    call getHplusMfr(env, lr, rpa, transChrg, denseDesc, species0, ovrXev, grndEigVecs, gammaMat,&
        & 3, t, vecHovT)

    rhs(iGlobal:fGlobal) = rhs(iGlobal:fGlobal) - vecHovT(iGlobal:fGlobal)

    ! Woo^+ += 0.5 * H^+_ij[T+Z] / Omega_mn, Z part computed later
    call getHplusMfr(env, lr, rpa, transChrg, denseDesc, species0, ovrXev, grndEigVecs, gammaMat,&
        & 1, t, vecHooT)

    do s = 1, nSpin
      call distributeRangeInChunks(env, 1, rpa%nxoo_ud(s), iGlobal, fGlobal)
      do ij = iGlobal, fGlobal
        ijs = ij + soo(s)
        woo(ij,s) = woo(ij,s) + vecHooT(ijs)
      end do
    end do

    ! Contributions due to range-separation
    if (rpa%tHybridXc) then
      allocate(vecHvvXpY(sum(rpa%nxvv_ud)))
      allocate(vecHvvXmY(sum(rpa%nxvv_ud)))
      allocate(vecHooXpY(sum(rpa%nxoo_ud)))
      allocate(vecHooXmY(sum(rpa%nxoo_ud)))

      ! Long-range part of H^+[(X+Y)^n] or H^-[(X-Y)^n] for occ-occ and vir-vir comp. of H
      call getHvvXY(env, orb, lr, rpa, transChrg, denseDesc, ovrXev, grndEigVecs, lrGamma,  1, xpyn,&
          & vecHvvXpY)

      call getHvvXY(env, orb, lr, rpa, transChrg, denseDesc, ovrXev, grndEigVecs, lrGamma, -1, xmyn,&
          & vecHvvXmY)

      call getHooXY(env, orb, lr, rpa, transChrg, denseDesc, ovrXev, grndEigVecs, lrGamma,  1, xpyn,&
          & vecHooXpY)

      call getHooXY(env, orb, lr, rpa, transChrg, denseDesc, ovrXev, grndEigVecs, lrGamma, -1, xmyn,&
          & vecHooXmY)

      call distributeRangeInChunks(env, 1, nxov, iGlobal, fGlobal)
      do ias = iGlobal, fGlobal

        call indXov(rpa%win, ias, rpa%getIA, i, a, s)
        do b = rpa%nocc_ud(s) + 1, nOrb
          ibs = rpa%iaTrans(i, b, s)
          abs = rpa%iaTrans(a, b, s)
          rhs(ias) = rhs(ias) - cExchange * 0.5_dp * xpym(ibs) * vecHvvXpY(abs)
          if (a >= b) then
            rhs(ias) = rhs(ias) - cExchange * 0.5_dp * xmym(ibs) * vecHvvXmY(abs)
          ! Only a>b is stored in vecHvvXmY, which is anti-symmetric
          else
            rhs(ias) = rhs(ias) + cExchange * 0.5_dp * xmym(ibs) * vecHvvXmY(abs)
          end if
        end do

        do j = 1, rpa%nocc_ud(s)
          jas = rpa%iaTrans(j, a, s)
          ijs = rpa%iaTrans(i, j, s)
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

      ! Long-range part of H^+[(X+Y)^n] or H^-[(X-Y)^n] for occ-occ and vir-vir comp. of H
      call getHvvXY(env, orb, lr, rpa, transChrg, denseDesc, ovrXev, grndEigVecs, lrGamma, 1, xpym,&
          & vecHvvXpY)

      call getHvvXY(env, orb, lr, rpa, transChrg, denseDesc, ovrXev, grndEigVecs, lrGamma, -1, xmym,&
          & vecHvvXmY)

      call getHooXY(env, orb, lr, rpa, transChrg, denseDesc, ovrXev, grndEigVecs, lrGamma,  1, xpym,&
          & vecHooXpY)

      call getHooXY(env, orb, lr, rpa, transChrg, denseDesc, ovrXev, grndEigVecs, lrGamma, -1, xmym,&
          & vecHooXmY)

      do ias = iGlobal, fGlobal

        call indXov(rpa%win, ias, rpa%getIA, i, a, s)
        do b = rpa%nocc_ud(s) + 1, nOrb
          ibs = rpa%iaTrans(i, b, s)
          abs = rpa%iaTrans(a, b, s)
          rhs(ias) = rhs(ias) - cExchange * 0.5_dp * xpyn(ibs) * vecHvvXpY(abs)
          if (a >= b) then
            rhs(ias) = rhs(ias) - cExchange * 0.5_dp * xmyn(ibs) * vecHvvXmY(abs)
          ! Only a>b is stored in vecHvvXmY, which is anti-symmetric
          else
            rhs(ias) = rhs(ias) + cExchange * 0.5_dp * xmyn(ibs) * vecHvvXmY(abs)
          end if
        end do

        do j = 1, rpa%nocc_ud(s)
          jas = rpa%iaTrans(j, a, s)
          ijs = rpa%iaTrans(i, j, s)
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
      call getHovT(env, orb, lr, rpa, transChrg, denseDesc, ovrXev, grndEigVecs, lrGamma, t, vecHovT)

      rhs(iGlobal:fGlobal) = rhs(iGlobal:fGlobal) - cExchange * vecHovT(iGlobal:fGlobal)

      ! Woo^+ += 0.5 * H^+_ij[T+Z] / Omega_mn, Z part computed later
      call getHooT(env, orb, lr, rpa, transChrg, denseDesc, ovrXev, grndEigVecs, lrGamma, t, vecHooT)

      do s = 1, nSpin
        call distributeRangeInChunks(env, 1, rpa%nxoo_ud(s), iGlobal, fGlobal)
        do ij = iGlobal, fGlobal
          ijs = ij + soo(s)
          woo(ij, s) = woo(ij, s) + cExchange * vecHooT(ijs)
        end do
      end do

    endif

    call assembleChunks(env, rhs)
    call assembleChunks(env, woo)
    call assembleChunks(env, wov)
    call assembleChunks(env, wvv)

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
  subroutine addNadiaGradients(env, orb, lr, rpa, transChrg, hybridXc, denseDesc, sym, species0,&
      & ovrXev, grndEigVecs, gammaMat, lrGamma, coord0, dq_ud, dqex, shift, xpyn, xmyn, xpym, xmym,&
      & woo, wov, wvv, skHamCont, skOverCont, derivator, rhoSqr, pc, nacv, deltaRho)

    !> Environment settings
    type(TEnvironment), intent(inout) :: env

    !> Data type for atomic orbital information
    type(TOrbitals), intent(in) :: orb

    !> Data structure for linear response
    type(TLinResp), intent(in) :: lr

    !> Run time parameters of the Casida routine
    type(TCasidaParameter), intent(in) :: rpa

    !> machinery for transition charges between single particle levels
    type(TTransCharges), intent(in) :: transChrg

    !> Data for range-separated calculation
    class(THybridXcFunc), allocatable, intent(inout) :: hybridXc

    !> Dense matrix descriptor
    type(TDenseDescr), intent(in) :: denseDesc

    !> Symmetry of the transition
    character, intent(in) :: sym

    !> Central cell chemical species
    integer, intent(in) :: species0(:)

    !> Overlap times ground state eigenvectors
    real(dp), intent(in) :: ovrXev(:,:,:)

    !> Ground state eigenvectors
    real(dp), intent(in) :: grndEigVecs(:,:,:)

    !> Softened coulomb matrix
    real(dp), intent(in) :: gammaMat(:,:)

    !> Electrostatic matrix, long-range corrected
    real(dp), allocatable, intent(in) :: lrGamma(:,:)

    !> Central cell atomic coordinates
    real(dp), intent(in) :: coord0(:,:)

    !> Ground state gross charges
    real(dp), intent(in) :: dq_ud(:,:)

    !> Charge differences from ground to excited state
    real(dp), intent(in) :: dqex(:,:)

    !> Ground state potentials (shift vector)
    real(dp), intent(in) :: shift(:)

    !> X+Y Furche term (state n)
    real(dp), intent(in) :: xpyn(:)

    !> X-Y Furche term (state n)
    real(dp), intent(in) :: xmyn(:)

    !> X+Y Furche term (state m)
    real(dp), intent(in) :: xpym(:)

    !> X-Y Furche term (state m)
    real(dp), intent(in) :: xmym(:)

    !> W vector occupied part
    real(dp), intent(in) :: woo(:,:)

    !> W vector occupied-virtual part
    real(dp), intent(in) :: wov(:)

    !> W vector virtual part
    real(dp), intent(in) :: wvv(:,:)

    !> H0 data
    type(TSlakoCont), intent(in) :: skHamCont

    !> Overlap data
    type(TSlakoCont), intent(in) :: skOverCont

    !> Differentiator for the non-scc matrices
    class(TNonSccDiff), intent(in) :: derivator

    !> Ground state density matrix
    real(dp), intent(in) :: rhoSqr(:,:,:)

    !> Transition density matrix
    real(dp), intent(in) :: pc(:,:,:)

    !> Resulting non-adiabatic coupling
    real(dp), intent(out) :: nacv(:,:)

    !> Difference density matrix (vs. uncharged atoms)
    real(dp), intent(inout), optional :: deltaRho(:,:,:)

    real(dp), allocatable :: shift_excited(:,:), xpyq(:,:), xpyqds(:,:)
    real(dp), allocatable :: shxpyq(:,:,:), xpycc(:,:,:,:), wcc(:,:,:), tmp5(:), tmp7(:,:), tmp11(:)
    real(dp), allocatable :: qTr(:), temp(:), dq(:), dm(:), dsigma(:)
    real(dp), allocatable :: dH0(:,:,:), dSo(:,:,:)
    real(dp), allocatable :: Dens(:,:), SpinDens(:,:)
    real(dp), allocatable :: xmycc(:,:,:,:), xpyas(:,:,:,:), xmyas(:,:,:,:)
    real(dp), allocatable :: overlap(:,:), lrGammaOrb(:,:), gammaLongRangePrime(:,:,:)
    real(dp), allocatable :: PS(:,:,:), DS(:,:,:), SPS(:,:,:), SDS(:,:,:), SX(:,:,:,:)
    real(dp), allocatable :: XS(:,:,:,:), SXS(:,:,:,:), SY(:,:,:,:), YS(:,:,:,:), SYS(:,:,:,:)
    real(dp), allocatable :: deltaRhoGlobal(:,:,:), grndEigVecsGlobal(:,:,:), DensGlobal(:,:,:)
    real(dp), allocatable :: xpy(:,:), xmy(:,:)
    real(dp) :: tmp1, tmp2, tmp3, tmp4, tmp6, tmp8, tmp9, tmp10, rab
    real(dp) :: diffvec(3), dgab(3), tmpVec(3), tmp3a, tmp3b, tmprs, tmprs2, tmps(2)
    integer, allocatable :: species(:)
    integer :: ia, i, j, a, b, ab, ij, m, n, mu, nu, xyz, iAt1, iAt2, ka
    integer :: indalpha, indalpha1, indbeta, indbeta1, soo(2), svv(2)
    integer :: iSp1, iSp2, iSpin, nSpin, iState, nOrb, iGlobal, fGlobal

    nSpin = size(grndEigVecs, dim=3)
    nOrb = orb%nOrb

    allocate(shift_excited(lr%nAtom, nSpin))
    allocate(xpyq(lr%nAtom, 2))
    allocate(shxpyq(lr%nAtom, nSpin, 2))
    allocate(xpycc(nOrb, nOrb, nSpin, 2))
    allocate(wcc(nOrb, nOrb, nSpin))
    allocate(qTr(lr%nAtom))
    allocate(temp(nOrb))
    allocate(tmp5(nSpin))
    allocate(tmp7(nSpin, 2))

    ! This should be changed to save memory
    allocate(xpy(rpa%nxov_rd, 2))
    allocate(xmy(rpa%nxov_rd, 2))
    xpy(:,1) = xpyn
    xpy(:,2) = xpym
    xmy(:,1) = xmyn
    xmy(:,2) = xmym

    allocate(Dens(nOrb, nOrb), DensGlobal(nOrb,nOrb,size(rhoSqr,dim=3)))
    do iSpin = 1, size(rhoSqr, dim=3)
      call distrib2replicated(env%blacs%orbitalGrid, denseDesc%blacsOrbSqr, &
                           &  rhoSqr(:,:,iSpin), DensGlobal(:,:,iSpin))
    enddo
    Dens(:,:) = sum(DensGlobal, dim=3)
    deallocate(DensGlobal)

    allocate(dH0(orb%mOrb, orb%mOrb, 3))
    allocate(dSo(orb%mOrb, orb%mOrb, 3))

    soo(:) = [0, rpa%nxoo_ud(1)]
    svv(:) = [0, rpa%nxvv_ud(1)]

    allocate(dq(lr%nAtom))
    dq(:) = dq_ud(:,1)

    if (lr%tSpin) then
      allocate(dm(lr%nAtom))
      allocate(xpyqds(lr%nAtom, 2))
      xpyqds = 0.0_dp
      allocate(tmp11(nSpin))

      ! FIXME: here nOrb is the global value but rhoSqr has dimension of nOrb local
      ! TODO: nOrb is the global but rhoSqr has the dimension of nOrb local
      !      The test NH forces does not run even for single process.
      allocate(SpinDens(nOrb, nOrb))
      SpinDens(:,:) = rhoSqr(:,:,1) - rhoSqr(:,:,2)

      allocate(dsigma(2))
      dsigma(1) = 1.0_dp
      dsigma(2) = -1.0_dp
      dm(:) = dq_ud(:,2)
    end if

    if (rpa%tHybridXc) then
      allocate(xmycc(nOrb, nOrb, nSpin, 2))
      allocate(xpyas(nOrb, nOrb, nSpin, 2))
      allocate(xmyas(nOrb, nOrb, nSpin, 2))
      allocate(PS(nOrb, nOrb, nSpin))
      allocate(DS(nOrb, nOrb, nSpin))
      allocate(SPS(nOrb, nOrb, nSpin))
      allocate(SDS(nOrb, nOrb, nSpin))
      allocate(SX(nOrb, nOrb, nSpin, 2))
      allocate(XS(nOrb, nOrb, nSpin, 2))
      allocate(SXS(nOrb, nOrb, nSpin, 2))
      allocate(SY(nOrb, nOrb, nSpin, 2))
      allocate(YS(nOrb, nOrb, nSpin, 2))
      allocate(SYS(nOrb, nOrb, nSpin, 2))
      allocate(overlap(nOrb, nOrb))
      allocate(lrGammaOrb(nOrb, nOrb))
      allocate(gammaLongRangePrime(3, lr%nAtom, lr%nAtom))

      ! Convert local arrays to global
      allocate(deltaRhoGlobal(norb,norb,size(deltaRho,dim=3)))
      allocate(grndEigVecsGlobal(norb,norb,size(grndEigVecs,dim=3)))
      do iSpin = 1, nSpin
        call distrib2replicated(env%blacs%orbitalGrid, denseDesc%blacsOrbSqr, &
                             &  deltaRho(:,:,iSpin), deltaRhoGlobal(:,:,iSpin))
        call distrib2replicated(env%blacs%orbitalGrid, denseDesc%blacsOrbSqr, &
                             &  grndEigVecs(:,:,iSpin), grndEigVecsGlobal(:,:,iSpin))
      enddo

      ! Symmetrize deltaRhoGlobal
      do mu = 1, size(deltaRhoGlobal, dim=1)
        do nu = mu + 1, size(deltaRhoGlobal, dim=2)
          deltaRhoGlobal(mu,nu,:) = deltaRhoGlobal(nu,mu,:)
        end do
      end do

      ! Compute long-range gamma derivative
      call distributeRangeInChunks(env, 1, lr%nAtom, iGlobal, fGlobal)
      gammaLongRangePrime(:,:,:) = 0.0_dp
      do iAt1 = iGlobal, fGlobal
        do iAt2 = 1, lr%nAtom
          if (iAt1 /= iAt2) then
            call getDirectionalCamGammaPrimeValue(hybridXc, tmpVec, iAt1, iAt2)
            gammaLongRangePrime(:, iAt1, iAt2) = tmpVec
          end if
        end do
      end do
      call assembleChunks(env,gammaLongRangePrime)

      ! Symmetrize S (can't we get S from caller?)
      call getSqrS(coord0, lr%nAtom, skOverCont, orb, denseDesc%iAtomStart, species0, overlap)
      call getSqrGamma(lr%nAtom, lrGamma, denseDesc%iAtomStart, lrGammaOrb)

    end if

    nacv(:,:) = 0.0_dp

    ! Excited state potentials at atomic sites
    do iSpin = 1, nSpin
      call hemv(shift_excited(:,iSpin), gammaMat, dqex(:,iSpin))
    end do

    ! xypq(alpha) = sum_ia (X+Y)_ia q^ia(alpha)
    ! Complexity nOrb * nOrb * nOrb
    xpyq = 0.0_dp
    do iState = 1, 2
      call transChrg%qMatVec(env, denseDesc, ovrXev, grndEigVecs, rpa%getIA, rpa%win,&
           & xpy(:,iState), xpyq(:,iState))
      ! complexity nOrb * nOrb
      shxpyq(:,:,iState) = 0.0_dp
      if (.not. lr%tSpin) then
        if (sym == "S") then
          call hemv(shxpyq(:,1,iState), gammaMat, xpyq(:,iState))
        else
          shxpyq(:,1,iState) = xpyq(:,iState) * lr%spinW(species0)
        end if
      else
        call transChrg%qMatVecDs(env, denseDesc, ovrXev, grndEigVecs, rpa%getIA, rpa%win,&
             & xpy(:,iState), xpyqds(:,iState))
        do iSpin = 1, nSpin
          call hemv(shxpyq(:,iSpin,iState), gammaMat, xpyq(:,iState))
          shxpyq(:,iSpin,iState) = shxpyq(:,iSpin,iState) + dsigma(iSpin)&
               & * lr%spinW(species0) * xpyqds(:,iState)
          shxpyq(:,iSpin,iState) = 0.5_dp * shxpyq(:,iSpin,iState)
        end do
      end if

      ! Calculate xpycc
      ! (xpycc)_{mu nu} = sum_{ia} (X + Y)_{ia} (grndEigVecs(mu,i)grndEigVecs(nu,a)
      ! + grndEigVecs(nu,i)grndEigVecs(mu,a))
      ! Complexity nOrb * nOrb * nOrb
      !
      ! xpycc(mu,nu) = sum_ia (X+Y)_ia grndEigVecs(mu,i) grndEigVecs(nu,a)
      ! xpycc(mu, nu) += sum_ia (X+Y)_ia grndEigVecs(mu,a) grndEigVecs(nu,i)
      call distributeRangeInChunks(env, 1, rpa%nxov_rd, iGlobal, fGlobal)
      xpycc(:,:,:,iState) = 0.0_dp
      do ia = iGlobal, fGlobal
        call indxov(rpa%win, ia, rpa%getIA, i, a, iSpin)
        ! Should replace with DSYR2 call:
        do nu = 1, nOrb
          do mu = 1, nOrb
            xpycc(mu,nu,iSpin,iState) = xpycc(mu,nu,iSpin,iState) + xpy(ia,iState)&
                & * (grndEigVecsGlobal(mu,i,iSpin)*grndEigVecsGlobal(nu,a,iSpin)&
                & + grndEigVecsGlobal(mu,a,iSpin)*grndEigVecsGlobal(nu,i,iSpin))
          end do
        end do
      end do
    end do
    call assembleChunks(env, xpycc)

    if (rpa%tHybridXc) then

      xmycc(:,:,:,:) = 0.0_dp
      xpyas(:,:,:,:) = 0.0_dp
      xmyas(:,:,:,:) = 0.0_dp

      call distributeRangeInChunks(env, 1, rpa%nxov_rd, iGlobal, fGlobal)
      do iState = 1,2
        ! Asymmetric contribution: xmycc_as = sum_ias (X-Y)_ias c_mas c_nis
        do ia = iGlobal, fGlobal
          call indxov(rpa%win, ia, rpa%getIA, i, a, iSpin)
          ! Should replace with DSYR2 call:
          do nu = 1, nOrb
            do mu = 1, nOrb
               xmycc(mu,nu,iSpin,iState) = xmycc(mu,nu,iSpin,iState) + xmy(ia,iState) *&
                & ( grndEigVecsGlobal(mu,i,iSpin) * grndEigVecsGlobal(nu,a,iSpin)&
                & + grndEigVecsGlobal(mu,a,iSpin) * grndEigVecsGlobal(nu,i,iSpin) )
               xpyas(mu,nu,iSpin,iState) = xpyas(mu,nu,iSpin,iState) + xpy(ia,iState) *&
                & grndEigVecsGlobal(mu,i,iSpin) * grndEigVecsGlobal(nu,a,iSpin)
               xmyas(mu,nu,iSpin,iState) = xmyas(mu,nu,iSpin,iState) + xmy(ia,iState) *&
                & grndEigVecsGlobal(mu,i,iSpin) * grndEigVecsGlobal(nu,a,iSpin)
            end do
          end do
        end do
        call assembleChunks(env, xmycc(:,:,:,iState))
        call assembleChunks(env, xpyas(:,:,:,iState))
        call assembleChunks(env, xmyas(:,:,:,iState))

        ! Account for normalization of S/T versus spin-polarized X+/-Y
        ! We have (X+Y)^S = 1/sqrt(2) [(X+Y)_up + (X+Y)_dn]
        if (lr%tSpin) then
          xmycc(:,:,:,:) = xmycc / sqrt(2._dp)
          xpyas(:,:,:,:) = xpyas / sqrt(2._dp)
          xmyas(:,:,:,:) = xmyas / sqrt(2._dp)
        end if

        do iSpin = 1, nSpin
          call symm(PS(:,:,iSpin), 'R', overlap, pc(:,:,iSpin), 'U', 1.0_dp, 0.0_dp, nOrb, nOrb)
          call symm(SPS(:,:,iSpin), 'L', overlap, PS(:,:,iSpin), 'U', 1.0_dp, 0.0_dp, nOrb, nOrb)
          call symm(DS(:,:,iSpin), 'R', overlap, deltaRhoGlobal(:,:,iSpin), 'U', 1.0_dp, 0.0_dp, nOrb,&
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
    ! Complexity nOrb * nOrb * nOrb

    ! Calculate the occ-occ part
    wcc(:,:,:) = 0.0_dp
    do iSpin = 1, nSpin
      call distributeRangeInChunks(env, 1, rpa%nxoo_ud(iSpin), iGlobal, fGlobal)
      do ij = iGlobal, fGlobal
        i = rpa%getIJ(ij + soo(iSpin), 1)
        j = rpa%getIJ(ij + soo(iSpin), 2)
        ! Replace with DSYR2 call:
        do mu = 1, nOrb
          do nu = 1, nOrb
            wcc(mu,nu,iSpin) = wcc(mu,nu,iSpin) + woo(ij,iSpin) *&
                & ( grndEigVecsGlobal(mu,i,iSpin)*grndEigVecsGlobal(nu,j,iSpin)&
                & + grndEigVecsGlobal(mu,j,iSpin)*grndEigVecsGlobal(nu,i,iSpin) )
          end do
        end do
      end do
    end do

    ! Calculate the occ-virt part: the same way as for xpycc
    call distributeRangeInChunks(env, 1, rpa%nxov_rd, iGlobal, fGlobal)
    do ia = iGlobal, fGlobal
      call indxov(rpa%win, ia, rpa%getIA, i, a, iSpin)
      ! Again replace with DSYR2 call:
      do nu = 1, nOrb
        do mu = 1, nOrb
          wcc(mu,nu,iSpin) = wcc(mu,nu,iSpin) + wov(ia) *&
              & ( grndEigVecsGlobal(mu,i,iSpin)*grndEigVecsGlobal(nu,a,iSpin)&
              & + grndEigVecsGlobal(mu,a,iSpin)*grndEigVecsGlobal(nu,i,iSpin) )
        end do
      end do
    end do

    ! Calculate the virt - virt part
    do iSpin = 1, nSpin
      call distributeRangeInChunks(env, 1, rpa%nxvv_ud(iSpin), iGlobal, fGlobal)
      do ab = iGlobal, fGlobal
        a = rpa%getAB(ab + svv(iSpin), 1)
        b = rpa%getAB(ab + svv(iSpin), 2)
        ! Replace with DSYR2 call:
        do mu = 1, nOrb
          do nu = 1, nOrb
            wcc(mu,nu,iSpin) = wcc(mu,nu,iSpin) + wvv(ab,iSpin) *&
                & ( grndEigVecsGlobal(mu,a,iSpin)*grndEigVecsGlobal(nu,b,iSpin)&
                & + grndEigVecsGlobal(mu,b,iSpin)*grndEigVecsGlobal(nu,a,iSpin) )
          end do
        end do
      end do
    end do
    call assembleChunks(env, wcc)

    ! Now calculating the force complexity : nOrb * nOrb * 3

    ! As have already performed nOrb**3 operation to get here, calculate for all atoms

    ! BA: only for non-periodic systems!
    do iAt1 = 1, lr%nAtom
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
        dgab(:) = diffvec * (-1.0_dp / rab**2 - expGammaPrime(rab, lr%HubbardU(iSp1),&
            & lr%HubbardU(iSp2)))

        tmp3a = 0.0_dp
        do iSpin = 1, nSpin
          tmp3a = tmp3a + dq(iAt1) * dqex(iAt2,iSpin) + dqex(iAt1,iSpin) * dq(iAt2)
        end do

        ! This term is symmetrized for NACV
        if (.not. lr%tSpin) then
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

        if (lr%tSpin) then
          tmp9 = lr%spinW(iSp1) * dm(iAt1) + lr%spinW(iSp2) * dm(iAt2)
          tmp11(:) = lr%spinW(iSp1) * dqex(iAt1,:) + lr%spinW(iSp2) * dqex(iAt2,:)
        end if

        if (rpa%tHybridXc) then
          tmprs = 0.0_dp
          tmps(:) = 0.0_dp
          do iSpin = 1, nSpin
            do mu = indAlpha, indAlpha1
              do nu = indBeta, indBeta1
                tmprs = tmprs +&
                    & ( 2.0_dp * (PS(mu,nu,iSpin) * DS(nu,mu,iSpin)&
                    & + PS(nu,mu,iSpin) * DS(mu,nu,iSpin))&
                    & + SPS(mu,nu,iSpin) * deltaRhoGlobal(mu,nu,iSpin)&
                    & + SPS(nu,mu,iSpin) * deltaRhoGlobal(nu,mu,iSpin)&
                    & + pc(mu,nu,iSpin) * SDS(mu,nu,iSpin)&
                    & + pc(nu,mu,iSpin) * SDS(nu,mu,iSpin) )

                tmprs = tmprs +&
                    & ( xpyas(mu,nu,iSpin,1) * SXS(mu,nu,iSpin,2)&
                    & + xpyas(nu,mu,iSpin,2) * SXS(nu,mu,iSpin,1)&
                    & + SX(mu,nu,iSpin,1) * XS(mu,nu,iSpin,2)&
                    & + SX(nu,mu,iSpin,2) * XS(nu,mu,iSpin,1) )
                tmprs = tmprs +&
                    & ( xpyas(mu,nu,iSpin,2) * SXS(mu,nu,iSpin,1)&
                    & + xpyas(nu,mu,iSpin,1) * SXS(nu,mu,iSpin,2)&
                    & + SX(mu,nu,iSpin,2) * XS(mu,nu,iSpin,1)&
                    & + SX(nu,mu,iSpin,1) * XS(nu,mu,iSpin,2) )

                tmprs = tmprs + 0.5_dp *&
                    & ( XS(mu,nu,iSpin,1) * XS(nu,mu,iSpin,2)&
                    & + XS(nu,mu,iSpin,2) * XS(mu,nu,iSpin,1)&
                    & + SXS(mu,nu,iSpin,1) * xpyas(nu,mu,iSpin,2)&
                    & + SXS(nu,mu,iSpin,2) * xpyas(mu,nu,iSpin,1)&
                    & + xpyas(mu,nu,iSpin,1) * SXS(nu,mu,iSpin,2)&
                    & + xpyas(nu,mu,iSpin,2) * SXS(mu,nu,iSpin,1)&
                    & + SX(mu,nu,iSpin,1) * SX(nu,mu,iSpin,2)&
                    & + SX(nu,mu,iSpin,2) * SX(mu,nu,iSpin,1) )
                tmprs = tmprs + 0.5_dp *&
                    & ( XS(mu,nu,iSpin,2) * XS(nu,mu,iSpin,1)&
                    & + XS(nu,mu,iSpin,1) * XS(mu,nu,iSpin,2)&
                    & + SXS(mu,nu,iSpin,2) * xpyas(nu,mu,iSpin,1)&
                    & + SXS(nu,mu,iSpin,1) * xpyas(mu,nu,iSpin,2)&
                    & + xpyas(mu,nu,iSpin,2) * SXS(nu,mu,iSpin,1)&
                    & + xpyas(nu,mu,iSpin,1) * SXS(mu,nu,iSpin,2)&
                    & + SX(mu,nu,iSpin,2) * SX(nu,mu,iSpin,1)&
                    & + SX(nu,mu,iSpin,1) * SX(mu,nu,iSpin,2) )

                tmprs = tmprs +&
                    & ( xmyas(mu,nu,iSpin,1) * SYS(mu,nu,iSpin,2)&
                    & + xmyas(nu,mu,iSpin,2) * SYS(nu,mu,iSpin,1)&
                    & + SY(mu,nu,iSpin,1) * YS(mu,nu,iSpin,2)&
                    & + SY(nu,mu,iSpin,2) * YS(nu,mu,iSpin,1) )
                tmprs = tmprs +&
                    & ( xmyas(mu,nu,iSpin,2) * SYS(mu,nu,iSpin,1)&
                    & + xmyas(nu,mu,iSpin,1) * SYS(nu,mu,iSpin,2)&
                    & + SY(mu,nu,iSpin,2) * YS(mu,nu,iSpin,1)&
                    & + SY(nu,mu,iSpin,1) * YS(nu,mu,iSpin,2) )

                tmprs = tmprs - 0.5_dp *&
                    & ( YS(mu,nu,iSpin,1) * YS(nu,mu,iSpin,2)&
                    & + YS(nu,mu,iSpin,2) * YS(mu,nu,iSpin,1)&
                    & + SYS(mu,nu,iSpin,1) * xmyas(nu,mu,iSpin,2)&
                    & + SYS(nu,mu,iSpin,2) * xmyas(mu,nu,iSpin,1)&
                    & + xmyas(mu,nu,iSpin,1) * SYS(nu,mu,iSpin,2)&
                    & + xmyas(nu,mu,iSpin,2) * SYS(mu,nu,iSpin,1)&
                    & + SY(mu,nu,iSpin,1) * SY(nu,mu,iSpin,2)&
                    & + SY(nu,mu,iSpin,2) * SY(mu,nu,iSpin,1) )
                tmprs = tmprs - 0.5_dp *&
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

                if (lr%tSpin) then
                  tmp8 = tmp8 + tmp9 * dSo(n,m,xyz) * dsigma(iSpin) * pc(mu,nu,iSpin)
                  tmp10 = tmp10 + tmp11(iSpin) * dSo(n,m,xyz) * dsigma(iSpin) * SpinDens(mu,nu)
                end if

                if (rpa%tHybridXc) then
                  tmprs = 0.0_dp
                  do ka = 1, nOrb
                    tmprs = tmprs + (lrGammaOrb(mu,ka) + lrGammaOrb(nu,ka)) *&
                        & ( PS(mu,ka,iSpin) * deltaRhoGlobal(nu,ka,iSpin)&
                        & + PS(nu,ka,iSpin) * deltaRhoGlobal(mu,ka,iSpin)&
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
  subroutine getHplusXYfr(env, lr, rpa, transChrg, denseDesc, sym, species0, ovrXev, grndEigVecs,&
      & frGamma, XorY, vecHoo, vecHvv)

    !> Environment settings
    type(TEnvironment), intent(inout) :: env

    !> Data structure for linear response
    type(TLinResp), intent(in) :: lr

    !> Run time parameters of the Casida routine
    type(TCasidaParameter), intent(in) :: rpa

    !> machinery for transition charges between single particle levels
    type(TTransCharges), intent(in) :: transChrg

    !> Dense matrix descriptor
    type(TDenseDescr), intent(in) :: denseDesc

    !> Symmetry of the transition
    character, intent(in) :: sym

    !> Central cell chemical species
    integer, intent(in) :: species0(:)

    !> overlap times eigenvector. (nOrb, nOrb) [distributed]
    real(dp), intent(in) :: ovrXev(:,:,:)

    !> eigenvectors (nOrb, nOrb) [distributed]
    real(dp), intent(in) :: grndEigVecs(:,:,:)

    !> full-range Gamma
    real(dp), intent(in) :: frGamma(:,:)

    !> RPA eigenvectors, either (X+Y) or (X-Y)
    real(dp), intent(in) :: XorY(:)

    !> Output vector H[V] occupied-occupied
    real(dp), optional, intent(out) :: vecHoo(:)

    !> Output vector H[V] virtual-virtual
    real(dp), optional, intent(out) :: vecHvv(:)

    integer :: nSpin, ab, s, abs, svv(2), ij, ijs, soo(2), iGlobal, fGlobal
    real(dp) :: fact
    real(dp), allocatable  :: xpyq(:), gamxpyq(:), qTr(:), xpyqds(:), gamxpyqds(:)

    nSpin = size(grndEigVecs, dim=3)

    allocate(xpyq(lr%nAtom))
    allocate(qTr(lr%nAtom))
    allocate(gamxpyq(lr%nAtom))
    soo(:) = [0, rpa%nxoo_ud(1)]
    svv(:) = [0, rpa%nxvv_ud(1)]

    if (present(vecHoo)) then
      vecHoo(:) = 0.0_dp
    end if
    if (present(vecHvv)) then
      vecHvv(:) = 0.0_dp
    end if

    xpyq(:) = 0.0_dp
    call transChrg%qMatVec(env, denseDesc, ovrXev, grndEigVecs, rpa%getIA, rpa%win, XorY, xpyq)

    if (.not. lr%tSpin) then  ! ---- spin-unpolarized case ----
      ! vecHvv(ab) = sum_jc K_ab,jc (X+Y)_jc
      if (sym == "S") then
        call hemv(gamxpyq, frGamma, xpyq)
        if (present(vecHvv)) then
          call distributeRangeInChunks(env, 1, rpa%nxvv_ud(1), iGlobal, fGlobal)
          do ab = iGlobal, fGlobal
            qTr(:) = transChrg%qTransAB(ab, env, denseDesc, ovrXev, grndEigVecs, rpa%getAB)
            vecHvv(ab) = 2.0_dp * sum(qTr * gamxpyq)
          end do
        end if
        if (present(vecHoo)) then
          call distributeRangeInChunks(env, 1, rpa%nxoo_ud(1), iGlobal, fGlobal)
          do ij = iGlobal, fGlobal
            qTr(:) = transChrg%qTransIJ(ij, env, denseDesc, ovrXev, grndEigVecs, rpa%getIJ)
            ! vecHoo(ij) = sum_kb K_ij,kb (X+Y)_kb
            vecHoo(ij) = 2.0_dp * sum(qTr * gamxpyq)
          end do
        end if
      else ! Triplet case
        if (present(vecHvv)) then
          call distributeRangeInChunks(env, 1, rpa%nxvv_ud(1), iGlobal, fGlobal)
          do ab = iGlobal, fGlobal
            qTr(:) = transChrg%qTransAB(ab, env, denseDesc, ovrXev, grndEigVecs, rpa%getAB)
            vecHvv(ab) = 2.0_dp * sum(qTr * xpyq * lr%spinW(species0))
          end do
        end if
        if (present(vecHoo)) then
          call distributeRangeInChunks(env, 1, rpa%nxoo_ud(1), iGlobal, fGlobal)
          do ij = iGlobal, fGlobal
            qTr(:) = transChrg%qTransIJ(ij, env, denseDesc, ovrXev, grndEigVecs, rpa%getIJ)
            vecHoo(ij) = 2.0_dp * sum(qTr * xpyq * lr%spinW(species0))
          end do
        end if
      end if

    else  ! ---- spin-polarized case -----

      allocate(xpyqds(lr%nAtom))
      allocate(gamxpyqds(lr%nAtom))
      xpyqds(:) = 0.0_dp
      call transChrg%qMatVecDs(env, denseDesc, ovrXev, grndEigVecs, rpa%getIA, rpa%win, XorY,&
          & xpyqds)

      call hemv(gamxpyq, frGamma,  xpyq)
      do s = 1, 2
        if (s == 1) then
          fact = 1.0_dp
        else
          fact = -1.0_dp
        end if
        if (present(vecHvv)) then
          call distributeRangeInChunks(env, 1, rpa%nxvv_ud(s), iGlobal, fGlobal)
          do ab = iGlobal, fGlobal
            abs = ab + svv(s)
            qTr(:) = transChrg%qTransAB(abs, env, denseDesc, ovrXev, grndEigVecs, rpa%getAB)
            vecHvv(abs) = sum(qTr * gamxpyq)
            ! Magnetization part
            vecHvv(abs) = vecHvv(abs) + fact * sum(qTr * xpyqds * lr%spinW(species0))
          end do
        end if
        if (present(vecHoo)) then
          call distributeRangeInChunks(env, 1, rpa%nxoo_ud(s), iGlobal, fGlobal)
          do ij = iGlobal, fGlobal
            ijs = ij + soo(s)
            qTr(:) = transChrg%qTransIJ(ijs, env, denseDesc, ovrXev, grndEigVecs, rpa%getIJ)
            vecHoo(ijs) = sum(qTr * gamxpyq)
            !magnetization part
            vecHoo(ijs) = vecHoo(ijs) + fact * sum(qTr * xpyqds * lr%spinW(species0))
          end do
        end if
      end do

    end if

    if (present(vecHoo)) then
      call assembleChunks(env, vecHoo)
    end if
    if (present(vecHvv)) then
      call assembleChunks(env, vecHvv)
    end if

  end subroutine getHplusXYfr


  !> Computes full range part of H^+_pq [M] as defined in Furche JCP 117 7433 (2002) eq. 20
  !! Here pq are arbitrary orbitals and M is a general matrix with ov,oo,vv components
  !! iMode = 1: returns oo components of H
  !! iMode = 2: returns vv components of H
  !! iMode = 3: returns ov components of H
  !! Note: The full range part of H^- is zero!
  !! Routine currently does not work for M_ia /= 0 on entry!
  subroutine getHplusMfr(env, lr, rpa, transChrg, denseDesc, species0, ovrXev, grndEigVecs,&
      & frGamma, iMode, matM, vecH)

    !> Environment settings
    type(TEnvironment), intent(inout) :: env

    !> Data structure for linear response
    type(TLinResp), intent(in) :: lr

    !> Run time parameters of the Casida routine
    type(TCasidaParameter), intent(in) :: rpa

    !> machinery for transition charges between single particle levels
    type(TTransCharges), intent(in) :: transChrg

    !> Dense matrix descriptor
    type(TDenseDescr), intent(in) :: denseDesc

    !> Central cell chemical species
    integer, intent(in) :: species0(:)

    !> overlap times eigenvector. (nOrb, nOrb) [distributed]
    real(dp), intent(in) :: ovrXev(:,:,:)

    !> eigenvectors (nOrb, nOrb) [distributed]
    real(dp), intent(in) :: grndEigVecs(:,:,:)

    !> full-range Gamma
    real(dp), intent(in) :: frGamma(:,:)

    !> Type of return vector (oo, vv, ov)
    integer , intent(in)  :: iMode

    !> Input Matrix spin-resolved
    real(dp), intent(in) :: matM(:,:,:)

    !> Output vector H[M]
    real(dp), intent(out) :: vecH(:)

    integer :: nSpin, ab, i, j, a, b, s, svv(2), ij, soo(2), ias, iGlobal, fGlobal
    real(dp), dimension(2) :: spinFactor = [1.0_dp, -1.0_dp]
    real(dp), allocatable  :: xpyq(:), gamxpyq(:), qTr(:), xpyqds(:), gamxpyqds(:)
    real(dp), allocatable  :: gamqt(:)

    if (iMode == 1) then
      @:ASSERT(size(vecH) == sum(rpa%nxoo_ud))
    else if (iMode == 2) then
      @:ASSERT(size(vecH) == sum(rpa%nxvv_ud))
    else
      @:ASSERT(size(vecH) == rpa%nxov_rd)
    end if

    nSpin = size(grndEigVecs, dim=3)
    if (lr%tSpin) then
      allocate(xpyqds(lr%nAtom))
      allocate(gamxpyqds(lr%nAtom))
    endif
    allocate(xpyq(lr%nAtom))
    allocate(qTr(lr%nAtom))
    allocate(gamxpyq(lr%nAtom))
    allocate(gamqt(lr%nAtom))
    soo(:) = [0, rpa%nxoo_ud(1)]
    svv(:) = [0, rpa%nxvv_ud(1)]

    vecH(:) = 0.0_dp
    ! gamxpyq(iAt2) = sum_ij q_ij(iAt2) M_ij
    gamxpyq(:) = 0.0_dp
    if (lr%tSpin) then
      gamxpyqds(:) = 0.0_dp
    end if

    do s = 1, nSpin
      call distributeRangeInChunks(env, 1, rpa%nxoo_ud(s), iGlobal, fGlobal)
      do ij = iGlobal, fGlobal
        i = rpa%getIJ(ij + soo(s), 1)
        j = rpa%getIJ(ij + soo(s), 2)
        qTr(:) = transChrg%qTransIJ(ij + soo(s), env, denseDesc, ovrXev, grndEigVecs, rpa%getIJ)
        if (i == j) then
          gamxpyq(:) = gamxpyq(:) + matM(i,j,s) * qTr(:)
          if (lr%tSpin) then
            gamxpyqds(:) = gamxpyqds(:) + matM(i,j,s) * qTr(:) * spinFactor(s)
          end if
        else
          gamxpyq(:) = gamxpyq(:) + (matM(i,j,s) + matM(j,i,s)) * qTr(:)
          if (lr%tSpin) then
            gamxpyqds(:) = gamxpyqds(:) + (matM(i,j,s) + matM(j,i,s)) * qTr(:) * spinFactor(s)
          end if
        end if
      end do

      call distributeRangeInChunks(env, 1, rpa%nxvv_ud(s), iGlobal, fGlobal)
      ! gamxpyq(iAt2) += sum_ab q_ab(iAt2) M_ab
      do ab = iGlobal, fGlobal
        a = rpa%getAB(ab + svv(s), 1)
        b = rpa%getAB(ab + svv(s), 2)
        qTr(:) = transChrg%qTransAB(ab + svv(s), env, denseDesc, ovrXev, grndEigVecs, rpa%getAB)
        if (a == b) then
          gamxpyq(:) = gamxpyq(:) + matM(a,b,s) * qTr(:)
          if (lr%tSpin) then
            gamxpyqds(:) = gamxpyqds(:) + matM(a,b,s) * qTr(:) * spinFactor(s)
          end if
        else
          ! Factor 2 because of symmetry of the matrix
          gamxpyq(:) = gamxpyq(:) + (matM(a,b,s) + matM(b,a,s)) * qTr(:)
          if (lr%tSpin) then
            gamxpyqds(:) = gamxpyqds(:) + (matM(a,b,s) + matM(b,a,s)) * qTr(:) * spinFactor(s)
          end if
        end if
      end do

    end do

    call distributeRangeInChunks(env, 1, rpa%nxov_rd, iGlobal, fGlobal)
    ! gamxpyq(iAt2) += sum_ab q_ab(iAt2) M_ia
    do ias = iGlobal, fGlobal
      i = rpa%getIA(ias, 1)
      a = rpa%getIA(ias, 2)
      s = rpa%getIA(ias, 3)
      qTr(:) = transChrg%qTransIA(ias, env, denseDesc, ovrXev, grndEigVecs, rpa%getIA, rpa%win)
      gamxpyq(:) = gamxpyq(:) + (matM(i,a,s) + matM(a,i,s)) * qTr(:)
      if (lr%tSpin) then
         gamxpyqds(:) = gamxpyqds(:) + (matM(i,a,s) + matM(a,i,s)) * qTr(:) * spinFactor(s)
      end if
    end do

    call assembleChunks(env, gamxpyq)
    if (lr%tSpin) then
      call assembleChunks(env, gamxpyqds)
    end if

    ! gamqt(iAt1) = sum_iAt2 gamma_iAt1,iAt2 gamxpyq(iAt2)
    call hemv(gamqt, frGamma, gamxpyq)

    if (iMode == 1) then
      do s = 1, nSpin
        call distributeRangeInChunks(env, 1, rpa%nxoo_ud(s), iGlobal, fGlobal)
        do ij = iGlobal, fGlobal
          i = rpa%getIJ(ij + soo(s), 1)
          j = rpa%getIJ(ij + soo(s), 2)
          qTr(:) = transChrg%qTransIJ(ij + soo(s), env, denseDesc, ovrXev, grndEigVecs, rpa%getIJ)
          if (.not. lr%tSpin) then
            vecH(ij + soo(s)) = 4.0_dp * dot_product(gamqt,qTr)
          else
            vecH(ij + soo(s)) = 2.0_dp * dot_product(gamqt,qTr)
            vecH(ij + soo(s)) = vecH(ij + soo(s)) +&
                & spinFactor(s) * dot_product(gamxpyqds * lr%spinW(species0), qTr)
          end if
        end do
      end do
    else if (iMode == 2) then
      do s = 1, nSpin
        call distributeRangeInChunks(env, 1, rpa%nxvv_ud(s), iGlobal, fGlobal)
        do ab = iGlobal, fGlobal
          a = rpa%getAB(ab + svv(s), 1)
          b = rpa%getAB(ab + svv(s), 2)
          qTr(:) = transChrg%qTransAB(ab + svv(s), env, denseDesc, ovrXev, grndEigVecs, rpa%getAB)
          if (.not. lr%tSpin) then
            vecH(ab + svv(s)) = 4.0_dp * dot_product(gamqt,qTr)
          else
            vecH(ab + svv(s)) = 2.0_dp * dot_product(gamqt,qTr)
            vecH(ab + svv(s)) = vecH(ab + svv(s)) +&
                & spinFactor(s) * dot_product(gamxpyqds * lr%spinW(species0), qTr)
          end if
        end do
      end do
    else
      call distributeRangeInChunks(env, 1, rpa%nxov_rd, iGlobal, fGlobal)
      do ias = iGlobal, fGlobal
        i = rpa%getIA(ias, 1)
        a = rpa%getIA(ias, 2)
        s = rpa%getIA(ias, 3)
        qTr(:) = transChrg%qTransIA(ias, env, denseDesc, ovrXev, grndEigVecs, rpa%getIA, rpa%win)
        if (.not. lr%tSpin) then
           vecH(ias) = 4.0_dp * dot_product(gamqt,qTr)
        else
          vecH(ias) =  2.0_dp * dot_product(gamqt,qTr)
          vecH(ias) = vecH(ias) +&
              & spinFactor(s) * dot_product(gamxpyqds * lr%spinW(species0), qTr)
        end if
      end do
    end if
    call assembleChunks(env, vecH)

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
      X1(:) = reshape(excDerivs(:,:,nexcGrad), [3 * nAtoms])
    else
      X1(:) = reshape(excDerivs(:,:,nexcGrad) - excDerivs(:,:,1), [3 * nAtoms])
    end if
    ! Last entry of array holds coupling between states indNACouplings(1/2)
    X2(:) = reshape(naCouplings(:,:,nCoupl), [3 * nAtoms])

    ! Original Bearbark suggestion would be dS1/dq, Harabuchi chooses 0.5 (dS0/dq + dS1/dq)
    dE2(:) = reshape(derivs + 0.5_dp * (excDerivs(:,:,nexcGrad) + excDerivs(:,:,1)), [3 * nAtoms])

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


  !> Computes transition charges for a given excited state according to
  !! JCP 140, 174101 (2014). Works for singlet states and spin polarized calculations.
  !! Initialization prevents entering this routine for triplets for which the charges
  !! would simply be zero.
  subroutine getAndWriteTransitionCharges(env, lr, rpa, transChrg, sym, denseDesc, ovrXev, &
      & grndEigVecs, xpy, fdTransQ, fdTagged, taggedWriter)

    !> Environment settings
    type(TEnvironment), intent(inout) :: env

    !> Data structure for linear response
    type(TLinResp), intent(in) :: lr

    !> Run time parameters of the Casida routine
    type(TCasidaParameter), intent(in) :: rpa

    !> machinery for transition charges between single particle levels
    type(TTransCharges), intent(in) :: transChrg

    !> symmetry flag (singlet or triplet)
    character, intent(in) :: sym

    !> Dense matrix descriptor
    type(TDenseDescr), intent(in) :: denseDesc

    !> overlap times eigenvector. (nOrb, nOrb) [distributed]
    real(dp), intent(in) :: ovrXev(:,:,:)

    !> eigenvectors (nOrb, nOrb) [distributed]
    real(dp), intent(in) :: grndEigVecs(:,:,:)

    !> (X+Y) RPA eigenvectors (global)
    real(dp), intent(in) :: xpy(:,:)

    !> File descriptor for ATQ.DAT
    type(TFileDescr) :: fdTransQ

    !> File unit for tagged output (> -1 for write out)
    type(TFileDescr), intent(in) :: fdTagged

    !> Tagged writer
    type(TTaggedWriter), intent(inout) :: taggedWriter

    integer :: i, ia
    real(dp), allocatable :: atomicTransQ(:), qia(:)
    real(dp) :: preFactor

    @:ASSERT(size(xpy(:, lr%nstat)) == rpa%nxov_rd)

    if (lr%tSpin) then
      preFactor = 1.0_dp
    else
      preFactor = sqrt(2.0_dp)
    endif

    allocate(atomicTransQ(lr%nAtom))
    allocate(qia(lr%nAtom))

    atomicTransQ(:) = 0.0_dp
    do ia = 1, rpa%nxov_rd
      qia(:) = transChrg%qTransIA(ia, env, denseDesc, ovrXev, grndEigVecs, rpa%getIA, rpa%win)
      atomicTransQ(:) = atomicTransQ + preFactor * qia * xpy(ia,lr%nstat)
    end do

    call openfile(fdTransQ, transChrgOut, mode="w")


    write(fdTransQ%unit, '(a)') "#"
    write(fdTransQ%unit, '(a,2x,a,5x,a)') "#", "atom", "transition charge"
    write(fdTransQ%unit, '(a)') "#"

    do i = 1, lr%nAtom
      write(fdTransQ%unit, '(2x,i5,5x,f12.9)') i, atomicTransQ(i)
    end do

    call closeFile(fdTransQ)

    if (fdTagged%isConnected()) then
      call taggedWriter%write(fdTagged%unit, tagLabels%transQ, atomicTransQ**2)
    end if

  end subroutine getAndWriteTransitionCharges

end module dftbp_timedep_linrespgrad
