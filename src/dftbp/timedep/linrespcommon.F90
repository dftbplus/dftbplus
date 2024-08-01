!--------------------------------------------------------------------------------------------------!
!  DFTB+: general package for performing fast atomistic simulations                                !
!  Copyright (C) 2006 - 2023  DFTB+ developers group                                               !
!                                                                                                  !
!  See the LICENSE file for terms of usage and distribution.                                       !
!--------------------------------------------------------------------------------------------------!

#:include 'common.fypp'

!> Helper routines for the linear response modules.
module dftbp_timedep_linrespcommon
  use dftbp_common_accuracy, only : dp, elecTolMax
  use dftbp_common_constants, only: Hartree__eV, au__Debye, cExchange
  use dftbp_common_schedule, only : assembleChunks, gatherChunks, distributeRangeInChunks
  use dftbp_common_file, only : TFileDescr, openFile, closeFile
  use dftbp_dftb_onsitecorrection, only : getOnsME
  use dftbp_io_message, only : error
  use dftbp_math_blasroutines, only : hemv, gemm
  use dftbp_math_eigensolver, only : heev
  use dftbp_math_sorting, only : index_heap_sort
  use dftbp_timedep_linresptypes, only : TLinResp, TCasidaParameter
  use dftbp_timedep_transcharges, only : TTransCharges, transq
  use dftbp_type_commontypes, only : TOrbitals
  use dftbp_type_densedescr, only : TDenseDescr
  use dftbp_common_environment, only : TEnvironment

#:if WITH_SCALAPACK
  use dftbp_math_scalafxext, only : distrib2replicated
#:endif
  
  implicit none

  public

  !> Prefactor of 2/3.
  real(dp), parameter :: twothird = 2.0_dp / 3.0_dp

  !> Names of output files
  character(*), parameter :: excitedQOut = "XCH.DAT"
  character(*), parameter :: excitedDipoleOut = "XREST.DAT"
  character(*), parameter :: singlePartOut = "SPX.DAT"

contains


  !> Find (possibly degenerate) transitions with stronger dipole transition strengths than a
  !! tolerance, count them and place at the start of the appropriate arrays.
  subroutine dipSelect(wij, sposz, win, transd, nxov_r, threshold, grndEigVal, getIA)

    !> Energies of transitions
    real(dp), intent(inout) :: wij(:)

    !> Transition strengths
    real(dp), intent(inout) :: sposz(:)

    !> Index array for transitions to single particle transitions
    integer, intent(inout) :: win(:)

    !> Transition dipole matrix elements for each atom
    real(dp), intent(inout) :: transd(:,:)

    !> Number of excitations in space
    integer, intent(out) :: nxov_r

    !> Threshold for dipole cutoff
    real(dp), intent(in) :: threshold

    !> Ground state eigenvalues
    real(dp), intent(in) :: grndEigVal(:,:)

    !> Index array of transitions
    integer, intent(in) :: getIA(:,:)

    integer :: nxov, ii, jj, iOcc, iVrt, iStart, iSpin
    real(dp) :: eOcc, eExc, mu

    nxov = size(wij)
    @:ASSERT(nxov == size(sposz))
    @:ASSERT(nxov == size(win))
    @:ASSERT(all(shape(transd) == [nxov, 3]))
    @:ASSERT(size(getIA, dim=1) >= nxov)
    @:ASSERT(size(getIA, dim=2) == 3)
    @:ASSERT(threshold >= 0.0_dp)

    call indxov(win, 1, getIA, iOcc, iVrt, iSpin)
    eOcc = grndEigVal(iOcc, 1)
    eExc = wij(1)
    iStart = 1
    nxov_r = 0
    ! Check, to single precision tolerance, for degeneracies when selecting bright transitions
    do ii = 2, nxov
      call indxov(win, ii, getIA, iOcc, iVrt, iSpin)
      ! Check if this is a still within a degenerate group, otherwise process the group
      if (abs(grndEigVal(iOcc, 1)-eOcc) > epsilon(0.0) .or. abs(wij(ii) - eExc) > epsilon(0.0)) then
        eOcc = grndEigVal(iOcc, 1)
        eExc = wij(ii)
        mu = 0.0_dp
        ! Loop over transitions in the group and check the oscillator strength
        do jj = iStart, ii - 1
          call indxov(win, jj, getIA, iOcc, iVrt, iSpin)
          mu = mu + sposz(jj)
        end do
        ! If something in the group is bright, so include them all
        if (mu > threshold) then
          do jj = iStart, ii - 1
            nxov_r = nxov_r + 1
            win(nxov_r) = win(jj)
            wij(nxov_r) = wij(jj)
            sposz(nxov_r) = sposz(jj)
            transd(nxov_r,:) = transd(jj,:)
          end do
        end if
        iStart = ii
      end if
    end do

    ! Last group in the transitions
    mu = 0.0_dp
    do jj = iStart, nxov
      call indxov(win, jj, getIA, iOcc, iVrt, iSpin)
      mu = mu + sposz(jj)
    end do
    if (mu > threshold) then
      do jj = iStart, nxov
        nxov_r = nxov_r + 1
        win(nxov_r) = win(jj)
        wij(nxov_r) = wij(jj)
        sposz(nxov_r) = sposz(jj)
        transd(nxov_r,:) = transd(jj,:)
      end do
    end if

  end subroutine dipSelect


  !> Computes individual indices from the compound occ-virt excitation index.
  pure subroutine indxov(win, indx, getIA, ii, jj, ss)

    !> Index array after sorting of eigenvalues
    integer, intent(in) :: win(:)

    !> Compound excitation index
    integer, intent(in) :: indx

    !> Index array of transitions after sorting of eigenvalues
    integer, intent(in) :: getIA(:,:)

    !> Initial (filled) state
    integer, intent(out) :: ii

    !> Final (empty) state
    integer, intent(out) :: jj

    !> Spin channel: 1 (up) or 2 (down)
    integer, intent(out) :: ss

    integer :: indo

    indo = win(indx)
    ii = getIA(indo, 1)
    jj = getIA(indo, 2)
    ss = getIA(indo, 3)

  end subroutine indxov


  !> Computes individual indices from the compound occ-occ excitation index.
  subroutine indxoo(indx, ii, jj)

    !> Compound excitation index
    integer, intent(in) :: indx

    !> Initial state
    integer, intent(out) :: ii

    !> Final state
    integer, intent(out) :: jj

    call indxvv(0, indx, ii, jj)

  end subroutine indxoo


  !> Computes individual indices from the compound virtual-virtual excitation index.
  subroutine indxvv(nocc, indx, ii, jj)

    !> Compund excitation index
    integer, intent(in) :: nocc

    !> Number of occupied states
    integer, intent(in) :: indx

    !> Initial state
    integer, intent(out) :: ii

    !> Final state
    integer, intent(out) :: jj

    real(dp) :: re

    @:ASSERT(indx > 0)

    ! Solve a quadratic for the row of a matrix between occupied and virtual states, given the
    ! number of the element in the matrix
    re = sqrt(real(indx, dp) * 2.0_dp - 1.75_dp) - 0.5_dp

    ! actual row
    ii  = floor(re) + 1 + nocc
    ! column
    jj  = indx - ((ii - 1 - nocc) * (ii - nocc)) / 2 + nocc

    @:ASSERT(ii > 0)
    @:ASSERT(jj > 0)

  end subroutine indxvv


  !> Builds array to convert from orbital pairs to a compound index, reverse of index generated by
  !! getSPExcitations (which were then potentially windowed afterwards).
  subroutine rindxov_array(win, nxov, nxoo, nxvv, getIA, getIJ, getAB, iatrans)

    !> Array for indexing excitations
    integer, intent(in) :: win(:)

    !> Number of transitions from occupied to virtual
    integer, intent(in) :: nxov

    !> Number of transitions from occupied to occupied
    integer, intent(in) :: nxoo

    !> Number of transitions from virtual to virtual
    integer, intent(in) :: nxvv

    !> Array of the occupied->virtual pairs
    integer, intent(in) :: getIA(:,:)

    !> Array of the occupied->occupied pairs
    integer, intent(in) :: getIJ(:,:)

    !> Array of the virtual->virtual pairs
    integer, intent(in) :: getAB(:,:)

    !> Resulting index array from orbital pairs to compound index
    integer, intent(out) :: iatrans(:,:,:)

    integer :: ia, ij, ab

    @:ASSERT(size(getIA, dim=2) == 3)
    @:ASSERT(size(getIJ, dim=2) == 3)
    @:ASSERT(size(getAB, dim=2) == 3)

    ! Store reverse indices

    iaTrans(:,:,:) = 0

    ! The transition charges q_rs are symmetrical in rs
    ! We only need one index per pair of orbitals
    !$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(ia) SCHEDULE(RUNTIME)
    do ia = 1, nxov
      iatrans(getIA(win(ia), 1), getIA(win(ia), 2), getIA(win(ia), 3)) = ia
      iatrans(getIA(win(ia), 2), getIA(win(ia), 1), getIA(win(ia), 3)) = ia
    end do
    !$OMP END PARALLEL DO

    !$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(ij) SCHEDULE(RUNTIME)
    do ij = 1, nxoo
      iatrans(getIJ(ij, 1), getIJ(ij, 2), getIJ(ij, 3)) = ij
      iatrans(getIJ(ij, 2), getIJ(ij, 1), getIJ(ij, 3)) = ij
    end do
    !$OMP END PARALLEL DO

    !$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(ab) SCHEDULE(RUNTIME)
    do ab = 1, nxvv
      iatrans(getAB(ab, 1), getAB(ab, 2), getAB(ab, 3)) = ab
      iatrans(getAB(ab, 2), getAB(ab, 1), getAB(ab, 3)) = ab
    end do
    !$OMP END PARALLEL DO

  end subroutine rindxov_array


  !> Calculate the transition block at a specific atom.
  subroutine transDens(ii, jj, iAt, iAtomStart, nOrb, updwn, ovrXev, grndEigVecs, qq_ij)

    !> Index of initial state
    integer, intent(in) :: ii

    !> Index of final state
    integer, intent(in) :: jj

    !> Atom number
    integer, intent(in) :: iAt

    !> Starting position of each atom in the list of orbitals
    integer, intent(in) :: iAtomStart(:)

    !> Up spin channel (T) or down spin channel (F)
    logical, intent(in) :: updwn

    !> Overlap times eigenvector: sum_m Smn cmi (nOrb, nOrb)
    real(dp), intent(in) :: ovrXev(:,:,:)

    !> Eigenvectors (nOrb, nOrb)
    real(dp), intent(in) :: grndEigVecs(:,:,:)

    !> Transition charge block
    real(dp), intent(out) :: qq_ij(:,:)

    integer :: nOrb, iOrb1, iOrb2
    integer :: mu, nu, ss

    ss = 1
    if (.not. updwn) ss = 2

    !mu = iAtomStart(iAt)
    qq_ij(:,:) = 0.0_dp

    do iOrb1 = 1, nOrb
      do iOrb2 = iOrb1, nOrb
        mu = iAtomStart(iAt) + iOrb1 - 1
        nu = iAtomStart(iAt) + iOrb2 - 1
        qq_ij(iOrb1, iOrb2) = 0.25_dp&
            & * (grndEigVecs(mu, ii, ss) * ovrXev(nu, jj, ss)&
            & + grndEigVecs(mu, jj, ss) * ovrXev(nu, ii, ss)&
            & + grndEigVecs(nu, ii, ss) * ovrXev(mu, jj, ss)&
            & + grndEigVecs(nu, jj, ss) * ovrXev(mu, ii, ss))
        if (iOrb1 /= iOrb2) then
          qq_ij(iOrb2, iOrb1) = qq_ij(iOrb1, iOrb2)
        end if
      end do
    end do

  end subroutine transDens


  !> Returns the (spatial) MO overlap between orbitals in different spin channels.
  function MOoverlap(pp, qq, ovrXev, grndEigVecs) result(S_pq)

    !> Orbital in alpha channel
    integer, intent(in) :: pp

    !> Orbital in  beta channel
    integer, intent(in) :: qq

    !> Overlap times ground single particle state wavefunctions
    real(dp), intent(in) :: ovrXev(:,:,:)

    !> Ground state single particle wavefunctions
    real(dp), intent(in) :: grndEigVecs(:,:,:)

    !> Resulting overlap between states
    real(dp):: S_pq

    @:ASSERT(all(shape(grndEigVecs) == shape(ovrXev)))
    @:ASSERT(size(grndEigVecs, dim=3) == 2)

    S_pq = sum(grndEigVecs(:, pp, 1) * ovrXev(:, qq, 2))

  end function MOoverlap


  !> Square root of differences in occupations between filled and empty states.
  !! We need occupations per spin channel ([0:1]) also for closed shell systems.
  subroutine getSqrOcc(occNr, win, nmatup, nmat, getIA, tSpin, n_ij)

    !> Occupation of states
    real(dp), intent(in) :: occNr(:,:)

    !> Index array for transitions to single particle transitions
    integer, intent(in) :: win(:)

    !> Number of up-up excitations
    integer, intent(in) :: nmatup

    !> Number of transitions to scale
    integer, intent(in) :: nmat

    !> Array of the occupied->virtual pairs
    integer, intent(in) :: getIA(:,:)

    !> Is system spin polarized?
    logical, intent(in) :: tSpin

    !> Resulting scaled matrix
    real(dp), intent(out) :: n_ij(:)

    integer :: ia, ii, jj, ss
    logical :: updwn
    real(dp) :: docc_ij

    !$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(ia, ii, jj, updwn, docc_ij) SCHEDULE(RUNTIME)
    do ia = 1, nmat
      call indxov(win, ia, getIA, ii, jj, ss)
      updwn = (win(ia) <= nmatup)
      if (updwn) then
        docc_ij = occNr(ii, 1) - occNr(jj, 1)
      else
        docc_ij = occNr(ii, 2) - occNr(jj, 2)
      end if
      if (.not. tSpin) then
         docc_ij = docc_ij / 2.0_dp
      end if
      n_ij(ia) = sqrt(docc_ij)
    end do
    !$OMP END PARALLEL DO

  end subroutine getSqrOcc

  !> Multiplies the excitation supermatrix with a supervector.
  !! For the hermitian RPA eigenvalue problem this corresponds to \Omega_ias,jbt * v_jbt
  !! (spin polarized case) or \Omega^{S/T}_ia,jb * v_jb (singlet/triplet)
  !!
  !! For the standard RPA, (A+B)_ias,jbt * v_jbt needs to be computed (similar for singlet/triplet)
  !! (see definitions in Marc Casida, in Recent Advances in Density Functional Methods,
  !!  World Scientific, 1995, Part I, p. 155.)
  !! Note: we actually compute sqrt(n_is-n_as) (A+B)_ias,jbt sqrt(n_jt-n_bt), with the
  !! occupations n_is, correct for finite T.
  !! See also Dominguez JCTC 9 4901 (2013), Kranz JCTC 13 1737 (2017) for DFTB specifics.
  !!
  !! Note: In order not to store the entire supermatrix (nmat, nmat), the various pieces are
  !! assembled individually and multiplied directly with the corresponding part of the supervector.
  !! Note MPI: The supervector is distributed: iGlobal and fGlobal mark the relevant indices 
  !! in global arrays
  subroutine actionAplusB(iGlobal, fGlobal, env, orb, lr, rpa, transChrg, sym, denseDesc, species0,&
    &  ovrXev, grndEigVecs, frGamma, tAplusB, vin, vout, lrGamma)

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

    !> machinery for transition charges between single particle levels
    type(TTransCharges), intent(in) :: transChrg
    
    !> symmetry flag (singlet or triplet)
    character, intent(in) :: sym

    !> Dense matrix descriptor
    type(TDenseDescr), intent(in) :: denseDesc

    !> chemical species of the atoms
    integer, intent(in) :: species0(:)

    !> overlap times eigenvector. (nOrb, nOrb) [distributed]
    real(dp), intent(in) :: ovrXev(:,:,:)
    
    !> eigenvectors (nOrb, nOrb) [distributed]
    real(dp), intent(in) :: grndEigVecs(:,:,:)

    !> DFTB gamma matrix (nAtm, nAtom)
    real(dp), intent(in) :: frGamma(:,:)

    !> long-range Gamma if in use
    real(dp), allocatable, optional, intent(in) :: lrGamma(:,:)

    !> whether (A+B)v or (Omega)v should be returned
    logical, intent(in) :: tAplusB

    !> vector to multiply with size(nMat)
    real(dp), intent(in) :: vin(:)

    !> vector containing the result on exit size(nMat)
    real(dp), intent(out) :: vout(:)

    integer :: nMat, nOrb, ia, nxvv_a
    integer :: ii, aa, ss, jj, bb, ias, jbs, abs, ijs, jas, ibs
    integer :: nLoc, myia, myab, myja, nLocAB, iGlobalAB, fGlobalAB
    integer, allocatable ::  getABasym(:,:)

    ! somewhat ugly, but fast small arrays on stack:
    real(dp) :: otmp(size(frGamma, dim=1)), gtmp(size(frGamma, dim=1))
    real(dp) :: qij(size(frGamma, dim=1)) ! qij Working array (used for excitation charges). (nAtom)
    real(dp) :: vout_ons(size(vin))
    real(dp), allocatable :: vLoc(:), vGlb(:), vGlb2(:)
    real(dp), parameter :: spinFactor(2) = [1.0_dp, -1.0_dp]
    ! for later use to change HFX contribution
    real(dp), allocatable :: qv(:,:)

    @:ASSERT(size(vin) == size(vout))

    ! The dimension is guessed from input vector
    nMat = size(vin)
    @:ASSERT(nMat <= rpa%nxov_rd)
    nOrb = orb%nOrb

    ! Local chunk of RPA vectors have this size under MPI
    nLoc = fGlobal - iGlobal + 1
    @:ASSERT(nMat == nLoc)
    allocate(vLoc, mold=vin)
    vLoc(:) = 0.0_dp
    vout(:) = 0.0_dp

    if(tAplusB) then
      vLoc(:) = vin
    else
      vLoc(:) = vin * sqrt(rpa%wij(iGlobal:fGlobal))
    endif

    ! Compute w_ia = q^ia_A G_AB q^jb v_jb 

    oTmp(:) = 0.0_dp
    do ia = iGlobal, fGlobal
      myia = ia - iGlobal + 1
      qij(:) = transChrg%qTransIA(ia, env, denseDesc, ovrXev, grndEigVecs, rpa%getIA, rpa%win)
      otmp(:) = otmp(:) + qij(:) * vLoc(myia) * rpa%sqrOccIA(ia)
    enddo
    
    call assembleChunks(env, otmp)

    if (.not. lr%tSpin) then !-----------spin-unpolarized systems--------------

      if (sym == 'S') then

        call hemv(gtmp, frGamma, otmp)

        ! 4 * wn * (g * Q)
        vOut(:) = 0.0_dp
        do ia = iGlobal, fGlobal
          myia = ia - iGlobal + 1
          qij(:) = transChrg%qTransIA(ia, env, denseDesc, ovrXev, grndEigVecs, rpa%getIA, rpa%win)
          vOut(myia) = 4.0_dp * rpa%sqrOccIA(ia) * dot_product(qij, gTmp)
        enddo  

      else
       
        otmp(:) = otmp * lr%spinW(species0)

        ! 4 * wn * (o * Q)
        vOut(:) = 0.0_dp
        call transChrg%qVecMat(env, denseDesc, ovrXev, grndEigVecs, rpa%getIA, rpa%win, oTmp, vOut,&
            & iGlobal-1)
        vOut(:) = 4.0_dp * rpa%sqrOccIA(iGlobal:fGlobal) * vOut


      end if

    else !--------------spin-polarized systems--------

      call hemv(gtmp, frGamma, otmp)

      otmp(:) = 0.0_dp

      do ia = iGlobal, fGlobal
        myia = ia - iGlobal + 1
        
        qij(:) = transChrg%qTransIA(ia, env, denseDesc, ovrXev, grndEigVecs, rpa%getIA, rpa%win)

        ! singlet gamma part (S)
        vOut(myia) = 2.0_dp * rpa%sqrOccIA(ia) * dot_product(qij, gtmp)

        ! magnetization part (T1)
        ss = rpa%getIA(rpa%win(ia), 3)

        otmp(:) = otmp + spinFactor(ss) * rpa%sqrOccIA(ia) * vLoc(myia) * qij

      end do

      call assembleChunks(env, otmp)

      otmp(:) = otmp * lr%spinW(species0)

      do ia = iGlobal, fGlobal
        myia = ia - iGlobal + 1
        
        qij(:) = transChrg%qTransIA(ia, env, denseDesc, ovrXev, grndEigVecs, rpa%getIA, rpa%win)

        ss = rpa%getIA(rpa%win(ia), 3)
      
        vOut(myia) = vOut(myia) + 2.0_dp * rpa%sqrOccIA(ia) * spinFactor(ss) * &
            & dot_product(qij, otmp)

      end do

    end if
    
    ! The largest matrices in this part are of dim nAtoms*nVir*nVir, seems logical to distribute 
    ! them. The local dimensions are different from the case of nOcc*nVir, we require new offsets.
    ! A new index array is required because the standard abs index is symmetrical in a and b  
    if (rpa%tHybridXc) then
      
      ! Number of vir-vir transitions a->b _and_ b->a, summed over spin channels
      nxvv_a = sum(rpa%nvir_ud**2)
      allocate(getABasym(nxvv_a, 3))
  
      do ss = 1, size(rpa%iaTrans, dim=3)
        do aa = rpa%nocc_ud(ss) + 1, nOrb
          do bb = rpa%nocc_ud(ss) + 1, nOrb
            abs = (bb - rpa%nocc_ud(ss) - 1) * rpa%nvir_ud(ss) + aa - rpa%nocc_ud(ss)
            if (ss==2) then
              abs = abs + rpa%nvir_ud(1)**2
            end if
            getABasym(abs, 1) = aa
            getABasym(abs, 2) = bb
            getABasym(abs, 3) = ss
          end do
        end do
      end do

      call distributeRangeInChunks(env, 1, nxvv_a, iGlobalAB, fGlobalAB)
      nLocAB = fGlobalAB - iGlobalAB + 1

      allocate(qv(lr%nAtom, max(nLoc, nLocAB)))
      allocate(vGlb(rpa%nxov_rd))
      allocate(vGlb2(rpa%nxov_rd))

      ! TD-LC-DFTB seems to require two additional global arrays of dim nOcc*nVir, qv is local 
      qv(:,:) = 0.0_dp
      vGlb(:) = 0.0_dp
      vGlb2(:) = 0.0_dp

      ! Gather local arrays in corresponding global array
      call gatherChunks(env, 1, rpa%nxov_rd, vLoc, vGlb)

      ! Compute w_ia = q^ij_A GLR_AB q^ab v_jb
      do jas = iGlobal, fGlobal
        myja = jas - iGlobal + 1
        jj = rpa%getIA(rpa%win(jas), 1)
        aa = rpa%getIA(rpa%win(jas), 2)
        ss = rpa%getIA(rpa%win(jas), 3)
        do bb = rpa%nocc_ud(ss) + 1, nOrb
          abs = rpa%iaTrans(aa, bb, ss)
          qij(:) = transChrg%qTransAB(abs, env, denseDesc, ovrXev, grndEigVecs, rpa%getAB)
          jbs = rpa%iaTrans(jj, bb, ss)
          qv(:,myja) = qv(:,myja) + qij(:) * vGlb(jbs) * rpa%sqrOccIA(jbs)
        end do
        
        otmp(:) = qv(:,myja)
        call hemv(qv(:,myja), lrGamma, otmp, uplo='U')
        
        do ii = 1, rpa%nocc_ud(ss)
          ijs = rpa%iaTrans(ii, jj, ss)
          qij(:) = transChrg%qTransIJ(ijs, env, denseDesc, ovrXev, grndEigVecs, rpa%getIJ)
          ias = rpa%iaTrans(ii, aa, ss)
          vGlb2(ias) = vGlb2(ias) - cExchange * rpa%sqrOccIA(ias) * dot_product(qij, qv(:,myja))
        end do
      end do

      ! Compute w_ia = q^ib_A GLR_AB q^ja v_jb 
      qv(:,:) = 0.0_dp
      do abs = iGlobalAB, fGlobalAB
        myab = abs - iGlobalAB + 1     
        aa = getABasym(abs, 1)
        bb = getABasym(abs, 2)
        ss = getABasym(abs, 3)
        do jj = 1, rpa%nocc_ud(ss)
          jas =  rpa%iaTrans(jj, aa, ss)
          jbs =  rpa%iaTrans(jj, bb, ss)
          qij(:) = transChrg%qTransIA(jas, env, denseDesc, ovrXev, grndEigVecs, rpa%getIA, rpa%win)
          qv(:,myab) = qv(:,myab) + qij(:) * vGlb(jbs) * rpa%sqrOccIA(jbs)
        end do

        otmp(:) = qv(:,myab)
        call hemv(qv(:,myab), lrGamma, otmp, uplo='U')

        do ii = 1, rpa%nocc_ud(ss)
          ibs = rpa%iaTrans(ii, bb, ss)
          qij(:) = transChrg%qTransIA(ibs, env, denseDesc, ovrXev, grndEigVecs, rpa%getIA, rpa%win)
          ias = rpa%iaTrans(ii, aa, ss)
          vGlb2(ias) = vGlb2(ias) - cExchange * rpa%sqrOccIA(ias) * dot_product(qij, qv(:,myab)) 
        end do
      end do

      ! Get contribution of all ranks to global array
      call assembleChunks(env, vGlb2) 

      do jas = iGlobal, fGlobal
        myja = jas - iGlobal + 1        
        vout(myja) = vout(myja) + vGlb2(jas) 
      end do

    end if
   
    if (allocated(lr%OnsiteMatrixElements)) then
      call onsiteEner(env, orb, lr, rpa, denseDesc, sym, species0, ovrXev, grndEigVecs, vLoc,&
          & vout_ons, iGlobal-1)
      
      vout(:) = vout + vout_ons
    end if
    
    vout(:) = vout + rpa%wij(iGlobal:fGlobal) * vLoc
    
    if (.not. tAplusB) then
      vOut(:) = vOut * sqrt(rpa%wij(iGlobal:fGlobal))
    endif

  end subroutine actionAplusB

  !> Multiplies the excitation supermatrix with a supervector.
  !! (A-B)_ias,jbt * v_jbt is computed (and similar for singlet/triplet)
  !! (see definitions in Marc Casida, in Recent Advances in Density Functional Methods,
  !!  World Scientific, 1995, Part I, p. 155.)
  !! See also Dominguez JCTC 9 4901 (2013), Kranz JCTC 13 1737 (2017) for DFTB specifics.
  !! Note MPI: The supervector is distributed, locSize gives the local size.
  subroutine actionAminusB(iGlobal, fGlobal, env, orb, lr, rpa, transChrg, denseDesc,  ovrXev,&
      & grndEigVecs, vin, vout, lrGamma)

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

    !> machinery for transition charges between single particle levels
    type(TTransCharges), intent(in) :: transChrg

    !> Dense matrix descriptor
    type(TDenseDescr), intent(in) :: denseDesc

    !> overlap times eigenvector. (nOrb, nOrb) [distributed]
    real(dp), intent(in) :: ovrXev(:,:,:)
    
    !> eigenvectors (nOrb, nOrb) [distributed]
    real(dp), intent(in) :: grndEigVecs(:,:,:)

    !> vector to multiply with size(nMat)
    real(dp), intent(in) :: vin(:)

    !> vector containing the result on exit size(nMat)
    real(dp), intent(out) :: vout(:)

    !> long-range Gamma if in use
    real(dp), allocatable, optional, intent(in) :: lrGamma(:,:)

    integer :: nMat, nOrb, nxvv_a
    integer :: ii, aa, ss, jj, bb, ias, jbs, abs, ijs, jas, ibs
    integer :: nLoc, myab, myja, nLocAB, iGlobalAB, fGlobalAB
    integer, allocatable ::  getABasym(:,:)
    real(dp), allocatable :: vGlb(:), vGlb2(:)
    real(dp), allocatable :: qv(:,:)
    real(dp), allocatable :: qij(:), otmp(:)

    @:ASSERT(size(vin) == size(vout))

    nMat = size(vin)
    @:ASSERT(nMat <= rpa%nxov_rd)
    vout(:) = 0.0_dp

    nOrb = orb%nOrb

    ! Local chunk of RPA vectors have this size under MPI   
    nLoc = fGlobal - iGlobal + 1
    @:ASSERT(nMat == nLoc)

    if (rpa%tHybridXc) then

      allocate(qij(lr%nAtom))
      allocate(otmp(lr%nAtom))
      ! Number of vir-vir transitions a->b _and_ b->a, summed over spin channels
      nxvv_a = sum(rpa%nvir_ud**2)
      allocate(getABasym(nxvv_a, 3))
  
      do ss = 1, size(rpa%iaTrans, dim=3)
        do aa = rpa%nocc_ud(ss) + 1, nOrb
          do bb = rpa%nocc_ud(ss) + 1, nOrb
            abs = (bb - rpa%nocc_ud(ss) - 1) * rpa%nvir_ud(ss) + aa - rpa%nocc_ud(ss)
            if (ss==2) then
              abs = abs + rpa%nvir_ud(1)**2
            end if
            getABasym(abs, 1) = aa
            getABasym(abs, 2) = bb
            getABasym(abs, 3) = ss
          end do
        end do
      end do

      call distributeRangeInChunks(env, 1, nxvv_a, iGlobalAB, fGlobalAB)
      nLocAB = fGlobalAB - iGlobalAB + 1

      allocate(qv(lr%nAtom, max(nLoc, nLocAB)))
      allocate(vGlb(rpa%nxov_rd))
      allocate(vGlb2(rpa%nxov_rd))

      ! TD-LC-DFTB seems to require two additional global arrays of dim nOcc*nVir, qv is local 
      qv(:,:) = 0.0_dp
      vGlb(:) = 0.0_dp
      vGlb2(:) = 0.0_dp

      ! Gather local arrays in corresponding global array
      call gatherChunks(env, 1, rpa%nxov_rd, vin, vGlb)

      ! Compute w_ia = q^ij_A GLR_AB q^ab v_jb
      do jas = iGlobal, fGlobal
        myja = jas - iGlobal + 1      
        jj = rpa%getIA(rpa%win(jas), 1)
        aa = rpa%getIA(rpa%win(jas), 2)
        ss = rpa%getIA(rpa%win(jas), 3)
        do bb = rpa%nocc_ud(ss) + 1, nOrb
          abs = rpa%iaTrans(aa, bb, ss)
          qij(:) = transChrg%qTransAB(abs, env, denseDesc, ovrXev, grndEigVecs, rpa%getAB)
          jbs = rpa%iaTrans(jj, bb, ss)
          qv(:,myja) = qv(:,myja) + qij(:) * vGlb(jbs) * rpa%sqrOccIA(jbs)
        end do
        
        otmp(:) = qv(:,myja)
        call hemv(qv(:,myja), lrGamma, otmp, uplo='U')
        
        do ii = 1, rpa%nocc_ud(ss)
          ijs = rpa%iaTrans(ii, jj, ss)
          qij(:) = transChrg%qTransIJ(ijs, env, denseDesc, ovrXev, grndEigVecs, rpa%getIJ)
          ias = rpa%iaTrans(ii, aa, ss)
          vGlb2(ias) = vGlb2(ias) - cExchange * rpa%sqrOccIA(ias) * dot_product(qij, qv(:,myja))
        end do
      end do

      ! Compute w_ia = q^ib_A GLR_AB q^ja v_jb 
      qv(:,:) = 0.0_dp
      do abs = iGlobalAB, fGlobalAB
        myab = abs - iGlobalAB + 1        
        aa = getABasym(abs, 1)
        bb = getABasym(abs, 2)
        ss = getABasym(abs, 3)
        do jj = 1, rpa%nocc_ud(ss)
          jas =  rpa%iaTrans(jj, aa, ss)
          jbs =  rpa%iaTrans(jj, bb, ss)
          qij(:) = transChrg%qTransIA(jas, env, denseDesc, ovrXev, grndEigVecs, rpa%getIA, rpa%win)
          qv(:,myab) = qv(:,myab) + qij(:) * vGlb(jbs) * rpa%sqrOccIA(jbs)
        end do

        otmp(:) = qv(:,myab)
        call hemv(qv(:,myab), lrGamma, otmp, uplo='U')

        do ii = 1, rpa%nocc_ud(ss)
          ibs = rpa%iaTrans(ii, bb, ss)
          qij(:) = transChrg%qTransIA(ibs, env, denseDesc, ovrXev, grndEigVecs, rpa%getIA, rpa%win)
          ias = rpa%iaTrans(ii, aa, ss)
          vGlb2(ias) = vGlb2(ias) + cExchange * rpa%sqrOccIA(ias) * dot_product(qij, qv(:,myab)) 
        end do
      end do

      ! Get contribution of all ranks to global array
      call assembleChunks(env, vGlb2) 

      do jas = iGlobal, fGlobal
        myja = jas - iGlobal + 1   
        vout(myja) = vout(myja) +  vGlb2(jas) 
      end do
      
    endif
    
    ! orb. energy difference diagonal contribution
    vout(:) = vout + rpa%wij(iGlobal:fGlobal) * vin
    
  end subroutine actionAminusB

  !> Generates initial matrices M+ and M- for the RPA algorithm by Stratmann
  !! (JCP 109 8218 (1998).
  !! M+/- = (A+/-B)_ias,jbt (spin polarized) (A+/-B)^{S/T}_ia,jb (closed shell)
  !! Here ias,jbt <= initDim
  !! Also computed is v+/- = (A+/-B)_ias,jbt with ias <= nMat, jbt <= initDim
  !! Note: Routine not set up to handle onsite corrections.
  !! Note: Not yet OpenMP parallelized
  !! Note MPI: The supervector is distributed: iGlobal and fGlobal mark the relevant indices 
  !! in global arrays
  subroutine initialSubSpaceMatrixApmB(iGlobal, fGlobal, env, lr, rpa, transChrg, sym,&
      & denseDesc, species0, ovrXev, grndEigVecs, frGamma, lrGamma, initDim, vP, vM, mP, mM) 

    !> Starting index of current rank in global RPA vectors
    integer, intent(in) :: iGlobal

    !> End index of current rank in global RPA vectors
    integer, intent(in) :: fGlobal

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

    !> chemical species of the atoms
    integer, intent(in) :: species0(:)

    !> overlap times eigenvector. (nOrb, nOrb) [distributed]
    real(dp), intent(in) :: ovrXev(:,:,:)
    
    !> eigenvectors (nOrb, nOrb) [distributed]
    real(dp), intent(in) :: grndEigVecs(:,:,:)

    !> DFTB gamma matrix (nAtm, nAtom)
    real(dp), intent(in) :: frGamma(:,:)

    !> long-range Gamma if in use
    real(dp), allocatable, intent(in) :: lrGamma(:,:)

    !> initial RPA subspace
    integer, intent(in) :: initDim

    !> output vector v+ (nMat, initDim)
    real(dp), intent(out) :: vP(:,:)

    !> output vector v+ (nMat, initDim)
    real(dp), intent(out) :: vM(:,:)

    !> output matrix M+ (initDim, initDim)
    real(dp), intent(out) :: mP(:,:)

    !> output matrix M- (initDim, initDim)
    real(dp), intent(out) :: mM(:,:)

    integer :: nMat, ia, jb, ii, jj, ss, tt
    real(dp), allocatable :: oTmp(:), gTmp(:), qTr(:)
    real(dp), parameter :: spinFactor(2) = [1.0_dp, -1.0_dp]
    real(dp) :: rTmp
    integer :: aa, bb, iat, jbs, abs, ijs, ibs, jas
    integer :: nLoc, myia, myii

    nMat = size(vP, dim=1) ! also known as nXov
    
    ! Local chunk of RPA vectors have this size under MPI
    nLoc = fGlobal - iGlobal + 1
    @:ASSERT(nMat == nLoc)

    allocate(gTmp(lr%nAtom))
    allocate(oTmp(lr%nAtom))
    allocate(qTr(lr%nAtom))

    vP(:,:) = 0.0_dp
    vM(:,:) = 0.0_dp
    mP(:,:) = 0.0_dp
    mM(:,:) = 0.0_dp

    oTmp(:) = 0.0_dp

    !-----------spin-unpolarized systems--------------
    if (.not. lr%tSpin) then

      if (sym == 'S') then

        ! Full range coupling matrix contribution: 4 * sum_A q^ia_A sum_B gamma_AB q^jb_B
        do jb = 1, initDim
          qTr(:) = transChrg%qTransIA(jb, env, denseDesc, ovrXev, grndEigVecs, rpa%getIA, rpa%win)
          qTr(:) = qTr * rpa%sqrOccIA(jb)

          call hemv(oTmp, frGamma, qTr)

          do ia = iGlobal, fGlobal
            myia = ia - iGlobal + 1
            qTr(:) = transChrg%qTransIA(ia, env, denseDesc, ovrXev, grndEigVecs,&
                & rpa%getIA, rpa%win)
            vP(myia,jb) = 4.0_dp * rpa%sqrOccIA(ia) * dot_product(qTr, oTmp)
          end do

        end do

      else

        ! Full range triplets contribution: 2 * sum_A q^ia_A M_A q^jb_A
        do jb = 1, initDim
          qTr(:) = transChrg%qTransIA(jb, env, denseDesc, ovrXev, grndEigVecs, rpa%getIA, rpa%win)
          oTmp(:) = qTr * rpa%sqrOccIA(jb) * lr%spinW(species0)

          do ia = iGlobal, fGlobal
            myia = ia - iGlobal + 1
            qTr(:) = transChrg%qTransIA(ia, env, denseDesc, ovrXev, grndEigVecs, rpa%getIA,&
                & rpa%win)
            vP(myia,jb) = vP(myia,jb) + 4.0_dp * rpa%sqrOccIA(ia) * dot_product(qTr, oTmp)
          end do

        end do

      end if
    !--------------spin-polarized systems--------
    else

      do jb = 1, initDim
        qTr(:) = transChrg%qTransIA(jb, env, denseDesc, ovrXev, grndEigVecs, rpa%getIA, rpa%win)
        qTr(:) = qTr * rpa%sqrOccIA(jb)

        call hemv(gTmp, frGamma, qTr)

        oTmp(:) = qTr

        do ia = iGlobal, fGlobal
          myia = ia - iGlobal + 1
          qTr(:) = transChrg%qTransIA(ia, env, denseDesc, ovrXev, grndEigVecs, rpa%getIA, rpa%win)
          vP(myia,jb) = 2.0_dp * rpa%sqrOccIA(ia) * dot_product(qTr, gTmp)
        end do

        ss = rpa%getIA(rpa%win(jb), 3)

        oTmp(:)  =  spinFactor(ss) * 2.0_dp * lr%spinW(species0) * oTmp

        do ia = iGlobal, fGlobal
          myia = ia - iGlobal + 1
          qTr(:) = transChrg%qTransIA(ia, env, denseDesc, ovrXev, grndEigVecs, rpa%getIA, rpa%win)
          tt = rpa%getIA(rpa%win(ia), 3)
          vP(myia,jb) = vP(myia,jb) + spinFactor(tt) * rpa%sqrOccIA(ia) * dot_product(qTr, oTmp)
        end do

      end do

    end if

    if (rpa%tHybridXc) then

      do jbs = 1, initDim
        jj = rpa%getIA(rpa%win(jbs), 1)
        bb = rpa%getIA(rpa%win(jbs), 2)
        ss = rpa%getIA(rpa%win(jbs), 3)

        do iat = iGlobal, fGlobal
          myia = iat - iGlobal + 1
          ii = rpa%getIA(rpa%win(iat), 1)
          aa = rpa%getIA(rpa%win(iat), 2)
          tt = rpa%getIA(rpa%win(iat), 3)

          if (ss /= tt) cycle

          abs = rpa%iaTrans(aa, bb, ss)
          qTr(:) = transChrg%qTransAB(abs, env, denseDesc, ovrXev, grndEigVecs, rpa%getAB)
          oTmp(:) = 0.0_dp
          call hemv(oTmp, lrGamma, qTr)

          ijs = rpa%iaTrans(ii, jj, ss)
          qTr(:) = transChrg%qTransIJ(ijs, env, denseDesc, ovrXev, grndEigVecs, rpa%getIJ)
          rTmp = cExchange * rpa%sqrOccIA(iat) * rpa%sqrOccIA(jbs) * dot_product(qTr, oTmp)
          vP(myia,jbs) = vP(myia,jbs) - rTmp
          vM(myia,jbs) = vM(myia,jbs) - rTmp

          ibs = rpa%iaTrans(ii, bb, ss)
          qTr(:) = transChrg%qTransIA(ibs, env, denseDesc, ovrXev, grndEigVecs, rpa%getIA, rpa%win)
          oTmp(:) = 0.0_dp
          call hemv(oTmp, lrGamma, qTr)

          jas = rpa%iaTrans(jj, aa, ss)
          qTr(:) = transChrg%qTransIA(jas, env, denseDesc, ovrXev, grndEigVecs, rpa%getIA, rpa%win)
          rTmp = cExchange * rpa%sqrOccIA(iat) * rpa%sqrOccIA(jbs) * dot_product(qTr, oTmp)
          vP(myia,jbs) = vP(myia,jbs) - rTmp
          vM(myia,jbs) = vM(myia,jbs) + rTmp
        end do

      end do
    end if

    do ia = iGlobal, fGlobal
      myia = ia - iGlobal + 1
      if(ia > initDim) cycle
      vP(myia,ia) = vP(myia,ia) + rpa%wij(ia)
      vM(myia,ia) = vM(myia,ia) + rpa%wij(ia)
    end do

    do ii = iGlobal, fGlobal
      myii = ii - iGlobal + 1
      if(ii > initDim) exit
      do jj = ii, initDim
        mP(ii,jj) = vP(myii,jj)
        mP(jj,ii) = mP(ii,jj)
        mM(ii,jj) = vM(myii,jj)
        mM(jj,ii) = mM(ii,jj)
      end do
    end do
    
    call assembleChunks(env, mP)
    call assembleChunks(env, mM)
    
  end subroutine initialSubSpaceMatrixApmB  

  !> Onsite energy corrections
  !! Routine also works for MPI if optional index offset is provided. The arrays ovrXev and
  !! grndEigVecs still need to be global!
  subroutine onsiteEner(env, orb, lr, rpa, denseDesc, sym, species0, ovrXev, eigVec, vin, vout,&
      & indexOffSet) 

    !> Environment settings
    type(TEnvironment), intent(inout) :: env

    !> Data type for atomic orbital information
    type(TOrbitals), intent(in) :: orb
     
    !> Data structure for linear response
    type(TLinResp), intent(in) :: lr

    !> Run time parameters of the Casida routine
    type(TCasidaParameter), intent(in) :: rpa

    !> Dense matrix descriptor
    type(TDenseDescr), intent(in) :: denseDesc

    !> Symmetry flag (singlet or triplet)
    character, intent(in) :: sym

    !> Chemical species of the atoms
    integer, intent(in) :: species0(:)   

    !> Overlap times ground state wavefunctions
    real(dp), intent(in) :: ovrXev(:,:,:)

    !> Ground state wavefunctions
    real(dp), intent(in) :: eigVec(:,:,:)

    !> Vector to multiply with size(nmat)
    real(dp), intent(in) :: vin(:)

    !> Vector containing the result on exit size(nmat)
    real(dp), intent(out) :: vout(:)

    !> Offset of vector index (optional) which determines orbital pair
    integer, intent(in), optional :: indexOffSet

    real(dp), allocatable :: eigVecGlb(:,:,:), ovrXevGlb(:,:,:)
    real(dp) :: otmp(orb%mOrb,orb%mOrb,size(species0),2)
    real(dp) :: fact
    real(dp) :: onsite(orb%mOrb,orb%mOrb,2)
    real(dp) :: qq_ij(orb%mOrb,orb%mOrb)
    logical :: updwn
    integer :: nmat, nOrb, iOff
    integer :: ia, iAGlb, iAt, iSp, iSh, iOrb, ii, jj, ss, sindx(2), iSpin, nSpin
    real(dp) :: degeneracy, partTrace

    if (present(indexOffSet)) then
      iOff = indexOffSet
    else
      iOff = 0
    end if

    nmat = size(vin)
    nSpin = size(ovrXev, dim=3)

    ! Global arrays required, because transdens is difficult to parallelize 
  #:if WITH_SCALAPACK
    
    allocate(eigVecGlb(orb%nOrb,orb%nOrb,nSpin))
    allocate(ovrXevGlb(orb%nOrb,orb%nOrb,nSpin))
    do ss = 1, nSpin
      call distrib2replicated(env%blacs%orbitalGrid, denseDesc%blacsOrbSqr, eigVec(:,:,ss),&
          & eigVecGlb(:,:,ss))
      call distrib2replicated(env%blacs%orbitalGrid, denseDesc%blacsOrbSqr, ovrXev(:,:,ss),&
          & ovrXevGlb(:,:,ss))
    end do
    
  #:endif  
    
    vout(:) = 0.0_dp
    otmp(:,:,:,:)  = 0.0_dp

    fact = 1.0_dp
    if (.not. lr%tSpin) then
      if (sym == 'T') then
        fact = -1.0_dp
      end if
    end if

    do iAt = 1, lr%nAtom
      iSp = species0(iAt)
      nOrb = orb%nOrbAtom(iAt)
      call getOnsME(orb, iSp, lr%OnsiteMatrixElements, nOrb, onsite)
      do ia = 1, nmat
        iaGlb = ia + iOff
        call indxov(rpa%win, iaGlb, rpa%getIA, ii, jj, ss)
        updwn = (rpa%win(iaGlb) <= rpa%nxov_ud(1))
        
      #:if WITH_SCALAPACK
        
        call transDens(ii, jj, iAt, denseDesc%iAtomStart, nOrb, updwn, ovrXevGlb, eigVecGlb, qq_ij)
        
      #:else
        
        call transDens(ii, jj, iAt, denseDesc%iAtomStart, nOrb, updwn, ovrXev, eigVec, qq_ij)
        
      #:endif
        
        if (lr%tSpin) then
          if (.not. updwn) then
            sindx = [2, 1]
          else
            sindx = [1, 2]
          end if
          do iSpin = 1, 2
            otmp(:nOrb, :nOrb, iAt, iSpin) = otmp(:nOrb, :nOrb, iAt, iSpin) + qq_ij(:nOrb, :nOrb)&
                & * rpa%sqrOccIA(iaGlb) * vin(ia) * onsite(:nOrb, :nOrb, sindx(iSpin))
          end do
        else
          ! closed shell
          otmp(:nOrb, :nOrb, iAt, 1) = otmp(:nOrb, :nOrb, iAt, 1) + &
              & qq_ij(:nOrb, :nOrb) * rpa%sqrOccIA(iaGlb) * vin(ia)&
              & * (onsite(:nOrb, :nOrb, 1) + fact * onsite(:nOrb, :nOrb, 2))
        end if
      end do

      ! Rotational invariance corection for diagonal part
      do iSpin = 1, nSpin
        do iSh = 1, orb%nShell(iSp)
          degeneracy = real(2 * orb%angShell(iSh, iSp) + 1, dp)
          partTrace = 0.0_dp
          do iOrb = orb%posShell(iSh, iSp), orb%posShell(iSh + 1, iSp) - 1
            partTrace = partTrace + otmp(iOrb, iOrb, iAt, iSpin)
          end do
          partTrace = partTrace / degeneracy
          do iOrb = orb%posShell(iSh, iSp), orb%posShell(iSh + 1, iSp) - 1
            otmp(iOrb, iOrb, iAt, iSpin) = otmp(iOrb, iOrb, iAt, iSpin) - partTrace
          end do
        end do
      end do
    end do

    call assembleChunks(env, otmp)
     
    do ia = 1, nmat
      iaGlb = ia + iOff
      call indxov(rpa%win, iaGlb, rpa%getIA, ii, jj, ss)
      updwn = (rpa%win(iaGlb) <= rpa%nxov_ud(1))
      do iAt = 1, lr%nAtom
        nOrb = orb%nOrbAtom(iAt)

      #:if WITH_SCALAPACK
        
        call transDens(ii, jj, iAt, denseDesc%iAtomStart, nOrb, updwn, ovrXevGlb, eigVecGlb, qq_ij)

      #:else

        call transDens(ii, jj, iAt, denseDesc%iAtomStart, nOrb, updwn, ovrXev, eigVec, qq_ij)

      #:endif        
        
        vout(ia) = vout(ia) + 4.0_dp * rpa%sqrOccIA(iaGlb) * &
            & sum(qq_ij(:nOrb, :nOrb) * otmp(:nOrb, :nOrb, iAt, ss))
      end do
    end do

  end subroutine onsiteEner  

  !> Calculating spin polarized excitations.
  !! Note: the subroutine is generalized to account for spin and partial occupancy
  subroutine getSPExcitations(nOcc, nVir, grndEigVal, filling, wij, getIA, getIJ, getAB)

    !> Number of occupied states per spin channel
    integer, intent(in) :: nOcc(:)

    !> Number of virtual states per spin channel
    integer, intent(in) :: nVir(:)

    !> Ground state eigenvalues
    real(dp), intent(in) :: grndEigVal(:,:)

    !> Occupation numbers for the ground state
    real(dp), intent(in) :: filling(:,:)

    !> Kohn-Sham energy differences between empty and filled states
    real(dp), intent(out) :: wij(:)

    !> Index of occ-vir pairs of KS states for the transitions in wij
    integer, intent(out) :: getIA(:,:)

    !> Index of occ-occ pairs of KS states for the transitions in wij
    integer, intent(out) :: getIJ(:,:)

    !> Index of vir-vir pairs of KS states for the transitions in wij
    integer, intent(out) :: getAB(:,:)

    integer :: ind, ii, jj, aa, bb, off
    integer :: norb, iSpin, nSpin

    @:ASSERT(all(shape(grndEigVal) == shape(filling)))

    norb = size(grndEigVal, dim=1)
    nSpin = size(grndEigVal, dim=2)

    ind = 0

    do iSpin = 1, nSpin
      do ii = 1, norb - 1
        do aa = ii, norb
          if (filling(ii, iSpin) > filling(aa, iSpin) + elecTolMax) then
            ind = ind + 1
            wij(ind) = grndEigVal(aa, iSpin) - grndEigVal(ii, iSpin)
            getIA(ind,:) = [ii, aa, iSpin]
          end if
        end do
      end do
    end do

    off = 0

    do iSpin = 1, nSpin
      do ii = 1, nOcc(iSpin)
        do jj = 1, ii
          ind = jj + ((ii - 1) * ii) / 2 + off
          getIJ(ind,:) = [ii, jj, iSpin]
        end do
      end do
      off = (nOcc(1) * (nOcc(1) + 1)) / 2
    end do

    off = 0

    do iSpin = 1, nSpin
      do aa = 1, norb - nOcc(iSpin)
        do bb = 1, aa
          ind = bb + ((aa  - 1) * aa) / 2 + off
          getAB(ind,:) = [aa + nOcc(iSpin), bb + nOcc(iSpin), iSpin]
        end do
      end do
      off = (nVir(1) * (nVir(1) + 1)) / 2
    end do

  end subroutine getSPExcitations


  !> Calculate dipole moment for transitions.
  function oscillatorStrength(tSpin, transd, omega, xpy, sqrOccIA) result(osz)

    !> Spin-polarized calculation?
    logical, intent(in) :: tSpin

    !> Transition dipole matrix elements for each atom
    real(dp), intent(in) :: transd(:,:)

    !> Excitation energy
    real(dp), intent(in) :: omega

    !> Excitation eigenvector (X+Y) for state in question
    real(dp), intent(in) :: xpy(:)

    !> Square root of KS occupation differences
    real(dp), intent(in) :: sqrOccIA(:)

    !> Oscillator strength
    real(dp) :: osz

    real(dp) :: rtmp
    integer :: ii

    osz = 0.0_dp
    do ii = 1, 3
      rtmp = sum(transd(:,ii) * sqrOccIA * xpy)
      osz = osz + rtmp * rtmp
    end do
    osz = twothird * omega * osz

    ! For spin-unpolarized systems
    ! (X+Y)_ia_up = (X+Y)_ia_dn = (X+Y)_ia^Singlet / sqrt(2)
    if (.not. tSpin) then
      osz = osz * 2.0_dp
    end if

  end function oscillatorStrength


  !> Transition dipoles between single particle states.
  subroutine transitionDipole(tSpin, transd, xpy, sqrOccIA, tdip)

    !> Spin-polarized calculation?
    logical, intent(in) :: tSpin

    !> Transition dipole for transitions between KS states
    real(dp), intent(in) :: transd(:,:)

    !> Eigenvectors RPA equations (X+Y)
    real(dp), intent(in) :: xpy(:,:)

    !> Square root of KS occupation differences
    real(dp), intent(in) :: sqrOccIA(:)

    !> Resulting transition dipoles
    real(dp), intent(out) :: tdip(:,:)

    integer :: ii, ll

    tdip(:,:) = 0.0_dp
    !$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(ii) SCHEDULE(RUNTIME)
    do ii = 1, size(xpy, dim=2)
      do ll = 1, 3
        tdip(ii, ll) = sum(transd(:,ll) * sqrOccIA * xpy(:,ii))
      end do
    end do
    !$OMP END PARALLEL DO

    ! For spin-unpolarized systems
    ! (X+Y)_ia_up = (X+Y)_ia_dn = (X+Y)_ia^Singlet / sqrt(2)
    if (.not. tSpin) then
      tdip(:,:) = tdip * sqrt(2.0_dp)
    end if

  end subroutine transitionDipole


  !> Calculate transition moments for transitions between Kohn-Sham states, including spin-flipping
  !! transitions
  subroutine calcTransitionDipoles(coord0, win, getIA, transChrg, env, denseDesc, ovrXev,&
      & grndEigVecs, snglPartTransDip)

    !> Atomic positions
    real(dp), intent(in) :: coord0(:,:)

    !> Single particle transition index
    integer, intent(in) :: win(:) 

    !> Environment settings
    type(TEnvironment), intent(inout) :: env

    !> Dense matrix descriptor
    type(TDenseDescr), intent(in) :: denseDesc 

    !> Index array for excitation pairs
    integer, intent(in) :: getIA(:,:)

    !> machinery for transition charges between single particle levels
    type(TTransCharges), intent(in) :: transChrg

    !> overlap times ground state wavefunctions
    real(dp), intent(in) :: ovrXev(:,:,:)

    !> Ground state wavefunctions
    real(dp), intent(in) :: grndEigVecs(:,:,:)

    !> Resulting transition dipoles
    real(dp), intent(out) :: snglPartTransDip(:,:)

    integer :: nxov, natom, indm
    real(dp), allocatable :: qij(:)

    nxov = size(win)
    natom = size(coord0, dim=2)

    allocate(qij(natom))

    ! Calculate transition dipole elements
    do indm = 1, nxov
      qij(:) = transChrg%qTransIA(indm, env, denseDesc, ovrXev, grndEigVecs, getIA, win)
      snglPartTransDip(indm, :) = matmul(coord0, qij)
    end do

  end subroutine calcTransitionDipoles

  !> Calculate <S^2> as a measure of spin contamination (smaller magnitudes are better, 0.5 is
  !! considered an upper threshold for reliability according to Garcia thesis).
  subroutine getExcSpin(env, orb, rpa, denseDesc, Ssq, xpy, filling, ovrXev, eigVec)

    !> Environment settings
    type(TEnvironment), intent(inout) :: env

    !> Data type for atomic orbital information
    type(TOrbitals), intent(in) :: orb
    
    !> Run time parameters of the Casida routine
    type(TCasidaParameter), intent(in) :: rpa

    !> Dense matrix descriptor
    type(TDenseDescr), intent(in) :: denseDesc
    
    !> Spin contamination
    real(dp), intent(out) :: Ssq(:)

    !> Casida excited eigenvectors (X+Y)
    real(dp), intent(in) :: xpy(:,:)

    !> Occupations in ground state
    real(dp), intent(in) :: filling(:,:)

    !> Overlap times ground state eigenvectors
    real(dp), intent(in) :: ovrXev(:,:,:)

    !> Ground state eigenvectors
    real(dp), intent(in) :: eigVec(:,:,:)

    integer, allocatable :: TDvin(:)
    integer:: i, k, l, m, ia, jb, ii, aa, jj, bb, ss
    integer:: nSpin, nmat, nexc, nup, ndwn
    real(dp), allocatable :: TDvec(:), TDvec_sq(:)
    real(dp), allocatable :: eigVecGlb(:,:,:), ovrXevGlb(:,:,:)
    real(dp) :: TDvnorm
    real(dp) :: s_iaja, s_iaib, s_iajb, tmp
    logical :: ud_ia, ud_jb
    
    nmat = size(xpy, dim=1)
    nexc = size(Ssq)
    nup = ceiling(sum(filling(:,1)))
    ndwn = ceiling(sum(filling(:,2)))
    nSpin = size(eigVec, dim=3)
    
    allocate(TDvec(nmat))
    allocate(TDvec_sq(nmat))
    allocate(TDvin(nmat))

    ! Global arrays required, routine is difficult to parallelize
    allocate(eigVecGlb(orb%nOrb,orb%nOrb,nSpin))
    allocate(ovrXevGlb(orb%nOrb,orb%nOrb,nSpin))
    
  #:if WITH_SCALAPACK
    
    do ss = 1, nSpin
      call distrib2replicated(env%blacs%orbitalGrid, denseDesc%blacsOrbSqr, eigVec(:,:,ss),&
          & eigVecGlb(:,:,ss))
      call distrib2replicated(env%blacs%orbitalGrid, denseDesc%blacsOrbSqr, ovrXev(:,:,ss),&
          & ovrXevGlb(:,:,ss))
    end do

  #:else   

    eigVecGlb(:,:,:) = eigVec
    ovrXevGlb(:,:,:) = ovrXev
    
  #:endif      

    do i = 1, nexc
      TDvec(:) = xpy(:,i)
      TDvnorm = 1.0_dp / sqrt(sum(TDvec**2))
      TDvec(:) = TDvec * TDvnorm
      TDvec_sq = TDvec**2

      ! Put these transition dipoles in order of descending magnitude
      call index_heap_sort(TDvin, TDvec_sq)
      TDvin = TDvin(nmat:1:-1)
      TDvec_sq = TDvec_sq(TDvin)

      ! S_{ia,ja}
      s_iaja = 0.0_dp
      do k = 1, nmat
        ia = TDvin(k)
        call indxov(rpa%win, ia, rpa%getIA, ii, aa, ss)
        ud_ia = (rpa%win(ia) <= rpa%nxov_ud(1))
        do l = 1, nmat
          jb = TDvin(l)
          call indxov(rpa%win, jb, rpa%getIA, jj, bb, ss)
          ud_jb = (rpa%win(jb) <= rpa%nxov_ud(1))

          if ((bb /= aa) .or. (ud_jb .neqv. ud_ia)) cycle

          tmp = 0.0_dp
          if (ud_ia) then
            do m = 1, ndwn
              tmp = tmp + MOoverlap(ii, m, ovrXevGlb, eigVecGlb)&
                  & * MOoverlap(jj, m, ovrXevGlb, eigVecGlb)
            end do
          else
            do m = 1, nup
              tmp = tmp + MOoverlap(m, ii, ovrXevGlb, eigVecGlb)&
                  & * MOoverlap(m, jj, ovrXevGlb, eigVecGlb)
            end do
          end if

          s_iaja = s_iaja + TDvec(ia) * TDvec(jb) * tmp
        end do
      end do

      ! S_{ia,ib}
      s_iaib = 0.0_dp
      do k = 1, nmat
        ia = TDvin(k)
        call indxov(rpa%win, ia, rpa%getIA, ii, aa, ss)
        ud_ia = (rpa%win(ia) <= rpa%nxov_ud(1))
        do l = 1, nmat
          jb = TDvin(l)
          call indxov(rpa%win, jb, rpa%getIA, jj, bb, ss)
          ud_jb = (rpa%win(jb) <= rpa%nxov_ud(1))

          if ((ii /= jj) .or. (ud_jb .neqv. ud_ia)) cycle

          tmp = 0.0_dp
          if (ud_ia) then
            do m = 1, ndwn
              tmp = tmp + MOoverlap(aa, m, ovrXevGlb, eigVecGlb)&
                  & * MOoverlap(bb, m, ovrXevGlb, eigVecGlb)
            end do
          else
            do m = 1, nup
              tmp = tmp + MOoverlap(m, aa, ovrXevGlb, eigVecGlb)&
                  & * MOoverlap(m, bb, ovrXevGlb, eigVecGlb)
            end do
          end if

         s_iaib = s_iaib + TDvec(ia) * TDvec(jb) * tmp
        end do
      end do

      ! S_{ia,jb}
      s_iajb = 0.0_dp
      do k = 1, nmat
        ia = TDvin(k)
        call indxov(rpa%win, ia, rpa%getIA, ii, aa, ss)
        ud_ia = (rpa%win(ia) <= rpa%nxov_ud(1))
        if (.not. ud_ia) cycle
        do l = 1, nmat
          jb = TDvin(l)
          call indxov(rpa%win, jb, rpa%getIA, jj, bb, ss)
          ud_jb = (rpa%win(jb) <= rpa%nxov_ud(1))

          if (ud_jb) cycle

          s_iajb = s_iajb + TDvec(ia) * TDvec(jb) * MOoverlap(aa, bb, ovrXevGlb, eigVecGlb)&
              & * MOoverlap(ii, jj, ovrXevGlb, eigVecGlb)
        end do
      end do

      Ssq(i) =  s_iaja - s_iaib - 2.0_dp * s_iajb

    end do

  end subroutine getExcSpin

  !> Write single particle excitations to a file.
  subroutine writeSPExcitations(lr, rpa, sposz)

    !> Data structure for linear response
    type(TLinResp), intent(in) :: lr

    !> Run time parameters of the Casida routine
    type(TCasidaParameter), intent(in) :: rpa

    !> Single particle oscilation strengths
    real(dp), intent(in) :: sposz(:)

    integer :: indm, m, n, s
    logical :: updwn
    character :: sign
    type(TFileDescr) :: fdSPTrans

    @:ASSERT(size(sposz) >= rpa%nxov_rd)

    if (lr%writeSPTrans) then
      ! Single particle excitations
      call openFile(fdSPTrans, singlePartOut, mode="w")
      write(fdSPTrans%unit, *)
      write(fdSPTrans%unit, '(7x,a,7x,a,8x,a)') '#      w [eV]', 'Osc.Str.', 'Transition'
      write(fdSPTrans%unit, *)
      write(fdSPTrans%unit, '(1x,58("="))')
      write(fdSPTrans%unit, *)
      do indm = 1, rpa%nxov_rd
        call indxov(rpa%win, indm, rpa%getIA, m, n, s)
        sign = " "
        if (lr%tSpin) then
          updwn = (rpa%win(indm) <= rpa%nxov_ud(1))
          if (updwn) then
            sign = "U"
          else
            sign = "D"
          end if
        end if
        write(fdSPTrans%unit,&
            & '(1x,i7,3x,f8.3,3x,f13.7,4x,i5,3x,a,1x,i5,1x,1a)')&
            & indm, Hartree__eV * rpa%wij(indm), sposz(indm), m, '->', n, sign
      end do
      write(fdSPTrans%unit, *)
      call closeFile(fdSPTrans)
    end if

  end subroutine writeSPExcitations


  !> Excited state Mulliken charges and dipole moments written to disc.
  subroutine writeExcMulliken(sym, nstat, dq, dqex, coord0)

    !> Symmetry label
    character, intent(in) :: sym

    !> State index
    integer, intent(in) :: nstat

    !> Ground state gross charge
    real(dp), intent(in) :: dq(:)

    !> Change in atomic charges from ground to excited state
    real(dp), intent(in) :: dqex(:)

    !> Central cell coordinates
    real(dp), intent(in) :: coord0(:,:)

    type(TFileDescr) :: fdMulliken
    integer :: natom, m
    real(dp) :: dipol(3), dipabs

    natom = size(dq)

    @:ASSERT(size(dq) == size(dqex))
    @:ASSERT(all(shape(coord0) == [3, nAtom]))

    ! Output of excited state Mulliken charges
    call openFile(fdMulliken, excitedQOut, mode="a")
    write(fdMulliken%unit, "(a,a,i2)") "# MULLIKEN CHARGES of excited state ", sym, nstat
    write(fdMulliken%unit, "(a,2x,A,i4)") "#", 'Natoms =',natom
    write(fdMulliken%unit, "('#',1X,A4,T15,A)")'Atom','netCharge'
    write(fdMulliken%unit,'("#",41("="))')
    do m = 1, natom
      write(fdMulliken%unit,"(i5,1x,f16.8)") m, -dq(m) - dqex(m)
    end do
    call closeFile(fdMulliken)

    ! Calculation of excited state dipole moment
    dipol(:) = -1.0_dp * matmul(coord0, dq + dqex)
    dipabs = sqrt(sum(dipol**2))

    call openFile(fdMulliken, excitedDipoleOut, mode="a")
    write(fdMulliken%unit, "(a,a,i2)") "Mulliken analysis of excited state ", sym, nstat
    write(fdMulliken%unit, '(42("="))')
    write(fdMulliken%unit, "(a)") " "
    write(fdMulliken%unit, "(a)") "Mulliken exc. state dipole moment [Debye]"
    write(fdMulliken%unit, '(42("="))')
    write(fdMulliken%unit, "(3f14.8)") (dipol(m) * au__Debye, m = 1, 3)
    write(fdMulliken%unit, "(a)") " "
    write(fdMulliken%unit, "(a)") "Norm of exc. state dipole moment [Debye]"
    write(fdMulliken%unit, '(42("="))')
    write(fdMulliken%unit, "(e20.12)") dipabs * au__Debye
    write(fdMulliken%unit, *)
    call closeFile(fdMulliken)

  end subroutine writeExcMulliken


  !> Increase dimension of vector 
  pure subroutine incSizeVec(oldDim, newDim, vec)

    !> Size of the input vector to copy over to resized vector
    integer, intent(in) :: oldDim

    !> New size of vector 
    integer, intent(in) :: newDim

    !> Vector to re-size, retaining initial elements in output
    real(dp), allocatable, intent(inout) :: vec(:)

    real(dp), allocatable :: temp(:)

    allocate(temp(newDim))
    temp(:) = 0.0_dp
    temp(1:oldDim) = vec
    call move_alloc(temp, vec)

  end subroutine incSizeVec


  !> Increase size of (oldDim, n) array to (newDim, n).
  pure subroutine incSizeMatDimOne(oldDim, newDim, mat)

    !> Size of the input matrix first dimension to copy over to resized matrix
    integer, intent(in) :: oldDim

    !> New dimension 
    integer, intent(in) :: newDim

    !> Matrix to re-size, retaining initial elements in output
    real(dp), allocatable, intent(inout) :: mat(:,:)

    integer :: dim2
    real(dp), allocatable :: temp(:,:)

    dim2 = size(mat, dim=2)
    allocate(temp(newDim, dim2))
    temp(:,:) = 0.0_dp
    temp(1:oldDim,:) = mat
    call move_alloc(temp, mat)

  end subroutine incSizeMatDimOne


  !> Increase size of (n, oldDim) array to (n, newDim).
  pure subroutine incSizeMatDimTwo(oldDim, newDim, mat)

    !> Size of the input matrix second dimension to copy over to resized matrix
    integer, intent(in) :: oldDim

    !> New dimension
    integer, intent(in) :: newDim

    !> Matrix to re-size, retaining initial elements in output
    real(dp), allocatable, intent(inout) :: mat(:,:)

    integer :: dim1
    real(dp), allocatable :: temp(:,:)

    dim1 = size(mat, dim=1)
    allocate(temp(dim1, newDim))
    temp = 0.0_dp
    temp(:,1:oldDim) = mat
    call move_alloc(temp, mat)

  end subroutine incSizeMatDimTwo


  !> Increase size of (oldDim, oldDim) square array to (newDim, newDim).
  pure subroutine incSizeMatBothDim(oldDim, newDim, mat)

    !> Size of the input matrix second dimension to copy over to resized matrix (square)
    integer, intent(in) :: oldDim

    !> New dimension
    integer, intent(in) :: newDim

    !> Matrix to re-size, retaining initial elements in output
    real(dp), allocatable, intent(inout) :: mat(:,:)

    real(dp), allocatable :: temp(:,:)

    allocate(temp(newDim, newDim))
    temp(:,:) = 0.0_dp
    temp(1:oldDim, 1:oldDim) = mat
    call move_alloc(temp, mat)

  end subroutine incSizeMatBothDim


  !> Calculate square root and inverse of sqrt of a real, symmetric positive definite matrix.
  subroutine calcMatrixSqrt(matIn, spaceDim, matOut, matInvOut)

    !> Matrix to operate on
    real(dp), intent(in) :: matIn(:,:)

    !> Dimensions of input matrix
    integer, intent(in) :: spaceDim

    !> Matrix square root
    real(dp), intent(out) :: matOut(:,:)

    !> Inverse of matrix square root
    real(dp), intent(out) :: matInvOut(:,:)

    real(dp) :: dummyEV(spaceDim)
    real(dp) :: dummyM(spaceDim, spaceDim), dummyM2(spaceDim, spaceDim)
    integer :: ii

    dummyM(:,:) = matIn

    call heev(dummyM, dummyEV, 'U', 'V') 

    ! Calc. sqrt
    do ii = 1, spaceDim
      dummyM2(:,ii) = sqrt(dummyEV(ii)) * dummyM(:,ii)
    end do

    call gemm(matOut, dummyM2, dummyM, transB='T')

    ! Calc. inv. of sqrt
    do ii = 1, spaceDim
      dummyM2(:,ii) = dummyM(:,ii) / sqrt(dummyEV(ii))
    end do

    call gemm(matInvOut, dummyM2, dummyM, transB='T')

  end subroutine calcMatrixSqrt


  !> Perform modified Gram-Schmidt orthonormalization of vectors in columns of vec(1:end). Assume
  !! vectors 1:(start-1) are already orthonormal
  subroutine orthonormalizeVectors(env, iStart, iEnd, vec)

    !> Environment settings
    type(TEnvironment), intent(in) :: env

    !> Starting place in vectors to work from
    integer, intent(in) :: iStart

    !> Ending place in vectors
    integer, intent(in) :: iEnd

    !> Vectors to be orthogonalized against 1:end vectors
    real(dp), intent(inout) :: vec(:,:)

    integer :: ii, jj
    real(dp) :: dummyReal

    ! Obviously, not optimal in terms of communication, can be optimized if necessary
    do ii = iStart, iEnd
      do jj = 1, ii - 1
        dummyReal = dot_product(vec(:,ii), vec(:,jj))
        call assembleChunks(env, dummyReal)
        vec(:,ii) = vec(:,ii) - dummyReal * vec(:,jj)
      end do
      dummyReal = dot_product(vec(:,ii), vec(:,ii))
      call assembleChunks(env, dummyReal)
      vec(:,ii) = vec(:,ii) / sqrt(dummyReal)
    end do
    
  end subroutine orthonormalizeVectors


  !> Encapsulate memory expansion for Stratmann solver.
  subroutine incMemStratmann(oldDim, newDim, vecB, vP, vM, mP, mM, mH, mMsqrt, mMsqrtInv, dummyM,&
      & evalInt, evecL, evecR)

    !> Previous size of subspace
    integer, intent(in) :: oldDim

    !> New size of subspace
    integer, intent(in) :: newDim

    !> Basis of subspace
    real(dp), allocatable, intent(inout) :: vecB(:,:)

    !> Vec. for (A+B)b_i
    real(dp), allocatable, intent(inout) :: vP(:,:)

    !> Vec. for (A-B)b_i
    real(dp), allocatable, intent(inout) :: vM(:,:)

    !> M_plus
    real(dp), allocatable, intent(inout) :: mP(:,:)

    !> M_minus
    real(dp), allocatable, intent(inout) :: mM(:,:)

    !> M_herm matrix
    real(dp), allocatable, intent(inout) :: mH(:,:)

    !> Matrix square root
    real(dp), allocatable, intent(inout) :: mMsqrt(:,:)

    !> Inverse of matrix sqrt
    real(dp), allocatable, intent(inout) :: mMsqrtInv(:,:)

    !> Work array
    real(dp), allocatable, intent(inout) :: dummyM(:,:)

    !> Internal eigenvector storage
    real(dp), allocatable, intent(inout) :: evalInt(:)

    !> Left eigenvectors
    real(dp), allocatable, intent(inout) :: evecL(:,:)

    !> Right eigenvectors
    real(dp), allocatable, intent(inout) :: evecR(:,:)

    call incSizeMatDimTwo(oldDim, newDim, vecB)
    call incSizeMatDimTwo(oldDim, newDim, vP)
    call incSizeMatDimTwo(oldDim, newDim, vM)
    call incSizeMatBothDim(oldDim, newDim, mP)
    call incSizeMatBothDim(oldDim, newDim, mM)
    call incSizeMatBothDim(oldDim, newDim, mH)
    call incSizeMatBothDim(oldDim, newDim, mMsqrt)
    call incSizeMatBothDim(oldDim, newDim, mMsqrtInv)
    call incSizeMatBothDim(oldDim, newDim, dummyM)
    call incSizeVec(oldDim, newDim, evalInt)
    call incSizeMatDimOne(oldDim, newDim, evecL)
    call incSizeMatDimOne(oldDim, newDim, evecR)

  end subroutine incMemStratmann

end module dftbp_timedep_linrespcommon
