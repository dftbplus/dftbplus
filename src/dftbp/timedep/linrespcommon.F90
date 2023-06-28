!--------------------------------------------------------------------------------------------------!
!  DFTB+: general package for performing fast atomistic simulations                                !
!  Copyright (C) 2006 - 2023  DFTB+ developers group                                               !
!                                                                                                  !
!  See the LICENSE file for terms of usage and distribution.                                       !
!--------------------------------------------------------------------------------------------------!

#:include 'common.fypp'

!> Helper routines for the linear response modules.
module dftbp_timedep_linrespcommon
  use mpi
  use dftbp_common_accuracy, only : dp, elecTolMax
  use dftbp_common_constants, only: Hartree__eV, au__Debye, cExchange
  use dftbp_common_file, only : TFileDescr, openFile, closeFile
  use dftbp_dftb_onsitecorrection, only : getOnsME
  use dftbp_io_message, only : error
  use dftbp_math_blasroutines, only : hemv
  use dftbp_math_sorting, only : index_heap_sort
  use dftbp_timedep_transcharges, only : TTransCharges, transq
  use dftbp_type_commontypes, only : TOrbitals
  use dftbp_type_densedescr, only : TDenseDescr
  use dftbp_common_environment, only : TEnvironment
  implicit none

  public

  !> prefactor of 2/3.
  real(dp), parameter :: twothird = 2.0_dp / 3.0_dp

  !> Names of output files
  character(*), parameter :: excitedQOut = "XCH.DAT"
  character(*), parameter :: excitedDipoleOut = "XREST.DAT"
  character(*), parameter :: singlePartOut = "SPX.DAT"

contains


  !> find (possibly degenerate) transitions with stronger dipole
  !> transition strengths than a tolerance, count them and place at
  !> the start of the appropriate arrays
  subroutine dipSelect(wij,sposz,win,transd,nxov_r,threshold,grndEigVal, getIA)

    !> Energies of transitions
    real(dp), intent(inout) :: wij(:)

    !> transition strengths
    real(dp), intent(inout) :: sposz(:)

    !> index array for transitions to single particle transitions
    integer, intent(inout) :: win(:)

    !> transition dipole matrix elements for each atom
    real(dp), intent(inout) :: transd(:,:)

    !> number of excitations in space
    integer, intent(out) :: nxov_r

    !> threshold for dipole cutoff
    real(dp), intent(in) :: threshold

    !> ground state eigenvalues
    real(dp), intent(in) :: grndEigVal(:,:)

    !> Index array of transitions
    integer, intent(in) :: getIA(:,:)

    integer :: nxov, ii, jj, iOcc, iVrt, iStart, iSpin
    real(dp) :: eOcc, eExc, mu

    nxov = size(wij)
    @:ASSERT(nxov == size(sposz))
    @:ASSERT(nxov == size(win))
    @:ASSERT(all(shape(transd) == [nxov,3]))
    @:ASSERT(size(getIA,dim=1) >= nxov)
    @:ASSERT(size(getIA,dim=2) == 3)
    @:ASSERT(threshold >= 0.0_dp)

    call indxov(win, 1, getIA, iOcc, iVrt, iSpin)
    eOcc = grndEigVal(iOcc,1)
    eExc = wij(1)
    iStart = 1
    nxov_r = 0
    ! Check, to single precision tolerance, for degeneracies when selecting bright transitions
    do ii = 2, nxov
      call indxov(win, ii, getIA, iOcc, iVrt, iSpin)
      ! check if this is a still within a degenerate group, otherwise process the group
      if ( abs(grndEigVal(iOcc,1)-eOcc) > epsilon(0.0) .or. &
          & abs(wij(ii)-eExc) > epsilon(0.0) ) then
        eOcc = grndEigVal(iOcc,1)
        eExc = wij(ii)
        mu = 0.0_dp
        ! loop over transitions in the group and check the oscillator strength
        do jj = iStart, ii-1
          call indxov(win, jj, getIA, iOcc, iVrt, iSpin)
          mu = mu + sposz(jj)
        end do
        ! if something in the group is bright, so include them all
        if (mu>threshold) then
          do jj = iStart, ii-1
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

    ! last group in the transitions
    mu = 0.0_dp
    do jj = iStart, nxov
      call indxov(win, jj, getIA, iOcc, iVrt, iSpin)
      mu = mu + sposz(jj)
    end do
    if (mu>threshold) then
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

    !> index array after sorting of eigenvalues.
    integer, intent(in) :: win(:)

    !> Compound excitation index.
    integer, intent(in) :: indx

    !> Index array of transitions after sorting of eigenvalues
    integer, intent(in) :: getIA(:,:)

    !> Initial (filled) state.
    integer, intent(out) :: ii

    !> Final (empty) state.
    integer, intent(out) :: jj

    !> Spin channel: 1 (up) or 2 (down)
    integer, intent(out) :: ss

    integer :: indo

    indo = win(indx)
    ii = getIA(indo,1)
    jj = getIA(indo,2)
    ss = getIA(indo,3)

  end subroutine indxov


  !> Computes individual indices from the compound occ-occ excitation index.
  subroutine indxoo(indx, ii, jj)

    !> Compound excitation index.
    integer, intent(in) :: indx

    !> Initial state.
    integer, intent(out) :: ii

    !> Final state.
    integer, intent(out) :: jj

    call indxvv(0, indx, ii, jj)

  end subroutine indxoo


  !> Computes individual indices from the compound virtual-virtual excitation index.
  subroutine indxvv(nocc, indx, ii, jj)

    !> Compund excitation index.
    integer, intent(in) :: nocc

    !> Number of occupied states.
    integer, intent(in) :: indx

    !> Initial state.
    integer, intent(out) :: ii

    !> Final state.
    integer, intent(out) :: jj

    real(dp) :: re

    @:ASSERT(indx > 0)

    ! solve a quadratic for the row of a matrix between occupied and
    ! virtual states, given the number of the element in the matrix
    re = sqrt(real(indx, dp) * 2.0_dp - 1.75_dp) - 0.5_dp

    ii  = floor(re) + 1 + nocc ! actual row
    jj  = indx - ((ii - 1 - nocc) * (ii - nocc)) / 2 + nocc ! column

    @:ASSERT(ii > 0)
    @:ASSERT(jj > 0)

  end subroutine indxvv


  !> Builds array to convert from orbital pairs to a compound index, reverse of index generated by
  !> getSPExcitations (which were then potentially windowed afterwards)
  subroutine rindxov_array(win, nxov, nxoo, nxvv, getIA, getIJ, getAB, iatrans)

    !> array for indexing excitations
    integer, intent(in) :: win(:)

    !> Number of transitions from occupied to virtual
    integer, intent(in) :: nxov

    !> Number of transitions from occupied to occupied
    integer, intent(in) :: nxoo

    !> Number of transitions from virtual to virtual
    integer, intent(in) :: nxvv

    !> array of the occupied->virtual pairs
    integer, intent(in) :: getIA(:,:)

    !> array of the occupied->occupied pairs
    integer, intent(in) :: getIJ(:,:)

    !> array of the virtual->virtual pairs
    integer, intent(in) :: getAB(:,:)

    !> resulting index array from orbital pairs to compound index
    integer, intent(out) :: iatrans(:,:,:)

    integer :: ia, ij, ab

    @:ASSERT(size(getIA,dim=2) == 3)
    @:ASSERT(size(getIJ,dim=2) == 3)
    @:ASSERT(size(getAB,dim=2) == 3)

    ! Store reverse indices

    iaTrans(:,:,:) = 0

    ! The transition charges q_rs are symmetrical in rs
    ! We only need one index per pair of orbitals
    !$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(ia) SCHEDULE(RUNTIME)
    do ia = 1, nxov
      iatrans( getIA(win(ia),1), getIA(win(ia),2), getIA(win(ia),3)) = ia
      iatrans( getIA(win(ia),2), getIA(win(ia),1), getIA(win(ia),3)) = ia
    end do
    !$OMP  END PARALLEL DO

    !$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(ij) SCHEDULE(RUNTIME)
    do ij = 1, nxoo
      iatrans( getIJ(ij,1), getIJ(ij,2), getIJ(ij,3)) = ij
      iatrans( getIJ(ij,2), getIJ(ij,1), getIJ(ij,3)) = ij
    end do
    !$OMP  END PARALLEL DO

    !$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(ab) SCHEDULE(RUNTIME)
    do ab = 1, nxvv
      iatrans( getAB(ab,1), getAB(ab,2), getAB(ab,3)) = ab
      iatrans( getAB(ab,2), getAB(ab,1), getAB(ab,3)) = ab
    end do
    !$OMP  END PARALLEL DO


  end subroutine rindxov_array


  !> calculate the transition block at a specific atom
  subroutine transDens(ii, jj, iAt, iAtomStart, nOrb, updwn, ovrXev, grndEigVecs, qq_ij)

    !> Index of initial state.
    integer, intent(in) :: ii

    !> Index of final state.
    integer, intent(in) :: jj

    !> Atom number
    integer, intent(in) :: iAt

    !> Starting position of each atom in the list of orbitals.
    integer, intent(in) :: iAtomStart(:)

    !> up spin channel (T) or down spin channel (F)
    logical, intent(in) :: updwn

    !> Overlap times eigenvector: sum_m Smn cmi (nOrb, nOrb).
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
    !call ger(qq_ij(:nOrb,:nOrb), 0.5_dp,grndEigVecs(mu:mu+nOrb-1,ii,ss),ovrXev(mu:mu+nOrb-1,jj,ss))
    !call ger(qq_ij(:nOrb,:nOrb), 0.5_dp,grndEigVecs(mu:mu+nOrb-1,jj,ss),ovrXev(mu:mu+nOrb-1,ii,ss))
    !qq_ij(:nOrb,:nOrb) = 0.5_dp * (qq_ij(:nOrb,:nOrb) + transpose(qq_ij(:nOrb,:nOrb)))

    do iOrb1 = 1, nOrb
      do iOrb2 = iOrb1, nOrb
        mu = iAtomStart(iAt) + iOrb1 - 1
        nu = iAtomStart(iAt) + iOrb2 - 1
        qq_ij(iOrb1,iOrb2) = 0.25_dp*( grndEigVecs(mu,ii,ss)*ovrXev(nu,jj,ss) &
             &                       + grndEigVecs(mu,jj,ss)*ovrXev(nu,ii,ss) &
             &                       + grndEigVecs(nu,ii,ss)*ovrXev(mu,jj,ss) &
             &                       + grndEigVecs(nu,jj,ss)*ovrXev(mu,ii,ss) )
        if (iOrb1 /= iOrb2) then
          qq_ij(iOrb2,iOrb1) = qq_ij(iOrb1,iOrb2)
        end if
      end do
    end do

  end subroutine transDens


  !> Returns the (spatial) MO overlap between orbitals in different spin channels
  function MOoverlap(pp, qq, ovrXev, grndEigVecs) result(S_pq)

    !> orbital in alpha channel
    integer, intent(in) :: pp

    !> orbital in  beta channel
    integer, intent(in) :: qq

    !> overlap times ground single particle state wavefunctiona
    real(dp), intent(in) :: ovrXev(:,:,:)

    !> ground state single particle wavefunctions
    real(dp), intent(in) :: grndEigVecs(:,:,:)

    !> Resulting overlap between states
    real(dp):: S_pq

    @:ASSERT(all(shape(grndEigVecs) == shape(ovrXev)))
    @:ASSERT(size(grndEigVecs, dim=3) == 2)

    S_pq = sum(grndEigVecs(:,pp,1)*ovrXev(:,qq,2))

  end function MOoverlap


  !> Square root of differences in occupations between filled and empty states
  !> We need occupations per spin channel ([0:1]) also for closed shell systems
  subroutine getSqrOcc(occNr, win, nmatup, nmat, getIA, tSpin, n_ij)

    !> occupation of states
    real(dp), intent(in) :: occNr(:,:)

    !> index array for transitions to single particle transitions
    integer, intent(in) :: win(:)

    !> number of up-up excitations
    integer, intent(in) :: nmatup

    !> number of transitions to scale
    integer, intent(in) :: nmat

    !> array of the occupied->virtual pairs
    integer, intent(in) :: getIA(:,:)

    !> is system spin polarized?
    logical, intent(in) :: tSpin

    !> resulting scaled matrix
    real(dp), intent(out) :: n_ij(:)

    integer :: ia, ii, jj, ss
    logical :: updwn
    real(dp) :: docc_ij

    !$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(ia, ii, jj, updwn, docc_ij) SCHEDULE(RUNTIME)
    do ia = 1, nmat
      call indxov(win, ia, getIA, ii, jj, ss)
      updwn = (win(ia) <= nmatup)
      if(updwn) then
        docc_ij = occNr(ii,1) - occNr(jj,1)
      else
        docc_ij = occNr(ii,2) - occNr(jj,2)
      end if
      if(.not. tSpin) then
         docc_ij = docc_ij / 2.0_dp
      endif
      n_ij(ia) = sqrt(docc_ij)
    end do
    !$OMP  END PARALLEL DO

  end subroutine getSqrOcc


  !> Multiplies the excitation supermatrix with a supervector.
  !> For the hermitian RPA eigenvalue problem this corresponds to \Omega_ias,jbt * v_jbt
  !> (spin polarized case) or \Omega^{S/T}_ia,jb * v_jb (singlet/triplet)
  !>
  !> For the standard RPA, (A+B)_ias,jbt * v_jbt needs to be computed (similar for singlet/triplet)
  !> (see definitions in Marc Casida, in Recent Advances in Density Functional Methods,
  !>  World Scientific, 1995, Part I, p. 155.)
  !> Note: we actually compute sqrt(n_is-n_as) (A+B)_ias,jbt sqrt(n_jt-n_bt), with the
  !> occupations n_is, correct for finite T.
  !> See also Dominguez JCTC 9 4901 (2013), Kranz JCTC 13 1737 (2017) for DFTB specifics.
  !>
  !> Note: In order not to store the entire supermatrix (nmat, nmat), the various pieces are
  !> assembled individually and multiplied directly with the corresponding part of the supervector.
  subroutine actionAplusB(spin, wij, sym, win, nocc_ud, nvir_ud, nxoo_ud, nxvv_ud, nxov_ud,&
      & nxov_rd, iaTrans, getIA, getIJ, getAB, env, denseDesc, ovrXev, grndEigVecs, occNr, sqrOccIA,&
      & gamma, species0, spinW, onsMEs, orb, tAplusB, transChrg, vin, vout, tRangeSep, lrGamma)

    !> logical spin polarization
    logical, intent(in) :: spin

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

    !> Environment settings
    type(TEnvironment), intent(in) :: env

    !> Dense matrix descriptor
    type(TDenseDescr), intent(in) :: denseDesc 

    !> overlap times eigenvector. (nOrb, nOrb)
    real(dp), intent(in) :: ovrXev(:,:,:)

    !> eigenvectors (nOrb, nOrb)
    real(dp), intent(in) :: grndEigVecs(:,:,:)

    !> occupation numbers
    real(dp), intent(in) :: occNr(:,:)

    !> Square root of occupation difference between vir and occ states
    real(dp), intent(in) :: sqrOccIA(:)

    !> DFTB gamma matrix (nAtm, nAtom)
    real(dp), intent(in) :: gamma(:,:)

    !> chemical species of the atoms
    integer, intent(in) :: species0(:)

    !> ground state spin constants for each species
    real(dp), intent(in) :: spinW(:)

    !> onsite matrix elements for shells (elements between s orbitals on the same shell are ignored)
    real(dp), intent(in), allocatable :: onsMEs(:,:,:,:)

    !> data type for atomic orbital information
    type(TOrbitals), intent(in) :: orb

    !> whether (A+B)v or (Omega)v should be returned
    logical, intent(in) :: tAplusB

    !> machinery for transition charges between single particle levels
    type(TTransCharges), intent(in) :: transChrg

    !> vector to multiply with size(nmat)
    real(dp), intent(in) :: vin(:)

    !> vector containing the result on exit size(nmat)
    real(dp), intent(out) :: vout(:)

    !> is calculation range-separated?
    logical, intent(in) :: tRangeSep

    !> long-range Gamma if in use
    real(dp), allocatable, optional, intent(in) :: lrGamma(:,:)

    integer :: nmat, natom, norb, ia, nxvv_a
    integer :: ii, aa, ss, jj, bb, ias, jbs, abs, ijs, jas, ibs

    ! somewhat ugly, but fast small arrays on stack:
    real(dp) :: otmp(size(gamma, dim=1)), gtmp(size(gamma, dim=1))
    real(dp) :: qij(size(gamma, dim=1)) ! qij Working array (used for excitation charges). (nAtom)
    real(dp) :: vout_ons(size(vin))
    real(dp) :: vTmp(size(vin))
    real(dp), dimension(2) :: spinFactor = (/ 1.0_dp, -1.0_dp /)
    ! for later use to change HFX contribution
    real(dp), allocatable :: qv(:,:)

    @:ASSERT(size(vin) == size(vout))

    !! dimension guessed from input vector
    nmat = size(vin)
    @:ASSERT(nmat <= nxov_rd)
    natom = size(gamma, dim=1)
    norb = size(ovrXev, dim=1)
    @:ASSERT(norb == nocc_ud(1)+nvir_ud(1))

    vout = 0.0_dp

    if(tAplusB) then
       vTmp(:) = vin
    else
       vTmp(:) = vin * sqrt(wij)
    endif

    ! product charges with the v*wn product, i.e. Q * v*wn
    oTmp(:) = 0.0_dp
    call transChrg%qMatVec(denseDesc, ovrXev, grndEigVecs, getIA, win, &
        & vTmp * sqrOccIA, oTmp)

    if (.not.spin) then !-----------spin-unpolarized systems--------------

      if (sym == 'S') then

        call hemv(gtmp, gamma, otmp)

        ! 4 * wn * (g * Q)
        vOut(:) = 0.0_dp
        call transChrg%qVecMat(denseDesc, ovrXev, grndEigVecs, getIA, win, gTmp, vOut)
        vOut(:) = 4.0_dp * sqrOccIA * vOut

      else

        otmp = otmp * spinW(species0)

        ! 4 * wn * (o * Q)
        vOut = 0.0_dp
        call transChrg%qVecMat(denseDesc, ovrXev, grndEigVecs, getIA, win, oTmp, vOut)
        vOut(:) = 4.0_dp * sqrOccIA * vOut


      end if

    else !--------------spin-polarized systems--------

      call hemv(gtmp, gamma, otmp)

      otmp(:) = 0.0_dp

      !$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(ia,ss,qij) &
      !$OMP& SCHEDULE(RUNTIME) REDUCTION(+:otmp)
      do ia = 1,nmat

        qij(:) = transChrg%qTransIA(ia, denseDesc, ovrXev, grndEigVecs, getIA, win)

        ! singlet gamma part (S)
        vout(ia) = 2.0_dp * sqrOccIA(ia) * dot_product(qij, gtmp)

        ! magnetization part (T1)
        ss = getIA(win(ia), 3)

        otmp(:) = otmp + spinFactor(ss) * sqrOccIA(ia) * vTmp(ia) * qij

      end do
      !$OMP  END PARALLEL DO

      otmp = otmp * spinW(species0)

      !$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(ia,ss,qij) &
      !$OMP& SCHEDULE(RUNTIME)
      do ia = 1,nmat

        qij(:) = transChrg%qTransIA(ia, denseDesc, ovrXev, grndEigVecs, getIA, win)

        ss = getIA(win(ia), 3)

        vout(ia) = vout(ia) + 2.0_dp * sqrOccIA(ia) * spinFactor(ss) * dot_product(qij, otmp)

      end do
      !$OMP  END PARALLEL DO

    end if

    if (tRangeSep) then
      !! Number of vir-vir transitions a->b _and_ b->a, summed over spin channels
      nxvv_a = sum(nvir_ud**2)
      allocate(qv(natom, max(sum(nxov_ud), nxvv_a)))
      qv(:,:) = 0.0_dp
      do jbs = 1, nmat
        jj = getIA(win(jbs), 1)
        bb = getIA(win(jbs), 2)
        ss = getIA(win(jbs), 3)
        do aa = nocc_ud(ss) + 1, norb
          abs = iaTrans(aa, bb, ss)
          qij(:) = transChrg%qTransAB(abs, denseDesc, ovrXev, grndEigVecs, getAB)
          jas = iaTrans(jj, aa, ss)
          qv(:,jas) = qv(:,jas) + qij * vTmp(jbs) * sqrOccIA(jbs)
        end do
      end do
      !! Attention: Necessary to diff. between nmat/nxov_rd and all possible exc
      do jas = 1, sum(nxov_ud)
        otmp(:) = qv(:,jas)
        call hemv(qv(:,jas), lrGamma, otmp)
      end do
      do ias = 1, nmat
        ii = getIA(win(ias), 1)
        aa = getIA(win(ias), 2)
        ss = getIA(win(ias), 3)
        do jj = 1, nocc_ud(ss)
          ijs = iaTrans(ii, jj, ss)
          jas = iaTrans(jj, aa, ss)
          qij(:) = transChrg%qTransIJ(ijs, denseDesc, ovrXev, grndEigVecs, getIJ)
          vout(ias) = vout(ias) - cExchange * sqrOccIA(ias) * dot_product(qij, qv(:, jas))
        end do
      end do

      qv(:,:) = 0.0_dp
      do jbs = 1, nmat
        jj = getIA(win(jbs), 1)
        bb = getIA(win(jbs), 2)
        ss = getIA(win(jbs), 3)
        !! Here, index abs is different for a,b and b,a.
        do aa = nocc_ud(ss) + 1, norb
          jas = iaTrans(jj, aa, ss)
          qij(:) = transChrg%qTransIA(jas, denseDesc, ovrXev, grndEigVecs, getIA, win)
          abs = (bb - nocc_ud(ss) - 1) * nvir_ud(ss) + aa - nocc_ud(ss)
          if (ss==2) then
            abs = abs + nvir_ud(1)**2
          end if
          !! qv is not symmetric in a and b
          qv(:,abs) = qv(:,abs) + qij * vTmp(jbs) * sqrOccIA(jbs)
        end do
      end do

      do abs = 1, nxvv_a
        otmp(:) = qv(:,abs)
        call hemv(qv(:,abs), lrGamma, otmp)
      end do

      do ias = 1, nmat
        ii = getIA(win(ias), 1)
        aa = getIA(win(ias), 2)
        ss = getIA(win(ias), 3)
        do bb = nocc_ud(ss) + 1, norb
          ibs = iaTrans(ii, bb, ss)
          qij(:) = transChrg%qTransIA(ibs, denseDesc, ovrXev, grndEigVecs, getIA, win)
          abs = (bb - nocc_ud(ss) - 1) * nvir_ud(ss) + aa - nocc_ud(ss)
          if (ss==2) then
            abs = abs + nvir_ud(1)**2
          end if
          vout(ias) = vout(ias) - cExchange * sqrOccIA(ias) * dot_product(qij, qv(:, abs))
       end do
      end do
    endif

    if (allocated(onsMEs)) then
      call onsiteEner(spin, sym, wij, sqrOccIA, win, nxov_ud(1), denseDesc%iAtomStart, getIA,&
          & species0, ovrXev, grndEigVecs, onsMEs, orb, vTmp, vout_ons)
      vout(:) = vout + vout_ons
    end if

    vout(:) = vout + wij * vTmp

    if(.not. tAplusB) then
       vout(:) = vout * sqrt(wij)
    endif

  end subroutine actionAplusB

    !> Multiplies the excitation supermatrix with a supervector.
  !> For the hermitian RPA eigenvalue problem this corresponds to \Omega_ias,jbt * v_jbt
  !> (spin polarized case) or \Omega^{S/T}_ia,jb * v_jb (singlet/triplet)
  !>
  !> For the standard RPA, (A+B)_ias,jbt * v_jbt needs to be computed (similar for singlet/triplet)
  !> (see definitions in Marc Casida, in Recent Advances in Density Functional Methods,
  !>  World Scientific, 1995, Part I, p. 155.)
  !> Note: we actually compute sqrt(n_is-n_as) (A+B)_ias,jbt sqrt(n_jt-n_bt), with the
  !> occupations n_is, correct for finite T.
  !> See also Dominguez JCTC 9 4901 (2013), Kranz JCTC 13 1737 (2017) for DFTB specifics.
  !>
  !> Note: In order not to store the entire supermatrix (nmat, nmat), the various pieces are
  !> assembled individually and multiplied directly with the corresponding part of the supervector.
  subroutine actionAplusB_MPI(comm, locSize, vOffset, spin, wij, sym, win, nocc_ud, nvir_ud,&
      & nxoo_ud, nxvv_ud, nxov_ud, nxov_rd, iaTrans, getIA, getIJ, getAB, env, denseDesc, ovrXev,&
      & grndEigVecs, occNr, sqrOccIA, gamma, species0, spinW, onsMEs, orb, tAplusB, transChrg, vin,&
      & vout, tRangeSep, lrGamma)

    !> MPI communicator
    integer, intent(in) :: comm

    !> Size of local vectors for each rank
    integer, intent(in) :: locSize(:)

    !> Vector offset for each rank
    integer, intent(in) :: vOffset(:)

    !> logical spin polarization
    logical, intent(in) :: spin

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

    !> Environment settings
    type(TEnvironment), intent(in) :: env

    !> Dense matrix descriptor
    type(TDenseDescr), intent(in) :: denseDesc 

    !> overlap times eigenvector. (nOrb, nOrb)
    real(dp), intent(in) :: ovrXev(:,:,:)

    !> eigenvectors (nOrb, nOrb)
    real(dp), intent(in) :: grndEigVecs(:,:,:)

    !> occupation numbers
    real(dp), intent(in) :: occNr(:,:)

    !> Square root of occupation difference between vir and occ states
    real(dp), intent(in) :: sqrOccIA(:)

    !> DFTB gamma matrix (nAtm, nAtom)
    real(dp), intent(in) :: gamma(:,:)

    !> chemical species of the atoms
    integer, intent(in) :: species0(:)

    !> ground state spin constants for each species
    real(dp), intent(in) :: spinW(:)

    !> onsite matrix elements for shells (elements between s orbitals on the same shell are ignored)
    real(dp), intent(in), allocatable :: onsMEs(:,:,:,:)

    !> data type for atomic orbital information
    type(TOrbitals), intent(in) :: orb

    !> whether (A+B)v or (Omega)v should be returned
    logical, intent(in) :: tAplusB

    !> machinery for transition charges between single particle levels
    type(TTransCharges), intent(in) :: transChrg

    !> vector to multiply with size(nmat)
    real(dp), intent(in) :: vin(:)

    !> vector containing the result on exit size(nmat)
    real(dp), intent(out) :: vout(:)

    !> is calculation range-separated?
    logical, intent(in) :: tRangeSep

    !> long-range Gamma if in use
    real(dp), allocatable, optional, intent(in) :: lrGamma(:,:)

    integer :: nmat, natom, norb, ia, nxvv_a
    integer :: ii, aa, ss, jj, bb, ias, jbs, abs, ijs, jas, ibs
    integer :: iam, nLoc, myia, ierr, iGlb, fGlb

    ! somewhat ugly, but fast small arrays on stack:
    real(dp) :: otmp(size(gamma, dim=1)), gtmp(size(gamma, dim=1))
    real(dp) :: qij(size(gamma, dim=1)) ! qij Working array (used for excitation charges). (nAtom)
    real(dp) :: vout_ons(size(vin))
    !!real(dp) :: vTmp(size(vin))
    real(dp), allocatable :: vTmp(:)
    real(dp), dimension(2) :: spinFactor = (/ 1.0_dp, -1.0_dp /)
    ! for later use to change HFX contribution
    real(dp), allocatable :: qv(:,:)
    external mpi_comm_rank, mpi_allreduce


    @:ASSERT(size(vin) == size(vout))

    !! dimension guessed from input vector
    nmat = size(vin)
    @:ASSERT(nmat <= nxov_rd)
    natom = size(gamma, dim=1)
    norb = nocc_ud(1)+nvir_ud(1)
    call mpi_comm_rank(comm, iam, ierr)
    nLoc = locSize(iam+1)
    @:ASSERT(nmat == nLoc)
    iGlb = vOffset(iam+1) + 1 
    fGlb = vOffset(iam+1) + nLoc
    allocate(vTmp, mold=vin)
    vTmp = 0.0_dp
    vout = 0.0_dp

    if(tAplusB) then
      vTmp(:) = vin
    else
      vTmp(:) = vin * sqrt(wij(iGlb:fGlb))
    endif

    vTmp(:) =  vTmp * sqrOccIA(iGlb:fGlb)
    oTmp(:) = 0.0_dp    
    do myia = 1, nLoc
      ia = vOffset(iam+1) + myia
      qij(:) = transChrg%qTransIA(ia, denseDesc, ovrXev, grndEigVecs, getIA, win)
      otmp(:) = otmp(:) + qij(:)*vTmp(myia)
    enddo

    call mpi_allreduce(MPI_IN_PLACE, otmp, natom, MPI_DOUBLE_PRECISION, MPI_SUM, comm, ierr)

    if (.not.spin) then !-----------spin-unpolarized systems--------------

      if (sym == 'S') then

        call hemv(gtmp, gamma, otmp)

        ! 4 * wn * (g * Q)
        vOut(:) = 0.0_dp
        do myia = 1, nLoc
          ia = vOffset(iam+1) + myia
          qij(:) = transChrg%qTransIA(ia, denseDesc, ovrXev, grndEigVecs, getIA, win)
          vOut(myia) = 4.0_dp * sqrOccIA(ia) * dot_product(qij, gTmp)
        enddo  

      else

        otmp = otmp * spinW(species0)

        ! 4 * wn * (o * Q)
        vOut = 0.0_dp
        call transChrg%qVecMat(denseDesc, ovrXev, grndEigVecs, getIA, win, oTmp, vOut, vOffset(iam+1))
        vOut(:) = 4.0_dp * sqrOccIA(iGlb:fGlb) * vOut


      end if

    else !--------------spin-polarized systems--------

      call hemv(gtmp, gamma, otmp)

      otmp(:) = 0.0_dp

      !$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(ia,ss,qij) &
      !$OMP& SCHEDULE(RUNTIME) REDUCTION(+:otmp)
      do myia = 1, nLoc
         
        ia = vOffset(iam+1) + myia
        qij(:) = transChrg%qTransIA(ia, denseDesc, ovrXev, grndEigVecs, getIA, win)

        ! singlet gamma part (S)
        vOut(myia) = 2.0_dp * sqrOccIA(ia) * dot_product(qij, gtmp)

        ! magnetization part (T1)
        ss = getIA(win(ia), 3)

        otmp(:) = otmp + spinFactor(ss) * sqrOccIA(ia) * vOut(myia) * qij

      end do
      !$OMP  END PARALLEL DO

      otmp = otmp * spinW(species0)

      !$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(ia,ss,qij) &
      !$OMP& SCHEDULE(RUNTIME)
      do myia = 1, nLoc

        ia = vOffset(iam+1) + myia        
        qij(:) = transChrg%qTransIA(ia, denseDesc, ovrXev, grndEigVecs, getIA, win)

        ss = getIA(win(ia), 3)

        vOut(myia) = vOut(myia) + 2.0_dp * sqrOccIA(ia) * spinFactor(ss) * dot_product(qij, otmp)

      end do
      !$OMP  END PARALLEL DO

    end if
    
    if (tRangeSep) then
      !! Number of vir-vir transitions a->b _and_ b->a, summed over spin channels
      nxvv_a = sum(nvir_ud**2)
      allocate(qv(natom, max(sum(nxov_ud), nxvv_a)))
      qv(:,:) = 0.0_dp
      do jbs = 1, nmat
        jj = getIA(win(jbs), 1)
        bb = getIA(win(jbs), 2)
        ss = getIA(win(jbs), 3)
        do aa = nocc_ud(ss) + 1, norb
          abs = iaTrans(aa, bb, ss)
          qij(:) = transChrg%qTransAB(abs, denseDesc, ovrXev, grndEigVecs, getAB)
          jas = iaTrans(jj, aa, ss)
          qv(:,jas) = qv(:,jas) + qij * vTmp(jbs) * sqrOccIA(jbs)
        end do
      end do
      !! Attention: Necessary to diff. between nmat/nxov_rd and all possible exc
      do jas = 1, sum(nxov_ud)
        otmp(:) = qv(:,jas)
        call hemv(qv(:,jas), lrGamma, otmp)
      end do
      do ias = 1, nmat
        ii = getIA(win(ias), 1)
        aa = getIA(win(ias), 2)
        ss = getIA(win(ias), 3)
        do jj = 1, nocc_ud(ss)
          ijs = iaTrans(ii, jj, ss)
          jas = iaTrans(jj, aa, ss)
          qij(:) = transChrg%qTransIJ(ijs, denseDesc, ovrXev, grndEigVecs, getIJ)
          vout(ias) = vout(ias) - cExchange * sqrOccIA(ias) * dot_product(qij, qv(:, jas))
        end do
      end do

      qv(:,:) = 0.0_dp
      do jbs = 1, nmat
        jj = getIA(win(jbs), 1)
        bb = getIA(win(jbs), 2)
        ss = getIA(win(jbs), 3)
        !! Here, index abs is different for a,b and b,a.
        do aa = nocc_ud(ss) + 1, norb
          jas = iaTrans(jj, aa, ss)
          qij(:) = transChrg%qTransIA(jas, denseDesc, ovrXev, grndEigVecs, getIA, win)
          abs = (bb - nocc_ud(ss) - 1) * nvir_ud(ss) + aa - nocc_ud(ss)
          if (ss==2) then
            abs = abs + nvir_ud(1)**2
          end if
          !! qv is not symmetric in a and b
          qv(:,abs) = qv(:,abs) + qij * vTmp(jbs) * sqrOccIA(jbs)
        end do
      end do

      do abs = 1, nxvv_a
        otmp(:) = qv(:,abs)
        call hemv(qv(:,abs), lrGamma, otmp)
      end do

      do ias = 1, nmat
        ii = getIA(win(ias), 1)
        aa = getIA(win(ias), 2)
        ss = getIA(win(ias), 3)
        do bb = nocc_ud(ss) + 1, norb
          ibs = iaTrans(ii, bb, ss)
          qij(:) = transChrg%qTransIA(ibs, denseDesc, ovrXev, grndEigVecs, getIA, win)
          abs = (bb - nocc_ud(ss) - 1) * nvir_ud(ss) + aa - nocc_ud(ss)
          if (ss==2) then
            abs = abs + nvir_ud(1)**2
          end if
          vout(ias) = vout(ias) - cExchange * sqrOccIA(ias) * dot_product(qij, qv(:, abs))
       end do
      end do
    endif

   #:if WITH_SCALAPACK
    
    if (allocated(onsMEs)) then
      call error('Onsite corrections not parallelized yet.')
    end if
    
   #:else
   
    if (allocated(onsMEs)) then
      call onsiteEner(spin, sym, wij, sqrOccIA, win, nxov_ud(1), denseDesc%iAtomStart, getIA, species0, &
          & ovrXev, grndEigVecs, onsMEs, orb, vTmp, vout_ons)      
      vout(:) = vout + vout_ons
    end if
     
   #:endif
    
    vOut(:) = vOut + wij(iGlb:fGlb) * vTmp
    
    if(.not. tAplusB) then
      vOut(:) = vOut * sqrt(wij(iGlb:fGlb))
    endif

  end subroutine actionAplusB_MPI


  !> Multiplies the excitation supermatrix with a supervector.
  !> (A-B)_ias,jbt * v_jbt is computed (and similar for singlet/triplet)
  !> (see definitions in Marc Casida, in Recent Advances in Density Functional Methods,
  !>  World Scientific, 1995, Part I, p. 155.)
  !> See also Dominguez JCTC 9 4901 (2013), Kranz JCTC 13 1737 (2017) for DFTB specifics.
  subroutine actionAminusB(spin, wij, win, nocc_ud, nvir_ud, nxoo_ud, nxvv_ud, nxov_ud, nxov_rd,&
      & iaTrans, getIA, getIJ, getAB, env, denseDesc, ovrXev, grndEigVecs, occNr, sqrOccIA, transChrg,&
      & vin, vout, tRangeSep, lrGamma)

    !> logical spin polarization
    logical, intent(in) :: spin

    !> excitation energies (wij = epsion_j - epsilon_i)
    real(dp), intent(in) :: wij(:)

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

    !> Environment settings
    type(TEnvironment), intent(in) :: env

    !> Dense matrix descriptor
    type(TDenseDescr), intent(in) :: denseDesc

    !> overlap times eigenvector. (nOrb, nOrb)
    real(dp), intent(in) :: ovrXev(:,:,:)

    !> eigenvectors (nOrb, nOrb)
    real(dp), intent(in) :: grndEigVecs(:,:,:)

    !> occupation numbers
    real(dp), intent(in) :: occNr(:,:)

    ! Square root of occupation difference between vir and occ states
    real(dp), intent(in) :: sqrOccIA(:)

    !> machinery for transition charges between single particle levels
    type(TTransCharges), intent(in) :: transChrg

    !> vector to multiply with size(nmat)
    real(dp), intent(in) :: vin(:)

    !> vector containing the result on exit size(nmat)
    real(dp), intent(out) :: vout(:)

    !> is calculation range-separated?
    logical, intent(in) :: tRangeSep

    !> long-range Gamma if in use
    real(dp), allocatable, optional, intent(in) :: lrGamma(:,:)

    integer :: nmat, natom, norb, nxvv_a
    integer :: ii, aa, ss, jj, bb, ias, jbs, abs, ijs, jas, ibs
    real(dp), allocatable :: qv(:,:)
    real(dp), allocatable :: qij(:), otmp(:)

    nmat = size(vin)
    norb = size(ovrXev, dim=1)
    vout(:) = 0.0_dp

    if (tRangeSep) then
      natom = size(lrGamma, dim=1)
      !! Number of vir-vir transitions a->b _and_ b->a, summed over spin channels
      nxvv_a = sum(nvir_ud**2)
      allocate(qij(natom))
      allocate(otmp(natom))
      allocate(qv(natom, max(sum(nxov_ud),nxvv_a)))

      qv(:,:) = 0.0_dp
      do jbs = 1, nmat
        jj = getIA(win(jbs), 1)
        bb = getIA(win(jbs), 2)
        ss = getIA(win(jbs), 3)
        do aa = nocc_ud(ss) + 1, norb
          abs = iaTrans(aa, bb, ss)
          qij(:) = transChrg%qTransAB(abs, denseDesc, ovrXev, grndEigVecs, getAB)
          jas = iaTrans(jj, aa, ss)
          qv(:,jas) = qv(:,jas) + qij * vin(jbs) * sqrOccIA(jbs)
        end do
      end do

      !! Attention: Necessary to diff. between nmat/nxov_rd and all possible exc
      do jas = 1, sum(nxov_ud)
        otmp(:) = qv(:,jas)
        call hemv(qv(:,jas), lrGamma, otmp)
      end do

      do ias = 1, nmat
        ii = getIA(win(ias), 1)
        aa = getIA(win(ias), 2)
        ss = getIA(win(ias), 3)
        do jj = 1, nocc_ud(ss)
          ijs = iaTrans(ii, jj, ss)
          jas = iaTrans(jj, aa, ss)
          qij(:) = transChrg%qTransIJ(ijs, denseDesc, ovrXev, grndEigVecs, getIJ)
          vout(ias) = vout(ias) - cExchange * sqrOccIA(ias) * dot_product(qij, qv(:, jas))
        end do
      end do

      qv(:,:) = 0.0_dp
      do jbs = 1, nmat
        jj = getIA(win(jbs), 1)
        bb = getIA(win(jbs), 2)
        ss = getIA(win(jbs), 3)
        !! Here, index abs is different for a,b and b,a.
        do aa = nocc_ud(ss) + 1, norb
          jas = iaTrans(jj, aa, ss)
          qij(:) = transChrg%qTransIA(jas, denseDesc, ovrXev, grndEigVecs, getIA, win)
          abs = (bb - nocc_ud(ss) - 1) * nvir_ud(ss) + aa - nocc_ud(ss)
          if (ss==2) then
            abs = abs + nvir_ud(1)**2
          end if
          !! qv is not symmetric in a and b
          qv(:,abs) = qv(:,abs) + qij * vin(jbs) * sqrOccIA(jbs)
        end do
      end do

      do abs = 1, nxvv_a
        otmp(:) = qv(:,abs)
        call hemv(qv(:,abs), lrGamma, otmp)
      end do

      do ias = 1, nmat
        ii = getIA(win(ias), 1)
        aa = getIA(win(ias), 2)
        ss = getIA(win(ias), 3)
        do bb = nocc_ud(ss) + 1, norb
          ibs = iaTrans(ii, bb, ss)
          qij(:) = transChrg%qTransIA(ibs, denseDesc, ovrXev, grndEigVecs, getIA, win)
          abs = (bb - nocc_ud(ss) - 1) * nvir_ud(ss) + aa - nocc_ud(ss)
          if (ss==2) then
            abs = abs + nvir_ud(1)**2
          end if
          vout(ias) = vout(ias) + cExchange * sqrOccIA(ias) * dot_product(qij, qv(:, abs))
       end do
      end do
    endif

    ! orb. energy difference diagonal contribution
    vout(:) = vout + wij * vin

  end subroutine actionAminusB


  !> Generates initial matrices M+ and M- for the RPA algorithm by Stratmann
  !> (JCP 109 8218 (1998).
  !> M+/- = (A+/-B)_ias,jbt (spin polarized) (A+/-B)^{S/T}_ia,jb (closed shell)
  !> Here ias,jbt <= initDim
  !> Also computed is v+/- = (A+/-B)_ias,jbt with ias <= nMat, jbt <= initDim
  !> Note: Routine not set up to handle onsite corrections.
  !> Note: Not yet OpenMP parallelized
  subroutine initialSubSpaceMatrixApmB(transChrg, initDim, wIJ, sym, win, nmatup, env, denseDesc,&
      & sTimesGrndEigVecs, grndEigVecs, occNr, sqrOccIA, getIA, getIJ, getAB, iaTrans, gamma,&
      & lrGamma, species0, spinW, tSpin, tRangeSep, vP, vM, mP, mM)

    !> machinery for transition charges between single particle levels
    type(TTransCharges), intent(in) :: transChrg

    !> initial RPA subspace
    integer, intent(in) :: initDim

    !> number of same-spin transitions
    integer, intent(in) :: nmatup

    !> excitation energies (wij = epsion_j - epsilon_i)
    real(dp), intent(in) :: wIJ(:)

    !> occupation numbers
    real(dp), intent(in) :: occNr(:,:)

    ! Square root of occupation difference between vir and occ states
    real(dp), intent(in) :: sqrOccIA(:)

    !> symmetry flag (singlet or triplet)
    character, intent(in) :: sym

    !> sorting index of the excitation energies.
    integer, intent(in) :: win(:)

    !> Environment settings
    type(TEnvironment), intent(in) :: env

    !> Dense matrix descriptor
    type(TDenseDescr), intent(in) :: denseDesc 

    !> overlap times eigenvector. (nOrb, nOrb)
    real(dp), intent(in) :: sTimesGrndEigVecs(:,:,:)

    !> eigenvectors (nOrb, nOrb)
    real(dp), intent(in) :: grndEigVecs(:,:,:)

    !> DFTB gamma matrix (nAtm, nAtom)
    real(dp), intent(in) :: gamma(:,:)

    !> DFTB gamma matrix (long-range corrected)
    real(dp), intent(in), allocatable :: lrGamma(:,:)

    !> index array for excitations
    integer, intent(in) :: getIA(:,:)

    !> index array for occ-occ single particle excitations
    integer, intent(in) :: getIJ(:,:)

    !> index array for vir-vir single particle excitations
    integer, intent(in) :: getAB(:,:)

    !> index array from orbital pairs to compound index
    integer, intent(in) :: iaTrans(:,:,:)

    !> chemical species of the atoms
    integer, intent(in) :: species0(:)

    !> ground state spin constants for each species
    real(dp), intent(in) :: spinW(:)

    !> logical spin polarization
    logical, intent(in) :: tSpin

    !> is calculation range-separated?
    logical, intent(in) :: tRangeSep

    !> output vector v+ (nMat, initDim)
    real(dp), intent(out) :: vP(:,:)

    !> output vector v+ (nMat, initDim)
    real(dp), intent(out) :: vM(:,:)

    !> output matrix M+ (initDim, initDim)
    real(dp), intent(out) :: mP(:,:)

    !> output matrix M- (initDim, initDim)
    real(dp), intent(out) :: mM(:,:)

    integer :: nMat, nAtom
    integer :: ia, jb, ii, jj, ss, tt
    real(dp), allocatable :: oTmp(:), gTmp(:), qTr(:)
    real(dp), dimension(2) :: spinFactor = (/ 1.0_dp, -1.0_dp /)
    real(dp) :: rTmp
    integer :: aa, bb, iat, jbs, abs, ijs, ibs, jas

    nMat = size(vP, dim=1) ! also known as nXov
    nAtom = size(gamma, dim=1)

    allocate(gTmp(nAtom))
    allocate(oTmp(nAtom))
    allocate(qTr(nAtom))

    vP(:,:) = 0.0_dp
    vM(:,:) = 0.0_dp
    mP(:,:) = 0.0_dp
    mM(:,:) = 0.0_dp

    oTmp(:) = 0.0_dp

    !-----------spin-unpolarized systems--------------
    if(.not. tSpin) then

      if (sym == 'S') then

        ! full range coupling matrix contribution: 4 * sum_A q^ia_A sum_B gamma_AB q^jb_B
        do jb = 1, initDim
          qTr(:) = transChrg%qTransIA(jb, denseDesc, sTimesGrndEigVecs, grndEigVecs, getIA, win)
          qTr(:) = qTr * sqrOccIA(jb)

          call hemv(oTmp, gamma, qTr)

          do ia = 1, nMat
            qTr(:) = transChrg%qTransIA(ia, denseDesc, sTimesGrndEigVecs, grndEigVecs, getIA, win)
            vP(ia,jb) = 4.0_dp * sqrOccIA(ia) * dot_product(qTr, oTmp)
          end do

        end do

      else

        ! full range triplets contribution: 2 * sum_A q^ia_A M_A q^jb_A
        do jb = 1, initDim
          qTr(:) = transChrg%qTransIA(jb, denseDesc, sTimesGrndEigVecs, grndEigVecs, getIA, win)
          oTmp(:) = qTr * sqrOccIA(jb) * spinW(species0)

          do ia = 1, nMat
            qTr(:) = transChrg%qTransIA(ia, denseDesc, sTimesGrndEigVecs, grndEigVecs, getIA, win)
            vP(ia,jb) = vP(ia,jb) + 4.0_dp * sqrOccIA(ia) * dot_product(qTr, oTmp)
          end do

        end do

      end if
    !--------------spin-polarized systems--------
    else

      do jb = 1, initDim
        qTr(:) = transChrg%qTransIA(jb, denseDesc, sTimesGrndEigVecs, grndEigVecs, getIA, win)
        qTr(:) = qTr * sqrOccIA(jb)

        call hemv(gTmp, gamma, qTr)

        oTmp(:) = qTr

        do ia = 1, nMat
          qTr(:) = transChrg%qTransIA(ia, denseDesc, sTimesGrndEigVecs, grndEigVecs, getIA, win)
          vP(ia,jb) = 2.0_dp * sqrOccIA(ia) * dot_product(qTr, gTmp)
        end do

        ss = getIA(win(jb), 3)

        oTmp(:)  =  spinFactor(ss) * 2.0_dp * spinW(species0) * oTmp

        do ia = 1, nMat
          qTr(:) = transChrg%qTransIA(ia, denseDesc, sTimesGrndEigVecs, grndEigVecs, getIA, win)
          tt = getIA(win(ia), 3)
          vP(ia,jb) = vP(ia,jb) + spinFactor(tt) * sqrOccIA(ia) * dot_product(qTr, oTmp)
        end do

      end do

    end if

    if (tRangeSep) then

      do jbs = 1, initDim
        jj = getIA(win(jbs), 1)
        bb = getIA(win(jbs), 2)
        ss = getIA(win(jbs), 3)

        do iat = 1, nMat
          ii = getIA(win(iat), 1)
          aa = getIA(win(iat), 2)
          tt = getIA(win(iat), 3)

          if (ss /= tt) cycle

          abs = iaTrans(aa, bb, ss)
          qTr(:) = transChrg%qTransAB(abs, denseDesc, sTimesGrndEigVecs, grndEigVecs, getAB)
          oTmp(:) = 0.0_dp
          call hemv(oTmp, lrGamma, qTr)

          ijs = iaTrans(ii, jj, ss)
          qTr(:) = transChrg%qTransIJ(ijs, denseDesc, sTimesGrndEigVecs, grndEigVecs, getIJ)
          rTmp = cExchange * sqrOccIA(iat) * sqrOccIA(jbs) * dot_product(qTr, oTmp)
          vP(iat,jbs) = vP(iat,jbs) - rTmp
          vM(iat,jbs) = vM(iat,jbs) - rTmp

          ibs = iaTrans(ii, bb, ss)
          qTr(:) = transChrg%qTransIA(ibs, denseDesc, sTimesGrndEigVecs, grndEigVecs, getIA, win)
          oTmp(:) = 0.0_dp
          call hemv(oTmp, lrGamma, qTr)

          jas = iaTrans(jj, aa, ss)
          qTr(:) = transChrg%qTransIA(jas, denseDesc, sTimesGrndEigVecs, grndEigVecs, getIA, win)
          rTmp = cExchange * sqrOccIA(iat) * sqrOccIA(jbs) * dot_product(qTr, oTmp)
          vP(iat,jbs) = vP(iat,jbs) - rTmp
          vM(iat,jbs) = vM(iat,jbs) + rTmp
        end do

      end do
    endif

    do jb = 1, initDim
      vP(jb,jb) = vP(jb,jb) + wIJ(jb)
      vM(jb,jb) = vM(jb,jb) + wIJ(jb)
    end do

    do ii = 1, initDim
      do jj = ii, initDim
        mP(ii,jj) = vP(ii,jj)
        mP(jj,ii) = mP(ii,jj)
        mM(ii,jj) = vM(ii,jj)
        mM(jj,ii) = mM(ii,jj)
      end do
    end do

  end subroutine initialSubSpaceMatrixApmB


  !> Onsite energy corrections
  subroutine onsiteEner(spin, sym, wij, sqrOccIA, win, nmatup, iAtomStart, getIA, species0, ovrXev,&
      & grndEigVecs, ons_en, orb, vin, vout)

    !> logical spin polarization
    logical, intent(in) :: spin

    !> symmetry flag (singlet or triplet)
    character, intent(in) :: sym

    !> excitation energies (wij = epsion_j - epsilon_i)
    real(dp), intent(in) :: wij(:)

    !> Square root of occupation difference between vir and occ states
    real(dp), intent(in) :: sqrOccIA(:)

    !> sorting index of the excitation energies.
    integer, intent(in) :: win(:)

    !> number of same-spin transitions
    integer, intent(in) :: nmatup

    !> starting position of each atom in the list of orbitals.
    integer, intent(in) :: iAtomStart(:)

    !> index array for excitations
    integer, intent(in) :: getIA(:,:)

    !> chemical species of the atoms
    integer, intent(in) :: species0(:)

    !> overlap times ground state wavefunctions
    real(dp), intent(in) :: ovrXev(:,:,:)

    !> ground state wavefunctions
    real(dp), intent(in) :: grndEigVecs(:,:,:)

    !> onsite matrix elements for shells (elements between s orbitals on the same shell are ignored)
    real(dp), intent(in) :: ons_en(:,:,:,:)

    !> data type for atomic orbital information
    type(TOrbitals), intent(in) :: orb

    !> vector to multiply with size(nmat)
    real(dp), intent(in) :: vin(:)

    !> vector containing the result on exit size(nmat)
    real(dp), intent(out) :: vout(:)

    real(dp) :: otmp(orb%mOrb,orb%mOrb,size(species0),2)
    real(dp) :: fact
    real(dp) :: onsite(orb%mOrb,orb%mOrb,2)
    real(dp) :: qq_ij(orb%mOrb,orb%mOrb)
    logical :: updwn
    integer :: nmat, nAtom, nOrb
    integer :: ia, iAt, iSp, iSh, iOrb, ii, jj, ss, sindx(2), iSpin, nSpin
    real(dp) :: degeneracy, partTrace

    nmat = size(vin)
    nAtom = size(species0)
    nSpin = size(ovrXev, dim=3)

    vout(:) = 0.0_dp
    otmp(:,:,:,:)  = 0.0_dp

    fact = 1.0_dp
    if (.not.Spin) then
      if (sym == 'T') then
        fact = -1.0_dp
      end if
    end if

    do iAt = 1, nAtom
      iSp = species0(iAt)
      nOrb = orb%nOrbAtom(iAt)
      call getOnsME(orb, iSp, ons_en, nOrb, onsite)
      do ia = 1, nmat
        call indxov(win, ia, getIA, ii, jj, ss)
        updwn = (win(ia) <= nmatup)
        call transDens(ii, jj, iAt, iAtomStart, nOrb, updwn, ovrXev, grndEigVecs, qq_ij)
        if (spin) then
          if (.not. updwn) then
            sindx = [2, 1]
          else
            sindx = [1, 2]
          end if
          do iSpin = 1, 2
            otmp(:nOrb, :nOrb, iAt, iSpin) = otmp(:nOrb, :nOrb, iAt, iSpin) + qq_ij(:nOrb, :nOrb)&
                & * sqrOccIA(ia) * vin(ia) * onsite(:nOrb, :nOrb, sindx(iSpin))
          end do
        else
          ! closed shell
          otmp(:nOrb, :nOrb, iAt, 1) = otmp(:nOrb, :nOrb, iAt, 1) + &
              & qq_ij(:nOrb, :nOrb) * sqrOccIA(ia) * vin(ia)&
              & * (onsite(:nOrb, :nOrb, 1) + fact * onsite(:nOrb, :nOrb, 2))
        end if
      end do
      ! rotational invariance corection for diagonal part
      do iSpin = 1, nSpin
        do iSh = 1, orb%nShell(iSp)
          degeneracy = real(2*orb%angShell(iSh, iSp) + 1, dp)
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

    do ia = 1, nmat
      call indxov(win, ia, getIA, ii, jj, ss)
      updwn = (win(ia) <= nmatup)
      do iAt = 1, nAtom
        nOrb = orb%nOrbAtom(iAt)
        call transDens(ii, jj, iAt, iAtomStart, nOrb, updwn, ovrXev, grndEigVecs, qq_ij)
        vout(ia) = vout(ia) + 4.0_dp * sqrOccIA(ia) * &
            & sum(qq_ij(:nOrb, :nOrb) * otmp(:nOrb, :nOrb, iAt, ss))
      end do
    end do

  end subroutine onsiteEner


  !> calculating spin polarized excitations
  !> Note: the subroutine is generalized to account for spin and partial occupancy
  subroutine getSPExcitations(nOcc, nVir, grndEigVal, filling, wij, getIA, getIJ, getAB)

    !> number of occupied states per spin channel
    integer, intent(in) :: nOcc(:)

    !> number of virtual states per spin channel
    integer, intent(in) :: nVir(:)

    !> ground state eigenvalues
    real(dp), intent(in) :: grndEigVal(:,:)

    !> occupation numbers for the ground state
    real(dp), intent(in) :: filling(:,:)

    !> Kohn-Sham energy differences between empty and filled states
    real(dp), intent(out) :: wij(:)

    !> index of occ-vir pairs of KS states for the transitions in wij
    integer, intent(out) :: getIA(:,:)

    !> index of occ-occ pairs of KS states for the transitions in wij
    integer, intent(out) :: getIJ(:,:)

    !> index of vir-vir pairs of KS states for the transitions in wij
    integer, intent(out) :: getAB(:,:)

    integer :: ind, ii, jj, aa, bb, off
    integer :: norb, iSpin, nSpin

    @:ASSERT(all(shape(grndEigVal)==shape(filling)))

    norb = size(grndEigVal, dim=1)
    nSpin = size(grndEigVal, dim=2)

    ind = 0

    do iSpin = 1, nSpin
      do ii = 1, norb - 1
        do aa = ii, norb
          if (filling(ii,iSpin) > filling(aa,iSpin) + elecTolMax) then
            ind = ind + 1
            wij(ind) = grndEigVal(aa,iSpin) - grndEigVal(ii,iSpin)
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


  !> Calculate dipole moment for transitions
  function oscillatorStrength(tSpin, transd, omega, xpy, sqrOccIA) result(osz)

    !> spin-polarized calculation ?
    logical, intent(in) :: tSpin

    !> transition dipole matrix elements for each atom
    real(dp), intent(in) :: transd(:,:)

    !> excitation energy
    real(dp), intent(in) :: omega

    !> excitation eigenvector (X+Y) for state in question
    real(dp), intent(in) :: xpy(:)

    !> square root of KS occupation differences
    real(dp), intent(in) :: sqrOccIA(:)

    !> oscillator strength
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


  !> Transition dipoles between single particle states
  subroutine transitionDipole(tSpin, transd, xpy, sqrOccIA, tdip)

    !> spin-polarized calculation ?
    logical, intent(in) :: tSpin

    !> transition dipole for transitions between KS states
    real(dp), intent(in) :: transd(:,:)

    !> eigenvectors RPA equations (X+Y)
    real(dp), intent(in) :: xpy(:,:)

    !> square root of KS occupation differences
    real(dp), intent(in) :: sqrOccIA(:)

    !> resulting transition dipoles
    real(dp), intent(out) :: tdip(:,:)

    integer :: ii, ll

    tdip(:,:) = 0.0_dp
    !$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(ii) SCHEDULE(RUNTIME)
    do ii = 1, size(xpy, dim=2)
      do ll = 1, 3
        tdip(ii,ll) = sum(transd(:,ll) * sqrOccIA * xpy(:,ii))
      end do
    end do
    !$OMP  END PARALLEL DO

    ! For spin-unpolarized systems
    ! (X+Y)_ia_up = (X+Y)_ia_dn = (X+Y)_ia^Singlet / sqrt(2)
    if (.not. tSpin) then
      tdip(:,:) = tdip * sqrt(2.0_dp)
    end if

  end subroutine transitionDipole


  !> Calculate transition moments for transitions between Kohn-Sham states, including spin-flipping
  !> transitions
  subroutine calcTransitionDipoles(coord0, win, nmatup, getIA, transChrg, env, denseDesc, ovrXev, &
      & grndEigVecs, snglPartTransDip)

    !> Atomic positions
    real(dp), intent(in) :: coord0(:,:)

    !> single particle transition index
    integer, intent(in) :: win(:)

    !> number of same-spin transitions
    integer, intent(in) :: nmatup

    !> Environment settings
    type(TEnvironment), intent(in) :: env

    !> Dense matrix descriptor
    type(TDenseDescr), intent(in) :: denseDesc 

    !> index array for excitation pairs
    integer, intent(in) :: getIA(:,:)

    !> machinery for transition charges between single particle levels
    type(TTransCharges), intent(in) :: transChrg

    !> overlap times ground state wavefunctions
    real(dp), intent(in) :: ovrXev(:,:,:)

    !> ground state wavefunctions
    real(dp), intent(in) :: grndEigVecs(:,:,:)

    !> resulting transition dipoles
    real(dp), intent(out) :: snglPartTransDip(:,:)

    integer :: nxov, natom, indm
    real(dp), allocatable :: qij(:)

    nxov = size(win)
    natom = size(coord0, dim=2)

    allocate(qij(natom))

    ! Calculate transition dipole elements
    do indm = 1, nxov
      qij(:) = transChrg%qTransIA(indm, denseDesc, ovrXev, grndEigVecs, getIA, win) 
      snglPartTransDip(indm, :) = matmul(coord0, qij)
    end do

  end subroutine calcTransitionDipoles


  !> Calculate <S^2> as a measure of spin contamination (smaller magnitudes are better, 0.5 is
  !> considered an upper threshold for reliability according to Garcia thesis)
  subroutine getExcSpin(Ssq, nmatup, getIA, win, eval, xpy, filling, ovrXev, grndEigVecs)

    !> spin contamination
    real(dp), intent(out) :: Ssq(:)

    !> number of spin up excitations
    integer, intent(in) :: nmatup

    !> index for composite excitations to specific occupied and empty states
    integer, intent(in) :: getIA(:,:)

    !> single particle excitation index
    integer, intent(in) :: win(:)

    !> Casida exitation energies
    real(dp), intent(in) :: eval(:)

    !> Casida excited eigenvectors (X+Y)
    real(dp), intent(in) :: xpy(:,:)

    !> occupations in ground state
    real(dp), intent(in) :: filling(:,:)

    !> Overlap times ground state eigenvectors
    real(dp), intent(in) :: ovrXev(:,:,:)

    !> Ground state eigenvectors
    real(dp), intent(in) :: grndEigVecs(:,:,:)

    integer:: i, k, l, m, ia, jb, ii, aa, jj, bb, ss
    integer:: nmat, nexc, nup, ndwn
    real(dp) :: TDvnorm
    real(dp), allocatable :: TDvec(:), TDvec_sq(:)
    integer, allocatable :: TDvin(:)
    logical :: ud_ia, ud_jb
    real(dp) :: s_iaja, s_iaib, s_iajb, tmp

    nmat = size(xpy, dim=1)
    nexc = size(Ssq)
    nup = ceiling(sum(filling(:,1)))
    ndwn = ceiling(sum(filling(:,2)))
    allocate(TDvec(nmat))
    allocate(TDvec_sq(nmat))
    allocate(TDvin(nmat))

    do i = 1, nexc
      TDvec(:) = xpy(:,i)
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
        call indxov(win, ia, getIA, ii, aa, ss)
        ud_ia = (win(ia) <= nmatup)
        do l = 1, nmat
          jb = TDvin(l)
          call indxov(win, jb, getIA, jj, bb, ss)
          ud_jb = (win(jb) <= nmatup)

          if ( (bb /= aa) .or. (ud_jb .neqv. ud_ia) ) then
            cycle
          end if

          tmp = 0.0_dp
          if (ud_ia) then
            do m = 1,ndwn
              tmp = tmp + MOoverlap(ii,m,ovrXev,grndEigVecs) * MOoverlap(jj,m,ovrXev,grndEigVecs)
            end do
          else
            do m = 1,nup
              tmp = tmp + MOoverlap(m,ii,ovrXev,grndEigVecs) * MOoverlap(m,jj,ovrXev,grndEigVecs)
            end do
          end if

          s_iaja = s_iaja + TDvec(ia)*TDvec(jb)*tmp

        end do
      end do

      ! S_{ia,ib}
      s_iaib = 0.0_dp
      do k = 1, nmat
        ia = TDvin(k)
        call indxov(win, ia, getIA, ii, aa, ss)
        ud_ia = (win(ia) <= nmatup)
        do l = 1, nmat
          jb = TDvin(l)
          call indxov(win, jb, getIA, jj, bb, ss)
          ud_jb = (win(jb) <= nmatup)

          if ( (ii /= jj) .or. (ud_jb .neqv. ud_ia) ) then
            cycle
          end if

          tmp = 0.0_dp
          if (ud_ia) then
            do m = 1,ndwn
              tmp = tmp + MOoverlap(aa,m,ovrXev,grndEigVecs) * MOoverlap(bb,m,ovrXev,grndEigVecs)
            end do
          else
            do m = 1,nup
              tmp = tmp + MOoverlap(m,aa,ovrXev,grndEigVecs) * MOoverlap(m,bb,ovrXev,grndEigVecs)
            end do
          end if

         s_iaib = s_iaib + TDvec(ia)*TDvec(jb)*tmp
        end do
      end do

      ! S_{ia,jb}
      s_iajb = 0.0_dp
      do k = 1, nmat
        ia = TDvin(k)
        call indxov(win, ia, getIA, ii, aa, ss)
        ud_ia = (win(ia) <= nmatup)
        if (.not. ud_ia ) then
          cycle
        end if
        do l = 1, nmat
          jb = TDvin(l)
          call indxov(win, jb, getIA, jj, bb, ss)
          ud_jb = (win(jb) <= nmatup)

          if ( ud_jb ) cycle

          s_iajb = s_iajb + TDvec(ia)*TDvec(jb) * MOoverlap(aa,bb,ovrXev,grndEigVecs)&
              & * MOoverlap(ii,jj,ovrXev,grndEigVecs)

        end do
      end do

      Ssq(i) =  s_iaja - s_iaib - 2.0_dp*s_iajb

    end do

  end subroutine getExcSpin


  !> Write single particle excitations to a file
  subroutine writeSPExcitations(wij, win, nmatup, getIA, writeSPTrans, sposz, nxov, tSpin)

    !> single particle excitation energies
    real(dp), intent(in) :: wij(:)

    !> index array for single particle transitions
    integer, intent(in) :: win(:)

    !> number of transitions within same spin channel
    integer, intent(in) :: nmatup

    !> index from composite index to occupied and virtual single particle states
    integer, intent(in) :: getIA(:,:)

    !> whether single particle excitation data should be written
    logical, intent(in) :: writeSPTrans

    !> single particle oscilation strengths
    real(dp), intent(in) :: sposz(:)

    !> Number of included single particle excitations to print out (assumes that win and wij are
    !> sorted so that the wanted transitions are first in the array)
    integer, intent(in) :: nxov

    !> is this a spin-polarized calculation?
    logical, intent(in) :: tSpin

    integer :: indm, m, n, s
    logical :: updwn
    character :: sign
    type(TFileDescr) :: fdSPTrans

    @:ASSERT(size(sposz)>=nxov)

    if (writeSPTrans) then
      ! single particle excitations
      call openFile(fdSPTrans, singlePartOut, mode="w")
      write(fdSPTrans%unit,*)
      write(fdSPTrans%unit,'(7x,a,7x,a,8x,a)') '#      w [eV]',&
          & 'Osc.Str.', 'Transition'
      write(fdSPTrans%unit,*)
      write(fdSPTrans%unit,'(1x,58("="))')
      write(fdSPTrans%unit,*)
      do indm = 1, nxov
        call indxov(win, indm, getIA, m, n, s)
        sign = " "
        if (tSpin) then
          updwn = (win(indm) <= nmatup)
          if (updwn) then
            sign = "U"
          else
            sign = "D"
          end if
        end if
        write(fdSPTrans%unit,&
            & '(1x,i7,3x,f8.3,3x,f13.7,4x,i5,3x,a,1x,i5,1x,1a)')&
            & indm, Hartree__eV * wij(indm), sposz(indm), m, '->', n, sign
      end do
      write(fdSPTrans%unit,*)
      call closeFile(fdSPTrans)
    end if

  end subroutine writeSPExcitations


  !> Excited state Mulliken charges and dipole moments written to disc
  subroutine writeExcMulliken(sym, nstat, dq, dqex, coord0)

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

    type(TFileDescr) :: fdMulliken
    integer :: natom, m
    real(dp) :: dipol(3), dipabs

    natom = size(dq)

    @:ASSERT(size(dq) == size(dqex))
    @:ASSERT(all(shape(coord0) == [3,nAtom]))

    ! Output of excited state Mulliken charges
    call openFile(fdMulliken, excitedQOut, mode="a")
    write(fdMulliken%unit, "(a,a,i2)") "# MULLIKEN CHARGES of excited state ", sym, nstat
    write(fdMulliken%unit, "(a,2x,A,i4)") "#", 'Natoms =',natom
    write(fdMulliken%unit, "('#',1X,A4,T15,A)")'Atom','netCharge'
    write(fdMulliken%unit,'("#",41("="))')
    do m = 1,  natom
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


  !> increase dimension of vector from sizeIn to fac*sizeIn
  pure subroutine incSizeVec(sizeIn, fac, vec)

    !> Size of the input vector to copy over to resized vector
    integer, intent(in) :: sizeIn

    !> Increment factor
    integer, intent(in) :: fac

    !> Vector to re-size, retaining initial elements in output
    real(dp), allocatable, intent(inout) :: vec(:)

    real(dp), allocatable :: temp(:)

    allocate(temp(fac * sizeIn))
    temp(:) = 0.0_dp
    temp(1:sizeIn) = vec
    call move_alloc(temp, vec)

  end subroutine incSizeVec


  !> increase size of (sizeIn,n) array to (fac*sizeIn,n)
  pure subroutine incSizeMatDimOne(sizeIn, fac, mat)

    !> Size of the input matrix first dimension to copy over to resized matrix
    integer, intent(inout) :: sizeIn

    !> Increment factor
    integer, intent(in) :: fac

    !> Matrix to re-size, retaining initial elements in output
    real(dp), allocatable, intent(inout) :: mat(:,:)

    integer :: dim2
    real(dp), allocatable :: temp(:,:)

    dim2 = size(mat, dim=2)
    allocate(temp(fac * sizeIn, dim2))
    temp(:,:) = 0.0_dp
    temp(1:sizeIn,:) = mat
    call move_alloc(temp, mat)

  end subroutine incSizeMatDimOne


  !> increase size of (n,sizeIn) array to (n, fac*sizeIn)
  pure subroutine incSizeMatDimTwo(sizeIn, fac, mat)

    !> Size of the input matrix second dimension to copy over to resized matrix
    integer, intent(inout) :: sizeIn

    !> Increment factor
    integer, intent(in) :: fac

    !> Matrix to re-size, retaining initial elements in output
    real(dp), allocatable, intent(inout) :: mat(:,:)

    integer :: dim1
    real(dp), allocatable :: temp(:,:)

    dim1 = size(mat, dim=1)
    allocate(temp(dim1, fac * sizeIn))
    temp(:,:) = 0.0_dp
    temp(:,1:sizeIn) = mat
    call move_alloc(temp, mat)

  end subroutine incSizeMatDimTwo


  !> increase size of (sizeIn,sizeIn) array to (fac1*sizeIn,fac2*sizeIn)
  pure subroutine incSizeMatBothDim(sizeIn, fac1, fac2, mat)

    !> Size of the input matrix second dimension to copy over to resized matrix (square)
    integer, intent(inout) :: sizeIn

    !> Increment factor for first dimension
    integer, intent(in) :: fac1

    !> Increment factor for second dimension
    integer, intent(in) :: fac2

    !> Matrix to re-size, retaining initial elements in output
    real(dp), allocatable, intent(inout) :: mat(:,:)

    real(dp), allocatable :: temp(:,:)

    allocate(temp(fac1 * sizeIn, fac2 * sizeIn))
    temp(:,:) = 0.0_dp
    temp(1:sizeIn, 1:sizeIn) = mat
    call move_alloc(temp, mat)

  end subroutine incSizeMatBothDim


  !> Calculate square root and inverse of sqrt of a real, symmetric positive definite matrix
  subroutine calcMatrixSqrt(matIn, spaceDim, memDim, workArray, workDim, matOut, matInvOut)

    !> Matrix to operate on
    real(dp), intent(in) :: matIn(:,:)

    !> Dimensions of input matrix
    integer, intent(in) :: spaceDim

    !> Leading dimension of resulting matrices
    integer, intent(in) :: memDim

    !> Matrix square root
    real(dp), intent(out) :: matOut(:,:)

    !> Inverse of matrix square root
    real(dp), intent(out) :: matInvOut(:,:)

    !> workspace array
    real(dp), intent(out) :: workArray(:)

    !> Size of work array
    integer :: workDim

    real(dp) :: dummyEV(spaceDim)
    real(dp) :: dummyM(spaceDim, spaceDim), dummyM2(spaceDim, spaceDim)
    integer :: info
    integer :: ii
    external dsyev, dgemm

    dummyM(:,:) = matIn(1:spaceDim,1:spaceDim)

    call dsyev('V', 'U', spaceDim, dummyM, spaceDim, dummyEV, workArray, workDim, info)

    !calc. sqrt
    do ii = 1, spaceDim
      dummyM2(:,ii) = sqrt(dummyEV(ii)) * dummyM(:,ii)
    end do

    call dgemm('N', 'T', spaceDim, spaceDim, spaceDim, 1.0_dp, dummyM2, spaceDim, dummyM,&
        & spaceDim, 0.0_dp, matOut, memDim)

    !calc. inv. of sqrt
    do ii = 1, spaceDim
      dummyM2(:,ii) = dummyM(:,ii) / sqrt(dummyEV(ii))
    end do

    call dgemm('N', 'T', spaceDim, spaceDim, spaceDim, 1.0_dp, dummyM2, spaceDim, dummyM,&
        & spaceDim, 0.0_dp, matInvOut, memDim)

  end subroutine calcMatrixSqrt


  !> Perform modified Gram-Schmidt orthonormalization of vectors in columns of vec(1:end). Assume
  !> vectors 1:(start-1) are already orthonormal
  pure subroutine orthonormalizeVectors(start, end, vec)

    !> Starting place in vectors to work from
    integer, intent(in) :: start

    !> Ending place in vectors
    integer, intent(in) :: end

    !> Vectors to be orthogonalized against 1:end vectors
    real(dp), intent(inout) :: vec(:,:)

    integer :: ii, jj

    do ii = start, end
      do jj = 1, ii - 1
        vec(:,ii) = vec(:,ii) - dot_product(vec(:,ii), vec(:,jj)) * vec(:,jj)
      end do
      vec(:,ii) = vec(:,ii) / sqrt(dot_product(vec(:,ii), vec(:,ii)))
    end do

  end subroutine orthonormalizeVectors


  !> Encapsulate memory expansion for Stratmann solver
  subroutine incMemStratmann(memDim, workDim, vecB, vP, vM, mP, mM, mH, mMsqrt, mMsqrtInv, &
       &  dummyM, evalInt, workArray, evecL, evecR, vecNorm)

    !> size of subspace
    integer, intent(inout) :: memDim

    !> Work-space large enough for the subspace
    integer, intent(inout) :: workDim

    !> basis of subspace
    real(dp), allocatable, intent(inout) :: vecB(:,:)

    !> vec. for (A+B)b_i
    real(dp), allocatable, intent(inout) :: vP(:,:)

    !> vec. for (A-B)b_i
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

    !> Workspace array
    real(dp), allocatable, intent(inout) :: workArray(:)

    !> Left eigenvectors
    real(dp), allocatable, intent(inout) :: evecL(:,:)

    !> Right eigenvectors
    real(dp), allocatable, intent(inout) :: evecR(:,:)

    !> Norm of eigen vectors
    real(dp), allocatable, intent(inout) :: vecNorm(:)

    call incSizeMatDimTwo(memDim, 3, vecB)
    call incSizeMatDimTwo(memDim, 3, vP)
    call incSizeMatDimTwo(memDim, 3, vM)
    call incSizeMatBothDim(memDim, 3, 2, mP)
    call incSizeMatBothDim(memDim, 3, 2, mM)
    call incSizeMatBothDim(memDim, 3, 2, mH)
    call incSizeMatBothDim(memDim, 3, 2, mMsqrt)
    call incSizeMatBothDim(memDim, 3, 2, mMsqrtInv)
    call incSizeMatBothDim(memDim, 3, 2, dummyM)
    call incSizeVec(memDim, 3, evalInt)
    call incSizeVec(workDim, 3, workArray)
    call incSizeMatDimOne(memDim, 3, evecL)
    call incSizeMatDimOne(memDim, 3, evecR)
    call incSizeVec(2 * memDim, 3, vecNorm)
    memDim = 3 * memDim
    workDim = 3 * workDim

  end subroutine incMemStratmann

end module dftbp_timedep_linrespcommon
