!--------------------------------------------------------------------------------------------------!
!  DFTB+: general package for performing fast atomistic simulations                                !
!  Copyright (C) 2018  DFTB+ developers group                                                      !
!                                                                                                  !
!  See the LICENSE file for terms of usage and distribution.                                       !
!--------------------------------------------------------------------------------------------------!

#:include 'common.fypp'

!> Helper routines for the linear response modules.
module linrespcommon
  use assert
  use accuracy
  use blasroutines
  use sorting
  use message
  use commontypes
  implicit none
  public


  !> prefactor of 2/3.
  real(dp), parameter :: twothird = 2.0_dp / 3.0_dp

contains


  !> Sort arrays in reverse order on basis of values in the sposz array
  subroutine dipsort(wij, sposz, win, transd)

    !> Energies of transitions
    real(dp), intent(inout) :: wij(:)

    !> array to be sorted containing single particle transition strengths
    real(dp), intent(inout) :: sposz(:)

    !> index array for transitions to single particle transitions
    integer, intent(inout) :: win(:)

    !> dipole moments (to be ordered on first index according to sposz sorting)
    real(dp), intent(inout) :: transd(:,:)

    integer, allocatable :: tmpIndx(:)

    allocate(tmpIndx(size(win)))

    @:ASSERT(size(wij) == size(sposz))
    @:ASSERT(size(wij) == size(win))
    @:ASSERT(size(wij) == size(transd,dim=1))

    call index_heap_sort(tmpIndx, sposz)
    tmpIndx = tmpIndx(size(win):1:-1)
    win = win(tmpIndx)
    wij = wij(tmpIndx)
    sposz = sposz(tmpIndx)
    transd(:,:) = transd(tmpIndx,:)

  end subroutine dipsort


  !> find (possibly degenerate) transitions with stronger dipole
  !> transition strengths than a tolerance, count them and place at
  !> the start of the appropriate arrays
  subroutine dipSelect(wij,sposz,win,transd,nxov_r,threshold,grndEigVal, getij)

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
    integer, intent(in) :: getij(:,:)

    integer :: nxov, ii, jj, iOcc, iVrt, iStart
    real(dp) :: eOcc, eExc, mu

    nxov = size(wij)
    @:ASSERT(nxov == size(sposz))
    @:ASSERT(nxov == size(win))
    @:ASSERT(all(shape(transd) == [nxov,3]))
    @:ASSERT(size(getij,dim=1) >= nxov)
    @:ASSERT(size(getij,dim=2) == 2)
    @:ASSERT(threshold >= 0.0_dp)

    call indxov(win, 1, getij, iOcc, iVrt)
    eOcc = grndEigVal(iOcc,1)
    eExc = wij(1)
    iStart = 1
    nxov_r = 0
    ! Check, to single precision tolerance, for degeneracies when selecting bright transitions
    do ii = 2, nxov
      call indxov(win, ii, getij, iOcc, iVrt)
      ! check if this is a still within a degenerate group, otherwise process the group
      if ( abs(grndEigVal(iOcc,1)-eOcc) > epsilon(0.0) .or. &
          & abs(wij(ii)-eExc) > epsilon(0.0) ) then
        eOcc = grndEigVal(iOcc,1)
        eExc = wij(ii)
        mu = 0.0_dp
        ! loop over transitions in the group and check the oscillator strength
        do jj = iStart, ii-1
          call indxov(win, jj, getij, iOcc, iVrt)
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
      call indxov(win, jj, getij, iOcc, iVrt)
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


  !> counts number of transition involving either the HOMO or LUMO levels.  Used to find number of
  !> states in the occupied and empty spaces.
  subroutine getNorb_r(nxov, win, getij, homo, no, nv)

    !> number of excitations
    integer, intent(in) :: nxov

    !> index array after sorting of eigenvalues
    integer, intent(in) :: win(:)

    !> Index array of transitions
    integer, intent(in) :: getij(:,:)

    !> index of highest occupied level
    integer, intent(in) :: homo

    !> number of occupied states with transitions into LUMO
    integer, intent(out) :: no

    !> number of virtual states involved in transitions from HOMO
    integer, intent(out) :: nv

    integer :: ia, ii, jj

    no = 0
    nv = 0

    !$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(ia,ii,jj) SCHEDULE(RUNTIME) REDUCTION(+:nv,no)
    do ia = 1, nxov
      call indxov(win, ia, getij, ii, jj)
      if (ii == homo) nv = nv +1
      if (jj == homo +1) no = no +1
    end do
    !$OMP  END PARALLEL DO

  end subroutine getNorb_r


  !> Computes individual indices from the compound occ-virt excitation index.
  subroutine indxov(win, indx, getij, ii, jj)

    !> index array after sorting of eigenvalues.
    integer, intent(in) :: win(:)

    !> Compound excitation index.
    integer, intent(in) :: indx

    !> Index array of transitions after sorting of eigenvalues
    integer, intent(in) :: getij(:,:)

    !> Initial (filled) state.
    integer, intent(out) :: ii

    !> Final (empty) state.
    integer, intent(out) :: jj

    integer :: indo

    indo = win(indx)
    ii = getij(indo,1)
    jj = getij(indo,2)

  end subroutine indxov


  !> Computes individual indices from the compound occ-occ excitation index.
  subroutine indxoo(nocc, nocc_r, indx, ii, jj)

    !> number of occupied states
    integer, intent(in) :: nocc

    !> number of required occupied-occupied transitions states
    integer, intent(in) :: nocc_r

    !> Compound excitation index.
    integer, intent(in) :: indx

    !> Initial state.
    integer, intent(out) :: ii

    !> Final state.
    integer, intent(out) :: jj

    call indxvv(nocc - nocc_r, indx, ii, jj)

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


  !> Computes the virtual-virtual compound index from two indices.
  subroutine rindxvv(nocc, ii, jj, indx)

    !> Number of occupied states.
    integer, intent(in) :: nocc

    !> First virtual state.
    integer, intent(in) :: ii

    !> second virtual state.
    integer, intent(in) :: jj

    !> Compound indx.
    integer, intent(out) :: indx

    indx = (jj - nocc) + (((ii - nocc) - 1) * (ii - nocc)) / 2

  end subroutine rindxvv


  !> Builds array to convert from orbital pairs to a compound index, reverse of index generated by
  !> getSPExcitations (which were then potentially windowed afterwards)
  subroutine rindxov_array(win, nocc, nxov, getij, iatrans)

    !> array for indexing excitations
    integer, intent(in) :: win(:)

    !> number of occupied states
    integer, intent(in) :: nocc

    !> Number of transitions from occupied to virtual
    integer, intent(in) :: nxov

    !> array of the occupied->virtual pairs
    integer, intent(in) :: getij(:,:)

    !> resulting index array from orbital pairs to compound index
    integer, intent(out) :: iatrans(:,nocc+1:)

    integer :: ia

    @:ASSERT(size(getij,dim=2) == 2)

    ! Store reverse indices

    ! If wij was not sorted, it would be a trivial transformation.
    !$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(ia) SCHEDULE(RUNTIME)
    do ia = 1, nxov
      iatrans( getij(win(ia),1), getij(win(ia),2) ) = ia
    end do
    !$OMP  END PARALLEL DO

  end subroutine rindxov_array


  !> Calculates atomic transition charges for a certain excitation.
  !> Calculates qij = 0.5 * (c_i S c_j + c_j S c_i) where c_i and c_j are selected eigenvectors, and
  !> S the overlap matrix. Since qij is atomic quantity (so far) the corresponding values are summed
  !> up.
  !> Note: the parameters 'updwn' were added for spin alpha and beta channels.
  subroutine transq(ii, jj, iAtomStart, updwn, stimc, grndEigVecs, qij)

    !> Index of inital state.
    integer, intent(in) :: ii

    !> Index of final state.
    integer, intent(in) :: jj

    !> Starting position of each atom in the list of orbitals.
    integer, intent(in) :: iAtomStart(:)

    !> up spin channel (T) or down spin channel (F)
    logical, intent(in) :: updwn

    !> Overlap times eigenvector: sum_m Smn cmi (nOrb, nOrb).
    real(dp), intent(in) :: stimc(:,:,:)

    !> Eigenvectors (nOrb, nOrb)
    real(dp), intent(in) :: grndEigVecs(:,:,:)

    !> Transition charge on exit. (nAtom)
    real(dp), intent(out) :: qij(:)

    integer :: kk, aa, bb, ss
    real(dp) :: qTmp(size(stimc,dim=1))

    @:ASSERT(all(shape(stimc) == shape(grndEigVecs)))

    ss = 1
    if (.not. updwn) ss = 2
    qTmp(:) =  grndEigVecs(:,ii,ss) * stimc(:,jj,ss) &
        & + grndEigVecs(:,jj,ss) * stimc(:,ii,ss)
    do kk = 1, size(qij)
      aa = iAtomStart(kk)
      bb = iAtomStart(kk + 1) -1
      qij(kk) = 0.5_dp * sum(qTmp(aa:bb))
    end do

  end subroutine transq


  !> Returns the (spatial) MO overlap between orbitals in different spin channels
  function MOoverlap(pp, qq, stimc, grndEigVecs) result(S_pq)

    !> orbital in alpha channel
    integer, intent(in) :: pp

    !> orbital in  beta channel
    integer, intent(in) :: qq

    !> overlap times ground single particle state wavefunctiona
    real(dp), intent(in) :: stimc(:,:,:)

    !> ground state single particle wavefunctions
    real(dp), intent(in) :: grndEigVecs(:,:,:)

    !> Resulting overlap between states
    real(dp):: S_pq

    @:ASSERT(all(shape(grndEigVecs) == shape(stimc)))
    @:ASSERT(size(grndEigVecs, dim=3) == 2)

    S_pq = sum(grndEigVecs(:,pp,1)*stimc(:,qq,2))

  end function MOoverlap


  !> Scale excitation energies by differences in occupations between filled and empty states
  subroutine wtdn(wij, occNr, win, nmatup, nmat, getij, wn_ij)

    !> K-S energies
    real(dp), intent(in) :: wij(:)

    !> occupation of states
    real(dp), intent(in) :: occNr(:,:)

    !> index array for transitions to single particle transitions
    integer, intent(in) :: win(:)

    !> number of up-up excitations
    integer, intent(in) :: nmatup

    !> number of transitions to scale
    integer, intent(in) :: nmat

    !> array of the occupied->virtual pairs
    integer, intent(in) :: getij(:,:)

    !> resulting scaled matrix
    real(dp), intent(out) :: wn_ij(:)

    integer :: ia, ii, jj
    logical :: updwn
    real(dp) :: docc_ij

    !$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(ia, ii, jj, updwn, docc_ij) SCHEDULE(RUNTIME)
    do ia = 1, nmat
      call indxov(win, ia, getij, ii, jj)
      updwn = (win(ia) <= nmatup)
      if(updwn) then
        docc_ij = occNr(ii,1) - occNr(jj,1)
      else
        docc_ij = occNr(ii,2) - occNr(jj,2)
      end if
      wn_ij(ia) = wij(ia) * docc_ij
    end do
    !$OMP  END PARALLEL DO

  end subroutine wtdn


  !> Multiplies the excitation supermatrix with a supervector.
  !>
  !> Note: In order not to store the entire supermatrix (nmat, nmat), the various pieces are
  !> assembled individually and multiplied directly with the corresponding part of the supervector.
  subroutine omegatvec(spin, vin, vout, wij, sym, win, nmatup, iAtomStart, stimc, grndEigVecs, &
      & occNr, getij, gamma, species0, spinW )

    !> logical spin polarization
    logical, intent(in) :: spin

    !> Vector to multiply with. (nmat)
    real(dp), intent(in) :: vin(:)

    !> Vector containing the result on exit. (nmat,)
    real(dp), intent(out) :: vout(:)

    !> Excitation energies (wij = epsion_j - epsilon_i)
    real(dp), intent(in) :: wij(:)

    !> Symmetry flag (singlet or triplet)
    character, intent(in) :: sym

    !> Sorting index of the excitation energies.
    integer, intent(in) :: win(:)

    !> number of same-spin transitions
    integer, intent(in) :: nmatup

    !> Starting position of each atom in the list of orbitals.
    integer, intent(in) :: iAtomStart(:)

    !> overlap times eigenvector. (nOrb, nOrb)
    real(dp), intent(in) :: stimc(:,:,:)

    !> Eigenvectors (nOrb, nOrb)
    real(dp), intent(in) :: grndEigVecs(:,:,:)

    !> Occupation numbers
    real(dp), intent(in) :: occNr(:,:)

    !> index array for excitations
    real(dp), intent(in) :: gamma(:,:)

    !> DFTB gamma matrix (nAtm, nAtom)
    integer, intent(in) :: getij(:,:)

    !> Chemical species of the atoms
    integer, intent(in) :: species0(:)

    !> ground state spin constants for each species
    real(dp), intent(in) :: spinW(:)

    integer :: nmat, natom
    integer :: ia, ii, jj
    real(dp) :: fact
    ! somewhat ugly, but fast small arrays on stack:
    real(dp) :: otmp(size(gamma, dim=1)), gtmp(size(gamma, dim=1))
    real(dp) :: qij(size(gamma, dim=1)) ! qij Working array (used for excitation charges). (nAtom)
    real(dp) :: wnij(size(wij))
    logical :: updwn

    @:ASSERT(size(vin) == size(vout))

    nmat = size(vin)
    natom = size(gamma, dim=1)

    vout = 0.0_dp

    call wtdn(wij, occNr, win, nmatup, nmat, getij, wnij)
    wnij = sqrt(wnij) ! always used as root(wnij) after this

    otmp = 0.0_dp
    !$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(ia,ii,jj,updwn,qij) &
    !$OMP& SCHEDULE(RUNTIME) REDUCTION(+:otmp)
    do ia = 1, nmat
      call indxov(win, ia, getij, ii, jj)
      updwn = (win(ia) <= nmatup)
      call transq(ii, jj, iAtomStart, updwn, stimc, grndEigVecs, qij)
      otmp = otmp + vin(ia) * wnij(ia) * qij
    end do
    !$OMP  END PARALLEL DO

    if (.not.spin) then !-----------spin-unpolarized systems--------------

      if (sym == 'S') then

        call hemv(gtmp, gamma, otmp)

        !$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(ia,ii,jj,updwn,qij) &
        !$OMP& SCHEDULE(RUNTIME)
        do ia = 1, nmat
          call indxov(win, ia, getij, ii, jj)
          updwn = (win(ia) <= nmatup)
          call transq(ii, jj, iAtomStart, updwn, stimc, grndEigVecs, qij)
          vout(ia) = 2.0_dp * wnij(ia) * dot_product(qij, gtmp)
        end do
        !$OMP  END PARALLEL DO
      else

        otmp = 2.0_dp * otmp * spinW(species0)

        !$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(ia,ii,jj,updwn,qij) &
        !$OMP& SCHEDULE(RUNTIME)
        do ia = 1, nmat
          call indxov(win, ia, getij, ii, jj)
          updwn = (win(ia) <= nmatup)
          call transq(ii, jj, iAtomStart, updwn, stimc, grndEigVecs, qij)
          ! Note: 2 times atomic magnetization m_A
          ! vout = sum_A q_A^ia m_A * otmp(A)
          vout(ia) = vout(ia) + wnij(ia) * dot_product(qij, otmp)
        end do
        !$OMP  END PARALLEL DO

      end if

    else !--------------spin-polarized systems--------

      call hemv(gtmp, gamma, otmp)

      otmp(:) = 0.0_dp

      !$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(ia,ii,jj,updwn,qij,fact) &
      !$OMP& SCHEDULE(RUNTIME) REDUCTION(+:otmp)
      do ia = 1,nmat

        call indxov(win, ia, getij, ii, jj)
        updwn = (win(ia) <= nmatup)

        call transq(ii, jj, iAtomStart, updwn, stimc, grndEigVecs, qij)

        !singlet gamma part (S)
        vout(ia) = 2.0_dp * wnij(ia) * dot_product(qij, gtmp)

        !magnetization part (T1)
        if (updwn) then
          fact = 1.0_dp
        else
          fact =-1.0_dp
        end if
        otmp(:) = otmp(:) + fact * vin(ia) * wnij(ia) * qij(:)

      end do
      !$OMP  END PARALLEL DO

      otmp = 2.0_dp * otmp * spinW(species0)

      !$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(ia,ii,jj,updwn,qij,fact) &
      !$OMP& SCHEDULE(RUNTIME)
      do ia = 1,nmat

        call indxov(win, ia, getij, ii, jj)
        updwn = (win(ia) <= nmatup)
        call transq(ii, jj, iAtomStart, updwn, stimc, grndEigVecs, qij)
        if (updwn) then
          fact = 1.0_dp
        else
          fact =-1.0_dp
        end if
        vout(ia) = vout(ia) + fact * wnij(ia) * dot_product(qij, otmp)

      end do
      !$OMP  END PARALLEL DO

    end if

    vout = vout + ( wij**2 ) * vin

  end subroutine omegatvec


  !> Multiplies the supermatrix (A+B) with a given vector.
  subroutine apbw(rkm1, rhs2, wij, nmat, natom, win, nmatup, getij, iAtomStart, stimc, grndEigVecs,&
      & gamma)

    !> Resulting vector on return.
    real(dp), intent(out) :: rkm1(:)

    !> Vector, which (A+B) should act on.
    real(dp), intent(in) :: rhs2(:)

    !> Excitation energies
    real(dp), intent(in) :: wij(:)

    !> matrix dimension
    integer, intent(in) :: nmat

    !> number of atoms
    integer, intent(in) :: natom

    !> index array
    integer, intent(in) :: win(:)

    !> number of same spin transitions
    integer, intent(in) :: nmatup

    !> index array
    integer, intent(in) :: getij(:,:)

    !> indexing array for atomic basis functions
    integer, intent(in) :: iAtomStart(:)

    !> overlap times ground state wavefunction
    real(dp), intent(in) :: stimc(:,:,:)

    !> ground state wave function coefficients
    real(dp), intent(in) :: grndEigVecs(:,:,:)

    !> transition charges
    real(dp), intent(in) :: gamma(:,:)

    !> gamma matrix
    integer :: ia, ii, jj
    real(dp) :: tmp(natom), gtmp(natom), qij(natom)
    logical :: updwn

    @:ASSERT(size(rkm1) == nmat)

    tmp(:) = 0.0_dp
    !$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(ia,ii,jj,updwn,qij)&
    !$OMP& SCHEDULE(RUNTIME) REDUCTION(+:tmp)
    do ia = 1, nmat
      call indxov(win, ia, getij, ii, jj)
      updwn = (win(ia) <= nmatup)
      call transq(ii, jj, iAtomStart, updwn, stimc, grndEigVecs, qij)
      tmp(:) = tmp + rhs2(ia) * qij
    end do
    !$OMP  END PARALLEL DO

    call hemv(gtmp, gamma, tmp)

    rkm1 = 0.0_dp
    !$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(ia,ii,jj,updwn,qij)&
    !$OMP& SCHEDULE(RUNTIME)
    do ia = 1, nmat
      call indxov(win, ia, getij, ii, jj)
      updwn = (win(ia) <= nmatup)
      call transq(ii, jj, iAtomStart, updwn, stimc, grndEigVecs, qij)
      rkm1(ia) = 4.0_dp * dot_product(gtmp, qij)
    end do
    !$OMP  END PARALLEL DO

    rkm1 = rkm1 + wij * rhs2

  end subroutine apbw


  !> calculating spin polarized excitations
  !> Note: the subroutine is generalized to account for spin and partial occupancy
  subroutine getSPExcitations(grndEigVal, filling, wij, getij)

    !> ground state eigenvalues
    real(dp), intent(in) :: grndEigVal(:,:)

    !> occupation numbers for the ground state
    real(dp), intent(in) :: filling(:,:)

    !> Kohn-Sham energy differences between empty and filled states
    real(dp), intent(out) :: wij(:)

    !> index of pairs of KS states for the transitions in wij
    integer, intent(out) :: getij(:,:)

    integer :: ind, ii, jj
    integer :: norb, iSpin, nSpin

    @:ASSERT(all(shape(grndEigVal)==shape(filling)))

    norb = size(grndEigVal, dim=1)
    nSpin = size(grndEigVal, dim=2)

    ind = 0

    do iSpin = 1, nSpin
      do ii = 1, norb - 1
        do jj = ii, norb
          if (filling(ii,iSpin) > filling(jj,iSpin) + elecTolMax) then
            ind = ind + 1
            wij(ind) = grndEigVal(jj,iSpin) - grndEigVal(ii,iSpin)
            getij(ind,:) = [ii,jj]
          end if
        end do
      end do
    end do

  end subroutine getSPExcitations


  !> Calculate dipole moment for transitions
  function oscillatorStrength(transd, wij, evec) result(osz)

    !> transition dipole matrix elements for each atom
    real(dp), intent(in) :: transd(:,:)

    !> energy differences between single particle states
    real(dp), intent(in) :: wij(:)

    !> excitation eigenvector
    real(dp), intent(in) :: evec(:)

    !> oscillator strength
    real(dp) :: osz

    real(dp) :: rtmp
    integer :: ii

    osz = 0.0_dp
    do ii = 1, 3
      rtmp = sum(transd(:,ii) * sqrt(wij) * evec)
      osz = osz + rtmp * rtmp
    end do
    osz = twothird * osz

  end function oscillatorStrength


  !> Transition dipoles between single particle states
  subroutine transitionDipole(transd, wij, eval, evec, tdip)

    !> transition dipole for transitions between KS states
    real(dp), intent(in) :: transd(:,:)

    !> energy differences between single particle states
    real(dp), intent(in) :: wij(:)

    !> single particle excitation eigenvalues
    real(dp), intent(in) :: eval(:)

    !> single particle eigenvectors
    real(dp), intent(in) :: evec(:,:)

    !> resulting transition dipoles
    real(dp), intent(out) :: tdip(:,:)

    real(dp) :: rtmp
    integer :: ii, ll

    tdip(:,:) = 0.0_dp
    !$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(ii,rtmp) SCHEDULE(RUNTIME)
    do ii = 1, size(evec, dim=2)
      rtmp = eval(ii)**(-4) ! 1._dp / sqrt(sqrt(eval(ii)))
      do ll = 1, 3
        tdip(ii,ll) = sum(transd(:,ll) * sqrt(wij) * evec(:,ii)) * rtmp
      end do
    end do
    !$OMP  END PARALLEL DO

  end subroutine transitionDipole

end module linrespcommon
