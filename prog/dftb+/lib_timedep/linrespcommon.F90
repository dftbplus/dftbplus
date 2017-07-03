!--------------------------------------------------------------------------------------------------!
!  DFTB+: general package for performing fast atomistic simulations                                !
!  Copyright (C) 2017  DFTB+ developers group                                                      !
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

  real(dp), parameter :: twothird = 2.0_dp / 3.0_dp ! prefactor of 2/3.

contains

  !> Sort arrays in reverse order for values in sposz array
  !! \param wij Energies of transitions
  !! \param sposz array to be sorted containing transition strengths
  !! \param win index array for transitions to single particle transitions
  !! \param transd dipole moments (to be ordered on first index
  !! according to sposz sorting)
  subroutine dipsort(wij, sposz, win, transd)
    real(dp), intent(inout) :: wij(:)
    real(dp), intent(inout) :: sposz(:)
    integer, intent(inout)  :: win(:)
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
  !> transition strength than a tollerance, count them and place at
  !> the start of the appropriate arrays
  subroutine dipSelect(wij,sposz,win,transd,nxov_r,threshold,grndEigVal, getij)
    real(dp), intent(inout) :: wij(:)
    real(dp), intent(inout) :: sposz(:)
    integer, intent(inout)  :: win(:)
    real(dp), intent(inout) :: transd(:,:)
    integer, intent(out)    :: nxov_r
    real(dp), intent(in)    :: threshold
    real(dp), intent(in)    :: grndEigVal(:,:)
    integer, intent(in)     :: getij(:,:)

    integer  :: nxov, ii, jj, iOcc, iVrt, iStart
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
    do ii = 2, nxov
      call indxov(win, ii, getij, iOcc, iVrt)
      if ( abs(grndEigVal(iOcc,1)-eOcc) > epsilon(0.0) .or. &
          & abs(wij(ii)-eExc) > epsilon(0.0) ) then
        eOcc = grndEigVal(iOcc,1)
        eExc = wij(ii)
        mu = 0.0_dp
        do jj = iStart, ii-1
          call indxov(win, jj, getij, iOcc, iVrt)
          mu = mu + sposz(jj)
        end do
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

  !> counts number of transition involving either the HOMO or LUMO levels.
  !! Used to find number of states in the occupied and empty spaces.
  !! \param nxov number of excitations
  !! \param win index array after sorting of eigenvalues
  !! \param getij Index array of transitions
  !! \param homo index of highest occupied level
  !! \param no number of occupied states with transitions into LUMO
  !! \param nv number of virtual states involved in transitions from HOMO
  subroutine getNorb_r(nxov, win, getij, homo, no, nv)
    integer, intent(in) :: nxov, win(:), getij(:,:), homo
    integer, intent(out) :: no, nv

    integer :: ia, ii, jj

    no = 0
    nv = 0

    do ia = 1, nxov
      call indxov(win, ia, getij, ii, jj)
      if (ii == homo) nv = nv +1
      if (jj == homo +1) no = no +1
    end do

  end subroutine getNorb_r

  !> Computes individual indices from the compound occ-virt excitation index.
  !! \param win index array after sorting of eigenvalues.
  !! \param indx  Compound excitation index.
  !! \param getij Index array of transitions after sorting of eigenvalues
  !! \param ii  Initial (filled) state.
  !! \param jj  Final (empty) state.
  subroutine indxov(win, indx, getij, ii, jj)
    integer, intent(in) :: win(:), indx
    integer, intent(in) :: getij(:,:)
    integer, intent(out) :: ii, jj

    integer :: indo

    indo = win(indx)
    ii = getij(indo,1)
    jj = getij(indo,2)

  end subroutine indxov

  !> Computes individual indices from the compound occ-occ excitation index.
  !! \param nocc number of occupied states
  !! \param nocc_r number of required occupied-occupied transitions states
  !! \param indx Compound excitation index.
  !! \param ii Initial state.
  !! \param jj Final state.
  subroutine indxoo(nocc, nocc_r, indx, ii, jj)
    integer, intent(in) ::  nocc, nocc_r, indx
    integer, intent(out) :: ii, jj

    call indxvv(nocc - nocc_r, indx, ii, jj)

  end subroutine indxoo

  !> Computes individual indices from the compound virtual-virtual
  !! excitation index.
  !! \param indx  Compund excitation index.
  !! \param nocc  Number of occupied states.
  !! \param ii  Initial state.
  !! \param jj  Final state.
  subroutine indxvv(nocc, indx, ii, jj)
    integer, intent(in) ::  nocc, indx
    integer, intent(out) :: ii, jj

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
  !! \param nocc  Number of occupied states.
  !! \param ii  First virtual state.
  !! \param jj  second virtual state.
  !! \param indx  Compound indx.
  subroutine rindxvv(nocc, ii, jj, indx)
    integer, intent(in) :: nocc, ii, jj
    integer, intent(out) :: indx

    indx = (jj - nocc) + (((ii - nocc) - 1) * (ii - nocc)) / 2

  end subroutine rindxvv

  !> Builds array to convert from orbital pairs to compound index
  !! \param win array for indexing excitations
  !! \param nocc number of occupied states
  !! \param nxov Number of transitions from occupied to virtual
  !! \param getij array of the occupied-virtual pairs
  !! \param iatrans resulting index array
  subroutine rindxov_array(win, nocc, nxov, getij, iatrans)
    integer, intent(in)  :: win(:)
    integer, intent(in)  :: nocc
    integer, intent(in)  :: nxov
    integer, intent(in)  :: getij(:,:)
    integer, intent(out) :: iatrans(:,nocc+1:)

    integer :: ia

    @:ASSERT(size(getij,dim=2) == 2)

    ! Store reverse indices.
    ! BA: If wij would not be sorted, it would be a trivial transformation.
    do ia = 1, nxov
      iatrans( getij(win(ia),1), getij(win(ia),2) ) = ia
    end do

  end subroutine rindxov_array

  !> Calculates atomic transition charges for a certain excitation.
  !! Calculates qij = 0.5 * (c_i S c_j + c_j S c_i) where c_i and c_j are
  !! selected eigenvectors, and S the overlap matrix. Since qij is atomic
  !! quantity (so far) the corresponding values are summed up.
  !! \param ii  Index of inital state.
  !! \param jj  Index of final state.
  !! \param iAtomStart Starting position of each atom in the list of orbitals.
  !! \param stimc  Overlap times eigenvector: sum_m Smn cmi (nOrb, nOrb).
  !! \param grndEigVecs  Eigenvectors (nOrb, nOrb)
  !! \param qij  Transition charge on exit. (nAtom)
  !! \note ADG: the parameters 'updwn' were added for spin alpha and beta
  !! channels.
  subroutine transq(ii, jj, iAtomStart, updwn, stimc, grndEigVecs, qij)
    integer, intent(in) :: ii, jj, iAtomStart(:)
    logical, intent(in) :: updwn
    real(dp), intent(in) :: stimc(:,:,:), grndEigVecs(:,:,:)
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

  !> Returns the (spatial) MO overlap between orbitals in different
  !> spin channels
  !! \param pp orbital in alpha channel
  !! \param qq orbital in  beta channel
  !! \param stimc S * ci
  !! \param grndEigVecs ci
  function MOoverlap(pp, qq, stimc, grndEigVecs) result(S_pq)
    integer, intent(in) :: pp, qq
    real(dp), intent(in) :: stimc(:,:,:), grndEigVecs(:,:,:)
    real(dp):: S_pq

    @:ASSERT(all(shape(grndEigVecs) == shape(stimc)))
    @:ASSERT(size(grndEigVecs, dim=3) == 2)

    S_pq = sum(grndEigVecs(:,pp,1)*stimc(:,qq,2))

  end function MOoverlap

  !> scale excitation energies by differences in occupations between
  !> filled and empty states
  !! \param wij K-S energies
  !! \param occNr occupation of states
  !! \param win index array
  !! \param nmatup number of up-up excitations
  !! \param nmat number of transitions to scale
  !! \param getij index array
  !! \param wn_ij resulting scaled matrix
  subroutine wtdn(wij, occNr, win, nmatup, nmat, getij, wn_ij)
    real(dp), intent(in) :: wij(:)
    real(dp), intent(in) :: occNr(:,:)
    integer, intent(in) :: win(:), nmatup, nmat, getij(:,:)
    real(dp), intent(out) :: wn_ij(:)

    integer :: ia, ii, jj
    logical :: updwn
    real(dp) :: docc_ij

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

  end subroutine wtdn

  !> Multiplies the excitation supermatrix with a supervector.
  !! \note In order not to store the entire supermatrix (nmat, nmat), the
  !! various pieces are assembled individually and multiplied directly
  !! with the corresponding part of the supervector.
  !! \param spin logical spin polarization
  !! \param vin  Vector to multiply with. (nmat)
  !! \param vout  Vector containing the result on exit. (nmat,)
  !! \param wij  Excitation energies (wij = epsion_j - epsilon_i)
  !! \param sym  Symmetry flag (singlet or triplet)
  !! \param win  Sorting index of the excitation energies.
  !! \param nmatup number of same-spin transitions
  !! \param iAtomStart Starting position of each atom in the list of orbitals.
  !! \param stimc overlap times eigenvector. (nOrb, nOrb)
  !! \param grndEigVecs  Eigenvectors (nOrb, nOrb)
  !! \param occNr Occupation numbers
  !! \param getij index array for excitations
  !! \param gamma DFTB gamma matrix (nAtm, nAtom)
  !! \param species0 Chemical species of the atoms
  !! \param spinW ground state spin derivatives for each species
  subroutine omegatvec(spin, vin, vout, wij, sym, win, nmatup, iAtomStart,&
      & stimc, grndEigVecs, occNr, getij, gamma, species0, spinW )
    logical, intent(in)   :: spin
    real(dp), intent(in)  :: vin(:)
    real(dp), intent(out) :: vout(:)
    real(dp), intent(in)  :: wij(:)
    character, intent(in) :: sym
    integer, intent(in)   :: win(:), nmatup, iAtomStart(:)
    real(dp), intent(in)  :: stimc(:,:,:), grndEigVecs(:,:,:)
    real(dp), intent(in)  :: occNr(:,:)
    real(dp), intent(in)  :: gamma(:,:)
    integer,  intent(in)  :: getij(:,:)
    integer, intent(in)   :: species0(:)
    real(dp), intent(in)  :: spinW(:)

    integer :: nmat, natom
    integer :: ia, ii, jj
    real(dp) :: fact
    ! somewhat ugly, but fast small arrays on stack:
    real(dp) :: otmp(size(gamma, dim=1)), gtmp(size(gamma, dim=1))
    real(dp) :: qij(size(gamma, dim=1)) ! qij Working array (used for
                                        ! excitation charges). (nAtom)
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
  !! \param rkm1  Resulting vector on return.
  !! \param rhs2  Vector, which (A+B) should act on.
  !! \param wij Excitation energies
  !! \param nmat matrix dimension
  !! \param nAtom number of atoms
  !! \param win index array
  !! \param nmatup number of same spin transitions
  !! \param getij index array
  !! \param iAtomStart indexing array for atomic basis functions
  !! \param stimc overlap times ground state wavefunction
  !! \param grndEigVecs ground state wave function coefficients
  !! \param qij transition charges
  !! \param gamma gamma matrix
  subroutine apbw(rkm1, rhs2, wij, nmat, natom, win, nmatup, getij, iAtomStart,&
      & stimc, grndEigVecs, gamma)
    real(dp), intent(out) :: rkm1(:)
    real(dp), intent(in) :: rhs2(:), wij(:)
    integer, intent(in) :: nmat, natom, win(:)
    integer, intent(in) :: nmatup, getij(:,:), iAtomStart(:)
    real(dp), intent(in) :: stimc(:,:,:), grndEigVecs(:,:,:)
    real(dp), intent(in) :: gamma(:,:)

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
  !! \param eigval ground state eigenvalues
  !! \param filling occupation numbers for the ground state
  !! \param wij Kohn-Sham energy differences between empty and filled states
  !! \param getij index of pairs of KS states for the transitions in wij
  !! \note ADG: the subroutine is generalized to account for spin and
  !! partial occupancy
  subroutine getSPExcitations(grndEigVal, filling, wij, getij)
    real(dp), intent(in) :: grndEigVal(:,:)
    real(dp), intent(in) :: filling(:,:)
    real(dp), intent(out) :: wij(:)
    integer,  intent(out) :: getij(:,:)

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
  !! \param transd transition dipole matrix elements for each atom
  !! \param wij energy differences
  !! \parm evec eigenvectors
  function oscillatorStrength(transd, wij, evec) result(osz)
    real(dp), intent(in) :: transd(:,:), wij(:), evec(:)
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

  !> Transition dipoles between states
  !! \param transd transition dipole for transitions between KS states
  !! \param eigenvalues
  !! \param eigenvectors
  !! \param tdip resulting transition dipoles
  subroutine transitionDipole(transd, wij, eval, evec, tdip)
    real(dp), intent(in) :: transd(:,:), wij(:), eval(:), evec(:,:)
    real(dp), intent(out) :: tdip(:,:)

    real(dp) :: rtmp
    integer :: ii, ll

    tdip(:,:) = 0.0_dp
    do ii = 1, size(evec, dim=2)
      rtmp = sqrt(sqrt(eval(ii)))
      do ll = 1, 3
        tdip(ii,ll) = sum(transd(:,ll) * sqrt(wij) * evec(:,ii)) / rtmp
      end do
    end do

  end subroutine transitionDipole

  !> count initial number of transitions from occupied to empty states
  !! \param nOrb Number of ground state orbitals
  !! \param nSpin Number of spin channels
  !! \param filling occupations of ground state SP levels
  function countSPexcitations(nOrb,nSpin,filling) result(nxov)
    integer, intent(in)  :: nOrb, nSpin
    real(dp), intent(in) :: filling(:,:)
    real(dp) :: nxov

    real(dp) :: nxov_ud(nSpin)
    integer :: ii, jj, iSpin

    @:ASSERT(size(filling,dim=1)==nOrb)
    @:ASSERT(size(filling,dim=2)==nSpin)

    nxov_ud = 0
    do iSpin = 1, nSpin
      do ii = 1, norb - 1
        do jj = ii, norb
          if (filling(ii,iSpin) > filling(jj,iSpin) + elecTolMax) then
            nxov_ud(iSpin) = nxov_ud(iSpin) + 1
          end if
        end do
      end do
    end do
    nxov = sum(nxov_ud)

  end function countSPexcitations

end module linrespcommon
