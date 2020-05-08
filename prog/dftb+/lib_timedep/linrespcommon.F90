!--------------------------------------------------------------------------------------------------!
!  DFTB+: general package for performing fast atomistic simulations                                !
!  Copyright (C) 2006 - 2020  DFTB+ developers group                                               !
!                                                                                                  !
!  See the LICENSE file for terms of usage and distribution.                                       !
!--------------------------------------------------------------------------------------------------!

#:include 'common.fypp'

!> Helper routines for the linear response modules.
module dftbp_linrespcommon
  use dftbp_assert
  use dftbp_accuracy
  use dftbp_blasroutines
  use dftbp_sorting
  use dftbp_message
  use dftbp_commontypes
  use dftbp_transcharges
  use dftbp_onsitecorrection, only : getOnsME
  implicit none
  public


  !> prefactor of 2/3.
  real(dp), parameter :: twothird = 2.0_dp / 3.0_dp

contains


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


  !> Computes individual indices from the compound occ-virt excitation index.
  pure subroutine indxov(win, indx, getij, ii, jj)

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


  !> calculate the transition block at a specific atom
  subroutine transDens(ii, jj, iAt, iAtomStart, nOrb, updwn, stimc, grndEigVecs, qq_ij)

    !> Index of inital state.
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
    real(dp), intent(in) :: stimc(:,:,:)

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
    !call ger(qq_ij(:nOrb,:nOrb), 0.5_dp, grndEigVecs(mu:mu+nOrb-1,ii,ss),stimc(mu:mu+nOrb-1,jj,ss))
    !call ger(qq_ij(:nOrb,:nOrb), 0.5_dp, grndEigVecs(mu:mu+nOrb-1,jj,ss),stimc(mu:mu+nOrb-1,ii,ss))
    !qq_ij(:nOrb,:nOrb) = 0.5_dp * (qq_ij(:nOrb,:nOrb) + transpose(qq_ij(:nOrb,:nOrb)))

    do iOrb1 = 1, nOrb
      do iOrb2 = iOrb1, nOrb
        mu = iAtomStart(iAt) + iOrb1 - 1
        nu = iAtomStart(iAt) + iOrb2 - 1
        qq_ij(iOrb1,iOrb2) = 0.25_dp*( grndEigVecs(mu,ii,ss)*stimc(nu,jj,ss) &
             &                       + grndEigVecs(mu,jj,ss)*stimc(nu,ii,ss) &
             &                       + grndEigVecs(nu,ii,ss)*stimc(mu,jj,ss) &
             &                       + grndEigVecs(nu,jj,ss)*stimc(mu,ii,ss) )
        if (iOrb1 /= iOrb2) then
          qq_ij(iOrb2,iOrb1) = qq_ij(iOrb1,iOrb2)
        end if
      end do
    end do

  end subroutine transDens


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
      & occNr, getij, gamma, species0, spinW, onsMEs, orb, transChrg)

    !> logical spin polarization
    logical, intent(in) :: spin

    !> Vector to multiply with size(nmat)
    real(dp), intent(in) :: vin(:)

    !> Vector containing the result on exit size(nmat)
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

    !> onsite matrix elements for shells (elements between s orbitals on the same shell are ignored)
    real(dp), intent(in), allocatable :: onsMEs(:,:,:,:)

    !> data type for atomic orbital information
    type(TOrbitals), intent(in) :: orb

    !> machinery for transition charges between single particle levels
    type(TTransCharges), intent(in) :: transChrg

    integer :: nmat, natom
    integer :: ia, ii, jj
    real(dp) :: fact
    ! somewhat ugly, but fast small arrays on stack:
    real(dp) :: otmp(size(gamma, dim=1)), gtmp(size(gamma, dim=1))
    real(dp) :: qij(size(gamma, dim=1)) ! qij Working array (used for excitation charges). (nAtom)
    real(dp) :: wnij(size(wij))
    logical :: updwn
    real(dp) :: vout_ons(size(vin))

    @:ASSERT(size(vin) == size(vout))

    nmat = size(vin)
    natom = size(gamma, dim=1)

    vout = 0.0_dp

    call wtdn(wij, occNr, win, nmatup, nmat, getij, wnij)
    wnij = sqrt(wnij) ! always used as root(wnij) after this

    ! product charges with the v*wn product, i.e. Q * v*wn
    oTmp(:) = 0.0_dp
    call transChrg%qMatVec(iAtomStart, stimc, grndEigVecs, getij, win, vin(:nmat) * wnij(:nmat),&
        & oTmp)

    if (.not.spin) then !-----------spin-unpolarized systems--------------

      if (sym == 'S') then

        call hemv(gtmp, gamma, otmp)

        ! 2 * wn * (g * Q)
        vOut(:) = 0.0_dp
        call transChrg%qVecMat(iAtomStart, stimc, grndEigVecs, getij, win, gTmp, vOut)
        vOut(:) = 2.0_dp * wnij(:) * vOut(:)

      else

        otmp = 2.0_dp * otmp * spinW(species0)

        ! wn * (o * Q)
        vOut = 0.0_dp
        call transChrg%qVecMat(iAtomStart, stimc, grndEigVecs, getij, win, oTmp, vOut)
        vOut(:) = wnij(:) * vOut(:)


      end if

    else !--------------spin-polarized systems--------

      call hemv(gtmp, gamma, otmp)

      otmp(:) = 0.0_dp

      !$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(ia,ii,jj,updwn,qij,fact) &
      !$OMP& SCHEDULE(RUNTIME) REDUCTION(+:otmp)
      do ia = 1,nmat

        updwn = (win(ia) <= nmatup)

        qij(:) = transChrg%qTransIJ(ia, iAtomStart, stimc, grndEigVecs, getij, win)

        ! singlet gamma part (S)
        vout(ia) = 2.0_dp * wnij(ia) * dot_product(qij, gtmp)

        ! magnetization part (T1)
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

        qij(:) = transChrg%qTransIJ(ia, iAtomStart, stimc, grndEigVecs, getij, win)

        updwn = (win(ia) <= nmatup)
        if (updwn) then
          fact = 1.0_dp
        else
          fact =-1.0_dp
        end if
        vout(ia) = vout(ia) + fact * wnij(ia) * dot_product(qij, otmp)

      end do
      !$OMP  END PARALLEL DO

    end if

    if (allocated(onsMEs)) then
      call onsiteEner(spin, sym, wnij, win, nmatup, iAtomStart, getij, species0, stimc,&
          & grndEigVecs, onsMEs, orb, vin, vout_ons)
      vout(:) = vout + vout_ons
    end if

    vout(:) = vout(:) + ( wij**2 ) * vin

  end subroutine omegatvec


  subroutine onsiteEner(spin, sym, wnij, win, nmatup, iAtomStart, getij, species0, stimc,&
      & grndEigVecs, ons_en, orb, vin, vout)

    logical, intent(in) :: spin
    character, intent(in) :: sym
    real(dp), intent(in) :: wnij(:)
    integer, intent(in) :: win(:)
    integer, intent(in) :: nmatup
    integer, intent(in) :: iAtomStart(:)
    integer, intent(in) :: getij(:,:)
    integer, intent(in) :: species0(:)
    real(dp), intent(in) :: stimc(:,:,:)
    real(dp), intent(in) :: grndEigVecs(:,:,:)
    real(dp), intent(in) :: ons_en(:,:,:,:)
    type(TOrbitals), intent(in) :: orb
    real(dp), intent(in) :: vin(:)
    real(dp), intent(out) :: vout(:)

    real(dp) :: otmp(orb%mOrb,orb%mOrb,size(species0),2)
    real(dp) :: fact
    real(dp) :: onsite(orb%mOrb,orb%mOrb,2)
    real(dp) :: qq_ij(orb%mOrb,orb%mOrb)
    logical :: updwn
    integer :: nmat, nAtom, nOrb
    integer :: ia, iAt, iSp, iSh, iOrb, iOrb1, iOrb2, ii, jj, mu, nu, ss, sindx(2), iSpin
    real(dp) :: degeneracy, partTrace

    nmat = size(vin)
    nAtom = size(species0)

    vout(:) = 0.0_dp
    otmp(:,:,:,:)  = 0.0_dp

    fact = 1.0_dp
    if (.not.Spin) then
      if (sym == 'T') then
        fact = -1.0_dp
      end if
    end if

    if (spin) then
      ss = 2
    else
      ss = 1
    end if

    do iAt = 1, nAtom
      iSp = species0(iAt)
      nOrb = orb%nOrbAtom(iAt)
      call getOnsME(orb, iSp, ons_en, nOrb, onsite)
      do ia = 1, nmat
        call indxov(win, ia, getij, ii, jj)
        updwn = (win(ia) <= nmatup)
        call transDens(ii, jj, iAt, iAtomStart, nOrb, updwn, stimc, grndEigVecs, qq_ij)
        if (spin) then
          if (.not. updwn) then
            sindx = [2, 1]
          else
            sindx = [1, 2]
          end if
          do iSpin = 1, 2
            otmp(:nOrb, :nOrb, iAt, iSpin) = otmp(:nOrb, :nOrb, iAt, iSpin) + qq_ij(:nOrb, :nOrb)&
                & * wnij(ia) * vin(ia) * onsite(:nOrb, :nOrb, sindx(iSpin))
          end do
        else
          ! closed shell
          otmp(:nOrb, :nOrb, iAt, 1) = otmp(:nOrb, :nOrb, iAt, 1) + 0.5_dp &
              & * qq_ij(:nOrb, :nOrb) * wnij(ia) * vin(ia)&
              & * (onsite(:nOrb, :nOrb, 1) + fact * onsite(:nOrb, :nOrb, 2))
        end if
      end do
      ! rotational invariance corection for diagonal part
      do iSpin = 1, ss
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
      call indxov(win, ia, getij, ii, jj)
      updwn = (win(ia) <= nmatup)
      ss = 1
      if (.not. updwn) then
        ss = 2
      end if
      do iAt = 1, nAtom
        nOrb = orb%nOrbAtom(iAt)
        call transDens(ii, jj, iAt, iAtomStart, nOrb, updwn, stimc, grndEigVecs, qq_ij)
        vout(ia) = vout(ia) + 4.0_dp * wnij(ia) *&
            & sum(qq_ij(:nOrb, :nOrb) * otmp(:nOrb, :nOrb, iAt, ss))
      end do
    end do

  end subroutine onsiteEner


  !> Multiplies the supermatrix (A+B) with a given vector.
  subroutine apbw(rkm1, rhs2, wij, nmat, natom, win, nmatup, getij, iAtomStart, stimc, grndEigVecs,&
      & gamma, transChrg)

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

    !> machinery for transition charges between single particle levels
    type(TTransCharges), intent(in) :: transChrg

    !> gamma matrix
    integer :: ia, ii, jj
    real(dp) :: tmp(natom), gtmp(natom), qij(natom)

    @:ASSERT(size(rkm1) == nmat)

    tmp(:) = 0.0_dp
    call transChrg%qMatVec(iAtomStart, stimc, grndEigVecs, getij, win, rhs2, tmp)

    call hemv(gtmp, gamma, tmp)

    rkm1(:) = 0.0_dp
    call transChrg%qVecMat(iAtomStart, stimc, grndEigVecs, getij, win, gTmp, rkm1)

    rkm1(:) = 4.0_dp * rkm1(:) + wij(:) * rhs2(:)

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

end module dftbp_linrespcommon
