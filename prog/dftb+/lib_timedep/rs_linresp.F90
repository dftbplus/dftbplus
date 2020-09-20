!--------------------------------------------------------------------------------------------------!
!  DFTB+: general package for performing fast atomistic simulations                                !
!  Copyright (C) 2006 - 2020  DFTB+ developers group                                               !
!                                                                                                  !
!  See the LICENSE file for terms of usage and distribution.                                       !
!--------------------------------------------------------------------------------------------------!

#:include 'common.fypp'

!> This module heavily relies on the linresp module by Dominguez et al.
!> and the rangesep module by Lusker et al.
!> Periodic systems and 3rd order calculations are not supported so far.
module dftbp_rs_linearresponse
  use dftbp_assert
  use dftbp_arpack
  use dftbp_commontypes
  use dftbp_slakoCont
  use dftbp_shortgamma
  use dftbp_accuracy
  use dftbp_constants, only: Hartree__eV, au__Debye
  use dftbp_nonscc, only: TNonSccDiff
  use dftbp_scc, only: TScc
  use dftbp_blasroutines
  use dftbp_eigensolver
  use dftbp_message
  use dftbp_taggedOutput, only: TTaggedWriter, tagLabels
  use dftbp_linresp
  use dftbp_linrespcommon
  use dftbp_rangeseparated
  use dftbp_sorting
  use dftbp_qm
  use dftbp_transcharges
  use dftbp_sk, only: rotateH0
  !!^-- Check: once done, see which modules are actually required
  implicit none

  private

  public :: linRespCalcExcitationsRS

  type :: rs_linresp
    !> type to hold data required only for rangesep linresp but not linresp itself
    !> Check: If they turn out to be feq, consider adding them to linresp itself

  end type rs_linresp

  interface incSize
    module procedure incSizeInt
    module procedure incSizeVec
    module procedure incSizeMat
  end interface incSize

  !> Names of output files
  character(*), parameter :: arpackOut = "ARPACK.DAT"
  character(*), parameter :: testArpackOut = "TEST_ARPACK.DAT"
  character(*), parameter :: transitionsOut = "TRA.DAT"
  character(*), parameter :: XplusYOut = "XplusY.DAT"
  character(*), parameter :: excitedCoefsOut = "COEF.DAT"
  character(*), parameter :: excitationsOut = "EXC.DAT"
  character(*), parameter :: transDipOut = "TDP.DAT"
  character(*), parameter :: transChargesOut = "ATQ.DAT"

  !> Tolerance for ARPACK solver.
  real(dp), parameter :: ARTOL = epsilon(1.0_rsp)

  !> Maximal allowed iteration in the ARPACK solver.
  integer, parameter :: MAX_AR_ITER = 300

  character(lc) :: tmpStr

  !> Communication with ARPACK for progress information
  integer :: ndigit
  integer :: msaupd, msaup2, msaitr, mseigt, msapps, msgets, mseupd

contains

  !!> Computes compound occ-occ excitation index from individual indices.
!!! \param ii  Initial state.
!!! \param jj  Final state.
!!! \param index  Compound excitation index.
  !subroutine rIndXoo(nOcc, nOcc_r, ii, jj, index)
  !  integer, intent(in) ::  nOcc, nOcc_r, ii, jj
  !  integer, intent(out) :: index

  !  if (jj > ii) then
  !     index = ii - nOcc + nOcc_r + ((jj - nOcc + nOcc_r - 1) * (jj - nOcc + nOcc_r)) / 2
  !  else
  !     index = jj - nOcc + nOcc_r + ((ii - nOcc + nOcc_r - 1) * (ii - nOcc + nOcc_r)) / 2
  !  end if

  !end subroutine rIndXoo

  !> Computes the occupied-occupied compound index from two individual indices; assume jj <= ii
  pure subroutine rIndXoo(nOcc, nOcc_r, ii, jj, indX)

    !> Number of occupied states.
    integer, intent(in) :: nOcc

    !> Number of "active" occupied states.
    integer, intent(in) :: nOcc_r

    !> First (initial) occupied state.
    integer, intent(in) :: ii

    !> Second (final) occupied state.
    integer, intent(in) :: jj

    !> Compound indx.
    integer, intent(out) :: indX

    indX = jj - nOcc + nOcc_r + ((ii - nOcc + nOcc_r - 1) * (ii - nOcc + nOcc_r)) / 2

  end subroutine rIndXoo

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
    vect(:) = vect(:) / dist
    call getSKIntegrals(skOverCont, interSKOver, dist, iSp1, iSp2)
    call rotateH0(Sblock, interSKOver, vect(1), vect(2), vect(3), iSp1, iSp2, orb)

  end subroutine getSOffsite


  !> Set value of sizeIn to new, three times as large, value
  pure subroutine incSizeInt(sizeIn)

    integer, intent(inout) :: sizeIn
    sizeIn = 3 * sizeIn

  end subroutine incSizeInt


  !> increase size of a NxsizeIn array to 3*sizeIn
  pure subroutine incSizeMat(sizeIn, vec)

    integer, intent(inout) :: sizeIn
    real(dp), allocatable, intent(inout) :: vec(:,:)

    integer :: dim1
    real(dp), allocatable :: temp(:,:)
    dim1 = size(vec, dim=1)
    allocate(temp(dim1, 3 * sizeIn))
    ! Check: would be nice if temp could become memory for vec, see if possible in fortran
    temp(:,1:sizeIn) = vec
    call move_alloc(temp, vec)

  end subroutine incSizeMat


  !>increase size of a sizeIn array to 3*sizeIn
  pure subroutine incSizeVec(sizeIn, vec)

    integer, intent(in) :: sizeIn
    real(dp), allocatable, intent(inout) :: vec(:)

    real(dp), allocatable :: temp(:)
    allocate(temp(3 * sizeIn))
    ! Check: would be nice if temp could become memory for vec, see if possible in fortran
    temp(1:sizeIn) = vec
    call move_alloc(temp, vec)

  end subroutine incSizeVec


  !>increase size of a sizeInxN array to 3*sizeInxN
  pure subroutine incSizeMatSwapped(sizeIn, vec)

    integer, intent(inout) :: sizeIn
    real(dp), allocatable, intent(inout) :: vec(:,:)

    integer :: dim1
    real(dp), allocatable :: temp(:,:)
    dim1 = size(vec, dim=2)
    allocate(temp(3 * sizeIn, dim1))
    ! Check: would be nice if temo could become memory for vec, see if possible in fortran
    temp(1:sizeIn,:) = vec
    call move_alloc(temp, vec)

  end subroutine incSizeMatSwapped


  !>increase size of a sizeInxsizeIn array 
  pure subroutine incSizeMatBoth(sizeIn, vec)

    integer, intent(inout) :: sizeIn
    real(dp), allocatable, intent(inout) :: vec(:,:)

    real(dp), allocatable :: temp(:,:)

    allocate(temp(3 * sizeIn, 2 * sizeIn))
    ! Check: would be nice if temo could become memory for vec, see if possible in fortran
    temp(1:sizeIn, 1:sizeIn) = vec
    call move_alloc(temp, vec)

  end subroutine incSizeMatBoth


  !> Calculate square root and inverse of sqrt of a real, symmetric positive definite matrix
  subroutine calcMatrixSqrt(matIn, spaceDim, memDim, workArray, workDim, matOut, matInvOut)

    real(dp), intent(in) :: matIn(:,:)
    integer, intent(in) :: spaceDim, memDim
    real(dp), intent(out) :: matOut(:,:), matInvOut(:,:)
    real(dp) :: workArray(:)
    integer :: workDim

    real(dp) :: dummyEV(spaceDim)
    real(dp) :: dummyM(spaceDim, spaceDim), dummyM2(spaceDim, spaceDim)
    integer :: info
    integer :: ii

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
  !> vectors 1:(start-1) are orthonormal yet
  pure subroutine orthonormalizeVectors(start, end, vec)
    integer, intent(in) :: start
    integer, intent(in) :: end
    real(dp), intent(inout) :: vec(:,:)

    integer :: ii, jj

    do ii = start, end
      do jj = 1, ii - 1
        vec(:,ii) = vec(:,ii) - dot_product(vec(:,ii), vec(:,jj)) * vec(:,jj)
      end do
      vec(:,ii) = vec(:,ii) / sqrt(dot_product(vec(:,ii), vec(:,ii)))
    end do

  end subroutine orthonormalizeVectors

  !> Calculate the product of the matrix A+B and a vector. A,B as usual in linear response TD-DFT.
  subroutine multApBVecFast(vin, wIJ, sym, win, nMatUp, homo, nOcc, nVir, occNr, getIJ, gamma,&
      & lrGamma, species0, spinW, iaTrans, gqvTmp, tQov, tQoo, tQvv, vOut)
    real(dp), intent(in) :: vin(:)
    real(dp), intent(in) :: wIJ(:)
    character, intent(in) :: sym
    integer, intent(in) :: win(:), nMatUp, homo, nOcc, nVir !, ind(:)
    real(dp), intent(in) :: occNr(:,:)
    real(dp), intent(in) :: gamma(:,:), lrGamma(:,:)
    integer,  intent(in) :: getIJ(:,:)
    integer, intent(in) :: species0(:)
    real(dp), intent(in) :: spinW(:)
    real(dp), intent(in) :: tQov(:,:), tQoo(:,:), tQvv(:,:)
    integer, intent(in) :: iaTrans(1:,homo+1:)
    real(dp), intent(out) :: vOut(:)
    real(dp), intent(out) :: gqvTmp(:,:)   !for faster, more memory intensive routine

    integer :: izpAlpha, nMat, nAtom
    integer :: ia, ib, ja, jb, alpha, ii, jj, ij, aa, bb, ab
    real(dp) :: tmp2
    real(dp) :: oTmp(size(gamma, dim=1)), gTmp(size(gamma, dim=1))
    real(dp) :: wnIJ(size(wIJ))

    nMat = size(vin) ! also known as nXov
    nAtom = size(gamma, dim=1)
    gTmp(:) = 0.0_dp  !used to be oTmp before optimization

    call wtdn(wIJ, occNr, win, nMatUp, nMat, getIJ, wnIJ)

    !----only spin unpolarized case for now-------------------------------

    if (sym == 'S') then
      !full range coupling matrix contribution
      do ia = 1, nMat
        !call hemv(gTmp, gamma, tQov(:,ia))
        !oTmp(:) = oTmp + vin(ia) * gTmp
        gTmp(:) = gTmp + tQov(:,ia) * vin(ia) !gTmp=sum_jb q_B^jb * vin_jb
      end do

      call hemv(oTmp, gamma, gTmp)

      do ia= 1, nMat
        !call indXov(win, ia, getIJ, ii, jj)
        !lUpDwn = (win(ia) <= nMatUp)
        !qIJ = transQ(ii, jj, ind, lUpDwn, sTimesGrndEigVecs, grndEigVecs)
        vOut(ia) = 4.0 * dot_product(tQov(:,ia), oTmp)
        ! Check: prefactor! need 4.0 * for symmetry reduced sum?
      end do

    else
      do ia = 1, nMat
        !call indXov(win, ia, getIJ, ii, jj)
        !lUpDwn = (win(ia) <= nMatUp)
        !qIJ = transQ(ii, jj, ind, lUpDwn, sTimesGrndEigVecs, grndEigVecs)
        oTmp(:) = oTmp + vin(ia) * tQov(:,ia)
      end do

      do ia = 1, nMat
        vOut(ia) = 0.0_dp
        !call indXov(win, ia, getIJ, ii, jj)
        !lUpDwn = (win(ia) <= nMatUp)
        !qIJ = transQ(ii, jj, ind, lUpDwn, sTimesGrndEigVecs, grndEigVecs)
        do alpha = 1, nAtom
          izpAlpha = species0(alpha)
          tmp2 = 2.0_dp * (spinW(izpAlpha)) !== 2*magnetization
          vOut(ia) = vOut(ia) + oTmp(alpha) * tmp2 * tQov(alpha,ia)
        end do
      end do

    end if

    !   Old implementation, probably slower, but need less memory
    !   Reintroduce as option later on

    !!long range coupling matrix contribution
    !do ia = 1, nMat !adds -sum_AB ( q_ij^A * q_ab^B * Gamma_lr_AB + q_ib^A * q_ja^B * Gamma_lr_AB )
    !  do jb = 1, nMat

    !    call indXov(win, ia, getIJ, ii, aa)
    !    call indXov(win, jb, getIJ, jj, bb)

    !    lUpDwn = (win(iaTrans(jj, aa)) <= nMatUp) ! UNTESTED
    !    qIJ = transQ(jj, aa, ind, lUpDwn, sTimesGrndEigVecs, grndEigVecs)
    !    call hemv(gTmp, lrGamma, qIJ)

    !    lUpDwn = (win(iaTrans(ii, bb)) <= nMatUp) ! UNTESTED
    !    qIJ = transQ(ii, bb, ind, lUpDwn, sTimesGrndEigVecs, grndEigVecs)
    !    vOut(ia) = vOut(ia) - dot_product(qIJ, gTmp) * vin(jb)

    !    lUpDwn = .true. ! ???
    !    qIJ = transQ(ii, jj, ind, lUpDwn, sTimesGrndEigVecs, grndEigVecs)
    !    call hemv(gTmp, lrGamma, qIJ)

    !    lUpDwn = .true. ! ???
    !    qIJ = transQ(aa, bb, ind, lUpDwn, sTimesGrndEigVecs, grndEigVecs)
    !    vOut(ia) = vOut(ia) - dot_product(qIJ, gTmp) * vin(jb)

    !  end do
    !end do

    ! End old implementation

    ! Faster but more memory intensive routine

    gqvTmp(:,:) = 0.0_dp
    do jb = 1, nMat
      call indXov(win, jb, getIJ, jj, bb)
      !lUpDwn = (win(jb) <= nMatUp)
      do aa = homo + 1, homo + nVir
        !lUpDwn = .true.
        !qIJ = transQ(aa, bb, ind, lUpDwn, sTimesGrndEigVecs, grndEigVecs)
        if (bb < aa) then
          call rIndXvv(homo, aa, bb, ab)
        else
          call rIndXvv(homo, bb, aa, ab)
        end if
        !call hemv(gTmp, lrGamma, tQvv(:,ab))
        ja = nVir * (jj - 1) + aa - nOcc
        gqvTmp(:,ja) = gqvTmp(:,ja) + tQvv(:,ab) * vin(jb)

      end do
    end do

    do ja = 1, nMat
      gTmp(:) = gqvTmp(:,ja)
      call hemv(gqvTmp(:,ja), lrGamma, gTmp)
    end do

    do jj = 1, nOcc
      do ia = 1, nMat
        !would prefer to swap loops, but want to avoid calculation ia as unsure how implemented
        call indXov(win, ia, getIJ, ii, aa)
        !lUpDwn = (win(ia) <= nMatUp)
        !lUpDwn = .true.
        !qIJ = transQ(ii, jj, ind, lUpDwn, sTimesGrndEigVecs, grndEigVecs)
        if (jj < ii) then
          ! inverted indXoo function
          ij = jj + ((ii - 1 - homo + nOcc) * (ii -homo + nOcc)) / 2 + homo - nOcc
        else
          ij = ii + ((jj - 1 - homo + nOcc) * (jj -homo + nOcc)) / 2 + homo - nOcc
        end if
        ja = nVir * (jj - 1) + aa - nOcc
        vOut(ia) = vOut(ia) - dot_product(tQoo(:,ij), gqvTmp(:,ja))

      end do
    end do

    gqvTmp(:,:) = 0.0_dp
    do jb = 1, nMat
      call indXov(win, jb, getIJ, jj, bb)
      do aa = homo + 1, homo + nVir
        ab = (bb - nOcc - 1) * nVir + aa - nOcc
        ja = iaTrans(jj, aa)
        gqvTmp(:,ab) = gqvTmp(:,ab) + tQov(:,ja) * vin(jb)
      end do
    end do

    do ab = 1, nVir * nVir
      gTmp(:) = gqvTmp(:,ab)
      call hemv(gqvTmp(:,ab), lrGamma, gTmp)
    end do

    do bb = homo + 1, homo + nVir
      do ia = 1, nMat
        call indXov(win, ia, getIJ, ii, aa)
        ib = iaTrans(ii,bb)
        ab = (bb - nOcc - 1) * nVir + aa - nOcc
        vOut(ia) = vOut(ia) - dot_product(tQov(:,ib), gqvTmp(:,ab))
      end do
    end do

    ! End faster implementation

    !----only spin unpolarized case for now------------------------------------
    !vOut(:) = vOut + 0.5_dp * wnIJ * vin !orb. energy difference diagonal contribution
    vOut(:) = vOut + 0.5_dp * wnIJ(1:nMat) * vin(1:nMat)
    !Check: What about the factor 0.5 introduced here? Check derivation!

  end subroutine multApBVecFast


  !> Calculate the product of A-B and a vector.
  subroutine multAmBVecFast(vin, wIJ, win, nMatUp, homo, nOcc, nVir,&
      & occNr, getIJ, gamma, lrGamma, iaTrans, gqvTmp, tQov, tQoo, tQvv, vOut)
    !logical, intent(in) :: spin !for now ignore spin and assume unpolarized system
    real(dp), intent(in) :: vin(:)
    real(dp), intent(in) :: wIJ(:)
    !character, intent(in) :: sym
    integer, intent(in) :: win(:), nMatUp, homo, nOcc, nVir !, ind(:)
    !real(dp), intent(in) :: sTimesGrndEigVecs(:,:,:), grndEigVecs(:,:,:)
    real(dp), intent(in) :: occNr(:,:)
    real(dp), intent(in) :: gamma(:,:), lrGamma(:,:)
    integer,  intent(in) :: getIJ(:,:)
    !integer, intent(in) :: species0(:)
    !type(TOrbitals), intent(in) :: orb
    real(dp), intent(in) :: tQov(:,:), tQoo(:,:), tQvv(:,:)
    integer, intent(in) :: iaTrans(1:,homo+1:)
    real(dp), intent(out) :: vOut(:)

    real(dp) :: gqvTmp(:,:)   !for faster, more memory intensive routine
    integer :: nMat, nAtom
    integer :: ia, ib, ja, jb, ii, jj, ij, aa, bb, ab
    !real(dp) :: oTmp(size(gamma, dim=1))
    real(dp) :: gTmp(size(gamma, dim=1))
    !real(dp) :: qIJ(size(gamma, dim=1))
    real(dp) :: wnIJ(size(wIJ))
    !logical :: lUpDwn

    nMat = size(vin) ! also known as nXov
    nAtom = size(gamma, dim=1)

    call wtdn(wIJ, occNr, win, nMatUp, nMat, getIJ, wnIJ)

    !----only spin unpolarized case and singlet for now-------------------------------

    ! Old routine, reintroduce as option later on, see A+B fct

    !! Long range coupling matrix contribution
    !do ia = 1, nMat !Calc -sum_AB ( q_ij^A * q_ab^B * Gamma_lr_AB - q_ib^A * q_ja^B * Gamma_lr_AB )
    !  vOut(ia) = 0.0_dp
    !  do jb = 1, nMat

    !    call indXov(win, ia, getIJ, ii, aa)
    !    call indXov(win, jb, getIJ, jj, bb)

    !    lUpDwn = (win(iaTrans(jj, aa)) <= nMatUp) ! UNTESTED
    !    qIJ = transQ(jj, aa, ind, lUpDwn, sTimesGrndEigVecs, grndEigVecs)
    !    call hemv(gTmp, lrGamma, qIJ)

    !    lUpDwn = (win(iaTrans(ii, bb)) <= nMatUp) ! UNTESTED
    !    qIJ = transQ(ii, bb, ind, lUpDwn, sTimesGrndEigVecs, grndEigVecs)
    !    vOut(ia) = vOut(ia) + dot_product(qIJ, gTmp) * vin(jb)

    !    lUpDwn = .true.
    !    qIJ = transQ(ii, jj, ind, lUpDwn, sTimesGrndEigVecs, grndEigVecs)
    !    call hemv(gTmp, lrGamma, qIJ)

    !    lUpDwn = .true.
    !    qIJ = transQ(aa, bb, ind, lUpDwn, sTimesGrndEigVecs, grndEigVecs)
    !    vOut(ia) = vOut(ia) - dot_product(qIJ, gTmp) * vin(jb)

    !  end do
    !end do

    ! End old routine

    ! Faster but more memory intensive routine

    gqvTmp(:,:) = 0.0_dp
    do jb = 1, nMat
      call indXov(win, jb, getIJ, jj, bb)
      !lUpDwn = (win(jb) <= nMatUp)
      do aa = homo + 1, homo + nVir
        !lUpDwn = .true.
        !qIJ = transQ(aa, bb, ind, lUpDwn, sTimesGrndEigVecs, grndEigVecs)
        if (bb < aa) then
          call rIndXvv(homo, aa, bb, ab)
        else
          call rIndXvv(homo, bb, aa, ab)
        end if
        !call hemv(gTmp, lrGamma, tQvv(:,ab))
        ja = nVir * (jj - 1) + aa - nOcc
        gqvTmp(:,ja) = gqvTmp(:,ja) + tQvv(:,ab) * vin(jb)
      end do
    end do

    do ja = 1, nMat
      gTmp(:) = gqvTmp(:,ja)
      call hemv(gqvTmp(:,ja), lrGamma, gTmp)
    end do

    do jj = 1, nOcc
      do ia = 1, nMat
        ! would prefer to swap loops, but want to avoid calculation ia as unsure how implemted
        call indXov(win, ia, getIJ, ii, aa)
        !lUpDwn = (win(ia) <= nMatUp)
        !lUpDwn = .true.
        !qIJ = transQ(ii, jj, ind, lUpDwn, sTimesGrndEigVecs, grndEigVecs)
        if (jj < ii) then
          ! inverted indXoo function
          ij = jj + ((ii - 1 - homo + nOcc) * (ii - homo + nOcc)) / 2 + homo - nOcc
        else
          ij = ii + ((jj - 1 - homo + nOcc) * (jj - homo + nOcc)) / 2 + homo - nOcc
        end if
        ja = nVir * (jj - 1) + aa - nOcc
        vOut(ia) = vOut(ia) - dot_product(tQoo(:,ij), gqvTmp(:,ja))

      end do
    end do

    gqvTmp(:,:) = 0.0_dp
    do jb = 1, nMat
      call indXov(win, jb, getIJ, jj, bb)
      do aa = homo + 1, homo + nVir
        ab = (bb - nOcc - 1) * nVir + aa - nOcc
        ja = iaTrans(jj, aa)
        gqvTmp(:,ab) = gqvTmp(:,ab) + tQov(:,ja) * vin(jb)
      end do
    end do

    do ab = 1, nVir * nVir
      gTmp(:) = gqvTmp(:,ab)
      call hemv(gqvTmp(:,ab), lrGamma, gTmp)
    end do

    do bb = homo + 1, homo + nVir
      do ia = 1, nMat
        call indXov(win, ia, getIJ, ii, aa)
        ib = iaTrans(ii,bb)
        ab = (bb - nOcc - 1) * nVir + aa - nOcc
        vOut(ia) = vOut(ia) + dot_product(tQov(:,ib), gqvTmp(:,ab))
      end do
    end do

    ! End faster implementation

    ! orb. energy difference diagonal contribution
    vOut(:) = vOut + 0.5_dp * wnIJ(1:nMat) * vin(1:nMat)


  end subroutine multAmBVecFast


  !> Generates initial matrices Mm and Mp without calculating full Mat Vec product,
  !> not required for init. space.
  !> Use precalculated transition charges
  subroutine setupInitMatFast(initDim, wIJ, sym, win, nMatUp, nOcc, homo, occNr, getIJ, iaTrans,&
      & gamma, lrGamma, species0, spinW, tQov, tQoo, tQvv, vP, vM, mP, mM)
    integer, intent(in) :: initDim
    !logical, intent(in) :: spin !for now ignore spin and assume unpolarized system
    real(dp), intent(in) :: wIJ(:)
    character, intent(in) :: sym
    integer, intent(in) :: win(:), nMatUp, nOcc, homo ! , ind(:)
    !real(dp), intent(in) :: sTimesGrndEigVecs(:,:,:), grndEigVecs(:,:,:)
    real(dp), intent(in) :: occNr(:,:)
    real(dp), intent(in) :: gamma(:,:), lrGamma(:,:)
    integer,  intent(in) :: getIJ(:,:)
    integer, intent(in) :: iaTrans(1:, homo+1:)
    integer, intent(in) :: species0(:)
    real(dp), intent(in) :: spinW(:)
    !type(TOrbitals), intent(in) :: orb
    real(dp), intent(in) :: tQov(:,:), tQoo(:,:), tQvv(:,:)
    real(dp), intent(out) :: vP(:,:), vM(:,:)
    real(dp), intent(out) :: mP(:,:), mM(:,:)

    integer :: izpAlpha, nMat, nAtom, nw
    integer :: ia, jb, alpha, ii, jj, aa, bb, aaOld, bbOld, jjOld, ij, ab
    real(dp) :: tmp2
    real(dp), allocatable :: oTmp(:), oTmp2(:), qIJ(:), wnIJ(:)
    real(dp) :: dtmp, dtmp2

    nMat = size(vP, dim=1) ! also known as nXov
    nAtom = size(gamma, dim=1)
    nw = size(wIJ)

    allocate(oTmp(nAtom))
    allocate(oTmp2(nAtom))
    allocate(qIJ(nAtom))
    allocate(wnIJ(nw))

    call wtdn(wIJ, occNr, win, nMatUp, nMat, getIJ, wnIJ)

    !calc. vP, vM

    vP(:,:) = 0.0_dp
    vM(:,:) = 0.0_dp
    mP(:,:) = 0.0_dp
    mM(:,:) = 0.0_dp

    oTmp(:) = 0.0_dp

    !----only spin unpolarized case for now-------------------------------

    if (sym == 'S') then
      ! full range coupling matrix contribution
      do jb = 1, initDim ! calc 4* sum_A q^ia_A sum_B gamma_AB q^jb_B
        call hemv(oTmp, gamma, tQov(:,jb))
        do ia = 1, nMat
          vP(ia,jb) = 4.0 * dot_product(tQov(:,ia), oTmp)
        end do
      end do
    else
      do jb = 1, initDim
        oTmp(:) = oTmp + tQov(:,jb)
        do ia = 1, nMat
          vP(ia,jb) = 0.0_dp
          do alpha = 1, nAtom
            izpAlpha = species0(alpha)
            tmp2 = 2.0_dp * (spinW(izpAlpha)) !== 2*magnetization
            vP(ia,jb) = vP(ia,jb) + oTmp(alpha) * tmp2 * tQov(alpha, jb)
          end do
        end do
      end do
    end if

    aaOld = -1
    bbOld = -1
    jjOld = -1

    ! long range coupling matrix contribution
    do jb = 1, initDim
      call indXov(win, jb, getIJ, jj, bb)
      do ia = 1, nMat
        call indXov(win, ia, getIJ, ii, aa)
        if ((aaOld .ne. aa) .or. (bbOld .ne. bb)) then
          if (bb < aa) then
            call rIndXvv(homo, aa, bb, ab)
          else
            call rIndXvv(homo, bb, aa, ab)
          end if
          call hemv(oTmp, lrGamma, tQvv(:,ab))
        end if

        if (jj < ii) then
          ij = jj + ((ii - 1 - homo + nOcc) * (ii - homo + nOcc)) / 2 + homo - nOcc
        else
          ij = ii + ((jj - 1 - homo + nOcc) * (jj - homo + nOcc)) / 2 + homo - nOcc
        end if

        dtmp = dot_product(tQoo(:,ij), oTmp)

        vP(ia,jb) = vP(ia,jb) - dtmp
        vM(ia,jb) = vM(ia,jb) - dtmp

        if ((jjOld .ne. jj) .or. (aaOld .ne. aa)) then
          call hemv(oTmp2, lrGamma, tQov(:,iaTrans(jj,aa)))
        end if

        dtmp2 = dot_product(tQov(:,iaTrans(ii,bb)), oTmp2)

        vP(ia,jb) = vP(ia,jb) - dtmp2
        vM(ia,jb) = vM(ia,jb) + dtmp2

        aaOld = aa
      end do
      jjOld = jj
      bbOld = bb
    end do

    do jb = 1, initDim
      vP(jb,jb) = vP(jb,jb) + 0.5_dp * wnIJ(jb)  !orb. energy difference diagonal contribution
      vM(jb,jb) = vM(jb,jb) + 0.5_dp * wnIJ(jb)
    end do
    ! end calc. vP, vM

    ! Calc. matrices
    do ii = 1, initDim
      do jj = ii, initDim
        mP(ii,jj) = vP(ii,jj)
        mP(jj,ii) = mP(ii,jj)
        mM(ii,jj) = vM(ii,jj)
        mM(jj,ii) = mM(ii,jj)
      end do
    end do

    deallocate(oTmp)
    deallocate(oTmp2)
    deallocate(qIJ)
    deallocate(wnIJ)

  end subroutine setupInitMatFast


  !> Run Linear Response calc. The actual diagonalization of the RPA equation
  !> is outsourced to rsLinRespCalc. Most code copied from non range-separated
  !> version by A. Dominguez.
  subroutine runRsLinRespCalc(spin, tOnsite, nAtom, iAtomStart, grndEigVecs, grndEigVal, sccCalc,&
      & dQ, coord0, nExc, nStat0, cSym, SSqr, filling, species0, nBeweg, HubbardU, spinW,&
      & rNel, iNeighbor, img2CentCell, orb, rsData, tWriteTagged, fdTagged, taggedWriter,&
      & fdMulliken, fdCoeffs, fdXplusY, fdTrans, fdSPTrans, fdTraDip, fdTransQ, tArnoldi,&
      & fdArnoldi, fdExc, tEnergyWindow, energyWindow, tOscillatorWindow, oscillatorWindow,&
      & tCacheCharges, omega, shift, skHamCont, skOverCont, derivator, deltaRho, excGrad, dQAtomEx)
    logical, intent(in) :: spin
    logical, intent(in) :: tOnsite
    integer, intent(in) :: nAtom, iAtomStart(:)
    real(dp), intent(in) :: grndEigVecs(:,:,:), grndEigVal(:,:), dQ(:), coord0(:,:)
    type(TScc), intent(in) :: sccCalc
    integer, intent(in) :: nExc, nStat0
    character, intent(in) :: cSym
    real(dp), intent(in) :: SSqr(:,:)
    real(dp), intent(in) :: filling(:,:)
    integer, intent(in) :: species0(:), nBeweg
    real(dp), intent(in) :: HubbardU(:), spinW(:), rNel ! maybe needed for gradients later
    integer, intent(in) :: iNeighbor(0:,:), img2CentCell(:)
    type(TOrbitals), intent(in) :: orb
    type(TRangeSepFunc), intent(inout) :: rsData ! contains long-range gamma
    !> print tag information
    logical, intent(in) :: tWriteTagged

    !> file id for tagging information
    integer, intent(in) :: fdTagged

    !> tagged writer
    type(TTaggedWriter), intent(inout) :: taggedWriter

    !logical :: tMulliken, tCoeffs, tXplusY, tTrans, tTraDip, tArnoldi

    !> file unit for excited Mulliken populations?
    integer, intent(in) :: fdMulliken
    !> file unit if the coefficients for the excited states should be written to disc
    integer, intent(in) :: fdCoeffs
    !> file for X+Y data
    integer, intent(in) :: fdXplusY
    !> File unit for single particle (KS) transitions if required
    integer, intent(in) :: fdTrans
    !> File unit for single particle transition dipole strengths
    integer, intent(in) :: fdSPTrans
    !> File unit for transition dipole data
    integer, intent(in) :: fdTraDip
    !> File unit for transition charges
    integer, intent(in) :: fdTransQ
    !> write state of Arnoldi solver to disc
    logical, intent(in) :: tArnoldi
    !> file unit for Arnoldi write out
    integer, intent(in) :: fdArnoldi
    !> file handle for excitation energies
    integer, intent(in) :: fdExc

    !real(dp), intent(in) :: ons_en(:,:), ons_dip(:,:)
    logical, intent(in) :: tEnergyWindow, tOscillatorWindow
    real(dp), intent(in) :: energyWindow, oscillatorWindow
    logical, intent(in) :: tCacheCharges
    real(dp), intent(out) :: omega
    real(dp), intent(in), optional :: shift(:)
    real(dp), intent(inout), optional :: deltaRho(:,:)
    type(TSlakoCont), intent(in), optional :: skHamCont, skOverCont
    class(TNonSccDiff), intent(in), optional :: derivator
    real(dp), intent(inout), optional :: excGrad(:,:)   !also may be needed for gradients later on
    real(dp), intent(inout), optional :: dQAtomEx(:)

    real(dp) :: Ssq(nExc)
    real(dp), allocatable :: gamma(:,:), lrGamma(:,:), snglPartTransDip(:,:)
    real(dp), allocatable :: sTimesGrndEigVecs(:,:,:), wIJ(:)
    real(dp), allocatable :: snglPartOscStrength(:), oscStrength(:), vecXmY(:), vecXpY(:), pc(:,:)
    real(dp), allocatable :: t(:,:), rhs(:), vWoo(:), vWvv(:), vWov(:)
    real(dp), allocatable :: evec(:,:), eval(:), allXpY(:,:), transitionDipoles(:,:)
    real(dp), allocatable :: atomicTransQ(:)
    !to hold precalculated transition charges, alloc and calc in rs calc
    real(dp), allocatable :: tQov(:,:), tQoo(:,:), tQvv(:,:)
    integer, allocatable :: win(:), iaTrans(:,:), getIJ(:,:)
    character, allocatable :: symmetries(:)

    integer :: nStat, nOcc, nOccR, nVirR, nXooR, nXvvR
    integer :: nXov, nXovUD(2), nXovR, nXovD, nXovRD, nOrb
    integer :: i, j, isym
    integer :: spinDim
    character :: cSym2
    logical :: tGrads, tZVector
    real(dp) :: energyThreshold
    integer :: logfil
    !real(dp) ::  t0,t1,t2,t3,t4,t5

    !> spin-polarized calculation?
    logical :: tSpin
    integer :: nSpin, iSpin
    !> control variables
    logical :: tCoeffs, tTraDip, tTransQ, tTrans, tXplusY !, tZVector
    !> printing data
    logical :: tMulliken

    !> transition charges, either cached or evaluated on demand
    type(TTransCharges) :: transChrg

    !> aux variables
    integer :: mu, nu

    ! For now, spin-polarized calculation not supported
    tSpin = .false.
    nSpin = 1
    ! For now, ons not supported with range separation
    if (tOnsite) then
      call error("Onsite corrections not currently supported for range-separated excited state&
          & calculations")
    end if

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
    end if
    ! End of ARPACK communication variables

    spinDim = size(grndEigVal, dim=2)
    nOrb = orb%nOrb

    @:ASSERT(present(excGrad) .eqv. present(shift))
    @:ASSERT(present(shift) .eqv. present(skHamCont))
    @:ASSERT(present(skHamCont) .eqv. present(skOverCont))

    ! work out which data files are required, based on whether they have valid file IDs (>0)
    tMulliken = (fdMulliken > 0)
    tCoeffs = (fdCoeffs > 0)
    tTraDip = (fdTraDip > 0)
    tTransQ = (fdTransQ > 0)
    tTrans = (fdTrans > 0)
    tXplusY = (fdXplusY > 0)

    if (tMulliken) then
      open(fdMulliken, file=excitedQOut, position="rewind", status="replace")
      close(fdMulliken)
      open(fdMulliken, file=excitedDipoleOut, position="rewind", status="replace")
      close(fdMulliken)
    end if

    @:ASSERT(fdArnoldi > 0)
    if (tArnoldi) then
      open(fdArnoldi, file=arpackOut, position="rewind", status="replace")
    end if

    ! Transition charges
    if (tTransQ) then
      allocate(atomicTransQ(nAtom))
    end if

    ! count initial number of transitions from occupied to empty states
    nXovUD = 0
    do iSpin = 1, nSpin
      !$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(i,j) SCHEDULE(RUNTIME) REDUCTION(+:nXovUD)
      do i = 1, nOrb - 1
        do j = i, nOrb
          if (filling(i, iSpin) > filling(j, iSpin) + elecTolMax) then
            nXovUD(iSpin) = nXovUD(iSpin) + 1
          end if
        end do
      end do
      !$OMP  END PARALLEL DO
    end do
    nXov = sum(nXovUD)

    if (nExc + 1 >= nXov) then
      write(tmpStr,"(' Insufficent single particle excitations, ', I0,&
          & ', for required number of excited states ', I0)") nXov, nExc
      call error(tmpStr)
    end if

    ! are gradients required?
    tGrads = present(excGrad)  ! For now no gradients, but keep for later!
    ! is a Z-vector required?
    tZVector = tGrads .or. tMulliken .or. tCoeffs .or. tTransQ

    ! Sanity checks
    nStat = nStat0
    if (nStat < 0 .and. cSym /= "S") then
      call error("Linresp: Brightest mode only for singlets.")
    elseif (nStat /= 0 .and. cSym == "B") then
      call error("Linresp: Both symmetries not allowed if specific state is excited")
    elseif (nStat == 0 .and. tZVector) then
      call error("Linresp: Gradient, charges, coefficients and charges only with selected&
          & excitation.")
    elseif (tGrads .and. nExc > nXov) then
      call error("Linresp: With gradients nExc can be max. the number of occ-virt excitations")
    end if

    ! Select symmetries to process
    !ADG: temporal solution for spin case.
    if (.not. spin) then
      select case (cSym)
      case ("B")
        allocate(symmetries(2))
        symmetries(:) = ["T", "S"]
      case ("S")
        allocate(symmetries(1))
        symmetries(:) = ["S"]
      case ("T")
        allocate(symmetries(1))
        symmetries(:) = ["T"]
      end select
    else
      allocate(symmetries(1))
      symmetries(:) = [" "]
    end if

    ! Allocation for general arrays
    allocate(gamma(nAtom, nAtom))
    allocate(lrGamma(nAtom, nAtom))
    allocate(snglPartTransDip(nXov, 3))
    allocate(sTimesGrndEigVecs(nOrb, nOrb, spinDim))
    allocate(wIJ(nXov))
    allocate(win(nXov))
    allocate(eval(nExc))
    allocate(getIJ(nXov, 2))
    allocate(transitionDipoles(nXov, 3))
    allocate(snglPartOscStrength(nXov))

    ! Arrays for gradients and Mulliken analysis !Again not required for now
    if (tGrads .or. tMulliken) then
      allocate(pc(nOrb, nOrb))
    end if

    ! Overlap times wave function coefficients - most routines in DFTB+ use lower triangle
    ! (would remove the need to symmetrize the overlap and ground state density matrix
    ! in the main code if this could be used everywhere in these routines)
    do iSpin = 1, nSpin
      call symm(sTimesGrndEigVecs(:,:,iSpin), "L", SSqr, grndEigVecs(:,:,iSpin))
    end do

    ! ground state Hubbard U softened coulombic interactions
    call sccCalc%getAtomicGammaMatrix(gamma, iNeighbor, img2CentCell)
    ! gamma is symmetric but still the upper triangle is not enough
    do mu = 2, nAtom
      do nu = 1, mu-1
        gamma(nu, mu) = gamma(mu, nu)
      end do
    end do
    lrGamma(:,:) = 0._dp
    call rsData%getLrGamma(lrGamma)

    ! Oscillator strengths for excited states, when needed.
    !if (nStat <= 0) then
    allocate(oscStrength(nExc))
    !end if

    ! Find all single particle transitions and KS energy differences
    !   for cases that go from filled to empty states
    call getSPExcitations(grndEigVal, filling, wIJ, getIJ)

    ! put them in ascending energy order
    if (tOscillatorWindow) then
      ! use a stable sort so that degenerate transitions from the same single particle state are
      ! grouped together in the results, allowing these to be selected together (since how intensity
      ! is shared out over degenerate transitions is arbitrary between eigensolvers/platforms).
      call merge_sort(win, wIJ, 1.0_dp * epsilon(1.0))
    else
      ! do not require stability, use the usual routine to sort, saving an O(N) workspace
      call index_heap_sort(win, wIJ)
    end if
    wIJ = wIJ(win)

    ! dipole strength of transitions between K-S states
    call calcTransitionDipoles(coord0, win, nXovUD(1), getIJ, iAtomStart, sTimesGrndEigVecs,&
        & grndEigVecs, snglPartTransDip)

    ! single particle excitation oscillator strengths
    snglPartOscStrength(:) = twothird * wIJ(:) * sum(snglPartTransDip**2, dim=2)

    if (tOscillatorWindow .and. tZVector ) then
      call error("Incompabilitity between excited state property evaluation&
          & and an oscillator strength window at the moment.")
    end if

    if (tOscillatorWindow .or. tEnergyWindow) then
      if (.not. tEnergyWindow) then
        ! find transitions that are strongly dipole allowed (> oscillatorWindow)
        call dipselect(wIJ, snglPartOscStrength, win, snglPartTransDip, nXovRD, oscillatorWindow,&
            & grndEigVal, getIJ)
      else
        ! energy window above the lowest nExc single particle transitions
        energyThreshold = wIJ(nExc) + energyWindow
        nXovR = count(wIJ <= energyThreshold)
        nXovD = 0
        if (tOscillatorWindow) then
          ! find transitions that are strongly dipole allowed (> oscillatorWindow)
          if (nXovR < nXov) then
            ! find transitions that are strongly dipole allowed (> oscillatorWindow)
            call dipselect(wIJ(nXovR+1:), snglPartOscStrength(nXovR+1:), win(nXovR+1:),&
                & snglPartTransDip(nXovR+1:,:), nXovD, oscillatorWindow, grndEigVal, getIJ)
          end if
        end if
        nXovRD = nXovR + nXovD
      end if
    else
      nXovRD = nXov
    end if

    ! just in case energy/dipole windows add no extra states, and is due to an arpack solver
    ! requirement combined with the need to get at least nExc states
    nXovRD = max(nXovRD, min(nExc+1, nXov))

    call TTransCharges_init(transChrg, iAtomStart, sTimesGrndEigVecs, grndEigVecs,&
        & nXovRD, nXovUD(1), getIJ, win, tCacheCharges)

    !if (nStat == 0) then
    !  if(tTrans) then
    !    call writeSPExcitations(wIJ, snglPartTransDip, win, nXovUD(1), getIJ,&
    !        & tWriteTagged, fdTagged, taggedWriter, snglPartOscStrength)
    !  end if
    !  call openExcitationFiles(spin, tXplusY, tTrans, tTraDip, tArnoldi, logfil)
    !end if

    if (tXplusY) then
      open(fdXplusY, file=XplusYOut, position="rewind", status="replace")
    end if

    if(tTrans) then
      open(fdTrans, file=transitionsOut, position="rewind", status="replace")
      write(fdTrans,*)
    end if

    ! single particle transition dipole file
    if (tTraDip) then
      open(fdTraDip, file=transDipOut, position="rewind", status="replace")
      write(fdTraDip,*)
      write(fdTraDip,'(5x,a,5x,a,2x,a)') "#", 'w [eV]', 'Transition dipole (x,y,z) [Debye]'
      write(fdTraDip,*)
      write(fdTraDip,'(1x,57("="))')
      write(fdTraDip,*)
    end if

    ! excitation energies
    open(fdExc, file=excitationsOut, position="rewind", status="replace")
    write(fdExc,*)
    if (tSpin) then
      write(fdExc,'(5x,a,7x,a,9x,a,9x,a,6x,a,4x,a)')&
          & 'w [eV]', 'Osc.Str.', 'Transition', 'Weight', 'KS [eV]', 'D<S*S>'
    else
      write(fdExc,'(5x,a,7x,a,9x,a,9x,a,6x,a,4x,a)')&
          & 'w [eV]', 'Osc.Str.', 'Transition', 'Weight', 'KS [eV]', 'Sym.'
    end if

    write(fdExc,*)
    write(fdExc,'(1x,80("="))')
    write(fdExc,*)

    ! single particle excitations (output file and tagged file if needed).  Was used for nXovRD =
    ! size(wIJ), but now for just states that are actually included in the excitation calculation.
    call writeSPExcitations(wIJ, win, nXovUD(1), getIJ, fdSPTrans, snglPartOscStrength, nXovRD,&
        & tSpin)

    ! redefine if needed (generalize it for spin-polarized and fractional occupancy)
    nOcc = nint(rNel) / 2
    nOccR = nOcc
    nVirR = nOrb - nOcc

    ! elements in a triangle plus the diagonal of the occ-occ and virt-virt blocks
    nXooR = (nOccR * (nOccR + 1)) / 2
    nXvvR = (nVirR * (nVirR + 1)) / 2

    allocate(evec(nXovRD, nExc))

    !allocate(allXpY(nXovUD(1), nExc)) ! allocate(allXpY(nXovRD, nExc))
    allocate(allXpY(nXovRD, nExc))

    if (nStat > 0) then
      !allocate(vecXmY(nXovUD(1)))
      !allocate(vecXpY(nXovUD(1)))
      allocate(vecXmY(nXovRD))
      allocate(vecXpY(nXovRD))
    end if

    ! Arrays needed for Z vector
    if (tZVector) then
      allocate(t(nOrb, nOrb))
      allocate(rhs(nXovRD))
      allocate(vWoo(nXooR))
      allocate(vWvv(nXvvR))
      allocate(vWov(nXovRD))
    end if

    ! for fast init. mat calc. need combined index
    allocate(iaTrans(1:nOcc, nOcc+1:nOrb))
    call rIndXov_array(win, nOcc, nXov, getIJ, iaTrans)

    ! run lin. resp. calculation:
    do isym = 1, size(symmetries)
      cSym2 = symmetries(isym)

      call rsLinRespCalc(tZVector, tTransQ, wIJ, nExc, cSym2, win, nXovUD(1), nXovRD, nOcc, nOccR,&
          & nVirR, iAtomStart, sTimesGrndEigVecs, grndEigVecs, filling, getIJ, iaTrans, gamma,&
          & lrGamma, species0, spinW, eval, evec, allXpY, nStat, vecXmY, tQov, tQoo, tQvv)

      call getOscillatorStrengthsRS(cSym2, snglPartTransDip(1:nXovRD,:), eval, allXpY, filling,&
          & nStat, oscStrength, tTraDip, transitionDipoles)

      if (spin) then
        call getExcSpin(Ssq, nXovUD(1), getIJ, win, eval, evec, wIJ, filling, sTimesGrndEigVecs,&
            & grndEigVecs)
        call writeExcitationsRS(cSym2, oscStrength, nExc, nXovUD(1), getIJ, win, eval, allXpY,&
            & wIJ(1:nXovRD), fdXplusY, fdTrans, fdTraDip, transitionDipoles, tWriteTagged,&
            & fdTagged, taggedWriter, fdExc, Ssq)
      else
        call writeExcitationsRS(cSym2, oscStrength, nExc, nXovUD(1), getIJ, win, eval, allXpY,&
            & wIJ(1:nXovRD), fdXplusY, fdTrans, fdTraDip, transitionDipoles, tWriteTagged,&
            & fdTagged, taggedWriter, fdExc)
      end if

      if (tTransQ) then
        call transitionChargesRS(allXpY(:,nStat), tQov, atomicTransQ)
        if (.not. tZVector) then
          deallocate(tQov)
        end if
        call writeTransitionChargesRS(fdTransQ, atomicTransQ)
      end if

    end do

    if (tArnoldi) then
      close(fdArnoldi)
    end if

    if (fdTrans > 0) close(fdTrans)
    if (fdXplusY > 0) close(fdXplusY)
    if (fdExc > 0) close(fdExc)
    if (fdTraDip > 0) close(fdTraDip)

    if (nStat == 0) then
      omega = 0.0_dp
      return
    end if

    ! Attention: right now I take sqrt in rsLinRespCalc. May not want to do this!!!
    omega = sqrt(eval(nStat))

    if (tZVector) then

      ! Furche terms: X+Y, X-Y
      !vecXmY(:) = sqrt(omega) / sqrt(wIJ) * evec(:,nStat)    ! vecXmY was set in rsLinRespCalc
      vecXpY(:) = allXpY(:,nStat)

      ! already performed above (for efficient init mat calc)
      !call rIndXov_array(win, nOcc, nXov, getIJ, iaTrans)

      call getZVectorEqRhsRS(vecXpY, vecXmY, win, iAtomStart, nOcc, nOccR, nXovUD(1), getIJ,&
          & iaTrans, nAtom, species0, grndEigVal(:,1), sTimesGrndEigVecs, grndEigVecs, gamma,&
          & lrGamma, spinW, omega, cSym2, rhs, t, vWov, vWoo, vWvv)

      call solveZVectorEqRS(rhs, win, nXovUD(1), getIJ, filling, nAtom, species0, spinW,&
          & gamma, lrGamma, wIJ, nOcc, nOccR, nVirR, iaTrans, tQov, tQoo, tQvv)

      call calcWVectorZRS(rhs, win, nOcc, nOrb, nXovUD(1), getIJ, iAtomStart, sTimesGrndEigVecs,&
          & grndEigVecs, gamma, lrGamma, grndEigVal(:,1), vWov, vWoo, vWvv, iaTrans)

      call calcPMatrix(t, rhs, win, getIJ, pc)

      call transformMO2AODense(pc, grndEigVecs(:,:,1))

      call getExcMulliken(iAtomStart, pc, SSqr, dQAtomEx)
      if (tMulliken) then
        call writeExcMulliken(cSym2, nStat, dQ, dQAtomEx, coord0, fdMulliken)
      end if

      if (tGrads) then
        call addGradientsRS(cSym2, nXovRD, nAtom, species0, iAtomStart, nOrb, nOcc, nXovUD(1),&
            & getIJ, win, grndEigVecs, pc, dQ, dQAtomEx, gamma, lrGamma, rsData, HubbardU,&
            & spinW, shift, vWoo, vWov, vWvv, vecXpY, vecXmY, nBeweg, coord0, orb,&
            & skHamCont, skOverCont, tQov, derivator, deltaRho, excGrad)
      end if
    end if
  end subroutine runRsLinRespCalc


  !> Perform linear response calculation with algorithm by Stratmann and Scuseria (JCP 1998)
  subroutine rsLinRespCalc(tZVector, tTransQ, wIJ, nExc, cSym, win, nMatUp, nXov, homo, nOcc,&
      & nVir, iAtomStart, sTimesGrndEigVecs, grndEigVecs, occNr, getIJ, iaTrans, gamma, lrGamma,&
      & species0, spinW, eval, evec, vecXpY, nStat, vecXmY, tQov, tQoo, tQvv)
    logical, intent(in) :: tZVector, tTransQ
    real(dp),intent(in) :: wIJ(:)
    character, intent(in) :: cSym
    integer, intent(in) :: nExc ! desired number of excitations
    integer, intent(in) :: win(:), nMatUp, nXov, homo, nOcc, nVir, iAtomStart(:)
    real(dp), intent(in) :: sTimesGrndEigVecs(:,:,:), grndEigVecs(:,:,:)
    real(dp), intent(in) :: occNr(:,:)
    real(dp), intent(in) :: gamma(:,:), lrGamma(:,:)
    integer, intent(in) :: getIJ(:,:), species0(:)
    integer, intent(in) :: iaTrans(1:,homo+1:) ! needed in fast setup of init mat
    real(dp), intent(in) :: spinW(:)
    integer, intent(in) :: nStat
    ! lin. resp. eigenvalues and eigenvectors (F_I and X_I+Y_I)
    real(dp), intent(out) :: eval(:), evec(:,:), vecXpY(:,:)
    real(dp), intent(out) :: vecXmY(:)
    ! pre-calculated transition charges
    real(dp), allocatable, intent(out) :: tQov(:,:), tQoo(:,:), tQvv(:,:)

    integer :: ii, jj, ia, ij, aa, bb, ab, nXovAct
    ! cur. size of subspace and memory for subspace vectors and matrices
    integer :: subSpaceDim, prevSubSpaceDim, memDim, workDim
    real(dp), allocatable :: vecB(:,:) ! basis of subspace
    real(dp), allocatable :: evecL(:,:), evecR(:,:) ! left and right eigenvectors of Mnh
    real(dp), allocatable :: vP(:,:), vM(:,:) ! vec. for (A+B)b_i, (A-B)b_i
    ! matrices M_plus, M_minus, M_minus^(1/2), M_minus^(-1/2) and M_herm~=resp. mat on subsapce
    real(dp), allocatable :: mP(:,:), mM(:,:), mMsqrt(:,:), mMsqrtInv(:,:), mH(:,:), vecXpYtest(:,:)
    real(dp), allocatable :: evalInt(:) ! store eigenvectors within routine
    real(dp), allocatable :: dummyM(:,:), workArray(:)
    real(dp), allocatable :: vecNorm(:) ! will hold norms of residual vectors
    real(dp), allocatable :: gqvTmp(:,:) ! work array used in matrix multiplications
    real(dp), allocatable :: qIJ(:)
    real(dp) :: dummyReal
    logical :: lUpdwn
    integer :: info
    integer :: dummyInt
    logical :: didConverge
    real(dp), parameter :: convThreshold = 0.0001_dp
    integer :: newVec
    integer :: initExc
    initExc = min(20 * nExc, nOcc * nVir)

    nXovAct = nXov ! or nOcc * nVir
    allocate(vecXpYtest(size(evec,dim=1), size(evec,dim=2)))

    ! initial allocations
    ! start with lowest excitations. Inital number somewhat arbritary.
    subSpaceDim = min(initExc, nXov)
    memDim = min(subSpaceDim + 6 * nExc, nXov) ! memory available for subspace calcs.
    workDim = 3 * memDim + 1
    allocate(vecB(nXovAct, memDim)) ! nXov, memDim))
    allocate(vP(nXovAct, memDim)) ! nXov, memDim))
    allocate(vM(nXovAct, memDim)) ! nXov, memDim))
    allocate(mP(memDim, memDim))
    allocate(mM(memDim, memDim))
    allocate(mMsqrt(memDim, memDim))
    allocate(mMsqrtInv(memDim, memDim))
    allocate(mH(memDim, memDim))
    allocate(dummyM(memDim, memDim))
    allocate(evalInt(memDim))
    allocate(evecL(memDim, nExc))
    allocate(evecR(memDim, nExc))
    allocate(workArray(3 * memDim + 1))
    allocate(vecNorm(2 * memDim))

    ! Work array for faster implementation of (A+-B)*v. Make optional later!
    ! for gwvTmp: could use symmetry to improve nVir * nVir case
    allocate(gqvTmp(size(gamma, dim=1), max(nVir, nOcc) * nVir))
    allocate(tQov(size(gamma, dim=1), nOcc * nVir))
    allocate(tQoo(size(gamma, dim=1), nOcc * (nOcc + 1) / 2))
    allocate(tQvv(size(gamma, dim=1), nVir * (nVir + 1) / 2))
    allocate(qIJ(size(gamma, dim=1)))

    ! Precalculate transition charges
    do ia = 1, nVir * nOcc
      call indXov(win, ia, getIJ, ii, aa)
      lUpdwn = (win(ia) <= nMatUp)
      qIJ = transQ(ii, aa, iAtomStart, lUpdwn, sTimesGrndEigVecs, grndEigVecs)
      tQov(:,ia) = qIJ
    end do
    do ij = 1, nOcc * (nOcc + 1) / 2
      call indXoo(ij, ii, jj)
      lUpdwn = .true. ! UNTESTED
      qIJ = transQ(ii, jj, iAtomStart, lUpdwn, sTimesGrndEigVecs, grndEigVecs)
      tQoo(:,ij) = qIJ
    end do
    do ab = 1, nVir * (nVir + 1) / 2
      call indXvv(homo, ab, aa, bb)
      lUpdwn = .true. ! UNTESTED
      qIJ = transQ(aa, bb, iAtomStart, lUpdwn, sTimesGrndEigVecs, grndEigVecs)
      tQvv(:,ab) = qIJ
    end do

    ! set initial bs
    vecB(:,:) = 0.0_dp
    do ii = 1, subSpaceDim
      vecB(ii, ii) = 1.0
    end do

    ! Could calc. init. matrix here, Due to special form of bs, matrix multiplication handles a lot
    ! of zeros.
    ! Add this later!

    prevSubSpaceDim = 0
    didConverge = .false.

    ! Solve the linear response problem. Iterative expansion of subspace:
    solveLinResp: do

      if (prevSubSpaceDim > 0) then

        !extend subspace matrices:
        do ii = prevSubSpaceDim + 1, subSpaceDim
          call multApBVecFast(vecB(:,ii), wIJ, cSym, win, nMatUp, homo, nOcc, nVir, occNr, getIJ,&
              & gamma, lrGamma, species0, spinW, iaTrans, gqvTmp, tQov, tQoo, tQvv,&
              & vP(:,ii))
          call multAmBVecFast(vecB(:,ii), wIJ, win, nMatUp, homo, nOcc, nVir, occNr, getIJ, gamma,&
              & lrGamma, iaTrans, gqvTmp, tQov, tQoo, tQvv, vM(:,ii))
        end do

        do ii = prevSubSpaceDim + 1, subSpaceDim
          do jj = 1, ii
            mP(ii,jj) = dot_product(vecB(:,jj), vP(:,ii))
            mP(jj,ii) = mP(ii,jj)
            mM(ii,jj) = dot_product(vecB(:,jj), vM(:,ii))
            mM(jj,ii) = mM(ii,jj)
          end do
        end do

      else

        call setupInitMatFast(subSpaceDim, wIJ, cSym, win, nMatUp, nOcc, homo, occNr, getIJ,&
            & iaTrans, gamma, lrGamma, species0, spinW, tQov, tQoo, tQvv, vP, vM, mP, mM)
      end if

      call calcMatrixSqrt(mM, subSpaceDim, memDim, workArray, workDim, mMsqrt, mMsqrtInv)
      call dsymm('L', 'U', subSpaceDim, subSpaceDim, 1.0_dp, mP, memDim, mMsqrt, memDim,&
          & 0.0_dp, dummyM, memDim)
      call dsymm('L', 'U', subSpaceDim, subSpaceDim, 1.0_dp, mMsqrt, memDim, dummyM, memDim,&
          & 0.0_dp, mH, memDim)

      ! diagonalize in subspace
      call dsyev('V', 'U', subSpaceDim, mH, memDim, evalInt, workArray, workDim, info)

      ! This yields T=(A-B)^(-1/2)|X+Y>.
      ! Calc. |R_n>=|X+Y>=(A-B)^(1/2)T and |L_n>=|X-Y>=(A-B)^(-1/2)T.
      ! Transformation preserves orthonormality.
      ! Only compute up to nExc index, because only that much needed.
      call dsymm('L', 'U', subSpaceDim, nExc, 1.0_dp, Mmsqrt, memDim, Mh, memDim, 0.0_dp,&
          & evecR, memDim)
      call dsymm('L', 'U', subSpaceDim, nExc, 1.0_dp, Mmsqrtinv, memDim, Mh, memDim, 0.0_dp,&
          & evecL, memDim)
      ! check if dsymm can actually be used

      ! Need |X-Y>=sqrt(w)(A-B)^(-1/2)T, |X+Y>=(A-B)^(1/2)T/sqrt(w) for proper solution to original
      ! EV problem
      do ii = 1, nExc ! only use first nExc vectors
        dummyReal = sqrt(sqrt(evalInt(ii)))
        evecR(:,ii) = evecR(:,ii) / dummyReal
        evecL(:,ii) = evecL(:,ii) * dummyReal
      end do

      !see if more memory is needed to save extended basis. If so increase amount of memory.
      if (subSpaceDim + 2 * nExc > memDim) then
        call incSize(memDim, vecB)
        call incSize(memDim, vP)
        call incSize(memDim, vM)
        call incSizeMatBoth(memDim, mP)
        call incSizeMatBoth(memDim, mM)
        call incSizeMatBoth(memDim, mH)
        call incSizeMatBoth(memDim, mMsqrt)
        call incSizeMatBoth(memDim, mMsqrtInv)
        call incSizeMatBoth(memDim, dummyM)
        call incSize(memDim, evalInt)
        call incSize(workDim, workArray)
        call incSizeMatSwapped(memDim, evecL)
        call incSizeMatSwapped(memDim, evecR)
        call incSize(2 * memDim, vecNorm)
        call incSize(memDim)
        call incSize(workDim)
      end if

      ! Calculate the residual vectors
      !   calcs. all |R_n>
      call dgemm('N', 'N', nXovAct, nExc, subSpaceDim, 1.0_dp, vecB, nXovAct, evecR, memDim,&
          & 0.0_dp, vecB(1,subSpaceDim+1), nXovAct)
      !   calcs. all |L_n>
      call dgemm('N', 'N', nXovAct, nExc, subSpaceDim, 1.0_dp, vecB, nXovAct, evecL, memDim,&
          & 0.0_dp, vecB(1,subSpaceDim+1+nExc), nXovAct)

      do ii = 1, nExc
        dummyReal = -sqrt(evalInt(ii))
        vecB(:,subSpaceDim + ii) = dummyReal * vecB(:, subSpaceDim + ii)
        vecB(:,subSpaceDim + nExc + ii) = dummyReal * vecB(:, subSpaceDim + nExc + ii)
      end do

      ! (A-B)|L_n> for all n=1,..,nExc
      call dgemm('N', 'N', nXovAct, nExc, subSpaceDim, 1.0_dp, vM, nXovAct, evecL, memDim, 1.0_dp,&
          & vecB(1, subSpaceDim + 1), nXovAct)
      ! (A+B)|R_n> for all n=1,..,nExc
      call dgemm('N', 'N', nXovAct, nExc, subSpaceDim, 1.0_dp, vP, nXovAct, evecR, memDim, 1.0_dp,&
          & vecB(1, subSpaceDim + 1 + nExc), nXovAct)

      ! calc. norms of residual vectors to check for convergence
      didConverge = .true.
      do ii = subSpaceDim + 1, subSpaceDim + nExc
        vecNorm(ii-subSpaceDim) = dot_product(vecB(:,ii), vecB(:,ii))
        if (vecNorm(ii-subSpaceDim) .gt. convThreshold) then
          didConverge = .false.
        end if
      end do

      if (didConverge) then
        do ii = subSpaceDim + nExc + 1, subSpaceDim + 2 * nExc
          vecNorm(ii-subSpaceDim) = dot_product(vecB(:,ii), vecB(:,ii))
          if (vecNorm(ii-subSpaceDim) .gt. convThreshold) then
            didConverge = .false.
          end if
        end do
      end if

      if ((.not. didConverge) .and. (subSpaceDim + 2 * nExc > nXov)) then
        print *, 'Error: Linear Response calculation in subspace did not converge!'
        stop
      end if

      ! if converged then exit loop:
      if (didConverge) then
        ! do some wrapup, e.g.:
        write (*,'(A,I4)')&
            & 'Linear Response calculation converged. Required no. of basis vectors: ', subSpaceDim
        eval(:) = evalInt(1:nExc)
        ! Calc. F_I
        ! vecXpY_test = matmul(vecB(:,1:subSpaceDim), Mh(1:subSpaceDim,1:nExc))
        evec = matmul(vecB, Mh(:,1:nExc))
        ! Calc. X+Y
        vecXpY = matmul(vecB(:,1:subSpaceDim), evecR(1:subSpaceDim,:))
        ! Calc. X-Y for exc. nStat, will be used in force calculation
        if (nStat > 0) then
          vecXmY(:) = 0.0_dp
          do ii = 1, subSpaceDim
            vecXmY(:) = vecXmY + evecL(ii,nStat) * vecB(:,ii)
          end do
        end if

        exit solveLinResp ! terminate diag. routine
      end if

      ! Otherwise calculate new basis vectors and extend subspace with them
      ! only include new vectors if they add meaningful residue component
      newVec = 0
      do ii = 1, nExc
        if (vecNorm(ii) .gt. convThreshold) then
          newVec = newVec + 1
          dummyReal = sqrt(evalInt(ii))
          info = subSpaceDim + ii
          dummyInt = subSpaceDim + newVec
          do jj = 1, nXov
            vecB(jj,dummyInt) = vecB(jj,info) / (dummyReal - wIJ(jj))
          end do
        end if
      end do

      do ii = 1, nExc
        if (vecNorm(nExc+ii) .gt. convThreshold) then
          newVec = newVec + 1
          info = subSpaceDim + nExc + ii
          dummyInt = subSpaceDim + newVec
          do jj = 1, nXov
            vecB(jj,dummyInt) = vecB(jj,info) / (dummyReal - wIJ(jj))
          end do
        end if
      end do

      prevSubSpaceDim = subSpaceDim
      subSpaceDim = subSpaceDim + newVec

      call orthonormalizeVectors(prevSubSpaceDim + 1, subSpaceDim, vecB) !create orthogonal basis

    end do solveLinResp

    if (.not. tZVector) then !not required anymore if gradient not computed
      if (.not. tTransQ) then
        deallocate(tQov)
      end if
      deallocate(tQoo)
      deallocate(tQvv)
    end if

  end subroutine rsLinRespCalc


  !> Write out transitions from ground to excited state along with single particle transitions
  !> and dipole strengths
  !> Modified routine because some expresions in original routine don't apply in the RS case
  subroutine writeExcitationsRS(cSym, oscStrength, nExc, nMatUp, getIJ, win, eval, mXpYall, wIJ,&
      & fdXplusY, fdTrans, fdTraDip, transitionDipoles, tWriteTagged, fdTagged, taggedWriter,&
      & fdExc, Ssq)

    !> Symmetry label for the type of transition
    character, intent(in) :: cSym

    !> oscillator strengths for transitions from ground to excited states
    real(dp), intent(in) :: oscStrength(:)

    !> number of excited states to solve for
    integer, intent(in) :: nExc

    !> number of same spin excitations
    integer, intent(in) :: nMatUp

    !> index array between transitions in square and 1D representations
    integer, intent(in) :: getIJ(:,:)

    !> index array for single particle excitions
    integer, intent(in) :: win(:)

    !> excitation energies
    real(dp), intent(in) :: eval(:)

    !> eigenvectors of excited states
    !real(dp), intent(in) :: evec(:,:)

    !> X+Y, all of them
    real(dp), intent(in) :: mXpYall(:,:)

    !> single particle excitation energies
    real(dp), intent(in) :: wIJ(:)

    !> single particle transition dipole moments
    real(dp), intent(in) :: transitionDipoles(:,:)

    !> should tagged information be written out
    logical, intent(in) :: tWriteTagged

    !> file unit for transition dipoles
    integer, intent(in) :: fdTraDip

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

    integer :: nMat
    integer :: i, j, iWeight, indO, m, n
    integer :: iDeg
    real(dp), allocatable :: vecW(:)
    real(dp), allocatable :: eDeg(:)
    real(dp), allocatable :: oDeg(:)
    integer, allocatable :: vecWi(:)
    real(dp) :: weight, vecWnorm
    logical :: lUpdwn, tSpin
    character :: cSign

    @:ASSERT(fdExc > 0)

    tSpin = present(Ssq)
    nMat = size(wIJ)
    allocate(vecW(nMat))
    allocate(vecWi(nMat))
    allocate(eDeg(nExc))
    allocate(oDeg(nExc))
    vecW = 0.0_dp
    vecWi = 0
    eDeg = 0.0_dp
    oDeg = 0.0_dp

    if(fdXplusY > 0) then
      write(fdXplusY,*) nMat, nExc
    end if

    do i = 1, nExc
      if (eval(i) > 0.0_dp) then

        vecW(:) = mXpYall(:,i)**2
        vecWnorm = 1.0_dp / sqrt(sum(vecW**2))
        vecW(:) = vecW(:) * vecWnorm

        ! find largest coefficient in CI - should use maxloc
        call index_heap_sort(vecWi, vecW)
        vecWi = vecWi(size(vecWi):1:-1)
        vecW = vecW(vecWi)

        weight = vecW(1)
        iWeight = vecWi(1)

        call indXov(win, iWeight, getIJ, m, n)
        cSign = cSym
        if (tSpin) then
          cSign = " "
          write (fdExc, '(1x,f10.3,4x,f14.8,2x,i5,3x,a,1x,i5,7x,f6.3,2x,f10.3,4x,f6.3)')&
              & Hartree__eV * sqrt(eval(i)), oscStrength(i), m, '->', n, weight,&
              & Hartree__eV * wIJ(iWeight), Ssq(i)
        else
          write (fdExc, '(1x,f10.3,4x,f14.8,5x,i5,3x,a,1x,i5,7x,f6.3,2x,f10.3,6x,a)')&
              & Hartree__eV * sqrt(eval(i)), oscStrength(i), m, '->', n, weight,&
              & Hartree__eV * wIJ(iWeight), cSign
        end if

        if (fdXplusY > 0) then
          if (tSpin) then
            lUpdwn = (win(iWeight) <= nMatUp)
            cSign = "D"
            if (lUpdwn) cSign = "U"
          end if
          write(fdXplusY, '(1x,i5,3x,a,3x,ES17.10)') i, cSign, sqrt(eval(i))
          write(fdXplusY, '(6(1x,ES17.10))') mXpYall(:,i)
        end if

        if (fdTrans > 0) then
          write(fdTrans, '(2x,a,T12,i5,T21,ES17.10,1x,a,2x,a)')&
              & 'Energy ', i, Hartree__eV * sqrt(eval(i)), 'eV', cSign
          write(fdTrans, *)
          write(fdTrans, '(2x,a,9x,a,8x,a)') 'Transition', 'Weight', 'KS [eV]'
          write(fdTrans, '(1x,45("="))')

          cSign = " "
          do j = 1, nMat
            indO = vecWi(j)
            call indXov(win, indO, getIJ, m, n)
            if (tSpin) then
              lUpdwn = (win(indO) <= nMatUp)
              cSign = "D"
              if (lUpdwn) cSign = "U"
            end if
            write(fdTrans, '(i5,3x,a,1x,i5,1x,1a,T22,f10.8,T33,f14.8)')&
                & m, '->', n, cSign, vecW(j), Hartree__eV * wIJ(vecWi(j))
          end do
          write(fdTrans, *)
        end if

        if(fdTraDip > 0) then
          write(fdTraDip, '(1x,i5,1x,f10.3,2x,3(ES13.6))')&
              & i, Hartree__eV * sqrt(eval(i)), (transitionDipoles(i,j) * au__Debye, j=1,3)
        end if
      else

        ! find largest coefficient in CI - should use maxloc
        call index_heap_sort(vecWi, vecW)
        vecWi = vecWi(size(vecWi):1:-1)
        vecW = vecW(vecWi)

        weight = vecW(1)
        iWeight = vecWi(1)
        call indXov(win, iWeight, getIJ, m, n)
        cSign = cSym

        if (tSpin) then
          cSign = " "
          write(fdExc, '(6x,A,T12,4x,f14.8,2x,i5,3x,a,1x,i5,7x,A,2x,f10.3,4x,f6.3)')&
              & '< 0', oscStrength(i), m, '->', n, '-', Hartree__eV * wIJ(iWeight), Ssq(i)
        else
          write(fdExc, '(6x,A,T12,4x,f14.8,2x,i5,3x,a,1x,i5,7x,f6.3,2x,f10.3,6x,a)')&
              & '< 0', oscStrength(i), m, '->', n, weight, Hartree__eV * wIJ(iWeight), cSign
        end if

        if (fdXplusY > 0) then
          if (tSpin) then
            lUpdwn = (win(iWeight) <= nMatUp)
            cSign = "D"
            if (lUpdwn) cSign = "U"
          end if
          write(fdXplusY, '(1x,i5,3x,a,3x,A)') i, cSign, '-'
        end if

        if (fdTrans > 0) then
          write (fdTrans, '(2x,a,1x,i5,5x,a,1x,a,3x,a)') 'Energy ', i,  '-', 'eV', cSign
          write (fdTrans,*)
        end if

        if(fdTraDip > 0) then
          write(fdTraDip, '(1x,i5,1x,A)') i, '-'
        end if

      end if
    end do

    ! Determine degenerate levels and sum oscillator strength over any degenerate levels
    iDeg = 1
    eDeg(1) = eval(1)
    oDeg(1) = oscStrength(1)
    do i = 2, nExc
      if(abs(eval(i) - eval(i-1)) < elecTolMax) then
        oDeg(iDeg) = oDeg(iDeg) + oscStrength(i)
      else
        iDeg = iDeg + 1
        eDeg(iDeg) = eval(i)
        oDeg(iDeg) = oscStrength(i)
      end if
    end do

    if (tWriteTagged) then
      call taggedWriter%write(fdTagged, tagLabels%excEgy, eDeg(:iDeg))
      call taggedWriter%write(fdTagged, tagLabels%excOsc, oDeg(:iDeg))
    end if

  end subroutine writeExcitationsRS


  !> Write atomic transition charges to file
  subroutine writeTransitionChargesRS(fdTransQ, atomicTransQ)
    !> file unit for transition charges
    integer, intent(in)  :: fdTransQ
    !> transition charges to write
    real(dp), intent(in) :: atomicTransQ(:)

    integer :: nAtom, i
    nAtom = size(atomicTransQ)

    open(fdTransQ, file=transChargesOut, action="write", status="replace")
    write(fdTransQ, '(a)') "#"
    write(fdTransQ, '(a,2x,a,5x,a)') "#", "atom", "transition charge"
    write(fdTransQ, '(a)') "#"

    do i = 1, nAtom
      write(fdTransQ, '(2x,i5,5x,f12.9)') i, atomicTransQ(i)
    end do

    close(fdTransQ)

  end subroutine writeTransitionChargesRS


  !> Computes excitation spectrum through range separated response calculation
  subroutine linRespCalcExcitationsRS(spin, tOnsite, this, iAtomStart, eigVec, eigVal, sccCalc,&
      & SSqrReal, filling, coords0, dqAt, specie0, hubbUAtom, iNeighbor, img2CentCell, orb, rsData,&
      & tWriteTagged, fdTagged, taggedWriter, excEnergy, skHamCont, skOverCont, derivator,&
      & deltaRho, excGrad, dQAtomEx)
    logical, intent(in) :: spin
    logical, intent(in) :: tOnsite
    type(TLinResp), intent(inout) :: this
    integer, intent(in) :: iAtomStart(:)
    real(dp), intent(in) :: eigVec(:,:,:)
    real(dp), intent(in) :: eigVal(:,:), SSqrReal(:,:)
    !> Self-consistent charge module settings
    type(TScc), intent(in) :: sccCalc
    real(dp), intent(in) :: filling(:,:)
    real(dp), intent(in) :: coords0(:,:)
    real(dp), intent(in) :: dqAt(:)
    integer, intent(in) :: specie0(:)
    real(dp), intent(in) :: hubbUAtom(:)
    integer, intent(in) :: iNeighbor(0:,:), img2CentCell(:)
    type(TOrbitals), intent(in) :: orb
    type(TRangeSepFunc), intent(inout) :: rsData
    !> print tag information
    logical, intent(in) :: tWriteTagged

    !> file id for tagging information
    integer, intent(in) :: fdTagged

    !> tagged writer
    type(TTaggedWriter), intent(inout) :: taggedWriter

    !real(dp), intent(in) :: ons_en(:,:), ons_dip(:,:)
    real(dp), intent(out) :: excEnergy
    !real(dp), intent(in), optional :: shift(:)
    type(TSlakoCont), intent(in), optional :: skHamCont, skOverCont
    !> Differentiatior for the non-scc matrices
    class(TNonSccDiff), intent(in), optional :: derivator
    real(dp), intent(inout), optional :: deltaRho(:,:)
    real(dp), intent(inout), optional :: excGrad(:,:)
    real(dp), intent(inout), optional :: dQAtomEx(:)

    real(dp), allocatable :: shiftPerAtom(:), shiftPerL(:,:)
    integer :: nAtom !, i
    real(dp), allocatable :: occNr(:,:)

    if (any(abs(mod(filling, 2.0_dp)) > epsilon(0.0_dp))) then
      call error("Fractionally occupied states not currently supported for range separated linear&
          & response excitations")
    end if

    nAtom = size(orb%nOrbAtom)

    allocate(occNr(size(filling,dim=1), size(filling,dim=2)))
    occNr = filling(:,:)

    if (.not. present(excGrad)) then
      call runRsLinRespCalc(spin, tOnsite, nAtom, iAtomStart, eigVec, eigVal, sccCalc, dqAt,&
          & coords0, this%nExc, this%nStat, this%symmetry, SSqrReal, occNr, specie0, this%nAtom,&
          & hubbUAtom, this%spinW, this%nEl, iNeighbor, img2CentCell, orb, rsData, tWriteTagged,&
          & fdTagged, taggedWriter, this%fdMulliken, this%fdCoeffs, this%fdXplusY, this%fdTrans,&
          & this%fdSPTrans, this%fdTraDip, this%fdTransQ, this%tArnoldi, this%fdArnoldi,&
          & this%fdExc, this%tEnergyWindow, this%energyWindow, this%tOscillatorWindow,&
          & this%oscillatorWindow, this%tCacheCharges, excEnergy)
    else
      allocate(shiftPerAtom(nAtom))
      allocate(shiftPerL(orb%mShell, nAtom))
      call sccCalc%getShiftPerAtom(shiftPerAtom)
      call sccCalc%getShiftPerL(shiftPerL)
      shiftPerAtom = shiftPerAtom + shiftPerL(1,:)

      call runRsLinRespCalc(spin, tOnsite, nAtom, iAtomStart, eigVec, eigVal, sccCalc, dqAt,&
          & coords0, this%nExc, this%nStat, this%symmetry, SSqrReal, occNr, specie0, this%nAtom,&
          & hubbUAtom, this%spinW, this%nEl, iNeighbor, img2CentCell, orb, rsData, tWriteTagged,&
          & fdTagged, taggedWriter, this%fdMulliken, this%fdCoeffs, this%fdXplusY, this%fdTrans,&
          & this%fdSPTrans, this%fdTraDip, this%fdTransQ, this%tArnoldi, this%fdArnoldi,&
          & this%fdExc, this%tEnergyWindow, this%energyWindow, this%tOscillatorWindow,&
          & this%oscillatorWindow, this%tCacheCharges, excEnergy, shiftPerAtom, skHamCont,&
          & skOverCont, derivator, deltaRho, excGrad, dQAtomEx)
    end if

  end subroutine linRespCalcExcitationsRS


  ! Functions for gradient calculations:

  !> Create P = T + 1/2 Z symmetric (paper has T + Z asymmetric)
  !> (Zab = Zij = 0, Tia = 0)
  subroutine calcPMatrix(t, rhs, win, getIJ, pc)
    real(dp), intent(in) :: t(:,:), rhs(:)
    integer, intent(in) :: win(:), getIJ(:,:)
    real(dp), intent(out) :: pc(:,:)

    integer :: ia, i, a

    pc(:,:) = t(:,:)
    do ia = 1, size(rhs)
      call indXov(win, ia, getIJ, i, a)
      ! shouldn't be pc = pc + 1/2 rhs?
      pc(i,a) = 0.5_dp * rhs(ia)
      pc(a,i) = 0.5_dp * rhs(ia)
    end do

  end subroutine calcPMatrix

  !> Transform dense matrix from MO-representation into AO representation
  subroutine transformMO2AODense(aa, cc)
    real(dp), intent(inout) :: aa(:,:)
    real(dp), intent(in) :: cc(:,:)

    real(dp), allocatable :: buffer(:,:)

    allocate(buffer(size(aa, dim=1), size(aa, dim=2)))

    !  Transform A from MO representation into AO repr:
    buffer(:,:) = 0.0_dp
    call gemm(buffer, aa, cc, transA="N", transB="T")
    call gemm(aa, cc, buffer, transA="N", transB="N")

  end subroutine transformMO2AODense

  subroutine getExcMulliken(iAtomStart, pc, s, dQAtomEx)
    integer, intent(in) :: iAtomStart(:)
    real(dp), intent(in) :: pc(:,:), s(:,:)
    real(dp), intent(out) :: dQAtomEx(:)

    integer :: alpha, beta, ia1, ia2, ib1, ib2
    real(dp) :: tmp

    dQAtomEx(:) = 0.0_dp
    do alpha = 1, size(dQAtomEx)
      ia1 = iAtomStart(alpha)
      ia2 = iAtomStart(alpha + 1) - 1
      do beta = 1, alpha
        ib1 = iAtomStart(beta)
        ib2 = iAtomStart(beta + 1) - 1
        tmp = sum(pc(ia1:ia2, ib1:ib2) * s(ia1:ia2, ib1:ib2))
        dQAtomEx(alpha) = dQAtomEx(alpha) + tmp
        if (alpha /= beta) then
          dQAtomEx(beta) = dQAtomEx(beta) + tmp
        end if
      end do
    end do

  end subroutine getExcMulliken


  !> Build right hand side of the equation for the Z-vector and those parts of the W-vectors which
  !> do not depend on Z.
  subroutine getZVectorEqRhsRS(vecXpY, vecXmY, win, iAtomStart, homo, nOcc, nMatUp, getIJ,&
      & iaTrans, nAtom, species0, ev, sTimesGrndEigVecs, grndEigVecs, gamma, lrGamma,&
      & spinW, omega, cSym, rhs, t, vWov, vWoo, vWvv)
    real(dp), intent(in) :: vecXpY(:), vecXmY(:)
    integer, intent(in) :: win(:), iAtomStart(:)
    integer, intent(in) :: homo, nOcc, nMatUp, getIJ(:,:), iaTrans(1:,homo+1:), nAtom, species0(:)
    real(dp), intent(in) :: ev(:), sTimesGrndEigVecs(:,:,:), grndEigVecs(:,:,:)
    real(dp), intent(in) :: gamma(:,:), lrGamma(:,:), spinW(:), omega
    character, intent(in) :: cSym
    real(dp), intent(out) :: rhs(:), t(:,:)
    real(dp), intent(out) :: vWov(:), vWoo(:), vWvv(:)

    real(dp), allocatable :: vecXpYQ(:), qIJ(:), gamXpYQ(:), qGamXPyQ(:), gamQt(:)
    real(dp), allocatable :: vecHvvX(:), vecHvvY(:), vecHooX(:), vecHooY(:), vecHovT(:), vecHooT(:)
    integer :: nVirt, nOrb, nXov, nXoo, nXvv
    integer :: i, j, a, b, ia, ib, ij, ab, ja
    real(dp) :: tmp1, tmp2
    logical :: lUpDwn

    nXov = size(rhs)
    nXoo = size(vWoo)
    nXvv = size(vWvv)
    ! the below does not work with energyWindow (nXovRD applies and not nXov, and nXovRD < nXov)
    !nVirt = nXov / nOcc
    !nOrb = nOcc + nVirt
    nOrb = iAtomStart(nAtom+1) - 1
    nVirt = nOrb - nOcc

    @:ASSERT(nXvv == nVirt * (nVirt + 1) / 2)
    @:ASSERT(nXoo == nOcc * (nOcc + 1) / 2)

    allocate(vecXpYQ(nAtom))
    allocate(qIJ(nAtom))
    allocate(gamXpYQ(nAtom))
    allocate(gamQt(nAtom))
    allocate(qGamXpYQ(max(nXoo, nXvv) * nOrb)) !check how much memory is really required!

    t(:,:) = 0.0_dp
    rhs(:) = 0.0_dp
    vWov(:) = 0.0_dp
    vWoo(:) = 0.0_dp
    vWvv(:) = 0.0_dp
    vecXpYQ(:) = 0.0_dp

    ! optim code
    ! Build t_ab = 0.5 * sum_i (X+Y)_ia (X+Y)_ib + (X-Y)_ia (X-Y)_ib
    ! and w_ab = Q_ab with Q_ab as in (B16) but with corrected sign.
    ! factor 1 / (1 + delta_ab) follows later
    do ia = 1, nXov
      call indXov(win, ia, getIJ, i, a)
      lUpDwn = (win(ia) <= nMatUp)
      ! BA: Is T_aa = 0?
      do b = homo + 1, a
        ib = iaTrans(i, b)
        ab = (b - homo) + (((a - homo) - 1) * (a - homo)) / 2
        tmp1 = vecXpY(ia) * vecXpY(ib) + vecXmY(ia) * vecXmY(ib)
        tmp2 = omega * (vecXpY(ia) * vecXmY(ib)+ vecXmY(ia) * vecXpY(ib))
        t(a,b) = t(a,b) + 0.5_dp * tmp1
        ! to prevent double counting
        if (a /= b) then
          t(b,a) = t(b,a) + 0.5_dp * tmp1
        end if
        ! Note: diagonal elements will be multiplied by 0.5 later.
        vWvv(ab) = vWvv(ab) + ev(i) * tmp1 + tmp2
      end do

      ! Build t_ij = 0.5 * sum_a (X+Y)_ia (X+Y)_ja + (X-Y)_ia (X-Y)_ja
      ! and 1 / (1 + delta_ij) Q_ij with Q_ij as in eq. (B9) (1st part of w_ij)
      do j = i, homo
        ja = iaTrans(j,a)
        ij = i - homo + nOcc + ((j - homo + nOcc - 1) * (j - homo + nOcc)) / 2
        tmp1 = vecXpY(ia) * vecXpY(ja) + vecXmY(ia) * vecXmY(ja)
        tmp2 = omega * (vecXpY(ia) * vecXmY(ja) + vecXmY(ia) * vecXpY(ja))
        ! Note, there is a typo in Heringer et al. J. Comp Chem 28, 2589.
        ! The sign must be negative see Furche, J. Chem. Phys, 117 7433 (2002).
        t(i,j) = t(i,j) - 0.5_dp * tmp1
        ! to prevent double counting
        if (i /= j) then
          t(j,i) = t(j,i) - 0.5_dp * tmp1
        end if
        vWoo(ij) = vWoo(ij) - ev(a) * tmp1 + tmp2
      end do
    end do

    ! Build xpyq = sum_ia (X+Y)_ia
    do ia = 1, nXov
      call indXov(win, ia, getIJ, i, a)
      lUpDwn = (win(ia) <= nMatUp)
      qIJ = transQ(i, a, iAtomStart, lUpDwn, sTimesGrndEigVecs, grndEigVecs)
      vecXpYQ(:) = vecXpYQ + vecXpY(ia) * qIJ
    end do

    gamXpYQ(:) = matmul(gamma, vecXpYQ)

    ! qgamxpyq(ab) = sum_jc K_ab,jc (X+Y)_jc
    do ab = 1, nXvv
      call indXvv(homo, ab, a, b)
      lUpDwn = .true. ! UNTESTED
      qIJ = transQ(a, b, iAtomStart, lUpDwn, sTimesGrndEigVecs, grndEigVecs)
      if (cSym == "S") then
        qGamXpYQ(ab) = 2.0_dp * sum(qIJ * gamXpYQ)
      else
        qGamXpYQ(ab) = sum(qIJ * vecXpYq * (spinW(species0)))
      end if
    end do

    ! rhs(ia) -= Qia = sum_b (X+Y)_ib * qgamxpyq(ab))
    do ia = 1, nXov
      call indXov(win, ia, getIJ, i, a)
      lUpDwn = (win(ia) <= nMatUp)
      do b = homo + 1, a
        if (b < a) then
          call rIndXvv(homo, a, b, ab)
        else
          call rIndXvv(homo, b, a, ab)
        end if
        ib = iaTrans(i,b)
        rhs(ia) = rhs(ia) - 2.0_dp * vecXpY(ib) * qGamXpYQ(ab)
        ! Since qgamxpyq has only upper triangle
        if (a /= b) then
          rhs(ib) = rhs(ib) - 2.0_dp * vecXpY(ia) * qGamXpYQ(ab)
        end if
      end do
    end do

    ! -rhs = -rhs - sum_j (X + Y)_ja H + _ij[X + Y]
    ! qgamxpyq(ij) = sum_kb K_ij,kb (X+Y)_kb
    do ij = 1, nXoo
      qGamXpYQ(ij) = 0.0_dp
      call indXoo(ij, i, j)
      lUpDwn = .true.
      qIJ = transQ(i, j, iAtomStart, lUpDwn, sTimesGrndEigVecs, grndEigVecs)
      if (cSym == "S") then
        qGamXpYQ(ij) = 2.0_dp * sum(qIJ * gamXpYQ)
      else
        qGamXpYQ(ij) = sum(qIJ * vecXpYQ * (spinW(species0)))
      end if
    end do

    ! rhs(ia) += Qai = sum_j (X+Y)_ja qgamxpyq(ij)
    ! add Qai to Wia as well.
    do ia = 1, nXov
      call indXov(win, ia, getIJ, i, a)
      do j = i, homo
        ja = iaTrans(j, a)
        ij = i - homo + nOcc + ((j - homo + nOcc - 1) * (j - homo + nOcc)) / 2
        tmp1 = 2.0_dp * vecXpY(ja) * qGamXpYQ(ij)
        rhs(ia) = rhs(ia) + tmp1
        vWov(ia) = vWov(ia) + tmp1
        if (i /= j) then
          tmp2 = 2.0_dp * vecXpY(ia) * qGamXpYQ(ij)
          rhs(ja) = rhs(ja) + tmp2
          vWov(ja) = vWov(ja) + tmp2
        end if
      end do
    end do

    ! gamxpyq(beta) = sum_ij q_ij(beta) T_ij
    gamXpYQ(:) = 0.0_dp
    do ij = 1, nXoo
      call indXoo(ij, i, j)
      lUpDwn = .true.
      qIJ = transQ(i, j, iAtomStart, lUpDwn, sTimesGrndEigVecs, grndEigVecs)
      if (i == j) then
        gamXpYQ(:) = gamXpYQ(:) + t(i,j) * qIJ(:)
      else
        ! factor 2 because of symmetry of the matrix
        gamXpYQ(:) = gamXpYQ(:) + 2.0_dp  * t(i,j) * qIJ(:)
      end if
    end do

    ! gamxpyq(beta) += sum_ab q_ab(beta) T_ab
    do ab = 1, nXvv
      call indXvv(homo, ab, a, b)
      lUpDwn = .true.
      qIJ = transQ(a, b, iAtomStart, lUpDwn, sTimesGrndEigVecs, grndEigVecs)
      if (a == b) then
        gamXpYQ(:) = gamXpYQ(:) + t(a,b) * qIJ(:)
      else
        ! factor 2 because of symmetry of the matrix
        gamXpYQ(:) = gamXpYQ(:) + 2.0_dp * t(a,b) * qIJ(:)
      end if
    end do

    ! gamqt(alpha) = sum_beta gamma_alpha,beta gamxpyq(beta)
    gamQt(:) = matmul(gamma, gamXpYQ)

    ! rhs -= sum_q^ia(alpha) gamxpyq(alpha)
    do ia = 1, nXov
      call indXov(win, ia, getIJ, i, a)
      lUpDwn = (win(ia) <= nMatUp)
      qIJ = transQ(i, a, iAtomStart, lUpDwn, sTimesGrndEigVecs, grndEigVecs)
      rhs(ia) = rhs(ia) - 4.0_dp * sum(qIJ * gamQt)
    end do

    ! Furche vectors
    do ij = 1, nXoo
      call indXoo(ij, i, j)
      lUpDwn = .true.
      qIJ = transQ(i, j, iAtomStart, lUpDwn, sTimesGrndEigVecs, grndEigVecs)
      vWoo(ij) = vWoo(ij) + 4.0_dp * sum(qIJ * gamQt)
    end do

    ! Contributions due to range-separation
    allocate(vecHvvX(nXvv))
    allocate(vecHvvY(nXvv))
    allocate(vecHooX(nXoo))
    allocate(vecHooY(nXoo))
    allocate(vecHovT(nXov))

    call getHvvXY(1, nXov, nXvv, nOrb, homo, nAtom, nMatUp, iaTrans, getIJ, win, iAtomStart,&
        & sTimesGrndEigVecs, grndEigVecs, lrGamma, vecXpY, vecHvvX)
    call getHvvXY(-1, nXov, nXvv, nOrb, homo, nAtom, nMatUp, iaTrans, getIJ, win, iAtomStart,&
        & sTimesGrndEigVecs, grndEigVecs, lrGamma, vecXmY, vecHvvY)
    call getHooXY(1, nXov, nXoo, nOrb, homo, nAtom, nMatUp, iaTrans, getIJ, win, iAtomStart,&
        & sTimesGrndEigVecs, grndEigVecs, lrGamma, vecXpY, vecHooX)
    call getHooXY(-1, nXov, nXoo, nOrb, homo, nAtom, nMatUp, iaTrans, getIJ, win, iAtomStart,&
        & sTimesGrndEigVecs, grndEigVecs, lrGamma, vecXmY, vecHooY)
    call getHovT(nOcc, nXov, nXoo, nXvv, nOrb, homo, nAtom, nMatUp, iaTrans, getIJ, win,&
        & iAtomStart, sTimesGrndEigVecs, grndEigVecs, lrGamma, t, vecHovT)

    do ia = 1, nXov
      call indXov(win, ia, getIJ, i, a)
      lUpDwn = (win(ia) <= nMatUp)
      do b = homo + 1, nOrb
        ib = iaTrans(i, b)
        if (ib <= nXov) then
          if (b < a) then
            call rIndXvv(homo, a, b, ab)
          else
            call rIndXvv(homo, b, a, ab)
          end if
          rhs(ia) = rhs(ia) - vecXpY(ib) * vecHvvX(ab)
          if (a >= b) then
            rhs(ia) = rhs(ia) - vecXmY(ib) * vecHvvY(ab)
          else
            rhs(ia) = rhs(ia) + vecXmY(ib) * vecHvvY(ab)
          end if
        end if
      end do

      do j = 1, homo
        ja = iaTrans(j, a)
        if (ja <= nXov) then
          if (j < i) then
            call rIndXoo(homo, nOcc, i, j, ij)
          else
            call rIndXoo(homo, nOcc, j, i, ij)
          end if
          rhs(ia) = rhs(ia) + vecXpY(ja)*vecHooX(ij)
          vWov(ia) = vWov(ia) + vecXpY(ja)*vecHooX(ij)
          if (i >= j) then
            rhs(ia) = rhs(ia) + vecXmY(ja)*vecHooY(ij)
            vWov(ia) = vWov(ia) + vecXmY(ja)*vecHooY(ij)
          else
            rhs(ia) = rhs(ia) - vecXmY(ja)*vecHooY(ij)
            vWov(ia) = vWov(ia) - vecXmY(ja)*vecHooY(ij)
          end if
        end if
      end do
      rhs(ia) = rhs(ia) - vecHovT(ia)
    end do

    deallocate(vecHvvX)
    deallocate(vecHvvY)
    deallocate(vecHooX)
    deallocate(vecHooY)
    deallocate(vecHovT)

    allocate(vecHooT(nXoo))
    call getHooT(nOcc, nXov, nXoo, nOrb, homo, nAtom, nMatUp, iaTrans, getIJ, win, iAtomStart,&
        & sTimesGrndEigVecs, grndEigVecs, lrGamma, t, vecHooT)
    vWoo(:) = vWoo + vecHooT
    deallocate(vecHooT)

    deallocate(vecXpYq)
    deallocate(qIJ)
    deallocate(gamXpYQ)
    deallocate(gamQt)
    deallocate(qGamXpYQ)

  end subroutine getZVectorEqRhsRS


  !> Solve the (A+B) Z = -R equation via conjugate gradient optimization.
  subroutine solveZVectorEqRS(rhs, win, nMatUp, getIJ, occNr, nAtom, species0, spinW, gamma,&
      & lrGamma, wIJ, homo, nOcc, nVir, iaTrans, tQov, tQoo, tQvv)
    real(dp), intent(inout) :: rhs(:)
    integer, intent(in) :: win(:), nMatUp, getIJ(:,:), nAtom
    real(dp), intent(in) :: occNr(:,:)
    real(dp), intent(in) :: gamma(:,:), lrGamma(:,:), wIJ(:)
    integer, intent(in) :: species0(:)
    real(dp), intent(in) :: spinW(:)
    integer, intent(in) :: homo, nOcc, nVir
    real(dp), intent(in) :: tQov(:,:), tQoo(:,:), tQvv(:,:)
    integer, intent(in) :: iaTrans(1:,homo+1:)

    integer :: nXov
    integer :: ia, i, a, k, ii, aa
    real(dp), allocatable :: qIJ(:), gXqIJ(:), rhs2(:), rkm1(:), pkm1(:), apk(:)
    real(dp) :: rs, tmp1, tmp2
    logical :: lUpDwn
    real(dp), allocatable :: gqvTmp(:,:)

    !dummy variables to cal. (A+B)*vec
    character :: cSym
    logical :: spin
    cSym = 'S'
    spin = .false.

    if (dot_product(rhs, rhs) <= 1.0E-6_dp) then ! May actually be zero, for example for H_2
      rhs(:) = 0.0_dp
      return
    end if

    nXov = size(rhs)
    allocate(qIJ(nAtom))
    allocate(gXqIJ(nAtom))
    allocate(rhs2(nXov))
    allocate(rkm1(nXov))
    allocate(pkm1(nXov))
    allocate(apk(nXov))
    allocate(gqvTmp(size(gamma, dim=1), max(nVir, nOcc) * nVir))

    ! Choosing a start value
    ! rhs2 = rhs / (A+B)_ia,ia (diagonal of the supermatrix sum A+B)
    do ia = 1, nXov
      call indXov(win, ia, getIJ, i, a)
      lUpDwn = (win(ia) <= nMatUp)
      ! qIJ = transQ(i, a, iAtomStart, lUpDwn, sTimesGrndEigVecs, grndEigVecs)
      rs = 4.0_dp * dot_product(tQov(:,ia), matmul(gamma, tQov(:,ia))) + wIJ(ia)
      rs = rs - dot_product(tQov(:,ia), matmul(lrGamma, tQov(:,ia))) ! part of rs contirb
      ! !Bug here! ii = i - homo + nOcc + ((i - homo + nOcc + 1) * (i - homo + nOcc))/2
      ii = i - homo + nOcc + ((i - homo + nOcc - 1) * (i - homo + nOcc))/2
      ! lUpDwn = .true.
      ! qIJ = transQ(i, i, iAtomStart, lUpDwn, sTimesGrndEigVecs, grndEigVecs)
      gXqIJ(:) = matmul(lrGamma, tQoo(:,ii))
      call rIndXvv(homo, a, a, aa)
      ! lUpDwn = .true.
      ! qIJ = transQ(a, a, iAtomStart, lUpDwn, sTimesGrndEigVecs, grndEigVecs)
      rs = rs - dot_product(tQvv(:,aa), gXqIJ)
      rhs2(ia) = rhs(ia) / rs
    end do

    call multApBVecFast(rhs2, wIJ, cSym, win, nMatUp, homo, nOcc, nVir, occNr, getIJ, gamma,&
        & lrGamma, species0, spinW, iaTrans, gqvTmp, tQov, tQoo, tQvv, rkm1)

    rkm1(:) = rhs(:) - rkm1(:)
    pkm1(:) = rkm1(:)

    ! Iteration:  should be convergent at most nXov steps
    do k = 1, nXov

      call multApBVecFast(pkm1, wIJ, cSym, win, nMatUp, homo, nOcc, nVir, occNr, getIJ, gamma,&
          & lrGamma, species0, spinW, iaTrans, gqvTmp, tQov, tQoo, tQvv, apk)
      tmp1 = dot_product(rkm1, rkm1)
      tmp2 = dot_product(pkm1, apk)

      rhs2(:) = rhs2(:) + (tmp1 / tmp2) * pkm1(:)
      rkm1(:) = rkm1(:) - (tmp1 / tmp2) * apk(:)

      ! Residual
      tmp2 = dot_product(rkm1, rkm1)

      if (tmp2 <= 1.0e-14_dp) then
        exit
      end if

      if (k == nXov) then
        call error("LrespoGrad : Z vector not converged!")
      end if

      pkm1(:) = (tmp2 / tmp1) * pkm1(:) + rkm1(:)

    end do

    rhs(:) = rhs2(:)

    deallocate(qIJ)
    deallocate(gXqIJ)
    deallocate(rhs2)
    deallocate(rkm1)
    deallocate(pkm1)
    deallocate(apk)
    deallocate(gqvTmp)
  end subroutine solveZVectorEqRS

  !> Calculate Z-dependent parts of the W-vectors including rs contributions and divide diagonal
  !> elements of W_ij and W_ab by 2.
  subroutine calcWVectorZRS(zz, win, homo, nOrb, nMatUp, getIJ, iAtomStart, sTimesGrndEigVecs,&
      & grndEigVecs, gamma, lrGamma, ev, vWov, vWoo, vWvv, iaTrans)
    real(dp), intent(in) :: zz(:)
    integer, intent(in) :: win(:), homo, nOrb, nMatUp, getIJ(:,:), iAtomStart(:)
    real(dp), intent(in) :: sTimesGrndEigVecs(:,:,:), grndEigVecs(:,:,:), gamma(:,:), lrGamma(:,:)
    real(dp), intent(in) :: ev(:)
    real(dp), intent(inout) :: vWov(:), vWoo(:), vWvv(:)
    integer, intent(in) :: iaTrans(1:,homo+1:)

    integer :: nXov, nXoo, nXvv, nAtom
    integer :: ij, ia, ab, i, j, a, b, alpha
    real(dp), allocatable :: qIJ(:), gamXpYQ(:), zQ(:), vecHooZ(:)
    logical :: lUpDwn

    nXov = size(zz)
    nAtom = size(gamma, dim=1)
    nXoo = size(vWoo)
    nXvv = size(vWvv)
    allocate(qIJ(nAtom))
    allocate(gamXpYQ(nAtom))
    allocate(zQ(nAtom))
    ! for rangesep contributions:
    allocate(vecHooZ(nXoo))

    ! Adding missing epsilon_i * Z_ia term to W_ia
    do ia = 1, nXov
      call indXov(win, ia, getIJ, i, a)
      lUpDwn = (win(ia) <= nMatUp)
      vWov(ia) = vWov(ia) + zz(ia) * ev(i)
    end do

    ! Missing sum_kb 4 K_ijkb Z_kb term in W_ij:
    ! zq(alpha) = sum_kb q^kb(alpha) Z_kb
    do alpha = 1, nAtom
      zQ(alpha) = 0.0_dp
      do ia = 1, nXov
        call indXov(win, ia, getIJ, i, a)
        lUpDwn = (win(ia) <= nMatUp)
        qIJ = transQ(i, a, iAtomStart, lUpDwn, sTimesGrndEigVecs, grndEigVecs)
        zQ(alpha) = zQ(alpha) + zz(ia) * qIJ(alpha)
      end do
    end do

    ! gamxpyq(alpha) = sum_beta gamma(alpha, beta) zq(beta)
    gamXpYQ = matmul(gamma, zQ)

    ! sum_alpha qIJ(alpha) gamxpyq(alpha)
    do ij = 1, nXoo
      call indXoo(ij, i, j)
      lUpDwn = .true.
      qIJ = transQ(i, j, iAtomStart, lUpDwn, sTimesGrndEigVecs, grndEigVecs)
      do alpha = 1, nAtom
        ! W contains 1/2 for i == j.
        vWoo(ij) = vWoo(ij) + 4.0_dp * qIJ(alpha) * gamXpYQ(alpha)
      end do
    end do

    ! Contributions due to range-separation
    call getHooXY(1, nXov, nXoo, nOrb, homo, nAtom, nMatUp, iaTrans, getIJ, win, iAtomStart,&
        & sTimesGrndEigVecs, grndEigVecs, lrGamma, zz, vecHooZ)
    vWoo(:) = vWoo(:) + vecHooZ(:)

    ! Divide diagonal elements of W_ij by 2.
    do ij = 1, nXoo
      call indXoo(ij, i, j)
      if (i == j) then
        vWoo(ij) = 0.5_dp * vWoo(ij)
      end if
    end do

    ! Divide diagonal elements of W_ab by 2.
    do ab = 1, nXvv
      call indXvv(homo, ab, a, b)
      if (a == b) then
        vWvv(ab) = 0.5_dp * vWvv(ab)
      end if
    end do

    deallocate(vecHooZ)
    deallocate(qIJ)
    deallocate(gamXpYQ)
    deallocate(zQ)
  end subroutine calcWVectorZRS


  subroutine addGradientsRS(cSym, nXov, nAtom, species0, iAtomStart, nOrb, homo, nMatUp, getIJ,&
      & win, grndEigVecs, pc, dQ, dQAtomEx, gamma, lrGamma, rsData, HubbardU, spinW, shift, vWoo,&
      & vWov, vWvv, vecXpY, vecXmY, nBeweg, coord0, orb, skHamCont, skOverCont, tQov, derivator,&
      & deltaRho, excGrad)
    character, intent(in) :: cSym
    integer, intent(in) :: nXov, nAtom, species0(:), iAtomStart(:), nOrb, homo, win(:)
    integer, intent(in) :: nMatUp, getIJ(:,:)
    real(dp), intent(in) :: grndEigVecs(:,:,:)
    real(dp), intent(in) :: pc(:,:), dQ(:), dQAtomEx(:), gamma(:,:), lrGamma(:,:)
    real(dp), intent(in) :: HubbardU(:), spinW(:), shift(:), vWoo(:), vWov(:), vWvv(:)
    real(dp), intent(in) :: vecXpY(:), vecXmY(:), coord0(:,:), tQov(:,:)
    type(TRangeSepFunc), intent(inout) :: rsData
    integer, intent(in) :: nBeweg
    type(TOrbitals), intent(in) :: orb
    type(TSlakoCont), intent(in) :: skHamCont, skOverCont
    !> Differentiatior for the non-scc matrices
    class(TNonSccDiff), intent(in) :: derivator
    real(dp), intent(inout) :: excGrad(:,:), deltaRho(:,:)

    real(dp), allocatable :: shEx(:), vecXpYq(:), shXpYQ(:), vecXpYcc(:,:), vecXmYcc(:,:), wcc(:,:)
    real(dp), allocatable :: temp(:), qIJ(:), vecXpYas(:,:), vecXmYas(:,:)
    real(dp), allocatable :: dmn(:,:)
    integer :: ia, i, j, a, b, ab, ij, m, n, mu, nu, xyz, alpha, beta
    integer :: indAlpha, indAlpha1, indBeta, indBeta1, izpAlpha, izpBeta
    real(dp) :: tmp1, tmp2, tmp3, tmp4, tmp5, tmp6, tmp7, tmprs, tmprs2, rab
    real(dp) :: diffvec(3), dgab(3), tmp3a(3), tmp3b(3)
    real(dp) :: dHMNdR, dSMNdR
    real(dp), parameter :: deltaX = epsilon(1.0_dp)**0.25_dp
    real(dp), parameter :: rcDX = 1.0_dp / deltaX
    real(dp), allocatable :: overlap(:,:), lrGammaOrb(:,:)
    real(dp), allocatable :: PS(:,:), DS(:,:), SPS(:,:), SDS(:,:), gammaLongRangePrime(:,:,:)
    real(dp), allocatable :: SX(:,:), XS(:,:), SXS(:,:), SY(:,:), YS(:,:), SYS(:,:)
    real(dp) tmp3vec(3)
    integer iAt1, iAt2, ka

    integer, allocatable :: species(:)

    real(dp), allocatable :: derH0(:,:,:), derS(:,:,:)

    logical :: lUpDwn

    integer :: nXoo, nXvv

    allocate(shEx(nAtom))
    allocate(vecXpYq(nAtom))
    allocate(shXpYQ(nAtom))
    allocate(vecXpYcc(nOrb, nOrb))
    allocate(vecXmYcc(nOrb, nOrb))
    allocate(vecXpYas(nOrb, nOrb))
    allocate(vecXmYas(nOrb, nOrb))
    allocate(wcc(nOrb, nOrb))
    allocate(qIJ(nAtom))
    allocate(temp(nOrb))
    allocate(dmn(nOrb, nOrb))

    !!  RS matrices
    allocate(PS(nOrb, nOrb))
    allocate(DS(nOrb, nOrb))
    allocate(SPS(nOrb, nOrb))
    allocate(SDS(nOrb, nOrb))
    allocate(SX(nOrb, nOrb))
    allocate(XS(nOrb, nOrb))
    allocate(SXS(nOrb, nOrb))
    allocate(SY(nOrb, nOrb))
    allocate(YS(nOrb, nOrb))
    allocate(SYS(nOrb, nOrb))
    allocate(overlap(nOrb, nOrb))
    allocate(lrGammaOrb(nOrb, nOrb))
    allocate(gammaLongRangePrime(3, nAtom, nAtom))

    allocate(derH0(orb%mOrb, orb%mOrb, 3))
    allocate(derS(orb%mOrb, orb%mOrb, 3))

    nXoo = size(vWoo)
    nXvv = size(vWvv)

    ! Density matrix by rank k-update
    ! BA: density matrix should be provided from outside
    dmn = 0._dp
    call herk(dmn, grndEigVecs(:,1:homo,1), alpha=2.0_dp)

    ! Symmetrize deltaRho
    do mu = 1, nOrb
      do nu = mu + 1, nOrb
        deltaRho(mu, nu) = deltaRho(nu, mu)
      end do
    end do

    ! Compute long-range gamma derivative
    gammaLongRangePrime(:,:,:) = 0._dp
    call rsData%getSpecies(species)
    do iAt1 = 1, nAtom
      do iAt2 = 1, nAtom
        if(iAt1 /= iAt2) then
          call getGammaPrimeValue(rsData, tmp3vec, iAt1, iAt2, coord0, species)
          gammaLongRangePrime(:, iAt1, iAt2) = tmp3vec(:)
        end if
      end do
    end do

    ! Symmetrize S
    call getSqrS(coord0, nAtom, skOverCont, orb, iAtomStart, species0, overlap)
    call getSqrGamma(nAtom, lrGamma, iAtomStart, lrGammaOrb)

    shEx(:) = matmul(dQAtomEx, gamma)

    ! xypq(alpha) = sum_ia (X+Y)_ia q^ia(alpha)
    ! complexity nOrb * nOrb * nOrb
    vecXpYq(:) = 0.0_dp
    do ia = 1, nXov
      vecXpYq(:) = vecXpYq(:) + vecXpY(ia) * tQov(:,ia)
    end do

    ! complexity nOrb * nOrb
    shXpYQ(:) = 0.0_dp
    if (cSym == "S") then
      shXpYQ(:) = matmul(vecXpYq, gamma)
    else
      shXpYQ(:) = 0.5_dp * vecXpYq(:) * (spinW(species0))
    end if

    ! calculate xpycc
    ! (xpycc)_{mu nu} = sum_{ia} (X + Y)_{ia} (c(mu,i)c(nu,a) + c(nu,i)c(mu,a))
    ! complexity nOrb * nOrb * nOrb
    vecXpYcc(:,:) = 0.0_dp
    vecXmYcc(:,:) = 0.0_dp
    vecXpYas(:,:) = 0.0_dp
    vecXmYas(:,:) = 0.0_dp
    ! xpycc(mu,nu) = sum_ia (X+Y)_ia c(mu,i) c(nu,a)
    ! xpycc(mu, nu) += sum_ia (X+Y)_ia c(mu,a) c(nu,i)
    do ia = 1, nXov
      call indXov(win, ia, getIJ, i, a)
      lUpDwn = (win(ia) <= nMatUp)
      do nu = 1, nOrb
        do mu = 1, nOrb
          vecXpYcc(mu, nu) = vecXpYcc(mu, nu) + vecXpY(ia) *&
              & (grndEigVecs(mu, i, 1) * grndEigVecs(nu, a, 1) +&
              &  grndEigVecs(mu, a, 1) * grndEigVecs(nu, i, 1))
          vecXmYcc(mu, nu) = vecXmYcc(mu, nu) + vecXmY(ia) *&
              & (grndEigVecs(mu, i, 1) * grndEigVecs(nu, a, 1) +&
              &  grndEigVecs(mu, a, 1) * grndEigVecs(nu, i, 1))
          vecXpYas(mu, nu) = vecXpYas(mu, nu) +&
              & vecXpY(ia) * grndEigVecs(mu, i, 1) * grndEigVecs(nu, a, 1)
          vecXmYas(mu, nu) = vecXmYas(mu, nu) +&
              & vecXmY(ia) * grndEigVecs(mu, i, 1) * grndEigVecs(nu, a, 1)
        end do
      end do
    end do

    ! calculate wcc = c_mu,i * W_ij * c_j,nu
    ! We have only W_ab b > a and W_ij j > i:
    !   wcc(m,n) = sum_{pq, p <= q} w_pq (c(mu,p)c(nu,q) + c(nu,p)c(mu,q))
    ! complexity nOrb * nOrb * nOrb

    ! calculate the occ-occ part
    wcc(:,:) = 0.0_dp
    do ij = 1, nXoo
      call indXoo(ij, i, j)
      do mu = 1, nOrb
        do nu = 1, nOrb
          wcc(mu, nu) = wcc(mu, nu) + vWoo(ij) *&
              & (grndEigVecs(mu, i, 1) * grndEigVecs(nu, j, 1) +&
              &  grndEigVecs(mu, j, 1)*grndEigVecs(nu, i, 1))
        end do
      end do
    end do

    ! calculate the occ-virt part : the same way as for xpycc
    do ia = 1, nXov
      call indXov(win, ia, getIJ, i, a)
      lUpDwn = (win(ia) <= nMatUp)
      do nu = 1, nOrb
        do mu = 1, nOrb
          wcc(mu, nu) = wcc(mu, nu) + vWov(ia) *&
              & (grndEigVecs(mu, i, 1) * grndEigVecs(nu, a, 1) +&
              &  grndEigVecs(mu, a, 1) * grndEigVecs(nu, i, 1))
        end do
      end do
    end do

    !calculate the virt - virt part
    do ab =1, nXvv
      call indXvv(homo, ab, a, b)
      do mu = 1, nOrb
        do nu = 1, nOrb
          wcc(mu, nu) = wcc(mu, nu) + vWvv(ab) *&
              & (grndEigVecs(mu, a, 1) * grndEigVecs(nu, b, 1) +&
              &  grndEigVecs(mu, b, 1) * grndEigVecs(nu, a, 1))
        end do
      end do
    end do

    ! now calculating the force !
    ! complexity : nOrb * nOrb * 3
    call symm(PS, 'R', overlap, pc, 'U', 1.0_dp, 0.0_dp, nOrb, nOrb)
    call symm(SPS, 'L', overlap, PS, 'U', 1.0_dp, 0.0_dp, nOrb, nOrb)
    call symm(DS, 'R', overlap, deltaRho, 'U', 1.0_dp, 0.0_dp, nOrb, nOrb)
    call symm(SDS, 'L', overlap, DS, 'U', 1.0_dp, 0.0_dp, nOrb, nOrb)
    call symm(XS, 'R', overlap, vecXpYas, 'U', 1.0_dp, 0.0_dp, nOrb, nOrb)
    call symm(SX, 'L', overlap, vecXpYas, 'U', 1.0_dp, 0.0_dp, nOrb, nOrb)
    call symm(SXS, 'L', overlap, XS, 'U', 1.0_dp, 0.0_dp, nOrb, nOrb)
    call symm(YS, 'R', overlap, vecXmYas, 'U', 1.0_dp, 0.0_dp, nOrb, nOrb)
    call symm(SY, 'L', overlap, vecXmYas, 'U', 1.0_dp, 0.0_dp, nOrb, nOrb)
    call symm(SYS, 'L', overlap, YS, 'U', 1.0_dp, 0.0_dp, nOrb, nOrb)

    !only for non-periodic systems!
    do alpha = 1, nBeweg
      indAlpha = iAtomStart(alpha)
      indAlpha1 = iAtomStart(alpha + 1) - 1
      izpAlpha = species0(alpha)
      do beta = 1, alpha - 1
        indBeta = iAtomStart(beta)
        indBeta1 = iAtomStart(beta + 1) - 1
        izpBeta = species0(beta)
        diffvec = coord0(:,alpha) - coord0(:,beta)
        rab = sqrt(sum(diffvec**2))
        ! Note: Here dgamma/dr * 1/r is needed as returned by old gam121
        tmp1 = (-1.0_dp / rab**2 - expGammaPrime(rab, HubbardU(izpAlpha), HubbardU(izpBeta))) / rab
        ! calculate the derivative of gamma
        dgab(:) = diffvec(:) * tmp1
        tmp3a(:) = dgab(:) * (dQ(alpha) * dQAtomEx(beta) + dQAtomEx(alpha) * dQ(beta))
        if (cSym == "S") then
          tmp3b(:) = 4.0_dp * dgab(:) * vecXpYq(alpha) * vecXpYq(beta)
        else
          tmp3b(:) = 0.0_dp
        end if
        excGrad(:,alpha) = excGrad(:,alpha) + tmp3a + tmp3b
        excGrad(:,beta) = excGrad(:,beta) - tmp3a - tmp3b
        tmp5 = shEx(beta) + shEx(alpha)
        tmp7 = 2.0_dp * shXpYQ(beta) + 2.0_dp * shXpYQ(alpha)

        tmprs = 0.0_dp
        do mu = indAlpha, indAlpha1
          do nu = indBeta, indBeta1
            tmprs = tmprs +&
                & (2.0_dp * (PS(mu,nu) * DS(nu,mu) + PS(nu,mu) * DS(mu,nu)) +&
                &   SPS(mu,nu) * deltaRho(mu,nu) + SPS(nu,mu) * deltaRho(nu,mu) +&
                &   pc(mu,nu) * SDS(mu,nu) + pc(nu,mu) * SDS(nu,mu))
            tmprs = tmprs + 2.0_dp *&
                & (vecXpYas(mu,nu) * SXS(mu,nu) + vecXpYas(nu,mu) * SXS(nu,mu) +&
                &   SX(mu,nu) * XS(mu,nu) + SX(nu,mu) * XS(nu,mu))
            tmprs = tmprs +&
                & (XS(mu,nu) * XS(nu,mu) + XS(nu,mu) * XS(mu,nu) +&
                &   SXS(mu,nu) * vecXpYas(nu,mu) + SXS(nu,mu) * vecXpYas(mu,nu) +&
                &   vecXpYas(mu,nu) * SXS(nu,mu) + vecXpYas(nu,mu) * SXS(mu,nu) +&
                &   SX(mu,nu) * SX(nu,mu) + SX(nu,mu) * SX(mu,nu))
            tmprs = tmprs + 2.0_dp *&
                & (vecXmYas(mu,nu) * SYS(mu,nu) + vecXmYas(nu,mu) * SYS(nu,mu) +&
                &   SY(mu,nu) * YS(mu,nu) + SY(nu,mu) * YS(nu,mu))
            tmprs = tmprs -&
                & (YS(mu,nu) * YS(nu,mu) + YS(nu,mu) * YS(mu,nu) +&
                &   SYS(mu,nu) * vecXmYas(nu,mu) + SYS(nu,mu) * vecXmYas(mu,nu) +&
                &   vecXmYas(mu,nu) * SYS(nu,mu) + vecXmYas(nu,mu) * SYS(mu,nu) +&
                &   SY(mu,nu) * SY(nu,mu) + SY(nu,mu) * SY(mu,nu))
          end do
        end do

        excGrad(:,alpha) = excGrad(:,alpha) - 0.125_dp * tmprs * gammaLongRangePrime(:,alpha,beta)
        excGrad(:,beta) = excGrad(:,beta) + 0.125_dp * tmprs * gammaLongRangePrime(:,alpha,beta)

        call derivator%getFirstDeriv(derH0, skHamCont, coord0, species0, alpha, beta, orb)
        call derivator%getFirstDeriv(derS, skOverCont, coord0, species0, alpha, beta, orb)

        do xyz = 1, 3
          tmp1 = 0.0_dp
          tmp2 = 0.0_dp
          tmp3 = 0.0_dp
          tmp4 = 0.0_dp
          tmp6 = 0.0_dp
          tmprs2 = 0.0_dp
          do mu = indAlpha, indAlpha1
            do nu = indBeta, indBeta1
              m = mu - indAlpha + 1
              n = nu - indBeta + 1

              dHMNdR = derH0(n,m,xyz)
              dSMNdR = derS(n,m,xyz)

              tmp1 = tmp1 + 2.0_dp * dHMNdR * pc(mu,nu)
              tmp2 = tmp2 + dSMNdR * pc(mu,nu) * (shift(alpha) + shift(beta))
              tmp3 = tmp3 - dSMNdR * wcc(mu,nu)
              tmp4 = tmp4 + tmp5 * dSMNdR * dmn(mu,nu)
              tmp6 = tmp6 + tmp7 * dSMNdR * vecXpYcc(mu,nu)

              tmprs = 0.0_dp
              do ka = 1, nOrb
                tmprs = tmprs +&
                    & (PS(mu,ka) * deltaRho(nu,ka) + PS(nu,ka) * deltaRho(mu,ka) +&
                    &  pc(mu,ka) * DS(nu,ka) + pc(nu,ka) * DS(mu,ka)) *&
                    & (lrGammaOrb(mu,ka) + lrGammaOrb(nu,ka))
                tmprs = tmprs +&
                    & (vecXpYas(mu,ka) * XS(nu,ka) + vecXpYas(ka,mu) * SX(ka,nu) +&
                    &  vecXpYas(nu,ka) * XS(mu,ka) + vecXpYas(ka,nu) * SX(ka,mu)) *&
                    & (lrGammaOrb(mu,ka) + lrGammaOrb(nu,ka))
                tmprs = tmprs +&
                    & (vecXmYas(mu,ka) * YS(nu,ka) + vecXmYas(ka,mu) * SY(ka,nu) +&
                    &  vecXmYas(nu,ka) * YS(mu,ka) + vecXmYas(ka,nu) * SY(ka,mu)) *&
                    & (lrGammaOrb(mu,ka) + lrGammaOrb(nu,ka))
                tmprs = tmprs +&
                    & (XS(mu,ka) * vecXpYas(ka,nu) + XS(nu,ka) * vecXpYas(ka,mu) +&
                    &  vecXpYas(mu,ka) * SX(ka,nu) + vecXpYas(nu,ka) * SX(ka,mu)) *&
                    & (lrGammaOrb(mu,ka) + lrGammaOrb(nu,ka))
                tmprs = tmprs -&
                    & (YS(mu,ka) * vecXmYas(ka,nu) + YS(nu,ka) * vecXmYas(ka,mu) +&
                    &  vecXmYas(mu,ka) * SY(ka,nu) + vecXmYas(nu,ka) * SY(ka,mu)) *&
                    & (lrGammaOrb(mu,ka) + lrGammaOrb(nu,ka))
              end do
              tmprs2 = tmprs2 + dSMNdR * tmprs
            end do
          end do

          excGrad(xyz,alpha) = excGrad(xyz,alpha) + tmp1 + tmp2 + tmp4 + tmp6 + tmp3&
              & - 0.25_dp * tmprs2
          excGrad(xyz,beta) = excGrad(xyz,beta) - tmp1 - tmp2 - tmp4 - tmp6 - tmp3&
              & + 0.25_dp * tmprs2
        end do
      end do
    end do

    deallocate(shEx)
    deallocate(vecXpYq)
    deallocate(shXpYQ)
    deallocate(wcc)
    deallocate(qIJ)
    deallocate(temp)
    deallocate(dmn)

    deallocate(vecXpYcc)
    deallocate(vecXmYcc)
    deallocate(vecXpYas)
    deallocate(vecXmYas)
    deallocate(gammaLongRangePrime)
    deallocate(overlap)
    deallocate(lrGammaOrb)
    deallocate(PS)
    deallocate(DS)
    deallocate(SPS)
    deallocate(SDS)
    deallocate(SX)
    deallocate(XS)
    deallocate(SXS)
    deallocate(SY)
    deallocate(YS)
    deallocate(SYS)

  end subroutine addGradientsRS


  !> Calculate oscillator strength for a given excitation range sep version
  subroutine getOscillatorStrengthsRS(cSym, snglPartTransDip, eval, vecXpY, occNr,&
      & istat, oscStrength, tTraDip, transDip)
    character, intent(in) :: cSym
    real(dp), intent(in) :: snglPartTransDip(:,:), eval(:), vecXpY(:,:)
    real(dp), intent(in) :: occNr(:,:)
    logical :: tTraDip
    integer, intent(inout) :: istat
    real(dp), intent(out) :: oscStrength(:), transDip(:,:)

    integer :: ii, nMat
    real(dp) :: oscStrength0
    logical :: spin

    nMat = size(vecXpY, dim=1)
    spin = .false.
    if (size(occNr, dim=2) == 2) spin = .true.

    ! Triplet oscillator strength and transition dipole is zero
    if ((.not. spin) .and. (cSym == "T")) then
      oscStrength = 0.0_dp
      if (tTraDip) transDip(:,:) = 0.0_dp
      return
    end if

    do ii = 1, size(vecXpY, dim=2)
      oscStrength(ii) = oscillatorStrengthRS(snglPartTransDip, sqrt(eval(ii)), vecXpY(:,ii))
    end do

    if (tTraDip) call transitionDipoleRS(snglPartTransDip, vecXpY, transDip)

    if (istat < 0) then
      istat = 1
      oscStrength0 = oscStrength(1)
      do ii = 2, size(oscStrength)
        if (oscStrength(ii) > oscStrength0) then
          oscStrength0 = oscStrength(ii)
          istat = ii
        end if
      end do
    end if

  contains

    pure function oscillatorStrengthRS(snglPartTransDip, omega, vecXpY) result(oscStrength)
      real(dp), intent(in) :: snglPartTransDip(:,:), omega, vecXpY(:)
      real(dp) :: oscStrength

      real(dp) :: rtmp
      integer :: ii

      oscStrength = 0.0_dp
      do ii = 1, 3
        rtmp = sum(snglPartTransDip(:,ii) * vecXpY)
        oscStrength = oscStrength + rtmp * rtmp
      end do
      oscStrength = twothird * 2.0_dp * omega * oscStrength

    end function oscillatorStrengthRS

    pure subroutine transitionDipoleRS(snglPartTransDip, vecXpY, transDip)
      real(dp), intent(in) :: snglPartTransDip(:,:), vecXpY(:,:)
      real(dp), intent(out) :: transDip(:,:)

      integer :: ii, ll

      transDip(:,:) = 0.0_dp
      do ii = 1, size(vecXpY, dim=2)
        do ll = 1, 3
          transDip(ii,ll) = sum(snglPartTransDip(:,ll) * sqrt(2.0_dp) * vecXpY(:,ii))
        end do
      end do

    end subroutine transitionDipoleRS

  end subroutine getOscillatorStrengthsRS


  !> Calculate transition charges for a given excitation
  !> (assumes singlet excitation for now, should set result to 0 for triplets)
  subroutine transitionChargesRS(vecXpY, tQov, atomicTransQ)
    real(dp), intent(in) :: vecXpY(:), tQov(:,:)
    real(dp), intent(out) :: atomicTransQ(:)
    real(dp) :: prefactor

    prefactor = sqrt(2.0_dp) ! Check sqrt(2.0) !-> should be fine
    atomicTransQ(:) = prefactor * matmul(tQov, vecXpY)

  end subroutine transitionChargesRS


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
      S(mu,mu) = 1.0_dp !Diagonal entries
    end do

  end subroutine getSqrS


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

  subroutine getHvvXY(ipm, nXov, nXvv, nOrb, homo, nAtom, nMatUp, iaTrans, getIJ, win,&
      & iAtomStart, sTimesGrndEigVecs, grndEigVecs, lrGamma, XorY, vecHvv)
    real(dp), intent(in) :: sTimesGrndEigVecs(:,:,:), grndEigVecs(:,:,:), lrGamma(:,:), XorY(:)
    real(dp), intent(out) :: vecHvv(:)
    real(dp), allocatable :: qIJ(:), gqIJ(:), qX(:,:), Gq(:,:)
    integer, intent(in) :: ipm, nXov, nXvv, nOrb, homo, nAtom, nMatUp, win(:), iAtomStart(:)
    integer, intent(in) :: iaTrans(1:,homo+1:), getIJ(:,:)

    integer :: i, a, b, ia, ib, ab, nXovAct, nOcc, nVirt
    logical :: lUpDwn

    nVirt = 0
    do while (nVirt*(nVirt+1)/2 < nXvv)
      nVirt = nVirt + 1
    end do
    nOcc = nOrb - nVirt

    nXovAct = nOcc * nVirt

    allocate(qIJ(nAtom))
    allocate(gqIJ(nAtom))
    allocate(qX(nAtom, nXovAct)) ! size(XorY))) ! allocate(qX(nAtom, nXov))
    allocate(Gq(nAtom, nXovAct)) ! size(XorY))) ! allocate(Gq(nAtom, nXov))

    qX(:,:) = 0.0_dp
    do ia = 1, nXovAct ! size(XorY) ! nXov
      call indXov(win, ia, getIJ, i, a)
      lUpDwn = (win(ia) <= nMatUp)
      do b = homo + 1, nOrb
        ib = iaTrans(i, b)
        if (ib <= nXov) then
          lUpDwn = .true.
          qIJ = transQ(a, b, iAtomStart, lUpDwn, sTimesGrndEigVecs, grndEigVecs)
          qX(:,ia) = qX(:,ia) + qIJ(:)*XorY(ib)
        end if
      end do
    end do

    Gq(:,:)  = 0.0_dp
    do ia = 1, nXovAct ! size(XorY)
      call indXov(win, ia, getIJ, i, a)
      lUpDwn = (win(ia) <= nMatUp)
      qIJ = transQ(i, a, iAtomStart, lUpDwn, sTimesGrndEigVecs, grndEigVecs)
      call dsymv('U', nAtom, 1.0_dp, lrGamma, nAtom, qIJ, 1, 0.0_dp, gqIJ, 1)
      Gq(:,ia) = gqIJ(:)
    end do

    vecHvv(:) = 0.0_dp
    do ab = 1, nXvv
      call indXvv(homo, ab, a, b)
      do i = 1, homo
        ia = iaTrans(i, a)
        ib = iaTrans(i, b)
        vecHvv(ab) = vecHvv(ab) - ipm * (dot_product(qX(:,ia), Gq(:,ib))&
            & + ipm * dot_product(Gq(:,ia), qX(:,ib)))
      end do
    end do

    deallocate(qIJ)
    deallocate(gqIJ)
    deallocate(qX)
    deallocate(Gq)

  end subroutine getHvvXY

  subroutine getHooXY(ipm, nXov, nXoo, nOrb, homo, nAtom, nMatUp, iaTrans, getIJ, win,&
      & iAtomStart, sTimesGrndEigVecs, grndEigVecs, lrGamma, XorY, vecHoo)
    real(dp), intent(in) :: sTimesGrndEigVecs(:,:,:), grndEigVecs(:,:,:), lrGamma(:,:), XorY(:)
    integer, intent(in) :: ipm, nXov, nXoo, nOrb, homo, nAtom, nMatUp, win(:), iAtomStart(:)
    integer, intent(in) :: iaTrans(1:,homo+1:), getIJ(:,:)
    real(dp), intent(out) :: vecHoo(:)

    real(dp), allocatable :: qIJ(:), gqIJ(:), qX(:,:), Gq(:,:)
    integer :: i, a, j, ia, ja, ij, nXovAct, nOcc, nVirt
    logical :: lUpDwn

    nOcc = 0
    do while (nOcc*(nOcc+1)/2 < nXoo)
      nOcc = nOcc + 1
    end do
    nVirt = nOrb - nOcc

    nXovAct = nOcc * nVirt

    allocate(qIJ(nAtom))
    allocate(gqIJ(nAtom))
    allocate(qX(nAtom, nXovAct)) ! size(XorY))) ! allocate(qX(nAtom, nXov))
    allocate(Gq(nAtom, nXovAct)) ! size(XorY))) ! allocate(Gq(nAtom, nXov))

    qX(:,:) = 0.0_dp
    do ia = 1, nXovAct ! size(XorY) ! nXov
      call indXov(win, ia, getIJ, i, a)
      lUpDwn = (win(ia) <= nMatUp)
      do j = 1, homo
        ja = iaTrans(j, a)
        if (ja <= nXov) then
          lUpDwn = .true.
          qIJ = transQ(i, j, iAtomStart, lUpDwn, sTimesGrndEigVecs, grndEigVecs)
          qX(:,ia) = qX(:,ia) + qIJ(:)*XorY(ja)
        end if
      end do
    end do

    Gq(:,:)  = 0.0_dp
    do ia = 1, nXovAct ! size(XorY) ! nXov
      call indXov(win, ia, getIJ, i, a)
      lUpDwn = (win(ia) <= nMatUp)
      qIJ = transQ(i, a, iAtomStart, lUpDwn, sTimesGrndEigVecs, grndEigVecs)
      call dsymv('U', nAtom, 1.0_dp, lrGamma, nAtom, qIJ, 1, 0.0_dp, gqIJ, 1)
      Gq(:,ia) = gqIJ(:)
    end do

    vecHoo(:) = 0.0_dp
    do ij = 1, nXoo
      call indXoo(ij, i, j)
      do a = homo + 1, nOrb
        ia = iaTrans(i, a)
        ja = iaTrans(j, a)
        vecHoo(ij) = vecHoo(ij) - ipm * (dot_product(qX(:,ia), Gq(:,ja))&
            & + ipm * dot_product(Gq(:,ia), qX(:,ja)))
      end do
    end do

    deallocate(qIJ)
    deallocate(gqIJ)
    deallocate(qX)
    deallocate(Gq)

  end subroutine getHooXY

  subroutine getHovT(nOcc, nXov, nXoo, nXvv, nOrb, homo, nAtom, nMatUp, iaTrans, getIJ, win,&
      & iAtomStart, sTimesGrndEigVecs, grndEigVecs, lrGamma, t, vecHov)
    real(dp), intent(in) :: sTimesGrndEigVecs(:,:,:), grndEigVecs(:,:,:), lrGamma(:,:), t(:,:)
    integer, intent(in) :: nOcc, nXov, nXoo, nXvv, nOrb, homo, nAtom, nMatUp, win(:), iAtomStart(:)
    integer, intent(in) :: iaTrans(1:,homo+1:), getIJ(:,:)
    real(dp), intent(out) :: vecHov(:)

    real(dp), allocatable :: qIJ(:), gqIJ(:), qX(:,:), Gq(:,:)
    integer :: i, a, b, j, ia, ib, ab, ja, ij
    logical :: lUpDwn

    allocate(qIJ(nAtom))
    allocate(gqIJ(nAtom))
    allocate(qX(nAtom, nOcc*(nOrb-nOcc))) ! allocate(qX(nAtom, nXov))
    allocate(Gq(nAtom, max(nXoo,nXvv)))

    qX(:,:) = 0.0_dp
    do ia = 1, nOcc*(nOrb-nOcc)
      call indXov(win, ia, getIJ, i, a)
      lUpDwn = (win(ia) <= nMatUp)
      do b = homo + 1, nOrb
        lUpDwn = (win(iaTrans(i, b)) <= nMatUp) ! UNTESTED
        qIJ = transQ(i, b, iAtomStart, lUpDwn, sTimesGrndEigVecs, grndEigVecs)
        qX(:,ia) = qX(:,ia) + qIJ(:)*t(a,b)
      end do
    end do

    Gq(:,:)  = 0.0_dp
    do ab = 1, nXvv
      call indXvv(homo, ab, a, b)
      lUpDwn = .true.
      qIJ = transQ(a, b, iAtomStart, lUpDwn, sTimesGrndEigVecs, grndEigVecs)
      call dsymv('U', nAtom, 1.0_dp, lrGamma, nAtom, qIJ, 1, 0.0_dp, gqIJ, 1)
      Gq(:,ab) = gqIJ(:)
    end do

    vecHov(:) = 0.0_dp
    do ia = 1, nXov ! nOcc*(nOrb-nOcc)
      call indXov(win, ia, getIJ, i, a)
      lUpDwn = (win(ia) <= nMatUp)
      do b = homo + 1, nOrb
        ib = iaTrans(i, b)
        if (b < a) then
          call rIndXvv(homo, a, b, ab)
        else
          call rIndXvv(homo, b, a, ab)
        end if
        vecHov(ia) = vecHov(ia) - 2.0_dp * dot_product(qX(:,ib), Gq(:,ab))
      end do
    end do

    qX(:,:) = 0.0_dp
    do ia = 1, nOcc*(nOrb-nOcc) ! nXov
      call indXov(win, ia, getIJ, i, a)
      lUpDwn = (win(ia) <= nMatUp)
      do j = 1, homo
        lUpDwn = (win(iaTrans(j, a)) <= nMatUp) ! UNTESTED
        qIJ = transQ(j, a, iAtomStart, lUpDwn, sTimesGrndEigVecs, grndEigVecs)
        qX(:,ia) = qX(:,ia) + qIJ(:)*t(i,j)
      end do
    end do

    Gq(:,:)  = 0.0_dp
    do ij = 1, nXoo
      call indXoo(ij, i, j)
      lUpDwn = .true.
      qIJ = transQ(i, j, iAtomStart, lUpDwn, sTimesGrndEigVecs, grndEigVecs)
      call dsymv('U', nAtom, 1.0_dp, lrGamma, nAtom, qIJ, 1, 0.0_dp, gqIJ, 1)
      Gq(:,ij) = gqIJ(:)
    end do

    do ia = 1, nXov ! nOcc*(nOrb-nOcc)
      call indXov(win, ia, getIJ, i, a)
      lUpDwn = (win(ia) <= nMatUp)
      do j = 1, homo
        ja = iaTrans(j, a)
        if (j < i) then
          call rIndXoo(homo, nOcc, i, j, ij)
        else
          call rIndXoo(homo, nOcc, j, i, ij)
        end if
        vecHov(ia) = vecHov(ia) - 2.0_dp * dot_product(qX(:,ja), Gq(:,ij))
      end do
    end do

    deallocate(qIJ)
    deallocate(gqIJ)
    deallocate(qX)
    deallocate(Gq)

  end subroutine getHovT

  subroutine getHooT(nOcc, nXov, nXoo, nOrb, homo, nAtom, nMatUp, iaTrans, getIJ, win,&
      & iAtomStart, sTimesGrndEigVecs, grndEigVecs, lrGamma, t, vecHoo)
    real(dp), intent(in) :: sTimesGrndEigVecs(:,:,:), grndEigVecs(:,:,:), lrGamma(:,:), t(:,:)
    integer, intent(in) :: nOcc, nXov, nXoo, nOrb, homo, nAtom, nMatUp, win(:), iAtomStart(:)
    integer, intent(in) :: iaTrans(1:,homo+1:), getIJ(:,:)
    real(dp), intent(out) :: vecHoo(:)

    real(dp), allocatable :: qIJ(:), gqIJ(:), qX(:,:), qXa(:,:,:), Gq(:,:)
    integer :: i, a, b, j, k, ia, ja, ij, jk, nXovAct
    logical :: lUpDwn

    nXovAct = nOcc*(nOrb-nOcc)

    allocate(qIJ(nAtom))
    allocate(gqIJ(nAtom))
    allocate(qX(nAtom, nXovAct)) ! nXov))
    allocate(Gq(nAtom, max(nXoo, nXovAct))) ! nXov)))

    qX(:,:) = 0.0_dp
    do ia = 1, nXovAct ! nXov
      call indXov(win, ia, getIJ, i, a)
      lUpDwn = (win(ia) <= nMatUp)
      do b = homo + 1, nOrb
        lUpDwn = (win(iaTrans(i, b)) <= nMatUp) ! UNTESTED
        qIJ = transQ(i, b, iAtomStart, lUpDwn, sTimesGrndEigVecs, grndEigVecs)
        qX(:,ia) = qX(:,ia) + qIJ(:)*t(a,b)
      end do
    end do

    Gq(:,:)  = 0.0_dp
    do ia = 1, nXovAct ! nXov
      call indXov(win, ia, getIJ, i, a)
      lUpDwn = (win(ia) <= nMatUp)
      qIJ = transQ(i, a, iAtomStart, lUpDwn, sTimesGrndEigVecs, grndEigVecs)
      call dsymv('U', nAtom, 1.0_dp, lrGamma, nAtom, qIJ, 1, 0.0_dp, gqIJ, 1)
      Gq(:,ia) = gqIJ(:)
    end do

    vecHoo(:) = 0.0_dp
    do ij = 1, nXoo
      call indXoo(ij, i, j)
      do a = homo + 1, nOrb
        ia = iaTrans(i, a)
        ja = iaTrans(j, a)
        vecHoo(ij) = vecHoo(ij) - 2.0_dp * dot_product(qX(:,ia), Gq(:,ja))
      end do
    end do
    deallocate(qX)

    allocate(qXa(nAtom, nOcc, nOcc))
    qXa(:,:,:) = 0.0_dp
    do i = 1, homo
      do k = 1, homo
        lUpDwn = .true.
        qIJ = transQ(i, k, iAtomStart, lUpDwn, sTimesGrndEigVecs, grndEigVecs)
        do j = 1, homo
          qXa(:,i,j) = qXa(:,i,j) + qIJ(:)*t(j,k)
        end do
      end do
    end do

    Gq(:,:)  = 0.0_dp
    do ij = 1, nXoo
      call indXoo(ij, i, j)
      lUpDwn = .true.
      qIJ = transQ(i, j, iAtomStart, lUpDwn, sTimesGrndEigVecs, grndEigVecs)
      call dsymv('U', nAtom, 1.0_dp, lrGamma, nAtom, qIJ, 1, 0.0_dp, gqIJ, 1)
      Gq(:,ij) = gqIJ(:)
    end do

    do ij = 1, nXoo
      call indXoo(ij, i, j)
      do k = 1, homo
        if (k < j) then
          call rIndXoo(homo, nOcc, j, k, jk)
        else
          call rIndXoo(homo, nOcc, k, j, jk)
        end if
        vecHoo(ij) = vecHoo(ij) - 2.0_dp * dot_product(qXa(:,i,k), Gq(:,jk))
      end do
    end do

    deallocate(qIJ)
    deallocate(gqIJ)
    deallocate(qXa)
    deallocate(Gq)

  end subroutine getHooT

end module dftbp_rs_linearresponse
