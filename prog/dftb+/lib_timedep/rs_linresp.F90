!!
!!This module heavily relies on the linresp module by Dominguez et al. 
!!and the rangesep module by Lusker et al.
!!Periodic systems and 3rd order calculations are not supported so far.

#:include 'common.fypp'

module dftbp_rs_linearresponse
  use dftbp_assert
  use dftbp_arpack
  use dftbp_commontypes
  use dftbp_slakoCont
  use dftbp_shortgamma
  use dftbp_accuracy
  use dftbp_constants, only : Hartree__eV, au__Debye
  use dftbp_nonscc, only : NonSccDiff
  use dftbp_scc, only: TScc
  use dftbp_blasroutines
  use dftbp_eigensolver
  use dftbp_message
  use dftbp_taggedOutput, only : TTaggedWriter, tagLabels
  use dftbp_linresp
  use dftbp_linrespcommon
  use dftbp_linrespgrad,!only: calcTransitionDipoles, getExcSpin, writeSPExcitations, writeExcMulliken
  use dftbp_rangeseparated
  use dftbp_sorting
  use dftbp_qm
  use dftbp_transcharges
  use dftbp_sk, only: rotateH0
  implicit none 
  private
  !!^-- Check: once done, see which modules are actually required

  public :: LinResp_calcExcitations_rs

  type :: rs_linresp
     !> type to hold data required only for rangesep linresp but not linresp itself
     !Check: If they turn out to be feq, consider adding them to linresp itself

  end type rs_linresp

  interface incSize
     module procedure incSizeInt
     module procedure incSizeVec
     module procedure incSizeMat
  end interface 

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
  character(*), parameter :: transChargesOut = "ATQ.DAT"
  character(*), parameter :: singlePartOut = "SPX.DAT"

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
 !subroutine rindxoo(nocc, nocc_r, ii, jj, index)
 !  integer, intent(in) ::  nocc, nocc_r, ii, jj
 !  integer, intent(out) :: index

 !  if (jj > ii) then  
 !     index = ii - nocc + nocc_r + ((jj - nocc + nocc_r - 1) * (jj - nocc + nocc_r)) / 2   
 !  else
 !     index = jj - nocc + nocc_r + ((ii - nocc + nocc_r - 1) * (ii - nocc + nocc_r)) / 2    
 !  end if

 !end subroutine rindxoo

  !> Computes the occupied-occupied compound index from two individual indices; assume jj <= ii
  pure subroutine rindxoo(nocc, nocc_r, ii, jj, indx)

    !> Number of occupied states.
    integer, intent(in) :: nocc

    !> Number of "active" occupied states.
    integer, intent(in) :: nocc_r

    !> First (initial) occupied state.
    integer, intent(in) :: ii

    !> Second (final) occupied state.
    integer, intent(in) :: jj

    !> Compound indx.
    integer, intent(out) :: indx

    indx = jj - nocc + nocc_r + ((ii - nocc + nocc_r - 1) * (ii - nocc + nocc_r)) / 2

  end subroutine rindxoo

  subroutine getSOffsite(coords1, coords2, iSp1, iSp2, orb, skOverCont, Sblock)
    implicit none
    real(dp), intent(in) :: coords1(:), coords2(:)
    integer, intent(in) :: iSp1, iSp2
    type(TOrbitals), intent(in) :: orb
    type(OSlakoCont), intent(in) :: skOverCont
    real(dp), intent(out) :: Sblock(:,:)

    real(dp) :: interSKOver(getMIntegrals(skOverCont))
    real(dp) :: vect(3), dist

    @:ASSERT(size(coords1) == 3)
    @:ASSERT(size(coords2) == 3)
    @:ASSERT(all(shape(Sblock) >= [ orb%mOrb, orb%mOrb ]))

    vect(:) = coords2 - coords1
    dist = sqrt(sum(vect**2))
    vect(:) = vect(:) / dist
    call getSKIntegrals(skOverCont, interSKOver, dist, iSp1, iSp2)
    call rotateH0(Sblock, interSKOver, vect(1), vect(2), vect(3), iSp1, iSp2, orb)

  end subroutine getSOffsite


  !> Set value of sizeIn to new, twice as large, value
  pure subroutine incSizeInt(sizeIn)
  
    implicit none
    integer, intent(inout) :: sizeIn
    sizeIn = 3 * sizeIn

  end subroutine incSizeInt


  !>increase size of a NxsizeIn array to 2*sizeIn
  pure subroutine incSizeMat(sizeIn, vec)
  
    implicit none
    integer, intent(inout) :: sizeIn
    real(dp), allocatable, intent(inout) :: vec(:,:)
    
    integer :: dim1
    real(dp), allocatable :: temp(:,:)
    dim1 = size(vec, dim=1)
    allocate(temp(dim1, 3 * sizeIn))
    temp(:,1:sizeIn) = vec   !Check: would be nice if temp could become memory for vec, see if possible in fortran
    call move_alloc(temp, vec)
    
  end subroutine incSizeMat


  !>increase size of a sizeIn array to 2*sizeIn
  pure subroutine incSizeVec(sizeIn, vec)
  
    implicit none
    integer, intent(in) :: sizeIn
    real(dp), allocatable, intent(inout) :: vec(:)
    
    real(dp), allocatable :: temp(:)
    allocate(temp(3 * sizeIn))
    temp(1:sizeIn) = vec   !Check: would be nice if temp could become memory for vec, see if possible in fortran
    call move_alloc(temp, vec)
    
  end subroutine incSizeVec


  !>increase size of a sizeInxN array to 2*sizeInxN
  pure subroutine incSizeMatSwapped(sizeIn, vec)
  
    implicit none
    integer, intent(inout) :: sizeIn
    real(dp), allocatable, intent(inout) :: vec(:,:)
    
    integer :: dim1
    real(dp), allocatable :: temp(:,:)
    dim1 = size(vec, dim=2)
    allocate(temp(3 * sizeIn, dim1))
    temp(1:sizeIn,:) = vec   !Check: would be nice if temo could become memory for vec, see if possible in fortran
    call move_alloc(temp, vec)
    
  end subroutine incSizeMatSwapped


  !>increase size of a sizeInxsizeIn array to 2*sizeIn
  pure subroutine incSizeMatBoth(sizeIn, vec)
  
    implicit none
    integer, intent(inout) :: sizeIn
    real(dp), allocatable, intent(inout) :: vec(:,:)
    
    real(dp), allocatable :: temp(:,:)
    
    allocate(temp(3 * sizeIn, 2 * sizeIn))
    temp(1:sizeIn,1:sizeIn) = vec   !Check: would be nice if temo could become memory for vec, see if possible in fortran
    call move_alloc(temp, vec)
    
  end subroutine incSizeMatBoth


  !>Calculate square root and inverse of sqrt of a real, symmetric positive definite matrix
  subroutine calcMatrixSqrt(MatIn, spaceDim, memDim, workArray, workDim, MatOut, MatInvOut)

    implicit none
    real(dp), intent(in) :: MatIn(:,:)
    integer, intent(in) :: spaceDim, memDim
    real(dp), intent(out) :: MatOut(:,:), MatInvOut(:,:)
    real(dp) :: workArray(:)
    integer :: workDim
    
    real(dp) :: dummyEV(spaceDim)
    real(dp) :: dummyM(spaceDim, spaceDim), dummyM2(spaceDim, spaceDim)
    integer :: info
    integer :: ii

    dummyM(:,:) = MatIn(1:spaceDim,1:spaceDim)
 !   write(*,*) ,'Debug mmm ', dummyM(:,:)

    call dsyev('V', 'U', spaceDim, dummyM, spaceDim, dummyEV, workArray, workDim, info)
    !calc. sqrt
    do ii=1, spaceDim
       dummyM2(:,ii)=sqrt(dummyEV(ii))*dummyM(:,ii)
    end do

!    write(*,*), 'debug ev ', dummyEV(:)
    call dgemm('N', 'T', spaceDim, spaceDim, spaceDim, 1.0_dp, dummyM2, spaceDim, dummyM, spaceDim, 0.0_dp, MatOut, memDim)
    !calc. inv. of sqrt
        do ii=1, spaceDim
       dummyM2(:,ii) = dummyM(:,ii) / sqrt(dummyEV(ii))
    end do

    call dgemm('N', 'T', spaceDim, spaceDim, spaceDim, 1.0_dp, dummyM2, spaceDim, dummyM, spaceDim, 0.0_dp, MatInvOut, memDim)

  end subroutine calcMatrixSqrt


  !>Perform modified Gram-Schmidt orthonormalization of vectors in columns of vec(1:end). Assume vectors 1:(start-1) are orthonormal yet
  pure subroutine orthonormalizeVectors(start, end, vec)
    implicit none
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

  !>Calculate the product of the matrix A+B and a vector. A,B as usual in linear response TD-DFT.
 !subroutine multApBVecFast(spin, vin, wij, sym, win, nMatUp, homo, nocc, nvir, ind, &
 !    & sTimesGrndEigVecs, grndEigVecs, occNr, getij, gamma, lrGamma, species0, uhbuu, uhbud, orb, iatrans, &
 !    & gqvtmp, tqov, tqoo, tqvv, vout)
  subroutine multApBVecFast(vin, wij, sym, win, nMatUp, homo, nocc, nvir, & ! sTimesGrndEigVecs,
      & occNr, getij, gamma, lrGamma, species0, uhbuu, uhbud, iatrans, &
      & gqvtmp, tqov, tqoo, tqvv, vout)
    implicit none
   !logical, intent(in)   :: spin !for now ignore spin and assume unpolarized system
    real(dp), intent(in)  :: vin(:)
    real(dp), intent(in)  :: wij(:)
    character, intent(in) :: sym
    integer, intent(in)   :: win(:), nMatUp, homo, nocc, nvir !, ind(:)
   !real(dp), intent(in)  :: sTimesGrndEigVecs(:,:,:), grndEigVecs(:,:,:)
    real(dp), intent(in)  :: occNr(:,:)
    real(dp), intent(in)  :: gamma(:,:), lrGamma(:,:)
    integer,  intent(in)  :: getij(:,:)     
    integer, intent(in)   :: species0(:)
    real(dp), intent(in)  :: uhbuu(:), uhbud(:)
   !type(TOrbitals), intent(in) :: orb    
    real(dp), intent(in)  :: tqov(:,:), tqoo(:,:), tqvv(:,:)
    integer, intent(in)   :: iatrans(1:,homo+1:)
    real(dp), intent(out) :: vout(:)
    real(dp), intent(out) :: gqvtmp(:,:)   !for faster, more memory intensive routine

    integer :: izpalpha, nmat, natom
    integer :: ia, ib, ja, jb, alpha, ii, jj, ij, aa, bb, ab
    real(dp) :: tmp2
    real(dp) :: otmp(size(gamma, dim=1)), gtmp(size(gamma, dim=1))
   !real(dp) :: qij(size(gamma, dim=1))
    real(dp) :: wnij(size(wij))
   !logical :: lUpDwn

    nmat = size(vin) ! also known as nxov
    natom = size(gamma, dim=1)
    gtmp(:) = 0.0_dp  !used to be otmp before optimization

    call wtdn(wij, occNr, win, nMatUp, nmat, getij, wnij)    

    !----only spin unpolarized case for now-------------------------------
    
    if (sym == 'S') then
       !full range coupling matrix contribution
       do ia = 1, nmat  
          !call hemv(gtmp, gamma, tqov(:,ia)) 
          !otmp(:) = otmp + vin(ia) * gtmp
          gtmp(:) = gtmp + tqov(:,ia) * vin(ia) !gtmp=sum_jb q_B^jb * vin_jb

       end do
       
       call hemv(otmp, gamma, gtmp)

       do ia= 1, nmat
          !call indxov(win, ia, getij, ii, jj)
          !lUpDwn = (win(ia) <= nMatUp)
          !qij = transq(ii, jj, ind, lUpDwn, sTimesGrndEigVecs, grndEigVecs)
          vout(ia) = 4.0 * dot_product(tqov(:,ia), otmp) !Check: prefactor! need 4.0 * for symmetry reduced sum? 
       end do

    else
       do ia = 1, nmat
          !call indxov(win, ia, getij, ii, jj)
          !lUpDwn = (win(ia) <= nMatUp)
          !qij = transq(ii, jj, ind, lUpDwn, sTimesGrndEigVecs, grndEigVecs)
          otmp(:) = otmp + vin(ia) * tqov(:,ia)
       end do
       
       do ia = 1, nmat
          vout(ia) = 0.0_dp
          !call indxov(win, ia, getij, ii, jj)
          !lUpDwn = (win(ia) <= nMatUp)
          !qij = transq(ii, jj, ind, lUpDwn, sTimesGrndEigVecs, grndEigVecs)
          do alpha = 1, natom
             izpalpha = species0(alpha)
             tmp2 = 2.0_dp * (uhbuu(izpalpha) - uhbud(izpalpha)) !== 2*magnetization
             vout(ia) = vout(ia) + otmp(alpha) * tmp2 * tqov(alpha,ia) 
          end do
       end do
  
    end if

!!!!Old implementation, probably slower, but need less memory 
!!!!Reintroduce as option later on

    ! !long range coupling matrix contribution
    ! do ia = 1, nmat   !adds -sum_AB ( q_ij^A * q_ab^B * Gamma_lr_AB + q_ib^A * q_ja^B * Gamma_lr_AB  )
    !    do jb = 1, nmat
         
    !       call indxov(win, ia, getij, ii, aa)
    !       call indxov(win, jb, getij, jj, bb)  

    !       lUpDwn = (win(iatrans(jj, aa)) <= nMatUp) ! UNTESTED
    !       qij = transq(jj, aa, ind, lUpDwn, sTimesGrndEigVecs, grndEigVecs)
    !       call hemv(gtmp, lrGamma, qij)

    !       lUpDwn = (win(iatrans(ii, bb)) <= nMatUp) ! UNTESTED
    !       qij = transq(ii, bb, ind, lUpDwn, sTimesGrndEigVecs, grndEigVecs)
    !       vout(ia) = vout(ia) - dot_product(qij, gtmp) * vin(jb)  

    !       lUpDwn = .true. ! ???
    !       qij = transq(ii, jj, ind, lUpDwn, sTimesGrndEigVecs, grndEigVecs)
    !       call hemv(gtmp, lrGamma, qij)

    !       lUpDwn = .true. ! ???
    !       qij = transq(aa, bb, ind, lUpDwn, sTimesGrndEigVecs, grndEigVecs)
    !       vout(ia) = vout(ia) - dot_product(qij, gtmp) * vin(jb)          

    !    end do
    ! end do

!!!!End old implementation

    !!!Faster but more memory intensive routine
   
    gqvtmp(:,:) = 0.0_dp
    do jb = 1, nmat
       call indxov(win, jb, getij, jj, bb)      
       !lUpDwn = (win(jb) <= nMatUp)
       do aa = homo + 1, homo + nvir 
          !lUpDwn = .true.
          !qij = transq(aa, bb, ind, lUpDwn, sTimesGrndEigVecs, grndEigVecs)
          if (bb < aa) then
             call rindxvv(homo, aa, bb, ab)
          else
             call rindxvv(homo, bb, aa, ab)
          end if
          !call hemv(gtmp, lrGamma, tqvv(:,ab))
          ja = nvir * (jj - 1) + aa - nocc
          gqvtmp(:,ja) = gqvtmp(:,ja) + tqvv(:,ab) * vin(jb)

       end do
    end do

    do ja = 1, nmat
       gtmp(:) = gqvtmp(:,ja)
       call hemv(gqvtmp(:,ja), lrGamma, gtmp)
    end do

    do jj = 1, nocc
       do ia = 1, nmat
          call indxov(win, ia, getij, ii, aa)  !would prefer to swap loops, but want to avoid calculation ia as unsure how implemted
          !lUpDwn = (win(ia) <= nMatUp)
          !lUpDwn = .true.
          !qij = transq(ii, jj, ind, lUpDwn, sTimesGrndEigVecs, grndEigVecs)
          if (jj < ii) then
             ij = jj + ((ii - 1 - homo + nocc) * (ii -homo + nocc)) / 2 + homo - nocc !inverted indxoo fct.
          else
             ij = ii + ((jj - 1 - homo + nocc) * (jj -homo + nocc)) / 2 + homo - nocc 
          end if
          ja = nvir * (jj - 1) + aa - nocc
          vout(ia) = vout(ia) - dot_product(tqoo(:,ij), gqvtmp(:,ja))

       end do
    end do
  
    gqvtmp(:,:) = 0.0_dp
    do jb = 1, nmat
       call indxov(win, jb, getij, jj, bb)      
       do aa = homo + 1, homo + nvir
          ab = (bb - nocc - 1) * nvir + aa - nocc
          ja = iatrans(jj, aa)
         !if (ja <= nmat) then ! TOMAS KUBAR ADDED
             gqvtmp(:,ab) = gqvtmp(:,ab) + tqov(:,ja) * vin(jb)
         !end if
       end do
    end do

    do ab = 1, nvir * nvir
       gtmp(:) = gqvtmp(:,ab)
       call hemv(gqvtmp(:,ab), lrGamma, gtmp) 
    end do

    do bb = homo + 1, homo + nvir
       do ia = 1, nmat
          call indxov(win, ia, getij, ii, aa)
          !lUpDwn = (win(ia) <= nMatUp)
          ib = iatrans(ii,bb)
         !if (ib <= nmat) then ! TOMAS KUBAR ADDED
             !lUpDwn = (win(jb) <= nMatUp)
             !qij = transq(ii, bb, ind, lUpDwn, sTimesGrndEigVecs, grndEigVecs)
             ab = (bb - nocc - 1) * nvir + aa - nocc 
             vout(ia) = vout(ia) - dot_product(tqov(:,ib), gqvtmp(:,ab))
         !end if
       end do
    end do
  
    !!!!end faster implementation

    !----only spin unpolarized case for now------------------------------------
   !vout(:) = vout + 0.5_dp * wnij * vin !orb. energy difference diagonal contribution
    vout(:) = vout + 0.5_dp * wnij(1:nmat) * vin(1:nmat)
    !Check: What about the factor 0.5 introduced here? Check derivation!

  end subroutine multApBVecFast


  !>Calculate the product of A-B and a vector.
 !subroutine multAmBVecFast(spin, vin, wij, sym, win, nMatUp, homo, nocc, nvir, ind, &
 !    & sTimesGrndEigVecs, grndEigVecs, occNr, getij, gamma, lrGamma, species0, orb, iatrans, gqvtmp, tqov, tqoo, tqvv, vout)
  subroutine multAmBVecFast(vin, wij, win, nMatUp, homo, nocc, nvir, &
      & occNr, getij, gamma, lrGamma, iatrans, gqvtmp, tqov, tqoo, tqvv, vout)
    implicit none
   !logical, intent(in) :: spin !for now ignore spin and assume unpolarized system
    real(dp), intent(in) :: vin(:)
    real(dp), intent(in) :: wij(:)
   !character, intent(in) :: sym
    integer, intent(in) :: win(:), nMatUp, homo, nocc, nvir !, ind(:)
   !real(dp), intent(in) :: sTimesGrndEigVecs(:,:,:), grndEigVecs(:,:,:)
    real(dp), intent(in) :: occNr(:,:)
    real(dp), intent(in) :: gamma(:,:), lrGamma(:,:)
    integer,  intent(in) :: getij(:,:)     
   !integer, intent(in) :: species0(:)
   !type(TOrbitals), intent(in) :: orb
    real(dp), intent(in) :: tqov(:,:), tqoo(:,:), tqvv(:,:)    
    integer, intent(in) :: iatrans(1:,homo+1:)
    real(dp), intent(out) :: vout(:)

    real(dp) :: gqvtmp(:,:)   !for faster, more memory intensive routine
    integer :: nmat, natom
    integer :: ia, ib, ja, jb, ii, jj, ij, aa, bb, ab
   !real(dp) :: otmp(size(gamma, dim=1))
    real(dp) :: gtmp(size(gamma, dim=1))    
   !real(dp) :: qij(size(gamma, dim=1))
    real(dp) :: wnij(size(wij))
   !logical :: lUpDwn

    nmat = size(vin) ! also known as nxov
    natom = size(gamma, dim=1)

    call wtdn(wij, occNr, win, nMatUp, nmat, getij, wnij)    
    
    !----only spin unpolarized case and singlet for now-------------------------------

!!!!Old routine, reintroduce as option later on, see A+B fct

    ! !long range coupling matrix contribution
    ! do ia = 1, nmat   !calcs -sum_AB ( q_ij^A * q_ab^B * Gamma_lr_AB - q_ib^A * q_ja^B * Gamma_lr_AB  )
    !    vout(ia) = 0.0_dp
    !    do jb = 1, nmat

    !       call indxov(win, ia, getij, ii, aa)
    !       call indxov(win, jb, getij, jj, bb)  

    !       lUpDwn = (win(iatrans(jj, aa)) <= nMatUp) ! UNTESTED
    !       qij = transq(jj, aa, ind, lUpDwn, sTimesGrndEigVecs, grndEigVecs)
    !       call hemv(gtmp, lrGamma, qij)

    !       lUpDwn = (win(iatrans(ii, bb)) <= nMatUp) ! UNTESTED
    !       qij = transq(ii, bb, ind, lUpDwn, sTimesGrndEigVecs, grndEigVecs)
    !       vout(ia) = vout(ia) + dot_product(qij, gtmp) * vin(jb)  

    !       lUpDwn = .true.
    !       qij = transq(ii, jj, ind, lUpDwn, sTimesGrndEigVecs, grndEigVecs)
    !       call hemv(gtmp, lrGamma, qij)

    !       lUpDwn = .true.
    !       qij = transq(aa, bb, ind, lUpDwn, sTimesGrndEigVecs, grndEigVecs)
    !       vout(ia) = vout(ia) - dot_product(qij, gtmp) * vin(jb)          

    !    end do
    ! end do

!!!End old routine

    !!!Faster but more memory intensive routine
   
    gqvtmp(:,:) = 0.0_dp
    do jb = 1, nmat
       call indxov(win, jb, getij, jj, bb)      
      !lUpDwn = (win(jb) <= nMatUp)
       do aa = homo + 1, homo + nvir 
          !lUpDwn = .true.
          !qij = transq(aa, bb, ind, lUpDwn, sTimesGrndEigVecs, grndEigVecs)
          if (bb < aa) then
             call rindxvv(homo, aa, bb, ab)
          else
             call rindxvv(homo, bb, aa, ab)
          end if
         ! call hemv(gtmp, lrGamma, tqvv(:,ab))
          ja = nvir * (jj - 1) + aa - nocc
          gqvtmp(:,ja) = gqvtmp(:,ja) + tqvv(:,ab) * vin(jb)
       end do
    end do

    do ja = 1, nmat
       gtmp(:) = gqvtmp(:,ja)
       call hemv(gqvtmp(:,ja), lrGamma, gtmp)
    end do

    do jj = 1, nocc
       do ia = 1, nmat
          call indxov(win, ia, getij, ii, aa)  !would prefer to swap loops, but want to avoid calculation ia as unsure how implemted
          !lUpDwn = (win(ia) <= nMatUp)
          !lUpDwn = .true.
          !qij = transq(ii, jj, ind, lUpDwn, sTimesGrndEigVecs, grndEigVecs)
          if (jj < ii) then
             ij = jj + ((ii - 1 - homo + nocc) * (ii -homo + nocc)) / 2 + homo - nocc !inverted indxoo fct.
          else
             ij = ii + ((jj - 1 - homo + nocc) * (jj -homo + nocc)) / 2 + homo - nocc 
          end if
          ja = nvir * (jj - 1) + aa - nocc
          vout(ia) = vout(ia) - dot_product(tqoo(:,ij), gqvtmp(:,ja))

       end do
    end do

    gqvtmp(:,:) = 0.0_dp
    do jb = 1, nmat
       call indxov(win, jb, getij, jj, bb)      
       do aa = homo + 1, homo + nvir
          ab = (bb - nocc - 1) * nvir + aa - nocc
          ja = iatrans(jj, aa)
         !if (ja <= nmat) then ! TOMAS KUBAR ADDED
             gqvtmp(:,ab) = gqvtmp(:,ab) + tqov(:,ja) * vin(jb)
         !end if
       end do
    end do

    do ab = 1, nvir * nvir
       gtmp(:) = gqvtmp(:,ab)
       call hemv(gqvtmp(:,ab), lrGamma, gtmp)
    end do

    do bb = homo + 1, homo + nvir
       do ia = 1, nmat
          call indxov(win, ia, getij, ii, aa)
          !lUpDwn = (win(ia) <= nMatUp)
          !lUpDwn = (win(iatrans(ii, bb)) <= nMatUp) ! UNTESTED
          !qij = transq(ii, bb, ind, lUpDwn, sTimesGrndEigVecs, grndEigVecs)
          ib = iatrans(ii,bb)
         !if (ib <= nmat) then ! TOMAS KUBAR ADDED
             ab = (bb - nocc - 1) * nvir + aa - nocc 
             vout(ia) = vout(ia) + dot_product(tqov(:,ib), gqvtmp(:,ab))
         !end if
       end do
    end do

    !!!!end faster implementation

   !vout(:) = vout + 0.5_dp * wnij * vin ! orb. energy difference diagonal contribution
    vout(:) = vout + 0.5_dp * wnij(1:nmat) * vin(1:nmat)
    !Check: does wnij really contain the right thing?

  end subroutine multAmBVecFast

!>Generates initial matrices Mm and Mp without calculating full Mat Vec product, not required for init. space.
! Use precalculated transition charges
  subroutine setupInitMatFast(initDim, & ! spin,
      & wij, sym, win, nMatUp, nocc, homo, & ! ind, sTimesGrndEigVecs, grndEigVecs,
      & occNr, getij, iatrans, gamma, lrGamma, species0, uhbuu, uhbud, & ! orb,
      & tqov, tqoo, tqvv, vp, vm, Mp,Mm)
    implicit none
    integer, intent(in) :: initDim
   !logical, intent(in) :: spin !for now ignore spin and assume unpolarized system
    real(dp), intent(in) :: wij(:)
    character, intent(in) :: sym
    integer, intent(in) :: win(:), nMatUp, nocc, homo ! , ind(:)
   !real(dp), intent(in) :: sTimesGrndEigVecs(:,:,:), grndEigVecs(:,:,:)
    real(dp), intent(in) :: occNr(:,:)
    real(dp), intent(in) :: gamma(:,:), lrGamma(:,:)
    integer,  intent(in) :: getij(:,:)     
    integer, intent(in) :: iatrans(1:,homo+1:)
    integer, intent(in) :: species0(:)
    real(dp), intent(in) :: uhbuu(:), uhbud(:)
   !type(TOrbitals), intent(in) :: orb    
    real(dp), intent(in) :: tqov(:,:), tqoo(:,:), tqvv(:,:)
    real(dp), intent(out) :: vp(:,:), vm(:,:)
    real(dp), intent(out) :: Mp(:,:), Mm(:,:) 

    integer :: izpalpha, nmat, natom, nw
    integer :: ia, jb, alpha, ii, jj, aa, bb, aa_old, bb_old, jj_old, ij, ab
    real(dp) :: tmp2
    real(dp), allocatable :: otmp(:), otmp2(:), qij(:), wnij(:)
    real(dp) :: dtmp, dtmp2    

    nmat = size(vp, dim=1) ! also known as nxov
    natom = size(gamma, dim=1)
    nw = size(wij)
    
    allocate(otmp(natom))
    allocate(otmp2(natom))
    allocate(qij(natom))
    allocate(wnij(nw))
   
    call wtdn(wij, occNr, win, nMatUp, nmat, getij, wnij)    

    !calc. vp, vm
    
    vp(:,:) = 0.0_dp 
    vm(:,:) = 0.0_dp
    Mp(:,:) = 0.0_dp
    Mm(:,:) = 0.0_dp
    
    otmp(:) = 0.0_dp
    
    !----only spin unpolarized case for now-------------------------------

    if (sym == 'S') then
       !full range coupling matrix contribution
       do jb = 1, initDim !calc 4* sum_A q^ia_A sum_B gamma_AB q^jb_B
          call hemv(otmp, gamma, tqov(:,jb))
          do ia = 1, nmat
             vp(ia,jb) = 4.0 * dot_product(tqov(:,ia), otmp) 
          end do
       end do
    else
       do jb = 1, initDim
          otmp(:) = otmp + tqov(:,jb)
          do ia = 1, nmat
             vp(ia,jb) = 0.0_dp
             do alpha = 1, natom
                izpalpha = species0(alpha)
                tmp2 = 2.0_dp * (uhbuu(izpalpha) - uhbud(izpalpha)) !== 2*magnetization
                vp(ia,jb) = vp(ia,jb) + otmp(alpha) * tmp2 * tqov(alpha,jb) 
             end do
          end do
       end do
    end if

    aa_old = -1
    bb_old = -1
    jj_old = -1

    !long range coupling matrix contribution
    do jb = 1, initDim
       call indxov(win, jb, getij, jj, bb)
       do ia = 1, nmat
          call indxov(win, ia, getij, ii, aa)
          if ((aa_old .ne. aa) .or. (bb_old .ne. bb)) then
             if (bb < aa) then
                call rindxvv(homo, aa, bb, ab)
             else
                call rindxvv(homo, bb, aa, ab)
             end if
             call hemv(otmp, lrGamma, tqvv(:,ab))
          end if

          if (jj < ii) then
             ij = jj + ((ii - 1 - homo + nocc) * (ii - homo + nocc)) / 2 + homo - nocc
          else
             ij = ii + ((jj - 1 - homo + nocc) * (jj - homo + nocc)) / 2 + homo - nocc
          end if

          dtmp = dot_product(tqoo(:,ij), otmp)

          vp(ia,jb) = vp(ia,jb) - dtmp
          vm(ia,jb) = vm(ia,jb) - dtmp

         !if (iatrans(jj,aa) <= nmat .and. iatrans(ii,bb) <= nmat) then ! TOMAS KUBAR ADDED

             if ((jj_old .ne. jj) .or. (aa_old .ne. aa)) then
                call hemv(otmp2, lrGamma, tqov(:,iatrans(jj,aa)))
             end if
             
             dtmp2 = dot_product(tqov(:,iatrans(ii,bb)), otmp2)
             
             vp(ia,jb) = vp(ia,jb) - dtmp2
             vm(ia,jb) = vm(ia,jb) + dtmp2

         !end if

          aa_old = aa
       end do
       jj_old = jj
       bb_old = bb
    end do

    do jb = 1, initDim
       vp(jb,jb) = vp(jb,jb) + 0.5_dp * wnij(jb)  !orb. energy difference diagonal contribution
       vm(jb,jb) = vm(jb,jb) + 0.5_dp * wnij(jb)
    end do
    !end calc. vp, vm

    !calc. vm

    !calc. matrices 
    do ii = 1, initDim
       do jj = ii, initDim
          Mp(ii,jj) = vp(ii,jj)
          Mp(jj,ii) = Mp(ii,jj)
          Mm(ii,jj) = vm(ii,jj)
          Mm(jj,ii) = Mm(ii,jj)
       end do
    end do

!    write(*,*), 'debug vm ', vm(:,1) ! seems to come out the same

    deallocate(otmp)
    deallocate(otmp2)
    deallocate(qij)
    deallocate(wnij)

  end subroutine setupInitMatFast


!>Run Linear Response calc. The actual diagonalization of the RPA equation is outsourced to Rs_LinRespCalc. Most code copied from non range-separated version by A. Dominguez.
  subroutine runRs_LinRespCalc(spin, tOnsite, natom, iAtomStart, grndEigVecs, grndEigVal, sccCalc, dQ, coord0,&
      & nexc, nstat0, symc, SSqr, filling, species0, nbeweg, uhubb, uhbuu, uhbud, rnel,&
      & iNeighbor, img2CentCell, orb, rs_data, tWriteTagged, fdTagged, taggedWriter,&
     !& tMulliken, tCoeffs, tXplusY, tTrans, tTradip, tArnoldi,&
      & fdMulliken, fdCoeffs, fdXplusY, fdTrans, fdSPTrans, fdTradip, fdTransQ, tArnoldi, fdArnoldi, fdExc,& !ons_en, ons_dip,
      & tEnergyWindow, energyWindow, tOscillatorWindow, oscillatorWindow, tCacheCharges,&
      & omega, shift, skHamCont,&
      & skOverCont, derivator, deltaRho, excgrad, dQAtomEx)
    implicit none
    logical, intent(in) :: spin
    logical, intent(inout) :: tOnsite !intent out so it can be forced to .false. for now
    integer, intent(in) :: natom, iAtomStart(:)
    real(dp), intent(in) :: grndEigVecs(:,:,:), grndEigVal(:,:), dQ(:), coord0(:,:)
    !> Self-consistent charge module settings
    type(TScc), intent(in) :: sccCalc
    integer, intent(in) :: nexc, nstat0
    character, intent(in) :: symc
    real(dp), intent(in) :: SSqr(:,:)
    real(dp), intent(in) :: filling(:,:)
    integer, intent(in) :: species0(:), nbeweg
    real(dp), intent(in) :: uhubb(:), uhbuu(:), uhbud(:), rnel !maybe needed for gradients later on
    integer, intent(in) :: iNeighbor(0:,:), img2CentCell(:)
    type(TOrbitals), intent(in) :: orb
    type(RangeSepFunc), intent(inout) :: rs_data !contains long-range gamma
    !> print tag information
    logical, intent(in) :: tWriteTagged

    !> file id for tagging information
    integer, intent(in) :: fdTagged

    !> tagged writer
    type(TTaggedWriter), intent(inout) :: taggedWriter

   !logical :: tMulliken, tCoeffs, tXplusY, tTrans, tTradip, tArnoldi

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
    integer, intent(in) :: fdTradip
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
    type(OSlakoCont), intent(in), optional :: skHamCont, skOverCont
    !> Differentiatior for the non-scc matrices
    class(NonSccDiff), intent(in), optional :: derivator
    real(dp), intent(inout), optional :: excgrad(:,:)   !also may be needed for gradients later on
    real(dp), intent(inout), optional :: dQAtomEx(:)

    real(dp) :: Ssq(nexc)
    real(dp), allocatable :: gamma(:,:), lrGamma(:,:), snglPartTransDip(:,:), sTimesGrndEigVecs(:,:,:), wij(:)
    real(dp), allocatable :: snglPartOscStrength(:), oscStrength(:), xmy(:), xpy(:), pc(:,:)
    real(dp), allocatable :: t(:,:), rhs(:), woo(:), wvv(:), wov(:)
    real(dp), allocatable :: evec(:,:), eval(:), xpyAll(:,:), transitionDipoles(:,:), atomicTransQ(:)
    real(dp), allocatable :: tqov(:,:), tqoo(:,:), tqvv(:,:)!to hold precalculated transition charges, alloc and calc in rs calc
    integer, allocatable :: win(:), iatrans(:,:), getij(:,:)
    character, allocatable :: symmetries(:)

    integer :: nstat, nocc, nocc_r, nvir_r, nxoo_r, nxvv_r
    integer :: nxov, nxov_ud(2), nxov_r, nxov_d, nxov_rd, norb
    integer :: i, j, isym
    integer :: spin_dim
    character :: sym
    logical :: tGrads, tZVector
    real(dp) :: energyThreshold
    integer :: logfil
   !real(dp) ::  t0,t1,t2,t3,t4,t5

    !> spin-polarized calculation?
    logical :: tSpin
    integer :: nSpin, iSpin
    !> control variables
    logical :: tCoeffs, tTradip, tTransQ, tTrans, tXplusY !, tZVector
    !> printing data
    logical :: tMulliken

    !> transition charges, either cached or evaluated on demand
    type(TTransCharges) :: transChrg

    !> aux variables
    integer :: mu, nu
    
    tSpin = .false. !For now, spin-polarized calculation not supported
    nSpin = 1
    tOnsite = .false. !For now, ons not supported with range separation
   !logfil = 46

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

    spin_dim = size(grndEigVal, dim=2)
    norb = orb%nOrb

    @:ASSERT(present(excgrad) .eqv. present(shift))
    @:ASSERT(present(shift) .eqv. present(skHamCont))
    @:ASSERT(present(skHamCont) .eqv. present(skOverCont))

    ! work out which data files are required, based on whether they have valid file IDs (>0)
    tMulliken = (fdMulliken > 0)
    tCoeffs = (fdCoeffs > 0)
    tTradip = (fdTradip > 0)
    tTransQ = (fdTransQ > 0)
    tTrans = (fdTrans > 0)
    tXplusY = (fdXplusY > 0)

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
    
    ! Transition charges
    if (tTransQ) then
       allocate(atomicTransQ(natom))
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

    if (nexc + 1 >= nxov) then
      write(tmpStr,"(' Insufficent single particle excitations, ',I0,&
          & ', for required number of excited states ',I0)")nxov, nexc
      call error(tmpStr)
    end if

    ! are gradients required?
    tGrads = present(excgrad)  ! For now no gradients, but keep for later!
    ! is a Z-vector required?
    tZVector = tGrads .or. tMulliken .or. tCoeffs

    ! Sanity checks
    nstat = nstat0
    if (nstat < 0 .and. symc /= "S") then
      call error("Linresp: Brightest mode only for singlets.")
    elseif (nstat /= 0 .and. symc == "B") then
      call error("Linresp: Both symmetries not allowed if specific state is&
          & excited")
    elseif (nstat == 0 .and. (tGrads .or. tMulliken .or. tCoeffs)) then
      call error("Linresp: Gradient, charges and coefficients only with&
          & selected excitation.")
    elseif (tGrads .and. nexc > nxov) then
      call error("Linresp: With gradients nexc can be max. the number of&
          & occupied-virtual excitations")
    end if

    ! Select symmetries to process
    !ADG: temporal solution for spin case.
    if (.not. spin) then
      select case (symc)
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
      allocate(symmetries(1))
      symmetries(:) = [ " " ]
    end if

    ! Allocation for general arrays
    allocate(gamma(natom, natom))
    allocate(lrGamma(natom, natom))
    allocate(snglPartTransDip(nxov, 3))
    !allocate(sTimesGrndEigVecs(norb, norb))
    allocate(sTimesGrndEigVecs(norb, norb, spin_dim))
    allocate(wij(nxov))
    allocate(win(nxov))
!    allocate(evec(nxov, nexc))
    allocate(eval(nexc))
    allocate(getij(nxov, 2))
    allocate(transitionDipoles(nxov, 3))
    allocate(snglPartOscStrength(nxov))
    
    ! Arrays for gradients and Mulliken analysis !Again not required for now
    if (tGrads .or. tMulliken) then
      allocate(pc(norb, norb))
    end if

    ! Overlap times wave function coefficients - most routines in DFTB+ use lower triangle (would
    ! remove the need to symmetrize the overlap and ground state density matrix in the main code if
    ! this could be used everywhere in these routines)
    do iSpin = 1, nSpin
      call symm(sTimesGrndEigVecs(:,:,iSpin), "L", SSqr, grndEigVecs(:,:,iSpin))
    end do

    ! TEST OK -- sTimesGrndEigVecs

    ! ground state Hubbard U softened coulombic interactions
    call sccCalc%getAtomicGammaMatrix(gamma, iNeighbor, img2CentCell)
    ! gamma is symmetric but still the upper triangle is not enough
    do mu = 2, natom
       do nu = 1, mu-1
          gamma(nu, mu) = gamma(mu, nu)
       enddo
    enddo
    call rs_data%getLrGammaEval(lrGamma)

    ! TEST OK -- gamma, lrGamma (gamma was triangular only, previously)

    ! Oscillator strengths for excited states, when needed.
   !if (nstat <= 0) then
      allocate(oscStrength(nexc))
   !end if

    ! Find all single particle transitions and KS energy differences
    !   for cases that go from filled to empty states
    call getSPExcitations(grndEigVal, filling, wij, getij)

    ! put them in ascending energy order
    if (tOscillatorWindow) then
      ! use a stable sort so that degenerate transitions from the same single particle state are
      ! grouped together in the results, allowing these to be selected together (since how intensity
      ! is shared out over degenerate transitions is arbitrary between eigensolvers/platforms).
      call merge_sort(win, wij, 1.0_dp*epsilon(1.0))
    else
      ! do not require stability, use the usual routine to sort, saving an O(N) workspace
      call index_heap_sort(win, wij)
    end if
    wij = wij(win)

    ! TEST OK -- wij, getij

    ! dipole strength of transitions between K-S states
    call calcTransitionDipoles(coord0, win, nxov_ud(1), getij, iAtomStart, sTimesGrndEigVecs, grndEigVecs, snglPartTransDip)

    ! single particle excitation oscillator strengths
    snglPartOscStrength(:) = twothird * wij(:) * sum(snglPartTransDip**2, dim=2)

    ! TEST OK -- snglPartOscStrength

    if (tOscillatorWindow .and. tZVector ) then
      call error("Incompabilitity between excited state property evaluation and an oscillator&
          & strength window at the moment.")
    end if

    if (tOscillatorWindow .or. tEnergyWindow) then

      if (.not. tEnergyWindow) then

        ! find transitions that are strongly dipole allowed (> oscillatorWindow)
        call dipselect(wij, snglPartOscStrength, win, snglPartTransDip, nxov_rd, oscillatorWindow, grndEigVal, getij)

      else

        ! energy window above the lowest nexc single particle transitions
        energyThreshold = wij(nexc) + energyWindow
        nxov_r = count(wij <= energyThreshold)
        nxov_d = 0

        if (tOscillatorWindow) then

          ! find transitions that are strongly dipole allowed (> oscillatorWindow)
          if (nxov_r < nxov) then
            ! find transitions that are strongly dipole allowed (> oscillatorWindow)
            call dipselect(wij(nxov_r+1:), snglPartOscStrength(nxov_r+1:), win(nxov_r+1:),&
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
    nxov_rd = max(nxov_rd, min(nexc+1, nxov))

    ! TEST OK -- nxov_rd

    call TTransCharges_init(transChrg, iAtomStart, sTimesGrndEigVecs, grndEigVecs, nxov_rd, nxov_ud(1), getij,&
        & win, tCacheCharges)

   !if (nstat == 0) then
   !   if(tTrans) then
   !      call writeSPExcitations(wij, snglPartTransDip, win, nxov_ud(1), getij, tWriteTagged, fdTagged, taggedWriter, snglPartOscStrength)
   !   endif
   !   call openExcitationFiles(spin, tXplusY, tTrans, tTradip, tArnoldi, logfil)
   !end if

    if (tXplusY) then
      open(fdXplusY, file=XplusYOut, position="rewind", status="replace")
    end if

    if(tTrans) then
      open(fdTrans, file=transitionsOut, position="rewind", status="replace")
      write(fdTrans,*)
    endif

    ! single particle transition dipole file
    if (tTradip) then
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
    call writeSPExcitations(wij, win, nxov_ud(1), getij, fdSPTrans, snglPartOscStrength, nxov_rd, tSpin)

    ! redefine if needed (generalize it for spin-polarized and fractional occupancy)
    nocc = nint(rnel) / 2
    nocc_r = nOcc
    nvir_r = nOrb - nOcc

    ! elements in a triangle plus the diagonal of the occ-occ and virt-virt blocks
    nxoo_r = (nocc_r * (nocc_r + 1)) / 2
    nxvv_r = (nvir_r * (nvir_r + 1)) / 2

    allocate(evec(nxov_rd, nexc))

   !allocate(xpyAll(nxov_ud(1), nexc)) ! allocate(xpyAll(nxov_rd, nexc))
    allocate(xpyAll(nxov_rd, nexc))

    if (nstat > 0) then
      !allocate(xmy(nxov_ud(1)))
      !allocate(xpy(nxov_ud(1)))
       allocate(xmy(nxov_rd))
       allocate(xpy(nxov_rd))
    end if

    ! Arrays needed for Z vector
    if (tZVector) then
      allocate(t(norb, norb))
      allocate(rhs(nxov_rd))
      allocate(woo(nxoo_r))
      allocate(wvv(nxvv_r))
      allocate(wov(nxov_rd))
    end if
    
    !for fast init. mat calc. need combined index
    allocate(iatrans(1:nocc, nocc+1:norb))
    call rindxov_array(win, nocc, nxov, getij, iatrans)

    !run lin. resp. calculation:
    do isym = 1, size(symmetries)
      sym = symmetries(isym)
      !call buildAndDiagExcMatrix(spin, tOnsite, wij, sym, win, nxov_ud(1), nxov_rd, iAtomStart,&
      !    & sTimesGrndEigVecs, grndEigVecs, filling, getij, gamma,&
      !    & species0, uhbuu, uhbud, ons_en, orb, eval, evec)

     !call Rs_LinRespCalc(spin, tZVector, wij, nexc, sym, win, nxov_ud(1), nxov_rd, nocc, nocc_r, nvir_r, iAtomStart,&
     !      & sTimesGrndEigVecs, grndEigVecs, filling, getij, iatrans, gamma, lrGamma,&
     !      & species0, uhbuu, uhbud, orb, eval, evec, xpyAll, nstat, xmy, tqov, tqoo, tqvv)
      call Rs_LinRespCalc(tZVector, tTransQ, wij, nexc, sym, win, nxov_ud(1), nxov_rd, nocc, nocc_r, nvir_r, iAtomStart,&
            & sTimesGrndEigVecs, grndEigVecs, filling, getij, iatrans, gamma, lrGamma,&
            & species0, uhbuu, uhbud, eval, evec, xpyAll, nstat, xmy, tqov, tqoo, tqvv)

    ! TEST OK -- eval, evec
    ! TEST OK -- tqov, tqoo, tqvv

     !do i=1,nexc
     !write (*,'(I3,3X,F7.3)') i, eval(i)
     !end do

      write (*,*) "1st eigenstate and X+Y with energy", eval(1), "is:"
       do i=1,size(evec,dim=1)
          write (*,'(I5,F9.5)') i, evec(i,1)
          write (*,'(I5,A,F9.5)') i, "          ", xpyAll(i,1)
       end do

     !do i=1,nxov_rd
     !write (*,'(I3,3X,3F7.3)') i, snglPartTransDip(i,:)
     !end do

      write (*,*) "xpyAll"
      do i=1,nexc
      write (*,'(I3,3X,200F7.3)') i, xpyAll(:,i)
      end do

     !if (nstat <= 0) then
        call getOscillatorStrengths_rs(sym, snglPartTransDip(1:nxov_rd,:), eval, xpyAll, filling,&
            & nstat, oscStrength, tTradip, transitionDipoles) 
      write (*,*) "osc strengths"
      do i=1,nexc
      write (*,'(I3,3X,F7.3)') i, oscStrength(i)
      end do
       !if (nstat == 0) then
          if (spin) then
            call getExcSpin(Ssq, nxov_ud(1), getij, win, eval, evec, wij,&
                & filling, sTimesGrndEigVecs, grndEigVecs)
            call writeExcitations_rs(sym, oscStrength, nexc, nxov_ud(1), getij, win, eval,& ! evec,&
                & xpyAll, wij(1:nxov_rd), fdXplusY, fdTrans, fdTradip, transitionDipoles,&
                & tWriteTagged, fdTagged, taggedWriter, fdExc, Ssq)
          else
            call writeExcitations_rs(sym, oscStrength, nexc, nxov_ud(1), getij, win, eval,& ! evec,&
                & xpyAll, wij(1:nxov_rd), fdXplusY, fdTrans, fdTradip, transitionDipoles,&
                & tWriteTagged, fdTagged, taggedWriter, fdExc)
          end if
       !end if
     !end if

      if (tTransQ) then
         call transitionCharges_rs(xpyAll(:,nstat), tqov, atomicTransQ)
         if (.not. tZVector) then
            deallocate(tqov)
         end if
         call writeTransitionCharges_rs(fdTransQ, atomicTransQ)
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
     !call closeExcitationFiles(tXplusY, tTrans, tTradip, tArnoldi, logfil)
      omega = 0.0_dp
      return
    end if

    omega = sqrt(eval(nstat))  !Attention: right now I take sqrt in Rs_LinRespCalc. May not want to do this!!!!!!!!!
    
    ! TEST OK -- omega

    if (tZVector) then

      ! Furche terms: X+Y, X-Y
     !xmy(:) = sqrt(omega) / sqrt(wij) * evec(:,nstat)    ! xmy was set in Rs_LinRespCalc
      xpy(:) = xpyAll(:,nstat)

      !call rindxov_array(win, nocc, nxov, getij, iatrans) ! already performed above (for efficient init mat calc)

      !!call cpu_time(t0)
      call getZVectorEqRHS_rs(xpy, xmy, win, iAtomStart, nocc, nocc_r, nxov_ud(1), getij, iatrans, natom, species0,&
           & grndEigVal(:,1), sTimesGrndEigVecs, grndEigVecs, gamma, lrGamma, uhbuu, uhbud, omega, sym,&
           & rhs, t, wov, woo, wvv)

    ! TEST OK -- RHS, T, WOV, WOO, WVV

      !!call cpu_time(t1)
      call solveZVectorEq_rs(rhs, win, nxov_ud(1), getij, filling, natom, & ! iAtomStart, sTimesGrndEigVecs,
           & species0, uhbuu, uhbud, gamma, lrGamma, wij, & ! grndEigVecs, orb,
           & nocc, nocc_r, nvir_r, iatrans, tqov, tqoo, tqvv)

    ! TEST OK -- RHS SOLVED

      !!call cpu_time(t2)
      call calcWVectorZ_rs(rhs, win, nocc, norb, nxov_ud(1), getij, iAtomStart, sTimesGrndEigVecs, grndEigVecs,&
           & gamma, lrGamma, grndEigVal(:,1), wov, woo, wvv, iatrans)

    ! TEST OK -- WOV, WOO, WVV MODIF

      !!call cpu_time(t3)
      call calcPMatrix(t, rhs, win, getij, pc)

    ! TEST OK -- PC

      ! if (tCoeffs) then
      ! !this has to be checked in case of orbital constraints 
      !   call writeCoeffs(pc, grndEigVecs, filling, nocc)
      ! end if
      call transformMO2AODense(pc, grndEigVecs(:,:,1))

    ! TEST OK -- PC TRANSFORMED

      call getExcMulliken(iAtomStart, pc, SSqr, dQAtomEx)
      if (tMulliken) then
        call writeExcMulliken(sym, nstat, dQ, dQAtomEx, coord0, fdMulliken)
      end if

      ! TEST OK -- DQEX

      if (tGrads) then
      !call cpu_time(t4)
         call addGradients_rs(sym, nxov_rd, natom, species0, iAtomStart, norb, nocc, nxov_ud(1), getij, win,& !iatrans,
              & grndEigVecs, pc, dQ, dQAtomEx, gamma, lrGamma, rs_data, uhubb, uhbuu, uhbud,&
              & shift, woo, wov, wvv, xpy, xmy, nbeweg, coord0, orb,&
              & skHamCont, skOverCont, tqov, derivator, deltaRho, excgrad)

      ! TEST OK -- EXGRAD

      !!call cpu_time(t5)
      !!write(*,'(a,2x,f18.10)') 'Time getZVectorEqRHS', t1-t0
      !!write(*,'(a,2x,f18.10)') 'Time solveZVectorEq', t2-t1
      !!write(*,'(a,2x,f18.10)') 'Time calcWVectorZ_rs', t3-t2
      !!write(*,'(a,2x,f18.10)') 'Time addGradients', t5-t4
      end if
   end if
  end subroutine runRs_LinRespCalc


!>Perform linear response calculation with algorithm by Stratmann and Scuseria (JCP 1998)
  subroutine Rs_LinRespCalc( & !spin,
       & tZVector, tTransQ, wij, nexc, sym, win, nMatUp, nxov, homo, nocc, nvir, iAtomStart,&
       & sTimesGrndEigVecs, grndEigVecs, occNr, getij, iatrans, gamma, lrGamma, species0, uhbuu, uhbud,& ! orb,
       & eval, evec, xpy, nstat, xmy, tqov, tqoo, tqvv)
    implicit none
    logical, intent(in) :: tZVector, tTransQ !, spin
    real(dp),intent(in) :: wij(:)
    character, intent(in) :: sym
    integer, intent(in) :: nexc !desired number of excitations
    integer, intent(in) :: win(:), nMatUp, nxov, homo, nocc, nvir, iAtomStart(:)
    real(dp), intent(in) :: sTimesGrndEigVecs(:,:,:), grndEigVecs(:,:,:)  
    real(dp), intent(in) :: occNr(:,:)
    real(dp), intent(in) :: gamma(:,:), lrGamma(:,:)
    integer, intent(in) :: getij(:,:), species0(:)
    integer, intent(in) :: iatrans(1:,homo+1:) !needed in fast setup of init mat
    real(dp), intent(in) :: uhbuu(:), uhbud(:)
   !type(TOrbitals), intent(in) :: orb
    integer, intent(in) :: nstat
    real(dp), intent(out) :: eval(:), evec(:,:), xpy(:,:) !lin. resp. eigenvalues and eigenvectors (F_I and X_I+Y_I)
    real(dp), intent(out) :: xmy(:)
    real(dp), allocatable, intent(out) :: tqov(:,:), tqoo(:,:), tqvv(:,:)!hold precalculated transition charges

    integer :: ii, jj, ia, ij, aa, bb, ab, nxov_act
    integer :: subSpaceDim, prevSubSpaceDim, memDim, workDim !cur. size of subspace and memory for subspace vectors and matrices
    real(dp), allocatable :: bvec(:,:) !basis of subspace
    real(dp), allocatable :: Lvec(:,:), Rvec(:,:) !left and right eigenvectors of Mnh
    real(dp), allocatable :: vp(:,:), vm(:,:) !vec. for (A+B)b_i, (A-B)b_i
    real(dp), allocatable :: Mp(:,:), Mm(:,:), Mmsqrt(:,:), Mmsqrtinv(:,:),  Mh(:,:), xpy_test(:,:)
           !Matrices M_plus, M_minus, M_minus^(1/2), M_minus^(-1/2) and M_herm~=resp. mat on subsapce
    real(dp), allocatable :: evalInt(:) !store eigenvectors within routine
    real(dp), allocatable :: dummyM(:,:), workArray(:)
    real(dp), allocatable :: vecNorm(:) !will hold norms of residual vectors
    real(dp), allocatable :: gqvtmp(:,:) !work array used in matrix multiplications
    real(dp), allocatable :: qij(:)
    real(dp) :: dummyReal
    logical :: lUpdwn
    integer :: info
    integer :: dummyInt
    logical :: didConverge
    real(dp), parameter :: convThreshold = 0.0001_dp 
    integer :: newVec
    integer :: initExc 
    initExc = min(20 * nexc, nocc * nvir)

    nxov_act = nxov ! or nocc * nvir
    allocate(xpy_test(size(evec,dim=1),size(evec,dim=2)))

    !initial allocations 
    subSpaceDim = min(initExc, nxov) !start with lowest excitations. Inital number somewhat arbritary.
    memDim = min(subSpaceDim + 6 * nexc, nxov) !memory available for subspace calcs.
    workDim = 3 * memDim + 1
    allocate(bvec(nxov_act, memDim)) ! nxov, memDim))
    allocate(vp(nxov_act, memDim)) ! nxov, memDim))
    allocate(vm(nxov_act, memDim)) ! nxov, memDim))
    allocate(Mp(memDim, memDim))
    allocate(Mm(memDim, memDim))
    allocate(Mmsqrt(memDim, memDim))
    allocate(Mmsqrtinv(memDim, memDim))
    allocate(Mh(memDim, memDim))
    allocate(dummyM(memDim, memDim))
    allocate(evalInt(memDim))
    allocate(Lvec(memDim, nexc))
    allocate(Rvec(memDim, nexc))
    allocate(workArray(3 * memDim + 1))
    allocate(vecNorm(2 * memDim))

 !!!Work array for faster implementation of (A+-B)*v. Make optional later!
    allocate(gqvtmp(size(gamma, dim=1), max(nvir, nocc) * nvir)) !could use symmetry to improve nvir * nvir case
    allocate(tqov(size(gamma, dim=1), nocc * nvir)) ! nxov -- TOMAS KUBAR CHANGED
    allocate(tqoo(size(gamma, dim=1), nocc * (nocc + 1) / 2))
    allocate(tqvv(size(gamma, dim=1), nvir * (nvir + 1) / 2))
    allocate(qij(size(gamma, dim=1)))
 !   allocate(iatrans, (1:homo,homo+1:homo + nvir)) !maybe need nvir_unrestricted
  !!precalculate transition charges
    do ia = 1, nvir * nocc ! nxov ! nMatUp -- TOMAS KUBAR CHANGED
       call indxov(win, ia, getij, ii, aa)
       lUpdwn = (win(ia) <= nMatUp)
       qij = transq(ii, aa, iAtomStart, lUpdwn, sTimesGrndEigVecs, grndEigVecs)
       tqov(:,ia) = qij 
    end do
    do ij = 1, nocc * (nocc + 1) / 2 
       call indxoo(ij, ii, jj)
       lUpdwn = .true. ! UNTESTED
       qij = transq(ii, jj, iAtomStart, lUpdwn, sTimesGrndEigVecs, grndEigVecs)
       tqoo(:,ij) = qij
    end do
    do ab = 1, nvir * (nvir + 1) / 2
       call indxvv(homo, ab, aa, bb)
       lUpdwn = .true. ! UNTESTED
       qij = transq(aa, bb, iAtomStart, lUpdwn, sTimesGrndEigVecs, grndEigVecs)
       tqvv(:,ab) = qij
    end do
!!!calc. array with ia indices 
!!   call rindxov_array(win, nMatUp, homo, nocc * nvir, getij, iatrans)
!    call rindxov_array(win, homo, nocc * nvir, getij, iatrans)

!    write(*,*) ' debug gamma  ', gamma(:,:)

    !set initial bs
    bvec(:,:) = 0.0_dp
    do ii = 1, subSpaceDim
       bvec(ii, ii) = 1.0
    end do
   
    !Could calc. init. matrix here, Due to special form of bs, matrix multiplication handles a lot of zeros.
    !Add this latter!

    prevSubSpaceDim = 0
    didConverge = .false.

    !solution of linear response problem. Iterative expansion of subspace:
    solveLinResp: do

       if (prevSubSpaceDim > 0) then

          !extend subspace matrices:
          do ii = prevSubSpaceDim + 1, subSpaceDim
            !call multApBVecFast(spin, bvec(:,ii), wij, sym, win, nMatUp, homo, nocc, nvir, iAtomStart, &
            !    & sTimesGrndEigVecs, grndEigVecs, occNr, getij, gamma, lrGamma, species0, uhbuu, uhbud, &
            !    & orb, iatrans, gqvtmp, tqov, tqoo, tqvv, vp(:,ii))
            !call multAmBVecFast(spin, bvec(:,ii), wij, sym, win, nMatUp, homo, nocc, nvir, iAtomStart, &
            !    & sTimesGrndEigVecs, grndEigVecs, occNr, getij, gamma, lrGamma, species0, &
            !    & orb, iatrans, gqvtmp, tqov, tqoo, tqvv, vm(:,ii))
             call multApBVecFast(bvec(:,ii), wij, sym, win, nMatUp, homo, nocc, nvir, &
                 & occNr, getij, gamma, lrGamma, species0, uhbuu, uhbud, iatrans, gqvtmp, tqov, tqoo, tqvv, vp(:,ii))
             call multAmBVecFast(bvec(:,ii), wij, win,      nMatUp, homo, nocc, nvir, &
                 & occNr, getij, gamma, lrGamma,                         iatrans, gqvtmp, tqov, tqoo, tqvv, vm(:,ii))
          end do

          do ii = prevSubSpaceDim + 1, subSpaceDim
             do jj = 1, ii
                Mp(ii,jj) = dot_product(bvec(:,jj), vp(:,ii))
                Mp(jj,ii) = Mp(ii,jj)
                Mm(ii,jj) = dot_product(bvec(:,jj), vm(:,ii))
                Mm(jj,ii) = Mm(ii,jj)
             end do
          end do

       else
         
          call setupInitMatFast(subSpaceDim, wij, sym, win, nMatUp, nocc, homo, &
              & occNr, getij, iatrans, gamma, lrGamma, species0, uhbuu, uhbud, tqov, tqoo, tqvv, vp, vm, Mp, Mm)
          !call setupInitMat(initExc, spin, wij, sym, win, nMatUp, iAtomStart, sTimesGrndEigVecs, &
          !    &grndEigVecs, occNr, getij, gamma, lrGamma, species0, uhbuu, uhbud, orb, vp, vm, Mp, Mm)
!          write(*,*), 'debug after init ', Mm 
       end if
!      write(*,*) 'debug mp ', mp(1:10,1:10)
      
      call calcMatrixSqrt(Mm, subSpaceDim, memDim, workArray, workDim,  Mmsqrt, Mmsqrtinv)
      call dsymm('L', 'U', subSpaceDim, subSpaceDim, 1.0_dp, Mp, memDim, Mmsqrt, memDim, 0.0_dp, dummyM, memDim)
      call dsymm('L', 'U', subSpaceDim, subSpaceDim, 1.0_dp, Mmsqrt, memDim, dummyM, memDim, 0.0_dp, Mh, memDim) 
!      write(*,*) 'debug msqrt ', mmsqrt(1:10,1:10)     

!      write(*,*) ,' debug Mh: ', Mh(1:subspacedim,1:subspacedim)
      !diagonalize in subspace
      call dsyev('V', 'U', subSpaceDim, Mh, memDim, evalInt, workArray, workDim, info)
     
      !This yields T=(A-B)^(-1/2)|X+Y>. Calc. |R_n>=|X+Y>=(A-B)^(1/2)T |L_n>=|X-Y>=(A-B)^(-1/2)T. Transformation preserves orthonormality. 
      !Only compute up to nexc index, because only that much needed.
      call dsymm('L', 'U', subSpaceDim, nexc, 1.0_dp, Mmsqrt,    memDim, Mh, memDim, 0.0_dp, Rvec, memDim) 
      call dsymm('L', 'U', subSpaceDim, nexc, 1.0_dp, Mmsqrtinv, memDim, Mh, memDim, 0.0_dp, Lvec, memDim)
        !check if dsymm can actually be used

      !Need |X-Y>=sqrt(w)(A-B)^(-1/2)T, |X+Y>=(A-B)^(1/2)T/sqrt(w) for proper solution to original EV problem
      do ii = 1, nexc !only use first nexc vectors
         dummyReal = sqrt(sqrt(evalInt(ii)))
         Rvec(:,ii) = Rvec(:,ii) / dummyReal
         Lvec(:,ii) = dummyReal * Lvec(:,ii) 
      end do
     
      !see if more memory is needed to save extended basis. If so increase amount of memory.
      if (subSpaceDim + 2 * nexc > memDim) then
         call incSize(memDim, bvec)
         call incSize(memDim, vp)
         call incSize(memDim, vm)
         call incSizeMatBoth(memDim, Mp)
         call incSizeMatBoth(memDim, Mm)
         call incSizeMatBoth(memDim, Mh)
         call incSizeMatBoth(memDim, Mmsqrt)
         call incSizeMatBoth(memDim, Mmsqrtinv)
         call incSizeMatBoth(memDim, dummyM)
         call incSize(memDim, evalInt)
         call incSize(workDim, workArray)
         call incSizeMatSwapped(memDim, Lvec)
         call incSizeMatSwapped(memDim, Rvec)
         call incSize(2 * memDim, vecNorm)
         call incSize(memDim)
         call incSize(workDim)
      end if

      !Calculate the residual vectors
        !calcs. all |R_n>
     !call dgemm('N', 'N', nxov, nexc, subSpaceDim, 1.0_dp, bvec, nxov, Rvec, memDim, 0.0_dp, bvec(1,subSpaceDim+1),      nxov) 
      call dgemm('N', 'N', nxov_act, nexc, subSpaceDim, 1.0_dp, bvec, nxov_act, &
                & Rvec, memDim, 0.0_dp, bvec(1,subSpaceDim+1),      nxov_act) 
        !calcs. all |L_n>
     !call dgemm('N', 'N', nxov, nexc, subSpaceDim, 1.0_dp, bvec, nxov, Lvec, memDim, 0.0_dp, bvec(1,subSpaceDim+1+nexc), nxov) 
      call dgemm('N', 'N', nxov_act, nexc, subSpaceDim, 1.0_dp, bvec, nxov_act, &
                & Lvec, memDim, 0.0_dp, bvec(1,subSpaceDim+1+nexc), nxov_act) 
   
      do ii = 1, nexc
         dummyReal = -sqrt(evalInt(ii))
         bvec(:,subSpaceDim+ii)      = dummyReal * bvec(:,subSpaceDim+ii)
         bvec(:,subSpaceDim+nexc+ii) = dummyReal * bvec(:,subSpaceDim+nexc+ii)
      end do

      ! (A-B)|L_n> for all n=1,..,nexc 
     !call dgemm('N', 'N', nxov, nexc, subSpaceDim, 1.0_dp, vm, nxov, Lvec, memDim, 1.0_dp, bvec(1,subSpaceDim+1),      nxov)
      call dgemm('N', 'N', nxov_act, nexc, subSpaceDim, 1.0_dp, vm, nxov_act, &
                & Lvec, memDim, 1.0_dp, bvec(1,subSpaceDim+1),      nxov_act)
      ! (A+B)|R_n> for all n=1,..,nexc
     !call dgemm('N', 'N', nxov, nexc, subSpaceDim, 1.0_dp, vp, nxov, Rvec, memDim, 1.0_dp, bvec(1,subSpaceDim+1+nexc), nxov)
      call dgemm('N', 'N', nxov_act, nexc, subSpaceDim, 1.0_dp, vp, nxov_act, &
                & Rvec, memDim, 1.0_dp, bvec(1,subSpaceDim+1+nexc), nxov_act)
     
      !calc. norms of residual vectors to check for convergence
      didConverge = .true.
     
      do ii = subSpaceDim + 1, subSpaceDim + nexc
         vecNorm(ii-subSpaceDim) = dot_product(bvec(:,ii),bvec(:,ii))
         if (vecNorm(ii-subSpaceDim) .gt. convThreshold) then
            didConverge = .false.
         end if
      end do

      if (didConverge) then
         do ii = subSpaceDim + nexc + 1, subSpaceDim + 2 * nexc
            vecNorm(ii-subSpaceDim) = dot_product(bvec(:,ii),bvec(:,ii))
            if (vecNorm(ii-subSpaceDim) .gt. convThreshold) then
               didConverge = .false.
            end if
         end do
      end if

      if ((.not. didConverge) .and. (subSpaceDim + 2 * nexc > nxov)) then
         print *, 'Error: Linear Response calculation in subspace did not converge!'
         stop
      end if
      
      !if converged then exit loop:
      if (didConverge) then
         !do some wrapup, e.g.:
         write (*,'(A,I4)') 'Linear Response calculation converged. Required number of basis vectors: ', subSpaceDim 
         eval(:) = evalInt(1:nexc) 
      
        write (*,*) "evec", size(evec, dim=1), size(evec, dim=2)
        write (*,*) "xpy",  size(xpy,  dim=1), size(xpy,  dim=2)
        write (*,*) "bvec", size(bvec, dim=1), size(bvec, dim=2)
        write (*,*) "Mh",   size(Mh,   dim=1), size(Mh,   dim=2)
        write (*,*) "Rvec", size(Rvec, dim=1), size(Rvec, dim=2)
        write (*,*) "subSpaceDim", subSpaceDim
        !call dgemm('N', 'N', nxov, nexc, subSpaceDim, 1.0_dp, bvec, nxov, Mh,   memDim, 0.0_dp, evec, nxov) !Calc. F_I
        !call dgemm('N', 'N', nxov_act, nexc, subSpaceDim, 1.0_dp, bvec, nxov_act, Mh,   memDim, 0.0_dp, evec, nxov_act) !Calc. F_I
        !xpy_test = matmul(bvec(:,1:subSpaceDim), Mh(1:subSpaceDim,1:nexc))
         evec = matmul(bvec, Mh(:,1:nexc))
        !call dgemm('N', 'N', nxov_act, nexc, subSpaceDim, 1.0_dp, bvec, nxov_act, Rvec, memDim, 0.0_dp, xpy,  nvir*nocc) !Calc. X+Y
        !call dgemm('N', 'N', nxov, nexc, subSpaceDim, 1.0_dp, bvec, nxov, Rvec, memDim, 0.0_dp, xpy,  nxov) !Calc. X+Y
         xpy = matmul(bvec(:,1:subSpaceDim), Rvec(1:subSpaceDim,:))

         if (nstat > 0) then !Calc. X-Y for exc. nstat, will be used in force calculation
            xmy(:) = 0.0_dp
            do ii = 1, subSpaceDim
               xmy(:) = xmy + Lvec(ii,nstat) * bvec(:,ii)
            end do
         end if

         exit solveLinResp !terminate diag. routine
      end if
    
      !Otherwise calculate new basis vectors and extend subspace with them
      !only include new vectors if they add meaningful residue component
      newVec = 0
      do ii = 1, nexc
         if (vecNorm(ii) .gt. convThreshold) then
            newVec = newVec + 1
            dummyReal = sqrt(evalInt(ii))
            info = subSpaceDim + ii
            dummyInt = subSpaceDim + newVec
            do jj = 1, nxov
               bvec(jj,dummyInt) = bvec(jj,info) / (dummyReal - wij(jj))       
            end do
         end if
      end do

      do ii = 1, nexc
         if (vecNorm(nexc+ii) .gt. convThreshold) then
            newVec = newVec + 1
            info = subSpaceDim + nexc + ii
            dummyInt = subSpaceDim + newVec
            do jj = 1, nxov
               bvec(jj,dummyInt) = bvec(jj,info) / (dummyReal - wij(jj))
            end do
         end if
      end do

      prevSubSpaceDim = subSpaceDim
      subSpaceDim = subSpaceDim + newVec
      
      call orthonormalizeVectors(prevSubSpaceDim + 1, subSpaceDim, bvec) !create orthogonal basis
     
    end do solveLinResp

    if (.not. tZVector) then !not required anymore if gradient not computed
       if (.not. tTransQ) then
          deallocate(tqov)
       end if
       deallocate(tqoo)
       deallocate(tqvv)
    end if

  end subroutine Rs_LinRespCalc


  !> Write out transitions from ground to excited state along with single particle transitions and
  !> dipole strengths
  !> Modified routine because some expresions in original routine don't apply in the RS case
  subroutine writeExcitations_rs(sym, oscStrength, nexc, nMatUp, getij, win, eval, & ! evec, 
      & xplyAll, wij, fdXplusY, fdTrans, fdTradip, transitionDipoles, tWriteTagged, fdTagged, taggedWriter, fdExc, Ssq)
    implicit none

    !> Symmetry label for the type of transition
    character, intent(in) :: sym

    !> oscillator strengths for transitions from ground to excited states
    real(dp), intent(in) :: oscStrength(:)

    !> number of excited states to solve for
    integer, intent(in) :: nexc

    !> number of same spin excitations
    integer, intent(in) :: nMatUp

    !> index array between transitions in square and 1D representations
    integer, intent(in) :: getij(:,:)

    !> index array for single particle excitions
    integer, intent(in) :: win(:)

    !> excitation energies
    real(dp), intent(in) :: eval(:)

    !> eigenvectors of excited states
    !real(dp), intent(in) :: evec(:,:)

    !> X+Y, all of them
    real(dp), intent(in) :: xplyAll(:,:)

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

    !> tagged writer
    type(TTaggedWriter), intent(inout) :: taggedWriter

    !> file unit for excitation energies
    integer, intent(in) :: fdExc

    !> For spin polarized systems, measure of spin
    real(dp), intent(in), optional :: Ssq(:)

    integer :: nmat
    integer :: i, j, iweight, indo, m, n
    integer :: iDeg
    real(dp), allocatable :: wvec(:)
    real(dp), allocatable :: eDeg(:)
    real(dp), allocatable :: oDeg(:)
    integer, allocatable :: wvin(:)
    real(dp) :: weight, wvnorm !, rsqw
    logical :: lUpdwn, tSpin
    character :: sign

    @:ASSERT(fdExc > 0)

    tSpin = present(Ssq)
    nmat = size(wij)
    ALLOCATE(wvec(nmat))
    ALLOCATE(wvin(nmat))
    ALLOCATE(eDeg(nexc))
    ALLOCATE(oDeg(nexc))
    wvec = 0.0_dp
    wvin = 0
    eDeg = 0.0_dp
    oDeg = 0.0_dp

    if(fdXplusY > 0) then
      write(fdXplusY,*) nmat, nexc
    end if

    do i = 1, nexc
      if (eval(i) > 0.0_dp) then

        wvec(:) = xplyAll(:,i)**2
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
              & Hartree__eV * sqrt(eval(i)), oscStrength(i), m, '->', n, weight,&
              & Hartree__eV * wij(iweight), Ssq(i)
        else
          write(fdExc,&
              & '(1x,f10.3,4x,f14.8,5x,i5,3x,a,1x,i5,7x,f6.3,2x,f10.3,6x,a)')&
              & Hartree__eV * sqrt(eval(i)), oscStrength(i), m, '->', n, weight,&
              & Hartree__eV * wij(iweight), sign
        end if

        if(fdXplusY > 0) then
          if (tSpin) then
            lUpdwn = (win(iweight) <= nMatUp)
            sign = "D"
            if (lUpdwn) sign = "U"
          end if
          write(fdXplusY,'(1x,i5,3x,a,3x,ES17.10)') i, sign, sqrt(eval(i))
          write(fdXplusY,'(6(1x,ES17.10))') xplyAll(:,i)
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
              lUpdwn = (win(indo) <= nMatUp)
              sign = "D"
              if (lUpdwn) sign = "U"
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
              & '< 0', oscStrength(i), m, '->', n, '-', Hartree__eV * wij(iweight),&
              & Ssq(i)
        else
          write(fdExc,&
              & '(6x,A,T12,4x,f14.8,2x,i5,3x,a,1x,i5,7x,f6.3,2x,f10.3,6x,a)')&
              & '< 0', oscStrength(i), m, '->', n, weight,&
              & Hartree__eV * wij(iweight), sign
        end if

        if(fdXplusY > 0) then
          if (tSpin) then
            lUpdwn = (win(iweight) <= nMatUp)
            sign = "D"
            if (lUpdwn) sign = "U"
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
    oDeg(1) = oscStrength(1)
    do i = 2, nexc
      if(abs(eval(i)-eval(i-1)) < elecTolMax) then
        oDeg(iDeg) = oDeg(iDeg) + oscStrength(i)
      else
        iDeg = iDeg + 1
        eDeg(iDeg) = eval(i)
        oDeg(iDeg) = oscStrength(i)
      endif
    end do
    if (tWriteTagged) then
      call taggedWriter%write(fdTagged, tagLabels%excEgy, eDeg(:iDeg))
      call taggedWriter%write(fdTagged, tagLabels%excOsc, oDeg(:iDeg))
    end if

  end subroutine writeExcitations_rs


  !> Write atomic transition charges to file
  subroutine writeTransitionCharges_rs(fdTransQ, atomicTransQ)
    !> file unit for transition charges
    integer, intent(in)  :: fdTransQ
    !> transition charges to write
    real(dp), intent(in) :: atomicTransQ(:)

    integer :: natom, i
    natom = size(atomicTransQ)

    open(fdTransQ, file=transChargesOut, action="write", status="replace")
    write(fdTransQ, '(a)') "#"
    write(fdTransQ, '(a,2x,a,5x,a)') "#", "atom", "transition charge"
    write(fdTransQ, '(a)') "#"

    do i = 1, natom
       write(fdTransQ, '(2x,i5,5x,f12.9)') i, atomicTransQ(i)
    end do

    close(fdTransQ)

  end subroutine writeTransitionCharges_rs


  !> Computes excitation spectrum through range separated response calculation
  subroutine LinResp_calcExcitations_rs(spin, tOnsite, self, iAtomStart, eigVec,&
      & eigVal, sccCalc, SSqrReal, filling, coords0, dqAt, specie0, hubbUAtom, iNeighbor, &
      & img2CentCell, orb, rs_data, tWriteTagged, fdTagged, taggedWriter, & ! ons_en, ons_dip,
      & excEnergy, skHamCont, skOverCont, derivator, deltaRho, excgrad, dQAtomEx)
    implicit none
    logical, intent(in) :: spin
    logical, intent(inout) :: tOnsite  
    type(linresp), intent(inout) :: self
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
    type(RangeSepFunc), intent(inout) :: rs_data
    !> print tag information
    logical, intent(in) :: tWriteTagged

    !> file id for tagging information
    integer, intent(in) :: fdTagged

    !> tagged writer
    type(TTaggedWriter), intent(inout) :: taggedWriter

   !real(dp), intent(in) :: ons_en(:,:), ons_dip(:,:)
    real(dp), intent(out) :: excEnergy
!    real(dp), intent(in), optional :: shift(:)
    type(OSlakoCont), intent(in), optional :: skHamCont, skOverCont
    !> Differentiatior for the non-scc matrices
    class(NonSccDiff), intent(in), optional :: derivator
    real(dp), intent(inout), optional :: deltaRho(:,:)
    real(dp), intent(inout), optional :: excgrad(:,:)
    real(dp), intent(inout), optional :: dQAtomEx(:)

    real(dp), allocatable :: shiftPerAtom(:), shiftPerL(:,:)
    integer :: nAtom !, i
    real(dp), allocatable :: occNr(:,:)

   !real(dp), allocatable :: hubbUDerivUp(:), hubbUDerivDn(:)

    nAtom = size(orb%nOrbAtom)
    
    allocate(occNr(size(filling,dim=1),size(filling,dim=2)))
    occNr = real(filling, 4)

   !allocate(hubbUDerivUp(8))
   !allocate(hubbUDerivDn(8))
   !hubbUDerivUp(:) = 0.
   !hubbUDerivDn(:) = 0.

    if (.not. present(excgrad)) then 
       call runRs_LinRespCalc(spin, tOnsite, nAtom, iAtomStart, eigVec,&
            & eigVal, sccCalc, dqAt, coords0, self%nExc, self%nStat, self%symmetry,&
            & SSqrReal, occNr, specie0, self%nMoved,&
            & hubbUAtom, self%hubbUDerivUp, self%hubbUDerivDn, self%nEl,&
       !    & hubbUAtom, hubbUDerivUp, hubbUDerivDn, self%nEl,&
            & iNeighbor, img2CentCell, orb, rs_data, tWriteTagged, fdTagged, taggedWriter, &
       !    & self%tMulliken, self%tCoeffs, self%tXplusY, self%tTrans, self%tTradip, self%tArnoldi&
            & self%fdMulliken, self%fdCoeffs, self%fdXplusY, self%fdTrans, self%fdSPTrans, &
            & self%fdTradip, self%fdTransQ, self%tArnoldi, self%fdArnoldi, self%fdExc,& !ons_en, ons_dip,
            & self%tEnergyWindow, self%energyWindow, self%tOscillatorWindow, self%oscillatorWindow, &
            & self%tCacheCharges, excEnergy)
 
    else

       allocate(shiftPerAtom(nAtom))
       allocate(shiftPerL(orb%mShell, nAtom))
       call sccCalc%getShiftPerAtom(shiftPerAtom)
       call sccCalc%getShiftPerL(shiftPerL)
       shiftPerAtom = shiftPerAtom + shiftPerL(1,:)

       call runRs_LinRespCalc(spin, tOnsite, nAtom, iAtomStart, eigVec,&
            & eigVal, sccCalc, dqAt, coords0, self%nExc, self%nStat, self%symmetry,&
            & SSqrReal, occNr, specie0, self%nMoved,&
            & hubbUAtom, self%hubbUDerivUp, self%hubbUDerivDn, self%nEl,&
       !    & hubbUAtom, hubbUDerivUp, hubbUDerivDn, self%nEl,&
            & iNeighbor, img2CentCell, orb, rs_data, tWriteTagged, fdTagged, taggedWriter, &
       !    & self%tMulliken, self%tCoeffs, self%tXplusY, self%tTrans, self%tTradip, self%tArnoldi,&
            & self%fdMulliken, self%fdCoeffs, self%fdXplusY, self%fdTrans, self%fdSPTrans, &
            & self%fdTradip, self%fdTransQ, self%tArnoldi, self%fdArnoldi, self%fdExc,& ! ons_en, ons_dip,
            & self%tEnergyWindow, self%energyWindow, self%tOscillatorWindow, self%oscillatorWindow, &
            & self%tCacheCharges, excEnergy,&
       !    & shiftPerAtom(:,1), skHamCont, skOverCont, derivator, deltaRho, excgrad)
            & shiftPerAtom, skHamCont, skOverCont, derivator, deltaRho, excgrad, dQAtomEx)

    end if

  end subroutine LinResp_calcExcitations_rs


!Functions for gradient calculations:

  ! Create P = T + 1/2 Z symmetric (paper has T + Z asymmetric)
  ! (Zab = Zij = 0, Tia = 0)
  subroutine calcPMatrix(t, rhs, win, getij, pc)
    implicit none
    real(dp), intent(in) :: t(:,:), rhs(:)
    integer, intent(in) :: win(:), getij(:,:)
    real(dp), intent(out) :: pc(:,:)

    integer :: ia, i, a

    pc(:,:) = t(:,:)
    do ia = 1, size(rhs)
      call indxov(win, ia, getij, i, a)
      ! shouldn't be pc = pc + 1/2 rhs?
      pc(i,a) = 0.5_dp * rhs(ia)
      pc(a,i) = 0.5_dp * rhs(ia)
    end do
      
  end subroutine calcPMatrix

  !> Transform dense matrix from MO-representation into AO representation
  subroutine transformMO2AODense(aa, cc)
    implicit none
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
    implicit none
    integer, intent(in) :: iAtomStart(:)
    real(dp), intent(in) :: pc(:,:), s(:,:)
    real(dp), intent(out) :: dQAtomEx(:)

    integer :: alpha, beta, ia1, ia2, ib1, ib2
    real(dp) :: tmp

    dQAtomEx(:) = 0.0_dp
    do alpha = 1, size(dQAtomEx)
      ia1 = iAtomStart(alpha) ! TOMAS KUBAR SHIFTED INDICES BY 1
      ia2 = iAtomStart(alpha + 1) - 1
      do beta = 1, alpha
        ib1 = iAtomStart(beta) ! TOMAS KUBAR SHIFTED INDICES BY 1
        ib2 = iAtomStart(beta + 1) - 1
        tmp = sum(pc(ia1:ia2, ib1:ib2) * s(ia1:ia2, ib1:ib2))
        dQAtomEx(alpha) = dQAtomEx(alpha) + tmp
        if (alpha /= beta) then
          dQAtomEx(beta) = dQAtomEx(beta) + tmp
        end if
      end do
    end do
  end subroutine getExcMulliken

      
  ! Build right hand side of the equation for the Z-vector and those parts
  ! of the W-vectors which do not depend on Z.
  subroutine getZVectorEqRHS_rs(xpy, xmy, win, iAtomStart, homo, nocc, nMatUp, getij, iatrans, natom,&
        & species0, ev, sTimesGrndEigVecs, grndEigVecs, gamma, lrGamma, uhbuu, uhbud, omega, sym,& ! tqov, tqoo, tqvv
        & rhs, t, wov, woo, wvv)
    implicit none
    real(dp), intent(in) :: xpy(:), xmy(:)
    integer, intent(in) :: win(:), iAtomStart(:)
    integer, intent(in) :: homo, nocc, nMatUp, getij(:,:), iatrans(1:,homo+1:), natom, species0(:)
    real(dp), intent(in) :: ev(:), sTimesGrndEigVecs(:,:,:), grndEigVecs(:,:,:), gamma(:,:), lrGamma(:,:)
    real(dp), intent(in) :: uhbuu(:), uhbud(:), omega
    character, intent(in) :: sym
   !real(dp), intent(in) :: tqov(:,:), tqoo(:,:), tqvv(:,:)
    real(dp), intent(out) :: rhs(:), t(:,:)
    real(dp), intent(out) :: wov(:), woo(:), wvv(:)

    real(dp), allocatable :: xpyq(:), qij(:), gamxpyq(:), qgamxpyq(:), gamqt(:)
    real(dp), allocatable :: HvvX(:), HvvY(:), HooX(:), HooY(:), HovT(:), HooT(:)
    integer :: nvirt, norb, nxov, nxoo, nxvv
    integer :: i, j, a, b, ia, ib, ij, ab, ja
    real(dp) :: tmp1, tmp2
    logical :: lUpDwn

    nxov = size(rhs)
    nxoo = size(woo)
    nxvv = size(wvv)
   !nvirt = nxov / nocc ! this does not work if there is energyWindow (nxov_rd applies and not nxov, and nxov_rd < nxov)
   !norb = nocc + nvirt
    norb = iAtomStart(natom+1) - 1
    nvirt = norb - nocc

    @:ASSERT(nxvv == nvirt * (nvirt + 1) / 2)
    @:ASSERT(nxoo == nocc * (nocc + 1) / 2)

    allocate(xpyq(natom))
    allocate(qij(natom))
    allocate(gamxpyq(natom))
    allocate(gamqt(natom))
    allocate(qgamxpyq(max(nxoo, nxvv) * norb)) !check how much memory is really required!
   
    t(:,:) = 0.0_dp
    rhs(:) = 0.0_dp
    wov(:) = 0.0_dp
    woo(:) = 0.0_dp
    wvv(:) = 0.0_dp
    xpyq(:) = 0.0_dp

    ! optim code
    ! Build t_ab = 0.5 * sum_i (X+Y)_ia (X+Y)_ib + (X-Y)_ia (X-Y)_ib
    ! and w_ab = Q_ab with Q_ab as in (B16) but with corrected sign.
    ! factor 1 / (1 + delta_ab) follows later
    do ia = 1, nxov
      call indxov(win, ia, getij, i, a)
      lUpDwn = (win(ia) <= nMatUp)
      ! BA: Is T_aa = 0?
      do b = homo + 1, a
        ib = iatrans(i, b) 
       !if (ib <= nxov) then ! TOMAS KUBAR ADDED
          ab = (b - homo) + (((a - homo) - 1) * (a - homo))/2
          tmp1 = xpy(ia) * xpy(ib) + xmy(ia) * xmy(ib)
          tmp2 = omega * (xpy(ia) * xmy(ib)+ xmy(ia) * xpy(ib))
          t(a,b) = t(a,b) + 0.5_dp * tmp1
          ! to prevent double counting            
          if (a /= b) then
            t(b,a) = t(b,a) + 0.5_dp * tmp1 
          end if
          ! Note: diagonal elements will be multiplied by 0.5 later.
          wvv(ab) = wvv(ab) + ev(i) * tmp1 + tmp2
       !end if
      end do

      ! Build t_ij = 0.5 * sum_a (X+Y)_ia (X+Y)_ja + (X-Y)_ia (X-Y)_ja
      ! and 1 / (1 + delta_ij) Q_ij with Q_ij as in eq. (B9) (1st part of w_ij)
      do j = i, homo
        ja = iatrans(j,a) 
       !if (ja <= nxov) then ! TOMAS KUBAR ADDED
          ij = i-homo+nocc + ((j-homo+nocc - 1) * (j-homo+nocc)) / 2      
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
       !end if
      end do
    end do

    ! Build xpyq = sum_ia (X+Y)_ia 
    do ia = 1, nxov
      call indxov(win, ia, getij, i, a)
      lUpDwn = (win(ia) <= nMatUp)
      qij = transq(i, a, iAtomStart, lUpDwn, sTimesGrndEigVecs, grndEigVecs)
      xpyq(:) = xpyq + xpy(ia) * qij
    end do

    gamxpyq(:) = matmul(gamma, xpyq)
    
    ! qgamxpyq(ab) = sum_jc K_ab,jc (X+Y)_jc
    do ab = 1, nxvv
      call indxvv(homo, ab, a, b)
      lUpDwn = .true. ! UNTESTED
      qij = transq(a, b, iAtomStart, lUpDwn, sTimesGrndEigVecs, grndEigVecs)
      if (sym == "S") then
        qgamxpyq(ab) = 2.0_dp * sum(qij * gamxpyq)
      else
        qgamxpyq(ab) = sum(qij * xpyq * (uhbuu(species0) - uhbud(species0)))
      end if
    end do

    ! rhs(ia) -= Qia = sum_b (X+Y)_ib * qgamxpyq(ab))
    do ia = 1, nxov
      call indxov(win, ia, getij, i, a)
      lUpDwn = (win(ia) <= nMatUp)
      do b = homo + 1, a
        ! TOMAS KUBAR -- CORRECTED:
        !call rindxvv(homo, a, b, ab)
        if (b < a) then
           call rindxvv(homo, a, b, ab)
        else
           call rindxvv(homo, b, a, ab)
        end if
        ib = iatrans(i,b)
       !if (ib <= nxov) then ! TOMAS KUBAR ADDED
          rhs(ia) = rhs(ia) - 2.0_dp * xpy(ib) * qgamxpyq(ab)
          ! Since qgamxpyq has only upper triangle
          if (a /= b) then
            rhs(ib) = rhs(ib) - 2.0_dp * xpy(ia) * qgamxpyq(ab)
          end if
       !end if
      end do
    end do

    ! -rhs = -rhs - sum_j (X + Y)_ja H + _ij[X + Y]
    ! qgamxpyq(ij) = sum_kb K_ij,kb (X+Y)_kb
    do ij = 1, nxoo
      qgamxpyq(ij) = 0.0_dp
      call indxoo(ij, i, j)
      lUpDwn = .true.
      qij = transq(i, j, iAtomStart, lUpDwn, sTimesGrndEigVecs, grndEigVecs)
      if (sym == "S") then
        qgamxpyq(ij) = 2.0_dp * sum(qij * gamxpyq)
      else
        qgamxpyq(ij) = sum(qij * xpyq * (uhbuu(species0) - uhbud(species0)))
      end if
    end do

    ! rhs(ia) += Qai = sum_j (X+Y)_ja qgamxpyq(ij)
    ! add Qai to Wia as well.
    do ia = 1, nxov
      call indxov(win, ia, getij, i, a)
     !lUpDwn = (win(ia) <= nMatUp)
      do j = i, homo
        ja = iatrans(j, a)
       !if (ja <= nxov) then ! TOMAS KUBAR ADDED
          ij = i-homo+nocc + ((j-homo+nocc - 1) * (j-homo+nocc)) / 2
          tmp1 = 2.0_dp * xpy(ja) * qgamxpyq(ij)
          rhs(ia) = rhs(ia) + tmp1
          wov(ia) = wov(ia) + tmp1
          if (i /= j) then
            tmp2 = 2.0_dp * xpy(ia) * qgamxpyq(ij)
            rhs(ja) = rhs(ja) + tmp2
            wov(ja) = wov(ja) + tmp2
          end if
       !end if
      end do
    end do

    ! gamxpyq(beta) = sum_ij q_ij(beta) T_ij 
    gamxpyq(:) = 0.0_dp
    do ij = 1, nxoo
      call indxoo(ij, i, j)
      lUpDwn = .true.
      qij = transq(i, j, iAtomStart, lUpDwn, sTimesGrndEigVecs, grndEigVecs)
      if (i == j) then
        gamxpyq(:) = gamxpyq(:) + t(i,j) * qij(:)
      else
        ! factor 2 because of symmetry of the matrix
        gamxpyq(:) = gamxpyq(:) + 2.0_dp  * t(i,j) * qij(:)
      end if
    end do

    ! gamxpyq(beta) += sum_ab q_ab(beta) T_ab 
    do ab = 1, nxvv
      call indxvv(homo, ab, a, b)
      lUpDwn = .true.
      qij = transq(a, b, iAtomStart, lUpDwn, sTimesGrndEigVecs, grndEigVecs)
      if (a == b) then
        gamxpyq(:) = gamxpyq(:) + t(a,b) * qij(:)
      else
        ! factor 2 because of symmetry of the matrix
        gamxpyq(:) = gamxpyq(:) + 2.0_dp * t(a,b) * qij(:)
      end if
    end do

    ! gamqt(alpha) = sum_beta gamma_alpha,beta gamxpyq(beta)
    gamqt(:) = matmul(gamma, gamxpyq)

    ! rhs -= sum_q^ia(alpha) gamxpyq(alpha)
    do ia = 1, nxov
      call indxov(win, ia, getij, i, a)
      lUpDwn = (win(ia) <= nMatUp)
      qij = transq(i, a, iAtomStart, lUpDwn, sTimesGrndEigVecs, grndEigVecs)
      rhs(ia) = rhs(ia) - 4.0_dp * sum(qij * gamqt)
    end do
    
    ! Furche vectors
    do ij = 1, nxoo
      call indxoo(ij, i, j)
      lUpDwn = .true.
      qij = transq(i, j, iAtomStart, lUpDwn, sTimesGrndEigVecs, grndEigVecs)
      woo(ij) = woo(ij) + 4.0_dp * sum(qij * gamqt)
    end do
    
    ! Contributions due to range-separation
    allocate(HvvX(nxvv))
    allocate(HvvY(nxvv))
    allocate(HooX(nxoo))
    allocate(HooY(nxoo))
    allocate(HovT(nxov))

    call GetHvvXY(  1, nxov, nxvv,       norb, homo, natom, nMatUp, iatrans, getij, win,&
                & iAtomStart, sTimesGrndEigVecs, grndEigVecs, lrGamma, xpy, HvvX)
    call GetHvvXY( -1, nxov, nxvv,       norb, homo, natom, nMatUp, iatrans, getij, win,&
                & iAtomStart, sTimesGrndEigVecs, grndEigVecs, lrGamma, xmy, HvvY)
    call GetHooXY(  1, nxov, nxoo,       norb, homo, natom, nMatUp, iatrans, getij, win,&
                & iAtomStart, sTimesGrndEigVecs, grndEigVecs, lrGamma, xpy, HooX)
    call GetHooXY( -1, nxov, nxoo,       norb, homo, natom, nMatUp, iatrans, getij, win,&
                & iAtomStart, sTimesGrndEigVecs, grndEigVecs, lrGamma, xmy, HooY)
    call GetHovT(nocc, nxov, nxoo, nxvv, norb, homo, natom, nMatUp, iatrans, getij, win,&
                & iAtomStart, sTimesGrndEigVecs, grndEigVecs, lrGamma, t,   HovT)

    do ia = 1, nxov
      call indxov(win, ia, getij, i, a)
      lUpDwn = (win(ia) <= nMatUp)
      do b = homo + 1, norb
         ib = iatrans(i, b)
         if (ib <= nxov) then ! TOMAS KUBAR ADDED
            ! TOMAS KUBAR -- CORRECTED
            !call rindxvv(nocc, a, b, ab)
            if (b < a) then
               call rindxvv(homo, a, b, ab)
            else
               call rindxvv(homo, b, a, ab)
            end if
            rhs(ia) = rhs(ia) - xpy(ib)*HvvX(ab)
            if (a >= b) then
               rhs(ia) = rhs(ia) - xmy(ib)*HvvY(ab)
            else
               rhs(ia) = rhs(ia) + xmy(ib)*HvvY(ab)
            endif
         end if
      enddo

      do j = 1, homo
         ja = iatrans(j, a)
         if (ja <= nxov) then ! TOMAS KUBAR ADDED
            !call rindxoo(homo, nocc, i, j, ij)
            if (j < i) then
               call rindxoo(homo, nocc, i, j, ij)
            else
               call rindxoo(homo, nocc, j, i, ij)
            end if
            rhs(ia) = rhs(ia) + xpy(ja)*HooX(ij)
            wov(ia) = wov(ia) + xpy(ja)*HooX(ij)
            if (i >= j) then
               rhs(ia) = rhs(ia) + xmy(ja)*HooY(ij)
               wov(ia) = wov(ia) + xmy(ja)*HooY(ij)
            else
               rhs(ia) = rhs(ia) - xmy(ja)*HooY(ij)
               wov(ia) = wov(ia) - xmy(ja)*HooY(ij)
            endif
         end if
      enddo
      rhs(ia) = rhs(ia) - HovT(ia)  
    enddo
    deallocate(HvvX)
    deallocate(HvvY)
    deallocate(HooX)
    deallocate(HooY)
    deallocate(HovT)

    allocate(HooT(nxoo))
    call GetHooT(nocc, nxov, nxoo, norb, homo, natom, nMatUp, iatrans, getij, win,&
               & iAtomStart, sTimesGrndEigVecs, grndEigVecs, lrGamma, t, HooT)
    woo(:) = woo(:) + HooT(:) 
    deallocate(HooT)

    deallocate(xpyq)
    deallocate(qij)
    deallocate(gamxpyq)
    deallocate(gamqt)
    deallocate(qgamxpyq) 
  end subroutine getZVectorEqRHS_rs


  ! Solving the (A+B) Z = -R equation via conjugate gradient optimization.
 !subroutine solveZVectorEq_rs(rhs, win, nMatUp, getij, occNr, natom, iAtomStart, sTimesGrndEigVecs,&
 !    & species0, uhbuu, uhbud, gamma, lrGamma, wij, grndEigVecs,&
 !    & orb, homo, nocc, nvir, iatrans, tqov, tqoo, tqvv)
  subroutine solveZVectorEq_rs(rhs, win, nMatUp, getij, occNr, natom, &
      & species0, uhbuu, uhbud, gamma, lrGamma, wij, &
      & homo, nocc, nvir, iatrans, tqov, tqoo, tqvv)
    implicit none
    real(dp), intent(inout) :: rhs(:)
    integer, intent(in) :: win(:), nMatUp, getij(:,:), natom !, iAtomStart(:)
    real(dp), intent(in) :: occNr(:,:)
    real(dp), intent(in) :: gamma(:,:), lrGamma(:,:), wij(:) !, grndEigVecs(:,:,:), sTimesGrndEigVecs(:,:,:)
    integer, intent(in) :: species0(:)
    real(dp), intent(in) :: uhbuu(:), uhbud(:)
   !type(TOrbitals), intent(in) :: orb
    integer, intent(in) :: homo, nocc, nvir
    real(dp), intent(in) :: tqov(:,:), tqoo(:,:), tqvv(:,:)
    integer, intent(in) :: iatrans(1:,homo+1:)

    integer :: nxov
    integer :: ia, i, a, k, ii, aa
    real(dp), allocatable :: qij(:), gxqij(:), rhs2(:), rkm1(:), pkm1(:), apk(:)
    real(dp) :: rs, tmp1, tmp2
    logical :: lUpDwn
    real(dp), allocatable :: gqvtmp(:,:)

    !dummy variables to cal. (A+B)*vec
    character :: sym
    logical :: spin
    sym = 'S'
    spin = .false.
   
    if (sum(rhs(:) * rhs(:)) .le. 0.000001) then !May actually be zero, for example for H_2
       rhs(:) = 0.0_dp
       return
    end if

    nxov = size(rhs)
    allocate(qij(natom))
    allocate(gxqij(natom))
    allocate(rhs2(nxov))
    allocate(rkm1(nxov))
    allocate(pkm1(nxov))
    allocate(apk(nxov))
    allocate(gqvtmp(size(gamma, dim=1), max(nvir, nocc) * nvir)) 
    
    ! Choosing a start value
    ! rhs2 = rhs / (A+B)_ia,ia (diagonal of the supermatrix sum A+B)
    do ia = 1, nxov
      call indxov(win, ia, getij, i, a)
      lUpDwn = (win(ia) <= nMatUp)
   !   qij = transq(i, a, iAtomStart, lUpDwn, sTimesGrndEigVecs, grndEigVecs)
      rs = 4.0_dp * dot_product(tqov(:,ia), matmul(gamma, tqov(:,ia))) + wij(ia)
      rs = rs - dot_product(tqov(:,ia), matmul(lrGamma, tqov(:,ia))) !part of rs contirb
      !!Bug here! ii = i - homo + nocc + ((i - homo + nocc + 1) * (i - homo + nocc))/2
      ii = i - homo + nocc + ((i - homo + nocc - 1) * (i - homo + nocc))/2
      !lUpDwn = .true.
      !qij = transq(i, i, iAtomStart, lUpDwn, sTimesGrndEigVecs, grndEigVecs)
      gxqij(:) = matmul(lrGamma, tqoo(:,ii))
      call rindxvv(homo, a, a, aa)
      !lUpDwn = .true.
      !qij = transq(a, a, iAtomStart, lUpDwn, sTimesGrndEigVecs, grndEigVecs)
      rs = rs - dot_product(tqvv(:,aa), gxqij)
      rhs2(ia) = rhs(ia) / rs
    end do

    call multApBVecFast(rhs2, wij, sym, win, nMatUp, homo, nocc, nvir, & ! sTimesGrndEigVecs, 
        & occNr, getij, gamma, lrGamma, species0, uhbuu, uhbud, iatrans, gqvtmp, tqov, tqoo, tqvv, rkm1)
    
    rkm1(:) = rhs(:) - rkm1(:)
    pkm1(:) = rkm1(:)
    
    ! Iteration:  should be convergent at most nxov steps
    do k = 1, nxov

      call multApBVecFast(pkm1, wij, sym, win, nMatUp, homo, nocc, nvir, & ! sTimesGrndEigVecs,
          & occNr, getij, gamma, lrGamma, species0, uhbuu, uhbud, iatrans, gqvtmp, tqov, tqoo, tqvv, apk)
      tmp1 = dot_product(rkm1, rkm1)
      tmp2 = dot_product(pkm1, apk)

      rhs2(:) = rhs2(:) + (tmp1 / tmp2) * pkm1(:)
      rkm1(:) = rkm1(:) - (tmp1 / tmp2) *  apk(:)

      tmp2 = dot_product(rkm1, rkm1)
!       print *, "residuo", tmp2 

      if (tmp2 <= 1.0e-14_dp) then
        exit
      end if

      if (k == nxov) then
        call error("LrespoGrad : Z vector not converged!")
      end if

      pkm1(:) = (tmp2 / tmp1) * pkm1(:) + rkm1(:)

    end do

    rhs(:) = rhs2(:)   

    deallocate(qij)
    deallocate(gxqij)
    deallocate(rhs2)
    deallocate(rkm1)
    deallocate(pkm1)
    deallocate(apk)
    deallocate(gqvtmp)     
  end subroutine solveZVectorEq_rs

  ! Calculate Z-dependent parts of the W-vectors including rs contributions 
  !  and divide diagonal elements of W_ij and W_ab by 2.
  !subroutine calcWVectorZ_rs(zz, win, homo, nocc, nvirt, norb, nMatUp, getij, iAtomStart, sTimesGrndEigVecs,&
  !    & grndEigVecs, gamma, lrGamma, ev, wov, woo, wvv, tqov, tqoo, iatrans)
  subroutine calcWVectorZ_rs(zz, win, homo, norb, nMatUp, getij, iAtomStart, sTimesGrndEigVecs,&
      & grndEigVecs, gamma, lrGamma, ev, wov, woo, wvv, iatrans)
    implicit none
    real(dp), intent(in) :: zz(:)
    integer, intent(in) :: win(:), homo, norb, nMatUp, getij(:,:), iAtomStart(:) !, nocc, nvirt
    real(dp), intent(in) :: sTimesGrndEigVecs(:,:,:), grndEigVecs(:,:,:), gamma(:,:), lrGamma(:,:), ev(:) !, tqov(:,:), tqoo(:,:)
    real(dp), intent(inout) :: wov(:), woo(:), wvv(:)
    integer, intent(in) :: iatrans(1:,homo+1:)

    integer :: nxov, nxoo, nxvv, natom
    integer :: ij, ia, ab, i, j, a, b, alpha
    real(dp), allocatable :: qij(:), gamxpyq(:), zq(:), HooZ(:)
    logical :: lUpDwn

    nxov = size(zz)
    natom = size(gamma, dim=1)   
    nxoo = size(woo)
    nxvv = size(wvv)
    allocate(qij(natom))
    allocate(gamxpyq(natom))
    allocate(zq(natom))
    !for rangesep contributions:
    allocate(HooZ(nxoo))

    ! Adding missing epsilon_i * Z_ia term to W_ia
    do ia = 1, nxov
      call indxov(win, ia, getij, i, a)
      lUpDwn = (win(ia) <= nMatUp)
      wov(ia) = wov(ia) + zz(ia) * ev(i)
    end do

    ! Missing sum_kb 4 K_ijkb Z_kb term in W_ij:
    ! zq(alpha) = sum_kb q^kb(alpha) Z_kb
    do alpha = 1, natom
      zq(alpha) = 0.0_dp
      do ia = 1, nxov
        call indxov(win, ia, getij, i, a)
        lUpDwn = (win(ia) <= nMatUp)
        qij = transq(i, a, iAtomStart, lUpDwn, sTimesGrndEigVecs, grndEigVecs)
        zq(alpha) = zq(alpha) + zz(ia) * qij(alpha)
      end do
    end do

    ! gamxpyq(alpha) = sum_beta gamma(alpha, beta) zq(beta)
    gamxpyq = matmul(gamma, zq)

    ! sum_alpha qij(alpha) gamxpyq(alpha)
    do ij = 1, nxoo
      call indxoo(ij, i, j)
      lUpDwn = .true.
      qij = transq(i, j, iAtomStart, lUpDwn, sTimesGrndEigVecs, grndEigVecs)
      do alpha = 1, natom
        ! W contains 1/2 for i == j.
        woo(ij) = woo(ij)&
            & + 4.0_dp * qij(alpha) * gamxpyq(alpha)
      end do
    end do

    ! Contributions due to range-separation
    call GetHooXY(1, nxov, nxoo, norb, homo, natom, nMatUp, iatrans, getij, win, iAtomStart,&
                & sTimesGrndEigVecs, grndEigVecs, lrGamma, zz, HooZ)
    woo(:) = woo(:) + HooZ(:)

    ! Divide diagonal elements of W_ij by 2.
    do ij = 1, nxoo
      call indxoo(ij, i, j)
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

    deallocate(HooZ)
    deallocate(qij)
    deallocate(gamxpyq)
    deallocate(zq)
  end subroutine calcWVectorZ_rs

 !subroutine addGradients_rs(sym, nxov, natom, species0, iAtomStart, norb, homo, nocc, nMatUp, getij, win,& ! iatrans,
 !    & grndEigVecs, pc, sTimesGrndEigVecs, dQ, dQAtomEx, gamma, lrGamma, rs_data, uhubb, uhbuu, uhbud, shift,&
 !    & woo, wov, wvv, xpy, xmy, nbeweg, coord0, orb, skHamCont,&
 !    & skOverCont, tqov, tqoo, tqvv, derivator, deltaRho, excgrad)
  subroutine addGradients_rs(sym, nxov, natom, species0, iAtomStart, norb, homo, nMatUp, getij, win,& ! iatrans,
      & grndEigVecs, pc, dQ, dQAtomEx, gamma, lrGamma, rs_data, uhubb, uhbuu, uhbud, shift,&
      & woo, wov, wvv, xpy, xmy, nbeweg, coord0, orb, skHamCont,&
      & skOverCont, tqov, derivator, deltaRho, excgrad)
    implicit none
    character, intent(in) :: sym
    integer, intent(in) :: nxov, natom, species0(:), iAtomStart(:), norb, homo, win(:) !, nocc
    integer, intent(in) :: nMatUp, getij(:,:) !, iatrans(1:,homo+1:)
    real(dp), intent(in) :: grndEigVecs(:,:,:), pc(:,:), dQ(:), dQAtomEx(:) !, sTimesGrndEigVecs(:,:,:)
    real(dp), intent(in) :: gamma(:,:), lrGamma(:,:), uhubb(:), uhbuu(:), uhbud(:), shift(:)
    type(RangeSepFunc), intent(inout) :: rs_data
    real(dp), intent(in) :: woo(:), wov(:), wvv(:), xpy(:), xmy(:)
    integer, intent(in) :: nbeweg
    real(dp), intent(in) :: coord0(:,:)
    type(TOrbitals), intent(in) :: orb
    type(OSlakoCont), intent(in) :: skHamCont, skOverCont
    real(dp), intent(in) :: tqov(:,:) !, tqoo(:,:), tqvv(:,:)
    !> Differentiatior for the non-scc matrices
    class(NonSccDiff), intent(in) :: derivator
    real(dp), intent(inout) :: excgrad(:,:), deltaRho(:,:)
    
    real(dp), allocatable :: shex(:), xpyq(:), shxpyq(:), xpycc(:,:), xmycc(:,:), wcc(:,:)
    real(dp), allocatable :: temp(:), qij(:), xpyas(:,:), xmyas(:,:)
   !real(dp), allocatable :: au(:,:), bu(:,:), auh(:,:), buh(:,:),
    real(dp), allocatable :: dmn(:,:)
    integer :: ia, i, j, a, b, ab, ij, m, n, mu, nu, xyz, alpha, beta
    integer :: indalpha, indalpha1, indbeta, indbeta1, izpalpha, izpbeta
    real(dp) :: tmp1, tmp2, tmp3, tmp4, tmp5, tmp6, tmp7, tmprs, tmprs2, rab !, xhelp
    real(dp) :: diffvec(3), dgab(3), tmp3a(3), tmp3b(3)
    real(dp) :: dhmndr, dsmndr
    real(dp), parameter :: deltax = epsilon(1.0_dp)**0.25_dp
    real(dp), parameter :: rcdx = 1.0_dp / deltax
    !Working arrays for the Hamilton/Overlap builder to make things faster.
   !real(dp) :: interSKHam(getMIntegrals(skHamCont))
   !real(dp) :: interSKOver(getMIntegrals(skOverCont))
    real(dp), allocatable :: overlap(:,:), lrGammaOrb(:,:)
    real(dp), allocatable :: PS(:,:), DS(:,:), SPS(:,:), SDS(:,:), gammaLongRangePrime(:,:,:)
    real(dp), allocatable :: SX(:,:), XS(:,:), SXS(:,:), SY(:,:), YS(:,:), SYS(:,:)
    real(dp) tmpura(3)
    integer iAt1, iAt2, ka

    integer, allocatable :: species(:)

    real(dp), allocatable :: derH0(:,:,:), derS(:,:,:)
    
    logical :: lUpDwn
    
    integer :: nxoo, nxvv
    allocate(shex(natom))
    allocate(xpyq(natom))
    allocate(shxpyq(natom))
    allocate(xpycc(norb, norb))
    allocate(xmycc(norb, norb))
    allocate(xpyas(norb, norb))
    allocate(xmyas(norb, norb))
    allocate(wcc(norb, norb))
    allocate(qij(natom))
    allocate(temp(norb))
    allocate(dmn(norb, norb))

!!  RS matrices
    allocate(PS(norb, norb))
    allocate(DS(norb, norb))
    allocate(SPS(norb, norb))    
    allocate(SDS(norb, norb))
    allocate(SX(norb, norb))
    allocate(XS(norb, norb)) 
    allocate(SXS(norb, norb))
    allocate(SY(norb, norb))
    allocate(YS(norb, norb))
    allocate(SYS(norb, norb))
    allocate(overlap(norb, norb))
    allocate(lrGammaOrb(norb, norb))
    allocate(gammaLongRangePrime(3,natom,natom))

    allocate(derH0(orb%mOrb, orb%mOrb, 3))
    allocate(derS(orb%mOrb, orb%mOrb, 3))

    nxoo = size(woo)
    nxvv = size(wvv)

     !write (*,*) "C"
     !write (*,'(10F12.8)') c
     !write (*,*)

    !Density matrix by rank k-update
    !BA: density matrix should be provided from outside
    dmn = 0._dp
    call herk(dmn, grndEigVecs(:,1:homo,1), alpha=2.0_dp)

     !write (*,*) "HERK DMN"
     !write (*,'(22F12.8)') dmn
     !write (*,*)

    ! TEST OK -- DMN (UPPER TRIANGLE ONLY, WHICH IS FINE)

    !! Symmetrize deltaRho
    do mu = 1, norb
       do nu = mu + 1, norb
          deltaRho(mu, nu) = deltaRho(nu, mu)
       enddo
    enddo

     !write (*,*) "DELTA RHO SYMM"
     !write (*,'(22F12.8)') deltaRho
     !write (*,*)

    ! TEST OK -- DELTA RHO SYMM

    !! Compute long-range gamma derivative
    gammaLongRangePrime(:,:,:) = 0._dp
    call rs_data%getSpecies(species)
    do iAt1 = 1, natom
       do iAt2 = 1, natom
          if(iAt1 /= iAt2) then
             call getGammaPrimeValue(rs_data, tmpura, iAt1, iAt2, coord0, species)
             gammaLongRangePrime(:,iAt1, iAt2) = tmpura(:)
          end if
       end do
    end do

     !write (*,*) "GAMMA LR PRIME"
     !write (*,'(10F14.10)') gammaLongRangePrime
     !write (*,*)

     ! Symmetrize S
   !call getSqrS(coord0, natom, skOverCont, interSKOver, orb, iAtomStart, species0, overlap)
    call getSqrS(coord0, natom, skOverCont, orb, iAtomStart, species0, overlap)
    call getSqrGamma(natom, lrGamma, iAtomStart, lrGammaOrb)

     !write (*,*) "OVERLAP"
     !write (*,'(10F12.8)') overlap
     !write (*,*)

     !write (*,*) "GAMMA LR ORB"
     !write (*,'(10F12.8)') lrGammaOrb
     !write (*,*)

    shex(:) = matmul(dQAtomEx, gamma)

     !write (*,*) "SHEX"
     !write (*,'(10F12.8)') shex
     !write (*,*)

    ! TEST OK -- OVERLAP, GAMMA LR ORB, SHEX

    !xypq(alpha) = sum_ia (X+Y)_ia q^ia(alpha)
    !complexity norb * norb * norb
    xpyq(:) = 0.0_dp
    do ia = 1, nxov
!      call indxov(win, ia, getij, i, a)
!      lUpDwn = (win(ia) <= nMatUp)
!      qij = transq(i, a, iAtomStart, lUpDwn, sTimesGrndEigVecs, c)
      xpyq(:) = xpyq(:) + xpy(ia) * tqov(:,ia)! * qij(:)
    end do

    !complexity norb * norb
    shxpyq(:) = 0.0_dp
    if (sym == "S") then
      shxpyq(:) = matmul(xpyq, gamma)
    else
      shxpyq(:) = 0.5_dp * xpyq(:) * (uhbuu(species0) - uhbud(species0))
    end if

    !calculate xpycc
    !(xpycc)_{mu nu} = 
    !=  sum_{ia} (X + Y)_{ia} (c(mu,i)c(nu,a) + c(nu,i)c(mu,a))
    !complexity norb * norb * norb
    xpycc(:,:) = 0.0_dp
    xmycc(:,:) = 0.0_dp
    xpyas(:,:) = 0.0_dp
    xmyas(:,:) = 0.0_dp   

    !xpycc(mu,nu) = sum_ia (X+Y)_ia c(mu,i) c(nu,a)
    !xpycc(mu, nu) += sum_ia (X+Y)_ia c(mu,a) c(nu,i)    
    do ia = 1, nxov
      call indxov(win, ia, getij, i, a)
      lUpDwn = (win(ia) <= nMatUp)
      do nu = 1, norb
        do mu = 1, norb
          xpycc(mu,nu) = xpycc(mu,nu) + xpy(ia) *&
           &( grndEigVecs(mu,i,1)*grndEigVecs(nu,a,1) + grndEigVecs(mu,a,1)*grndEigVecs(nu,i,1) )
          xmycc(mu,nu) = xmycc(mu,nu) + xmy(ia) *&
           &( grndEigVecs(mu,i,1)*grndEigVecs(nu,a,1) + grndEigVecs(mu,a,1)*grndEigVecs(nu,i,1) )
          xpyas(mu,nu) = xpyas(mu,nu) + xpy(ia) * grndEigVecs(mu,i,1) * grndEigVecs(nu,a,1) 
          xmyas(mu,nu) = xmyas(mu,nu) + xmy(ia) * grndEigVecs(mu,i,1) * grndEigVecs(nu,a,1) 
        end do
      end do        
    end do

    !calculate wcc = c_mu,i * W_ij * c_j,nu
    !We have only W_ab b > a and W_ij j > i:
    !wcc(m,n) = sum_{pq, p <= q} w_pq (c(mu,p)c(nu,q) + c(nu,p)c(mu,q))
    !complexity norb * norb * norb

    !calculate the occ-occ part
    !BA: Does not the diagonal contain twice as much as needed?
    wcc(:,:) = 0.0_dp
    
    do ij = 1, nxoo
      call indxoo(ij, i, j)
      do mu = 1, norb
        do nu = 1, norb
          wcc(mu,nu) = wcc(mu,nu) + woo(ij) *&
           &( grndEigVecs(mu,i,1)*grndEigVecs(nu,j,1) + grndEigVecs(mu,j,1)*grndEigVecs(nu,i,1) )
        end do
      end do  
    end do

    !calculate the occ-virt part : the same way as for xpycc
    
    do ia = 1, nxov
      call indxov(win, ia, getij, i, a)
      lUpDwn = (win(ia) <= nMatUp)
      do nu = 1, norb
        do mu = 1, norb
          wcc(mu,nu) = wcc(mu,nu) + wov(ia) *&
           &( grndEigVecs(mu,i,1)*grndEigVecs(nu,a,1) + grndEigVecs(mu,a,1)*grndEigVecs(nu,i,1) )
        end do
      end do        
    end do
    
    !calculate the virt - virt part
    
    do ab =1, nxvv
      call indxvv(homo, ab, a, b)
      do mu = 1, norb
        do nu = 1, norb
          wcc(mu,nu) = wcc(mu,nu) + wvv(ab) *&
           &( grndEigVecs(mu,a,1)*grndEigVecs(nu,b,1) + grndEigVecs(mu,b,1)*grndEigVecs(nu,a,1) )
        end do
      end do 
    end do  

    !now calculating the force !
    !complexity : norb * norb * 3

    call symm( PS, 'R', overlap, pc,       'U', 1.0_dp, 0.0_dp, norb, norb)
    call symm(SPS, 'L', overlap, PS,       'U', 1.0_dp, 0.0_dp, norb, norb)
    call symm( DS, 'R', overlap, deltaRho, 'U', 1.0_dp, 0.0_dp, norb, norb)
    call symm(SDS, 'L', overlap, DS,       'U', 1.0_dp, 0.0_dp, norb, norb)
    call symm( XS, 'R', overlap, xpyas,    'U', 1.0_dp, 0.0_dp, norb, norb)
    call symm(SX,  'L', overlap, xpyas,    'U', 1.0_dp, 0.0_dp, norb, norb)
    call symm(SXS, 'L', overlap, XS,       'U', 1.0_dp, 0.0_dp, norb, norb)
    call symm( YS, 'R', overlap, xmyas,    'U', 1.0_dp, 0.0_dp, norb, norb)
    call symm(SY,  'L', overlap, xmyas,    'U', 1.0_dp, 0.0_dp, norb, norb)
    call symm(SYS, 'L', overlap, YS,       'U', 1.0_dp, 0.0_dp, norb, norb)

    !BA: only for non-periodic systems!
    do alpha = 1, nbeweg
      ! TOMAS KUBAR MODIFIED
      !indalpha = iAtomStart(alpha) + 1
      !indalpha1 = iAtomStart(alpha + 1) 
      indalpha = iAtomStart(alpha)
      indalpha1 = iAtomStart(alpha + 1) - 1
      izpalpha = species0(alpha)
      do beta = 1, alpha - 1
        ! TOMAS KUBAR MODIFIED
        !indbeta = iAtomStart(beta) + 1
        !indbeta1 = iAtomStart(beta + 1) 
        indbeta = iAtomStart(beta)
        indbeta1 = iAtomStart(beta + 1) - 1
        izpbeta = species0(beta)
        diffvec = coord0(:,alpha) - coord0(:,beta)
        rab = sqrt(sum(diffvec**2))
        !Note: Here dgamma/dr * 1/r is needed as returned by old gam121
        tmp1 = (-1.0_dp / rab**2&
            & - expGammaPrime(rab, uhubb(izpalpha), uhubb(izpbeta)))&
            & / rab
        !calculate the derivative of gamma 
        dgab(:) = diffvec(:) * tmp1
        tmp3a(:) = dgab(:) * (dQ(alpha) * dQAtomEx(beta) + dQAtomEx(alpha) * dQ(beta))
        if (sym == "S") then
          tmp3b(:) = 4.0_dp * dgab(:) * xpyq(alpha) * xpyq(beta)
        else
          tmp3b(:) = 0.0_dp
        end if
        excgrad(:,alpha) = excgrad(:,alpha) + tmp3a + tmp3b
        excgrad(:,beta) = excgrad(:,beta) - tmp3a - tmp3b
        tmp5 = shex(beta) + shex(alpha)
        tmp7 = 2.0_dp * shxpyq(beta) + 2.0_dp * shxpyq(alpha)

        tmprs = 0.0_dp
        do mu = indalpha, indalpha1
            do nu = indbeta, indbeta1
               tmprs  = tmprs +           &
                          &   ( 2.0_dp*(PS(mu,nu)*DS(nu,mu) + PS(nu,mu)*DS(mu,nu))      + &
                          &     SPS(mu,nu)*deltaRho(mu,nu) + SPS(nu,mu)*deltaRho(nu,mu) + &
                          &     pc(mu,nu)*SDS(mu,nu) + pc(nu,mu)*SDS(nu,mu)             )
               tmprs  = tmprs +  2.0_dp * &
                          &   ( xpyas(mu,nu)*SXS(mu,nu) + xpyas(nu,mu)*SXS(nu,mu)       + &
                          &     SX(mu,nu)*XS(mu,nu) + SX(nu,mu)*XS(nu,mu)               ) 
               tmprs  = tmprs +           &
                          &   ( XS(mu,nu)*XS(nu,mu) + XS(nu,mu)*XS(mu,nu)               + &
                          &     SXS(mu,nu)*xpyas(nu,mu) + SXS(nu,mu)*xpyas(mu,nu)       + &
                          &     xpyas(mu,nu)*SXS(nu,mu) + xpyas(nu,mu)*SXS(mu,nu)       + &
                          &     SX(mu,nu)*SX(nu,mu) + SX(nu,mu)*SX(mu,nu)               ) 
               tmprs  = tmprs +  2.0_dp * &
                          &   ( xmyas(mu,nu)*SYS(mu,nu) + xmyas(nu,mu)*SYS(nu,mu)       + &
                          &     SY(mu,nu)*YS(mu,nu) + SY(nu,mu)*YS(nu,mu)               )
               tmprs  = tmprs -           &
                          &   ( YS(mu,nu)*YS(nu,mu) + YS(nu,mu)*YS(mu,nu)               + &
                          &     SYS(mu,nu)*xmyas(nu,mu) + SYS(nu,mu)*xmyas(mu,nu)       + &
                          &     xmyas(mu,nu)*SYS(nu,mu) + xmyas(nu,mu)*SYS(mu,nu)       + &
                          &     SY(mu,nu)*SY(nu,mu) + SY(nu,mu)*SY(mu,nu) )                  
            enddo 
        enddo

        excgrad(:,alpha) = excgrad(:,alpha) - 0.125_dp * tmprs * gammaLongRangePrime(:,alpha,beta)
        excgrad(:,beta)  = excgrad(:,beta)  + 0.125_dp * tmprs * gammaLongRangePrime(:,alpha,beta)

        call derivator%getFirstDeriv(derH0, skHamCont, coord0, species0, alpha, beta, orb)
        call derivator%getFirstDeriv(derS, skOverCont, coord0, species0, alpha, beta, orb)

        do xyz = 1, 3

          tmp1 = 0.0_dp
          tmp2 = 0.0_dp
          tmp3 = 0.0_dp
          tmp4 = 0.0_dp
          tmp6 = 0.0_dp
          tmprs2 = 0.0_dp
          do mu = indalpha, indalpha1
            do nu = indbeta, indbeta1
              m = mu - indalpha + 1
              n = nu - indbeta + 1

              dhmndr = derH0(n,m,xyz)
              dsmndr = derS(n,m,xyz)

              tmp1 = tmp1 + 2.0_dp * dhmndr * pc(mu,nu)
              tmp2 = tmp2 + dsmndr * pc(mu,nu) *&
                  & (shift(alpha) + shift(beta))
              tmp3 = tmp3 - dsmndr * wcc(mu,nu)
              tmp4 = tmp4 + tmp5 * dsmndr * dmn(mu,nu)
              tmp6 = tmp6 + tmp7 * dsmndr * xpycc(mu,nu)

              tmprs = 0.0_dp
              do ka = 1, norb
                 tmprs = tmprs +   &
                      &  ( PS(mu,ka)*deltaRho(nu,ka) + PS(nu,ka)*deltaRho(mu,ka)   + &
                      &    pc(mu,ka)*DS(nu,ka) + pc(nu,ka)*DS(mu,ka)             ) * &
                      &             (lrGammaOrb(mu,ka) + lrGammaOrb(nu,ka)) 
                 tmprs = tmprs +   &
                      &  ( xpyas(mu,ka)*XS(nu,ka) + xpyas(ka,mu)*SX(ka,nu)         + &
                      &    xpyas(nu,ka)*XS(mu,ka) + xpyas(ka,nu)*SX(ka,mu)       ) * &
                      &             (lrGammaOrb(mu,ka) + lrGammaOrb(nu,ka)) 
                 tmprs = tmprs +   &
                      &  ( xmyas(mu,ka)*YS(nu,ka) + xmyas(ka,mu)*SY(ka,nu)         + &
                      &    xmyas(nu,ka)*YS(mu,ka) + xmyas(ka,nu)*SY(ka,mu)       ) * &
                      &  (lrGammaOrb(mu,ka) + lrGammaOrb(nu,ka))     
                 tmprs = tmprs +   &
                      &  ( XS(mu,ka)*xpyas(ka,nu) + XS(nu,ka)*xpyas(ka,mu)         + &
                      &    xpyas(mu,ka)*SX(ka,nu) + xpyas(nu,ka)*SX(ka,mu)       ) * &
                      &  (lrGammaOrb(mu,ka) + lrGammaOrb(nu,ka)) 
                 tmprs = tmprs -   &
                      &  ( YS(mu,ka)*xmyas(ka,nu) + YS(nu,ka)*xmyas(ka,mu)         + &
                      &    xmyas(mu,ka)*SY(ka,nu) + xmyas(nu,ka)*SY(ka,mu)       ) * &
                      &  (lrGammaOrb(mu,ka) + lrGammaOrb(nu,ka)) 
              enddo  
              tmprs2 = tmprs2 + dsmndr * tmprs 
            end do
          end do

          excgrad(xyz,alpha) = excgrad(xyz,alpha) + tmp1 + tmp2 + tmp4 + tmp6 + tmp3 - 0.25_dp * tmprs2
          excgrad(xyz,beta)  = excgrad(xyz,beta)  - tmp1 - tmp2 - tmp4 - tmp6 - tmp3 + 0.25_dp * tmprs2
        end do
      end do
    end do

    deallocate(shex)
    deallocate(xpyq)
    deallocate(shxpyq)
    deallocate(wcc)
    deallocate(qij)
    deallocate(temp)
    deallocate(dmn)

    deallocate(xpycc)
    deallocate(xmycc)
    deallocate(xpyas)
    deallocate(xmyas)
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

  end subroutine addGradients_rs


  ! Calculate oscillator strength for a given excitation range sep version
  subroutine getOscillatorStrengths_rs(sym, snglPartTransDip, eval, xpy, occNr,&
      & istat, oscStrength, tTradip, transDip)
    implicit none
    character, intent(in) :: sym
    real(dp), intent(in) :: snglPartTransDip(:,:), eval(:), xpy(:,:)
    real(dp), intent(in) :: occNr(:,:)
    logical :: tTradip
    integer, intent(inout) :: istat
    real(dp), intent(out) :: oscStrength(:), transDip(:,:)

    integer :: ii, nmat
    real(dp) :: oscStrength0
    logical :: spin
    
    nmat = size(xpy, dim=1)
    spin = .false. 
    if (size(occNr, dim=2) == 2) spin = .true. 

    ! Triplet oscillator strength and transition dipole is zero
    if ((.not. spin) .and. (sym == "T")) then
      oscStrength = 0.0_dp
      if (tTradip) transDip(:,:) = 0.0_dp
      return
    end if

    do ii = 1, size(xpy, dim=2)
      oscStrength(ii) = oscillatorStrength_rs(snglPartTransDip, sqrt(eval(ii)), xpy(:,ii))
    end do
    
    if (tTradip) call transitionDipole_rs(snglPartTransDip, xpy, transDip)

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
    
    pure function oscillatorStrength_rs(snglPartTransDip, omega, xpy) result(oscStrength)
     implicit none
      real(dp), intent(in) :: snglPartTransDip(:,:), omega, xpy(:)
      real(dp) :: oscStrength
      
      real(dp) :: rtmp
      integer :: ii
      
      oscStrength = 0.0_dp
      do ii = 1, 3
        rtmp = sum(snglPartTransDip(:,ii) * xpy)
        oscStrength = oscStrength + rtmp * rtmp
      end do
      oscStrength = twothird * 2.0_dp * omega * oscStrength
      
    end function oscillatorStrength_rs
    
    pure subroutine transitionDipole_rs(snglPartTransDip, xpy, transDip) 
      implicit none
      real(dp), intent(in) :: snglPartTransDip(:,:), xpy(:,:)
      real(dp), intent(out) :: transDip(:,:)
      
      integer :: ii, ll
      
      transDip(:,:) = 0.0_dp
      do ii = 1, size(xpy, dim=2)
        do ll = 1, 3
          transDip(ii,ll) = sum(snglPartTransDip(:,ll) * sqrt(2.0_dp) * xpy(:,ii))
        enddo
      end do
      
    end subroutine transitionDipole_rs
    
  end subroutine getOscillatorStrengths_rs


  !> Calculate transition charges for a given excitation
  !> (assumes singlet excitation for now, should set result to 0 for triplets)
  subroutine transitionCharges_rs(xpy, tqov, atomicTransQ)
    real(dp), intent(in) :: xpy(:), tqov(:,:)
    real(dp), intent(out) :: atomicTransQ(:)
    real(dp) :: prefactor 

    prefactor = sqrt(2.0_dp) !Check sqrt(2.0) !-> should be fine
    atomicTransQ(:) = prefactor * matmul(tqov, xpy) 

  end subroutine transitionCharges_rs


!subroutine getSqrS(coord, natom, skOverCont, interSKOver, orb, iAtomStart, species0, S)
subroutine getSqrS(coord, natom, skOverCont, orb, iAtomStart, species0, S)
  real(dp), intent(in) :: coord(:,:)
  integer,intent(in) :: natom, iAtomStart(:), species0(:)
  type(OSlakoCont), intent(in) :: skOverCont
  type(TOrbitals), intent(in) :: orb
  real(dp), intent(out) :: S(:,:)

 !real(dp) :: interSKOver(:)
  real(dp) :: SBlock(9,9)
  integer :: iAt1, iAt2, mu, nu, m, n
  
  S(:,:) = 0.0_dp

  do iAt1 = 1, natom
     do iAt2 = 1, iAt1-1  
           
        call getSOffsite(coord(:,iAt1), coord(:,iAt2), species0(iAt1), species0(iAt2),&
            & orb, skOverCont, SBlock)
           !& orb, skOverCont, SBlock, interSKOver)
        
        do mu = iAtomStart(iAt1), iAtomStart(iAt1+1)-1
           m = mu - iAtomStart(iAt1) + 1
           do nu = iAtomStart(iAt2), iAtomStart(iAt2+1)-1
              n = nu - iAtomStart(iAt2) + 1
              S(mu,nu) = SBlock(n,m) 
              S(nu,mu) = S(mu,nu)
           end do
        end do

     end do
  end do

  do mu = 1, size(S,dim=1)
     S(mu,mu) = 1.0_dp !Diagonal entries
  end do

end subroutine getSqrS


subroutine getSqrGamma(natom, lrGamma, iAtomStart, lrGammaOrb)
  implicit none
  real(dp), intent(in) :: lrGamma(:,:)
  integer,intent(in) :: natom, iAtomStart(:)
!Hblock, SBlock pass au, bu
!Hblock1, SBlock1 pass auh, buh
  real(dp), intent(out) :: lrGammaOrb(:,:)
  integer :: at1, at2, mu, nu, indat1, indat1p1, indat2, indat2p1
  
  lrGammaOrb(:,:) = 0.0_dp

  do at1 = 1, natom
     ! TOMAS KUBAR MODIFIED
     !indat1 = iAtomStart(at1) + 1
     !indat1p1 =iAtomStart(at1+1)
     indat1 = iAtomStart(at1)
     indat1p1 =iAtomStart(at1+1) - 1
     do at2 = 1, at1
        ! TOMAS KUBAR MODIFIED
        !indat2 = iAtomStart(at2) + 1
        !indat2p1 = iAtomStart(at2+1)
        indat2 = iAtomStart(at2)
        indat2p1 = iAtomStart(at2+1) - 1
        do mu = indat1, indat1p1
           do nu = indat2, indat2p1
              lrGammaOrb(mu,nu) = lrGamma(at1,at2) 
              lrGammaOrb(nu,mu) = lrGammaOrb(mu,nu)
           end do
        end do
     end do
  end do

end subroutine getSqrGamma

subroutine GetHvvXY(ipm, nxov, nxvv, norb, homo, natom, nMatUp, iatrans, getij, win, &
          & iAtomStart, sTimesGrndEigVecs, grndEigVecs, lrGamma, XorY, Hvv)
implicit none
real(dp), intent(in)  :: sTimesGrndEigVecs(:,:,:), grndEigVecs(:,:,:), lrGamma(:,:), XorY(:)
real(dp), intent(out) :: Hvv(:) 
real(dp), allocatable :: qij(:), gqij(:), qX(:,:), Gq(:,:)
integer, intent(in)   :: ipm, nxov, nxvv, norb, homo, natom, nMatUp, win(:), iAtomStart(:)
integer, intent(in)   :: iatrans(1:,homo+1:), getij(:,:)

integer               :: i, a, b, ia, ib, ab, nxov_act, nocc, nvirt
logical               :: lUpDwn

nvirt = 0
do while (nvirt*(nvirt+1)/2 < nxvv)
  nvirt = nvirt + 1
end do
nocc = norb - nvirt

nxov_act = nocc * nvirt

allocate(qij(natom))
allocate(gqij(natom))
allocate(qX(natom, nxov_act)) ! size(XorY))) ! allocate(qX(natom, nxov))
allocate(Gq(natom, nxov_act)) ! size(XorY))) ! allocate(Gq(natom, nxov))

qX(:,:) = 0.0_dp
do ia = 1, nxov_act ! size(XorY) ! nxov
   call indxov(win, ia, getij, i, a)
   lUpDwn = (win(ia) <= nMatUp)
   do b = homo + 1, norb
      ib = iatrans(i, b)
      if (ib <= nxov) then ! TOMAS KUBAR ADDED
         lUpDwn = .true.
         qij = transq(a, b, iAtomStart, lUpDwn, sTimesGrndEigVecs, grndEigVecs)
         qX(:,ia) = qX(:,ia) + qij(:)*XorY(ib)
      end if
   enddo
enddo

Gq(:,:)  = 0.0_dp
do ia = 1, nxov_act ! size(XorY)
   call indxov(win, ia, getij, i, a)
   lUpDwn = (win(ia) <= nMatUp)
   qij = transq(i, a, iAtomStart, lUpDwn, sTimesGrndEigVecs, grndEigVecs)
   call dsymv('U', natom, 1.0_dp, lrGamma, natom, qij, 1, 0.0_dp, gqij, 1)
   Gq(:,ia) = gqij(:)
enddo   

Hvv(:) = 0.0_dp
do ab = 1, nxvv
   call indxvv(homo, ab, a, b)
   do i = 1, homo
      ia = iatrans(i, a)
      ib = iatrans(i, b)
     !if ((ia <= nxov) .and. (ib <= nxov)) then ! TOMAS KUBAR ADDED
         Hvv(ab) = Hvv(ab) - ipm * ( dot_product(qX(:,ia), Gq(:,ib)) + ipm * dot_product(Gq(:,ia), qX(:,ib)) )  
     !end if
   enddo
enddo

deallocate(qij)
deallocate(gqij)
deallocate(qX)
deallocate(Gq)

end subroutine GetHvvXY

subroutine GetHooXY(ipm,& !nocc,
          & nxov, nxoo, norb, homo, natom, nMatUp, iatrans, getij, win, &
          & iAtomStart, sTimesGrndEigVecs, grndEigVecs, lrGamma, XorY, Hoo)
implicit none
real(dp), intent(in)  :: sTimesGrndEigVecs(:,:,:), grndEigVecs(:,:,:), lrGamma(:,:), XorY(:)
integer, intent(in)   :: ipm, nxov, nxoo, norb, homo, natom, nMatUp, win(:), iAtomStart(:) !, nocc
integer, intent(in)   :: iatrans(1:,homo+1:), getij(:,:)
real(dp), intent(out) :: Hoo(:) 

real(dp), allocatable :: qij(:), gqij(:), qX(:,:), Gq(:,:)
integer               :: i, a, j, ia, ja, ij, nxov_act, nocc, nvirt
logical               :: lUpDwn

nocc = 0
do while (nocc*(nocc+1)/2 < nxoo)
  nocc = nocc + 1
end do
nvirt = norb - nocc

nxov_act = nocc * nvirt

allocate(qij(natom))
allocate(gqij(natom))
allocate(qX(natom, nxov_act)) ! size(XorY))) ! allocate(qX(natom, nxov))
allocate(Gq(natom, nxov_act)) ! size(XorY))) ! allocate(Gq(natom, nxov))

qX(:,:) = 0.0_dp
do ia = 1, nxov_act ! size(XorY) ! nxov
   call indxov(win, ia, getij, i, a)
   lUpDwn = (win(ia) <= nMatUp)
   do j = 1, homo
      ja = iatrans(j, a)
      if (ja <= nxov) then ! TOMAS KUBAR ADDED
         lUpDwn = .true.
         qij = transq(i, j, iAtomStart, lUpDwn, sTimesGrndEigVecs, grndEigVecs)
         qX(:,ia) = qX(:,ia) + qij(:)*XorY(ja)
      end if
   enddo
enddo

Gq(:,:)  = 0.0_dp
do ia = 1, nxov_act ! size(XorY) ! nxov
   call indxov(win, ia, getij, i, a)
   lUpDwn = (win(ia) <= nMatUp)
   qij = transq(i, a, iAtomStart, lUpDwn, sTimesGrndEigVecs, grndEigVecs)
   call dsymv('U', natom, 1.0_dp, lrGamma, natom, qij, 1, 0.0_dp, gqij, 1)
   Gq(:,ia) = gqij(:)
enddo   

Hoo(:) = 0.0_dp
do ij = 1, nxoo
   call indxoo(ij, i, j)
   do a = homo + 1, norb
      ia = iatrans(i, a)
      ja = iatrans(j, a)
     !if ((ia <= nxov) .and. (ja <= nxov)) then ! TOMAS KUBAR ADDED
         Hoo(ij) = Hoo(ij) - ipm * ( dot_product(qX(:,ia), Gq(:,ja)) + ipm * dot_product(Gq(:,ia), qX(:,ja)) )  
     !end if
   enddo
enddo

deallocate(qij)
deallocate(gqij)
deallocate(qX)
deallocate(Gq)

end subroutine GetHooXY

subroutine GetHovT(nocc, nxov, nxoo, nxvv, norb, homo, natom, nMatUp, iatrans, getij, win, &
          & iAtomStart, sTimesGrndEigVecs, grndEigVecs, lrGamma, t, Hov)
implicit none
real(dp), intent(in)  :: sTimesGrndEigVecs(:,:,:), grndEigVecs(:,:,:), lrGamma(:,:), t(:,:)
integer, intent(in)   :: nocc, nxov, nxoo, nxvv, norb, homo, natom, nMatUp, win(:), iAtomStart(:)
integer, intent(in)   :: iatrans(1:,homo+1:), getij(:,:)
real(dp), intent(out) :: Hov(:) 

real(dp), allocatable :: qij(:), gqij(:), qX(:,:), Gq(:,:)
integer               :: i, a, b, j, ia, ib, ab, ja, ij
logical               :: lUpDwn

allocate(qij(natom))
allocate(gqij(natom))
allocate(qX(natom, nocc*(norb-nocc))) ! allocate(qX(natom, nxov))
allocate(Gq(natom, max(nxoo,nxvv)))

qX(:,:) = 0.0_dp
do ia = 1, nocc*(norb-nocc) ! nxov
   call indxov(win, ia, getij, i, a)
   lUpDwn = (win(ia) <= nMatUp)
   do b = homo + 1, norb
      ! TOMAS KUBAR -- SHOULD THE FOLLOWING BE if (iatrans(i, b) <= nxov) then ???
     !if (iatrans(i,b) <= nxov) then
         lUpDwn = (win(iatrans(i, b)) <= nMatUp) ! UNTESTED
         qij = transq(i, b, iAtomStart, lUpDwn, sTimesGrndEigVecs, grndEigVecs)
         qX(:,ia) = qX(:,ia) + qij(:)*t(a,b)
     !end if
   enddo
enddo

Gq(:,:)  = 0.0_dp
do ab = 1, nxvv
   call indxvv(homo, ab, a, b) 
   lUpDwn = .true.
   qij = transq(a, b, iAtomStart, lUpDwn, sTimesGrndEigVecs, grndEigVecs)
   call dsymv('U', natom, 1.0_dp, lrGamma, natom, qij, 1, 0.0_dp, gqij, 1)
   Gq(:,ab) = gqij(:)
enddo   

Hov(:) = 0.0_dp
do ia = 1, nxov ! nocc*(norb-nocc)
   call indxov(win, ia, getij, i, a)
   lUpDwn = (win(ia) <= nMatUp)
   do b = homo + 1, norb
      ib = iatrans(i, b)
     !if (ib <= nxov) then ! TOMAS KUBAR ADDED
         ! TOMAS KUBAR -- CORRECTED:
         !call rindxvv(nocc, a, b, ab)
         if (b < a) then
            call rindxvv(homo, a, b, ab)
         else
            call rindxvv(homo, b, a, ab)
         end if
         Hov(ia) = Hov(ia) - 2.0_dp * dot_product(qX(:,ib), Gq(:,ab)) 
     !end if
   enddo
enddo

qX(:,:) = 0.0_dp
do ia = 1, nocc*(norb-nocc) ! nxov
   call indxov(win, ia, getij, i, a)
   lUpDwn = (win(ia) <= nMatUp)
   do j = 1, homo
      ! TOMAS KUBAR -- SHOULD THE FOLLOWING BE if (iatrans(j, a) <= nxov) then ???
     !if (iatrans(j,a) <= nxov) then
         lUpDwn = (win(iatrans(j, a)) <= nMatUp) ! UNTESTED
         qij = transq(j, a, iAtomStart, lUpDwn, sTimesGrndEigVecs, grndEigVecs)
         qX(:,ia) = qX(:,ia) + qij(:)*t(i,j)
     !end if
   enddo
enddo

Gq(:,:)  = 0.0_dp
do ij = 1, nxoo
   call indxoo(ij, i, j)
   lUpDwn = .true.
   qij = transq(i, j, iAtomStart, lUpDwn, sTimesGrndEigVecs, grndEigVecs)
   call dsymv('U', natom, 1.0_dp, lrGamma, natom, qij, 1, 0.0_dp, gqij, 1)
   Gq(:,ij) = gqij(:)
enddo   

do ia = 1, nxov ! nocc*(norb-nocc)
   call indxov(win, ia, getij, i, a)
   lUpDwn = (win(ia) <= nMatUp)
   do j = 1, homo
      ja = iatrans(j, a)
     !if (ja <= nxov) then ! TOMAS KUBAR ADDED
        !call rindxoo(homo, nocc, i, j, ij)
         if (j < i) then
            call rindxoo(homo, nocc, i, j, ij)
         else
            call rindxoo(homo, nocc, j, i, ij)
         end if
         Hov(ia) = Hov(ia) - 2.0_dp * dot_product(qX(:,ja), Gq(:,ij)) 
     !end if
   enddo
enddo

deallocate(qij)
deallocate(gqij)
deallocate(qX)
deallocate(Gq)

end subroutine GetHovT

subroutine GetHooT(nocc, nxov, nxoo, norb, homo, natom, nMatUp, iatrans, getij, win, &
          & iAtomStart, sTimesGrndEigVecs, grndEigVecs, lrGamma, t, Hoo)
implicit none
real(dp), intent(in)  :: sTimesGrndEigVecs(:,:,:), grndEigVecs(:,:,:), lrGamma(:,:), t(:,:)
integer, intent(in)   :: nocc, nxov, nxoo, norb, homo, natom, nMatUp, win(:), iAtomStart(:)
integer, intent(in)   :: iatrans(1:,homo+1:), getij(:,:)
real(dp), intent(out) :: Hoo(:) 

real(dp), allocatable :: qij(:), gqij(:), qX(:,:), qXa(:,:,:), Gq(:,:)
integer               :: i, a, b, j, k, ia, ja, ij, jk, nxov_act !, ik
logical               :: lUpDwn

nxov_act = nocc*(norb-nocc)

allocate(qij(natom))
allocate(gqij(natom))
allocate(qX(natom, nxov_act)) ! nxov))
allocate(Gq(natom, max(nxoo, nxov_act))) ! nxov)))

qX(:,:) = 0.0_dp
do ia = 1, nxov_act ! nxov
   call indxov(win, ia, getij, i, a)
   lUpDwn = (win(ia) <= nMatUp)
   do b = homo + 1, norb
      ! TOMAS KUBAR -- SHOULD THE FOLLOWING BE if (iatrans(i, b) <= nxov) then ???
     !if (iatrans(i,b) <= nxov) then
         lUpDwn = (win(iatrans(i, b)) <= nMatUp) ! UNTESTED
         qij = transq(i, b, iAtomStart, lUpDwn, sTimesGrndEigVecs, grndEigVecs)
         qX(:,ia) = qX(:,ia) + qij(:)*t(a,b)
     !end if
   enddo
enddo

Gq(:,:)  = 0.0_dp
do ia = 1, nxov_act ! nxov
   call indxov(win, ia, getij, i, a)
   lUpDwn = (win(ia) <= nMatUp)
   qij = transq(i, a, iAtomStart, lUpDwn, sTimesGrndEigVecs, grndEigVecs)
   call dsymv('U', natom, 1.0_dp, lrGamma, natom, qij, 1, 0.0_dp, gqij, 1)
   Gq(:,ia) = gqij(:)
enddo   

Hoo(:) = 0.0_dp
do ij = 1, nxoo
   call indxoo(ij, i, j)
   do a = homo + 1, norb
      ia = iatrans(i, a)
      ja = iatrans(j, a)
     !if (ia <= nxov .and. ja <= nxov) then ! TOMAS KUBAR ADDED
         Hoo(ij) = Hoo(ij) - 2.0_dp * dot_product(qX(:,ia), Gq(:,ja)) 
     !end if
   enddo
enddo
deallocate(qX)

allocate(qXa(natom, nocc, nocc))
qXa(:,:,:) = 0.0_dp
do i = 1, homo
   do k = 1, homo
      lUpDwn = .true.
      qij = transq(i, k, iAtomStart, lUpDwn, sTimesGrndEigVecs, grndEigVecs)
      do j = 1, homo
         qXa(:,i,j) = qXa(:,i,j) + qij(:)*t(j,k)
      enddo
   enddo
enddo

Gq(:,:)  = 0.0_dp
do ij = 1, nxoo
   call indxoo(ij, i, j)
   lUpDwn = .true.
   qij = transq(i, j, iAtomStart, lUpDwn, sTimesGrndEigVecs, grndEigVecs)
   call dsymv('U', natom, 1.0_dp, lrGamma, natom, qij, 1, 0.0_dp, gqij, 1)
   Gq(:,ij) = gqij(:)
enddo   

do ij = 1, nxoo
   call indxoo(ij, i, j)
   do k = 1, homo
     !call rindxoo(homo, nocc, i, k, ik) -- NOT NEEDED
     !call rindxoo(homo, nocc, j, k, jk)
      if (k < j) then
         call rindxoo(homo, nocc, j, k, jk)
      else
         call rindxoo(homo, nocc, k, j, jk)
      end if
      Hoo(ij) = Hoo(ij) - 2.0_dp * dot_product(qXa(:,i,k), Gq(:,jk))
   enddo
enddo

deallocate(qij)
deallocate(gqij)
deallocate(qXa)
deallocate(Gq)

end subroutine GetHooT

end module dftbp_rs_linearresponse
