!--------------------------------------------------------------------------------------------------!
!  DFTB+: general package for performing fast atomistic simulations                                !
!  Copyright (C) 2006 - 2023  DFTB+ developers group                                               !
!                                                                                                  !
!  See the LICENSE file for terms of usage and distribution.                                       !
!--------------------------------------------------------------------------------------------------!

#:include 'common.fypp'

!> Helper routines for transition charges between levels.
module dftbp_timedep_transcharges
  use mpi
  use dftbp_common_accuracy, only : dp
  use dftbp_type_densedescr, only : TDenseDescr
  use dftbp_common_environment, only : TEnvironment
  use dftbp_extlibs_scalapackfx, only : DLEN_, NB_, CSRC_, MB_, RSRC_, scalafx_indxl2g
  use dftbp_extlibs_mpifx, only : mpifx_allreduceip
  use dftbp_io_message, only : error
  use dftbp_extlibs_scalapackfx, only : pblasfx_pgemm
  implicit none

  private
  public :: TTransCharges, TTransCharges_init
  public :: transq

  !> Internal data of the transition charges
  type :: TTransCharges
    private

    !> should transition charges be cached in memory or evaluated when needed?
    logical :: tCacheChargesOccVir

    !> same for occ-occ/vir-vir transitions
    logical :: tCacheChargesSame

    !> storage if caching the occupied -> virtual transition charges
    real(dp), allocatable :: qCacheOccVir(:,:)

    !> storage if caching the occupied -> occupied transition charges
    real(dp), allocatable :: qCacheOccOcc(:,:)

    !> storage if caching the virtual -> virtual transition charges
    real(dp), allocatable :: qCacheVirVir(:,:)

    !> Number of occ-vir transitions in the cache
    integer :: nTransitions

    !> Number of atoms in the system
    integer :: nAtom

    !> Number of occ-vir transitions within the spin-up block
    integer :: nMatUp

    !> Number of occ-occ transitions within the spin-up block
    integer :: nMatUpOccOcc

    !> Number of vir-vir transitions within the spin-up block
    integer :: nMatUpVirVir

  contains

    procedure :: qTransIA =>TTransCharges_qTransIA
    procedure :: qTransIJ =>TTransCharges_qTransIJ
    procedure :: qTransAB =>TTransCharges_qTransAB
    procedure :: qMatVec => TTransCharges_qMatVec
    procedure :: qVecMat => TTransCharges_qVecMat
    procedure :: qMatVecDs => TTransCharges_qMatVecDs
    procedure :: qVecMatDs => TTransCharges_qVecMatDs

  end type TTransCharges


contains

  !> initialise the cache/on-the fly transition charge evaluator
  subroutine TTransCharges_init(this, env, denseDesc, sTimesGrndEigVecs, grndEigVecs, nTrans, nMatUp,&
      & nXooUD, nXvvUD, getia, getij, getab, iatrans, win, tStoreOccVir, tStoreSame)

    !> Instance
    type(TTransCharges), intent(out) :: this

    !> Environment settings
    type(TEnvironment), intent(in) :: env

    !> Dense matrix descriptor
    type(TDenseDescr), intent(in) :: denseDesc 

    !> Overlap times eigenvector: sum_m Smn cmi (nOrb, nOrb)
    real(dp), intent(in) :: sTimesGrndEigVecs(:,:,:)

    !> Eigenvectors (nOrb, nOrb)
    real(dp), intent(in) :: grndEigVecs(:,:,:)

    !> number of transitions in the system
    integer, intent(in) :: nTrans

    !> number of up-up occ-vir excitations
    integer, intent(in) :: nMatUp

    !> number of occ-occ excitations per spin channel
    integer, intent(in) :: nXooUD(:)

    !> number of vir-vir excitations per spin channel
    integer, intent(in) :: nXvvUD(:)

    !> should occ-vir transitions be stored?
    !> this is sufficient for single-point TD-DFTB
    logical, intent(in) :: tStoreOccVir

    !> should also occ-occ and vir-vir transitions be stored?
    !> required for excited state forces and TD-LC-DFTB
    logical, intent(in) :: tStoreSame

    !> index array for occ-vir single particle excitations
    integer, intent(in) :: getia(:,:)

    !> index array for occ-occ single particle excitations
    integer, intent(in) :: getij(:,:)

    !> index array for vir-vir single particle excitations
    integer, intent(in) :: getab(:,:)
    
    !> array from pairs of single particles states to compound index 
    integer, intent(in) :: iatrans(:,:,:)

    !> index array for single particle excitations that are included
    integer, intent(in) :: win(:)

    integer :: ia, ij, ii, jj, kk, ab, aa, bb, iam, nProc, iSpin, iAtom, nSpin
    integer :: ss, mm, desc(DLEN_), iLoc, aLoc, mLoc, iGlb, aGlb, mGlb, jLoc, jGlb, nOrbs
    integer :: nLocRow, nLocCol
     real(dp) :: rSum
    real(dp), allocatable :: maskMat(:,:), workLocal(:,:), workGlobal(:,:,:)
    logical :: updwn

    this%nTransitions = nTrans
    this%nAtom = size(denseDesc%iAtomStart) - 1
    this%nMatUp = nMatUp
    this%nMatUpOccOcc = nXooUD(1)
    this%nMatUpVirVir = nXvvUD(1)
    nSpin = size(grndEigVecs, dim=3)
    nOrbs = size(iatrans, dim=1)

    if (tStoreOccVir) then

      @:ASSERT(.not.allocated(this%qCacheOccVir))
      allocate(this%qCacheOccVir(this%nAtom, nTrans))

    #:if WITH_SCALAPACK
      call blacs_pinfo(iam, nProc)
      desc(:) = denseDesc%blacsOrbSqr
      allocate(workGlobal(nOrbs,nOrbs,nSpin))
      nLocRow = size(grndEigVecs,dim=1)
      nLocCol = size(grndEigVecs,dim=2)
      allocate(workLocal(nLocRow,nLocCol))
      allocate(maskMat(nLocRow,nLocCol))

      do iSpin = 1, nSpin
        do iAtom = 1, this%nAtom
          workGlobal = 0.0_dp
          workLocal = 0.0_dp
          maskMat = 0.0_dp
          do mLoc = 1, size(grndEigVecs, dim=1)
            mGlb = scalafx_indxl2g(mLoc, desc(MB_), env%blacs%orbitalGrid%myrow, desc(RSRC_), &
                & env%blacs%orbitalGrid%nrow)
            if (mGlb >= denseDesc%iAtomStart(iAtom) .and. mGlb < denseDesc%iAtomStart(iAtom + 1)) then
              maskMat(mLoc,:) = sTimesGrndEigVecs(mLoc,:,iSpin)
            end if
          end do
          call pblasfx_pgemm(grndEigVecs(:,:,iSpin), denseDesc%blacsOrbSqr, maskMat, &
              & denseDesc%blacsOrbSqr, workLocal, denseDesc%blacsOrbSqr, transa="T")
          do jLoc = 1, size(grndEigVecs, dim=2)
             jGlb = scalafx_indxl2g(jLoc, desc(NB_), env%blacs%orbitalGrid%mycol, desc(CSRC_), &
                 & env%blacs%orbitalGrid%ncol)
            do iLoc = 1, size(grndEigVecs, dim=1)
             iGlb = scalafx_indxl2g(iLoc, desc(MB_), env%blacs%orbitalGrid%myrow, desc(RSRC_), &
                 & env%blacs%orbitalGrid%nrow)
              workGlobal(iGlb,jGlb,iSpin) =  workLocal(iLoc,jLoc)
            enddo
          enddo
          call mpifx_allreduceip(env%mpi%groupComm, workGlobal, MPI_SUM)
          do ia = 1, nTrans
            kk = win(ia)
            ii = getia(kk,1)
            aa = getia(kk,2)
            if(kk <= nMatUp) then
              this%qCacheOccVir(iAtom, ia) = 0.5_dp * (workGlobal(ii,aa,1) + workGlobal(aa,ii,1))
            else
              this%qCacheOccVir(iAtom, ia) = 0.5_dp * (workGlobal(ii,aa,2) + workGlobal(aa,ii,2))
            end if
          enddo
        enddo
      enddo
                       
!!$      this%qCacheOccVir = 0.0_dp
!!$      do ss = 1, size(grndEigVecs, dim=3)
!!$        do iLoc = 1, size(grndEigVecs, dim=2)
!!$          iGlb = scalafx_indxl2g(iLoc, desc(NB_), env%blacs%orbitalGrid%mycol, desc(CSRC_), &
!!$              & env%blacs%orbitalGrid%ncol)
!!$          do aLoc = 1, size(grndEigVecs, dim=2)
!!$            aGlb = scalafx_indxl2g(aLoc, desc(NB_), env%blacs%orbitalGrid%mycol, desc(CSRC_), &
!!$                & env%blacs%orbitalGrid%ncol)
!!$            iasGlb = iaTrans(iGlb, aGlb, ss)
!!$            print *,'the prev',iGlb,aGlb,iasGlb
!!$            qTmp(:) =  0.5_dp * (grndEigVecs(:,iLoc,ss) * sTimesGrndEigVecs(:,aLoc,ss)&
!!$                & + grndEigVecs(:,aLoc,ss) * sTimesGrndEigVecs(:,iLoc,ss))
!!$            do kk = 1, this%nAtom
!!$              rSum = 0.0_dp
!!$              do mLoc = 1, size(grndEigVecs, dim=1)
!!$                mGlb = scalafx_indxl2g(mLoc, desc(MB_), env%blacs%orbitalGrid%myrow, desc(RSRC_), &
!!$                    & env%blacs%orbitalGrid%nrow)
!!$                if(iasGlb == 1081) then
!!$                  print *,'indixes',mGlb,mLoc,denseDesc%iAtomStart(kk),denseDesc%iAtomStart(kk + 1)
!!$                endif
!!$                if (mGlb >= denseDesc%iAtomStart(kk) .and. mGlb < denseDesc%iAtomStart(kk + 1)) then
!!$                  rSum = rSum + qTmp(mLoc)
!!$                end if
!!$              end do
!!$              this%qCacheOccVir(kk, iasGlb) = rSum
!!$            end do
!!$          end do
!!$        end do
!!$      end do
!!$      call mpifx_allreduceip(env%mpi%groupComm, this%qCacheOccVir, MPI_SUM)
!!$      print *,'afterwards'
          
      ! allocate(aIbJ(desc(M_), nTrans), source=1.0_dp)
      ! allocate(aJbI(desc(M_), nTrans), source=1.0_dp)

      ! do ia = 1, nTrans
      !   kk = win(ia)
      !   ii = getia(kk,1)
      !   aa = getia(kk,2)
        
      !   if(kk <= this%nMatUp) then
      !     ss = 1
      !   else
      !     ss = 2
      !   endif
          
      !   do mm = 1, desc(M_)
      !     call scalafx_infog2l(env%blacs%orbitalGrid, desc, mm, ii, mLoc, iLoc, rSrc, cSrc)
      !     if(env%blacs%orbitalGrid%myrow == rSrc .and. env%blacs%orbitalGrid%mycol == cSrc) then
      !       aIbJ(mm, ia) = aIbJ(mm, ia) * grndEigVecs(mLoc, iLoc, ss)
      !       aJbI(mm, ia) = aJbI(mm, ia) * sTimesGrndEigVecs(mLoc, iLoc, ss)
      !     end if
      !     call scalafx_infog2l(env%blacs%orbitalGrid, desc, mm, aa, mLoc, aLoc, rSrc, cSrc)
      !     if(env%blacs%orbitalGrid%myrow == rSrc .and. env%blacs%orbitalGrid%mycol == cSrc) then
      !       aIbJ(mm, ia) = aIbJ(mm, ia) * sTimesGrndEigVecs(mLoc,aLoc,ss)
      !       aJbI(mm, ia) = aJbI(mm, ia) * grndEigVecs(mLoc, aLoc, ss)
      !     endif
      !   end do
      ! end do
   
      ! ! I think, this shold be actually the group communicator mpiGroupComm
      ! ! For the moment, we should ensure, that only one group is used, and later check carefully...
      ! call mpifx_allreduceip(env%mpi%globalComm, aIbJ, MPI_PROD)
      ! call mpifx_allreduceip(env%mpi%globalComm, aJbI, MPI_PROD)
      ! do ia = 1, nTrans
      !   do kk = 1, this%nAtom
      !     aa = denseDesc%iAtomStart(kk)
      !     bb = denseDesc%iAtomStart(kk + 1) -1
      !     this%qCacheOccVir(kk, ia) =  0.5_dp * sum(aIbJ(aa:bb, ia) + aJbI(aa:bb, ia))      
      !   end do
      ! end do

   #:else

      !!$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(ia,ii,aa,kk,updwn) SCHEDULE(RUNTIME)
      do ia = 1, nTrans
        kk = win(ia)
        ii = getia(kk,1)
        aa = getia(kk,2)
        updwn = (kk <= this%nMatUp)
        this%qCacheOccVir(:,ia) = transq(ii, aa, env, denseDesc, updwn,  sTimesGrndEigVecs,&
            & grndEigVecs)
      end do
      !!$OMP  END PARALLEL DO

   #:endif
      
    else

      this%tCacheChargesOccVir = .false.

    end if

    if (tStoreSame) then
      @:ASSERT(.not.allocated(this%qCacheOccOcc))
      allocate(this%qCacheOccOcc(this%nAtom, sum(nXooUD)))
      @:ASSERT(.not.allocated(this%qCacheVirVir))
      allocate(this%qCacheVirVir(this%nAtom, sum(nXvvUD)))

      !!$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(ij,ii,jj,updwn) SCHEDULE(RUNTIME)
      do ij = 1, sum(nXooUD)
        ii = getij(ij,1)
        jj = getij(ij,2)
        updwn = (ij <= nXooUD(1))
        this%qCacheOccOcc(:,ij) = transq(ii, jj, env, denseDesc, updwn,  sTimesGrndEigVecs,&
             & grndEigVecs)
      enddo
      !!$OMP  END PARALLEL DO

      !!$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(ab,aa,bb,updwn) SCHEDULE(RUNTIME)
      do ab = 1, sum(nXvvUD)
        aa = getab(ab,1)
        bb = getab(ab,2)
        updwn = (ab <= nXvvUD(1))
        this%qCacheVirVir(:,ab) = transq(aa, bb, env, denseDesc, updwn,  sTimesGrndEigVecs,&
            & grndEigVecs)
      end do
      !!$OMP  END PARALLEL DO

      this%tCacheChargesSame = .true.

    else

      this%tCacheChargesSame = .false.

    end if

  end subroutine TTransCharges_init


  !> returns transition charges between occ-vir single particle levels
  function TTransCharges_qTransIA(this, ij, env, denseDesc, sTimesGrndEigVecs, grndEigVecs, getia,&
      & win) result(q)

    !> instance of the transition charge object
    class(TTransCharges), intent(in) :: this

    !> Index of transition
    integer, intent(in) :: ij

    !> Environment settings
    type(TEnvironment), intent(in) :: env

    !> Dense matrix descriptor
    type(TDenseDescr), intent(in) :: denseDesc 

    real(dp), intent(in) :: sTimesGrndEigVecs(:,:,:)

    !> Eigenvectors (nOrb, nOrb)
    real(dp), intent(in) :: grndEigVecs(:,:,:)

    !> index array for for single particle excitations
    integer, intent(in) :: getia(:,:)

    !> index array for single particle excitations that are included
    integer, intent(in) :: win(:)

    !> Transition charge on exit. (nAtom)
    real(dp), dimension(size(denseDesc%iAtomStart)-1) :: q

    logical :: updwn
    integer :: ii, jj, kk

    if (allocated(this%qCacheOccVir)) then
      q(:) = this%qCacheOccVir(:, ij)
    else
      kk = win(ij)
      ii = getia(kk,1)
      jj = getia(kk,2)
      updwn = (kk <= this%nMatUp)
      q(:) = transq(ii, jj, env, denseDesc, updwn, sTimesgrndEigVecs, grndEigVecs)
    end if

  end function TTransCharges_qTransIA

  !> returns transition charges between occ-occ single particle levels
  function TTransCharges_qTransIJ(this, ij, env, denseDesc, sTimesGrndEigVecs, grndEigVecs, getij)&
      & result(q)

    !> instance of the transition charge object
    class(TTransCharges), intent(in) :: this

    !> Index of transition
    integer, intent(in) :: ij

    !> Environment settings
    type(TEnvironment), intent(in) :: env

    !> Dense matrix descriptor
    type(TDenseDescr), intent(in) :: denseDesc 

    !> Overlap times eigenvector: sum_m Smn cmi (nOrb, nOrb)
    real(dp), intent(in) :: sTimesGrndEigVecs(:,:,:)

    !> Eigenvectors (nOrb, nOrb)
    real(dp), intent(in) :: grndEigVecs(:,:,:)

    !> index array for occ-occ single particle excitations
    integer, intent(in) :: getij(:,:)

    !> Transition charge on exit. (nAtom)
    real(dp), dimension(size(denseDesc%iAtomStart)-1) :: q

    logical :: updwn
    integer :: ii, jj

    if (allocated(this%qCacheOccOcc)) then
      q(:) = this%qCacheOccOcc(:, ij)
    else
      ii = getij(ij,1)
      jj = getij(ij,2)
      updwn = (ij <= this%nMatUpOccOcc)
      q(:) = transq(ii, jj, env, denseDesc, updwn, sTimesgrndEigVecs, grndEigVecs)
    end if

  end function TTransCharges_qTransIJ

  !> returns transition charges between vir-vir single particle levels
  function TTransCharges_qTransAB(this, ab, env, denseDesc, sTimesGrndEigVecs, grndEigVecs, getab)&
      & result(q)

    !> instance of the transition charge object
    class(TTransCharges), intent(in) :: this

    !> Index of transition
    integer, intent(in) :: ab

    !> Environment settings
    type(TEnvironment), intent(in) :: env

    !> Dense matrix descriptor
    type(TDenseDescr), intent(in) :: denseDesc 

    !> Overlap times eigenvector: sum_m Smn cmi (nOrb, nOrb)
    real(dp), intent(in) :: sTimesGrndEigVecs(:,:,:)

    !> Eigenvectors (nOrb, nOrb)
    real(dp), intent(in) :: grndEigVecs(:,:,:)

    !> index array for occ-occ single particle excitations
    integer, intent(in) :: getab(:,:)

    !> Transition charge on exit. (nAtom)
    real(dp), dimension(size(denseDesc%iAtomStart)-1) :: q

    logical :: updwn
    integer :: aa, bb

    if (allocated(this%qCacheVirVir)) then
      q(:) = this%qCacheVirVir(:, ab)
    else
      aa = getab(ab,1)
      bb = getab(ab,2)
      updwn = (ab <= this%nMatUpVirVir)
      q(:) = transq(aa, bb, env, denseDesc, updwn, sTimesgrndEigVecs, grndEigVecs)
    end if

  end function TTransCharges_qTransAB


  !> Transition charges left produced with a vector Q * v
  !> qProduct has dimension nAtom
  subroutine TTransCharges_qMatVec(this, env, denseDesc, sTimesGrndEigVecs, grndEigVecs, getia,&
      & win, vector, qProduct, indexOffSet)

    !> instance of the transition charge object
    class(TTransCharges), intent(in) :: this

    !> Environment settings
    type(TEnvironment), intent(in) :: env

    !> Dense matrix descriptor
    type(TDenseDescr), intent(in) :: denseDesc 

    !> Overlap times eigenvector: sum_m Smn cmi (nOrb, nOrb)
    real(dp), intent(in) :: sTimesGrndEigVecs(:,:,:)

    !> Eigenvectors (nOrb, nOrb)
    real(dp), intent(in) :: grndEigVecs(:,:,:)

    !> index array for for single particle excitations
    integer, intent(in) :: getia(:,:)

    !> index array for single particle excitations that are included
    integer, intent(in) :: win(:)

    !> vector to product with the transition charges
    real(dp), intent(in) :: vector(:)

    !> Product on exit
    real(dp), intent(inout) :: qProduct(:)

    !> Offset of vector index (optional) which determines orbital pair
    integer, intent(in), optional :: indexOffSet

    real(dp), allocatable :: qij(:)
    integer :: ii, jj, ij, kk, iOff, iGlb, fGlb
    logical :: updwn

    ! RPA vectors are distributed for MPI: size(vector) < nOcc*nVir
    if (present(indexOffSet)) then
      iOff = indexOffSet
    else
      iOff = 0
    end if
    iGlb = iOff + 1
    fGlb = iOff + size(vector)

    if (this%tCacheChargesOccVir) then

      qProduct(:) = qProduct + matmul(this%qCacheOccVir(:,iGlb:fGlb), vector)

    else
      
      allocate(qij(this%nAtom))

      do ij = 1, size(vector)
        kk = win(iOff+ij)
        ii = getia(kk,1)
        jj = getia(kk,2)
        
        updwn = (kk <= this%nMatUp)
        qij(:) = transq(ii, jj, env, denseDesc, updwn, sTimesGrndEigVecs, grndEigVecs)
        qProduct(:) = qProduct + qij * vector(ij)
      end do

      deallocate(qij)

    end if

  end subroutine TTransCharges_qMatVec


  !> Transition charges right produced with a vector v * Q
  !> qProduct has dimension nOcc*nVir
  subroutine TTransCharges_qVecMat(this, env, denseDesc, sTimesGrndEigVecs, grndEigVecs, getia,&
      & win, vector, qProduct, indexOffSet)

    !> instance of the transition charge object
    class(TTransCharges), intent(in) :: this

    !> Environment settings
    type(TEnvironment), intent(in) :: env

    !> Dense matrix descriptor
    type(TDenseDescr), intent(in) :: denseDesc 

    !> Overlap times eigenvector: sum_m Smn cmi (nOrb, nOrb)
    real(dp), intent(in) :: sTimesGrndEigVecs(:,:,:)

    !> Eigenvectors (nOrb, nOrb)
    real(dp), intent(in) :: grndEigVecs(:,:,:)

    !> index array for for single particle excitations
    integer, intent(in) :: getia(:,:)

    !> index array for single particle excitations that are included
    integer, intent(in) :: win(:)

    !> vector to product with the transition charges
    real(dp), intent(in) :: vector(:)

    !> Product on exit
    real(dp), intent(inout) :: qProduct(:)

    !> Offset of vector index (optional) which determines orbital pair
    integer, intent(in), optional :: indexOffSet

    real(dp), allocatable :: qij(:)
    integer :: ii, jj, ij, kk, iOff, iGlb, fGlb
    logical :: updwn

    ! RPA vectors are distributed for MPI: size(vector) < nOcc*nVir
    if (present(indexOffSet)) then
      iOff = indexOffSet
    else
      iOff = 0
    end if
    iGlb = iOff + 1
    fGlb = iOff + size(qProduct)
    
    if (this%tCacheChargesOccVir) then

      qProduct(:) = qProduct + matmul(vector, this%qCacheOccVir(:,iGlb:fGlb))

    else
      
      allocate(qij(this%nAtom))
      !!$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(ij,ii,jj,kk,updwn,qij)&
      !!$OMP& SCHEDULE(RUNTIME)
      do ij = 1, size(qProduct)
        kk = win(iOff+ij)
        ii = getia(kk,1)
        jj = getia(kk,2)
        updwn = (kk <= this%nMatUp)
        qij(:) = transq(ii, jj, env, denseDesc, updwn, sTimesGrndEigVecs, grndEigVecs)
        qProduct(ij) = qProduct(ij) + dot_product(qij, vector)
      end do
      !!$OMP  END PARALLEL DO

      deallocate(qij)

    end if

  end subroutine TTransCharges_qVecMat


  !> Transition charges left produced with a vector Q * v for spin up
  !> minus Transition charges left produced with a vector Q * v for spin down
  !> sum_ias q_ias V_ias delta_s,  where delta_s = 1 for spin up and delta_s = -1 for spin down
  subroutine TTransCharges_qMatVecDs(this, env, denseDesc, sTimesGrndEigVecs, grndEigVecs, getia,&
      & win, vector, qProduct)

    !> instance of the transition charge object
    class(TTransCharges), intent(in) :: this

    !> Environment settings
    type(TEnvironment), intent(in) :: env

    !> Dense matrix descriptor
    type(TDenseDescr), intent(in) :: denseDesc

    !> Overlap times eigenvector: sum_m Smn cmi (nOrb, nOrb)
    real(dp), intent(in) :: sTimesGrndEigVecs(:,:,:)

    !> Eigenvectors (nOrb, nOrb)
    real(dp), intent(in) :: grndEigVecs(:,:,:)

    !> index array for for single particle excitations
    integer, intent(in) :: getia(:,:)

    !> index array for single particle excitations that are included
    integer, intent(in) :: win(:)

    !> vector to product with the transition charges
    real(dp), intent(in) :: vector(:)

    !> Product on exit
    real(dp), intent(inout) :: qProduct(:)

    real(dp), allocatable :: qij(:)
    integer :: ii, jj, ij, kk
    logical :: updwn

    allocate(qij(this%nAtom))

    !!$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(ij,ii,jj,kk,updwn,qij)&
    !!$OMP& SCHEDULE(RUNTIME) REDUCTION(+:qProduct)
    do ij = 1, this%nTransitions
      kk = win(ij)
      ii = getia(kk,1)
      jj = getia(kk,2)
      updwn = (kk <= this%nMatUp)
      qij(:) = transq(ii, jj, env, denseDesc, updwn, sTimesGrndEigVecs, grndEigVecs)
      if (updwn) then
        qProduct(:) = qProduct + qij * vector(ij)
      else
        qProduct(:) = qProduct - qij * vector(ij)
      end if
    end do
    !!$OMP  END PARALLEL DO

    deallocate(qij)

  end subroutine TTransCharges_qMatVecDs


  !> Transition charges right produced with a vector v * Q for spin up
  !> and negative transition charges right produced with a vector v * Q for spin down
  !> R_ias = delta_s sum_A q_A^(ias) V_A,  where delta_s = 1 for spin up and delta_s = -1 for spin down
  subroutine TTransCharges_qVecMatDs(this, env, denseDesc, sTimesGrndEigVecs, grndEigVecs, getia,&
      & win, vector, qProduct)

    !> instance of the transition charge object
    class(TTransCharges), intent(in) :: this

    !> Environment settings
    type(TEnvironment), intent(in) :: env

    !> Dense matrix descriptor
    type(TDenseDescr), intent(in) :: denseDesc 

    !> Overlap times eigenvector: sum_m Smn cmi (nOrb, nOrb)
    real(dp), intent(in) :: sTimesGrndEigVecs(:,:,:)

    !> Eigenvectors (nOrb, nOrb)
    real(dp), intent(in) :: grndEigVecs(:,:,:)

    !> index array for for single particle excitations
    integer, intent(in) :: getia(:,:)

    !> index array for single particle excitations that are included
    integer, intent(in) :: win(:)

    !> vector to product with the transition charges
    real(dp), intent(in) :: vector(:)

    !> Product on exit
    real(dp), intent(inout) :: qProduct(:)

    real(dp), allocatable :: qij(:)
    integer :: ii, jj, ij, kk
    logical :: updwn

    allocate(qij(this%nAtom))

    !!$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(ij,ii,jj,kk,updwn,qij)&
    !!$OMP& SCHEDULE(RUNTIME)
    do ij = 1, this%nTransitions
      kk = win(ij)
      ii = getia(kk,1)
      jj = getia(kk,2)
      updwn = (kk <= this%nMatUp)
      qij(:) = transq(ii, jj, env, denseDesc, updwn, sTimesGrndEigVecs, grndEigVecs)
      if (updwn) then
        qProduct(ij) = qProduct(ij) + dot_product(qij, vector)
      else
        qProduct(ij) = qProduct(ij) - dot_product(qij, vector)
      end if
    end do
    !!$OMP  END PARALLEL DO

    deallocate(qij)

  end subroutine TTransCharges_qVecMatDs


  !> Calculates atomic transition charges for a specified excitation.
  !> Calculates qij = 0.5 * (c_i S c_j + c_j S c_i) where c_i and c_j are selected eigenvectors, and
  !> S the overlap matrix.
  !> Since qij is atomic quantity (so far) the corresponding values for the atom are summed up.
  !> Note: the parameters 'updwn' were added for spin alpha and beta channels.
  function transq(ii, jj, env, denseDesc, updwn, sTimesGrndEigVecs, grndEigVecs) result(qij)

    !> Index of initial state.
    integer, intent(in) :: ii

    !> Index of final state.
    integer, intent(in) :: jj

    !> Environment settings
    type(TEnvironment), intent(in) :: env

    !> Dense matrix descriptor
    type(TDenseDescr), intent(in) :: denseDesc 

    !> up spin channel (T) or down spin channel (F)
    logical, intent(in) :: updwn

    !> Overlap times eigenvector: sum_m Smn cmi (nOrb, nOrb).
    real(dp), intent(in) :: sTimesGrndEigVecs(:,:,:)

    !> Eigenvectors (nOrb, nOrb)
    real(dp), intent(in) :: grndEigVecs(:,:,:)

    !> Transition charge on exit. (nAtom)
    real(dp) :: qij(size(denseDesc%iAtomStart)-1)

    integer :: kk, aa, bb, ss, mm
    real(dp), allocatable :: qTmp(:)

    ss = 1
    if (.not. updwn) then
      ss = 2
    end if

  #:if WITH_SCALAPACK

    call error('Direct call to transq not possible with MPI.')
    
   ! Dimension of qTmp is nOrb also for MPI 
   !  allocate(qTmp(desc(M_)))
   !  allocate(aIbJ(desc(M_)), source=1.0_dp)
   !  allocate(aJbI(desc(M_)), source=1.0_dp)

   !  do mm = 1, desc(M_)
   !    call scalafx_infog2l(env%blacs%orbitalGrid, desc, mm, ii, mLoc, iLoc, rSrc, cSrc)
   !    if(env%blacs%orbitalGrid%myrow == rSrc .and. env%blacs%orbitalGrid%mycol == cSrc) then
   !      aIbJ(mm) = aIbJ(mm) * grndEigVecs(mLoc, iLoc, ss)
   !      aJbI(mm) = aJbI(mm) * sTimesGrndEigVecs(mLoc, iLoc, ss)
   !    end if
   !    call scalafx_infog2l(env%blacs%orbitalGrid, desc, mm, jj, mLoc, jLoc, rSrc, cSrc)
   !    if(env%blacs%orbitalGrid%myrow == rSrc .and. env%blacs%orbitalGrid%mycol == cSrc) then
   !      aIbJ(mm) = aIbJ(mm) * sTimesGrndEigVecs(mLoc,jLoc,ss)
   !      aJbI(mm) = aJbI(mm) * grndEigVecs(mLoc, jLoc, ss)
   !    endif
   ! end do
   
   ! ! I think, this shold be actually the group communicator mpiGroupComm
   ! ! For the moment, we should ensure, that only one group is used, and later check carefully...
   ! call mpifx_allreduceip(env%mpi%globalComm, aIbJ, MPI_PROD)
   ! call mpifx_allreduceip(env%mpi%globalComm, aJbI, MPI_PROD)
   ! qTmp = aIbJ + aJbI
   
   #:else
    
    allocate(qTmp, size(grndEigVecs,dim=1))
    qTmp(:) =  grndEigVecs(:,ii,ss) * sTimesGrndEigVecs(:,jj,ss)&
          & + grndEigVecs(:,jj,ss) * sTimesGrndEigVecs(:,ii,ss)

   #:endif
    
    do kk = 1, size(qij)
      aa = denseDesc%iAtomStart(kk)
      bb = denseDesc%iAtomStart(kk + 1) -1
      qij(kk) = 0.5_dp * sum(qTmp(aa:bb))
    end do

  end function transq


end module dftbp_timedep_transcharges
