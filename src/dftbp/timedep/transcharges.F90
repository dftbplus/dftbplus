!--------------------------------------------------------------------------------------------------!
!  DFTB+: general package for performing fast atomistic simulations                                !
!  Copyright (C) 2006 - 2023  DFTB+ developers group                                               !
!                                                                                                  !
!  See the LICENSE file for terms of usage and distribution.                                       !
!--------------------------------------------------------------------------------------------------!

#:include 'common.fypp'

!> Helper routines for transition charges between levels.
module dftbp_timedep_transcharges
  use dftbp_common_accuracy, only : dp
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
  subroutine TTransCharges_init(this, iAtomStart, sTimesGrndEigVecs, grndEigVecs, nTrans, nMatUp,&
      & nXooUD, nXvvUD, getia, getij, getab, win, tStoreOccVir, tStoreSame)

    !> Instance
    type(TTransCharges), intent(out) :: this

    !> Starting position of each atom in the list of orbitals
    integer, intent(in) :: iAtomStart(:)

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

    !> index array for single particle excitations that are included
    integer, intent(in) :: win(:)

    integer :: ia, ij, ii, jj, kk, ab, aa, bb
    logical :: updwn

    this%nTransitions = nTrans
    this%nAtom = size(iAtomStart) - 1
    this%nMatUp = nMatUp
    this%nMatUpOccOcc = nXooUD(1)
    this%nMatUpVirVir = nXvvUD(1)

    if (tStoreOccVir) then

      @:ASSERT(.not.allocated(this%qCacheOccVir))
      allocate(this%qCacheOccVir(this%nAtom, nTrans))

      !!$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(ia,ii,aa,kk,updwn) SCHEDULE(RUNTIME)
      do ia = 1, nTrans
        kk = win(ia)
        ii = getia(kk,1)
        aa = getia(kk,2)
        updwn = (kk <= this%nMatUp)
        this%qCacheOccVir(:,ia) = transq(ii, aa, iAtomStart, updwn,  sTimesGrndEigVecs,&
            & grndEigVecs)
      end do
      !!$OMP  END PARALLEL DO

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
        this%qCacheOccOcc(:,ij) = transq(ii, jj, iAtomStart, updwn,  sTimesGrndEigVecs,&
             & grndEigVecs)
      enddo
      !!$OMP  END PARALLEL DO

      !!$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(ab,aa,bb,updwn) SCHEDULE(RUNTIME)
      do ab = 1, sum(nXvvUD)
        aa = getab(ab,1)
        bb = getab(ab,2)
        updwn = (ab <= nXvvUD(1))
        this%qCacheVirVir(:,ab) = transq(aa, bb, iAtomStart, updwn,  sTimesGrndEigVecs,&
            & grndEigVecs)
      end do
      !!$OMP  END PARALLEL DO

      this%tCacheChargesSame = .true.

    else

      this%tCacheChargesSame = .false.

    end if

  end subroutine TTransCharges_init


  !> returns transition charges between occ-vir single particle levels
  pure function TTransCharges_qTransIA(this, ij, iAtomStart, sTimesGrndEigVecs, grndEigVecs, getia,&
      & win) result(q)

    !> instance of the transition charge object
    class(TTransCharges), intent(in) :: this

    !> Index of transition
    integer, intent(in) :: ij

    !> Starting position of each atom in the list of orbitals
    integer, intent(in) :: iAtomStart(:)

    !> Overlap times eigenvector: sum_m Smn cmi (nOrb, nOrb)
    real(dp), intent(in) :: sTimesGrndEigVecs(:,:,:)

    !> Eigenvectors (nOrb, nOrb)
    real(dp), intent(in) :: grndEigVecs(:,:,:)

    !> index array for for single particle excitations
    integer, intent(in) :: getia(:,:)

    !> index array for single particle excitations that are included
    integer, intent(in) :: win(:)

    !> Transition charge on exit. (nAtom)
    real(dp), dimension(size(iAtomStart)-1) :: q

    logical :: updwn
    integer :: ii, jj, kk

    if (allocated(this%qCacheOccVir)) then
      q(:) = this%qCacheOccVir(:, ij)
    else
      kk = win(ij)
      ii = getia(kk,1)
      jj = getia(kk,2)
      updwn = (kk <= this%nMatUp)
      q(:) = transq(ii, jj, iAtomStart, updwn, sTimesgrndEigVecs, grndEigVecs)
    end if

  end function TTransCharges_qTransIA

  !> returns transition charges between occ-occ single particle levels
  pure function TTransCharges_qTransIJ(this, ij, iAtomStart, sTimesGrndEigVecs, grndEigVecs, getij)&
      & result(q)

    !> instance of the transition charge object
    class(TTransCharges), intent(in) :: this

    !> Index of transition
    integer, intent(in) :: ij

    !> Starting position of each atom in the list of orbitals
    integer, intent(in) :: iAtomStart(:)

    !> Overlap times eigenvector: sum_m Smn cmi (nOrb, nOrb)
    real(dp), intent(in) :: sTimesGrndEigVecs(:,:,:)

    !> Eigenvectors (nOrb, nOrb)
    real(dp), intent(in) :: grndEigVecs(:,:,:)

    !> index array for occ-occ single particle excitations
    integer, intent(in) :: getij(:,:)

    !> Transition charge on exit. (nAtom)
    real(dp), dimension(size(iAtomStart)-1) :: q

    logical :: updwn
    integer :: ii, jj

    if (allocated(this%qCacheOccOcc)) then
      q(:) = this%qCacheOccOcc(:, ij)
    else
      ii = getij(ij,1)
      jj = getij(ij,2)
      updwn = (ij <= this%nMatUpOccOcc)
      q(:) = transq(ii, jj, iAtomStart, updwn, sTimesgrndEigVecs, grndEigVecs)
    end if

  end function TTransCharges_qTransIJ

  !> returns transition charges between vir-vir single particle levels
  pure function TTransCharges_qTransAB(this, ab, iAtomStart, sTimesGrndEigVecs, grndEigVecs, getab)&
      & result(q)

    !> instance of the transition charge object
    class(TTransCharges), intent(in) :: this

    !> Index of transition
    integer, intent(in) :: ab

    !> Starting position of each atom in the list of orbitals
    integer, intent(in) :: iAtomStart(:)

    !> Overlap times eigenvector: sum_m Smn cmi (nOrb, nOrb)
    real(dp), intent(in) :: sTimesGrndEigVecs(:,:,:)

    !> Eigenvectors (nOrb, nOrb)
    real(dp), intent(in) :: grndEigVecs(:,:,:)

    !> index array for occ-occ single particle excitations
    integer, intent(in) :: getab(:,:)

    !> Transition charge on exit. (nAtom)
    real(dp), dimension(size(iAtomStart)-1) :: q

    logical :: updwn
    integer :: aa, bb

    if (allocated(this%qCacheVirVir)) then
      q(:) = this%qCacheVirVir(:, ab)
    else
      aa = getab(ab,1)
      bb = getab(ab,2)
      updwn = (ab <= this%nMatUpVirVir)
      q(:) = transq(aa, bb, iAtomStart, updwn, sTimesgrndEigVecs, grndEigVecs)
    end if

  end function TTransCharges_qTransAB


  !> Transition charges left produced with a vector Q * v
  pure subroutine TTransCharges_qMatVec(this, iAtomStart, sTimesGrndEigVecs, grndEigVecs, getia,&
      & win, vector, qProduct)

    !> instance of the transition charge object
    class(TTransCharges), intent(in) :: this

    !> Starting position of each atom in the list of orbitals
    integer, intent(in) :: iAtomStart(:)

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

    if (this%tCacheChargesOccVir) then

      qProduct(:) = qProduct + matmul(this%qCacheOccVir, vector)

    else

      allocate(qij(this%nAtom))

      !!$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(ij,ii,jj,kk,updwn,qij)&
      !!$OMP& SCHEDULE(RUNTIME) REDUCTION(+:qProduct)
      do ij = 1, this%nTransitions
        kk = win(ij)
        ii = getia(kk,1)
        jj = getia(kk,2)
        updwn = (kk <= this%nMatUp)
        qij(:) = transq(ii, jj, iAtomStart, updwn, sTimesGrndEigVecs, grndEigVecs)
        qProduct(:) = qProduct + qij * vector(ij)
      end do
      !!$OMP  END PARALLEL DO

      deallocate(qij)

    end if

  end subroutine TTransCharges_qMatVec


  !> Transition charges right produced with a vector v * Q
  pure subroutine TTransCharges_qVecMat(this, iAtomStart, sTimesGrndEigVecs, grndEigVecs, getia,&
      & win, vector, qProduct)

    !> instance of the transition charge object
    class(TTransCharges), intent(in) :: this

    !> Starting position of each atom in the list of orbitals
    integer, intent(in) :: iAtomStart(:)

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

    if (this%tCacheChargesOccVir) then

      qProduct(:) = qProduct + matmul(vector, this%qCacheOccVir)

    else

      allocate(qij(this%nAtom))

      !!$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(ij,ii,jj,kk,updwn,qij)&
      !!$OMP& SCHEDULE(RUNTIME)
      do ij = 1, this%nTransitions
        kk = win(ij)
        ii = getia(kk,1)
        jj = getia(kk,2)
        updwn = (kk <= this%nMatUp)
        qij(:) = transq(ii, jj, iAtomStart, updwn, sTimesGrndEigVecs, grndEigVecs)
        qProduct(ij) = qProduct(ij) + dot_product(qij, vector)
      end do
      !!$OMP  END PARALLEL DO

      deallocate(qij)

    end if

  end subroutine TTransCharges_qVecMat


  !> Transition charges left produced with a vector Q * v for spin up
  !> minus Transition charges left produced with a vector Q * v for spin down
  !> sum_ias q_ias V_ias delta_s,  where delta_s = 1 for spin up and delta_s = -1 for spin down
  pure subroutine TTransCharges_qMatVecDs(this, iAtomStart, sTimesGrndEigVecs, grndEigVecs, getia,&
      & win, vector, qProduct)

    !> instance of the transition charge object
    class(TTransCharges), intent(in) :: this

    !> Starting position of each atom in the list of orbitals
    integer, intent(in) :: iAtomStart(:)

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
      qij(:) = transq(ii, jj, iAtomStart, updwn, sTimesGrndEigVecs, grndEigVecs)
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
  pure subroutine TTransCharges_qVecMatDs(this, iAtomStart, sTimesGrndEigVecs, grndEigVecs, getia,&
      & win, vector, qProduct)

    !> instance of the transition charge object
    class(TTransCharges), intent(in) :: this

    !> Starting position of each atom in the list of orbitals
    integer, intent(in) :: iAtomStart(:)

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
      qij(:) = transq(ii, jj, iAtomStart, updwn, sTimesGrndEigVecs, grndEigVecs)
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
  pure function transq(ii, jj, iAtomStart, updwn, sTimesGrndEigVecs, grndEigVecs) result(qij)

    !> Index of initial state.
    integer, intent(in) :: ii

    !> Index of final state.
    integer, intent(in) :: jj

    !> Starting position of each atom in the list of orbitals.
    integer, intent(in) :: iAtomStart(:)

    !> up spin channel (T) or down spin channel (F)
    logical, intent(in) :: updwn

    !> Overlap times eigenvector: sum_m Smn cmi (nOrb, nOrb).
    real(dp), intent(in) :: sTimesGrndEigVecs(:,:,:)

    !> Eigenvectors (nOrb, nOrb)
    real(dp), intent(in) :: grndEigVecs(:,:,:)

    !> Transition charge on exit. (nAtom)
    real(dp) :: qij(size(iAtomStart)-1)

    integer :: kk, aa, bb, ss
    real(dp) :: qTmp(size(grndEigVecs,dim=1))

    ss = 1
    if (.not. updwn) then
      ss = 2
    end if

    qTmp(:) =  grndEigVecs(:,ii,ss) * sTimesGrndEigVecs(:,jj,ss)&
        & + grndEigVecs(:,jj,ss) * sTimesGrndEigVecs(:,ii,ss)
    do kk = 1, size(qij)
      aa = iAtomStart(kk)
      bb = iAtomStart(kk + 1) -1
      qij(kk) = 0.5_dp * sum(qTmp(aa:bb))
    end do

  end function transq


end module dftbp_timedep_transcharges
