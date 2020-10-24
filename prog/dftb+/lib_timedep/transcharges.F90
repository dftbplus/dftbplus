!--------------------------------------------------------------------------------------------------!
!  DFTB+: general package for performing fast atomistic simulations                                !
!  Copyright (C) 2006 - 2020  DFTB+ developers group                                               !
!                                                                                                  !
!  See the LICENSE file for terms of usage and distribution.                                       !
!--------------------------------------------------------------------------------------------------!

#:include 'common.fypp'

!> Helper routines for transition charges between levels.
module dftbp_transcharges
  use dftbp_assert
  use dftbp_accuracy
  implicit none
  private

  public :: TTransCharges, TTransCharges_init
  public :: transq

  !> Internal data of the transition charges
  type :: TTransCharges
    private

    !> should transition charges be cached in memory or evaluated when needed?
    logical :: tCacheCharges

    !> storage if caching the occupied -> virtual transition charges
    real(dp), allocatable :: qCacheOccVirt(:,:)

    !> Number of transitions in the cache
    integer :: nTransitions

    !> Number of atoms in the system
    integer :: nAtom

    !> Number of transitions within the spin-up block
    integer :: nMatUp

  contains

    procedure :: qTransIJ =>TTransCharges_qTransIJ
    procedure :: qMatVec => TTransCharges_qMatVec
    procedure :: qVecMat => TTransCharges_qVecMat
    procedure :: qMatVecDs => TTransCharges_qMatVecDs
    procedure :: qVecMatDs => TTransCharges_qVecMatDs

  end type TTransCharges


contains

  !> initialise the cache/on-the fly transition charge evaluator
  subroutine TTransCharges_init(this, iAtomStart, sTimesGrndEigVecs, grndEigVecs, nTrans, nMatUp,&
      & getij, win, tStore)

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

    !> number of up-up excitations
    integer, intent(in) :: nMatUp

    !> should transitions be stored?
    logical, intent(in) :: tStore

    !> index array for for single particle excitations
    integer, intent(in) :: getij(:,:)

    !> index array for single particle excitions that are included
    integer, intent(in) :: win(:)

    integer :: ij, ii, jj, kk
    logical :: updwn

    this%nTransitions = nTrans
    this%nAtom = size(iAtomStart) - 1
    this%nMatUp = nMatUp

    if (tStore) then

      @:ASSERT(.not.allocated(this%qCacheOccVirt))
      allocate(this%qCacheOccVirt(this%nAtom, nTrans))

      !!$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(ij,ii,jj,kk,updwn) SCHEDULE(RUNTIME)
      do ij = 1, nTrans
        kk = win(ij)
        ii = getij(kk,1)
        jj = getij(kk,2)
        !call indxov(win, ij, getij, ii, jj)
        updwn = (kk <= this%nMatUp)
        this%qCacheOccVirt(:,ij) = transq(ii, jj, iAtomStart, updwn,  sTimesGrndEigVecs,&
            & grndEigVecs)
      end do
      !!$OMP  END PARALLEL DO

      this%tCacheCharges = .true.

    else

      this%tCacheCharges = .false.

    end if

  end subroutine TTransCharges_init


  !> returns transtion charges between single particle levels
  pure function TTransCharges_qTransIJ(this, ij, iAtomStart, sTimesGrndEigVecs, grndEigVecs, getij,&
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
    integer, intent(in) :: getij(:,:)

    !> index array for single particle excitions that are included
    integer, intent(in) :: win(:)

    !> Transition charge on exit. (nAtom)
    real(dp), dimension(size(iAtomStart)-1) :: q

    logical :: updwn
    integer :: ii, jj, kk

    if (allocated(this%qCacheOccVirt)) then
      q(:) = this%qCacheOccVirt(:, ij)
    else
      kk = win(ij)
      ii = getij(kk,1)
      jj = getij(kk,2)
      updwn = (kk <= this%nMatUp)
      q(:) = transq(ii, jj, iAtomStart, updwn, sTimesgrndEigVecs, grndEigVecs)
    end if

  end function TTransCharges_qTransIJ


  !> Transition charges left producted with a vector Q * v
  pure subroutine TTransCharges_qMatVec(this, iAtomStart, sTimesGrndEigVecs, grndEigVecs, getij,&
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
    integer, intent(in) :: getij(:,:)

    !> index array for single particle excitions that are included
    integer, intent(in) :: win(:)

    !> vector to product with the transition charges
    real(dp), intent(in) :: vector(:)

    !> Product on exit
    real(dp), intent(inout) :: qProduct(:)

    real(dp), allocatable :: qij(:)
    integer :: ii, jj, ij, kk
    logical :: updwn

    if (this%tCacheCharges) then

      qProduct(:) = qProduct + matmul(this%qCacheOccVirt, vector)

    else

      allocate(qij(this%nAtom))

      !!$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(ij,ii,jj,kk,updwn,qij)&
      !!$OMP& SCHEDULE(RUNTIME) REDUCTION(+:qProduct)
      do ij = 1, this%nTransitions
        kk = win(ij)
        ii = getij(kk,1)
        jj = getij(kk,2)
        updwn = (kk <= this%nMatUp)
        qij(:) = transq(ii, jj, iAtomStart, updwn, sTimesGrndEigVecs, grndEigVecs)
        qProduct(:) = qProduct + qij * vector(ij)
      end do
      !!$OMP  END PARALLEL DO

      deallocate(qij)

    end if

  end subroutine TTransCharges_qMatVec


  !> Transition charges right producted with a vector v * Q
  pure subroutine TTransCharges_qVecMat(this, iAtomStart, sTimesGrndEigVecs, grndEigVecs, getij,&
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
    integer, intent(in) :: getij(:,:)

    !> index array for single particle excitions that are included
    integer, intent(in) :: win(:)

    !> vector to product with the transition charges
    real(dp), intent(in) :: vector(:)

    !> Product on exit
    real(dp), intent(inout) :: qProduct(:)

    real(dp), allocatable :: qij(:)
    integer :: ii, jj, ij, kk
    logical :: updwn

    if (this%tCacheCharges) then

      qProduct(:) = qProduct + matmul(vector, this%qCacheOccVirt)

    else

      allocate(qij(this%nAtom))

      !!$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(ij,ii,jj,kk,updwn,qij)&
      !!$OMP& SCHEDULE(RUNTIME)
      do ij = 1, this%nTransitions
        kk = win(ij)
        ii = getij(kk,1)
        jj = getij(kk,2)
        updwn = (kk <= this%nMatUp)
        qij(:) = transq(ii, jj, iAtomStart, updwn, sTimesGrndEigVecs, grndEigVecs)
        qProduct(ij) = qProduct(ij) + dot_product(qij, vector)
      end do
      !!$OMP  END PARALLEL DO

      deallocate(qij)

    end if

  end subroutine TTransCharges_qVecMat


  !> Transition charges left producted with a vector Q * v for spin up
  !> minus Transition charges left producted with a vector Q * v for spin down
  !> sum_ias q_ias V_ias delta_s,  where delta_s = 1 for spin up and delta_s = -1 for spin down
  pure subroutine TTransCharges_qMatVecDs(this, iAtomStart, sTimesGrndEigVecs, grndEigVecs, getij,&
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
    integer, intent(in) :: getij(:,:)

    !> index array for single particle excitions that are included
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
      ii = getij(kk,1)
      jj = getij(kk,2)
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


  !> Transition charges right producted with a vector v * Q for spin up
  !> and negative transition charges right producted with a vector v * Q for spin down
  !> R_ias = delta_s sum_A q_A^(ias) V_A,  where delta_s = 1 for spin up and delta_s = -1 for spin down
  pure subroutine TTransCharges_qVecMatDs(this, iAtomStart, sTimesGrndEigVecs, grndEigVecs, getij,&
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
    integer, intent(in) :: getij(:,:)

    !> index array for single particle excitions that are included
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
      ii = getij(kk,1)
      jj = getij(kk,2)
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

    !> Index of inital state.
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


end module dftbp_transcharges
