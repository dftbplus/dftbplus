!--------------------------------------------------------------------------------------------------!
!  DFTB+: general package for performing fast atomistic simulations                                !
!  Copyright (C) 2018  DFTB+ developers group                                                      !
!                                                                                                  !
!  See the LICENSE file for terms of usage and distribution.                                       !
!--------------------------------------------------------------------------------------------------!

#:include 'common.fypp'

!> Helper routines for transition charges between levels.
module transitionCharges
  use assert
  use accuracy
  implicit none

  public :: qTransition, qTransitionInit

  !> Internal data of the transition charges
  type :: qTransition
    private

    !> should transition charges be cached in memory or evaluated when needed?
    logical :: tCacheCharges = .false.

    !> storage if caching the occupied -> virtual transition charges
    real(dp), allocatable :: qCacheOccVirt(:,:)

    !> Number of transitions in the cache
    integer :: nTransitions

    !> Number of atoms in the system
    integer :: nAtom

    !> Number of transitions within the spin-up block
    integer :: nMatUp

  contains

    procedure, nopass :: qTransIJ
    procedure, nopass :: qMatVec

  end type qTransition

contains

  subroutine qTransitionInit(this, iAtomStart, sTimesGrndEigVecs, grndEigVecs, nxov_rd, nMatUp,&
      & getij, win, tStore)

    type(qTransition), intent(inout) :: this

    !> Starting position of each atom in the list of orbitals
    integer, intent(in) :: iAtomStart(:)

    !> Overlap times eigenvector: sum_m Smn cmi (nOrb, nOrb)
    real(dp), intent(in) :: sTimesGrndEigVecs(:,:,:)

    !> Eigenvectors (nOrb, nOrb)
    real(dp), intent(in) :: grndEigVecs(:,:,:)

    !> number of transitions in the system
    integer, intent(in) :: nxov_rd

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

    this%nTransitions = nxov_rd
    this%nAtom = size(iAtomStart) - 1
    this%nMatUp = nMatUp

    if (tStore) then

    @:ASSERT(.not.allocated(this%qCacheOccVirt))
      allocate(this%qCacheOccVirt(this%nAtom,nxov_rd))

      !!$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(ij,ii,jj,updwn) SCHEDULE(RUNTIME)
      do ij = 1, nxov_rd
        kk = win(ij)
        ii = getij(kk,1)
        jj = getij(kk,2)
        !call indxov(win, ij, getij, ii, jj)
        updwn = (kk <= this%nMatUp)
        call transq(ii, jj, iAtomStart, updwn,  sTimesGrndEigVecs, grndEigVecs,&
            & this%qCacheOccVirt(:,ij))
      end do
      !!$OMP  END PARALLEL DO

      this%tCacheCharges = .true.

    else

      this%tCacheCharges = .false.

    end if

  end subroutine qTransitionInit


  !> returns transtion charges between single particle levels
  pure function qTransIJ(this, ij, iAtomStart, sTimesGrndEigVecs, grndEigVecs, getij, win) result(q)

    !> instance of the transition charge object
    class(qTransition), intent(in) :: this

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
      !call indxov(win, ij, getij, ii, jj)
      updwn = (kk <= this%nMatUp)
      call transq(ii, jj, iAtomStart, updwn, sTimesgrndEigVecs, grndEigVecs, q(:))
    end if

  end function qTransIJ


  !> Transition charges producted with a matrix
  pure subroutine qMatVec(this, iAtomStart, sTimesGrndEigVecs, grndEigVecs, getij, win, vector,&
        & qProduct)

    !> instance of the transition charge object
    class(qTransition), intent(in) :: this

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
    real(dp), intent(out) :: qProduct(:)

    real(dp), allocatable :: qij(:)
    integer :: ii, jj, ij, kk
    logical :: updwn

    if (this%tCacheCharges) then

      qProduct(:) = matmul(this%qCacheOccVirt(:, :), vector(:))

    else

      allocate(qij(this%nAtom))

      !!$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(ij,ii,jj,updwn,qij)&
      !!$OMP& SCHEDULE(RUNTIME)
      do ij = 1, this%nTransitions
        kk = win(ij)
        ii = getij(kk,1)
        jj = getij(kk,2)
        !call indxov(win, ij, getij, ii, jj)
        updwn = (kk <= this%nMatUp)
        call transq(ii, jj, iAtomStart, updwn, sTimesGrndEigVecs, grndEigVecs, qij(:))
        qProduct(ij) = dot_product(qij(:), vector)
      end do
      !!$OMP  END PARALLEL DO

      deallocate(qij)

    end if

  end subroutine qMatVec


  !> Calculates atomic transition charges for a specified excitation.
  !> Calculates qij = 0.5 * (c_i S c_j + c_j S c_i) where c_i and c_j are selected eigenvectors, and
  !> S the overlap matrix.
  !> Since qij is atomic quantity (so far) the corresponding values for the atom are summed up.
  !> Note: the parameters 'updwn' were added for spin alpha and beta channels.
  pure subroutine transq(ii, jj, iAtomStart, updwn, sTimesGrndEigVecs, grndEigVecs, qij)

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
    real(dp), intent(out) :: qij(:)

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

  end subroutine transq


end module transitionCharges
