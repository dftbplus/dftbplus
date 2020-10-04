!--------------------------------------------------------------------------------------------------!
!  DFTB+: general package for performing fast atomistic simulations                                !
!  Copyright (C) 2006 - 2020  DFTB+ developers group                                               !
!                                                                                                  !
!  See the LICENSE file for terms of usage and distribution.                                       !
!--------------------------------------------------------------------------------------------------!

#:include 'common.fypp'

!> Calculate either the whole single particle density matrix and energy
!> weighted density matrix or only the elements dictated by a neighbour map.
!> Calculation of the whole matrix scales as O(N**3), the sparse form as
!> O(N**2) but with a larger pre-factor.
!> Note: Dense code based on suggestions from Thomas Heine
!> Caveat: The routines create the transposed and complex conjugated of the density matrices! (cc*
!> instead of the conventional c*c)
module dftbp_densitymatrix
  use dftbp_assert
  use dftbp_accuracy
  use dftbp_blasroutines
  use dftbp_sorting
  use dftbp_commontypes
#:if WITH_SCALAPACK
  use dftbp_scalapackfx
  use dftbp_blacsenv
#:endif
  implicit none
  private

  public :: makeDensityMatrix

#:if WITH_SCALAPACK
  public :: makeDensityMtxRealBlacs, makeDensityMtxCplxBlacs
#:endif


  !> Provides an interface to calculate the two types of dm - regular and
  !> weighted and put them into packed storage
  interface makeDensityMatrix
    ! dense cases
    module procedure fullDensityMatrix_real
    module procedure fullDensityMatrix_cmplx
    module procedure fullEnergyDensityMatrix_real
    module procedure fullEnergyDensityMatrix_cmplx
    ! now the sparse routines
    module procedure sp_density_matrix_real
    module procedure sp_density_matrix_cmplx
    module procedure sp_energy_density_matrix_real
    module procedure sp_energy_density_matrix_cmplx
  end interface makeDensityMatrix

  real(dp), parameter :: arbitraryConstant = 0.1_dp

contains


  !> Make a regular density matrix for the real wave-function case
  !> Note: In order to save memory, the eigenvectors (which should be intent(in) parameters) are
  !> overwritten and then restored again
  subroutine fullDensityMatrix_real(dm, eigenvecs, filling)

    !> the resulting nOrb*nOrb density matrix
    real(dp), intent(out) :: dm(:,:)

    !> the eigenvectors of the system
    real(dp), intent(inout) :: eigenvecs(:,:)

    !> the occupation numbers of the orbitals
    real(dp), intent(in) :: filling(:)

    integer :: ii, nLevels
    real(dp) :: shift

    @:ASSERT(all(shape(eigenvecs) == shape(dm)))
    @:ASSERT(size(eigenvecs,dim=1) == size(eigenvecs,dim=2))
    @:ASSERT(size(eigenvecs,dim=1) == size(filling))

    dm(:,:) = 0.0_dp
    do ii =  size(filling), 1, -1
      nLevels = ii
      if (abs(filling(ii)) >= epsilon(1.0_dp)) then
        exit
      end if
    end do
    shift = minval(filling(1:nLevels))
    if (shift > epsilon(1.0_dp)) then
      ! all fillings are definitely positive

      !$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(ii) SCHEDULE(RUNTIME)
      do ii = 1, nLevels
        eigenvecs(:,ii) = sqrt(filling(ii)) * eigenvecs(:,ii)
      end do
      !$OMP  END PARALLEL DO

      call herk(dm, eigenvecs(:,1:nLevels))
      !$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(ii) SCHEDULE(RUNTIME)
      do ii = 1, nLevels
        eigenvecs(:,ii) = eigenvecs(:,ii) / sqrt(filling(ii))
      end do
      !$OMP  END PARALLEL DO

    else

      ! shift matrix so that filling operations are positive
      call herk(dm, eigenvecs(:,1:nLevels))
      shift = shift - arbitraryConstant
      !$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(ii) SCHEDULE(RUNTIME)
      do ii = 1, nLevels
        eigenvecs(:,ii) = sqrt(filling(ii)-shift) * eigenvecs(:,ii)
      end do
      !$OMP  END PARALLEL DO

      call herk(dm, eigenvecs(:,1:nLevels), beta=shift)
      !$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(ii) SCHEDULE(RUNTIME)
      do ii = 1, nLevels
        eigenvecs(:,ii) = eigenvecs(:,ii) / sqrt(filling(ii)-shift)
      end do
      !$OMP  END PARALLEL DO

    end if
  end subroutine fullDensityMatrix_real


  !> Make a regular density matrix for the complex wave-function case
  !> Note: in order to save memory, the eigenvectors (which should be intent(in) parameters) are
  !> overwritten and then restored again
  subroutine fullDensityMatrix_cmplx(dm, eigenvecs, filling)

    !> the resulting nOrb*nOrb density matrix
    complex(dp), intent(out) :: dm(:,:)

    !> the eigenvectors of the system
    complex(dp), intent(inout) :: eigenvecs(:,:)

    !> the occupation numbers of the orbitals
    real(dp), intent(in) :: filling(:)

    integer :: ii, nLevels
    real(dp) :: shift

    @:ASSERT(all(shape(eigenvecs) == shape(dm)))
    @:ASSERT(size(eigenvecs,dim=1) == size(eigenvecs,dim=2))
    @:ASSERT(size(eigenvecs,dim=1) == size(filling))

    dm(:,:) = cmplx(0.0_dp,0.0_dp,dp)

    do ii =  size(filling), 1, -1
      nLevels = ii
      if (abs(filling(ii)) >= epsilon(1.0_dp)) then
        exit
      end if
    end do
    shift = minval(filling(1:nLevels))
    if (shift > epsilon(1.0_dp)) then
      ! all fillings are definitely positive
      !$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(ii) SCHEDULE(RUNTIME)
      do ii = 1, nLevels
        eigenvecs(:,ii) = sqrt(filling(ii))*eigenvecs(:,ii)
      end do
      !$OMP  END PARALLEL DO

      call herk(dm, eigenvecs(:,1:nLevels))
      !$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(ii) SCHEDULE(RUNTIME)
      do ii = 1, nLevels
        eigenvecs(:,ii) = eigenvecs(:,ii) / sqrt(filling(ii))
      end do
      !$OMP  END PARALLEL DO

    else
      ! shift matrix so that filling operations are positive
      call herk(dm, eigenvecs(:,1:nLevels))
      shift = shift - arbitraryConstant
      !$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(ii) SCHEDULE(RUNTIME)
      do ii = 1, nLevels
        eigenvecs(:,ii) = sqrt(filling(ii)-shift)*eigenvecs(:,ii)
      end do
      !$OMP  END PARALLEL DO

      call herk(dm, eigenvecs(:,1:nLevels), beta=shift)
      !$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(ii) SCHEDULE(RUNTIME)
      do ii = 1, nLevels
        eigenvecs(:,ii) = eigenvecs(:,ii) / sqrt(filling(ii)-shift)
      end do
      !$OMP  END PARALLEL DO

    end if
  end subroutine fullDensityMatrix_cmplx


  !> Make an energy weighted density matrix for the real wave-function case
  !> Note: in order to save memory, the eigenvectors (which should be intent(in) parameters) are
  !> overwritten and then restored again
  subroutine fullEnergyDensityMatrix_real(dm, eigenvecs, filling, eigen)

    !> the resulting nOrb*nOrb density matrix
    real(dp), intent(out) :: dm(:,:)

    !> the eigenvectors of the system
    real(dp), intent(inout) :: eigenvecs(:,:)

    !> the occupation numbers of the orbitals
    real(dp), intent(in) :: filling(:)

    !> eigenvalues of the system
    real(dp), intent(in) :: eigen(:)

    integer :: ii, nLevels
    real(dp) :: shift
    real(dp) :: fillProduct(size(filling))

    @:ASSERT(all(shape(eigenvecs) == shape(dm)))
    @:ASSERT(size(eigenvecs,dim=1) == size(eigenvecs,dim=2))
    @:ASSERT(size(eigenvecs,dim=1) == size(filling))
    @:ASSERT(size(eigen) == size(filling))

    dm(:,:) = 0.0_dp
    do ii =  size(filling), 1, -1
      nLevels = ii
      if (abs(filling(ii)) >= epsilon(1.0_dp)) then
        exit
      end if
    end do
    fillProduct(1:nLevels) = filling(1:nLevels) * eigen(1:nLevels)
    if ((minval(fillProduct(1:nLevels)) < 0.0_dp&
        & .eqv. maxval(fillProduct(1:nLevels)) < 0.0_dp)&
        & .and. abs(minval(fillProduct(1:nLevels))) > epsilon(1.0_dp)&
        & .and. abs(maxval(fillProduct(1:nLevels))) > epsilon(1.0_dp)) then
      ! all fillings the same sign, and fairly large

      !$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(ii) SCHEDULE(RUNTIME)
      do ii = 1, nLevels
        eigenvecs(:,ii) = sqrt(abs(fillProduct(ii)))*eigenvecs(:,ii)
      end do
      !$OMP  END PARALLEL DO

      call herk(dm, eigenvecs(:,1:nLevels), alpha=sign(1.0_dp, maxval(fillProduct(1:nLevels))))
      !$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(ii) SCHEDULE(RUNTIME)
      do ii = 1, nLevels
        eigenvecs(:,ii) = eigenvecs(:,ii) / sqrt(abs(fillProduct(ii)))
      end do
      !$OMP  END PARALLEL DO

    else

      ! shift matrix so that filling operations are positive
      call herk(dm, eigenvecs(:,1:nLevels))
      shift = minval(fillProduct(1:nLevels)) - arbitraryConstant
      !$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(ii) SCHEDULE(RUNTIME)
      do ii = 1, nLevels
        eigenvecs(:,ii) = sqrt(fillProduct(ii)-shift) * eigenvecs(:,ii)
      end do
      !$OMP  END PARALLEL DO

      call herk(dm, eigenvecs(:,1:nLevels), beta=shift)
      !$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(ii) SCHEDULE(RUNTIME)
      do ii = 1, nLevels
        eigenvecs(:,ii) = eigenvecs(:,ii) / sqrt(fillProduct(ii)-shift)
      end do
      !$OMP  END PARALLEL DO

    end if
  end subroutine fullEnergyDensityMatrix_real


  !> Make an energy weighted density matrix for the complex wave-function case
  !> Note: in order to save memory, the eigenvectors (which should be intent(in) parameters) are
  !> overwritten and then restored again
  subroutine fullEnergyDensityMatrix_cmplx(dm, eigenvecs, filling, eigen)

    !> the resulting nOrb*nOrb density matrix
    complex(dp), intent(out) :: dm(:,:)

    !> the eigenvectors of the system
    complex(dp), intent(inout) :: eigenvecs(:,:)

    !> the occupation numbers of the orbitals
    real(dp), intent(in) :: filling(:)

    !> eigenvalues of the system
    real(dp), intent(in) :: eigen(:)

    integer :: ii, nLevels
    real(dp) :: shift
    real(dp) :: fillProduct(size(filling))

    @:ASSERT(all(shape(eigenvecs) == shape(dm)))
    @:ASSERT(size(eigenvecs,dim=1) == size(eigenvecs,dim=2))
    @:ASSERT(size(eigenvecs,dim=1) == size(filling))
    @:ASSERT(size(eigen) == size(filling))

    dm(:,:) = cmplx(0.0_dp,0.0_dp,dp)

    do ii =  size(filling), 1, -1
      nLevels = ii
      if (abs(filling(ii)) >= epsilon(1.0_dp)) then
        exit
      end if
    end do

    fillProduct(1:nLevels) = filling(1:nLevels) * eigen(1:nLevels)
    if ((minval(fillProduct(1:nLevels)) < 0.0_dp&
        & .eqv. maxval(fillProduct(1:nLevels)) < 0.0_dp)&
        & .and. abs(minval(fillProduct(1:nLevels))) > epsilon(1.0_dp)&
        & .and. abs(maxval(fillProduct(1:nLevels))) > epsilon(1.0_dp)) then
      ! all fillings the same sign, and fairly large
      !$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(ii) SCHEDULE(RUNTIME)
      do ii = 1, nLevels
        eigenvecs(:,ii) = sqrt(abs(fillProduct(ii))) * eigenvecs(:,ii)
      end do
      !$OMP  END PARALLEL DO

      call herk(dm, eigenvecs(:,1:nLevels), alpha=sign(1.0_dp, maxval(fillProduct(1:nLevels))))
      !$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(ii) SCHEDULE(RUNTIME)
      do ii = 1, nLevels
        eigenvecs(:,ii) = eigenvecs(:,ii) / sqrt(abs(fillProduct(ii)))
      end do
      !$OMP  END PARALLEL DO

    else

      ! shift matrix so that filling operations are positive
      call herk(dm, eigenvecs(:,1:nLevels))
      shift = minval(fillProduct(1:nLevels)) - arbitraryConstant
      !$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(ii) SCHEDULE(RUNTIME)
      do ii = 1, nLevels
        eigenvecs(:,ii) = sqrt(fillProduct(ii)-shift)*eigenvecs(:,ii)
      end do
      !$OMP  END PARALLEL DO

      call herk(dm, eigenvecs(:,1:nLevels), beta=shift)
      !$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(ii) SCHEDULE(RUNTIME)
      do ii = 1, nLevels
        eigenvecs(:,ii) = eigenvecs(:,ii) / sqrt(fillProduct(ii)-shift)
      end do
      !$OMP  END PARALLEL DO

    end if
  end subroutine fullEnergyDensityMatrix_cmplx


  !> Make a regular density matrix for the real wave-function case
  subroutine sp_density_matrix_real(dm, eigenvecs, filling, iNeighbour, nNeighbourSK, orb,&
      & iAtomStart, img2CentCell)

    !> the resulting nOrb*nOrb density matrix with only the elements of interest
    !> calculated, instead of the whole dm
    real(dp), intent(out) :: dm(:,:)

    !> the eigenvectors of the system
    real(dp), intent(in) :: eigenvecs(:,:)

    !> the occupation numbers of the orbitals
    real(dp), intent(in) :: filling(:)

    !> Neighbour list for each atom (First index from 0!)
    integer, intent(in) :: iNeighbour(0:,:)

    !> Nr. of neighbours for each atom (incl. itself).
    integer, intent(in) :: nNeighbourSK(:)

    !> Information about the orbitals of the atoms
    type(TOrbitals), intent(in) :: orb

    !> Atom offset for the squared Hamiltonian
    integer, intent(in) :: iAtomStart(:)

    !> Atomic mapping indexes.
    integer, intent(in) :: img2CentCell(:)

    integer :: iAt1, iNeigh1, nOrb1, nOrb2, jj
    integer :: nAtom, nLevels, start1, start2, ii
    integer, allocatable :: inCellNeighbour(:,:)
    integer, allocatable :: nInCellNeighbour(:)
    real(dp), allocatable :: tmpEigen(:,:)

    @:ASSERT(all(shape(eigenvecs) == shape(dm)))
    @:ASSERT(size(eigenvecs,dim=1) == size(eigenvecs,dim=2))
    @:ASSERT(size(eigenvecs,dim=1) == size(filling))

    allocate(inCellNeighbour(0:size(iNeighbour,dim=1),size(iNeighbour,dim=2)))
    allocate(nInCellNeighbour(size(iNeighbour,dim=2)))

    nAtom = size(orb%nOrbAtom)
    dm(:,:) = 0.0_dp

    inCellNeighbour(:,:) = 0
    nInCellNeighbour(:) = 0

    do ii =  size(filling), 1, -1
      nLevels = ii
      if (abs(filling(ii)) >= epsilon(1.0_dp)) then
        exit
      end if
    end do

    allocate(tmpEigen(nLevels,orb%mOrb))
    !$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(iAt1) SCHEDULE(RUNTIME)
    do iAt1 = 1, nAtom
      inCellNeighbour(0:nNeighbourSK(iAt1),iAt1) =&
          & img2CentCell(iNeighbour(0:nNeighbourSK(iAt1),iAt1))
      call heap_sort(inCellNeighbour(:nNeighbourSK(iAt1),iAt1))
      nInCellNeighbour(iAt1) = unique(inCellNeighbour(:,iAt1), nNeighbourSK(iAt1)+1) - 1
    end do
    !$OMP  END PARALLEL DO

    do iAt1 = 1, nAtom
      nOrb1 = orb%nOrbAtom(iAt1)
      start1 = iAtomStart(iAt1)
      tmpEigen(1:nLevels,1:nOrb1) = transpose(eigenvecs(start1:start1+nOrb1-1,1:nLevels))
      !$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(jj) SCHEDULE(RUNTIME)
      do jj = 1, nLevels
        tmpEigen(jj,1:nOrb1) = filling(jj) * tmpEigen(jj,1:nOrb1)
      end do
      !$OMP  END PARALLEL DO
      !$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(iNeigh1,start2,nOrb2) SCHEDULE(RUNTIME)
      do iNeigh1 = 0, nInCellNeighbour(iAt1)
        start2 = iAtomStart(inCellNeighbour(iNeigh1,iAt1))
        nOrb2 = orb%nOrbAtom(inCellNeighbour(iNeigh1,iAt1))
        dm(start2:start2+nOrb2-1, start1:start1+nOrb1-1) =&
            & matmul(eigenvecs(start2:start2+nOrb2-1,1:nLevels), tmpEigen(1:nLevels,1:nOrb1))
      end do
      !$OMP  END PARALLEL DO
    end do

  end subroutine sp_density_matrix_real


  !> Make a regular density matrix for the complex wave-function case
  subroutine sp_density_matrix_cmplx(dm, eigenvecs, filling, iNeighbour, nNeighbourSK, orb,&
      & iAtomStart, img2CentCell)

    !> the resulting nOrb*nOrb density matrix with only the elements of interest
    !> calculated, instead of the whole dm
    complex(dp), intent(out) :: dm(:,:)

    !> the eigenvectors of the system
    complex(dp), intent(in) :: eigenvecs(:,:)

    !> the occupation numbers of the orbitals
    real(dp), intent(in) :: filling(:)

    !> Neighbour list for each atom (First index from 0!)
    integer, intent(in) :: iNeighbour(0:,:)

    !> Nr. of neighbours for each atom (incl. itself).
    integer, intent(in) :: nNeighbourSK(:)

    !> Informatio about the orbitals.
    type(TOrbitals), intent(in) :: orb

    !> Atom offset for the squared Hamiltonian
    integer, intent(in) :: iAtomStart(:)

    !> Atomic mapping indexes.
    integer, intent(in) :: img2CentCell(:)

    integer :: iAt1, iNeigh1, nOrb1, nOrb2, jj, nAtom
    integer :: nLevels, start1, start2, ii
    integer, allocatable :: inCellNeighbour(:,:)
    integer, allocatable :: nInCellNeighbour(:)
    complex(dp), allocatable :: tmpEigen(:,:)

    @:ASSERT(all(shape(eigenvecs) == shape(dm)))
    @:ASSERT(size(eigenvecs,dim=1) == size(eigenvecs,dim=2))
    @:ASSERT(size(eigenvecs,dim=1) == size(filling))

    allocate(inCellNeighbour(0:size(iNeighbour,dim=1),size(iNeighbour,dim=2)))
    allocate(nInCellNeighbour(size(iNeighbour,dim=2)))

    nAtom = size(orb%nOrbAtom)

    dm(:,:) = cmplx(0.0_dp,0.0_dp,dp)
    inCellNeighbour(:,:) = 0
    nInCellNeighbour(:) = 0

    do ii =  size(filling), 1, -1
      nLevels = ii
      if (abs(filling(ii)) >= epsilon(1.0_dp)) then
        exit
      end if
    end do

    allocate(tmpEigen(nLevels, orb%mOrb))
    !$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(iAt1) SCHEDULE(RUNTIME)
    do iAt1 = 1, nAtom
      inCellNeighbour(0:nNeighbourSK(iAt1),iAt1) =&
          & img2CentCell(iNeighbour(0:nNeighbourSK(iAt1),iAt1))
      call heap_sort(inCellNeighbour(:nNeighbourSK(iAt1),iAt1))
      nInCellNeighbour(iAt1) = unique(inCellNeighbour(:,iAt1), nNeighbourSK(iAt1)+1) - 1
    end do
    !$OMP  END PARALLEL DO

    do iAt1 = 1, nAtom
      nOrb1 = orb%nOrbAtom(iAt1)
      start1 = iAtomStart(iAt1)
      tmpEigen(1:nLevels,1:nOrb1) = transpose(eigenvecs(start1:start1+nOrb1-1,1:nLevels))

      !$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(jj) SCHEDULE(RUNTIME)
      do jj = 1, nLevels
        tmpEigen(jj,1:nOrb1) = filling(jj)*conjg(tmpEigen(jj,1:nOrb1))
      end do
      !$OMP  END PARALLEL DO

      !$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(iNeigh1,start2,nOrb2) SCHEDULE(RUNTIME)
      do iNeigh1 = 0, nInCellNeighbour(iAt1)
        start2 = iAtomStart(inCellNeighbour(iNeigh1,iAt1))
        nOrb2 = orb%nOrbAtom(inCellNeighbour(iNeigh1,iAt1))
        dm(start2:start2+nOrb2-1, start1:start1+nOrb1-1) =&
            & matmul(eigenvecs(start2:start2+nOrb2-1,1:nLevels),&
            & tmpEigen(1:nLevels,1:nOrb1))
      end do
      !$OMP  END PARALLEL DO
    end do

  end subroutine sp_density_matrix_cmplx


  !> Make an energy weighted density matrix for the real wave-function case
  subroutine sp_energy_density_matrix_real(dm, eigenvecs, filling, eigen, iNeighbour, nNeighbourSK,&
      & orb, iAtomStart, img2CentCell)

    !> the resulting nOrb*nOrb density matrix with only the elements of interest
    !> calculated, instead of the whole dm
    real(dp), intent(out) :: dm(:,:)

    !> the eigenvectors of the system
    real(dp), intent(in) :: eigenvecs(:,:)

    !> the occupation numbers of the orbitals
    real(dp), intent(in) :: filling(:)

    !> eigenvalues of the system
    real(dp), intent(in) :: eigen(:)

    !> Neighbour list for each atom (First index from 0!)
    integer, intent(in) :: iNeighbour(0:,:)

    !> Nr. of neighbours for each atom (incl. itself).
    integer, intent(in) :: nNeighbourSK(:)

    !> Information about the orbitals.
    type(TOrbitals), intent(in) :: orb

    !> Atom offset for the squared Hamiltonian
    integer, intent(in) :: iAtomStart(:)

    !> Atomic mapping indexes.
    integer, intent(in) :: img2CentCell(:)

    integer :: iAt1, iNeigh1, jj, nOrb1, nOrb2, nAtom
    integer :: nLevels, start1, start2, ii
    integer, allocatable :: inCellNeighbour(:,:)
    integer, allocatable :: nInCellNeighbour(:)
    real(dp), allocatable :: tmpEigen(:,:)

    @:ASSERT(all(shape(eigenvecs) == shape(dm)))
    @:ASSERT(size(eigenvecs,dim=1) == size(eigenvecs,dim=2))
    @:ASSERT(size(eigenvecs,dim=1) == size(filling))
    @:ASSERT(size(eigen) == size(filling))

    allocate(inCellNeighbour(0:size(iNeighbour,dim=1),size(iNeighbour,dim=2)))
    allocate(nInCellNeighbour(size(iNeighbour,dim=2)))

    nAtom = size(orb%nOrbAtom)
    dm(:,:) = 0.0_dp

    inCellNeighbour(:,:) = 0
    nInCellNeighbour(:) = 0

    do ii =  size(filling), 1, -1
      nLevels = ii
      if (abs(filling(ii)) >= epsilon(1.0_dp)) then
        exit
      end if
    end do

    allocate(tmpEigen(nLevels,orb%mOrb))
    !$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(iAt1) SCHEDULE(RUNTIME)
    do iAt1 = 1, nAtom
      inCellNeighbour(0:nNeighbourSK(iAt1),iAt1) =&
          & img2CentCell(iNeighbour(0:nNeighbourSK(iAt1),iAt1))
      call heap_sort(inCellNeighbour(:nNeighbourSK(iAt1),iAt1))
      nInCellNeighbour(iAt1) = unique(inCellNeighbour(:,iAt1), nNeighbourSK(iAt1)+1) - 1
    end do
    !$OMP  END PARALLEL DO

    do iAt1 = 1, nAtom
      nOrb1 = orb%nOrbAtom(iAt1)
      start1 = iAtomStart(iAt1)
      tmpEigen(1:nLevels,1:nOrb1) = transpose(eigenvecs(start1:start1+nOrb1-1,1:nLevels))

      !$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(jj) SCHEDULE(RUNTIME)
      do jj = 1, nLevels
        tmpEigen(jj,1:nOrb1) = filling(jj)*eigen(jj)*tmpEigen(jj,1:nOrb1)
      end do
      !$OMP  END PARALLEL DO
      !$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(iNeigh1,start2,nOrb2) SCHEDULE(RUNTIME)
      do iNeigh1 = 0, nInCellNeighbour(iAt1)
        start2 = iAtomStart(inCellNeighbour(iNeigh1,iAt1))
        nOrb2 = orb%nOrbAtom(inCellNeighbour(iNeigh1,iAt1))
        dm(start2:start2+nOrb2-1, start1:start1+nOrb1-1) =&
            & matmul(eigenvecs(start2:start2+nOrb2-1,1:nLevels), tmpEigen(1:nLevels,1:nOrb1))
      end do
      !$OMP  END PARALLEL DO
    end do

  end subroutine sp_energy_density_matrix_real


  !> Make an energy weighted density matrix for the complex wave-function case
  subroutine sp_energy_density_matrix_cmplx(dm, eigenvecs, filling, eigen, iNeighbour,&
      & nNeighbourSK, orb, iAtomStart, img2CentCell)

    !> the resulting nOrb*nOrb density matrix with only the elements of interest
    !> calculated, instead of the whole dm
    complex(dp), intent(out) :: dm(:,:)

    !> the eigenvectors of the system
    complex(dp), intent(in) :: eigenvecs(:,:)

    !> the occupation numbers of the orbitals
    real(dp), intent(in) :: filling(:)

    !> Neighbour list for each atom (First index from 0!)
    real(dp), intent(in) :: eigen(:)

    !> Nr. of neighbours for each atom (incl. itself).
    integer, intent(in) :: iNeighbour(0:,:)

    !> Information about the orbitals.
    integer, intent(in) :: nNeighbourSK(:)

    !> Atom offset for the squared Hamiltonian
    type(TOrbitals), intent(in) :: orb

    !> Atomic mapping indexes.
    integer, intent(in) :: iAtomStart(:)

    integer, intent(in) :: img2CentCell(:)

    integer :: iAt1, iNeigh1, jj, nOrb1, nOrb2
    integer :: nAtom, nLevels, start1, start2, ii
    integer, allocatable :: inCellNeighbour(:,:)
    integer, allocatable :: nInCellNeighbour(:)
    complex(dp), allocatable :: tmpEigen(:,:)

    @:ASSERT(all(shape(eigenvecs) == shape(dm)))
    @:ASSERT(size(eigenvecs,dim=1) == size(eigenvecs,dim=2))
    @:ASSERT(size(eigenvecs,dim=1) == size(filling))
    @:ASSERT(size(eigen) == size(filling))

    allocate(inCellNeighbour(0:size(iNeighbour,dim=1),size(iNeighbour,dim=2)))
    allocate(nInCellNeighbour(size(iNeighbour,dim=2)))

    nAtom = size(orb%nOrbAtom)
    dm(:,:) = cmplx(0.0_dp,0.0_dp,dp)

    inCellNeighbour(:,:) = 0
    nInCellNeighbour(:) = 0

    do ii =  size(filling), 1, -1
      nLevels = ii
      if (abs(filling(ii)) >= epsilon(1.0_dp)) then
        exit
      end if
    end do

    allocate(tmpEigen(nLevels, orb%mOrb))
    !$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(iAt1) SCHEDULE(RUNTIME)
    do iAt1 = 1, nAtom
      inCellNeighbour(0:nNeighbourSK(iAt1),iAt1) =&
          & img2CentCell(iNeighbour(0:nNeighbourSK(iAt1),iAt1))
      call heap_sort(inCellNeighbour(:nNeighbourSK(iAt1),iAt1))
      nInCellNeighbour(iAt1) = unique(inCellNeighbour(:,iAt1), nNeighbourSK(iAt1)+1) - 1
    end do
    !$OMP  END PARALLEL DO

    do iAt1 = 1, nAtom
      nOrb1 = orb%nOrbAtom(iAt1)
      start1 = iAtomStart(iAt1)
      tmpEigen(1:nLevels,1:nOrb1) = transpose(eigenvecs(start1:start1+nOrb1-1,1:nLevels))

      !$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(jj) SCHEDULE(RUNTIME)
      do jj = 1, nLevels
        tmpEigen(jj,1:nOrb1) = filling(jj)*eigen(jj)*conjg(tmpEigen(jj,1:nOrb1))
      end do
      !$OMP  END PARALLEL DO
      !$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(iNeigh1,start2,nOrb2) SCHEDULE(RUNTIME)
      do iNeigh1 = 0, nInCellNeighbour(iAt1)
        start2 = iAtomStart(inCellNeighbour(iNeigh1,iAt1))
        nOrb2 = orb%nOrbAtom(inCellNeighbour(iNeigh1,iAt1))
        dm(start2:start2+nOrb2-1, start1:start1+nOrb1-1) =&
            &matmul(eigenvecs(start2:start2+nOrb2-1,1:nLevels), tmpEigen(1:nLevels,1:nOrb1))
      end do
      !$OMP  END PARALLEL DO
    end do

  end subroutine sp_energy_density_matrix_cmplx


#:if WITH_SCALAPACK
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!! Scalapack routines
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Create density or energy weighted density matrix (real).
  !!
  subroutine makeDensityMtxRealBlacs(myBlacs, desc, filling, eigenVecs, densityMtx, eigenVals)
    type(blacsgrid), intent(in) :: myBlacs
    integer, intent(in) :: desc(:)
    real(dp), target, intent(in) :: filling(:)
    real(dp), intent(inout) :: eigenVecs(:,:)
    real(dp), intent(out) :: densityMtx(:,:)
    real(dp), intent(in), optional :: eigenVals(:)

    integer  :: ii, jj, iGlob, iLoc, blockSize, nLevel
    real(dp) :: minFill, maxFill, alpha, beta
    real(dp), allocatable, target :: eFilling(:)
    real(dp), pointer :: myFilling(:)
    type(blocklist) :: blocks

    densityMtx(:,:) = 0.0_dp

    ! Make non-zero fillings positive definite
    nLevel = size(filling)
    do while (abs(filling(nLevel)) < epsilon(1.0_dp) .and. nLevel > 1)
      nLevel = nLevel - 1
    end do
    if (present(eigenVals)) then
      allocate(eFilling(nLevel))
      eFilling(:) = filling(1:nLevel) * eigenVals(1:nLevel)
      myFilling => eFilling
    else
      myFilling => filling
    end if
    minFill = minval(myFilling)
    maxFill = maxval(myFilling)
    if ((minFill < 0.0_dp .eqv. maxFill < 0.0_dp) .and. abs(minFill) >= epsilon(1.0_dp)&
        & .and. abs(maxFill) >= epsilon(1.0_dp)) then
      alpha = sign(1.0_dp, maxFill)
      beta = 0.0_dp
    else
      alpha = 1.0_dp
      beta = minFill - arbitraryConstant
      call pblasfx_psyrk(eigenVecs, desc, densityMtx, desc, kk=nLevel)
    end if

    ! Scale eigenvectors
    call blocks%init(myBlacs, desc, "c")
    do ii = 1, size(blocks)
      call blocks%getblock(ii, iGlob, iLoc, blockSize)
      do jj = 0, min(blockSize - 1, nLevel - iGlob)
        eigenVecs(:,iLoc + jj) = eigenVecs(:,iLoc + jj) * sqrt(abs(myFilling(iGlob + jj) - beta))
      end do
    end do

    ! Create matrix by rank-k update
    call pblasfx_psyrk(eigenVecs, desc, densityMtx, desc, kk=nLevel,&
        & alpha=alpha, beta=beta)

    ! Revert eigenvectors to their original value
    do ii = 1, size(blocks)
      call blocks%getblock(ii, iGlob, iLoc, blockSize)
      do jj = 0, min(blockSize - 1, nLevel - iGlob)
        eigenVecs(:,iLoc + jj) = eigenVecs(:,iLoc + jj) / sqrt(abs(myFilling(iGlob + jj) - beta))
      end do
    end do

  end subroutine makeDensityMtxRealBlacs


  !> Create density or energy weighted density matrix (complex).
  !!
  subroutine makeDensityMtxCplxBlacs(myBlacs, desc, filling, eigenVecs, densityMtx, eigenVals)
    type(blacsgrid), intent(in) :: myBlacs
    integer, intent(in) :: desc(:)
    real(dp), target, intent(in) :: filling(:)
    complex(dp), intent(inout) :: eigenVecs(:,:)
    complex(dp), intent(out) :: densityMtx(:,:)
    real(dp), intent(in), optional :: eigenVals(:)

    integer  :: ii, jj, iGlob, iLoc, blockSize, nLevel
    real(dp) :: minFill, maxFill, alpha, beta
    real(dp), allocatable, target :: eFilling(:)
    real(dp), pointer :: myFilling(:)
    type(blocklist) :: blocks

    densityMtx(:,:) = 0.0_dp

    ! Make non-zero fillings positive definite
    nLevel = size(filling)
    do while (abs(filling(nLevel)) < epsilon(1.0_dp) .and. nLevel > 1)
      nLevel = nLevel - 1
    end do
    if (present(eigenVals)) then
      allocate(eFilling(nLevel))
      eFilling(:) = filling(1:nLevel) * eigenVals(1:nLevel)
      myFilling => eFilling
    else
      myFilling => filling
    end if
    minFill = minval(myFilling)
    maxFill = maxval(myFilling)
    if ((minFill < 0.0_dp .eqv. maxFill < 0.0_dp)&
        & .and. abs(minFill) >= epsilon(1.0_dp)&
        & .and. abs(maxFill) >= epsilon(1.0_dp)) then
      alpha = sign(1.0_dp, maxFill)
      beta = 0.0_dp
    else
      alpha = 1.0_dp
      beta = minFill - arbitraryConstant
      call pblasfx_pherk(eigenVecs, desc, densityMtx, desc, kk=nLevel)
    end if

    ! Scale eigenvectors
    call blocks%init(myBlacs, desc, "c")
    do ii = 1, size(blocks)
      call blocks%getblock(ii, iGlob, iLoc, blockSize)
      do jj = 0, min(blockSize - 1, nLevel - iGlob)
        eigenVecs(:,iLoc + jj) = eigenVecs(:,iLoc + jj) * sqrt(abs(myFilling(iGlob + jj) - beta))
      end do
    end do

    ! Create matrix by rank-k update
    call pblasfx_pherk(eigenVecs, desc, densityMtx, desc, kk=nLevel,&
        & alpha=alpha, beta=beta)

    ! Revert eigenvectors to their original value
    do ii = 1, size(blocks)
      call blocks%getblock(ii, iGlob, iLoc, blockSize)
      do jj = 0, min(blockSize - 1, nLevel - iGlob)
        eigenVecs(:,iLoc + jj) = eigenVecs(:,iLoc + jj) / sqrt(abs(myFilling(iGlob + jj) - beta))
      end do
    end do

  end subroutine makeDensityMtxCplxBlacs

#:endif

end module dftbp_densitymatrix
