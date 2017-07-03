!--------------------------------------------------------------------------------------------------!
!  DFTB+: general package for performing fast atomistic simulations                                !
!  Copyright (C) 2017  DFTB+ developers group                                                      !
!                                                                                                  !
!  See the LICENSE file for terms of usage and distribution.                                       !
!--------------------------------------------------------------------------------------------------!

#:include 'common.fypp'

!!* Calculate either the whole single particle density matrix and energy
!!* weighted density matrix or only the elements dictated by a neighbor map.
!!* Calculation of the whole matrix scales as O(N**3), the sparse form as
!!* O(N**2) but with a larger pre-factor.
!!* @author Ben Hourahine
!!* @note Dense code based on implementation by Thomas Heine
!!* @caveat The routines create the transposed and complex conjugated of
!!*   the density matrices! (cc* instead of the conventional c*c)
module densitymatrix
  use assert
  use accuracy
  use blasroutines
  use sorting
  use commontypes
  implicit none
  private

  public :: makeDensityMatrix


  !!* Provides an interface to calculate the two types of dm - regular and
  !!* weighted and put them into packed storage
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
  end interface

  real(dp), parameter :: arbitraryConstant = 0.1_dp

contains


  !!* Make a regular density matrix for the real wave-function case
  !!* @param dm the resulting nOrb*nOrb density matrix
  !!* @param eigenvecs the eigenvectors of the system
  !!* @param filling the occupation numbers of the orbitals
  !!* @note In order to save memory, the eigenvectors (which should be
  !!* intent(in) parameters) are overwritten and then restored again
  subroutine fullDensityMatrix_real(dm, eigenvecs, filling)
    real(dp), intent(out) :: dm(:,:)
    real(dp), intent(inout) :: eigenvecs(:,:)
    real(dp), intent(in) :: filling(:)

    integer  :: ii, nLevels
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




  !!* Make a regular density matrix for the complex wave-function case
  !!* @param dm the resulting nOrb*nOrb density matrix
  !!* @param eigenvecs the eigenvectors of the system
  !!* @param filling the occupation numbers of the orbitals
  !!* @note In order to save memory, the eigenvectors (which should be
  !!* intent(in) parameters) are overwritten and then restored again
  subroutine fullDensityMatrix_cmplx(dm, eigenvecs, filling)
    complex(dp), intent(out) :: dm(:,:)
    complex(dp), intent(inout) :: eigenvecs(:,:)
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




  !!* Make an energy weighted density matrix for the real
  !!* wave-function case
  !!* @param dm the resulting nOrb*nOrb density matrix
  !!* @param eigenvecs the eigenvectors of the system
  !!* @param filling the occupation numbers of the orbitals
  !!* @param eigen eigenvalues of the system
  !!* @note In order to save memory, the eigenvectors (which should be
  !!* intent(in) parameters) are overwritten and then restored again
  subroutine fullEnergyDensityMatrix_real(dm, eigenvecs, filling, eigen)
    real(dp), intent(out) :: dm(:,:)
    real(dp), intent(inout) :: eigenvecs(:,:)
    real(dp), intent(in) :: filling(:)
    real(dp), intent(in) :: eigen(:)

    integer  :: ii, nLevels
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
    if ((minval(fillProduct(1:nLevels)) < 0.0_dp &
        &.eqv. maxval(fillProduct(1:nLevels)) < 0.0_dp) &
        &.and. abs(minval(fillProduct(1:nLevels))) > epsilon(1.0_dp) &
        &.and. abs(maxval(fillProduct(1:nLevels))) > epsilon(1.0_dp)) then
      ! all fillings the same sign, and fairly large

!$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(ii) SCHEDULE(RUNTIME)
      do ii = 1, nLevels
        eigenvecs(:,ii) = sqrt(abs(fillProduct(ii)))*eigenvecs(:,ii)
      end do
!$OMP  END PARALLEL DO

      call herk(dm, eigenvecs(:,1:nLevels), &
          &alpha=sign(1.0_dp, maxval(fillProduct(1:nLevels))))
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




  !!* Make an energy weighted density matrix for the complex
  !!* wave-function case
  !!* @param dm the resulting nOrb*nOrb density matrix
  !!* @param eigenvecs the eigenvectors of the system
  !!* @param filling the occupation numbers of the orbitals
  !!* @param eigen eigenvalues of the system
  !!* @note In order to save memory, the eigenvectors (which should be
  !!* intent(in) parameters) are overwritten and then restored again
  subroutine fullEnergyDensityMatrix_cmplx(dm, eigenvecs, filling, eigen)
    complex(dp), intent(out) :: dm(:,:)
    complex(dp), intent(inout) :: eigenvecs(:,:)
    real(dp), intent(in) :: filling(:)
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
        &.eqv. maxval(fillProduct(1:nLevels)) < 0.0_dp)&
        &.and. abs(minval(fillProduct(1:nLevels))) > epsilon(1.0_dp) &
        &.and. abs(maxval(fillProduct(1:nLevels))) > epsilon(1.0_dp)) then
      ! all fillings the same sign, and fairly large
!$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(ii) SCHEDULE(RUNTIME)
      do ii = 1, nLevels
        eigenvecs(:,ii) = sqrt(abs(fillProduct(ii))) * eigenvecs(:,ii)
      end do
!$OMP  END PARALLEL DO

      call herk(dm, eigenvecs(:,1:nLevels), &
          & alpha=sign(1.0_dp, maxval(fillProduct(1:nLevels))))
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




  !!* Make a regular density matrix for the real wave-function case
  !!* @param dm the resulting nOrb*nOrb density matrix with only the elements
  !!* of interest calculated, instead of the whole dm
  !!* @param eigenvecs the eigenvectors of the system
  !!* @param filling the occupation numbers of the orbitals
  !!* @param iNeighbor Neighbor list for each atom (First index from 0!)
  !!* @param nNeighbor Nr. of neighbors for each atom (incl. itself).
  !!* @param orb Information about the orbitals of the atoms
  !!* @param iAtomStart Atom offset for the squared Hamiltonian
  !!* @param img2CentCell Atomic mapping indexes.
  subroutine sp_density_matrix_real(dm, eigenvecs, filling, iNeighbor, &
      &nNeighbor, orb, iAtomStart, img2CentCell)
    real(dp), intent(out) :: dm(:,:)
    real(dp), intent(in) :: eigenvecs(:,:)
    real(dp), intent(in) :: filling(:)
    integer, intent(in) :: iNeighbor(0:,:)
    integer, intent(in) :: nNeighbor(:)
    type(TOrbitals), intent(in) :: orb
    integer, intent(in) :: iAtomStart(:)
    integer, intent(in) :: img2CentCell(:)

    integer                  :: iAt1, iNeigh1, nOrb1, nOrb2, jj
    integer                  :: nAtom, nLevels, start1, start2, ii
    integer, allocatable     :: inCellNeighbor(:,:)
    integer, allocatable     :: nInCellNeighbor(:)
    real(dp), allocatable    :: tmpEigen(:,:)

    @:ASSERT(all(shape(eigenvecs) == shape(dm)))
    @:ASSERT(size(eigenvecs,dim=1) == size(eigenvecs,dim=2))
    @:ASSERT(size(eigenvecs,dim=1) == size(filling))

    allocate(inCellNeighbor(0:size(iNeighbor,dim=1),size(iNeighbor,dim=2)))
    allocate(nInCellNeighbor(size(iNeighbor,dim=2)))

    nAtom = size(orb%nOrbAtom)
    dm(:,:) = 0.0_dp

    inCellNeighbor(:,:) = 0
    nInCellNeighbor(:) = 0

    do ii =  size(filling), 1, -1
      nLevels = ii
      if (abs(filling(ii)) >= epsilon(1.0_dp)) then
        exit
      end if
    end do

    allocate(tmpEigen(nLevels,orb%mOrb))
!$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(iAt1) SCHEDULE(RUNTIME)
    do iAt1 = 1, nAtom
      inCellNeighbor(0:nNeighbor(iAt1),iAt1) = &
          &img2CentCell(iNeighbor(0:nNeighbor(iAt1),iAt1))
      call heap_sort(inCellNeighbor(:nNeighbor(iAt1),iAt1))
      nInCellNeighbor(iAt1) = &
          &unique(inCellNeighbor(:,iAt1), nNeighbor(iAt1)+1) - 1
    end do
!$OMP  END PARALLEL DO

    do iAt1 = 1, nAtom
      nOrb1 = orb%nOrbAtom(iAt1)
      start1 = iAtomStart(iAt1)
      tmpEigen(1:nLevels,1:nOrb1) = &
          & transpose(eigenvecs(start1:start1+nOrb1-1,1:nLevels))
!$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(jj) SCHEDULE(RUNTIME)
      do jj = 1, nLevels
        tmpEigen(jj,1:nOrb1) = filling(jj) * tmpEigen(jj,1:nOrb1)
      end do
!$OMP  END PARALLEL DO
!$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(iNeigh1,start2,nOrb2) &
!$OMP& SCHEDULE(RUNTIME)
      do iNeigh1 = 0, nInCellNeighbor(iAt1)
        start2 = iAtomStart(inCellNeighbor(iNeigh1,iAt1))
        nOrb2 = orb%nOrbAtom(inCellNeighbor(iNeigh1,iAt1))
        dm(start2:start2+nOrb2-1, start1:start1+nOrb1-1) = &
            &matmul(eigenvecs(start2:start2+nOrb2-1,1:nLevels), &
            &tmpEigen(1:nLevels,1:nOrb1))
      end do
!$OMP  END PARALLEL DO
    end do

  end subroutine sp_density_matrix_real




  !!* Make a regular density matrix for the complex wave-function case
  !!* @param dm the resulting nOrb*nOrb density matrix with only the elements
  !!* of interest calculated, instead of the whole dm
  !!* @param eigenvecs the eigenvectors of the system
  !!* @param filling the occupation numbers of the orbitals
  !!* @param iNeighbor Neighbor list for each atom (First index from 0!)
  !!* @param nNeighbor Nr. of neighbors for each atom (incl. itself).
  !!* @param orb Informatio about the orbitals.
  !!* @param iAtomStart Atom offset for the squared Hamiltonian
  !!* @param img2CentCell Atomic mapping indexes.
  subroutine sp_density_matrix_cmplx(dm, eigenvecs, filling, iNeighbor, &
      &nNeighbor, orb, iAtomStart, img2CentCell)
    complex(dp), intent(out) :: dm(:,:)
    complex(dp), intent(in) :: eigenvecs(:,:)
    real(dp), intent(in) :: filling(:)
    integer, intent(in) :: iNeighbor(0:,:)
    integer, intent(in) :: nNeighbor(:)
    type(TOrbitals), intent(in) :: orb
    integer, intent(in) :: iAtomStart(:)
    integer, intent(in) :: img2CentCell(:)

    integer                  :: iAt1, iNeigh1, nOrb1, nOrb2, jj, nAtom
    integer                  :: nLevels, start1, start2, ii
    integer, allocatable     :: inCellNeighbor(:,:)
    integer, allocatable     :: nInCellNeighbor(:)
    complex(dp), allocatable :: tmpEigen(:,:)

    @:ASSERT(all(shape(eigenvecs) == shape(dm)))
    @:ASSERT(size(eigenvecs,dim=1) == size(eigenvecs,dim=2))
    @:ASSERT(size(eigenvecs,dim=1) == size(filling))

    allocate(inCellNeighbor(0:size(iNeighbor,dim=1),size(iNeighbor,dim=2)))
    allocate(nInCellNeighbor(size(iNeighbor,dim=2)))

    nAtom = size(orb%nOrbAtom)

    dm(:,:) = cmplx(0.0_dp,0.0_dp,dp)
    inCellNeighbor(:,:) = 0
    nInCellNeighbor(:) = 0

    do ii =  size(filling), 1, -1
      nLevels = ii
      if (abs(filling(ii)) >= epsilon(1.0_dp)) then
        exit
      end if
    end do

    allocate(tmpEigen(nLevels, orb%mOrb))
!$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(iAt1) SCHEDULE(RUNTIME)
    do iAt1 = 1, nAtom
      inCellNeighbor(0:nNeighbor(iAt1),iAt1) = &
          &img2CentCell(iNeighbor(0:nNeighbor(iAt1),iAt1))
      call heap_sort(inCellNeighbor(:nNeighbor(iAt1),iAt1))
      nInCellNeighbor(iAt1) = &
          &unique(inCellNeighbor(:,iAt1), nNeighbor(iAt1)+1) - 1
    end do
!$OMP  END PARALLEL DO

    do iAt1 = 1, nAtom
      nOrb1 = orb%nOrbAtom(iAt1)
      start1 = iAtomStart(iAt1)
      tmpEigen(1:nLevels,1:nOrb1) = &
          &transpose(eigenvecs(start1:start1+nOrb1-1,1:nLevels))

!$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(jj) SCHEDULE(RUNTIME)
      do jj = 1, nLevels
        tmpEigen(jj,1:nOrb1) = filling(jj)*conjg(tmpEigen(jj,1:nOrb1))
      end do
!$OMP  END PARALLEL DO

!$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(iNeigh1,start2,nOrb2) &
!$OMP& SCHEDULE(RUNTIME)
      do iNeigh1 = 0, nInCellNeighbor(iAt1)
        start2 = iAtomStart(inCellNeighbor(iNeigh1,iAt1))
        nOrb2 = orb%nOrbAtom(inCellNeighbor(iNeigh1,iAt1))
        dm(start2:start2+nOrb2-1, start1:start1+nOrb1-1) = &
            &matmul(eigenvecs(start2:start2+nOrb2-1,1:nLevels), &
            &tmpEigen(1:nLevels,1:nOrb1))
      end do
!$OMP  END PARALLEL DO
    end do

  end subroutine sp_density_matrix_cmplx




  !!* Make an energy weighted density matrix for the real wave-function case
  !!* @param dm the resulting nOrb*nOrb density matrix with only the elements
  !!* of interest calculated, instead of the whole dm
  !!* @param eigenvecs the eigenvectors of the system
  !!* @param filling the occupation numbers of the orbitals
  !!* @param eigen eigenvalues of the system
  !!* @param iNeighbor Neighbor list for each atom (First index from 0!)
  !!* @param nNeighbor Nr. of neighbors for each atom (incl. itself).
  !!* @param orb Information about the orbitals.
  !!* @param iAtomStart Atom offset for the squared Hamiltonian
  !!* @param img2CentCell Atomic mapping indexes.
  subroutine sp_energy_density_matrix_real(dm, eigenvecs, filling, eigen, &
      &iNeighbor, nNeighbor, orb, iAtomStart, img2CentCell)
    real(dp), intent(out) :: dm(:,:)
    real(dp), intent(in) :: eigenvecs(:,:)
    real(dp), intent(in) :: filling(:)
    real(dp), intent(in) :: eigen(:)
    integer, intent(in) :: iNeighbor(0:,:)
    integer, intent(in) :: nNeighbor(:)
    type(TOrbitals), intent(in) :: orb
    integer, intent(in) :: iAtomStart(:)
    integer, intent(in) :: img2CentCell(:)

    integer                  :: iAt1, iNeigh1, jj, nOrb1, nOrb2, nAtom
    integer                  :: nLevels, start1, start2, ii
    integer, allocatable     :: inCellNeighbor(:,:)
    integer, allocatable     :: nInCellNeighbor(:)
    real(dp), allocatable    :: tmpEigen(:,:)


    @:ASSERT(all(shape(eigenvecs) == shape(dm)))
    @:ASSERT(size(eigenvecs,dim=1) == size(eigenvecs,dim=2))
    @:ASSERT(size(eigenvecs,dim=1) == size(filling))
    @:ASSERT(size(eigen) == size(filling))

    allocate(inCellNeighbor(0:size(iNeighbor,dim=1),size(iNeighbor,dim=2)))
    allocate(nInCellNeighbor(size(iNeighbor,dim=2)))

    nAtom = size(orb%nOrbAtom)
    dm(:,:) = 0.0_dp

    inCellNeighbor(:,:) = 0
    nInCellNeighbor(:) = 0

    do ii =  size(filling), 1, -1
      nLevels = ii
      if (abs(filling(ii)) >= epsilon(1.0_dp)) then
        exit
      end if
    end do

    allocate(tmpEigen(nLevels,orb%mOrb))
!$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(iAt1) SCHEDULE(RUNTIME)
    do iAt1 = 1, nAtom
      inCellNeighbor(0:nNeighbor(iAt1),iAt1) = &
          & img2CentCell(iNeighbor(0:nNeighbor(iAt1),iAt1))
      call heap_sort(inCellNeighbor(:nNeighbor(iAt1),iAt1))
      nInCellNeighbor(iAt1) = &
          &unique(inCellNeighbor(:,iAt1), nNeighbor(iAt1)+1) - 1
    end do
!$OMP  END PARALLEL DO

    do iAt1 = 1, nAtom
      nOrb1 = orb%nOrbAtom(iAt1)
      start1 = iAtomStart(iAt1)
      tmpEigen(1:nLevels,1:nOrb1) = &
          & transpose(eigenvecs(start1:start1+nOrb1-1,1:nLevels))

!$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(jj) SCHEDULE(RUNTIME)
      do jj = 1, nLevels
        tmpEigen(jj,1:nOrb1) = filling(jj)*eigen(jj)*tmpEigen(jj,1:nOrb1)
      end do
!$OMP  END PARALLEL DO
!$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(iNeigh1,start2,nOrb2) &
!$OMP& SCHEDULE(RUNTIME)
      do iNeigh1 = 0, nInCellNeighbor(iAt1)
        start2 = iAtomStart(inCellNeighbor(iNeigh1,iAt1))
        nOrb2 = orb%nOrbAtom(inCellNeighbor(iNeigh1,iAt1))
        dm(start2:start2+nOrb2-1, start1:start1+nOrb1-1) = &
            &matmul(eigenvecs(start2:start2+nOrb2-1,1:nLevels), &
            &tmpEigen(1:nLevels,1:nOrb1))
      end do
!$OMP  END PARALLEL DO
    end do

  end subroutine sp_energy_density_matrix_real



  !!* Make an energy weighted density matrix for the complex wave-function
  !!* case
  !!* @param dm the resulting nOrb*nOrb density matrix with only the elements
  !!* of interest calculated, instead of the whole dm
  !!* @param eigenvecs the eigenvectors of the system
  !!* @param filling the occupation numbers of the orbitals
  !!* @param iNeighbor Neighbor list for each atom (First index from 0!)
  !!* @param nNeighbor Nr. of neighbors for each atom (incl. itself).
  !!* @param orb Information about the orbitals.
  !!* @param iAtomStart Atom offset for the squared Hamiltonian
  !!* @param img2CentCell Atomic mapping indexes.
  subroutine sp_energy_density_matrix_cmplx(dm, eigenvecs, filling, eigen, &
      &iNeighbor, nNeighbor, orb, iAtomStart, img2CentCell)
    complex(dp), intent(out) :: dm(:,:)
    complex(dp), intent(in) :: eigenvecs(:,:)
    real(dp), intent(in) :: filling(:)
    real(dp), intent(in) :: eigen(:)
    integer, intent(in) :: iNeighbor(0:,:)
    integer, intent(in) :: nNeighbor(:)
    type(TOrbitals), intent(in) :: orb
    integer, intent(in) :: iAtomStart(:)
    integer, intent(in) :: img2CentCell(:)

    integer                  :: iAt1, iNeigh1, jj, nOrb1, nOrb2
    integer                  :: nAtom, nLevels, start1, start2, ii
    integer, allocatable     :: inCellNeighbor(:,:)
    integer, allocatable     :: nInCellNeighbor(:)
    complex(dp), allocatable :: tmpEigen(:,:)

    @:ASSERT(all(shape(eigenvecs) == shape(dm)))
    @:ASSERT(size(eigenvecs,dim=1) == size(eigenvecs,dim=2))
    @:ASSERT(size(eigenvecs,dim=1) == size(filling))
    @:ASSERT(size(eigen) == size(filling))

    allocate(inCellNeighbor(0:size(iNeighbor,dim=1),size(iNeighbor,dim=2)))
    allocate(nInCellNeighbor(size(iNeighbor,dim=2)))

    nAtom = size(orb%nOrbAtom)
    dm(:,:) = cmplx(0.0_dp,0.0_dp,dp)

    inCellNeighbor(:,:) = 0
    nInCellNeighbor(:) = 0

    do ii =  size(filling), 1, -1
      nLevels = ii
      if (abs(filling(ii)) >= epsilon(1.0_dp)) then
        exit
      end if
    end do

    allocate(tmpEigen(nLevels, orb%mOrb))
!$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(iAt1) SCHEDULE(RUNTIME)
    do iAt1 = 1, nAtom
      inCellNeighbor(0:nNeighbor(iAt1),iAt1) = &
          &img2CentCell(iNeighbor(0:nNeighbor(iAt1),iAt1))
      call heap_sort(inCellNeighbor(:nNeighbor(iAt1),iAt1))
      nInCellNeighbor(iAt1) = &
          &unique(inCellNeighbor(:,iAt1), nNeighbor(iAt1)+1) - 1
    end do
!$OMP  END PARALLEL DO

    do iAt1 = 1, nAtom
      nOrb1 = orb%nOrbAtom(iAt1)
      start1 = iAtomStart(iAt1)
      tmpEigen(1:nLevels,1:nOrb1) = &
          & transpose(eigenvecs(start1:start1+nOrb1-1,1:nLevels))

!$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(jj) SCHEDULE(RUNTIME)
      do jj = 1, nLevels
        tmpEigen(jj,1:nOrb1) = filling(jj)*eigen(jj)*conjg(tmpEigen(jj,1:nOrb1))
      end do
!$OMP  END PARALLEL DO
!$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(iNeigh1,start2,nOrb2) &
!$OMP& SCHEDULE(RUNTIME)
      do iNeigh1 = 0, nInCellNeighbor(iAt1)
        start2 = iAtomStart(inCellNeighbor(iNeigh1,iAt1))
        nOrb2 = orb%nOrbAtom(inCellNeighbor(iNeigh1,iAt1))
        dm(start2:start2+nOrb2-1, start1:start1+nOrb1-1) = &
            &matmul(eigenvecs(start2:start2+nOrb2-1,1:nLevels), &
            &tmpEigen(1:nLevels,1:nOrb1))
      end do
!$OMP  END PARALLEL DO
    end do

  end subroutine sp_energy_density_matrix_cmplx


end module densitymatrix
