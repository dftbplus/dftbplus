!--------------------------------------------------------------------------------------------------!
!  DFTB+: general package for performing fast atomistic simulations                                !
!  Copyright (C) 2017  DFTB+ developers group                                                      !
!                                                                                                  !
!  See the LICENSE file for terms of usage and distribution.                                       !
!--------------------------------------------------------------------------------------------------!

#:include 'common.fypp'

!> Contains routines to calculate the coulombic interaction in non periodic and periodic systems.
module coulomb
#:if WITH_SCALAPACK
  use scalapackfx
#:endif
  use assert
  use accuracy
  use message
  use errorfunction
  use constants, only : pi
  implicit none

  private

  public :: invR, sumInvR, addInvRPrime, getOptimalAlphaEwald, getMaxGEwald
  public :: getMaxREwald, ewald, invR_stress
  public :: addInvRPrimeXlbomd
#:if WITH_SCALAPACK
  public :: getInvRClusterBlacs, getInvRPeriodicBlacs
  public :: getDInvRClusterBlacs, getDInvRPeriodicBlacs
  public :: getDInvRXlbomdClusterBlacs, getDInvRXlbomdPeriodicBlacs
#:endif


  !> 1/r interaction for all atoms with another group
  interface sumInvR
    module procedure sumInvR_cluster_asymm
    module procedure sumInvR_periodic_asymm
  end interface sumInvR


  !> 1/r interaction
  interface invR
    module procedure invR_cluster
    module procedure invR_periodic
  end interface invR


  !> 1/r^2
  interface addInvRPrime
    module procedure addInvRPrime_cluster
    module procedure addInvRPrime_cluster_asymm
    module procedure addInvRPrime_periodic
    module procedure addInvRPrime_periodic_asymm
  end interface addInvRPrime


  !> 1/r^2 term for extended lagrangian
  interface addInvRPrimeXlbomd
    module procedure addInvRPrimeXlbomd_cluster
    module procedure addInvRPrimeXlbomd_periodic
  end interface addInvRPrimeXlbomd


  character(len=*), parameter :: ftTooClose = &
      &"('The objects with the following indexes are too close to each other',I5,I5)"

  !> Maximal argument value of erf, after which it is constant
  real(dp), parameter :: erfArgLimit = 10.0_dp

contains


  !> Calculates the 1/R Matrix for all atoms for the non-periodic case.  Only the lower triangle is
  !> constructed.
  subroutine invR_cluster(invRMat, nAtom, coord)

    !> Matrix of 1/R values for each atom pair.
    real(dp), intent(out) :: invRMat(:,:)

    !> number of atoms
    integer, intent(in) :: nAtom

    !> List of atomic coordinates.
    real(dp), intent(in) :: coord(:,:)

    integer :: ii, jj
    real(dp) :: dist, vect(3)

    @:ASSERT(all(shape(invRMat) == (/ nAtom, nAtom /)))
    @:ASSERT(size(coord, dim=1) == 3)
    @:ASSERT(size(coord, dim=2) >= nAtom)

    invRMat(:,:) = 0.0_dp
    !$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(ii,jj,vect,dist) SCHEDULE(RUNTIME)
    do ii = 1, nAtom
      do jj = ii + 1, nAtom
        vect(:) = coord(:,ii) - coord(:,jj)
        dist = sqrt(sum(vect(:)**2))
        invRMat(jj,ii) = 1.0_dp/dist
      end do
    end do
    !$OMP  END PARALLEL DO

  end subroutine invR_cluster

#:if WITH_SCALAPACK

  !> Calculates the 1/R Matrix for all atoms for the non-periodic case.  Only the lower triangle is
  !> constructed.
  subroutine getInvRClusterBlacs(grid, coord, descInvRMat, invRMat)

    !> BLACS grid involved in 1/R calculation
    type(blacsgrid), intent(in) :: grid

    !> List of atomic coordinates.
    real(dp), intent(in) :: coord(:,:)

    !> Matrix descriptor for 1/R matrix
    integer, intent(in) :: descInvRMat(DLEN_)

    !> Matrix of 1/R values for each atom pair.
    real(dp), intent(out) :: invRMat(:,:)

    integer :: ii, jj, iAt2, iAt1
    real(dp) :: dist, vect(3)

    invRMat(:,:) = 0.0_dp
    do jj = 1, size(invRMat, dim=2)
      iAt1 = scalafx_indxl2g(jj, descInvRMat(NB_), grid%mycol, descInvRMat(CSRC_), grid%ncol)
      do ii = 1, size(invRMat, dim=1)
        iAt2 = scalafx_indxl2g(ii, descInvRMat(MB_), grid%myrow, descInvRMat(RSRC_), grid%nrow)
        if (iAt2 <= iAt1) then
          cycle
        end if
        vect(:) = coord(:,iAt1) - coord(:,iAt2)
        dist = sqrt(sum(vect**2))
        invRMat(ii, jj) = 1.0_dp / dist
      end do
    end do

  end subroutine getInvRClusterBlacs

#:endif


  !> Calculates the summed 1/R vector for all atoms for the non-periodic case asymmmetric case (like
  !> interaction of atoms with point charges).
  subroutine sumInvR_cluster_asymm(invRVec, nAtom0, nAtom1, coord0, coord1, charges1, blurWidths1)

    !> Vector of sum_i q_i/|R_atom - R_i] values for each atom
    real(dp), intent(out) :: invRVec(:)

    !> Number of atoms in the first group
    integer, intent(in) :: nAtom0

    !> Number of atoms in the second group
    integer, intent(in) :: nAtom1

    !> Coordinates of the first group of objects (atoms)
    real(dp), intent(in) :: coord0(:,:)

    !> Coordinates of the 2nd group of objects (point charges)
    real(dp), intent(in) :: coord1(:,:)

    !> Charges of the 2nd group of objects
    real(dp), intent(in) :: charges1(:)

    !> Gaussian blur widht of the charges in the 2nd group
    real(dp), intent(in), optional :: blurWidths1(:)

    integer :: iAt0, iAt1
    real(dp) :: dist, vect(3), fTmp
    character(len=100) :: errorString

    @:ASSERT(size(invRVec) == nAtom0)
    @:ASSERT(size(coord0, dim=2) >= nAtom0)
    @:ASSERT(size(coord0, dim=1) == 3)
    @:ASSERT(size(coord1, dim=2) >= nAtom1)
    @:ASSERT(size(coord1, dim=1) == 3)
    @:ASSERT(size(charges1) == nAtom1)
#:call ASSERT_CODE
    if (present(blurWidths1)) then
      @:ASSERT(size(blurWidths1) == nAtom1)
    end if
#:endcall ASSERT_CODE

    ! Doing blurring and non blurring case separately in order to avoid
    ! the if branch deep in the loop
    invRVec(:) = 0.0_dp
    if (present(blurWidths1)) then
      !$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(iAt0,iAt1,vect,dist,fTmp,errorString) &
      !$OMP& SCHEDULE(RUNTIME)
      do iAt0 = 1, nAtom0
        do iAt1 = 1, nAtom1
          vect(:) = coord0(:,iAt0) - coord1(:,iAt1)
          dist = sqrt(sum(vect(:)**2))
          if (dist > epsilon(0.0_dp)) then
            fTmp = charges1(iAt1) / dist
            if (dist < erfArgLimit * blurWidths1(iAt1)) then
              fTmp = fTmp * erfwrap(dist/blurWidths1(iAt1))
            end if
            invRVec(iAt0) = invRVec(iAt0) + fTmp
          else
            write (errorString, ftTooClose) iAt0, iAt1
            call error(trim(errorString))
          end if
        end do
      end do
      !$OMP  END PARALLEL DO
    else
      !$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(iAt0,iAt1,vect,dist,errorString) &
      !$OMP& SCHEDULE(RUNTIME)
      do iAt0 = 1, nAtom0
        do iAt1 = 1, nAtom1
          vect(:) = coord0(:,iAt0) - coord1(:,iAt1)
          dist = sqrt(sum(vect(:)**2))
          if (dist > epsilon(0.0_dp)) then
            invRVec(iAt0) = invRVec(iAt0) + charges1(iAt1) / dist
          else
            write (errorString, ftTooClose) iAt0, iAt1
            call error(trim(errorString))
          end if
        end do
      end do
      !$OMP  END PARALLEL DO
    end if

  end subroutine sumInvR_cluster_asymm


  !> Calculates the 1/R Matrix for all atoms for the periodic case.  Only the lower triangle is
  !> constructed.
  subroutine invR_periodic(invRMat, nAtom, coord, nNeighborEwald, iNeighbor, img2CentCell,&
      & recPoint, alpha, volume)

    !> Matrix of 1/R values for each atom pair.
    real(dp), intent(out) :: invRMat(:,:)

    !> Number of atoms.
    integer, intent(in) :: nAtom

    !> List of atomic coordinates (all atoms).
    real(dp), intent(in) :: coord(:,:)

    !> Nr. of neighbors for each atom for real part of Ewald.
    integer, intent(in) :: nNeighborEwald(:)

    !> List of neighbors for the real space part of Ewald.
    integer, intent(in) :: iNeighbor(0:,:)

    !> Image index for each atom in the central cell.
    integer, intent(in) :: img2CentCell(:)

    !> Contains the points included in the reciprocal sum.  The set should not include the origin or
    !> inversion related points.
    real(dp), intent(in) :: recPoint(:,:)

    !> Parameter for Ewald summation.
    real(dp), intent(in) :: alpha

    !> Volume of the real space unit cell.
    real(dp), intent(in) :: volume

    integer :: iAtom1, iAtom2, iAtom2f, iNeigh

    @:ASSERT(all(shape(invRMat) == (/ nAtom, nAtom /)))
    @:ASSERT(size(coord, dim=1) == 3)
    @:ASSERT(size(coord, dim=2) >= nAtom)
    @:ASSERT(size(nNeighborEwald) == nAtom)
    @:ASSERT(size(iNeighbor, dim=2) == nAtom)
    @:ASSERT(volume > 0.0_dp)

    invRMat(:,:) = 0.0_dp

    ! Real space part of the Ewald sum.
    !$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(iAtom1,iNeigh,iAtom2,iAtom2f) SCHEDULE(RUNTIME)
    do iAtom1 = 1, nAtom
      do iNeigh = 1, nNeighborEwald(iAtom1)
        iAtom2 = iNeighbor(iNeigh, iAtom1)
        iAtom2f = img2CentCell(iAtom2)
        invRMat(iAtom2f, iAtom1) = invRMat(iAtom2f, iAtom1)&
            & + rTerm(sqrt(sum((coord(:,iAtom1)-coord(:,iAtom2))**2)), alpha)
      end do
    end do
    !$OMP  END PARALLEL DO

    ! Reciprocal space part of the Ewald sum.
    !$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(iAtom1,iAtom2) SCHEDULE(RUNTIME)
    do iAtom1 = 1, nAtom
      do iAtom2 = iAtom1, nAtom
        invRMat(iAtom2, iAtom1) = invRMat(iAtom2, iAtom1)&
            & + ewaldReciprocal(coord(:,iAtom1)-coord(:,iAtom2), recPoint, alpha, volume)&
            & - pi / (volume * alpha**2)
      end do
    end do
    !$OMP  END PARALLEL DO

    ! Extra contribution for self interaction.
    !$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(iAtom1) SCHEDULE(RUNTIME)
    do iAtom1 = 1, nAtom
      invRMat(iAtom1, iAtom1) = invRMat(iAtom1, iAtom1) - 2.0_dp * alpha / sqrt(pi)
    end do
    !$OMP  END PARALLEL DO

  end subroutine invR_periodic


#:if WITH_SCALAPACK

  !> Calculates the 1/R Matrix for all atoms for the periodic case.  Only the lower triangle is
  !> constructed.
  subroutine getInvRPeriodicBlacs(grid, coord, nNeighborEwald, iNeighbor, img2CentCell, recPoint,&
      & alpha, volume, descInvRMat, invRMat)

    !> BLACS grid involved in 1/R calculation
    type(blacsgrid), intent(in) :: grid

    !> List of atomic coordinates (all atoms).
    real(dp), intent(in) :: coord(:,:)

    !> Nr. of neighbors for each atom for real part of Ewald.
    integer, intent(in) :: nNeighborEwald(:)

    !> List of neighbors for the real space part of Ewald.
    integer, intent(in) :: iNeighbor(0:,:)

    !> Image index for each atom in the central cell.
    integer, intent(in) :: img2CentCell(:)

    !> Contains the points included in the reciprocal sum.  The set should not include the origin or
    !> inversion related points.
    real(dp), intent(in) :: recPoint(:,:)

    !> Parameter for Ewald summation.
    real(dp), intent(in) :: alpha

    !> Volume of the real space unit cell.
    real(dp), intent(in) :: volume

    !> Descriptor for the 1/R matrix
    integer, intent(in) :: descInvRMat(DLEN_)

    !> Matrix of 1/R values for each atom pair.
    real(dp), intent(out) :: invRMat(:,:)

    integer :: iAt1, iAt2, iAt2f, iNeigh, jj, ii, iLoc, jLoc
    logical :: tLocal

    invRMat(:,:) = 0.0_dp

    ! Real space part of the Ewald sum.
    do jj = 1, size(invRMat, dim=2)
      iAt1 = scalafx_indxl2g(jj, descInvRMat(NB_), grid%mycol, descInvRMat(CSRC_), grid%ncol)
      do iNeigh = 1, nNeighborEwald(iAt1)
        iAt2 = iNeighbor(iNeigh, iAt1)
        iAt2f = img2CentCell(iAt2)
        call scalafx_islocal(grid, descInvRMat, iAt2f, iAt1, tLocal, iLoc, jLoc)
        if (tLocal) then
          invRMat(iLoc, jLoc) = invRMat(iLoc, jLoc)&
              &  + rTerm(sqrt(sum((coord(:,iAt1) - coord(:,iAt2))**2)), alpha)
        end if
      end do
    end do

    ! Reciprocal space part of the Ewald sum.
    do jj = 1, size(invRMat, dim=2)
      iAt1 = scalafx_indxl2g(jj, descInvRMat(NB_), grid%mycol, descInvRMat(CSRC_), grid%ncol)
      do ii = 1, size(invRMat, dim=1)
        iAt2 = scalafx_indxl2g(ii, descInvRMat(MB_), grid%myrow, descInvRMat(RSRC_),&
            & grid%nrow)
        if (iAt2 < iAt1) then
          cycle
        end if
        invRMat(ii, jj) = invRMat(ii, jj)&
            & + ewaldReciprocal(coord(:,iAt1) - coord(:,iAt2), recPoint, alpha, volume)&
            & - pi / (volume * alpha**2)
      end do
    end do

    ! Extra contribution for self interaction.
    do jj = 1, size(invRMat, dim=2)
      iAt1 = scalafx_indxl2g(jj, descInvRMat(NB_), grid%mycol, descInvRMat(CSRC_), grid%ncol)
      call scalafx_islocal(grid, descInvRMat, iAt1, iAt1, tLocal, iLoc, jLoc)
      if (tLocal) then
        invRMat(iLoc, jLoc) = invRMat(iLoc, jLoc) - 2.0_dp * alpha / sqrt(pi)
      end if
    end do

  end subroutine getInvRPeriodicBlacs

#:endif

  !> Calculates summed 1/R vector for two groups of objects for the periodic case.
  subroutine sumInvR_periodic_asymm(invRVec, nAtom0, nAtom1, coord0, coord1, charges1, rLat, gLat,&
      & alpha, volume)

    !> Vector of sum_i q_i/|R_atom - R_i] values for each atom
    real(dp), intent(out) :: invRVec(:)

    !> Number of atoms in the first group
    integer, intent(in) :: nAtom0

    !> Number of atoms in the second group
    integer, intent(in) :: nAtom1

    !> Coordinates of the first group of objects (atoms)
    real(dp), intent(in) :: coord0(:,:)

    !> Coordinates of the 2nd group of objects (point charges)
    real(dp), intent(in) :: coord1(:,:)

    !> Charges of the 2nd group of objects
    real(dp), intent(in) :: charges1(:)

    !> Lattice vectors to be used for the real Ewald summation
    real(dp), intent(in) :: rLat(:,:)

    !> Lattice vectors to be used for the reciprocal Ewald summation.
    real(dp), intent(in) :: gLat(:,:)

    !> Parameter of the Ewald summation
    real(dp), intent(in) :: alpha

    !> Volume of the supercell.
    real(dp), intent(in) :: volume

    integer :: iAt0, iAt1
    real(dp) :: rTmp, rr(3)

    @:ASSERT(size(invRVec) == nAtom0)
    @:ASSERT(size(coord0, dim=2) >= nAtom0)
    @:ASSERT(size(coord0, dim=1) == 3)
    @:ASSERT(size(coord1, dim=2) >= nAtom1)
    @:ASSERT(size(coord1, dim=1) == 3)
    @:ASSERT(size(charges1) == nAtom1)
    @:ASSERT(size(rLat, dim=1) == 3)
    @:ASSERT(size(gLat, dim=1) == 3)
    @:ASSERT(volume > 0.0_dp)

    invRVec(:) = 0.0_dp
    !$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(iAt0,iAt1,rr,rTmp) SCHEDULE(RUNTIME)
    do iAt0 = 1, nAtom0
      do iAt1 = 1, nAtom1
        rr = coord0(:,iAt0) - coord1(:,iAt1)

        rTmp = ewald(rr, rLat, gLat, alpha, volume)
        invRVec(iAt0) = invRVec(iAt0) + rTmp * charges1(iAt1)
      end do
    end do
    !$OMP  END PARALLEL DO

  end subroutine sumInvR_periodic_asymm


  !> Calculates the -1/R**2 deriv contribution for all atoms for the non-periodic case, without
  !> storing anything.
  subroutine addInvRPrime_cluster(deriv, nAtom, coord, deltaQAtom)

    !> Contains the derivative on exit.
    real(dp), intent(inout) :: deriv(:,:)

    !> Number of atoms.
    integer, intent(in) :: nAtom

    !> List of atomic coordinates.
    real(dp), intent(in) :: coord(:,:)

    !> List of charges on each atom.
    real(dp), intent(in) :: deltaQAtom(:)

    integer :: ii, jj
    real(dp) :: dist, vect(3), fTmp

    @:ASSERT(size(deriv, dim=1) == 3)
    @:ASSERT(size(deriv, dim=2) >= nAtom)
    @:ASSERT(size(coord, dim=1) == 3)
    @:ASSERT(size(coord, dim=2) >= nAtom)
    @:ASSERT(size(deltaQAtom) == nAtom)

    !$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(ii,jj,vect,dist,ftmp) &
    !$OMP& SCHEDULE(RUNTIME) REDUCTION(+:deriv)
    do ii = 1, nAtom
      do jj = ii + 1, nAtom
        vect(:) = coord(:,ii) - coord(:,jj)
        dist = sqrt(sum(vect(:)**2))
        fTmp = -deltaQAtom(ii) * deltaQAtom(jj) / (dist**3)
        deriv(:,ii) = deriv(:,ii) + vect(:)*fTmp
        ! Skew-symmetric 1/r2 interaction, so the other triangle is calculated :
        deriv(:,jj) = deriv(:,jj) - vect(:)*fTmp
      end do
    end do
    !$OMP  END PARALLEL DO

  end subroutine addInvRPrime_cluster

#:if WITH_SCALAPACK

  !> Calculates the -1/R**2 deriv contribution for all atoms for the non-periodic case, without
  !> storing anything.
  subroutine getDInvRClusterBlacs(grid, descAtomSqr, localShape, coord, deltaQAtom, deriv)

    !> BLACS grid of the distributed derivative vector
    type(blacsgrid), intent(in) :: grid

    !> Descriptor for the nAtom x nAtom electrostatic matrix distributed on the grid
    integer, intent(in) :: descAtomSqr(DLEN_)

    !> Local shape of the distributed nAtom x nAtom electrostatic matrix
    integer, intent(in) :: localShape(:)

    !> List of atomic coordinates.
    real(dp), intent(in) :: coord(:,:)

    !> List of charges on each atom.
    real(dp), intent(in) :: deltaQAtom(:)

    !> Contains the derivative on exit.
    real(dp), intent(out) :: deriv(:,:)

    integer :: ii, jj, iAt1, iAt2
    real(dp) :: dist, vect(3), fTmp

    deriv(:,:) = 0.0_dp
    do jj = 1, localShape(2)
      iAt1 = scalafx_indxl2g(jj, descAtomSqr(NB_), grid%mycol, descAtomSqr(CSRC_), grid%ncol)
      do ii = 1, localShape(1)
        iAt2 = scalafx_indxl2g(ii, descAtomSqr(MB_), grid%myrow, descAtomSqr(RSRC_), grid%nrow)
        if (iAt1 /= iAt2) then
          vect(:) = coord(:,iAt1) - coord(:,iAt2)
          dist = sqrt(sum(vect**2))
          fTmp = -deltaQAtom(iAt1) * deltaQAtom(iAt2) / (dist**3)
          deriv(:,iAt1) = deriv(:,iAt1) + vect * fTmp
        end if
      end do
    end do

  end subroutine getDInvRClusterBlacs

#:endif


  !> Calculates the -1/R**2 deriv contribution for extended lagrangian dynamics forces in a periodic
  !> geometry
  subroutine addInvRPrimeXlbomd_cluster(nAtom, coord, dQInAtom, dQOutAtom, deriv)

    !> number of atoms
    integer, intent(in) :: nAtom

    !> coordinates of atoms
    real(dp), intent(in) :: coord(:,:)

    !> input charge fluctuations
    real(dp), intent(in) :: dQInAtom(:)

    !> output charge fluctuations
    real(dp), intent(in) :: dQOutAtom(:)

    !> energy derivative to add contribution to
    real(dp), intent(inout) :: deriv(:,:)

    integer :: iAt1, iAt2
    real(dp) :: dist, vect(3), fTmp, prefac

    @:ASSERT(size(deriv, dim=1) == 3)
    @:ASSERT(size(deriv, dim=2) >= nAtom)
    @:ASSERT(size(coord, dim=1) == 3)
    @:ASSERT(size(coord, dim=2) >= nAtom)
    @:ASSERT(size(dQInAtom) == nAtom)
    @:ASSERT(size(dQOutAtom) == nAtom)

    !$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(iAt1,iAt2,vect,dist,prefac,ftmp) &
    !$OMP& SCHEDULE(RUNTIME) REDUCTION(+:deriv)
    do iAt1 = 1, nAtom
      do iAt2 = iAt1 + 1, nAtom
        vect(:) = coord(:,iAt1) - coord(:,iAt2)
        dist = sqrt(sum(vect(:)**2))
        prefac = dQOutAtom(iAt1) * dQInAtom(iAt2) + dQInAtom(iAt1) * dQOutAtom(iAt2) &
            & - dQInAtom(iAt1) * dQInAtom(iAt2)
        fTmp = -prefac / (dist**3)
        deriv(:,iAt1) = deriv(:,iAt1) + vect * fTmp
        ! Skew-symmetric 1/r2 interaction, so the other triangle is calculated
        deriv(:,iAt2) = deriv(:,iAt2) - vect *fTmp
      end do
    end do
    !$OMP  END PARALLEL DO

  end subroutine addInvRPrimeXlbomd_cluster


#:if WITH_SCALAPACK

  !> Calculates the -1/R**2 deriv contribution for extended lagrangian dynamics forces in a periodic
  !> geometry
  subroutine getDInvRXlbomdClusterBlacs(grid, descAtomSqr, localShape, coord, dQInAtom, dQOutAtom,&
      & deriv)

    !> BLACS grid of the distributed derivative vector
    type(blacsgrid), intent(in) :: grid

    !> Descriptor for an nAtom x nAtom matrix distributed on the grid
    integer, intent(in) :: descAtomSqr(DLEN_)

    !> Local shape of the distributed nAtom x nAtom matrix
    integer, intent(in) :: localShape(:)

    !> coordinates of atoms
    real(dp), intent(in) :: coord(:,:)

    !> input charge fluctuations
    real(dp), intent(in) :: dQInAtom(:)

    !> output charge fluctuations
    real(dp), intent(in) :: dQOutAtom(:)

    !> energy derivative to add contribution to
    real(dp), intent(out) :: deriv(:,:)

    integer :: ii, jj, iAt1, iAt2
    real(dp) :: dist, vect(3), fTmp, prefac

    deriv(:,:) = 0.0_dp
    do jj = 1, localShape(2)
      iAt1 = scalafx_indxl2g(jj, descAtomSqr(NB_), grid%mycol, descAtomSqr(CSRC_), grid%ncol)
      do ii = 1, localShape(1)
        iAt2 = scalafx_indxl2g(ii, descAtomSqr(MB_), grid%myrow, descAtomSqr(RSRC_), grid%nrow)
        if (iAt1 /= iAt2) then
          vect(:) = coord(:,iAt1) - coord(:,iAt2)
          dist = sqrt(sum(vect**2))
          prefac = dQOutAtom(iAt1) * dQInAtom(iAt2) + dQInAtom(iAt1) * dQOutAtom(iAt2) &
              & - dQInAtom(iAt1) * dQInAtom(iAt2)
          fTmp = -prefac / (dist**3)
          deriv(:,iAt1) = deriv(:,iAt1) + vect * fTmp
        end if
      end do
    end do

  end subroutine getDInvRXlbomdClusterBlacs

#:endif


  !> Calculates the -1/R**2 deriv contribution for charged atoms interacting with a group of charged
  !> objects (like point charges) for the non-periodic case, without storing anything.
  subroutine addInvRPrime_cluster_asymm(deriv0, deriv1, nAtom0, nAtom1, coord0, coord1, charge0,&
      & charge1, blurWidths1)

    !> Contains the derivative for the first group
    real(dp), intent(inout) :: deriv0(:,:)

    !> Contains the derivative for the second group
    real(dp), intent(inout) :: deriv1(:,:)

    !> Number of atoms in the first group
    integer, intent(in) :: nAtom0

    !> Number of atoms in the second group
    integer, intent(in) :: nAtom1

    !> List of atomic coordinates.
    real(dp), intent(in) :: coord0(:,:)

    !> List of the point charge coordinates
    real(dp), intent(in) :: coord1(:,:)

    !> Charge of the atoms.
    real(dp), intent(in) :: charge0(:)

    !> Charge of the point charges.
    real(dp), intent(in) :: charge1(:)

    !> if gaussian distribution for the charge
    real(dp), intent(in), optional :: blurWidths1(:)

    integer :: iAt0, iAt1
    real(dp) :: dist, vect(3), fTmp(3), sigma, rs

    @:ASSERT(size(deriv0, dim=1) == 3)
    @:ASSERT(size(deriv0, dim=2) >= nAtom0)
    @:ASSERT(size(deriv1, dim=1) == 3)
    @:ASSERT(size(deriv1, dim=2) >= nAtom1)
    @:ASSERT(size(coord0, dim=1) == 3)
    @:ASSERT(size(coord0, dim=2) >= nAtom0)
    @:ASSERT(size(coord1, dim=1) == 3)
    @:ASSERT(size(coord1, dim=2) >= nAtom1)
    @:ASSERT(size(charge0) == nAtom0)
    @:ASSERT(size(charge1) == nAtom1)
#:call ASSERT_CODE
    if (present(blurWidths1)) then
      @:ASSERT(size(blurWidths1) == nAtom1)
    end if
#:endcall ASSERT_CODE

    ! Doing blured and unblured cases separately to avoid ifs in the loop
    if (present(blurWidths1)) then
      !$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(iAt0,iAt1,vect,dist,ftmp,sigma,rs) &
      !$OMP& SCHEDULE(RUNTIME) REDUCTION(+:deriv0,deriv1)
      do iAt0 = 1, nAtom0
        do iAt1 = 1, nAtom1
          vect(:) = coord0(:,iAt0) - coord1(:,iAt1)
          dist = sqrt(sum(vect(:)**2))
          fTmp = -vect(:) / (dist**3)
          if (dist < erfArgLimit * blurWidths1(iAt1)) then
            sigma = blurWidths1(iAt1)
            rs = dist / sigma
            fTmp = fTmp * (erfwrap(rs) - 2.0_dp/(sqrt(pi)*sigma) * dist * exp(-(rs**2)))
          end if
          fTmp = charge0(iAt0) * charge1(iAt1) * fTmp
          deriv0(:,iAt0) = deriv0(:,iAt0) + fTmp(:)
          deriv1(:,iAt1) = deriv1(:,iAt1) - fTmp(:)
        end do
      end do
      !$OMP END PARALLEL DO
    else
      !$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(iAt0,iAt1,vect,dist,ftmp) &
      !$OMP& SCHEDULE(RUNTIME) REDUCTION(+:deriv0,deriv1)
      do iAt0 = 1, nAtom0
        do iAt1 = 1, nAtom1
          vect(:) = coord0(:,iAt0) - coord1(:,iAt1)
          dist = sqrt(sum(vect(:)**2))
          fTmp = -charge0(iAt0) * charge1(iAt1) / (dist**3) * vect(:)
          deriv0(:,iAt0) = deriv0(:,iAt0) + fTmp(:)
          deriv1(:,iAt1) = deriv1(:,iAt1) - fTmp(:)
        end do
      end do
      !$OMP  END PARALLEL DO
    end if

  end subroutine addInvRPrime_cluster_asymm


  !> Calculates the -1/R**2 deriv contribution for the periodic case, without storing anything.
  subroutine addInvRPrime_periodic(deriv, nAtom, coord, nNeighborEwald, iNeighbor, img2CentCell,&
      & recPoint, alpha, volume,deltaQAtom)

    !> Derivative on exit
    real(dp), intent(inout) :: deriv(:,:)

    !> Number of atoms
    integer, intent(in) :: nAtom

    !> List of atomic coordinates (all atoms).
    real(dp), intent(in) :: coord(:,:)

    !> Nr. of neighbors for each atom for real part ofEwald.
    integer, intent(in) :: nNeighborEwald(:)

    !> list of neighbours for each atom
    integer, intent(in) :: iNeighbor(0:,:)

    !> mapping from image atoms to central cell
    integer, intent(in) :: img2CentCell(:)

    !> Contains the points included in the reciprocal sum. The set should not include the origin or
    !> inversion related points.
    real(dp), intent(in) :: recPoint(:,:)

    !> Parameter for Ewald summation.
    real(dp), intent(in) :: alpha

    !> Volume of the real space unit cell.
    real(dp), intent(in) :: volume

    !> List of charges on each atom
    real(dp), intent(in) :: deltaQAtom(:)

    integer :: iAtom1, iAtom2, iAtom2f, iNeigh
    real(dp) :: r(3)

    @:ASSERT(size(deriv, dim=1) == 3)
    @:ASSERT(size(deriv, dim=2) >= nAtom)
    @:ASSERT(size(coord, dim=1) == 3)
    @:ASSERT(size(coord, dim=2) >= nAtom)
    @:ASSERT(size(nNeighborEwald) == nAtom)
    @:ASSERT(size(iNeighbor, dim=2) == nAtom)
    @:ASSERT(volume > 0.0_dp)
    @:ASSERT(size(deltaQAtom) == nAtom)

    ! d(1/R)/dr real space
    !$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(iAtom1,iNeigh,iAtom2,iAtom2f,r) &
    !$OMP& SCHEDULE(RUNTIME) REDUCTION(+:deriv)
    do iAtom1 = 1, nAtom
      do iNeigh = 1, nNeighborEwald(iAtom1)
        iAtom2 = iNeighbor(iNeigh, iAtom1)
        iAtom2f = img2CentCell(iAtom2)
        if (iAtom2f /= iAtom1) then
          r(:) = coord(:,iAtom1)-coord(:,iAtom2)
          deriv(:,iAtom1) = deriv(:,iAtom1)&
              & + derivRTerm(r,alpha) * deltaQAtom(iAtom1) * deltaQAtom(iAtom2f)
          deriv(:,iAtom2f) = deriv(:,iAtom2f)&
              & - derivRTerm(r,alpha) * deltaQAtom(iAtom1) * deltaQAtom(iAtom2f)
        end if
      end do
    end do
    !$OMP  END PARALLEL DO

    ! d(1/R)/dr reciprocal space
    !$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(iAtom1,iAtom2,r) &
    !$OMP& SCHEDULE(RUNTIME) REDUCTION(+:deriv)
    do iAtom1 = 1, nAtom
      do iAtom2 = iAtom1+1, nAtom
        r(:) = coord(:,iAtom1)-coord(:,iAtom2)
        deriv(:,iAtom1) = deriv(:,iAtom1)&
            & + derivEwaldReciprocal(r,recPoint,alpha,volume)*deltaQAtom(iAtom1)*deltaQAtom(iAtom2)
        deriv(:,iAtom2) = deriv(:,iAtom2)&
            & - derivEwaldReciprocal(r,recPoint,alpha,volume)*deltaQAtom(iAtom1)*deltaQAtom(iAtom2)
      end do
    end do
    !$OMP  END PARALLEL DO

  end subroutine addInvRPrime_periodic


#:if WITH_SCALAPACK

  subroutine getDInvRPeriodicBlacs(grid, descAtomSqr, localShape, coord, nNeighborEwald,&
      & iNeighbor, img2CentCell, recPoint, alpha, volume, deltaQAtom, deriv)

    !> BLACS grid of the distributed derivative vector
    type(blacsgrid), intent(in) :: grid

    !> Descriptor for an nAtom x nAtom matrix distributed on the grid
    integer, intent(in) :: descAtomSqr(DLEN_)

    !> Local shape of the distributed nAtom x nAtom matrix
    integer, intent(in) :: localShape(:)

    !> List of atomic coordinates (all atoms).
    real(dp), intent(in) :: coord(:,:)

    !> Nr. of neighbors for each atom for real part ofEwald.
    integer, intent(in) :: nNeighborEwald(:)

    !> list of neighbours for each atom
    integer, intent(in) :: iNeighbor(0:,:)

    !> mapping from image atoms to central cell
    integer, intent(in) :: img2CentCell(:)

    !> Contains the points included in the reciprocal sum. The set should not include the origin or
    !> inversion related points.
    real(dp), intent(in) :: recPoint(:,:)

    !> Parameter for Ewald summation.
    real(dp), intent(in) :: alpha

    !> Volume of the real space unit cell.
    real(dp), intent(in) :: volume

    !> List of charges on each atom
    real(dp), intent(in) :: deltaQAtom(:)

    !> Derivative on exit
    real(dp), intent(out) :: deriv(:,:)

    integer :: ii, jj, iAt1, iAt2, iAt2f, iNeigh, iLoc, jLoc
    real(dp) :: rr(3), contrib(3)
    logical :: tLocal

    deriv(:,:) = 0.0_dp

    do jj = 1, localShape(2)
      iAt1 = scalafx_indxl2g(jj, descAtomSqr(NB_), grid%mycol, descAtomSqr(CSRC_), grid%ncol)
      do iNeigh = 1, nNeighborEwald(iAt1)
        iAt2 = iNeighbor(iNeigh, iAt1)
        iAt2f = img2CentCell(iAt2)
        call scalafx_islocal(grid, descAtomSqr, iAt2f, iAt1, tLocal, iLoc, jLoc)
        if (tLocal .and. iAt2f /= iAt1) then
          rr(:) = coord(:,iAt1) - coord(:,iAt2)
          contrib = derivRTerm(rr, alpha) * deltaQAtom(iAt1) * deltaQAtom(iAt2f)
          deriv(:,iAt1) = deriv(:,iAt1) + contrib
          ! Neighbor list only contains lower triange: add also to screw symmetric equivalent
          deriv(:,iAt2f) = deriv(:,iAt2f) - contrib
        end if
      end do
    end do

    do jj = 1, localShape(2)
      iAt1 = scalafx_indxl2g(jj, descAtomSqr(NB_), grid%mycol, descAtomSqr(CSRC_), grid%ncol)
      do ii = 1, localShape(1)
        iAt2 = scalafx_indxl2g(ii, descAtomSqr(MB_), grid%myrow, descAtomSqr(RSRC_), grid%nrow)
        if (iAt2 /= iAt1) then
          rr(:) = coord(:,iAt1) - coord(:,iAt2)
          deriv(:,iAt1) = deriv(:,iAt1)&
              & + derivEwaldReciprocal(rr, recPoint, alpha, volume) * deltaQAtom(iAt1)&
              & * deltaQAtom(iAt2)
        end if
      end do
    end do

  end subroutine getDInvRPeriodicBlacs

#:endif


  !> Calculates the -1/R**2 deriv contribution for extended lagrangian dynamics forces
  subroutine addInvRPrimeXlbomd_periodic(nAtom, coord, nNeighborEwald, iNeighbor, img2CentCell,&
      & recPoint, alpha, volume, dQInAtom, dQOutAtom, deriv)

    !> number of atoms
    integer, intent(in) :: nAtom

    !> coordinates of atoms
    real(dp), intent(in) :: coord(:,:)

    !> Nr. of neighbors for each atom for real part of Ewald
    integer, intent(in) :: nNeighborEwald(:)

    !> List of neighbors for the real space part of Ewald.
    integer, intent(in) :: iNeighbor(0:,:)

    !> Image of each atom in the central cell.
    integer, intent(in) :: img2CentCell(:)

    !> Contains the points included in the reciprocal sum. The set should not include the origin or
    !> inversion related points.
    real(dp), intent(in) :: recPoint(:,:)

    !> Ewald parameter
    real(dp), intent(in) :: alpha

    !> cell volume
    real(dp), intent(in) :: volume

    !> input charge fluctuations
    real(dp), intent(in) :: dQInAtom(:)

    !> output charge fluctuations
    real(dp), intent(in) :: dQOutAtom(:)

    !> energy derivative to add contribution to
    real(dp), intent(inout) :: deriv(:,:)

    integer :: iAt1, iAt2, iAt2f, iNeigh
    real(dp) :: rr(3), contrib(3), prefac

    @:ASSERT(size(deriv, dim=1) == 3)
    @:ASSERT(size(deriv, dim=2) >= nAtom)
    @:ASSERT(size(coord, dim=1) == 3)
    @:ASSERT(size(coord, dim=2) >= nAtom)
    @:ASSERT(size(nNeighborEwald) == nAtom)
    @:ASSERT(size(iNeighbor, dim=2) == nAtom)
    @:ASSERT(volume > 0.0_dp)
    @:ASSERT(size(dQOutAtom) == nAtom)
    @:ASSERT(size(dQInAtom) == nAtom)

    ! real space
    !$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(iAt1,iNeigh,iAt2,iAt2f,rr,prefac,contrib) &
    !$OMP& SCHEDULE(RUNTIME) REDUCTION(+:deriv)
    do iAt1 = 1, nAtom
      do iNeigh = 1, nNeighborEwald(iAt1)
        iAt2 = iNeighbor(iNeigh, iAt1)
        iAt2f = img2CentCell(iAt2)
        if (iAt2f == iAt1) then
          cycle
        end if
        rr(:) = coord(:,iAt1) - coord(:,iAt2)
        prefac = dQOutAtom(iAt1) * dQInAtom(iAt2f) + dQInAtom(iAt1) * dQOutAtom(iAt2f)&
            & - dQInAtom(iAt1) * dQInAtom(iAt2f)
        contrib(:) = prefac * derivRTerm(rr, alpha)
        deriv(:,iAt1) = deriv(:,iAt1) + contrib
        deriv(:,iAt2f) = deriv(:,iAt2f) - contrib
      end do
    end do
    !$OMP  END PARALLEL DO

    ! reciprocal space
    !$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(iAt1,iAt2,rr,prefac,contrib) &
    !$OMP& SCHEDULE(RUNTIME) REDUCTION(+:deriv)
    do iAt1 = 1, nAtom
      do iAt2 = iAt1 + 1, nAtom
        rr(:) = coord(:,iAt1) - coord(:,iAt2)
        prefac = dQOutAtom(iAt1) * dQInAtom(iAt2) + dQInAtom(iAt1) * dQOutAtom(iAt2)&
            & - dQInAtom(iAt1) * dQInAtom(iAt2)
        contrib(:) = prefac * derivEwaldReciprocal(rr, recPoint, alpha, volume)
        deriv(:,iAt1) = deriv(:,iAt1) + contrib
        deriv(:,iAt2) = deriv(:,iAt2) - contrib
      end do
    end do
    !$OMP  END PARALLEL DO

  end subroutine addInvRPrimeXlbomd_periodic


#:if WITH_SCALAPACK

    !> Calculates the -1/R**2 deriv contribution for extended lagrangian dynamics forces
  subroutine getDInvRXlbomdPeriodicBlacs(grid, descAtomSqr, localShape, coord, nNeighborEwald,&
      & iNeighbor, img2CentCell, recPoint, alpha, volume, dQInAtom, dQOutAtom, deriv)

    !> BLACS grid of the distributed derivative vector
    type(blacsgrid), intent(in) :: grid

    !> Descriptor for an nAtom x nAtom matrix distributed on the grid
    integer, intent(in) :: descAtomSqr(DLEN_)

    !> Local shape of the distributed nAtom x nAtom matrix
    integer, intent(in) :: localShape(:)

    !> coordinates of atoms
    real(dp), intent(in) :: coord(:,:)

    !> Nr. of neighbors for each atom for real part of Ewald
    integer, intent(in) :: nNeighborEwald(:)

    !> List of neighbors for the real space part of Ewald.
    integer, intent(in) :: iNeighbor(0:,:)

    !> Image of each atom in the central cell.
    integer, intent(in) :: img2CentCell(:)

    !> Contains the points included in the reciprocal sum. The set should not include the origin or
    !> inversion related points.
    real(dp), intent(in) :: recPoint(:,:)

    !> Ewald parameter
    real(dp), intent(in) :: alpha

    !> cell volume
    real(dp), intent(in) :: volume

    !> input charge fluctuations
    real(dp), intent(in) :: dQInAtom(:)

    !> output charge fluctuations
    real(dp), intent(in) :: dQOutAtom(:)

    !> energy derivative to add contribution to
    real(dp), intent(out) :: deriv(:,:)

    integer :: ii, jj, iAt1, iAt2, iAt2f, iNeigh, iLoc, jLoc
    real(dp) :: rr(3), contrib(3), prefac
    logical :: tLocal

    deriv(:,:) = 0.0_dp

    ! Real space contribution
    do jj = 1, localShape(2)
      iAt1 = scalafx_indxl2g(jj, descAtomSqr(NB_), grid%mycol, descAtomSqr(CSRC_), grid%ncol)
      do iNeigh = 1, nNeighborEwald(iAt1)
        iAt2 = iNeighbor(iNeigh, iAt1)
        iAt2f = img2CentCell(iAt2)
        call scalafx_islocal(grid, descAtomSqr, iAt2f, iAt1, tLocal, iLoc, jLoc)
        if (tLocal .and. iAt2f /= iAt1) then
          rr(:) = coord(:,iAt1) - coord(:,iAt2)
          prefac = dQOutAtom(iAt1) * dQInAtom(iAt2f) + dQInAtom(iAt1) * dQOutAtom(iAt2f) &
              & - dQInAtom(iAt1) * dQInAtom(iAt2f)
          contrib(:) = prefac * derivRTerm(rr, alpha)
          deriv(:,iAt1) = deriv(:,iAt1) + contrib
          ! Neighbor list only contains lower triange: add also to screw symmetric equivalent
          deriv(:,iAt2f) = deriv(:,iAt2f) - contrib
        end if
      end do
    end do

    ! Reciprocal space contribution
    do jj = 1, localShape(2)
      iAt1 = scalafx_indxl2g(jj, descAtomSqr(NB_), grid%mycol, descAtomSqr(CSRC_), grid%ncol)
      do ii= 1, localShape(1)
        iAt2 = scalafx_indxl2g(ii, descAtomSqr(MB_), grid%myrow, descAtomSqr(RSRC_), grid%nrow)
        if (iAt2 /= iAt1) then
          rr(:) = coord(:,iAt1) - coord(:,iAt2)
          prefac = dQOutAtom(iAt1) * dQInAtom(iAt2) + dQInAtom(iAt1) * dQOutAtom(iAt2) &
              & - dQInAtom(iAt1) * dQInAtom(iAt2)
          contrib(:) = prefac * derivEwaldReciprocal(rr, recPoint, alpha, volume)
          deriv(:,iAt1) = deriv(:,iAt1) + contrib
        end if
      end do
    end do

  end subroutine getDInvRXlbomdPeriodicBlacs

#:endif


  !> Calculates the -1/R**2 deriv contribution for charged atoms interacting with a group of charged
  !> objects (like point charges) for the periodic case, without storing anything.
  subroutine addInvRPrime_periodic_asymm(deriv0, deriv1, nAtom0, nAtom1, coord0, coord1, charge0,&
      & charge1, rVec, gVec, alpha, vol)

    !> Contains the derivative for the first group on exit
    real(dp), intent(inout) :: deriv0(:,:)

    !> Contains the derivative for the second group on exit
    real(dp), intent(inout) :: deriv1(:,:)

    !> Number of atoms in the first group
    integer, intent(in) :: nAtom0

    !> Number of atoms in the second group
    integer, intent(in) :: nAtom1

    !> List of atomic coordinates (first group)
    real(dp), intent(in) :: coord0(:,:)

    !> List of the point charge coordinates (second group)
    real(dp), intent(in) :: coord1(:,:)

    !> Charge of the atoms in group 1.
    real(dp), intent(in) :: charge0(:)

    !> Charge of the point charges.
    real(dp), intent(in) :: charge1(:)

    !> Lattice vectors to be used for the real Ewald summation
    real(dp), intent(in) :: rVec(:,:)

    !> Lattice vectors to be used for the reciprocal Ewald summation.
    real(dp), intent(in) :: gVec(:,:)

    !> Parameter of the Ewald summation
    real(dp), intent(in) :: alpha

    !> Volume of the supercell.
    real(dp), intent(in) :: vol

    integer :: iAt0, iAt1
    real(dp) :: dist, vect(3), fTmp(3)

    @:ASSERT(size(deriv0, dim=1) == 3)
    @:ASSERT(size(deriv0, dim=2) == nAtom0)
    @:ASSERT(size(deriv1, dim=1) == 3)
    @:ASSERT(size(deriv1, dim=2) >= nAtom1)
    @:ASSERT(size(coord0, dim=1) == 3)
    @:ASSERT(size(coord0, dim=2) >= nAtom0)
    @:ASSERT(size(coord1, dim=1) == 3)
    @:ASSERT(size(coord1, dim=2) >= nAtom1)
    @:ASSERT(size(charge0) == nAtom0)
    @:ASSERT(size(charge1) == nAtom1)
    @:ASSERT(size(rVec, dim=1) == 3)
    @:ASSERT(size(gVec, dim=1) == 3)
    @:ASSERT(vol > 0.0_dp)

    !$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(iAt0,iAt1,vect,dist,fTmp) &
    !$OMP& SCHEDULE(RUNTIME) REDUCTION(+:deriv0,deriv1)
    do iAt0 = 1, nAtom0
      do iAt1 = 1, nAtom1
        vect(:) = coord0(:,iAt0) - coord1(:,iAt1)
        dist = sqrt(sum(vect(:)**2))
        fTmp(:) = (derivEwaldReal(vect, rVec, alpha)&
            & + derivEwaldReciprocal(vect, gVec, alpha, vol)) * charge0(iAt0) * charge1(iAt1)
        deriv0(:,iAt0) = deriv0(:,iAt0) + fTmp(:)
        deriv1(:,iAt1) = deriv1(:,iAt1) - fTmp(:)
      end do
    end do
    !$OMP  END PARALLEL DO

  end subroutine addInvRPrime_periodic_asymm


  !> Get optimal alpha-parameter for the Ewald summation by finding alpha, where decline of real and
  !> reciprocal part of Ewald are equal.
  !> The function stops, if the optimal alpha cannot be found.
  function getOptimalAlphaEwald(latVec, recVec, volume, tolerance) result(alpha)

    !> Lattice vectors.
    real(dp), intent(in) :: latVec(:,:)

    !> Reciprocal vectors.
    real(dp), intent(in) :: recVec(:,:)

    !> Volume of the unit cell.
    real(dp), intent(in) :: volume

    !> Tolerance for difference in real and rec. part.
    real(dp), intent(in) :: tolerance

    !> Optimal alpha.
    real(dp) :: alpha

    real(dp) :: alphaLeft, alphaRight
    real(dp), parameter :: alphaInit = 1.0e-8_dp

    real(dp) :: minG, minR, diff
    integer :: iIter
    integer :: iError
    character(len=100) :: errorString

    @:ASSERT(all(shape(latVec) == (/3, 3/)))
    @:ASSERT(all(shape(recVec) == (/3, 3/)))
    @:ASSERT(volume > 0.0_dp)
    @:ASSERT(tolerance > 0.0_dp)

    minG = sqrt(minval(sum(recVec(:,:)**2, dim=1)))
    minR = sqrt(minval(sum(latVec(:,:)**2, dim=1)))

    iError = 0
    alpha = alphaInit
    diff = diffRecReal(alpha, minG, minR, volume)
    do while (diff < -tolerance .and. alpha <= huge(1.0_dp))
      alpha = 2.0_dp * alpha
      diff = diffRecReal(alpha, minG, minR, volume)
    end do
    if (alpha > huge(1.0_dp)) then
      iError = 1
    elseif (alpha == alphaInit) then
      iError = 2
    end if

    if (iError == 0) then
      alphaLeft = 0.5_dp * alpha
      do while (diff < tolerance .and. alpha <= huge(1.0_dp))
        alpha = 2.0_dp * alpha
        diff = diffRecReal(alpha, minG, minR, volume)
      end do
      if (alpha > huge(1.0_dp)) then
        iError = 3
      end if
    end if

    if (iError == 0) then
      alphaRight = alpha
      alpha = (alphaLeft + alphaRight) / 2.0
      iIter = 0
      diff = diffRecReal(alpha, minG, minR, volume)
      do while (abs(diff) > tolerance .and. iIter <= nSearchIter)
        if (diff < 0) then
          alphaLeft = alpha
        else
          alphaRight = alpha
        end if
        alpha = (alphaLeft + alphaRight) / 2.0
        diff = diffRecReal(alpha, minG, minR, volume)
        iIter = iIter + 1
      end do
      if (iIter > nSearchIter) then
        iError = 4
      end if
    end if

    if (iError /= 0) then
      !alpha = exp(-0.310104 * log(volume) + 0.786382) / 2.0
99000 format ('Failure in determining optimal alpha for Ewaldsum.', ' Error code: ',I3)
      write(errorString, 99000) iError
      call error(errorString)
    end if

  end function getOptimalAlphaEwald


  !> Returns the longest reciprocal vector which gives a bigger contribution to the Ewald sum than a
  !> certain tolerance.
  function getMaxGEwald(alpha, volume, minValue) result(xx)

    !> Parameter of the ewald summation.
    real(dp), intent(in) :: alpha

    !> Volume of the unit cell.
    real(dp), intent(in) :: volume

    !> Tolerance value.
    real(dp), intent(in) :: minValue

    !> magnitude of reciprocal vector
    real(dp) :: xx

    real(dp), parameter :: gInit = 1.0e-8_dp
    real(dp) :: xLeft, xRight, yLeft, yRight, yy
    integer :: iError, iIter
    character(len=100) :: errorString

    iError = 0
    xx = gInit
    yy = gTerm(xx, alpha, volume)
    do while (yy > minValue .and. xx <= huge(1.0_dp))
      xx = 2.0_dp * xx
      yy = gTerm(xx, alpha, volume)
    end do
    if (xx > huge(1.0_dp)) then
      iError = 1
    elseif (xx == gInit) then
      iError = 2
    end if

    if (iError == 0) then
      xLeft = 0.5_dp * xx
      yLeft = gTerm(xLeft, alpha, volume)
      xRight = xx
      yRight = yy

      iIter = 1
      do while (yLeft - yRight > minValue .and. iIter <= nSearchIter)
        xx = 0.5_dp * (xLeft + xRight)
        yy = gTerm(xx, alpha, volume)
        if (yy >= minValue) then
          xLeft = xx
          yLeft = yy
        else
          xRight = xx
          yRight = yy
        end if
        iIter = iIter + 1
      end do
      if (iIter > nSearchIter) then
        iError = 3
      end if
    end if

    if (iError /= 0) then
99010 format ('Failure in getMaxGEwald.', ' Error nr: ',I3)
      write(errorString, 99010) iError
      call error(errorString)
    end if

  end function getMaxGEwald


  !> Returns the longest real space vector which gives a bigger contribution to the Ewald sum than a
  !> certain tolerance.
  function getMaxREwald(alpha, minValue) result(xx)

    !> Parameter of the ewald summation.
    real(dp), intent(in) :: alpha

    !> Tolerance value.
    real(dp), intent(in) :: minValue

    !> Magnitude of real space vector
    real(dp) :: xx

    real(dp), parameter :: rInit = 1.0e-8_dp
    real(dp) :: xLeft, xRight, yLeft, yRight, yy
    integer :: iError, iIter
    character(len=100) :: errorString

    iError = 0
    xx = rInit
    yy = rTerm(xx, alpha)
    do while (yy > minValue .and. xx <= huge(1.0_dp))
      xx = 2.0_dp * xx
      yy = rTerm(xx, alpha)
    end do
    if (xx > huge(1.0_dp)) then
      iError = 1
    elseif (xx == rInit) then
      iError = 2
    end if

    if (iError == 0) then
      xLeft = 0.5_dp * xx
      yLeft = rTerm(xLeft, alpha)
      xRight = xx
      yRight = yy

      iIter = 1
      do while (yLeft - yRight > minValue .and. iIter <= nSearchIter)
        xx = 0.5_dp * (xLeft + xRight)
        yy = rTerm(xx, alpha)
        if (yy >= minValue) then
          xLeft = xx
          yLeft = yy
        else
          xRight = xx
          yRight = yy
        end if
        iIter = iIter + 1
      end do
      if (iIter > nSearchIter) then
        iError = 3
      end if
    end if

    if (iError /= 0) then
99020 format ('Failure in getMaxREwald.', ' Error nr: ',I3)
      write(errorString, 99020) iError
      call error(errorString)
    end if

  end function getMaxREwald


  !> Returns the Ewald sum for a given lattice in a given point.
  function ewald(rr, rVec, gVec, alpha, vol)

    !> Vector where to calculate the Ewald sum.
    real(dp), intent(in) :: rr(:)

    !> Real space vectors to sum over. (Should contain origin).
    real(dp), intent(in) :: rVec(:,:)

    !> Reciprocal space vectors to sum over (Should not contain either origin nor inversion related
    !> points).
    real(dp), intent(in) :: gVec(:,:)

    !> Parameter for the Ewald summation.
    real(dp), intent(in) :: alpha

    !> Volume of the real space unit cell.
    real(dp), intent(in) :: vol

    !> Result
    real(dp) :: ewald

    ewald = ewaldReciprocal(rr, gVec, alpha, vol) + ewaldReal(rr, rVec, alpha) - pi / (vol*alpha**2)
    if (sum(rr(:)**2) < tolSameDist2) then
      ewald = ewald - 2.0_dp * alpha / sqrt(pi)
    end if

  end function ewald


  !> Returns the reciprocal part of the Ewald sum.
  function ewaldReciprocal(rr, gVec, alpha, vol) result(recSum)

    !> Vector where to calculate the Ewald sum.
    real(dp), intent(in) :: rr(:)

    !> Reciprocal space vectors to sum over (Should not contain either origin nor inversion related
    !> points).
    real(dp), intent(in) :: gVec(:,:)

    !> Parameter for the Ewald summation.
    real(dp), intent(in) :: alpha

    !> Volume of the real space unit cell.
    real(dp), intent(in) :: vol

    !> contribution to the sum
    real(dp) :: recSum

    real(dp) :: gg(3), g2
    integer :: iG

    @:ASSERT(size(gVec, dim=1) == 3)
    @:ASSERT(size(rr) == 3)
    @:ASSERT(vol > 0.0_dp)

    recSum = 0.0_dp
    do iG = 1, size(gVec, dim=2)
      gg = gVec(:,iG)
      g2 = sum(gg(:)**2)
      recSum = recSum + exp(-g2/(4.0_dp*alpha**2))/g2 * cos(dot_product(gg,rr))
    end do
    ! note factor of 2 as only summing half of reciprocal space
    recSum = 2.0_dp * recSum * 4.0_dp * pi / vol

  end function ewaldReciprocal


  !> Returns the derivative of the reciprocal part of the Ewald sum.
  function derivEwaldReciprocal(rr, gVec, alpha, vol) result(recSum)

    !> Vector where to calculate the Ewald sum.
    real(dp), intent(in) :: rr(:)

    !> Reciprocal space vectors to sum over (Should not contain either origin nor inversion related
    !> points).
    real(dp), intent(in) :: gVec(:,:)

    !> Parameter for the Ewald summation.
    real(dp), intent(in) :: alpha

    !> Volume of the real space unit cell.
    real(dp), intent(in) :: vol

    !> contribution to the derivative value
    real(dp) :: recSum(3)

    real(dp) :: gg(3), g2
    integer :: iG

    @:ASSERT(size(gVec, dim=1) == 3)
    @:ASSERT(size(rr) == 3)
    @:ASSERT(vol > 0.0_dp)

    recSum(:) = 0.0_dp
    do iG = 1, size(gVec, dim=2)
      gg(:) = gVec(:,iG)
      g2 = sum(gg(:)**2)
      recSum(:) = recSum(:) - gg(:)*sin(dot_product(gg,rr))*exp(-g2/(4.0_dp*alpha**2))/g2
    end do
    ! note factor of 2 as only summing over half of reciprocal space
    recSum(:) = 2.0_dp * recSum(:) * 4.0_dp * pi / vol

  end function derivEwaldReciprocal


  !> Returns the real space part of the Ewald sum.
  function ewaldReal(rr, rVec, alpha) result(realSum)

    !> Real space vectors to sum over. (Should contain origin).
    real(dp), intent(in) :: rVec(:,:)

    !> Parameter for the Ewald summation.
    real(dp), intent(in) :: alpha

    !> Vector where to calculate the Ewald sum.
    real(dp), intent(in) :: rr(:)

    !> contribution to sum
    real(dp) :: realSum

    real(dp) :: absRR
    integer :: iR

    @:ASSERT(size(rVec, dim=1) == 3)
    @:ASSERT(size(rr) == 3)

    realSum = 0.0_dp;
    do iR = 1, size(rVec, dim=2)
      absRR = sqrt(sum((rr(:) + rVec(:,iR))**2))
      if (absRR < tolSameDist) then
        cycle
      end if
      realSum = realSum + erfcwrap(alpha*absRR)/absRR
    end do

  end function ewaldReal


  !> Returns the derivative of the real space part of the Ewald sum.
  function derivEwaldReal(rdiff, rVec, alpha) result(dewr)

    !> Vector where to calculate the Ewald sum.
    real(dp), intent(in) :: rdiff(:)

    !> Real space vectors to sum over. (Should contain origin).
    real(dp), intent(in) :: rVec(:,:)

    !> Parameter for the Ewald summation.
    real(dp), intent(in) :: alpha

    !> contribution to derivative
    real(dp) :: dewr(3)

    real(dp) :: rNew(3)
    real(dp) :: rr
    integer :: iR

    @:ASSERT(size(rVec, dim=1) == 3)
    @:ASSERT(size(rdiff) == 3)

    dewr = 0.0_dp;
    do iR = 1, size(rVec, dim=2)
      rNew(:) = rdiff(:) + rVec(:,iR)
      rr = sqrt(sum(rNew**2))
      if (rr < tolSameDist2) then
        cycle
      end if
      dewr(:) = dewr + rNew(:) * (-2.0_dp/sqrt(pi)*exp(-alpha*alpha*rr*rr)* alpha*rr&
          & - erfcwrap(alpha*rr))/(rr*rr*rr)
    end do

  end function derivEwaldReal


  !> Returns the difference in the decrease of the real and reciprocal parts of the Ewald sum.
  !> In order to make the real space part shorter as the reciprocal space part, the inclinations are
  !> taken at different points for the the real space and the reciprocal space part.
  function diffRecReal(alpha, minG, minR, volume) result(diff)

    !> Parameter for the Ewald summation.
    real(dp), intent(in) :: alpha

    !> Length of the shortest reciprocal space vector in the sum.
    real(dp), intent(in) :: minG

    !> Length of the shortest real space vector in the sum.
    real(dp), intent(in) :: minR

    !> Volume of the real space unit cell.
    real(dp), intent(in) :: volume

    !> difference between changes in the two terms
    real(dp) :: diff

    @:ASSERT(volume > 0.0_dp)

    diff = ((gTerm(4.0_dp*minG, alpha, volume) &
        &- gTerm(5.0_dp*minG, alpha, volume))) &
        &- (rTerm(2.0_dp*minR, alpha) - rTerm(3.0_dp*minR, alpha))

  end function diffRecReal


  !> Returns the max. value of a term in the reciprocal space part of the Ewald summation for a
  !> given vector length.
  function gTerm(gg, alpha, vol)

    !> Length of the reciprocal space vector.
    real(dp), intent(in) :: gg

    !> Parameter of the Ewald summation.
    real(dp), intent(in) :: alpha

    !> Volume of the real space unit cell.
    real(dp), intent(in) :: vol

    !> reciprocal term
    real(dp) :: gTerm

    gTerm = 4.0_dp*pi*(exp(-0.25_dp*gg**2/(alpha**2))/(vol*gg**2))

  end function gTerm


  !> Returns the max. value of a term in the real space part of the Ewald summation for a given
  !> vector length.
  function rTerm(rr, alpha)

    !> Length of the real space vector.
    real(dp), intent(in) :: rr

    !> Parameter of the Ewald summation.
    real(dp), intent(in) :: alpha

    !> real space term
    real(dp) :: rTerm

    @:ASSERT(rr >= epsilon(1.0_dp))

    rTerm = erfcwrap(alpha*rr)/rr

  end function rTerm


  !> Returns the derivative of a term in the real space part of the Ewald summation for a given
  !> vector length.
  function derivRTerm(r, alpha)

    !> Length of the real space vector.
    real(dp), intent(in) :: r(3)

    !> Parameter of the Ewald summation.
    real(dp), intent(in) :: alpha

    !> real space derivative term
    real(dp) :: derivRTerm(3)

    real(dp) :: rr
    rr = sqrt(sum(r(:)**2))

    @:ASSERT(rr >= epsilon(1.0_dp))

    derivRTerm (:) = r(:)*(-2.0_dp/sqrt(pi)*exp(-alpha*alpha*rr*rr)* &
        & alpha*rr - erfcwrap(alpha*rr))/(rr*rr*rr)

  end function derivRTerm


  !> Calculates the stress tensor derivatives of the Ewald electrostatics
  !> Aguard and Madden J Chem Phys 119 7471 (2003)
  subroutine invR_stress(stress, nAtom, coord, nNeighborEwald, iNeighbor, &
      & img2CentCell, recPoint, alpha, volume, q)

    !> Stress tensor
    real(dp), intent(out) :: stress(:,:)

    !> Number of atoms.
    integer, intent(in) :: nAtom

    !> List of atomic coordinates (all atoms).
    real(dp), intent(in) :: coord(:,:)

    !> Nr. of neighbors for each atom for real part of Ewald.
    integer, intent(in) :: nNeighborEwald(:)

    !> List of neighbors for the real space part of Ewald.
    integer, intent(in) :: iNeighbor(0:,:)

    !> Image of each atom in the central cell.
    integer, intent(in) :: img2CentCell(:)

    !> Contains the points included in the reciprocal sum. The set should not include the origin or
    !> inversion related points.
    real(dp), intent(in) :: recPoint(:,:)

    !> Parameter for Ewald summation.
    real(dp), intent(in) :: alpha

    !> Volume of the real space unit cell.
    real(dp), intent(in) :: volume

    !> charges in the cell
    real(dp), intent(in) :: q(:)

    integer :: iAtom1, iAtom2, iAtom2f, iNeigh, iInv, ii, jj, kk
    real(dp) :: r(3), f(3), g(3), g2, intermed, intermed2
    real(dp) :: stressTmp(3,3)

    @:ASSERT(all(shape(stress)==(/3,3/)))
    @:ASSERT(size(coord, dim=1) == 3)
    @:ASSERT(size(coord, dim=2) >= nAtom)
    @:ASSERT(size(nNeighborEwald) == nAtom)
    @:ASSERT(size(iNeighbor, dim=2) == nAtom)
    @:ASSERT(size(q) == nAtom)
    @:ASSERT(volume > 0.0_dp)

    stress(:,:) = 0.0_dp

    ! Reciprocal space part of the Ewald sum.
    do ii = 1, size(recpoint, dim=2)
      do iInv = -1, 1, 2
        g(:) = real(iInv,dp)*recpoint(:,ii)
        intermed = 0.0_dp
        intermed2 = 0.0_dp
        !$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(iAtom1) &
        !$OMP& SCHEDULE(RUNTIME) REDUCTION(+:intermed,intermed2)
        do iAtom1 = 1, nAtom
          intermed = intermed &
              & + q(iAtom1)*cos(dot_product(g(:),coord(:,iAtom1)))
          intermed2 = intermed2 &
              & + q(iAtom1)*sin(dot_product(g(:),coord(:,iAtom1)))
        end do
        !$OMP  END PARALLEL DO
        intermed = intermed**2 + intermed2**2
        g2 = sum(g(:)**2)
        intermed = intermed*exp(-g2/(4.0_dp*alpha*alpha))/g2
        stressTmp(:,:) = 0.0_dp
        do jj = 1, 3
          stressTmp(jj,jj) = 1.0_dp
          do kk = 1,3
            stressTmp(kk,jj) = stressTmp(kk,jj) &
                & -2.0_dp*(1.0_dp/(4.0_dp*alpha*alpha) + 1.0_dp/g2) &
                & *g(kk)*g(jj)
          end do
        end do
        stress = stress + stressTmp * intermed
      end do
    end do
    stress = -stress * 4.0_dp * pi / volume

    ! Real space part of the Ewald sum.
    !$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(iAtom1,iNeigh,iAtom2,iAtom2f,r,f,ii,jj) &
    !$OMP& SCHEDULE(RUNTIME) REDUCTION(+:stress)
    do iAtom1 = 1, nAtom
      do iNeigh = 1, nNeighborEwald(iAtom1)
        iAtom2 = iNeighbor(iNeigh, iAtom1)
        iAtom2f = img2CentCell(iAtom2)
        r(:) = coord(:,iAtom1)-coord(:,iAtom2)
        f(:) = derivRTerm(r,alpha) * Q(iAtom1) * Q(iAtom2f)
        if (iAtom2f /= iAtom1) then
          do ii = 1, 3
            do jj = 1, 3
              stress(jj,ii) = stress(jj,ii) + (r(jj)*f(ii) + f(jj)*r(ii))
            end do
          end do
        else
          do ii = 1, 3
            do jj = 1, 3
              stress(jj,ii) = stress(jj,ii) + 0.5_dp*(r(jj)*f(ii) + f(jj)*r(ii))
            end do
          end do
        end if
      end do
    end do
    !$OMP  END PARALLEL DO

    stress = stress / volume

  end subroutine invR_stress


end module coulomb
