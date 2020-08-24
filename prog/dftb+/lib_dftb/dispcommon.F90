!--------------------------------------------------------------------------------------------------!
!  DFTB+: general package for performing fast atomistic simulations                                !
!  Copyright (C) 2006 - 2020  DFTB+ developers group                                               !
!                                                                                                  !
!  See the LICENSE file for terms of usage and distribution.                                       !
!--------------------------------------------------------------------------------------------------!

#:include 'common.fypp'

!> Common routines for the dispersion modules.
!> Periodic summation from the following references:
!> N. Karasawa et al., J. Phys. Chem. 93, 7320-7327 (1989)
!> Zhou-Min Chen et al., J. Comp. Chem. 18, 1365 (1997)
module dftbp_dispcommon
  use dftbp_assert
  use dftbp_accuracy
  use dftbp_constants, only : pi
  use dftbp_environment, only : TEnvironment
  use dftbp_errorfunction
  use dftbp_message
  use dftbp_simplealgebra, only : cross3
  use dftbp_schedule, only : distributeRangeInChunks, assembleChunks
  use dftbp_sorting
  implicit none
  private

  public :: addDispEGr_per_species, addDispEGr_per_atom
  public :: getOptimalEta, getMaxRDispersion, getMaxGDispersion

contains


  !> Adds the energy per atom and the gradients for periodic 1/r^6 summation.
  !> Fast converging Ewald like summation on 1/r^6 type interactions.  For details see the
  !> references in the module description.
  !> Note: the interaction parameter C6 is specified atomwise.
  subroutine addDispEGr_per_atom(env, nAtom, coords, nNeighbourSK, iNeighbour, neighDist2,&
      & img2CentCell, c6, eta, vol, gLatVecs, energies, gradients, stress)

    !> Computational environment settings
    type(TEnvironment), intent(in) :: env

    !> Nr. of atoms (without periodic images)
    integer, intent(in) :: nAtom

    !> Coordinates of the atoms (including images)
    real(dp), intent(in) :: coords(:,:)

    !> Nr. of neighbours for each atom
    integer, intent(in) :: nNeighbourSK(:)

    !> Neighbourlist.
    integer, intent(in) :: iNeighbour(0:,:)

    !> Square distances of the neighbours.
    real(dp), intent(in) :: neighDist2(0:,:)

    !> Mapping of periodic images onto the central cell
    integer, intent(in) :: img2CentCell(:)

    !> Van der Waals coefficients (nAtom, nAtom)
    real(dp), intent(in) :: c6(:,:)

    !> Controling partitioning between real and reciprocal space sum
    real(dp), intent(in) :: eta

    !> Volume of the unit cell
    real(dp), intent(in) :: vol

    !> Set of reciprocal space vectors to include in the sum
    real(dp), intent(in) :: gLatVecs(:,:)

    !> Updated energy vector at return
    real(dp), intent(inout) :: energies(:)

    !> Updated gradient vector at return
    real(dp), intent(inout) :: gradients(:,:)

    !> Updated stress tensor
    real(dp), intent(inout) :: stress(:,:)

    integer :: iAtFirst, iAtLast
    integer :: iAt1, iNeigh, iAt2, iAt2f, iG, ii
    real(dp) :: rSum, rSum3(3), gSum, gSum3(3), gg(3), ggAbs, vec(3)
    real(dp) :: aam2, bb, bbm2, rTmp, rTmp2, rTmp3, rc, r3c, gc, g3c, ddp
    real(dp) :: etam3, rTmp33, gsum33(3,3)
    real(dp), allocatable :: localDeriv(:,:), localStress(:, :), localEnergies(:)

    @:ASSERT(size(energies) == nAtom)
    @:ASSERT(all(shape(gradients) == (/ 3, nAtom /)))
    @:ASSERT(all(shape(stress) == (/ 3, 3 /)))

    etam3 = eta**(-3)
    rc = 0.5_dp * etam3 * etam3
    r3c = 2.0_dp * rc / (eta * eta)
    gc = pi**1.5_dp / (12.0_dp * vol)
    g3c = pi**1.5_dp / (6.0_dp * vol)

    call distributeRangeInChunks(env, 1, nAtom, iAtFirst, iAtLast)

    allocate(localEnergies(nAtom), localDeriv(3, nAtom), localStress(3, 3))
    localEnergies(:) = 0.0_dp
    localDeriv(:,:) = 0.0_dp
    localStress(:,:) = 0.0_dp

    !$omp parallel do default(none) schedule(runtime) &
    !$omp reduction(+:localEnergies, localDeriv, localStress) shared(etam3, nAtom) &
    !$omp shared(iAtFirst, iAtLast, nNeighbourSK, iNeighbour, coords, c6, eta) &
    !$omp shared(img2CentCell, neighDist2, vol, gLatVecs, rc, r3c, gc, g3c) &
    !$omp private(iAt1, iNeigh, iAt2, iAt2f, iG, ii, rSum, rSum3, gSum, gSum3) &
    !$omp private(gg, ggAbs, vec, aam2, bb, bbm2, rTmp, rTmp2, rTmp3, ddp) &
    !$omp private(rTmp33, gsum33)
    do iAt1 = iAtFirst, iAtLast
      do iNeigh = 1, nNeighbourSK(iAt1)
        iAt2 = iNeighbour(iNeigh, iAt1)
        vec = coords(:,iAt1) - coords(:,iAt2)
        iAt2f = img2CentCell(iAt2)
        aam2 = (sqrt(neighDist2(iNeigh, iAt1))/eta)**(-2)
        rSum =  rc * c6(iAt2f, iAt1) * ((aam2 + 1.0_dp)*aam2 + 0.5_dp)*aam2 * exp(-1.0_dp / aam2)
        rSum3(:) = r3c * c6(iAt2f, iAt1) * exp(-1.0_dp / aam2)&
            & * (((6.0_dp*aam2 + 6.0_dp)*aam2 + 3.0_dp)*aam2 + 1.0_dp)*aam2 * vec
        localEnergies(iAt1) = localEnergies(iAt1) - rSum
        localDeriv(:,iAt1) = localDeriv(:,iAt1) + rSum3
        if (iAt1 /= iAt2f) then
          localEnergies(iAt2f) = localEnergies(iAt2f) - rSum
          localDeriv(:,iAt2f) = localDeriv(:,iAt2f) - rSum3
          do ii = 1, 3
            localStress(:,ii) = localStress(:,ii) - rSum3 * vec(ii) / vol
          end do
        else
          do ii = 1, 3
            localStress(:,ii) = localStress(:,ii) - 0.5_dp * rSum3 * vec(ii) / vol
          end do
        end if
      end do

      gSum = 0.0_dp
      gSum3(:) = 0.0_dp
      gSum33(:,:) = 0.0_dp
      do iG = 1, size(gLatVecs, dim=2)
        gg = gLatVecs(:, iG)
        ggAbs = sqrt(sum(gg**2))
        bb = 0.5_dp * ggAbs * eta
        bbm2 = bb**(-2)
        rTmp = 0.0_dp
        rTmp2 = 0.0_dp
        do iAt2 = 1, nAtom
          ddp = dot_product(gg, coords(:,iAt1) - coords(:,iAt2))
          rTmp = rTmp + cos(ddp) * c6(iAt1, iAt2)
          rTmp2 = rTmp2 + sin(ddp) * c6(iAt1, iAt2)
        end do
        rTmp3 = sqrt(pi) * erfcwrap(bb) + (0.5_dp * bbm2 - 1.0_dp) / bb * exp(-1.0_dp / bbm2)
        rTmp33 = sqrt(pi) * erfcwrap(bb) - exp(-bb*bb)/ bb
        gSum = gSum + rTmp * ggAbs**3 * rTmp3
        gSum3(:) = gSum3 + gg * rTmp2 * ggAbs**3 * rTmp3
        do ii = 1, 3
          gSum33(:,ii) = gSum33(:,ii) + rTmp * 3.0_dp * ggAbs * rTmp33 * gg * gg(ii)
        end do
      end do
      gSum = gSum * gc
      gSum3(:) = gSum3 * g3c
      gSum33 = gSum33 * gc
      localEnergies(iAt1) = localEnergies(iAt1) - gSum + (c6(iAt1,iAt1)/12.0_dp * etam3&
          & - pi**1.5 * sum(c6(:,iAt1))/(6.0_dp * vol)) * etam3
      do ii = 1, 3
        localStress(ii,ii) = localStress(ii,ii) - gSum/vol&
            & - (pi**1.5 * sum(c6(:,iAt1))/(6.0_dp * vol*vol)) * etam3
      end do
      localStress = localStress - gSum33 / vol
      localDeriv(:,iAt1) =  localDeriv(:,iAt1) +  gSum3
    end do
    !$omp end parallel do

    call assembleChunks(env, localEnergies)
    call assembleChunks(env, localDeriv)
    call assembleChunks(env, localStress)

    energies(:) = energies + localEnergies
    gradients(:,:) = gradients + localDeriv
    stress(:,:) = stress + localStress

  end subroutine addDispEGr_per_atom


  !> Adds the energy per atom and the gradients for periodic 1/r^6 summation.
  !> Fast converging Ewald like summation on 1/r^6 type interactions.  The 1/r^12 term is
  !> summed in direct space.
  !> Note: the interaction coefficients (c6) are specified specieswise.
  subroutine addDispEGr_per_species(env, nAtom, coords, species0, nNeighbourSK, iNeighbour,&
      & neighDist2, img2CentCell, c6, eta, vol, gLatVecs, energies, gradients, stress)

    !> Computational environment settings
    type(TEnvironment), intent(in) :: env

    !> Nr. of atoms (without periodic images)
    integer, intent(in) :: nAtom

    !> Coordinates of the atoms (including images)
    real(dp), intent(in) :: coords(:,:)

    !> chemical species of atoms
    integer, intent(in) :: species0(:)

    !> Nr. of neighbours for each atom
    integer, intent(in) :: nNeighbourSK(:)

    !> Neighbourlist.
    integer, intent(in) :: iNeighbour(0:,:)

    !> Square distances of the neighbours.
    real(dp), intent(in) :: neighDist2(0:,:)

    !> Mapping of periodic images onto the central cell
    integer, intent(in) :: img2CentCell(:)

    !> Van der Waals coefficients (nSpecies, nSpecies)
    real(dp), intent(in) :: c6(:,:)

    !> Controling partitioning between real and reciprocal space sum
    real(dp), intent(in) :: eta

    !> Volume of the unit cell
    real(dp), intent(in) :: vol

    !> Set of reciprocal space vectors to include in the sum
    real(dp), intent(in) :: gLatVecs(:,:)

    !> Updated energy vector at return
    real(dp), intent(inout) :: energies(:)

    !> Updated gradient vector at return
    real(dp), intent(inout) :: gradients(:,:)

    !> Updated stress tensor
    real(dp), intent(inout) :: stress(:,:)

    integer :: iAtFirst, iAtLast
    integer :: iAt1, iNeigh, iAt2, iAt2f, iG, iSp1, iSp2, ii
    real(dp) :: rSum, rSum3(3), gSum, gSum3(3),gsum33(3,3),gg(3), ggAbs, vec(3)
    real(dp) :: aam2, bb, bbm2, rTmp, rTmp2, rTmp3, rc, r3c, gc, g3c, ddp, etam3
    real(dp) :: rTmp33
    real(dp), allocatable :: localDeriv(:,:), localStress(:, :), localEnergies(:)

    @:ASSERT(size(energies) == nAtom)
    @:ASSERT(all(shape(gradients) == (/ 3, nAtom /)))
    @:ASSERT(all(shape(stress) == (/ 3, 3 /)))

    etam3 = eta**(-3)
    rc = 0.5_dp * etam3 * etam3
    r3c = 2.0_dp * rc / (eta * eta)
    gc = pi**1.5_dp / (12.0_dp * vol)
    g3c = pi**1.5_dp / (6.0_dp * vol)

    call distributeRangeInChunks(env, 1, nAtom, iAtFirst, iAtLast)

    allocate(localEnergies(nAtom), localDeriv(3, nAtom), localStress(3, 3))
    localEnergies(:) = 0.0_dp
    localDeriv(:,:) = 0.0_dp
    localStress(:,:) = 0.0_dp

    !$omp parallel do default(none) schedule(runtime) &
    !$omp reduction(+:localEnergies, localDeriv, localStress) &
    !$omp shared(iAtFirst, iAtLast, nNeighbourSK, iNeighbour, coords, c6, eta) &
    !$omp shared(img2CentCell, neighDist2, vol, gLatVecs, rc, r3c, gc, g3c) &
    !$omp shared(species0, etam3, nAtom) &
    !$omp private(iAt1, iNeigh, iAt2, iAt2f, iG, ii, rSum, rSum3, gSum, gSum3) &
    !$omp private(gg, ggAbs, vec, aam2, bb, bbm2, rTmp, rTmp2, rTmp3, ddp) &
    !$omp private(rTmp33, gsum33, iSp1, iSp2)
    do iAt1 = iAtFirst, iAtLast
      iSp1 = species0(iAt1)
      do iNeigh = 1, nNeighbourSK(iAt1)
        iAt2 = iNeighbour(iNeigh, iAt1)
        vec = coords(:,iAt1) - coords(:,iAt2)
        iAt2f = img2CentCell(iAt2)
        iSp2 = species0(iAt2f)
        aam2 = (sqrt(neighDist2(iNeigh, iAt1))/eta)**(-2)
        rSum =  rc * c6(iSp2, iSp1) * ((aam2 + 1.0_dp)*aam2 + 0.5_dp)*aam2 * exp(-1.0_dp / aam2)
        rSum3(:) = r3c * c6(iSp2, iSp1) * exp(-1.0_dp / aam2)&
            & * (((6.0_dp*aam2 + 6.0_dp)*aam2 + 3.0_dp)*aam2 + 1.0_dp)*aam2 * vec
        localEnergies(iAt1) = localEnergies(iAt1) - rSum
        localDeriv(:,iAt1) = localDeriv(:,iAt1) + rSum3
        if (iAt1 /= iAt2f) then
          localEnergies(iAt2f) = localEnergies(iAt2f) - rSum
          localDeriv(:,iAt2f) = localDeriv(:,iAt2f) - rSum3
          do ii = 1, 3
            localStress(:,ii) = localStress(:,ii) - rSum3 * vec(ii) / vol
          end do
        else
          do ii = 1, 3
            localStress(:,ii) = localStress(:,ii) - 0.5_dp * rSum3 * vec(ii) / vol
          end do
        end if
      end do

      gSum = 0.0_dp
      gSum3(:) = 0.0_dp
      gSum33(:,:) = 0.0_dp
      do iG = 1, size(gLatVecs, dim=2)
        gg = gLatVecs(:, iG)
        ggAbs = sqrt(sum(gg**2))
        bb = 0.5_dp * ggAbs * eta
        bbm2 = bb**(-2)
        rTmp = 0.0_dp
        rTmp2 = 0.0_dp
        do iAt2 = 1, nAtom
          iSp2 = species0(iAt2)
          ddp = dot_product(gg, coords(:,iAt1) - coords(:,iAt2))
          rTmp = rTmp + cos(ddp) * c6(iSp1, iSp2)
          rTmp2 = rTmp2 + sin(ddp) * c6(iSp1, iSp2)
        end do
        rTmp3 = sqrt(pi) * erfcwrap(bb) + (0.5_dp * bbm2 - 1.0_dp) / bb * exp(-1.0_dp / bbm2)
        rTmp33 = sqrt(pi) * erfcwrap(bb) - exp(-bb*bb)/ bb
        gSum = gSum + rTmp * ggAbs**3 * rTmp3
        gSum3(:) = gSum3 + gg * rTmp2 * ggAbs**3 * rTmp3
        do ii = 1, 3
          gSum33(:,ii) = gSum33(:,ii) + rTmp * 3.0_dp * ggAbs * rTmp33 * gg * gg(ii)
        end do
      end do
      gSum = gSum * gc
      gSum3(:) = gSum3 * g3c
      gSum33 = gSum33 * gc
      localEnergies(iAt1) = localEnergies(iAt1) - gSum + (c6(iSp1,iSp1)/12.0_dp * etam3&
          & - pi**1.5 * sum(c6(species0(1:nAtom),iSp1))/(6.0_dp * vol)) * etam3
      do ii = 1, 3
        localStress(ii,ii) = localStress(ii,ii) - gSum/vol&
            & - (pi**1.5 * sum(c6(species0(1:nAtom),iSp1))/(6.0_dp * vol*vol)) * etam3
      end do
      localStress = localStress - gSum33 / vol
      localDeriv(:,iAt1) =  localDeriv(:,iAt1) +  gSum3
    end do
    !$omp end parallel do

    call assembleChunks(env, localEnergies)
    call assembleChunks(env, localDeriv)
    call assembleChunks(env, localStress)

    energies(:) = energies + localEnergies
    gradients(:,:) = gradients + localDeriv
    stress(:,:) = stress + localStress

  end subroutine addDispEGr_per_species


  !> Delivers the optimal paramater eta for the Ewald-like summation
  function getOptimalEta(latVecs, vol) result(eta)

    !> Lattice vectors
    real(dp), intent(in) :: latVecs(:,:)

    !> unit cell volume
    real(dp), intent(in) :: vol

    !> Optimal parameter.
    real(dp) :: eta

    real(dp) :: vecLens(3), tmp(3)
    integer :: indx(3)

    vecLens(:) = sqrt(sum(latVecs**2, dim=1))
    call index_heap_sort(indx, vecLens)
    tmp(:) = cross3(latVecs(:,indx(1)), latVecs(:,indx(2)))
    eta = sqrt(vecLens(indx(1)) * vol / (pi * sqrt(sum(tmp**2))))

  end function getOptimalEta


  !> Returns the longest real space vector needed to achive a given accuracy in the Ewald summation
  !> for the dispersion.
  function getMaxRDispersion(eta, c6sum, vol, minValue) result(xx)

    !> Parameter of the ewald summation.
    real(dp), intent(in) :: eta

    !> Sum of the Van der Waals coefficients (for every possible interaction in the system).
    real(dp), intent(in) :: c6sum

    !> Volume of the unit cell
    real(dp), intent(in) :: vol

    !> Tolerance value.
    real(dp), intent(in) :: minValue

    !> Cutoff radius.
    real(dp) :: xx

    real(dp), parameter :: rInit = 1.0e-8_dp
    real(dp) :: xLeft, xRight, yLeft, yRight, yy
    integer :: iError, iIter
    character(lc) :: strError

    iError = 0
    xx = rInit
    yy = getDispRealError(xx, c6sum, vol, eta)
    do while (yy > minValue .and. xx <= huge(1.0_dp))
      xx = 2.0_dp * xx
      yy = getDispRealError(xx, c6sum, vol, eta)
    end do
    if (xx > huge(1.0_dp)) then
      iError = 1
    elseif (xx == rInit) then
      iError = 2
    end if

    if (iError == 0) then
      xLeft = 0.5_dp * xx
      yLeft = getDispRealError(xLeft, c6sum, vol, eta)
      xRight = xx
      yRight = yy

      iIter = 1
      do while (yLeft - yRight > minValue .and. iIter <= nSearchIter)
        xx = 0.5_dp * (xLeft + xRight)
        yy = getDispRealError(xx, c6sum, vol, eta)
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
99020 format ('Failure in getMaxRDispersion_.', ' Error nr: ',I3)
      write(strError, 99020) iError
      call error(strError)
    end if

  end function getMaxRDispersion


  !> Returns the longest reciprocal space vector needed to achive a given accuracy in the Ewald
  !> summation for the dispersion.
  function getMaxGDispersion(eta, c6sum, minValue) result(xx)

    !> Parameter of the ewald summation.
    real(dp), intent(in) :: eta

    !> Sum of the absolute values of the c6 coeffs for every atom pair.
    real(dp), intent(in) :: c6sum

    !> Tolerance value.
    real(dp), intent(in) :: minValue

    !> Cutoff radius.
    real(dp) :: xx

    real(dp), parameter :: gInit = 1.0e-8_dp
    real(dp) :: xLeft, xRight, yLeft, yRight, yy
    integer :: iError, iIter
    character(lc) :: strError

    iError = 0
    xx = gInit
    yy = getDispReciprocalError(xx, c6sum, eta)
    do while (yy > minValue .and. xx <= huge(1.0_dp))
      xx = 2.0_dp * xx
      yy = getDispReciprocalError(xx, c6sum, eta)
    end do
    if (xx > huge(1.0_dp)) then
      iError = 1
    elseif (xx == gInit) then
      iError = 2
    end if

    if (iError == 0) then
      xLeft = 0.5_dp * xx
      yLeft = getDispReciprocalError(xLeft, c6sum, eta)
      xRight = xx
      yRight = yy

      iIter = 1
      do while (yLeft - yRight > minValue .and. iIter <= nSearchIter)
        xx = 0.5_dp * (xLeft + xRight)
        yy = getDispReciprocalError(xx, c6sum, eta)
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
99010 format ('Failure in getMaxGDispersion_.', ' Error nr: ',I3)
      write(strError, 99010) iError
      call error(strError)
    end if

  end function getMaxGDispersion


  !> Returns an estimate of the error of the real space summation for a certain cutoff
  function getDispRealError(rr, c6sum, vol, eta) result(err)

    !> Cutoff radius
    real(dp), intent(in) :: rr

    !> VdW coefficient sum.
    real(dp), intent(in) :: c6sum

    !> Volume of the supercell
    real(dp), intent(in) :: vol

    !> Summation parameter
    real(dp), intent(in) :: eta

    !> Estimated error of the summation.
    real(dp) :: err

    err = (pi**1.5_dp * c6sum / vol) * eta * (1.0_dp/rr**4&
        & + 1.0_dp/(rr**2 * eta**2) + 1.0_dp / (2.0_dp * eta**4)) * erfcwrap(rr/eta)

  end function getDispRealError


  !> Returns the error of the reciprocal space summation for a certain cutoff
  function getDispReciprocalError(gg, c6sum, eta) result(err)

    !> Cutoff radius
    real(dp), intent(in) :: gg

    !> VdW coefficient sum.
    real(dp), intent(in) :: c6sum

    !> Summation parameter
    real(dp), intent(in) :: eta

    !> Estimated error of the summation.
    real(dp) :: err

    err = c6sum/(6.0_dp * sqrt(pi)) * (1.0_dp / eta**6)&
        & * (gg * eta * exp(-1.0_dp * (0.5_dp*gg*eta)**2) + sqrt(pi) * erfcwrap(0.5_dp * gg * eta))

  end function getDispReciprocalError

end module dftbp_dispcommon
