!--------------------------------------------------------------------------------------------------!
!  DFTB+: general package for performing fast atomistic simulations                                !
!  Copyright (C) 2006 - 2020  DFTB+ developers group                                               !
!                                                                                                  !
!  See the LICENSE file for terms of usage and distribution.                                       !
!--------------------------------------------------------------------------------------------------!

#:include 'common.fypp'

!> Implementation of the electronegativity equilibration charge model used
!  for the charge scaling in DFT-D4.
!
!  This implementation is general enough to be used outside of DFT-D4.
!
module dftbp_encharges
  use dftbp_assert
  use dftbp_accuracy, only : dp
  use dftbp_constants, only : pi
  use dftbp_errorfunction, only : erfwrap
  use dftbp_coulomb, only : ewaldReal, ewaldReciprocal, derivStressEwaldRec
  use dftbp_blasroutines, only : hemv, gemv, gemm
  use dftbp_lapackroutines, only : symmatinv
  implicit none
  private

  public :: getEEQcharges

  real(dp), parameter :: sqrtpi = sqrt(pi)
  real(dp), parameter :: sqrt2pi = sqrt(2.0_dp/pi)

contains

  !> Generates full interaction matrix for Gaussian charge distributions.
  subroutine getCoulombMatrixCluster(nAtom, coords, species, gam, rad, aMat)

    !> number of atoms
    integer, intent(in) :: nAtom

    !> List of atomic coordinates.
    real(dp), intent(in) :: coords(:, :)

    !> Species of every atom.
    integer, intent(in) :: species(:)

    !> element-specific chemical hardnesses
    real(dp), intent(in) :: gam(:)

    !> element-specific charge widths / atomic radii
    real(dp), intent(in) :: rad(:)

    !> Interaction Matrix for each atom pair.
    real(dp), intent(out) :: aMat(:, :)

    integer :: iAt1, iSp1, iAt2f, iSp2
    real(dp) :: dist, vec(3), eta12

    aMat(:, :) = 0.0_dp

    !$OMP PARALLEL DEFAULT(NONE) &
    !$OMP SHARED(nAtom, species, gam, rad, coords, aMat) &
    !$OMP PRIVATE(iSp1, iAt2f, iSp2, vec, dist, eta12)
    !$OMP DO SCHEDULE(RUNTIME)
    do iAt1 = 1, nAtom
      iSp1 = species(iAt1)
      aMat(iAt1, iAt1) = gam(iSp1) + sqrt2pi / rad(iSp1)
      do iAt2f = iAt1 + 1, nAtom
        iSp2 = species(iAt2f)
        vec(:) = coords(:, iAt1) - coords(:, iAt2f)
        dist = sqrt(sum(vec**2))
        eta12 = 1.0_dp / sqrt(rad(iSp1)**2 + rad(iSp2)**2)
        aMat(iAt2f, iAt1) = erfwrap(eta12 * dist) / dist
        aMat(iAt1, iAt2f) = aMat(iAt2f, iAt1)
      end do
    end do
    !$OMP END DO
    !$OMP END PARALLEL

  end subroutine getCoulombMatrixCluster


  !> Derivatives of interaction matrix for Gaussian charge distributions.
  subroutine getCoulombDerivsCluster(nAtom, coords, species, rad, qVec, aMatdr, aTrace)

    !> Number of atoms.
    integer, intent(in) :: nAtom

    real(dp), intent(in) :: coords(:, :)

    !> Species of every atom.
    integer, intent(in) :: species(:)

    !> element-specific charge widths / atomic radii
    real(dp), intent(in) :: rad(:)

    !> List of charges on each atom.
    real(dp), intent(in) :: qVec(:)

    !> Contains the derivative on exit.
    real(dp), intent(out) :: aMatdr(:, :, :)

    !> Contains the `trace' derivative on exit.
    real(dp), intent(out) :: aTrace(:, :)

    integer :: iAt1, iAt2f, iSp1, iSp2
    real(dp) :: dist, vec(3), aTmp(3), arg, eta12

    aMatdr(:,:,:) = 0.0_dp
    aTrace(:,:) = 0.0_dp

    !$OMP PARALLEL DEFAULT(NONE) REDUCTION(+:aTrace, aMatdr) &
    !$OMP SHARED(nAtom, species, rad, coords, qVec) &
    !$OMP PRIVATE(iAt2f, iSp1, iSp2, vec, dist, aTmp, arg, eta12)
    !$OMP DO SCHEDULE(RUNTIME)
    do iAt1 = 1, nAtom
      iSp1 = species(iAt1)
      do iAt2f = 1, iAt1 - 1
        iSp2 = species(iAt2f)
        vec(:) = coords(:, iAt1) - coords(:, iAt2f)
        dist = sqrt(sum(vec**2))
        eta12 = 1.0_dp / sqrt(rad(iSp1)**2 + rad(iSp2)**2)
        arg = dist * eta12
        aTmp = vec * (2.0_dp * eta12 / sqrtpi * exp(-arg*arg) - erfwrap(arg)/dist) / dist**2
        aTrace(:, iAt1) = aTrace(:, iAt1) + aTmp * qVec(iAt2f)
        aTrace(:, iAt2f) = aTrace(:, iAt2f) - aTmp * qVec(iAt1)
        aMatdr(:, iAt1, iAt2f) = aMatdr(:, iAt1, iAt2f) + aTmp * qVec(iAt1)
        aMatdr(:, iAt2f, iAt1) = aMatdr(:, iAt2f, iAt1) - aTmp * qVec(iAt2f)
      end do
    end do
    !$OMP END DO
    !$OMP END PARALLEL

  end subroutine getCoulombDerivsCluster


  !> Generates full interaction matrix for Gaussian charge distributions.
  subroutine getCoulombMatrixPeriodic(nAtom, coords, species, nNeighbour, iNeighbour, neighDist2,&
      & img2CentCell, recPoint, gam, rad, alpha, volume, aMat)

    !> Nr. of atoms (without periodic images)
    integer, intent(in) :: nAtom

    !> Species of every atom.
    integer, intent(in) :: species(:)

    !> List of atomic coordinates (all atoms).
    real(dp), intent(in) :: coords(:, :)

    !> Nr. of neighbours for each atom
    integer, intent(in) :: nNeighbour(:)

    !> Neighbourlist.
    integer, intent(in) :: iNeighbour(0:, :)

    !> Square distances of the neighbours.
    real(dp), intent(in) :: neighDist2(0:, :)

    !> Mapping into the central cell.
    integer, intent(in) :: img2CentCell(:)

    !> Contains the points included in the reciprocal sum.
    !  The set should not include the origin or inversion related points.
    real(dp), intent(in) :: recPoint(:, :)

    !> element-specific chemical hardnesses
    real(dp), intent(in) :: gam(:)

    !> element-specific charge widths / atomic radii
    real(dp), intent(in) :: rad(:)

    !> Parameter for Ewald summation.
    real(dp), intent(in) :: alpha

    !> Volume of the real space unit cell.
    real(dp), intent(in) :: volume

    !> Matrix of 1/R values for each atom pair.
    real(dp), intent(out) :: aMat(:, :)

    aMat(:, :) = 0.0_dp

    ! Real space part of the Ewald sum.
    call addRealSpaceContribs(nAtom, coords, species, nNeighbour, iNeighbour, neighDist2,&
        & img2CentCell, gam, rad, alpha, aMat)

    ! Reciprocal space part of the Ewald sum.
    call addEwaldContribs(nAtom, coords, recPoint, alpha, volume, aMat)

  end subroutine getCoulombMatrixPeriodic


  !> Derivatives of interaction matrix for Gaussian charge distributions.
  subroutine getCoulombDerivsPeriodic(nAtom, coords, species, nNeighbour, iNeighbour, neighDist2,&
      & img2CentCell, recPoint, alpha, volume, rad, qVec, aMatdr, aMatdL, aTrace)

    !> Number of atoms
    integer, intent(in) :: nAtom

    !> List of atomic coordinates (all atoms).
    real(dp), intent(in) :: coords(:, :)

    !> Species of every atom.
    integer, intent(in) :: species(:)

    !> Nr. of neighbours for each atom
    integer, intent(in) :: nNeighbour(:)

    !> Neighbourlist.
    integer, intent(in) :: iNeighbour(0:, :)

    !> Square distances of the neighbours.
    real(dp), intent(in) :: neighDist2(0:, :)

    !> Mapping into the central cell.
    integer, intent(in) :: img2CentCell(:)

    !> Contains the points included in the reciprocal sum.
    !  The set should not include the origin or inversion related points.
    real(dp), intent(in) :: recPoint(:, :)

    !> Parameter for Ewald summation.
    real(dp), intent(in) :: alpha

    !> Volume of the real space unit cell.
    real(dp), intent(in) :: volume

    !> element-specific charge widths / atomic radii
    real(dp), intent(in) :: rad(:)

    !> List of charges on each atom
    real(dp), intent(in) :: qVec(:)

    !> Contains the derivative on exit.
    real(dp), intent(out) :: aMatdr(:, :, :)

    !> Contains the strain derivative on exit.
    real(dp), intent(out) :: aMatdL(:, :, :)

    !> Contains the `trace' derivative on exit.
    real(dp), intent(out) :: aTrace(:, :)

    @:ASSERT(volume > 0.0_dp)

    aMatdr(:, :, :) = 0.0_dp
    aMatdL(:, :, :) = 0.0_dp
    aTrace(:, :) = 0.0_dp

    ! d(1/R)/dr real space
    call addRealSpaceDerivs(nAtom, coords, species, nNeighbour, iNeighbour, neighDist2,&
        & img2CentCell, alpha, rad, qVec, aMatdr, aMatdL, aTrace)

    ! d(1/R)/dr reciprocal space
    call addEwaldDerivs(nAtom, coords, recPoint, alpha, volume, qVec, aMatdr, aMatdL, aTrace)

  end subroutine getCoulombDerivsPeriodic


  !> Reciprocal space contributions to interaction matrix.
  subroutine addEwaldContribs(nAtom, coords, recPoint, alpha, volume, aMat)

    !> Nr. of atoms (without periodic images)
    integer, intent(in) :: nAtom

    !> List of atomic coordinates (all atoms).
    real(dp), intent(in) :: coords(:, :)

    !> Contains the points included in the reciprocal sum.
    !  The set should not include the origin or inversion related points.
    real(dp), intent(in) :: recPoint(:, :)

    !> Parameter for Ewald summation.
    real(dp), intent(in) :: alpha

    !> Volume of the real space unit cell.
    real(dp), intent(in) :: volume

    !> Interaction Matrix for each atom pair.
    real(dp), intent(inout) :: aMat(:, :)

    integer :: iAt1, iAt2f
    real(dp) :: vec(3), rTerm

    @:ASSERT(volume > 0.0_dp)

    ! Reciprocal space part of the Ewald sum.
    !$OMP PARALLEL DEFAULT(NONE) REDUCTION(+:aMat) PRIVATE(iAt2f, vec, rTerm) &
    !$OMP SHARED(nAtom, alpha, coords, recPoint, volume)
    !$OMP DO SCHEDULE(RUNTIME)
    do iAt1 = 1, nAtom
      aMat(iAt1, iAt1) = aMat(iAt1, iAt1) - alpha / sqrt(pi) + pi / (volume * alpha**2)
      do iAt2f = iAt1, nAtom
        vec(:) = coords(:, iAt1)-coords(:, iAt2f)
        rTerm = ewaldReciprocal(vec, recPoint, alpha, volume) - pi / (volume * alpha**2)
        aMat(iAt2f, iAt1) = aMat(iAt2f, iAt1) + rTerm
        aMat(iAt1, iAt2f) = aMat(iAt1, iAt2f) + rTerm
      end do
    end do
    !$OMP END DO
    !$OMP END PARALLEL

  end subroutine addEwaldContribs


  !> Reciprocal space contributions to derivative of interaction matrix.
  subroutine addEwaldDerivs(nAtom, coords, recPoint, alpha, volume, qVec, aMatdr, aMatdL, aTrace)

    !> Number of atoms
    integer, intent(in) :: nAtom

    !> List of atomic coordinates (all atoms).
    real(dp), intent(in) :: coords(:, :)

    !> Contains the points included in the reciprocal sum.
    !  The set should not include the origin or inversion related points.
    real(dp), intent(in) :: recPoint(:, :)

    !> Parameter for Ewald summation.
    real(dp), intent(in) :: alpha

    !> Volume of the real space unit cell.
    real(dp), intent(in) :: volume

    !> List of charges on each atom
    real(dp), intent(in) :: qVec(:)

    !> Contains the derivative on exit.
    real(dp), intent(inout) :: aMatdr(:, :, :)

    !> Contains the strain derivative on exit.
    real(dp), intent(inout) :: aMatdL(:, :, :)

    !> Contains the `trace' derivative on exit.
    real(dp), intent(inout) :: aTrace(:, :)

    integer :: iAt1, iAt2f
    real(dp) :: vec(3), aTmp(3), sigma(3, 3)

    !$OMP PARALLEL DEFAULT(NONE) REDUCTION(+:aTrace, aMatdr, aMatdL) &
    !$OMP SHARED(nAtom, coords, recPoint, alpha, volume, qVec) &
    !$OMP PRIVATE(iAt2f, vec, aTmp, sigma)
    !$OMP DO SCHEDULE(RUNTIME)
    do iAt1 = 1, nAtom
      do iAt2f = iAt1, nAtom
        vec(:) = coords(:, iAt1) - coords(:, iAt2f)
        call derivStressEwaldRec(vec, recPoint, alpha, volume, aTmp, sigma)
        aTrace(:, iAt1) = aTrace(:, iAt1) + aTmp * qVec(iAt2f)
        aTrace(:, iAt2f) = aTrace(:, iAt2f) - aTmp * qVec(iAt1)
        aMatdr(:, iAt1, iAt2f) = aMatdr(:, iAt1, iAt2f) + aTmp * qVec(iAt1)
        aMatdr(:, iAt2f, iAt1) = aMatdr(:, iAt2f, iAt1) - aTmp * qVec(iAt2f)
        aMatdL(:, :, iAt1) = aMatdL(:, :, iAt1) + sigma * qVec(iAt2f)
        if (iAt1 /= iAt2f) then
          aMatdL(:, :, iAt2f) = aMatdL(:, :, iAt2f) + sigma * qVec(iAt1)
        end if
      end do
    end do
    !$OMP END DO
    !$OMP END PARALLEL

  end subroutine addEwaldDerivs


  !> Real space contributions to interaction matrix.
  subroutine addRealSpaceContribs(nAtom, coords, species, nNeighbour, iNeighbour, neighDist2,&
      & img2CentCell, gam, rad, alpha, aMat)

    !> Nr. of atoms (without periodic images)
    integer, intent(in) :: nAtom

    !> Species of every atom.
    integer, intent(in) :: species(:)

    !> List of atomic coordinates (all atoms).
    real(dp), intent(in) :: coords(:, :)

    !> Nr. of neighbours for each atom
    integer, intent(in) :: nNeighbour(:)

    !> Neighbourlist.
    integer, intent(in) :: iNeighbour(0:, :)

    !> Square distances of the neighbours.
    real(dp), intent(in) :: neighDist2(0:, :)

    !> Mapping into the central cell.
    integer, intent(in) :: img2CentCell(:)

    !> Parameter for Ewald summation.
    real(dp), intent(in) :: alpha

    !> element-specific chemical hardnesses
    real(dp), intent(in) :: gam(:)

    !> element-specific charge widths / atomic radii
    real(dp), intent(in) :: rad(:)

    !> Interaction Matrix for each atom pair.
    real(dp), intent(inout) :: aMat(:, :)

    real(dp) :: dist, eta12, rTerm
    integer :: iAt1, iAt2, iAt2f, iSp1, iSp2, iNeigh

    !$OMP PARALLEL DEFAULT(NONE) REDUCTION(+:aMat) &
    !$OMP SHARED(nAtom, species, gam, rad, nNeighbour, iNeighbour, img2CentCell) &
    !$OMP SHARED(neighDist2, alpha) &
    !$OMP PRIVATE(iNeigh, iAt2, iAt2f, iSp1, iSp2, dist, rTerm, eta12)
    !$OMP DO SCHEDULE(RUNTIME)
    do iAt1 = 1, nAtom
      iSp1 = species(iAt1)
      aMat(iAt1, iAt1) = aMat(iAt1, iAt1) + gam(iSp1) + sqrt2pi/rad(iSp1)
      do iNeigh = 1, nNeighbour(iAt1)
        iAt2 = iNeighbour(iNeigh, iAt1)
        iAt2f = img2CentCell(iAt2)
        iSp2 = species(iAt2f)
        dist = sqrt(neighDist2(iNeigh, iAt1))
        eta12 = 1.0_dp/sqrt(rad(iSp1)**2 + rad(iSp2)**2)
        rTerm = (erfwrap(eta12 * dist) - erfwrap(alpha * dist)) / dist
        aMat(iAt2f, iAt1) = aMat(iAt2f, iAt1) + rTerm
        aMat(iAt1, iAt2f) = aMat(iAt1, iAt2f) + rTerm
      end do
    end do
    !$OMP END DO
    !$OMP END PARALLEL

  end subroutine addRealSpaceContribs


  !> Real space contributions to derivative of interaction matrix.
  subroutine addRealSpaceDerivs(nAtom, coords, species, nNeighbour, iNeighbour, neighDist2,&
      & img2CentCell, alpha, rad, qVec, aMatdr, aMatdL, aTrace)

    !> Nr. of atoms (without periodic images)
    integer, intent(in) :: nAtom

    !> Species of every atom.
    integer, intent(in) :: species(:)

    !> List of atomic coordinates (all atoms).
    real(dp), intent(in) :: coords(:, :)

    !> Nr. of neighbours for each atom
    integer, intent(in) :: nNeighbour(:)

    !> Neighbourlist.
    integer, intent(in) :: iNeighbour(0:, :)

    !> Square distances of the neighbours.
    real(dp), intent(in) :: neighDist2(0:, :)

    !> Mapping into the central cell.
    integer, intent(in) :: img2CentCell(:)

    !> Parameter for Ewald summation.
    real(dp), intent(in) :: alpha

    !> element-specific charge widths / atomic radii
    real(dp), intent(in) :: rad(:)

    !> List of charges on each atom.
    real(dp), intent(in) :: qVec(:)

    !> Contains the derivative on exit.
    real(dp), intent(inout) :: aMatdr(:, :, :)

    !> Contains the strain derivative on exit.
    real(dp), intent(inout) :: aMatdL(:, :, :)

    !> Contains the `trace' derivative on exit.
    real(dp), intent(inout) :: aTrace(:, :)

    real(dp) :: dist, vec(3), aTmp(3), arg, eta12, ewl, sigma(3, 3)
    integer :: iAt1, iAt2, iAt2f, iSp1, iSp2, iNeigh

    !$OMP PARALLEL DEFAULT(NONE) REDUCTION(+:aTrace, aMatdr, aMatdL) &
    !$OMP SHARED(nAtom, species, nNeighbour, iNeighbour, coords, img2CentCell) &
    !$OMP SHARED(neighDist2, rad, qVec, alpha) &
    !$OMP PRIVATE(iAt2, iAt2f, iSp1, iSp2, vec, dist, aTmp, arg, ewl, eta12, sigma)
    !$OMP DO SCHEDULE(RUNTIME)
    do iAt1 = 1, nAtom
      iSp1 = species(iAt1)
      do iNeigh = 1, nNeighbour(iAt1)
        iAt2 = iNeighbour(iNeigh, iAt1)
        vec(:) = coords(:, iAt1) - coords(:, iAt2)
        iAt2f = img2CentCell(iAt2)
        iSp2 = species(iAt2f)
        dist = sqrt(neighDist2(iNeigh, iAt1))
        eta12 = 1.0_dp/sqrt(rad(iSp1)**2 + rad(iSp2)**2)
        arg = dist * eta12
        ewl = dist * alpha
        aTmp = ((2*eta12 / sqrtpi * exp(-arg *arg) - erfwrap(arg) / dist) &
            & - (2*alpha / sqrtpi * exp(-ewl *ewl) - erfwrap(ewl) / dist)) * vec / dist**2
        aTrace(:, iAt1) = aTrace(:, iAt1) + aTmp * qVec(iAt2f)
        aTrace(:, iAt2f) = aTrace(:, iAt2f) - aTmp * qVec(iAt1)
        aMatdr(:, iAt1, iAt2f) = aMatdr(:, iAt1, iAt2f) + aTmp * qVec(iAt1)
        aMatdr(:, iAt2f, iAt1) = aMatdr(:, iAt2f, iAt1) - aTmp * qVec(iAt2f)
        sigma(:,:) = spread(aTmp, 1, 3) * spread(vec, 2, 3)
        aMatdL(:, :, iAt1) = aMatdL(:, :, iAt1) + sigma * qVec(iAt2f)
        if (iAt1 /= iAt2f) then
          aMatdL(:, :, iAt2f) = aMatdL(:, :, iAt2f) + sigma * qVec(iAt1)
        end if
      end do
    end do
    !$OMP END DO
    !$OMP END PARALLEL

  end subroutine addRealSpaceDerivs


  !> Electronegativity equilibration charge model
  subroutine getEEQCharges(nAtom, coords, species, charge, nNeighbour, iNeighbour, neighDist2,&
      & img2CentCell, recPoint, alpha, volume, chi, kcn, gam, rad, cn, dcndr, dcndL, energies,&
      & gradients, stress, qAtom, dqdr, dqdL)

    !> number of atoms
    integer, intent(in) :: nAtom

    !> List of atomic coordinates.
    real(dp), intent(in) :: coords(:, :)

    !> Species of every atom.
    integer, intent(in) :: species(:)

    !> Total molecular charge.
    real(dp), intent(in) :: charge

    !> Nr. of neighbours for each atom
    integer, intent(in) :: nNeighbour(:)

    !> Neighbourlist.
    integer, intent(in) :: iNeighbour(0:, :)

    !> Square distances of the neighbours.
    real(dp), intent(in) :: neighDist2(0:, :)

    !> Mapping into the central cell.
    integer, intent(in) :: img2CentCell(:)

    !> Contains the points included in the reciprocal sum.
    !  The set should not include the origin or inversion related points.
    real(dp), allocatable, intent(in) :: recPoint(:, :)

    !> Parameter for Ewald summation.
    real(dp), intent(in) :: alpha

    !> Volume of the real space unit cell.
    real(dp), intent(in) :: volume

    !> element-specific electronegativity
    real(dp), intent(in) :: chi(:)

    !> element-specific chemical hardnesses
    real(dp), intent(in) :: gam(:)

    !> element-specific CN scaling constant
    real(dp), intent(in) :: kcn(:)

    !> element-specific charge widths / atomic radii
    real(dp), intent(in) :: rad(:)

    !> Error function coordination number.
    real(dp), intent(in) :: cn(:)

    !> Derivative of the CN with respect to the Cartesian coordinates.
    real(dp), intent(in) :: dcndr(:, :, :)

    !> Derivative of the CN with respect to strain deformations.
    real(dp), intent(in) :: dcndL(:, :, :)

    !> Updated energy vector at return
    real(dp), intent(inout), optional :: energies(:)

    !> Updated gradient vector at return
    real(dp), intent(inout), optional :: gradients(:, :)

    !> Updated stress tensor at return
    real(dp), intent(inout), optional :: stress(:, :)

    !> Atomic partial charges.
    real(dp), intent(out), optional :: qAtom(:)

    !> Derivative of partial charges w.r.t. cartesian coordinates.
    real(dp), intent(out), optional :: dqdr(:, :, :)

    !> Derivative of partial charges w.r.t. strain deformations.
    real(dp), intent(out), optional :: dqdL(:, :, :)

    real(dp), parameter :: small = 1.0e-14_dp

    !> Matrix of 1/R values for each atom pair.
    real(dp), allocatable :: aMat(:, :)
    real(dp), allocatable :: xVec(:)
    real(dp), allocatable :: xFac(:)
    real(dp), allocatable :: qVec(:)

    real(dp), allocatable :: aInv(:, :)

    real(dp), allocatable :: aMatdr(:, :, :)
    real(dp), allocatable :: aMatdL(:, :, :)
    real(dp), allocatable :: aTrace(:, :)

    logical :: tPeriodic
    integer :: nDim
    integer :: iAt1, iSp1, ii, jj
    real(dp) :: tmp

    tPeriodic = allocated(recPoint)

    nDim = nAtom + 1
    allocate(xVec(nDim), xFac(nAtom), aMatdr(3, nAtom, nDim), aMat(nDim, nDim), &
        & aInv(nDim, nDim), qVec(nDim), aMatdL(3, 3, nDim), aTrace(3, nAtom))

    ! Step 1: contruct RHS for linear system
    do iAt1 = 1, nAtom
      iSp1 = species(iAt1)
      tmp = kcn(iSp1) / (sqrt(cn(iAt1)) + small)
      xVec(iAt1) = -chi(iSp1) + tmp*cn(iAt1)
      xFac(iAt1) = -0.5_dp*tmp
    end do
    xVec(nDim) = charge

    ! Step 2: construct interaction matrix
    if (tPeriodic) then
      call getCoulombMatrixPeriodic(nAtom, coords, species, nNeighbour, iNeighbour, neighDist2,&
          & img2CentCell, recPoint, gam, rad, alpha, volume, aMat)
    else
      call getCoulombMatrixCluster(nAtom, coords, species, gam, rad, aMat)
    end if
    aMat(nDim, 1:nAtom) = 1.0_dp
    aMat(1:nAtom, nDim) = 1.0_dp
    aMat(nDim, nDim) = 0.0_dp

    ! Step 3: invert linear system
    aInv(:, :) = aMat
    call symmatinv(aInv)
    qVec(:) = 0.0_dp
    call hemv(qVec, aInv, xVec)

    ! Step 4: return atom resolved energies if requested, xVec is scratch
    if (present(energies)) then
      call hemv(xVec, aMat, qVec, alpha=0.5_dp, beta=-1.0_dp)
      energies(:) = energies(:) + xVec(:nAtom) * qVec(:nAtom)
    end if
    deallocate(xVec) ! free xVec to avoid using it later

    ! Step 5: get derivative of interaction matrix
    if (tPeriodic) then
      call getCoulombDerivsPeriodic(nAtom, coords, species, nNeighbour, iNeighbour, neighDist2,&
          & img2CentCell, recPoint, alpha, volume, rad, qVec, aMatdr, aMatdL, aTrace)
    else
      call getCoulombDerivsCluster(nAtom, coords, species, rad, qVec, aMatdr, aTrace)
    end if
    do iAt1 = 1, nAtom
      aMatdr(:, :, iAt1) = aMatdr(:, :, iAt1) + dcndr(:, :, iAt1) * xFac(iAt1)
      aMatdL(:, :, iAt1) = aMatdL(:, :, iAt1) + dcndL(:, :, iAt1) * xFac(iAt1)
    end do

    if (present(gradients)) then
      call gemv(gradients, aMatdr, qVec, beta=1.0_dp)
    end if

    if (present(stress)) then
      call gemv(stress, aMatdL, qVec, beta=1.0_dp)
    end if

    do iAt1 = 1, nAtom
      aMatdr(:, iAt1, iAt1) = aMatdr(:, iAt1, iAt1) + aTrace(:, iAt1)
    end do

    if (present(dqdr) .or. present(dqdL)) then
      ! Symmetrise inverted matrix
      do ii = 1, nDim
        aInv(ii, ii+1:nDim) = aInv(ii+1:nDim, ii)
      end do
    end if

    if (present(dqdr)) then
      dqdr(:, :, :) = 0.0_dp
      call gemm(dqdr, aMatdr, aInv, alpha=-1.0_dp)
    end if

    if (present(dqdL)) then
      dqdL(:, :, :) = 0.0_dp
      call gemm(dqdL, aMatdL, aInv, alpha=-1.0_dp)
    end if

    if (present(qAtom)) then
      qAtom(:) = qVec(:nAtom)
    end if

  end subroutine getEEQCharges


end module dftbp_encharges
