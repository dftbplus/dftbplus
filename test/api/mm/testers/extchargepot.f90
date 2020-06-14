!--------------------------------------------------------------------------------------------------!
!  DFTB+: general package for performing fast atomistic simulations                                !
!  Copyright (C) 2006 - 2020  DFTB+ developers group                                               !
!                                                                                                  !
!  See the LICENSE file for terms of usage and distribution.                                       !
!--------------------------------------------------------------------------------------------------!

!> Module containing some routines to calculate external point charges related quantities, as used
!> by a few API tests.
!>
!> NOTE: This is not a production module, it should only used to drive the API tests.
!> The tasks implemented here should be provided by the MM-driver in a production code.
!>
module extchargepot
  implicit none
  private

  public :: getPointChargePotential, getPointChargeGradients


  !> double precision kind
  integer, parameter :: dp = kind(1.0d0)


contains

  !> Calculate the potential and its first derivatives at DFTB atoms due to external electrostatic
  !> charges (perhaps from a surrounding MM region)
  subroutine getPointChargePotential(coordsMm, chargesMm, coordsQm, extPot, extPotGrad)

    !> Coordinates of the external charges (xyz,:nAtomMm) in atomic units
    real(dp), intent(in) :: coordsMm(:,:)

    !> Charges of MM region atoms, in atomic units (:nAtomMm)
    real(dp), intent(in) :: chargesMm(:)

    !> Coordinates of DFTB QM atoms (xyz,:nAtomQm) in atomic units
    real(dp), intent(in) :: coordsQm(:,:)

    !> Potentials at DFTB atomic sites (:nAtomQm)
    real(dp), intent(out) :: extPot(:)

    !> Gradient of potentials with respect to DFTB atom displacement (xyz,:nAtomQm)
    real(dp), intent(out) :: extPotGrad(:,:)

    real(dp) :: atomPosQm(3), atomPosMm(3)
    real(dp) :: chargeMm, dist
    integer :: nAtomQm, nAtomMm
    integer :: iAtQm, iAtMm

    nAtomQm = size(coordsQm, dim=2)
    nAtomMm = size(coordsMm, dim=2)
    extPot(:) = 0.0_dp
    extPotGrad(:,:) = 0.0_dp
    do iAtQm = 1, nAtomQm
      atomPosQm(:) = coordsQm(:, iAtQm)
      do iAtMm = 1, nAtomMm
        atomPosMm(:) = coordsMm(1:3, iAtMm)
        chargeMm = chargesMm(iAtMm)
        dist = sqrt(sum((atomPosQm - atomPosMm)**2))
        extPot(iAtQm) = extPot(iAtQm) - chargeMm / dist
        extPotGrad(:, iAtQm) = extPotGrad(:, iAtQm) + chargeMm * (atomPosQm - atomPosMm) / dist**3
      end do
    end do

  end subroutine getPointChargePotential


  !> Calculate the gradient of the electrostatic energy wrt to external charges in the field from
  !> the charges of DFTB atoms
  subroutine getPointChargeGradients(coordsQm, chargesQm, coordsMm, chargesMm, gradients)

    !> Coordinates of the DFTB QM atoms (xyz,:nAtomQm) in atomic units
    real(dp), intent(in) :: coordsQm(:,:)

    !> Charges of QM DFTB region atoms, in atomic units (:nAtomMm)
    real(dp), intent(in) :: chargesQm(:)

    !> Coordinates of the external charges (xyz,:nAtomMm) in atomic units
    real(dp), intent(in) :: coordsMm(:,:)

    !> Charges of MM region atoms, in atomic units (:nAtomMm)
    real(dp), intent(in) :: chargesMm(:)

    !> Gradient of potentials with respect to MM atom displacement (xyz,:nAtomMm)
    real(dp), intent(out) :: gradients(:,:)

    real(dp) :: atomPosQm(3), atomPosMm(3)
    real(dp) :: chargeQm, chargeMm, dist
    integer :: nAtomQm, nAtomMm
    integer :: iAtQm, iAtMm

    gradients(:,:) = 0.0_dp
    nAtomQm = size(coordsQm, dim=2)
    nAtomMm = size(coordsMm, dim=2)
    do iAtMm = 1, nAtomMm
      atomPosMm(:) = coordsMm(:, iAtMm)
      chargeMm = chargesMm(iAtMm)
      do iAtQm = 1, nAtomQm
        atomPosQm(:) = coordsQm(:, iAtQm)
        chargeQm = chargesQm(iAtQm)
        dist = sqrt(sum((atomPosQm - atomPosMm)**2))
        gradients(:, iAtMm) = gradients(:, iAtMm) &
            & + chargeQm * chargeMm * (atomPosMm - atomPosQm) / dist**3
      end do
    end do

  end subroutine getPointChargeGradients

end module extchargepot
