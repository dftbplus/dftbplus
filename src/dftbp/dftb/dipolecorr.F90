!--------------------------------------------------------------------------------------------------!
!  DFTB+: general package for performing fast atomistic simulations                                !
!  Copyright (C) 2006 - 2023  DFTB+ developers group                                               !
!                                                                                                  !
!  See the LICENSE file for terms of usage and distribution.                                       !
!--------------------------------------------------------------------------------------------------!

#:include 'common.fypp'

!> Implements a dipole correction for neutral supercell calculations
!>
!> Unchecked constraints:
!>
!> * Cell must be orthogonal
!> * Dipole correction layer is orthogonal to the z-axis
!> * Cell is neutral
!>
module dftbp_dftb_dipolecorr
  use dftbp_common_accuracy, only : dp
  use dftbp_common_constants, only : pi
  use dftbp_math_simplealgebra, only : determinant33
  implicit none


  type :: TDipoleCorrInput
    real(dp) :: z0
  end type


  type :: TDipoleCorr
   private

   real(dp) :: z0_ = 0.0_dp
   real(dp) :: cellHeight_ = 0.0_dp
   real(dp) :: cellVol_ = 0.0_dp
   real(dp), allocatable :: zCoords_(:)
   real(dp) :: dipoleZ_
  contains
    procedure :: updateLatVecs => TDipoleCorr_updateLatVecs
    procedure :: updateCoords => TDipoleCorr_updateCoords
    procedure :: updateCharges => TDipoleCorr_updateCharges
    procedure :: getShiftPerAtom => TDipoleCorr_getShiftPerAtom
    procedure :: getPotential => TDipoleCorr_getPotential
    procedure :: addEnergyPerAtom => TDipoleCorr_addEnergyPerAtom
    procedure :: addForceDc => TDipoleCorr_addForceDc
  end type


contains

  subroutine TDipoleCorr_init(this, inp)
    type(TDipoleCorr), intent(out) :: this
    type(TDipoleCorrInput), intent(in) :: inp

    this%z0_ = inp%z0

  end subroutine TDipoleCorr_init


  subroutine TDipoleCorr_updateLatVecs(this, latVecs)
    class(TDipoleCorr), intent(inout) :: this
    real(dp), intent(in) :: latVecs(:,:)

    this%cellHeight_ = latVecs(3, 3)
    this%cellVol_ = determinant33(latVecs)

  end subroutine TDipoleCorr_updateLatVecs


  subroutine TDipoleCorr_updateCoords(this, coords0)
    class(TDipoleCorr), intent(inout) :: this
    real(dp), intent(in) :: coords0(:,:)

    this%zCoords_ = foldedZCoords_(coords0, this%z0_, this%cellHeight_)

  end subroutine TDipoleCorr_updateCoords


  subroutine TDipoleCorr_updateCharges(this, dQAtom)
    class(TDipoleCorr), intent(inout) :: this
    real(dp), intent(in) :: dQAtom(:)

    integer :: nAtom

    nAtom = size(dQAtom)
    this%dipoleZ_ = sum((this%zCoords_ - (this%z0_ + this%cellHeight_ / 2.0_dp)) * dQAtom)

  end subroutine TDipoleCorr_updateCharges


  subroutine TDipoleCorr_getShiftPerAtom(this, shiftPerAtom)
    class(TDipoleCorr), intent(in) :: this
    real(dp), intent(out) :: shiftPerAtom(:)

    shiftPerAtom(:) = dipoleCorrectionPot_(this%zCoords_, this%z0_, this%dipoleZ_, this%cellVol_,&
        & this%cellHeight_)

  end subroutine TDipoleCorr_getShiftPerAtom


  subroutine TDipoleCorr_getPotential(this, coords, potentials)
    class(TDipoleCorr), intent(in) :: this
    real(dp), intent(in) :: coords(:,:)
    real(dp), intent(out) :: potentials(:)

    real(dp), allocatable :: zCoords(:)

    zCoords = foldedZCoords_(coords, this%z0_, this%cellHeight_)
    potentials(:) = dipoleCorrectionPot_(zCoords, this%z0_, this%dipoleZ_, this%cellVol_,&
        & this%cellHeight_)

  end subroutine TDipoleCorr_getPotential


  subroutine TDipoleCorr_addEnergyPerAtom(this, deltaQAtom, energyPerAtom)
    class(TDipoleCorr), intent(in) :: this
    real(dp), intent(in) :: deltaQAtom(:)
    real(dp), intent(inout) :: energyPerAtom(:)

    energyPerAtom(:) = energyPerAtom &
        & + 0.5_dp * deltaQAtom * dipoleCorrectionPot_(this%zCoords_, this%z0_, this%dipoleZ_,&
        & this%cellVol_, this%cellHeight_)

  end subroutine TDipoleCorr_addEnergyPerAtom


  subroutine TDipoleCorr_addForceDc(this, forces, deltaQAtom)
    class(TDipoleCorr), intent(in) :: this
    real(dp), intent(inout) :: forces(:,:)
    real(dp), intent(in) :: deltaQAtom(:)

    forces(3, :) = forces(3, :) + 4.0_dp * pi / this%cellVol_ * this%dipoleZ_ * deltaQAtom

  end subroutine TDipoleCorr_addForceDc


  pure function foldedZCoords_(coords, z0, cellHeight) result(zCoords)
    real(dp), intent(in) :: coords(:,:)
    real(dp), intent(in) :: z0, cellHeight
    real(dp) :: zCoords(size(coords, dim=2))

    zCoords(:) = modulo(coords(3, :) - z0, cellHeight) + z0

  end function foldedZCoords_


  pure function dipoleCorrectionPot_(zCoords, z0, dipoleZ, cellVol, cellHeight) result(pot)
    real(dp), intent(in) :: zCoords(:)
    real(dp), intent(in) :: z0, dipoleZ, cellVol, cellHeight
    real(dp) :: pot(size(zCoords))

    pot(:) = 4.0_dp * pi / cellVol * dipoleZ * (zCoords - (z0 + cellHeight / 2.0_dp))

  end function dipoleCorrectionPot_

end module dftbp_dftb_dipolecorr
