!--------------------------------------------------------------------------------------------------!
!  DFTB+: general package for performing fast atomistic simulations                                !
!  Copyright (C) 2006 - 2025  DFTB+ developers group                                               !
!                                                                                                  !
!  See the LICENSE file for terms of usage and distribution.                                       !
!--------------------------------------------------------------------------------------------------!

#:include 'common.fypp'
#:include 'error.fypp'

!> Implements a dipole correction for neutral supercell calculations
module dftbp_dftb_dipolecorr
  use dftbp_common_accuracy, only : dp, elecTolMax
  use dftbp_common_constants, only : pi
  use dftbp_common_status, only : TStatus
  use dftbp_math_simplealgebra, only : determinant33
  implicit none

  private
  public :: TDipoleCorr, TDipoleCorr_init, TDipoleCorrInput, ensureCorrRequirements

  !> Input settings for correction
  type :: TDipoleCorrInput

    !> Index of the lattice vector, which is normal to the slab
    integer :: iNormalVec = 3

    !> Index of the vector component, which represents the normal direction
    integer :: iNormalComp = 3

    !> Position of the dipole correction layer
    real(dp) :: z0 = 0.0_dp

  end type TDipoleCorrInput


  !> Type containing correction parameters
  type, extends(TDipoleCorrInput) :: TDipoleCorr
    private
    !> Cell length along normal direction
    real(dp) :: cellHeight_ = 0.0_dp
    !> Cell volume
    real(dp) :: cellVol_ = 0.0_dp
    !> Z components of atomic coodinates
    real(dp), allocatable :: zCoords_(:)
    !> Z component of dipole moment on plane
    real(dp) :: dipoleZ_ = 0.0_dp
  contains
    procedure :: updateLatVecs => TDipoleCorr_updateLatVecs
    procedure :: updateCoords => TDipoleCorr_updateCoords
    procedure :: updateCharges => TDipoleCorr_updateCharges
    procedure :: getShiftPerAtom => TDipoleCorr_getShiftPerAtom
    procedure :: getPotential => TDipoleCorr_getPotential
    procedure :: addEnergyPerAtom => TDipoleCorr_addEnergyPerAtom
    procedure :: addForceDc => TDipoleCorr_addForceDc
  end type TDipoleCorr

contains


  !> Check fulfillment of the requirements for using the dipole correction as implemented
  subroutine ensureCorrRequirements(this, errStatus, latVecs, netCharge)

    !> Instance
    class(TDipoleCorrInput), intent(in) :: this

    !> Error status
    type(TStatus), intent(out) :: errStatus

    !> Lattice vectors
    real(dp), intent(in), optional :: latVecs(:,:)

    !> System charge
    real(dp), intent(in), optional :: netCharge

    real(dp) :: proj(3), vec(3)

    @:ASSERT(present(latVecs) .or. present(netCharge))

    if (present(latVecs)) then
      vec(:) = 0.0_dp
      vec(this%iNormalComp) = 1.0_dp
      proj(:) = abs(matmul(latVecs, vec))
      if (count(proj > 1E-12_dp) /= 1 .or. abs(proj(this%iNormalComp)) < 1E-12_dp) then
        @:RAISE_ERROR(errStatus, -1, "Dipole correction only applicable if only one lattice vector&
            & (the slab normal vector) has non-zero z-component")
      end if
    end if

    if (present(netCharge)) then
      if (abs(netCharge) > elecTolMax) then
        @:RAISE_ERROR(errStatus, -1, "Dipole correction only applicable for charge neutral systems")
      end if
    end if

  end subroutine ensureCorrRequirements


  !> Initialise type
  subroutine TDipoleCorr_init(this, inp)

    !> Instance
    type(TDipoleCorr), intent(out) :: this

    !> Input settings
    type(TDipoleCorrInput), intent(in) :: inp

    this%iNormalVec = inp%iNormalVec
    this%iNormalComp = inp%iNormalComp
    this%z0 = inp%z0

  end subroutine TDipoleCorr_init


  !> Modify lattice vectors of the system
  subroutine TDipoleCorr_updateLatVecs(this, latVecs, errStatus)

    !> Instance
    class(TDipoleCorr), intent(inout) :: this

    !> Lattice vectors
    real(dp), intent(in) :: latVecs(:,:)

    !> Error status
    type(TStatus), intent(out) :: errStatus

    call ensureCorrRequirements(this, errStatus, latvecs=latVecs)
    @:PROPAGATE_ERROR(errStatus)
    this%cellHeight_ = latVecs(this%iNormalComp, this%iNormalVec)
    this%cellVol_ = determinant33(latVecs)

  end subroutine TDipoleCorr_updateLatVecs


  !> Modify atomic coordinates of the system
  subroutine TDipoleCorr_updateCoords(this, coords0)

    !> Instance
    class(TDipoleCorr), intent(inout) :: this

    !> Central cell coordinates
    real(dp), intent(in) :: coords0(:,:)

    this%zCoords_ = foldedZCoords_(coords0, this%iNormalComp, this%z0, this%cellHeight_)

  end subroutine TDipoleCorr_updateCoords


  !> Modify atomic charges of the system
  subroutine TDipoleCorr_updateCharges(this, dQAtom, errStatus)

    !> Instance
    class(TDipoleCorr), intent(inout) :: this

    !> Atomic charges
    real(dp), intent(in) :: dQAtom(:)

    !> Error status
    type(TStatus), intent(out) :: errStatus

    integer :: nAtom

    call ensureCorrRequirements(this, errStatus, netCharge=sum(dQAtom))
    @:PROPAGATE_ERROR(errStatus)
    nAtom = size(dQAtom)
    this%dipoleZ_ = sum((this%zCoords_ - (this%z0 + this%cellHeight_ / 2.0_dp)) * dQAtom)

  end subroutine TDipoleCorr_updateCharges


  !> Return atomic site potential shifts
  subroutine TDipoleCorr_getShiftPerAtom(this, shiftPerAtom)

    !> Instance
    class(TDipoleCorr), intent(in) :: this

    !> Resulting shifts
    real(dp), intent(out) :: shiftPerAtom(:)

    shiftPerAtom(:) = dipoleCorrectionPot_(this%zCoords_, this%z0, this%dipoleZ_, this%cellVol_,&
        & this%cellHeight_)

  end subroutine TDipoleCorr_getShiftPerAtom


  !> Get z-component of electrostatic potential at specified coordinates
  subroutine TDipoleCorr_getPotential(this, coords0, potentials)

    !> Instance
    class(TDipoleCorr), intent(in) :: this

    !> Sites to obtain potential
    real(dp), intent(in) :: coords0(:,:)

    !> Resulting potentials at sites
    real(dp), intent(out) :: potentials(:)

    real(dp), allocatable :: zCoords_(:)

    zCoords_ = foldedZCoords_(coords0, this%iNormalComp, this%z0, this%cellHeight_)
    potentials(:) = dipoleCorrectionPot_(zCoords_, this%z0, this%dipoleZ_, this%cellVol_,&
        & this%cellHeight_)

  end subroutine TDipoleCorr_getPotential


  !> Adds energy contributions to each atom from the correction
  subroutine TDipoleCorr_addEnergyPerAtom(this, deltaQAtom, energyPerAtom)

    !> Instance
    class(TDipoleCorr), intent(in) :: this

    !> Atomic charges
    real(dp), intent(in) :: deltaQAtom(:)

    !> Energy associated with each atom
    real(dp), intent(inout) :: energyPerAtom(:)

    energyPerAtom(:) = energyPerAtom &
        & + 0.5_dp * deltaQAtom * dipoleCorrectionPot_(this%zCoords_, this%z0, this%dipoleZ_,&
        & this%cellVol_, this%cellHeight_)

  end subroutine TDipoleCorr_addEnergyPerAtom


  !> Adds force contributions to each atom from the correction
  subroutine TDipoleCorr_addForceDc(this, forces, deltaQAtom)

    !> Instance
    class(TDipoleCorr), intent(in) :: this

    !> Forces acting on each atom
    real(dp), intent(inout) :: forces(:,:)

    !> Atomic charges
    real(dp), intent(in) :: deltaQAtom(:)

    forces(this%iNormalComp, :) = &
        & forces(this%iNormalComp, :) + 4.0_dp * pi / this%cellVol_ * this%dipoleZ_ * deltaQAtom

  end subroutine TDipoleCorr_addForceDc


  ! Internal routines

  !> Folds z component of atomic coordinates with respect to the correction plane
  pure function foldedZCoords_(coords0, iNormalComp, z0, cellHeight_) result(zCoords_)

    !> Atomic coordinates for central cell atoms
    real(dp), intent(in) :: coords0(:,:)

    !> Cartesian component to fold
    integer, intent(in) :: iNormalComp

    !> Plane position to fold with respect to
    real(dp), intent(in) :: z0

    !> Length of cell normal to plane
    real(dp), intent(in) :: cellHeight_

    real(dp) :: zCoords_(size(coords0, dim=2))

    zCoords_(:) = modulo(coords0(iNormalComp, :) - z0, cellHeight_) + z0

  end function foldedZCoords_


  !> Potential from dipole correction
  pure function dipoleCorrectionPot_(zCoords_, z0, dipoleZ_, cellVol_, cellHeight_) result(pot)

    !> Z part of charge coordinates
    real(dp), intent(in) :: zCoords_(:)

    !> Plane distance along z
    real(dp), intent(in) :: z0

    !> Dipole moment along z
    real(dp), intent(in) :: dipoleZ_

    !> Cell volume
    real(dp), intent(in) :: cellVol_

    !> Cell length
    real(dp), intent(in) :: cellHeight_

    !> Resulting potentials at atom sites
    real(dp) :: pot(size(zCoords_))

    pot(:) = 4.0_dp * pi / cellVol_ * dipoleZ_ * (zCoords_ - (z0 + cellHeight_ / 2.0_dp))

  end function dipoleCorrectionPot_

end module dftbp_dftb_dipolecorr
