!--------------------------------------------------------------------------------------------------!
!  DFTB+: general package for performing fast atomistic simulations                                !
!  Copyright (C) 2018  DFTB+ developers group                                                      !
!                                                                                                  !
!  See the LICENSE file for terms of usage and distribution.                                       !
!--------------------------------------------------------------------------------------------------!

!> Dummy generator for an external potential given by point charges
module extchargepotgen
  use dftbp_qdepextpotgen
  use extchargepot, only : getPointChargePotential, getPointChargeGradients
  implicit none
  private

  public :: TExtChargePotGen, TExtChargePotGen_init

  integer, parameter :: dp = kind(1.0d0)


  type, extends(TQDepExtPotGen) :: TExtChargePotGen
    private
    real(dp), allocatable :: atomCoords(:,:)
    real(dp), allocatable :: extChargeCoords(:,:)
    real(dp), allocatable :: extCharges(:)
    real(dp), allocatable :: extPotGrad(:,:)
  contains
    procedure :: getExternalPot => TExtChargePotGen_getExternalPot
    procedure :: getExternalPotGrad => TExtChargePotGen_getExternalPotGrad
  end type TExtChargePotGen

contains

  subroutine TExtChargePotGen_init(this, atomCoords, extChargeCoords, extCharges)
    type(TExtChargePotGen), intent(out) :: this
    real(dp), intent(in) :: atomCoords(:,:)
    real(dp), intent(in) :: extChargeCoords(:,:)
    real(dp), intent(in) :: extCharges(:)

    this%atomCoords = atomCoords
    this%extChargeCoords = extChargeCoords
    this%extCharges = extCharges
    allocate(this%extPotGrad(3, size(atomCoords, dim=2)))

  end subroutine TExtChargePotGen_init


  subroutine TExtChargePotGen_getExternalPot(this, chargePerAtom, chargePerShell, extPotAtom,&
      & extPotShell)
    class(TExtChargePotGen), intent(inout) :: this
    real(dp), intent(in) :: chargePerAtom(:)
    real(dp), intent(in) :: chargePerShell(:,:)
    real(dp), intent(out) :: extPotAtom(:)
    real(dp), intent(out) :: extPotShell(:,:)

    call getPointChargePotential(this%extChargeCoords, -this%extCharges, this%atomCoords,&
        & extPotAtom, this%extPotGrad)
    extPotShell(:,:) = 0.0_dp

  end subroutine TExtChargePotGen_getExternalPot


  subroutine TExtChargePotGen_getExternalPotGrad(this, chargePerAtom, chargePerShell, extPotGrad)
    class(TExtChargePotGen), intent(inout) :: this
    real(dp), intent(in) :: chargePerAtom(:)
    real(dp), intent(in) :: chargePerShell(:,:)
    real(dp), intent(out) :: extPotGrad(:,:)

    extPotGrad = this%extPotGrad

  end subroutine TExtChargePotGen_getExternalPotGrad


end module extchargepotgen
