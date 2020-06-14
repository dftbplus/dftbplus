!--------------------------------------------------------------------------------------------------!
!  DFTB+: general package for performing fast atomistic simulations                                !
!  Copyright (C) 2006 - 2020  DFTB+ developers group                                               !
!                                                                                                  !
!  See the LICENSE file for terms of usage and distribution.                                       !
!--------------------------------------------------------------------------------------------------!

!> Dummy generator for an external potential given by surrounding point charges
module extchargepotgen
  use dftbp_qdepextpotgen, only : TQDepExtPotGen
  use extchargepot, only : getPointChargePotential
  implicit none
  private

  public :: TExtChargePotGen, TExtChargePotGen_init

  !> double precision
  integer, parameter :: dp = kind(1.0d0)


  !> type extension to add this generator to the internal DFTB structure
  type, extends(TQDepExtPotGen) :: TExtChargePotGen
    private
    !> coordinates of DFTB atoms
    real(dp), allocatable :: atomCoords(:,:)
    !> external charge coordinates
    real(dp), allocatable :: extChargeCoords(:,:)
    !> values of external charges
    real(dp), allocatable :: extCharges(:)
    !> gradients of the external potential wrt to DFTB atom positions
    real(dp), allocatable :: extPotGrad(:,:)
  contains
    !> method for evaluating external potential
    procedure :: getExternalPot => TExtChargePotGen_getExternalPot
    !> method for evaluating gradient of external potential
    procedure :: getExternalPotGrad => TExtChargePotGen_getExternalPotGrad
  end type TExtChargePotGen

contains

  !> Initialise a source of external potential
  subroutine TExtChargePotGen_init(this, atomCoords, extChargeCoords, extCharges)

    !> instance
    type(TExtChargePotGen), intent(out) :: this

    !> coordinates of DFTB atoms
    real(dp), intent(in) :: atomCoords(:,:)

    !> coordinates of external charges
    real(dp), intent(in) :: extChargeCoords(:,:)

    !> values of external charges
    real(dp), intent(in) :: extCharges(:)

    this%atomCoords = atomCoords
    this%extChargeCoords = extChargeCoords
    this%extCharges = extCharges
    allocate(this%extPotGrad(3, size(atomCoords, dim=2)))

  end subroutine TExtChargePotGen_init


  !> evaluate external potential at DFTB atom locations
  subroutine TExtChargePotGen_getExternalPot(this, chargePerAtom, chargePerShell, extPotAtom,&
      & extPotShell)

    !> instance
    class(TExtChargePotGen), intent(inout) :: this

    !> charges for DFTB atoms
    real(dp), intent(in) :: chargePerAtom(:)

    !> charges for each angular shell of the DFTB atoms
    real(dp), intent(in) :: chargePerShell(:,:)

    !> external charges
    real(dp), intent(out) :: extPotAtom(:)

    !> external charges shell charges
    real(dp), intent(out) :: extPotShell(:,:)

    ! call to simple routine to evaluate external potential from this data at the DFTB atom
    ! coordinates
    call getPointChargePotential(this%extChargeCoords, this%extCharges, this%atomCoords,&
        & extPotAtom, this%extPotGrad)
    extPotShell(:,:) = 0.0_dp

  end subroutine TExtChargePotGen_getExternalPot


  !> Generates the derivative of the potential wrt the atom positions
  subroutine TExtChargePotGen_getExternalPotGrad(this, chargePerAtom, chargePerShell, extPotGrad)

    !> instance
    class(TExtChargePotGen), intent(inout) :: this

    !> charges for each DFTB atom
    real(dp), intent(in) :: chargePerAtom(:)

    !> shell resolved atomic charges
    real(dp), intent(in) :: chargePerShell(:,:)

    !> gradient of the potential wrt to atom positions
    real(dp), intent(out) :: extPotGrad(:,:)

    extPotGrad = this%extPotGrad

  end subroutine TExtChargePotGen_getExternalPotGrad


end module extchargepotgen
