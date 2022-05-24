!--------------------------------------------------------------------------------------------------!
!  DFTB+: general package for performing fast atomistic simulations                                !
!  Copyright (C) 2006 - 2022  DFTB+ developers group                                               !
!                                                                                                  !
!  See the LICENSE file for terms of usage and distribution.                                       !
!--------------------------------------------------------------------------------------------------!

#:include 'common.fypp'

!> Routines for calculating the interaction with external charges
module dftbp_dftb_extcharges
  use dftbp_common_accuracy, only : dp
  use dftbp_common_constants, only : pi
  use dftbp_common_environment, only : TEnvironment
  use dftbp_dftb_boundarycond, only : TBoundaryConditions
  use dftbp_dftb_coulomb, only : TCoulomb
  implicit none

  private
  public :: TExtCharges, TExtCharges_init


  !> Private module variables
  type TExtCharges
    private

    !> Number of point charges
    integer :: nChrg

    !> Number of atoms
    integer :: nAtom

    !> Coordinates of the point charges
    real(dp), allocatable :: coords(:,:)

    !> Charge of the point charges
    real(dp), allocatable :: charges(:)

    !> Shift vector
    real(dp), allocatable :: shift(:)

    !> If charges should be blured
    logical :: tBlur

    !> Blur width for the charges.
    real(dp), allocatable :: blurWidths(:)

    !> System periodic?
    logical :: tPeriodic

    !> Calculator initialised?
    logical :: tInitialized = .false.

    !> First coordinates received?
    logical :: tUpdated = .false.

  contains

    !> Sets the coordinates of the charges
    procedure :: setCoordinates

    !> Updates the lattice vectors
    procedure :: setLatticeVectors

    !> Sets the external charges
    procedure :: setExternalCharges

    !> Adds energy contributions per atom
    procedure :: addEnergyPerAtom

    !> Adds the potential from the point charges
    procedure :: addShiftPerAtom

    !> Adds force double counting component
    procedure :: addForceDc

    !> Returns the electrostatic potential on a grid
    procedure :: getElStatPotential

    !> Copy Q * inverse R contribution for the point charges
    procedure :: copyInvRvec

  end type TExtCharges


contains


  !> Initializes the calculator for external charges
  subroutine TExtCharges_init(this, nAtom, nExtCharge, tPeriodic)

    !> External charge object
    type(TExtCharges), intent(out) :: this

    !> Nr. of atoms
    integer, intent(in) :: nAtom

    !> Nr. of external charges
    integer, intent(in) :: nExtCharge

    !> Whether the system is periodic
    logical, intent(in) :: tPeriodic

    this%nAtom = nAtom
    this%nChrg = nExtCharge
    allocate(this%coords(3, this%nChrg))
    allocate(this%charges(this%nChrg))
    allocate(this%shift(nAtom))
    this%tPeriodic = tPeriodic
    this%tBlur = .false.

    this%tUpdated = .false.
    this%tInitialized = .false.

  end subroutine TExtCharges_init


  !> Sets the external point charges.
  !>
  !> This call should be executed before setLattticeVectors()
  !>
  subroutine setExternalCharges(this, chargeCoords, chargeQs, blurWidths)

    !> Instance
    class(TExtCharges), intent(inout) :: this

    !> Coordinates of the external charges as (3, nExtCharges) array
    real(dp), intent(in) :: chargeCoords(:,:)

    !> Charges of the point charges (sign: eletron is negative)
    real(dp), intent(in) :: chargeQs(:)

    !> Width of the Gaussians if the charges are blurred
    real(dp), intent(in), optional :: blurWidths(:)

    this%coords(:,:) = chargeCoords

    ! Adapt to internal sign convention for potentials (electron has positive charge)
    this%charges(:) = -chargeQs

    if (present(blurWidths)) then
      this%tBlur = any(blurWidths > 1.0e-7_dp)
    else
      this%tBlur = .false.
    end if

    if (this%tBlur) then
      if (.not. allocated(this%blurWidths)) then
        allocate(this%blurWidths(this%nChrg))
      end if
      this%blurWidths(:) = blurWidths
    end if

    this%tInitialized = .true.

  end subroutine setExternalCharges


  !> Updates the module, if the lattice vectors had been changed
  subroutine setLatticeVectors(this, latVecs, boundaryConds)

    !> External charges structure
    class(TExtCharges), intent(inout) :: this

    !> New lattice vectors
    real(dp), intent(in) :: latVecs(:,:)

    !> Boundary conditions on the calculation
    type(TBoundaryConditions), intent(in) :: boundaryConds

    @:ASSERT(this%tInitialized .and. this%tPeriodic)

    !! Fold charges back to unit cell
    call boundaryConds%foldCoordsToCell(this%coords, latVecs)

  end subroutine setLatticeVectors


  !> Builds the new shift vectors for new atom coordinates
  subroutine setCoordinates(this, env, coords, coulomb)

    !> External charges structure
    class(TExtCharges), intent(inout) :: this

    !> Environment settings
    type(TEnvironment), intent(in) :: env

    !> Coordinates of the atoms
    real(dp), intent(in) :: coords(:,:)

    !> Coulomb calculator
    type(TCoulomb), intent(inout) :: coulomb

    @:ASSERT(this%tInitialized)

    call coulomb%getPotential(env, coords, this%coords, this%charges, this%shift,&
        & blurWidths=this%blurWidths)
    this%tUpdated = .true.

  end subroutine setCoordinates


  !> Adds the contribution of the external charges to the shift vector
  subroutine addShiftPerAtom(this, shift)

    !> External charges structure
    class(TExtCharges), intent(in) :: this

    !> Shift vector to add the contribution to.
    real(dp), intent(inout) :: shift(:)

    @:ASSERT(this%tInitialized .and. this%tUpdated)
    @:ASSERT(size(shift) == this%nAtom)

    shift(:) = shift + this%shift

  end subroutine addShiftPerAtom


  !> Adds the atomic energy contribution due to the external charges.
  subroutine addEnergyPerAtom(this, atomCharges, energy)

    !> External charges structure
    class(TExtCharges), intent(in) :: this

    !> Charge of the atoms
    real(dp), intent(in) :: atomCharges(:)

    !> Vector containing the energy per atom values.
    real(dp), intent(inout) :: energy(:)

    @:ASSERT(this%tInitialized .and. this%tUpdated)
    @:ASSERT(size(atomCharges) == this%nAtom)
    @:ASSERT(size(energy) == this%nAtom)

    energy(:) = energy + this%shift * atomCharges

  end subroutine addEnergyPerAtom


  !> Adds that part of force contribution due to the external charges, which is not contained in the
  !> term with the shift vectors.
  subroutine addForceDc(this, env, atomForces, chrgForces, atomCoords, atomCharges, coulomb)

    !> External charges structure
    class(TExtCharges), intent(in) :: this

    !> Environment settings
    type(TEnvironment), intent(in) :: env

    !> Force vectors on the atoms
    real(dp), intent(inout) :: atomForces(:,:)

    !> Force vectors on the external charges
    real(dp), intent(inout) :: chrgForces(:,:)

    !> Coordinates of the atoms in the central cell
    real(dp), intent(in) :: atomCoords(:,:)

    !> Charges of the atoms.
    real(dp), intent(in) :: atomCharges(:)

    !> Coulomb calculator
    type(TCoulomb), intent(inout) :: coulomb

    @:ASSERT(size(atomForces, dim=1) == 3)
    @:ASSERT(size(atomForces, dim=2) == this%nAtom)
    @:ASSERT(size(atomCoords, dim=1) == 3)
    @:ASSERT(size(atomCoords, dim=2) == this%nAtom)

    call coulomb%addExternalPotGrad(env, atomCoords, this%coords, atomCharges, this%charges,&
        & atomForces, chrgForces, tHamDeriv=.false., extChargeBlurWidths=this%blurWidths)

  end subroutine addForceDc


  !> Returns potential from external charges (periodic case)
  subroutine getElStatPotential(this, env, locations, V, coulomb, epsSoften)

    !> Instance of SCC calculation
    class(TExtCharges), intent(in) :: this

    !> Computational environment settings
    type(TEnvironment), intent(in) :: env

    !> sites to calculate potential
    real(dp), intent(in) :: locations(:,:)

    !> Resulting potentials
    real(dp), intent(out) :: V(:)

    !> Coulomb calculator
    type(TCoulomb), intent(inout) :: coulomb

    !> optional potential softening
    real(dp), optional, intent(in) :: epsSoften

    @:ASSERT(all(shape(locations) == [3, size(V)]))

    call coulomb%getPotential(env, locations, this%coords, -this%charges, V,&
        & blurWidths=this%blurWidths, epsSoften=epsSoften)

  end subroutine getElStatPotential


  !> Copy Q * inverse R contribution for the point charges
  subroutine copyInvRvec(this, qInvR)

    !> Instance of SCC calculation
    class(TExtCharges), intent(in) :: this

    !> (Q * invR) contribution
    real(dp), intent(out) :: qInvR(:)

    qInvR(:) = this%shift

  end subroutine copyInvRvec


end module dftbp_dftb_extcharges
