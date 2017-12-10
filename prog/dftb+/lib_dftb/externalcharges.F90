!--------------------------------------------------------------------------------------------------!
!  DFTB+: general package for performing fast atomistic simulations                                !
!  Copyright (C) 2017  DFTB+ developers group                                                      !
!                                                                                                  !
!  See the LICENSE file for terms of usage and distribution.                                       !
!--------------------------------------------------------------------------------------------------!

#:include 'common.fypp'

!> Routines for calculating the interaction with external charges
module ExternalCharges
  use assert
  use accuracy
  use blasroutines
  use coulomb
  use constants
  use periodic, only : getCellTranslations, foldCoordToUnitCell
  use environment
  implicit none

  private

  public :: TExtCharge, TExtCharge_init


  !> Private module variables
  type TExtCharge
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
    real(dp), allocatable :: invRVec(:)

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

    !> Real lattice points for Ewald-sum.
    real(dp), allocatable :: rCellVec(:,:)

  contains
    private

    procedure :: updateCoordsCluster
    procedure :: updateCoordsPeriodic
    procedure :: addForceDcCluster
    procedure :: addForceDcPeriodic

    !> Updates the stored coordinates for point charges
    generic, public :: updateCoords => updateCoordsCluster, updateCoordsPeriodic

    !> Updates the lattice vectors
    procedure, public :: updateLatVecs

    !> Adds energy contributions per atom
    procedure, public :: addEnergyPerAtom

    !> Adds the potential from the point charges
    procedure, public :: addShiftPerAtom

    !> Adds force double counting component
    generic, public :: addForceDc => addForceDcCluster, addForceDcPeriodic

  end type TExtCharge


contains


  !> Initializes the calculator for external charges
  subroutine TExtCharge_init(this, coordsAndCharges, nAtom, latVecs, recVecs, ewaldCutoff,&
      & blurWidths)

    !> External charge object
    type(TExtCharge), intent(out) :: this

    !> (4, nAtom) array with coordinates and charges
    real(dp), intent(in) :: coordsAndCharges(:,:)

    !> Number of atoms
    integer, intent(in) :: nAtom

    !> Blurring for the external charges
    real(dp), intent(in), optional :: latVecs(:,:)

    !> Lattice vectors of the supercell
    real(dp), intent(in), optional :: recVecs(:,:)

    !> Reciprocal vectors
    real(dp), intent(in), optional :: ewaldCutoff

    !> Cutoff for the ewald summation
    real(dp), intent(in), optional :: blurWidths(:)

    real(dp), allocatable :: dummy(:,:)

    this%nChrg = size(coordsAndCharges, dim=2)

    @:ASSERT(size(coordsAndCharges, dim=1) == 4)
    @:ASSERT(this%nChrg > 0)
    @:ASSERT(present(latVecs) .eqv. present(recVecs))
    @:ASSERT(present(latVecs) .eqv. present(ewaldCutoff))
#:call ASSERT_CODE
    if (present(blurWidths)) then
      @:ASSERT(size(blurWidths) == this%nChrg)
    end if
#:endcall ASSERT_CODE

    this%nAtom = nAtom
    allocate(this%coords(3, this%nChrg))
    this%coords = coordsAndCharges(1:3,:)
    allocate(this%charges(this%nChrg))
    this%charges = -1.0_dp * coordsAndCharges(4,:)
    allocate(this%invRVec(nAtom))
    this%tPeriodic = present(latVecs)
    if (this%tPeriodic) then
      !! Fold charges back to unit cell
      call foldCoordToUnitCell(this%coords, latVecs, recVecs / (2.0_dp * pi))

      !! Creating the real lattice for the Ewald summation (no neighbor list) The reciprocal part
      !! will be passed from the SCC module, since it is also needed there.
      call getCellTranslations(dummy, this%rCellVec, latVecs, recVecs/(2.0_dp*pi), ewaldCutoff)
    end if

    !! Create blurring array
    if (present(blurWidths)) then
      this%tBlur = any(abs(blurWidths) > 1.0e-7_dp)
    else
      this%tBlur = .false.
    end if
    if (this%tBlur) then
      allocate(this%blurWidths(this%nChrg))
      this%blurWidths = blurWidths
    end if

    this%tUpdated = .false.
    this%tInitialized = .true.

  end subroutine TExtCharge_init


  !> Updates the module, if the lattice vectors had been changed
  subroutine updateLatVecs(this, latVecs, recVecs, ewaldCutoff)

    !> External charges structure
    class(TExtCharge), intent(inout) :: this

    !> New lattice vectors
    real(dp), intent(in) :: latVecs(:,:)

    !> New reciprocal lattice vectors.
    real(dp), intent(in) :: recVecs(:,:)

    !> Cutoff for the Ewald function.
    real(dp), intent(in) :: ewaldCutoff

    real(dp), allocatable :: dummy(:,:)

    @:ASSERT(this%tInitialized .and. this%tPeriodic)

    !! Fold charges back to unit cell
    call foldCoordToUnitCell(this%coords, latVecs, recVecs / (2.0_dp * pi))
    call getCellTranslations(dummy, this%rCellVec, latVecs, recVecs / (2.0_dp * pi), ewaldCutoff)

  end subroutine updateLatVecs


  !> Builds the new shift vectors for new atom coordinates
  subroutine updateCoordsCluster(this, env, atomCoords)

    !> External charges structure
    class(TExtCharge), intent(inout) :: this

    !> Environment settings
    type(TEnvironment), intent(in) :: env

    !> Coordinates of the atoms (not the point charges!)
    real(dp), intent(in) :: atomCoords(:,:)

    @:ASSERT(this%tInitialized)
    @:ASSERT(.not. this%tPeriodic)
    @:ASSERT(size(atomCoords, dim=1) == 3)
    @:ASSERT(size(atomCoords, dim=2) >= this%nAtom)

    if (this%tBlur) then
      call sumInvR(this%invRVec, env, this%nAtom, this%nChrg, atomCoords, this%coords,&
          & this%charges, blurWidths1=this%blurWidths)
    else
      call sumInvR(this%invRVec, env, this%nAtom, this%nChrg, atomCoords, this%coords, this%charges)
    end if

    this%tUpdated = .true.

  end subroutine updateCoordsCluster


  !> Builds the new shift vectors for new atom coordinates
  subroutine updateCoordsPeriodic(this, env, atomCoords, gLat, alpha, volume)

    !> External charges structure
    class(TExtCharge), intent(inout) :: this

    !> Environment settings
    type(TEnvironment), intent(in) :: env

    !> Coordinates of the atoms (not the point charges!)
    real(dp), intent(in) :: atomCoords(:,:)

    !> Reciprocal lattice vectors
    real(dp), intent(in) :: gLat(:,:)

    !> Ewald parameter
    real(dp), intent(in) :: alpha

    !> Cell volume
    real(dp), intent(in) :: volume

    @:ASSERT(this%tInitialized)
    @:ASSERT(this%tPeriodic)
    @:ASSERT(size(atomCoords, dim=1) == 3)
    @:ASSERT(size(atomCoords, dim=2) >= this%nAtom)
    @:ASSERT(size(gLat, dim=1) == 3)

    if (this%tBlur) then
      call sumInvR(this%invRVec, env, this%nAtom, this%nChrg, atomCoords, this%coords,&
          & this%charges, this%rCellVec, gLat, alpha, volume, blurWidths1=this%blurWidths)
    else
      call sumInvR(this%invRVec, env, this%nAtom, this%nChrg, atomCoords, this%coords,&
          & this%charges, this%rCellVec, gLat, alpha, volume)
    end if

    this%tUpdated = .true.

  end subroutine updateCoordsPeriodic


  !> Adds the contribution of the external charges to the shift vector
  subroutine addShiftPerAtom(this, shift)

    !> External charges structure
    class(TExtCharge), intent(in) :: this

    !> Shift vector to add the contribution to.
    real(dp), intent(inout) :: shift(:)

    @:ASSERT(this%tInitialized .and. this%tUpdated)
    @:ASSERT(size(shift) == this%nAtom)

    shift(:) = shift + this%invRVec

  end subroutine addShiftPerAtom


  !> Adds the atomic energy contribution due to the external charges.
  subroutine addEnergyPerAtom(this, atomCharges, energy)

    !> External charges structure
    class(TExtCharge), intent(in) :: this

    !> Charge of the atoms
    real(dp), intent(in) :: atomCharges(:)

    !> Vector containing the energy per atom values.
    real(dp), intent(inout) :: energy(:)

    @:ASSERT(this%tInitialized .and. this%tUpdated)
    @:ASSERT(size(atomCharges) == this%nAtom)
    @:ASSERT(size(energy) == this%nAtom)

    energy(:) = energy(:) + this%invRVec(:) * atomCharges(:)

  end subroutine addEnergyPerAtom


  !> Adds that part of force contribution due to the external charges, which is not contained in the
  !> term with the shift vectors.
  subroutine addForceDcCluster(this, atomForces, chrgForces, env, atomCoords, atomCharges)

    !> External charges structure
    class(TExtCharge), intent(in) :: this

    !> Force vectors on the atoms
    real(dp), intent(inout) :: atomForces(:,:)

    !> Force vectors on the external charges
    real(dp), intent(inout) :: chrgForces(:,:)

    !> Environment settings
    type(TEnvironment), intent(in) :: env

    !> Coordinates of the atoms.
    real(dp), intent(in) :: atomCoords(:,:)

    !> Charges of the atoms.
    real(dp), intent(in) :: atomCharges(:)

    @:ASSERT(size(atomForces, dim=1) == 3)
    @:ASSERT(size(atomForces, dim=2) == this%nAtom)
    @:ASSERT(size(atomCoords, dim=1) == 3)
    @:ASSERT(size(atomCoords, dim=2) == this%nAtom)

    if (this%tBlur) then
      call addInvRPrime(atomForces, chrgForces, env, this%nAtom, this%nChrg, atomCoords,&
          & this%coords, atomCharges, this%charges, blurWidths1=this%blurWidths)
    else
      call addInvRPrime(atomForces, chrgForces, env, this%nAtom, this%nChrg, atomCoords,&
          & this%coords, atomCharges, this%charges)
    end if

  end subroutine addForceDcCluster


  !> Adds that part of force contribution due to the external charges, which is not contained in the
  !> term with the shift vectors.
  subroutine addForceDcPeriodic(this, atomForces, chrgForces, env, atomCoords, atomCharges, gVec,&
      & alpha, vol)

    !> External charges structure
    class(TExtCharge), intent(in) :: this

    !> Force vectors on the atoms
    real(dp), intent(inout) :: atomForces(:,:)

    !> Force vectors on the external charges
    real(dp), intent(inout) :: chrgForces(:,:)

    !> Environment settings
    type(TEnvironment), intent(in) :: env

    !> Coordinates of the atoms.
    real(dp), intent(in) :: atomCoords(:,:)

    !> Charges of the atoms.
    real(dp), intent(in) :: atomCharges(:)

    !> Reciprocal lattice vectors for the Ewald summation
    real(dp), intent(in) :: gVec(:,:)

    !> Parameter of the Ewald summation
    real(dp), intent(in) :: alpha

    !> Cell volume
    real(dp), intent(in) :: vol

    @:ASSERT(size(atomForces, dim=1) == 3)
    @:ASSERT(size(atomForces, dim=2) == this%nAtom)
    @:ASSERT(size(atomCoords, dim=1) == 3)
    @:ASSERT(size(atomCoords, dim=2) >= this%nAtom)

    if (this%tBlur) then
      call addInvRPrime(atomForces, chrgForces, env, this%nAtom, this%nChrg, atomCoords,&
          & this%coords, atomCharges, this%charges, this%rCellVec, gVec, alpha, vol,&
          & blurWidths1=this%blurWidths)
    else
      call addInvRPrime(atomForces, chrgForces, env, this%nAtom, this%nChrg, atomCoords,&
          & this%coords, atomCharges, this%charges, this%rCellVec, gVec, alpha, vol)
    end if

  end subroutine addForceDcPeriodic

end module ExternalCharges
