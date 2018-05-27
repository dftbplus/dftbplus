!--------------------------------------------------------------------------------------------------!
!  DFTB+: general package for performing fast atomistic simulations                                !
!  Copyright (C) 2018  DFTB+ developers group                                                      !
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

  contains
    private

    procedure :: updateCoordsCluster
    procedure :: updateCoordsPeriodic
    procedure :: addForceDcCluster
    procedure :: addForceDcPeriodic
    procedure :: getElStatPotentialCluster
    procedure :: getElStatPotentialPeriodic
    
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

    !> Returns the electrostatic potential on a grid
    generic, public :: getElStatPotential => getElStatPotentialCluster, getElStatPotentialPeriodic
    
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
    end if

    !! Create blurring array
    if (present(blurWidths)) then
      this%tBlur = any(blurWidths > 1.0e-7_dp)
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
  subroutine updateLatVecs(this, latVecs, recVecs)

    !> External charges structure
    class(TExtCharge), intent(inout) :: this

    !> New lattice vectors
    real(dp), intent(in) :: latVecs(:,:)

    !> New reciprocal lattice vectors.
    real(dp), intent(in) :: recVecs(:,:)

    @:ASSERT(this%tInitialized .and. this%tPeriodic)

    !! Fold charges back to unit cell
    call foldCoordToUnitCell(this%coords, latVecs, recVecs / (2.0_dp * pi))

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
      call sumInvR(env, this%nAtom, this%nChrg, atomCoords, this%coords, this%charges,&
          & this%invRVec, blurWidths1=this%blurWidths)
    else
      call sumInvR(env, this%nAtom, this%nChrg, atomCoords, this%coords, this%charges, this%invRVec)
    end if

    this%tUpdated = .true.

  end subroutine updateCoordsCluster


  !> Builds the new shift vectors for new atom coordinates
  subroutine updateCoordsPeriodic(this, env, atomCoords, rCellVec, gLat, alpha, volume)

    !> External charges structure
    class(TExtCharge), intent(inout) :: this

    !> Environment settings
    type(TEnvironment), intent(in) :: env

    !> Coordinates of the atoms (not the point charges!)
    real(dp), intent(in) :: atomCoords(:,:)

    !> Real lattice points for asymmetric Ewald sum
    real(dp), intent(in) :: rCellVec(:,:)

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
      call sumInvR(env, this%nAtom, this%nChrg, atomCoords, this%coords, this%charges,&
          & rCellVec, gLat, alpha, volume, this%invRVec, blurWidths1=this%blurWidths)
    else
      call sumInvR(env, this%nAtom, this%nChrg, atomCoords, this%coords, this%charges,&
          & rCellVec, gLat, alpha, volume, this%invRVec)
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
  subroutine addForceDcCluster(this, env, atomForces, chrgForces, atomCoords, atomCharges)

    !> External charges structure
    class(TExtCharge), intent(in) :: this

    !> Environment settings
    type(TEnvironment), intent(in) :: env

    !> Force vectors on the atoms
    real(dp), intent(inout) :: atomForces(:,:)

    !> Force vectors on the external charges
    real(dp), intent(inout) :: chrgForces(:,:)

    !> Coordinates of the atoms.
    real(dp), intent(in) :: atomCoords(:,:)

    !> Charges of the atoms.
    real(dp), intent(in) :: atomCharges(:)

    @:ASSERT(size(atomForces, dim=1) == 3)
    @:ASSERT(size(atomForces, dim=2) == this%nAtom)
    @:ASSERT(size(atomCoords, dim=1) == 3)
    @:ASSERT(size(atomCoords, dim=2) == this%nAtom)

    if (this%tBlur) then
      call addInvRPrime(env, this%nAtom, this%nChrg, atomCoords, this%coords, atomCharges,&
          & this%charges, atomForces, chrgForces, blurWidths1=this%blurWidths)
    else
      call addInvRPrime(env, this%nAtom, this%nChrg, atomCoords, this%coords, atomCharges,&
          & this%charges, atomForces, chrgForces)
    end if

  end subroutine addForceDcCluster


  !> Adds that part of force contribution due to the external charges, which is not contained in the
  !> term with the shift vectors.
  subroutine addForceDcPeriodic(this, env, atomForces, chrgForces, atomCoords, atomCharges,&
      & rCellVec, gVec, alpha, vol)

    !> External charges structure
    class(TExtCharge), intent(in) :: this

    !> Environment settings
    type(TEnvironment), intent(in) :: env

    !> Force vectors on the atoms
    real(dp), intent(inout) :: atomForces(:,:)

    !> Force vectors on the external charges
    real(dp), intent(inout) :: chrgForces(:,:)

    !> Coordinates of the atoms.
    real(dp), intent(in) :: atomCoords(:,:)

    !> Charges of the atoms.
    real(dp), intent(in) :: atomCharges(:)

    !> Real lattice points for asymmetric Ewald sum
    real(dp), intent(in) :: rCellVec(:,:)

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
      call addInvRPrime(env, this%nAtom, this%nChrg, atomCoords, this%coords, atomCharges,&
          & this%charges, rCellVec, gVec, alpha, vol, atomForces, chrgForces,&
          & blurWidths1=this%blurWidths)
    else
      call addInvRPrime(env, this%nAtom, this%nChrg, atomCoords, this%coords, atomCharges,&
          & this%charges, rCellVec, gVec, alpha, vol, atomForces, chrgForces)
    end if

  end subroutine addForceDcPeriodic


  !> Returns potential from external charges (periodic case)
  subroutine getElStatPotentialPeriodic(this, env, locations, rCellVec, gLatPoint, alpha, volume,&
      & V, epsSoften)

    !> Instance of SCC calculation
    class(TExtCharge), intent(in) :: this

    !> Computational environment settings
    type(TEnvironment), intent(in) :: env

    !> sites to calculate potential
    real(dp), intent(in) :: locations(:,:)

    !> Real lattice points for Ewald-sum.
    real(dp), intent(in) :: rCellVec(:,:)

    !> Lattice points for reciprocal Ewald
    real(dp), intent(in) :: gLatPoint(:,:)

    !> Parameter for Ewald
    real(dp), intent(in) :: alpha

    !> Cell volume
    real(dp), intent(in) :: volume

    !> Resulting potentials
    real(dp), intent(out) :: V(:)

    !> optional potential softening
    real(dp), optional, intent(in) :: epsSoften

    @:ASSERT(all(shape(locations) == [3, size(V)]))

    V(:) = 0.0_dp

    if (allocated(this%blurWidths)) then
      call sumInvR(env, size(V), size(this%charges), locations, this%coords, -this%charges,&
          & rCellVec, gLatPoint, alpha, volume, V, this%blurWidths, epsSoften=epsSoften)
    else
      call sumInvR(env, size(V), size(this%charges), locations, this%coords, -this%charges,&
          & rCellVec, gLatPoint, alpha, volume, V, epsSoften=epsSoften)
    end if

  end subroutine getElStatPotentialPeriodic


  !> Returns potential from external charges (cluster case)
  subroutine getElStatPotentialCluster(this, env, locations, V, epsSoften)

    !> Instance of SCC calculation
    class(TExtCharge), intent(in) :: this

    !> Computational environment settings
    type(TEnvironment), intent(in) :: env

    !> sites to calculate potential
    real(dp), intent(in) :: locations(:,:)

    !> Resulting potentials
    real(dp), intent(out) :: V(:)

    !> optional potential softening
    real(dp), optional, intent(in) :: epsSoften

    @:ASSERT(all(shape(locations) == [3, size(V)]))

    V(:) = 0.0_dp

    if (allocated(this%blurWidths)) then
      call sumInvR(env, size(V), size(this%charges), locations, this%coords, -this%charges, V,&
          & this%blurWidths, epsSoften=epsSoften)
    else
      call sumInvR(env, size(V), size(this%charges), locations, this%coords, -this%charges, V,&
          & epsSoften=epsSoften)
    end if

  end subroutine getElStatPotentialCluster


end module ExternalCharges
