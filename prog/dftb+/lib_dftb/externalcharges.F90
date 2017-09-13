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
  implicit none

  private

  public :: init_ExtChrg
  public :: updateCoords_ExtChrg, updateLatVecs_ExtChrg, addShiftPerAtom_ExtChrg
  public :: addEnergyPerAtom_ExtChrg, addForceDCSCC_ExtChrg


  !> Updates the stored coordinates for point charges
  interface updateCoords_ExtChrg
    module procedure updateCoords_ExtChrg_cluster
    module procedure updateCoords_ExtChrg_periodic
  end interface updateCoords_ExtChrg


  !> Add the potential from the point charges
  interface addShiftPerAtom_ExtChrg
    module procedure addShift1
    module procedure addShift2
  end interface addShiftPerAtom_ExtChrg


  !> force contributions from external charges
  interface addForceDCSCC_ExtChrg
    module procedure addForceDCSCC_ExtChrg_cluster
    module procedure addForceDCSCC_ExtChrg_periodic
  end interface addForceDCSCC_ExtChrg

  ! Private module variables


  !> Number of point charges
  integer :: nChrg_

  !> Number of atoms
  integer :: nAtom_

  !> Coordinates of the point charges
  real(dp), allocatable :: coords_(:,:)

  !> Charge of the point charges
  real(dp), allocatable :: charges_(:)

  !> Shift vector
  real(dp), allocatable :: invRVec_(:)

  !> If charges should be blured
  logical :: tBlur_

  !> Blur width for the charges.
  real(dp), allocatable :: blurWidths_(:)

  !> System periodic?
  logical :: tPeriodic_

  !> Calculator initialised?
  logical :: tInitialized_ = .false.

  !> First coordinates received?
  logical :: tUpdated_ = .false.

  !> Real lattice points for Ewald-sum.
  real(dp), allocatable :: rCellVec_(:,:)

contains


  !> Initializes the calculator for external charges
  !>
  !> Note: Blurring of point charges is currently not possible with periodic boundary conditions.
  subroutine init_ExtChrg(coordsAndCharges, nAtom, latVecs, recVecs, ewaldCutoff, blurWidths)

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

    nChrg_ = size(coordsAndCharges, dim=2)

    @:ASSERT(size(coordsAndCharges, dim=1) == 4)
    @:ASSERT(nChrg_ > 0)
    @:ASSERT(present(latVecs) .eqv. present(recVecs))
    @:ASSERT(present(latVecs) .eqv. present(ewaldCutoff))
    @:ASSERT(present(latVecs) .neqv. present(blurWidths))
#:call ASSERT_CODE
    if (present(blurWidths)) then
      @:ASSERT(size(blurWidths) == nChrg_)
    end if
#:endcall ASSERT_CODE

    nAtom_ = nAtom
    allocate(coords_(3, nChrg_))
    coords_ = coordsAndCharges(1:3,:)
    allocate(charges_(nChrg_))
    charges_ = -1.0_dp * coordsAndCharges(4,:)
    allocate(invRVec_(nAtom))
    tPeriodic_ = present(latVecs)
    if (tPeriodic_) then
      !! Fold charges back to unit cell
      call foldCoordToUnitCell(coords_, latVecs, recVecs / (2.0_dp * pi))

      !! Creating the real lattice for the Ewald summation (no neighbor list) The reciprocal part
      !! will be passed from the SCC module, since it is also needed there.
      call getCellTranslations(dummy, rCellVec_, latVecs, recVecs/(2.0_dp*pi), &
          & ewaldCutoff)
    else
      !! Create blurring array for the cluster modell
      if (present(blurWidths)) then
        tBlur_ = any(abs(blurWidths) > 1.0e-7_dp)
      else
        tBlur_ = .false.
      end if
      if (tBlur_) then
        allocate(blurWidths_(nChrg_))
        blurWidths_ = blurWidths
      end if
    end if

    tUpdated_ = .false.
    tInitialized_ = .true.

  end subroutine init_ExtChrg


  !> Updates the module, if the lattice vectors had been changed
  subroutine updateLatVecs_ExtChrg(latVecs, recVecs, ewaldCutoff)

    !> New lattice vectors
    real(dp), intent(in) :: latVecs(:,:)

    !> New reciprocal lattice vectors.
    real(dp), intent(in) :: recVecs(:,:)

    !> Cutoff for the Ewald function.
    real(dp), intent(in) :: ewaldCutoff

    real(dp), allocatable :: dummy(:,:)

    @:ASSERT(tInitialized_ .and. tPeriodic_)

    !! Fold charges back to unit cell
    call foldCoordToUnitCell(coords_, latVecs, recVecs / (2.0_dp * pi))
    call getCellTranslations(dummy, rCellVec_, latVecs, &
        & recVecs / (2.0_dp * pi), ewaldCutoff)

  end subroutine updateLatVecs_ExtChrg


  !> Builds the new shift vectors for new atom coordinates
  subroutine updateCoords_ExtChrg_cluster(atomCoords)

    !> Coordinates of the atoms (not the point charges!)
    real(dp), intent(in) :: atomCoords(:,:)

    @:ASSERT(tInitialized_)
    @:ASSERT(.not. tPeriodic_)
    @:ASSERT(size(atomCoords, dim=1) == 3)
    @:ASSERT(size(atomCoords, dim=2) >= nAtom_)

    if (tBlur_) then
      call sumInvR(invRVec_, nAtom_, nChrg_, atomCoords, coords_, charges_, &
          &blurWidths1=blurWidths_)
    else
      call sumInvR(invRVec_, nAtom_, nChrg_, atomCoords, coords_, charges_)
    end if

    tUpdated_ = .true.

  end subroutine updateCoords_ExtChrg_cluster


  !> Builds the new shift vectors for new atom coordinates
  subroutine updateCoords_ExtChrg_periodic(atomCoords, gLat, alpha, volume)

    !> Coordinates of the atoms (not the point charges!)
    real(dp), intent(in) :: atomCoords(:,:)

    !> Reciprocal lattice vectors
    real(dp), intent(in) :: gLat(:,:)

    !> Ewald parameter
    real(dp), intent(in) :: alpha

    !> Cell volume
    real(dp), intent(in) :: volume

    @:ASSERT(tInitialized_)
    @:ASSERT(tPeriodic_)
    @:ASSERT(size(atomCoords, dim=1) == 3)
    @:ASSERT(size(atomCoords, dim=2) >= nAtom_)
    @:ASSERT(size(gLat, dim=1) == 3)

    call sumInvR(invRVec_, nAtom_, nChrg_, atomCoords, coords_, charges_, &
        &rCellVec_, gLat, alpha, volume)

    tUpdated_ = .true.

  end subroutine updateCoords_ExtChrg_periodic


  !> Adds the contribution of the external charges to the shift vector
  subroutine addShift1(shift)

    !> Shift vector to add the contribution to.
    real(dp), intent(inout) :: shift(:)

    @:ASSERT(tInitialized_ .and. tUpdated_)
    @:ASSERT(size(shift) == nAtom_)

    shift(:) = shift(:) + invRVec_(:)

  end subroutine addShift1


  !> Adds the contribution of the external charges to the shift vector
  subroutine addShift2(shift)

    !> Shift vector to add the contribution to.
    real(dp), intent(inout) :: shift(:,:)

    @:ASSERT(tInitialized_ .and. tUpdated_)
    @:ASSERT(size(shift) == nAtom_)

    shift(:,1) = shift(:,1) + invRVec_(:)

  end subroutine addShift2


  !> Adds the atomic energy contribution do to the external charges.
  subroutine addEnergyPerAtom_ExtChrg(atomCharges, energy)

    !> Charge of the atoms
    real(dp), intent(in) :: atomCharges(:)

    !> Vector containing the energy per atom values.
    real(dp), intent(inout) :: energy(:)

    @:ASSERT(tInitialized_ .and. tUpdated_)
    @:ASSERT(size(atomCharges) == nAtom_)
    @:ASSERT(size(energy) == nAtom_)

    energy(:) = energy(:) + invRVec_(:) * atomCharges(:)

  end subroutine addEnergyPerAtom_ExtChrg


  !> Adds that part of force contribution due to the external charges, which is not contained in the
  !> term with the shift vectors.
  subroutine addForceDCSCC_ExtChrg_cluster(atomForces, chrgForces, atomCoords, atomCharges)

    !> Force vectors on the atoms
    real(dp), intent(inout) :: atomForces(:,:)

    !> Force vectors on the external charges
    real(dp), intent(inout) :: chrgForces(:,:)

    !> Coordinates of the atoms.
    real(dp), intent(in) :: atomCoords(:,:)

    !> Charges of the atoms.
    real(dp), intent(in) :: atomCharges(:)

    @:ASSERT(size(atomForces, dim=1) == 3)
    @:ASSERT(size(atomForces, dim=2) == nAtom_)
    @:ASSERT(size(atomCoords, dim=1) == 3)
    @:ASSERT(size(atomCoords, dim=2) == nAtom_)

    if (tBlur_) then
      call addInvRPrime(atomForces, chrgForces, nAtom_, nChrg_, atomCoords, &
          &coords_, atomCharges, charges_, blurWidths1=blurWidths_)
    else
      call addInvRPrime(atomForces, chrgForces, nAtom_, nChrg_, atomCoords, &
          &coords_, atomCharges, charges_)
    end if

  end subroutine addForceDCSCC_ExtChrg_cluster


  !> Adds that part of force contribution due to the external charges, which is not contained in the
  !> term with the shift vectors.
  subroutine addForceDCSCC_ExtChrg_periodic(atomForces, chrgForces, atomCoords, atomCharges, gVec, &
      & alpha, vol)

    !> Force vectors on the atoms
    real(dp), intent(inout) :: atomForces(:,:)

    !> Force vectors on the external charges
    real(dp), intent(inout) :: chrgForces(:,:)

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
    @:ASSERT(size(atomForces, dim=2) == nAtom_)
    @:ASSERT(size(atomCoords, dim=1) == 3)
    @:ASSERT(size(atomCoords, dim=2) >= nAtom_)

    call addInvRPrime(atomForces, chrgForces, nAtom_, nChrg_, &
        &atomCoords, coords_, atomCharges, charges_, rCellVec_, gVec, alpha, &
        &vol)

  end subroutine addForceDCSCC_ExtChrg_periodic

end module ExternalCharges
