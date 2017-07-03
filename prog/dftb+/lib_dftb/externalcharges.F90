!--------------------------------------------------------------------------------------------------!
!  DFTB+: general package for performing fast atomistic simulations                                !
!  Copyright (C) 2017  DFTB+ developers group                                                      !
!                                                                                                  !
!  See the LICENSE file for terms of usage and distribution.                                       !
!--------------------------------------------------------------------------------------------------!

#:include 'common.fypp'

!!* Routines for calculating the interaction with external charges
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

  interface updateCoords_ExtChrg
    module procedure updateCoords_ExtChrg_cluster
    module procedure updateCoords_ExtChrg_periodic
  end interface

  interface addShiftPerAtom_ExtChrg
    module procedure addShift1
    module procedure addShift2
  end interface

  interface addForceDCSCC_ExtChrg
    module procedure addForceDCSCC_ExtChrg_cluster
    module procedure addForceDCSCC_ExtChrg_periodic
  end interface

  !! Private module variables
  integer :: nChrg_        ! Number of point charges
  integer :: nAtom_        ! Number of atoms
  real(dp), allocatable :: coords_(:,:)  ! Coordinates of the point charges
  real(dp), allocatable :: charges_(:)   ! Charge of the point charges
  real(dp), allocatable :: invRVec_(:)   ! Shift vector
  logical :: tBlur_                      ! If charges should be blured
  real(dp), allocatable :: blurWidths_(:)! Blur width for the charges.
  logical :: tPeriodic_                  ! System periodic?
  logical :: tInitialized_ = .false.     ! Calculator initialised?
  logical :: tUpdated_ = .false.         ! First coordinates received?
  real(dp), allocatable :: rCellVec_(:,:)  ! Real lattice points for Ewald-sum.



contains

  !!* Initializes the calculator for external charges
  !!* @param coordsAndCharges (4, nAtom) array with coordinates and charges
  !!* @param nAtom Number of atoms
  !!* @param blurWidths Blurring for the external charges
  !!* @param latVecs Lattice vectors of the supercell
  !!* @param recVecs Reciprocal vectors
  !!* @param ewaldCutoff Cutoff for the ewald summation
  !!* @note Blurring of point charges is currently not possible with periodic
  !!*  boundary conditions.
  subroutine init_ExtChrg(coordsAndCharges, nAtom, latVecs, recVecs, &
      &ewaldCutoff, blurWidths)
    real(dp), intent(in) :: coordsAndCharges(:,:)
    integer, intent(in) :: nAtom
    real(dp), intent(in), optional :: latVecs(:,:)
    real(dp), intent(in), optional :: recVecs(:,:)
    real(dp), intent(in), optional :: ewaldCutoff
    real(dp), intent(in), optional :: blurWidths(:)

    real(dp), allocatable :: dummy(:,:)

    nChrg_ = size(coordsAndCharges, dim=2)

    @:ASSERT(size(coordsAndCharges, dim=1) == 4)
    @:ASSERT(nChrg_ > 0)
    @:ASSERT(present(latVecs) .eqv. present(recVecs))
    print *, present(latVecs), present(ewaldCutoff)
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

      !! Creating the real lattice for the Ewald summation (no neighbor list)
      !! The reciprocal part will be passed from the SCC module, since it is
      !! also needed there.
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


  !!* Updates the module, if the lattice vectors had been changed
  !!* @param latVecs  New lattice vectors
  !!* @param recVecs  New reciprocal lattice vectors.
  !!* @param ewaldCutoff  Cutoff for the Ewald function.
  subroutine updateLatVecs_ExtChrg(latVecs, recVecs, ewaldCutoff)
    real(dp), intent(in) :: latVecs(:,:), recVecs(:,:)
    real(dp), intent(in) :: ewaldCutoff

    real(dp), allocatable :: dummy(:,:)

    @:ASSERT(tInitialized_ .and. tPeriodic_)

    !! Fold charges back to unit cell
    call foldCoordToUnitCell(coords_, latVecs, recVecs / (2.0_dp * pi))
    call getCellTranslations(dummy, rCellVec_, latVecs, &
        & recVecs / (2.0_dp * pi), ewaldCutoff)

  end subroutine updateLatVecs_ExtChrg


  !!* Builds the new shift vectors for new atom coordinates
  !!* @param atomCoords Coordinates of the atoms (not the point charges!)
  subroutine updateCoords_ExtChrg_cluster(atomCoords)
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



  !!* Builds the new shift vectors for new atom coordinates
  !!* @param atomCoords Coordinates of the atoms (not the point charges!)
  subroutine updateCoords_ExtChrg_periodic(atomCoords, gLat, alpha, &
      &volume)
    real(dp), intent(in) :: atomCoords(:,:)
    real(dp), intent(in) :: gLat(:,:)
    real(dp), intent(in) :: alpha
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



  !!* Adds the contribution of the external charges to the shift vector
  !!* @param shift Shift vector to add the contribution to.
  subroutine addShift1(shift)
    real(dp), intent(inout) :: shift(:)

    @:ASSERT(tInitialized_ .and. tUpdated_)
    @:ASSERT(size(shift) == nAtom_)

    shift(:) = shift(:) + invRVec_(:)

  end subroutine addShift1

  !!* Adds the contribution of the external charges to the shift vector
  !!* @param shift Shift vector to add the contribution to.
  subroutine addShift2(shift)
    real(dp), intent(inout) :: shift(:,:)

    @:ASSERT(tInitialized_ .and. tUpdated_)
    @:ASSERT(size(shift) == nAtom_)

    shift(:,1) = shift(:,1) + invRVec_(:)

  end subroutine addShift2



  !!* Adds the atomic energy contribution do to the external charges.
  !!* @param atomCharges Charge of the atoms
  !!* @param energy Vector containing the energy per atom values.
  subroutine addEnergyPerAtom_ExtChrg(atomCharges, energy)
    real(dp), intent(in) :: atomCharges(:)
    real(dp), intent(inout) :: energy(:)

    @:ASSERT(tInitialized_ .and. tUpdated_)
    @:ASSERT(size(atomCharges) == nAtom_)
    @:ASSERT(size(energy) == nAtom_)

    energy(:) = energy(:) + invRVec_(:) * atomCharges(:)

  end subroutine addEnergyPerAtom_ExtChrg



  !!* Adds that part of force contribution due to the external charges, which is
  !!* not contained in the term with the shift vectors.
  !!* @param atomForces Force vectors on the atoms
  !!* @param chrgForces Force vectors on the external charges
  !!* @param atomCoords Coordinates of the atoms.
  !!* @param atomCharges Charges of the atoms.
  subroutine addForceDCSCC_ExtChrg_cluster(atomForces, chrgForces, &
      &atomCoords, atomCharges)
    real(dp), intent(inout) :: atomForces(:,:)
    real(dp), intent(inout) :: chrgForces(:,:)
    real(dp), intent(in) :: atomCoords(:,:)
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



  !!* Adds that part of force contribution due to the external charges, which is
  !!* not contained in the term with the shift vectors.
  !!* @param atomForces Force vectors on the atoms
  !!* @param chrgForces Force vectors on the external charges
  !!* @param atomCoords Coordinates of the atoms.
  !!* @param atomCharges Charges of the atoms.
  !!* @param gVec Reciprocal lattice vectors for the Ewald summation
  !!* @param alpha Parameter of the Ewald summation
  !!* @param vol Cell volume
  subroutine addForceDCSCC_ExtChrg_periodic(atomForces, chrgForces, &
      &atomCoords, atomCharges, gVec, alpha, vol)
    real(dp), intent(inout) :: atomForces(:,:)
    real(dp), intent(inout) :: chrgForces(:,:)
    real(dp), intent(in) :: atomCoords(:,:)
    real(dp), intent(in) :: atomCharges(:)
    real(dp), intent(in) :: gVec(:,:)
    real(dp), intent(in) :: alpha
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
