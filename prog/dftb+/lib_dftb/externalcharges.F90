!--------------------------------------------------------------------------------------------------!
!  DFTB+: general package for performing fast atomistic simulations                                !
!  Copyright (C) 2006 - 2020  DFTB+ developers group                                               !
!                                                                                                  !
!  See the LICENSE file for terms of usage and distribution.                                       !
!--------------------------------------------------------------------------------------------------!

#:include 'common.fypp'

!> Routines for calculating the interaction with external charges
module dftbp_externalcharges
  use dftbp_assert
  use dftbp_accuracy
  use dftbp_blasroutines
  use dftbp_coulomb
  use dftbp_constants
  use dftbp_periodic, only : getCellTranslations, foldCoordToUnitCell
  use dftbp_environment
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
    procedure :: setCoordsCluster
    procedure :: setCoordsPeriodic
    procedure :: addForceDcCluster
    procedure :: addForceDcPeriodic
    procedure :: getElStatPotentialCluster
    procedure :: getElStatPotentialPeriodic
    
    !> Updates the stored coordinates for point charges
    generic, public :: setCoordinates => setCoordsCluster, setCoordsPeriodic

    !> Updates the lattice vectors
    procedure, public :: setLatticeVectors

    !> Sets the external charges
    procedure, public :: setExternalCharges

    !> Adds energy contributions per atom
    procedure, public :: addEnergyPerAtom

    !> Adds the potential from the point charges
    procedure, public :: addShiftPerAtom

    !> Adds force double counting component
    generic, public :: addForceDc => addForceDcCluster, addForceDcPeriodic

    !> Returns the electrostatic potential on a grid
    generic, public :: getElStatPotential => getElStatPotentialCluster, getElStatPotentialPeriodic

    !> Copy Q * inverse R contribution for the point charges
    procedure, public :: copyInvRvec
    
  end type TExtCharge


contains


  !> Initializes the calculator for external charges
  subroutine TExtCharge_init(this, nAtom, nExtCharge, tPeriodic)

    !> External charge object
    type(TExtCharge), intent(out) :: this

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
    allocate(this%invRVec(nAtom))
    this%tPeriodic = tPeriodic
    this%tBlur = .false.

    this%tUpdated = .false.
    this%tInitialized = .false.

  end subroutine TExtCharge_init


  !> Sets the external point charges.
  !>
  !> This call should be executed before setLattticeVectors()
  !>
  subroutine setExternalCharges(this, chargeCoords, chargeQs, blurWidths)

    !> Instance
    class(TExtCharge), intent(inout) :: this

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
  subroutine setLatticeVectors(this, latVecs, recVecs)

    !> External charges structure
    class(TExtCharge), intent(inout) :: this

    !> New lattice vectors
    real(dp), intent(in) :: latVecs(:,:)

    !> New reciprocal lattice vectors.
    real(dp), intent(in) :: recVecs(:,:)

    @:ASSERT(this%tInitialized .and. this%tPeriodic)

    !! Fold charges back to unit cell
    call foldCoordToUnitCell(this%coords, latVecs, recVecs / (2.0_dp * pi))

  end subroutine setLatticeVectors


  !> Builds the new shift vectors for new atom coordinates
  subroutine setCoordsCluster(this, env, atomCoords)

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

  end subroutine setCoordsCluster


  !> Builds the new shift vectors for new atom coordinates
  subroutine setCoordsPeriodic(this, env, atomCoords, rCellVec, gLat, alpha, volume)

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

  end subroutine setCoordsPeriodic


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
          & this%charges, atomForces, chrgForces, tHamDeriv=.false., blurWidths1=this%blurWidths)
    else
      call addInvRPrime(env, this%nAtom, this%nChrg, atomCoords, this%coords, atomCharges,&
          & this%charges, atomForces, chrgForces, tHamDeriv=.false.)
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
          & this%charges, rCellVec, gVec, alpha, vol, atomForces, chrgForces, tHamDeriv=.false.,&
          & blurWidths1=this%blurWidths)
    else
      call addInvRPrime(env, this%nAtom, this%nChrg, atomCoords, this%coords, atomCharges,&
          & this%charges, rCellVec, gVec, alpha, vol, atomForces, chrgForces, tHamDeriv=.false.)
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


  !> Copy Q * inverse R contribution for the point charges
  subroutine copyInvRvec(this, qInvR)

    !> Instance of SCC calculation
    class(TExtCharge), intent(in) :: this

    !> (Q * invR) contribution
    real(dp), intent(out) :: qInvR(:)

    qInvR(:) = this%invRVec

  end subroutine copyInvRvec

  
end module dftbp_externalcharges
