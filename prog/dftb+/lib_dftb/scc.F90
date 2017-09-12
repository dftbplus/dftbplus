!--------------------------------------------------------------------------------------------------!
!  DFTB+: general package for performing fast atomistic simulations                                !
!  Copyright (C) 2017  DFTB+ developers group                                                      !
!                                                                                                  !
!  See the LICENSE file for terms of usage and distribution.                                       !
!--------------------------------------------------------------------------------------------------!

#:include 'common.fypp'

!> Functions and local variables for the SCC calculation.
module scc
  use assert
  use accuracy
  use message
  use coulomb
  use shortgamma
  use fileid
  use constants
  use periodic
  use externalcharges
  use blasroutines
  use commontypes
  use chargeconstr
  use shift
  implicit none

  private

  public :: TSCCInit, typSCC, init_SCC
  public :: SCC_getOrbitalEquiv, updateCoords_SCC, updateLatVecs_SCC, updateCharges_SCC
  public :: getSCCCutoff, getEnergyPerAtom_SCC, getEnergyPerAtom_SCC_Xlbomd
  public :: addForceDCSCC, addForceDCSCC_Xlbomd, addStressDCSCC
  public ::  getShiftPerAtom, getShiftPerL, getSCCEwaldPar
  public :: getAtomicGammaMatrix

  !> tolerance for Ewald - should be changed to be possible to override it
  real(dp), parameter :: tolEwald = 1.0e-9_dp


  !> Data necessary to initialize the SCC module
  type TSCCInit

    !> atomic orbital information
    type(TOrbitals), pointer :: orb

    !> hubbard U values
    real(dp), allocatable :: hubbU(:,:)

    !> use damping of short range part
    logical, allocatable :: tDampedShort(:)

    !> exponent if damping used
    real(dp) :: dampExp = 0.0_dp

    !> lattice vectors
    real(dp), allocatable :: latVecs(:,:)

    !> reciprocal lattice vectors
    real(dp), allocatable :: recVecs(:,:)

    !> cell volume
    real(dp) :: volume = 0.0_dp

    !> external charges
    real(dp), allocatable :: extCharges(:,:)

    !> if broadened out
    real(dp), allocatable :: blurWidths(:)

    !> any constraints on atomic charges
    real(dp), allocatable :: chrgConstraints(:,:)

    !> third order energy contributions
    real(dp), allocatable :: thirdOrderOn(:,:)

    !> if > 0 -> manual setting for alpha
    real(dp) :: ewaldAlpha = 0.0_dp
  end type TSCCInit

  !> private module variables for SCC
  type typSCC
    private

    !> If module is initialised
    logical :: tInitialised_ = .false.

    !> If coordinates updated at least once
    logical :: tCoordUp_

    !> If charges updated at least once
    logical :: tChrgUp_

    !> Nr. of atoms
    integer :: nAtom_

    !> Nr. of species
    integer :: nSpecies_

    !> Max nr. of orbitals per atom
    integer :: mOrb_

    !> Max nr. of shells per atom
    integer :: mShell_

    !> Maximal number of unique U values
    integer :: mHubbU_

    !> Stores 1/r between atom pairs
    real(dp), allocatable :: invRMat_(:,:)

    !> Shift vector per atom
    real(dp), allocatable :: shiftPerAtom_(:)

    !> Shift vector per l-shell
    real(dp), allocatable :: shiftPerL_(:,:)

    !> Short range interaction
    real(dp), allocatable :: shortGamma_(:,:,:,:)

    !> Cutoff for short range int.
    real(dp), allocatable :: shortCutoff_(:,:,:,:)

    !> Maximal cutoff
    real(dp) :: cutoff_

    !> Lattice points for reciprocal Ewald
    real(dp), allocatable :: gLatPoint_(:,:)

    !> Nr. of neighbors for short range interaction
    integer, allocatable :: nNeighShort_(:,:,:,:)

    !> Nr. of neigh for real Ewald
    integer, allocatable :: nNeighEwald_(:)

    !> Cutoff for real Ewald
    real(dp) :: maxREwald_

    !> Parameter for Ewald
    real(dp) :: alpha_

    !> Cell volume
    real(dp) :: volume_

    !> Uniq Us per species
    real(dp), allocatable :: uniqHubbU_(:,:)

    !> Nr. of uniq Us per species
    integer, allocatable :: nHubbU_(:)

    !> Mapping L-shell -> uniq U
    integer, allocatable :: iHubbU_(:,:)

    !> Are external charges present?
    logical :: tExtChrg_

    !> Negative net charge
    real(dp), allocatable :: deltaQ_(:,:)

    !> Negative net charge per shell
    real(dp), allocatable :: deltaQPerLShell_(:,:)

    !> Negative net charge per atom
    real(dp), allocatable :: deltaQAtom_(:)

    !> Negative net charge per U
    real(dp), allocatable :: deltaQUniqU_(:,:)

    !> Damped short range? (nSpecies)
    logical, allocatable :: tDampedShort_(:)

    !> Damping exponent
    real(dp) :: dampExp_

    !> Is the system periodic?
    logical :: tPeriodic_

    !> Shifts due charge constrains?
    logical :: tChrgConstr_

    !> Object for charge constraints
    type(OChrgConstr), allocatable :: chrgConstr_

    !> use third order contributions
    logical :: tThirdOrder_

    !> Shifts due to 3rd order
    type(OChrgConstr), allocatable :: thirdOrder_

    !> evaluate Ewald parameter
    logical :: tAutoEwald_

  end type typSCC


contains


  !> Initialises the SCC module
  subroutine init_SCC(inp, OSCC)

    !> SCC ADT
    type(TSCCInit), intent(inout) :: inp

    !> Resulting module variables
    type(typSCC), intent(inout) :: OSCC

    integer :: iSp1, iSp2, iU1, iU2, iL
    real(dp) :: maxGEwald

    @:ASSERT(.not. OSCC%tInitialised_)

    OSCC%nSpecies_ = size(inp%orb%nOrbSpecies)
    OSCC%nAtom_ = size(inp%orb%nOrbAtom)
    OSCC%mShell_ = inp%orb%mShell
    OSCC%mOrb_ = inp%orb%mOrb

    @:ASSERT(allocated(inp%latVecs) .eqv. allocated(inp%recVecs))
    @:ASSERT(allocated(inp%latVecs) .eqv. (inp%volume > 0.0_dp))
    @:ASSERT(size(inp%hubbU, dim=1) == OSCC%mShell_)
    @:ASSERT(size(inp%hubbU, dim=2) == OSCC%nSpecies_)
    @:ASSERT(size(inp%tDampedShort) == OSCC%nSpecies_)
    @:ASSERT(allocated(inp%extCharges) .or. .not. allocated(inp%blurWidths))
#:call ASSERT_CODE
    if (allocated(inp%extCharges)) then
      @:ASSERT(size(inp%extCharges, dim=1) == 4)
      @:ASSERT(size(inp%extCharges, dim=2) > 0)
    end if
#:endcall ASSERT_CODE

    allocate(OSCC%invRMat_(OSCC%nAtom_, OSCC%nAtom_))
    allocate(OSCC%shiftPerAtom_(OSCC%nAtom_))
    allocate(OSCC%shiftPerL_(OSCC%mShell_, OSCC%nAtom_))
    allocate(OSCC%shortGamma_(0, 0, 0, 0))
    OSCC%tPeriodic_ = allocated(inp%latVecs)
    OSCC%tExtChrg_ = allocated(inp%extCharges)

    ! Initialize Hubbard U's
    allocate(OSCC%uniqHubbU_(OSCC%mShell_, OSCC%nSpecies_))
    allocate(OSCC%nHubbU_(OSCC%nSpecies_))
    allocate(OSCC%iHubbU_(OSCC%mShell_, OSCC%nSpecies_))
    OSCC%iHubbU_(:,:) = 0
    OSCC%iHubbU_(1,:) = 1
    OSCC%nHubbU_(:) = 1
    OSCC%uniqHubbU_(:,:) = 0.0_dp
    OSCC%uniqHubbU_(1,:) = inp%hubbU(1,:)
    do iSp1 = 1, OSCC%nSpecies_
      do iL = 2, inp%orb%nShell(iSp1)
        do iU1 = 1, OSCC%nHubbU_(iSp1)
          if (abs(inp%hubbU(iL,iSp1) - OSCC%uniqHubbU_(iU1,iSp1)) < MinHubDiff) then
            OSCC%iHubbU_(iL,iSp1) = iU1
            exit
          end if
        end do
        if (OSCC%iHubbU_(iL,iSp1) == 0) then
          OSCC%nHubbU_(iSp1) = OSCC%nHubbU_(iSp1) + 1
          OSCC%uniqHubbU_(OSCC%nHubbU_(iSp1),iSp1) = inp%hubbU(iL,iSp1)
          OSCC%iHubbU_(iL,iSp1) = OSCC%nHubbU_(iSp1)
        end if
      end do
    end do
    OSCC%mHubbU_ = maxval(OSCC%nHubbU_)

    ! Get cutoff for short range coulomb
    allocate(OSCC%shortCutoff_(OSCC%mHubbU_, OSCC%mHubbU_, OSCC%nSpecies_, OSCC%nSpecies_))
    OSCC%shortCutoff_(:,:,:,:) = 0.0_dp
    do iSp1 = 1, OSCC%nSpecies_
      do iSp2 = iSp1, OSCC%nSpecies_
        do iU1 = 1, OSCC%nHubbU_(iSp1)
          do iU2 = 1, OSCC%nHubbU_(iSp2)
            OSCC%shortCutoff_(iU2, iU1, iSp2, iSp1) = &
                & expGammaCutoff(OSCC%uniqHubbU_(iU2, iSp2), OSCC%uniqHubbU_(iU1, iSp1))
            OSCC%shortCutoff_(iU1, iU2, iSp1, iSp2) = OSCC%shortCutoff_(iU2, iU1, iSp2, iSp1)
          end do
        end do
      end do
    end do
    OSCC%cutoff_ = maxval(OSCC%shortCutoff_)

    ! Initialize Ewald summation for the periodic case
    if (OSCC%tPeriodic_) then
      OSCC%volume_ = inp%volume
      OSCC%tAutoEwald_ = inp%ewaldAlpha <= 0.0_dp
      if (OSCC%tAutoEwald_) then
        OSCC%alpha_ = getOptimalAlphaEwald(inp%latVecs, inp%recVecs, OSCC%volume_, tolEwald)
      else
        OSCC%alpha_ = inp%ewaldAlpha
      end if
      OSCC%maxREwald_ = getMaxREwald(OSCC%alpha_, tolEwald)
      maxGEwald = getMaxGEwald(OSCC%alpha_, OSCC%volume_, tolEwald)
      call getLatticePoints(OSCC%gLatPoint_, inp%recVecs, inp%latVecs/(2.0_dp*pi), maxGEwald, &
          & onlyInside=.true., reduceByInversion=.true., withoutOrigin=.true.)
      OSCC%gLatPoint_(:,:) = matmul(inp%recVecs, OSCC%gLatPoint_)
      OSCC%cutoff_ = max(OSCC%cutoff_, OSCC%maxREwald_)
    end if

    ! Number of neighbors for short range cutoff and real part of Ewald
    allocate(OSCC%nNeighShort_(OSCC%mHubbU_, OSCC%mHubbU_, OSCC%nSpecies_, OSCC%nAtom_))
    if (OSCC%tPeriodic_) then
      allocate(OSCC%nNeighEwald_(OSCC%nAtom_))
    end if

    ! Initialise external charges
    if (OSCC%tExtChrg_) then
      if (OSCC%tPeriodic_) then
        call init_ExtChrg(inp%extCharges, OSCC%nAtom_, inp%latVecs, inp%recVecs, OSCC%maxREwald_)
      else
        @:ASSERT(allocated(inp%blurWidths))
        call init_ExtChrg(inp%extCharges, OSCC%nAtom_, blurWidths=inp%blurWidths)
      end if
    end if

    OSCC%tChrgConstr_ = allocated(inp%chrgConstraints)
    if (OSCC%tChrgConstr_) then
      allocate(OSCC%chrgConstr_)
      call init(OSCC%chrgConstr_, inp%chrgConstraints, 2)
    end if
    OSCC%tThirdOrder_ = allocated(inp%thirdOrderOn)
    if (OSCC%tThirdOrder_) then
      allocate(OSCC%thirdOrder_)
      ! Factor 1/6 in the energy is put into the Hubbard derivatives
      call init(OSCC%thirdOrder_, inp%thirdOrderOn / 6.0_dp, 3)
    end if

    ! Initialise arrays for charge differences
    allocate(OSCC%deltaQ_(OSCC%mOrb_, OSCC%nAtom_))
    allocate(OSCC%deltaQPerLShell_(OSCC%mShell_, OSCC%nAtom_))
    allocate(OSCC%deltaQAtom_(OSCC%nAtom_))
    allocate(OSCC%deltaQUniqU_(OSCC%mHubbU_, OSCC%nAtom_))

    ! Initialise short range damping
    allocate(OSCC%tDampedShort_(OSCC%nSpecies_))
    OSCC%tDampedShort_(:) = inp%tDampedShort(:)
    OSCC%dampExp_ = inp%dampExp

    OSCC%tCoordUp_ = .false.
    OSCC%tChrgUp_ = .false.

    OSCC%tInitialised_ = .true.

  end subroutine init_SCC


  !> Returns a minimal cutoff for the neighborlist, which must be passed to various functions in
  !> this module.
  function getSCCCutoff(OSCC) result(cutoff)

    !> Module variables
    type(typSCC), intent(in) :: OSCC

    !> cutoff The neighborlists, passed to scc routines, should contain neighbour information at
    !> least up to that cutoff.
    real(dp) :: cutoff

    @:ASSERT(OSCC%tInitialised_)
    cutoff = OSCC%cutoff_

  end function getSCCCutoff


  !> Returns the currenty used alpha parameter of the Ewald-summation
  function getSCCEwaldPar(OSCC) result(alpha)

    !> Module variables
    type(typSCC), intent(in) :: OSCC

    !> Parameter in the Ewald summation.
    real(dp) :: alpha

    @:ASSERT(OSCC%tInitialised_)
    alpha = OSCC%alpha_

  end function getSCCEwaldPar


  !> Updates the number of neighbors for the SCC module (local).
  subroutine updateNNeigh_(species, neighList, OSCC)

    !> Species for each atom
    integer, intent(in) :: species(:)

    !> Neighbor list for the atoms in the system.
    type(TNeighborList), intent(in) :: neighList

    !> Module variables
    type(typSCC), intent(inout) :: OSCC

    integer :: iAt1, iSp2, iU1, iU2

    OSCC%nNeighShort_(:,:,:,:) = 0
    do iAt1 = 1, OSCC%nAtom_
      do iSp2 = 1, OSCC%nSpecies_
        do iU1 = 1, OSCC%nHubbU_(species(iAt1))
          do iU2 = 1, OSCC%nHubbU_(iSp2)
            OSCC%nNeighShort_(iU2, iU1, iSp2, iAt1) = &
                &getNrOfNeighbors(neighList, OSCC%shortCutoff_(iU2, iU1, iSp2, species(iAt1)), iAt1)
          end do
        end do
      end do
    end do
    if (OSCC%tPeriodic_) then
      call getNrOfNeighborsForAll(OSCC%nNeighEwald_, neighList, OSCC%maxREwald_)
    end if

  end subroutine updateNNeigh_


  !> Updates the atom coordinates for the SCC module.
  subroutine updateCoords_SCC(coord, species, neighList, img2CentCell, OSCC)

    !> New coordinates of the atoms
    real(dp), intent(in) :: coord(:,:)

    !> Species of the atoms (should not change during run)
    integer, intent(in) :: species(:)

    !> Neighbor list for the atoms.
    type(TNeighborList), intent(in) :: neighList

    !> Mapping to the central cell for the atoms
    integer, intent(in) :: img2CentCell(:)

    !> Module variables
    type(typSCC), intent(inout) :: OSCC

    @:ASSERT(OSCC%tInitialised_)

    call updateNNeigh_(species, neighList, OSCC)
    if (OSCC%tPeriodic_) then
      call invR(OSCC%invRMat_, OSCC%nAtom_, coord, OSCC%nNeighEwald_, neighList%iNeighbor, &
          &img2CentCell, OSCC%gLatPoint_, OSCC%alpha_, OSCC%volume_)
    else
      call invR(OSCC%invRMat_, OSCC%nAtom_, coord)
    end if
    call initGamma_(OSCC, coord, species, neighList%iNeighbor)

    if (OSCC%tExtChrg_) then
      if (OSCC%tPeriodic_) then
        call updateCoords_ExtChrg(coord, OSCC%gLatPoint_, OSCC%alpha_, OSCC%volume_)
      else
        call updateCoords_ExtChrg(coord)
      end if
    end if

    OSCC%tCoordUp_ = .true.
    OSCC%tChrgUp_ = .false.

  end subroutine updateCoords_SCC


  !> Updates the SCC module, if the lattice vectors had been changed
  subroutine updateLatVecs_SCC(latVec, recVec, vol, OSCC)

    !> New lattice vectors
    real(dp), intent(in) :: latVec(:,:)

    !> New reciprocal lattice vectors
    real(dp), intent(in) :: recVec(:,:)

    !> New volume
    real(dp), intent(in) :: vol

    !> Module variables
    type(typSCC), intent(inout) :: OSCC

    real(dp) :: maxGEwald

    @:ASSERT(OSCC%tInitialised_)
    @:ASSERT(OSCC%tPeriodic_)

    OSCC%volume_ = vol
    if (OSCC%tAutoEwald_) then
      OSCC%alpha_ = getOptimalAlphaEwald(latVec, recVec, OSCC%volume_, tolEwald)
      OSCC%maxREwald_ = getMaxREwald(OSCC%alpha_, tolEwald)
    end if
    maxGEwald = getMaxGEwald(OSCC%alpha_, OSCC%volume_, tolEwald)
    call getLatticePoints(OSCC%gLatPoint_, recVec, latVec/(2.0_dp*pi), maxGEwald, &
        &onlyInside=.true., reduceByInversion=.true., withoutOrigin=.true.)
    OSCC%gLatPoint_ = matmul(recVec, OSCC%gLatPoint_)
    OSCC%cutoff_ = max(OSCC%cutoff_, OSCC%maxREwald_)

    if (OSCC%tExtChrg_) then
      call updateLatVecs_extChrg(latVec, recVec, OSCC%maxREwald_)
    end if

    OSCC%tCoordUp_ = .false.
    OSCC%tChrgUp_ = .false.

  end subroutine updateLatVecs_SCC


  !> Updates the SCC module, if the charges have been changed
  subroutine updateCharges_SCC(qOrbital, q0, orb, species, iNeighbor, img2CentCell, OSCC)

    !> Orbital resolved charges
    real(dp), intent(in) :: qOrbital(:,:,:)

    !> Reference charge distribution (neutral atoms)
    real(dp), intent(in) :: q0(:,:,:)

    !> Contains information about the atomic orbitals in the system
    type(TOrbitals), intent(in) :: orb

    !> Species of the atoms (should not change during run)
    integer, intent(in) :: species(:)

    !> Neighbor indexes
    integer, intent(in) :: iNeighbor(0:,:)

    !> Mapping on atoms in the central cell
    integer, intent(in) :: img2CentCell(:)

    !> Resulting module variables
    type(typSCC), intent(inout) :: OSCC

    @:ASSERT(OSCC%tInitialised_)
    @:ASSERT(OSCC%tCoordUp_)

    call getNetCharges_(OSCC%nAtom_, OSCC%iHubbU_, species, orb, qOrbital, q0, OSCC%deltaQ_, &
        & OSCC%deltaQAtom_, OSCC%deltaQPerLShell_, OSCC%deltaQUniqU_)
    call buildShifts_(OSCC, orb, species, iNeighbor, img2CentCell, OSCC%deltaQAtom_, &
        & OSCC%deltaQPerLShell_, OSCC%shiftPerAtom_, OSCC%shiftPerL_)
    if (OSCC%tChrgConstr_) then
      call buildShift(OSCC%chrgConstr_, OSCC%deltaQAtom_)
    end if
    if (OSCC%tThirdOrder_) then
      call buildShift(OSCC%thirdOrder_, OSCC%deltaQAtom_)
    end if

    OSCC%tChrgUp_ = .true.

  end subroutine updateCharges_SCC


  !> Calculates various net charges needed by the SCC module.
  subroutine getNetCharges_(nAtom, iHubbU, species, orb, qOrbital, q0, dQ, dQAtom, dQShell, dQUniqU)

    !> Number of atoms in the system
    integer, intent(in) :: nAtom

    !> Mapping L-shell -> uniq U
    integer, intent(in) :: iHubbU(:,:)

    !> chemical species for atoms
    integer, intent(in) :: species(:)

    !> Contains information about the atomic orbitals in the system
    type(TOrbitals), intent(in) :: orb

    !> Orbital resolved charges
    real(dp), intent(in) :: qOrbital(:,:,:)

    !> Reference charge distribution (neutral atoms)
    real(dp), intent(in) :: q0(:,:,:)

    !> net charge for each orbital
    real(dp), intent(out) :: dQ(:,:)

    !> net charge for each atom
    real(dp), intent(out) :: dQAtom(:)

    !> net charge for each atomic shell
    real(dp), intent(out) :: dQShell(:,:)

    !> net charge for shells with the same U value on atoms
    real(dp), intent(out), optional :: dQUniqU(:,:)

    call getNetChargesPerOrbital_(qOrbital(:,:,1), q0(:,:,1), dQ)
    call getNetChargesPerAtom_(dQ, dQAtom)
    call getNetChargesPerLShell_(nAtom,species, orb, dQ, dQShell)
    if (present(dQUniqU)) then
      call getNetChargesPerUniqU_(nAtom, iHubbU, species, orb, dQShell, dQUniqU)
    end if

  end subroutine getNetCharges_


  !> net charges for each atomic orbital
  subroutine getNetChargesPerOrbital_(qOrbital, q0, deltaQ)

    !> orbital charges
    real(dp), intent(in) :: qOrbital(:,:)

    !> reference charges
    real(dp), intent(in) :: q0(:,:)

    !> resulting net charges
    real(dp), intent(out) :: deltaQ(:,:)

    deltaQ(:,:) = qOrbital(:,:) - q0(:,:)

  end subroutine getNetChargesPerOrbital_


  !> net charges per atom
  subroutine getNetChargesPerAtom_(deltaQ, deltaQAtom)

    !> net charges
    real(dp), intent(in) :: deltaQ(:,:)

    !> net charge per atom
    real(dp), intent(out) :: deltaQAtom(:)

    deltaQAtom(:) = sum(deltaQ(:,:), dim=1)

  end subroutine getNetChargesPerAtom_


  !> net charge per atomic shell
  subroutine getNetChargesPerLShell_(nAtom,species, orb, deltaQ, deltaQPerLShell)

    !> species of atom
    integer, intent(in) :: nAtom

    !> species of atom
    integer, intent(in) :: species(:)

    !> orbital information
    type(TOrbitals), intent(in) :: orb

    !> net charge for orbitals
    real(dp), intent(in) :: deltaQ(:,:)

    !> net charge per atomic shell
    real(dp), intent(out) :: deltaQPerLShell(:,:)

    integer :: iAt, iSp, iSh, iStart, iend

    deltaQPerLShell(:,:) = 0.0_dp
    do iAt = 1, nAtom
      iSp = species(iAt)
      do iSh = 1, orb%nShell(iSp)
        iStart = orb%posShell(iSh,iSp)
        iEnd = orb%posShell(iSh+1,iSp)-1
        deltaQPerLShell(iSh, iAt) = sum(deltaQ(iStart:iEnd, iAt))
      end do
    end do

  end subroutine getNetChargesPerLShell_


  !> net charges for orbitals with the same hubard U value on an atom
  subroutine getNetChargesPerUniqU_(nAtom, iHubbU, species, orb, deltaQPerLShell, deltaQUniqU)

    !> species of atom
    integer, intent(in) :: nAtom

    !> Mapping L-shell -> uniq U
    integer, intent(in) :: iHubbU(:,:)
    
    !> chemical species
    integer, intent(in) :: species(:)

    !> orbital information
    type(TOrbitals), intent(in) :: orb

    !> Net charge per l-shell of atom
    real(dp), intent(in) :: deltaQPerLShell(:,:)

    !> Net charge for unique groups of orbitals
    real(dp), intent(out) :: deltaQUniqU(:,:)

    integer :: iAt, iSp, iSh

    deltaQUniqU(:,:) = 0.0_dp
    do iAt = 1, nAtom
      iSp = species(iAt)
      do iSh = 1, orb%nShell(iSp)
        deltaQUniqU(iHubbU(iSh, iSp), iAt) =  deltaQUniqU(iHubbU(iSh, iSp), iAt) &
            & + deltaQPerLShell(iSh, iAt)
      end do
    end do

  end subroutine getNetChargesPerUniqU_


  !> Set up the storage and internal values for the short range part of Gamma.
  subroutine initGamma_(OSCC, coord, species, iNeighbor)

    !> Resulting module variables
    type(typSCC), intent(inout) :: OSCC

    !> List of coordinates
    real(dp), intent(in) :: coord(:,:)

    !> List of the species for each atom.
    integer, intent(in) :: species(:)

    !> Index of neighboring atoms for each atom.
    integer, intent(in) :: iNeighbor(0:,:)

    integer :: iAt1, iAt2, iU1, iU2, iNeigh, iSp1, iSp2
    real(dp) :: rab, u1, u2

    @:ASSERT(OSCC%tInitialised_)

    ! Reallocate shortGamma, if it does not contain enough neighbors
    if (size(OSCC%shortGamma_, dim=3) < maxval(OSCC%nNeighShort_)+1) then
      deallocate(OSCC%shortGamma_)
      allocate(OSCC%shortGamma_(OSCC%mHubbU_, OSCC%mHubbU_, 0:maxval(OSCC%nNeighShort_), &
          & OSCC%nAtom_))
    end if
    OSCC%shortGamma_(:,:,:,:) = 0.0_dp

    ! some additional symmetry not used, as the value of gamma for atoms
    ! interacting with themselves is the same for all atoms of the same species
    do iAt1 = 1, OSCC%nAtom_
      iSp1 = species(iAt1)
      do iNeigh = 0, maxval(OSCC%nNeighShort_(:,:,:, iAt1))
        iAt2 = iNeighbor(iNeigh, iAt1)
        iSp2 = species(iAt2)
        rab = sqrt(sum((coord(:,iAt1) - coord(:,iAt2))**2))
        do iU1 = 1, OSCC%nHubbU_(species(iAt1))
          u1 = OSCC%uniqHubbU_(iU1, iSp1)
          do iU2 = 1, OSCC%nHubbU_(species(iAt2))
            u2 = OSCC%uniqHubbU_(iU2, iSp2)
            if (iNeigh <= OSCC%nNeighShort_(iU2,iU1,iSp2,iAt1)) then
              if (OSCC%tDampedShort_(iSp1) .or. OSCC%tDampedShort_(iSp2)) then
                OSCC%shortGamma_(iU2 ,iU1, iNeigh, iAt1) = expGammaDamped(rab, u2, u1, &
                    & OSCC%dampExp_)
              else
                OSCC%shortGamma_(iU2 ,iU1, iNeigh, iAt1) = expGamma(rab, u2, u1)
              end if
            end if
          end do
        end do
      end do
    end do

  end subroutine initGamma_


  !> Routine for returning lower triangle of atomic resolved gamma as a matrix
  subroutine getAtomicGammaMatrix(gammamat, OSCC, iNeighbor, img2CentCell)

    !> Atom resolved gamma
    real(dp), intent(out) :: gammamat(:,:)

    !> Resulting module variables
    type(typSCC), intent(in) :: OSCC

    !> neighbours of atoms
    integer, intent(in) :: iNeighbor(0:,:)

    !> index array to fold images to central cell
    integer, intent(in) :: img2CentCell(:)

    integer :: iAt1, iAt2, iAt2f, iNeigh

    @:ASSERT(OSCC%tInitialised_)
    @:ASSERT(OSCC%tCoordUp_)
    @:ASSERT(all(shape(gammamat) == [ OSCC%nAtom_, OSCC%nAtom_ ]))
    @:ASSERT(all(OSCC%nHubbU_ == 1))

    gammamat(:,:) = OSCC%invRMat_
    do iAt1 = 1, OSCC%nAtom_
      do iNeigh = 0, maxval(OSCC%nNeighShort_(:,:,:, iAt1))
        iAt2 = iNeighbor(iNeigh, iAt1)
        iAt2f = img2CentCell(iAt2)
        gammamat(iAt2f, iAt1) = gammamat(iAt2f, iAt1) - OSCC%shortGamma_(1, 1, iNeigh, iAt1)
      end do
    end do

  end subroutine getAtomicGammaMatrix


  !> Calculate  the derivative of the short range part of Gamma.
  subroutine addGammaPrime_(force, OSCC, coord, species, iNeighbor, img2CentCell)

    !> force vector to add the short-range part of gamma contribution
    real(dp), intent(inout) :: force(:,:)

    !> Resulting module variables
    type(typSCC), intent(in) :: OSCC

    !> list of coordinates
    real(dp), intent(in) :: coord(:,:)

    !> List of the species for each atom.
    integer, intent(in) :: species(:)

    !> Index of neighboring atoms for each atom.
    integer, intent(in) :: iNeighbor(0:,:)

    !> Image of each atom in the central cell.
    integer, intent(in) :: img2CentCell(:)

    integer :: iAt1, iAt2, iAt2f, iU1, iU2, iNeigh, ii, iSp1, iSp2
    real(dp) :: rab, tmpGammaPrime, u1, u2

    @:ASSERT(size(force,dim=1) == 3)
    @:ASSERT(size(force,dim=2) == OSCC%nAtom_)
    @:ASSERT(OSCC%tInitialised_)

    ! some additional symmetry not used
    do iAt1 = 1, OSCC%nAtom_
      iSp1 = species(iAt1)
      do iNeigh = 1, maxval(OSCC%nNeighShort_(:,:,:, iAt1))
        iAt2 = iNeighbor(iNeigh, iAt1)
        iAt2f = img2CentCell(iAt2)
        iSp2 = species(iAt2f)
        rab = sqrt(sum((coord(:,iAt1) - coord(:,iAt2))**2))
        do iU1 = 1, OSCC%nHubbU_(species(iAt1))
          u1 = OSCC%uniqHubbU_(iU1, iSp1)
          do iU2 = 1, OSCC%nHubbU_(species(iAt2f))
            u2 = OSCC%uniqHubbU_(iU2, species(iAt2f))
            if (iNeigh <= OSCC%nNeighShort_(iU2,iU1,species(iAt2f),iAt1)) then
              if (OSCC%tDampedShort_(iSp1) .or. OSCC%tDampedShort_(iSp2)) then
                tmpGammaPrime = expGammaDampedPrime(rab, u2, u1, OSCC%dampExp_)
              else
                tmpGammaPrime = expGammaPrime(rab, u2, u1)
              end if
              do ii = 1,3
                force(ii,iAt1) = force(ii,iAt1) - OSCC%deltaQUniqU_(iU1,iAt1) * &
                    & OSCC%deltaQUniqU_(iU2,iAt2f)*tmpGammaPrime*(coord(ii,iAt1) &
                    & - coord(ii,iAt2))/rab
                force(ii,iAt2f) = force(ii,iAt2f) + OSCC%deltaQUniqU_(iU1,iAt1) * &
                    & OSCC%deltaQUniqU_(iU2,iAt2f)*tmpGammaPrime*(coord(ii,iAt1) &
                    & - coord(ii,iAt2))/rab
              end do
            end if
          end do
        end do
      end do
    end do

  end subroutine addGammaPrime_


  !> Calculate  the derivative of the short range part of Gamma.
  subroutine addSTGammaPrime_(st, OSCC, coord, species, iNeighbor, img2CentCell)

    !> Stress tensor component to add the short-range part of the gamma contribution
    real(dp), intent(out) :: st(:,:)

    !> Resulting module variables
    type(typSCC), intent(in) :: OSCC

    !> list of coordinates
    real(dp), intent(in) :: coord(:,:)

    !> List of the species for each atom.
    integer, intent(in) :: species(:)

    !> Index of neighboring atoms for each atom.
    integer, intent(in) :: iNeighbor(0:,:)

    !> Image of each atom in the central cell.
    integer, intent(in) :: img2CentCell(:)

    integer :: iAt1, iAt2, iAt2f, iU1, iU2, iNeigh, ii, jj, iSp1, iSp2
    real(dp) :: rab, tmpGammaPrime, u1, u2
    real(dp) :: intermed(3), vect(3)

    @:ASSERT(all(shape(st)==(/3,3/)))
    @:ASSERT(OSCC%tInitialised_)

    st(:,:) = 0.0_dp
    ! some additional symmetry not used
    do iAt1 = 1, OSCC%nAtom_
      iSp1 = species(iAt1)
      do iNeigh = 1, maxval(OSCC%nNeighShort_(:,:,:, iAt1))
        iAt2 = iNeighbor(iNeigh, iAt1)
        iAt2f = img2CentCell(iAt2)
        iSp2 = species(iAt2f)
        vect(:) = coord(:,iAt1) - coord(:,iAt2)
        rab = sqrt(sum((vect)**2))
        intermed(:) = 0.0_dp
        do iU1 = 1, OSCC%nHubbU_(species(iAt1))
          u1 = OSCC%uniqHubbU_(iU1, iSp1)
          do iU2 = 1, OSCC%nHubbU_(species(iAt2f))
            u2 = OSCC%uniqHubbU_(iU2, species(iAt2f))
            if (iNeigh <= OSCC%nNeighShort_(iU2,iU1,species(iAt2f),iAt1)) then
              if (OSCC%tDampedShort_(iSp1) .or. OSCC%tDampedShort_(iSp2)) then
                tmpGammaPrime = expGammaDampedPrime(rab, u2, u1, OSCC%dampExp_)
              else
                tmpGammaPrime = expGammaPrime(rab, u2, u1)
              end if
              do ii = 1,3
                intermed(ii) = intermed(ii) &
                    & - OSCC%deltaQUniqU_(iU1,iAt1) * OSCC%deltaQUniqU_(iU2,iAt2f) &
                    & *tmpGammaPrime*vect(ii)/rab
              end do
            end if
          end do
        end do
        if (iAt2f /= iAt1) then
          do ii = 1, 3
            do jj = 1, 3
              st(jj,ii) = st(jj,ii) + (vect(jj)*intermed(ii) + intermed(jj)*vect(ii))
            end do
          end do
        else
          do ii = 1, 3
            do jj = 1, 3
              st(jj,ii) = st(jj,ii) + 0.5_dp * (vect(jj)*intermed(ii) + intermed(jj)*vect(ii))
            end do
          end do
        end if
      end do
    end do

    st(:,:) = st(:,:) / OSCC%volume_

  end subroutine addSTGammaPrime_


  !> Calculates the contribution of the charge consistent part to the energy per atom.
  subroutine getEnergyPerAtom_SCC(eSCC, OSCC)

    !> The SCC contribution to the energy
    real(dp), intent(out) :: eSCC(:)

    !> Resulting module variables
    type(typSCC), intent(in) :: OSCC

    @:ASSERT(OSCC%tInitialised_)
    @:ASSERT(size(eSCC) == OSCC%nAtom_)

    eSCC(:) = 0.5_dp * (OSCC%shiftPerAtom_ * OSCC%deltaQAtom_ &
        & + sum(OSCC%shiftPerL_ * OSCC%deltaQPerLShell_, dim=1))
    if (OSCC%tExtChrg_) then
      call addEnergyPerAtom_ExtChrg(OSCC%deltaQAtom_, eSCC)
    end if
    if (OSCC%tChrgConstr_) then
      call addEnergyPerAtom(OSCC%chrgConstr_, eSCC, OSCC%deltaQAtom_)
    end if
    if (OSCC%tThirdOrder_) then
      call addEnergyPerAtom(OSCC%thirdOrder_, eSCC, OSCC%deltaQAtom_)
    end if

  end subroutine getEnergyPerAtom_SCC


  !> Calculates SCC energy contribution using the linearized XLBOMD form.
  !> Note: When SCC is driven in XLBOMD mode, the charges should NOT be updated after diagonalizing
  !> the Hamiltonian, so the charge stored in the module are the input (auxiliary) charges, used to
  !> build the Hamiltonian.  However, the linearized energy expession needs also the output charges,
  !> therefore these are passed in as an extra variable.
  subroutine getEnergyPerAtom_SCC_Xlbomd(OSCC, species, orb, qOut, q0, eSCC)

    !> Resulting module variables
    type(typSCC), intent(in) :: OSCC

    !> atomic species
    integer, intent(in) :: species(:)

    !> orbital information
    type(TOrbitals), intent(in) :: orb

    !> output charges
    real(dp), intent(in) :: qOut(:,:,:)

    !> reference charges
    real(dp), intent(in) :: q0(:,:,:)

    !> energy contributions
    real(dp), intent(out) :: eSCC(:)

    real(dp), allocatable :: dQOut(:,:), dQOutAtom(:), dQOutShell(:,:)

    @:ASSERT(OSCC%tInitialised_)
    @:ASSERT(size(eSCC) == OSCC%nAtom_)

    allocate(dQOut(orb%mOrb, OSCC%nAtom_))
    allocate(dQOutAtom(OSCC%nAtom_))
    allocate(dQOutShell(OSCC%mShell_, OSCC%nAtom_))

    call getNetCharges_(OSCC%nAtom_, OSCC%iHubbU_, species, orb, qOut, q0, dQOut, dQOutAtom, &
        & dQOutShell)

    ! 1/2 sum_A (2 q_A - n_A) * shift(n_A)
    eSCC(:) = 0.5_dp * (OSCC%shiftPerAtom_ * (2.0_dp * dQOutAtom - OSCC%deltaQAtom_) &
        & + sum(OSCC%shiftPerL_ * (2.0_dp * dQOutShell - OSCC%deltaQPerLShell_), dim=1))

    if (OSCC%tExtChrg_) then
      call error("XLBOMD not working with external charges yet")
      !call addEnergyPerAtom_ExtChrg(OSCC%deltaQAtom_, eSCC)
    end if
    if (OSCC%tChrgConstr_) then
      call error("XLBOMD not working with charge constraints yet")
      !call addEnergyPerAtom(OSCC%chrgConstr_, eSCC, OSCC%deltaQAtom_)
    end if
    if (OSCC%tThirdOrder_) then
      call error("XLBOMD not working with third order yet")
      !call addEnergyPerAtom(OSCC%thirdOrder_, eSCC, OSCC%deltaQAtom_)
    end if

  end subroutine getEnergyPerAtom_SCC_Xlbomd


  !> Calculates the contribution of the charge consistent part to the forces for molecules/clusters,
  !> which is not covered in the term with the shift vectors.
  subroutine addForceDCSCC(force, OSCC, species, iNeighbor, img2CentCell, coord,chrgForce)

    !> has force contribution added
    real(dp), intent(inout) :: force(:,:)

    !> Resulting module variables
    type(typSCC), intent(in) :: OSCC

    !> Species for each atom.
    integer, intent(in) :: species(:)

    !> List of neighbors for each atom.
    integer, intent(in) :: iNeighbor(0:,:)

    !> Indexing of images of the atoms in the central cell.
    integer, intent(in) :: img2CentCell(:)

    !> List of coordinates
    real(dp), intent(in) :: coord(:,:)

    !> Force contribution due to the external charges, which is not contained in the term with the
    !> shift vectors.
    real(dp), intent(inout), optional :: chrgForce(:,:)

    @:ASSERT(OSCC%tInitialised_)
    @:ASSERT(size(force,dim=1) == 3)
    @:ASSERT(size(force,dim=2) == OSCC%nAtom_)
    @:ASSERT(present(chrgForce) .eqv. OSCC%tExtChrg_)

    ! Short-range part of gamma contribution
    call addGammaPrime_(force,OSCC,coord,species,iNeighbor,img2CentCell)

    ! 1/R contribution
    if (OSCC%tPeriodic_) then
      call addInvRPrime(force, OSCC%nAtom_, coord, OSCC%nNeighEwald_, iNeighbor, &
          & img2CentCell, OSCC%gLatPoint_, OSCC%alpha_, OSCC%volume_, OSCC%deltaQAtom_)
      if (OSCC%tExtChrg_) then
        call addForceDCSCC_ExtChrg(force, chrgForce, coord, OSCC%deltaQAtom_, OSCC%gLatPoint_, &
            & OSCC%alpha_, OSCC%volume_)
      end if
    else
      call addInvRPrime(force, OSCC%nAtom_, coord, OSCC%deltaQAtom_)
      if (OSCC%tExtChrg_) then
        call addForceDCSCC_ExtChrg(force, chrgForce, coord, OSCC%deltaQAtom_)
      end if
    end if

  end subroutine addForceDCSCC


  !> Calculates the contribution of the stress tensor which is not covered in the term with the
  !> shift vectors.
  subroutine addStressDCSCC(st, OSCC, species, iNeighbor, img2CentCell, coord) !,chrgForce)

    !> Add stress tensor contribution to this
    real(dp), intent(inout) :: st(:,:)

    !> Resulting module variables
    type(typSCC), intent(in) :: OSCC

    !> Species for each atom.
    integer, intent(in) :: species(:)

    !> List of neighbors for each atom.
    integer, intent(in) :: iNeighbor(0:,:)

    !> Indexing of images of the atoms in the central cell.
    integer, intent(in) :: img2CentCell(:)

    !> List of coordinates
    real(dp), intent(in) :: coord(:,:)

    real(dp) :: stTmp(3,3)

    @:ASSERT(OSCC%tInitialised_)
    @:ASSERT(OSCC%tPeriodic_)
    @:ASSERT(all(shape(st)==(/3,3/)))

    stTmp = 0.0_dp

    ! Short-range part of gamma contribution
    call addSTGammaPrime_(stTmp,OSCC,coord,species,iNeighbor,img2CentCell)

    st(:,:) = st(:,:) - 0.5_dp * stTmp(:,:)

    ! 1/R contribution
    ! call invRstress

    stTmp = 0.0_dp
    call invR_stress(stTmp, OSCC%nAtom_, coord, OSCC%nNeighEwald_, iNeighbor,img2CentCell, &
        & OSCC%gLatPoint_, OSCC%alpha_, OSCC%volume_, OSCC%deltaQAtom_)

    st(:,:) = st(:,:) - 0.5_dp * stTmp(:,:)

    ! if (tExtChrg_) then
    ! ????
    ! end if

  end subroutine addStressDCSCC


  !> Constructs the shift vectors for the SCC contributions.
  !> The full shift vector must be constructed by adding shiftAtom and shiftShell accordingly.
  subroutine buildShifts_(OSCC, orb, species, iNeighbor, img2CentCell, dQAtom, dQShell, shiftAtom, &
      & shiftShell)

    !> Resulting module variables
    type(typSCC), intent(in) :: OSCC

    !> Contains information about the atomic orbitals in the system
    type(TOrbitals), intent(in) :: orb

    !> List of the species for each atom.
    integer, intent(in) :: species(:)

    !> List of surrounding neighbours for each atom.
    integer, intent(in) :: iNeighbor(0:,:)

    !> Image of each atom in the central cell.
    integer, intent(in) :: img2CentCell(:)

    !> Net population per atom.
    real(dp), intent(in) :: dQAtom(:)

    !> Net population per l-shell.
    real(dp), intent(in) :: dQShell(:,:)

    !> Shift vector for each atom.
    real(dp), intent(out) :: shiftAtom(:)

    !> Shift vector per l-shell.
    real(dp), intent(out) :: shiftShell(:,:)

    integer :: iAt1, iSp1, iSh1, iU1, iNeigh, iAt2f, iSp2, iSh2, iU2

    @:ASSERT(OSCC%tInitialised_)

    ! 1/R contribution [shiftAtom(A) = \sum_B 1/R_AB * (Q_B - Q0_B)]
    shiftAtom(:) = 0.0_dp
    call hemv(shiftAtom, OSCC%invRMat_, dQAtom,'L')

    ! Contribution of the short range part of gamma to the shift
    ! sgamma'_{A,l} = sum_B sum_{u\in B} S(U(A,l),u)*q_u
    shiftShell(:,:) = 0.0_dp
    do iAt1 = 1, OSCC%nAtom_
      iSp1 = species(iAt1)
      do iSh1 = 1, orb%nShell(iSp1)
        iU1 = OSCC%iHubbU_(iSh1, iSp1)
        do iNeigh = 0, maxval(OSCC%nNeighShort_(:,:,:,iAt1))
          iAt2f = img2CentCell(iNeighbor(iNeigh, iAt1))
          iSp2 = species(iAt2f)
          do iSh2 = 1, orb%nShell(iSp2)
            iU2 = OSCC%iHubbU_(iSh2, iSp2)
            shiftShell(iSh1, iAt1) = shiftShell(iSh1, iAt1) &
                &- OSCC%shortGamma_(iU2, iU1, iNeigh, iAt1) * dQShell(iSh2, iAt2f)
            if (iAt2f /= iAt1) then
              shiftShell(iSh2, iAt2f) = shiftShell(iSh2, iAt2f) &
                  &- OSCC%shortGamma_(iU2, iU1, iNeigh, iAt1) * dQShell(iSh1, iAt1)
            end if
          end do
        end do
      end do
    end do

  end subroutine buildShifts_


  !> Returns the shift per atom coming from the SCC part
  subroutine getShiftPerAtom(shift, OSCC)

    !> Contains the shift on exit.
    real(dp), intent(out) :: shift(:)

    !> Module variables
    type(typSCC), intent(in) :: OSCC

    @:ASSERT(OSCC%tInitialised_)
    @:ASSERT(size(shift) == size(OSCC%shiftPerAtom_))

    shift = 0.0_dp
    shift = OSCC%shiftPerAtom_
    if (OSCC%tExtChrg_) then
      call addShiftPerAtom_ExtChrg(shift)
    end if
    if (OSCC%tChrgConstr_) then
      call addShiftPerAtom(OSCC%chrgConstr_, shift)
    end if
    if (OSCC%tThirdOrder_) then
      call addShiftPerAtom(OSCC%thirdOrder_, shift)
    end if
    
  end subroutine getShiftPerAtom


  !> Returns the shift per L contribution of the SCC.
  subroutine getShiftPerL(shift, OSCC)

    !> Contains the shift on exit.
    real(dp), intent(out) :: shift(:,:)

    !> Module variables
    type(typSCC), intent(in) :: OSCC

    @:ASSERT(OSCC%tInitialised_)
    @:ASSERT(size(shift,dim=1) == size(OSCC%shiftPerL_,dim=1))
    @:ASSERT(size(shift,dim=2) == size(OSCC%shiftPerL_,dim=2))
    shift = 0.0_dp
    shift = OSCC%shiftPerL_

  end subroutine getShiftPerL


  !> Returns the equivalency relations between orbitals of the atoms. If transfering charge between
  !> the orbitals does not change the electrostatic energy, they are considered equivalent.
  subroutine SCC_getOrbitalEquiv(OSCC, orb, species, equiv)

    !> Resulting module variables
    type(typSCC), intent(in) :: OSCC

    !> Contains information about the atomic orbitals in the system
    type(TOrbitals), intent(in) :: orb

    !> Type of each atom (nAtom).
    integer, intent(in) :: species(:)

    !> The vector describing the equivalence on return.
    integer, intent(out) :: equiv(:,:,:)

    integer :: iAt, iOrb, iS, iSp
    integer :: nSpin, shift

    nSpin = size(equiv, dim=3)

    @:ASSERT(OSCC%tInitialised_)
    @:ASSERT(size(species) == OSCC%nAtom_)
    @:ASSERT(size(equiv, dim=1) == orb%mOrb)
    @:ASSERT(all(shape(equiv) == (/ orb%mOrb, OSCC%nAtom_, nSpin /)))

    equiv(:,:,:) = 0
    shift = 0
    do iAt = 1, OSCC%nAtom_
      iSp = species(iAt)
      do iOrb = 1, orb%nOrbSpecies(iSp)
        equiv(iOrb, iAt, 1) = OSCC%iHubbU_(orb%iShellOrb(iOrb, iSp), iSp) + shift
      end do
      shift = shift + maxval(OSCC%iHubbU_(:, iSp))
    end do
    do iS = 2, nSpin
      equiv(:,:,iS) = equiv(:,:,1)
    end do

  end subroutine SCC_getOrbitalEquiv


  !> Calculate the "double counting" force term using linearized XLBOMD form.
  !> Note: When SCC is driven in XLBOMD mode, the charges should NOT be updated after diagonalizing
  !> the Hamiltonian, so the charge stored in the module are the input (auxiliary) charges, used to
  !> build the Hamiltonian.  However, the linearized energy expession needs also the output charges,
  !> therefore these are passed in as an extra variable.
  subroutine addForceDCSCC_Xlbomd(OSCC, species, orb, iNeighbor, img2CentCell, coord, qOrbitalOut, &
      & q0, force)

    !> Resulting module variables
    type(typSCC), intent(in) :: OSCC

    !> atomic species
    integer, intent(in) :: species(:)

    !> orbital information
    type(TOrbitals), intent(in) :: orb

    !> neighbours surrounding each atom
    integer, intent(in) :: iNeighbor(0:,:)

    !> index from image atoms to central cell
    integer, intent(in) :: img2CentCell(:)

    !> coordinates of all atoms
    real(dp), intent(in) :: coord(:,:)

    !> output charges
    real(dp), intent(in) :: qOrbitalOut(:,:,:)

    !> reference charges
    real(dp), intent(in) :: q0(:,:,:)

    !> Force terms are added to this
    real(dp), intent(inout) :: force(:,:)

    real(dp), allocatable :: dQOut(:,:), dQOutAtom(:)
    real(dp), allocatable :: dQOutLShell(:,:), dQOutUniqU(:,:)

    allocate(dQOut(OSCC%mOrb_, OSCC%nAtom_))
    allocate(dQOutAtom(OSCC%nAtom_))
    allocate(dQOutLShell(OSCC%mShell_, OSCC%nAtom_))
    allocate(dQOutUniqU(OSCC%mHubbU_, OSCC%nAtom_))

    call getNetCharges_(OSCC%nAtom_, OSCC%iHubbU_, species, orb, qOrbitalOut, q0, dQOut, &
        & dQOutAtom, dQOutLShell, dQOutUniqU)

    ! Short-range part of gamma contribution
    call addGammaPrimeXlbomd_(OSCC, OSCC%deltaQUniqU_, dQOutUniqU, coord, species, iNeighbor, &
        & img2CentCell, force)

    ! 1/R contribution
    if (OSCC%tPeriodic_) then
      call addInvRPrimeXlbomd(OSCC%nAtom_, coord, OSCC%nNeighEwald_, iNeighbor, img2CentCell, &
          & OSCC%gLatPoint_, OSCC%alpha_, OSCC%volume_, OSCC%deltaQAtom_, dQOutAtom, force)
      if (OSCC%tExtChrg_) then
        call error("XLBOMD with external charges does not work yet!")
        !call addForceDCSCC_ExtChrg(force, chrgForce, coord, deltaQAtom_, gLatPoint_, alpha_, &
        !& volume_)
      end if
    else
      call addInvRPrimeXlbomd(OSCC%nAtom_, coord, OSCC%deltaQAtom_, dQOutAtom, force)
      if (OSCC%tExtChrg_) then
        call error("XLBOMD with external charges does not work yet!")
        !call addForceDCSCC_ExtChrg(force, chrgForce, coord, deltaQAtom_)
      end if
    end if

  end subroutine addForceDCSCC_Xlbomd


  !> Calculate the derivative of the short range contributions using the linearized XLBOMD
  !> formulation with auxiliary charges.
  subroutine addGammaPrimeXlbomd_(OSCC, dQInUniqU, dQOutUniqU, coord, species, iNeighbor, &
      & img2CentCell, force)

    !> Resulting module variables
    type(typSCC), intent(in) :: OSCC

    !> Input charges
    real(dp), intent(in) :: dQInUniqU(:,:)

    !> output charges
    real(dp), intent(in) :: dQOutUniqU(:,:)

    !> atomic coordinates
    real(dp), intent(in) :: coord(:,:)

    !> chemical species
    integer, intent(in) :: species(:)

    !> neighbours around atoms
    integer, intent(in) :: iNeighbor(0:,:)

    !> image to real atom indexing
    integer, intent(in) :: img2CentCell(:)

    !> term to add force contributions to
    real(dp), intent(inout) :: force(:,:)

    integer :: iAt1, iAt2, iAt2f, iU1, iU2, iNeigh, iSp1, iSp2
    real(dp) :: rab, tmpGammaPrime, u1, u2, prefac
    real(dp) :: contrib(3)

    @:ASSERT(OSCC%tInitialised_)
    @:ASSERT(size(force,dim=1) == 3)
    @:ASSERT(size(force,dim=2) == OSCC%nAtom_)

    do iAt1 = 1, OSCC%nAtom_
      iSp1 = species(iAt1)
      do iNeigh = 1, maxval(OSCC%nNeighShort_(:,:,:, iAt1))
        iAt2 = iNeighbor(iNeigh, iAt1)
        iAt2f = img2CentCell(iAt2)
        iSp2 = species(iAt2f)
        rab = sqrt(sum((coord(:,iAt1) - coord(:,iAt2))**2))
        do iU1 = 1, OSCC%nHubbU_(species(iAt1))
          u1 = OSCC%uniqHubbU_(iU1, iSp1)
          do iU2 = 1, OSCC%nHubbU_(species(iAt2f))
            u2 = OSCC%uniqHubbU_(iU2, species(iAt2f))
            if (iNeigh <= OSCC%nNeighShort_(iU2,iU1,species(iAt2f),iAt1)) then
              if (OSCC%tDampedShort_(iSp1) .or. OSCC%tDampedShort_(iSp2)) then
                tmpGammaPrime = expGammaDampedPrime(rab, u2, u1, OSCC%dampExp_)
              else
                tmpGammaPrime = expGammaPrime(rab, u2, u1)
              end if
              prefac = dQOutUniqU(iU1, iAt1) * dQInUniqU(iU2, iAt2f) &
                  & + dQInUniqU(iU1, iAt1) * dQOutUniqU(iU2, iAt2f) &
                  & - dQInUniqU(iU1, iAt1) * dQInUniqU(iU2, iAt2f)
              contrib(:) = prefac * tmpGammaPrime / rab  * (coord(:,iAt2) - coord(:,iAt1))
              force(:,iAt1) = force(:,iAt1) + contrib
              force(:,iAt2f) = force(:,iAt2f) - contrib
            end if
          end do
        end do
      end do
    end do

  end subroutine addGammaPrimeXlbomd_

end module scc
