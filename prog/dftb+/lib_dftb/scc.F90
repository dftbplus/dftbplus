!--------------------------------------------------------------------------------------------------!
!  DFTB+: general package for performing fast atomistic simulations                                !
!  Copyright (C) 2017  DFTB+ developers group                                                      !
!                                                                                                  !
!  See the LICENSE file for terms of usage and distribution.                                       !
!--------------------------------------------------------------------------------------------------!

#:include 'common.fypp'

!!* Functions and local variables for the SCC calculation.
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

  public :: TSCCInit, init_SCC
  public :: updateCoords_SCC, updateLatVecs_SCC, updateCharges_SCC
  public :: getSCCCutoff, getEnergyPerAtom_SCC, getEnergyPerAtom_SCC_Xlbomd
  public :: addForceDCSCC, getShiftPerAtom, getShiftPerL, getSCCEwaldPar
  public :: SCC_getOrbitalEquiv, addStressDCSCC
  public :: getAtomicGammaMatrix
  public :: addForceDCSCC_Xlbomd

  !!* Data necessary for initialize the SCC module
  type TSCCInit
    type(TOrbitals), pointer :: orb
    real(dp), allocatable :: hubbU(:,:)
    logical, allocatable :: tDampedShort(:)
    real(dp) :: dampExp = 0.0_dp
    real(dp), allocatable :: latVecs(:,:)
    real(dp), allocatable :: recVecs(:,:)
    real(dp) :: volume = 0.0_dp
    real(dp), allocatable :: extCharges(:,:)
    real(dp), allocatable :: blurWidths(:)
    real(dp), allocatable :: chrgConstraints(:,:)
    real(dp), allocatable :: thirdOrderOn(:,:)
    real(dp) :: ewaldAlpha = 0.0_dp    ! if > 0 -> manual setting for alpha
  end type TSCCInit

  real(dp), parameter :: tolEwald = 1.0e-9_dp !* tolerance for Ewald


  !! Private module variables (suffixed with "_" for clarity)
  logical :: tInitialised_ = .false.      ! If module is initialised
  logical :: tCoordUp_                    ! If coordinates updated at least once
  logical :: tChrgUp_                     ! If charges updated at least once
  integer :: nAtom_, nSpecies_             ! Nr. of atoms and species
  integer :: mOrb_, mShell_       ! Max nr. of orbitals/shells per atom
  integer :: mHubbU_    ! Maximal number of unique U values
  real(dp), allocatable :: invRMat_(:,:)  ! Stores 1/r between atom pairs
  real(dp), allocatable :: shiftPerAtom_(:)      ! Shift vector per atom
  real(dp), allocatable :: shiftPerL_(:,:)       ! Shift vector per l-shell
  real(dp), allocatable :: shortGamma_(:,:,:,:)  ! Short range interaction
  real(dp), allocatable :: shortCutoff_(:,:,:,:) ! Cutoff for short range int.
  real(dp) :: cutoff_                            ! Maximal cutoff
  real(dp), allocatable :: gLatPoint_(:,:) ! Lattice points for reciprocal Ewald
  integer, allocatable :: nNeighShort_(:,:,:,:) ! Nr. of neighbors for short
                                                ! range interaction
  integer, allocatable :: nNeighEwald_(:)       ! Nr. of neigh for real Ewald
  real(dp) :: maxREwald_                        ! Cutoff for real Ewald
  real(dp) :: alpha_                            ! Parameter for Ewald
  real(dp) :: volume_                           ! Cell volume
  real(dp), allocatable :: uniqHubbU_(:,:)      ! Uniq Us per species
  integer, allocatable :: nHubbU_(:)            ! Nr. of uniq Us per species
  integer, allocatable :: iHubbU_(:,:)          ! Mapping L-shell -> uniq U
  logical :: tExtChrg_                          ! Are external charges present?

  real(dp), allocatable :: deltaQ_(:,:)          ! Negative net charge
  real(dp), allocatable :: deltaQPerLShell_(:,:) ! Negative net charge per shell
  real(dp), allocatable :: deltaQAtom_(:)        ! Negative net charge per atom
  real(dp), allocatable :: deltaQUniqU_(:,:)     ! Negative net charge per U

  logical, allocatable :: tDampedShort_(:)      ! Damped short range? (nSpecies)
  real(dp) :: dampExp_                          ! Damping exponent

  logical :: tPeriodic_                          ! Is the system periodic?

  logical :: tChrgConstr_                       ! Shifts due charge constrains?
  type(OChrgConstr), allocatable, save :: chrgConstr_        ! Object for charge constraints
  logical :: tThirdOrder_                       ! Shifts due to 3rd order
  type(OChrgConstr), allocatable, save :: thirdOrder_
  logical :: tAutoEwald_


contains

  !!* Initialises the SCC module
  !!* @param nAtom Number of atoms in the system
  !!* @param nSpecies Number of species in the system
  !!* @param nShell Nr. of shells for each species
  !!* @param mShell Maximal nr. of shells for any species
  !!* @param mOrb   Maximal nr. of orbitals for any species
  !!* @param hubbU The Hubbard Us read from the SK files
  !!* @param tDampedShort True for species which needs damped interaction
  !!* @param latVec Lattice vectors (if the system is periodic)
  !!* @param recVec Reciprocal vectors (if the system is periodic)
  !!* @param volume Volume of the supercell (if the system is periodic)
  !!* @param extCharges Coordinate and charge of external point charges.
  !!* @param blurWidths widths of Gaussian blurring on charges (non-periodic
  !!* only)
  subroutine init_SCC(inp)
    type(TSCCInit), intent(inout) :: inp

    integer :: iSp1, iSp2, iU1, iU2, iL
    real(dp) :: maxGEwald

    nSpecies_ = size(inp%orb%nOrbSpecies)
    nAtom_ = size(inp%orb%nOrbAtom)
    mShell_ = inp%orb%mShell
    mOrb_ = inp%orb%mOrb

    @:ASSERT(.not. tInitialised_)
    @:ASSERT(allocated(inp%latVecs) .eqv. allocated(inp%recVecs))
    @:ASSERT(allocated(inp%latVecs) .eqv. (inp%volume > 0.0_dp))
    @:ASSERT(size(inp%hubbU, dim=1) == mShell_)
    @:ASSERT(size(inp%hubbU, dim=2) == nSpecies_)
    @:ASSERT(size(inp%tDampedShort) == nSpecies_)
    @:ASSERT(allocated(inp%extCharges) .or. .not. allocated(inp%blurWidths))
  #:call ASSERT_CODE
    if (allocated(inp%extCharges)) then
      @:ASSERT(size(inp%extCharges, dim=1) == 4)
      @:ASSERT(size(inp%extCharges, dim=2) > 0)
    end if
  #:endcall ASSERT_CODE

    allocate(invRMat_(nAtom_, nAtom_))
    allocate(shiftPerAtom_(nAtom_))
    allocate(shiftPerL_(mShell_, nAtom_))
    allocate(shortGamma_(0, 0, 0, 0))
    tPeriodic_ = allocated(inp%latVecs)
    tExtChrg_ = allocated(inp%extCharges)

    !! Initialize Hubbard U's
    allocate(uniqHubbU_(mShell_, nSpecies_))
    allocate(nHubbU_(nSpecies_))
    allocate(iHubbU_(mShell_, nSpecies_))
    iHubbU_(:,:) = 0
    iHubbU_(1,:) = 1
    nHubbU_(:) = 1
    uniqHubbU_(:,:) = 0.0_dp
    uniqHubbU_(1,:) = inp%hubbU(1,:)
    do iSp1 = 1, nSpecies_
      do iL = 2, inp%orb%nShell(iSp1)
        do iU1 = 1, nHubbU_(iSp1)
          if (abs(inp%hubbU(iL,iSp1) - uniqHubbU_(iU1,iSp1)) < MinHubDiff) then
            iHubbU_(iL,iSp1) = iU1
            exit
          end if
        end do
        if (iHubbU_(iL,iSp1) == 0) then
          nHubbU_(iSp1) = nHubbU_(iSp1) + 1
          uniqHubbU_(nHubbU_(iSp1),iSp1) = inp%hubbU(iL,iSp1)
          iHubbU_(iL,iSp1) = nHubbU_(iSp1)
        end if
      end do
    end do
    mHubbU_ = maxval(nHubbU_)

    !! Get cutoff for short range coulomb
    allocate(shortCutoff_(mHubbU_, mHubbU_, nSpecies_, nSpecies_))
    shortCutoff_(:,:,:,:) = 0.0_dp
    do iSp1 = 1, nSpecies_
      do iSp2 = iSp1, nSpecies_
        do iU1 = 1, nHubbU_(iSp1)
          do iU2 = 1, nHubbU_(iSp2)
            shortCutoff_(iU2, iU1, iSp2, iSp1) = &
                & expGammaCutoff(uniqHubbU_(iU2, iSp2), uniqHubbU_(iU1, iSp1))
            shortCutoff_(iU1, iU2, iSp1, iSp2) = &
                & shortCutoff_(iU2, iU1, iSp2, iSp1)
          end do
        end do
      end do
    end do
    cutoff_ = maxval(shortCutoff_)

    !! Initialize Ewald summation for the periodic case
    if (tPeriodic_) then
      volume_ = inp%volume
      tAutoEwald_ = inp%ewaldAlpha <= 0.0_dp
      if (tAutoEwald_) then
        alpha_ = getOptimalAlphaEwald(inp%latVecs, inp%recVecs, volume_, &
            &tolEwald)
      else
        alpha_ = inp%ewaldAlpha
      end if
      maxREwald_ = getMaxREwald(alpha_, tolEwald)
      maxGEwald = getMaxGEwald(alpha_, volume_, tolEwald)
      call getLatticePoints(gLatPoint_, inp%recVecs, inp%latVecs/(2.0_dp*pi), &
          &maxGEwald, onlyInside=.true., reduceByInversion=.true., &
          &withoutOrigin=.true.)
      gLatPoint_(:,:) = matmul(inp%recVecs, gLatPoint_)
      cutoff_ = max(cutoff_, maxREwald_)
    end if

    !! Number of neighbors for short range cutoff and real part of Ewald
    allocate(nNeighShort_(mHubbU_, mHubbU_, nSpecies_, nAtom_))
    if (tPeriodic_) then
      allocate(nNeighEwald_(nAtom_))
    end if

    !! Initialise external charges
    if (tExtChrg_) then
      if (tPeriodic_) then
        call init_ExtChrg(inp%extCharges, nAtom_, inp%latVecs, inp%recVecs, &
            &maxREwald_)
      else
        @:ASSERT(allocated(inp%blurWidths))
        call init_ExtChrg(inp%extCharges, nAtom_, blurWidths=inp%blurWidths)
      end if
    end if

    tChrgConstr_ = allocated(inp%chrgConstraints)
    if (tChrgConstr_) then
      allocate(chrgConstr_)
      call init(chrgConstr_, inp%chrgConstraints, 2)
    end if
    tThirdOrder_ = allocated(inp%thirdOrderOn)
    if (tThirdOrder_) then
      allocate(thirdOrder_)
      ! Factor 1/6 in the energy is put into the Hubbard derivatives
      call init(thirdOrder_, inp%thirdOrderOn / 6.0_dp, 3)
    end if

    !! Initialise arrays for charge differences
    allocate(deltaQ_(mOrb_, nAtom_))
    allocate(deltaQPerLShell_(mShell_, nAtom_))
    allocate(deltaQAtom_(nAtom_))
    allocate(deltaQUniqU_(mHubbU_, nAtom_))

    !! Initialise short range damping
    allocate(tDampedShort_(nSpecies_))
    tDampedShort_(:) = inp%tDampedShort(:)
    dampExp_ = inp%dampExp

    tInitialised_ = .true.
    tCoordUp_ = .false.
    tChrgUp_ = .false.

  end subroutine init_SCC


  !!* Returns a minimal cutoff for the neighborlist, which must be passed
  !!* to various functions in this module.
  !!* @return cutoff The neighborlists, passed to scc routines, should contain
  !!*   neigbhor information at least up to that cutoff.
  function getSCCCutoff() result(cutoff)
    real(dp) :: cutoff

    @:ASSERT(tInitialised_)
    cutoff = cutoff_

  end function getSCCCutoff


  !!* Returns the currenty used alpha parameter of the Ewald-summation
  !!* @return Parameter in the Ewald summation.
  function getSCCEwaldPar() result(alpha)
    real(dp) :: alpha

    @:ASSERT(tInitialised_)
    alpha = alpha_

  end function getSCCEwaldPar


  !!* Updates the number of neighbors for the SCC module (local).
  !!* @param species Species for each atom
  !!* @param neighList Neighbor list for the atoms in the system.
  subroutine updateNNeigh_(species, neighList)
    integer, intent(in) :: species(:)
    type(TNeighborList), intent(in) :: neighList

    integer :: iAt1, iSp2, iU1, iU2

    nNeighShort_(:,:,:,:) = 0
    do iAt1 = 1, nAtom_
      do iSp2 = 1, nSpecies_
        do iU1 = 1, nHubbU_(species(iAt1))
          do iU2 = 1, nHubbU_(iSp2)
            nNeighShort_(iU2, iU1, iSp2, iAt1) = &
                & getNrOfNeighbors(neighList, &
                & shortCutoff_(iU2, iU1, iSp2, species(iAt1)), iAt1)
          end do
        end do
      end do
    end do
    if (tPeriodic_) then
      call getNrOfNeighborsForAll(nNeighEwald_, neighList, maxREwald_)
    end if

  end subroutine updateNNeigh_



  !!* Updates the atom coordinates for the SCC module.
  !!* @param coord New coordinates of the atoms
  !!* @param species Species of the atoms (should not change during run)
  !!* @param mAngSpecies Maximal angular momentum of the species (should not
  !!* change during run)
  !!* @param neighList Neighbor list for the atoms.
  !!* @param img2CentCell Mapping to the central cell for the atoms
  subroutine updateCoords_SCC(coord, species, neighList, img2CentCell)
    real(dp), intent(in) :: coord(:,:)
    integer,  intent(in) :: species(:)
    type(TNeighborList), intent(in) :: neighList
    integer,  intent(in) :: img2CentCell(:)

    @:ASSERT(tInitialised_)

    call updateNNeigh_(species, neighList)
    if (tPeriodic_) then
      call invR(invRMat_, nAtom_, coord, nNeighEwald_, neighList%iNeighbor, &
          &img2CentCell, gLatPoint_, alpha_, volume_)
    else
      call invR(invRMat_, nAtom_, coord)
    end if
    call initGamma_(coord, species, neighList%iNeighbor)

    if (tExtChrg_) then
      if (tPeriodic_) then
        call updateCoords_ExtChrg(coord, gLatPoint_, alpha_, volume_)
      else
        call updateCoords_ExtChrg(coord)
      end if
    end if

    tCoordUp_ = .true.
    tChrgUp_ = .false.

  end subroutine updateCoords_SCC


  !!* Updates the SCC module, if the lattice vectors had been changed
  !!* @param latVec  New lattice vectors
  !!* @param recVec  New reciprocal lattice vectors
  !!* @param vol  New volume
  subroutine updateLatVecs_SCC(latVec, recVec, vol)
    real(dp), intent(in) :: latVec(:,:), recVec(:,:)
    real(dp), intent(in) :: vol

    real(dp) :: maxGEwald

    @:ASSERT(tInitialised_)
    @:ASSERT(tPeriodic_)

    volume_ = vol
    if (tAutoEwald_) then
      alpha_ = getOptimalAlphaEwald(latVec, recVec, volume_, tolEwald)
      maxREwald_ = getMaxREwald(alpha_, tolEwald)
    end if
    maxGEwald = getMaxGEwald(alpha_, volume_, tolEwald)
    call getLatticePoints(gLatPoint_, recVec, latVec/(2.0_dp*pi), maxGEwald, &
        &onlyInside=.true., reduceByInversion=.true., withoutOrigin=.true.)
    gLatPoint_ = matmul(recVec, gLatPoint_)
    cutoff_ = max(cutoff_, maxREwald_)

    if (tExtChrg_) then
      call updateLatVecs_extChrg(latVec, recVec, maxREwald_)
    end if

    tCoordUp_ = .false.
    tChrgUp_ = .false.

  end subroutine updateLatVecs_SCC


  !!* Updates the SCC module, if the charges have been changed
  !!* @param qOrbital Orbital resolved charges
  !!* @param q0 Reference charge distribution (neutral atoms)
  !!* @param mAngSpecies Maximal angular momentum of the species (should not
  !!* change during run)
  !!* @param orb Contains information about the atomic orbitals in the system
  !!* @param species Species of the atoms (should not change during run)
  !!* @param iNeighbor Neighbor indexes
  !!* @param img2CentCell Mapping on atoms in the centrall cell
  subroutine updateCharges_SCC(qOrbital, q0, orb, species, iNeighbor, &
      &img2CentCell)
    real(dp), intent(in) :: qOrbital(:,:,:)
    real(dp), intent(in) :: q0(:,:,:)
    type(TOrbitals), intent(in) :: orb
    integer, intent(in) :: species(:)
    integer, intent(in) :: iNeighbor(0:,:)
    integer, intent(in) :: img2CentCell(:)

    @:ASSERT(tInitialised_)
    @:ASSERT(tCoordUp_)

    call getNetCharges_(species, orb, qOrbital, q0, deltaQ_, deltaQAtom_, &
        & deltaQPerLShell_, deltaQUniqU_)
    call buildShifts_(orb, species, iNeighbor, img2CentCell, deltaQAtom_, &
        & deltaQPerLShell_, shiftPerAtom_, shiftPerL_)
    if (tChrgConstr_) then
      call buildShift(chrgConstr_, deltaQAtom_)
    end if
    if (tThirdOrder_) then
      call buildShift(thirdOrder_, deltaQAtom_)
    end if

    tChrgUp_ = .true.

  end subroutine updateCharges_SCC


  !> Calculates various net charges needed by the SCC module.
  !!
  subroutine getNetCharges_(species, orb, qOrbital, q0, dQ, dQAtom, dQShell, &
      & dQUniqU)
    integer, intent(in) :: species(:)
    type(TOrbitals), intent(in) :: orb
    real(dp), intent(in) :: qOrbital(:,:,:), q0(:,:,:)
    real(dp), intent(out) :: dQ(:,:), dQAtom(:), dQShell(:,:)
    real(dp), intent(out), optional :: dQUniqU(:,:)

    call getNetChargesPerOrbital_(qOrbital(:,:,1), q0(:,:,1), dQ)
    call getNetChargesPerAtom_(dQ, dQAtom)
    call getNetChargesPerLShell_(species, orb, dQ, dQShell)
    if (present(dQUniqU)) then
      call getNetChargesPerUniqU_(species, orb, dQShell, dQUniqU)
    end if

  end subroutine getNetCharges_



  subroutine getNetChargesPerOrbital_(qOrbital, q0, deltaQ)
    real(dp), intent(in) :: qOrbital(:,:), q0(:,:)
    real(dp), intent(out) :: deltaQ(:,:)

    deltaQ(:,:) = qOrbital(:,:) - q0(:,:)

  end subroutine getNetChargesPerOrbital_



  subroutine getNetChargesPerAtom_(deltaQ, deltaQAtom)
    real(dp), intent(in) :: deltaQ(:,:)
    real(dp), intent(out) :: deltaQAtom(:)

    deltaQAtom(:) = sum(deltaQ(:,:), dim=1)

  end subroutine getNetChargesPerAtom_



  subroutine getNetChargesPerLShell_(species, orb, deltaQ, deltaQPerLShell)
    integer, intent(in) :: species(:)
    type(TOrbitals), intent(in) :: orb
    real(dp), intent(in) :: deltaQ(:,:)
    real(dp), intent(out) :: deltaQPerLShell(:,:)

    integer :: iAt, iSp, iSh, iStart, iend

    deltaQPerLShell(:,:) = 0.0_dp
    do iAt = 1, nAtom_
      iSp = species(iAt)
      do iSh = 1, orb%nShell(iSp)
        iStart = orb%posShell(iSh,iSp)
        iEnd = orb%posShell(iSh+1,iSp)-1
        deltaQPerLShell(iSh, iAt) = sum(deltaQ(iStart:iEnd, iAt))
      end do
    end do

  end subroutine getNetChargesPerLShell_



  subroutine getNetChargesPerUniqU_(species, orb, deltaQPerLShell, deltaQUniqU)
    integer, intent(in) :: species(:)
    type(TOrbitals), intent(in) :: orb
    real(dp), intent(in) :: deltaQPerLShell(:,:)
    real(dp), intent(out) :: deltaQUniqU(:,:)

    integer :: iAt, iSp, iSh

    deltaQUniqU(:,:) = 0.0_dp
    do iAt = 1, nAtom_
      iSp = species(iAt)
      do iSh = 1, orb%nShell(iSp)
        deltaQUniqU(iHubbU_(iSh, iSp), iAt) = &
            & deltaQUniqU(iHubbU_(iSh, iSp), iAt) &
            &+ deltaQPerLShell(iSh, iAt)
      end do
    end do

  end subroutine getNetChargesPerUniqU_



  !!* Set up the storage and internal values for the short range part of Gamma.
  !!* @param coord        List of coordinates
  !!* @param species       List of the species for each atom.
  !!* @param iNeighbor    Index of neighboring atoms for each atom.
  subroutine initGamma_(coord, species, iNeighbor)
    real(dp), intent(in) :: coord(:,:)
    integer,  intent(in) :: species(:)
    integer,  intent(in) :: iNeighbor(0:,:)

    integer :: iAt1, iAt2, iU1, iU2, iNeigh, iSp1, iSp2
    real(dp) :: rab, u1, u2

    ! Reallocate shortGamma, if it does not contain enough neighbors
    if (size(shortGamma_, dim=3) < maxval(nNeighShort_)+1) then
      deallocate(shortGamma_)
      allocate(shortGamma_(mHubbU_, mHubbU_, 0:maxval(nNeighShort_), nAtom_))
    end if
    shortGamma_(:,:,:,:) = 0.0_dp

    ! some additional symmetry not used, as the value of gamma for atoms
    ! interacting with themselves is the same for all atoms of the same species
    do iAt1 = 1, nAtom_
      iSp1 = species(iAt1)
      do iNeigh = 0, maxval(nNeighShort_(:,:,:, iAt1))
        iAt2 = iNeighbor(iNeigh, iAt1)
        iSp2 = species(iAt2)
        rab = sqrt(sum((coord(:,iAt1) - coord(:,iAt2))**2))
        do iU1 = 1, nHubbU_(species(iAt1))
          u1 = uniqHubbU_(iU1, iSp1)
          do iU2 = 1, nHubbU_(species(iAt2))
            u2 = uniqHubbU_(iU2, iSp2)
            if (iNeigh <= nNeighShort_(iU2,iU1,iSp2,iAt1)) then
              if (tDampedShort_(iSp1) .or. tDampedShort_(iSp2)) then
                shortGamma_(iU2 ,iU1, iNeigh, iAt1) = &
                    &expGammaDamped(rab, u2, u1, dampExp_)
              else
                shortGamma_(iU2 ,iU1, iNeigh, iAt1) = expGamma(rab, u2, u1)
              end if
            end if
          end do
        end do
      end do
    end do

  end subroutine initGamma_

  !!* Routine for lower triangle of atomic resolved gamma as a matrix
  !!* @param gammamat Atom resolved gamma
  !!* @param iNeighbor neighbours of atoms
  !!* @param img2CentCell index array to fold images to central cell
  subroutine getAtomicGammaMatrix(gammamat, iNeighbor, img2CentCell)
    real(dp), intent(out) :: gammamat(:,:)
    integer, intent(in) :: iNeighbor(0:,:), img2CentCell(:)

    integer :: iAt1, iAt2, iAt2f, iNeigh

    @:ASSERT(tInitialised_)
    @:ASSERT(tCoordUp_)
    @:ASSERT(all(shape(gammamat) == [ nAtom_, nAtom_ ]))
    @:ASSERT(all(nHubbU_ == 1))

    gammamat(:,:) = invRMat_
    do iAt1 = 1, nAtom_
      do iNeigh = 0, maxval(nNeighShort_(:,:,:, iAt1))
        iAt2 = iNeighbor(iNeigh, iAt1)
        iAt2f = img2CentCell(iAt2)
        gammamat(iAt2f, iAt1) = gammamat(iAt2f, iAt1)&
            & - shortGamma_(1, 1, iNeigh, iAt1)
      end do
    end do

  end subroutine getAtomicGammaMatrix


  !!* Calculate  the derivative of the short range part of Gamma.
  !!* @param force        force vector to add the short-range part of gamma
  !!* contribution
  !!* @param coord        list of coordinates
  !!* @param species       List of the species for each atom.
  !!* @param iNeighbor    Index of neighboring atoms for each atom.
  !!* @param img2CentCell Image of each atom in the central cell.
  subroutine addGammaPrime_(force, coord, species, iNeighbor, img2CentCell)
    real(dp), intent(inout) :: force(:,:)
    real(dp), intent(in)    :: coord(:,:)
    integer,  intent(in)    :: species(:)
    integer,  intent(in)    :: iNeighbor(0:,:)
    integer,  intent(in)    :: img2CentCell(:)

    integer :: iAt1, iAt2, iAt2f, iU1, iU2, iNeigh, ii, iSp1, iSp2
    real(dp) :: rab, tmpGammaPrime, u1, u2

    @:ASSERT(size(force,dim=1) == 3)
    @:ASSERT(size(force,dim=2) == nAtom_)

    ! some additional symmetry not used
    do iAt1 = 1, nAtom_
      iSp1 = species(iAt1)
      do iNeigh = 1, maxval(nNeighShort_(:,:,:, iAt1))
        iAt2 = iNeighbor(iNeigh, iAt1)
        iAt2f = img2CentCell(iAt2)
        iSp2 = species(iAt2f)
        rab = sqrt(sum((coord(:,iAt1) - coord(:,iAt2))**2))
        do iU1 = 1, nHubbU_(species(iAt1))
          u1 = uniqHubbU_(iU1, iSp1)
          do iU2 = 1, nHubbU_(species(iAt2f))
            u2 = uniqHubbU_(iU2, species(iAt2f))
            if (iNeigh <= nNeighShort_(iU2,iU1,species(iAt2f),iAt1)) then
              if (tDampedShort_(iSp1) .or. tDampedShort_(iSp2)) then
                tmpGammaPrime = expGammaDampedPrime(rab, u2, u1, dampExp_)
              else
                tmpGammaPrime = expGammaPrime(rab, u2, u1)
              end if
              do ii = 1,3
                force(ii,iAt1) = force(ii,iAt1) - &
                    & deltaQUniqU_(iU1,iAt1) * &
                    & deltaQUniqU_(iU2,iAt2f) &
                    & *tmpGammaPrime*(coord(ii,iAt1) - coord(ii,iAt2))/rab
                force(ii,iAt2f) = force(ii,iAt2f)&
                    & +deltaQUniqU_(iU1,iAt1) * &
                    & deltaQUniqU_(iU2,iAt2f) &
                    & *tmpGammaPrime*(coord(ii,iAt1) - coord(ii,iAt2))/rab
              end do
            end if
          end do
        end do
      end do
    end do

  end subroutine addGammaPrime_


  !!* Calculate  the derivative of the short range part of Gamma.
  !!* @param st           Stress tensor component to add the short-range
  !!*  part of the gamma contribution
  !!* @param coord        list of coordinates
  !!* @param species       List of the species for each atom.
  !!* @param iNeighbor    Index of neighboring atoms for each atom.
  !!* @param img2CentCell Image of each atom in the central cell.
  subroutine addSTGammaPrime_(st, coord, species, iNeighbor, img2CentCell)
    real(dp), intent(out)   :: st(:,:)
    real(dp), intent(in)    :: coord(:,:)
    integer,  intent(in)    :: species(:)
    integer,  intent(in)    :: iNeighbor(0:,:)
    integer,  intent(in)    :: img2CentCell(:)

    integer  :: iAt1, iAt2, iAt2f, iU1, iU2, iNeigh, ii, jj, iSp1, iSp2
    real(dp) :: rab, tmpGammaPrime, u1, u2
    real(dp) :: intermed(3), vect(3)

    @:ASSERT(all(shape(st)==(/3,3/)))

    st(:,:) = 0.0_dp
    ! some additional symmetry not used
    do iAt1 = 1, nAtom_
      iSp1 = species(iAt1)
      do iNeigh = 1, maxval(nNeighShort_(:,:,:, iAt1))
        iAt2 = iNeighbor(iNeigh, iAt1)
        iAt2f = img2CentCell(iAt2)
        iSp2 = species(iAt2f)
        vect(:) = coord(:,iAt1) - coord(:,iAt2)
        rab = sqrt(sum((vect)**2))
        intermed(:) = 0.0_dp
        do iU1 = 1, nHubbU_(species(iAt1))
          u1 = uniqHubbU_(iU1, iSp1)
          do iU2 = 1, nHubbU_(species(iAt2f))
            u2 = uniqHubbU_(iU2, species(iAt2f))
            if (iNeigh <= nNeighShort_(iU2,iU1,species(iAt2f),iAt1)) then
              if (tDampedShort_(iSp1) .or. tDampedShort_(iSp2)) then
                tmpGammaPrime = expGammaDampedPrime(rab, u2, u1, dampExp_)
              else
                tmpGammaPrime = expGammaPrime(rab, u2, u1)
              end if
              do ii = 1,3
                intermed(ii) = intermed(ii) - &
                    & deltaQUniqU_(iU1,iAt1) * &
                    & deltaQUniqU_(iU2,iAt2f) &
                    & *tmpGammaPrime*vect(ii)/rab
              end do
            end if
          end do
        end do
        if (iAt2f /= iAt1) then
          do ii = 1, 3
            do jj = 1, 3
              st(jj,ii) = st(jj,ii) &
                  & + (vect(jj)*intermed(ii) + intermed(jj)*vect(ii))
            end do
          end do
        else
          do ii = 1, 3
            do jj = 1, 3
              st(jj,ii) = st(jj,ii) + 0.5_dp * &
                  & (vect(jj)*intermed(ii) + intermed(jj)*vect(ii))
            end do
          end do
        end if
      end do
    end do

    st(:,:) = st(:,:) / volume_

  end subroutine addSTGammaPrime_


  !!* Calculates the contribution of the charge consistent part to the energy
  !!* per atom.
  !!* @param eSCC         The SCC contribution to the energy
  !!* @param species       Species for each atom.
  !!* @param iNeighbor    List of neighbors for each atom.
  !!* @param img2CentCell Image of the atoms in the central cell.
  subroutine getEnergyPerAtom_SCC(eSCC)
    real(dp), intent(out) :: eSCC(:)

    @:ASSERT(tInitialised_)
    @:ASSERT(size(eSCC) == nAtom_)

    eSCC(:) = 0.5_dp * (shiftPerAtom_ * deltaQAtom_ &
        & + sum(shiftPerL_ * deltaQPerLShell_, dim=1))
    if (tExtChrg_) then
      call addEnergyPerAtom_ExtChrg(deltaQAtom_, eSCC)
    end if
    if (tChrgConstr_) then
      call addEnergyPerAtom(chrgConstr_, eSCC, deltaQAtom_)
    end if
    if (tThirdOrder_) then
      call addEnergyPerAtom(thirdOrder_, eSCC, deltaQAtom_)
    end if

  end subroutine getEnergyPerAtom_SCC


  !> Calculates SCC energy contribution using the linearized XLBOMD form.
  !!
  !! Note: When SCC is driven in XLBOMD mode, the charges should NOT be updated
  !! after diagonalizing the Hamiltonian, so the charge stored in the module are
  !! the ingoing (auxiliary) charges, used to build the Hamiltonian.  However,
  !! the linearized energy expession needs also the outgoing ones,
  !! therefore, there must be passed extra.
  !!
  subroutine getEnergyPerAtom_SCC_Xlbomd(species, orb, qOut, q0, eSCC)
    integer, intent(in) :: species(:)
    type(TOrbitals), intent(in) :: orb
    real(dp), intent(in) :: qOut(:,:,:), q0(:,:,:)
    real(dp), intent(out) :: eSCC(:)

    real(dp), allocatable :: dQOut(:,:), dQOutAtom(:), dQOutShell(:,:)

    @:ASSERT(tInitialised_)
    @:ASSERT(size(eSCC) == nAtom_)

    allocate(dQOut(orb%mOrb, nAtom_))
    allocate(dQOutAtom(nAtom_))
    allocate(dQOutShell(mShell_, nAtom_))

    call getNetCharges_(species, orb, qOut, q0, dQOut, dQOutAtom, dQOutShell)

    ! 1/2 sum_A (2 q_A - n_A) * shift(n_A)
    eSCC(:) = 0.5_dp * (shiftPerAtom_ * (2.0_dp * dQOutAtom - deltaQAtom_) &
        & + sum(shiftPerL_ * (2.0_dp * dQOutShell - deltaQPerLShell_), dim=1))

    if (tExtChrg_) then
      call error("XLBOMD not working with external charges yet")
      call addEnergyPerAtom_ExtChrg(deltaQAtom_, eSCC)
    end if
    if (tChrgConstr_) then
      call error("XLBOMD not working with charge constraints yet")
      call addEnergyPerAtom(chrgConstr_, eSCC, deltaQAtom_)
    end if
    if (tThirdOrder_) then
      call error("XLBOMD not working with third order yet")
      call addEnergyPerAtom(thirdOrder_, eSCC, deltaQAtom_)
    end if

  end subroutine getEnergyPerAtom_SCC_Xlbomd



  !!* Calculates the contribution of the charge consitent part to the forces
  !!* for molecules/clusters, which is not covered in the term with the shift
  !!* vectors.
  !!* @param force        Force contribution
  !!* @param mAngSpecies   Maximal angular momentum per species.
  !!* @param species       Species for each atom.
  !!* @param iNeighbor    List of neighbors for each atom.
  !!* @param img2CentCell Image of the atoms in the central cell.
  !!* @param coord        List of coordinates
  !!* @param chrgForce    Force contribution due to the external charges, which
  !!* is not contained in the term with the shift vectors.
  subroutine addForceDCSCC(force, species, iNeighbor, img2CentCell, &
      & coord,chrgForce)
    real(dp), intent(inout) :: force(:,:)
    integer,  intent(in)    :: species(:)
    integer,  intent(in)    :: iNeighbor(0:,:)
    integer,  intent(in)    :: img2CentCell(:)
    real(dp), intent(in)    :: coord(:,:)
    real(dp), intent(inout), optional :: chrgForce(:,:)


    @:ASSERT(size(force,dim=1) == 3)
    @:ASSERT(size(force,dim=2) == nAtom_)
    @:ASSERT(present(chrgForce) .eqv. tExtChrg_)

    ! Short-range part of gamma contribution
    call addGammaPrime_(force,coord,species,iNeighbor,img2CentCell)

    ! 1/R contribution
    if (tPeriodic_) then
      call addInvRPrime(force, nAtom_, coord, nNeighEwald_, iNeighbor, &
          & img2CentCell, gLatPoint_, alpha_, volume_, deltaQAtom_)
      if (tExtChrg_) then
        call addForceDCSCC_ExtChrg(force, chrgForce, coord, deltaQAtom_, &
            &gLatPoint_, alpha_, volume_)
      end if
    else
      call addInvRPrime(force, nAtom_, coord, deltaQAtom_)
      if (tExtChrg_) then
        call addForceDCSCC_ExtChrg(force, chrgForce, coord, deltaQAtom_)
      end if
    end if

  end subroutine addForceDCSCC



  !!* Calculates the contribution of the charge consitent part to the forces
  !!* for molecules/clusters, which is not covered in the term with the shift
  !!* vectors.
  !!* @param force        Force contribution
  !!* @param mAngSpecies   Maximal angular momentum per species.
  !!* @param species       Species for each atom.
  !!* @param iNeighbor    List of neighbors for each atom.
  !!* @param img2CentCell Image of the atoms in the central cell.
  !!* @param coord        List of coordinates
  subroutine addStressDCSCC(st, species, iNeighbor, img2CentCell, &
      & coord)!,chrgForce)
    real(dp), intent(inout) :: st(:,:)
    integer,  intent(in)    :: species(:)
    integer,  intent(in)    :: iNeighbor(0:,:)
    integer,  intent(in)    :: img2CentCell(:)
    real(dp), intent(in)    :: coord(:,:)

    real(dp) :: stTmp(3,3)

    @:ASSERT(tPeriodic_)
    @:ASSERT(all(shape(st)==(/3,3/)))

    stTmp = 0.0_dp

    ! Short-range part of gamma contribution
    call addSTGammaPrime_(stTmp,coord,species,iNeighbor,img2CentCell)

    st(:,:) = st(:,:) - 0.5_dp * stTmp(:,:)

    ! 1/R contribution
    ! call invRstress

    stTmp = 0.0_dp
    call invR_stress(stTmp, nAtom_, coord, nNeighEwald_, &
        & iNeighbor,img2CentCell, gLatPoint_, &
        & alpha_, volume_, deltaQAtom_)

    st(:,:) = st(:,:) - 0.5_dp * stTmp(:,:)

    ! if (tExtChrg_) then
    ! ????
    ! end if

  end subroutine addStressDCSCC


  !> Constructs the shift vectors for the SCC contributions
  !! \param orb  Contains information about the atomic orbitals in the system
  !! \param species  List of the species for each atom.
  !! \param iNeighbor  List of surrounding neighbours for each atom.
  !! \param img2CentCell  Image of each atom in the central cell.
  !! \param dQAtom  Net population per atom.
  !! \param dQShell  Net population per l-shell.
  !! \param shiftAtom  Shift vector for each atom.
  !! \param shiftShell  Shift vector per l-shell.
  !!
  !! \note The full shift vector must be constructed by adding shiftAtom
  !! and shiftShell accordingly.
  !!
  subroutine buildShifts_(orb, species, iNeighbor, img2CentCell, dQAtom, &
      & dQShell, shiftAtom, shiftShell)
    type(TOrbitals), intent(in) :: orb
    integer,  intent(in) :: species(:), iNeighbor(0:,:), img2CentCell(:)
    real(dp), intent(in) :: dQAtom(:), dQShell(:,:)
    real(dp), intent(out) :: shiftAtom(:), shiftShell(:,:)

    integer :: iAt1, iSp1, iSh1, iU1, iNeigh, iAt2f, iSp2, iSh2, iU2

    @:ASSERT(tInitialised_)

    !! 1/R contribution [shiftAtom(A) = \sum_B 1/R_AB * (Q_B - Q0_B)]
    shiftAtom(:) = 0.0_dp
    call hemv(shiftAtom, invRMat_, dQAtom,'L')

    !! Contribution of the short range part of gamma to the shift
    !! sgamma'_{A,l} = sum_B sum_{u\in B} S(U(A,l),u)*q_u
    shiftShell(:,:) = 0.0_dp
    do iAt1 = 1, nAtom_
      iSp1 = species(iAt1)
      do iSh1 = 1, orb%nShell(iSp1)
        iU1 = iHubbU_(iSh1, iSp1)
        do iNeigh = 0, maxval(nNeighShort_(:,:,:,iAt1))
          iAt2f = img2CentCell(iNeighbor(iNeigh, iAt1))
          iSp2 = species(iAt2f)
          do iSh2 = 1, orb%nShell(iSp2)
            iU2 = iHubbU_(iSh2, iSp2)
            shiftShell(iSh1, iAt1) = shiftShell(iSh1, iAt1) &
                &- shortGamma_(iU2, iU1, iNeigh, iAt1) * dQShell(iSh2, iAt2f)
            if (iAt2f /= iAt1) then
              shiftShell(iSh2, iAt2f) = shiftShell(iSh2, iAt2f) &
                  &- shortGamma_(iU2, iU1, iNeigh, iAt1) * dQShell(iSh1, iAt1)
            end if
          end do
        end do
      end do
    end do

  end subroutine buildShifts_



  !!* Returns the shift per atom coming from the SCC part (with a spin
  !!* index)
  !!* @param shift Contains the shift on exit.
  subroutine getShiftPerAtom(shift)
    real(dp), intent(out) :: shift(:,:)

    @:ASSERT(size(shift,dim=1) == size(shiftPerAtom_,dim=1))

    shift(:,:) = 0.0_dp
    shift(:,1) = shiftPerAtom_
    if (tExtChrg_) then
      call addShiftPerAtom_ExtChrg(shift(:,1))
    end if
    if (tChrgConstr_) then
      call addShiftPerAtom(chrgConstr_, shift(:,1))
    end if
    if (tThirdOrder_) then
      call addShiftPerAtom(thirdOrder_, shift(:,1))
    end if

  end subroutine getShiftPerAtom



  !!* Returns the shift per L contribution of the SCC. (with a spin index)
  !!* @param shift Contains the shift on exit.
  subroutine getShiftPerL(shift)
    real(dp), intent(out) :: shift(:,:,:)

    @:ASSERT(size(shift,dim=1) == size(shiftPerL_,dim=1))
    @:ASSERT(size(shift,dim=2) == size(shiftPerL_,dim=2))
    @:ASSERT(size(shift,dim=3) > 0)
    shift(:,:,:) = 0.0_dp
    shift(:,:,1) = shiftPerL_(:,:)

  end subroutine getShiftPerL



  !!* Returns the equivalency relations between the orbitals of the atoms.
  !!* @param orb Contains information about the atomic orbitals in the system
  !!* @param species Type of each atom (nAtom).
  !!* @param equiv The vector describing the equivalence on return.
  subroutine SCC_getOrbitalEquiv(orb, species, equiv)
    type(TOrbitals), intent(in) :: orb
    integer, intent(in) :: species(:)
    integer, intent(out) :: equiv(:,:,:)

    integer :: iAt, iOrb, iS, iSp
    integer :: nSpin, shift

    nSpin = size(equiv, dim=3)

    @:ASSERT(tInitialised_)
    @:ASSERT(size(species) == nAtom_)
    @:ASSERT(size(equiv, dim=1) == orb%mOrb)
    @:ASSERT(all(shape(equiv) == (/ orb%mOrb, nAtom_, nSpin /)))

    equiv(:,:,:) = 0
    shift = 0
    do iAt = 1, nAtom_
      iSp = species(iAt)
      do iOrb = 1, orb%nOrbSpecies(iSp)
        equiv(iOrb, iAt, 1) = iHubbU_(orb%iShellOrb(iOrb, iSp), iSp) + shift
      end do
      shift = shift + maxval(iHubbU_(:, iSp))
    end do
    do iS = 2, nSpin
      equiv(:,:,iS) = equiv(:,:,1)
    end do

  end subroutine SCC_getOrbitalEquiv


  !> Calculate the "double counting" force term using linearized XLBOMD form.
  !!
  !! Note: When SCC is driven in XLBOMD mode, the charges should NOT be updated
  !! after diagonalizing the Hamiltonian, so the charge stored in the module are
  !! the ingoing (auxiliary) charges, used to build the Hamiltonian.  However,
  !! the force double counting expession needs also the outgoing ones,
  !! therefore, there must be passed extra.
  !!
  subroutine addForceDCSCC_Xlbomd(species, orb, iNeighbor, img2CentCell, &
      & coord, qOrbitalOut, q0, force)
    integer,  intent(in) :: species(:)
    type(TOrbitals), intent(in) :: orb
    integer,  intent(in) :: iNeighbor(0:,:)
    integer,  intent(in) :: img2CentCell(:)
    real(dp), intent(in) :: coord(:,:)
    real(dp), intent(in) :: qOrbitalOut(:,:,:)
    real(dp), intent(in) :: q0(:,:,:)
    real(dp), intent(inout) :: force(:,:)

    real(dp), allocatable :: dQOut(:,:), dQOutAtom(:)
    real(dp), allocatable :: dQOutLShell(:,:), dQOutUniqU(:,:)

    allocate(dQOut(mOrb_, nAtom_))
    allocate(dQOutAtom(nAtom_))
    allocate(dQOutLShell(mShell_, nAtom_))
    allocate(dQOutUniqU(mHubbU_, nAtom_))

    call getNetCharges_(species, orb, qOrbitalOut, q0, dQOut, dQOutAtom, &
        & dQOutLShell, dQOutUniqU)

    ! Short-range part of gamma contribution
    call addGammaPrimeXlbomd_(deltaQUniqU_, dQOutUniqU, coord, species, &
        & iNeighbor, img2CentCell, force)

    ! 1/R contribution
    if (tPeriodic_) then
      call addInvRPrimeXlbomd(nAtom_, coord, nNeighEwald_, iNeighbor, &
          & img2CentCell, gLatPoint_, alpha_, volume_, deltaQAtom_, &
          & dQOutAtom, force)
      if (tExtChrg_) then
        call error("XLBOMD with external charges does not work yet!")
        !call addForceDCSCC_ExtChrg(force, chrgForce, coord, deltaQAtom_, &
        !    & gLatPoint_, alpha_, volume_)
      end if
    else
      call addInvRPrimeXlbomd(nAtom_, coord, deltaQAtom_, dQOutAtom, force)
      if (tExtChrg_) then
        call error("XLBOMD with external charges does not work yet!")
        !call addForceDCSCC_ExtChrg(force, chrgForce, coord, deltaQAtom_)
      end if
    end if

  end subroutine addForceDCSCC_Xlbomd


  !> Calculate  the derivative of the short range using the
  !! linearized XLBOMD formulation with auxiliary charges.
  !!
  subroutine addGammaPrimeXlbomd_(dQInUniqU, dQOutUniqU, coord, species, &
      & iNeighbor, img2CentCell, force)
    real(dp), intent(in) :: dQInUniqU(:,:), dQOutUniqU(:,:)
    real(dp), intent(in) :: coord(:,:)
    integer,  intent(in) :: species(:), iNeighbor(0:,:), img2CentCell(:)
    real(dp), intent(inout) :: force(:,:)

    integer :: iAt1, iAt2, iAt2f, iU1, iU2, iNeigh, iSp1, iSp2
    real(dp) :: rab, tmpGammaPrime, u1, u2, prefac
    real(dp) :: contrib(3)

    @:ASSERT(size(force,dim=1) == 3)
    @:ASSERT(size(force,dim=2) == nAtom_)

    do iAt1 = 1, nAtom_
      iSp1 = species(iAt1)
      do iNeigh = 1, maxval(nNeighShort_(:,:,:, iAt1))
        iAt2 = iNeighbor(iNeigh, iAt1)
        iAt2f = img2CentCell(iAt2)
        iSp2 = species(iAt2f)
        rab = sqrt(sum((coord(:,iAt1) - coord(:,iAt2))**2))
        do iU1 = 1, nHubbU_(species(iAt1))
          u1 = uniqHubbU_(iU1, iSp1)
          do iU2 = 1, nHubbU_(species(iAt2f))
            u2 = uniqHubbU_(iU2, species(iAt2f))
            if (iNeigh <= nNeighShort_(iU2,iU1,species(iAt2f),iAt1)) then
              if (tDampedShort_(iSp1) .or. tDampedShort_(iSp2)) then
                tmpGammaPrime = expGammaDampedPrime(rab, u2, u1, dampExp_)
              else
                tmpGammaPrime = expGammaPrime(rab, u2, u1)
              end if
              prefac = dQOutUniqU(iU1, iAt1) * dQInUniqU(iU2, iAt2f) &
                  & + dQInUniqU(iU1, iAt1) * dQOutUniqU(iU2, iAt2f) &
                  & - dQInUniqU(iU1, iAt1) * dQInUniqU(iU2, iAt2f)
              contrib(:) = prefac * tmpGammaPrime / rab  &
                  & * (coord(:,iAt2) - coord(:,iAt1))
              force(:,iAt1) = force(:,iAt1) + contrib
              force(:,iAt2f) = force(:,iAt2f) - contrib
            end if
          end do
        end do
      end do
    end do

  end subroutine addGammaPrimeXlbomd_


end module scc
