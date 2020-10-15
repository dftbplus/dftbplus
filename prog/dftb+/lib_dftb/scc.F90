!--------------------------------------------------------------------------------------------------!
!  DFTB+: general package for performing fast atomistic simulations                                !
!  Copyright (C) 2006 - 2020  DFTB+ developers group                                               !
!                                                                                                  !
!  See the LICENSE file for terms of usage and distribution.                                       !
!--------------------------------------------------------------------------------------------------!

#:include 'common.fypp'

!> Functions and local variables for the SCC calculation.
module dftbp_scc
  use dftbp_assert
  use dftbp_accuracy
  use dftbp_blasroutines
  use dftbp_charges, only : getSummedCharges
  use dftbp_chargeconstr
  use dftbp_commontypes
  use dftbp_constants
  use dftbp_coulomb, only : TCoulombCont, TCoulombInput, init, sumInvR
  use dftbp_dynneighlist
  use dftbp_environment
  use dftbp_fileid
  use dftbp_externalcharges
  use dftbp_h5correction
  use dftbp_message
#:if WITH_MPI
  use dftbp_mpifx
#:endif
  use dftbp_periodic
#:if WITH_SCALAPACK
  use dftbp_scalapackfx
#:endif
  use dftbp_shift
  use dftbp_shortgamma
  implicit none

  private

  public :: TSccInp, TScc, initialize


  !> Data necessary to initialize the SCC module
  type TSccInp

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

    !> if broadened external charges
    real(dp), allocatable :: blurWidths(:)

    !> any constraints on atomic charges
    real(dp), allocatable :: chrgConstraints(:,:)

    !> third order energy contributions
    real(dp), allocatable :: thirdOrderOn(:,:)

    !> H5 correction object
    type(TH5Corr), allocatable :: h5Correction

    !> Coulomb input
    type(TCoulombInput) :: coulombInput

    !> Whether shift vector is set externally -> skip internal shift calculation
    logical :: hasExternalShifts

  end type TSccInp


  !> private module variables for SCC
  type TScc
    private

    !> If module is initialised
    logical :: tInitialised = .false.

    !> Nr. of atoms
    integer :: nAtom

    !> Nr. of species
    integer :: nSpecies

    !> Max nr. of orbitals per atom
    integer :: mOrb

    !> Max nr. of shells per atom
    integer :: mShell

    !> Maximal number of unique U values
    integer :: mHubbU

    !> Shift vector per atom
    real(dp), allocatable :: shiftPerAtom(:)

    !> Shift vector per l-shell
    real(dp), allocatable :: shiftPerL(:,:)

    !> Short range interaction
    real(dp), allocatable :: shortGamma(:,:,:,:)

    !> CutOff for short range int.
    real(dp), allocatable :: shortCutOff(:,:,:,:)

    !> Maximal cutoff
    real(dp) :: cutoff

    !> Nr. of neighbours for short range part of gamma
    integer, allocatable :: nNeighShort(:,:,:,:)

    !> Atomic coordinates
    real(dp), allocatable :: coord(:,:)

    !> lattice vectors
    real(dp), allocatable :: latVecs(:,:)

    !> reciprocal lattice vectors
    real(dp), allocatable :: recVecs(:,:)

    !> Cell volume
    real(dp) :: volume

    !> Uniq Us per species
    real(dp), allocatable :: uniqHubbU(:,:)

    !> Nr. of uniq Us per species
    integer, allocatable :: nHubbU(:)

    !> Mapping L-shell -> uniq U
    integer, allocatable :: iHubbU(:,:)

    !> Negative gross charge
    real(dp), allocatable :: deltaQ(:,:)

    !> Negative gross charge per shell
    real(dp), allocatable :: deltaQPerLShell(:,:)

    !> Negative gross charge per atom
    real(dp), allocatable :: deltaQAtom(:)

    !> Negative gross charge per U
    real(dp), allocatable :: deltaQUniqU(:,:)

    !> Damped short range? (nSpecies)
    logical, allocatable :: tDampedShort(:)

    !> Damping exponent
    real(dp) :: dampExp

    !> Is the system periodic?
    logical :: tPeriodic

    !> Shifts due charge constrains?
    logical :: tChrgConstr

    !> Object for charge constraints
    type(TChrgConstr), allocatable :: chrgConstr

    !> use third order contributions
    logical :: tThirdOrder

    !> Shifts due to 3rd order
    type(TChrgConstr), allocatable :: thirdOrder

    !> External charges
    type(TExtCharge), allocatable :: extCharge

    !> Flag for activating the H5 H-bond correction
    logical :: tH5

    !> H5 correction object
    type(TH5Corr), allocatable :: h5Correction

    !> Coulombic interaction container
    type(TCoulombCont), public :: coulombCont

    !> Whether shift vector is set externally -> skip internal shift calculation
    logical :: hasExternalShifts

  contains

    !> Returns a minimal cutoff for the neighbourlist
    procedure :: getCutOff

    !> Returns the currenty used alpha parameter of the Ewald-summation
    procedure :: getEwaldPar

    !> Updates the atom coordinates for the SCC module
    procedure :: updateCoords

    !> Updates the SCC module, if the lattice vectors had been changed
    procedure :: updateLatVecs

    !> Updates the SCC module, if the charges have been changed
    procedure :: updateCharges

    !> Set external charge field
    procedure :: setExternalCharges

    !> Update potential shifts
    procedure :: updateShifts

    !> Routine for returning lower triangle of atomic resolved gamma as a matrix
    procedure :: getAtomicGammaMatrix

    !> Routine for returning lower triangle of atomic resolved gamma for specified U values
    procedure :: getAtomicGammaMatU

    !> Calculates the contribution of the SCC to the energy per atom
    procedure :: getEnergyPerAtom

    !> Calculates SCC energy contribution using the linearized XLBOMD form
    procedure :: getEnergyPerAtomXlbomd

    !> Calculates the contribution of the SCC to the forces
    procedure :: addForceDc

    !> Calculates the contribution of the stress tensor
    procedure :: addStressDc

    !> Returns the shift per atom coming from the SCC part
    procedure :: getShiftPerAtom

    !> Returns the shift per L contribution of the SCC.
    procedure :: getShiftPerL

    !> set electrostatic shifts (e.g. Poisson solver)
    procedure :: setShiftPerAtom

    !> set the shifts from outside (e.g. Poisson solver)
    procedure :: setShiftPerL

    !> Returns the equivalency relations between orbitals of the atoms
    procedure :: getOrbitalEquiv

    !> Calculate the "double counting" force term using linearized XLBOMD form
    procedure :: addForceDcXlbomd

    !> Returns potential from DFTB charges
    procedure :: getInternalElStatPotential

    !> Returns potential from external charges
    procedure :: getExternalElStatPotential

    !> Calculate gamma integral derivatives in SCC part
    procedure :: getGammaDeriv

    !> Get Q * inverse R contribution for the point charges
    procedure :: getShiftOfPC

  end type TScc


  !> Initialize SCC container from input data
  interface initialize
    module procedure Scc_initialize
  end interface initialize


contains


  !> Initialize SCC container from input data
  subroutine Scc_initialize(this, env, inp)

    !> Resulting instance
    type(TScc), intent(out) :: this

    !> Environment settings
    type(TEnvironment), intent(in) :: env

    !> Scc input
    type(TSccInp), intent(in) :: inp

    integer :: iSp1, iSp2, iU1, iU2, iL

    @:ASSERT(.not. this%tInitialised)

    this%nSpecies = size(inp%orb%nOrbSpecies)
    this%nAtom = size(inp%orb%nOrbAtom)
    this%mShell = inp%orb%mShell
    this%mOrb = inp%orb%mOrb

    @:ASSERT(allocated(inp%latVecs) .eqv. allocated(inp%recVecs))
    @:ASSERT(allocated(inp%latVecs) .eqv. (inp%volume > 0.0_dp))
    @:ASSERT(allocated(inp%extCharges) .or. .not. allocated(inp%blurWidths))

    if (allocated(inp%latVecs)) then
      call init(this%coulombCont, inp%coulombInput, env, this%nAtom, &
          & inp%latVecs, inp%recVecs, inp%volume)
    else
      call init(this%coulombCont, inp%coulombInput, env, this%nAtom)
    end if

    allocate(this%shiftPerAtom(this%nAtom))
    allocate(this%shiftPerL(this%mShell, this%nAtom))
    allocate(this%shortGamma(0, 0, 0, 0))
    this%tPeriodic = allocated(inp%latVecs)

    ! Initialize Hubbard U's
    allocate(this%uniqHubbU(this%mShell, this%nSpecies))
    allocate(this%nHubbU(this%nSpecies))
    allocate(this%iHubbU(this%mShell, this%nSpecies))
    this%iHubbU(:,:) = 0
    this%iHubbU(1,:) = 1
    this%nHubbU(:) = 1
    this%uniqHubbU(:,:) = 0.0_dp
    this%uniqHubbU(1,:) = inp%hubbU(1,:)
    do iSp1 = 1, this%nSpecies
      do iL = 2, inp%orb%nShell(iSp1)
        do iU1 = 1, this%nHubbU(iSp1)
          if (abs(inp%hubbU(iL,iSp1) - this%uniqHubbU(iU1,iSp1)) < MinHubDiff) then
            this%iHubbU(iL,iSp1) = iU1
            exit
          end if
        end do
        if (this%iHubbU(iL,iSp1) == 0) then
          this%nHubbU(iSp1) = this%nHubbU(iSp1) + 1
          this%uniqHubbU(this%nHubbU(iSp1),iSp1) = inp%hubbU(iL,iSp1)
          this%iHubbU(iL,iSp1) = this%nHubbU(iSp1)
        end if
      end do
    end do
    this%mHubbU = maxval(this%nHubbU)

    ! Get cutoff for short range coulomb
    allocate(this%shortCutOff(this%mHubbU, this%mHubbU, this%nSpecies, this%nSpecies))
    this%shortCutOff(:,:,:,:) = 0.0_dp
    do iSp1 = 1, this%nSpecies
      do iSp2 = iSp1, this%nSpecies
        do iU1 = 1, this%nHubbU(iSp1)
          do iU2 = 1, this%nHubbU(iSp2)
            this%shortCutOff(iU2, iU1, iSp2, iSp1) =&
                & expGammaCutOff(this%uniqHubbU(iU2, iSp2), this%uniqHubbU(iU1, iSp1))
            this%shortCutOff(iU1, iU2, iSp1, iSp2) = this%shortCutOff(iU2, iU1, iSp2, iSp1)
          end do
        end do
      end do
    end do
    this%cutoff = maxval(this%shortCutOff)

    if (this%tPeriodic) then
      this%latVecs = inp%latVecs
      this%recVecs = inp%recVecs
      this%volume = inp%volume
    end if

    ! Number of neighbours for short range cutoff and real part of Ewald
    allocate(this%nNeighShort(this%mHubbU, this%mHubbU, this%nSpecies, this%nAtom))

    ! Initialise external charges
    if (allocated(inp%extCharges)) then
      call this%setExternalCharges(inp%extCharges(1:3,:), inp%extCharges(4,:),&
          & blurWidths=inp%blurWidths)
    end if

    this%tChrgConstr = allocated(inp%chrgConstraints)
    if (this%tChrgConstr) then
      allocate(this%chrgConstr)
      call init(this%chrgConstr, inp%chrgConstraints, 2)
    end if
    this%tThirdOrder = allocated(inp%thirdOrderOn)
    if (this%tThirdOrder) then
      allocate(this%thirdOrder)
      ! Factor 1/6 in the energy is put into the Hubbard derivatives
      call init(this%thirdOrder, inp%thirdOrderOn / 6.0_dp, 3)
    end if

    ! Initialise arrays for charge differences
    allocate(this%deltaQ(this%mOrb, this%nAtom))
    allocate(this%deltaQPerLShell(this%mShell, this%nAtom))
    allocate(this%deltaQAtom(this%nAtom))
    allocate(this%deltaQUniqU(this%mHubbU, this%nAtom))

    ! Initialise short range damping
    allocate(this%tDampedShort(this%nSpecies))
    this%tDampedShort(:) = inp%tDampedShort(:)
    this%dampExp = inp%dampExp

    ! H5 correction
    this%tH5 = allocated(inp%h5Correction)
    if (this%tH5) then
      this%h5Correction = inp%h5Correction
    end if

    this%hasExternalShifts = inp%hasExternalShifts

    this%tInitialised = .true.

  end subroutine Scc_initialize


  !> Returns a minimal cutoff for the neighbourlist, which must be passed to various functions in
  !> this module.
  function getCutOff(this) result(cutoff)

    !> Instance
    class(TScc), intent(in) :: this

    !> The neighbourlists, passed to scc routines, should contain neighbour information at
    !> least up to that cutoff.
    real(dp) :: cutoff

    @:ASSERT(this%tInitialised)
    cutoff = this%cutoff

  end function getCutOff


  !> Returns the currenty used alpha parameter of the Ewald-summation
  function getEwaldPar(this) result(alpha)

    !> Instance
    class(TScc), intent(in) :: this

    !> Parameter in the Ewald summation.
    real(dp) :: alpha

    @:ASSERT(this%tInitialised)
    alpha = this%coulombCont%alpha

  end function getEwaldPar


  !> Updates the atom coordinates for the SCC module.
  subroutine updateCoords(this, env, coord, species, neighList)

    !> Instance
    class(TScc), intent(inout) :: this

    !> Environment settings
    type(TEnvironment), intent(in) :: env

    !> New coordinates of the atoms
    real(dp), intent(in) :: coord(:,:)

    !> Species of the atoms (should not change during run)
    integer, intent(in) :: species(:)

    !> Neighbour list for the atoms.
    type(TNeighbourList), intent(in) :: neighList

    @:ASSERT(this%tInitialised)

    this%coord = coord

    call this%coulombCont%updateCoords(env, neighList, coord, species)

    call updateNNeigh_(this, species, neighList)

    call initGamma_(this, species, neighList%iNeighbour)

    if (allocated(this%extCharge)) then
      if (this%tPeriodic) then
        call this%extCharge%setCoordinates(env, this%coord, this%coulombCont%rCellVec,&
            & this%coulombCont%gLatPoint, this%coulombCont%alpha, this%volume)
      else
        call this%extCharge%setCoordinates(env, this%coord)
      end if
    end if

  end subroutine updateCoords


  !> Updates the SCC module, if the lattice vectors had been changed
  subroutine updateLatVecs(this, latVec, recVec, vol)

    !> Instance
    class(TScc), intent(inout) :: this

    !> New lattice vectors
    real(dp), intent(in) :: latVec(:,:)

    !> New reciprocal lattice vectors
    real(dp), intent(in) :: recVec(:,:)

    !> New volume
    real(dp), intent(in) :: vol

    @:ASSERT(this%tInitialised)
    @:ASSERT(this%tPeriodic)

    call this%coulombCont%updateLatVecs(latVec, recVec, vol)

    this%volume = vol
    this%latVecs(:,:) = latVec
    this%recVecs(:,:) = recVec

    if (allocated(this%extCharge)) then
      call this%extCharge%setLatticeVectors(latVec, recVec)
    end if

  end subroutine updateLatVecs


  !> Updates the SCC module, if the charges have been changed
  subroutine updateCharges(this, env, qOrbital, q0, orb, species)

    !> Resulting module variables
    class(TScc), intent(inout) :: this

    !> Environment settings
    type(TEnvironment), intent(in) :: env

    !> Orbital resolved charges
    real(dp), intent(in) :: qOrbital(:,:,:)

    !> Reference charge distribution (neutral atoms)
    real(dp), intent(in) :: q0(:,:,:)

    !> Contains information about the atomic orbitals in the system
    type(TOrbitals), intent(in) :: orb

    !> Species of the atoms (should not change during run)
    integer, intent(in) :: species(:)

    @:ASSERT(this%tInitialised)

    call getSummedCharges(species, orb, qOrbital, q0, iHubbU=this%iHubbU, dQ=this%deltaQ, &
        & dQAtom=this%deltaQAtom, dQShell=this%deltaQPerLShell, dQUniqU=this%deltaQUniqU)

    call this%coulombCont%updateCharges(env, qOrbital, q0, orb, species, &
        & this%deltaQ, this%deltaQAtom, this%deltaQPerLShell, this%deltaQUniqU)

  end subroutine updateCharges


  !> Update potential shifts. Call after updateCharges
  subroutine updateShifts(this, env, orb, species, iNeighbour, img2CentCell)

    !> Resulting module variables
    class(TScc), intent(inout) :: this

    !> Environment settings
    type(TEnvironment), intent(in) :: env

    !> Contains information about the atomic orbitals in the system
    type(TOrbitals), intent(in) :: orb

    !> Species of the atoms (should not change during run)
    integer, intent(in) :: species(:)

    !> Neighbour indexes
    integer, intent(in) :: iNeighbour(0:,:)

    !> Mapping on atoms in the central cell
    integer, intent(in) :: img2CentCell(:)

    @:ASSERT(this%tInitialised)
    @:ASSERT(.not. this%hasExternalShifts)

    call buildShifts_(this, env, orb, species, iNeighbour, img2CentCell)

    call this%coulombCont%updateShifts(env, orb, species, iNeighbour, img2CentCell)

    if (this%tChrgConstr) then
      call buildShift(this%chrgConstr, this%deltaQAtom)
    end if
    if (this%tThirdOrder) then
      call buildShift(this%thirdOrder, this%deltaQAtom)
    end if

  end subroutine updateShifts


  !> set external charge field
  subroutine setExternalCharges(this, chargeCoords, chargeQs, blurWidths)

    !> Instance
    class(TScc), intent(inout) :: this

    !> Coordinates of external charges
    real(dp), intent(in) :: chargeCoords(:,:)

    !> Magnitude of external charges
    real(dp), intent(in) :: chargeQs(:)

    !> Spatial extension of external charge distribution
    real(dp), intent(in), optional :: blurWidths(:)

    if (.not. allocated(this%extCharge)) then
      allocate(this%extCharge)
      call TExtCharge_init(this%extCharge, this%nAtom, size(chargeQs), this%tPeriodic)
    end if

    if (present(blurWidths)) then
      if (this%tPeriodic .and. any(blurWidths > 1.0e-7_dp)) then
        if (1.0_dp / maxval(blurWidths) < this%coulombCont%alpha) then
          call error("Charge blur widths are too wide compared to the Ewald real space sum")
        end if
      end if
      call this%extCharge%setExternalCharges(chargeCoords, chargeQs, blurWidths=blurWidths)
    else
      call this%extCharge%setExternalCharges(chargeCoords, chargeQs)
    end if

  end subroutine setExternalCharges


  !> Routine for returning lower triangle of atomic resolved gamma as a matrix
  subroutine getAtomicGammaMatrix(this, gammamat, iNeighbour, img2CentCell)

    !> Instance
    class(TScc), intent(in) :: this

    !> Atom resolved gamma
    real(dp), intent(out) :: gammamat(:,:)

    !> neighbours of atoms
    integer, intent(in) :: iNeighbour(0:,:)

    !> index array between images and central cell
    integer, intent(in) :: img2CentCell(:)

    integer :: iAt1, iAt2, iAt2f, iNeigh

    @:ASSERT(this%tInitialised)
    @:ASSERT(all(shape(gammamat) == [ this%nAtom, this%nAtom ]))
    @:ASSERT(all(this%nHubbU == 1))

  #:if WITH_SCALAPACK
    call error("scc:getAtomicGammaMatrix does not work with MPI yet")
  #:endif
    gammamat(:,:) = this%coulombCont%invRMat
    do iAt1 = 1, this%nAtom
      do iNeigh = 0, maxval(this%nNeighShort(:,:,:, iAt1))
        iAt2 = iNeighbour(iNeigh, iAt1)
        iAt2f = img2CentCell(iAt2)
        gammamat(iAt2f, iAt1) = gammamat(iAt2f, iAt1) - this%shortGamma(1, 1, iNeigh, iAt1)
      end do
    end do

  end subroutine getAtomicGammaMatrix


  !> Routine for returning lower triangle of atomic resolved Coulomb matrix
  subroutine getAtomicGammaMatU(this, gammamat, U_h, species, iNeighbour, img2CentCell)

    !> Instance
    class(TScc), intent(in) :: this

    !> Atom resolved gamma
    real(dp), intent(out) :: gammamat(:,:)

    !> ppRPA Hubbard parameters
    real(dp), intent(in) :: U_h(:)

    !> List of the species for each atom.
    integer,  intent(in) :: species(:)

    !> neighbours of atoms
    integer, intent(in) :: iNeighbour(0:,:)

    !> index array between images and central cell
    integer, intent(in) :: img2CentCell(:)

    integer  :: iAt1, iAt2, iSp1, iSp2, iAt2f, iNeigh
    real(dp) :: R_ab

    @:ASSERT(this%tInitialised)
    @:ASSERT(all(shape(gammamat) == [ this%nAtom, this%nAtom ]))
    @:ASSERT(all(this%nHubbU == 1))

  #:if WITH_SCALAPACK
    call error("scc:getAtomicGammaMatU does not work with MPI yet")
  #:endif
    gammamat(:,:) = this%coulombCont%invRMat
    do iAt1 = 1, this%nAtom
      iSp1 = species(iAt1)
      do iNeigh = 0, maxval(this%nNeighShort(:,:,:, iAt1))
        iAt2 = iNeighbour(iNeigh, iAt1)
        iAt2f = img2CentCell(iAt2)
        iSp2 = species(iAt2f)
        R_ab = sqrt(sum((this%coord(:,iAt1) - this%coord(:,iAt2))**2))
        gammamat(iAt2f, iAt1) = gammamat(iAt2f, iAt1) - expGamma(R_ab, U_h(iSp2), U_h(iSp1))
      end do
    end do

    do iAt1 = 1, this%nAtom
      do iAt2 = 1, iAt1 - 1
        gammamat(iAt2, iAt1) = gammamat(iAt1, iAt2)
      end do
    end do

  end subroutine getAtomicGammaMatU


  !> Calculates the contribution of the charge consistent part to the energy per atom.
  subroutine getEnergyPerAtom(this, eScc)

    !> Resulting module variables
    class(TScc), intent(in) :: this

    !> The SCC contribution to the energy
    real(dp), intent(out) :: eScc(:)

    @:ASSERT(this%tInitialised)
    @:ASSERT(size(eScc) == this%nAtom)

    eScc(:) = 0.5_dp * (this%shiftPerAtom * this%deltaQAtom&
        & + sum(this%shiftPerL * this%deltaQPerLShell, dim=1))

    if (.not. this%hasExternalShifts) then
      call this%coulombCont%addEnergy(eScc)
    end if

    if (allocated(this%extCharge)) then
      call this%extCharge%addEnergyPerAtom(this%deltaQAtom, eScc)
    end if
    if (this%tChrgConstr) then
      call addEnergyPerAtom(this%chrgConstr, eScc, this%deltaQAtom)
    end if
    if (this%tThirdOrder) then
      call addEnergyPerAtom(this%thirdOrder, eScc, this%deltaQAtom)
    end if

  end subroutine getEnergyPerAtom


  !> Calculates SCC energy contribution using the linearized XLBOMD form.
  !> Note: When SCC is driven in XLBOMD mode, the charges should NOT be updated after diagonalizing
  !> the Hamiltonian, so the charge stored in the module are the input (auxiliary) charges, used to
  !> build the Hamiltonian.  However, the linearized energy expession needs also the output charges,
  !> therefore these are passed in as an extra variable.
  subroutine getEnergyPerAtomXlbomd(this, species, orb, qOut, q0, eScc)

    !> Resulting module variables
    class(TScc), intent(in) :: this

    !> atomic species
    integer, intent(in) :: species(:)

    !> orbital information
    type(TOrbitals), intent(in) :: orb

    !> output charges
    real(dp), intent(in) :: qOut(:,:,:)

    !> reference charges
    real(dp), intent(in) :: q0(:,:,:)

    !> energy contributions
    real(dp), intent(out) :: eScc(:)

    real(dp), allocatable :: dQOut(:,:), dQOutAtom(:), dQOutShell(:,:)

    @:ASSERT(this%tInitialised)
    @:ASSERT(size(eScc) == this%nAtom)

    allocate(dQOut(orb%mOrb, this%nAtom))
    allocate(dQOutAtom(this%nAtom))
    allocate(dQOutShell(this%mShell, this%nAtom))

    call getSummedCharges(species, orb, qOut, q0, dQ=dQOut, dQAtom=dQOutAtom, dQShell=dQOutShell)

    ! 1/2 sum_A (2 q_A - n_A) * shift(n_A)
    eScc(:) = 0.5_dp * (this%shiftPerAtom * (2.0_dp * dQOutAtom - this%deltaQAtom)&
        & + sum(this%shiftPerL * (2.0_dp * dQOutShell - this%deltaQPerLShell), dim=1))

    if (.not. this%hasExternalShifts) then
      call this%coulombCont%addEnergy(eScc, dQOut, dQOutAtom, dQOutShell)
    end if

    if (allocated(this%extCharge)) then
      call error("XLBOMD not working with external charges yet")
      !call addEnergyPerAtom_ExtChrg(this%deltaQAtom, eScc)
    end if
    if (this%tChrgConstr) then
      call error("XLBOMD not working with charge constraints yet")
      !call addEnergyPerAtom(this%chrgConstr, eScc, this%deltaQAtom)
    end if
    if (this%tThirdOrder) then
      call error("XLBOMD not working with third order yet")
      !call addEnergyPerAtom(this%thirdOrder, eScc, this%deltaQAtom)
    end if

  end subroutine getEnergyPerAtomXlbomd


  !> Calculates the contribution of the charge consistent part to the forces for molecules/clusters,
  !> which is not covered in the term with the shift vectors.
  subroutine addForceDc(this, env, force, species, iNeighbour, img2CentCell, chrgForce)

    !> Resulting module variables
    class(TScc), intent(in) :: this

    !> Environment settings
    type(TEnvironment), intent(in) :: env

    !> has force contribution added
    real(dp), intent(inout) :: force(:,:)

    !> Species for each atom.
    integer, intent(in) :: species(:)

    !> List of neighbours for each atom.
    integer, intent(in) :: iNeighbour(0:,:)

    !> Indexing of images of the atoms in the central cell.
    integer, intent(in) :: img2CentCell(:)

    !> Force contribution due to the external charges, which is not contained in the term with the
    !> shift vectors.
    real(dp), intent(inout), optional :: chrgForce(:,:)

    @:ASSERT(this%tInitialised)
    @:ASSERT(size(force,dim=1) == 3)
    @:ASSERT(size(force,dim=2) == this%nAtom)
    @:ASSERT(present(chrgForce) .eqv. allocated(this%extCharge))

    ! Short-range part of gamma contribution
    call addGammaPrime_(this, force, species, iNeighbour, img2CentCell)

    call this%coulombCont%addGradients(env, this%coord, species, iNeighbour, &
        & img2CentCell, force)

    if (allocated(this%extCharge)) then
      if (this%tPeriodic) then
        call this%extCharge%addForceDc(env, force, chrgForce, this%coord, this%deltaQAtom,&
            & this%coulombCont%rCellVec, this%coulombCont%gLatPoint, this%coulombCont%alpha,&
            & this%volume)
      else
        call this%extCharge%addForceDc(env, force, chrgForce, this%coord, this%deltaQAtom)
      end if
    end if

  end subroutine addForceDc


  !> Calculates the contribution of the stress tensor which is not covered in the term with the
  !> shift vectors.
  subroutine addStressDc(this, st, env, species, iNeighbour, img2CentCell)

    !> Resulting module variables
    class(TScc), intent(in) :: this

    !> Add stress tensor contribution to this
    real(dp), intent(inout) :: st(:,:)

    !> Computational environment settings
    type(TEnvironment), intent(in) :: env

    !> Species for each atom.
    integer, intent(in) :: species(:)

    !> List of neighbours for each atom.
    integer, intent(in) :: iNeighbour(0:,:)

    !> Indexing of images of the atoms in the central cell.
    integer, intent(in) :: img2CentCell(:)

    real(dp) :: stTmp(3,3)

    @:ASSERT(this%tInitialised)
    @:ASSERT(this%tPeriodic)
    @:ASSERT(all(shape(st)==(/3,3/)))

    stTmp(:,:) = 0.0_dp

    ! Short-range part of gamma contribution
    call addSTGammaPrime_(stTmp,this,species,iNeighbour,img2CentCell)

    st(:,:) = st(:,:) - 0.5_dp * stTmp(:,:)

    ! 1/R contribution
    call this%coulombCont%addStress(env, this%coord, species, iNeighbour, &
        & img2CentCell, st)

    ! if (tExtChrg_) then
    ! ????
    ! end if

  end subroutine addStressDc


  !> Returns the shift per atom coming from the SCC part
  subroutine getShiftPerAtom(this, shift)

    !> Instance
    class(TScc), intent(in) :: this

    !> Contains the shift on exit.
    real(dp), intent(out) :: shift(:)

    @:ASSERT(this%tInitialised)
    @:ASSERT(size(shift) == size(this%shiftPerAtom))

    shift(:) = this%shiftPerAtom
    if (.not. this%hasExternalShifts) then
      call this%coulombCont%addShiftPerAtom(shift)
    end if
    if (allocated(this%extCharge)) then
      call this%extCharge%addShiftPerAtom(shift)
    end if
    if (this%tChrgConstr) then
      call addShiftPerAtom(this%chrgConstr, shift)
    end if
    if (this%tThirdOrder) then
      call addShiftPerAtom(this%thirdOrder, shift)
    end if

  end subroutine getShiftPerAtom


  !> Returns the shift per L contribution of the SCC.
  subroutine getShiftPerL(this, shift)

    !> Instance
    class(TScc), intent(in) :: this

    !> Contains the shift on exit.
    real(dp), intent(out) :: shift(:,:)


    @:ASSERT(this%tInitialised)
    @:ASSERT(size(shift,dim=1) == size(this%shiftPerL,dim=1))
    @:ASSERT(size(shift,dim=2) == size(this%shiftPerL,dim=2))

    shift(:, :) = this%shiftPerL
    if (.not. this%hasExternalShifts) then
      call this%coulombCont%addShiftPerShell(shift)
    end if

  end subroutine getShiftPerL


  !> set electrostatic shifts (e.g. Poisson solver)
  subroutine setShiftPerAtom(this, shift)

    !> Instance
    class(TScc), intent(inout) :: this

    !> Contains the input shifts: shift(atom).
    real(dp), intent(in) :: shift(:)

    @:ASSERT(this%tInitialised)
    @:ASSERT(size(shift) == size(this%shiftPerAtom,dim=1))
    @:ASSERT(this%hasExternalShifts)

    this%shiftPerAtom(:) = shift

  end subroutine setShiftPerAtom


  !> set the shifts from outside (e.g. Poisson solver)
  subroutine setShiftPerL(this, shift)

    !> Instance
    class(TScc), intent(inout) :: this

    !> Contains the input shifts (shell, Atom)
    real(dp), intent(in) :: shift(:,:)

    @:ASSERT(this%tInitialised)
    @:ASSERT(size(shift,dim=1) == size(this%shiftPerL,dim=1))
    @:ASSERT(size(shift,dim=2) == size(this%shiftPerL,dim=2))
    @:ASSERT(this%hasExternalShifts)

    this%shiftPerL(:, :) = shift

  end subroutine setShiftPerL

  !> Returns the equivalency relations between orbitals of the atoms. If transfering charge between
  !> the orbitals does not change the electrostatic energy, they are considered equivalent.
  subroutine getOrbitalEquiv(this, orb, species, equiv)

    !> Resulting module variables
    class(TScc), intent(in) :: this

    !> Contains information about the atomic orbitals in the system
    type(TOrbitals), intent(in) :: orb

    !> Type of each atom (nAtom).
    integer, intent(in) :: species(:)

    !> The vector describing the equivalence on return.
    integer, intent(out) :: equiv(:,:,:)

    integer :: iAt, iOrb, iS, iSp
    integer :: nSpin, shift

    nSpin = size(equiv, dim=3)

    @:ASSERT(this%tInitialised)
    @:ASSERT(size(species) == this%nAtom)
    @:ASSERT(size(equiv, dim=1) == orb%mOrb)
    @:ASSERT(all(shape(equiv) == (/ orb%mOrb, this%nAtom, nSpin /)))

    equiv(:,:,:) = 0
    shift = 0
    do iAt = 1, this%nAtom
      iSp = species(iAt)
      do iOrb = 1, orb%nOrbSpecies(iSp)
        equiv(iOrb, iAt, 1) = this%iHubbU(orb%iShellOrb(iOrb, iSp), iSp) + shift
      end do
      shift = shift + maxval(this%iHubbU(:, iSp))
    end do
    do iS = 2, nSpin
      equiv(:,:,iS) = equiv(:,:,1)
    end do

  end subroutine getOrbitalEquiv


  !> Calculate the "double counting" force term using linearized XLBOMD form.
  !> Note: When SCC is driven in XLBOMD mode, the charges should NOT be updated after diagonalizing
  !> the Hamiltonian, so the charge stored in the module are the input (auxiliary) charges, used to
  !> build the Hamiltonian.  However, the linearized energy expession needs also the output charges,
  !> therefore these are passed in as an extra variable.
  subroutine addForceDcXlbomd(this, env, species, orb, iNeighbour, img2CentCell, qOrbitalOut,&
      & q0, force)

    !> Resulting module variables
    class(TScc), intent(in) :: this

    !> Environment settings
    type(TEnvironment), intent(in) :: env

    !> atomic species
    integer, intent(in) :: species(:)

    !> orbital information
    type(TOrbitals), intent(in) :: orb

    !> neighbours surrounding each atom
    integer, intent(in) :: iNeighbour(0:,:)

    !> index from image atoms to central cell
    integer, intent(in) :: img2CentCell(:)

    !> output charges
    real(dp), intent(in) :: qOrbitalOut(:,:,:)

    !> reference charges
    real(dp), intent(in) :: q0(:,:,:)

    !> Force terms are added to this
    real(dp), intent(inout) :: force(:,:)

    real(dp), allocatable :: dQOut(:,:), dQOutAtom(:)
    real(dp), allocatable :: dQOutLShell(:,:), dQOutUniqU(:,:)

    allocate(dQOut(this%mOrb, this%nAtom))
    allocate(dQOutAtom(this%nAtom))
    allocate(dQOutLShell(this%mShell, this%nAtom))
    allocate(dQOutUniqU(this%mHubbU, this%nAtom))

    call getSummedCharges(species, orb, qOrbitalOut, q0, iHubbU=this%iHubbU, dQ=dQOut,&
        & dQAtom=dQOutAtom, dQShell=dQOutLShell, dQUniqU=dQOutUniqU)

    ! Short-range part of gamma contribution
    call addGammaPrimeXlbomd_(this, this%deltaQUniqU, dQOutUniqU, species, iNeighbour,&
        & img2CentCell, force)

    call this%coulombCont%addGradients(env, this%coord, species, iNeighbour, &
        & img2CentCell, force, dQOut, dQOutAtom, dQOutLShell)

    if (allocated(this%extCharge)) then
      call error("XLBOMD with external charges does not work yet!")
    end if

  end subroutine addForceDcXlbomd


  !> Returns potential from DFTB charges
  subroutine getInternalElStatPotential(this, V, env, locations, epsSoften)

    !> Instance of SCC calculation
    class(TScc), intent(in) :: this

    !> Resulting potentials
    real(dp), intent(out) :: V(:)

    !> Computational environment settings
    type(TEnvironment), intent(in) :: env

    !> sites to calculate potential
    real(dp), intent(in) :: locations(:,:)

    !> optional potential softening
    real(dp), optional, intent(in) :: epsSoften

    @:ASSERT(this%tInitialised)
    @:ASSERT(all(shape(locations) == [3,size(V)]))

    V(:) = 0.0_dp

    if (this%tPeriodic) then
      call sumInvR(env, size(V), this%nAtom, locations, this%coord, this%deltaQAtom,&
          & this%coulombCont%rCellVec, this%coulombCont%gLatPoint, this%coulombCont%alpha,&
          & this%volume, V, epsSoften=epsSoften)
    else
      call sumInvR(env, size(V), this%nAtom, locations, this%coord, this%deltaQAtom, V,&
          & epsSoften=epsSoften)
    end if

  end subroutine getInternalElStatPotential


  !> Returns potential from external charges
  subroutine getExternalElStatPotential(this, V, env, locations, epsSoften)

    !> Instance
    class(TScc), intent(in) :: this

    !> Resulting potentials
    real(dp), intent(out) :: V(:)

    !> Computational environment settings
    type(TEnvironment), intent(in) :: env

    !> sites to calculate potential
    real(dp), intent(in) :: locations(:,:)

    !> optional potential softening
    real(dp), optional, intent(in) :: epsSoften

    @:ASSERT(this%tInitialised)

    if (allocated(this%extCharge)) then
      if (this%tPeriodic) then
        call this%extCharge%getElStatPotential(env, locations, this%coulombCont%rCellVec,&
            & this%coulombCont%gLatPoint, this%coulombCont%alpha, this%volume, V,&
            & epsSoften=epsSoften)
      else
        call this%extCharge%getElStatPotential(env, locations, V, epsSoften=epsSoften)
      end if
    else
      V(:) = 0.0_dp
    end if

  end subroutine getExternalElStatPotential


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!! Private routines
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Updates the number of neighbours for the SCC module (local).
  subroutine updateNNeigh_(this, species, neighList)

    !> Instance
    type(TScc), intent(inout) :: this

    !> Species for each atom
    integer, intent(in) :: species(:)

    !> Neighbour list for the atoms in the system.
    type(TNeighbourList), intent(in) :: neighList


    integer :: iAt1, iSp2, iU1, iU2

    this%nNeighShort(:,:,:,:) = 0
    do iAt1 = 1, this%nAtom
      do iSp2 = 1, this%nSpecies
        do iU1 = 1, this%nHubbU(species(iAt1))
          do iU2 = 1, this%nHubbU(iSp2)
            this%nNeighShort(iU2, iU1, iSp2, iAt1) =&
                & getNrOfNeighbours(neighList, this%shortCutOff(iU2, iU1, iSp2, species(iAt1)),&
                & iAt1)
          end do
        end do
      end do
    end do

  end subroutine updateNNeigh_


  !> Constructs the shift vectors for the SCC contributions.
  !> The full shift vector must be constructed by adding shiftAtom and shiftShell accordingly.
  subroutine buildShifts_(this, env, orb, species, iNeighbour, img2CentCell)

    !> Resulting module variables
    type(TScc), intent(inout) :: this

    !> Environment settings
    type(TEnvironment), intent(in) :: env

    !> Contains information about the atomic orbitals in the system
    type(TOrbitals), intent(in) :: orb

    !> List of the species for each atom.
    integer, intent(in) :: species(:)

    !> List of surrounding neighbours for each atom.
    integer, intent(in) :: iNeighbour(0:,:)

    !> Image of each atom in the central cell.
    integer, intent(in) :: img2CentCell(:)

    @:ASSERT(this%tInitialised)

    this%shiftPerAtom(:) = 0.0_dp

    call buildShiftPerShell_(this, env, orb, species, iNeighbour, img2CentCell)

  end subroutine buildShifts_


  !> Builds the short range shell resolved part of the shift vector
  subroutine buildShiftPerShell_(this, env, orb, species, iNeighbour, img2CentCell)

    !> Resulting module variables
    type(TScc), intent(inout), target :: this

    !> Environment settings
    type(TEnvironment), intent(in) :: env

    !> Contains information about the atomic orbitals in the system
    type(TOrbitals), intent(in) :: orb

    !> List of the species for each atom.
    integer, intent(in) :: species(:)

    !> List of surrounding neighbours for each atom.
    integer, intent(in) :: iNeighbour(0:,:)

    !> Image of each atom in the central cell.
    integer, intent(in) :: img2CentCell(:)

    integer :: iAt1, iAt2f, iSp1, iSp2, iSh1, iSh2, iU1, iU2, iNeigh, iAt1Start, iAt1End


  #:if WITH_SCALAPACK
    if (env%blacs%atomGrid%iProc /= -1) then
      iAt1Start = env%blacs%atomGrid%iproc * this%nAtom / env%blacs%atomGrid%nproc + 1
      iAt1End = (env%blacs%atomGrid%iproc + 1) * this%nAtom / env%blacs%atomGrid%nproc
    else
      ! Do not calculate anything if process is not part of the atomic grid
      iAt1Start = 0
      iAt1End = -1
    end if
  #:else
    iAt1Start = 1
    iAt1End = this%nAtom
  #:endif

    this%shiftPerL(:,:) = 0.0_dp
    do iAt1 = iAt1Start, iAt1End
      iSp1 = species(iAt1)
      do iSh1 = 1, orb%nShell(iSp1)
        iU1 = this%iHubbU(iSh1, iSp1)
        do iNeigh = 0, maxval(this%nNeighShort(:,:,:,iAt1))
          iAt2f = img2CentCell(iNeighbour(iNeigh, iAt1))
          iSp2 = species(iAt2f)
          do iSh2 = 1, orb%nShell(iSp2)
            iU2 = this%iHubbU(iSh2, iSp2)
            this%shiftPerL(iSh1, iAt1) = this%shiftPerL(iSh1, iAt1)&
                & - this%shortGamma(iU2, iU1, iNeigh, iAt1) * this%deltaQPerLShell(iSh2, iAt2f)
            if (iAt2f /= iAt1) then
              this%shiftPerL(iSh2, iAt2f) = this%shiftPerL(iSh2, iAt2f)&
                  & - this%shortGamma(iU2, iU1, iNeigh, iAt1) * this%deltaQPerLShell(iSh1, iAt1)
            end if
          end do
        end do
      end do
    end do

  #:if WITH_SCALAPACK
    call mpifx_allreduceip(env%mpi%groupComm, this%shiftPerL, MPI_SUM)
  #:endif

  end subroutine buildShiftPerShell_


  !> Calculate the derivative of the short range contributions using the linearized XLBOMD
  !> formulation with auxiliary charges.
  subroutine addGammaPrimeXlbomd_(this, dQInUniqU, dQOutUniqU, species, iNeighbour,&
      & img2CentCell, force)

    !> Resulting module variables
    type(TScc), intent(in) :: this

    !> Input charges
    real(dp), intent(in) :: dQInUniqU(:,:)

    !> output charges
    real(dp), intent(in) :: dQOutUniqU(:,:)

    !> chemical species
    integer, intent(in) :: species(:)

    !> neighbours around atoms
    integer, intent(in) :: iNeighbour(0:,:)

    !> image to real atom indexing
    integer, intent(in) :: img2CentCell(:)

    !> term to add force contributions to
    real(dp), intent(inout) :: force(:,:)

    integer :: iAt1, iAt2, iAt2f, iU1, iU2, iNeigh, iSp1, iSp2
    real(dp) :: rab, tmpGammaPrime, u1, u2, prefac
    real(dp) :: contrib(3)

    @:ASSERT(this%tInitialised)
    @:ASSERT(size(force,dim=1) == 3)
    @:ASSERT(size(force,dim=2) == this%nAtom)

    do iAt1 = 1, this%nAtom
      iSp1 = species(iAt1)
      do iNeigh = 1, maxval(this%nNeighShort(:,:,:, iAt1))
        iAt2 = iNeighbour(iNeigh, iAt1)
        iAt2f = img2CentCell(iAt2)
        iSp2 = species(iAt2f)
        rab = sqrt(sum((this%coord(:,iAt1) - this%coord(:,iAt2))**2))
        do iU1 = 1, this%nHubbU(species(iAt1))
          u1 = this%uniqHubbU(iU1, iSp1)
          do iU2 = 1, this%nHubbU(species(iAt2f))
            u2 = this%uniqHubbU(iU2, species(iAt2f))
            if (iNeigh <= this%nNeighShort(iU2,iU1,species(iAt2f),iAt1)) then
              if (this%tDampedShort(iSp1) .or. this%tDampedShort(iSp2)) then
                tmpGammaPrime = expGammaDampedPrime(rab, u2, u1, this%dampExp)
              else
                tmpGammaPrime = expGammaPrime(rab, u2, u1)
              end if
              prefac = dQOutUniqU(iU1, iAt1) * dQInUniqU(iU2, iAt2f)&
                  & + dQInUniqU(iU1, iAt1) * dQOutUniqU(iU2, iAt2f)&
                  & - dQInUniqU(iU1, iAt1) * dQInUniqU(iU2, iAt2f)
              contrib(:) = prefac * tmpGammaPrime / rab  * (this%coord(:,iAt2) - this%coord(:,iAt1))
              force(:,iAt1) = force(:,iAt1) + contrib
              force(:,iAt2f) = force(:,iAt2f) - contrib
            end if
          end do
        end do
      end do
    end do

  end subroutine addGammaPrimeXlbomd_


  !> Set up the storage and internal values for the short range part of Gamma.
  subroutine initGamma_(this, species, iNeighbour)

    !> Resulting module variables
    type(TScc), intent(inout) :: this

    !> List of the species for each atom.
    integer, intent(in) :: species(:)

    !> Index of neighbouring atoms for each atom.
    integer, intent(in) :: iNeighbour(0:,:)

    integer :: iAt1, iAt2, iU1, iU2, iNeigh, iSp1, iSp2
    real(dp) :: rab, u1, u2

    @:ASSERT(this%tInitialised)

    ! Reallocate shortGamma, if it does not contain enough neighbours
    if (size(this%shortGamma, dim=3) < maxval(this%nNeighShort)+1) then
      deallocate(this%shortGamma)
      allocate(this%shortGamma(this%mHubbU, this%mHubbU, 0:maxval(this%nNeighShort),&
          & this%nAtom))
    end if
    this%shortGamma(:,:,:,:) = 0.0_dp

    ! some additional symmetry not used, as the value of gamma for atoms
    ! interacting with themselves is the same for all atoms of the same species
    do iAt1 = 1, this%nAtom
      iSp1 = species(iAt1)
      do iNeigh = 0, maxval(this%nNeighShort(:,:,:, iAt1))
        iAt2 = iNeighbour(iNeigh, iAt1)
        iSp2 = species(iAt2)
        rab = sqrt(sum((this%coord(:,iAt1) - this%coord(:,iAt2))**2))
        do iU1 = 1, this%nHubbU(species(iAt1))
          u1 = this%uniqHubbU(iU1, iSp1)
          do iU2 = 1, this%nHubbU(species(iAt2))
            u2 = this%uniqHubbU(iU2, iSp2)
            if (iNeigh <= this%nNeighShort(iU2,iU1,iSp2,iAt1)) then
              if (this%tDampedShort(iSp1) .or. this%tDampedShort(iSp2)) then
                this%shortGamma(iU2 ,iU1, iNeigh, iAt1) = expGammaDamped(rab, u2, u1, this%dampExp)
              else
                this%shortGamma(iU2 ,iU1, iNeigh, iAt1) = expGamma(rab, u2, u1)
                if (this%tH5) then
                  call this%h5Correction%scaleShortGamma(&
                      & this%shortGamma(iU2 ,iU1, iNeigh, iAt1), iSp1, iSp2, rab)
                end if
              end if
            end if
          end do
        end do
      end do
    end do

  end subroutine initGamma_


  !> Calculate  the derivative of the short range part of Gamma.
  subroutine addGammaPrime_(this, force, species, iNeighbour, img2CentCell)

    !> Resulting module variables
    type(TScc), intent(in) :: this

    !> force vector to add the short-range part of gamma contribution
    real(dp), intent(inout) :: force(:,:)

    !> List of the species for each atom.
    integer, intent(in) :: species(:)

    !> Index of neighbouring atoms for each atom.
    integer, intent(in) :: iNeighbour(0:,:)

    !> Image of each atom in the central cell.
    integer, intent(in) :: img2CentCell(:)

    integer :: iAt1, iAt2, iAt2f, iU1, iU2, iNeigh, ii, iSp1, iSp2
    real(dp) :: rab, tmpGammaPrime, u1, u2
    real(dp) :: tmpGamma

    @:ASSERT(size(force,dim=1) == 3)
    @:ASSERT(size(force,dim=2) == this%nAtom)
    @:ASSERT(this%tInitialised)

    ! some additional symmetry not used
    do iAt1 = 1, this%nAtom
      iSp1 = species(iAt1)
      do iNeigh = 1, maxval(this%nNeighShort(:,:,:, iAt1))
        iAt2 = iNeighbour(iNeigh, iAt1)
        iAt2f = img2CentCell(iAt2)
        iSp2 = species(iAt2f)
        rab = sqrt(sum((this%coord(:,iAt1) - this%coord(:,iAt2))**2))
        do iU1 = 1, this%nHubbU(species(iAt1))
          u1 = this%uniqHubbU(iU1, iSp1)
          do iU2 = 1, this%nHubbU(species(iAt2f))
            u2 = this%uniqHubbU(iU2, species(iAt2f))
            if (iNeigh <= this%nNeighShort(iU2,iU1,species(iAt2f),iAt1)) then
              if (this%tDampedShort(iSp1) .or. this%tDampedShort(iSp2)) then
                tmpGammaPrime = expGammaDampedPrime(rab, u2, u1, this%dampExp)
              else
                tmpGammaPrime = expGammaPrime(rab, u2, u1)
                if (this%tH5) then
                  tmpGamma = expGamma(rab, u2, u1)
                  call this%h5Correction%scaleShortGammaDeriv(tmpGamma, tmpGammaPrime, iSp1, iSp2,&
                      & rab)
                end if
              end if
              do ii = 1,3
                force(ii,iAt1) = force(ii,iAt1) - this%deltaQUniqU(iU1,iAt1) *&
                    & this%deltaQUniqU(iU2,iAt2f)*tmpGammaPrime*(this%coord(ii,iAt1)&
                    & - this%coord(ii,iAt2))/rab
                force(ii,iAt2f) = force(ii,iAt2f) + this%deltaQUniqU(iU1,iAt1) *&
                    & this%deltaQUniqU(iU2,iAt2f)*tmpGammaPrime*(this%coord(ii,iAt1)&
                    & - this%coord(ii,iAt2))/rab
              end do
            end if
          end do
        end do
      end do
    end do

  end subroutine addGammaPrime_


  !> Calculate  the derivative of the short range part of Gamma.
  subroutine addSTGammaPrime_(st, this, species, iNeighbour, img2CentCell)

    !> Stress tensor component to add the short-range part of the gamma contribution
    real(dp), intent(out) :: st(:,:)

    !> Resulting module variables
    type(TScc), intent(in) :: this

    !> List of the species for each atom.
    integer, intent(in) :: species(:)

    !> Index of neighbouring atoms for each atom.
    integer, intent(in) :: iNeighbour(0:,:)

    !> Image of each atom in the central cell.
    integer, intent(in) :: img2CentCell(:)

    integer :: iAt1, iAt2, iAt2f, iU1, iU2, iNeigh, ii, jj, iSp1, iSp2
    real(dp) :: rab, tmpGammaPrime, u1, u2
    real(dp) :: intermed(3), vect(3)
    ! H5 correction temp. vars
    real(dp) :: tmpGamma
    ! H5 correction end

    @:ASSERT(all(shape(st)==(/3,3/)))
    @:ASSERT(this%tInitialised)

    st(:,:) = 0.0_dp
    ! some additional symmetry not used
    do iAt1 = 1, this%nAtom
      iSp1 = species(iAt1)
      do iNeigh = 1, maxval(this%nNeighShort(:,:,:, iAt1))
        iAt2 = iNeighbour(iNeigh, iAt1)
        iAt2f = img2CentCell(iAt2)
        iSp2 = species(iAt2f)
        vect(:) = this%coord(:,iAt1) - this%coord(:,iAt2)
        rab = sqrt(sum((vect)**2))
        intermed(:) = 0.0_dp
        do iU1 = 1, this%nHubbU(species(iAt1))
          u1 = this%uniqHubbU(iU1, iSp1)
          do iU2 = 1, this%nHubbU(species(iAt2f))
            u2 = this%uniqHubbU(iU2, species(iAt2f))
            if (iNeigh <= this%nNeighShort(iU2,iU1,species(iAt2f),iAt1)) then
              if (this%tDampedShort(iSp1) .or. this%tDampedShort(iSp2)) then
                tmpGammaPrime = expGammaDampedPrime(rab, u2, u1, this%dampExp)
              else
                tmpGammaPrime = expGammaPrime(rab, u2, u1)
                if (this%tH5) then
                  tmpGamma = expGamma(rab, u2, u1)
                  call this%h5Correction%scaleShortGammaDeriv(tmpGamma, tmpGammaPrime, iSp1, iSp2,&
                      & rab)
                end if
              end if
              do ii = 1,3
                intermed(ii) = intermed(ii) &
                    & - this%deltaQUniqU(iU1,iAt1) * this%deltaQUniqU(iU2,iAt2f)&
                    & * tmpGammaPrime*vect(ii)/rab
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

    st(:,:) = st(:,:) / this%volume

  end subroutine addSTGammaPrime_


  !> Calculate gamma integral derivatives in SCC part
  subroutine getGammaDeriv(this, env, species, iNeighbour, img2CentCell, GammaDeriv)

    !> Instance
    class(TScc), intent(in) :: this

    !> Computational environment settings
    type(TEnvironment), intent(in) :: env

    !> list of all atomic species
    integer, intent(in) :: species(:)

    !> neighbour list for atoms
    integer, intent(in) :: iNeighbour(0:,:)

    !> indexing array for periodic image atoms
    integer, intent(in) :: img2CentCell(:)

    !> atom resolved scc gamma derivative, \gamma_{A,B}
    !> gamma_deriv = (-1/R^2 - S')*((x or y,z)/R)
    real(dp), intent(out) :: GammaDeriv(:,:,:)

    real(dp), allocatable :: shortGammaDeriv(:,:,:)
    real(dp), allocatable :: invRDeriv(:,:,:)

    real(dp) :: rab, u1, u2, tmpGamma, tmpGammaPrime
    integer :: iAt1, iSp1, iNeigh, iAt2, iAt2f, iSp2, iU1, iU2, ii

    allocate(shortGammaDeriv(this%nAtom,this%nAtom,3))
    allocate(invRDeriv(this%nAtom,this%nAtom,3))

    ! shortGamma contribution to gamma derivative
    shortGammaDeriv(:,:,:) = 0.0_dp
    do iAt1 = 1, this%nAtom
      iSp1 = species(iAt1)
      do iNeigh = 1, maxval(this%nNeighShort(:,:,:, iAt1))
        iAt2 = iNeighbour(iNeigh, iAt1)
        iAt2f = img2CentCell(iAt2)
        iSp2 = species(iAt2f)
        rab = sqrt(sum((this%coord(:,iAt1) - this%coord(:,iAt2))**2))
        do iU1 = 1, this%nHubbU(species(iAt1))
          u1 = this%uniqHubbU(iU1, iSp1)
          do iU2 = 1, this%nHubbU(species(iAt2f))
            u2 = this%uniqHubbU(iU2, species(iAt2f))
            if (iNeigh <= this%nNeighShort(iU2,iU1,species(iAt2f),iAt1)) then
              if (this%tDampedShort(iSp1) .or. this%tDampedShort(iSp2)) then
                tmpGammaPrime = expGammaDampedPrime(rab, u2, u1, this%dampExp)
              else
                tmpGammaPrime = expGammaPrime(rab, u2, u1)
                if (this%tH5) then
                  tmpGamma = expGamma(rab, u2, u1)
                  call this%h5Correction%scaleShortGammaDeriv(tmpGamma, &
                      & tmpGammaPrime, iSp1, iSp2, rab)
                end if
              end if
              do ii = 1,3
                shortGammaDeriv(iAt2f,iAt1,ii) = -tmpGammaPrime * &
                    & ( this%coord(ii,iAt1) - this%coord(ii,iAt2) ) / rab
              end do
            end if
          end do
        end do
      end do
    end do

    ! 1/R contribution to gamma derivative
    invRDeriv(:,:,:) = 0.0_dp
    if (this%tPeriodic) then
      call this%coulombCont%addInvRPrimePeriodicMat(env, this%coord,&
          & this%coulombCont%neighListGen, this%coulombCont%gLatPoint, this%coulombCont%alpha,&
          & this%volume, invRDeriv)
    else
      call this%coulombCont%addInvRPrimeClusterMat(env, this%coord, invRDeriv)
    end if

    GammaDeriv(:,:,:) = invRDeriv + shortGammaDeriv

  end subroutine getGammaDeriv


  !> Get Q * inverse R contribution for the point charges
  subroutine getShiftOfPC(this, QinvR)

    !> Instance
    class(TScc), intent(in) :: this

    !> (Q * invR) contribution
    real(dp), intent(out) :: QinvR(:)

    call this%extCharge%copyInvRvec(QinvR)

  end subroutine getShiftOfPC


end module dftbp_scc
