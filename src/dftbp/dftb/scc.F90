!--------------------------------------------------------------------------------------------------!
!  DFTB+: general package for performing fast atomistic simulations                                !
!  Copyright (C) 2006 - 2022  DFTB+ developers group                                               !
!                                                                                                  !
!  See the LICENSE file for terms of usage and distribution.                                       !
!--------------------------------------------------------------------------------------------------!

#:include 'common.fypp'

!> Functions and local variables for the SCC calculation.
module dftbp_dftb_scc
  use dftbp_common_accuracy, only : dp
  use dftbp_common_environment, only : TEnvironment
  use dftbp_dftb_boundarycond, only : boundaryConditions, TBoundaryConditions
  use dftbp_dftb_chargeconstr, only : TChrgConstr, TChrgConstr_init
  use dftbp_dftb_charges, only : getSummedCharges
  use dftbp_dftb_coulomb, only : TCoulombInput, TCoulomb, TCoulomb_init
  use dftbp_dftb_extcharges, only : TExtCharges, TExtCharges_init
  use dftbp_dftb_periodic, only : TNeighbourList
  use dftbp_dftb_shortgamma, only : TShortGammaInput, TShortGamma, TShortGamma_init
  use dftbp_dftbplus_elstattypes, only : elstatTypes
  use dftbp_extlibs_poisson, only : TPoissonInput, TPoisson, TPoisson_init
  use dftbp_io_message, only : error
  use dftbp_type_commontypes, only : TOrbitals
  implicit none

  private
  public :: TSccInput, TScc, TScc_init


  !> Data necessary to initialize the SCC module
  type TSccInput

    !> external charges
    real(dp), allocatable :: extCharges(:,:)

    !> if broadened external charges
    real(dp), allocatable :: blurWidths(:)

    !> any constraints on atomic charges
    real(dp), allocatable :: chrgConstraints(:,:)

    !> third order energy contributions
    real(dp), allocatable :: thirdOrderOn(:,:)

    !> Calculator for short gamma
    type(TShortGammaInput), allocatable :: shortGammaInput

    !> Coulomb calculator
    type(TCoulombInput), allocatable :: coulombInput

    !> Poisson solver for calculating electrostatics (instead of shortGamma + coulombCalc)
    type(TPoissonInput), allocatable :: poissonInput

    !> Boundary condition of the system
    integer :: boundaryCond = boundaryConditions%unknown

  end type TSccInput


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

    !> Shift vector per atom
    real(dp), allocatable :: shiftPerAtom(:)

    !> Shift vector per l-shell
    real(dp), allocatable :: shiftPerL(:,:)

    !> Atomic coordinates
    real(dp), allocatable :: coord(:,:)

    !> Cell volume
    real(dp) :: volume

    !> Negative gross charge
    real(dp), allocatable :: deltaQ(:,:)

    !> Negative gross charge per shell
    real(dp), allocatable :: deltaQShell(:,:)

    !> Negative gross charge per atom
    real(dp), allocatable :: deltaQAtom(:)

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
    type(TExtCharges), allocatable :: extCharges

    !> Coulombic interaction container
    ! TODO(BA): Prevent direct access from REKS and make it private
    type(TCoulomb), allocatable, public :: coulomb

    !> Calculator for short range gamma interaction
    type(TShortGamma), allocatable :: shortGamma

    !> Poisson solver for calculating electrostatics (instead of shortGamma + coulomb)
    type(TPoisson), allocatable :: poisson

    !> Which electrostatic solver should be used?
    integer :: elstatType

  contains

    !> Returns a minimal cutoff for the neighbourlist
    procedure :: getCutOff

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

    !> Triggers all instructions which must be done once the SCC-loop had been finished
    procedure :: finishSccLoop

  end type TScc


  ! Boundary conditions the module can handle
  integer, parameter :: implementedBoundaryConds_(*) =&
      & [boundaryConditions%cluster, boundaryConditions%pbc3d]


contains


  !> Initialize SCC container from input data
  subroutine TScc_init(this, input, env, orb)

    !> Resulting instance
    type(TScc), intent(out) :: this

    !> Scc input
    type(TSccInput), intent(inout) :: input

    !> Environment
    type(TEnvironment), intent(inout) :: env

    !> Orbital information
    type(TOrbitals), intent(in) :: orb

    @:ASSERT(.not. this%tInitialised)
    @:ASSERT(allocated(input%extCharges) .or. .not. allocated(input%blurWidths))
    @:ASSERT(any(implementedBoundaryConds_ == input%boundaryCond))
    @:ASSERT(allocated(input%shortGammaInput) .eqv. allocated(input%coulombInput))
    @:ASSERT(allocated(input%shortGammaInput) .neqv. allocated(input%poissonInput))

    this%nSpecies = size(orb%nOrbSpecies)
    this%nAtom = size(orb%nOrbAtom)
    this%mShell = orb%mShell
    this%mOrb = orb%mOrb

    allocate(this%shiftPerAtom(this%nAtom))
    allocate(this%shiftPerL(this%mShell, this%nAtom))

    if (allocated(input%shortGammaInput)) then
      this%elstatType = elstatTypes%gammaFunc
    else
      this%elstatType = elstatTypes%poisson
    end if

    select case (this%elstatType)
    case (elstatTypes%gammaFunc)
      call initShortGamma_(input%shortGammaInput, orb, this%shortGamma)
      call initCoulomb_(input%coulombInput, env, this%nAtom, this%coulomb)
    case (elstatTypes%poisson)
      #:block REQUIRES_COMPONENT('Poisson-solver', WITH_POISSON)
        call initPoisson_(input%poissonInput, env, orb, this%poisson)
      #:endblock
    end select

    this%tPeriodic = (input%boundaryCond == boundaryConditions%pbc3d)

    ! Initialise external charges
    if (allocated(input%extCharges)) then
      if (.not. allocated(this%coulomb)) then
        call error("External charges require Coulomb calculator (Gamma-electrostatics)")
      end if
      call this%setExternalCharges(input%extCharges(1:3,:), input%extCharges(4,:),&
          & blurWidths=input%blurWidths)
    end if

    this%tChrgConstr = allocated(input%chrgConstraints)
    if (this%tChrgConstr) then
      allocate(this%chrgConstr)
      call TChrgConstr_init(this%chrgConstr, input%chrgConstraints, 2)
    end if
    this%tThirdOrder = allocated(input%thirdOrderOn)
    if (this%tThirdOrder) then
      allocate(this%thirdOrder)
      ! Factor 1/6 in the energy is put into the Hubbard derivatives
      call TChrgConstr_init(this%thirdOrder, input%thirdOrderOn / 6.0_dp, 3)
    end if

    ! Initialise arrays for charge differences
    allocate(this%deltaQ(this%mOrb, this%nAtom))
    allocate(this%deltaQShell(this%mShell, this%nAtom))
    allocate(this%deltaQAtom(this%nAtom))

    this%tInitialised = .true.

  end subroutine TScc_init


  !> Returns a minimal cutoff for the neighbourlist, which must be passed to various functions in
  !> this module.
  function getCutOff(this) result(cutoff)

    !> Instance
    class(TScc), intent(in) :: this

    !> The neighbourlists, passed to scc routines, should contain neighbour information at
    !> least up to that cutoff.
    real(dp) :: cutoff

    @:ASSERT(this%tInitialised)

    select case (this%elstatType)
    case (elstatTypes%gammaFunc)
      cutoff = this%shortGamma%getCutOff()
    case default
      cutoff = 0.0_dp
    end select

  end function getCutOff


  !> Updates the atom coordinates for the SCC module.
  subroutine updateCoords(this, env, coord0, coord, species, neighList)

    !> Instance
    class(TScc), intent(inout) :: this

    !> Environment settings
    type(TEnvironment), intent(in) :: env

    !> New unfolded coordinates of the atoms in the central cell
    real(dp), intent(in) :: coord0(:,:)

    !> New coordinates of the atoms
    real(dp), intent(in) :: coord(:,:)

    !> Species of the atoms (should not change during run)
    integer, intent(in) :: species(:)

    !> Neighbour list for the atoms.
    type(TNeighbourList), intent(in) :: neighList

    @:ASSERT(this%tInitialised)

    this%coord = coord

    select case (this%elstatType)
    case (elstatTypes%gammaFunc)
      call this%coulomb%updateCoords(env, neighList, coord, species)
      call this%shortGamma%updateCoords(coord, species, neighList)
    case (elstatTypes%poisson)
      #:block REQUIRES_COMPONENT('Poisson-solver', WITH_POISSON)
        ! Poisson solver needs currently the unfolded central cell coordinates (coord0), because the
        ! folding can mess up the contact position.
        call this%poisson%updateCoords(coord0)
      #:endblock
    end select

    if (allocated(this%extCharges)) then
      call this%extCharges%setCoordinates(env, coord(:, 1:this%nAtom), this%coulomb)
    end if

  end subroutine updateCoords


  !> Updates the SCC module, if the lattice vectors had been changed
  subroutine updateLatVecs(this, latVec, recVec, boundaryConds, vol)

    !> Instance
    class(TScc), intent(inout) :: this

    !> New lattice vectors
    real(dp), intent(in) :: latVec(:,:)

    !> New reciprocal lattice vectors
    real(dp), intent(in) :: recVec(:,:)

    !> Boundary conditions on the calculation
    type(TBoundaryConditions), intent(in) :: boundaryConds

    !> New volume
    real(dp), intent(in) :: vol

    @:ASSERT(this%tInitialised)
    @:ASSERT(this%tPeriodic)

    this%volume = vol

    select case (this%elstatType)
    case (elstatTypes%gammaFunc)
      call this%coulomb%updateLatVecs(latVec, recVec, vol)
    case (elstatTypes%poisson)
      #:block REQUIRES_COMPONENT('Poisson-solver', WITH_POISSON)
        call this%poisson%updateLatVecs(latVec)
      #:endblock
    end select

    if (allocated(this%extCharges)) then
      call this%extCharges%setLatticeVectors(latVec, boundaryConds)
    end if

  end subroutine updateLatVecs


  !> Updates the SCC module, if the charges have been changed
  subroutine updateCharges(this, env, qOrbital, orb, species, q0)

    !> Resulting module variables
    class(TScc), intent(inout) :: this

    !> Environment settings
    type(TEnvironment), intent(inout) :: env

    !> Orbital resolved charges. Shape: [mOrb, nAtom, nSpin].
    real(dp), intent(in) :: qOrbital(:,:,:)

    !> Contains information about the atomic orbitals in the system
    type(TOrbitals), intent(in) :: orb

    !> Species of the atoms (should not change during run). Shape: [nSpecies]
    integer, intent(in) :: species(:)

    !> Reference charge distribution (neutral atoms). Shape: [mOrb, nAtom, nSpin]
    real(dp), intent(in), optional :: q0(:,:,:)

    @:ASSERT(this%tInitialised)

    call getSummedCharges(species, orb, qOrbital, q0, dQ=this%deltaQ, dQAtom=this%deltaQAtom,&
        & dQShell=this%deltaQShell)

    select case (this%elstatType)
    case (elstatTypes%gammaFunc)
      call this%shortGamma%updateCharges(orb, species, this%deltaQShell)
      call this%coulomb%updateCharges(env, qOrbital, orb, species, this%deltaQ, this%deltaQAtom,&
          & this%deltaQShell)
    case (elstatTypes%poisson)
      #:block REQUIRES_COMPONENT('Poisson-solver', WITH_POISSON)
        call this%poisson%updateCharges(env, qOrbital(:,:,1), q0)
      #:endblock
    end select

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

    select case (this%elstatType)
    case (elstatTypes%gammaFunc)
      call this%coulomb%updateShifts(env, orb, species, iNeighbour, img2CentCell)
      this%shiftPerAtom(:) = 0.0_dp
      call this%shortGamma%updateShifts(env, orb, species, iNeighbour, img2CentCell)
      call this%shortGamma%getShiftPerShell(this%shiftPerL)
    case (elstatTypes%poisson)
      #:block REQUIRES_COMPONENT('Poisson-solver', WITH_POISSON)
        this%shiftPerAtom(:) = 0.0_dp
        this%shiftPerL(:,:) = 0.0_dp
        call this%poisson%addPotentials(this%shiftPerL)
      #:endblock
    end select

    if (this%tChrgConstr) then
      call this%chrgConstr%buildShift(this%deltaQAtom)
    end if
    if (this%tThirdOrder) then
      call this%thirdOrder%buildShift(this%deltaQAtom)
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

    if (.not. allocated(this%extCharges)) then
      allocate(this%extCharges)
      call TExtCharges_init(this%extCharges, this%nAtom, size(chargeQs), this%tPeriodic)
    end if

    call this%extCharges%setExternalCharges(chargeCoords, chargeQs, blurWidths=blurWidths)

    if (present(blurWidths)) then
      if (this%tPeriodic .and. any(blurWidths > 1.0e-7_dp)) then
        if (1.0_dp / maxval(blurWidths) < this%coulomb%alpha) then
          call error("Charge blur widths are too wide compared to the Ewald real space sum")
        end if
      end if
    else
      call this%extCharges%setExternalCharges(chargeCoords, chargeQs)
    end if

  end subroutine setExternalCharges


  !> Routine for returning lower triangle of atomic resolved gamma as a matrix
  !>
  !> Works only, if SCC-instance uses Gamma-electrostatics.
  !>
  subroutine getAtomicGammaMatrix(this, gammamat, iNeighbour, img2CentCell)

    !> Instance
    class(TScc), intent(in) :: this

    !> Atom resolved gamma
    real(dp), intent(out) :: gammamat(:,:)

    !> neighbours of atoms
    integer, intent(in) :: iNeighbour(0:,:)

    !> index array between images and central cell
    integer, intent(in) :: img2CentCell(:)

    @:ASSERT(this%tInitialised)
    @:ASSERT(all(shape(gammamat) == [ this%nAtom, this%nAtom ]))
    @:ASSERT(this%elstatType == elstatTypes%gammaFunc)

  #:if WITH_SCALAPACK
    call error("scc:getAtomicGammaMatrix does not work with MPI yet")
  #:endif

    gammamat(:,:) = this%coulomb%invRMat
    call this%shortGamma%addAtomicMatrix(gammamat, iNeighbour, img2CentCell)

  end subroutine getAtomicGammaMatrix


  !> Routine for returning lower triangle of atomic resolved Coulomb matrix
  !>
  !> Works only, if SCC-instance uses Gamma-electrostatics.
  !>
  subroutine getAtomicGammaMatU(this, gammamat, hubbU, species, iNeighbour, img2CentCell)

    !> Instance
    class(TScc), intent(in) :: this

    !> Atom resolved gamma
    real(dp), intent(out) :: gammamat(:,:)

    !> ppRPA Hubbard parameters
    real(dp), intent(in) :: hubbU(:)

    !> List of the species for each atom.
    integer,  intent(in) :: species(:)

    !> neighbours of atoms
    integer, intent(in) :: iNeighbour(0:,:)

    !> index array between images and central cell
    integer, intent(in) :: img2CentCell(:)

    integer  :: iAt1, iAt2

    @:ASSERT(this%tInitialised)
    @:ASSERT(all(shape(gammamat) == [this%nAtom, this%nAtom]))
    @:ASSERT(this%elstatType == elstatTypes%gammaFunc)

  #:if WITH_SCALAPACK
    call error("scc:getAtomicGammaMatU does not work with MPI yet")
  #:endif

    gammamat(:,:) = this%coulomb%invRMat
    call this%shortGamma%addAtomicMatrixCustomU(gammamat, hubbU, species, this%coord, iNeighbour,&
        & img2CentCell)
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
        & + sum(this%shiftPerL * this%deltaQShell, dim=1))
    if (this%elstatType == elstatTypes%gammaFunc) then
      call this%coulomb%addEnergy(eScc)
    end if

    if (allocated(this%extCharges)) then
      call this%extCharges%addEnergyPerAtom(this%deltaQAtom, eScc)
    end if

    if (this%tChrgConstr) then
      call this%chrgConstr%addEnergyPerAtom(eScc, this%deltaQAtom)
    end if

    if (this%tThirdOrder) then
      call this%thirdOrder%addEnergyPerAtom(eScc, this%deltaQAtom)
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
        & + sum(this%shiftPerL * (2.0_dp * dQOutShell - this%deltaQShell), dim=1))
    if (this%elstatType == elstatTypes%gammaFunc) then
      call this%coulomb%addEnergy(eScc, dQOut, dQOutAtom, dQOutShell)
    end if

    if (allocated(this%extCharges)) then
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
    class(TScc), intent(inout) :: this

    !> Environment settings
    type(TEnvironment), intent(inout) :: env

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

    real(dp), allocatable :: tmpDerivs(:,:)

    @:ASSERT(this%tInitialised)
    @:ASSERT(size(force,dim=1) == 3)
    @:ASSERT(size(force,dim=2) == this%nAtom)
    @:ASSERT(present(chrgForce) .eqv. allocated(this%extCharges))

    select case (this%elstatType)
    case (elstatTypes%gammaFunc)
      call this%shortGamma%addGradientsDc(force, species, this%coord, iNeighbour, img2CentCell)
      call this%coulomb%addGradients(env, this%coord, species, iNeighbour, img2CentCell, force)
    case (elstatTypes%poisson)
      #:block REQUIRES_COMPONENT('Poisson-solver', WITH_POISSON)
        allocate(tmpDerivs, mold=force)
        call this%poisson%getGradients(env, tmpDerivs)
        force(:,:) = force + tmpDerivs
      #:endblock
    end select

    if (allocated(this%extCharges)) then
      call this%extCharges%addForceDc(env, force, chrgForce, this%coord(:, 1:this%nAtom),&
          & this%deltaQAtom, this%coulomb)
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

    @:ASSERT(this%tInitialised)
    @:ASSERT(this%tPeriodic)
    @:ASSERT(all(shape(st)==(/3,3/)))

    select case (this%elstatType)
    case (elstatTypes%gammaFunc)
      call this%shortGamma%addStressDc(st, this%coord, species, this%volume, iNeighbour,&
          & img2CentCell)
      call this%coulomb%addStress(env, this%coord, species, iNeighbour, img2CentCell, st)
    case (elstatTypes%poisson)
      #:block REQUIRES_COMPONENT('Poisson-solver', WITH_POISSON)
        st(:,:) = 0.0_dp
      #:endblock
    end select

    !! NOTE: no stress contribution from external charges

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
    if (this%elstatType == elstatTypes%gammaFunc) then
      call this%coulomb%addShiftPerAtom(shift)
    end if
    if (allocated(this%extCharges)) then
      call this%extCharges%addShiftPerAtom(shift)
    end if
    if (this%tChrgConstr) then
      call this%chrgConstr%addShiftPerAtom(shift)
    end if
    if (this%tThirdOrder) then
      call this%thirdOrder%addShiftPerAtom(shift)
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
    if (this%elstatType == elstatTypes%gammaFunc) then
      call this%coulomb%addShiftPerShell(shift)
    end if

  end subroutine getShiftPerL


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
    real(dp), allocatable :: dQOutLShell(:,:)

    if (this%elstatType /= elstatTypes%gammaFunc) then
      call error("XLBOMD only works with gamma-electrostatics")
    end if

    allocate(dQOut(this%mOrb, this%nAtom))
    allocate(dQOutAtom(this%nAtom))
    allocate(dQOutLShell(this%mShell, this%nAtom))

    call getSummedCharges(species, orb, qOrbitalOut, q0, dQ=dQOut, dQAtom=dQOutAtom,&
        & dQShell=dQOutLShell)

    call this%shortGamma%addGradientsDcXlbomd(this%deltaQShell, dQOutLShell, this%coord, species,&
        & orb, iNeighbour, img2CentCell, force)
    call this%coulomb%addGradients(env, this%coord, species, iNeighbour, img2CentCell, force,&
        & dQOut, dQOutAtom, dQOutLShell)

    if (allocated(this%extCharges)) then
      call error("XLBOMD with external charges does not work yet!")
    end if

  end subroutine addForceDcXlbomd


  !> Returns potential from DFTB charges
  subroutine getInternalElStatPotential(this, pot, env, locations, epsSoften)

    !> Instance of SCC calculation
    class(TScc), intent(in) :: this

    !> Resulting potentials
    real(dp), intent(out) :: pot(:)

    !> Computational environment settings
    type(TEnvironment), intent(in) :: env

    !> sites to calculate potential
    real(dp), intent(in) :: locations(:,:)

    !> optional potential softening
    real(dp), optional, intent(in) :: epsSoften

    @:ASSERT(this%tInitialised)
    @:ASSERT(all(shape(locations) == [3,size(pot)]))

    if (this%elstatType /= elstatTypes%gammaFunc) then
      call error("getInternalElStatPotential only works with gamma-electrostatics")
    end if

    pot(:) = 0.0_dp
    call this%coulomb%getPotential(env, locations, this%coord, this%deltaQAtom, pot,&
        & epsSoften=epsSoften)

  end subroutine getInternalElStatPotential


  !> Returns potential from external charges
  subroutine getExternalElStatPotential(this, V, env, locations, epsSoften)

    !> Instance
    class(TScc), intent(inout) :: this

    !> Resulting potentials
    real(dp), intent(out) :: V(:)

    !> Computational environment settings
    type(TEnvironment), intent(in) :: env

    !> sites to calculate potential
    real(dp), intent(in) :: locations(:,:)

    !> optional potential softening
    real(dp), optional, intent(in) :: epsSoften

    @:ASSERT(this%tInitialised)

    if (allocated(this%extCharges)) then
      call this%extCharges%getElStatPotential(env, locations, V, this%coulomb,&
          & epsSoften=epsSoften)
    else
      V(:) = 0.0_dp
    end if

  end subroutine getExternalElStatPotential


  !> Calculate gamma integral derivatives in SCC part
  subroutine getGammaDeriv(this, env, species, iNeighbour, img2CentCell, gammaDeriv)

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
    real(dp), intent(out) :: gammaDeriv(:,:,:)

    if (this%elstatType /= elstatTypes%gammaFunc) then
      call error("getGammaDeriv only works with gamma-electrostatics")
    end if

    gammaDeriv(:,:,:) = 0.0_dp
    call this%shortGamma%addDerivativeMatrix(this%coord, species, iNeighbour, img2CentCell,&
        & gammaDeriv)
    if (this%tPeriodic) then
      call this%coulomb%addInvRPrimePeriodicMat(env, this%coord,gammaDeriv)
    else
      call this%coulomb%addInvRPrimeClusterMat(env, this%coord, gammaDeriv)
    end if

  end subroutine getGammaDeriv


  !> Get Q * inverse R contribution for the point charges
  subroutine getShiftOfPC(this, QinvR)

    !> Instance
    class(TScc), intent(in) :: this

    !> (Q * invR) contribution
    real(dp), intent(out) :: QinvR(:)

    call this%extCharges%copyInvRvec(QinvR)

  end subroutine getShiftOfPC


  !> Triggers all instructions which must be done once the SCC-loop had been finished
  subroutine finishSccLoop(this, env)

    !> Instance
    class(TScc), intent(inout) :: this

    !> Computational environment settings
    type(TEnvironment), intent(inout) :: env

    if (this%elstatType == elstatTypes%poisson) then
      #:block REQUIRES_COMPONENT('Poisson-solver', WITH_POISSON)
        call this%poisson%savePotential(env)
      #:endblock
    end if

  end subroutine finishSccLoop


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!  Private routines
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

#:if WITH_POISSON

  ! Initializes the Poisson solver
  subroutine initPoisson_(poissonInput, env, orb, poisson)
    type(TPoissonInput), intent(inout) :: poissonInput
    type(TEnvironment), intent(inout) :: env
    type(TOrbitals), intent(in) :: orb
    type(TPoisson), allocatable, intent(out) :: poisson

    logical :: success

    allocate(poisson)
    call TPoisson_init(poisson, poissonInput, env, orb, success)
    if (.not. success) then
      call error("Poisson solver not initialized")
    end if

  end subroutine initPoisson_

#:endif


  ! Initializes the short gamma calculator
  subroutine initShortGamma_(shortGammaInput, orb, shortGamma)
    type(TShortGammaInput), intent(inout) :: shortGammaInput
    type(TOrbitals), intent(in) :: orb
    type(TShortGamma), allocatable, intent(out) :: shortGamma

    allocate(shortGamma)
    call TShortGamma_init(shortGamma, shortGammaInput, orb)

  end subroutine initShortGamma_


  ! Initializes the Coulomb calculator
  subroutine initCoulomb_(coulombInput, env, nAtom, coulomb)
    type(TCoulombInput), intent(inout) :: coulombInput
    type(TEnvironment), intent(inout) :: env
    integer, intent(in) :: nAtom
    type(TCoulomb), allocatable, intent(out) :: coulomb

    allocate(coulomb)
    call TCoulomb_init(coulomb, coulombInput, env, nAtom)

  end subroutine initCoulomb_


end module dftbp_dftb_scc
