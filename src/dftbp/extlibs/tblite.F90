!--------------------------------------------------------------------------------------------------!
!  DFTB+: general package for performing fast atomistic simulations                                !
!  Copyright (C) 2006 - 2022  DFTB+ developers group                                               !
!                                                                                                  !
!  See the LICENSE file for terms of usage and distribution.                                       !
!--------------------------------------------------------------------------------------------------!
#:include "common.fypp"
#:include "error.fypp"

!> Proxy module for interfacing with the tblite library.
!>
!> @note
!> The library calculates energies, gradients (∂E/∂R) and strain derivatives (∂E/∂ε = –V·σ),
!> while DFTB+ works with atom-resolved energies, gradients and the stress tensor.
!> The atom-resolved energy partitioning is archived by distributing the energy equivalently
!> over all atoms. The strain derivative is saved as such in the container and only
!> transformed to the stress tensor on demand using σ = –1/V·∂E/∂ε.
!>
!> @warning
!> This module has to account for changing between sign conventions between DFTB+ and tblite.
!> Generally, all intermediate quantities passed from DFTB+ to the library are using the
!> conventions of DFTB+, while all intermediate quantities passed from the library to DFTB+
!> will follow tblite's conventions (usually encapsulated in derived types already).
!>
!> Both tblite and DFTB+ use consistent ordering of spherical harmonics
!> in the standard sorting, *i.e.* [-l, ..., 0, ..., l].
module dftbp_extlibs_tblite
  use dftbp_common_accuracy, only : dp
  use dftbp_common_environment, only : TEnvironment
  use dftbp_common_schedule, only : distributeRangeInChunks, assembleChunks
  use dftbp_dftb_charges, only : getSummedCharges
  use dftbp_dftb_periodic, only : TNeighbourList
  use dftbp_io_message, only : error
  use dftbp_math_blasroutines, only : gemv
  use dftbp_math_simplealgebra, only : determinant33
  use dftbp_type_commontypes, only : TOrbitals
  use dftbp_type_integral, only : TIntegral
#:if WITH_TBLITE
  use mctc_env, only : error_type
  use mctc_io, only : structure_type, new
  use mctc_io_symbols, only : symbol_length
  use tblite_basis_type, only : get_cutoff, basis_type
  use tblite_context_type, only : context_type
  use tblite_container, only : container_cache
  use tblite_cutoff, only : get_lattice_points
  use tblite_integral_multipole, only : multipole_cgto, multipole_grad_cgto, maxl, msao
  use tblite_param, only : param_record
  use tblite_scf_info, only : scf_info, atom_resolved, shell_resolved, orbital_resolved, &
      & not_used
  use tblite_scf_potential, only : potential_type, new_potential
  use tblite_version, only : get_tblite_version
  use tblite_wavefunction_type, only : wavefunction_type, new_wavefunction
  use tblite_xtb_calculator, only : xtb_calculator, new_xtb_calculator
  use tblite_xtb_gfn1, only : new_gfn1_calculator
  use tblite_xtb_gfn2, only : new_gfn2_calculator
  use tblite_xtb_h0, only : get_selfenergy, get_hamiltonian, get_occupation, &
      & get_hamiltonian_gradient, tb_hamiltonian
  use tblite_xtb_ipea1, only : new_ipea1_calculator
  use tblite_xtb_singlepoint, only : xtb_singlepoint
#:endif
  implicit none
  private

  public :: TTBLite, TTBLiteInput, TTBLite_init, writeTBLiteInfo
  public :: tbliteMethod


  !> Possible methods available in this library
  type :: EnumMethod

    !> Selector for GFN2-xTB Hamiltonian
    integer :: gfn2xtb = 2

    !> Selector for GFN1-xTB Hamiltonian
    integer :: gfn1xtb = 1

    !> Selector for IPEA1-xTB Hamiltonian
    integer :: ipea1xtb = 11

  end type EnumMethod

  !> Actual numerated method selector
  type(EnumMethod), parameter :: tbliteMethod = EnumMethod()


  !> Information on the library setup
  type :: TBLiteInfo

    !> Parametrization name
    character(len=:), allocatable :: name

  end type TBLiteInfo


  !> Input for the library
  type :: TTBLiteInput

  #:if WITH_TBLITE
    !> Molecular structure data
    type(structure_type) :: mol

    !> Parametrisation data
    type(xtb_calculator) :: calc
  #:endif

    !> Parametrization info
    type(TBLiteInfo) :: info

  contains

    !> Create geometry data for library
    procedure :: setupGeometry

    !> Create parametrization data for library
    generic :: setupCalculator => setupCalculatorFromEnum, setupCalculatorFromFile
    procedure :: setupCalculatorFromEnum
    procedure :: setupCalculatorFromFile

    !> Create orbital information from input data
    procedure :: setupOrbitals

  end type TTBLiteInput


  !> Library interface handler
  type :: TTBLite
  #:if WITH_TBLITE
    !> Molecular structure data
    type(structure_type) :: mol

    !> Calculation context
    type(context_type) :: ctx

    !> Parametrisation data
    type(xtb_calculator) :: calc

    !> Wavefunction data
    type(wavefunction_type) :: wfn

    !> Density-dependent potential data
    type(potential_type) :: pot

    !> Reuseable data for Coulombic interactions
    type(container_cache) :: cache

    !> Reuseable data for Dispersion interactions
    type(container_cache) :: dcache
  #:endif

    !> Parametrization info
    type(TBLiteInfo) :: info

    !> Mapping between species and identifiers
    integer, allocatable :: sp2id(:)

    !> Coordination number
    real(dp), allocatable :: cn(:)

    !> Derivative of the coordination number w.r.t. atomic displacements
    real(dp), allocatable :: dcndr(:, :, :)

    !> Derivative of the coordination number w.r.t. strain deformations
    real(dp), allocatable :: dcndL(:, :, :)

    !> Diagonal elements of the Hamiltonian
    real(dp), allocatable :: selfenergy(:)

    !> Derivatives of the diagonal elements w.r.t. the coordination number
    real(dp), allocatable :: dsedcn(:)

    !> Repulsion energy
    real(dp), allocatable :: erep(:)

    !> Halogen bonding energy
    real(dp), allocatable :: ehal(:)

    !> Non-self consistent dispersion energy
    real(dp), allocatable :: edisp(:)

    !> Self-consistent dispersion energy
    real(dp), allocatable :: escd(:)

    !> Electrostatic energy
    real(dp), allocatable :: ees(:)

    !> Contributions to the gradient
    real(dp), allocatable :: gradient(:, :)

    !> Contributions to the virial
    real(dp) :: sigma(3, 3)

  contains

    !> Update internal copy of coordinates
    procedure :: updateCoords

    !> Update internal copy of lattice vectors
    procedure :: updateLatVecs

    !> Get real space cutoff
    procedure :: getRCutoff

    !> Get energy contributions
    procedure :: getEnergies

    !> Get force contributions
    procedure :: addGradients

    !> Get stress tensor contributions
    procedure :: getStress

    !> Updates with changed charges for the instance.
    procedure :: updateCharges

    !> Returns shifts per atom
    procedure :: getShifts

    !> Get orbital information
    procedure :: getOrbitalInfo

    !> Get information about required multipolar contributions
    procedure :: getMultipoleInfo

    !> Get reference occupation
    procedure :: getReferenceN0

    !> Returns the equivalence to get the correct mixing of charge dependent contributions
    procedure :: getOrbitalEquiv

    !> Get Hubbard parameters from second order electrostatic
    procedure :: getHubbardU

    !> Remove second order electrostatics
    procedure :: removeES2

    !> Construct Hamiltonian and overlap related integrals
    procedure :: buildSH0

    !> Evaluate shift related derivatives from Hamiltonian and overlap related integrals
    procedure :: buildDerivativeShift

    !> Calculates nonadiabatic matrix: overlap gradient (Sprime) times velocities (Rdot)
    procedure :: buildRdotSprime

  end type TTBLite


  !> Number of dipole components used in tblite library (x, y, z)
  integer, parameter :: dimDipole = 3

  !> Number of quadrupole components used in tblite library (xx, xy, yy, xz, yz, zz)
  integer, parameter :: dimQuadrupole = 6


contains


  !> Setup geometry information for input data
  subroutine setupGeometry(this, nAtom, species0, coords0, speciesNames, latVecs)

    !> Input data
    class(TTBLiteInput), intent(inout) :: this

    !> Nr. of atoms in the system
    integer, intent(in) :: nAtom

    !> Species of every atom in the unit cell
    integer, intent(in) :: species0(:)

    !> Atomic coordinates in the unit cell
    real(dp), intent(in) :: coords0(:,:)

    !> Symbols of the species
    character(len=*), intent(in) :: speciesNames(:)

    !> Lattice vectors, if the system is periodic
    real(dp), intent(in), optional :: latVecs(:,:)

  #:if WITH_TBLITE
    character(len=symbol_length), allocatable :: symbol(:)

    symbol = speciesNames(species0)
    call new(this%mol, symbol, coords0, lattice=latVecs)

    if (any(this%mol%num <= 0)) then
      call error("Unidentified species present in species list")
    end if
  #:else
    call notImplementedError
  #:endif
  end subroutine setupGeometry


  !> Setup calculator for input data
  subroutine setupCalculatorFromEnum(this, method)

    !> Input data
    class(TTBLiteInput), intent(inout) :: this

    !> Selected method
    integer, intent(in) :: method

  #:if WITH_TBLITE
    call getCalculator(method, this%mol, this%calc)
  #:else
    call notImplementedError
  #:endif
  end subroutine setupCalculatorFromEnum


  !> Setup calculator for input data
  subroutine setupCalculatorFromFile(this, method)

    !> Input data
    class(TTBLiteInput), intent(inout) :: this

    !> Selected method
    character(len=*), intent(in) :: method

  #:if WITH_TBLITE
    type(param_record) :: param
    type(error_type), allocatable :: err

    call param%load(method, err)
    if (allocated(err)) then
      call error(err%message)
    end if
    call new_xtb_calculator(this%calc, this%mol, param, err)
    if (allocated(err)) then
      call error(err%message)
    end if
    this%info%name = param%name
  #:else
    call notImplementedError
  #:endif
  end subroutine setupCalculatorFromFile


  !> Setup orbital information from input data
  subroutine setupOrbitals(this, species0, orb)

    !> Input data
    class(TTBLiteInput), intent(in) :: this

    !> Species of every atom in the unit cell
    integer, intent(in) :: species0(:)

    !> Orbital information
    type(TOrbitals), intent(out) :: orb

  #:if WITH_TBLITE
    call setupOrbitalInfo(this%calc%bas, this%mol%id, species0, orb)
  #:else
    call notImplementedError
  #:endif
  end subroutine setupOrbitals


  !> Constructor for the library interface
  subroutine TTBLite_init(this, input, nAtom, species0, speciesNames, coords0, latVecs)

    !> Instance of the library interface
    type(TTBLite), intent(out) :: this

    !> Input to construct the library interface from
    type(TTBLiteInput), intent(in) :: input

    !> Nr. of atoms in the system
    integer, intent(in) :: nAtom

    ! Spin channels in the system
    integer, parameter :: nSpin = 1

    !> Species of every atom in the unit cell
    integer, intent(in) :: species0(:)

    !> Atomic coordinates in the unit cell
    real(dp), intent(in) :: coords0(:,:)

    !> Symbols of the species
    character(len=*), intent(in) :: speciesNames(:)

    !> Lattice vectors, if the system is periodic
    real(dp), intent(in), optional :: latVecs(:,:)

  #:if WITH_TBLITE
    type(scf_info) :: info

    this%mol = input%mol
    this%calc = input%calc
    this%info = input%info

    info = this%calc%variable_info()
    if (info%charge > shell_resolved) then
      call error("Library interface does not support orbital-resolved charge communication")
    end if
    if (info%dipole > atom_resolved) then
      call error("Library interface does not support shell-resolved dipole moment communication")
    end if
    if (info%quadrupole > atom_resolved) then
      call error("Library interface does not support shell-resolved quadrupole moment communication")
    end if

    call new_wavefunction(this%wfn, this%mol%nat, this%calc%bas%nsh, this%calc%bas%nao, &
        & nSpin, 0.0_dp)

    call new_potential(this%pot, this%mol, this%calc%bas, this%wfn%nspin)

    if (allocated(this%calc%ncoord)) then
      allocate(this%cn(this%mol%nat))
      allocate(this%dcndr(3, this%mol%nat, this%mol%nat), this%dcndL(3, 3, this%mol%nat))
    end if

    allocate(this%selfenergy(this%calc%bas%nsh), this%dsedcn(this%calc%bas%nsh))

    allocate(this%erep(this%mol%nat), this%ehal(this%mol%nat), this%edisp(this%mol%nat), &
        & this%escd(this%mol%nat), this%ees(this%mol%nat))
    allocate(this%gradient(3, this%mol%nat))

    allocate(this%sp2id(maxval(species0)))
    call getSpeciesIdentifierMap(this%sp2id, species0, this%mol%id)
  #:else
    call notImplementedError
  #:endif
  end subroutine TTBLite_init


#:if WITH_TBLITE
  subroutine getCalculator(method, mol, calc)

    !> Selected method
    integer, intent(in) :: method

    !> Molecular structure data
    type(structure_type), intent(in) :: mol

    !> Parametrisation data
    type(xtb_calculator), intent(out) :: calc

    select case(method)
    case default
      call error("Unknown method selector")
    case(tbliteMethod%gfn2xtb)
      call new_gfn2_calculator(calc, mol)
    case(tbliteMethod%gfn1xtb)
      call new_gfn1_calculator(calc, mol)
    case(tbliteMethod%ipea1xtb)
      call new_ipea1_calculator(calc, mol)
    end select
  end subroutine getCalculator
#:endif


  subroutine getSpeciesIdentifierMap(sp2id, species, id)

    !> Mapping from species to identifiers
    integer, intent(out) :: sp2id(:)

    !> Element species used in DFTB+
    integer, intent(in) :: species(:)

    !> Element identifiers used in tblite
    integer, intent(in) :: id(:)

    integer :: nSpecies, iAt, iSp, iId
    logical, allocatable :: done(:)

    nSpecies = maxval(species)
    allocate(done(nSpecies))
    done(:) = .false.
    do iAt = 1, size(species)
      iId = id(iAt)
      iSp = species(iAt)
      if (done(iSp)) cycle
      sp2id(iSp) = iId
      done(iSp) = .true.
    end do
  end subroutine getSpeciesIdentifierMap


  !> Write information about library setup
  subroutine writeTBLiteInfo(unit, this)

    !> Formatted unit for output
    integer, intent(in) :: unit

    !> Data structure
    type(TTBLite), intent(in) :: this

  #:if WITH_TBLITE
    character(len=:), allocatable :: version_string
    character(len=*), parameter :: fmt = '(a, ":", t30, a)'

    call get_tblite_version(string=version_string)
    write(unit, fmt) "tblite library version", version_string
    if (allocated(this%info%name)) then
      write(unit, fmt) "-> parametrization", this%info%name
    end if
    write(unit, fmt) "-> repulsion", abool(allocated(this%calc%repulsion))
    write(unit, fmt) "-> dispersion", abool(allocated(this%calc%dispersion))
    write(unit, fmt) "-> halogen bonding", abool(allocated(this%calc%halogen))
    write(unit, fmt) "-> electrostatics", abool(allocated(this%calc%coulomb))
    if (allocated(this%calc%coulomb)) then
      write(unit, fmt) "   -> isotropic", abool(allocated(this%calc%coulomb%es2))
      write(unit, fmt) "   -> anisotropic", abool(allocated(this%calc%coulomb%aes2))
      write(unit, fmt) "   -> third-order", abool(allocated(this%calc%coulomb%es3))
    end if
  contains
    pure function abool(cond) result(str)
      logical, intent(in) :: cond
      character(len=merge(3, 2, cond)) :: str
      if (cond) then
        str = "Yes"
      else
        str = "No"
      end if
    end function abool
  #:else
    call notImplementedError
  #:endif
  end subroutine writeTBLiteInfo


  !> Update internal stored coordinates
  subroutine updateCoords(this, env, neighList, img2CentCell, coords, species0)

    !> Data structure
    class(TTBLite), intent(inout) :: this

    !> Computational environment settings
    type(TEnvironment), intent(in) :: env

    !> List of neighbours to atoms
    type(TNeighbourList), intent(in) :: neighList

    !> Image to central cell atom index
    integer, intent(in) :: img2CentCell(:)

    !> Atomic coordinates
    real(dp), intent(in) :: coords(:,:)

    !> Central cell chemical species
    integer, intent(in) :: species0(:)

  #:if WITH_TBLITE
    type(container_cache) :: hcache, rcache

    this%mol%xyz(:, :) = coords(:, :this%mol%nat)
    this%ehal(:) = 0.0_dp
    this%erep(:) = 0.0_dp
    this%edisp(:) = 0.0_dp
    this%gradient(:, :) = 0.0_dp
    this%sigma(:, :) = 0.0_dp

    if (allocated(this%calc%halogen)) then
      call this%calc%halogen%update(this%mol, hcache)
      call this%calc%halogen%get_engrad(this%mol, hcache, this%ehal, &
          & this%gradient, this%sigma)
    end if

    if (allocated(this%calc%repulsion)) then
      call this%calc%repulsion%update(this%mol, rcache)
      call this%calc%repulsion%get_engrad(this%mol, rcache, this%erep, &
          & this%gradient, this%sigma)
    end if

    if (allocated(this%calc%dispersion)) then
      call this%calc%dispersion%update(this%mol, this%dcache)
      call this%calc%dispersion%get_engrad(this%mol, this%dcache, this%edisp, &
          & this%gradient, this%sigma)
    end if

    call new_potential(this%pot, this%mol, this%calc%bas, this%wfn%nspin)
    if (allocated(this%calc%coulomb)) then
      call this%calc%coulomb%update(this%mol, this%cache)
    end if

    if (allocated(this%calc%ncoord)) then
      call this%calc%ncoord%get_cn(this%mol, this%cn, this%dcndr, this%dcndL)
    end if

    call get_selfenergy(this%calc%h0, this%mol%id, this%calc%bas%ish_at, &
        & this%calc%bas%nsh_id, cn=this%cn, selfenergy=this%selfenergy, dsedcn=this%dsedcn)
  #:else
    call notImplementedError
  #:endif
  end subroutine updateCoords


  !> Update internal copy of lattice vectors
  subroutine updateLatVecs(this, latVecs)

    !> Data structure
    class(TTBLite), intent(inout) :: this

    !> Lattice vectors
    real(dp), intent(in) :: latVecs(:,:)

  #:if WITH_TBLITE
    this%mol%lattice(:, :) = latVecs
  #:else
    call notImplementedError
  #:endif
  end subroutine updateLatVecs


  !> Get energy contributions
  subroutine getEnergies(this, energies)

    !> Data structure
    class(TTBLite), intent(inout) :: this

    !> Energy contributions for each atom
    real(dp), intent(out) :: energies(:)

  #:if WITH_TBLITE
    energies(:) = this%ehal + this%erep + this%edisp + this%escd + this%ees
  #:else
    call notImplementedError
  #:endif
  end subroutine getEnergies


  !> Get force contributions
  subroutine addGradients(this, env, neighList, species, coords, img2CentCell, gradients)

    !> Data structure
    class(TTBLite), intent(inout) :: this

    !> Computational environment settings
    type(TEnvironment), intent(in) :: env

    !> Neighbour list.
    type(TNeighbourList), intent(in) :: neighList

    !> Specie for each atom.
    integer, intent(in) :: species(:)

    !> Coordinate of each atom.
    real(dp), intent(in) :: coords(:,:)

    !> Mapping of atoms to cetnral cell.
    integer, intent(in) :: img2CentCell(:)

    !> Gradient contributions for each atom
    real(dp), intent(inout) :: gradients(:,:)

  #:if WITH_TBLITE
    gradients(:, :) = gradients + this%gradient
  #:else
    call notImplementedError
  #:endif
  end subroutine addGradients


  !> Get stress tensor contributions, by converting the saved strain derivatives.
  !> Calculating the stress tensor includes a sign change from the strain derivatives
  !> and a normalization with the cell volume
  subroutine getStress(this, stress)

    !> Data structure
    class(TTBLite), intent(inout) :: this

    !> Stress tensor contributions
    real(dp), intent(out) :: stress(:,:)

  #:if WITH_TBLITE
    stress(:, :) = -this%sigma / abs(determinant33(this%mol%lattice))
  #:else
    call notImplementedError
  #:endif
  end subroutine getStress


  !> Distance cut off for dispersion interactions
  function getRCutoff(this) result(cutoff)

    !> Data structure
    class(TTBLite), intent(inout) :: this

    !> Resulting cutoff
    real(dp) :: cutoff

  #:if WITH_TBLITE
    cutoff = get_cutoff(this%calc%bas)
  #:else
    call notImplementedError
  #:endif
  end function getRCutoff


  !> Updates with changed charges for the instance.
  !>
  !> This routine will be called for both potential and energy calculations.
  !> In case of the energy calculations the orbital charges qq are the actual
  !> output charges, while in case of the potential calculation the orbital
  !> charges in qq are an incomplete output from the mixer and are only accurate
  !> up to the requested charges from the variable_info routine.
  !>
  !> Also the orbital charges have the opposite sign of the ones requested from
  !> the library.
  subroutine updateCharges(this, env, species, neighList, qq, q0, dipAtom, quadAtom, &
      & img2CentCell, orb)

    !> Data structure
    class(TTBLite), intent(inout) :: this

    !> Computational environment settings
    type(TEnvironment), intent(in) :: env

    !> Species, shape: [nAtom]
    integer, intent(in) :: species(:)

    !> Neighbour list.
    type(TNeighbourList), intent(in) :: neighList

    !> Orbital populations
    real(dp), intent(in) :: qq(:,:,:)

    !> Reference orbital populations
    real(dp), intent(in) :: q0(:,:,:)

    !> Cumulative atomic dipole populations
    real(dp), intent(in), optional :: dipAtom(:,:,:)

    !> Cumulative atomic quadrupole populations
    real(dp), intent(in), optional :: quadAtom(:,:,:)

    !> Mapping on atoms in central cell.
    integer, intent(in) :: img2CentCell(:)

    !> Orbital information
    type(TOrbitals), intent(in) :: orb

  #:if WITH_TBLITE
    integer :: iAt, iSh, ii
    real(dp), allocatable :: dQAtom(:), dQShell(:, :)

    call this%pot%reset
    this%escd(:) = 0.0_dp
    this%ees(:) = 0.0_dp

    allocate(dQAtom(this%mol%nat), dQShell(orb%mShell, this%mol%nat))
    call getSummedCharges(species, orb, qq, q0, dQAtom=dQAtom, dQShell=dQShell)

    this%wfn%qat(:, 1) = -dQAtom
    do iAt = 1, size(dQShell, 2)
      ii = this%calc%bas%ish_at(iAt)
      do iSh = 1, this%calc%bas%nsh_at(iAt)
        this%wfn%qsh(ii+iSh, 1) = -dQShell(iSh, iAt)
      end do
    end do

    if (present(dipAtom)) then
      this%wfn%dpat(:, :, 1) = -dipAtom(:, :, 1)
    end if

    if (present(quadAtom)) then
      this%wfn%qpat(:, :, 1) = -quadAtom(:, :, 1)
    end if

    if (allocated(this%calc%coulomb)) then
      call this%calc%coulomb%get_potential(this%mol, this%cache, this%wfn, this%pot)
    end if
    if (allocated(this%calc%dispersion)) then
      call this%calc%dispersion%get_potential(this%mol, this%dcache, this%wfn, this%pot)
    end if

    if (allocated(this%calc%coulomb)) then
      call this%calc%coulomb%get_energy(this%mol, this%cache, this%wfn, this%ees)
    end if
    if (allocated(this%calc%dispersion)) then
      call this%calc%dispersion%get_energy(this%mol, this%dcache, this%wfn, this%escd)
    end if
  #:else
    call notImplementedError
  #:endif
  end subroutine updateCharges


  !> Returns atom-resolved and shell-resolved shifts from library.
  !>
  !> Potential shifts in tblite are calculated as derivatives w.r.t. the partial
  !> charges / multipole moments, while in DFTB+ they are calculated as derivative
  !> w.r.t. the population. This means we have to flip the sign here.
  subroutine getShifts(this, shiftPerAtom, shiftPerShell, dipShift, quadShift)

    !> Data structure
    class(TTBLite), intent(inout) :: this

    !> Shift per atom
    real(dp), intent(out) :: shiftPerAtom(:)

    !> Shift per shell
    real(dp), intent(out) :: shiftPerShell(:,:)

    !> Dipolar shift per atom
    real(dp), intent(out), optional :: dipShift(:,:)

    !> Quadrupolar shift per atom
    real(dp), intent(out), optional :: quadShift(:,:)

  #:if WITH_TBLITE
    integer :: iAt, iSh, ii

    shiftPerAtom(:) = -this%pot%vat(:, 1)

    shiftPerShell(:,:) = 0.0_dp
    do iAt = 1, size(shiftPerShell, 2)
      ii = this%calc%bas%ish_at(iAt)
      do iSh = 1, this%calc%bas%nsh_at(iAt)
        shiftPerShell(iSh, iAt) = -this%pot%vsh(ii+iSh, 1)
      end do
    end do

    if (present(dipShift)) then
      dipShift(:,:) = -this%pot%vdp(:, :, 1)
    end if

    if (present(quadShift)) then
      quadshift(:,:) = -this%pot%vqp(:, :, 1)
    end if
  #:else
    call notImplementedError
  #:endif
  end subroutine getShifts


  !> Create orbital information. The orbital information is already generated in
  !> this%calc%bas, but might be incomplete w.r.t. the information required here.
  subroutine getOrbitalInfo(this, species0, orb)

    !> Data structure
    class(TTBLite), intent(in) :: this

    !> Species of each atom, shape: [nAtom]
    integer, intent(in) :: species0(:)

    !> Orbital information
    type(TOrbitals), intent(out) :: orb

  #:if WITH_TBLITE
    call setupOrbitalInfo(this%calc%bas, this%mol%id, species0, orb)
  #:else
    call notImplementedError
  #:endif
  end subroutine getOrbitalInfo


#:if WITH_TBLITE
  !> Create orbital information from tight-binding basis set.
  !>
  !> Both this%mol%id and species are identifiers for unique groups of elements,
  !> since we cannot guarantee that the species in DFTB+ are identical to the ones
  !> found in the library we will reconstruct the species information from the atoms.
  subroutine setupOrbitalInfo(bas, id, species0, orb)

    !> Basis set data
    class(basis_type), intent(in) :: bas

    !> Identifier of each atom, shape: [nAtom]
    integer, intent(in) :: id(:)

    !> Species of each atom, shape: [nAtom]
    integer, intent(in) :: species0(:)

    !> Orbital information
    type(TOrbitals), intent(out) :: orb

    integer :: nShell, mShell, nSpecies, iAt, iSp, iId, iSh, ind, ii
    logical, allocatable :: done(:)

    nSpecies = maxval(species0)
    mShell = maxval(bas%nsh_id)
    allocate(done(nSpecies))
    allocate(orb%angShell(mShell, nSpecies))
    done(:) = .false.
    do iAt = 1, size(species0)
      iId = id(iAt)
      iSp = species0(iAt)
      if (done(iSp)) cycle
      do iSh = 1, bas%nsh_at(iAt)
        orb%angShell(iSh, iSp) = bas%cgto(iSh, iId)%ang
      end do
      done(iSp) = .true.
    end do

    allocate(orb%nShell(nSpecies))
    allocate(orb%nOrbSpecies(nSpecies))
    allocate(orb%nOrbAtom(size(species0)))

    done(:) = .false.
    do iAt = 1, size(species0)
      iSp = species0(iAt)
      ii = bas%ish_at(iAt)
      nShell = bas%nsh_at(iAt)
      orb%nOrbAtom(iAt) = sum(bas%nao_sh(ii+1:ii+nShell))
      if (done(iSp)) cycle
      orb%nOrbSpecies(iSp) = orb%nOrbAtom(iAt)
      orb%nShell(iSp) = nShell
      done(iSp) = .true.
    end do

    orb%mShell = mShell
    orb%mOrb = maxval(orb%nOrbSpecies)
    orb%nOrb = sum(orb%nOrbAtom)

    allocate(orb%iShellOrb(orb%mOrb, nSpecies))
    allocate(orb%posShell(orb%mShell+1, nSpecies))
    do iSp = 1, nSpecies
      ind = 1
      do iSh = 1, orb%nShell(iSp)
        orb%posShell(iSh, iSp) = ind
        orb%iShellOrb(ind:ind+2*orb%angShell(iSh, iSp), iSp) = iSh
        ind = ind + 2 * orb%angShell(iSh, iSp) + 1
      end do
      orb%posShell(orb%nShell(iSp)+1, iSp) = ind
    end do
  end subroutine setupOrbitalInfo
#:endif


  !> Get information on required multipolar contributions
  subroutine getMultipoleInfo(this, nDipole, nQuadrupole)

    !> Data structure
    class(TTBLite), intent(in) :: this

    !> Number of dipole moment components
    integer, intent(out) :: nDipole

    !> Number of quadrupole moment components
    integer, intent(out) :: nQuadrupole

  #:if WITH_TBLITE
    nDipole = dimDipole
    nQuadrupole = dimQuadrupole
  #:else
    call notImplementedError
  #:endif
  end subroutine getMultipoleInfo


  !> Get reference occupation numbers.
  subroutine getReferenceN0(this, species0, referenceN0)

    !> Data structure
    class(TTBLite), intent(in) :: this

    !> Species of each atom, shape: [nAtom]
    integer, intent(in) :: species0(:)

    !> Reference occupation numbers
    real(dp), intent(out) :: referenceN0(:, :)

  #:if WITH_TBLITE
    integer :: iAt, iSp, iId, iSh
    logical, allocatable :: done(:)

    referenceN0(:,:) = 0.0_dp
    allocate(done(maxval(species0)))
    done(:) = .false.
    do iAt = 1, size(species0)
      iId = this%mol%id(iAt)
      iSp = species0(iAt)
      if (done(iSp)) cycle
      do iSh = 1, this%calc%bas%nsh_at(iAt)
        referenceN0(iSh, iSp) = this%calc%h0%refocc(iSh, iId)
      end do
      done(iSp) = .true.
    end do
  #:else
    call notImplementedError
  #:endif
  end subroutine getReferenceN0


  !> Returns the equivalence to get the correct mixing of charge dependent contributions
  subroutine getOrbitalEquiv(this, equivOrb, equivDip, equivQuad, orb, species)

    !> Data structure
    class(TTBLite), intent(inout) :: this

    !> The equivalence vector for orbital populations
    integer, intent(out) :: equivOrb(:,:,:)

    !> The equivalence vector for cumulative atomic dipole populations
    integer, intent(out) :: equivDip(:,:)

    !> The equivalence vector for cumulative atomic quadrupole populations
    integer, intent(out) :: equivQuad(:,:)

    !> Information about the orbitals and their angular momenta
    type(TOrbitals), intent(in) :: orb

    !> Species of each atom
    integer, intent(in) :: species(:)

  #:if WITH_TBLITE
    type(scf_info) :: info
    integer :: nAtom, iCount, iSpin, nSpin, iAt, iSp, ii

    nAtom = size(equivOrb, dim=2)
    nSpin = size(equivOrb, dim=3)

    equivOrb(:,:,:) = 0
    equivDip(:,:) = 0
    equivQuad(:,:) = 0

    info = this%calc%variable_info()
    select case(info%charge)
    case(atom_resolved)
      iCount = 0
      do iSpin = 1, nSpin
        do iAt = 1, nAtom
          iSp = species(iAt)
          iCount = iCount + 1
          do ii = 1, orb%nOrbSpecies(iSp)
            equivOrb(ii, iAt, iSpin) = iCount
          end do
        end do
      end do
    case(shell_resolved)
      iCount = 0
      do iSpin = 1, nSpin
        do iAt = 1, nAtom
          iSp = species(iAt)
          do ii = 1, orb%nOrbSpecies(iSp)
            equivOrb(ii, iAt, iSpin) = iCount + orb%iShellOrb(ii, iSp)
          end do
          iCount = iCount + orb%nShell(iSp)
        end do
      end do
    case(orbital_resolved)
      iCount = 0
      do iSpin = 1, nSpin
        do iAt = 1, nAtom
          iSp = species(iAt)
          do ii = 1, orb%nOrbSpecies(iSp)
            iCount = iCount + 1
            equivOrb(ii, iAt, iSpin) = iCount
          end do
        end do
      end do
    end select

    select case(info%dipole)
    case(atom_resolved)
      iCount = 0
      do iAt = 1, nAtom
        do ii = 1, dimDipole
          iCount = iCount + 1
          equivDip(ii, iAt) = iCount
        end do
      end do
    end select

    select case(info%quadrupole)
    case(atom_resolved)
      iCount = 0
      do iAt = 1, nAtom
        do ii = 1, dimQuadrupole
          iCount = iCount + 1
          equivQuad(ii, iAt) = iCount
        end do
      end do
    end select
  #:else
    call notImplementedError
  #:endif
  end subroutine getOrbitalEquiv


  !> Return Hubbard parameters from extended tight binding model
  subroutine getHubbardU(this, hubbU)

    !> Data structure
    class(TTBLite), intent(in) :: this

    !> Hubbard parameters for
    real(dp), intent(out) :: hubbU(:, :)

  #:if WITH_TBLITE
    hubbU(:, :) = 0.0_dp
    if (.not.allocated(this%calc%coulomb)) return
    if (.not.allocated(this%calc%coulomb%es2)) return

    block
      use tblite_coulomb_charge, only : gamma_coulomb, effective_coulomb
      integer :: iSp, iSh
      select type(es2 => this%calc%coulomb%es2)
      class is(gamma_coulomb)
        do iSp = 1, size(hubbU, 2)
          hubbU(:, iSp) = es2%hubbard(:, this%sp2id(iSp))
        end do
      class is(effective_coulomb)
        do iSp = 1, size(hubbU, 2)
          do iSh = 1, size(hubbU, 1)
            hubbU(iSh, iSp) = es2%hubbard(iSh, iSh, this%sp2id(iSp), this%sp2id(iSp))
          end do
        end do
      end select
    end block
  #:else
    call notImplementedError
  #:endif
  end subroutine getHubbardU

  !> Remove second order electrostatics from container
  subroutine removeES2(this)

    !> Data structure
    class(TTBLite), intent(inout) :: this

  #:if WITH_TBLITE
    if (.not.allocated(this%calc%coulomb)) return
    if (.not.allocated(this%calc%coulomb%es2)) return
    deallocate(this%calc%coulomb%es2)
  #:else
    call notImplementedError
  #:endif
  end subroutine removeES2

  !> Build atomic block sparse compressed Hamiltonian and overlap related integrals
  subroutine buildSH0(this, env, species, coords, nNeighbour, iNeighbours, img2centCell, &
      & iPair, orb, hamiltonian, overlap, dpintBra, dpintKet, qpintBra, qpintKet)

    !> Data structure
    class(TTBLite), intent(inout) :: this

    !> Computational environment settings
    type(TEnvironment), intent(in) :: env

    !> Returns the non-self-consistent Hamiltonian
    real(dp), intent(out) :: hamiltonian(:)

    !> Returns the non-self-consistent Hamiltonian
    real(dp), intent(out) :: overlap(:)

    !> Dipole moment integral matrix, operator on the bra function
    real(dp), intent(inout) :: dpintBra(:, :)

    !> Dipole moment integral matrix, operator on the ket function
    real(dp), intent(inout) :: dpintKet(:, :)

    !> Quadrupole moment integral matrix, operator on the bra function
    real(dp), intent(inout) :: qpintBra(:, :)

    !> Quadrupole moment integral matrix, operator on the ket function
    real(dp), intent(inout) :: qpintKet(:, :)

    !> Atomic coordinates
    real(dp), intent(in) :: coords(:,:)

    !> Number of surrounding neighbours for each atom
    integer, intent(in) :: nNeighbour(:)

    !> List of surrounding neighbours for each atom
    integer, intent(in) :: iNeighbours(0:,:)

    !> Mapping of images back to atoms in the central cell
    integer, intent(in) :: img2centCell(:)

    !> Chemical species of each atom
    integer, intent(in) :: species(:)

    !> Shift vector, where the interaction between two atoms
    integer, intent(in) :: iPair(0:,:)

    !> Information about the orbitals in the system
    type(TOrbitals), intent(in) :: orb

  #:if WITH_TBLITE
    integer :: nAtom, iAtFirst, iAtLast

    nAtom = size(nNeighbour)
    hamiltonian(:) = 0.0_dp
    overlap(:) = 0.0_dp
    dpintBra(:, :) = 0.0_dp
    dpintKet(:, :) = 0.0_dp
    qpintBra(:, :) = 0.0_dp
    qpintKet(:, :) = 0.0_dp

    call distributeRangeInChunks(env, 1, nAtom, iAtFirst, iAtLast)

    call buildDiagonalBlocks(iAtFirst, iAtLast, this%mol%id, iPair, orb%nOrbAtom, &
        & this%calc%bas, this%selfenergy, hamiltonian, overlap, dpintBra, dpintKet, &
        & qpintBra, qpintKet)

    call buildDiatomicBlocks(iAtFirst, iAtLast, this%mol%id, coords, &
        & nNeighbour, iNeighbours, img2centCell, iPair, orb%nOrbAtom, &
        & this%calc%bas, this%calc%h0, this%selfenergy, hamiltonian, overlap, &
        & dpintBra, dpintKet, qpintBra, qpintKet)

    call assembleChunks(env, hamiltonian)
    call assembleChunks(env, overlap)
    call assembleChunks(env, dpintBra)
    call assembleChunks(env, dpintKet)
    call assembleChunks(env, qpintBra)
    call assembleChunks(env, qpintKet)
  #:else
    call notImplementedError
  #:endif
  end subroutine buildSH0


#:if WITH_TBLITE
  !> Put the on-site energies into the Hamiltonian,
  !> and <lm|l'm'> = delta_l,l' * delta_m,m' for the overlap
  !>
  !> This routine requires access to the internals of the tblite library,
  !> breaking the abstraction layer of the library interface.
  !> Candidate for (partial) upstreaming in tblite library.
  subroutine buildDiagonalBlocks(iAtFirst, iAtLast, species, iPair, nOrbAtom, &
      & bas, selfenergy, hamiltonian, overlap, dpintBra, dpintKet, qpintBra, qpintKet)

    !> Atom range for this processor to evaluate
    integer, intent(in) :: iAtFirst, iAtLast

    !> Chemical species of each atom
    integer, intent(in) :: species(:)

    !> Shift vector, where the interaction between two atoms starts in the sparse format.
    integer, intent(in) :: iPair(0:,:)

    !> Size of the block in spare format for each atom
    integer, intent(in) :: nOrbAtom(:)

    !> Basis set information
    type(basis_type), intent(in) :: bas

    !> Diagonal elements of the Hamiltonian
    real(dp), intent(in) :: selfenergy(:)

    !> Overlap integral matrix
    real(dp), intent(inout) :: overlap(:)

    !> Dipole moment integral matrix, operator on the bra function
    real(dp), intent(inout) :: dpintBra(:, :)

    !> Dipole moment integral matrix, operator on the ket function
    real(dp), intent(inout) :: dpintKet(:, :)

    !> Quadrupole moment integral matrix, operator on the bra function
    real(dp), intent(inout) :: qpintBra(:, :)

    !> Quadrupole moment integral matrix, operator on the ket function
    real(dp), intent(inout) :: qpintKet(:, :)

    !> Effective Hamiltonian
    real(dp), intent(inout) :: hamiltonian(:)

    integer :: iAt, iZp, ind, iOrb, ii, jj, nBlk, is, io, ish, jsh, iao, jao, ij, iblk, nao
    integer :: li, lj
    real(dp), parameter :: r2 = 0.0_dp, vec(3) = 0.0_dp
    real(dp), allocatable :: stmp(:)
    real(dp), allocatable :: dtmp(:, :), qtmp(:, :)

    allocate(stmp(msao(bas%maxl)**2))
    allocate(dtmp(dimDipole, msao(bas%maxl)**2), qtmp(dimQuadrupole, msao(bas%maxl)**2))

    !$omp parallel do schedule(runtime) default(none) &
    !$omp shared(iAtFirst, iAtLast, species, bas, iPair, nOrbAtom, selfenergy, &
    !$omp& overlap, hamiltonian, dpintBra, dpintKet, qpintBra, qpintKet) &
    !$omp private(iAt, iZp, is, io, ind, ish, jsh, ii, jj, iao, jao, nBlk, ij, iblk, nao, &
    !$omp& stmp, dtmp, qtmp, li, lj)
    do iAt = iAtFirst, iAtLast
      iZp = species(iAt)
      is = bas%ish_at(iAt)
      io = bas%iao_sh(is+1)
      ind = iPair(0, iAt)
      nBlk = nOrbAtom(iAt)
      do iSh = 1, bas%nsh_id(iZp)
        ii = bas%iao_sh(is+iSh) - io
        do iao = 1, msao(bas%cgto(iSh, iZp)%ang)
          overlap(ind + ii+iao + nBlk*(ii+iao-1)) = & !overlap(ii+iao, ii+iao) &
            1.0_dp

          hamiltonian(ind + ii+iao + nblk*(ii+iao-1)) = & !hamiltonian(ii+iao, ii+iao) &
            selfenergy(is+iSh)
        end do
      end do

      do ish = 1, bas%nsh_id(iZp)
        ii = bas%iao_sh(is+iSh) - io
        do jsh = 1, bas%nsh_id(iZp)
          jj = bas%iao_sh(is+jSh) - io
          call multipole_cgto(bas%cgto(jSh, iZp), bas%cgto(iSh, iZp), &
              & r2, vec, bas%intcut, stmp, dtmp, qtmp)
          !call overlap_cgto(bas%cgto(jSh, iZp), bas%cgto(iSh, iZp), &
          !    & r2, vec, bas%intcut, stmp)

          li = bas%cgto(iSh, iZp)%ang
          lj = bas%cgto(jSh, iZp)%ang
          nao = msao(lj)
          do iao = 1, msao(li)
            do jao = 1, nao
              ij = jao + nao*(iao-1)
              iblk = ind + jj+jao + nBlk*(ii+iao-1)

              dpintBra(:, iblk) = dtmp(:, ij)
              dpintKet(:, iblk) = dtmp(:, ij)

              qpintBra(:, iblk) = qtmp(:, ij)
              qpintKet(:, iblk) = qtmp(:, ij)
            end do
          end do

        end do
      end do

    end do
  end subroutine buildDiagonalBlocks


  !> The Hamiltonian is saved in an atomic block sparse compressed format.
  !> We calculate a shell pair as a contiguous blocks and spread it to the
  !> contiguous atomic block.
  !>
  !> Candidate for (partial) upstreaming in tblite library.
  subroutine buildDiatomicBlocks(iAtFirst, iAtLast, species, coords, &
      & nNeighbour, iNeighbours, img2centCell, iPair, nOrbAtom, &
      & bas, h0, selfenergy, hamiltonian, overlap, dpintBra, dpintKet, &
      & qpintBra, qpintKet)

    !> Atom range for this processor to evaluate
    integer, intent(in) :: iAtFirst, iAtLast

    !> Chemical species of each atom
    integer, intent(in) :: species(:)

    !> Atomic coordinates
    real(dp), intent(in) :: coords(:,:)

    !> Number of surrounding neighbours for each atom
    integer, intent(in) :: nNeighbour(:)

    !> List of surrounding neighbours for each atom
    integer, intent(in) :: iNeighbours(0:,:)

    !> Shift vector, where the interaction between two atoms starts in the sparse format.
    integer, intent(in) :: iPair(0:,:)

    !> Size of the block in spare format for each atom
    integer, intent(in) :: nOrbAtom(:)

    !> Mapping of images back to atoms in the central cell
    integer, intent(in) :: img2centCell(:)

    !> Basis set information
    type(basis_type), intent(in) :: bas

    !> Hamiltonian interaction data
    type(tb_hamiltonian), intent(in) :: h0

    !> Diagonal elements of the Hamiltonian
    real(dp), intent(in) :: selfenergy(:)

    !> Overlap integral matrix
    real(dp), intent(inout) :: overlap(:)

    !> Dipole moment integral matrix, operator on the bra function
    real(dp), intent(inout) :: dpintBra(:, :)

    !> Dipole moment integral matrix, operator on the ket function
    real(dp), intent(inout) :: dpintKet(:, :)

    !> Quadrupole moment integral matrix, operator on the bra function
    real(dp), intent(inout) :: qpintBra(:, :)

    !> Quadrupole moment integral matrix, operator on the ket function
    real(dp), intent(inout) :: qpintKet(:, :)

    !> Effective Hamiltonian
    real(dp), intent(inout) :: hamiltonian(:)

    integer :: iAt, jAt, iZp, jZp, iNeigh, img, ind, io, jo, iblk, ij, li, lj
    integer :: iSh, jSh, is, js, ii, jj, iao, jao, nao, nBlk
    real(dp) :: rr, r2, vec(3), hij, shpoly, dtmpj(dimDipole), qtmpj(dimQuadrupole)
    real(dp), allocatable :: stmp(:)
    real(dp), allocatable :: dtmpi(:, :), qtmpi(:, :)

    allocate(stmp(msao(bas%maxl)**2))
    allocate(dtmpi(dimDipole, msao(bas%maxl)**2), qtmpi(dimQuadrupole, msao(bas%maxl)**2))

    !$omp parallel do schedule(runtime) default(none) &
    !$omp shared(iatfirst, iatlast, species, bas, overlap, hamiltonian, h0, selfenergy, &
    !$omp& nNeighbour, ineighbours, img2centCell, ipair, nOrbAtom, coords, &
    !$omp& dpintBra, dpintKet, qpintBra, qpintKet) &
    !$omp private(iAt, jAt, iZp, jZp, is, js, iSh, jSh, ii, jj, iao, jao, nao, nBlk, &
    !$omp& io, jo, iNeigh, img, ind, r2, vec, stmp, dtmpj, dtmpi, qtmpj, qtmpi, hij, &
    !$omp& shpoly, rr, iblk, ij, li, lj)
    do iAt = iAtFirst, iAtLast
      iZp = species(iAt)
      is = bas%ish_at(iAt)
      io = bas%iao_sh(is+1)
      do iNeigh = 1, nNeighbour(iAt)
        img = iNeighbours(iNeigh, iAt)
        jAt = img2centCell(img)
        jZp = species(jAt)
        js = bas%ish_at(jAt)
        ind = iPair(iNeigh, iAt)
        jo = bas%iao_sh(js+1)
        nBlk = nOrbAtom(jAt)

        vec(:) = coords(:, iAt) - coords(:, img)
        r2 = vec(1)**2 + vec(2)**2 + vec(3)**2
        rr = sqrt(sqrt(r2) / (h0%rad(jZp) + h0%rad(iZp)))
        do ish = 1, bas%nsh_id(iZp)
          ii = bas%iao_sh(is+iSh) - io
          do jsh = 1, bas%nsh_id(jZp)
            jj = bas%iao_sh(js+jSh) - jo
            call multipole_cgto(bas%cgto(jSh, jZp), bas%cgto(iSh, iZp), &
                & r2, vec, bas%intcut, stmp, dtmpi, qtmpi)
            !call overlap_cgto(bas%cgto(jSh, jZp), bas%cgto(iSh, iZp), &
            !    & r2, vec, bas%intcut, stmp)

            shpoly = (1.0_dp + h0%shpoly(iSh, iZp)*rr) &
              & * (1.0_dp + h0%shpoly(jSh, jZp)*rr)

            hij = 0.5_dp * (selfenergy(is+iSh) + selfenergy(js+jSh)) &
              * h0%hscale(jSh, iSh, jZp, iZp) * shpoly


            li = bas%cgto(iSh, iZp)%ang
            lj = bas%cgto(jSh, jZp)%ang
            nao = msao(lj)
            do iao = 1, msao(li)
              do jao = 1, nao
                ij = jao + nao*(iao-1)
                iblk = ind + jj+jao + nBlk*(ii+iao-1)
                call shiftOperator(vec, stmp(ij), dtmpi(:, ij), qtmpi(:, ij), dtmpj, qtmpj)

                overlap(iblk) = stmp(ij)

                dpintBra(:, iblk) = dtmpj
                dpintKet(:, iblk) = dtmpi(:, ij)

                qpintBra(:, iblk) = qtmpj
                qpintKet(:, iblk) = qtmpi(:, ij)

                hamiltonian(iblk) = stmp(ij) * hij
              end do
            end do

          end do
        end do

      end do
    end do
  end subroutine buildDiatomicBlocks
#:endif


  !> Shift multipole operator from Ket function (center i) to Bra function (center j),
  !> the multipole operator on the Bra function can be assembled from the lower moments
  !> on the Ket function and the displacement vector using horizontal shift rules.
  !>
  !> This is usually done inside the tblite library, but since we want to have both
  !> Bra and Ket contributions at once and do not want to iterate over both triangles
  !> of the multipole integral matrix we perform the shift of the moment operator here.
  !>
  !> Candidate for (partial) upstreaming in tblite library.
  pure subroutine shiftOperator(vec, s, di, qi, dj, qj)

    !> Displacement vector of center i and j
    real(dp),intent(in) :: vec(:)

    !> Overlap integral between basis functions
    real(dp),intent(in) :: s

    !> Dipole integral with operator on Ket function (center i)
    real(dp),intent(in) :: di(:)

    !> Quadrupole integral with operator on Ket function (center i)
    real(dp),intent(in) :: qi(:)

    !> Dipole integral with operator on Bra function (center j)
    real(dp),intent(out) :: dj(:)

    !> Quadrupole integral with operator on Bra function (center j)
    real(dp),intent(out) :: qj(:)

    real(dp) :: tr

    ! Create dipole operator on Bra function from Ket function and shift contribution
    ! due to monopol displacement
    dj(1) = di(1) + vec(1)*s
    dj(2) = di(2) + vec(2)*s
    dj(3) = di(3) + vec(3)*s

    ! For the quadrupole operator on the Bra function we first construct the shift
    ! contribution from the dipole and monopol displacement, since we have to remove
    ! the trace contribution from the shift and the moment integral on the Ket function
    ! is already traceless
    qj(1) = 2*vec(1)*di(1) + vec(1)**2*s
    qj(3) = 2*vec(2)*di(2) + vec(2)**2*s
    qj(6) = 2*vec(3)*di(3) + vec(3)**2*s
    qj(2) = vec(1)*di(2) + vec(2)*di(1) + vec(1)*vec(2)*s
    qj(4) = vec(1)*di(3) + vec(3)*di(1) + vec(1)*vec(3)*s
    qj(5) = vec(2)*di(3) + vec(3)*di(2) + vec(2)*vec(3)*s
    ! Now collect the trace of the shift contribution
    tr = 0.5_dp * (qj(1) + qj(3) + qj(6))

    ! Finally, assemble the quadrupole operator on the Bra function from the operator
    ! on the Ket function and the traceless shift contribution
    qj(1) = qi(1) + 1.5_dp * qj(1) - tr
    qj(2) = qi(2) + 1.5_dp * qj(2)
    qj(3) = qi(3) + 1.5_dp * qj(3) - tr
    qj(4) = qi(4) + 1.5_dp * qj(4)
    qj(5) = qi(5) + 1.5_dp * qj(5)
    qj(6) = qi(6) + 1.5_dp * qj(6) - tr
  end subroutine shiftOperator


  !> Evaluate derivatives of potential shifts from Hamiltonian and overlap related
  !> integrals.
  subroutine buildDerivativeShift(this, env, DM, EDM, coords, species, &
      & nNeighbour, iNeighbours, img2CentCell, iPair, orb, shift, dipShift, quadShift)

    !> Data structure
    class(TTBLite), intent(inout) :: this

    !> Computational environment settings
    type(TEnvironment), intent(in) :: env

    !> density matrix in packed format
    real(dp), intent(in) :: DM(:,:)

    !> energy-weighted density matrix in packed format
    real(dp), intent(in) :: EDM(:)

    !> list of all atomic coordinates
    real(dp), intent(in) :: coords(:,:)

    !> list of all atomic species
    integer, intent(in) :: species(:)

    !> number of neighbours of each atom
    integer, intent(in) :: nNeighbour(:)

    !> neighbour list for atoms
    integer, intent(in) :: iNeighbours(0:,:)

    !> indexing array for periodic image atoms
    integer, intent(in) :: img2CentCell(:)

    !> indexing array for the Hamiltonian
    integer, intent(in) :: iPair(0:,:)

    !> Information about the shells and orbitals in the system.
    type(TOrbitals), intent(in) :: orb

    !> block shift from the potential
    real(dp), intent(in) :: shift(:,:,:,:)

    !> Dipole potential shift, shape: [nDipole, nAtom]
    real(dp), intent(in) :: dipShift(:, :)

    !> Quadrupole potential shift, shape: [nQuadrupole, nAtom]
    real(dp), intent(in) :: quadShift(:, :)

  #:if WITH_TBLITE
    integer :: nAtom, iAtFirst, iAtLast
    real(dp), allocatable :: gradient(:, :), sigma(:, :), dEdcn(:)

    if (allocated(this%calc%coulomb)) then
      call this%calc%coulomb%get_gradient(this%mol, this%cache, this%wfn, this%gradient, &
          & this%sigma)
    end if

    if (allocated(this%calc%dispersion)) then
      call this%calc%dispersion%get_gradient(this%mol, this%dcache, this%wfn, this%gradient, &
          & this%sigma)
    end if

    nAtom = size(nNeighbour)
    allocate(gradient(3, nAtom), sigma(3, 3), dEdcn(nAtom))
    gradient(:, :) = 0.0_dp
    sigma(:, :) = 0.0_dp
    dEdcn(:) = 0.0_dp

    call distributeRangeInChunks(env, 1, nAtom, iAtFirst, iAtLast)

    call buildDiagonalDerivs(iAtFirst, iAtLast, this%mol%id, iPair, orb%nOrbAtom, &
        & this%calc%bas, this%dsedcn, DM, dEdcn)
    call buildDiatomicDerivs(iAtFirst, iAtLast, this%mol%id, coords, nNeighbour, &
        & iNeighbours, img2centCell, iPair, orb%nOrbAtom, this%calc%bas, this%calc%h0, &
        & this%selfenergy, this%dsedcn, DM, EDM, shift, dipShift, quadShift, &
        & dEdcn, gradient, sigma)

    call assembleChunks(env, dEdcn)
    call assembleChunks(env, gradient)
    call assembleChunks(env, sigma)

    call gemv(gradient, this%dcndr, dEdcn, beta=1.0_dp)
    call gemv(sigma, this%dcndL, dEdcn, beta=1.0_dp)

    this%gradient(:, :) = this%gradient + gradient
    this%sigma(:, :) = this%sigma + sigma
  #:else
    call notImplementedError
  #:endif
  end subroutine buildDerivativeShift


#:if WITH_TBLITE
  !> Calculate derivatives for all diagonal elements.
  !>
  !> Candidate for (partial) upstreaming in tblite library.
  subroutine buildDiagonalDerivs(iAtFirst, iAtLast, species, iPair, nOrbAtom, &
      & bas, dsedcn, pmat, dEdcn)

    !> Atom range for this processor to evaluate
    integer, intent(in) :: iAtFirst, iAtLast

    !> Chemical species of each atom
    integer, intent(in) :: species(:)

    !> Shift vector, where the interaction between two atoms starts in the sparse format.
    integer, intent(in) :: iPair(0:,:)

    !> Size of the block in spare format for each atom
    integer, intent(in) :: nOrbAtom(:)

    !> Basis set information
    type(basis_type), intent(in) :: bas

    !> Derivative of diagonal elements of the Hamiltonian w.r.t. coordination number
    real(dp), intent(in) :: dsedcn(:)

    !> Density matrix
    real(dp), intent(in) :: pmat(:, :)

    !> Derivative of energy w.r.t. coordination number
    real(dp), intent(inout) :: dEdcn(:)

    integer :: iat, izp, is, ish, ii, iao, io, ind, nBlk
    real(dp) :: dhdcni, dcni

    !$omp parallel do schedule(runtime) default(none) reduction(+:dEdcn) &
    !$omp shared(iAtFirst, iAtLast, species, bas, dsedcn, pmat, iPair, nOrbAtom) &
    !$omp private(iat, izp, is, ish, ii, iao, io, ind, nBlk, dcni, dhdcni)
    do iat = iAtFirst, iAtLast
      izp = species(iat)
      is = bas%ish_at(iat)
      io = bas%iao_sh(is+1)
      ind = iPair(0, iAt)
      nBlk = nOrbAtom(iAt)
      do ish = 1, bas%nsh_id(izp)
        ii = bas%iao_sh(is+ish) - io
        dhdcni = dsedcn(is+ish)
        dcni = 0.0_dp
        do iao = 1, msao(bas%cgto(ish, izp)%ang)
          dcni = dcni + dhdcni * pmat(ind + ii+iao + nBlk*(ii+iao-1), 1)
        end do
        dEdcn(iat) = dEdcn(iat) + dcni
      end do
    end do
  end subroutine buildDiagonalDerivs


  !> Calculate derivatives for all diatomic pairs.
  !>
  !> Todo: Handle actual block potential shifts
  !>
  !> Candidate for (partial) upstreaming in tblite library.
  subroutine buildDiatomicDerivs(iAtFirst, iAtLast, species, coords, &
      & nNeighbour, iNeighbours, img2centCell, iPair, nOrbAtom, bas, h0, &
      & selfenergy, dsedcn, pmat, xmat, shift, dipShift, quadShift, &
      & dEdcn, gradient, sigma)

    !> Atom range for this processor to evaluate
    integer, intent(in) :: iAtFirst, iAtLast

    !> Chemical species of each atom
    integer, intent(in) :: species(:)

    !> Atomic coordinates
    real(dp), intent(in) :: coords(:,:)

    !> Number of surrounding neighbours for each atom
    integer, intent(in) :: nNeighbour(:)

    !> List of surrounding neighbours for each atom
    integer, intent(in) :: iNeighbours(0:,:)

    !> Shift vector, where the interaction between two atoms starts in the sparse format.
    integer, intent(in) :: iPair(0:,:)

    !> Size of the block in spare format for each atom
    integer, intent(in) :: nOrbAtom(:)

    !> Mapping of images back to atoms in the central cell
    integer, intent(in) :: img2centCell(:)

    !> Basis set information
    type(basis_type), intent(in) :: bas

    !> Hamiltonian interaction data
    type(tb_hamiltonian), intent(in) :: h0

    !> Diagonal elements of the Hamiltonian
    real(dp), intent(in) :: selfenergy(:)

    !> Derivative of diagonal elements of the Hamiltonian w.r.t. coordination number
    real(dp), intent(in) :: dsedcn(:)

    !> Density matrix in packed form with spin information
    real(dp), intent(in) :: pmat(:, :)

    !> Energy weighted density matrix in packed form
    real(dp), intent(in) :: xmat(:)

    !> Block potential shift, shape: [nOrb, nOrb, nAtom, nSpin]
    real(dp), intent(in) :: shift(:, :, :, :)

    !> Dipole potential shift, shape: [nDipole, nAtom]
    real(dp), intent(in) :: dipShift(:, :)

    !> Quadrupole potential shift, shape: [nQuadrupole, nAtom]
    real(dp), intent(in) :: quadShift(:, :)

    !> Derivative of energy w.r.t. coordination number
    real(dp), intent(inout) :: dEdcn(:)

    !> Energy derivative w.r.t. atomic displacements
    real(dp), intent(inout) :: gradient(:, :)

    !> Energy derivative w.r.t. strain deformations
    real(dp), intent(inout) :: sigma(:, :)

    integer :: iat, jat, izp, jzp, itr, io, jo, nblk, iNeigh, img, ind, li, lj
    integer :: ish, jsh, is, js, ii, jj, iao, jao, nao, ij, iSpin, nSpin, iblk
    real(dp) :: rr, r2, vec(3), cutoff2, hij, shpoly, dshpoly, dG(3), hscale
    real(dp) :: sval, dcni, dcnj, dhdcni, dhdcnj, hpij, pij
    real(dp), allocatable :: stmp(:), dtmp(:, :), qtmp(:, :)
    real(dp), allocatable :: dstmp(:, :), ddtmpi(:, :, :), dqtmpi(:, :, :)
    real(dp), allocatable :: ddtmpj(:, :, :), dqtmpj(:, :, :)

    allocate(stmp(msao(bas%maxl)**2), dstmp(3, msao(bas%maxl)**2), &
      & dtmp(dimDipole, msao(bas%maxl)**2), ddtmpi(3, dimDipole, msao(bas%maxl)**2), &
      & qtmp(dimQuadrupole, msao(bas%maxl)**2), dqtmpi(3, dimQuadrupole, msao(bas%maxl)**2), &
      & ddtmpj(3, dimDipole, msao(bas%maxl)**2), dqtmpj(3, dimQuadrupole, msao(bas%maxl)**2))

    nSpin = size(shift, 4)

    !$omp parallel do schedule(runtime) default(none) reduction(+:dEdcn, gradient, sigma) &
    !$omp shared(iAtFirst, iAtLast, species, bas, h0, selfenergy, dsedcn, coords, &
    !$omp& nNeighbour, iNeighbours, img2centCell, iPair, nOrbAtom, pmat, xmat, shift, &
    !$omp& dipShift, quadShift, nSpin) &
    !$omp private(iat, jat, izp, jzp, is, js, ish, jsh, ii, jj, iao, jao, nao, ij, ind, &
    !$omp& iNeigh, io, jo, nblk, img, r2, vec, stmp, dtmp, qtmp, dstmp, ddtmpi, ddtmpj, &
    !$omp& dqtmpi, dqtmpj, hij, shpoly, dshpoly, dG, dcni, dcnj, dhdcni, dhdcnj, hpij, rr, &
    !$omp& sval, hscale, pij, iSpin, iblk, li, lj)
    do iat = iAtFirst, iAtLast
      izp = species(iat)
      is = bas%ish_at(iat)
      io = bas%iao_sh(is+1)
      do iNeigh = 1, nNeighbour(iAt)
        img = iNeighbours(iNeigh, iAt)
        jAt = img2centCell(img)
        jZp = species(jAt)
        js = bas%ish_at(jAt)
        ind = iPair(iNeigh, iAt)
        jo = bas%iao_sh(js+1)
        nBlk = nOrbAtom(jAt)

        vec(:) = coords(:, iAt) - coords(:, img)
        r2 = vec(1)**2 + vec(2)**2 + vec(3)**2
        rr = sqrt(sqrt(r2) / (h0%rad(jzp) + h0%rad(izp)))
        do ish = 1, bas%nsh_id(izp)
          ii = bas%iao_sh(is+ish) - io
          do jsh = 1, bas%nsh_id(jzp)
            jj = bas%iao_sh(js+jsh) - jo
            call multipole_grad_cgto(bas%cgto(jsh, jzp), bas%cgto(ish, izp), &
                & r2, vec, bas%intcut, stmp, dtmp, qtmp, dstmp, ddtmpj, dqtmpj, &
                & ddtmpi, dqtmpi)
            !call overlap_grad_cgto(bas%cgto(jsh, jzp), bas%cgto(ish, izp), &
            !    & r2, vec, bas%intcut, stmp, dstmp)

            shpoly = (1.0_dp + h0%shpoly(ish, izp)*rr) &
              & * (1.0_dp + h0%shpoly(jsh, jzp)*rr)
            dshpoly = ((1.0_dp + h0%shpoly(ish, izp)*rr)*h0%shpoly(jsh, jzp)*rr &
              & + (1.0_dp + h0%shpoly(jsh, jzp)*rr)*h0%shpoly(ish, izp)*rr) &
              & * 0.5_dp / r2

            hscale = h0%hscale(jsh, ish, jzp, izp)
            hij = 0.5_dp * (selfenergy(is+ish) + selfenergy(js+jsh)) * hscale
            dhdcni = dsedcn(is+ish) * shpoly * hscale
            dhdcnj = dsedcn(js+jsh) * shpoly * hscale

            dG(:) = 0.0_dp
            dcni = 0.0_dp
            dcnj = 0.0_dp
            li = bas%cgto(iSh, iZp)%ang
            lj = bas%cgto(jSh, jZp)%ang
            nao = msao(lj)
            do iao = 1, msao(li)
              do jao = 1, nao
                ij = jao + nao*(iao-1)
                iblk = ind + jj+jao + nBlk*(ii+iao-1)

                pij = pmat(iblk, 1)
                hpij = pij * hij * shpoly
                sval = 2*hpij - 2*xmat(iblk)
                do iSpin = 1, nSpin
                   sval = sval + pmat(iblk, iSpin) * (shift(jj+jao, jj+jao, jat, iSpin) &
                      & + shift(ii+iao, ii+iao, iat, iSpin))
                end do

                dG(:) = dG + sval * dstmp(:, ij) &
                  & + 2*hpij*stmp(ij) * dshpoly / shpoly * vec &
                  & + pij * matmul(ddtmpi(:, :, ij), dipShift(:, iat)) &
                  & + pij * matmul(ddtmpj(:, :, ij), dipShift(:, jat)) &
                  & + pij * matmul(dqtmpi(:, :, ij), quadShift(:, iat)) &
                  & + pij * matmul(dqtmpj(:, :, ij), quadShift(:, jat))

                dcni = dcni + dhdcni * pij * stmp(ij)
                dcnj = dcnj + dhdcnj * pij * stmp(ij)
              end do
            end do
            dEdcn(iat) = dEdcn(iat) + dcni
            dEdcn(jat) = dEdcn(jat) + dcnj
            gradient(:, iat) = gradient(:, iat) + dG
            gradient(:, jat) = gradient(:, jat) - dG
            if (iat == jat) then
              sigma(:, :) = sigma + 0.25_dp * (spread(vec, 1, 3) * spread(dG, 2, 3) &
                & + spread(dG, 1, 3) * spread(vec, 2, 3))
            else
              sigma(:, :) = sigma + 0.50_dp * (spread(vec, 1, 3) * spread(dG, 2, 3) &
                & + spread(dG, 1, 3) * spread(vec, 2, 3))
            end if

          end do
        end do

      end do
    end do
  end subroutine buildDiatomicDerivs
#:endif


  !> Calculates nonadiabatic matrix: overlap gradient (Sprime) times velocities (Rdot)
  subroutine buildRdotSprime(this, env, RdotSprime, coords, dcoord, species, nNeighbour, &
      & iNeighbour, img2CentCell, iSquare, orb)

    !> Data structure
    class(TTBLite), intent(inout) :: this

    !> Computational environment settings
    type(TEnvironment), intent(in) :: env

    !> Nonadiabatic coupling matrix elements
    complex(dp), intent(out) :: RdotSprime(:,:)

    !> Coords of the atoms (3, nAllAtom)
    real(dp), intent(in) :: coords(:,:)

    !> Change in coords of the atoms (3, nAtom)
    real(dp), intent(in) :: dcoord(:,:)

    !> List of all atomic species
    integer, intent(in) :: species(:)

    !> Number of neighbours for atoms out to max interaction distance (excluding Ewald terms)
    integer, intent(in) :: nNeighbour(:)

    !> Neighbour list for atoms
    integer, intent(in) :: iNeighbour(0:,:)

    !> Image atoms to their equivalent in the central cell
    integer, intent(in) :: img2CentCell(:)

    !> Index array for start of atomic block in dense matrices
    integer, intent(in) :: iSquare(:)

    !> Data type for atomic orbital information
    type(TOrbitals), intent(in) :: orb

    call error("Forces currently not available in Ehrenfest dynamic with this Hamiltonian")
  end subroutine buildRdotSprime


#:if not WITH_TBLITE
  subroutine notImplementedError

    call error("DFTB+ compiled without support for tblite library")
  end subroutine notImplementedError
#:endif


end module dftbp_extlibs_tblite
