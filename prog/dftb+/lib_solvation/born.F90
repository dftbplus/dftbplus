!--------------------------------------------------------------------------------------------------!
!  DFTB+: general package for performing fast atomistic simulations                                !
!  Copyright (C) 2006 - 2020  DFTB+ developers group                                               !
!                                                                                                  !
!  See the LICENSE file for terms of usage and distribution.                                       !
!--------------------------------------------------------------------------------------------------!

#:include 'common.fypp'

!> Generalized Born solvation model.
module dftbp_born
  use dftbp_assert
  use dftbp_accuracy, only : dp
  use dftbp_blasroutines, only : hemv, gemv
  use dftbp_charges, only : getSummedCharges
  use dftbp_cm5, only : TChargeModel5, TCM5Input, TChargeModel5_init
  use dftbp_commontypes, only : TOrbitals
  use dftbp_constants, only : Hartree__eV
  use dftbp_environment, only : TEnvironment
  use dftbp_periodic, only : TNeighbourList, getNrOfNeighboursForAll
  use dftbp_sasa, only : TSASACont, TSASAInput, TSASACont_init, writeSASAContInfo
  use dftbp_schedule, only : distributeRangeInChunks, assembleChunks
  use dftbp_simplealgebra, only : determinant33
  use dftbp_solvation, only : TSolvation
  implicit none
  private

  public :: TGeneralizedBorn, TGBInput, TGeneralizedBorn_init
  public :: writeGeneralizedBornInfo, fgbKernel


  !> Possible interaction kernel
  type :: TFGBKernelEnum

    !> Canonical Still interaction kernel
    integer :: still = 1

    !> P16 interaction kernel
    integer :: p16 = 2

  end type TFGBKernelEnum

  !> Actual enumerator for available interaction kernel
  type(TFGBKernelEnum), parameter :: fgbKernel = TFGBKernelEnum()


  !> Global parameters for the solvation
  type :: TGBParameters

    !> Energy shift to the reference system
    real(dp) :: freeEnergyShift = 0.0_dp

    !> Dielectric screening
    real(dp) :: keps = 1.0_dp

    !> Scaling factor for Born radii
    real(dp) :: bornScale = 1.0_dp

    !> Offset parameter for Born radii
    real(dp) :: bornOffset = 0.0_dp

    !> Onufriev--Bashford--Case correction to Born radii
    real(dp) :: obc(3) = [1.00_dp, 0.80_dp, 4.85_dp]

    !> Van-der-Waals radii
    real(dp), allocatable :: vdwRad(:)

    !> Analytical linearized Poission-Boltzmann parameter alpha
    real(dp) :: alpbet = 0.0_dp

    !> Used interaction kernel
    integer :: kernel = fgbKernel%still

  end type TGBParameters


  !> Input parameters to initialize generalized Born model
  type, extends(TGBParameters) :: TGBInput

    !> Dielectric descreening parameter
    real(dp), allocatable :: descreening(:)

    !> Real space cutoff
    real(dp) :: rCutoff = 0.0_dp

    !> Use charge model 5
    type(TCM5Input), allocatable :: cm5Input

    !> Input for solvent accessible surface area model
    type(TSASAInput), allocatable :: sasaInput

    !> Parameter for H-bond correction
    real(dp), allocatable :: hBondPar(:)

  end type TGBInput


  !> Data for the Generalized Born solvation model
  type, extends(TSolvation) :: TGeneralizedBorn
    private

    !> number of atoms
    integer :: nAtom = 0

    !> solvation free energy
    real(dp), allocatable :: energies(:)

    !> lattice vectors if periodic
    real(dp) :: latVecs(3, 3) = 0.0_dp

    !> Volume of the unit cell
    real(dp) :: volume = 0.0_dp

    !> Strain derivatives
    real(dp) :: sigma(3, 3) = 0.0_dp

    !> is this periodic
    logical :: tPeriodic

    !> are the coordinates current?
    logical :: tCoordsUpdated = .false.

    !> are the charges current?
    logical :: tChargesUpdated = .false.

    !> Real space cutoff
    real(dp) :: rCutoff = 0.0_dp

    !> Model parameters
    type(TGBParameters) :: param

    !> Correction to charges with CM5
    type(TChargeModel5), allocatable :: cm5

    !> Born shifts to the hamiltonian
    real(dp), allocatable :: shift(:)

    !> Charges
    real(dp), allocatable :: chargesPerAtom(:)

    !> Born radii
    real(dp), allocatable :: bornRad(:)

    !> Born matrix
    real(dp), allocatable :: bornMat(:, :)

    !> Pair descreening approximation radii
    real(dp), allocatable :: rho(:)

    !> Gradient of the Born radii
    real(dp), allocatable :: dbrdr(:, :, :)

    !> Strain derivative of the Born radii
    real(dp), allocatable :: dbrdL(:, :, :)

    !> Solvent accessible surface area model
    type(TSASACont), allocatable :: sasaCont

    !> Parameter for H-bond correction
    real(dp), allocatable :: hBondStrength(:)

  contains

    !> update internal copy of coordinates
    procedure :: updateCoords

    !> update internal copy of lattice vectors
    procedure :: updateLatVecs

    !> get real space cutoff
    procedure :: getRCutoff

    !> get energy contributions
    procedure :: getEnergies

    !> get force contributions
    procedure :: addGradients

    !> get stress tensor contributions
    procedure :: getStress

    !> Updates with changed charges for the instance
    procedure :: updateCharges

    !> Returns shifts per atom
    procedure :: getShifts

    !> Query if object is actually an analytical linearized Poisson Boltzmann model
    procedure :: isALPB
  end type TGeneralizedBorn


  !> P16 zeta parameter
  real(dp), parameter :: zetaP16 = 1.028_dp

  !> P16 zeta parameter over 16
  real(dp), parameter :: zetaP16o16 = zetaP16 / 16.0_dp


contains


  !> Initialize generalized Born model from input data
  subroutine TGeneralizedBorn_init(this, input, nAtom, species0, speciesNames, &
      & latVecs)

    !> Initialised instance at return
    type(TGeneralizedBorn), intent(out) :: this

    !> Specific input parameters for generalized Born
    type(TGBinput), intent(in) :: input

    !> Nr. of atoms in the system
    integer, intent(in) :: nAtom

    !> Species of every atom in the unit cell
    integer, intent(in) :: species0(:)

    !> Symbols of the species
    character(len=*), intent(in) :: speciesNames(:)

    !> Lattice vectors, if the system is periodic
    real(dp), intent(in), optional :: latVecs(:,:)

    integer :: nSpecies
    integer :: iAt1, iSp1

    nSpecies = size(speciesNames)
    this%tPeriodic = present(latVecs)

    if (allocated(input%sasaInput)) then
       allocate(this%sasaCont)
       if (this%tPeriodic) then
         call TSASACont_init(this%sasaCont, input%sasaInput, nAtom, species0, &
             & speciesNames, latVecs)
       else
         call TSASACont_init(this%sasaCont, input%sasaInput, nAtom, species0, &
             & speciesNames)
       end if
    end if

    if (this%tPeriodic) then
      call this%updateLatVecs(LatVecs)
    end if
    this%nAtom = nAtom

    allocate(this%energies(nAtom))
    allocate(this%shift(nAtom))
    allocate(this%chargesPerAtom(nAtom))
    allocate(this%bornRad(nAtom))
    allocate(this%bornMat(nAtom, nAtom))
    allocate(this%rho(nSpecies))
    allocate(this%dbrdr(3, nAtom, nAtom))
    allocate(this%dbrdL(3, 3, nAtom))

    this%param = input%TGBParameters
    this%rho(:) = input%vdwRad(:) * input%descreening(:)

    if (allocated(this%sasaCont) .and. allocated(input%hBondPar)) then
      if (any(input%hBondPar /= 0.0_dp)) then
        allocate(this%hBondStrength(nAtom))
        do iAt1 = 1, nAtom
          iSp1 = species0(iAt1)
          this%hBondStrength(iAt1) = input%hBondPar(iSp1) / this%sasaCont%probeRad(iSp1)**2
        end do
      end if
    end if

    this%rCutoff = input%rCutoff

    if (allocated(input%cm5Input)) then
      allocate(this%cm5)
      if (this%tPeriodic) then
        call TChargeModel5_init(this%cm5, input%cm5Input, nAtom, speciesNames, &
           & .true., latVecs)
      else
        call TChargeModel5_init(this%cm5, input%cm5Input, nAtom, speciesNames, &
           & .true.)
      end if
    end if

    this%tCoordsUpdated = .false.
    this%tChargesUpdated = .false.

  end subroutine TGeneralizedBorn_init


  !> Print the solvation model used
  subroutine writeGeneralizedBornInfo(unit, solvation)

    !> Formatted unit for IO
    integer, intent(in) :: unit

    !> Solvation model
    type(TGeneralizedBorn), intent(in) :: solvation

    write(unit, '(a, ":", t30, es14.6)') "Dielectric constant", &
        & 1.0_dp/(solvation%param%keps * (1.0_dp + solvation%param%alpbet) + 1.0_dp)
    write(unit, '(a, ":", t30, es14.6, 1x, a, t50, es14.6, 1x, a)') &
        & "Free energy shift", solvation%param%freeEnergyShift, "H", &
        & Hartree__eV * solvation%param%freeEnergyShift, "eV"

    write(unit, '(a, ":", t30)', advance='no') "Born interaction kernel"
    select case(solvation%param%kernel)
    case default
      write(unit, '(a)') "unknown (internal error)"
    case(fgbKernel%still)
      write(unit, '(a)') "Still"
    case(fgbKernel%p16)
      write(unit, '(a)') "P16"
    end select

    write(unit, '(a, ":", t30, a)') "Born radii integrator", "GBOBC"

    write(unit, '(a, ":", t30)', advance='no') "SASA model"
    if (allocated(solvation%sasaCont)) then
      write(unit, '(a)') "Yes"
      call writeSASAContInfo(unit, solvation%sasaCont)
    else
      write(unit, '(a)') "No"
    end if

    write(unit, '(a, ":", t30)', advance='no') "CM5 correction"
    if (allocated(solvation%cm5)) then
      write(unit, '(a)') "Yes"
    else
      write(unit, '(a)') "No"
    end if

    write(unit, '(a, ":", t30)', advance='no') "Hydrogen bond correction"
    if (allocated(solvation%hBondStrength)) then
      write(unit, '(a)') "Yes"
    else
      write(unit, '(a)') "No"
    end if
  end subroutine writeGeneralizedBornInfo


  !> Check if this is actually an analyical linearized Poisson-Boltzmann model
  !  masquerading as a generalized Born one
  pure function isALPB(this) result(alpb)

    !> Data structure
    class(TGeneralizedBorn), intent(in) :: this

    !> Analytical linearized Poisson-Boltzmann model used
    logical :: alpb

    alpb = this%param%alpbet > 0.0_dp

  end function isALPB


  !> Update internal stored coordinates
  subroutine updateCoords(this, env, neighList, img2CentCell, coords, species0)

    !> Data structure
    class(TGeneralizedBorn), intent(inout) :: this

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

    integer, allocatable :: nNeigh(:)
    real(dp) :: aDet

    if (allocated(this%sasaCont)) then
      call this%sasaCont%updateCoords(env, neighList, img2CentCell, coords, species0)
    end if

    allocate(nNeigh(this%nAtom))
    call getNrOfNeighboursForAll(nNeigh, neighList, this%rCutoff)
    call getBornRadii(env, nNeigh, neighList%iNeighbour, img2CentCell, &
        & neighList%neighDist2, species0, coords, this%param%vdwRad, &
        & this%rho, this%param%bornOffset, this%param%bornScale, this%param%obc, &
        & this%bornRad, this%dbrdr)
    call getBornMatrixCluster(env, this%nAtom, coords, this%param%kernel, &
        & this%param%keps, this%bornRad, this%bornMat)

    ! Analytical linearized Poission-Boltzmann contribution for charged systems
    if (this%param%alpbet > 0.0_dp) then
      call getADet(this%nAtom, coords, species0, this%param%vdwRad, aDet)
      this%bornMat(:, :) = this%bornMat + this%param%kEps * this%param%alpbet / aDet
    end if

    if (allocated(this%cm5)) then
      call this%cm5%updateCoords(neighList, img2CentCell, coords, species0)
    end if

    this%tCoordsUpdated = .true.
    this%tChargesUpdated = .false.

  end subroutine updateCoords


  !> Update internal copy of lattice vectors
  subroutine updateLatVecs(this, latVecs)

    !> Data structure
    class(TGeneralizedBorn), intent(inout) :: this

    !> Lattice vectors
    real(dp), intent(in) :: latVecs(:,:)

    @:ASSERT(this%tPeriodic)
    @:ASSERT(all(shape(latvecs) == shape(this%latvecs)))

    if (allocated(this%sasaCont)) then
      call this%sasaCont%updateLatVecs(latVecs)
    end if

    this%volume = abs(determinant33(latVecs))
    this%latVecs(:,:) = latVecs

    if (allocated(this%cm5)) then
      call this%cm5%updateLatVecs(LatVecs)
    end if

    this%tCoordsUpdated = .false.
    this%tChargesUpdated = .false.

  end subroutine updateLatVecs


  !> Get energy contributions
  subroutine getEnergies(this, energies)

    !> data structure
    class(TGeneralizedBorn), intent(inout) :: this

    !> energy contributions for each atom
    real(dp), intent(out) :: energies(:)

    @:ASSERT(this%tCoordsUpdated)
    @:ASSERT(this%tChargesUpdated)
    @:ASSERT(size(energies) == this%nAtom)

    if (allocated(this%sasaCont)) then
      call this%sasaCont%getEnergies(energies)
    else
      energies(:) = 0.0_dp
    end if

    energies(:) = energies + 0.5_dp * (this%shift * this%chargesPerAtom) &
       & + this%param%freeEnergyShift / real(this%nAtom, dp)

  end subroutine getEnergies


  !> Get force contributions
  subroutine addGradients(this, env, neighList, species, coords, img2CentCell, gradients)

    !> Data structure
    class(TGeneralizedBorn), intent(inout) :: this

    !> Computational environment settings
    type(TEnvironment), intent(in) :: env

    !> Neighbour list
    type(TNeighbourList), intent(in) :: neighList

    !> Specie for each atom
    integer, intent(in) :: species(:)

    !> Coordinate of each atom
    real(dp), intent(in) :: coords(:,:)

    !> Mapping of atoms to cetnral cell
    integer, intent(in) :: img2CentCell(:)

    !> Gradient contributions for each atom
    real(dp), intent(inout) :: gradients(:,:)

    real(dp) :: sigma(3, 3)
    real(dp), allocatable :: dEdcm5(:)
    integer, allocatable :: nNeigh(:)
    real(dp), allocatable :: dhbds(:)

    @:ASSERT(this%tCoordsUpdated)
    @:ASSERT(this%tChargesUpdated)
    @:ASSERT(all(shape(gradients) == [3, this%nAtom]))

    if (allocated(this%sasaCont)) then
      call this%sasaCont%addGradients(env, neighList, species, coords, img2CentCell, gradients)
      if (allocated(this%hBondStrength)) then
        allocate(dhbds(this%nAtom))
        dhbds(:) = this%hBondStrength * this%chargesPerAtom**2
        call gemv(gradients, this%sasaCont%dsdr, dhbds, beta=1.0_dp)
        deallocate(dhbds)
      end if
    end if

    allocate(nNeigh(this%nAtom))
    sigma(:, :) = 0.0_dp
    this%energies(:) = 0.0_dp

    call getNrOfNeighboursForAll(nNeigh, neighList, this%rCutoff)
    call getBornEGCluster(env, this%nAtom, coords, this%chargesPerAtom, &
      & this%bornRad, this%dbrdr, this%dbrdL, this%param%kernel, this%param%keps, &
      & this%energies, gradients, sigma)

    ! Analytical linearized Poission-Boltzmann contribution for charged systems
    if (this%param%alpbet > 0.0_dp) then
      call getADetDeriv(this%nAtom, coords, species, this%param%vdwRad, &
          & this%param%kEps*this%param%alpbet, this%chargesPerAtom, gradients)
    end if

    if (allocated(this%cm5)) then
      allocate(dEdcm5(this%nAtom))
      dEdcm5(:) = 0.0_dp
      call hemv(dEdcm5, this%bornMat, this%chargesPerAtom)
      call this%cm5%addGradients(dEdcm5, gradients)
      call this%cm5%addSigma(dEdcm5, sigma)
    end if

    this%energies = this%energies + this%param%freeEnergyShift / real(this%nAtom, dp)

    if (this%tPeriodic) then
      this%sigma(:, :) = sigma
    end if

  end subroutine addGradients


  !> get stress tensor contributions
  subroutine getStress(this, stress)

    !> data structure
    class(TGeneralizedBorn), intent(inout) :: this

    !> Stress tensor contributions
    real(dp), intent(out) :: stress(:,:)

    @:ASSERT(this%tCoordsUpdated)
    @:ASSERT(this%tChargesUpdated)
    @:ASSERT(all(shape(stress) == [3, 3]))
    @:ASSERT(this%tPeriodic)
    @:ASSERT(this%volume > 0.0_dp)

    if (allocated(this%sasaCont)) then
      call this%sasaCont%getStress(stress)
    else
      stress(:, :) = 0.0_dp
    end if

    stress(:,:) = stress + this%sigma / this%volume

  end subroutine getStress


  !> Distance cut off for generalized Born calculations
  function getRCutoff(this) result(cutoff)

    !> data structure
    class(TGeneralizedBorn), intent(inout) :: this

    !> resulting cutoff
    real(dp) :: cutoff

    cutoff = this%rCutoff
    if (allocated(this%cm5)) then
      cutoff = max(cutoff, this%cm5%getRCutoff())
    end if

    if (allocated(this%sasaCont)) then
      cutoff = max(cutoff, this%sasaCont%getRCutoff())
    end if

  end function getRCutoff


  !> Updates with changed charges for the instance.
  subroutine updateCharges(this, env, species, neighList, qq, q0, img2CentCell, orb)

    !> Data structure
    class(TGeneralizedBorn), intent(inout) :: this

    !> Computational environment settings
    type(TEnvironment), intent(in) :: env

    !> Species, shape: [nAtom]
    integer, intent(in) :: species(:)

    !> Neighbour list.
    type(TNeighbourList), intent(in) :: neighList

    !> Orbital charges.
    real(dp), intent(in) :: qq(:,:,:)

    !> Reference orbital charges.
    real(dp), intent(in) :: q0(:,:,:)

    !> Mapping on atoms in central cell.
    integer, intent(in) :: img2CentCell(:)

    !> Orbital information
    type(TOrbitals), intent(in) :: orb

    @:ASSERT(this%tCoordsUpdated)

    if (allocated(this%sasaCont)) then
      call this%sasaCont%updateCharges(env, species, neighList, qq, q0, img2CentCell, orb)
    end if

    call getSummedCharges(species, orb, qq, q0, dQAtom=this%chargesPerAtom)
    if (allocated(this%cm5)) then
      call this%cm5%addCharges(this%chargesPerAtom)
    end if

    if (allocated(this%sasaCont) .and. allocated(this%hBondStrength)) then
      this%shift(:) = 2.0_dp * this%sasaCont%sasa * this%hBondStrength * this%chargesPerAtom
    else
      this%shift(:) = 0.0_dp
    end if
    call hemv(this%shift, this%bornMat, this%chargesPerAtom, beta=1.0_dp)

    this%tChargesUpdated = .true.

  end subroutine updateCharges


  !> Returns shifts per atom
  subroutine getShifts(this, shiftPerAtom, shiftPerShell)

    !> Data structure
    class(TGeneralizedBorn), intent(inout) :: this

    !> Shift per atom
    real(dp), intent(out) :: shiftPerAtom(:)

    !> Shift per shell
    real(dp), intent(out) :: shiftPerShell(:,:)

    @:ASSERT(this%tCoordsUpdated)
    @:ASSERT(this%tChargesUpdated)
    @:ASSERT(size(shiftPerAtom) == this%nAtom)
    @:ASSERT(size(shiftPerShell, dim=2) == this%nAtom)

    if (allocated(this%sasaCont)) then
      call this%sasaCont%getShifts(shiftPerAtom, shiftPerShell)
    else
      shiftPerAtom(:) = 0.0_dp
      shiftPerShell(:,:) = 0.0_dp
    end if

    shiftPerAtom(:) = shiftPerAtom + this%shift

  end subroutine getShifts


  !> Calculate Born radii for a given geometry
  subroutine getBornRadii(env, nNeighbour, iNeighbour, img2CentCell, &
      & neighDist2, species, coords, vdwRad, rho, bornOffset, bornScale, obcPar, &
      & bornRad, dbrdr)

    !> Computational environment settings
    type(TEnvironment), intent(in) :: env

    !> Nr. of neighbours for each atom
    integer, intent(in) :: nNeighbour(:)

    !> Neighbourlist
    integer, intent(in) :: iNeighbour(0:, :)

    !> Mapping into the central cell
    integer, intent(in) :: img2CentCell(:)

    !> Square distances of the neighbours
    real(dp), intent(in) :: neighDist2(0:, :)

    !> Central cell chemical species
    integer, intent(in) :: species(:)

    !> current atomic positions
    real(dp), intent(in) :: coords(:, :)

    !> Van-der-Waals radii
    real(dp), intent(in) :: vdwRad(:)

    !> Descreened radii
    real(dp), intent(in) :: rho(:)

    !> Offset parameter for the Born radii
    real(dp), intent(in) :: bornOffset

    !> Scaling parameter for the Born radii
    real(dp), intent(in) :: bornScale

    !> Onufriev-Bashford-Case correction parameters
    real(dp), intent(in) :: obcPar(3)

    !> Born radii
    real(dp), intent(out) :: bornRad(:)

    !> Derivatives of the Born radii w.r.t. the cartesian coordinates
    real(dp), intent(out) :: dbrdr(:,:,:)

    integer :: nAtom, iAt1, iSp1
    real(dp) :: br, dpsi, svdwi,vdwri, s1, v1, s2, arg, arg2, th, ch

    nAtom = size(nNeighbour)

    call getPsi(env, nNeighbour, iNeighbour, img2CentCell, neighDist2, &
        & species, coords, vdwRad, rho, bornRad, dbrdr)

    !$omp parallel do default(none) schedule(runtime) shared(bornRad, dbrdr) &
    !$omp shared(nAtom, species, vdwRad, bornOffset, bornScale, obcPar) &
    !$omp private(iAt1, iSp1, br, dpsi, svdwi,vdwri, s1, v1, s2, arg, arg2, th, ch)
    do iAt1 = 1, nAtom
      iSp1 = species(iAt1)

      br = bornRad(iAt1)

      svdwi = vdwRad(iSp1) - bornOffset
      vdwri = vdwRad(iSp1)
      s1 = 1.0_dp/svdwi
      v1 = 1.0_dp/vdwri
      s2 = 0.5_dp*svdwi

      br = br*s2

      arg2 = br*(obcPar(3)*br-obcPar(2))
      arg = br*(obcPar(1)+arg2)
      arg2 = 2.0_dp*arg2+obcPar(1)+obcPar(3)*br*br

      th = tanh(arg)
      ch = cosh(arg)

      br = 1.0_dp/(s1-v1*th)
      ! Include GBMV2-like scaling
      br = bornScale*br

      dpsi = ch*(s1-v1*th)
      dpsi = s2*v1*arg2/(dpsi*dpsi)
      dpsi = bornScale*dpsi

      bornRad(iAt1) = br
      dbrdr(:, :, iAt1) = dbrdr(:, :, iAt1) * dpsi
      !dbrdL(:, :, iAt1) = dbrdL(:, :, iAt1) * dpsi

    end do
    !$omp end parallel do

  end subroutine getBornRadii


  !> Evaluate volume integrals, intermediate values are stored in Born radii fields
  subroutine getPsi(env, nNeighbour, iNeighbour, img2CentCell, &
      & neighDist2, species, coords, vdwRad, rho, psi, dpsidr)

    !> Computational environment settings
    type(TEnvironment), intent(in) :: env

    !> Nr. of neighbours for each atom
    integer, intent(in) :: nNeighbour(:)

    !> Neighbourlist
    integer, intent(in) :: iNeighbour(0:, :)

    !> Mapping into the central cell
    integer, intent(in) :: img2CentCell(:)

    !> Square distances of the neighbours
    real(dp), intent(in) :: neighDist2(0:, :)

    !> Central cell chemical species
    integer, intent(in) :: species(:)

    !> current atomic positions
    real(dp), intent(in) :: coords(:, :)

    !> Van-der-Waals radii
    real(dp), intent(in) :: vdwRad(:)

    !> Descreened radii
    real(dp), intent(in) :: rho(:)

    !> Atomic volumes
    real(dp), intent(out) :: psi(:)

    !> Derivative of atomic volumes
    real(dp), intent(out) :: dpsidr(:,:,:)

    integer :: nAtom, iAtFirst, iAtLast, iAt1, iNeigh, iAt2, iAt2f, iSp1, iSp2
    logical :: tOvij, tOvji
    real(dp) :: vec(3), dist, rhoi, rhoj
    real(dp) :: gi, gj, ap, am, lnab, rhab, ab, dgi, dgj
    real(dp) :: dGr(3)
    real(dp) :: rh1, rhr1, r24, r1, aprh1, r12
    real(dp) :: rvdwi, rvdwj
    real(dp), allocatable :: dpsitr(:,:)

    nAtom = size(nNeighbour)

    call distributeRangeInChunks(env, 1, nAtom, iAtFirst, iAtLast)

    allocate(dpsitr(3, nAtom))
    psi(:) = 0.0_dp
    dpsidr(:, :, :) = 0.0_dp
    dpsitr(:, :) = 0.0_dp

    !$omp parallel do default(none) schedule(runtime) &
    !$omp reduction(+:psi, dpsidr, dpsitr) shared(iAtFirst, iAtLast, species) &
    !$omp shared(nNeighbour, iNeighbour, img2CentCell, coords, neighDist2, rho) &
    !$omp shared(vdwRad) private(iAt1, iSp1, iNeigh, iAt2, iAt2f, iSp2, dist) &
    !$omp private(tOvij, tOvji, vec, rhoi, rhoj, gi, gj, ap, am, lnab, rhab) &
    !$omp private(ab, dgi, dgj, dGr, rh1, rhr1, r24, r1, aprh1, r12, rvdwi, rvdwj)
    do iAt1 = iAtFirst, iAtLast
      iSp1 = species(iAt1)
      do iNeigh = 1, nNeighbour(iAt1)
        iAt2 = iNeighbour(iNeigh, iAt1)
        iAt2f = img2CentCell(iAt2)
        iSp2 = species(iAt2f)
        vec(:) = coords(:, iAt1) - coords(:, iAt2)
        dist = sqrt(neighDist2(iNeigh, iAt1))

        rhoi = rho(iSp1)
        rhoj = rho(iSp2)
        rvdwi = vdwRad(iSp1)
        rvdwj = vdwRad(iSp2)

        tOvij = dist < (rvdwi + rhoj)
        tOvji = dist < (rhoi + rvdwj)

        tOverlap: if (.not. tOvij .and. .not. tOvji) then ! ij do not overlap; ji do not overlap
          ! nonoverlaping spheres
          if(abs(rhoi-rhoj) < 1.e-8_dp) then
            ! equal reduced radiAt1
            r1 = 1.0_dp/dist
            ap = dist+rhoj
            am = dist-rhoj
            ab = ap*am
            rhab = rhoj/ab
            lnab = 0.5_dp*log(am/ap)*r1
            gi = rhab+lnab
            dgi = -2.0_dp*rhab/ab+(rhab-lnab)*r1*r1
            ! accumulate psi
            psi(iAt1) = psi(iAt1)+gi
            psi(iAt2f) = psi(iAt2f)+gi
            ! accumulate psi gradient
            dGr(:) = dgi*vec(:)
            dpsitr(:,iAt1) = dpsitr(:,iAt1)+dGr(:)
            dpsidr(:,iAt2f,iAt1) = dpsidr(:,iAt2f,iAt1)-dGr(:)
            dpsitr(:,iAt2f) = dpsitr(:,iAt2f)-dGr(:)
            dpsidr(:,iAt1,iAt2f) = dpsidr(:,iAt1,iAt2f)+dGr(:)
          else
            ! unequal reduced radiAt1
            ! ij contribution
            r1 = 1.0_dp/dist
            ap = dist+rhoj
            am = dist-rhoj
            ab = ap*am
            rhab = rhoj/ab
            lnab = 0.5_dp*log(am/ap)*r1
            gi = rhab+lnab
            dgi = -2.0_dp*rhab/ab+(rhab-lnab)*r1*r1
            ! ji contribution
            ap = dist+rhoi
            am = dist-rhoi
            ab = ap*am
            rhab = rhoi/ab
            lnab = 0.5_dp*log(am/ap)*r1
            gj = rhab+lnab
            dgj = -2.0_dp*rhab/ab+(rhab-lnab)*r1*r1
            ! accumulate psi
            psi(iAt1) = psi(iAt1)+gi
            psi(iAt2f) = psi(iAt2f)+gj
            ! accumulate psi gradient
            dGr(:) = dgi*vec(:)
            dpsitr(:,iAt1) = dpsitr(:,iAt1)+dGr(:)
            dpsidr(:,iAt2f,iAt1) = dpsidr(:,iAt2f,iAt1)-dGr(:)

            dGr(:) = dgj*vec(:)
            dpsitr(:,iAt2f) = dpsitr(:,iAt2f)-dGr(:)
            dpsidr(:,iAt1,iAt2f) = dpsidr(:,iAt1,iAt2f)+dGr(:)
          end if

        else if (.not. tOvij .and. tOvji) then tOverlap ! ij do not overlap; ji overlap

          ! ij contribution
          r1 = 1.0_dp/dist
          ap = dist+rhoj
          am = dist-rhoj
          ab = ap*am
          rhab = rhoj/ab
          lnab = 0.5_dp*log(am/ap)*r1
          gi = rhab+lnab
          dgi = -2.0_dp*rhab/ab+(rhab-lnab)*r1*r1
          ! accumulate psi
          psi(iAt1) = psi(iAt1)+gi
          ! accumulate psi gradient
          dGr(:) = dgi*vec(:)
          dpsitr(:,iAt1) = dpsitr(:,iAt1)+dGr(:)
          dpsidr(:,iAt2f,iAt1) = dpsidr(:,iAt2f,iAt1)-dGr(:)

          if((dist+rhoi) > rvdwj) then
            ! ji contribution
            r1 = 1.0_dp/dist
            r12 = 0.5_dp*r1
            r24 = r12*r12

            ap = dist+rhoi
            am = dist-rhoi
            rh1 = 1.0_dp/rvdwj
            rhr1 = 1.0_dp/ap
            aprh1 = ap*rh1
            lnab = log(aprh1)

            gj = rh1-rhr1+r12*(0.5_dp*am*(rhr1-rh1*aprh1)-lnab)

            dgj = rhr1*rhr1*(1.0_dp-0.25_dp*am*r1*(1.0_dp+aprh1*aprh1))+ &
              & rhoi*r24*(rhr1-rh1*aprh1)+r12*(r1*lnab-rhr1)
            dgj = dgj*r1
            ! accumulate psi
            psi(iAt2f) = psi(iAt2f)+gj
            ! accumulate psi gradient
            dGr(:) = dgj*vec(:)
            dpsitr(:,iAt2f) = dpsitr(:,iAt2f)-dGr(:)
            dpsidr(:,iAt1,iAt2f) = dpsidr(:,iAt1,iAt2f)+dGr(:)
          end if

        else if (tOvij .and. .not. tOvji) then ! ij overlap; ji do not overlap

          if((dist+rhoj) > rvdwi) then
            ! ij contribution
            r1 = 1.0_dp/dist
            r12 = 0.5_dp*r1
            r24 = r12*r12

            ap = dist+rhoj
            am = dist-rhoj
            rh1 = 1.0_dp/rvdwi
            rhr1 = 1.0_dp/ap
            aprh1 = ap*rh1
            lnab = log(aprh1)

            gi = rh1-rhr1+r12*(0.5_dp*am*(rhr1-rh1*aprh1)-lnab)

            dgi = rhr1*rhr1*(1.0_dp-0.25_dp*am*r1*(1.0_dp+aprh1*aprh1))+ &
              & rhoj*r24*(rhr1-rh1*aprh1)+r12*(r1*lnab-rhr1)
            dgi = dgi*r1
            ! accumulate psi
            psi(iAt1) = psi(iAt1)+gi
            ! accumulate psi gradient
            dGr(:) = dgi*vec(:)
            dpsitr(:,iAt1) = dpsitr(:,iAt1)+dGr(:)
            dpsidr(:,iAt2f,iAt1) = dpsidr(:,iAt2f,iAt1)-dGr(:)
          end if

          ! ji contribution
          ap = dist+rhoi
          am = dist-rhoi
          ab = ap*am
          rhab = rhoi/ab
          lnab = 0.5_dp*log(am/ap)*r1
          gj = rhab+lnab
          dgj = -2.0_dp*rhab/ab+(rhab-lnab)*r1*r1
          ! accumulate psi
          psi(iAt2f) = psi(iAt2f)+gj
          ! accumulate psi gradient
          dGr(:) = dgj*vec(:)
          dpsitr(:,iAt2f) = dpsitr(:,iAt2f)-dGr(:)
          dpsidr(:,iAt1,iAt2f) = dpsidr(:,iAt1,iAt2f)+dGr(:)

        else if (tOvij .and. tOvji) then tOverlap ! ij and ji overlap
          ! overlaping spheres
          if((dist+rhoj) > rvdwi) then
            ! ij contribution
            r1 = 1.0_dp/dist
            r12 = 0.5_dp*r1
            r24 = r12*r12

            ap = dist+rhoj
            am = dist-rhoj
            rh1 = 1.0_dp/rvdwi
            rhr1 = 1.0_dp/ap
            aprh1 = ap*rh1
            lnab = log(aprh1)

            gi = rh1-rhr1+r12*(0.5_dp*am*(rhr1-rh1*aprh1)-lnab)

            dgi = rhr1*rhr1*(1.0_dp-0.25_dp*am*r1*(1.0_dp+aprh1*aprh1))+ &
              & rhoj*r24*(rhr1-rh1*aprh1)+r12*(r1*lnab-rhr1)
            dgi = dgi*r1
            ! accumulate psi
            psi(iAt1) = psi(iAt1)+gi
            ! accumulate psi gradient
            dGr(:) = dgi*vec(:)
            dpsitr(:,iAt1) = dpsitr(:,iAt1)+dGr(:)
            dpsidr(:,iAt2f,iAt1) = dpsidr(:,iAt2f,iAt1)-dGr(:)
          end if

          if((dist+rhoi) > rvdwj) then
            ! ji contribution
            r1 = 1.0_dp/dist
            r12 = 0.5_dp*r1
            r24 = r12*r12

            ap = dist+rhoi
            am = dist-rhoi
            rh1 = 1.0_dp/rvdwj
            rhr1 = 1.0_dp/ap
            aprh1 = ap*rh1
            lnab = log(aprh1)

            gj = rh1-rhr1+r12*(0.5_dp*am*(rhr1-rh1*aprh1)-lnab)

            dgj = rhr1*rhr1*(1.0_dp-0.25_dp*am*r1*(1.0_dp+aprh1*aprh1))+ &
              & rhoi*r24*(rhr1-rh1*aprh1)+r12*(r1*lnab-rhr1)
            dgj = dgj*r1
            ! accumulate psi
            psi(iAt2f) = psi(iAt2f)+gj
            ! accumulate psi gradient
            dGr(:) = dgj*vec(:)
            dpsitr(:,iAt2f) = dpsitr(:,iAt2f)-dGr(:)
            dpsidr(:,iAt1,iAt2f) = dpsidr(:,iAt1,iAt2f)+dGr(:)
          end if

        end if tOverlap

      end do
    end do
    !$omp end parallel do

    ! save one-center terms
    do iAt1 = 1, nAtom
      dpsidr(:,iAt1,iAt1) = dpsidr(:,iAt1,iAt1) + dpsitr(:,iAt1)
    end do

    call assembleChunks(env, psi)
    call assembleChunks(env, dpsidr)

  end subroutine getPsi


  !> compute Born matrix
  subroutine getBornMatrixCluster(env, nAtom, coords0, kernel, keps, bornRad, &
      & bornMat)

    !> Computational environment settings
    type(TEnvironment), intent(in) :: env

    !> Number of atoms
    integer, intent(in) :: nAtom

    !> coordinates in the central cell
    real(dp), intent(in) :: coords0(:, :)

    !> Born radii
    real(dp), intent(in) :: bornRad(:)

    !> Dielectric screening
    real(dp), intent(in) :: keps

    !> Interaction kernel
    integer, intent(in) :: kernel

    !> Born matrix
    real(dp), intent(out) :: bornMat(:, :)

    integer :: iAtFirst, iAtLast, iAt1

    call distributeRangeInChunks(env, 1, nAtom, iAtFirst, iAtLast)

    bornMat(:, :) = 0.0_dp

    select case(kernel)
    case(fgbKernel%still)
      call getBornMatrixStillCluster(iAtFirst, iAtLast, bornRad, coords0, &
          & keps, bornMat)
    case(fgbKernel%p16)
      call getBornMatrixP16Cluster(iAtFirst, iAtLast, bornRad, coords0, &
          & keps, bornMat)
    end select

    !> self-energy part
    do iAt1 = iAtFirst, iAtLast
      bornMat(iAt1, iAt1) = keps/bornRad(iAt1)
    end do

    call assembleChunks(env, bornMat)

  end subroutine getBornMatrixCluster


  !> compute Born matrix using Still interaction kernel
  subroutine getBornMatrixStillCluster(iAtFirst, iAtLast, bornRad, coords0, keps, &
      & bornMat)

    !> Number of atoms
    integer, intent(in) :: iAtFirst

    !> Number of atoms
    integer, intent(in) :: iAtLast

    !> Born radii for each atom
    real(dp), intent(in) :: bornRad(:)

    !> coordinates in the central cell
    real(dp), intent(in) :: coords0(:, :)

    !> Dielectric scaling
    real(dp), intent(in) :: keps

    !> Born matrix
    real(dp), intent(inout) :: bornMat(:, :)

    integer :: iAt1, iAt2
    real(dp) :: aa, dist2, dd, expd, dfgb, fgb

    !$omp parallel do default(none) &
    !$omp shared(bornMat, iAtFirst, iAtLast, coords0, bornRad, kEps) &
    !$omp private(iAt1, iAt2, dist2, aa, dd, expd, fgb, dfgb)
    do iAt1 = iAtFirst, iAtLast
      do iAt2 = 1, iAt1-1
        dist2 = sum((coords0(:, iAt1) - coords0(:, iAt2))**2)

        aa = bornRad(iAt1)*bornRad(iAt2)
        dd = 0.25_dp*dist2/aa
        expd = exp(-dd)
        dfgb = 1.0_dp/(dist2+aa*expd)
        fgb = keps*sqrt(dfgb)
        bornMat(iAt1, iAt2) = bornMat(iAt1, iAt2) + fgb
        bornMat(iAt2, iAt1) = bornMat(iAt2, iAt1) + fgb
      end do
    end do
    !$omp end parallel do

  end subroutine getBornMatrixStillCluster


  !> compute Born matrix using Still interaction kernel
  subroutine getBornMatrixP16Cluster(iAtFirst, iAtLast, bornRad, coords0, keps, &
      & bornMat)

    !> Number of atoms
    integer, intent(in) :: iAtFirst

    !> Number of atoms
    integer, intent(in) :: iAtLast

    !> Born radii for each atom
    real(dp), intent(in) :: bornRad(:)

    !> coordinates in the central cell
    real(dp), intent(in) :: coords0(:, :)

    !> Dielectric scaling
    real(dp), intent(in) :: keps

    !> Born matrix
    real(dp), intent(inout) :: bornMat(:, :)

    integer :: iAt1, iAt2
    real(dp) :: r1, ab, arg, fgb, dfgb

    !$omp parallel do default(none) &
    !$omp shared(bornMat, iAtFirst, iAtLast, coords0, bornRad, kEps) &
    !$omp private(iAt1, iAt2, r1, ab, arg, fgb, dfgb)
    do iAt1 = iAtFirst, iAtLast
      do iAt2 = 1, iAt1-1
        r1 = sqrt(sum((coords0(:, iAt1) - coords0(:, iAt2))**2))
        ab = sqrt(bornRad(iAt1) * bornRad(iAt2))
        arg = ab / (ab + zetaP16o16*r1) ! ab / (1 + ζR/(16·ab))
        arg = arg * arg ! ab / (1 + ζR/(16·ab))²
        arg = arg * arg ! ab / (1 + ζR/(16·ab))⁴
        arg = arg * arg ! ab / (1 + ζR/(16·ab))⁸
        arg = arg * arg ! ab / (1 + ζR/(16·ab))¹⁶
        fgb = r1 + ab*arg
        dfgb = 1.0_dp / fgb
        bornMat(iAt2, iAt1) = bornMat(iAt2, iAt1) + dfgb * kEps
        bornMat(iAt1, iAt2) = bornMat(iAt1, iAt2) + dfgb * kEps
      end do
    end do
    !$omp end parallel do

  end subroutine getBornMatrixP16Cluster


  !> GB energy and gradient
  subroutine getBornEGCluster(env, nAtom, coords, chargesPerAtom, bornRad, &
      & dbrdr, dbrdL, kernel, keps, energies, gradients, sigma)

    !> Computational environment settings
    type(TEnvironment), intent(in) :: env

    !> Number of atoms
    integer, intent(in) :: nAtom

    !> Current atomic positions
    real(dp), intent(in) :: coords(:, :)

    !> Charges
    real(dp), intent(in) :: chargesPerAtom(:)

    !> Born radii
    real(dp), intent(in) :: bornRad(:)

    !> Gradient of the Born radii
    real(dp), intent(in) :: dbrdr(:, :, :)

    !> Strain derivative of the Born radii
    real(dp), intent(in) :: dbrdL(:, :, :)

    !> Dielectric screening
    real(dp), intent(in) :: keps

    !> Interaction kernel
    integer, intent(in) :: kernel

    !> Atom resolved energies
    real(dp), intent(out) :: energies(:)

    !> Molecular gradient
    real(dp), intent(inout) :: gradients(:, :)

    !> Strain derivative
    real(dp), intent(inout) :: sigma(:, :)

    integer :: iAt1, iAtFirst, iAtLast
    real(dp) :: qq, bp
    real(dp) :: grddbi
    real(dp), allocatable :: dEdbr(:), localSigma(:, :), derivs(:, :)

    call distributeRangeInChunks(env, 1, nAtom, iAtFirst, iAtLast)

    allocate(dEdbr(nAtom), derivs(3, nAtom), localSigma(3, 3))

    localSigma(:, :) = 0.0_dp
    derivs(:, :) = 0.0_dp
    energies(:) = 0.0_dp
    dEdbr(:) = 0.0_dp

    select case(kernel)
    case(fgbKernel%still)
      call getBornEGStillCluster(iAtFirst, iAtLast, coords, chargesPerAtom, &
          & bornRad, keps, energies, derivs, localSigma, dEdbr)
    case(fgbKernel%p16)
      call getBornEGP16Cluster(iAtFirst, iAtLast, coords, chargesPerAtom, &
          & bornRad, keps, energies, derivs, localSigma, dEdbr)
    end select

    !> self-energy part
    do iAt1 = iAtFirst, iAtLast
      bp = 1.0_dp/bornRad(iAt1)
      qq = chargesPerAtom(iAt1)*bp
      energies(iAt1) = energies(iAt1) + 0.5_dp*chargesPerAtom(iAt1)*qq*keps
      grddbi = -0.5_dp*keps*qq*bp
      dEdbr(iAt1) = dEdbr(iAt1) + grddbi*chargesPerAtom(iAt1)
    end do

    call assembleChunks(env, energies)
    call assembleChunks(env, derivs)
    call assembleChunks(env, localSigma)
    call assembleChunks(env, dEdbr)

    gradients(:, :) = gradients + derivs
    sigma(:, :) = sigma + localSigma

    !> contract with the Born radii derivatives
    call gemv(gradients, dbrdr, dEdbr, beta=1.0_dp)
    call gemv(sigma, dbrdL, dEdbr, beta=1.0_dp)

  end subroutine getBornEGCluster


  !> GB energy and gradient using Still interaction kernel
  subroutine getBornEGStillCluster(iAtFirst, iAtLast, coords, chargesPerAtom, &
      & bornRad, keps, energies, derivs, sigma, dEdbr)

    !> Number of atoms
    integer, intent(in) :: iAtFirst

    !> Number of atoms
    integer, intent(in) :: iAtLast

    !> Current atomic positions
    real(dp), intent(in) :: coords(:, :)

    !> Charges
    real(dp), intent(in) :: chargesPerAtom(:)

    !> Born radii
    real(dp), intent(in) :: bornRad(:)

    !> Dielectric screening
    real(dp), intent(in) :: keps

    !> Atom resolved energies
    real(dp), intent(inout) :: energies(:)

    !> Molecular gradient
    real(dp), intent(inout) :: derivs(:, :)

    !> Strain derivative
    real(dp), intent(inout) :: sigma(:, :)

    !> Strain derivative
    real(dp), intent(inout) :: dEdbr(:)

    integer :: iAt1, iAt2
    real(dp) :: aa, dist2, fgb2, qq, dd, expd, dfgb, dfgb2, dfgb3, ap, bp
    real(dp) :: grddbi, grddbj, vec(3), dGr(3), dSr(3, 3)

    !$omp parallel do default(none) reduction(+:energies, derivs, dEdbr, sigma) &
    !$omp shared(iAtFirst, iAtLast, coords, chargesPerAtom, bornRad, kEps) &
    !$omp private(iAt1, iAt2, aa, dist2, fgb2, qq, dd, expd, dfgb, dfgb2, dfgb3) &
    !$omp private(ap, bp, grddbi, grddbj, vec, dGr, dSr)
    do iAt1 = iAtFirst, iAtLast
      do iAt2 = 1, iAt1-1
        vec(:) = coords(:, iAt1) - coords(:, iAt2)
        dist2 = sum(vec**2)

        ! dielectric scaling of the charges
        qq = chargesPerAtom(iAt1)*chargesPerAtom(iAt2)
        aa = bornRad(iAt1)*bornRad(iAt2)
        dd = 0.25_dp*dist2/aa
        expd = exp(-dd)
        fgb2 = dist2+aa*expd
        dfgb2 = 1.0_dp/fgb2
        dfgb = sqrt(dfgb2)
        dfgb3 = dfgb2*dfgb*keps

        energies(iAt1) = energies(iAt1) + qq*keps*dfgb/2
        if (iAt1 /= iAt2) then
          energies(iAt2) = energies(iAt2) + qq*keps*dfgb/2
        end if

        ap = (1.0_dp-0.25_dp*expd)*dfgb3
        dGr = ap*vec
        derivs(:,iAt1) = derivs(:,iAt1) - dGr*qq
        derivs(:,iAt2) = derivs(:,iAt2) + dGr*qq

        dSr = spread(dGr, 1, 3) * spread(vec, 2, 3)
        if (iAt1 /= iAt2) then
          sigma = sigma + dSr
        else
          sigma = sigma + dSr/2
        end if

        bp = -0.5_dp*expd*(1.0_dp+dd)*dfgb3
        grddbi = bornRad(iAt2)*bp
        grddbj = bornRad(iAt1)*bp
        dEdbr(iAt1) = dEdbr(iAt1) + grddbi*qq
        if (iAt1 /= iAt2) then
          dEdbr(iAt2) = dEdbr(iAt2) + grddbj*qq
        end if

      end do
    end do
    !$omp end parallel do

  end subroutine getBornEGStillCluster


  !> GB energy and gradient using P16 interaction kernel
  subroutine getBornEGP16Cluster(iAtFirst, iAtLast, coords, chargesPerAtom, &
      & bornRad, keps, energies, derivs, sigma, dEdbr)

    !> Number of atoms
    integer, intent(in) :: iAtFirst

    !> Number of atoms
    integer, intent(in) :: iAtLast

    !> Current atomic positions
    real(dp), intent(in) :: coords(:, :)

    !> Charges
    real(dp), intent(in) :: chargesPerAtom(:)

    !> Born radii
    real(dp), intent(in) :: bornRad(:)

    !> Dielectric screening
    real(dp), intent(in) :: keps

    !> Atom resolved energies
    real(dp), intent(inout) :: energies(:)

    !> Molecular gradient
    real(dp), intent(inout) :: derivs(:, :)

    !> Strain derivative
    real(dp), intent(inout) :: sigma(:, :)

    !> Strain derivative
    real(dp), intent(inout) :: dEdbr(:)

    integer :: iAt1, iAt2
    real(dp) :: vec(3), r2, r1, ab, arg1, arg16, qq, fgb, dfgb, dfgb2
    real(dp) :: dEdbr1, dEdbr2, dG(3), ap, bp, dS(3, 3)

    !$omp parallel do default(none) reduction(+:energies, derivs, dEdbr, sigma) &
    !$omp shared(iAtFirst, iAtLast, coords, chargesPerAtom, bornRad, kEps) &
    !$omp private(iAt1, iAt2, vec, r1, r2, ab, arg1, arg16, fgb, dfgb, dfgb2, ap, &
    !$omp& bp, qq, dEdbr1, dEdbr2, dG, dS)
    do iAt1 = iAtFirst, iAtLast
      do iAt2 = 1, iAt1-1
        vec(:) = coords(:, iAt1) - coords(:, iAt2)
        r2 = sum(vec**2)
        r1 = sqrt(r2)
        qq = chargesPerAtom(iAt1)*chargesPerAtom(iAt2)

        ab = sqrt(bornRad(iAt1) * bornRad(iAt2))
        arg1 = ab / (ab + zetaP16o16*r1) ! 1 / (1 + ζR/(16·ab))
        arg16 = arg1 * arg1 ! 1 / (1 + ζR/(16·ab))²
        arg16 = arg16 * arg16 ! 1 / (1 + ζR/(16·ab))⁴
        arg16 = arg16 * arg16 ! 1 / (1 + ζR/(16·ab))⁸
        arg16 = arg16 * arg16 ! 1 / (1 + ζR/(16·ab))¹⁶

        fgb = r1 + ab*arg16
        dfgb = 1.0_dp / fgb
        dfgb2 = dfgb * dfgb

        energies(iAt1) = energies(iAt1) + qq*keps*dfgb/2
        if (iAt1 /= iAt2) then
          energies(iAt2) = energies(iAt2) + qq*keps*dfgb/2
        end if

        ! (1 - ζ/(1 + Rζ/(16 ab))^17)/(R + ab/(1 + Rζ/(16 ab))¹⁶)²
        ap = (1.0_dp - zetaP16 * arg1 * arg16) * dfgb2
        dG(:) = ap * vec * kEps / r1 * qq
        derivs(:, iAt1) = derivs(:, iAt1) - dG
        derivs(:, iAt2) = derivs(:, iAt2) + dG

        dS = spread(dG, 1, 3) * spread(vec, 2, 3)
        if (iAt1 /= iAt2) then
          sigma = sigma + dS
        else
          sigma = sigma + dS/2
        end if

        ! -(Rζ/(2·ab²·(1 + Rζ/(16·ab))¹⁷) + 1/(2·ab·(1 + Rζ/(16·ab))¹⁶))/(R + ab/(1 + Rζ/(16·ab))¹⁶)²
        bp = -0.5_dp*(r1 * zetaP16 / ab * arg1 + 1.0_dp) / ab * arg16 * dfgb2
        dEdbr1 = bornRad(iAt2) * bp * kEps * qq
        dEdbr2 = bornRad(iAt1) * bp * kEps * qq
        dEdbr(iAt1) = dEdbr(iAt1) + dEdbr1
        dEdbr(iAt2) = dEdbr(iAt2) + dEdbr2

      end do
    end do
    !$omp end parallel do

  end subroutine getBornEGP16Cluster


  !> Evaluate inertia tensor for solid spheres with mass rad**3
  pure subroutine getInertia(nAtom, coord, species, rad, center, inertia)

    !> Number of atoms
    integer, intent(in) :: nAtom

    !> Cartesian coordinates
    real(dp), intent(in) :: coord(:, :)

    !> Species identifiers for each atom
    integer, intent(in) :: species(:)

    !> Atomic radii
    real(dp), intent(in) :: rad(:)

    !> Center of mass
    real(dp), intent(in) :: center(:)

    !> Inertia tensor
    real(dp), intent(out) :: inertia(:, :)

    integer :: iAt, iSp
    real(dp) :: r2, rad2, rad3, vec(3)
    real(dp), parameter :: tof = 2.0_dp/5.0_dp, unity(3, 3) = reshape(&
        & [1.0_dp, 0.0_dp, 0.0_dp, 0.0_dp, 1.0_dp, 0.0_dp, 0.0_dp, 0.0_dp, 1.0_dp], &
        & [3, 3])

    inertia(:, :) = 0.0_dp
    do iAt = 1, nAtom
      iSp = species(iAt)
      rad2 = rad(iSp) * rad(iSp)
      rad3 = rad2 * rad(iSp)
      vec(:) = coord(:, iAt) - center
      r2 = sum(vec**2)
      inertia(:, :) = inertia + rad3 * ((r2 + tof*rad2) * unity &
          & - spread(vec, 1, 3) * spread(vec, 2, 3))
    end do

  end subroutine getInertia


  !> Molecular shape descriptor
  subroutine getADet(nAtom, coord, species, rad, aDet)

    !> Number of atoms
    integer, intent(in) :: nAtom

    !> Cartesian coordinates
    real(dp), intent(in) :: coord(:, :)

    !> Species identifiers for each atom
    integer, intent(in) :: species(:)

    !> Atomic radii
    real(dp), intent(in) :: rad(:)

    !> Shape descriptor of the structure
    real(dp), intent(out) :: aDet

    integer :: iAt, iSp
    real(dp) :: rad2, rad3, totRad3, center(3), inertia(3, 3)

    totRad3 = 0.0_dp
    center(:) = 0.0_dp
    do iAt = 1, nAtom
      iSp = species(iAt)
      rad2 = rad(iSp) * rad(iSp)
      rad3 = rad2 * rad(iSp)
      totRad3 = totRad3 + rad3
      center(:) = center + coord(:, iAt) * rad3
    end do
    center = center / totRad3

    call getInertia(nAtom, coord, species, rad, center, inertia)

    aDet = sqrt(determinant33(inertia)**(1.0_dp/3.0_dp)/(2.0_dp*totRad3)*5.0_dp)

  end subroutine getADet


  !> Derivative of the molecular shape descriptor
  subroutine getADetDeriv(nAtom, coord, species, rad, kEps, qvec, gradient)

    !> Number of atoms
    integer, intent(in) :: nAtom

    !> Cartesian coordinates
    real(dp), intent(in) :: coord(:, :)

    !> Species identifiers for each atom
    integer, intent(in) :: species(:)

    !> Atomic radii
    real(dp), intent(in) :: rad(:)

    !> Dielectric constant, including alpha times beta
    real(dp), intent(in) :: kEps

    !> Atomic gross charges
    real(dp), intent(in) :: qvec(:)

    !> Molecular gradient
    real(dp), intent(inout) :: gradient(:, :)

    integer :: iAt, iSp
    real(dp) :: rad2, rad3, totRad3, vec(3), center(3), inertia(3, 3), aDet
    real(dp) :: aDeriv(3, 3), qtotal

    qtotal = 0.0_dp
    totRad3 = 0.0_dp
    center(:) = 0.0_dp
    do iAt = 1, nAtom
      iSp = species(iAt)
      rad2 = rad(iSp) * rad(iSp)
      rad3 = rad2 * rad(iSp)
      totRad3 = totRad3 + rad3
      center(:) = center + coord(:, iAt) * rad3
      qtotal = qtotal + qvec(iAt)
    end do
    center = center / totRad3

    call getInertia(nAtom, coord, species, rad, center, inertia)

    aDet = sqrt(determinant33(inertia)**(1.0_dp/3.0_dp)/(2.0_dp*totRad3)*5.0_dp)

    aDeriv(:, :) = reshape([&
        & inertia(1,1)*(inertia(2,2)+inertia(3,3))-inertia(1,2)**2-inertia(1,3)**2, &
        & inertia(1,2)*inertia(3,3)-inertia(1,3)*inertia(2,3), & ! xy
        & inertia(1,3)*inertia(2,2)-inertia(1,2)*inertia(3,2), & ! xz
        & inertia(1,2)*inertia(3,3)-inertia(1,3)*inertia(2,3), & ! xy
        & inertia(2,2)*(inertia(1,1)+inertia(3,3))-inertia(1,2)**2-inertia(2,3)**2, &
        & inertia(1,1)*inertia(2,3)-inertia(1,2)*inertia(1,3), & ! yz
        & inertia(1,3)*inertia(2,2)-inertia(1,2)*inertia(3,2), & ! xz
        & inertia(1,1)*inertia(2,3)-inertia(1,2)*inertia(1,3), & ! yz
        & inertia(3,3)*(inertia(1,1)+inertia(2,2))-inertia(1,3)**2-inertia(2,3)**2],&
        & shape=[3, 3]) * (250.0_dp / (48.0_dp * totRad3**3 * aDet**5)) &
        & * (-0.5_dp * kEps * qtotal**2 / aDet**2)

    do iAt = 1, nAtom
      iSp = species(iAt)
      rad2 = rad(iSp) * rad(iSp)
      rad3 = rad2 * rad(iSp)
      vec(:) = coord(:, iAt) - center
      gradient(:, iAt) = gradient(:, iAt) + rad3 * matmul(aDeriv, vec)
    end do

  end subroutine getADetDeriv


end module dftbp_born
