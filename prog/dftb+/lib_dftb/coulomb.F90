!--------------------------------------------------------------------------------------------------!
!  DFTB+: general package for performing fast atomistic simulations                                !
!  Copyright (C) 2006 - 2020  DFTB+ developers group                                               !
!                                                                                                  !
!  See the LICENSE file for terms of usage and distribution.                                       !
!--------------------------------------------------------------------------------------------------!

#:include 'common.fypp'

!> Contains routines to calculate the coulombic interaction in non periodic and periodic systems.
module dftbp_coulomb
  use dftbp_assert
  use dftbp_accuracy
  use dftbp_blasroutines, only : hemv
  use dftbp_commontypes, only : TOrbitals
  use dftbp_constants, only : pi
  use dftbp_dynneighlist
  use dftbp_environment
  use dftbp_errorfunction
  use dftbp_message
#:if WITH_MPI
  use dftbp_mpifx
#:endif
  use dftbp_periodic, only : TNeighbourList, getLatticePoints, getCellTranslations
#:if WITH_SCALAPACK
  use dftbp_scalapackfx
#:endif
  use dftbp_schedule
  implicit none

  private

  public :: boundaryCondition, TCoulombInput, TCoulombCont, init
  public :: invRCluster, invRPeriodic, sumInvR, addInvRPrime, getOptimalAlphaEwald, getMaxGEwald
  public :: getMaxREwald, invRStress
  public :: addInvRPrimeXlbomd
  public :: ewaldReal, ewaldReciprocal, derivEwaldReal, derivEwaldReciprocal, derivStressEwaldRec


  !> Possible boundary conditions for evaluating coulombic interactions.
  type :: TBoundaryConditionEnum

    !> Finite molecular cluster
    integer :: cluster = 0

    !> Three dimensional infinite periodic boundary conditions
    integer :: pbc3d = 3

  end type TBoundaryConditionEnum

  !> Actual instance of the boundary condition enumerator
  type(TBoundaryConditionEnum), parameter :: boundaryCondition = TBoundaryConditionEnum()


  !> Input data for coulombic interaction container
  type :: TCoulombInput

    !> if > 0 -> manual setting for alpha
    real(dp) :: ewaldAlpha = 0.0_dp

    !> Ewald tolerance
    real(dp) :: tolEwald = 0.0_dp

  end type TCoulombInput


  !> Calculates the coulombic interaction in non-periodic and periodic systems.
  type :: TCoulombCont

    !> number of atoms
    integer :: nAtom

    !> Boundary condition for coulombic interaction evaluation
    integer :: boundaryCondition

    !> lattice vectors if periodic
    real(dp) :: latVecs(3, 3)

    !> reciprocal lattice vectors
    real(dp) :: recVecs(3, 3)

    !> Cell volume
    real(dp) :: volume

    !> Lattice points for reciprocal Ewald
    real(dp), allocatable :: gLatPoint(:,:)

    !> Real lattice points for asymmetric Ewald sum
    real(dp), allocatable :: rCellVec(:,:)

    !> Dynamic neighbour list for the real space Ewald summation
    type(TDynNeighList), allocatable :: neighListGen

    !> evaluate Ewald parameter
    logical :: tAutoEwald

    !> Parameter for Ewald
    real(dp) :: alpha

    !> Ewald tolerance
    real(dp) :: tolEwald

    !> are the coordinates current?
    logical :: tCoordsUpdated

    !> are the charges current?
    logical :: tChargesUpdated

    !> Stores 1/r between atom pairs
    real(dp), allocatable :: invRMat(:,:)

    !> Shift vector per atom
    real(dp), allocatable :: shiftPerAtom(:)

    !> Negative gross charge per atom
    real(dp), allocatable :: deltaQAtom(:)

#:if WITH_SCALAPACK
    !> Descriptor for 1/R matrix
    integer :: descInvRMat(DLEN_)

    !> Descriptor for charge vector
    integer :: descQVec(DLEN_)

    !> Distributed shift per atom vector
    real(dp), allocatable :: shiftPerAtomGlobal(:,:)

    !> Distributed charge vector
    real(dp), allocatable :: qGlobal(:,:)
#:endif

  contains

    !> update internal copy of coordinates
    procedure :: updateCoords

    !> update internal copy of lattice vectors
    procedure :: updateLatVecs

    !> get energy contributions
    procedure :: addEnergy

    !> get force contributions
    procedure :: addGradients

    !> get stress tensor contributions
    procedure :: addStress

    !> Updates with changed charges for the instance
    procedure :: updateCharges

    !> Update potential shifts for the instance
    procedure :: updateShifts

    !> Returns shifts per atom
    procedure :: addShiftPerAtom

    !> Returns shifts per shell
    procedure :: addShiftPerShell

    !> Get the variables relate to periodic information
    procedure :: getPeriodicInfo

    !> Calculates the -1/R**2 deriv contribution for all atoms for the non-periodic case, without
    !> storing anything.
    procedure :: addInvRPrimeClusterMat

    !> Calculates the -1/R**2 deriv contribution for the periodic case, without storing anything.
    procedure :: addInvRPrimePeriodicMat

  end type TCoulombCont


  !> Initialize coulombic interaction container from input data
  interface init
    module procedure :: initialize
  end interface init


  !> 1/r interaction for all atoms with another group
  interface sumInvR
    module procedure sumInvRClusterAsymm
    module procedure sumInvRPeriodicAsymm
  end interface sumInvR


  !> 1/r^2
  interface addInvRPrime
    module procedure addInvRPrimeCluster
    module procedure addInvRPrimeClusterAsymm
    module procedure addInvRPrimePeriodic
    module procedure addInvRPrimePeriodicAsymm
  end interface addInvRPrime


  !> 1/r^2 term for extended lagrangian
  interface addInvRPrimeXlbomd
    module procedure addInvRPrimeXlbomdCluster
    module procedure addInvRPrimeXlbomdPeriodic
  end interface addInvRPrimeXlbomd

  !> Maximal argument value of erf, after which it is constant
  real(dp), parameter :: erfArgLimit = 10.0_dp

  !> Chunk size to use when obtaining neighbours dynamically via an iterator
  integer, parameter :: iterChunkSize = 1000


contains


  subroutine initialize(this, input, env, nAtom, latVecs, recVecs, volume)

    !> Data structure
    class(TCoulombCont), intent(out) :: this

    !> Input data for coulombic interaction container
    class(TCoulombInput), intent(in) :: input

    !> Environment settings
    type(TEnvironment), intent(in) :: env

    !> Nr. of atoms in the system
    integer, intent(in) :: nAtom

    !> Lattice vectors, if the system is periodic
    real(dp), intent(in), optional :: latVecs(:,:)

    !> reciprocal lattice vectors
    real(dp), intent(in), optional :: recVecs(:,:)

    !> Cell volume
    real(dp), intent(in), optional :: volume

    real(dp) :: maxREwald, maxGEwald
  #:if WITH_SCALAPACK
    integer :: nRowLoc, nColLoc
  #:endif

    @:ASSERT(present(latVecs) .eqv. present(recVecs))
    @:ASSERT(present(recVecs) .eqv. present(volume))

    this%nAtom = nAtom

    if (present(latVecs) .and. present(recVecs) .and. present(volume)) then
      this%boundaryCondition = boundaryCondition%pbc3d
      this%latVecs(:, :) = latVecs
      this%recVecs(:, :) = recVecs
      this%volume = volume

      ! Initialize Ewald summation for the periodic case
      this%tAutoEwald = input%ewaldAlpha <= 0.0_dp
      this%tolEwald = input%tolEwald
      if (this%tAutoEwald) then
        this%alpha = getOptimalAlphaEwald(latVecs, recVecs, volume, this%tolEwald)
      else
        this%alpha = input%ewaldAlpha
      end if
      maxREwald = getMaxREwald(this%alpha, this%tolEwald)
      maxGEwald = getMaxGEwald(this%alpha, volume, this%tolEwald)
      call getLatticePoints(this%gLatPoint, recVecs, latVecs/(2.0_dp*pi), maxGEwald,&
          & onlyInside=.true., reduceByInversion=.true., withoutOrigin=.true.)
      this%gLatPoint(:,:) = matmul(recVecs, this%gLatPoint)

      allocate(this%neighListGen)
      call TDynNeighList_init(this%neighListGen, maxREwald, nAtom, .true.)
    else
      this%boundaryCondition = boundaryCondition%cluster
    end if

  #:if WITH_SCALAPACK
    if (env%blacs%atomGrid%iproc /= -1) then
      call scalafx_getdescriptor(env%blacs%atomGrid, nAtom, nAtom,&
          & env%blacs%rowBlockSize, env%blacs%columnBlockSize, this%descInvRMat)
      call scalafx_getlocalshape(env%blacs%atomGrid, this%descInvRMat, nRowLoc, nColLoc)
      allocate(this%invRMat(nRowLoc, nColLoc))
      call scalafx_getdescriptor(env%blacs%atomGrid, 1, nAtom, env%blacs%rowBlockSize,&
          & env%blacs%columnBlockSize, this%descQVec)
      call scalafx_getlocalshape(env%blacs%atomGrid, this%descQVec, nRowLoc, nColLoc)
      allocate(this%shiftPerAtomGlobal(nRowLoc, nColLoc))
      allocate(this%qGlobal(nRowLoc, nColLoc))
    end if
  #:else
    allocate(this%invRMat(nAtom, nAtom))
  #:endif

    ! Initialise arrays for charge differences
    allocate(this%deltaQAtom(nAtom))

    ! Initialise arrays for potential shifts
    allocate(this%shiftPerAtom(nAtom))

  end subroutine initialize


  !> Update internal stored coordinates
  subroutine updateCoords(this, env, neighList, coords, species)

    !> Data structure
    class(TCoulombCont), intent(inout) :: this

    !> Computational environment settings
    type(TEnvironment), intent(in) :: env

    !> List of neighbours to atoms
    type(TNeighbourList), intent(in) :: neighList

    !> Atomic coordinates
    real(dp), intent(in) :: coords(:,:)

    !> Central cell chemical species
    integer, intent(in) :: species(:)

    if (this%boundaryCondition == boundaryCondition%pbc3d) then
      call this%neighListGen%updateCoords(coords(:, 1:this%nAtom))
    end if

    ! If process is outside of atom grid, skip invRMat calculation
    if (allocated(this%invRMat)) then
      if (this%boundaryCondition == boundaryCondition%pbc3d) then
        call invRPeriodic(env, this%nAtom, coords, this%neighListGen, this%gLatPoint,&
            & this%alpha, this%volume, this%invRMat)
      else
        call invRCluster(env, this%nAtom, coords, this%invRMat)
      end if
    end if

    this%tCoordsUpdated = .true.
    this%tChargesUpdated = .false.

  end subroutine updateCoords


  !> Update internal copy of lattice vectors
  subroutine updateLatVecs(this, latVecs, recVecs, volume)

    !> Data structure
    class(TCoulombCont), intent(inout) :: this

    !> New lattice vectors
    real(dp), intent(in) :: latVecs(:,:)

    !> New reciprocal lattice vectors
    real(dp), intent(in) :: recVecs(:,:)

    !> New volume
    real(dp), intent(in) :: volume

    real(dp) :: maxREwald, maxGEwald

    real(dp), allocatable :: dummy(:,:)

    @:ASSERT(all(shape(latVecs) == shape(this%latVecs)))

    if (this%tAutoEwald) then
      this%alpha = getOptimalAlphaEwald(latVecs, recVecs, volume, this%tolEwald)
      maxREwald = getMaxREwald(this%alpha, this%tolEwald)
    end if
    maxGEwald = getMaxGEwald(this%alpha, volume, this%tolEwald)
    call getLatticePoints(this%gLatPoint, recVecs, latVecs/(2.0_dp*pi), maxGEwald,&
        &onlyInside=.true., reduceByInversion=.true., withoutOrigin=.true.)
    this%gLatPoint = matmul(recVecs, this%gLatPoint)

    this%latVecs(:, :) = latVecs
    this%recVecs(:, :) = recVecs
    this%volume = volume

    ! Fold charges back to unit cell
    call getCellTranslations(dummy, this%rCellVec, latVecs, recVecs / (2.0_dp * pi), maxREwald)

    call this%neighListGen%updateLatVecs(latVecs, recVecs / (2.0_dp * pi))

    this%tCoordsUpdated = .false.
    this%tChargesUpdated = .false.

  end subroutine updateLatVecs


  !> Get energy contributions from coulombic interactions
  subroutine addEnergy(this, energies, dQOut, dQOutAtom, dQOutShell)

    !> Data structure
    class(TCoulombCont), intent(in) :: this

    !> Energy contributions for each atom
    real(dp), intent(inout) :: energies(:)

    !> Negative gross charge (present for XLBOMD)
    real(dp), intent(in), optional :: dQOut(:,:)

    !> Negative gross charge per atom (present for XLBOMD)
    real(dp), intent(in), optional :: dQOutAtom(:)

    !> Negative gross charge per shell (present for XLBOMD)
    real(dp), intent(in), optional :: dQOutShell(:,:)

    @:ASSERT(this%tCoordsUpdated)
    @:ASSERT(this%tChargesUpdated)
    @:ASSERT(present(dQOut) .eqv. present(dQOutAtom))
    @:ASSERT(present(dQOut) .eqv. present(dQOutShell))
    @:ASSERT(size(energies) == this%nAtom)

    if (present(dQOutAtom)) then
      ! XLBOMD: 1/2 sum_A (2 q_A - n_A) * shift(n_A)
      energies(:) = energies + 0.5_dp * this%shiftPerAtom * (2.0_dp * dQOutAtom &
          & - this%deltaQAtom)
    else
      energies(:) = energies + 0.5_dp * this%shiftPerAtom * this%deltaQAtom
    end if

  end subroutine addEnergy


  !> Get force contributions
  subroutine addGradients(this, env, coords, species, iNeighbour, img2CentCell, &
      & gradients, dQOut, dQOutAtom, dQOutShell)

    !> Data structure
    class(TCoulombCont), intent(in) :: this

    !> Environment settings
    type(TEnvironment), intent(in) :: env

    !> Atomic coordinates
    real(dp), intent(in) :: coords(:,:)

    !> Species for each atom
    integer, intent(in) :: species(:)

    !> List of neighbours for each atom
    integer, intent(in) :: iNeighbour(0:,:)

    !> Indexing of images of the atoms in the central cell
    integer, intent(in) :: img2CentCell(:)

    !> Gradient contributions
    real(dp), intent(inout) :: gradients(:,:)

    !> Negative gross charge (present for XLBOMD)
    real(dp), intent(in), optional :: dQOut(:,:)

    !> Negative gross charge per atom (present for XLBOMD)
    real(dp), intent(in), optional :: dQOutAtom(:)

    !> Negative gross charge per shell (present for XLBOMD)
    real(dp), intent(in), optional :: dQOutShell(:,:)

    @:ASSERT(this%tCoordsUpdated)
    @:ASSERT(this%tChargesUpdated)
    @:ASSERT(present(dQOut) .eqv. present(dQOutAtom))
    @:ASSERT(present(dQOut) .eqv. present(dQOutShell))
    @:ASSERT(all(shape(gradients) == [3, this%nAtom]))

    ! 1/R contribution
    if (present(dQOutAtom)) then
      if (this%boundaryCondition == boundaryCondition%pbc3d) then
        call addInvRPrimeXlbomd(env, this%nAtom, coords, this%neighListGen, &
            & this%gLatPoint, this%alpha, this%volume, this%deltaQAtom, &
            & dQOutAtom, gradients)
      else
        call addInvRPrimeXlbomd(env, this%nAtom, coords, this%deltaQAtom, &
            & dQOutAtom, gradients)
      end if
    else
      if (this%boundaryCondition == boundaryCondition%pbc3d) then
        call addInvRPrime(env, this%nAtom, coords, this%neighListGen, &
            & this%gLatPoint, this%alpha, this%volume, this%deltaQAtom, gradients)
      else
        call addInvRPrime(env, this%nAtom, coords, this%deltaQAtom, gradients)
      end if
    end if

  end subroutine addGradients


  !> Get stress tensor contributions
  subroutine addStress(this, env, coords, species, iNeighbour, img2CentCell, &
      & stress)

    !> Data structure
    class(TCoulombCont), intent(in) :: this

    !> Computational environment settings
    type(TEnvironment), intent(in) :: env

    !> Atomic coordinates
    real(dp), intent(in) :: coords(:,:)

    !> Species for each atom.
    integer, intent(in) :: species(:)

    !> List of neighbours for each atom.
    integer, intent(in) :: iNeighbour(0:,:)

    !> Indexing of images of the atoms in the central cell.
    integer, intent(in) :: img2CentCell(:)

    !> Stress tensor contributions
    real(dp), intent(inout) :: stress(:,:)

    real(dp) :: stTmp(3,3)

    @:ASSERT(this%tCoordsUpdated)
    @:ASSERT(this%tChargesUpdated)
    @:ASSERT(all(shape(stress) == [3, 3]))

    ! 1/R contribution
    stTmp = 0.0_dp
    call invRStress(env, this%nAtom, coords, this%neighListGen, this%gLatPoint, this%alpha,&
        & this%volume, this%deltaQAtom, stTmp)

    stress(:,:) = stress(:,:) - 0.5_dp * stTmp(:,:)

  end subroutine addStress


  !> Updates with changed charges for the instance.
  subroutine updateCharges(this, env, qOrbital, q0, orb, species, deltaQ, &
        & deltaQAtom, deltaQPerLShell, deltaQUniqU)

    !> Data structure
    class(TCoulombCont), intent(inout) :: this

    !> Environment settings
    type(TEnvironment), intent(in) :: env

    !> Orbital resolved charges
    real(dp), intent(in) :: qOrbital(:,:,:)

    !> Reference charge distribution (neutral atoms)
    real(dp), intent(in) :: q0(:,:,:)

    !> Contains information about the atomic orbitals in the system
    type(TOrbitals), intent(in) :: orb

    !> Species, shape: [nAtom]
    integer, intent(in) :: species(:)

    !> Negative gross charge
    real(dp), intent(in) :: deltaQ(:,:)

    !> Negative gross charge per shell
    real(dp), intent(in) :: deltaQPerLShell(:,:)

    !> Negative gross charge per atom
    real(dp), intent(in) :: deltaQAtom(:)

    !> Negative gross charge per U
    real(dp), intent(in) :: deltaQUniqU(:,:)

    @:ASSERT(this%tCoordsUpdated)

    this%deltaQAtom(:) = deltaQAtom

    this%tChargesUpdated = .true.

  end subroutine updateCharges


  !> Update potential shifts. Call after updateCharges
  subroutine updateShifts(this, env, orb, species, iNeighbour, img2CentCell)

    !> Data structure
    class(TCoulombCont), intent(inout), target :: this

    !> Environment settings
    type(TEnvironment), intent(in) :: env

    !> Contains information about the atomic orbitals in the system
    type(TOrbitals), intent(in) :: orb

    !> Species of the atoms (should not change during run)
    integer, intent(in) :: species(:)

    !> Neighbor indexes
    integer, intent(in) :: iNeighbour(0:,:)

    !> Mapping on atoms in the central cell
    integer, intent(in) :: img2CentCell(:)

  #:if WITH_SCALAPACK
    real(dp), pointer :: deltaQAtom2D(:,:), shiftPerAtom2D(:,:)
    integer :: ll
  #:endif

    @:ASSERT(this%tCoordsUpdated)
    @:ASSERT(this%tChargesUpdated)

    this%shiftPerAtom(:) = 0.0_dp

  #:if WITH_SCALAPACK
    if (env%blacs%atomGrid%iproc /= -1) then
      ll = size(this%deltaQAtom)
      deltaQAtom2D(1:1, 1:ll) => this%deltaQAtom
      ll = size(this%shiftPerAtom)
      shiftPerAtom2D(1:1, 1:ll) => this%shiftPerAtom
      call scalafx_cpl2g(env%blacs%atomGrid, deltaQAtom2D, this%descQVec, 1, 1, &
          & this%qGlobal)
      call pblasfx_psymv(this%invRMat, this%descInvRMat, this%qGlobal, &
          & this%descQVec, this%shiftPerAtomGlobal, this%descQVec)
      call scalafx_cpg2l(env%blacs%atomGrid, this%descQVec, 1, 1, &
          & this%shiftPerAtomGlobal, shiftPerAtom2D)
    end if
    call mpifx_allreduceip(env%mpi%groupComm, this%shiftPerAtom, MPI_SUM)
  #:else
    call hemv(this%shiftPerAtom, this%invRMat, this%deltaQAtom, 'L')
  #:endif

  end subroutine updateShifts


  !> Returns shifts per atom
  subroutine addShiftPerAtom(this, shiftPerAtom)

    !> Data structure
    class(TCoulombCont), intent(in) :: this

    !> Shift per atom
    real(dp), intent(inout) :: shiftPerAtom(:)

    @:ASSERT(this%tCoordsUpdated)
    @:ASSERT(this%tChargesUpdated)
    @:ASSERT(size(shiftPerAtom) == this%nAtom)

    shiftPerAtom(:) = shiftPerAtom + this%shiftPerAtom

  end subroutine addShiftPerAtom


  !> Returns shifts per atom
  subroutine addShiftPerShell(this, shiftPerShell)

    !> Data structure
    class(TCoulombCont), intent(in) :: this

    !> Shift per shell
    real(dp), intent(inout) :: shiftPerShell(:,:)

    @:ASSERT(this%tCoordsUpdated)
    @:ASSERT(this%tChargesUpdated)
    @:ASSERT(size(shiftPerShell, dim=2) == this%nAtom)

  end subroutine addShiftPerShell


  !> Calculates the 1/R Matrix for all atoms for the non-periodic case.  Only the lower triangle is
  !> constructed.
  subroutine invRCluster(env, nAtom, coord, invRMat)

    !> Computational environment settings
    type(TEnvironment), intent(in) :: env

    !> number of atoms
    integer, intent(in) :: nAtom

    !> List of atomic coordinates.
    real(dp), intent(in) :: coord(:,:)

    !> Matrix of 1/R values for each atom pair.
    real(dp), intent(out) :: invRMat(:,:)

  #:if WITH_SCALAPACK
    integer :: descInvRMat(DLEN_)
  #:endif

   #:if WITH_SCALAPACK
    if (env%blacs%atomGrid%iproc == -1) then
      return
    end if
    call scalafx_getdescriptor(env%blacs%atomGrid, nAtom, nAtom, env%blacs%rowBlockSize,&
        & env%blacs%columnBlockSize, descInvRMat)
    call invRClusterBlacs(env%blacs%atomGrid, coord, descInvRMat, invRMat)
  #:else
    call invRClusterSerial(coord, invRMat)
  #:endif

  end subroutine invRCluster


#:if WITH_SCALAPACK

  !> Calculates the 1/R Matrix for all atoms for the non-periodic case (BLACS version)
  !>
  !> Note: Only the lower triangle is constructed.
  !>
  subroutine invRClusterBlacs(grid, coord, descInvRMat, invRMat)

    !> Grid to use for the computation
    type(blacsgrid), intent(in) :: grid

    !> List of atomic coordinates.
    real(dp), intent(in) :: coord(:,:)

    !> Descriptor of the distributed matrix
    integer, intent(in) :: descInvRMat(DLEN_)

    !> Matrix of 1/R values for each atom pair.
    real(dp), intent(out) :: invRMat(:,:)

    integer :: ii, jj, iAt1, iAt2
    real(dp) :: dist, vect(3)

    invRMat(:,:) = 0.0_dp

    !$OMP PARALLEL DO&
    !$OMP& DEFAULT(SHARED) PRIVATE(ii, iAt1, iAt2, vect, dist) SCHEDULE(RUNTIME)
    do jj = 1, size(invRMat, dim=2)
      iAt1 = scalafx_indxl2g(jj, descInvRMat(NB_), grid%mycol, descInvRMat(CSRC_), grid%ncol)
      do ii = 1, size(invRMat, dim=1)
        iAt2 = scalafx_indxl2g(ii, descInvRMat(MB_), grid%myrow, descInvRMat(RSRC_), grid%nrow)
        if (iAt2 <= iAt1) then
          ! wrong triangle
          cycle
        end if
        vect(:) = coord(:,iAt1) - coord(:,iAt2)
        dist = sqrt(sum(vect**2))
        invRMat(ii, jj) = 1.0_dp / dist
      end do
    end do
    !$OMP END PARALLEL DO

  end subroutine invRClusterBlacs


#:else


  !> Calculates the 1/R Matrix for all atoms for the non-periodic case (serial version)
  !>
  !> Note: Only the lower triangle is constructed.
  !>
  subroutine invRClusterSerial(coord, invRMat)

    !> List of atomic coordinates.
    real(dp), intent(in) :: coord(:,:)

    !> Matrix of 1/R values for each atom pair.
    real(dp), intent(out) :: invRMat(:,:)

    integer :: ii, jj, iAt1, iAt2
    real(dp) :: dist, vect(3)

    invRMat(:,:) = 0.0_dp

    !$OMP PARALLEL DO&
    !$OMP& DEFAULT(SHARED) PRIVATE(iAt1, ii, iAt2, vect, dist) SCHEDULE(RUNTIME)
    do jj = 1, size(invRMat, dim=2)
      iAt1 = jj
      do ii = jj + 1, size(invRMat, dim=1)
        iAt2 = ii
        vect(:) = coord(:,iAt1) - coord(:,iAt2)
        dist = sqrt(sum(vect**2))
        invRMat(ii, jj) = 1.0_dp / dist
      end do
    end do
    !$OMP END PARALLEL DO

  end subroutine invRClusterSerial

#:endif


  !> Calculates the summed 1/R vector for all atoms for the non-periodic case asymmmetric case (like
  !> interaction of atoms with point charges).
  subroutine sumInvRClusterAsymm(env, nAtom0, nAtom1, coord0, coord1, charges1, invRVec,&
      & blurWidths1, epsSoften)

    !> Computational environment settings
    type(TEnvironment), intent(in) :: env

    !> Number of atoms in the first group
    integer, intent(in) :: nAtom0

    !> Number of atoms in the second group
    integer, intent(in) :: nAtom1

    !> Coordinates of the first group of objects (atoms)
    real(dp), intent(in) :: coord0(:,:)

    !> Coordinates of the 2nd group of objects (point charges)
    real(dp), intent(in) :: coord1(:,:)

    !> Charges of the 2nd group of objects
    real(dp), intent(in) :: charges1(:)

    !> Vector of sum_i q_i/|R_atom - R_i] values for each atom
    real(dp), intent(out) :: invRVec(:)

    !> Gaussian blur width of the charges in the 2nd group
    real(dp), intent(in), optional :: blurWidths1(:)

    !> Short distance softening
    real(dp), intent(in), optional :: epsSoften

    integer :: iAt0, iAt1
    real(dp) :: dist, vect(3), fTmp, epsSoften2
    integer :: iAtFirst0, iAtLast0, iAtFirst1, iAtLast1

    if (present(epsSoften)) then
      epsSoften2 = epsSoften**2
    else
      epsSoften2 = 0.0_dp
    end if

    call distributeRangeInChunks2(env, 1, nAtom0, 1, nAtom1, iAtFirst0, iAtLast0, iAtFirst1,&
        & iAtLast1)

    invRVec(:) = 0.0_dp

    ! Doing blurring and non blurring case separately in order to avoid the if branch in the loop
    if (present(blurWidths1)) then
      !$OMP PARALLEL DO&
      !$OMP& DEFAULT(SHARED) PRIVATE(iAt1, vect, dist, fTmp) SCHEDULE(RUNTIME)
      do iAt0 = iAtFirst0, iAtLast0
        do iAt1 = iAtFirst1, iAtLast1
          vect(:) = coord0(:,iAt0) - coord1(:,iAt1)
          dist = sum(vect**2)
          fTmp = charges1(iAt1) / sqrt(dist + epsSoften2)
          if (dist < (erfArgLimit * blurWidths1(iAt1))**2) then
            fTmp = fTmp * erfwrap(sqrt(dist) / blurWidths1(iAt1))
          end if
          invRVec(iAt0) = invRVec(iAt0) + fTmp
        end do
      end do
      !$OMP END PARALLEL DO
    else
      !$OMP PARALLEL DO&
      !$OMP& DEFAULT(SHARED) PRIVATE(iAt1, vect, dist) SCHEDULE(RUNTIME)
      do iAt0 = iAtFirst0, iAtLast0
        do iAt1 = iAtFirst1, iAtLast1
          vect(:) = coord0(:,iAt0) - coord1(:,iAt1)
          dist = sqrt(sum(vect**2) + epsSoften2)
          invRVec(iAt0) = invRVec(iAt0) + charges1(iAt1) / dist
        end do
      end do
      !$OMP END PARALLEL DO
    end if

    call assembleChunks(env, invRVec)

  end subroutine sumInvRClusterAsymm


  !> Calculates the 1/R Matrix for all atoms for the periodic case.  Only the lower triangle is
  !> constructed.
  subroutine invRPeriodic(env, nAtom, coord, neighList, recPoint, alpha, volume, invRMat)

    !> Computational environment settings
    type(TEnvironment), intent(in) :: env

    !> Number of atoms.
    integer, intent(in) :: nAtom

    !> List of atomic coordinates (all atoms).
    real(dp), intent(in) :: coord(:,:)

    !> Neighbour list for the real space summation of the Ewald.
    type(TDynNeighList), intent(in) :: neighList

    !> Contains the points included in the reciprocal sum.  The set should not include the origin or
    !> inversion related points.
    real(dp), intent(in) :: recPoint(:,:)

    !> Parameter for Ewald summation.
    real(dp), intent(in) :: alpha

    !> Volume of the real space unit cell.
    real(dp), intent(in) :: volume

    !> Matrix of 1/R values for each atom pair.
    real(dp), intent(out) :: invRMat(:,:)

  #:if WITH_SCALAPACK
    integer :: descInvRMat(DLEN_)
  #:endif

    @:ASSERT(volume > 0.0_dp)

    ! Somewhat redundant test, but stops compiler complaints as for serial case nAtom not referenced
    @:ASSERT(size(coord,dim=2) >= nAtom)

  #:if WITH_SCALAPACK
    if (env%blacs%atomGrid%iproc == -1) then
      return
    end if
    call scalafx_getdescriptor(env%blacs%atomGrid, nAtom, nAtom, env%blacs%rowBlockSize,&
        & env%blacs%columnBlockSize, descInvRMat)
    call invRPeriodicBlacs(env%blacs%atomGrid, coord, neighList, recPoint, alpha, volume,&
        & descInvRMat, invRMat)
  #:else
    call invRPeriodicSerial(coord, neighList, recPoint, alpha, volume, invRMat)
  #:endif

  end subroutine invRPeriodic


#:if WITH_SCALAPACK

  !> Calculates the 1/R Matrix for all atoms for the periodic case (BLACS version)
  !>
  !> Note: Only the lower triangle is constructed.
  !>
  subroutine invRPeriodicBlacs(grid, coord, neighList, recPoint, alpha, volume, descInvRMat,&
      & invRMat)

    !> Grid to use for the computation
    type(blacsgrid), intent(in) :: grid

    !> List of atomic coordinates (all atoms).
    real(dp), intent(in) :: coord(:,:)

    !> Neighbour list for the real space summation of the Ewald.
    type(TDynNeighList), target, intent(in) :: neighList

    !> Contains the points included in the reciprocal sum.  The set should not include the origin or
    !> inversion related points.
    real(dp), intent(in) :: recPoint(:,:)

    !> Parameter for Ewald summation.
    real(dp), intent(in) :: alpha

    !> Volume of the real space unit cell.
    real(dp), intent(in) :: volume

    !> Descriptor of the distributed matrix
    integer, intent(in) :: descInvRMat(DLEN_)

    !> Matrix of 1/R values for each atom pair.
    real(dp), intent(out) :: invRMat(:,:)

    type(TNeighIterator) :: neighIter
    type(TDynNeighList), pointer :: pNeighList
    real(dp) :: neighCoords(3, iterChunkSize)
    integer :: neighImages(iterChunkSize)
    integer :: iAt1, iAt2, iAt2f, iNeigh, nNeigh, jj, ii, iLoc, jLoc
    logical :: tLocal

    @:ASSERT(volume > 0.0_dp)

    invRMat(:,:) = 0.0_dp
    pNeighList => neighList

    ! Real space part of the Ewald sum.
    do jj = 1, size(invRMat, dim=2)
      iAt1 = scalafx_indxl2g(jj, descInvRMat(NB_), grid%mycol, descInvRMat(CSRC_), grid%ncol)
      call TNeighIterator_init(neighIter, pNeighList, iAt1)
      nNeigh = iterChunkSize
      do while (nNeigh == iterChunkSize)
        call neighIter%getNextNeighbours(nNeigh, coords=neighCoords, img2CentCell=neighImages)
        do iNeigh = 1, nNeigh
          iAt2f = neighImages(iNeigh)
          call scalafx_islocal(grid, descInvRMat, iAt2f, iAt1, tLocal, iLoc, jLoc)
          if (tLocal) then
            invRMat(iLoc, jLoc) = invRMat(iLoc, jLoc)&
                &  + rTerm(sqrt(sum((coord(:,iAt1) - neighCoords(:,iNeigh))**2)), alpha)
          end if
        end do
      end do
    end do

    ! Reciprocal space part of the Ewald sum.
    !$OMP PARALLEL DO&
    !$OMP& DEFAULT(SHARED) PRIVATE(iAt1, ii, iAt2) SCHEDULE(RUNTIME)
    do jj = 1, size(invRMat, dim=2)
      iAt1 = scalafx_indxl2g(jj, descInvRMat(NB_), grid%mycol, descInvRMat(CSRC_), grid%ncol)
      do ii = 1, size(invRMat, dim=1)
        iAt2 = scalafx_indxl2g(ii, descInvRMat(MB_), grid%myrow, descInvRMat(RSRC_), grid%nrow)
        if (iAt2 < iAt1) then
          cycle
        end if
        invRMat(ii, jj) = invRMat(ii, jj)&
            & + ewaldReciprocal(coord(:,iAt1) - coord(:,iAt2), recPoint, alpha, volume)&
            & - pi / (volume * alpha**2)
      end do
    end do
    !$OMP END PARALLEL DO

    ! Extra contribution for self interaction.
    !$OMP PARALLEL DO&
    !$OMP& DEFAULT(SHARED) PRIVATE(iAt1, tLocal, iLoc, jLoc) SCHEDULE(RUNTIME)
    do jj = 1, size(invRMat, dim=2)
      iAt1 = scalafx_indxl2g(jj, descInvRMat(NB_), grid%mycol, descInvRMat(CSRC_), grid%ncol)
      call scalafx_islocal(grid, descInvRMat, iAt1, iAt1, tLocal, iLoc, jLoc)
      if (tLocal) then
        invRMat(iLoc, jLoc) = invRMat(iLoc, jLoc) - 2.0_dp * alpha / sqrt(pi)
      end if
    end do
    !$OMP END PARALLEL DO

  end subroutine invRPeriodicBlacs


#:else


  !> Calculates the 1/R Matrix for all atoms for the periodic case (serial version)
  !>
  !> Note: Only the lower triangle is constructed.
  !>
  subroutine invRPeriodicSerial(coord, neighList, recPoint, alpha, volume, invRMat)

    !> List of atomic coordinates (all atoms).
    real(dp), intent(in) :: coord(:,:)

    !> Neighbour list for the real space summation of the Ewald.
    type(TDynNeighList), target, intent(in) :: neighList

    !> Contains the points included in the reciprocal sum.  The set should not include the origin or
    !> inversion related points.
    real(dp), intent(in) :: recPoint(:,:)

    !> Parameter for Ewald summation.
    real(dp), intent(in) :: alpha

    !> Volume of the real space unit cell.
    real(dp), intent(in) :: volume

    !> Matrix of 1/R values for each atom pair.
    real(dp), intent(out) :: invRMat(:,:)

    type(TDynNeighList), pointer :: pNeighList
    integer :: nAtom
    integer :: iAt1, iAt2

    nAtom = size(invRMat, dim=1)

    invRMat(:,:) = 0.0_dp
    pNeighList => neighList

    ! Real space part of the Ewald sum.
    !$OMP PARALLEL DO&
    !$OMP& DEFAULT(SHARED) SCHEDULE(RUNTIME)
    do iAt1 = 1, nAtom
      call addNeighbourContribsIS(iAt1, pNeighList, coord, alpha, invRMat)
    end do
    !$OMP END PARALLEL DO

    ! Reciprocal space part of the Ewald sum.
    !$OMP PARALLEL DO&
    !$OMP& DEFAULT(SHARED) PRIVATE(iAt2) SCHEDULE(RUNTIME)
    do iAt1 = 1, nAtom
      do iAt2 = iAt1, nAtom
        invRMat(iAt2, iAt1) = invRMat(iAt2, iAt1)&
            & + ewaldReciprocal(coord(:,iAt1)-coord(:,iAt2), recPoint, alpha, volume)&
            & - pi / (volume * alpha**2)
      end do
    end do
    !$OMP END PARALLEL DO

    ! Extra contribution for self interaction.
    !$OMP PARALLEL DO&
    !$OMP& DEFAULT(SHARED) SCHEDULE(RUNTIME)
    do iAt1 = 1, nAtom
      invRMat(iAt1, iAt1) = invRMat(iAt1, iAt1) - 2.0_dp * alpha / sqrt(pi)
    end do
    !$OMP END PARALLEL DO

  end subroutine invRPeriodicSerial

  !> Neighbour summation with local scope for predictable OMP <= 4.0 behaviour
  subroutine addNeighbourContribsIS(iAt1, pNeighList, coords, alpha, invRMat)
    integer, intent(in) :: iAt1
    type(TDynNeighList), pointer, intent(in) :: pNeighList
    real(dp), intent(in) :: coords(:,:)
    real(dp), intent(in) :: alpha
    real(dp), intent(inout) :: invRMat(:,:)

    type(TNeighIterator) :: neighIter
    real(dp) :: neighCoords(3, iterChunkSize)
    integer :: neighImages(iterChunkSize)
    integer :: iAt2f, iNeigh, nNeigh

    call TNeighIterator_init(neighIter, pNeighList, iAt1)
    nNeigh = iterChunkSize
    do while (nNeigh == iterChunkSize)
      call neighIter%getNextNeighbours(nNeigh, coords=neighCoords, img2CentCell=neighImages)
      do iNeigh = 1, nNeigh
        iAt2f = neighImages(iNeigh)
        invRMat(iAt2f, iAt1) = invRMat(iAt2f, iAt1)&
            & + rTerm(sqrt(sum((coords(:,iAt1) - neighCoords(:,iNeigh))**2)), alpha)
      end do
    end do

  end subroutine addNeighbourContribsIS

#:endif


  !> Calculates summed 1/R vector for two groups of objects for the periodic case.
  subroutine sumInvRPeriodicAsymm(env, nAtom0, nAtom1, coord0, coord1, charges1, rLat, gLat, alpha,&
      & volume, invRVec, blurwidths1, epsSoften)

    !> Computational environment settings
    type(TEnvironment), intent(in) :: env

    !> Number of atoms in the first group
    integer, intent(in) :: nAtom0

    !> Number of atoms in the second group
    integer, intent(in) :: nAtom1

    !> Coordinates of the first group of objects (atoms)
    real(dp), intent(in) :: coord0(:,:)

    !> Coordinates of the 2nd group of objects (point charges)
    real(dp), intent(in) :: coord1(:,:)

    !> Charges of the 2nd group of objects
    real(dp), intent(in) :: charges1(:)

    !> Lattice vectors to be used for the real Ewald summation
    real(dp), intent(in) :: rLat(:,:)

    !> Lattice vectors to be used for the reciprocal Ewald summation.
    real(dp), intent(in) :: gLat(:,:)

    !> Parameter of the Ewald summation
    real(dp), intent(in) :: alpha

    !> Volume of the supercell.
    real(dp), intent(in) :: volume

    !> Vector of sum_i q_i/|R_atom - R_i] values for each atom
    real(dp), intent(out) :: invRVec(:)

    !> Gaussian blur width of the charges in the 2nd group
    real(dp), intent(in), optional :: blurWidths1(:)

    !> Short distance softening
    real(dp), intent(in), optional :: epsSoften

    integer :: iAt0, iAt1
    real(dp) :: rTmp, rr(3)
    integer :: iAtFirst0, iAtLast0, iAtFirst1, iAtLast1

    @:ASSERT(volume > 0.0_dp)

    call distributeRangeInChunks2(env, 1, nAtom0, 1, nAtom1, iAtFirst0, iAtLast0, iAtFirst1,&
        & iAtLast1)

    invRVec(:) = 0.0_dp

    if (present(blurWidths1)) then
      !$OMP PARALLEL DO&
      !$OMP& DEFAULT(SHARED) PRIVATE(iAt1, rr, rTmp) SCHEDULE(RUNTIME)
      do iAt0 = iAtFirst0, iAtLast0
        do iAt1 = iAtFirst1, iAtLast1
          rr = coord0(:,iAt0) - coord1(:,iAt1)
          rTmp = ewald(rr, rLat, gLat, alpha, volume, blurWidth=blurWidths1(iAt1),&
              & epsSoften=epsSoften)
          invRVec(iAt0) = invRVec(iAt0) + rTmp * charges1(iAt1)
        end do
      end do
      !$OMP END PARALLEL DO
    else
      !$OMP PARALLEL DO&
      !$OMP& DEFAULT(SHARED) PRIVATE(iAt1, rr, rTmp) SCHEDULE(RUNTIME)
      do iAt0 = iAtFirst0, iAtLast0
        do iAt1 = iAtFirst1, iAtLast1
          rr = coord0(:,iAt0) - coord1(:,iAt1)
          rTmp = ewald(rr, rLat, gLat, alpha, volume, epsSoften=epsSoften)
          invRVec(iAt0) = invRVec(iAt0) + rTmp * charges1(iAt1)
        end do
      end do
      !$OMP END PARALLEL DO
    end if

    call assembleChunks(env, invRVec)

  end subroutine sumInvRPeriodicAsymm


  !> Calculates the -1/R**2 deriv contribution for all atoms for the non-periodic case, without
  !> storing anything.
  subroutine addInvRPrimeCluster(env, nAtom, coord, deltaQAtom, deriv)

    !> Computational environment settings
    type(TEnvironment), intent(in) :: env

    !> Number of atoms.
    integer, intent(in) :: nAtom

    !> List of atomic coordinates.
    real(dp), intent(in) :: coord(:,:)

    !> List of charges on each atom.
    real(dp), intent(in) :: deltaQAtom(:)

    !> Contains the derivative on exit.
    real(dp), intent(inout) :: deriv(:,:)

    integer :: ii, jj
    real(dp) :: dist, vect(3), fTmp
    real(dp), allocatable :: localDeriv(:,:)
    integer :: iAtFirst, iAtLast

    call distributeRangeInChunks(env, 1, nAtom, iAtFirst, iAtLast)

    allocate(localDeriv(3, nAtom))
    localDeriv(:,:) = 0.0_dp

    !$OMP PARALLEL DO&
    !$OMP& DEFAULT(SHARED) PRIVATE(jj, vect, dist, ftmp) REDUCTION(+:localDeriv) SCHEDULE(RUNTIME)
    do ii = iAtFirst, iAtLast
      do jj = ii + 1, nAtom
        vect(:) = coord(:,ii) - coord(:,jj)
        dist = sqrt(sum(vect(:)**2))
        fTmp = -deltaQAtom(ii) * deltaQAtom(jj) / (dist**3)
        localDeriv(:,ii) = localderiv(:,ii) + vect(:)*fTmp
        ! Skew-symmetric 1/r2 interaction, so the other triangle is calculated :
        localDeriv(:,jj) = localderiv(:,jj) - vect(:)*fTmp
      end do
    end do
    !$OMP END PARALLEL DO

    call assembleChunks(env, localDeriv)

    deriv(:,:) = deriv + localDeriv

  end subroutine addInvRPrimeCluster


  !> Calculates the -1/R**2 deriv contribution for extended lagrangian dynamics forces in a periodic
  !> geometry
  subroutine addInvRPrimeXlbomdCluster(env, nAtom, coord, dQInAtom, dQOutAtom, deriv)

    !> Computational environment settings
    type(TEnvironment), intent(in) :: env

    !> number of atoms
    integer, intent(in) :: nAtom

    !> coordinates of atoms
    real(dp), intent(in) :: coord(:,:)

    !> input charge fluctuations
    real(dp), intent(in) :: dQInAtom(:)

    !> output charge fluctuations
    real(dp), intent(in) :: dQOutAtom(:)

    !> energy derivative to add contribution to
    real(dp), intent(inout) :: deriv(:,:)

    integer :: iAt1, iAt2
    real(dp) :: dist, vect(3), fTmp, prefac
    real(dp), allocatable :: localDeriv(:,:)
    integer :: iAtFirst, iAtLast

    allocate(localDeriv(3, nAtom))
    localDeriv(:,:) = 0.0_dp

    call distributeRangeInChunks(env, 1, nAtom, iAtFirst, iAtLast)

    !$OMP PARALLEL DO&
    !$OMP& DEFAULT(SHARED) PRIVATE(iAt2, vect, dist, prefac, ftmp) REDUCTION(+:localDeriv)&
    !$OMP& SCHEDULE(RUNTIME)
    do iAt1 = iAtFirst, iAtLast
      do iAt2 = iAt1 + 1, nAtom
        vect(:) = coord(:,iAt1) - coord(:,iAt2)
        dist = sqrt(sum(vect(:)**2))
        prefac = dQOutAtom(iAt1) * dQInAtom(iAt2) + dQInAtom(iAt1) * dQOutAtom(iAt2) &
            & - dQInAtom(iAt1) * dQInAtom(iAt2)
        fTmp = -prefac / (dist**3)
        localDeriv(:,iAt1) = localDeriv(:,iAt1) + vect * fTmp
        ! Skew-symmetric 1/r2 interaction, so the other triangle is calculated
        localDeriv(:,iAt2) = localDeriv(:,iAt2) - vect *fTmp
      end do
    end do
    !$OMP END PARALLEL DO

    call assembleChunks(env, localDeriv)

    deriv(:,:) = deriv + localDeriv

  end subroutine addInvRPrimeXlbomdCluster


  !> Calculates the -1/R**2 deriv contribution for charged atoms interacting with a group of charged
  !> objects (like point charges) for the non-periodic case, without storing anything.
  subroutine addInvRPrimeClusterAsymm(env, nAtom0, nAtom1, coord0, coord1, charge0, charge1,&
      & deriv0, deriv1, tHamDeriv, blurWidths1)

    !> Computational environment settings
    type(TEnvironment), intent(in) :: env

    !> Number of atoms in the first group
    integer, intent(in) :: nAtom0

    !> Number of atoms in the second group
    integer, intent(in) :: nAtom1

    !> List of atomic coordinates.
    real(dp), intent(in) :: coord0(:,:)

    !> List of the point charge coordinates
    real(dp), intent(in) :: coord1(:,:)

    !> Charge of the atoms.
    real(dp), intent(in) :: charge0(:)

    !> Charge of the point charges.
    real(dp), intent(in) :: charge1(:)

    !> Contains the derivative for the first group
    real(dp), intent(inout) :: deriv0(:,:)

    !> Contains the derivative for the second group
    real(dp), intent(inout) :: deriv1(:,:)

    !> Compute the derivative of Hamiltonians? Otherwise, compute the force
    logical, intent(in) :: tHamDeriv

    !> if gaussian distribution for the charge
    real(dp), intent(in), optional :: blurWidths1(:)

    integer :: iAt0, iAt1
    real(dp) :: dist, vect(3), fTmp(3), sigma, rs
    integer :: iAtFirst0, iAtLast0, iAtFirst1, iAtLast1
    real(dp), allocatable :: localDeriv0(:,:), localDeriv1(:,:)

    allocate(localDeriv0(3, nAtom0))
    localDeriv0(:,:) = 0.0_dp

    if (.not. tHamDeriv) then
      allocate(localDeriv1(3, nAtom1))
      localDeriv1(:,:) = 0.0_dp
    end if

    call distributeRangeInChunks2(env, 1, nAtom0, 1, nAtom1, iAtFirst0, iAtLast0, iAtFirst1,&
        & iAtLast1)

    ! Doing blured and unblured cases separately to avoid ifs in the loop
    if (present(blurWidths1)) then
      if (tHamDeriv) then
        !$OMP PARALLEL DO&
        !$OMP& DEFAULT(SHARED) PRIVATE(iAt1, vect, dist, ftmp, sigma, rs)&
        !$OMP& REDUCTION(+:localDeriv0) SCHEDULE(RUNTIME)
        do iAt0 = iAtFirst0, iAtLast0
          do iAt1 = iAtFirst1, iAtLast1
            vect(:) = coord0(:,iAt0) - coord1(:,iAt1)
            dist = sqrt(sum(vect(:)**2))
            fTmp = -vect(:) / (dist**3)
            if (dist < erfArgLimit * blurWidths1(iAt1)) then
              sigma = blurWidths1(iAt1)
              rs = dist / sigma
              fTmp = fTmp * (erfwrap(rs) - 2.0_dp/(sqrt(pi)*sigma) * dist * exp(-(rs**2)))
            end if
            fTmp = charge1(iAt1) * fTmp
            localDeriv0(:,iAt0) = localDeriv0(:,iAt0) + fTmp(:)
          end do
        end do
        !$OMP END PARALLEL DO
      else
        !$OMP PARALLEL DO&
        !$OMP& DEFAULT(SHARED) PRIVATE(iAt1, vect, dist, ftmp, sigma, rs)&
        !$OMP& REDUCTION(+:localDeriv0, localDeriv1) SCHEDULE(RUNTIME)
        do iAt0 = iAtFirst0, iAtLast0
          do iAt1 = iAtFirst1, iAtLast1
            vect(:) = coord0(:,iAt0) - coord1(:,iAt1)
            dist = sqrt(sum(vect(:)**2))
            fTmp = -vect(:) / (dist**3)
            if (dist < erfArgLimit * blurWidths1(iAt1)) then
              sigma = blurWidths1(iAt1)
              rs = dist / sigma
              fTmp = fTmp * (erfwrap(rs) - 2.0_dp/(sqrt(pi)*sigma) * dist * exp(-(rs**2)))
            end if
            fTmp = charge0(iAt0) * charge1(iAt1) * fTmp
            localDeriv0(:,iAt0) = localDeriv0(:,iAt0) + fTmp(:)
            localDeriv1(:,iAt1) = localDeriv1(:,iAt1) - fTmp(:)
          end do
        end do
        !$OMP END PARALLEL DO
      end if
    else
      if (tHamDeriv) then
        !$OMP PARALLEL DO&
        !$OMP& DEFAULT(SHARED) PRIVATE(iAt1, vect, dist, ftmp)&
        !$OMP& REDUCTION(+:localDeriv0) SCHEDULE(RUNTIME)
        do iAt0 = iAtFirst0, iAtLast0
          do iAt1 = iAtFirst1, iAtLast1
            vect(:) = coord0(:,iAt0) - coord1(:,iAt1)
            dist = sqrt(sum(vect(:)**2))
            fTmp = -charge1(iAt1) / (dist**3) * vect(:)
            localDeriv0(:,iAt0) = localDeriv0(:,iAt0) + fTmp(:)
          end do
        end do
        !$OMP END PARALLEL DO
      else
        !$OMP PARALLEL DO&
        !$OMP& DEFAULT(SHARED) PRIVATE(iAt1, vect, dist, ftmp)&
        !$OMP& REDUCTION(+:localDeriv0, localDeriv1) SCHEDULE(RUNTIME)
        do iAt0 = iAtFirst0, iAtLast0
          do iAt1 = iAtFirst1, iAtLast1
            vect(:) = coord0(:,iAt0) - coord1(:,iAt1)
            dist = sqrt(sum(vect(:)**2))
            fTmp = -charge0(iAt0) * charge1(iAt1) / (dist**3) * vect(:)
            localDeriv0(:,iAt0) = localDeriv0(:,iAt0) + fTmp(:)
            localDeriv1(:,iAt1) = localDeriv1(:,iAt1) - fTmp(:)
          end do
        end do
        !$OMP END PARALLEL DO
      end if
    end if

    call assembleChunks(env, localDeriv0)
    deriv0(:,:) = deriv0 + localDeriv0

    if (.not. tHamDeriv) then
      call assembleChunks(env, localDeriv1)
      deriv1(:,:) = deriv1 + localDeriv1
    end if

  end subroutine addInvRPrimeClusterAsymm


  !> Calculates the -1/R**2 deriv contribution for the periodic case, without storing anything.
  subroutine addInvRPrimePeriodic(env, nAtom, coord, neighList, recPoint, alpha, volume,&
      & deltaQAtom, deriv)

    !> Computational environment settings
    type(TEnvironment), intent(in) :: env

    !> Number of atoms
    integer, intent(in) :: nAtom

    !> List of atomic coordinates (all atoms).
    real(dp), intent(in) :: coord(:,:)

    !> Dynamic neighbour list to be used in the real part of Ewald
    type(TDynNeighList), target, intent(in) :: neighList

    !> Contains the points included in the reciprocal sum. The set should not include the origin or
    !> inversion related points.
    real(dp), intent(in) :: recPoint(:,:)

    !> Parameter for Ewald summation.
    real(dp), intent(in) :: alpha

    !> Volume of the real space unit cell.
    real(dp), intent(in) :: volume

    !> List of charges on each atom
    real(dp), intent(in) :: deltaQAtom(:)

    !> Derivative on exit
    real(dp), intent(inout) :: deriv(:,:)

    type(TDynNeighList), pointer :: pNeighList
    integer :: iAtom1, iAtom2
    real(dp) :: r(3)
    real(dp), allocatable :: localDeriv(:,:)
    integer :: iAtFirst, iAtLast

    pNeighList => neighList
    allocate(localDeriv(3, nAtom))
    localDeriv(:,:) = 0.0_dp

    call distributeRangeInChunks(env, 1, nAtom, iAtFirst, iAtLast)

    ! d(1/R)/dr real space
    !$OMP PARALLEL DO&
    !$OMP& DEFAULT(SHARED) REDUCTION(+:localDeriv) SCHEDULE(RUNTIME)
    do iAtom1 = iAtFirst, iAtLast
      call addNeighbourContribsInvRP(iAtom1, pNeighList, coord, deltaQAtom, alpha, localDeriv)
    end do
    !$OMP END PARALLEL DO

    ! d(1/R)/dr reciprocal space
    !$OMP PARALLEL DO&
    !$OMP& DEFAULT(SHARED) PRIVATE(iAtom2, r) REDUCTION(+:localDeriv) SCHEDULE(RUNTIME)
    do iAtom1 = iAtFirst, iAtLast
      do iAtom2 = iAtom1+1, nAtom
        r(:) = coord(:,iAtom1)-coord(:,iAtom2)
        localDeriv(:,iAtom1) = localDeriv(:,iAtom1)&
            & + derivEwaldReciprocal(r,recPoint,alpha,volume)*deltaQAtom(iAtom1)*deltaQAtom(iAtom2)
        localDeriv(:,iAtom2) = localDeriv(:,iAtom2)&
            & - derivEwaldReciprocal(r,recPoint,alpha,volume)*deltaQAtom(iAtom1)*deltaQAtom(iAtom2)
      end do
    end do
    !$OMP END PARALLEL DO

    call assembleChunks(env, localDeriv)

    deriv(:,:) = deriv + localDeriv

  end subroutine addInvRPrimePeriodic

  !> Neighbour summation with local scope for predictable OMP <= 4.0 behaviour
  subroutine addNeighbourContribsInvRP(iAtom1, pNeighList, coords, deltaQAtom, alpha, deriv)
    integer, intent(in) :: iAtom1
    type(TDynNeighList), pointer, intent(in) :: pNeighList
    real(dp), intent(in) :: coords(:,:)
    real(dp), intent(in) :: deltaQAtom(:)
    real(dp), intent(in) :: alpha
    real(dp), intent(inout) :: deriv(:,:)

    type(TNeighIterator) :: neighIter
    real(dp) :: neighCoords(3, iterChunkSize)
    integer :: neighImages(iterChunkSize)
    integer :: iAtom2f, iNeigh, nNeigh
    real(dp) :: rr(3)

    call TNeighIterator_init(neighIter, pNeighList, iAtom1)
    nNeigh = iterChunkSize
    do while (nNeigh == iterChunkSize)
      call neighIter%getNextNeighbours(nNeigh, coords=neighCoords, img2CentCell=neighImages)
      do iNeigh = 1, nNeigh
        iAtom2f = neighImages(iNeigh)
        if (iAtom2f /= iAtom1) then
          rr(:) = coords(:,iAtom1) - neighCoords(:,iNeigh)
          deriv(:,iAtom1) = deriv(:,iAtom1)&
              & + derivRTerm(rr, alpha) * deltaQAtom(iAtom1) * deltaQAtom(iAtom2f)
          deriv(:,iAtom2f) = deriv(:,iAtom2f)&
              & - derivRTerm(rr, alpha) * deltaQAtom(iAtom1) * deltaQAtom(iAtom2f)
        end if
      end do
    end do

  end subroutine addNeighbourContribsInvRP


  !> Calculates the -1/R**2 deriv contribution for extended lagrangian dynamics forces
  subroutine addInvRPrimeXlbomdPeriodic(env, nAtom, coord, neighList, recPoint, alpha, volume,&
      & dQInAtom, dQOutAtom, deriv)

    !> Computational environment settings
    type(TEnvironment), intent(in) :: env

    !> number of atoms
    integer, intent(in) :: nAtom

    !> coordinates of atoms
    real(dp), intent(in) :: coord(:,:)

    !> Dynamic neighbour list to be used in the real part of Ewald
    type(TDynNeighList), target, intent(in) :: neighList

    !> Contains the points included in the reciprocal sum. The set should not include the origin or
    !> inversion related points.
    real(dp), intent(in) :: recPoint(:,:)

    !> Ewald parameter
    real(dp), intent(in) :: alpha

    !> cell volume
    real(dp), intent(in) :: volume

    !> input charge fluctuations
    real(dp), intent(in) :: dQInAtom(:)

    !> output charge fluctuations
    real(dp), intent(in) :: dQOutAtom(:)

    !> energy derivative to add contribution to
    real(dp), intent(inout) :: deriv(:,:)

    type(TDynNeighList), pointer :: pNeighList
    integer :: iAt1, iAt2
    real(dp) :: rr(3), contrib(3), prefac
    real(dp), allocatable :: localDeriv(:,:)
    integer :: iAtFirst, iAtLast

    pNeighList => neighList
    allocate(localDeriv(3, nAtom))
    localDeriv(:,:) = 0.0_dp

    call distributeRangeInChunks(env, 1, nAtom, iAtFirst, iAtLast)

    ! real space
    !$OMP PARALLEL DO&
    !$OMP& DEFAULT(SHARED) REDUCTION(+:localDeriv) SCHEDULE(RUNTIME)
    do iAt1 = iAtFirst, iAtLast
      call addNeighbourContribsXl(iAt1, pNeighList, coord, dQInAtom, dQOutAtom, alpha, localDeriv)
    end do
    !$OMP END PARALLEL DO

    ! reciprocal space
    !$OMP PARALLEL DO&
    !$OMP& DEFAULT(SHARED) PRIVATE(iAt2, rr, prefac, contrib)  REDUCTION(+:localDeriv)
    do iAt1 = iAtFirst, iAtLast
      do iAt2 = iAt1 + 1, nAtom
        rr(:) = coord(:,iAt1) - coord(:,iAt2)
        prefac = dQOutAtom(iAt1) * dQInAtom(iAt2) + dQInAtom(iAt1) * dQOutAtom(iAt2)&
            & - dQInAtom(iAt1) * dQInAtom(iAt2)
        contrib(:) = prefac * derivEwaldReciprocal(rr, recPoint, alpha, volume)
        localDeriv(:,iAt1) = localDeriv(:,iAt1) + contrib
        localDeriv(:,iAt2) = localDeriv(:,iAt2) - contrib
      end do
    end do
    !$OMP END PARALLEL DO

    call assembleChunks(env, localDeriv)
    deriv(:,:) = deriv + localDeriv

  end subroutine addInvRPrimeXlbomdPeriodic

  !> Neighbour summation with local scope for predictable OMP <= 4.0 behaviour
  subroutine addNeighbourContribsXl(iAt1, pNeighList, coords, dQInAtom, dQOutAtom, alpha, deriv)
    integer, intent(in) :: iAt1
    type(TDynNeighList), pointer, intent(in) :: pNeighList
    real(dp), intent(in) :: coords(:,:)
    real(dp), intent(in) :: dQInAtom(:)
    real(dp), intent(in) :: dQOutAtom(:)
    real(dp), intent(in) :: alpha
    real(dp), intent(inout) :: deriv(:,:)

    type(TNeighIterator) :: neighIter
    real(dp) :: neighCoords(3, iterChunkSize)
    real(dp) :: rr(3), contrib(3)
    real(dp) :: prefac
    integer :: neighImages(iterChunkSize)
    integer :: iAt2f, iNeigh, nNeigh

    call TNeighIterator_init(neighIter, pNeighList, iAt1)
    nNeigh = iterChunkSize
    do while (nNeigh == iterChunkSize)
      call neighIter%getNextNeighbours(nNeigh, coords=neighCoords, img2CentCell=neighImages)
      do iNeigh = 1, nNeigh
        iAt2f = neighImages(iNeigh)
        if (iAt2f == iAt1) then
          cycle
        end if
        rr(:) = coords(:,iAt1) - neighCoords(:,iNeigh)
        prefac = dQOutAtom(iAt1) * dQInAtom(iAt2f) + dQInAtom(iAt1) * dQOutAtom(iAt2f)&
            & - dQInAtom(iAt1) * dQInAtom(iAt2f)
        contrib(:) = prefac * derivRTerm(rr, alpha)
        deriv(:,iAt1) = deriv(:,iAt1) + contrib
        deriv(:,iAt2f) = deriv(:,iAt2f) - contrib
      end do
    end do

  end subroutine addNeighbourContribsXl


  !> Calculates the -1/R**2 deriv contribution for charged atoms interacting with a group of charged
  !> objects (like point charges) for the periodic case, without storing anything.
  subroutine addInvRPrimePeriodicAsymm(env, nAtom0, nAtom1, coord0, coord1, charge0, charge1, rVec,&
      & gVec, alpha, vol, deriv0, deriv1, tHamDeriv, blurWidths1)

    !> Computational environment settings
    type(TEnvironment), intent(in) :: env

    !> Number of atoms in the first group
    integer, intent(in) :: nAtom0

    !> Number of atoms in the second group
    integer, intent(in) :: nAtom1

    !> List of atomic coordinates (first group)
    real(dp), intent(in) :: coord0(:,:)

    !> List of the point charge coordinates (second group)
    real(dp), intent(in) :: coord1(:,:)

    !> Charge of the atoms in group 1.
    real(dp), intent(in) :: charge0(:)

    !> Charge of the point charges.
    real(dp), intent(in) :: charge1(:)

    !> Lattice vectors to be used for the real Ewald summation
    real(dp), intent(in) :: rVec(:,:)

    !> Lattice vectors to be used for the reciprocal Ewald summation.
    real(dp), intent(in) :: gVec(:,:)

    !> Parameter of the Ewald summation
    real(dp), intent(in) :: alpha

    !> Volume of the supercell.
    real(dp), intent(in) :: vol

    !> Contains the derivative for the first group on exit
    real(dp), intent(inout) :: deriv0(:,:)

    !> Contains the derivative for the second group on exit
    real(dp), intent(inout) :: deriv1(:,:)

    !> Compute the derivative of Hamiltonians? Otherwise, compute the force
    logical, intent(in) :: tHamDeriv

    !> Gaussian blur width of the charges in the 2nd group
    real(dp), intent(in), optional :: blurWidths1(:)

    integer :: iAt0, iAt1
    real(dp) :: vect(3), fTmp(3)
    integer :: iAtFirst0, iAtLast0, iAtFirst1, iAtLast1
    real(dp), allocatable :: localDeriv0(:,:), localDeriv1(:,:)

    @:ASSERT(vol > 0.0_dp)

    allocate(localDeriv0(3, nAtom0))
    localDeriv0(:,:) = 0.0_dp

    if (.not. tHamDeriv) then
      allocate(localDeriv1(3, nAtom1))
      localDeriv1(:,:) = 0.0_dp
    end if

    call distributeRangeInChunks2(env, 1, nAtom0, 1, nAtom1, iAtFirst0, iAtLast0, iAtFirst1,&
        & iAtLast1)

    ! real space part
    if (present(blurwidths1)) then
      if (tHamDeriv) then
        !$OMP PARALLEL DO&
        !$OMP& DEFAULT(SHARED) PRIVATE(iAt1, vect, fTmp) REDUCTION(+:localDeriv0)&
        !$OMP& SCHEDULE(RUNTIME)
        do iAt0 = iAtFirst0, iAtLast0
          do iAt1 = iAtFirst1, iAtLast1
            vect(:) = coord0(:,iAt0) - coord1(:,iAt1)
            fTmp(:) = derivEwaldReal(vect, rVec, alpha, blurWidth=blurWidths1(iAt1))&
                & * charge1(iAt1)
            localDeriv0(:,iAt0) = localDeriv0(:,iAt0) + fTmp(:)
          end do
        end do
        !$OMP END PARALLEL DO
      else
        !$OMP PARALLEL DO&
        !$OMP& DEFAULT(SHARED) PRIVATE(iAt1, vect, fTmp) REDUCTION(+:localDeriv0, localDeriv1)&
        !$OMP& SCHEDULE(RUNTIME)
        do iAt0 = iAtFirst0, iAtLast0
          do iAt1 = iAtFirst1, iAtLast1
            vect(:) = coord0(:,iAt0) - coord1(:,iAt1)
            fTmp(:) = derivEwaldReal(vect, rVec, alpha, blurWidth=blurWidths1(iAt1))&
                & * charge0(iAt0) * charge1(iAt1)
            localDeriv0(:,iAt0) = localDeriv0(:,iAt0) + fTmp(:)
            localDeriv1(:,iAt1) = localDeriv1(:,iAt1) - fTmp(:)
          end do
        end do
        !$OMP END PARALLEL DO
      end if
    else
      if (tHamDeriv) then
        !$OMP PARALLEL DO&
        !$OMP& DEFAULT(SHARED) PRIVATE(iAt1, vect, fTmp) REDUCTION(+:localDeriv0)&
        !$OMP& SCHEDULE(RUNTIME)
        do iAt0 = iAtFirst0, iAtLast0
          do iAt1 = iAtFirst1, iAtLast1
            vect(:) = coord0(:,iAt0) - coord1(:,iAt1)
            fTmp(:) = derivEwaldReal(vect, rVec, alpha) * charge1(iAt1)
            localDeriv0(:,iAt0) = localDeriv0(:,iAt0) + fTmp(:)
          end do
        end do
        !$OMP END PARALLEL DO
      else
        !$OMP PARALLEL DO&
        !$OMP& DEFAULT(SHARED) PRIVATE(iAt1, vect, fTmp) REDUCTION(+:localDeriv0, localDeriv1)&
        !$OMP& SCHEDULE(RUNTIME)
        do iAt0 = iAtFirst0, iAtLast0
          do iAt1 = iAtFirst1, iAtLast1
            vect(:) = coord0(:,iAt0) - coord1(:,iAt1)
            fTmp(:) = derivEwaldReal(vect, rVec, alpha) * charge0(iAt0) * charge1(iAt1)
            localDeriv0(:,iAt0) = localDeriv0(:,iAt0) + fTmp(:)
            localDeriv1(:,iAt1) = localDeriv1(:,iAt1) - fTmp(:)
          end do
        end do
        !$OMP END PARALLEL DO
      end if
    end if

    if (tHamDeriv) then
      ! reciprocal space part
      !$OMP PARALLEL DO&
      !$OMP& DEFAULT(SHARED) PRIVATE(iAt1, vect, fTmp) REDUCTION(+:localDeriv0)&
      !$OMP& SCHEDULE(RUNTIME)
      do iAt0 = iAtFirst0, iAtLast0
        do iAt1 = iAtFirst1, iAtLast1
          vect(:) = coord0(:,iAt0) - coord1(:,iAt1)
          fTmp(:) = derivEwaldReciprocal(vect, gVec, alpha, vol) * charge1(iAt1)
          localDeriv0(:,iAt0) = localDeriv0(:,iAt0) + fTmp(:)
        end do
      end do
      !$OMP END PARALLEL DO
    else
      ! reciprocal space part
      !$OMP PARALLEL DO&
      !$OMP& DEFAULT(SHARED) PRIVATE(iAt1, vect, fTmp) REDUCTION(+:localDeriv0, localDeriv1)&
      !$OMP& SCHEDULE(RUNTIME)
      do iAt0 = iAtFirst0, iAtLast0
        do iAt1 = iAtFirst1, iAtLast1
          vect(:) = coord0(:,iAt0) - coord1(:,iAt1)
          fTmp(:) = derivEwaldReciprocal(vect, gVec, alpha, vol) * charge0(iAt0) * charge1(iAt1)
          localDeriv0(:,iAt0) = localDeriv0(:,iAt0) + fTmp(:)
          localDeriv1(:,iAt1) = localDeriv1(:,iAt1) - fTmp(:)
        end do
      end do
      !$OMP END PARALLEL DO
    end if

    call assembleChunks(env, localDeriv0)
    deriv0(:,:) = deriv0 + localDeriv0

    if (.not. tHamDeriv) then
      call assembleChunks(env, localDeriv1)
      deriv1(:,:) = deriv1 + localDeriv1
    end if

  end subroutine addInvRPrimePeriodicAsymm


  !> Get optimal alpha-parameter for the Ewald summation by finding alpha, where decline of real and
  !> reciprocal part of Ewald are equal.
  !> The function stops, if the optimal alpha cannot be found.
  function getOptimalAlphaEwald(latVec, recVec, volume, tolerance) result(alpha)

    !> Lattice vectors.
    real(dp), intent(in) :: latVec(:,:)

    !> Reciprocal vectors.
    real(dp), intent(in) :: recVec(:,:)

    !> Volume of the unit cell.
    real(dp), intent(in) :: volume

    !> Tolerance for difference in real and rec. part.
    real(dp), intent(in) :: tolerance

    !> Optimal alpha.
    real(dp) :: alpha

    real(dp) :: alphaLeft, alphaRight
    real(dp), parameter :: alphaInit = 1.0e-8_dp

    real(dp) :: minG, minR, diff
    integer :: iIter
    integer :: iError
    character(len=100) :: errorString

    @:ASSERT(volume > 0.0_dp)
    @:ASSERT(tolerance > 0.0_dp)

    minG = sqrt(minval(sum(recVec(:,:)**2, dim=1)))
    minR = sqrt(minval(sum(latVec(:,:)**2, dim=1)))

    iError = 0
    alpha = alphaInit
    diff = diffRecReal(alpha, minG, minR, volume)
    do while (diff < -tolerance .and. alpha <= huge(1.0_dp))
      alpha = 2.0_dp * alpha
      diff = diffRecReal(alpha, minG, minR, volume)
    end do
    if (alpha > huge(1.0_dp)) then
      iError = 1
    elseif (alpha == alphaInit) then
      iError = 2
    end if

    if (iError == 0) then
      alphaLeft = 0.5_dp * alpha
      do while (diff < tolerance .and. alpha <= huge(1.0_dp))
        alpha = 2.0_dp * alpha
        diff = diffRecReal(alpha, minG, minR, volume)
      end do
      if (alpha > huge(1.0_dp)) then
        iError = 3
      end if
    end if

    if (iError == 0) then
      alphaRight = alpha
      alpha = (alphaLeft + alphaRight) / 2.0
      iIter = 0
      diff = diffRecReal(alpha, minG, minR, volume)
      do while (abs(diff) > tolerance .and. iIter <= nSearchIter)
        if (diff < 0) then
          alphaLeft = alpha
        else
          alphaRight = alpha
        end if
        alpha = (alphaLeft + alphaRight) / 2.0
        diff = diffRecReal(alpha, minG, minR, volume)
        iIter = iIter + 1
      end do
      if (iIter > nSearchIter) then
        iError = 4
      end if
    end if

    if (iError /= 0) then
      !alpha = exp(-0.310104 * log(volume) + 0.786382) / 2.0
99000 format ('Failure in determining optimal alpha for Ewaldsum.', ' Error code: ',I3)
      write(errorString, 99000) iError
      call error(errorString)
    end if

  end function getOptimalAlphaEwald


  !> Returns the longest reciprocal vector which gives a bigger contribution to the Ewald sum than a
  !> certain tolerance.
  function getMaxGEwald(alpha, volume, minValue) result(xx)

    !> Parameter of the ewald summation.
    real(dp), intent(in) :: alpha

    !> Volume of the unit cell.
    real(dp), intent(in) :: volume

    !> Tolerance value.
    real(dp), intent(in) :: minValue

    !> magnitude of reciprocal vector
    real(dp) :: xx

    real(dp), parameter :: gInit = 1.0e-8_dp
    real(dp) :: xLeft, xRight, yLeft, yRight, yy
    integer :: iError, iIter
    character(len=100) :: errorString

    iError = 0
    xx = gInit
    yy = gTerm(xx, alpha, volume)
    do while (yy > minValue .and. xx <= huge(1.0_dp))
      xx = 2.0_dp * xx
      yy = gTerm(xx, alpha, volume)
    end do
    if (xx > huge(1.0_dp)) then
      iError = 1
    elseif (xx == gInit) then
      iError = 2
    end if

    if (iError == 0) then
      xLeft = 0.5_dp * xx
      yLeft = gTerm(xLeft, alpha, volume)
      xRight = xx
      yRight = yy

      iIter = 1
      do while (yLeft - yRight > minValue .and. iIter <= nSearchIter)
        xx = 0.5_dp * (xLeft + xRight)
        yy = gTerm(xx, alpha, volume)
        if (yy >= minValue) then
          xLeft = xx
          yLeft = yy
        else
          xRight = xx
          yRight = yy
        end if
        iIter = iIter + 1
      end do
      if (iIter > nSearchIter) then
        iError = 3
      end if
    end if

    if (iError /= 0) then
99010 format ('Failure in getMaxGEwald.', ' Error nr: ',I3)
      write(errorString, 99010) iError
      call error(errorString)
    end if

  end function getMaxGEwald


  !> Returns the longest real space vector which gives a bigger contribution to the Ewald sum than a
  !> certain tolerance.
  function getMaxREwald(alpha, minValue) result(xx)

    !> Parameter of the ewald summation.
    real(dp), intent(in) :: alpha

    !> Tolerance value.
    real(dp), intent(in) :: minValue

    !> Magnitude of real space vector
    real(dp) :: xx

    real(dp), parameter :: rInit = 1.0e-8_dp
    real(dp) :: xLeft, xRight, yLeft, yRight, yy
    integer :: iError, iIter
    character(len=100) :: errorString

    iError = 0
    xx = rInit
    yy = rTerm(xx, alpha)
    do while (yy > minValue .and. xx <= huge(1.0_dp))
      xx = 2.0_dp * xx
      yy = rTerm(xx, alpha)
    end do
    if (xx > huge(1.0_dp)) then
      iError = 1
    elseif (xx == rInit) then
      iError = 2
    end if

    if (iError == 0) then
      xLeft = 0.5_dp * xx
      yLeft = rTerm(xLeft, alpha)
      xRight = xx
      yRight = yy

      iIter = 1
      do while (yLeft - yRight > minValue .and. iIter <= nSearchIter)
        xx = 0.5_dp * (xLeft + xRight)
        yy = rTerm(xx, alpha)
        if (yy >= minValue) then
          xLeft = xx
          yLeft = yy
        else
          xRight = xx
          yRight = yy
        end if
        iIter = iIter + 1
      end do
      if (iIter > nSearchIter) then
        iError = 3
      end if
    end if

    if (iError /= 0) then
99020 format ('Failure in getMaxREwald.', ' Error nr: ',I3)
      write(errorString, 99020) iError
      call error(errorString)
    end if

  end function getMaxREwald


  !> Returns the Ewald sum for a given lattice in a given point.
  function ewald(rr, rVec, gVec, alpha, vol, blurWidth, epsSoften)

    !> Vector where to calculate the Ewald sum.
    real(dp), intent(in) :: rr(:)

    !> Real space vectors to sum over. (Should contain origin).
    real(dp), intent(in) :: rVec(:,:)

    !> Reciprocal space vectors to sum over (Should not contain either origin nor inversion related
    !> points).
    real(dp), intent(in) :: gVec(:,:)

    !> Parameter for the Ewald summation.
    real(dp), intent(in) :: alpha

    !> Volume of the real space unit cell.
    real(dp), intent(in) :: vol

    !> Gaussian blur width of the charges in the 2nd group
    real(dp), intent(in), optional :: blurWidth

    !> Short distance softening
    real(dp), intent(in), optional :: epsSoften

    !> Result
    real(dp) :: ewald

    ewald = ewaldReciprocal(rr, gVec, alpha, vol)&
        & + ewaldReal(rr, rVec, alpha, blurWidth=blurWidth, epsSoften=epsSoften)&
        & - pi / (vol*alpha**2)
    if (sum(rr(:)**2) < tolSameDist2) then
      ewald = ewald - 2.0_dp * alpha / sqrt(pi)
    end if

  end function ewald


  !> Returns the reciprocal part of the Ewald sum.
  function ewaldReciprocal(rr, gVec, alpha, vol) result(recSum)

    !> Vector where to calculate the Ewald sum.
    real(dp), intent(in) :: rr(:)

    !> Reciprocal space vectors to sum over (Should not contain either origin nor inversion related
    !> points).
    real(dp), intent(in) :: gVec(:,:)

    !> Parameter for the Ewald summation.
    real(dp), intent(in) :: alpha

    !> Volume of the real space unit cell.
    real(dp), intent(in) :: vol

    !> contribution to the sum
    real(dp) :: recSum

    real(dp) :: gg(3), g2
    integer :: iG

    @:ASSERT(vol > 0.0_dp)

    recSum = 0.0_dp
    do iG = 1, size(gVec, dim=2)
      gg = gVec(:,iG)
      g2 = sum(gg(:)**2)
      recSum = recSum + exp(-g2/(4.0_dp*alpha**2))/g2 * cos(dot_product(gg,rr))
    end do
    ! note factor of 2 as only summing half of reciprocal space
    recSum = 2.0_dp * recSum * 4.0_dp * pi / vol

  end function ewaldReciprocal


  !> Returns the derivative of the reciprocal part of the Ewald sum.
  function derivEwaldReciprocal(rr, gVec, alpha, vol) result(recSum)

    !> Vector where to calculate the Ewald sum.
    real(dp), intent(in) :: rr(:)

    !> Reciprocal space vectors to sum over (Should not contain either origin nor inversion related
    !> points).
    real(dp), intent(in) :: gVec(:,:)

    !> Parameter for the Ewald summation.
    real(dp), intent(in) :: alpha

    !> Volume of the real space unit cell.
    real(dp), intent(in) :: vol

    !> contribution to the derivative value
    real(dp) :: recSum(3)

    real(dp) :: gg(3), g2
    integer :: iG

    @:ASSERT(vol > 0.0_dp)

    recSum(:) = 0.0_dp
    do iG = 1, size(gVec, dim=2)
      gg(:) = gVec(:,iG)
      g2 = sum(gg(:)**2)
      recSum(:) = recSum(:) - gg(:)*sin(dot_product(gg,rr))*exp(-g2/(4.0_dp*alpha**2))/g2
    end do
    ! note factor of 2 as only summing over half of reciprocal space
    recSum(:) = 2.0_dp * recSum(:) * 4.0_dp * pi / vol

  end function derivEwaldReciprocal


  !> Returns the derivative and stress of the reciprocal part of the Ewald sum
  subroutine derivStressEwaldRec(rr, gVec, alpha, vol, recSum, sigma)

    !> Vector where to calculate the Ewald sum.
    real(dp), intent(in) :: rr(:)

    !> Reciprocal space vectors to sum over.
    !  Should not contain either origin nor inversion related points.
    real(dp), intent(in) :: gVec(:, :)

    !> Parameter for the Ewald summation.
    real(dp), intent(in) :: alpha

    !> Volume of the real space unit cell.
    real(dp), intent(in) :: vol

    !> contribution to the derivative value
    real(dp), intent(out) :: recSum(3)

    !> contribution to the derivative value
    real(dp), intent(out) :: sigma(3, 3)

    real(dp), parameter :: unity(3, 3) = reshape(&
        & [1.0_dp, 0.0_dp, 0.0_dp, 0.0_dp, 1.0_dp, 0.0_dp, 0.0_dp, 0.0_dp, 1.0_dp], [3, 3])

    integer :: iG
    real(dp) :: gg(3), g2, rg, eTerm, sTmp

    recSum(:) = 0.0_dp
    sigma(:,:) = 0.0_dp
    do iG = 1, size(gVec, dim=2)
      gg = gVec(:, iG)
      g2 = sum(gg**2)
      rg = dot_product(gg, rr)
      eTerm = exp(-g2 / (4.0_dp * alpha**2)) / g2
      recSum(:) = recSum - gg * sin(rg) * eTerm
      sTmp = 2.0_dp * (1.0_dp / (4.0_dp * alpha * alpha) + 1.0_dp / g2)
      sigma(:,:) = sigma + (-unity + sTmp * spread(gg, 1, 3)*spread(gg, 2, 3)) * cos(rg) * eTerm
    end do
    ! note factor of 2 as only summing over half of reciprocal space
    recSum(:) = 2.0_dp * recSum * 4.0_dp * pi / vol
    sigma(:,:) = 2.0_dp * sigma * 4.0_dp * pi / vol

  end subroutine derivStressEwaldRec


  !> Returns the real space part of the Ewald sum.
  function ewaldReal(rr, rVec, alpha, blurWidth, epsSoften) result(realSum)

    !> Real space vectors to sum over. (Should contain origin).
    real(dp), intent(in) :: rVec(:,:)

    !> Parameter for the Ewald summation.
    real(dp), intent(in) :: alpha

    !> Vector where to calculate the Ewald sum.
    real(dp), intent(in) :: rr(:)

    !> Gaussian blur width of the 2nd charge
    real(dp), intent(in), optional :: blurWidth

    !> Short distance softening
    real(dp), intent(in), optional :: epsSoften

    !> contribution to sum
    real(dp) :: realSum

    real(dp) :: absRR, epsSoften2
    integer :: iR

    realSum = 0.0_dp

     if (present(epsSoften)) then
      epsSoften2 = epsSoften**2
    else
      epsSoften2 = 0.0_dp
    end if

    if (present(blurWidth)) then
      do iR = 1, size(rVec, dim=2)
        absRR = sum((rr(:) + rVec(:,iR))**2)
        if (absRR < tolSameDist**2) then
          cycle
        end if
        if (absRR < (erfArgLimit * blurWidth)**2) then
          realSum = realSum + erfwrap(sqrt(absRR) / blurWidth)/sqrt(absRR+epsSoften2)
        else
          realSum = realSum + 1.0_dp/sqrt(absRR+epsSoften2)
        end if
        realSum = realSum -erfwrap(alpha*sqrt(absRR))/sqrt(absRR+epsSoften2)
      end do
    else
      do iR = 1, size(rVec, dim=2)
        absRR = sum((rr(:) + rVec(:,iR))**2)
        if (absRR < tolSameDist**2) then
          cycle
        end if
        realSum = realSum + erfcwrap(alpha*sqrt(absRR))/sqrt(absRR+epsSoften2)
      end do
    end if

  end function ewaldReal


  !> Returns the derivative of the real space part of the Ewald sum.
  function derivEwaldReal(rdiff, rVec, alpha, blurWidth) result(dewr)

    !> Vector where to calculate the Ewald sum.
    real(dp), intent(in) :: rdiff(:)

    !> Real space vectors to sum over. (Should contain origin).
    real(dp), intent(in) :: rVec(:,:)

    !> Parameter for the Ewald summation.
    real(dp), intent(in) :: alpha

    !> Gaussian blur width of the second charge
    real(dp), intent(in), optional :: blurWidth

    !> contribution to derivative
    real(dp) :: dewr(3)

    real(dp) :: rNew(3)
    real(dp) :: rr
    integer :: iR

    dewr = 0.0_dp

    if (present(blurWidth)) then
      do iR = 1, size(rVec, dim=2)
        rNew(:) = rdiff(:) + rVec(:,iR)
        rr = sqrt(sum(rNew**2))
        if (rr < tolSameDist2) then
          cycle
        end if
        ! derivative of -erf(alpha*r)/r
        dewr(:) = dewr + rNew(:) * (-2.0_dp/sqrt(pi)*exp(-alpha*alpha*rr*rr)* alpha*rr&
            & - erfcwrap(alpha*rr))/(rr*rr*rr)
        ! deriv of erf(r/blur)/r
        if (rr < erfArgLimit * blurWidth) then
          dewr(:) = dewr + rNew(:) * (2.0_dp/sqrt(pi)*exp(-rr*rr/(blurWidth**2))*rr/blurWidth&
              & + erfcwrap(rr/blurWidth))/(rr*rr*rr)
        end if
      end do
    else
      do iR = 1, size(rVec, dim=2)
        rNew(:) = rdiff(:) + rVec(:,iR)
        rr = sqrt(sum(rNew**2))
        if (rr < tolSameDist2) then
          cycle
        end if
        dewr(:) = dewr + rNew(:) * (-2.0_dp/sqrt(pi)*exp(-alpha*alpha*rr*rr)* alpha*rr&
            & - erfcwrap(alpha*rr))/(rr*rr*rr)
      end do
    end if

  end function derivEwaldReal


  !> Returns the difference in the decrease of the real and reciprocal parts of the Ewald sum.  In
  !> order to make the real space part shorter than the reciprocal space part, the values are taken
  !> at different distances for the real and the reciprocal space parts.
  function diffRecReal(alpha, minG, minR, volume) result(diff)

    !> Parameter for the Ewald summation.
    real(dp), intent(in) :: alpha

    !> Length of the shortest reciprocal space vector in the sum.
    real(dp), intent(in) :: minG

    !> Length of the shortest real space vector in the sum.
    real(dp), intent(in) :: minR

    !> Volume of the real space unit cell.
    real(dp), intent(in) :: volume

    !> difference between changes in the two terms
    real(dp) :: diff

    @:ASSERT(volume > 0.0_dp)

    diff = ((gTerm(4.0_dp*minG, alpha, volume) &
        &- gTerm(5.0_dp*minG, alpha, volume))) &
        &- (rTerm(2.0_dp*minR, alpha) - rTerm(3.0_dp*minR, alpha))

  end function diffRecReal


  !> Returns the max. value of a term in the reciprocal space part of the Ewald summation for a
  !> given vector length.
  function gTerm(gg, alpha, vol)

    !> Length of the reciprocal space vector.
    real(dp), intent(in) :: gg

    !> Parameter of the Ewald summation.
    real(dp), intent(in) :: alpha

    !> Volume of the real space unit cell.
    real(dp), intent(in) :: vol

    !> reciprocal term
    real(dp) :: gTerm

    gTerm = 4.0_dp*pi*(exp(-0.25_dp*gg**2/(alpha**2))/(vol*gg**2))

  end function gTerm


  !> Returns the max. value of a term in the real space part of the Ewald summation for a given
  !> vector length.
  function rTerm(rr, alpha)

    !> Length of the real space vector.
    real(dp), intent(in) :: rr

    !> Parameter of the Ewald summation.
    real(dp), intent(in) :: alpha

    !> real space term
    real(dp) :: rTerm

    @:ASSERT(rr >= epsilon(1.0_dp))

    rTerm = erfcwrap(alpha*rr)/rr

  end function rTerm


  !> Returns the derivative of a term in the real space part of the Ewald summation for a given
  !> vector length.
  function derivRTerm(r, alpha)

    !> Length of the real space vector.
    real(dp), intent(in) :: r(3)

    !> Parameter of the Ewald summation.
    real(dp), intent(in) :: alpha

    !> real space derivative term
    real(dp) :: derivRTerm(3)

    real(dp) :: rr
    rr = sqrt(sum(r(:)**2))

    @:ASSERT(rr >= epsilon(1.0_dp))

    derivRTerm (:) = r(:)*(-2.0_dp/sqrt(pi)*exp(-alpha*alpha*rr*rr)* &
        & alpha*rr - erfcwrap(alpha*rr))/(rr*rr*rr)

  end function derivRTerm


  !> Calculates the stress tensor derivatives of the Ewald electrostatics
  !> Aguard and Madden J Chem Phys 119 7471 (2003)
  subroutine invRStress(env, nAtom, coord, neighList, recPoint, alpha, volume, q, stress)

    !> Computational environment settings
    type(TEnvironment), intent(in) :: env

    !> Number of atoms.
    integer, intent(in) :: nAtom

    !> List of atomic coordinates (all atoms).
    real(dp), intent(in) :: coord(:,:)

    !> Dynamic neighbour list to be used in the real part of Ewald
    type(TDynNeighList), target, intent(in) :: neighList

    !> Contains the points included in the reciprocal sum. The set should not include the origin or
    !> inversion related points.
    real(dp), intent(in) :: recPoint(:,:)

    !> Parameter for Ewald summation.
    real(dp), intent(in) :: alpha

    !> Volume of the real space unit cell.
    real(dp), intent(in) :: volume

    !> charges in the cell
    real(dp), intent(in) :: q(:)

    !> Stress tensor
    real(dp), intent(out) :: stress(:,:)

    type(TDynNeighList), pointer :: pNeighList
    integer :: iAtom1, iInv, ii, jj, kk
    real(dp) :: g(3), g2, intermed, intermed2
    real(dp) :: stressTmp(3,3), localStress(3,3)
    integer :: iFirst, iLast

    @:ASSERT(volume > 0.0_dp)

    ! Reciprocal space part of the Ewald sum.
    call distributeRangeInChunks(env, 1, size(recpoint, dim=2), iFirst, iLast)
    localStress(:,:) = 0.0_dp
    do ii = iFirst, iLast
      do iInv = -1, 1, 2
        g(:) = real(iInv,dp)*recpoint(:,ii)
        intermed = 0.0_dp
        intermed2 = 0.0_dp
        !$OMP PARALLEL DO&
        !$OMP& DEFAULT(SHARED) REDUCTION(+:intermed, intermed2) SCHEDULE(RUNTIME)
        do iAtom1 = 1, nAtom
          intermed = intermed + q(iAtom1) * cos(dot_product(g, coord(:,iAtom1)))
          intermed2 = intermed2 + q(iAtom1) * sin(dot_product(g, coord(:,iAtom1)))
        end do
        !$OMP END PARALLEL DO
        intermed = intermed**2 + intermed2**2
        g2 = sum(g(:)**2)
        intermed = intermed*exp(-g2/(4.0_dp*alpha*alpha))/g2
        stressTmp(:,:) = 0.0_dp
        do jj = 1, 3
          stressTmp(jj,jj) = 1.0_dp
          do kk = 1,3
            stressTmp(kk,jj) = stressTmp(kk,jj) &
                & -2.0_dp*(1.0_dp/(4.0_dp*alpha*alpha) + 1.0_dp/g2) &
                & *g(kk)*g(jj)
          end do
        end do
        localStress(:,:) = localStress + stressTmp * intermed
      end do
    end do
    localStress(:,:) = -localStress * 4.0_dp * pi / volume

    call assembleChunks(env, localStress)
    stress(:,:) = localStress


    ! Real space part of the Ewald sum.
    pNeighList => neighList
    call distributeRangeInChunks(env, 1, nAtom, iFirst, iLast)
    localStress = 0.0_dp
    !$OMP PARALLEL DO&
    !$OMP& DEFAULT(SHARED) REDUCTION(+:localStress) SCHEDULE(RUNTIME)
    do iAtom1 = iFirst, iLast
      call addNeighbourContribsStress(iAtom1, pNeighList, coord, alpha, Q, localStress)
    end do
    !$OMP END PARALLEL DO

    call assembleChunks(env, localStress)
    stress(:,:) = stress + localStress

    stress(:,:) = stress / volume

  end subroutine invRStress

  !> Neighbour summation with local scope for predictable OMP <= 4.0 behaviour
  subroutine addNeighbourContribsStress(iAtom1, pNeighList, coords, alpha, dQAtom, stress)
    integer, intent(in) :: iAtom1
    type(TDynNeighList), pointer, intent(in) :: pNeighList
    real(dp), intent(in) :: coords(:,:)
    real(dp), intent(in) :: dQAtom(:)
    real(dp), intent(in) :: alpha
    real(dp), intent(inout) :: stress(:,:)

    type(TNeighIterator) :: neighIter
    real(dp) :: neighCoords(3, iterChunkSize)
    integer :: neighImages(iterChunkSize)
    integer :: iAtom2f, iNeigh, nNeigh, ii, jj
    real(dp) :: r(3), f(3)

    call TNeighIterator_init(neighIter, pNeighList, iAtom1)
    nNeigh = iterChunkSize
    do while (nNeigh == iterChunkSize)
      call neighIter%getNextNeighbours(nNeigh, coords=neighCoords, img2CentCell=neighImages)
      do iNeigh = 1, nNeigh
        iAtom2f = neighImages(iNeigh)
        r(:) = coords(:,iAtom1) - neighCoords(:,iNeigh)
        f(:) = derivRTerm(r, alpha) * dQAtom(iAtom1) * dQAtom(iAtom2f)
        if (iAtom2f /= iAtom1) then
          do ii = 1, 3
            do jj = 1, 3
              stress(jj,ii) = stress(jj,ii) + (r(jj) * f(ii) + f(jj) * r(ii))
            end do
          end do
        else
          do ii = 1, 3
            do jj = 1, 3
              stress(jj,ii) = stress(jj,ii) + 0.5_dp * (r(jj) * f(ii) + f(jj) * r(ii))
            end do
          end do
        end if
      end do
    end do

  end subroutine addNeighbourContribsStress


  !> Calculates the -1/R**2 deriv contribution for all atoms for the non-periodic case, without
  !> storing anything.
  subroutine addInvRPrimeClusterMat(this, env, coord, invRDeriv)

    class(TCoulombCont), intent(in) :: this

    !> Computational environment settings
    type(TEnvironment), intent(in) :: env

    !> List of atomic coordinates.
    real(dp), intent(in) :: coord(:,:)

    !> derivative of inverse R matrix
    real(dp), intent(inout) :: invRDeriv(:,:,:)

    real(dp) :: dist, vect(3)
    integer :: nAtom, iAtFirst, iAtLast
    integer :: ii, jj

    nAtom = size(invRDeriv,dim=1)

    call distributeRangeInChunks(env, 1, nAtom, iAtFirst, iAtLast)

    !$OMP PARALLEL DO&
    !$OMP& DEFAULT(SHARED) PRIVATE(jj, vect, dist) REDUCTION(+:invRDeriv) SCHEDULE(RUNTIME)
    do ii = iAtFirst, iAtLast
      do jj = ii + 1, nAtom
        vect(:) = coord(:,ii) - coord(:,jj)
        dist = sqrt(sum(vect(:)**2))
        invRDeriv(jj,ii,:) = -vect(:) / (dist**3)
      end do
    end do
    !$OMP END PARALLEL DO

    call assembleChunks(env, invRDeriv)

  end subroutine addInvRPrimeClusterMat


  !> Calculates the -1/R**2 deriv contribution for the periodic case, without storing anything.
  subroutine addInvRPrimePeriodicMat(this, env, coord, neighList, recPoint, alpha, volume,&
      & invRDeriv)

    class(TCoulombCont), intent(in) :: this

    !> Computational environment settings
    type(TEnvironment), intent(in) :: env

    !> List of atomic coordinates (all atoms).
    real(dp), intent(in) :: coord(:,:)

    !> Dynamic neighbour list to be used in the real part of Ewald
    type(TDynNeighList), target, intent(in) :: neighList

    !> Contains the points included in the reciprocal sum. The set should not include the origin or
    !> inversion related points.
    real(dp), intent(in) :: recPoint(:,:)

    !> Parameter for Ewald summation.
    real(dp), intent(in) :: alpha

    !> Volume of the real space unit cell.
    real(dp), intent(in) :: volume

    !> derivative of inverse R matrix
    real(dp), intent(inout) :: invRDeriv(:,:,:)

    type(TDynNeighList), pointer :: pNeighList
    real(dp) :: r(3)
    integer :: nAtom, iAtom1, iAtom2
    integer :: iAtFirst, iAtLast

    nAtom = size(invRDeriv,dim=1)
    pNeighList => neighList

    call distributeRangeInChunks(env, 1, nAtom, iAtFirst, iAtLast)

    ! d(1/R)/dr real space
    !$OMP PARALLEL DO&
    !$OMP& DEFAULT(SHARED) REDUCTION(+:invRDeriv) SCHEDULE(RUNTIME)
    do iAtom1 = iAtFirst, iAtLast
      call addNeighbourContribsInvRPMat(iAtom1, pNeighList, &
          & coord, alpha, invRDeriv)
    end do
    !$OMP END PARALLEL DO

    ! d(1/R)/dr reciprocal space
    !$OMP PARALLEL DO&
    !$OMP& DEFAULT(SHARED) PRIVATE(iAtom2,r) REDUCTION(+:invRDeriv) SCHEDULE(RUNTIME)
    do iAtom1 = iAtFirst, iAtLast
      do iAtom2 = iAtom1+1, nAtom
        r(:) = coord(:,iAtom1) - coord(:,iAtom2)
        invRDeriv(iAtom2,iAtom1,:) = invRDeriv(iAtom2,iAtom1,:) &
            & + derivEwaldReciprocal(r,recPoint,alpha,volume)
      end do
    end do
    !$OMP END PARALLEL DO

    call assembleChunks(env, invRDeriv)

  end subroutine addInvRPrimePeriodicMat


  !> Neighbour summation with local scope for predictable OMP <= 4.0 behaviour
  subroutine addNeighbourContribsInvRPMat(iAtom1, pNeighList, coords, alpha, invRDeriv)
    integer, intent(in) :: iAtom1
    type(TDynNeighList), pointer, intent(in) :: pNeighList
    real(dp), intent(in) :: coords(:,:)
    real(dp), intent(in) :: alpha
    real(dp), intent(inout) :: invRDeriv(:,:,:)

    type(TNeighIterator) :: neighIter
    real(dp) :: neighCoords(3, iterChunkSize)
    integer :: neighImages(iterChunkSize)
    integer :: iAtom2f, iNeigh, nNeigh
    real(dp) :: rr(3)

    call TNeighIterator_init(neighIter, pNeighList, iAtom1)
    nNeigh = iterChunkSize
    do while (nNeigh == iterChunkSize)
      call neighIter%getNextNeighbours(nNeigh, coords=neighCoords, img2CentCell=neighImages)
      do iNeigh = 1, nNeigh
        iAtom2f = neighImages(iNeigh)
        if (iAtom2f /= iAtom1) then
          rr(:) = coords(:,iAtom1) - neighCoords(:,iNeigh)
          invRDeriv(iAtom2f,iAtom1,:) = &
              & invRDeriv(iAtom2f,iAtom1,:) + derivRTerm(rr,alpha)
        end if
      end do
    end do

  end subroutine addNeighbourContribsInvRPMat


  !> Get the variables relate to periodic information
  subroutine getPeriodicInfo(this, rVec, gVec, alpha, vol)

    !> Instance
    class(TCoulombCont), intent(in) :: this

    !> real lattice points for Ewald-sum
    real(dp), allocatable, intent(out) :: rVec(:,:)

    !> lattice points for reciprocal Ewald
    real(dp), allocatable, intent(out) :: gVec(:,:)

    !> parameter for Ewald
    real(dp), intent(out) :: alpha

    !> parameter for cell volume
    real(dp), intent(out) :: vol

    ! Use of automatic reallocation for gVec and rVec
    gVec = this%gLatPoint
    rVec = this%rCellVec

    alpha = this%alpha
    vol = this%volume

  end subroutine getPeriodicInfo


end module dftbp_coulomb
