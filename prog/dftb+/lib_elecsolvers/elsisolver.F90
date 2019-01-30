!--------------------------------------------------------------------------------------------------!
!  DFTB+: general package for performing fast atomistic simulations                                !
!  Copyright (C) 2018  DFTB+ developers group                                                      !
!                                                                                                  !
!  See the LICENSE file for terms of usage and distribution.                                       !
!--------------------------------------------------------------------------------------------------!

#:include 'common.fypp'

!> Contains the interface to the ELSI solvers
module elsisolver
#:if WITH_MPI
  use mpifx
#:endif
  use accuracy, only : dp, lc
  use environment, only : TEnvironment, globalTimers
  use elecsolvertypes, only : electronicSolverTypes
  use elsiiface
  use densedescr
  use periodic
  use orbitals
  use message, only : error
  use commontypes, only : TParallelKS, TOrbitals
  use energies, only : TEnergies
  use sparse2dense
  use constants, only : pi
  use assert
  use message
  use spin, only : ud2qm
  use angmomentum, only : getLOnsite
  use spinorbit, only : addOnsiteSpinOrbitHam, getOnsiteSpinOrbitEnergy
  use potentials, only : TPotentials
  implicit none
  private

  public :: TElsiSolverInp
  public :: TElsiSolver
  public :: TElsiSolver_init, TElsiSolver_final
  public :: TElsiSolver_getDensity
  public :: initPexsiDeltaVRanges
  public :: updatePexsiDeltaVRanges
  public :: getEDensityMtxElsiPauli, getEDensityMtxElsiReal, getEDensityMtxElsiCmplx
  public :: getEDensityRealElsi, getEDensityComplexElsi
  public :: TSparse2Sparse, initSparse2Sparse


  !> Input data for the ELSI solvers
  type :: TElsiSolverInp

    !> Choice of the solver
    integer :: iSolver

    !> Choice of ELPA solver
    integer :: ELPA_Solver = 2

    !> Iterations of ELPA solver before OMM minimization
    integer :: OMM_IterationsELPA = 5

    !> Halting tolerance for OMM iterations
    real(dp) :: OMM_Tolerance = 1.0E-10_dp

    !> Should the overlap be factorized before minimization
    logical :: OMM_Choleskii = .true.

    !> number of poles for PEXSI expansion
    integer :: PEXSI_n_pole = 20

    !> number of processors per pole for PEXSI
    integer :: PEXSI_np_per_pole = 1

    !> number of interpolation points for mu (Fermi) search
    integer :: PEXSI_n_mu = 2

    !> number of processors for symbolic factorisation
    integer :: PEXSI_np_symbo = 1

    !> spectral radius (range of eigenvalues) if available
    real(dp) :: PEXSI_delta_e = 10.0_dp

    !> density matrix purification algorithm
    integer :: NTPoly_method = 2

    !> truncation threshold for sparse matrix multiplication
    real(dp) :: NTPoly_truncation = 1.0E-10_dp

    !> convergence tolerance for density matrix purification
    real(dp) :: NTPoly_tolerance = 1.0E-5_dp

    !> Use sparse CSR format
    logical :: ELSI_CSR = .false.

  end type TElsiSolverInp


  !> Contains settings for the solvers of the ELSI library
  type :: TElsiSolver

    !> Solver id (using the central DFTB+ solver type)
    integer :: iSolver

    !> Handle for the ELSI library
    type(elsi_handle), public :: handle

    !> Use sparse CSR format
    logical, public :: CSR = .false.

    !> solver number (according to ELSI enumeration)
    integer :: SOLVER

    !> level of output from solver
    integer :: outputLevel

    !> should the code write matrices and stop
    logical :: tWriteHS

    !> handle for matrix IO
    type(elsi_rw_handle) :: rwHandle

    !> parallelisation strategy
    integer :: parallel

    !> Dense BLACS
    integer :: BLACS_DENSE

    integer :: n_basis
    real(dp) :: n_electron
    real(dp) :: spin_degeneracy
    integer :: n_state
    integer :: MPI_COMM_WORLD
    integer :: my_COMM_WORLD
    integer :: my_BLACS_Ctxt
    integer :: BLACS_blockSize
    integer :: CSR_blockSize

    integer :: n_spin
    integer :: n_kpoint

    integer :: mu_broaden_scheme
    integer :: mu_mp_order

    ! ELPA settings
    integer, public :: ELPA_SOLVER_Option

    ! OMM settings
    integer, public :: OMM_iter
    real(dp), public :: OMM_Tolerance
    logical, public :: OMM_Choleskii

    ! PEXSI settings
    real(dp), public :: PEXSI_mu_min
    real(dp), public :: PEXSI_mu_max
    real(dp), public :: PEXSI_DeltaVmin
    real(dp), public :: PEXSI_DeltaVmax
    real(dp), allocatable, public :: PEXSI_VOld(:)
    integer, public :: PEXSI_n_pole
    integer, public :: PEXSI_np_per_pole
    integer, public :: PEXSI_n_mu
    integer, public :: PEXSI_np_symbo
    real(dp) :: PEXSI_delta_e

    ! NTPoly settings
    integer, public :: NTPoly_method
    real(dp), public :: NTPoly_truncation
    real(dp), public :: NTPoly_tolerance

    !> count of the number of times ELSI has been reset (usually every geometry step)
    integer :: nResets = 0

    logical :: tCholeskyDecomposed

  contains
    procedure :: reset => TElsiSolver_reset
    procedure :: getSolverName => TElsiSolver_getSolverName
    procedure, private :: ensurePauliCompatibility => TElsiSolver_ensurePauliCompatibility
  end type TElsiSolver


  !> data type for indexing transformation from DFTB+ compressed block to CSR
  type :: TSparse2Sparse

    logical :: tInit = .false.

    !> number of non-zero matrix elements on this processor
    integer :: nnzLocal

    !>  numberof non-zero matrix elements in the whole sparse matrix
    integer :: nnzGlobal

    !> On which column does this processor start its matrix
    integer :: colStartLocal

    !> On which column does this processor end its matrix
    integer :: colEndLocal

    !> Number of local columns on this processor
    integer :: numColLocal

    !> Local column pointer in the CSC format on this processor
    integer, allocatable :: colPtrLocal(:)

    !> Index for starting row of blocks in nzValLocal
    integer, allocatable :: blockRow(:,:)

    !> Local row index for non-zero elements
    integer, allocatable :: rowIndLocal(:)

    !> List of atoms with elements in the columns held locally
    integer, allocatable :: atomsInColumns(:)

    !> Count of the atoms with elements in the columns held locally
    integer :: nAtomsInColumns

  end type TSparse2Sparse


  #:set FLAVOURS = [('real'), ('complex')]

  interface addBlock2Elsi
  #:for SUFFIX in FLAVOURS
    module procedure addBlock2Elsi${SUFFIX}$
  #:endfor
  end interface addBlock2Elsi

  interface cpElsi2Block
  #:for SUFFIX in FLAVOURS
    module procedure cpElsi2Block${SUFFIX}$
  #:endfor
  end interface cpElsi2Block

contains


  !> Initialise ELSI solver
  subroutine TElsiSolver_init(this, inp, env, nBasisFn, nEl, iDistribFn, nSpin, nKPoint, tWriteHS)

    !> control structure for solvers, including ELSI data
    class(TElsiSolver), intent(out) :: this

    !> input structure for ELSI
    type(TElsiSolverInp), intent(in) :: inp

    !> Environment settings
    type(TEnvironment), intent(in) :: env

    !> number of orbitals in the system
    integer, intent(in) :: nBasisFn

    !> number of electrons
    real(dp), intent(in) :: nEl(:)

    !> filling function
    integer, intent(in) :: iDistribFn

    !> total number of spin channels.
    integer, intent(in) :: nSpin

    !> total number of k-points
    integer, intent(in) :: nKPoint

    !> Should the matrices be written out
    logical, intent(in) :: tWriteHS

    this%iSolver = inp%iSolver

    select case(this%iSolver)

    case (electronicSolverTypes%elpa)
      this%SOLVER = 1
      ! ELPA is asked for all states
      this%n_state = nBasisFn

    case (electronicSolverTypes%omm)
      this%SOLVER = 2
      ! OMM solves only over occupied space
      ! spin degeneracies for closed shell
      if (nSpin == 1) then
        this%n_state = nint(sum(nEl)*0.5_dp)
      else
        this%n_state = nint(sum(nEl))
      end if

    case (electronicSolverTypes%pexsi)
      this%SOLVER = 3
      ! ignored by PEXSI:
      this%n_state = nBasisFn

    case (electronicSolverTypes%ntpoly)
      this%SOLVER = 6
      ! ignored by NTPoly:
      this%n_state = nBasisFn

    end select

    ! parallelism with multiple processes
    this%parallel = 1

    ! number of basis functions
    this%n_basis = nBasisFn

    ! number of electrons in the problem
    this%n_electron = sum(nEl)

    ! should the code write the hamiltonian and overlap and stop
    this%tWriteHS = tWriteHS

    if (nSpin == 2) then
      this%spin_degeneracy = nEl(1) - nEl(2)
    else
      this%spin_degeneracy = 0.0_dp
    end if

    this%n_spin = nSpin
    this%n_kpoint = nKPoint

    this%MPI_COMM_WORLD = env%mpi%globalComm%id
    this%my_COMM_WORLD = env%mpi%groupComm%id
    this%my_BLACS_Ctxt = env%blacs%orbitalGrid%ctxt

    ! assumes row and column sizes the same
    this%BLACS_blockSize = env%blacs%rowBlockSize

    this%CSR_blockSize = this%n_basis / env%mpi%groupComm%size

    if (this%CSR_blockSize * env%mpi%groupComm%size < this%n_basis) then
      this%CSR_blockSize = this%CSR_blockSize + 1
    end if

    this%mu_broaden_scheme = min(iDistribFn,2)
    if (iDistribFn > 1) then
      ! set Meth-Pax order
      this%mu_mp_order = iDistribFn - 2
    else
      this%mu_mp_order = 0
    end if

    ! ELPA settings
    this%ELPA_SOLVER_Option = inp%ELPA_Solver

    ! OMM settings
    this%OMM_iter = inp%OMM_IterationsELPA
    this%OMM_Tolerance = inp%OMM_Tolerance
    this%OMM_Choleskii = inp%OMM_Choleskii

    ! PEXSI settings
    this%PEXSI_n_pole = inp%PEXSI_n_pole
    this%PEXSI_np_per_pole = inp%PEXSI_np_per_pole
    this%PEXSI_n_mu = inp%PEXSI_n_mu
    this%PEXSI_np_symbo = inp%PEXSI_np_symbo
    this%PEXSI_delta_e = inp%PEXSI_delta_e

    ! NTPoly settings
    this%NTPoly_method = inp%NTPoly_method
    this%NTPoly_truncation = inp%NTPoly_truncation
    this%NTPoly_tolerance = inp%NTPoly_tolerance

    this%CSR = inp%ELSI_CSR

    ! data as dense BLACS blocks
    if (this%CSR) then
      ! CSR format
      this%BLACS_DENSE = 2
    else
      this%BLACS_DENSE = 0
    end if

    ! customize output level, note there are levels 0..3 not DFTB+ 0..2
    this%OutputLevel = 0
  #:call DEBUG_CODE
    this%OutputLevel = 3
  #:endcall DEBUG_CODE

    if (nSpin == 4) then
      call this%ensurePauliCompatibility()
    end if

    this%tCholeskyDecomposed = .false.

  end subroutine TElsiSolver_init


  !> Finalizes the ELSI solver.
  subroutine TELsiSolver_final(this)

    !> Instance
    type(TElsiSolver), intent(inout) :: this

    call elsi_finalize(this%handle)

  end subroutine TELsiSolver_final


  !> reset the ELSI solver - safer to do this on geometry change, due to the lack of a Choleskii
  !> refactorization option
  subroutine TElsiSolver_reset(this, tempElec, iSpin, iKPoint, kWeight)

    !> Instance
    class(TElsiSolver), intent(inout) :: this

    !> electron temperature
    real(dp), intent(in) :: tempElec

    !> current spin value, set to be 1 if non-collinear
    integer, intent(in) :: iSpin

    !> current k-point value
    integer, intent(in) :: iKPoint

    !> weight for current k-point
    real(dp), intent(in) :: kWeight

    if (this%nResets > 0) then
      ! destroy previous instance of solver if called before
      call elsi_finalize(this%handle)
    end if
    this%nResets = this%nResets + 1

    this%tCholeskyDecomposed = .false.

    call elsi_init(this%handle, this%SOLVER, this%parallel, this%BLACS_DENSE, this%n_basis,&
        & this%n_electron, this%n_state)

    call elsi_set_mpi_global(this%handle, this%MPI_COMM_WORLD)
    call elsi_set_sing_check(this%handle, 0) ! disable singularity check
    call elsi_set_mpi(this%handle, this%my_COMM_WORLD)

    if (this%CSR) then
      call elsi_set_csc_blk(this%handle, this%CSR_blockSize)
    end if
    call elsi_set_blacs(this%handle, this%my_BLACS_Ctxt, this%BLACS_blockSize)

    if (this%tWriteHS) then
      ! setup to write a matrix
      call elsi_init_rw(this%rwHandle, 1, this%parallel, this%n_basis, this%n_electron)
      ! MPI comm
      call elsi_set_rw_mpi(this%rwHandle, this%MPI_COMM_WORLD)
      if (.not.this%CSR) then
        ! dense matrices
        call elsi_set_rw_blacs(this%rwHandle, this%my_BLACS_Ctxt, this%BLACS_blockSize)
      end if
    end if

    call elsi_set_mu_broaden_scheme(this%handle, this%mu_broaden_scheme)
    if (this%mu_broaden_scheme == 2) then
      call elsi_set_mu_mp_order(this%handle, this%mu_mp_order)
    end if

    ! set filling temperature/width
    call elsi_set_mu_broaden_width(this%handle, tempElec)

    select case(this%iSolver)
    case(electronicSolverTypes%elpa)

      select case(this%ELPA_SOLVER_Option)
      case(1)
        call elsi_set_elpa_solver(this%handle, 1)
      case(2)
        call elsi_set_elpa_solver(this%handle, 2)
      case default
        call error("Unknown ELPA solver modes")
      end select

    case(electronicSolverTypes%omm)
      ! libOMM
      if (this%OMM_Choleskii) then
        call elsi_set_omm_flavor(this%handle, 2)
      else
        call elsi_set_omm_flavor(this%handle, 0)
      end if
      call elsi_set_omm_n_elpa(this%handle, this%OMM_iter)
      call elsi_set_omm_tol(this%handle, this%OMM_Tolerance)

    case(electronicSolverTypes%pexsi)

      this%PEXSI_mu_min = -10.0_dp
      this%PEXSI_mu_max = 10.0_dp
      this%PEXSI_DeltaVmin = 0.0_dp
      this%PEXSI_DeltaVmax = 0.0_dp

      ! processors per pole to invert for
      call elsi_set_pexsi_np_per_pole(this%handle, this%PEXSI_np_per_pole)

      call elsi_set_pexsi_mu_min(this%handle, this%PEXSI_mu_min +this%PEXSI_DeltaVmin)
      call elsi_set_pexsi_mu_max(this%handle, this%PEXSI_mu_max +this%PEXSI_DeltaVmax)

      ! set filling temperature/width
      call elsi_set_pexsi_temp(this%handle, tempElec)

      ! number of poles for the expansion
      call elsi_set_pexsi_n_pole(this%handle, this%PEXSI_n_pole)

      ! number of interpolation points for mu
      call elsi_set_pexsi_n_mu(this%handle, this%PEXSI_n_mu)

      ! number of processors for symbolic factorisation task
      call elsi_set_pexsi_np_symbo(this%handle, this%PEXSI_np_symbo)

      ! spectral radius (range of eigenspectrum, if known, otherwise defaul usually fine)
      call elsi_set_pexsi_delta_e(this%handle, this%PEXSI_delta_e)

    case(electronicSolverTypes%ntpoly)

      ! NTPoly
      ! set purification method
      call elsi_set_ntpoly_method(this%handle, this%NTPoly_method)

      ! set truncation tolerance for sparse matrix multiplications
      call elsi_set_ntpoly_filter(this%handle, this%NTPoly_truncation)

      ! set purification convergence threshold
      call elsi_set_ntpoly_tol(this%handle, this%NTPoly_tolerance)

    end select

    if (any(this%iSolver == [electronicSolverTypes%omm,&
        & electronicSolverTypes%pexsi, electronicSolverTypes%ntpoly])) then
      ! density matrix build needs to know the number of spin channels to normalize against
      select case(this%n_spin)
      case(1)
        call elsi_set_spin(this%handle, 1, iSpin)
      case(2)
        call elsi_set_spin(this%handle, 2, iSpin)
      case(4)
        call elsi_set_spin(this%handle, 2, iSpin)
      end select
      if (this%n_kpoint > 1) then
        call elsi_set_kpoint(this%handle, this%n_kpoint, iKPoint, kWeight)
      end if
    end if

    call elsi_set_output(this%handle, this%OutputLevel)
    if (this%OutputLevel == 3) then
      call elsi_set_output_log(this%handle, 1)
    end if

  end subroutine TElsiSolver_reset


  function TElsiSolver_getSolverName(this) result(solverName)

    !> Instance.
    class(TElsiSolver), intent(in) :: this

    !> Name of the solver.
    character(:), allocatable :: solverName

    character(lc) :: buffer

    select case (this%iSolver)

    case(electronicSolverTypes%elpa)
      if (this%ELPA_SOLVER_Option == 1) then
        write(buffer, "(A)") "ELSI interface to the 1 stage ELPA solver"
      else
        write(buffer, "(A)") "ELSI interface to the 2 stage ELPA solver"
      end if

    case(electronicSolverTypes%omm)
      write(buffer, "(A,I0,A,E8.2)") "ELSI solver libOMM with ",&
          & this%OMM_iter, " ELPA iterations",this%OMM_Tolerance
      if (this%CSR) then
        write(buffer, "(A)") "ELSI solver libOMM Sparse"
      else
        write(buffer, "(A)") "ELSI solver libOMM Dense"
      end if

    case(electronicSolverTypes%pexsi)
      if (this%CSR) then
        write(buffer, "(A)") "ELSI solver PEXSI Sparse"
      else
        write(buffer, "(A)") "ELSI solver PEXSI Dense"
      end if

    case(electronicSolverTypes%ntpoly)
      if (this%CSR) then
        write(buffer, "(A)") "ELSI solver NTPoly Sparse"
      else
        write(buffer, "(A)") "ELSI solver NTPoly Dense"
      end if

    case default
      write(buffer, "(A,I0,A)") "Invalid electronic solver! (iSolver = ", this%iSolver, ")"

    end select

    solverName = trim(buffer)

  end function TElsiSolver_getSolverName


  !> Ensures Pauli compatibility
  subroutine TElsiSolver_ensurePauliCompatibility(this)
    class(TElsiSolver), intent(in) :: this

    if (this%CSR .and. this%solver /= 1) then
      call error("Current solver configuration not avaible for two component complex&
          & hamiltonians")
    end if

  end subroutine TElsiSolver_ensurePauliCompatibility



  !> Returns the density matrix using ELSI non-diagonalisation routines.
  subroutine TElsiSolver_getDensity(this, env, denseDesc, ham, over, neighbourList, nNeighbourSK,&
      & iSparseStart, img2CentCell, iCellVec, cellVec, kPoint, kWeight, orb, species,&
      & tRealHS, tSpinSharedEf, tSpinOrbit, tDualSpinOrbit, tFillKSep, tFixEf, tMulliken,&
      & iDistribFn, tempElec, nEl, parallelKS, Ef, mu, energy, eigen, filling, rhoPrim, Eband, TS,&
      & E0, iHam, xi, orbitalL, HSqrReal, SSqrReal, eigvecsReal, iRhoPrim, HSqrCplx, SSqrCplx,&
      & eigvecsCplx, tLargeDenseMatrices, sparseIndexing)

    !> Electronic solver information
    type(TElsiSolver), intent(inout) :: this

    !> Environment settings
    type(TEnvironment), intent(inout) :: env

    !> Dense matrix descriptor
    type(TDenseDescr), intent(in) :: denseDesc

    !> hamiltonian in sparse storage
    real(dp), intent(in) :: ham(:,:)

    !> sparse overlap matrix
    real(dp), intent(in) :: over(:)

    !> list of neighbours for each atom
    type(TNeighbourList), intent(in) :: neighbourList

    !> Number of neighbours for each of the atoms
    integer, intent(in) :: nNeighbourSK(:)

    !> Index array for the start of atomic blocks in sparse arrays
    integer, intent(in) :: iSparseStart(:,:)

    !> map from image atoms to the original unique atom
    integer, intent(in) :: img2CentCell(:)

    !> Index for which unit cell atoms are associated with
    integer, intent(in) :: iCellVec(:)

    !> Vectors (in units of the lattice constants) to cells of the lattice
    real(dp), intent(in) :: cellVec(:,:)

    !> k-points
    real(dp), intent(in) :: kPoint(:,:)

    !> Weights for k-points
    real(dp), intent(in) :: kWeight(:)

    !> Atomic orbital information
    type(TOrbitals), intent(in) :: orb

    !> species of all atoms in the system
    integer, intent(in) :: species(:)

    !> Is the hamitonian real (no k-points/molecule/gamma point)?
    logical, intent(in) :: tRealHS

    !> Is the Fermi level common accross spin channels?
    logical, intent(in) :: tSpinSharedEf

    !> Are spin orbit interactions present
    logical, intent(in) :: tSpinOrbit

    !> Are block population spin orbit interactions present
    logical, intent(in) :: tDualSpinOrbit

    !> Fill k-points separately if true (no charge transfer accross the BZ)
    logical, intent(in) :: tFillKSep

    !> Whether fixed Fermi level(s) should be used. (No charge conservation!)
    logical, intent(in) :: tFixEf

    !> Should Mulliken populations be generated/output
    logical, intent(in) :: tMulliken

    !> occupation function for electronic states
    integer, intent(in) :: iDistribFn

    !> Electronic temperature
    real(dp), intent(in) :: tempElec

    !> Number of electrons
    real(dp), intent(in) :: nEl(:)

    !> K-points and spins to process
    type(TParallelKS), intent(in) :: parallelKS

    !> Fermi level(s)
    real(dp), intent(inout) :: Ef(:)

    !> Electrochemical potentials (contact, spin)
    real(dp), allocatable, intent(in) :: mu(:,:)

    !> Energy contributions and total
    type(TEnergies), intent(inout) :: energy

    !> eigenvalues (level, kpoint, spin)
    real(dp), intent(out) :: eigen(:,:,:)

    !> occupations (level, kpoint, spin)
    real(dp), intent(out) :: filling(:,:,:)

    !> sparse density matrix
    real(dp), intent(out) :: rhoPrim(:,:)

    !> band structure energy
    real(dp), intent(out) :: Eband(:)

    !> electronic entropy times temperature
    real(dp), intent(out) :: TS(:)

    !> extrapolated 0 temperature band energy
    real(dp), intent(out) :: E0(:)

    !> imaginary part of hamitonian
    real(dp), intent(in), allocatable :: iHam(:,:)

    !> spin orbit constants
    real(dp), intent(in), allocatable :: xi(:,:)

    !> orbital moments of atomic shells
    real(dp), intent(inout), allocatable :: orbitalL(:,:,:)

    !> imaginary part of density matrix
    real(dp), intent(inout), allocatable :: iRhoPrim(:,:)

    !> dense real hamiltonian storage
    real(dp), intent(inout), allocatable :: HSqrReal(:,:)

    !> dense real overlap storage
    real(dp), intent(inout), allocatable :: SSqrReal(:,:)

    !> real eigenvectors on exit
    real(dp), intent(inout), allocatable :: eigvecsReal(:,:,:)

    !> dense complex (k-points) hamiltonian storage
    complex(dp), intent(inout), allocatable :: HSqrCplx(:,:)

    !> dense complex (k-points) overlap storage
    complex(dp), intent(inout), allocatable :: SSqrCplx(:,:)

    !> complex eigenvectors on exit
    complex(dp), intent(inout), allocatable :: eigvecsCplx(:,:,:)

    !> Are dense matrices for H, S, etc. being used
    logical, intent(in) :: tLargeDenseMatrices

    !> sparse matrices indexing data structure
    type(TSparse2Sparse), intent(inout) :: sparseIndexing

    integer :: nSpin
    integer :: iKS, iK, iS

    ! as each spin and k-point combination forms a separate group for this solver, then iKS = 1
    iKS = 1
    iK = parallelKS%localKS(1, iKS)
    iS = parallelKS%localKS(2, iKS)

    nSpin = size(ham, dim=2)
    if (nSpin == 2 .and. .not. tSpinSharedEf) then
      call error("ELSI currently requires shared Fermi levels over spin")
    end if

    rhoPrim(:,:) = 0.0_dp
    if (allocated(iRhoPrim)) then
      iRhoPrim(:,:) = 0.0_dp
    end if

    if (this%iSolver == electronicSolverTypes%pexsi) then
      call elsi_set_pexsi_mu_min(this%handle, this%PEXSI_mu_min + this%PEXSI_DeltaVmin)
      call elsi_set_pexsi_mu_max(this%handle, this%PEXSI_mu_max + this%PEXSI_DeltaVmax)
    end if

    if (nSpin /= 4) then
      if (tRealHS) then
        if (tLargeDenseMatrices) then
          call getDensityFromElsiRealDense(env, denseDesc, ham, over, neighbourList, nNeighbourSK,&
              & iSparseStart, img2CentCell, orb, this, parallelKS, rhoPrim, Eband,&
              & HSqrReal, SSqrReal, eigvecsReal)
        else
          !call getDensityFromElsiRealSparse()
          call calcDensityRealElsi(sparseIndexing, parallelKS, this, ham, over,&
              & neighbourList%iNeighbour, nNeighbourSK, denseDesc%iAtomStart, iSparseStart,&
              & img2CentCell, orb, rhoPrim, Eband)
        end if
      else
        if (tLargeDenseMatrices) then
          call getDensityFromElsiCplxDense(env, denseDesc, ham, over, neighbourList, nNeighbourSK,&
              & iSparseStart, img2CentCell, iCellVec, cellVec, kPoint, kWeight, orb,&
              & this, parallelKS, rhoPrim, Eband, HSqrCplx, SSqrCplx, eigvecsCplx)
        else
          !call getDensityFromElsiCplxSparse()
          call calcDensityComplexElsi(sparseIndexing, parallelKS, this, kPoint(:,iK),&
              & kWeight(iK), iCellVec, cellVec, ham, over, neighbourList%iNeighbour,&
              & nNeighbourSK, denseDesc%iAtomStart, iSparseStart, img2CentCell, orb, rhoPrim,&
              & Eband)
        end if
      end if
      call ud2qm(rhoPrim)
    else
      if (tLargeDenseMatrices) then
        call getDensityFromElsiPauliDense(env, denseDesc, ham, over, neighbourList, nNeighbourSK,&
            & iSparseStart, img2CentCell, iCellVec, cellVec, kPoint, kWeight, orb, species,&
            & this, tSpinOrbit, tDualSpinOrbit, tMulliken, parallelKS, energy, rhoPrim,&
            & Eband, iHam, xi, orbitalL, iRhoPrim, HSqrCplx, SSqrCplx)
      else
        call error("Internal error: getDensityFromElsiSolver : sparse pauli not yet implemented")
      end if
    end if

    Ef(:) = 0.0_dp
    TS(:) = 0.0_dp
    if (env%mpi%tGroupMaster) then
      call elsi_get_mu(this%handle, Ef(iS))
      call elsi_get_entropy(this%handle, TS(iS))
    end if
    call mpifx_allreduceip(env%mpi%globalComm, Ef, MPI_SUM)
    call mpifx_allreduceip(env%mpi%globalComm, TS, MPI_SUM)

    if (this%iSolver == electronicSolverTypes%pexsi) then
      call elsi_get_pexsi_mu_min(this%handle, this%PEXSI_mu_min)
      call elsi_get_pexsi_mu_max(this%handle, this%PEXSI_mu_max)
    end if

    ! Add up and distribute density matrix contribution from each group
    call mpifx_allreduceip(env%mpi%globalComm, rhoPrim, MPI_SUM)
    if (allocated(iRhoPrim)) then
      call mpifx_allreduceip(env%mpi%globalComm, iRhoPrim, MPI_SUM)
    end if

    call env%globalTimer%stopTimer(globalTimers%densityMatrix)

  end subroutine TElsiSolver_getDensity


  !> Returns the density matrix using ELSI non-diagonalisation routines (real dense case).
  subroutine getDensityFromElsiRealDense(env, denseDesc, ham, over, neighbourList, nNeighbourSK,&
      & iSparseStart, img2CentCell, orb, elsiSolver, parallelKS, rhoPrim, Eband, HSqrReal,&
      & SSqrReal, eigvecsReal)

    !> Environment settings
    type(TEnvironment), intent(inout) :: env

    !> Dense matrix descriptor
    type(TDenseDescr), intent(in) :: denseDesc

    !> hamiltonian in sparse storage
    real(dp), intent(in) :: ham(:,:)

    !> sparse overlap matrix
    real(dp), intent(in) :: over(:)

    !> list of neighbours for each atom
    type(TNeighbourList), intent(in) :: neighbourList

    !> Number of neighbours for each of the atoms
    integer, intent(in) :: nNeighbourSK(:)

    !> Index array for the start of atomic blocks in sparse arrays
    integer, intent(in) :: iSparseStart(:,:)

    !> map from image atoms to the original unique atom
    integer, intent(in) :: img2CentCell(:)

    !> Atomic orbital information
    type(TOrbitals), intent(in) :: orb

    !> Electronic solver information
    type(TElsiSolver), intent(inout) :: elsiSolver

    !> K-points and spins to process
    type(TParallelKS), intent(in) :: parallelKS

    !> sparse density matrix
    real(dp), intent(inout) :: rhoPrim(:,:)

    !> band structure energy
    real(dp), intent(out) :: Eband(:)

    !> dense real hamiltonian storage
    real(dp), intent(inout), allocatable :: HSqrReal(:,:)

    !> dense real overlap storage
    real(dp), intent(inout), allocatable :: SSqrReal(:,:)

    !> real eigenvectors on exit
    real(dp), intent(inout), allocatable :: eigvecsReal(:,:,:)

    real(dp), allocatable :: rhoSqrReal(:,:)
    integer :: iKS, iS

    ! as each spin and k-point combination forms a separate group for this solver, then iKS = 1
    iKS = 1
    iS = parallelKS%localKS(2, iKS)

    call unpackHSRealBlacs(env%blacs, ham(:,iS), neighbourList%iNeighbour, nNeighbourSK,&
        & iSparseStart, img2CentCell, denseDesc, HSqrReal)
    if (.not. elsiSolver%tCholeskyDecomposed) then
      call unpackHSRealBlacs(env%blacs, over, neighbourList%iNeighbour, nNeighbourSK,&
          & iSparseStart, img2CentCell, denseDesc, SSqrReal)
      if (elsiSolver%OMM_Choleskii .or. elsiSolver%iSolver == electronicSolverTypes%pexsi) then
        elsiSolver%tCholeskyDecomposed = .true.
      end if
    end if

    allocate(rhoSqrReal(size(HSqrReal, dim=1), size(HSqrReal, dim=2)))
    rhoSqrReal(:,:) = 0.0_dp
    Eband(iS) = 0.0_dp
    if (elsiSolver%tWriteHS) then
      call elsi_write_mat_real(elsiSolver%rwHandle, "ELSI_Hreal.bin", HSqrReal)
      call elsi_write_mat_real(elsiSolver%rwHandle, "ELSI_Sreal.bin", SSqrReal)
      call elsi_finalize_rw(elsiSolver%rwHandle)
      call cleanShutdown("Finished dense matrix write")
    end if
    call elsi_dm_real(elsiSolver%handle, HSqrReal, SSqrReal, rhoSqrReal, Eband(iS))

    call packRhoRealBlacs(env%blacs, denseDesc, rhoSqrReal, neighbourList%iNeighbour,&
        & nNeighbourSK, orb%mOrb, iSparseStart, img2CentCell, rhoPrim(:,iS))

  end subroutine getDensityFromElsiRealDense


  !> Returns the density matrix using ELSI non-diagonalisation routines (complex dense case).
  subroutine getDensityFromElsiCplxDense(env, denseDesc, ham, over, neighbourList, nNeighbourSK,&
      & iSparseStart, img2CentCell, iCellVec, cellVec, kPoint, kWeight, orb, elsiSolver,&
      & parallelKS, rhoPrim, Eband, HSqrCplx, SSqrCplx, eigvecsCplx)

    !> Environment settings
    type(TEnvironment), intent(inout) :: env

    !> Dense matrix descriptor
    type(TDenseDescr), intent(in) :: denseDesc

    !> hamiltonian in sparse storage
    real(dp), intent(in) :: ham(:,:)

    !> sparse overlap matrix
    real(dp), intent(in) :: over(:)

    !> list of neighbours for each atom
    type(TNeighbourList), intent(in) :: neighbourList

    !> Number of neighbours for each of the atoms
    integer, intent(in) :: nNeighbourSK(:)

    !> Index array for the start of atomic blocks in sparse arrays
    integer, intent(in) :: iSparseStart(:,:)

    !> map from image atoms to the original unique atom
    integer, intent(in) :: img2CentCell(:)

    !> Index for which unit cell atoms are associated with
    integer, intent(in) :: iCellVec(:)

    !> Vectors (in units of the lattice constants) to cells of the lattice
    real(dp), intent(in) :: cellVec(:,:)

    !> k-points
    real(dp), intent(in) :: kPoint(:,:)

    !> Weights for k-points
    real(dp), intent(in) :: kWeight(:)

    !> Atomic orbital information
    type(TOrbitals), intent(in) :: orb

    !> Electronic solver information
    type(TElsiSolver), intent(inout) :: elsiSolver

    !> K-points and spins to process
    type(TParallelKS), intent(in) :: parallelKS

    !> sparse density matrix
    real(dp), intent(inout) :: rhoPrim(:,:)

    !> band structure energy
    real(dp), intent(out) :: Eband(:)

    !> electronic entropy times temperature
    !> dense complex (k-points) hamiltonian storage
    complex(dp), intent(inout), allocatable :: HSqrCplx(:,:)

    !> dense complex (k-points) overlap storage
    complex(dp), intent(inout), allocatable :: SSqrCplx(:,:)

    !> complex eigenvectors on exit
    complex(dp), intent(inout), allocatable :: eigvecsCplx(:,:,:)

    complex(dp), allocatable :: rhoSqrCplx(:,:)
    integer :: iKS, iK, iS

    ! as each spin and k-point combination forms a separate group for this solver, then iKS = 1
    iKS = 1
    iK = parallelKS%localKS(1, iKS)
    iS = parallelKS%localKS(2, iKS)

    HSqrCplx(:,:) = 0.0_dp
    call unpackHSCplxBlacs(env%blacs, ham(:,iS), kPoint(:,iK), neighbourList%iNeighbour,&
        & nNeighbourSK, iCellVec, cellVec, iSparseStart, img2CentCell, denseDesc, HSqrCplx)
    if (.not. elsiSolver%tCholeskyDecomposed) then
      SSqrCplx(:,:) = 0.0_dp
      call unpackHSCplxBlacs(env%blacs, over, kPoint(:,iK), neighbourList%iNeighbour,&
          & nNeighbourSK, iCellVec, cellVec, iSparseStart, img2CentCell, denseDesc,&
          & SSqrCplx)
      if (elsiSolver%OMM_Choleskii .or. elsiSolver%iSolver == electronicSolverTypes%pexsi) then
        elsiSolver%tCholeskyDecomposed = .true.
      end if
    end if

    allocate(rhoSqrCplx(size(HSqrCplx,dim=1),size(HSqrCplx,dim=2)))
    rhoSqrCplx(:,:) = 0.0_dp
    if (elsiSolver%tWriteHS) then
      call elsi_write_mat_complex(elsiSolver%rwHandle, "ELSI_Hcmplex.bin",&
          & HSqrCplx)
      call elsi_write_mat_complex(elsiSolver%rwHandle, "ELSI_Scmplx.bin",&
          & SSqrCplx)
      call elsi_finalize_rw(elsiSolver%rwHandle)
      call cleanShutdown("Finished dense matrix write")
    end if
    call elsi_dm_complex(elsiSolver%handle, HSqrCplx, SSqrCplx, rhoSqrCplx,&
        & Eband(iS))
    call packRhoCplxBlacs(env%blacs, denseDesc, rhoSqrCplx, kPoint(:,iK), kWeight(iK),&
        & neighbourList%iNeighbour, nNeighbourSK, orb%mOrb, iCellVec, cellVec,&
        & iSparseStart, img2CentCell, rhoPrim(:,iS))

  end subroutine getDensityFromElsiCplxDense


  !> Returns the density matrix using ELSI non-diagonalisation routines (Pauli dense case).
  subroutine getDensityFromElsiPauliDense(env, denseDesc, ham, over, neighbourList, nNeighbourSK,&
      & iSparseStart, img2CentCell, iCellVec, cellVec, kPoint, kWeight, orb, species,&
      & elsiSolver, tSpinOrbit, tDualSpinOrbit, tMulliken, parallelKS, energy, rhoPrim,&
      & Eband, iHam, xi, orbitalL, iRhoPrim, HSqrCplx, SSqrCplx)

    !> Environment settings
    type(TEnvironment), intent(inout) :: env

    !> Dense matrix descriptor
    type(TDenseDescr), intent(in) :: denseDesc

    !> hamiltonian in sparse storage
    real(dp), intent(in) :: ham(:,:)

    !> sparse overlap matrix
    real(dp), intent(in) :: over(:)

    !> list of neighbours for each atom
    type(TNeighbourList), intent(in) :: neighbourList

    !> Number of neighbours for each of the atoms
    integer, intent(in) :: nNeighbourSK(:)

    !> Index array for the start of atomic blocks in sparse arrays
    integer, intent(in) :: iSparseStart(:,:)

    !> map from image atoms to the original unique atom
    integer, intent(in) :: img2CentCell(:)

    !> Index for which unit cell atoms are associated with
    integer, intent(in) :: iCellVec(:)

    !> Vectors (in units of the lattice constants) to cells of the lattice
    real(dp), intent(in) :: cellVec(:,:)

    !> k-points
    real(dp), intent(in) :: kPoint(:,:)

    !> Weights for k-points
    real(dp), intent(in) :: kWeight(:)

    !> Atomic orbital information
    type(TOrbitals), intent(in) :: orb

    !> species of all atoms in the system
    integer, intent(in) :: species(:)

    !> Electronic solver information
    type(TElsiSolver), intent(inout) :: elsiSolver

    !> Are spin orbit interactions present
    logical, intent(in) :: tSpinOrbit

    !> Are block population spin orbit interactions present
    logical, intent(in) :: tDualSpinOrbit

    !> Should Mulliken populations be generated/output
    logical, intent(in) :: tMulliken

    !> K-points and spins to process
    type(TParallelKS), intent(in) :: parallelKS

    !> Energy contributions and total
    type(TEnergies), intent(inout) :: energy

    !> sparse density matrix
    real(dp), intent(inout) :: rhoPrim(:,:)

    !> band structure energy
    real(dp), intent(out) :: Eband(:)

    !> imaginary part of hamitonian
    real(dp), intent(in), allocatable :: iHam(:,:)

    !> spin orbit constants
    real(dp), intent(in), allocatable :: xi(:,:)

    !> orbital moments of atomic shells
    real(dp), intent(inout), allocatable :: orbitalL(:,:,:)

    !> imaginary part of density matrix
    real(dp), intent(inout), allocatable :: iRhoPrim(:,:)

    !> dense complex (k-points) hamiltonian storage
    complex(dp), intent(inout), allocatable :: HSqrCplx(:,:)

    !> dense complex (k-points) overlap storage
    complex(dp), intent(inout), allocatable :: SSqrCplx(:,:)

    complex(dp), allocatable :: rhoSqrCplx(:,:)
    real(dp), allocatable :: rVecTemp(:), orbitalLPart(:,:,:)
    integer :: iKS, iK, iS
    integer :: nAtom
    logical :: tImHam

    ! as each spin and k-point combination forms a separate group for this solver, then iKS = 1
    iKS = 1
    iK = parallelKS%localKS(1, iKS)
    iS = parallelKS%localKS(2, iKS)

    nAtom = size(nNeighbourSK)

    if (tSpinOrbit .and. .not. tDualSpinOrbit) then
      energy%atomLS(:) = 0.0_dp
      allocate(rVecTemp(nAtom))
    end if

    if (tSpinOrbit .and. .not. tDualSpinOrbit) then
      allocate(orbitalLPart(3, orb%mShell, nAtom))
      orbitalL(:,:,:) = 0.0_dp
    end if

    if (allocated(iHam)) then
      call unpackHPauliBlacs(env%blacs, ham, kPoint(:,iK), neighbourList%iNeighbour,&
          & nNeighbourSK, iCellVec, cellVec, iSparseStart, img2CentCell, orb%mOrb, denseDesc,&
          & HSqrCplx, iorig=iHam)
    else
      call unpackHPauliBlacs(env%blacs, ham, kPoint(:,iK), neighbourList%iNeighbour,&
          & nNeighbourSK, iCellVec, cellVec, iSparseStart, img2CentCell, orb%mOrb, denseDesc,&
          & HSqrCplx)
    end if
    if (.not. elsiSolver%tCholeskyDecomposed) then
      SSqrCplx(:,:) = 0.0_dp
      call unpackSPauliBlacs(env%blacs, over, kPoint(:,iK), neighbourList%iNeighbour,&
          & nNeighbourSK, iCellVec, cellVec, iSparseStart, img2CentCell, orb%mOrb, denseDesc,&
          & SSqrCplx)
      if (elsiSolver%OMM_Choleskii .or. elsiSolver%iSolver == electronicSolverTypes%pexsi) then
        elsiSolver%tCholeskyDecomposed = .true.
      end if
    end if
    if (allocated(xi) .and. .not. allocated(iHam)) then
      call addOnsiteSpinOrbitHam(env, xi, species, orb, denseDesc, HSqrCplx)
    end if
    allocate(rhoSqrCplx(size(HSqrCplx,dim=1),size(HSqrCplx,dim=2)))
    rhoSqrCplx(:,:) = 0.0_dp
    if (elsiSolver%tWriteHS) then
      call elsi_write_mat_complex(elsiSolver%rwHandle, "ELSI_Hcmplex.bin",&
          & HSqrCplx)
      call elsi_write_mat_complex(elsiSolver%rwHandle, "ELSI_Scmplx.bin",&
          & SSqrCplx)
      call elsi_finalize_rw(elsiSolver%rwHandle)
      call cleanShutdown("Finished dense matrix write")
    end if
    call elsi_dm_complex(elsiSolver%handle, HSqrCplx, SSqrCplx, rhoSqrCplx,&
        & Eband(iS))
    if (tSpinOrbit .and. .not. tDualSpinOrbit) then
      call getOnsiteSpinOrbitEnergy(env, rVecTemp, rhoSqrCplx, denseDesc, xi, orb, species)
      energy%atomLS = energy%atomLS + kWeight(iK) * rVecTemp
      if (tMulliken) then
        orbitalLPart(:,:,:) = 0.0_dp
        call getLOnsite(env, orbitalLPart, rhoSqrCplx, denseDesc, orb, species)
        orbitalL(:,:,:) = orbitalL + kWeight(iK) * orbitalLPart
      end if
    end if
    if (tImHam) then
      call packRhoPauliBlacs(env%blacs, denseDesc, rhoSqrCplx, kPoint(:,iK), kWeight(iK),&
          & neighbourList%iNeighbour, nNeighbourSK, orb%mOrb, iCellVec, cellVec,&
          & iSparseStart, img2CentCell, rhoPrim, iRhoPrim)
      iRhoPrim(:,:) = 2.0_dp * iRhoPrim
    else
      call packRhoPauliBlacs(env%blacs, denseDesc, rhoSqrCplx, kPoint(:,iK), kWeight(iK),&
          & neighbourList%iNeighbour, nNeighbourSK, orb%mOrb, iCellVec, cellVec,&
          & iSparseStart, img2CentCell, rhoPrim)
    end if
    RhoPrim(:,:) = 2.0_dp * RhoPrim

  end subroutine getDensityFromElsiPauliDense


  subroutine initSparse2Sparse(self, env, elsiSolver, iNeighbour, nNeighbourSK, iAtomStart,&
      & iSparseStart, img2CentCell)

    !> Sparse conversion instance
    type(TSparse2Sparse), intent(inout) :: self

    !> Environment settings
    type(TEnvironment), intent(in) :: env

    !> Electronic solver information
    type(TElsiSolver), intent(inout) :: elsiSolver

    !> Neighbour list for the atoms (First index from 0!)
    integer, intent(in) :: iNeighbour(0:,:)

    !> Nr. of neighbours for the atoms.
    integer, intent(in) :: nNeighbourSK(:)

    !> Atom offset for the squared matrix
    integer, intent(in) :: iAtomStart(:)

    !> indexing array for the sparse Hamiltonian
    integer, intent(in) :: iSparseStart(0:,:)

    !> Mapping between image atoms and corresponding atom in the central cell.
    integer, intent(in) :: img2CentCell(:)

    integer :: numCol, iAtom1, iAtom2, iAtom2f, nAtom, ii, jj, iNeigh, iNext, nOrb1, nOrb2
    integer :: iOrig
    integer, allocatable :: blockList(:,:)
    logical, allocatable :: tRowTrans(:)

    if (self%tInit) then
      ! this has been called before, clean out the storage
      deallocate(self%colPtrLocal)
      deallocate(self%rowIndLocal)
      deallocate(self%blockRow)
      deallocate(self%atomsInColumns)
    end if

    numCol = elsiSolver%n_basis

    self%colStartLocal = env%mpi%groupComm%rank * elsiSolver%CSR_blockSize + 1
    if (env%mpi%groupComm%rank /= env%mpi%groupComm%size - 1) then
      self%numColLocal = elsiSolver%CSR_blockSize
    else
      self%numColLocal = numCol - (env%mpi%groupComm%size - 1) * elsiSolver%CSR_blockSize
    end if
    self%colEndLocal = self%colStartLocal + self%numColLocal - 1

    allocate(self%colPtrLocal(self%numColLocal + 1))

    call pack2colptr_parallel(iNeighbour, nNeighbourSK, iAtomStart, iSparseStart, img2CentCell,&
        & self%colStartLocal, self%colEndLocal, self%nnzLocal, self%colPtrLocal)

    self%nnzGlobal = 0
    call mpifx_allreduce(env%mpi%groupComm, self%nnzLocal, self%nnzGlobal, MPI_SUM)

    nAtom = size(nNeighbourSK)

    allocate(self%blockRow(0:size(iNeighbour,dim=1)-1,nAtom))
    self%blockRow(:,:) = 0

    allocate(self%atomsInColumns(nAtom))
    self%atomsInColumns(nAtom) = 0
    self%nAtomsInColumns = 0

    allocate(blockList(nAtom,2))
    allocate(tRowTrans(nAtom))
    blockList = 0

    ! Offset in column belonging to transposed triangle
    blockList(:,2) = 1

    do iAtom1 = 1, nAtom
      ii = iAtomStart(iAtom1)
      nOrb1 = iAtomStart(iAtom1+1) - ii
      ! Offset in current column
      blockList(:,1) = 0
      tRowTrans = .false.
      ! Starting index for column in DFTB+ sparse structure, because column probaly already contains
      ! blocks coming from transposing previously processed elements.
      iNext = blockList(iAtom1, 2)
      do iNeigh = 0, nNeighbourSK(iAtom1)
        iOrig = iSparseStart(iNeigh,iAtom1) + 1
        iAtom2 = iNeighbour(iNeigh, iAtom1)
        iAtom2f = img2CentCell(iAtom2)

        jj = iAtomStart(iAtom2f)
        nOrb2 = iAtomStart(iAtom2f+1) - jj

        if (blockList(iAtom2f,1) == 0) then
          blockList(iAtom2f,1) = iNext
          iNext = iNext + nOrb2
        end if

        self%blockRow(iNeigh, iAtom1) = blockList(iAtom2f,1)

        if ( .not. isBlockInLocal(jj,jj+nOrb2-1, ii, ii+nOrb1-1, self%colStartLocal,&
            & self%colEndLocal)) then
          cycle
        end if

        if (iAtom1 /= self%atomsInColumns(max(self%nAtomsInColumns,1))) then
          ! this atom is required for the local columns
          self%nAtomsInColumns = self%nAtomsInColumns + 1
          self%atomsInColumns(self%nAtomsInColumns) = iAtom1
        end if

        if (ii == jj) then
          ! on the diagonal, can not be in other triangle
          cycle
        end if

        ! Because of folding of periodic images, it can happen that the transposed block has already
        ! been registered.
        if (.not. tRowTrans(iAtom2f)) then
          blockList(iAtom2f,2) = blockList(iAtom2f,2) + nOrb1
          tRowTrans(iAtom2f) = .true.
        end if

      end do
    end do

    self%tInit = .true.

  end subroutine initSparse2Sparse


  !> Calculates density matrix using the elsi routine.
  subroutine calcDensityRealElsi(self, parallelKS, elsiSolver, ham, over,&
      & iNeighbour, nNeighbourSK, iAtomStart, iSparseStart, img2CentCell, orb, rho, Eband)

    !> Sparse conversion instance
    type(TSparse2Sparse), intent(inout) :: self

    !> Contains (iK, iS) tuples to be processed in parallel by various processor groups
    type(TParallelKS), intent(in) :: parallelKS

    !> Electronic solver information
    type(TElsiSolver), intent(inout) :: elsiSolver

    !> hamiltonian in sparse storage
    real(dp), intent(in) :: ham(:,:)

    !> overlap matrix in sparse storage
    real(dp), intent(in) :: over(:)

    !> Neighbour list for the atoms (First index from 0!)
    integer, intent(in) :: iNeighbour(0:,:)

    !> Nr. of neighbours for the atoms.
    integer, intent(in) :: nNeighbourSK(:)

    !> Atom offset for the squared matrix
    integer, intent(in) :: iAtomStart(:)

    !> indexing array for the sparse Hamiltonian
    integer, intent(in) :: iSparseStart(0:,:)

    !> Mapping between image atoms and corresponding atom in the central cell.
    integer, intent(in) :: img2CentCell(:)

    !> data structure with atomic orbital information
    type(TOrbitals), intent(in) :: orb

    !> Density matrix in DFTB+ sparse format
    real(dp), intent(out) :: rho(:,:)

    !> Band energy
    real(dp), intent(out) :: Eband(:)

    integer :: iS

    real(dp), allocatable :: HnzValLocal(:), SnzValLocal(:)
    real(dp), allocatable :: DMnzValLocal(:)
    logical :: tFirstCall

    tFirstCall = .false.
    if (.not. allocated(self%rowIndLocal)) then
      allocate(self%rowIndLocal(self%nnzLocal))
      tFirstCall = .true.
    end if

    allocate(HnzValLocal(self%nnzLocal))
    allocate(SnzValLocal(self%nnzLocal))
    allocate(DMnzValLocal(self%nnzLocal))

    if (tFirstCall) then
      ! also generate rowIndLocal for the new structure
      call pack2elsi_real(self, over, iNeighbour, nNeighbourSK, iAtomStart, iSparseStart,&
          & img2CentCell, self%colStartLocal, self%colEndLocal, self%colPtrLocal, SnzValLocal,&
          & self%rowIndLocal)
      call elsi_set_csc(elsiSolver%handle, self%nnzGlobal, self%nnzLocal,&
          & self%numColLocal, self%rowIndLocal, self%colPtrLocal)
    else
      call pack2elsi_real(self, over, iNeighbour, nNeighbourSK, iAtomStart, iSparseStart,&
          & img2CentCell, self%colStartLocal, self%colEndLocal, self%colPtrLocal, SnzValLocal)
    end if

    iS = parallelKS%localKS(2, 1)
    call pack2elsi_real(self, ham(:,iS), iNeighbour, nNeighbourSK, iAtomStart,&
        & iSparseStart, img2CentCell, self%colStartLocal, self%colEndLocal, self%colPtrLocal,&
        & HnzValLocal)

    if (elsiSolver%tWriteHS) then
      ! set up for sparse matrix writing
      call elsi_set_rw_csc(elsiSolver%rwHandle, self%nnzGlobal, self%nnzLocal,&
          & self%numColLocal)

      call elsi_write_mat_real_sparse(elsiSolver%rwHandle, "ELSI_HrealSparse.bin",&
          & self%rowIndLocal, self%colPtrLocal, HnzValLocal)
      call elsi_write_mat_real_sparse(elsiSolver%rwHandle, "ELSI_SrealSparse.bin",&
          & self%rowIndLocal, self%colPtrLocal, SnzValLocal)
      call elsi_finalize_rw(elsiSolver%rwHandle)
      call cleanShutdown("Finished matrix write")
    end if

    ! Load the matrix into ELSI and solve DM
    call elsi_dm_real_sparse(elsiSolver%handle, HnzValLocal, SnzValLocal, DMnzValLocal,&
        & Eband(iS))

    ! get DM back into DFTB+ format
    rho(:,:) = 0.0_dp
    call elsi2pack_real(self, self%colStartLocal, self%colEndLocal, iNeighbour,&
        & nNeighbourSK, orb%mOrb, iAtomStart, iSparseStart, img2CentCell, self%colPtrLocal,&
        & DMnzValLocal, self%blockRow, rho(:,iS))

  end subroutine calcDensityRealElsi


  subroutine initPexsiDeltaVRanges(tSccCalc, potential, elsiSolver)

    !> Whether we have an SCC calculation
    logical, intent(in) :: tSccCalc

    !> Potentials acting
    type(TPotentials), intent(in) :: potential

    !> Electronic solver information
    type(TElsiSolver), intent(inout) :: elsiSolver

    if (tSccCalc) then
      if (.not.allocated(elsiSolver%PEXSI_VOld)) then
        allocate(elsiSolver%PEXSI_VOld(size(potential%intBlock(:,:,:,1))))
      end if
      elsiSolver%PEXSI_VOld = reshape(potential%intBlock(:,:,:,1) + potential%extBlock(:,:,:,1),&
          & [size(potential%extBlock(:,:,:,1))])
    end if
    elsiSolver%PEXSI_DeltaVmin = 0.0_dp
    elsiSolver%PEXSI_DeltaVmax = 0.0_dp

  end subroutine initPexsiDeltaVRanges


  subroutine updatePexsiDeltaVRanges(potential, elsiSolver)

    !> Potentials acting
    type(TPotentials), intent(in) :: potential

    !> Electronic solver information
    type(TElsiSolver), intent(inout) :: elsiSolver


    elsiSolver%PEXSI_DeltaVmin = minval(elsiSolver%PEXSI_VOld&
        & - reshape(potential%intBlock(:,:,:,1) + potential%extBlock(:,:,:,1),&
        & [size(potential%extBlock(:,:,:,1))]))
    elsiSolver%PEXSI_DeltaVmax = maxval(elsiSolver%PEXSI_VOld&
        & - reshape(potential%intBlock(:,:,:,1) + potential%extBlock(:,:,:,1),&
        & [size(potential%extBlock(:,:,:,1))]))
    elsiSolver%PEXSI_VOld = reshape(potential%intBlock(:,:,:,1) + potential%extBlock(:,:,:,1),&
        & [size(potential%extBlock(:,:,:,1))])

  end subroutine updatePexsiDeltaVRanges


  subroutine getEDensityMtxElsiPauli(env, elsiSolver, denseDesc, kPoint, kWeight, neighbourList,&
      & nNeighbourSK, orb, iSparseStart, img2CentCell, iCellVec, cellVec, parallelKS, ERhoPrim,&
      & SSqrCplx)

    !> Environment settings
    type(TEnvironment), intent(inout) :: env

    !> Electronic solver information
    type(TElsiSolver), intent(inout) :: elsiSolver

    !> Dense matrix descriptor
    type(TDenseDescr), intent(in) :: denseDesc

    !> K-points
    real(dp), intent(in) :: kPoint(:,:)

    !> Weights for k-points
    real(dp), intent(in) :: kWeight(:)

    !> list of neighbours for each atom
    type(TNeighbourList), intent(in) :: neighbourList

    !> Number of neighbours for each of the atoms
    integer, intent(in) :: nNeighbourSK(:)

    !> Atomic orbital information
    type(TOrbitals), intent(in) :: orb

    !> Index array for the start of atomic blocks in sparse arrays
    integer, intent(in) :: iSparseStart(:,:)

    !> map from image atoms to the original unique atom
    integer, intent(in) :: img2CentCell(:)

    !> Index for which unit cell atoms are associated with
    integer, intent(in) :: iCellVec(:)

    !> Vectors (in units of the lattice constants) to cells of the lattice
    real(dp), intent(in) :: cellVec(:,:)

    !> K-points and spins to process
    type(TParallelKS), intent(in) :: parallelKS

    !> Energy weighted sparse matrix
    real(dp), intent(out) :: ERhoPrim(:)

    !> Storage for dense overlap matrix
    complex(dp), intent(inout), allocatable :: SSqrCplx(:,:)

    integer :: iK

    if (elsiSolver%CSR) then
      call error("EDensity via sparse ELSI solver not supported for Pauli-matrices")
    else
      ! iKS always 1, as number of groups matches the number of k-points
      iK = parallelKS%localKS(1, 1)
      call elsi_get_edm_complex(elsiSolver%handle, SSqrCplx)
      call packERhoPauliBlacs(env%blacs, denseDesc, SSqrCplx, kPoint(:,iK), kWeight(iK),&
          & neighbourList%iNeighbour, nNeighbourSK, orb%mOrb, iCellVec, cellVec, iSparseStart,&
          & img2CentCell, ERhoPrim)
      call mpifx_allreduceip(env%mpi%globalComm, ERhoPrim, MPI_SUM)
    end if

  end subroutine getEDensityMtxElsiPauli


  subroutine getEDensityMtxElsiReal(env, elsiSolver, sparseIndexing, denseDesc, neighbourList,&
      & nNeighbourSK, orb, iSparseStart, img2CentCell, ERhoPrim, SSqrReal)

    !> Environment settings
    type(TEnvironment), intent(inout) :: env

    !> Electronic solver information
    type(TElsiSolver), intent(inout) :: elsiSolver

    !> Structure for sparse data representions, if using ELSI
    type(TSparse2Sparse), intent(inout) :: sparseIndexing

    !> Dense matrix descriptor
    type(TDenseDescr), intent(in) :: denseDesc

    !> list of neighbours for each atom
    type(TNeighbourList), intent(in) :: neighbourList

    !> Number of neighbours for each of the atoms
    integer, intent(in) :: nNeighbourSK(:)

    !> Atomic orbital information
    type(TOrbitals), intent(in) :: orb

    !> Index array for the start of atomic blocks in sparse arrays
    integer, intent(in) :: iSparseStart(:,:)

    !> map from image atoms to the original unique atom
    integer, intent(in) :: img2CentCell(:)

    !> Energy weighted sparse matrix
    real(dp), intent(out) :: ERhoPrim(:)

    !> Storage for dense overlap matrix
    real(dp), intent(inout), allocatable :: SSqrReal(:,:)

    if (elsiSolver%CSR) then
      call getEDensityRealElsi(sparseIndexing, elsiSolver,&
          & neighbourList%iNeighbour, nNeighbourSK, denseDesc%iAtomStart, iSparseStart,&
          & img2CentCell, orb, erhoPrim)
    else
      call elsi_get_edm_real(elsiSolver%handle, SSqrReal)
      ERhoPrim(:) = 0.0_dp
      call packRhoRealBlacs(env%blacs, denseDesc, SSqrReal, neighbourList%iNeighbour, nNeighbourSK,&
          & orb%mOrb, iSparseStart, img2CentCell, ERhoPrim)
    end if

    ! add contributions from different spin channels together if necessary
    call mpifx_allreduceip(env%mpi%globalComm, ERhoPrim, MPI_SUM)

  end subroutine getEDensityMtxElsiReal


  subroutine getEDensityMtxElsiCmplx(env, elsiSolver, sparseIndexing, denseDesc, kPoint, kWeight,&
      & neighbourList, nNeighbourSK, orb, iSparseStart, img2CentCell, iCellVec, cellVec,&
      & parallelKS, ERhoPrim, SSqrCplx)

    !> Environment settings
    type(TEnvironment), intent(inout) :: env

    !> Electronic solver information
    type(TElsiSolver), intent(inout) :: elsiSolver

    !> Structure for sparse data representions, if using ELSI
    type(TSparse2Sparse), intent(inout) :: sparseIndexing

    !> Dense matrix descriptor
    type(TDenseDescr), intent(in) :: denseDesc

    !> K-points
    real(dp), intent(in) :: kPoint(:,:)

    !> Weights for k-points
    real(dp), intent(in) :: kWeight(:)

    !> list of neighbours for each atom
    type(TNeighbourList), intent(in) :: neighbourList

    !> Number of neighbours for each of the atoms
    integer, intent(in) :: nNeighbourSK(:)

    !> Atomic orbital information
    type(TOrbitals), intent(in) :: orb

    !> Index array for the start of atomic blocks in sparse arrays
    integer, intent(in) :: iSparseStart(:,:)

    !> map from image atoms to the original unique atom
    integer, intent(in) :: img2CentCell(:)

    !> Index for which unit cell atoms are associated with
    integer, intent(in) :: iCellVec(:)

    !> Vectors (in units of the lattice constants) to cells of the lattice
    real(dp), intent(in) :: cellVec(:,:)

    !> K-points and spins to process
    type(TParallelKS), intent(in) :: parallelKS

    !> Energy weighted sparse matrix
    real(dp), intent(out) :: ERhoPrim(:)

    !> Storage for dense overlap matrix (complex case)
    complex(dp), intent(inout), allocatable :: SSqrCplx(:,:)

    integer :: iK

    ERhoPrim(:) = 0.0_dp
    ! iKS always 1
    iK = parallelKS%localKS(1, 1)
    if (elsiSolver%CSR) then
      call getEDensityComplexElsi(sparseIndexing, elsiSolver, kPoint(:,iK), kWeight(iK), iCellVec,&
          & cellVec, neighbourList%iNeighbour, nNeighbourSK, denseDesc%iAtomStart, iSparseStart,&
          & img2CentCell, orb, erhoPrim)
    else
      call elsi_get_edm_complex(elsiSolver%handle, SSqrCplx)
      call packRhoCplxBlacs(env%blacs, denseDesc, SSqrCplx, kPoint(:,iK), kWeight(iK),&
          & neighbourList%iNeighbour, nNeighbourSK, orb%mOrb, iCellVec, cellVec, iSparseStart,&
          & img2CentCell, ERhoPrim)
    end if
    call mpifx_allreduceip(env%mpi%globalComm, ERhoPrim, MPI_SUM)

  end subroutine getEDensityMtxElsiCmplx


  !> Gets energy density matrix using the elsi routine.
  subroutine getEDensityRealElsi(self, elsiSolver, iNeighbour,&
      & nNeighbourSK, iAtomStart, iSparseStart, img2CentCell, orb, Erho)

    !> Sparse conversion instance
    type(TSparse2Sparse), intent(inout) :: self

    !> Electronic solver information
    type(TElsiSolver), intent(inout) :: elsiSolver

    !> Neighbour list for the atoms (First index from 0!)
    integer, intent(in) :: iNeighbour(0:,:)

    !> Nr. of neighbours for the atoms.
    integer, intent(in) :: nNeighbourSK(:)

    !> Atom offset for the squared matrix
    integer, intent(in) :: iAtomStart(:)

    !> indexing array for the sparse Hamiltonian
    integer, intent(in) :: iSparseStart(0:,:)

    !> Mapping between image atoms and corresponding atom in the central cell.
    integer, intent(in) :: img2CentCell(:)

    !> data structure with atomic orbital information
    type(TOrbitals), intent(in) :: orb

    !> Density matrix in DFTB+ sparse format
    real(dp), intent(out) :: Erho(:)

    integer :: iS
    real(dp), allocatable :: EDMnzValLocal(:)

    allocate(EDMnzValLocal(self%nnzLocal))

    ! get the energy weighted density matrix from ELSI
    call elsi_get_edm_real_sparse(elsiSolver%handle, EDMnzValLocal)

    ! get EDM back into DFTB+ format
    erho(:) = 0.0_dp
    call elsi2pack_real(self, self%colStartLocal, self%colEndLocal, iNeighbour,&
        & nNeighbourSK, orb%mOrb, iAtomStart, iSparseStart, img2CentCell, self%colPtrLocal,&
        & EDMnzValLocal, self%blockRow, erho)

  end subroutine getEDensityRealElsi


  !> Calculates density matrix using the elsi routine.
  subroutine calcDensityComplexElsi(self, parallelKS, elsiSolver, kPoint, kWeight, iCellVec,&
      & cellVec, ham, over, iNeighbour, nNeighbourSK, iAtomStart, iSparseStart, img2CentCell, orb,&
      & rho, Eband)

    !> Sparse conversion instance
    type(TSparse2Sparse), intent(inout) :: self

    !> Contains (iK, iS) tuples to be processed in parallel by various processor groups
    type(TParallelKS), intent(in) :: parallelKS

    !> Electronic solver information
    type(TElsiSolver), intent(inout) :: elsiSolver

    !> Current k-point
    real(dp), intent(in) :: kPoint(:)

    !> Weight for current k-points
    real(dp), intent(in) :: kWeight

    !> Index for which unit cell atoms are associated with
    integer, intent(in) :: iCellVec(:)

    !> Vectors (in units of the lattice constants) to cells of the lattice
    real(dp), intent(in) :: cellVec(:,:)

    !> hamiltonian in sparse storage
    real(dp), intent(in) :: ham(:,:)

    !> overlap matrix in sparse storage
    real(dp), intent(in) :: over(:)

    !> Neighbour list for the atoms (First index from 0!)
    integer, intent(in) :: iNeighbour(0:,:)

    !> Nr. of neighbours for the atoms.
    integer, intent(in) :: nNeighbourSK(:)

    !> Atom offset for the squared matrix
    integer, intent(in) :: iAtomStart(:)

    !> indexing array for the sparse Hamiltonian
    integer, intent(in) :: iSparseStart(0:,:)

    !> Mapping between image atoms and corresponding atom in the central cell.
    integer, intent(in) :: img2CentCell(:)

    !> data structure with atomic orbital information
    type(TOrbitals), intent(in) :: orb

    !> Density matrix in DFTB+ sparse format
    real(dp), intent(out) :: rho(:,:)

    !> Band energy
    real(dp), intent(out) :: Eband(:)

    integer :: iS

    complex(dp), allocatable :: HnzValLocal(:), SnzValLocal(:)
    complex(dp), allocatable :: DMnzValLocal(:)
    logical :: tFirstCall

    tFirstCall = .false.
    if (.not. allocated(self%rowIndLocal)) then
      allocate(self%rowIndLocal(self%nnzLocal))
      tFirstCall = .true.
    end if

    allocate(HnzValLocal(self%nnzLocal))
    allocate(SnzValLocal(self%nnzLocal))
    allocate(DMnzValLocal(self%nnzLocal))

    if (tFirstCall) then
      ! also generate rowIndLocal for the new structure
      call pack2elsi_cmplx(self, over, iNeighbour, nNeighbourSK, iAtomStart, iSparseStart,&
          & img2CentCell, kPoint, kWeight, iCellVec, cellVec, self%colStartLocal, self%colEndLocal,&
          & self%colPtrLocal, SnzValLocal, self%rowIndLocal)
      call elsi_set_csc(elsiSolver%handle, self%nnzGlobal, self%nnzLocal,&
          & self%numColLocal, self%rowIndLocal, self%colPtrLocal)
    else
      call pack2elsi_cmplx(self, over, iNeighbour, nNeighbourSK, iAtomStart, iSparseStart,&
          & img2CentCell, kPoint, kWeight, iCellVec, cellVec, self%colStartLocal, self%colEndLocal,&
          & self%colPtrLocal, SnzValLocal)
    end if

    iS = parallelKS%localKS(2, 1)
    call pack2elsi_cmplx(self, ham(:,iS), iNeighbour, nNeighbourSK, iAtomStart,&
        & iSparseStart, img2CentCell, kPoint, kWeight, iCellVec, cellVec, self%colStartLocal,&
        & self%colEndLocal, self%colPtrLocal, HnzValLocal)

    if (elsiSolver%tWriteHS) then
      ! set up for sparse matrix writing
      call elsi_set_rw_csc(elsiSolver%rwHandle, self%nnzGlobal, self%nnzLocal,&
          & self%numColLocal)

      call elsi_write_mat_complex_sparse(elsiSolver%rwHandle, "ELSI_HcmplxSparse.bin",&
          & self%rowIndLocal, self%colPtrLocal, HnzValLocal)
      call elsi_write_mat_complex_sparse(elsiSolver%rwHandle, "ELSI_ScmplxSparse.bin",&
          & self%rowIndLocal, self%colPtrLocal, SnzValLocal)
      call elsi_finalize_rw(elsiSolver%rwHandle)
      call cleanShutdown("Finished matrix write")
    end if

    ! Load the matrix into ELSI and solve DM
    call elsi_dm_complex_sparse(elsiSolver%handle, HnzValLocal, SnzValLocal,&
        & DMnzValLocal, Eband(iS))

    ! get DM back into DFTB+ format
    rho(:,:) = 0.0_dp
    call elsi2pack_cmplx(self, self%colStartLocal, self%colEndLocal, iNeighbour,&
        & nNeighbourSK, orb%mOrb, iAtomStart, iSparseStart, img2CentCell, kPoint, kWeight,&
        & iCellVec, cellVec, self%colPtrLocal, DMnzValLocal, self%blockRow, rho(:,iS))

  end subroutine calcDensityComplexElsi


  !> Gets energy density matrix using the elsi routine.
  subroutine getEDensityComplexElsi(self, elsiSolver, kPoint, kWeight, iCellVec, cellVec,&
      & iNeighbour, nNeighbourSK, iAtomStart, iSparseStart, img2CentCell, orb, erho)

    !> Sparse conversion instance
    type(TSparse2Sparse), intent(inout) :: self

    !> Electronic solver information
    type(TElsiSolver), intent(inout) :: elsiSolver

    !> Current k-point
    real(dp), intent(in) :: kPoint(:)

    !> Weight for current k-points
    real(dp), intent(in) :: kWeight

    !> Index for which unit cell atoms are associated with
    integer, intent(in) :: iCellVec(:)

    !> Vectors (in units of the lattice constants) to cells of the lattice
    real(dp), intent(in) :: cellVec(:,:)

    !> Neighbour list for the atoms (First index from 0!)
    integer, intent(in) :: iNeighbour(0:,:)

    !> Nr. of neighbours for the atoms.
    integer, intent(in) :: nNeighbourSK(:)

    !> Atom offset for the squared matrix
    integer, intent(in) :: iAtomStart(:)

    !> indexing array for the sparse Hamiltonian
    integer, intent(in) :: iSparseStart(0:,:)

    !> Mapping between image atoms and corresponding atom in the central cell.
    integer, intent(in) :: img2CentCell(:)

    !> data structure with atomic orbital information
    type(TOrbitals), intent(in) :: orb

    !> Density matrix in DFTB+ sparse format
    real(dp), intent(out) :: erho(:)

    integer :: iS
    complex(dp), allocatable :: EDMnzValLocal(:)

    allocate(EDMnzValLocal(self%nnzLocal))

    ! get the energy weighted density matrix from ELSI
    call elsi_get_edm_complex_sparse(elsiSolver%handle, EDMnzValLocal)

    ! get EDM back into DFTB+ format
    erho(:) = 0.0_dp
    call elsi2pack_cmplx(self, self%colStartLocal, self%colEndLocal, iNeighbour,&
        & nNeighbourSK, orb%mOrb, iAtomStart, iSparseStart, img2CentCell, kPoint, kWeight,&
        & iCellVec, cellVec, self%colPtrLocal, EDMnzValLocal, self%blockRow, erho)

  end subroutine getEDensityComplexElsi


  !> Creating colptr and nnz for CSC matrix format from packed format
  subroutine pack2colptr_parallel(iNeighbour, nNeighbourSK, iAtomStart, iSparseStart,&
      & img2CentCell, colStartLocal, colEndLocal, nnzLocal, colPtrLocal)

    !> Neighbour list for each atom (first index from 0!).
    integer, intent(in) :: iNeighbour(0:,:)

    !> Nr. of neighbours for each atom (incl. itself).
    integer, intent(in) :: nNeighbourSK(:)

    !> Atom offset for the squared Hamiltonian
    integer, intent(in) :: iAtomStart(:)

    !> Indexing array for the sparse Hamiltonian
    integer, intent(in) :: iSparseStart(0:,:)

    !> Atomic mapping indexes
    integer, intent(in) :: img2CentCell(:)

    !> Column of global matrix where local matrix starts
    integer, intent(in) :: colStartLocal

    !> Column of global matrix where local matrix ends
    integer, intent(in) :: colEndLocal

    !> Local number of  non-zero elements
    integer, intent(out) :: nnzLocal

    !> Local column pointer
    integer, intent(out) :: colPtrLocal(:)

    integer :: nAtom
    integer :: iOrig, ii, jj, nOrb1, nOrb2
    integer :: iNeigh
    integer :: iAtom1, iAtom2, iAtom2f, kk

    logical, allocatable :: blockList(:)

    nAtom = size(iNeighbour, dim=2)

    colPtrLocal(:) = 0
    nnzLocal = 0

    allocate(blockList(nAtom))

    ! Loop over all atom blocks
    do iAtom1 = 1, nAtom
      ii = iAtomStart(iAtom1)
      nOrb1 = iAtomStart(iAtom1+1) - ii
      blockList(:) = .false.
      do iNeigh = 0, nNeighbourSK(iAtom1)
        iOrig = iSparseStart(iNeigh,iAtom1) + 1
        iAtom2 = iNeighbour(iNeigh, iAtom1)
        iAtom2f = img2CentCell(iAtom2)
        jj = iAtomStart(iAtom2f)
        nOrb2 = iAtomStart(iAtom2f+1) - jj

        if (blockList(iAtom2f)) then
          cycle
        else
          blockList(iAtom2f) = .true.
        end if

        if (isBlockInLocal(jj, jj+nOrb2-1, ii, ii+nOrb1-1, colStartLocal, colEndLocal)) then
          do kk = ii, ii + nOrb1 - 1
            if (isColumnInLocal(kk, colStartLocal, colEndLocal)) then
              nnzLocal = nnzLocal + nOrb2
              colPtrLocal(kk - colStartLocal + 2) = colPtrLocal(kk - colStartLocal + 2) + nOrb2
            end if
          end do
        end if

        if (ii == jj) then
          ! on the diagonal, can not be in other triangle
          cycle
        end if

        do kk = jj, jj + nOrb2 - 1
          if (isColumnInLocal(kk, colStartLocal, colEndLocal)) then
            nnzLocal = nnzLocal + nOrb1
            colPtrLocal(kk -colStartLocal + 2) = colPtrLocal(kk - colStartLocal + 2) + nOrb1
          end if
        end do
      end do
    end do
    colPtrLocal(1) = 1
    do ii = 2, size(colPtrLocal)
      colPtrLocal(ii) = colPtrLocal(ii-1) + colPtrLocal(ii)
    end do
    colPtrLocal(size(colPtrLocal)) = nnzLocal + 1

  end subroutine pack2colptr_parallel


  !> Convert sparse DFTB+ matrix to distributed CSC matrix format for ELSI calculations for real
  !> matrices
  !>
  !> NOTE: ELSI needs the full matrix (both triangles)
  subroutine pack2elsi_real(self, orig, iNeighbour, nNeighbourSK, iAtomStart, iSparseStart,&
      & img2CentCell, colStartLocal, colEndLocal, colPtrLocal, nzValLocal, rowIndLocal)

    !> Sparse conversion instance
    type(TSparse2Sparse), intent(in) :: self

    !> Sparse Hamiltonian
    real(dp), intent(in) :: orig(:)

    !> Neighbour list for each atom (first index from 0!).
    integer, intent(in) :: iNeighbour(0:,:)

    !> Nr. of neighbours for each atom (incl. itself).
    integer, intent(in) :: nNeighbourSK(:)

    !> Atom offset for the squared Hamiltonian
    integer, intent(in) :: iAtomStart(:)

    !> Indexing array for the sparse Hamiltonian
    integer, intent(in) :: iSparseStart(0:,:)

    !> Atomic mapping indexes.
    integer, intent(in) :: img2CentCell(:)

    !> Column of global matrix where local matrix starts
    integer, intent(in) :: colStartLocal

    !> Column of global matrix where local matrix ends
    integer, intent(in) :: colEndLocal

    !> Local column pointer
    integer, intent(in) :: colPtrLocal(:)

    !> Local non-zero elements
    real(dp), intent(out) :: nzValLocal(:)

    !> Local row index pointer
    integer, intent(out), optional :: rowIndLocal(:)

    integer :: nAtom
    integer :: iOrig, ii, jj, nOrb1, nOrb2
    integer :: iNeigh
    integer :: iAt, iAtom1, iAtom2, iAtom2f

    integer :: iNext
    integer, allocatable :: blockList(:,:)
    logical, allocatable :: tRowTrans(:)

    nAtom = size(iNeighbour, dim=2)

    @:ASSERT(nAtom > 0)
    @:ASSERT(size(nNeighbourSK) == nAtom)
    @:ASSERT(size(iAtomStart) == nAtom + 1)

    allocate(blockList(nAtom,2))
    allocate(tRowTrans(nAtom))
    if (present(rowIndLocal)) then
      rowIndLocal(:) = 0
    end if
    nzValLocal(:) = 0.0_dp

    ! Offset in column belonging to transposed triangle
    blockList(:,2) = 1

    ! loop over atoms relevant to this processor
    do iAt = 1, self%nAtomsInColumns
      iAtom1 = self%atomsInColumns(iAt)
      ii = iAtomStart(iAtom1)
      nOrb1 = iAtomStart(iAtom1+1) - ii
      ! Offset in current column
      blockList(:,1) = 0
      tRowTrans = .false.
      ! Starting index for column in DFTB+ sparse structure, because column probaly already contains
      ! blocks coming from transposing previously processed elements.
      iNext = blockList(iAtom1, 2)
      do iNeigh = 0, nNeighbourSK(iAtom1)
        iOrig = iSparseStart(iNeigh,iAtom1) + 1
        iAtom2 = iNeighbour(iNeigh, iAtom1)
        iAtom2f = img2CentCell(iAtom2)

        jj = iAtomStart(iAtom2f)
        @:ASSERT(jj >= ii)
        nOrb2 = iAtomStart(iAtom2f+1) - jj

        if (blockList(iAtom2f,1) == 0) then
          blockList(iAtom2f,1) = iNext
          iNext = iNext + nOrb2
        end if

        if ( .not. isBlockInLocal(jj,jj+nOrb2-1, ii, ii+nOrb1-1, colStartLocal, colEndLocal)) then
          cycle
        end if
        call addBlock2Elsi(reshape(orig(iOrig:iOrig+nOrb1*nOrb2-1), [ nOrb2, nOrb1 ]),&
            & colStartLocal, colEndLocal, jj, ii, colPtrLocal, blockList(iAtom2f,1)-1, nzValLocal,&
            & rowIndLocal)

        if (ii == jj) then
          cycle
        end if

        ! Because of folding of periodic images, it can happen that the transposed block has already
        ! been registered.
        if (.not. tRowTrans(iAtom2f)) then
          blockList(iAtom2f,2) = blockList(iAtom2f,2) + nOrb1
          tRowTrans(iAtom2f) = .true.
        end if

        call addBlock2Elsi(transpose(reshape(orig(iOrig:iOrig+nOrb1*nOrb2-1), [ nOrb2, nOrb1 ])),&
            & colStartLocal, colEndLocal, ii, jj, colPtrLocal, blockList(iAtom2f,2)-nOrb1-1,&
            & nzValLocal, rowIndLocal)
      end do
    end do

  end subroutine pack2elsi_real


  !> Convert sparse DFTB+ matrix to distributed CSC matrix format for ELSI calculations for complex
  !> matrices
  !>
  !> NOTE: ELSI needs the full matrix (both triangles)
  subroutine pack2elsi_cmplx(self, orig, iNeighbour, nNeighbourSK, iAtomStart,&
      & iSparseStart, img2CentCell, kPoint, kWeight, iCellVec, cellVec, colStartLocal, colEndLocal,&
      & colPtrLocal, nzValLocal, rowIndLocal)

    !> Sparse conversion instance
    type(TSparse2Sparse), intent(in) :: self

    !> Sparse Hamiltonian
    real(dp), intent(in) :: orig(:)

    !> Neighbour list for each atom (first index from 0!).
    integer, intent(in) :: iNeighbour(0:,:)

    !> Nr. of neighbours for each atom (incl. itself).
    integer, intent(in) :: nNeighbourSK(:)

    !> Atom offset for the squared Hamiltonian
    integer, intent(in) :: iAtomStart(:)

    !> Indexing array for the sparse Hamiltonian
    integer, intent(in) :: iSparseStart(0:,:)

    !> Atomic mapping indexes.
    integer, intent(in) :: img2CentCell(:)

    !> Current k-point
    real(dp), intent(in) :: kPoint(:)

    !> Weight for current k-points
    real(dp), intent(in) :: kWeight

    !> Index for which unit cell atoms are associated with
    integer, intent(in) :: iCellVec(:)

    !> Vectors (in units of the lattice constants) to cells of the lattice
    real(dp), intent(in) :: cellVec(:,:)

    !> Column of global matrix where local matrix starts
    integer, intent(in) :: colStartLocal

    !> Column of global matrix where local matrix ends
    integer, intent(in) :: colEndLocal

    !> Local column pointer
    integer, intent(in) :: colPtrLocal(:)

    !> Local non-zero elements
    complex(dp), intent(out) :: nzValLocal(:)

    !> Local row index pointer
    integer, intent(out), optional :: rowIndLocal(:)

    integer :: nAtom
    integer :: iOrig, ii, jj, nOrb1, nOrb2
    integer :: iNeigh
    integer :: iAt, iAtom1, iAtom2, iAtom2f

    integer :: iNext, iVec
    integer, allocatable :: blockList(:,:)
    logical, allocatable :: tRowTrans(:)
    real(dp) :: kPoint2p(3)
    complex(dp) :: phase

    kPoint2p(:) = 2.0_dp * pi * kPoint(:)
    nAtom = size(iNeighbour, dim=2)

    @:ASSERT(nAtom > 0)
    @:ASSERT(size(nNeighbourSK) == nAtom)
    @:ASSERT(size(iAtomStart) == nAtom + 1)

    allocate(blockList(nAtom,2))
    allocate(tRowTrans(nAtom))
    if (present(rowIndLocal)) then
      rowIndLocal(:) = 0
    end if
    nzValLocal(:) = cmplx(0,0,dp)

    ! Offset in column belonging to transposed triangle
    blockList(:,2) = 1

    ! loop over atoms relevant to this processor
    do iAt = 1, self%nAtomsInColumns
      iAtom1 = self%atomsInColumns(iAt)
      ii = iAtomStart(iAtom1)
      nOrb1 = iAtomStart(iAtom1+1) - ii
      ! Offset in current column
      blockList(:,1) = 0
      tRowTrans = .false.
      ! Starting index for column in DFTB+ sparse structure, because column probaly already contains
      ! blocks coming from transposing previously processed elements.
      iNext = blockList(iAtom1, 2)
      do iNeigh = 0, nNeighbourSK(iAtom1)
        iOrig = iSparseStart(iNeigh,iAtom1) + 1
        iAtom2 = iNeighbour(iNeigh, iAtom1)
        iAtom2f = img2CentCell(iAtom2)

        jj = iAtomStart(iAtom2f)
        @:ASSERT(jj >= ii)
        nOrb2 = iAtomStart(iAtom2f+1) - jj

        if (blockList(iAtom2f,1) == 0) then
          blockList(iAtom2f,1) = iNext
          iNext = iNext + nOrb2
        end if

        if ( .not. isBlockInLocal(jj,jj+nOrb2-1, ii, ii+nOrb1-1, colStartLocal, colEndLocal)) then
          cycle
        end if

        iVec = iCellVec(iAtom2)
        phase = exp(cmplx(0, 1, dp) * dot_product(kPoint2p, cellVec(:, iVec)))

        call addBlock2Elsi(phase*reshape(orig(iOrig:iOrig+nOrb1*nOrb2-1), [ nOrb2, nOrb1 ]),&
            & colStartLocal, colEndLocal, jj, ii, colPtrLocal, blockList(iAtom2f,1)-1, nzValLocal,&
            & rowIndLocal)

        if (ii == jj) then
          cycle
        end if

        ! other triangle, so Hermitian symmetry
        phase = conjg(phase)

        ! Because of folding of periodic images, it can happen that the transposed block has already
        ! been registered.
        if (.not. tRowTrans(iAtom2f)) then
          blockList(iAtom2f,2) = blockList(iAtom2f,2) + nOrb1
          tRowTrans(iAtom2f) = .true.
        end if

        call addBlock2Elsi(&
            & phase * transpose(reshape(orig(iOrig:iOrig+nOrb1*nOrb2-1), [ nOrb2, nOrb1 ])),&
            & colStartLocal, colEndLocal, ii, jj, colPtrLocal, blockList(iAtom2f,2)-nOrb1-1,&
            & nzValLocal, rowIndLocal)
      end do
    end do

  end subroutine pack2elsi_cmplx


  !> Checks if atom block is part of local matrix (Elsi)
  pure logical function isBlockInLocal(rowStartBlock, rowEndBlock, colStartBlock, colEndBlock,&
      & colStartLocal, colEndLocal)

    !> Row of global matrix where block starts
    integer, intent(in) :: rowStartBlock

    !> Row of global matrix where block ends
    integer, intent(in) :: rowEndBlock

    !> Column of global matrix where block starts
    integer, intent(in) :: colStartBlock

    !> Column of global matrix where block ends
    integer, intent(in) :: colEndBlock

    !> Column of global matrix where local part starts
    integer, intent(in) :: colStartLocal

    !> Column of global matrix where local part ends
    integer, intent(in) :: colEndLocal

    isBlockInLocal = .false.
    if ( (colStartLocal <= colStartBlock) .and. (colStartBlock <= colEndLocal) ) then
      isBlockInLocal = .true.
    else if ( (colStartLocal <= colEndBlock) .and. (colEndBlock <= colEndLocal) ) then
      isBlockInLocal = .true.
    else if ( (colStartBlock <= colStartLocal) .and. (colStartLocal <= colEndBlock)) then
      isBlockInLocal = .true.
    else if ( (colStartBlock <= colEndLocal) .and. (colEndLocal <= colEndBlock)) then
      isBlockInLocal = .true.
    else if ( (colStartLocal <= rowStartBlock) .and. (rowStartBlock <= colEndLocal) ) then
      isBlockInLocal = .true.
    else if ( (colStartLocal <= rowEndBlock) .and. (rowEndBlock <= colEndLocal) ) then
      isBlockInLocal = .true.
    else if ( (rowStartBlock <= colStartLocal) .and. (colStartLocal <= rowEndBlock)) then
      isBlockInLocal = .true.
    else if ( (rowStartBlock <= colEndLocal) .and. (colEndLocal <= rowEndBlock)) then
      isBlockInLocal = .true.
    end if

  end function isBlockInLocal


  !> Checks whether global column is in local matrix
  pure logical function isColumnInLocal(col, colStartLocal, colEndLocal)

    !> global column
    integer, intent(in) :: col

    !> Column of global matrix where local matrix starts
    integer, intent(in) :: colStartLocal

    !> Column of global matrix where local matrix ends
    integer, intent(in) :: colEndLocal

    if ( (colStartLocal <= col) .and. (col <= colEndLocal) ) then
      isColumnInLocal = .true.
    else
      isColumnInLocal = .false.
    end if

  end function isColumnInLocal

#:for SUFFIX in FLAVOURS

  !> Add the content of a local matrix block to ELSI CSC format
  subroutine addBlock2Elsi${SUFFIX}$(loc, colStart, colEnd, ii, jj, colptr, rowOffset, val,&
      & rowIndLocal)

    !> Local matrix.
    ${SUFFIX}$(dp), intent(in) :: loc(:,:)

    !> Column of global matrix where local matrix starts
    integer, intent(in) :: colStart

    !> Column of global matrix where local matrix ends
    integer, intent(in) :: colEnd

    !> Starting row in the global matrix.
    integer, intent(in) :: ii

    !> Starting column in the global matrix
    integer, intent(in) :: jj

    !> column pointer
    integer, intent(in) :: colptr(:)

    !> index of the next element per column
    integer, intent(in) :: rowOffset

    !> values in CSC format
    ${SUFFIX}$(dp), intent(inout) :: val(:)

    !> row index pointer
    integer, intent(inout), optional :: rowIndLocal(:)

    integer :: j2, iloc, jloc
    integer :: jStart, jEnd

    jStart = max(jj, colStart) - colStart + 1
    jEnd = min(jj+size(loc,dim=2)-1, colEnd) - colStart + 1
    jloc = max(jj, colStart) - jj + 1

    if (present(rowIndLocal)) then
      do j2 = jStart, jEnd
        do iloc = 1, size(loc,1)
          val(colptr(j2)+rowOffset+iloc-1) = val(colptr(j2) + rowOffset+iloc-1) + loc(iloc,jloc)
          rowIndLocal(colptr(j2) + rowOffset + iloc-1) = iloc + ii - 1
        end do
        jloc = jloc + 1
      end do
    else
      do j2 = jStart, jEnd
        do iloc = 1, size(loc,1)
          val(colptr(j2)+rowOffset+iloc-1) = val(colptr(j2) + rowOffset+iloc-1) + loc(iloc,jloc)
        end do
        jloc = jloc + 1
      end do
    end if

  end subroutine addBlock2Elsi${SUFFIX}$
#:endfor

  !> Convert CSC matrix format into DFTB+ sparse format
  !>
  !> Note: primitive will not be set to zero on startup, and values are added to enable addition of
  !> spin components. Make sure, you set it to zero before invoking this routine the first time.
  subroutine elsi2pack_real(self, colStart, colEnd, iNeighbour, nNeighbourSK, mOrb,&
      & iAtomStart, iSparseStart, img2CentCell, colptr, nzval, blockRow, primitive)

    !> Sparse conversion instance
    type(TSparse2Sparse), intent(in) :: self

    !> Column of global matrix where local matrix starts
    integer, intent(in) :: colStart

    !> Column of global matrix where local matrix ends
    integer, intent(in) :: colEnd

    !> Neighbour list for each atom (first index from 0!).
    integer, intent(in) :: iNeighbour(0:,:)

    !> Nr. of neighbours for each atom (incl. itself).
    integer, intent(in) :: nNeighbourSK(:)

    !> Maximal number of orbitals on an atom.
    integer, intent(in) :: mOrb

    !> Atom offset for the squared Hamiltonian
    integer, intent(in) :: iAtomStart(:)

    !> Indexing array for the sparse Hamiltonian
    integer, intent(in) :: iSparseStart(0:,:)

    !> Atomic mapping indexes.
    integer, intent(in) :: img2CentCell(:)

    !> Local column pointer
    integer, intent(in) :: colptr(:)

    !> Local non-zero elements
    real(dp), intent(in) :: nzval(:)

    !> Saves starting row of blocks in nzValLocal
    integer, intent(in) :: blockRow(0:,:)

    !> Sparse Hamiltonian
    real(dp), intent(inout) :: primitive(:)

    integer :: nAtom
    integer :: iOrig, ii, jj, kk
    integer :: iNeigh
    integer :: iAt, iAtom1, iAtom2, iAtom2f
    integer :: nOrb1, nOrb2
    real(dp) :: tmpSqr(mOrb,mOrb)

    nAtom = size(iNeighbour, dim=2)

    tmpSqr(:,:) = 0.0_dp

    ! loop over relevant atoms to pack back
    do iAt = 1, self%nAtomsInColumns
      iAtom1 = self%atomsInColumns(iAt)
      ii = iAtomStart(iAtom1)
      nOrb1 = iAtomStart(iAtom1+1) - ii
      do iNeigh = 0, nNeighbourSK(iAtom1)
        iOrig = iSparseStart(iNeigh,iAtom1) + 1
        iAtom2 = iNeighbour(iNeigh, iAtom1)
        iAtom2f = img2CentCell(iAtom2)
        jj = iAtomStart(iAtom2f)
        nOrb2 = iAtomStart(iAtom2f+1) - jj

        if ( .not. isBlockInLocal(jj,jj+nOrb2-1, ii, ii+nOrb1-1, colStart, colEnd)) then
          cycle
        end if

        call cpElsi2Block(colStart, colEnd, colptr, nzval, ii, blockRow(iNeigh, iAtom1),&
            & tmpSqr(1:nOrb2,1:nOrb1))

        ! Symmetrize the on-site block before packing, just in case
        if (iAtom1 == iAtom2f) then
          do kk = 1, nOrb2
            tmpSqr(kk, kk+1:nOrb1) = tmpSqr(kk+1:nOrb1, kk)
          end do
        end if

        primitive(iOrig : iOrig + nOrb1 * nOrb2 - 1) = primitive(iOrig : iOrig + nOrb1 * nOrb2 - 1)&
            & + reshape(tmpSqr(1:nOrb2,1:nOrb1), [ nOrb1*nOrb2 ])

      end do
    end do

  end subroutine elsi2pack_real


  !> Convert CSC matrix format into DFTB+ sparse format
  !>
  !> Note: primitive will not be set to zero on startup, and values are added to enable addition of
  !> spin components. Make sure, you set it to zero before invoking this routine the first time.
  subroutine elsi2pack_cmplx(self, colStart, colEnd, iNeighbour, nNeighbourSK, mOrb,&
      & iAtomStart, iSparseStart, img2CentCell, kPoint, kWeight, iCellVec, cellVec, colptr, nzval,&
      & blockRow, primitive)

    !> Sparse conversion instance
    type(TSparse2Sparse), intent(in) :: self

    !> Column of global matrix where local matrix starts
    integer, intent(in) :: colStart

    !> Column of global matrix where local matrix ends
    integer, intent(in) :: colEnd

    !> Neighbour list for each atom (first index from 0!).
    integer, intent(in) :: iNeighbour(0:,:)

    !> Nr. of neighbours for each atom (incl. itself).
    integer, intent(in) :: nNeighbourSK(:)

    !> Maximal number of orbitals on an atom.
    integer, intent(in) :: mOrb

    !> Atom offset for the squared Hamiltonian
    integer, intent(in) :: iAtomStart(:)

    !> Indexing array for the sparse Hamiltonian
    integer, intent(in) :: iSparseStart(0:,:)

    !> Atomic mapping indexes.
    integer, intent(in) :: img2CentCell(:)

    !> Current k-point
    real(dp), intent(in) :: kPoint(:)

    !> Weight for current k-points
    real(dp), intent(in) :: kWeight

    !> Index for which unit cell atoms are associated with
    integer, intent(in) :: iCellVec(:)

    !> Vectors (in units of the lattice constants) to cells of the lattice
    real(dp), intent(in) :: cellVec(:,:)

    !> Local column pointer
    integer, intent(in) :: colptr(:)

    !> Local non-zero elements
    complex(dp), intent(in) :: nzval(:)

    !> Saves starting row of blocks in nzValLocal
    integer, intent(in) :: blockRow(0:,:)

    !> Sparse Hamiltonian
    real(dp), intent(inout) :: primitive(:)

    integer :: nAtom, iVec
    integer :: iOrig, ii, jj, kk
    integer :: iNeigh
    integer :: iAt, iAtom1, iAtom2, iAtom2f
    integer :: nOrb1, nOrb2
    complex(dp) :: tmpSqr(mOrb,mOrb)
    real(dp) :: kPoint2p(3)
    complex(dp) :: phase

    kPoint2p(:) = 2.0_dp * pi * kPoint(:)
    nAtom = size(iNeighbour, dim=2)

    tmpSqr(:,:) = cmplx(0,0,dp)

    ! loop over relevant atoms to pack back
    do iAt = 1, self%nAtomsInColumns
      iAtom1 = self%atomsInColumns(iAt)
      ii = iAtomStart(iAtom1)
      nOrb1 = iAtomStart(iAtom1+1) - ii
      do iNeigh = 0, nNeighbourSK(iAtom1)
        iOrig = iSparseStart(iNeigh,iAtom1) + 1
        iAtom2 = iNeighbour(iNeigh, iAtom1)
        iAtom2f = img2CentCell(iAtom2)
        jj = iAtomStart(iAtom2f)
        nOrb2 = iAtomStart(iAtom2f+1) - jj

        if ( .not. isBlockInLocal(jj,jj+nOrb2-1, ii, ii+nOrb1-1, colStart, colEnd)) then
          cycle
        end if

        call cpElsi2Block(colStart, colEnd, colptr, nzval, ii, blockRow(iNeigh, iAtom1),&
            & tmpSqr(1:nOrb2,1:nOrb1))

        iVec = iCellVec(iAtom2)
        phase = exp(cmplx(0, -1, dp) * dot_product(kPoint2p, cellVec(:, iVec)))

        tmpSqr(1:nOrb2,1:nOrb1) = phase * tmpSqr(1:nOrb2,1:nOrb1)

        ! Hermitian conjugate the on-site block before packing, just in case
        if (iAtom1 == iAtom2f) then
          do kk = 1, nOrb2
            tmpSqr(kk, kk+1:nOrb1) = conjg(tmpSqr(kk+1:nOrb1, kk))
          end do
        end if

        primitive(iOrig : iOrig + nOrb1 * nOrb2 - 1) = primitive(iOrig : iOrig + nOrb1 * nOrb2 - 1)&
            & + kWeight * real(reshape(tmpSqr(1:nOrb2,1:nOrb1), [ nOrb1*nOrb2 ]))

      end do
    end do

  end subroutine elsi2pack_cmplx


#:for SUFFIX in FLAVOURS

  !> Copies the content from the ELSI structure to block
  subroutine cpElsi2Block${SUFFIX}$(colStart, colEnd, colptr, nzval, jj, blockRow, loc)

    !> Column of global matrix where local matrix starts
    integer, intent(in) :: colStart

    !> Column of global matrix where local matrix ends
    integer, intent(in) :: colEnd

    !> column pointer
    integer, intent(in) :: colptr(:)

    !> non-zero values
    ${SUFFIX}$(dp), intent(in) :: nzval(:)

    !> Starting column in the global matrix
    integer, intent(in) :: jj

    !> Saves starting row of blocks in nzValLocal
    integer, intent(in) :: blockRow

    !> Local block of matrix.
    ${SUFFIX}$(dp), intent(out) :: loc(:,:)

    integer :: j2, iloc

    loc(:,:) = 0.0_dp

    do j2 = 1, size(loc, dim=2)
      if ( isColumnInLocal(jj + j2 - 1, colStart, colEnd )) then
        iloc = blockRow + colptr(jj-colStart+j2) - 1
        loc(:, j2) = nzval(iloc:iloc+size(loc,dim=1)-1)
      end if
    end do

  end subroutine cpElsi2Block${SUFFIX}$

#:endfor

end module elsisolver
