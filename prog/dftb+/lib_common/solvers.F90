!--------------------------------------------------------------------------------------------------!
!  DFTB+: general package for performing fast atomistic simulations                                !
!  Copyright (C) 2018  DFTB+ developers group                                                      !
!                                                                                                  !
!  See the LICENSE file for terms of usage and distribution.                                       !
!--------------------------------------------------------------------------------------------------!

#:include 'common.fypp'

!> Contains computer environment settings
module solvers
  use accuracy
  use message
  use environment
  use elsiiface
  implicit none

  private
  public :: TElectronicSolverInp
  public :: TElectronicSolver, TElectronicSolver_init
  public :: electronicSolverTypes


  type :: TElsiSolverInp

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


  !> Input for electronic/eigen solver block
  type :: TElectronicSolverInp

    !> Solver type
    integer :: iSolver

    type(TElsiSolverInp), allocatable :: elsi

  end type TElectronicSolverInp


  !> Contains settings for the solvers of the ELSI library
  type :: TElsiSolver

    !> Handle for the ELSI library
    type(elsi_handle), public :: handle

    !> Use sparse CSR format
    logical, public :: CSR = .false.

    !> solver number
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

  end type TElsiSolver


  !> Eigensolver state and settings
  type :: TElectronicSolver

    !> Electronic solver number
    integer, public :: iSolver

    !> Whether it is an ELSI solver
    logical, public :: isElsiSolver

    !> Whether the solver provides eigenvalues
    logical, public :: providesEigenvals

    !> Are Choleskii factors already available for the overlap matrix
    logical, public, allocatable :: tCholeskiiDecomposed(:)

    type(TElsiSolver), allocatable :: elsi

  contains

    procedure :: initElsi => TElectronicSolver_initElsi
    procedure :: resetElsi => TElectronicSolver_resetElsi
    procedure :: finalElsi => TElectronicSolver_finalElsi

    procedure :: getSolverName => TElectronicSolver_getSolverName
    procedure, private :: ensurePauliCompatibility => TElectronicSolver_ensurePauliCompatibility

  end type TElectronicSolver


  !> Namespace for possible solver methods
  type :: TElectronicSolverTypesEnum
    integer :: qr = 1
    integer :: divideandconquer = 2
    integer :: relativelyrobust = 3
    ! elsi provided solvers
    integer :: elpa = 4
    integer :: omm = 5
    integer :: pexsi = 6
    integer :: dummy1 = 7
    integer :: dummy2 = 8
    integer :: ntpoly = 9
    ! transport related
    integer :: gf = 10
    integer :: onlyTransport = 11
  end type TElectronicSolverTypesEnum

  !> Actual values for electronicSolverTypes.
  type(TElectronicSolverTypesEnum), parameter :: electronicSolverTypes =&
      & TElectronicSolverTypesEnum()

contains

  subroutine TElectronicSolver_init(this, iSolver)
    type(TElectronicSolver), intent(out) :: this
    integer, intent(in) :: iSolver

    this%iSolver = iSolver
    this%isElsiSolver = any(this%iSolver ==&
        & [electronicSolverTypes%elpa, electronicSolverTypes%omm, electronicSolverTypes%pexsi,&
        & electronicSolverTypes%ntpoly])
    this%providesEigenvals = any(this%iSolver ==&
        & [electronicSolverTypes%qr, electronicSolverTypes%divideandconquer,&
        & electronicSolverTypes%relativelyrobust, electronicSolverTypes%elpa])

  end subroutine TElectronicSolver_init


  !> Initialise extra settings relevant to ELSI in the solver data structure
  subroutine TElectronicSolver_initElsi(this, inp, env, nBasisFn, nEl, iDistribFn,&
      & nSpin, nKPoint, tWriteHS)

    !> control structure for solvers, including ELSI data
    class(TElectronicSolver), intent(inout) :: this

    !> input structure for ELSI
    type(TElectronicSolverInp), intent(in) :: inp

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

    allocate(this%elsi)

    ! DFTB+ to ELSI solver index
    select case(this%iSolver)
    case (electronicSolverTypes%elpa)
      this%elsi%SOLVER = 1
      ! ELPA is asked for all states
      this%elsi%n_state = nBasisFn
    case (electronicSolverTypes%omm)
      this%elsi%SOLVER = 2
      ! OMM solves only over occupied space
      ! spin degeneracies for closed shell
      if (nSpin == 1) then
        this%elsi%n_state = nint(sum(nEl)*0.5_dp)
      else
        this%elsi%n_state = nint(sum(nEl))
      end if
    case (electronicSolverTypes%pexsi)
      this%elsi%SOLVER = 3
      ! ignored by PEXSI:
      this%elsi%n_state = nBasisFn
    case (electronicSolverTypes%ntpoly)
      this%elsi%SOLVER = 6
      ! ignored by NTPoly:
      this%elsi%n_state = nBasisFn
    end select

    ! parallelism with multiple processes
    this%elsi%parallel = 1

    ! number of basis functions
    this%elsi%n_basis = nBasisFn

    ! number of electrons in the problem
    this%elsi%n_electron = sum(nEl)

    ! should the code write the hamiltonian and overlap and stop
    this%elsi%tWriteHS = tWriteHS

    if (nSpin == 2) then
      this%elsi%spin_degeneracy = nEl(1) - nEl(2)
    else
      this%elsi%spin_degeneracy = 0.0_dp
    end if

    this%elsi%n_spin = nSpin
    this%elsi%n_kpoint = nKPoint

    this%elsi%MPI_COMM_WORLD = env%mpi%globalComm%id
    this%elsi%my_COMM_WORLD = env%mpi%groupComm%id
    this%elsi%my_BLACS_Ctxt = env%blacs%orbitalGrid%ctxt

    ! assumes row and column sizes the same
    this%elsi%BLACS_blockSize = env%blacs%rowBlockSize

    this%elsi%CSR_blockSize = this%elsi%n_basis / env%mpi%groupComm%size

    if (this%elsi%CSR_blockSize * env%mpi%groupComm%size < this%elsi%n_basis) then
      this%elsi%CSR_blockSize = this%elsi%CSR_blockSize + 1
    end if

    this%elsi%mu_broaden_scheme = min(iDistribFn,2)
    if (iDistribFn > 1) then
      ! set Meth-Pax order
      this%elsi%mu_mp_order = iDistribFn - 2
    else
      this%elsi%mu_mp_order = 0
    end if

    ! ELPA settings
    this%elsi%ELPA_SOLVER_Option = inp%elsi%ELPA_Solver

    ! OMM settings
    this%elsi%OMM_iter = inp%elsi%OMM_IterationsELPA
    this%elsi%OMM_Tolerance = inp%elsi%OMM_Tolerance
    this%elsi%OMM_Choleskii = inp%elsi%OMM_Choleskii

    ! PEXSI settings
    this%elsi%PEXSI_n_pole = inp%elsi%PEXSI_n_pole
    this%elsi%PEXSI_np_per_pole = inp%elsi%PEXSI_np_per_pole
    this%elsi%PEXSI_n_mu = inp%elsi%PEXSI_n_mu
    this%elsi%PEXSI_np_symbo = inp%elsi%PEXSI_np_symbo
    this%elsi%PEXSI_delta_e = inp%elsi%PEXSI_delta_e

    ! NTPoly settings
    this%elsi%NTPoly_method = inp%elsi%NTPoly_method
    this%elsi%NTPoly_truncation = inp%elsi%NTPoly_truncation
    this%elsi%NTPoly_tolerance = inp%elsi%NTPoly_tolerance

    this%elsi%CSR = inp%elsi%ELSI_CSR

    ! data as dense BLACS blocks
    if (this%elsi%CSR) then
      ! CSR format
      this%elsi%BLACS_DENSE = 2
    else
      this%elsi%BLACS_DENSE = 0
    end if

    ! customize output level, note there are levels 0..3 not DFTB+ 0..2
    this%elsi%OutputLevel = 0
  #:call DEBUG_CODE
    this%elsi%OutputLevel = 3
  #:endcall DEBUG_CODE

    if (nSpin == 4) then
      call this%ensurePauliCompatibility()
    end if

  end subroutine TElectronicSolver_initElsi


  !> reset the ELSI solver - safer to do this on geometry change, due to the lack of a Choleskii
  !> refactorization option
  subroutine TElectronicSolver_resetElsi(this, tempElec, iSpin, iKPoint, kWeight)

    !> Instance
    class(TElectronicSolver), intent(inout) :: this

    !> electron temperature
    real(dp), intent(in) :: tempElec

    !> current spin value, set to be 1 if non-collinear
    integer, intent(in) :: iSpin

    !> current k-point value
    integer, intent(in) :: iKPoint

    !> weight for current k-point
    real(dp), intent(in) :: kWeight

    if (this%elsi%nResets > 0) then
      ! destroy previous instance of solver if called before
      call elsi_finalize(this%elsi%handle)
    end if
    this%elsi%nResets = this%elsi%nResets + 1

    call elsi_init(this%elsi%handle, this%elsi%SOLVER, this%elsi%parallel, this%elsi%BLACS_DENSE,&
        & this%elsi%n_basis, this%elsi%n_electron, this%elsi%n_state)

    call elsi_set_mpi_global(this%elsi%handle, this%elsi%MPI_COMM_WORLD)
    call elsi_set_sing_check(this%elsi%handle, 0) ! disable singularity check
    call elsi_set_mpi(this%elsi%handle, this%elsi%my_COMM_WORLD)

    if (this%elsi%CSR) then
      call elsi_set_csc_blk(this%elsi%handle, this%elsi%CSR_blockSize)
    end if
    call elsi_set_blacs(this%elsi%handle, this%elsi%my_BLACS_Ctxt, this%elsi%BLACS_blockSize)

    if (this%elsi%tWriteHS) then
      ! setup to write a matrix
      call elsi_init_rw(this%elsi%rwHandle, 1, this%elsi%parallel, this%elsi%n_basis,&
          & this%elsi%n_electron)
      ! MPI comm
      call elsi_set_rw_mpi(this%elsi%rwHandle, this%elsi%MPI_COMM_WORLD)
      if (.not.this%elsi%CSR) then
        ! dense matrices
        call elsi_set_rw_blacs(this%elsi%rwHandle, this%elsi%my_BLACS_Ctxt,&
            & this%elsi%BLACS_blockSize)
      end if
    end if


    call elsi_set_mu_broaden_scheme(this%elsi%handle, this%elsi%mu_broaden_scheme)
    if (this%elsi%mu_broaden_scheme == 2) then
      call elsi_set_mu_mp_order(this%elsi%handle, this%elsi%mu_mp_order)
    end if

    ! set filling temperature/width
    call elsi_set_mu_broaden_width(this%elsi%handle, tempElec)

    select case(this%iSolver)
    case(electronicSolverTypes%elpa)

      select case(this%elsi%ELPA_SOLVER_Option)
      case(1)
        call elsi_set_elpa_solver(this%elsi%handle, 1)
      case(2)
        call elsi_set_elpa_solver(this%elsi%handle, 2)
      case default
        call error("Unknown ELPA solver modes")
      end select

    case(electronicSolverTypes%omm)
      ! libOMM
      if (this%elsi%OMM_Choleskii) then
        call elsi_set_omm_flavor(this%elsi%handle, 2)
      else
        call elsi_set_omm_flavor(this%elsi%handle, 0)
      end if
      call elsi_set_omm_n_elpa(this%elsi%handle, this%elsi%OMM_iter)
      call elsi_set_omm_tol(this%elsi%handle, this%elsi%OMM_Tolerance)

    case(electronicSolverTypes%pexsi)

      this%elsi%PEXSI_mu_min = -10.0_dp
      this%elsi%PEXSI_mu_max = 10.0_dp
      this%elsi%PEXSI_DeltaVmin = 0.0_dp
      this%elsi%PEXSI_DeltaVmax = 0.0_dp

      ! processors per pole to invert for
      call elsi_set_pexsi_np_per_pole(this%elsi%handle, this%elsi%PEXSI_np_per_pole)

      call elsi_set_pexsi_mu_min(this%elsi%handle, this%elsi%PEXSI_mu_min +this%elsi%PEXSI_DeltaVmin)
      call elsi_set_pexsi_mu_max(this%elsi%handle, this%elsi%PEXSI_mu_max +this%elsi%PEXSI_DeltaVmax)

      ! set filling temperature/width
      call elsi_set_pexsi_temp(this%elsi%handle, tempElec)

      ! number of poles for the expansion
      call elsi_set_pexsi_n_pole(this%elsi%handle, this%elsi%PEXSI_n_pole)

      ! number of interpolation points for mu
      call elsi_set_pexsi_n_mu(this%elsi%handle, this%elsi%PEXSI_n_mu)

      ! number of processors for symbolic factorisation task
      call elsi_set_pexsi_np_symbo(this%elsi%handle, this%elsi%PEXSI_np_symbo)

      ! spectral radius (range of eigenspectrum, if known, otherwise defaul usually fine)
      call elsi_set_pexsi_delta_e(this%elsi%handle, this%elsi%PEXSI_delta_e)

    case(electronicSolverTypes%ntpoly)

      ! NTPoly
      ! set purification method
      call elsi_set_ntpoly_method(this%elsi%handle, this%elsi%NTPoly_method)

      ! set truncation tolerance for sparse matrix multiplications
      call elsi_set_ntpoly_filter(this%elsi%handle, this%elsi%NTPoly_truncation)

      ! set purification convergence threshold
      call elsi_set_ntpoly_tol(this%elsi%handle, this%elsi%NTPoly_tolerance)

    end select

    if (any(this%iSolver == [electronicSolverTypes%omm,&
        & electronicSolverTypes%pexsi, electronicSolverTypes%ntpoly])) then
      ! density matrix build needs to know the number of spin channels to normalize against
      select case(this%elsi%n_spin)
      case(1)
        call elsi_set_spin(this%elsi%handle, 1, iSpin)
      case(2)
        call elsi_set_spin(this%elsi%handle, 2, iSpin)
      case(4)
        call elsi_set_spin(this%elsi%handle, 2, iSpin)
      end select
      if (this%elsi%n_kpoint > 1) then
        call elsi_set_kpoint(this%elsi%handle, this%elsi%n_kpoint, iKPoint, kWeight)
      end if
    end if

    call elsi_set_output(this%elsi%handle, this%elsi%OutputLevel)
    if (this%elsi%OutputLevel == 3) then
      call elsi_set_output_log(this%elsi%handle, 1)
    end if

  end subroutine TElectronicSolver_resetElsi


  subroutine TElectronicSolver_finalElsi(this)
    class(TElectronicSolver), intent(inout) :: this

    call elsi_finalize(this%elsi%handle)

  end subroutine TElectronicSolver_finalElsi


  function TElectronicSolver_getSolverName(this) result(solverName)
    class(TElectronicSolver), intent(in) :: this
    character(:), allocatable :: solverName

    character(lc) :: buffer

    select case (this%iSolver)

    case(electronicSolverTypes%qr)
      write(buffer, "(A)") "Standard"

    case(electronicSolverTypes%divideandconquer)
      write(buffer, "(A)") "Divide and Conquer"

    case(electronicSolverTypes%relativelyrobust)
      write(buffer, "(A)") "Relatively robust"

    case(electronicSolverTypes%elpa)
      if (this%elsi%ELPA_SOLVER_Option == 1) then
        write(buffer, "(A)") "ELSI interface to the 1 stage ELPA solver"
      else
        write(buffer, "(A)") "ELSI interface to the 2 stage ELPA solver"
      end if

    case(electronicSolverTypes%omm)
      write(buffer, "(A,I0,A,E8.2)") "ELSI solver libOMM with ",&
          & this%elsi%OMM_iter, " ELPA iterations",this%elsi%OMM_Tolerance
      if (this%elsi%CSR) then
        write(buffer, "(A)") "ELSI solver libOMM Sparse"
      else
        write(buffer, "(A)") "ELSI solver libOMM Dense"
      end if

    case(electronicSolverTypes%pexsi)
      if (this%elsi%CSR) then
        write(buffer, "(A)") "ELSI solver PEXSI Sparse"
      else
        write(buffer, "(A)") "ELSI solver PEXSI Dense"
      end if

    case(electronicSolverTypes%ntpoly)
      if (this%elsi%CSR) then
        write(buffer, "(A)") "ELSI solver NTPoly Sparse"
      else
        write(buffer, "(A)") "ELSI solver NTPoly Dense"
      end if

    case(electronicSolverTypes%gf)
      write(buffer, "(A)") "Green's functions"

    case(electronicSolverTypes%onlyTransport)
      write(buffer, "(A)") "Transport Only (no energies)"

    case default
      write(buffer, "(A,I0,A)") "Invalid electronic solver! (iSolver = ", this%iSolver, ")"

    end select
    solverName = trim(buffer)

  end function TElectronicSolver_getSolverName


  subroutine TElectronicSolver_ensurePauliCompatibility(this)
    class(TElectronicSolver), intent(in) :: this

    logical :: tPauliIncompat

    tPauliIncompat = this%elsi%CSR&
        & .and. any(this%iSolver == [electronicSolverTypes%omm, electronicSolverTypes%pexsi,&
        & electronicSolverTypes%ntpoly])
    if (tPauliIncompat) then
      call error("Current solver configuration not avaible for two component complex&
          & hamiltonians")
    end if

  end subroutine TElectronicSolver_ensurePauliCompatibility


end module solvers
