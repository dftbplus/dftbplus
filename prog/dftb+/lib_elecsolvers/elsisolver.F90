!--------------------------------------------------------------------------------------------------!
!  DFTB+: general package for performing fast atomistic simulations                                !
!  Copyright (C) 2018  DFTB+ developers group                                                      !
!                                                                                                  !
!  See the LICENSE file for terms of usage and distribution.                                       !
!--------------------------------------------------------------------------------------------------!

#:include 'common.fypp'

!> Contains the interface to the ELSI solvers
module elsisolver
  use accuracy, only : dp, lc
  use environment, only : TEnvironment
  use elecsolvertypes, only : electronicSolverTypes
  use elsiiface
  use message, only : error
  implicit none
  private

  public :: TElsiSolverInp
  public :: TElsiSolver
  public :: TElsiSolver_init, TElsiSolver_final


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

  contains
    procedure :: reset => TElsiSolver_reset
    procedure :: getSolverName => TElsiSolver_getSolverName
    procedure, private :: ensurePauliCompatibility => TElsiSolver_ensurePauliCompatibility
  end type TElsiSolver


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

end module elsisolver
