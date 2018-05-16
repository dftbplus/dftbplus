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
#:if WITH_ELSI
  use ELSI
#:endif
  implicit none

  private
  public :: TElectronicSolverInp, TElectronicSolver
#:if WITH_ELSI
  public :: init
#:endif

#:if WITH_ELSI
  interface init
    module procedure init_ELSI
  end interface init
#:endif

  !> Input for electronic/eigen solver block
  type :: TElectronicSolverInp

    !> Solver type
    integer :: iSolver

  #:if WITH_ELSI

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

  #:endif

  end type TElectronicSolverInp

  !> Eigensolver state and settings
  type :: TElectronicSolver

    !> Is the ELSI solver being used
    logical :: tUsingELSI = .false.

    !> Are Choleskii factors already available for the overlap matrix
    logical, public, allocatable :: tCholeskiiDecomposed(:)

    !> Electronic solver number
    integer :: iSolver

  #:if WITH_ELSI

    !> Handle for the ELSI library
    type(elsi_handle), public :: elsiHandle

    !> solver number
    integer :: ELSI_SOLVER

    !> level of output from solver
    integer :: ELSI_OutputLevel

    !> parallelisation strategy
    integer :: ELSI_parallel

    !> Dense BLACS
    integer :: ELSI_BLACS_DENSE

    integer :: ELSI_n_basis
    real(dp) :: ELSI_n_electron
    real(dp) :: ELSI_spin_degeneracy
    integer :: ELSI_n_state
    integer :: ELSI_MPI_COMM_WORLD
    integer :: ELSI_my_COMM_WORLD
    integer :: ELSI_my_BLACS_Ctxt
    integer :: ELSI_blockSize

    integer :: ELSI_n_spin
    integer :: ELSI_n_kpoint

    integer :: ELSI_mu_broaden_scheme
    integer :: ELSI_mu_mp_order

    ! ELPA settings
    integer, public :: ELSI_ELPA_SOLVER_Option

    ! OMM settings
    integer, public :: ELSI_OMM_iter
    real(dp), public :: ELSI_OMM_Tolerance
    logical, public :: ELSI_OMM_Choleskii

    ! PEXSI settings
    real(dp), public :: ELSI_PEXSI_mu_min
    real(dp), public :: ELSI_PEXSI_mu_max
    real(dp), public :: ELSI_PEXSI_DeltaVmin
    real(dp), public :: ELSI_PEXSI_DeltaVmax
    real(dp), allocatable, public :: ELSI_PEXSI_VOld(:)
    integer, public :: ELSI_PEXSI_n_pole
    integer, public :: ELSI_PEXSI_np_per_pole
    integer, public :: ELSI_PEXSI_n_mu
    integer, public :: ELSI_PEXSI_np_symbo
    real(dp) :: ELSI_PEXSI_delta_e

    !> count of the number of times ELSI has been reset (usually every geometry step)
    integer :: nELSI_resets = 0

  contains

    procedure :: resetELSI

  #:endif

  end type TElectronicSolver

  !> Namespace for possible solver methods
  type :: electronicSolverTypesEnum
    integer :: qr
    integer :: divideandconquer
    integer :: relativelyrobust
    integer :: elpa
    integer :: omm
    integer :: pexsi
  end type electronicSolverTypesEnum

  !> Actual values for electronicSolverTypes.
  type(electronicSolverTypesEnum), parameter :: electronicSolverTypes =&
      & electronicSolverTypesEnum(1, 2, 3, 4, 5, 6)

contains

#: if WITH_ELSI

  !> Initialise extra settings relevant to ELSI in the solver data structure
  subroutine init_ELSI(inp, this, env, nBasisFn, nEl, iDistribFn, tWriteDetailedOutBands,&
      & nIndepHam, nSpin, nKPoint)

    !> input structure for ELSI
    type(TElectronicSolverInp), intent(in) :: inp

    !> control structure for solvers, including ELSI data
    type(TElectronicSolver), intent(inout) :: this

    !> Environment settings
    type(TEnvironment), intent(in) :: env

    !> number of orbitals in the system
    integer, intent(in) :: nBasisFn

    !> number of electrons
    real(dp), intent(in) :: nEl(:)

    !> filling function
    integer, intent(in) :: iDistribFn

    !> Should bands be produced?
    logical, intent(inout) :: tWriteDetailedOutBands

    !> total number of independent spin channels. In the case of non-collinear, should pass 1
    integer, intent(in) :: nIndepHam

    !> total number of spin channels.
    integer, intent(in) :: nSpin

    !> total number of k-points
    integer, intent(in) :: nKPoint

    ! use of ELSI to solve electronic states
    this%tUsingELSI = .true.

    ! DFTB+ to ELSI solver index
    this%ELSI_SOLVER = this%iSolver -3

    select case(this%ELSI_SOLVER)
    case (1)
      ! ELPA is asked for all states
      this%ELSI_n_state = nBasisFn
    case (2)
      ! OMM solves only over occupied space
      ! spin degeneracies for closed shell
      if (nSpin == 1) then
        this%ELSI_n_state = nint(sum(nEl)*0.5_dp)
      else
        this%ELSI_n_state = nint(sum(nEl))
      end if
    case (3)
      ! PEXSI ignores this
      this%ELSI_n_state = nBasisFn
    end select

    ! bands only available for ELPA
    if (this%ELSI_SOLVER > 1) then
      tWriteDetailedOutBands = .false.
    end if

    ! data as dense BLACS blocks
    this%ELSI_BLACS_DENSE = 0

    ! parallelism with multiple processes
    this%ELSI_parallel = 1

    ! number of basis functions
    this%ELSI_n_basis = nBasisFn

    ! number of electrons in the problem
    this%ELSI_n_electron = sum(nEl)

    if (nSpin == 2) then
      this%ELSI_spin_degeneracy = nEl(1) - nEl(2)
    else
      this%ELSI_spin_degeneracy = 0.0_dp
    end if

    this%ELSI_n_spin = nIndepHam
    this%ELSI_n_kpoint = nKPoint

    this%ELSI_MPI_COMM_WORLD = env%mpi%globalComm%id
    this%ELSI_my_COMM_WORLD = env%mpi%groupComm%id
    this%ELSI_my_BLACS_Ctxt = env%blacs%orbitalGrid%ctxt

    ! assumes row and column sizes the same
    this%ELSI_blockSize = env%blacs%rowBlockSize

    this%ELSI_mu_broaden_scheme = min(iDistribFn,2)
    if (iDistribFn > 1) then
      ! set Meth-Pax order
      this%ELSI_mu_mp_order = iDistribFn - 2
    else
      this%ELSI_mu_mp_order = 0
    end if

    ! ELPA settings
    this%ELSI_ELPA_SOLVER_Option = inp%ELPA_Solver

    ! OMM settings
    this%ELSI_OMM_iter = inp%OMM_IterationsELPA
    this%ELSI_OMM_Tolerance = inp%OMM_Tolerance
    this%ELSI_OMM_Choleskii = inp%OMM_Choleskii

    ! PEXSI settings
    this%ELSI_PEXSI_n_pole = inp%PEXSI_n_pole
    this%ELSI_PEXSI_np_per_pole = inp%PEXSI_np_per_pole
    this%ELSI_PEXSI_n_mu = inp%PEXSI_n_mu
    this%ELSI_PEXSI_np_symbo = inp%PEXSI_np_symbo
    this%ELSI_PEXSI_delta_e = inp%PEXSI_delta_e

    ! customize output level, note there are levels 0..3 not DFTB+ 0..2
    this%ELSI_OutputLevel = 0
  #:call DEBUG_CODE
    this%ELSI_OutputLevel = 3
  #:endcall DEBUG_CODE

  end subroutine init_ELSI

  !> reset the ELSI solver - safer to do this on geometry change, due to the lack of a Choleskii
  !> refactorization option
  subroutine resetELSI(this, tempElec, iSpin, iKPoint, kWeight)

    !> Instance
    class(TElectronicSolver), intent(inout) :: this

    !> electron temperature
    real(dp), intent(in) :: tempElec

    !> current spin value, se to be 1 if non-collinear
    integer, intent(in) :: iSpin

    !> current k-point value
    integer, intent(in) :: iKPoint

    !> weight for current k-point
    real(dp), intent(in) :: kWeight

    if (this%nELSI_resets > 0) then
      ! destroy previous instance of solver if called before
      call elsi_finalize(this%elsiHandle)
    end if
    this%nELSI_resets = this%nELSI_resets + 1

    call elsi_init(this%elsiHandle, this%ELSI_SOLVER, this%ELSI_parallel, this%ELSI_BLACS_DENSE,&
        & this%ELSI_n_basis, this%ELSI_n_electron, this%ELSI_n_state)
    call elsi_set_mpi_global(this%elsiHandle, this%ELSI_MPI_COMM_WORLD)
    call elsi_set_sing_check(this%elsiHandle, 0) ! disable singularity check
    call elsi_set_mpi(this%elsiHandle, this%ELSI_my_COMM_WORLD)
    call elsi_set_blacs(this%elsiHandle, this%ELSI_my_BLACS_Ctxt, this%ELSI_blockSize)

    call elsi_set_mu_broaden_scheme(this%elsiHandle, this%ELSI_mu_broaden_scheme)
    if (this%ELSI_mu_broaden_scheme == 2) then
      call elsi_set_mu_mp_order(this%elsiHandle, this%ELSI_mu_mp_order)
    end if

    ! set filling temperature/width
    call elsi_set_mu_broaden_width(this%elsiHandle, tempElec)

    select case(this%ELSI_SOLVER)
    case(1)
      ! ELPA
      select case(this%ELSI_ELPA_SOLVER_Option)
      case(1)
        call elsi_set_elpa_solver(this%elsiHandle, 1)
      case(2)
        call elsi_set_elpa_solver(this%elsiHandle, 2)
      case default
        call error("Unknown ELPA solver modes")
      end select
    case(2)
      ! libOMM
      if (this%ELSI_OMM_Choleskii) then
        call elsi_set_omm_flavor(this%elsiHandle, 2)
      else
        call elsi_set_omm_flavor(this%elsiHandle, 0)
      end if
      call elsi_set_omm_n_elpa(this%elsiHandle, this%ELSI_OMM_iter)
      call elsi_set_omm_tol(this%elsiHandle, this%ELSI_OMM_Tolerance)

    case(3)
      ! PEXSI
      this%ELSI_PEXSI_mu_min = -10.0_dp
      this%ELSI_PEXSI_mu_max = 10.0_dp
      this%ELSI_PEXSI_DeltaVmin = 0.0_dp
      this%ELSI_PEXSI_DeltaVmax = 0.0_dp

      ! processors per pole to invert for
      call elsi_set_pexsi_np_per_pole(this%elsiHandle, this%ELSI_PEXSI_np_per_pole)

      call elsi_set_pexsi_mu_min(this%elsiHandle, this%ELSI_PEXSI_mu_min +this%ELSI_PEXSI_DeltaVmin)
      call elsi_set_pexsi_mu_max(this%elsiHandle, this%ELSI_PEXSI_mu_max +this%ELSI_PEXSI_DeltaVmax)

      ! set filling temperature/width
      call elsi_set_pexsi_temp(this%elsiHandle, tempElec)

      ! number of poles for the expansion
      call elsi_set_pexsi_n_pole(this%elsiHandle, this%ELSI_PEXSI_n_pole)

      ! number of interpolation points for mu
      call elsi_set_pexsi_n_mu(this%elsiHandle, this%ELSI_PEXSI_n_mu)

      ! number of processors for symbolic factorisation task
      call elsi_set_pexsi_np_symbo(this%elsiHandle, this%ELSI_PEXSI_np_symbo)

      ! spectral radius (range of eigenspectrum, if known, otherwise defaul usually fine)
      call elsi_set_pexsi_delta_e(this%elsiHandle, this%ELSI_PEXSI_delta_e)

    end select

    if (this%ELSI_SOLVER > 1) then
      if (this%ELSI_n_spin > 1) then
        call elsi_set_spin(this%elsiHandle, this%ELSI_n_spin, iSpin)
        call elsi_set_mu_spin_degen(this%elsiHandle, this%ELSI_spin_degeneracy)
      end if
      if (this%ELSI_n_kpoint > 1) then
        call elsi_set_kpoint(this%elsiHandle, this%ELSI_n_kpoint, iKPoint, kWeight)
      end if
    end if

    call elsi_set_output(this%elsiHandle, this%ELSI_OutputLevel)
    if (this%ELSI_OutputLevel == 3) then
      call elsi_set_output_log(this%elsiHandle, 1)
    end if

  end subroutine resetELSI

#:endif

end module solvers
