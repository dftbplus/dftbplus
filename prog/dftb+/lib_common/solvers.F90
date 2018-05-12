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
#:if WITH_ELSI
  use ELSI
#:endif
  implicit none

  private
  public :: TElectronicSolverInp, TElectronicSolver

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

    integer, public :: ELSI_SOLVER
    integer, public :: ELSI_OutputLevel
    integer, public :: ELSI_parallel
    integer, public :: ELSI_BLACS_DENSE
    integer, public :: ELSI_n_basis
    real(dp), public :: ELSI_n_electron
    integer, public :: ELSI_n_state
    integer, public :: ELSI_MPI_COMM_WORLD
    integer, public :: ELSI_my_COMM_WORLD
    integer, public :: ELSI_my_BLACS_Ctxt
    integer, public :: ELSI_blockSize

    integer, public :: ELSI_mu_broaden_scheme
    integer, public :: ELSI_mu_mp_order

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
    real(dp), public, allocatable :: ELSI_PEXSI_VOld(:)
    integer, public :: ELSI_PEXSI_n_pole
    integer, public :: ELSI_PEXSI_n_mu
    integer, public :: ELSI_PEXSI_np_symbo
    real(dp), public :: ELSI_PEXSI_delta_e

    !> count of the number of times ELSI has been reset (usually every geometry step)
    integer :: nELSI_reset = 0

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

  subroutine resetELSI(this, tempElec)

    !> Instance
    class(TElectronicSolver), intent(inout) :: this

    !> electron temperature
    real(dp), intent(in) :: tempElec

    if (this%nELSI_reset > 0) then
      ! destroy previous instance of solver
      call elsi_finalize(this%elsiHandle)
    end if
    this%nELSI_reset = this%nELSI_reset + 1

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

      call elsi_set_pexsi_np_per_pole(this%elsiHandle, 1)

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
    call elsi_set_output(this%elsiHandle, this%ELSI_OutputLevel)
    if (this%ELSI_OutputLevel == 3) then
      call elsi_set_output_log(this%elsiHandle, 1)
    end if

  end subroutine resetELSI

#:endif

end module solvers
