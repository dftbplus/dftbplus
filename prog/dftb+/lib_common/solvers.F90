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
    integer, public :: ELSI_blockSize

    integer, public :: ELSI_ELPA_SOLVER_Option

    integer, public :: ELSI_OMM_iter
    real(dp), public :: ELSI_OMM_Tolerance
    logical, public :: ELSI_OMM_Choleskii

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

  subroutine resetELSI(this)
    class(TElectronicSolver), intent(inout) :: this

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
    call elsi_set_blacs(this%elsiHandle,this%ELSI_MPI_COMM_WORLD, this%ELSI_blockSize)
    select case(this%ELSI_SOLVER)
    case(1)
      select case(this%ELSI_ELPA_SOLVER_Option)
      case(1)
        call elsi_set_elpa_solver(this%elsiHandle, 1)
      case(2)
        call elsi_set_elpa_solver(this%elsiHandle, 2)
      case default
        call error("Unknown ELPA solver modes")
      end select
    case(2)
      if (this%ELSI_OMM_Choleskii) then
        call elsi_set_omm_flavor(this%elsiHandle, 2)
      else
        call elsi_set_omm_flavor(this%elsiHandle, 0)
      end if
      call elsi_set_omm_n_elpa(this%elsiHandle, this%ELSI_OMM_iter)
      call elsi_set_omm_tol(this%elsiHandle, this%ELSI_OMM_Tolerance)
    end select
    call elsi_set_output(this%elsiHandle, this%ELSI_OutputLevel)

  end subroutine resetELSI

#:endif

end module solvers
