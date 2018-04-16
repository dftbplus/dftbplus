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

  !> Eigensolver state and settings
  type :: TElectronicSolver

    !> Is the ELSI solver being used
    logical :: tUsingELSI = .false.

    !> Are Choleskii factors already available for the overlap matrix
    logical, public, allocatable :: tCholeskiiDecomposed(:)

    integer, public :: ELSI_SOLVER_Option

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

    integer :: nELSI_reset = 0

  contains

    procedure :: resetELSI

  #:endif

  end type TElectronicSolver

#:if WITH_ELSI

contains

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
    if (this%ELSI_SOLVER == 1) then
      select case(this%ELSI_SOLVER_option)
      case(1)
        call elsi_set_elpa_solver(this%elsiHandle, 1)
      case(2)
        call elsi_set_elpa_solver(this%elsiHandle, 2)
      case default
        call error("Unknown ELPA solver modes")
      end select
    end if
    call elsi_set_output(this%elsiHandle, this%ELSI_OutputLevel)

  end subroutine resetELSI

#:endif

end module solvers
