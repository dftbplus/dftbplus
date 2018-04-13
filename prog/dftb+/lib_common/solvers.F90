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
#:if WITH_ELSI
  use ELSI
#:endif
  implicit none

  !> Eigensolver state and settings
  type :: TElectronicSolver

    !> Is the ELSI solver being used
    logical :: tUsingELSI = .false.

    !> Are Chileskii factors already available for the overlap matrix
    logical, public, allocatable :: tCholeskiiDecomposed(:)

  #:if WITH_ELSI
    !> Handle for the ELSI library
    type(elsi_handle), public :: elsiHandle
    integer, public :: ELSI_SOLVER
    integer, public :: ELSI_SOLVER_Option
    integer, public :: ELSI_OutputLevel
    integer, public :: ELSI_parallel
    integer, public :: ELSI_BLACS_DENSE
    integer, public :: ELSI_n_basis
    real(dp), public :: ELSI_n_electron
    integer, public :: ELSI_n_state
    integer, public :: ELSI_MPI_COMM_WORLD
    integer, public :: ELSI_my_COMM_WORLD
    integer, public :: ELSI_blockSize
  #:endif

  end type TElectronicSolver

end module solvers
