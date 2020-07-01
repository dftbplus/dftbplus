!--------------------------------------------------------------------------------------------------!
!  DFTB+: general package for performing fast atomistic simulations                                !
!  Copyright (C) 2018  DFTB+ developers group                                                      !
!                                                                                                  !
!  See the LICENSE file for terms of usage and distribution.                                       !
!--------------------------------------------------------------------------------------------------!

#:include 'common.fypp'

!> Contains eigensolver related data types.
module dftbp_solvers
#:if WITH_PROGRESS
  use sp2progress
#:endif
  implicit none
  private

  public :: solverTypes
  public :: TSolver

  !> Helper types listing possible solver types as fields with appropriate names
  type :: TSolverTypesHelper
    integer :: lapackQr
    integer :: lapackDivAndConq
    integer :: lapackRelRobust
    integer :: progressSp2
  end type TSolverTypesHelper

  !> Contains possible solver types
  type(TSolverTypesHelper), parameter :: solverTypes = TSolverTypesHelper(1, 2, 3, 4)


  !> Contains an actual solver
  type :: TSolver

    !> Type number of the solver (one of the fields of solverTypes)
    integer :: solverType
    
  #:if WITH_PROGRESS
    !> An SP2 solver instance
    type(TSp2Solver), allocatable :: sp2Solver
  #:endif

  end type TSolver


end module dftbp_solvers
