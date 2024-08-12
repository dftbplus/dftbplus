!--------------------------------------------------------------------------------------------------!
!  DFTB+: general package for performing fast atomistic simulations                                !
!  Copyright (C) 2006 - 2023  DFTB+ developers group                                               !
!                                                                                                  !
!  See the LICENSE file for terms of usage and distribution.                                       !
!--------------------------------------------------------------------------------------------------!

#:include 'common.fypp'

!> Contains data type representing the input data for DFTB.
module modes_inputdata
  use dftbp_common_accuracy, only : dp
  use dftbp_type_typegeometry, only : TGeometry
  implicit none

  private
  public :: solverTypes
  public :: TInputData, TInputData_init
  public :: TControl, TParallelOpts, TBlacsOpts


  !> Contains Blacs specific options.
  type :: TBlacsOpts

    !> Block size for matrix rows and columns.
    integer :: blockSize

  end type TBlacsOpts


  !> Contains the parallel options.
  type :: TParallelOpts

    !> Number of processor groups
    integer :: nGroup

    !> Blacs options
    type(TBlacsOpts) :: blacsOpts

    !> Whether hybrid parallelisation is enable
    logical :: tOmpThreads

  end type TParallelOpts


  !> Main control data for program as extracted by the parser.
  type TControl

    !> Geometry
    type(TGeometry) :: geo

    !> Atomic masses to build dynamical matrix
    real(dp), allocatable :: atomicMasses(:)

    !> File name of Hessian (direct read)
    character(len=:), allocatable :: hessianFile

    !> Hessian matrix
    real(dp), allocatable :: hessian(:,:)

    !> Born charges matrix
    real(dp), allocatable :: bornMatrix(:)

    !> Derivatives of Born charges matrix with respect to electric field, i.e. polarizability
    !! derivatives with respect to atom locations
    real(dp), allocatable :: bornDerivsMatrix(:)

    !> Produce plots of modes
    logical :: tPlotModes

    !> Modes to produce xyz file for
    integer, allocatable :: modesToPlot(:)

    !> If animating, number of cycles to show in an animation
    integer :: nCycles = 3

    !> Steps in an animation cycle
    integer :: nSteps = 10

    !> List of atoms in dynamical matrix
    integer, allocatable :: iMovedAtoms(:)

    !> Animate mode  or as vectors
    logical :: tAnimateModes

    !> Remove translation modes
    logical :: tRemoveTranslate

    !> Remove rotation modes
    logical :: tRemoveRotate

    !> File access type to use when opening binary files for reading and writing
    character(20) :: binaryAccessTypes(2)

    !> Parallelization related data
    type(TParallelOpts), allocatable :: parallelOpts

    !> Eigensolver choice
    integer :: iSolver

  end type TControl


  !> Container for input data constituents.
  type TInputData

    !> Initialized?
    logical :: tInitialized = .false.

    !> Main control data
    type(TControl) :: ctrl

    !> Geometry data
    type(TGeometry) :: geo

  end type TInputData


  !> Namespace for possible eigensolver methods
  type :: TSolverTypesEnum

    !> Lapack/Scalapack solvers
    integer :: qr = 1
    integer :: divideAndConquer = 2
    integer :: relativelyRobust = 3

    !> GPU accelerated solver using MAGMA
    integer :: magmaEvd = 4

  end type TSolverTypesEnum


  !> Actual values for solverTypes.
  type(TSolverTypesEnum), parameter :: solverTypes = TSolverTypesEnum()


contains

  !> Mark data structure as initialised.
  subroutine TInputData_init(this)

    !> Instance
    type(TInputData), intent(out) :: this

    this%tInitialized = .true.

  end subroutine TInputData_init

end module modes_inputdata
