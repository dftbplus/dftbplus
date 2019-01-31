!--------------------------------------------------------------------------------------------------!
!  DFTB+: general package for performing fast atomistic simulations                                !
!  Copyright (C) 2018  DFTB+ developers group                                                      !
!                                                                                                  !
!  See the LICENSE file for terms of usage and distribution.                                       !
!--------------------------------------------------------------------------------------------------!

#:include 'common.fypp'

!> Contains computer environment settings
module elecsolvers
  use accuracy, only : lc
  use elecsolvertypes, only : electronicSolverTypes
  use elsisolver
  implicit none

  private
  public :: TElectronicSolverInp
  public :: TElectronicSolver, TElectronicSolver_init
  public :: TElsiSolverInp
  public :: electronicSolverTypes


  !> Input for electronic/eigen solver block
  type :: TElectronicSolverInp

    !> Solver type
    integer :: iSolver

    !> Input for the ELSI solver
    type(TElsiSolverInp), allocatable :: elsi

  end type TElectronicSolverInp


  !> Eigensolver state and settings
  type :: TElectronicSolver
    private

    !> Electronic solver number
    integer, public :: iSolver

    !> Whether it is an ELSI solver
    logical, public :: isElsiSolver

    !> Whether the solver provides eigenvalues
    logical, public :: providesEigenvals

    !> Are Choleskii factors already available for the overlap matrix
    logical, public, allocatable :: tCholeskiiDecomposed(:)

    type(TElsiSolver), public, allocatable :: elsi

  contains
    procedure :: getSolverName => TElectronicSolver_getSolverName

  end type TElectronicSolver


contains

  !> Initializes an electronic solver
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

    if (this%isElsiSolver) then
      allocate(this%elsi)
    end if

  end subroutine TElectronicSolver_init


  !> Returns the name of the solver used.
  function TElectronicSolver_getSolverName(this) result(solverName)

    !> Instance.
    class(TElectronicSolver), intent(in) :: this

    !> Name of the solver.
    character(:), allocatable :: solverName

    character(lc) :: buffer

    if (this%isElsiSolver) then
      solverName = this%elsi%getSolverName()
      return
    end if

    select case (this%iSolver)

    case(electronicSolverTypes%qr)
      write(buffer, "(A)") "Standard"

    case(electronicSolverTypes%divideandconquer)
      write(buffer, "(A)") "Divide and Conquer"

    case(electronicSolverTypes%relativelyrobust)
      write(buffer, "(A)") "Relatively robust"

    case(electronicSolverTypes%gf)
      write(buffer, "(A)") "Green's functions"

    case(electronicSolverTypes%onlyTransport)
      write(buffer, "(A)") "Transport Only (no energies)"

    case default
      write(buffer, "(A,I0,A)") "Invalid electronic solver! (iSolver = ", this%iSolver, ")"

    end select

    solverName = trim(buffer)

  end function TElectronicSolver_getSolverName


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!! TElsiSolver
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


end module elecsolvers
