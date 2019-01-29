!--------------------------------------------------------------------------------------------------!
!  DFTB+: general package for performing fast atomistic simulations                                !
!  Copyright (C) 2018  DFTB+ developers group                                                      !
!                                                                                                  !
!  See the LICENSE file for terms of usage and distribution.                                       !
!--------------------------------------------------------------------------------------------------!

#:include 'common.fypp'

!> Contains computer environment settings
module elecsolvers
  use accuracy, only : dp, lc
  use environment, only : TEnvironment
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

    procedure :: initElsi => TElectronicSolver_initElsi
    procedure :: resetElsi => TElectronicSolver_resetElsi
    procedure :: finalElsi => TElectronicSolver_finalElsi

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
    call TElsiSolver_init(this%elsi, inp%elsi, env, nBasisFn, nEl, iDistribFn, nSpin, nKPoint,&
        & tWriteHS)

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

    call this%elsi%reset(tempElec, iSpin, iKPoint, kWeight)

  end subroutine TElectronicSolver_resetElsi


  !> Finalizes the ELSI module.
  subroutine TElectronicSolver_finalElsi(this)
    class(TElectronicSolver), intent(inout) :: this

    call TElsiSolver_final(this%elsi)

  end subroutine TElectronicSolver_finalElsi


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
