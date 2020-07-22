!--------------------------------------------------------------------------------------------------!
!  DFTB+: general package for performing fast atomistic simulations                                !
!  Copyright (C) 2006 - 2020  DFTB+ developers group                                               !
!                                                                                                  !
!  See the LICENSE file for terms of usage and distribution.                                       !
!--------------------------------------------------------------------------------------------------!

#:include 'common.fypp'

!> Contains computer environment settings
module dftbp_elecsolvers
  use dftbp_accuracy, only : dp, lc
  use dftbp_assert
  use dftbp_elecsolvertypes, only : electronicSolverTypes
  use dftbp_elsisolver
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

    !> Whether the solver provides a band energy (sum of the eigenvalues)
    logical, public :: providesBandEnergy

    !> Whether the solver provides the TS term for electrons
    logical, public :: providesElectronEntropy

    !> Whether the solver provides electronic free energy (or is consistent with its evaluation)
    logical, public :: providesFreeEnergy

    !> Whether the solver provides the electron chemical potential (or is consistent with its
    !> evaluation or it being suplied externally)
    logical, public :: elecChemPotAvailable

    !> Data for ELSI solvers
    type(TElsiSolver), public, allocatable :: elsi

    !> Are Choleskii factors already available for the overlap matrix
    logical, public, allocatable :: hasCholesky(:)

    !> Buffer for storing real overlap matrices between calls
    real(dp), allocatable :: choleskyBufferReal(:,:,:)

    !> Buffer for storing complex overlap matrices between calls
    complex(dp), allocatable :: choleskyBufferCmplx(:,:,:)

    !> Number of buffered overlap matrices
    integer :: nCholesky

  contains
    procedure :: getSolverName => TElectronicSolver_getSolverName
    procedure :: reset => TElectronicSolver_reset
    procedure, private :: getCholeskyReal => TElectronicSolver_getCholeskyReal
    procedure, private :: getCholeskyCmplx => TElectronicSolver_getCholeskyCmplx
    procedure, private :: storeCholeskyReal => TElectronicSolver_storeCholeskyReal
    procedure, private :: storeCholeskyCmplx => TElectronicSolver_storeCholeskyCmplx
    generic :: storeCholesky => storeCholeskyReal, storeCholeskyCmplx
    generic :: getCholesky => getCholeskyReal, getCholeskyCmplx
    procedure :: updateElectronicTemp => TElectronicSolver_updateElectronicTemp

  end type TElectronicSolver


contains

  !> Initializes an electronic solver
  subroutine TElectronicSolver_init(this, iSolver, nCholesky)

    !> Instance.
    type(TElectronicSolver), intent(out) :: this

    !> Solver type to be used.
    integer, intent(in) :: iSolver

    !> Number of Cholesky-decompositions which will be buffered.
    integer, intent(in) :: nCholesky

    this%iSolver = iSolver

    this%isElsiSolver = any(this%iSolver ==&
        & [electronicSolverTypes%elpa, electronicSolverTypes%omm, electronicSolverTypes%pexsi,&
        & electronicSolverTypes%ntpoly, electronicSolverTypes%elpadm])

    !> Eigenvalues for hamiltonian available
    this%providesEigenvals = any(this%iSolver ==&
        & [electronicSolverTypes%qr, electronicSolverTypes%divideandconquer,&
        & electronicSolverTypes%relativelyrobust, electronicSolverTypes%elpa,&
        & electronicSolverTypes%magma_gvd])

    !> Band energy for electrons available
    this%providesBandEnergy = any(this%iSolver ==&
        & [electronicSolverTypes%qr, electronicSolverTypes%divideandconquer,&
        & electronicSolverTypes%relativelyrobust, electronicSolverTypes%elpa,&
        & electronicSolverTypes%elpadm, electronicSolverTypes%ntpoly,&
        & electronicSolverTypes%magma_gvd, electronicSolverTypes%pexsi])

    !> TS term for electrons is available
    this%providesElectronEntropy = any(this%iSolver ==&
        & [electronicSolverTypes%qr, electronicSolverTypes%divideandconquer,&
        & electronicSolverTypes%relativelyrobust, electronicSolverTypes%elpa,&
        & electronicSolverTypes%elpadm, electronicSolverTypes%magma_gvd,&
        & electronicSolverTypes%pexsi])

    !> Electron chemical potential is either available or provided externally. Note this can get
    !> over-riden in initprogram (e.g. if the boundary conditions are contacts with different
    !> chemical potentials which prevents evaluation of mu * N_elec terms)
    this%elecChemPotAvailable = any(this%iSolver ==&
        & [electronicSolverTypes%qr, electronicSolverTypes%divideandconquer,&
        & electronicSolverTypes%relativelyrobust, electronicSolverTypes%elpa,&
        & electronicSolverTypes%elpadm, electronicSolverTypes%ntpoly,&
        & electronicSolverTypes%magma_gvd, electronicSolverTypes%pexsi, electronicSolverTypes%gf])

    !> The electronic Helmholtz free energy of the system is available (U - TS + mu N_elec). Note
    !> that chemical potential logical can be re-defined elsewhere.
    this%providesFreeEnergy =  this%providesElectronEntropy .and. this%elecChemPotAvailable

    this%nCholesky = nCholesky
    allocate(this%hasCholesky(this%nCholesky))
    this%hasCholesky(:) = .false.

    if (this%isElsiSolver) then
      allocate(this%elsi)
    end if

  end subroutine TElectronicSolver_init


  !> Resets the electronic solver for the next geometry step.
  subroutine TElectronicSolver_reset(this)

    !> Instance.
    class(TElectronicSolver), intent(inout) :: this

    this%hasCholesky(:) = .false.
    if (this%isElsiSolver) then
      call this%elsi%reset()
    end if

    ! The electronic Helmholtz free energy of the system is available (U - TS + mu N_elec). Note
    ! that chemical potential logical can be re-defined.
    this%providesFreeEnergy =  this%providesElectronEntropy .and. this%elecChemPotAvailable

  end subroutine TElectronicSolver_reset


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

    case(electronicSolverTypes%magma_gvd)
      write(buffer, "(A)") "Divide and Conquer (MAGMA GPU version)"

    case default
      write(buffer, "(A,I0,A)") "Invalid electronic solver! (iSolver = ", this%iSolver, ")"

    end select

    solverName = trim(buffer)

  end function TElectronicSolver_getSolverName


#:for SUFFIX, TYPE in [('Real', 'real'), ('Cmplx', 'complex')]

  !> Stores a Cholesky-decomposition.
  subroutine TElectronicSolver_storeCholesky${SUFFIX}$(this, iCholesky, cholesky)

    !> Instance.
    class(TElectronicSolver), intent(inout) :: this

    !> Serial number of the cholesky decomposition to store.
    integer, intent(in) :: iCholesky

    !> Cholesky decomposition of the overlap
    ${TYPE}$(dp), intent(in) :: cholesky(:,:)

    @:ASSERT(iCholesky >= 1 .and. iCholesky <= this%nCholesky)

    if (.not. allocated(this%choleskyBuffer${SUFFIX}$)) then
      allocate(this%choleskyBuffer${SUFFIX}$&
          &(size(cholesky, dim=1), size(cholesky, dim=2), this%nCholesky))
    end if
    this%choleskyBuffer${SUFFIX}$(:,:,iCholesky) = cholesky
    this%hasCholesky(iCholesky) = .true.

  end subroutine TElectronicSolver_storeCholesky${SUFFIX}$


  !> Returns a stored Cholesky decomposition
  subroutine TElectronicSolver_getCholesky${SUFFIX}$(this, iCholesky, cholesky)

    !> Instance.
    class(TElectronicSolver), intent(inout) :: this

    !> Serial number of the cholesky decomposition to store.
    integer, intent(in) :: iCholesky

    !> Cholesky decomposition of the overlap
    ${TYPE}$(dp), intent(out) :: cholesky(:,:)

    @:ASSERT(iCholesky >= 1 .and. iCholesky <= this%nCholesky)
    @:ASSERT(allocated(this%choleskyBuffer${SUFFIX}$))

    cholesky(:,:) = this%choleskyBuffer${SUFFIX}$(:,:,iCholesky)

  end subroutine TElectronicSolver_getCholesky${SUFFIX}$

#:endfor

  !> Updates the electronic temperatures for the solvers
  subroutine TElectronicSolver_updateElectronicTemp(this, tempElec)

    !> Instance.
    class(TElectronicSolver), intent(inout) :: this

    !> Electronic temperature
    real(dp), intent(in) :: tempElec

    if (this%isElsiSolver) then
      call this%elsi%updateElectronicTemp(tempElec)
    end if

  end subroutine TElectronicSolver_updateElectronicTemp


end module dftbp_elecsolvers
