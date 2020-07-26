!--------------------------------------------------------------------------------------------------!
!  DFTB+: general package for performing fast atomistic simulations                                !
!  Copyright (C) 2006 - 2020  DFTB+ developers group                                               !
!                                                                                                  !
!  See the LICENSE file for terms of usage and distribution.                                       !
!--------------------------------------------------------------------------------------------------!

#:include 'common.fypp'

!> Contains the interface to the ELSI solvers
module dftbp_elsisolver
#:if WITH_MPI
  use dftbp_mpifx
#:endif
  use dftbp_accuracy, only : dp, lc
  use dftbp_environment, only : TEnvironment, globalTimers
  use dftbp_globalenv, only : stdOut
  use dftbp_elecsolvertypes, only : electronicSolverTypes
  use dftbp_elsiiface
  use dftbp_elsicsc
  use dftbp_densedescr
  use dftbp_periodic
  use dftbp_orbitals
  use dftbp_message, only : error, warning, cleanshutdown
  use dftbp_commontypes, only : TParallelKS, TOrbitals
  use dftbp_energytypes, only : TEnergies
  use dftbp_etemp, only : fillingTypes
  use dftbp_sparse2dense
  use dftbp_assert
  use dftbp_spin, only : ud2qm
  use dftbp_angmomentum, only : getLOnsite
  use dftbp_spinorbit, only : addOnsiteSpinOrbitHam, getOnsiteSpinOrbitEnergy
  use dftbp_potentials, only : TPotentials
  implicit none
  private

  public :: TElsiSolverInp
  public :: TElsiSolver, TElsiSolver_init, TElsiSolver_final


  !> Input data for the ELSI solvers
  type :: TElsiSolverInp

    !> Choice of the solver
    integer :: iSolver

    !> Choice of ELPA solver
    integer :: elpaSolver = 2

    !> Iterations of ELPA solver before OMM minimization
    integer :: ommIterationsElpa = 5

    !> Halting tolerance for OMM iterations
    real(dp) :: ommTolerance = 1.0E-10_dp

    !> Should the overlap be factorized before minimization
    logical :: ommCholesky = .true.

    !> PEXSI pole expansion method
    integer :: pexsiMethod = 3

    !> number of poles for PEXSI expansion
    integer :: pexsiNPole = 30

    !> number of processors per pole for PEXSI
    integer :: pexsiNpPerPole = 1

    !> number of interpolation points for mu (Fermi) search
    integer :: pexsiNMu = 2

    !> number of processors for symbolic factorisation
    integer :: pexsiNpSymbo = 1

    !> spectral radius (range of eigenvalues) if available
    real(dp) :: pexsiDeltaE = 10.0_dp

    !> density matrix purification algorithm
    integer :: ntpolyMethod = 2

    !> truncation threshold for sparse matrix multiplication
    real(dp) :: ntpolyTruncation = 1.0E-10_dp

    !> convergence tolerance for density matrix purification
    real(dp) :: ntpolyTolerance = 1.0E-5_dp

    !> Use sparse CSR format
    logical :: elsiCsr = .false.

    !> Tolerance for converting from dense matrices to internal sparse storage for libOMM, PEXSI and
    !> NTPoly.
    real(dp) :: elsi_zero_def

  end type TElsiSolverInp


  !> Contains settings for the solvers of the ELSI library. See ELSI manual for detailed meaning of
  !> variables
  type :: TElsiSolver

    private

    !> should the code write matrices and stop
    logical, public :: tWriteHS

    !> handle for matrix IO
    type(elsi_rw_handle), public :: rwHandle

    !> Solver id (using the central DFTB+ solver type)
    integer :: iSolver

    !> Handle for the ELSI library
    type(elsi_handle), public :: handle

    !> Use sparse CSR format
    logical, public :: isSparse = .false.

    !> Tolerance for converting from dense matrices to internal sparse storage for libOMM, PEXSI and
    !> NTPoly.
    real(dp) :: elsi_zero_def

    !> ELSI Solver choice
    integer :: solver

    !> Level of solver information output
    integer :: outputLevel

    !> See ELSI manual - parallelisation strategy
    integer :: parallel

    !> See ELSI manual - type of parallel data decomposition
    integer :: denseBlacs

    !> Number of basis functions
    integer :: nBasis

    !> Number of electrons in the system
    real(dp) :: nElectron

    !> Maximum spin occupation for levels
    real(dp) :: spinDegeneracy

    !> Number of states to solve when relevant
    integer :: nState

    !> Global comm
    integer :: mpiCommWorld

    !> Group comm
    integer :: myCommWorld

    !> BLACS matrix context in use
    integer :: myBlacsCtx

    !> BLACS block sizes
    integer :: BlacsBlockSize

    !> CSC blocksize
    integer :: csrBlockSize

    !> Number of independent spins
    integer :: nSpin

    !> Number of independent k-points
    integer :: nKPoint

    !> Index for current spin processed here
    integer :: iSpin

    !> Index for current k-point processed here
    integer :: iKPoint

    !> Weighting for current k-point
    real(dp) :: kWeight

    !> Choice of broadening function
    integer :: muBroadenScheme

    !> If Meth-Paxton, order of scheme
    integer :: muMpOrder

    !> Whether solver had been already initialised
    logical :: tSolverInitialised = .false.

    !> Has overlap been factorized
    logical :: tCholeskyDecomposed

    !> Is this the first calls for this geometry
    logical :: tFirstCalc = .true.

    !> Sparse format data structure
    type(TElsiCsc), allocatable :: elsiCsc


    !> Version number for ELSI
    integer, public :: major, minor, patch

    !! ELPA settings

    !> ELPA solver choice
    integer :: elpaSolverOption

    !! OMM settings

    !> Starting iterations with ELPA
    integer :: ommIter

    !> Tolerance for minimisation
    real(dp) :: ommTolerance

    !> Whether to Cholesky factorize and transform or work with general
    logical :: ommCholesky

    !! PEXSI settings

    !> Minimum contour range
    real(dp) :: pexsiMuMin

    !> Maximum contour range
    real(dp) :: pexsiMuMax

    !> Most negative potential
    real(dp) :: pexsiDeltaVMin

    !> Most positive potential
    real(dp) :: pexsiDeltaVMax

    !> Previous potentials
    real(dp), allocatable :: pexsiVOld(:)

    !> PEXSI pole expansion method
    integer :: pexsiMethod

    !> Number of poles for expansion
    integer :: pexsiNPole

    !> Processors per pole
    integer :: pexsiNpPerPole

    !> Processes for chemical potential search
    integer :: pexsiNMu

    !> Processors used for symbolic factorization stage
    integer :: pexsiNpSymbo

    !> Spectral radius
    real(dp) :: pexsiDeltaE

    !! NTPoly settings

    !> Choice of minimisation method
    integer :: ntpolyMethod

    !> Truncation threshold for elements
    real(dp) :: ntpolyTruncation

    !> Convergence tolerance
    real(dp) :: ntpolyTolerance

  contains

    procedure :: reset => TElsiSolver_reset
    procedure :: updateGeometry => TElsiSolver_updateGeometry
    procedure :: updateElectronicTemp => TElsiSolver_updateElectronicTemp
    procedure :: getDensity => TElsiSolver_getDensity
    procedure :: getEDensity => TElsiSolver_getEDensity
    procedure :: initPexsiDeltaVRanges => TElsiSolver_initPexsiDeltaVRanges
    procedure :: updatePexsiDeltaVRanges => TElsiSolver_updatePexsiDeltaVRanges
    procedure :: getSolverName => TElsiSolver_getSolverName

  end type TElsiSolver


contains


  !> Initialise ELSI solver
  subroutine TElsiSolver_init(this, inp, env, nBasisFn, nEl, iDistribFn, nSpin, iSpin, nKPoint,&
      & iKPoint, kWeight, tWriteHS, providesElectronEntropy)

    !> control structure for solvers, including ELSI data
    class(TElsiSolver), intent(out) :: this

    !> input structure for ELSI
    type(TElsiSolverInp), intent(in) :: inp

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

    !> spin channel processed by current process
    integer, intent(in) :: iSpin

    !> total number of k-points
    integer, intent(in) :: nKPoint

    !> K-point processed by current process.
    integer, intent(in) :: iKPoint

    !> Weight of current k-point
    real(dp), intent(in) :: kWeight

    !> Should the matrices be written out
    logical, intent(in) :: tWriteHS

    !> Whether the solver provides the TS term for electrons
    logical, intent(inout) :: providesElectronEntropy

    integer :: dateStamp

  #:if WITH_ELSI

    character(lc) :: buffer

    call elsi_get_version(this%major, this%minor, this%patch)
    call elsi_get_datestamp(dateStamp)

    call supportedVersionNumber(this, dateStamp)

    this%iSolver = inp%iSolver

    select case(this%iSolver)

    case (electronicSolverTypes%elpa)
      this%solver = 1
      ! ELPA is asked for all states
      this%nState = nBasisFn

    case (electronicSolverTypes%omm)
      this%solver = 2
      ! OMM solves only over occupied space
      ! spin degeneracies for closed shell
      if (nSpin == 1) then
        this%nState = nint(sum(nEl)*0.5_dp)
      else
        this%nState = nint(sum(nEl))
      end if

    case (electronicSolverTypes%pexsi)
      this%solver = 3
      ! ignored by PEXSI:
      this%nState = nBasisFn

    case (electronicSolverTypes%ntpoly)
      this%solver = 6
      ! ignored by NTPoly, but set anyway:
      this%nState = nBasisFn

    case (electronicSolverTypes%elpadm)
      this%solver = 1
      ! ignored by density matrix from ELPA, but set anyway:
      this%nState = nBasisFn

    end select

    ! parallelism with multiple processes
    this%parallel = 1

    ! number of basis functions
    this%nBasis = nBasisFn

    ! number of electrons in the problem
    this%nElectron = sum(nEl)

    ! should the code write the hamiltonian and overlap then stop
    this%tWriteHS = tWriteHS

    if (nSpin == 2) then
      this%spinDegeneracy = nEl(1) - nEl(2)
    else
      this%spinDegeneracy = 0.0_dp
    end if

    ! Number of spin channels passed to the ELSI library
    this%nSpin = min(nSpin, 2)
    this%iSpin = iSpin
    if (any(this%iSolver == [electronicSolverTypes%ntpoly, electronicSolverTypes%omm])) then
      if (nSpin == 4) then
        this%nSpin = 2
        this%iSpin = 1
        this%nElectron = 2.0_dp * this%nElectron
      end if
    end if
    this%nKPoint = nKPoint
    this%iKPoint = iKPoint
    this%kWeight = kWeight

    this%mpiCommWorld = env%mpi%globalComm%id
    this%myCommWorld = env%mpi%groupComm%id
    this%myBlacsCtx = env%blacs%orbitalGrid%ctxt

    ! assumes row and column sizes the same
    this%BlacsBlockSize = env%blacs%rowBlockSize

    this%csrBlockSize = this%nBasis / env%mpi%groupComm%size

    if (this%csrBlockSize * env%mpi%groupComm%size < this%nBasis) then
      this%csrBlockSize = this%csrBlockSize + 1
    end if

    this%muMpOrder = 0
    if (iDistribFn == fillingTypes%Fermi) then
      this%muBroadenScheme = 1
    else if (iDistribFn == fillingTypes%Gaussian) then
      this%muBroadenScheme = 0
    else if (iDistribFn >= fillingTypes%Methfessel) then
      this%muBroadenScheme = 2
      ! set Meth-Pax order
      this%muMpOrder = iDistribFn - fillingTypes%Methfessel - 1
    else
      call error("Unknown electronic filling type")
    end if

    if (iDistribFn /= fillingTypes%Fermi .and.&
        & any([electronicSolverTypes%pexsi,electronicSolverTypes%ntpoly] == this%iSolver)) then
      call error("This electronic solver can only be used for Fermi function distributed electrons")
    end if

    ! ELPA settings
    this%elpaSolverOption = inp%elpaSolver

    ! OMM settings
    this%ommIter = inp%ommIterationsElpa
    this%ommTolerance = inp%ommTolerance
    this%ommCholesky = inp%ommCholesky

    ! PEXSI settings
    this%pexsiMethod = inp%pexsiMethod
    this%pexsiNPole = inp%pexsiNPole

    if (this%pexsiNPole < 10) then
      call error("Too few PEXSI poles")
    end if
    select case(this%pexsiMethod)
    case(1)
      if (mod(this%pexsiNPole,10) /= 0 .or. this%pexsiNPole > 120) then
        call error("Unsupported number of PEXSI poles for method 1")
      end if
    case(2,3)
      if (mod(this%pexsiNPole,5) /= 0 .or. this%pexsiNPole > 40) then
        write(buffer,"(A,I0)")"Unsupported number of PEXSI poles for method ", this%pexsiMethod
        call error(trim(buffer))
      end if
    end select

    this%pexsiNpPerPole = inp%pexsiNpPerPole
    this%pexsiNMu = inp%pexsiNMu
    this%pexsiNpSymbo = inp%pexsiNpSymbo
    this%pexsiDeltaE = inp%pexsiDeltaE

    ! NTPoly settings
    this%ntpolyMethod = inp%ntpolyMethod
    this%ntpolyTruncation = inp%ntpolyTruncation
    this%ntpolyTolerance = inp%ntpolyTolerance

    this%isSparse = inp%elsiCsr

    this%elsi_zero_def = inp%elsi_zero_def

    ! data as dense BLACS blocks
    if (this%isSparse) then
      ! CSR format
      this%denseBlacs = 2
    else
      this%denseBlacs = 0
    end if

    if (this%isSparse) then
      allocate(this%elsiCsc)
    end if

    ! customize output level, note there are levels 0..3 not DFTB+ 0..2
    this%OutputLevel = 0
  #:call DEBUG_CODE
    this%OutputLevel = 3
  #:endcall DEBUG_CODE

    if (nSpin == 4 .and. this%isSparse) then
      call error("Sparse solver not currently avaible for two component complex hamiltonians")
    end if

    this%tCholeskyDecomposed = .false.

    if (this%pexsiMethod == 2) then
      providesElectronEntropy = .false.
    end if

  #:else

    call error("Internal error: TElsiSolver_init() called despite missing ELSI support")

  #:endif

  end subroutine TElsiSolver_init


  !> Checks for supported ELSI api version, ideally 2.6.0 or 2.6.1 (correct electronic entropy
  !> return and # PEXSI poles change between 2.5.0 and 2.6.0), but 2.5.0 can also be used with
  !> warnings.
  subroutine supportedVersionNumber(this, dateStamp)

    !> Version value components inside the structure
    class(TElsiSolver), intent(in) :: this

    !> git commit date for the library
    integer, intent(in) :: dateStamp

    logical :: isSupported, isPartSupported

    isSupported = .true.
    if (this%major < 2) then
      isSupported = .false.
    elseif (this%major == 2 .and. this%minor < 6) then
      isSupported = .false.
    elseif (this%major == 2 .and. this%minor == 6 .and. all(this%patch /= [0,1])) then
      ! library must be 2.6.{0,1}
      isSupported = .false.
    end if

    isPartSupported = isSupported
    if (.not.isPartSupported) then
      if (all([this%major, this%minor, this%patch] == [2,5,0])) then
        isPartSupported = .true.
      end if
    end if

    if (.not. (isSupported .or. isPartSupported)) then
      call error("Unsuported ELSI version for DFTB+, requires release >= 2.5.0")
    end if

    if (.not.isSupported .and. isPartSupported) then
      call warning("ELSI version 2.5.0 is only partially supported due to changes in default solver&
          & behaviour for PEXSI at 2.6.0")
    end if

    write(stdOut,"(A,T30,I0,'.',I0,'.',I0)")'ELSI library version :', this%major, this%minor,&
        & this%patch

    if (all([this%major,this%minor,this%patch] == [2,5,0]) .and. dateStamp /= 20200204) then
      call warning("ELSI 2.5.0 library version is between releases")
    end if
    if (all([this%major,this%minor,this%patch] == [2,6,0]) .and. dateStamp /= 20200617) then
      call warning("ELSI 2.6.0 library version is between releases")
    end if
    if (all([this%major,this%minor,this%patch] == [2,6,1]) .and. dateStamp /= 20200625) then
      call warning("ELSI 2.6.1 library version is between releases")
    end if


  end subroutine supportedVersionNumber


  !> Finalizes the ELSI solver.
  subroutine TELsiSolver_final(this)

    !> Instance
    type(TElsiSolver), intent(inout) :: this

  #:if WITH_ELSI

    call elsi_finalize(this%handle)

  #:else

    call error("Internal error: TELsiSolver_final() called despite missing ELSI&
        & support")

  #:endif

  end subroutine TELsiSolver_final


  !> reset the ELSI solver - safer to do this on geometry change, due to the lack of a Choleskii
  !> refactorization option
  subroutine TElsiSolver_reset(this)

    !> Instance
    class(TElsiSolver), intent(inout) :: this


  #:if WITH_ELSI

    if (this%tSolverInitialised) then

      ! reset previous instance of solver
      call elsi_reinit(this%handle)

      if (this%iSolver == electronicSolverTypes%pexsi) then
        ! reset PEXSI chemical potential search
        this%pexsiMuMin = -10.0_dp
        this%pexsiMuMax = 10.0_dp
        this%pexsiDeltaVMin = 0.0_dp
        this%pexsiDeltaVMax = 0.0_dp
      end if

    else

      ! initialise solver

      call elsi_init(this%handle, this%solver, this%parallel, this%denseBlacs, this%nBasis,&
          & this%nElectron, this%nState)

      call elsi_set_mpi_global(this%handle, this%mpiCommWorld)
      call elsi_set_sing_check(this%handle, 0) ! disable singularity check
      call elsi_set_mpi(this%handle, this%myCommWorld)

      if (this%isSparse) then
        call elsi_set_csc_blk(this%handle, this%csrBlockSize)
      else
        call elsi_set_zero_def(this%handle, this%elsi_zero_def)
        ! sparsity of both H and S used as the pattern
        call elsi_set_sparsity_mask(this%handle, 0)
      end if
      call elsi_set_blacs(this%handle, this%myBlacsCtx, this%BlacsBlockSize)

      if (this%tWriteHS) then
        ! setup to write a matrix
        call elsi_init_rw(this%rwHandle, 1, this%parallel, this%nBasis, this%nElectron)
        ! MPI comm
        call elsi_set_rw_mpi(this%rwHandle, this%mpiCommWorld)
        if (.not.this%isSparse) then
          ! dense matrices
          call elsi_set_rw_blacs(this%rwHandle, this%myBlacsCtx, this%BlacsBlockSize)
          ! use default of 1E-15 for the moment, but could be over-ridden if needed, probably with a
          ! different variable name:
          !call elsi_set_rw_zero_def(this%rwHandle, this%elsi_zero_def)
        end if
      end if


      select case(this%iSolver)
      case(electronicSolverTypes%elpa, electronicSolverTypes%elpadm)

        select case(this%elpaSolverOption)
        case(1)
          ! single stage
          call elsi_set_elpa_solver(this%handle, 1)
        case(2)
          ! two stage
          call elsi_set_elpa_solver(this%handle, 2)
        case default
          call error("Unknown ELPA solver modes")
        end select

      case(electronicSolverTypes%omm)
        ! libOMM
        if (this%ommCholesky) then
          call elsi_set_omm_flavor(this%handle, 2)
        else
          call elsi_set_omm_flavor(this%handle, 0)
        end if
        call elsi_set_omm_n_elpa(this%handle, this%ommIter)
        call elsi_set_omm_tol(this%handle, this%ommTolerance)

      case(electronicSolverTypes%pexsi)

        this%pexsiMuMin = -10.0_dp
        this%pexsiMuMax = 10.0_dp
        this%pexsiDeltaVMin = 0.0_dp
        this%pexsiDeltaVMax = 0.0_dp

        ! processors per pole to invert for
        call elsi_set_pexsi_np_per_pole(this%handle, this%pexsiNpPerPole)

        call elsi_set_pexsi_mu_min(this%handle, this%pexsiMuMin +this%pexsiDeltaVMin)
        call elsi_set_pexsi_mu_max(this%handle, this%pexsiMuMax +this%pexsiDeltaVMax)

        ! set the pole expansion method for PEXSI
        call elsi_set_pexsi_method(this%handle, this%pexsiMethod)

        ! number of poles for the expansion
        call elsi_set_pexsi_n_pole(this%handle, this%pexsiNPole)

        ! number of interpolation points for mu
        call elsi_set_pexsi_n_mu(this%handle, this%pexsiNMu)

        ! number of processors for symbolic factorisation task
        call elsi_set_pexsi_np_symbo(this%handle, this%pexsiNpSymbo)

        ! spectral radius (range of eigenspectrum, if known, otherwise default usually fine)
        call elsi_set_pexsi_delta_e(this%handle, this%pexsiDeltaE)

      case(electronicSolverTypes%ntpoly)

        ! NTPoly
        ! set purification method
        call elsi_set_ntpoly_method(this%handle, this%ntpolyMethod)

        ! set truncation tolerance for sparse matrix multiplications
        call elsi_set_ntpoly_filter(this%handle, this%ntpolyTruncation)

        ! set purification convergence threshold
        call elsi_set_ntpoly_tol(this%handle, this%ntpolyTolerance)

      end select

      if (any(this%iSolver == [electronicSolverTypes%omm, electronicSolverTypes%pexsi,&
          & electronicSolverTypes%ntpoly, electronicSolverTypes%elpadm])) then
        ! density matrix build needs to know the number of spin channels to normalize against
        select case(this%nSpin)
        case(1)
          call elsi_set_spin(this%handle, 1, this%iSpin)
        case(2)
          call elsi_set_spin(this%handle, 2, this%iSpin)
        case(4)
          call elsi_set_spin(this%handle, 2, this%iSpin)
        end select
        if (this%nKPoint > 1) then
          call elsi_set_kpoint(this%handle, this%nKPoint, this%iKPoint, this%kWeight)
        end if
      end if

      call elsi_set_output(this%handle, this%OutputLevel)
      if (this%OutputLevel == 3) then
        call elsi_set_output_log(this%handle, 1)
      end if
      this%tSolverInitialised = .true.
    end if

    this%tCholeskyDecomposed = .false.
    this%tFirstCalc = .true.

  #:else

    call error("Internal error: TElsiSolver_reset() called despite missing ELSI support")

  #:endif

  end subroutine TElsiSolver_reset


  !> Resets solver due to geometry changes
  subroutine TElsiSolver_updateGeometry(this, env, neighList, nNeighbourSK, iAtomStart,&
      & iSparseStart, img2CentCell)

    !> Instance
    class(TElsiSolver), intent(inout) :: this

    !> Environment settings
    type(TEnvironment), intent(inout) :: env

    !> Neighbour list
    type(TNeighbourList), intent(in) :: neighList

    !> Number of neighbours for each of the atoms
    integer, intent(in) :: nNeighbourSK(:)

    !> Atom offset for the squared matrix
    integer, intent(in) :: iAtomStart(:)

    !> indexing array for the sparse Hamiltonian
    integer, intent(in) :: iSparseStart(0:,:)

    !> Mapping between image atoms and corresponding atom in the central cell.
    integer, intent(in) :: img2CentCell(:)

  #:if WITH_ELSI

    call TElsiCsc_init(this%elsiCsc, env, this%nBasis, this%csrBlockSize, neighList,&
        & nNeighbourSK, iAtomStart, iSparseStart, img2CentCell)

  #:else

    call error("Internal error: TELsiSolver_updateGeometry() called despite missing ELSI&
        & support")

  #:endif

  end subroutine TElsiSolver_updateGeometry


  !> Updated the electronic temperature
  subroutine TElsiSolver_updateElectronicTemp(this, tempElec)

    !> Instance
    class(TElsiSolver), intent(inout) :: this

    !> Electronic temperature
    real(dp), intent(in) :: tempElec

  #:if WITH_ELSI

    call elsi_set_mu_broaden_scheme(this%handle, this%muBroadenScheme)
    if (this%muBroadenScheme == 2) then
      call elsi_set_mu_mp_order(this%handle, this%muMpOrder)
    end if
    call elsi_set_mu_broaden_width(this%handle, tempElec)
    if (this%iSolver == electronicSolverTypes%pexsi) then
      call elsi_set_pexsi_temp(this%handle, tempElec)
    end if

  #:else

    call error("Internal error: TELsiSolver_updateElectronicTemp() called despite missing ELSI&
        & support")

  #:endif

  end subroutine TElsiSolver_updateElectronicTemp


  !> Returns the density matrix using ELSI non-diagonalisation routines.
  subroutine TElsiSolver_getDensity(this, env, denseDesc, ham, over, neighbourList, nNeighbourSK,&
      & iSparseStart, img2CentCell, iCellVec, cellVec, kPoint, kWeight, tHelical, orb, species,&
      & coord, tRealHS, tSpinSharedEf, tSpinOrbit, tDualSpinOrbit, tMulliken, parallelKS, Ef,&
      & energy, rhoPrim, Eband, TS, iHam, xi, orbitalL, HSqrReal, SSqrReal, iRhoPrim, HSqrCplx,&
      & SSqrCplx)

    !> Electronic solver information
    class(TElsiSolver), intent(inout) :: this

    !> Environment settings
    type(TEnvironment), intent(inout) :: env

    !> Dense matrix descriptor
    type(TDenseDescr), intent(in) :: denseDesc

    !> hamiltonian in sparse storage
    real(dp), intent(in) :: ham(:,:)

    !> sparse overlap matrix
    real(dp), intent(in) :: over(:)

    !> list of neighbours for each atom
    type(TNeighbourList), intent(in) :: neighbourList

    !> Number of neighbours for each of the atoms
    integer, intent(in) :: nNeighbourSK(:)

    !> Index array for the start of atomic blocks in sparse arrays
    integer, intent(in) :: iSparseStart(:,:)

    !> map from image atoms to the original unique atom
    integer, intent(in) :: img2CentCell(:)

    !> Index for which unit cell atoms are associated with
    integer, intent(in) :: iCellVec(:)

    !> Vectors (in units of the lattice constants) to cells of the lattice
    real(dp), intent(in) :: cellVec(:,:)

    !> k-points
    real(dp), intent(in) :: kPoint(:,:)

    !> Weights for k-points
    real(dp), intent(in) :: kWeight(:)

    !> Is the geometry helical
    logical, intent(in) :: tHelical

    !> Atomic orbital information
    type(TOrbitals), intent(in) :: orb

    !> species of all atoms in the system
    integer, intent(in) :: species(:)

    !> atomic coordinates
    real(dp), intent(in) :: coord(:,:)

    !> Is the hamiltonian real (no k-points/molecule/gamma point)?
    logical, intent(in) :: tRealHS

    !> Is the Fermi level common accross spin channels?
    logical, intent(in) :: tSpinSharedEf

    !> Are spin orbit interactions present
    logical, intent(in) :: tSpinOrbit

    !> Are block population spin orbit interactions present
    logical, intent(in) :: tDualSpinOrbit

    !> Should Mulliken populations be generated/output
    logical, intent(in) :: tMulliken

    !> K-points and spins to process
    type(TParallelKS), intent(in) :: parallelKS

    !> Fermi level(s)
    real(dp), intent(inout) :: Ef(:)

    !> Energy contributions and total
    type(TEnergies), intent(inout) :: energy

    !> sparse density matrix
    real(dp), intent(out) :: rhoPrim(:,:)

    !> band structure energy
    real(dp), intent(out) :: Eband(:)

    !> electronic entropy times temperature
    real(dp), intent(out) :: TS(:)

    !> imaginary part of hamiltonian
    real(dp), intent(in), allocatable :: iHam(:,:)

    !> spin orbit constants
    real(dp), intent(in), allocatable :: xi(:,:)

    !> orbital moments of atomic shells
    real(dp), intent(inout), allocatable :: orbitalL(:,:,:)

    !> imaginary part of density matrix
    real(dp), intent(inout), allocatable :: iRhoPrim(:,:)

    !> dense real hamiltonian storage
    real(dp), intent(inout), allocatable :: HSqrReal(:,:)

    !> dense real overlap storage
    real(dp), intent(inout), allocatable :: SSqrReal(:,:)

    !> dense complex (k-points) hamiltonian storage
    complex(dp), intent(inout), allocatable :: HSqrCplx(:,:)

    !> dense complex (k-points) overlap storage
    complex(dp), intent(inout), allocatable :: SSqrCplx(:,:)

  #:if WITH_ELSI

    integer :: nSpin
    integer :: iKS, iK, iS

    ! as each spin and k-point combination forms a separate group for this solver, then iKS = 1
    iKS = 1
    iK = parallelKS%localKS(1, iKS)
    iS = parallelKS%localKS(2, iKS)

    nSpin = size(ham, dim=2)
    if (nSpin == 2 .and. .not. tSpinSharedEf) then
      call error("ELSI currently requires shared Fermi levels over spin")
    end if

    rhoPrim(:,:) = 0.0_dp
    if (allocated(iRhoPrim)) then
      iRhoPrim(:,:) = 0.0_dp
    end if

    if (this%iSolver == electronicSolverTypes%pexsi) then
      call elsi_set_pexsi_mu_min(this%handle, this%pexsiMuMin + this%pexsiDeltaVMin)
      call elsi_set_pexsi_mu_max(this%handle, this%pexsiMuMax + this%pexsiDeltaVMax)
    end if

    if (nSpin /= 4) then
      if (tRealHS) then
        if (this%isSparse) then
          call getDensityRealSparse(this, parallelKS, ham, over, neighbourList%iNeighbour,&
              & nNeighbourSK, denseDesc%iAtomStart, iSparseStart, img2CentCell, tHelical, orb,&
              & species, coord, rhoPrim, Eband)
        else
          call getDensityRealDense(this, env, denseDesc, ham, over, neighbourList,&
              & nNeighbourSK, iSparseStart, img2CentCell, tHelical, orb, species, coord,&
              & parallelKS, rhoPrim, Eband, HSqrReal, SSqrReal)
        end if
      else
        if (this%isSparse) then
          call getDensityCmplxSparse(this, parallelKS, kPoint(:,iK), kWeight(iK), iCellVec,&
              & cellVec, ham, over, neighbourList%iNeighbour, nNeighbourSK, denseDesc%iAtomStart,&
              & iSparseStart, img2CentCell, tHelical, orb, species, coord, rhoPrim, Eband)
        else
          call getDensityCmplxDense(this, env, denseDesc, ham, over, neighbourList,&
              & nNeighbourSK, iSparseStart, img2CentCell, iCellVec, cellVec, kPoint, kWeight,&
              & tHelical, orb, species, coord, parallelKS, rhoPrim, Eband, HSqrCplx, SSqrCplx)
        end if
      end if
      call ud2qm(rhoPrim)
    else
      if (this%isSparse) then
        call error("Internal error: getDensityFromElsiSolver : sparse pauli not yet implemented")
      else
        call getDensityPauliDense(this, env, denseDesc, ham, over, neighbourList,&
            & nNeighbourSK, iSparseStart, img2CentCell, iCellVec, cellVec, kPoint, kWeight, orb,&
            & species, tSpinOrbit, tDualSpinOrbit, tMulliken, parallelKS, energy, rhoPrim, Eband,&
            & iHam, xi, orbitalL, iRhoPrim, HSqrCplx, SSqrCplx)
      end if
    end if

    TS(:) = 0.0_dp
    if (any(this%iSolver == [electronicSolverTypes%pexsi, electronicSolverTypes%elpadm])) then
      call elsi_get_entropy(this%handle, TS(iS))
    end if

    Ef(:) = 0.0_dp
    if (any(this%iSolver == [electronicSolverTypes%pexsi, electronicSolverTypes%ntpoly,&
        & electronicSolverTypes%elpadm])) then
      call elsi_get_mu(this%handle, Ef(iS))
    end if

    if (this%iSolver == electronicSolverTypes%pexsi) then
      call elsi_get_pexsi_mu_min(this%handle, this%pexsiMuMin)
      call elsi_get_pexsi_mu_max(this%handle, this%pexsiMuMax)
    end if

    ! Add up and distribute density matrix contribution from each group
    call mpifx_allreduceip(env%mpi%globalComm, rhoPrim, MPI_SUM)
    if (allocated(iRhoPrim)) then
      call mpifx_allreduceip(env%mpi%globalComm, iRhoPrim, MPI_SUM)
    end if

    if (tSpinOrbit) then
      call mpifx_allreduceip(env%mpi%globalComm, energy%atomLS, MPI_SUM)
      energy%atomLS(:) = 2.0_dp * energy%atomLS
    end if
    if (tMulliken .and. tSpinOrbit .and. .not. tDualSpinOrbit) then
      call mpifx_allreduceip(env%mpi%globalComm, orbitalL, MPI_SUM)
      orbitalL(:,:,:) = 2.0_dp * orbitalL
    end if
    if (tSpinOrbit .and. .not. tDualSpinOrbit) then
      energy%ELS = sum(energy%atomLS)
    end if

    call env%globalTimer%stopTimer(globalTimers%densityMatrix)

  #:else

    call error("Internal error: TELsiSolver_getDensity() called despite missing ELSI support")

  #:endif

  end subroutine TElsiSolver_getDensity


  ! Returns the energy weighted density matrix using ELSI non-diagonalisation routines.
  subroutine TElsiSolver_getEDensity(this, env, denseDesc, nSpin, kPoint, kWeight, neighbourList,&
      & nNeighbourSK, tHelical, orb, species, coord, iSparseStart, img2CentCell, iCellVec, cellVec,&
      & tRealHS, parallelKS, ERhoPrim, SSqrReal, SSqrCplx)

    !> Electronic solver information
    class(TElsiSolver), intent(inout) :: this

    !> Environment settings
    type(TEnvironment), intent(inout) :: env

    !> Dense matrix descriptor
    type(TDenseDescr), intent(in) :: denseDesc

    !> Number of spin channels
    integer, intent(in) :: nSpin

    !> K-points
    real(dp), intent(in) :: kPoint(:,:)

    !> Weights for k-points
    real(dp), intent(in) :: kWeight(:)

    !> list of neighbours for each atom
    type(TNeighbourList), intent(in) :: neighbourList

    !> Number of neighbours for each of the atoms
    integer, intent(in) :: nNeighbourSK(:)

    !> Is the geometry helical
    logical, intent(in) :: tHelical

    !> Atomic orbital information
    type(TOrbitals), intent(in) :: orb

    !> species of all atoms in the system
    integer, intent(in) :: species(:)

    !> atomic coordinates
    real(dp), intent(in) :: coord(:,:)

    !> Index array for the start of atomic blocks in sparse arrays
    integer, intent(in) :: iSparseStart(:,:)

    !> map from image atoms to the original unique atom
    integer, intent(in) :: img2CentCell(:)

    !> Index for which unit cell atoms are associated with
    integer, intent(in) :: iCellVec(:)

    !> Vectors (in units of the lattice constants) to cells of the lattice
    real(dp), intent(in) :: cellVec(:,:)

    !> Is the hamiltonian real (no k-points/molecule/gamma point)?
    logical, intent(in) :: tRealHS

    !> K-points and spins to process
    type(TParallelKS), intent(in) :: parallelKS

    !> Energy weighted sparse matrix
    real(dp), intent(out) :: ERhoPrim(:)

    !> Storage for dense overlap matrix
    real(dp), intent(inout), allocatable :: SSqrReal(:,:)

    !> Storage for dense overlap matrix (complex case)
    complex(dp), intent(inout), allocatable :: SSqrCplx(:,:)

  #:if WITH_ELSI

    if (nSpin == 4) then
      call getEDensityMtxPauli(this, env, denseDesc, kPoint, kWeight,&
          & neighbourList, nNeighbourSK, orb, iSparseStart, img2CentCell, iCellVec, cellVec,&
          & parallelKS, ERhoPrim, SSqrCplx)
    else
      if (tRealHS) then
        call getEDensityMtxReal(this, env, denseDesc, neighbourList, nNeighbourSK, tHelical, orb,&
            & species, coord, iSparseStart, img2CentCell, ERhoPrim, SSqrReal)
      else
        call getEDensityMtxCmplx(this, env, denseDesc, kPoint, kWeight, neighbourList,&
            & nNeighbourSK, tHelical, orb, species, coord, iSparseStart, img2CentCell, iCellVec,&
            & cellVec, parallelKS, ERhoPrim, SSqrCplx)
      end if
    end if

  #:else

    call error("Internal error: TELsiSolver_getDensity() called despite missing ELSI support")

  #:endif

  end subroutine TElsiSolver_getEDensity


  !> Initializes the PEXSI potentials.
  subroutine TElsiSolver_initPexsiDeltaVRanges(this, tSccCalc, potential)

    !> Electronic solver information
    class(TElsiSolver), intent(inout) :: this

    !> Whether we have an SCC calculation
    logical, intent(in) :: tSccCalc

    !> Potentials acting
    type(TPotentials), intent(in) :: potential

  #:if WITH_ELSI

    if (tSccCalc) then
      if (.not.allocated(this%pexsiVOld)) then
        allocate(this%pexsiVOld(size(potential%intBlock(:,:,:,1))))
      end if
      this%pexsiVOld = reshape(potential%intBlock(:,:,:,1) + potential%extBlock(:,:,:,1),&
          & [size(potential%extBlock(:,:,:,1))])
    end if
    this%pexsiDeltaVMin = 0.0_dp
    this%pexsiDeltaVMax = 0.0_dp

  #:else

    call error("Internal error: TELsiSolver_initPexsiDeltaVRanges() called despite missing ELSI&
        & support")

  #:endif

  end subroutine TElsiSolver_initPexsiDeltaVRanges


  !> Updates the PEXSI potentials.
  subroutine TElsiSolver_updatePexsiDeltaVRanges(this, potential)

    !> Electronic solver information
    class(TElsiSolver), intent(inout) :: this

    !> Potentials acting
    type(TPotentials), intent(in) :: potential

  #:if WITH_ELSI

    this%pexsiDeltaVMin = minval(this%pexsiVOld&
        & - reshape(potential%intBlock(:,:,:,1) + potential%extBlock(:,:,:,1),&
        & [size(potential%extBlock(:,:,:,1))]))
    this%pexsiDeltaVMax = maxval(this%pexsiVOld&
        & - reshape(potential%intBlock(:,:,:,1) + potential%extBlock(:,:,:,1),&
        & [size(potential%extBlock(:,:,:,1))]))
    this%pexsiVOld = reshape(potential%intBlock(:,:,:,1) + potential%extBlock(:,:,:,1),&
        & [size(potential%extBlock(:,:,:,1))])

  #:else

    call error("Internal error: TElsiSolver_updatePexsiDeltaVRanges() called despite missing ELSI&
        & support")

  #:endif

  end subroutine TElsiSolver_updatePexsiDeltaVRanges


  !> Returns the name of the solver
  function TElsiSolver_getSolverName(this) result(solverName)

    !> Instance.
    class(TElsiSolver), intent(in) :: this

    !> Name of the solver.
    character(:), allocatable :: solverName

  #:if WITH_ELSI

    character(lc) :: buffer

    select case (this%iSolver)

    case(electronicSolverTypes%elpa)
      if (this%elpaSolverOption == 1) then
        write(buffer, "(A)") "ELSI interface to the 1 stage ELPA solver"
      else
        write(buffer, "(A)") "ELSI interface to the 2 stage ELPA solver"
      end if

    case(electronicSolverTypes%omm)
      write(buffer, "(A,I0,A,E9.2)") "ELSI solver libOMM with ",&
          & this%ommIter, " ELPA iterations",this%ommTolerance
      if (this%isSparse) then
        write(buffer, "(A)") "ELSI solver libOMM Sparse"
      else
        write(buffer, "(A)") "ELSI solver libOMM Dense"
      end if

    case(electronicSolverTypes%pexsi)
      if (this%isSparse) then
        write(buffer, "(A)") "ELSI solver PEXSI Sparse"
      else
        write(buffer, "(A)") "ELSI solver PEXSI Dense"
      end if

    case(electronicSolverTypes%ntpoly)
      if (this%isSparse) then
        write(buffer, "(A)") "ELSI solver NTPoly Sparse"
      else
        write(buffer, "(A)") "ELSI solver NTPoly Dense"
      end if

    case default
      write(buffer, "(A,I0,A)") "Invalid electronic solver! (iSolver = ", this%iSolver, ")"

    end select

    solverName = trim(buffer)

  #:else

    solverName = ""

    call error("Internal error: TElsiSolver_getSolverName() called despite missing ELSI support")

  #:endif

  end function TElsiSolver_getSolverName


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!  Private routines
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

#:if WITH_ELSI

  !> Returns the density matrix using ELSI non-diagonalisation routines (real dense case).
  subroutine getDensityRealDense(this, env, denseDesc, ham, over, neighbourList, nNeighbourSK,&
      & iSparseStart, img2CentCell, tHelical, orb, species, coord, parallelKS, rhoPrim, Eband,&
      & HSqrReal, SSqrReal)

    !> Electronic solver information
    type(TElsiSolver), intent(inout) :: this

    !> Environment settings
    type(TEnvironment), intent(inout) :: env

    !> Dense matrix descriptor
    type(TDenseDescr), intent(in) :: denseDesc

    !> hamiltonian in sparse storage
    real(dp), intent(in) :: ham(:,:)

    !> sparse overlap matrix
    real(dp), intent(in) :: over(:)

    !> list of neighbours for each atom
    type(TNeighbourList), intent(in) :: neighbourList

    !> Number of neighbours for each of the atoms
    integer, intent(in) :: nNeighbourSK(:)

    !> Index array for the start of atomic blocks in sparse arrays
    integer, intent(in) :: iSparseStart(:,:)

    !> map from image atoms to the original unique atom
    integer, intent(in) :: img2CentCell(:)

    !> Is the geometry helicalxs
    logical, intent(in) :: tHelical

    !> Atomic orbital information
    type(TOrbitals), intent(in) :: orb

    !> species of all atoms in the system
    integer, intent(in) :: species(:)

    !> atomic coordinates
    real(dp), intent(in) :: coord(:,:)

    !> K-points and spins to process
    type(TParallelKS), intent(in) :: parallelKS

    !> sparse density matrix
    real(dp), intent(inout) :: rhoPrim(:,:)

    !> band structure energy
    real(dp), intent(out) :: Eband(:)

    !> dense real hamiltonian storage
    real(dp), intent(inout), allocatable :: HSqrReal(:,:)

    !> dense real overlap storage
    real(dp), intent(inout), allocatable :: SSqrReal(:,:)

    real(dp), allocatable :: rhoSqrReal(:,:)
    integer :: iKS, iS

    ! as each spin and k-point combination forms a separate group for this solver, then iKS = 1
    iKS = 1
    iS = parallelKS%localKS(2, iKS)

    if (tHelical) then
      call unpackHSHelicalRealBlacs(env%blacs, ham(:,iS), neighbourList%iNeighbour, nNeighbourSK,&
          & iSparseStart, img2CentCell, orb, species, coord, denseDesc, HSqrReal)
      if (.not. this%tCholeskyDecomposed) then
        call unpackHSHelicalRealBlacs(env%blacs, over, neighbourList%iNeighbour, nNeighbourSK,&
            & iSparseStart, img2CentCell, orb, species, coord, denseDesc, SSqrReal)
      end if
    else
      call unpackHSRealBlacs(env%blacs, ham(:,iS), neighbourList%iNeighbour, nNeighbourSK,&
          & iSparseStart, img2CentCell, denseDesc, HSqrReal)
      if (.not. this%tCholeskyDecomposed) then
        call unpackHSRealBlacs(env%blacs, over, neighbourList%iNeighbour, nNeighbourSK,&
            & iSparseStart, img2CentCell, denseDesc, SSqrReal)
      end if
    end if

    if (.not. this%tCholeskyDecomposed) then
      if (this%ommCholesky .or. this%iSolver == electronicSolverTypes%pexsi) then
        this%tCholeskyDecomposed = .true.
      end if
    end if

    allocate(rhoSqrReal(size(HSqrReal, dim=1), size(HSqrReal, dim=2)))
    rhoSqrReal(:,:) = 0.0_dp
    Eband(iS) = 0.0_dp
    if (this%tWriteHS) then
      call elsi_write_mat_real(this%rwHandle, "ELSI_Hreal.bin", HSqrReal)
      call elsi_write_mat_real(this%rwHandle, "ELSI_Sreal.bin", SSqrReal)
      call elsi_finalize_rw(this%rwHandle)
      call cleanShutdown("Finished dense matrix write")
    end if
    call elsi_dm_real(this%handle, HSqrReal, SSqrReal, rhoSqrReal, Eband(iS))

    if (tHelical) then
      call packRhoHelicalRealBlacs(env%blacs, denseDesc, rhoSqrReal, neighbourList%iNeighbour,&
          & nNeighbourSK, iSparseStart, img2CentCell, orb, species, coord, rhoPrim(:,iS))
    else
      call packRhoRealBlacs(env%blacs, denseDesc, rhoSqrReal, neighbourList%iNeighbour,&
          & nNeighbourSK, orb%mOrb, iSparseStart, img2CentCell, rhoPrim(:,iS))
    end if

  end subroutine getDensityRealDense


  !> Returns the density matrix using ELSI non-diagonalisation routines (complex dense case).
  subroutine getDensityCmplxDense(this, env, denseDesc, ham, over, neighbourList,&
      & nNeighbourSK, iSparseStart, img2CentCell, iCellVec, cellVec, kPoint, kWeight, tHelical,&
      & orb, species, coord, parallelKS, rhoPrim, Eband, HSqrCplx, SSqrCplx)

    !> Electronic solver information
    type(TElsiSolver), intent(inout) :: this

    !> Environment settings
    type(TEnvironment), intent(inout) :: env

    !> Dense matrix descriptor
    type(TDenseDescr), intent(in) :: denseDesc

    !> hamiltonian in sparse storage
    real(dp), intent(in) :: ham(:,:)

    !> sparse overlap matrix
    real(dp), intent(in) :: over(:)

    !> list of neighbours for each atom
    type(TNeighbourList), intent(in) :: neighbourList

    !> Number of neighbours for each of the atoms
    integer, intent(in) :: nNeighbourSK(:)

    !> Index array for the start of atomic blocks in sparse arrays
    integer, intent(in) :: iSparseStart(:,:)

    !> map from image atoms to the original unique atom
    integer, intent(in) :: img2CentCell(:)

    !> Index for which unit cell atoms are associated with
    integer, intent(in) :: iCellVec(:)

    !> Vectors (in units of the lattice constants) to cells of the lattice
    real(dp), intent(in) :: cellVec(:,:)

    !> k-points
    real(dp), intent(in) :: kPoint(:,:)

    !> Weights for k-points
    real(dp), intent(in) :: kWeight(:)

    !> Is the geometry helical
    logical, intent(in) :: tHelical

    !> Atomic orbital information
    type(TOrbitals), intent(in) :: orb

    !> species of all atoms in the system
    integer, intent(in) :: species(:)

    !> atomic coordinates
    real(dp), intent(in) :: coord(:,:)

    !> K-points and spins to process
    type(TParallelKS), intent(in) :: parallelKS

    !> sparse density matrix
    real(dp), intent(inout) :: rhoPrim(:,:)

    !> band structure energy
    real(dp), intent(out) :: Eband(:)

    !> electronic entropy times temperature
    !> dense complex (k-points) hamiltonian storage
    complex(dp), intent(inout), allocatable :: HSqrCplx(:,:)

    !> dense complex (k-points) overlap storage
    complex(dp), intent(inout), allocatable :: SSqrCplx(:,:)

    complex(dp), allocatable :: rhoSqrCplx(:,:)
    integer :: iKS, iK, iS

    ! as each spin and k-point combination forms a separate group for this solver, then iKS = 1
    iKS = 1
    iK = parallelKS%localKS(1, iKS)
    iS = parallelKS%localKS(2, iKS)

    HSqrCplx(:,:) = 0.0_dp
    if (tHelical) then
      call unpackHSHelicalCplxBlacs(env%blacs, ham(:,iS), kPoint(:,iK), neighbourList%iNeighbour,&
          & nNeighbourSK, iCellVec, cellVec, iSparseStart, img2CentCell, orb, species, coord,&
          & denseDesc, HSqrCplx)
      if (.not. this%tCholeskyDecomposed) then
        SSqrCplx(:,:) = 0.0_dp
        call unpackHSHelicalCplxBlacs(env%blacs, over, kPoint(:,iK), neighbourList%iNeighbour,&
            & nNeighbourSK, iCellVec, cellVec, iSparseStart, img2CentCell, orb, species, coord,&
            & denseDesc, SSqrCplx)
      end if
    else
      call unpackHSCplxBlacs(env%blacs, ham(:,iS), kPoint(:,iK), neighbourList%iNeighbour,&
          & nNeighbourSK, iCellVec, cellVec, iSparseStart, img2CentCell, denseDesc, HSqrCplx)
      if (.not. this%tCholeskyDecomposed) then
        SSqrCplx(:,:) = 0.0_dp
        call unpackHSCplxBlacs(env%blacs, over, kPoint(:,iK), neighbourList%iNeighbour,&
            & nNeighbourSK, iCellVec, cellVec, iSparseStart, img2CentCell, denseDesc,&
            & SSqrCplx)
      end if
    end if
    if (.not. this%tCholeskyDecomposed) then
      if (this%ommCholesky .or. this%iSolver == electronicSolverTypes%pexsi) then
        this%tCholeskyDecomposed = .true.
      end if
    end if

    allocate(rhoSqrCplx(size(HSqrCplx,dim=1),size(HSqrCplx,dim=2)))
    rhoSqrCplx(:,:) = 0.0_dp
    if (this%tWriteHS) then
      call elsi_write_mat_complex(this%rwHandle, "ELSI_Hcmplex.bin",&
          & HSqrCplx)
      call elsi_write_mat_complex(this%rwHandle, "ELSI_Scmplx.bin",&
          & SSqrCplx)
      call elsi_finalize_rw(this%rwHandle)
      call cleanShutdown("Finished dense matrix write")
    end if
    call elsi_dm_complex(this%handle, HSqrCplx, SSqrCplx, rhoSqrCplx,&
        & Eband(iS))
    if (tHelical) then
      call packRhoHelicalCplxBlacs(env%blacs, denseDesc, rhoSqrCplx, kPoint(:,iK), kWeight(iK),&
          & neighbourList%iNeighbour, nNeighbourSK, iCellVec, cellVec, iSparseStart, img2CentCell,&
          & orb, species, coord, rhoPrim(:,iS))
    else
      call packRhoCplxBlacs(env%blacs, denseDesc, rhoSqrCplx, kPoint(:,iK), kWeight(iK),&
          & neighbourList%iNeighbour, nNeighbourSK, orb%mOrb, iCellVec, cellVec,&
          & iSparseStart, img2CentCell, rhoPrim(:,iS))
    end if

  end subroutine getDensityCmplxDense


  !> Returns the density matrix using ELSI non-diagonalisation routines (Pauli dense case).
  subroutine getDensityPauliDense(this, env, denseDesc, ham, over, neighbourList,&
      & nNeighbourSK, iSparseStart, img2CentCell, iCellVec, cellVec, kPoint, kWeight, orb, species,&
      & tSpinOrbit, tDualSpinOrbit, tMulliken, parallelKS, energy, rhoPrim, Eband, iHam, xi,&
      & orbitalL, iRhoPrim, HSqrCplx, SSqrCplx)

    !> Electronic solver information
    type(TElsiSolver), intent(inout) :: this

    !> Environment settings
    type(TEnvironment), intent(inout) :: env

    !> Dense matrix descriptor
    type(TDenseDescr), intent(in) :: denseDesc

    !> hamiltonian in sparse storage
    real(dp), intent(in) :: ham(:,:)

    !> sparse overlap matrix
    real(dp), intent(in) :: over(:)

    !> list of neighbours for each atom
    type(TNeighbourList), intent(in) :: neighbourList

    !> Number of neighbours for each of the atoms
    integer, intent(in) :: nNeighbourSK(:)

    !> Index array for the start of atomic blocks in sparse arrays
    integer, intent(in) :: iSparseStart(:,:)

    !> map from image atoms to the original unique atom
    integer, intent(in) :: img2CentCell(:)

    !> Index for which unit cell atoms are associated with
    integer, intent(in) :: iCellVec(:)

    !> Vectors (in units of the lattice constants) to cells of the lattice
    real(dp), intent(in) :: cellVec(:,:)

    !> k-points
    real(dp), intent(in) :: kPoint(:,:)

    !> Weights for k-points
    real(dp), intent(in) :: kWeight(:)

    !> Atomic orbital information
    type(TOrbitals), intent(in) :: orb

    !> species of all atoms in the system
    integer, intent(in) :: species(:)

    !> Are spin orbit interactions present
    logical, intent(in) :: tSpinOrbit

    !> Are block population spin orbit interactions present
    logical, intent(in) :: tDualSpinOrbit

    !> Should Mulliken populations be generated/output
    logical, intent(in) :: tMulliken

    !> K-points and spins to process
    type(TParallelKS), intent(in) :: parallelKS

    !> Energy contributions and total
    type(TEnergies), intent(inout) :: energy

    !> sparse density matrix
    real(dp), intent(inout) :: rhoPrim(:,:)

    !> band structure energy
    real(dp), intent(out) :: Eband(:)

    !> imaginary part of hamiltonian
    real(dp), intent(in), allocatable :: iHam(:,:)

    !> spin orbit constants
    real(dp), intent(in), allocatable :: xi(:,:)

    !> orbital moments of atomic shells
    real(dp), intent(inout), allocatable :: orbitalL(:,:,:)

    !> imaginary part of density matrix
    real(dp), intent(inout), allocatable :: iRhoPrim(:,:)

    !> dense complex (k-points) hamiltonian storage
    complex(dp), intent(inout), allocatable :: HSqrCplx(:,:)

    !> dense complex (k-points) overlap storage
    complex(dp), intent(inout), allocatable :: SSqrCplx(:,:)

    complex(dp), allocatable :: rhoSqrCplx(:,:)
    real(dp), allocatable :: rVecTemp(:), orbitalLPart(:,:,:)
    integer :: iKS, iK, iS
    integer :: nAtom

    ! as each spin and k-point combination forms a separate group for this solver, then iKS = 1
    iKS = 1
    iK = parallelKS%localKS(1, iKS)
    iS = parallelKS%localKS(2, iKS)

    nAtom = size(nNeighbourSK)

    if (tSpinOrbit .and. .not. tDualSpinOrbit) then
      energy%atomLS(:) = 0.0_dp
      allocate(rVecTemp(nAtom))
    end if

    if (tSpinOrbit .and. .not. tDualSpinOrbit) then
      allocate(orbitalLPart(3, orb%mShell, nAtom))
      orbitalL(:,:,:) = 0.0_dp
    end if

    if (allocated(iHam)) then
      call unpackHPauliBlacs(env%blacs, ham, kPoint(:,iK), neighbourList%iNeighbour,&
          & nNeighbourSK, iCellVec, cellVec, iSparseStart, img2CentCell, orb%mOrb, denseDesc,&
          & HSqrCplx, iorig=iHam)
    else
      call unpackHPauliBlacs(env%blacs, ham, kPoint(:,iK), neighbourList%iNeighbour,&
          & nNeighbourSK, iCellVec, cellVec, iSparseStart, img2CentCell, orb%mOrb, denseDesc,&
          & HSqrCplx)
    end if
    if (.not. this%tCholeskyDecomposed) then
      SSqrCplx(:,:) = 0.0_dp
      call unpackSPauliBlacs(env%blacs, over, kPoint(:,iK), neighbourList%iNeighbour,&
          & nNeighbourSK, iCellVec, cellVec, iSparseStart, img2CentCell, orb%mOrb, denseDesc,&
          & SSqrCplx)
      if (this%ommCholesky .or. this%iSolver == electronicSolverTypes%pexsi) then
        this%tCholeskyDecomposed = .true.
      end if
    end if
    if (allocated(xi) .and. .not. allocated(iHam)) then
      call addOnsiteSpinOrbitHam(env, xi, species, orb, denseDesc, HSqrCplx)
    end if
    allocate(rhoSqrCplx(size(HSqrCplx,dim=1),size(HSqrCplx,dim=2)))
    rhoSqrCplx(:,:) = 0.0_dp
    if (this%tWriteHS) then
      call elsi_write_mat_complex(this%rwHandle, "ELSI_Hcmplex.bin",&
          & HSqrCplx)
      call elsi_write_mat_complex(this%rwHandle, "ELSI_Scmplx.bin",&
          & SSqrCplx)
      call elsi_finalize_rw(this%rwHandle)
      call cleanShutdown("Finished dense matrix write")
    end if
    call elsi_dm_complex(this%handle, HSqrCplx, SSqrCplx, rhoSqrCplx,&
        & Eband(iS))
    if (tSpinOrbit .and. .not. tDualSpinOrbit) then
      call getOnsiteSpinOrbitEnergy(env, rVecTemp, rhoSqrCplx, denseDesc, xi, orb, species)
      energy%atomLS = energy%atomLS + kWeight(iK) * rVecTemp
      if (tMulliken) then
        orbitalLPart(:,:,:) = 0.0_dp
        call getLOnsite(env, orbitalLPart, rhoSqrCplx, denseDesc, orb, species)
        orbitalL(:,:,:) = orbitalL + kWeight(iK) * orbitalLPart
      end if
    end if
    if (allocated(iHam)) then
      call packRhoPauliBlacs(env%blacs, denseDesc, rhoSqrCplx, kPoint(:,iK), kWeight(iK),&
          & neighbourList%iNeighbour, nNeighbourSK, orb%mOrb, iCellVec, cellVec,&
          & iSparseStart, img2CentCell, rhoPrim, iRhoPrim)
      iRhoPrim(:,:) = 2.0_dp * iRhoPrim
    else
      call packRhoPauliBlacs(env%blacs, denseDesc, rhoSqrCplx, kPoint(:,iK), kWeight(iK),&
          & neighbourList%iNeighbour, nNeighbourSK, orb%mOrb, iCellVec, cellVec,&
          & iSparseStart, img2CentCell, rhoPrim)
    end if
    RhoPrim(:,:) = 2.0_dp * RhoPrim

  end subroutine getDensityPauliDense


  !> Calculates density matrix using the elsi routine.
  subroutine getDensityRealSparse(this, parallelKS, ham, over, iNeighbour, nNeighbourSK,&
      & iAtomStart, iSparseStart, img2CentCell, tHelical, orb, species, coord, rho, Eband)

    !> Electronic solver information
    type(TElsiSolver), intent(inout) :: this

    !> Contains (iK, iS) tuples to be processed in parallel by various processor groups
    type(TParallelKS), intent(in) :: parallelKS

    !> hamiltonian in sparse storage
    real(dp), intent(in) :: ham(:,:)

    !> overlap matrix in sparse storage
    real(dp), intent(in) :: over(:)

    !> Neighbour list for the atoms (First index from 0!)
    integer, intent(in) :: iNeighbour(0:,:)

    !> Nr. of neighbours for the atoms.
    integer, intent(in) :: nNeighbourSK(:)

    !> Atom offset for the squared matrix
    integer, intent(in) :: iAtomStart(:)

    !> indexing array for the sparse Hamiltonian
    integer, intent(in) :: iSparseStart(0:,:)

    !> Mapping between image atoms and corresponding atom in the central cell.
    integer, intent(in) :: img2CentCell(:)

    !> Is the geometry helical
    logical, intent(in) :: tHelical

    !> data structure with atomic orbital information
    type(TOrbitals), intent(in) :: orb

    !> species of all atoms in the system
    integer, intent(in) :: species(:)

    !> atomic coordinates
    real(dp), intent(in) :: coord(:,:)

    !> Density matrix in DFTB+ sparse format
    real(dp), intent(out) :: rho(:,:)

    !> Band energy
    real(dp), intent(out) :: Eband(:)

    integer :: iS

    real(dp), allocatable :: HnzValLocal(:), SnzValLocal(:)
    real(dp), allocatable :: DMnzValLocal(:)

    allocate(HnzValLocal(this%elsiCsc%nnzLocal))
    allocate(SnzValLocal(this%elsiCsc%nnzLocal))
    allocate(DMnzValLocal(this%elsiCsc%nnzLocal))

    if (tHelical) then
      call this%elsiCsc%convertPackedToElsiReal(over, iNeighbour, nNeighbourSK, iAtomStart,&
          & iSparseStart, img2CentCell, SnzValLocal, orb, species, coord)
    else
      call this%elsiCsc%convertPackedToElsiReal(over, iNeighbour, nNeighbourSK, iAtomStart,&
          & iSparseStart, img2CentCell, SnzValLocal, orb)
    end if

    if (this%tFirstCalc) then
      call elsi_set_csc(this%handle, this%elsiCsc%nnzGlobal, this%elsiCsc%nnzLocal,&
          & this%elsiCsc%numColLocal, this%elsiCsc%rowIndLocal, this%elsiCsc%colPtrLocal)
      this%tFirstCalc = .false.
    end if

    iS = parallelKS%localKS(2, 1)

    if (tHelical) then
      call this%elsiCsc%convertPackedToElsiReal(ham(:,iS), iNeighbour, nNeighbourSK, iAtomStart,&
          & iSparseStart, img2CentCell, HnzValLocal, orb, species, coord)
    else
      call this%elsiCsc%convertPackedToElsiReal(ham(:,iS), iNeighbour, nNeighbourSK, iAtomStart,&
          & iSparseStart, img2CentCell, HnzValLocal, orb)
    end if

    if (this%tWriteHS) then
      call elsi_set_rw_csc(this%rwHandle, this%elsiCsc%nnzGlobal,&
          & this%elsiCsc%nnzLocal, this%elsiCsc%numColLocal)
      call elsi_write_mat_real_sparse(this%rwHandle, "ELSI_HrealSparse.bin",&
          & this%elsiCsc%rowIndLocal, this%elsiCsc%colPtrLocal, HnzValLocal)
      call elsi_write_mat_real_sparse(this%rwHandle, "ELSI_SrealSparse.bin",&
          & this%elsiCsc%rowIndLocal, this%elsiCsc%colPtrLocal, SnzValLocal)
      call elsi_finalize_rw(this%rwHandle)
      call cleanShutdown("Finished matrix write")
    end if

    call elsi_dm_real_sparse(this%handle, HnzValLocal, SnzValLocal, DMnzValLocal,&
        & Eband(iS))

    rho(:,:) = 0.0_dp

    if (tHelical) then
      call this%elsiCsc%convertElsiToPackedReal(iNeighbour, nNeighbourSK, orb, species, coord,&
          & iAtomStart, iSparseStart, img2CentCell, DMnzValLocal, rho(:,iS))
    else
      call this%elsiCsc%convertElsiToPackedReal(iNeighbour, nNeighbourSK, orb, iAtomStart,&
          & iSparseStart, img2CentCell, DMnzValLocal, rho(:,iS))
    end if

  end subroutine getDensityRealSparse


  !> Calculates density matrix using the elsi routine.
  subroutine getDensityCmplxSparse(this, parallelKS, kPoint, kWeight, iCellVec, cellVec, ham,&
      & over, iNeighbour, nNeighbourSK, iAtomStart, iSparseStart, img2CentCell, tHelical, orb,&
      & species, coord, rho, Eband)

    !> Electronic solver information
    type(TElsiSolver), intent(inout) :: this

    !> Contains (iK, iS) tuples to be processed in parallel by various processor groups
    type(TParallelKS), intent(in) :: parallelKS

    !> Current k-point
    real(dp), intent(in) :: kPoint(:)

    !> Weight for current k-points
    real(dp), intent(in) :: kWeight

    !> Index for which unit cell atoms are associated with
    integer, intent(in) :: iCellVec(:)

    !> Vectors (in units of the lattice constants) to cells of the lattice
    real(dp), intent(in) :: cellVec(:,:)

    !> hamiltonian in sparse storage
    real(dp), intent(in) :: ham(:,:)

    !> overlap matrix in sparse storage
    real(dp), intent(in) :: over(:)

    !> Neighbour list for the atoms (First index from 0!)
    integer, intent(in) :: iNeighbour(0:,:)

    !> Nr. of neighbours for the atoms.
    integer, intent(in) :: nNeighbourSK(:)

    !> Atom offset for the squared matrix
    integer, intent(in) :: iAtomStart(:)

    !> indexing array for the sparse Hamiltonian
    integer, intent(in) :: iSparseStart(0:,:)

    !> Mapping between image atoms and corresponding atom in the central cell.
    integer, intent(in) :: img2CentCell(:)

    !> Is the geometry helical
    logical, intent(in) :: tHelical

    !> data structure with atomic orbital information
    type(TOrbitals), intent(in) :: orb

    !> species of all atoms in the system
    integer, intent(in) :: species(:)

    !> atomic coordinates
    real(dp), intent(in) :: coord(:,:)

    !> Density matrix in DFTB+ sparse format
    real(dp), intent(out) :: rho(:,:)

    !> Band energy
    real(dp), intent(out) :: Eband(:)

    integer :: iS

    complex(dp), allocatable :: HnzValLocal(:), SnzValLocal(:)
    complex(dp), allocatable :: DMnzValLocal(:)

    allocate(HnzValLocal(this%elsiCsc%nnzLocal))
    allocate(SnzValLocal(this%elsiCsc%nnzLocal))
    allocate(DMnzValLocal(this%elsiCsc%nnzLocal))

    if (tHelical) then
      call this%elsiCsc%convertPackedToElsiCmplx(over, iNeighbour, nNeighbourSK, iAtomStart,&
          & iSparseStart, img2CentCell, kPoint, iCellVec, cellVec, SnzValLocal, orb, species, coord)
    else
      call this%elsiCsc%convertPackedToElsiCmplx(over, iNeighbour, nNeighbourSK, iAtomStart,&
          & iSparseStart, img2CentCell, kPoint, iCellVec, cellVec, SnzValLocal, orb)
    end if
    if (this%tFirstCalc) then
      call elsi_set_csc(this%handle, this%elsiCsc%nnzGlobal, this%elsiCsc%nnzLocal,&
          & this%elsiCsc%numColLocal, this%elsiCsc%rowIndLocal, this%elsiCsc%colPtrLocal)
      this%tFirstCalc = .false.
    end if

    iS = parallelKS%localKS(2, 1)

    if (tHelical) then
      call this%elsiCsc%convertPackedToElsiCmplx(ham(:,iS), iNeighbour, nNeighbourSK, iAtomStart,&
          & iSparseStart, img2CentCell, kPoint, iCellVec, cellVec, HnzValLocal, orb, species, coord)
    else
      call this%elsiCsc%convertPackedToElsiCmplx(ham(:,iS), iNeighbour, nNeighbourSK, iAtomStart,&
          & iSparseStart, img2CentCell, kPoint, iCellVec, cellVec, HnzValLocal, orb)
    end if

    if (this%tWriteHS) then
      call elsi_set_rw_csc(this%rwHandle, this%elsiCsc%nnzGlobal, this%elsiCsc%nnzLocal,&
          & this%elsiCsc%numColLocal)
      call elsi_write_mat_complex_sparse(this%rwHandle, "ELSI_HcmplxSparse.bin",&
          & this%elsiCsc%rowIndLocal, this%elsiCsc%colPtrLocal, HnzValLocal)
      call elsi_write_mat_complex_sparse(this%rwHandle, "ELSI_ScmplxSparse.bin",&
          & this%elsiCsc%rowIndLocal, this%elsiCsc%colPtrLocal, SnzValLocal)
      call elsi_finalize_rw(this%rwHandle)
      call cleanShutdown("Finished matrix write")
    end if

    call elsi_dm_complex_sparse(this%handle, HnzValLocal, SnzValLocal, DMnzValLocal, Eband(iS))

    rho(:,:) = 0.0_dp
    if (tHelical) then
      call this%elsiCsc%convertElsiToPackedCmplx(iNeighbour, nNeighbourSK, orb, species, coord,&
          & iAtomStart, iSparseStart, img2CentCell, kPoint, kWeight, iCellVec, cellVec,&
          & DMnzValLocal, rho(:,iS))
    else
      call this%elsiCsc%convertElsiToPackedCmplx(iNeighbour, nNeighbourSK, orb, iAtomStart,&
          & iSparseStart, img2CentCell, kPoint, kWeight, iCellVec, cellVec, DMnzValLocal, rho(:,iS))
    end if

  end subroutine getDensityCmplxSparse


  !> Returns the energy weighted density matrix using ELSI non-diagonalisation routines.
  subroutine getEDensityMtxReal(this, env, denseDesc, neighbourList, nNeighbourSK, tHelical, orb,&
      & species, coord, iSparseStart, img2CentCell, ERhoPrim, SSqrReal)

    !> Electronic solver information
    type(TElsiSolver), intent(inout) :: this

    !> Environment settings
    type(TEnvironment), intent(inout) :: env

    !> Dense matrix descriptor
    type(TDenseDescr), intent(in) :: denseDesc

    !> list of neighbours for each atom
    type(TNeighbourList), intent(in) :: neighbourList

    !> Number of neighbours for each of the atoms
    integer, intent(in) :: nNeighbourSK(:)

    !> Is the geometry helical
    logical, intent(in) :: tHelical

    !> Atomic orbital information
    type(TOrbitals), intent(in) :: orb

    !> species of all atoms in the system
    integer, intent(in) :: species(:)

    !> atomic coordinates
    real(dp), intent(in) :: coord(:,:)

    !> Index array for the start of atomic blocks in sparse arrays
    integer, intent(in) :: iSparseStart(:,:)

    !> map from image atoms to the original unique atom
    integer, intent(in) :: img2CentCell(:)

    !> Energy weighted sparse matrix
    real(dp), intent(out) :: ERhoPrim(:)

    !> Storage for dense overlap matrix
    real(dp), intent(inout), allocatable :: SSqrReal(:,:)

    real(dp), allocatable :: EDMnzValLocal(:)

    ERhoPrim(:) = 0.0_dp
    if (this%isSparse) then
      allocate(EDMnzValLocal(this%elsiCsc%nnzLocal))
      call elsi_get_edm_real_sparse(this%handle, EDMnzValLocal)
      if (tHelical) then
        call this%elsiCsc%convertElsiToPackedReal(neighbourList%iNeighbour, nNeighbourSK, orb,&
            & species, coord, denseDesc%iAtomStart, iSparseStart, img2CentCell, EDMnzValLocal,&
            & ErhoPrim)
      else
        call this%elsiCsc%convertElsiToPackedReal(neighbourList%iNeighbour, nNeighbourSK, orb,&
            & denseDesc%iAtomStart, iSparseStart, img2CentCell, EDMnzValLocal, ErhoPrim)
      end if

    else
      call elsi_get_edm_real(this%handle, SSqrReal)
      if (tHelical) then
        call packRhoHelicalRealBlacs(env%blacs, denseDesc, SSqrReal, neighbourList%iNeighbour,&
            & nNeighbourSK, iSparseStart, img2CentCell, orb, species, coord, ERhoPrim)
      else
        call packRhoRealBlacs(env%blacs, denseDesc, SSqrReal, neighbourList%iNeighbour,&
            & nNeighbourSK, orb%mOrb, iSparseStart, img2CentCell, ERhoPrim)
      end if
    end if

    ! add contributions from different spin channels together if necessary
    call mpifx_allreduceip(env%mpi%globalComm, ERhoPrim, MPI_SUM)

  end subroutine getEDensityMtxReal


  !> Returns the energy weighted density matrix using ELSI non-diagonalisation routines.
  subroutine getEDensityMtxCmplx(this, env, denseDesc, kPoint, kWeight, neighbourList,&
      & nNeighbourSK, tHelical, orb, species, coord, iSparseStart, img2CentCell, iCellVec, cellVec,&
      & parallelKS, ERhoPrim, SSqrCplx)

    !> Electronic solver information
    type(TElsiSolver), intent(inout) :: this

    !> Environment settings
    type(TEnvironment), intent(inout) :: env

    !> Dense matrix descriptor
    type(TDenseDescr), intent(in) :: denseDesc

    !> K-points
    real(dp), intent(in) :: kPoint(:,:)

    !> Weights for k-points
    real(dp), intent(in) :: kWeight(:)

    !> list of neighbours for each atom
    type(TNeighbourList), intent(in) :: neighbourList

    !> Number of neighbours for each of the atoms
    integer, intent(in) :: nNeighbourSK(:)

    !> Is the geometry helical
    logical, intent(in) :: tHelical

    !> Atomic orbital information
    type(TOrbitals), intent(in) :: orb

    !> species of all atoms in the system
    integer, intent(in) :: species(:)

    !> atomic coordinates
    real(dp), intent(in) :: coord(:,:)

    !> Index array for the start of atomic blocks in sparse arrays
    integer, intent(in) :: iSparseStart(:,:)

    !> map from image atoms to the original unique atom
    integer, intent(in) :: img2CentCell(:)

    !> Index for which unit cell atoms are associated with
    integer, intent(in) :: iCellVec(:)

    !> Vectors (in units of the lattice constants) to cells of the lattice
    real(dp), intent(in) :: cellVec(:,:)

    !> K-points and spins to process
    type(TParallelKS), intent(in) :: parallelKS

    !> Energy weighted sparse matrix
    real(dp), intent(out) :: ERhoPrim(:)

    !> Storage for dense overlap matrix (complex case)
    complex(dp), intent(inout), allocatable :: SSqrCplx(:,:)

    complex(dp), allocatable :: EDMnzValLocal(:)
    integer :: iK

    ERhoPrim(:) = 0.0_dp
    ! iKS always 1
    iK = parallelKS%localKS(1, 1)
    if (this%isSparse) then
      allocate(EDMnzValLocal(this%elsiCsc%nnzLocal))
      call elsi_get_edm_complex_sparse(this%handle, EDMnzValLocal)
      if (tHelical) then
        call this%elsiCsc%convertElsiToPackedCmplx(neighbourList%iNeighbour, nNeighbourSK, orb,&
            & species, coord, denseDesc%iAtomStart, iSparseStart, img2CentCell, kPoint(:,iK),&
            & kWeight(iK), iCellVec, cellVec, EDMnzValLocal, ERhoPrim)
      else
        call this%elsiCsc%convertElsiToPackedCmplx(neighbourList%iNeighbour, nNeighbourSK, orb,&
            & denseDesc%iAtomStart, iSparseStart, img2CentCell, kPoint(:,iK), kWeight(iK),&
            & iCellVec, cellVec, EDMnzValLocal, ERhoPrim)

      end if
    else
      call elsi_get_edm_complex(this%handle, SSqrCplx)
      if (tHelical) then
        call packRhoHelicalCplxBlacs(env%blacs, denseDesc, SSqrCplx, kPoint(:,iK), kWeight(iK),&
            & neighbourList%iNeighbour, nNeighbourSK, iCellVec, cellVec, iSparseStart,&
            & img2CentCell, orb, species, coord, ERhoPrim)
      else
        call packRhoCplxBlacs(env%blacs, denseDesc, SSqrCplx, kPoint(:,iK), kWeight(iK),&
            & neighbourList%iNeighbour, nNeighbourSK, orb%mOrb, iCellVec, cellVec, iSparseStart,&
            & img2CentCell, ERhoPrim)
      end if
    end if
    call mpifx_allreduceip(env%mpi%globalComm, ERhoPrim, MPI_SUM)

  end subroutine getEDensityMtxCmplx


  !> Returns the energy weighted density matrix using ELSI non-diagonalisation routines.
  subroutine getEDensityMtxPauli(this, env, denseDesc, kPoint, kWeight, neighbourList,&
      & nNeighbourSK, orb, iSparseStart, img2CentCell, iCellVec, cellVec, parallelKS, ERhoPrim,&
      & SSqrCplx)

    !> Electronic solver information
    type(TElsiSolver), intent(inout) :: this

    !> Environment settings
    type(TEnvironment), intent(inout) :: env

    !> Dense matrix descriptor
    type(TDenseDescr), intent(in) :: denseDesc

    !> K-points
    real(dp), intent(in) :: kPoint(:,:)

    !> Weights for k-points
    real(dp), intent(in) :: kWeight(:)

    !> list of neighbours for each atom
    type(TNeighbourList), intent(in) :: neighbourList

    !> Number of neighbours for each of the atoms
    integer, intent(in) :: nNeighbourSK(:)

    !> Atomic orbital information
    type(TOrbitals), intent(in) :: orb

    !> Index array for the start of atomic blocks in sparse arrays
    integer, intent(in) :: iSparseStart(:,:)

    !> map from image atoms to the original unique atom
    integer, intent(in) :: img2CentCell(:)

    !> Index for which unit cell atoms are associated with
    integer, intent(in) :: iCellVec(:)

    !> Vectors (in units of the lattice constants) to cells of the lattice
    real(dp), intent(in) :: cellVec(:,:)

    !> K-points and spins to process
    type(TParallelKS), intent(in) :: parallelKS

    !> Energy weighted sparse matrix
    real(dp), intent(out) :: ERhoPrim(:)

    !> Storage for dense overlap matrix
    complex(dp), intent(inout), allocatable :: SSqrCplx(:,:)

    integer :: iK

    if (this%isSparse) then
      call error("EDensity via sparse ELSI solver not supported for Pauli-matrices")
    else
      ! iKS always 1, as number of groups matches the number of k-points
      iK = parallelKS%localKS(1, 1)
      call elsi_get_edm_complex(this%handle, SSqrCplx)
      call packERhoPauliBlacs(env%blacs, denseDesc, SSqrCplx, kPoint(:,iK), kWeight(iK),&
          & neighbourList%iNeighbour, nNeighbourSK, orb%mOrb, iCellVec, cellVec, iSparseStart,&
          & img2CentCell, ERhoPrim)
      call mpifx_allreduceip(env%mpi%globalComm, ERhoPrim, MPI_SUM)
    end if

  end subroutine getEDensityMtxPauli

#:endif

end module dftbp_elsisolver
