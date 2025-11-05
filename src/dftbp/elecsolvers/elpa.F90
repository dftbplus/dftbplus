!--------------------------------------------------------------------------------------------------!
!  DFTB+: general package for performing fast atomistic simulations                                !
!  Copyright (C) 2006 - 2025  DFTB+ developers group                                               !
!                                                                                                  !
!  See the LICENSE file for terms of usage and distribution.                                       !
!--------------------------------------------------------------------------------------------------!

#:include 'common.fypp'

!> Contains the interface to the ELPA solver
module dftbp_elecsolvers_elpa
  use, intrinsic :: iso_c_binding, only : c_int
  use dftbp_common_accuracy, only : dp
  use dftbp_common_environment, only : TEnvironment
#:if WITH_ELPA
  use dftbp_extlibs_elpa, only : elpa_allocate, elpa_autotune_deallocate, elpa_autotune_t,&
      & elpa_deallocate, elpa_init, elpa_t, elpa_uninit, ELPA_AUTOTUNE_DOMAIN_COMPLEX,&
      & ELPA_AUTOTUNE_DOMAIN_REAL, ELPA_AUTOTUNE_FAST, ELPA_AUTOTUNE_MEDIUM, ELPA_OK,&
      & ELPA_SOLVER_1STAGE, ELPA_SOLVER_2STAGE
#:else
  use dftbp_extlibs_elpa, only : elpa_autotune_t, elpa_t
#:endif
#:if WITH_MPI
  use dftbp_extlibs_mpifx, only : mpifx_bcast, mpifx_comm
#:endif
#:if WITH_SCALAPACK
  use dftbp_extlibs_scalapackfx, only : blacsfx_gemr2d, blacsgrid, scalafx_getdescriptor,&
      & scalafx_numroc, CTXT_, DLEN_
#:endif
  use dftbp_io_message, only : error
  implicit none

  private
  public :: TElpaInp
  public :: TElpa, TElpa_init, TElpa_final


  !> Input data for the ELPA solver
  type :: TElpaInp

    !> Choice of ELPA solver
    integer :: solver = 2

    !> Enable ELPA autotuning
    logical :: autotune = .false.

    !> Enable GPU usage in ELPA
    logical :: gpu = .false.

    !> On what fraction of the original number of ranks to redistribute the matrix
    integer :: redistributeFactor = 1

  end type TElpaInp


  !> ELPA solver state
  type :: TElpa
    private

    !> Handle of the ELPA instance
    class(elpa_t), pointer :: handle

    !> Handle of the ELPA autotuning
    class(elpa_autotune_t), pointer :: autotune

    !> Whether GPUs are enabled
    logical :: gpu = .false.

    !> Whether the solver is called for the first time for the current geometry
    logical :: firstCall = .true.

    !> Whether to print ELPA's timing information an the end
    logical :: printTimings = .false.

    !> Whether we are currently autotuning
    logical :: autotuning = .false.

    !> Whether we should redistribute the matrix each call
    logical :: redistributing = .false.

    !> Whether the current process holds parts of the redistributed matrix and joins ELPA calls
    logical :: joinElpaCalls = .true.

    !> Global size of the square matrix
    integer :: matrixSize

    !> Number of rows in the local matrix
    integer :: matrixLocalRows = 1

    !> Number of columns in the local matrix
    integer :: matrixLocalColumns = 1

  #:if WITH_SCALAPACK
    !> BLACS grid to use for redistribution
    type(blacsgrid) :: redistributeGrid
  #:endif

  #:if WITH_MPI
    !> MPI communicator to use for redistribution
    type(mpifx_comm) :: redistributeComm

    !> MPI communicator of all ranks in the current group
    type(mpifx_comm) :: groupComm
  #:endif

    !> BLACS context
    integer :: contextOrig

  #:if WITH_SCALAPACK
    !> Original  descriptor of the matrix
    integer :: descOrig(DLEN_)

    !> Descriptor to be used in ELPA, possibly redistributed
    integer :: desc(DLEN_)
  #:endif

    !> First matrix used for redistribution
    real(dp), allocatable :: matrixReal1(:,:)

    !> Second matrix used for redistribution
    real(dp), allocatable :: matrixReal2(:,:)

    !> Eigenvector storage used for redistribution
    real(dp), allocatable :: eigenvectorsReal(:,:)

    !> First matrix used for redistribution
    complex(dp), allocatable :: matrixComplex1(:,:)

    !> Second matrix used for redistribution
    complex(dp), allocatable :: matrixComplex2(:,:)

    !> Eigenvector storage used for redistribution
    complex(dp), allocatable :: eigenvectorsComplex(:,:)

  contains

    procedure, private :: TElpa_solveReal
    procedure, private :: TElpa_solveComplex
    generic :: solve => TElpa_solveReal, TElpa_solveComplex
    procedure :: reset => TElpa_reset
  #:if WITH_ELPA
    procedure, private :: setConfig => Telpa_setConfig
    procedure, private :: initConfig => Telpa_initConfig
    procedure, private :: initRedistribute => Telpa_initRedistribute
  #:endif

  end type TElpa

contains


  !> Initialise the ELPA solver
  subroutine TElpa_init(this, env, inp, nBasisFn, timingLevel)

    !> Instance
    type(TElpa), intent(out) :: this

    !> Environment settings
    type(TEnvironment), intent(in) :: env

    !> Input structure read from config file
    type(TElpaInp), intent(in) :: inp

    !> Number of orbitals in the system
    integer, intent(in) :: nBasisFn

    !> Timing level read from config file
    integer, intent(in) :: timingLevel

    integer :: na_rows, na_cols
    integer(kind=c_int) :: status

  #:if WITH_ELPA

    this%matrixSize = nBasisFn
    this%contextOrig = env%blacs%orbitalGrid%ctxt

    call scalafx_getdescriptor(env%blacs%orbitalGrid, this%matrixSize, this%matrixSize,&
        & env%blacs%rowBlockSize, env%blacs%columnBlockSize, this%descOrig)

    if (env%blacs%rowBlockSize /= env%blacs%columnBlockSize) then
      call error("Error during ELPA initialization: different block sizes for rows and columns")
    end if

    status = elpa_init(20250131) ! corresponds to ELPA version 2025.01.001
    if (status /= ELPA_OK) then
      call error("ELPA error: elpa_init failed")
    end if

    if (inp%redistributeFactor /= 1) then
      call this%initRedistribute(inp%redistributeFactor, env%mpi%groupComm)

      if (.not. this%joinElpaCalls) then
        this%desc(:) = 0
        this%desc(CTXT_) = -1
        return
      end if
    end if

    this%handle => elpa_allocate(status)
    if (status /= ELPA_OK) then
      call error("ELPA error: elpa_allocate failed")
    end if

    if (this%redistributing) then
      call this%initConfig(this%redistributeGrid, this%redistributeComm, env%blacs%rowBlockSize)
    else
      call this%initConfig(env%blacs%orbitalGrid, env%mpi%groupComm, env%blacs%rowBlockSize)
    end if

    if (timingLevel < 0) then
      call this%setConfig("timings", 1)
      if (this%redistributing) then
        this%printTimings = this%redistributeComm%lead
      else
        this%printTimings = env%mpi%groupComm%lead
      end if
    end if

    status = this%handle%setup()
    if (status /= elpa_ok) then
      call error("ELPA error: elpa_setup failed")
    end if

    if (inp%gpu) then
      call this%setConfig("nvidia-gpu", 1)

      status = this%handle%setup_gpu()
      if (status /= ELPA_OK) then
        call error("ELPA error: elpa_setup_gpu failed")
      end if
      this%gpu = .true.
    end if

    if (inp%autotune) then
      this%autotuning = .true.
    else
      select case (inp%solver)
        case (1)
          call this%setConfig("solver", ELPA_SOLVER_1STAGE)
        case (2)
          call this%setConfig("solver", ELPA_SOLVER_2STAGE)
        case default
          call error("Invalid solver choice for ELPA")
      end select
    end if

  #:else
    call error("Internal error: TElpa_init() called despite missing ELPA support")
  #:endif

  end subroutine TElpa_init


#:if WITH_ELPA
  !> Initialize ELPA settings
  subroutine TElpa_initConfig(this, grid, comm, nblk)

    !> Instance
    class(TElpa), intent(inout) :: this

    !> BLACS grid to use in ELPA
    type(blacsgrid), intent(in) :: grid

    !> MPI communicator to use in ELPA
    type(mpifx_comm), intent(in) :: comm

    !> BLACS block size
    integer, intent(in) :: nblk

    call scalafx_getdescriptor(grid, this%matrixSize, this%matrixSize, nblk, nblk, this%desc)

    this%matrixLocalRows = scalafx_numroc(this%matrixSize, nblk, grid%myrow, grid%leadrow,&
        & grid%nrow)
    this%matrixLocalColumns = scalafx_numroc(this%matrixSize, nblk, grid%mycol, grid%leadcol,&
        & grid%ncol)

    call this%setConfig("na", this%matrixSize)
    call this%setConfig("nev", this%matrixSize)
    call this%setConfig("nblk", nblk)
    call this%setConfig("local_nrows", this%matrixLocalRows)
    call this%setConfig("local_ncols", this%matrixLocalColumns)
    call this%setConfig("mpi_comm_parent", comm%id)
    call this%setConfig("blacs_context", grid%ctxt)
    call this%setConfig("process_row", grid%myrow)
    call this%setConfig("process_col", grid%mycol)

  end subroutine TElpa_initConfig


  !> Initialize matrix redistribution
  subroutine TElpa_initRedistribute(this, nprocStep, groupComm)

    !> Instance
    class(TElpa), intent(inout) :: this

    !> How many processes to choose from the group communicator
    !> nprocStep == 2 means: select every second rank
    integer, intent(in) :: nprocStep

    !> Group communicator from which a new communicator will be created
    type(mpifx_comm), intent(in) :: groupComm

    integer :: np_rows, np_cols, row, col
    integer :: nprocs, proc, splitkey, rankkey
    integer :: status
    integer, allocatable :: gridmap(:,:)

    if (nprocStep < 1 .or. modulo(groupComm%size, nprocStep) /= 0) then
      call error("Invalid value for number of processes in ELPA redistribution")
    end if
    nprocs = groupComm%size / nprocStep

    do np_rows = nint(sqrt(real(nprocs))), 2, -1
      if (mod(nprocs, np_rows) == 0) exit
    enddo
    np_cols = nprocs / np_rows

    allocate(gridmap(np_rows, np_cols))

    splitkey = 0
    rankkey = 0
    proc = 0
    do row = 1, np_rows
      do col = 1, np_cols
        gridmap(row, col) = proc
        if (proc == groupComm%rank) then
          splitkey = 1
        end if
        proc = proc + nprocStep
      end do
    end do

    call this%redistributeGrid%initmappedgrids(gridmap)

    if (splitkey == 1) then
      rankkey = this%redistributeGrid%iproc
    end if

    call groupComm%split(splitkey, rankkey, this%redistributeComm, status)
    if (status /= 0) then
      call error("Error during communicator setup for ELPA redistribution")
    end if

    this%joinElpaCalls = splitkey == 1
    this%redistributing = .true.
    this%groupComm = groupComm

  end subroutine TElpa_initRedistribute
#:endif


  !> Reset the solver state when the geometry has changed
  subroutine TElpa_reset(this)

    !> Instance
    class(TElpa), intent(inout) :: this

    this%firstCall = .true.

  end subroutine TElpa_reset


#:if WITH_ELPA
  !> Set ELPA flags with error handling
  subroutine TElpa_setConfig(this, name, val)

    !> Instance
    class(TElpa), intent(inout) :: this

    !> Name of the flag
    character(*), intent(in) :: name

    !> Value of the flag
    integer, intent(in) :: val

    integer(kind=c_int) :: status

    call this%handle%set(name, int(val, kind=c_int), status)
    if (status /= ELPA_OK) then
      call error("Error during ELPA initialization: setting " // name // " failed")
    end if

  end subroutine TElpa_setConfig
#:endif


  !> Finalize the ELPA solver
  subroutine TElpa_final(this)

    !> Instance
    type(TElpa), intent(inout) :: this

    integer(kind=c_int) :: status

  #:if WITH_ELPA

    if (this%redistributing) then
      call this%redistributeGrid%destruct()
      call this%redistributeComm%free()

      if (allocated(this%matrixReal1)) then
        deallocate(this%matrixReal1)
      end if
      if (allocated(this%matrixReal2)) then
        deallocate(this%matrixReal2)
      end if
      if (allocated(this%eigenvectorsReal)) then
        deallocate(this%eigenvectorsReal)
      end if
      if (allocated(this%matrixComplex1)) then
        deallocate(this%matrixComplex1)
      end if
      if (allocated(this%matrixComplex2)) then
        deallocate(this%matrixComplex2)
      end if
      if (allocated(this%eigenvectorsComplex)) then
        deallocate(this%eigenvectorsComplex)
      end if
    end if

    if (associated(this%autotune)) then
      call elpa_autotune_deallocate(this%autotune, status)
      if (status /= ELPA_OK) then
        call error("ELPA error: elpa_autotune_deallocate failed")
      end if
    end if

    if (associated(this%handle)) then
      if (this%printTimings) then
        call this%handle%print_times("ELPA generalized_eigenvectors")
        if (this%redistributing) then
          call this%handle%print_times("ELPA redistribute")
        end if
      end if

      call elpa_deallocate(this%handle, status)
      if (status /= ELPA_OK) then
        call error("ELPA error: elpa_deallocate failed")
      end if
    end if

    call elpa_uninit(status)
    if (status /= ELPA_OK) then
      call error("ELPA error: elpa_uninit failed")
    end if

  #:else
    call error("Internal error: TElpa_final() called despite missing ELPA support")
  #:endif

  end subroutine TElpa_final


#:for DTYPE, NAME in [('complex', 'Complex'), ('real', 'Real')]
  !> Solve the eigenproblem (${DTYPE}$ case)
  subroutine TElpa_solve${NAME}$(this, HSqr, SSqr, eigenVals, eigenVecs)

    !> Instance
    class(TElpa), intent(inout) :: this

    !> Hamiltonian
    ${DTYPE}$(dp), intent(inout) :: HSqr(:,:)

    !> Overlap matrix
    ${DTYPE}$(dp), intent(inout) :: SSqr(:,:)

    !> Eigenvalues
    real(dp), intent(out) :: eigenVals(:)

    !> Eigenvectors
    ${DTYPE}$(dp), intent(out) :: eigenVecs(:,:)

    integer(kind=c_int) :: status
    logical :: unfinished

  #:if WITH_ELPA

    if (this%autotuning .and. this%joinElpaCalls) then
      if (.not. associated(this%autotune)) then
        if (this%gpu) then
          this%autotune => this%handle%autotune_setup(ELPA_AUTOTUNE_MEDIUM,&
              & ELPA_AUTOTUNE_DOMAIN_${DTYPE.upper()}$, status)
        else
          this%autotune => this%handle%autotune_setup(ELPA_AUTOTUNE_FAST,&
              & ELPA_AUTOTUNE_DOMAIN_${DTYPE.upper()}$, status)
        end if
        if (status /= elpa_ok) then
          call error("elpa error: elpa_autotune_setup failed")
        end if
      end if

      unfinished = this%handle%autotune_step(this%autotune, status)
      if (status /= ELPA_OK) then
        call error("ELPA error: elpa_autotune_step failed")
      end if

      if (.not. unfinished) then
        call this%handle%autotune_print_state(this%autotune, status)
        if (status /= ELPA_OK) then
          call error("ELPA error: elpa_autotune_print_state failed")
        end if

        call this%handle%autotune_set_best(this%autotune, status)
        if (status /= ELPA_OK) then
          call error("ELPA error: elpa_autotune_set_best failed")
        end if

        this%autotuning = .false.
      end if
    end if

    if (this%redistributing) then
      if (.not. allocated(this%matrix${NAME}$1)) then
        allocate(this%matrix${NAME}$1(this%matrixLocalRows, this%matrixLocalColumns))
      end if
      if (.not. allocated(this%matrix${NAME}$2)) then
        allocate(this%matrix${NAME}$2(this%matrixLocalRows, this%matrixLocalColumns))
      end if
      if (.not. allocated(this%eigenvectors${NAME}$)) then
        allocate(this%eigenvectors${NAME}$(this%matrixLocalRows, this%matrixLocalColumns))
      end if

      if (this%joinElpaCalls) then
        call this%handle%timer_start("ELPA redistribute")
      end if
      call blacsfx_gemr2d(this%matrixSize, this%matrixSize, HSqr, 1, 1, this%descOrig,&
          & this%matrix${NAME}$1, 1, 1, this%desc, this%contextOrig)
      call blacsfx_gemr2d(this%matrixSize, this%matrixSize, SSqr, 1, 1, this%descOrig,&
          & this%matrix${NAME}$2, 1, 1, this%desc, this%contextOrig)
      if (this%joinElpaCalls) then
        call this%handle%timer_stop("ELPA redistribute")
      end if

      if (this%joinElpaCalls) then
        call this%handle%timer_start("ELPA generalized_eigenvectors")
        call this%handle%generalized_eigenvectors(this%matrix${NAME}$1, this%matrix${NAME}$2,&
            & eigenVals, this%eigenvectors${NAME}$, .not. this%firstCall, status)
        call this%handle%timer_stop("ELPA generalized_eigenvectors")

        if (status /= ELPA_OK) then
          call error("ELPA error: generalized_eigenvectors failed")
        end if
      end if

      if (this%joinElpaCalls) then
        call this%handle%timer_start("ELPA redistribute")
      end if
      if (this%firstCall) then
        call blacsfx_gemr2d(this%matrixSize, this%matrixSize, this%matrix${NAME}$2, 1, 1,&
            & this%desc, SSqr, 1, 1, this%descOrig, this%contextOrig)
      end if
      call blacsfx_gemr2d(this%matrixSize, this%matrixSize, this%eigenvectors${NAME}$, 1, 1,&
          & this%desc, eigenVecs, 1, 1, this%descOrig, this%contextOrig)
      call mpifx_bcast(this%groupComm, eigenVals)
      if (this%joinElpaCalls) then
        call this%handle%timer_stop("ELPA redistribute")
      end if
    else
      call this%handle%timer_start("ELPA generalized_eigenvectors")
      call this%handle%generalized_eigenvectors(HSqr, SSqr, eigenVals, eigenVecs,&
          & .not. this%firstCall, status)
      call this%handle%timer_stop("ELPA generalized_eigenvectors")

      if (status /= ELPA_OK) then
        call error("ELPA error: generalized_eigenvectors failed")
      end if
    end if

    this%firstCall = .false.

  #:else
    eigenVals(:) = 0._dp
    eigenVecs(:,:) = 0._dp
    call error("Internal error: TElpa_solve${NAME}$() called despite missing ELPA support")
  #:endif

  end subroutine TElpa_solve${NAME}$
#:endfor


end module dftbp_elecsolvers_elpa
