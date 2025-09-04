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
#:if WITH_SCALAPACK
  use dftbp_extlibs_scalapackfx, only : scalafx_numroc
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

  contains

    procedure, private :: TElpa_solveReal
    procedure, private :: TElpa_solveComplex
    generic :: solve => TElpa_solveReal, TElpa_solveComplex
    procedure :: reset => TElpa_reset
    procedure, private :: setConfig => Telpa_setConfig

  end type TElpa

contains


  !> Initialise the ELPA solver
  subroutine TElpa_init(this, env, inp, nBasisFn, timingLevel)

    !> Instance
    class(TElpa), intent(out) :: this

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

    na_rows = scalafx_numroc(nBasisFn, env%blacs%rowBlockSize, env%blacs%orbitalGrid%myrow,&
        & env%blacs%orbitalGrid%leadrow, env%blacs%orbitalGrid%nrow)
    na_cols = scalafx_numroc(nBasisFn, env%blacs%columnBlockSize, env%blacs%orbitalGrid%mycol,&
        & env%blacs%orbitalGrid%leadcol, env%blacs%orbitalGrid%ncol)

    status = elpa_init(20250131) ! corresponds to ELPA version 2025.01.001
    if (status /= ELPA_OK) then
      call error("ELPA error: elpa_init failed")
    end if

    this%handle => elpa_allocate(status)
    if (status /= ELPA_OK) then
      call error("ELPA error: elpa_allocate failed")
    end if

    call this%setConfig("na", nBasisFn)
    call this%setConfig("nev", nBasisFn)
    call this%setConfig("nblk", env%blacs%rowBlockSize)
    call this%setConfig("local_nrows", na_rows)
    call this%setConfig("local_ncols", na_cols)
    call this%setConfig("mpi_comm_parent", env%mpi%groupComm%id)
    call this%setConfig("blacs_context", env%blacs%orbitalGrid%ctxt)
    call this%setConfig("process_row", env%blacs%orbitalGrid%myrow)
    call this%setConfig("process_col", env%blacs%orbitalGrid%mycol)

    if (timingLevel < 0) then
      call this%setConfig("timings", 1)
      this%printTimings = env%mpi%globalComm%lead
    end if

    status = this%handle%setup()
    if (status /= elpa_ok) then
      call error("elpa error: elpa_setup failed")
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


  !> Reset the solver state when the geometry has changed
  subroutine TElpa_reset(this)

    !> Instance
    class(TElpa), intent(inout) :: this

    this%firstCall = .true.

  end subroutine TElpa_reset


  !> Set ELPA flags with error handling
  subroutine TElpa_setConfig(this, name, val)

    !> Instance
    class(TElpa), intent(inout) :: this

    !> Name of the flag
    character(*), intent(in) :: name

    !> Value of the flag
    integer, intent(in) :: val

    integer(kind=c_int) :: status

  #:if WITH_ELPA
    call this%handle%set(name, int(val, kind=c_int), status)
    if (status /= ELPA_OK) then
      call error("Error during ELPA initialization: setting " // name // " failed")
    end if
  #:endif

  end subroutine TElpa_setConfig


  !> Finalize the ELPA solver
  subroutine TElpa_final(this)

    !> Instance
    type(TElpa), intent(inout) :: this

    integer(kind=c_int) :: status

  #:if WITH_ELPA

    if (this%printTimings) then
      call this%handle%print_times("ELPA")
    end if

    if (associated(this%autotune)) then
      call elpa_autotune_deallocate(this%autotune, status)
      if (status /= ELPA_OK) then
        call error("ELPA error: elpa_autotune_deallocate failed")
      end if
    end if

    if (associated(this%handle)) then
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

    if (this%autotuning) then
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

    call this%handle%timer_start("ELPA")
    call this%handle%generalized_eigenvectors(HSqr, SSqr, eigenVals, eigenVecs,&
        & .not. this%firstCall, status)
    call this%handle%timer_stop("ELPA")

    this%firstCall = .false.

    if (status /= ELPA_OK) then
      call error("ELPA error: generalized_eigenvectors failed")
    end if

  #:else
    eigenVals(:) = 0._dp
    eigenVecs(:,:) = 0._dp
    call error("Internal error: TElpa_solve${NAME}$() called despite missing ELPA support")
  #:endif

  end subroutine TElpa_solve${NAME}$
#:endfor

end module dftbp_elecsolvers_elpa
