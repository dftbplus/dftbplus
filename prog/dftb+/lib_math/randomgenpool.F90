!--------------------------------------------------------------------------------------------------!
!  DFTB+: general package for performing fast atomistic simulations                                !
!  Copyright (C) 2006 - 2020  DFTB+ developers group                                               !
!                                                                                                  !
!  See the LICENSE file for terms of usage and distribution.                                       !
!--------------------------------------------------------------------------------------------------!

#:include 'common.fypp'

!> Implements a random generator pool, returning random generators on request. The status of the
!> subsequently returned random generators (and hence the random numbers they will produce) is
!> uniquely determined by the seed value used to initialise the random generator pool itself.
module dftbp_randomgenpool
#:if WITH_MPI
  use dftbp_mpifx
#:endif
  use dftbp_environment
  use dftbp_accuracy, only : dp
  use dftbp_ranlux
  use dftbp_assert
  implicit none
  private

  public :: TRandomGenPool, init


  !> Random generator pool
  !>
  !> Note: To ensure random numbers being independent from the nr. of processes being used,
  !> all random generator pool methods must always be called collectively by all processes.
  !>
  type :: TRandomGenPool
    private

    !> random number generator for pool
    type(TRanlux), allocatable :: generator

    !> Random values that have been generated
    integer :: served = -1

    !> Compatibility to old behaviour
    logical :: oldCompat = .false.
  contains

    !> returns a random generator
    procedure :: getGenerator
  end type TRandomGenPool


  !> initialise the generator
  interface init
    module procedure RandomGenPool_init
  end interface init


  !> Size of random pool necessary in order to emulate the old global random generator.
  integer, parameter :: OLDCOMPAT_POOL_SIZE = 10

contains


  !> Intialises a random generator pool.
  subroutine RandomGenPool_init(this, env, seed, oldCompat)

    !> Instance.
    class(TRandomGenPool), intent(out) :: this

    !> Environment settings
    type(TEnvironment), intent(in) :: env

    !> Seed to use for initialisation of the random generator pool.
    !> If value is less than one, a random seed will be chosen (and passed back to the calling
    !> routine).
    integer, intent(inout) :: seed


    !> Whether the first random generator returned should deliver the same random number sequence
    !> as the old global random generator in DFTB+ (default: .false.)
    logical, intent(in), optional :: oldCompat


    !> system time if available
    integer :: timeValues(8)

    !> real temporary
    real(dp) :: rTmp

    if (env%tGlobalLead) then
      if (seed < 1) then
        call system_clock(seed)
      end if

      if (seed < 1) then
        call date_and_time(values=timeValues)
        if (timeValues(5) >= 0) then
          seed = 1000 * (60 * (60 * timeValues(5) + timeValues(6)) + timeValues(7)) + timeValues(8)
        end if
      end if

      if (seed < 1) then
        call random_seed()
        call random_number(rTmp)
        ! Make sure seed > 1
        seed = int(real(huge(seed) - 1, dp) * rTmp) + 1
      end if
    end if

  #:if WITH_MPI
    call mpifx_bcast(env%mpi%globalComm, seed)
  #:endif

    allocate(this%generator)
    call init(this%generator, 3, initSeed=seed)
    this%served = 0

    if (present(oldCompat)) then
      this%oldCompat = oldCompat
    end if

  end subroutine RandomGenPool_init


  !> Returns a random generator.
  subroutine getGenerator(this, env, randomGenerator)

    !> Instance.
    class(TRandomGenPool), intent(inout) :: this

    !> Environment settings.
    type(TEnvironment), intent(in) :: env

    !> Initialised random generator.
    type(TRanlux), allocatable, intent(out) :: randomGenerator

    integer :: seed
    real(dp) :: randomPool(OLDCOMPAT_POOL_SIZE)
    real(dp) :: rTmp

    @:ASSERT(this%served >= 0)

    ! First random generator returned needs special treatment to yield the same random numbers
    ! as the previous global random generator.
    if (this%served == 0 .and. this%oldCompat) then
      call getRandom(this%generator, randompool)
      seed = int(real(huge(seed) - 1, dp) * randompool(1)) + 1
    #:if WITH_MPI
      call mpifx_bcast(env%mpi%globalComm, seed)
    #:endif
      call move_alloc(this%generator, randomGenerator)
      allocate(this%generator)
      call init(this%generator, initSeed=seed)
    else
      call getRandom(this%generator, rTmp)
      seed = int(real(huge(seed) - 1, dp) * rTmp) + 1
    #:if WITH_MPI
      call mpifx_bcast(env%mpi%globalComm, seed)
    #:endif
      allocate(randomGenerator)
      call init(randomGenerator, initSeed=seed)
    end if
    this%served = this%served + 1

    @:ASSERT(this%served < huge(this%served))

  end subroutine getGenerator

end module dftbp_randomgenpool
