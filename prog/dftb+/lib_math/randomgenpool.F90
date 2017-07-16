!--------------------------------------------------------------------------------------------------!
!  DFTB+: general package for performing fast atomistic simulations                                !
!  Copyright (C) 2017  DFTB+ developers group                                                      !
!                                                                                                  !
!  See the LICENSE file for terms of usage and distribution.                                       !
!--------------------------------------------------------------------------------------------------!

#:include 'common.fypp'

!> Implements a random generator pool.
!!
!! A random generator pool returns random generators on request. The status of the subsequently
!! returned random generators (and hence the random numbers they will produce) is uniquely
!! determined by the seed value used to initialise the random generator pool itself.
!!
module randomgenpool
  use accuracy, only : dp
  use ranlux
  implicit none
  private

  public :: ORandomGenPool, init

  !> Random generator pool
  type :: ORandomGenPool
    private
    type(ORanlux), allocatable :: generator
    integer :: served = -1
    logical :: oldCompat = .false.
  contains
    procedure :: getGenerator
  end type ORandomGenPool


  interface init
    module procedure RandomGenPool_init
  end interface init

  ! Size of random pool necessary in order to emulate the old global random generator.
  integer, parameter :: OLDCOMPAT_POOL_SIZE = 10

contains


  !> Intialises a random generator pool.
  subroutine RandomGenPool_init(this, seed, oldCompat)

    !> Instance.
    class(ORandomGenPool), intent(out) :: this

    !> Seed to use for initialisation of the random generator pool.
    !!
    !! If value is less than one, a random seed will be chosen (and passed back to the calling
    !! routine).
    integer, intent(inout) :: seed

    !> Whether the first random generator returned should deliver the same random number sequence
    !! as the old global random generator in DFTB+ (default: .false.)
    logical, intent(in), optional :: oldCompat

    integer :: timeValues(8)
    real(dp) :: rTmp

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

    allocate(this%generator)
    call init(this%generator, 3, initSeed=seed)
    this%served = 0

    if (present(oldCompat)) then
      this%oldCompat = oldCompat
    end if

  end subroutine RandomGenPool_init


  !> Returns a random generator.
  !!
  subroutine getGenerator(this, randomGenerator)

    !> Instance.
    class(ORandomGenPool), intent(inout) :: this

    !> Initialised random generator.
    type(ORanlux), allocatable, intent(out) :: randomGenerator

    integer :: seed
    real(dp) :: randomPool(OLDCOMPAT_POOL_SIZE)
    real(dp) :: rTmp

    @:ASSERT(this%served >= 0)

    ! First random generator returned needs special treatment to yield the same random numbers
    ! as the previous global random generator.
    if (this%served == 0 .and. this%oldCompat) then
      call getRandom(this%generator, randompool)
      seed = int(real(huge(seed) - 1, dp) * randompool(1)) + 1
      call move_alloc(this%generator, randomGenerator)
      allocate(this%generator)
      call init(this%generator, seed)
    else
      call getRandom(this%generator, rTmp)
      seed = int(real(huge(seed) - 1, dp) * rTmp) + 1
      allocate(randomGenerator)
      call init(randomGenerator, initSeed=seed)
    end if
    this%served = this%served + 1

  end subroutine getGenerator


end module randomgenpool
