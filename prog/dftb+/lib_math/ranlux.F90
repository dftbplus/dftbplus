!--------------------------------------------------------------------------------------------------!
!  DFTB+: general package for performing fast atomistic simulations                                !
!  Copyright (C) 2006 - 2020  DFTB+ developers group                                               !
!                                                                                                  !
!  See the LICENSE file for terms of usage and distribution.                                       !
!--------------------------------------------------------------------------------------------------!

#:include 'common.fypp'

!> High quality pseudo random generator for "luxury pseudorandom numbers".
!>
!>   This is a subtract-and-borrow random generator proposed by Masaglia and
!>   Zaman, implemented by F. James with the name RCARRY in 1991, and later
!>   improved by M. Luescher in 1993. Fortran 77 coded by F. James 1993.
!>   The current version is a repackaging of the integer based version made
!>   by K.G. Hamilton and F. James.
!>
!>  The following luxury levels are available:
!>
!> Level 0 p=24 Equivalent to the original RCARRY of Marsaglia and Zaman, very long period, but
!> fails many tests.
!>
!> Level 1 p=48 Considerable improvement in quality over level 0, now passes the gap test, but still
!> fails spectral test.
!>
!> Level 2 p=97 Passes all known tests, but theoretically still defective
!>
!> Level 3 p=223 DEFAULT VALUE. Any theoretically possible correlations have very small chance of
!> being observed.
!>
!> Level 4 p=389 Highest possible luxury, all 24 bits chaotic.
!>
!>
!> The validation was made by obtaining the difference between the F90 version of the original code
!> and the current module for 1e5 calls each filling a vector with 1e6 random
!> numbers. (i686-linux-ifort81, DEBUG=0) Luxury level was 3, the initial seed 123456. Since the
!> original code uses single precision, while the current code uses double precision, differences
!> less than 1e-11 occur in the generated numbers. The integers describing the inner state of the
!> generators had been compared after each call and had been found to be identical every time.
!>
!> See M. Luscher, Computer Physics Communications 79 (1994) 100 and F. James, Computer Physics
!> Communications 79 (1994) 111
module dftbp_ranlux
  use dftbp_assert
  use dftbp_accuracy, only : dp
  implicit none

  private


  !> Internal variables for the luxury pseudorandom generator
  type TRanlux
    integer :: next(24)
    integer :: luxlev
    integer :: nskip
    integer :: in24
    integer :: i24
    integer :: j24
    integer :: iseeds(24)
    integer :: icarry
    real(dp) :: twom24
    real(dp) :: twom12
  end type TRanlux


  !> Creates a ranlux random number generator
  interface init
    module procedure Ranlux_init_default
    !module procedure Ranlux_init_restart
  end interface init


  !> Fills a vector with random numbers
  interface getRandom
    module procedure Ranlux_getRandomVector
    module procedure Ranlux_getRandom2DArray
    module procedure Ranlux_getRandomNumber
  end interface getRandom


  !> Return the state of the generator
  interface getState
    module procedure Ranlux_getState
  end interface getState

  public :: TRanlux
  public :: init, getRandom, getState


  !> Maximal luxury level
  integer, parameter :: maxlev = 4


  !> Default luxury level
  integer, parameter :: lxdflt = 3


  !> Default seed
  integer, parameter :: jsdflt = 314159265


  !> Nr. of random numbers to throw away to destroy coherence
  integer, parameter :: ndskip(0:maxlev) = (/ 0, 24, 73, 199, 365 /)


  !> 2**24 as integer
  integer, parameter :: itwo24 = 2**24


  !> Auxiliary constant
  integer, parameter :: icons = 2147483563


  !> Mask for the lowest 24 bits
  integer, parameter :: masklo = itwo24 - 1


  !> Mask for all but the lowest 24 bits
  integer, parameter :: maskhi = not(masklo)

contains


  !> Creates and initializes a random generator
  subroutine Ranlux_init_default(this, luxlev, initSeed)

    !> Initialized random generator on exit
    type(TRanlux), intent(out) :: this

    !> Luxury level. Possible values: 0, 1, 2, 3, 4. (Default: 3)
    integer, intent(in), optional :: luxlev

    !> Initial seed value. (Default: 314159265)
    integer, intent(in), optional :: initSeed

    integer :: jseed
    integer :: ii, kk

  #:block DEBUG_CODE
    if (present(luxlev)) then
      @:ASSERT(luxlev >= 0 .and. luxlev <= maxlev)
    end if
    if (present(initSeed)) then
      @:ASSERT(initSeed > 0)
    end if
  #:endblock DEBUG_CODE

    !! Set luxury level
    this%luxlev = lxdflt
    if (present(luxlev)) then
      if (luxlev >= 0 .and. luxlev <= maxlev) then
        this%luxlev = luxlev
      end if
    end if

    !! Set initial seed
    jseed = jsdflt
    if (present(initSeed)) then
      if (initSeed > 0) then
        jseed = initSeed
      end if
    end if

    this%nskip = ndskip(this%luxlev)
    this%in24 = 0
    this%twom24 = 1.0_dp

    !! Calculate seeds
    do ii = 1, 24
      this%twom24 = this%twom24 * 0.5_dp
      kk = jseed / 53668
      jseed = 40014 * (jseed-kk*53668) - kk * 12211
      if (jseed < 0) then
        jseed = jseed + icons
      end if
      this%iseeds(ii) = mod(jseed,itwo24)
      this%next(ii) = ii - 1
    end do

    this%twom12 = this%twom24 * 4096.0_dp
    this%next(1) = 24
    this%i24 = 24
    this%j24 = 10
    this%icarry = 0
    if (iand(this%iseeds(24), maskhi) /= 0) then
      this%icarry = 1
    end if

  end subroutine Ranlux_init_default


  !> Creates and initializes a random generator with previously saved values.
  subroutine Ranlux_init_restart(this, isdext)

    !> Initialized random generator instance on exit
    type(TRanlux), intent(out) :: this

    !> Contains the state of a saved generator as produced by Ranlux_getState.
    integer, intent(in) :: isdext(:)

    integer :: ii, isd

    @:ASSERT(size(isdext) == 25)

    this%twom24 = 1.0_dp
    do ii = 1, 24
      this%next(ii) = ii - 1
      this%twom24 = this%twom24 * 0.5_dp
    end do
    this%next(1) = 24
    this%twom12 = this%twom24 * 4096.0_dp
    this%iseeds(1:24) = isdext(1:24)
    this%icarry = 0
    if (isdext(25) < 0) then
      this%icarry = 1
    end if

    isd = iabs(isdext(25))
    this%i24 = mod(isd,100)
    isd = isd / 100
    this%j24 = mod(isd,100)
    isd = isd / 100
    this%in24 = mod(isd,100)
    isd = isd / 100
    this%luxlev = isd
    if (this%luxlev <= maxlev) then
      this%nskip = ndskip(this%luxlev)
    else if (this%luxlev >= 24) then
      this%nskip = this%luxlev - 24
    else
      this%nskip = ndskip(maxlev)
      this%luxlev = maxlev
    end if

  end subroutine Ranlux_init_restart


  !> Fills a given vector with random numbers.
  subroutine Ranlux_getRandomVector(this, rvec)

    !> Ranlux instance
    type(TRanlux), intent(inout) :: this

    !> Vector containing the random numbers on exit.
    real(dp), intent(out) :: rvec(:)

    call getRandomVector_local(rvec, this%iseeds, this%icarry, this%in24, &
        &this%i24, this%j24, this%next, this%nskip, this%twom24, this%twom12)

  end subroutine Ranlux_getRandomVector


  !> Fills a given 2D array with random numbers.
  subroutine Ranlux_getRandom2DArray(this, r2Darray)

    !> Ranlux instance
    type(TRanlux), intent(inout) :: this

    !> Vector containing the random numbers on exit.
    real(dp), intent(out) :: r2Darray(:,:)

    real(dp), allocatable :: rvec(:)

    allocate(rvec(size(r2Darray,dim=1)*size(r2Darray,dim=2)))
    call getRandomVector_local(rvec, this%iseeds, this%icarry, this%in24, &
        &this%i24, this%j24, this%next, this%nskip, this%twom24, this%twom12)
    r2Darray = reshape(rvec,shape(r2Darray))

  end subroutine Ranlux_getRandom2DArray


  !> Returns a random number
  subroutine Ranlux_getRandomNumber(this, rnum)

    !> Ranlux instance
    type(TRanlux), intent(inout) :: this

    !> Contains the random number on exit.
    real(dp), intent(out) :: rnum

    real(dp) :: rvec(1)

    call getRandomVector_local(rvec, this%iseeds, this%icarry, this%in24, &
        &this%i24, this%j24, this%next, this%nskip, this%twom24, this%twom12)
    rnum = rvec(1)

  end subroutine Ranlux_getRandomNumber


  !> Workhorse for the Ranlux_getRandom* methods.
  subroutine getRandomVector_local(rvec,iseeds,icarry,in24,i24,j24,next,nskip,twom24,twom12)

    !> Vector containing the random numbers on exit
    real(dp), intent(out) :: rvec(:)

    !> Stored seeds
    integer, intent(inout) :: iseeds(:)

    !> Carry bit
    integer, intent(inout) :: icarry

    !> Auxiliary variable
    integer, intent(inout) :: in24

    !> Auxiliary variable
    integer, intent(inout) :: i24

    !> Auxiliary variable
    integer, intent(inout) :: j24

    !> Auxiliary variable
    integer, intent(in) :: next(:)

    !> Nr. of numbers to throw away to destroy coherence
    real(dp), intent(in) :: twom24

    !> 2**-24 as real
    real(dp), intent(in) :: twom12

    !> 2**-12 as real
    integer, intent(in) :: nskip

    integer :: lenv, ivec, iuni
    real(dp) :: uni
    integer :: isk

    @:ASSERT(size(iseeds) == 24)
    @:ASSERT(size(next) == 24)

    lenv = size(rvec)
    do ivec = 1, lenv
      iuni = iseeds(j24) - iseeds(i24) - icarry
      if (iand(iuni, maskhi) /= 0) then
        iuni = iand(iuni, masklo)
        icarry = 1
      else
        icarry = 0
      end if
      iseeds(i24) = iuni
      i24 = next(i24)
      j24 = next(j24)
      uni = real(iuni, dp) * twom24
      if (uni < twom12) then
        uni = uni + (real(iseeds(j24), dp)*twom24)*twom24
        if (uni <= 0.0_dp) then
          uni = twom24 * twom24
        end if
      end if
      rvec(ivec) = real(uni, dp)

      in24 = in24 + 1
      if (in24 == 24) then
        in24 = 0
        do isk = 1, nskip
          iuni = iseeds(j24) - iseeds(i24) - icarry
          if (iand(iuni, maskhi) /= 0) then
            iuni = iand(iuni, masklo)
            icarry = 1
          else
            icarry = 0
          end if
          iseeds(i24) = iuni
          i24 = next(i24)
          j24 = next(j24)
        end do
      end if
    end do

  end subroutine getRandomVector_local


  !> Saves the state of the random generator in an integer array
  subroutine Ranlux_getState(this, isdext)

    !> Ranlux instance.
    type(TRanlux), intent(in) :: this

    !> Contains the state of the generator as integer array.
    integer, intent(out) :: isdext(:)

    @:ASSERT(size(isdext) == 25)

    isdext(1:24) = this%iseeds(1:24)
    isdext(25) = this%i24 + 100 * this%j24 + 10000 * this%in24 &
        &+ 1000000 * this%luxlev
    if (this%icarry /= 0) then
      isdext(25) = -isdext(25)
    end if

  end subroutine Ranlux_getState

end module dftbp_ranlux
