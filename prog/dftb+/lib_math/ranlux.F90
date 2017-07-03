!--------------------------------------------------------------------------------------------------!
!  DFTB+: general package for performing fast atomistic simulations                                !
!  Copyright (C) 2017  DFTB+ developers group                                                      !
!                                                                                                  !
!  See the LICENSE file for terms of usage and distribution.                                       !
!--------------------------------------------------------------------------------------------------!

#:include 'common.fypp'

!!* High quality pseudo random generator for "luxury pseudorandom numbers".
!!* @desc
!!* <p>
!!*   This is a subtract-and-borrow random generator proposed by Masaglia and
!!*   Zaman, implemented by F. James with the name RCARRY in 1991, and later
!!*   improved by M. Luescher in 1993. Fortran 77 coded by F. James 1993.
!!*   The current version is a repackaging of the integer based version made
!!*   by K.G. Hamilton and F. James.
!!* </p>
!!* <p>
!!*  The following luxury levels are available:
!!*  <table border="1">
!!*    <tr><th>Level</th><th>p</th><th>Description</th></tr>
!!*    <tr><td>0</td><td>24</td><td>Equivalent to the original RCARRY of
!!*      Marsaglia and Zaman, very long period, but fails many tests.</td></tr>
!!*   <tr><td>1</td><td>48</td><td>Considerable improvement in quality over
!!*     level 0, now passes the gap test, but still fails spectral test.
!!*   </td></tr>
!!*   <tr><td>2</td><td>97</td><td>Passes all known tests, but theoretically
!!*     still defective</td></tr>
!!*   <tr><td>3</td><td>223</td><td>DEFAULT VALUE. Any theoretically possible
!!*     correlations have very small chance of being observed.</td></tr>
!!*   <tr><td>4</td><td>389</td><td>Highest possible luxury, all 24 bits
!!*     chaotic.</td></tr>
!!*   </table>
!!* </p>
!!* <p>
!!*   The validation was made by obtaining the difference between the F90
!!*   version of the original code and the current module for 1e5 calls
!!*   each filling a vector with 1e6 random numbers. (i686-linux-ifort81,
!!*   DEBUG=0) Luxury level was 3, the initial seed 123456. Since the original
!!*   code uses single precision, while the current code uses double precision,
!!*   differences less than 1e-11 occur in the generated numbers. The integers
!!*   describing the inner state of the generators had been compared after each
!!*   call and had been found to be identical every time.
!!* </p>
!!* @see M. Luscher, Computer Physics Communications  79 (1994) 100
!!* @see F. James, Computer Physics Communications 79 (1994) 111
module ranlux
  use assert
  use accuracy, only : dp
  implicit none

  private

  !!* Internal variables for the luxury pseudorandom generator
  type ORanlux
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
  end type ORanlux


  !!* Creates a ranlux random number generator
  interface init
    module procedure Ranlux_init_default
    module procedure Ranlux_init_restart
  end interface

  !!* Fills a vector with random numbers
  interface getRandom
    module procedure Ranlux_getRandomVector
    module procedure Ranlux_getRandom2DArray
    module procedure Ranlux_getRandomNumber
  end interface

  !!* Return the state of the generator
  interface getState
    module procedure Ranlux_getState
  end interface

  public :: ORanlux
  public :: init, getRandom, getState


  !!* Maximal luxury level
  integer, parameter :: maxlev = 4

  !!* Default luxury level
  integer, parameter :: lxdflt = 3

  !!* Default seed
  integer, parameter :: jsdflt = 314159265

  !!* Nr. of random numbers to throw away to destroy coherence
  integer, parameter :: ndskip(0:maxlev) = (/ 0, 24, 73, 199, 365 /)

  !! 2**24 as integer
  integer, parameter :: itwo24 = 2**24

  !!* Auxiliary constant
  integer, parameter :: icons = 2147483563

  !!* Mask for the lowest 24 bits
  integer, parameter :: masklo = itwo24 - 1

  !!* Mask for all but the lowest 24 bits
  integer, parameter :: maskhi = not(masklo)


contains

  !!* Creates and initializes a random generator
  !!* @param self     Initialized random generator on exit
  !!* @param luxlev   Luxury level. Possible values: 0, 1, 2, 3, 4. (Default: 3)
  !!* @param initSeed Initial seed value. (Default: 314159265)
  subroutine Ranlux_init_default(self, luxlev, initSeed)
    type(ORanlux), intent(out) :: self
    integer, intent(in), optional :: luxlev
    integer, intent(in), optional :: initSeed

    integer :: jseed
    integer :: ii, kk

  #:call ASSERT_CODE
    if (present(luxlev)) then
      @:ASSERT(luxlev >= 0 .and. luxlev <= maxlev)
    end if
    if (present(initSeed)) then
      @:ASSERT(initSeed > 0)
    end if
  #:endcall ASSERT_CODE

    !! Set luxury level
    self%luxlev = lxdflt
    if (present(luxlev)) then
      if (luxlev >= 0 .and. luxlev <= maxlev) then
        self%luxlev = luxlev
      end if
    end if

    !! Set initial seed
    jseed = jsdflt
    if (present(initSeed)) then
      if (initSeed > 0) then
        jseed = initSeed
      end if
    end if

    self%nskip = ndskip(self%luxlev)
    self%in24 = 0
    self%twom24 = 1.0_dp

    !! Calculate seeds
    do ii = 1, 24
      self%twom24 = self%twom24 * 0.5_dp
      kk = jseed / 53668
      jseed = 40014 * (jseed-kk*53668) - kk * 12211
      if (jseed < 0) then
        jseed = jseed + icons
      end if
      self%iseeds(ii) = mod(jseed,itwo24)
      self%next(ii) = ii - 1
    end do

    self%twom12 = self%twom24 * 4096.0_dp
    self%next(1) = 24
    self%i24 = 24
    self%j24 = 10
    self%icarry = 0
    if (iand(self%iseeds(24), maskhi) /= 0) then
      self%icarry = 1
    end if

  end subroutine Ranlux_init_default



  !!* Creates and initializes a random generator with previously saved
  !!* values.
  !!* @param self   Initialized random generator instance on exit
  !!* @param isdext Contains the state of a saved generator as
  !!*   produced by Ranlux_getState.
  subroutine Ranlux_init_restart(self, isdext)
    type(ORanlux), intent(out) :: self
    integer, intent(in) :: isdext(:)

    integer :: ii, isd

    @:ASSERT(size(isdext) == 25)

    self%twom24 = 1.0_dp
    do ii = 1, 24
      self%next(ii) = ii - 1
      self%twom24 = self%twom24 * 0.5_dp
    end do
    self%next(1) = 24
    self%twom12 = self%twom24 * 4096.0_dp
    self%iseeds(1:24) = isdext(1:24)
    self%icarry = 0
    if (isdext(25) < 0) then
      self%icarry = 1
    end if

    isd = iabs(isdext(25))
    self%i24 = mod(isd,100)
    isd = isd / 100
    self%j24 = mod(isd,100)
    isd = isd / 100
    self%in24 = mod(isd,100)
    isd = isd / 100
    self%luxlev = isd
    if (self%luxlev <= maxlev) then
      self%nskip = ndskip(self%luxlev)
    else if (self%luxlev >= 24) then
      self%nskip = self%luxlev - 24
    else
      self%nskip = ndskip(maxlev)
      self%luxlev = maxlev
    end if

  end subroutine Ranlux_init_restart


  !!* Fills a given vector with random numbers.
  !!* @param self Ranlux instance
  !!* @param rvec Vector containing the random numbers on exit.
  subroutine Ranlux_getRandomVector(self, rvec)
    type(ORanlux), intent(inout) :: self
    real(dp), intent(out) :: rvec(:)

    call getRandomVector_local(rvec, self%iseeds, self%icarry, self%in24, &
        &self%i24, self%j24, self%next, self%nskip, self%twom24, self%twom12)

  end subroutine Ranlux_getRandomVector


  !!* Fills a given 2D array with random numbers.
  !!* @param self Ranlux instance
  !!* @param r2Darray Vector containing the random numbers on exit.
  subroutine Ranlux_getRandom2DArray(self, r2Darray)
    type(ORanlux), intent(inout) :: self
    real(dp), intent(out) :: r2Darray(:,:)

    real(dp), allocatable :: rvec(:)

    allocate(rvec(size(r2Darray,dim=1)*size(r2Darray,dim=2)))
    call getRandomVector_local(rvec, self%iseeds, self%icarry, self%in24, &
        &self%i24, self%j24, self%next, self%nskip, self%twom24, self%twom12)
    r2Darray = reshape(rvec,shape(r2Darray))

  end subroutine Ranlux_getRandom2DArray


  !!* Returns a random number
  !!* @param self Ranlux instance
  !!* @param rnum Contains the random number on exit.
  subroutine Ranlux_getRandomNumber(self, rnum)
    type(ORanlux), intent(inout) :: self
    real(dp), intent(out) :: rnum

    real(dp) :: rvec(1)

    call getRandomVector_local(rvec, self%iseeds, self%icarry, self%in24, &
        &self%i24, self%j24, self%next, self%nskip, self%twom24, self%twom12)
    rnum = rvec(1)

  end subroutine Ranlux_getRandomNumber



  !!* Workhorse for the Ranlux_getRandom* methods.
  !!* @param rvec    Vector containing the random numbers on exit
  !!* @param iseeds  Stored seeds
  !!* @param icarry  Carry bit
  !!* @param in24    Auxiliary variable
  !!* @param i24     Auxiliary variable
  !!* @param j24     Auxiliary variable
  !!* @param next    Auxiliary variable
  !!* @param nskip   Nr. of numbers to throw away to destroy coherence
  !!* @param twom24  2**-24 as real
  !!* @param twom12  2**-12 as real
  subroutine getRandomVector_local(rvec, iseeds, icarry, in24, i24, j24, next, &
      & nskip, twom24, twom12)
    real(dp), intent(out) :: rvec(:)
    integer, intent(inout) :: iseeds(:)
    integer, intent(inout) :: icarry
    integer, intent(inout) :: in24
    integer, intent(inout) :: i24, j24
    integer, intent(in) :: next(:)
    real(dp), intent(in) :: twom24, twom12
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



  !!* Saves the state of the random generator in an integer array
  !!* @param self   Ranlux instance.
  !!* @param isdext Contains the state of the generator as integer array.
  subroutine Ranlux_getState(self, isdext)
    type(ORanlux), intent(in) :: self
    integer, intent(out) :: isdext(:)

    @:ASSERT(size(isdext) == 25)

    isdext(1:24) = self%iseeds(1:24)
    isdext(25) = self%i24 + 100 * self%j24 + 10000 * self%in24 &
        &+ 1000000 * self%luxlev
    if (self%icarry /= 0) then
      isdext(25) = -isdext(25)
    end if

  end subroutine Ranlux_getState


end module ranlux
