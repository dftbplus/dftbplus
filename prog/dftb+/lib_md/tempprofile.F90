!--------------------------------------------------------------------------------------------------!
!  DFTB+: general package for performing fast atomistic simulations                                !
!  Copyright (C) 2017  DFTB+ developers group                                                      !
!                                                                                                  !
!  See the LICENSE file for terms of usage and distribution.                                       !
!--------------------------------------------------------------------------------------------------!

#:include 'common.fypp'

!!* Contains simple temperature profiles for molecular dynamics.
module tempprofile
  use assert
  use accuracy
  implicit none
  private

  public :: OTempProfile, init, next, getTemperature
  public :: constProf, linProf, expProf

  !!* Data for the temperature profile.
  type OTempProfile
    integer, allocatable :: tempMethods(:)
    integer, allocatable :: tempInts(:)
    real(dp), allocatable :: tempValues(:)
    real(dp) :: curTemp
    integer :: iInt
    integer :: nInt
    integer :: iStep
    real(dp) :: incr
  end type OTempProfile

  interface init
    module procedure TempProfile_init
  end interface

  interface next
    module procedure TempProfile_next
  end interface

  interface getTemperature
    module procedure TempProfile_getTemperature
  end interface

  !! Constants for the different profile options
  integer, parameter :: constProf = 1
  integer, parameter :: linProf = 2
  integer, parameter :: expProf = 3

  !! Default starting temperature
  real(dp), parameter :: startingTemp_ = minTemp


contains


  !!* Creates a TempProfile instance.
  !!* @param self TempProfile instane on return.
  !!* @param tempMethods The annealing method for each intervall.
  !!* @param tempInts The length of the intervalls (in steps)
  !!* @param tempValues Target temperature for each intervall. This temperature
  !!*   will be reached after the specified number of steps, using the
  !!*   specified profile (constant, linear, exponential)
  subroutine TempProfile_init(self, tempMethods, tempInts, tempValues)
    type(OTempProfile), intent(out) :: self
    integer, intent(in) :: tempMethods(:)
    integer, intent(in) :: tempInts(:)
    real(dp), intent(in) :: tempValues(:)

    integer :: ii, iTmp

    @:ASSERT(all(tempMethods == constProf .or. tempMethods == linProf&
        & .or. tempMethods == expProf))
    @:ASSERT(size(tempInts) > 0)
    @:ASSERT(size(tempInts) == size(tempValues) .and. size(tempInts) == size(tempMethods))
    @:ASSERT(all(tempInts >= 0))
    @:ASSERT(all(tempValues >= 0.0_dp))

    self%nInt = size(tempInts)
    allocate(self%tempInts(0:self%nInt))
    allocate(self%tempValues(0:self%nInt))
    allocate(self%tempMethods(self%nInt))
    self%tempInts(0) = 0
    self%tempInts(1:) = tempInts(:)
    self%tempValues(0) = startingTemp_
    self%tempValues(1:) = tempValues(:)
    self%tempMethods(:) = tempMethods(:)
    iTmp = self%tempInts(1)
    do ii = 2, self%nInt
      iTmp = iTmp + self%tempInts(ii)
      self%tempInts(ii) = iTmp
    end do
    self%incr = 0
    self%iStep = 1
    self%iInt = 1
    do while (self%tempInts(self%iInt) == 0)
      self%iInt = self%iInt + 1
    end do
    self%curTemp = self%tempValues(self%iInt)

  end subroutine TempProfile_init


  !!* Changes the temperature to the next value.
  !!* @param self Pointer to the TempProfile object.
  subroutine TempProfile_next(self)
    type(OTempProfile), intent(inout) :: self

    real(dp) :: subVal, supVal
    integer :: sub, sup
    logical :: tChanged

    self%iStep = self%iStep + 1
    if (self%iStep > self%tempInts(self%nInt)) then
      return
    end if
    !! Looking for the next interval which contains the relevant information
    tChanged = .false.
    do while (self%tempInts(self%iInt) < self%iStep)
      self%iInt = self%iInt + 1
      tChanged = .true.
    end do
    sup = self%tempInts(self%iInt)
    sub = self%tempInts(self%iInt-1)
    supVal = self%tempValues(self%iInt)
    subVal = self%tempValues(self%iInt-1)

    select case (self%tempMethods(self%iInt))
    case (constProf)
      self%curTemp = self%tempValues(self%iInt)
    case (linProf)
      if (tChanged) then
        self%incr = (supVal - subVal) / real(sup - sub, dp)
      end if
      self%curTemp = subVal + self%incr * real(self%iStep - sub, dp)
    case (expProf)
      if (tChanged) then
        self%tempValues(self%iInt) = supVal
        self%tempValues(self%iInt-1) = subVal
        self%incr = log(supVal/subVal) / real(sup - sub, dp)
      end if
      self%curTemp = subVal * exp(self%incr * real(self%iStep - sub, dp))
    end select

  end subroutine TempProfile_next



  !!* Returns the current temperature.
  !!* @param self Pointer to the TempProfile object.
  !!* @param temp Temperature on return.
  subroutine TempProfile_getTemperature(self, temp)
    type(OTempProfile), intent(in) :: self
    real(dp), intent(out) :: temp

    temp = self%curTemp

  end subroutine TempProfile_getTemperature


end module tempprofile
