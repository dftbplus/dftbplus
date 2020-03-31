!--------------------------------------------------------------------------------------------------!
!  DFTB+: general package for performing fast atomistic simulations                                !
!  Copyright (C) 2006 - 2020  DFTB+ developers group                                               !
!                                                                                                  !
!  See the LICENSE file for terms of usage and distribution.                                       !
!--------------------------------------------------------------------------------------------------!

#:include 'common.fypp'

!> Contains simple temperature profiles for molecular dynamics.
module dftbp_tempprofile
  use dftbp_assert
  use dftbp_accuracy
  implicit none
  private

  public :: TTempProfile, init, next, getTemperature
  public :: constProf, linProf, expProf


  !> Data for the temperature profile.
  type TTempProfile

    !> The annealing method for each interval.
    integer, allocatable :: tempMethods(:)

    !> The length of the intervals (in MD steps)
    integer, allocatable :: tempInts(:)

    !> Target temperature for each interval
    real(dp), allocatable :: tempValues(:)

    !> Current kinetic temperature
    real(dp) :: curTemp

    !> Index for current temperature value
    integer :: iInt

    !> Number of intervals in total
    integer :: nInt

    !> Current interval
    integer :: iStep

    !> Temperature increment to next step
    real(dp) :: incr
  end type TTempProfile


  !> Initialise the profile
  interface init
    module procedure TempProfile_init
  end interface init


  !> Next temperature in the profile
  interface next
    module procedure TempProfile_next
  end interface next


  !> Get the current temperature in the profile
  interface getTemperature
    module procedure TempProfile_getTemperature
  end interface getTemperature

  ! Constants for the different profile options

  !> Constant temperature
  integer, parameter :: constProf = 1

  !> linear change in profile
  integer, parameter :: linProf = 2

  !> exponentially changing profile
  integer, parameter :: expProf = 3


  !> Default starting temperature
  real(dp), parameter :: startingTemp_ = minTemp

contains


  !> Creates a TempProfile instance.
  subroutine TempProfile_init(self, tempMethods, tempInts, tempValues)

    !> TempProfile instane on return.
    type(TTempProfile), intent(out) :: self

    !> The annealing method for each interval.
    integer, intent(in) :: tempMethods(:)

    !> The length of the intervals (in steps)
    integer, intent(in) :: tempInts(:)

    !> Target temperature for each interval. This temperature will be reached after the specified
    !> number of steps, using the specified profile (constant, linear, exponential)
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


  !> Changes the temperature to the next value.
  subroutine TempProfile_next(self)

    !> The TempProfile object.
    type(TTempProfile), intent(inout) :: self

    real(dp) :: subVal, supVal
    integer :: sub, sup
    logical :: tChanged

    self%iStep = self%iStep + 1
    if (self%iStep > self%tempInts(self%nInt)) then
      return
    end if
    ! Looking for the next interval which contains the relevant information
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


  !> Returns the current temperature.
  subroutine TempProfile_getTemperature(self, temp)

    !> Pointer to the TempProfile object.
    type(TTempProfile), intent(in) :: self

    !> Temperature on return.
    real(dp), intent(out) :: temp

    temp = self%curTemp

  end subroutine TempProfile_getTemperature

end module dftbp_tempprofile
