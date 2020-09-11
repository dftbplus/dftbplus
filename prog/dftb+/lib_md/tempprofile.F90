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
  use dftbp_charmanip, only : tolower
  implicit none
  private

  public :: TTempProfile, TempProfile_init, identifyTempProfile

  ! Internal constants for the different profiles:
  !> Constant temperature
  integer, parameter :: constProf = 1
  !> linear change in profile
  integer, parameter :: linProf = 2
  !> exponentially changing profile
  integer, parameter :: expProf = 3


  !> Data for the temperature profile.
  type TTempProfile

    private

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

  contains

    procedure :: next
    procedure :: getTemperature

  end type TTempProfile

  !> Default starting temperature
  real(dp), parameter :: startingTemp_ = minTemp

contains


  !> Creates a TempProfile instance.
  subroutine TempProfile_init(this, tempMethods, tempInts, tempValues)

    !> TempProfile instane on return.
    type(TTempProfile), intent(out) :: this

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

    this%nInt = size(tempInts)
    allocate(this%tempInts(0:this%nInt))
    allocate(this%tempValues(0:this%nInt))
    allocate(this%tempMethods(this%nInt))
    this%tempInts(0) = 0
    this%tempInts(1:) = tempInts(:)
    this%tempValues(1:) = tempValues(:)
    this%tempMethods(:) = tempMethods(:)
    iTmp = this%tempInts(1)
    do ii = 2, this%nInt
      iTmp = iTmp + this%tempInts(ii)
      this%tempInts(ii) = iTmp
    end do
    this%incr = 0
    this%iStep = 1
    this%iInt = 1
    do while (this%tempInts(this%iInt) == 0)
      this%iInt = this%iInt + 1
    end do
    if (this%tempMethods(this%iInt) == constProf) then
      this%tempValues(0) = this%tempValues(1)
    else
      this%tempValues(0) = startingTemp_
    end if
    this%curTemp = this%tempValues(0)

  end subroutine TempProfile_init


  !> Changes the temperature to the next value.
  subroutine next(this)

    !> The TempProfile object.
    class(TTempProfile), intent(inout) :: this

    real(dp) :: subVal, supVal
    integer :: sub, sup

    this%iStep = this%iStep + 1
    if (this%iStep > this%tempInts(this%nInt)) then
      return
    end if
    ! Looking for the next interval which contains the relevant information
    do while (this%tempInts(this%iInt) < this%iStep)
      this%iInt = this%iInt + 1
    end do
    sup = this%tempInts(this%iInt)
    sub = this%tempInts(this%iInt-1)
    supVal = this%tempValues(this%iInt)
    subVal = this%tempValues(this%iInt-1)

    select case (this%tempMethods(this%iInt))
    case (constProf)
      this%curTemp = this%tempValues(this%iInt)
    case (linProf)
      this%incr = (supVal - subVal) / real(sup - sub, dp)
      this%curTemp = subVal + this%incr * real(this%iStep - sub, dp)
    case (expProf)
      this%tempValues(this%iInt) = supVal
      this%tempValues(this%iInt-1) = subVal
      this%incr = log(supVal/subVal) / real(sup - sub, dp)
      this%curTemp = subVal * exp(this%incr * real(this%iStep - sub, dp))
    end select

  end subroutine next


  !> Returns the current temperature.
  subroutine getTemperature(this, temp)

    !> Pointer to the TempProfile object.
    class(TTempProfile), intent(in) :: this

    !> Temperature on return.
    real(dp), intent(out) :: temp

    temp = this%curTemp

  end subroutine getTemperature


  !> Maps a (supported) profile name onto integer identifier
  subroutine identifyTempProfile(iProfile, profileName, success)

    !> Internal profile identifying number
    integer, intent(out) :: iProfile

    !> Possible profile name
    character(*), intent(in) :: profileName

    !> was the profile correctly identified
    logical, intent(out) :: success

    success = .true.
    select case (tolower(trim(profileName)))
    case ("constant")
      iProfile = constProf
    case ("linear")
      iProfile = linProf
    case ("exponential")
      iProfile = expProf
    case default
      success = .false.
    end select

  end subroutine identifyTempProfile

end module dftbp_tempprofile
