!--------------------------------------------------------------------------------------------------!
!  DFTB+: general package for performing fast atomistic simulations                                !
!  Copyright (C) 2006 - 2025  DFTB+ developers group                                               !
!                                                                                                  !
!  See the LICENSE file for terms of usage and distribution.                                       !
!--------------------------------------------------------------------------------------------------!

#:include 'common.fypp'

!> Contains simple temperature profiles for molecular dynamics.
module dftbp_md_tempprofile
  use dftbp_common_accuracy, only : dp, minTemp
  use dftbp_io_charmanip, only : tolower
  implicit none

  private
  public :: TTempProfile, TempProfile_init, identifyTempProfile
  public :: TTempProfileInput
  public :: tempProfileTypes


  !> Helper type for defining temperature profile types
  type :: TTempProfileTypes

    !> Constant temperature
    integer :: constant = 1

    !> Linear temperature change
    integer :: linear = 2

    !> Exponential temperature change
    integer :: exponential = 3

  end type TTempProfileTypes


  !> Temperature profile types
  type(TTempProfileTypes), parameter :: tempProfileTypes = TTempProfileTypes()


  !> Input for defining a temperature profiles
  type :: TTempProfileInput

    !> The annealing method for each interval.
    integer, allocatable :: tempMethods(:)

    !> The length of the intervals (in MD steps)
    integer, allocatable :: tempInts(:)

    !> Target temperature for each interval (zeroth entry is dummy with staring temperature)
    real(dp), allocatable :: tempValues(:)

  end type TTempProfileInput


  !> Data for the temperature profile.
  type :: TTempProfile
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
  subroutine TempProfile_init(this, input)

    !> TempProfile instane on return.
    type(TTempProfile), intent(out) :: this

    !> Input data
    type(TTempProfileInput), intent(in) :: input

    integer :: ii, iTmp

    @:ASSERT(all(input%tempMethods == tempProfileTypes%constant&
        & .or. input%tempMethods == tempProfileTypes%linear&
        & .or. input%tempMethods == tempProfileTypes%exponential))
    @:ASSERT(size(input%tempInts) > 0)
    @:ASSERT(size(input%tempInts) == size(input%tempValues))
    @:ASSERT(size(input%tempInts) == size(input%tempMethods))
    @:ASSERT(all(input%tempInts >= 0))
    @:ASSERT(all(input%tempValues >= 0.0_dp))


    this%tempInts = [0, input%tempInts]
    this%tempValues = [startingTemp_, input%tempValues]
    this%tempMethods = [tempProfileTypes%constant, input%tempMethods]
    ! Convert temperature interval lengths to absolute step numbers
    iTmp = this%tempInts(1)
    do ii = 2, size(this%tempInts)
      iTmp = iTmp + this%tempInts(ii)
      this%tempInts(ii) = iTmp
    end do
    this%incr = 0
    this%iStep = 1
    this%iInt = 1
    do while (this%tempInts(this%iInt) == 0)
      this%iInt = this%iInt + 1
    end do
    if (this%tempMethods(this%iInt) == tempProfileTypes%constant) then
      this%tempValues(this%iInt - 1) = this%tempValues(this%iInt)
    end if
    this%curTemp = this%tempValues(this%iInt - 1)

  end subroutine TempProfile_init


  !> Changes the temperature to the next value.
  subroutine next(this)

    !> The TempProfile object.
    class(TTempProfile), intent(inout) :: this

    real(dp) :: subVal, supVal
    integer :: sub, sup

    this%iStep = this%iStep + 1
    if (this%iStep > this%tempInts(size(this%tempInts))) then
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
    case (tempProfileTypes%constant)
      this%curTemp = this%tempValues(this%iInt)
    case (tempProfileTypes%linear)
      this%incr = (supVal - subVal) / real(sup - sub, dp)
      this%curTemp = subVal + this%incr * real(this%iStep - sub, dp)
    case (tempProfileTypes%exponential)
      this%tempValues(this%iInt) = supVal
      this%tempValues(this%iInt-1) = subVal
      this%incr = log(supVal/subVal) / real(sup - sub, dp)
      this%curTemp = subVal * exp(this%incr * real(this%iStep - sub, dp))
    end select
    print *, "TempProfile:", this%iStep, this%curTemp

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
      iProfile = tempProfileTypes%constant
    case ("linear")
      iProfile = tempProfileTypes%linear
    case ("exponential")
      iProfile = tempProfileTypes%exponential
    case default
      success = .false.
    end select

  end subroutine identifyTempProfile

end module dftbp_md_tempprofile
