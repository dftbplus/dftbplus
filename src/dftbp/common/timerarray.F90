!--------------------------------------------------------------------------------------------------!
!  DFTB+: general package for performing fast atomistic simulations                                !
!  Copyright (C) 2006 - 2025  DFTB+ developers group                                               !
!                                                                                                  !
!  See the LICENSE file for terms of usage and distribution.                                       !
!--------------------------------------------------------------------------------------------------!

#:include 'common.fypp'

module dftbp_common_timerarray
  use dftbp_common_accuracy, only : dp
  use dftbp_common_timer, only : TTimer
  implicit none

  private
  public :: TTimerItem
  public :: TTimerArray, TTimerArray_init

  !> Initialisation data for a single timer in the array
  type :: TTimerItem
    character(40) :: name
    integer :: level
  end type TTimerItem

  !> Implements global timer with subtimers
  type :: TTimerArray
    private
    type(TTimer) :: myTimer
    character(40), allocatable :: timerNames(:)
    integer, allocatable :: timerLevels(:)
    type(TTimer), allocatable :: timers(:)
    real(dp), allocatable :: cpuTimes(:)
    real(dp), allocatable :: wallClockTimes(:)
    integer :: maxLevel
    character(:), allocatable :: header
    integer :: unit
  contains
    procedure :: startTimer
    procedure :: stopTimer
    procedure :: writeTimings
    procedure :: reset
  end type TTimerArray

contains

  !> Initializes a global timer
  subroutine TTimerArray_init(this, timerItems, unit, maxLevel, header)

    !> Instance
    type(TTimerArray), intent(out) :: this

    !> Names of the sub-timers to use
    type(TTimerItem), intent(in) :: timerItems(:)

    !> File unit to write the statistics to
    integer, intent(in) :: unit

    !> Last timer level to be included in printing (default: all)
    integer, intent(in), optional :: maxLevel

    !> Optional header message for the timings
    character(*), intent(in), optional :: header

    integer :: nTimer

    if (present(maxLevel)) then
      this%maxLevel = maxLevel
    else
      this%maxLevel = -1
    end if

    if (present(header)) then
      this%header = header
    else
      this%header = "Timing"
    end if

    this%unit = unit

    this%timerNames = timerItems(:)%name
    this%timerLevels = timerItems(:)%level
    nTimer = size(timerItems)
    allocate(this%timers(nTimer))
    allocate(this%cpuTimes(nTimer))
    allocate(this%wallClockTimes(nTimer))
    call this%reset()

  end subroutine TTimerArray_init


  !> Resets the timers in the global timer
  subroutine reset(this)

    !> Instance
    class(TTimerArray), intent(inout) :: this

    call this%myTimer%start()
    this%cpuTimes(:) = 0.0_dp
    this%wallClockTimes(:) = 0.0_dp

  end subroutine reset


  !> Starts a given sub-timer.
  subroutine startTimer(this, timerIndex)

    !> Instance.
    class(TTimerArray), intent(inout) :: this

    !> Index of the sub-timer.
    integer, intent(in) :: timerIndex

    @:ASSERT(timerIndex >= 1 .and. timerIndex <= size(this%timers))

    call this%timers(timerIndex)%start()

  end subroutine startTimer


  !> Stops a given sub-timer.
  subroutine stopTimer(this, timerIndex)

    !> Instance.
    class(TTimerArray), intent(inout) :: this

    !> Index of the timer.
    integer, intent(in) :: timerIndex

    @:ASSERT(timerIndex >= 1 .and. timerIndex <= size(this%timers))

    call this%timers(timerIndex)%stop()
    this%cpuTimes(timerIndex) = this%cpuTimes(timerIndex) + this%timers(timerIndex)%getCpuTime()
    this%wallClockTimes(timerIndex) = this%wallClockTimes(timerIndex)&
        & + this%timers(timerIndex)%getWallClockTime()

  end subroutine stopTimer


  !> Writes the current timing values
  subroutine writeTimings(this)

    !> Instance
    class(TTimerArray), intent(in) :: this

    integer :: fp
    real(dp) :: totalCpu, totalWall, cpuTime, wallTime, allCpu, allWall
    integer :: iTimer, level, maxLevel
    character :: operation
    character(:), allocatable :: prefix

    if (this%maxLevel < 0) then
      maxLevel = maxval(this%timerLevels)
    else
      maxLevel = this%maxLevel
    end if
    if (maxLevel < 1) then
      return
    end if

    totalCpu = this%myTimer%getCpuTime()
    totalWall = this%myTimer%getWallClockTime()

    fp = this%unit
    write(fp, *)
    write(fp, "(A)") repeat("-", 80)
    write(fp, "(A,T46,A,T66,A)") this%header, 'cpu [s]', 'wall clock [s]'
    write(fp, "(A)") repeat("-", 80)
    allCpu = 0.0
    allWall = 0.0
    do iTimer = 1, size(this%timers)
      level = this%timerLevels(iTimer)
      if (level > maxLevel) then
        cycle
      end if
      cpuTime = this%cpuTimes(iTimer)
      wallTime = this%wallClockTimes(iTimer)
      if (abs(cpuTime) < 1e-2_dp .and. abs(wallTime) < 1e-2) then
        cycle
      end if
      prefix = repeat(" ", 2 * (level - 1))
      if (level == 1) then
        operation = "+"
      else
        operation = " "
      end if
      write(fp, "(A,A,T40,A,T41,F11.2,1X,'(',F5.1,'%)',T61,F11.2,1X,'(',F5.1,'%)')")&
          & prefix, trim(this%timerNames(iTimer)), operation, cpuTime,&
          & relative_(cpuTime, totalCpu) * 100.0_dp, wallTime,&
          & relative_(wallTime, totalWall) * 100.0_dp
      if (this%timerLevels(iTimer) == 1) then
        allCpu = allCpu + cpuTime
        allWall = allWall + wallTime
      end if
    end do
    write(fp, "(A)") repeat("-", 80)
    write(fp, "(A,T40,A,T41,F11.2,1X,'(',F5.1,'%)',T61,F11.2,1X,'(',F5.1,'%)')")&
        & "Missing", "+", abs(totalCpu - allCpu),&
        & relative_(abs(totalCpu - allCpu), totalCpu) * 100.0_dp,&
        & abs(totalWall - allWall), relative_(abs(totalWall - allWall), totalWall) * 100.0_dp
    write(fp, "(A,T40,A,T41,F11.2,1X,'(',F5.1,'%)',T61,F11.2,1X,'(',F5.1,'%)')")&
        & "Total", "=", totalCpu, 100.0_dp, totalWall, 100.0_dp
    write(fp, "(A)") repeat("-", 80)

  end subroutine writeTimings


  ! Gives the relative share, or 1.0 if the total is too small
  pure function relative_(item, total) result(rel)
    real(dp), intent(in) :: item, total
    real(dp) :: rel

    ! Omit total share calculations for total times shorter than one microsecond
    real(dp), parameter :: tolerance = 1e-6_dp

    if (total >= tolerance) then
      rel = item / total
    else
      rel = 1.0_dp
    end if

  end function relative_

end module dftbp_common_timerarray
