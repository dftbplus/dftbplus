!--------------------------------------------------------------------------------------------------!
!  DFTB+: general package for performing fast atomistic simulations                                !
!  Copyright (C) 2018  DFTB+ developers group                                                      !
!                                                                                                  !
!  See the LICENSE file for terms of usage and distribution.                                       !
!--------------------------------------------------------------------------------------------------!

!> Time stages in the code
module timer
  use globalenv, only : stdOut
  implicit none
  private

  public :: TTimer

  !> Simple timer object
  type :: TTimer
    private
    real :: startTime = 0.0
    real :: endTime = -1.0
    integer :: startCount = 0
    integer :: endCount = -1
    integer :: countRate = 1
  contains
    procedure :: start
    procedure :: stop => stopTimer
    procedure :: getCpuTime
    procedure :: getWallClockTime
    procedure :: writeTimes
  end type TTimer

contains

  !> Starts the timer.
  subroutine start(this)

    !> Instance.
    class(TTimer), intent(inout)   :: this

    call cpu_time(this%startTime)
    call system_clock(count=this%startCount, count_rate=this%countRate)

  end subroutine start


  !> Stops the timer.
  subroutine stopTimer(this)

    !> Instance.
    class(TTimer), intent(inout) :: this

    call cpu_time(this%endTime)
    call system_clock(count=this%endCount)

  end subroutine stopTimer


  !> Returns the measured CPU time
  function getCpuTime(this) result(cpuTime)

    !> Instance.
    class(TTimer), intent(in) :: this

    !> Cpu time evolved between the last start() and stop() calls.
    real :: cpuTime

    real :: endTime

    if (this%endTime < 0.0) then
      call cpu_time(endTime)
    else
      endTime = this%endTime
    end if
    cpuTime = endTime - this%startTime

  end function getCpuTime


  !> Returns the measured wall clock time.
  function getWallClockTime(this) result(wallClockTime)

    !> Instance.
    class(TTimer), intent(in) :: this

    !> Wall clock time evolved between the last start() and stop() calls.
    real :: wallClockTime

    integer :: endCount

    if (this%endCount < 0) then
      call system_clock(count=endCount)
    else
      endCount = this%endCount
    end if
    if (this%countRate == 0) then
      wallClockTime = 0.0
    else
      wallClockTime = real(endCount - this%startCount) / real(this%countRate)
    end if

  end function getWallClockTime


  !> Writes the current measured times.
  subroutine writeTimes(this, msg, node, fp)

    !> Instance.
    class(TTimer), intent(in) :: this

    !> Message to print along the timings
    character(*), intent(in) :: msg

    !> Node information (default: no node information is printed)
    integer, intent(in), optional :: node

    !> File in which to write the times (default: actual standard out)
    integer, intent(in), optional :: fp

    integer :: fp0

    if (present(fp)) then
      fp0 = fp
    else
      fp0 = stdOut
    end if

    if (present(node)) then
      write(fp0, "(A,1X,I5.5,A,1X,A,T48,A,1X,F8.2,5X,A,1X,F8.2)") 'NODE', node, '|TIME',&
          & trim(msg), 'CPU:', this%getCpuTime(), 'WALL:', this%getWallClockTime()
    else
      write(fp0, "(A,1X,A,T48,A,1X,F8.2,5X,A,1X,F8.2)") 'TIME', trim(msg), 'CPU:',&
          & this%getCpuTime(), 'WALL:', this%getWallClockTime()
    end if

  end subroutine writeTimes


end module timer
