!--------------------------------------------------------------------------------------------------!
!  DFTB+: general package for performing fast atomistic simulations                                !
!  Copyright (C) 2006 - 2021  DFTB+ developers group                                               !
!                                                                                                  !
!  See the LICENSE file for terms of usage and distribution.                                       !
!--------------------------------------------------------------------------------------------------!

#:include 'common.fypp'

!> File system utilities
module dftbp_common_filesystem
  implicit none
  private

  public :: getEnvVar, isAbsolutePath


contains


  !> Check whether a file path is absolute
  function isAbsolutePath(path) result(absolute)

    !> Name of variable
    character(len=*), intent(in) :: path

    !> Absolute file path
    logical :: absolute

    character(len=*), parameter :: chars = "ABCDEFGHIJKLMNOPQRSTUVWXYZ"

    absolute = index(path, "/") == 1 .or. index(path, ":\") == 2 .and. scan(path, chars) == 1
  end function isAbsolutePath


  !> Obtain the value of an environment variable
  subroutine getEnvVar(var, val)

    !> Name of variable
    character(len=*), intent(in) :: var

    !> Value of variable
    character(len=:), allocatable, intent(out) :: val

    integer :: length, stat

    call get_environment_variable(var, length=length, status=stat)
    if (stat /= 0) then
      return
    endif

    allocate(character(len=length) :: val, stat=stat)
    if (stat /= 0) then
      return
    endif

    if (length > 0) then
      call get_environment_variable(var, val, status=stat)
      if (stat /= 0) then
        deallocate(val)
        return
      end if
    end if
  end subroutine getEnvVar


end module dftbp_common_filesystem
