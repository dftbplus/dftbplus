!--------------------------------------------------------------------------------------------------!
!  DFTB+: general package for performing fast atomistic simulations                                !
!  Copyright (C) 2006 - 2023  DFTB+ developers group                                               !
!                                                                                                  !
!  See the LICENSE file for terms of usage and distribution.                                       !
!--------------------------------------------------------------------------------------------------!

#:include 'common.fypp'

!> File system utilities
module dftbp_common_filesystem
  use dftbp_io_charmanip, only : endsWith
  use dftbp_extlibs_xmlf90, only : string, char, assignment(=)
  implicit none
  private

  public :: getEnvVar, isAbsolutePath, joinPaths, findFile, getParamSearchPaths, joinPathsPrettyErr


  !> Whether the operating system we are using is UNIX (FIXME: use preprocessor here instead)
  logical, parameter :: is_unix = .true.

  !> Path separator for this platform
  character(len=*), parameter :: separator = merge("/", "\", is_unix)

  !> Environment variable containing parameter directory for DFTB+
  character(len=*), parameter :: paramEnv = "DFTBPLUS_PARAM_DIR"


contains


  !> Check whether a file path is absolute
  pure function isAbsolutePath(path) result(absolute)

    !> Name of variable
    character(len=*), intent(in) :: path

    !> Absolute file path
    logical :: absolute

    character(len=*), parameter :: chars = &
        & "ABCDEFGHIJKLMNOPQRSTUVWXYZabcdefghijklmnopqrstuvwxyz"

    absolute = index(path, separator) == 1
    if (.not.(is_unix.or.absolute)) then
      absolute = scan(path, chars) == 1 .and. index(path, ":") == 2 .and. scan(path, "\/") == 3
    end if
  end function isAbsolutePath


  !> Obtain the value of an environment variable
  subroutine getEnvVar(var, val)

    !> Name of variable
    character(len=*), intent(in) :: var

    !> Value of variable. Unallocated, if variable could not be retrieved.
    character(len=:), allocatable, intent(out) :: val

    integer :: length, stat

    call get_environment_variable(var, length=length, status=stat)
    if (stat /= 0) then
      return
    end if

    allocate(character(len=length) :: val, stat=stat)
    if (stat /= 0) then
      return
    end if

    if (length > 0) then
      call get_environment_variable(var, val, status=stat)
      if (stat /= 0) then
        deallocate(val)
        return
      end if
    end if
  end subroutine getEnvVar


  !> Join two paths together
  pure function joinPaths(prefix, suffix) result(path)

    !> First path component
    character(len=*), intent(in) :: prefix

    !> Next path component
    character(len=*), intent(in) :: suffix

    !> Joined path
    character(len=:), allocatable :: path

    if (isAbsolutePath(suffix)) then
      path = suffix
    else if (endsWith(prefix, separator)) then
      path = prefix // suffix
    else
      path = prefix // separator // suffix
    end if

  end function joinPaths


  !> Find a file in a PATH-like context
  subroutine findFile(searchPath, inName, outName)

    !> Paths to search in, the current working directory is always searched first
    type(string), intent(in) :: searchPath(:)

    !> File to search for
    character(len=*), intent(in) :: inName

    !> Contains first match on output, unallocated if no file was found
    character(len=:), intent(out), allocatable :: outName

    integer :: iPath
    logical :: exists

    if (isAbsolutePath(inName)) then
      inquire(file=trim(inName), exist=exists)
      if (exists) outName = trim(inName)
      return
    end if

    do iPath = 1, size(searchPath)
      outName = joinPaths(char(searchPath(iPath)), trim(inName))
      inquire(file=outName, exist=exists)
      if (exists) exit
    end do

    if (.not. exists .and. allocated(outName)) deallocate(outName)

  end subroutine findFile


  !> Get the environment variable describing the search path for DFTB+. Additionally, the current
  !! working directory ('.') is pre-pended.
  subroutine getParamSearchPaths(path)

    !> Instance
    type(string), intent(out), allocatable :: path(:)

    character(len=:), allocatable :: var

    call getEnvVar(paramEnv, var)

    if (allocated(var)) then
      allocate(path(2))
      path(2) = var
    else
      allocate(path(1))
    end if

    path(1) = "."

  end subroutine getParamSearchPaths


  !> Joins search paths for pretty error printing.
  function joinPathsPrettyErr(path) result(prettyStr)

    !> Instance
    type(string), intent(in), allocatable :: path(:)

    !! Joined paths for pretty error printing
    character(len=:), allocatable :: prettyStr

    !! Auxiliary variables
    integer :: iPath

    if (allocated(path)) then
      if (size(path) > 0) then
        prettyStr = "'" // char(path(1)) // "'"
      end if
    else
      return
    end if

    do iPath = 2, size(path)
      prettyStr = prettyStr // ", '" // char(path(iPath)) // "'"
    end do

  end function joinPathsPrettyErr

end module dftbp_common_filesystem
