!--------------------------------------------------------------------------------------------------!
!  DFTB+: general package for performing fast atomistic simulations                                !
!  Copyright (C) 2006 - 2023  DFTB+ developers group                                               !
!                                                                                                  !
!  See the LICENSE file for terms of usage and distribution.                                       !
!--------------------------------------------------------------------------------------------------!

#:include "common.fypp"

!> Contains a file descriptor and methods to set the default file access type globally.
module dftbp_common_file
  use dftbp_io_charmanip, only : i2c
  use dftbp_io_message, only : error
#:block DEBUG_CODE
  use dftbp_common_globalenv, only : stdOut0
#:endblock
  implicit none

  private
  public :: TFileDescr, TOpenOptions, openFile, closeFile, clearFile, fileExists
  public :: defaultFileAccess, setDefaultFileAccess, fileAccessValues


  !> Length of the character string holding a file option value
  integer, parameter :: openOptionCharLen_ = 20


  !> File opening options
  type :: TOpenOptions

    !> File access ("sequential", "direct", "stream" or "default"), default: "default"
    !> "default" indicates to use the global default value defined in defaultFileAccess
    character(openOptionCharLen_) :: access = "default"

    !> File action ("read", "write", "readwrite"), default: "read"
    character(openOptionCharLen_) :: action = "read"

    !> File format ("formatted", "unformatted"), default: "formatted"
    character(openOptionCharLen_) :: form = "formatted"

    !> File status ("unknown", "old", "replace", "oldnew"), default: "unknown"
    !> Oldnew is an extension of the standard values. If the file exists, it will be opened,
    !> if it does not, it will be newly created.
    character(openOptionCharLen_) :: status = "unknown"

    !> File position ("asis", "rewind", "append"), default: "asis"
    character(openOptionCharLen_) :: position = "asis"

  contains

    procedure :: setMode => TOpenOptions_setMode

  end type TOpenOptions


  !> Implements a file descriptor object.
  type :: TFileDescr

    !> File unit
    !>
    !> The file unit should be treated as read-only field and never changed manually.
    !> Also, never call any Fortran statements, which change the unit association status, such
    !> as open() or close(). If you wish to close the file, deallocate the file descriptor instance.
    !>
    integer :: unit = -1

    !> The file the descriptor is associated with (unallocated if not associated)
    character(:), allocatable :: file

    ! Whether unit must be closed on finalization
    logical, private :: needsClosing_ = .false.

  contains
    procedure :: connectToFile => TFileDescr_connectToFile
    procedure :: connectToUnit => TFileDescr_connectToUnit
    procedure :: disconnect => TFileDescr_disconnect
    procedure :: isConnected => TFileDescr_isConnected
    procedure, pass(rhs) :: TFileDescr_assign
    generic :: assignment(=) => TFileDescr_assign
    final :: TFileDescr_final_
  end type TFileDescr


  !> Valid values for file acess
  character(*), parameter :: fileAccessValues(*) =&
      & [character(openOptionCharLen_) :: "sequential", "stream"]


  !> Current default file access type for read, write and readwrite actions
  character(openOptionCharLen_), protected :: defaultFileAccess(3) = &
      & [character(openOptionCharLen_) :: "stream", "stream", "stream"]

contains


  !> Sets the file opening options according to a C-style mode string
  subroutine TOpenOptions_setMode(this, mode, ioStat, ioMsg)

    !> Instance.
    class(TOpenOptions), intent(inout) :: this

    !> C-like mode string ("r", "r+", "w", "w+", "a", "a+" followed by optional "t" or "b")
    character(*), intent(in) :: mode

    !> Zero if the mode string could be interpreted or a negative value otherwise
    integer, intent(out) :: ioStat

    !> Error string generated if some error occured
    character(:), allocatable, optional, intent(out) :: ioMsg

    integer :: pos
    character(2) :: action
    character(1) :: form
    character(:), allocatable :: modestr

    ioStat = 0
    modestr = trim(mode)
    action = ""
    pos = 1
    if (len(modestr) >= pos) then
      action = modestr(pos:pos)
      pos = pos + 1
    end if
    if (len(modestr) >= pos) then
      if (modestr(pos:pos) == "+") then
        action = trim(action) // modestr(pos:pos)
        pos = pos + 1
      end if
    end if

    select case (action)
    case ("r")
      this%action = "read"
      this%status = "old"
      this%position = "rewind"
    case ("r+")
      this%action = "readwrite"
      this%status = "old"
      this%position = "rewind"
    case ("w")
      this%action = "write"
      this%status = "replace"
      this%position = "rewind"
    case ("w+")
      this%action = "readwrite"
      this%status = "replace"
      this%position = "rewind"
    case ("a")
      this%action = "write"
      this%status = "oldnew"
      this%position = "append"
    case ("a+")
      this%action = "readwrite"
      this%status = "oldnew"
      this%position = "append"
    case default
      ioStat = -1
      if (present(ioMsg)) then
        ioMsg = "Invalid action '" // trim(action) // "' in mode specification '" // modestr // "'"
      end if
      return
    end select

    if (len(modestr) >= pos) then
      form = modestr(pos:pos)
      pos = pos + 1
    else
      form = "t"
    end if

    select case (form)
    case ("t")
      this%form = "formatted"
    case ("b")
      this%form = "unformatted"
    case default
      ioStat = -2
      if (present(ioMsg)) then
        ioMsg = "Invalid format '" // trim(form) // "' in mode specification '" // modestr // "'"
      end if
      return
    end select

    if (len(modestr) >= pos) then
      ioStat = -3
      if (present(ioMsg)) then
        ioMsg = "Superfluous characters '" // modestr(pos:) // "' in mode specification '" &
            & // modestr // "'"
      end if
      return
    end if

  end subroutine TOpenOptions_setMode


  !> Opens a file and returns the descriptor.
  subroutine TFileDescr_connectToFile(this, file, options, mode, ioStat, ioMsg)

    !> File descriptor on exit (invalid if the file could not be opened and ioStat /= 0)
    class(TFileDescr), intent(out) :: this

    !> Name of the file to open
    character(*), intent(in) :: file

    !> File opening options
    type(TOpenOptions), optional, intent(in) :: options

    !> C-style mode specification for file opening options with following possible values:
    !> * "r": read (file must exist, positioned at start),
    !> * "r+": read/write (file must exist, positioned at start),
    !> * "w": write (file created or truncated if it already exists)
    !> * "w+": readwrite (file created or truncated if it already exists)
    !> * "a": append-write (file opened if exists, otherwise created, positioned at the end)
    !> * "a+": append-read/write (file opened if exists, otherwise created, positioned at the end)
    !> The values above can be followed by "t" for text/unformatted mode (default) or "b" for
    !> binary/unformatted mode.
    !>
    !> When arguments options and mode are both specified, the resulting options are determined by
    !> applying options first and then set the fields which mode manipulates.
    !>
    character(*), optional, intent(in) :: mode

    !> I/O stat error generated during open, zero on exit, if no error occured.
    integer, optional, intent(out) :: ioStat

    !> I/O stat message generated during open, unallocated on exit, if no error occured.
    character(:), allocatable, optional, intent(out) :: ioMsg

    type(TOpenOptions) :: opts
    integer :: ioStat_
    character(1024) :: ioMsg_

    if (present(options)) then
      opts = options
    end if

    if (present(mode)) then
      call opts%setMode(mode, ioStat_, ioMsg)
      if (ioStat_ /= 0) then
        if (present(ioStat)) then
          ioStat = ioStat_
          return
        else
          call error("Invalid mode specification '" // trim(mode) // "'")
        end if
      end if
    end if

    if (opts%access == "default") then
      select case (opts%action)
      case ("read")
        opts%access = defaultFileAccess(1)
      case ("write")
        opts%access = defaultFileAccess(2)
      case ("readwrite")
        opts%access = defaultFileAccess(3)
      end select
    end if

    if (opts%status == "oldnew") then
      if (fileExists(file)) then
        opts%status = "old"
      else
        opts%status = "new"
      end if
    end if

    open(newunit=this%unit, file=file, access=opts%access, action=opts%action, form=opts%form,&
        & status=opts%status, position=opts%position, iostat=ioStat_, iomsg=ioMsg_)

    if (ioStat_ /= 0) then
      this%unit = -1
      if (present(ioStat) .or. present(ioMsg)) then
        if (present(ioStat)) then
          ioStat = ioStat_
        end if
        if (present(ioMsg)) then
          ioMsg = trim(ioMsg_)
        end if
        return
      else
        call error("Failed to open file '" // trim(file) // "' [(" // i2c(ioStat_) // ") "&
            &  // trim(ioMsg_) // "]")
      end if
    end if
    if (present(ioStat)) then
      ioStat = 0
    end if
    this%file = trim(file)
    this%needsClosing_ = .true.

  #:block DEBUG_CODE
    write(stdOut0, "(5a)") "[Debug] File '", trim(file), "' opened (action='", trim(opts%action),&
        & "')"
    flush(stdOut0)
  #:endblock DEBUG_CODE

  end subroutine TFileDescr_connectToFile


  !> Connects the descriptor to an existing file unit
  subroutine TFileDescr_connectToUnit(this, unit)

    !> File descriptor on exit
    class(TFileDescr), intent(out) :: this

    !> Unit number to wrap as file (unit will NOT be closed on finalization)
    integer, intent(in) :: unit

    this%unit = unit
    this%file = ""

  end subroutine TFileDescr_connectToUnit


  !> Finalizes (closes) a file
  elemental impure subroutine TFileDescr_disconnect(this)

    !> Instance
    class(TFileDescr), intent(inout) :: this

    if (this%unit /= -1 .and. this%needsClosing_) then
      close(this%unit)
    #:block DEBUG_CODE
      write(stdOut0, "(3a)") "[Debug] File '", trim(this%file), "' closed"
      flush(stdOut0)
    #:endblock DEBUG_CODE
    end if
    this%unit = -1
    if (allocated(this%file)) then
      deallocate(this%file)
    end if
    this%needsClosing_ = .false.

  end subroutine TFileDescr_disconnect


  !> Queries, whether the descriptor is connected to a file or unit.
  elemental function TFileDescr_isConnected(this) result(connected)

    !> Instance
    class(TFileDescr), intent(in) :: this

    !> Whether the descriptor is connected to a file/unit
    logical :: connected

    connected = this%unit /= -1

  end function TFileDescr_isConnected


  !> Assignment (should never be called as it stops the code)
  elemental impure subroutine TFileDescr_assign(lhs, rhs)

    !> Right hand side of the assignment
    class(*), intent(inout) :: lhs

    !> Left hand side of the assignment
    class(TFileDescr), intent(in) :: rhs

    call error("Internal error: TFileDescr object is not permitted on the RHS of an assignment")

  end subroutine TFileDescr_assign


  ! Finalizes the instance
  elemental impure subroutine TFileDescr_final_(this)
    type(TFileDescr), intent(inout) :: this

    call this%disconnect()

  end subroutine TFileDescr_final_


  !> Convenience wrapper connecting a file descriptor to a file
  subroutine openFile(fileDescr, file, options, mode, ioStat, ioMsg)

    !> File descriptor connected to the open file
    type(TFileDescr), intent(out) :: fileDescr

    !> Name of the file to open
    character(*), intent(in) :: file

    !> Fortran style file opening option specification
    type(TOpenOptions), optional, intent(in) :: options

    !> C-style mode specification for file opening options with following possible values:
    !> * "r": read (file must exist, positioned at start),
    !> * "r+": read/write (file must exist, positioned at start),
    !> * "w": write (file created or truncated if it already exists)
    !> * "w+": readwrite (file created or truncated if it already exists)
    !> * "a": append-write (file opened if exists, otherwise created, positioned at the end)
    !> * "a+": append-read/write (file opened if exists, otherwise created, positioned at the end)
    !> The values above can be followed by "t" for text/unformatted mode (default) or "b" for
    !> binary/unformatted mode.
    !>
    !> When arguments options and mode are both specified, the resulting options are determined by
    !> applying options first and then set the fields which mode manipulates.
    !>
    character(*), optional, intent(in) :: mode

    !> I/O stat error generated during open
    integer, optional, intent(out) :: ioStat

    !> I/O message if error occured during open, unallocated otherwise
    character(:), allocatable, optional, intent(out) :: ioMsg

    call fileDescr%connectToFile(file, options=options, mode=mode, ioStat=ioStat, ioMsg=ioMsg)

  end subroutine openFile


  !> Convenience wrapper for disconnecting a file descriptor (which closes the file)
  elemental impure subroutine closeFile(fileDescr)

    !> File descriptor, will be unconnected on exit
    type(TFileDescr), intent(inout) :: fileDescr

    call fileDescr%disconnect()

  end subroutine closeFile


  !> Sets the default access type for file opening operations
  subroutine setDefaultFileAccess(readAccess, writeAccess, readwriteAccess)

    !> Access type to use for read access ("sequential", "direct", "stream")
    character(*), intent(in) :: readAccess

    !> Access type to use for write access ("sequential", "direct", "stream")
    character(*), intent(in) :: writeAccess

    !> Access type to use for readwrite access ("sequential", "direct", "stream")
    character(*), intent(in) :: readwriteAccess

    @:ASSERT(any(readAccess == fileAccessValues))
    @:ASSERT(any(writeAccess == fileAccessValues))
    @:ASSERT(any(readwriteAccess == fileAccessValues))

    defaultFileAccess(:) = [readAccess, writeAccess, readwriteAccess]

  end subroutine setDefaultFileAccess


  !> Checks, whether a certain file exists in the file system.
  function fileExists(file) result(exists)

    !> File to check for existence
    character(*), intent(in) :: file

    !> True, if the file exists.
    logical :: exists

    inquire(file=file, exist=exists)

  end function fileExists


  !> Creates empty file without content or truncates files to zero content if it already exists
  subroutine clearFile(file, ioStat, ioMsg)

    !> Name of the file to create
    character(*), intent(in) :: file

    !> I/O stat error generated during open, zero on exit, if no error occured.
    integer, optional, intent(out) :: ioStat

    !> I/O stat message generated during open, unallocated on exit, if no error occured.
    character(:), allocatable, optional, intent(out) :: ioMsg

    type(TFileDescr) :: fd

    call openFile(fd, file, mode="w", ioStat=ioStat, ioMsg=ioMsg)
    call closeFile(fd)

  end subroutine clearFile

end module dftbp_common_file
