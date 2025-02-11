!--------------------------------------------------------------------------------------------------!
!  DFTB+: general package for performing fast atomistic simulations                                !
!  Copyright (C) 2006 - 2025  DFTB+ developers group                                               !
!                                                                                                  !
!  See the LICENSE file for terms of usage and distribution.                                       !
!--------------------------------------------------------------------------------------------------!

#:include "fytest.fypp"

#:block TEST_SUITE("file")
  use dftbp_common_file, only : TFileDescr, TOpenOptions, fileExists, openFile, closeFile,&
      & clearFile, defaultBinaryAccess, defaultTextAccess
  use dftbp_io_charmanip, only : tolower
  implicit none

  ! Equality for TOpenOptions instances for easier testing
  interface operator(==)
    module procedure TOpenOptions_isEqual_
  end interface

  ! Unequality for TOpenOptions instances for easier testing
  interface operator(/=)
    module procedure TOpenOptions_isNotEqual_
  end interface

#:contains

  #:block TEST_FIXTURE("openOptions")

    type(TOpenOptions) :: opts
    integer :: ioStat
    character(:), allocatable :: ioMsg

  #:contains

    #:block TEST("equal")
      @:ASSERT(opts == TOpenOptions())
      @:ASSERT(opts /= TOpenOptions(access="stream"))
      @:ASSERT(opts /= TOpenOptions(action="write"))
      @:ASSERT(opts /= TOpenOptions(form="unformatted"))
      @:ASSERT(opts /= TOpenOptions(status="old"))
      @:ASSERT(opts /= TOpenOptions(position="rewind"))
    #:endblock


    #:block TEST("mode_r")
      call opts%setMode("r", ioStat, ioMsg)
      @:ASSERT(ioStat == 0 .and. .not. allocated(ioMsg))
      @:ASSERT(opts == TOpenOptions(action="read", form="formatted", position="rewind",&
          & status="old"))
    #:endblock


    #:block TEST("mode_rt")
      call opts%setMode("rt", ioStat, ioMsg)
      @:ASSERT(ioStat == 0 .and. .not. allocated(ioMsg))
      @:ASSERT(opts == TOpenOptions(action="read", form="formatted", position="rewind",&
          & status="old"))
    #:endblock


    #:block TEST("mode_rb")
      call opts%setMode("rb", ioStat, ioMsg)
      @:ASSERT(ioStat == 0 .and. .not. allocated(ioMsg))
      @:ASSERT(opts == TOpenOptions(action="read", form="unformatted", position="rewind",&
          & status="old"))
    #:endblock


    #:block TEST("mode_rp")
      call opts%setMode("r+", ioStat, ioMsg)
      @:ASSERT(ioStat == 0 .and. .not. allocated(ioMsg))
      @:ASSERT(opts == TOpenOptions(action="readwrite", form="formatted", position="rewind",&
          & status="old"))
    #:endblock


    #:block TEST("mode_rpt")
      call opts%setMode("r+t", ioStat, ioMsg)
      @:ASSERT(ioStat == 0 .and. .not. allocated(ioMsg))
      @:ASSERT(opts == TOpenOptions(action="readwrite", form="formatted", position="rewind",&
          & status="old"))
    #:endblock


    #:block TEST("mode_rpb")
      call opts%setMode("r+b", ioStat, ioMsg)
      @:ASSERT(ioStat == 0 .and. .not. allocated(ioMsg))
      @:ASSERT(opts == TOpenOptions(action="readwrite", form="unformatted", position="rewind",&
          & status="old"))
    #:endblock


    #:block TEST("mode_w")
      call opts%setMode("w", ioStat, ioMsg)
      @:ASSERT(ioStat == 0 .and. .not. allocated(ioMsg))
      @:ASSERT(opts == TOpenOptions(action="write", form="formatted", position="rewind",&
          & status="replace"))
    #:endblock


    #:block TEST("mode_wt")
      call opts%setMode("wt", ioStat, ioMsg)
      @:ASSERT(ioStat == 0 .and. .not. allocated(ioMsg))
      @:ASSERT(opts == TOpenOptions(action="write", form="formatted", position="rewind",&
          & status="replace"))
    #:endblock


    #:block TEST("mode_wb")
      call opts%setMode("wb", ioStat, ioMsg)
      @:ASSERT(ioStat == 0 .and. .not. allocated(ioMsg))
      @:ASSERT(opts == TOpenOptions(action="write", form="unformatted", position="rewind",&
          & status="replace"))
    #:endblock


    #:block TEST("mode_wp")
      call opts%setMode("w+", ioStat, ioMsg)
      @:ASSERT(ioStat == 0 .and. .not. allocated(ioMsg))
      @:ASSERT(opts == TOpenOptions(action="readwrite", form="formatted", position="rewind",&
          & status="replace"))
    #:endblock


    #:block TEST("mode_wpt")
      call opts%setMode("w+t", ioStat, ioMsg)
      @:ASSERT(ioStat == 0 .and. .not. allocated(ioMsg))
      @:ASSERT(opts == TOpenOptions(action="readwrite", form="formatted", position="rewind",&
          & status="replace"))
    #:endblock


    #:block TEST("mode_wpb")
      call opts%setMode("w+b", ioStat, ioMsg)
      @:ASSERT(ioStat == 0 .and. .not. allocated(ioMsg))
      @:ASSERT(opts == TOpenOptions(action="readwrite", form="unformatted", position="rewind",&
          & status="replace"))
    #:endblock


    #:block TEST("mode_a")
      call opts%setMode("a", ioStat, ioMsg)
      @:ASSERT(ioStat == 0 .and. .not. allocated(ioMsg))
      @:ASSERT(opts == TOpenOptions(action="write", form="formatted", position="append",&
          & status="oldnew"))
    #:endblock


    #:block TEST("mode_at")
      call opts%setMode("at", ioStat, ioMsg)
      @:ASSERT(ioStat == 0 .and. .not. allocated(ioMsg))
      @:ASSERT(opts == TOpenOptions(action="write", form="formatted", position="append",&
          & status="oldnew"))
    #:endblock


    #:block TEST("mode_ab")
      call opts%setMode("ab", ioStat, ioMsg)
      @:ASSERT(ioStat == 0 .and. .not. allocated(ioMsg))
      @:ASSERT(opts == TOpenOptions(action="write", form="unformatted", position="append",&
          & status="oldnew"))
    #:endblock


    #:block TEST("mode_ap")
      call opts%setMode("a+", ioStat, ioMsg)
      @:ASSERT(ioStat == 0 .and. .not. allocated(ioMsg))
      @:ASSERT(opts == TOpenOptions(action="readwrite", form="formatted", position="append",&
          & status="oldnew"))
    #:endblock


    #:block TEST("mode_apt")
      call opts%setMode("a+t", ioStat, ioMsg)
      @:ASSERT(ioStat == 0 .and. .not. allocated(ioMsg))
      @:ASSERT(opts == TOpenOptions(action="readwrite", form="formatted", position="append",&
          & status="oldnew"))
    #:endblock


    #:block TEST("mode_apb")
      call opts%setMode("a+b", ioStat, ioMsg)
      @:ASSERT(ioStat == 0 .and. .not. allocated(ioMsg))
      @:ASSERT(opts == TOpenOptions(action="readwrite", form="unformatted", position="append",&
          & status="oldnew"))
    #:endblock


    #:block TEST("mode_invalid_action1")
      call opts%setMode("g", ioStat, ioMsg)
      @:ASSERT(ioStat == -1 .and. allocated(ioMsg))
      @:ASSERT(ioMsg(1:14) == "Invalid action")
    #:endblock


    #:block TEST("mode_invalid_action2")
      call opts%setMode("gt", ioStat, ioMsg)
      @:ASSERT(ioStat == -1 .and. allocated(ioMsg))
      @:ASSERT(ioMsg(1:14) == "Invalid action")
    #:endblock


    #:block TEST("mode_invalid_action3")
      call opts%setMode("g+t", ioStat, ioMsg)
      @:ASSERT(ioStat == -1 .and. allocated(ioMsg))
      @:ASSERT(ioMsg(1:14) == "Invalid action")
    #:endblock


    #:block TEST("mode_invalid_format1")
      call opts%setMode("rp", ioStat, ioMsg)
      @:ASSERT(ioStat == -2 .and. allocated(ioMsg))
      @:ASSERT(ioMsg(1:14) == "Invalid format")
    #:endblock


    #:block TEST("mode_invalid_format2")
      call opts%setMode("r+p", ioStat, ioMsg)
      @:ASSERT(ioStat == -2 .and. allocated(ioMsg))
      @:ASSERT(ioMsg(1:14) == "Invalid format")
    #:endblock


    #:block TEST("mode_superfluous_chars1")
      call opts%setMode("rtq", ioStat, ioMsg)
      @:ASSERT(ioStat == -3 .and. allocated(ioMsg))
      @:ASSERT(ioMsg(1:22) == "Superfluous characters")
    #:endblock


    #:block TEST("mode_superfluous_chars2")
      call opts%setMode("r+tp", ioStat, ioMsg)
      @:ASSERT(ioStat == -3 .and. allocated(ioMsg))
      @:ASSERT(ioMsg(1:22) == "Superfluous characters")
    #:endblock

  #:endblock TEST_FIXTURE


  #:block TEST_FIXTURE("fileDescr")

    integer, parameter :: lineLength = 10
    character(lineLength), parameter :: dummyContent(*) = &
        & [character(lineLength) :: "1st line", "2nd line"]
    type(TFileDescr) :: fd
    integer :: ioStat

  #:contains

    #:block TEST("close_on_exit")
      character(*), parameter :: fname = "file_close_on_exit.tmp"

      call deleteFile_(fname)
      @:ASSERT(.not. fileExists(fname))
      block
        type(TFileDescr) :: fd
        @:ASSERT(.not. connected_(fname))
        call openFile(fd, fname, mode="w", ioStat=ioStat)
        @:ASSERT(ioStat == 0 .and. fd%unit /= -1 .and. connected_(fname))
      end block
      @:ASSERT(fileExists(fname) .and. .not. connected_(fname))
    #:endblock


    #:block TEST("close_on_exit_rank_1")
      character(31), parameter :: fnames(*) = [&
        & "file_close_on_exit_rank_1.a.tmp", "file_close_on_exit_rank_1.b.tmp"]
      integer :: iFile

      do iFile = 1, size(fnames)
        call deleteFile_(fnames(iFile))
        @:ASSERT(.not. fileExists(fnames(iFile)))
      end do
      block
        type(TFileDescr) :: fds(size(fnames))
        do iFile = 1, size(fnames)
          @:ASSERT(.not. connected_(fnames(iFile)))
          call openFile(fds(iFile), fnames(iFile), mode="w", ioStat=ioStat)
          @:ASSERT(ioStat == 0)
        end do
        do iFile = 1, size(fnames)
          @:ASSERT(fds(iFile)%unit /= -1 .and. connected_(fnames(iFile)))
        end do
      end block
      do iFile = 1, size(fnames)
        @:ASSERT(fileExists(fnames(iFile)) .and. .not. connected_(fnames(iFile)))
      end do
    #:endblock


    #:block TEST("close_on_close")
      character(*), parameter :: fname = "file_close_on_exit.tmp"

      call deleteFile_(fname)
      @:ASSERT(.not. fileExists(fname) .and. .not. connected_(fname))
      call openFile(fd, fname, mode="w", ioStat=ioStat)
      @:ASSERT(ioStat == 0 .and. fd%unit /= -1 .and. connected_(fname))
      call closeFile(fd)
      @:ASSERT(fileExists(fname) .and. .not. connected_(fname))
    #:endblock


    #:block TEST("close_on_close_rank_1")
      character(32), parameter :: fnames(*) = [&
        & "file_close_on_exit_rank_1.a.tmp", "file_close_on_exit_rank_1.b.tmp"]
      type(TFileDescr) :: fds(size(fnames))
      integer :: iFile

      do iFile = 1, size(fnames)
        call deleteFile_(fnames(iFile))
        @:ASSERT(.not. fileExists(fnames(iFile)))
      end do
      do iFile = 1, size(fnames)
        @:ASSERT(.not. connected_(fnames(iFile)))
        call openFile(fds(iFile), fnames(iFile), mode="w", ioStat=ioStat)
        @:ASSERT(ioStat == 0)
      end do
      do iFile = 1, size(fnames)
        @:ASSERT(fds(iFile)%unit /= -1 .and. connected_(fnames(iFile)))
      end do
      call closeFile(fds)
      do iFile = 1, size(fnames)
        @:ASSERT(fileExists(fnames(iFile)) .and. .not. connected_(fnames(iFile)))
      end do
    #:endblock


    #:block TEST("open_r")
      character(*), parameter :: fname = "file_open_r.tmp"

      call createTextFile_(fname, dummyContent)
      @:ASSERT(fileExists(fname))
      call openFile(fd, fname, mode="r", ioStat=ioStat)
      @:ASSERT(ioStat == 0 .and. connected_(fname))
      @:ASSERT(checkUnitProperties_(fd%unit, action="read", form="formatted", position="rewind"))
    #:endblock


    #:block TEST("open_rb")
      character(*), parameter :: fname = "file_open_rb.tmp"

      call createTextFile_(fname, dummyContent)
      @:ASSERT(fileExists(fname))
      call openFile(fd, fname, mode="rb", ioStat=ioStat)
      @:ASSERT(ioStat == 0 .and. connected_(fname))
      @:ASSERT(checkUnitProperties_(fd%unit, action="read", form="unformatted", position="rewind"))
    #:endblock


    #:block TEST("open_rp")
      character(*), parameter :: fname = "file_open_rp.tmp"

      call createTextFile_(fname, dummyContent)
      @:ASSERT(fileExists(fname))
      call openFile(fd, fname, mode="r+")
      @:ASSERT(connected_(fname))
      @:ASSERT(checkUnitProperties_(fd%unit, action="readwrite", form="formatted",&
          & position="rewind"))
    #:endblock


    #:block TEST("open_w")
      character(*), parameter :: fname = "file_open_w.tmp"

      call openFile(fd, fname, mode="w", ioStat=ioStat)
      @:ASSERT(ioStat == 0 .and. connected_(fname))
      @:ASSERT(checkUnitProperties_(fd%unit, action="write", form="formatted", position="rewind"))
    #:endblock


    #:block TEST("open_wb")
      character(*), parameter :: fname = "file_open_wb.tmp"

      call openFile(fd, fname, mode="wb", ioStat=ioStat)
      @:ASSERT(ioStat == 0 .and. connected_(fname))
      @:ASSERT(checkUnitProperties_(fd%unit, action="write", form="unformatted", position="rewind"))
    #:endblock


    #:block TEST("open_wp")
      character(*), parameter :: fname = "file_open_wp.tmp"

      call openFile(fd, fname, mode="w+")
      @:ASSERT(connected_(fname))
      @:ASSERT(checkUnitProperties_(fd%unit, action="readwrite", form="formatted",&
          & position="rewind"))
    #:endblock


    #:block TEST("open_a_new")
      character(*), parameter :: fname = "file_open_a_existing.tmp"

      call deleteFile_(fname)
      call openFile(fd, fname, mode="a", ioStat=ioStat)
      @:ASSERT(ioStat == 0 .and. connected_(fname))
      @:ASSERT(checkUnitProperties_(fd%unit, action="write", form="formatted", position="append"))
    #:endblock


    #:block TEST("open_a_existing")
      character(*), parameter :: fname = "file_open_a_existing.tmp"

      call clearFile(fname)
      call openFile(fd, fname, mode="a", ioStat=ioStat)
      @:ASSERT(ioStat == 0 .and. connected_(fname))
      @:ASSERT(checkUnitProperties_(fd%unit, action="write", form="formatted", position="append"))
    #:endblock


    #:block TEST("open_ab")
      character(*), parameter :: fname = "file_open_ab.tmp"

      call openFile(fd, fname, mode="ab", ioStat=ioStat)
      @:ASSERT(ioStat == 0 .and. connected_(fname))
      @:ASSERT(checkUnitProperties_(fd%unit, action="write", form="unformatted", position="append"))
    #:endblock


    #:block TEST("open_ap")
      character(*), parameter :: fname = "file_open_ap.tmp"

      call openFile(fd, fname, mode="a+")
      @:ASSERT(connected_(fname))
      @:ASSERT(checkUnitProperties_(fd%unit, action="readwrite", form="formatted",&
          & position="append"))
    #:endblock


    #:block TEST("open_r_fail")
      character(*), parameter :: fname = "file_open_r_fail.tmp"

      @:ASSERT(.not. fileExists(fname))
      call openFile(fd, fname, mode="r", ioStat=ioStat)
      @:ASSERT((ioStat > 0 .and. fd%unit == -1 .and. .not. connected_(fname)))

    #:endblock

  #:endblock TEST_FIXTURE


  #:block TEST_FIXTURE("clearFile")

  #:contains

    #:block TEST("non_existing")
      character(*), parameter :: fname = "file_non_existing.tmp"

      call deleteFile_(fname)
      @:ASSERT(.not. fileExists(fname))
      call clearFile(fname)
      @:ASSERT(fileExists(fname) .and. fileSize_(fname) == 0)

    #:endblock


    #:block TEST("existing")
      character(*), parameter :: fname = "file_existing.tmp"
      character(1), parameter :: dummyContent(*) = ["a", "b"]

      call createTextFile_(fname, dummyContent)
      @:ASSERT(fileExists(fname) .and. fileSize_(fname) > 0)
      call clearFile(fname)
      @:ASSERT(fileExists(fname) .and. fileSize_(fname) == 0)

    #:endblock

  #:endblock


  #:block TEST_FIXTURE("fileAccess")

    type(TFileDescr) :: fd
    type(TOpenOptions) :: opts
    integer :: ioStat

  #:contains

    #:block TEST("default_text_access")
      character(*), parameter :: fname = "default_text_access"

      ! Note: openFile() closes the file connected to the descriptor (if any)
      call openFile(fd, fname, mode="w", ioStat=ioStat)
      @:ASSERT(ioStat == 0)
      @:ASSERT(checkUnitProperties_(fd%unit, access=defaultTextAccess(2)))
      call openFile(fd, fname, mode="r", ioStat=ioStat)
      @:ASSERT(ioStat == 0)
      @:ASSERT(checkUnitProperties_(fd%unit, access=defaultTextAccess(1)))
      call openFile(fd, fname, mode="r+", ioStat=ioStat)
      @:ASSERT(ioStat == 0)
      @:ASSERT(checkUnitProperties_(fd%unit, access=defaultTextAccess(3)))
    #:endblock


    #:block TEST("default_binary_access")
      character(*), parameter :: fname = "default_binary_access"

      ! Note: openFile() closes the file connected to the descriptor (if any)
      call openFile(fd, fname, mode="wb", ioStat=ioStat)
      @:ASSERT(ioStat == 0)
      @:ASSERT(checkUnitProperties_(fd%unit, access=defaultBinaryAccess(2)))
      call openFile(fd, fname, mode="rb", ioStat=ioStat)
      @:ASSERT(ioStat == 0)
      @:ASSERT(checkUnitProperties_(fd%unit, access=defaultBinaryAccess(1)))
      call openFile(fd, fname, mode="r+b", ioStat=ioStat)
      @:ASSERT(ioStat == 0)
      @:ASSERT(checkUnitProperties_(fd%unit, access=defaultBinaryAccess(3)))
    #:endblock


    #:block TEST("stream_access")
      character(*), parameter :: fname = "stream_access"
      character(*), parameter :: access = "stream"

      ! Note: openFile() closes the file connected to the descriptor (if any)
      opts%access = access
      call openFile(fd, fname, options=opts, mode="w", ioStat=ioStat)
      @:ASSERT(ioStat == 0)
      @:ASSERT(checkUnitProperties_(fd%unit, access=access))
      call openFile(fd, fname, options=opts, mode="r", ioStat=ioStat)
      @:ASSERT(ioStat == 0)
      @:ASSERT(checkUnitProperties_(fd%unit, access=access))
      call openFile(fd, fname, options=opts, mode="r+", ioStat=ioStat)
      @:ASSERT(ioStat == 0)
      @:ASSERT(checkUnitProperties_(fd%unit, access=access))
    #:endblock


    #:block TEST("sequential_access")
      character(*), parameter :: fname = "sequential_access"
      character(*), parameter :: access = "sequential"

      ! Note: openFile() closes the file connected to the descriptor (if any)
      opts%access = access
      call openFile(fd, fname, options=opts, mode="w", ioStat=ioStat)
      @:ASSERT(ioStat == 0)
      @:ASSERT(checkUnitProperties_(fd%unit, access=access))
      call openFile(fd, fname, options=opts, mode="r", ioStat=ioStat)
      @:ASSERT(ioStat == 0)
      @:ASSERT(checkUnitProperties_(fd%unit, access=access))
      call openFile(fd, fname, options=opts, mode="r+", ioStat=ioStat)
      @:ASSERT(ioStat == 0)
      @:ASSERT(checkUnitProperties_(fd%unit, access=access))
    #:endblock

  #:endblock


  ! Compares, whether two options are equal
  function TOpenOptions_isEqual_(this, other) result(isEqual)
    type(TOpenOptions), intent(in) :: this, other
    logical :: isEqual

    isEqual = this%access == other%access .and. this%action == other%action&
        & .and. this%form == other%form .and. this%position == other%position&
        & .and. this%status == other%status

  end function TOpenOptions_isEqual_


  ! Compares, whether two options are not equal
  function TOpenOptions_isNotEqual_(this, other) result(isNotEqual)
    type(TOpenOptions), intent(in) :: this, other
    logical :: isNotEqual

    isNotEqual = .not. (this == other)

  end function TOpenOptions_isNotEqual_


  ! Checks, whether a file is already connected to a unit
  function connected_(fname) result(isConnected)
    character(*), intent(in) :: fname
    logical :: isConnected

    inquire(file=fname, opened=isConnected)

  end function connected_


  ! Creates a test file with specified content
  subroutine createTextFile_(fname, lines)
    character(*), intent(in) :: fname
    character(*), intent(in) :: lines(:)

    integer :: unit
    integer :: iLine

    open(newunit=unit, file=fname, form="formatted", action="write", status="replace")
    do iLine = 1, size(lines)
      write(unit, "(a)") trim(lines(iLine))
    end do
    close(unit)

  end subroutine createTextFile_


  ! Deletes a file if present.
  subroutine deleteFile_(fname)
    character(*), intent(in) :: fname

    integer :: unit
    logical :: exists

    inquire(file=fname, exist=exists)
    if (exists) then
      open(newunit=unit, file=fname)
      close(unit, status="delete")
    end if

  end subroutine deleteFile_


  ! Checks whether unit possesses certain properties
  function checkUnitProperties_(unit, access, action, form, position) result(isOk)
    integer, intent(in) :: unit
    character(*), optional, intent(in) :: access, action, form, position
    logical :: isOk

    character(20) :: access_, action_, form_, position_

    inquire(unit=unit, access=access_, action=action_, form=form_, position=position_)
    isOK = .true.
    if (present(access)) isOK = isOK .and. tolower(access) == tolower(access_)
    if (present(action)) isOK = isOK .and. tolower(action) == tolower(action_)
    if (present(form)) isOK = isOK .and. tolower(form) == tolower(form_)
    if (present(position)) isOK = isOK .and. tolower(position) == tolower(position_)

  end function checkUnitProperties_


  ! Checks whether file possesses certain properties
  function fileSize_(file) result(fileSize)
    character(*), intent(in) :: file
    integer :: fileSize

    inquire(file=file, size=fileSize)

  end function fileSize_

#:endblock TEST_SUITE


@:TEST_DRIVER()
