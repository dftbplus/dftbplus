!--------------------------------------------------------------------------------------------------!
!  DFTB+: general package for performing fast atomistic simulations                                !
!  Copyright (C) 2006 - 2025  DFTB+ developers group                                               !
!                                                                                                  !
!  See the LICENSE file for terms of usage and distribution.                                       !
!--------------------------------------------------------------------------------------------------!

#:include "fortuno_serial.fypp"

module test_common_file
  use dftbp_common_file, only : TFileDescr, TOpenOptions, fileExists, openFile, closeFile,&
      & clearFile, defaultBinaryAccess, defaultTextAccess
  use dftbp_io_charmanip, only : tolower
  use fortuno_serial, only : suite => serial_suite_item, test_list
  $:FORTUNO_SERIAL_IMPORTS()
  implicit none

  private
  public :: tests


  ! Equality for TOpenOptions instances for easier testing
  interface operator(==)
    module procedure TOpenOptions_isEqual_
  end interface

  ! Unequality for TOpenOptions instances for easier testing
  interface operator(/=)
    module procedure TOpenOptions_isNotEqual_
  end interface


  type :: TOpenOptionsFx
    type(TOpenOptions) :: opts
    integer :: ioStat
    character(:), allocatable :: ioMsg
  end type TOpenOptionsFx


  type :: TFileDescrFx
    integer :: lineLength = 10
    character(10) :: dummyContent(2) = [character(10) :: "1st line", "2nd line"]
    type(TFileDescr) :: fd
    integer :: ioStat
  end type TFileDescrFx


  type :: TFileAccessFx
    type(TFileDescr) :: fd
    type(TOpenOptions) :: opts
    integer :: ioStat
  end type TFileAccessFx

contains


  $:TEST("equal", label="openOptions")
    type(TOpenOptionsFx) :: fx
    @:ASSERT(fx%opts == TOpenOptions())
    @:ASSERT(fx%opts /= TOpenOptions(access="stream"))
    @:ASSERT(fx%opts /= TOpenOptions(action="write"))
    @:ASSERT(fx%opts /= TOpenOptions(form="unformatted"))
    @:ASSERT(fx%opts /= TOpenOptions(status="old"))
    @:ASSERT(fx%opts /= TOpenOptions(position="rewind"))
  $:END_TEST()


  $:TEST("mode_r", label="openOptions")
    type(TOpenOptionsFx) :: fx
    call fx%opts%setMode("r", fx%ioStat, fx%ioMsg)
    @:ASSERT(fx%ioStat == 0 .and. .not. allocated(fx%ioMsg))
    @:ASSERT(fx%opts == TOpenOptions(action="read", form="formatted", position="rewind",&
        & status="old"))
  $:END_TEST()


  $:TEST("mode_rt", label="openOptions")
    type(TOpenOptionsFx) :: fx
    call fx%opts%setMode("rt", fx%ioStat, fx%ioMsg)
    @:ASSERT(fx%ioStat == 0 .and. .not. allocated(fx%ioMsg))
    @:ASSERT(fx%opts == TOpenOptions(action="read", form="formatted", position="rewind",&
        & status="old"))
  $:END_TEST()


  $:TEST("mode_rb", label="openOptions")
    type(TOpenOptionsFx) :: fx
    call fx%opts%setMode("rb", fx%ioStat, fx%ioMsg)
    @:ASSERT(fx%ioStat == 0 .and. .not. allocated(fx%ioMsg))
    @:ASSERT(fx%opts == TOpenOptions(action="read", form="unformatted", position="rewind",&
        & status="old"))
  $:END_TEST()


  $:TEST("mode_rp", label="openOptions")
    type(TOpenOptionsFx) :: fx
    call fx%opts%setMode("r+", fx%ioStat, fx%ioMsg)
    @:ASSERT(fx%ioStat == 0 .and. .not. allocated(fx%ioMsg))
    @:ASSERT(fx%opts == TOpenOptions(action="readwrite", form="formatted", position="rewind",&
        & status="old"))
  $:END_TEST()


  $:TEST("mode_rpt", label="openOptions")
    type(TOpenOptionsFx) :: fx
    call fx%opts%setMode("r+t", fx%ioStat, fx%ioMsg)
    @:ASSERT(fx%ioStat == 0 .and. .not. allocated(fx%ioMsg))
    @:ASSERT(fx%opts == TOpenOptions(action="readwrite", form="formatted", position="rewind",&
        & status="old"))
  $:END_TEST()


  $:TEST("mode_rpb", label="openOptions")
    type(TOpenOptionsFx) :: fx
    call fx%opts%setMode("r+b", fx%ioStat, fx%ioMsg)
    @:ASSERT(fx%ioStat == 0 .and. .not. allocated(fx%ioMsg))
    @:ASSERT(fx%opts == TOpenOptions(action="readwrite", form="unformatted", position="rewind",&
        & status="old"))
  $:END_TEST()


  $:TEST("mode_w", label="openOptions")
    type(TOpenOptionsFx) :: fx
    call fx%opts%setMode("w", fx%ioStat, fx%ioMsg)
    @:ASSERT(fx%ioStat == 0 .and. .not. allocated(fx%ioMsg))
    @:ASSERT(fx%opts == TOpenOptions(action="write", form="formatted", position="rewind",&
        & status="replace"))
  $:END_TEST()


  $:TEST("mode_wt", label="openOptions")
    type(TOpenOptionsFx) :: fx
    call fx%opts%setMode("wt", fx%ioStat, fx%ioMsg)
    @:ASSERT(fx%ioStat == 0 .and. .not. allocated(fx%ioMsg))
    @:ASSERT(fx%opts == TOpenOptions(action="write", form="formatted", position="rewind",&
        & status="replace"))
  $:END_TEST()


  $:TEST("mode_wb", label="openOptions")
    type(TOpenOptionsFx) :: fx
    call fx%opts%setMode("wb", fx%ioStat, fx%ioMsg)
    @:ASSERT(fx%ioStat == 0 .and. .not. allocated(fx%ioMsg))
    @:ASSERT(fx%opts == TOpenOptions(action="write", form="unformatted", position="rewind",&
        & status="replace"))
  $:END_TEST()


  $:TEST("mode_wp", label="openOptions")
    type(TOpenOptionsFx) :: fx
    call fx%opts%setMode("w+", fx%ioStat, fx%ioMsg)
    @:ASSERT(fx%ioStat == 0 .and. .not. allocated(fx%ioMsg))
    @:ASSERT(fx%opts == TOpenOptions(action="readwrite", form="formatted", position="rewind",&
        & status="replace"))
  $:END_TEST()


  $:TEST("mode_wpt", label="openOptions")
    type(TOpenOptionsFx) :: fx
    call fx%opts%setMode("w+t", fx%ioStat, fx%ioMsg)
    @:ASSERT(fx%ioStat == 0 .and. .not. allocated(fx%ioMsg))
    @:ASSERT(fx%opts == TOpenOptions(action="readwrite", form="formatted", position="rewind",&
        & status="replace"))
  $:END_TEST()


  $:TEST("mode_wpb", label="openOptions")
    type(TOpenOptionsFx) :: fx
    call fx%opts%setMode("w+b", fx%ioStat, fx%ioMsg)
    @:ASSERT(fx%ioStat == 0 .and. .not. allocated(fx%ioMsg))
    @:ASSERT(fx%opts == TOpenOptions(action="readwrite", form="unformatted", position="rewind",&
        & status="replace"))
  $:END_TEST()


  $:TEST("mode_a", label="openOptions")
    type(TOpenOptionsFx) :: fx
    call fx%opts%setMode("a", fx%ioStat, fx%ioMsg)
    @:ASSERT(fx%ioStat == 0 .and. .not. allocated(fx%ioMsg))
    @:ASSERT(fx%opts == TOpenOptions(action="write", form="formatted", position="append",&
        & status="oldnew"))
  $:END_TEST()


  $:TEST("mode_at", label="openOptions")
    type(TOpenOptionsFx) :: fx
    call fx%opts%setMode("at", fx%ioStat, fx%ioMsg)
    @:ASSERT(fx%ioStat == 0 .and. .not. allocated(fx%ioMsg))
    @:ASSERT(fx%opts == TOpenOptions(action="write", form="formatted", position="append",&
        & status="oldnew"))
  $:END_TEST()


  $:TEST("mode_ab", label="openOptions")
    type(TOpenOptionsFx) :: fx
    call fx%opts%setMode("ab", fx%ioStat, fx%ioMsg)
    @:ASSERT(fx%ioStat == 0 .and. .not. allocated(fx%ioMsg))
    @:ASSERT(fx%opts == TOpenOptions(action="write", form="unformatted", position="append",&
        & status="oldnew"))
  $:END_TEST()


  $:TEST("mode_ap", label="openOptions")
    type(TOpenOptionsFx) :: fx
    call fx%opts%setMode("a+", fx%ioStat, fx%ioMsg)
    @:ASSERT(fx%ioStat == 0 .and. .not. allocated(fx%ioMsg))
    @:ASSERT(fx%opts == TOpenOptions(action="readwrite", form="formatted", position="append",&
        & status="oldnew"))
  $:END_TEST()


  $:TEST("mode_apt", label="openOptions")
    type(TOpenOptionsFx) :: fx
    call fx%opts%setMode("a+t", fx%ioStat, fx%ioMsg)
    @:ASSERT(fx%ioStat == 0 .and. .not. allocated(fx%ioMsg))
    @:ASSERT(fx%opts == TOpenOptions(action="readwrite", form="formatted", position="append",&
        & status="oldnew"))
  $:END_TEST()


  $:TEST("mode_apb", label="openOptions")
    type(TOpenOptionsFx) :: fx
    call fx%opts%setMode("a+b", fx%ioStat, fx%ioMsg)
    @:ASSERT(fx%ioStat == 0 .and. .not. allocated(fx%ioMsg))
    @:ASSERT(fx%opts == TOpenOptions(action="readwrite", form="unformatted", position="append",&
        & status="oldnew"))
  $:END_TEST()


  $:TEST("mode_invalid_action1", label="openOptions")
    type(TOpenOptionsFx) :: fx
    call fx%opts%setMode("g", fx%ioStat, fx%ioMsg)
    @:ASSERT(fx%ioStat == -1 .and. allocated(fx%ioMsg))
    @:ASSERT(fx%ioMsg(1:14) == "Invalid action")
  $:END_TEST()


  $:TEST("mode_invalid_action2", label="openOptions")
    type(TOpenOptionsFx) :: fx
    call fx%opts%setMode("gt", fx%ioStat, fx%ioMsg)
    @:ASSERT(fx%ioStat == -1 .and. allocated(fx%ioMsg))
    @:ASSERT(fx%ioMsg(1:14) == "Invalid action")
  $:END_TEST()


  $:TEST("mode_invalid_action3", label="openOptions")
    type(TOpenOptionsFx) :: fx
    call fx%opts%setMode("g+t", fx%ioStat, fx%ioMsg)
    @:ASSERT(fx%ioStat == -1 .and. allocated(fx%ioMsg))
    @:ASSERT(fx%ioMsg(1:14) == "Invalid action")
  $:END_TEST()


  $:TEST("mode_invalid_format1", label="openOptions")
    type(TOpenOptionsFx) :: fx
    call fx%opts%setMode("rp", fx%ioStat, fx%ioMsg)
    @:ASSERT(fx%ioStat == -2 .and. allocated(fx%ioMsg))
    @:ASSERT(fx%ioMsg(1:14) == "Invalid format")
  $:END_TEST()


  $:TEST("mode_invalid_format2", label="openOptions")
    type(TOpenOptionsFx) :: fx
    call fx%opts%setMode("r+p", fx%ioStat, fx%ioMsg)
    @:ASSERT(fx%ioStat == -2 .and. allocated(fx%ioMsg))
    @:ASSERT(fx%ioMsg(1:14) == "Invalid format")
  $:END_TEST()


  $:TEST("mode_superfluous_chars1", label="openOptions")
    type(TOpenOptionsFx) :: fx
    call fx%opts%setMode("rtq", fx%ioStat, fx%ioMsg)
    @:ASSERT(fx%ioStat == -3 .and. allocated(fx%ioMsg))
    @:ASSERT(fx%ioMsg(1:22) == "Superfluous characters")
  $:END_TEST()


  $:TEST("mode_superfluous_chars2", label="openOptions")
    type(TOpenOptionsFx) :: fx
    call fx%opts%setMode("r+tp", fx%ioStat, fx%ioMsg)
    @:ASSERT(fx%ioStat == -3 .and. allocated(fx%ioMsg))
    @:ASSERT(fx%ioMsg(1:22) == "Superfluous characters")
  $:END_TEST()


  $:TEST("close_on_exit", label="fileDescr")
    character(*), parameter :: fname = "file_close_on_exit.tmp"
    type(TFileDescrFx) :: fx

    call deleteFile_(fname)
    @:ASSERT(.not. fileExists(fname))
    block
      type(TFileDescr) :: fd
      @:ASSERT(.not. connected_(fname))
      call openFile(fd, fname, mode="w", ioStat=fx%ioStat)
      @:ASSERT(fx%ioStat == 0 .and. fd%unit /= -1 .and. connected_(fname))
    end block
    @:ASSERT(fileExists(fname) .and. .not. connected_(fname))
  $:END_TEST()


  $:TEST("close_on_exit_rank_1", label="fileDescr")
    character(31), parameter :: fnames(*) = [&
        & "file_close_on_exit_rank_1.a.tmp", "file_close_on_exit_rank_1.b.tmp"]
    type(TFileDescrFx) :: fx
    integer :: iFile

    do iFile = 1, size(fnames)
      call deleteFile_(fnames(iFile))
      @:ASSERT(.not. fileExists(fnames(iFile)))
    end do
    block
      type(TFileDescr) :: fds(size(fnames))
      do iFile = 1, size(fnames)
        @:ASSERT(.not. connected_(fnames(iFile)))
        call openFile(fds(iFile), fnames(iFile), mode="w", ioStat=fx%ioStat)
        @:ASSERT(fx%ioStat == 0)
      end do
      do iFile = 1, size(fnames)
        @:ASSERT(fds(iFile)%unit /= -1 .and. connected_(fnames(iFile)))
      end do
    end block
    do iFile = 1, size(fnames)
      @:ASSERT(fileExists(fnames(iFile)) .and. .not. connected_(fnames(iFile)))
    end do
  $:END_TEST()


  $:TEST("close_on_close", label="fileDescr")
    character(*), parameter :: fname = "file_close_on_exit.tmp"
    type(TFileDescrFx) :: fx

    call deleteFile_(fname)
    @:ASSERT(.not. fileExists(fname) .and. .not. connected_(fname))
    call openFile(fx%fd, fname, mode="w", ioStat=fx%ioStat)
    @:ASSERT(fx%ioStat == 0 .and. fx%fd%unit /= -1 .and. connected_(fname))
    call closeFile(fx%fd)
    @:ASSERT(fileExists(fname) .and. .not. connected_(fname))
  $:END_TEST()


  $:TEST("close_on_close_rank_1", label="fileDescr")
    character(32), parameter :: fnames(*) = [&
      & "file_close_on_exit_rank_1.a.tmp", "file_close_on_exit_rank_1.b.tmp"]
    type(TFileDescrFx) :: fx
    type(TFileDescr) :: fds(size(fnames))
    integer :: iFile

    do iFile = 1, size(fnames)
      call deleteFile_(fnames(iFile))
      @:ASSERT(.not. fileExists(fnames(iFile)))
    end do
    do iFile = 1, size(fnames)
      @:ASSERT(.not. connected_(fnames(iFile)))
      call openFile(fds(iFile), fnames(iFile), mode="w", ioStat=fx%ioStat)
      @:ASSERT(fx%ioStat == 0)
    end do
    do iFile = 1, size(fnames)
      @:ASSERT(fds(iFile)%unit /= -1 .and. connected_(fnames(iFile)))
    end do
    call closeFile(fds)
    do iFile = 1, size(fnames)
      @:ASSERT(fileExists(fnames(iFile)) .and. .not. connected_(fnames(iFile)))
    end do
  $:END_TEST()


  $:TEST("open_r", label="fileDescr")
    character(*), parameter :: fname = "file_open_r.tmp"
    type(TFileDescrFx) :: fx

    call createTextFile_(fname, fx%dummyContent)
    @:ASSERT(fileExists(fname))
    call openFile(fx%fd, fname, mode="r", ioStat=fx%ioStat)
    @:ASSERT(fx%ioStat == 0 .and. connected_(fname))
    @:ASSERT(checkUnitProperties_(fx%fd%unit, action="read", form="formatted", position="rewind"))
  $:END_TEST()


  $:TEST("open_rb", label="fileDescr")
    character(*), parameter :: fname = "file_open_rb.tmp"
    type(TFileDescrFx) :: fx

    call createTextFile_(fname, fx%dummyContent)
    @:ASSERT(fileExists(fname))
    call openFile(fx%fd, fname, mode="rb", ioStat=fx%ioStat)
    @:ASSERT(fx%ioStat == 0 .and. connected_(fname))
    @:ASSERT(checkUnitProperties_(fx%fd%unit, action="read", form="unformatted", position="rewind"))
  $:END_TEST()


  $:TEST("open_rp", label="fileDescr")
    character(*), parameter :: fname = "file_open_rp.tmp"
    type(TFileDescrFx) :: fx

    call createTextFile_(fname, fx%dummyContent)
    @:ASSERT(fileExists(fname))
    call openFile(fx%fd, fname, mode="r+")
    @:ASSERT(connected_(fname))
    @:ASSERT(checkUnitProperties_(fx%fd%unit, action="readwrite", form="formatted",&
        & position="rewind"))
  $:END_TEST()


  $:TEST("open_w", label="fileDescr")
    character(*), parameter :: fname = "file_open_w.tmp"
    type(TFileDescrFx) :: fx

    call openFile(fx%fd, fname, mode="w", ioStat=fx%ioStat)
    @:ASSERT(fx%ioStat == 0 .and. connected_(fname))
    @:ASSERT(checkUnitProperties_(fx%fd%unit, action="write", form="formatted", position="rewind"))
  $:END_TEST()


  $:TEST("open_wb", label="fileDescr")
    character(*), parameter :: fname = "file_open_wb.tmp"
    type(TFileDescrFx) :: fx

    call openFile(fx%fd, fname, mode="wb", ioStat=fx%ioStat)
    @:ASSERT(fx%ioStat == 0 .and. connected_(fname))
    @:ASSERT(checkUnitProperties_(fx%fd%unit, action="write", form="unformatted", position="rewind"))
  $:END_TEST()


  $:TEST("open_wp", label="fileDescr")
    character(*), parameter :: fname = "file_open_wp.tmp"
    type(TFileDescrFx) :: fx

    call openFile(fx%fd, fname, mode="w+")
    @:ASSERT(connected_(fname))
    @:ASSERT(checkUnitProperties_(fx%fd%unit, action="readwrite", form="formatted",&
        & position="rewind"))
  $:END_TEST()


  $:TEST("open_a_new", label="fileDescr")
    character(*), parameter :: fname = "file_open_a_existing.tmp"
    type(TFileDescrFx) :: fx

    call deleteFile_(fname)
    call openFile(fx%fd, fname, mode="a", ioStat=fx%ioStat)
    @:ASSERT(fx%ioStat == 0 .and. connected_(fname))
    @:ASSERT(checkUnitProperties_(fx%fd%unit, action="write", form="formatted", position="append"))
  $:END_TEST()


  $:TEST("open_a_existing", label="fileDescr")
    character(*), parameter :: fname = "file_open_a_existing.tmp"
    type(TFileDescrFx) :: fx

    call clearFile(fname)
    call openFile(fx%fd, fname, mode="a", ioStat=fx%ioStat)
    @:ASSERT(fx%ioStat == 0 .and. connected_(fname))
    @:ASSERT(checkUnitProperties_(fx%fd%unit, action="write", form="formatted", position="append"))
  $:END_TEST()


  $:TEST("open_ab", label="fileDescr")
    character(*), parameter :: fname = "file_open_ab.tmp"
    type(TFileDescrFx) :: fx

    call openFile(fx%fd, fname, mode="ab", ioStat=fx%ioStat)
    @:ASSERT(fx%ioStat == 0 .and. connected_(fname))
    @:ASSERT(checkUnitProperties_(fx%fd%unit, action="write", form="unformatted",&
        & position="append"))
  $:END_TEST()


  $:TEST("open_ap", label="fileDescr")
    character(*), parameter :: fname = "file_open_ap.tmp"
    type(TFileDescrFx) :: fx

    call openFile(fx%fd, fname, mode="a+")
    @:ASSERT(connected_(fname))
    @:ASSERT(checkUnitProperties_(fx%fd%unit, action="readwrite", form="formatted",&
        & position="append"))
  $:END_TEST()


  $:TEST("open_r_fail", label="fileDescr")
    character(*), parameter :: fname = "file_open_r_fail.tmp"
    type(TFileDescrFx) :: fx

    @:ASSERT(.not. fileExists(fname))
    call openFile(fx%fd, fname, mode="r", ioStat=fx%ioStat)
    @:ASSERT((fx%ioStat > 0 .and. fx%fd%unit == -1 .and. .not. connected_(fname)))

  $:END_TEST()


  $:TEST("non_existing", label="clearFile")
    character(*), parameter :: fname = "file_non_existing.tmp"

    call deleteFile_(fname)
    @:ASSERT(.not. fileExists(fname))
    call clearFile(fname)
    @:ASSERT(fileExists(fname) .and. fileSize_(fname) == 0)

  $:END_TEST()


  $:TEST("existing", label="clearFile")
    character(*), parameter :: fname = "file_existing.tmp"
    character(1), parameter :: dummyContent(*) = ["a", "b"]

    call createTextFile_(fname, dummyContent)
    @:ASSERT(fileExists(fname) .and. fileSize_(fname) > 0)
    call clearFile(fname)
    @:ASSERT(fileExists(fname) .and. fileSize_(fname) == 0)

  $:END_TEST()


  $:TEST("default_text_access", label="fileAccess")
    character(*), parameter :: fname = "default_text_access"
    type(TFileAccessFx) :: fx

    ! Note: openFile() closes the file connected to the descriptor (if any)
    call openFile(fx%fd, fname, mode="w", ioStat=fx%ioStat)
    @:ASSERT(fx%ioStat == 0)
    @:ASSERT(checkUnitProperties_(fx%fd%unit, access=defaultTextAccess(2)))
    call openFile(fx%fd, fname, mode="r", ioStat=fx%ioStat)
    @:ASSERT(fx%ioStat == 0)
    @:ASSERT(checkUnitProperties_(fx%fd%unit, access=defaultTextAccess(1)))
    call openFile(fx%fd, fname, mode="r+", ioStat=fx%ioStat)
    @:ASSERT(fx%ioStat == 0)
    @:ASSERT(checkUnitProperties_(fx%fd%unit, access=defaultTextAccess(3)))
  $:END_TEST()


  $:TEST("default_binary_access", label="fileAccess")
    character(*), parameter :: fname = "default_binary_access"
    type(TFileAccessFx) :: fx

    ! Note: openFile() closes the file connected to the descriptor (if any)
    call openFile(fx%fd, fname, mode="wb", ioStat=fx%ioStat)
    @:ASSERT(fx%ioStat == 0)
    @:ASSERT(checkUnitProperties_(fx%fd%unit, access=defaultBinaryAccess(2)))
    call openFile(fx%fd, fname, mode="rb", ioStat=fx%ioStat)
    @:ASSERT(fx%ioStat == 0)
    @:ASSERT(checkUnitProperties_(fx%fd%unit, access=defaultBinaryAccess(1)))
    call openFile(fx%fd, fname, mode="r+b", ioStat=fx%ioStat)
    @:ASSERT(fx%ioStat == 0)
    @:ASSERT(checkUnitProperties_(fx%fd%unit, access=defaultBinaryAccess(3)))
  $:END_TEST()


  $:TEST("stream_access", label="fileAccess")
    character(*), parameter :: fname = "stream_access"
    character(*), parameter :: access = "stream"
    type(TFileAccessFx) :: fx

    ! Note: openFile() closes the file connected to the descriptor (if any)
    fx%opts%access = access
    call openFile(fx%fd, fname, options=fx%opts, mode="w", ioStat=fx%ioStat)
    @:ASSERT(fx%ioStat == 0)
    @:ASSERT(checkUnitProperties_(fx%fd%unit, access=access))
    call openFile(fx%fd, fname, options=fx%opts, mode="r", ioStat=fx%ioStat)
    @:ASSERT(fx%ioStat == 0)
    @:ASSERT(checkUnitProperties_(fx%fd%unit, access=access))
    call openFile(fx%fd, fname, options=fx%opts, mode="r+", ioStat=fx%ioStat)
    @:ASSERT(fx%ioStat == 0)
    @:ASSERT(checkUnitProperties_(fx%fd%unit, access=access))
  $:END_TEST()


  $:TEST("sequential_access", label="fileAccess")
    character(*), parameter :: fname = "sequential_access"
    character(*), parameter :: access = "sequential"
    type(TFileAccessFx) :: fx

    ! Note: openFile() closes the file connected to the descriptor (if any)
    fx%opts%access = access
    call openFile(fx%fd, fname, options=fx%opts, mode="w", ioStat=fx%ioStat)
    @:ASSERT(fx%ioStat == 0)
    @:ASSERT(checkUnitProperties_(fx%fd%unit, access=access))
    call openFile(fx%fd, fname, options=fx%opts, mode="r", ioStat=fx%ioStat)
    @:ASSERT(fx%ioStat == 0)
    @:ASSERT(checkUnitProperties_(fx%fd%unit, access=access))
    call openFile(fx%fd, fname, options=fx%opts, mode="r+", ioStat=fx%ioStat)
    @:ASSERT(fx%ioStat == 0)
    @:ASSERT(checkUnitProperties_(fx%fd%unit, access=access))
  $:END_TEST()


  function tests()
    type(test_list) :: tests

    tests = test_list([&
        suite("file", test_list([&
            suite("openOptions", test_list([&
                $:TEST_ITEMS(label="openOptions")
            ])),&
            suite("fileDescr", test_list([&
                $:TEST_ITEMS(label="fileDescr")
            ])),&
            suite("clearFile", test_list([&
                $:TEST_ITEMS(label="clearFile")
            ])),&
            suite("fileAccess", test_list([&
                $:TEST_ITEMS(label="fileAccess")
            ]))&
        ]))&
    ])
    $:STOP_ON_MISSING_TEST_ITEMS()

  end function tests


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

end module test_common_file
