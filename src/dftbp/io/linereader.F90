!--------------------------------------------------------------------------------------------------!
!  DFTB+: general package for performing fast atomistic simulations                                !
!  Copyright (C) 2006 - 2023  DFTB+ developers group                                               !
!                                                                                                  !
!  See the LICENSE file for terms of usage and distribution.                                       !
!--------------------------------------------------------------------------------------------------!

!> Contains a line reader implementation
module dftbp_io_linereader
  implicit none

  private
  public :: TLineReader


  !> Simple line reader returning content of a formatted file line by line
  type :: TLineReader
    private
    integer :: unit = -1
    integer :: ioStat = -1
    integer :: bufferSize = 4096
  contains
    procedure :: readLine => TLineReader_readLine
  end type TLineReader

  !> Constructor interface
  interface TLineReader
    module procedure TLineReader_construct
  end interface TLineReader


contains

  !> Constructs a line reader instance.
  function TLineReader_construct(unit, bufferSize) result(this)

    !> File unit to read the lines from, file must be opened with formatted I/O and PAD=YES.
    integer, intent(in) :: unit

    !> Size of the buffer to use (default: 4096 bytes)
    integer, optional, intent(in) :: bufferSize

    !> Initialized instance
    type(TLineReader) :: this

    this%unit = unit
    if (present(bufferSize)) this%bufferSize = bufferSize
    this%ioStat = 0

  end function TLineReader_construct


  !> Reads a line form the file.
  !>
  subroutine TLineReader_readLine(this, line, ioStat)

    !> Instance
    class(TLineReader), intent(inout) :: this

    !> Contains the line read
    !>
    !> Note: argument may be unallocated on return if the returned ioStat value is non-zero.
    !>
    character(:), allocatable, intent(out) :: line

    !> I/O status value
    !>
    !> If non-zero iostat is reported, subsequent calls to this routine will always report the same
    !> value (and return no line content any more).
    !>
    !> Note: Whenever valid content is returned in line, the returned ioStat will be zero. When no
    !> content can be returned any more (line is unallocated on return), either end-of-file will be
    !> returned or a positive value indicating the I/O error occured during read.
    !>
    integer, intent(out) :: ioStat

    character(this%bufferSize) :: buffer
    integer :: nReadChars

    if (this%ioStat /= 0) then
      ioStat = this%ioStat
      return
    end if
    do
      read(this%unit, "(a)", advance="no", iostat=ioStat, size=nReadChars) buffer
      if (ioStat > 0) then
        this%ioStat = ioStat
        return
      end if
      if (nReadChars > 0) then
        if (allocated(line)) then
          line = line // buffer(1:nReadChars)
        else
          line = buffer(1:nReadChars)
        end if
      end if
      if (ioStat < 0) then
        if (is_iostat_eor(ioStat)) then
          ! End of record signalizes end of line, so we report no-error to the caller
          this%ioStat = 0
          ioStat = 0
        else if (is_iostat_end(ioStat) .and. nReadChars > 0) then
          ! There are compilers, where non-advancing I/O reads valid characters and signalizes
          ! at the same time end-of-file (typically at last line of files with missing trailing
          ! newline char) -> Report OK to caller, but remember eof for the next call.
          this%ioStat = ioStat
          ioStat = 0
        else
          this%ioStat = ioStat
        end if
        exit
      end if
    end do
    ! Make sure, line is always allocated if returned ioStat is zero
    if (ioStat == 0 .and. .not. allocated(line)) line = ""

  end subroutine TLineReader_readLine


end module dftbp_io_linereader
