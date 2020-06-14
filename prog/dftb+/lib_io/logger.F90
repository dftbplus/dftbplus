!--------------------------------------------------------------------------------------------------!
!  DFTB+: general package for performing fast atomistic simulations                                !
!  Copyright (C) 2006 - 2020  DFTB+ developers group                                               !
!                                                                                                  !
!  See the LICENSE file for terms of usage and distribution.                                       !
!--------------------------------------------------------------------------------------------------!

!> Contains a simple logger which helps to avoid direct write statements.
module dftbp_logger
  use dftbp_accuracy, only : dp
  use dftbp_globalenv, only : stdOut
  use dftbp_optarg
  implicit none
  private

  public :: LogWriter, LogWriter_init


  !> Logger
  type :: LogWriter
    private
    integer :: verbosity
  contains
    ! internal write procedures

    !> write a string
    procedure :: writeStr

    !> write an integer
    procedure :: writeInt

    !> write a real
    procedure :: writeReal

    !> write a real vector
    procedure :: writeReal1

    !> write a real array
    procedure :: writeReal2

    !> generic for all of these
    generic :: write => writeStr, writeInt, writeReal, writeReal1, writeReal2
  end type LogWriter


  !> Constructor for LogWriter
  interface LogWriter
    module procedure construct
  end interface LogWriter


  !> Default verbosity level
  integer, parameter :: DEFAULT_VERBOSITY = 1


  !> Maximum line length
  integer, parameter :: MAX_LEN = 80

contains


  !> Initialises a logger
  subroutine LogWriter_init(this, verbosity)


    !> Initialised instance on exit
    type(LogWriter), intent(out) :: this


    !> verbosity level
    integer, intent(in), optional :: verbosity

    if (present(verbosity)) then
      this%verbosity = verbosity
    else
      this%verbosity = DEFAULT_VERBOSITY
    end if

  end subroutine LogWriter_init


  !> Constructs a logger.
  function construct(verbosity) result(this)


    !> verbosity level
    integer, intent(in), optional :: verbosity


    !> Initialised instance on exit
    type(LogWriter) :: this

    call LogWriter_init(this, verbosity)

  end function construct


  !> Writes a message into the log (string).
  subroutine writeStr(this, msg, verbosity, formStr)


    !> Instance
    class(LogWriter), intent(in) :: this


    !> Message to write
    character(*), intent(in) :: msg


    !> Verbosity level
    integer, intent(in), optional :: verbosity


    !> Format string for a single item
    character(*), intent(in), optional :: formStr

    character(*), parameter :: DEFAULT_FORM_STR = '(A)'
    character(:), allocatable :: formStr0
    integer :: verbosity0

    call getOptionalArg(DEFAULT_VERBOSITY, verbosity0, verbosity)
    call getOptionalArg(DEFAULT_FORM_STR, formStr0, formStr)
    if (verbosity0 <= this%verbosity) then
      write(stdout, formStr0) msg
    end if

  end subroutine writeStr


  !> Writes a message into the log (int).
  subroutine writeInt(this, msg, verbosity, formStr)


    !> Instance
    class(LogWriter), intent(in) :: this


    !> Message to write
    integer, intent(in) :: msg


    !> Verbosity level
    integer, intent(in), optional :: verbosity


    !> Format string for a single item
    character(*), intent(in), optional :: formStr

    character(*), parameter :: DEFAULT_FORM_STR = '(I0)'
    character(:), allocatable :: formStr0
    integer :: verbosity0

    call getOptionalArg(DEFAULT_VERBOSITY, verbosity0, verbosity)
    call getOptionalArg(DEFAULT_FORM_STR, formStr0, formStr)
    if (verbosity0 <= this%verbosity) then
      write(stdout, formStr0) msg
    end if

  end subroutine writeInt


  !> Writes a message into the log (real).
  subroutine writeReal(this, msg, verbosity, formStr)


    !> Instance
    class(LogWriter), intent(in) :: this


    !> Message to write
    real(dp), intent(in) :: msg


    !> Verbosity level
    integer, intent(in), optional :: verbosity


    !> Format string for a single item
    character(*), intent(in), optional :: formStr

    character(*), parameter :: DEFAULT_FORM_STR = '(ES23.15)'
    character(:), allocatable :: formStr0
    integer :: verbosity0

    call getOptionalArg(DEFAULT_VERBOSITY, verbosity0, verbosity)
    call getOptionalArg(DEFAULT_FORM_STR, formStr0, formStr)
    if (verbosity0 <= this%verbosity) then
      write(stdout, formStr0) msg
    end if

  end subroutine writeReal


  !> Writes a message into the log (real1).
  subroutine writeReal1(this, msg, verbosity, formStr, columnwise)


    !> Instance
    class(LogWriter), intent(in) :: this


    !> Message to write
    real(dp), intent(in) :: msg(:)


    !> Verbosity level
    integer, intent(in), optional :: verbosity


    !> Format string for a single item
    character(*), intent(in), optional :: formStr


    !> Whether column vectors should be written columnwise (default: rowwise)
    logical, intent(in), optional :: columnwise

    character(*), parameter :: DEFAULT_FORM_STR = '(ES23.15)'
    character(:), allocatable :: formStr0, formStrRow
    integer :: verbosity0
    logical :: columnwise0

    call getOptionalArg(DEFAULT_VERBOSITY, verbosity0, verbosity)
    call getOptionalArg(DEFAULT_FORM_STR, formStr0, formStr)
    call getOptionalArg(.false., columnwise0, columnwise)

    if (verbosity0 <= this%verbosity) then
      if (columnwise0) then
        call getRowFormat(formStr0, 1, formStrRow)
        write(stdout, formStrRow) msg
      else
        call getRowFormat(formStr0, size(msg), formStrRow)
        write(stdout, formStrRow) msg
      end if
    end if

  end subroutine writeReal1


  !> Writes a message into the log (real2).
  subroutine writeReal2(this, msg, verbosity, formStr, columnwise)


    !> Instance
    class(LogWriter), intent(in) :: this


    !> Message to write
    real(dp), intent(in) :: msg(:,:)


    !> Verbosity level
    integer, intent(in), optional :: verbosity


    !> Format string for a single item
    character(*), intent(in), optional :: formStr


    !> Whether column vectors should be written columnwise (default: rowwise)
    logical, intent(in), optional :: columnwise

    integer :: verbosity0
    logical :: columnwise0
    integer :: ii

    call getOptionalArg(DEFAULT_VERBOSITY, verbosity0, verbosity)
    call getOptionalArg(.false., columnwise0, columnwise)

    if (verbosity0 <= this%verbosity) then
      if (columnwise0) then
        do ii = 1, size(msg, dim=1)
          call this%write(msg(ii,:), verbosity, formStr, .false.)
        end do
      else
        do ii = 1, size(msg, dim=2)
          call this%write(msg(:,ii), verbosity, formStr, .false.)
        end do
      end if
    end if

  end subroutine writeReal2


  !> Returns the format string for an entire row.
  subroutine getRowFormat(formStr, nItems, formStrRow)
    character(*), intent(in) :: formStr
    integer, intent(in) :: nItems
    character(:), allocatable, intent(out) :: formStrRow

    character(200) :: dummy
    integer :: nn

    write(dummy, '(I0)') nItems
    nn = len_trim(formStr) + 2 + len_trim(dummy)
    allocate(character(nn) :: formStrRow)
    write(formStrRow, '(A,I0,A,A)') '(', nItems, trim(formStr), ')'

  end subroutine getRowFormat

end module dftbp_logger
