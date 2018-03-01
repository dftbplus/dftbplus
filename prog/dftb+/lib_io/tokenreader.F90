!--------------------------------------------------------------------------------------------------!
!  DFTB+: general package for performing fast atomistic simulations                                !
!  Copyright (C) 2018  DFTB+ developers group                                                      !
!                                                                                                  !
!  See the LICENSE file for terms of usage and distribution.                                       !
!--------------------------------------------------------------------------------------------------!

#:include 'common.fypp'

!> Contains the low level utilities for parsing xml-data into intrinsic Fortran types.
!>
!> This module contains the utilities which can parse a strings into Fortran intrinsic types. Tokens
!> are assumed to be separated by white space, therefore strings with spaces inside can not
!> currently be handled.
module tokenreader
  use assert
  use charmanip
  use message, only : error
  use accuracy, only : dp
  use xmlf90
  implicit none

  private


  !> Flag for signals successfull token reading
  integer, parameter :: TOKEN_OK = 0


  !> Flag for signals end of string
  integer, parameter :: TOKEN_EOS = -1


  !> Flag for signals reading error
  integer, parameter :: TOKEN_ERROR = -2


  !> Contains procedures which read the next token from a string
  interface getNextToken
    module procedure getNextToken_string
    module procedure getNextToken_integer
    module procedure getNextToken_integerR1
    module procedure getNextToken_real
    module procedure getNextToken_realR1
    module procedure getNextToken_logical
    module procedure getNextToken_logicalR1
  end interface getNextToken


  !> Character representation of the logical true value
  character(len=*), parameter :: LOGICAL_TRUE = "Yes"


  !> Lower cased character representation of the logical true value
  character(len=*), parameter :: LOGICAL_TRUE_LO = "yes"


  !> Character representation of the logical false value
  character(len=*), parameter :: LOGICAL_FALSE = "No"


  !> Lower cased character representation of the logical false value
  character(len=*), parameter :: LOGICAL_FALSE_LO = "no"

  public :: getNextToken, TOKEN_OK, TOKEN_EOS, TOKEN_ERROR
  public :: LOGICAL_TRUE, LOGICAL_FALSE, LOGICAL_TRUE_LO, LOGICAL_FALSE_LO

contains


  !> Returns the next token from the provided string as integer
  subroutine getNextToken_integer(str, tokenValue, start, iostat)

    !> String to parse
    character(len=*), intent(in) :: str

    !> Contains the value of the token on return
    integer, intent(out) :: tokenValue

    !> Starting position for the parsing on call, first position after the end of the token on
    !> return.
    integer, intent(inout) :: start

    !> Token reader i/o status flag on return
    integer, intent(out), optional :: iostat

    integer :: iStart, iError, tokStart, tokLen, tokEnd

    tokenValue = 0
    iStart = start
    call getNextToken_local(str, tokStart, tokEnd, tokLen, iStart)
    if (tokLen == 0) then
      iError = TOKEN_EOS
    else if (.not. validIntegerStart(str(tokStart:tokStart))) then
      ! Workaround:PGI 17.10 -> strings starting with certain letters are interpreted as integers
      ! therefore we check explicitly
      iError = TOKEN_ERROR
    else
      read (str(tokStart:tokEnd), *, iostat=iError) tokenValue
      if (iError /= 0) then
        iError = TOKEN_ERROR
      else
        iError = TOKEN_OK
        start = iStart
      end if
    end if
    if (present(iostat)) then
      iostat = iError
    elseif (iError /= TOKEN_OK) then
      call error("Integer reading error")
    end if

  end subroutine getNextToken_integer


  !> Returns the next token from the provided string as rank one integer array
  subroutine getNextToken_integerR1(str, tokenValue, start, iostat, nItem)

    !> String to parse
    character(len=*), intent(in) :: str

    !> Contains the value of the token on return
    integer, intent(out) :: tokenValue(:)

    !> Starting position for the parsing on call, first position after the end of the token on
    !> return.
    integer, intent(inout) :: start

    !> Token reader i/o status flag on return
    integer, intent(out), optional :: iostat

    !> Nr. of read items
    integer, intent(out), optional :: nItem

    integer :: iStart, iError, nReadItem
    integer :: ii
    integer :: tmp

    tokenValue(:) = 0
    iError = TOKEN_OK
    iStart = start
    do ii = 1, size(tokenValue)
      call getNextToken(str, tmp, iStart, iError)
      if (iError /= TOKEN_OK) then
        exit
      end if
      tokenValue(ii) = tmp
    end do

    if (iError == TOKEN_OK) then
      start = iStart
      nReadItem = size(tokenValue)
    else
      nReadItem = ii - 1
    end if

    if (present(nItem)) then
      nItem = nReadItem
    end if

    if (present(iostat)) then
      iostat = iError
    elseif (iError /= TOKEN_OK) then
      call error("Integer reading error")
    end if

  end subroutine getNextToken_integerR1


  !> Returns the next token from the provided string as string
  subroutine getNextToken_string(str, tokenValue, start, iostat)

    !> String to parse
    character(len=*), intent(in) :: str

    !> Contains the value of the token on return
    type(string), intent(inout) :: tokenValue

    !> Starting position for the parsing on call, first position after the end of the token on
    !> return.
    integer, intent(inout) :: start

    !> Token reader i/o status flag on return
    integer, intent(out), optional :: iostat

    integer :: iError, tokStart, tokEnd, tokLen

    call getNextToken_local(str, tokStart, tokEnd, tokLen, start, &
        &ignoreQuotation=.false.)
    if (tokLen == 0) then
      iError = TOKEN_EOS
    else
      iError = TOKEN_OK
      tokenValue = str(tokStart:tokEnd)
    end if

    if (present(iostat)) then
      iostat = iError
    end if

  end subroutine getNextToken_string


  !> Returns the next token from the provided string as a real value.
  subroutine getNextToken_real(str, tokenValue, start, iostat)

    !> String to parse
    character(len=*), intent(in) :: str

    !> Contains the value of the token on return
    real(dp), intent(out) :: tokenValue

    !> Starting position for the parsing on call, first position after the end of the token on
    !> return.
    integer, intent(inout) :: start

    !> Token reader i/o status flag on return
    integer, intent(out), optional :: iostat

    integer :: iStart, iError, tokStart, tokEnd, tokLen

    tokenValue = 0.0_dp
    iStart = start
    call getNextToken_local(str, tokStart, tokEnd, tokLen, iStart)
    if (tokLen == 0) then
      iError = TOKEN_EOS
    else
      read (str(tokStart:tokEnd), *, iostat=iError) tokenValue
      if (iError /= 0) then
        iError = TOKEN_ERROR
      else
        iError = TOKEN_OK
        start = iStart
      end if
    end if
    if (present(iostat)) then
      iostat = iError
    else
      if (iError == TOKEN_ERROR) then
        call error("Real reading error")
      end if
    end if

  end subroutine getNextToken_real


  !> Returns the next token from the provided string as rank one real array
  subroutine getNextToken_realR1(str, tokenValue, start, iostat, nItem)

    !> String to parse
    character(len=*), intent(in) :: str

    !> Contains the value of the token on return
    real(dp), intent(out) :: tokenValue(:)

    !> Starting position for the parsing on call, first position after the end of the token on
    !> return.
    integer, intent(inout) :: start

    !> Token reader i/o status flag on return
    integer, intent(out), optional :: iostat

    !> Nr. of read items
    integer, intent(out), optional :: nItem

    integer :: iStart, iError, nReadItem
    integer :: ii
    real(dp) :: tmp

    tokenValue(:) = 0.0_dp
    iError = TOKEN_OK
    iStart = start
    do ii = 1, size(tokenValue)
      call getNextToken(str, tmp, iStart, iError)
      if (iError /= TOKEN_OK) then
        exit
      end if
      tokenValue(ii) = tmp
    end do

    if (iError == TOKEN_OK) then
      start = iStart
      nReadItem = size(tokenValue)
    else
      nReadItem = ii - 1
    end if

    if (present(nItem)) then
      nItem = nReadItem
    end if

    if (present(iostat)) then
      iostat = iError
    elseif (iError /= TOKEN_OK) then
      call error("Real reading error")
    end if

  end subroutine getNextToken_realR1


  !> Returns the next token from the provided string as logical
  subroutine getNextToken_logical(str, tokenValue, start, iostat)

    !> String to parse
    character(len=*), intent(in) :: str

    !> Contains the value of the token on return
    logical, intent(out) :: tokenValue

    !> Starting position for the parsing on call, first position after the end of the token on
    !> return.
    integer, intent(inout) :: start

    !> Token reader i/o status flag on return
    integer, intent(out), optional :: iostat

    integer :: iStart, iError, tokStart, tokEnd, tokLen
    character(len=len(str)) :: buffer

    tokenValue = .false.
    iStart = start
    iError = TOKEN_OK
    call getNextToken_local(str, tokStart, tokEnd, tokLen, iStart)
    if (tokLen == 0) then
      iError = TOKEN_EOS
    else
      buffer = tolower(str(tokStart:tokEnd))
      if (trim(buffer) == LOGICAL_TRUE_LO) then
        tokenValue = .true.
      elseif (trim(buffer) == LOGICAL_FALSE_LO) then
        tokenValue = .false.
      else
        iError = TOKEN_ERROR
      end if
    end if
    start = iStart

    if (present(iostat)) then
      iostat = iError
    else
      if (iError == TOKEN_ERROR) then
        call error("Token reading error")
      end if
    end if

  end subroutine getNextToken_logical


  !> Returns the next token from the provided string as logical
  subroutine getNextToken_logicalR1(str, tokenValue, start, iostat, nItem)

    !> String to parse
    character(len=*), intent(in) :: str

    !> Contains the value of the token on return
    logical, intent(out) :: tokenValue(:)

    !> Starting position for the parsing on call, first position after the end of the token on
    !> return.
    integer, intent(inout) :: start

    !> Token reader i/o status flag on return
    integer, intent(out), optional :: iostat

    !> Nr. of read items
    integer, intent(out), optional :: nItem

    integer :: iStart, iError, nReadItem
    integer :: ii
    logical :: tmp

    tokenValue = .false.
    iStart = start
    iError = TOKEN_OK
    do ii = 1, size(tokenValue)
      call getNextToken(str, tmp, iStart, iError)
      if (iError /= TOKEN_OK) then
        exit
      end if
      tokenValue(ii) = tmp
    end do

    if (iError == TOKEN_OK) then
      start = iStart
      nReadItem = size(tokenValue)
    else
      nReadItem = ii - 1
    end if

    if (present(nItem)) then
      nItem = nReadItem
    end if

    if (present(iostat)) then
      iostat = iError
    elseif (iError /= TOKEN_OK) then
      call error("Logical reading error")
    end if

  end subroutine getNextToken_logicalR1


  !> Returns the next token from the provided string
  !>
  !> If the string does not contain any tokens, the empty string is returned.
  subroutine getNextToken_local(str, tokStart, tokEnd, tokLen, start, ignoreQuotation)

    !> String to parse
    character(len=*), intent(in) :: str

    !> stating location of token in string
    integer, intent(out) :: tokStart

    !> end position of token in string
    integer, intent(out) :: tokEnd

    !> length of the token
    integer, intent(out) :: tokLen

    !> Starting position for the parsing on call, first position after the end of the token on
    !> return
    integer, intent(inout) :: start

    !> ignore quotation marks (Default: yes)
    logical, intent(in), optional :: ignoreQuotation

    integer :: lenStr
    logical :: tIgnoreQuotation

    @:ASSERT(start > 0)

    if (present(ignoreQuotation)) then
      tIgnoreQuotation = ignoreQuotation
    else
      tIgnoreQuotation = .true.
    end if

    lenStr = len(str)

    tokStart = complementaryScan(str(start:), whiteSpaces)
    if (tokStart == 0) then
      tokLen = 0
      tokEnd = 0
      return
    else
      tokStart = tokStart + start - 1
    end if
    if (tIgnoreQuotation) then
      tokEnd = scan(str(tokStart:), whiteSpaces)
    else
      tokEnd = unquotedScan(str(tokStart:), whiteSpaces)
    end if
    if (tokEnd == 0) then
      tokEnd = lenStr
    else
      tokEnd = tokEnd + tokStart - 2
    end if

    tokLen = tokEnd - tokStart + 1
    start = tokEnd + 2

  end subroutine getNextToken_local


  !> Cheks whether a given character represents a valid staring character for an integer.
  pure function validIntegerStart(char) result(tValid)

    !> Character to check
    character, intent(in) :: char

    !> Whether it can be a starting character of a string representing an integer.
    logical :: tValid

    integer :: ind

    ind = iachar(char)
    ! Either a digit 0-9 or - or +
    tValid = ((ind >= 48 .and. ind <= 57) .or. ind == 45 .or. ind == 43)

  end function validIntegerStart

end module tokenreader
