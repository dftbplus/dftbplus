!--------------------------------------------------------------------------------------------------!
!  DFTB+: general package for performing fast atomistic simulations                                !
!  Copyright (C) 2017  DFTB+ developers group                                                      !
!                                                                                                  !
!  See the LICENSE file for terms of usage and distribution.                                       !
!--------------------------------------------------------------------------------------------------!

#:include 'common.fypp'

!!* Contains the low level utilities for parsing xml-data into intrinsic Fortran
!!* types.
!!* @desc This module contains the utilities which can parse a strings into
!!*   Fortran intrinsic types. Tokens are assumed being separated by space,
!!*   therefore strings with spaces inside can not be handled yet.
module tokenreader
  use assert
  use charmanip
  use message, only : error
  use accuracy, only : dp
  use xmlf90
  implicit none

  private

  !!* Flag for signalising successfull token reading
  integer, parameter :: TOKEN_OK = 0

  !!* Flag for signalising end of string
  integer, parameter :: TOKEN_EOS = -1

  !!* Flag for signalising reading error
  integer, parameter :: TOKEN_ERROR = -2


  !!* Contains procedures which read the next token from a string
  interface getNextToken
    module procedure getNextToken_string
    module procedure getNextToken_integer
    module procedure getNextToken_integerR1
    module procedure getNextToken_real
    module procedure getNextToken_realR1
    module procedure getNextToken_logical
    module procedure getNextToken_logicalR1
  end interface

  !!* Character representation of the logical true value
  character(len=*), parameter :: LOGICAL_TRUE = "Yes"

  !!* Lower cased character representation of the logical true value
  character(len=*), parameter :: LOGICAL_TRUE_LO = "yes"

  !!* Character representation of the logical false value
  character(len=*), parameter :: LOGICAL_FALSE = "No"

  !!* Lower cased character representation of the logical false value
  character(len=*), parameter :: LOGICAL_FALSE_LO = "no"


  public :: getNextToken, TOKEN_OK, TOKEN_EOS, TOKEN_ERROR
  public :: LOGICAL_TRUE, LOGICAL_FALSE, LOGICAL_TRUE_LO, LOGICAL_FALSE_LO


contains


  !!* Returns the next token from the provided string as integer
  !!* @param str    String to parse
  !!* @param value  Contains the value of the token on return
  !!* @param start  Starting position for the parsing on call, first position
  !!*   after the end of the token on return.
  !!* @param iostat Token reader i/o status flag on return
  subroutine getNextToken_integer(str, value, start, iostat)
    character(len=*), intent(in) :: str
    integer, intent(out) :: value
    integer, intent(inout) :: start
    integer, intent(out), optional :: iostat

    integer :: iStart, iError, tokStart, tokLen, tokEnd
    integer :: iTmp

    value = 0
    iStart = start
    call getNextToken_local(str, tokStart, tokEnd, tokLen, iStart)
    if (tokLen == 0) then
      iError = TOKEN_EOS
    else
      !read (str(tokStart:tokEnd), *, iostat=iError) value
      !! Ugly hack for intel, because it interprets F, T as integers...
      iTmp = tokEnd - tokStart + 1
      read (str(tokStart:tokEnd), "(I"// i2c(iTmp) // ")", iostat=iError) value
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



  !!* Returns the next token from the provided string as rank one integer array
  !!* @param str    String to parse
  !!* @param value  Contains the value of the token on return
  !!* @param start  Starting position for the parsing on call, first position
  !!*   after the end of the token on return.
  !!* @param iostat Token reader i/o status flag on return
  !!* @param nItem    Nr. of read items
  subroutine getNextToken_integerR1(str, value, start, iostat, nItem)
    character(len=*), intent(in) :: str
    integer, intent(out) :: value(:)
    integer, intent(inout) :: start
    integer, intent(out), optional :: iostat
    integer, intent(out), optional :: nItem

    integer :: iStart, iError, nReadItem
    integer :: ii
    integer :: tmp

    value(:) = 0
    iError = TOKEN_OK
    iStart = start
    do ii = 1, size(value)
      call getNextToken(str, tmp, iStart, iError)
      if (iError /= TOKEN_OK) then
        exit
      end if
      value(ii) = tmp
    end do

    if (iError == TOKEN_OK) then
      start = iStart
      nReadItem = size(value)
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



  !!* Returns the next token from the provided string as string
  !!* @param str    String to parse
  !!* @param value  Contains the value of the token on return
  !!* @param start  Starting position for the parsing on call, first position
  !!*   after the end of the token on return.
  !!* @param iostat Token reader i/o status flag on return
  subroutine getNextToken_string(str, value, start, iostat)
    character(len=*), intent(in) :: str
    type(string), intent(inout) :: value
    integer, intent(inout) :: start
    integer, intent(out), optional :: iostat

    integer :: iError, tokStart, tokEnd, tokLen

    call getNextToken_local(str, tokStart, tokEnd, tokLen, start, &
        &ignoreQuotation=.false.)
    if (tokLen == 0) then
      iError = TOKEN_EOS
    else
      iError = TOKEN_OK
      value = str(tokStart:tokEnd)
    end if

    if (present(iostat)) then
      iostat = iError
    end if

  end subroutine getNextToken_string



  !!* Returns the next token from the provided string as real.
  !!* @param str    String to parse
  !!* @param value  Contains the value of the token on return
  !!* @param start  Starting position for the parsing on call, first position
  !!*   after the end of the token on return.
  !!* @param iostat Token reader i/o status flag on return
  subroutine getNextToken_real(str, value, start, iostat)
    character(len=*), intent(in) :: str
    real(dp), intent(out) :: value
    integer, intent(inout) :: start
    integer, intent(out), optional :: iostat

    integer :: iStart, iError, tokStart, tokEnd, tokLen

    value = 0.0_dp
    iStart = start
    call getNextToken_local(str, tokStart, tokEnd, tokLen, iStart)
    if (tokLen == 0) then
      iError = TOKEN_EOS
    else
      read (str(tokStart:tokEnd), *, iostat=iError) value
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



  !!* Returns the next token from the provided string as rank one real array
  !!* @param str    String to parse
  !!* @param value  Contains the value of the token on return
  !!* @param start  Starting position for the parsing on call, first position
  !!*   after the end of the token on return.
  !!* @param iostat Token reader i/o status flag on return
  !!* @param nItem    Nr. of read items
  subroutine getNextToken_realR1(str, value, start, iostat, nItem)
    character(len=*), intent(in) :: str
    real(dp), intent(out) :: value(:)
    integer, intent(inout) :: start
    integer, intent(out), optional :: iostat
    integer, intent(out), optional :: nItem

    integer :: iStart, iError, nReadItem
    integer :: ii
    real(dp) :: tmp

    value(:) = 0.0_dp
    iError = TOKEN_OK
    iStart = start
    do ii = 1, size(value)
      call getNextToken(str, tmp, iStart, iError)
      if (iError /= TOKEN_OK) then
        exit
      end if
      value(ii) = tmp
    end do

    if (iError == TOKEN_OK) then
      start = iStart
      nReadItem = size(value)
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



  !!* Returns the next token from the provided string as logical
  !!* @param str    String to parse
  !!* @param value  Contains the value of the token on return
  !!* @param start  Starting position for the parsing on call, first position
  !!*   after the end of the token on return.
  !!* @param iostat Token reader i/o status flag on return
  subroutine getNextToken_logical(str, value, start, iostat)
    character(len=*), intent(in) :: str
    logical, intent(out) :: value
    integer, intent(inout) :: start
    integer, intent(out), optional :: iostat

    integer :: iStart, iError, tokStart, tokEnd, tokLen
    character(len=len(str)) :: buffer

    value = .false.
    iStart = start
    iError = TOKEN_OK
    call getNextToken_local(str, tokStart, tokEnd, tokLen, iStart)
    if (tokLen == 0) then
      iError = TOKEN_EOS
    else
      buffer = tolower(str(tokStart:tokEnd))
      if (trim(buffer) == LOGICAL_TRUE_LO) then
        value = .true.
      elseif (trim(buffer) == LOGICAL_FALSE_LO) then
        value = .false.
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


  !!* Returns the next token from the provided string as logical
  !!* @param str    String to parse
  !!* @param value  Contains the value of the token on return
  !!* @param start  Starting position for the parsing on call, first position
  !!*   after the end of the token on return.
  !!* @param iostat Token reader i/o status flag on return
  !!* @param nItem    Nr. of read items
  subroutine getNextToken_logicalR1(str, value, start, iostat, nItem)
    character(len=*), intent(in) :: str
    logical, intent(out) :: value(:)
    integer, intent(inout) :: start
    integer, intent(out), optional :: iostat
    integer, intent(out), optional :: nItem

    integer :: iStart, iError, nReadItem
    integer :: ii
    logical :: tmp

    value = .false.
    iStart = start
    iError = TOKEN_OK
    do ii = 1, size(value)
      call getNextToken(str, tmp, iStart, iError)
      if (iError /= TOKEN_OK) then
        exit
      end if
      value(ii) = tmp
    end do

    if (iError == TOKEN_OK) then
      start = iStart
      nReadItem = size(value)
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

  !!* Returns the next token from the provided string
  !!* @param str             String to parse
  !!* @param token           Contains the value of the token on return
  !!* @param start           Starting position for the parsing on call, first
  !!*   position after the end of the token on return
  !!* @param ignoreQuotation Should quotation be ignored? (Default: yes)
  !!* @note If the string does not contain any tokens, the empty string is
  !!*   returned.
  subroutine getNextToken_local(str, tokStart, tokEnd, tokLen, start, &
      &ignoreQuotation)
    character(len=*), intent(in) :: str
    integer, intent(out) :: tokStart
    integer, intent(out) :: tokEnd
    integer, intent(out) :: tokLen
    integer, intent(inout) :: start
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


end module tokenreader
