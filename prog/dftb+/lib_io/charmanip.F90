!--------------------------------------------------------------------------------------------------!
!  DFTB+: general package for performing fast atomistic simulations                                !
!  Copyright (C) 2018  DFTB+ developers group                                                      !
!                                                                                                  !
!  See the LICENSE file for terms of usage and distribution.                                       !
!--------------------------------------------------------------------------------------------------!

#:include 'common.fypp'

!> Contains character manipulation routines
module charmanip
  use assert
  implicit none

  private

  ! Quotation related quantities

  !> count of these
  integer, parameter :: nQuoteChar = 2

  !> types of quote mark
  character(len=1), parameter :: quoteChars(nQuoteChar) = (/ "'", '"' /)

  ! Whitespace like characters

  !> blank space
  character(len=1), parameter :: space = " "

  !> line feed character
  character(len=1), parameter :: lineFeed = achar(10)

  !> carriage return
  character(len=1), parameter :: carriageReturn = achar(13)

  !> tab symbol
  character(len=1), parameter :: tabulator = achar(9)


  !> collected whitespace characters
  character(len=*), parameter :: whiteSpaces = space // lineFeed &
      &// carriageReturn // tabulator


  !> Whitespace characters not treated as such by fortran
  character(len=*), parameter :: unhandledWhiteSpaces = lineFeed // tabulator


  !> Newline character
  character(len=*), parameter :: newline = lineFeed


  !> Maximal character length for integers (including sign)
  integer, parameter :: maxIntLen = range(1) + 2

  public :: unquotedIndex, unquote, trim2, len_trim2, tolower, i2c
  public :: getNextQuotationPos, getFirstOccurance, complementaryScan
  public :: unquotedScan
  public :: space, lineFeed, carriageReturn, tabulator, whiteSpaces, newline
  public :: convertWhitespaces

contains


  !> Returns the first unquoted occurance of a substring in a string
  function unquotedIndex(string, substring) result(unqIndex)

    !> String to investigate (hay)
    character(len=*), intent(in) :: string

    !> Substring to look for (needle)
    character(len=*), intent(in) :: substring

    !> Position of the first occurance or zero if substring was not found
    integer :: unqIndex

    integer :: strPos, quoteStart, quoteEnd
    logical :: tFinished
    integer :: shift, lenStr

    quoteStart = 0
    quoteEnd = 0
    strPos = index(string, substring)
    tFinished = .false.
    lenStr = len_trim(string)

    do while (.not. tFinished)


      ! Occurance of substr after last quotation end -> get next quotation
      if (strPos > quoteEnd) then
        shift = quoteEnd
        call getNextQuotationPos(string(shift+1:lenStr), quoteStart, quoteEnd)
        quoteStart = quoteStart + shift
        quoteEnd = quoteEnd + shift


        ! Substring occurs after quotation start but before quotation end -> Look for next
        ! occurance. (If not found, no unquoted occurance exists.)
      elseif (strPos > quoteStart) then
        shift = strPos
        strPos = index(string(strPos+1:lenStr), substring)
        if (strPos == 0) then
          tFinished = .true.
        else
          strPos = strPos + shift
        end if


        ! Substring before quotation start -> unquoted occurance found
      else
        tFinished = .true.
      end if
    end do

    unqIndex = strPos

  end function unquotedIndex


  !> Unquotes a string by removing the paired quotation marks
  function unquote(string, optLower) result(unquoted)

    !> String to remove the quotation marks from
    character(len=*), intent(in) :: string

    !> Should unquoted part of string be converted to lowercase?
    logical, intent(in), optional :: optLower

    !> Unquoted string
    character(len=len(string)) :: unquoted

    integer :: quoteStart, quoteEnd, shift
    integer :: tmp
    integer :: lastPos
    integer :: lenStr
    logical :: tFinished, lower

    lower = .false.
    if (present(optLower)) then
      lower = optLower
    end if

    unquoted = repeat(" ", len(unquoted))
    quoteStart = 0
    quoteEnd = 0
    lastPos = 1
    tFinished = .false.
    lenStr = len(string)

    do while (.not. tFinished)
      shift = quoteEnd
      call getNextQuotationPos(string(shift+1:), quoteStart, quoteEnd)
      quoteStart = quoteStart + shift
      quoteEnd = quoteEnd + shift

      tmp = lastPos + quoteStart - shift - 2
      if (tmp >= lastPos) then
        unquoted(lastPos:tmp) = string(shift+1:quoteStart-1)
        lastPos = tmp + 1
      end if
      if (quoteStart < lenStr) then
        tmp = lastPos + quoteEnd - quoteStart - 2
        if (tmp >= lastPos) then
          if (lower) then
            unquoted(lastPos:tmp) = tolower(string(quoteStart+1:quoteEnd-1))
          else
            unquoted(lastPos:tmp) = string(quoteStart+1:quoteEnd-1)
          end if
          lastPos = tmp + 1
        end if
      else
        tFinished = .true.
      end if
    end do

  end function unquote


  !> Returns the starting and ending position of the next quotation
  !> Note: Starting and ending positions are retuned as integer greater than the string length if
  !> there is not a match for the quotes string
  subroutine getNextQuotationPos(str, qStart, qEnd)

    !> String to investigate
    character(len=*), intent(in) :: str

    !> Starting position of the quotation on return
    integer, intent(out) :: qStart

    !> Ending position of the quotation on return
    integer, intent(out) :: qEnd

    integer :: iType, tmp, lenStr
    integer :: ii

    iType = 0
    lenStr = len(str)
    qStart = lenStr + 1
    do ii = 1, nQuoteChar
      tmp = index(str, quoteChars(ii))
      if (tmp /= 0 .and. tmp < qStart) then
        qStart = tmp
        iType = ii
      end if
    end do

    ! If quotation start was found, look for the appropriate quotation end
    if (qStart < lenStr) then
      qEnd = index(str(qStart+1:), quoteChars(iType))
      if (qEnd == 0) then
        qEnd = lenStr + 1
      else
        qEnd = qStart + qEnd
      end if
    else
      qEnd = qStart + 1
    end if

  end subroutine getNextQuotationPos


  !> Returns the first occurance of any of the passed substrings in a string
  subroutine getFirstOccurance(string, substrs, masks, iSubstr, pos)

    !> String to investigate
    character(len=*), intent(in) :: string

    !> substrings to search for
    character(len=*), intent(in) :: substrs(:)

    !> mask for supplied substrings
    logical, intent(in) :: masks(:)

    !> which substring was found
    integer, intent(out) :: iSubstr

    !> its position in the string
    integer, intent(out) :: pos

    integer :: ii, iTmp, nSubstr

    nSubstr = size(substrs)

    @:ASSERT(size(masks) == nSubstr)


    !> Get first occurance of a separator
    iSubstr = 0
    pos = len(string) + 1
    do ii = 1, nSubstr
      if (masks(ii)) then
        iTmp = unquotedIndex(string, trim(substrs(ii)))
        if (iTmp /= 0 .and. iTmp < pos) then
          pos = iTmp
          iSubstr = ii
        end if
      end if
    end do

  end subroutine getFirstOccurance


  !> Scans a string for the first character not part of a supplied set.
  pure function complementaryScan(string, set, back) result(ind)

    !> String to investigate
    character(len=*), intent(in) :: string

    !> Set containing the non-interesting characters
    character(len=*), intent(in) :: set

    !> If search should be made backwards
    logical, intent(in), optional :: back

    !> Index of the first character not contained in the set or zero if not found.
    integer :: ind

    character(len=1) :: cc
    integer :: iStart, iEnd, iStep, lenSet
    logical :: tFound
    integer :: ii, jj

    lenSet = len(set)

    if (present(back)) then
      if (back) then
        iStart = len(string)
        iEnd = 1
        iStep = -1
      end if
    else
      iStart = 1
      iEnd = len(string)
      iStep = 1
    end if

    ind = 0
    do ii = iStart, iEnd, iStep
      cc = string(ii:ii)
      tFound = .false.
      do jj = 1, lenSet
        tFound = tFound .or. (cc == set(jj:jj))
      end do
      if (.not. tFound) then
        ind = ii
        exit
      end if
    end do

  end function complementaryScan


  !> Returns the first unquoted occurance of a substring in a string
  function unquotedScan(string, set) result(unqIndex)

    !> String to investigate (hay)
    character(len=*), intent(in) :: string

    !> Substring to look for (needle)
    character(len=*), intent(in) :: set

    !> Position of the first occurance or zero if substring was not found
    integer :: unqIndex

    integer :: strPos, quoteStart, quoteEnd
    logical :: tFinished
    integer :: shift

    quoteStart = 0
    quoteEnd = 0
    strPos = scan(string, set)
    tFinished = .false.

    do while (.not. tFinished)


      !> Occurance after last quotation end -> get next quotation
      if (strPos > quoteEnd) then
        shift = quoteEnd
        call getNextQuotationPos(string(shift+1:), quoteStart, quoteEnd)
        quoteStart = quoteStart + shift
        quoteEnd = quoteEnd + shift

        ! Char occurs after quotation start but before quotation end ->
        ! Look for next occurance. (If not found, no unquoted occurance exists.)
      elseif (strPos > quoteStart) then
        shift = strPos
        strPos = scan(string(shift+1:), set)
        if (strPos == 0) then
          tFinished = .true.
        else
          strPos = strPos + shift
        end if

        ! Substring before quotation start -> unquoted occurance found
      else
        tFinished = .true.
      end if
    end do

    unqIndex = strPos

  end function unquotedScan


  !> Length of a trimmed string if CR, LF and TAB count as trimmed characters.
  pure function len_trim2(string)

    !> String to investigate
    character(len=*), intent(in) :: string

    !> Length of the string
    integer :: len_trim2

    len_trim2 = complementaryScan(string, whiteSpaces, back=.true.)

  end function len_trim2


  !> Returns a trimmed string if CR, LF and TAB count as trimmed characters.
  function trim2(string)

    !> String to trim
    character(len=*), intent(in) :: string

    character(:), allocatable :: trim2

    trim2 = string(:len_trim2(string))

  end function trim2


  !> Returns a lowercase string
  pure function tolower(str) result(lower)

    !> String to convert to lowercase
    character(len=*), intent(in) :: str

    !> Lowercase string
    character(len=len(str)) :: lower

    integer :: ii, iTmp

    do ii = 1, len(str)
      iTmp = iachar(str(ii:ii))
      if (65 <= iTmp .and. iTmp <= 90) then
        lower(ii:ii) = achar(iTmp + 32)
      else
        lower(ii:ii) = str(ii:ii)
      end if
    end do

  end function tolower



  !> Converts an integer to a character string
  pure function i2c(number)

    !> Integer to convert
    integer, intent(in) :: number

    !> converted string
    character(:), allocatable :: i2c

    character(len=maxIntLen) :: buffer

    write(buffer, "(I0)") number
    i2c = trim(buffer)

  end function i2c


  !> Replaces whitespace characters not recognised by Fortran as such by spaces.
  subroutine convertWhitespaces(str)

    !> String to process.
    character(len=*), intent(inout) :: str

    integer :: ii

    ii = scan(str, unhandledWhiteSpaces)
    do while (ii > 0)
      str(ii:ii) = space
      ii = scan(str, unhandledWhiteSpaces)
    end do

  end subroutine convertWhitespaces

end module charmanip
