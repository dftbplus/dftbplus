!--------------------------------------------------------------------------------------------------!
!  DFTB+: general package for performing fast atomistic simulations                                !
!  Copyright (C) 2006 - 2020  DFTB+ developers group                                               !
!                                                                                                  !
!  See the LICENSE file for terms of usage and distribution.                                       !
!--------------------------------------------------------------------------------------------------!

#:include 'common.fypp'
#:include 'error.fypp'

!> Implements simple index selection with grammar containing trivial logical operators.
!>
!> Grammar inspired by https://gist.github.com/mlabbe/81d667bac36aa60787fee60e3647a0a8/
!>
module dftbp_io_indexselection
  use dftbp_assert
  use dftbp_accuracy, only : dp
  use dftbp_charmanip, only : tolower
  use dftbp_status, only : TStatus
  implicit none

  private
  public :: getIndexSelection, errors


  !> Enumeration of possible token types
  type :: TTokenTypeEnum_

    !> No token
    integer :: empty = 0

    !> Atom selector (index, index range and eventually species name)
    integer :: selector = 1

    ! Logical NOT operator is implicit without token representation

    !> Logical AND operator
    integer :: and = 3

    !> Logical NOT operator
    integer :: not = 4

    !> Opening paranthesis
    integer :: open = 5

    !> Closing paranthesis
    integer :: close = 6

  end type  TTokenTypeEnum_

  !> Available token types
  type(TTokenTypeEnum_), parameter :: tokenTypes_ = TTokenTypeEnum_()


  !> Container for token information
  type :: TToken_

    !> Character content of the token
    character(:), allocatable :: content

    !> Type of the content
    integer :: type = tokenTypes_%empty

  end type TToken_


  !> Tokenizer returning tokens of a selection expression
  type :: TTokenizer_
    private
    character(:), allocatable :: expression
    integer :: pos
  contains
    procedure :: getNextToken => TTokenizer_getNextToken
  end type TTokenizer_


  !> Permanent infos needed during processing the selection expression
  type :: TSelection_

    !> Start and end index of the selectable range.
    integer :: selectionRange(2)

    !> Start and end index of the index range
    integer :: indexRange(2)

    !> Names of the species, in case it is an atom index selection. Shape: [nSpecies]
    character(:), allocatable :: speciesNames(:)

    !> Species of each atom, in case it is an atom index selection. Shape: [endInd - startInd + 1]
    integer, allocatable :: species(:)

    !> Whether backward indexing via negative index values is allowed
    logical :: backwardIndexing

  end type TSelection_


  !> Error types
  type :: TErrorEnum

    !> Syntax error in the selection expression
    integer :: syntaxError = 1

    !> Invalid value of an atom selector (wrong index value or non-existing species)
    integer :: invalidSelector = 2

  end type TErrorEnum

  !> Container for the errors this module can generate
  type(TErrorEnum), parameter :: errors = TErrorEnum()


contains

  !> Selects indices based on a selection expression.
  !>
  !> The selection expression must statisfy following grammar:
  !>
  !> expr := addTerm { addTerm }
  !> addTerm := mulTerm { "&" mulTerm }
  !> mulTerm := term |  "!"term
  !> term := selector | "(" expr ")"
  !> selector := index{":"index} | speciesName
  !> index := {-}[0-9]+
  !> speciesName := [a-zA-Z][0-9a-zA-Z_]*
  !>
  !> If no species names are provided, only numerical indices can be selected.
  !>
  subroutine getIndexSelection(expression, selectionRange, selected, status, indexRange,&
      & backwardIndexing, speciesNames, species)

    !> Selection expression
    character(*), intent(in) :: expression

    !> The range of indices [from, to] offered for selection.
    integer, intent(in) :: selectionRange(:)

    !> Logical indicating for each index, whether it is part of the selection or not. First element
    !> of the array corresponds to index selectionRange(1).
    logical, intent(out) :: selected(:)

    !> Success status with following error codes:
    !>
    !>  errors%syntaxError: Syntax error in the expression.
    !>  errors%invalidSelector: Invalid selector (e.g. index out of bonds or invalid species name)
    !>
    type(TStatus), intent(out) :: status

    !> The range of indices [from, to] available in general. It should at least contain the range
    !> specified in selectionRange. Default: the same range as in selectionRange.
    integer, optional, intent(in) :: indexRange(:)

    !> Whether backward indexing through negative index values is allowed. Applies only when both,
    !> the start and the end index of indexRange are positive. Default: .true.
    logical, optional, intent(in) :: backwardIndexing

    !> Names of atom species. Shape: [nSpecies].
    character(*), optional, intent(in) :: speciesNames(:)

    !> Species of each atom. Shape: [nAtom].
    integer, optional, intent(in) :: species(:)

    type(TSelection_) :: selection
    integer :: nSelectable
    type(TTokenizer_) :: tokenizer
    type(TToken_) :: token

    @:ASSERT(size(selectionRange) == 2)
    @:ASSERT(selectionRange(2) >= selectionRange(1))
    @:ASSERT(present(speciesNames) .eqv. present(species))

    selection%selectionRange(:) = selectionRange
    nSelectable = selection%selectionRange(2) - selection%selectionRange(1) + 1
    @:ASSERT(size(selected) == nSelectable)
    if (present(indexRange)) then
      @:ASSERT(size(indexRange) == 2)
      @:ASSERT(indexRange(1) <= selectionRange(1) .and. indexRange(2) >= selectionRange(2))
      selection%indexRange(:) = indexRange
    else
      selection%indexRange(:) = selection%selectionRange
    end if
    if (present(speciesNames)) then
      @:ASSERT(size(species) == nSelectable)
      selection%species = species
      selection%speciesNames = tolower(speciesNames)
    end if

    selection%backwardIndexing = all(selection%indexRange >= 0)
    if (selection%backwardIndexing .and. present(backwardIndexing)) then
      selection%backwardIndexing = backwardIndexing
    end if

    call TTokenizer_init(tokenizer, expression)
    call tokenizer%getNextToken(token)
    selected(:) = .false.
    call evalExpression_(selection, 1, tokenizer, token, selected, status)
    @:PROPAGATE_ERROR(status)

  end subroutine getIndexSelection


  !> Evaluates a selection expression.
  subroutine evalExpression_(selection, level, tokenizer, token, selected, status)
    type(TSelection_), intent(in) :: selection
    integer, intent(in) :: level
    type(TTokenizer_), intent(inout) :: tokenizer
    type(TToken_), intent(inout) :: token
    logical, intent(out) :: selected(:)
    type(TStatus), intent(inout) :: status

    logical, allocatable :: buffer(:)

    call evalAdditiveTerm_(selection, level, tokenizer, token, selected, status)
    @:PROPAGATE_ERROR(status)
    ! Apply implicit OR operator if next token is of appropriate type
    if (any(token%type == [tokenTypes_%not, tokenTypes_%selector, tokenTypes_%open])) then
      allocate(buffer(size(selected)))
      do while (any(token%type == [tokenTypes_%not, tokenTypes_%selector, tokenTypes_%open]))
        call evalAdditiveTerm_(selection, level, tokenizer, token, buffer, status)
        selected(:) = selected .or. buffer
        @:PROPAGATE_ERROR(status)
      end do
    end if
    if (level == 1 .and. token%type /= tokenTypes_%empty) then
      @:RAISE_FORMATTED_ERROR(status, errors%syntaxError,&
          & "('Unexpected token ''', A, ''' found')", token%content)
    end if

  end subroutine evalExpression_


  !> Evaluetes an additive term.
  subroutine evalAdditiveTerm_(selection, level, tokenizer, token, selected, status)
    type(TSelection_), intent(in) :: selection
    integer, intent(in) :: level
    type(TTokenizer_), intent(inout) :: tokenizer
    type(TToken_), intent(inout) :: token
    logical, intent(out) :: selected(:)
    type(TStatus), intent(inout) :: status

    logical, allocatable :: buffer(:)

    call evalMultiplicativeTerm_(selection, level, tokenizer, token, selected, status)
    @:PROPAGATE_ERROR(status)
    if (token%type == tokenTypes_%and) then
      allocate(buffer(size(selected)))
      do while (token%type == tokenTypes_%and)
        call tokenizer%getNextToken(token)
        call evalMultiplicativeTerm_(selection, level, tokenizer, token, buffer, status)
        selected(:) = selected .and. buffer
        @:PROPAGATE_ERROR(status)
      end do
    end if

  end subroutine evalAdditiveTerm_


  !> Evaluates a multiplicative term.
  subroutine evalMultiplicativeTerm_(selection, level, tokenizer, token, selected, status)
    type(TSelection_), intent(in) :: selection
    integer, intent(in) :: level
    type(TTokenizer_), intent(inout) :: tokenizer
    type(TToken_), intent(inout) :: token
    logical, intent(out) :: selected(:)
    type(TStatus), intent(inout) :: status

    logical, allocatable :: buffer(:)

    if (token%type == tokenTypes_%not) then
      call tokenizer%getNextToken(token)
      allocate(buffer(size(selected)))
      call evalTerm_(selection, level, tokenizer, token, buffer, status)
      @:PROPAGATE_ERROR(status)
      selected(:) = .not. buffer
    else
      call evalTerm_(selection, level, tokenizer, token, selected, status)
      @:PROPAGATE_ERROR(status)
    end if

  end subroutine evalMultiplicativeTerm_


  !> Evaluates a simple term.
  subroutine evalTerm_(selection, level, tokenizer, token, selected, status)
    type(TSelection_), intent(in) :: selection
    integer, intent(in) :: level
    type(TTokenizer_), intent(inout) :: tokenizer
    type(TToken_), intent(inout) :: token
    logical, intent(out) :: selected(:)
    type(TStatus), intent(inout) :: status

    if (token%type == tokenTypes_%empty) then
      selected(:) = .false.
      return
    else if (token%type == tokenTypes_%selector) then
      call evalSelector_(token%content, selection, selected, status)
      @:PROPAGATE_ERROR(status)
    else if (token%type == tokenTypes_%open) then
      call tokenizer%getNextToken(token)
      call evalExpression_(selection, level + 1, tokenizer, token, selected, status)
      @:PROPAGATE_ERROR(status)
      if (token%type /= tokenTypes_%close) then
        @:RAISE_ERROR(status, errors%syntaxError, "Missing closing paranthesis ')'")
      end if
    else
      @:RAISE_FORMATTED_ERROR(status, errors%syntaxError,&
          & "('Unexpected token ''', A, ''' found')", token%content)
    end if
    call tokenizer%getNextToken(token)

  end subroutine evalTerm_


  !> Evaluates an atom selector.
  subroutine evalSelector_(selector, selection, selected, status)
    character(*), intent(in) :: selector
    type(TSelection_), intent(in) :: selection
    logical, intent(inout) :: selected(:)
    type(TStatus), intent(inout) :: status

    integer :: seppos, iFirst, iLast, iFirstNorm, iLastNorm, iSpecies, iFirstChar
    integer :: iostat
    character :: firstChar

    firstChar = selector(1:1)
    iFirstChar = ichar(firstChar)

    selected(:) = .false.

    ! Numerical index value
    if (iFirstChar  >= ichar('0') .and. iFirstChar <= ichar('9') .or. firstChar == '-') then
      seppos = index(selector, ':')
      if (seppos > 0) then
        read(selector(: seppos - 1), *, iostat=iostat) iFirst
        if (iostat /= 0) then
          @:RAISE_FORMATTED_ERROR(status, errors%syntaxError,&
              & "('Non-integer index ''', A,  '''')", selector(: seppos - 1))
        end if
        read(selector(seppos + 1 :), *, iostat=iostat) iLast
        if (iostat /= 0) then
          @:RAISE_FORMATTED_ERROR(status, errors%syntaxError,&
              & "('Non-integer index ''', A,  '''')", selector(seppos + 1 :))
        end if
      else
        read(selector, *, iostat=iostat) iFirst
        if (iostat /= 0) then
          @:RAISE_FORMATTED_ERROR(status, errors%syntaxError,&
              & "('Non-integer index ''', A,  '''')", selector)
        end if
        iLast = iFirst
      end if
      iFirstNorm = normalizedIndex_(iFirst, selection%selectionRange, selection%indexRange,&
          & selection%backwardIndexing)
      if (iFirstNorm == 0) then
        @:RAISE_FORMATTED_ERROR(status, errors%invalidSelector,&
            & "('Out of bounds index value ', I0)", iFirst)
      end if
      iLastNorm = normalizedIndex_(iLast, selection%selectionRange, selection%indexRange,&
          & selection%backwardIndexing)
      if (iLastNorm == 0) then
        @:RAISE_FORMATTED_ERROR(status, errors%invalidSelector,&
            & "('Out of bounds index value ', I0)", iLast)
      end if
      selected(iFirstNorm : iLastNorm) = .true.

    ! Non-numerical value. If this is an atom selection, try to interprete it as species name.
    else if (allocated(selection%speciesNames)) then
      iSpecies = findloc(selection%speciesNames, tolower(selector), dim=1)
      if (iSpecies == 0) then
        @:RAISE_FORMATTED_ERROR(status, errors%invalidSelector,&
            & "('Invalid species name ''', A, '''')", trim(selector))
      end if
      where (selection%species == iSpecies)
        selected = .true.
      end where

    ! Non-numerical value in index selection
    else
      @:RAISE_FORMATTED_ERROR(status, errors%syntaxError,&
          & "('Invalid index selector ''', A, '''')", trim(selector))
    end if

  end subroutine evalSelector_


  !> Returns normalized atom index or 0 if index is invalid.
  pure function normalizedIndex_(ind, selectionRange, indexRange, backwardIndexing) &
        & result(normInd)
    integer, intent(in) :: ind, selectionRange(:), indexRange(:)
    logical, intent(in) :: backwardIndexing
    integer :: normInd

    normInd = ind
    if (ind < 0 .and. backwardIndexing) then
      normInd = indexRange(2) + ind + 1
    end if
    if (normInd >= selectionRange(1) .and. normInd <= selectionRange(2)) then
      normInd = normInd - selectionRange(1) + 1
    else
      normInd = 0
    end if

  end function normalizedIndex_


  !> Initializes a tokenizer.
  subroutine TTokenizer_init(this, expression)
    type(TTokenizer_), intent(out) :: this
    character(*), intent(in) :: expression

    this%expression = expression
    this%pos = 1

  end subroutine TTokenizer_init


  !> Delivers the next token. If no more tokens are available, the token type will be 'empty'.
  subroutine TTokenizer_getNextToken(this, token)
    class(TTokenizer_), intent(inout) :: this
    type(TToken_), intent(out) :: token

    call getNextToken_(this%expression, this%pos, token)

  end subroutine TTokenizer_getNextToken


  !> Helper function for TTokenizer_getNextToken().
  subroutine getNextToken_(expression, pos, token)
    character(*), intent(in) :: expression
    integer, intent(inout) :: pos
    type(TToken_), intent(out) :: token

    integer :: startPos
    character :: curChar

    do while (pos <= len(expression))
      if (expression(pos:pos) /= " ") then
        exit
      end if
      pos = pos + 1
    end do

    if (pos > len(expression)) then
      return
    end if

    if (scan(expression(pos:pos), "!()|&") /= 0) then
      token%content = expression(pos:pos)
      pos = pos + 1
    else
      startPos = pos
      pos = pos + 1
      do while (pos <= len(expression))
        curChar = expression(pos:pos)
        if (scan(curChar, " )|&") /= 0) then
          exit
        end if
        pos = pos + 1
      end do
      token%content = expression(startPos : pos - 1)
    end if

    token%type = getTokenType_(token%content)

  end subroutine getNextToken_


  !> Determines the type of a token.
  function getTokenType_(token) result(tokenType)
    character(*), intent(in) :: token
    integer :: tokenType

    if (token == "&") then
      tokenType = tokenTypes_%and
    else if (token == "!") then
      tokenType = tokenTypes_%not
    else if (token == "(") then
      tokenType = tokenTypes_%open
    else if (token == ")") then
      tokenType = tokenTypes_%close
    else
      tokenType = tokenTypes_%selector
    end if

  end function getTokenType_


end module dftbp_io_indexselection
