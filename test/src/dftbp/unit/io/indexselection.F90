!--------------------------------------------------------------------------------------------------!
!  DFTB+: general package for performing fast atomistic simulations                                !
!  Copyright (C) 2006 - 2025  DFTB+ developers group                                               !
!                                                                                                  !
!  See the LICENSE file for terms of usage and distribution.                                       !
!--------------------------------------------------------------------------------------------------!

#:include "fortuno_serial.fypp"

module test_io_indexselection
  use dftbp_common_status, only : TStatus
  use dftbp_io_indexselection, only : getIndexSelection, errors
  use fortuno_serial, only : suite => serial_suite_item, test_list
  $:FORTUNO_SERIAL_IMPORTS()
  implicit none

  private
  public :: tests


  ! Default initialized fixture for atomSelection tests
  type :: TAtomSelectionFx
    character(2) :: speciesNames(2) = [character(2) :: "Si", "C"]
    integer :: nAtom = 16
    integer :: species(16) = [1, 2, 1, 2, 1, 2, 1, 2, 1, 2, 1, 2, 1, 2, 1, 2]
    integer :: fullRange(2) = [1, 16]
    logical :: selected(16)
    type(TStatus) ::errStatus
  end type TAtomSelectionFx

  ! Default initialized fixture for numeric tests
  type :: TNumericFx
    integer :: fullRange(2) = [1, 10]
    integer :: nItem = 10
    logical :: selected(10)
    type(TStatus) :: errStatus
  end type TNumericFx


  ! Default initialized fixture for negativeNumeric tests
  type :: TNegativeNumericFx
    integer :: fullRange(2) = [-2, 2]
    integer :: nItem = 5
    logical :: selected(5)
    type(TStatus) :: errStatus
  end type TNegativeNumericFx

contains


  $:TEST("explicit", label="atomSelection")
    type(TAtomSelectionFx) :: fx
    call getIndexSelection("1 2 4", fx%fullRange, fx%selected, fx%errStatus,&
        & speciesNames=fx%speciesNames, species=fx%species)
    @:ASSERT(fx%errStatus%isOk())
    @:ASSERT(checkSelection(fx%selected, fx%fullRange, [1, 2, 4]))
  $:END_TEST()


  $:TEST("allAsRange", label="atomSelection")
    type(TAtomSelectionFx) :: fx
    call getIndexSelection("1:-1", fx%fullRange, fx%selected, fx%errStatus,&
        & speciesNames=fx%speciesNames, species=fx%species)
    @:ASSERT(fx%errStatus%isOk())
    @:ASSERT(all(fx%selected))
  $:END_TEST()


  $:TEST("range", label="atomSelection")
    type(TAtomSelectionFx) :: fx
    call getIndexSelection("2:6",fx%fullRange,fx%selected, fx%errStatus,&
        & speciesNames=fx%speciesNames, species=fx%species)
    @:ASSERT(fx%errStatus%isOk())
    @:ASSERT(checkSelection(fx%selected,fx%fullRange, [2, 3, 4, 5, 6]))
  $:END_TEST()


  $:TEST("negativeSubrange", label="atomSelection")
    type(TAtomSelectionFx) :: fx
    call getIndexSelection("-6:-2",fx%fullRange,fx%selected, fx%errStatus,&
        & speciesNames=fx%speciesNames, species=fx%species)
    @:ASSERT(fx%errStatus%isOk())
    @:ASSERT(checkSelection(fx%selected,fx%fullRange, [11, 12, 13, 14, 15]))
  $:END_TEST()


  $:TEST("speciesName", label="atomSelection")
    type(TAtomSelectionFx) :: fx
    call getIndexSelection("Si",fx%fullRange,fx%selected, fx%errStatus,&
        & speciesNames=fx%speciesNames, species=fx%species)
    @:ASSERT(fx%errStatus%isOk())
    @:ASSERT(checkSelection(fx%selected,fx%fullRange, [1, 3, 5, 7, 9, 11, 13, 15]))
  $:END_TEST()


  $:TEST("speciesNameDiffCase", label="atomSelection")
    type(TAtomSelectionFx) :: fx
    call getIndexSelection("si",fx%fullRange,fx%selected, fx%errStatus,&
        & speciesNames=fx%speciesNames, species=fx%species)
    @:ASSERT(fx%errStatus%isOk())
    @:ASSERT(checkSelection(fx%selected,fx%fullRange, [1, 3, 5, 7, 9, 11, 13, 15]))
  $:END_TEST()


  $:TEST("excludeIndices", label="atomSelection")
    type(TAtomSelectionFx) :: fx
    call getIndexSelection("!4:14",fx%fullRange,fx%selected, fx%errStatus,&
        & speciesNames=fx%speciesNames, species=fx%species)
    @:ASSERT(fx%errStatus%isOk())
    @:ASSERT(checkSelection(fx%selected,fx%fullRange, [1, 2, 3, 15, 16]))
  $:END_TEST()


  $:TEST("excludeSpecies", label="atomSelection")
    type(TAtomSelectionFx) :: fx
    call getIndexSelection("!C",fx%fullRange,fx%selected, fx%errStatus,&
        & speciesNames=fx%speciesNames, species=fx%species)
    @:ASSERT(fx%errStatus%isOk())
    @:ASSERT(checkSelection(fx%selected,fx%fullRange, [1, 3, 5, 7, 9, 11, 13, 15]))
  $:END_TEST()


  $:TEST("excludesWithAnd", label="atomSelection")
    type(TAtomSelectionFx) :: fx
    call getIndexSelection("!Si & !1:4 & !12:16",fx%fullRange,fx%selected, fx%errStatus,&
        & speciesNames=fx%speciesNames, species=fx%species)
    @:ASSERT(fx%errStatus%isOk())
    @:ASSERT(checkSelection(fx%selected,fx%fullRange, [6, 8, 10]))
  $:END_TEST()


  $:TEST("notWithParenthesis", label="atomSelection")
    type(TAtomSelectionFx) :: fx
    call getIndexSelection("!(Si 1:4 12:16)",fx%fullRange,fx%selected, fx%errStatus,&
        & speciesNames=fx%speciesNames, species=fx%species)
    @:ASSERT(fx%errStatus%isOk())
    @:ASSERT(checkSelection(fx%selected,fx%fullRange, [6, 8, 10]))
  $:END_TEST()


  $:TEST("andOrPrecedence", label="atomSelection")
    type(TAtomSelectionFx) :: fx
    call getIndexSelection("1:3 6:8 & Si 9:12",fx%fullRange,fx%selected, fx%errStatus,&
        & speciesNames=fx%speciesNames, species=fx%species)
    @:ASSERT(fx%errStatus%isOk())
    @:ASSERT(checkSelection(fx%selected,fx%fullRange, [1, 2, 3, 7, 9, 10, 11, 12]))
  $:END_TEST()


  $:TEST("parenthesisPrecedence", label="atomSelection")
    type(TAtomSelectionFx) :: fx
    call getIndexSelection("(1:3 6:8) & Si 9:12",fx%fullRange,fx%selected, fx%errStatus,&
        & speciesNames=fx%speciesNames, species=fx%species)
    @:ASSERT(fx%errStatus%isOk())
    @:ASSERT(checkSelection(fx%selected,fx%fullRange, [1, 3, 7, 9, 10, 11, 12]))
  $:END_TEST()


  $:TEST("parenthesisPrecedence2", label="atomSelection")
    type(TAtomSelectionFx) :: fx
    call getIndexSelection("1 (1:3 6:8) & !Si 9:12",fx%fullRange,fx%selected, fx%errStatus,&
        & speciesNames=fx%speciesNames, species=fx%species)
    @:ASSERT(fx%errStatus%isOk())
    @:ASSERT(checkSelection(fx%selected,fx%fullRange, [1, 2, 6, 8, 9, 10, 11, 12]))
  $:END_TEST()


  $:TEST("parenthesisPrecedence3", label="atomSelection")
    type(TAtomSelectionFx) :: fx
    call getIndexSelection("1 & (1:3 6:8) & !Si 9:12",fx%fullRange,fx%selected, fx%errStatus,&
        & speciesNames=fx%speciesNames, species=fx%species)
    @:ASSERT(fx%errStatus%isOk())
    @:ASSERT(checkSelection(fx%selected,fx%fullRange, [9, 10, 11, 12]))
  $:END_TEST()


  $:TEST("zeroIndexingPositive", label="atomSelection")
    type(TAtomSelectionFx) :: fx
    call getIndexSelection("0 2 3", [0, fx%nAtom - 1],fx%selected, fx%errStatus,&
        & speciesNames=fx%speciesNames, species=fx%species)
    @:ASSERT(fx%errStatus%isOk())
    @:ASSERT(checkSelection(fx%selected,fx%fullRange, [1, 3, 4]))
  $:END_TEST()


  $:TEST("zeroIndexingNegative", label="atomSelection")
    type(TAtomSelectionFx) :: fx
    call getIndexSelection("-16:-14", [0, fx%nAtom - 1],fx%selected, fx%errStatus,&
        & speciesNames=fx%speciesNames, species=fx%species)
    @:ASSERT(fx%errStatus%isOk())
    @:ASSERT(checkSelection(fx%selected,fx%fullRange, [1, 2, 3]))
  $:END_TEST()


  $:TEST("subrange", label="atomSelection")
    type(TAtomSelectionFx) :: fx
    call getIndexSelection("11:13", [9, 16],fx%selected(9:16), fx%errStatus,&
        & speciesNames=fx%speciesNames, species=fx%species(9:16), indexRange=fx%fullRange)
    @:ASSERT(fx%errStatus%isOk())
    @:ASSERT(checkSelection(fx%selected(9:16), [9, 16], [11, 12, 13]))
  $:END_TEST()

  ! Failing tests

  $:TEST("selectWrongSpeciesName", label="atomSelection")
    type(TAtomSelectionFx) :: fx
    call getIndexSelection("Mg",fx%fullRange,fx%selected, fx%errStatus,&
        & speciesNames=fx%speciesNames, species=fx%species)
    @:ASSERT(fx%errStatus%hasError())
    @:ASSERT(fx%errStatus%code == errors%invalidSelector)
  $:END_TEST()


  $:TEST("outOfBounds", label="atomSelection")
    type(TAtomSelectionFx) :: fx
    call getIndexSelection("17",fx%fullRange,fx%selected, fx%errStatus,&
        & speciesNames=fx%speciesNames, species=fx%species)
    @:ASSERT(fx%errStatus%hasError())
    @:ASSERT(fx%errStatus%code == errors%invalidSelector)
  $:END_TEST()


  $:TEST("outOfBoundsPlusZeroInd", label="atomSelection")
    type(TAtomSelectionFx) :: fx
    call getIndexSelection("16", [0, fx%nAtom - 1],fx%selected, fx%errStatus,&
        & speciesNames=fx%speciesNames, species=fx%species)
    @:ASSERT(fx%errStatus%hasError())
    @:ASSERT(fx%errStatus%code == errors%invalidSelector)
  $:END_TEST()


  $:TEST("outOfBoundsMinus", label="atomSelection")
    type(TAtomSelectionFx) :: fx
    call getIndexSelection("-17",fx%fullRange,fx%selected, fx%errStatus,&
        & speciesNames=fx%speciesNames, species=fx%species)
    @:ASSERT(fx%errStatus%hasError())
    @:ASSERT(fx%errStatus%code == errors%invalidSelector)
  $:END_TEST()


  $:TEST("outOfBoundsMinusZeroInd", label="atomSelection")
    type(TAtomSelectionFx) :: fx
    call getIndexSelection("-17",fx%fullRange,fx%selected, fx%errStatus,&
        & speciesNames=fx%speciesNames, species=fx%species)
    @:ASSERT(fx%errStatus%hasError())
    @:ASSERT(fx%errStatus%code == errors%invalidSelector)
  $:END_TEST()


  $:TEST("negativeIndexWithoutBackward", label="atomSelection")
    type(TAtomSelectionFx) :: fx
    call getIndexSelection("-6:-2",fx%fullRange,fx%selected, fx%errStatus, &
      & speciesNames=fx%speciesNames, species=fx%species, backwardIndexing=.false.)
    @:ASSERT(fx%errStatus%hasError())
    @:ASSERT(fx%errStatus%code == errors%invalidSelector)
  $:END_TEST()


  $:TEST("outOfSubrangeBounds", label="atomSelection")
    type(TAtomSelectionFx) :: fx
    call getIndexSelection("1:10", [9, 16],fx%selected(9:16), fx%errStatus,&
        & speciesNames=fx%speciesNames, species=fx%species(9:16), indexRange=fx%fullRange)
    @:ASSERT(fx%errStatus%hasError())
    @:ASSERT(fx%errStatus%code == errors%invalidSelector)
  $:END_TEST()


  $:TEST("doubleOperator", label="atomSelection")
    type(TAtomSelectionFx) :: fx
    call getIndexSelection("1 && 2",fx%fullRange,fx%selected, fx%errStatus,&
        & speciesNames=fx%speciesNames, species=fx%species)
    @:ASSERT(fx%errStatus%hasError())
    @:ASSERT(fx%errStatus%code == errors%syntaxError)
  $:END_TEST()


  $:TEST("unclosedParenthesis", label="atomSelection")
    type(TAtomSelectionFx) :: fx
    call getIndexSelection("(1 & 2",fx%fullRange,fx%selected, fx%errStatus,&
        & speciesNames=fx%speciesNames, species=fx%species)
    @:ASSERT(fx%errStatus%hasError())
    @:ASSERT(fx%errStatus%code == errors%syntaxError)
  $:END_TEST()


  $:TEST("unopenedParenthesis", label="atomSelection")
    type(TAtomSelectionFx) :: fx
    call getIndexSelection("1 & 2)",fx%fullRange,fx%selected, fx%errStatus,&
        & speciesNames=fx%speciesNames, species=fx%species)
    @:ASSERT(fx%errStatus%hasError())
    @:ASSERT(fx%errStatus%code == errors%syntaxError)
  $:END_TEST()


  $:TEST("explicit", label="numeric")
    type(TNumericFx) :: fx
    call getIndexSelection("5 6 7", fx%fullRange, fx%selected, fx%errStatus)
    @:ASSERT(fx%errStatus%isOk())
    @:ASSERT(checkSelection(fx%selected, fx%fullRange, [5, 6, 7]))
  $:END_TEST()


  $:TEST("range", label="numeric")
    type(TNumericFx) :: fx
    call getIndexSelection("-4:-2", fx%fullRange, fx%selected, fx%errStatus)
    @:ASSERT(fx%errStatus%isOk())
    @:ASSERT(checkSelection(fx%selected, fx%fullRange, [7, 8, 9]))
  $:END_TEST()


  $:TEST("explicit", label="negativeNumeric")
    type(TNegativeNumericFx) :: fx
    call getIndexSelection("-1 0 1", fx%fullRange, fx%selected, fx%errStatus)
    @:ASSERT(fx%errStatus%isOk())
    @:ASSERT(checkSelection(fx%selected, fx%fullRange, [-1, 0, 1]))
  $:END_TEST()


  $:TEST("range", label="negativeNumeric")
    type(TNegativeNumericFx) :: fx
    call getIndexSelection("-1:1", fx%fullRange, fx%selected, fx%errStatus)
    @:ASSERT(fx%errStatus%isOk())
    @:ASSERT(checkSelection(fx%selected, fx%fullRange, [-1, 0, 1]))
  $:END_TEST()


  function tests()
    type(test_list) :: tests

    tests = test_list([&
        suite("indexselection", test_list([&
            $:TEST_ITEMS(label="atomSelection", suffix=",")
            $:TEST_ITEMS(label="numeric", suffix=",")
            $:TEST_ITEMS(label="negativeNumeric")
        ]))&
    ])
    $:STOP_ON_MISSING_TEST_ITEMS()

  end function tests


  !! Check whether the selected indices in a selection array match the expected ones
  function checkSelection(selected,fullRange, expectedIndices) result(matches)
    logical, intent(in) ::selected(:)
    integer, intent(in) ::fullRange(:)
    integer, intent(in) :: expectedIndices(:)
    logical :: matches

    integer, allocatable :: indices(:)
    integer :: ii

    matches = .false.
    indices = pack([(ii, ii =fullRange(1),fullRange(2))],selected)
    if (size(indices) /= size(expectedIndices)) return
    matches = all(indices == expectedIndices)

  end function checkSelection

end module test_io_indexselection
