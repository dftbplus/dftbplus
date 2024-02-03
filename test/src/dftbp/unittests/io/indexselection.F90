!--------------------------------------------------------------------------------------------------!
!  DFTB+: general package for performing fast atomistic simulations                                !
!  Copyright (C) 2006 - 2024  DFTB+ developers group                                               !
!                                                                                                  !
!  See the LICENSE file for terms of usage and distribution.                                       !
!--------------------------------------------------------------------------------------------------!

#:include "fytest.fypp"

#:block TEST_SUITE("indexselection")
  use dftbp_io_indexselection, only : getIndexSelection, errors
  use dftbp_common_status, only : TStatus
  implicit none

#:contains

  #:block TEST_FIXTURE("atomSelection")

    character(2), parameter :: speciesNames(*) = [character(2) :: "Si", "C"]
    integer, parameter :: species(*) = [1, 2, 1, 2, 1, 2, 1, 2, 1, 2, 1, 2, 1, 2, 1, 2]
    integer, parameter :: nAtom = size(species)
    integer, parameter :: fullRange(2) = [1, nAtom]
    logical :: selected(nAtom)
    type(TStatus) :: errStatus

  #:contains

    #:block TEST("explicit")
      call getIndexSelection("1 2 4", fullRange, selected, errStatus, speciesNames=speciesNames,&
          & species=species)
      @:ASSERT(errStatus%isOk())
      @:ASSERT(checkSelection(selected, fullRange, [1, 2, 4]))
    #:endblock

    #:block TEST("allAsRange")
      call getIndexSelection("1:-1", fullRange, selected, errStatus, speciesNames=speciesNames,&
          & species=species)
      @:ASSERT(errStatus%isOk())
      @:ASSERT(all(selected))
    #:endblock

    #:block TEST("range")
      call getIndexSelection("2:6", fullRange, selected, errStatus, speciesNames=speciesNames,&
          & species=species)
      @:ASSERT(errStatus%isOk())
      @:ASSERT(checkSelection(selected, fullRange, [2, 3, 4, 5, 6]))
    #:endblock

    #:block TEST("negativeSubrange")
      call getIndexSelection("-6:-2", fullRange, selected, errStatus, speciesNames=speciesNames,&
          & species=species)
      @:ASSERT(errStatus%isOk())
      @:ASSERT(checkSelection(selected, fullRange, [11, 12, 13, 14, 15]))
    #:endblock

    #:block TEST("speciesName")
      call getIndexSelection("Si", fullRange, selected, errStatus, speciesNames=speciesNames,&
          & species=species)
      @:ASSERT(errStatus%isOk())
      @:ASSERT(checkSelection(selected, fullRange, [1, 3, 5, 7, 9, 11, 13, 15]))
    #:endblock

    #:block TEST("speciesNameDiffCase")
      call getIndexSelection("si", fullRange, selected, errStatus, speciesNames=speciesNames,&
          & species=species)
      @:ASSERT(errStatus%isOk())
      @:ASSERT(checkSelection(selected, fullRange, [1, 3, 5, 7, 9, 11, 13, 15]))
    #:endblock

    #:block TEST("excludeIndices")
      call getIndexSelection("!4:14", fullRange, selected, errStatus, speciesNames=speciesNames,&
          & species=species)
      @:ASSERT(errStatus%isOk())
      @:ASSERT(checkSelection(selected, fullRange, [1, 2, 3, 15, 16]))
    #:endblock

    #:block TEST("excludeSpecies")
      call getIndexSelection("!C", fullRange, selected, errStatus, speciesNames=speciesNames,&
          & species=species)
      @:ASSERT(errStatus%isOk())
      @:ASSERT(checkSelection(selected, fullRange, [1, 3, 5, 7, 9, 11, 13, 15]))
    #:endblock

    #:block TEST("excludesWithAnd")
      call getIndexSelection("!Si & !1:4 & !12:16", fullRange, selected, errStatus,&
           speciesNames=speciesNames, species=species)
      @:ASSERT(errStatus%isOk())
      @:ASSERT(checkSelection(selected, fullRange, [6, 8, 10]))
    #:endblock

    #:block TEST("notWithParenthesis")
      call getIndexSelection("!(Si 1:4 12:16)", fullRange, selected, errStatus,&
           speciesNames=speciesNames, species=species)
      @:ASSERT(errStatus%isOk())
      @:ASSERT(checkSelection(selected, fullRange, [6, 8, 10]))
    #:endblock

    #:block TEST("andOrPrecedence")
      call getIndexSelection("1:3 6:8 & Si 9:12", fullRange, selected, errStatus,&
           speciesNames=speciesNames, species=species)
      @:ASSERT(errStatus%isOk())
      @:ASSERT(checkSelection(selected, fullRange, [1, 2, 3, 7, 9, 10, 11, 12]))
    #:endblock

    #:block TEST("parenthesisPrecedence")
      call getIndexSelection("(1:3 6:8) & Si 9:12", fullRange, selected, errStatus,&
          speciesNames=speciesNames, species=species)
      @:ASSERT(errStatus%isOk())
      @:ASSERT(checkSelection(selected, fullRange, [1, 3, 7, 9, 10, 11, 12]))
    #:endblock

    #:block TEST("parenthesisPrecedence2")
      call getIndexSelection("1 (1:3 6:8) & !Si 9:12", fullRange, selected, errStatus,&
          speciesNames=speciesNames, species=species)
      @:ASSERT(errStatus%isOk())
      @:ASSERT(checkSelection(selected, fullRange, [1, 2, 6, 8, 9, 10, 11, 12]))
    #:endblock

    #:block TEST("parenthesisPrecedence3")
      call getIndexSelection("1 & (1:3 6:8) & !Si 9:12", fullRange, selected, errStatus,&
          speciesNames=speciesNames, species=species)
      @:ASSERT(errStatus%isOk())
      @:ASSERT(checkSelection(selected, fullRange, [9, 10, 11, 12]))
    #:endblock

    #:block TEST("zeroIndexingPositive")
      call getIndexSelection("0 2 3", [0, nAtom - 1], selected, errStatus,&
          speciesNames=speciesNames, species=species)
      @:ASSERT(errStatus%isOk())
      @:ASSERT(checkSelection(selected, fullRange, [1, 3, 4]))
    #:endblock

    #:block TEST("zeroIndexingNegative")
      call getIndexSelection("-16:-14", [0, nAtom - 1], selected, errStatus,&
          speciesNames=speciesNames, species=species)
      @:ASSERT(errStatus%isOk())
      @:ASSERT(checkSelection(selected, fullRange, [1, 2, 3]))
    #:endblock

    #:block TEST("subrange")
      call getIndexSelection("11:13", [9, 16], selected(9:16), errStatus,&
          & speciesNames=speciesNames, species=species(9:16), indexRange=fullRange)
      @:ASSERT(errStatus%isOk())
      @:ASSERT(checkSelection(selected(9:16), [9, 16], [11, 12, 13]))
    #:endblock

    #! Failing tests

    #:block TEST("selectWrongSpeciesName")
      call getIndexSelection("Mg", fullRange, selected, errStatus, speciesNames=speciesNames,&
          & species=species)
      @:ASSERT(errStatus%hasError())
      @:ASSERT(errStatus%code == errors%invalidSelector)
    #:endblock

    #:block TEST("outOfBounds")
      call getIndexSelection("17", fullRange, selected, errStatus, speciesNames=speciesNames,&
          & species=species)
      @:ASSERT(errStatus%hasError())
      @:ASSERT(errStatus%code == errors%invalidSelector)
    #:endblock

    #:block TEST("outOfBoundsPlusZeroInd")
      call getIndexSelection("16", [0, nAtom - 1], selected, errStatus, speciesNames=speciesNames,&
          & species=species)
      @:ASSERT(errStatus%hasError())
      @:ASSERT(errStatus%code == errors%invalidSelector)
    #:endblock

    #:block TEST("outOfBoundsMinus")
      call getIndexSelection("-17", fullRange, selected, errStatus, speciesNames=speciesNames,&
          & species=species)
      @:ASSERT(errStatus%hasError())
      @:ASSERT(errStatus%code == errors%invalidSelector)
    #:endblock

    #:block TEST("outOfBoundsMinusZeroInd")
      call getIndexSelection("-17", fullRange, selected, errStatus, speciesNames=speciesNames,&
          & species=species)
      @:ASSERT(errStatus%hasError())
      @:ASSERT(errStatus%code == errors%invalidSelector)
    #:endblock

    #:block TEST("negativeIndexWithoutBackward")
      call getIndexSelection("-6:-2", fullRange, selected, errStatus, speciesNames=speciesNames,&
        & species=species, backwardIndexing=.false.)
      @:ASSERT(errStatus%hasError())
      @:ASSERT(errStatus%code == errors%invalidSelector)
    #:endblock

    #:block TEST("outOfSubrangeBounds")
      call getIndexSelection("1:10", [9, 16], selected(9:16), errStatus,&
          & speciesNames=speciesNames, species=species(9:16), indexRange=fullRange)
      @:ASSERT(errStatus%hasError())
      @:ASSERT(errStatus%code == errors%invalidSelector)
    #:endblock

    #:block TEST("doubleOperator")
      call getIndexSelection("1 && 2", fullRange, selected, errStatus, speciesNames=speciesNames,&
          & species=species)
      @:ASSERT(errStatus%hasError())
      @:ASSERT(errStatus%code == errors%syntaxError)
    #:endblock

    #:block TEST("unclosedParenthesis")
      call getIndexSelection("(1 & 2", fullRange, selected, errStatus, speciesNames=speciesNames,&
          & species=species)
      @:ASSERT(errStatus%hasError())
      @:ASSERT(errStatus%code == errors%syntaxError)
    #:endblock

    #:block TEST("unopenedParenthesis")
      call getIndexSelection("1 & 2)", fullRange, selected, errStatus, speciesNames=speciesNames,&
          & species=species)
      @:ASSERT(errStatus%hasError())
      @:ASSERT(errStatus%code == errors%syntaxError)
    #:endblock

  #:endblock TEST_FIXTURE


  #:block TEST_FIXTURE("numeric")

    integer, parameter :: fullRange(2) = [1, 10]
    integer, parameter :: nItem = 10
    logical :: selected(nItem)
    type(TStatus) :: errStatus

  #:contains

     #:block TEST("explicit")
       call getIndexSelection("5 6 7", fullRange, selected, errStatus)
       @:ASSERT(errStatus%isOk())
       @:ASSERT(checkSelection(selected, fullRange, [5, 6, 7]))
      #:endblock

      #:block TEST("range")
        call getIndexSelection("-4:-2", fullRange, selected, errStatus)
        @:ASSERT(errStatus%isOk())
        @:ASSERT(checkSelection(selected, fullRange, [7, 8, 9]))
       #:endblock

  #:endblock TEST_FIXTURE


  #:block TEST_FIXTURE("negativeNumeric")

    integer, parameter :: fullRange(2) = [-2, 2]
    integer, parameter :: nItem = 5
    logical :: selected(nItem)
    type(TStatus) :: errStatus

  #:contains

     #:block TEST("explicit")
       call getIndexSelection("-1 0 1", fullRange, selected, errStatus)
       @:ASSERT(errStatus%isOk())
       @:ASSERT(checkSelection(selected, fullRange, [-1, 0, 1]))
      #:endblock

      #:block TEST("range")
        call getIndexSelection("-1:1", fullRange, selected, errStatus)
        @:ASSERT(errStatus%isOk())
        @:ASSERT(checkSelection(selected, fullRange, [-1, 0, 1]))
       #:endblock

  #:endblock TEST_FIXTURE


  !> Check whether the selected indices in a selection array match the expected ones
  function checkSelection(selected, fullRange, expectedIndices) result(matches)
    logical, intent(in) :: selected(:)
    integer, intent(in) :: fullRange(:)
    integer, intent(in) :: expectedIndices(:)
    logical :: matches

    integer, allocatable :: indices(:)
    integer :: ii

    matches = .false.
    indices = pack([(ii, ii = fullRange(1), fullRange(2))], selected)
    if (size(indices) /= size(expectedIndices)) return
    matches = all(indices == expectedIndices)

  end function checkSelection

#:endblock TEST_SUITE


@:TEST_DRIVER()
