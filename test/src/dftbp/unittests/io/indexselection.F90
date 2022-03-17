!--------------------------------------------------------------------------------------------------!
!  DFTB+: general package for performing fast atomistic simulations                                !
!  Copyright (C) 2006 - 2022  DFTB+ developers group                                               !
!                                                                                                  !
!  See the LICENSE file for terms of usage and distribution.                                       !
!--------------------------------------------------------------------------------------------------!

#:include "fytest.fypp"

#:block TEST_SUITE("indexselection")
  use dftbp_common_exception, only : TException
  use dftbp_io_indexselection, only : getIndexSelection, errors
  implicit none

#:contains

  #:block TEST_FIXTURE("atomSelection")

    character(2), parameter :: speciesNames(*) = [character(2) :: "Si", "C"]
    integer, parameter :: species(*) = [1, 2, 1, 2, 1, 2, 1, 2, 1, 2, 1, 2, 1, 2, 1, 2]
    integer, parameter :: nAtom = size(species)
    integer, parameter :: fullRange(2) = [1, nAtom]
    logical :: selected(nAtom)
    type(TException), allocatable :: exc

  #:contains

    #:block TEST("explicit")
      call getIndexSelection(exc, "1 2 4", fullRange, selected, speciesNames=speciesNames,&
          & species=species)
      @:ASSERT(.not. allocated(exc))
      @:ASSERT(checkSelection(selected, fullRange, [1, 2, 4]))
    #:endblock

    #:block TEST("allAsRange")
      call getIndexSelection(exc, "1:-1", fullRange, selected, speciesNames=speciesNames,&
          & species=species)
      @:ASSERT(.not. allocated(exc))
      @:ASSERT(all(selected))
    #:endblock

    #:block TEST("range")
      call getIndexSelection(exc, "2:6", fullRange, selected, speciesNames=speciesNames,&
          & species=species)
      @:ASSERT(.not. allocated(exc))
      @:ASSERT(checkSelection(selected, fullRange, [2, 3, 4, 5, 6]))
    #:endblock

    #:block TEST("negativeSubrange")
      call getIndexSelection(exc, "-6:-2", fullRange, selected, speciesNames=speciesNames,&
          & species=species)
      @:ASSERT(.not. allocated(exc))
      @:ASSERT(checkSelection(selected, fullRange, [11, 12, 13, 14, 15]))
    #:endblock

    #:block TEST("speciesName")
      call getIndexSelection(exc, "Si", fullRange, selected, speciesNames=speciesNames,&
          & species=species)
      @:ASSERT(.not. allocated(exc))
      @:ASSERT(checkSelection(selected, fullRange, [1, 3, 5, 7, 9, 11, 13, 15]))
    #:endblock

    #:block TEST("speciesNameDiffCase")
      call getIndexSelection(exc, "si", fullRange, selected, speciesNames=speciesNames,&
          & species=species)
      @:ASSERT(.not. allocated(exc))
      @:ASSERT(checkSelection(selected, fullRange, [1, 3, 5, 7, 9, 11, 13, 15]))
    #:endblock

    #:block TEST("excludeIndices")
      call getIndexSelection(exc, "!4:14", fullRange, selected, speciesNames=speciesNames,&
          & species=species)
      @:ASSERT(.not. allocated(exc))
      @:ASSERT(checkSelection(selected, fullRange, [1, 2, 3, 15, 16]))
    #:endblock

    #:block TEST("excludeSpecies")
      call getIndexSelection(exc, "!C", fullRange, selected, speciesNames=speciesNames,&
          & species=species)
      @:ASSERT(.not. allocated(exc))
      @:ASSERT(checkSelection(selected, fullRange, [1, 3, 5, 7, 9, 11, 13, 15]))
    #:endblock

    #:block TEST("excludesWithAnd")
      call getIndexSelection(exc, "!Si & !1:4 & !12:16", fullRange, selected,&
           speciesNames=speciesNames, species=species)
      @:ASSERT(.not. allocated(exc))
      @:ASSERT(checkSelection(selected, fullRange, [6, 8, 10]))
    #:endblock

    #:block TEST("notWithParenthesis")
      call getIndexSelection(exc, "!(Si 1:4 12:16)", fullRange, selected,&
           speciesNames=speciesNames, species=species)
      @:ASSERT(.not. allocated(exc))
      @:ASSERT(checkSelection(selected, fullRange, [6, 8, 10]))
    #:endblock

    #:block TEST("andOrPrecedence")
      call getIndexSelection(exc, "1:3 6:8 & Si 9:12", fullRange, selected,&
           speciesNames=speciesNames, species=species)
      @:ASSERT(.not. allocated(exc))
      @:ASSERT(checkSelection(selected, fullRange, [1, 2, 3, 7, 9, 10, 11, 12]))
    #:endblock

    #:block TEST("parenthesisPrecedence")
      call getIndexSelection(exc, "(1:3 6:8) & Si 9:12", fullRange, selected,&
          speciesNames=speciesNames, species=species)
      @:ASSERT(.not. allocated(exc))
      @:ASSERT(checkSelection(selected, fullRange, [1, 3, 7, 9, 10, 11, 12]))
    #:endblock

    #:block TEST("parenthesisPrecedence2")
      call getIndexSelection(exc, "1 (1:3 6:8) & !Si 9:12", fullRange, selected,&
          speciesNames=speciesNames, species=species)
      @:ASSERT(.not. allocated(exc))
      @:ASSERT(checkSelection(selected, fullRange, [1, 2, 6, 8, 9, 10, 11, 12]))
    #:endblock

    #:block TEST("parenthesisPrecedence3")
      call getIndexSelection(exc, "1 & (1:3 6:8) & !Si 9:12", fullRange, selected,&
          speciesNames=speciesNames, species=species)
      @:ASSERT(.not. allocated(exc))
      @:ASSERT(checkSelection(selected, fullRange, [9, 10, 11, 12]))
    #:endblock

    #:block TEST("zeroIndexingPositive")
      call getIndexSelection(exc, "0 2 3", [0, nAtom - 1], selected,&
          speciesNames=speciesNames, species=species)
      @:ASSERT(.not. allocated(exc))
      @:ASSERT(checkSelection(selected, fullRange, [1, 3, 4]))
    #:endblock

    #:block TEST("zeroIndexingNegative")
      call getIndexSelection(exc, "-16:-14", [0, nAtom - 1], selected,&
          speciesNames=speciesNames, species=species)
      @:ASSERT(.not. allocated(exc))
      @:ASSERT(checkSelection(selected, fullRange, [1, 2, 3]))
    #:endblock

    #:block TEST("subrange")
      call getIndexSelection(exc, "11:13", [9, 16], selected(9:16),&
          & speciesNames=speciesNames, species=species(9:16), indexRange=fullRange)
      @:ASSERT(.not. allocated(exc))
      @:ASSERT(checkSelection(selected(9:16), [9, 16], [11, 12, 13]))
    #:endblock

    #! Failing tests

    #:block TEST("selectWrongSpeciesName")
      call getIndexSelection(exc, "Mg", fullRange, selected, speciesNames=speciesNames,&
          & species=species)
      @:ASSERT(allocated(exc))
      call exc%deactivate()
      @:ASSERT(exc%code == errors%invalidSelector)
    #:endblock

    #:block TEST("outOfBounds")
      call getIndexSelection(exc, "17", fullRange, selected, speciesNames=speciesNames,&
          & species=species)
      @:ASSERT(allocated(exc))
      call exc%deactivate()
      @:ASSERT(exc%code == errors%invalidSelector)
    #:endblock

    #:block TEST("outOfBoundsPlusZeroInd")
      call getIndexSelection(exc, "16", [0, nAtom - 1], selected, speciesNames=speciesNames,&
          & species=species)
      @:ASSERT(allocated(exc))
      call exc%deactivate()
      @:ASSERT(exc%code == errors%invalidSelector)
    #:endblock

    #:block TEST("outOfBoundsMinus")
      call getIndexSelection(exc, "-17", fullRange, selected, speciesNames=speciesNames,&
          & species=species)
      @:ASSERT(allocated(exc))
      call exc%deactivate()
      @:ASSERT(exc%code == errors%invalidSelector)
    #:endblock

    #:block TEST("outOfBoundsMinusZeroInd")
      call getIndexSelection(exc, "-17", fullRange, selected, speciesNames=speciesNames,&
          & species=species)
      @:ASSERT(allocated(exc))
      call exc%deactivate()
      @:ASSERT(exc%code == errors%invalidSelector)
    #:endblock

    #:block TEST("negativeIndexWithoutBackward")
      call getIndexSelection(exc, "-6:-2", fullRange, selected, speciesNames=speciesNames,&
        & species=species, backwardIndexing=.false.)
      @:ASSERT(allocated(exc))
      call exc%deactivate()
      @:ASSERT(exc%code == errors%invalidSelector)
    #:endblock

    #:block TEST("outOfSubrangeBounds")
      call getIndexSelection(exc, "1:10", [9, 16], selected(9:16),&
          & speciesNames=speciesNames, species=species(9:16), indexRange=fullRange)
      @:ASSERT(allocated(exc))
      call exc%deactivate()
      @:ASSERT(exc%code == errors%invalidSelector)
    #:endblock

    #:block TEST("doubleOperator")
      call getIndexSelection(exc, "1 && 2", fullRange, selected, speciesNames=speciesNames,&
          & species=species)
      @:ASSERT(allocated(exc))
      call exc%deactivate()
      @:ASSERT(exc%code == errors%syntaxError)
    #:endblock

    #:block TEST("unclosedParenthesis")
      call getIndexSelection(exc, "(1 & 2", fullRange, selected, speciesNames=speciesNames,&
          & species=species)
      @:ASSERT(allocated(exc))
      call exc%deactivate()
      @:ASSERT(exc%code == errors%syntaxError)
    #:endblock

    #:block TEST("unopenedParenthesis")
      call getIndexSelection(exc, "1 & 2)", fullRange, selected, speciesNames=speciesNames,&
          & species=species)
      @:ASSERT(allocated(exc))
      call exc%deactivate()
      @:ASSERT(exc%code == errors%syntaxError)
    #:endblock

  #:endblock TEST_FIXTURE


  #:block TEST_FIXTURE("numeric")

    integer, parameter :: fullRange(2) = [1, 10]
    integer, parameter :: nItem = 10
    logical :: selected(nItem)
    type(TException), allocatable :: exc

  #:contains

     #:block TEST("explicit")
       call getIndexSelection(exc, "5 6 7", fullRange, selected)
       @:ASSERT(.not. allocated(exc))
       @:ASSERT(checkSelection(selected, fullRange, [5, 6, 7]))
      #:endblock

      #:block TEST("range")
        call getIndexSelection(exc, "-4:-2", fullRange, selected)
        @:ASSERT(.not. allocated(exc))
        @:ASSERT(checkSelection(selected, fullRange, [7, 8, 9]))
       #:endblock

  #:endblock TEST_FIXTURE


  #:block TEST_FIXTURE("negativeNumeric")

    integer, parameter :: fullRange(2) = [-2, 2]
    integer, parameter :: nItem = 5
    logical :: selected(nItem)
    type(TException), allocatable :: exc

  #:contains

     #:block TEST("explicit")
       call getIndexSelection(exc, "-1 0 1", fullRange, selected)
       @:ASSERT(.not. allocated(exc))
       @:ASSERT(checkSelection(selected, fullRange, [-1, 0, 1]))
      #:endblock

      #:block TEST("range")
        call getIndexSelection(exc, "-1:1", fullRange, selected)
        @:ASSERT(.not. allocated(exc))
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
