!--------------------------------------------------------------------------------------------------!
!  DFTB+: general package for performing fast atomistic simulations                                !
!  Copyright (C) 2006 - 2021  DFTB+ developers group                                               !
!                                                                                                  !
!  See the LICENSE file for terms of usage and distribution.                                       !
!--------------------------------------------------------------------------------------------------!

#:include "common.fypp"

!> Provides a neighour list iterator over a pre-generated static list of neighours.
module dftbp_dftb_staticneighiter
  use dftbp_common_accuracy, only : dp
  use dftbp_dftb_neighbouriter, only : TNeighbourIter, TNeighbourIterFact
  use dftbp_dftb_periodic, only : TNeighbourList
  implicit none

  private
  public :: TStaticNeighIter, TStaticNeighIter_init
  public :: TStaticNeighIterFact, TStaticNeighIterFact_init


  type, extends(TNeighbourIter) :: TStaticNeighIter
    private
    type(TNeighbourList), pointer :: pNeighList => null()
    integer :: iAtom = 0
    real(dp) :: cutoff2 = 0.0_dp
    integer :: lastServed = -1
    integer :: chunkSize = 0
  contains
    procedure :: start => TStaticNeighIter_start
    procedure :: get => TStaticNeighIter_get
  end type TStaticNeighIter


  type, extends(TNeighbourIterFact) :: TStaticNeighIterFact
    type(TNeighbourList), pointer :: pNeighList => null()
  contains
    procedure :: getIterator => TStaticNeighIterFact_getIterator
  end type TStaticNeighIterFact


contains

  subroutine TStaticNeighIter_init(this, pNeighList)
    type(TStaticNeighIter), intent(out) :: this
    type(TNeighbourList), pointer, intent(in) :: pNeighList

    this%pNeighList => pNeighList

  end subroutine TStaticNeighIter_init


  subroutine TStaticNeighIter_start(this, iAtom, cutoff, includeSelf, chunkSize)
    class(TStaticNeighIter), intent(inout) :: this
    integer, intent(in) :: iAtom
    real(dp), intent(in) :: cutoff
    logical, intent(in) :: includeSelf
    integer, intent(in) :: chunkSize

    this%iAtom = iAtom
    this%cutoff2 = cutoff**2
    if (includeSelf) then
      this%lastServed = -1
    else
      this%lastServed = 0
    end if
    this%chunkSize = chunkSize

  end subroutine TStaticNeighIter_start


  subroutine TStaticNeighIter_get(this, nNeighs, iNeighs, distances2)
    class(TStaticNeighIter), intent(inout) :: this
    integer, intent(out) :: nNeighs
    integer, optional, intent(out) :: iNeighs(:)
    real(dp), optional, intent(out) :: distances2(:)

    integer :: iFirst, iLast

    do nNeighs = 1, min(this%chunkSize, this%pNeighList%nNeighbour(this%iAtom) - this%lastServed)
      if (this%pNeighList%neighDist2(this%lastServed + nNeighs, this%iAtom) > this%cutoff2) exit
    end do
    nNeighs = nNeighs - 1
    iFirst = this%lastServed + 1
    iLast = this%lastServed + nNeighs
    if (present(iNeighs)) then
      iNeighs(1 : nNeighs) = this%pNeighList%iNeighbour(iFirst : iLast, this%iAtom)
    end if
    if (present(distances2)) then
      distances2(1 : nNeighs) = this%pNeighList%neighDist2(iFirst : iLast, this%iAtom)
    end if
    this%lastServed = iLast

  end subroutine TStaticNeighIter_get


  subroutine TStaticNeighIterFact_init(this, pNeighList)
    type(TStaticNeighIterFact), intent(out) :: this
    type(TNeighbourList), pointer, intent(in) :: pNeighList

    this%pNeighList => pNeighList

  end subroutine TStaticNeighIterFact_init


  function TStaticNeighIterFact_getIterator(this, maxCutoff) result(neighIter)
    class(TStaticNeighIterFact), intent(in) :: this
    real(dp), intent(in) :: maxCutoff
    class(TNeighbourIter), allocatable :: neighIter

    type(TStaticNeighIter), allocatable :: staticNeighIter

    @:ASSERT(associated(this%pNeighList))
    @:ASSERT(this%pNeighList%cutoff >= maxCutoff)
    allocate(staticNeighIter)
    call TStaticNeighIter_init(staticNeighIter, this%pNeighList)
    call move_alloc(staticNeighIter, neighIter)

  end function TStaticNeighIterFact_getIterator


end module dftbp_dftb_staticneighiter
