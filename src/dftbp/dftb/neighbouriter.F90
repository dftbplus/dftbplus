!--------------------------------------------------------------------------------------------------!
!  DFTB+: general package for performing fast atomistic simulations                                !
!  Copyright (C) 2006 - 2021  DFTB+ developers group                                               !
!                                                                                                  !
!  See the LICENSE file for terms of usage and distribution.                                       !
!--------------------------------------------------------------------------------------------------!
module dftbp_dftb_neighbouriter
  use dftbp_common_accuracy, only : dp
  implicit none

  private
  public :: TNeighbourIter, TNeighbourIterFact


  type, abstract :: TNeighbourIter
  contains
    procedure(TNeighbourIter_start), deferred :: start
    procedure(TNeighbourIter_get), deferred :: get
  end type TNeighbourIter


  type, abstract :: TNeighbourIterFact
  contains
    procedure(TNeighbourIterFact_getIterator), deferred :: getIterator
  end type TNeighbourIterFact


  abstract interface

    subroutine TNeighbourIter_start(this, iAtom, cutoff, includeSelf, chunkSize)
      import :: TNeighbourIter, dp
      implicit none
      class(TNeighbourIter), intent(inout) :: this
      integer, intent(in) :: iAtom
      real(dp), intent(in) :: cutoff
      logical, intent(in) :: includeSelf
      integer, intent(in) :: chunkSize
    end subroutine TNeighbourIter_start


    subroutine TNeighbourIter_get(this, nNeighs, iNeighs, distances2)
      import :: TNeighbourIter, dp
      implicit none
      class(TNeighbourIter), intent(inout) :: this
      integer, intent(out) :: nNeighs
      integer, optional, intent(out) :: iNeighs(:)
      real(dp), optional, intent(out) :: distances2(:)
    end subroutine TNeighbourIter_get


    subroutine TNeighbourIterFact_getIterator(this, maxCutoff, neighIter)
      import :: TNeighbourIterFact, TNeighbourIter, dp
      implicit none
      class(TNeighbourIterFact), intent(in) :: this
      real(dp), intent(in) :: maxCutoff
      class(TNeighbourIter), allocatable, intent(out) :: neighIter
    end subroutine TNeighbourIterFact_getIterator

  end interface

end module dftbp_dftb_neighbouriter
