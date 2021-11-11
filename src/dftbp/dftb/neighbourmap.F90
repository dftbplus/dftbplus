!--------------------------------------------------------------------------------------------------!
!  DFTB+: general package for performing fast atomistic simulations                                !
!  Copyright (C) 2006 - 2021  DFTB+ developers group                                               !
!                                                                                                  !
!  See the LICENSE file for terms of usage and distribution.                                       !
!--------------------------------------------------------------------------------------------------!
module dftbp_dftb_neighbourmap
  use dftbp_common_accuracy, only : dp
  implicit none

  private
  public :: TNeighbourMap


  type, abstract :: TNeighbourMap
  contains
    procedure(TNeighbourMap_getNeighbours), deferred :: getNeighbours
  end type TNeighbourMap

  abstract interface

    subroutine TNeighbourMap_getNeighbours(this, iAtom, cutoff, includeSelf, iNeighs, distances2)
      import :: TNeighbourMap, dp
      implicit none
      class(TNeighbourMap), intent(in) :: this
      integer, intent(in) :: iAtom
      real(dp), intent(in) :: cutoff
      logical, optional, intent(in) :: includeSelf
      integer, allocatable, optional, intent(out) :: iNeighs(:)
      real(dp), allocatable, optional, intent(out) :: distances2(:)
    end subroutine TNeighbourMap_getNeighbours

  end interface

end module dftbp_dftb_neighbourmap