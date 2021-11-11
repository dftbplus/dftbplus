!--------------------------------------------------------------------------------------------------!
!  DFTB+: general package for performing fast atomistic simulations                                !
!  Copyright (C) 2006 - 2021  DFTB+ developers group                                               !
!                                                                                                  !
!  See the LICENSE file for terms of usage and distribution.                                       !
!--------------------------------------------------------------------------------------------------!
module dftbp_dftb_staticneighmap
  use dftbp_common_accuracy, only : dp, tolSameDist2
  use dftbp_dftb_neighbourmap, only : TNeighbourMap
  use dftbp_dftb_periodic, only : TNeighbourList
  use dftbp_math_bisect, only : bisection
  implicit none

  private
  public :: TStaticNeighMap, TStaticNeighMap_init


  type, extends(TNeighbourMap) :: TStaticNeighMap
    private
    type(TNeighbourList), pointer :: pNeighList => null()
  contains
    procedure :: getNeighbours => TStaticNeighMap_getNeighbours
  end type TStaticNeighMap


contains


  subroutine TStaticNeighMap_init(this, pNeighList)
    type(TStaticNeighMap), intent(out) :: this
    type(TNeighbourList), pointer, intent(in) :: pNeighList

    this%pNeighList => pNeighList

  end subroutine TStaticNeighMap_init


  subroutine TStaticNeighMap_getNeighbours(this, iAtom, cutoff, includeSelf, iNeighs, distances2)
    class(TStaticNeighMap), intent(in) :: this
    integer, intent(in) :: iAtom
    real(dp), intent(in) :: cutoff
    logical, optional, intent(in) :: includeSelf
    integer, allocatable, optional, intent(out) :: iNeighs(:)
    real(dp), allocatable, optional, intent(out) :: distances2(:)

    integer :: nNeighs
    integer :: iFirst, iLast
    logical :: includeSelf_

    call bisection(nNeighs,&
        & this%pNeighList%neighDist2(1 : this%pNeighList%nNeighbour(iAtom), iAtom), cutoff**2,&
        & tolSameDist2)
    includeSelf_ = .false.
    if (present(includeSelf)) includeSelf_ = includeSelf
    iFirst = 1
    iLast = nNeighs
    if (includeSelf_) then
      iFirst = 0
      nNeighs = nNeighs + 1
    end if
    if (present(iNeighs)) iNeighs = this%pNeighList%iNeighbour(iFirst : iLast, iAtom)
    if (present(distances2)) distances2 = this%pNeighList%neighDist2(iFirst : iLast, iAtom)

  end subroutine TStaticNeighMap_getNeighbours


end module dftbp_dftb_staticneighmap