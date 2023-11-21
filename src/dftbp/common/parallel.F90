!--------------------------------------------------------------------------------------------------!
!  DFTB+: general package for performing fast atomistic simulations                                !
!  Copyright (C) 2006 - 2023  DFTB+ developers group                                               !
!                                                                                                  !
!  See the LICENSE file for terms of usage and distribution.                                       !
!--------------------------------------------------------------------------------------------------!

#:include 'common.fypp'

!> Provides auxiliary routines for parallelization.
module dftbp_common_parallel

  implicit none

  private

  public :: getStartAndEndIndex


contains

  !> Returns the start and end index of an MPI process that calculates parts of a loop.
  subroutine getStartAndEndIndex(nElements, nProcs, iProc, iStart, iEnd)

    !> Array size to split
    integer, intent(in) :: nElements

    !> Number of available processes
    integer, intent(in) :: nProcs

    !> Current process index
    integer, intent(in) :: iProc

    !> Start and end index of current element range
    integer, intent(out) :: iStart, iEnd

    !! size of split index regions
    integer :: splitSize

    !! number of elements that exceed integer times nProcs
    integer :: offset

    @:ASSERT(iProc + 1 <= nProcs)

    splitSize = nElements / nProcs

    ! start and end indices assuming equal split sizes
    iStart = iProc * splitSize + 1
    iEnd = iStart + splitSize - 1

    ! distribute possible remainder to the ranges at the end
    offset = nProcs - mod(nElements, nProcs)
    if (iProc + 1 > offset) then
      iStart = iStart + iProc - offset
      iEnd = iEnd + iProc - offset + 1
    end if

  end subroutine getStartAndEndIndex

end module dftbp_common_parallel
