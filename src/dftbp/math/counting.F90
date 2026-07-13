!--------------------------------------------------------------------------------------------------!
!  DFTB+: general package for performing fast atomistic simulations                                !
!  Copyright (C) 2006 - 2025  DFTB+ developers group                                               !
!                                                                                                  !
!  See the LICENSE file for terms of usage and distribution.                                       !
!--------------------------------------------------------------------------------------------------!

#:include 'common.fypp'

!> Various types of value counting routines
module dftbp_math_counting
  use dftbp_common_accuracy, only : dp
  implicit none

  private
  public :: unique


contains


  !> Function to count number of unique elements in a sorted array of value greater than 0 and place
  !! them at the start of the array in order
  !! To do: check that the elements are in sorted order, and generalise for decreasing order as well
  !! as increasing
  function unique(array, arraySize) result(nUnique)

    !> Array to make unique
    integer, intent(inout) :: array(:)

    !> Constrains the effect of the subroutine on the first n elements, where n is the value for
    !! arraySize (default: size(array))
    integer, intent(in), optional :: arraySize

    !> Number of unique elements
    integer :: nUnique

    integer :: ii, ij, nn

    if (present(arraySize)) then
      nn = arraySize
    else
      nn = size(array)
    end if

    @:ASSERT(nn >= 1 )
    @:ASSERT(nn <= size(array))
    @:ASSERT(all(array(:nn) > 0))

    ii = 1
    do ij = 2, nn
      if (array(ij) /= array(ii)) then
        ii = ii + 1
        array(ii) = array(ij)
      end if
    end do
    nUnique = ii

  end function unique

end module dftbp_math_counting
