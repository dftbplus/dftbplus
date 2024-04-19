!--------------------------------------------------------------------------------------------------!
!  DFTB+: general package for performing fast atomistic simulations                                !
!  Copyright (C) 2006 - 2023  DFTB+ developers group                                               !
!                                                                                                  !
!  See the LICENSE file for terms of usage and distribution.                                       !
!--------------------------------------------------------------------------------------------------!

!> Contains routines relating to duplicated numbers
module dftbp_math_duplicate
  use dftbp_math_sorting, only : heap_sort, unique
  implicit none

  private
  public :: isRepeated

contains

  !> Is there a repeated entry in an array
  function isRepeated(array)

    !> Array to check
    integer, intent(in) :: array(:)

    !> Are there repeating elements in the array
    logical :: isRepeated

    integer :: work(size(array))
    integer :: nUnique

    work(:) = array
    call heap_sort(work)
    nUnique = unique(work)

    isRepeated = .not.(nUnique == size(array))

  end function isRepeated

end module dftbp_math_duplicate
