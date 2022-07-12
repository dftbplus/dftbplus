!--------------------------------------------------------------------------------------------------!
!  DFTB+: general package for performing fast atomistic simulations                                !
!  Copyright (C) 2006 - 2022  DFTB+ developers group                                               !
!                                                                                                  !
!  See the LICENSE file for terms of usage and distribution.                                       !
!--------------------------------------------------------------------------------------------------!

!> Contains various constants for memory management
module dftbp_common_memman
  use dftbp_common_accuracy, only : dp

  implicit none

  private
  public :: incrmntOfArray


contains


  !> figures out how much larger an array should be to minimize reallocations in future if the array
  !> grows more
  pure function incrmntOfArray(currentSize)
    integer :: incrmntOfArray

    !> current array size
    integer, intent(in) :: currentSize

    incrmntOfArray = currentSize + currentSize  / 2 + 1

  end function incrmntOfArray

end module dftbp_common_memman
