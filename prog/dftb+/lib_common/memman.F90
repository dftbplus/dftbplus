!--------------------------------------------------------------------------------------------------!
!  DFTB+: general package for performing fast atomistic simulations                                !
!  Copyright (C) 2006 - 2020  DFTB+ developers group                                               !
!                                                                                                  !
!  See the LICENSE file for terms of usage and distribution.                                       !
!--------------------------------------------------------------------------------------------------!

!> Contains various constants for memory management
module dftbp_memman
  use dftbp_accuracy

  implicit none

  private

  public :: incrmntOfArray


  !> If the space for an array gets to small and has to be reallocated,
  !> new_size = arrayIncrement * old_size.
  !> Setting it too low causes a lot of realloc operations to occur!.
  real(dp), parameter :: arrayIncrement = 2.0_dp

contains


  !> figures out how much larger an array should be to minimize reallocations in future if the array
  !> grows more
  function incrmntOfArray(currentSize)
    integer :: incrmntOfArray

    !> current array size
    integer, intent(in) :: currentSize

    incrmntOfArray = int(real(currentSize, dp) * arrayIncrement)

  end function incrmntOfArray

end module dftbp_memman
