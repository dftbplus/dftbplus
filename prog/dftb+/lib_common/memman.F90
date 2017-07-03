!--------------------------------------------------------------------------------------------------!
!  DFTB+: general package for performing fast atomistic simulations                                !
!  Copyright (C) 2017  DFTB+ developers group                                                      !
!                                                                                                  !
!  See the LICENSE file for terms of usage and distribution.                                       !
!--------------------------------------------------------------------------------------------------!

!!* Contains various constants for memory management
module memman
  use accuracy

  implicit none

  private

  public :: incrmntOfArray

  real(dp), parameter :: arrayIncrement = 2.0_dp !!* If the space for an
  !!* array gets to small and has to be reallocated,
  !!* new_size = arrayIncrement * old_size. Setting it too low causes
  !!* a lot of realloc operations!.

contains

  !!* figures out how much larger an array should be to minimize reallocations
  !!* in future if the array grows more
  !!* @param currentSize current array size
  function incrmntOfArray(currentSize)
    integer             :: incrmntOfArray
    integer, intent(in) :: currentSize

    incrmntOfArray = int(real(currentSize, dp) * arrayIncrement)

  end function incrmntOfArray

end module memman
