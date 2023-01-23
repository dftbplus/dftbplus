!--------------------------------------------------------------------------------------------------!
!  DFTB+: general package for performing fast atomistic simulations                                !
!  Copyright (C) 2006 - 2023  DFTB+ developers group                                               !
!                                                                                                  !
!  See the LICENSE file for terms of usage and distribution.                                       !
!--------------------------------------------------------------------------------------------------!

#:include "fytest.fypp"

#:block TEST_SUITE("memman")
  use dftbp_common_accuracy, only : dp
  use dftbp_common_memman, only : TAlignedArray
  use, intrinsic :: iso_c_binding, only : c_associated
  implicit none

#:contains

  #:block TEST_FIXTURE("allocateAligned")

    type(TAlignedArray) :: aligned

  #:contains

    #:block TEST("allocate32Bytes")
      call aligned%allocateAligned(7, 32)

      @:ASSERT(aligned%isAllocated)
      @:ASSERT(aligned%alignment == 32)

      @:ASSERT(aligned%size == 7)
      @:ASSERT(size(aligned%array, dim=1) == 7)
      @:ASSERT(lbound(aligned%array, dim=1) == 1)
      @:ASSERT(ubound(aligned%array, dim=1) == 7)

      aligned%array(:) = 1.0_dp
      aligned%array(7) = 42.0_dp
      @:ASSERT(aligned%array(7) == 42.0_dp)
    #:endblock

    #:block TEST("allocateDefault64Bytes")
      call aligned%allocateAligned(7)

      @:ASSERT(c_associated(aligned%memoryPointer))
      @:ASSERT(aligned%isAllocated)
      @:ASSERT(aligned%alignment == 64)
    #:endblock

    #:block TEST("deallocate")
      @:ASSERT(.not. c_associated(aligned%memoryPointer))
      @:ASSERT(.not. aligned%isAllocated)

      call aligned%allocateAligned(7)
      @:ASSERT(c_associated(aligned%memoryPointer))
      @:ASSERT(aligned%isAllocated)

      call triggerDeallocation(aligned)
      @:ASSERT(.not. c_associated(aligned%memoryPointer))
      @:ASSERT(.not. aligned%isAllocated)
    #:endblock


  #:endblock TEST_FIXTURE

  subroutine triggerDeallocation(instance)

    type(TAlignedArray), intent(out) :: instance

  end subroutine triggerDeallocation

#:endblock TEST_SUITE


@:TEST_DRIVER()
