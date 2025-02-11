!--------------------------------------------------------------------------------------------------!
!  DFTB+: general package for performing fast atomistic simulations                                !
!  Copyright (C) 2006 - 2025  DFTB+ developers group                                               !
!                                                                                                  !
!  See the LICENSE file for terms of usage and distribution.                                       !
!--------------------------------------------------------------------------------------------------!

#:include "fytest.fypp"

#:block TEST_SUITE("memman")
  use dftbp_common_accuracy, only : dp
  use dftbp_common_memman, only : TAlignedArray
  implicit none

#:contains

  #:block TEST_FIXTURE("allocate")

    type(TAlignedArray) :: aligned, aligned2
    real(dp), pointer :: array(:)

  #:contains

    #:block TEST("allocate32Bytes")
      call aligned%allocate(7, 32)
      call aligned%getArray(array)

      @:ASSERT(associated(array))
      @:ASSERT(aligned%alignment == 32)

      @:ASSERT(aligned%size == 7)
      @:ASSERT(size(array) == 7)
      @:ASSERT(size(array, dim=1) == 7)
      @:ASSERT(lbound(array, dim=1) == 1)
      @:ASSERT(ubound(array, dim=1) == 7)

      array(:) = 1.0_dp
      array(7) = 42.0_dp
      nullify(array)
      call aligned%getArray(array)
      @:ASSERT(array(7) == 42.0_dp)
    #:endblock

    #:block TEST("allocateDefault64Bytes")
      call aligned%allocate(7)
      call aligned%getArray(array)

      @:ASSERT(associated(array))
      @:ASSERT(aligned%alignment == 64)
    #:endblock

    #:block TEST("deallocate")
      call aligned%getArray(array)
      @:ASSERT(.not. associated(array))

      call aligned%allocate(7)
      call aligned%getArray(array)
      @:ASSERT(associated(array))

      call triggerDeallocation(aligned)
      call aligned%getArray(array)
      @:ASSERT(.not. associated(array))
    #:endblock

    #:block TEST("reallocate")
      call aligned%allocate(3)
      call aligned%getArray(array)
      array(1) = 3.0_dp

      call aligned2%allocate(4)
      call aligned2%getArray(array)
      array(1) = 4.0_dp

      aligned = aligned2
      @:ASSERT(aligned%size == 4)
      call aligned%getArray(array)
      @:ASSERT(array(1) == 4.0_dp)
    #:endblock


  #:endblock TEST_FIXTURE

  subroutine triggerDeallocation(instance)

    type(TAlignedArray), intent(out) :: instance

  end subroutine triggerDeallocation

#:endblock TEST_SUITE


@:TEST_DRIVER()
