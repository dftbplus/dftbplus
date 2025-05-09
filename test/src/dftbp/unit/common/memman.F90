!--------------------------------------------------------------------------------------------------!
!  DFTB+: general package for performing fast atomistic simulations                                !
!  Copyright (C) 2006 - 2025  DFTB+ developers group                                               !
!                                                                                                  !
!  See the LICENSE file for terms of usage and distribution.                                       !
!--------------------------------------------------------------------------------------------------!

#:include "fortuno_serial.fypp"

module test_common_memman
  use fortuno_serial, only : suite => serial_suite_item, test_list
  use dftbp_common_accuracy, only : dp
  use dftbp_common_memman, only : TAlignedArray
  $:FORTUNO_SERIAL_IMPORTS()
  implicit none

  private
  public :: tests

contains


  $:TEST("allocate32Bytes")
    type(TAlignedArray) :: aligned
    real(dp), pointer :: array(:)

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
  $:END_TEST()


  $:TEST("allocateDefault64Bytes")
    type(TAlignedArray) :: aligned
    real(dp), pointer :: array(:)

    call aligned%allocate(7)
    call aligned%getArray(array)

    @:ASSERT(associated(array))
    @:ASSERT(aligned%alignment == 64)
  $:END_TEST()


  $:TEST("deallocate")
    type(TAlignedArray) :: aligned
    real(dp), pointer :: array(:)

    call aligned%getArray(array)
    @:ASSERT(.not. associated(array))

    call aligned%allocate(7)
    call aligned%getArray(array)
    @:ASSERT(associated(array))

    call triggerDeallocation_(aligned)
    call aligned%getArray(array)
    @:ASSERT(.not. associated(array))
  $:END_TEST()


  $:TEST("reallocate")
    type(TAlignedArray) :: aligned, aligned2
    real(dp), pointer :: array(:)

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
  $:END_TEST()


  function tests()
    type(test_list) :: tests

    tests = test_list([&
        suite("memman", test_list([&
            $:TEST_ITEMS()
        ]))&
    ])
    $:STOP_ON_MISSING_TEST_ITEMS()

  end function tests


  subroutine triggerDeallocation_(instance)

    type(TAlignedArray), intent(out) :: instance

  end subroutine triggerDeallocation_

end module test_common_memman
