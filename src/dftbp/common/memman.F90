!--------------------------------------------------------------------------------------------------!
!  DFTB+: general package for performing fast atomistic simulations                                !
!  Copyright (C) 2006 - 2025  DFTB+ developers group                                               !
!                                                                                                  !
!  See the LICENSE file for terms of usage and distribution.                                       !
!--------------------------------------------------------------------------------------------------!

#:include 'common.fypp'

!> Contains various constants for memory management
module dftbp_common_memman
  use, intrinsic :: iso_c_binding, only : c_associated, c_f_pointer, c_int, c_null_ptr, c_ptr,&
      & c_size_t, c_sizeof
  use dftbp_common_accuracy, only : dp
  use dftbp_io_message, only : error

  implicit none

  private
  public :: incrmntOfArray, TAlignedArray

  !> Holds an allocated array being aligned in memory
  type TAlignedArray

    !> Number of array elements
    integer :: size = 0

    !> Chosen alignment in bytes (default: 64 bytes)
    integer :: alignment = 64

    !> Pointer to the allocated memory
    type(c_ptr), private :: memoryPointer_ = c_null_ptr

  contains

    procedure :: allocate => TAlignedArray_allocate
    procedure :: getArray => TAlignedArray_getArray
    procedure, pass(this) :: TAlignedArray_assign
    generic :: assignment(=) => TAlignedArray_assign
    final :: TAlignedArray_finalize

  end type TAlignedArray


  !> Bound to 'posix_memalign' to allocate aligned memory
  interface
    function posixMemalign(ptr, alignment, size) result(error) bind(C, name="posix_memalign")
      import c_ptr, c_size_t, c_int
      type(c_ptr), intent(inout) :: ptr
      integer(c_size_t), intent(in), value :: alignment, size
      integer(c_int) :: error
    end function
  end interface


  !> Bound to 'free' to deallocate memory again being previously allocated using 'posix_memalign'
  interface
    subroutine free(ptr) bind(C, name="free")
      import c_ptr
      type(c_ptr), value :: ptr
    end subroutine
  end interface


contains


  !> Figures out how much larger an array should be to minimize reallocations in future if the array
  !! grows more
  pure function incrmntOfArray(currentSize)
    integer :: incrmntOfArray

    !> Current array size
    integer, intent(in) :: currentSize

    incrmntOfArray = currentSize + currentSize  / 2 + 1

  end function incrmntOfArray


  !> Allocate an array aligned to the n byte boundary (typically used for aligning to 64 bytes)
  subroutine TAlignedArray_allocate(this, size, alignment)

    !> Instance of the aligned array type
    class(TAlignedArray), intent(inout) :: this

    !> Lower and upper bound for alignment
    integer, intent(in) :: size

    !> Alignment (default: 64 bytes)
    integer, intent(in), optional :: alignment

    integer :: dp_size

    if (c_associated(this%memoryPointer_)) then
      call error("Aligned array is already allocated")
    end if

    @:ASSERT(size > 0)

    dp_size = int(c_sizeof(0._dp))
    if (present(alignment)) then
      @:ASSERT(alignment >= dp_size)
      @:ASSERT(mod(alignment, dp_size) == 0)

      this%alignment = alignment
    end if

    if (posixMemalign(this%memoryPointer_, int(this%alignment, kind=c_size_t),&
          & int(size * dp_size, kind=c_size_t)) /= 0) then
      call error("Error during allocation of aligned array")
    end if

    @:ASSERT(c_associated(this%memoryPointer_))

    this%size = size

  end subroutine TAlignedArray_allocate


  !> Fetch the pointer to the array data
  subroutine TAlignedArray_getArray(this, array)

    !> Instance of the aligned array type
    class(TAlignedArray), intent(in) :: this

    !> Pointer to the array data
    real(dp), pointer, intent(out) :: array(:)

    if (c_associated(this%memoryPointer_)) then
      call c_f_pointer(this%memoryPointer_, array, [this%size])
    else
      array => null()
    end if

  end subroutine TAlignedArray_getArray


  !> On Assignment, deallocate, allocate, and copy array data
  elemental impure subroutine TAlignedArray_assign(this, other)

    !> Instance of the lhs of the assignment
    class(TAlignedArray), intent(out) :: this

    !> Instance of the rhs of the assignment
    type(TAlignedArray), intent(in) :: other

    real(dp), pointer :: thisArray(:), otherArray(:)

    call this%allocate(other%size, other%alignment)
    call this%getArray(thisArray)
    call other%getArray(otherArray)
    thisArray = otherArray

  end subroutine TAlignedArray_assign


  !> Clean-up an aligned array by deallocating it
  elemental impure subroutine TAlignedArray_finalize(this)

    !> Instance of the aligned array type
    type(TAlignedArray), intent(inout) :: this

    if (c_associated(this%memoryPointer_)) then
      call free(this%memoryPointer_)
      this%memoryPointer_ = c_null_ptr
    end if

  end subroutine TAlignedArray_finalize


end module dftbp_common_memman
