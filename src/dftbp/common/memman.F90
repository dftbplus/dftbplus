!--------------------------------------------------------------------------------------------------!
!  DFTB+: general package for performing fast atomistic simulations                                !
!  Copyright (C) 2006 - 2023  DFTB+ developers group                                               !
!                                                                                                  !
!  See the LICENSE file for terms of usage and distribution.                                       !
!--------------------------------------------------------------------------------------------------!

#:include 'common.fypp'

!> Contains various constants for memory management
module dftbp_common_memman
  use dftbp_common_accuracy, only : dp
  use dftbp_io_message, only : error, warning
  use, intrinsic :: iso_c_binding, only : c_ptr, c_null_ptr, c_f_pointer,&
      & c_int, c_intptr_t, c_size_t

  implicit none

  private
  public :: incrmntOfArray, TAlignedArray

  !> Holds an allocated array being aligned in memory
  type TAlignedArray

    !> Actual array data
    real(dp), pointer :: array(:)

    !> Number of array elements
    integer :: size = 0

    !> Chosen alignment in bytes (default: 64 byted)
    integer :: alignment = 64

    !> Whether the array is allocated
    logical :: isAllocated = .false.

    !> Pointer to the allocated memory
    type(c_ptr) :: memoryPointer = c_null_ptr

  contains

    procedure :: allocateAligned => TAlignedArray_allocateAligned
    final :: TAlignedArray_finalize

  end type TAlignedArray


  !> Bound to 'posix_memalign' to allocate aligned memory
  interface
    function posixMemalign(ptr, alignment, size) result(error) bind(C, name="posix_memalign")
      import c_ptr, c_intptr_t, c_int
      type(c_ptr), intent(inout) :: ptr
      integer(c_intptr_t), intent(in), value :: alignment, size
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


  !> figures out how much larger an array should be to minimize reallocations in future if the array
  !> grows more
  pure function incrmntOfArray(currentSize)
    integer :: incrmntOfArray

    !> current array size
    integer, intent(in) :: currentSize

    incrmntOfArray = currentSize + currentSize  / 2 + 1

  end function incrmntOfArray


  !> Allocate an array aligned to the n byte boundary (typically used for aligning to 64 byte)
  subroutine TAlignedArray_allocateAligned(this, size, alignment)

    !> Instance of the aligned array type
    class(TAlignedArray), intent(inout) :: this

    !> Lower and upper bound for alignment
    integer, intent(in) :: size

    !> Alignment (default: 64 bytes)
    integer, intent(in), optional :: alignment

    if (this%isAllocated) then
      call error("Aligned array is already allocated")
    end if

    @:ASSERT(size > 0)

    if (present(alignment)) then
      @:ASSERT(alignment >= dp)
      @:ASSERT(mod(alignment, dp) == 0)

      this%alignment = alignment
    end if

    if (posixMemalign(this%memoryPointer, int(this%alignment, kind=c_size_t),&
          & int(size * dp, kind=c_size_t)) /= 0) then
      call error("Error during allocation of aligned array")
    end if

    this%size = size
    call c_f_pointer(this%memoryPointer, this%array, [size])

    this%isAllocated = .true.

  end subroutine TAlignedArray_allocateAligned


  !> Clean-up an aligned array by deallocating it
  subroutine TAlignedArray_finalize(this)

    !> Instance of the aligned array type
    type(TAlignedArray), intent(inout) :: this

    if (this%isAllocated) then
      nullify(this%array)
      call free(this%memoryPointer)
      this%memoryPointer = c_null_ptr
      this%isAllocated = .false.
    end if

  end subroutine TAlignedArray_finalize


end module dftbp_common_memman
