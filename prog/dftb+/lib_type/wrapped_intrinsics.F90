!> Implementes various ragged arrays.
module WrappedIntrinsics
  use Accuracy
  implicit none
  private

  public :: WrappedInt1, WrappedReal1, WrappedLogical1

  type :: WrappedInt1
    integer, allocatable :: data(:)
  end type WrappedInt1

  type :: WrappedReal1
    real(dp), allocatable :: data(:)
  end type WrappedReal1

  type :: WrappedLogical1
    logical, allocatable :: data(:)
  end type WrappedLogical1

end module WrappedIntrinsics
