!--------------------------------------------------------------------------------------------------!
!  DFTB+: general package for performing fast atomistic simulations                                !
!  Copyright (C) 2006 - 2020  DFTB+ developers group                                               !
!                                                                                                  !
!  See the LICENSE file for terms of usage and distribution.                                       !
!--------------------------------------------------------------------------------------------------!

!> Implements various wrapped data types for use in creating ragged multi-dimensional arrays.
module dftbp_wrappedintr
  use dftbp_accuracy
  implicit none
  private

  public :: WrappedInt1, WrappedReal1, WrappedLogical1

  !> 1 dimensional integers
  type :: WrappedInt1
    integer, allocatable :: data(:)
  end type WrappedInt1

  !> 1 dimensional reals
  type :: WrappedReal1
    real(dp), allocatable :: data(:)
  end type WrappedReal1

  !> 1 dimensional logicals
  type :: WrappedLogical1
    logical, allocatable :: data(:)
  end type WrappedLogical1

end module dftbp_wrappedintr
