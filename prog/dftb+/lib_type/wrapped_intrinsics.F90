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

  public :: TWrappedInt1, TWrappedReal1, TWrappedLogical1

  !> 1 dimensional integers
  type :: TWrappedInt1
    integer, allocatable :: data(:)
  end type TWrappedInt1

  !> 1 dimensional reals
  type :: TWrappedReal1
    real(dp), allocatable :: data(:)
  end type TWrappedReal1

  !> 1 dimensional logicals
  type :: TWrappedLogical1
    logical, allocatable :: data(:)
  end type TWrappedLogical1

end module dftbp_wrappedintr
