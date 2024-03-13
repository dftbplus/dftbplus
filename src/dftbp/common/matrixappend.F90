!--------------------------------------------------------------------------------------------------!
!  DFTB+: general package for performing fast atomistic simulations                                !
!  Copyright (C) 2006 - 2024  DFTB+ developers group                                               !
!                                                                                                  !
!  See the LICENSE file for terms of usage and distribution.                                       !
!--------------------------------------------------------------------------------------------------!

#:include 'common.fypp'

#:set FLAVOURS = [('logical', 'logical', ''), ('integer', 'int', ''), ('real', 'real', '(dp)'),&
  & ('complex', 'cmplx', '(dp)')]

!> Implements various routines to append elements to an array by copying and reallocation.
module dftbp_common_matrixappend
  use dftbp_common_accuracy, only : dp
  implicit none

  private
#:for _, NAME, _ in FLAVOURS

  public :: appendToArray1d_${NAME}$
  public :: appendToArray2d_${NAME}$

#:endfor

contains

#:for TYPE, NAME, PREC in FLAVOURS

  !> Appends data to one-dimensional, ${TYPE}$-valued array.
  subroutine appendToArray1d_${NAME}$(array, data)

    !> Array to extend
    ${TYPE}$${PREC}$, intent(inout), allocatable :: array(:)

    !> Data to add
    ${TYPE}$${PREC}$, intent(in) :: data

    ! Temporary storage
    ${TYPE}$${PREC}$, allocatable :: tmp(:)

    ! Original number of entries in array
    integer :: nElem

    if (allocated(array)) then
      nElem = size(array)
      allocate(tmp(nElem + 1))
      tmp(1:nElem) = array
      tmp(nElem + 1) = data
      deallocate(array)
      call move_alloc(tmp, array)
    else
      allocate(array(1))
      array(1) = data
    end if

  end subroutine appendToArray1d_${NAME}$


  !> Appends data to two-dimensional, ${TYPE}$-valued array.
  subroutine appendToArray2d_${NAME}$(array, data)

    !> Array to extend
    ${TYPE}$${PREC}$, intent(inout), allocatable :: array(:,:)

    !> Data to add
    ${TYPE}$${PREC}$, intent(in) :: data(:)

    !! Temporary storage
    ${TYPE}$${PREC}$, allocatable :: tmp(:,:)

    !! Original number of entries in array
    integer :: nElem

    if (allocated(array)) then
      @:ASSERT(size(data, dim=1) == size(array, dim=1))
      nElem = size(array, dim=2)
      allocate(tmp(size(array, dim=1), nElem + 1))
      tmp(:, :nElem) = array
      tmp(:, nElem + 1) = data
      deallocate(array)
      call move_alloc(tmp, array)
    else
      allocate(array(size(data, dim=1), 1))
      array(:, 1) = data
    end if

  end subroutine appendToArray2d_${NAME}$

#:endfor

end module dftbp_common_matrixappend
