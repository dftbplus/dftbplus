!--------------------------------------------------------------------------------------------------!
!  DFTB+: general package for performing fast atomistic simulations                                !
!  Copyright (C) 2006 - 2025  DFTB+ developers group                                               !
!                                                                                                  !
!  See the LICENSE file for terms of usage and distribution.                                       !
!--------------------------------------------------------------------------------------------------!

#:include 'common.fypp'

!> Various types of value counting routines
module dftbp_math_counting
  use dftbp_common_accuracy, only : dp, elecTolMax, rsp
  use dftbp_dftb_etemp, only : fillingTypes
  implicit none

  private
  public :: unique, filledStates, emptyStates


contains


  !> Function to count number of unique elements in a sorted array of value greater than 0 and place
  !! them at the start of the array in order
  !! To do: check that the elements are in sorted order, and generalise for decreasing order as well
  !! as increasing
  function unique(array, arraySize) result(nUnique)

    !> Array to make unique
    integer, intent(inout) :: array(:)

    !> Constrains the effect of the subroutine on the first n elements, where n is the value for
    !! arraySize (default: size(array))
    integer, intent(in), optional :: arraySize

    !> Number of unique elements
    integer :: nUnique

    integer :: ii, ij, nn

    if (present(arraySize)) then
      nn = arraySize
    else
      nn = size(array)
    end if

    @:ASSERT(nn >= 1 )
    @:ASSERT(nn <= size(array))
    @:ASSERT(all(array(:nn) > 0))

    ii = 1
    do ij = 2, nn
      if (array(ij) /= array(ii)) then
        ii = ii + 1
        array(ii) = array(ij)
      end if
    end do
    nUnique = ii

  end function unique


  !> Finds last electronic state with at least some occupation
  subroutine filledStates(nFilled, fillings, iFilling, tol)

    !> Number of (partially) filled states for each k-point and spin channel
    integer, intent(out), allocatable :: nFilled(:,:)

    !> Fillings for each level at k-points and spin channels
    real(dp), intent(in) :: fillings(:,:,:)

    !> Choice of electron distribution function, defaults to Fermi if not present
    integer, intent(in), optional :: iFilling

    !> Tolerance to detect empty states, defaults to elecTolMax if not present
    real(dp), intent(in), optional :: tol

    integer :: iS, iK, iLev, iFilling_
    real(dp) :: tol_

    tol_ = elecTolMax
    if (present(tol)) tol_ = tol
    iFilling_ = fillingTypes%Fermi
    if (present(iFilling)) iFilling_ = iFilling

    allocate(nFilled(size(fillings, dim=3), size(fillings, dim=2)), source=0)

    do iS = 1, size(fillings, dim=3)
      do iK = 1, size(fillings, dim=2)
        lpLevels: do iLev = size(fillings, dim=1), 1, -1
          if ( fillings(iLev, iK, iS) > tol_ ) then
            nFilled(iS, iK) = iLev
            exit lpLevels
          end if
        end do lpLevels
      end do
    end do

  end subroutine filledStates


  !> First level of each k/spin channel; which is not filled to capacity
  subroutine emptyStates(nEmpty, fillings, maxFill, iFilling, tol)

    !> First (partially) empty state for each k-point and spin channel
    integer, intent(out), allocatable :: nEmpty(:,:)

    !> Fillings of levels
    real(dp), intent(in) :: fillings(:,:,:)

    !> Maximum occupation for levels
    real(dp), intent(in) :: maxFill

    !> Choice of electron distribution function, defaults to Fermi if not present
    integer, intent(in), optional :: iFilling

    !> Tolerance to detect empty
    real(dp), intent(in), optional :: tol

    integer :: iS, iK, iLev, iFilling_
    real(dp) :: tol_

    tol_ = elecTolMax
    if (present(tol)) tol_ = tol
    iFilling_ = fillingTypes%Fermi
    if (present(iFilling)) iFilling_ = iFilling

    ! Start by assuming each channel is empty, so first level is unoccupied
    allocate(nEmpty(size(fillings, dim=3), size(fillings, dim=2)), source = 1)
    ! Check where there is slightly empty states
    do iS = 1, size(fillings, dim=3)
      do iK = 1, size(fillings, dim=2)
        lpLevel: do iLev = 1, size(fillings, dim=1)
          if ( abs( fillings(iLev, iK, iS) - maxFill ) > epsilon(1.0)) then
            ! this is a partially filled level
            nEmpty(iS, iK) = iLev
            exit lpLevel
          end if
        end do lpLevel
      end do
    end do

  end subroutine emptyStates


end module dftbp_math_counting
