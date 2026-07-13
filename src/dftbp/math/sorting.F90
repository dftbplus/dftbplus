!--------------------------------------------------------------------------------------------------!
!  DFTB+: general package for performing fast atomistic simulations                                !
!  Copyright (C) 2006 - 2025  DFTB+ developers group                                               !
!                                                                                                  !
!  See the LICENSE file for terms of usage and distribution.                                       !
!--------------------------------------------------------------------------------------------------!

#:include 'common.fypp'

!> Various types of sorting routines
module dftbp_math_sorting
  use dftbp_common_accuracy, only : dp
  implicit none

  private
  public :: heap_sort, index_heap_sort, merge_sort


  !> Heap sort algorithm - O(N log(N)) time performance and in place, but not 'stable' in order of
  !! sorting
  interface heap_sort
    module procedure heap_sort_real
    module procedure heap_sort_int
  end interface heap_sort


contains


  !> Real case in-place heap sort
  !! Based on Numerical Recipes Software 1986-92
  pure subroutine heap_sort_real(array, tolerance)

    !> Array of values to be sorted
    real(dp), intent(inout) :: array(:)

    !> Tolerance for equality of two elements
    real(dp), intent(in), optional :: tolerance

    integer :: n, ir, ij, il, ii, ik
    real(dp) :: tmpReal
    real(dp) :: tol

    if (present(tolerance)) then
      tol = tolerance
    else
      tol = epsilon(1.0_dp)
    end if

    n = size(array)
    if (n <= 1) return
    il = n/2 + 1
    ir = n
    ik = 1
    do while (ik == 1)
      if (il > 1) then
        il = il - 1
        tmpReal = array(il)
      else
        tmpReal = array(ir)
        array(ir) = array(1)
        ir = ir - 1
        if(ir == 1)then
          array(1) = tmpReal
          return
        end if
      end if
      ii = il
      ij = 2 * il
      do while (ij <= ir)
        if (ij < ir) then
          if (array(ij) < array(ij+1) - tol) then
            ij = ij + 1
          end if
        end if
        if(tmpReal < array(ij) - tol) then
          array(ii) = array(ij)
          ii = ij
          ij = 2 * ij
        else
          ij = ir + 1
        end if
      end do
      array(ii) = tmpReal
    end do

  end subroutine heap_sort_real


  !> Integer case in-place heap sort
  !! based on Numerical Recipes Software 1986-92
  pure subroutine heap_sort_int(array)

    !> Array of values to be sorted
    integer, intent(inout) :: array(:)

    integer :: n, ii, ir, ij, il, ik
    integer :: tmpInt

    n = size(array)
    if (n <= 1) return

    il = n/2 + 1
    ir = n
    ik = 1
    do while (ik == 1)
      if (il > 1) then
        il = il - 1
        tmpInt = array(il)
      else
        tmpInt = array(ir)
        array(ir) = array(1)
        ir = ir - 1
        if(ir == 1)then
          array(1) = tmpInt
          return
        end if
      end if
      ii = il
      ij = 2 * il
      do while (ij <= ir)
        if (ij < ir) then
          if (array(ij) < array(ij+1)) then
            ij = ij + 1
          end if
        end if
        if(tmpInt < array(ij)) then
          array(ii) = array(ij)
          ii = ij
          ij = 2 * ij
        else
          ij = ir + 1
        end if
      end do
      array(ii) = tmpInt
    end do

  end subroutine  heap_sort_int


  !> Heap sort algorithm - O(N log(N)) performance, provides an index vector instead of re-ordering
  !! values, not a stable sort(!)  Comparisions are to within the provided tolerance. Routine based
  !! on Numerical Recipes Software 1986-92
  subroutine index_heap_sort(indx, array, tolerance)

    !> Indexing array on return
    integer, intent(out) :: indx(:)

    !> Array of values to be sorted
    real(dp), intent(in) :: array(:)

    !> Tolerance for equality of two elements
    real(dp), intent(in), optional :: tolerance

    integer :: n, ir, ij, il, ii, ik
    integer :: indxTmp
    real(dp) :: arrayTmp, tol

    @:ASSERT(size(array)==size(indx))

    if (present(tolerance)) then
      tol = tolerance
    else
      tol = epsilon(0.0_dp)
    end if

    do ii=1,size(indx)
      indx(ii) = ii
    end do
    n = size(array)
    if (n <= 1) return
    il=n/2+1
    ir=n
    ik = 1
    do while (ik == 1)
      if (il.gt.1) then
        il=il-1
        indxTmp=indx(il)
        arrayTmp=array(indxTmp)
      else
        indxTmp=indx(ir)
        arrayTmp=array(indxTmp)
        indx(ir)=indx(1)
        ir=ir-1
        if (ir.lt.1) then
          indx(1)=indxTmp
          return
        end if
      end if
      ii=il
      ij=2 * il
      do while (ij <= ir)
        if (ij < ir) then
          if (array(indx(ij)) < array(indx(ij+1)) - tol) then
            ij = ij + 1
          end if
        end if
        if(arrayTmp < array(indx(ij)) - tol) then
          indx(ii)=indx(ij)
          ii=ij
          ij=2*ij
        else
          ij = ir + 1
        end if
      end do
      indx(ii)=indxTmp
    end do

  end subroutine index_heap_sort


  !> Merge sort algorithm wrapper - O(N log(N)) performance, stable ordering but requires an O(N)
  !! workspace
  subroutine merge_sort(indx, arr, tolerance)

    !> Output index array
    integer, intent(out) :: indx(:)

    !> Data array to sort
    real(dp), intent(in) :: arr(:)

    !> Comparison equality tolerance between values
    real(dp), intent(in) :: tolerance

    integer :: n, ii
    integer, allocatable :: work(:)
    real(dp) :: tol

    n = size(arr)
    @:ASSERT(size(indx) >= n)

    ! Initialize an index array
    forall (ii = 1:n) indx(ii) = ii

    if (n > 1) then
      allocate(work(n))
      call merge_sort_index_track(indx, arr, work, 1, n, tolerance)
      deallocate(work)
    end if

  end subroutine merge_sort


  !> Recursive splitting of the index array
  recursive subroutine merge_sort_index_track(indx, arr, work, left, right, tol)

    !> Indexing array for data
    integer, intent(inout) :: indx(:)

    !> Data array to sort
    real(dp), intent(in) :: arr(:)

    !> Work array
    integer, intent(inout) :: work(:)

    !> Start of range
    integer, intent(in) :: left

    !> End of range
    integer, intent(in) :: right

    !> Tolerance for numerical comparisions
    real(dp), intent(in) :: tol

    integer :: midpoint

    if (left < right) then

      midpoint = left + (right - left) / 2 ! Overflow protected

      ! Recursively sort left and right parts
      call merge_sort_index_track(indx, arr, work, left, midpoint, tol)
      call merge_sort_index_track(indx, arr, work, midpoint + 1, right, tol)

      ! Merge the two sorted index segments
      call merge_index(arr, indx, work, left, midpoint, right, tol)

    end if

  end subroutine merge_sort_index_track


  !> Merges two sorted sub-segments using tolerance logic safely
  subroutine merge_index(arr, indx, work, left, midpoint, right, tol)

    !> Data array to sort
    real(dp), intent(in) :: arr(:)

    !> Indexing array for data
    integer, intent(inout) :: indx(:)

    !> Work array
    integer, intent(inout) :: work(:)

    !> Start of range
    integer, intent(in) :: left

    !> Middle of range
    integer, intent(in) :: midpoint

    !> End of range
    integer, intent(in) :: right

    !> Tolerance for numerical comparisions
    real(dp), intent(in) :: tol

    integer :: ii, jj, kk, nn
    real(dp) :: val_i, val_j

    ! Cache the current segment of array index into temporary workspace
    work(left:right) = indx(left:right)

    ii = left
    jj = midpoint + 1
    kk = left

    ! Merge back into indx()
    do while (ii <= midpoint .and. jj <= right)

      val_i = arr(work(ii))
      val_j = arr(work(jj))

      if ((val_j - val_i) > tol) then
        ! element ii is significantly smaller than element jj

        indx(kk) = work(ii)
        ii = ii + 1

      else if ((val_i - val_j) > tol) then
        ! element jj is significantly smaller than element i

        indx(kk) = work(jj)
        jj = jj + 1

      else
        ! elements are equal within tolerance -> enforce stability

        if (work(ii) <= work(jj)) then
          indx(kk) = work(ii)
          ii = ii + 1
        else
          indx(kk) = work(jj)
          jj = jj + 1
        end if

      end if

      kk = kk + 1

    end do

    ! Copy remaining elements from the left segment
    nn = midpoint - ii
    indx(kk:kk+nn) = work(ii:ii+nn)
    kk = kk + nn + 1

    ! Copy remaining elements from the right segment
    nn = right - jj
    indx(kk:kk+nn) = work(jj:jj+nn)

  end subroutine merge_index


end module dftbp_math_sorting
