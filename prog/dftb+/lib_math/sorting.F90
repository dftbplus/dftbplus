!--------------------------------------------------------------------------------------------------!
!  DFTB+: general package for performing fast atomistic simulations                                !
!  Copyright (C) 2018  DFTB+ developers group                                                      !
!                                                                                                  !
!  See the LICENSE file for terms of usage and distribution.                                       !
!--------------------------------------------------------------------------------------------------!

#:include 'common.fypp'

!> Various types of sorting routines, and related stuff
!> To do: add other algorithms, radix? definitely not quicksort though,
!> but adaptive heap sorts?
module sorting
  use assert
  use accuracy, only : dp
  implicit none
  private

  public :: heap_sort, index_heap_sort, merge_sort, unique


  !> Heap sort algorithm - O(N log(N)) time performance and in place, but not 'stable' in order of
  !> sorting
  interface heap_sort
    module procedure heap_sort_real
    module procedure heap_sort_int
  end interface heap_sort


  !> Heap sort algorithm - O(N log(N)) performance, provides an index vector instead of re-ordering
  !> values, again not stable
  interface index_heap_sort
    module procedure index_heap_sort_real
    module procedure index_heap_sort_int
  end interface index_heap_sort


  !> Merge sort algorithm - O(N log(N)) performance, stable ordering but requires an O(N)
  !> workspace. Versions with and without index array supplied.
  interface merge_sort
    module procedure merge_sort_int
    module procedure merge_sort_indx_int
    module procedure merge_sort_real
    module procedure merge_sort_indx_real
  end interface merge_sort


  !> Function to count number of unique elements in a sorted array of value greater than 0 and place
  !> them at the start of the array in order
  interface unique
    module procedure unique_int
  end interface unique

  ! non-public interfaces


  !> internal workhorse for merge sorts
  interface MergeSort
    module procedure MergeSort_int
    module procedure MergeSort_indx_int
    module procedure MergeSort_real
    module procedure MergeSort_indx_real
  end interface MergeSort


  !> merge two arrays together in order onto a third
  interface Merge
    module procedure Merge_int
    module procedure Merge_indx_int
    module procedure Merge_real
    module procedure Merge_indx_real
  end interface Merge

contains


  !> real case in-place heap sort
  !> Based on Numerical Recipes Software 1986-92
  subroutine heap_sort_real(array, tolerance)

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


  !> integer case in-place heap sort
  !> based on Numerical Recipes Software 1986-92
  subroutine heap_sort_int(array)

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


  !> Real case heap sort returning an index.
  !> based on Numerical Recipes Software 1986-92
  subroutine index_heap_sort_real(indx, array, tolerance)

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

  end subroutine index_heap_sort_real


  !> real case heap sort returning an index
  !> based on Numerical Recipes Software 1986-92
  subroutine index_heap_sort_int(indx, array)

    !> Indexing array on return
    integer, intent(out) :: indx(:)

    !> Array of values to be sorted
    integer, intent(in) :: array(:)

    integer :: n, ir, ij, il, ii, ik
    integer :: indxTmp, arrayTmp

    @:ASSERT(size(array)==size(indx))

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
          if (array(indx(ij)) < array(indx(ij+1))) then
            ij = ij + 1
          end if
        end if
        if(arrayTmp < array(indx(ij))) then
          indx(ii)=indx(ij)
          ii=ij
          ij=2*ij
        else
          ij = ir + 1
        end if
      end do
      indx(ii)=indxTmp
    end do

  end subroutine index_heap_sort_int


  !> Merge sort of integers
  subroutine merge_sort_int(array)

    !> vector to sort
    integer, intent(inout) :: array(:)

    integer, allocatable :: work(:)
    integer :: n

    n = size(array)

    allocate(work((n+1)/2))
    call mergeSort(array,n,work)

  end subroutine merge_sort_int


  !> Merge two arrays together in order onto a third
  subroutine merge_int(NA,NB,NC,A,B,C)

    !> first array of values
    integer, intent(in) :: NA

    !> elements in A
    integer, intent(in) :: NB

    !> second array of values
    integer, intent(in) :: NC

    !> elements in A
    integer, intent(in) :: A(NA)

    !> array to merge onto
    integer, intent(in) :: B(NB)

    !> elements in C
    integer, intent(inout) :: C(NC)

    integer :: I, J, K

    @:ASSERT((na+nb)==nc)

    I = 1; J = 1; K = 1;
    do while(I <= NA .and. J <= NB)
      if (A(I) <= B(J)) then
        C(K) = A(I)
        I = I+1
      else
        C(K) = B(J)
        J = J+1
      endif
      K = K + 1
    enddo
    do while (I <= NA)
      C(K) = A(I)
      I = I + 1
      K = K + 1
    enddo

  end subroutine merge_int


  !> Integer merge sort
  recursive subroutine mergeSort_int(A,N,T)

    !> array to sort
    integer, intent(inout) :: A(:)

    !> number of elements in array
    integer, intent(in) :: N

    !> workspace of at least (N+1)/2 size
    integer, intent (out) :: T(:)

    integer :: NA, NB, V

    if (N < 2) return

    if (N == 2) then
      if (A(1) > A(2)) then
        V = A(1)
        A(1) = A(2)
        A(2) = V
      endif
      return
    endif

    NA=(N+1)/2
    NB=N-NA

    call MergeSort(A,NA,T)
    call MergeSort(A(NA+1:),NB,T)

    if (A(NA) > A(NA+1)) then
      T(1:NA)=A(1:NA)
      call merge(NA,NB,N,T,A(NA+1:),A)
    endif

  end subroutine MergeSort_Int


  !> Merge sort of integers, using an array index instead of re-ordering
  subroutine merge_sort_indx_int(indx,array)

    !> index array for sort order
    integer, intent(out) :: indx(:)

    !> vector to sort
    integer, intent(in) :: array(:)

    integer, allocatable :: work(:,:), tmp(:,:)
    integer :: ii, n

    @:ASSERT(size(indx)==size(array))

    n = size(array)
    allocate(tmp(n,2))
    tmp(:,2) = array
    do ii = 1, n
      tmp(ii,1) = ii
    end do
    allocate(work((n+1)/2,2))
    call mergeSort(tmp,n,work)
    indx = tmp(:,1)

  end subroutine merge_sort_indx_int


  !> Merge two arrays together in order onto a third, where first dimension of both is index for
  !> original order and also value
  subroutine merge_indx_int(NA,NB,NC,A,B,C)

    !> first array of values
    integer, intent(in) :: NA

    !> elements in A
    integer, intent(in) :: NB

    !> second array of values
    integer, intent(in) :: NC

    !> elements in A
    integer, intent(in) :: A(NA,2)

    !> array to merge onto
    integer, intent(in) :: B(NB,2)

    !> elements in C
    integer, intent(inout) :: C(NC,2)

    integer :: I, J, K

    @:ASSERT((na+nb)==nc)

    I = 1; J = 1; K = 1;
    do while(I <= NA .and. J <= NB)
      if (A(I,2) <= B(J,2)) then
        C(K,:) = A(I,:)
        I = I+1
      else
        C(K,:) = B(J,:)
        J = J+1
      endif
      K = K + 1
    enddo
    do while (I <= NA)
      C(K,:) = A(I,:)
      I = I + 1
      K = K + 1
    enddo

  end subroutine merge_indx_int


  !> Integer merge sort, using an index
  recursive subroutine mergeSort_indx_int(A,N,T)

    !> array to sort, first element of first dimension is an index array, second element is actual
    !> value
    integer, intent(inout) :: A(:,:)

    !> N number of elements in array
    integer, intent(in) :: N

    !> workspace of at least (N+1)/2 size
    integer, intent (out) :: T(:,:)

    integer :: NA, NB, V(2)

    @:ASSERT(size(A,dim=2) == 2)
    @:ASSERT(size(T,dim=2) == 2)

    if (N < 2) return

    if (N == 2) then
      if (A(1,2) > A(2,2)) then
        V = A(1,:)
        A(1,:) = A(2,:)
        A(2,:) = V
      endif
      return
    endif

    NA=(N+1)/2
    NB=N-NA

    call MergeSort(A(:NA,:),NA,T)
    call MergeSort(A(NA+1:,:),NB,T)

    if (A(NA,2) > A(NA+1,2)) then
      T(1:NA,:)=A(1:NA,:)
      call merge(NA,NB,N,T(:NA,:),A(NA+1:N,:),A(:N,:))
    endif

  end subroutine mergeSort_indx_int


  !> Merge sort of reals
  subroutine merge_sort_real(array)

    !> vector to sort
    real(dp), intent(inout) :: array(:)

    real(dp), allocatable :: work(:)
    integer :: n

    n = size(array)

    allocate(work((n+1)/2))
    call mergeSort(array,n,work)

  end subroutine merge_sort_real


  !> Merge two arrays together in order onto a third
  subroutine merge_real(NA,NB,NC,A,B,C)

    !> first array of values
    integer, intent(in) :: NA

    !> elements in A
    integer, intent(in) :: NB

    !> second array of values
    integer, intent(in) :: NC

    !> elements in A
    real(dp), intent(in) :: A(NA)

    !> array to merge onto
    real(dp), intent(in) :: B(NB)

    !> elements in C
    real(dp), intent(inout) :: C(NC)

    integer :: I, J, K

    @:ASSERT((na+nb)==nc)

    I = 1; J = 1; K = 1;
    do while(I <= NA .and. J <= NB)
      if (A(I) <= B(J)) then
        C(K) = A(I)
        I = I+1
      else
        C(K) = B(J)
        J = J+1
      endif
      K = K + 1
    enddo
    do while (I <= NA)
      C(K) = A(I)
      I = I + 1
      K = K + 1
    enddo

  end subroutine merge_real


  !> Real merge sort

  recursive subroutine mergeSort_real(A,N,T)

    !> array to sort
    real(dp), intent(inout) :: A(:)

    !> number of elements in array
    integer, intent(in) :: N

    !> workspace of at least (N+1)/2 size
    real(dp), intent (out) :: T(:)

    integer :: NA, NB
    real(dp) :: V

    if (N < 2) return

    if (N == 2) then
      if (A(1) > A(2)) then
        V = A(1)
        A(1) = A(2)
        A(2) = V
      endif
      return
    endif

    NA=(N+1)/2
    NB=N-NA

    call MergeSort(A,NA,T)
    call MergeSort(A(NA+1:),NB,T)

    if (A(NA) > A(NA+1)) then
      T(1:NA)=A(1:NA)
      call merge(NA,NB,N,T,A(NA+1:),A)
    endif

  end subroutine mergeSort_real


  !> Merge sort of reals, using an array index instead of re-ordering

  subroutine merge_sort_indx_real(indx,array, tol)

    !> array of sorted order
    integer, intent(out) :: indx(:)

    !> vector to sort
    real(dp), intent(in) :: array(:)

    !> tolerance for comparisions
    real(dp), intent(in) :: tol

    real(dp), allocatable :: work(:,:), tmp(:,:)
    integer :: ii, n

    @:ASSERT(size(indx)==size(array))

    n = size(array)
    allocate(tmp(n,2))
    tmp(:,2) = array
    do ii = 1, n
      tmp(ii,1) = real(ii)
    end do
    allocate(work((n+1)/2,2))
    call mergeSort(tmp,n,work,tol)
    indx = nint(tmp(:,1))

  end subroutine merge_sort_indx_real


  !> Merge two arrays together in order onto a third, where first dimension of both is index for
  !> original order and also value
  subroutine merge_indx_real(NA,NB,NC,A,B,C, tol)

    !> first array of values
    integer, intent(in) :: NA

    !> elements in A
    integer, intent(in) :: NB

    !> second array of values
    integer, intent(in) :: NC

    !> elements in A
    real(dp), intent(in) :: A(NA,2)

    !> array to merge onto
    real(dp), intent(in) :: B(NB,2)

    !> elements in C
    real(dp), intent(inout) :: C(NC,2)

    real(dp), intent(in) :: tol

    integer :: I, J, K

    @:ASSERT((na+nb)==nc)

    I = 1; J = 1; K = 1;
    do while(I <= NA .and. J <= NB)
      if (A(I,2) <= B(J,2) .and. abs(A(I,2)-B(J,2)) > tol) then
        C(K,:) = A(I,:)
        I = I+1
      else
        C(K,:) = B(J,:)
        J = J+1
      endif
      K = K + 1
    enddo
    do while (I <= NA)
      C(K,:) = A(I,:)
      I = I + 1
      K = K + 1
    enddo

  end subroutine merge_indx_real


  !> Real merge sort, using an index
  recursive subroutine mergeSort_indx_real(A,N,T,tol)

    !> array to sort, first element of first dimension is an index array, second element is actual
    !> value
    real(dp), intent(inout) :: A(:,:)

    !> number of elements in array
    integer, intent(in) :: N

    !> workspace of at least (N+1)/2 size
    real(dp), intent (out) :: T(:,:)

    !> tolerance for comparisions
    real(dp), intent(in) :: tol

    integer :: NA, NB
    real(dp) :: V(2)

    @:ASSERT(size(A,dim=2) == 2)
    @:ASSERT(size(T,dim=2) == 2)

    if (N < 2) return

    if (N == 2) then
      if (A(1,2) > A(2,2) .and. abs(A(1,2) - A(2,2)) > tol) then
        V = A(1,:)
        A(1,:) = A(2,:)
        A(2,:) = V
      endif
      return
    endif

    NA=(N+1)/2
    NB=N-NA

    call MergeSort(A(:NA,:),NA,T,tol)
    call MergeSort(A(NA+1:,:),NB,T,tol)

    if (A(NA,2) > A(NA+1,2) .and. (A(NA,2) - A(NA+1,2)) > tol) then
      T(1:NA,:)=A(1:NA,:)
      call merge(NA,NB,N,T(:NA,:),A(NA+1:N,:),A(:N,:),tol)
    endif

  end subroutine mergeSort_indx_real


  !> Function to count number of unique elements in a sorted array of value greater than 0 and place
  !> them at the start of the array in order.
  !> To do: check that the elements are in sorted order, and generalise for
  !> decreasing order as well as increasing
  function unique_int(array, arraySize) result(nUnique)

    !> Array to make unique.
    integer, intent(inout) :: array(:)

    !> Constrains the effect of the subroutine on the first n elements, where n is the value for
    !> arraySize. (default: size(array))
    integer, intent(in), optional :: arraySize

    !> Number of unique elements.
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

  end function unique_int

end module sorting
