!--------------------------------------------------------------------------------------------------!
!  DFTB+: general package for performing fast atomistic simulations                                !
!  Copyright (C) 2006 - 2022  DFTB+ developers group                                               !
!                                                                                                  !
!  See the LICENSE file for terms of usage and distribution.                                       !
!--------------------------------------------------------------------------------------------------!

!> Simple geometry calculations
module dftbp_math_simplegeometry
  use dftbp_common_accuracy, only : dp
  use dftbp_math_sorting, only : merge_sort
  implicit none

  private
  public :: dist2Segment, radialSort, uniqueCoords

contains

  !> Square of the distance of a point to a line segment between points p1 and p2
  pure function dist2Segment(p1, p2, q) result(dist)

    !> Start of segment
    real(dp), intent(in) :: p1(:)

    !> End of segment
    real(dp), intent(in) :: p2(:)

    !> Point to measure distance from
    real(dp), intent(in) :: q(:)

    !> Resulting distance squared to segment
    real(dp) :: dist

    real(dp) :: u(size(p1)), v(size(p1)), p(size(p1)), d

    u(:) = p2 - p1
    v(:) = q - p1

    d = dot_product(u, v) / dot_product(u,u) ! projection onto line segment
    if (d < 0.0_dp) then ! projection is before start of segment
      dist = dot_product(q-p1,q-p1)
    else if (d > 1.0_dp) then ! projection is after end of seqment
      dist = dot_product(q-p2,q-p2)
    else ! projection is inside segment
      p(:) = p1 + d * u ! nearest point inside segment
      dist = dot_product(q-p,q-p)
    end if

  end function dist2Segment


  !> Sort a set of points with respect to the distance from the origin (using a stable sort), and
  !> return an index array, where the first nPoints are unique
  subroutine radialSort(points, indx, nPoints)

    !> Points to sort
    real(dp), intent(in) :: points(:,:)

    integer, intent(out) :: indx(:)

    !> Count of unique points
    integer, intent(out) :: nPoints

    real(dp), allocatable :: dist2(:)
    integer, allocatable :: indxR(:), iUniq(:)
    integer :: iStart, iEnd

    dist2 = sum(points**2, dim=1)
    allocate(indxR(size(points, dim=2)))
    call merge_sort(indxR, dist2, epsilon(0.0_dp))
    dist2(:) = dist2(indxR)
    nPoints = 0
    iStart = 1
    do iEnd = 1, size(points, dim=2)
      if (abs(dist2(iStart)-dist2(iEnd))<=epsilon(0.0_dp)) then
        cycle
      end if
      call uniqueCoords(points(:,indxR(iStart:iEnd-1)) ,iUniq)
      iUniq = iUniq + iStart - 1
      indx(nPoints+1:nPoints+size(iUniq)) = indxR(iUniq)
      nPoints = nPoints + size(iUniq)
      deallocate(iUniq)
      iStart = iEnd
    end do
    call uniqueCoords(points(:,indxR(iStart:iEnd-1)) ,iUniq)
    iUniq = iUniq + iStart - 1
    indx(nPoints+1:nPoints+size(iUniq)) = indxR(iUniq)
    nPoints = nPoints + size(iUniq)

  end subroutine radialSort


  !> Find first instances of unique coordinates in a set, returning an index to them
  subroutine uniqueCoords(coords, iUnique)

    !> Coordinates to check
    real(dp) :: coords(:,:)

    !> Index for unique atoms
    integer, allocatable, intent(out) :: iUnique(:)

    integer, allocatable :: work(:)
    integer :: nUnique, ii, jj, n

    n = size(coords,dim=2)
    allocate(work(n), source = 0)
    nUnique = 0
    do ii = 1, n
      if (work(ii) < 0) then
        cycle
      else
        nUnique = nUnique + 1
        work(nUnique) = ii
      end if
      do jj = ii+1, n
        if (all(abs(coords(:,ii) - coords(:,jj)) < 1024_dp*epsilon(0.0_dp))) then
          work(jj) = -1
        end if
      end do
    end do
    iUnique = work(:nUnique)

  end subroutine uniqueCoords

end module dftbp_math_simplegeometry
