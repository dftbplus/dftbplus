!--------------------------------------------------------------------------------------------------!
!  DFTB+: general package for performing fast atomistic simulations                                !
!  Copyright (C) 2006 - 2020  DFTB+ developers group                                               !
!                                                                                                  !
!  See the LICENSE file for terms of usage and distribution.                                       !
!--------------------------------------------------------------------------------------------------!

#:include 'common.fypp'

!> Implements an interator over lattice points up to a certain cutoff
module dftbp_latpointiter
  use dftbp_accuracy, only : dp
  use dftbp_assert
  implicit none
  private

  public :: TLatPointIter, TLatPointIter_init


  !> Lattice point iterator
  type :: TLatPointIter
    private

    !> Square of the cutoff
    real(dp) :: cutoff2

    !> Lattice vectors
    real(dp) :: latVecs(3, 3)

    !> Index ranges necessary to obtain all lattcie points within given cutoff
    integer :: ranges(2, 3)

    !> Indices of the current lattice point
    integer :: curPoint(3)

    !> Increase upper limit of necessary index ranges by this amount
    integer :: posExt = 0

    !> Decrease lower limit of necessary index ranges by this amount
    integer :: negExt = 0

    !> Return all lattice vector in the index ranges, not just those shorter than cutoff
    logical :: tAll = .true.

    !> Whether lattice vector (0, 0, 0) should be excluded when iterationg
    logical :: tExcludeOrig = .true.

    !> Whether to exclude points connected by inversion
    logical :: tExcludeInv = .false.

    !> Whether iterator has returned last point
    logical :: tFinished = .false.
  contains
    procedure :: getNextPoint => TLatPointIter_getNextPoint
    procedure :: getAllPoints => TLatPointIter_getAllPoints
  end type TLatPointIter

contains


  !> Initialises the lattice point iterator.
  subroutine TLatPointIter_init(this, latVecs, invLatVecs, cutoff, negExtension, posExtension,&
      & onlyInside, reduceByInversion, excludeOrigin)

    !> Instance
    type(TLatPointIter), intent(out) :: this

    !> Lattice vectors
    real(dp), intent(in) :: latVecs(:,:)

    !> Inverse lattice vectors
    real(dp), intent(in) :: invLatVecs(:,:)

    !> Cutoff radius for lattice points
    real(dp), intent(in) :: cutoff

    !> Same as posExtension for negative lattice vectors
    integer, intent(in), optional :: negExtension

    !> Extend the set along the positive lattice vectors with that many additional lattice vectors.
    integer, intent(in), optional :: posExtension

    !> Include only lattice points being inside the cutoff radius
    logical, intent(in), optional :: onlyInside

    !> whether to include time reversal symmetry when generating k-points
    logical,  intent(in), optional :: reduceByInversion

    !> whether to exclude the (0,0,0) point
    logical,  intent(in), optional :: excludeOrigin

    integer :: ranges(2, 3)

    if (present(negExtension)) then
      this%negExt = negExtension
    end if
    if (present(posExtension)) then
      this%posExt = posExtension
    end if
    if (present(onlyInside)) then
      this%tAll = .not. onlyInside
    end if
    if (present(reduceByInversion)) then
      this%tExcludeInv = reduceByInversion
    end if
    if (present(excludeOrigin)) then
      this%tExcludeOrig = excludeOrigin
    end if

    this%latVecs(:,:) = latVecs
    call getRanges(cutoff, invLatVecs, this%posExt, this%negExt, ranges)
    this%ranges(1,:) = minval(ranges, dim=1)
    this%ranges(2,:) = maxval(ranges, dim=1)
    this%curPoint(:) = this%ranges(1,:)
    this%cutoff2 = cutoff**2

  end subroutine TLatPointIter_init


  !> Delivers the next lattice point
  subroutine TLatPointIter_getNextPoint(this, latticePoint, tFinished)

    !> Instance
    class(TLatPointIter), intent(inout) :: this

    !> Next lattice point
    real(dp), intent(out) :: latticePoint(:)

    !> Whether the returned point was the last one.
    logical, intent(out) :: tFinished

    real(dp) :: rr(3)
    integer :: curPoint(3)

    curPoint(:) = this%curPoint
    do
      tFinished = (curPoint(1) > this%ranges(2, 1))
      if (tFinished) then
        exit
      end if
      if (this%tExcludeInv) then
        if (curPoint(1) < 0) then
          curPoint(1) = 0
        end if
        if (curPoint(2) < 0 .and. curPoint(1) == 0) then
          curPoint(2) = 0
        end if
        if (curPoint(3) < 0 .and. curPoint(1) == 0 .and. curPoint(2) == 0) then
          curPoint(3) = 0
        end if
      end if
      if (this%tExcludeOrig .and. all(curPoint == 0)) then
        call getNextPoint(this%ranges, curPoint)
        cycle
      end if
      latticePoint(:) = real(curPoint, dp)
      if (this%tAll) then
        exit
      end if
      rr(:) = latticePoint(1) * this%latVecs(:,1) + latticePoint(2) * this%latVecs(:,2)&
          & + latticePoint(3) * this%latVecs(:,3)
      if (sum(rr**2) <= this%cutoff2) then
        exit
      end if
      call getNextPoint(this%ranges, curPoint)
    end do

    if (.not. tFinished) then
      ! generate image for next call
      call getNextPoint(this%ranges, curPoint)
    end if
    this%curPoint(:) = curPoint

  end subroutine TLatPointIter_getNextPoint


  !> Returns all lattice points within a cutoff at the same time
  subroutine TLatPointIter_getAllPoints(this, latticePoints)

    !> Instance
    class(TLatPointIter), intent(inout) :: this

    !> Lattice points
    real(dp), allocatable, intent(out) :: latticePoints(:,:)

    real(dp), allocatable :: tmpLatPoints(:,:)
    integer :: maxLatPoints, iLatPoint
    logical :: tFinished

    maxLatPoints = product(this%ranges(2,:) - this%ranges(1,:) + 1)
    allocate(tmpLatPoints(3, maxLatPoints + 1))
    do iLatPoint = 1, maxLatPoints + 1
      call this%getNextPoint(tmpLatPoints(:,iLatPoint), tFinished)
      if (tFinished) then
        exit
      end if
    end do
    latticePoints = tmpLatPoints(:, 1:iLatPoint - 1)

  end subroutine TLatPointIter_getAllPoints


  !> Helper function to increase a tuple of 3 indices by one
  subroutine getNextPoint(ranges, point)

    !> Lower and upper bounds for the lattice point indices. Shape: (2, 3)
    integer, intent(in) :: ranges(:,:)

    !> Current lattice point, next one on exit.
    integer, intent(inout) :: point(3)
    
    point(3) = point(3) + 1
    if (point(3) > ranges(2, 3)) then
      point(3) = ranges(1, 3)
      point(2) = point(2) + 1
      if (point(2) > ranges(2, 2)) then
        point(2) = ranges(1, 2)
        point(1) = point(1) + 1
      end if
    end if

  end subroutine getNextPoint


  !> Calculate the range of images of the central cell that interact
  subroutine getRanges(dist, recVec2p, posExt, negExt, ranges)

    !> distance of interaction
    real(dp), intent(in) :: dist

    !> reciprocal lattice vector
    real(dp), intent(in) :: recVec2p(:,:)

    !> Extend the set along the positive lattice vectors with that many additional lattice vectors.
    integer, intent(in) :: posExt

    !> Same as posExtension for negative lattice vectors
    integer, intent(in) :: negExt

    !> Array of the two extremal points
    integer, intent(out) :: ranges(:,:)

    integer :: ii, iTmp

    @:ASSERT(dist >= 0.0_dp)

    do ii = 1, 3
      iTmp = floor(dist * sqrt(sum(recVec2p(:, ii)**2)))
      ranges(1, ii) = -(iTmp + negExt)
      ranges(2, ii) = iTmp + posExt
    end do

  end subroutine getRanges

end module dftbp_latpointiter
