!--------------------------------------------------------------------------------------------------!
!  DFTB+: general package for performing fast atomistic simulations                                !
!  Copyright (C) 2006 - 2020  DFTB+ developers group                                               !
!                                                                                                  !
!  See the LICENSE file for terms of usage and distribution.                                       !
!--------------------------------------------------------------------------------------------------!

#:include 'common.fypp'

!> Contains Types and subroutine to build up and query a Slater-Koster table where the integrals are
!> specified on an equidistant grid.
module dftbp_slakoeqgrid
  use dftbp_assert
  use dftbp_accuracy
  use dftbp_interpolation
  use dftbp_message
  implicit none
  private

  public :: TSlakoEqGrid, init
  public :: getSKIntegrals, getNIntegrals, getCutoff
  public :: skEqGridOld, skEqGridNew


  !> Represents an equally spaced Slater-Koster grid
  type TSlakoEqGrid
    private
    integer :: nGrid
    integer :: nInteg
    real(dp) :: dist
    real(dp), allocatable :: skTab(:,:)
    integer :: skIntMethod
    logical :: tInit = .false.
  end type TSlakoEqGrid


  !> Initialises SlakoEqGrid.
  interface init
    module procedure SlakoEqGrid_init
  end interface init


  !> Returns the integrals for a given distance.
  interface getSKIntegrals
    module procedure SlakoEqGrid_getSKIntegrals
  end interface getSKIntegrals


  !> Returns the number of integrals the table contains
  interface getNIntegrals
    module procedure SlakoEqGrid_getNIntegrals
  end interface getNIntegrals


  !> Returns the cutoff of the interaction.
  interface getCutoff
    module procedure SlakoEqGrid_getCutoff
  end interface getCutoff

  ! Interpolation methods

  !> Historical method
  integer, parameter :: skEqGridOld = 1

  !> Current method
  integer, parameter :: skEqGridNew = 2

  ! Nr. of grid points to use for the polynomial interpolation

  !> Historical choice
  integer, parameter :: nInterOld_ = 3

  !> Present choice
  integer, parameter :: nInterNew_ = 8

  ! Nr. of grid points on the right of the interpolated point.

  ! For an odd number of intervals, the number of right points should be bigger than the number of
  ! left points, to remain compatible with the old code.


  !> value nRightInterOld: floor(real(nInterOld_, dp) / 2.0_dp + 0.6_dp)
  integer, parameter :: nRightInterOld_ = 2

  !> value nRightInterNew: floor(real(nInterNew_, dp) / 2.0_dp + 0.6_dp)
  integer, parameter :: nRightInterNew_ = 4


  !> Displacement for deriving interpolated polynomials
  real(dp), parameter :: deltaR_ = 1e-5_dp

contains


  !> Initialises SlakoEqGrid.
  subroutine SlakoEqGrid_init(this, dist, table, skIntMethod)

    !> SlakoEqGrid instance.
    type(TSlakoEqGrid), intent(out) :: this

    !> Distance between the grid points.
    real(dp), intent(in) :: dist

    !> Slater-Koster table (first entry belongs to first grid point)
    real(dp), intent(in) :: table(:,:)

    !> Method for the interpolation between the entries.
    integer, intent(in) :: skintMethod

    @:ASSERT(.not. this%tInit)
    @:ASSERT(dist >= 0.0_dp)
    @:ASSERT(skIntMethod == skEqGridOld .or. skIntMethod == skEqGridNew)

    this%dist = dist
    this%nGrid = size(table, dim=1)
    this%nInteg = size(table, dim=2)
    allocate(this%skTab(this%nGrid, this%nInteg))
    this%skTab(:,:) = table(:,:)
    this%skIntMethod = skIntMethod
    this%tInit = .true.

  end subroutine SlakoEqGrid_init


  !> Returns the integrals for a given distance.
  subroutine SlakoEqGrid_getSKIntegrals(this, sk, dist)

    !> SlakoEqGrid instance.
    type(TSlakoEqGrid), intent(in) :: this

    !> Contains the interpolated integrals on exit
    real(dp), intent(out) :: sk(:)

    !> Distance for which the integrals should be interpolated.
    real(dp), intent(in) :: dist

    @:ASSERT(this%tInit)
    @:ASSERT(size(sk) >= this%nInteg)
    @:ASSERT(dist >= 0.0_dp)

    if (this%skIntMethod == skEqGridOld) then
      call SlakoEqGrid_interOld_(this, sk, dist)
    else
      call SlakoEqGrid_interNew_(this, sk, dist)
    end if

  end subroutine SlakoEqGrid_getSKIntegrals


  !> Returns the number of intgrals the table contains
  function SlakoEqGrid_getNIntegrals(this) result(nInt)

    !> SlakoEqGrid instance.
    type(TSlakoEqGrid), intent(in) :: this

    !> Number of integrals.
    integer :: nInt

    nInt = this%nInteg

  end function SlakoEqGrid_getNIntegrals


  !> Returns the cutoff of the interaction.
  function SlakoEqGrid_getCutoff(this) result(cutoff)

    !>  SlakoEqGrid instance.
    type(TSlakoEqGrid), intent(in) :: this

    !> grid cutoff
    real(dp) :: cutoff

    cutoff = real(this%nGrid, dp) * this%dist
    if (this%skIntMethod == skEqGridOld) then
      cutoff = cutoff + distFudgeOld
    else
      cutoff = cutoff + distFudge
    end if

  end function SlakoEqGrid_getCutoff


  !> Inter- and extrapolation for SK-tables, new method.
  subroutine SlakoEqGrid_interNew_(this, dd, rr)

    !> SlakoEqGrid table on equiv. grid
    type(TSlakoEqGrid), intent(in) :: this

    !> Output table of interpolated values.
    real(dp), intent(out) :: dd(:)

    !> distance bewteen two atoms of interest
    real(dp), intent(in) :: rr

    real(dp) :: xa(nInterNew_), ya(nInterNew_), yb(this%nInteg,nInterNew_), y1, y1p, y1pp
    real(dp) :: incr, dr, rMax, y0(this%nInteg), y2(this%nInteg)
    integer :: leng, ind, iLast
    integer :: ii

    real(dp), parameter :: invdistFudge = -1.0_dp / distFudge

    leng = this%nGrid
    incr = this%dist
    rMax = real(leng, dp) * incr + distFudge
    ind = floor(rr / incr)

    !! Sanity check, if SK-table contains enough entries
    if (leng < nInterNew_ + 1) then
      call error("SlakoEqGrid: Not enough points in the SK-table for &
          &interpolation!")
    end if

    dd(:) = 0.0_dp
    if (rr >= rMax) then
      !! Beyond last grid point + distFudge => no interaction
      dd(:) = 0.0_dp
    elseif (ind < leng) then
      !! Closer to origin than last grid point => polynomial fit
      iLast = min(leng, ind + nRightInterNew_)
      iLast = max(iLast, nInterNew_)
      do ii = 1, nInterNew_
        xa(ii) = real(iLast - nInterNew_ + ii, dp) * incr
      end do
      yb = transpose(this%skTab(iLast-nInterNew_+1:iLast,:this%nInteg))
      dd(:this%nInteg) = polyInterUniform(xa, yb, rr)
    else
      !! Beyond the grid => extrapolation with polynomial of 5th order
      dr = rr - rMax
      iLast = leng
      do ii = 1, nInterNew_
        xa(ii) = real(iLast - nInterNew_ + ii, dp) * incr
      end do
      yb = transpose(this%skTab(iLast-nInterNew_+1:iLast,:this%nInteg))
      y0 = polyInterUniform(xa, yb, xa(nInterNew_) - deltaR_)
      y2 = polyInterUniform(xa, yb, xa(nInterNew_) + deltaR_)
      do ii = 1, this%nInteg
        ya(:) = this%skTab(iLast-nInterNew_+1:iLast, ii)
        y1 = ya(nInterNew_)
        y1p = (y2(ii) - y0(ii)) / (2.0_dp * deltaR_)
        y1pp = (y2(ii) + y0(ii) - 2.0_dp * y1) / (deltaR_ * deltaR_)
        dd(ii) = poly5ToZero(y1, y1p, y1pp, dr, -1.0_dp * distFudge, invDistFudge)
      end do
    end if

  end subroutine SlakoEqGrid_interNew_


  !> Inter- and extra-polation for SK-tables equivalent to the old DFTB code.
  subroutine SlakoEqGrid_interOld_(this, dd, rr)

    !> Data structure for SK interpolation
    type(TSlakoEqGrid), intent(in) :: this

    !> Output table of interpolated values.
    real(dp), intent(out) :: dd(:)

    !> distance bewteen two atoms of interest
    real(dp), intent(in) :: rr

    real(dp) :: xa(nInterOld_), yb(this%nInteg,nInterOld_),y0, y1, y2, y1p, y1pp
    real(dp) :: incr, dr
    integer :: leng, ind, mInd, iLast
    integer :: ii
    real(dp) :: r1, r2
    real(dp) :: invdistFudge

    leng = this%nGrid
    incr = this%dist
    mInd = leng + floor(distFudgeOld/incr)
    ind = floor(rr / incr)

    invdistFudge = -1.0_dp / (real(mInd - leng -1, dp) * incr)

    !! Sanity check, if SK-table contains enough entries
    if (leng < nInterOld_ + 1) then
      call error("skspar: Not enough points in the SK-table for interpolation!")
    end if

    dd(:) = 0.0_dp
    if (ind < leng-1) then
      !! Distance closer than penultimate grid point => polynomial fit
      iLast = min(leng, ind + nRightInterOld_)
      iLast = max(iLast, nInterOld_)
      do ii = 1, nInterOld_
        xa(ii) = real(iLast - nInterOld_ + ii, dp) * incr
      end do
      yb = transpose(this%skTab(iLast-nInterOld_+1:iLast,:this%nInteg))
      dd(:this%nInteg) = polyInterUniform(xa, yb, rr)
    elseif (ind < leng) then
      !! Distance between penultimate and last grid point => free cubic spline
      dr = rr - real(leng - 1, dp) * incr
      do ii = 1, this%nInteg
        y0 = this%skTab(leng-2, ii)
        y1 = this%skTab(leng-1, ii)
        y2 = this%skTab(leng, ii)
        y1p = (y2 - y0) / (2.0_dp * incr)
        y1pp = (y2 + y0 - 2.0_dp * y1) / incr**2
        call freeCubicSpline(y1, y1p, y1pp, incr, y2, dr, dd(ii))
      end do
    elseif (ind < mInd - 1) then
      !! Extrapolation
      dr = rr - real(mInd - 1, dp) * incr
      do ii = 1, this%nInteg
        y0 = this%skTab(leng-2, ii)
        y1 = this%skTab(leng-1, ii)
        y2 = this%skTab(leng, ii)
        r1 = (y2 - y0) / (2.0_dp * incr)
        r2 = (y2 + y0 - 2.0_dp * y1) / incr**2
        call freeCubicSpline(y1, r1, r2, incr, y2, incr, yp=y1p, ypp=y1pp)
        dd(ii) = poly5ToZero(y2, y1p, y1pp, dr,&
            & -1.0_dp * real(mInd - leng -1, dp)*incr, invdistFudge)
      end do
    else
      !! Dist. greater than tabulated sk range + distFudge => no interaction
      dd(:) = 0.0_dp
    end if

  end subroutine SlakoEqGrid_interOld_

end module dftbp_slakoeqgrid
