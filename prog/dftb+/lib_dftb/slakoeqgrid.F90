!--------------------------------------------------------------------------------------------------!
!  DFTB+: general package for performing fast atomistic simulations                                !
!  Copyright (C) 2017  DFTB+ developers group                                                      !
!                                                                                                  !
!  See the LICENSE file for terms of usage and distribution.                                       !
!--------------------------------------------------------------------------------------------------!

#:include 'common.fypp'

!!* Contains Types and subroutine to build up and query a Slater-Koster table
!!* where the integrals are specified on an equidistant grid.
module slakoeqgrid
  use assert
  use accuracy
  use interpolation
  use message
  implicit none
  private

  public :: OSlakoEqGrid, init
  public :: getSKIntegrals, getNIntegrals, getCutoff
  public :: skEqGridOld, skEqGridNew

  !!* Represents an equally spaced Slater-Koster grid
  type OSlakoEqGrid
    private
    integer :: nGrid
    integer :: nInteg
    real(dp) :: dist
    real(dp), allocatable :: skTab(:,:)
    integer :: skIntMethod
    logical :: tInit = .false.
  end type OSlakoEqGrid

  !!* Initialises SlakoEqGrid.
  interface init
    module procedure SlakoEqGrid_init
  end interface

  !!* Returns the integrals for a given distance.
  interface getSKIntegrals
    module procedure SlakoEqGrid_getSKIntegrals
  end interface

  !!* Returns the number of integrals the table contains
  interface getNIntegrals
    module procedure SlakoEqGrid_getNIntegrals
  end interface

  !!* Returns the cutoff of the interaction.
  interface getCutoff
    module procedure SlakoEqGrid_getCutoff
  end interface


  !! Interpolation methods
  integer, parameter :: skEqGridOld = 1
  integer, parameter :: skEqGridNew = 2


  !! Nr. of grid points to use for the polynomial interpolation
  integer, parameter :: nInterOld_ = 3
  integer, parameter :: nInterNew_ = 8

  !! Nr. of grid points on the right of the interpolated point
  !! For odd nr. of intervals, nr. of right points should be bigger than
  !! nr. of left points, to remain compatible with the old code.

  !! value nRightInterOld: floor(real(nInterOld_, dp) / 2.0_dp + 0.6_dp)
  integer, parameter :: nRightInterOld_ = 2
  !! value nRightInterNew: floor(real(nInterNew_, dp) / 2.0_dp + 0.6_dp)
  integer, parameter :: nRightInterNew_ = 4

  !! Displacement for deriving interpolated polynomials
  real(dp), parameter :: deltaR_ = 1e-5_dp



contains

  !!* Initialises SlakoEqGrid.
  !!* @param self SlakoEqGrid instance.
  !!* @param dist Distance between the grid points.
  !!* @param table Slater-Koster table (first entry belongs to first grid point)
  !!* @param skIntMethod Method for the interpolation between the entries.
  subroutine SlakoEqGrid_init(self, dist, table, skIntMethod)
    type(OSlakoEqGrid), intent(out) :: self
    real(dp), intent(in) :: dist
    real(dp), intent(in) :: table(:,:)
    integer, intent(in) :: skintMethod

    @:ASSERT(.not. self%tInit)
    @:ASSERT(dist >= 0.0_dp)
    @:ASSERT(skIntMethod == skEqGridOld .or. skIntMethod == skEqGridNew)

    self%dist = dist
    self%nGrid = size(table, dim=1)
    self%nInteg = size(table, dim=2)
    allocate(self%skTab(self%nGrid, self%nInteg))
    self%skTab(:,:) = table(:,:)
    self%skIntMethod = skIntMethod
    self%tInit = .true.

  end subroutine SlakoEqGrid_init



  !!* Returns the integrals for a given distance.
  !!* @param self SlakoEqGrid instance.
  !!* @param sk Contains the interpolated integrals on exit
  !!* @param dist Distance for which the integrals should be interpolated.
  subroutine SlakoEqGrid_getSKIntegrals(self, sk, dist)
    type(OSlakoEqGrid), intent(in) :: self
    real(dp), intent(out) :: sk(:)
    real(dp), intent(in) :: dist

    @:ASSERT(self%tInit)
    @:ASSERT(size(sk) >= self%nInteg)
    @:ASSERT(dist >= 0.0_dp)

    if (self%skIntMethod == skEqGridOld) then
      call SlakoEqGrid_interOld_(self, sk, dist)
    else
      call SlakoEqGrid_interNew_(self, sk, dist)
    end if

  end subroutine SlakoEqGrid_getSKIntegrals



  !!* Returns the number of intgrals the table contains
  !!* @param self SlakoEqGrid instance.
  !!* @return Number of integrals.
  function SlakoEqGrid_getNIntegrals(self) result(nInt)
    type(OSlakoEqGrid), intent(in) :: self
    integer :: nInt

    nInt = self%nInteg

  end function SlakoEqGrid_getNIntegrals



  !!* Returns the cutoff of the interaction.
  !!* @param self  SlakoEqGrid instance.
  !!* @return Cutoff.
  function SlakoEqGrid_getCutoff(self) result(cutoff)
    type(OSlakoEqGrid), intent(in) :: self
    real(dp) :: cutoff

    cutoff = real(self%nGrid, dp) * self%dist
    if (self%skIntMethod == skEqGridOld) then
      cutoff = cutoff + distFudgeOld
    else
      cutoff = cutoff + distFudge
    end if

  end function SlakoEqGrid_getCutoff



  !!* Inter- and extrapolation for SK-tables, new method.
  !!* @param self SlakoEqGrid table on equiv. grid
  !!* @param dd Output table of interpolated values.
  !!* @param rr distance bewteen two atoms of interest
  subroutine SlakoEqGrid_interNew_(self, dd, rr)
    type(OSlakoEqGrid), intent(in) :: self
    real(dp), intent(out) :: dd(:)
    real(dp), intent(in) :: rr

    real(dp) :: xa(nInterNew_), ya(nInterNew_), y0, y1, y2, y1p, y1pp
    real(dp) :: incr, dr, rMax
    integer :: leng, ind, iLast
    integer :: ii


    leng = self%nGrid
    incr = self%dist
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
      do ii = 1, self%nInteg
        ya(:) = self%skTab(iLast-nInterNew_+1:iLast, ii)
        dd(ii) = polyInter(xa, ya, rr)
      end do
    else
      !! Beyond the grid => extrapolation with polynomial of 5th order
      dr = rr - rMax
      iLast = leng
      do ii = 1, nInterNew_
        xa(ii) = real(iLast - nInterNew_ + ii, dp) * incr
      end do
      do ii = 1, self%nInteg
        ya(:) = self%skTab(iLast-nInterNew_+1:iLast, ii)
        y1 = ya(nInterNew_)
        y0 = polyInter(xa, ya, xa(nInterNew_) - deltaR_)
        y2 = polyInter(xa, ya, xa(nInterNew_) + deltaR_)
        y1p = (y2 - y0) / (2.0_dp * deltaR_)
        y1pp = (y2 + y0 - 2.0_dp * y1) / (deltaR_ * deltaR_)
        dd(ii) = poly5ToZero(y1, y1p, y1pp, dr, -1.0_dp * distFudge)
      end do
    end if

  end subroutine SlakoEqGrid_interNew_



  !!* Inter- and extrapolation for SK-tables like in the old DFTB code.
  !!* @param iSp1 Atom species for the first atom of interest.
  !!* @param iSp2 Atom species for the second atom of interest.
  !!* @param rr distance bewteen two atoms of interest
  !!* @param dd Output table of interpolated values.
  !!* @param skIncr Increment between elements in the SK table for species pairs
  !!* @param skLen Number of elements in the SK table for species pairs
  !!* @param mAngSpecies Maximum values of l for the species of atoms
  !!* @param skTab Input SK table
  !!* @param llmIndx index array to convert from l1,l2,m value into column in
  !!*   the SK table
  !!* @param off ofset for counting in llm interpolation =0 for all table
  !!*   =1 for heteronuclear second call - then only l1!=l2 case is calculated
  !!* @desc If the point is between the 1st
  subroutine SlakoEqGrid_interOld_(self, dd, rr)
    type(OSlakoEqGrid), intent(in) :: self
    real(dp), intent(out) :: dd(:)
    real(dp), intent(in) :: rr

    real(dp) :: xa(nInterOld_), ya(nInterOld_), y0, y1, y2, y1p, y1pp
    real(dp) :: incr, dr
    integer :: leng, ind, mInd, iLast
    integer :: ii
    real(dp) :: r1, r2

    leng = self%nGrid
    incr = self%dist
    mInd = leng + floor(distFudgeOld/incr)
    ind = floor(rr / incr)

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
      do ii = 1, self%nInteg
        ya(:) = self%skTab(iLast-nInterOld_+1:iLast, ii)
        dd(ii) = polyInter(xa, ya, rr)
      end do
    elseif (ind < leng) then
      !! Distance between penultimate and last grid point => free cubic spline
      dr = rr - real(leng - 1, dp) * incr
      do ii = 1, self%nInteg
        y0 = self%skTab(leng-2, ii)
        y1 = self%skTab(leng-1, ii)
        y2 = self%skTab(leng, ii)
        y1p = (y2 - y0) / (2.0_dp * incr)
        y1pp = (y2 + y0 - 2.0_dp * y1) / incr**2
        call freeCubicSpline(y1, y1p, y1pp, incr, y2, dr, dd(ii))
      end do
    elseif (ind < mInd - 1) then
      !! Extrapolation
      dr = rr - real(mInd - 1, dp) * incr
      do ii = 1, self%nInteg
        y0 = self%skTab(leng-2, ii)
        y1 = self%skTab(leng-1, ii)
        y2 = self%skTab(leng, ii)
        r1 = (y2 - y0) / (2.0_dp * incr)
        r2 = (y2 + y0 - 2.0_dp * y1) / incr**2
        call freeCubicSpline(y1, r1, r2, incr, y2, incr, yp=y1p, ypp=y1pp)
        dd(ii) = poly5ToZero(y2, y1p, y1pp, dr, &
            &-1.0_dp * real(mInd - leng -1, dp)*incr)
      end do
    else
      !! Dist. greater than tabulated sk range + distFudge => no interaction
      dd(:) = 0.0_dp
    end if

  end subroutine SlakoEqGrid_interOld_


end module slakoeqgrid
