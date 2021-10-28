!--------------------------------------------------------------------------------------------------!
!  DFTB+: general package for performing fast atomistic simulations                                !
!  Copyright (C) 2006 - 2021  DFTB+ developers group                                               !
!                                                                                                  !
!  See the LICENSE file for terms of usage and distribution.                                       !
!--------------------------------------------------------------------------------------------------!

#:include 'common.fypp'

!> Module for creating and modifying three-dimensional, equidistant grids.
module waveplot_grids

  use dftbp_common_accuracy, only : dp
  use dftbp_math_simplealgebra, only : invert33
  use dftbp_io_message, only : error

  use waveplot_interp, only : linearInterpolation, trilinearInterpolation, splineInterpolation,&
      & locateByBisection
  use waveplot_parallel, only : getStartAndEndIndices
  use waveplot_slater, only: RealTessY

#:if WITH_OMP
  use omp_lib
#:endif

  implicit none
  private

  public :: TGrid, TGrid_init
  public :: TGridData, TGridData_init
  public :: TRealTessY, TRealTessY_init
  public :: subgridsToGlobalGrid


  type :: TGrid

    !> cartesian coordinates of grid origin
    real(dp) :: origin(3)

    !> array of column vectors spanning dense grid basis
    real(dp) :: basis(3,3)

    !> inverted basis
    real(dp) :: invBasis(3,3)

    !> lower range bound for every axis
    integer :: lowerBounds(3)

    !> upper range bound for every axis
    integer :: upperBounds(3)

    !> number of grid points along each dimension
    integer :: nPoints(3)

  contains

    procedure :: setOrigin => TGrid_setOrigin
    procedure :: getCorners => TGrid_getCorners
    procedure :: cartesianToGridcoord => TGrid_cartesianToGridcoord
    procedure :: cartesianToRealGridcoord => TGrid_cartesianToRealGridcoord
    procedure :: gridcoordToCartesian => TGrid_gridcoordToCartesian
    procedure :: realGridcoordToCartesian => TGrid_realGridcoordToCartesian
    procedure :: getSubgridRanges => TGrid_getSubgridRanges
    procedure :: hasGridcoords => TGrid_hasGridcoords
    procedure :: hasSubgrid => TGrid_hasSubgrid
    procedure :: getIntersecGridPoints => TGrid_getIntersecGridPoints

  end type TGrid


  type :: TRealTessY

    !> representation of parallelepiped grid
    class(TGrid), allocatable :: grid

    !> tabulated real tesseral spherical harmonics
    real(dp), allocatable :: data(:,:,:,:)

    !> possible combinations of angular momenta l and magnetic quantum numbers m
    integer, allocatable :: lmCombs(:,:)

    !> maximum angular momentum of the system
    integer :: maxAng

  contains

    procedure :: tabulateRealTessY => TRealTessY_tabulateRealTessY

  end type TRealTessY


  type :: TGridData

    !> representation of parallelepiped grid
    class(TGrid), allocatable :: grid

    !> tabulated real tesseral spherical harmonics for all l,m-combinations
    class(TRealTessY), allocatable :: rty

    !> volumetric data pointer
    real(dp), pointer :: data(:,:,:)

    procedure(add), pointer, pass(this) :: add => null()
    procedure(tabulateBasis), pointer, pass(this) :: tabulateBasis => null()

  contains

    procedure :: setGridOrigin => TGridData_setGridOrigin
    procedure :: setRwInterpolationType => TGridData_setRwInterpolationType
    procedure :: setGridInterpolationType => TGridData_setGridInterpolationType
    procedure :: getValue => TGridData_getValue
    procedure :: getInterpolatedValueGc => TGridData_getInterpolatedValueGc
    procedure :: getSubcubeCornsAndVals => TGridData_getSubcubeCornsAndVals

  end type TGridData


  abstract interface

    subroutine add(this, griddata, eigvec, position, square)

      import :: TGridData, dp
      implicit none

      !> representation of volumetric grid and data
      class(TGridData), intent(inout) :: this

      !> volumetric data to add to current instance
      class(TGridData), intent(in) :: griddata

      !> eigenvector of current basis
      real(dp), intent(in) :: eigvec

      !> position where to place basis grid onto total grid, shape: [3]
      real(dp), intent(in) :: position(:)

      !> True, if data should be squared, i.e. for densities
      logical, intent(in), optional :: square

    end subroutine add


    subroutine tabulateBasis(this, rwf, rty, ll, mm)

      import :: TGridData, TRealTessY, dp
      implicit none

      !> representation of volumetric grid and data
      class(TGridData), intent(inout) :: this

      !> tabulated radial wavefunction, shape: [nDistances, 3]
      real(dp), intent(in) :: rwf(:,:)

      !> tabulated real tesseral spherical harmonics for all l,m-combinations
      class(TRealTessY), intent(in) :: rty

      !> angular momentum of the species
      integer, intent(in) :: ll

      !> magnetic quantum number
      integer, intent(in) :: mm

    end subroutine tabulateBasis

  end interface


  interface subgridsToGlobalGrid
    module procedure :: subgridsToUncachedGlobalGrid, subgridsToCachedGlobalGrids
  end interface subgridsToGlobalGrid


contains

  subroutine TGrid_init(this, origin, basis, nPoints)

    !> representation of parallelepiped grid
    class(TGrid), intent(inout) :: this

    !> cartesian coordinates of grid origin
    real(dp), intent(in) :: origin(:)

    !> array of column vectors spanning dense grid basis
    real(dp), intent(in) :: basis(:,:)

    !> number of points along each axis
    integer, intent(in) :: nPoints(:)

    this%origin = origin
    this%nPoints = nPoints
    this%basis = basis

    this%lowerBounds = [0, 0, 0]
    this%upperBounds = nPoints - 1

    call invert33(this%invBasis, this%basis)

  end subroutine TGrid_init


  subroutine TGrid_setOrigin(this, origin)

    !> representation of parallelepiped grid
    class(TGrid), intent(inout) :: this

    !> cartesian coordinates of grid origin
    real(dp), intent(in) :: origin(:)

    this%origin = origin

  end subroutine TGrid_setOrigin


  subroutine TGrid_gridcoordToCartesian(this, gridcoords, cartcoords)

    !> representation of parallelepiped grid
    class(TGrid), intent(in) :: this

    !> integer grid coordinates
    integer, intent(in) :: gridcoords(:,:)

    !> corresponding cartesian grid coordinates
    real(dp), intent(out), allocatable :: cartcoords(:,:)

    !> auxiliary variable
    integer :: ii

    allocate(cartcoords(3, size(gridcoords, dim=2)))

    do ii = 1, size(gridcoords, dim=2)
      cartcoords(:, ii) = matmul(this%basis, gridcoords(:, ii)) + this%origin
    end do

  end subroutine TGrid_gridcoordToCartesian


  subroutine TGrid_realGridcoordToCartesian(this, gridcoords, cartcoords)

    !> representation of parallelepiped grid
    class(TGrid), intent(in) :: this

    !> real grid coordinates
    real(dp), intent(in) :: gridcoords(:,:)

    !> corresponding cartesian grid coordinates
    real(dp), intent(out), allocatable :: cartcoords(:,:)

    !> auxiliary variable
    integer :: ii

    allocate(cartcoords(3, size(gridcoords, dim=2)))

    do ii = 1, size(gridcoords, dim=2)
      cartcoords(:, ii) = matmul(this%basis, gridcoords(:, ii)) + this%origin
    end do

  end subroutine TGrid_realGridcoordToCartesian


  subroutine TGrid_cartesianToGridcoord(this, cartcoords, gridcoords)

    !> representation of parallelepiped grid
    class(TGrid), intent(in) :: this

    !> cartesian coordinates
    real(dp), intent(in) :: cartcoords(:,:)

    !> corresponding integer grid coordinates
    integer, intent(out), allocatable :: gridcoords(:,:)

    !> auxiliary variable
    integer :: ii

    allocate(gridcoords(size(cartcoords, dim=1), size(cartcoords, dim=2)))

    do ii = 1, size(cartcoords, dim=2)
      gridcoords(:, ii) = floor(matmul(cartcoords(:, ii) - this%origin, this%invBasis))
    end do

  end subroutine TGrid_cartesianToGridcoord


  subroutine TGrid_cartesianToRealGridcoord(this, cartcoords, gridcoords)

    !> representation of parallelepiped grid
    class(TGrid), intent(in) :: this

    !> cartesian coordinates
    real(dp), intent(in) :: cartcoords(:,:)

    !> corresponding real grid coordinates
    real(dp), intent(out), allocatable :: gridcoords(:,:)

    !> auxiliary variable
    integer :: ii

    allocate(gridcoords(size(cartcoords, dim=1), size(cartcoords, dim=2)))

    do ii = 1, size(cartcoords, dim=2)
      gridcoords(:, ii) = matmul(cartcoords(:, ii) - this%origin, this%invBasis)
    end do

  end subroutine TGrid_cartesianToRealGridcoord


  subroutine TGrid_getCorners(this, corncoords)

    !> representation of parallelepiped grid
    class(TGrid), intent(in) :: this

    !> integer grid coordinates of the 2**3=8 corners of the parallelepipedon
    integer, intent(out) :: corncoords(3,8)

    !> indices of corners
    integer, allocatable :: cornerInds(:,:)

    !> grid integer ranges along each axis
    integer :: ranges(2,3)

    !> auxiliary variables
    integer :: ii, jj, kk, idx, tmp(3)

    cornerInds = reshape([&
        & 1, 1, 1,&
        & 1, 1, 2,&
        & 1, 2, 1,&
        & 1, 2, 2,&
        & 2, 1, 1,&
        & 2, 1, 2,&
        & 2, 2, 1,&
        & 2, 2, 2&
        &], [3, 8])

    do ii = 1, 3
      ranges(:, ii) = [this%lowerBounds(ii), this%upperBounds(ii)]
    end do

    do jj = 1, size(cornerInds, dim=2)
      do kk = 1, size(cornerInds, dim=1)
        idx = cornerInds(kk, jj)
        tmp(kk) = ranges(idx, kk)
      end do
      corncoords(:, jj) = tmp(:)
    end do

  end subroutine TGrid_getCorners


  subroutine getGridcoords(grid, gridcoords)

    !> representation of parallelepiped grid
    type(TGrid), intent(in) :: grid

    !> integer grid coordinates
    integer, intent(out), allocatable :: gridcoords(:,:)

    !> auxiliary variables
    integer :: ii, jj, kk, ind

    allocate(gridcoords(3,&
        & (grid%upperBounds(1) - grid%lowerBounds(1) + 1) *&
        & (grid%upperBounds(2) - grid%lowerBounds(2) + 1) *&
        & (grid%upperBounds(3) - grid%lowerBounds(3) + 1)))

    ind = 1

    do kk = grid%lowerBounds(3), grid%upperBounds(3)
      do jj = grid%lowerBounds(2), grid%upperBounds(2)
        do ii = grid%lowerBounds(1), grid%upperBounds(1)
          gridcoords(1, ind) = ii
          gridcoords(2, ind) = jj
          gridcoords(3, ind) = kk
          ind = ind + 1
        end do
      end do
    end do

  end subroutine getGridcoords


  subroutine TGrid_getSubgridRanges(this, subgrid, position, luGc, luSubGc)

    !> representation of parallelepiped grid
    class(TGrid), intent(in) :: this

    !> subgrid to calculate intersection with
    class(TGrid), intent(in) :: subgrid

    !> position where to place subgrid onto grid, expected shape: [3]
    real(dp), intent(in) :: position(:)

    !> lower and upper bounds of intersection in integer grid coordinates
    integer, intent(out), allocatable :: luGc(:,:)

    !> lower and upper bounds of intersection in real subgrid coordinates
    real(dp), intent(out), allocatable :: luSubGc(:,:)

    !> lower and upper bounds of intersection in real grid coordinates
    real(dp), allocatable :: realLuGc(:,:)

    !> lower and upper bounds of intersection in real subgrid coordinates
    real(dp), allocatable :: realLuSubGc(:,:)

    !> lower and upper bounds of intersection in cartesian coordinates
    real(dp) :: luCc(3, 2)

    !> auxiliary array to containing cartesian coordinates of processed intersection
    real(dp), allocatable :: roundedLuCc(:,:)

    !> lower and upper cartesian bounds of grid and (partly) contained subgrid
    real(dp), allocatable :: lCcTotGrid(:,:), uCcTotGrid(:,:), lCcSubGrid(:,:), uCcSubGrid(:,:)

    !> auxiliary variable
    integer :: ii

    call this%gridcoordToCartesian(reshape(this%lowerBounds, [3, 1]), lCcTotGrid)
    call this%gridcoordToCartesian(reshape(this%upperBounds + 1, [3, 1]), uCcTotGrid)
    call subgrid%gridcoordToCartesian(reshape(subgrid%lowerBounds, [3, 1]), lCcSubGrid)
    call subgrid%gridcoordToCartesian(reshape(subgrid%upperBounds, [3, 1]), uCcSubGrid)

    luCc(:, 1) = max(lCcSubGrid(:, 1) + position, lCcTotGrid(:, 1))
    luCc(:, 2) = min(uCcSubGrid(:, 1) + position, uCcTotGrid(:, 1))

    call this%cartesianToGridcoord(luCc, luGc)

    luGc(:, 1) = max(luGc(:, 1), this%lowerBounds)
    luGc(:, 1) = min(luGc(:, 1), this%upperBounds)
    luGc(:, 2) = max(luGc(:, 2), this%lowerBounds)
    luGc(:, 2) = min(luGc(:, 2), this%upperBounds)

    call this%gridcoordToCartesian(luGc, roundedLuCc)
    call subgrid%cartesianToRealGridcoord(roundedLuCc - spread(position, 2, 2), luSubGc)

  end subroutine TGrid_getSubgridRanges


  subroutine TGrid_hasSubgrid(this, subgrid, position, isInside)

    !> representation of parallelepiped grid
    class(TGrid), intent(in) :: this

    !> grid to calculate intersection with
    class(TGrid), intent(in) :: subgrid

    !> position of subgrid, expected shape: [3]
    real(dp), intent(in) :: position(:)

    !> true, if subgrid is a true subgrid of current grid, otherwise false
    logical, intent(out) :: isInside

    !> lower and upper cartesian bounds of grid and (partly) contained subgrid
    real(dp), allocatable :: lCcTotGrid(:,:), uCcTotGrid(:,:), lCcSubGrid(:,:), uCcSubGrid(:,:)

    isInside = .true.

    call this%gridcoordToCartesian(reshape(this%lowerBounds, [3, 1]), lCcTotGrid)
    call this%gridcoordToCartesian(reshape(this%upperBounds, [3, 1]), uCcTotGrid)
    call subgrid%gridcoordToCartesian(reshape(subgrid%lowerBounds, [3, 1]), lCcSubGrid)
    call subgrid%gridcoordToCartesian(reshape(subgrid%upperBounds, [3, 1]), uCcSubGrid)

    if (any(lCcSubGrid(:, 1) + position < lCcTotGrid(:, 1))) then
      isInside = .false.
    elseif (any(uCcSubGrid(:, 1) + position > uCcTotGrid(:, 1))) then
      isInside = .false.
    end if

  end subroutine TGrid_hasSubgrid


  subroutine TGrid_hasGridcoords(this, gridcoords, isInside)

    !> representation of parallelepiped grid
    class(TGrid), intent(in) :: this

    !> integer grid coordinates, shape: [3, nGridPoints]
    integer, intent(in) :: gridcoords(:,:)

    !> true, if the grid coordinate is within the grid region, shape: [nGridPoints]
    logical, intent(out), allocatable :: isInside(:)

    !> auxiliary variable
    integer :: ii

    allocate(isInside(size(gridcoords, dim=2)))

    do ii = 1, size(gridcoords, dim=2)
      if (any(gridcoords(:, ii) < this%lowerBounds) .or. &
          & any(gridcoords(:, ii) > this%upperBounds)) then
        isInside(ii) = .false.
      else
        isInside(ii) = .true.
      end if
    end do

  end subroutine TGrid_hasGridcoords


  subroutine getDistances(grid, position, distances)

    !> representation of a grid
    type(TGrid), intent(in) :: grid

    !> reference point for calculating the distances
    real(dp), intent(in) :: position(:)

    !> contains calculated distances on exit
    real(dp), intent(out), allocatable :: distances(:,:,:)

    !> cartesian coordinates of current grid point
    real(dp), allocatable :: cartcoords(:,:)

    !> data indices, corresponding to a grid point
    integer :: datainds(3)

    !> auxiliary variables
    integer :: ii, jj, kk, idx(3, 1)

    allocate(distances(grid%nPoints(1), grid%nPoints(2), grid%nPoints(3)))

    do kk = grid%lowerBounds(3), grid%upperBounds(3)
      idx(3, 1) = kk
      do jj = grid%lowerBounds(2), grid%upperBounds(2)
        idx(2, 1) = jj
        do ii = grid%lowerBounds(1), grid%upperBounds(1)
          idx(1, 1) = ii
          datainds(:) = idx(:, 1) - grid%lowerBounds + 1
          call grid%gridcoordToCartesian(idx, cartcoords)
          distances(datainds(1), datainds(2), datainds(3)) = norm2(cartcoords(:, 1) - position)
        end do
      end do
    end do

  end subroutine getDistances


  subroutine TGridData_init(this, grid, data, rwInterType, gridInterType)

    !> representation of volumetric grid and data
    type(TGridData), intent(out) :: this

    !> representation of parallelepiped grid
    type(TGrid), intent(in) :: grid

    !> volumetric data pointer
    real(dp), intent(in), pointer :: data(:,:,:)

    !> interpolation type to use for obtaining radial wavefunctions at arbitrary distances
    character(len=*), intent(in), optional :: rwInterType

    !> interpolation type to use for adding subgrid data
    character(len=*), intent(in), optional :: gridInterType

    if (present(rwInterType)) then
      call this%setRwInterpolationType(rwInterType)
    else
      call this%setRwInterpolationType('spline')
    end if

    if (present(gridInterType)) then
      call this%setGridInterpolationType(gridInterType)
    else
      call this%setGridInterpolationType('linear')
    end if

    this%grid = grid
    this%data => data

  end subroutine TGridData_init


  subroutine TGridData_setRwInterpolationType(this, interType)

    !> representation of volumetric grid and data
    class(TGridData), intent(inout) :: this

    !> interpolation type to use for obtaining radial wavefunctions at arbitrary distances
    character(len=*), intent(in) :: interType

    select case(trim(interType))
      case('trivial')
        this%tabulateBasis => TGridData_tabulateBasisTrivial
      case('linear')
        this%tabulateBasis => TGridData_tabulateBasisLinear
      case('spline')
        this%tabulateBasis => TGridData_tabulateBasisSpline
      case default
        this%tabulateBasis => TGridData_tabulateBasisSpline
    end select

  end subroutine TGridData_setRwInterpolationType


  subroutine TGridData_setGridInterpolationType(this, interType)

    !> representation of volumetric grid and data
    class(TGridData), intent(inout) :: this

    !> interpolation type to use for adding subgrid data
    character(len=*), intent(in) :: interType

    select case(trim(interType))

      case('trivial')
        this%add => TGridData_addTrivial

      case('linear')
        this%add => TGridData_addLinear

      case default
        this%add => TGridData_addLinear

    end select

  end subroutine TGridData_setGridInterpolationType


  subroutine TGridData_setGridOrigin(this, origin)

    !> representation of volumetric grid and data
    class(TGridData), intent(inout) :: this

    !> cartesian grid coordinates of grid origin
    real(dp), intent(in) :: origin(:)

    call this%grid%setOrigin(origin)

  end subroutine TGridData_setGridOrigin


  ! subroutine TGridData_getSubgridDataview(this, subgrid, position, dataview)

  !   !> representation of volumetric grid and data
  !   class(TGridData), intent(in), target :: this

  !   !> subgrid, which must be fully contained
  !   type(TGrid), intent(in) :: subgrid

  !   !> position where to place subgrid onto grid, shape: [3]
  !   real(dp), intent(in) :: position(:)

  !   !> part of the data array, corresponding to the subgrid
  !   real(dp), intent(out), pointer :: dataview(:,:,:)

  !   !> lower and upper bounds of contained subgrid in integer grid coordinates
  !   integer, allocatable :: luGc(:,:)

  !   call this%grid%getSubgridRanges(subgrid, position, luGc)

  !   dataview => this%data(&
  !       & luGc(1, 1) + 1:luGc(1, 2) + 1,&
  !       & luGc(2, 1) + 1:luGc(2, 2) + 1,&
  !       & luGc(3, 1) + 1:luGc(3, 2) + 1)

  ! end subroutine TGridData_getSubgridDataview


  subroutine TGridData_getSubcubeCornsAndVals(this, gridcoords, data)

    !> representation of volumetric grid and data
    class(TGridData), intent(in) :: this

    !> integer grid points to calculate subcubes for
    !> expected shape: [3, nxPoints, nyPoints, nzPoints]
    integer, intent(in) :: gridcoords(:,:,:,:)

    !> data corresponding to the eight corners of the subcube
    !> expected shape: [2, 2, 2, nxPoints, nyPoints, nzPoints]
    real(dp), intent(out), allocatable :: data(:,:,:,:,:,:)

    !> data indices, corresponding to an integer grid point
    integer :: inds(3)

    !> auxiliary variables
    integer :: ii, jj, kk

    allocate(data(2, 2, 2, size(gridcoords, dim=2), size(gridcoords, dim=3),&
        & size(gridcoords, dim=4)))

    do kk = 1, size(gridcoords, dim=4)
      do jj = 1, size(gridcoords, dim=3)
        do ii = 1, size(gridcoords, dim=2)
          inds(:) = gridcoords(:, ii, jj, kk) - this%grid%lowerBounds + 1
          inds(:) = max(inds, this%grid%lowerBounds)
          inds(:) = min(inds, this%grid%upperBounds - 1)
          data(1, 1, 1, ii, jj, kk) = this%data(inds(1) + 0, inds(2) + 0, inds(3) + 0)
          data(2, 1, 1, ii, jj, kk) = this%data(inds(1) + 1, inds(2) + 0, inds(3) + 0)
          data(1, 2, 1, ii, jj, kk) = this%data(inds(1) + 0, inds(2) + 1, inds(3) + 0)
          data(1, 1, 2, ii, jj, kk) = this%data(inds(1) + 0, inds(2) + 0, inds(3) + 1)
          data(2, 2, 1, ii, jj, kk) = this%data(inds(1) + 1, inds(2) + 1, inds(3) + 0)
          data(2, 1, 2, ii, jj, kk) = this%data(inds(1) + 1, inds(2) + 0, inds(3) + 1)
          data(1, 2, 2, ii, jj, kk) = this%data(inds(1) + 0, inds(2) + 1, inds(3) + 1)
          data(2, 2, 2, ii, jj, kk) = this%data(inds(1) + 1, inds(2) + 1, inds(3) + 1)
        end do
      end do
    end do

  end subroutine TGridData_getSubcubeCornsAndVals


  subroutine TGridData_tabulateBasisTrivial(this, rwf, rty, ll, mm)

    !> representation of volumetric grid and data
    class(TGridData), intent(inout) :: this

    !> tabulated radial wavefunction, shape: [nDistances, 3]
    real(dp), intent(in) :: rwf(:,:)

    !> tabulated real tesseral spherical harmonics for all l,m-combinations
    class(TRealTessY), intent(in) :: rty

    !> angular momentum of the species
    integer, intent(in) :: ll

    !> magnetic quantum number
    integer, intent(in) :: mm

    !> cartesian distances of grid points to origin
    real(dp), allocatable :: distances(:,:,:)

    !> data indices, corresponding to a grid point
    integer :: datainds(3)

    !> index of nearest distance
    integer :: loc

    !> auxiliary variables
    integer :: ii, jj, kk, rtyIdx, endIdx, idx(3)

    rtyIdx = ll * (ll + 1) + 1 + mm
    endIdx = this%grid%nPoints(1) * this%grid%nPoints(2) * this%grid%nPoints(3) - 1

    call getDistances(this%grid, [0.0_dp, 0.0_dp, 0.0_dp], distances)

    do kk = this%grid%lowerBounds(3), this%grid%upperBounds(3)
      idx(3) = kk
      do jj = this%grid%lowerBounds(2), this%grid%upperBounds(2)
        idx(2) = jj
        do ii = this%grid%lowerBounds(1), this%grid%upperBounds(1)
          idx(1) = ii
          datainds(:) = idx - this%grid%lowerBounds + 1
          loc = max(min(locateByBisection(rwf(:, 1),&
              & distances(datainds(1), datainds(2), datainds(3))), endIdx), 1)
          this%data(datainds(1), datainds(2), datainds(3)) =&
              & this%data(datainds(1), datainds(2), datainds(3))&
              & + rwf(loc, 2) * rty%data(datainds(1), datainds(2), datainds(3), rtyIdx)
        end do
      end do
    end do

  end subroutine TGridData_tabulateBasisTrivial


  subroutine TGridData_tabulateBasisLinear(this, rwf, rty, ll, mm)

    !> representation of volumetric grid and data
    class(TGridData), intent(inout) :: this

    !> tabulated radial wavefunction, shape: [nDistances, 3]
    real(dp), intent(in) :: rwf(:,:)

    !> tabulated real tesseral spherical harmonics for all l,m-combinations
    class(TRealTessY), intent(in) :: rty

    !> angular momentum of the species
    integer, intent(in) :: ll

    !> magnetic quantum number
    integer, intent(in) :: mm

    !> cartesian distances of grid points to origin
    real(dp), allocatable :: distances(:,:,:)

    !> data indices, corresponding to a grid point
    integer :: datainds(3)

    !> index of nearest distance
    integer :: loc

    !> auxiliary variables
    integer :: ii, jj, kk, rtyIdx, endIdx, idx(3)

    rtyIdx = ll * (ll + 1) + 1 + mm
    endIdx = this%grid%nPoints(1) * this%grid%nPoints(2) * this%grid%nPoints(3) - 1

    call getDistances(this%grid, [0.0_dp, 0.0_dp, 0.0_dp], distances)

    do kk = this%grid%lowerBounds(3), this%grid%upperBounds(3)
      idx(3) = kk
      do jj = this%grid%lowerBounds(2), this%grid%upperBounds(2)
        idx(2) = jj
        do ii = this%grid%lowerBounds(1), this%grid%upperBounds(1)
          idx(1) = ii
          datainds(:) = idx - this%grid%lowerBounds + 1
          loc = max(min(locateByBisection(rwf(:, 1),&
              & distances(datainds(1), datainds(2), datainds(3))), endIdx), 1)
          this%data(datainds(1), datainds(2), datainds(3)) =&
              & this%data(datainds(1), datainds(2), datainds(3))&
              & + linearInterpolation(rwf(loc, 1), rwf(loc + 1, 1), rwf(loc, 2), rwf(loc + 1, 2),&
              & distances(datainds(1), datainds(2), datainds(3)))&
              & * rty%data(datainds(1), datainds(2), datainds(3), rtyIdx)
        end do
      end do
    end do

  end subroutine TGridData_tabulateBasisLinear


  subroutine TGridData_tabulateBasisSpline(this, rwf, rty, ll, mm)

    !> representation of volumetric grid and data
    class(TGridData), intent(inout) :: this

    !> tabulated radial wavefunction, shape: [nDistances, 3]
    real(dp), intent(in) :: rwf(:,:)

    !> tabulated real tesseral spherical harmonics for all l,m-combinations
    class(TRealTessY), intent(in) :: rty

    !> angular momentum of the species
    integer, intent(in) :: ll

    !> magnetic quantum number
    integer, intent(in) :: mm

    !> cartesian distances of grid points to origin
    real(dp), allocatable :: distances(:,:,:)

    !> data indices, corresponding to a grid point
    integer :: datainds(3)

    !> auxiliary variables
    integer :: ii, jj, kk, rtyIdx, endIdx, idx(3)

    rtyIdx = ll * (ll + 1) + 1 + mm
    endIdx = this%grid%nPoints(1) * this%grid%nPoints(2) * this%grid%nPoints(3) - 1

    call getDistances(this%grid, [0.0_dp, 0.0_dp, 0.0_dp], distances)

    do kk = this%grid%lowerBounds(3), this%grid%upperBounds(3)
      idx(3) = kk
      do jj = this%grid%lowerBounds(2), this%grid%upperBounds(2)
        idx(2) = jj
        do ii = this%grid%lowerBounds(1), this%grid%upperBounds(1)
          idx(1) = ii
          datainds(:) = idx - this%grid%lowerBounds + 1
          this%data(datainds(1), datainds(2), datainds(3)) =&
              & this%data(datainds(1), datainds(2), datainds(3)) +&
              & splineInterpolation(rwf(:, 1), rwf(:, 2), rwf(:, 3), distances(datainds(1),&
              & datainds(2), datainds(3))) * rty%data(datainds(1), datainds(2), datainds(3), rtyIdx)
        end do
      end do
    end do

  end subroutine TGridData_tabulateBasisSpline


  subroutine TGrid_getIntersecGridPoints(this, subgrid, luGc, luSubGc, subgridGc)

    !> representation of parallelepiped grid
    class(TGrid), intent(in) :: this

    !> subgrid to calculate intersection with
    class(TGrid), intent(in) :: subgrid

    !> lower and upper bounds of intersection in integer grid coordinates
    integer, intent(in) :: luGc(:,:)

    !> lower and upper bounds of intersection in real subgrid coordinates
    real(dp), intent(in) :: luSubGc(:,:)

    !> real grid coordinates of subgrid that intersect with grid
    real(dp), intent(out), allocatable :: subgridGc(:,:,:,:)

    !> auxiliary variable
    real(dp) :: step(3), gc(3)

    !> auxiliary variables
    integer :: ii, jj, kk, idx(3)

    allocate(subgridGc(3,&
        & luGc(1, 2) - luGc(1, 1) + 1,&
        & luGc(2, 2) - luGc(2, 1) + 1,&
        & luGc(3, 2) - luGc(3, 1) + 1))

    step = norm2(this%basis, dim=1) / norm2(subgrid%basis, dim=1)

    do kk = 1, luGc(3, 2) - luGc(3, 1) + 1
      idx(3) = kk
      do jj = 1, luGc(2, 2) - luGc(2, 1) + 1
        idx(2) = jj
        do ii = 1, luGc(1, 2) - luGc(1, 1) + 1
          idx(1) = ii
          gc(:) = luSubGc(:, 1) + (idx - 1.0_dp) * step
          gc(:) = max(gc, real(subgrid%lowerBounds, dp))
          gc(:) = min(gc, real(subgrid%upperBounds, dp))
          subgridGc(:, ii, jj, kk) = gc
        end do
      end do
    end do

  end subroutine TGrid_getIntersecGridPoints


  subroutine TGridData_addTrivial(this, griddata, eigvec, position, square)

    !> representation of volumetric grid and data
    class(TGridData), intent(inout) :: this

    !> volumetric data to add to current instance
    class(TGridData), intent(in) :: griddata

    !> eigenvector of current basis
    real(dp), intent(in) :: eigvec

    !> position where to place subgrid onto grid, expected shape: [3]
    real(dp), intent(in) :: position(:)

    !> True, if data should be squared, i.e. for densities
    logical, intent(in), optional :: square

    !> subgrid points that intersect with surrounding grid
    real(dp), allocatable :: intersecSubGc(:,:,:,:)

    !> lower and upper bounds of intersection in integer grid coordinates
    integer, allocatable :: luGc(:,:)

    !> lower and upper bounds of intersection in real subgrid coordinates
    real(dp), allocatable :: luSubGc(:,:)

    !> interpolated data in intersection area
    real(dp), allocatable :: interp(:,:,:)

    !> equals square if present, otherwise false
    logical :: tSquare

    !> data indices, corresponding to a grid point
    integer :: inds(3)

    !> auxiliary variables
    integer :: ii, jj, kk, idx(3)

    if (present(square)) then
      tSquare = square
    else
      tSquare = .false.
    end if

    call checkParallelBasis(this%grid%basis, griddata%grid%basis)

    call this%grid%getSubgridRanges(griddata%grid, position, luGc, luSubGc)
    call this%grid%getIntersecGridPoints(griddata%grid, luGc, luSubGc, intersecSubGc)
    call griddata%getInterpolatedValueGc(intersecSubGc, interp)

    if (tSquare) then
      interp = interp**2
    end if

    do kk = 1, luGc(3, 2) - luGc(3, 1) + 1
      idx(3) = kk
      do jj = 1, luGc(2, 2) - luGc(2, 1) + 1
        idx(2) = jj
        do ii = 1, luGc(1, 2) - luGc(1, 1) + 1
          idx(1) = ii
          inds(:) = idx + luGc(:, 1) - this%grid%lowerBounds
          this%data(inds(1), inds(2), inds(3)) = this%data(inds(1), inds(2), inds(3))&
              & + eigvec * interp(ii, jj, kk)
        end do
      end do
    end do

  end subroutine TGridData_addTrivial


  subroutine TGridData_getInterpolatedValueGc(this, gridcoords, data)

    !> representation of volumetric grid and data
    class(TGridData), intent(in) :: this

    !> real subgrid coordinates to calculate interpolation for
    real(dp), intent(in) :: gridcoords(:,:,:,:)

    !> corresponding (interpolated) data
    real(dp), intent(out), allocatable :: data(:,:,:)

    call this%getValue(floor(gridcoords), data)

  end subroutine TGridData_getInterpolatedValueGc


  subroutine TGridData_getValue(this, gridcoords, vals)

    !> representation of volumetric grid and data
    class(TGridData), intent(in) :: this

    !> integer grid coordinates to return data for
    integer, intent(in) :: gridcoords(:,:,:,:)

    !> data points corresponding to integer grid coordinates
    real(dp), intent(out), allocatable :: vals(:,:,:)

    !> data indices, specifying a data access point
    integer :: inds(3)

    !> auxiliary variables
    integer :: ii, jj, kk

    allocate(vals(size(gridcoords, dim=2), size(gridcoords, dim=3), size(gridcoords, dim=4)))

    do kk = 1, size(gridcoords, dim=4)
      do jj = 1, size(gridcoords, dim=3)
        do ii = 1, size(gridcoords, dim=2)
          inds(:) = gridcoords(:, ii, jj, kk) - this%grid%lowerBounds + 1
          vals(ii, jj, kk) = this%data(inds(1), inds(2), inds(3))
        end do
      end do
    end do

  end subroutine TGridData_getValue


  subroutine TGridData_addLinear(this, griddata, eigvec, position, square)

    !> representation of volumetric grid and data
    class(TGridData), intent(inout) :: this

    !> volumetric data to add to current instance
    class(TGridData), intent(in) :: griddata

    !> eigenvector of current basis
    real(dp), intent(in) :: eigvec

    !> position where to place subgrid onto grid, expected shape: [3]
    real(dp), intent(in) :: position(:)

    !> True, if data should be squared, i.e. for densities
    logical, intent(in), optional :: square

    !> lower and upper bounds of intersection in integer grid coordinates
    integer, allocatable :: luGc(:,:)

    !> lower and upper bounds of intersection in real subgrid coordinates
    real(dp), allocatable :: luSubGc(:,:)

    !> lower subcube corners of grid points
    real(dp), allocatable :: intersecSubGc(:,:,:,:)

    !> floor rounded lower subcube corners of grid points
    integer, allocatable :: roundedIntersecSubGc(:,:,:,:)

    !> data corresponding to the eight corners of the subcubes
    real(dp), allocatable :: scData(:,:,:,:,:,:)

    !> equals square if present, otherwise false
    logical :: tSquare

    !> auxiliary variables
    integer :: ii, jj, kk, inds(3), idxval(3)

    if (present(square)) then
      tSquare = square
    else
      tSquare = .false.
    end if

    call checkParallelBasis(this%grid%basis, griddata%grid%basis)

    call this%grid%getSubgridRanges(griddata%grid, position, luGc, luSubGc)

    call this%grid%getIntersecGridPoints(griddata%grid, luGc, luSubGc, intersecSubGc)
    roundedIntersecSubGc = floor(intersecSubGc)
    call griddata%getSubcubeCornsAndVals(roundedIntersecSubGc, scData)

    if (tSquare) then
      scData(:,:,:,:,:,:) = scData * scData
    end if

    do kk = 1, luGc(3, 2) - luGc(3, 1) + 1
      idxval(3) = kk
      do jj = 1, luGc(2, 2) - luGc(2, 1) + 1
        idxval(2) = jj
        do ii = 1, luGc(1, 2) - luGc(1, 1) + 1
          idxval(1) = ii
          inds(:) = idxval + luGc(:, 1) - this%grid%lowerBounds
          this%data(inds(1), inds(2), inds(3)) = this%data(inds(1), inds(2), inds(3))&
              & + eigvec * trilinearInterpolation(real(roundedIntersecSubGc(:, ii, jj, kk), dp),&
              & real(roundedIntersecSubGc(:, ii, jj, kk) + 1, dp), scData(:,:,:, ii, jj, kk),&
              & intersecSubGc(:, ii, jj, kk))
        end do
      end do
    end do

  end subroutine TGridData_addLinear


  subroutine subgridsToUncachedGlobalGrid(globalGridDat, subgridsDat, coords, coeffs,&
      & orbitalToAtom, orbitalToSubgrid, nRegions, addDensities)

    !> representation of volumetric grid and data
    type(TGridData), intent(inout) :: globalGridDat

    !> subgrids to add to global grid instance, exp. shape: [nOrbitals]
    type(TGridData), intent(in) :: subgridsDat(:)

    !> coordinates, defining the cartesian position of the subgrids, exp. shape: [3, nAtom, nCell]
    real(dp), intent(in) :: coords(:,:,:)

    !> orbital-resolved coefficients for volumetric data, exp. shape: [nOrbitals]
    real(dp), intent(in) :: coeffs(:)

    !> index mapping from orbital to atom, exp. shape: [nOrbitals]
    integer, intent(in) :: orbitalToAtom(:)

    !> index mapping from orbital to subgrid, exp. shape: [nOrbitals]
    integer, intent(in) :: orbitalToSubgrid(:)

    !> number of parallel global regions
    integer, intent(in) :: nRegions

    !> add densities instead of wave functions
    logical, intent(in), optional :: addDensities

    !> grid parallel regions, holding global grid sections
    type(TGridData), allocatable :: regionGridDat(:)

    !> start and end indices of all parallel tiles, exp. shape: [2, nRegions]
    integer, allocatable :: tiling(:,:)

    !> true, if optional argument provided and .true., otherwise False
    logical :: tAddDensities

    !> temporary grid instance
    type(TGrid) :: tmpGrid

    !> temporary data containers
    real(dp), allocatable, target :: regionData(:,:,:)
    real(dp), pointer :: pRegionData(:,:,:), pTmpData(:,:,:)

    !> auxiliary variables
    integer :: iGrid, iOrb, iAtom, iCell, iSubgrid

    ! try to ensure a smooth runtime
    @:ASSERT(size(coeffs) == size(orbitalToAtom))
    @:ASSERT(size(coeffs) == size(orbitalToSubgrid))
    @:ASSERT(size(subgridsDat) == maxval(orbitalToSubgrid))
    @:ASSERT(size(coords, dim=1) == 3)
    @:ASSERT(size(coords, dim=2) == maxval(orbitalToAtom))
    @:ASSERT(nRegions >= 1)

    if (present(addDensities)) then
      tAddDensities = addDensities
    else
      tAddDensities = .false.
    end if

    ! get start and end index of every parallel region of the global grid instance
    call getStartAndEndIndices(size(globalGridDat%data, dim=3), nRegions, tiling)

    ! allocate regional global grid tiles
    allocate(regionGridDat(nRegions))

    ! assign subgrids to parallel regions of the total grid
    do iGrid = 1, nRegions
      tmpGrid = globalGridDat%grid
      tmpGrid%lowerBounds(3) = tiling(1, iGrid)
      tmpGrid%upperBounds(3) = tiling(2, iGrid)
      tmpGrid%nPoints(3) = tiling(2, iGrid) - tiling(1, iGrid) + 1
      pTmpData => globalGridDat%data(:,:, tiling(1, iGrid):tiling(2, iGrid))
      call TGridData_init(regionGridDat(iGrid), tmpGrid, pTmpData, rwInterType='spline',&
          & gridInterType='linear')
    end do

    !$omp parallel do default(none)&
    !$omp& shared(nRegions, regionGridDat, subgridsDat, orbitalToAtom, orbitalToSubgrid,&
    !$omp& coords, coeffs, tAddDensities)&
    !$omp& private(iGrid, iOrb, iAtom, iSubgrid, iCell)
    lpGrid: do iGrid = 1, nRegions
      lpOrb: do iOrb = 1, size(coeffs)

        iAtom = orbitalToAtom(iOrb)
        iSubgrid = orbitalToSubgrid(iOrb)

        lpCell: do iCell = 1, size(coords, dim=3)

          call regionGridDat(iGrid)%add(subgridsDat(iSubgrid), coeffs(iOrb),&
              & coords(:, iAtom, iCell), square=tAddDensities)

        end do lpCell

      end do lpOrb
    end do lpGrid
    !$omp end parallel do

  end subroutine subgridsToUncachedGlobalGrid


  subroutine subgridsToCachedGlobalGrids(globalGridDat, subgridsDat, coords, coeffs, levelIndex,&
      & orbitalToAtom, orbitalToSubgrid, nRegions, pCopyBuffers, globalGridsDat, addDensities)

    !> prototype of global grid instance
    type(TGridData), intent(in) :: globalGridDat

    !> subgrids to add to global grid instances, exp. shape: [nOrbitals]
    type(TGridData), intent(in) :: subgridsDat(:)

    !> coordinates, defining the cartesian positions of the subgrids, exp. shape: [3, nAtom, nCell]
    real(dp), intent(in) :: coords(:,:,:)

    !> orbital-resolved coefficients for volumetric data, exp. shape: [nOrbitals, nCache]
    real(dp), intent(in) :: coeffs(:,:)

    !> Index Mapping from state index to [Level, KPoint, Spin]
    integer, intent(in) :: levelIndex(:,:)

    !> index mapping from orbital to atom, exp. shape: [nOrbitals]
    integer, intent(in) :: orbitalToAtom(:)

    !> index mapping from orbital to subgrid, exp. shape: [nOrbitals]
    integer, intent(in) :: orbitalToSubgrid(:)

    !> number of parallel global regions
    integer, intent(in) :: nRegions

    !> pointer to target array that holds the volumetric data of grid caches
    real(dp), intent(in), pointer :: pCopyBuffers(:,:,:,:)

    !> representation of volumetric grids and data, exp. shape: [nCache]
    type(TGridData), intent(out), allocatable :: globalGridsDat(:)

    !> add densities instead of wave functions
    logical, intent(in), optional :: addDensities

    !> grid parallel regions, holding global grid sections
    type(TGridData), allocatable :: regionGridDat(:)

    !> start and end indices of all parallel tiles, exp. shape: [2, nRegions]
    integer, allocatable :: tiling(:,:)

    !> true, if optional argument provided and .true., otherwise False
    logical :: tAddDensities

    !> temporary grid instance
    type(TGrid) :: tmpGrid

    !> temporary data containers
    real(dp), allocatable, target :: regionData(:,:,:)
    real(dp), pointer :: pRegionData(:,:,:), pTmpData(:,:,:)

    !> auxiliary variables
    integer :: iGrid, iOrb, iCell, iAtom, iSubgrid, iState

    ! try to ensure a smooth runtime
    @:ASSERT(size(coeffs, dim=1) == size(orbitalToAtom))
    @:ASSERT(size(coeffs, dim=1) == size(orbitalToSubgrid))
    @:ASSERT(size(subgridsDat) == maxval(orbitalToSubgrid))
    @:ASSERT(size(coords, dim=1) == 3)
    @:ASSERT(size(coords, dim=2) == maxval(orbitalToAtom))
    @:ASSERT(size(coords, dim=3) >= 1)
    @:ASSERT(nRegions >= 1)

    if (present(addDensities)) then
      tAddDensities = addDensities
    else
      tAddDensities = .false.
    end if

    ! generate global grid caches from prototype
    allocate(globalGridsDat(size(coeffs, dim=2)))
    globalGridsDat(:) = globalGridDat

    ! separate data pointers for grid caches
    do iState = 1, size(levelIndex, dim=2)
      globalGridsDat(iState)%data => pCopyBuffers(:,:,:, iState)
    end do

    ! get start and end index of every parallel region of the global grid instance
    call getStartAndEndIndices(size(globalGridDat%data, dim=3), nRegions, tiling)

    ! allocate regional global grid tiles
    allocate(regionGridDat(nRegions))

    ! assign subgrids to parallel regions of the total grid
    do iGrid = 1, nRegions
      tmpGrid = globalGridDat%grid
      tmpGrid%lowerBounds(3) = tiling(1, iGrid)
      tmpGrid%upperBounds(3) = tiling(2, iGrid)
      tmpGrid%nPoints(3) = tiling(2, iGrid) - tiling(1, iGrid) + 1
      pTmpData => globalGridDat%data(:,:, tiling(1, iGrid):tiling(2, iGrid))
      call TGridData_init(regionGridDat(iGrid), tmpGrid, pTmpData, rwInterType='spline',&
          & gridInterType='linear')
    end do

    !$omp parallel do default(none)&
    !$omp& shared(nRegions, regionGridDat, subgridsDat, orbitalToAtom, orbitalToSubgrid,&
    !$omp& coords, coeffs, levelIndex, globalGridsDat, tiling, tAddDensities)&
    !$omp& private(iGrid, iOrb, iAtom, iSubgrid, iCell, iState)
    lpGrid: do iGrid = 1, nRegions
      lpOrb: do iOrb = 1, size(coeffs, dim=1)

        ! reset volumetric grid values
        regionGridDat(iGrid)%data(:,:,:) = 0.0_dp

        iAtom = orbitalToAtom(iOrb)
        iSubgrid = orbitalToSubgrid(iOrb)

        lpCell: do iCell = 1, size(coords, dim=3)

          call regionGridDat(iGrid)%add(subgridsDat(iSubgrid), 1.0_dp, coords(:, iAtom, iCell),&
              & square=tAddDensities)

        end do lpCell

        ! distribute molecular orbitals to state-resolved global grid caches
        lpStates: do iState = 1, size(levelIndex, dim=2)

          globalGridsDat(iState)%data(:,:, tiling(1, iGrid):tiling(2, iGrid)) =&
              & globalGridsDat(iState)%data(:,:, tiling(1, iGrid):tiling(2, iGrid))&
              & + coeffs(iOrb, iState) * regionGridDat(iGrid)%data

        end do lpStates

      end do lpOrb
    end do lpGrid
    !$omp end parallel do

  end subroutine subgridsToCachedGlobalGrids


  subroutine TRealTessY_init(this, grid, maxAng)

    !> representation of real tesseral spherical harmonics
    class(TRealTessY), intent(inout) :: this

    !> representation of parallelepiped grid
    class(TGrid), intent(in) :: grid

    !> maximum angular momentum of the system
    integer, intent(in) :: maxAng

    !> angular momentum and magnetic quantum number
    integer :: ll, mm

    this%grid = grid
    this%maxAng = maxAng

    allocate(this%lmCombs((this%maxAng + 1)**2, 2))

    do ll = 0, this%maxAng
      do mm = - ll, ll
        this%lmCombs(ll * (ll + 1) + 1 + mm, 1) = ll
        this%lmCombs(ll * (ll + 1) + 1 + mm, 2) = mm
      end do
    end do

    allocate(this%data(this%grid%nPoints(1), this%grid%nPoints(2), this%grid%nPoints(3),&
        & size(this%lmCombs, dim=1)))

    call this%tabulateRealTessY()

  end subroutine TRealTessY_init


  subroutine TRealTessY_tabulateRealTessY(this)

    !> representation of real tesseral spherical harmonics
    class(TRealTessY), intent(inout) :: this

    !> cartesian coordinates of current grid point
    real(dp), allocatable :: cartcoords(:,:)

    !> auxiliary variables
    integer :: ii, jj, kk, lmComb, idx(3, 1), datainds(3)

    do kk = this%grid%lowerBounds(3), this%grid%upperBounds(3)
      idx(3, 1) = kk
      do jj = this%grid%lowerBounds(2), this%grid%upperBounds(2)
        idx(2, 1) = jj
        do ii = this%grid%lowerBounds(1), this%grid%upperBounds(1)
          idx(1, 1) = ii
          call this%grid%gridcoordToCartesian(idx, cartcoords)
          do lmComb = 1, size(this%lmCombs, dim=1)
            datainds(:) = idx(:, 1) - this%grid%lowerBounds + 1
            this%data(datainds(1), datainds(2), datainds(3), lmComb) =&
                & RealTessY(this%lmCombs(lmComb, 1), this%lmCombs(lmComb, 2), cartcoords(:, 1))
          end do
        end do
      end do
    end do

  end subroutine TRealTessY_tabulateRealTessY


  subroutine checkParallelBasis(basis1, basis2)

    !> basis vectors of grids to compare
    real(dp), intent(in) :: basis1(:,:), basis2(:,:)

    !> cross product of two vectors
    real(dp) :: cross(3)

    !> auxiliary variable
    integer :: ii

    do ii = 1, 3
      cross = crossProduct(basis1(:, ii), basis2(:, ii))
      if (norm2(cross) > 1e-10_dp) then
        call error('Non-parallel basis vectors not supported.')
      end if
    end do

  end subroutine checkParallelBasis


  pure function crossProduct(vec1, vec2) result(cross)

    !> vectors to calculate the cross product for
    real(dp), intent(in) :: vec1(:), vec2(:)

    !> resulting cross product
    real(dp) :: cross(3)

    cross(1) = vec1(2) * vec2(3) - vec1(3) * vec2(2)
    cross(2) = vec1(3) * vec2(1) - vec1(1) * vec2(3)
    cross(3) = vec1(1) * vec2(2) - vec1(2) * vec2(1)

  end function crossProduct

end module waveplot_grids
