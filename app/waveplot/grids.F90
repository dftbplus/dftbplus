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
  use dftbp_common_constants, only : pi
  use waveplot_interp, only : linearInterpolation, trilinearInterpolation, polynomialInterpolation,&
      & locateByBisection
  use waveplot_parallel, only : getStartAndEndIndices
  use waveplot_slater, only: RealTessY, SlaterOrbital_getValue_explicit, TSlaterOrbital, &
      & SlaterOrbital_getValue

#:if WITH_OMP
  use omp_lib
#:endif

  implicit none
  private

  real(dp), parameter :: roundingEpsilon = 1.0e-12_dp

  public :: TGrid, TGrid_init
  public :: TGridData, TGridData_init, TTabulationTypesEnum, TGridInterpolationTypesEnum
  public :: TRealTessY, TRealTessY_init
  public :: subgridsToGlobalGrid
  public :: modifyEigenvecs
  public :: parallelGridSlices
  public :: calcCartCoords


  ! Enumerates the different types which are used for the tabulation of the radial wf.
  type :: TTabulationTypesEnum

  ! Trivial interpolation
  integer :: trivial = 0

  ! Linear interpolation
  integer :: linear = 1

  ! Polynomial interpolation
  integer :: polynomial = 2

  ! Explicit calculation
  integer :: explicit = 3

  end type TTabulationTypesEnum

  !> Container for enumerated radial wf tabulation types
  type(TTabulationTypesEnum), parameter :: rwTabulationTypes = TTabulationTypesEnum()


  ! Enumerates the two different types which are used for the grid interpolation.
  type :: TGridInterpolationTypesEnum

  ! Trivial interpolation
  integer :: trivial = 0

  ! Linear interpolation
  integer :: linear = 1

  ! Explicit calculation, no interpolation
  integer :: explicit = 2

  end type TGridInterpolationTypesEnum

  !> Container for enumerated grid interpolation types
  type(TGridInterpolationTypesEnum), parameter :: gridInterpolTypes = TGridInterpolationTypesEnum()


  !> Contains properties of a three-dimensional, equidistant (not necessarily cubic) grid
  type :: TGrid

    !> cartesian coordinates of grid origin
    real(dp) :: origin(3)

    !> array of column vectors spanning dense grid basis
    real(dp) :: basis(3,3)

    !> inverted basis
    real(dp) :: invBasis(3,3)

    !> lower range bound for every axis, starting with 0 and ends at n-1
    integer :: lowerBounds(3)

    !> upper range bound for every axis, starting with 0 and ends at n-1
    integer :: upperBounds(3)

    !> number of grid points along each dimension
    integer :: nPoints(3)

  contains

    procedure :: setOrigin => TGrid_setOrigin
    procedure :: getCorners => TGrid_getCorners
    procedure :: cartesianToGridcoord => TGrid_cartesianToGridCoord
    procedure :: cartesianToRealGridcoord => TGrid_cartesianToRealGridCoord
    procedure :: gridcoordToCartesian => TGrid_gridCoordToCartesian
    procedure :: realGridcoordToCartesian => TGrid_realGridCoordToCartesian
    procedure :: getSubgridRanges => TGrid_getSubgridRanges
    procedure :: hasGridcoords => TGrid_hasGridCoords
    procedure :: hasSubgrid => TGrid_hasSubgrid
    procedure :: getIntersecGridPoints => TGrid_getIntersecGridPoints

  end type TGrid


  !> Data of the real tesseral spherical harmonics on a given grid
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


  !> Data on a given grid
  type :: TGridData

    !> representation of parallelepiped grid
    class(TGrid), allocatable :: grid

    !> tabulated real tesseral spherical harmonics for all l,m-combinations
    class(TRealTessY), allocatable :: rty

    !> volumetric data pointer, shape: [nxPoints, nyPoints, nzPoints]
    real(dp), pointer :: data(:,:,:)

    procedure(tabulateBasis), pointer, pass(this) :: tabulateBasis => null()

  contains

    procedure :: setGridOrigin => TGridData_setGridOrigin
    procedure :: setRwTabulationType => TGridData_setRwTabulationType
    procedure :: getValue => TGridData_getValue
    procedure :: getInterpolatedValueGc => TGridData_getInterpolatedValueGc
    procedure :: getSubcubeCornsAndVals => TGridData_getSubcubeCornsAndVals

  end type TGridData


  interface subgridsToGlobalGrid
    module procedure :: subgridsToUncachedGlobalGrid
    module procedure :: subgridsToCachedGlobalGrids
    module procedure :: subgridsToCachedGlobalGridsCplx
  end interface subgridsToGlobalGrid


  interface gridInterpolationLinear
    module procedure gridInterpolationRealLinear
    module procedure gridInterpolationCplxLinear
  end interface gridInterpolationLinear


  interface gridInterpolationTrivial
    module procedure gridInterpolationRealTrivial
    module procedure gridInterpolationCplxTrivial
  end interface gridInterpolationTrivial


  interface explicitSTOcalculation
    module procedure explicitSTOcalculationReal
    module procedure explicitSTOcalculationCplx
  end interface explicitSTOcalculation


  interface modifyEigenvecs
    module procedure modifyRealEigenvecs
    module procedure modifyCplxEigenvecs
  end interface modifyEigenvecs


contains

  !> Pretabulates the slater type orbitals
  subroutine tabulateBasis(this, rty, ll, mm, rwf, alpha, aa)

    implicit none

    !> representation of volumetric grid and data
    class(TGridData), intent(inout) :: this

    !> tabulated real tesseral spherical harmonics for all l,m-combinations
    class(TRealTessY), intent(in) :: rty

    !> angular momentum of the species
    integer, intent(in) :: ll

    !> magnetic quantum number
    integer, intent(in) :: mm

    !> tabulated radial wavefunction, shape: [nDistances, 3], it contains the radius, the value of
    !> the wf at this radius and the second derivative of the radial wf
    real(dp), intent(in) :: rwf(:,:)

    !> exponential coefficients for explicit calculation of the slater type orbitals,
    !> shape: [nAlpha]
    real(dp), intent(in) :: alpha(:)

    !> summation coefficients, shape: [nCoeffPerAlpha, nAlpha]
    real(dp), intent(in) :: aa(:,:)

  end subroutine tabulateBasis


  !> Initializes grid instance
  subroutine TGrid_init(this, origin, basis, nPoints)

    !> representation of parallelepiped grid
    class(TGrid), intent(inout) :: this

    !> cartesian coordinates of grid origin, shape [3]
    real(dp), intent(in) :: origin(:)

    !> array of column vectors spanning dense grid basis, shape [3,3]
    real(dp), intent(in) :: basis(:,:)

    !> number of points along each axis, shape [3]
    integer, intent(in) :: nPoints(:)

    this%origin = origin
    this%nPoints = nPoints
    this%basis = basis

    this%lowerBounds = [0, 0, 0]
    this%upperBounds = nPoints - 1

    call invert33(this%invBasis, this%basis)

  end subroutine TGrid_init


  !> Initializes the grid data instance
  subroutine TGridData_init(this, grid, data, rwTabulationType)

    !> representation of volumetric grid and data
    type(TGridData), intent(out) :: this

    !> representation of parallelepiped grid
    type(TGrid), intent(in) :: grid

    !> volumetric data pointer, shape: [nxPoints, nyPoints, nzPoints]
    real(dp), intent(in), pointer :: data(:,:,:)

    !> one dimensional tabulation type to use for obtaining radial wavefunctions at arbitrary
    !> distances
    integer, intent(in), optional :: rwTabulationType

    if (present(rwTabulationType)) then
      call this%setRwTabulationType(rwTabulationType)
    else
      call this%setRwTabulationType(rwTabulationTypes%explicit)
    end if

    this%grid = grid
    this%data => data

  end subroutine TGridData_init


  !> Intializes the grid instance of the real tesseral spherical harmonics
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


  !> Tabulates the values of the real tesseral spherical harmonics
  subroutine TRealTessY_tabulateRealTessY(this)

    !> representation of real tesseral spherical harmonics
    class(TRealTessY), intent(inout) :: this

    !> cartesian coordinates of current grid point, shape: [3, 1]
    real(dp), allocatable :: cartCoords(:,:)

    !> auxiliary variables
    integer :: ii, jj, kk, lmComb, idx(3, 1), datainds(3)

    do kk = this%grid%lowerBounds(3), this%grid%upperBounds(3)
      idx(3, 1) = kk
      do jj = this%grid%lowerBounds(2), this%grid%upperBounds(2)
        idx(2, 1) = jj
        do ii = this%grid%lowerBounds(1), this%grid%upperBounds(1)
          idx(1, 1) = ii
          call this%grid%gridcoordToCartesian(idx, cartCoords)
          do lmComb = 1, size(this%lmCombs, dim=1)
            datainds(:) = idx(:, 1) - this%grid%lowerBounds + 1
            this%data(datainds(1), datainds(2), datainds(3), lmComb) =&
                & RealTessY(this%lmCombs(lmComb, 1), this%lmCombs(lmComb, 2), cartCoords(:, 1))
          end do
        end do
      end do
    end do

  end subroutine TRealTessY_tabulateRealTessY


  !> Sets origin af a grid instance
  subroutine TGrid_setOrigin(this, origin)

    !> representation of parallelepiped grid
    class(TGrid), intent(inout) :: this

    !> cartesian coordinates of grid origin, shape [3]
    real(dp), intent(in) :: origin(:)

    this%origin = origin

  end subroutine TGrid_setOrigin


  !> Transformation from integer grid coordinates to cartesian coordinates
  subroutine TGrid_gridCoordToCartesian(this, gridcoords, cartcoords)

    !> representation of parallelepiped grid
    class(TGrid), intent(in) :: this

    !> integer grid coordinates, shape: [3, nCoords]
    integer, intent(in) :: gridcoords(:,:)

    !> corresponding cartesian grid coordinates, shape: [3, nCoords]
    real(dp), intent(out), allocatable :: cartcoords(:,:)

    !> auxiliary variable
    integer :: ii

    allocate(cartcoords(3, size(gridcoords, dim=2)))

    do ii = 1, size(gridcoords, dim=2)
      cartcoords(:, ii) = matmul(this%basis, gridcoords(:, ii)) + this%origin
    end do

  end subroutine TGrid_gridCoordToCartesian


  !> Transformation from real grid coordinates to cartesian coordinates
  subroutine TGrid_realGridCoordToCartesian(this, gridcoords, cartcoords)

    !> representation of parallelepiped grid
    class(TGrid), intent(in) :: this

    !> real grid coordinates, shape: [3, nCoords]
    real(dp), intent(in) :: gridcoords(:,:)

    !> corresponding cartesian grid coordinates, shape: [3, nCoords]
    real(dp), intent(out), allocatable :: cartcoords(:,:)

    !> auxiliary variable
    integer :: ii

    allocate(cartcoords(3, size(gridcoords, dim=2)))

    do ii = 1, size(gridcoords, dim=2)
      cartcoords(:, ii) = matmul(this%basis, gridcoords(:, ii)) + this%origin
    end do

  end subroutine TGrid_realGridCoordToCartesian


  !> Transformation from cartesian coordinates to integer grid coordinates
  subroutine TGrid_cartesianToGridCoord(this, cartcoords, gridcoords)

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
      gridcoords(:, ii) = floor(matmul(cartcoords(:, ii) - this%origin + roundingEpsilon,&
          & this%invBasis))
    end do

  end subroutine TGrid_cartesianToGridCoord


  !> Transformation from cartesian coordinates to real grid coordinates
  subroutine TGrid_cartesianToRealGridCoord(this, cartcoords, gridcoords)

    !> representation of parallelepiped grid
    class(TGrid), intent(in) :: this

    !> cartesian coordinates, shape: [3, nCoords]
    real(dp), intent(in) :: cartcoords(:,:)

    !> corresponding real grid coordinates, shape: [3, nCoords]
    real(dp), intent(out), allocatable :: gridcoords(:,:)

    !> auxiliary variable
    integer :: ii

    allocate(gridcoords(size(cartcoords, dim=1), size(cartcoords, dim=2)))

    do ii = 1, size(cartcoords, dim=2)
      gridcoords(:, ii) = matmul(this%invBasis, cartcoords(:, ii) - this%origin)
    end do

  end subroutine TGrid_cartesianToRealGridCoord


  subroutine TGrid_getCorners(this, cornCoords)

    !> representation of parallelepiped grid
    class(TGrid), intent(in) :: this

    !> integer grid coordinates of the 2**3=8 corners of the parallelepipedon
    integer, intent(out) :: cornCoords(3,8)

    !> indices of corners, shape: [3, 8]
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
      cornCoords(:, jj) = tmp(:)
    end do

  end subroutine TGrid_getCorners


  subroutine getGridCoords(grid, gridcoords)

    !> representation of parallelepiped grid
    type(TGrid), intent(in) :: grid

    !> integer grid coordinates
    integer, intent(out), allocatable :: gridCoords(:,:)

    !> auxiliary variables
    integer :: ii, jj, kk, ind

    allocate(gridCoords(3,&
        & (grid%upperBounds(1) - grid%lowerBounds(1) + 1) *&
        & (grid%upperBounds(2) - grid%lowerBounds(2) + 1) *&
        & (grid%upperBounds(3) - grid%lowerBounds(3) + 1)))

    ind = 1

    do kk = grid%lowerBounds(3), grid%upperBounds(3)
      do jj = grid%lowerBounds(2), grid%upperBounds(2)
        do ii = grid%lowerBounds(1), grid%upperBounds(1)
          gridCoords(1, ind) = ii
          gridCoords(2, ind) = jj
          gridCoords(3, ind) = kk
          ind = ind + 1
        end do
      end do
    end do

  end subroutine getGridCoords


  !> Calculates lower and upper bounds of the intersection of the total and the subgrid in integer
  !> total grid coordinates and cartesian subgrid coordinates.
  subroutine TGrid_getSubgridRanges(this, subgrid, position, luGc, luSubGc, intersection)

    !> representation of parallelepiped grid
    class(TGrid), intent(in) :: this

    !> subgrid to calculate intersection with
    class(TGrid), intent(in) :: subgrid

    !> position where to place subgrid onto grid, expected shape: [3]
    real(dp), intent(in) :: position(:)

    !> lower and upper bounds of intersection in integer grid coordinates, shape: [3,2]
    integer, intent(out), allocatable :: luGc(:,:)

    !> lower and upper bounds of intersection in real subgrid coordinates, shape: [3,2]
    real(dp), intent(out), allocatable :: luSubGc(:,:)

    !> true if there is an intersection between both grids, false otherwise
    logical, intent(out) :: intersection

    !> lower and upper bounds of intersection in real grid coordinates, shape: [3,2]
    real(dp), allocatable :: realLuGc(:,:)

    !> lower and upper bounds of intersection in real subgrid coordinates, shape: [3,2]
    real(dp), allocatable :: realLuSubGc(:,:)

    !> lower and upper bounds of intersection in real grid coordinates, shape: [3,2]
    real(dp) :: luRGc(3,2)

    !> lower and upper bounds of intersection in cartesian coordinates, shape: [3,2]
    real(dp) :: lCc(3, 2), uCc(3,2)

    !> auxiliary array for containing cartesian coordinates of processed intersection, shape: [3,2]
    real(dp), allocatable :: roundedLuCc(:,:)

    !> lower and upper cartesian bounds of grid and (partly) contained subgrid, shape: [3,1]
    real(dp), allocatable :: lCcTotGrid(:,:), uCcTotGrid(:,:), lCcSubGrid(:,:), uCcSubGrid(:,:)

    !> lower and upper real grid coordinates of intersection, shape: [3,1]
    real(dp), allocatable :: lRGc(:,:), uRGc(:,:)

    !> auxiliary variable
    integer :: ii

    allocate(luGc(3,2))

    call this%gridcoordToCartesian(reshape(this%lowerBounds, [3,1]), lCcTotGrid)
    call this%gridcoordToCartesian(reshape(this%upperBounds, [3,1]), uCcTotGrid)
    call subgrid%gridcoordToCartesian(reshape(subgrid%lowerBounds, [3,1]), lCcSubGrid)
    call subgrid%gridcoordToCartesian(reshape(subgrid%upperBounds, [3,1]), uCcSubGrid)

    lCc(:,1) = lCcSubGrid(:,1) + position
    uCc(:,1) = uCcSubGrid(:,1) + position

    call this%cartesianToRealGridcoord(lCc, lRGc)
    call this%cartesianToRealGridcoord(uCc, uRGc)

    luRGc(:,1) = max(lRGc(:,1), real(this%lowerBounds, dp))
    luRGc(:,2) = min(uRGc(:,1), real(this%upperBounds, dp))

    luGc(:,:) = floor(luRGc(:,:))

    luGc(:,1) = max(luGc(:,1), this%lowerBounds)
    luGc(:,1) = min(luGc(:,1), this%upperBounds)
    luGc(:,2) = max(luGc(:,2), this%lowerBounds)
    luGc(:,2) = min(luGc(:,2), this%upperBounds)

    call this%gridcoordToCartesian(luGc, roundedLuCc)
    call subgrid%cartesianToRealGridcoord(roundedLuCc - spread(position, 2, 2), luSubGc)

    ! Check if there is an intersection of both grids
    if ((any(min(uCc(:,1), uCcTotGrid(:,1)) - max(lCc(:,1), lCcTotGrid(:,1)) .lt. 0.0)) .or.&
        & any(luRGc(:,2) - luRGc(:,1) .lt. 0)) then
      intersection = .false.
    else
      intersection = .true.
    end if

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

    !> lower and upper cartesian bounds of grid and (partly) contained subgrid, shape: [3,1]
    real(dp), allocatable :: lCcTotGrid(:,:), uCcTotGrid(:,:), lCcSubGrid(:,:), uCcSubGrid(:,:)

    isInside = .true.

    call this%gridcoordToCartesian(reshape(this%lowerBounds, [3,1]), lCcTotGrid)
    call this%gridcoordToCartesian(reshape(this%upperBounds, [3,1]), uCcTotGrid)
    call subgrid%gridcoordToCartesian(reshape(subgrid%lowerBounds, [3,1]), lCcSubGrid)
    call subgrid%gridcoordToCartesian(reshape(subgrid%upperBounds, [3,1]), uCcSubGrid)

    if (any(lCcSubGrid(:,1) + position < lCcTotGrid(:,1))) then
      isInside = .false.
    elseif (any(uCcSubGrid(:,1) + position > uCcTotGrid(:,1))) then
      isInside = .false.
    end if

  end subroutine TGrid_hasSubgrid


  subroutine TGrid_hasGridCoords(this, gridcoords, isInside)

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

  end subroutine TGrid_hasGridCoords


  !> Calculates the cartesian distances between a reference position and all grid points of a given
  !> grid
  subroutine getDistances(grid, position, distances)

    !> representation of a grid
    type(TGrid), intent(in) :: grid

    !> reference point for calculating the distances, shape: [3]
    real(dp), intent(in) :: position(:)

    !> contains calculated distances on exit, shape: [nxPoints, nyPoints, nzPoints]
    real(dp), intent(out), allocatable :: distances(:,:,:)

    !> cartesian coordinates of current grid point, shape: [3,1]
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


  !> Set the subroutine wich should be used for the tabulation of the rwfs.
  subroutine TGridData_setRwTabulationType(this, interType)

    !> representation of volumetric grid and data
    class(TGridData), intent(inout) :: this

    !> interpolation type to use for obtaining radial wavefunctions at arbitrary distances
    integer, intent(in) :: interType

    select case (interType)
    case (rwTabulationTypes%trivial)
      this%tabulateBasis => TGridData_tabulateBasisTrivial
    case (rwTabulationTypes%linear)
      this%tabulateBasis => TGridData_tabulateBasisLinear
    case (rwTabulationTypes%polynomial)
      this%tabulateBasis => TGridData_tabulateBasisPolynomial
    case (rwTabulationTypes%explicit)
      this%tabulateBasis => TGridData_tabulateBasisExplicit
    case default
      this%tabulateBasis => TGridData_tabulateBasisExplicit
  end select

  end subroutine TGridData_setRwTabulationType


  !> Sets the origin of a grid instance in cartesian coordinates
  subroutine TGridData_setGridOrigin(this, origin)

    !> representation of volumetric grid and data
    class(TGridData), intent(inout) :: this

    !> cartesian grid coordinates of grid origin, shape: [3]
    real(dp), intent(in) :: origin(:)

    call this%grid%setOrigin(origin)

  end subroutine TGridData_setGridOrigin


  !> Calculates the data at the eight corners of the subcube around a subgrid point.
  subroutine TGridData_getSubcubeCornsAndVals(this, gridcoords, data)

    !> representation of volumetric grid and data
    class(TGridData), intent(in) :: this

    !> integer subgrid grid coordinates to calculate subcubes for
    !> expected shape: [3, nxIntersecPoints, nyIntersecPoints, nzIntersecPoints]
    integer, intent(in) :: gridCoords(:,:,:,:)

    !> data corresponding to the eight corners of the subcube
    !> expected shape: [2, 2, 2, nxIntersecPoints, nyIntersecPoints, nzIntersecPoints]
    real(dp), intent(out), allocatable :: data(:,:,:,:,:,:)

    !> data indices, corresponding to an integer grid point
    integer :: inds(3)

    !> auxiliary variables
    integer :: ii, jj, kk

    allocate(data(2, 2, 2, size(gridCoords, dim=2), size(gridCoords, dim=3),&
        & size(gridCoords, dim=4)))

    do kk = 1, size(gridCoords, dim=4)
      do jj = 1, size(gridCoords, dim=3)
        do ii = 1, size(gridCoords, dim=2)
          inds(:) = gridCoords(:, ii, jj, kk) - this%grid%lowerBounds + 1
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


  !> Tabulates the slater type orbitals, with an explicit calculation of the sto values.
  subroutine TGridData_tabulateBasisExplicit(this, rty, ll, mm, rwf, alpha, aa)

    !> representation of volumetric grid and data
    class(TGridData), intent(inout) :: this

    !> tabulated real tesseral spherical harmonics for all l,m-combinations
    class(TRealTessY), intent(in) :: rty

    !> angular momentum of the species
    integer, intent(in) :: ll

    !> magnetic quantum number
    integer, intent(in) :: mm

    !> tabulated radial wavefunction, shape: [nDistances, 3], (unused in this routine)
    real(dp), intent(in) :: rwf(:,:)

    !> exponential coefficients, shape: [nAlpha]
    real(dp), intent(in) :: alpha(:)

    !> summation coefficients, shape: [nCoeffPerAlpha, nAlpha]
    real(dp), intent(in) :: aa(:,:)

    !> cartesian distances of grid points to origin, shape depends on intersection of the grids,
    !> general shape: [nxIntersecPoints, nyIntersecPoints, nzIntersecPoints]
    real(dp), allocatable :: distances(:,:,:)

    !> radial part of the slater type orbital
    real(dp) :: sto_rwf

    !> data indices, corresponding to a grid point
    integer :: datainds(3)

    !> auxiliary variables
    integer :: ii, jj, kk, rtyIdx, idx(3)
    integer :: nAlpha
    integer :: nPow

    nAlpha = size(alpha)
    nPow = size(aa, dim=1)

    rtyIdx = ll * (ll + 1) + 1 + mm

    call getDistances(this%grid, [0.0_dp, 0.0_dp, 0.0_dp], distances)

    do kk = this%grid%lowerBounds(3), this%grid%upperBounds(3)
      idx(3) = kk
      do jj = this%grid%lowerBounds(2), this%grid%upperBounds(2)
        idx(2) = jj
        do ii = this%grid%lowerBounds(1), this%grid%upperBounds(1)
          idx(1) = ii
          datainds(:) = idx - this%grid%lowerBounds + 1
          call SlaterOrbital_getValue_explicit(ll, nPow, nAlpha, aa, alpha, &
              & distances(datainds(1), datainds(2), datainds(3)), sto_rwf)
          this%data(datainds(1), datainds(2), datainds(3)) =&
              & this%data(datainds(1), datainds(2), datainds(3))&
              & + sto_rwf * rty%data(datainds(1), datainds(2), datainds(3), rtyIdx)
        end do
      end do
    end do

  end subroutine TGridData_tabulateBasisExplicit


  !> Tabulates the slater type orbitals, with a trivial interpolation of the radial WF values.
  subroutine TGridData_tabulateBasisTrivial(this, rty, ll, mm, rwf, alpha, aa)

    !> representation of volumetric grid and data
    class(TGridData), intent(inout) :: this

    !> tabulated real tesseral spherical harmonics for all l,m-combinations
    class(TRealTessY), intent(in) :: rty

    !> angular momentum of the species
    integer, intent(in) :: ll

    !> magnetic quantum number
    integer, intent(in) :: mm

    !> tabulated radial wavefunction, shape: [nDistances, 3]
    real(dp), intent(in) :: rwf(:,:)

    !> exponential coefficients, (unused in this routine)
    real(dp), intent(in) :: alpha(:)

    !> summation coefficients, shape: [nCoeffPerAlpha, nAlpha], (unused in this routine)
    real(dp), intent(in) :: aa(:,:)

    !> cartesian distances of grid points to origin, shape depends on intersection of the grids,
    !> general shape: [nxIntersecPoints, nyIntersecPoints, nzIntersecPoints]
    real(dp), allocatable :: distances(:,:,:)

    !> data indices, corresponding to a grid point
    integer :: datainds(3)

    !> index of nearest distance
    integer :: loc

    !> auxiliary variables
    integer :: ii, jj, kk, rtyIdx, endIdx, idx(3)

    rtyIdx = ll * (ll + 1) + 1 + mm
    endIdx = size(rwf, dim=1)

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


  !> Tabulates the slater type orbitals, with a linear interpolation of the radial WF values.
  subroutine TGridData_tabulateBasisLinear(this, rty, ll, mm, rwf, alpha, aa)

    !> representation of volumetric grid and data
    class(TGridData), intent(inout) :: this

    !> tabulated real tesseral spherical harmonics for all l,m-combinations
    class(TRealTessY), intent(in) :: rty

    !> angular momentum of the species
    integer, intent(in) :: ll

    !> magnetic quantum number
    integer, intent(in) :: mm

    !> tabulated radial wavefunction, shape: [nDistances, 3]
    real(dp), intent(in) :: rwf(:,:)

    !> exponential coefficients, (unused in this routine)
    real(dp), intent(in) :: alpha(:)

    !> summation coefficients, shape: [nCoeffPerAlpha, nAlpha], (unused in this routine)
    real(dp), intent(in) :: aa(:,:)

    !> cartesian distances of grid points to origin, shape depends on intersection of the grids,
    !> general shape: [nxIntersecPoints, nyIntersecPoints, nzIntersecPoints]
    real(dp), allocatable :: distances(:,:,:)

    !> data indices, corresponding to a grid point
    integer :: datainds(3)

    !> index of nearest distance
    integer :: loc

    !> auxiliary variables
    integer :: ii, jj, kk, rtyIdx, endIdx, idx(3)

    rtyIdx = ll * (ll + 1) + 1 + mm
    endIdx = size(rwf, dim=1)

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


  !> Tabulates the slater type orbitals, with a polynomial interpolation of the radial WF values.
  subroutine TGridData_tabulateBasisPolynomial(this, rty, ll, mm, rwf, alpha, aa)

    !> representation of volumetric grid and data
    class(TGridData), intent(inout) :: this

    !> tabulated real tesseral spherical harmonics for all l,m-combinations
    class(TRealTessY), intent(in) :: rty

    !> angular momentum of the species
    integer, intent(in) :: ll

    !> magnetic quantum number
    integer, intent(in) :: mm

    !> tabulated radial wavefunction, shape: [nDistances, 3]
    real(dp), intent(in) :: rwf(:,:)

    !> exponential coefficients, (unused in this routine)
    real(dp), intent(in) :: alpha(:)

    !> summation coefficients, shape: [nCoeffPerAlpha, nAlpha], (unused in this routine)
    real(dp), intent(in) :: aa(:,:)

    !> cartesian distances of grid points to origin, shape depends on intersection of the grids,
    !> general shape: [nxIntersecPoints, nyIntersecPoints, nzIntersecPoints]
    real(dp), allocatable :: distances(:,:,:)

    !> data indices, corresponding to a grid point
    integer :: datainds(3)

    !> auxiliary variables
    integer :: ii, jj, kk, rtyIdx, idx(3)

    rtyIdx = ll * (ll + 1) + 1 + mm

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
              & polynomialInterpolation(rwf(:, 1), rwf(:, 2), rwf(:, 3), distances(datainds(1),&
              & datainds(2), datainds(3))) * rty%data(datainds(1), datainds(2), datainds(3), rtyIdx)
        end do
      end do
    end do

  end subroutine TGridData_tabulateBasisPolynomial


  !> Calculates real grid coordinates of the subgrid points which intersect with the total grid.
  subroutine TGrid_getIntersecGridPoints(this, subgrid, luGc, luSubGc, subgridGc)

    !> representation of parallelepiped grid
    class(TGrid), intent(in) :: this

    !> subgrid to calculate intersection with
    class(TGrid), intent(in) :: subgrid

    !> lower and upper bounds of intersection in integer grid coordinates, shape: [3,2]
    integer, intent(in) :: luGc(:,:)

    !> lower and upper bounds of intersection in real subgrid coordinates, shape: [3,2]
    real(dp), intent(in) :: luSubGc(:,:)

    !> real subgrid grid coordinates that intersect with grid,
    !> shape: [3, nxPoints, nyPoints, nzPoints]
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


  subroutine TGridData_getInterpolatedValueGc(this, gridCoords, data)

    !> representation of volumetric grid and data
    class(TGridData), intent(in) :: this

    !> real subgrid coordinates to calculate interpolation for
    real(dp), intent(in) :: gridCoords(:,:,:,:)

    !> corresponding (interpolated) data
    real(dp), intent(out), allocatable :: data(:,:,:)

    call this%getValue(floor(gridCoords + roundingEpsilon), data)

  end subroutine TGridData_getInterpolatedValueGc


  subroutine TGridData_getValue(this, gridcoords, vals)

    !> representation of volumetric grid and data
    class(TGridData), intent(in) :: this

    !> integer grid coordinates to return data for
    integer, intent(in) :: gridCoords(:,:,:,:)

    !> data points corresponding to integer grid coordinates
    real(dp), intent(out), allocatable :: vals(:,:,:)

    !> data indices, specifying a data access point
    integer :: inds(3)

    !> auxiliary variables
    integer :: ii, jj, kk

    allocate(vals(size(gridCoords, dim=2), size(gridCoords, dim=3), size(gridCoords, dim=4)))

    do kk = 1, size(gridCoords, dim=4)
      do jj = 1, size(gridCoords, dim=3)
        do ii = 1, size(gridCoords, dim=2)
          inds(:) = gridCoords(:, ii, jj, kk) - this%grid%lowerBounds + 1
          vals(ii, jj, kk) = this%data(inds(1), inds(2), inds(3))
        end do
      end do
    end do

  end subroutine TGridData_getValue


  !> Collects and adds all contributions from the subgrids onto the global grid weighted with the
  !> coefficients.
  subroutine subgridsToUncachedGlobalGrid(globalGridDat, subgridsDat, atomicGridDat, coords, &
      & coeffs, orbitalToAtom, orbitalToSubgrid, nRegions, stos, orbitalToAngs, orbitalToM, &
      & orbitalToStos, gridInterType, addDensities)

    !> representation of volumetric grid and data
    type(TGridData), intent(in) :: globalGridDat

    !> subgrids to add to global grid instance, exp. shape: [nOrbitals]
    type(TGridData), intent(in) :: subgridsDat(:)

    !> array containing values of the atomic charge, shape: [nxPoints, nyPoints, nzPoints]
    real(dp), allocatable, intent(out) :: atomicGridDat(:,:,:)

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

    !> Tabulated radial WF of the slater type orbitals
    type(TSlaterOrbital), intent(in) :: stos(:)

    !> Index mapping from orbitals to angular momentum quantum numbers
    integer, intent(in) :: orbitalToAngs(:)

    !> Index mapping from orbitals to magnetic quantum numbers
    integer, intent(in) :: orbitalToM(:)

    !> Index mapping from orbitals to slater type orbital
    integer, intent(in) :: orbitalToStos(:)

    !> Interpolation type for the grid interpolation
    integer, intent(in) :: gridInterType

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
    real(dp), pointer :: pTmpData(:,:,:)

    !> auxiliary variables
    integer :: iGrid, iOrb, iAtom, iCell, iSubgrid, iAng, iM, idx(3,1), ii, jj, kk, aa(3), iStos

    !> Temporary stored cartesian coordinates of one grid point, shape: [3,1]
    real(dp), allocatable :: cartcoordsTmp(:,:)

    !> Cartesian coordinates of the total grid points, shape: [3, nxPoints, nyPoints, nzPoints]
    real(dp), allocatable :: cartCoords(:,:,:,:)

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
      tmpGrid%lowerBounds(3) = tiling(1, iGrid) - 1
      tmpGrid%upperBounds(3) = tiling(2, iGrid) - 1
      tmpGrid%nPoints(3) = tiling(2, iGrid) - tiling(1, iGrid) + 1
      pTmpData => globalGridDat%data(:,:, tiling(1, iGrid):tiling(2, iGrid))
      call TGridData_init(regionGridDat(iGrid), tmpGrid, pTmpData,&
          & rwTabulationType=rwTabulationTypes%explicit)
    end do

    allocate(atomicGridDat(globalGridDat%grid%nPoints(1), globalGridDat%grid%nPoints(2), &
        & globalGridDat%grid%nPoints(3)))
    atomicGridDat(:,:,:) = 0.0_dp

    ! Calculate the cartesian coordinates of the grid points of the total grid.
    if (gridInterType == gridInterpolTypes%explicit) then
      allocate(cartCoords(3, globalGridDat%grid%nPoints(1), &
          & globalGridDat%grid%nPoints(2), globalGridDat%grid%nPoints(3)))

      do iGrid = 1, nRegions
        do kk = regionGridDat(iGrid)%grid%lowerBounds(3), regionGridDat(iGrid)%grid%upperBounds(3)
          idx(3, 1) = kk
          do jj = regionGridDat(iGrid)%grid%lowerBounds(2), regionGridDat(iGrid)%grid%upperBounds(2)
            idx(2, 1) = jj
            do ii = regionGridDat(iGrid)%grid%lowerBounds(1),&
                  & regionGridDat(iGrid)%grid%upperBounds(1)
              idx(1, 1) = ii
              aa(:) = idx(:, 1) - regionGridDat(iGrid)%grid%lowerBounds + 1
              call regionGridDat(iGrid)%grid%gridcoordToCartesian(idx, cartCoordsTmp)
              cartCoords(:, aa(1), aa(2), tiling(1, iGrid) + aa(3) - 1) = cartcoordsTmp(:,1)
            end do
          end do
        end do
      end do
    end if

    !$omp parallel do default(none)&
    !$omp& shared(nRegions, regionGridDat, subgridsDat, orbitalToAtom, orbitalToSubgrid,&
    !$omp& coords, coeffs, tAddDensities, atomicGridDat, gridInterType, tiling, stos,&
    !$omp& orbitalToM, cartCoords, orbitalToStos, orbitalToAngs)&
    !$omp& private(iGrid, iOrb, iAtom, iSubgrid, iCell, iAng, iM, idx, ii, jj, kk, cartCoordsTmp,&
    !$omp& aa, iStos)
    lpGrid: do iGrid = 1, nRegions
      lpOrb: do iOrb = 1, size(coeffs)

        iAtom = orbitalToAtom(iOrb)
        iSubgrid = orbitalToSubgrid(iOrb)
        iAng = orbitalToAngs(iOrb)
        iM = orbitalToM(iOrb)
        iStos = orbitalToStos(iOrb)

        lpCell: do iCell = 1, size(coords, dim=3)

          if (gridInterType == gridInterpolTypes%explicit) then
            atomicGridDat(:,:, tiling(1, iGrid):tiling(2, iGrid)) = &
                & atomicGridDat(:,:, tiling(1, iGrid):tiling(2, iGrid)) + &
                & explicitSTOcalculation(regionGridDat(iGrid), subgridsDat(iSubgrid), &
                & coeffs(iOrb), coords(:, iAtom, iCell), &
                & cartCoords(:,:,:,tiling(1, iGrid):tiling(2, iGrid)), &
                & stos(iStos), iAng, iM, square=tAddDensities)
          else if (gridInterType .eq. gridInterpolTypes%trivial) then
            atomicGridDat(:,:, tiling(1, iGrid):tiling(2, iGrid)) = &
                & atomicGridDat(:,:, tiling(1, iGrid):tiling(2, iGrid)) + &
                & gridInterpolationTrivial(regionGridDat(iGrid), subgridsDat(iSubgrid), &
                & coeffs(iOrb), coords(:, iAtom, iCell), square=tAddDensities)
          else if (gridInterType .eq. gridInterpolTypes%linear) then
            atomicGridDat(:,:, tiling(1, iGrid):tiling(2, iGrid)) = &
                & atomicGridDat(:,:, tiling(1, iGrid):tiling(2, iGrid)) + &
                & gridInterpolationLinear(regionGridDat(iGrid), subgridsDat(iSubgrid), &
                & coeffs(iOrb), coords(:, iAtom, iCell), square=tAddDensities)
          end if

        end do lpCell

      end do lpOrb
    end do lpGrid
    !$omp end parallel do

  end subroutine subgridsToUncachedGlobalGrid


  !> Collects and add all contributions from the subgrids onto the global grid weighted with the
  !> coefficients.
  subroutine subgridsToCachedGlobalGrids(globalGridDat, subgridsDat, coords, coeffs, levelIndex,&
      & requiredLevels, requiredSpins, orbitalToAtom, orbitalToSubgrid, nRegions, regionGridDat, &
      & cartCoords, tiling, statesTiling, stos, orbitalToAngs, orbitalToM, orbitalToStos, &
      & globalGridsDat, gridInterType, addDensities)

    !> prototype of global grid instance
    type(TGridData), intent(in) :: globalGridDat

    !> subgrids to add to global grid instances, exp. shape: [nOrbitals]
    type(TGridData), intent(in) :: subgridsDat(:)

    !> coordinates, defining the cartesian positions of the subgrids, exp. shape: [3, nAtom, nCell]
    real(dp), intent(in) :: coords(:,:,:)

    !> orbital-resolved coefficients for volumetric data, shape: [nOrbital, nLevel, nKPoint, nSpin]
    real(dp), intent(in) :: coeffs(:,:,:,:)

    !> Index Mapping from state index to [Level, KPoint, Spin]
    integer, intent(in) :: levelIndex(:,:)

    !> All levels which has to be calculated
    integer, intent(in) :: requiredLevels(:)

    !> All spins which has to be calculated
    integer, intent(in) :: requiredSpins(:)

    !> index mapping from orbital to atom, exp. shape: [nOrbitals]
    integer, intent(in) :: orbitalToAtom(:)

    !> index mapping from orbital to subgrid, exp. shape: [nOrbitals]
    integer, intent(in) :: orbitalToSubgrid(:)

    !> number of parallel global regions
    integer, intent(in) :: nRegions

    !> grid parallel regions, holding global grid sections
    type(TGridData), intent(in) :: regionGridDat(:)

    !> Cartesian coordinates of the total grid points, shape: [3, nxPoints, nyPoints, nzPoints]
    real(dp), intent(in) :: cartCoords(:,:,:,:)

    !> start and end indices of all parallel tiles, exp. shape: [2, nRegions]
    integer, intent(in) :: tiling(:,:)

    !> start and end index of states which need to be calculated
    integer, intent(in) :: statesTiling(:)

    !> Tabulated radial WF of the slater type orbitals
    type(TSlaterOrbital), intent(in) :: stos(:)

    !> Index mapping from orbitals to angular momentum quantum numbers
    integer, intent(in) :: orbitalToAngs(:)

    !> Index mapping from orbitals to magnetic quantum numbers
    integer, intent(in) :: orbitalToM(:)

    !> Index mapping from orbitals to slater type orbital
    integer, intent(in) :: orbitalToStos(:)

    !> array holding the data of the grid for every level which is needed,
    !> shape: [nxPoints, nyPoints, nzPoints, nLevel]
    real(dp), intent(out), allocatable :: globalGridsDat(:,:,:,:,:)

    !> Interpolation type for the grid interpolation
    integer, intent(in) :: gridInterType

    !> add densities instead of wave functions
    logical, intent(in), optional :: addDensities

    !> true, if optional argument provided and .true., otherwise False
    logical :: tAddDensities

    !> array holding the values of the STOs for every periodic Cell,
    !> shape: [nxPoints, nyPoints, nzPoints, nCells]
    real(dp), allocatable :: basis(:,:,:)

    !> auxiliary variables
    integer :: iGrid, iOrb, iCell, iAtom, iSubgrid, iSpin, iL, iAng, iM, iStos, ind
    integer :: currentLevel, currentSpin
    integer :: levelInd

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

    allocate(globalGridsDat(globalGridDat%grid%nPoints(1), globalGridDat%grid%nPoints(2), &
        & globalGridDat%grid%nPoints(3), statesTiling(2) - statesTiling(1) + 1, &
        & size(requiredSpins)))
    globalGridsDat(:,:,:,:,:) = 0.0_dp

    allocate(basis(globalGridDat%grid%nPoints(1), globalGridDat%grid%nPoints(2), &
        & globalGridDat%grid%nPoints(3)))
    basis(:,:,:) = 0.0_dp

    !$omp parallel do default(none)&
    !$omp& shared(nRegions, regionGridDat, subgridsDat, orbitalToAtom, orbitalToSubgrid, stos,&
    !$omp& orbitalToAngs, orbitalToM, cartCoords, orbitalToStos, coords, coeffs, levelIndex, &
    !$omp& globalGridsDat, tiling, tAddDensities, gridInterType, basis, statesTiling,&
    !$omp& requiredLevels, requiredSpins)&
    !$omp& private(iGrid, iOrb, iAtom, iSubgrid, iCell, iAng, iM, iStos, ind, currentLevel, &
    !$omp& currentSpin, iL, iSpin, levelInd)
    lpGrid: do iGrid = 1, nRegions
      lpOrb: do iOrb = 1, size(coeffs, dim=1)

        iAtom = orbitalToAtom(iOrb)
        iSubgrid = orbitalToSubgrid(iOrb)
        iAng = orbitalToAngs(iOrb)
        iM = orbitalToM(iOrb)
        iStos = orbitalToStos(iOrb)

        lpCell: do iCell = 1, size(coords, dim=3)

          if (gridInterType == gridInterpolTypes%explicit) then
            basis(:,:, tiling(1, iGrid):tiling(2, iGrid)) = &
                & basis(:,:, tiling(1, iGrid):tiling(2, iGrid)) + &
                & explicitSTOcalculation(regionGridDat(iGrid), subgridsDat(iSubgrid), &
                & 1.0_dp, coords(:, iAtom, iCell), &
                & cartCoords(:,:,:,tiling(1, iGrid):tiling(2, iGrid)), stos(iStos), iAng, iM, &
                & square=tAddDensities)
          else if (gridInterType .eq. gridInterpolTypes%trivial) then
            basis(:,:, tiling(1, iGrid):tiling(2, iGrid)) = &
                & basis(:,:, tiling(1, iGrid):tiling(2, iGrid)) + &
                & gridInterpolationTrivial(regionGridDat(iGrid), subgridsDat(iSubgrid), &
                & 1.0_dp, coords(:, iAtom, iCell), square=tAddDensities)
          else if (gridInterType .eq. gridInterpolTypes%linear) then
            basis(:,:, tiling(1, iGrid):tiling(2, iGrid)) = &
                & basis(:,:, tiling(1, iGrid):tiling(2, iGrid)) + &
                & gridInterpolationLinear(regionGridDat(iGrid), subgridsDat(iSubgrid), &
                & 1.0_dp, coords(:, iAtom, iCell), square=tAddDensities)
          end if

        end do lpCell

        levelInd = 1
        do iL = statesTiling(1), statesTiling(2)
          currentLevel = requiredLevels(iL)

            lpSpin: do iSpin = 1, size(requiredSpins)
              currentSpin = requiredSpins(iSpin)

              globalGridsDat(:,:, tiling(1, iGrid):tiling(2, iGrid), levelInd, iSpin) =&
                  & globalGridsDat(:,:, tiling(1, iGrid):tiling(2, iGrid), levelInd, iSpin)&
                  & + coeffs(iOrb, currentLevel, 1, currentSpin) * &
                  & basis(:,:, tiling(1, iGrid):tiling(2, iGrid))

            end do lpSpin

          levelInd = levelInd + 1

        end do

        basis(:,:, tiling(1, iGrid):tiling(2, iGrid)) = 0.0_dp

      end do lpOrb
    end do lpGrid
    !$omp end parallel do

  end subroutine subgridsToCachedGlobalGrids


  !> Collects and adds all contributions from the subgrids onto the global grid weighted with the
  !> coefficients.
  subroutine subgridsToCachedGlobalGridsCplx(globalGridDat, subgridsDat, coords, coeffs,&
      & levelIndex, requiredLevels, requiredKPoints, requiredSpins, requiredKPointsForLevel, &
      & orbitalToAtom, orbitalToSubgrid, nRegions, regionGridDat, cartCoords, tiling, &
      & statesTiling, stos, orbitalToAngs, orbitalToM, orbitalToStos, globalGridsDat,&
      & gridInterType, phases, addDensities)

    !> prototype of global grid instance
    type(TGridData), intent(in) :: globalGridDat

    !> subgrids to add to global grid instances, exp. shape: [nOrbitals]
    type(TGridData), intent(in) :: subgridsDat(:)

    !> coordinates, defining the cartesian positions of the subgrids, exp. shape: [3, nAtom, nCell]
    real(dp), intent(in) :: coords(:,:,:)

    !> orbital-resolved coefficients for volumetric data,
    !> shape: [totOrbital, totLevel, totKPoint, totSpin]
    complex(dp), intent(in) :: coeffs(:,:,:,:)

    !> Index Mapping from state index to [Level, KPoint, Spin]
    integer, intent(in) :: levelIndex(:,:)

    !> All levels which has to be calculated
    integer, intent(in) :: requiredLevels(:)

    !> All k-points which has to be calculated
    integer, intent(in) :: requiredKPoints(:)

    !> All Spins which has to be calculated
    integer, intent(in) :: requiredSpins(:)

    !> The required k-points for each level. Shape: [nRequiredLevels, nRequiredKPoints]
    integer, intent(in) :: requiredKPointsForLevel(:,:)

    !> index mapping from orbital to atom, exp. shape: [nOrbitals]
    integer, intent(in) :: orbitalToAtom(:)

    !> index mapping from orbital to subgrid, exp. shape: [nOrbitals]
    integer, intent(in) :: orbitalToSubgrid(:)

    !> number of parallel global regions
    integer, intent(in) :: nRegions

    !> grid parallel regions, holding global grid sections
    type(TGridData), intent(in) :: regionGridDat(:)

    !> Cartesian coordinates of the total grid points, shape: [3, nxPoints, nyPoints, nzPoints]
    real(dp), intent(in) :: cartCoords(:,:,:,:)

    !> start and end indices of all parallel tiles, exp. shape: [2, nRegions]
    integer, intent(in) :: tiling(:,:)

    !> start and end index of states which need to be calculated
    integer, intent(in) :: statesTiling(:)

    !> Tabulated radial WF of the slater type orbitals
    type(TSlaterOrbital), intent(in) :: stos(:)

    !> Index mapping from orbitals to angular momentum quantum numbers
    integer, intent(in) :: orbitalToAngs(:)

    !> Index mapping from orbitals to magnetic quantum numbers
    integer, intent(in) :: orbitalToM(:)

    !> Index mapping from orbitals to slater type orbital
    integer, intent(in) :: orbitalToStos(:)

    !> array holding the data of the grid for every level which is needed,
    !> shape: [nxPoints, nyPoints, nzPoints, nLevel, nKPoints, nSpin]
    complex(dp), intent(out), allocatable :: globalGridsDat(:,:,:,:,:,:)

    !> Interpolation type for the grid interpolation
    integer, intent(in) :: gridInterType

    !> Complex phase factor. It has the shape: [nCell, nKPoints]
    complex(dp), intent(in) :: phases(:,:)

    !> add densities instead of wave functions
    logical, intent(in), optional :: addDensities

    !> true, if optional argument provided and .true., otherwise False
    logical :: tAddDensities

    !> array holding the values of the STOs for every periodic Cell,
    !> shape: [nxPoints, nyPoints, nzPoints, nCells]
    complex(dp), allocatable :: basis(:,:,:,:)
    complex(dp), allocatable :: basisTmp(:,:,:)

    !> auxiliary variables
    integer :: iGrid, iOrb, iCell, iAtom, iSubgrid, iSpin, iKPoint, iL, iAng, iM, &
        & iStos, ind, levelInd
    integer :: currentKPoint, currentLevel, currentSpin

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

    allocate(globalGridsDat(globalGridDat%grid%nPoints(1), globalGridDat%grid%nPoints(2), &
        & globalGridDat%grid%nPoints(3), statesTiling(2) - statesTiling(1) + 1, &
        & size(requiredKPoints), size(requiredSpins)))
    globalGridsDat(:,:,:,:,:,:) = 0.0_dp

    !$omp parallel do default(none)&
    !$omp& shared(nRegions, regionGridDat, subgridsDat, orbitalToAtom, orbitalToSubgrid,&
    !$omp& coords, coeffs, globalGridsDat, tiling, tAddDensities, gridInterType,&
    !$omp& phases, stos, orbitalToAngs, orbitalToM, cartcoords, orbitalToStos, statesTiling, &
    !$omp& requiredLevels, requiredKPoints, requiredSpins, requiredKPointsForLevel, globalGridDat)&
    !$omp& private(iGrid, iOrb, iAtom, iSubgrid, iCell, ind, iAng, iM, iStos, basisTmp,&
    !$omp& currentKPoint, currentLevel, currentSpin, iL, iSpin, iKPoint, levelInd, basis)
    lpGrid: do iGrid = 1, nRegions
      allocate(basisTmp(globalGridDat%grid%nPoints(1), globalGridDat%grid%nPoints(2), &
          & tiling(2, iGrid) - tiling(1, iGrid) + 1))
      basisTmp(:,:,:) = 0.0_dp
      allocate(basis(globalGridDat%grid%nPoints(1), globalGridDat%grid%nPoints(2), &
          & tiling(2, iGrid) - tiling(1, iGrid) + 1, size(requiredKPoints)))
      basis(:,:,:,:) = 0.0_dp

      lpOrb: do iOrb = 1, size(coeffs, dim=1)
        basis(:,:,:,:) = 0.0_dp

        iAtom = orbitalToAtom(iOrb)
        iSubgrid = orbitalToSubgrid(iOrb)
        iAng = orbitalToAngs(iOrb)
        iM = orbitalToM(iOrb)
        iStos = orbitalToStos(iOrb)

        lpCell: do iCell = 1, size(coords, dim=3)

          if (gridInterType == gridInterpolTypes%explicit) then
            basisTmp(:,:,:) = &
                & explicitSTOcalculation(regionGridDat(iGrid), subgridsDat(iSubgrid), &
                & (1.0_dp, 0.0_dp), coords(:, iAtom, iCell), &
                & cartCoords(:,:,:,tiling(1, iGrid):tiling(2, iGrid)), stos(iStos), iAng, iM, &
                & square=tAddDensities)
          else if (gridInterType .eq. gridInterpolTypes%trivial) then
            basisTmp(:,:,:) = &
                & gridInterpolationTrivial(regionGridDat(iGrid), subgridsDat(iSubgrid), &
                & (1.0_dp, 0.0_dp), coords(:, iAtom, iCell), square=tAddDensities)
          else if (gridInterType .eq. gridInterpolTypes%linear) then
            basisTmp(:,:,:) = &
                & gridInterpolationLinear(regionGridDat(iGrid), subgridsDat(iSubgrid), &
                & (1.0_dp, 0.0_dp), coords(:, iAtom, iCell), square=tAddDensities)
          end if

          do iKPoint = 1, size(requiredKPoints)
            basis(:,:,:, iKPoint) = basis(:,:,:, iKPoint) + basisTmp(:,:,:) * phases(iCell, iKPoint)
          end do

        end do lpCell

        levelInd = 1
        lpLevel: do iL = statesTiling(1), statesTiling(2)
          currentLevel = requiredLevels(iL)
          lpKPoint: do iKPoint = 1, size(requiredKPointsForLevel(levelInd, :))
            currentKPoint = requiredKPointsForLevel(levelInd, iKPoint)

            if (currentKPoint .ne. 0) then

              lpSpin: do iSpin = 1, size(requiredSpins)
                currentSpin = requiredSpins(iSpin)

                globalGridsDat(:,:, tiling(1, iGrid):tiling(2, iGrid), levelInd, iKPoint, iSpin) = &
                    & globalGridsDat(:,:, tiling(1, iGrid):tiling(2, iGrid), levelInd, &
                    & iKPoint, iSpin) + coeffs(iOrb, currentLevel, currentKPoint, currentSpin) * &
                    & basis(:,:,:, currentKPoint)

              end do lpSpin
            end if

          end do lpKPoint

          levelInd = levelInd + 1

        end do lpLevel

      end do lpOrb
      deallocate(basisTmp)
      deallocate(basis)
    end do lpGrid
    !$omp end parallel do

  end subroutine subgridsToCachedGlobalGridsCplx


  subroutine checkParallelBasis(basis1, basis2)

    !> basis vectors of grids to compare, shapes: [3,3]
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


  ! Divides a grid in multiple parallel grid slices.
  subroutine parallelGridSlices(totGridDat, regionGridDat, parallelRegionNum, tiling)

    !> Grid data instance to divide in slices
    type(TGridData), intent(inout) :: totGridDat

    !> Grid data instances representing the slices of the total grid, shape: [nParallelRegions]
    type(TGridData), allocatable, intent(inout) :: regionGridDat(:)

    !> Number of parallel regions of the total grid
    integer, intent(in) :: parallelRegionNum

    !> start and end indices of all parallel tiles, exp. shape: [2, nRegions]
    integer, intent(in) :: tiling(:,:)

    !> temporary grid instance
    type(TGrid) :: tmpGrid
    real(dp), pointer :: pTmpData(:,:,:)

    !> auxiliary variable
    integer :: iGrid

    ! allocate regional global grid tiles
    allocate(regionGridDat(parallelRegionNum))

    ! assign subgrids to parallel regions of the total grid
    do iGrid = 1, parallelRegionNum
      tmpGrid = totGridDat%grid
      tmpGrid%lowerBounds(3) = tiling(1, iGrid) - 1
      tmpGrid%upperBounds(3) = tiling(2, iGrid) - 1
      tmpGrid%nPoints(3) = tiling(2, iGrid) - tiling(1, iGrid) + 1
      pTmpData => totGridDat%data(:,:, tiling(1, iGrid):tiling(2, iGrid))
      call TGridData_init(regionGridDat(iGrid), tmpGrid, pTmpData,&
          & rwTabulationType=rwTabulationTypes%explicit)
    end do

  end subroutine parallelGridSlices


  ! Calculates the cartesian coordinates of the grid points.
  subroutine calcCartCoords(totGridDat, cartCoords)

    !> Grid data instance to calculate the cartesian coordinates for
    type(TGridData), intent(inout) :: totGridDat

    !> Cartesian coordinates of the grid points, shape: [3, nxPoints, nyPoints, nzPoints]
    real(dp), allocatable, intent(inout) :: cartCoords(:,:,:,:)

    !> Temporary stored cartesian coordinates of one grid point, shape: [3,1]
    real(dp), allocatable :: cartcoordsTmp(:,:)

    !> Auxiliary variables
    integer :: iGrid, kk, jj, ii, idx(3, 1), aa(3)

    allocate(cartCoords(3, totGridDat%grid%nPoints(1), &
        & totGridDat%grid%nPoints(2), totGridDat%grid%nPoints(3)))

    do kk = totGridDat%grid%lowerBounds(3), totGridDat%grid%upperBounds(3)
      idx(3, 1) = kk
      do jj = totGridDat%grid%lowerBounds(2), totGridDat%grid%upperBounds(2)
        idx(2, 1) = jj
        do ii = totGridDat%grid%lowerBounds(1), &
              & totGridDat%grid%upperBounds(1)
          idx(1, 1) = ii
          aa(:) = idx(:, 1) - totGridDat%grid%lowerBounds + 1
          call totGridDat%grid%gridcoordToCartesian(idx, cartCoordsTmp)
          cartCoords(:, aa(1), aa(2), aa(3)) = cartcoordsTmp(:,1)
        end do
      end do
    end do

  end subroutine calcCartCoords


#:for DTYPE, NAME in [('complex', 'Cplx'), ('real', 'Real')]

  !> Takes a 2d array of eigenvectors where the first dimension runs over the orbitals and the
  !> second dimension runs over the eigenvectors for each level, k-point and spin. It converts this
  !> array to a 4d array, where the dimensions are the following:
  !> [orbitals, eigenvectors, k-points, spin].
  subroutine modify${NAME}$Eigenvecs(eigenvecsOld, eigenvecsNew, kPointNum, spinNum)

    !> eigenvectors with old shape: [nOrb, nStates]
    ${DTYPE}$(dp), intent(in) :: eigenvecsOld(:,:)

    !> eigenvectors with new shape: [orbitals, eigenvectors, k-points, spin]
    ${DTYPE}$(dp), intent(out) :: eigenvecsNew(:,:,:,:)

    integer, intent(in) :: spinNum
    integer, intent(in) :: kPointNum

    integer :: nBlock

    !> auxiliary variables for do loops
    integer :: iKP, iS, iL

    nBlock = 0
    do iS = 1, spinNum
      do iKP = 1, kPointNum
        do iL = 1, size(eigenvecsOld, dim=1)
          eigenvecsNew(:, iL, iKP, iS) = eigenvecsOld(:, iL + size(eigenvecsOld, dim=1) * nBlock)
        end do
        nBlock = nBlock + 1
      end do
    end do

  end subroutine modify${NAME}$Eigenvecs

#:endfor


  pure function crossProduct(vec1, vec2) result(cross)

    !> vectors to calculate the cross product for
    real(dp), intent(in) :: vec1(:), vec2(:)

    !> resulting cross product
    real(dp) :: cross(3)

    cross(1) = vec1(2) * vec2(3) - vec1(3) * vec2(2)
    cross(2) = vec1(3) * vec2(1) - vec1(1) * vec2(3)
    cross(3) = vec1(1) * vec2(2) - vec1(2) * vec2(1)

  end function crossProduct


#:for DTYPE, NAME in [('complex', 'Cplx'), ('real', 'Real')]

  !> Calculates the values of the total grid with an explicit calculation of the sherical harmonics.
  function explicitSTOcalculation${NAME}$(regionGridDat, griddata, eigvec, position, &
      & cartCoords, stos, ll, mm, square) result(stoValue)

    !> representation of volumetric grid and data
    class(TGridData), intent(in) :: regionGridDat

    !> volumetric data to add to current instance
    class(TGridData), intent(in) :: griddata

    !> eigenvector of current basis
    ${DTYPE}$(dp), intent(in) :: eigvec

    !> position where to place subgrid onto grid, expected shape: [3]
    real(dp), intent(in) :: position(:)

    !> Cartesian coordinates of the total grid, shape: [3, nxPoints, nyPoints, nzPoints]
    real(dp), intent(in) :: cartCoords(:,:,:,:)

    !> tabulated radial wavefunction, shape: [nDistances, 3]
    type(TSlaterOrbital), intent(in) :: stos

    !> Angular momentum quantum number for current STO
    integer, intent(in) :: ll

    !> Magnetic quantum number for current STO
    integer, intent(in) :: mm

    !> True, if data should be squared, i.e. for densities
    logical, intent(in), optional :: square

    !> lower and upper bounds of intersection in integer grid coordinates, shape: [3,2]
    integer, allocatable :: luGc(:,:)

    !> lower and upper bounds of intersection in real grid coordinates, shape: [3,2]
    real(dp), allocatable :: luSubGc(:,:)

    !> equals square if present, otherwise false
    logical :: tSquare

    !> true if there is an intersection between both grids, false otherwise
    logical :: intersec

    !> auxiliary variables
    integer :: ii, jj, kk, inds(3), endIdx, loc, idx(3, 1)

    !> interpolated data, shape: [nxRegionPoints, nyRegionPoints, nzRegionPoints]
    real(dp), allocatable :: stoValue(:,:,:)

    !> Distance from origin of subgrid to coordinate of grid point to calculate
    real(dp) :: distance
    real(dp) :: rwfValue

    if (present(square)) then
      tSquare = square
    else
      tSquare = .false.
    end if

    allocate(stoValue(regionGridDat%grid%upperBounds(1) + 1, regionGridDat%grid%upperBounds(2) + 1,&
        & regionGridDat%grid%upperBounds(3) - regionGridDat%grid%lowerBounds(3) + 1))

    call checkParallelBasis(regionGridDat%grid%basis, griddata%grid%basis)

    ! Get lower and upper bounds of intersection in grid and subgrid coordinates.
    call regionGridDat%grid%getSubgridRanges(griddata%grid, position, luGc, luSubGc, intersec)

    stoValue(:,:,:) = 0.0_dp
    endIdx = size(stos%gridValue, dim=1)

    if (intersec) then
      if (tSquare) then
        do kk = 1, luGc(3, 2) - luGc(3, 1) + 1
          idx(3, 1) = kk
          do jj = 1, luGc(2, 2) - luGc(2, 1) + 1
            idx(2, 1) = jj
            do ii = 1, luGc(1, 2) - luGc(1, 1) + 1
              idx(1, 1) = ii
              inds(:) = idx(:, 1) + luGc(:, 1) - regionGridDat%grid%lowerBounds
              distance = norm2(cartCoords(:, inds(1), inds(2), inds(3)) - position)
              loc = max(min(locateByBisection(stos%gridValue(:, 1), distance), endIdx), 1)
              if (loc >= endIdx) then
                stoValue(inds(1), inds(2), inds(3)) = 0.0_dp
              else
                stoValue(inds(1), inds(2), inds(3)) = (linearInterpolation(stos%gridValue(loc, 1),&
                    & stos%gridValue(loc + 1, 1), stos%gridValue(loc, 2), &
                    & stos%gridValue(loc + 1, 2), distance) * &
                    & RealTessY(ll, mm, cartCoords(:, inds(1), inds(2), inds(3)) - position))**2 &
                    & * eigvec
              end if
            end do
          end do
        end do
      else
        do kk = 1, luGc(3, 2) - luGc(3, 1) + 1
          idx(3, 1) = kk
          do jj = 1, luGc(2, 2) - luGc(2, 1) + 1
            idx(2, 1) = jj
            do ii = 1, luGc(1, 2) - luGc(1, 1) + 1
              idx(1, 1) = ii
              inds(:) = idx(:, 1) + luGc(:, 1) - regionGridDat%grid%lowerBounds
              distance = norm2(cartCoords(:, inds(1), inds(2), inds(3)) - position)
              loc = max(min(locateByBisection(stos%gridValue(:, 1),&
                  & distance), endIdx), 1)
              if (loc >= endIdx) then
                stoValue(inds(1), inds(2), inds(3)) = 0.0_dp
              else
                stoValue(inds(1), inds(2), inds(3)) = linearInterpolation(stos%gridValue(loc, 1),&
                    & stos%gridValue(loc + 1, 1), stos%gridValue(loc, 2), &
                    & stos%gridValue(loc + 1, 2), distance) * &
                    & RealTessY(ll, mm, cartCoords(:, inds(1), inds(2), inds(3)) - position) &
                    & * eigvec
              end if
            end do
          end do
        end do
      end if
    end if

  end function explicitSTOcalculation${NAME}$

#:endfor


#:for DTYPE, NAME in [('complex', 'Cplx'), ('real', 'Real')]

  !> Executes the three dimensional linear grid interpolation with real or complex values.
  function gridInterpolation${NAME}$linear(regionGridDat, griddata, eigvec, position, square) &
      & result(interp)

    !> representation of volumetric grid and data
    class(TGridData), intent(in) :: regionGridDat

    !> volumetric data to add to current instance
    class(TGridData), intent(in) :: griddata

    !> eigenvector of current basis
    ${DTYPE}$(dp), intent(in) :: eigvec

    !> position where to place subgrid onto grid, expected shape: [3]
    real(dp), intent(in) :: position(:)

    !> True, if data should be squared, i.e. for densities
    logical, intent(in), optional :: square

    !> lower and upper bounds of intersection in integer grid coordinates, shape: [3,2]
    integer, allocatable :: luGc(:,:)

    !> lower and upper bounds of intersection in real grid coordinates, shape: [3,2]
    real(dp), allocatable :: luSubGc(:,:)

    !> real grid coordinates of subgrid that intersect with grid,
    !> shape: [3, nxIntersecPoints, nyIntersecPoints, nzIntersecPoints]
    real(dp), allocatable :: intersecSubGc(:,:,:,:)

    !> floor rounded lower subcube corners of grid points,
    !> shape: [3, nxIntersecPoints, nyIntersecPoints, nzIntersecPoints]
    integer, allocatable :: roundedIntersecSubGc(:,:,:,:)

    !> data corresponding to the eight corners of the subcube,
    !> expected shape: [2, 2, 2, nxIntersecPoints, nyIntersecPoints, nzIntersecPoints]
    real(dp), allocatable :: scData(:,:,:,:,:,:)

    !> equals square if present, otherwise false
    logical :: tSquare

    !> true if there is an intersection between both grids, false otherwise
    logical :: intersec

    !> data indices, corresponding to a grid point
    integer :: inds(3)

    !> auxiliary variables
    integer :: ii, jj, kk, idxval(3)

    !> interpolated data, shape: [nxRegionPoints, nyRegionPoints, nzRegionPoints]
    ${DTYPE}$(dp), allocatable :: interp(:,:,:)

    if (present(square)) then
      tSquare = square
    else
      tSquare = .false.
    end if

    allocate(interp(regionGridDat%grid%upperBounds(1) + 1, regionGridDat%grid%upperBounds(2) + 1, &
        & regionGridDat%grid%upperBounds(3) - regionGridDat%grid%lowerBounds(3) + 1))

    call checkParallelBasis(regionGridDat%grid%basis, griddata%grid%basis)

    ! Get lower and upper bounds of intersection in grid and subgrid coordinates.
    call regionGridDat%grid%getSubgridRanges(griddata%grid, position, luGc, luSubGc, intersec)

    ! Calculates for any grid point which intersects with the subgrid the subgrids coordinates of
    ! the eight subgrid points which surround the grid point.
    call regionGridDat%grid%getIntersecGridPoints(griddata%grid, luGc, luSubGc, intersecSubGc)
    roundedIntersecSubGc = floor(intersecSubGc)

    ! Calculates the values of the subgrid cubes around total grid points in intersection area.
    call griddata%getSubcubeCornsAndVals(roundedIntersecSubGc, scData)

    if (tSquare) then
      scData(:,:,:,:,:,:) = scData * scData
    end if

    interp(:,:,:) = 0.0_dp

    if (intersec) then
      do kk = 1, luGc(3, 2) - luGc(3, 1) + 1
        idxval(3) = kk
        do jj = 1, luGc(2, 2) - luGc(2, 1) + 1
          idxval(2) = jj
          do ii = 1, luGc(1, 2) - luGc(1, 1) + 1
            idxval(1) = ii
            inds(:) = idxval + luGc(:, 1) - regionGridDat%grid%lowerBounds
            interp(inds(1), inds(2), inds(3)) = eigvec * &
                & trilinearInterpolation(real(roundedIntersecSubGc(:, ii, jj, kk), dp),&
                & real(roundedIntersecSubGc(:, ii, jj, kk) + 1, dp), scData(:,:,:, ii, jj, kk),&
                & intersecSubGc(:, ii, jj, kk))
          end do
        end do
      end do
    end if

  end function gridInterpolation${NAME}$linear

#:endfor


#:for DTYPE, NAME in [('complex', 'Cplx'), ('real', 'Real')]

  !> Executes the three dimensional trivial grid interpolation with real or complex values.
  function gridInterpolation${NAME}$trivial(regionGridDat, griddata, eigvec, position, square) &
      & result(grid_interp)

    !> representation of volumetric grid and data
    class(TGridData), intent(in) :: regionGridDat

    !> volumetric data to add to current instance
    class(TGridData), intent(in) :: griddata

    !> eigenvector of current basis
    ${DTYPE}$(dp), intent(in) :: eigvec

    !> position where to place subgrid onto grid, expected shape: [3]
    real(dp), intent(in) :: position(:)

    !> True, if data should be squared, i.e. for densities
    logical, intent(in), optional :: square

    !> real grid coordinates of subgrid that intersect with grid,
    !> shape: [3, nxIntersecPoints, nyIntersecPoints, nzIntersecPoints]
    real(dp), allocatable :: intersecSubGc(:,:,:,:)

    !> lower and upper bounds of intersection in integer grid coordinates, shape: [3,2]
    integer, allocatable :: luGc(:,:)

    !> lower and upper bounds of intersection in real grid coordinates, shape: [3,2]
    real(dp), allocatable :: luSubGc(:,:)

    !> interpolated data, shape: [nxRegionPoints, nyRegionPoints, nzRegionPoints]
    real(dp), allocatable :: interp(:,:,:)

    !> equals square if present, otherwise false
    logical :: tSquare

    !> true if there is an intersection between both grids, false otherwise
    logical :: intersec

    !> data indices, corresponding to a grid point
    integer :: inds(3)

    !> auxiliary variables
    integer :: ii, jj, kk, idx(3)

    ${DTYPE}$(dp), allocatable :: grid_interp(:,:,:)

    allocate(grid_interp(regionGridDat%grid%upperBounds(1) + 1, &
        & regionGridDat%grid%upperBounds(2) + 1, &
        & regionGridDat%grid%upperBounds(3) - regionGridDat%grid%lowerBounds(3) + 1))

    grid_interp(:,:,:) = 0.0_dp

    if (present(square)) then
      tSquare = square
    else
      tSquare = .false.
    end if

    call checkParallelBasis(regionGridDat%grid%basis, griddata%grid%basis)

    ! Get lower and upper bounds of intersection in grid and subgrid coordinates.
    call regionGridDat%grid%getSubgridRanges(griddata%grid, position, luGc, luSubGc, intersec)

    ! Calculates for any grid point which intersects with the subgrid the subgrids coordinates of
    ! the eight subgrid points which surround the grid point.
    call regionGridDat%grid%getIntersecGridPoints(griddata%grid, luGc, luSubGc, intersecSubGc)
    call griddata%getInterpolatedValueGc(intersecSubGc, interp)

    if (tSquare) then
      interp = interp**2
    end if

    if (intersec) then
      do kk = 1, luGc(3, 2) - luGc(3, 1) + 1
        idx(3) = kk
        do jj = 1, luGc(2, 2) - luGc(2, 1) + 1
          idx(2) = jj
          do ii = 1, luGc(1, 2) - luGc(1, 1) + 1
            idx(1) = ii
            inds(:) = idx + luGc(:, 1) - regionGridDat%grid%lowerBounds
            grid_interp(inds(1), inds(2), inds(3)) = grid_interp(inds(1), inds(2), inds(3))&
                & + eigvec * cmplx(interp(ii, jj, kk), 0.0_dp, dp)
          end do
        end do
      end do
    end if

  end function gridInterpolation${NAME}$trivial

#:endfor

end module waveplot_grids
