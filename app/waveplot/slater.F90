!--------------------------------------------------------------------------------------------------!
!  DFTB+: general package for performing fast atomistic simulations                                !
!  Copyright (C) 2006 - 2022  DFTB+ developers group                                               !
!                                                                                                  !
!  See the LICENSE file for terms of usage and distribution.                                       !
!--------------------------------------------------------------------------------------------------!

#:include 'common.fypp'

!> Routines to calculate a Slater type orbital (STO)
module waveplot_slater

  use dftbp_common_accuracy, only : dp
  implicit none

  private
  public :: RealTessY, SlaterOrbital_getValue_explicit, TSlaterOrbital,&
      & initSto, SlaterOrbital_getValue
  save


  !> Data for STOs
  type TSlaterOrbital

    !> Number of coeffs in aa
    integer :: nPow
    
    !> Number of exponentials
    integer :: nAlpha

    !> Angular momentum
    integer :: ll

    !> Summation coefficients
    real(dp), allocatable :: aa(:,:)

    !> Exponential coefficients
    real(dp), allocatable :: alpha(:)

    !> Values on the radial grid for the rwf, shape: [nGrid, 2] or [nGrid, 3]. Last dimension
    !> contains the distance and the value of the rwf at that distance if the rwf is explicitly
    !> calculated from the coefficients. If the rwf is read in from an wfc.hsd file, then the
    !> the last dimension contains also the second derivative of the rwf.
    real(dp), allocatable :: gridValue(:,:)

    !> Distance between two grid points
    real(dp) :: gridDist

    !> Number of grid points
    integer :: nGrid
  end type TSlaterOrbital


  !> Initialises a SlaterOrbital
  interface initSto
    module procedure SlaterOrbital_init
  end interface


contains


  !> Returns the real tesseral spherical harmonics in a given point
  !> This function only works for angular momenta between 0 and 3 (s-f).
  function RealTessY(ll, mm, coord, rrOpt) result (rty)

    !> Angular momentum of the spherical harmonics (0 <= ll <= 3)
    integer, intent(in) :: ll

    !> Magnetic quantum number
    integer, intent(in) :: mm

    !> Coordinate where the value should be calculated
    real(dp), intent(in) :: coord(:)

    !> Length of the coordinate vector, if known in advance
    real(dp), intent(in), optional :: rrOpt

    real(dp) :: rty

    real(dp) :: rr, xx, yy, zz

    if (present(rrOpt)) then
      rr = rrOpt
    else
      rr = sqrt(sum(coord**2))
    end if

    @:ASSERT(ll >= 0 .and. ll <= 3)
    @:ASSERT(abs(mm) <= ll)
    @:ASSERT(size(coord) == 3)
    @:ASSERT(rr >= 0.0_dp)

    xx = coord(1)
    yy = coord(2)
    zz = coord(3)

    if (rr < epsilon(1.0_dp) .and. ll /= 0) then
      rty = 0.0_dp
      return
    end if

    select case (ll)
    case(0)
      rty = 0.2820947917738782_dp
    case(1)
      select case(mm)
      case(-1)
        ! y
        rty = 0.4886025119029198_dp * yy / rr
      case(0)
        ! z
        rty = 0.4886025119029198_dp * zz / rr
      case(1)
        ! x
        rty = 0.4886025119029198_dp * xx / rr
      end select
    case(2)
      select case(mm)
      case(-2)
        ! xy
        rty = 1.092548430592079_dp * xx * yy / rr**2
      case(-1)
        ! yz
        rty = 1.092548430592079_dp * yy * zz / rr**2
      case(0)
        ! z**2
        rty = -0.3153915652525200_dp * (-2.0_dp * zz**2 + xx**2 + yy**2) / rr**2
      case(1)
        ! xz
        rty = 1.092548430592079_dp * xx * zz / rr**2
      case(2)
        ! x**2-y**2
        rty = 0.5462742152960395_dp * (xx**2 - yy**2) / rr**2
      end select
    case(3)

      !> general set for f orbitals (not cubic), see
      !> http://winter.group.shef.ac.uk/orbitron/AOs/4f/equations.html
      select case (mm)
      case(-3)
        ! y(3x**2-y**2)
        rty = 0.5900435899266435_dp * yy * (3.0_dp * xx**2 - yy**2) &
            &/ rr**3
      case(-2)
        ! x**2+y**2+z**2
        rty = 2.890611442640554_dp * xx * yy *zz / rr**3
      case(-1)
        ! yz**2
        rty = -0.4570457994644658_dp * (-4.0_dp * zz**2 + xx**2 + yy**2) * yy &
            &/ rr**3
      case(0)
        ! z**3
        rty = -0.3731763325901155_dp * zz &
            &*(-2.0_dp * zz**2 + 3.0_dp * xx**2 + 3.0_dp * yy**2)/ rr**3
      case(1)
        ! xz**2
        rty = -0.4570457994644658_dp * (-4.0_dp * zz**2 + xx**2 + yy**2) * xx &
            &/ rr**3
      case(2)
        ! z(x**2-y**2)
        rty = 1.445305721320277_dp * zz * (xx**2 - yy**2) / rr**3
      case(3)
        ! x(x**2-3y**2)
        rty = 0.5900435899266435_dp * xx * (xx**2 - 3.0_dp * yy**2) / rr**3
      end select
    end select

  end function RealTessY


  !> Initialises a SlaterOrbital.
  subroutine SlaterOrbital_init(this, aa, alpha, ll, cutoff, resolution, rwf)

    !> SlaterOrbital instance to initialise
    type(TSlaterOrbital), intent(inout) :: this

    !> Summation coefficients (nCoeffPerAlpha, nAlpha)
    real(dp), intent(in) :: aa(:,:)

    !> Exponential coefficients
    real(dp), intent(in) :: alpha(:)

    !> Angular momentum of the orbital
    integer, intent(in) :: ll

    !> Cutoff, after which orbital is assumed to be zero
    real(dp), intent(in) :: cutoff

    !> Grid distance for the orbital
    real(dp), intent(in), optional :: resolution

    !>  Tabulated values of the radial wf: [distance, value, 2nd derivative]
    real(dp), intent(in), optional :: rwf(:,:)

    integer :: nAlpha
    integer :: nPow
    integer :: iGrid
    real(dp) :: rr

    nAlpha = size(alpha)
    nPow = size(aa, dim=1)

    @:ASSERT(cutoff > 0.0_dp)
    if (present(resolution)) then
      @:ASSERT(resolution > 0.0_dp)
    end if

    allocate(this%aa(nPow, nAlpha))
    allocate(this%alpha(nAlpha))

    this%aa(:,:) = aa
    this%alpha(:) = -1.0_dp * alpha
    this%nPow = nPow
    this%nAlpha = nAlpha
    this%ll = ll

    if (present(rwf)) then
      allocate(this%gridValue(size(rwf, dim=1), 3))
      this%gridValue = rwf
      this%nGrid = size(rwf, dim=1)
    else
      ! Obtain STO on a grid
      this%nGrid = floor(cutoff / resolution) + 2
      this%gridDist = resolution
      allocate(this%gridValue(this%nGrid, 2))
      do iGrid = 1, this%nGrid
        rr = real(iGrid - 1, dp) * resolution
        this%gridValue(iGrid, 1) = rr
        call SlaterOrbital_getValue_explicit(ll, nPow, nAlpha, aa, this%alpha, rr,&
            & this%gridValue(iGrid, 2))
      end do
    end if

  end subroutine SlaterOrbital_init


  !> Retunrns the value of the SlaterOrbital in a given point
  subroutine SlaterOrbital_getValue(this, rr, sto)

    !> SlaterOrbital instance
    type(TSlaterOrbital), intent(in) :: this

    !> Distance, where STO should be calculated
    real(dp), intent(in) :: rr

    !> Contains the value of the function on return
    real(dp), intent(out) :: sto

    integer :: ind
    real(dp) :: frac

    @:ASSERT(rr >= 0.0_dp)

    ! ind = 1 means zero distance as rr = (ind - 1) * gridDist
    ind = floor(rr / this%gridDist) + 1
    if (ind < this%nGrid) then
      frac = mod(rr, this%gridDist) / this%gridDist
      sto = (1.0_dp - frac) * this%gridValue(ind, 2) + frac * this%gridValue(ind+1, 2)
    else
      sto = 0.0_dp
    end if

  end subroutine SlaterOrbital_getValue


  !> Calculates the value of an STO analytically
  subroutine SlaterOrbital_getValue_explicit(ll, nPow, nAlpha, aa, alpha, rr, sto)

    !> Angular momentum of the STO
    integer, intent(in) :: ll

    !> Maximal power of the distance in the STO
    integer, intent(in) :: nPow

    !> Number of exponential coefficients
    integer, intent(in) :: nAlpha

    !> summation coefficients, shape: [nCoeffPerAlpha, nAlpha]
    real(dp), intent(in) :: aa(:,:)

    !> exponential coefficients
    real(dp), intent(in) :: alpha(:)

    !> Distance, where the STO should be calculated
    real(dp), intent(in) :: rr

    !> Value of the STO on return
    real(dp), intent(out) :: sto

    real(dp) :: pows(nPow)
    real(dp) :: rTmp
    integer :: ii, jj

    ! Avoid 0.0**0 as it may lead to arithmetic exception
    if (ll == 0 .and. rr < epsilon(1.0_dp)) then
      rTmp = 1.0_dp
    else
      rTmp = rr**ll
    end if
    do ii = 1, nPow
      pows(ii) = rTmp
      rTmp = rTmp * rr
    end do
    sto = 0.0_dp
    do ii = 1, nAlpha
      rTmp = 0.0_dp
      do jj = 1, nPow
        rTmp = rTmp + aa(jj, ii) * pows(jj)
      end do
      sto = sto + rTmp * exp(alpha(ii) * rr)
    end do

  end subroutine SlaterOrbital_getValue_explicit

end module waveplot_slater
