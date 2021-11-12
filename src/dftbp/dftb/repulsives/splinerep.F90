!--------------------------------------------------------------------------------------------------!
!  DFTB+: general package for performing fast atomistic simulations                                !
!  Copyright (C) 2006 - 2021  DFTB+ developers group                                               !
!                                                                                                  !
!  See the LICENSE file for terms of usage and distribution.                                       !
!--------------------------------------------------------------------------------------------------!

#:include 'common.fypp'

!> Implements a repulsive potential between two atoms represented by cubic splines.
module dftbp_dftb_repulsives_splinerep
  use dftbp_common_accuracy, only : dp, minNeighDist
  use dftbp_dftb_repulsives_pairrepulsive, only : TPairRepulsive
  use dftbp_math_bisect, only : bisection
  implicit none

  private
  public :: TSplineRepInp, TSplineRep, TSplineRep_init


  !> Initialisation type for ORepSpline
  type TSplineRepInp

    !> Starting pos. for each spline
    real(dp), allocatable :: xStart(:)

    !> Spline coeffs (3, nSpline-1)
    real(dp), allocatable :: spCoeffs(:,:)

    !> Coeffs. of the last spline
    real(dp) :: spLastCoeffs(6)

    !> Coeffs for the exponential head
    real(dp) :: expCoeffs(3)

    !> Cutoff for the last spline
    real(dp) :: cutoff

  end type TSplineRepInp


  !> Contains the spline representation of a repulsive.
  type, extends(TPairRepulsive) :: TSplineRep
    private

    !> Nr. of splines.
    integer :: nSpline

    !> Starting point for each spline
    real(dp), allocatable :: xStart(:)

    !> Spline coeffs (3, nSpline-1)
    real(dp), allocatable :: spCoeffs(:,:)

    !> Coeffs of the last spline
    real(dp) :: spLastCoeffs(6)

    !> Exponential head
    real(dp) :: expCoeffs(3)

    !> Cutoff of the last spline
    real(dp) :: cutoff

  contains

    procedure :: getCutoff => TSplineRep_getCutoff
    procedure :: getValue => TSplineRep_getValue

  end type TSplineRep


  !> Structure constructor for TSplineRep
  interface TSplineRep
    module procedure TSplineRep_construct
  end interface TSplineRep


contains


  !> Initialises spline repulsive.
  subroutine TSplineRep_init(this, inp)

    !> Spline repulsive.
    type(TSplineRep), intent(out) :: this

    !> Input parameters for the spline repulsive.
    type(TSplineRepInp), intent(in) :: inp

    @:ASSERT(size(inp%xStart) > 0)
    @:ASSERT(size(inp%spCoeffs, dim=1) == 4)
    @:ASSERT(size(inp%spCoeffs, dim=2) == size(inp%xStart) - 1)
    @:ASSERT(inp%cutoff >= 0.0_dp)

    this%nSpline = size(inp%xStart)
    allocate(this%xStart(this%nSpline))
    allocate(this%spCoeffs(4, this%nSpline - 1))
    this%xStart(:) = inp%xStart(:)
    this%spCoeffs(:,:) = inp%spCoeffs(:,:)
    this%spLastCoeffs(:) = inp%spLastCoeffs(:)
    this%expCoeffs(:) = inp%expCoeffs
    this%cutoff = inp%cutoff

  end subroutine TSplineRep_init


  !> Constructor for the spline repulsive
  function TSplineRep_construct(inp) result(this)

    !> Input parameters for the spline repulsive.
    type(TSplineRepInp), intent(in) :: inp
    
    !> Spline repulsive.
    type(TSplineRep) :: this

    call TSplineRep_init(this, inp)

  end function TSplineRep_construct


  !> Returns cutoff of the repulsive.
  function TSplineRep_getCutoff(this) result(cutoff)

    !> Spline repulsive.
    class(TSplineRep), intent(in) :: this

    !> Cutoff.
    real(dp) :: cutoff

    cutoff = this%cutoff

  end function TSplineRep_getCutoff


  !> Returns energy of the repulsive for a given distance.
  subroutine TSplineRep_getValue(this, rr, energy, dEnergy, d2Energy)

    !> Spline repulsive.
    class(TSplineRep), intent(in) :: this

    !> Distance between interacting atoms.
    real(dp), intent(in) :: rr

    !> Repulsive energy
    real(dp), optional, intent(out) :: energy

    !> First derivative of the energy
    real(dp), optional, intent(out) :: dEnergy

    !> Second derivative of the energy
    real(dp), optional, intent(out) :: d2Energy

    integer :: imatch
    real(dp) :: dr

    @:ASSERT(rr >= 0.0_dp)

    if (rr < minNeighDist .or. rr >= this%cutoff) then
      if (present(energy)) energy = 0.0_dp
      if (present(dEnergy)) dEnergy = 0.0_dp
      if (present(d2Energy)) d2Energy = 0.0_dp
      return
    end if

    if (rr < this%xStart(1)) then
      call getExponentialHead(this%expCoeffs, rr, energy, dEnergy, d2Energy)
      return
    end if

    call bisection(iMatch, this%xStart, rr)
    dr = rr - this%xStart(iMatch)
    if (iMatch < this%nSpline) then
      call getSpline(this%spCoeffs(:, iMatch), dr, energy, dEnergy, d2Energy)
    else
      call getPolynomialTail(this%spLastCoeffs, dr, energy, dEnergy, d2Energy)
    end if

  end subroutine TSplineRep_getValue


  !> Calculates the exponential head of the repulsive
  subroutine getExponentialHead(coeffs, rr, energy, dEnergy, d2Energy)
    real(dp), intent(in) :: coeffs(:), rr
    real(dp), optional, intent(out) :: energy, dEnergy, d2Energy

    if (present(energy)) energy = exp(-coeffs(1) * rr + coeffs(2)) + coeffs(3)
    if (present(dEnergy)) dEnergy = -coeffs(1) * exp(-coeffs(1) * rr + coeffs(2))
    if (present(d2Energy)) d2Energy = coeffs(1)**2 * exp(-coeffs(1) * rr + coeffs(2))

  end subroutine getExponentialHead


  !> Calculates the 5ht order polynomial tail
  subroutine getPolynomialTail(coeffs, dr, energy, dEnergy, d2Energy)
    real(dp), intent(in) :: coeffs(:), dr
    real(dp), optional, intent(out) :: energy, dEnergy, d2Energy

    real(dp) :: xh
    integer :: ii

    if (present(energy)) then
      energy = coeffs(1)
      xh = dr
      do ii = 2, 6
        energy = energy + coeffs(ii) * xh
        xh = xh * dr
      end do
    end if

    if (present(dEnergy)) then
      dEnergy = 0.0_dp
      xh = 1.0_dp
      do ii = 2, 6
        dEnergy = dEnergy + real(ii - 1, dp) * coeffs(ii) * xh
        xh = xh * dr
      end do
    end if

    if (present(d2Energy)) then
      d2Energy = 0.0_dp
      xh = 1.0_dp
      do ii = 3, 6
        d2Energy = d2Energy + real(ii - 2, dp) * real(ii - 1, dp) * xh * coeffs(ii)
        xh = xh * dr
      end do
    end if

  end subroutine getPolynomialTail


  !> Calculates the cubic spline
  subroutine getSpline(coeffs, dr, energy, dEnergy, d2Energy)
    real(dp), intent(in) :: coeffs(:), dr
    real(dp), optional, intent(out) :: energy, dEnergy, d2Energy

    real(dp) :: xh
    integer :: ii

    if (present(energy)) then
      xh = dr
      energy = coeffs(1)
      do ii = 2, 4
        energy = energy + coeffs(ii) * xh
        xh = xh * dr
      end do
    end if

    if (present(dEnergy)) then
      xh = 1.0_dp
      dEnergy = 0.0_dp
      do ii = 2, 4
        dEnergy = dEnergy + real(ii - 1, dp) * coeffs(ii) * xh
        xh = xh * dr
      end do
    end if

    if (present(d2Energy)) then
      xh = 1.0_dp
      d2Energy = 0.0_dp
      do ii = 3, 4
        d2Energy = d2Energy + real(ii - 2, dp) * real(ii - 1, dp) * coeffs(ii) * xh
        xh = xh * dr
      end do
    end if

  end subroutine getSpline

end module dftbp_dftb_repulsives_splinerep
