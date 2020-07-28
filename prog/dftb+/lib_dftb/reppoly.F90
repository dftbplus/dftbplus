!--------------------------------------------------------------------------------------------------!
!  DFTB+: general package for performing fast atomistic simulations                                !
!  Copyright (C) 2006 - 2020  DFTB+ developers group                                               !
!                                                                                                  !
!  See the LICENSE file for terms of usage and distribution.                                       !
!--------------------------------------------------------------------------------------------------!

#:include 'common.fypp'

!> Implements a repulsive potential between two atoms represented by a polynomial of 9th degree
module dftbp_reppoly
  use dftbp_assert
  use dftbp_accuracy
  use dftbp_bisect
  implicit none
  private

  public :: powMin, powMax
  public :: TRepPolyIn, TRepPoly, init
  public :: getCutoff, getEnergy, getEnergyDeriv


  !> Minimal power appearing in the polynomial
  integer, parameter :: powMin = 2

  !> Maximal power appearing in the polynomial
  integer, parameter :: powMax = 9

  !> smallest power between 1 and powMin
  integer, parameter :: powMin1 = max(powMin, 1)


  !> Initialisation type for TRepPoly
  type TRepPolyIn

    !> Polynomial coefficients
    real(dp) :: polyCoeffs(powMin:powMax)

    !> Cutoff distance
    real(dp) :: cutoff
  end type TRepPolyIn


  !> Contains the polynomial representation of a repulsive.
  type TRepPoly
    private

    !> coefficients of the polynomial
    real(dp) :: polyCoeffs(powMin:powMax)

    !> Cutoff distance
    real(dp) :: cutoff

    !> initialised the repulsive
    logical :: tInit = .false.
  end type TRepPoly


  !> Initialises polynomial repulsive.
  interface init
    module procedure RepPoly_init
  end interface


  !> Returns cutoff of the repulsive.
  interface getCutoff
    module procedure RepPoly_getCutoff
  end interface


  !> Returns energy of the repulsive for a given distance.
  interface getEnergy
    module procedure RepPoly_getEnergy
  end interface


  !> Returns gradient of the repulsive for a given distance.
  interface getEnergyDeriv
    module procedure RepPoly_getEnergyDeriv
  end interface

contains


  !> Initialises polynomial repulsive
  subroutine RepPoly_init(this, inp)

    !> Polynomial repulsive
    type(TRepPoly), intent(out) :: this

    !> Input parameters for the polynomial repulsive
    type(TRepPolyIn), intent(in) :: inp

    @:ASSERT(.not. this%tInit)
    @:ASSERT(inp%cutoff >= 0.0_dp)

    this%polyCoeffs(:) = inp%polyCoeffs(:)
    this%cutoff = inp%cutoff
    this%tInit = .true.

  end subroutine RepPoly_init


  !> Returns cutoff of the repulsive
  function RepPoly_getCutoff(this) result(cutoff)

    !> self Polynomial repulsive
    type(TRepPoly), intent(in) :: this

    !> Cutoff
    real(dp) :: cutoff

    @:ASSERT(this%tInit)
    cutoff = this%cutoff

  end function RepPoly_getCutoff


  !> Returns energy of the repulsive for a given distance
  subroutine RepPoly_getEnergy(this, res, rr)

    !> Polynomial repulsive
    type(TRepPoly), intent(in) :: this

    !> Energy contribution
    real(dp), intent(out) :: res

    !> Distance between interacting atoms
    real(dp), intent(in) :: rr

    real(dp) :: rrr
    integer :: ii

    @:ASSERT(this%tInit)

    if (rr >= this%cutoff) then
      res = 0.0_dp
    else
      rrr = this%cutoff - rr
      res = this%polyCoeffs(powMax)
      do ii = powMax-1, powMin, -1
        res = res * rrr + this%polyCoeffs(ii)
      end do
      do ii = powMin, 1, -1
        res = res * rrr
      end do
    end if

  end subroutine RepPoly_getEnergy


  !> Returns gradient of the repulsive for a given distance.
  subroutine RepPoly_getEnergyDeriv(this, res, xx, d2)

    !> Polynomial repulsive.
    type(TRepPoly), intent(in) :: this

    !> Resulting contribution
    real(dp), intent(out) :: res(3)

    !> Actual vector between atoms
    real(dp), intent(in) :: xx(3)

    !> second derivative in direction of repulsive derivative
    real(dp), intent(out), optional :: d2

    integer :: ii
    real(dp) :: rr, rrr, xh

    @:ASSERT(this%tInit)

    res(:) = 0.0_dp
    if (present(d2)) then
      d2 = 0.0_dp
    end if

    rr = sqrt(sum(xx**2))
    if (rr >= this%cutoff) then
      return
    end if

    rrr = this%cutoff - rr
    ! ((n * c_n * x + (n-1) * c_{n-1}) * x + (n-2) * c_{n-2}) * x + ... c_2
    xh = real(powMax, dp) * this%polyCoeffs(powMax)
    do ii = powMax - 1, powMin1, -1
      xh = xh * rrr + real(ii, dp) * this%polyCoeffs(ii)
    end do
    do ii = powMin1, 2, -1
      xh = xh * rrr
    end do
    res(:) = -xh * xx(:) / rr

    if (present(d2)) then
      xh = real(powMax * (powMax -1), dp) * this%polyCoeffs(powMax)
      do ii = powMax - 1, powMin1, -1
        xh = xh * rrr + real(ii * (ii - 1), dp) * this%polyCoeffs(ii)
      end do
      ! reserved for future expansion for repulsive potentials using higher
      ! order terms only:
      ! Commented as it gives a warning
      !do ii = powMin1, 3, -1
      !  xh = xh * rrr
      !end do
      d2 = xh
    end if

  end subroutine RepPoly_getEnergyDeriv

end module dftbp_reppoly
