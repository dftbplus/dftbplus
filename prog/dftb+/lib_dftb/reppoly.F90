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
  subroutine RepPoly_init(self, inp)

    !> Polynomial repulsive
    type(TRepPoly), intent(out) :: self

    !> Input parameters for the polynomial repulsive
    type(TRepPolyIn), intent(in) :: inp

    @:ASSERT(.not. self%tInit)
    @:ASSERT(inp%cutoff >= 0.0_dp)

    self%polyCoeffs(:) = inp%polyCoeffs(:)
    self%cutoff = inp%cutoff
    self%tInit = .true.

  end subroutine RepPoly_init


  !> Returns cutoff of the repulsive
  function RepPoly_getCutoff(self) result(cutoff)

    !> self Polynomial repulsive
    type(TRepPoly), intent(in) :: self

    !> Cutoff
    real(dp) :: cutoff

    @:ASSERT(self%tInit)
    cutoff = self%cutoff

  end function RepPoly_getCutoff


  !> Returns energy of the repulsive for a given distance
  subroutine RepPoly_getEnergy(self, res, rr)

    !> Polynomial repulsive
    type(TRepPoly), intent(in) :: self

    !> Energy contribution
    real(dp), intent(out) :: res

    !> Distance between interacting atoms
    real(dp), intent(in) :: rr

    real(dp) :: rrr
    integer :: ii

    @:ASSERT(self%tInit)

    if (rr >= self%cutoff) then
      res = 0.0_dp
    else
      rrr = self%cutoff - rr
      res = self%polyCoeffs(powMax)
      do ii = powMax-1, powMin, -1
        res = res * rrr + self%polyCoeffs(ii)
      end do
      do ii = powMin, 1, -1
        res = res * rrr
      end do
    end if

  end subroutine RepPoly_getEnergy


  !> Returns gradient of the repulsive for a given distance.
  subroutine RepPoly_getEnergyDeriv(self, xx, res, res2)

    !> Polynomial repulsive.
    type(TRepPoly), intent(in) :: self

    !> Actual vector between atoms
    real(dp), intent(in) :: xx(3)

    !> Resulting contribution
    real(dp), intent(out) :: res(3)

    !> second derivatives
    real(dp), intent(out), optional :: res2(3,3)

    integer :: ii, jj
    real(dp) :: r2, rr, rrr, d1, d2

  @:ASSERT(self%tInit)
    r2 = sum(xx**2)
    rr = sqrt(r2)

    res(:) = 0.0_dp
    if (present(res2)) then
      res2 = 0.0_dp
    end if

    if (rr < self%cutoff) then
      rrr = self%cutoff - rr

      !do ii = powMin1, powMax
      !  d1 = d1 + real(ii, dp) * self%polyCoeffs(ii) * rrr**(ii-1)
      !end do
      d1 = real(powMax, dp) * self%polyCoeffs(powMax)
      do ii = powMax-1, powMin1, -1
        d1 = d1 * rrr + real(ii, dp) * self%polyCoeffs(ii)
      end do
      do ii = powMin1, 2, -1
        d1 = d1 * rrr
      end do

      res(:) = -d1 * xx(:) / rr

      if (present(res2)) then
        do ii = powMin1, powMax
          d2 = d2 + real(ii*(ii-1), dp) * self%polyCoeffs(ii) * rrr**(ii-2)
        end do
        do ii = 1, 3
          do jj = 1, 3
            res2(ii,jj) = xx(ii)*xx(jj)
          end do
        end do
        res2 = res2 * ( d2 / r2 + d1 / (rr*r2) )
        do ii = 1, 3
          res2(ii,ii) = res2(ii,ii) - d1 / rr
        end do
      end if

    end if

  end subroutine RepPoly_getEnergyDeriv

end module dftbp_reppoly
