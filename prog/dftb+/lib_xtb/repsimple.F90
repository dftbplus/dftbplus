!--------------------------------------------------------------------------------------------------!
!  DFTB+: general package for performing fast atomistic simulations                                !
!  Copyright (C) 2006 - 2020  DFTB+ developers group                                               !
!                                                                                                  !
!  See the LICENSE file for terms of usage and distribution.                                       !
!--------------------------------------------------------------------------------------------------!

#:include 'common.fypp'

!> Implements a simple repulsive potential between two atoms.
module dftbp_repsimple
  use dftbp_assert
  use dftbp_accuracy
  implicit none

  public :: TRepSimpleIn, ORepSimple, init
  public :: getCutoff, getEnergy, getEnergyDeriv
  private


  type :: TRepSimpleIn
    real(dp) :: zeff
    real(dp) :: alpha
    real(dp) :: kexp
    real(dp) :: rexp
    !> Cutoff of the last spline
    real(dp) :: cutoff
  end type TRepSimpleIn

  !> Contains the simple representation of a repulsive.
  type, extends(TRepSimpleIn) :: ORepSimple
    private
    !> Initialisation status
    logical :: tInit = .false.
  end type ORepSimple


  !> Initialises simple repulsive.
  interface init
    module procedure RepSimple_init
  end interface


  !> Returns cutoff of the repulsive.
  interface getCutoff
    module procedure RepSimple_getCutoff
  end interface


  !> Returns energy of the repulsive for a given distance.
  interface getEnergy
    module procedure RepSimple_getEnergy
  end interface


  !> Returns gradient of the repulsive for a given distance.
  interface getEnergyDeriv
    module procedure RepSimple_getEnergyDeriv
  end interface

contains


  !> Initialises spline repulsive.
  subroutine RepSimple_init(self, inp)
    !> Simple repulsive.
    type(ORepSimple), intent(out) :: self
    !> Input parameters for the spline repulsive.
    type(TRepSimpleIn), intent(in) :: inp

    @:ASSERT(.not. self%tInit)
    @:ASSERT(inp%cutoff >= 0.0_dp)

    self%TRepSimpleIn = inp

    self%tInit = .true.

  end subroutine RepSimple_init


  !> Returns cutoff of the repulsive.
  function RepSimple_getCutoff(self) result(cutoff)
    !> Simple repulsive.
    type(ORepSimple), intent(in) :: self
    !> Cutoff.
    real(dp) :: cutoff

    cutoff = self%cutoff

  end function RepSimple_getCutoff


  !> Returns energy of the repulsive for a given distance.
  subroutine RepSimple_getEnergy(self, res, rr)
    !> Simple repulsive.
    type(ORepSimple), intent(in) :: self
    !> repulsive contribution
    real(dp), intent(out) :: res
    !> Distance between interacting atoms.
    real(dp), intent(in) :: rr

    real(dp) :: t1, t2

    if (rr >= self%cutoff) then
      res = 0.0_dp
    else
      t1 = rr**self%kexp
      t2 = exp(-self%alpha*t1)/rr**self%rexp
      res = self%zeff * t2
    end if

  end subroutine RepSimple_getEnergy


  !> Returns gradient of the repulsive for a given distance.
  subroutine RepSimple_getEnergyDeriv(self, grad, xx, d2)
    !> Simple repulsive.
    type(ORepSimple), intent(in) :: self
    !> Resulting contribution
    real(dp), intent(out) :: grad(3)
    !> Actual vector between atoms
    real(dp), intent(in) :: xx(3)
    !> Second derivative in direction of xx, if needed.
    real(dp), intent(out), optional :: d2

    real(dp) :: rr, t1, t2, d1

    rr = sqrt(sum(xx**2))
    if (rr >= self%cutoff) then
      d1 = 0.0_dp
    else
      t1 = rr**self%kexp
      t2 = exp(-self%alpha*t1)/rr**self%rexp
      d1 = -self%zeff*t2/rr*(self%alpha*t1*self%kexp + self%rexp)
    end if
    grad = d1 * xx / rr

  end subroutine RepSimple_getEnergyDeriv

end module dftbp_repsimple
