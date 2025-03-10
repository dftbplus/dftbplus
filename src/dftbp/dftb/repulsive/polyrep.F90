!--------------------------------------------------------------------------------------------------!
!  DFTB+: general package for performing fast atomistic simulations                                !
!  Copyright (C) 2006 - 2025  DFTB+ developers group                                               !
!                                                                                                  !
!  See the LICENSE file for terms of usage and distribution.                                       !
!--------------------------------------------------------------------------------------------------!

#:include 'common.fypp'

!> Implements a repulsive potential between two atoms represented by a polynomial of 9th degree
module dftbp_dftb_repulsive_polyrep
  use dftbp_common_accuracy, only : dp
  use dftbp_dftb_repulsive_pairrepulsive, only : TPairRepulsive
  implicit none

  private
  public :: TPolyRepInp, TPolyRep, TPolyRep_init


  !> Minimal power appearing in the polynomial
  integer, parameter :: powMin = 2

  !> Maximal power appearing in the polynomial
  integer, parameter :: powMax = 9

  !> smallest power between 1 and powMin
  integer, parameter :: powMin1 = max(powMin, 1)


  !> Initialisation type for TPolyRep
  type TPolyRepInp

    !> Polynomial coefficients
    real(dp) :: polyCoeffs(powMin:powMax)

    !> Cutoff distance
    real(dp) :: cutoff

  end type TPolyRepInp


  !> Contains the polynomial representation of a repulsive.
  type, extends(TPairRepulsive) :: TPolyRep
    private

    !> coefficients of the polynomial
    real(dp) :: polyCoeffs(powMin:powMax)

    !> Cutoff distance
    real(dp) :: cutoff

    !> initialised the repulsive
    logical :: tInit = .false.

  contains

    procedure :: getCutoff => TPolyRep_getCutoff
    procedure :: getValue => TPolyRep_getValue

  end type TPolyRep


  !> Structure constructor for TPolyRep
  interface TPolyRep
    module procedure TPolyRep_construct
  end interface TPolyRep


contains


  !> Initialises polynomial repulsive
  subroutine TPolyRep_init(this, inp)

    !> Polynomial repulsive
    type(TPolyRep), intent(out) :: this

    !> Input parameters for the polynomial repulsive
    type(TPolyRepInp), intent(in) :: inp

    @:ASSERT(inp%cutoff >= 0.0_dp)

    this%polyCoeffs(:) = inp%polyCoeffs(:)
    this%cutoff = inp%cutoff

  end subroutine TPolyRep_init


  !> Constructs polynomial repulsive
  function TPolyRep_construct(inp) result(this)

    !> Input parameters for the polynomial repulsive
    type(TPolyRepInp), intent(in) :: inp

    !> Polynomial repulsive
    type(TPolyRep) :: this

    call TPolyRep_init(this, inp)

  end function TPolyRep_construct


  !> Returns cutoff of the repulsive
  function TPolyRep_getCutoff(this) result(cutoff)

    !> self Polynomial repulsive
    class(TPolyRep), intent(in) :: this

    !> Cutoff
    real(dp) :: cutoff

    cutoff = this%cutoff

  end function TPolyRep_getCutoff


  !> Returns energy of the repulsive for a given distance
  subroutine TPolyRep_getValue(this, rr, energy, dEnergy, d2Energy)

    !> Polynomial repulsive
    class(TPolyRep), intent(in) :: this

    !> Distance between interacting atoms
    real(dp), intent(in) :: rr

    !> Energy
    real(dp), optional, intent(out) :: energy

    !> First derivative of the energy
    real(dp), optional, intent(out) :: dEnergy

    !> Second derivative of the energy
    real(dp), optional, intent(out) :: d2Energy

    real(dp) :: dr, xh
    integer :: ii

    @:ASSERT(rr >= 0.0_dp)

    if (rr >= this%cutoff) then
      if (present(energy)) energy = 0.0_dp
      if (present(dEnergy)) dEnergy = 0.0_dp
      if (present(d2Energy)) d2Energy = 0.0_dp
      return
    end if

    dr = this%cutoff - rr
    if (present(energy)) then
      energy = this%polyCoeffs(powMax)
      do ii = powMax - 1, powMin, -1
        energy = energy * dr + this%polyCoeffs(ii)
      end do
      do ii = powMin, 1, -1
        energy = energy * dr
      end do
    end if

    if (present(dEnergy)) then
      ! ((n * c_n * x + (n-1) * c_{n-1}) * x + (n-2) * c_{n-2}) * x + ... c_2
      xh = real(powMax, dp) * this%polyCoeffs(powMax)
      do ii = powMax - 1, powMin1, -1
        xh = xh * dr + real(ii, dp) * this%polyCoeffs(ii)
      end do
      do ii = powMin1, 2, -1
        xh = xh * dr
      end do
      ! polynomial is written as p(xcutoff - x), negative sign in first deriv
      dEnergy = -xh
    end if

    if (present(d2Energy)) then
      xh = real(powMax * (powMax -1), dp) * this%polyCoeffs(powMax)
      do ii = powMax - 1, powMin1, -1
        xh = xh * dr + real(ii * (ii - 1), dp) * this%polyCoeffs(ii)
      end do
      do ii = powMin1, 3, -1
        xh = xh * dr
      end do
      d2Energy = xh
    end if

  end subroutine TPolyRep_getValue

end module dftbp_dftb_repulsive_polyrep
