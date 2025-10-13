!--------------------------------------------------------------------------------------------------!
!  DFTB+: general package for performing fast atomistic simulations                                !
!  Copyright (C) 2006 - 2025  DFTB+ developers group                                               !
!                                                                                                  !
!  See the LICENSE file for terms of usage and distribution.                                       !
!--------------------------------------------------------------------------------------------------!

#:include 'common.fypp'
#:set pure = "" if defined('WITH_ASSERT') else "pure"

!> Real spherical harmonics up to l=4.
!! See also: spharmonics.cuh
module dftbp_wavegrid_basis_spharmonics
  use dftbp_common_accuracy, only : dp
  implicit none

  private

  public :: realTessY

  ! Normalisation Constants labeled C_l_m
  ! l = 0
  real(dp), parameter :: C0_0 = 0.28209479177387814_dp ! 1/2 * sqrt(1/pi)

  ! l = 1
  real(dp), parameter :: C1 = 0.4886025119029199_dp ! m =-1,0,1  1/2 * sqrt(3/pi)

  ! l = 2
  real(dp), parameter :: C2_0         = 0.31539156525252005_dp  ! m=0        1/4 * sqrt(5/pi)
  real(dp), parameter :: C2_abs1_neg2 = 1.0925484305920792_dp   ! m=-2,-1,1  1/2 * sqrt(15/pi)
  real(dp), parameter :: C2_2         = 0.5462742152960396_dp   ! m=2        1/4 * sqrt(15/pi)

  ! l = 3
  real(dp), parameter :: C3_0    = 0.3731763325901154_dp   ! m=0     1/4 * sqrt(7/pi)
  real(dp), parameter :: C3_abs1 = 0.4570457994644658_dp   ! m=-1,1  1/4 * sqrt(21/(2*pi))
  real(dp), parameter :: C3_neg2 = 2.890611442640554_dp    ! m=-2    1/2 * sqrt(105/pi)
  real(dp), parameter :: C3_2    = 1.445305721320277_dp    ! m=2     1/4 * sqrt(105/pi)
  real(dp), parameter :: C3_abs3 = 0.5900435899266435_dp   ! m=-3,3  1/4 * sqrt(35/(2*pi))

  ! l = 4
  real(dp), parameter :: C4_0    = 0.10578554691520431_dp  ! m=0     3/16 * sqrt(1/pi)
  real(dp), parameter :: C4_abs1 = 0.6690465435572892_dp   ! m=-1,1  3/4 * sqrt(5/(2*pi))
  real(dp), parameter :: C4_neg2 = 0.9461746957575601_dp   ! m=-2    3/4 * sqrt(5/pi)
  real(dp), parameter :: C4_2    = 0.47308734787878004_dp  ! m=2     3/8 * sqrt(5/pi)
  real(dp), parameter :: C4_abs3 = 1.7701307697799304_dp   ! m=-3,3  3/4 * sqrt(35/(2*pi))
  real(dp), parameter :: C4_neg4 = 2.5033429417967046_dp   ! m=-4    3/4 * sqrt(35/pi)
  real(dp), parameter :: C4_4    = 0.6258357354491761_dp   ! m=4     3/16 * sqrt(35/pi)

contains

  !> Computes real tesseral spherical harmonics Y_lm(r) up to l=4 without performing divisions.
  !> All Definitions in accordance with:
  !> https://en.wikipedia.org/wiki/Table_of_spherical_harmonics#Real_spherical_harmonics
  !> See also (general set):
  !> https://winter.group.shef.ac.uk/orbitron/atomic_orbitals/4f/4f_equations.html
  ${pure}$ function realTessY(l, m, coord, inv_r) result(rty)

    !> Orbital quantum number (0 <= l <= 4)
    integer, intent(in) :: l
    !> Magnetic quantum number
    integer, intent(in) :: m
    !> Coordinate vector (x, y, z) where the value should be calculated
    real(dp), intent(in) :: coord(:)
    !> Pre-calculated 1/r. (At origin, pass 0 to avoid NaNs)
    real(dp), intent(in) :: inv_r

    !> The value of the real spherical harmonic.
    real(dp) :: rty

    real(dp) :: x_r, y_r, z_r
    real(dp) :: x_r2, y_r2, z_r2

    @:ASSERT(l >= 0 .and. l <= 4)
    @:ASSERT(abs(m) <= l)
    @:ASSERT(size(coord) == 3)

    x_r = coord(1) * inv_r
    y_r = coord(2) * inv_r
    z_r = coord(3) * inv_r

    rty = 0.0_dp

    select case (l)
      case (0) ! s orbital
        rty = C0_0
      case (1) ! p orbitals
        select case (m)
          case (-1); rty = C1 * y_r ! p_y
          case (0);  rty = C1 * z_r ! p_z
          case (1);  rty = C1 * x_r ! p_x
        end select
      case (2) ! d orbitals
        x_r2 = x_r * x_r
        y_r2 = y_r * y_r
        z_r2 = z_r * z_r
        select case (m)
          case (-2); rty = C2_abs1_neg2 * (x_r * y_r)         ! d_xy
          case (-1); rty = C2_abs1_neg2 * (y_r * z_r)         ! d_yz
          case (0);  rty = C2_0 * (3.0_dp * z_r2 - 1.0_dp)    ! d_z^2
          case (1);  rty = C2_abs1_neg2 * (x_r * z_r)         ! d_xz
          case (2);  rty = C2_2 * (x_r2 - y_r2)               ! d_x^2-y^2
        end select
      case (3) ! f orbitals
        x_r2 = x_r * x_r
        y_r2 = y_r * y_r
        z_r2 = z_r * z_r
        select case (m)
          case (-3); rty = C3_abs3 * y_r * (3.0_dp * x_r2 - y_r2)    ! f_y(3x^2-y^2)
          case (-2); rty = C3_neg2 * (x_r * y_r * z_r)               ! f_xyz
          case (-1); rty = C3_abs1 * y_r * (5.0_dp * z_r2 - 1.0_dp)  ! f_y(5z^2-r^2) -> f_yz^2
          case (0);  rty = C3_0 * z_r * (5.0_dp * z_r2 - 3.0_dp)     ! f_z(5z^2-3r^2) -> f_z^3
          case (1);  rty = C3_abs1 * x_r * (5.0_dp * z_r2 - 1.0_dp)  ! f_x(5z^2-r^2) -> f_xz^2
          case (2);  rty = C3_2 * z_r * (x_r2 - y_r2)                ! f_z(x^2-y^2)
          case (3);  rty = C3_abs3 * x_r * (x_r2 - 3.0_dp * y_r2)    ! f_x(x^2-3y^2)
        end select
      case (4) ! g orbitals
        x_r2 = x_r * x_r
        y_r2 = y_r * y_r
        z_r2 = z_r * z_r
        select case (m)
          case (-4); rty = C4_neg4 * (x_r * y_r * (x_r2 - y_r2))
          case (-3); rty = C4_abs3 * (z_r * y_r * (3.0_dp * x_r2 - y_r2))
          case (-2); rty = C4_neg2 * (x_r * y_r * (7.0_dp * z_r2 - 1.0_dp))
          case (-1); rty = C4_abs1 * (z_r * y_r * (7.0_dp * z_r2 - 3.0_dp))
          case (0);  rty = C4_0 * (35.0_dp * z_r2 * z_r2 - 30.0_dp * z_r2 + 3.0_dp)
          case (1);  rty = C4_abs1 * (z_r * x_r * (7.0_dp * z_r2 - 3.0_dp))
          case (2);  rty = C4_2 * ((x_r2 - y_r2) * (7.0_dp * z_r2 - 1.0_dp))
          case (3);  rty = C4_abs3 * (z_r * x_r * (x_r2 - 3.0_dp * y_r2))
          case (4);  rty = C4_4 * (x_r2 * x_r2 - 6.0_dp * x_r2 * y_r2 + y_r2 * y_r2)
        end select
    end select

  end function realTessY


end module dftbp_wavegrid_basis_spharmonics
