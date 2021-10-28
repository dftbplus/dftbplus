!--------------------------------------------------------------------------------------------------!
!  DFTB+: general package for performing fast atomistic simulations                                !
!  Copyright (C) 2006 - 2021  DFTB+ developers group                                               !
!                                                                                                  !
!  See the LICENSE file for terms of usage and distribution.                                       !
!--------------------------------------------------------------------------------------------------!

#:include 'common.fypp'

!> Routines to calculate a Slater type orbital (STO)
module waveplot_slater

  use dftbp_common_accuracy, only : dp
  implicit none

  private
  public :: RealTessY, TRadialWavefunc, init
  save


  !> Data type containing radial wavefunction components
  type TRadialWavefunc

    !> Radial wavefunction components
    real(dp), allocatable :: rwf(:,:)

  end type TRadialWavefunc


  !> Initialises a RadialWavefunc
  interface init
    module procedure RadialWavefunc_init
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


  !> Initialises a RadialWavefunc.
  subroutine RadialWavefunc_init(this, rwf)

    !> RadialWavefunc instance to initialise
    type(TRadialWavefunc), intent(inout) :: this

    !> Radial wavefunction component
    real(dp), intent(in) :: rwf(:,:)

    this%rwf = rwf

  end subroutine RadialWavefunc_init

end module waveplot_slater
