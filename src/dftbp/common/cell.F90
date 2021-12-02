!--------------------------------------------------------------------------------------------------!
!  DFTB+: general package for performing fast atomistic simulations                                !
!  Copyright (C) 2006 - 2021  DFTB+ developers group                                               !
!                                                                                                  !
!  See the LICENSE file for terms of usage and distribution.                                       !
!--------------------------------------------------------------------------------------------------!

!> Calculation of cell parameters
module dftbp_common_cell
  use dftbp_common_accuracy, only : dp
  use dftbp_common_constants, only : pi
  implicit none
  private

  public :: TCell, latticeToCell, hexagonalCAxes

  !> Cell parameters
  type :: TCell
    !> Cell lengths
    real(dp) :: length(3)
    !> Cell angles
    real(dp) :: angle(3)
  end type TCell

contains

  !> Convert lattice vectors to cell parameters
  pure function latticeToCell(lattice) result(cell)
    !> Lattice parameters
    real(dp),intent(in)  :: lattice(3,3)
    !> Cell parameters
    type(TCell) :: cell

    cell%length(1) = norm2(lattice(:, 1))
    cell%length(2) = norm2(lattice(:, 2))
    cell%length(3) = norm2(lattice(:, 3))

    cell%angle(1) = acos(dot_product(lattice(:, 2), lattice(:, 3)) &
        & / (cell%length(2)*cell%length(3)))
    cell%angle(2) = acos(dot_product(lattice(:, 1), lattice(:, 3)) &
        & / (cell%length(1)*cell%length(3)))
    cell%angle(3) = acos(dot_product(lattice(:, 1), lattice(:, 2)) &
        & / (cell%length(1)*cell%length(2)))

  end function latticeToCell

  !> Check whether we are dealing with a hexagonal lattice
  pure function hexagonalCAxes(lattice) result(axes)
    !> Lattice parameters
    real(dp), intent(in)  :: lattice(3,3)
    !> Is c axis in hexagonal lattice
    logical :: axes(3)

    type(TCell) :: cell
    real(dp), parameter :: eps = sqrt(epsilon(0.0_dp))
    real(dp), parameter :: a60 = 60.0_dp / 180.0_dp * pi, a120 = 120.0_dp / 180.0_dp * pi

    cell = latticeToCell(lattice)

    axes(1) = abs(cell%length(2) - cell%length(3)) < eps &
      & .and. abs(cell%angle(1) - a120) < eps &
      & .and. abs(cell%angle(2) - a60) < eps &
      & .and. abs(cell%angle(3) - a60) < eps

    axes(2) = abs(cell%length(1) - cell%length(3)) < eps &
      & .and. abs(cell%angle(1) - a60) < eps &
      & .and. abs(cell%angle(2) - a120) < eps &
      & .and. abs(cell%angle(3) - a60) < eps

    axes(3) = abs(cell%length(2) - cell%length(1)) < eps &
      & .and. abs(cell%angle(1) - a60) < eps &
      & .and. abs(cell%angle(2) - a60) < eps &
      & .and. abs(cell%angle(3) - a120) < eps

  end function hexagonalCAxes

end module dftbp_common_cell
