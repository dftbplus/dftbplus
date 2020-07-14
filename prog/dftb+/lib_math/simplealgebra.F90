!--------------------------------------------------------------------------------------------------!
!  DFTB+: general package for performing fast atomistic simulations                                !
!  Copyright (C) 2006 - 2020  DFTB+ developers group                                               !
!                                                                                                  !
!  See the LICENSE file for terms of usage and distribution.                                       !
!--------------------------------------------------------------------------------------------------!

#:include 'common.fypp'

!> Simple algebraic stuff for cases where LAPACK would be overkill
module dftbp_simplealgebra
  use dftbp_assert
  use dftbp_accuracy
  implicit none
  private

  public :: cross3, determinant33, derivDeterminant33, invert33

contains

  !> Cross product
  pure function cross3(v1, v2) result(res)

    !> Resulting vector.
    real(dp) :: res(3)

    !> First vector.
    real(dp), intent(in) :: v1(:)

    !> Second vector.
    real(dp), intent(in) :: v2(:)

    res(1) = v1(2) * v2(3) - v1(3) * v2(2)
    res(2) = v1(3) * v2(1) - v1(1) * v2(3)
    res(3) = v1(1) * v2(2) - v1(2) * v2(1)

  end function cross3


  !> Signed determinant of a 3x3 matrix
  function  determinant33(matrix)

    !> The matrix for which to calculate the determinant.
    real(dp), intent(in) :: matrix(:,:)

    !> Resulting det(matrix)
    real(dp) :: determinant33

    real(dp) :: tmp

    @:ASSERT(all(shape(matrix) == [3, 3]))

    tmp = matrix(1, 1) &
        &* (matrix(2, 2) * matrix(3, 3) - matrix(3, 2) * matrix(2, 3))
    tmp = tmp - matrix(1, 2) &
        &* (matrix(2, 1) * matrix(3, 3) - matrix(3, 1) * matrix(2, 3))
    tmp = tmp + matrix(1, 3) &
        &* (matrix(2, 1) * matrix(3, 2) - matrix(3, 1) * matrix(2, 2))

    determinant33 = tmp

  end function determinant33


  !> Derivative of determinant of a 3x3 matrix
  subroutine  derivDeterminant33(deriv,matrix)

    !> derivative of the determinant
    real(dp), intent(out) :: deriv(:, :)

    !> The matrix from which to calculate the determinant.
    real(dp), intent(in) :: matrix(:, :)

    deriv(1,1) = matrix(2, 2) * matrix(3, 3) - matrix(3, 2) * matrix(2, 3)
    deriv(1,2) = matrix(2, 3) * matrix(3, 1) - matrix(3, 3) * matrix(2, 1)
    deriv(1,3) = matrix(2, 1) * matrix(3, 2) - matrix(3, 1) * matrix(2, 2)
    deriv(2,1) = matrix(1, 3) * matrix(3, 2) - matrix(1, 2) * matrix(3, 3)
    deriv(2,2) = matrix(1, 1) * matrix(3, 3) - matrix(1, 3) * matrix(3, 1)
    deriv(2,3) = matrix(1, 2) * matrix(3, 1) - matrix(1, 1) * matrix(3, 2)
    deriv(3,1) = matrix(1, 2) * matrix(2, 3) - matrix(1, 3) * matrix(2, 2)
    deriv(3,2) = matrix(1, 3) * matrix(2, 1) - matrix(1, 1) * matrix(2, 3)
    deriv(3,3) = matrix(1, 1) * matrix(2, 2) - matrix(1, 2) * matrix(2, 1)

    deriv = deriv * sign(1.0_dp,determinant33(matrix))

  end subroutine derivDeterminant33


  !> Inverts a 3x3 matrix
  subroutine invert33(inverted, orig, optDet)

    !> Contains the inverted matrix on return.
    real(dp), intent(out) :: inverted(:, :)

    !> Matrix to invert.
    real(dp), intent(in) :: orig(:, :)

    !> Determinant of the matrix, if already known.
    real(dp), intent(in), optional :: optDet

    real(dp) :: det

    @:ASSERT(all(shape(inverted) == [3, 3]))
    @:ASSERT(all(shape(orig) == [3, 3]))

    if (present(optDet)) then
      det = optDet
    else
      det = determinant33(orig)
    end if

    inverted(1, 1) = -orig(2, 3) * orig(3, 2) + orig(2, 2) * orig(3, 3)
    inverted(2, 1) =  orig(2, 3) * orig(3, 1) - orig(2, 1) * orig(3, 3)
    inverted(3, 1) = -orig(2, 2) * orig(3, 1) + orig(2, 1) * orig(3, 2)

    inverted(1, 2) =  orig(1, 3) * orig(3, 2) - orig(1, 2) * orig(3, 3)
    inverted(2, 2) = -orig(1, 3) * orig(3, 1) + orig(1, 1) * orig(3, 3)
    inverted(3, 2) =  orig(1, 2) * orig(3, 1) - orig(1, 1) * orig(3, 2)

    inverted(1, 3) = -orig(1, 3) * orig(2, 2) + orig(1, 2) * orig(2, 3)
    inverted(2, 3) =  orig(1, 3) * orig(2, 1) - orig(1, 1) * orig(2, 3)
    inverted(3, 3) = -orig(1, 2) * orig(2, 1) + orig(1, 1) * orig(2, 2)
    inverted = inverted / det

  end subroutine invert33

end module dftbp_simplealgebra
