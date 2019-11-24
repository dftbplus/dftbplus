!--------------------------------------------------------------------------------------------------!
!  DFTB+: general package for performing fast atomistic simulations                                !
!  Copyright (C) 2006 - 2019  DFTB+ developers group                                               !
!                                                                                                  !
!  See the LICENSE file for terms of usage and distribution.                                       !
!--------------------------------------------------------------------------------------------------!

!> Contains routines relating to quaternion rotation in 3D
module dftbp_quaternions
  use dftbp_accuracy, only : dp
  use dftbp_simplealgebra, only : cross3
  implicit none

  private

  public :: quaternionProduct, quaternionInvert, quaternionRotate, quaternionConstruct, rotate3

  !> Product between two quaternions
  interface quaternionProduct
    module procedure quaternionProduct_
  end interface

  !> Inverse of a quaternion
  interface quaternionInvert
    module procedure quaternionInvert_
  end interface

  !> Rotation of a vector by a quaternion
  interface quaternionRotate
    module procedure quaternionRotate_
  end interface

  !> build a quaternion from angle + axis
  interface quaternionConstruct
    module procedure quaternionConstruct_
  end interface

  !> rotate a 3 vector by an angle around an axis
  interface rotate3
    module procedure rotate3_
  end interface

contains


  !> Product between two quaternions
  pure subroutine quaternionProduct_(res,a,b)

    !> Result of product
    real(dp), intent(out) :: res(4)

    !> First quaternion
    real(dp), intent(in)  :: a(4)

    !> second quaternion
    real(dp), intent(in)  :: b(4)

    real(dp) axb(3)

    axb(:) = cross3(a(2:4), b(2:4))

    res(1) = a(1)*b(1) - dot_product(a(2:),b(2:))
    res(2:) = a(1)*b(2:) + b(1)*a(2:) + axb

  end subroutine quaternionProduct_


  !> Form q^-1
  pure subroutine quaternionInvert_(res)

    !> Initial quaternion, overwritten by result
    real(dp), intent(inout) :: res(4)

    real(dp) :: mag2

    mag2 = sum(res**2)

    res(2:) = -res(2:)
    res = res / mag2

  end subroutine quaternionInvert_


  !> Rotates a 3 vector using a quaternion
  pure subroutine quaternionRotate_(vec,quat)

    !> Initial 3 vector, overwritten by result
    real(dp), intent(inout) :: vec(3)

    !> Quaternion
    real(dp), intent(in)  :: quat(4)

    real(dp) :: qTmp(4,3)

    qTmp = 0.0_dp
    qTmp(2:,1) = vec

    qTmp(:,2) = quat
    call  quaternionInvert(qTmp(:,2))

    call quaternionProduct(qTmp(:,3),qTmp(:,1),qTmp(:,2))
    call quaternionProduct(qTmp(:,1),quat,qTmp(:,3))

    vec = qTmp(2:,1)

  end subroutine quaternionRotate_


  !> Constructs a quaternion equivalent to rotation by an angle around an axis
  pure subroutine quaternionConstruct_(res,angle,axis)

    !> The quaternion
    real(dp), intent(out) :: res(4)

    !> Angle of rotation
    real(dp), intent(in)  :: angle

    !> 3 axis, does not require normalisation
    real(dp), intent(in)  :: axis(3)

    real(dp) :: norm(3)

    norm = axis / sqrt(sum(axis**2))

    res(1) = cos(0.5_dp * angle)
    res(2:) = sin(0.5_dp * angle) * norm

  end subroutine quaternionConstruct_


  !> Rotates a 3 vector by a specified angle around an axis
  pure subroutine rotate3_(x,angle,axis)

    !> Vector, overwritten when rotated
    real(dp), intent(inout) :: x(3)

    !> Angle in radians
    real(dp), intent(in) :: angle

    !> Specified axis (does not require normalisation, but must be non-zero)
    real(dp), intent(in) :: axis(3)

    real(dp) :: quat(4)

    call quaternionConstruct(quat,angle,axis)
    call quaternionRotate(x,quat)

  end subroutine rotate3_

end module dftbp_quaternions
