!--------------------------------------------------------------------------------------------------!
!  DFTB+: general package for performing fast atomistic simulations                                !
!  Copyright (C) 2017  DFTB+ developers group                                                      !
!                                                                                                  !
!  See the LICENSE file for terms of usage and distribution.                                       !
!--------------------------------------------------------------------------------------------------!

#:include 'common.fypp'

!> Contains routines to calculate matrix determinants
module determinant
  use accuracy
  use lapackroutines
  implicit none

  private

  interface det
    module procedure det_real
    module procedure det_cmplx
  end interface det

  public :: det

contains

  !> Determinant of a real matrix, matrix destroyed in process
  function det_real(A) result(det)

    !> The matrix
    real(dp), intent(inout) :: A(:,:)

    !> resulting determinant
    real(dp) :: det

    integer, allocatable  :: ipiv(:)
    integer :: ii, n, exponent

    n = minval(shape(A))
    allocate(ipiv(n))

    call getrf(A,ipiv)

    det = 1.0_dp
    exponent = 0
    do ii = 1, n
      if (ipiv(ii).ne.ii) then
        det = -det * A(ii,ii)
      else
        det = det * A(ii,ii)
      end if
      if (det == 0.0_dp) then
        return
      end if
      do while (abs(det) > 2.0_dp)
        det = det / 2.0_dp
        exponent = exponent + 1
      end do
      do while (abs(det) < 0.5_dp)
        det = det * 2.0_dp
        exponent = exponent - 1
      end do
    end do
    det = det * 2.0_dp ** exponent

  end function det_real

  !> Determinant of a complex matrix, matrix destroyed in process
  function det_cmplx(A) result(det)

    !> The matrix
    complex(dp), intent(inout) :: A(:,:)

    !> resulting determinant
    complex(dp) :: det

    integer, allocatable  :: ipiv(:)
    integer :: ii, n, exponent

    n = minval(shape(A))
    allocate(ipiv(n))

    call getrf(A,ipiv)

    det = cmplx(1,0,dp)
    exponent = 0
    do ii = 1, n
      if (ipiv(ii).ne.ii) then
        det = -det * A(ii,ii)
      else
        det = det * A(ii,ii)
      end if
      if (det == 0.0_dp) then
        return
      end if
      do while (abs(det) > 2.0_dp)
        det = det / 2.0_dp
        exponent = exponent + 1
      end do
      do while (abs(det) < 0.5_dp)
        det = det * 2.0_dp
        exponent = exponent - 1
      end do
    end do
    det = det * 2.0_dp ** exponent

  end function det_cmplx

end module determinant
