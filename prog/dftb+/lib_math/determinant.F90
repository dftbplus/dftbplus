!--------------------------------------------------------------------------------------------------!
!  DFTB+: general package for performing fast atomistic simulations                                !
!  Copyright (C) 2017  DFTB+ developers group                                                      !
!                                                                                                  !
!  See the LICENSE file for terms of usage and distribution.                                       !
!--------------------------------------------------------------------------------------------------!

#:include 'common.fypp'

#:set FLAVORS = [('Real', 'real', 'real'), ('Cmplx', 'complex', 'cmplx')]

!> Contains routines to calculate matrix determinants
module determinant
  use accuracy
  use lapackroutines
  implicit none
  private

  public :: det

  interface det
  #:for SUFFIX, _, _ in FLAVORS
    module procedure det${SUFFIX}$
  #:endfor
  end interface det


contains

#:for SUFFIX, TYPE, CONVERT in

  !> Determinant of a real matrix, matrix destroyed in process
  function det${SUFFIX}$(A) result(det)

    !> The matrix
    ${TYPE}$(dp), intent(inout) :: A(:,:)

    !> resulting determinant
    ${TYPE}$(dp) :: det

    integer, allocatable  :: ipiv(:)
    integer :: ii, n, exponent

    n = minval(shape(A))
    allocate(ipiv(n))

    call getrf(A,ipiv)

    det = ${CONVERT}$(1, kind=dp)
    exponent = 0
    do ii = 1, n
      if (ipiv(ii) /= ii) then
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

  end function det${SUFFIX}$

#:endfor

end module determinant
