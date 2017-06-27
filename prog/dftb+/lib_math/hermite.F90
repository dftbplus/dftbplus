!--------------------------------------------------------------------------------------------------!
!  DFTB+: general package for performing fast atomistic simulations                                !
!  Copyright (C) 2017  DFTB+ developers group                                                      !
!                                                                                                  !
!  See the LICENSE file for terms of usage and distribution.                                       !
!--------------------------------------------------------------------------------------------------!

#:include 'common.fypp'

!!* Contains routines relating to Hermite polynomials
!!* @todo Proper documentation, and trap overflow and underflows
module hermite
  use assert
  use accuracy, only : dp
  implicit none
contains
  !!* Calculate the Hermite polynomials up to a given order
  !!* @param h Resulting polynomials $H_n(x)$, starting from $H_0$ upwards
  !!* @param n Order of the polynomials to calculate up to
  !!* @param x Value to calculate the polynomials for
  !!* @desc Calculates all of the Hermite polynomials of x from
  !!* order 0 up to order n, using the recurrence relation, $H_{n+1}(x) =
  !!* 2\left( xH_{n}(x) - nH_{n-1}(x) \right)$, with $H_0(x) = 1$ and $H_0(x) =
  !!* 2x$.
  !!* @author B. Hourahine
  subroutine hx(h,n,x)
    real(dp), intent(out) :: h(0:)
    integer, intent(in) :: n
    real(dp), intent(in) :: x
    integer :: i
    @:ASSERT(n >=0 )
    @:ASSERT(size(H) > 0)
    @:ASSERT(size(H) <= n+1)
    h(:) = 0.0_dp
    h(0) = 1.0_dp
    if (n > 0) then
       h(1) = 2.0_dp * x
       do i = 2,n
          h(i) = 2.0_dp*(x*h(i-1) - real((i-1),dp)*h(i-2))
       end do
    end if
  end subroutine hx
end module hermite
