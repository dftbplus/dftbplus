!--------------------------------------------------------------------------------------------------!
!  DFTB+: general package for performing fast atomistic simulations                                !
!  Copyright (C) 2006 - 2024  DFTB+ developers group                                               !
!                                                                                                  !
!  See the LICENSE file for terms of usage and distribution.                                       !
!--------------------------------------------------------------------------------------------------!

#:include 'common.fypp'

!> Contains routines for Kramers-Kronig transformations
module dftbp_math_kramerskronig
  use dftbp_common_accuracy, only : dp
  use dftbp_common_constants, only : pi
  implicit none

  private
  public :: kki2r

contains

  !> Kramers-Kronig transformation of the imaginary part of a frequency dependent quantity into its
  !! real part
  subroutine kki2r(delta, lowerLimit, im, re, eta)

    !> Grid spacing for data
    real(dp), intent(in) :: delta

    !> Lower limit for integral (i.e. lowest value grid point)
    real(dp), intent(in) :: lowerLimit

    !> Imaginary part
    real(dp), intent(in) :: im(:,:)

    !> Real part
    real(dp), intent(out) :: re(:,:)

    !> Width of broadening
    complex(dp), intent(in) :: eta

    integer :: ii, jj, n
    !! Real values of grid points
    real(dp) :: er, ei

    @:ASSERT(all(shape(re) == shape(im)))

    n = size(im, dim=2)
    re(:,:) = 0.0_dp

    !$OMP PARALLEL DO&
    !$OMP& DEFAULT(SHARED) PRIVATE(jj, er, ei) SCHEDULE(RUNTIME) REDUCTION(+:re)
    do ii = 1, n
      er = lowerLimit + real(ii-1, dp) * delta
      do jj = 1, n
        ei = lowerLimit + real(jj-1, dp) * delta
        re(:,ii) = re(:,ii) + im(:,jj) * (2.0_dp*delta/pi)*real(ei/(ei**2 - er**2 + eta), dp)
      end do
    end do
    !$OMP END PARALLEL DO

  end subroutine kki2r

end module dftbp_math_kramerskronig
