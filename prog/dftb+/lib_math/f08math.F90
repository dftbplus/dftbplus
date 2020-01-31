!--------------------------------------------------------------------------------------------------!
!  DFTB+: general package for performing fast atomistic simulations                                !
!  Copyright (C) 2006 - 2020  DFTB+ developers group                                               !
!                                                                                                  !
!  See the LICENSE file for terms of usage and distribution.                                       !
!--------------------------------------------------------------------------------------------------!

#:include 'common.fypp'

!> Contains replacements for some mathematical routines introduced in Fortran 2008.
!>
!> If the compiler does not provide these routines, preprocess the module with the
!> -DEMULATE_F08_MATH option.
!>
module dftbp_f08math
  use dftbp_accuracy, only : dp
  implicit none
  private

#:if EMULATE_F08_MATH

  public :: norm2

contains

  !> Calculate the l2 norm of a 1 dimensional array
  pure function norm2(array)

    !> Array
    real(dp), intent(in) :: array(:)

    !> Resulting norm for the array
    real(dp) :: norm2

    norm2 = sqrt(sum(array**2))

  end function norm2

#:endif

end module dftbp_f08math
