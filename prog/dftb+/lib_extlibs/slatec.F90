!--------------------------------------------------------------------------------------------------!
!  DFTB+: general package for performing fast atomistic simulations                                !
!  Copyright (C) 2006 - 2019  DFTB+ developers group                                               !
!                                                                                                  !
!  See the LICENSE file for terms of usage and distribution.                                       !
!--------------------------------------------------------------------------------------------------!

#:include 'common.fypp'

!> Interfaces for the SLATEC routines needed in DFTB+
module dftbp_slatec
  use dftbp_accuracy, only : rdp
  implicit none
  private

  public :: zbesi, dgamit, dgamic


  interface

    !> Compute a sequence of the Bessel functions I(a,z) for complex argument z and real nonnegative
    !> orders a=b,b+1, b+2,... where b>0.  A scaling option is available to help avoid overflow.
    subroutine zbesi(zr, zi, fnu, kode, n, cyr, cyi, nz, ierr)
      import rdp
      real(rdp), intent(in) :: zr
      real(rdp), intent(in) :: zi
      real(rdp), intent(in) :: fnu
      integer, intent(in) :: kode
      integer, intent(in) :: n
      real(rdp), dimension (n), intent(out) :: cyr
      real(rdp), dimension (n), intent(out) :: cyi
      integer, intent(out) :: nz
      integer, intent(out) :: ierr
    end subroutine zbesi


    !> Calculate Tricomi's form of the incomplete Gamma function.
    function dgamit (a, x)
      import rdp
      real(rdp) dgamit
      real(rdp), intent(in) :: x
      real(rdp), intent(in) :: a
    end function dgamit


    !> Calculate the complementary incomplete Gamma function.
    function dgamic (a, x)
      import rdp
      real(rdp) dgamic
      real(rdp), intent(in) :: x
      real(rdp), intent(in) :: a
    end function dgamic

  end interface

end module dftbp_slatec
