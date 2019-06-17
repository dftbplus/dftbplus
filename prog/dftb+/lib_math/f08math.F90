#:include 'common.fypp'

!> Contains replacements for some mathematical routines introduced in Fortran 2008.
!>
!> If the compiler does not provide those mathematical routines, preprocess the module
!> with the -DEMULATE_F08_MATH option.
!>
module dftbp_f08math
  use dftbp_accuracy, only : dp
  implicit none
  private

#:if EMULATE_F08_MATH  

  public :: norm2

contains

  pure function norm2(array)
    real(dp), intent(in) :: array(:)
    real(dp) :: norm2

    norm2 = sqrt(sum(array**2))

  end function norm2

#:endif
  
end module dftbp_f08math
