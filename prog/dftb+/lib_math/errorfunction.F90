!--------------------------------------------------------------------------------------------------!
!  DFTB+: general package for performing fast atomistic simulations                                !
!  Copyright (C) 2017  DFTB+ developers group                                                      !
!                                                                                                  !
!  See the LICENSE file for terms of usage and distribution.                                       !
!--------------------------------------------------------------------------------------------------!

!> Wrappers for the functions erf(x) and erfc(x).
!!
!! \details Based on the preprocessor settings, the error function is wrapped
!! differently:
!!
!! - no special definitions: the intrinsic error function is used (officially
!!   first available in the Fortran 2008 standard, but most F95/2003 compilers
!!   already implements this).
!! - EXTERNALERFC is defined: single precision and double precision external
!!   routines are expected (erf(x), erfc(x), derf(x), derfc(x)).
!! - INTERNALERFC is defined: erf(x) and erfc(x) are internally calculated
!!   by the code.
module errorfunction
  use accuracy
#ifdef INTERNALERFC
  use erfcalc, only: erf, erfc
#endif
  implicit none
  private

  public :: erfwrap, erfcwrap

#ifdef EXTERNALERFC
  interface exterf
    module procedure exterf_real, exterf_double
  end interface exterf
  
  interface exterfc
    module procedure exterfc_real, exterfc_double
  end interface exterfc

  interface
    function erf(xx)
      import rsp
      real(rsp), intent(in) :: xx
      real(rsp) :: erf
    end function erf

    function derf(xx)
      import rdp
      real(rdp), intent(in) :: xx
      real(rdp) :: derf
    end function derf
    
    function erfc(xx)
      import rsp
      real(rsp), intent(in) :: xx
      real(rsp) :: erfc
    end function erfc

    function derfc(xx)
      import rdp
      real(rdp), intent(in) :: xx
      real(rdp) :: derfc
    end function derfc
  end interface
#endif


contains

  !> Calculates the value of the error function.
  !! \param x  Function argument.
  !! \return erf(x)
  function erfwrap(xx) result(res)
    real(dp), intent(in) :: xx
    real(dp) :: res
    
#ifdef EXTERNALERFC    
    res = exterf(xx)
#else
    res = erf(xx)
#endif

  end function erfwrap


  !> Calculates the value of the complementary error function.
  !! \param x  Function argument.
  !! \return erf(x)
  function erfcwrap(xx) result(res)
    real(dp), intent(in) :: xx
    real(dp) :: res

#ifdef EXTERNALERFC    
    res = exterfc(xx)
#else
    res = erfc(xx)
#endif
    
  end function erfcwrap
  

#ifdef EXTERNALERFC
  
  !> Single precision external complementary error function routine.
  !! \param xx value to calculate erfc(x) of
  !! \return erfc(x)
  function exterfc_real(xx) result(res)
    real(rsp), intent(in) :: xx
    real(rsp) :: res

    res = erfc(xx)

  end function exterfc_real


  !> Double precision external complementary error function routine.
  !! \param xx value to calculate erfc(x) of
  !! \return erfc(x)
  function exterfc_double(xx) result(res)
    real(rdp), intent(in) :: xx
    real(rdp) :: res

    res = derfc(xx)

  end function exterfc_double


  !> Single precision external error function routine.
  !! \param xx value to calculate erf(x) of
  !! \return erf(x)
  function exterf_real(xx) result(res)
    real(rsp), intent(in) :: xx
    real(rsp) :: res

    res = erf(xx)

  end function exterf_real


  !> Double precision external error function routine.
  !! \param xx value to calculate erfc(x) of
  !! \return erf(x)
  function exterf_double(xx) result(res)
    real(rdp), intent(in) :: xx
    real(rdp) :: res

    res = derf(xx)

  end function exterf_double

#endif

end module errorfunction
