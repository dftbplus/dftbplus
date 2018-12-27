!--------------------------------------------------------------------------------------------------!
!  DFTB+: general package for performing fast atomistic simulations                                !
!  Copyright (C) 2018  DFTB+ developers group                                                      !
!                                                                                                  !
!  See the LICENSE file for terms of usage and distribution.                                       !
!--------------------------------------------------------------------------------------------------!

!> Implementing the special functions erf(x), erfc(x) and exp(x*x) * erfc(x).
!>
!> The special functions erf(x), erfc(x) and exp(x * x) * erfc(x) are
!> implemented using the appropriate routine in NETLIB/SPECFUNC. The routines
!> have been converted to Fortran 2003. They can handle single and double
!> precision calls.
!>
!> Compared to iforts built in erf routine, the max. deviation is 4e-16 for
!> the double precision implementation.
module erfcalc

  !> wp: working precision, sp: real single, dp: real double
  use accuracy,  only : wp => dp, sp => rsp, dp => rdp
  implicit none
  private

  public :: erf, erfc, erfcx


  !> Evaluate erf()
  interface erfcalc_calc
    module procedure erfcalc_calcsingle, erfcalc_calcdouble
  end interface

contains


  !> Calculates the value of the error function.
  function erf(x)

    !> Function argument.
    real(wp), intent(in) :: x

    !> erf(x)
    real(wp) :: erf

    call erfcalc_calc(x, erf, 0)

  end function erf


  !> Calculates the value of the complementary error function.
  function erfc(x)

    !> Function argument.
    real(wp), intent(in) :: x

    !> erfc(x)
    real(wp) :: erfc

    call erfcalc_calc(x, erfc, 1)

  end function erfc


  !> Calculates the value of the function exp(x**2) * erfc(x)
  function erfcx(x)

    !> Function argument.
    real(wp), intent(in) :: x

    !> exp(x**2) * erfc(x)
    real(wp) :: erfcx

    call erfcalc_calc(x, erfcx, 2)

  end function erfcx


  !> Calculates the appropriate function in double precision.
  subroutine erfcalc_calcdouble(arg, res, jint)

    !> Where to evaluate the function (x).
    real(dp), intent(in) :: arg

    !> Result.
    real(dp), intent(out) :: res

    !> Function type: 1 - erf(x), 2 - erfc(x), 3 - exp(x**2)*erfc(x).
    integer, intent(in) :: jint

    real(dp), parameter :: xsmall = 1.11E-16_dp
    real(dp), parameter :: xbig = 26.543_dp
    real(dp), parameter :: xhuge = 6.71E7_dp
    real(dp), parameter :: xmax = 2.53E307_dp

    real(dp), parameter :: four = 4.0_dp
    real(dp), parameter :: one = 1.0_dp
    real(dp), parameter :: zero = 0.0_dp
    real(dp), parameter :: sqrpi = 5.6418958354775628695E-1_dp
    real(dp), parameter :: thresh = 0.46875E0_dp
    real(dp), parameter :: sixteen = 16.0_dp
    real(dp), parameter :: aa(5) = [&
        & 3.16112374387056560E00_dp, 1.13864154151050156E02_dp,&
        & 3.77485237685302021E02_dp, 3.20937758913846947E03_dp,&
        & 1.85777706184603153E-1_dp ]
    real(dp), parameter :: bb(4) = [&
        & 2.36012909523441209E01_dp, 2.44024637934444173E02_dp,&
        & 1.28261652607737228E03_dp, 2.84423683343917062E03_dp ]
    real(dp), parameter :: cc(9) = [&
        & 5.64188496988670089E-1_dp, 8.88314979438837594E00_dp,&
        & 6.61191906371416295E01_dp, 2.98635138197400131E02_dp,&
        & 8.81952221241769090E02_dp, 1.71204761263407058E03_dp,&
        & 2.05107837782607147E03_dp, 1.23033935479799725E03_dp,&
        & 2.15311535474403846E-8_dp ]
    real(dp), parameter :: dd(8) = [&
        & 1.57449261107098347E01_dp, 1.17693950891312499E02_dp,&
        & 5.37181101862009858E02_dp, 1.62138957456669019E03_dp,&
        & 3.29079923573345963E03_dp, 4.36261909014324716E03_dp,&
        & 3.43936767414372164E03_dp, 1.23033935480374942E03_dp ]
    real(dp), parameter :: pp(6) = [&
        & 3.05326634961232344E-1_dp, 3.60344899949804439E-1_dp,&
        & 1.25781726111229246E-1_dp, 1.60837851487422766E-2_dp,&
        & 6.58749161529837803E-4_dp, 1.63153871373020978E-2_dp ]
    real(dp), parameter :: qq(5) = [&
        & 2.56852019228982242E00_dp, 1.87295284992346047E00_dp,&
        & 5.27905102951428412E-1_dp, 6.05183413124413191E-2_dp,&
        & 2.33520497626869185E-3_dp ]

    integer :: ii
    real(dp) :: del, xx, xden, xnum, yy, ysq

    xx = arg
    yy = abs(xx)

    if (yy <= thresh) then
      !------------------------------------------------------------------
      !  evaluate  erf  for  |x| <= 0.46875
      !------------------------------------------------------------------
      ysq = zero
      if (yy > xsmall) ysq = yy * yy
      xnum = aa(5)*ysq
      xden = ysq
      do ii = 1, 3
        xnum = (xnum + aa(ii)) * ysq
        xden = (xden + bb(ii)) * ysq
      end do
      res = xx * (xnum + aa(4)) / (xden + bb(4))
      if (jint /= 0) res = one - res
      if (jint == 2) res = exp(ysq) * res
      return

    elseif (yy <= four) then
      !------------------------------------------------------------------
      !  evaluate  erfc  for 0.46875 <= |x| <= 4.0
      !------------------------------------------------------------------
      xnum = cc(9) * yy
      xden = yy
      do ii = 1, 7
        xnum = (xnum + cc(ii)) * yy
        xden = (xden + dd(ii)) * yy
      end do
      res = (xnum + cc(8)) / (xden + dd(8))
      if (jint /= 2) then
        ysq = aint(yy * sixteen) / sixteen
        del = (yy - ysq) * (yy + ysq)
        res = exp(-ysq * ysq) * exp(-del) * res
      end if
    else
      !------------------------------------------------------------------
      !  evaluate  erfc  for |x| > 4.0
      !------------------------------------------------------------------
      res = zero
      if (yy >= xbig) then
        if ((jint /= 2) .or. (yy >= xmax)) then
          call fixnegf_dp(res, ysq, del, yy, jint, xx)
          return
        end if
        if (yy >= xhuge) then
          res = sqrpi / yy
          call fixnegf_dp(res, ysq, del, yy, jint, xx)
          return
        end if
      end if
      ysq = one / (yy * yy)
      xnum = pp(6)*ysq
      xden = ysq
      do ii = 1, 4
        xnum = (xnum + pp(ii)) * ysq
        xden = (xden + qq(ii)) * ysq
      end do
      res = ysq *(xnum + pp(5)) / (xden + qq(5))
      res = (sqrpi -  res) / yy
      if (jint /= 2) then
        ysq = aint(yy * sixteen) / sixteen
        del = (yy - ysq) * (yy + ysq)
        res = exp(-ysq*ysq) * exp(-del) * res
      end if
    end if
    call fixnegf_dp(res, ysq, del, yy, jint, xx)

  end subroutine erfcalc_calcdouble

  !> fix up for negative argument, erf, etc.
  subroutine fixnegf_dp(res, ysq, del, yy, jint, xx)
    real(dp), intent(inout) :: res
    real(dp), intent(inout) :: ysq
    real(dp), intent(inout) :: del
    real(dp), intent(inout) :: yy
    integer, intent(in) :: jint
    real(dp), intent(in) :: xx

    real(dp), parameter :: zero = 0.0_dp
    real(dp), parameter :: half = 0.5_dp
    real(dp), parameter :: two = 2.0_dp
    real(dp), parameter :: sixteen = 16.0_dp
    real(dp), parameter :: xneg = -26.628_dp
    real(dp), parameter :: xinf = 1.79E308_dp

    if (jint == 0) then
      res = (half - res) + half
      if (xx < zero) res = -res
    else if (jint == 1) then
      if (xx < zero) res = two - res
    else
      if (xx < zero) then
        if (xx < xneg) then
          res = xinf
        else
          ysq = aint(xx * sixteen) / sixteen
          del = (xx - ysq) * (xx + ysq)
          yy = exp(ysq * ysq) * exp(del)
          res = (yy + yy) - res
        end if
      end if
    end if

  end subroutine fixnegf_dp

  !> Calculates the appropriate function in single precision.
  subroutine erfcalc_calcsingle(arg, res, jint)

    !> Where to evaluate the function (x).
    real(sp), intent(in) :: arg

    !> Result.
    real(sp), intent(out) :: res

    !> Function type: 1 - erf(x), 2 - erfc(x), 3 - exp(x**2)*erfc(x).
    integer, intent(in) :: jint

    real(sp), parameter :: xsmall = 5.96E-8_sp
    real(sp), parameter :: xbig = 9.194E0_sp
    real(sp), parameter :: xhuge = 2.90E3_sp
    real(sp), parameter :: xmax = 4.79E37_dp

    real(sp), parameter :: four = 4.0_sp
    real(sp), parameter :: one = 1.0_sp
    real(sp), parameter :: zero = 0.0_sp
    real(sp), parameter :: sqrpi = 5.6418958354775628695E-1_sp
    real(sp), parameter :: thresh = 0.46875E0_sp
    real(sp), parameter :: sixteen = 16.0_sp
    real(sp), parameter :: aa(5) = [&
        & 3.16112374387056560E00_sp, 1.13864154151050156E02_sp,&
        & 3.77485237685302021E02_sp, 3.20937758913846947E03_sp,&
        & 1.85777706184603153E-1_sp ]
    real(sp), parameter :: bb(4) = [&
        & 2.36012909523441209E01_sp, 2.44024637934444173E02_sp,&
        & 1.28261652607737228E03_sp, 2.84423683343917062E03_sp ]
    real(sp), parameter :: cc(9) = [&
        & 5.64188496988670089E-1_sp, 8.88314979438837594E00_sp,&
        & 6.61191906371416295E01_sp, 2.98635138197400131E02_sp,&
        & 8.81952221241769090E02_sp, 1.71204761263407058E03_sp,&
        & 2.05107837782607147E03_sp, 1.23033935479799725E03_sp,&
        & 2.15311535474403846E-8_sp ]
    real(sp), parameter :: dd(8) = [&
        & 1.57449261107098347E01_sp, 1.17693950891312499E02_sp,&
        & 5.37181101862009858E02_sp, 1.62138957456669019E03_sp,&
        & 3.29079923573345963E03_sp, 4.36261909014324716E03_sp,&
        & 3.43936767414372164E03_sp, 1.23033935480374942E03_sp ]
    real(sp), parameter :: pp(6) = [&
        & 3.05326634961232344E-1_sp, 3.60344899949804439E-1_sp,&
        & 1.25781726111229246E-1_sp, 1.60837851487422766E-2_sp,&
        & 6.58749161529837803E-4_sp, 1.63153871373020978E-2_sp ]
    real(sp), parameter :: qq(5) = [&
        & 2.56852019228982242E00_sp, 1.87295284992346047E00_sp,&
        & 5.27905102951428412E-1_sp, 6.05183413124413191E-2_sp,&
        & 2.33520497626869185E-3_sp ]

    integer :: ii
    real(sp) :: del, xx, xden, xnum, yy, ysq

    xx = arg
    yy = abs(xx)

    if (yy <= thresh) then
      !------------------------------------------------------------------
      !  evaluate  erf  for  |x| <= 0.46875
      !------------------------------------------------------------------
      ysq = zero
      if (yy > xsmall) ysq = yy * yy
      xnum = aa(5)*ysq
      xden = ysq
      do ii = 1, 3
        xnum = (xnum + aa(ii)) * ysq
        xden = (xden + bb(ii)) * ysq
      end do
      res = xx * (xnum + aa(4)) / (xden + bb(4))
      if (jint /= 0) res = one - res
      if (jint == 2) res = exp(ysq) * res
      return

    elseif (yy <= four) then
      !------------------------------------------------------------------
      !  evaluate  erfc  for 0.46875 <= |x| <= 4.0
      !------------------------------------------------------------------
      xnum = cc(9) * yy
      xden = yy
      do ii = 1, 7
        xnum = (xnum + cc(ii)) * yy
        xden = (xden + dd(ii)) * yy
      end do
      res = (xnum + cc(8)) / (xden + dd(8))
      if (jint /= 2) then
        ysq = aint(yy * sixteen) / sixteen
        del = (yy - ysq) * (yy + ysq)
        res = exp(-ysq * ysq) * exp(-del) * res
      end if
    else
      !------------------------------------------------------------------
      !  evaluate  erfc  for |x| > 4.0
      !------------------------------------------------------------------
      res = zero
      if (yy >= xbig) then
        if ((jint /= 2) .or. (yy >= xmax)) then
          call fixnegf_sp(res, ysq, del, yy, jint, xx)
          return
        end if
        if (yy >= xhuge) then
          res = sqrpi / yy
          call fixnegf_sp(res, ysq, del, yy, jint, xx)
          return
        end if
      end if
      ysq = one / (yy * yy)
      xnum = pp(6)*ysq
      xden = ysq
      do ii = 1, 4
        xnum = (xnum + pp(ii)) * ysq
        xden = (xden + qq(ii)) * ysq
      end do
      res = ysq *(xnum + pp(5)) / (xden + qq(5))
      res = (sqrpi -  res) / yy
      if (jint /= 2) then
        ysq = aint(yy * sixteen) / sixteen
        del = (yy - ysq) * (yy + ysq)
        res = exp(-ysq*ysq) * exp(-del) * res
      end if
    end if
    call fixnegf_sp(res, ysq, del, yy, jint, xx)

  end subroutine erfcalc_calcsingle

  !> fix up for negative argument, erf, etc.
  subroutine fixnegf_sp(res, ysq, del, yy, jint, xx)
    real(sp), intent(inout) :: res
    real(sp), intent(inout) :: ysq
    real(sp), intent(inout) :: del
    real(sp), intent(inout) :: yy
    integer, intent(in) :: jint
    real(sp), intent(in) :: xx

    real(sp), parameter :: zero = 0.0_sp
    real(sp), parameter :: half = 0.5_sp
    real(sp), parameter :: two = 2.0_sp
    real(sp), parameter :: sixteen = 16.0_sp
    real(sp), parameter :: xneg = -9.382E0_sp
    real(sp), parameter :: xinf = 3.40E+38_sp

    if (jint == 0) then
      res = (half - res) + half
      if (xx < zero) res = -res
    else if (jint == 1) then
      if (xx < zero) res = two - res
    else
      if (xx < zero) then
        if (xx < xneg) then
          res = xinf
        else
          ysq = aint(xx * sixteen) / sixteen
          del = (xx - ysq) * (xx + ysq)
          yy = exp(ysq * ysq) * exp(del)
          res = (yy + yy) - res
        end if
      end if
    end if
  end subroutine fixnegf_sp

end module erfcalc
