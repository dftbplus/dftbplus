!--------------------------------------------------------------------------------------------------!
!  DFTB+: general package for performing fast atomistic simulations                                !
!  Copyright (C) 2006 - 2020  DFTB+ developers group                                               !
!                                                                                                  !
!  See the LICENSE file for terms of usage and distribution.                                       !
!--------------------------------------------------------------------------------------------------!

#:include 'common.fypp'

!> Implementing the special functions erf(x), erfc(x) and exp(x*x) * erfc(x).
!>
!> The special functions erf(x), erfc(x) and exp(x * x) * erfc(x) are
!> implemented using the appropriate routine in NETLIB/SPECFUNC. The routines
!> have been converted to Fortran 2003. They can handle single and double
!> precision calls.
!>
!> Compared to iforts built in erf routine, the max. deviation is 4e-16 for
!> the double precision implementation.
module dftbp_erfcalc

  !> wp: working precision, sp: real single, dp: real double
  use dftbp_accuracy,  only : wp => dp, sp => rsp, dp => rdp
  implicit none
  private

#:if INTERNAL_ERFC

  public :: erf, erfc

  !> Evaluate erf()
  interface erfcalc_calc
    module procedure erfcalc_calcsingle, erfcalc_calcdouble
  end interface

contains


  !> Calculates the value of the error function.
  elemental function erf(x)

    !> Function argument.
    real(wp), intent(in) :: x

    !> erf(x)
    real(wp) :: erf

    erf = erfcalc_calc(x, 0)

  end function erf


  !> Calculates the value of the complementary error function.
  elemental function erfc(x)

    !> Function argument.
    real(wp), intent(in) :: x

    !> erfc(x)
    real(wp) :: erfc

    erfc = erfcalc_calc(x, 1)

  end function erfc

#:for VC, LABEL in [('sp', 'single'), ('dp', 'double')]

  !> Calculates the appropriate function in double precision.
  elemental function erfcalc_calc${LABEL}$(arg, jint) result(res)

    !> Where to evaluate the function (x).
    real(${VC}$), intent(in) :: arg

    !> Function type: 1 - erf(x), 2 - erfc(x)
    integer, intent(in) :: jint

    !> Result
    real(${VC}$) :: res

    real(${VC}$), parameter :: xsmall = 1.11E-16_${VC}$
    real(${VC}$), parameter :: xbig = 26.543_${VC}$
    real(${VC}$), parameter :: xhuge = 6.71E7_${VC}$

    ! see https://www.netlib.org/specfun/erf :
  #:if VC == 'sp'
    real(sp), parameter :: xmax = 4.79E37_sp
  #:else
    real(dp), parameter :: xmax = 2.53E307_dp
  #:endif

    real(${VC}$), parameter :: four = 4.0_${VC}$
    real(${VC}$), parameter :: one = 1.0_${VC}$
    real(${VC}$), parameter :: zero = 0.0_${VC}$
    real(${VC}$), parameter :: sqrpi = 5.6418958354775628695E-1_${VC}$
    real(${VC}$), parameter :: thresh = 0.46875E0_${VC}$
    real(${VC}$), parameter :: sixteen = 16.0_${VC}$
    real(${VC}$), parameter :: aa(5) = [&
        & 3.16112374387056560E00_${VC}$, 1.13864154151050156E02_${VC}$,&
        & 3.77485237685302021E02_${VC}$, 3.20937758913846947E03_${VC}$,&
        & 1.85777706184603153E-1_${VC}$ ]
    real(${VC}$), parameter :: bb(4) = [&
        & 2.36012909523441209E01_${VC}$, 2.44024637934444173E02_${VC}$,&
        & 1.28261652607737228E03_${VC}$, 2.84423683343917062E03_${VC}$ ]
    real(${VC}$), parameter :: cc(9) = [&
        & 5.64188496988670089E-1_${VC}$, 8.88314979438837594E00_${VC}$,&
        & 6.61191906371416295E01_${VC}$, 2.98635138197400131E02_${VC}$,&
        & 8.81952221241769090E02_${VC}$, 1.71204761263407058E03_${VC}$,&
        & 2.05107837782607147E03_${VC}$, 1.23033935479799725E03_${VC}$,&
        & 2.15311535474403846E-8_${VC}$ ]
    real(${VC}$), parameter :: dd(8) = [&
        & 1.57449261107098347E01_${VC}$, 1.17693950891312499E02_${VC}$,&
        & 5.37181101862009858E02_${VC}$, 1.62138957456669019E03_${VC}$,&
        & 3.29079923573345963E03_${VC}$, 4.36261909014324716E03_${VC}$,&
        & 3.43936767414372164E03_${VC}$, 1.23033935480374942E03_${VC}$ ]
    real(${VC}$), parameter :: pp(6) = [&
        & 3.05326634961232344E-1_${VC}$, 3.60344899949804439E-1_${VC}$,&
        & 1.25781726111229246E-1_${VC}$, 1.60837851487422766E-2_${VC}$,&
        & 6.58749161529837803E-4_${VC}$, 1.63153871373020978E-2_${VC}$ ]
    real(${VC}$), parameter :: qq(5) = [&
        & 2.56852019228982242E00_${VC}$, 1.87295284992346047E00_${VC}$,&
        & 5.27905102951428412E-1_${VC}$, 6.05183413124413191E-2_${VC}$,&
        & 2.33520497626869185E-3_${VC}$ ]

    integer :: ii
    real(${VC}$) :: del, xx, xden, xnum, yy, ysq

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
          call fixnegf_${VC}$(res, ysq, del, yy, jint, xx)
          return
        end if
        if (yy >= xhuge) then
          res = sqrpi / yy
          call fixnegf_${VC}$(res, ysq, del, yy, jint, xx)
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
    call fixnegf_${VC}$(res, ysq, del, yy, jint, xx)

  end function erfcalc_calc${LABEL}$


  !> fix up for negative argument, erf, etc.
  pure subroutine fixnegf_${VC}$(res, ysq, del, yy, jint, xx)
    real(${VC}$), intent(inout) :: res
    real(${VC}$), intent(inout) :: ysq
    real(${VC}$), intent(inout) :: del
    real(${VC}$), intent(inout) :: yy
    integer, intent(in) :: jint
    real(${VC}$), intent(in) :: xx

    real(${VC}$), parameter :: zero = 0.0_${VC}$
    real(${VC}$), parameter :: half = 0.5_${VC}$
    real(${VC}$), parameter :: two = 2.0_${VC}$
    real(${VC}$), parameter :: sixteen = 16.0_${VC}$
    real(${VC}$), parameter :: xneg = -26.628_${VC}$
    real(${VC}$), parameter :: xinf = 1.79E308_${VC}$

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

  end subroutine fixnegf_${VC}$

#:endfor

#:endif

end module dftbp_erfcalc
