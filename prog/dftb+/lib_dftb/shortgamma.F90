!--------------------------------------------------------------------------------------------------!
!  DFTB+: general package for performing fast atomistic simulations                                !
!  Copyright (C) 2006 - 2020  DFTB+ developers group                                               !
!                                                                                                  !
!  See the LICENSE file for terms of usage and distribution.                                       !
!--------------------------------------------------------------------------------------------------!

#:include "error.fypp"

!> Contains functions to calculate the short-ranged part of gamma, and the distance beyond which it
!> becomes negligible.
!>
!> Note, assumes Hubbard U values are above minHubTol and shell resolved U differ by more than
!> minHubDiff
module dftbp_shortgamma
  use dftbp_accuracy, only : dp, minHubTol, minHubDiff, tolShortGamma, minShortGamma, tolSameDist
  implicit none

  private

  public :: expGamma, expGammaDamped, expGammaPrime, expGammaDampedPrime
  public :: expGammaCutoff

contains


  !> Determines the cut off where the short range part goes to zero (a small constant)
  subroutine expGammaCutoff(cutOff, U2, U1, minValue, iErr)

    !> Returned cut off distance
    real(dp), intent(out) :: cutOff

    !> Hubbard U value in a.u.
    real(dp), intent(in) :: U2

    !> Hubbard U value in a.u.
    real(dp), intent(in) :: U1

    !> value below which the short range contribution is considered negligible. If not set this
    !> comes from a constant in the precision module.
    real(dp), intent(in), optional :: minValue

    !> Error return (if used)
    integer, intent(out), optional :: iErr

    real(dp) :: cutValue
    real(dp) :: rab, MaxGamma, MinGamma, lowerGamma, gamma
    real(dp) :: cut, MaxCutOff, MinCutOff

    cutoff = 0.0_dp

    if (present(minValue)) then
      cutValue = minValue
    else
      cutValue = minShortGamma
    end if

    if (cutValue < tolShortGamma) then
      @:FORMATTED_ERROR_HANDLING(iErr, -1, "(A,F12.6)",&
          & "Failure in determining short-range cut-off -ve cutoff negative :", cutValue)
    else if (U1 < minHubTol .or. U2 < minHubTol) then
      @:FORMATTED_ERROR_HANDLING(iErr, -1, "(A,2F12.6)",&
          & "Failure in short-range gamma, U too small :", U1, U2)
    end if

    rab = 1.0_dp
    do while(expGamma(rab,U2,U1) > cutValue)
      rab = 2.0_dp*rab
    end do
    if (rab < 2.0_dp) then
      @:FORMATTED_ERROR_HANDLING(iErr, -1, "(A,F10.6)",&
          & "Failure in short-range gamma cut-off : requested tolerance too large :", cutValue)
    end if
    ! bisection search for the value where the contribution drops below cutValue
    MinCutOff = rab + 0.1_dp
    MaxCutOff = 0.5_dp * rab - 0.1_dp
    maxGamma = expGamma(MaxCutOff,U2,U1)
    minGamma = expGamma(MinCutOff,U2,U1)
    lowerGamma =  expGamma(MinCutOff,U2,U1)
    cut = MaxCutOff + 0.1_dp
    gamma =  expGamma(cut,U2,U1)
    do while ((gamma-lowerGamma) > tolShortGamma)
      MaxCutOff = 0.5_dp*(cut + MinCutOff)
      if ((maxGamma >= minGamma) .eqv. (cutValue >= &
          & expGamma(MaxCutOff,U2,U1))) then
        MinCutOff = MaxCutOff
        lowerGamma =  expGamma(MinCutOff,U2,U1)
      else
        cut = MaxCutOff
        gamma =  expGamma(cut,U2,U1)
      end if
    end do
    cutoff = MinCutOff

  end subroutine expGammaCutoff


  !> Determines the value of the short range contribution to gamma with the exponential form
  pure function expGamma(rab, Ua, Ub)

    !> separation of sites a and b
    real(dp), intent(in) :: rab

    !> Hubbard U for site a
    real(dp), intent(in) :: Ua

    !> Hubbard U for site b
    real(dp), intent(in) :: Ub

    !> contribution
    real(dp) :: expGamma

    real(dp) :: tauA, tauB, tauMean

    ! 16/5 * U, see review papers / theses
    tauA = 3.2_dp*Ua
    tauB = 3.2_dp*Ub
    if (rab < tolSameDist) then
      ! on-site case with R~0
      if (abs(Ua - Ub) < MinHubDiff) then
        ! same Hubbard U values, onsite , NOTE SIGN CHANGE!
        expGamma = -0.5_dp*(Ua + Ub)
      else
        ! Ua /= Ub Hubbard U values - limiting case, NOTE SIGN CHANGE!
        expGamma = &
            & -0.5_dp*((tauA*tauB)/(tauA+tauB) + (tauA*tauB)**2/(tauA+tauB)**3)
      end if
    else if (abs(Ua - Ub) < MinHubDiff) then
      ! R > 0 and same Hubbard U values
      tauMean = 0.5_dp*(tauA + tauB)
      expGamma = &
          & exp(-tauMean*rab) * (1.0_dp/rab + 0.6875_dp*tauMean &
          & + 0.1875_dp*rab*(tauMean**2) &
          & + 0.02083333333333333333_dp*(rab**2)*(tauMean**3))
    else
      ! using the sign convention in the review articles, not Joachim Elstner's thesis -- there's a
      ! typo there
      expGamma = gammaSubExprn(rab,tauA,tauB) + gammaSubExprn(rab,tauB,tauA)
    end if

  end function expGamma


  !> Determines the value of the derivative of the short range contribution to gamma with the
  !> exponential form
  pure function expGammaPrime(rab, Ua, Ub)

    !> separation of sites a and b
    real(dp), intent(in) :: rab

    !> Hubbard U for site a
    real(dp), intent(in) :: Ua

    !> Hubbard U for site b
    real(dp), intent(in) :: Ub

    !> returned contribution
    real(dp) :: expGammaPrime

    real(dp) :: tauA, tauB, tauMean

    ! on-site case with R~0
    if (rab < tolSameDist) then
      expGammaPrime = 0.0_dp
    else if (abs(Ua - Ub) < MinHubDiff) then
      ! R > 0 and same Hubbard U values
      ! 16/5 * U, see review papers
      tauMean = 3.2_dp * 0.5_dp * (Ua + Ub)
      expGammaPrime = &
          & -tauMean * exp(-tauMean*rab) * &
          & ( 1.0_dp/rab + 0.6875_dp*tauMean &
          & + 0.1875_dp*rab*(tauMean**2) &
          & + 0.02083333333333333333_dp*(rab**2)*(tauMean**3) ) &
          & + exp(-tauMean*rab) * &
          & ( -1.0_dp/rab**2 + 0.1875_dp*(tauMean**2) &
          & + 2.0_dp*0.02083333333333333333_dp*rab*(tauMean**3) )
    else
      ! 16/5 * U, see review papers
      tauA = 3.2_dp*Ua
      tauB = 3.2_dp*Ub
      ! using the sign convention in the review articles, not Joachim Elstner's thesis -- there's a
      ! typo there
      expGammaPrime = gammaSubExprnPrime(rab, tauA, tauB) + gammaSubExprnPrime(rab, tauB, tauA)
    end if

  end function expGammaPrime


  !> Determines the value of the short range contribution to gamma with the exponential form with
  !> damping.
  !> See J. Phys. Chem. A, 111, 10865 (2007).
  pure function expGammaDamped(rab, Ua, Ub, dampExp)

    !> separation of sites a and b
    real(dp), intent(in) :: rab

    !> Hubbard U for site a
    real(dp), intent(in) :: Ua

    !> Hubbard U for site b
    real(dp), intent(in) :: Ub

    !> Damping exponent
    real(dp), intent(in) :: dampExp

    !> returned contribution
    real(dp) :: expGammaDamped

    real(dp) :: rTmp

    rTmp = -1.0_dp * (0.5_dp * (Ua + Ub))**dampExp
    expGammaDamped = expGamma(rab, Ua, Ub) * exp(rTmp * rab**2)

  end function expGammaDamped


  !> Determines the value of the derivative of the short range contribution to gamma with the
  !> exponential form with damping
  pure function expGammaDampedPrime(rab, Ua, Ub, dampExp)

    !> separation of sites a and b
    real(dp), intent(in) :: rab

    !> Hubbard U for site a
    real(dp), intent(in) :: Ua

    !> Hubbard U for site b
    real(dp), intent(in) :: Ub

    !> Damping exponent
    real(dp), intent(in) :: dampExp

    !> returned contribution
    real(dp) :: expGammaDampedPrime

    real(dp) :: rTmp

    rTmp = -1.0_dp * (0.5_dp *(Ua + Ub))**dampExp
    expGammaDampedPrime = expGammaPrime(rab, Ua, Ub) * exp(rTmp * rab**2) &
        &+ 2.0_dp * expGamma(rab, Ua, Ub) * exp(rTmp * rab**2) * rab * rTmp

  end function expGammaDampedPrime


  !> Determines the value of the short range contribution to gamma using the old Ohno/Klopman form
  !> Caveat: This is too long ranged to use in a periodic calculation
  pure function OhnoKlopman(rab, Ua, Ub)

    !> separation of sites a and b
    real(dp), intent(in) :: rab

    !> Hubbard U for site a
    real(dp), intent(in) :: Ua

    !> Hubbard U for site b
    real(dp), intent(in) :: Ub

    !> contribution
    real(dp) :: OhnoKlopman

    OhnoKlopman = 1.0_dp/sqrt(rab**2 + 0.25_dp*(1.0_dp/Ua + 1.0_dp/Ub)**2)

  end function OhnoKlopman


  !> Determines the value of the derivative of the short range contribution to gamma using the old
  !> Ohno/Klopman form. Caveat: This is too long ranged to use in a periodic calculation
  pure function OhnoKlopmanPrime(rab, Ua, Ub)

    !> separation of sites a and b
    real(dp), intent(in) :: rab

    !> Hubbard U for site a
    real(dp), intent(in) :: Ua

    !> Hubbard U for site b
    real(dp), intent(in) :: Ub

    !> contribution
    real(dp) :: OhnoKlopmanPrime

    OhnoKlopmanPrime = -rab / (sqrt(rab**2 + 0.25_dp*(1.0_dp/Ua + 1.0_dp/Ub)**2)**3)

  end function OhnoKlopmanPrime


  !> Determines the value of a part of the short range contribution to the exponential gamma, when
  !> Ua /= Ub and R > 0
  pure function gammaSubExprn(rab,tau1,tau2)

    !> separation of sites a and b
    real(dp), intent(in) :: rab

    !> Charge fluctuation for site a
    real(dp), intent(in) :: tau1

    !> Charge fluctuation U for site b
    real(dp), intent(in) :: tau2

    !> contribution
    real(dp) :: gammaSubExprn

    gammaSubExprn = exp(-tau1 * rab) * &
        & ( (0.5_dp*tau2**4*tau1/(tau1**2-tau2**2)**2) - &
        & (tau2**6-3.0_dp*tau2**4*tau1**2)/(rab*(tau1**2-tau2**2)**3) )

  end function gammaSubExprn


  !> Determines the derivative of the value of a part of the short range contribution to the
  !> exponential gamma, when Ua /= Ub and R > 0
  pure function gammaSubExprnPrime(rab,tau1,tau2)

    !> separation of sites a and b
    real(dp), intent(in) :: rab

    !> Charge fluctuation for site a
    real(dp), intent(in) :: tau1

    !> Charge fluctuation U for site b
    real(dp), intent(in) :: tau2

    !> contribution
    real(dp) :: gammaSubExprnPrime

    gammaSubExprnPrime = -tau1 * exp(- tau1 * rab) * &
        &( (0.5_dp*tau2**4*tau1/(tau1**2-tau2**2)**2) - &
        & (tau2**6-3.0_dp*tau2**4*tau1**2)/(rab*(tau1**2-tau2**2)**3) ) + &
        & exp(- tau1 * rab) * (tau2**6-3.0_dp*tau2**4*tau1**2) &
        & / (rab**2 *(tau1**2-tau2**2)**3)

  end function gammaSubExprnPrime

end module dftbp_shortgamma
