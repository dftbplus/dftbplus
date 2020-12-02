!--------------------------------------------------------------------------------------------------!
!  DFTB+: general package for performing fast atomistic simulations                                !
!  Copyright (C) 2006 - 2020  DFTB+ developers group                                               !
!                                                                                                  !
!  See the LICENSE file for terms of usage and distribution.                                       !
!--------------------------------------------------------------------------------------------------!

!> Contains functions to calculate the short-ranged part of gamma, and the distance beyond which it
!> becomes negligible.
module dftbp_shortgamma
  use dftbp_accuracy
  use dftbp_message
  implicit none

  private

  public :: expGamma, expGammaDamped, expGammaPrime, expGammaDampedPrime
  public :: expGammaCutoff

contains


  !> Determines the cut off where the short range part goes to zero (a small constant)
  function expGammaCutoff(U2, U1, minValue)

    !> Hubbard U value in a.u.
    real(dp), intent(in) :: U2

    !> Hubbard U value in a.u.
    real(dp), intent(in) :: U1

    !> value below which the short range contribution is considered negligible. If not set this
    !> comes from a constant in the precision module.
    real(dp), intent(in), optional :: minValue

    !> Returned cut off distance
    real(dp) :: expGammaCutoff

    real(dp) :: cutValue
    real(dp) :: rab, MaxGamma, MinGamma, lowerGamma, gamma
    real(dp) :: cut, MaxCutOff, MinCutOff
    character(len=100) :: error_string

    expGammaCutoff = 0.0_dp

    if (present(minValue)) then
      cutValue = minValue
    else
      cutValue = minShortGamma
    end if

    if (cutValue < tolShortGamma) then
99000 format ('Failure in determining short-range cut-off,', &
          & ' -ve cutoff negative :',f12.6)
      write(error_string, 99000) cutValue
      call error(error_string)
    else if (U1 < minHubTol .or. U2 < minHubTol) then
99010 format ('Failure in short-range gamma, U too small :',f12.6,f12.6)
      write(error_string, 99010) U1, U2
      call error(error_string)
    end if

    rab = 1.0_dp
    do while(expGamma(rab,U2,U1) > cutValue)
      rab = 2.0_dp*rab
    end do
    if (rab < 2.0_dp) then
99020 format ('Failure in short-range gamma cut-off : ', &
          & 'requested tolerance too large : ',f10.6)
      write(error_string, 99020) cutValue
      call error(error_string)
    end if
    ! bisection search for the value where the contribution drops below cutValue
    MinCutOff = rab + 0.1_dp
    MaxCutOff = 0.5_dp * rab - 0.1_dp
    maxGamma = expGamma(MaxCutOff,U2,U1)
    minGamma = expGamma(MinCutOff,U2,U1)
    lowerGamma =  expGamma(MinCutOff,U2,U1)
    cut = MaxCutOff + 0.1_dp
    gamma =  expGamma(cut,U2,U1)
    do While ((gamma-lowerGamma) > tolShortGamma)
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
    expGammaCutoff = MinCutOff

  end function expGammaCutoff


  !> Determines the value of the short range contribution to gamma with the exponential form
  function expGamma(rab,Ua,Ub)

    !> separation of sites a and b
    real(dp), intent(in) :: rab

    !> Hubbard U for site a
    real(dp), intent(in) :: Ua

    !> Hubbard U for site b
    real(dp), intent(in) :: Ub

    !> contribution
    real(dp) :: expGamma

    real(dp) :: tauA, tauB, tauMean
    character(len=100) :: error_string

    if (rab < 0.0_dp) then
99030 format ('Failure in short-range gamma, r_ab negative :',f12.6)
      write(error_string, 99030) rab
      call error(error_string)
    else if (Ua < MinHubTol) then
99040 format ('Failure in short-range gamma, U too small :',f12.6)
      write(error_string, 99040) Ua
      call error(error_string)
    else if (Ub < MinHubTol) then
99050 format ('Failure in short-range gamma, U too small : ',f12.6)
      write(error_string, 99050) Ub
      call error(error_string)
    end if

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
  function expGammaPrime(rab,Ua,Ub)

    !> separation of sites a and b
    real(dp), intent(in) :: rab

    !> Hubbard U for site a
    real(dp), intent(in) :: Ua

    !> Hubbard U for site b
    real(dp), intent(in) :: Ub

    !> returned contribution
    real(dp) :: expGammaPrime

    real(dp) :: tauA, tauB, tauMean
    character(len=100) :: error_string

    if (rab < 0.0_dp) then
99060 format ('Failure in short-range gamma, r_ab negative :',f12.6)
      write(error_string, 99060) rab
      call error(error_string)
    else if (Ua < MinHubTol) then
99070 format ('Failure in short-range gamma, U too small :',f12.6)
      write(error_string, 99070) Ua
      call error(error_string)
    else if (Ub < MinHubTol) then
99080 format ('Failure in short-range gamma, U too small : ',f12.6)
      write(error_string, 99080) Ub
      call error(error_string)
    end if

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
      expGammaPrime = gammaSubExprnPrime(rab,tauA,tauB) + &
          & gammaSubExprnPrime(rab,tauB,tauA)
    end if
  end function expGammaPrime


  !> Determines the value of the short range contribution to gamma with the exponential form with
  !> damping.
  !> See J. Phys. Chem. A, 111, 10865 (2007).
  function expGammaDamped(rab, Ua, Ub, dampExp)

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
  function expGammaDampedPrime(rab, Ua, Ub, dampExp)

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
  function OhnoKlopman(rab,Ua,Ub)

    !> separation of sites a and b
    real(dp), intent(in) :: rab

    !> Hubbard U for site a
    real(dp), intent(in) :: Ua

    !> Hubbard U for site b
    real(dp), intent(in) :: Ub

    !> contribution
    real(dp) :: OhnoKlopman

    character(len=100) :: error_string

    if (Ua < MinHubTol) then
99090 format ('Failure in short-range gamma, U too small : ',f12.6)
      write(error_string, 99090) Ua
      call error(error_string)
    else if (Ub < MinHubTol) then
99100 format ('Failure in short-range gamma, U too small : ',f12.6)
      write(error_string, 99100) Ub
      call error(error_string)
    end if

    OhnoKlopman = 1.0_dp/sqrt(rab**2 + 0.25_dp*(1.0_dp/Ua + 1.0_dp/Ub)**2)

  end function OhnoKlopman


  !> Determines the value of the derivative of the short range contribution to gamma using the old
  !> Ohno/Klopman form. Caveat: This is too long ranged to use in a periodic calculation
  function OhnoKlopmanPrime(rab,Ua,Ub)

    !> separation of sites a and b
    real(dp), intent(in) :: rab

    !> Hubbard U for site a
    real(dp), intent(in) :: Ua

    !> Hubbard U for site b
    real(dp), intent(in) :: Ub

    !> contribution
    real(dp) :: OhnoKlopmanPrime

    character(len=100) :: error_string

    if (Ua < MinHubTol) then
99110 format ('Failure in short-range gamma, U too small : ',f12.6)
      write(error_string, 99110) Ua
      call error(error_string)
    else if (Ub < MinHubTol) then
99120 format ('Failure in short-range gamma, U too small : ',f12.6)
      write(error_string, 99120) Ub
      call error(error_string)
    end if

    OhnoKlopmanPrime = -rab / &
        & (sqrt(rab**2 + 0.25_dp*(1.0_dp/Ua + 1.0_dp/Ub)**2)**3)

  end function OhnoKlopmanPrime


  !> Determines the value of a part of the short range contribution to the exponential gamma, when
  !> Ua /= Ub and R > 0
  function gammaSubExprn(rab,tau1,tau2)

    !> separation of sites a and b
    real(dp), intent(in) :: rab

    !> Charge fluctuation for site a
    real(dp), intent(in) :: tau1

    !> Charge fluctuation U for site b
    real(dp), intent(in) :: tau2

    !> contribution
    real(dp) :: gammaSubExprn

    character(len=100) :: error_string

    if (abs(tau1 - tau2) < 3.2_dp*MinHubDiff) then
99130 format ('Failure in gammaSubExprn, both tau degenerate ',f12.6,f12.6)
      write(error_string, 99130) tau1,tau2
      call error(error_string)
    else if (rab < tolSameDist) then
99140 format ('Atoms on top of each other in gammaSubExprn')
      write(error_string, 99140)
      call error(error_string)
    end if

    gammaSubExprn = exp(-tau1 * rab) * &
        & ( (0.5_dp*tau2**4*tau1/(tau1**2-tau2**2)**2) - &
        & (tau2**6-3.0_dp*tau2**4*tau1**2)/(rab*(tau1**2-tau2**2)**3) )

  end function gammaSubExprn


  !> Determines the derivative of the value of a part of the short range contribution to the
  !> exponential gamma, when Ua /= Ub and R > 0
  function gammaSubExprnPrime(rab,tau1,tau2)

    !> separation of sites a and b
    real(dp), intent(in) :: rab

    !> Charge fluctuation for site a
    real(dp), intent(in) :: tau1

    !> Charge fluctuation U for site b
    real(dp), intent(in) :: tau2

    !> contribution
    real(dp) :: gammaSubExprnPrime

    character(len=100) :: error_string

    if (abs(tau1 - tau2) < 3.2_dp*MinHubDiff) then
99150 format ('Failure in gammaSubExprn, both tau degenerate ',f12.6,f12.6)
      write(error_string, 99150) tau1,tau2
      call error(error_string)
    else if (rab < tolSameDist) then
99160 format ('Atoms on top of each other in gammaSubExprn')
      write(error_string, 99160)
      call error(error_string)
    end if

    gammaSubExprnPrime = -tau1 * exp(- tau1 * rab) * &
        &( (0.5_dp*tau2**4*tau1/(tau1**2-tau2**2)**2) - &
        & (tau2**6-3.0_dp*tau2**4*tau1**2)/(rab*(tau1**2-tau2**2)**3) ) + &
        & exp(- tau1 * rab) * (tau2**6-3.0_dp*tau2**4*tau1**2) &
        & / (rab**2 *(tau1**2-tau2**2)**3)

  end function gammaSubExprnPrime

end module dftbp_shortgamma
