!--------------------------------------------------------------------------------------------------!
!  DFTB+: general package for performing fast atomistic simulations                                !
!  Copyright (C) 2006 - 2023  DFTB+ developers group                                               !
!                                                                                                  !
!  See the LICENSE file for terms of usage and distribution.                                       !
!--------------------------------------------------------------------------------------------------!

!> Contains mathematical functions to calculate the short-ranged part of gamma, and the distance
!> beyond which it becomes negligible.
module dftbp_dftb_shortgammafuncs
  use dftbp_common_accuracy, only : dp, lc, minHubTol, minHubDiff, tolSameDist, minShortGamma,&
      & tolShortGamma
  use dftbp_io_message, only : error
  implicit none

  private
  public :: expGammaCutoff, expGamma, expGammaPrime, expGammaDamped, expGammaDampedPrime
  public :: expGammaDoublePrime, expGammaTriplePrime
  public :: expGammaQuadruplePrime, expGammaQuintuplePrime

  !> Error handling string
  character(lc) :: errorString

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
    real(dp) :: rab, maxGamma, MinGamma, lowerGamma, gamma
    real(dp) :: cut, maxCutOff, MinCutOff

    expGammaCutoff = 0.0_dp

    if (present(minValue)) then
      cutValue = minValue
    else
      cutValue = minShortGamma
    end if

    if (cutValue < tolShortGamma) then
      write(errorString,&
          & "('Failure in determining short-range cut-off, cutoff negative :',f12.6)")&
          & cutValue
      call error(errorString)
    else if (U1 < minHubTol .or. U2 < minHubTol) then
      write(errorString, "('Failure in short-range gamma, U too small :',f12.6,f12.6)") U1, U2
      call error(errorString)
    end if

    rab = 1.0_dp
    do while(expGamma(rab,U2,U1) > cutValue)
      rab = 2.0_dp*rab
    end do
    if (rab < 2.0_dp) then
      write(errorString,&
          & "('Failure in short-range gamma cut-off: requested tolerance too large: ',f10.6)")&
          & cutValue
      call error(errorString)
    end if
    ! bisection search for the value where the contribution drops below cutValue
    MinCutOff = rab + 0.1_dp
    maxCutOff = 0.5_dp * rab - 0.1_dp
    maxGamma = expGamma(maxCutOff,U2,U1)
    minGamma = expGamma(MinCutOff,U2,U1)
    lowerGamma =  expGamma(MinCutOff,U2,U1)
    cut = maxCutOff + 0.1_dp
    gamma =  expGamma(cut,U2,U1)
    do While ((gamma-lowerGamma) > tolShortGamma)
      maxCutOff = 0.5_dp*(cut + MinCutOff)
      if ((maxGamma >= minGamma) .eqv. (cutValue >= &
          & expGamma(maxCutOff,U2,U1))) then
        MinCutOff = maxCutOff
        lowerGamma =  expGamma(MinCutOff,U2,U1)
      else
        cut = maxCutOff
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

    if (rab < 0.0_dp) then
      write(errorString, "('Failure in short-range gamma, r_ab negative :', f12.6)") rab
      call error(errorString)
    else if (Ua < MinHubTol) then
      write(errorString, "('Failure in short-range gamma, U too small :', f12.6)") Ua
      call error(errorString)
    else if (Ub < MinHubTol) then
      write(errorString, "('Failure in short-range gamma, U too small : ', f12.6)") Ub
      call error(errorString)
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
      expGamma = gammaSubExprn_(rab,tauA,tauB) + gammaSubExprn_(rab,tauB,tauA)
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

    if (rab < 0.0_dp) then
      write(errorString, "('Failure in short-range gamma, r_ab negative :', f12.6)") rab
      call error(errorString)
    else if (Ua < MinHubTol) then
      write(errorString, "('Failure in short-range gamma, U too small :', f12.6)") Ua
      call error(errorString)
    else if (Ub < MinHubTol) then
      write(errorString, "('Failure in short-range gamma, U too small : ',f12.6)") Ub
      call error(errorString)
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
      expGammaPrime = gammaSubExprnPrime_(rab,tauA,tauB)&
          & + gammaSubExprnPrime_(rab,tauB,tauA)
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


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!  Private routines
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


  !> Determines the value of a part of the short range contribution to the exponential gamma, when
  !> Ua /= Ub and R > 0
  function gammaSubExprn_(rab, tau1, tau2) result(gammaSub)

    !> separation of sites a and b
    real(dp), intent(in) :: rab

    !> Charge fluctuation for site a
    real(dp), intent(in) :: tau1

    !> Charge fluctuation U for site b
    real(dp), intent(in) :: tau2

    !> contribution
    real(dp) :: gammaSub

    if (abs(tau1 - tau2) < 3.2_dp*MinHubDiff) then
      write(errorString, "('Failure in gammaSubExprn, both tau degenerate ', f12.6, f12.6)")&
          & tau1, tau2
      call error(errorString)
    else if (rab < tolSameDist) then
      write(errorString, "('Atoms on top of each other in gammaSubExprn')")
      call error(errorString)
    end if

    gammaSub = exp(-tau1 * rab)&
        & * ((0.5_dp*tau2**4*tau1/(tau1**2-tau2**2)**2)&
        & - (tau2**6-3.0_dp*tau2**4*tau1**2)/(rab*(tau1**2-tau2**2)**3))

  end function gammaSubExprn_


  !> Determines the derivative of the value of a part of the short range contribution to the
  !> exponential gamma, when Ua /= Ub and R > 0
  function gammaSubExprnPrime_(rab, tau1, tau2) result(gammaSubPrime)

    !> separation of sites a and b
    real(dp), intent(in) :: rab

    !> Charge fluctuation for site a
    real(dp), intent(in) :: tau1

    !> Charge fluctuation U for site b
    real(dp), intent(in) :: tau2

    !> contribution
    real(dp) :: gammaSubPrime

    if (abs(tau1 - tau2) < 3.2_dp * MinHubDiff) then
      write(errorString, "('Failure in gammaSubExprn, both tau degenerate ', f12.6, f12.6)")&
          & tau1, tau2
      call error(errorString)
    else if (rab < tolSameDist) then
      write(errorString, "('Atoms on top of each other in gammaSubExprn')")
      call error(errorString)
    end if

    gammaSubPrime = -tau1 * exp(- tau1 * rab)&
        & * ((0.5_dp*tau2**4*tau1/(tau1**2-tau2**2)**2)&
        & - (tau2**6-3.0_dp*tau2**4*tau1**2)/(rab*(tau1**2-tau2**2)**3) )&
        & + exp(- tau1 * rab) * (tau2**6-3.0_dp*tau2**4*tau1**2) &
        & / (rab**2 *(tau1**2-tau2**2)**3)

  end function gammaSubExprnPrime_


  !> Determines the value of a part of the short range contribution to the exponential gamma, when
  !> Ua /= Ub and R > 0
  function gammaSubfExprn(rab,tau1,tau2)

    !> separation of sites a and b
    real(dp), intent(in) :: rab

    !> Charge fluctuation for site a
    real(dp), intent(in) :: tau1

    !> Charge fluctuation U for site b
    real(dp), intent(in) :: tau2

    !> contribution
    real(dp) :: gammaSubfExprn

    if (abs(tau1 - tau2) < 3.2_dp*MinHubDiff) then
99130 format ('Failure in gammaSubfExprn, both tau degenerate ',f12.6,f12.6)
      write(errorString, 99130) tau1,tau2
      call error(errorString)
    else if (rab < tolSameDist) then
99140 format ('Atoms on top of each other in gammaSubfExprn')
      write(errorString, 99140)
      call error(errorString)
    end if

    gammaSubfExprn = &
        & (0.5_dp * tau2**4 * tau1 / (tau1**2 - tau2**2)**2) &
        & - ((tau2**6 - 3.0_dp * tau2**4 * tau1**2) / (rab * (tau1**2 - tau2**2)**3))

  end function gammaSubfExprn


  !> Determines the value of a part of the short range contribution to the exponential gamma, when
  !> Ua /= Ub and R > 0
  function gammaSubgExprn(rab,tauMean)

    !> separation of sites a and b
    real(dp), intent(in) :: rab

    !> Charge fluctuation for site a
    real(dp), intent(in) :: tauMean

    !> contribution
    real(dp) :: gammaSubgExprn

    if (rab < tolSameDist) then
99140 format ('Atoms on top of each other in gammaSubgExprn')
      write(errorString, 99140)
      call error(errorString)
    end if

    gammaSubgExprn = &
       & 1.0_dp / rab + 0.6875_dp * tauMean &
       & + 0.1875_dp * rab * (tauMean**2) &
       & + 0.02083333333333333333_dp * (rab**2) * (tauMean**3)

  end function gammaSubgExprn


  !> Determines the derivative of the value of a part of the short range contribution to the
  !> exponential gamma, when Ua /= Ub and R > 0
  function gammaSubfExprnPrime(rab, tau1, tau2)

    !> separation of sites a and b
    real(dp), intent(in) :: rab

    !> Charge fluctuation for site a
    real(dp), intent(in) :: tau1

    !> Charge fluctuation U for site b
    real(dp), intent(in) :: tau2

    !> contribution
    real(dp) :: gammaSubfExprnPrime

    if (abs(tau1 - tau2) < 3.2_dp*MinHubDiff) then
99150 format ('Failure in gammaSubfExprnPrime, both tau degenerate ',f12.6,f12.6)
      write(errorString, 99150) tau1,tau2
      call error(errorString)
    else if (rab < tolSameDist) then
99160 format ('Atoms on top of each other in gammaSubfExprnPrime')
      write(errorString, 99160)
      call error(errorString)
    end if

    gammaSubfExprnPrime = &
        & (tau2**6 - 3.0_dp * tau2**4 * tau1**2) &
        & / (rab**2 *(tau1**2 - tau2**2)**3)

  end function gammaSubfExprnPrime


  !> Determines the value of a part of the short range contribution to the exponential gamma, when
  !> Ua /= Ub and R > 0
  function gammaSubgExprnPrime(rab, tauMean)

    !> separation of sites a and b
    real(dp), intent(in) :: rab

    !> Charge fluctuation for site a
    real(dp), intent(in) :: tauMean

    !> contribution
    real(dp) :: gammaSubgExprnPrime

    if (rab < tolSameDist) then
99140 format ('Atoms on top of each other in gammaSubgExprnPrime')
      write(errorString, 99140)
      call error(errorString)
    end if

    gammaSubgExprnPrime = &
        & - 1.0_dp / rab**2 + 0.1875_dp * (tauMean**2) &
        & + 2.0_dp * 0.02083333333333333333_dp * rab * (tauMean**3)

  end function gammaSubgExprnPrime



  !> Determines the seconde derivative of the value of a part of the short range contribution to the
  !> exponential gamma, when Ua /= Ub and R > 0
  function gammaSubfExprnDoublePrime(rab, tau1, tau2)

    !> separation of sites a and b
    real(dp), intent(in) :: rab

    !> Charge fluctuation for site a
    real(dp), intent(in) :: tau1

    !> Charge fluctuation U for site b
    real(dp), intent(in) :: tau2

    !> contribution
    real(dp) :: gammaSubfExprnDoublePrime

    if (abs(tau1 - tau2) < 3.2_dp*MinHubDiff) then
99150 format ('Failure in gammaSubfExprnDoublePrime, both tau degenerate ',f12.6,f12.6)
      write(errorString, 99150) tau1,tau2
      call error(errorString)
    else if (rab < tolSameDist) then
99160 format ('Atoms on top of each other in gammaSubfExprnDoublePrime')
      write(errorString, 99160)
      call error(errorString)
    end if

    gammaSubfExprnDoublePrime = &
        & -2.0_dp * (tau2**6 - 3.0_dp * tau2**4 * tau1**2) &
        & / (rab**3 * (tau1**2 - tau2**2)**3)
  end function gammaSubfExprnDoublePrime

  !> Determines the value of a part of the short range contribution to the exponential gamma, when
  !> Ua /= Ub and R > 0
  function gammaSubgExprnDoublePrime(rab, tauMean)

    !> separation of sites a and b
    real(dp), intent(in) :: rab

    !> Charge fluctuation for site a
    real(dp), intent(in) :: tauMean

    !> contribution
    real(dp) :: gammaSubgExprnDoublePrime

    if (rab < tolSameDist) then
99140 format ('Atoms on top of each other in gammaSubgExprnDoublePrime')
      write(errorString, 99140)
      call error(errorString)
    end if

    gammaSubgExprnDoublePrime = &
        & 2.0_dp / rab**3 + 2.0_dp * 0.02083333333333333333_dp * (tauMean**3)

  end function gammaSubgExprnDoublePrime


  !> Determines the seconde derivative of the value of a part of the short range contribution to the
  !> exponential gamma, when Ua /= Ub and R > 0
  function gammaSubfExprnTriplePrime(rab, tau1, tau2)

    !> separation of sites a and b
    real(dp), intent(in) :: rab

    !> Charge fluctuation for site a
    real(dp), intent(in) :: tau1

    !> Charge fluctuation U for site b
    real(dp), intent(in) :: tau2

    !> contribution
    real(dp) :: gammaSubfExprnTriplePrime

    if (abs(tau1 - tau2) < 3.2_dp*MinHubDiff) then
99150 format ('Failure in gammaSubfExprnTriplePrime, both tau degenerate ',f12.6,f12.6)
      write(errorString, 99150) tau1,tau2
      call error(errorString)
    else if (rab < tolSameDist) then
99160 format ('Atoms on top of each other in gammaSubfExprnTriplePrime')
      write(errorString, 99160)
      call error(errorString)
    end if

    gammaSubfExprnTriplePrime = &
        & 6.0_dp * (tau2**6 - 3.0_dp * tau2**4 * tau1**2) &
        & / (rab**4 * (tau1**2 - tau2**2)**3)
  end function gammaSubfExprnTriplePrime

  !> Determines the value of a part of the short range contribution to the exponential gamma, when
  !> Ua /= Ub and R > 0
  function gammaSubgExprnTriplePrime(rab, tauMean)

    !> separation of sites a and b
    real(dp), intent(in) :: rab

    !> Charge fluctuation for site a
    real(dp), intent(in) :: tauMean

    !> contribution
    real(dp) :: gammaSubgExprnTriplePrime

    if (rab < tolSameDist) then
99140 format ('Atoms on top of each other in gammaSubgExprnTriplePrime')
      write(errorString, 99140)
      call error(errorString)
    end if

    gammaSubgExprnTriplePrime = -6.0_dp / rab**4

  end function gammaSubgExprnTriplePrime


  !> Determines the seconde derivative of the value of a part of the short range contribution to the
  !> exponential gamma, when Ua /= Ub and R > 0
  function gammaSubfExprnQuadruplePrime(rab, tau1, tau2)

    !> separation of sites a and b
    real(dp), intent(in) :: rab

    !> Charge fluctuation for site a
    real(dp), intent(in) :: tau1

    !> Charge fluctuation U for site b
    real(dp), intent(in) :: tau2

    !> contribution
    real(dp) :: gammaSubfExprnQuadruplePrime

    if (abs(tau1 - tau2) < 3.2_dp*MinHubDiff) then
99150 format ('Failure in gammaSubfExprnQuadruplePrime, both tau degenerate ',f12.6,f12.6)
      write(errorString, 99150) tau1,tau2
      call error(errorString)
    else if (rab < tolSameDist) then
99160 format ('Atoms on top of each other in gammaSubfExprnQuadruplePrime')
      write(errorString, 99160)
      call error(errorString)
    end if

    gammaSubfExprnQuadruplePrime = &
        & -24.0_dp * (tau2**6 - 3.0_dp * tau2**4 * tau1**2) &
        & / (rab**5 * (tau1**2 - tau2**2)**3)
  end function gammaSubfExprnQuadruplePrime

  !> Determines the value of a part of the short range contribution to the exponential gamma, when
  !> Ua /= Ub and R > 0
  function gammaSubgExprnQuadruplePrime(rab, tauMean)

    !> separation of sites a and b
    real(dp), intent(in) :: rab

    !> Charge fluctuation for site a
    real(dp), intent(in) :: tauMean

    !> contribution
    real(dp) :: gammaSubgExprnQuadruplePrime

    if (rab < tolSameDist) then
99140 format ('Atoms on top of each other in gammaSubgExprnQuadruplePrime')
      write(errorString, 99140)
      call error(errorString)
    end if

    gammaSubgExprnQuadruplePrime = 24.0_dp / rab**5

  end function gammaSubgExprnQuadruplePrime


  !> Determines the seconde derivative of the value of a part of the short range contribution to the
  !> exponential gamma, when Ua /= Ub and R > 0
  function gammaSubfExprnQuintuplePrime(rab, tau1, tau2)

    !> separation of sites a and b
    real(dp), intent(in) :: rab

    !> Charge fluctuation for site a
    real(dp), intent(in) :: tau1

    !> Charge fluctuation U for site b
    real(dp), intent(in) :: tau2

    !> contribution
    real(dp) :: gammaSubfExprnQuintuplePrime

    if (abs(tau1 - tau2) < 3.2_dp*MinHubDiff) then
99150 format ('Failure in gammaSubfExprnQuintuplePrime, both tau degenerate ',f12.6,f12.6)
      write(errorString, 99150) tau1,tau2
      call error(errorString)
    else if (rab < tolSameDist) then
99160 format ('Atoms on top of each other in gammaSubfExprnQuintuplePrime')
      write(errorString, 99160)
      call error(errorString)
    end if

    gammaSubfExprnQuintuplePrime = &
        & 120.0_dp * (tau2**6 - 3.0_dp * tau2**4 * tau1**2) &
        & / (rab**6 * (tau1**2 - tau2**2)**3)
  end function gammaSubfExprnQuintuplePrime

  !> Determines the value of a part of the short range contribution to the exponential gamma, when
  !> Ua /= Ub and R > 0
  function gammaSubgExprnQuintuplePrime(rab, tauMean)

    !> separation of sites a and b
    real(dp), intent(in) :: rab

    !> Charge fluctuation for site a
    real(dp), intent(in) :: tauMean

    !> contribution
    real(dp) :: gammaSubgExprnQuintuplePrime

    if (rab < tolSameDist) then
99140 format ('Atoms on top of each other in gammaSubgExprnQuintuplePrime')
      write(errorString, 99140)
      call error(errorString)
    end if

    gammaSubgExprnQuintuplePrime = -120.0_dp / rab**6

  end function gammaSubgExprnQuintuplePrime


  !> Determines the value of the second derivative of the short range contribution to gamma with the
  !> exponential form
  function expGammaDoublePrime(rab, Ua, Ub)

    !> separation of sites a and b
    real(dp), intent(in) :: rab

    !> Hubbard U for site a
    real(dp), intent(in) :: Ua

    !> Hubbard U for site b
    real(dp), intent(in) :: Ub

    !> returned contribution
    real(dp) :: expGammaDoublePrime

    real(dp) :: tauA, tauB, tauMean

    if (rab < 0.0_dp) then
99060 format ('Failure in short-range gamma, r_ab negative :',f12.6)
      write(errorString, 99060) rab
      call error(errorString)
    else if (Ua < MinHubTol) then
99070 format ('Failure in short-range gamma, U too small :',f12.6)
      write(errorString, 99070) Ua
      call error(errorString)
    else if (Ub < MinHubTol) then
99080 format ('Failure in short-range gamma, U too small : ',f12.6)
      write(errorString, 99080) Ub
      call error(errorString)
    end if

    ! on-site case with R~0
    if (rab < tolSameDist) then
      if (abs(Ua - Ub) < MinHubDiff) then
        ! same Hubbard U values, onsite , NOTE SIGN CHANGE!
        tauMean = 3.2_dp * 0.5_dp * (Ua + Ub)
        expGammaDoublePrime = -tauMean**3 / 48.0_dp
      else
        ! Ua /= Ub Hubbard U values - limiting case, NOTE SIGN CHANGE!
        expGammaDoublePrime = 0.0_dp
      end if
    else if (abs(Ua - Ub) < MinHubDiff) then
      ! R > 0 and same Hubbard U values
      ! 16/5 * U, see review papers
      tauMean = 3.2_dp * 0.5_dp * (Ua + Ub)
      expGammaDoublePrime = &
          & exp(-tauMean * rab) * gammaSubgExprnDoublePrime(rab, tauMean) &
          & - 2.0_dp * tauMean * exp(-tauMean * rab) * gammaSubgExprnPrime(rab, tauMean) &
          & + tauMean**2 * exp(-tauMean * rab) * gammaSubgExprn(rab, tauMean)
    else
      ! 16/5 * U, see review papers
      tauA = 3.2_dp*Ua
      tauB = 3.2_dp*Ub
      ! using the sign convention in the review articles, not Joachim Elstner's thesis -- there's a
      ! typo there
      expGammaDoublePrime = &
          & exp(-tauA * rab) * gammaSubfExprnDoublePrime(rab, tauA, tauB) &
          & - 2.0_dp * tauA * exp(-tauA * rab) * gammaSubfExprnPrime(rab, tauA, tauB) &
          & + tauA**2 * exp(-tauA * rab) * gammaSubfExprn(rab, tauA, tauB) &
          & + exp(-tauB * rab) * gammaSubfExprnDoublePrime(rab, tauB, tauA) &
          & - 2.0_dp * tauB * exp(-tauB * rab) * gammaSubfExprnPrime(rab, tauB, tauA) &
          & + tauB**2 * exp(-tauB * rab) * gammaSubfExprn(rab, tauB, tauA)
    end if
  end function expGammaDoublePrime


  !> Determines the value of the third derivative of the short range contribution to gamma with the
  !> exponential form
  function expGammaTriplePrime(rab, Ua, Ub)

    !> separation of sites a and b
    real(dp), intent(in) :: rab

    !> Hubbard U for site a
    real(dp), intent(in) :: Ua

    !> Hubbard U for site b
    real(dp), intent(in) :: Ub

    !> returned contribution
    real(dp) :: expGammaTriplePrime

    real(dp) :: tauA, tauB, tauMean

    if (rab < 0.0_dp) then
99060 format ('Failure in short-range gamma, r_ab negative :',f12.6)
      write(errorString, 99060) rab
      call error(errorString)
    else if (Ua < MinHubTol) then
99070 format ('Failure in short-range gamma, U too small :',f12.6)
      write(errorString, 99070) Ua
      call error(errorString)
    else if (Ub < MinHubTol) then
99080 format ('Failure in short-range gamma, U too small : ',f12.6)
      write(errorString, 99080) Ub
      call error(errorString)
    end if

    ! on-site case with R~0
    if (rab < tolSameDist) then
      expGammaTriplePrime = 0.0_dp
    else if (abs(Ua - Ub) < MinHubDiff) then
      ! R > 0 and same Hubbard U values
      ! 16/5 * U, see review papers
      tauMean = 3.2_dp * 0.5_dp * (Ua + Ub)
      expGammaTriplePrime = &
          & exp(-tauMean * rab) * gammaSubgExprnTriplePrime(rab, tauMean) &
          & - 3.0_dp * tauMean * exp(-tauMean * rab) * gammaSubgExprnDoublePrime(rab, tauMean) &
          & + 3.0_dp * tauMean**2 * exp(-tauMean * rab) * gammaSubgExprnPrime(rab, tauMean) &
          & - tauMean**3 * exp(-tauMean * rab) * gammaSubgExprn(rab, tauMean)
    else
      ! 16/5 * U, see review papers
      tauA = 3.2_dp*Ua
      tauB = 3.2_dp*Ub
      ! using the sign convention in the review articles, not Joachim Elstner's thesis -- there's a
      ! typo there
      expGammaTriplePrime = &
          & exp(-tauA * rab) * gammaSubfExprnTriplePrime(rab, tauA, tauB) &
          & - 3.0_dp * tauA * exp(-tauA * rab) * gammaSubfExprnDoublePrime(rab, tauA, tauB) &
          & + 3.0_dp * tauA**2 * exp(-tauA * rab) * gammaSubfExprnPrime(rab, tauA, tauB) &
          & - tauA**3 * exp(-tauA * rab) * gammaSubfExprn(rab, tauA, tauB) &
          & + exp(-tauB * rab) * gammaSubfExprnTriplePrime(rab, tauB, tauA) &
          & - 3.0_dp * tauB * exp(-tauB * rab) * gammaSubfExprnDoublePrime(rab, tauB, tauA) &
          & + 3.0_dp * tauB**2 * exp(-tauB * rab) * gammaSubfExprnPrime(rab, tauB, tauA) &
          & - tauB**3 * exp(-tauB * rab) * gammaSubfExprn(rab, tauB, tauA)
    end if
  end function expGammaTriplePrime


  !> Determines the value of the fourth derivative of the short range contribution to gamma with the
  !> exponential form
  function expGammaQuadruplePrime(rab, Ua, Ub)

    !> separation of sites a and b
    real(dp), intent(in) :: rab

    !> Hubbard U for site a
    real(dp), intent(in) :: Ua

    !> Hubbard U for site b
    real(dp), intent(in) :: Ub

    !> returned contribution
    real(dp) :: expGammaQuadruplePrime

    real(dp) :: tauA, tauB, tauMean

    if (rab < 0.0_dp) then
99060 format ('Failure in short-range gamma, r_ab negative :',f12.6)
      write(errorString, 99060) rab
      call error(errorString)
    else if (Ua < MinHubTol) then
99070 format ('Failure in short-range gamma, U too small :',f12.6)
      write(errorString, 99070) Ua
      call error(errorString)
    else if (Ub < MinHubTol) then
99080 format ('Failure in short-range gamma, U too small : ',f12.6)
      write(errorString, 99080) Ub
      call error(errorString)
    end if

    ! on-site case with R~0
    if (rab < tolSameDist) then
      if (abs(Ua - Ub) < MinHubDiff) then
        ! same Hubbard U values, onsite , NOTE SIGN CHANGE!
        tauMean = 3.2_dp * 0.5_dp * (Ua + Ub)
        expGammaQuadruplePrime = tauMean**5 / 80.0_dp
      else
        ! Ua /= Ub Hubbard U values - limiting case, NOTE SIGN CHANGE!
        expGammaQuadruplePrime = 0.0_dp
      end if
    else if (abs(Ua - Ub) < MinHubDiff) then
      ! R > 0 and same Hubbard U values
      ! 16/5 * U, see review papers
      tauMean = 3.2_dp * 0.5_dp * (Ua + Ub)
      expGammaQuadruplePrime = &
          & exp(-tauMean * rab) * gammaSubgExprnQuadruplePrime(rab, tauMean) &
          & - 4.0_dp * tauMean * exp(-tauMean * rab) * gammaSubgExprnTriplePrime(rab, tauMean) &
          & + 6.0_dp * tauMean**2 * exp(-tauMean * rab) * gammaSubgExprnDoublePrime(rab, tauMean) &
          & - 4.0_dp * tauMean**3 * exp(-tauMean * rab) * gammaSubgExprnPrime(rab, tauMean) &
          & + tauMean**4 * exp(-tauMean * rab) * gammaSubgExprn(rab, tauMean)
    else
      ! 16/5 * U, see review papers
      tauA = 3.2_dp*Ua
      tauB = 3.2_dp*Ub
      ! using the sign convention in the review articles, not Joachim Elstner's thesis -- there's a
      ! typo there
      expGammaQuadruplePrime = &
          & exp(-tauA * rab) * gammaSubfExprnQuadruplePrime(rab, tauA, tauB) &
          & - 4.0_dp * tauA * exp(-tauA * rab) * gammaSubfExprnTriplePrime(rab, tauA, tauB) &
          & + 6.0_dp * tauA**2 * exp(-tauA * rab) * gammaSubfExprnDoublePrime(rab, tauA, tauB) &
          & - 4.0_dp * tauA**3 * exp(-tauA * rab) * gammaSubfExprnPrime(rab, tauA, tauB) &
          & + tauA**4 * exp(-tauA * rab) * gammaSubfExprn(rab, tauA, tauB) &
          & + exp(-tauB * rab) * gammaSubfExprnQuadruplePrime(rab, tauB, tauA) &
          & - 4.0_dp * tauB * exp(-tauB * rab) * gammaSubfExprnTriplePrime(rab, tauB, tauA) &
          & + 6.0_dp * tauB**2 * exp(-tauB * rab) * gammaSubfExprnDoublePrime(rab, tauB, tauA) &
          & - 4.0_dp * tauB**3 * exp(-tauB * rab) * gammaSubfExprnPrime(rab, tauB, tauA) &
          & + tauB**4 * exp(-tauB * rab) * gammaSubfExprn(rab, tauB, tauA)
    end if
  end function expGammaQuadruplePrime


  !> Determines the value of the fifth derivative of the short range contribution to gamma with the
  !> exponential form
  function expGammaQuintuplePrime(rab, Ua, Ub)

    !> separation of sites a and b
    real(dp), intent(in) :: rab

    !> Hubbard U for site a
    real(dp), intent(in) :: Ua

    !> Hubbard U for site b
    real(dp), intent(in) :: Ub

    !> returned contribution
    real(dp) :: expGammaQuintuplePrime

    real(dp) :: tauA, tauB, tauMean

    if (rab < 0.0_dp) then
99060 format ('Failure in short-range gamma, r_ab negative :',f12.6)
      write(errorString, 99060) rab
      call error(errorString)
    else if (Ua < MinHubTol) then
99070 format ('Failure in short-range gamma, U too small :',f12.6)
      write(errorString, 99070) Ua
      call error(errorString)
    else if (Ub < MinHubTol) then
99080 format ('Failure in short-range gamma, U too small : ',f12.6)
      write(errorString, 99080) Ub
      call error(errorString)
    end if

    ! on-site case with R~0
    if (rab < tolSameDist) then
      expGammaQuintuplePrime = 0.0_dp
    else if (abs(Ua - Ub) < MinHubDiff) then
      ! R > 0 and same Hubbard U values
      ! 16/5 * U, see review papers
      tauMean = 3.2_dp * 0.5_dp * (Ua + Ub)
      expGammaQuintuplePrime = &
          & exp(-tauMean * rab) * gammaSubgExprnQuintuplePrime(rab, tauMean) &
          & - 5.0_dp * tauMean * exp(-tauMean * rab) * gammaSubgExprnQuadruplePrime(rab, tauMean) &
          & + 10.0_dp * tauMean**2 * exp(-tauMean * rab) * gammaSubgExprnTriplePrime(rab, tauMean) &
          & - 10.0_dp * tauMean**3 * exp(-tauMean * rab) * gammaSubgExprnDoublePrime(rab, tauMean) &
          & + 5.0_dp * tauMean**4 * exp(-tauMean * rab) * gammaSubgExprnPrime(rab, tauMean) &
          & - tauMean**5 * exp(-tauMean * rab) * gammaSubgExprn(rab, tauMean)
    else
      ! 16/5 * U, see review papers
      tauA = 3.2_dp*Ua
      tauB = 3.2_dp*Ub
      ! using the sign convention in the review articles, not Joachim Elstner's thesis -- there's a
      ! typo there
      expGammaQuintuplePrime = &
          & exp(-tauA * rab) * gammaSubfExprnQuintuplePrime(rab, tauA, tauB) &
          & - 5.0_dp * tauA * exp(-tauA * rab) * gammaSubfExprnQuadruplePrime(rab, tauA, tauB) &
          & + 10.0_dp * tauA**2 * exp(-tauA * rab) * gammaSubfExprnTriplePrime(rab, tauA, tauB) &
          & - 10.0_dp * tauA**3 * exp(-tauA * rab) * gammaSubfExprnDoublePrime(rab, tauA, tauB) &
          & + 5.0_dp * tauA**4 * exp(-tauA * rab) * gammaSubfExprnPrime(rab, tauA, tauB) &
          & - tauA**5 * exp(-tauA * rab) * gammaSubfExprn(rab, tauA, tauB) &
          & + exp(-tauB * rab) * gammaSubfExprnQuintuplePrime(rab, tauB, tauA) &
          & - 5.0_dp * tauB * exp(-tauB * rab) * gammaSubfExprnQuadruplePrime(rab, tauB, tauA) &
          & + 10.0_dp * tauB**2 * exp(-tauB * rab) * gammaSubfExprnTriplePrime(rab, tauB, tauA) &
          & - 10.0_dp * tauB**3 * exp(-tauB * rab) * gammaSubfExprnDoublePrime(rab, tauB, tauA) &
          & + 5.0_dp * tauB**4 * exp(-tauB * rab) * gammaSubfExprnPrime(rab, tauB, tauA) &
          & - tauB**5 * exp(-tauB * rab) * gammaSubfExprn(rab, tauB, tauA)
    end if
  end function expGammaQuintuplePrime

end module dftbp_dftb_shortgammafuncs
