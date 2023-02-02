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

    character(len=100) :: errorString

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
    character(len=100) :: errorString

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
    character(lc) :: errorString

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

    character(lc) :: errorString

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

    character(lc) :: errorString

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


end module dftbp_dftb_shortgammafuncs
