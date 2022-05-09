!--------------------------------------------------------------------------------------------------!
!  DFTB+: general package for performing fast atomistic simulations                                !
!  Copyright (C) 2006 - 2022  DFTB+ developers group                                               !
!                                                                                                  !
!  See the LICENSE file for terms of usage and distribution.                                       !
!--------------------------------------------------------------------------------------------------!

#:include 'common.fypp'

!> Module containing finite temperature function distributions and derivatives for the Fermi-Dirac
!> distribution
module dftbp_derivs_fermihelper
  use dftbp_common_accuracy, only : dp, mExpArg
  implicit none

  private
  public :: theta, invDiff, deltamn, dEfda, dEida

  interface invDiff
    module procedure :: invDiffStatic
    module procedure :: invDiffDynamic
  end interface invDiff

contains

  !> Fermi function at energy En for a level at energy Em
  pure function theta(En, Em, sigma)

    !> Virtual energy
    real(dp), intent(in) :: En

    !> Occupied energy
    real(dp), intent(in) :: Em

    !> distribution width
    real(dp), intent(in) :: sigma

    !> Resulting occupation value
    real(dp) :: theta

    real(dp) :: x

    x = ( En - Em ) / sigma
  #:if EXP_TRAP
    ! Where the compiler does not handle inf gracefully, trap the exponential function for small
    ! values
    if (x > mExpArg) then
      theta = 0.0_dp
    else
      theta = 1.0_dp / (1.0_dp + exp(x))
    endif
  #:else
    theta = 1.0_dp / (1.0_dp + exp(x))
  #:endif

  end function theta


  !> Expression limit aware evaluation of (Theta(occ) - Theta(empty))/(Eempty-Efilled)
  pure function invDiffStatic(En, Em, Ef, sigma)

    !> Empty level energy
    real(dp), intent(in) :: En

    !> Filled level energy
    real(dp), intent(in) :: Em

    !> Fermi energy
    real(dp), intent(in) :: Ef

    !> Width of distribution
    real(dp), intent(in) :: sigma

    !> resulting value
    real(dp) :: invDiffStatic

    if ( abs(En - Em) > 10.0_dp * epsilon(1.0_dp) ) then
      invDiffStatic = ( theta(En,Ef,sigma) - theta(Em,Ef,sigma) ) / (En - Em)
    else
      invDiffStatic = -deltamn(En,Ef,sigma)
    end if

  end function invDiffStatic


  !> Expression limit aware evaluation of (Theta(occ) - Theta(empty))/(Eempty-Efilled + omega)
  pure function invDiffDynamic(En, Em, Ef, sigma, omega)

    !> Empty level energy
    real(dp), intent(in) :: En

    !> Filled level energy
    real(dp), intent(in) :: Em

    !> Fermi energy
    real(dp), intent(in) :: Ef

    !> Width of distribution
    real(dp), intent(in) :: sigma

    !> Complex frequency
    complex(dp), intent(in) :: omega

    !> resulting value
    complex(dp) :: invDiffDynamic

    invDiffDynamic = ( theta(En,Ef,sigma) - theta(Em,Ef,sigma) ) / (En - Em + omega)

  end function invDiffDynamic


  !> Limit of function occurring in perturbation expression sum over states
  pure function deltamn(En,Em,sigma)

    !> Energy at which to measure the smeared level
    real(dp), intent(in) :: En

    !> Level centre
    real(dp), intent(in) :: Em

    !> Distribution width
    real(dp), intent(in) :: sigma

    !> Resulting contribution
    real(dp) :: deltamn

    real(dp) :: invSigma

    invSigma = 1.0_dp / sigma

    deltamn = deltatilde( (En-Em)*invSigma ) * invSigma

  end function deltamn


  !> Broadening function underlying the Fermi distribution
  pure function deltaTilde(x) result(delta)

    !> Distance from centre of level
    real(dp), intent(in) :: x

    !> Broadened value
    real(dp) :: delta

    delta  = 1.0_dp / (2.0_dp * (1.0_dp + cosh(x)))

  end function deltaTilde


  !> Evaluates the derivative of the Fermi energy given derivatives of eigenvalues with respect to a
  !> perturbation
  pure function dEfda(f, dei)

    !> Occupations
    real(dp), intent(in) :: f(:)

    !> Derivatives of eigenvalues
    real(dp), intent(in) :: dei(:)

    !> Derivative of Fermi energy
    real(dp) :: dEfda

    dEfda = sum(f * (1.0_dp-f) * dei) / sum(f * (1.0_dp-f))

  end function dEfda


  !> Evaluates the derivative of level fillings, given the derivative of eigenvalues
  pure subroutine dEida(dfill, f, dEi, kT)

    !> Derivatives of fillings
    real(dp), intent(out) :: dfill(:)

    !> Occupations
    real(dp), intent(in) :: f(:)

    !> Derivatives of eigenvalues
    real(dp), intent(in) :: dEi(:)

    !> Electronic temperature
    real(dp), intent(in) :: kT

    integer :: ii
    real(dp) :: dEf

    dEf = dEfda(f, dEi)

    do ii = 1, size(f)

      dfill(ii) = -(dEi(ii) - deF) * f(ii) * (1.0_dp - f(ii)) / kT

    end do

  end subroutine dEida

end module dftbp_derivs_fermihelper
