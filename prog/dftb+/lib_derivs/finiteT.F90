!--------------------------------------------------------------------------------------------------!
!  DFTB+: general package for performing fast atomistic simulations                                !
!  Copyright (C) 2018  DFTB+ developers group                                                      !
!                                                                                                  !
!  See the LICENSE file for terms of usage and distribution.                                       !
!--------------------------------------------------------------------------------------------------!

#:include 'common.fypp'

!> Module containing finite temperature function distributions for perturbation-like expressions and
!> similar
module dftbp_finiteTHelper
  use dftbp_accuracy
  use dftbp_errorfunction
  implicit none

  public

contains
  
  !> Fermi function at En for level at Em
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
  pure function invDiff(En, Em, Ef, sigma) ! En == filled,  Em == empty

    !> Empty level energy
    real(dp), intent(in) :: En

    !> Filled level energy
    real(dp), intent(in) :: Em

    !> Fermi energy
    real(dp), intent(in) :: Ef

    !> Width of distribution
    real(dp), intent(in) :: sigma

    !> resulting value
    real(dp) :: invDiff
    
    if ( abs(En - Em) > 10.0_dp * epsilon(1.0_dp) ) then
      invDiff = ( theta(En,Ef,sigma) - theta(Em,Ef,sigma) ) / (En - Em)
    else
      invDiff = -deltamn(En,Ef,sigma)
    end if

  end function invDiff


  !> Limit function of perturbation expression in sum over states
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
  
end module dftbp_finiteTHelper
