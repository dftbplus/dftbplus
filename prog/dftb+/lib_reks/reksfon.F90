!--------------------------------------------------------------------------------------------------!
!  DFTB+: general package for performing fast atomistic simulations                                !
!  Copyright (C) 2006 - 2020  DFTB+ developers group                                               !
!                                                                                                  !
!  See the LICENSE file for terms of usage and distribution.                                       !
!--------------------------------------------------------------------------------------------------!

!> REKS and SI-SA-REKS formulation in DFTB as developed by Lee et al.
!>
!> The functionality of the module has some limitation:
!> * Third order does not work.
!> * Periodic system do not work yet appart from Gamma point.
!> * Orbital potentials or spin-orbit or external E-field does not work yet.
!> * Only for closed shell system.
!> * Onsite corrections are not included in this version
module dftbp_reksfon

  use dftbp_accuracy
  use dftbp_globalenv
  use dftbp_message
  use dftbp_reksvar, only : TReksCalc, reksTypes

  implicit none

  private

  public :: optimizeFons

  !> Parameter to distinguish between two asymptotic regimes
  !> in the behavior of the derivative of the REKS coefficient, fx
  real(dp) :: Threshold = 5.0_dp

  !> Convergence for NR solver
  real(dp) :: ConvergeLimit = 1.0E-10_dp
!  real(dp) :: ConvergeLimit = 1.0E-8_dp

  contains

  !> Optimize the fractional occupation numbers (FONs) in REKS
  subroutine optimizeFons(this)

    !> data type for REKS
    type(TReksCalc), intent(inout) :: this

    real(dp) :: x

    select case (this%reksAlg)
    case (reksTypes%noReks)
    case (reksTypes%ssr22)

      call getFONs22_(x, this%hess, this%enLtot, this%delta, this%FonMaxIter, this%Plevel)
      ! FONs(1,1) = n_a, FONs(2,1) = n_b
      this%FONs(1,1) = 2.0_dp * x
      this%FONs(2,1) = 2.0_dp - this%FONs(1,1)

    case (reksTypes%ssr44)

      call error("SSR(4,4) is not implemented yet")

    end select

  end subroutine optimizeFons


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!! Private routines
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Optimize FONs in REKS(2,2) case with Newton-Raphson method
  subroutine getFONs22_(x, hess0, enLtot, delta, maxIter, opt)

    !> converged x (= n_a/2)
    real(dp), intent(out) :: x

    !> converged Hessian of FONs
    real(dp), intent(out) :: hess0

    !> total energy for each microstate
    real(dp), intent(in) :: enLtot(:)

    !> Smoothing factor used in FON optimization
    real(dp), intent(in) :: delta

    !> Maximum iteration used in FON optimization
    integer, intent(in) :: maxIter

    !> Print level in standard output file
    integer, intent(in) :: opt

    real(dp) :: fac, fac1, fac2, fac3, fac4, fac5, fac6
    real(dp) :: Const, ConDeno
    real(dp) :: xUpLim, xDownLim
    real(dp) :: root, x0, x1, y, grad, hess
    real(dp) :: eps
    integer :: iter

    ! Calculate Const in equation 12.c
    ! Reference : JCP, 147, 034113 (2017) and its supporting information 
    ConDeno = enLtot(5) - enLtot(3)
    ! ConDeno should be negative
    ! In general, E2 - E1 > 0 due to the MO swap (always, n_a > n_b)
    ! So, negative ConDeno generates negative Const
    if (ConDeno > 0.0_dp) ConDeno = -ConDeno
    Const = (enLtot(2) - enLtot(1)) / ConDeno

    ! Set up the starting value for x, x0
    !
    ! if |Const| > threshold solve the equation
    ! abs[1-2*x] = (2*(2 + d)/(1 + d)) * (4*x*(1 - x))**(-d/(2*(1 + d)))
    ! In the equation abs[1-2*x] is actually abs[Const]
    !
    ! else solve the equation:
    ! Const = (-2)*(2*x - 1)*(4*x*(1 - x))**(-1/2)
    !
    if (abs(Const) > Threshold) then
      fac = (2.0_dp*(2.0_dp+delta)/((1.0_dp+delta)*abs(Const))) &
         & **(2.0_dp*(1.0_dp+delta)/delta)
      if (Const > 0.0_dp) then
        ! x ~= 0 case
        root = 0.5_dp * (1.0_dp - sqrt(1.0_dp - fac))
      else
        ! x ~= 1 case
        root = 0.5_dp * (1.0_dp + sqrt(1.0_dp - fac))
      end if
    else
      ! x ~= near 1/2 case
      root = 0.5_dp * (1.0_dp - Const/sqrt(4.0_dp+Const**2))
    end if

    xDownLim = 1.0E-14_dp
    xUpLim = 1.0_dp - xDownLim
    x0 = root
    eps = 0.0_dp

    NRsolver: do iter = 1, maxIter

      ! Update x1 value
      x1 = x0 + eps
      ! Move to the inside of limit ( 0 < x < 1 )
      if(x1 > xUpLim)   x1 = xUpLim
      if(x1 < xDownLim) x1 = xDownLim
      ! Decide y = 4*x*(1-x) where x = n_r/2
      y = 4.0_dp*x1*(1.0_dp-x1)
      ! Calculate gradient of f(x)
      fac1 = 2.0_dp * (1.0_dp - 2.0_dp*x1) / (1.0_dp + delta)
      fac2 = -0.5_dp * (y + delta) / (1.0_dp + delta)
      fac3 = 2.0_dp + delta - y - y*log(y)
      grad = fac1 * fac3 * y**fac2
      ! Calculate hessian of f(x)
      fac1 = 4.0_dp / (1.0_dp + delta)**2
      fac2 = -1.0_dp * (2.0_dp + 3.0_dp*delta + y) / (2.0_dp + 2.0_dp*delta)
      fac4 = -1.0_dp * delta**2 - y * ( 8.0_dp + (-8.0_dp+y)*y )
      fac5 = delta * ( -2.0_dp + 5.0_dp*(-1.0_dp+y)*y )
      fac6 = 4.0_dp + delta*(2.0_dp-3.0_dp*y) + y*(-7.0_dp+2.0_dp*y) + (-1.0_dp+y)*y*log(y)
      fac3 = fac4 + fac5 - y*log(y)*fac6
      hess = fac1 * fac3 * y**fac2
      ! Update eps value
      eps = (Const - grad) / hess
      if (opt >= 2) then
        write(stdOut,'(2x,a,1x,i4,4x,a,F18.14,1x,a,F18.14)') &
       & 'NR solver: Iteration', iter, 'X =', x1, 'Eps =', eps
      end if

      ! Convergence check
      if (abs(eps) > ConvergeLimit) then
        x0 = x1
        if (iter == maxIter) then
          write(stdOut,'(2x,a,i4,a)') &
         & 'Warning! Maximum number of iterations (', maxIter, &
         & ') is exceeded in NR solver'
        end if
      else
        if (opt >= 2) then
          write(stdOut,'(2x,a,1x,i4,1x,a)') &
         & 'Convergence reached in NR solver after', iter, 'iterations'
        end if
        exit NRsolver
      end if

    end do NRsolver

    ! Converged x and hessian value
    x = x1
    hess0 = hess

  end subroutine getFONs22_


end module dftbp_reksfon
