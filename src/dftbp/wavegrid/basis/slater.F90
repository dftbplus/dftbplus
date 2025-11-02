!--------------------------------------------------------------------------------------------------!
!  DFTB+: general package for performing fast atomistic simulations                                !
!  Copyright (C) 2006 - 2025  DFTB+ developers group                                               !
!                                                                                                  !
!  See the LICENSE file for terms of usage and distribution.                                       !
!--------------------------------------------------------------------------------------------------!

#:include 'common.fypp'
#:set pure = "" if defined('WITH_ASSERT') else "pure"
!> Holds TSlaterOrbital, a concrete implementation of TOrbital.
module dftbp_wavegrid_basis_slater
  use dftbp_common_accuracy, only : dp
  use dftbp_io_message, only : error
  use dftbp_wavegrid_basis_orbital, only : TOrbital
  implicit none

  private

  public :: TSlaterOrbital, TSlaterOrbital_init
  
  !> Represents a Slater-Type Orbital (STO):
  !! R(r) = r^l * Sum_i Sum_j [ aa(i,j) * r^(i-1) * exp(-alpha(j) * r) ]
  type, extends(TOrbital) :: TSlaterOrbital
    !> Maximum power of the radial distance
    integer :: nPow

    !> Summation coefficients. Shape: [nPow, nAlpha]
    real(dp), allocatable :: aa(:,:)

    !> Exponential coefficients
    real(dp), allocatable :: alpha(:)
  contains
    procedure :: getRadial => TSlaterOrbital_getRadial
    procedure, pass(lhs) :: assign => TSlaterOrbital_assign
  end type TSlaterOrbital

contains

  !> Initialises using STO parameters.
  subroutine TSlaterOrbital_init(this, aa, alpha, angMom, cutoff)

    !> TSlaterOrbital instance to initialise
    type(TSlaterOrbital), intent(out) :: this

    !> Summation coefficients. Shape: [nCoeffPerAlpha, nAlpha]
    real(dp), intent(in) :: aa(:,:)

    !> Exponential coefficients
    real(dp), intent(in) :: alpha(:)

    !> Angular momentum of the orbital (l)
    integer, intent(in) :: angMom

    !> Cutoff, after which orbital is assumed to be zero
    real(dp), intent(in) :: cutoff

    @:ASSERT(cutoff > 0.0_dp)

    this%angMom = angMom
    this%cutoffSq = cutoff ** 2

    @:ASSERT(size(aa, dim=2) == size(alpha))
    this%nPow = size(aa, dim=1)

    allocate(this%aa, source=aa)
    allocate(this%alpha, source=alpha)

  end subroutine TSlaterOrbital_init


  !> Calculates the value of an STO analytically.
  ${pure}$ function TSlaterOrbital_getRadial(this, r) result(sto)

    !> SlaterOrbital instance
    class(TSlaterOrbital), intent(in) :: this

    !> Distance, where the STO should be calculated
    real(dp), intent(in) :: r

    !> Value of the STO on return
    real(dp) :: sto

    real(dp) :: pows(this%nPow)
    real(dp) :: rTmp
    integer :: ii, jj

    ! Avoid 0.0**0 as it may lead to arithmetic exception on pre-2008 compilers
    if (this%angMom == 0 .and. r < epsilon(1.0_dp)) then
      rTmp = 1.0_dp
    else
      rTmp = r**this%angMom
    end if

    do ii = 1, this%nPow
      pows(ii) = rTmp
      rTmp = rTmp * r
    end do
    sto = 0.0_dp
    do ii = 1, size(this%alpha)
      rTmp = 0.0_dp
      do jj = 1, this%nPow
        rTmp = rTmp + this%aa(jj, ii) * pows(jj)
      end do
      sto = sto + rTmp * exp(-this%alpha(ii) * r)
    end do

  end function TSlaterOrbital_getRadial
  
  !> Assignment operator
  subroutine TSlaterOrbital_assign(lhs, rhs)
    class(TSlaterOrbital), intent(out) :: lhs
    class(TOrbital), intent(in) :: rhs


    select type (rhs)
      type is (TSlaterOrbital)
        lhs%angMom = rhs%angMom
        lhs%cutoffSq = rhs%cutoffSq
        lhs%nPow = rhs%nPow

        if (allocated(lhs%alpha)) deallocate(lhs%alpha)
        if (allocated(lhs%aa)) deallocate(lhs%aa)

        allocate(lhs%alpha, source=rhs%alpha)
        allocate(lhs%aa, source=rhs%aa)
      class default
        call error("Cannot assign non-Slater orbital to Slater orbital.")
    end select
  end subroutine TSlaterOrbital_assign

end module dftbp_wavegrid_basis_slater
