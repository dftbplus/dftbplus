!--------------------------------------------------------------------------------------------------!
!  DFTB+: general package for performing fast atomistic simulations                                !
!  Copyright (C) 2006 - 2023  DFTB+ developers group                                               !
!                                                                                                  !
!  See the LICENSE file for terms of usage and distribution.                                       !
!--------------------------------------------------------------------------------------------------!

#:include 'common.fypp'

!> Module for determining filled and empty bands
module dftbp_derivs_fillings
  use dftbp_common_accuracy, only : dp, rsp
  implicit none

  public

contains

  !> Counts number of decoupled matrix problems
  pure function nIndependentHam(nSpin)

    !> Number of spin channels (1=no-spin, 2=colinear, 4=non-collinear)
    integer, intent(in) :: nSpin

    !> Resulting number of separate matrix problems
    integer :: nIndependentHam

    select case(nSpin)
    case(1,4)
      nIndependentHam = 1
    case(2)
      nIndependentHam = 2
    end select

  end function nIndependentHam


  !> Maxium occupation for single particle states
  pure function maximumFillings(nSpin)

    !> Number of spin channels (1=no-spin, 2=colinear, 4=non-collinear)
    integer, intent(in) :: nSpin

    !> Resulting upper limit on occupation for a level
    real(dp) :: maximumFillings

    select case(nSpin)
    case(1)
      maximumFillings = 2.0_dp
    case(2,4)
      maximumFillings = 1.0_dp
    end select

  end function maximumFillings


  !> Determines number of filled and empty states (allowing for fractional occupation)
  subroutine filledOrEmpty(nFilled, nEmpty, nIndepHam, nKpts, filling, nOrbs, maxFill)

    !> Number of (at least partially) filled levels for each spin and k-point
    integer, allocatable, intent(out) :: nFilled(:,:)

    !> Number of (at least partially) empty levels for each spin and k-point
    integer, allocatable, intent(out) :: nEmpty(:,:)

    !> Number of independent hamiltonians
    integer, intent(in) :: nIndepHam

    !> Number of k-points
    integer, intent(in) :: nKpts

    !> Number of atomic orbitals
    integer, intent(in) :: nOrbs

    !> Occupations
    real(dp), intent(in) :: filling(:,:,:)

    !> Upper limit on occupation for a level
    real(dp), intent(in) :: maxFill

    integer :: ik, iLev, iS

    allocate(nFilled(nIndepHam, nKpts), source=-1)
    do iS = 1, nIndepHam
      do iK = 1, nKPts
        do iLev = 1, nOrbs
          if ( filling(iLev, iK, iS) < epsilon(1.0_rsp) ) then
            ! assumes Fermi filling, so above this is empty
            nFilled(iS, iK) = iLev - 1
            exit
          end if
        end do
        ! check if channel is fully filled, if the loop above found nothing
        if (nFilled(iS, iK) < 0) then
          nFilled(iS, iK) = nOrbs
        end if
      end do
    end do

    allocate(nEmpty(nIndepHam, nKpts), source=-1)
    do iS = 1, nIndepHam
      do iK = 1, nKpts
        do iLev = 1, nOrbs
          if ( abs( filling(iLev, iK, iS) - maxFill ) > epsilon(1.0_rsp)) then
            ! assumes Fermi filling, so this is filled
            nEmpty(iS, iK) = iLev
            exit
          end if
        end do
        ! check is channel is empty, if the loop above found nothing
        if (nEmpty(iS, iK) < 0) then
          nEmpty(iS, iK) = 1
        end if
      end do
    end do

  end subroutine filledOrEmpty

end module dftbp_derivs_fillings
