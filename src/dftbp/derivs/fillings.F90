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

  private
  public :: nIndependentHam, maximumFillings

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


end module dftbp_derivs_fillings
