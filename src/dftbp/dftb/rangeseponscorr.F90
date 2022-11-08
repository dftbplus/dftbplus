!--------------------------------------------------------------------------------------------------!
!  DFTB+: general package for performing fast atomistic simulations                                !
!  Copyright (C) 2006 - 2022  DFTB+ developers group                                               !
!                                                                                                  !
!  See the LICENSE file for terms of usage and distribution.                                       !
!--------------------------------------------------------------------------------------------------!

#:include 'common.fypp'

!> Contains subroutines to add contribution of range-separated hybrid functional to onsite correction
!> from doi: 10.1021/acs.jctc.2c00037
module dftbp_dftb_rangeseponscorr
  implicit none

  private
  public :: TRangeSepOnsCorrFunc


  !> Onsite correction with range-separated hybrid functional module
  type :: TRangeSepOnsCorrFunc
    private


  end type TRangeSepOnsCorrFunc



end module dftbp_dftb_rangeseponscorr
