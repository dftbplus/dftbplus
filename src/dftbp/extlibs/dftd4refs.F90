!--------------------------------------------------------------------------------------------------!
!  DFTB+: general package for performing fast atomistic simulations                                !
!  Copyright (C) 2006 - 2023  DFTB+ developers group                                               !
!                                                                                                  !
!  See the LICENSE file for terms of usage and distribution.                                       !
!--------------------------------------------------------------------------------------------------!

!> Interface wrapper for the DFT-D4 parameters
module dftbp_extlibs_dftd4refs
  use dftbp_common_accuracy, only : dp
  implicit none
  private

  public :: sscale, clsq, secaiw, alphaiw, clsh, refn, refsys, refcn, refcovcn
  public :: hcount, ascale

  include 'dftd4_references.fh'

end module dftbp_extlibs_dftd4refs
