!--------------------------------------------------------------------------------------------------!
!  DFTB+: general package for performing fast atomistic simulations                                !
!  Copyright (C) 2006 - 2020  DFTB+ developers group                                               !
!                                                                                                  !
!  See the LICENSE file for terms of usage and distribution.                                       !
!--------------------------------------------------------------------------------------------------!

!> Interface wrapper for the DFT-D4 parameters
module dftbp_dftd4refs
  use dftbp_accuracy, only : dp
  implicit none
  private

  public :: sscale, clsq, secaiw, alphaiw, clsh, refn, refsys, refcn, refcovcn
  public :: hcount, ascale

  include 'dftd4_references.fh'

end module dftbp_dftd4refs
