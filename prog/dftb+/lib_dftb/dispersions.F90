!--------------------------------------------------------------------------------------------------!
!  DFTB+: general package for performing fast atomistic simulations                                !
!  Copyright (C) 2006 - 2019  DFTB+ developers group                                               !
!                                                                                                  !
!  See the LICENSE file for terms of usage and distribution.                                       !
!--------------------------------------------------------------------------------------------------!

#:include 'common.fypp'

!> Offers everything which is publicly available when dealing with dispersions.
!!
module dftbp_dispersions
  use dftbp_dispiface
  use dftbp_dispuff
  use dftbp_dispuffdata
  use dftbp_dispslaterkirkw
#:if WITH_DFTD3
  use dftbp_dispdftd3
#:endif
  implicit none
  public


  !> Types of dispersion model
  type :: DispersionInp

    !> Based on universal force-field
    type(DispUffInp), allocatable :: uff

    !> Slater-Kirkwood
    type(DispSlaKirkInp), allocatable :: slakirk

  #:if WITH_DFTD3
    !> Grimme DFT-D3
    type(DispDftD3Inp), allocatable :: dftd3
  #:endif
  
  end type DispersionInp

end module dftbp_dispersions
