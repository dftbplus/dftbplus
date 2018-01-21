!--------------------------------------------------------------------------------------------------!
!  DFTB+: general package for performing fast atomistic simulations                                !
!  Copyright (C) 2017  DFTB+ developers group                                                      !
!                                                                                                  !
!  See the LICENSE file for terms of usage and distribution.                                       !
!--------------------------------------------------------------------------------------------------!

#:include 'common.fypp'

!> Offers everything which is publicly available when dealing with dispersions.
!!
module dispersions
  use dispiface
  use dispuff_module
  use dispuffdata
  use dispslaterkirkw
#:if WITH_DFTD3
  use dispdftd3_module
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

end module dispersions
