!--------------------------------------------------------------------------------------------------!
!  DFTB+: general package for performing fast atomistic simulations                                !
!  Copyright (C) 2006 - 2020  DFTB+ developers group                                               !
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
  use dftbp_dispdftd4
  implicit none
  public


  !> Types of dispersion model
  type :: TDispersionInp

    !> Based on universal force-field
    type(TDispUffInp), allocatable :: uff

    !> Slater-Kirkwood
    type(TDispSlaKirkInp), allocatable :: slakirk

  #:if WITH_DFTD3
    !> Grimme DFT-D3
    type(TDispDftD3Inp), allocatable :: dftd3
  #:endif

    !> D4 dispersion model.
    type(TDispDftD4Inp), allocatable :: dftd4
  
  end type TDispersionInp

end module dftbp_dispersions
