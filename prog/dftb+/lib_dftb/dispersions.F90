!--------------------------------------------------------------------------------------------------!
!  DFTB+: general package for performing fast atomistic simulations                                !
!  Copyright (C) 2006 - 2021  DFTB+ developers group                                               !
!                                                                                                  !
!  See the LICENSE file for terms of usage and distribution.                                       !
!--------------------------------------------------------------------------------------------------!

#:include 'common.fypp'

!> Offers everything which is publicly available when dealing with dispersions.
!!
module dftbp_dispersions
  use dftbp_dispiface, only : TDispersionIface
  use dftbp_dispuff, only : TDispUffInp, TDispUFF, dispuff_init
  use dftbp_dispuffdata, only : getuffvalues
  use dftbp_dispslaterkirkw, only : TDispSlaKirkInp, TDispSlaKirk, dispslakirk_init
#:if WITH_DFTD3
  use dftbp_dispdftd3, only : TDispDftD3Inp
#:endif
  use dftbp_simpledftd3, only : TSimpleDftD3Input, TSimpleDftD3, init
  use dftbp_dispdftd4, only : TDispDftD4Inp, TDispDftD4, init
#:if WITH_MBD
  use dftbp_dispmbd, only : TDispMbdInp
#:endif
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

    !> Simple D3 dispersion model.
    type(TSimpleDftD3Input), allocatable :: sdftd3

    !> D4 dispersion model.
    type(TDispDftD4Inp), allocatable :: dftd4
  
  #:if WITH_MBD
    !> Many-body dispersion
    type(TDispMbdInp), allocatable :: mbd
  #:endif

  end type TDispersionInp

end module dftbp_dispersions
