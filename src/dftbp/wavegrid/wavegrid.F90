!--------------------------------------------------------------------------------------------------!
!  DFTB+: general package for performing fast atomistic simulations                                !
!  Copyright (C) 2006 - 2025  DFTB+ developers group                                               !
!                                                                                                  !
!  See the LICENSE file for terms of usage and distribution.                                       !
!--------------------------------------------------------------------------------------------------!

module dftbp_wavegrid
  use dftbp_wavegrid_molorb, only: TSpeciesBasis, TMolecularOrbital, TMolecularOrbital_init
  use dftbp_wavegrid_molorb, only: getValue, getAtomicDensities, getTotalChrg
  use dftbp_wavegrid_molorb_spharmonics, only: realTessY
  use dftbp_wavegrid_slater, only: TSlaterOrbital
  implicit none

  public :: TSpeciesBasis, TMolecularOrbital, TSlaterOrbital
  public :: TMolecularOrbital_init
  public :: realTessY
  public :: getValue, getAtomicDensities, getTotalChrg

end module dftbp_wavegrid

