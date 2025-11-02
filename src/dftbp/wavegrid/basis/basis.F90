!--------------------------------------------------------------------------------------------------!
!  DFTB+: general package for performing fast atomistic simulations                                !
!  Copyright (C) 2006 - 2025  DFTB+ developers group                                               !
!                                                                                                  !
!  See the LICENSE file for terms of usage and distribution.                                       !
!--------------------------------------------------------------------------------------------------!

module dftbp_wavegrid_basis
  use dftbp_common_accuracy, only : dp
  use dftbp_wavegrid_basis_spharmonics, only: realTessY
  use dftbp_wavegrid_basis_orbital, only: TOrbital, TOrbitalWrapper
  use dftbp_wavegrid_basis_slater, only: TSlaterOrbital, TSlaterOrbital_init
  use dftbp_wavegrid_basis_lut, only: TRadialTableOrbital, TRadialTableOrbital_initFromArray, &
       & TRadialTableOrbital_initFromOrbital
  implicit none

  private

  public :: realTessY
  public :: TOrbital, TOrbitalWrapper
  public :: TSlaterOrbital, TRadialTableOrbital

  public :: TSlaterOrbital_init
  public :: TRadialTableOrbital_initFromArray, TRadialTableOrbital_initFromOrbital

end module dftbp_wavegrid_basis

