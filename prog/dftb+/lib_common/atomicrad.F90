!--------------------------------------------------------------------------------------------------!
!  DFTB+: general package for performing fast atomistic simulations                                !
!  Copyright (C) 2006 - 2020  DFTB+ developers group                                               !
!                                                                                                  !
!  See the LICENSE file for terms of usage and distribution.                                       !
!--------------------------------------------------------------------------------------------------!

#:include 'common.fypp'

!> Atomic Radii of the Elements
!>
!> M. Mantina, R. Valero, C. J. Cramer, and D. G. Truhlar,
!> in CRC Handbook of Chemistry and Physics, 91st Edition (2010-2011),
!> edited by W. M. Haynes (CRC Press, Boca Raton, FL, 2010), pages 9-49-9-50;
!> corrected Nov. 17, 2010 for the 92nd edition.
module dftbp_atomicrad
  use dftbp_accuracy, only : dp
  use dftbp_constants, only : AA__Bohr, symbolToNumber
  implicit none
  private

  public :: getAtomicRad


  !> Get atomic radius for a species
  interface getAtomicRad
    module procedure :: getAtomicRadSymbol
    module procedure :: getAtomicRadNumber
  end interface getAtomicRad


  !> Atomic radii
  real(dp), parameter :: atomicRadii(1:118) = AA__Bohr * [ &
      & 0.32_dp, 0.37_dp, 1.30_dp, 0.99_dp, 0.84_dp, 0.75_dp, 0.71_dp, 0.64_dp, &
      & 0.60_dp, 0.62_dp, 1.60_dp, 1.40_dp, 1.24_dp, 1.14_dp, 1.09_dp, 1.04_dp, &
      & 1.00_dp, 1.01_dp, 2.00_dp, 1.74_dp, 1.59_dp, 1.48_dp, 1.44_dp, 1.30_dp, &
      & 1.29_dp, 1.24_dp, 1.18_dp, 1.17_dp, 1.22_dp, 1.20_dp, 1.23_dp, 1.20_dp, &
      & 1.20_dp, 1.18_dp, 1.17_dp, 1.16_dp, 2.15_dp, 1.90_dp, 1.76_dp, 1.64_dp, &
      & 1.56_dp, 1.46_dp, 1.38_dp, 1.36_dp, 1.34_dp, 1.30_dp, 1.36_dp, 1.40_dp, &
      & 1.42_dp, 1.40_dp, 1.40_dp, 1.37_dp, 1.36_dp, 1.36_dp, 2.38_dp, 2.06_dp, &
      & 1.94_dp, 1.84_dp, 1.90_dp, 1.88_dp, 1.86_dp, 1.85_dp, 1.83_dp, 1.82_dp, &
      & 1.81_dp, 1.80_dp, 1.79_dp, 1.77_dp, 1.77_dp, 1.78_dp, 1.74_dp, 1.64_dp, &
      & 1.58_dp, 1.50_dp, 1.41_dp, 1.36_dp, 1.32_dp, 1.30_dp, 1.30_dp, 1.32_dp, &
      & 1.44_dp, 1.45_dp, 1.50_dp, 1.42_dp, 1.48_dp, 1.46_dp, 2.42_dp, 2.11_dp, &
      & 2.01_dp, 1.90_dp, 1.84_dp, 1.83_dp, 1.80_dp, 1.80_dp, 1.73_dp, 1.68_dp, &
      & 1.68_dp, 1.68_dp, 1.65_dp, 1.67_dp, 1.73_dp, 1.76_dp, 1.61_dp, 1.57_dp, &
      & 1.49_dp, 1.43_dp, 1.41_dp, 1.34_dp, 1.29_dp, 1.28_dp, 1.21_dp, 1.22_dp, &
      & 1.36_dp, 1.43_dp, 1.62_dp, 1.75_dp, 1.65_dp, 1.57_dp]


contains


  !> Get atomic radius for species with a given symbol
  elemental function getAtomicRadSymbol(symbol) result(radius)

    !> Element symbol
    character(len=*), intent(in) :: symbol

    !> atomic radius
    real(dp) :: radius

    radius = getAtomicRad(symbolToNumber(symbol))

  end function getAtomicRadSymbol


  !> Get atomic radius for species with a given atomic number
  elemental function getAtomicRadNumber(number) result(radius)

    !> Atomic number
    integer, intent(in) :: number

    !> atomic radius
    real(dp) :: radius

    if (number > 0 .and. number <= size(atomicRadii, dim=1)) then
      radius = AtomicRadii(number)
    else
      radius = -1.0_dp
    end if

  end function getAtomicRadNumber


end module dftbp_atomicrad
