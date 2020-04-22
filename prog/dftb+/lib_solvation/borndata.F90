!--------------------------------------------------------------------------------------------------!
!  DFTB+: general package for performing fast atomistic simulations                                !
!  Copyright (C) 2006 - 2020  DFTB+ developers group                                               !
!                                                                                                  !
!  See the LICENSE file for terms of usage and distribution.                                       !
!--------------------------------------------------------------------------------------------------!

#:include 'common.fypp'

!> Parameters for the generalized Born model
module dftbp_borndata
  use dftbp_accuracy, only : dp
  use dftbp_constants, only : AA__Bohr, symbolToNumber
  implicit none
  private

  public :: getVanDerWaalsRadiusD3


  !> Get van-der-Waals radius for a species
  interface getVanDerWaalsRadiusD3
    module procedure :: getVanDerWaalsRadiusD3Symbol
    module procedure :: getVanDerWaalsRadiusD3Number
  end interface getVanDerWaalsRadiusD3

  !> D3 pairwise van-der-Waals radii (only homoatomic pairs present here)
  real(dp), parameter :: vanDerWaalsRadiiD3(1:94) = AA__Bohr * [&
      & 1.09155_dp, 0.86735_dp, 1.74780_dp, 1.54910_dp, &
      & 1.60800_dp, 1.45515_dp, 1.31125_dp, 1.24085_dp, &
      & 1.14980_dp, 1.06870_dp, 1.85410_dp, 1.74195_dp, &
      & 2.00530_dp, 1.89585_dp, 1.75085_dp, 1.65535_dp, &
      & 1.55230_dp, 1.45740_dp, 2.12055_dp, 2.05175_dp, &
      & 1.94515_dp, 1.88210_dp, 1.86055_dp, 1.72070_dp, &
      & 1.77310_dp, 1.72105_dp, 1.71635_dp, 1.67310_dp, &
      & 1.65040_dp, 1.61545_dp, 1.97895_dp, 1.93095_dp, &
      & 1.83125_dp, 1.76340_dp, 1.68310_dp, 1.60480_dp, &
      & 2.30880_dp, 2.23820_dp, 2.10980_dp, 2.02985_dp, &
      & 1.92980_dp, 1.87715_dp, 1.78450_dp, 1.73115_dp, &
      & 1.69875_dp, 1.67625_dp, 1.66540_dp, 1.73100_dp, &
      & 2.13115_dp, 2.09370_dp, 2.00750_dp, 1.94505_dp, &
      & 1.86900_dp, 1.79445_dp, 2.52835_dp, 2.59070_dp, &
      & 2.31305_dp, 2.31005_dp, 2.28510_dp, 2.26355_dp, &
      & 2.24480_dp, 2.22575_dp, 2.21170_dp, 2.06215_dp, &
      & 2.12135_dp, 2.07705_dp, 2.13970_dp, 2.12250_dp, &
      & 2.11040_dp, 2.09930_dp, 2.00650_dp, 2.12250_dp, &
      & 2.04900_dp, 1.99275_dp, 1.94775_dp, 1.87450_dp, &
      & 1.72280_dp, 1.67625_dp, 1.62820_dp, 1.67995_dp, &
      & 2.15635_dp, 2.13820_dp, 2.05875_dp, 2.00270_dp, &
      & 1.93220_dp, 1.86080_dp, 2.53980_dp, 2.46470_dp, &
      & 2.35215_dp, 2.21260_dp, 2.22970_dp, 2.19785_dp, &
      & 2.17695_dp, 2.21705_dp]

contains

  !> Get van-der-Waals radius for species with a given symbol
  elemental function getVanDerWaalsRadiusD3Symbol(symbol) result(radius)

    !> Element symbol
    character(len=*), intent(in) :: symbol

    !> van-der-Waals radius
    real(dp) :: radius

    radius = getVanDerWaalsRadiusD3(symbolToNumber(symbol))

  end function getVanDerWaalsRadiusD3Symbol

  !> Get van-der-Waals radius for species with a given atomic number
  elemental function getVanDerWaalsRadiusD3Number(number) result(radius)

    !> Atomic number
    integer, intent(in) :: number

    !> van-der-Waals radius
    real(dp) :: radius

    if (number > 0 .and. number <= size(vanDerWaalsRadiiD3, dim=1)) then
      radius = vanDerWaalsRadiiD3(number)
    else
      radius = -1.0_dp
    end if

  end function getVanDerWaalsRadiusD3Number

end module dftbp_borndata
