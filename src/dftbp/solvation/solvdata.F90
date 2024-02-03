!--------------------------------------------------------------------------------------------------!
!  DFTB+: general package for performing fast atomistic simulations                                !
!  Copyright (C) 2006 - 2024  DFTB+ developers group                                               !
!                                                                                                  !
!  See the LICENSE file for terms of usage and distribution.                                       !
!--------------------------------------------------------------------------------------------------!

#:include 'common.fypp'

!> Parameters for the generalized Born model
module dftbp_solvation_solvdata
  use dftbp_common_accuracy, only : dp
  use dftbp_common_constants, only : AA__Bohr, symbolToNumber
  implicit none
  private

  public :: getVanDerWaalsRadiusD3, getVanDerWaalsRadiusCosmo, getVanDerWaalsRadiusBondi


  !> In case no van-der-Waals value is provided
  real(dp), parameter :: missing = -1.0_dp

  !> Get van-der-Waals radius for a species
  interface getVanDerWaalsRadiusD3
    module procedure :: getVanDerWaalsRadiusD3Symbol
    module procedure :: getVanDerWaalsRadiusD3Number
  end interface getVanDerWaalsRadiusD3

  !> D3 pairwise van-der-Waals radii (only homoatomic pairs present here)
  real(dp), parameter :: vanDerWaalsRadiiD3(1:94) = AA__Bohr * [&
      & 1.09155_dp, 0.86735_dp, 1.74780_dp, 1.54910_dp, &  ! H-Be
      & 1.60800_dp, 1.45515_dp, 1.31125_dp, 1.24085_dp, &  ! B-O
      & 1.14980_dp, 1.06870_dp, 1.85410_dp, 1.74195_dp, &  ! F-Mg
      & 2.00530_dp, 1.89585_dp, 1.75085_dp, 1.65535_dp, &  ! Al-S
      & 1.55230_dp, 1.45740_dp, 2.12055_dp, 2.05175_dp, &  ! Cl-Ca
      & 1.94515_dp, 1.88210_dp, 1.86055_dp, 1.72070_dp, &  ! Sc-Cr
      & 1.77310_dp, 1.72105_dp, 1.71635_dp, 1.67310_dp, &  ! Mn-Ni
      & 1.65040_dp, 1.61545_dp, 1.97895_dp, 1.93095_dp, &  ! Cu-Ge
      & 1.83125_dp, 1.76340_dp, 1.68310_dp, 1.60480_dp, &  ! As-Kr
      & 2.30880_dp, 2.23820_dp, 2.10980_dp, 2.02985_dp, &  ! Rb-Zr
      & 1.92980_dp, 1.87715_dp, 1.78450_dp, 1.73115_dp, &  ! Nb-Ru
      & 1.69875_dp, 1.67625_dp, 1.66540_dp, 1.73100_dp, &  ! Rh-Cd
      & 2.13115_dp, 2.09370_dp, 2.00750_dp, 1.94505_dp, &  ! In-Te
      & 1.86900_dp, 1.79445_dp, 2.52835_dp, 2.59070_dp, &  ! I-Ba
      & 2.31305_dp, 2.31005_dp, 2.28510_dp, 2.26355_dp, &  ! La-Nd
      & 2.24480_dp, 2.22575_dp, 2.21170_dp, 2.06215_dp, &  ! Pm-Gd
      & 2.12135_dp, 2.07705_dp, 2.13970_dp, 2.12250_dp, &  ! Tb-Er
      & 2.11040_dp, 2.09930_dp, 2.00650_dp, 2.12250_dp, &  ! Tm-Hf
      & 2.04900_dp, 1.99275_dp, 1.94775_dp, 1.87450_dp, &  ! Ta-Os
      & 1.72280_dp, 1.67625_dp, 1.62820_dp, 1.67995_dp, &  ! Ir-Hg
      & 2.15635_dp, 2.13820_dp, 2.05875_dp, 2.00270_dp, &  ! Tl-Po
      & 1.93220_dp, 1.86080_dp, 2.53980_dp, 2.46470_dp, &  ! At-Ra
      & 2.35215_dp, 2.21260_dp, 2.22970_dp, 2.19785_dp, &  ! Ac-U
      & 2.17695_dp, 2.21705_dp]                            ! Np-Pu


  !> Get van-der-Waals radius for a species
  interface getVanDerWaalsRadiusCosmo
    module procedure :: getVanDerWaalsRadiusCosmoSymbol
    module procedure :: getVanDerWaalsRadiusCosmoNumber
  end interface getVanDerWaalsRadiusCosmo


   !> Default value for unoptimized van-der-Waals radii
   real(dp), parameter :: cosmoStub = 2.223_dp

   !> COSMO optimized van-der-Waals radii
   real(dp), parameter :: vanDerWaalsRadiiCosmo(1:94) = AA__Bohr * [ &
       & 1.3000_dp, 1.6380_dp, 1.5700_dp, 1.0530_dp, &   ! H-Be
       & 2.0480_dp, 2.0000_dp, 1.8300_dp, 1.7200_dp, &   ! B-O
       & 1.7200_dp, 1.8018_dp, 1.8000_dp, 1.6380_dp, &   ! F-Mg
       & 2.1530_dp, 2.2000_dp, 2.1060_dp, 2.1600_dp, &   ! Al-S
       & 2.0500_dp, 2.2000_dp, 2.2230_dp, cosmoStub, &   ! Cl-Ca
       & cosmoStub, 2.2930_dp, cosmoStub, cosmoStub, &   ! Sc-Cr
       & cosmoStub, cosmoStub, cosmoStub, cosmoStub, &   ! Mn-Ni
       & cosmoStub, 1.6260_dp, cosmoStub, 2.7000_dp, &   ! Cu-Ge
       & 2.3500_dp, 2.2000_dp, 2.1600_dp, 2.3630_dp, &   ! As-Kr
       & cosmoStub, cosmoStub, cosmoStub, cosmoStub, &   ! Rb-Zr
       & cosmoStub, cosmoStub, cosmoStub, cosmoStub, &   ! Nb-Ru
       & cosmoStub, cosmoStub, cosmoStub, cosmoStub, &   ! Rh-Cd
       & 2.2580_dp, 2.5500_dp, 2.4100_dp, 2.4100_dp, &   ! In-Te
       & 2.3200_dp, 2.5270_dp, cosmoStub, cosmoStub, &   ! I-Ba
       & cosmoStub, cosmoStub, cosmoStub, cosmoStub, &   ! La-Nd
       & cosmoStub, cosmoStub, cosmoStub, cosmoStub, &   ! Pm-Gd
       & cosmoStub, cosmoStub, cosmoStub, cosmoStub, &   ! Tb-Er
       & cosmoStub, cosmoStub, cosmoStub, cosmoStub, &   ! Tm-Hf
       & cosmoStub, cosmoStub, cosmoStub, cosmoStub, &   ! Ta-Os
       & cosmoStub, cosmoStub, cosmoStub, cosmoStub, &   ! Ir-Hg
       & cosmoStub, 2.3600_dp, 2.4220_dp, 2.3050_dp, &   ! Tl-Po
       & 2.3630_dp, 2.5740_dp, cosmoStub, cosmoStub, &   ! At-Ra
       & cosmoStub, cosmoStub, cosmoStub, cosmoStub, &   ! Ac-U
       & cosmoStub, cosmoStub]                           ! Np-Pu


  !> Get van-der-Waals radius for a species
  interface getVanDerWaalsRadiusBondi
    module procedure :: getVanDerWaalsRadiusBondiSymbol
    module procedure :: getVanDerWaalsRadiusBondiNumber
  end interface getVanDerWaalsRadiusBondi

  !> Van-der-Waals radii from
  !> Manjeera Mantina, Adam C. Chamberlin, Rosendo Valero, Christopher J. Cramer,
  !> and Donald G. Truhlar, Consistent van der Waals Radii for the Whole Main Group,
  !> J. Phys. Chem. A 2009, 113, 19, 5806–5812. https://doi.org/10.1021/jp8111556
  real(dp), parameter :: vanDerWaalsRadiiBondi(1:88) = AA__Bohr * [ &
      & 1.10_dp, 1.40_dp, 1.81_dp, 1.53_dp, 1.92_dp, 1.70_dp, 1.55_dp, 1.52_dp, &  ! H-O
      & 1.47_dp, 1.54_dp, 2.27_dp, 1.73_dp, 1.84_dp, 2.10_dp, 1.80_dp, 1.80_dp, &  ! F-S
      & 1.75_dp, 1.88_dp, 2.75_dp, 2.31_dp, missing, missing, missing, missing, &  ! Cl-Cr
      & missing, missing, missing, missing, missing, missing, 1.87_dp, 2.11_dp, &  ! Mn-Ge
      & 1.85_dp, 1.90_dp, 1.83_dp, 2.02_dp, 3.03_dp, 2.49_dp, missing, missing, &  ! As-Zr
      & missing, missing, missing, missing, missing, missing, missing, missing, &  ! Nb-Cd
      & 1.93_dp, 2.17_dp, 2.06_dp, 2.06_dp, 1.98_dp, 2.16_dp, 3.43_dp, 2.68_dp, &  ! I-Ba
      & missing, missing, missing, missing, missing, missing, missing, missing, &  ! La-Gd
      & missing, missing, missing, missing, missing, missing, missing, missing, &  ! Tb-Hf
      & missing, missing, missing, missing, missing, missing, missing, missing, &  ! Ta-Hg
      & 1.96_dp, 2.02_dp, 2.07_dp, 1.97_dp, 2.02_dp, 2.20_dp, 3.48_dp, 2.83_dp]    ! Tl-Ra

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
      radius = missing
    end if

  end function getVanDerWaalsRadiusD3Number

  !> Get van-der-Waals radius for species with a given symbol
  elemental function getVanDerWaalsRadiusCosmoSymbol(symbol) result(radius)

    !> Element symbol
    character(len=*), intent(in) :: symbol

    !> van-der-Waals radius
    real(dp) :: radius

    radius = getVanDerWaalsRadiusCosmo(symbolToNumber(symbol))

  end function getVanDerWaalsRadiusCosmoSymbol

  !> Get van-der-Waals radius for species with a given atomic number
  elemental function getVanDerWaalsRadiusCosmoNumber(number) result(radius)

    !> Atomic number
    integer, intent(in) :: number

    !> van-der-Waals radius
    real(dp) :: radius

    if (number > 0 .and. number <= size(vanDerWaalsRadiiCosmo, dim=1)) then
      radius = vanDerWaalsRadiiCosmo(number)
    else
      radius = missing
    end if

  end function getVanDerWaalsRadiusCosmoNumber

  !> Get van-der-Waals radius for species with a given symbol
  elemental function getVanDerWaalsRadiusBondiSymbol(symbol) result(radius)

    !> Element symbol
    character(len=*), intent(in) :: symbol

    !> van-der-Waals radius
    real(dp) :: radius

    radius = getVanDerWaalsRadiusBondi(symbolToNumber(symbol))

  end function getVanDerWaalsRadiusBondiSymbol

  !> Get van-der-Waals radius for species with a given atomic number
  elemental function getVanDerWaalsRadiusBondiNumber(number) result(radius)

    !> Atomic number
    integer, intent(in) :: number

    !> van-der-Waals radius
    real(dp) :: radius

    if (number > 0 .and. number <= size(vanDerWaalsRadiiBondi, dim=1)) then
      radius = vanDerWaalsRadiiBondi(number)
    else
      radius = missing
    end if

  end function getVanDerWaalsRadiusBondiNumber

end module dftbp_solvation_solvdata
