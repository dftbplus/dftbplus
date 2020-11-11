!--------------------------------------------------------------------------------------------------!
!  DFTB+: general package for performing fast atomistic simulations                                !
!  Copyright (C) 2006 - 2020  DFTB+ developers group                                               !
!                                                                                                  !
!  See the LICENSE file for terms of usage and distribution.                                       !
!--------------------------------------------------------------------------------------------------!

#:include 'common.fypp'

!> Atomic masses for each element known in the PSE
module dftbp_atomicmass
  use dftbp_accuracy, only : dp
  use dftbp_constants, only : amu__au, symbolToNumber
  implicit none
  private

  public :: getAtomicMass


  !> Get atomic mass for a species
  interface getAtomicMass
    module procedure :: getAtomicMassSymbol
    module procedure :: getAtomicMassNumber
  end interface getAtomicMass


  !> Atomic masses
  real(dp), parameter :: atomicMassNist(1:118) = [ &
      &   1.00794075_dp,   4.00260193_dp,   6.94003660_dp,   9.01218307_dp,&
      &  10.81102805_dp,  12.01073590_dp,  14.00670321_dp,  15.99940492_dp,&
      &  18.99840316_dp,  20.18004638_dp,  22.98976928_dp,  24.30505162_dp,&
      &  26.98153853_dp,  28.08549871_dp,  30.97376200_dp,  32.06478741_dp,&
      &  35.45293758_dp,  39.94779856_dp,  39.09830091_dp,  40.07802251_dp,&
      &  44.95590828_dp,  47.86674496_dp,  50.94146504_dp,  51.99613176_dp,&
      &  54.93804391_dp,  55.84514443_dp,  58.93319429_dp,  58.69334711_dp,&
      &  63.54603995_dp,  65.37778253_dp,  69.72306607_dp,  72.62755016_dp,&
      &  74.92159457_dp,  78.95938856_dp,  79.90352778_dp,  83.79800000_dp,&
      &  85.46766360_dp,  87.61664447_dp,  88.90584030_dp,  91.22364160_dp,&
      &  92.90637300_dp,  95.95978854_dp,  97.90721240_dp, 101.06494014_dp,&
      & 102.90549800_dp, 106.41532751_dp, 107.86814963_dp, 112.41155782_dp,&
      & 114.81808663_dp, 118.71011259_dp, 121.75978367_dp, 127.60312648_dp,&
      & 126.90447190_dp, 131.29276145_dp, 132.90545196_dp, 137.32689163_dp,&
      & 138.90546887_dp, 140.11573074_dp, 140.90765760_dp, 144.24159603_dp,&
      & 144.91275590_dp, 150.36635571_dp, 151.96437813_dp, 157.25213065_dp,&
      & 158.92535470_dp, 162.49947282_dp, 164.93032880_dp, 167.25908265_dp,&
      & 168.93421790_dp, 173.05415017_dp, 174.96681496_dp, 178.48497872_dp,&
      & 180.94787564_dp, 183.84177755_dp, 186.20670455_dp, 190.22485963_dp,&
      & 192.21605165_dp, 195.08445686_dp, 196.96656879_dp, 200.59916703_dp,&
      & 204.38341284_dp, 207.21690806_dp, 208.98039910_dp, 208.98243080_dp,&
      & 209.98714790_dp, 222.01757820_dp, 223.01973600_dp, 226.02541030_dp,&
      & 227.02775230_dp, 232.03805580_dp, 231.03588420_dp, 238.02891046_dp,&
      & 237.04817360_dp, 244.06420530_dp, 243.06138130_dp, 247.07035410_dp,&
      & 247.07030730_dp, 251.07958860_dp, 252.08298000_dp, 257.09510610_dp,&
      & 258.09843150_dp, 259.10103000_dp, 262.10961000_dp, 267.12179000_dp,&
      & 269.12791000_dp, 271.13393000_dp, 270.13336000_dp, 276.14846000_dp,&
      & 276.15159000_dp, 280.16131000_dp, 282.16912000_dp, 284.17416000_dp,&
      & 284.17873000_dp, 289.19042000_dp, 288.19274000_dp, 293.20449000_dp,&
      & 292.20746000_dp, 294.21392000_dp]


contains


  !> Get atomic mass for species with a given symbol
  elemental function getAtomicMassSymbol(symbol) result(mass)

    !> Element symbol
    character(len=*), intent(in) :: symbol

    !> atomic mass
    real(dp) :: mass

    mass = getAtomicMass(symbolToNumber(symbol))

  end function getAtomicMassSymbol


  !> Get atomic mass for species with a given atomic number
  elemental function getAtomicMassNumber(number) result(mass)

    !> Atomic number
    integer, intent(in) :: number

    !> atomic mass
    real(dp) :: mass

    if (number > 0 .and. number <= size(atomicMassNist, dim=1)) then
      mass = atomicMassNist(number)
    else
      mass = -1.0_dp
    end if

  end function getAtomicMassNumber


end module dftbp_atomicmass
