!--------------------------------------------------------------------------------------------------!
!  DFTB+: general package for performing fast atomistic simulations                                !
!  Copyright (C) 2006 - 2020  DFTB+ developers group                                               !
!                                                                                                  !
!  See the LICENSE file for terms of usage and distribution.                                       !
!--------------------------------------------------------------------------------------------------!

#:include 'common.fypp'

!> Constants for commonly known solvents
module dftbp_solventdata
  use dftbp_accuracy, only : dp
  use dftbp_charmanip, only : tolower
  use dftbp_constants, only : kg__au, AA__Bohr, amu__au
  implicit none
  private

  public :: TSolventData
  public :: SolventFromName


  type :: TSolventData
    real(dp) :: dielectricConstant
    real(dp) :: molecularMass
    real(dp) :: density
  end type TSolventData

  real(dp), parameter :: kgPerL__au = 1.0e+3_dp*kg__au/(1.0e10_dp*AA__Bohr)**3


contains


  subroutine solventFromName(this, name, found)
    type(TSolventData), intent(out) :: this
    character(len=*), intent(in) :: name
    logical, intent(out) :: found

    found = .false.

    select case(tolower(name))
    case default
      return
    case("acetone")
      this = TSolventData(20.7_dp, 58.08_dp*amu__au, 0.79_dp*kgPerL__au)
    case("acetonitrile")
      this = TSolventData(37.5_dp, 41.05_dp*amu__au, 0.786_dp*kgPerL__au)
    case("toluene")
      this = TSolventData(7.0_dp, 92.14_dp*amu__au, 0.867_dp*kgPerL__au)
    case("benzene", "c6h6")
      this = TSolventData(7.0_dp, 78.11_dp*amu__au, 0.867_dp*kgPerL__au)
    case("chloroform", "trichloromethane", "chcl3")
      this = TSolventData(7.0_dp, 119.38_dp*amu__au, 1.49_dp*kgPerL__au)
    case("dichloromethane", "ch2cl2")
      this = TSolventData(7.0_dp, 84.93_dp*amu__au, 1.33_dp*kgPerL__au)
    case("cs2")
      this = TSolventData(2.6_dp, 76.13_dp*amu__au, 1.266_dp*kgPerL__au)
    case("dmf")
      this = TSolventData(37.0_dp, 73.1_dp*amu__au, 0.95_dp*kgPerL__au)
    case("dmso")
      this = TSolventData(47.2_dp, 78.13_dp*amu__au, 1.1_dp*kgPerL__au)
    case("ether")
      this = TSolventData(7.3_dp, 74.12_dp*amu__au, 0.713_dp*kgPerL__au)
    case("water", "h2o")
      this = TSolventData(78.5_dp, 18.0_dp*amu__au, 0.998_dp*kgPerL__au)
    case("methanole")
      this = TSolventData(33.6_dp, 32.04_dp*amu__au, 0.792_dp*kgPerL__au)
    case("nhexane", "n-hexane")
      this = TSolventData(1.88_dp, 86.18_dp*amu__au, 0.66_dp*kgPerL__au)
    case("thf")
      this = TSolventData(10.0_dp, 72.1061_dp*amu__au, 0.883_dp*kgPerL__au)
    end select

    found = .true.

  end subroutine solventFromName


end module dftbp_solventdata
