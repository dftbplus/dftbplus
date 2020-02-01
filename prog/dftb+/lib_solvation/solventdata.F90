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
  implicit none
  private

  public :: TSolventData
  public :: SolventFromName


  type :: TSolventData
    real(dp) :: dielectricConstant
    real(dp) :: molecularMass
    real(dp) :: density
  end type TSolventData


contains


  subroutine solventFromName(self, name, found)
    type(TSolventData), intent(out) :: self
    character(len=*), intent(in) :: name
    logical, intent(out) :: found

    found = .false.

    select case(tolower(name))
    case default
      return
    case("acetone")
      self = TSolventData(20.7_dp, 58.08_dp, 0.79_dp)
    case("acetonitrile")
      self = TSolventData(37.5_dp, 41.05_dp, 0.786_dp)
    case("toluene")
      self = TSolventData(7.0_dp, 92.14_dp, 0.867_dp)
    case("benzene", "c6h6")
      self = TSolventData(7.0_dp, 78.11_dp, 0.867_dp)
    case("chloroform", "trichloromethane", "chcl3")
      self = TSolventData(7.0_dp, 119.38_dp, 1.49_dp)
    case("dichloromethane", "ch2cl2")
      self = TSolventData(7.0_dp, 84.93_dp, 1.33_dp)
    case("cs2")
      self = TSolventData(2.6_dp, 76.13_dp, 1.266_dp)
    case("dmf")
      self = TSolventData(37.0_dp, 73.1_dp, 0.95_dp)
    case("dmso")
      self = TSolventData(47.2_dp, 78.13_dp, 1.1_dp)
    case("ether")
      self = TSolventData(7.3_dp, 74.12_dp, 0.713_dp)
    case("water", "h2o")
      self = TSolventData(78.5_dp, 18.0_dp, 0.998_dp)
    case("methanole")
      self = TSolventData(33.6_dp, 32.04_dp, 0.792_dp)
    case("nhexane", "n-hexane")
      self = TSolventData(1.88_dp, 86.18_dp, 0.66_dp)
    case("thf")
      self = TSolventData(10.0_dp, 72.1061_dp, 0.883_dp)
    end select

    found = .true.

  end subroutine solventFromName

end module dftbp_solventdata
