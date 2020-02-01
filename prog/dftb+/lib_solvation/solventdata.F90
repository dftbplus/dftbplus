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
  implicit none
  private

  public :: TSolventData

  type :: TSolventData
    real(dp) :: dielectricConstant
    real(dp) :: molecularMass
    real(dp) :: density
  contains
    procedure :: fromName
  end type TSolventData

contains

  subroutine fromName(self, name, found)
    class(TSolventData), intent(out) :: self
    character(len=*), intent(in) :: name
    logical, intent(out) :: found

    found = .false.

    select case(name)
    case default
      return
    end select

    found = .true.

  end subroutine fromName

end module dftbp_solventdata
