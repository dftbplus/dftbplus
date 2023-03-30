!--------------------------------------------------------------------------------------------------!
!  DFTB+: general package for performing fast atomistic simulations                                !
!  Copyright (C) 2006 - 2023  DFTB+ developers group                                               !
!                                                                                                  !
!  See the LICENSE file for terms of usage and distribution.                                       !
!--------------------------------------------------------------------------------------------------!

!> Common format strings
module dftbp_io_commonformats
  use dftbp_common_accuracy, only : dp
  implicit none


  private
  public :: formatHessian, formatBorn, formatdBorn, formatGeoOut, format1U, format1Ue,&
      & format2Ue, format1U1e
  public :: strFormat2U


  !> Format string for energy second derivative matrix
  character(len=*), parameter :: formatHessian = '(4f16.10)'

  !> Format string for Born charges
  character(len=*), parameter :: formatBorn = '(3f16.10)'

  !> Format string for derivative of Born charges
  character(len=*), parameter :: formatdBorn = '(9E16.8)'

  !> Atomic geometries format
  character(len=*), parameter :: formatGeoOut = "(I5, F16.8, F16.8, F16.8)"

  !> Format for a single value with units
  character(len=*), parameter :: format1U = "(A, ':', T32, F18.10, T51, A)"

  !> Format for a single value using exponential notation with units
  character(len=*), parameter :: format1Ue = "(A, ':', T37, E13.6, T51, A)"

  !> Format for two using exponential notation values with units
  character(len=*), parameter :: format2Ue = "(A, ':', T37, E13.6, T51, A, T57, E13.6, T71, A)"

  !> Format for mixed decimal and exponential values with units
  character(len=*), parameter :: format1U1e =&
      & "(' ', A, ':', T32, F18.10, T51, A, T57, E13.6, T71, A)"

contains


  !> Formats a string using the "(a, ':', t32, f18.10, t51, a, t54, f16.4, t71, a)" format.
  pure function strFormat2U(quantity, value1, unit1, value2, unit2) result(str)

    !> Name of the quantity to print (first character in the output)
    character(*), intent(in) :: quantity

    !> First value
    real(dp), intent(in) :: value1

    !> First unit
    character(*), intent(in) :: unit1

    !> Second value
    real(dp), intent(in) :: value2

    !> Second unit
    character(*), intent(in) :: unit2

    !> Formatted string
    character(70 + len(unit2)) :: str

    integer :: itemLen

    str = repeat(" ", len(str))
    itemLen = min(len(quantity), 30)
    str(1 : itemLen) = quantity(: itemLen)
    str(itemLen + 1 : itemLen + 1) = ":"
    write(str(32 : 49), "(F18.10)") value1
    itemLen = min(len(unit1), 3)
    str(51 : 51 + itemLen - 1) = unit1(: itemLen)
    write(str(54 : 69), "(F16.4)") value2
    str(71 : len(str)) = unit2

  end function strformat2U

end module dftbp_io_commonformats
