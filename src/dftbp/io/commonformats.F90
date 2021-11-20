!--------------------------------------------------------------------------------------------------!
!  DFTB+: general package for performing fast atomistic simulations                                !
!  Copyright (C) 2006 - 2022  DFTB+ developers group                                               !
!                                                                                                  !
!  See the LICENSE file for terms of usage and distribution.                                       !
!--------------------------------------------------------------------------------------------------!

!> Common format strings
module dftbp_io_commonformats
  implicit none

  private
  public :: formatHessian, formatGeoOut, format1U, format2U, format1Ue, format2Ue, format1U1e

  !> Format string for energy second derivative matrix
  character(len=*), parameter :: formatHessian = '(4f16.10)'

  !> Atomic geometries format
  character(len=*), parameter :: formatGeoOut = "(I5, F16.8, F16.8, F16.8)"

  !> Format for a single value with units
  character(len=*), parameter :: format1U = "(A, ':', T32, F18.10, T51, A)"

  !> Format for two values with units
  character(len=*), parameter :: format2U = "(A, ':', T32, F18.10, T51, A, T54, F16.4, T71, A)"

  !> Format for a single value using exponential notation with units
  character(len=*), parameter :: format1Ue = "(A, ':', T37, E13.6, T51, A)"

  !> Format for two using exponential notation values with units
  character(len=*), parameter :: format2Ue = "(A, ':', T37, E13.6, T51, A, T57, E13.6, T71, A)"

  !> Format for mixed decimal and exponential values with units
  character(len=*), parameter :: format1U1e =&
      & "(' ', A, ':', T32, F18.10, T51, A, T57, E13.6, T71, A)"

end module dftbp_io_commonformats
