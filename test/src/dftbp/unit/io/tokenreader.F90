!--------------------------------------------------------------------------------------------------!
!  DFTB+: general package for performing fast atomistic simulations                                !
!  Copyright (C) 2006 - 2025  DFTB+ developers group                                               !
!                                                                                                  !
!  See the LICENSE file for terms of usage and distribution.                                       !
!--------------------------------------------------------------------------------------------------!

#:include "fortuno_serial.fypp"


module test_io_tokenreader
  use fortuno_serial, only : suite => serial_suite_item, test_list
  use dftbp_common_accuracy, only : dp
  use dftbp_io_tokenreader, only : getNextToken
  $:FORTUNO_SERIAL_IMPORTS()
  implicit none

  private
  public :: tests


  ! Tolerance for recognizing floating point nr. when parsing
  real(dp), parameter :: tol = 10.0_dp * epsilon(1.0_dp)

contains


  $:TEST("int_re_pls_int_im", label="complex")
    character(*), parameter :: input = "1+1i"
    complex(dp), parameter :: expected = (1, 1)

    complex(dp) :: tokenValue
    integer :: ioStat, iStart

    iStart = 1
    call getNextToken(input, tokenValue, iStart, iostat=ioStat)
    @:ASSERT(ioStat == 0)
    @:ASSERT(iStart == len(input) + 2)
    @:ASSERT(abs(expected - tokenValue) <= tol)
  $:END_TEST()


  $:TEST("int_re_min_int_im", label="complex")
    character(*), parameter :: input = "1-1i"
    complex(dp), parameter :: expected = (1, -1)

    complex(dp) :: tokenValue
    integer :: ioStat, iStart

    iStart = 1
    call getNextToken(input, tokenValue, iStart, iostat=ioStat)
    @:ASSERT(ioStat == 0)
    @:ASSERT(iStart == len(input) + 2)
    @:ASSERT(abs(expected - tokenValue) <= tol)
  $:END_TEST()


  $:TEST("fixed_re_pls_fixed_im", label="complex")
    character(*), parameter :: input = "4.321+5.425i"
    complex(dp), parameter :: expected = (4.321_dp, +5.425_dp)

    complex(dp) :: tokenValue
    integer :: ioStat, iStart

    iStart = 1
    call getNextToken(input, tokenValue, iStart, iostat=ioStat)
    @:ASSERT(ioStat == 0)
    @:ASSERT(iStart == len(input) + 2)
    @:ASSERT(abs(expected - tokenValue) <= tol)
  $:END_TEST()


  $:TEST("fixed_re_min_fixed_im", label="complex")
    character(*), parameter :: input = "4.321-5.425i"
    complex(dp), parameter :: expected = (4.321_dp, -5.425_dp)

    complex(dp) :: tokenValue
    integer :: ioStat, iStart

    iStart = 1
    call getNextToken(input, tokenValue, iStart, iostat=ioStat)
    @:ASSERT(ioStat == 0)
    @:ASSERT(iStart == len(input) + 2)
    @:ASSERT(abs(expected - tokenValue) <= tol)
  $:END_TEST()


  $:TEST("sci_re_pls_sci_im", label="complex")
    character(*), parameter :: input = "4.321e-1+5.425e+2i"
    complex(dp), parameter :: expected = (4.321e-1_dp, 5.425e2_dp)

    complex(dp) :: tokenValue
    integer :: ioStat, iStart

    iStart = 1
    call getNextToken(input, tokenValue, iStart, iostat=ioStat)
    @:ASSERT(ioStat == 0)
    @:ASSERT(iStart == len(input) + 2)
    @:ASSERT(abs(expected - tokenValue) <= tol)
  $:END_TEST()


  $:TEST("sci_re_pls_sci_im_capitalised", label="complex")
    character(*), parameter :: input = "4.321E-1+5.425E+2i"
    complex(dp), parameter :: expected = (4.321e-1_dp, 5.425e2_dp)

    complex(dp) :: tokenValue
    integer :: ioStat, iStart

    iStart = 1
    call getNextToken(input, tokenValue, iStart, iostat=ioStat)
    @:ASSERT(ioStat == 0)
    @:ASSERT(iStart == len(input) + 2)
    @:ASSERT(abs(expected - tokenValue) <= tol)
  $:END_TEST()


  $:TEST("re_only", label="complex")
    character(*), parameter :: input = "4.321e-1"
    complex(dp), parameter :: expected = (4.321e-1_dp, 0.0_dp)

    complex(dp) :: tokenValue
    integer :: ioStat, iStart

    iStart = 1
    call getNextToken(input, tokenValue, iStart, iostat=ioStat)
    @:ASSERT(ioStat == 0)
    @:ASSERT(iStart == len(input) + 2)
    @:ASSERT(abs(expected - tokenValue) <= tol)
  $:END_TEST()


  $:TEST("pos_re_only", label="complex")
    character(*), parameter :: input = "+4.321e-1"
    complex(dp), parameter :: expected = (4.321e-1_dp, 0.0_dp)

    complex(dp) :: tokenValue
    integer :: ioStat, iStart

    iStart = 1
    call getNextToken(input, tokenValue, iStart, iostat=ioStat)
    @:ASSERT(ioStat == 0)
    @:ASSERT(iStart == len(input) + 2)
    @:ASSERT(abs(expected - tokenValue) <= tol)
  $:END_TEST()


  $:TEST("neg_re_only", label="complex")
    character(*), parameter :: input = "-4.321e-1"
    complex(dp), parameter :: expected = (-4.321e-1_dp, 0.0_dp)

    complex(dp) :: tokenValue
    integer :: ioStat, iStart

    iStart = 1
    call getNextToken(input, tokenValue, iStart, iostat=ioStat)
    @:ASSERT(ioStat == 0)
    @:ASSERT(iStart == len(input) + 2)
    @:ASSERT(abs(expected - tokenValue) <= tol)
  $:END_TEST()


  $:TEST("im_only", label="complex")
    character(*), parameter :: input = "4.321e-1i"
    complex(dp), parameter :: expected = (0.0_dp, 4.321e-1_dp)

    complex(dp) :: tokenValue
    integer :: ioStat, iStart

    iStart = 1
    call getNextToken(input, tokenValue, iStart, iostat=ioStat)
    @:ASSERT(ioStat == 0)
    @:ASSERT(iStart == len(input) + 2)
    @:ASSERT(abs(expected - tokenValue) <= tol)
  $:END_TEST()


  $:TEST("fail_jform", label="complex")
    character(*), parameter :: input = "4.321+5.425j"

    complex(dp) :: tokenValue
    integer :: ioStat, iStart

    iStart = 1
    call getNextToken(input, tokenValue, iStart, iostat=ioStat)
    @:ASSERT(ioStat /= 0)
  $:END_TEST()


  $:TEST("fail_with_parentheses", label="complex")
    character(*), parameter :: input = "(3.219-4.321e-1i)"

    complex(dp) :: tokenValue
    integer :: ioStat, iStart

    iStart = 1
    call getNextToken(input, tokenValue, iStart, iostat=ioStat)
    @:ASSERT(ioStat /= 0)
  $:END_TEST()


  $:TEST("fail_plus_minus", label="complex")
    character(*), parameter :: input = "1.235e-1+-4.321e+2i"

    complex(dp) :: tokenValue
    integer :: ioStat, iStart

    iStart = 1
    call getNextToken(input, tokenValue, iStart, iostat=ioStat)
    @:ASSERT(ioStat /= 0)
  $:END_TEST()


  $:TEST("fail_missing_exp_re", label="complex")
    character(*), parameter :: input = "1.235-1+4.321e+2i"

    complex(dp) :: tokenValue
    integer :: ioStat, iStart

    iStart = 1
    call getNextToken(input, tokenValue, iStart, iostat=ioStat)
    @:ASSERT(ioStat /= 0)
  $:END_TEST()


  $:TEST("fail_missing_exp_im", label="complex")
    character(*), parameter :: input = "1.235e-1+4.321+2i"

    complex(dp) :: tokenValue
    integer :: ioStat, iStart

    iStart = 1
    call getNextToken(input, tokenValue, iStart, iostat=ioStat)
    @:ASSERT(ioStat /= 0)
  $:END_TEST()


  function tests()
    type(test_list) :: tests

    tests = test_list([&
        suite("tokenreader", test_list([&
            $:TEST_ITEMS(label="complex")
        ]))&
    ])
    $:STOP_ON_MISSING_TEST_ITEMS()

  end function tests

end module test_io_tokenreader
