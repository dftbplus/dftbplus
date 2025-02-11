!--------------------------------------------------------------------------------------------------!
!  DFTB+: general package for performing fast atomistic simulations                                !
!  Copyright (C) 2006 - 2025  DFTB+ developers group                                               !
!                                                                                                  !
!  See the LICENSE file for terms of usage and distribution.                                       !
!--------------------------------------------------------------------------------------------------!

#:include "fytest.fypp"

#:block TEST_SUITE("tokenreader")
  use dftbp_common_accuracy, only : dp
  use dftbp_io_tokenreader, only : getNextToken
  implicit none

  integer, parameter :: tol = 10.0_dp * epsilon(1.0_dp)

#:contains

  #:block TEST_FIXTURE("complex")
    complex(dp) :: tokenValue
    integer :: ioStat, iStart

  #:contains

    #:block TEST("int_re_pls_int_im")
      character(*), parameter :: input = "1+1i"
      complex(dp), parameter :: expected = (1, 1)
      iStart = 1
      call getNextToken(input, tokenValue, iStart, iostat=ioStat)
      @:ASSERT(ioStat == 0)
      @:ASSERT(iStart == len(input) + 2)
      @:ASSERT(abs(expected - tokenValue) <= tol)
    #:endblock

    #:block TEST("int_re_min_int_im")
      character(*), parameter :: input = "1-1i"
      complex(dp), parameter :: expected = (1, -1)
      iStart = 1
      call getNextToken(input, tokenValue, iStart, iostat=ioStat)
      @:ASSERT(ioStat == 0)
      @:ASSERT(iStart == len(input) + 2)
      @:ASSERT(abs(expected - tokenValue) <= tol)
    #:endblock

    #:block TEST("fixed_re_pls_fixed_im")
      character(*), parameter :: input = "4.321+5.425i"
      complex(dp), parameter :: expected = (4.321_dp, +5.425_dp)
      iStart = 1
      call getNextToken(input, tokenValue, iStart, iostat=ioStat)
      @:ASSERT(ioStat == 0)
      @:ASSERT(iStart == len(input) + 2)
      @:ASSERT(abs(expected - tokenValue) <= tol)
    #:endblock

    #:block TEST("fixed_re_min_fixed_im")
      character(*), parameter :: input = "4.321-5.425i"
      complex(dp), parameter :: expected = (4.321_dp, -5.425_dp)
      iStart = 1
      call getNextToken(input, tokenValue, iStart, iostat=ioStat)
      @:ASSERT(ioStat == 0)
      @:ASSERT(iStart == len(input) + 2)
      @:ASSERT(abs(expected - tokenValue) <= tol)
    #:endblock

    #:block TEST("sci_re_pls_sci_im")
      character(*), parameter :: input = "4.321e-1+5.425e+2i"
      complex(dp), parameter :: expected = (4.321e-1_dp, 5.425e2_dp)
      iStart = 1
      call getNextToken(input, tokenValue, iStart, iostat=ioStat)
      @:ASSERT(ioStat == 0)
      @:ASSERT(iStart == len(input) + 2)
      @:ASSERT(abs(expected - tokenValue) <= tol)
    #:endblock

    #:block TEST("sci_re_pls_sci_im_capitalised")
      character(*), parameter :: input = "4.321E-1+5.425E+2i"
      complex(dp), parameter :: expected = (4.321e-1_dp, 5.425e2_dp)
      iStart = 1
      call getNextToken(input, tokenValue, iStart, iostat=ioStat)
      @:ASSERT(ioStat == 0)
      @:ASSERT(iStart == len(input) + 2)
      @:ASSERT(abs(expected - tokenValue) <= tol)
    #:endblock

    #:block TEST("re_only")
      character(*), parameter :: input = "4.321e-1"
      complex(dp), parameter :: expected = (4.321e-1_dp, 0.0_dp)
      iStart = 1
      call getNextToken(input, tokenValue, iStart, iostat=ioStat)
      @:ASSERT(ioStat == 0)
      @:ASSERT(iStart == len(input) + 2)
      @:ASSERT(abs(expected - tokenValue) <= tol)
    #:endblock

    #:block TEST("pos_re_only")
      character(*), parameter :: input = "+4.321e-1"
      complex(dp), parameter :: expected = (4.321e-1_dp, 0.0_dp)
      iStart = 1
      call getNextToken(input, tokenValue, iStart, iostat=ioStat)
      @:ASSERT(ioStat == 0)
      @:ASSERT(iStart == len(input) + 2)
      @:ASSERT(abs(expected - tokenValue) <= tol)
    #:endblock

    #:block TEST("neg_re_only")
      character(*), parameter :: input = "-4.321e-1"
      complex(dp), parameter :: expected = (-4.321e-1_dp, 0.0_dp)
      iStart = 1
      call getNextToken(input, tokenValue, iStart, iostat=ioStat)
      @:ASSERT(ioStat == 0)
      @:ASSERT(iStart == len(input) + 2)
      @:ASSERT(abs(expected - tokenValue) <= tol)
    #:endblock

    #:block TEST("im_only")
      character(*), parameter :: input = "4.321e-1i"
      complex(dp), parameter :: expected = (0.0_dp, 4.321e-1_dp)
      iStart = 1
      call getNextToken(input, tokenValue, iStart, iostat=ioStat)
      @:ASSERT(ioStat == 0)
      @:ASSERT(iStart == len(input) + 2)
      @:ASSERT(abs(expected - tokenValue) <= tol)
    #:endblock

    #:block TEST("fail_jform")
      character(*), parameter :: input = "4.321+5.425j"
      iStart = 1
      call getNextToken(input, tokenValue, iStart, iostat=ioStat)
      @:ASSERT(ioStat /= 0)
    #:endblock

    #:block TEST("fail_with_parentheses")
      character(*), parameter :: input = "(3.219-4.321e-1i)"
      iStart = 1
      call getNextToken(input, tokenValue, iStart, iostat=ioStat)
      @:ASSERT(ioStat /= 0)
    #:endblock

    #:block TEST("fail_plus_minus")
      character(*), parameter :: input = "1.235e-1+-4.321e+2i"
      iStart = 1
      call getNextToken(input, tokenValue, iStart, iostat=ioStat)
      @:ASSERT(ioStat /= 0)
    #:endblock

    #:block TEST("fail_missing_exp_re")
      character(*), parameter :: input = "1.235-1+4.321e+2i"
      iStart = 1
      call getNextToken(input, tokenValue, iStart, iostat=ioStat)
      @:ASSERT(ioStat /= 0)
    #:endblock

    #:block TEST("fail_missing_exp_im")
      character(*), parameter :: input = "1.235e-1+4.321+2i"
      iStart = 1
      call getNextToken(input, tokenValue, iStart, iostat=ioStat)
      @:ASSERT(ioStat /= 0)
    #:endblock

  #:endblock TEST_FIXTURE

#:endblock TEST_SUITE


@:TEST_DRIVER()
