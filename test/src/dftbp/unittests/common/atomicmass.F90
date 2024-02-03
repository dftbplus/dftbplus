!--------------------------------------------------------------------------------------------------!
!  DFTB+: general package for performing fast atomistic simulations                                !
!  Copyright (C) 2006 - 2024  DFTB+ developers group                                               !
!                                                                                                  !
!  See the LICENSE file for terms of usage and distribution.                                       !
!--------------------------------------------------------------------------------------------------!

#:include "fytest.fypp"

#:block TEST_SUITE("atomicmass")
  use dftbp_common_accuracy, only : dp
  use dftbp_common_atomicmass, only : getAtomicSymbol
  implicit none

#:contains

  #:block TEST_FIXTURE("getAtomicSymbol")

  #:contains

    #:block TEST("testSomeElements")
      @:ASSERT(getAtomicSymbol(1.0_dp) == 'H')
      @:ASSERT(getAtomicSymbol(12.0_dp) == 'C')
      @:ASSERT(getAtomicSymbol(16.0_dp) == 'O')
    #:endblock

    #:block TEST("testTooLargeDifference")
      @:ASSERT(getAtomicSymbol(500.0_dp) == '??')
    #:endblock

    #:block TEST("testNegativeMass")
      @:ASSERT(getAtomicSymbol(-1.0_dp) == '??')
    #:endblock

  #:endblock TEST_FIXTURE

#:endblock TEST_SUITE


@:TEST_DRIVER()
