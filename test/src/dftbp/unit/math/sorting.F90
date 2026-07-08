!--------------------------------------------------------------------------------------------------!
!  DFTB+: general package for performing fast atomistic simulations                                !
!  Copyright (C) 2006 - 2025  DFTB+ developers group                                               !
!                                                                                                  !
!  See the LICENSE file for terms of usage and distribution.                                       !
!--------------------------------------------------------------------------------------------------!

#:include "fortuno_serial.fypp"
module test_math_sorting
  use dftbp_common_accuracy, only : dp, rsp
  use dftbp_math_sorting, only : merge_sort
  use fortuno_serial, only : suite => serial_suite_item, test_list, all_close
  $:FORTUNO_SERIAL_IMPORTS()
  implicit none

  private
  public :: tests

contains


  $:TEST("mergesort_real")
    integer :: nn
    integer, allocatable :: indx(:)
    integer :: ii
    logical :: isSorted, isStable
    real(dp), parameter :: tol = real(epsilon(0.0_rsp),dp)

    real(dp), parameter :: data(42) = [& ! Data from C60_OscWindow test
        & 1.0009256850962192E-1_dp, 1.0009256850980475E-1_dp, 1.0009256850996931E-1_dp,&
        & 1.0009256824902155E-1_dp, 1.0009256824920437E-1_dp, 1.0009256824936894E-1_dp,&
        & 1.0009256824896715E-1_dp, 1.0009256824914997E-1_dp, 1.0009256824931453E-1_dp,&
        & 1.0009256824891125E-1_dp, 1.0009256824909407E-1_dp, 1.0009256824925863E-1_dp,&
        & 6.6042780967860620E-2_dp, 6.6042780968043446E-2_dp, 6.6042780968208009E-2_dp,&
        & 9.4447988867716701E-2_dp, 9.4447988867874991E-2_dp, 9.4447988868041677E-2_dp,&
        & 6.6042780967824899E-2_dp, 6.6042780968007725E-2_dp, 6.6042780968172288E-2_dp,&
        & 9.4447988867680979E-2_dp, 9.4447988867839269E-2_dp, 9.4447988868005955E-2_dp,&
        & 6.6042780967797060E-2_dp, 6.6042780967979886E-2_dp, 6.6042780968144449E-2_dp,&
        & 9.4447988867653140E-2_dp, 9.4447988867811430E-2_dp, 9.4447988867978117E-2_dp,&
        & 6.6042780825941200E-2_dp, 6.6042780826124026E-2_dp, 6.6042780826288588E-2_dp,&
        & 9.4447988725797280E-2_dp, 9.4447988725955570E-2_dp, 9.4447988726122256E-2_dp,&
        & 6.6042780825940506E-2_dp, 6.6042780826123332E-2_dp, 6.6042780826287895E-2_dp,&
        & 9.4447988725796586E-2_dp, 9.4447988725954876E-2_dp, 9.4447988726121562E-2_dp ]

    nn = size(data)

    allocate(indx(nn), source=0)
    call merge_sort(indx, data, tol)

    isSorted = .true.
    isStable = .true.
    lpCheck: do ii = 2, nn

      if (data(indx(ii-1)) - data(indx(ii)) >= tol) then
        write(*,*)
        write(*,*)'Data is not sorted'
        isSorted = .false.
        exit lpCheck
      else if (abs(data(indx(ii)) - data(indx(ii-1))) < tol .and. indx(ii) < indx(ii-1)) then
        write(*,*)
        write(*,*)'Sort is not stable'
        isStable = .false.
        exit lpCheck
      end if

    end do lpCheck

    if (.not.(isSorted .or. isStable)) then
      do ii = 1, nn
        write(*,*)indx(ii), data(indx(ii))
      end do
    end if

    @:ASSERT(isSorted)
    @:ASSERT(isStable)

  $:END_TEST()


  function tests()
    type(test_list) :: tests

    tests = test_list([&
        suite("sorting", test_list([&
            $:TEST_ITEMS()
        ]))&
    ])
    $:STOP_ON_MISSING_TEST_ITEMS()

  end function tests

end module test_math_sorting
