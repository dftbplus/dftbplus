!--------------------------------------------------------------------------------------------------!
!  DFTB+: general package for performing fast atomistic simulations                                !
!  Copyright (C) 2006 - 2020  DFTB+ developers group                                               !
!                                                                                                  !
!  See the LICENSE file for terms of usage and distribution.                                       !
!--------------------------------------------------------------------------------------------------!

module dftbp_elecsolvertypes
  implicit none
  private

  public :: electronicSolverTypes


  !> Namespace for possible solver methods
  type :: TElecSolverTypesEnum

    ! lapack/scalapack solvers
    integer :: qr = 1
    integer :: divideandconquer = 2
    integer :: relativelyrobust = 3

    ! elsi provided solvers
    integer :: elpa = 4
    integer :: omm = 5
    integer :: pexsi = 6
    integer :: dummy1 = 7
    integer :: dummy2 = 8
    integer :: ntpoly = 9
    integer :: elpadm = 10

    ! transport related
    integer :: gf = 11
    integer :: onlyTransport = 12

    ! GPU accellerated solvers using MAGMA
    integer :: magma_gvd = 13

  end type TElecSolverTypesEnum


  !> Actual values for elecSolverTypes.
  type(TElecSolverTypesEnum), parameter :: electronicSolverTypes = TElecSolverTypesEnum()

end module dftbp_elecsolvertypes
