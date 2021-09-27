!--------------------------------------------------------------------------------------------------!
!  DFTB+: general package for performing fast atomistic simulations                                !
!  Copyright (C) 2006 - 2021  DFTB+ developers group                                               !
!                                                                                                  !
!  See the LICENSE file for terms of usage and distribution.                                       !
!--------------------------------------------------------------------------------------------------!

!> Input for geometry optimizations
module dftbp_geoopt_input
  use dftbp_common_accuracy, only : dp
  use dftbp_geoopt_filter, only : TFilterInput
  use dftbp_geoopt_fire, only : TFireInput
  use dftbp_geoopt_lbfgs2, only : TLBFGSInput
  use dftbp_geoopt_rf, only : TRationalFunctionInput
  implicit none
  private
  public :: TGeoOptInput, TOptConv
  public :: TFilterInput, TRationalFunctionInput, TLBFGSInput, TFireInput


  !> Thresholds for optimization
  type :: TOptConv

    !> Convergence threshold for energy
    real(dp) :: ethr = huge(1.0_dp)

    !> Convergence threshold for gradient norm
    real(dp) :: gthr = huge(1.0_dp)

    !> Convergence threshold for gradient norm
    real(dp) :: gmax = huge(1.0_dp)

    !> Convergence threshold for displacement norm
    real(dp) :: dthr = huge(1.0_dp)

    !> Convergence threshold for displacement norm
    real(dp) :: dmax = huge(1.0_dp)

  end type TOptConv


  !> General input wrapper for geometry optimization
  type :: TGeoOptInput

    !> Input for coordinate transformation and filter step
    type(TFilterInput) :: filter

    !> Input for rational function optimization driver
    type(TRationalFunctionInput), allocatable :: rf

    !> Input for limited memory BFGS optimizer
    type(TLBFGSInput), allocatable :: lbfgs

    !> Input for fast inertial relaxation engine
    type(TFireInput), allocatable :: fire

    !> Number of allowed geometry optimization steps
    integer :: nGeoSteps = huge(1) - 1

    !> Prefix of the output file name
    character(len=:), allocatable :: outFile

    !> Accuracy to determine convergence thresholds
    real(dp) :: accuracy = 1.0_dp

    !> Threshold for optimization
    type(TOptConv) :: conv

  end type TGeoOptInput


end module dftbp_geoopt_input
