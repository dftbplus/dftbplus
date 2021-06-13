!--------------------------------------------------------------------------------------------------!
!  DFTB+: general package for performing fast atomistic simulations                                !
!  Copyright (C) 2006 - 2021  DFTB+ developers group                                               !
!                                                                                                  !
!  See the LICENSE file for terms of usage and distribution.                                       !
!--------------------------------------------------------------------------------------------------!

#:include 'common.fypp'

!> Contains data type representing the input data for setupgeom
module dftbp_inputsetup
  use dftbp_common_assert
  use dftbp_common_accuracy
  use dftbp_type_typegeometry
  use dftbp_io_message
  use dftbp_dftb_slakocont
  use dftbp_type_commontypes
  use dftbp_dftb_repcont
  use dftbp_type_linkedlist
  use dftbp_type_wrappedintr
  use dftbp_transport_negfvars

  implicit none
  private

  public :: TGeometry, TSlater, TInputData
  public :: init, destruct

  !> Slater-Koster data
  type TSlater
    real(dp), allocatable :: skSelf(:, :)
    real(dp), allocatable :: skHubbU(:, :)
    real(dp), allocatable :: skOcc(:, :)
    real(dp), allocatable :: mass(:)

    type(TSlakoCont), allocatable :: skHamCont
    type(TSlakoCont), allocatable :: skOverCont
    type(TRepCont), allocatable :: repCont
    type(TOrbitals), allocatable :: orb
  end type TSlater

  !> container for input data constituents
  type TInputData
    logical :: tInitialized = .false.
    type(TGeometry) :: geom
    type(TSlater) :: slako
    type(TTransPar) :: transpar
  end type TInputData


  !> Initialise the input data
  interface init
    module procedure InputData_init
  end interface init


  !> destroy input data for variables that do not go out of scope
  interface destruct
    module procedure InputData_destruct
  end interface destruct

contains


  !> Mark data structure as initialised
  subroutine InputData_init(this)
    type(TInputData), intent(out) :: this

    this%tInitialized = .true.

  end subroutine InputData_init


  !> destructor for parts that are not cleaned up when going out of scope
  subroutine InputData_destruct(this)
    type(TInputData), intent(inout) :: this

  end subroutine InputData_destruct

end module dftbp_inputsetup
