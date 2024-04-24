!--------------------------------------------------------------------------------------------------!
!  DFTB+: general package for performing fast atomistic simulations                                !
!  Copyright (C) 2006 - 2024  DFTB+ developers group                                               !
!                                                                                                  !
!  See the LICENSE file for terms of usage and distribution.                                       !
!--------------------------------------------------------------------------------------------------!

#:include 'common.fypp'

!> Plugin loading mechanism
module dftbp_plugins_plugin
  use, intrinsic :: iso_c_binding, only : c_associated, c_char, c_double, c_int, c_ptr
  use dftbp_common_accuracy, only : dp
  use dftbp_io_message, only : error
  implicit none

  private
  public :: TPlugin

  !> Keep the handle to the plugin
  type, bind(C) :: c_handle
    type(c_ptr) :: cptr
  end type

  !> Type to manage a plugin
  type TPlugin

    !> Handle to the plugin
    type(c_handle), private :: handle

    !> Whether the plugin provides certain features
    logical :: provides_getSKIntegrals = .false.

    !> Whether the plugin is initialized
    logical :: initialized = .false.

  contains

    procedure :: init => TPlugin_init
    procedure :: getSKIntegrals => TPlugin_getSKIntegrals
    final :: TPlugin_final

  end type TPlugin

  interface

    !> Load and initialize the plugin
    function init_plugin_c(filename) result(handle) bind(C, name='init_plugin')
      import c_char, c_handle
      character(kind=c_char), intent(in) :: filename(*)
      type(c_handle) :: handle
    end function init_plugin_c

    !> Close the plugin
    subroutine final_plugin_c(handle) bind(C, name='final_plugin')
      import c_handle
      type(c_handle), value :: handle
    end subroutine final_plugin_c

    !> Check if the plugin implements a certain function
    function provides_plugin_c(handle, func) result(success) bind(C, name='provides_plugin')
      import c_handle, c_char, c_int
      type(c_handle), value, intent(in) :: handle
      character(kind=c_char), intent(in) :: func(*)
      integer(c_int) :: success
    end function provides_plugin_c

    !> Call the implemented function
    function call_getSKIntegrals_c(handle, sk, dist, sp1, sp2) result(success)&
        & bind(C, name='call_getSKIntegrals')
      import c_handle, c_double, c_int
      type(c_handle), value, intent(in) :: handle
      real(c_double), intent(out) :: sk(:)
      real(c_double), intent(in) :: dist
      integer(c_int), intent(in) :: sp1, sp2
      integer(c_int) :: success
    end function call_getSKIntegrals_c

  end interface

contains

  !> Sets the handle to the plugin library
  function TPlugin_init(this, filename) result(success)

    !> Instance
    class(TPlugin), intent(inout) :: this

    !> Plugin file name
    character(len=*), intent(in) :: filename

    !> Success of initialization
    logical :: success

    this%handle = init_plugin_c(trim(filename) // char(0))
    this%initialized = c_associated(this%handle%cptr)

    if (this%initialized) then
      this%provides_getSKIntegrals = provides_plugin_c(this%handle, "getSKIntegrals" // char(0))&
          & == 1
    end if

    success = this%initialized

  end function TPlugin_init

  !> Frees all references to the plugin
  subroutine TPlugin_final(this)

    !> Instance
    type(TPlugin), intent(inout) :: this

    if (this%initialized) then
      call final_plugin_c(this%handle)
      this%initialized = .false.
    end if

  end subroutine TPlugin_final

  !> Returns the Slater-Koster integrals for a given distance for a given atom pair
  subroutine TPlugin_getSKIntegrals(this, sk, dist, sp1, sp2)

    !> Instance
    class(TPlugin), intent(in) :: this

    !> Contains the integrals on exit
    real(dp), intent(out) :: sk(:)

    !> Distance of the two atoms
    real(dp), intent(in) :: dist

    !> Index of the first interacting species
    integer, intent(in) :: sp1

    !> Index of the second interacting species
    integer, intent(in) :: sp2

    integer :: success

    if (.not. this%initialized) then
      call error("Trying to call a function in an uninitialized plugin")
    end if
    if (.not. this%provides_getSKIntegrals) then
      call error("Trying to call a function not provided by the plugin")
    end if

    success = call_getSKIntegrals_c(this%handle, sk, dist, sp1, sp2)

  end subroutine TPlugin_getSKIntegrals

end module dftbp_plugins_plugin
