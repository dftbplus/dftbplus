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

    !> Whether the plugin provides SK data
    logical :: provides_getSKIntegrals = .false.

    !> Whether the plugin needs the neighbour list
    logical :: provides_setNeighbourList = .false.

    !> Whether the plugin is initialized
    logical :: initialized = .false.

  contains

    procedure :: init => TPlugin_init
    procedure :: getSKIntegrals => TPlugin_getSKIntegrals
    procedure :: setNeighbourList => TPlugin_setNeighbourList
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

    !> Call the implemented function for SK integrals
    function call_getSKIntegrals_c(handle, nSkgrid, nSkIntg, skTab, dist, atom1, atom2, species1,&
          & species2, HorS, interdist) result(success) bind(C, name='call_getSKIntegrals')
      import c_handle, c_double, c_int
      type(c_handle), value, intent(in) :: handle
      integer(c_int), value, intent(in) :: nSkgrid, nSkIntg
      real(c_double), intent(in) :: skTab(*)
      real(c_double), value, intent(in) :: dist
      integer(c_int), value, intent(in) :: atom1, atom2, species1, species2, HorS
      real(c_double), value, intent(in) :: interdist
      integer(c_int) :: success
    end function call_getSKIntegrals_c

    !> Call the implemented function for setting the neighbour list
    subroutine call_setNeighbourList_c(handle, nAtoms, nAtomsCent, coords, img2CentCell,&
          & iNeighbour, neightDist2) bind(C, name='call_setNeighbourList')
      import c_handle, c_double, c_int
      type(c_handle), value, intent(in) :: handle
      integer(c_int), value, intent(in) :: nAtoms, nAtomsCent
      real(c_double), intent(in) :: coords(*)
      integer(c_int), intent(in) :: img2CentCell(*), iNeighbour(*)
      real(c_double), intent(in) :: neightDist2(*)
    end subroutine call_setNeighbourList_c

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
      this%provides_getSKIntegrals = provides_plugin_c(this%handle, "getSKIntegrals"//char(0))&
          & == 1
      this%provides_setNeighbourList = provides_plugin_c(this%handle, "setNeighbourList"//char(0))&
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
  function TPlugin_getSKIntegrals(this, skTab, dist, atom1, atom2, species1, species2, isH,&
        & interdist) result(success)

    !> Instance
    class(TPlugin), intent(in) :: this

    !> Contains the integrals on exit
    real(dp), intent(in) :: skTab(:,:)

    !> Distance of the two atoms
    real(dp), intent(in) :: dist

    !> Index of the first atom
    integer, intent(in) :: atom1

    !> Index of the second atom
    integer, intent(in) :: atom2

    !> Index of the first interacting species
    integer, intent(in) :: species1

    !> Index of the second interacting species
    integer, intent(in) :: species2

    !> Specified the container for Hamitonian (==.true.) or Overlap (==.false.).
    logical, intent(in) :: isH

    !> Distance between the two interpolation in skTab.
    real(dp), intent(in) :: interdist

    logical :: success
    integer :: HorS

    if (.not. this%initialized) then
      call error("Trying to call a function in an uninitialized plugin")
    end if
    if (.not. this%provides_getSKIntegrals) then
      call error("Trying to call a function not provided by the plugin")
    end if

    HorS = 0
    if (.not. isH) then
      HorS = 1
    end if
    success = call_getSKIntegrals_c(this%handle, size(skTab,1), size(skTab,2), skTab, dist, atom1,&
        & atom2, species1, species2, HorS, interdist) == 1

  end function TPlugin_getSKIntegrals

  !> Sets the neighbour list
  subroutine TPlugin_setNeighbourList(this, coords, img2CentCell, iNeighbour, neightDist2)

    !> Instance
    class(TPlugin), intent(in) :: this

    !> Atom coordinates
    real(dp), intent(in) :: coords(:,:)

    !> Mapping of atom number to central cell atom number
    integer, intent(in) :: img2CentCell(:)

    !> index of neighbour atoms
    integer, intent(in) :: iNeighbour(:,:)

    !> neighbour distances
    real(dp), intent(in) :: neightDist2(:,:)

    if (.not. this%initialized) then
      call error("Trying to call a function in an uninitialized plugin")
    end if
    if (.not. this%provides_setNeighbourList) then
      call error("Trying to call a function not provided by the plugin")
    end if

    call call_setNeighbourList_c(this%handle, size(img2CentCell, dim=1), size(iNeighbour, dim=2),&
        & coords, img2CentCell, iNeighbour, neightDist2)

  end subroutine TPlugin_setNeighbourList

end module dftbp_plugins_plugin
