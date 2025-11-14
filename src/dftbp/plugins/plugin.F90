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

  integer, parameter :: major = 1, minor = 0

  private
  public :: TPlugin

  !> Keep the handle to the plugin
  type, bind(C) :: c_handle
    type(c_ptr) :: cptr
  end type

  !> Type to manage plugin capabilities, i.e. what function the plugin actually implements
  type TPluginCapabilities

    !> Whether the plugin updates SK data
    logical :: provides_updateSKIntegrals = .false.

    !> Whether the plugin needs the neighbour list
    logical :: provides_readNeighbourList = .false.

    !> Whether the plugin needs the atomic self energy
    logical :: provides_readAtomSelfEnergy = .false.

    !> Whether the plugin needs the Hubbard Us
    logical :: provides_readHubbardU = .false.

  end type TPluginCapabilities

  !> Type to manage a plugin
  type TPlugin

    !> Handle to the plugin
    type(c_handle), private :: handle

    !> Plugin capabiltites
    type(TPluginCapabilities) :: capabilities

    !> Whether the plugin is initialized
    logical :: initialized = .false.

  contains

    procedure :: init => TPlugin_init
    procedure :: updateSKIntegrals => TPlugin_updateSKIntegrals
    procedure :: readNeighbourList => TPlugin_readNeighbourList
    procedure :: readAtomSelfEnergy => TPlugin_readAtomSelfEnergy
    procedure :: readHubbardU => TPlugin_readHubbardU
    final :: TPlugin_final

  end type TPlugin

  !> Type bound to C for fetching the plugin capabilities
  type, bind(C) :: TPluginCapabilities_c
    integer(c_int) :: provides_updateSKIntegrals = 0
    integer(c_int) :: provides_readNeighbourList = 0
    integer(c_int) :: provides_readAtomSelfEnergy = 0
    integer(c_int) :: provides_readHubbardU = 0
  end type TPluginCapabilities_c

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

    !> Send the current API version to the plugin and retrieve a success state
    function version_plugin_c(handle, major, minor) result(success) bind(C, name='version_plugin')
      import c_handle, c_int
      type(c_handle), value, intent(in) :: handle
      integer(c_int), value, intent(in) :: major, minor
      integer(c_int) :: success
    end function version_plugin_c

    !> Check if the plugin implements a certain function
    function provides_plugin_c(handle, capabilities) result(success) bind(C, name='provides_plugin')
      import TPluginCapabilities_c
      import c_handle, c_int
      type(c_handle), value, intent(in) :: handle
      type(TPluginCapabilities_c), intent(inout) :: capabilities
      integer(c_int) :: success
    end function provides_plugin_c

    !> Call the implemented function for SK integrals
    function call_updateSKIntegrals_c(handle, nSkgrid, nSkIntg, skTab, dist, atom1, atom2,&
          & species1, species2, HorS, interdist) result(updated)&
          & bind(C, name='call_updateSKIntegrals')
      import c_handle, c_double, c_int
      type(c_handle), value, intent(in) :: handle
      integer(c_int), value, intent(in) :: nSkgrid, nSkIntg
      real(c_double), intent(in) :: skTab(*)
      real(c_double), value, intent(in) :: dist
      integer(c_int), value, intent(in) :: atom1, atom2, species1, species2, HorS
      real(c_double), value, intent(in) :: interdist
      integer(c_int) :: updated
    end function call_updateSKIntegrals_c

    !> Call the implemented function with the neighbour list
    subroutine call_readNeighbourList_c(handle, nAtoms, nNeighbours, nAtomsCent, coords,&
          & img2CentCell, iNeighbour, neighDist2) bind(C, name='call_readNeighbourList')
      import c_handle, c_double, c_int
      type(c_handle), value, intent(in) :: handle
      integer(c_int), value, intent(in) :: nAtoms, nNeighbours, nAtomsCent
      real(c_double), intent(in) :: coords(*)
      integer(c_int), intent(in) :: img2CentCell(*), iNeighbour(*)
      real(c_double), intent(in) :: neighDist2(*)
    end subroutine call_readNeighbourList_c

    !> Call the implemented function with the atomic self energy
    subroutine call_readAtomSelfEnergy_c(handle, nOrbitals, nAtoms, atomEigVal)&
          & bind(C, name='call_readAtomSelfEnergy')
      import c_handle, c_double, c_int
      type(c_handle), value, intent(in) :: handle
      integer(c_int), value, intent(in) :: nOrbitals, nAtoms
      real(c_double), intent(in) :: atomEigVal(*)
    end subroutine call_readAtomSelfEnergy_c

    !> Call the implemented function with the Hubbard Us
    subroutine call_readHubbardU_c(handle, nShells, nSpecies, nHubbU, uniqHubbU)&
          & bind(C, name='call_readHubbardU')
      import c_handle, c_double, c_int
      type(c_handle), value, intent(in) :: handle
      integer(c_int), value, intent(in) :: nShells, nSpecies
      integer(c_int), intent(in) :: nHubbU(*)
      real(c_double), intent(in) :: uniqHubbU(*)
    end subroutine call_readHubbardU_c

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

    type(TPluginCapabilities_c) :: capabilities_c

    this%handle = init_plugin_c(trim(filename) // char(0))
    this%initialized = c_associated(this%handle%cptr)

    if (this%initialized) then
      if (version_plugin_c(this%handle, major, minor) /= 1) then
        call error("The plugin does not support the current API version.")
        this%initialized = .false.
      end if
    end if

    if (this%initialized) then
      if (provides_plugin_c(this%handle, capabilities_c) /= 1) then
        call error("Unable to determine the capabilities of the plugin.")
        this%initialized = .false.
      else
        this%capabilities%provides_updateSKIntegrals =&
            & capabilities_c%provides_updateSKIntegrals == 1
        this%capabilities%provides_readNeighbourList =&
            & capabilities_c%provides_readNeighbourList == 1
        this%capabilities%provides_readAtomSelfEnergy =&
            & capabilities_c%provides_readAtomSelfEnergy == 1
        this%capabilities%provides_readHubbardU =&
            & capabilities_c%provides_readHubbardU == 1
      end if
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

  !> Updates the Slater-Koster integrals for a given distance for a given atom pair
  function TPlugin_updateSKIntegrals(this, skTab, dist, atom1, atom2, species1, species2, isH,&
        & interdist) result(updated)

    !> Instance
    class(TPlugin), intent(in) :: this

    !> Contains the integrals on exit
    real(dp), intent(inout) :: skTab(:,:)

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

    logical :: updated
    integer :: HorS

    if (.not. this%initialized) then
      call error("Trying to call a function in an uninitialized plugin")
    end if
    if (.not. this%capabilities%provides_updateSKIntegrals) then
      call error("Trying to call a function not provided by the plugin")
    end if

    HorS = 0
    if (.not. isH) then
      HorS = 1
    end if
    updated = call_updateSKIntegrals_c(this%handle, size(skTab,1), size(skTab,2), skTab, dist,&
        & atom1, atom2, species1, species2, HorS, interdist) == 1

  end function TPlugin_updateSKIntegrals

  !> Sets the neighbour list
  subroutine TPlugin_readNeighbourList(this, coords, img2CentCell, iNeighbour, neighDist2)

    !> Instance
    class(TPlugin), intent(in) :: this

    !> Atom coordinates
    real(dp), intent(in) :: coords(:,:)

    !> Mapping of atom number to central cell atom number
    integer, intent(in) :: img2CentCell(:)

    !> index of neighbour atoms
    integer, intent(in) :: iNeighbour(:,:)

    !> neighbour distances
    real(dp), intent(in) :: neighDist2(:,:)

    if (.not. this%initialized) then
      call error("Trying to call a function in an uninitialized plugin")
    end if
    if (.not. this%capabilities%provides_readNeighbourList) then
      call error("Trying to call a function not provided by the plugin")
    end if

    call call_readNeighbourList_c(this%handle, size(img2CentCell, dim=1), size(iNeighbour, dim=1),&
        & size(iNeighbour, dim=2), coords, img2CentCell, iNeighbour, neighDist2)

  end subroutine TPlugin_readNeighbourList

  !> Sets the atomic self energy
  subroutine TPlugin_readAtomSelfEnergy(this, atomEigVal)

    !> Instance
    class(TPlugin), intent(in) :: this

    !> Self energy (orbital, atom)
    real(dp), intent(inout) :: atomEigVal(:,:)

    if (.not. this%initialized) then
      call error("Trying to call a function in an uninitialized plugin")
    end if
    if (.not. this%capabilities%provides_readAtomSelfEnergy) then
      call error("Trying to call a function not provided by the plugin")
    end if

    call call_readAtomSelfEnergy_c(this%handle, size(atomEigVal, dim=1), size(atomEigVal, dim=2),&
        & atomEigVal)

  end subroutine TPlugin_readAtomSelfEnergy

  !> Sets the Hubbard Us
  subroutine TPlugin_readHubbardU(this, nHubbU, uniqHubbU, iHubbU)

    !> Instance
    class(TPlugin), intent(in) :: this

    !> Nr. of uniq Us per species. Shape: [nSpecies]
    integer, intent(in) :: nHubbU(:)

    !> Uniq Us per species. Shape: [mShell, nSpecies]
    real(dp), intent(inout) :: uniqHubbU(:,:)

    !> Mapping L-shell to unique Hubbard U. Shape: [mShell, nSpecies]
    integer, intent(in) :: iHubbU(:, :)

    if (.not. this%initialized) then
      call error("Trying to call a function in an uninitialized plugin")
    end if
    if (.not. this%capabilities%provides_readHubbardU) then
      call error("Trying to call a function not provided by the plugin")
    end if

    call call_readHubbardU_c(this%handle, size(uniqHubbU, dim=1), size(uniqHubbU, dim=2),&
        & nHubbU, uniqHubbU)

  end subroutine TPlugin_readHubbardU

end module dftbp_plugins_plugin
