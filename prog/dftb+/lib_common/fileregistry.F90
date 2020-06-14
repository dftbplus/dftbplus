!--------------------------------------------------------------------------------------------------!
!  DFTB+: general package for performing fast atomistic simulations                                !
!  Copyright (C) 2006 - 2020  DFTB+ developers group                                               !
!                                                                                                  !
!  See the LICENSE file for terms of usage and distribution.                                       !
!--------------------------------------------------------------------------------------------------!

!> Contains a file registry which keeps track of eventually open files and closes them.
module dftbp_fileregistry
  implicit none
  private

  public :: TFileRegistry_init, TFileRegistry


  !> Implements a file registry
  type :: TFileRegistry
    private
    integer, allocatable :: registeredFiles(:)
    integer :: nRegistered = -1
  contains
    procedure :: register => TFileRegistry_register
    procedure :: closeAll => TFileRegistry_closeAll
  end type TFileRegistry


  !> Default registry size (nr. of files to be stored before reallocation)
  integer, parameter :: defaultSize = 20

contains


  !> Initializes the file registry.
  subroutine TFileRegistry_init(this)

    !> Initialized instance on exit.
    type(TFileRegistry), intent(out) :: this

    allocate(this%registeredFiles(defaultSize))
    this%nRegistered = 0

  end subroutine TFileRegistry_init


  !> Registers a file into the registry.
  subroutine TFileRegistry_register(this, unit)

    !> Instance
    class(TFileRegistry), intent(inout) :: this

    !> File unit to register
    integer, intent(in) :: unit

    integer, allocatable :: buffer(:)

    if (this%nRegistered == size(this%registeredFiles)) then
      allocate(buffer(2 * size(this%registeredFiles)))
      buffer(1:this%nRegistered) = this%registeredFiles
      call move_alloc(buffer, this%registeredFiles)
    end if
    this%nRegistered = this%nRegistered + 1
    this%registeredFiles(this%nRegistered) = unit

  end subroutine TFileRegistry_register


  !> Checks the status of all registered files and closes those which are still opened.
  subroutine TFileRegistry_closeAll(this)

    !> Instance
    class(TFileRegistry), intent(in) :: this

    integer :: unit, ii
    logical :: tOpen

    do ii = 1, this%nRegistered
      unit = this%registeredFiles(ii)
      inquire(unit=unit, opened=tOpen)
      if (tOpen) then
        close(unit)
      end if
    end do

  end subroutine TFileRegistry_closeAll
  

end module dftbp_fileregistry
