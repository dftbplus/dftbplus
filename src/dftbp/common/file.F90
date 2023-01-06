!--------------------------------------------------------------------------------------------------!
!  DFTB+: general package for performing fast atomistic simulations                                !
!  Copyright (C) 2006 - 2023  DFTB+ developers group                                               !
!                                                                                                  !
!  See the LICENSE file for terms of usage and distribution.                                       !
!--------------------------------------------------------------------------------------------------!

!> Contains a file descriptor implementation.
module dftbp_common_file
  use dftbp_io_message, only : error
  implicit none

  private
  public :: TFile, TFile_create, TFileOptions


  !> Default file mode settings
  type :: TFileOptions

    !> File access ("sequential", "direct", "stream")
    character(20) :: access = "sequential"

    !> File action ("read", "write", "readwrite")
    character(20) :: action = "read"

    !> File format ("formatted", "unformatted")
    character(20) :: form = "formatted"

    !> File status ("unknown", "old", "replace")
    character(20) :: status = "unknown"

  end type TFileOptions


  !> Implements a file resource
  type :: TFile

    !> File unit
    !>
    !> The file unit should be treated as read-only field and never changed manually.
    !> Also, never call any Fortran statements, which change the unit association status, such
    !> as open or close. If you wish to close the file, deallocate the file instance.
    !>
    integer :: unit = -1

  contains
    final :: TFile_final
    procedure, pass(rhs) :: TFile_assign
    generic :: assignment(=) => TFile_assign
  end type TFile


contains

  !> Opens a file and returns the descriptor.
  !>
  !> If the file could not be opened, the descriptor will be unallocated on return
  !>
  subroutine TFile_create(this, name, options, ioStat)

    !> File descriptor on exit, or unallocated if the file could not be opened
    type(TFile), allocatable, intent(out) :: this

    !> Name of the file to open
    character(*), intent(in) :: name

    !> File options
    type(TFileOptions), optional, intent(in) :: options

    !> I/O stat error generated during open
    integer, optional, intent(out) :: ioStat

    type(TFileOptions) :: opts

    if (present(options)) then
      opts = options
    end if
    allocate(this)
    if (present(ioStat)) then
      open(newunit=this%unit, file=name, access=opts%access, action=opts%action, form=opts%form,&
          & status=opts%status, iostat=ioStat)
      if (ioStat /= 0) deallocate(this)
    else
      open(newunit=this%unit, file=name, access=opts%access, action=opts%action, form=opts%form,&
          & status=opts%status)
    end if

  end subroutine TFile_create


  !> Finalizes (closes) a file
  subroutine TFile_final(this)

    !> Instance
    type(TFile), intent(inout) :: this

    if (this%unit /= -1) then
      close(this%unit)
    end if

  end subroutine TFile_final


  !> Assignment (should never be called as it stops the code)
  subroutine TFile_assign(lhs, rhs)

    !> Right hand side of the assignment
    class(*), intent(inout) :: lhs

    !> Left hand side of the assignment
    class(TFile), intent(in) :: rhs

    call error("Internal error: TFile object is not permitted on the RHS of an assignment")

  end subroutine TFile_assign


end module dftbp_common_file
