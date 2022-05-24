! NAG 7.1 can not use user-defined derived type I/O, when derived type is part of an array.
! Workaround: Make a local copy and use that in the I/O statement.
module unitconversion
  implicit none

  private
  public :: TUnit
  public :: lengthUnits

  type TUnit
    private
    character(20) :: name
    real :: conversionFact = 1.0
  contains
    procedure :: writeFormatted => TUnit_writeFormatted
    generic :: write(formatted) => writeFormatted
  end type TUnit

  type(TUnit), parameter :: lengthUnits(*) = [&
      & TUnit("angstrom            ", 1.0),&
      & TUnit("aa                  ", 1.0)&
      & ]

contains

  subroutine TUnit_writeFormatted(this, unit, iotype, vlist, iostat, iomsg)
    class(TUnit), intent(in) :: this
    integer, intent(in) :: unit
    character(*), intent(in) :: iotype
    integer, intent(in) :: vlist(:)
    integer, intent(out) :: iostat
    character(*), intent(inout) :: iomsg

    write(unit,"(a20, a, e24.14, a)") this%name, ':', this%conversionFact, " * x"
    iostat = 0

  end subroutine TUnit_writeFormatted

end module unitconversion


program printunits
  use unitconversion, only : lengthUnits, TUnit
  implicit none

  integer :: ii
  type(TUnit) :: localUnit
  do ii = 1, size(lengthUnits)
    ! Not working:
    ! write(*,"(1x,dt)") lengthUnits(ii)
    ! Workaround
    localUnit = lengthUnits(ii)
    write(*,"(1x,dt)") localUnit
  end do

end program printunits
