!--------------------------------------------------------------------------------------------------!
!  DFTB+: general package for performing fast atomistic simulations                                !
!  Copyright (C) 2006 - 2020  DFTB+ developers group                                               !
!                                                                                                  !
!  See the LICENSE file for terms of usage and distribution.                                       !
!--------------------------------------------------------------------------------------------------!

!> Provides functions to check relations between versions.
module dftbp_version
  implicit none

  private
  public :: TVersion


  integer, parameter :: maxVersionNumbers = 3


  !> Contains version numbers (major, minor, patch) and provides relational operators
  type :: TVersion

    !> Various version number components
    integer :: numbers(maxVersionNumbers)

  contains

    generic :: operator(>) => greater
    procedure, private :: greater

    generic :: operator(<) => less
    procedure, private :: less

    generic :: operator(==) => equal
    procedure, private :: equal

    generic :: operator(>=) => greater_equal
    procedure, private :: greater_equal

    generic :: operator(<=) => less_equal
    procedure, private :: less_equal

    generic :: operator(/=) => unequal
    procedure, private :: unequal

  end type TVersion


  ! Constructor for TVersion
  interface TVersion
    module procedure TVersion_construct
  end interface TVersion


contains

  !> Constructs a version instance.
  function TVersion_construct(major, minor, patch) result(this)

    !> Major version number
    integer, intent(in) :: major

    !> Minor version number
    integer, intent(in), optional :: minor

    !> Patch version number
    integer, intent(in), optional :: patch

    !> Created instance.
    type(TVersion) :: this

    this%numbers(1) = major
    this%numbers(2:) = 0
    if (present(minor)) then
      this%numbers(2) = minor
    end if
    if (present(patch)) then
      this%numbers(3) = patch
    end if

  end function TVersion_construct


  !> Checks whether two versions are greater.
  elemental function greater(version, refVersion) result(compResult)

    !> Version (LHS of the comparison)
    class(TVersion), intent(in) :: version

    !> Reference version (RHS of the comparison)
    class(TVersion), intent(in) :: refVersion

    !> Result of the comparison
    logical :: compResult

    integer :: ii

    do ii = 1, maxVersionNumbers
      if (version%numbers(ii) /= refVersion%numbers(ii)) then
        compResult = version%numbers(ii) > refVersion%numbers(ii)
        return
      end if
    end do
    compResult = .false.

  end function greater


  !> Checks whether two versions are less
  elemental function less(version, refVersion) result(compResult)

    !> Version (LHS of the comparison)
    class(TVersion), intent(in) :: version

    !> Reference version (RHS of the comparison)
    class(TVersion), intent(in) :: refVersion

    !> Result of the comparison
    logical :: compResult

    compResult = refVersion > version

  end function less


  !> Checks whether two versions are equal.
  elemental function equal(version, refVersion) result(compResult)

    !> Version (LHS of the comparison)
    class(TVersion), intent(in) :: version

    !> Reference version (RHS of the comparison)
    class(TVersion), intent(in) :: refVersion

    !> Result of the comparison
    logical :: compResult

    compResult = .not. (version > refVersion)
    if (compResult) then
      compResult = .not. (refVersion > version)
    end if

  end function equal


  !> Checks whether two versions are greater equal.
  elemental function greater_equal(version, refVersion) result(compResult)

    !> Version (LHS of the comparison)
    class(TVersion), intent(in) :: version

    !> Reference version (RHS of the comparison)
    class(TVersion), intent(in) :: refVersion

    !> Result of the comparison
    logical :: compResult

    compResult = .not. (refVersion > version)

  end function greater_equal


  !> Checks whether two version are less equal.
  elemental function less_equal(version, refVersion) result(compResult)

    !> Version (LHS of the comparison)
    class(TVersion), intent(in) :: version

    !> Reference version (RHS of the comparison)
    class(TVersion), intent(in) :: refVersion

    !> Result of the comparison
    logical :: compResult

    compResult = .not. (version > refVersion)

  end function less_equal


  !> Checks whether two versions are unequal.
  elemental function unequal(version, refVersion) result(compResult)

    !> Version (LHS of the comparison)
    class(TVersion), intent(in) :: version

    !> Reference version (RHS of the comparison)
    class(TVersion), intent(in) :: refVersion

    !> Result of the comparison
    logical :: compResult

    compResult = version > refVersion
    if (.not. compResult) then
      compResult = refVersion > version
    end if

  end function unequal


end module dftbp_version
