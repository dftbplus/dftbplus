!--------------------------------------------------------------------------------------------------!
!  DFTB+: general package for performing fast atomistic simulations                                !
!  Copyright (C) 2006 - 2020  DFTB+ developers group                                               !
!                                                                                                  !
!  See the LICENSE file for terms of usage and distribution.                                       !
!--------------------------------------------------------------------------------------------------!

!> Provides functions to check relations between versions.
module dftbp_versioncheck
  implicit none

  private
  public :: version_equal, version_greater, version_greater_equal, version_less, version_less_equal

contains

  #:for name, op in [('greater', '>'), ('less', '<')]

    !> Checks whether two versions are greater/less.
    pure function version_${name}$(version, refVersion) result(compResult)

      !> Version (LHS of the comparison)
      integer, intent(in) :: version(:)

      !> Reference version (RHS of the comparison)
      integer, intent(in) :: refVersion(:)

      !> Result of the comparison
      logical :: compResult

      integer :: versionFull(max(size(version), size(refVersion)))
      integer :: refVersionFull(size(versionFull))
      integer :: ii

      versionFull(:) = 0
      versionFull(1 : size(version)) = version
      refVersionFull(:) = 0
      refVersionFull(1 : size(refVersion)) = refVersion

      do ii = 1, size(versionFull)
        if (versionFull(ii) /= refVersionFull(ii)) then
          compResult = versionFull(ii) ${op}$ refVersionFull(ii)
          return
        end if
      end do
      compResult = .false.

    end function version_${name}$

  #:endfor


  !> Checks whether two versions are equal.
  pure function version_equal(version, refVersion) result(compResult)

    !> Version (LHS of the comparison)
    integer, intent(in) :: version(:)

    !> Reference version (RHS of the comparison)
    integer, intent(in) :: refVersion(:)

    !> Result of the comparison
    logical :: compResult

    integer :: versionFull(max(size(version), size(refVersion)))
    integer :: refVersionFull(size(versionFull))

    versionFull(:) = 0
    versionFull(1 : size(version)) = version
    refVersionFull(:) = 0
    refVersionFull(1 : size(refVersion)) = refVersion
    compResult = all(versionFull == refVersionFull)

  end function version_equal


  !> Checks whether two versions are greater equal.
  pure function version_greater_equal(version, refVersion) result(compResult)

    !> Version (LHS of the comparison)
    integer, intent(in) :: version(:)

    !> Reference version (RHS of the comparison)
    integer, intent(in) :: refVersion(:)

    !> Result of the comparison
    logical :: compResult

    compResult = version_equal(version, refVersion) .or. version_greater(version, refVersion)

  end function version_greater_equal


  !> Checks whether two version are less equal.
  pure function version_less_equal(version, refVersion) result(compResult)

    !> Version (LHS of the comparison)
    integer, intent(in) :: version(:)

    !> Reference version (RHS of the comparison)
    integer, intent(in) :: refVersion(:)

    !> Result of the comparison
    logical :: compResult

    compResult = version_equal(version, refVersion) .or. version_less(version, refVersion)

  end function version_less_equal


end module dftbp_versioncheck
