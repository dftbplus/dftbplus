!--------------------------------------------------------------------------------------------------!
!  DFTB+: general package for performing fast atomistic simulations                                !
!  Copyright (C) 2006 - 2025 DFTB+ developers group                                                !
!                                                                                                  !
!  See the LICENSE file for terms of usage and distribution.                                       !
!--------------------------------------------------------------------------------------------------!

!> Contains routines to symmetrize parts of a calculation
module dftbp_math_symmetrize
  use dftbp_common_accuracy, only : dp
  implicit none

  !> Atoms that are equivalent in some sense
  type TEquivalentAtoms

    !> Number of atoms that are equivalents
    integer, allocatable :: nMappedAtoms(:)

    !> List of equivalent atoms
    integer, allocatable :: iMappedAtoms(:,:)

  contains

    procedure, private :: map0D
    procedure, private :: map1D
    procedure, private :: map2D
    generic :: mapAtomicQuantities => map0D, map1D, map2D

  end type TEquivalentAtoms


  private
  public :: initTEquivalentAtoms, TEquivalentAtoms

contains


  !> Processing for transport contact calculations only at the moment, mapping atoms in two
  !! principle layers to be identical
  subroutine initTEquivalentAtoms(this, nAtom)

    !> Instance
    class(TEquivalentAtoms), intent(out) :: this

    !> Number of atoms in system
    integer, intent(in) :: nAtom

    integer :: iAt

    allocate(this%nMappedAtoms(nAtom/2), source=2)
    allocate(this%iMappedAtoms(2, nAtom/2), source=0)
    do iAt = 1, nAtom/2
      this%iMappedAtoms(1, iAt) = iAt
    end do
    this%iMappedAtoms(2, :nAtom/2) = this%iMappedAtoms(1, :nAtom/2) + nAtom/2

  end subroutine initTEquivalentAtoms


  !> Map (spin channel resolved) atom charges
  subroutine map0D(this, q)

    !> Instance
    class(TEquivalentAtoms), intent(in) :: this

    !> Charges
    real(dp), intent(inout) :: q(:, :)

    integer :: iGrp, nAtoms, iAtoms
    real(dp) :: work(size(q, dim=2))

    do iGrp = 1, size(this%nMappedAtoms)
      nAtoms = this%nMappedAtoms(iGrp)
      work(:) = sum(q(this%iMappedAtoms(:nAtoms, iGrp), :), dim=1)
      work(:) = work / real(nAtoms, dp)
      do iAtoms = 1, nAtoms
        q(this%iMappedAtoms(iAtoms, iGrp), :) = work
      end do
    end do

  end subroutine map0D


  !> Map (spin channel resolved) atomic shell or orbital charges
  subroutine map1D(this, q)

    !> Instance
    class(TEquivalentAtoms), intent(in) :: this

    !> Charges
    real(dp), intent(inout) :: q(:, :, :)

    integer :: iGrp, nAtoms, iAtoms
    real(dp) :: work(size(q, dim=1), size(q, dim=3))

    do iGrp = 1, size(this%nMappedAtoms)
      nAtoms = this%nMappedAtoms(iGrp)
      work(:,:) = sum(q(:, this%iMappedAtoms(:nAtoms, iGrp), :), dim=2)
      work(:,:) = work / real(nAtoms, dp)
      do iAtoms = 1, nAtoms
        q(:, this%iMappedAtoms(iAtoms, iGrp), :) = work
      end do
    end do

  end subroutine map1D


  !> Map (spin channel resolved) atomic block charges
  subroutine map2D(this, q)

    !> Instance
    class(TEquivalentAtoms), intent(in) :: this

    !> Charges
    real(dp), intent(inout) :: q(:, :, :, :)

    integer :: iGrp, nAtoms, iAtoms
    real(dp) :: work(size(q, dim=1), size(q, dim=2), size(q, dim=4))

    do iGrp = 1, size(this%nMappedAtoms)
      nAtoms = this%nMappedAtoms(iGrp)
      work(:,:,:) = sum(q(:, :, this%iMappedAtoms(:nAtoms, iGrp), :), dim=3)
      work(:,:,:) = work / real(nAtoms, dp)
      do iAtoms = 1, nAtoms
        q(:, :, this%iMappedAtoms(iAtoms, iGrp), :) = work
      end do
    end do

  end subroutine map2D

end module dftbp_math_symmetrize
