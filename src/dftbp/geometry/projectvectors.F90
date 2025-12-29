!--------------------------------------------------------------------------------------------------!
!  DFTB+: general package for performing fast atomistic simulations                                !
!  Copyright (C) 2006 - 2025  DFTB+ developers group                                               !
!                                                                                                  !
!  See the LICENSE file for terms of usage and distribution.                                       !
!--------------------------------------------------------------------------------------------------!

#:include 'common.fypp'


!> Contains routines for projecting out components of vectors attached to multiple atoms
module dftbp_geometry_projectvectors
  use dftbp_common_accuracy, only : dp
  use dftbp_math_eigensolver, only : heev
  use dftbp_math_matrixops, only : adjointLowerTriangle
  implicit none

  private
  public :: getCentreOfMass, getPrincipleAxes, inertiaTensor

contains

  !> Evaluate centre of mass
  pure function getCentreOfMass(coords, masses, atoms) result(centreOfMass)

    !> Coordinates
    real(dp), intent(in) :: coords(:,:)

    !> Masses
    real(dp), intent(in) :: masses(:)

    !> Atoms over which the moment of inertia tensor is formed
    integer, intent(in), optional :: atoms(:)

    !> Resulting centre of mass
    real(dp) :: centreOfMass(3)

    integer :: ii, iAt, nAtom

    nAtom = size(masses)
    centreOfMass(:) = 0.0_dp
    if (present(atoms)) then
      do ii = 1, nAtom
        iAt = atoms(ii)
        centreOfMass(:) = centreOfMass + coords(:, iAt) * masses(iAt)
      end do
      centreOfMass(:) = centreOfMass / sum(masses(atoms))
    else
      do iAt = 1, nAtom
        centreOfMass(:) = centreOfMass + coords(:, iAt) * masses(iAt)
      end do
      centreOfMass(:) = centreOfMass / sum(masses)
    end if

  end function getCentreOfMass


  !> Principle axes of the moment of inertia tensor
  subroutine getPrincipleAxes(inertia, ei, coords, masses, centreOfMass, atoms)

    !> Inertia axes
    real(dp), intent(out) :: inertia(3,3)

    !> Moments
    real(dp), intent(out) :: ei(3)

    !> Coordinates
    real(dp), intent(in) :: coords(:,:)

    !> Masses
    real(dp), intent(in) :: masses(:)

    !> The centre of mass
    real(dp), intent(in) :: centreOfMass(3)

    !> Atoms over which the moment of inertia tensor is formed
    integer, intent(in), optional :: atoms(:)

    integer :: ii, jj, kk, iAt, nAtom

    if (present(atoms)) then
      nAtom = size(atoms)
    else
      nAtom = size(masses)
    end if

    ei(:) = 0.0_dp
    if (nAtom < 1) return

    call inertiaTensor(inertia, coords, masses, centreOfMass, atoms)

    call heev(inertia, ei, 'L', 'V')

  end subroutine getPrincipleAxes


  !> Calculate the inertia tensor with respect to an origin point
  pure subroutine inertiaTensor(inertia, coords, masses, origin, atoms)

    !> Inertia tensor
    real(dp), intent(out) :: inertia(3,3)

    !> Coordinates
    real(dp), intent(in) :: coords(:,:)

    !> Masses
    real(dp), intent(in) :: masses(:)

    !> The origin of the coordinate system
    real(dp), intent(in) :: origin(3)

    !> Atoms over which the moment of inertia tensor is formed
    integer, intent(in), optional :: atoms(:)

    integer :: iAt, ii, jj, kk

    inertia(:,:) = 0.0_dp
    if (present(atoms)) then

      do kk = 1, size(atoms)
        iAt = atoms(kk)
        do ii = 1, 3
          inertia(ii, ii) = inertia(ii, ii) + masses(iAt)*sum((coords(:,iAt) - origin(:))**2)
          do jj = ii, 3
            inertia(jj, ii) = inertia(jj, ii) - masses(iAt) * (coords(jj,iAt) - origin(jj))&
                & * (coords(ii,iAt) - origin(ii))
          end do
        end do
      end do

    else

      do iAt = 1, size(masses)
        do ii = 1, 3
          inertia(ii, ii) = inertia(ii, ii) + masses(iAt) * sum((coords(:,iAt) - origin(:))**2)
          do jj = ii, 3
            inertia(jj, ii) = inertia(jj, ii) - masses(iAt) * (coords(jj,iAt) - origin(jj))&
                & * (coords(ii,iAt) - origin(ii))
          end do
        end do
      end do

    end if

  end subroutine inertiaTensor

end module dftbp_geometry_projectvectors
