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
  public :: getCentreOfMass, getPrincipleAxes


contains

  !> Evaluate centre of mass
  function getCentreOfMass(coords, masses, atoms) result(centreOfMass)

    !> Coordinates
    real(dp), intent(in) :: coords(:,:)

    !> Masses
    real(dp), intent(in) :: masses(:)

    !> Atoms over which the moment of intertia tensor is formed
    integer, intent(in) :: atoms(:)

    !> Resulting centre of mass
    real(dp) :: centreOfMass(3)

    integer :: ii, iAt, nAtom

    nAtom = size(atoms)
    centreOfMass(:) = 0.0_dp

    !$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(iAt) SCHEDULE(RUNTIME) REDUCTION(+:centreOfMass)
    do ii = 1, nAtom
      iAt = atoms(ii)
      centreOfMass(:) = centreOfMass + coords(:, iAt) * masses(iAt)
    end do
    !$OMP END PARALLEL DO
    centreOfMass(:) = centreOfMass / sum(masses(atoms))

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

    !> Atoms over which the moment of intertia tensor is formed
    integer, intent(in) :: atoms(:)

    integer :: ii, jj, kk, iAt, nAtom

    nAtom = size(atoms)
    inertia(:,:) = 0.0_dp
    ei(:) = 0.0_dp
    if (nAtom < 1) return

    !$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(iAt, ii, jj) SCHEDULE(RUNTIME) REDUCTION(+:inertia)
    do kk = 1, nAtom
      iAt = atoms(kk)
      do ii = 1, 3
        inertia(ii, ii) = inertia(ii, ii) + masses(iAt) * sum((coords(:,iAt) - centreOfMass(:))**2)
        do jj = ii, 3
          inertia(jj, ii) = inertia(jj, ii) - masses(iAt) * (coords(jj,iAt) - centreOfMass(jj))&
              & * (coords(ii,iAt) - centreOfMass(ii))
        end do
      end do
    end do
    !$OMP END PARALLEL DO

    call heev(inertia, ei, 'L', 'V')

  end subroutine getPrincipleAxes

end module dftbp_geometry_projectvectors
