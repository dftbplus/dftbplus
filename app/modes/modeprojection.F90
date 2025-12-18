!--------------------------------------------------------------------------------------------------!
!  DFTB+: general package for performing fast atomistic simulations                                !
!  Copyright (C) 2006 - 2025  DFTB+ developers group                                               !
!                                                                                                  !
!  See the LICENSE file for terms of usage and distribution.                                       !
!--------------------------------------------------------------------------------------------------!

#:include 'common.fypp'

!> Removal of translation or rotation related modes.
module modes_modeprojection
  use dftbp_common_accuracy, only : dp
  use dftbp_io_message, only : warning
  use dftbp_geometry_projectvectors, only : getCentreOfMass, getPrincipleAxes
  use dftbp_math_blasroutines, only : herk
  use dftbp_math_matrixops, only : adjointLowerTriangle, makeSimilarityTrans
  use dftbp_math_simplealgebra, only : cross3
  use dftbp_type_typegeometry, only : TGeometry
  implicit none

  private
  public :: project


contains

  !> Projection out of the space of the dynamical matrix.
  subroutine project(dynMatrix, tRemoveTranslate, tRemoveRotate, nDerivs, nMovedAtom, geo,&
      & atomicMasses)

    !> Dynamical matrix
    real(dp), intent(inout) :: dynMatrix(:,:)

    !> Whether translation is removed
    logical, intent(in) :: tRemoveTranslate

    !> Whether rotation is removed
    logical, intent(in) :: tRemoveRotate

    !> Number of derivatives (dimension of dynMatrix)
    integer, intent(in) :: nDerivs

    !> Number of moving atoms
    integer, intent(in) :: nMovedAtom

    !> Geometry information on the system
    type(TGeometry), intent(in) :: geo

    !> Masses of the atoms in the system
    real(dp), intent(in) :: atomicMasses(:)

    real(dp), allocatable :: vectorsToNull(:,:), projector(:,:)
    real(dp) :: centreOfMass(3), rTmp(3), vTmp(3), inertia(3,3), moments(3)
    integer :: nToNull, ii, jj, iAt
    integer, allocatable :: movedAtoms(:)

    nToNull = 0
    if (tRemoveTranslate) nToNull = nToNull + 3
    if (tRemoveRotate) then
      nToNull = nToNull + 3
      allocate(movedAtoms(nMovedAtom))
      do ii = 1, nMovedAtom
        movedAtoms(ii) = ii
      end do
    end if
    if (nToNull == 0) return

    allocate(vectorsToNull(nDerivs, nToNull), source=0.0_dp)

    ! symmetrize dynamical matrix
    call adjointLowerTriangle(dynMatrix)

    if (tRemoveTranslate) then
      do iAt = 1, nMovedAtom
        do ii = 1, 3
          vectorsToNull((iAt - 1) * 3 + ii, ii) = 1.0_dp
        end do
      end do
    end if

    if (tRemoveRotate) then
      if (geo%tPeriodic) then
        call warning("Rotational modes were requested to be removed for a periodic geometry -&
            & results probably unphysical!")
      end if

      centreOfMass(:) = getCentreOfMass(geo%coords, atomicMasses, movedAtoms)

      call getPrincipleAxes(inertia, moments, geo%coords, atomicMasses, centreOfMass, movedAtoms)

      ! axis to project with respect to
      do ii = 1, 3
        if (moments(ii) < epsilon(0.0_dp)) then
          ! zero moment of inertia - linear molecule, and this direction is along its axis
          cycle
        end if
        vTmp(:) = inertia(:,ii)
        do iAt = 1, nMovedAtom
          rTmp = cross3(vTmp, geo%coords(:,iAt) - centreOfMass)
          vectorsToNull((iAt - 1) * 3 + 1 : iAt * 3, nToNull - ii + 1) = rTmp
        end do
      end do
    end if

    ! Change from displacements to weighted displacements basis of the Hessian
    do iAt = 1, nMovedAtom
      vectorsToNull((iAt - 1) * 3 + 1 : iAt * 3, :) = vectorsToNull((iAt - 1) * 3 + 1 : iAt * 3, :)&
          & * sqrt(atomicMasses(iAt))
    end do

    ! normalise vectors
    do ii = 1, nToNull
      if (sum(vectorsToNull(:,ii)**2) > epsilon(1.0_dp)) then
        vectorsToNull(:,ii) = vectorsToNull(:,ii) / sqrt(sum(vectorsToNull(:,ii)**2))
      else
        vectorsToNull(:,ii) = 0.0_dp
      end if
    end do

    ! Form projector matrix
    allocate(projector(nDerivs, nDerivs), source=0.0_dp)
    !$OMP PARALLEL DO DEFAULT(SHARED) SCHEDULE(RUNTIME)
    do jj = 1, nDerivs
      projector(jj, jj) = 1.0_dp
    end do
    !$OMP END PARALLEL DO
    call herk(projector, vectorsToNull, alpha=-1.0_dp, beta=1.0_dp)

    deallocate(vectorsToNull)

    ! copy to other triangle
    call adjointLowerTriangle(projector)

    ! project out removed degrees of freedom
    call makeSimilarityTrans(dynMatrix, projector)

    deallocate(projector)

  end subroutine project

end module modes_modeprojection
