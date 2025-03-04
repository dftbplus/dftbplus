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
  use dftbp_math_blasroutines, only : herk
  use dftbp_math_eigensolver, only : heev
  use dftbp_math_simplealgebra, only : cross3
  use dftbp_type_typegeometry, only : TGeometry
  use dftbp_type_densedescr, only : TDenseDescr
  use dftbp_math_matrixops, only : adjointLowerTriangle
#:if WITH_MPI
  use dftbp_common_environment, only : TEnvironment
  use dftbp_math_matrixops, only : adjointLowerTriangle_BLACS
#:endif
  implicit none

  private
  public :: project
#:if WITH_MPI
  public :: projectBlacs
#:endif


contains

#:if WITH_MPI
  !> Projection out of the space of the dynamical matrix.
  subroutine projectBlacs(env, denseDesc, dynMatrix, tRemoveTranslate, tRemoveRotate, nDerivs,&
      & nMovedAtom, geo, atomicMasses)

    !> Environment settings
    type(TEnvironment), intent(in) :: env

    !> Dense matrix descriptor
    type(TDenseDescr), intent(in) :: denseDesc

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
    real(dp) :: centerOfMass(3), rTmp(3), vTmp(3), inertia(3,3), moments(3)
    integer :: nToNull, ii, jj, iAt

    nToNull = 0
    if (tRemoveTranslate) nToNull = nToNull + 3
    if (tRemoveRotate) nToNull = nToNull + 3
    if (nToNull == 0) return

    allocate(vectorsToNull(nDerivs, nToNull), source=0.0_dp)
    allocate(projector(nDerivs, nDerivs), source=0.0_dp)

    ! symmetrize dynamical matrix
    call adjointLowerTriangle_BLACS(denseDesc%blacsOrbSqr, env%blacs%orbitalGrid%myCol,&
        & env%blacs%orbitalGrid%myRow, env%blacs%orbitalGrid%nCol, env%blacs%orbitalGrid%nRow,&
        & dynMatrix)

    do jj = 1, nDerivs
      projector(jj, jj) = 1.0_dp
    end do

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

      call getCenterOfMass(nMovedAtom, atomicMasses, geo%coords, centerOfMass)
      call getPrincipleAxes(geo%coords, atomicMasses, centerOfMass, nMovedAtom, inertia, moments)

      ! axis to project with respect to
      do ii = 1, 3
        if (moments(ii) < epsilon(0.0_dp)) then
          ! zero moment of inertia - linear molecule, and this direction is along its axis
          cycle
        end if
        vTmp(:) = inertia(:,ii)
        do iAt = 1, nMovedAtom
          rTmp = cross3(vTmp, geo%coords(:,iAt) - centerOfMass)
          vectorsToNull((iAt - 1) * 3 + 1 : iAt * 3, nToNull - ii + 1) = rTmp
        end do
      end do
    end if

    ! change from displacements to weighted displacements basis of the Hessian
    do iAt = 1, nMovedAtom
      vectorsToNull((iAt - 1) * 3 + 1 : iAt * 3, :) = vectorsToNull((iAt - 1) * 3 + 1 : iAt * 3, :)&
          & * sqrt(atomicMasses(iAt))
    end do

    ! normalise non-null vectors
    do ii = 1, nToNull
      if (sum(vectorsToNull(:, ii)**2) > epsilon(1.0_dp)) then
        vectorsToNull(:, ii) = vectorsToNull(:, ii) / sqrt(sum(vectorsToNull(:, ii)**2))
      end if
    end do

    ! projector = projector - vectorsToNull * vectorsToNull^T
    call herk(projector, vectorsToNull, alpha=-1.0_dp, beta=1.0_dp)

    ! copy to other triangle
    call adjointLowerTriangle_BLACS(denseDesc%blacsOrbSqr, env%blacs%orbitalGrid%myCol,&
        & env%blacs%orbitalGrid%myRow, env%blacs%orbitalGrid%nCol, env%blacs%orbitalGrid%nRow,&
        & projector)

    ! project out removed degrees of freedom
    dynMatrix(:,:) = matmul(projector, matmul(dynMatrix, projector))

  end subroutine projectBlacs

#:endif

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
    real(dp) :: centerOfMass(3), rTmp(3), vTmp(3), inertia(3,3), moments(3)
    integer :: nToNull, ii, jj, iAt

    nToNull = 0
    if (tRemoveTranslate) nToNull = nToNull + 3
    if (tRemoveRotate) nToNull = nToNull + 3
    if (nToNull == 0) return

    allocate(vectorsToNull(nDerivs, nToNull), source=0.0_dp)
    allocate(projector(nDerivs, nDerivs), source=0.0_dp)

    ! symmetrize dynamical matrix
    call adjointLowerTriangle(dynMatrix)

    do jj = 1, nDerivs
      projector(jj, jj) = 1.0_dp
    end do

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

      call getCenterOfMass(nMovedAtom, atomicMasses, geo%coords, centerOfMass)
      call getPrincipleAxes(geo%coords, atomicMasses, centerOfMass, nMovedAtom, inertia, moments)

      ! axis to project with respect to
      do ii = 1, 3
        if (moments(ii) < epsilon(0.0_dp)) then
          ! zero moment of inertia - linear molecule, and this direction is along its axis
          cycle
        end if
        vTmp(:) = inertia(:,ii)
        do iAt = 1, nMovedAtom
          rTmp = cross3(vTmp, geo%coords(:,iAt) - centerOfMass)
          vectorsToNull((iAt - 1) * 3 + 1 : iAt * 3, nToNull - ii + 1) = rTmp
        end do
      end do
    end if

    ! change from displacements to weighted displacements basis of the Hessian
    do iAt = 1, nMovedAtom
      vectorsToNull((iAt - 1) * 3 + 1 : iAt * 3, :) = vectorsToNull((iAt - 1) * 3 + 1 : iAt * 3, :)&
          & * sqrt(atomicMasses(iAt))
    end do

    ! normalise non-null vectors
    do ii = 1, nToNull
      if (sum(vectorsToNull(:, ii)**2) > epsilon(1.0_dp)) then
        vectorsToNull(:, ii) = vectorsToNull(:, ii) / sqrt(sum(vectorsToNull(:, ii)**2))
      end if
    end do

    call herk(projector, vectorsToNull, alpha=-1.0_dp, beta=1.0_dp)

    ! copy to other triangle
    call adjointLowerTriangle(projector)

    ! project out removed degrees of freedom
    dynMatrix(:,:) = matmul(projector, matmul(dynMatrix, projector))

  end subroutine project


  !> Calculates the center of mass of a given geometry.
  subroutine getCenterOfMass(nMovedAtom, atomicMasses, coords, centerOfMass)

    !> Number of moving atoms
    integer, intent(in) :: nMovedAtom

    !> Masses
    real(dp), intent(in) :: atomicMasses(:)

    !> Coordinates
    real(dp), intent(in) :: coords(:,:)

    !> The center of mass
    real(dp), intent(out) :: centerOfMass(3)

    !! Auxiliary variables
    integer :: iAt

    centerOfMass(:) = 0.0_dp
    do iAt = 1, nMovedAtom
      centerOfMass(:) = centerOfMass + coords(:, iAt) * atomicMasses(iAt)
    end do
    centerOfMass(:) = centerOfMass / sum(atomicMasses(:nMovedAtom))

  end subroutine getCenterOfMass


  !> Determines principle moment of inertia axes.
  subroutine getPrincipleAxes(coords, atomicMasses, centerOfMass, nMovedAtom, inertia, ei)

    !> Coordinates
    real(dp), intent(in) :: coords(:,:)

    !> Masses
    real(dp), intent(in) :: atomicMasses(:)

    !> The center of mass
    real(dp), intent(in) :: centerOfMass(3)

    !> Number of moving atoms
    integer, intent(in) :: nMovedAtom

    !> Intertia axes
    real(dp), intent(out) :: inertia(3, 3)

    !> Moments
    real(dp), intent(out) :: ei(3)

    integer :: ii, jj, iAt

    inertia(:,:) = 0.0_dp
    ei(:) = 0.0_dp

    do iAt = 1, nMovedAtom
      do ii = 1, 3
        inertia(ii, ii) = inertia(ii, ii)&
            & + atomicMasses(iAt) * sum((coords(:,iAt) - centerOfMass)**2)
      end do
    end do
    do iAt = 1, nMovedAtom
      do ii = 1, 3
        do jj = 1, 3
          inertia(jj, ii) = inertia(jj, ii)&
              & - atomicMasses(iAt) * (coords(jj,iAt) - centerOfMass(jj))&
              & * (coords(ii,iAt) - centerOfMass(ii))
        end do
      end do
    end do

    call heev(inertia, ei, 'U', 'V')

  end subroutine getPrincipleAxes

end module modes_modeprojection
