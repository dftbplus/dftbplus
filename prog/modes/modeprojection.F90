!--------------------------------------------------------------------------------------------------!
!  DFTB+: general package for performing fast atomistic simulations                                !
!  Copyright (C) 2006 - 2020  DFTB+ developers group                                               !
!                                                                                                  !
!  See the LICENSE file for terms of usage and distribution.                                       !
!--------------------------------------------------------------------------------------------------!

#:include 'common.fypp'

!> Removal of translation or rotation related modes
module dftbp_modeprojection
  use dftbp_accuracy, only : dp
  use dftbp_message
  use dftbp_typegeometry
  use dftbp_blasroutines, only : herk
  use dftbp_simplealgebra, only : cross3
  use dftbp_eigensolver, only : heev
  implicit none

  private
  public :: project

contains

  !> Projection out of the space of the dynamical matrix
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

    nToNull = 0
    if (tRemoveTranslate) then
      nToNull = nToNull + 3
    end if
    if (tRemoveRotate) then
      nToNull = nToNull + 3
    end if

    if (nToNull == 0) then
      return
    end if

    allocate(vectorsToNull(nDerivs, nToNull))
    allocate(projector(nDerivs, nDerivs))
    projector(:,:) = 0.0_dp
    vectorsToNull(:,:) = 0.0_dp

    ! symmetrize dynamical matrix
    do jj = 1, nDerivs
      dynMatrix(jj, jj + 1 :) = dynMatrix(jj + 1 :, jj)
    end do

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
      centreOfMass(:) = 0.0_dp
      do iAt = 1, nMovedAtom
        centreOfMass(:) = centreOfMass + geo%coords(:,iAt) * atomicMasses(iAt)
      end do
      centreOfMass(:) = centreOfMass / sum(atomicMasses(:nMovedAtom))

      call getPrincipleAxes(inertia, moments, geo%coords, atomicMasses, centreOfMass, nMovedAtom)

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
      vectorsToNull((iAt - 1) * 3 + 1 : iAt * 3, :) =&
          & vectorsToNull((iAt - 1) * 3 + 1 : iAt * 3, :) * sqrt(atomicMasses(iAt))
    end do

    ! normalise non-null vectors
    do ii = 1, nToNull
      if (sum(vectorsToNull(:,ii)**2) > epsilon(1.0_dp)) then
        vectorsToNull(:,ii) = vectorsToNull(:,ii) / sqrt(sum(vectorsToNull(:,ii)**2))
      end if
    end do

    call herk(projector, vectorsToNull, alpha=-1.0_dp, beta=1.0_dp)

    deallocate(vectorsToNull)

    ! copy to other triangle
    do jj = 1, nDerivs
      projector(jj,jj+1:) = projector(jj+1:,jj)
    end do

    ! project out removed degrees of freedom
    dynMatrix(:,:) = matmul(projector, matmul(dynMatrix, projector))

    deallocate(projector)

  end subroutine project


  !> Principle moment of inertia axes
  subroutine getPrincipleAxes(inertia, ei, coords, masses, centreOfMass, nMovedAtom)

    !> Intertia axes
    real(dp), intent(out) :: inertia(3,3)

    !> Moments
    real(dp), intent(out) :: ei(3)

    !> coordinates
    real(dp), intent(in) :: coords(:,:)

    !> masses
    real(dp), intent(in) :: masses(:)

    !> the centre of mass
    real(dp), intent(in) :: centreOfMass(3)

    !> number of moving atoms
    integer, intent(in) :: nMovedAtom

    integer :: ii, jj, iAt

    inertia(:,:) = 0.0_dp
    ei(:) = 0.0_dp

    do iAt = 1, nMovedAtom
      do ii = 1, 3
        inertia(ii, ii) = inertia(ii, ii) + masses(iAt) * sum((coords(:,iAt) - centreOfMass(:))**2)
      end do
    end do
    do iAt = 1, nMovedAtom
      do ii = 1, 3
        do jj = 1, 3
          inertia(jj, ii) = inertia(jj, ii) - masses(iAt) * (coords(jj,iAt) - centreOfMass(jj))&
              & * (coords(ii,iAt) - centreOfMass(ii))
        end do
      end do
    end do

    call heev(inertia, ei, 'U', 'V')

  end subroutine getPrincipleAxes

end module dftbp_modeprojection
