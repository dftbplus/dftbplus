!--------------------------------------------------------------------------------------------------!
!  DFTB+: general package for performing fast atomistic simulations                                !
!  Copyright (C) 2006 - 2022  DFTB+ developers group                                               !
!                                                                                                  !
!  See the LICENSE file for terms of usage and distribution.                                       !
!--------------------------------------------------------------------------------------------------!

#:include 'common.fypp'
#:include 'error.fypp'

!> Contains geometrical boundary condition information on the calculation
module dftbp_dftb_boundarycond
  use dftbp_common_accuracy, only : dp
  use dftbp_common_constants, only : pi
  use dftbp_common_status, only : TStatus
  use dftbp_math_angmomentum, only : rotateZ
  use dftbp_math_quaternions, only : rotate3
  use dftbp_math_simplealgebra, only : invert33
  use dftbp_type_commontypes, only : TOrbitals
  implicit none

  private
  public :: zAxis, boundaryConditions, TBoundaryConditions, TBoundaryConditions_init

#:set FLAVOURS = [('complex', 'cplx'), ('real', 'real')]

  !> z direction vector for screw rotation
  real(dp), parameter :: zAxis(3) = [0.0_dp, 0.0_dp, 1.0_dp]


  !> Enumerator containing possible boundary conditions
  type :: TBoundaryConditionEnum_

    !> Unknown/undefined boundary condition
    integer :: unknown = 0

    !> Finite molecular cluster
    integer :: cluster = 1

    !> Three dimensional infinite periodic boundary conditions
    integer :: pbc3d = 2

    !> Helical boundary conditions
    integer :: helical = 3

  end type TBoundaryConditionEnum_


  !> Actual instance of the boundary condition enumerator
  type(TBoundaryConditionEnum_), parameter :: boundaryConditions = TBoundaryConditionEnum_()


  type TBoundaryConditions

    integer :: iBoundaryCondition = boundaryConditions%unknown

  contains

    procedure :: foldInDiatomicBlock_real
    procedure :: foldInDiatomicBlock_cplx
    generic :: foldInDiatomicBlock => foldInDiatomicBlock_real, foldInDiatomicBlock_cplx

    procedure :: foldOutDiatomicBlock_real
    procedure :: foldOutDiatomicBlock_cplx
    generic :: foldOutDiatomicBlock => foldOutDiatomicBlock_real, foldOutDiatomicBlock_cplx

    procedure :: alignVectorCentralCell

    procedure :: foldCoordsToCell

  end type TBoundaryConditions

contains

  subroutine TBoundaryConditions_init(this, iBoundaryCondition, errStatus)

    !> Instance
    type(TBoundaryConditions), intent(out) :: this

    !> Boundary condition choice
    integer, intent(in) :: iBoundaryCondition

    !> Status of routine
    type(TStatus), intent(out) :: errStatus

    if (.not. any([boundaryConditions%cluster, boundaryConditions%pbc3d,&
        & boundaryConditions%helical] == iBoundaryCondition)) then
      @:RAISE_ERROR(errStatus, -1, "Unknown boundary condition specified")
    end if

    this%iBoundaryCondition = iBoundaryCondition

  end subroutine TBoundaryConditions_init


#:for TYPE, LABEL in FLAVOURS

  !> convert matrix elements from image to central cell atoms
  pure subroutine foldInDiatomicBlock_${LABEL}$(this, sqrBlock, iAt1, iAt2, coords, species,&
      & img2CentCell, orb)

    !> Instance
    class(TBoundaryConditions), intent(in) :: this

    !> Basis set expanded quantity to manipulate as per boundary conditions
    ${TYPE}$(dp), intent(inout) :: sqrBlock(:,:)

    !> Atom already in central cell
    integer, intent(in) :: iAt1

    !> Atom to be folded in
    integer, intent(in) :: iAt2

    !> Coordinates of all atoms
    real(dp), intent(in) :: coords(:,:)

    !> Species of each atom
    integer, intent(in) :: species(:)

    !> Map from images of atoms to central cell atoms
    integer, intent(in) :: img2CentCell(:)

    !> data type for atomic orbital information
    type(TOrbitals), intent(in) :: orb

    integer :: iAt2f, iSp, iSh, nOrb1, nOrb2, lShellVals(orb%mShell)
    real(dp) :: theta, rotZ(orb%mOrb,orb%mOrb)

    select case(this%iBoundaryCondition)

    case(boundaryConditions%cluster, boundaryConditions%pbc3d)

      ! orbitals all in the same orientation in space

      return

    case(boundaryConditions%helical)

      if (iAt1 == iAt2) then
        return
      end if
      iAt2f = img2CentCell(iAt2)
      iSp = species(iAt2f)
      iSh = orb%nShell(iSp)
      nOrb1 = orb%nOrbAtom(iAt1)
      nOrb2 = orb%nOrbAtom(iAt2f)
      lShellVals(:iSh) = orb%angShell(:iSh,iSp)
      theta = atan2(coords(2,iAt2f),coords(1,iAt2f)) - atan2(coords(2,iAt2),coords(1,iAt2))
      theta = mod(theta,2.0_dp*pi)
      call rotateZ(rotZ, lShellVals(:iSh), theta)
      sqrBlock(:nOrb2, :nOrb1) = matmul(rotZ(:nOrb2, :nOrb2), sqrBlock(:nOrb2, :nOrb1))

    end select

  end subroutine foldInDiatomicBlock_${LABEL}$


  !> convert matrix elements from central cell to image atoms
  pure subroutine foldOutDiatomicBlock_${LABEL}$(this, sqrBlock, iAt1, iAt2, coords, species,&
      & img2CentCell, orb)

    !> Instance
    class(TBoundaryConditions), intent(in) :: this

    !> Basis set expanded quantity to manipulate as per boundary conditions
    ${TYPE}$(dp), intent(inout) :: sqrBlock(:,:)

    !> Atom already in central cell
    integer, intent(in) :: iAt1

    !> Atom to be folded in
    integer, intent(in) :: iAt2

    !> Coordinates of all atoms
    real(dp), intent(in) :: coords(:,:)

    !> Species of each atom
    integer, intent(in) :: species(:)

    !> Map from images of atoms to central cell atoms
    integer, intent(in) :: img2CentCell(:)

    !> data type for atomic orbital information
    type(TOrbitals), intent(in) :: orb

    integer :: iAt2f, iSp, iSh, nOrb1, nOrb2, lShellVals(orb%mShell)
    real(dp) :: theta, rotZ(orb%mOrb,orb%mOrb)

    select case(this%iBoundaryCondition)

    case(boundaryConditions%cluster, boundaryConditions%pbc3d)

      ! orbitals all in the same orientation in space

      return

    case(boundaryConditions%helical)

      if (iAt1 == iAt2) then
        return
      end if
      iAt2f = img2CentCell(iAt2)
      iSp = species(iAt2f)
      iSh = orb%nShell(iSp)
      nOrb1 = orb%nOrbAtom(iAt1)
      nOrb2 = orb%nOrbAtom(iAt2f)
      lShellVals(:iSh) = orb%angShell(:iSh,iSp)
      theta = atan2(coords(2,iAt2),coords(1,iAt2)) - atan2(coords(2,iAt2f),coords(1,iAt2f))
      theta = mod(theta,2.0_dp*pi)
      call rotateZ(rotZ, lShellVals(:iSh), theta)
      sqrBlock(:nOrb2, :nOrb1) = matmul(rotZ(:nOrb2, :nOrb2), sqrBlock(:nOrb2, :nOrb1))

    end select

  end subroutine foldOutDiatomicBlock_${LABEL}$

#:endfor


  !> Transform alignment of vector array between alignment of coord0 and coord (central cell only)
  subroutine alignVectorCentralCell(this, vectors, coords, coords0, nAtom)

    !> Instance
    class(TBoundaryConditions), intent(in) :: this

    !> vectors to transform
    real(dp), intent(inout) :: vectors(:,:)

    !> Coordinates in central/objective cell
    real(dp), intent(in) :: coords(:,:)

    !> Coordinates in (potentially) extended structure
    real(dp), intent(in) :: coords0(:,:)

    !> Number of atoms in central/objective cell
    integer, intent(in) :: nAtom

    integer :: iAt
    real(dp) :: deltaTheta

    select case(this%iBoundaryCondition)

    case(boundaryConditions%helical)

      do iAt = 1, nAtom
        deltaTheta = atan2(coords0(2,iAt),coords0(1,iAt)) - atan2(coords(2,iAt),coords(1,iAt))
        call rotate3(vectors(:,iAt), deltaTheta, zAxis)
      end do

    case default

      return

    end select

  end subroutine alignVectorCentralCell


  !> Fold coordinates back in the central cell.
  !>
  !> Throw away the integer part of the relative coordinates of every atom. If the resulting
  !> coordinate is very near to 1.0 (closer than 1e-12 in absolute length), fold it to 0.0 to make
  !> the algorithm more predictable and independent of numerical noise.
  subroutine foldCoordsToCell(this, coord, latVec)

    !> Instance
    class(TBoundaryConditions), intent(in) :: this

    !> Contains the original coordinates on call and the folded ones on return.
    real(dp), intent(inout) :: coord(:,:)

    !> Lattice vectors (column format).
    real(dp), intent(in) :: latVec(:,:)

    integer :: nAtom, iAt, jj
    real(dp) :: frac(3), frac2(3), tmp3(3), vecLen(3), thetaNew, thetaOld, invLatVecs(3,3)

    nAtom = size(coord, dim=2)

    select case(this%iBoundaryCondition)

    case(boundaryConditions%cluster)

      ! No unit cell to fold into

      return

    case(boundaryConditions%pbc3d)

      call invert33(invLatVecs, latVec)
      vecLen(:) = sqrt(sum(latVec**2, dim=1))
      do iAt = 1, nAtom
        frac(:) = matmul(invLatVecs, coord(:,iAt))
        tmp3(:) = coord(:,iAt)
        frac2(:) = frac(:) - real(floor(frac(:)), dp)
        where (abs(vecLen*(1.0_dp - frac2)) < epsilon(0.0_dp)) frac2 = 0.0_dp
        coord(:, iAt) = matmul(latVec, frac2)
      end do

    case(boundaryConditions%helical)

      do iAt = 1, nAtom
        jj = floor(coord(3,iAt)/latVec(1,1))
        ! want coordinate in eventual range 0..latVec(1,1) hence floor
        tmp3(:) = coord(:,iAt)
        coord(3,iAt) = coord(3,iAt) - jj * latVec(1,1)
        call rotate3(coord(:,iAt), -jj*latVec(2,1),zAxis)
        thetaOld = atan2(coord(2,iAt), coord(1,iAt))
        thetaNew = mod(thetaOld+2.0_dp*pi, 2.0_dp*pi/latvec(3,1))
        call rotate3(coord(:,iAt),-thetaOld+thetaNew,zAxis)
      end do

    end select

  end subroutine foldCoordsToCell

end module dftbp_dftb_boundarycond
