!--------------------------------------------------------------------------------------------------!
!  DFTB+: general package for performing fast atomistic simulations                                !
!  Copyright (C) 2006 - 2025  DFTB+ developers group                                               !
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
  use dftbp_math_matrixops, only : pseudoInv
  use dftbp_math_quaternions, only : rotate3
  use dftbp_math_simplealgebra, only : invert33, determinant33
  use dftbp_type_commontypes, only : TOrbitals
  implicit none

  private
  public :: zAxis, boundaryCondsEnum, TBoundaryConds, TBoundaryConds_init

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
  type(TBoundaryConditionEnum_), parameter :: boundaryCondsEnum = TBoundaryConditionEnum_()


  !> Boundary condition type with information on geometry and methods for applying geometric
  !! boundary conditions
  type TBoundaryConds

    !> Which type of boundary conditions are present
    integer :: iBoundaryCondition = boundaryCondsEnum%unknown

    !> Is this a transport calculation with semi-infinite contacts?
    logical :: areSIContactsPresent

    !> Cartesian/lattice vector indices unaffected by transport contacts, hence actually periodic if
    !> suitable boundary conditions
    integer, allocatable :: nonContactDirs(:)

  contains

    !> Map diatomic hamiltonian/overlap block in extended structure into central cell
    procedure :: foldInDiatomicBlock_real
    procedure :: foldInDiatomicBlock_cplx
    generic :: foldInDiatomicBlock => foldInDiatomicBlock_real, foldInDiatomicBlock_cplx

    !> Map diatomic hamiltonian/overlap block in central cell out into extended structure
    procedure :: foldOutDiatomicBlock_real
    procedure :: foldOutDiatomicBlock_cplx
    generic :: foldOutDiatomicBlock => foldOutDiatomicBlock_real, foldOutDiatomicBlock_cplx

    !> Map vectors from alignment of central and extended coordinates in the central cell
    procedure :: alignVectorCentralCell

    !> Folds geometry into central cell
    procedure :: foldCoordsToCell

    !> Setup for any semi-infinite contacts that could interact with boundary conditions
    procedure :: setTransportDirections

    !> Update quantities dependent on boundary conditions, based on central cell definition
    procedure :: handleBoundaryChanges

  end type TBoundaryConds


contains


  !> Initialize instance of the boundary conditions type
  subroutine TBoundaryConds_init(this, iBoundaryCondition, errStatus)

    !> Instance
    type(TBoundaryConds), intent(out) :: this

    !> Boundary condition choice
    integer, intent(in) :: iBoundaryCondition

    !> Status of routine
    type(TStatus), intent(out) :: errStatus

    if (.not. any([boundaryCondsEnum%cluster, boundaryCondsEnum%pbc3d,&
        & boundaryCondsEnum%helical] == iBoundaryCondition)) then
      @:RAISE_ERROR(errStatus, -1, "Unknown boundary condition specified")
    end if

    this%iBoundaryCondition = iBoundaryCondition

    this%areSIContactsPresent = .false.

  end subroutine TBoundaryConds_init


  !> In the case of calculations with semi-infinite contacts, set the directions which are impacted
  !> by contacts for things like periodic folding.
  subroutine setTransportDirections(this, contactDirs)

    !> Instance
    class(TBoundaryConds), intent(inout) :: this

    !> Which, if any, of the three directions (either cartesian, if molecular, or lattice vectors if
    !> periodic) are impacted by presence of semi-infinite contacts. Elements are true if there is a
    !> contact impinging on that direction.
    logical, intent(in) :: contactDirs(:)

    integer :: iCart, iDir

    @:ASSERT(this%iBoundaryCondition == boundaryCondsEnum%pbc3d)
    @:ASSERT(size(contactDirs) >= 1 .and. size(contactDirs) <= 3)

    allocate(this%nonContactDirs(count(.not.contactDirs)))
    this%areSIContactsPresent = .false.
    iDir = 0
    do iCart = 1, 3
      if (.not.contactDirs(iCart)) then
        iDir = iDir + 1
        this%nonContactDirs(iDir) = iCart
      else
        this%areSIContactsPresent = .true.
      end if
    end do

  end subroutine setTransportDirections


#:for TYPE, LABEL in FLAVOURS

  !> Convert matrix elements from image to central cell atoms
  pure subroutine foldInDiatomicBlock_${LABEL}$(this, sqrBlock, iAt1, iAt2, coords, species,&
      & img2CentCell, orb)

    !> Instance
    class(TBoundaryConds), intent(in) :: this

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

    !> Data type for atomic orbital information
    type(TOrbitals), intent(in) :: orb

    integer :: iAt2f, iSp, iSh, nOrb1, nOrb2, lShellVals(orb%mShell)
    real(dp) :: theta, rotZ(orb%mOrb,orb%mOrb)

    select case(this%iBoundaryCondition)

    case(boundaryCondsEnum%cluster, boundaryCondsEnum%pbc3d)

      ! orbitals all in the same orientation in space

      return

    case(boundaryCondsEnum%helical)

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


  !> Convert matrix elements from central cell to image atoms
  pure subroutine foldOutDiatomicBlock_${LABEL}$(this, sqrBlock, iAt1, iAt2, coords, species,&
      & img2CentCell, orb)

    !> Instance
    class(TBoundaryConds), intent(in) :: this

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

    !> Data type for atomic orbital information
    type(TOrbitals), intent(in) :: orb

    integer :: iAt2f, iSp, iSh, nOrb1, nOrb2, lShellVals(orb%mShell)
    real(dp) :: theta, rotZ(orb%mOrb,orb%mOrb)

    select case(this%iBoundaryCondition)

    case(boundaryCondsEnum%cluster, boundaryCondsEnum%pbc3d)

      ! orbitals all in the same orientation in space

      return

    case(boundaryCondsEnum%helical)

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
    class(TBoundaryConds), intent(in) :: this

    !> Vectors to transform
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

    case(boundaryCondsEnum%helical)

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
    class(TBoundaryConds), intent(in) :: this

    !> Contains the original coordinates on call and the folded ones on return.
    real(dp), intent(inout) :: coord(:,:)

    !> Lattice vectors (column format).
    real(dp), intent(in) :: latVec(:,:)

    integer :: nAtom, iAt, jj, nDir
    real(dp) :: frac(3), tmp3(3), vecLen(3), thetaNew, thetaOld, invLatVecs(3,3)
    real(dp), allocatable :: work(:,:)

    nAtom = size(coord, dim=2)

    select case(this%iBoundaryCondition)

    case(boundaryCondsEnum%cluster)

      ! No unit cell to fold into

      return

    case(boundaryCondsEnum%pbc3d)

      if (this%areSIContactsPresent) then

        nDir = size(this%nonContactDirs)
        work = latVec(:,this%nonContactDirs)
        vecLen(:nDir) = sqrt(sum(work**2, dim=1))
        call pseudoInv(work, invLatVecs(:,:nDir))
        !$OMP PARALLEL DO&
        !$OMP& DEFAULT(SHARED) PRIVATE(frac, tmp3) SCHEDULE(RUNTIME)
        do iAt = 1, nAtom
          frac(:nDir) = matmul(transpose(invLatVecs(:,:nDir)), coord(:,iAt))
          tmp3(:) = coord(:,iAt) - matmul(latVec(:,:nDir), frac(:nDir))
          frac(:nDir) = frac(:nDir) - real(int(frac(:nDir)), dp)
          where (abs(vecLen(:nDir)*(1.0_dp - frac(:nDir))) < epsilon(0.0_dp)) frac(:nDir) = 0.0_dp
          coord(:, iAt) = tmp3 + matmul(latVec(:,:nDir), frac(:nDir))
        end do
        !$OMP END PARALLEL DO

      else

        call invert33(invLatVecs, latVec)
        vecLen(:) = sqrt(sum(latVec**2, dim=1))
        !$OMP PARALLEL DO&
        !$OMP& DEFAULT(SHARED) PRIVATE(frac) SCHEDULE(RUNTIME)
        do iAt = 1, nAtom
          frac(:) = matmul(invLatVecs, coord(:,iAt))
          frac(:) = frac - real(int(frac), dp)
          where (abs(vecLen*(1.0_dp - frac)) < epsilon(0.0_dp)) frac = 0.0_dp
          coord(:, iAt) = matmul(latVec, frac)
        end do
        !$OMP END PARALLEL DO

      end if

    case(boundaryCondsEnum%helical)

      do iAt = 1, nAtom
        jj = floor(coord(3,iAt)/latVec(1,1))
        ! want coordinate in eventual range 0..latVec(1,1) hence floor
        coord(3,iAt) = coord(3,iAt) - jj * latVec(1,1)
        call rotate3(coord(:,iAt), -jj*latVec(2,1),zAxis)
        thetaOld = atan2(coord(2,iAt), coord(1,iAt))
        thetaNew = mod(thetaOld+2.0_dp*pi, 2.0_dp*pi/latvec(3,1))
        call rotate3(coord(:,iAt),-thetaOld+thetaNew,zAxis)
      end do

    end select

  end subroutine foldCoordsToCell


  !> Updates quantities, like reciprocal lattice, that are dependant on the central cell boundary
  !> condition definitions
  subroutine handleBoundaryChanges(this, latVecs, invLatVecs, recVecs, cellVol, recCellVol)

    !> Instance
    class(TBoundaryConds), intent(in) :: this

    !> New lattice vector data
    real(dp), intent(in) :: latVecs(:,:)

    !> Inverse of lattice vectors / reciprocal lattice vectors in units of 2 pi
    real(dp), intent(out) :: invLatVecs(:,:)

    !> Reciprocal lattice vectors
    real(dp), intent(out) :: recVecs(:,:)

    !> Unit cell volume if relevant
    real(dp), intent(out) :: cellVol

    !> Reciprocal lattice unit cell volume
    real(dp), intent(out) :: recCellVol

    cellVol = 0.0_dp
    recCellVol = 0.0_dp
    select case(this%iBoundaryCondition)
    case(boundaryCondsEnum%pbc3d)
      if (this%areSIContactsPresent) then
        ! semi-infinite contacts block periodicity in some direction(s)
        invLatVecs(:,:size(this%nonContactDirs)) = latVecs(:, this%nonContactDirs)
        call pseudoInv(invLatVecs(:,:size(this%nonContactDirs)),&
            & recVecs(:,:size(this%nonContactDirs)))
        invLatVecs(:,:) = 0.0_dp
        invLatVecs(:, this%nonContactDirs) = recVecs(:,:size(this%nonContactDirs))
        recVecs(:,:) = 2.0_dp * pi * invLatVecs
      else
        invLatVecs(:,:) = latVecs
        call invert33(invLatVecs)
        invLatVecs(:,:) = transpose(invLatVecs)
        recVecs(:,:) = 2.0_dp * pi * invLatVecs
        cellVol = abs(determinant33(latVecs))
        recCellVol = abs(determinant33(recVecs))
      end if
    end select

  end subroutine handleBoundaryChanges

end module dftbp_dftb_boundarycond
