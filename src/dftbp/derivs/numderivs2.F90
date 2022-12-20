!--------------------------------------------------------------------------------------------------!
!  DFTB+: general package for performing fast atomistic simulations                                !
!  Copyright (C) 2006 - 2022  DFTB+ developers group                                               !
!                                                                                                  !
!  See the LICENSE file for terms of usage and distribution.                                       !
!--------------------------------------------------------------------------------------------------!

#:include 'common.fypp'

!> Module containing routines for numerical second derivs of energy using central finite difference.
!> To Do: Option to restart the calculation
module dftbp_derivs_numderivs2
  use dftbp_common_accuracy, only : dp
  implicit none

  private
  public :: TNumDerivs, create, next, getHessianMatrix, dipoleAdd, polAdd


  !> Contains necessary data for the derivs
  type :: TNumDerivs
    private

    !> Internal matrix to hold derivative and intermediate values for their construction
    !>
    !> Must be pointer, so that the type can safely return a pointer to it.
    !>
    real(dp), pointer :: forceDerivs(:,:) => null()

    !> Coordinates at x=0 to differentiate at
    real(dp), allocatable :: x0(:,:)

    !> How many derivates are moved
    integer :: nMovedAtoms

    !> Which atom are we currently differentiating with respect to?
    integer :: iAtom

    !> Which component, x,y,z are we currently differentiating with respect to?
    integer :: iComponent

    !> displacement along + or - for central difference
    real(dp) :: iDelta

    !> Step size for derivative
    real(dp) :: delta

    !> Dipole moment derivatives
    real(dp), pointer :: dipoleDerivs(:,:) => null()

    !> Polarizability derivatives
    real(dp), pointer :: polDerivs(:,:,:) => null()

  contains

    final :: TNumDerivs_final

  end type TNumDerivs


  !> Create numerical second derivatives instance
  interface create
    module procedure derivs_create
  end interface


  !> Delivers the next set of coordinates for evaluation of forces
  interface next
    module procedure derivs_next
  end interface


  !> Get the Hessian matrix of derivatives for the system
  interface getHessianMatrix
    module procedure getDerivMatrixPtr
  end interface

contains


  !> Create new instance of derivative object
  !> Note: Use pre-relaxed coordinates when starting this, as the truncation at second
  !> derivatives is only valid at the minimum position.
  !> The subroutine can allocate a rectangular matrix with parameter nDerivAtoms,
  !> Useful for distributed calculations of the Hessian
  subroutine derivs_create(this, xInit, nDerivAtoms, delta, isDipoleDiff, isPolDiff)

    !> Pointer to the initialised object on exit.
    type(TNumDerivs), allocatable, intent(out) :: this

    !> initial atomic coordinates (3, nMovedAtom)
    real(dp), intent(inout) :: xInit(:,:)

    !> number of atoms for which derivatives should be calculated (>= nMovedAtom)
    integer, intent(in) :: nDerivAtoms

    !> step size for numerical derivative
    real(dp), intent(in) :: delta

    !> Are dipole derivatives accumulated
    logical, intent(in), optional :: isDipoleDiff

    !> Are polarisability derivatives accumulated
    logical, intent(in), optional :: isPolDiff

    integer :: nMovedAtoms

    @:ASSERT(size(xInit,dim=1)==3)
    nMovedAtoms = size(xInit,dim=2)

    allocate(this)
    allocate(this%x0(3, nMovedAtoms))
    this%x0(:,:) = xInit(:,:)
    allocate(this%forceDerivs(3 * nDerivAtoms, 3 * nMovedAtoms), source=0.0_dp)
    if (present(isDipoleDiff)) then
      if (isDipoleDiff) then
        allocate(this%dipoleDerivs(3, 3 * nMovedAtoms), source=0.0_dp)
      end if
    end if
    if (present(isPolDiff)) then
      if (isPolDiff) then
        allocate(this%polDerivs(3, 3, 3 * nMovedAtoms), source=0.0_dp)
      end if
    end if
    this%nMovedAtoms = nMovedAtoms
    this%delta = delta

    this%iAtom = 1
    this%iComponent = 1
    this%iDelta = -1.0_dp

    xInit(this%iComponent,this%iAtom) = &
        & xInit(this%iComponent,this%iAtom) + this%iDelta*this%delta

  end subroutine derivs_create


  !> Takes the next step for derivatives using the central difference formula to choose the new
  !> coordinates for differentiation of the forces with respect to atomic coordinates
  subroutine derivs_next(this, xNew, fOld, tGeomEnd)

    !> Derivatives instance to propagate
    type(TNumDerivs), intent(inout) :: this

    !> New coordinates for the next step
    real(dp), intent(out) :: xNew(:,:)

    !> Forces for the previous geometry
    real(dp), intent(in) :: fOld(:,:)

    !> Has the process terminated? If so internally calculate the Hessian matrix.
    logical, intent(out) :: tGeomEnd

    integer :: ii, jj, nDerivAtoms

    @:ASSERT(all(shape(xNew)==(/3,this%nMovedAtoms/)))
    nDerivAtoms = size(this%forceDerivs, dim=1)/3
    @:ASSERT(size(fOld,1)==3)
    @:ASSERT(size(fOld,2)==nDerivAtoms)

    tGeomEnd = (this%iAtom == this%nMovedAtoms .and. this%iComponent == 3&
        & .and. this%iDelta > 0.0_dp)

    do ii = 1, nDerivAtoms
      do jj = 1, 3
        this%forceDerivs((ii-1)*3+jj,(this%iAtom-1)*3+this%iComponent) = &
            & this%forceDerivs((ii-1)*3+jj,(this%iAtom-1)*3+this%iComponent) &
            & + this%iDelta * fOld(jj,ii)
      end do
    end do

    if (.not.tGeomEnd) then

      if (this%iDelta < 0.0_dp) then
        this%iDelta = 1.0_dp
      else
        this%iDelta = -1.0_dp
        if (this%iComponent == 3) then
          this%iAtom = this%iAtom + 1
        end if
        this%iComponent = mod(this%iComponent,3) + 1
      end if

      xNew(:,:) = this%x0(:,:)
      xNew(this%iComponent,this%iAtom) = xNew(this%iComponent,this%iAtom) + &
          & this%iDelta * this%delta
    else
      ! assemble actual derivatives
      this%forceDerivs(:,:) = 0.5_dp * this%forceDerivs / this%delta
      if (associated(this%dipoleDerivs)) then
        this%dipoleDerivs(:, :) = 0.5_dp * this%dipoleDerivs / this%delta
      end if
      if (associated(this%polDerivs)) then
        this%polDerivs(:, :, :) = 0.5_dp * this%polDerivs / this%delta
      end if
      ! set xnew to an arbitrary value
      xNew(:,:) = this%x0
    end if

  end subroutine derivs_next


  !> Append dipole data
  subroutine dipoleAdd(this, dipole)

    !> Derivatives instance to propagate
    type(TNumDerivs), intent(inout) :: this

    !> Dipole moment
    real(dp), intent(in) :: dipole(3)

    integer :: jj

    do jj = 1, 3
      this%dipoleDerivs(jj,(this%iAtom-1)*3+this%iComponent) = &
          & this%dipoleDerivs(jj,(this%iAtom-1)*3+this%iComponent) + this%iDelta * dipole(jj)
    end do

  end subroutine dipoleAdd


  !> Append polarisation data
  subroutine polAdd(this, pol)

    !> Derivatives instance to propagate
    type(TNumDerivs), intent(inout) :: this

    !> Dipole moment
    real(dp), intent(in) :: pol(3,3)

    integer :: ii, jj

    do jj = 1, 3
      do ii = 1, 3
        this%polDerivs(ii,jj,(this%iAtom-1)*3+this%iComponent) = &
            & this%polDerivs(ii,jj,(this%iAtom-1)*3+this%iComponent) + this%iDelta * pol(ii,jj)
      end do
    end do

  end subroutine polAdd


  !> Routine to return pointer to internal matrix of derivative elements.
  subroutine getDerivMatrixPtr(this, d2, dip, pol)

    !> Derivatives instance including the Hessian internally
    type(TNumDerivs), intent(in) :: this

    !> Pointer to the Hessian matrix to allow retrieval
    real(dp), pointer, intent(out) :: d2(:, :)

    !> Pointer to the dipole derivative matrix to allow retrieval
    real(dp), pointer, intent(out), optional :: dip(:, :)

    !> Pointer to the polarisability derivative matrix to allow retrieval
    real(dp), pointer, intent(out), optional :: pol(:, :, :)

    d2 => this%forceDerivs

    if (present(dip)) then
      dip => this%dipoleDerivs
    end if

    if (present(pol)) then
      pol => this%polDerivs
    end if

  end subroutine getDerivMatrixPtr


  !> Finalizes TNumDerivs instance.
  subroutine TNumDerivs_final(this)

    !> Instance
    type(TNumDerivs), intent(inout) :: this

    if (associated(this%forceDerivs)) then
      deallocate(this%forceDerivs)
    end if
    if (associated(this%dipoleDerivs)) then
      deallocate(this%dipoleDerivs)
    end if
    if (associated(this%polDerivs)) then
      deallocate(this%polDerivs)
    end if

  end subroutine TNumDerivs_final

end module dftbp_derivs_numderivs2
