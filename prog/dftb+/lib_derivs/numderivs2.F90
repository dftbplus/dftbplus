!--------------------------------------------------------------------------------------------------!
!  DFTB+: general package for performing fast atomistic simulations                                !
!  Copyright (C) 2006 - 2020  DFTB+ developers group                                               !
!                                                                                                  !
!  See the LICENSE file for terms of usage and distribution.                                       !
!--------------------------------------------------------------------------------------------------!

#:include 'common.fypp'

!> Module containing routines for numerical second derivs of energy using central finite difference.
!> To Do: Option to restart the calculation
module dftbp_numderivs2
  use dftbp_assert
  use dftbp_accuracy, only : dp
  implicit none
  private

  public :: TNumDerivs, create, next, getHessianMatrix


  !> Contains necessary data for the derivs
  type TNumDerivs
    private

    !> Internal matrix to hold derivative and intermediate values for their construction
    real(dp), allocatable :: derivs(:,:)


    !> Coordinates at x=0 to differentiate at
    real(dp), allocatable :: x0(:,:)


    !> How many derivates are needed
    integer :: nDerivs


    !> Which atom are we currently differentiating with respect to?
    integer :: iAtom


    !> Which component, x,y,z are we currently differentiating with respect to?
    integer :: iComponent


    !> displacement along + or - for central difference
    real(dp) :: iDelta


    !> Step size for derivative
    real(dp) :: Delta
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
  !> Note: Use pre-relaxed coordinates when starting this, as the the truncation at second
  !> derivatives is only valid at the minimum position.
  subroutine derivs_create(self,xInit,Delta)

    !> Pointer to the initialised object on exit.
    type(TNumDerivs), allocatable, intent(out) :: self

    !> initial atomic coordinates (3,:)
    real(dp), intent(inout) :: xInit(:,:)

    !> step size for numerical derivative
    real(dp), intent(in) :: Delta

    integer :: nDerivs

    @:ASSERT(size(xInit,dim=1)==3)
    nDerivs = size(xInit,dim=2)

    allocate(self)
    allocate(self%x0(3, nDerivs))
    self%x0(:,:) = xInit(:,:)
    allocate(self%derivs(3*nDerivs,3*nDerivs))
    self%derivs(:,:) = 0.0_dp
    self%nDerivs = nDerivs
    self%Delta = Delta

    self%iAtom = 1
    self%iComponent = 1
    self%iDelta = -1.0_dp

    xInit(self%iComponent,self%iAtom) = &
        & xInit(self%iComponent,self%iAtom) + self%iDelta*self%Delta

  end subroutine derivs_create


  !> Takes the next step for derivatives using the central difference formula to choose the new
  !> coordinates for differentiation of the forces with respect to atomic coordinates
  subroutine derivs_next(self,xNew,fOld,tGeomEnd)

    !> Derivatives instance to propogate
    type(TNumDerivs), intent(inout) :: self

    !> New coordinates for the next step
    real(dp), intent(out) :: xNew(:,:)

    !> Forces for the previous geometry
    real(dp), intent(in) :: fOld(:,:)

    !> Has the process terminated? If so internally calculate the Hessian matrix.
    logical, intent(out) :: tGeomEnd

    integer :: ii, jj

    @:ASSERT(all(shape(xNew)==shape(fOld)))
    @:ASSERT(all(shape(xNew)==(/3,self%nDerivs/)))

    if (self%iAtom==self%nDerivs .and. self%iComponent == 3 .and. &
        & self%iDelta > 0.0_dp) then
      tGeomEnd = .true.
    else
      tGeomEnd = .false.
    end if

    do ii = 1, self%nDerivs
      do jj = 1, 3
        self%derivs((ii-1)*3+jj,(self%iAtom-1)*3+self%iComponent) = &
            & self%derivs((ii-1)*3+jj,(self%iAtom-1)*3+self%iComponent) &
            & + self%iDelta * fOld(jj,ii)
      end do
    end do

    if (.not.tGeomEnd) then

      if (self%iDelta < 0.0_dp) then
        self%iDelta = 1.0_dp
      else
        self%iDelta = -1.0_dp
        if (self%iComponent == 3) self%iAtom = self%iAtom + 1
        self%iComponent = mod(self%iComponent,3) + 1
      end if

      xNew(:,:) = self%x0(:,:)
      xNew(self%iComponent,self%iAtom) = xNew(self%iComponent,self%iAtom) + &
          & self%iDelta * self%Delta
    else
      ! get actual derivatives
      self%derivs(:,:) = 0.5_dp*self%derivs(:,:)/(self%Delta)
      ! set xnew to an arbitrary value
      xNew(:,:) = self%x0
    end if

  end subroutine derivs_next


  !> Routine to return pointer to internal matrix of derivative elements.
  subroutine getDerivMatrixPtr(self,d)

    !> Derivatives instance including the Hessian internally
    type(TNumDerivs), intent(in), target :: self

    !> Pointer to the Hessian matrix to allow retrieval
    real(dp), pointer :: d(:,:)

    d => self%derivs

  end subroutine getDerivMatrixPtr

end module dftbp_numderivs2
