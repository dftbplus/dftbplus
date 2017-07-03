!--------------------------------------------------------------------------------------------------!
!  DFTB+: general package for performing fast atomistic simulations                                !
!  Copyright (C) 2017  DFTB+ developers group                                                      !
!                                                                                                  !
!  See the LICENSE file for terms of usage and distribution.                                       !
!--------------------------------------------------------------------------------------------------!

#:include 'common.fypp'

!!* Module containing routines for numerical second derivs of energy using
!!* central finite difference.
!!* @todo Add option to restart the calculation
module numderivs2
  use assert
  use accuracy, only : dp
  implicit none
  private

  public :: OnumDerivs, create, next, getHessianMatrix

  !!* Contains necessary data for the derivs
  type OnumDerivs
    private
    ! Internal matrix to hold derivative and
    ! intermediate values for their construction
    real(dp), allocatable :: derivs(:,:)

    ! Coordinates at x=0 to differentiate at
    real(dp), allocatable :: x0(:,:)

    integer           :: nDerivs     !!* How many derivates are needed
    integer           :: iAtom       !!* Which atom are we currently
    !!* differentiating with respect to?
    integer           :: iComponent  !!* Which component, x,y,z are we currently
    !!* differentiating with respect to?
    real(dp)          :: iDelta      !!* displacement along + or - for central
    !!* difference
    real(dp)          :: Delta       !!* Step size for derivative
  end type OnumDerivs

  !!* Create numerical second derivatives instance
  interface create
    module procedure derivs_create
  end interface

  !!* Delivers the next set of coordinates for evaluation of forces
  interface next
    module procedure derivs_next
  end interface

  !!* Get pointer to the Hessian matrix of derivatives for the system
  interface getHessianMatrix
    module procedure getDerivMatrixPtr
  end interface

contains

  !!* Create new instance of derivative object
  !!* @param self Pointer to the initialised object on exit.
  !!* @param xInit initial atomic coordinates (3,:)
  !!* @param Delta step size for numerical derivative
  !!* @note Use pre-relaxed coordinates when starting this, as the the
  !!* truncation at second derivatives is only valid at the minimum position.
  subroutine derivs_create(self,xInit,Delta)
    type(OnumDerivs), allocatable, intent(out) :: self
    real(dp), intent(inout) :: xInit(:,:)
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


  !!* Takes the next step for derivatives using the central difference
  !!* formula to choose the new coordinates for differentiation of the
  !!* forces with respect to atomic coordinates
  !!* @param self Derivatives instance to propogate
  !!* @param xNew New coordinates for the next step
  !!* @param fOld Forces for the previous geometry
  !!* @param tGeomEnd Has the process terminated? If so internally calculate
  !!* the Hessian matrix.
  subroutine derivs_next(self,xNew,fOld,tGeomEnd)
    type(OnumDerivs), intent(inout) :: self
    real(dp), intent(out)     :: xNew(:,:)
    real(dp), intent(in)      :: fOld(:,:)
    logical, intent(out)      :: tGeomEnd

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

  !!* Routine to return pointer to internal matrix of derivative elements.
  !!* @param self Derivatives instance including the Hessian internally
  !!* @param d Pointer to the Hessian matrix to allow retrieval
  subroutine getDerivMatrixPtr(self,d)
    type(OnumDerivs), intent(in), target :: self
    real(dp), pointer :: d(:,:)

    d => self%derivs

  end subroutine getDerivMatrixPtr

end module numderivs2
