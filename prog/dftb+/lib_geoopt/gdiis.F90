!--------------------------------------------------------------------------------------------------!
!  DFTB+: general package for performing fast atomistic simulations                                !
!  Copyright (C) 2017  DFTB+ developers group                                                      !
!                                                                                                  !
!  See the LICENSE file for terms of usage and distribution.                                       !
!--------------------------------------------------------------------------------------------------!

#:include 'common.fypp'

!!* Contains a geometry DIIS optimizer interface.
module gdiis
  use assert
  use accuracy
  use diismixer
  implicit none

  private

  !!* Contains data for the DIIS mimimizer
  type ODIIS
    private
    type(ODIISMixer) :: pDIIS
    real(dp), allocatable :: x(:)
    integer  :: nElem
    real(dp) :: tolerance    !* Tolerance criteria for convergence
    logical  :: tInitialized !* If object is initialized
  end type ODIIS

  !!* Creates gDIIS instance
  interface init
    module procedure gDIIS_init
  end interface

  !!* Resets the gDIIS instance
  interface reset
    module procedure gDIIS_reset
  end interface

  !!* Passes calculated gradient to the minimizer and returns a new point
  interface next
    module procedure gDIIS_next
  end interface

  public :: ODIIS
  public :: init, reset, next


contains

  subroutine gDIIS_init(self, nElem, tol, alpha, nGens)
    type(ODIIS), intent(out) :: self
    integer, intent(in)  :: nElem
    real(dp), intent(in) :: tol
    real(dp), intent(in) :: alpha
    integer, intent(in)  :: nGens

    self%nElem = nElem
    self%tolerance = tol
    allocate(self%x(self%nElem))
    call init(self%pDIIS,nGens,alpha,.true.,alpha)
    self%tInitialized = .true.

  end subroutine gDIIS_init


  subroutine gDIIS_reset(self,x)
    type(ODIIS), intent(inout) :: self
    real(dp) :: x(:)

    call reset(self%pDIIS, self%nElem)
    self%x(:) = x(:)

  end subroutine gDIIS_reset


  subroutine gDIIS_next(self,dx, xNew, tConverged)
    type(ODIIS), intent(inout) :: self
    real(dp), intent(in)  :: dx(:)
    real(dp), intent(out) :: xNew(:)
    logical,  intent(out) :: tConverged

    @:ASSERT(self%tInitialized)
    @:ASSERT(size(xNew) == self%nElem)
    @:ASSERT(size(dx) == self%nElem)

    xNew = self%x
    call mix(self%pDIIS,self%x,dx)
    if (maxval(abs(xNew-self%x)) < self%tolerance &
        & .or. (maxval(abs(dx)) < self%tolerance)) then
      tConverged = .true.
    else
      tConverged = .false.
    end if
    xNew = self%x

  end subroutine gDIIS_next

end module gdiis
