!--------------------------------------------------------------------------------------------------!
!  DFTB+: general package for performing fast atomistic simulations                                !
!  Copyright (C) 2017  DFTB+ developers group                                                      !
!                                                                                                  !
!  See the LICENSE file for terms of usage and distribution.                                       !
!--------------------------------------------------------------------------------------------------!

#:include 'common.fypp'

!!* Function minimization with the standard conjugate gradient technique.
!!* @see Numerical Recipes
module conjgrad
  use assert
  use accuracy
  use linemin
  implicit none

  private


  !!* Contains data for the conjugate gradient minimizer
  type OConjGrad
    private
    integer :: state                   !* State of the minimizer
    integer :: nElem                   !* Nr. of variables
    real(dp), allocatable :: gg(:)         !* Gradient in previous cycle
    real(dp), allocatable :: hh(:)         !* Conjugate gradient
    real(dp), allocatable :: uu(:)         !* Last calculated point
    real(dp) :: tolerance              !* Tolerance criteria for convergence
    real(dp) :: maxDisp                !* Maximal displacement along one
                                       !* coordinate in one step
    logical :: tConverged              !* If CG converged
    logical :: tInitialized            !* If object is initialized
    type(OLineMin) :: pLinMin  !* Line minimizer
  end type OConjGrad


  !!* Initialises CG instance
  interface init
    module procedure ConjGrad_init
  end interface

  !!* Resets CG
  interface reset
    module procedure ConjGrad_reset
  end interface

  !!* Passes calculated function value and gradient to the minimizer and gives
  !!* back a new coordinate.
  interface next
    module procedure ConjGrad_next
  end interface

  !!* Coordinate of the minimum.
  interface getMinX
    module procedure ConjGrad_getMinX
  end interface

  !!* Function value in the minimum.
  interface getMinY
    module procedure ConjGrad_getMinY
  end interface

  !!* Gradient in the minimum.
  interface getMinGrad
    module procedure ConjGrad_getMinGrad
  end interface


  public :: OConjGrad
  public :: init, reset, next, getMinX, getMinY, getMinGrad

  integer, parameter  :: st_1 = 1, st_2 = 2

contains

  !!* Creates a conjugate gradient instance
  !!* @param self    Conjugate gradient instance on exit
  !!* @param nElem   Nr. of elements in the vectors
  !!* @param tol     Tolerance for the gradient
  !!* @param maxDisp Maximal displacement in one element in one step
  subroutine ConjGrad_init(self, nElem, tol, maxDisp)
    type(OConjGrad), intent(out) :: self
    integer, intent(in) :: nElem
    real(dp), intent(in) :: tol
    real(dp), intent(in) :: maxDisp

    @:ASSERT(nElem > 0)
    @:ASSERT(tol > 0.0_dp)
    @:ASSERT(maxDisp > 0.0_dp)

    self%nElem = nElem
    self%tolerance = tol
    self%maxDisp = maxDisp
    allocate(self%gg(nElem))
    allocate(self%hh(nElem))
    allocate(self%uu(nElem))
    !! Line minimizer is created with an extrem big tolerance: it just brackets
    !! the minimum and returns an approximative minimum between them. Seems
    !! to give in most cases better results as making many line min. steps.
    call init(self%pLinMin, nElem, 10, 10000.0_dp, self%maxDisp)
    self%tInitialized = .false.

  end subroutine ConjGrad_init


  !!* Resets CG minimizer
  !!* @param self CG minimizer
  !!* @param x0   Point to start from
  subroutine ConjGrad_reset(self, x0)
    type(OConjGrad), intent(inout) :: self
    real(dp), intent(in) :: x0(:)

    @:ASSERT(size(x0) == self%nElem)

    self%uu(:) = x0(:)
    self%state = st_1
    self%tConverged = .false.
    self%tInitialized = .true.

  end subroutine ConjGrad_reset



  !!* Passes calculated function value and gradient to the minimizare and
  !!* gives a new coordinate back.
  !!* @param self       CG minimizer
  !!* @param fx         Function value for last point returned by this routine
  !!* @param dx         Gradient in the last point
  !!* @param xNew       New proposed point
  !!* @param tConverged True, if gradient got below the specified tolerance.
  !!* @note When calling the first time, funciton value and gradient for the
  !!*   starting point of the minimization should be passed.
  subroutine ConjGrad_next(self, fx, dx, xNew, tConverged)
    type(OConjGrad), intent(inout) :: self
    real(dp), intent(in)  :: fx
    real(dp), intent(in)  :: dx(:)
    real(dp), intent(out) :: xNew(:)
    logical,  intent(out) :: tConverged

    @:ASSERT(self%tInitialized)
    @:ASSERT(size(xNew) == self%nElem)
    @:ASSERT(size(dx) == self%nElem)

    if (.not. self%tConverged) then
      call next_local(self%state, self%gg, self%hh, self%uu, self%tConverged, &
          &self%tolerance, self%pLinMin, fx, dx)
    end if
    xNew(:) = self%uu(:)
    tConverged = self%tConverged
  end subroutine ConjGrad_next



  !!* Working horse for the CG minimizer
  !!* @param fu Function value in the last point
  !!* @param du Gradient in the last point
  !!* @param uu New calculated point
  !!* @note Most parameters are identical with the fields of the derived type
  !!*   for the CG minimizer, documentation is to be found there.
  subroutine next_local(state, gg, hh, uu, tConverged, tolerance, pLinMin, fu,&
      &du)
    integer, intent(inout) :: state
    real(dp), intent(inout) :: gg(:)
    real(dp), intent(inout) :: hh(:)
    real(dp), intent(inout) :: uu(:)
    logical,  intent(inout) :: tConverged
    real(dp), intent(in)    :: tolerance
    type(OLineMin), intent(inout) :: pLinMin
    real(dp), intent(in) :: fu
    real(dp), intent(in) :: du(:)

    real(dp) :: ggAbs, dgAbs, rTmp
    logical :: tConvLine
    real(dp), allocatable  :: xi(:)

    if (state == st_1) then
      !! If first gradient converged: Reuse internal variables to store results
      if (maxval(abs(du)) < tolerance) then
        tConverged = .true.
        gg(:) = du(:)
        hh(:) = 0.0_dp
        hh(1) = fu
      else
        gg(:) = -du(:)
        hh(:) = gg(:)
      end if
      !! First step x = F/k, where F is the acting force and
      !! k is a spring constant in the magnitude of a C-C vibration mode
      call reset(pLinMin, uu, hh, 5.0_dp * sqrt(sum(gg**2)))
      state = st_2
    end if

    tConvLine = .true.
    do while ((.not. tConverged) .and. tConvLine)
      call next(pLinMin, fu, du, uu, tConvLine)
      if (tConvLine) then
        allocate(xi(size(gg)))
        call getMinGrad(pLinMin, xi)
        if (maxval(abs(xi)) < tolerance) then
          tConverged = .true.
        else
          ggAbs = dot_product(gg, gg)
          dgAbs = dot_product(xi+gg, xi)
          if (ggAbs < epsilon(0.0_dp)) then
            tConverged = .true.
          else
            gg(:) = -xi(:)
            hh(:) = gg(:) + (dgAbs / ggAbs) * hh(:)
            call getMinX(pLinMin, uu)
            call getMinLambda(pLinMin, rTmp)
            call reset(pLinMin, uu, hh, rTmp)
          end if
        end if
      end if
    end do

    !! If converged, reuse internal variables to store the result
    if (tConverged .and. state /= st_1) then
      call getMinX(pLinMin, uu)
      call getMinGrad(pLinMin, gg)
      call getMinY(pLinMin, hh(1))
    end if

  end subroutine next_local



  !!* Gives the coordinate of the minimal point back
  !!* @param self CG minimizer
  !!* @param minX Coordinate of the minimal point
  !!* @note The passed back value is meaningless if the subroutine is called
  !!*   before the CG minimizer signalizes convergence.
  subroutine ConjGrad_getMinX(self, minX)
    type(OConjGrad), intent(in) :: self
    real(dp), intent(out) :: minX(:)

    @:ASSERT(self%tInitialized .and. self%tConverged)
    @:ASSERT(size(minX) == self%nElem)
    call getMinX(self%pLinMin, minX)

  end subroutine ConjGrad_getMinX



  !!* Gives the function value in the minimal point back
  !!* @param self CG minimizer
  !!* @param minY Coordinate of the minimal point
  !!* @note The passed back value is meaningless if the subroutine is called
  !!*   before the CG minimizer signalizes convergence.
  subroutine ConjGrad_getMinY(self, minY)
    type(OConjGrad), intent(in) :: self
    real(dp), intent(out) :: minY

    @:ASSERT(self%tInitialized .and. self%tConverged)
    call getMinY(self%pLinMin, minY)

  end subroutine ConjGrad_getMinY



  !!* Gives the gradient in the minimal point back
  !!* @param self    CG minimizer
  !!* @param minGrad Coordinate of the minimal point
  !!* @note The passed back value is meaningless if the subroutine is called
  !!*   before the CG minimizer signalizes convergence.
  subroutine ConjGrad_getMinGrad(self, minGrad)
    type(OConjGrad), intent(in) :: self
    real(dp), intent(out) :: minGrad(:)

    @:ASSERT(self%tInitialized .and. self%tConverged)
    @:ASSERT(size(minGrad) == self%nElem)
    call getMinGrad(self%pLinMin, minGrad)

  end subroutine ConjGrad_getMinGrad


end module conjgrad
