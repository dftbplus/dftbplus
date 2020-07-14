!--------------------------------------------------------------------------------------------------!
!  DFTB+: general package for performing fast atomistic simulations                                !
!  Copyright (C) 2006 - 2020  DFTB+ developers group                                               !
!                                                                                                  !
!  See the LICENSE file for terms of usage and distribution.                                       !
!--------------------------------------------------------------------------------------------------!

#:include 'common.fypp'

!> Extended Lagrangian dynamics
module dftbp_extlagrangian
  use dftbp_assert
  use dftbp_accuracy, only : dp
  use dftbp_message
  implicit none
  private

  public ::ExtLagrangian, ExtLagrangian_init
  public :: ExtLagrangianInp


  !> Input for an extended Lagrangian integrator.
  !!
  type :: ExtLagrangianInp


    !> Nr. of timesteps to consider for the integration.
    integer :: nTimeSteps


    !> Nr. of elements in the vectors which should be propagated.
    integer :: nElems

  end type ExtLagrangianInp


  !> Represents an extended Lagrangian integrator.
  type :: ExtLagrangian
    private
    integer :: phase
    integer :: nElems
    integer :: nTimeSteps
    integer :: iStep, nSteps, nTransientSteps
    integer :: ind
    real(dp) :: scale
    real(dp), allocatable :: precondMtx(:,:)
    real(dp), allocatable :: auxVectors(:,:)
    real(dp), allocatable :: auxCoeffs(:)
    real(dp) :: alpha, kappa
  contains

    !> turn on the integrator
    procedure :: turnOn

    !> input for next time steo
    procedure :: getNextInput

    !> are converged values required for next step
    procedure :: needsConvergedValues

    !> Preconditioner for integrator
    procedure :: setPreconditioner

    !> Internal helper function
    procedure, private :: updatePhaseAndSteps
  end type ExtLagrangian


  !> Internal type for enumerating different states of the integrator.
  type :: ExtLagrangianPhases
    integer :: off
    integer :: fillingUp
    integer :: interpolating
    integer :: on
  end type ExtLagrangianPhases


  !> Internal states of the integrator
  type(ExtLagrangianPhases), parameter :: phases =&
      & ExtLagrangianPhases(1, 2, 3, 4)


  !> Various integrator parameters for integration with 5 time steps
  real(dp), parameter :: auxCoeffs5(0:5) = &
      & [-6.0_dp, 14.0_dp, -8.0_dp, -3.0_dp, 4.0_dp, -1.0_dp]
  real(dp), parameter :: alpha5 = 18e-3_dp
  real(dp), parameter :: kappa5 = 1.82_dp


  !> Various integrator parameters for integration with 6 time steps
  real(dp), parameter :: auxCoeffs6(0:6) = &
      & [-14.0_dp, 36.0_dp, -27.0_dp, -2.0_dp, 12.0_dp, -6.0_dp, 1.0_dp]
  real(dp), parameter :: alpha6 = 5.5e-3_dp
  real(dp), parameter :: kappa6 = 1.84_dp


  !> Various integrator parameters for integration with 7 time steps
  real(dp), parameter :: auxCoeffs7(0:7) = &
      & [-36.0_dp, 99.0_dp, -88.0_dp, 11.0_dp, 32.0_dp, -25.0_dp, 8.0_dp,&
      & -1.0_dp]
  real(dp), parameter :: alpha7 = 1.6e-3
  real(dp), parameter :: kappa7 = 1.86_dp

contains


  !> Initializes an extended Lagrangian integrator.
  subroutine ExtLagrangian_init(this, input)


    !> Initialized instance at exit.
    class(ExtLagrangian), intent(inout) :: this


    !> Input container.
    class(ExtLagrangianInp), intent(in) :: input

    this%nElems = input%nElems
    this%scale = 1.0_dp
    this%phase = phases%off
    this%iStep = 1
    this%nSteps = 0
    this%nTransientSteps = 0

    this%nTimeSteps = input%nTimeSteps
    allocate(this%auxCoeffs(0:this%nTimeSteps))
    allocate(this%auxVectors(this%nElems, this%nTimeSteps + 1))
    select case (this%nTimeSteps)
    case (5)
      this%auxCoeffs(:) = auxCoeffs5
      this%alpha = alpha5
      this%kappa = kappa5
    case (6)
      this%auxCoeffs(:) = auxCoeffs6
      this%alpha = alpha6
      this%kappa = kappa6
    case(7)
      this%auxCoeffs(:) = auxCoeffs7
      this%alpha = alpha7
      this%kappa = kappa7
    case default
      call error("Invalid number of generations for ExtLagrangian")
    end select

  end subroutine ExtLagrangian_init


  !> Turns on the integrator.
  !!
  !! The integrator will start of filling up its database with the subsequent
  !! vectors passed to it. Whether the database filling is complete, can be
  !! queried by the `needsConvergedValues()` method. When it returns `.false.`
  !! the integrator is ready.
  !!
  subroutine turnOn(this, nTransientSteps)


    !> Instance variable.
    class(ExtLagrangian), intent(inout) :: this


    !> Nr. of transient steps to do *additional* to the ones needed to fill up
    !! the integrator. During those additional steps, the integrator still needs
    !! converged quantities and will use them to correct its internal database
    !! with the predicted input charges to enable a smoother change between the
    !! fully converged calculations and XL predicted ones.
    integer, intent(in), optional :: nTransientSteps

    integer :: nTransientSteps0

    if (present(nTransientSteps)) then
      nTransientSteps0 = nTransientSteps
    else
      nTransientSteps0 = 0
    end if
    this%phase = phases%fillingUp
    this%iStep = 1
    this%nSteps = this%nTimeSteps + 1
    this%nTransientSteps = nTransientSteps0

  end subroutine turnOn


  !> Reads the last output and provides the input for the next timestep.
  !!
  subroutine getNextInput(this, outLast, inNext)


    !> Instance.
    class(ExtLagrangian), intent(inout) :: this


    !> Output quantity of the last iteration.
    real(dp), intent(in) :: outLast(:)


    !> Input quantity for the next iteration
    real(dp), intent(out) :: inNext(:)

    real(dp), allocatable :: diff(:)
    integer :: ind, ind0, ind1
    integer :: ii
    real(dp) :: cc

    select case (this%phase)

    case (phases%off)

      inNext(:) = outLast

    case (phases%fillingUp)

      this%ind = modIndex(this%ind + 1, this%nTimeSteps)
      this%auxVectors(:, this%ind) = outLast
      inNext(:) = outLast

    case (phases%interpolating, phases%on)

      ind0 = this%ind

      ! Override last prediction in our database with an interpolation between
      ! it and the exact result, which just have been passed.
      if (this%phase == phases%interpolating) then
        cc = real(this%iStep, dp) / real(this%nSteps + 1, dp)
        this%auxVectors(:,ind0) = cc * this%auxVectors(:,ind0)&
            & + (1.0_dp - cc) * outLast
      end if

      ! Build new auxiliary vector by integration
      allocate(diff(size(inNext)))
      ind1 = modIndex(ind0 - 1, this%nTimeSteps)
      inNext(:) = 2.0_dp * this%auxVectors(:,ind0) - this%auxVectors(:,ind1)
      diff(:) = outLast - this%auxVectors(:,ind0)
      if (allocated(this%precondMtx)) then
        diff(:) = matmul(this%precondMtx, diff)
      end if
      inNext(:) = inNext + this%kappa * this%scale * diff

      ! Add dissipation
      do ii = 0, this%nTimeSteps
        ind = modIndex(ind0 - ii, this%nTimeSteps)
        inNext(:) = inNext&
            & + this%alpha * this%auxCoeffs(ii) * this%auxVectors(:,ind)
      end do

      ! Store predicted vector in the database
      this%ind = modIndex(this%ind + 1, this%nTimeSteps)
      this%auxVectors(:,this%ind) = inNext
    end select

    call this%updatePhaseAndSteps()

  end subroutine getNextInput

  !> helper function
  pure function modIndex(ind, nTimeSteps)
    integer, intent(in) :: ind
    integer, intent(in) :: nTimeSteps
    integer :: modIndex

    modIndex = modulo(ind - 1, nTimeSteps + 1) + 1

  end function modIndex

  !> Whether next output quantity passed to the integrator should still contain
  !! fully converged values.
  !!
  function needsConvergedValues(this) result(needsConverged)


    !> Instance.
    class(ExtLagrangian), intent(in) :: this


    !> Whether converged values are needed.
    logical :: needsConverged

    needsConverged = any(this%phase ==&
        & [phases%off, phases%fillingUp, phases%interpolating])

  end function needsConvergedValues


  !> Sets a preconditioner for the integrator.
  !!
  subroutine setPreconditioner(this, scale, precondMtx)


    !> Instance variable.
    class(ExtLagrangian), intent(inout) :: this


    !> Scaling factor for the difference vector (e.g. scaling factor for
    !! SCF-free XLBOMD). Default: 1.0.
    real(dp), intent(in), optional :: scale


    !> Preconditioning matrix for the difference vector (e.g. inverse Jacobian)
    !! Default: identity matrix.
    real(dp), intent(in), optional :: precondMtx(:,:)

  #:block DEBUG_CODE
    if (present(precondMtx)) then
      @:ASSERT(all(shape(precondMtx) == [this%nElems, this%nElems]))
    end if
  #:endblock DEBUG_CODE

    if (present(scale)) then
      this%scale = scale
    end if
    if (present(precondMtx)) then
      if (allocated(this%precondMtx)) then
        deallocate(this%precondMtx)
        allocate(this%precondMtx(this%nElems, this%nElems))
      end if
      this%precondMtx(:,:) = precondMtx
    end if

  end subroutine setPreconditioner

  ! Private methods


  !> helper function
  subroutine updatePhaseAndSteps(this)
    class(ExtLagrangian), intent(inout) :: this

    if (any(this%phase == [phases%off, phases%on])) then
      return
    end if

    this%iStep = this%iStep + 1
    if (this%phase == phases%fillingUp .and. this%iStep > this%nSteps) then
      this%phase = phases%interpolating
      this%iStep = 1
      this%nSteps = this%nTransientSteps
    end if
    if (this%phase == phases%interpolating .and. this%iStep > this%nSteps) then
      this%phase = phases%on
    end if

  end subroutine updatePhaseAndSteps

end module dftbp_extlagrangian
