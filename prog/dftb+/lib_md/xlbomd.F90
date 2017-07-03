!--------------------------------------------------------------------------------------------------!
!  DFTB+: general package for performing fast atomistic simulations                                !
!  Copyright (C) 2017  DFTB+ developers group                                                      !
!                                                                                                  !
!  See the LICENSE file for terms of usage and distribution.                                       !
!--------------------------------------------------------------------------------------------------!

#:include 'common.fypp'

!> Implements the Extended Lagrangian Born-Oppenheimer MD.
!!
!! \see Aradi et al. Extended lagrangian density functional
!!  tight-binding molecular dynamics for molecules and
!!  solids. J. Chem. Theory Comput. 11:3357-3363, 2015
!!
module xlbomd_module
  use assert
  use accuracy
  use message
  use extlagrangian_module
  implicit none
  private

  public :: XlbomdInp, Xlbomd, Xlbomd_init


  character(*), parameter :: JacobianKernelFile = "neginvjac.dat"


  !> Input for the Xlbomd driver.
  !!
  type :: XlbomdInp

    !> Number of generation to consider during the integration (5, 6, 7)
    integer :: nKappa

    !> Scaling factor (only for fast Xlbomd, otherwise 1.0)
    real(dp) :: scale

    !> Scc parameter override: minimal Scc iterations per timestep
    integer :: minSccIter

    !> Scc parameter override: maximal Scc iterations per timestep
    integer :: maxSccIter

    !> Scc parameter override: scc tolerance
    real(dp) :: sccTol

    !> Number of time steps to precede the actual start of the XL integrator
    integer :: nPreSteps

    !> Number of full SCC steps after the XL integrator has been started. Those
    !! steps are used to fill up the integrator and to average the Jacobian.
    integer :: nFullSccSteps

    !> Number of transient steps (during which prediction is corrected)
    integer :: nTransientSteps

    !> Whether Xlbomd should use inverse Jacobian instead of scaling
    logical :: useInverseJacobian

    !> Whether Jacobian should be read from disk.
    logical :: readInverseJacobian

  end type XlbomdInp


  !> Contains the data for the Xlbomd driver.
  !!
  type :: Xlbomd
    private
    type(ExtLagrangian) :: extLagr
    integer :: nKappa
    integer :: minSccIter, maxSccIter, minSccIter0, maxSccIter0
    real(dp) :: sccTol, sccTol0
    integer :: iStep
    integer :: iStartXl, nPreSteps, nTransientSteps, nFullSccSteps
    logical :: useInverseJacobian, readInverseJacobian
    real(dp), allocatable :: invJacobian(:,:)
  contains
    procedure :: setDefaultSccParameters
    procedure :: isActive
    procedure :: getSccParameters
    procedure :: getNextCharges
    procedure :: needsInverseJacobian
    procedure :: setInverseJacobian
    procedure, private :: readJacobianKernel
  end type Xlbomd


contains

  !> Initializes the Xlbomd instance.
  !!
  subroutine Xlbomd_init(this, input, nElems)

    !> Instance.
    type(Xlbomd), intent(out) :: this

    !> Basic input parameters.
    type(XlbomdInp), intent(in) :: input

    !> Nr. of elements in the charge vector
    integer, intent(in) :: nElems


    type(ExtLagrangianInp) :: extLagrInp

    @:ASSERT(input%scale >= 0.0_dp)
    @:ASSERT(input%minSccIter >= 1)
    @:ASSERT(input%maxSccIter >= 1 .and. input%maxSccIter >= input%minSccIter)
    @:ASSERT(input%sccTol > 0.0_dp)

    extLagrInp%nTimeSteps = input%nKappa
    extLagrInp%nElems = nElems
    call ExtLagrangian_init(this%extLagr, extLagrInp)

    call this%extLagr%setPreconditioner(scale=input%scale)

    this%nKappa = input%nKappa
    this%iStep = 1
    this%nTransientSteps = input%nTransientSteps
    this%nPreSteps = input%nPreSteps
    this%nFullSccSteps = input%nFullSccSteps
    this%iStartXl = this%nPreSteps + input%nFullSccSteps - input%nKappa

    this%minSccIter = input%minSccIter
    this%maxSccIter = input%maxSccIter
    this%sccTol = input%sccTol

    this%useInverseJacobian = input%useInverseJacobian
    this%readInverseJacobian = input%readInverseJacobian
    if (this%useInverseJacobian) then
      allocate(this%invJacobian(nElems, nElems))
      if (this%readInverseJacobian) then
        call this%readJacobianKernel()
      else
        this%invJacobian(:,:) = 0.0_dp
      end if
    end if

  end subroutine Xlbomd_init


  !> Delivers charges for the next time step.
  !!
  subroutine getNextCharges(this, qCurrent, qNext)

    !> Instance.
    class(Xlbomd), intent(inout) :: this

    !> Charges for current time step.
    real(dp), intent(in) :: qCurrent(:)

    !> Charges for next time step.
    real(dp), intent(out) :: qNext(:)

    if (this%iStep == this%iStartXl) then
      call this%extLagr%turnOn(this%nTransientSteps)
    end if
    call this%extLagr%getNextInput(qCurrent, qNext)
    this%iStep = this%iStep + 1

  end subroutine getNextCharges


  !> Sets default Scc parameters to return, if no override is done.
  !!
  subroutine setDefaultSccParameters(this, minSccIter, maxSccIter, sccTol)
    !> Instance.
    class(Xlbomd), intent(inout) :: this

    !> Minimal number of Scc iterations.
    integer, intent(in) :: minSccIter

    !> Maximal number of Scc iterations.
    integer, intent(in) :: maxSccIter

    !> Scc tolerance.
    real(dp), intent(in) :: sccTol


    this%minSccIter0 = minSccIter
    this%maxSccIter0 = maxSccIter
    this%sccTol0 = sccTol

  end subroutine setDefaultSccParameters


  !> Signalizes whether XLBOMD integration is active.
  !!
  function isActive(this)

    !> Instance.
    class(Xlbomd), intent(in) :: this

    !> True if XLBOMD integration is active (no SCC convergence needed)
    logical :: isActive

    isActive = this%extLagr%needsConvergedValues()

  end function isActive


  !> Returns the Scc parameters to be used when the integrator is active.
  !!
  subroutine getSccParameters(this, minSccIter, maxSccIter, sccTol)

    !> Instance.
    class(Xlbomd), intent(in) :: this

    !> Minimal number of Scc cycles.
    integer, intent(out) :: minSccIter

    !> Maximal number of Scc cycles
    integer, intent(out) :: maxSccIter

    !> Tolerance for Scc convergence.
    real(dp), intent(out) :: sccTol

    if (this%extLagr%needsConvergedValues()) then
      minSccIter = this%minSccIter0
      maxSccIter = this%maxSccIter0
      sccTol = this%sccTol0
    else
      minSccIter = this%minSccIter
      maxSccIter = this%maxSccIter
      sccTol = this%sccTol
    end if

  end subroutine getSccParameters


  !> Whether integrator needs the inverse Jacobian at the current step.
  !!
  function needsInverseJacobian(this)

    !> Instance.
    class(Xlbomd), intent(in) :: this

    !> True, if a inverse Jacobian is needed. It should be passed via the
    !! setInverseJacobian() procedure.
    logical :: needsInverseJacobian

    integer :: iStart, iEnd

    needsInverseJacobian = .false.
    if (this%useInverseJacobian .and. .not. this%readInverseJacobian) then
      iStart = this%nPreSteps + 1
      iEnd = this%nPreSteps + this%nFullSccSteps
      needsInverseJacobian = (this%iStep >= iStart .and. this%iStep <= iEnd)
    end if

  end function needsInverseJacobian


  !> Sets the inverse Jacobian for driving the Xlbomd simulation.
  !!
  !! \param this  Instance.
  !! \param invJacobian  Inverse Jacobian.
  !!
  subroutine setInverseJacobian(this, invJacobian)
    class(Xlbomd), intent(inout) :: this
    real(dp), intent(in) :: invJacobian(:,:)

    real(dp) :: coeffOld, coeffNew, normFactor
    integer :: nn

    if (this%needsInverseJacobian()) then
      nn = this%iStep - this%nPreSteps
      normFactor = 1.0_dp / sum(invJacobian(:,1))
      coeffNew = 1.0_dp / real(nn, dp)
      coeffOld = real(nn - 1, dp) * coeffNew
      this%invJacobian(:,:) = this%invJacobian * coeffOld &
          & + normFactor * invJacobian * coeffNew
      call this%extLagr%setPreconditioner(precondMtx=this%invJacobian)
    end if

  end subroutine setInverseJacobian


  subroutine readJacobianKernel(this)
    class(Xlbomd), intent(inout) :: this

    open(12, file=JacobianKernelFile, status="old", action="read")
    read(12, *) this%invJacobian
    close(12)
    this%invJacobian = transpose(this%invJacobian)
    write(*, "(A,A,A)") "Negative inverse Jacobian read from '", &
        & JacobianKernelFile, "'"

  end subroutine readJacobianKernel


end module xlbomd_module
