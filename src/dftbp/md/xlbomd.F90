!--------------------------------------------------------------------------------------------------!
!  DFTB+: general package for performing fast atomistic simulations                                !
!  Copyright (C) 2006 - 2025  DFTB+ developers group                                               !
!                                                                                                  !
!  See the LICENSE file for terms of usage and distribution.                                       !
!--------------------------------------------------------------------------------------------------!

#:include 'common.fypp'

!> Implements the Extended Lagrangian Born-Oppenheimer MD.
!> Aradi et al. Extended lagrangian density functional tight-binding molecular dynamics for
!> molecules and solids. J. Chem. Theory Comput. 11:3357-3363, 2015
module dftbp_md_xlbomd
  use dftbp_common_accuracy, only : dp
  use dftbp_md_extlagrangian, only : ExtLagrangian, ExtLagrangian_init, ExtLagrangianInp
  implicit none

  private
  public :: TXLBOMDInp, TXLBOMD, Xlbomd_init


  !> Input for the Xlbomd driver.
  type :: TXLBOMDInp

    !> Number of generation to consider during the integration (5, 6, 7)
    integer :: nKappa

    !> Scaling factor (only for fast Xlbomd, otherwise 1.0)
    real(dp) :: scale

    !> SCC parameter override: minimal SCC iterations per timestep
    integer :: minSCCIter

    !> SCC parameter override: maximal SCC iterations per timestep
    integer :: maxSCCIter

    !> SCC parameter override: scc tolerance
    real(dp) :: sccTol

    !> Number of time steps to precede the actual start of the XL integrator
    integer :: nPreSteps

    !> Number of full SCC steps after the XL integrator has been started. Those
    !> steps are used to fill up the integrator.
    integer :: nFullSCCSteps

    !> Number of transient steps (during which prediction is corrected)
    integer :: nTransientSteps

  end type TXLBOMDInp


  !> Contains the data for the Xlbomd driver.
  type :: TXLBOMD
    private
    type(ExtLagrangian) :: extLagr
    integer :: nKappa
    integer :: minSCCIter, maxSCCIter, minSCCIter0, maxSCCIter0
    real(dp) :: sccTol, sccTol0
    integer :: iStep
    integer :: iStartXl, nPreSteps, nTransientSteps, nFullSCCSteps
  contains
    procedure :: setDefaultSCCParameters
    procedure :: isActive
    procedure :: getSCCParameters
    procedure :: getNextCharges
  end type TXLBOMD

contains


  !> Initializes the Xlbomd instance.
  subroutine Xlbomd_init(this, input, nElems)

    !> Instance.
    type(TXLBOMD), intent(out) :: this

    !> Basic input parameters.
    type(TXLBOMDInp), intent(in) :: input

    !> Nr. of elements in the charge vector
    integer, intent(in) :: nElems

    type(ExtLagrangianInp) :: extLagrInp

    @:ASSERT(input%scale >= 0.0_dp)
    @:ASSERT(input%minSCCIter >= 1)
    @:ASSERT(input%maxSCCIter >= 1 .and. input%maxSCCIter >= input%minSCCIter)
    @:ASSERT(input%sccTol > 0.0_dp)

    extLagrInp%nTimeSteps = input%nKappa
    extLagrInp%nElems = nElems
    call ExtLagrangian_init(this%extLagr, extLagrInp)

    call this%extLagr%setPreconditioner(scale=input%scale)

    this%nKappa = input%nKappa
    this%iStep = 1
    this%nTransientSteps = input%nTransientSteps
    this%nPreSteps = input%nPreSteps
    this%nFullSCCSteps = input%nFullSCCSteps
    this%iStartXl = this%nPreSteps + input%nFullSCCSteps - input%nKappa

    this%minSCCIter = input%minSCCIter
    this%maxSCCIter = input%maxSCCIter
    this%sccTol = input%sccTol

  end subroutine Xlbomd_init


  !> Delivers charges for the next time step.
  subroutine getNextCharges(this, qCurrent, qNext)

    !> Instance.
    class(TXLBOMD), intent(inout) :: this

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


  !> Sets default SCC parameters to return, if no override is done.
  subroutine setDefaultSCCParameters(this, minSCCIter, maxSCCIter, sccTol)

    !> Instance.
    class(TXLBOMD), intent(inout) :: this

    !> Minimal number of SCC iterations.
    integer, intent(in) :: minSCCIter

    !> Maximal number of SCC iterations.
    integer, intent(in) :: maxSCCIter

    !> SCC tolerance.
    real(dp), intent(in) :: sccTol

    this%minSCCIter0 = minSCCIter
    this%maxSCCIter0 = maxSCCIter
    this%sccTol0 = sccTol

  end subroutine setDefaultSCCParameters


  !> Signals whether XLBOMD integration is active.
  function isActive(this)

    !> Instance.
    class(TXLBOMD), intent(in) :: this

    !> True if XLBOMD integration is active (no SCC convergence needed)
    logical :: isActive

    isActive = this%extLagr%needsConvergedValues()

  end function isActive


  !> Returns the SCC parameters to be used when the integrator is active.
  subroutine getSCCParameters(this, minSCCIter, maxSCCIter, sccTol)

    !> Instance.
    class(TXLBOMD), intent(in) :: this

    !> Minimal number of SCC cycles.
    integer, intent(out) :: minSCCIter

    !> Maximal number of SCC cycles
    integer, intent(out) :: maxSCCIter

    !> Tolerance for SCC convergence.
    real(dp), intent(out) :: sccTol

    if (this%extLagr%needsConvergedValues()) then
      minSCCIter = this%minSCCIter0
      maxSCCIter = this%maxSCCIter0
      sccTol = this%sccTol0
    else
      minSCCIter = this%minSCCIter
      maxSCCIter = this%maxSCCIter
      sccTol = this%sccTol
    end if

  end subroutine getSCCParameters

end module dftbp_md_xlbomd
