!--------------------------------------------------------------------------------------------------!
!  DFTB+: general package for performing fast atomistic simulations                                !
!  Copyright (C) 2018  DFTB+ developers group                                                      !
!                                                                                                  !
!  See the LICENSE file for terms of usage and distribution.                                       !
!--------------------------------------------------------------------------------------------------!

#:include 'common.fypp'

!> Provides a general mixer which contains the desired actual mixer.
module mixer
  use assert
  use accuracy
  use simplemixer
  use andersonmixer
  use broydenmixer
  use diismixer
  use message
  implicit none

  private


  !> Interface type for various mixers.
  type OMixer
    private

    !> numerical type of mixer 1:4
    integer :: mixerType

    !> simple mixer instance
    type(OSimpleMixer),   allocatable :: pSimpleMixer

    !> Anderson mixer instance
    type(OAndersonMixer), allocatable :: pAndersonMixer

    !> Broyden mixer instance
    type(OBroydenMixer),  allocatable :: pBroydenMixer

    !> modified DIIS mixer instance
    type(ODIISMixer),  allocatable :: pDIISMixer
  end type OMixer


  !> Initialises specific mixer in use
  interface init
    module procedure Mixer_initSimple
    module procedure Mixer_initAnderson
    module procedure Mixer_initBroyden
    module procedure Mixer_initDIIS
  end interface


  !> Resets mixer
  interface reset
    module procedure Mixer_reset
  end interface reset


  !> Does the actual mixing
  interface mix
    module procedure Mixer_mix
  end interface mix


  !> Is J^-1 available?
  interface hasInverseJacobian
    module procedure Mixer_hasInverseJacobian
  end interface hasInverseJacobian


  !> Return J^-1 if possible
  interface getInverseJacobian
    module procedure Mixer_getInverseJacobian
  end interface getInverseJacobian

  public :: OMixer
  public :: init, reset, mix
  public :: hasInverseJacobian, getInverseJacobian


  !> Identifying constant for each of the different mixers.
  integer, parameter :: iSimpleMixer = 1
  integer, parameter :: iAndersonMixer = 2
  integer, parameter :: iBroydenMixer = 3
  integer, parameter :: iDIISMixer = 4

contains


  !> Initializes a simple mixer.
  subroutine Mixer_initSimple(self, pSimple)

    !> Mixer instance
    type(OMixer), intent(out) :: self

    !> A valid simple mixer instance on exit.
    type(OSimpleMixer), allocatable, intent(inout) :: pSimple

    self%mixerType = iSimpleMixer
    call move_alloc(pSimple, self%pSimpleMixer)

  end subroutine Mixer_initSimple


  !> Initializes an Anderson mixer.
  subroutine Mixer_initAnderson(self, pAnderson)

    !> Mixer instance
    type(OMixer), intent(out) :: self

    !> A valid Anderson mixer instance on exit.
    type(OAndersonMixer), allocatable, intent(inout) :: pAnderson

    self%mixerType = iAndersonMixer
    call move_alloc(pAnderson, self%pAndersonMixer)

  end subroutine Mixer_initAnderson


  !> Initializes a Broyden mixer
  subroutine Mixer_initBroyden(self, pBroyden)

    !> Mixer instance
    type(OMixer), intent(out) :: self

    !> A valid Broyden mixer instance on exit.
    type(OBroydenMixer), allocatable, intent(inout) :: pBroyden

    self%mixerType = iBroydenMixer
    call move_alloc(pBroyden, self%pBroydenMixer)

  end subroutine Mixer_initBroyden


  !> Initializes a DIIS mixer
  subroutine Mixer_initDIIS(self, pDIIS)

    !> Mixer instance
    type(OMixer), intent(out) :: self

    !> A valid DIIS mixer instance on exit.
    type(ODIISMixer), allocatable, intent(inout) :: pDIIS

    self%mixerType = iDIISMixer
    call move_alloc(pDIIS, self%pDIISMixer)

  end subroutine Mixer_initDIIS


  !> Resets the mixer
  subroutine Mixer_reset(self, nElem)

    !> Mixer instance.
    type(OMixer), intent(inout) :: self

    !> Size of the vectors to mix.
    integer, intent(in) :: nElem

    select case (self%mixerType)
    case(iSimpleMixer)
      call reset(self%pSimpleMixer, nElem)
    case (iAndersonMixer)
      call reset(self%pAndersonMixer, nElem)
    case (iBroydenMixer)
      call reset(self%pBroydenMixer, nElem)
    case (iDIISMixer)
      call reset(self%pDIISMixer, nElem)
    end select

  end subroutine Mixer_reset


  !> Mixes vectors together
  subroutine Mixer_mix(self, qInpRes, qDiff)

    !> Mixer instance.
    type(OMixer), intent(inout) :: self

    !> Input vector on entry, result vector on exit.
    real(dp),      intent(inout) :: qInpRes(:)

    !> Difference between input and output vectors (measure of lack of convergence)
    real(dp),      intent(in) :: qDiff(:)

    select case (self%mixerType)
    case (iSimpleMixer)
      call mix(self%pSimpleMixer, qInpRes, qDiff)
    case (iAndersonMixer)
      call mix(self%pAndersonMixer, qInpRes, qDiff)
    case (iBroydenMixer)
      call mix(self%pBroydenMixer, qInpRes, qDiff)
    case (iDIISMixer)
      call mix(self%pDIISMixer, qInpRes, qDiff)
    end select

  end subroutine Mixer_mix


  !> Tells whether the mixer is able to provide the inverse Jacobian.
  function Mixer_hasInverseJacobian(self) result(has)

    !> Mixer instance.
    type(OMixer), intent(inout) :: self

    !> Size of the vectors to mix.
    logical :: has

    select case (self%mixerType)
    case(iSimpleMixer)
      has = .false.
    case (iAndersonMixer)
      has = .false.
    case (iBroydenMixer)
      has = .true.
    case (iDIISMixer)
      has = .false.
    end select

  end function Mixer_hasInverseJacobian


  !> Return an inverse Jacobian if possible, halting if not
  subroutine Mixer_getInverseJacobian(self, invJac)

    !> Mixer instance.
    type(OMixer), intent(inout) :: self

    !> Inverse Jacobian matrix if available
    real(dp), intent(out) :: invJac(:,:)

    select case (self%mixerType)
    case(iSimpleMixer)
      call error("Simple mixer does not provide inverse Jacobian")
    case (iAndersonMixer)
      call error("Anderson mixer does not provide inverse Jacobian")
    case (iBroydenMixer)
      call getInverseJacobian(self%pBroydenMixer, invJac)
    case (iDIISMixer)
      call error("DIIS mixer does not provide inverse Jacobian")
    end select

  end subroutine Mixer_getInverseJacobian

end module mixer
