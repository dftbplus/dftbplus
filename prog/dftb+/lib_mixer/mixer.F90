!--------------------------------------------------------------------------------------------------!
!  DFTB+: general package for performing fast atomistic simulations                                !
!  Copyright (C) 2017  DFTB+ developers group                                                      !
!                                                                                                  !
!  See the LICENSE file for terms of usage and distribution.                                       !
!--------------------------------------------------------------------------------------------------!

#:include 'common.fypp'

!!* Contains a general mixer which calls the desired real mixers.
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

  !!* Interface type for various mixers.
  type OMixer
    private
    integer :: mixerType
    type(OSimpleMixer),   allocatable :: pSimpleMixer
    type(OAndersonMixer), allocatable :: pAndersonMixer
    type(OBroydenMixer),  allocatable :: pBroydenMixer
    type(ODIISMixer),  allocatable :: pDIISMixer
  end type OMixer

  !!* Inits mixer
  interface init
    module procedure Mixer_initSimple
    module procedure Mixer_initAnderson
    module procedure Mixer_initBroyden
    module procedure Mixer_initDIIS
  end interface

  !!* Resets mixer
  interface reset
    module procedure Mixer_reset
  end interface reset

  !!* Does actual mixing
  interface mix
    module procedure Mixer_mix
  end interface mix

  !!* Is J^-1 available?
  interface hasInverseJacobian
    module procedure Mixer_hasInverseJacobian
  end interface hasInverseJacobian

  !!* Return J^-1 if possible
  interface getInverseJacobian
    module procedure Mixer_getInverseJacobian
  end interface getInverseJacobian


  public :: OMixer
  public :: init, reset, mix
  public :: hasInverseJacobian, getInverseJacobian


  !! Constants for the different mixers.
  integer, parameter :: iSimpleMixer = 1
  integer, parameter :: iAndersonMixer = 2
  integer, parameter :: iBroydenMixer = 3
  integer, parameter :: iDIISMixer = 4

contains

  !!* Initializes a mixer as simple mixer.
  !!* @param self Mixer instance
  !!* @param pSimple Pointer to a valid simple mixer instance.
  subroutine Mixer_initSimple(self, pSimple)
    type(OMixer), intent(out) :: self
    type(OSimpleMixer), allocatable, intent(inout) :: pSimple

    self%mixerType = iSimpleMixer
    call move_alloc(pSimple, self%pSimpleMixer)

  end subroutine Mixer_initSimple


  !!* Initializes a mixer as Anderson mixer.
  !!* @param self Mixer instance
  !!* @param pAnderson Pointer to a valid Anderson mixer instance.
  subroutine Mixer_initAnderson(self, pAnderson)
    type(OMixer), intent(out) :: self
    type(OAndersonMixer), allocatable, intent(inout) :: pAnderson

    self%mixerType = iAndersonMixer
    call move_alloc(pAnderson, self%pAndersonMixer)

  end subroutine Mixer_initAnderson


  !!* Initializes a mixer as Broyden mixer
  !!* @param self Mixer instance
  !!* @param pBroyed Pointer to a valid Broyden mixer instance.
  subroutine Mixer_initBroyden(self, pBroyden)
    type(OMixer), intent(out) :: self
    type(OBroydenMixer), allocatable, intent(inout) :: pBroyden

    self%mixerType = iBroydenMixer
    call move_alloc(pBroyden, self%pBroydenMixer)

  end subroutine Mixer_initBroyden


  !!* Initializes a mixer as DIIS mixer
  !!* @param self Mixer instance
  !!* @param pDIIS Pointer to a valid DIIS mixer instance.
  subroutine Mixer_initDIIS(self, pDIIS)
    type(OMixer), intent(out) :: self
    type(ODIISMixer), allocatable, intent(inout) :: pDIIS

    self%mixerType = iDIISMixer
    call move_alloc(pDIIS, self%pDIISMixer)

  end subroutine Mixer_initDIIS


  !!* Resets the mixer
  !!* @param self  Mixer instance.
  !!* @param nElem Size of the vectors to mix.
  subroutine Mixer_reset(self, nElem)
    type(OMixer), intent(inout) :: self
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


  !!* Mixes to vectors together
  !!* @param self Mixer instance.
  !!* @param qInpRes Input vector on entry, result vector on exit.
  !!* @param qDiff   Difference vector
  subroutine Mixer_mix(self, qInpRes, qDiff)
    type(OMixer), intent(inout) :: self
    real(dp),      intent(inout) :: qInpRes(:)
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



  !!* Tells whether the mixer is able to provide the inverse Jacobian.
  !!* @param self  Mixer instance.
  !!* @param nElem Size of the vectors to mix.
  function Mixer_hasInverseJacobian(self) result(has)
    type(OMixer), intent(inout) :: self
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


  subroutine Mixer_getInverseJacobian(self, invJac)
    type(OMixer), intent(inout) :: self
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
