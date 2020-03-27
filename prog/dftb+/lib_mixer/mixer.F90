!--------------------------------------------------------------------------------------------------!
!  DFTB+: general package for performing fast atomistic simulations                                !
!  Copyright (C) 2006 - 2020  DFTB+ developers group                                               !
!                                                                                                  !
!  See the LICENSE file for terms of usage and distribution.                                       !
!--------------------------------------------------------------------------------------------------!

#:include 'common.fypp'

!> Provides a general mixer which contains the desired actual mixer.
module dftbp_mixer
  use dftbp_assert
  use dftbp_accuracy
  use dftbp_simplemixer
  use dftbp_andersonmixer
  use dftbp_broydenmixer
  use dftbp_diismixer
  use dftbp_message
  implicit none
  private

  public :: TMixer
  public :: init, reset, mix
  public :: hasInverseJacobian, getInverseJacobian
  public :: mixerTypes


  !> Interface type for various mixers.
  type TMixer
    private

    !> numerical type of mixer 1:4
    integer :: mixerType

    !> simple mixer instance
    type(TSimpleMixer),   allocatable :: pSimpleMixer

    !> Anderson mixer instance
    type(TAndersonMixer), allocatable :: pAndersonMixer

    !> Broyden mixer instance
    type(TBroydenMixer),  allocatable :: pBroydenMixer

    !> modified DIIS mixer instance
    type(TDIISMixer),  allocatable :: pDIISMixer
  end type TMixer


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


  type :: TMixerTypesEnum
    integer :: simple = 1
    integer :: anderson = 2
    integer :: broyden = 3
    integer :: diis = 4
  end type TMixerTypesEnum

  !> Contains mixer types
  type(TMixerTypesEnum), parameter :: mixerTypes = TMixerTypesEnum()

contains


  !> Initializes a simple mixer.
  subroutine Mixer_initSimple(self, pSimple)

    !> Mixer instance
    type(TMixer), intent(out) :: self

    !> A valid simple mixer instance on exit.
    type(TSimpleMixer), allocatable, intent(inout) :: pSimple

    self%mixerType = mixerTypes%simple
    call move_alloc(pSimple, self%pSimpleMixer)

  end subroutine Mixer_initSimple


  !> Initializes an Anderson mixer.
  subroutine Mixer_initAnderson(self, pAnderson)

    !> Mixer instance
    type(TMixer), intent(out) :: self

    !> A valid Anderson mixer instance on exit.
    type(TAndersonMixer), allocatable, intent(inout) :: pAnderson

    self%mixerType = mixerTypes%anderson
    call move_alloc(pAnderson, self%pAndersonMixer)

  end subroutine Mixer_initAnderson


  !> Initializes a Broyden mixer
  subroutine Mixer_initBroyden(self, pBroyden)

    !> Mixer instance
    type(TMixer), intent(out) :: self

    !> A valid Broyden mixer instance on exit.
    type(TBroydenMixer), allocatable, intent(inout) :: pBroyden

    self%mixerType = mixerTypes%broyden
    call move_alloc(pBroyden, self%pBroydenMixer)

  end subroutine Mixer_initBroyden


  !> Initializes a DIIS mixer
  subroutine Mixer_initDIIS(self, pDIIS)

    !> Mixer instance
    type(TMixer), intent(out) :: self

    !> A valid DIIS mixer instance on exit.
    type(TDIISMixer), allocatable, intent(inout) :: pDIIS

    self%mixerType = mixerTypes%diis
    call move_alloc(pDIIS, self%pDIISMixer)

  end subroutine Mixer_initDIIS


  !> Resets the mixer
  subroutine Mixer_reset(self, nElem)

    !> Mixer instance.
    type(TMixer), intent(inout) :: self

    !> Size of the vectors to mix.
    integer, intent(in) :: nElem

    select case (self%mixerType)
    case(mixerTypes%simple)
      call reset(self%pSimpleMixer, nElem)
    case (mixerTypes%anderson)
      call reset(self%pAndersonMixer, nElem)
    case (mixerTypes%broyden)
      call reset(self%pBroydenMixer, nElem)
    case (mixerTypes%diis)
      call reset(self%pDIISMixer, nElem)
    end select

  end subroutine Mixer_reset


  !> Mixes vectors together
  subroutine Mixer_mix(self, qInpRes, qDiff)

    !> Mixer instance.
    type(TMixer), intent(inout) :: self

    !> Input vector on entry, result vector on exit.
    real(dp),      intent(inout) :: qInpRes(:)

    !> Difference between input and output vectors (measure of lack of convergence)
    real(dp),      intent(in) :: qDiff(:)

    select case (self%mixerType)
    case (mixerTypes%simple)
      call mix(self%pSimpleMixer, qInpRes, qDiff)
    case (mixerTypes%anderson)
      call mix(self%pAndersonMixer, qInpRes, qDiff)
    case (mixerTypes%broyden)
      call mix(self%pBroydenMixer, qInpRes, qDiff)
    case (mixerTypes%diis)
      call mix(self%pDIISMixer, qInpRes, qDiff)
    end select

  end subroutine Mixer_mix


  !> Tells whether the mixer is able to provide the inverse Jacobian.
  function Mixer_hasInverseJacobian(self) result(has)

    !> Mixer instance.
    type(TMixer), intent(inout) :: self

    !> Size of the vectors to mix.
    logical :: has

    select case (self%mixerType)
    case(mixerTypes%simple)
      has = .false.
    case (mixerTypes%anderson)
      has = .false.
    case (mixerTypes%broyden)
      has = .true.
    case (mixerTypes%diis)
      has = .false.
    end select

  end function Mixer_hasInverseJacobian


  !> Return an inverse Jacobian if possible, halting if not
  subroutine Mixer_getInverseJacobian(self, invJac)

    !> Mixer instance.
    type(TMixer), intent(inout) :: self

    !> Inverse Jacobian matrix if available
    real(dp), intent(out) :: invJac(:,:)

    select case (self%mixerType)
    case(mixerTypes%simple)
      call error("Simple mixer does not provide inverse Jacobian")
    case (mixerTypes%anderson)
      call error("Anderson mixer does not provide inverse Jacobian")
    case (mixerTypes%broyden)
      call getInverseJacobian(self%pBroydenMixer, invJac)
    case (mixerTypes%diis)
      call error("DIIS mixer does not provide inverse Jacobian")
    end select

  end subroutine Mixer_getInverseJacobian

end module dftbp_mixer
