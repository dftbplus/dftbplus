!--------------------------------------------------------------------------------------------------!
!  DFTB+: general package for performing fast atomistic simulations                                !
!  Copyright (C) 2006 - 2024  DFTB+ developers group                                               !
!                                                                                                  !
!  See the LICENSE file for terms of usage and distribution.                                       !
!--------------------------------------------------------------------------------------------------!

#:include 'common.fypp'
#:set FLAVOURS = [('cmplx', 'complex', 'Cmplx'), ('real', 'real', 'Real')]

!> Provides a general mixer which contains the desired actual mixer.
module dftbp_mixer_mixer
  use dftbp_common_accuracy, only : dp
  use dftbp_io_message, only : error
#:for NAME, TYPE, LABEL in FLAVOURS
  use dftbp_mixer_andersonmixer, only : TAndersonMixer${LABEL}$, TAndersonMixer${LABEL}$_mix,&
      & TAndersonMixer${LABEL}$_reset
  use dftbp_mixer_broydenmixer, only : TBroydenMixer${LABEL}$, TBroydenMixer${LABEL}$_mix,&
      & TBroydenMixer${LABEL}$_reset
  use dftbp_mixer_diismixer, only : TDiisMixer${LABEL}$, TDiisMixer${LABEL}$_mix,&
      & TDiisMixer${LABEL}$_reset
  use dftbp_mixer_simplemixer, only : TSimpleMixer${LABEL}$, TSimpleMixer${LABEL}$_mix,&
      & TSimpleMixer${LABEL}$_reset
#:endfor
  use dftbp_mixer_broydenmixer, only : TBroydenMixerReal_getInverseJacobian
  implicit none

  private
#:for NAME, TYPE, LABEL in FLAVOURS
  public :: TMixer${LABEL}$
  public :: TMixer${LABEL}$_init, TMixer${LABEL}$_reset, TMixer${LABEL}$_mix
#:endfor
  public :: TMixerReal_hasInverseJacobian, TMixerReal_getInverseJacobian
  public :: mixerTypes


#:for NAME, TYPE, LABEL in FLAVOURS
  !> Interface type for various mixers.
  type TMixer${LABEL}$
    private

    !> Numerical type of mixer 1:4
    integer :: mixerType

    !> Simple mixer instance
    type(TSimpleMixer${LABEL}$), allocatable :: pSimpleMixer

    !> Anderson mixer instance
    type(TAndersonMixer${LABEL}$), allocatable :: pAndersonMixer

    !> Broyden mixer instance
    type(TBroydenMixer${LABEL}$), allocatable :: pBroydenMixer

    !> Modified DIIS mixer instance
    type(TDiisMixer${LABEL}$), allocatable :: pDIISMixer

  end type TMixer${LABEL}$


  !> Initialises specific mixer in use
  interface TMixer${LABEL}$_init
    module procedure TMixer${LABEL}$_initSimple
    module procedure TMixer${LABEL}$_initAnderson
    module procedure TMixer${LABEL}$_initBroyden
    module procedure TMixer${LABEL}$_initDiis
  end interface TMixer${LABEL}$_init


  !> Does the actual mixing
  interface TMixer${LABEL}$_mix
    module procedure TMixer${LABEL}$_mix1D
    module procedure TMixer${LABEL}$_mix3D
    module procedure TMixer${LABEL}$_mix6D
  end interface TMixer${LABEL}$_mix
#:endfor


  type :: TMixerTypesEnum
    integer :: simple = 1
    integer :: anderson = 2
    integer :: broyden = 3
    integer :: diis = 4
  end type TMixerTypesEnum

  !> Contains mixer types
  type(TMixerTypesEnum), parameter :: mixerTypes = TMixerTypesEnum()


contains

#:for NAME, TYPE, LABEL in FLAVOURS
  !> Initializes a simple mixer.
  subroutine TMixer${LABEL}$_initSimple(this, pSimple)

    !> Mixer instance
    type(TMixer${LABEL}$), intent(out) :: this

    !> A valid Simple mixer instance on exit
    type(TSimpleMixer${LABEL}$), allocatable, intent(inout) :: pSimple

    this%mixerType = mixerTypes%simple
    call move_alloc(pSimple, this%pSimpleMixer)

  end subroutine TMixer${LABEL}$_initSimple


  !> Initializes an Anderson mixer.
  subroutine TMixer${LABEL}$_initAnderson(this, pAnderson)

    !> Mixer instance
    type(TMixer${LABEL}$), intent(out) :: this

    !> A valid Anderson mixer instance on exit
    type(TAndersonMixer${LABEL}$), allocatable, intent(inout) :: pAnderson

    this%mixerType = mixerTypes%anderson
    call move_alloc(pAnderson, this%pAndersonMixer)

  end subroutine TMixer${LABEL}$_initAnderson


  !> Initializes a Broyden mixer.
  subroutine TMixer${LABEL}$_initBroyden(this, pBroyden)

    !> Mixer instance
    type(TMixer${LABEL}$), intent(out) :: this

    !> A valid Broyden mixer instance on exit
    type(TBroydenMixer${LABEL}$), allocatable, intent(inout) :: pBroyden

    this%mixerType = mixerTypes%broyden
    call move_alloc(pBroyden, this%pBroydenMixer)

  end subroutine TMixer${LABEL}$_initBroyden


  !> Initializes a DIIS mixer.
  subroutine TMixer${LABEL}$_initDiis(this, pDIIS)

    !> Mixer instance
    type(TMixer${LABEL}$), intent(out) :: this

    !> A valid DIIS mixer instance on exit
    type(TDiisMixer${LABEL}$), allocatable, intent(inout) :: pDIIS

    this%mixerType = mixerTypes%diis
    call move_alloc(pDIIS, this%pDIISMixer)

  end subroutine TMixer${LABEL}$_initDiis


  !> Resets the mixer.
  subroutine TMixer${LABEL}$_reset(this, nElem)

    !> Mixer instance
    type(TMixer${LABEL}$), intent(inout) :: this

    !> Size of the vectors to mix
    integer, intent(in) :: nElem

    select case (this%mixerType)
    case(mixerTypes%simple)
      call TSimpleMixer${LABEL}$_reset(this%pSimpleMixer, nElem)
    case(mixerTypes%anderson)
      call TAndersonMixer${LABEL}$_reset(this%pAndersonMixer, nElem)
    case(mixerTypes%broyden)
      call TBroydenMixer${LABEL}$_reset(this%pBroydenMixer, nElem)
    case(mixerTypes%diis)
      call TDiisMixer${LABEL}$_reset(this%pDIISMixer, nElem)
    end select

  end subroutine TMixer${LABEL}$_reset


  !> Mixes two vectors.
  subroutine TMixer${LABEL}$_mix1D(this, qInpRes, qDiff)

    !> Mixer instance
    type(TMixer${LABEL}$), intent(inout) :: this

    !> Input vector on entry, result vector on exit
    ${TYPE}$(dp), intent(inout) :: qInpRes(:)

    !> Difference between input and output vectors (measure of lack of convergence)
    ${TYPE}$(dp), intent(in) :: qDiff(:)

    select case (this%mixerType)
    case(mixerTypes%simple)
      call TSimpleMixer${LABEL}$_mix(this%pSimpleMixer, qInpRes, qDiff)
    case(mixerTypes%anderson)
      call TAndersonMixer${LABEL}$_mix(this%pAndersonMixer, qInpRes, qDiff)
    case(mixerTypes%broyden)
      call TBroydenMixer${LABEL}$_mix(this%pBroydenMixer, qInpRes, qDiff)
    case(mixerTypes%diis)
      call TDiisMixer${LABEL}$_mix(this%pDIISMixer, qInpRes, qDiff)
    end select

  end subroutine TMixer${LABEL}$_mix1D


  !> Mixes two 3D matrices.
  subroutine TMixer${LABEL}$_mix3D(this, qInpResSqr, qDiffSqr)

    !> Mixer instance
    type(TMixer${LABEL}$), intent(inout) :: this

    !> Input vector on entry, result vector on exit
    ${TYPE}$(dp), intent(inout), contiguous, target :: qInpResSqr(:,:,:)

    !> Difference between input and output vectors (measure of lack of convergence)
    ${TYPE}$(dp), intent(in), contiguous, target :: qDiffSqr(:,:,:)

    !! Difference between input and output vectors (1D pointer)
    ${TYPE}$(dp), pointer :: qDiff(:)

    !! Input vector on entry, result vector on exit (1D pointer)
    ${TYPE}$(dp), pointer :: qInpRes(:)

    qDiff(1:size(qDiffSqr)) => qDiffSqr
    qInpRes(1:size(qInpResSqr)) => qInpResSqr

    call TMixer${LABEL}$_mix1D(this, qInpRes, qDiff)

  end subroutine TMixer${LABEL}$_mix3D


  !> Mixes two 6D matrices.
  subroutine TMixer${LABEL}$_mix6D(this, qInpResSqr, qDiffSqr)

    !> Mixer instance
    type(TMixer${LABEL}$), intent(inout) :: this

    !> Input vector on entry, result vector on exit
    ${TYPE}$(dp), intent(inout), contiguous, target :: qInpResSqr(:,:,:,:,:,:)

    !> Difference between input and output vectors (measure of lack of convergence)
    ${TYPE}$(dp), intent(in), contiguous, target :: qDiffSqr(:,:,:,:,:,:)

    !!Difference between input and output vectors (1D pointer)
    ${TYPE}$(dp), pointer :: qDiff(:)

    !! Input vector on entry, result vector on exit (1D pointer)
    ${TYPE}$(dp), pointer :: qInpRes(:)

    qDiff(1:size(qDiffSqr)) => qDiffSqr
    qInpRes(1:size(qInpResSqr)) => qInpResSqr

    call TMixer${LABEL}$_mix1D(this, qInpRes, qDiff)

  end subroutine TMixer${LABEL}$_mix6D
#:endfor


  !> Tells whether the mixer is able to provide the inverse Jacobian.
  function TMixerReal_hasInverseJacobian(this) result(has)

    !> Mixer instance
    type(TMixerReal), intent(inout) :: this

    !> Size of the vectors to mix
    logical :: has

    select case (this%mixerType)
    case(mixerTypes%simple)
      has = .false.
    case (mixerTypes%anderson)
      has = .false.
    case (mixerTypes%broyden)
      has = .true.
    case (mixerTypes%diis)
      has = .false.
    end select

  end function TMixerReal_hasInverseJacobian


  !> Returns an inverse Jacobian if possible, halting if not.
  subroutine TMixerReal_getInverseJacobian(this, invJac)

    !> Mixer instance
    type(TMixerReal), intent(inout) :: this

    !> Inverse Jacobian matrix if available
    real(dp), intent(out) :: invJac(:,:)

    select case (this%mixerType)
    case(mixerTypes%simple)
      call error("Simple mixer does not provide inverse Jacobian")
    case (mixerTypes%anderson)
      call error("Anderson mixer does not provide inverse Jacobian")
    case (mixerTypes%broyden)
      call TBroydenMixerReal_getInverseJacobian(this%pBroydenMixer, invJac)
    case (mixerTypes%diis)
      call error("DIIS mixer does not provide inverse Jacobian")
    end select

  end subroutine TMixerReal_getInverseJacobian

end module dftbp_mixer_mixer
