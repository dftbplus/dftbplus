!--------------------------------------------------------------------------------------------------!
!  DFTB+: general package for performing fast atomistic simulations                                !
!  Copyright (C) 2006 - 2023  DFTB+ developers group                                               !
!                                                                                                  !
!  See the LICENSE file for terms of usage and distribution.                                       !
!--------------------------------------------------------------------------------------------------!

#:include 'common.fypp'

!> Provides a general mixer which contains the desired actual mixer.
module dftbp_mixer_mixer
  use dftbp_common_accuracy, only : dp
  use dftbp_io_message, only : error
  use dftbp_mixer_andersonmixer, only : TAndersonMixer, TAndersonMixer_mix, TAndersonMixer_reset
  use dftbp_mixer_broydenmixer, only : TBroydenMixerReal, TBroydenMixerCmplx,&
      & TBroydenMixerReal_mix, TBroydenMixerCmplx_mix, TBroydenMixerReal_reset,&
      & TBroydenMixerCmplx_reset, TBroydenMixerReal_getInverseJacobian
  use dftbp_mixer_diismixer, only : TDiisMixer, TDiisMixer_mix, TDiisMixer_reset
  use dftbp_mixer_simplemixer, only : TSimpleMixerReal, TSimpleMixerCmplx, TSimpleMixerReal_mix,&
      & TSimpleMixerCmplx_mix, TSimpleMixerReal_reset, TSimpleMixerCmplx_reset
  implicit none

#:set FLAVOURS = [('cmplx', 'complex', 'Cmplx'), ('real', 'real', 'Real')]

  private
#:for NAME, TYPE, LABEL in FLAVOURS
  public :: TMixer${LABEL}$
  public :: TMixer${LABEL}$_init, TMixer${LABEL}$_reset, TMixer${LABEL}$_mix
#:endfor
  public :: TMixerReal_hasInverseJacobian, TMixerReal_getInverseJacobian
  public :: mixerTypes


  !> Interface type for various mixers.
  type TMixerReal
    private

    !> Numerical type of mixer 1:4
    integer :: mixerType

    !> Simple mixer instance
    type(TSimpleMixerReal), allocatable :: pSimpleMixer

    !> Anderson mixer instance
    type(TAndersonMixer), allocatable :: pAndersonMixer

    !> Broyden mixer instance
    type(TBroydenMixerReal), allocatable :: pBroydenMixer

    !> Modified DIIS mixer instance
    type(TDiisMixer), allocatable :: pDIISMixer

  end type TMixerReal


  !> Interface type for various mixers.
  type TMixerCmplx
    private

    !> Numerical type of mixer 1,3
    integer :: mixerType

    !> Simple mixer instance
    type(TSimpleMixerCmplx), allocatable :: pSimpleMixer

    !> Broyden mixer instance
    type(TBroydenMixerCmplx), allocatable :: pBroydenMixer

  end type TMixerCmplx


  !> Initialises specific mixer in use
  interface TMixerReal_init
    module procedure TMixerReal_initSimple
    module procedure TMixerReal_initAnderson
    module procedure TMixerReal_initBroyden
    module procedure TMixerReal_initDIIS
  end interface


  !> Initialises specific mixer in use
  interface TMixerCmplx_init
    module procedure TMixerCmplx_initSimple
    module procedure TMixerCmplx_initBroyden
  end interface TMixerCmplx_init


#:for NAME, TYPE, LABEL in FLAVOURS
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

    !> A valid simple mixer instance on exit.
    type(TSimpleMixer${LABEL}$), allocatable, intent(inout) :: pSimple

    this%mixerType = mixerTypes%simple
    call move_alloc(pSimple, this%pSimpleMixer)

  end subroutine TMixer${LABEL}$_initSimple


  !> Initializes a Broyden mixer
  subroutine TMixer${LABEL}$_initBroyden(this, pBroyden)

    !> Mixer instance
    type(TMixer${LABEL}$), intent(out) :: this

    !> A valid Broyden mixer instance on exit.
    type(TBroydenMixer${LABEL}$), allocatable, intent(inout) :: pBroyden

    this%mixerType = mixerTypes%broyden
    call move_alloc(pBroyden, this%pBroydenMixer)

  end subroutine TMixer${LABEL}$_initBroyden
#:endfor


  !> Initializes an Anderson mixer.
  subroutine TMixerReal_initAnderson(this, pAnderson)

    !> Mixer instance
    type(TMixerReal), intent(out) :: this

    !> A valid Anderson mixer instance on exit.
    type(TAndersonMixer), allocatable, intent(inout) :: pAnderson

    this%mixerType = mixerTypes%anderson
    call move_alloc(pAnderson, this%pAndersonMixer)

  end subroutine TMixerReal_initAnderson


  !> Initializes a DIIS mixer
  subroutine TMixerReal_initDIIS(this, pDIIS)

    !> Mixer instance
    type(TMixerReal), intent(out) :: this

    !> A valid DIIS mixer instance on exit.
    type(TDiisMixer), allocatable, intent(inout) :: pDIIS

    this%mixerType = mixerTypes%diis
    call move_alloc(pDIIS, this%pDIISMixer)

  end subroutine TMixerReal_initDIIS


  !> Resets the mixer
  subroutine TMixerReal_reset(this, nElem)

    !> Mixer instance.
    type(TMixerReal), intent(inout) :: this

    !> Size of the vectors to mix.
    integer, intent(in) :: nElem

    select case (this%mixerType)
    case(mixerTypes%simple)
      call TSimpleMixerReal_reset(this%pSimpleMixer, nElem)
    case(mixerTypes%anderson)
      call TAndersonMixer_reset(this%pAndersonMixer, nElem)
    case(mixerTypes%broyden)
      call TBroydenMixerReal_reset(this%pBroydenMixer, nElem)
    case(mixerTypes%diis)
      call TDiisMixer_reset(this%pDIISMixer, nElem)
    end select

  end subroutine TMixerReal_reset


  !> Resets the mixer
  subroutine TMixerCmplx_reset(this, nElem)

    !> Mixer instance.
    type(TMixerCmplx), intent(inout) :: this

    !> Size of the vectors to mix.
    integer, intent(in) :: nElem

    select case (this%mixerType)
    case(mixerTypes%simple)
      call TSimpleMixerCmplx_reset(this%pSimpleMixer, nElem)
    case(mixerTypes%broyden)
      call TBroydenMixerCmplx_reset(this%pBroydenMixer, nElem)
    end select

  end subroutine TMixerCmplx_reset


  !> Mixes two vectors.
  subroutine TMixerReal_mix1D(this, qInpRes, qDiff)

    !> Mixer instance.
    type(TMixerReal), intent(inout) :: this

    !> Input vector on entry, result vector on exit.
    real(dp), intent(inout) :: qInpRes(:)

    !> Difference between input and output vectors (measure of lack of convergence)
    real(dp), intent(in) :: qDiff(:)

    select case (this%mixerType)
    case(mixerTypes%simple)
      call TSimpleMixerReal_mix(this%pSimpleMixer, qInpRes, qDiff)
    case(mixerTypes%anderson)
      call TAndersonMixer_mix(this%pAndersonMixer, qInpRes, qDiff)
    case(mixerTypes%broyden)
      call TBroydenMixerReal_mix(this%pBroydenMixer, qInpRes, qDiff)
    case(mixerTypes%diis)
      call TDiisMixer_mix(this%pDIISMixer, qInpRes, qDiff)
    end select

  end subroutine TMixerReal_mix1D


  !> Mixes two vectors.
  subroutine TMixerCmplx_mix1D(this, qInpRes, qDiff)

    !> Mixer instance
    type(TMixerCmplx), intent(inout) :: this

    !> Input vector on entry, result vector on exit
    complex(dp), intent(inout) :: qInpRes(:)

    !> Difference between input and output vectors (measure of lack of convergence)
    complex(dp), intent(in) :: qDiff(:)

    select case (this%mixerType)
    case(mixerTypes%simple)
      call TSimpleMixerCmplx_mix(this%pSimpleMixer, qInpRes, qDiff)
    case(mixerTypes%broyden)
      call TBroydenMixerCmplx_mix(this%pBroydenMixer, qInpRes, qDiff)
    end select

  end subroutine TMixerCmplx_mix1D


#:for NAME, TYPE, LABEL in FLAVOURS
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

    !> Mixer instance.
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

    !> Mixer instance.
    type(TMixerReal), intent(inout) :: this

    !> Size of the vectors to mix.
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


  !> Return an inverse Jacobian if possible, halting if not
  subroutine TMixerReal_getInverseJacobian(this, invJac)

    !> Mixer instance.
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
