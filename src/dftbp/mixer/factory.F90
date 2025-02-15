!--------------------------------------------------------------------------------------------------!
!  DFTB+: general package for performing fast atomistic simulations                                !
!  Copyright (C) 2006 - 2023  DFTB+ developers group                                               !
!                                                                                                  !
!  See the LICENSE file for terms of usage and distribution.                                       !
!--------------------------------------------------------------------------------------------------!

#:include 'common.fypp'
#:set FLAVOURS = [('real', 'real', 'Real'),('cmplx', 'complex', 'Cmplx')]
module dftbp_mixer_factory
  use dftbp_io_message, only : error
  use dftbp_mixer_mixer, only : TMixerReal, TMixerCmplx
  ! Mixer Input
  use dftbp_mixer_andersonmixer, only : TAndersonMixerInp
  use dftbp_mixer_broydenmixer, only : TBroydenMixerInp
  use dftbp_mixer_diismixer, only : TDiisMixerInp
  use dftbp_mixer_simplemixer, only : TSimpleMixerInp
  ! Real Mixer Variants
  use dftbp_mixer_andersonmixer, only : TAndersonMixerReal, TAndersonMixerReal_init
  use dftbp_mixer_broydenmixer, only : TBroydenMixerReal, TBroydenMixerReal_init
  use dftbp_mixer_diismixer, only : TDiisMixerReal, TDiisMixerReal_init
  use dftbp_mixer_simplemixer, only : TSimpleMixerReal, TSimpleMixerReal_init
  ! Complex Mixer Variants
  use dftbp_mixer_andersonmixer, only : TAndersonMixerCmplx, TAndersonMixerCmplx_init
  use dftbp_mixer_broydenmixer, only : TBroydenMixerCmplx, TBroydenMixerCmplx_init
  use dftbp_mixer_diismixer, only : TDiisMixerCmplx, TDiisMixerCmplx_init
  use dftbp_mixer_simplemixer, only : TSimpleMixerCmplx, TSimpleMixerCmplx_init

  implicit none
  private
  public :: TMixerInput, TMixerFactoryReal, TMixerFactoryCmplx


  type :: TMixerInput
    type(TAndersonMixerInp), allocatable :: andersonMixerInp
    type(TBroydenMixerInp), allocatable :: broydenMixerInp
    type(TDiisMixerInp), allocatable :: diisMixerInp
    type(TSimpleMixerInp), allocatable :: simpleMixerInp
  end type TMixerInput


contains

  #:for NAME, TYPE, LABEL in FLAVOURS
  subroutine TMixerFactory${LABEL}$(pChrgMixer${LABEL}$, mixerInp)
    class(TMixer${LABEL}$), allocatable, intent(out) :: pChrgMixer${LABEL}$
    type(TMixerInput), intent(in) :: mixerInp

    !> Anderson mixer (if used)
    type(TAndersonMixer${LABEL}$), allocatable :: pAndersonMixer${LABEL}$

    !> Broyden mixer (if used)
    type(TBroydenMixer${LABEL}$), allocatable :: pBroydenMixer${LABEL}$

    !> DIIS mixer (if used)
    type(TDiisMixer${LABEL}$), allocatable :: pDiisMixer${LABEL}$

    !> Simple mixer (if used)
    type(TSimpleMixer${LABEL}$), allocatable :: pSimpleMixer${LABEL}$

    if (allocated(mixerInp%andersonMixerInp)) then
      allocate(pAndersonMixer${LABEL}$)
      call TAndersonMixer${LABEL}$_init(pAndersonMixer${LABEL}$, mixerInp%andersonMixerInp)
      call move_alloc(pAndersonMixer${LABEL}$, pChrgMixer${LABEL}$)
    else if (allocated(mixerInp%broydenMixerInp)) then
      allocate(pBroydenMixer${LABEL}$)
      call TBroydenMixer${LABEL}$_init(pBroydenMixer${LABEL}$, mixerInp%broydenMixerInp)
      call move_alloc(pBroydenMixer${LABEL}$, pChrgMixer${LABEL}$)
    else if (allocated(mixerInp%diisMixerInp)) then
      allocate(pDiisMixer${LABEL}$)
      call TDiisMixer${LABEL}$_init(pDiisMixer${LABEL}$, mixerInp%diisMixerInp)
      call move_alloc(pDiisMixer${LABEL}$, pChrgMixer${LABEL}$)
    else if (allocated(mixerInp%simpleMixerInp)) then
      allocate(pSimpleMixer${LABEL}$)
      call TSimpleMixer${LABEL}$_init(pSimpleMixer${LABEL}$, mixerInp%simpleMixerInp)
      call move_alloc(pSimpleMixer${LABEL}$, pChrgMixer${LABEL}$)
    else
      call error("TMixerFactory${LABEL}$ : No mixer input provided")
    end if
  end subroutine TMixerFactory${LABEL}$
  #:endfor

end module dftbp_mixer_factory
