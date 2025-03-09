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
  use dftbp_mixer_andersonmixer, only : TAndersonMixerCmplx, TAndersonMixerCmplx_init, &
      & TAndersonMixerInp, TAndersonMixerReal, TAndersonMixerReal_init
  use dftbp_mixer_broydenmixer, only : TBroydenMixerCmplx, TBroydenMixerCmplx_init, &
      & TBroydenMixerInp, TBroydenMixerReal, TBroydenMixerReal_init
  use dftbp_mixer_diismixer, only : TDiisMixerCmplx, TDiisMixerCmplx_init, &
      & TDiisMixerInp, TDiisMixerReal, TDiisMixerReal_init
  use dftbp_mixer_mixer, only : TMixerCmplx, TMixerReal
  use dftbp_mixer_simplemixer, only : TSimpleMixerCmplx, TSimpleMixerCmplx_init, &
      & TSimpleMixerInp, TSimpleMixerReal, TSimpleMixerReal_init

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
  subroutine TMixerFactory${LABEL}$(chrgMixer${LABEL}$, mixerInp)
    !> The resulting allocated and initialised mixer
    class(TMixer${LABEL}$), allocatable, intent(out) :: chrgMixer${LABEL}$

    !> Mixer Input data structure, containing one of the specific mixer inputs
    type(TMixerInput), intent(in) :: mixerInp

    type(TAndersonMixer${LABEL}$), allocatable :: andersonMixer${LABEL}$
    type(TBroydenMixer${LABEL}$), allocatable :: broydenMixer${LABEL}$
    type(TDiisMixer${LABEL}$), allocatable :: diisMixer${LABEL}$
    type(TSimpleMixer${LABEL}$), allocatable :: simpleMixer${LABEL}$

    if (allocated(mixerInp%andersonMixerInp)) then
      allocate(andersonMixer${LABEL}$)
      call TAndersonMixer${LABEL}$_init(andersonMixer${LABEL}$, mixerInp%andersonMixerInp)
      call move_alloc(andersonMixer${LABEL}$, chrgMixer${LABEL}$)
    else if (allocated(mixerInp%broydenMixerInp)) then
      allocate(broydenMixer${LABEL}$)
      call TBroydenMixer${LABEL}$_init(broydenMixer${LABEL}$, mixerInp%broydenMixerInp)
      call move_alloc(broydenMixer${LABEL}$, chrgMixer${LABEL}$)
    else if (allocated(mixerInp%diisMixerInp)) then
      allocate(diisMixer${LABEL}$)
      call TDiisMixer${LABEL}$_init(diisMixer${LABEL}$, mixerInp%diisMixerInp)
      call move_alloc(diisMixer${LABEL}$, chrgMixer${LABEL}$)
    else if (allocated(mixerInp%simpleMixerInp)) then
      allocate(simpleMixer${LABEL}$)
      call TSimpleMixer${LABEL}$_init(simpleMixer${LABEL}$, mixerInp%simpleMixerInp)
      call move_alloc(simpleMixer${LABEL}$, chrgMixer${LABEL}$)
    else
      call error("TMixerFactory${LABEL}$ : No mixer input provided")
    end if
  end subroutine TMixerFactory${LABEL}$
#:endfor

end module dftbp_mixer_factory
