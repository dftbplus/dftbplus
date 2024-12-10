!--------------------------------------------------------------------------------------------------!
!  DFTB+: general package for performing fast atomistic simulations                                !
!  Copyright (C) 2006 - 2023  DFTB+ developers group                                               !
!                                                                                                  !
!  See the LICENSE file for terms of usage and distribution.                                       !
!--------------------------------------------------------------------------------------------------!

#:include 'common.fypp'
#:set FLAVOURS = [('real', 'real', 'Real'),('cmplx', 'complex', 'Cmplx')]
#:set VARIANTS = ['Anderson', 'Broyden', 'Diis', 'Simple']
module dftbp_mixer_factory
  use dftbp_io_message, only : error

  use dftbp_mixer_mixer, only : TMixerReal, TMixerCmplx
  #:for VARIANT in VARIANTS
  use dftbp_mixer_${VARIANT}$mixer, only : T${VARIANT}$MixerInp
    #:for NAME, TYPE, LABEL in FLAVOURS
    use dftbp_mixer_${VARIANT}$mixer, only : T${VARIANT}$Mixer${LABEL}$, T${VARIANT}$Mixer${LABEL}$_init
    #:endfor
  #:endfor

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
        call TAndersonMixer${LABEL}$_init(pAndersonMixer${LABEL}$, mixerInp%andersonMixerInp)
        call move_alloc(pAndersonMixer${LABEL}$, pChrgMixer${LABEL}$)
      else if (allocated(mixerInp%broydenMixerInp)) then
        call TBroydenMixer${LABEL}$_init(pBroydenMixer${LABEL}$, mixerInp%broydenMixerInp)
        call move_alloc(pBroydenMixer${LABEL}$, pChrgMixer${LABEL}$)
      else if (allocated(mixerInp%diisMixerInp)) then
        call TDiisMixer${LABEL}$_init(pDiisMixer${LABEL}$, mixerInp%diisMixerInp)
        call move_alloc(pDiisMixer${LABEL}$, pChrgMixer${LABEL}$)
      else if (allocated(mixerInp%simpleMixerInp)) then
        call TSimpleMixer${LABEL}$_init(pSimpleMixer${LABEL}$, mixerInp%simpleMixerInp)
        call move_alloc(pSimpleMixer${LABEL}$, pChrgMixer${LABEL}$)
      else
        call error("TMixerFactory${LABEL}$ : No mixer input provided")
      end if
    end subroutine TMixerFactory${LABEL}$
    #:endfor

end module dftbp_mixer_factory
