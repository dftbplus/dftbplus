!--------------------------------------------------------------------------------------------------!
!  DFTB+: general package for performing fast atomistic simulations                                !
!  Copyright (C) 2006 - 2023  DFTB+ developers group                                               !
!                                                                                                  !
!  See the LICENSE file for terms of usage and distribution.                                       !
!--------------------------------------------------------------------------------------------------!

#:include 'common.fypp'
#:set FLAVOURS = [('real', 'real', 'Real'),('cmplx', 'complex', 'Cmplx')]

module dftbp_mixer_factory
  use dftbp_common_accuracy, only : dp
  use dftbp_io_message, only : error
  !> Input Data structure
  use dftbp_mixer_mixer, only : TMixerReal, TMixerCmplx, mixerTypes, TMixerInput
  !> Mixer Implementations
  use dftbp_mixer_andersonmixer, only : TAndersonMixerReal, TAndersonMixerCmplx
  use dftbp_mixer_diismixer, only : TDiisMixerReal, TDiisMixerCmplx
  use dftbp_mixer_broydenmixer, only : TBroydenMixerReal,  TBroydenMixerCmplx
  use dftbp_mixer_simplemixer, only : TSimpleMixerReal, TSimpleMixerCmplx
  implicit none 

  contains

#:for NAME, TYPE, LABEL in FLAVOURS
    subroutine TMixerFactory${LABEL}$(pChrgMixer${LABEL}$, mixerInp)
      class(TMixer${LABEL}$), allocatable, intent(out) :: pChrgMixer${LABEL}$
      type(TMixerInput), intent(in) :: mixerInp

      !> Simple mixer (if used)
      type(TSimpleMixer${LABEL}$), allocatable :: pSimpleMixer${LABEL}$

      !> Anderson mixer (if used)
      type(TAndersonMixer${LABEL}$), allocatable :: pAndersonMixer${LABEL}$

      !> Broyden mixer (if used)
      type(TBroydenMixer${LABEL}$), allocatable :: pBroydenMixer${LABEL}$

      !> DIIS mixer (if used)
      type(TDiisMixer${LABEL}$), allocatable :: pDiisMixer${LABEL}$

      select case (mixerInp%iMixSwitch)
        case (mixerTypes%simple)
          allocate(pSimplemixer${LABEL}$)
          call move_alloc(pSimpleMixer${LABEL}$, pChrgMixer${LABEL}$)
        case (mixerTypes%anderson)
          allocate(pAndersonMixer${LABEL}$)
          call move_alloc(pAndersonMixer${LABEL}$, pChrgMixer${LABEL}$)
        case (mixerTypes%broyden)
          allocate(pBroydenMixer${LABEL}$)
          call move_alloc(pBroydenMixer${LABEL}$, pChrgMixer${LABEL}$)
        case (mixerTypes%diis)
          allocate(pDiisMixer${LABEL}$)
          call move_alloc(pDiisMixer${LABEL}$, pChrgMixer${LABEL}$)
        case default
          call error("Unknown charge/density mixer type.")
        end select
        call pChrgMixer${LABEL}$%init(mixerInp)
    end subroutine TMixerFactory${LABEL}$
#:endfor

end module dftbp_mixer_factory
