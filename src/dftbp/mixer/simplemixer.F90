!--------------------------------------------------------------------------------------------------!
!  DFTB+: general package for performing fast atomistic simulations                                !
!  Copyright (C) 2006 - 2025  DFTB+ developers group                                               !
!                                                                                                  !
!  See the LICENSE file for terms of usage and distribution.                                       !
!--------------------------------------------------------------------------------------------------!

#:include 'common.fypp'
#:set FLAVOURS = [('cmplx', 'complex', 'Cmplx'), ('real', 'real', 'Real')]

!> Simple mixer for mixing charges
module dftbp_mixer_simplemixer
  use dftbp_common_accuracy, only : dp
  use dftbp_mixer_mixer, only : TMixerCmplx, TMixerReal
  implicit none


  private
  public :: TSimpleMixerInp
  public :: TSimpleMixerReal, TSimpleMixerReal_init
  public :: TSimpleMixerCmplx, TSimpleMixerCmplx_init

  type :: TSimpleMixerInp
    !> Mixing parameter
    real(dp) :: mixParam = 0.0_dp
  end type TSimpleMixerInp

#:for NAME, TYPE, LABEL in FLAVOURS
  !> Contains data for a simple mixer
  type, extends (TMixer${LABEL}$) :: TSimpleMixer${LABEL}$
    private
    !> Mixing parameter
    real(dp) :: mixParam

  contains
    procedure :: reset => TSimpleMixer${LABEL}$_reset
    procedure :: mix1D => TSimpleMixer${LABEL}$_mix
  end type TSimpleMixer${LABEL}$
#:endfor


contains

#:for NAME, TYPE, LABEL in FLAVOURS

  !> Initializes a simple mixer.
  subroutine TSimpleMixer${LABEL}$_init(this, mixerInp)

    !> Simple mixer instance on exit
    type(TSimpleMixer${LABEL}$), intent(out) :: this

    !> TSimpleMixerInp input data struct
    type(TSimpleMixerInp), intent(in) :: mixerInp

    this%mixParam = mixerInp%mixParam

  end subroutine TSimpleMixer${LABEL}$_init


  !> Resets the mixer.
  subroutine TSimpleMixer${LABEL}$_reset(this, nElem)

    !> Simple mixer instance
    class(TSimpleMixer${LABEL}$), intent(inout) :: this

    !> Length of the vectors to mix
    integer, intent(in) :: nElem

    @:ASSERT(nElem > 0)

    continue

  end subroutine TSimpleMixer${LABEL}$_reset


  !> Does the actual mixing.
  subroutine TSimpleMixer${LABEL}$_mix(this, qInpResult, qDiff)

    !> SimpleMixer instance
    class(TSimpleMixer${LABEL}$), intent(inout) :: this

    !> Input charge on entry, mixed charge on exit
    ${TYPE}$(dp), intent(inout) :: qInpResult(:)

    !> Charge difference
    ${TYPE}$(dp), intent(in) :: qDiff(:)

    @:ASSERT(size(qInpResult) == size(qDiff))

    qInpResult(:) = qInpResult + this%mixParam * qDiff

  end subroutine TSimpleMixer${LABEL}$_mix
#:endfor

end module dftbp_mixer_simplemixer
