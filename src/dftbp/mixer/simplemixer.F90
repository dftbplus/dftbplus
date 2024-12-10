!--------------------------------------------------------------------------------------------------!
!  DFTB+: general package for performing fast atomistic simulations                                !
!  Copyright (C) 2006 - 2023  DFTB+ developers group                                               !
!                                                                                                  !
!  See the LICENSE file for terms of usage and distribution.                                       !
!--------------------------------------------------------------------------------------------------!

#:include 'common.fypp'
#:set FLAVOURS = [('cmplx', 'complex', 'Cmplx'), ('real', 'real', 'Real')]

!> Simple mixer for mixing charges
module dftbp_mixer_simplemixer
  use dftbp_common_accuracy, only : dp
  #:for NAME, TYPE, LABEL in FLAVOURS
  use dftbp_mixer_mixer, only: TMixer${LABEL}$
  #:endfor
  implicit none


  private
  public :: TSimpleMixerInp
#:for NAME, TYPE, LABEL in FLAVOURS
  public :: TSimpleMixer${LABEL}$, TSimpleMixer${LABEL}$_init
#:endfor

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

  !> Creates a simple mixer.
  subroutine TSimpleMixer${LABEL}$_init(this, mixerInp)

    !> Simple mixer instance on exit
    type(TSimpleMixer${LABEL}$), allocatable, intent(out) :: this
    
    !> TSimpleMixerInp input data struct
    type(TSimpleMixerInp), intent(in) :: mixerInp

    allocate(TSimpleMixer${LABEL}$ :: this)

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
