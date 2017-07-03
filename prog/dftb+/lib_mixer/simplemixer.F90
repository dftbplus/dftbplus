!--------------------------------------------------------------------------------------------------!
!  DFTB+: general package for performing fast atomistic simulations                                !
!  Copyright (C) 2017  DFTB+ developers group                                                      !
!                                                                                                  !
!  See the LICENSE file for terms of usage and distribution.                                       !
!--------------------------------------------------------------------------------------------------!

#:include 'common.fypp'

!!* Simple mixer for mixing charges
module simplemixer
  use assert
  use accuracy
  implicit none

  private

  !!* Contains data for a simple mixer
  type OSimpleMixer
    private
    real(dp) :: mixParam       !* Mixing parameter
  end type OSimpleMixer


  !!* Creates a SimpleMixer instance
  interface init
    module procedure SimpleMixer_init
  end interface

  !!* Resets a SimpleMixer
  interface reset
    module procedure SimpleMixer_reset
  end interface

  !!* Does the simple mixing
  interface mix
    module procedure SimpleMixer_mix
  end interface


  public :: OSimpleMixer
  public :: init, reset, mix


contains

  !!* Creates a simple mixer
  !!* @param self     Simple mixer instance on exit
  !!* @param mixParam Mixing parameter
  subroutine SimpleMixer_init(self, mixParam)
    type(OSimpleMixer), intent(out) :: self
    real(dp), intent(in) :: mixParam

    self%mixParam = mixParam

  end subroutine SimpleMixer_init


  !!* Resets the mixer
  !!* @param self Simple mixer instance
  !!* @param nElem Length of the vectors to mix
  subroutine SimpleMixer_reset(self, nElem)
    type(OSimpleMixer), intent(inout) :: self
    integer, intent(in) :: nElem

    @:ASSERT(nElem > 0)

    continue

  end subroutine SimpleMixer_reset



  !!* Does the actual mixing
  !!* @param self       SimpleMixer instance
  !!* @param qInpResult Input charge on entry, mixed charge on exit
  !!* @param qDiff      Charge difference
  subroutine SimpleMixer_mix(self, qInpResult, qDiff)
    type(OSimpleMixer), intent(inout) :: self
    real(dp), intent(inout) :: qInpResult(:)
    real(dp), intent(in)    :: qDiff(:)

    @:ASSERT(size(qInpResult) == size(qDiff))

    qInpResult(:) = qInpResult(:) + self%mixParam * qDiff(:)

  end subroutine SimpleMixer_mix


end module simplemixer
