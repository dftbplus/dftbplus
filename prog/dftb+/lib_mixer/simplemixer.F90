!--------------------------------------------------------------------------------------------------!
!  DFTB+: general package for performing fast atomistic simulations                                !
!  Copyright (C) 2017  DFTB+ developers group                                                      !
!                                                                                                  !
!  See the LICENSE file for terms of usage and distribution.                                       !
!--------------------------------------------------------------------------------------------------!

!!* Simple mixer for mixing charges
module simplemixer
#include "assert.h"  
#include "allocate.h"  
  use accuracy
  implicit none
  
  private

  !!* Contains data for a simple mixer
  type OSimpleMixer
    private
    real(dp) :: mixParam       !* Mixing parameter
  end type OSimpleMixer


  !!* Creates a SimpleMixer instance
  interface create
    module procedure SimpleMixer_create
  end interface

  !!* Destroys a SimpleMixer instance
  interface destroy
    module procedure SimpleMixer_destroy
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
  public :: create, destroy, reset, mix


contains

  !!* Creates a simple mixer
  !!* @param self     Simple mixer instance on exit
  !!* @param mixParam Mixing parameter
  subroutine SimpleMixer_create(self, mixParam)
    type(OSimpleMixer), pointer :: self
    real(dp), intent(in) :: mixParam

    INITALLOCATE_P(self)
    self%mixParam = mixParam

  end subroutine SimpleMixer_create



  !!* Destroys the simple mixer
  !!* @param self Simple mixer to destroy.
  subroutine SimpleMixer_destroy(self)
    type(OSimpleMixer), pointer :: self

    DEALLOCATE_P(self)

  end subroutine SimpleMixer_destroy



  !!* Resets the mixer
  !!* @param self Simple mixer instance
  !!* @param nElem Length of the vectors to mix
  subroutine SimpleMixer_reset(self, nElem)
    type(OSimpleMixer), pointer :: self
    integer, intent(in) :: nElem

    ASSERT(nElem > 0)

    continue
    
  end subroutine SimpleMixer_reset

  

  !!* Does the actual mixing
  !!* @param self       SimpleMixer instance
  !!* @param qInpResult Input charge on entry, mixed charge on exit
  !!* @param qDiff      Charge difference
  subroutine SimpleMixer_mix(self, qInpResult, qDiff)
    type(OSimpleMixer), pointer :: self
    real(dp), intent(inout) :: qInpResult(:)
    real(dp), intent(in)    :: qDiff(:)

    ASSERT(size(qInpResult) == size(qDiff))

    !print *, "SIMPLE MIXER:"
    !print *, "QINPRESULT:"
    !print *, qInpResult
    !print *, "QDIFF:"
    !print *, qDiff
    qInpResult(:) = qInpResult(:) + self%mixParam * qDiff(:)
    !print *, "MIXED CHARGES:"
    !print *, qInpResult
    

  end subroutine SimpleMixer_mix


end module simplemixer
