!--------------------------------------------------------------------------------------------------!
!  DFTB+: general package for performing fast atomistic simulations                                !
!  Copyright (C) 2018  DFTB+ developers group                                                      !
!                                                                                                  !
!  See the LICENSE file for terms of usage and distribution.                                       !
!--------------------------------------------------------------------------------------------------!

!> Module containing routines for the automatic testing of the API functionality

module testhelpers
  use dftbp_accuracy, only : dp
  use dftbp_taggedoutput, only : tagLabels, TTaggedWriter, TTaggedWriter_init
  implicit none
  private

  public :: writeAutotestTag


contains

  !> Writes an autotest.tag file with the basic quantities
  subroutine writeAutotestTag(merminEnergy, gradients, grossCharges, extChargeGradients)

    real(dp), optional, intent(in) :: merminEnergy
    real(dp), optional, intent(in) :: gradients(:,:)
    real(dp), optional, intent(in) :: grossCharges(:)
    real(dp), optional, intent(in) :: extChargeGradients(:,:)

    type(TTaggedWriter) :: taggedWriter
    integer :: autotestTag

    ! Write out quantities in tagged form for the internal testing system
    call TTaggedWriter_init(taggedWriter)
    open(newunit=autotestTag, file="autotest.tag", action="write")
    if (present(merminEnergy)) then
      call taggedWriter%write(autotestTag, tagLabels%freeEgy, merminEnergy)
    end if
    if (present(gradients)) then
      call taggedWriter%write(autotestTag, tagLabels%forceTot, -gradients)
    end if
    if (present(grossCharges)) then
      call taggedWriter%write(autotestTag, tagLabels%qOutAtGross, grossCharges)
    end if
    if (present(extChargeGradients)) then
      call taggedWriter%write(autotestTag, tagLabels%chrgForces, -extChargeGradients)
    end if
    close(autotestTag)

  end subroutine writeAutotestTag


end module testhelpers
