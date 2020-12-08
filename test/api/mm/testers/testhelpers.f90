!--------------------------------------------------------------------------------------------------!
!  DFTB+: general package for performing fast atomistic simulations                                !
!  Copyright (C) 2006 - 2020  DFTB+ developers group                                               !
!                                                                                                  !
!  See the LICENSE file for terms of usage and distribution.                                       !
!--------------------------------------------------------------------------------------------------!

!> Module containing routines for the automatic testing of the API functionality

module testhelpers
  use, intrinsic :: iso_c_binding
  use dftbp_accuracy, only : dp
  use dftbp_taggedoutput, only : tagLabels, TTaggedWriter, TTaggedWriter_init
  implicit none
  private

  public :: writeAutotestTag, c_writeAutotestTag

contains

  !> Writes an autotest.tag file with the basic quantities
  subroutine writeAutotestTag(merminEnergy, gradients, stressTensor, &
      & grossCharges, extChargeGradients, atomMasses)

    !> Mermin energy
    real(dp), optional, intent(in) :: merminEnergy

    !> Gradients
    real(dp), optional, intent(in) :: gradients(:,:)

    !> Stress tensor of the periodic system.
    real(dp), optional, intent(in) :: stressTensor(:,:)

    !> Gross mulliken charges
    real(dp), optional, intent(in) :: grossCharges(:)

    !> Gradients on external charges.
    real(dp), optional, intent(in) :: extChargeGradients(:,:)

    !> Atomic masses
    real(dp), optional, intent(in) :: atomMasses(:)

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
    if(present(stressTensor)) then
      call taggedWriter%write(autotestTag, tagLabels%stressTot, stressTensor)
    end if
    if(present(atomMasses)) then
      call taggedWriter%write(autotestTag, tagLabels%atomMass, atomMasses)
    end if
    close(autotestTag)

  end subroutine writeAutotestTag


  !> C wrapper for the write autotest tag routine.
  subroutine c_writeAutotestTag(nAtom, nExtCharge, merminEnergy, gradients, stressTensor, &
      & grossCharges, extChargeGradients) bind(C, name='dftbp_write_autotest_tag')

    !> Number of atoms
    integer(c_int), intent(in), value :: nAtom

    !> Number of external charges (set it to zero, if none)
    integer(c_int), intent(in), value :: nExtCharge

    !> Mermin energy
    real(c_double), intent(in), value :: merminEnergy

    !> Gradients or null pointer, if not avaialable.
    type(c_ptr), intent(in), value :: gradients

    !> Stress tensor or null pointer, if not avaialable.
    type(c_ptr), intent(in), value :: stressTensor

    !> Gross Mulliken charges or null pointer, if not avaialable.
    type(c_ptr), intent(in), value :: grossCharges

    !> Gradients on the external charges or null pointer, if not avaialable.
    type(c_ptr), intent(in), value :: extChargeGradients

    real(dp), pointer :: pGradients(:,:), pGrossCharges(:), &
      & pExtChargeGradients(:,:), pStressTensor(:,:)

    if (c_associated(gradients)) then
      call c_f_pointer(gradients, pGradients, [3, nAtom])
    else
      pGradients => null()
    end if
    
    if (c_associated(stressTensor)) then
      call c_f_pointer(stressTensor, pStressTensor, [3, 3])
    else
      pStressTensor => null()
    end if

    if (c_associated(grossCharges)) then
      call c_f_pointer(grossCharges, pGrossCharges, [nAtom])
    else
      pGrossCharges => null()
    end if

    if (nExtCharge > 0 .and. c_associated(extChargeGradients)) then
      call c_f_pointer(extChargeGradients, pExtChargeGradients, [3, nExtCharge])
    else
      pExtChargeGradients => null()
    end if

    call writeAutotestTag(merminEnergy=merminEnergy, gradients=pGradients, &
        & stressTensor=pStressTensor, grossCharges=pGrossCharges, &
        & extChargeGradients=pExtChargeGradients)

  end subroutine c_writeAutotestTag


end module testhelpers
