!--------------------------------------------------------------------------------------------------!
!  DFTB+: general package for performing fast atomistic simulations                                !
!  Copyright (C) 2006 - 2023  DFTB+ developers group                                               !
!                                                                                                  !
!  See the LICENSE file for terms of usage and distribution.                                       !
!--------------------------------------------------------------------------------------------------!

!> Module containing routines for the automatic testing of the API functionality

module testhelpers
  use, intrinsic :: iso_c_binding
  use dftbp_common_accuracy, only : dp
  use dftbp_common_file, only : TFileDescr, openFile, closeFile
  use dftbp_io_taggedoutput, only : tagLabels, TTaggedWriter, TTaggedWriter_init
  implicit none
  private

  public :: writeAutotestTag, c_writeAutotestTag

contains

  !> Writes an autotest.tag file with the basic quantities
  subroutine writeAutotestTag(merminEnergy, gradients, stressTensor, &
      & grossCharges, extChargeGradients, tdDipole, tdEnergy, tdCharges, &
      & tdCoords, tdForces, atomMasses, potential, cm5Charges, cutOff)

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

    !> Time-dependent dipole moment
    real(dp), optional, intent(in) :: tdDipole(:,:)

    !> Time-dependent total energy
    real(dp), optional, intent(in) :: tdEnergy

    !> Time-dependent negative gross charge
    real(dp), optional, intent(in) :: tdCharges(:,:)

    !> Time-dependent coordinates in Ehrenfest dynamics
    real(dp), optional, intent(in) :: tdCoords(:,:)

    !> Time-dependent Ehrenfest forces
    real(dp), optional, intent(in) :: tdForces(:,:)

    !> Atomic masses
    real(dp), optional, intent(in) :: atomMasses(:)

    !> Electrostatic potential in points
    real(dp), optional, intent(in) :: potential(:)

    !> Gross CM5 charges
    real(dp), optional, intent(in) :: cm5Charges(:)

    !> Cutoff distance
    real(dp), optional, intent(in) :: cutOff

    type(TTaggedWriter) :: taggedWriter
    type(TFileDescr) :: autotestTag

    ! Write out quantities in tagged form for the internal testing system
    call TTaggedWriter_init(taggedWriter)
    call openFile(autotestTag, "autotest.tag", mode="w")
    if (present(merminEnergy)) then
      call taggedWriter%write(autotestTag%unit, tagLabels%freeEgy, merminEnergy)
    end if
    if (present(gradients)) then
      call taggedWriter%write(autotestTag%unit, tagLabels%forceTot, -gradients)
    end if
    if (present(grossCharges)) then
      call taggedWriter%write(autotestTag%unit, tagLabels%qOutAtGross, grossCharges)
    end if
    if (present(extChargeGradients)) then
      call taggedWriter%write(autotestTag%unit, tagLabels%chrgForces, -extChargeGradients)
    end if
    if(present(stressTensor)) then
      call taggedWriter%write(autotestTag%unit, tagLabels%stressTot, stressTensor)
    end if
    if (present(tdEnergy)) then
      call taggedWriter%write(autotestTag%unit, tagLabels%tdenergy, tdEnergy)
    end if
    if (present(tdDipole)) then
      call taggedWriter%write(autotestTag%unit, tagLabels%tddipole, tdDipole)
    end if
    if (present(tdCharges)) then
      call taggedWriter%write(autotestTag%unit, tagLabels%tdcharges, tdCharges)
    end if
    if (present(tdCoords)) then
      call taggedWriter%write(autotestTag%unit, tagLabels%ehrencoords, tdCoords)
    end if
    if (present(tdForces)) then
      call taggedWriter%write(autotestTag%unit, tagLabels%ehrenforces, tdForces)
    end if
    if(present(atomMasses)) then
      call taggedWriter%write(autotestTag%unit, tagLabels%atomMass, atomMasses)
    end if
    if(present(potential)) then
      call taggedWriter%write(autotestTag%unit, tagLabels%internField, potential)
    end if
    if (present(cm5Charges)) then
      call taggedWriter%write(autotestTag%unit, tagLabels%qOutAtCM5, cm5Charges)
    end if
    if (present(cutOff)) then
      call taggedWriter%write(autotestTag%unit, "cutoff", cutOff)
    end if
    call closeFile(autotestTag)

  end subroutine writeAutotestTag


  !> C wrapper for the write autotest tag routine.
  subroutine c_writeAutotestTag(nAtom, nExtCharge, nPotLocations, merminEnergy, gradients,&
      & stressTensor, grossCharges, extChargeGradients, potential, cm5Charges)&
      & bind(C, name='dftbp_write_autotest_tag')

    !> Number of atoms
    integer(c_int), intent(in), value :: nAtom

    !> Number of external charges (set it to zero, if none)
    integer(c_int), intent(in), value :: nExtCharge

    !> Number of requeseted potential points
    integer(c_int), intent(in), value :: nPotLocations

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

    !> Electrostatic potential in nPotLocations points
    type(c_ptr), intent(in), value :: potential

    !> Gross CM5 charges or null pointer, if not avaialable.
    type(c_ptr), intent(in), value :: cm5Charges

    real(dp), pointer :: pGradients(:,:), pGrossCharges(:), &
      & pExtChargeGradients(:,:), pStressTensor(:,:), pPotential(:), &
      & pCM5Charges(:)

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

    if (nPotLocations > 0 .and. c_associated(potential)) then
      call c_f_pointer(potential, pPotential, [nPotLocations])
    else
      pPotential => null()
    end if

    if (c_associated(cm5Charges)) then
      call c_f_pointer(cm5Charges, pCM5Charges, [nAtom])
    else
      pCM5Charges => null()
    end if

    call writeAutotestTag(merminEnergy=merminEnergy, gradients=pGradients, &
        & stressTensor=pStressTensor, grossCharges=pGrossCharges, &
        & extChargeGradients=pExtChargeGradients, potential=pPotential, &
        & cm5Charges=pCM5charges)

  end subroutine c_writeAutotestTag


end module testhelpers
