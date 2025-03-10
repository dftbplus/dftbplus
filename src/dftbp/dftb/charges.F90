!--------------------------------------------------------------------------------------------------!
!  DFTB+: general package for performing fast atomistic simulations                                !
!  Copyright (C) 2006 - 2025  DFTB+ developers group                                               !
!                                                                                                  !
!  See the LICENSE file for terms of usage and distribution.                                       !
!--------------------------------------------------------------------------------------------------!

#:include 'common.fypp'

!> Module to calculate atomic charges
module dftbp_dftb_charges
  use dftbp_common_accuracy, only : dp
  use dftbp_dftb_uniquehubbard, only : TUniqueHubbard
  use dftbp_type_commontypes, only : TOrbitals
  implicit none

  private
  public :: getSummedCharges, getSummedChargesPerOrbital, getSummedChargesPerAtom
  public :: getSummedChargesPerLShell, getSummedChargesPerUniqueU


contains


  !> Calculates various gross charges.
  subroutine getSummedCharges(species, orb, qOrbital, q0, dQ, dQAtom, dQShell)

    !> Species of each atom.
    integer, intent(in) :: species(:)

    !> Orbital information
    type(TOrbitals), intent(in) :: orb

    !> Orbital populations. Shape [iChannels, mOrb, nAtom].
    real(dp), intent(in) :: qOrbital(:,:,:)

    !> Reference populations.
    real(dp), intent(in), optional :: q0(:,:,:)

    !> Charge per atomic orbital. Shape [qOrb, nAtom].
    real(dp), target, intent(out), optional :: dQ(:,:)

    !> Summed charge per atom.
    real(dp), intent(out), optional :: dQAtom(:)

    !> Summed charge per shell.
    real(dp), target, intent(out), optional :: dQShell(:,:)

    real(dp), allocatable, target :: dQLocal(:,:)
    real(dp), pointer :: dQWork(:,:)
    integer :: nAtom

    nAtom = size(orb%nOrbAtom)
    if (present(dQ)) then
      dQWork => dQ
    else
      allocate(dQLocal(orb%mOrb, nAtom))
      dQWork => dQLocal
    end if

    if (present(q0)) then
      call getSummedChargesPerOrbital(qOrbital(:,:,1), q0(:,:,1), dQWork)
    else
      dQWork(:,:) = qOrbital(:,:,1)
    end if
    if (present(dQAtom)) then
      call getSummedChargesPerAtom(dQWork, dQAtom)
    end if
    if (present(dQShell)) then
      call getSummedChargesPerLShell(species, orb, dQWork, dQShell)
    end if

  end subroutine getSummedCharges


  !> Orbital resolved charges
  pure subroutine getSummedChargesPerOrbital(qOrbital, q0, deltaQ)

    !> Charges per orbital
    real(dp), intent(in) :: qOrbital(:,:)

    !> Reference atomic charges
    real(dp), intent(in) :: q0(:,:)

    !> Summed charges (q - q0)
    real(dp), intent(out) :: deltaQ(:,:)

    deltaQ(:,:) = qOrbital - q0

  end subroutine getSummedChargesPerOrbital


  !> Atom resolved charges
  subroutine getSummedChargesPerAtom(deltaQ, deltaQAtom)

    !> Gross charge for all atomic orbitals on atoms
    real(dp), intent(in) :: deltaQ(:,:)

    !> Gross charge for each atom
    real(dp), intent(out) :: deltaQAtom(:)

    deltaQAtom(:) = sum(deltaQ, dim=1)

  end subroutine getSummedChargesPerAtom


  !> Shell resolved charges
  subroutine getSummedChargesPerLShell(species, orb, deltaQ, deltaQPerLShell)

    !> Chemical species of each atom
    integer, intent(in) :: species(:)

    !> Species resolved atomic orbital information
    type(TOrbitals), intent(in) :: orb

    !> Gross charge for each orbital
    real(dp), intent(in) :: deltaQ(:,:)

    !> Gross charge for each atomic shell
    real(dp), intent(out) :: deltaQPerLShell(:,:)

    integer :: iAt, iSp, iSh, iStart, iEnd

    deltaQPerLShell(:,:) = 0.0_dp
    do iAt = 1, size(orb%nOrbAtom)
      iSp = species(iAt)
      do iSh = 1, orb%nShell(iSp)
        iStart = orb%posShell(iSh,iSp)
        iEnd = orb%posShell(iSh+1,iSp) - 1
        deltaQPerLShell(iSh, iAt) = sum(deltaQ(iStart:iEnd, iAt))
      end do
    end do

  end subroutine getSummedChargesPerLShell


  !> Charges for regions with common U values
  subroutine getSummedChargesPerUniqueU(species, orb, hubbU, deltaQPerLShell, deltaQUniqU)

    !> Chemical species of each atom
    integer, intent(in) :: species(:)

    !> Species resolved atomic orbital information
    type(TOrbitals), intent(in) :: orb

    !> Contracted Hubbard U parameters
    type(TUniqueHubbard), intent(in) :: hubbU

    !> Charge per atomic shell
    real(dp), intent(in) :: deltaQPerLShell(:,:)

    !> Charges for atomic shells with a common U value
    real(dp), intent(out) :: deltaQUniqU(:,:)

    call hubbU%sumOverUniqueU(deltaQPerLShell, species, orb, deltaQUniqU)

  end subroutine getSummedChargesPerUniqueU

end module dftbp_dftb_charges
