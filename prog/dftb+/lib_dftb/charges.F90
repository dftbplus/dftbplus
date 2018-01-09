!--------------------------------------------------------------------------------------------------!
!  DFTB+: general package for performing fast atomistic simulations                                !
!  Copyright (C) 2017  DFTB+ developers group                                                      !
!                                                                                                  !
!  See the LICENSE file for terms of usage and distribution.                                       !
!--------------------------------------------------------------------------------------------------!

#:include 'common.fypp'

!> Module to calculate atomic charges
module charges
  use assert
  use accuracy
  use commontypes, only : TOrbitals
  implicit none
  private

  public :: getSummedCharges

contains


  !> Calculates various gross charges.
  subroutine getSummedCharges(species, orb, qOrbital, q0, iHubbU, dQ, dQAtom, dQShell, dQUniqU)


    !> Species of each atom.
    integer, intent(in) :: species(:)


    !> Orbital information
    type(TOrbitals), intent(in) :: orb


    !> Orbital populations.
    real(dp), intent(in) :: qOrbital(:,:,:)


    !> Reference populations.
    real(dp), intent(in) :: q0(:,:,:)


    !> Unique Hubbard U indices (only needed if dQUniqU is passed)
    integer, intent(in), optional :: iHubbU(:,:)


    !> charge per orbital.
    real(dp), target, intent(out), optional :: dQ(:,:)


    !> Summed charge per atom.
    real(dp), target, intent(out), optional :: dQAtom(:)


    !> Summed charge per shell.
    real(dp), target, intent(out), optional :: dQShell(:,:)


    !> Summed charges per unique Hubbard U
    real(dp), target, intent(out), optional :: dQUniqU(:,:)

    real(dp), allocatable, target :: dQLocal(:,:), dQShellLocal(:,:)
    real(dp), pointer :: dQWork(:,:), dQShellWork(:,:)
    integer :: nAtom

    @:ASSERT(present(iHubbU) .eqv. present(dQUniqU))

    nAtom = size(orb%nOrbAtom)
    if (present(dQ)) then
      dQWork => dQ
    else
      allocate(dQLocal(orb%mOrb, nAtom))
      dQWork => dQLocal
    end if
    if (present(dQUniqU) .or. present(dQShell)) then
      if (present(dQShell)) then
        dQShellWork => dQShell
      else
        allocate(dQShellLocal(orb%mShell, nAtom))
        dQShellWork => dQShellLocal
      end if
    end if

    call getSummedChargesPerOrbital(qOrbital(:,:,1), q0(:,:,1), dQWork)
    if (present(dQAtom)) then
      call getSummedChargesPerAtom(dQWork, dQAtom)
    end if
    if (present(dQShell) .or. present(dQUniqU)) then
      call getSummedChargesPerLShell(species, orb, dQWork, dQShellWork)
    end if
    if (present(dQUniqU)) then
      call getSummedChargesPerUniqU(species, orb, dQShellWork, iHubbU, dQUniqU)
    end if

  end subroutine getSummedCharges

  ! Private routines


  !> orbital resolved charges
  subroutine getSummedChargesPerOrbital(qOrbital, q0, deltaQ)

    !> charges per orbital
    real(dp), intent(in) :: qOrbital(:,:)

    !> reference atomic charges
    real(dp), intent(in) :: q0(:,:)

    !> Summed charges (q - q0)
    real(dp), intent(out) :: deltaQ(:,:)

    deltaQ(:,:) = qOrbital - q0

  end subroutine getSummedChargesPerOrbital


  !> atom resolved charges
  subroutine getSummedChargesPerAtom(deltaQ, deltaQAtom)

    !> gross charge for all atomic orbitals on atoms
    real(dp), intent(in) :: deltaQ(:,:)

    !> gross charge for each atom
    real(dp), intent(out) :: deltaQAtom(:)

    deltaQAtom(:) = sum(deltaQ, dim=1)

  end subroutine getSummedChargesPerAtom


  !> shell resolved charges
  subroutine getSummedChargesPerLShell(species, orb, deltaQ, deltaQPerLShell)

    !> chemical species of each atom
    integer, intent(in) :: species(:)

    !> species resolved atomic orbital information
    type(TOrbitals), intent(in) :: orb

    !> gross charge for each orbital
    real(dp), intent(in) :: deltaQ(:,:)

    !> gross charge for each atomic shell
    real(dp), intent(out) :: deltaQPerLShell(:,:)

    integer :: iAt, iSp, iSh, iStart, iend

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


  !> charges for regions with common U values
  subroutine getSummedChargesPerUniqU(species, orb, deltaQPerLShell, iHubbU, deltaQUniqU)

    !> chemical species of each atom
    integer, intent(in) :: species(:)

    !> species resolved atomic orbital information
    type(TOrbitals), intent(in) :: orb

    !> charge per atomic shell
    real(dp), intent(in) :: deltaQPerLShell(:,:)

    !> index for unique Hubbard U values
    integer, intent(in) :: iHubbU(:,:)

    !> charges for atomic shells with a common U value
    real(dp), intent(out) :: deltaQUniqU(:,:)

    integer :: iAt, iSp, iSh

    deltaQUniqU(:,:) = 0.0_dp
    do iAt = 1, size(orb%nOrbAtom)
      iSp = species(iAt)
      do iSh = 1, orb%nShell(iSp)
        deltaQUniqU(iHubbU(iSh, iSp), iAt) = deltaQUniqU(iHubbU(iSh, iSp), iAt)&
            & + deltaQPerLShell(iSh, iAt)
      end do
    end do

  end subroutine getSummedChargesPerUniqU

end module charges
