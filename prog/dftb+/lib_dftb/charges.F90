!--------------------------------------------------------------------------------------------------!
!  DFTB+: general package for performing fast atomistic simulations                                !
!  Copyright (C) 2017  DFTB+ developers group                                                      !
!                                                                                                  !
!  See the LICENSE file for terms of usage and distribution.                                       !
!--------------------------------------------------------------------------------------------------!

#:include 'common.fypp'

module charges
  use assert
  use accuracy
  use commontypes, only : TOrbitals
  implicit none
  private

  public :: getNetCharges

contains

  !> Calculates various net charges.
  !!
  subroutine getNetCharges(species, orb, qOrbital, q0, iHubbU, dQ, dQAtom,&
      & dQShell, dQUniqU)

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

    !> Net charge per orbital.
    real(dp), target, intent(out), optional :: dQ(:,:)

    !> Net charge per atom.
    real(dp), target, intent(out), optional :: dQAtom(:)

    !> Net charge per shell.
    real(dp), target, intent(out), optional :: dQShell(:,:)

    !> Net charges per unique Hubbard U
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

    call getNetChargesPerOrbital(qOrbital(:,:,1), q0(:,:,1), dQWork)
    if (present(dQAtom)) then
      call getNetChargesPerAtom(dQWork, dQAtom)
    end if
    if (present(dQShell) .or. present(dQUniqU)) then
      call getNetChargesPerLShell(species, orb, dQWork, dQShellWork)
    end if
    if (present(dQUniqU)) then
      call getNetChargesPerUniqU(species, orb, dQShellWork, iHubbU, dQUniqU)
    end if

  end subroutine getNetCharges


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!  Private routines
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine getNetChargesPerOrbital(qOrbital, q0, deltaQ)
    real(dp), intent(in) :: qOrbital(:,:), q0(:,:)
    real(dp), intent(out) :: deltaQ(:,:)

    deltaQ(:,:) = qOrbital - q0

  end subroutine getNetChargesPerOrbital


  subroutine getNetChargesPerAtom(deltaQ, deltaQAtom)
    real(dp), intent(in) :: deltaQ(:,:)
    real(dp), intent(out) :: deltaQAtom(:)

    deltaQAtom(:) = sum(deltaQ, dim=1)

  end subroutine getNetChargesPerAtom


  subroutine getNetChargesPerLShell(species, orb, deltaQ, deltaQPerLShell)
    integer, intent(in) :: species(:)
    type(TOrbitals), intent(in) :: orb
    real(dp), intent(in) :: deltaQ(:,:)
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

  end subroutine getNetChargesPerLShell


  subroutine getNetChargesPerUniqU(species, orb, deltaQPerLShell, iHubbU,&
      & deltaQUniqU)
    integer, intent(in) :: species(:)
    type(TOrbitals), intent(in) :: orb
    real(dp), intent(in) :: deltaQPerLShell(:,:)
    integer, intent(in) :: iHubbU(:,:)
    real(dp), intent(out) :: deltaQUniqU(:,:)

    integer :: iAt, iSp, iSh

    deltaQUniqU(:,:) = 0.0_dp
    do iAt = 1, size(orb%nOrbAtom)
      iSp = species(iAt)
      do iSh = 1, orb%nShell(iSp)
        deltaQUniqU(iHubbU(iSh, iSp), iAt) = &
            & deltaQUniqU(iHubbU(iSh, iSp), iAt) &
            &+ deltaQPerLShell(iSh, iAt)
      end do
    end do

  end subroutine getNetChargesPerUniqU


end module charges
