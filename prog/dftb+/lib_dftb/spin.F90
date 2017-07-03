!--------------------------------------------------------------------------------------------------!
!  DFTB+: general package for performing fast atomistic simulations                                !
!  Copyright (C) 2017  DFTB+ developers group                                                      !
!                                                                                                  !
!  See the LICENSE file for terms of usage and distribution.                                       !
!--------------------------------------------------------------------------------------------------!

#:include 'common.fypp'

!!* Module containing various routines for spin polarised calculations
!!* Intended to be used with SCC switched on !
module spin
  use assert
  use accuracy
  use message
  use commontypes
  use shift
  implicit none
  private

  public :: getEnergySpin, addSpinShift
  public :: Spin_getOrbitalEquiv, ud2qm, qm2ud

  interface getEnergySpin
    module procedure getEnergySpin_total
    module procedure getEnergySpin_atom
  end interface

  !!* swap from up/down to charge/magnetisation
  interface ud2qm
    module procedure ud2qm2
    module procedure ud2qm3
    module procedure ud2qm4
  end interface

  !!* swap from charge/magnetisation to up/down
  interface qm2ud
    module procedure qm2ud2
    module procedure qm2ud3
    module procedure qm2ud4
  end interface

contains

  !!* Constructs the spin-polarised shell shift from
  !!* shift_l = sum_l' W_ll' p_l'
  !!* @param shift resulting shell-shifts for the system
  !!* @param chargePerShell spin resolved charges for each shell
  !!* @param species Species of each atom
  !!* @param orb  Information about the orbitals and their angular momenta
  !!* @param spinW Spin coupling constants.
  !!* @todo Add more asserts
  subroutine addSpinShift(shift,chargePerShell,species,orb,spinW)
    real(dp), intent(inout)     :: shift(:,:,0:)
    real(dp), intent(in)        :: chargePerShell(:,:,0:)
    integer, intent(in)         :: species(:)
    type(TOrbitals), intent(in) :: orb
    real(dp), intent(in)        :: spinW(:,:,:)

    integer  :: nAtom, iAtom, iSpecies, iShell, iShell2, nSpin, iSpin

    nAtom = size(chargePerShell,dim=2)
    @:ASSERT(nAtom > 0)
    @:ASSERT(size(shift,dim=2)==nAtom)
    @:ASSERT(all(shape(chargePerShell)==shape(shift)))
    nSpin = size(chargePerShell,dim=3) - 1 ! counts from 0
    @:ASSERT(nSpin == 1 .or. nSpin == 3)

    do iSpin = 1, nSpin
      do iAtom = 1, nAtom
        iSpecies = species(iAtom)
        do iShell = 1, orb%nShell(iSpecies)
          do iShell2 = 1, orb%nShell(iSpecies)
            shift(iShell,iAtom,iSpin) =  shift(iShell,iAtom,iSpin) + &
                & spinW(iShell,iShell2,iSpecies) * &
                & chargePerShell(iShell2,iAtom,iSpin)
          end do
        end do
      end do
    end do

  end subroutine addSpinShift

  !!* Returns the total energy contribution of the spin polarisation
  !!* @param rslt Contains the atomic contributions on exit
  !!* @param chargePerShell spin resolved charges for each shell
  !!* @param shiftPerShell spin shift for each shell
  subroutine getEnergySpin_total(rslt, chargePerShell, shiftPerShell)
    real(dp), intent(out) :: rslt
    real(dp), intent(in)  :: chargePerShell(:,:,:)
    real(dp), intent(in)  :: shiftPerShell(:,:,:)

    @:ASSERT(all(shape(chargePerShell)==shape(shiftPerShell)))
    @:ASSERT(size(chargePerShell,dim=3)>1 .and. size(chargePerShell,dim=3)<5)

    ! safe as the shift for the spin=0 component is 0 at the moment
    rslt = sum(chargePerShell(:,:,:)*shiftPerShell(:,:,:))

  end subroutine getEnergySpin_total

  !!* Atom resolved part of the spin energy
  !!* @param rslt Contains the atomic contributions on exit
  !!* @param chargePerShell spin resolved charges for each shell
  !!* @param shiftPerShell spin shift for each shell
  subroutine getEnergySpin_atom(rslt, chargePerShell, shiftPerShell)
    real(dp), intent(out) :: rslt(:)
    real(dp), intent(in)  :: chargePerShell(:,:,:)
    real(dp), intent(in)  :: shiftPerShell(:,:,:)

    @:ASSERT(size(rslt)==size(chargePerShell,dim=2))
    @:ASSERT(all(shape(chargePerShell)==shape(shiftPerShell)))
    @:ASSERT(size(chargePerShell,dim=3)>1 .and. size(chargePerShell,dim=3)<5)

    ! safe as the shift for the spin=0 component is 0 at the moment
    rslt(:) = sum(sum(chargePerShell(:,:,:)*shiftPerShell(:,:,:),dim=3),dim=1)

  end subroutine getEnergySpin_atom

  !!* Returns the equivalence between the orbitals in the spin interaction.
  !!* @param orb  Information about the orbitals and their angular momenta
  !!* @param species Species of each atom
  !!* @param equiv The equivalence vector on return.
  !!* @note The current version assumes, that no shells only the orbitals
  !!* inside each shell are equivalent, which is in most cases true anyway.
  !!* @todo Proper analysis of the spin coupling constants to watch for
  !!*   eventual equivalence.
  subroutine Spin_getOrbitalEquiv(orb, species, equiv)
    type(TOrbitals), intent(in) :: orb
    integer, intent(in) :: species(:)
    integer, intent(out) :: equiv(:,:,:)

    integer :: nAtom, nSpin
    integer :: iAt, iOrb, iS, ind, iSp

    nAtom = size(equiv, dim=2)
    nSpin = size(equiv, dim=3)

    @:ASSERT(size(equiv, dim=1) == orb%mOrb)
    @:ASSERT(nSpin == 1 .or. nSpin == 2 .or. nSpin == 4)
    @:ASSERT(nAtom > 0)

    equiv(:,:,:) = 0
    ind = 1
    do iAt = 1, nAtom
      iSp = species(iAt)
      do iOrb = 1, orb%nOrbSpecies(iSp)
        equiv(iOrb, iAt, 1) = ind + orb%iShellOrb(iOrb, iSp) - 1
      end do
      ind = ind + orb%nShell(iSp)
    end do
    do iS = 2, nSpin
      ind = maxval(equiv)
      where (equiv(:,:,1) /= 0)
        equiv(:,:,iS) = equiv(:,:,1) + ind
      end where
    end do

  end subroutine Spin_getOrbitalEquiv

  !!* converts a charge/magnetization set into a up/down
  !!* @param x array of data, last index spin
  subroutine qm2ud2(x)
    real(dp), intent(inout) :: x(:,:)

    integer :: nSpin, nElements, ii

    nElements = size(x,dim=1)
    nSpin = size(x,dim=2)
    @:ASSERT( nSpin == 1 .or. nSpin == 2 .or. nSpin == 4 )

    select case(nSpin)
    case (1)
      ! nothing to do
    case (2)
      do ii = 1, nElements
        x(ii,1) = 0.5_dp*(x(ii,1) + x(ii,2))
        x(ii,2) = x(ii,1) - x(ii,2)
      end do
    case (4)
      ! nothing to do
    end select

  end subroutine qm2ud2

  !!* converts a charge/magnetization set into a up/down
  !!* @param x array of data, last index spin
  subroutine qm2ud3(x)
    real(dp), intent(inout) :: x(:,:,:)

    integer :: nSpin, ii, jj
    integer :: nElements(2)

    nElements(1) = size(x,dim=1)
    nElements(2) = size(x,dim=2)
    nSpin = size(x,dim=3)
    @:ASSERT( nSpin == 1 .or. nSpin == 2 .or. nSpin == 4)

    select case(nSpin)
    case (1)
      ! nothing to do
    case (2)
      do jj = 1, nElements(2)
        do ii = 1, nElements(1)
          x(ii,jj,1) = 0.5_dp*(x(ii,jj,1) + x(ii,jj,2))
          x(ii,jj,2) = x(ii,jj,1) - x(ii,jj,2)
        end do
      end do
    case (4)
      ! nothing to do
    end select

  end subroutine qm2ud3

  !!* converts a charge/magnetization set into a up/down
  !!* @param x array of data, last index spin
  subroutine qm2ud4(x)
    real(dp), intent(inout) :: x(:,:,:,:)

    integer :: nSpin, ii, jj, kk
    integer :: nElements(3)

    nElements(1) = size(x,dim=1)
    nElements(2) = size(x,dim=2)
    nElements(3) = size(x,dim=3)
    nSpin = size(x,dim=4)
    @:ASSERT( nSpin == 1 .or. nSpin == 2 .or. nSpin == 4)

    select case(nSpin)
    case (1)
      ! nothing to do
    case(2)
      do kk = 1, nElements(3)
        do jj = 1, nElements(2)
          do ii = 1, nElements(1)
            x(ii,jj,kk,1) = 0.5_dp*(x(ii,jj,kk,1) + x(ii,jj,kk,2))
            x(ii,jj,kk,2) = x(ii,jj,kk,1) - x(ii,jj,kk,2)
          end do
        end do
      end do
    case (4)
      ! nothing to do
    end select

  end subroutine qm2ud4

  !!* converts a up/down set into a charge/magnetization
  !!* @param x array of data, last index spin
  subroutine ud2qm2(x)
    real(dp), intent(inout) :: x(:,:)

    integer :: nSpin, nElements, ii

    nElements = size(x,dim=1)
    nSpin = size(x,dim=2)
    @:ASSERT( nSpin == 1 .or. nSpin == 2 .or. nSpin == 4)

    select case(nSpin)
    case (1)
      ! nothing to do
    case (2)
      do ii = 1, nElements
        x(ii,1) = x(ii,1) + x(ii,2)
        x(ii,2) = x(ii,1) - 2.0_dp * x(ii,2)
      end do
    case (4)
      ! nothing to do
    end select

  end subroutine ud2qm2

  !!* converts a up/down set into a charge/magnetization
  !!* @param x array of data, last index spin
  subroutine ud2qm3(x)
    real(dp), intent(inout) :: x(:,:,:)

    integer :: nSpin, ii, jj
    integer :: nElements(2)

    nElements(1) = size(x,dim=1)
    nElements(2) = size(x,dim=2)
    nSpin = size(x,dim=3)
    @:ASSERT( nSpin == 1 .or. nSpin == 2 .or. nSpin == 4 )

    select case(nSpin)
    case (1)
      ! nothing to do
    case (2)
      do jj = 1, nElements(2)
        do ii = 1, nElements(1)
          x(ii,jj,1) = x(ii,jj,1) + x(ii,jj,2)
          x(ii,jj,2) = x(ii,jj,1) - 2.0_dp * x(ii,jj,2)
        end do
      end do
    case (4)
      ! nothing to do
    end select

  end subroutine ud2qm3

  !!* converts a charge/magnetization set into a up/down
  !!* @param x array of data, last index spin
  subroutine ud2qm4(x)
    real(dp), intent(inout) :: x(:,:,:,:)

    integer :: nSpin, ii, jj, kk
    integer :: nElements(3)

    nElements(1) = size(x,dim=1)
    nElements(2) = size(x,dim=2)
    nElements(3) = size(x,dim=3)
    nSpin = size(x,dim=4)
    @:ASSERT( nSpin == 1 .or. nSpin == 2 .or. nSpin == 4)

    select case(nSpin)
    case (1)
      ! nothing to do
    case (2)
      do kk = 1, nElements(3)
        do jj = 1, nElements(2)
          do ii = 1, nElements(1)
            x(ii,jj,kk,1) = x(ii,jj,kk,1) + x(ii,jj,kk,2)
            x(ii,jj,kk,2) = x(ii,jj,kk,1) - 2.0_dp * x(ii,jj,kk,2)
          end do
        end do
      end do
    case (4)
      ! nothing to do
    end select

  end subroutine ud2qm4

end module spin
