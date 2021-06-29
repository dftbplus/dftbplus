!--------------------------------------------------------------------------------------------------!
!  DFTB+: general package for performing fast atomistic simulations                                !
!  Copyright (C) 2006 - 2021  DFTB+ developers group                                               !
!                                                                                                  !
!  See the LICENSE file for terms of usage and distribution.                                       !
!--------------------------------------------------------------------------------------------------!

#:include 'common.fypp'

!> Contains unique Hubbard U values for a given system
module dftbp_dftb_uniquehubbard
  use dftbp_common_accuracy, only : dp, minHubDiff
  use dftbp_type_commontypes, only : TOrbitals
  implicit none

  private
  public :: TUniqueHubbard, TUniqueHubbard_init


  !> Contains information about species (eventually contracted unique) Hubbard U values
  type :: TUniqueHubbard

    !> Nr. of uniq Us per species. Shape: [nSpecies]
    integer, allocatable :: nHubbU(:)

    !> Uniq Us per species. Shape: [mShell, nSpecies]
    real(dp), allocatable :: uniqHubbU(:,:)

    !> Mapping L-shell to unique Hubbard U. Shape: [mShell, nSpecies]
    integer, allocatable :: iHubbU(:, :)

    !> Maximal nr. of Hubbard Us per species
    integer :: mHubbU

  contains

    procedure :: sumOverUniqueU
    procedure :: getOrbitalEquiv

  end type TUniqueHubbard


contains


  !> Initialises a TUniqueHubbard instance.
  subroutine TUniqueHubbard_init(this, hubbU, orb)

    !> Instance.
    type(TUniqueHubbard), intent(out) :: this

    !> Hubbard U value for each species. Shape: [mShell, nSpecies].
    real(dp), intent(in) :: hubbU(:,:)

    !> Orbital information
    type(TOrbitals), intent(in) :: orb

    integer :: nSpecies, mShell
    integer :: iSp, iL, iU

    nSpecies = size(orb%nOrbSpecies)
    mShell = orb%mShell

    allocate(this%nHubbU(nSpecies))
    allocate(this%uniqHubbU(mShell, nSpecies))
    allocate(this%iHubbU(mShell, nSpecies))

    this%iHubbU(:, :) = 0
    this%iHubbU(1, :) = 1
    this%nHubbU(:) = 1
    this%uniqHubbU(:, :) = 0.0_dp
    this%uniqHubbU(1, :) = hubbU(1, :)
    do iSp = 1, nSpecies
      lpShell: do iL = 2, orb%nShell(iSp)
        do iU = 1, this%nHubbU(iSp)
          if (abs(hubbU(iL, iSp) - this%uniqHubbU(iU, iSp)) < minHubDiff) then
            this%iHubbU(iL, iSp) = iU
            cycle lpShell
          end if
        end do
        @:ASSERT(this%iHubbU(iL, iSp) == 0)
        this%nHubbU(iSp) = this%nHubbU(iSp) + 1
        this%uniqHubbU(this%nHubbU(iSp), iSp) = hubbU(iL, iSp)
        this%iHubbU(iL, iSp) = this%nHubbU(iSp)
      end do lpShell
    end do

    this%mHubbU = maxval(this%nHubbU)

  end subroutine TUniqueHubbard_init


  !> Performs a sum over unqiue U values
  subroutine sumOverUniqueU(this, valueShell, species, orb, valueUniqueU)

    !> Instance
    class(TUniqueHubbard), intent(in) :: this

    !> Value to sum up. Shape [maxShellPerSpecies, nAtom].
    real(dp), intent(in) :: valueShell(:,:)

    !> Species index of each atom. Shape: [nAtom]
    integer, intent(in) :: species(:)

    !> Basis information
    type(TOrbitals), intent(in) :: orb

    !> Values summed up over unique Hubbard U values. Shape: [maxUniqueUPerSpecies, nAtom].
    real(dp), intent(out) :: valueUniqueU(:,:)

    integer :: iAt, iSp, iSh

    valueUniqueU(:,:) = 0.0_dp
    do iAt = 1, size(orb%nOrbAtom)
      iSp = species(iAt)
      do iSh = 1, orb%nShell(iSp)
        valueUniqueU(this%iHubbU(iSh, iSp), iAt) = valueUniqueU(this%iHubbU(iSh, iSp), iAt)&
            & + valueShell(iSh, iAt)
      end do
    end do

  end subroutine sumOverUniqueU


  !> Return an orbital equivalency relation with the given unique U values
  subroutine getOrbitalEquiv(this, orb, species, equiv)

    !> Resulting module variables
    class(TUniqueHubbard), intent(in) :: this

    !> Contains information about the atomic orbitals in the system
    type(TOrbitals), intent(in) :: orb

    !> Type of each atom (nAtom).
    integer, intent(in) :: species(:)

    !> The vector describing the equivalence on return.
    integer, intent(out) :: equiv(:,:,:)

    integer :: nAtom, nSpin, nUniqueOrb
    integer :: iAt, iSp, iOrb, iS

    nAtom = size(equiv, dim=2)
    nSpin = size(equiv, dim=3)

    equiv(:,:,1) = 0
    nUniqueOrb = 0
    do iAt = 1, nAtom
      iSp = species(iAt)
      do iOrb = 1, orb%nOrbSpecies(iSp)
        equiv(iOrb, iAt, 1) = this%iHubbU(orb%iShellOrb(iOrb, iSp), iSp) + nUniqueOrb
      end do
      nUniqueOrb = nUniqueOrb + maxval(this%iHubbU(:, iSp))
    end do

    do iS = 2, size(equiv, dim=3)
      equiv(:,:,iS) = equiv(:,:,1)
    end do

  end subroutine getOrbitalEquiv


end module dftbp_dftb_uniquehubbard
