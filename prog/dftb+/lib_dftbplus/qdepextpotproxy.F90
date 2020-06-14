!--------------------------------------------------------------------------------------------------!
!  DFTB+: general package for performing fast atomistic simulations                                !
!  Copyright (C) 2006 - 2020  DFTB+ developers group                                               !
!                                                                                                  !
!  See the LICENSE file for terms of usage and distribution.                                       !
!--------------------------------------------------------------------------------------------------!

!> Contains a proxy communicating with external generators of population dependent potentials.
module dftbp_qdepextpotproxy
  use dftbp_accuracy, only : dp
  use dftbp_qdepextpotgen, only : TQDepExtPotGen, TQDepExtPotGenWrapper
  use dftbp_commontypes, only : TOrbitals
  use dftbp_shift, only : total_shift
  implicit none
  private


  public :: TQDepExtPotProxy, TQDepExtPotProxy_init


  !> Collection of external q-dependent potentials queried during the SCC-cycle
  type :: TQDepExtPotProxy
    private
    !> collection of external potentials
    type(TQDepExtPotGenWrapper), allocatable :: generators(:)
    !> energy contributions to DFTB atoms due to potentials
    real(dp), allocatable :: energyAtom(:)
  contains
    !> add potential contribution
    procedure :: addPotential => TQDepExtPotProxy_addPotential
    !> add energy contribution
    procedure :: addEnergy => TQDepExtPotProxy_addEnergy
    !> add force contribution
    procedure :: addGradientDc => TQDepExtPotProxy_addGradientDc
  end type TQDepExtPotProxy


contains

  !> Initializes proxy for querying a collection of q-dependent external potential generators.
  subroutine TQDepExtPotProxy_init(this, extPotGenerators)

    !> Instance.
    type(TQDepExtPotProxy), intent(out) :: this

    !> External potential generators to consider.
    type(TQDepExtPotGenWrapper), intent(in) :: extPotGenerators(:)

    this%generators = extPotGenerators

  end subroutine TQDepExtPotProxy_init


  !> Adds the external potential from the external potential generators.
  subroutine TQDepExtPotProxy_addPotential(this, deltaQAtom, deltaQShell, orb, species, potential)

    !> Instance.
    class(TQDepExtPotProxy), intent(inout) :: this

    !> Spin unpolarised net population per atom.
    real(dp), intent(in) :: deltaQAtom(:)

    !> Spin unpolarised net population per shell.
    real(dp), intent(in) :: deltaQShell(:,:)

    !> Orbital information
    type(TOrbitals), intent(in) :: orb

    !> Species
    integer, intent(in) :: species(:)

    !> Shell resolved potential to update.
    real(dp), intent(inout) :: potential(:,:,:,:)

    real(dp), allocatable :: potAtom(:,:), potShell(:,:,:), potAtomTmp(:), potShellTmp(:,:)
    integer :: mShell, nAtom, nSpin
    integer :: iGen

    mShell = size(deltaQShell, dim=1)
    nAtom = size(potential, dim=3)
    nSpin = size(potential, dim=4)
    allocate(potAtom(nAtom, nSpin))
    allocate(potShell(mShell, nAtom, nSpin))
    allocate(potAtomTmp(nAtom))
    allocate(potShellTmp(mShell, nAtom))
    potAtom(:,:) = 0.0_dp
    potShell(:,:,:) = 0.0_dp
    if (.not. allocated(this%energyAtom)) then
      allocate(this%energyAtom(nAtom))
    end if
    this%energyAtom(:) = 0.0_dp
    do iGen = 1, size(this%generators)
      call this%generators(iGen)%instance%getExternalPot(deltaQAtom, deltaQShell, potAtomTmp,&
          & potShellTmp)
      potAtom(:,1) = potAtom(:,1) + potAtomTmp
      potShell(:,:,1) = potShell(:,:,1) + potShellTmp
      this%energyAtom(:) = this%energyAtom + deltaQAtom * potAtomTmp
      this%energyAtom(:) = this%energyAtom + sum(deltaQShell * potShellTmp, dim=1)
    end do
    call total_shift(potShell, potAtom, orb, species)
    call total_shift(potential, potShell, orb, species)

  end subroutine TQDepExtPotProxy_addPotential


  !> Adds the energy contribution of the q-dependent external potentials.
  !>
  !> Note: This should only called after the potential had been queried via the
  !> addPotential() procedure.
  !>
  subroutine TQDepExtPotProxy_addEnergy(this, energies)

    !> Instance.
    class(TQDepExtPotProxy), intent(inout) :: this

    !> Energy per atoms, to which the energies from the external potentials should be added to.
    real(dp), intent(inout) :: energies(:)

    energies(:) = energies + this%energyAtom

  end subroutine TQDepExtPotProxy_addEnergy


  !> Adds the "double counting" gradient contribution of the q-dependent external potentials.
  !>
  !> The double counting part of the gradients is the one, which is not obtained by the derivatives
  !> of the shift-vectors.
  !>
  subroutine TQDepExtPotProxy_addGradientDc(this, deltaQAtom, deltaQShell, gradients)

    !> Instance.
    class(TQDepExtPotProxy), intent(inout) :: this

    !> Net population per atom.
    real(dp), intent(in) :: deltaQAtom(:)

    !> Net population per shell.
    real(dp), intent(in) :: deltaQShell(:,:)

    !> Gradients to upgrade.
    real(dp), intent(inout) :: gradients(:,:)

    real(dp), allocatable :: extPotGrad(:,:), deltaQAtomSpread(:,:)
    integer :: iGen

    allocate(extPotGrad(size(gradients, dim=1), size(gradients, dim=2)))
    deltaQAtomSpread = spread(deltaQAtom, 1, 3)
    do iGen = 1, size(this%generators)
      call this%generators(iGen)%instance%getExternalPotGrad(deltaQAtom, deltaQShell, extPotGrad)
      gradients(:,:) = gradients + deltaQAtomSpread * extPotGrad
    end do

  end subroutine TQDepExtPotProxy_addGradientDc


end module dftbp_qdepextpotproxy
