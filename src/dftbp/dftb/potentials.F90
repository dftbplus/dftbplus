!--------------------------------------------------------------------------------------------------!
!  DFTB+: general package for performing fast atomistic simulations                                !
!  Copyright (C) 2006 - 2022  DFTB+ developers group                                               !
!                                                                                                  !
!  See the LICENSE file for terms of usage and distribution.                                       !
!--------------------------------------------------------------------------------------------------!

#:include 'common.fypp'

!> Module to wrap around the different shift contributions in the DFTB energy expressions
module dftbp_dftb_potentials
  use dftbp_common_accuracy, only : dp
  use dftbp_type_commontypes, only : TOrbitals
  implicit none

  private
  public :: TPotentials, TPotentials_init, TAtomExtPotInput


  !> data type to store components of the internal and external potential as named variables - makes
  !> extending expressions easier.
  !>
  !> Note: the reason for splitting internal and external potentials is that internals require 0.5
  !> scaling for double counting when evaluating total energies of the form V * n
  type TPotentials

    !> Is data structure initialised
    logical :: tInitialised = .false.

    !> internal atom and spin resolved potential
    real(dp), allocatable :: intAtom(:,:)

    !> internal shell and spin resolved potential
    real(dp), allocatable :: intShell(:,:,:)

    !> internal block and spin resolved potential
    real(dp), allocatable :: intBlock(:,:,:,:)

    !> external atom and spin resolved potential
    real(dp), allocatable :: extAtom(:,:)

    !> external shell and spin resolved potential
    real(dp), allocatable :: extShell(:,:,:)

    !> external block and spin resolved potential
    real(dp), allocatable :: extBlock(:,:,:,:)

    !> gradient of the external potential with respect of nucleus coordinates
    real(dp), allocatable :: extGrad(:,:)

    !> pSIC/DFTB+U etc. potential
    real(dp), allocatable :: orbitalBlock(:,:,:,:)

    !> L.S etc where these are imaginary coefficients
    real(dp), allocatable :: iorbitalBlock(:,:,:,:)

    !> If performing a contact calculation, variable for retaining the shell resolved electrostatics
    !> for later storage. Ony the charge related potential is stored, so last index will be
    !> allocated as 1 in most cases
    real(dp), allocatable :: coulombShell(:,:,:)

    !> internal atom and spin resolved onsite only potential (i.e. relating to net, not gross
    !> populations)
    real(dp), allocatable :: intOnSiteAtom(:,:)

    !> external atom and spin resolved onsite only potential (i.e. relating to net, not gross
    !> populations)
    real(dp), allocatable :: extOnSiteAtom(:,:)

    !> Atomic dipolar contribution to the Hamiltonian
    real(dp), allocatable :: dipoleAtom(:,:)

    !> Atomic quadrupolar contribution to the Hamiltonian
    real(dp), allocatable :: quadrupoleAtom(:,:)

    !> External dipolar contribution to the Hamiltonian
    real(dp), allocatable :: extDipoleAtom(:,:)

  end type TPotentials


  !> External potentials at atomic sites
  type TAtomExtPotInput

    !> Potentials experienced by gross population charge on atom
    real(dp), allocatable :: Vext(:)

    !> Atom(s) at which potential is present
    integer, allocatable :: iAt(:)

    !> Potentials experienced by net population charge on atom
    real(dp), allocatable :: VextOnSite(:)

    !> Atom(s) at which on-site potential is present
    integer, allocatable :: iAtOnSite(:)

  end type TAtomExtPotInput

contains


  !> Allocates storage for the potential components
  subroutine TPotentials_init(this, orb, nAtom, nSpin, nDipole, nQuadrupole, extAtPotentials)

    !> data structure to allocate
    type(TPotentials), intent(out) :: this

    !> information about the orbitals and their angular momenta
    type(TOrbitals), intent(in) :: orb

    !> number of atoms needed for atom resolved arrays
    integer, intent(in) :: nAtom

    !> number of spins
    integer, intent(in) :: nSpin

    !> Number of dipole moment components
    integer, intent(in) :: nDipole

    !> Number of quadrupole moment components
    integer, intent(in) :: nQuadrupole

    !> Should the on site potentials be allocated?
    type(TAtomExtPotInput), intent(in), optional :: extAtPotentials

    @:ASSERT(.not. this%tInitialised)
    @:ASSERT(nSpin == 1 .or. nSpin == 2 .or. nSpin == 4)
    @:ASSERT(nAtom > 0)
    @:ASSERT(nDipole >= 0)
    @:ASSERT(nQuadrupole >= 0)
    @:ASSERT(orb%mShell > 0)
    @:ASSERT(orb%mOrb > 0)

    allocate(this%intAtom(nAtom,nSpin))
    allocate(this%intShell(orb%mShell,nAtom,nSpin))
    allocate(this%intBlock(orb%mOrb,orb%mOrb,nAtom,nSpin))
    allocate(this%extAtom(nAtom,nSpin))
    allocate(this%extShell(orb%mShell,nAtom,nSpin))
    allocate(this%extBlock(orb%mOrb,orb%mOrb,nAtom,nSpin))
    allocate(this%extGrad(3, nAtom))
    allocate(this%orbitalBlock(orb%mOrb,orb%mOrb,nAtom,nSpin))
    allocate(this%iorbitalBlock(orb%mOrb,orb%mOrb,nAtom,nSpin))
    this%intAtom = 0.0_dp
    this%intShell = 0.0_dp
    this%intBlock = 0.0_dp
    this%extAtom = 0.0_dp
    this%extShell = 0.0_dp
    this%extBlock = 0.0_dp
    this%extGrad(:,:) = 0.0_dp
    this%orbitalBlock = 0.0_dp
    this%iorbitalBlock = 0.0_dp
    if (nDipole > 0) then
      allocate(this%dipoleAtom(nDipole, nAtom))
      this%dipoleAtom(:,:) = 0.0_dp
      allocate(this%extDipoleAtom(nDipole, nAtom))
      this%extDipoleAtom(:,:) = 0.0_dp
    end if
    if (nQuadrupole > 0) then
      allocate(this%quadrupoleAtom(nQuadrupole, nAtom))
      this%quadrupoleAtom(:,:) = 0.0_dp
    end if

    if (present(extAtPotentials)) then

      if (allocated(extAtPotentials%VextOnSite)) then

        allocate(this%extOnSiteAtom(nAtom,nSpin))
        this%extOnSiteAtom(:,:) = 0.0_dp

        allocate(this%intOnSiteAtom(nAtom,nSpin))
        this%intOnSiteAtom(:,:) = 0.0_dp

      end if

    end if

    this%tInitialised = .true.

  end subroutine TPotentials_init

end module dftbp_dftb_potentials
