!--------------------------------------------------------------------------------------------------!
!  DFTB+: general package for performing fast atomistic simulations                                !
!  Copyright (C) 2018  DFTB+ developers group                                                      !
!                                                                                                  !
!  See the LICENSE file for terms of usage and distribution.                                       !
!--------------------------------------------------------------------------------------------------!

#:include 'common.fypp'

!> Module to wrap around the different shift contributions in the DFTB energy expressions
module potentials
  use assert
  use accuracy, only : dp
  use commontypes
  implicit none

  public :: TPotentials, init

  private


  !> data type to store components of the internal and external potential as named variables - makes
  !> extending expressions easier.
  !>
  !> Note: the reason for splitting internal and external potentials is that internals require 0.5
  !> scaling for double counting when evaluating total energies of the form V * n
  type TPotentials
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

    !> pSIC/DFTB+U etc. potential
    real(dp), allocatable :: orbitalBlock(:,:,:,:)

    !> L.S etc where these are imaginary coefficients
    real(dp), allocatable :: iorbitalBlock(:,:,:,:)
  end type TPotentials


  !> Initialise the structure
  interface init
    module procedure Potentials_init
  end interface

contains


  !> Allocates storage for the potential components
  subroutine Potentials_init(self,orb,nAtom,nSpin)

    !> data structure to allocate
    type(TPotentials), intent(out) :: self

    !> information about the orbitals and their angular momenta
    type(TOrbitals), intent(in) :: orb

    !> number of atoms needed for atom resolved arrays
    integer, intent(in) :: nAtom

    !> number of spins
    integer, intent(in) :: nSpin

    @:ASSERT(.not. self%tInitialised)
    @:ASSERT(nSpin == 1 .or. nSpin == 2 .or. nSpin == 4)
    @:ASSERT(nAtom > 0)
    @:ASSERT(orb%mShell > 0)
    @:ASSERT(orb%mOrb > 0)

    allocate(self%intAtom(nAtom,nSpin))
    allocate(self%intShell(orb%mShell,nAtom,nSpin))
    allocate(self%intBlock(orb%mOrb,orb%mOrb,nAtom,nSpin))
    allocate(self%extAtom(nAtom,nSpin))
    allocate(self%extShell(orb%mShell,nAtom,nSpin))
    allocate(self%extBlock(orb%mOrb,orb%mOrb,nAtom,nSpin))
    allocate(self%orbitalBlock(orb%mOrb,orb%mOrb,nAtom,nSpin))
    allocate(self%iorbitalBlock(orb%mOrb,orb%mOrb,nAtom,nSpin))
    self%intAtom = 0.0_dp
    self%intShell = 0.0_dp
    self%intBlock = 0.0_dp
    self%extAtom = 0.0_dp
    self%extShell = 0.0_dp
    self%extBlock = 0.0_dp
    self%orbitalBlock = 0.0_dp
    self%iorbitalBlock = 0.0_dp
    self%tInitialised = .true.

  end subroutine Potentials_init

end module potentials
