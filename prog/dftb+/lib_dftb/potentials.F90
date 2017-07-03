!--------------------------------------------------------------------------------------------------!
!  DFTB+: general package for performing fast atomistic simulations                                !
!  Copyright (C) 2017  DFTB+ developers group                                                      !
!                                                                                                  !
!  See the LICENSE file for terms of usage and distribution.                                       !
!--------------------------------------------------------------------------------------------------!

#:include 'common.fypp'

!!* Module to wrap around the different shift contributions in the DFTB
!!* energy expressions
module potentials
  use assert
  use accuracy, only : dp
  use commontypes
  implicit none

  public :: TPotentials, init

  private

  !!* data type to store components of the internal and external potential as
  !!* named variables - makes extending expressions easier.
  !!* @note the reason for splitting internal and external potentials is that
  !!* internals require 0.5 scaling for double counting when evaluating total
  !!* energies
  type TPotentials
    logical :: tInitialised = .false.
    real(dp), allocatable :: intAtom(:,:) ! internal atom and spin resolved pot.
    real(dp), allocatable :: intShell(:,:,:)
    real(dp), allocatable :: intBlock(:,:,:,:)
    real(dp), allocatable :: extAtom(:,:) ! external atom and spin resolved pot.
    real(dp), allocatable :: extShell(:,:,:)
    real(dp), allocatable :: extBlock(:,:,:,:)
    real(dp), allocatable :: orbitalBlock(:,:,:,:) ! pSIC/DFTB+U etc.
    ! L.S etc where these are imaginary coefficients
    real(dp), allocatable :: iorbitalBlock(:,:,:,:)
  end type TPotentials

  interface init
    module procedure Potentials_init
  end interface


contains

  !!* Allocates storage for the potential components
  !!* @param self data structure to allocate
  !!* @param orb information about the orbitals and their angular momenta
  !!* @param nAtom number of atoms needed for atom resolved arrays
  !!* @param nSpin number of spins
  subroutine Potentials_init(self,orb,nAtom,nSpin)
    type(TPotentials), intent(out) :: self
    type(TOrbitals), intent(in) :: orb
    integer, intent(in) :: nAtom
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
