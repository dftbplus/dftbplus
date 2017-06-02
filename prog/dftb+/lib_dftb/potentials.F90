!--------------------------------------------------------------------------------------------------!
!  DFTB+: general package for performing fast atomistic simulations                                !
!  Copyright (C) 2017  DFTB+ developers group                                                      !
!                                                                                                  !
!  See the LICENSE file for terms of usage and distribution.                                       !
!--------------------------------------------------------------------------------------------------!

!!* Module to wrap around the different shift contributions in the DFTB 
!!* energy expressions
module potentials
#include "allocate.h"
#include "assert.h"
  use accuracy, only : dp
  use commontypes
  implicit none
  
  public :: TPotentials, create, destroy
  
  private
  
  !!* data type to store components of the internal and external potential as
  !!* named variables - makes extending expressions easier.
  !!* @note the reason for splitting internal and external potentials is that
  !!* internals require 0.5 scaling for double counting when evaluating total
  !!* energies
  type TPotentials
    logical :: tInitialised = .false.
    real(dp), pointer :: intAtom(:,:) ! internal atom and spin resolved pot.
    real(dp), pointer :: intShell(:,:,:)
    real(dp), pointer :: intBlock(:,:,:,:)
    real(dp), pointer :: extAtom(:,:) ! external atom and spin resolved pot.
    real(dp), pointer :: extShell(:,:,:)
    real(dp), pointer :: extBlock(:,:,:,:)
    real(dp), pointer :: orbitalBlock(:,:,:,:) ! pSIC/DFTB+U etc.
    real(dp), pointer :: iorbitalBlock(:,:,:,:) ! L.S etc where these are
                                                ! imaginary coefficients
  end type TPotentials
  
  interface create
    module procedure Potentials_create
  end interface
  
  interface destroy
    module procedure Potentials_destroy
  end interface
  
contains
  
  !!* Allocates storage for the potential components
  !!* @param self data structure to allocate
  !!* @param orb information about the orbitals and their angular momenta
  !!* @param nAtom number of atoms needed for atom resolved arrays
  !!* @param nSpin number of spins
  subroutine Potentials_create(self,orb,nAtom,nSpin)
    type(TPotentials), intent(out) :: self
    type(TOrbitals), pointer :: orb
    integer, intent(in) :: nAtom
    integer, intent(in) :: nSpin
    
    ASSERT(.not. self%tInitialised)
    ASSERT(nSpin == 1 .or. nSpin == 2 .or. nSpin == 4)
    ASSERT(nAtom > 0)
    ASSERT(orb%mShell > 0)
    ASSERT(orb%mOrb > 0)
    
    INITALLOCATE_PARR(self%intAtom,(nAtom,nSpin))
    INITALLOCATE_PARR(self%intShell,(orb%mShell,nAtom,nSpin))
    INITALLOCATE_PARR(self%intBlock,(orb%mOrb,orb%mOrb,nAtom,nSpin))
    INITALLOCATE_PARR(self%extAtom,(nAtom,nSpin))
    INITALLOCATE_PARR(self%extShell,(orb%mShell,nAtom,nSpin))
    INITALLOCATE_PARR(self%extBlock,(orb%mOrb,orb%mOrb,nAtom,nSpin))
    INITALLOCATE_PARR(self%orbitalBlock,(orb%mOrb,orb%mOrb,nAtom,nSpin))
    INITALLOCATE_PARR(self%iorbitalBlock,(orb%mOrb,orb%mOrb,nAtom,nSpin))
    self%intAtom = 0.0_dp
    self%intShell = 0.0_dp
    self%intBlock = 0.0_dp
    self%extAtom = 0.0_dp
    self%extShell = 0.0_dp
    self%extBlock = 0.0_dp    
    self%orbitalBlock = 0.0_dp
    self%iorbitalBlock = 0.0_dp
    self%tInitialised = .true.
    
  end subroutine Potentials_create
  
  !!* De-allocates storage for the potential components
  subroutine Potentials_destroy(self)
    type(TPotentials), intent(inout) :: self
    
    ASSERT(self%tInitialised)
    
    DEALLOCATE_PARR(self%intAtom)
    DEALLOCATE_PARR(self%intShell)
    DEALLOCATE_PARR(self%intBlock)
    DEALLOCATE_PARR(self%extAtom)
    DEALLOCATE_PARR(self%extShell)    
    DEALLOCATE_PARR(self%extBlock)
    DEALLOCATE_PARR(self%orbitalBlock)
    DEALLOCATE_PARR(self%iorbitalBlock)

    self%tInitialised = .false.
    
  end subroutine Potentials_destroy
  
end module potentials
