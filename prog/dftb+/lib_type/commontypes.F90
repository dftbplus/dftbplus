!--------------------------------------------------------------------------------------------------!
!  DFTB+: general package for performing fast atomistic simulations                                !
!  Copyright (C) 2017  DFTB+ developers group                                                      !
!                                                                                                  !
!  See the LICENSE file for terms of usage and distribution.                                       !
!--------------------------------------------------------------------------------------------------!

!!* Contains some widely used types (at the moment only TOrbitals)
module commontypes
#include "allocate.h"  
  use accuracy
  implicit none
  private

  public :: destruct
  public :: TOrbitals


  !!* Contains information about the orbitals of the species/atoms in the system
  type TOrbitals
    !!* Nr. of shells for each atomic species (nSpecies)
    integer, pointer :: nShell(:) => null()      
    !!* Nr. of orbitals for each atomic species (nSpecies)
    integer, pointer :: nOrbSpecies(:) => null()  
    !!* Nr. of orbitals for each atom (nAtom)
    integer, pointer :: nOrbAtom(:) => null()
    !!* Ang. momentum of the a particular l-shell on a particular species
    !!* (maxval(nShell), nSpecies)
    integer, pointer :: angShell(:,:) => null() 
    !!* The shell which contains the given orbital on an atom
    !!* (maxval(nOrbSpecies), nSpecies)
    integer, pointer :: iShellOrb(:,:) => null()
    !!* Starting pos. within the atomic block of the each of the shells of
    !!* each species (maxval(nShell)+1, nSpecies)
    integer, pointer :: posShell(:,:) => null()
    integer :: mShell                   !* Max. nr. of shells for any species
    integer :: mOrb                     !* Max. nr. of orbitals for any species
    integer :: nOrb                     !* Total number of orbitals in system.
  end type TOrbitals


  !!* Destroys the components of a TOrbitals instance
  interface destruct
    module procedure Orbitals_destruct
  end interface


contains

  !!* Destroys the components of a TOrbitals instance
  !!* @param self Instance.
  subroutine Orbitals_destruct(self)
    type(TOrbitals), intent(inout) :: self

    DEALLOCATE_PARR(self%nShell)
    DEALLOCATE_PARR(self%nOrbSpecies)
    DEALLOCATE_PARR(self%nOrbAtom)
    DEALLOCATE_PARR(self%angShell)
    DEALLOCATE_PARR(self%iShellOrb)
    DEALLOCATE_PARR(self%posShell)

  end subroutine Orbitals_destruct
  

end module commontypes
