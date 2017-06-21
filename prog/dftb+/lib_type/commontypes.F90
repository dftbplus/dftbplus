!--------------------------------------------------------------------------------------------------!
!  DFTB+: general package for performing fast atomistic simulations                                !
!  Copyright (C) 2017  DFTB+ developers group                                                      !
!                                                                                                  !
!  See the LICENSE file for terms of usage and distribution.                                       !
!--------------------------------------------------------------------------------------------------!

!!* Contains some widely used types (at the moment only TOrbitals)
module commontypes
  use accuracy
  implicit none
  private

  public :: TOrbitals


  !!* Contains information about the orbitals of the species/atoms in the system
  type TOrbitals
    !!* Nr. of shells for each atomic species (nSpecies)
    integer, allocatable :: nShell(:)
    !!* Nr. of orbitals for each atomic species (nSpecies)
    integer, allocatable :: nOrbSpecies(:)
    !!* Nr. of orbitals for each atom (nAtom)
    integer, allocatable :: nOrbAtom(:)
    !!* Ang. momentum of the a particular l-shell on a particular species
    !!* (maxval(nShell), nSpecies)
    integer, allocatable :: angShell(:,:)
    !!* The shell which contains the given orbital on an atom
    !!* (maxval(nOrbSpecies), nSpecies)
    integer, allocatable :: iShellOrb(:,:)
    !!* Starting pos. within the atomic block of the each of the shells of
    !!* each species (maxval(nShell)+1, nSpecies)
    integer, allocatable :: posShell(:,:)
    integer :: mShell                   !* Max. nr. of shells for any species
    integer :: mOrb                     !* Max. nr. of orbitals for any species
    integer :: nOrb                     !* Total number of orbitals in system.
  end type TOrbitals

end module commontypes
