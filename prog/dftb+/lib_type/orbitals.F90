!--------------------------------------------------------------------------------------------------!
!  DFTB+: general package for performing fast atomistic simulations                                !
!  Copyright (C) 2018  DFTB+ developers group                                                      !
!                                                                                                  !
!  See the LICENSE file for terms of usage and distribution.                                       !
!--------------------------------------------------------------------------------------------------!

!> Contains a type with basis information
module orbitals
  use accuracy
  implicit none
  private

  public :: TOrbitals


  !> Contains information about the orbitals of the species/atoms in the system
  type TOrbitals

    !> Nr. of shells for each atomic species (nSpecies)
    integer, allocatable :: nShell(:)

    !> Nr. of orbitals for each atomic species (nSpecies)
    integer, allocatable :: nOrbSpecies(:)

    !> Nr. of orbitals for each atom (nAtom)
    integer, allocatable :: nOrbAtom(:)

    !> Ang. momentum of the a particular l-shell on a particular species (maxval(nShell), nSpecies)
    integer, allocatable :: angShell(:,:)

    !> The shell which contains the given orbital on an atom
    !> (maxval(nOrbSpecies), nSpecies)
    integer, allocatable :: iShellOrb(:,:)

    !> Starting pos. within the atomic block of the each of the shells of each species
    !> (maxval(nShell)+1, nSpecies)
    integer, allocatable :: posShell(:,:)

    !> Max. nr. of shells for any species
    integer :: mShell

    !> Max. nr. of orbitals for any species
    integer :: mOrb

    !> Total number of orbitals in system.
    integer :: nOrb
  end type TOrbitals

end module orbitals
