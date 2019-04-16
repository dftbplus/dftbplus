!--------------------------------------------------------------------------------------------------!
!  DFTB+: general package for performing fast atomistic simulations                                !
!  Copyright (C) 2018  DFTB+ developers group                                                      !
!                                                                                                  !
!  See the LICENSE file for terms of usage and distribution.                                       !
!--------------------------------------------------------------------------------------------------!

!> Contains a type with basis information
module dftbp_orbitals
  use dftbp_accuracy
  use dftbp_message
  use dftbp_constants, only : orbitalNames
  implicit none
  private

  public :: TOrbitals, getOrbitalNames


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


contains

  !> Builds a unique names for the atomic orbitals
  !> Assign 's', 'p', 'd' to first occurring shells, then 's2', 'p2', ...
  subroutine getOrbitalNames(iSpecie, orb, shellnames)

    !> atomic specie
    integer, intent(in) :: iSpecie

    !> orbital info
    type(TOrbitals), intent(in) :: orb

    !> output string naming the atomic orbital
    character(sc), intent(out), allocatable :: shellnames(:)

    integer :: ii, jj
    integer, allocatable :: ind(:)
    character(sc) :: sindx

    !allocate(names(orb%nShell(iSpecie)))
    allocate(shellnames(orb%nShell(iSpecie)))
    allocate(ind(orb%nShell(iSpecie)))
    ind = 1

    do ii = 1, orb%nShell(iSpecie)
      write(shellnames(ii), "(A)") orbitalNames(orb%angShell(ii, iSpecie) + 1)
      if (any(shellnames(ii) == shellnames(1:ii-1))) then
        ind(ii) = ind(ii) + 1
        if (ind(ii) == 1) then
          sindx = ""
        else
          write(sindx,'(I0)') ind(ii)
        end if
        if (len(trim(adjustl(shellnames(ii)))) + len(trim(sindx)) > sc) then
          call error("Shell labels are too long: "//trim(adjustl(shellnames(ii)))//trim(sindx))
        else
          shellnames(ii) = trim(adjustl(shellnames(ii)))//trim(sindx)
        end if
      end if
    end do
    deallocate(ind)

  end subroutine getOrbitalNames


end module dftbp_orbitals
