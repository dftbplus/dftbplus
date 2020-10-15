!--------------------------------------------------------------------------------------------------!
!  DFTB+: general package for performing fast atomistic simulations                                !
!  Copyright (C) 2006 - 2020  DFTB+ developers group                                               !
!                                                                                                  !
!  See the LICENSE file for terms of usage and distribution.                                       !
!--------------------------------------------------------------------------------------------------!

!> Contains a type with basis information
module dftbp_orbitals
  use dftbp_accuracy
  use dftbp_message
  use dftbp_constants, only : shellNames
  implicit none
  private

  public :: TOrbitals, getShellNames, orbitalNames


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

  !> Length of labels for atomic orbitals
  integer, parameter :: lenOrbitalNames = 9

  !> Names of the atomic orbitals in tesseral basis
  !> general set for f orbitals (not cubic), see
  !> http://winter.group.shef.ac.uk/orbitron/AOs/4f/equations.html
  character(lenOrbitalNames), parameter :: orbitalNames(-3:3,0:3) = reshape([&
      & '         ','         ','         ','         ','         ','         ','         ',&
      & '         ','         ','y        ','z        ','x        ','         ','         ',&
      & '         ','xy       ','yz       ','z2       ','xz       ','x2-y2    ','         ',&
      & 'y(3x2-y2)','x2+y2+z2 ','yz2      ','z3       ','xz2      ','z(x2-y2) ','x(x2-3y2)'&
      &], [7,4])

contains

  !> Builds a unique names for the atomic orbitals
  !> Assign 's', 'p', 'd' to first occurring shells, then 's2', 'p2', ...
  subroutine getShellNames(iSpecie, orb, shellNamesTmp)

    !> atomic specie
    integer, intent(in) :: iSpecie

    !> orbital info
    type(TOrbitals), intent(in) :: orb

    !> output string naming the atomic orbital
    character(sc), intent(out), allocatable :: shellNamesTmp(:)

    integer :: ii
    integer, allocatable :: ind(:)
    character(sc) :: sindx

    !allocate(names(orb%nShell(iSpecie)))
    allocate(shellNamesTmp(orb%nShell(iSpecie)))
    allocate(ind(orb%nShell(iSpecie)))
    ind = 1

    do ii = 1, orb%nShell(iSpecie)
      write(shellNamesTmp(ii), "(A)") shellNames(orb%angShell(ii, iSpecie) + 1)
      if (any(shellNamesTmp(ii) == shellNamesTmp(1:ii-1))) then
        ! at least one example of this shell already
        ind(ii) = ind(ii) + 1
        write(sindx,'(I0)') ind(ii)
        if (len(trim(adjustl(shellNamesTmp(ii)))) + len(trim(sindx)) > sc) then
          call error("Shell labels are too long: "//trim(adjustl(shellNamesTmp(ii)))//trim(sindx))
        else
          shellNamesTmp(ii) = trim(adjustl(shellNamesTmp(ii)))//trim(sindx)
        end if
      end if
    end do
    deallocate(ind)

  end subroutine getShellNames


end module dftbp_orbitals
