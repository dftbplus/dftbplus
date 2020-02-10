!--------------------------------------------------------------------------------------------------!
!  DFTB+: general package for performing fast atomistic simulations                                !
!  Copyright (C) 2006 - 2020  DFTB+ developers group                                               !
!                                                                                                  !
!  See the LICENSE file for terms of usage and distribution.                                       !
!--------------------------------------------------------------------------------------------------!

#:include 'common.fypp'

!> Non-selfconsistent part of the extended tight binding Hamiltonian
module dftbp_xtbh0
  use dftbp_accuracy, only : dp
  use dftbp_commontypes, only : TOrbitals
  use dftbp_environment, only : TEnvironment
  use dftbp_gtocont, only : TGTOCont
  implicit none
  private

  public :: buildSH0, xtbSelfEnergy

  !> Maximal angular momentum, for which rotations are present
  integer, parameter :: mAngRot_ = 3

contains

  !> Calculate coordination number dependency of xTB self-energy
  subroutine xtbSelfEnergy(selfEnergy, species, orb, atomEnergy, kcn, cn)

    !> On-site energies for each atom
    real(dp), intent(out) :: selfEnergy(:, :)

    !> Chemical species of each atom
    integer, intent(in) :: species(:)

    !> Information about the orbitals in the system
    type(TOrbitals), intent(in) :: orb

    !> On-site energies for each species
    real(dp), intent(in) :: atomEnergy(:, :)

    !> Enhancement factor for the coordination number dependence
    real(dp), intent(in) :: kcn(:, :)

    !> Coordination number of every atom
    real(dp), intent(in) :: cn(:)

    integer :: nAtom
    integer :: iAt1, iSp1, iSh1

    nAtom = size(selfEnergy, dim=2)

    @:ASSERT(size(atomEnergy, dim=1) == size(selfEnergy, dim=1))
    @:ASSERT(size(atomEnergy, dim=2) == maxval(species))

    !$omp parallel do schedule(runtime) default(none) &
    !$omp shared(nAtom, species, orb, selfEnergy, atomEnergy, cn, kcn) &
    !$omp private(iAt1, iSp1, iSh1)
    do iAt1 = 1, nAtom
      iSp1 = species(iAt1)
      do iSh1 = 1, orb%nShell(iSp1)
        selfEnergy(iSh1, iAt1) = atomEnergy(iSh1, iSp1) + cn(iAt1) * kcn(iSh1, iSp1)
      end do
    end do
    !$omp end parallel do

  end subroutine xtbSelfEnergy
 

  subroutine buildSH0(env, ovlp, ham, gtoCont, selfEnergy, coords, nNeighbour, &
      & iNeighbours, species, iPair, orb)

    !> Computational environment settings
    type(TEnvironment), intent(in) :: env

    !> Returns the non-self-consistent Hamiltonian
    real(dp), intent(out) :: ovlp(:)

    !> Returns the non-self-consistent Hamiltonian
    real(dp), intent(out) :: ham(:)

    !> Container for the SlaKo Hamiltonian integrals
    type(TGTOCont), intent(in) :: gtoCont

    !> On-site energies for each atom
    real(dp), intent(in) :: selfEnergy(:, :)

    !> List of all coordinates, including possible periodic images of atoms
    real(dp), intent(in) :: coords(:, :)

    !> Number of surrounding neighbours for each atom
    integer, intent(in) :: nNeighbour(:)

    !> List of surrounding neighbours for each atom
    integer, intent(in) :: iNeighbours(0:, :)

    !> Chemical species of each atom
    integer, intent(in) :: species(:)

    !> Shift vector, where the interaction between two atoms
    integer, intent(in) :: iPair(0:, :)

    !> Information about the orbitals in the system
    type(TOrbitals), intent(in) :: orb

    integer :: nAtom, iAt1, iNeigh1, iSp1, iAt2, iSp2, ind
    integer :: iOrb1, nOrb1, nOrb2
    real(dp) :: hDiag
    real(dp) :: tmpH(orb%mOrb, orb%mOrb), tmpS(orb%mOrb, orb%mOrb)
    real(dp) :: vec(3), dist

    nAtom = size(nNeighbour)
    ham(:) = 0.0_dp

    @:ASSERT(size(selfEnergy, dim=2) == nAtom)

    !$omp parallel do schedule(runtime) default(none) &
    !$omp shared(nAtom, species, iPair, orb, selfEnergy, ham) &
    !$omp private(iAt1, iSp1, ind, iOrb1)
    do iAt1 = 1, nAtom
      iSp1 = species(iAt1)
      ind = iPair(0, iAt1) + 1
      do iOrb1 = 1, orb%nOrbAtom(iAt1)
        ham(ind) = selfEnergy(orb%iShellOrb(iOrb1, iSp1), iAt1)
        ind = ind + orb%nOrbAtom(iAt1) + 1
      end do
    end do
    !$omp end parallel do

    ! Do the diatomic blocks for each of the atoms with its nNeighbour
    !$omp parallel do schedule(runtime) default(none) shared(ham, ovlp) &
    !$omp shared(nAtom, species, orb, nNeighbour, iNeighbours, iPair, coords) &
    !$omp private(iAt1, iSp1, nOrb1, iNeigh1, iAt2, iSp2, nOrb2, ind, vec, dist) &
    !$omp private(tmpH, tmpS)
    do iAt1 = 1, nAtom
      iSp1 = species(iAt1)
      nOrb1 = orb%nOrbSpecies(iSp1)
      do iNeigh1 = 1, nNeighbour(iAt1)
        iAt2 = iNeighbours(iNeigh1, iAt1)
        iSp2 = species(iAt2)
        nOrb2 = orb%nOrbSpecies(iSp2)
        ind = iPair(iNeigh1, iAt1)
        vec(:) = coords(:,iAt2) - coords(:,iAt1)
        call getDiatomicShellBlock(tmpS, tmpH, iSp1, iSp2, orb, vec)
        ham(ind + 1 : ind + nOrb2 * nOrb1) = reshape(tmpH(1:nOrb2, 1:nOrb1), [nOrb2 * nOrb1])
        ovlp(ind + 1 : ind + nOrb2 * nOrb1) = reshape(tmpS(1:nOrb2, 1:nOrb1), [nOrb2 * nOrb1])
      end do
    end do
    !$omp end parallel do

  end subroutine buildSH0

  subroutine getDiatomicShellBlock(ss, hh, iSp1, iSp2, orb, vec)

    !> the rectangular matrix containing the resulting diatomic matrix elements
    real(dp), intent(out) :: hh(:,:)

    !> the rectangular matrix containing the resulting diatomic matrix elements
    real(dp), intent(out) :: ss(:,:)

    !> Chemical species of atom i
    integer, intent(in) :: iSp1

    !> chemical species of atom j
    integer, intent(in) :: iSp2

    !> Information about the orbitals of chemical species in the system.
    type(TOrbitals), intent(in) :: orb

    !> distance vector
    real(dp), intent(in) :: vec(3)

    real(dp) :: dist

    integer :: iCol, iRow, ind, iSh1, iSh2
    integer :: ang1, ang2, nOrb1, nOrb2
    real(dp) :: tmpS(2*mAngRot_+1,2*mAngRot_+1)
    real(dp) :: tmpH(2*mAngRot_+1,2*mAngRot_+1)
    dist = norm2(vec)

    ind = 1
    iCol = 1
    do iSh1 = 1, orb%nShell(iSp1)
      ang1 = orb%angShell(iSh1, iSp1)
      nOrb1 = 2 * ang1 + 1
      iRow = 1
      do iSh2 = 1, orb%nShell(iSp2)
        ang2 = orb%angShell(iSh2, iSp2)
        nOrb2 = 2 * ang2 + 1

        call hamiltonianFromOverlap(tmpH, tmpS)

        ss(iRow:iRow+nOrb2-1,iCol:iCol+nOrb1-1) = tmpS(1:nOrb2,1:nOrb1)
        hh(iRow:iRow+nOrb2-1,iCol:iCol+nOrb1-1) = tmpH(1:nOrb2,1:nOrb1)
      end do
    end do
  end subroutine getDiatomicShellBlock

  subroutine hamiltonianFromOverlap(hh, ss)

    !> the rectangular matrix containing the resulting diatomic matrix elements
    real(dp), intent(out) :: hh(:,:)

    !> the rectangular matrix containing the resulting diatomic matrix elements
    real(dp), intent(in) :: ss(:,:)
  end subroutine hamiltonianFromOverlap

end module dftbp_xtbh0
