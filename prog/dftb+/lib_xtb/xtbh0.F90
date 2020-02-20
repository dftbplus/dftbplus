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
  use dftbp_schedule, only : distributeRangeInChunks, assembleChunks
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
 

  subroutine buildSH0(env, ovlp, ham, gtoCont, coords, nNeighbour, &
      & iNeighbours, species, iPair, orb)
    use dftbp_constants

    !> Computational environment settings
    type(TEnvironment), intent(in) :: env

    !> Returns the non-self-consistent Hamiltonian
    real(dp), intent(out) :: ovlp(:)

    !> Returns the non-self-consistent Hamiltonian
    real(dp), intent(out) :: ham(:)

    !> Container for the SlaKo Hamiltonian integrals
    type(TGTOCont), intent(in) :: gtoCont

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
    integer :: iAtFirst, iAtLast
    integer :: iOrb1, nOrb1, nOrb2
    real(dp) :: hDiag
    real(dp) :: tmpH(orb%mOrb, orb%mOrb), tmpS(orb%mOrb, orb%mOrb)
    real(dp) :: vec(3), dist

    nAtom = size(nNeighbour)
    ham(:) = 0.0_dp
    ovlp(:) = 0.0_dp

    @:ASSERT(allocated(gtoCont%selfEnergy))
    @:ASSERT(size(gtoCont%selfEnergy, dim=2) == nAtom)

    call distributeRangeInChunks(env, 1, nAtom, iAtFirst, iAtLast)

    ! Put the on-site energies into the Hamiltonian,
    ! and <lm|l'm'> = delta_l,l' * delta_m',m' for the overlap
    ! omp parallel do schedule(runtime) default(none) &
    ! omp shared(iAtFirst, iAtLast, species, iPair, orb, gtoCont, ham, ovlp) &
    ! omp private(iAt1, iSp1, ind, iOrb1)
    do iAt1 = iAtFirst, iAtLast
      iSp1 = species(iAt1)
      do iOrb1 = 1, orb%nOrbAtom(iAt1)
        ind = iPair(0, iAt1) + 1 + (iOrb1 - 1) * (orb%nOrbAtom(iAt1) + 1)
        print*, iSp1, gtoCont%selfEnergy(orb%iShellOrb(iOrb1, iSp1), iAt1) * Hartree__eV
        ham(ind) = gtoCont%selfEnergy(orb%iShellOrb(iOrb1, iSp1), iAt1)
        ovlp(ind) = 1.0_dp
      end do
    end do
    ! omp end parallel do

    ! Do the diatomic blocks for each of the atoms with its nNeighbour
    !$omp parallel do schedule(runtime) default(none) &
    !$omp shared(iAtFirst, iAtLast, species, orb, nNeighbour, iNeighbours, iPair) &
    !$omp shared(coords, ham, ovlp, gtoCont) &
    !$omp private(iAt1, iSp1, nOrb1, iNeigh1, iAt2, iSp2, nOrb2, ind, vec, dist) &
    !$omp private(tmpH, tmpS)
    do iAt1 = iAtFirst, iAtLast
      iSp1 = species(iAt1)
      nOrb1 = orb%nOrbSpecies(iSp1)
      do iNeigh1 = 1, nNeighbour(iAt1)
        iAt2 = iNeighbours(iNeigh1, iAt1)
        iSp2 = species(iAt2)
        nOrb2 = orb%nOrbSpecies(iSp2)
        ind = iPair(iNeigh1, iAt1)
        vec(:) = coords(:,iAt2) - coords(:,iAt1)
        call getDiatomicShellBlock(tmpS, tmpH, gtoCont, iSp1, iSp2, orb, vec)
        ham(ind+1:ind+nOrb2*nOrb1) = reshape(tmpH(1:nOrb2, 1:nOrb1), [nOrb2*nOrb1])
        ovlp(ind+1:ind+nOrb2*nOrb1) = reshape(tmpS(1:nOrb2, 1:nOrb1), [nOrb2*nOrb1])
      end do
    end do
    !$omp end parallel do

    call assembleChunks(env, ham)
    call assembleChunks(env, ovlp)

  end subroutine buildSH0

  subroutine getDiatomicShellBlock(ss, hh, gtoCont, iSp1, iSp2, orb, vec)

    !> the rectangular matrix containing the resulting diatomic matrix elements
    real(dp), intent(out) :: hh(:,:)

    !> the rectangular matrix containing the resulting diatomic matrix elements
    real(dp), intent(out) :: ss(:,:)

    !> Container for the SlaKo Hamiltonian integrals
    type(TGTOCont), intent(in) :: gtoCont

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
    real(dp) :: h11, h22, h12, poly1, poly2, poly12, rad12, dEN, kEht
    real(dp) :: tmpS(2*mAngRot_+1,2*mAngRot_+1)
    real(dp) :: tmpH(2*mAngRot_+1,2*mAngRot_+1)

    dist = norm2(vec)
    rad12 = gtoCont%atomicRad(iSp1) + gtoCont%atomicRad(iSp2)
    dEN = (gtoCont%electronegativity(iSp1) - gtoCont%electronegativity(iSp2))**2
    dEN = 1.0_dp - gtoCont%kENScale * dEN

    ind = 1
    iCol = 1
    do iSh1 = 1, orb%nShell(iSp1)
      ang1 = orb%angShell(iSh1, iSp1)
      nOrb1 = 2 * ang1 + 1
      h11 = gtoCont%selfEnergy(iSh1, iSp1)
      poly1 = gtoCont%shellPoly(iSh1, iSp1)
      iRow = 1
      do iSh2 = 1, orb%nShell(iSp2)
        ang2 = orb%angShell(iSh2, iSp2)
        nOrb2 = 2 * ang2 + 1
        h22 = gtoCont%selfEnergy(iSh2, iSp2)
        h12 = 0.5_dp * (h11 + h22)
        poly2 = gtoCont%shellPoly(iSh1, iSp1)
        call shellPoly(poly12, dist, rad12, poly1, poly2)

        call gtoCont%getOverlapIntegrals(tmpS, vec, dist, iSp1, iSh1, iSp2, iSh2)
        tmpH(:, :) = tmpS * h12 * poly12 * gtoCont%h0Scale(iSh2, iSh1, iSp2, iSp1)

        ss(iRow:iRow+nOrb2-1,iCol:iCol+nOrb1-1) = tmpS(1:nOrb2,1:nOrb1)
        hh(iRow:iRow+nOrb2-1,iCol:iCol+nOrb1-1) = tmpH(1:nOrb2,1:nOrb1)
        iRow = iRow + nOrb2
      end do
      iCol = iCol + nOrb1
    end do
  end subroutine getDiatomicShellBlock


  !> Enhancement factor for Hamiltonian integrals
  elemental subroutine shellPoly(poly, dist, rad12, poly1, poly2)

    !> Enhancement factor
    real(dp), intent(out) :: poly

    !> Actual distance between atom 1 and 2
    real(dp), intent(in) :: dist

    !> Sum of atomic radii between atom 1 and 2
    real(dp), intent(in) :: rad12

    !> Polynomial coefficient for atom 1
    real(dp), intent(in) :: poly1

    !> Polynomial coefficient for atom 2
    real(dp), intent(in) :: poly2

    real(dp) :: rr, rf1, rf2

    rr = sqrt(dist/rad12)
    rf1 = 1.0_dp + 0.01_dp * poly1 * rr
    rf2 = 1.0_dp + 0.01_dp * poly1 * rr
    poly = rf1 * rf2

  end subroutine shellPoly

end module dftbp_xtbh0
