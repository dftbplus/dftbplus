!--------------------------------------------------------------------------------------------------!
!  DFTB+: general package for performing fast atomistic simulations                                !
!  Copyright (C) 2018  DFTB+ developers group                                                      !
!                                                                                                  !
!  See the LICENSE file for terms of usage and distribution.                                       !
!--------------------------------------------------------------------------------------------------!

#:include 'common.fypp'

!> Calculates various types of charge populations
!> To do: extend to other populations than Mulliken
module populations
  use assert
  use accuracy
  use constants
  use periodic
  use commontypes
  implicit none
  private

  public :: mulliken, skewMulliken, getChargePerShell


  !> Provides an interface to calculate Mulliken populations, either dual basis atomic block,
  !> orbitally resolved or atom resolved
  interface mulliken
    module procedure mullikenPerBlock
    module procedure mullikenPerOrbital
    module procedure mullikenPerAtom
  end interface mulliken


  !> Provides an interface to calculate Mulliken populations for anti-symmetric density matrices
  interface skewMulliken
    module procedure skewMullikenPerBlock
  end interface skewMulliken

contains


  !> Calculate the Mulliken population for each atom in the system, by summing
  !> the individual orbital contriutions on each atom
  subroutine mullikenPerAtom(qq, over, rho, orb, iNeighbour, nNeighbourSK, img2CentCell, iPair)

    !> The charge per atom
    real(dp), intent(inout) :: qq(:)

    !> Overlap matrix in packed format
    real(dp), intent(in) :: over(:)

    !> Density matrix in Packed format
    real(dp), intent(in) :: rho(:)

    !> Information about the orbitals.
    type(TOrbitals), intent(in) :: orb

    !> Number of neighbours of each real atom (central cell)
    integer, intent(in) :: iNeighbour(0:,:)

    !> List of neighbours for each atom, starting at 0 for itself
    integer, intent(in) :: nNeighbourSK(:)

    !> indexing array to convert images of atoms back into their number in the central cell
    integer, intent(in) :: img2CentCell(:)

    !> indexing array for the Hamiltonian
    integer, intent(in) :: iPair(0:,:)

    real(dp), allocatable :: qPerOrbital(:,:)
    integer :: nAtom

    nAtom = size(orb%nOrbAtom)
    @:ASSERT(size(qq) == nAtom)
    @:ASSERT(size(over) == size(rho))

    allocate(qPerOrbital(orb%mOrb, nAtom))
    qPerOrbital(:,:) = 0.0_dp

    call mullikenPerOrbital( qPerOrbital,over,rho,orb,iNeighbour,nNeighbourSK, img2CentCell,iPair )

    qq(:) = qq(:) + sum(qPerOrbital, dim=1)
    deallocate(qPerOrbital)

  end subroutine mullikenPerAtom


  !> Calculate the Mulliken population for each orbital in the system using purely real-space
  !> overlap and density matrix values.  Currently Mulliken is transformed into real space sums over
  !> one triangle of real space extended matrices
  !>
  !> To do: add description of algorithm to programer manual / documentation.
  subroutine mullikenPerOrbital(qq, over, rho, orb, iNeighbour, nNeighbourSK, img2CentCell, iPair)

    !> The charge per orbital
    real(dp), intent(inout) :: qq(:,:)

    !> Overlap matrix in packed format
    real(dp), intent(in) :: over(:)

    !> Density matrix in Packed format
    real(dp), intent(in) :: rho(:)

    !> Information about the orbitals.
    type(TOrbitals), intent(in) :: orb

    !> Number of neighbours of each real atom (central cell)
    integer, intent(in) :: iNeighbour(0:,:)

    !> List of neighbours for each atom, starting at 0 for itself
    integer, intent(in) :: nNeighbourSK(:)

    !> indexing array to convert images of atoms back into their number in the central cell
    integer, intent(in) :: img2CentCell(:)

    !> indexing array for the Hamiltonian
    integer, intent(in) :: iPair(0:,:)

    integer :: iOrig
    integer :: iNeigh
    integer :: nAtom, iAtom1, iAtom2, iAtom2f
    integer :: nOrb1, nOrb2
    real(dp) :: sqrTmp(orb%mOrb,orb%mOrb)
    real(dp) :: mulTmp(orb%mOrb**2)

    nAtom = size(orb%nOrbAtom)

    @:ASSERT(all(shape(qq) == (/orb%mOrb, nAtom/)))
    @:ASSERT(size(over) == size(rho))

    do iAtom1 = 1, nAtom
      nOrb1 = orb%nOrbAtom(iAtom1)
      do iNeigh = 0, nNeighbourSK(iAtom1)
        sqrTmp(:,:) = 0.0_dp
        mulTmp(:) = 0.0_dp
        iAtom2 = iNeighbour(iNeigh, iAtom1)
        iAtom2f = img2CentCell(iAtom2)
        nOrb2 = orb%nOrbAtom(iAtom2f)
        iOrig = iPair(iNeigh,iAtom1) + 1
        mulTmp(1:nOrb1*nOrb2) = over(iOrig:iOrig+nOrb1*nOrb2-1) * rho(iOrig:iOrig+nOrb1*nOrb2-1)
        sqrTmp(1:nOrb2,1:nOrb1) = reshape(mulTmp(1:nOrb1*nOrb2), (/nOrb2,nOrb1/))
        qq(1:nOrb1,iAtom1) = qq(1:nOrb1,iAtom1) + sum(sqrTmp(1:nOrb2,1:nOrb1), dim=1)
        ! Add contribution to the other triangle sum, using the symmetry
        ! but only when off diagonal
        if (iAtom1 /= iAtom2f) then
          qq(1:nOrb2,iAtom2f) = qq(1:nOrb2,iAtom2f) + sum(sqrTmp(1:nOrb2,1:nOrb1), dim=2)
        end if
      end do
    end do

  end subroutine mullikenPerOrbital


  !> Calculate the Mulliken population for each element of the dual atomic blocks in the system
  !> using purely real-space overlap and density matrix values.
  !>
  !> To do: add description of algorithm to programer manual / documentation.
  subroutine mullikenPerBlock(qq, over, rho, orb, iNeighbour, nNeighbourSK, img2CentCell, iPair)

    !> The charge per atom block
    real(dp), intent(inout) :: qq(:,:,:)

    !> Overlap matrix in packed format
    real(dp), intent(in) :: over(:)

    !> Density matrix in Packed format
    real(dp), intent(in) :: rho(:)

    !> Information about the orbitals.
    type(TOrbitals), intent(in) :: orb

    !> Number of neighbours of each real atom (central cell)
    integer, intent(in) :: iNeighbour(0:,:)

    !> List of neighbours for each atom, starting at 0 for itself
    integer, intent(in) :: nNeighbourSK(:)

    !> indexing array to convert images of atoms back into their number in the central cell
    integer, intent(in) :: img2CentCell(:)

    !> indexing array for the Hamiltonian
    integer, intent(in) :: iPair(0:,:)

    integer :: iOrig
    integer :: iNeigh
    integer :: nAtom, iAtom1, iAtom2, iAtom2f
    integer :: nOrb1, nOrb2
    real(dp) :: STmp(orb%mOrb,orb%mOrb)
    real(dp) :: rhoTmp(orb%mOrb,orb%mOrb)

    nAtom = size(orb%nOrbAtom)

    @:ASSERT(all(shape(qq) == (/orb%mOrb,orb%mOrb,nAtom/)))
    @:ASSERT(size(over) == size(rho))

    do iAtom1 = 1, nAtom
      nOrb1 = orb%nOrbAtom(iAtom1)
      do iNeigh = 0, nNeighbourSK(iAtom1)
        iAtom2 = iNeighbour(iNeigh, iAtom1)
        iAtom2f = img2CentCell(iAtom2)
        nOrb2 = orb%nOrbAtom(iAtom2f)
        iOrig = iPair(iNeigh,iAtom1) + 1
        sTmp(1:nOrb2,1:nOrb1) = reshape(over(iOrig:iOrig+nOrb1*nOrb2-1), (/nOrb2,nOrb1/))
        rhoTmp(1:nOrb2,1:nOrb1) = reshape(rho(iOrig:iOrig+nOrb1*nOrb2-1), (/nOrb2,nOrb1/))
        qq(1:nOrb1,1:nOrb1,iAtom1) = qq(1:nOrb1,1:nOrb1,iAtom1)&
            & + 0.5_dp*( matmul(transpose(rhoTmp(1:nOrb2,1:nOrb1)), sTmp(1:nOrb2,1:nOrb1) )&
            & + matmul(transpose(sTmp(1:nOrb2,1:nOrb1)), rhoTmp(1:nOrb2,1:nOrb1) ) )
        ! Add contribution to the other triangle sum, using the symmetry but only when off diagonal
        if (iAtom1 /= iAtom2f) then
          qq(1:nOrb2,1:nOrb2,iAtom2f) = qq(1:nOrb2,1:nOrb2,iAtom2f)&
              & + 0.5_dp*( matmul( rhoTmp(1:nOrb2,1:nOrb1), transpose(sTmp(1:nOrb2,1:nOrb1)) )&
              & + matmul( sTmp(1:nOrb2,1:nOrb1), transpose(rhoTmp(1:nOrb2,1:nOrb1)) ) )
        end if
      end do
    end do

  end subroutine mullikenPerBlock


  !> Calculate the Mulliken population for each element of the dual atomic block orbital in the
  !> system using purely real-space overlap and density matrix values.
  !>
  !> To do: add description of algorithm to programer manual / documentation.
  subroutine skewMullikenPerBlock(qq, over, rho, orb, iNeighbour, nNeighbourSK, img2CentCell, iPair)

    !> The charge per atom block
    real(dp), intent(inout) :: qq(:,:,:)

    !> Overlap matrix in packed format
    real(dp), intent(in) :: over(:)

    !> Density matrix in Packed format
    real(dp), intent(in) :: rho(:)

    !> Information about the orbitals.
    type(TOrbitals), intent(in) :: orb

    !> Number of neighbours of each real atom (central cell)
    integer, intent(in) :: iNeighbour(0:,:)

    !> List of neighbours for each atom, starting at 0 for itself
    integer, intent(in) :: nNeighbourSK(:)

    !> indexing array to convert images of atoms back into their number in the central cell
    integer, intent(in) :: img2CentCell(:)

    !> indexing array for the Hamiltonian
    integer, intent(in) :: iPair(0:,:)

    integer :: iOrig
    integer :: iNeigh
    integer :: nAtom, iAtom1, iAtom2, iAtom2f
    integer :: nOrb1, nOrb2
    real(dp) :: STmp(orb%mOrb,orb%mOrb)
    real(dp) :: rhoTmp(orb%mOrb,orb%mOrb)

    nAtom = size(orb%nOrbAtom)

    @:ASSERT(all(shape(qq) == (/orb%mOrb,orb%mOrb,nAtom/)))
    @:ASSERT(size(over) == size(rho))

    do iAtom1 = 1, nAtom
      nOrb1 = orb%nOrbAtom(iAtom1)
      do iNeigh = 0, nNeighbourSK(iAtom1)
        iAtom2 = iNeighbour(iNeigh, iAtom1)
        iAtom2f = img2CentCell(iAtom2)
        nOrb2 = orb%nOrbAtom(iAtom2f)
        iOrig = iPair(iNeigh,iAtom1) + 1
        sTmp(1:nOrb2,1:nOrb1) = reshape(over(iOrig:iOrig+nOrb1*nOrb2-1), (/nOrb2,nOrb1/))
        rhoTmp(1:nOrb2,1:nOrb1) = reshape(rho(iOrig:iOrig+nOrb1*nOrb2-1), (/nOrb2,nOrb1/))
        qq(1:nOrb1,1:nOrb1,iAtom1) = qq(1:nOrb1,1:nOrb1,iAtom1)&
            & - 0.5_dp*( matmul(transpose(rhoTmp(1:nOrb2,1:nOrb1)), sTmp(1:nOrb2,1:nOrb1) )&
            & - matmul(transpose(sTmp(1:nOrb2,1:nOrb1)), rhoTmp(1:nOrb2,1:nOrb1) ) )
        ! Add contribution to the other triangle sum, using the symmetry but only when off diagonal
        if (iAtom1 /= iAtom2f) then
          qq(1:nOrb2,1:nOrb2,iAtom2f) = qq(1:nOrb2,1:nOrb2,iAtom2f)&
              & + 0.5_dp*( matmul( rhoTmp(1:nOrb2,1:nOrb1), transpose(sTmp(1:nOrb2,1:nOrb1)) )&
              & - matmul( sTmp(1:nOrb2,1:nOrb1), transpose(rhoTmp(1:nOrb2,1:nOrb1)) ) )
        end if
      end do
    end do

  end subroutine skewMullikenPerBlock


  !> Calculate the number of charges per shell given the orbital charges.
  subroutine getChargePerShell(qq, orb, species, chargePerShell)

    !> charges in each orbital, for each atom and spin channel
    real(dp), intent(in) :: qq(:,:,:)

    !> orbital information
    type(TOrbitals), intent(in) :: orb

    !> species of each atom
    integer, intent(in) :: species(:)

    !> Resulting charges in atomic shells
    real(dp), intent(out) :: chargePerShell(:,:,:)

    integer :: iAt, iSp, iSh
    integer :: nAtom, nSpin

    nAtom = size(chargePerShell, dim=2)
    nSpin = size(chargePerShell, dim=3)
    chargePerShell(:,:,:) = 0.0_dp
    do iAt = 1, nAtom
      iSp = species(iAt)
      do iSh = 1, orb%nShell(iSp)
        chargePerShell(iSh, iAt, 1:nSpin) = chargePerShell(iSh, iAt, 1:nSpin)&
            & + sum(qq(orb%posShell(iSh, iSp) : orb%posShell(iSh + 1, iSp) - 1, iAt, 1:nSpin),&
            & dim=1)
      end do
    end do

  end subroutine getChargePerShell



end module populations
