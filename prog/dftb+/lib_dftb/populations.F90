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

  public :: mulliken, skewMulliken, denseSubtractDensityOfAtoms, getChargePerShell


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

  
  !> Interface to subtract superposition of atomic densities from dense density matrix.
  !> Required for rangeseparated calculations
  interface denseSubtractDensityOfAtoms
     module procedure denseSubtractDensityOfAtoms_nospin
     module procedure denseSubtractDensityOfAtoms_spin
  end interface denseSubtractDensityOfAtoms

contains


  !> Calculate the Mulliken population for each atom in the system, by summing
  !> the individual orbital contriutions on each atom
  subroutine mullikenPerAtom(qq, over, rho, orb, iNeighbor, nNeighbor, img2CentCell, iPair)

    !> The charge per atom
    real(dp), intent(inout) :: qq(:)

    !> Overlap matrix in packed format
    real(dp), intent(in) :: over(:)

    !> Density matrix in Packed format
    real(dp), intent(in) :: rho(:)

    !> Information about the orbitals.
    type(TOrbitals), intent(in) :: orb

    !> Number of neighbours of each real atom (central cell)
    integer, intent(in) :: iNeighbor(0:,:)

    !> List of neighbours for each atom, starting at 0 for itself
    integer, intent(in) :: nNeighbor(:)

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

    call mullikenPerOrbital( qPerOrbital,over,rho,orb,iNeighbor,nNeighbor, &
        &img2CentCell,iPair )

    qq(:) = qq(:) + sum(qPerOrbital, dim=1)
    deallocate(qPerOrbital)

  end subroutine mullikenPerAtom


  !> Calculate the Mulliken population for each orbital in the system using purely real-space
  !> overlap and density matrix values.  Currently Mulliken is transformed into real space sums over
  !> one triangle of real space extended matrices
  !>
  !> To do: add description of algorithm to programer manual / documentation.
  subroutine mullikenPerOrbital(qq, over, rho, orb, iNeighbor, nNeighbor, img2CentCell, iPair)

    !> The charge per orbital
    real(dp), intent(inout) :: qq(:,:)

    !> Overlap matrix in packed format
    real(dp), intent(in) :: over(:)

    !> Density matrix in Packed format
    real(dp), intent(in) :: rho(:)

    !> Information about the orbitals.
    type(TOrbitals), intent(in) :: orb

    !> Number of neighbours of each real atom (central cell)
    integer, intent(in) :: iNeighbor(0:,:)

    !> List of neighbours for each atom, starting at 0 for itself
    integer, intent(in) :: nNeighbor(:)

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
      do iNeigh = 0, nNeighbor(iAtom1)
        sqrTmp(:,:) = 0.0_dp
        mulTmp(:) = 0.0_dp
        iAtom2 = iNeighbor(iNeigh, iAtom1)
        iAtom2f = img2CentCell(iAtom2)
        nOrb2 = orb%nOrbAtom(iAtom2f)
        iOrig = iPair(iNeigh,iAtom1) + 1
        mulTmp(1:nOrb1*nOrb2) = over(iOrig:iOrig+nOrb1*nOrb2-1) &
            &* rho(iOrig:iOrig+nOrb1*nOrb2-1)
        sqrTmp(1:nOrb2,1:nOrb1) = &
            & reshape(mulTmp(1:nOrb1*nOrb2), (/nOrb2,nOrb1/))
        qq(1:nOrb1,iAtom1) = qq(1:nOrb1,iAtom1) &
            &+ sum(sqrTmp(1:nOrb2,1:nOrb1), dim=1)
        ! Add contribution to the other triangle sum, using the symmetry
        ! but only when off diagonal
        if (iAtom1 /= iAtom2f) then
          qq(1:nOrb2,iAtom2f) = qq(1:nOrb2,iAtom2f) &
              &+ sum(sqrTmp(1:nOrb2,1:nOrb1), dim=2)
        end if
      end do
    end do

  end subroutine mullikenPerOrbital


  !> Calculate the Mulliken population for each element of the dual atomic blocks in the system
  !> using purely real-space overlap and density matrix values.
  !>
  !> To do: add description of algorithm to programer manual / documentation.
  subroutine mullikenPerBlock(qq, over, rho, orb, iNeighbor, nNeighbor, img2CentCell, iPair)

    !> The charge per atom block
    real(dp), intent(inout) :: qq(:,:,:)

    !> Overlap matrix in packed format
    real(dp), intent(in) :: over(:)

    !> Density matrix in Packed format
    real(dp), intent(in) :: rho(:)

    !> Information about the orbitals.
    type(TOrbitals), intent(in) :: orb

    !> Number of neighbours of each real atom (central cell)
    integer, intent(in) :: iNeighbor(0:,:)

    !> List of neighbours for each atom, starting at 0 for itself
    integer, intent(in) :: nNeighbor(:)

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
      do iNeigh = 0, nNeighbor(iAtom1)
        iAtom2 = iNeighbor(iNeigh, iAtom1)
        iAtom2f = img2CentCell(iAtom2)
        nOrb2 = orb%nOrbAtom(iAtom2f)
        iOrig = iPair(iNeigh,iAtom1) + 1
        sTmp(1:nOrb2,1:nOrb1) = &
            & reshape(over(iOrig:iOrig+nOrb1*nOrb2-1), (/nOrb2,nOrb1/))
        rhoTmp(1:nOrb2,1:nOrb1) = &
            & reshape(rho(iOrig:iOrig+nOrb1*nOrb2-1), (/nOrb2,nOrb1/))
        qq(1:nOrb1,1:nOrb1,iAtom1) = qq(1:nOrb1,1:nOrb1,iAtom1) &
            & + 0.5_dp*( matmul(transpose(rhoTmp(1:nOrb2,1:nOrb1)), &
            & sTmp(1:nOrb2,1:nOrb1) ) &
            & + matmul(transpose(sTmp(1:nOrb2,1:nOrb1)), &
            & rhoTmp(1:nOrb2,1:nOrb1) ) )
        ! Add contribution to the other triangle sum, using the symmetry but only when off diagonal
        if (iAtom1 /= iAtom2f) then
          qq(1:nOrb2,1:nOrb2,iAtom2f) = qq(1:nOrb2,1:nOrb2,iAtom2f) &
              & + 0.5_dp*( matmul( rhoTmp(1:nOrb2,1:nOrb1), &
              & transpose(sTmp(1:nOrb2,1:nOrb1)) ) &
              & + matmul( sTmp(1:nOrb2,1:nOrb1), &
              & transpose(rhoTmp(1:nOrb2,1:nOrb1)) ) )
        end if
      end do
    end do

  end subroutine mullikenPerBlock


  !> Calculate the Mulliken population for each element of the dual atomic block orbital in the
  !> system using purely real-space overlap and density matrix values.
  !>
  !> To do: add description of algorithm to programer manual / documentation.
  subroutine skewMullikenPerBlock(qq, over, rho, orb, iNeighbor, nNeighbor, img2CentCell, iPair)

    !> The charge per atom block
    real(dp), intent(inout) :: qq(:,:,:)

    !> Overlap matrix in packed format
    real(dp), intent(in) :: over(:)

    !> Density matrix in Packed format
    real(dp), intent(in) :: rho(:)

    !> Information about the orbitals.
    type(TOrbitals), intent(in) :: orb

    !> Number of neighbours of each real atom (central cell)
    integer, intent(in) :: iNeighbor(0:,:)

    !> List of neighbours for each atom, starting at 0 for itself
    integer, intent(in) :: nNeighbor(:)

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
      do iNeigh = 0, nNeighbor(iAtom1)
        iAtom2 = iNeighbor(iNeigh, iAtom1)
        iAtom2f = img2CentCell(iAtom2)
        nOrb2 = orb%nOrbAtom(iAtom2f)
        iOrig = iPair(iNeigh,iAtom1) + 1
        sTmp(1:nOrb2,1:nOrb1) = &
            & reshape(over(iOrig:iOrig+nOrb1*nOrb2-1), (/nOrb2,nOrb1/))
        rhoTmp(1:nOrb2,1:nOrb1) = &
            & reshape(rho(iOrig:iOrig+nOrb1*nOrb2-1), (/nOrb2,nOrb1/))
        qq(1:nOrb1,1:nOrb1,iAtom1) = qq(1:nOrb1,1:nOrb1,iAtom1) &
            & - 0.5_dp*( matmul(transpose(rhoTmp(1:nOrb2,1:nOrb1)), &
            & sTmp(1:nOrb2,1:nOrb1) ) &
            & - matmul(transpose(sTmp(1:nOrb2,1:nOrb1)), &
            & rhoTmp(1:nOrb2,1:nOrb1) ) )
        ! Add contribution to the other triangle sum, using the symmetry but only when off diagonal
        if (iAtom1 /= iAtom2f) then
          qq(1:nOrb2,1:nOrb2,iAtom2f) = qq(1:nOrb2,1:nOrb2,iAtom2f) &
              & + 0.5_dp*( matmul( rhoTmp(1:nOrb2,1:nOrb1), &
              & transpose(sTmp(1:nOrb2,1:nOrb1)) ) &
              & - matmul( sTmp(1:nOrb2,1:nOrb1), &
              & transpose(rhoTmp(1:nOrb2,1:nOrb1)) ) )
        end if
      end do
    end do

  end subroutine skewMullikenPerBlock


  !> Subtracts superposition of atomic densities from dense density matrix.
  !> Works only for closed shell!
  subroutine denseSubtractDensityOfAtoms_nospin(q0, iSquare, rho)

    !> Reference atom populations
    real(dp), intent(in) :: q0(:,:,:)

    !> Atom positions in the row/column of square matrices
    integer, intent(in) :: iSquare(:)

    !>Spin polarized (lower triangular) density matrix
    real(dp), intent(inout) :: rho(:,:,:)

    integer :: nAtom, iAtom, nSpin, iStart, iEnd, iOrb, iSpin

    nAtom = size(iSquare) - 1
    nSpin = size(rho, dim=3)
    do iSpin = 1, nSpin
       do iAtom = 1, nAtom
          iStart = iSquare(iAtom)
          iEnd = iSquare(iAtom + 1) - 1
          do iOrb = 1, iEnd - iStart + 1
             rho(iStart+iOrb-1, iStart+iOrb-1, iSpin) = &
                  & rho(iStart+iOrb-1, iStart+iOrb-1, iSpin)&
                  & - q0(iOrb, iAtom, iSpin)
          end do
       end do
    end do

  end subroutine denseSubtractDensityOfAtoms_nospin


  !> Subtracts superposition of atomic densities from dense density matrix.
  !> The spin unrestricted version
  !> RangeSep: for spin-unrestricted calculation 
  !> the initial guess should be equally distributed to
  !> alpha and beta density matrices 
  subroutine denseSubtractDensityOfAtoms_spin(q0, iSquare, rho, iSpin)
 
    !> Rerence atom populations
    real(dp), intent(in) :: q0(:,:,:)

    !> Atom positions in the row/colum of square matrix
    integer, intent(in) :: iSquare(:)

    !> Spin polarized (lower triangular) matrix
    real(dp), intent(inout) :: rho(:,:,:)

    !> Spin index
    integer, intent(in) :: iSpin

    integer :: nAtom, iAtom, nSpin, iStart, iEnd, iOrb

    nAtom = size(iSquare) - 1
    nSpin = size(rho, dim=3)
    
    do iAtom = 1, nAtom
       iStart = iSquare(iAtom)
       iEnd = iSquare(iAtom + 1) - 1
       do iOrb = 1, iEnd - iStart + 1
          rho(iStart+iOrb-1, iStart+iOrb-1, iSpin) = &
               & rho(iStart+iOrb-1, iStart+iOrb-1, iSpin)&
               & - q0(iOrb, iAtom, 1)*0.5_dp
       end do
    end do


  end subroutine denseSubtractDensityOfAtoms_spin


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
